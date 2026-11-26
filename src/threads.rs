//! Work spread over the machine's threads.

/// Runs `work` on as many threads as the machine offers, at most `most`, and returns what each run
/// returned. `work` takes its items from a queue its runs share until none is left, so however
/// many runs there are, between them they do all of it.
///
/// One run is always the calling thread's, so nothing is spawned where there is one thread to use.
/// A thread that cannot be spawned is done without: a count of processors does not say that a
/// thread can be had, as on a WebAssembly host without threads, and the calling thread's run then
/// does the work alone. A panic in `work` is the caller's, whichever thread it struck.
pub(crate) fn spread<T: Send>(most: usize, work: impl Fn() -> T + Sync) -> Vec<T> {
    let offered = std::thread::available_parallelism().map_or(1, |n| n.get());
    on(offered.min(most), work)
}

/// `f` of each of `0..n`, in order, with the items taken `chunk` at a time by as many runs of
/// [`spread`] as there are chunks. Each value is what `f` gives for its item on whichever thread,
/// so the result is the one a loop over `0..n` gives; one chunk is that loop.
pub(crate) fn map<T: Send>(n: usize, chunk: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let chunk = chunk.max(1);
    let chunks = n.div_ceil(chunk);
    if chunks <= 1 {
        return (0..n).map(f).collect();
    }
    let next = AtomicUsize::new(0);
    let runs = spread(chunks, || {
        let mut mine = Vec::new();
        loop {
            let c = next.fetch_add(1, Ordering::Relaxed);
            if c >= chunks {
                break mine;
            }
            let items = c * chunk..((c + 1) * chunk).min(n);
            mine.push((c, items.map(&f).collect::<Vec<T>>()));
        }
    });
    let mut parts: Vec<(usize, Vec<T>)> = runs.into_iter().flatten().collect();
    parts.sort_unstable_by_key(|p| p.0);
    parts.into_iter().flat_map(|p| p.1).collect()
}

/// [`spread`] over `threads` threads, the calling thread among them, or over as many of them as
/// can be spawned.
fn on<T: Send>(threads: usize, work: impl Fn() -> T + Sync) -> Vec<T> {
    std::thread::scope(|scope| {
        // Spawning stops at the first thread refused: the next would be refused too.
        let spawned: Vec<_> = (1..threads)
            .map_while(|_| std::thread::Builder::new().spawn_scoped(scope, &work).ok())
            .collect();
        let mut out = vec![work()];
        for worker in spawned {
            match worker.join() {
                Ok(done) => out.push(done),
                Err(panic) => std::panic::resume_unwind(panic),
            }
        }
        out
    })
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering};

    use super::on;

    /// The squares of `0..n` in order, each run of the work taking the next number until none is
    /// left, as the kernel's workers take their items; and how many runs there were.
    fn squares(threads: usize, n: usize) -> (Vec<usize>, usize) {
        let next = AtomicUsize::new(0);
        let runs = on(threads, || {
            let mut mine = Vec::new();
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= n {
                    break mine;
                }
                mine.push((i, i * i));
            }
        });
        let count = runs.len();
        let mut all: Vec<(usize, usize)> = runs.into_iter().flatten().collect();
        all.sort_unstable();
        (all.into_iter().map(|p| p.1).collect(), count)
    }

    #[test]
    fn work_gives_the_same_results_on_one_thread_and_on_several() {
        let want: Vec<usize> = (0..1000).map(|i| i * i).collect();
        for threads in [0, 1, 2, 8] {
            let (got, runs) = squares(threads, 1000);
            assert_eq!(got, want, "on {threads} threads");
            // The calling thread's run, and one for each thread spawned: fewer than asked for
            // where a thread is refused, as every one is on a target without threads.
            assert!(
                (1..=threads.max(1)).contains(&runs),
                "{runs} runs on {threads} threads"
            );
        }
    }

    #[test]
    fn a_map_gives_each_item_in_order_however_it_is_chunked() {
        let want: Vec<usize> = (0..1000).map(|i| 3 * i + 1).collect();
        for chunk in [0, 1, 7, 64, 1000, 5000] {
            assert_eq!(
                super::map(1000, chunk, |i| 3 * i + 1),
                want,
                "chunks of {chunk}"
            );
        }
        assert!(super::map(0, 4, |i| i).is_empty());
    }

    #[test]
    fn the_calling_thread_always_takes_a_run() {
        let here = std::thread::current().id();
        assert_eq!(on(1, || std::thread::current().id()), [here]);
        assert!(on(4, || std::thread::current().id()).contains(&here));
    }

    #[cfg(panic = "unwind")]
    #[test]
    fn a_panic_in_the_work_is_the_callers_whichever_thread_it_struck() {
        use std::panic::{AssertUnwindSafe, catch_unwind};
        use std::sync::atomic::AtomicBool;

        let message = |caught: std::thread::Result<Vec<()>>| {
            caught
                .err()
                .and_then(|panic| panic.downcast_ref::<&str>().copied())
        };
        // On the calling thread.
        let caught = catch_unwind(|| on(1, || panic!("the work failed")));
        assert_eq!(message(caught), Some("the work failed"));
        // On every other thread, the calling thread's own run ending well: the caller sees a
        // spawned run's panic, where a thread could be spawned to make one.
        let here = std::thread::current().id();
        let struck = AtomicBool::new(false);
        let caught = catch_unwind(AssertUnwindSafe(|| {
            on(4, || {
                if std::thread::current().id() != here {
                    struck.store(true, Ordering::Relaxed);
                    panic!("the work failed");
                }
            })
        }));
        let want = struck.load(Ordering::Relaxed).then_some("the work failed");
        assert_eq!(message(caught), want);
    }
}
