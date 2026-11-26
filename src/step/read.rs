//! Reading a STEP file: the product structure gives each solid's placements, and each solid is
//! built once and copied to every placement.

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

use super::{
    WRITER,
    geometry::{Geometry, Placing},
    heal::Healed,
    part21::{self, File, Param},
    topology,
};
use crate::{
    Error, Result, check,
    math::Point3,
    threads, tol,
    topo::{Brep, FaceId, FaceName, SolidId},
};

/// How a STEP file is read. The default makes each solid agree with itself and keeps what it can,
/// reporting every change in the [`Import`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct ReadOptions {
    /// Keep only solids whose geometry the file already gives agreeing with itself within
    /// [`tol::LINEAR`], and that do not pass through themselves. A solid the default read would
    /// heal, keep with local tolerances, or keep while it passes through itself is rejected
    /// instead, with the reason: so `healed` and `tolerant` stay empty, every solid kept is the
    /// file's own geometry, and `crossed` holds only places where a solid touches itself. Coincident
    /// vertices and edges are still joined, a shell written inside out is still turned, a bound
    /// flagged against the faces round it is still turned ([`TurnedBound`]), and units and
    /// placements still apply: none of those moves geometry.
    pub strict: bool,
}

impl ReadOptions {
    #[deprecated(since = "2.1.0", note = "use `ReadOptions::default().strict(..)`")]
    pub const fn new(strict: bool) -> ReadOptions {
        ReadOptions { strict }
    }

    /// These options with [`strict`](ReadOptions::strict) set to `on`.
    pub const fn strict(self, on: bool) -> ReadOptions {
        ReadOptions { strict: on, ..self }
    }
}

/// What a read produced: every solid that became a valid kernel solid, how far healing moved
/// each one's geometry, and why any other was rejected.
#[derive(Debug)]
#[non_exhaustive]
pub struct Import {
    pub brep: Brep,
    /// Rejected solids and placements, each as `#id: error`.
    pub rejected: Vec<String>,
    /// One entry per placed solid whose geometry healing moved, in placement order.
    pub healed: Vec<Healing>,
    /// One entry per placed solid kept as the file has it, with local tolerances, in placement
    /// order.
    pub tolerant: Vec<Tolerant>,
    /// Every place a placed solid passes through or touches itself, in placement order. Such a
    /// solid is kept as the file has it: an imported mesh stitched with folds microns deep is still
    /// the part, and other systems keep it too.
    pub crossed: Vec<Crossed>,
    /// One entry per placed solid written with voids outside its outer shell, each of which was
    /// read as a solid of its own, in placement order. A strict read rejects such a solid instead.
    pub separated: Vec<Separated>,
    /// One entry per face of a placed solid whose bounds the read takes in the sense the faces
    /// across their edges agree with rather than the sense the file's flags give them, in
    /// placement order. A strict read turns them too: turning a bound moves no geometry.
    pub turned_bounds: Vec<TurnedBound>,
    /// One entry per solid of `brep`, in its order: what the file calls it and how it draws it.
    pub bodies: Vec<Body>,
    /// The faces the file gives a colour of their own, apart from their body's.
    pub face_colours: Vec<(FaceId, [f64; 3])>,
}

/// A solid of an import as its file describes it.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Body {
    pub solid: SolidId,
    /// The file's id of the solid it came from; an assembly places one solid many times.
    pub step_id: u64,
    /// The name of the product whose shape holds the solid, where the file has products.
    pub name: Option<String>,
    /// The uses that place it, from the root assembly down: each use's name, or its id where it
    /// has none. Empty for a solid of a root product.
    pub instance: Vec<String>,
    /// The body's surface colour, red, green and blue in 0 to 1, where the file gives one.
    pub colour: Option<[f64; 3]>,
}

/// Two faces of one placed solid that share no edge but meet.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Crossed {
    /// The solid's position among `Import::brep`'s solids.
    pub index: usize,
    /// The solid's STEP instance.
    pub step_id: u64,
    /// The two faces, in `Import::brep`.
    pub faces: (FaceId, FaceId),
    /// Where they meet, placed as the solid is.
    pub at: Point3,
    /// How far one passes through the other, in mm: zero where they only touch.
    pub depth: f64,
}

/// A placed solid whose file lists, among its voids, shells that lie outside its outer shell. A
/// void is a hole in the material the outer shell bounds (ISO 10303-42, `brep_with_voids`), so such
/// a file writes a compound of bodies as one solid; the read keeps each of those shells as a solid
/// of its own, bounded by it as the file has it, however the file turns it.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Separated {
    /// The position among `Import::brep`'s solids of the solid the outer shell bounds, with any
    /// voids inside it; the solids read from the shells outside it follow it, in `shells`' order.
    pub index: usize,
    /// The solid's STEP instance.
    pub step_id: u64,
    /// The shells read as solids of their own, by the instances the solid lists them as.
    pub shells: Vec<u64>,
}

/// A face whose bounds the read turned against the file's flags. On a surface closed in both
/// directions, such as a torus or a sphere, either sense of a loop bounds a region of it, so a
/// `FACE_BOUND`'s orientation flag alone picks which one the face is (ISO 10303-42, `face_bound`).
/// Where the flag turns a loop the same way as the faces across every one of its shared edges,
/// so that the shell would not close, the read takes the sense those faces agree with.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct TurnedBound {
    /// The solid's position among `Import::brep`'s solids.
    pub index: usize,
    /// The solid's STEP instance.
    pub step_id: u64,
    /// The face, in `Import::brep`.
    pub face: FaceId,
    /// The face's STEP instance.
    pub face_step_id: u64,
    /// How many of the face's bounds were turned.
    pub bounds: usize,
}

/// A placed solid whose geometry could not be healed to meet within [`crate::tol::LINEAR`]
/// without moving it further than the file says it meets, so its edges and vertices carry the
/// tolerances it does meet within. Those are measured, and may be wider than the file states.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Tolerant {
    /// The solid's position among `Import::brep`'s solids. A void read as a solid of its own (see
    /// [`Separated`]) is reported apart, at its own position.
    pub index: usize,
    /// The solid's STEP instance.
    pub step_id: u64,
    /// How many of its vertices and edges carry a tolerance above the kernel's.
    pub vertices: usize,
    pub edges: usize,
    /// The largest tolerance of any of them, in mm.
    pub largest: f64,
    /// The uncertainty the file states for the solid, in mm, where it states one.
    pub stated: Option<f64>,
    /// How many of them carry a tolerance above `stated`: where the file misses its own word.
    pub past_stated: usize,
}

/// How far healing moved one placed solid's geometry to close it within the kernel tolerance.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Healing {
    /// The solid's position among `Import::brep`'s solids. A void read as a solid of its own (see
    /// [`Separated`]) is reported apart, at its own position.
    pub index: usize,
    /// The solid's STEP instance.
    pub step_id: u64,
    /// The largest distance any vertex, edge or surface moved, in mm.
    pub moved: f64,
    /// Where that movement was, placed as the solid is.
    pub at: Point3,
    /// What moved most, such as `vertex #12` or `the surface of face #34`.
    pub what: String,
    /// The largest distance any of its surfaces moved, which bounds how far its faces moved.
    pub surfaces: f64,
}

const SOLIDS: [&str; 4] = [
    "MANIFOLD_SOLID_BREP",
    "BREP_WITH_VOIDS",
    "FACETED_BREP",
    "TESSELLATED_SOLID",
];

/// The most solid placements one file may expand to; an assembly that reuses its parts
/// exponentially deep is refused rather than expanded.
const MAX_PLACEMENTS: usize = 1_000_000;

fn cyclic<T>(id: u64, what: &str) -> Result<T> {
    Err(Error::InvalidInput(format!("STEP #{id}: {what}")))
}

#[derive(Default)]
struct Placements {
    solids: Vec<Placement>,
    rejected: Vec<String>,
}

/// A solid to build, the representation context its units come from, and where it goes.
struct Placement {
    solid: u64,
    context: Option<u64>,
    at: Placing,
    /// The product definition whose shape holds the solid.
    product: Option<u64>,
    /// The assembly uses that place it, root first.
    uses: Vec<u64>,
}

impl Placement {
    /// What the solid is built from: the same solid in the same units at the same scale is built
    /// once.
    fn key(&self) -> Key {
        (self.solid, self.context, self.at.scale.to_bits())
    }
}

type Key = (u64, Option<u64>, u64);

fn params(f: &File, id: u64, kind: &str) -> Option<Vec<Param>> {
    f.get(id)?.part(kind).map(<[Param]>::to_vec)
}

fn refs(p: Option<&Param>) -> Vec<u64> {
    p.and_then(Param::as_list)
        .map(|l| l.iter().filter_map(Param::as_ref).collect())
        .unwrap_or_default()
}

/// The representation an instance is (any entity with a `REPRESENTATION` part or subtype that
/// carries `name, items, context`).
fn representation(f: &File, id: u64) -> Option<(Vec<u64>, Option<u64>)> {
    let e = f.get(id)?;
    let (_, p) = e.parts.iter().find(|(k, p)| {
        k.ends_with("REPRESENTATION") && !k.contains("RELATIONSHIP") && p.len() >= 3
    })?;
    Some((refs(p.get(1)), p.get(2).and_then(Param::as_ref)))
}

struct Assembly<'a> {
    file: &'a File,
    /// Representations joined without a transformation (a part's shape and its B-rep).
    links: HashMap<u64, Vec<u64>>,
    /// Each product definition's shape representations.
    reps: HashMap<u64, Vec<u64>>,
    /// Each product definition's children: the child, the relation that places it, and the use.
    children: HashMap<u64, Vec<(u64, u64, u64)>>,
    related: HashSet<u64>,
}

impl<'a> Assembly<'a> {
    fn new(file: &'a File) -> Assembly<'a> {
        let mut a = Assembly {
            file,
            links: HashMap::new(),
            reps: HashMap::new(),
            children: HashMap::new(),
            related: HashSet::new(),
        };
        let mut nauo_place: HashMap<u64, u64> = HashMap::new();
        for e in &file.entities {
            if let Some(p) = e.part("SHAPE_DEFINITION_REPRESENTATION") {
                let (Some(pds), Some(rep)) = (
                    p.first().and_then(Param::as_ref),
                    p.get(1).and_then(Param::as_ref),
                ) else {
                    continue;
                };
                if let Some(pd) = params(file, pds, "PRODUCT_DEFINITION_SHAPE")
                    .and_then(|q| q.get(2).and_then(Param::as_ref))
                {
                    a.reps.entry(pd).or_default().push(rep);
                }
            } else if let Some(p) = e.part("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION") {
                let (Some(rel), Some(pds)) = (
                    p.first().and_then(Param::as_ref),
                    p.get(1).and_then(Param::as_ref),
                ) else {
                    continue;
                };
                if let Some(nauo) = params(file, pds, "PRODUCT_DEFINITION_SHAPE")
                    .and_then(|q| q.get(2).and_then(Param::as_ref))
                {
                    nauo_place.insert(nauo, rel);
                }
            } else if !e.has("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
                && let Some(p) = e.part("REPRESENTATION_RELATIONSHIP").or_else(|| {
                    e.part("SHAPE_REPRESENTATION_RELATIONSHIP")
                        .filter(|p| p.len() >= 4)
                })
                && let (Some(r1), Some(r2)) = (
                    p.get(2).and_then(Param::as_ref),
                    p.get(3).and_then(Param::as_ref),
                )
            {
                a.links.entry(r1).or_default().push(r2);
                a.links.entry(r2).or_default().push(r1);
            }
        }
        for (i, e) in file.entities.iter().enumerate() {
            if let Some(p) = e.part("NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
                let id = file.ids[i];
                let (Some(parent), Some(child)) = (
                    p.get(3).and_then(Param::as_ref),
                    p.get(4).and_then(Param::as_ref),
                ) else {
                    continue;
                };
                if let Some(&rel) = nauo_place.get(&id) {
                    a.children.entry(parent).or_default().push((child, rel, id));
                    a.related.insert(child);
                }
            }
        }
        a
    }

    /// Every representation joined to `rep` without a transformation, `rep` included.
    fn linked(&self, rep: u64) -> Vec<u64> {
        let mut seen = vec![rep];
        let mut i = 0;
        while i < seen.len() {
            for &r in self.links.get(&seen[i]).into_iter().flatten() {
                if !seen.contains(&r) {
                    seen.push(r);
                }
            }
            i += 1;
        }
        seen
    }

    /// The solids in representation `rep` and those joined to it, at `at`. `path` holds the
    /// representations mapped into this one, so a map that contains itself is an error.
    fn rep_solids(
        &self,
        geo: &mut Geometry,
        rep: u64,
        at: &Placing,
        path: &mut Vec<u64>,
        out: &mut Placements,
    ) -> Result<()> {
        if path.contains(&rep) {
            return cyclic(rep, "a mapped representation contains itself");
        }
        path.push(rep);
        let r = self.rep_solids_in(geo, rep, at, path, out);
        path.pop();
        r
    }

    fn rep_solids_in(
        &self,
        geo: &mut Geometry,
        rep: u64,
        at: &Placing,
        path: &mut Vec<u64>,
        out: &mut Placements,
    ) -> Result<()> {
        for r in self.linked(rep) {
            let Some((items, context)) = representation(self.file, r) else {
                continue;
            };
            for it in items {
                let Some(e) = self.file.get(it) else { continue };
                if SOLIDS.contains(&e.kind()) {
                    if out.solids.len() >= MAX_PLACEMENTS {
                        return cyclic(it, "the assembly places too many solids");
                    }
                    out.solids.push(Placement {
                        solid: it,
                        context,
                        at: *at,
                        product: None,
                        uses: Vec::new(),
                    });
                } else if let Some(p) = e.part("MAPPED_ITEM") {
                    match self.mapped(geo, it, p, context) {
                        Ok((src, place)) => {
                            self.rep_solids(geo, src, &place.then(at), path, out)?;
                        }
                        Err(e) => out.rejected.push(format!("#{it}: {e}")),
                    }
                }
            }
        }
        Ok(())
    }

    /// A `MAPPED_ITEM`'s source representation and its placement in the target.
    fn mapped(
        &self,
        geo: &mut Geometry,
        item: u64,
        p: &[Param],
        context: Option<u64>,
    ) -> Result<(u64, Placing)> {
        let map_id = geo.reference(item, p.get(1))?;
        let map = self
            .file
            .get(map_id)
            .and_then(|e| e.part("REPRESENTATION_MAP"))
            .ok_or_else(|| {
                Error::InvalidInput(format!("STEP #{map_id}: expected a representation map"))
            })?;
        let src = geo.reference(map_id, map.get(1))?;
        geo.use_context(context)?;
        let origin = geo.transform_of(geo.reference(map_id, map.first())?)?;
        let target = geo.transform_of(geo.reference(item, p.get(2))?)?;
        Ok((src, origin.inverse().then(&target)))
    }

    /// The motion that places `child`'s shape in its parent's, from the placing relation.
    fn placement(&self, geo: &mut Geometry, child: u64, rel: u64) -> Result<Placing> {
        let e = self
            .file
            .get(rel)
            .ok_or_else(|| Error::InvalidInput(format!("STEP #{rel}: missing instance")))?;
        let rr = e.part("REPRESENTATION_RELATIONSHIP").ok_or_else(|| {
            Error::InvalidInput(format!(
                "STEP #{rel}: expected a representation relationship"
            ))
        })?;
        let r1 = geo.reference(rel, rr.get(2))?;
        let r2 = geo.reference(rel, rr.get(3))?;
        let transformation = e
            .part("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
            .ok_or_else(|| {
                Error::InvalidInput(format!("STEP #{rel}: expected a placing transformation"))
            })?;
        let op = geo.reference(rel, transformation.first())?;
        let child_reps: Vec<u64> = self
            .reps
            .get(&child)
            .into_iter()
            .flatten()
            .flat_map(|&r| self.linked(r))
            .collect();
        let child_is_r1 = child_reps.contains(&r1) || !child_reps.contains(&r2);
        let context = |r: u64| representation(self.file, r).and_then(|x| x.1);
        if let Some(it) = self
            .file
            .get(op)
            .and_then(|e| e.part("ITEM_DEFINED_TRANSFORMATION"))
        {
            geo.use_context(context(r1))?;
            let a = geo.transform_of(geo.reference(op, it.get(2))?)?;
            geo.use_context(context(r2))?;
            let b = geo.transform_of(geo.reference(op, it.get(3))?)?;
            let (from, to) = if child_is_r1 { (a, b) } else { (b, a) };
            Ok(from.inverse().then(&to))
        } else {
            geo.use_context(context(if child_is_r1 { r2 } else { r1 }))?;
            let t = geo.transform_of(op)?;
            Ok(if child_is_r1 { t } else { t.inverse() })
        }
    }

    /// Places `pd` and its children at `at`. `path` holds the products above it, so an assembly
    /// that contains itself is an error. With no placement, the walk only checks product links.
    fn walk(
        &self,
        geo: &mut Geometry,
        pd: u64,
        at: Option<&Placing>,
        path: &mut Vec<u64>,
        seen: &mut HashSet<u64>,
        out: &mut Placements,
    ) -> Result<()> {
        if path.contains(&pd) {
            return cyclic(pd, "the assembly contains itself");
        }
        if !seen.insert(pd) && at.is_none() {
            return Ok(());
        }
        path.push(pd);
        if let Some(at) = at {
            let own = out.solids.len();
            for &rep in self.reps.get(&pd).into_iter().flatten() {
                self.rep_solids(geo, rep, at, &mut Vec::new(), out)?;
            }
            for p in &mut out.solids[own..] {
                p.product = Some(pd);
            }
        }
        for &(child, rel, used) in self.children.get(&pd).into_iter().flatten() {
            let child_at = at.and_then(|at| match self.placement(geo, child, rel) {
                Ok(t) => Some(t.then(at)),
                Err(e) => {
                    out.rejected.push(format!("#{used}: {e}"));
                    None
                }
            });
            let below = out.solids.len();
            self.walk(geo, child, child_at.as_ref(), path, seen, out)?;
            for p in &mut out.solids[below..] {
                p.uses.insert(0, used);
            }
        }
        path.pop();
        Ok(())
    }

    fn placements(&self, geo: &mut Geometry) -> Result<Placements> {
        let mut out = Placements::default();
        let mut roots: Vec<u64> = self
            .reps
            .keys()
            .chain(self.children.keys())
            .copied()
            .filter(|pd| !self.related.contains(pd))
            .collect();
        roots.sort_unstable();
        roots.dedup();
        let mut seen = HashSet::new();
        for pd in roots {
            self.walk(
                geo,
                pd,
                Some(&Placing::IDENTITY),
                &mut Vec::new(),
                &mut seen,
                &mut out,
            )?;
        }
        // A product no root reaches sits on a cycle of assembly uses.
        let mut unreached: Vec<u64> = self
            .children
            .keys()
            .copied()
            .filter(|pd| !seen.contains(pd))
            .collect();
        unreached.sort_unstable();
        match unreached.first() {
            Some(&pd) => cyclic(pd, "the assembly contains itself"),
            None => Ok(out),
        }
    }
}

/// Reads a STEP file as one part: the union of its solids, as the pieces of a solid that touches
/// itself are written one by one and join again. A file with a solid this kernel rejects is an
/// error, since the part would lack it.
pub fn read_part(text: &[u8], feature: &str) -> Result<Brep> {
    read_part_with(text, feature, &ReadOptions::default())
}

/// [`read_part`], read as `options` say.
pub fn read_part_with(text: &[u8], feature: &str, options: &ReadOptions) -> Result<Brep> {
    let import = read_with(text, feature, options)?;
    if !import.rejected.is_empty() {
        return Err(Error::InvalidResult(import.rejected));
    }
    let b = &import.brep;
    let pieces: Vec<Brep> = b
        .solid_ids()
        .map(|s| crate::boolean::solids_of(b, &[s]))
        .collect();
    match &pieces[..] {
        [one] => Ok(one.clone()),
        _ => crate::boolean::union_all(&pieces.iter().collect::<Vec<_>>()),
    }
}

/// Reads a STEP file into kernel solids. Faces carry the names a file this kernel wrote
/// gave them; otherwise each is `<feature>.face_<i>`, numbered in the order the solids are placed.
pub fn read(text: &[u8], feature: &str) -> Result<Import> {
    read_with(text, feature, &ReadOptions::default())
}

/// [`read`], read as `options` say.
pub fn read_with(text: &[u8], feature: &str, options: &ReadOptions) -> Result<Import> {
    let file = part21::parse(text)?;
    let ours = file.originating_system() == Some(WRITER);
    let mut geo = Geometry::new(&file);
    let Placements {
        solids: mut placements,
        mut rejected,
    } = Assembly::new(&file).placements(&mut geo)?;
    if placements.is_empty() && rejected.is_empty() {
        // Without a product structure, a solid takes its units and uncertainty from a
        // representation that holds it, if one does.
        let mut contexts: HashMap<u64, u64> = HashMap::new();
        for &id in &file.ids {
            if let Some((items, Some(ctx))) = representation(&file, id) {
                for it in items {
                    contexts.entry(it).or_insert(ctx);
                }
            }
        }
        for (i, e) in file.entities.iter().enumerate() {
            if SOLIDS.contains(&e.kind()) {
                placements.push(Placement {
                    solid: file.ids[i],
                    context: contexts.get(&file.ids[i]).copied(),
                    at: Placing::IDENTITY,
                    product: None,
                    uses: Vec::new(),
                });
            }
        }
        if placements.is_empty() {
            return Err(Error::InvalidInput("STEP: the file holds no solids".into()));
        }
    }
    let mut keys: Vec<Key> = placements.iter().map(Placement::key).collect();
    keys.sort_unstable();
    keys.dedup();
    let built = build_all(&file, &keys, options.strict);
    let mut out = Brep::new();
    let mut healed = Vec::new();
    let mut tolerant = Vec::new();
    let mut crossed = Vec::new();
    let mut separated = Vec::new();
    let mut turned_bounds = Vec::new();
    let mut crossings: HashMap<Key, Vec<check::SelfCrossing>> = HashMap::new();
    let mut counter = 0usize;
    let mut index = 0usize;
    let colours = colours(&file);
    let mut bodies = Vec::new();
    let mut face_colours = Vec::new();
    for pl in placements {
        match &built[&pl.key()] {
            Err(e) => rejected.push(format!("#{}: {e}", pl.solid)),
            Ok((b, names, heal, stated, apart, turned)) => {
                let found = crossings
                    .entry(pl.key())
                    .or_insert_with(|| check::self_crossings(b))
                    .clone();
                if options.strict {
                    let through = found.iter().filter(|x| x.depth > tol::LINEAR);
                    if let Some(deepest) = through.clone().map(|x| x.depth).reduce(f64::max) {
                        rejected.push(format!(
                            "#{}: the solid passes through itself in {} places, {deepest:.1e} mm at \
                             the deepest, and a strict read keeps only solids that do not",
                            pl.solid,
                            through.count()
                        ));
                        continue;
                    }
                }
                let mut copy = b.clone();
                if let Err(e) = pl.at.place_scaled_body(&mut copy) {
                    rejected.push(format!("#{}: {e}", pl.solid));
                    continue;
                }
                for (k, h) in heal {
                    healed.push(Healing {
                        index: index + k,
                        step_id: pl.solid,
                        moved: h.moved,
                        at: pl.at.place_scaled_point(h.at),
                        what: h.what.clone(),
                        surfaces: h.surfaces,
                    });
                }
                if b.has_tolerances() {
                    let (vs, es) = b.tolerant();
                    // A solid read with voids outside it as solids of their own is several, each
                    // reported apart.
                    let (vertex_in, edge_in) = solid_of_each(b);
                    for (k, _) in b.solid_ids().enumerate() {
                        let vs: Vec<f64> = vs
                            .iter()
                            .filter(|x| vertex_in[x.0.index()] == k)
                            .map(|x| x.1)
                            .collect();
                        let es: Vec<f64> = es
                            .iter()
                            .filter(|x| edge_in[x.0.index()] == k)
                            .map(|x| x.1)
                            .collect();
                        if vs.is_empty() && es.is_empty() {
                            continue;
                        }
                        let past = |t: &&f64| stated.is_some_and(|s| **t > s);
                        tolerant.push(Tolerant {
                            index: index + k,
                            step_id: pl.solid,
                            vertices: vs.len(),
                            edges: es.len(),
                            largest: vs.iter().copied().fold(0.0, f64::max),
                            stated: *stated,
                            past_stated: vs.iter().chain(&es).filter(past).count(),
                        });
                    }
                }
                if !apart.is_empty() {
                    separated.push(Separated {
                        index,
                        step_id: pl.solid,
                        shells: apart.clone(),
                    });
                }
                let first = out.face_count();
                let solid_index = index;
                index += b.solid_ids().count();
                let mut own_colours = Vec::new();
                let faces = copy.face_ids().collect::<Vec<_>>();
                for (k, (f, (step, name))) in faces.into_iter().zip(names).enumerate() {
                    if let Some(&c) = colours.get(step) {
                        own_colours.push((k, c));
                    }
                    let parsed = if ours { parse_names(name) } else { Vec::new() };
                    copy.face_mut(f).names = if parsed.is_empty() {
                        [FaceName::new(feature, format!("face_{counter}"))].into()
                    } else {
                        parsed.into_iter().collect()
                    };
                    counter += 1;
                }
                let solids = out.append(&copy);
                let ids: Vec<FaceId> = out.face_ids().skip(first).collect();
                face_colours.extend(own_colours.into_iter().map(|(k, c)| (ids[k], c)));
                bodies.extend(solids.into_iter().map(|solid| Body {
                    solid,
                    step_id: pl.solid,
                    name: pl.product.and_then(|pd| product_name(&file, pd)),
                    instance: pl.uses.iter().map(|&u| use_name(&file, u)).collect(),
                    colour: colours.get(&pl.solid).copied(),
                }));
                for &(k, bounds) in turned {
                    // A solid read with voids outside it as solids of their own is several.
                    let within = b.face_ids().nth(k).and_then(|f| b.face(f).shell);
                    let within = within.and_then(|sh| b.shell(sh).solid);
                    turned_bounds.push(TurnedBound {
                        index: solid_index + within.map_or(0, |s| s.index()),
                        step_id: pl.solid,
                        face: ids[k],
                        face_step_id: names[k].0,
                        bounds,
                    });
                }
                for x in found {
                    // A solid read with voids outside it as solids of their own is several.
                    let within = b.face(x.a).shell.and_then(|sh| b.shell(sh).solid);
                    crossed.push(Crossed {
                        index: solid_index + within.map_or(0, |s| s.index()),
                        step_id: pl.solid,
                        faces: (ids[x.a.index()], ids[x.b.index()]),
                        at: pl.at.place_scaled_point(x.at),
                        depth: x.depth,
                    });
                }
            }
        }
    }
    Ok(Import {
        brep: out,
        rejected,
        healed,
        tolerant,
        crossed,
        separated,
        turned_bounds,
        bodies,
        face_colours,
    })
}

fn text(p: Option<&Param>) -> Option<String> {
    p.and_then(Param::as_str)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
}

/// The name of the product a product definition defines, or its id where it has no name.
fn product_name(file: &File, pd: u64) -> Option<String> {
    let formation = params(file, pd, "PRODUCT_DEFINITION")?.get(2)?.as_ref()?;
    let formed = file.get(formation)?;
    let product = formed
        .parts
        .iter()
        .find(|(k, _)| k.starts_with("PRODUCT_DEFINITION_FORMATION"))?
        .1
        .get(2)?
        .as_ref()?;
    let p = params(file, product, "PRODUCT")?;
    text(p.get(1)).or_else(|| text(p.first()))
}

/// What an assembly use is called: its name, else its id, else its number in the file.
fn use_name(file: &File, used: u64) -> String {
    params(file, used, "NEXT_ASSEMBLY_USAGE_OCCURRENCE")
        .and_then(|p| text(p.get(1)).or_else(|| text(p.first())))
        .unwrap_or_else(|| format!("#{used}"))
}

fn all_refs(p: &Param, out: &mut Vec<u64>) {
    match p {
        Param::Ref(r) => out.push(*r),
        Param::List(l) => l.iter().for_each(|x| all_refs(x, out)),
        Param::Typed(_, inner) => all_refs(inner, out),
        _ => {}
    }
}

/// The first colour reached from entity `id` within `depth` references.
fn colour_under(file: &File, id: u64, depth: usize) -> Option<[f64; 3]> {
    let e = file.get(id)?;
    if let Some(p) = e.part("COLOUR_RGB") {
        let at = |i: usize| p.get(i).and_then(Param::as_f64);
        return Some([at(1)?, at(2)?, at(3)?]);
    }
    if let Some(p) = e
        .part("DRAUGHTING_PRE_DEFINED_COLOUR")
        .or_else(|| e.part("PRE_DEFINED_COLOUR"))
    {
        return Some(match p.first().and_then(Param::as_str)? {
            "red" => [1.0, 0.0, 0.0],
            "green" => [0.0, 1.0, 0.0],
            "blue" => [0.0, 0.0, 1.0],
            "yellow" => [1.0, 1.0, 0.0],
            "magenta" => [1.0, 0.0, 1.0],
            "cyan" => [0.0, 1.0, 1.0],
            "black" => [0.0, 0.0, 0.0],
            "white" => [1.0, 1.0, 1.0],
            _ => return None,
        });
    }
    if depth == 0 {
        return None;
    }
    let mut below = Vec::new();
    for (_, ps) in &e.parts {
        ps.iter().for_each(|p| all_refs(p, &mut below));
    }
    below
        .into_iter()
        .find_map(|r| colour_under(file, r, depth - 1))
}

/// The surface colour the file's presentation styles give each solid or face they name, by its
/// id. A curve's or a point's colour is not a body's.
fn colours(file: &File) -> HashMap<u64, [f64; 3]> {
    // Usage, side style, fill area, fill style, fill colour, colour.
    const DEPTH: usize = 6;
    let mut out = HashMap::new();
    for e in &file.entities {
        let Some(p) = e
            .part("STYLED_ITEM")
            .or_else(|| e.part("OVER_RIDING_STYLED_ITEM"))
        else {
            continue;
        };
        let Some(item) = p.get(2).and_then(Param::as_ref) else {
            continue;
        };
        let mut styles = Vec::new();
        for assignment in refs(p.get(1)) {
            if let Some(a) = file.get(assignment) {
                for (_, ps) in &a.parts {
                    ps.iter().for_each(|x| all_refs(x, &mut styles));
                }
            }
        }
        let colour = styles
            .into_iter()
            .filter(|&s| file.get(s).is_some_and(|x| x.has("SURFACE_STYLE_USAGE")))
            .find_map(|s| colour_under(file, s, DEPTH));
        if let Some(c) = colour {
            out.entry(item).or_insert(c);
        }
    }
    out
}

/// A solid, its faces' ids and names, how far healing moved each of its solids it moved, by the
/// solid's position, the uncertainty its file states, the void shells outside its outer shell
/// read as solids of their own, and the faces whose bounds the read turned, as
/// [`topology::ReadSolid`].
type Built = std::result::Result<
    (
        Brep,
        Vec<(u64, String)>,
        Vec<(usize, Healed)>,
        Option<f64>,
        Vec<u64>,
        Vec<(usize, usize)>,
    ),
    String,
>;

/// The position among `b`'s solids of the solid each vertex and each edge of `b` bounds, by id;
/// `usize::MAX` for one that bounds none.
fn solid_of_each(b: &Brep) -> (Vec<usize>, Vec<usize>) {
    let mut vertices = vec![usize::MAX; b.vertex_ids().count()];
    let mut edges = vec![usize::MAX; b.edge_ids().count()];
    for (k, s) in b.solid_ids().enumerate() {
        for &sh in &b.solid(s).shells {
            for &f in &b.shell(sh).faces {
                for &lp in &b.face(f).loops {
                    for &c in &b.lp(lp).coedges {
                        let e = b.coedge(c).edge;
                        edges[e.index()] = k;
                        let edge = b.edge(e);
                        vertices[edge.start.index()] = k;
                        vertices[edge.end.index()] = k;
                    }
                }
            }
        }
    }
    (vertices, edges)
}

/// A rough measure of how long a solid takes to build: its faces, and the control points of its
/// spline surfaces, which healing fits.
fn weight(file: &File, solid: u64) -> usize {
    let refs_of = |id: u64, at: usize| -> Vec<u64> {
        let Some(e) = file.get(id) else {
            return Vec::new();
        };
        let Some(p) = e.parts.first().map(|p| &p.1) else {
            return Vec::new();
        };
        match p.get(at) {
            Some(Param::List(l)) => l.iter().filter_map(Param::as_ref).collect(),
            Some(x) => x.as_ref().into_iter().collect(),
            None => Vec::new(),
        }
    };
    let mut shells = refs_of(solid, 1);
    shells.extend(refs_of(solid, 2));
    let mut total = 0;
    for sh in shells {
        let sh = match file.get(sh).and_then(|e| e.part("ORIENTED_CLOSED_SHELL")) {
            Some(p) => p.get(2).and_then(Param::as_ref).unwrap_or(sh),
            None => sh,
        };
        for f in refs_of(sh, 1) {
            total += 1;
            let Some(face) = file.get(f) else { continue };
            let Some(s) = face.parts[0].1.get(2).and_then(Param::as_ref) else {
                continue;
            };
            let Some(e) = file.get(s) else { continue };
            // A complex instance's B_SPLINE_SURFACE part has no name before its degrees.
            let simple = super::geometry::SPLINE_SURFACES
                .iter()
                .find_map(|k| e.part(k));
            let net = match (e.part("B_SPLINE_SURFACE"), simple) {
                (Some(p), _) => p.get(2),
                (None, Some(p)) => p.get(3),
                _ => None,
            };
            total += net.and_then(Param::as_list).map_or(0, |rows| {
                rows.iter().filter_map(Param::as_list).map(<[_]>::len).sum()
            });
        }
    }
    total
}

/// What `build` makes with `cache`, or a rejection where it panics, so that a solid whose build
/// panics is rejected alone. A panic may leave `cache` half-written: an instance it was converting
/// stays marked as open, and the next solid to use that instance would be refused as referring to
/// itself. So `cache` is then replaced by a `fresh` one.
///
/// This holds only where panics unwind. Where they abort, as on wasm32-wasip1, a solid whose
/// build panics ends the whole read.
fn guarded<C, T>(
    cache: &mut C,
    fresh: impl FnOnce() -> C,
    build: impl FnOnce(&mut C) -> std::result::Result<T, String>,
) -> std::result::Result<T, String> {
    let built = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| build(&mut *cache)));
    built.unwrap_or_else(|_| {
        *cache = fresh();
        Err("its build panicked".to_string())
    })
}

/// Builds each distinct solid once. Solids share nothing but the file, so they are built on as
/// many threads as the machine offers, each with its own geometry cache, heaviest first so that
/// no long build starts last. Each result is filed under its own key, whatever finished first.
/// A solid whose build panics is rejected and the rest are still built, where panics unwind
/// ([`guarded`]).
fn build_all(file: &File, keys: &[Key], strict: bool) -> HashMap<Key, Built> {
    let mut order: Vec<(usize, usize)> = keys
        .iter()
        .enumerate()
        .map(|(k, key)| (weight(file, key.0), k))
        .collect();
    order.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let next = AtomicUsize::new(0);
    let results: Vec<Vec<(usize, Built)>> = threads::spread(keys.len(), || {
        let mut geo = Geometry::new(file);
        let mut mine = Vec::new();
        while let Some(&(_, k)) = order.get(next.fetch_add(1, Ordering::Relaxed)) {
            let (solid, context, scale) = keys[k];
            let r = guarded(
                &mut geo,
                || Geometry::new(file),
                |geo| {
                    geo.use_scaled(context, f64::from_bits(scale))
                        .and_then(|_| topology::solid(geo, solid, strict))
                        .map(|(b, names, healed, apart, turned)| {
                            (b, names, healed, geo.uncertainty(), apart, turned)
                        })
                        .map_err(|e| e.to_string())
                },
            );
            mine.push((k, r));
        }
        mine
    });
    // Every key was taken by one run or another, and each one taken was filed.
    let built = results.into_iter().flatten();
    built.map(|(k, r)| (keys[k], r)).collect()
}

/// The names a face this kernel wrote carries: `feature.face` entries joined by `;`.
fn parse_names(s: &str) -> Vec<FaceName> {
    s.split(';')
        .filter_map(|n| {
            let (f, face) = n.split_once('.')?;
            Some(FaceName::new(f, face))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::guarded;

    #[test]
    fn a_build_keeps_its_result_and_its_cache() {
        // The cache is a list a build adds to, and its result the list's length.
        let add = |c: &mut Vec<usize>| {
            c.push(c.len());
            Ok::<usize, String>(c.len())
        };
        let mut cache = Vec::new();
        assert_eq!(guarded(&mut cache, Vec::new, add), Ok(1));
        assert_eq!(guarded(&mut cache, Vec::new, add), Ok(2));
        assert_eq!(cache, [0, 1]);
        // A build's own error is its result, and keeps the cache too.
        let refuse = |_: &mut Vec<usize>| Err::<usize, String>("no".into());
        assert_eq!(guarded(&mut cache, Vec::new, refuse), Err("no".into()));
        assert_eq!(cache, [0, 1]);
    }

    #[cfg(panic = "unwind")]
    #[test]
    fn a_build_that_panics_is_rejected_alone_and_its_cache_replaced() {
        // Five solids built in turn with one cache, as one worker builds them. Each adds itself to
        // the cache and returns what the cache then holds; the third panics once it has.
        let mut cache: Vec<usize> = Vec::new();
        let built: Vec<Result<Vec<usize>, String>> = (0..5)
            .map(|k| {
                guarded(&mut cache, Vec::new, |c| {
                    c.push(k);
                    if k == 2 {
                        panic!("solid {k} cannot be built");
                    }
                    Ok(c.clone())
                })
            })
            .collect();
        // The first two share a cache; the third is rejected; the last two share a new cache,
        // without the third's mark in it.
        let want = [
            Ok(vec![0]),
            Ok(vec![0, 1]),
            Err("its build panicked".to_string()),
            Ok(vec![3]),
            Ok(vec![3, 4]),
        ];
        assert_eq!(built, want);
    }
}
