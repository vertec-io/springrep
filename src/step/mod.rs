//! STEP AP214 import and export: [`write()`] turns an arena's solids into an exchange file whose
//! faces carry their names, and [`write_assembly`] several parts placed by their instances;
//! [`read`] turns a file's advanced B-reps, placed by its product structure, into checked kernel
//! solids, each with its product's name and its colour.

mod adjust;
mod geometry;
mod heal;
mod part21;
mod read;
mod topology;
mod write;

pub use heal::{Healed, MAX_HEAL};
pub use read::{
    Body, Crossed, Healing, Import, ReadOptions, Separated, Tolerant, TurnedBound, read, read_part,
    read_part_with, read_with,
};
pub(crate) use topology::{PartFace, Parts, parts_of, rebuilt, sewn, solid_of};
pub use write::{
    Instance, Part, WriteOptions, write, write_assembly, write_assembly_with, write_faces,
    write_with,
};

/// The faces of a STEP file, each its own with its own vertices and edges, sewn into one solid
/// where their edges lie within `tolerance` of each other, then healed onto their surfaces.
pub fn sew(text: &[u8], tolerance: f64) -> crate::Result<crate::topo::Brep> {
    let file = part21::parse(text)?;
    let mut geo = geometry::Geometry::new(&file);
    topology::sew(&mut geo, tolerance)
}

/// The `originating_system` this kernel writes; a file that names it carries the kernel's face names.
pub const WRITER: &str = "springrep";

#[cfg(test)]
mod tests;
