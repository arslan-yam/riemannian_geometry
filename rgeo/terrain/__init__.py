"""Triangulated terrains as Riemannian manifolds: exact polyhedral lengths, exact labels
(MMP via pygeodesic), arbitrary-point queries."""
from .io import read_off, write_off, projected_areas, repair_height_field
from .mesh import TerrainMesh
from .exact import ExactTerrain
from .a2a import TerrainA2A

__all__ = ["read_off", "write_off", "projected_areas", "repair_height_field", "TerrainMesh", "ExactTerrain", "TerrainA2A"]
