# Masters Thesis: Improved Normal Estimation from Cross-Section Drawings

Codebase for masters thesis, *provided as-is*, includes all external dependencies. Built using Qt; UI created with ease of extension in mind so please feel free to expand on the codebase.

Thesis can be found [here](https://dx.doi.org/10.14288/1.0395345).

*This thesis is original, unpublished, work by the author S. Burla, and builds on the original CrossShade implementation by C. Shao.*

## ABSTRACT
Concept sketches are common in the early stages of production design, as they require little time to produce and are effective at communicating both shape and volume information via carefully placed strokes. However, they are limited in their ability to convey material properties, and they do not help to envision a fully shaded product. Affording artists and designers the ability for fast sketching of shapes with rich visual shading can save a significant amount of time during the design refinement process, where shading is often done manually across multiple variations of lighting and material properties per iteration of ideas. To meet this need, CrossShade is an existing tool designed by Shao et al. for automatically shading line drawings. By leveraging specialized curves, in particular cross-section curves and the intersections between them, surface information of the curve-described geometry can be solved for via optimization of cross-section curve plane normals. Cross-section curves enforce a set of geometric constraints that aid viewers in lifting and interpreting concept sketches into 3D space, and so readily supplement regularity cues and perceptual constraints in a 3D reconstruction framework. We adopt the same regularity cues employed by CrossShade to estimate normals along cross-section curves in a new optimization routine, individually target curve networks constituting the whole of a line drawing sketch, and include a simple solution selection heuristic yielding more efficient, and robust results than previously reported. Accompanying our work is an updated sketch interface to make the drawing process more enjoyable for artists and designers. 

## KEYWORDS
interactive tool, sketches, curve-network, normal reconstruction, 2D to 3D
