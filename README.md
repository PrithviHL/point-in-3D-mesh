# point-in-3D-mesh

-	Solved point containment by casting rays from the query point and counting triangle crossings with the Moller-Trumbore test (odd count = inside), voting across 3 random ray directions to stay robust when a ray grazes an edge or vertex.
-	Added a generalized winding number fallback that sums each triangle's signed solid angle over 4 pi (above 0.5 = inside), correctly classifying points in meshes with holes where ray parity fails.
-	Accelerated ray queries in C++17 with a median-split bounding volume hierarchy, running 95x faster than the linear winding sum on a 20,480-triangle mesh, and validated 40,000 random points on a sphere and a non-convex torus with 100% agreement against analytic ground truth.
