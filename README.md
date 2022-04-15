# 3D-Modeling

## Base features
- Compute Normal Vectors
- Compute Angle Statistics
- Write a `.obj` file
- Implementation for `enext()`, `sym()`, `org()`, `dest()`, and `fnext()`
- Compute Number of Components
- Implementation for `orientTriangles()`
- Compute Vertex Normal Vectors for Smooth Shading
- Visualise boundary edges
![teapot](https://github.com/AnuKritiW/3D-Modeling/blob/main/images/teapot.gif)

## Highlights
- Barycentric Subdivision
    - Create a new vertex in the middle of every edge
    - Create a new vertex in the centroid of every triangle
    - Split each triangle into 6 triangles using the newly added vertices

![smallcase](https://github.com/AnuKritiW/3D-Modeling/blob/main/images/smallobj_subdiv.gif)

![subdiv_summary](https://github.com/AnuKritiW/3D-Modeling/blob/main/images/subdiv_summary.png)

- Mesh Relaxation
    - Performed Laplacian smoothing 
    - Each vertex is moved to the average of all the neighbouring vertices
![pumpkin](https://github.com/AnuKritiW/3D-Modeling/blob/main/images/pumpkin_meshrelaxation.gif)

![relaxation_summary](https://github.com/AnuKritiW/3D-Modeling/blob/main/images/relaxation_summary.png)


