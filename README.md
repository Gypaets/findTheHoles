# README #

## About  
 findTheHoles is a 2D mesh reconstruction tool which automatically identifies holes in a point cloud.   

![Image](images/Titelbild.png?raw=true)

## Usage  
* Input:  
  + XY= Nx2 matrix with point coordinates.  
* Optional arguments:  
  + S = Critical area ratio (real and positive number). Knots with a value higher than the ratio of maximum to minimum neighbouring polygon area are identified as hole-border points. Default value is 3.  
  + M = Flag for manual editing. Default value is 0.  
      + 0: No manual editing.  
	  + 1: Manual editing of hole-border points. After the first calculation a figure with the point cloud is opened. Automatically recognized points are marked with an 'x'. The user can add/remove not/wrongly recognized border points selecting them with the brush and clicking the appropiate button.
	  ![Screenshot manual-editing mode](images/ManualMode.png?raw=true)
   + T = Triangulation (Nx3 matrix with polygon vertices). If not given the delaunay triangulation of X,Y is used.  
  
* Output:  
  * triangulation  
     + .Points: nx2 matrix with points' coordinates.  
     + .ConnectivityList: nx3 matrix with the triangulation polygons.   
  * holes: Cell array. Each cell element is an nx1 matrix with ordered hole border point indices.  
     
## Method  
 To identify the holes a triangulation of the points cloud is needed. If none is given, the delaunay triangulation is used. Points where the ratio between maximum and minimum neighbouring polygon area exceeds a critical value are identified as hole-border points. Polygons for which all vertices are identified as hole-border points are deleted. The adjacency matrix of the resulting mesh is clustered to get the holes.   
  
## Issues  
 The clustering algorithm isn't fully developed yet and has issues with 
 interconnected holes among others.  
  
## Examples  
 See file examples_findTheHoles.m  
  
## Author  
 https://github.com/Gypaets
