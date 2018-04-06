# engineering_drawing_software

It contains one code file with the name generate2D_v2.0.cpp file.
Run it through command line - g++ generate2D_v2.0.cpp -std=c++11
It will ask for one number.
Number 0 represents conversion from 2D to 3D and number 1 represents converting 2D views to 3D object.
After entering either 1 or 0, we need to enter the viewing direction where it is represented by vector (a,b,c) and (0,0,0).
Keep the files in the same directory as your code.

For 2D to 3D conversion - 
Input - 4 files (for any 2 views)
TopViewVertices.txt, TopViewEdges.txt
SideViewVertices.txt, SideViewEdges.txt
Output - 2 files
edge.txt and vertex.txt

For 3D to 2D conversion - 
Input - 2 files 
vertices.txt, edges.txt
Output - 6 files (for 3 views)
TopViewVertices.txt, TopViewEdges.txt
SideViewVertices.txt, SideViewEdges.txt
FrontViewVertices.txt, FrontViewEdges.txt
Also, a 2D visualisation on screen.

Note - Coordinates should be all positive number within the range [0, 3] if 2D visualisation is required.

Format of any file - 
Vertex file - 
v1 0 0 0
v2 1 1 1
v1 and v2 are the names of the vertices while (0,0,0) and (1,1,1) are the coordinates of the corresponding vertices.

Edge File - 
v1 v2
v2 v5
v1 v6
v1 v2 denotes an edge between the vertices v1 and v2. Similarly, v2 v5 and v1 v6 denotes another 2 edges.

Other Notes -
refman.pdf contains the doxygen documentation for the code.
Class diagram.jpg/xml/pdf are the class diagrams.
Functional_Specification.jpg contains the functional specification and mathematical analysis is the initial maths used to design the software.
cad_v1_0.cpp contains the same code using openGL which may or may not be working and is under process.
Makefile couldn't be completed due to incomplete graphics part.
