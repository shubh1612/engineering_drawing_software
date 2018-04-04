#include<bits/stdc++.h>
using namespace std;

#define least 0.0001

/** This class defines a vertex of a 2D 
  * object which contains two coordinates 
  * and the vector of edges from
  * that vertex. */

class Vertex2D
{
	public :
	double x,y;

	/// Name of vertex
	string name;

	/// vector of names of vertices which have edge to particular vertex
	vector<Vertex2D> vertex2DEdgeVector;

	/// Method to add neighbouring vertices to the vector
	void vertex2DAddNeighbour(Vertex2D v){ vertex2DEdgeVector.push_back(v);} 
};


/** This class defines an edge of a 2D  
  * object which contains two 2D vertices  
  * of the edge. */

class Edge2D
{
	public:
	Vertex2D v1, v2;
};


/** This class defines a model of a 2D
  * object. */

class Model2D
{
	public:

	/// vector of vertices in the model
	vector<Vertex2D> vertexVector1, vertexVector2, vertexVector3;

	/// vector of edges in the model
	vector<Edge2D> edgesVector1, edgesVector2, edgesVector3;

	int edge1Length(){
		return edgesVector1.size();
	}

	int vertex1Length(){
		return vertexVector1.size();
	}

	int edge2Length(){
		return edgesVector2.size();
	}

	int vertex2Length(){
		return vertexVector2.size();
	}

	int edge3Length(){
		return edgesVector3.size();
	}

	int vertex3Length(){
		return vertexVector3.size();
	}
};


/** This class defines a vertex of a 3D 
  * object which contains three coordinates 
  * and the vector of edges from
  * that vertex. */

class Vertex3D
{
	public :

	/// 3 coordinates of vertex
	double x,y,z;

	/// Name of vertex
	string name;

	/// vector of edges 
	vector<Vertex3D> vertex3DEdgeVector;

	/// Method to add neighbouring vertices to vertex
	void vertex3DAddNeighbour(Vertex3D v){ vertex3DEdgeVector.push_back(v);} 
};


/** This class defines an edge of a 3D  
  * object which contains two 3D vertices  
  * of the edge. */

class Edge3D
{
	public:
	Vertex3D v1, v2;
};


/** This class defines a model of a 3D
  * object. */

class Model3D
{
	public:

	/// vector of vertices in the model	
	vector<Vertex3D> vertexVector;

	/// vector of edges in the model
	vector<Edge3D> edgesVector;

	/// Method to add edge to the model
	int addEdge(Vertex3D v1, Vertex3D v2, int length){

		// Check if the corresponding edge already exists inside the edge list or not
		for (int i=0;i<length;i++){
			if((edgesVector[i].v1.name == v1.name && edgesVector[i].v2.name == v2.name) ||
			  (edgesVector[i].v1.name == v2.name && edgesVector[i].v2.name == v1.name)) return 0;
		}
		Edge3D e;
		e.v1 = v1; e.v2 = v2;
		edgesVector.push_back(e);
		return 1;
	}

	int edgeLength(){
		return edgesVector.size();
	}

	int vertexLength(){
		return vertexVector.size();
	}

};


/** Function to obtain the rotation matrix using
  * the viewing direction obtained from the 
  * user. */

void getRotaionMatrixz(float x, float y, float z, double a[][3])
{
	if(x == 0 && y == 0)
	{
		a[2][0] = a[2][1] = a[0][1] = a[0][2] = a[1][0] = a[1][2] = 0;
		a[0][0] = a[1][1] = a[2][2] = 1;
		return;
	}
	float value = z/(pow(x*x+y*y+z*z, 0.5));
	float theta = acos(value);
	float qr, q1, q2, q3 = 0;
	qr = cos(theta/2);
	q1 = -y*sin(theta/2)/(pow(x*x+y*y, 0.5));
	q2 = x*sin(theta/2)/(pow(x*x+y*y, 0.5));
	a[0][0] = 1 - 2*(q2*q2 + q3*q3);
	a[0][1] = 2*(q1*q2 - q3*qr);
	a[0][2] = 2*(q1*q3 + q2*qr);
	a[1][0] = 2*(q1*q2 + q3*qr);
	a[1][1] = 1 - 2*(q1*q1 + q3*q3);
	a[1][2] = 2*(q2*q3 - q1*qr);
	a[2][0] = 2*(q1*q3 - q2*qr);
	a[2][1] = 2*(q2*q3 + q1*qr);
	a[2][2] = 1 - 2*(q1*q1 + q2*q2);
}

void getRotaionMatrixx(float x, float y, float z, double a[][3])
{
	if(y == 0 && z == 0)
	{
		a[2][0] = a[2][1] = a[0][1] = a[0][2] = a[1][0] = a[1][2] = 0;
		a[0][0] = a[1][1] = a[2][2] = 1;
		return;
	}
	float value = x/(pow(x*x+y*y+z*z, 0.5));
	float theta = acos(value);
	float qr, q1 = 0, q2, q3;
	qr = cos(theta/2);
	q2 = y*sin(theta/2)/(pow(z*z+y*y, 0.5));
	q3 = -z*sin(theta/2)/(pow(z*z+y*y, 0.5));
	a[0][0] = 1 - 2*(q2*q2 + q3*q3);
	a[0][1] = 2*(q1*q2 - q3*qr);
	a[0][2] = 2*(q1*q3 + q2*qr);
	a[1][0] = 2*(q1*q2 + q3*qr);
	a[1][1] = 1 - 2*(q1*q1 + q3*q3);
	a[1][2] = 2*(q2*q3 - q1*qr);
	a[2][0] = 2*(q1*q3 - q2*qr);
	a[2][1] = 2*(q2*q3 + q1*qr);
	a[2][2] = 1 - 2*(q1*q1 + q2*q2);
}

void getRotaionMatrixy(float x, float y, float z, double a[][3])
{
	if(x == 0 && z == 0)
	{
		a[2][0] = a[2][1] = a[0][1] = a[0][2] = a[1][0] = a[1][2] = 0;
		a[0][0] = a[1][1] = a[2][2] = 1;
		return;
	}
	float value = y/(pow(x*x+y*y+z*z, 0.5));
	float theta = acos(value);
	float qr, q1, q2 = 0, q3;
	qr = cos(theta/2);
	q1 = z*sin(theta/2)/(pow(x*x+z*z, 0.5));
	q3 = -x*sin(theta/2)/(pow(x*x+z*z, 0.5));
	a[0][0] = 1 - 2*(q2*q2 + q3*q3);
	a[0][1] = 2*(q1*q2 - q3*qr);
	a[0][2] = 2*(q1*q3 + q2*qr);
	a[1][0] = 2*(q1*q2 + q3*qr);
	a[1][1] = 1 - 2*(q1*q1 + q3*q3);
	a[1][2] = 2*(q2*q3 - q1*qr);
	a[2][0] = 2*(q1*q3 - q2*qr);
	a[2][1] = 2*(q2*q3 + q1*qr);
	a[2][2] = 1 - 2*(q1*q1 + q2*q2);
}

/** Function to get the projections of each vertex of the model 
  * when rotated using the specified rotation
  * matrix. */

Vertex2D getProjection(Vertex3D v, double rotationMatrix[][3], int count)
{
	Vertex2D v_new;
	v_new.name = v.name;
	double newVertices[3] = {0}, oldVertices[3] = {v.x, v.y, v.z};
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			newVertices[i] += rotationMatrix[i][j]*oldVertices[j];

	v_new.x = newVertices[1];
	v_new.y = newVertices[2];
	
	/*if(count == 0)
	{
		v_new.x = newVertices[1];
		v_new.y = newVertices[2];
	}
	else if(count == 1)
	{
		v_new.x = newVertices[1];
		v_new.y = newVertices[2];
	}
	else if(count == 2)
	{
		v_new.x = newVertices[1];
		v_new.y = newVertices[2];
	}*/
	return v_new;
}


/** Function to output the complete transformed model */

Model2D rotatedModel3D(Model3D model, Model2D outputModel, double rotationMatrix[][3], int count)
{
	Vertex3D v;
	Vertex2D v_new;
	int len = model.vertexLength();
	for(int i=0;i<len;i++)
	{
		v = model.vertexVector[i];
		v_new = getProjection(v, rotationMatrix, count);
		if(count == 0)
			outputModel.vertexVector1.push_back(v_new);
		else if(count == 1)
			outputModel.vertexVector2.push_back(v_new);
		else if(count == 2)
			outputModel.vertexVector3.push_back(v_new);
	}

	return outputModel;
}


/** Function to obtain a graphical interface to view the model 
  * after applying all transformations. */

void graphicalOutput(Model2D outputModel)
{
	ofstream file_v1("FrontViewVertices.txt"), file_v2("TopViewVertices.txt"), file_v3("SideViewVertices.txt");
	ofstream file_e1("FrontViewEdges.txt"), file_e2("TopViewEdges.txt"), file_e3("SideViewEdges.txt");
	
	int vertexLength = outputModel.vertex1Length(), edgeLength = outputModel.edge1Length();
	Vertex2D v; Edge2D e;
	
	for(int i=0;i<vertexLength;i++)
	{
		v = outputModel.vertexVector1[i];
		file_v1<<v.name<<" "<<v.x<<" "<<v.y<<"\n";
	}
	for(int i=0;i<edgeLength;i++)
	{
		e = outputModel.edgesVector1[i];
		file_e1<<e.v1.name<<" "<<e.v2.name<<"\n";
	}

	vertexLength = outputModel.vertex2Length();
	edgeLength = outputModel.edge2Length(); 
	
	for(int i=0;i<vertexLength;i++)
	{
		v = outputModel.vertexVector2[i];
		file_v2<<v.name<<" "<<v.x<<" "<<v.y<<"\n";
	}

	for(int i=0;i<edgeLength;i++)
	{
		e = outputModel.edgesVector2[i];
		file_e2<<e.v1.name<<" "<<e.v2.name<<"\n";
	}

	vertexLength = outputModel.vertex3Length();
	edgeLength = outputModel.edge3Length(); 
	
	for(int i=0;i<vertexLength;i++)
	{
		v = outputModel.vertexVector3[i];
		file_v3<<v.name<<" "<<v.x<<" "<<v.y<<"\n";
	}
	for(int i=0;i<edgeLength;i++)
	{
		e = outputModel.edgesVector3[i];
		file_e3<<e.v1.name<<" "<<e.v2.name<<"\n";
	}
}

int findIndex(Model2D outputModel, string n)
{
	int len = outputModel.vertex1Length();
	for(int i=0; i<len; i++)
	{
		if(outputModel.vertexVector1[i].name == n)
			return i;
	}
}

bool sameVertex(Vertex2D v1, Vertex2D v2)
{
	if(abs(v1.x - v2.x) < least && abs(v1.y - v2.y) < least)
		return 1;
	return 0;
}

/**Function to get 2D model from the 3D Model after
  * removing one of the corrdinates in each of the view. */

Model2D projectedModel(Model3D model, Model2D outputModel)
{
	Vertex2D v11, v12;
	Edge2D e1;
	int len = model.edgeLength();

	for(int i=0;i<len;i++)
	{
		Edge3D e = model.edgesVector[i];
		Vertex3D v1 = e.v1, v2 = e.v2;
		string n1 = v1.name, n2 = v2.name;
		
		int index1 = findIndex(outputModel, n1);
		int index2 = findIndex(outputModel, n2);

		v11 = outputModel.vertexVector1[index1];
		v12 = outputModel.vertexVector1[index2];
		if(!sameVertex(v11, v12))
		{	
			e1.v1 = v11; e1.v2 = v12;
			v11.vertex2DAddNeighbour(v12);
			v12.vertex2DAddNeighbour(v11);
			outputModel.edgesVector1.push_back(e1);
		}

		v11 = outputModel.vertexVector2[index1];
		v12 = outputModel.vertexVector2[index2];
		if(!sameVertex(v11, v12))
		{
			e1.v1 = v11; e1.v2 = v12;
			v11.vertex2DAddNeighbour(v12);
			v12.vertex2DAddNeighbour(v11);
			outputModel.edgesVector2.push_back(e1);
		}

		v11 = outputModel.vertexVector3[index1];
		v12 = outputModel.vertexVector3[index2];
		if(!sameVertex(v11, v12))
		{
			e1.v1 = v11; e1.v2 = v12;
			v11.vertex2DAddNeighbour(v12);
			v12.vertex2DAddNeighbour(v11);
			outputModel.edgesVector3.push_back(e1);
		}
	}

	return outputModel;
}

/** Function to generate the 2D representation of 
  * the vector of edges and vertices parsed as a file
  * along with viewing direction coordinates. */

void generate2D(Model3D model, float x, float y, float z)
{
	double rotationMatrix[3][3];
	Model2D outputModel;

	/// Obtains the rotation matrix from the viewing direction coordinates
	getRotaionMatrixx(x, y, z, rotationMatrix);

	/// Generate the rotated 3D model
	outputModel = rotatedModel3D(model, outputModel, rotationMatrix, 0);

	/// Obtains the rotation matrix from the viewing direction coordinates
	getRotaionMatrixy(x, y, z, rotationMatrix);

	/// Generate the rotated 3D model
	outputModel = rotatedModel3D(model, outputModel, rotationMatrix, 1);

	/// Obtains the rotation matrix from the viewing direction coordinates
	getRotaionMatrixz(x, y, z, rotationMatrix);

	/// Generate the rotated 3D model
	outputModel = rotatedModel3D(model, outputModel, rotationMatrix, 2);

	/// Creates a 2D model from rotated 3D Model
	outputModel = projectedModel(model, outputModel);

	/// Final graphical output in form of files
	graphicalOutput(outputModel);
}


/** Function to create a 3D model using the vertices and 
  * edges vector of both the 2D faces and store it in the 
  * model object. */
vector<string> tokenizeString(string str)
{
	string out = "";
	vector <string> tokens;
	
	for(int i=0;i<str.length();i++)
	{
		if(str[i]!=' ') out=out+str[i];
		else{
			tokens.push_back(out);
			out = "";
		}
	}
	tokens.push_back(out);
	return tokens;
}


/** Create a 3D Model from input .txt files of
  * vertices and edges for 2D file. */
Model3D createModel3D(ifstream& v1, ifstream& v2, ifstream& e1,  ifstream& e2)
{
	Model3D model;
	float x,y,z;
	string str,name;
	vector <Vertex3D> vertex;
	Vertex3D u,v;
	vector<string> vertexNames;
	vector<string> tokens;
	
	getline(v1, str);
	tokens = tokenizeString(str);
	string view1 = tokens[0],	 view2 = tokens[1];
	while(getline(v1, str)){
		
		tokens = tokenizeString(str);
		v.name = tokens[0];
		vertexNames.push_back(tokens[0]);
		if(view1 == "x" && view2 == "y")
		{ 
			v.x = stoi(tokens[1]);
			v.y = stoi(tokens[2]);
		}
		else if(view1 == "y" && view2 == "z")
		{ 
			v.y = stoi(tokens[1]);
			v.z = stoi(tokens[2]);
		}
		else if(view1 == "x" && view2 == "z")
		{
			v.x = stoi(tokens[1]); 
			v.z = stoi(tokens[2]);
		}

		vertex.push_back(v);
	}
	
	getline(v2, str);
	tokens = tokenizeString(str);
	view1 = tokens[0];
	view2 = tokens[1];
	
	while(getline(v2, str)){
		
		tokens = tokenizeString(str);
		int index = find(vertexNames.begin(), vertexNames.end(), tokens[0]) - vertexNames.begin();
		if(view1 == "x" && view2 == "y")
		{ 
			vertex[index].x = stoi(tokens[1]);
			vertex[index].y = stoi(tokens[2]);
		}
		else if(view1 == "y" && view2 == "z")
		{ 
			vertex[index].y = stoi(tokens[1]); 
			vertex[index].z = stoi(tokens[2]);
		}
		else if(view1 == "x" && view2 == "z")
		{
			vertex[index].x = stoi(tokens[1]);
			vertex[index].z = stoi(tokens[2]);
		}
	}
	model.vertexVector = vertex;
	while(getline(e1, str)){
	
		tokens = tokenizeString(str);
		int index1 = find(vertexNames.begin(), vertexNames.end(), tokens[0]) - vertexNames.begin();
		int index2 = find(vertexNames.begin(), vertexNames.end(), tokens[1]) - vertexNames.begin();
		u = vertex[index1];
		v = vertex[index2];
		int len = model.edgeLength();
		vertex[index1].vertex3DAddNeighbour(v);
		vertex[index2].vertex3DAddNeighbour(u);
		model.addEdge(u, v, len);
	
	}
	
	int length = model.edgeLength();
	while(getline(e2, str)){
	
		tokens = tokenizeString(str);
		int index1 = find(vertexNames.begin(), vertexNames.end(), tokens[0]) - vertexNames.begin();
		int index2 = find(vertexNames.begin(), vertexNames.end(), tokens[1]) - vertexNames.begin();
		u = vertex[index1];
		v = vertex[index2];
		if(model.addEdge(u, v, length))
		{
			vertex[index1].vertex3DAddNeighbour(v);
			vertex[index2].vertex3DAddNeighbour(u);
		}
	}
	return model;
}


/** Function to get vertices from a 3D model in a file. */

void getVertex(Model3D model, ofstream& vertex)
{
	int len = model.vertexLength();
	Vertex3D v;	
	for(int i=0;i<len;i++)
	{	
		v = model.vertexVector[i];
		vertex<<v.name<<" "<<v.x<<" "<<v.y<<" "<<v.z<<"\n";
	}
}


/** Function to get edges from a 3D model in a file. */

void getEdges(Model3D model, ofstream& edge)
{
	int len = model.edgeLength();
	Edge3D e;
	for(int i=0;i<len;i++)
	{
		e = model.edgesVector[i];
		edge<<e.v1.name<<" "<<e.v2.name<<"\n";
	}
}


/** Function to generate 3D model from given vertex
  * and edges vector of both faces along with viewing direction.*/

void generate3D(ifstream& v1, ifstream& v2, ifstream& e1,  ifstream& e2, float x, float y, float z)
{
	Model3D model;
	ofstream vertex("vertex.txt");
	ofstream edge("edge.txt");

	/// Function to create the 3D model
	model = createModel3D(v1, v2, e1, e2);
	getEdges(model, edge);
	getVertex(model, vertex);

	/// Reduce the problem to projecting from 3D to 2D views
	generate2D(model, x, y, z);	
}


/** Create a 3D Model from input set of files */

Model3D model3DFromInputFiles(ifstream& v, ifstream& e)
{
	Model3D model;
	Vertex3D v1,v2;
	Edge3D e1;
	string str;
	vector<string> tokens, vertexNames;

	while(getline(v, str))
	{
		tokens = tokenizeString(str);
		v1.name = tokens[0];
		v1.x = stoi(tokens[1]);
		v1.y = stoi(tokens[2]);
		v1.z = stoi(tokens[3]);
		vertexNames.push_back(tokens[0]);
		model.vertexVector.push_back(v1);
	}
	while(getline(e, str)){
		tokens = tokenizeString(str);
		int index1 = find(vertexNames.begin(), vertexNames.end(), tokens[0]) - vertexNames.begin();
		int index2 = find(vertexNames.begin(), vertexNames.end(), tokens[1]) - vertexNames.begin();
		v1 = model.vertexVector[index1];
		v2 = model.vertexVector[index2];
		model.vertexVector[index1].vertex3DAddNeighbour(v2);
		model.vertexVector[index2].vertex3DAddNeighbour(v1);
		model.addEdge(v1, v2, 0);
	}

	return model;
}


/** Function for interface in case of 3D to 2D
  * contains buttons to upload vertex and edge files.
  * Generate button begins processing to display the model. */

void interfaceFor3Dto2D()
{
	float x,y,z;
	cin>>x>>y>>z;

	ifstream v, e;
	v.open("vertices.txt");
	e.open("edges.txt");

	if (!v || !e) {
	    cerr << "Unable to open file datafile.txt";
	    exit(1);   // call system to stop
	}
	
	Model3D model = model3DFromInputFiles(v, e);
	generate2D(model, x, y, z);
}


/** Function for interface in case of 2D to 3D
  * contains buttons to upload vertex and edge files.
  * Generate button begins processing to display the model. */

void interfaceFor2Dto3D()
{
	float x,y,z;
	cin>>x>>y>>z;

	ifstream v1,e1,v2,e2;
	v1.open("TopViewVertices.txt");
	e1.open("TopViewEdges.txt");
	v2.open("SideViewVertices.txt");
	e2.open("SideViewEdges.txt");

	if (!v1 || !e1 || !v2 || !e2) {
	    cerr << "Unable to open file datafile.txt";
	    exit(1);   // call system to stop
	}
	generate3D(v1, v2, e1, e2, x, y, z);
}


/** Function for initial interface showing 2 buttons
  * for 3D to 2D and 2D to 3D conversion. */

void interfaceforCAD(int num)
{
	if(num == 0) interfaceFor2Dto3D();
	else if(num == 1) interfaceFor3Dto2D();
	else cout<<"Wrong Input\n";
}


/** Main function initializing buttons for 2D to 3D 
  * and 3D to 2D. */

int main()
{
	int num;
	cin>>num;
	interfaceforCAD(num);
}
