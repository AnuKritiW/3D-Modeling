#include "mesh.h"

#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include "math.h"
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "mesh.h"
#include "Vector3d.h"
#include <map>
#include <queue>
#include <iomanip>
using namespace std;

void myObjType::draw(bool m_smooth, bool m_edges) {

	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glPushMatrix();
	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	glScalef(4.0 / longestSide, 4.0 / longestSide, 4.0 / longestSide);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);
	
	if (m_edges) {
		if (hasBoundary()) {
			drawEdges();
			glDisable(GL_LIGHTING);
			glPopMatrix();	
			return;
		} else {
			cout << "There are no edges to draw" << endl;
		}
		
	}
	
	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);
// uncomment the following after you computed the normals
		if (!m_smooth) {
			glNormal3dv(nlist[i]);
		}
		
		for (int j = 0; j < 3; j++) {
			if (m_smooth) {
				glNormal3dv(vnlist[tlist[i][j]]);
			}
			glVertex3dv(vlist[tlist[i][j]]);
		}
		glEnd();
	
	}
	glDisable(GL_LIGHTING);

	glPopMatrix();
}

void myObjType::writeFile(char* filename) {
	ofstream outFile;
	outFile.open(filename);
	for (int i = 1; i <= vcount; i++) {
		outFile << "v " << vlist[i][0] << " " << vlist[i][1] << " " << vlist[i][2] << "\n";
	}
	for (int i = 1; i <= tcount; i++) {
		outFile << "f " << tlist[i][0] << " " << tlist[i][1] << " " << tlist[i][2] << "\n";
	}
}

void myObjType::readFile(char* filename)
{
	cout << "Opening " << filename << endl;
	ifstream inFile;
	inFile.open(filename);
	if (!inFile.is_open()) {
		cout << "We cannot find your file " << filename << endl;
		exit(1);
	}

	string line;
	int i, j;
	bool firstVertex = 1;
	double currCood;

	while (getline(inFile, line))
	{
		if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
		{
			if (line[0] == 'v')
			{
				vcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ') j++;
					currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
					if (firstVertex) 
						lmin[k] = lmax[k] = currCood;
					else {
						if (lmin[k] > currCood)
							lmin[k] = currCood;
						if (lmax[k] < currCood)
							lmax[k] = currCood;
					}
					i = j;
				}

				firstVertex = 0;
			}
			if (line[0] == 'f')
			{
				tcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) {
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ' && linec[j] != '\\') j++;
					tlist[tcount][k] = atof(line.substr(i, j - i).c_str());
					i = j;
					// fnlist[tcount][k] = 0;
					while (linec[j] != ' ') j++;

				}

			}
		}
	}
	

	// We suggest you to compute the normals here
	populate_nlist();
	populate_fnlist();
	populate_vnlist();
	//numTriOriented = 0;
	orientTriangles();

    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();

	compute_components();
}

void myObjType::populate_nlist() {
	for (int i = 1; i <= tcount; i++) {
		int v1_idx = tlist[i][0];
		int v2_idx = tlist[i][1];
		int v3_idx = tlist[i][2];

		Vector3d v1(vlist[v1_idx]);
		Vector3d v2(vlist[v2_idx]);
		Vector3d v3(vlist[v3_idx]);

		Vector3d normal_v = triNormal(v1, v2, v3);

		double x = normal_v.x();
		double y = normal_v.y();
		double z = normal_v.z();

		nlist[i][0] = x;
		nlist[i][1] = y;
		nlist[i][2] = z;
	}
}

bool myObjType::hasBoundary() {
	std::set<std::pair<int, int> > edges;
	for (int t = 1; t <= tcount; t++) {
		for (int v = 0; v < 3; v++) {
			int neighbor = fnvlist[t][v];
			if (neighbor == 0) {
				std::pair<int, int> edge = getVerts(t, v);
				edges.insert(std::make_pair(edge.first, edge.second));
			}
		}
	}
	if (edges.empty()) {
		return false;
	}
	return true;
}

int compute_orientation(OrTri orTri, map<int, int> v_map) {
	int tVersion = ver(orTri); // && (0x0000111b);
	tIdx t_idx = idx(orTri); // >> 3;
	return makeOrTri(t_idx, tVersion); // (tIdx << 3) | v_map[tVersion];
}

int myObjType::enext(OrTri orTri) { //given the orTri number

	// allows versions to be mapped correctly
	map<int, int> v_map;
	v_map[0] = 1;
	v_map[1] = 2;
	v_map[2] = 0;
	v_map[3] = 5;
	v_map[4] = 3;
	v_map[5] = 4;

	return compute_orientation(orTri, v_map);
}

int myObjType::sym(OrTri orTri) { //given the orTri number

	// allows versions to be mapped correctly
	map<int, int> v_map;
	v_map[0] = 3;
	v_map[1] = 4;
	v_map[2] = 5;
	v_map[3] = 0;
	v_map[4] = 1;
	v_map[5] = 2;
	
	return compute_orientation(orTri, v_map);
}

int myObjType::org(OrTri orTri) {

	// allows versions to be mapped correctly
	map<int, int> v_map;
	v_map[0] = 0;
	v_map[1] = 1;
	v_map[2] = 2;
	v_map[3] = 1;
	v_map[4] = 2;
	v_map[5] = 0;

	return compute_orientation(orTri, v_map);
}

int myObjType::dest(OrTri orTri) {
	return org(sym(orTri));
}

void myObjType::relax_mesh() {
	double x_sum, y_sum, z_sum, num_neighbors;
	std::set<int>::iterator it;
	for (int i = 0; i < adjVtoV.size(); i++) {
		/*if (boundary_verts.count(i)) {
			continue;
		}*/
		x_sum = 0;
		y_sum = 0;
		z_sum = 0;
		num_neighbors = 0;
		for (it = adjVtoV[i].begin(); it != adjVtoV[i].end(); ++it) {
			int neighbor_v = *it;
			x_sum += vlist[neighbor_v][0];
			y_sum += vlist[neighbor_v][1];
			z_sum += vlist[neighbor_v][2];
			num_neighbors++;
		}
		relaxed_vlist[i][0] = x_sum / num_neighbors;
		relaxed_vlist[i][1] = y_sum / num_neighbors;
		relaxed_vlist[i][2] = z_sum / num_neighbors;
	}

	for (int i = 1; i < adjVtoV.size(); i++) {
		vlist[i][0] = relaxed_vlist[i][0];
		vlist[i][1] = relaxed_vlist[i][1];
		vlist[i][2] = relaxed_vlist[i][2];
	}

	populate_nlist();
	populate_vnlist();

	cout << "Completed Mesh Relaxation!" << endl;
}

void myObjType::populate_fnlist() {
	// create an ortri with each triangle in tlist => key in fnlist	
	adjVtoV = {};
	adjVtoE = {};
	adjFtoE = {};
	adjFtoV = {};
	adjFtoF = {};
    // for every triangle
    for (int t = 1; t <= tcount; t++) {
        //for every vertex from 0-2
        for (int v = 0; v < 3; v++) { //only considering first 3 versions; opposite faces are found later
            int v1 = tlist[t][v];
			int v2 = tlist[t][(v+1)%3]; //mod 3 so that the computation works for 2nd and 3rd indices as well
			int v3 = tlist[t][(v+2)%3];

			adjVtoV[v1].insert(v2);
			adjVtoV[v1].insert(v3);
			
			std::set<int> key = {v1, v2}; //key is the edge formed by the first two vertices
			tIdx t_idx = t;
			OrTri orTri = makeOrTri(t_idx, v); //convert current triangle into an oriented triangle

			adjVtoE[key].insert(v3);
			adjFtoE[key].insert(orTri);

			adjFtoV[v1].insert(t); //strong indices to average vertex normals for smoothing
        }
    }

	// for every triangle
    for (int t = 1; t <= tcount; t++) {
        //for every vertex from 0-2
        for (int v = 0; v < 3; v++) { //only considering first 3 versions; opposite faces are found later
			int v1 = tlist[t][v];
			int v2 = tlist[t][(v+1)%3]; //mod 3 so that the computation works for 2nd and 3rd indices as well
			int v3 = tlist[t][(v+2)%3];
			
			set<int> key = {v1, v2}; //key is the edge formed by the first two vertices
			
			tIdx t_idx = t;
			OrTri orTri = makeOrTri(t_idx, v); //convert current triangle into an oriented triangle

			std::set<OrTri> adj_faces = adjFtoE[key];
			
			OrTri opp_face = 0; //will be 0 if no opposite face is found (i.e. current face is a boundary)

			// set opp_face to be the face opposite to the current face.
			set<OrTri>::iterator it;
			for (it = adj_faces.begin(); it != adj_faces.end(); ++it) {
				OrTri face = *it;
				if (face != orTri) {
					opp_face = face;
				}

			}

			fnvlist[t][v] = opp_face; //this stores the versions for each tri
			fnlist[orTri] = opp_face; //opp_face is also an orTri

			if (opp_face != 0) {
				adjFtoF[t].insert(opp_face >> 3);
			}

		}
	}
}

OrTri myObjType::fnext(OrTri orTri) {
	return fnlist[orTri];
}

// for smoothing
void myObjType::populate_vnlist() {
	//for each vertex
	for (int v_idx = 1; v_idx <= vcount; v_idx++) {
		std::set<int> adj_faces = adjFtoV[v_idx];
		Vector3d sum_vec(0, 0, 0);
		
		// for each adjacent face
		set<int>::iterator it;
		for (it = adj_faces.begin(); it != adj_faces.end(); ++it) {
			int face_idx = *it;
			sum_vec[0] += nlist[face_idx][0];
			sum_vec[1] += nlist[face_idx][1];
			sum_vec[2] += nlist[face_idx][2];
		}

		vnlist[v_idx][0] = sum_vec[0] / sum_vec.length();
		vnlist[v_idx][1] = sum_vec[1] / sum_vec.length();
		vnlist[v_idx][2] = sum_vec[2] / sum_vec.length();
	}
}

std::pair<int, int> myObjType::getVerts(int t, int v) {
	int v0;
	int v1;
	if (v == 0) {
		v0 = tlist[t][0];
		v1 = tlist[t][1];
	} else if (v == 1) {
		v0 = tlist[t][1];
		v1 = tlist[t][2];
	} else { //v==2
		v0 = tlist[t][2];
		v1 = tlist[t][0];
	}
	std::pair<int, int> edge(v0, v1);
	return edge;
}

void myObjType::drawEdges() {
	//boundary_verts = {};
	glDisable(GL_LIGHTING);
	for (int t = 1; t <= tcount; t++) {
		for (int v = 0; v < 3; v++) {
			int neighbor = fnvlist[t][v]; //get index of tri
			if (neighbor == 0) { //i.e. boundary --> only want to draw edges for boundary edges

				std::pair<int, int> edge = getVerts(t, v);
				//boundary_verts.insert(edge.first);
				//boundary_verts.insert(edge.second);

				glBegin(GL_LINES);
				glColor3f(1.0f, 0.0f, 0.0f);
				glVertex3dv(vlist[edge.first]);
				glVertex3dv(vlist[edge.second]);
				glEnd();
			}
		}
		
	}
}

int myObjType::getUnvisitedFace(std::set<int> visited_faces) {
	for (int i = 1; i <= tcount; i++) {
		if (visited_faces.find(i) == visited_faces.end()) { //didn't find face i to be visited
			return i;
		}
	}
	return -1;
}

void myObjType::compute_components() {
	int num_comp = 0;
	std::set<int> visited_faces;
	std::queue<int> unvisited_faces;

	while (visited_faces.size() < tcount) { //for each root triangle in each unique comp
		int unvisited_face = getUnvisitedFace(visited_faces);
		unvisited_faces.push(unvisited_face); //push all neighbouring triangles

		while (!unvisited_faces.empty()) { //while there are still unvisited faces
			int visited_face = unvisited_faces.front();
			unvisited_faces.pop();

			if (visited_faces.find(visited_face) == visited_faces.end()) {
				visited_faces.insert(visited_face); //if not already in set, add

				std::set<int>::iterator it;
				for (it = adjFtoF[visited_face].begin(); it != adjFtoF[visited_face].end(); ++it) {
					int neighboring_face = *it;
					unvisited_faces.push(neighboring_face);
				}
			}
		}
		unvisited_faces = {};
		num_comp++;
	}
	cout << "No. of components: " << num_comp << endl;
}


bool myObjType::orientTriangles() {
	numTriOriented = 0;
	unvisited_faces_g = {}; //this holds verts that still need to be traversed.
	visited_faces_g = {};
	while (visited_faces_g.size() < tcount) { //for each root triangle in each unique comp
		int unvisited_face = getUnvisitedFace(visited_faces_g); //get the first index that has not been visited yet
		unvisited_faces_g.push(unvisited_face);
		visited_faces_g.insert(unvisited_face);

		while (!unvisited_faces_g.empty()) { //while there are still unvisited faces
			int visited_face = unvisited_faces_g.front();
			unvisited_faces_g.pop();

			bool neighbor_face = orientTriangles_helper(visited_face);
			if (!neighbor_face) {
				return false;
			}
		}
		unvisited_faces_g = {};
	}
	if (numTriOriented > 0) {
		cout << "No. of reoriented triangles: " << numTriOriented << endl;
		populate_fnlist();
		populate_nlist();
		populate_vnlist();
	}
	return true;
}

bool myObjType::orientTriangles_helper(int visited_face) {
	bool sameOrientation = false;
	int neighbor_face;
	for (int v = 0; v < 3; v++) {
		sameOrientation = false;
		int neighbor = fnvlist[visited_face][v]; //get each version from 0 to 2
		if (neighbor != 0) {
			neighbor_face = neighbor >> 3;
			int neighbor_v = neighbor & 0b111;
			std::pair<int, int> visited_face_verts = getVerts(visited_face, v);
			std::pair<int, int> neighbor_verts = getVerts(neighbor_face, neighbor_v);
			if (visited_face_verts == neighbor_verts) { //they are oriented in the same way
				sameOrientation = true;
			}
			if ((visited_faces_g.find(neighbor_face) != visited_faces_g.end()) && sameOrientation) { //neighbour face has been visited and it is in the same orientation
				cout << "cannot orient triangles!" << endl;
				return false;
			}
			else if (sameOrientation) { //i.e. neighbor face is still unvisited, but is in the same orientation
				//orient (exchange vertices 1 and 2)
				int tmp = tlist[neighbor_face][1];
				tlist[neighbor_face][1] = tlist[neighbor_face][2];
				tlist[neighbor_face][2] = tmp;

				//update fnvlist (exchange versions 0 and 2)
				int tmp_f = fnvlist[neighbor_face][0];
				fnvlist[neighbor_face][0] = fnvlist[neighbor_face][2];
				fnvlist[neighbor_face][2] = tmp_f;
				numTriOriented++;
			}
			if (!(visited_faces_g.find(neighbor_face) != visited_faces_g.end())) {
				unvisited_faces_g.push(neighbor_face);
				visited_faces_g.insert(neighbor_face);
			}
		}
	}
}

std::vector<Vector3d> myObjType::getUpdatedVerts(int i) { //i is the triangle index
	std::vector<Vector3d> all_verts; //vector of the all vertices in each triangle, incl new
	Vector3d new_v;

	for (int v = 0; v < 3; v++) { //for every pair of vertices in each triangle
		int v1 = tlist[i][v];
		int next_v = (v + 1) % 3;
		int v2 = tlist[i][next_v];

		// vector at each vertex of the edge
		Vector3d ec1(vlist[v1][0], vlist[v1][1], vlist[v1][2]); //get coordinates and make vector
		Vector3d ec2(vlist[v2][0], vlist[v2][1], vlist[v2][2]);

		new_v = ((ec1 + ec2) / 2.0);
		all_verts.push_back(new_v);
	}
	int v1 = tlist[i][0];
	int v2 = tlist[i][1];
	int v3 = tlist[i][2];

	//get vector of each vertex of triangle using coordinates
	Vector3d c1(vlist[v1][0], vlist[v1][1], vlist[v1][2]);
	Vector3d c2(vlist[v2][0], vlist[v2][1], vlist[v2][2]);
	Vector3d c3(vlist[v3][0], vlist[v3][1], vlist[v3][2]);

	all_verts.push_back(c1);
	all_verts.push_back(c2);
	all_verts.push_back(c3);

	//get centroid
	Vector3d centroid((c1 + c2 + c3) / 3.0);
	all_verts.push_back(centroid);

	return all_verts;
}

void myObjType::update_subdiv_tlist(Vector3d t) {
	subdiv_tcount++;
	subdiv_tlist[subdiv_tcount][0] = t[0];
	subdiv_tlist[subdiv_tcount][1] = t[1];
	subdiv_tlist[subdiv_tcount][2] = t[2];
}

bool myObjType::are_same(double a, double b, double epsilon) {
	return fabs(a - b) < epsilon;
}

void myObjType::barycentric_subdivision() {
	for (int i = 1; i <= tcount; i++) { //for each triangle
		std::vector<Vector3d> all_verts = getUpdatedVerts(i); // will have size == 7
		
		std::vector<int> all_vert_idxs; //will store all vertex indices
		// for each of the vertices, need to update vlist
		for (int j = 0; j < all_verts.size(); j++) {

			bool tmp_bool;
			int tmp_int;
			bool almost_equal_bool = false;
			for (int k = 1; k <= subdiv_vcount; k++) {
				//Vector3d v = all_ve0rts[j];
				double epsilon = 0.0002;
				if (myObjType::are_same(subdiv_vlist[k][0], all_verts[j][0], epsilon) && myObjType::are_same(subdiv_vlist[k][1], all_verts[j][1], epsilon) && myObjType::are_same(subdiv_vlist[k][2], all_verts[j][2], epsilon)) {
					tmp_bool = false;
					tmp_int = k;
					almost_equal_bool = true;
				}
			}
			if (!almost_equal_bool) {
				subdiv_vlist[subdiv_vcount + 1][0] = all_verts[j][0];
				subdiv_vlist[subdiv_vcount + 1][1] = all_verts[j][1];
				subdiv_vlist[subdiv_vcount + 1][2] = all_verts[j][2];
				tmp_bool = true;
				tmp_int = subdiv_vcount + 1;
			}
			
			all_vert_idxs.push_back(tmp_int);
			if (tmp_bool) {
				subdiv_vcount++;
			}
		}

		Vector3d t1(all_vert_idxs[6], all_vert_idxs[3], all_vert_idxs[0]);
		Vector3d t2(all_vert_idxs[6], all_vert_idxs[0], all_vert_idxs[4]);
		Vector3d t3(all_vert_idxs[6], all_vert_idxs[4], all_vert_idxs[1]);
		Vector3d t4(all_vert_idxs[6], all_vert_idxs[1], all_vert_idxs[5]);
		Vector3d t5(all_vert_idxs[6], all_vert_idxs[5], all_vert_idxs[2]);
		Vector3d t6(all_vert_idxs[6], all_vert_idxs[2], all_vert_idxs[3]);
		
		update_subdiv_tlist(t1);
		update_subdiv_tlist(t2);
		update_subdiv_tlist(t3);
		update_subdiv_tlist(t4);
		update_subdiv_tlist(t5);
		update_subdiv_tlist(t6);
		
	}
	//update global counts and lists with new information
	
	std::copy(&subdiv_tlist[0][0], &subdiv_tlist[0][0] + subdiv_tcount * 3 + 3, &tlist[0][0]);
	std::copy(&subdiv_vlist[0][0], &subdiv_vlist[0][0] + subdiv_vcount * 3 + 3, &vlist[0][0]);
	tcount = subdiv_tcount;
	vcount = subdiv_vcount;
	subdiv_tcount = 0;
	subdiv_vcount = 0;
	populate_fnlist();
	populate_nlist();
	populate_vnlist();
	compute_components();

	cout << "Completed subdivision!" << endl;
}



void myObjType::computeStat()
{
	int i;
    double minAngle = 0;
    double maxAngle = 0;
	double globalMin = 200;
	double globalMax = -1;

	for (int i = 0; i < 18; i++) {
		statMinAngle[i] = 0;
		statMaxAngle[i] = 0;
	}
	 
	for (int i = 1; i <= tcount; i++) {
		int v1_idx = tlist[i][0];
		int v2_idx = tlist[i][1];
		int v3_idx = tlist[i][2];

		Vector3d v1(vlist[v1_idx]);
		Vector3d v2(vlist[v2_idx]);
		Vector3d v3(vlist[v3_idx]);

		Vector3d v1_v2(v2 - v1);
		Vector3d v1_v3(v3 - v1);
		Vector3d v2_v3(v3 - v2);
		
		/*double v1_dot_v2 = dot(v1, v2);
		double v2_dot_v3 = dot(v2, v3);*/

		double v1_dot = dot(v1_v2, v1_v3);
		double v2_dot = dot(-v1_v2, v2_v3);
		//double v3_dot = dot(-v1_v3, -v2_v3);

		double v1_v2_len = v1_v2.length();
		double v1_v3_len = v1_v3.length();
		double v2_v3_len = v2_v3.length();

		//double temp = v1_dot_v2 / (v1_len * v2_len);
		double v1_angle = acos(v1_dot / (v1_v2_len * v1_v3_len)) * 180.0 / M_PI;
		double v2_angle = acos(v2_dot / (v1_v2_len * v2_v3_len)) * 180.0 / M_PI;
		double v3_angle = 180 - v1_angle - v2_angle;

		minAngle = min(v1_angle, v2_angle);
		minAngle = min(minAngle, v3_angle);
		globalMin = min(globalMin, minAngle);

		maxAngle = max(v1_angle, v2_angle);
		maxAngle = max(maxAngle, v3_angle);
		globalMax = max(globalMax, maxAngle);

		int min_bin = int(floor(minAngle / 10));
		statMinAngle[min_bin] += 1;

		int max_bin = int(floor(maxAngle / 10));
		statMaxAngle[max_bin] += 1;

	}
    
    cout << "Min. angle = " << globalMin << endl;
    cout << "Max. angle = " << globalMax << endl;

	cout << "Statistics for Maximum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMaxAngle[i] << " ";
	cout << endl;
	cout << "Statistics for Minimum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMinAngle[i] << " ";
	cout << endl;
}
