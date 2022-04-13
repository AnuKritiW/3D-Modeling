#pragma once

// maximum number of vertices and triangles
#define MAXV 1000000
#define MAXT 1000000
#include <map>
#include <set>
#include <queue>
#include <cmath>
#include <limits>
#include <iomanip>
#include <iostream>
#include <type_traits>
#include <algorithm>
#include "Vector3d.h"

typedef int OrTri;
typedef int tIdx;

inline OrTri makeOrTri(tIdx t, int version) { return (t << 3) | version; };
inline tIdx idx(OrTri ot) { return ot >> 3; };
inline int ver(OrTri ot) { return ot & 0b111; };

class myObjType {
	int vcount = 0;
	int tcount = 0;

	double vlist[MAXV][3];   // vertices list
	int tlist[MAXT][3];      // triangle list
	int fnvlist[MAXT][3];     // fnext list with versions
	OrTri fnlist[MAXT];		//fnext list with ortri only
	double nlist[MAXT][3];   // storing triangle normals
	double vnlist[MAXT][3];  //storing vertex normals for smoothing

	double lmax[3];          // the maximum coordinates of x,y,z
	double lmin[3];          // the minimum coordinates of x,y,z

	int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
	int statMaxAngle[18]; 
	int numTriOriented;

	double subdiv_vlist[MAXV][3];
	int subdiv_tlist[MAXT][3];
	int subdiv_vcount = 0;
	int subdiv_tcount = 0;
	
	double relaxed_vlist[MAXV][3];

public:
	myObjType() { vcount = 0; tcount = 0; };
	void readFile(char* filename);  // assumming file contains a manifold
	void writeFile(char* filename);  
	void draw(bool m_smooth, bool m_edges);
    void computeStat();
	OrTri enext(OrTri orTri);
	OrTri sym(OrTri orTri);
	OrTri org(OrTri orTri);
	OrTri dest(OrTri orTri);
	void myObjType::relax_mesh();
	void populate_fnlist();
	OrTri fnext(OrTri orTri);
	void populate_vnlist();
	void compute_components();
	int myObjType::getUnvisitedFace(std::set<int> visited_faces);
	bool myObjType::hasBoundary();
	std::pair<int, int> myObjType::getVerts(int t, int v);
	void myObjType::drawEdges();
	bool myObjType::orientTriangles();
	bool myObjType::orientTriangles_helper(int unvisited_face);
	void myObjType::populate_nlist();
	void myObjType::barycentric_subdivision();
	std::vector<Vector3d> myObjType::getUpdatedVerts(int i); //i is the triangle index
	void myObjType::update_subdiv_tlist(Vector3d t);
	bool myObjType::are_same(double a, double b, double epsilon);
	
	
private:
	//int compute_orientation(int orTri, map<int, int> v_map);
	std::map<std::set<int>, std::set<OrTri> > adjFtoE;
	std::map<std::set<int>, std::set<OrTri> > adjVtoE;
	std::map<int, std::set<int> > adjFtoV;
	std::map<int, std::set<int> > adjVtoV;
	std::map<int, std::set<int> > adjFtoF;
	std::set<int> visited_faces_g;
	std::queue<int> unvisited_faces_g;
	//std::set<int> boundary_verts;
};


