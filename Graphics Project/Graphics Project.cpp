// Graphics Project.cpp : Defines the entry point for the application.
//

#include "stdafx.h"
#include "math.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "Graphics Project.h"

#define MAX_LOADSTRING 100

using namespace std;

#define ID_Save 1
#define ID_Load 2
#define ID_Exit 3
#define ID_LineDDA 4
#define ID_LineMidPoint 5
#define ID_LineParametric 6 
#define ID_CircleCartesian 8
#define ID_CirclePolar 9
#define ID_CircleIterativePolar 10
#define ID_CircleMidPoint 11
#define ID_CurveFirstDegree 12
#define ID_CurveSecondDegree 13
#define ID_CurveThirdDegreeHermite 14
#define ID_CurveThirdDegreeBezier 15
#define ID_CurveSplines 16
#define ID_ConvixFilling 17
#define ID_RectangleClippingPoint 18
#define ID_RectangleClippingLine 19
#define ID_CircleClippingPoint 20
#define ID_CircleClippingLine 21
#define ID_RedColor 22
#define ID_GreenColor 23
#define ID_BlueColor 24
#define ID_BlackColor 25

bool OK = false;
bool _LineDDA = false;
bool _LineMidPoint = false;
bool _LineParametric = false;
bool _CircleCartesian = false;
bool _CirclePolar = false;
bool _CircleIterativePolar = false;
bool _CircleMidPoint = false;
bool _CurveFirstDegree = false;
bool _CurveSecondDegree = false;
bool _CurveThirdDegreeHermite = false;
bool _CurveThirdDegreeBezier = false;
bool _CurveSplines = false;
bool _ConvixFilling = false;
bool _RectangleClippingPoint = false;
bool _RectangleClippingLine = false;
bool _CircleClippingPoint = false;
bool _CircleClippingLine = false;
bool _TwoPoints = false;
bool _OnePoint = false;

void init() {
	_LineDDA = false;
	_LineMidPoint = false;
	_LineParametric = false;
	_CircleCartesian = false;
	_CirclePolar = false;
	_CircleIterativePolar = false;
	_CircleMidPoint = false;
	_CurveFirstDegree = false;
	_CurveSecondDegree = false;
	_CurveThirdDegreeHermite = false;
	_CurveThirdDegreeBezier = false;
	_CurveSplines = false;
	_ConvixFilling = false;
	_RectangleClippingPoint = false;
	_RectangleClippingLine = false;
	_CircleClippingPoint = false;
	_CircleClippingLine = false;
	_TwoPoints = false;
	_OnePoint = false;
	OK = false;
}
// Global Variables:
HINSTANCE hInst;								// current instance
TCHAR szTitle[MAX_LOADSTRING];					// The title bar text
TCHAR szWindowClass[MAX_LOADSTRING];			// the main window class name

vector<string> shapes;
HDC hdc;
HWND hWnd;
COLORREF color = RGB(255, 255, 255);

// Forward declarations of functions included in this code module:
ATOM				MyRegisterClass(HINSTANCE hInstance);
BOOL				InitInstance(HINSTANCE, int);
LRESULT CALLBACK	WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK	About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY _tWinMain(_In_ HINSTANCE hInstance,
	_In_opt_ HINSTANCE hPrevInstance,
	_In_ LPTSTR    lpCmdLine,
	_In_ int       nCmdShow)
{
	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(lpCmdLine);

	// TODO: Place code here.
	MSG msg;
	HACCEL hAccelTable;

	// Initialize global strings
	LoadString(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
	LoadString(hInstance, IDC_GRAPHICSPROJECT, szWindowClass, MAX_LOADSTRING);
	MyRegisterClass(hInstance);

	// Perform application initialization:
	if (!InitInstance(hInstance, nCmdShow))
	{
		return FALSE;
	}

	hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_GRAPHICSPROJECT));

	// Main message loop:
	while (GetMessage(&msg, NULL, 0, 0))
	{
		if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}

	return (int)msg.wParam;
}


/// ========================================================= Other Functions and Class Pointt =================================================
class Pointt {
public:
	int x;
	int y;
	Pointt() {}
	Pointt(int a, int b) {
		x = a;
		y = b;
	}
};

string ConvertToString(int n) {
	string Result;
	ostringstream convert;
	convert << n;
	Result = convert.str();
	return Result;
}
int ConvertToInteger(string s) {
	istringstream buffer(s);
	int value;
	buffer >> value;
	return value;
}

string FillString(string str, int xs, int ys, int xe, int ye) {
	string s = str + "-";
	string n = ConvertToString(xs);
	s += n;
	n = ConvertToString(ys);
	s += "-" + n;
	n = ConvertToString(xe);
	s += "-" + n;
	n = ConvertToString(ye);
	s += "-" + n;
	s += "-";
	return s;
}
string FillString2(string str, Pointt *p, int size) {
	string s = str + "-";
	for (int i = 0; i < size; i++) {
		string n = ConvertToString(p[i].x);
		string n1 = ConvertToString(p[i].y);
		s += (n + "-" + n1 + "-");
	}
	return s;
}
string GetType(string shape) {
	string s = "";
	for (int i = 0; i < shape.length(); i++) {
		if (shape[i] != '-')
			s += shape[i];
		else
			return s;
	}
}
vector<int> GetPoints(string shape) {
	vector<int> result;
	string s = "";
	for (int i = 0; i < shape.length(); i++) {
		if (shape[i] != '-')
			s += shape[i];
		else {
			result.push_back(ConvertToInteger(s));
			s = "";
		}
	}
	return result;
}
void swap(int &x, int &y) {
	int tmp = x;
	x = y;
	y = tmp;
}
int round(int x) {
	return (int)(x + 0.5);
}
/// ========================================================= Background Color ========================================================
void ChangeBackground(HWND hWnd, COLORREF color) {
	RECT rect;
	HDC hdc = GetDC(hWnd);
	HBRUSH brush = CreateSolidBrush(color);
	GetWindowRect(hWnd, &rect);
	FillRect(hdc, &rect, brush);
	DeleteObject(brush);
	ReleaseDC(hWnd, hdc);
}
/// ========================================================= Line Algorithms ========================================================
void DDADrawLine(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {
	int dx = xe - xs;
	int dy = ye - ys;
	double slope = (double)dy / dx;
	if (abs(dy) < abs(dx)) { // if slope < 1
		if (xs > xe) {
			swap(xs, xe);
			swap(ys, ye);
		}
		double x = xs;
		double y = ys;
		SetPixel(hdc, xs, ys, color);
		while (x < xe) {
			x = x + 1;
			y = y + slope;
			SetPixel(hdc, x, round(y), color);
		}
	}
	else {
		int islope = 1;
		if (dy != 0)
			islope = dx / dy;
		if (ys > ye) {
			swap(xs, xe);
			swap(ys, ye);
		}
		double x = xs;
		double y = ys;
		SetPixel(hdc, xs, ys, color);
		while (y < ye) {
			y = y + 1;
			x = x + 1 / slope;
			SetPixel(hdc, round(x), y, color);
		}
	}
}
void ParametricDrawLine(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {
	int B1 = xs;
	int B2 = ys;
	int A1 = xe - xs;
	int A2 = ye - ys;
	int n = max(abs(A1), abs(A2));
	double dt = 1.0 / n;
	double x = (double)xs;
	double y = (double)ys;
	double t = 0;
	for (int i = 0; i < n; i++) {
		SetPixel(hdc, round(x), round(y), color);
		x = A1 *t + B1;
		y = A2 *t + B2;
		t += dt;
	}
}
void MidPointLineDraw(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {
	int dx = xe - xs;
	int dy = ye - ys;
	if (abs(dy) <= abs(dx)) { // slope < 1
		if (xs > xe) {
			swap(xs, xe);
			swap(ys, ye);
			dx *= -1;
			dy *= -1;
		}
		int incY = (ys > ye ? -1 : 1);
		int d = dx - 2 * abs(dy);
		int change1 = 2 * (dx - abs(dy));
		int change2 = -2 * abs(dy);
		int x = xs, y = ys;
		SetPixel(hdc, x, y, color);
		while (x < xe) {
			if (d < 0) {
				y = y + incY;
				d = d + change1;
			}
			else {
				d = d + change2;
			}
			x++;
			SetPixel(hdc, x, y, color);
		}
	}
	else {
		if (ys > ye) {
			swap(xs, xe);
			swap(ys, ye);
			dx *= -1;
			dy *= -1;
		}
		int incX = (xs > xe ? -1 : 1);
		int d = 2 * abs(dx) - dy;
		int change1 = 2 * (abs(dx) - dy);
		int change2 = 2 * abs(dx);
		int x = xs, y = ys;
		SetPixel(hdc, x, y, color);
		while (y < ye) {
			if (d > 0) {
				x = x + incX;
				d = d + change1;
			}
			else {
				d = d + change2;
			}
			y++;
			SetPixel(hdc, x, y, color);
		}
	}
}
// ====================================== Circle Drawing Algorithms =============================================
void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color) {
	SetPixel(hdc, xc + a, yc + b, color);
	SetPixel(hdc, xc - a, yc + b, color);
	SetPixel(hdc, xc - a, yc - b, color);
	SetPixel(hdc, xc + a, yc - b, color);
	SetPixel(hdc, xc + b, yc + a, color);
	SetPixel(hdc, xc - b, yc + a, color);
	SetPixel(hdc, xc - b, yc - a, color);
	SetPixel(hdc, xc + b, yc - a, color);
}
void CartesianDraw(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
	double R = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	int x = 0;
	double y = R;

	while (x < y) {
		Draw8Points(hdc, x1, y1, x, round(y), color);
		x++;
		y = sqrt((R*R) - (x*x));
	}
}
void PolarDraw(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
	double R = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	double x = R;
	double y = 0;

	double theta = 0.0;
	double deltaTheta = 1.0 / R;

	while (x > y) {
		theta += deltaTheta;
		x = round(R*cos(theta));
		y = round(R*sin(theta));
		Draw8Points(hdc, x1, y1, round(x), round(y), color);
	}
}
void IterativePolarDrawCircle(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
	double R = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	double x = R, y = 0;
	double dtheta = 1.0 / R;
	double cosdtheta = cos(dtheta), sindtheta = sin(dtheta);
	Draw8Points(hdc, x1, y1, R, 0, color);

	while (x > y)
	{
		double x11 = x*cosdtheta - y*sindtheta;
		y = x*sindtheta + y*cosdtheta;
		x = x11;
		Draw8Points(hdc, x1, y1, round(x), round(y), color);
	}
}
void MidpointDrawCircle(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
	double R = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	int x = 0;
	int y = R;
	int d = 1 - R;
	int change1 = 3;
	int change2 = 5 - 2 * R;
	Draw8Points(hdc, x1, y1, x, y, color);
	while (x < y) {
		if (d < 0) {
			d += change1;
			change2 += 2;
		}
		else {
			d += change2;
			change2 += 4;
			y--;
		}
		change1 += 2;
		x++;
		Draw8Points(hdc, x1, y1, x, y, color);
	}
}

// ========================================= Clipping According to Rectangle =========================================
void Rectangle_PointClipping(HDC hdc, int x, int y, int xleft, int ytop, int xright, int ybottom, COLORREF color)
{
	if (x >= xleft && x <= xright && y >= ytop && y <= ybottom)
		SetPixel(hdc, x, y, color);
}
union outCode
{
	unsigned All : 4;
	struct { unsigned left : 1, right : 1, top : 1, bottom : 1; };
};
outCode GetOutCode(double x, double y, int xleft, int xright, int ytop, int ybottom)
{
	outCode res;
	res.All = 0;
	if (x<xleft)
		res.left = 1;
	else if (x>xright)
		res.right = 1;
	if (y<ytop)
		res.top = 1;
	else if (y>ybottom)
		res.bottom = 1;
	return res;
}
struct Point {
	double x, y;
	Point() {}
	Point(double a, double b) {
		x = a;
		y = b;
	}
};
Point VIntersect(double xs, double ys, double xe, double ye, int xEdge) {
	Point p;
	p.x = xEdge;
	p.y = ys + (xEdge - xs)*(ye - ys) / (xe - xs);
	return p;
}
Point HIntersect(double xs, double ys, double xe, double ye, int yEdge) {
	Point p;
	p.y = yEdge;
	p.x = xs + (yEdge - ys)*(xe - xs) / (ye - ys);
	return p;
}
void CohenSuth(HDC hdc, int xs, int ys, int xe, int ye, int xleft, int xright, int ytop, int ybottom, COLORREF color) {
	outCode out1 = GetOutCode(xs, ys, xleft, xright, ytop, ybottom);
	outCode out2 = GetOutCode(xe, ye, xleft, xright, ytop, ybottom);
	while ((out1.All || out2.All) && !(out1.All & out2.All)) {
		Point v;
		if (out1.All) {
			if (out1.left) {
				v = VIntersect(xs, ys, xe, ye, xleft);
			}
			else if (out1.right) {
				v = VIntersect(xs, ys, xe, ye, xright);
			}
			else if (out1.top) {
				v = HIntersect(xs, ys, xe, ye, ytop);
			}
			else if (out1.bottom) {
				v = HIntersect(xs, ys, xe, ye, ybottom);
			}
			xs = v.x, ys = v.y;
			out1 = GetOutCode(xs, ys, xleft, xright, ytop, ybottom);
		}
		else {
			if (out2.left) {
				v = VIntersect(xs, ys, xe, ye, xleft);
			}
			else if (out2.right) {
				v = VIntersect(xs, ys, xe, ye, xright);
			}
			else if (out2.top) {
				v = HIntersect(xs, ys, xe, ye, ytop);
			}
			else if (out2.bottom) {
				v = HIntersect(xs, ys, xe, ye, ybottom);
			}
			xe = v.x, ye = v.y;
			out2 = GetOutCode(xe, ye, xleft, xright, ytop, ybottom);
		}
	}
	if (!out1.All && !out2.All) {
		DDADrawLine(hdc, xs, ys, xe, ye, color);
	}
}
// =========================================  Clipping According to Circle =========================================
void Circle_PointClipping(HDC hdc, int X, int Y, int xc, int yc, int x, int y, COLORREF color)
{
	double R = sqrt((xc - x)*(xc - x) + (yc - y)*(yc - y));
	double distanceToCenter = sqrt((xc - X)*(xc - X) + (yc - Y)*(yc - Y));
	if (distanceToCenter <= R)
		SetPixel(hdc, x, y, color);
}
void Circle_LineClipping(HDC hdc, int xs, int ys, int xe, int ye, int xc, int yc, int x, int y, COLORREF color)
{
	double R = sqrt((xc - x)*(xc - x) + (yc - y)*(yc - y));
	int dx = xe - xs;
	int dy = ye - ys;
	double slope = (double)dy / dx;
	if (abs(dy) < abs(dx)) {
		if (xs > xe) {
			swap(xs, xe);
			swap(ys, ye);
		}
		double x = xs;
		double y = ys;
		double distanceToCenter = sqrt((xc - x)*(xc - x) + (yc - y)*(yc - y));
		if (distanceToCenter <= R)
			SetPixel(hdc, x, y, color);
		while (x < xe) {
			x = x + 1;
			y = y + slope;
			int Y = round(y);
			double distanceToCenter = sqrt((xc - x)*(xc - x) + (yc - Y)*(yc - Y));
			if (distanceToCenter <= R)
				SetPixel(hdc, x, Y, color);
		}
	}
	else {
		int islope = 1;
		if (dy != 0)
			islope = dx / dy;
		if (ys > ye) {
			swap(xs, xe);
			swap(ys, ye);
		}
		double x = xs;
		double y = ys;
		double distanceToCenter = sqrt((xc - x)*(xc - x) + (yc - y)*(yc - y));
		if (distanceToCenter <= R)
			SetPixel(hdc, x, y, color);
		while (y < ye) {
			y = y + 1;
			x = x + 1 / slope;
			int Y = round(y);
			double distanceToCenter = sqrt((xc - x)*(xc - x) + (yc - Y)*(yc - Y));
			if (distanceToCenter <= R)
				SetPixel(hdc, x, Y, color);
		}
	}

}
// ========================================= Convex Filling =========================================

class Edge {
public:
	int xleft;
	int xright;

	Edge() {}
	Edge(int a, int b) {
		xleft = a;
		xright = b;
	}
};

void InitEdgeTable(Edge* table, int n) {
	for (int i = 0; i < n; i++) {
		table[i].xleft = 1000000;
		table[i].xright = 0;
	}
}

void FromEdgeToTable(Pointt p1, Pointt p2, Edge * table) {
	if (p1.y == p2.y){ 
		return; 
	}
	if (p1.y > p2.y) {
		Pointt temp = p1;
		p1 = p2;
		p2 = temp;
	}

	int y = p1.y;
	double x = p1.x;
	double mi = (p2.x - p1.x) / (double)(p2.y - p1.y);
	while (y < p2.y) {
		if (x < table[y].xleft) {
			table[y].xleft = ceil(x);
		}
		if (x > table[y].xright) {
			table[y].xright = floor(x);
		}
		x += mi;
		y++;
	}
}

void FromPolygonToTable(Pointt* p, int n, Edge* table) {
	Pointt v1 = p[n - 1];
	for (int i = 0; i < n; i++) {
		Pointt v2 = p[i];
		FromEdgeToTable(v1, v2, table);
		v1 = p[i];
	}
}

void ConvexFilling(HDC hdc, Pointt* p, int n, COLORREF color) {

	Edge table[800];

	InitEdgeTable(table, 800);

	FromPolygonToTable(p, n, table);

	for (int i = 0; i < 800; i++) {

		if (table[i].xleft <= table[i].xright) {

			// Any Line Algorithm
			DDADrawLine(hdc, table[i].xleft, i, table[i].xright, i, color);
		}
	}
}
//
// ========================================= Curve =========================================
typedef double vec[4];
typedef double matrix[4][4];
void mul(matrix A, vec B, vec C) {
	for (int i = 0; i < 4; i++) {
		C[i] = 0;
		for (int j = 0; j < 4; j++) {
			C[i] += A[i][j] * B[j];
		}
	}
}

double dot(vec a, vec b) {
	double sum = 0.0;
	for (int i = 0; i < 4; i++) {
		sum += (a[i] * b[i]);
	}
	return sum;
}

void ThirdDegreeHermit(HDC hdc, Pointt p1, Pointt t1, Pointt p2, Pointt t2, COLORREF color) {
	static matrix H = { { 2, 1, -2, 1 },{ -3, -2, 3, -1 },{ 0, 1, 0, 0 },{ 1, 0, 0, 0 } };
	vec Vx = { p1.x, t1.x, p2.x, t2.x };
	vec Vy = { p1.y, t1.y, p2.y, t2.y };
	vec Gx, Gy;
	mul(H, Vx, Gx);
	mul(H, Vy, Gy);
	int n = 1000;
	double dt = 1.0 / n;
	double t = 0;
	for (int i = 0; i < n; i++) {
		double t2 = t*t;
		double t3 = t2*t;
		vec Vt = { t3, t2, t, 1 };
		double x = dot(Gx, Vt);
		double y = dot(Gy, Vt);
		SetPixel(hdc, round(x), round(y), color);
		t += dt;
	}
}
void FirstDegree(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {
	int B1 = xs;
	int B2 = ys;
	int A1 = xe - xs;
	int A2 = ye - ys;
	int n = max(abs(A1), abs(A2));
	double dt = 1.0 / n;
	double x = (double)xs;
	double y = (double)ys;
	double t = 0;
	for (int i = 0; i < n; i++) {
		SetPixel(hdc, round(x), round(y), color);
		x = A1 *t + B1;
		y = A2 *t + B2;
		t += dt;
	}
}
void SecondDegree(HDC hdc, Pointt p1, Pointt T1, Pointt p2, COLORREF color) {
	/*int B1, B2;
	B1 = T1.x;
	B2 = T1.y;
	int A1 = p2.x - p1.x - T1.x;
	int A2 = p2.y - p1.y - T1.y;
	int G1 = p1.x;
	int G2 = p1.y;
	int n = max(abs(p2.x = p1.x), abs(p2.y - p1.y));
	double dt = 1.0 / n;
	int x = p1.x;
	int y = p1.y;
	double t = 0;
	for (int i = 0; i < n; i++) {
		SetPixel(hdc, round(x), round(y), color);
		double t2 = t*t;
		x = A1 * t2 + B1 * t + G1;
		y = A2 * t2 + B2 * t + G2;
		t += dt;
	}*/
}
void ThirdDegreeBezeir(HDC hdc, Pointt p1, Pointt p2, Pointt p3, Pointt p4, COLORREF color) {
	Pointt T1(3 * (p2.x - p1.x), 3 * (p2.y - p1.y));
	Pointt T2(3 * (p4.x - p3.x), 3 * (p4.y - p3.y));
	ThirdDegreeHermit(hdc, p1, T1, p2, T2, RGB(0, 0, 0));
}
void DrawCardinalSpline(HDC hdc, Pointt P[], int n, double c)
{
	double c1 = 1 - c;
	Pointt T0(c1*(P[2].x - P[0].x), c1*(P[2].y - P[0].y));
	for (int i = 2; i < n - 1; i++)
	{
		Pointt T1(c1*(P[i + 1].x - P[i - 1].x), c1*(P[i + 1].y - P[i - 1].y));
		ThirdDegreeHermit(hdc, P[i - 1], T0, P[i], T1, RGB(200, 0, 0));
		T0 = T1;
	}
}
/// ========================================================= Save ========================================================
void Save(vector<string> v) {
	ofstream out("SaveData.txt");
	out << v.size() << endl;
	for (int i = 0; i < v.size(); i++) {
		out << v[i] << endl;
	}
	out << color << endl;
	out.close();
}
/// ========================================================= Load ========================================================
void Load(COLORREF color) {
	HDC hdc = GetDC(hWnd);
	COLORREF c;
	vector<string> v;
	ifstream in("SaveData.txt");
	
	string s;
	in >> s;
	int size = ConvertToInteger(s);
	for (int i = 0; i < size; i++) {
		string s;
		in >> s;
		v.push_back(s);
	}
	in >> color;
	ChangeBackground(hWnd, color);
	in.close();
	for (int i = 0; i < v.size(); i++) {
		string type = GetType(v[i]);
		vector<int> shapePoints = GetPoints(v[i]);
		if (type.compare("LineMidPoint") == 0) {
			MidPointLineDraw(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("LineDDA") == 0) {
			DDADrawLine(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("LineParametric") == 0) {
			ParametricDrawLine(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("CircleCartesian") == 0) {
			CartesianDraw(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("CirclePolar") == 0) {
			PolarDraw(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("CircleIterativePolar") == 0) {
			IterativePolarDrawCircle(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("CircleMidPoint") == 0) {
			MidpointDrawCircle(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("CurveFirstDegree") == 0) {
			FirstDegree(hdc, shapePoints[1], shapePoints[2], shapePoints[3], shapePoints[4], RGB(0, 0, 255));
		}
		else if (type.compare("CurveSecondDegree") == 0) {
			Pointt p[3];
			p[0].x = shapePoints[1]; p[0].y = shapePoints[2];
			p[1].x = shapePoints[3]; p[1].y = shapePoints[4];
			p[2].x = shapePoints[5]; p[2].y = shapePoints[6];
			SecondDegree(hdc, p[0], p[1], p[2], RGB(0, 0, 255));
		}
		else if (type.compare("ConvexFilling") == 0) {
			Pointt p[5];
			p[0].x = shapePoints[1]; p[0].y = shapePoints[2];
			p[1].x = shapePoints[3]; p[1].y = shapePoints[4];
			p[2].x = shapePoints[5]; p[2].y = shapePoints[6];
			p[3].x = shapePoints[7]; p[3].y = shapePoints[8];
			p[4].x = shapePoints[9]; p[4].y = shapePoints[10];
			ConvexFilling(hdc, p, 5, RGB(0, 0, 255));
		}
		else if (type.compare("CurveThirdDegreeHermite") == 0) {
			Pointt p[4];
			p[0].x = shapePoints[1]; p[0].y = shapePoints[2];
			p[1].x = shapePoints[3]; p[1].y = shapePoints[4];
			p[2].x = shapePoints[5]; p[2].y = shapePoints[6];
			p[3].x = shapePoints[7]; p[3].y = shapePoints[8];
			ThirdDegreeHermit(hdc, p[0], p[1], p[2], p[3], RGB(0, 0, 255));
		}
		else if (type.compare("CurveSecondDegreeBezier") == 0) {
			Pointt p[4];
			p[0].x = shapePoints[1]; p[0].y = shapePoints[2];
			p[1].x = shapePoints[3]; p[1].y = shapePoints[4];
			p[2].x = shapePoints[5]; p[2].y = shapePoints[6];
			p[3].x = shapePoints[7]; p[3].y = shapePoints[8];
			ThirdDegreeBezeir(hdc, p[0], p[1], p[2], p[3], RGB(0, 0, 255));
		}
		else if (type.compare("curveSpline") == 0) {
			Pointt p[4];
			p[0].x = shapePoints[1]; p[0].y = shapePoints[2];
			p[1].x = shapePoints[3]; p[1].y = shapePoints[4];
			p[2].x = shapePoints[5]; p[2].y = shapePoints[6];
			p[3].x = shapePoints[7]; p[3].y = shapePoints[8];
			DrawCardinalSpline(hdc, p, 4, 3);
		}
	}
}

//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
	WNDCLASSEX wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style = CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc = WndProc;
	wcex.cbClsExtra = 0;
	wcex.cbWndExtra = 0;
	wcex.hInstance = hInstance;
	wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_GRAPHICSPROJECT));
	wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
	wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	wcex.lpszMenuName = MAKEINTRESOURCE(IDC_GRAPHICSPROJECT);
	wcex.lpszClassName = szWindowClass;
	wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassEx(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{

	hInst = hInstance; // Store instance handle in our global variable

	hWnd = CreateWindow(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, NULL, NULL, hInstance, NULL);

	if (!hWnd)
	{
		return FALSE;
	}

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	return TRUE;
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_COMMAND	- process the application menu
//  WM_PAINT	- Paint the main window
//  WM_DESTROY	- post a quit message and return
//
//

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	int wmId, wmEvent;
	PAINTSTRUCT ps;
	static int xs, ys, xe, ye, cnt = 0;
	static bool flag = false;
	static Pointt points[100];
	switch (message)
	{
	case WM_CREATE: {
		ShowWindow(hWnd, SW_MAXIMIZE);

		HMENU hMenubar = CreateMenu();
		HMENU lineMenu = CreateMenu();
		HMENU circleMenu = CreateMenu();
		HMENU curveMenu = CreateMenu();
		HMENU curveTypeMenu = CreateMenu();
		HMENU fillingMenu = CreateMenu();
		HMENU clippingMenu1 = CreateMenu();
		HMENU clippingMenu2 = CreateMenu();
		HMENU fileMenu = CreateMenu();
		HMENU colorMenu = CreateMenu();

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)fileMenu, L"File");
		AppendMenu(fileMenu, MF_STRING, ID_Save, L"Save");
		AppendMenu(fileMenu, MF_STRING, ID_Load, L"Load");
		AppendMenu(fileMenu, MF_STRING, ID_Exit, L"Exit");

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)colorMenu, L"Background Color");
		AppendMenu(colorMenu, MF_STRING, ID_RedColor, L"Red");
		AppendMenu(colorMenu, MF_STRING, ID_GreenColor, L"Green");
		AppendMenu(colorMenu, MF_STRING, ID_BlueColor, L"Blue");
		AppendMenu(colorMenu, MF_STRING, ID_BlackColor, L"Black");

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)lineMenu, L"Line Algorithms");
		AppendMenu(lineMenu, MF_STRING, ID_LineDDA, L"DDA");
		AppendMenu(lineMenu, MF_STRING, ID_LineMidPoint, L"Midpoint");
		AppendMenu(lineMenu, MF_STRING, ID_LineParametric, L"Parametric");

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)circleMenu, L"Circle Algorithms");
		AppendMenu(circleMenu, MF_STRING, ID_CircleCartesian, L"Cartesian");
		AppendMenu(circleMenu, MF_STRING, ID_CirclePolar, L"Polar");
		AppendMenu(circleMenu, MF_STRING, ID_CircleIterativePolar, L"Iterative Polar");
		AppendMenu(circleMenu, MF_STRING, ID_CircleMidPoint, L"Midpoint");

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)curveMenu, L"Curve Algorithms");
		AppendMenu(curveMenu, MF_STRING, ID_CurveFirstDegree, L"First Degree");
		AppendMenu(curveMenu, MF_STRING, ID_CurveSecondDegree, L"Second Degree");

		AppendMenu(curveMenu, MF_POPUP, (UINT_PTR)curveTypeMenu, L"Third Degree");
		AppendMenu(curveTypeMenu, MF_STRING, ID_CurveThirdDegreeHermite, L"Hermite");
		AppendMenu(curveTypeMenu, MF_STRING, ID_CurveThirdDegreeBezier, L"Bezier");

		AppendMenu(curveMenu, MF_STRING, ID_CurveSplines, L"Splines");

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, ID_ConvixFilling, L"Convix Filling Algorithm");

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)clippingMenu1, L"Clipping according to Rectangle Window");
		AppendMenu(clippingMenu1, MF_STRING, ID_RectangleClippingPoint, L"Point Clipping");
		AppendMenu(clippingMenu1, MF_STRING, ID_RectangleClippingLine, L"Line Clipping");

		/**********************************************************************************************************/

		AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)clippingMenu2, L"Clipping according to Circle");
		AppendMenu(clippingMenu2, MF_STRING, ID_CircleClippingPoint, L"Point Clipping");
		AppendMenu(clippingMenu2, MF_STRING, ID_CircleClippingLine, L"Line Clipping");

		/**********************************************************************************************************/

		SetMenu(hWnd, hMenubar);

		break;
	}
	case WM_COMMAND:
		wmId = LOWORD(wParam);
		wmEvent = HIWORD(wParam);
		// Parse the menu selections:
		switch (wmId)
		{
		case IDM_ABOUT:
			DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
			break;
		case ID_Exit:
			DestroyWindow(hWnd);
			break;
		case ID_Save:
			Save(shapes);
			break;
		case ID_Load:
			Load(color);
			break;
		case ID_RedColor:
			init();
			ChangeBackground(hWnd, RGB(255, 0, 0));
			color = RGB(255, 0, 0);
			break;
		case ID_GreenColor:
			init();
			ChangeBackground(hWnd, RGB(0, 255, 0));
			color = RGB(0, 255, 0);
			break;
		case ID_BlueColor:
			init();
			ChangeBackground(hWnd, RGB(0, 0, 255));
			color = RGB(0, 0, 255);
			break;
		case ID_BlackColor:
			init();
			ChangeBackground(hWnd, RGB(0, 0, 0));
			color = RGB(0, 0, 255);
			break;
		case ID_LineDDA:
			init();
			_LineDDA = _TwoPoints = true;
			break;
		case ID_LineMidPoint:
			init();
			_LineMidPoint = _TwoPoints = true;
			break;
		case ID_LineParametric:
			init();
			_LineParametric = _TwoPoints = true;
			break;
		case ID_CircleCartesian:
			init();
			_CircleCartesian = _TwoPoints = true;
			break;
		case ID_CirclePolar:
			init();
			_CirclePolar = _TwoPoints = true;
			break;
		case ID_CircleIterativePolar:
			init();
			_CircleIterativePolar = _TwoPoints = true;
			break;
		case ID_CircleMidPoint:
			init();
			_CircleMidPoint = _TwoPoints = true;
			break;
		case ID_CurveFirstDegree:
			init();
			_CurveFirstDegree = _TwoPoints = true;
			break;
		case ID_CurveSecondDegree:
			init();
			_CurveSecondDegree = true;
			break;
		case ID_CurveThirdDegreeHermite:
			init();
			_CurveThirdDegreeHermite = true;
			break;
		case ID_CurveThirdDegreeBezier:
			init();
			_CurveThirdDegreeBezier = true;
			break;
		case ID_CurveSplines:
			init();
			_CurveSplines = true;
			break;
		case ID_ConvixFilling:
			init();
			_ConvixFilling = true;
			break;
		case ID_RectangleClippingPoint:
			init();
			_RectangleClippingPoint = _OnePoint = true;
			break;
		case ID_RectangleClippingLine:
			init();
			_RectangleClippingLine =  true;
			break;
		case ID_CircleClippingPoint:
			init();
			_CircleClippingPoint = _OnePoint = true;
			break;
		case ID_CircleClippingLine:
			init();
			_CircleClippingLine = _TwoPoints = true;
			break;

		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
		}
		break;
	case WM_LBUTTONDOWN:
		if (_OnePoint == true) {
			xs = LOWORD(lParam);
			ys = HIWORD(lParam);
			hdc = GetDC(hWnd);
			if (_CircleClippingPoint == true) {
				Circle_PointClipping(hdc, xs, ys, 200, 200, 250, 250, RGB(0, 0, 255));
			}
			else if (_RectangleClippingPoint == true) {
				Rectangle_PointClipping(hdc, xs, ys, 1, 400, 1, 400, RGB(0, 0, 255));
			}
			ReleaseDC(hWnd, hdc);
		}
		else if (_TwoPoints == true) {
			if (!flag) {
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
				flag = true;
			}
			else {
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);
				hdc = GetDC(hWnd);
				//  ------------ Line ------------
				if (_LineMidPoint == true) {
					MidPointLineDraw(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "LineMidPoint";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}
				else if (_LineDDA == true) {
					DDADrawLine(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "LineDDA";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}
				else if (_LineParametric == true) {
					ParametricDrawLine(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "LineParametric";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}
				//  ------------ Circle ------------
				else if (_CircleCartesian == true) {
					CartesianDraw(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "CircleCartesian";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}
				else if (_CirclePolar == true) {
					PolarDraw(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "CirclePolar";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}
				else if (_CircleIterativePolar == true) {
					IterativePolarDrawCircle(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "CircleIterativePolar";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}
				else if (_CircleMidPoint == true) {
					MidpointDrawCircle(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "CircleMidPoint";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}
				//  ------------ Clipping ------------
				else if (_CircleClippingLine == true)
					CartesianDraw(hdc, 200, 200, 250, 250, RGB(0, 0, 255)),
					Circle_LineClipping(hdc, xs, ys, xe, ye, 200, 200, 250, 250, RGB(0, 0, 255));
				// ------------ Curve 1st degree ------------
				else if (_CurveFirstDegree == true) {
					FirstDegree(hdc, xs, ys, xe, ye, RGB(0, 0, 255));
					string str = "CurveFirstDegree";
					string s = FillString(str, xs, ys, xe, ye);
					shapes.push_back(s);
				}

				ReleaseDC(hWnd, hdc);
				flag = false;
			}
		}
		else {
			if (_ConvixFilling == true) {
				if (cnt < 5) {
					points[cnt].x = LOWORD(lParam);
					points[cnt].y = HIWORD(lParam);
					cnt++;
				}
				else {
					hdc = GetDC(hWnd);
					ConvexFilling(hdc, points, 5, RGB(150, 0, 200));
					string s = FillString2("ConvexFilling", points, 5);
					shapes.push_back(s);

					ReleaseDC(hWnd, hdc);
					cnt = 0;
				}
			}
			else if (_CurveSecondDegree == true) {
				if (cnt < 3) {
					points[cnt].x = LOWORD(lParam);
					points[cnt].y = HIWORD(lParam);
					cnt++;
				}
				else {
					hdc = GetDC(hWnd);
					SecondDegree(hdc, points[0], points[1], points[2], RGB(150, 0, 200));
					string s = FillString2("CurveSecondDegree", points, 3);
					shapes.push_back(s);

					ReleaseDC(hWnd, hdc);
					cnt = 0;
				}
			}
			else if (_CurveThirdDegreeHermite == true) {
				if (cnt < 4) {
					points[cnt].x = LOWORD(lParam);
					points[cnt].y = HIWORD(lParam);
					cnt++;
				}
				else {
					hdc = GetDC(hWnd);
					ThirdDegreeHermit(hdc, points[0], points[1], points[2], points[3], RGB(150, 0, 200));
					string s = FillString2("CurveThirdDegreeHermite", points, 4);
					shapes.push_back(s);

					ReleaseDC(hWnd, hdc);
					cnt = 0;
				}
			}
			else if (_CurveThirdDegreeBezier == true) {
				if (cnt < 4) {
					points[cnt].x = LOWORD(lParam);
					points[cnt].y = HIWORD(lParam);
					cnt++;
				}
				else {
					hdc = GetDC(hWnd);
					ThirdDegreeBezeir(hdc, points[0], points[1], points[2], points[3], RGB(150, 0, 200));
					string s = FillString2("CurveSecondDegreeBezier", points, 4);
					shapes.push_back(s);

					ReleaseDC(hWnd, hdc);
					cnt = 0;
				}
			}
			else if (_CurveSplines == true) {
				if (cnt < 4) {
					points[cnt].x = LOWORD(lParam);
					points[cnt].y = HIWORD(lParam);
					cnt++;
				}
				else {
					hdc = GetDC(hWnd);
					DrawCardinalSpline(hdc, points, 4, 3);
					string s = FillString2("curveSpline", points, 4);
					shapes.push_back(s);

					ReleaseDC(hWnd, hdc);
					cnt = 0;
				}
			}
			else if (_RectangleClippingLine == true) {
				if (!OK) {
					if (cnt < 4) {
						points[cnt].x = LOWORD(lParam);
						points[cnt].y = HIWORD(lParam);
						cnt++;
					}
					else {
						hdc = GetDC(hWnd);
						OK = true;
						DDADrawLine(hdc, points[0].x, points[1].y, points[0].x, points[0].y, RGB(255, 0, 0));
						DDADrawLine(hdc, points[0].x, points[0].y, points[1].x, points[0].y, RGB(255, 0, 0));
						DDADrawLine(hdc, points[1].x, points[0].y, points[1].x, points[1].y, RGB(255, 0, 0));
						DDADrawLine(hdc, points[0].x, points[1].y, points[1].x, points[1].y, RGB(255, 0, 0));
						CohenSuth(hdc, points[2].x, points[2].y, points[3].x, points[3].y, points[0].x, points[1].x, points[1].y, points[0].y, RGB(255, 0, 0));
						ReleaseDC(hWnd, hdc);
						cnt = 0;
					}
				}
				else {
					if (!flag) {
						xs = LOWORD(lParam);
						ys = HIWORD(lParam);
						flag = true;
					}
					else {
						xe = LOWORD(lParam);
						ye = HIWORD(lParam);
						hdc = GetDC(hWnd);
						CohenSuth(hdc, xs, ys, xe, ye, points[0].x, points[1].x, points[1].y, points[0].y, RGB(255, 0, 0));
						ReleaseDC(hWnd, hdc);
						flag = false;
					}

				}
				
			}
		}
		break;
	case WM_PAINT:
		hdc = BeginPaint(hWnd, &ps);
		// TODO: Add any drawing code here...
		EndPaint(hWnd, &ps);
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
	switch (message)
	{
	case WM_INITDIALOG:
		return (INT_PTR)TRUE;

	case WM_COMMAND:
		if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
			return (INT_PTR)TRUE;
		}
		break;
	}
	return (INT_PTR)FALSE;
}
