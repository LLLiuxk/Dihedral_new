#pragma once
#include <opencv2/opencv.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <highgui/highgui.hpp>
#include <imgproc/imgproc.hpp>
#include <ctime>
#include <io.h>
#include <iostream>
#include <string>
#include <set>
#include <fstream>
#include <Eigen/Dense>


#define INF 1e20
#define eps 1e-8
#define PI  acos(-1.0)
#define TAR_num 202
#define OP Point2f(0,0)
#define FOR(i,s,t)  for(int i=(s); i<(t); i++)



using namespace std;
using namespace cv;
using namespace Eigen;

extern vector<string> frame_type;
extern vector<pair<string, Scalar>> colorbar;

//draw tools
void draw_poly(Mat &drawing_, vector<Point2f> contour_s, Point2f shift, int color = 0);
void draw_contour(Mat &drawing_, vector<Point2f> contour_s, Point2f shift, int color = 0, int thickness = 1);
void draw_contour_points(Mat &drawing_, vector<Point2f> contour_s, Point2f shift, int color = 0, int radius = 1);
void draw_pair(Mat &drawing_, vector<Point2f> contour1, vector<Point2f> contour2, vector<pair<int, int>> p, Point2f shift, int color = 0, int thickness = 1);
void draw_overlap(Mat &drawing_, vector<Point2f> contour1, vector<Point2f> contour2, Point2f shift);
void draw_path(vector<Point2f> contour1, vector<Point2f> contour2, vector<int> fea1, vector<int> fea2, vector<pair<int, int>> path_fea);
void progress_bar(double index, double total);

//geometry tools
//---------polygon process-----------
Point2f center_p(vector<Point2f> contour_);
double area_polygon(vector<Point2f> poly);
vector<Point2f> Trans_contour(vector<Point2f> c1, Point2f trans_shift);
vector<Point2f> Rotate_contour(vector<Point2f> src, Point2f center, double angle);
vector<Point2f> Flip_contour(vector<Point2f> cont_s);
vector<Point2f> Rotate_contour(vector<Point2f> src, Point2f center, double angle);
Point2f rotate_p(Point2f cen, Point2f rop, double angle);
double conotour_align(vector<Point2f>& cont1, vector<Point2f>& cont2, vector<pair<int, int>> path_min);

vector<int>  cal_feature(vector<Point2f> contour_, int  n_min, int n_max, double angle_cos, bool show_result = false); //n_min, n_max 表示计算cos度数选取的点个数
vector<Point2f> con_sample(vector<Point2f> contour_, vector<int> &feature_, int sam_num, bool show_result = false);
vector<Point2f> con_resample(vector<Point2f> contour_, int sam_num, bool show_result = false);
vector<Point2f> sampling_ave(vector<Point2f> contour_, int points_num);  //points_num是采样点的个数
vector<Point2f> sampling_seg(vector<Point2f> &segment, int points_num);
vector<int> relocate(vector<Point2f> anc_points, vector<Point2f> c2);

//bbx
vector<Point2f> bbx(vector<Point2f> &cont);
bool bbx_intersection(vector<Point2f> bbx1, vector<Point2f> bbx2, Point2f &bbx_ins);

//translate
vector<Point2f> PtoP2f(vector<Point> cont);
vector<Point> P2ftoP(vector<Point2f> cont);

// ---------- vector process------------
double length_2p(Point2f &u, Point2f &v);
double cos_2v(Point2f &v0, Point2f &v1);
double sin_2v(Point2f &v0, Point2f &v1);
double angle_2v(Point2f &v0, Point2f &v1);
double cos_3edges(double l1, double l2, double l3);
double crossProduct_2v(Point2f &v0, Point2f &v1);
double cross(Point2f a, Point2f b, Point2f c); //三点叉积(a-c)*(b-c)
Point2f unit_vec(Point2f vec);
Point2f vertical_vec(Point2f vec);
vector<Point2f> base_frame(vector<Point2f> frame, int type);  //two fixed point
vector<Point2f> base_frame2(vector<Point2f> frame, int type); //one fixed point
Point2f Polar_Car(Point2f origin, Point2f axis_p, double angle, double length);

//----------cross points------------
Point2f intersection(Point2f a, Point2f b, Point2f c, Point2f d); //两条直线的交点
int line_intersection(Point2f start1, Point2f end1, Point2f start2, Point2f end2, Point2f &cross_p);
vector<Point2f> line_polygon(Point2f start1, Point2f end1, vector<Point2f> contour, bool closed = true);
vector<Point2f> poly_poly(vector<Point2f> contour, vector<Point2f> contour_);
bool self_intersect(vector<Point2f> &contour_, int &first, int &second);
bool coll_detec(vector<Point2f> contour1, vector<Point2f> contour2, int threshold);
//--------geometry tool----------


//image and file tools
vector<Point2f> ima2contour(string imapath, bool show_result = true);
vector<Point2f> load_point_file(string filepath);
vector<vector<Point2f>> load_texture(string filepath);
void fileout(string filepath, vector<Point2f> contour_);
void save_svg(string svg_path, vector<Point2f> contour, Scalar color, Point2f shift, double zoom_scale);
void write_avi(vector<Mat> images, string filename, double fps);
void write_obj(string filepath, MatrixXd V, MatrixXi F);
void contour2obj(string savepath, vector<Point2f> contour, Point2f shift, double height);
void write_para(string filepath, vector<int> indexs, vector<Point2f> new_places);
void write_twoCon(string filepath, vector<int> in1, vector<Point2f> c1, vector<int> in2, vector<Point2f> c2);

//compare two contours by TAR
vector<vector<double>> compute_TAR(vector<Point2f> &contour_, double &shape_complexity, double frac = 0.5);
double tar_length_2p(vector<double> &p1, vector<double> &p2);
double tar_mismatch(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<pair<int, int>>& path, int &sec_shift, int width); //点对应匹配的筛选框宽度
double tar_mismatch_fea(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<int> first_fea, vector<int> sec_fea, vector<pair<int, int>>& path, int &sec_shift, int width); //增加特征点的权重
void print_TAR_Path(double d[][TAR_num], double dp[][TAR_num], int i, int j, vector<pair<int, int>>& path);
double isoperimetric_inequality(vector<Point2f> contour);
double getCurvature(Point2f p1, Point2f p2, Point2f p3);

//edge evaluation and optimization
void bound_recover(vector<Point2f>& old_edge, vector<Point2f>& new_edge);

double bound_collision(vector<Point2f> cont, vector<int> indexes, int type = -1);

double edge_nd_degree(vector<Point2f> edge, int type); // non-developable degree
double edge_nd_opt(vector<Point2f>& edge, int type);
double whole_con_opt(vector<Point2f> &cont, vector<int> &indexes, int type);

void contour_de_crossing(vector<Point2f> &contour_);  //交叉的优化
void contour_fine_tuning(vector<Point2f> &contour_);  //过近的优化

//contour triangulation
vector<Point2f> triangulateContour(vector<Point2f>& con_ori, MatrixXd& V, MatrixXi& F, double ratio = 0.08);
vector<Point2f> triangulate_2Contours(vector<Point2f>& cont1, vector<Point2f>& cont2, MatrixXd& V, MatrixXi& F);
vector<Point2f> triangulate_Contours_bbx(vector<Point2f>& cont1, vector<int> anc1, MatrixXd& V, MatrixXi& F);
vector<Point2f> triangulate_bbx(vector<Point2f>& cont1, MatrixXd& V, MatrixXi& F);
int add_points(vector<Point2f>& contour, double sparse_ratio);
int point_locate(vector<Point2f> con, Point2f p);

//calculate the area of shape
inline int dcmp(double x)
{
	if (x > eps) return 1;
	return x < -eps ? -1 : 0;
}
double CPIA(vector<Point2f> a, vector<Point2f> b, int na, int nb);//ConvexPolygonIntersectArea
double SPIA(vector<Point2f> a, vector<Point2f> b);//SimplePolygonIntersectArea 调用此函数
double evaluation_area(vector<Point2f> c1, vector<Point2f> c2);
double evaluation_area_pixels(vector<Point2f> c1, vector<Point2f> c2);

//calculate the Poisson’s ratios
double cal_Poisson_ratio(vector<Point2f> con_ori, vector<Point2f> con_new);
void drawLineGraph(vector<double>& x, vector<double>& y, string save_path);
vector<Point2f> mergedVector(vector<vector<Point2f>> vecs, vector<Point2f> shifts);
string d2str(double value, int tail_num);

template<typename T>

T delete_vector(vector<T> &vec, int index_p)
{
	vector<T> vec1;
	int vsize = vec.size();
	//cout << vsize << "  " << index_p << endl;
	for (int i = 0; i < index_p; i++)
	{
		vec1.push_back(vec[i]);
	}
	T delete_p = vec[index_p];
	for (int i = index_p + 1; i <vsize; i++)
	{
		vec1.push_back(vec[i]);
	}
	vec = vec1;
	return delete_p;
}


