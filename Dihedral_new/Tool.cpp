#include "Tool.h" 

vector<pair<string, Scalar>> colorbar = {
	{ "black", Scalar(0, 0, 0) },
	{ "white", Scalar(255, 255, 255) },
	{ "blue", Scalar(255, 0, 0) },
	{ "green", Scalar(0, 255, 0) },
	{ "red", Scalar(0, 0, 255) },
	{ "red1", Scalar(127, 110, 174) },
	{ "gray", Scalar(194, 194, 194) },
	{ "blue1", Scalar(251, 228, 169) },
	{ "blue2", Scalar(251, 204, 176) },
	{ "blue3", Scalar(251, 181, 105) },
	{ "green1", Scalar(130, 174, 89) },
	{ "green2", Scalar(222, 250, 167) },
	{ "yellow", Scalar(0, 255, 255) },
	{ "orange1", Scalar(70, 124, 217) },
	{ "orange2", Scalar(3, 142, 249) },
	{ "pink", Scalar(237, 171, 245) },
	{ "lightorange2", Scalar(175, 211, 249) },
	{ "lightpink", Scalar(244, 211, 247) } };


//draw tools
void draw_poly(Mat &drawing_, vector<Point2f> contour_s, Point2f shift, int color)
{
	Scalar col_sca;
	if (color < colorbar.size())
		col_sca = colorbar[color].second;
	else col_sca = Scalar(0, 0, 0);
	int n = contour_s.size();
	//cout << "n: " << n << endl;
	Point rook_points[1][2000];
	for (int t = 0; t < n; t++)
	{
		rook_points[0][t] = contour_s[t] + shift;
	}
	const Point* ppt[1] = { rook_points[0] };
	int npt[] = { n };
	fillPoly(drawing_,
		ppt,
		npt,
		1,
		col_sca //黑色
				//Scalar(255, 255, 255) //白色
		);
	//circle(drawing_, contour_s[0] + shift, 4, Scalar(255), 3);
	//circle(drawing_, center, 4, Scalar(0, 255, 255), -1);
}


void draw_contour(Mat &drawing_, vector<Point2f> contour_s, Point2f shift, int color, int thickness)
{
	Scalar col_sca;
	if (color < colorbar.size())
		col_sca = colorbar[color].second;
	else col_sca = Scalar(0, 0, 0);
	int n = contour_s.size();
	for (int t = 0; t < n; t++)
	{
		int lineType = 8;
		line(drawing_,
			contour_s[t] + shift,
			contour_s[(t + 1) % n] + shift,
			col_sca,
			thickness,
			lineType);
	}
	//cout << "n: " << n << endl;
}

void draw_contour_points(Mat &drawing_, vector<Point2f> contour_s, Point2f shift, int color, int radius)
{
	Scalar col_sca;
	if (color < colorbar.size())
		col_sca = colorbar[color].second;
	else col_sca = Scalar(0, 0, 0);
	int n = contour_s.size();
	for (int t = 0; t < n; t++)
		circle(drawing_, contour_s[t] + shift, radius, col_sca, -1);
	//cout << "n: " << n << endl;
}


//geometry tool
Point2f center_p(vector<Point2f> contour_)
{
	//利用轮廓的矩
	Moments mu = moments(contour_);
	return Point2f(mu.m10 / mu.m00, mu.m01 / mu.m00);
}

vector<Point2f> Trans_contour(vector<Point2f> c1, Point2f trans_shift)
{
	vector<Point2f> c2;
	for (int i = 0; i < c1.size(); i++)
	{
		c2.push_back(c1[i] + trans_shift);
	}
	return c2;
}

vector<Point2f> Rotate_contour(vector<Point2f> src, Point2f center, double angle)
{
	vector<Point2f> dst;
	cv::Mat rot_mat = cv::getRotationMatrix2D(center, angle, 1.0);
	/*cv::warpAffine(src, dst, rot_mat, src.size());*/
	cv::transform(src, dst, rot_mat);
	return dst;
}
vector<Point2f> Flip_contour(vector<Point2f> cont_s)
{
	Point2f ccen = center_p(cont_s);
	int cont_size = cont_s.size();
	//flip horizontal
	for (int i = 0; i < cont_size; i++)
	{
		cont_s[i].x = 2 * ccen.x - cont_s[i].x;
	}
	for (int i = 0; i < cont_size / 2; i++)
	{
		Point2f mid = cont_s[i];
		cont_s[i] = cont_s[cont_size - 1 - i];
		cont_s[cont_size - 1 - i] = mid;
	}
	return cont_s;
}

vector<int> cal_feature(vector<Point2f> contour_, int  n_min, int n_max, double angle_cos, bool show_result)
{
	vector<int> index_num;
	double angle_back = INF; //顶部对应点的值
	double angle_start = INF; //起始点对应的值
	int contoursize = contour_.size();
	//cout << "contoursize: " << contoursize << endl;
	//double arl = arcLength(contour_, true);
	//cout<<"arl: "<<arl<<"  "<<dmin<<" "<<dmax<<endl;
	int  dmid = 0.5* n_min + 0.5*n_max;
	//cout << dmid << endl;
	FOR(i, 0, contoursize)
	{
		int f = 0;
		double angle = 0;
		for (int k = n_min; k < n_max; k++)
		{
			double length_l = length_2p(contour_[i], contour_[(i + contoursize - k) % contoursize]);
			double length_r = length_2p(contour_[i], contour_[(i + k) % contoursize]);
			double length_op = length_2p(contour_[(i + contoursize - k) % contoursize], contour_[(i + k) % contoursize]);
			double angle1 = cos_3edges(length_l, length_r, length_op);
			//cout << length_l << "   " << length_r << "   " << length_op << "    " << angle1 << endl;
			if (angle1 < angle_cos)  //角度大于angle_cos的度数
			{
				f = 1;
				break;
			}
			else
			{
				angle += angle1; //计算平均值
								 //if (angle1 > angle) angle = angle1;    //angle1越小角度越大，满足条件的所有角里度数最大的那个
			}
		}
		angle = angle / (n_max - n_min);
		//cout <<i<< "   angle: " << angle << endl;
		if (f == 0 && angle < 2)
		{
			//cout << "ok~" << endl;
			if (index_num.empty())
			{
				angle_start = angle;
				angle_back = angle;
				index_num.push_back(i);
				//cout << "input fea: " << i << endl;
			}
			else
			{
				if ((i - index_num.back())> dmid)
				{
					angle_back = angle;
					index_num.push_back(i);
				}
				else
				{
					double multicross1 = crossProduct_2v(contour_[(i + contoursize - 1) % contoursize] - contour_[i], contour_[(i + 1) % contoursize] - contour_[i]);
					int back = index_num.back();
					double multicross2 = crossProduct_2v(contour_[(back + contoursize - 1) % contoursize] - contour_[back], contour_[(back + 1) % contoursize] - contour_[back]);
					//判断是凸点还是凹点，凸点优先，如果凹凸性一致，比较角度
					if (multicross1 > 0)
					{
						if (((multicross2 > 0) && (angle > angle_back)) || multicross2 < 0)
						{
							index_num.pop_back();
							angle_back = angle;
							if (index_num.empty())
								angle_start = angle;
							index_num.push_back(i);
							//cout << "index_num.pop   input fea: " << i << endl;
						}
					}
					else
					{
						if ((multicross2 < 0) && (angle > angle_back))
						{
							index_num.pop_back();
							angle_back = angle;
							if (index_num.empty())
								angle_start = angle;
							index_num.push_back(i);
						}
					}
				}
			}
		}
	}
	//判断最后一个点和第一个点，防止两点过于临近
	if ((!index_num.empty()) && ((index_num[0] + contoursize - index_num.back()) < dmid))
	{
		if (angle_back > angle_start)
		{
			vector<int> index2;
			index2.assign(index_num.begin() + 1, index_num.end());
			index_num.swap(index2);
			//index_num[0] = index_num.back();
		}
		else index_num.pop_back();
	}
	//show feature points
	cout << "Feature points:" << index_num.size() << endl;
	if (show_result)
	{
		Mat drwa = Mat::zeros(600, 600, CV_8UC3);
		int final_c_size = contour_.size();
		circle(drwa, contour_[0], 3, Scalar(0, 0, 255), -1);
		for (int i = 0, j = 0; i < final_c_size; i++)
		{
			if (j < index_num.size() && i == index_num[j])
			{
				circle(drwa, contour_[i], 2, Scalar(0, 255, 0), -1);
				j++;
				//cout << "j: " << i<<"  "<<j << endl;
			}
			else
				circle(drwa, contour_[i], 1, Scalar(255, 0, 0), -1);
		}

		imshow("feature test", drwa);
		//for (int g = 0; g < index_num.size(); g++)
		//	cout << index_num[g] << "   " << contour_[index_num[g]] << endl;
	}
	//
	return index_num;
}

vector<Point2f> con_sample(vector<Point2f> contour_, vector<int> &feature_, int sam_num, bool show_result)
{
	double length = arcLength(contour_, true);
	double Lambda = length / sam_num;

	vector<Point2f> contour_sam;
	vector<int> feature_temp;
	int csize = contour_.size();
	int fsize = feature_.size();
	if (fsize == 0)
	{
		contour_sam = sampling_ave(contour_, sam_num);
	}
	else
	{
		int f_max = 0;
		double l_max = 0;
		int p_num = 0;
		vector<int> actual_sam_num(fsize, 0);
		for (int i = 0; i < fsize; i++)  //依次计算相邻两个feature点的距离然后重采样,重采样轮廓的第一个点必是第一个特征点
		{
			double length_ = 0;
			int p_s = feature_[i];
			int p_e = feature_[(i + 1) % fsize];
			p_e = (p_e > p_s) ? p_e : p_e + csize;
			p_num = p_e - p_s;
			for (int j = p_s; j < p_e; j++)
			{
				length_ += length_2p(contour_[j% csize], contour_[(j + 1) % csize]);
			}
			if (length_ > l_max)
			{
				l_max = length_;
				f_max = i;
			}
			//确定重采样的长度和个数
			int num = length_ / Lambda;
			int inter_num = length_ >(0.5 + num)*Lambda ? num + 1 : num;
			double length2 = length_ / inter_num;
			//cout <<"length_: "<< length_<<"   Lambda: "<< Lambda<<"   inter_num: "<<inter_num <<"   "<< length2<<"  p_num:"<< p_num<<endl;
			Point2f sample = contour_[feature_[i]];
			contour_sam.push_back(sample);
			feature_temp.push_back(contour_sam.size() - 1);
			actual_sam_num[i] = 1;
			//cout << "fea sample: " << sample << endl;
			length_ = 0;
			for (int j = p_s; j < p_e; j++)
			{
				if (length_ == 0) length_ = length_2p(sample, contour_[(j + 1) % csize]);
				else length_ += length_2p(contour_[j%csize], contour_[(j + 1) % csize]);
				if (length_ > length2)
				{
					Point2f vec = unit_vec(contour_[(j + 1) % csize] - contour_[j % csize]);
					sample = contour_[(j + 1) % csize] - (length_ - length2) * vec;
					contour_sam.push_back(sample);
					actual_sam_num[i] ++;
					//cout << "sample: " << sample << endl;
					length_ = 0;
					j--;
				}
				else if (j == (p_e - 1) && length_ < 0.5* length2)
				{
					contour_sam.pop_back();
					actual_sam_num[i] --;
				}
			}
			//cout << p_num << "  actual_sam_num:   " << actual_sam_num[i] << endl;
		}
		//cout << "feature_.size(): " << feature_temp.size() << endl;
		//for (int g = 0; g < feature_temp.size(); g++)
		//	cout << "g: "<< feature_temp[g] << endl;
		//post process
		int cs_size = contour_sam.size();
		//cout << "cs_size: " << cs_size << "  f_max:"<< f_max<<endl;
		if (cs_size != sam_num) //need process
		{
			vector<Point2f> contour_sam_post;
			vector<Point2f> contour_sam_temp;
			int p_s = feature_[f_max];
			int p_e = feature_[(f_max + 1) % fsize];
			p_e = (p_e > p_s) ? p_e : p_e + csize;
			int new_num = actual_sam_num[f_max];
			int target_num = new_num + sam_num - cs_size;
			int count = 0;
			//int add_num = sam_num - cs_size;		
			while (new_num != target_num && count < 5)
			{
				contour_sam_post.swap(vector<Point2f>());
				count++;
				double length_ = 0;
				for (int j = p_s; j < p_e; j++)
				{
					length_ += length_2p(contour_[j% csize], contour_[(j + 1) % csize]);
				}
				int inter_num = target_num > new_num ? new_num + 1 : new_num - 1;
				//cout << "new_num: " << new_num << "  target:  " << target_num << "  inter:" << inter_num << endl;
				double length2 = length_ / inter_num;
				//cout << inter_num << endl;
				Point2f sample = contour_[feature_[f_max]];
				contour_sam_post.push_back(sample);
				new_num = 1;
				//cout << "fea sample: " << sample << endl;
				length_ = 0;
				for (int j = p_s; j < p_e; j++)
				{
					if (length_ == 0) length_ = length_2p(sample, contour_[(j + 1) % csize]);
					else length_ += length_2p(contour_[j%csize], contour_[(j + 1) % csize]);
					if (length_ > length2)
					{
						Point2f vec = unit_vec(contour_[(j + 1) % csize] - contour_[j % csize]);
						sample = contour_[(j + 1) % csize] - (length_ - length2) * vec;
						contour_sam_post.push_back(sample);
						//cout << "sample: " << sample << endl;
						new_num++;
						length_ = 0;
						j--;
					}
					else if (j == (p_e - 1) && length_ < 0.5* length2)
					{
						contour_sam_post.pop_back();
						new_num--;
					}
				}
				//cout << "new_num: " << new_num << "  target:  " << target_num << "  contour_sam_post:"<< contour_sam_post.size()<<endl;
			}
			for (int t = 0; t < feature_temp[f_max]; t++)
			{
				contour_sam_temp.push_back(contour_sam[t]);
			}
			contour_sam_temp.insert(contour_sam_temp.end(), contour_sam_post.begin(), contour_sam_post.end());
			for (int t = (feature_temp[(f_max + 1) % fsize] + cs_size - 1) % cs_size + 1; t < cs_size; t++)
			{
				contour_sam_temp.push_back(contour_sam[t]);
			}
			for (int g = 0; g < fsize; g++)
			{
				if (g>f_max) feature_temp[g] += (sam_num - cs_size);
			}
			//cout << "    hahahah:    " << contour_sam_temp.size() <<"   : "<< feature_temp.size()<< endl;
			contour_sam = contour_sam_temp;
		}
	}
	//post process finish
	//cout << "feature_.size(): " << feature_temp.size() << endl;
	//for (int g = 0; g < feature_temp.size(); g++)
	//	cout << "g: " <<g<<"  "<< feature_temp[g] << endl;
	feature_ = feature_temp;
	int final_c_size = contour_sam.size();
	cout << "Final contour num: " << csize << " -> " << final_c_size << "    Feature points num: " << feature_.size() << endl;
	//show feature points
	if (show_result)
	{
		Mat drwa = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		Mat drwa2 = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		FOR(i, 0, csize) circle(drwa, contour_[i], 2, Scalar(255, 0, 0), -1);
		FOR(i, 0, final_c_size)
		{
			circle(drwa, contour_sam[i], 2, Scalar(0, 0, 255), -1);
			circle(drwa2, contour_sam[i], 2, Scalar(0, 0, 0), -1);
		}
		FOR(i, 0, feature_.size()) circle(drwa2, contour_sam[feature_[i]], 2, Scalar(200, 200, 0), -1);
		circle(drwa, contour_sam[0], 3, Scalar(0, 255, 0), -1);
		circle(drwa2, contour_sam[0], 3, Scalar(0, 255, 0), -1);
		imshow("sample test", drwa);
		imshow("sample test2", drwa2);
	}
	//
	return contour_sam;
}

vector<Point2f> sampling_ave(vector<Point2f> contour_, int points_num)  //points_num是采样点的个数
{
	double length = arcLength(contour_, true);
	double Lambda = length / points_num;

	/*int sam_num = points_num * 100;
	double Lambda = length / sam_num;*/

	vector<Point2f> contour_sam;
	Point2f sample;
	//contour_sam_index记录初步采样点是哪些点
	contour_sam.push_back(contour_[0]);
	sample = contour_[0];
	int csize = contour_.size();
	for (int t = 1; t <= csize; t++)
	{
		double length_ = length_2p(sample, contour_[t%csize]);
		if (length_ > Lambda)
		{
			Point2f vec = unit_vec(contour_[t%csize] - sample);
			sample = sample + Lambda * vec;
			contour_sam.push_back(sample);
			t = t - 1;
		}
		else if (t < csize)
		{
			while ((length_ + length_2p(contour_[t], contour_[(t + 1) % csize])) < Lambda)
			{
				length_ = length_ + length_2p(contour_[t], contour_[(t + 1) % csize]);
				t++;
				if (t > (contour_.size() - 1)) break;
			}
			if (t >(contour_.size() - 1)) break;
			Point2f vec = unit_vec(contour_[(t + 1) % csize] - contour_[t]);
			sample = contour_[t] + (Lambda - length_) * vec;
			contour_sam.push_back(sample);
		}
	}
	//if (length_2p(contour_sam[0], contour_sam[contour_sam.size() - 1]) < 1)
	//{
	//	contour_sam.pop_back();
	//}
	//cout << "contour_sam: " << contour_sam.size() << endl;
	return contour_sam;
}

//bbx
vector<Point2f> bbx(vector<Point2f> &cont)
{
	vector<Point2f> four_cor;
	double bbx_max_x = -10000;
	double bbx_max_y = -10000;
	double bbx_min_x = 10000;
	double bbx_min_y = 10000;
	for (int i = 0; i < cont.size(); i++)
	{
		if (cont[i].x < bbx_min_x) bbx_min_x = cont[i].x;
		if (cont[i].x > bbx_max_x) bbx_max_x = cont[i].x;
		if (cont[i].y < bbx_min_y) bbx_min_y = cont[i].y;
		if (cont[i].y > bbx_max_y) bbx_max_y = cont[i].y;
	}
	four_cor.push_back(Point2f(bbx_min_x, bbx_max_y));
	four_cor.push_back(Point2f(bbx_min_x, bbx_min_y));
	four_cor.push_back(Point2f(bbx_max_x, bbx_min_y));
	four_cor.push_back(Point2f(bbx_max_x, bbx_max_y));
	return four_cor;
}


//translate
vector<Point2f> PtoP2f(vector<Point> cont)
{
	vector<Point2f> c1;
	int csize = cont.size();
	for (int i = 0; i < csize; i++)
	{
		c1.push_back(cont[i]);
	}
	return c1;
}
vector<Point> P2ftoP(vector<Point2f> cont)
{
	vector<Point> c1;
	int csize = cont.size();
	for (int i = 0; i < csize; i++)
	{
		c1.push_back(cont[i]);
	}
	return c1;
}


// ---------- vector process------------
double length_2p(Point2f &u, Point2f &v)
{
	return sqrt((u.x - v.x)*(u.x - v.x) + (u.y - v.y)*(u.y - v.y));
}

double cos_2v(Point2f &v0, Point2f &v1)
{
	Point2f v0u = unit_vec(v0);
	Point2f v1u = unit_vec(v1);
	return v0u.x*v1u.x + v0u.y*v1u.y;
}

double sin_2v(Point2f &v0, Point2f &v1)
{
	Point2f v0u = unit_vec(v0);
	Point2f v1u = unit_vec(v1);
	return v0u.x*v1u.y - v0u.y*v1u.x;
}

double cos_3edges(double l1, double l2, double l3)
{
	//l1,l2分别为左右两边,l3为对边
	return (l1 * l1 + l2 * l2 - l3 * l3) / (2 * l1 * l2);
}

double  crossProduct_2v(Point2f &v0, Point2f &v1)
{
	return v0.x*v1.y - v0.y*v1.x;
}

double cross(Point2f a, Point2f b, Point2f c) //三点叉积(a-c)*(b-c)
{
	return (a.x - c.x)*(b.y - c.y) - (b.x - c.x)*(a.y - c.y);
}

Point2f unit_vec(Point2f vec)
{
	double fenmu = sqrt(vec.x*vec.x + vec.y*vec.y);
	Point2f unit = Point2f(vec.x / fenmu, vec.y / fenmu);
	return unit;
}

Point2f vertical_vec(Point2f vec)
{
	Point2f vvec = Point2f(-vec.y, vec.x);
	return unit_vec(vvec);
}

int line_intersection(Point2f start1, Point2f end1, Point2f start2, Point2f end2, Point2f &cross_p)
{
	Point2f s10 = end1 - start1;
	Point2f s32 = end2 - start2;
	Point2f s02 = start1 - start2;
	float s_numer, t_numer, denom, t;
	denom = s10.x * s32.y - s32.x * s10.y;
	s_numer = s10.x * s02.y - s10.y * s02.x;
	t_numer = s32.x * s02.y - s32.y * s02.x;

	if (denom == 0)//平行或共线
	{
		if (s_numer == 0)//Collinear,返回离end1最近的点
		{
			double dis1 = sqrt((start2.x - end1.x)*(start2.x - end1.x) + (start2.y - end1.y)*(start2.y - end1.y));
			double dis2 = sqrt((end2.x - end1.x)*(end2.x - end1.x) + (end2.y - end1.y)*(end2.y - end1.y));
			if (dis1 > dis2)
			{
				cross_p = end2;
			}
			else {
				cross_p = start2;
			}
			return 2;
		}
		else return 0; // parallel
		cout << "denom == 0" << endl;
	}
	bool denomPositive = denom > 0;
	if ((s_numer < 0) == denomPositive)//参数是大于等于0且小于等于1的，分子分母必须同号且分子小于等于分母
		return 0; // No collision

	if ((t_numer < 0) == denomPositive)
		return 0; // No collision

	if (fabs(s_numer) > fabs(denom) || fabs(t_numer) > fabs(denom))
		return 0; // No collision
				  // Collision detected
	t = t_numer / denom;
	//cout << "t:" << t << endl;
	cross_p.x = start1.x + t * s10.x;
	cross_p.y = start1.y + t * s10.y;
	return 1;
}


bool self_intersect(vector<Point2f> &contour_, int &first, int &second)
{
	int sizec = contour_.size();
	for (int i = 0; i < sizec - 2; i++)
	{
		int start_n;
		int end_n;
		if (i == 0)
		{
			start_n = 2;
			end_n = sizec - 1;
		}
		else {
			start_n = i + 2;
			end_n = sizec;
		}
		FOR(j, start_n, end_n)
		{
			Point2f crossp;
			if (line_intersection(contour_[i], contour_[i + 1], contour_[j], contour_[(j + 1) % sizec], crossp) == 1)
			{
				first = i;
				second = j;
				return true;
			}
		}
	}
	first = -1;
	second = -1;
	return false;
}
//--------geometry tool----------

//image and file tools
vector<Point2f> ima2contour(string imapath, bool show_result)
{
	vector<Point2f> final_c;
	Mat src;
	Mat src_gray;
	Mat src_2;
	int thresh = 100;
	int max_thresh = 255;

	//read image
	cout << "Read "<<imapath << endl;
	src = imread(imapath, IMREAD_COLOR);
	//src = imread(imageName, CV_LOAD_IMAGE_UNCHANGED);
	if (src.empty())
	{
		cerr << "No image supplied ..." << endl;
		return final_c;
	}
	/*这里的处理过程：
	1.将任意一张图像转化成灰度图
	2.模糊，利于下一步进行筛选过滤
	3.转化二值图
	4.模糊方便提取轮廓*/
	bool valid_contour = false;
	Mat canny_output;
	//考虑到可能有多个轮廓
	vector<vector<cv::Point> > contours;
	vector<Vec4i> hierarchy;
	int count = 0;
	cvtColor(src, src_gray, COLOR_BGR2GRAY);
	while (!valid_contour)
	{
		if (count>0)  //提升处理强度
		{
			GaussianBlur(src_gray, src_gray, Size(3, 3), 10, 10);
			//未经处理的图需要四步，之前处理好的用四步会处理过头
			threshold(src_gray, src_gray, 128, 255, cv::THRESH_BINARY); //255 white
		}
		blur(src_gray, src_gray, Size(3, 3));
		//用candy法由灰度图求出掩码图
		Canny(src_gray, canny_output, thresh, thresh * 2, 3);
		//imshow("canny_output", canny_output);
		imwrite(imapath, src_gray);
		//由掩码图求出有序轮廓点
		findContours(canny_output, contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_NONE, cv::Point(0, 0));
		cout << "contours num:" << contours.size() <<"  1st points num:" << contours[0].size() <<endl;
		if (contours[0].size() < 1700 && contours[0].size() > 100)
		{
			for (int i = 0; i < contours[0].size(); i++)
				final_c.push_back(contours[0][i]);
			valid_contour = true;
		}
		count++;
		if (count > 5) return final_c;
	}
	cout << "It is the " << count << "th extracting" << endl;
	string file_suf = "C:/Users/liuxk/OneDrive/Recent/DualNPR/dataset/";
	if (show_result)
	{
		vector<Point2f> bbb = bbx(PtoP2f(contours[0]));
		int raw = abs(bbb[0].y - bbb[2].y) + 200;
		int col = abs(bbb[0].x - bbb[2].x) + 200;

		Mat drwa = Mat::zeros(raw, col, CV_8UC3);
		int i = 0;
		int final_c_size = final_c.size();
		for (; i < final_c_size / 4; i++)
		{
			circle(drwa, final_c[i], 1, Scalar(255, 0, 0), -1);
		}
		for (; i < final_c_size / 2; i++)
		{
			circle(drwa, final_c[i], 1, Scalar(0, 255, 0), -1);
		}
		for (; i < final_c_size; i++)
		{
			circle(drwa, final_c[i], 1, Scalar(0, 0, 255), -1);
		}
		//cout << "path length" << file_suf.length() << endl;
		imshow(imapath.substr(file_suf.length()), drwa);
	}
	//cout << imapath <<endl<< file_suf<<endl<< imapath.substr(file_suf.length(), imapath.length() - file_suf.length() - 3) << endl;
	string filepath = file_suf.substr(0, file_suf.length() - 8) + "contour/" + imapath.substr(file_suf.length(), imapath.length() - file_suf.length() - 3) + "txt";
	//string filepath = imapath.substr(0, imapath.length()-3) + "txt";
	fileout(filepath, final_c);
	return final_c;
}

vector<Point2f> load_point_file(string filepath)
{
	vector<Point2f> con_point;
	//读取一个存有轮廓点的文件，格式对应上一步计算轮廓点保存的文件
	ifstream in(filepath);
	if (!in.is_open())
	{
		cout << filepath << endl;
		cout << "Error opening file" << endl;
		return con_point;
	}
	//挨个处理每个字符
	//cout << "Opening file!!!" << endl;
	vector<char> each_point;
	int p_num;
	in >> p_num;
	cout << "points num: " << p_num << endl;
	while (!in.eof())
	{
		double aa;
		in >> aa;
		char bb;
		in >> bb;
		double cc;
		in >> cc;
		if (in.eof())
			break;
		//cout << "aa: " << aa <<" "<<bb<<" "<<cc<< endl;
		con_point.push_back(Point2f(aa, cc));
	}
	in.close();
	return con_point;
}

void fileout(string filepath, vector<Point2f> contour_)
{
	ofstream out(filepath);
	if (out.is_open())
	{
		out << contour_.size() << endl;//contours[0].size()
		for (int j = 0; j < contour_.size(); j++)
			out << contour_[j].x << "," << contour_[j].y << endl;
		//out << contour_[0].x << "," << contour_[0].y << endl;  //首尾连起来
	}
	cout << "contour has : " << contour_.size() << " points" << endl;
	out.close();
}

//compare two contours by TAR
vector<vector<double>> compute_TAR(vector<Point2f> &contour_, double &shape_complexity, double frac)
{
	vector<vector<double>> all_tar;
	int consize = contour_.size();
	//cout << "consize: " << consize << endl;
	int tar_num = frac * consize - 1;
	shape_complexity = 0;
	//cout << "consize: " << consize << " tar_num: " << tar_num << endl;
	vector<double> maxtar(tar_num, 0);// 记录最大值来进行归一化
	for (int i = 0; i < consize; i++)
	{
		vector<double> one_p_tar;
		for (int j = 0; j < tar_num; j++)
		{
			Point2f vpsubts_vp = contour_[(i - j - 1 + consize) % consize] - contour_[i];
			Point2f vpplusts_vp = contour_[(i + j + 1) % consize] - contour_[i];
			double tar = 0.5 * crossProduct_2v(vpplusts_vp, vpsubts_vp);
			//cout << vpsubts_vp << "  " << vpplusts_vp << endl;
			one_p_tar.push_back(tar);
			if (abs(tar) > maxtar[j]) maxtar[j] = abs(tar);

		}
		all_tar.push_back(one_p_tar);
	}
	//cout << "what happened??" << endl;
	for (int i = 0; i < consize; i++)
	{
		double max_tar_one = 0;
		double min_tar_one = 10000;
		for (int j = 0; j < tar_num; j++)
		{
			all_tar[i][j] = all_tar[i][j] / maxtar[j];
			if (all_tar[i][j] > max_tar_one) max_tar_one = all_tar[i][j];
			if (all_tar[i][j] < min_tar_one) min_tar_one = all_tar[i][j];
		}
		shape_complexity += abs(max_tar_one - min_tar_one);
	}
	shape_complexity = shape_complexity / consize;
	//cout << all_tar[0].size() << "    shape_com: " << shape_complexity << endl;
	return all_tar;
}

double tar_length_2p(vector<double> &p1, vector<double> &p2)
{
	int Ts = p1.size();
	if (Ts != p2.size())
	{
		cout << "point number not equal" << endl;
		return 0;
	}

	double result = 0;
	for (int i = 0; i < Ts; i++)
	{
		result += abs(p1[i] - p2[i]);
	}
	//cout << "test: "<<result << endl;
	return result / Ts;
}

double tar_mismatch(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<pair<int, int>>& path, int &sec_shift, int width) //点对应匹配的筛选框宽度
{
	double dis[TAR_num][TAR_num];//两组点之间的坐标差异
	double distance[TAR_num][TAR_num];
	//int step[TAR_num][TAR_num];//记录总的步数
	int first_num = first_arr.size();
	int second_num = second_arr.size(); //first 作为y轴 ,second为x轴
	if (first_num != second_num)
	{
		std::cout << "The sampling points of two contours are not equal: " << first_num << " - " << second_num << endl;
		//return 0;
	}
	if (first_arr[0].size() != second_arr[0].size())
	{
		std::cout << "The tar num of each point is not equal: " << first_arr[0].size() << " - " << second_arr[0].size() << endl;
		return 0;
	}
	double min_mis = 10000;
	vector<pair<int, int>> path_min;
	//double distance[202][202];
	//int step[202][202];//记录总的步数
	//std::cout << "first: " << first_num << "second: " << second_num << endl;
	for (int shift = 0; shift < second_num; shift++) //将first固定，分别对齐second的起点
	{
		//int ccc = 0;
		//int ddd = 0;
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				//ddd++;
				dis[i][j] = 0;
				if (max(0, i - width) <= j && j <= min(second_num - 1, i + width))
				{
					//ccc++;
					distance[i][j] = 0;
				}
				else distance[i][j] = 100000;

			}
		}
		//std::cout << "ccc :  " << ccc <<"ddd: "<<ddd<< endl;
		//distance[0][0]记录的是对齐的第一个点
		dis[0][0] = tar_length_2p(first_arr[0], second_arr[shift]);//
		distance[0][0] = dis[0][0];

		for (int i = 1; i < first_num; i++)
		{
			if (distance[i][0] == 100000) continue;
			dis[i][0] = tar_length_2p(first_arr[i], second_arr[shift]);
			distance[i][0] = distance[i - 1][0] + dis[i][0];
		}
		for (int i = 1; i < second_num; i++)
		{
			if (distance[0][i] == 100000) continue;
			dis[0][i] = tar_length_2p(first_arr[0], second_arr[(i + shift) % second_num]);
			distance[0][i] = distance[0][i - 1] + dis[0][i];
		}

		for (int i = 1; i < first_num; i++)
		{
			for (int j = 1; j < second_num; j++)
				//(int i = istart; i <= imax; i++)
			{
				if (distance[i][j] == 100000) continue;
				dis[i][j] = tar_length_2p(first_arr[i], second_arr[(j + shift) % second_num]);
				double g1 = distance[i - 1][j] + dis[i][j];
				double g2 = distance[i - 1][j - 1] + dis[i][j];
				double g3 = distance[i][j - 1] + dis[i][j];
				if (g1 < g2)
				{
					if (g1 < g3) distance[i][j] = g1;
					else distance[i][j] = g3;
				}
				else
				{
					if (g2 < g3) distance[i][j] = g2;
					else distance[i][j] = g3;
				}
			}
		}
		if (distance[first_num - 1][second_num - 1] < min_mis)
		{
			path_min.swap(vector<pair<int, int>>());
			min_mis = distance[first_num - 1][second_num - 1];
			print_TAR_Path(dis, distance, first_num - 1, second_num - 1, path_min);
			sec_shift = shift;
		}
	}
	path = path_min;
	return min_mis;
}

void print_TAR_Path(double d[][TAR_num], double dp[][TAR_num], int i, int j, vector<pair<int, int>>& path)
{
	if (i == 0 && j == 0) {
		//std::cout << first_arr[i] << " - " << second_arr[j] << endl;
		path.push_back(make_pair(i, j));
		return;
	}

	if (abs(dp[i][j] - (dp[i - 1][j - 1] + d[i][j])) < 0.001) {
		print_TAR_Path(d, dp, i - 1, j - 1, path);

	}
	else if (abs(dp[i][j] - (dp[i][j - 1] + d[i][j])) < 0.001) {
		print_TAR_Path(d, dp, i, j - 1, path);

	}
	else {
		print_TAR_Path(d, dp, i - 1, j, path);
	}
	path.push_back(make_pair(i, j));
}



double isoperimetric_inequality(vector<Point2f> contour)
{
	double arcl = arcLength(contour, true);
	double ratio = 4 * PI * contourArea(contour) / (arcl * arcl);
	return ratio;
}



