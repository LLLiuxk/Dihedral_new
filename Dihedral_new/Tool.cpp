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

void draw_pair(Mat &drawing_, vector<Point2f> contour1, vector<Point2f> contour2, vector<pair<int, int>> p, Point2f shift, int color, int thickness)
{
	draw_contour_points(drawing_, contour1, shift, 5, 2);
	draw_contour_points(drawing_, contour2, shift, 7, 2);
	FOR(i, 0, p.size())
	{
		line(drawing_, contour1[p[i].first] + shift, contour2[p[i].second] + shift, colorbar[6].second);
	}
}

void progress_bar(double index, double total)
{
	cout << "Loading: ";
	int show_num = index * 20 / total;
	for (int j = 0; j <= show_num; j++) cout << "█";
	cout << "  " << fixed << setprecision(2) << index*100.0 / (total - 1) << "%" << endl;
	if (index ==total - 1)
	{
		cout << endl << "======LOADING OVER======" << endl;
	}
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

double conotour_align(vector<Point2f>& cont1, vector<Point2f>& cont2, vector<pair<int, int>> path_min)
{
	double match_error = 0;
	int method = 1;
	int csize = cont1.size();
	if (cont2.size() != csize)
	{
		cout << "cannot align two contours with different size!" << endl;
		return match_error;
	}
	Point2f cen1 = center_p(cont1);
	Point2f cen2 = center_p(cont2);
	Point2f shift2 = cen1 - cen2;
	for (int i = 0; i < csize; i++)
	{
		cont2[i] = cont2[i] + shift2;
	}
	//two types of calculating the scale: i.based on the arclength  ii.based in uniform scale
	double scale = arcLength(cont1, true) / arcLength(cont2, true);
	std::cout << "scale: " << scale << endl;
	double s_1 = 0, s_2 = 0;
	FOR(i, 0, csize)
	{
		s_1 += pow(length_2p(cont1[i], cen1), 2);
		s_2 += pow(length_2p(cont2[2], cen1), 2);
	}
	double scale2 = sqrt(s_1 / s_2);
	cout << "scale2: " << scale2 << endl;
	if (method == 1)  
		//scale = scale2;
	for (int i = 0; i < cont2.size(); i++)
	{
		cont2[i] = cont2[i] * scale;
	}

	double length_min = 1000000;
	int angle_min = 0;
	int times = 0;
	double angl_al = 0;
	Point2f sh_al = Point2f(0, 0);
	Point2f shift_t = Point2f(1000, 1000);
	vector<Point2f> contour_tem;
	while (times < 3 || length_2p(shift_t, Point2f(0, 0)) > 10 || angle_min > 10)// 
	{
		length_min = 1000000;
		angle_min = 0;
		cen2 = center_p(cont2);
		//找到距离最近的角度
		for (int angle = 0; angle < 360; angle = angle + 2)
		{
			double leng = 0;
			Mat rot_mat(2, 3, CV_32FC1);
			rot_mat = getRotationMatrix2D(cen2, angle, 1);
			cv::transform(cont2, contour_tem, rot_mat);
			for (int m = 0; m < path_min.size(); m++)
			{
				leng += pow(length_2p(cont1[path_min[m].first], contour_tem[path_min[m].second]), 2);
				/*leng += length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]);*/
			}
			leng = sqrt(leng);
			if (leng < length_min)
			{
				angle_min = angle;
				length_min = leng;
			}
		}
		//cout << "angle_min: " << angle_min << "   length_min:  " << length_min << endl;
		Mat rot_mat1(2, 3, CV_32FC1);
		rot_mat1 = getRotationMatrix2D(cen2, angle_min, 1);
		cv::transform(cont2, contour_tem, rot_mat1);
		angl_al += angle_min;
		cont2 = contour_tem;
		//移动重心
		shift_t = Point2f(0, 0);
		for (int m = 0; m < path_min.size(); m++)
		{
			shift_t += cont1[path_min[m].first] - cont2[path_min[m].second];
		}
		shift_t = Point2f(shift_t.x / path_min.size(), shift_t.y / path_min.size());
		for (int i = 0; i < cont2.size(); i++)
		{
			cont2[i] = cont2[i] + shift_t;
		}
		sh_al += shift_t;
		times++;
		cout << "rotate angle_all: " << angl_al << "   length_min:  " << length_min << "    shift all: " << sh_al << "   shift_t: " << shift_t << endl;
	}
	return length_min;
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
	int  dmid = n_max;// 0.5* n_min + 0.5*n_max;
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
		Point2f shift_ = Point2f(300, 300) - center_p(contour_);
		Mat drwa = Mat::zeros(600, 600, CV_8UC3);
		int final_c_size = contour_.size();
		circle(drwa, contour_[0]+ shift_, 3, Scalar(0, 0, 255), -1);
		for (int i = 0, j = 0; i < final_c_size; i++)
		{
			int insize = index_num.size();
			if (j <insize && i == index_num[j])
			{
				circle(drwa, contour_[i]+ shift_, 2, Scalar(0, 255, 0), -1);
				j++;
				//cout << "j: " << i<<"  "<<j << endl;
			}
			else
				circle(drwa, contour_[i]+ shift_, 1, Scalar(255, 0, 0), -1);
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
		Point2f shift_ = Point2f(300, 300) - center_p(contour_sam);
		Mat drwa = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		Mat drwa2 = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		FOR(i, 0, csize) circle(drwa, contour_[i]+ shift_, 2, Scalar(255, 0, 0), -1);
		FOR(i, 0, final_c_size)
		{
			circle(drwa, contour_sam[i]+ shift_, 2, Scalar(0, 0, 255), -1);
			circle(drwa2, contour_sam[i]+ shift_, 2, Scalar(0, 0, 0), -1);
		}
		FOR(i, 0, feature_.size()) circle(drwa2, contour_sam[feature_[i]]+ shift_, 2, Scalar(200, 200, 0), -1);
		circle(drwa, contour_sam[0]+ shift_, 3, Scalar(0, 255, 0), -1);
		circle(drwa2, contour_sam[0]+ shift_, 3, Scalar(0, 255, 0), -1);
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

vector<Point2f> sampling_seg(vector<Point2f> &segment, int points_num)
{
	if (points_num <= 2) return segment;
	double length = 0; //contour_length(segment);
	for (int i = 0; i < segment.size() - 1; i++)
	{
		length += length_2p(segment[i], segment[i + 1]);
	}
	//cout << length << endl;
	double Lambda = length / (points_num - 1);

	/*int sam_num = points_num * 100;
	double Lambda = length / sam_num;*/

	vector<Point2f> contour_sam;
	Point2f sample;
	//contour_sam_index记录初步采样点是哪些点
	//contour_sam.push_back(contour_[0]);

	sample = segment[0];
	contour_sam.push_back(sample);
	int csize = segment.size();
	for (int t = 1; t <= csize; t++)
	{
		if (contour_sam.size() == points_num - 1) break;
		double length_ = length_2p(sample, segment[t%csize]);
		if (length_ > Lambda)
		{
			Point2f vec = unit_vec(segment[t%csize] - sample);
			sample = sample + Lambda * vec;
			contour_sam.push_back(sample);
			t = t - 1;
		}
		else
		{
			while ((length_ + length_2p(segment[t], segment[(t + 1) % csize])) < Lambda)
			{
				length_ = length_ + length_2p(segment[t], segment[(t + 1) % csize]);
				t++;
			}
			Point2f vec = unit_vec(segment[(t + 1) % csize] - segment[t]);
			sample = segment[t] + (Lambda - length_) * vec;
			contour_sam.push_back(sample);
		}
	}
	/*if (length_two_point2f(contour_sam[0], contour_sam[contour_sam.size() - 1]) < 1)
	{
	contour_sam.pop_back();
	}*/
	contour_sam.push_back(segment.back());
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

bool bbx_intersection(vector<Point2f> bbx1, vector<Point2f> bbx2, Point2f &bbx_ins)
{
	double left = max(bbx1[0].x, bbx2[0].x);
	double right = min(bbx1[2].x, bbx2[2].x);
	double bottom = max(bbx1[2].y, bbx2[2].y);
	double top = min(bbx1[0].y, bbx2[0].y);
	if (left < right && top < bottom) {
		bbx_ins = Point2f(0.5*(left + right), 0.5*(top + bottom));
		return 1;
	}
	else {
		return 0; //相交区域为空
	}
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

vector<Point2f> base_frame(vector<Point2f> frame, int type)  //two fixed point
{
	vector<Point2f> base_points = frame;
	//vector<double> four_edges = { length_2p(frame[0],frame[1]), length_2p(frame[1],frame[2]), length_2p(frame[2],frame[3]), length_2p(frame[3],frame[0]) };
	Point2f cen1 = 0.5* frame[0] + 0.5*frame[2];
	Point2f cen2 = 0.5* frame[1] + 0.5*frame[3];
	Point2f shift_ = cen1 - cen2;
	frame[1] += shift_;
	frame[3] += shift_;
	// parallelogram to the end
	if (type == 1) //rhombus
	{
		double length1 = 0.5*length_2p(frame[1], frame[3]);
		Point2f vec_v = vertical_vec(frame[2] - frame[0]);
		Point2f p1 = cen1 + length1*vec_v;
		Point2f p2 = cen1 - length1*vec_v;
		if (cos_2v(p1 - cen1, frame[1] - cen1) > 0)
		{
			frame[1] = p1;
			frame[3] = p2;
		}
		else
		{
			frame[1] = p2;
			frame[3] = p1;
		}
	}
	else if (type == 2) //rectangle 
	{
		double length1 = 0.5*length_2p(frame[0], frame[2]);
		frame[1] = cen1 + length1*unit_vec(frame[1] - cen1);
		frame[3] = cen1 + length1*unit_vec(frame[3] - cen1);
	}
	else if (type == 3) //square
	{
		double length1 = 0.5*length_2p(frame[1], frame[3]);
		Point2f vec_v = vertical_vec(frame[2] - frame[0]);
		Point2f p1 = cen1 + length1*vec_v;
		Point2f p2 = cen1 - length1*vec_v;
		if (cos_2v(p1 - cen1, frame[1] - cen1) > 0)
		{
			frame[1] = p1;
			frame[3] = p2;
		}
		else
		{
			frame[1] = p2;
			frame[3] = p1;
		}
		length1 = 0.5*length_2p(frame[0], frame[2]);
		frame[1] = cen1 + length1*unit_vec(frame[1] - cen1);
		frame[3] = cen1 + length1*unit_vec(frame[3] - cen1);
	}
	//for (int i = 0; i < four_edges.size(); i++)
	//	cout << "six_edges:  " << four_edges[i] << endl;
	//for (int i = 0; i < four_angles.size(); i++) 
	//	cout << "four_angles: " << four_angles[i] << endl;
	bool show = true;
	if (show)
	{
		Point2f shift = Point2f(300, 300) - center_p(base_points);
		Mat draw = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		for (int i = 0; i < base_points.size(); i++)
		{
			//cout << "base_points: " << base_points[i] << endl;
			circle(draw, base_points[i] + shift, 3, Scalar(0, 255, 0), -1);
			line(draw, base_points[i] + shift, base_points[(i + 1) % 4] + shift, Scalar(200, 200, 200), 2);
			circle(draw, frame[i] + shift, 3, Scalar(0, 0, 255), -1);
			line(draw, frame[i] + shift, frame[(i + 1) % 4] + shift, Scalar(0, 0, 0), 2);
		}
		imshow("frame: " + frame_type[type], draw);
	}
	/*cout << "length: " << length_2p(base_points[0], base_points[1]) << "   " << length_2p(base_points[1], base_points[2])
	<< "   " << length_2p(base_points[2], base_points[3]) << "   " << length_2p(base_points[3], base_points[0]) << endl;
	cout << "angle: " << 180 / PI * acos(cos_2v(base_points[3] - base_points[0], base_points[1] - base_points[0])) << "   " << 180 / PI * acos(cos_2v(base_points[0] - base_points[1], base_points[2] - base_points[1])) << endl;*/
	return frame;
}

vector<Point2f> base_frame2(vector<Point2f> frame, int type) //one fixed point
{
	vector<Point2f> base_points;
	vector<double> four_edges = { length_2p(frame[0],frame[1]), length_2p(frame[1],frame[2]), length_2p(frame[2],frame[3]), length_2p(frame[3],frame[0]) };
	vector<double> four_angles = { acos(cos_2v(frame[3] - frame[0], frame[1] - frame[0])), acos(cos_2v(frame[0] - frame[1], frame[2] - frame[1])),
		acos(cos_2v(frame[1] - frame[2], frame[3] - frame[2])), acos(cos_2v(frame[2] - frame[3], frame[0] - frame[3])) };
	if (type == 0) // parallelogram
	{
		double length1 = (four_edges[0] + four_edges[2]) / 2;
		double length2 = (four_edges[1] + four_edges[3]) / 2;
		double angle1 = (four_angles[0] + four_angles[2]) / 2;
		double angle2 = (four_angles[1] + four_angles[3]) / 2;
		base_points.push_back(frame[0]);
		base_points.push_back(frame[0] + length1*unit_vec(frame[1] - frame[0]));
		base_points.push_back(Polar_Car(base_points[1], base_points[0], angle2, length2));
		base_points.push_back(Polar_Car(base_points[2], base_points[1], angle1, length1));
	}
	else if (type == 1) //rhombus
	{
		double length1 = (four_edges[0] + four_edges[1] + four_edges[2] + four_edges[3]) / 4;
		double angle1 = (four_angles[0] + four_angles[2]) / 2;
		double angle2 = (four_angles[1] + four_angles[3]) / 2;
		base_points.push_back(frame[0]);
		base_points.push_back(frame[0] + length1*unit_vec(frame[1] - frame[0]));
		base_points.push_back(Polar_Car(base_points[1], base_points[0], angle2, length1));
		base_points.push_back(Polar_Car(base_points[2], base_points[1], angle1, length1));
	}
	else if (type == 2) //rectangle
	{
		double length1 = (four_edges[0] + four_edges[2]) / 2;
		double length2 = (four_edges[1] + four_edges[3]) / 2;
		double angle1 = acos(0);
		base_points.push_back(frame[0]);
		base_points.push_back(frame[0] + length1*unit_vec(frame[1] - frame[0]));
		base_points.push_back(Polar_Car(base_points[1], base_points[0], angle1, length2));
		base_points.push_back(Polar_Car(base_points[2], base_points[1], angle1, length1));
	}
	else if (type == 3) //square
	{
		double length1 = (four_edges[0] + four_edges[1] + four_edges[2] + four_edges[3]) / 4;
		double angle1 = acos(0);
		base_points.push_back(frame[0]);
		base_points.push_back(frame[0] + length1*unit_vec(frame[1] - frame[0]));
		base_points.push_back(Polar_Car(base_points[1], base_points[0], angle1, length1));
		base_points.push_back(Polar_Car(base_points[2], base_points[1], angle1, length1));
	}
	//for (int i = 0; i < four_edges.size(); i++)
	//	cout << "six_edges:  " << four_edges[i] << endl;
	//for (int i = 0; i < four_angles.size(); i++) 
	//	cout << "four_angles: " << four_angles[i] << endl;
	bool show = true;
	if (show)
	{
		Point2f shift = Point2f(300, 300) - center_p(base_points);
		Mat draw = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		for (int i = 0; i < base_points.size(); i++)
		{
			//cout << "base_points: " << base_points[i] << endl;
			circle(draw, base_points[i] + shift, 3, Scalar(0, 255, 0), -1);
			line(draw, base_points[i] + shift, base_points[(i + 1) % 4] + shift, Scalar(200, 200, 200), 2);
			circle(draw, frame[i] + shift, 3, Scalar(0, 0, 255), -1);
			line(draw, frame[i] + shift, frame[(i + 1) % 4] + shift, Scalar(0, 0, 0), 2);
		}
		imshow("frame: " + frame_type[type], draw);
	}
	/*cout << "length: " << length_2p(base_points[0], base_points[1]) << "   " << length_2p(base_points[1], base_points[2])
	<< "   " << length_2p(base_points[2], base_points[3]) << "   " << length_2p(base_points[3], base_points[0]) << endl;
	cout << "angle: " << 180 / PI * acos(cos_2v(base_points[3] - base_points[0], base_points[1] - base_points[0])) << "   " << 180 / PI * acos(cos_2v(base_points[0] - base_points[1], base_points[2] - base_points[1])) << endl;*/
	return base_points;
}

Point2f Polar_Car(Point2f origin, Point2f axis_p, double angle, double length)
{
	double angle_ori;
	Point2f axis = axis_p - origin;
	double angle_o_cos = cos_2v(Point2f(1, 0), axis);
	double angle_o_sin = sin_2v(Point2f(1, 0), axis);
	//cout << "axis:  "<<axis << "   " << angle_o_cos << "   " << angle_o_sin << endl;
	if (angle_o_sin > 0) angle_ori = acos(angle_o_cos);
	else angle_ori = -acos(angle_o_cos);
	//cout << angle_ori << endl;
	Point2f target = Point2f(cos(angle_ori + angle), sin(angle_ori + angle));
	target = length*  target + origin;
	return target;
}

//----------cross points------------
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

vector<Point2f> 
line_polygon(Point2f start1, Point2f end1, vector<Point2f> contour, bool closed)
{
	vector<Point2f> all_inter;
	if (closed) contour.push_back(contour[0]);
	int contsize = contour.size();
	for (int i = 0; i < contsize - 1; i++)
	{
		Point2f cen;
		int f = line_intersection(start1, end1, contour[i], contour[i + 1], cen);
		if (f == 1)
		{
			if (all_inter.empty() || length_2p(all_inter.back(), cen)>0.01)
			{
				all_inter.push_back(cen);
				//cout << "f==1" << cen << "  " << start1<<"  "<< end1 << "  " << contour[i] << "  " << contour[i + 1] << "  " << endl;
			}
		}
		else if (f == 2)
		{
			vector<Point2f> bbx1 = bbx(vector<Point2f>{ start1, end1 });
			vector<Point2f> bbx2 = bbx(vector<Point2f>{ contour[i], contour[i + 1] });
			if (bbx_intersection(bbx1, bbx2, cen))
			{
				all_inter.push_back(cen);
			}
			//all_inter.push_back(0.5 * (contour[i] + contour[i + 1]));
			//cout << "f==2" << cen << "  " << start1 << "  " << end1 << "  " << contour[i] << "  " << contour[i + 1] << "  " << endl; 
		}
	}
	//cout << "intersize :  "<<all_inter.size() << all_inter [0]<<" "<<all_inter[1]<< endl;
	return all_inter;
}

vector<Point2f> poly_poly(vector<Point2f> contour, vector<Point2f> contour_)
{
	vector<Point2f> all_inter;
	int contsize = contour.size();
	int consize2 = contour_.size();
	for (int i = 0; i < contsize; i++)
	{
		Point2f cen;
		vector<Point2f> interp = line_polygon(contour[i], contour[(i + 1) % contsize], contour_);
		if (!interp.empty())
		{
			cout << "interp: "<<interp << "   " << i << "   " << contour[i] << "   " << contour[(i + 1) % contsize] << endl;
			for (int j = 0; j < interp.size(); j++)
			{
				int f = 0;
				for (int t = 0; t < all_inter.size(); t++)
				{
					if (length_2p(all_inter[t], interp[j]) < 0.01)
						f = 1;
				}
				if (f == 0)
				{
					all_inter.push_back(interp[j]);
				}
			}
		}
	}
	return all_inter;
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

bool coll_detec(vector<Point2f> contour1, vector<Point2f> contour2, int threshold)
{
	//Mat drawing1 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
	//draw_contour_points(drawing1, contour1, Point2f(300,300), 4, 2);
	//draw_contour_points(drawing1, contour2, Point2f(300, 300), 5, 2);
	vector<Point2f> all_colli_points = poly_poly(contour1, contour2);
	//cout << "all_colli_points.size: " << all_colli_points.size()<<endl;
	//cout << contour1.size() << "   " << contour2.size() << endl;
	//for (auto p : all_colli_points) circle(drawing1, p+ Point2f(300, 300), 3, Scalar(0, 0, 255), -1);//cout << p << endl;
	//imshow("collision",drawing1);
	if (all_colli_points.size() > threshold) return true;
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
		//imwrite(imapath, src_gray);
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

//将tile保存为svg格式
void save_svg(string svg_path, vector<Point2f> contour, Scalar color, Point2f shift, double zoom_scale)
{
	//write head
	if (_access(svg_path.c_str(), 0) == -1) return;
	ofstream outfile(svg_path);
	outfile << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl
		<< "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl
		<< "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl << endl;
	outfile << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;

	//write body
	outfile << " <polygon points=\" ";
	//cout << " polygon points: " << contour[0] << endl;
	for (int i = 0; i < contour.size(); i++)
	{
		Point2f one = zoom_scale * contour[i] + shift;
		outfile << one.x << "," << one.y << " ";
	}
	outfile << "\"" << endl << "style=\"fill:rgb(" << to_string(int(color[2])) << "," << to_string(int(color[1])) << "," << to_string(int(color[0])) << ")\"/>" << endl;
	//write tail
	outfile << "</svg>" << endl;
	outfile.close();
}

void write_obj(string filepath, MatrixXd V, MatrixXi F)
{
	ofstream outfile(filepath, ios::out);
	//注意在生成.obj文件的时候，面的信息是由顶点下标+1组成的，是从1开始的！！！并不是由0开始！！！
	if (!outfile.is_open())
	{
		cerr << "open error";
		exit(1);
	}

	outfile << "#List of geometric vertices, with (x,y,z) coordinates" << endl;
	if (V.cols() == 2)
	{
		for (int i = 0; i < V.rows(); ++i) {
			outfile << "v " << V(i, 0) << " " << V(i, 1) << " " << 0.0 << endl;
		}
	}
	else if (V.cols() == 3)
	{
		for (int i = 0; i < V.rows(); ++i) {
			outfile << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << endl;
		}
	}

	for (int i = 0; i < F.rows(); ++i) {
		outfile << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1 << endl;

	}
	outfile.close();
}

void write_para(string filepath, vector<int> indexs, vector<Point2f> new_places)
{
	//para 文件的格式为：第一行 锚点个数，第二行锚点的下标，之后每行为每个锚点的目标坐标 
	ofstream outfile(filepath, ios::out);
	if (!outfile.is_open())
	{
		cerr << "open error";
		exit(1);
	}
	int anc_num = indexs.size();
	outfile << anc_num << endl;
	FOR(i,0, anc_num)
		outfile << indexs[i] << " ";
	outfile << endl;
	FOR(i, 0, anc_num)
		outfile << new_places[i].x<<" " << new_places[i].y << " " <<0.0<< endl;
	outfile.close();

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
		cout << "point number not equal: " <<Ts<<"   "<<p2.size()<< endl;
		return 0;
	}
	//cout << "p1.size(): " << p1.size()<<"  "<< p2.size() << endl;
	double result = 0;
	for (int i = 0; i < Ts; i++)
	{
		result += abs(p1[i] - p2[i]);
		//cout << "test: " << result << endl;
	}
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
				//cout << "dis[i][j] : " << dis[i][j] << endl;
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
		cout << shift<<"shift  min_mis: " << distance[first_num - 1][second_num - 1] << endl;
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

double tar_mismatch_fea(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<int> first_fea, vector<int> sec_fea, vector<pair<int, int>>& path, int &sec_shift, int width) //增加特征点的权重
{
	double dis[TAR_num][TAR_num];//两组点之间的坐标差异
	double distance[TAR_num][TAR_num];
	//int step[TAR_num][TAR_num];//记录总的步数
	int first_num = first_arr.size();
	int second_num = second_arr.size(); //first 作为y轴 ,second为x轴
	double fea_ratio = 0.5;
	double half_fea_ratio = 1;
	double com_ratio = 1;
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
	//std::cout << "first: " << first_num << "   second: " << second_num << endl;
	for (int shift = 0; shift < second_num; shift++) //将first固定，分别对齐second的起点
	{
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
		//distance[0][0]记录的是对齐的第一个点
		if (first_fea[0] > 1 && sec_fea[shift] > 1)
		{
			dis[0][0] = fea_ratio*tar_length_2p(first_arr[0], second_arr[shift]);//
		}
		else if (first_fea[0] > 1 || sec_fea[shift] > 1)  dis[0][0] = half_fea_ratio*tar_length_2p(first_arr[0], second_arr[shift]);
		else dis[0][0] = com_ratio*tar_length_2p(first_arr[0], second_arr[shift]);

		distance[0][0] = dis[0][0];
		for (int i = 1; i < first_num; i++)
		{
			if (distance[i][0] == 100000) continue;
			if (first_fea[i] > 1 && sec_fea[shift] > 1)
			{
				dis[i][0] = fea_ratio*tar_length_2p(first_arr[i], second_arr[shift]); //0;
			}
			else if (first_fea[i] > 1 || sec_fea[shift] > 1) dis[i][0] = half_fea_ratio* tar_length_2p(first_arr[i], second_arr[shift]);
			else  dis[i][0] = com_ratio*tar_length_2p(first_arr[i], second_arr[shift]);
			distance[i][0] = distance[i - 1][0] + dis[i][0];
		}
		for (int i = 1; i < second_num; i++)
		{
			if (distance[0][i] == 100000) continue;
			if (first_fea[0] > 1 && sec_fea[(i + shift) % second_num] > 1)
			{
				dis[0][i] = fea_ratio*tar_length_2p(first_arr[0], second_arr[(i + shift) % second_num]); //0;
			}
			else if (first_fea[0] > 1 || sec_fea[(i + shift) % second_num] > 1) dis[0][i] = half_fea_ratio* tar_length_2p(first_arr[0], second_arr[(i + shift) % second_num]);
			else	dis[0][i] = com_ratio*tar_length_2p(first_arr[0], second_arr[(i + shift) % second_num]);
			distance[0][i] = distance[0][i - 1] + dis[0][i];
		}
		for (int i = 1; i < first_num; i++)
		{
			for (int j = 1; j < second_num; j++)
				//(int i = istart; i <= imax; i++)
			{
				if (distance[i][j] == 100000) continue;
				if (first_fea[i] > 1 && sec_fea[(j + shift) % second_num] > 1)
				{
					dis[i][j] = fea_ratio*tar_length_2p(first_arr[i], second_arr[(j + shift) % second_num]); //0;
				}
				else if (first_fea[i] > 1 || sec_fea[(j + shift) % second_num] > 1)  dis[i][j] = half_fea_ratio*tar_length_2p(first_arr[i], second_arr[(j + shift) % second_num]);
				else dis[i][j] = com_ratio*tar_length_2p(first_arr[i], second_arr[(j + shift) % second_num]);
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
			//cout << "min_mis: " << distance[first_num - 1][second_num - 1] << endl;
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
	double thres = 0.001;
	if (abs(dp[i][j] - (dp[i - 1][j - 1] + d[i][j])) < thres) {
		print_TAR_Path(d, dp, i - 1, j - 1, path);

	}
	else if (abs(dp[i][j] - (dp[i][j - 1] + d[i][j])) < thres) {
		print_TAR_Path(d, dp, i, j - 1, path);

	}
	else  if (abs(dp[i][j] - (dp[i - 1][j] + d[i][j])) < thres) {
		print_TAR_Path(d, dp, i - 1, j, path);
	}
	else cout << "Out control!" << endl;
	path.push_back(make_pair(i, j));
}



double isoperimetric_inequality(vector<Point2f> contour)
{
	double arcl = arcLength(contour, true);
	double ratio = 4 * PI * contourArea(contour) / (arcl * arcl);
	return ratio;
}

//edge evaluation and optimization
double bound_collision(vector<Point2f> cont, vector<int> indexes, int type)
{
	int index_size = indexes.size();
	int csize = cont.size();
	vector<vector<Point2f>> four_edges;
	FOR(i, 0, index_size - 1)
	{
		vector<Point2f> edge;
		for (int j = indexes[i]; j <= indexes[i + 1]; j++)  edge.push_back(cont[j]);
		four_edges.push_back(edge);
	}
	vector<Point2f> edge;
	for (int j = indexes[index_size - 1]; j <= indexes[0] + csize; j++)  edge.push_back(cont[j%csize]);
	four_edges.push_back(edge);
	//提取四边
	double score = 0;
	double score1 = 0;
	if (type == -1)
	{
		FOR(i, 0, four_edges.size())
		{
			score += edge_nd_degree(four_edges[i], 0);
			score1 += edge_nd_degree(four_edges[i], 1);
		}
		if (score1 < score) score = score1;
	}
	else if (type == 0)
	{
		FOR(i, 0, four_edges.size()) score += edge_nd_degree(four_edges[i], 0);
	}
	else if (type == 1)
	{
		FOR(i, 0, four_edges.size()) score += edge_nd_degree(four_edges[i], 1);
	}
	//cout << "CCW: " << score << "   CW:" << score1 << endl;
	return score;
}

double edge_nd_degree(vector<Point2f> edge, int type)
{
	int esize = edge.size();
	if (type == 1)
	{
		vector<Point2f> edge_r;
		for (int t = esize - 1; t >= 0; t--) edge_r.push_back(edge[t]);
		edge.swap(edge_r);
	}
	int count = 1;
	Point2f origin_p = edge[0];
	Point2f end_p = edge[esize - 1];
	vector<double> Lengths;
	FOR(i, 0, esize)  Lengths.push_back(length_2p(origin_p, edge[i]));
	double total_angle = 0;
	FOR(j, 1, Lengths.size())
	{
		int col_flag = 0;
		double length_ = Lengths[j];
		for (int m = j + 1; m < Lengths.size(); m++)
		{
			if (Lengths[m] < length_) //计数，计算角度
			{
				col_flag = 1;
				Point2f tan_unit = vertical_vec(origin_p - edge[j]);
				Point2f vec_unit = unit_vec(edge[m] - edge[j]);
				double cos_ = cos_2v(tan_unit, vec_unit);
				//cout << tan_unit << "   " << vec_unit << endl;
				if (cos_ < 0) tan_unit = -tan_unit;
				double sin_ = sin_2v(tan_unit, vec_unit);
				double angle_ = asin(sin_);
				//cout << tan_unit << "   " << vec_unit <<"  "<< angle_<< endl;
				total_angle += abs(angle_);
				//if(sin_)
				break;
			}
		}
		if (col_flag == 1) count++;
	}
	//double colli_possi = count	*1.0 / (esize - 2);
	//cout << count << endl;
	return total_angle / (count*0.5*PI);
}

double edge_nd_opt(vector<Point2f>& edge, int type)
{
	int esize = edge.size();
	if (type == 1)
	{
		vector<Point2f> edge_r;
		for (int t = esize - 1; t >= 0; t--) edge_r.push_back(edge[t]);
		edge.swap(edge_r);
	}
	Point2f origin_p = edge[0];
	Point2f end_p = edge[esize - 1];
	double edge_colli_score = edge_nd_degree(edge, type);
	int count = 0;
	int max_times = 15;
	while (edge_colli_score > 0 && count<max_times)
	{
		count++;
		FOR(i, 1, esize - 1)
		{
			double len1 = length_2p(origin_p, edge[i]);
			FOR(j, i + 1, esize)
			{
				double len2 = length_2p(origin_p, edge[j]);
				if (len2 < len1) //计算角度
				{
					Point2f tan_unit = vertical_vec(origin_p - edge[i]);
					Point2f vec_unit = unit_vec(edge[j] - edge[i]);
					double cos_ = cos_2v(tan_unit, vec_unit);
					//cout << tan_unit << "   " << vec_unit << endl;
					if (cos_ < 0) tan_unit = -tan_unit;
					double sin_ = sin_2v(tan_unit, vec_unit);
					double angle_ = asin(sin_);
					//cout << tan_unit << "   " << vec_unit <<"  "<< angle_<< endl;
					int win_width = 0.5*max_times + 1;
					double max_r = 1.0 / max_times;
					double ratio = max_r / win_width;
					FOR(g, 0, win_width)
					{
						if (i - g > 0 && j + g < esize)
						{
							Point2f rot_center = 0.5*(edge[i - g] + edge[j + g]);
							Mat rot_mat = cv::getRotationMatrix2D(rot_center, (max_r - ratio*g)*angle_ / PI * 180, 1.0);
							vector<Point2f> src = { edge[i - g] ,edge[j + g] };
							vector<Point2f> dst;
							transform(src, dst, rot_mat);
							edge[i - g] = dst[0];
							edge[j + g] = dst[1];
						}
						//cout << "dst: " << dst.size() << "   " << dst[0] << "   " << dst[1] << endl;					
					}
				}
			}
		}
		edge_colli_score = edge_nd_degree(edge, type);
		//cout << "count: " << count << endl;
	}
	return edge_colli_score;
}

double whole_con_opt(vector<Point2f>& cont, vector<int>& indexes, int type)
{
	double degree_opt = 0;
	int index_size = indexes.size();
	int csize = cont.size();
	vector<vector<Point2f>> four_edges;
	FOR(i, 0, index_size - 1)
	{
		vector<Point2f> edge;
		for (int j = indexes[i]; j <= indexes[i + 1]; j++)  edge.push_back(cont[j]);
		four_edges.push_back(edge);
	}
	vector<Point2f> edge;
	for (int j = indexes[index_size - 1]; j <= indexes[0] + csize; j++)  edge.push_back(cont[j%csize]);
	four_edges.push_back(edge);

	FOR(i, 0, four_edges.size())  degree_opt+=edge_nd_opt(four_edges[i], type);
	vector<Point2f> new_c;
	vector<int> new_index;
	for (int m = 0; m < 4; m++)
	{
		new_index.push_back(new_c.size());
		for (int n = 0; n < four_edges[m].size() - 1; n++)
		{
			new_c.push_back(four_edges[m][n]);
		}
	}
	cont = new_c;
	indexes = new_index;
	return degree_opt;
}

void contour_de_crossing(vector<Point2f> &contour_)
{
	int first = 0, second = 0, time = 0;
	int iter_times = 5;
	while (self_intersect(contour_, first, second) && time < iter_times)
	{
		time++;
		vector<Point2f> mid_con;
		int csize = contour_.size();
		vector<Point2f> pre;
		vector<Point2f> post;
		std::cout << "Contour de-crossing: " << csize << " " << first << " " << second << endl;
		if (second - first < 4)
		{
			Point2f tem = contour_[first + 1];
			contour_[first + 1] = contour_[second];
			contour_[second] = tem;
		}
		else
		{
			int mid = (first + second) / 2;
			for (int t = first; t < mid; t++)
			{
				pre.push_back(contour_[t]);
			}
			for (int t = mid; t <= second + 1; t++)
			{
				post.push_back(contour_[t]);
			}
			int tt = 0;
			if ((second + 1) % csize == 0)
				tt = (second + 1) % csize;
			for (; tt < first; tt++)
			{
				mid_con.push_back(contour_[tt]);
			}
			if (mid_con.empty()) mid_con.push_back(0.5*contour_[csize - 1] + 0.5*post[post.size() - 2]);
			else mid_con.push_back(0.5*mid_con.back() + 0.5*post[post.size() - 2]);
			for (int t = post.size() - 2; t >= 0; t--)
			{
				mid_con.push_back(post[t]);
			}
			for (int t = pre.size() - 1; t > 0; t--)
			{
				mid_con.push_back(pre[t]);
			}
			mid_con.push_back(0.5*mid_con.back() + 0.5*contour_[(second + 2) % csize]);
			for (int t = second + 2; t < csize; t++)
			{
				mid_con.push_back(contour_[t]);
			}
			contour_.swap(mid_con);
		}
	}
	if (self_intersect(contour_, first, second))
	{
		cout << "Waring: still self_intersection! " << endl;
	}
	else
		cout << "No self_intersection! " << time<<" iters!"<< endl;
}

void contour_fine_tuning(vector<Point2f> &contour_)  //过近的优化
{
	Mat imagefine = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
	Point2f shift = Point2f(600, 600) - center_p(contour_);
	draw_contour_points(imagefine, contour_, shift, 8, 2);
	int dis_thres = 8;
	int win_width = 4;
	int csize = contour_.size();
	double con_length = arcLength(contour_, true);
	double len_thres = 2 * con_length / csize;
	cout << "length threshold: " << len_thres << endl;
	int exist_close_p = true;
	int iter_times = 0;
	while (exist_close_p)
	{
		exist_close_p = false;
		iter_times++;
		FOR(i, 0, csize)
		{
			FOR(j, i, csize)
			{
				if (abs(j - i) < dis_thres || abs(i + csize - j) < dis_thres)
					continue;
				else
				{
					if (length_2p(contour_[i], contour_[j]) < len_thres) //两个点过近
					{
						exist_close_p = true;
						cout << "before: " << length_2p(contour_[i], contour_[j]) << "   " << contour_[i] << "   " << contour_[j] << endl;
						circle(imagefine, contour_[i] + shift, 3, Scalar(0, 255, 0), -1);
						circle(imagefine, contour_[j] + shift, 3, Scalar(0, 255, 0), -1);
						Point2f unit_v = unit_vec(contour_[i] - contour_[j]);
						cout << unit_v << endl;
						double max_r = 0.25;
						contour_[i] = contour_[i] + max_r*unit_v;
						contour_[j] = contour_[j] - max_r*unit_v;
						FOR(g, 1, win_width)
						{
							double ratio = max_r - max_r/ win_width*g;
							contour_[(i + g) % csize] += ratio*unit_v;
							contour_[(i - g + csize)% csize] += ratio*unit_v;
							contour_[(j + g) % csize] -= ratio*unit_v;
							contour_[(j - g + csize) % csize] -= ratio*unit_v;
							//cout << "dst: " << dst.size() << "   " << dst[0] << "   " << dst[1] << endl;					
						}
						circle(imagefine, contour_[i] + shift, 3, Scalar(255, 0, 0), -1);
						circle(imagefine, contour_[j] + shift, 3, Scalar(255, 0, 0), -1);
						cout << "after: " << length_2p(contour_[i], contour_[j]) << contour_[i] << "   " << contour_[j] << endl;
					}
				}
			}
		}
	}
	draw_contour_points(imagefine, contour_, shift+Point2f(300,300), 9, 2);
	cout << "fine tuning times: " << iter_times << endl;
	imshow("contour_dst", imagefine);
}

int triangulateContour(vector<Point2f>& contour, MatrixXd& V, MatrixXi& F)
{
	vector<Point2f> con_ori = contour;
	Point2f cen = center_p(con_ori);
	int ori_num = add_points(contour, 0.08);
	// create a Subdiv2D object from the contour
	cv::Subdiv2D subdiv(Rect(cen.x - 500, cen.y - 500, 1000, 1000)); // change the rectangle size according to your contour
	for (int i = 0; i < contour.size(); i++) {
		subdiv.insert(contour[i]);
	}
	//calculate V
	V.resize(contour.size(), 2);
	for (int m = 0; m < contour.size(); m++)
	{
		V.row(m) << contour[m].x, contour[m].y;
		//cout << "rest: " << V.row(m) << endl;
	}

	// get the triangle list
	std::vector<Vec6f> triangleList;
	subdiv.getTriangleList(triangleList);
	// convert the triangle list to vertex and face matrices
	//V.resize(3 * triangleList.size(), 2);
	F.resize(triangleList.size(), 3);
	int i = 0;
	for (auto& t : triangleList) {
		for (int j = 0; j < 3; j++)
		{
			Point2f p(t[2 * j], t[2 * j + 1]);
			//cout << "tri: " << p.x << "   " << p.y << endl;
			//uset.emplace(p);
			auto it = find(contour.begin(), contour.end(), p);//vertexs.find(p);
			if (it != contour.end()) {
				//V.row(3 * i + j) << t[2 * j], t[2 * j + 1];
				F(i, j) = distance(contour.begin(), it);
			}
			else
			{
				//V.row(3 * i + j) << t[2 * j], t[2 * j + 1];
				F(i, j) = contour.size() - 1;
			}
		}
		i++;
	}
	cout << "contour.size " << contour.size() << "   triangleList.size " << triangleList.size() << "  V and F: " << V.size() << "  " << F.size() << endl;
	//for (int n = 0; n < F.size() / 3; n++) cout << "n row: "<<F.row(n) << endl;
	//cout << "triangleList.size(): " << triangleList.size() << "   " << F.size() << endl;
	// remove triangles outside the contour
	std::vector<bool> keep(F.rows(), false);
	for (int i = 0; i < F.rows(); i++) {
		Vector2d p1 = V.row(F(i, 0));
		Vector2d p2 = V.row(F(i, 1));
		Vector2d p3 = V.row(F(i, 2));
		if (pointPolygonTest(con_ori, 1.0 / 3 * (Point2f(p1[0], p1[1]) + Point2f(p2[0], p2[1]) + Point2f(p3[0], p3[1])), false) >= 0) {
			keep[i] = true;
		}
		//if (keep[i]) cout << "true!" << endl;
	}
	int numFaces = keep.size() - count(keep.begin(), keep.end(), false);
	cout << "numFaces: " << numFaces << endl;
	MatrixXi newF(numFaces, 3);
	int j = 0;
	for (int i = 0; i < F.rows(); i++) {
		if (keep[i]) {
			newF.row(j) = F.row(i);
			j++;
		}
	}
	F = newF;
	//V.conservativeResize(F.maxCoeff() + 1, 2);
	/*cout << "  V : " << V.size() << "  F: " << F.size() << endl;
	for (int n = 0; n < V.rows(); n++) cout << n << " V row: " << V.row(n) << endl;
	for (int n = 0; n < F.rows(); n++) cout << n << " F row: " << F.row(n) << endl;*/
	return ori_num;
}


vector<Point2f> triangulate_2Contours(vector<Point2f>& cont1, vector<Point2f>& cont2, MatrixXd& V, MatrixXi& F)
{
	vector<Point2f> con_tri;
	vector<Point2f> con_union;
	vector<Point2f> con_intersection;
	int csize = cont1.size();
	if (csize != cont2.size())
		cout << "different size of these contours in the triangulate_2Contours!" << endl;
	vector<int> endp1 = { 1000, -1 };
	vector<int> endp2 = { 1000, -1 };
	FOR(i, 0, csize)
	{
		FOR(j, 0, csize)
		{
			if (cont1[i] == cont2[j])
			{
				if (i < endp1[0])  endp1[0] = i;
				if (i > endp1[1])  endp1[1] = i;
				if (j < endp2[0])  endp2[0] = j;
				if (j > endp2[1])  endp2[1] = j;
			}
		}
	}
	cout << endp1[0] << "  " << endp1[1] << "  " << endp2[0] << "  "<<endp2[1] << endl;
	FOR(i, endp1[1], endp1[0] + csize) con_union.push_back(cont1[i%csize]);
	int size1 = con_union.size();
	//cout << "con_union size: " << con_union[0]<< "   " << con_union.back() << endl;
	FOR(i, endp2[1], endp2[0] + csize) con_union.push_back(cont2[i%csize]);
	//cout << "con_union2 size: " << con_union[size1] << "   " << con_union.back() << endl;
	//FOR(i, endp1[0] + 1, endp1[1]) con_intersection.push_back(cont1[i]);
	//cout << "con_intersection size: " << con_intersection[0] << "   " << con_intersection.back() << endl;
	//int ori_num1 = add_points(cont1, 0.1);
	//int ori_num2 = add_points(cont2, 0.1);
	con_tri = con_union;
	int ori_num = add_points(con_tri, 0.04);

	//con_tri.insert(con_tri.end(), con_intersection.begin(), con_intersection.end());
	//FOR(i, ori_num1, cont1.size())  con_tri.push_back(cont1[i]);
	//FOR(i, ori_num2, cont2.size())  con_tri.push_back(cont2[i]);
	Point2f cen = center_p(con_tri);
	// create a Subdiv2D object from the contour
	cv::Subdiv2D subdiv(Rect(cen.x - 600, cen.y - 600, 1200, 1200)); // change the rectangle size according to your contour
	for (int i = 0; i < con_tri.size(); i++) {
		subdiv.insert(con_tri[i]);
	}
	//calculate V
	V.resize(con_tri.size(), 2);
	for (int m = 0; m < con_tri.size(); m++)
	{
		V.row(m) << con_tri[m].x, con_tri[m].y;
		//cout << "rest: " << V.row(m) << endl;
	}

	// get the triangle list
	std::vector<Vec6f> triangleList;
	subdiv.getTriangleList(triangleList);
	// convert the triangle list to vertex and face matrices
	//V.resize(3 * triangleList.size(), 2);
	F.resize(triangleList.size(), 3);
	int i = 0;
	for (auto& t : triangleList) {
		for (int j = 0; j < 3; j++)
		{
			Point2f p(t[2 * j], t[2 * j + 1]);
			//cout << "tri: " << p.x << "   " << p.y << endl;
			//uset.emplace(p);
			auto it = find(con_tri.begin(), con_tri.end(), p);//vertexs.find(p);
			if (it != con_tri.end()) {
				//V.row(3 * i + j) << t[2 * j], t[2 * j + 1];
				F(i, j) = distance(con_tri.begin(), it);
			}
			else
			{
				//V.row(3 * i + j) << t[2 * j], t[2 * j + 1];
				F(i, j) = con_tri.size() - 1;
			}
		}
		i++;
	}
	cout << "contour.size " << con_tri.size() << "   triangleList.size " << triangleList.size() << "  V and F: " << V.size() << "  " << F.size() << endl;
	//for (int n = 0; n < F.size() / 3; n++) cout << "n row: "<<F.row(n) << endl;
	//cout << "triangleList.size(): " << triangleList.size() << "   " << F.size() << endl;
	// remove triangles outside the contour
	std::vector<bool> keep(F.rows(), false);
	for (int i = 0; i < F.rows(); i++) {
		Vector2d p1 = V.row(F(i, 0));
		Vector2d p2 = V.row(F(i, 1));
		Vector2d p3 = V.row(F(i, 2));
		if (pointPolygonTest(con_union, 1.0 / 3 * (Point2f(p1[0], p1[1]) + Point2f(p2[0], p2[1]) + Point2f(p3[0], p3[1])), false) >= 0) {
			keep[i] = true;
		}
		//if (keep[i]) cout << "true!" << endl;
	}
	int numFaces = keep.size() - count(keep.begin(), keep.end(), false);
	cout << "numFaces: " << numFaces << endl;
	MatrixXi newF(numFaces, 3);
	int j = 0;
	for (int i = 0; i < F.rows(); i++) {
		if (keep[i]) {
			newF.row(j) = F.row(i);
			j++;
		}
	}
	F = newF;
	//V.conservativeResize(F.maxCoeff() + 1, 2);
	/*cout << "  V : " << V.size() << "  F: " << F.size() << endl;
	for (int n = 0; n < V.rows(); n++) cout << n << " V row: " << V.row(n) << endl;
	for (int n = 0; n < F.rows(); n++) cout << n << " F row: " << F.row(n) << endl;*/
	Mat image_o = Mat(1600, 1600, CV_8UC3, Scalar(255, 255, 255));
	draw_contour_points(image_o, con_union, Point2f(600,600)-center_p(con_union), 5, 2);
	draw_contour_points(image_o, con_intersection, Point2f(600, 600) - center_p(con_union), 6, 2);
	draw_contour_points(image_o, con_tri, Point2f(1200, 600) - center_p(con_union), 5, 2);
	imshow("two triangle: ", image_o);
	//cout << con_union.size() << "  " << con_intersection.size() << "  " << ori_num1 << "   " << ori_num2 << "  " << con_tri.size() << endl;
	return con_tri;
}


vector<Point2f> triangulate_Contours_bbx(vector<Point2f>& cont1, vector<int> anc1, MatrixXd& V, MatrixXi& F)
{
	vector<Point2f> con_tri;
	vector<Point2f> con_union;
	vector<Point2f> con_intersection;
	vector<int> anc_tri;
	int csize = cont1.size();
	Point2f line1 = cont1[anc1[2]] - cont1[anc1[0]];
	Point2f line2 = cont1[anc1[3]] - cont1[anc1[1]];
	vector<Point2f> shifting;
	shifting.push_back(line1);
	shifting.push_back(line1 + line2);
	shifting.push_back(line2);

	// 提取平移围成的轮廓
	vector<vector<Point2f>> four_place;
	four_place.push_back(cont1);
	for (int i = 0; i < 3; i++)
	{
		vector<Point2f> one_loca;
		for (int j = 0; j < csize; j++)
		{
			one_loca.push_back(cont1[j]+ shifting[i]);
		}
		four_place.push_back(one_loca);
	}

	int total_num = 0;
	for (int i = 0; i < four_place.size(); i++)
	{
		anc_tri.push_back(total_num);
		int s_i = (i + 3) % 4;
		int e_i = (i + 2) % 4;
		int t_end = anc1[e_i];
		if (s_i>e_i)  t_end += csize;
		for (int t = anc1[s_i]; t < t_end; t++)                   //translation的提取规律为1的4-3, 2的1-4, 3的2-1, 4的3-2
		{
			con_tri.push_back(four_place[i][t % csize]);
			total_num++;
			//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
		}
	}
	con_union = con_tri;
	FOR(g, 0, 4)
	{
		cout << anc_tri[g] << "   "<< con_tri[anc_tri[g]]<<endl;
	}
	
	cout << "total_num: " << total_num << endl;

	for (int i = 0; i < four_place.size(); i++)
	{
		int s_i = (i + 3) % 4;
		int e_i = (i + 2) % 4;
		int ts = anc1[s_i];
		if (s_i<e_i)  ts += csize;
		for (int t = ts - 1; t > anc1[e_i]; t--)                   //translation的提取规律为1的4-3, 2的1-4, 3的2-1, 4的3-2
		{
			con_intersection.push_back(four_place[i][t % csize]);
		}
	}
	con_tri.insert(con_tri.end(), con_intersection.begin(), con_intersection.end());
	cout << "total_num2: " << con_tri.size() << endl;
	for (int i = 0; i < four_place.size(); i++)
	{
		int ori_num = add_points(four_place[i], 0.08);
		FOR(j, ori_num, four_place[i].size())  con_tri.push_back(four_place[i][j]);
	}
	int ori_num2 = add_points(con_intersection, 0.08);
	FOR(j, ori_num2, con_intersection.size())  con_tri.push_back(con_intersection[j]);
	cout << "total_num3: " << con_tri.size() << endl;

	Point2f cen = center_p(con_tri);
	// create a Subdiv2D object from the contour
	cv::Subdiv2D subdiv(Rect(cen.x - 800, cen.y - 800, 1600, 1600)); // change the rectangle size according to your contour
	for (int i = 0; i < con_tri.size(); i++) {
		subdiv.insert(con_tri[i]);
	}
	//calculate V
	V.resize(con_tri.size(), 2);
	for (int m = 0; m < con_tri.size(); m++)
	{
		V.row(m) << con_tri[m].x, con_tri[m].y;
		//cout << "rest: " << V.row(m) << endl;
	}

	// get the triangle list
	std::vector<Vec6f> triangleList;
	subdiv.getTriangleList(triangleList);
	// convert the triangle list to vertex and face matrices
	//V.resize(3 * triangleList.size(), 2);
	F.resize(triangleList.size(), 3);
	int i = 0;
	for (auto& t : triangleList) {
		for (int j = 0; j < 3; j++)
		{
			Point2f p(t[2 * j], t[2 * j + 1]);
			//cout << "tri: " << p.x << "   " << p.y << endl;
			//uset.emplace(p);
			auto it = find(con_tri.begin(), con_tri.end(), p);//vertexs.find(p);
			if (it != con_tri.end()) {
				//V.row(3 * i + j) << t[2 * j], t[2 * j + 1];
				F(i, j) = distance(con_tri.begin(), it);
			}
			else
			{
				//V.row(3 * i + j) << t[2 * j], t[2 * j + 1];
				F(i, j) = con_tri.size() - 1;
			}
		}
		i++;
	}
	cout << "contour.size " << con_tri.size() << "   triangleList.size " << triangleList.size() << "  V and F: " << V.size() << "  " << F.size() << endl;
	//for (int n = 0; n < F.size() / 3; n++) cout << "n row: "<<F.row(n) << endl;
	//cout << "triangleList.size(): " << triangleList.size() << "   " << F.size() << endl;
	// remove triangles outside the contour
	std::vector<bool> keep(F.rows(), false);
	for (int i = 0; i < F.rows(); i++) {
		Vector2d p1 = V.row(F(i, 0));
		Vector2d p2 = V.row(F(i, 1));
		Vector2d p3 = V.row(F(i, 2));
		if (pointPolygonTest(con_union, 1.0 / 3 * (Point2f(p1[0], p1[1]) + Point2f(p2[0], p2[1]) + Point2f(p3[0], p3[1])), false) >= 0) {
			keep[i] = true;
		}
		//if (keep[i]) cout << "true!" << endl;
	}
	int numFaces = keep.size() - count(keep.begin(), keep.end(), false);
	cout << "numFaces: " << numFaces << endl;
	MatrixXi newF(numFaces, 3);
	int j = 0;
	for (int i = 0; i < F.rows(); i++) {
		if (keep[i]) {
			newF.row(j) = F.row(i);
			j++;
		}
	}
	F = newF;

	Mat image_o = Mat(1600, 1600, CV_8UC3, Scalar(255, 255, 255));
	draw_contour_points(image_o, con_union, Point2f(600, 600) - center_p(con_union), 5, 2);
	draw_contour_points(image_o, con_intersection, Point2f(600, 600) - center_p(con_union), 6, 2);
	draw_contour_points(image_o, con_tri, Point2f(1200, 600) - center_p(con_union), 5, 2);
	imshow("two triangle: ", image_o);
	return con_tri;
}


int add_points(vector<Point2f>& contour, double sparse_ratio)
{
	int add_index = contour.size();
	vector<Point2f> bbx_p = bbx(contour);
	double minx = bbx_p[1].x;
	double miny = bbx_p[1].y;
	double maxx = bbx_p[3].x;
	double maxy = bbx_p[3].y;
	double interval_x = sparse_ratio*(maxx - minx);
	double interval_y = sparse_ratio*(maxy - miny);
	//cout << "bbx:"<<maxx << "   " << maxy<<"  "<< minx<<"   "<< miny << endl;
	double margin = 0.5*min(interval_x, interval_y);
	vector<Point2f> new_con;
	for (double nx = minx + interval_x; nx < maxx; nx += interval_x)
	{
		for (double ny = miny + interval_y; ny < maxy; ny += interval_y)
		{
			Point2f new_p(nx, ny);
			double min_dis = pointPolygonTest(contour, new_p, true);
			if (min_dis > margin)
			{
				new_con.push_back(new_p);
			}
		}
	}
	contour.insert(contour.end(), new_con.begin(), new_con.end());
	bool show = true;
	if (show)
	{
		Mat image_o = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
		draw_contour_points(image_o, contour, OP, 5, 2);
		imshow("add points: ", image_o);
		cout << "Add " << contour.size() - add_index << " new points!" << endl;
	}
	return add_index;
}

int point_locate(vector<Point2f> con, Point2f p)
{
	double thres_len = 0.01;
	//double min_length = 1000;
	for (int i = 0; i < con.size(); i++) 
	{
		double length = length_2p(con[i], p);
		if (length < thres_len)
			return i;
	}
	cout << "Not in the contour!" << endl;
	return -1;
}





