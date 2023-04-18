#include "tilingOpt.h"
vector<string> frame_type={ "parallelogram", "rhombus", "rectangle", "square" };
bool pers_trans = 1;
bool coll_opt = 0;
bool deve_opt = 0;

namespace Tiling_tiles {

	void Tiling_opt::tiliing_generation(string nameid)
	{
		load_dataset(true);
		string filepath = DefaultPath;
		filepath = filepath + "contour/" + nameid + ".txt";
		prototile_first = protoTile(filepath);
		string savepath = SavePath;
		savepath += nameid;
		const char *na = savepath.c_str();
		if (_access(na, 0) != -1) printf("The  file/dir had been Exisit \n");
		else	_mkdir(na);
		vector<int> cand_points = prototile_first.getCandPoints(prototile_first.contour_f);
		int trans = Tanslation_rule(cand_points, prototile_first.contour_f, savepath);
		//int flips = Flipping_rule(p_p_index, cont_orig, rootname);
		int all_inner_num = all_inner_conts.size();
		if (all_inner_num == 0)
		{
			std::cout << "no right placement" << endl;
			return;
		}
		FOR(i, 0, 0)
		//FOR(i, 0, all_inner_num)
		{
			match_candidate(i);
			for (int j = 0; j < 1; j++)
			{
				Mat draw = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
				vector<pair<int, int>> path = cand_paths[j];
				vector<Point_f> contour1 = prototile_mid.contour_f;
				vector<Point_f> contour2 = candidate_contours[j];
				Point2f sh = Point2f(300, 500) - center_p(conf_trans(contour1));
				draw_pair(draw, conf_trans(contour1), conf_trans(contour2), path, sh);
				
				double con_sc;
				vector<vector<double>> contour2_tar;
				//求contour_对应的point_f和tar值
				//contour2_tar = computeTAR(contour_, con_sc, 0.5);
				double angle = 165;
				double ratio_max = 100;
				double result_score = 0;
				double ratio = 0.5;
				vector<Point_f> contour_2 = morphing_dir(contour1, contour2, path, ratio);
				Point2f sh2 = Point2f(700, 500) - center_p(conf_trans(contour_2));
				draw_contour_points(draw, conf_trans(contour_2), sh2, 3, 2);
				imshow("correspond path", draw);
				vector<Point_f> con_re;
				vector<int> anc_mid;
				vector<int> anc_re;
				FOR(t, 0, contour_2.size())
				{
					if (contour_2[t].type == fixed_p)
						anc_mid.push_back(t);
				}
				//vector<Point2f> frame = { contour_2[anc_mid[0]].point, contour_2[anc_mid[1]].point, contour_2[anc_mid[2]].point, contour_2[anc_mid[3]].point };
				//vector<Point2f> frame_b = base_frame(frame, 3);
				//Mat rot_mat = getPerspectiveTransform(frame, frame_b);// getAffineTransform(frame_1, frame_b_1);
				//													  //cout << rot_mat << endl;
				//vector<Point2f> contour_dst;
				//perspectiveTransform(conf_trans(contour_2), contour_dst, rot_mat);
				//
				//whole_con_opt(contour_dst, anc_mid, 0);
				//contour_2 = set_flags(contour_dst, contour_2);

				Mat draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
				if (translation_placement(contour_2, con_re, anc_mid, anc_re, draw2))
				{
					cout << "OK!" << endl;
				}

				protoTile c1, c2;
				c1.show_contour(conf_trans(contour_2), anc_mid);
				c2.show_contour(conf_trans(con_re), anc_re);
				RotationVis(c1,c2, AntiClockWise);
				//for (; ratio <= 90; ratio += 5)
				//{
				//	//ratio越高，中间形状保持越多，对应的morphed_A的形状保持越多，eve2越高
				//	double ratio_d = ratio / 100.0;
				//	cout << "ratio_d: " << ratio_d << endl;
				//	//这里final_pettern的点数不一定为200点，但是后续的计算中会重新采样
				//	vector<Point2f> final_pettern = morphing(contour1, contour2, path, ratio_d);
				//	vector<int> return_p;
				//	vector<vector<Point2f>> four_place;
				//	vector<Point2f> morphed_A = extract_contour(final_pettern, mid_inter_morphed, return_p, four_place, all_inner_conts[i].type);
				//	return_p.swap(vector<int>());
				//	four_place.swap(vector<vector<Point2f>>());
				//	vector<Point2f> morphed_A_ori = extract_contour(prototile_mid->contour, mid_inter, return_p, four_place, all_inner_conts[i].type);

				//	double evaluation1 = evalua_deformation(final_pettern, contour_);
				//	double evaluation2 = evalua_deformation(morphed_A, morphed_A_ori);
				//	double score_ave = (1 - ratio_d)*evaluation2 + ratio_d*evaluation1;
				//	if (score_ave > result_score)
				//	{
				//		result_score = score_ave;
				//		ratio_max = ratio_d;
				//	}
				//}

				//vector<Point_f> final_pettern = morphing(contour1, contour2, path, ratio_max);

				//vector<int> return_p;
				//vector<vector<Point2f>> four_place;
				//vector<Point2f> morphed_A = extract_contour(final_pettern, mid_inter_morphed, return_p, four_place, all_inner_conts[i].type);
				////vector<int> return_p2;
				////vector<vector<Point2f>> four_place2;
				////vector<Point2f> morphed_A_ori = extract_contour(prototile_mid->contour, mid_inter, return_p2, four_place2, all_inner_conts[i].type);

				//cout << "ratio_max: " << ratio_max << "  --  result_score:  " << result_score << endl;
				//Mat drawing_pro = Mat(800, 3200, CV_8UC3, Scalar(255, 255, 255));
				//draw_poly(drawing_pro, prototile_mid->contour, Point2f(400, 400));
				//draw_poly(drawing_pro, contour_, Point2f(1200, 400));
				//draw_poly(drawing_pro, final_pettern, Point2f(2000, 400));
				//prototile_tem->loadPoints(final_pettern);
				//draw_poly(drawing_pro, morphed_A, Point2f(2800, 400));
				//string filename = rootname + "\\" + int2string(i) + "_Candidate_" + int2string(j) + "_" + double2string(ratio_max) + "_" + double2string(result_score) + ".png";
				//imwrite(filename, drawing_pro);

				////halftone generate
				//if (!contour_is_simple(morphed_A) || !contour_is_simple(final_pettern))
				//{
				//	cout << "--------------Warning!!!-------------"
				//		<< endl << "The " << j << " candidate patterns are not simple! No results!" << endl << endl;
				//	continue;
				//}
				//int halftone_num = 5000;
				//double l_offset = 2;
				//double scale_t = 0.5;
				//Mat drawing_tesse = Mat(halftone_num, halftone_num, CV_8UC3, Scalar(255, 255, 255));
				//vector<vector<Point2f>> all_tiles = tesse_all(morphed_A, return_p, all_inner_conts[i].type, halftone_num, l_offset, scale_t);
				//int alltilesnum = all_tiles.size();
				//cout << "The num of tiles in final tesse result: " << alltilesnum << endl;
				//for (int m = 0; m < alltilesnum; m++)
				//{
				//	//框架图
				//	//vector<Point2f> tile_dilate = contour_dilate(all_tiles[m], 5);
				//	//vector<Point2f> tile_erode = contour_erode(all_tiles[m], 8);
				//	////draw_poly(drawing_tesse, tile_dilate, center_p(tile_dilate));
				//	////draw_poly(drawing_tesse, tile_erode, center_p(tile_erode), 1);
				//	draw_poly(drawing_tesse, all_tiles[m], center_p(all_tiles[m]));
				//}
				//string filename2 = rootname + "\\" + int2string(i) + "_TessellationResult_" + int2string(j) + ".png";
				//imwrite(filename2, drawing_tesse);
			}
		}


	}

	void Tiling_opt::tiliing_gen_specify(string nameid, vector<int> anc_points)
	{
		load_dataset(false);
		string filepath = DefaultPath;
		filepath = filepath + "contour/" + nameid + ".txt";
		prototile_first = protoTile(filepath);
		string savepath = SavePath;
		savepath = savepath.substr(0, savepath.length() - 7)+"result_spe/";
		savepath += nameid;
		const char *na = savepath.c_str();
		if (_access(na, 0) != -1) printf("The  file/dir had been Exisit \n");
		else	_mkdir(na);
		vector<int> cand_points = prototile_first.getCandPoints(prototile_first.contour_f);
		int trans = Tanslation_rule_spec(cand_points, prototile_first.contour_f, savepath, anc_points);
		//int flips = Flipping_rule(p_p_index, cont_orig, rootname);
		int all_inner_num = all_inner_conts.size();
		if (all_inner_num == 0)
		{
			std::cout << "no right placement" << endl;
			return;
		}
		FOR(i, 0, 1)
			//FOR(i, 0, all_inner_num)
		{
			match_candidate_(i);
			//feature_match();
			for (int j = 0; j < 1; j++)
			{
				vector<pair<int, int>> path = cand_paths[j];
				cout << "path num: " << path.size() << endl;
				vector<Point_f> contour1 = prototile_mid.contour_f;
				vector<Point_f> contour2 = candidate_contours[j];
				prototile_second = protoTile(contour2);
				vector<pair<int, int>> path_fea = cand_fea_paths[j];

				//show the path
				Mat draw = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
				Point2f sh = Point2f(300, 500) - center_p(conf_trans(contour1));
				draw_pair(draw, conf_trans(contour1), conf_trans(contour2), path, sh);
				imshow("pair  match", draw);
				Point2f shh2 = sh + Point2f(400, 0);
				Mat draw22 = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
				draw_contour_points(draw22, prototile_mid.contour, sh, 5, 2);
				draw_contour_points(draw22, prototile_second.contour, shh2, 7, 2);
				FOR(i, 0, prototile_mid.feature_points.size())
					circle(draw22, prototile_mid.contour[prototile_mid.feature_points[i]]+ sh, 3, Scalar(0, 0, 255), -1);
				FOR(i, 0, prototile_second.feature_points.size())
					circle(draw22, prototile_second.contour[prototile_second.feature_points[i]] + shh2, 3, Scalar(0, 125, 255), -1);
				circle(draw22, prototile_mid.contour[prototile_mid.feature_points[0]]+ sh, 5, Scalar(120, 0, 255), -1);
				circle(draw22, prototile_second.contour[prototile_second.feature_points[0]] + shh2, 5, Scalar(0, 125, 255), -1);
				circle(draw22, prototile_mid.contour[0] + sh, 6, Scalar(25, 200, 25), 2);
				circle(draw22, prototile_second.contour[0] + shh2, 6, Scalar(25, 200, 25), 2);
				FOR(i, 0, path_fea.size())
				{
					//int tsfsize = prototile_second.feature_points.size();
					Point2f f1 = prototile_mid.contour[path_fea[i].first] + sh;
					Point2f f2 = prototile_second.contour[path_fea[i].second] + shh2;
					line(draw22, f1, f2, colorbar[6].second);
				}
				imshow("fea match", draw22);
				//----------show the feature------------

				double con_sc;
				vector<vector<double>> contour2_tar;
				//求contour_对应的point_f和tar值
				//contour2_tar = computeTAR(contour_, con_sc, 0.5);
				double angle = 165;
				double ratio_max = 100;
				double result_score = 0;
				double ratio = 0.5;
				vector<Point_f> contour_2 = morphing(contour1, contour2, path_fea, ratio);
				cout << "contour_2.size:  " << contour_2.size() << "   " << contour1.size() << endl;
				//FOR(mm, 0, contour_2.size()) cout << contour_2[mm].type << "   " << contour1[mm].type << endl;
				//vector<Point_f> contour_2 = morphing_dir(contour1, contour2, path, ratio);
				Point2f sh2 = Point2f(700, 500) - center_p(conf_trans(contour_2));
				draw_contour_points(draw, conf_trans(contour_2), sh2, 3, 2);
				imshow("correspond path", draw);

				//calculate two deformed shapes
				vector<Point_f> con_re;
				vector<int> anc_mid;
				vector<int> anc_re;
				int csize = contour_2.size();
				FOR(t, 0, csize)
				{
					if (contour_2[t].type == fixed_p)
						anc_mid.push_back(t);
				}
				Mat draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
				if (translation_spec(contour_2, con_re, anc_mid, anc_re, draw2))
				{
					cout << "OK! No intersection!" << endl;
				}

				vector<Point2f> contour_dst = conf_trans(contour_2);
				vector<Point2f> cont_re = conf_trans(con_re);
				Mat draw3 = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
				draw_contour_points(draw3, contour_dst, Point2f(600, 600) - center_p(contour_dst), 5, 2);
				draw_contour_points(draw3, cont_re, Point2f(600, 600) - center_p(contour_dst), 6, 2);
				imshow("2 contours:", draw3);
				MatrixXd V;
				MatrixXi F;
				//vector<Point2f> tt = triangulate_2Contours(contour_dst, cont_re,V,F);
				//vector<Point2f> tt = triangulate_Contours_bbx(cont_re, anc_re, V, F);
				//vector<Point2f> tt = triangulate_Contours_bbx(contour_dst, anc_mid, V, F);
				vector<Point2f> tt = triangulate_bbx(contour_dst, V, F);
				Mat image = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
				Point stf(0, 600);
				for (size_t i = 0; i < F.rows(); i++)
				{
					Point pt1(V.row(F(i, 0)).x(), V.row(F(i, 0)).y());
					Point pt2(V.row(F(i, 1)).x(), V.row(F(i, 1)).y());
					Point pt3(V.row(F(i, 2)).x(), V.row(F(i, 2)).y());
					//cout << F(i, 0) << "   " << F(i, 1) << "   " << F(i, 2) << endl;
					//cout << pt1 << "    " << pt2 << "   " << pt3 << endl;
					line(image, pt1+ stf, pt2+ stf, Scalar(0, 255, 0), 1, LINE_AA);
					line(image, pt2+ stf, pt3+ stf, Scalar(0, 255, 0), 1, LINE_AA);
					line(image, pt3+ stf, pt1+ stf, Scalar(0, 255, 0), 1, LINE_AA);
				}
				//imwrite("Triangulation.png", image);
				imshow("Triangulation", image);
				double degree_after_opt = 0;
				if (pers_trans)
				{
					vector<Point2f> frame = { contour_dst[anc_mid[0]], contour_dst[anc_mid[1]], contour_dst[anc_mid[2]], contour_dst[anc_mid[3]] };
					//vector<Point2f> frame = { cont_re[anc_re[0]], cont_re[anc_re[1]], cont_re[anc_re[2]], cont_re[anc_re[3]] };
					int min_type = 0;
					double min_l = 10000;
					vector<Point2f> frame_b;// = base_frame(frame, 0);
					FOR(f_type, 0, 4)
					{
						vector<Point2f> frame_bb = base_frame(frame, f_type);
						vector<pair<int, int>> path_min = { make_pair(0,0),make_pair(1,1),make_pair(2,2),make_pair(3,3) };
						double alige_e = conotour_align(frame, frame_bb, path_min);
						if (alige_e < min_l)
						{
							min_l = alige_e;
							min_type = f_type;
							frame_b = frame_bb;
						}
					}
					cout << "frame_type: " << frame_type[min_type] << "   " << min_l << endl;
					vector<int> anc_mid_;
					FOR(g, 0, 4)
					{
						cout << "frame[g]: " << frame[g] << endl;
						anc_mid_.push_back(point_locate(tt, frame[g]));
					}
						
					// 构建handle area
					//vector<int> handle_area;
					//vector<Point2f> handle_points;
					//FOR(index_num, 0, anc_mid.size())
					//{
					//	Point2f shift_new = frame_b[index_num] - frame[index_num];
					//	int area_width = handle_area_width;
					//	vector<Point2f> handle_one;
					//	FOR(width_n, -area_width, area_width + 1) 
					//	//for (int width_n = -area_width; width_n < area_width + 1; width_n += 2)
					//	{
					//		int index_ = (anc_mid[index_num] + width_n + csize) % csize;
					//		handle_area.push_back(index_);
					//		//handle_points.push_back(contour_2[index_].point + shift_new);
					//		handle_one.push_back(contour_2[index_].point + shift_new);
					//	}
					//	//cout << "handle_one before: " << handle_one[0] << "  " << handle_one[1] << "  " << handle_one[1] << endl;
					//	int cen_index = handle_one.size() / 2;
					//	Point2f cen_line = handle_one[0] + handle_one[handle_one.size() - 1] - 2*handle_one[cen_index];
					//	double cos1 = cos_2v(cen_line, shift_new);
					//	double sin1 = sin_2v(cen_line, shift_new);
					//	double angle = acos(cos1) / PI * 180;
					//	if (sin1 < 0) angle = -angle;		
					//	if (abs(angle)>90) angle = 0.5*(angle - angle / abs(angle) * 180);
					//	else angle = 0.5*angle;
					//	//cout << cen_line << " " << shift_new << cos1 << "  " << angle << endl;
					//	handle_one = Rotate_contour(handle_one, handle_one[cen_index], angle);
					//	//cout << "handle_one after: " << handle_one[0] << "  " << handle_one[1] << "  " << handle_one[1] << endl;
					//	for(auto p: handle_one)
					//		handle_points.push_back(p);
					//}
					Mat draw_ = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
					Point2f sh1=Point2f(600, 600)-center_p(cont_re);
					draw_contour(draw_, cont_re, sh1, 5);
					FOR(m, 0, 4) circle(draw_, frame[m] + sh1, 3, Scalar(125, 0, 0));
					//
					//MatrixXd V;
					//MatrixXi F;
					//vector<Point2f> con = contour_dst;
					////fileout("D:/vs2015project/Dihedral_new/Dihedral_new/mid_shape.txt", con);
					//int consize = con.size();
					//triangulateContour(con, V, F);

					//vector<Point2f> ttt = { tt[0],tt[1],tt[2],tt[3] };
					//Mat rot_mat = getPerspectiveTransform(frame, frame_b);// getAffineTransform(frame_1, frame_b_1);
					//perspectiveTransform(ttt, frame_b, rot_mat);
					/*vector<int> amid_ = { anc_mid_[0] + 4,anc_mid_[1] + 4,anc_mid_[2] + 4,anc_mid_[3] + 4 };
					anc_mid_ = amid_;*/

					string parap = ParaPath;
					string obj_path = parap +"mid_result.obj";
					string para_path = parap + "deform_para.txt";
					string deformed_c = parap + "deformed_c.txt";
					write_obj(obj_path, V, F);
					write_para(para_path, anc_mid_, frame_b);
					//write_para(para_path, anc_mid, frame_b);
					//write_para(para_path, handle_area, handle_points);
					string command = "D:/vs2015project/Dihedral_new/Dihedral_new/ARAP_Deform.exe  " + obj_path + "  " + para_path + " " + deformed_c; 
					//string command = "D:/vs2015project/ARAP_Deform/x64/Debug/ARAP_Deform.exe  " + obj_path + "  " + para_path + " " + deformed_c;
					cout << command << endl;
					system(command.c_str());
					contour_dst = load_point_file(deformed_c);
					vector<Point2f> con_tem;
					FOR(m, 0, csize) con_tem.push_back(contour_dst[m]);
					contour_dst = con_tem;

					draw_contour(draw_, contour_dst, sh1, 8);
					FOR(m, 0, 4) circle(draw_, frame_b[m] + sh1, 3, Scalar(0, 255, 0));
					imshow("123123", draw_);
					con_re = set_flags(contour_dst, con_re);
					//Mat draw_ = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
					//vector<Point2f> frame = { contour_2[anc_mid[0]].point, contour_2[anc_mid[1]].point, contour_2[anc_mid[2]].point, contour_2[anc_mid[3]].point };
					//vector<Point2f> frame_b = base_frame(frame, 0);
					//Point2f sh1(0, 400);
					//draw_contour(draw_, contour_dst, sh1, 5);
					//FOR(m, 0, 4) circle(draw_, frame[m] + sh1, 3, Scalar(125, 0, 0));
					//Mat rot_mat = getPerspectiveTransform(frame, frame_b);// getAffineTransform(frame_1, frame_b_1);
					//perspectiveTransform(contour_dst, contour_dst, rot_mat);
					//draw_contour(draw_, contour_dst, sh1, 8);
					//FOR(m, 0, 4) circle(draw_, frame_b[m] + sh1, 3, Scalar(0, 255, 0));
					//imshow("123123", draw_);
				}
				/*contour_2 = set_flags(contour_dst, contour_2);
				Mat draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
				int first = 0, second = 0;
				if (self_intersect(contour_dst, first, second)) cout << "Bad  intersection!" << endl;
				if (translation_spec(contour_2, con_re, anc_mid, anc_re, draw2))
				{
					cout << "OK! No intersection!" << endl;
				}
				contour_dst = conf_trans(con_re);*/
				//con_re = set_flags(contour_dst, con_re);
				if (coll_opt)
				{
					//contour_de_crossing(contour_dst);
					contour_fine_tuning(contour_dst);
					con_re = set_flags(contour_dst, con_re);
				}
				if (deve_opt)
				{
					vector<Point2f> resam_dst = sampling_ave(contour_dst, contour_dst.size());
					vector<Point2f> resam_ = { contour_dst[anc_re[0]],contour_dst[anc_re[1]], contour_dst[anc_re[2]], contour_dst[anc_re[3]] };
					anc_re = relocate(resam_, resam_dst);
					FOR(ii, 0, 4) 
					{
						resam_dst[anc_re[ii]] = resam_[ii];
						cout << resam_[ii] << "   " << resam_dst[anc_re[ii]] << endl;
					}

					degree_after_opt = whole_con_opt(resam_dst, anc_re, 0);
					//degree_after_opt = whole_con_opt(contour_dst, anc_re, 0);
					//degree_after_opt = whole_con_opt(contour_dst, anc_mid, 0);
					cout << "After developable optimation, the collision degree: " << degree_after_opt << endl;
					vector<Point_f> con_resam;
					FOR(ii, 0, resam_dst.size()) con_resam.push_back(Point_f(resam_dst[ii], general_p));
					FOR(jj, 0, anc_re.size()) con_resam[anc_re[jj]].type = fixed_p;
					con_re = con_resam;
				}
				
				draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
				if (translation_spec(con_re, contour_2, anc_re, anc_mid, draw2))
				{
					cout << "OK! No intersection!" << endl;
				}			
				else cout << "Bad result with intersection!" << endl;
				protoTile c1, c2;
				c1.show_contour(conf_trans(con_re), anc_re);
				c2.show_contour(conf_trans(contour_2), anc_mid);
				RotationVis(c1, c2, AntiClockWise);

			}
		}


	}

	void Tiling_opt::load_dataset(bool input_images)
	{
		if (!contour_dataset.empty() || !all_con_tars.empty())
		{
			std::cout << "the Contours OR TARs have been exist" << endl;
			contour_dataset.swap(vector<vector<Point_f>>());
			all_con_tars.swap(vector<vector<vector<double>>>());
			all_con_tars_flip.swap(vector<vector<vector<double>>>());
			all_shape_complexity.swap(vector<double>());
		}
		if (input_images)
		{
			//vector<int> t = { 26 };
			//int tin = 0;
			FOR(i, 192, 193)
			{			
				//if (i != t[tin]) continue;
				//else tin++;
				string imagepath = DefaultPath;
				imagepath = imagepath + "dataset/" + to_string(i) + ".png";
				ima2contour(imagepath, input_images);
			}
		}
		FOR(i, 0, AllTypes)
		{
			string filepath = DefaultPath;
			filepath = filepath + "contour/" + to_string(i) + ".txt";
			cout << "Image " << to_string(i) << endl;
			protoTile read_tile(filepath, false);
			if (read_tile.contour.empty()) continue;
			contour_dataset.push_back(read_tile.contour_f);
			double shape_com;
			vector<vector<double>> tar_all = compute_TAR(read_tile.contour, shape_com);
			vector<vector<double>> tar_all_flip = compute_TAR(Flip_contour(read_tile.contour), shape_com);
			all_con_tars.push_back(tar_all);
			all_con_tars_flip.push_back(tar_all_flip);
			all_shape_complexity.push_back(shape_com);
		}
		all_types = contour_dataset.size();
		std::cout << "All types： " << all_types << "  load over!" << endl;
		std::cout << "All TARs of contours have been computed" << endl;
	}

	void Tiling_opt::load_para(std::string filename)
	{
		ifstream fin(filename);
		if (fin.fail() == true)
		{
			cout << "Cannot open load_pare file." << endl;
		}
		string param;
		while (fin >> param)
		{
			if (param == "image_id")
				fin >> image_id;
		}
		cout << "Paramaters Loaded Successfully." << endl << endl;
	}


	void Tiling_opt::RotationVis(protoTile c1, protoTile c2, int clockorder)
	{
		//       0_______3 0_______3  0________3
		//       |                |  |                 |   |                | 
		//       |      c2      |  |     c1       |   |      c3      | clockwise describes c1's rotating order in opencv
		//       |_______|  |_______|   |_______|
		//       1                2  1               2  1                2
		//
		cout << "c1: "<<c1.contour.size() << c1.contour[0] << " " << c1.contour[1] << " " << c1.contour[2] << " " << c1.contour[3]<<  endl;
		cout << "c2: "<<c2.contour.size() << c2.contour[0] << " " << c2.contour[1] << " " << c2.contour[2] << " " << c2.contour[3] << endl;
		int rot_deg = 46;
		int count = 10;
		while (count > 0)
		{
			Clock_order = clockorder;
			int draw_col = 1000;
			int draw_row = 1000;
			double scale_ratio;
			if (count == 10) scale_ratio = 0.65;
			else scale_ratio = 1.0;
			c1.scale(scale_ratio);
			c2.scale(scale_ratio);
			c1.contour = Trans_contour(c1.contour, Point2f(draw_row / 2, draw_col / 2) - center_p(c1.contour));
			c2.contour = Trans_contour(c2.contour, c1.contour[c1.anchor_points[0]] - c2.contour[c2.anchor_points[3]]);
			//c2.contour = Rotate_contour(c2.contour, center_p(c2.contour), -1);
			if (Clock_order == ClockWise)
			{
				cout << "ClockWise" << endl;
				for (int degree = 0; degree < rot_deg; degree += 2)
				{
					Mat drawing = Mat(draw_row, draw_col, CV_8UC3, Scalar(255, 255, 255));
					vector<Point2f> c1_r = Rotate_contour(c1.contour, c1.contour[c1.anchor_points[1]], -degree);
					vector<Point2f> c2_r = Rotate_contour(c2.contour, c2.contour[c2.anchor_points[2]], degree);
					vector<Point2f> c3_r = Trans_contour(c2_r, c1_r[c1.anchor_points[3]] - c2_r[c2.anchor_points[0]]);
					//from top to bottom, from left to right
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[3]] - c1_r[c1.anchor_points[2]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[0]] - c2_r[c2.anchor_points[1]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[3]] - c1_r[c1.anchor_points[2]], 0);
					draw_contour(drawing, c2_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, Point2f(0, 0), 0);
					draw_contour(drawing, c3_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[1]] - c1_r[c1.anchor_points[0]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[2]] - c2_r[c2.anchor_points[3]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[1]] - c1_r[c1.anchor_points[0]], 0);
					imshow("Rotation Visualization", drawing);
					if (degree == 0) waitKey(1000);
					else waitKey(200);
				}
			}
			else if (Clock_order == AntiClockWise)
			{
				cout << "AntiClockWise" << endl;
				for (int degree = 0; degree < rot_deg; degree += 2)
				{
					Mat drawing = Mat(draw_row, draw_col, CV_8UC3, Scalar(255, 255, 255));
					vector<Point2f> c1_r = Rotate_contour(c1.contour, c1.contour[c1.anchor_points[0]], degree);
					vector<Point2f> c2_r = Rotate_contour(c2.contour, c2.contour[c2.anchor_points[3]], -degree);
					vector<Point2f> c3_r = Trans_contour(c2_r, c1_r[c1.anchor_points[2]] - c2_r[c2.anchor_points[1]]);
					//from top to bottom, from left to right
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[0]] - c1_r[c1.anchor_points[1]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[3]] - c2_r[c2.anchor_points[2]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[0]] - c1_r[c1.anchor_points[1]], 0);
					draw_contour(drawing, c2_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, Point2f(0, 0), 0);
					draw_contour(drawing, c3_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[2]] - c1_r[c1.anchor_points[3]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[1]] - c2_r[c2.anchor_points[0]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[2]] - c1_r[c1.anchor_points[3]], 0);
					//draw_contour(drawing, c2_r, shift2, 3);
					imshow("Rotation Visualization", drawing);
					if (degree == 0) waitKey(1000);
					waitKey(200);
				}
			}
			count--;
		}
		//imwrite("D://Rotation Visualization.png", drawing);
	}


	int Tiling_opt::Tanslation_rule(vector<int> cand_points, vector<Point_f> &contour_s, string rootname)
	{
		int trans = 0;
		int ppindex = cand_points.size();
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(conf_trans(contour_s));
		//std::cout << "contsize: " << cent_cont << endl;
		int times = 0;
		for (int i = 0; i < ppindex; i++)
		{
			progress_bar(i, ppindex);
			for (int j = i + 1; j < ppindex; j++)
			{
				for (int m = j + 1; m < ppindex; m++)
				{
					for (int n = m + 1; n < ppindex; n++)
					{					
						times++;
						//std::cout << i<<" "<<j<<" "<<m<<" "<<n << endl;
						//if (abs(part_points_index[n] - part_points_index[m]) < margin) continue;
						vector<Point_f> inner_contour;
						vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
						vector<int> indexes; //1.translati on：以下所有该类摆放都以1-3,2-4为轴摆放
						indexes.push_back(cand_points[i]);
						indexes.push_back(cand_points[j]);
						indexes.push_back(cand_points[m]);
						indexes.push_back(cand_points[n]);
						Mat drawing1 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));

						if (translation_placement( contour_s, inner_contour, indexes, mid_interval, drawing1))
						{
							//std::cout << ++trans << " Translation succeed" << endl;
							inPat one_situation(inner_contour, mid_interval, 0);
							all_inner_conts.push_back(one_situation);

							Point2f shift2 = Point2f(300, 600) - cent_cont;
							FOR(jj, 0, contsize)
							{
								if (contour_s[jj].type == general_p)  circle(drawing1, contour_s[jj].point+ shift2, 2, Scalar(0, 0, 0), -1);
								else if (contour_s[jj].type != general_p)  circle(drawing1, contour_s[jj].point + shift2, 2, Scalar(0, 255, 0), -1);
							}
							FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
							circle(drawing1, contour_s[0].point + shift2, 4, Scalar(255, 0, 55), -1);

							string filename = rootname + "/" + to_string(all_inner_conts.size() - 1) + "transPlacingResult.png";
							cv::imwrite(filename, drawing1);
						}
					}
				}
			}
		}
		std::cout << "Trans times: " << times << endl;
		return trans;
	}

	int Tiling_opt::Tanslation_rule_spec(vector<int> cand_points, vector<Point_f> &contour_s, string rootname, vector<int> anc_points)
	{
		int trans = 0;
		int ppindex = cand_points.size();
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(conf_trans(contour_s));
		//std::cout << "contsize: " << cent_cont << endl;
		int times = 0;
		int i = anc_points[0];
		int j = anc_points[1];
		int m = anc_points[2];
		int n = anc_points[3];
		times++;
		//std::cout << i<<" "<<j<<" "<<m<<" "<<n << endl;
		//if (abs(part_points_index[n] - part_points_index[m]) < margin) continue;
		vector<Point_f> inner_contour;
		vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
		vector<int> indexes; //1.translati on：以下所有该类摆放都以1-3,2-4为轴摆放
		indexes.push_back(cand_points[i]);
		indexes.push_back(cand_points[j]);
		indexes.push_back(cand_points[m]);
		indexes.push_back(cand_points[n]);
		Mat drawing1 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));

		if (translation_placement(contour_s, inner_contour, indexes, mid_interval, drawing1))
		{
			//std::cout << ++trans << " Translation succeed" << endl;
			inPat one_situation(inner_contour, mid_interval, 0);
			all_inner_conts.push_back(one_situation);

			Point2f shift2 = Point2f(300, 600) - cent_cont;
			FOR(jj, 0, contsize)
			{
				if (contour_s[jj].type == general_p)  circle(drawing1, contour_s[jj].point + shift2, 2, Scalar(0, 0, 0), -1);
				else if (contour_s[jj].type != general_p)  circle(drawing1, contour_s[jj].point + shift2, 2, Scalar(0, 255, 0), -1);
			}
			FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
			circle(drawing1, contour_s[0].point + shift2, 4, Scalar(255, 0, 55), -1);

			string filename = rootname + "/" + to_string(all_inner_conts.size() - 1) + "transPlacingResult.png";
			cv::imwrite(filename, drawing1);
		}
		std::cout << "Trans times: " << times << endl;
		return trans;
	}


	int Tiling_opt::Flipping_rule(vector<int> part_points_index, vector<Point_f> &contour_s, string rootname)
	{
		return 0;
	}

	bool Tiling_opt::translation_placement(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname)
	{
		bool test_coll = false;
		int csize = contour_s.size();
		Point2f line1 = contour_s[indexes[2]].point - contour_s[indexes[0]].point;
		Point2f line2 = contour_s[indexes[3]].point - contour_s[indexes[1]].point;
		vector<Point2f> shifting;
		shifting.push_back(line1);
		shifting.push_back(line2);
		shifting.push_back(line1 + line2);
		// 提取平移围成的轮廓
		vector<vector<Point_f>> four_place;
		four_place.push_back(contour_s);
		for (int i = 0; i < 3; i++)
		{
			vector<Point_f> one_loca;
			for (int j = 0; j < csize; j++)
			{
				one_loca.push_back(Point_f(contour_s[j].point + shifting[i], contour_s[j].type));
			}
			four_place.push_back(one_loca);
		}
		//std::cout << "compute four_place" << endl;
		//translation的提取规律为1的4-3,2的1-4,4的2-1,3的3-2
		//              ①     ③            共有①②③④四个图案，
		//               /\      /\			   检测①中的顶点3与②中顶点1，
		//             /1 \   /1 \           以及①中顶点4与③中顶点2的角度之和
		//            /      \/      \         
		//          2\   4/\2    /4
		//              \  /    \  /
		//             3\/      \/3
		//               /\      /\     
		//             /1 \   /1 \ 
		//            /     \/       \          为保证在一定范围内可靠，
		//          2\   4/\2  4/          比较时计算∠Pt-iPtPt+i的角度，
		//              \  /   \  /           i从1取到3
		//              3\/    \/3
		//               ②    ④
		if (test_coll)
		{
			Mat drawing1 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
			Point2f cent_cont = center_p(conf_trans(contour_s));
			Point2f shift2 = Point2f(400, 600) - cent_cont;
			draw_contour_points(drawing1, conf_trans(contour_s), shift2);
			FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
			Point2f shift1 = Point2f(1000, 600) - (center_p(conf_trans(contour_s)) + 0.5*(line1 + line2));
			FOR(i, 0, 4) draw_poly(drawing1, conf_trans(four_place[i]), shift1);
			string spath = "[" + to_string(indexes[0]) + "," + to_string(indexes[1]) + "," + to_string(indexes[2]) + "," + to_string(indexes[3]) + "]_placement";
			imshow(spath, drawing1);
		}

		if (coll_detec(conf_trans(four_place[0]), conf_trans(four_place[1]), 1) || coll_detec(conf_trans(four_place[0]), conf_trans(four_place[2]), 1)
			|| coll_detec(conf_trans(four_place[0]), conf_trans(four_place[3]), 0) || coll_detec(conf_trans(four_place[1]), conf_trans(four_place[2]), 0))
		{
			//std::cout << "polygon collision" << endl;
			return false;
		}
		int total_num = 0;
		ex_indexes.swap(vector<int>());
		extracted.swap(vector<Point_f>());
		ex_indexes.push_back(0);
		for (int t = indexes[3]; t > indexes[2]; t--)                   //translation的提取规律为1的4-3,2的1-4,4的2-1,3的3-2
		{
			extracted.push_back(four_place[0][t]);
			total_num++;
			//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
		}
		ex_indexes.push_back(total_num);
		for (int t = indexes[0] + csize; t > indexes[3]; t--)
		{
			extracted.push_back(four_place[1][t % csize]);
			total_num++;
		}
		ex_indexes.push_back(total_num);
		for (int t = indexes[1]; t > indexes[0]; t--)
		{
			extracted.push_back(four_place[3][t]);
			total_num++;
		}
		ex_indexes.push_back(total_num);
		for (int t = indexes[2]; t > indexes[1]; t--)
		{
			extracted.push_back(four_place[2][t]);
		}
		//提取后检测是否自交
		int index_1, index_2;
		bool self_inter = self_intersect(conf_trans(extracted), index_1, index_2);
		Point2f shift1 = Point2f(1100, 600) - (center_p(conf_trans(contour_s)) + 0.5*(line1 + line2));
		if (!self_inter)
		{
			FOR(i, 0, 4) draw_poly(countname, conf_trans(four_place[i]), shift1);
			draw_poly(countname, conf_trans(extracted), shift1, 5);
		}
		return !self_inter;
	}

	bool Tiling_opt::translation_spec(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname)
	{
		int csize = contour_s.size();
		Point2f line1 = contour_s[indexes[2]].point - contour_s[indexes[0]].point;
		Point2f line2 = contour_s[indexes[3]].point - contour_s[indexes[1]].point;
		vector<Point2f> shifting;
		shifting.push_back(line1);
		shifting.push_back(line1 + line2);
		shifting.push_back(line2);
		
		// 提取平移围成的轮廓
		vector<vector<Point_f>> four_place;
		four_place.push_back(contour_s);
		for (int i = 0; i < 3; i++)
		{
			vector<Point_f> one_loca;
			for (int j = 0; j < csize; j++)
			{
				one_loca.push_back(Point_f(contour_s[j].point + shifting[i], contour_s[j].type));
			}
			four_place.push_back(one_loca);
		}

		Point2f cent_cont = center_p(conf_trans(contour_s));
		Point2f shift2 = Point2f(400, 600) - cent_cont;
		draw_contour_points(countname, conf_trans(contour_s), shift2);
		FOR(jj, 0, 4) circle(countname, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
		Point2f shift1 = Point2f(1000, 600) - (center_p(conf_trans(contour_s)) + 0.5*(line1 + line2));
		FOR(i, 0, 4) draw_poly(countname, conf_trans(four_place[i]), shift1);
		string spath = "[" + to_string(indexes[0]) + "," + to_string(indexes[1]) + "," + to_string(indexes[2]) + "," + to_string(indexes[3]) + "]_placement";

		int total_num = 0;
		ex_indexes.swap(vector<int>());
		extracted.swap(vector<Point_f>());
		//ex_indexes.push_back(0);
		for (int i = 0; i < four_place.size(); i++)
		{
			ex_indexes.push_back(total_num);
			int s_i = (i + 3) % 4;
			int e_i = (i + 2) % 4;
			int t = indexes[s_i];
			if (s_i<e_i)  t += csize;
			for (; t > indexes[e_i]; t--)                   //translation的提取规律为1的4-3, 2的1-4, 3的2-1, 4的3-2
			{
				extracted.push_back(four_place[i][t % csize]);
				total_num++;
				//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
			}
		}

		//for (int t = indexes[3]; t > indexes[2]; t--)                   //translation的提取规律为1的4-3, 2的1-4, 3的2-1, 4的3-2
		//{
		//	extracted.push_back(four_place[0][t]);
		//	total_num++;
		//	//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
		//}
		//ex_indexes.push_back(total_num);
		//for (int t = indexes[0] + csize; t > indexes[3]; t--)
		//{
		//	extracted.push_back(four_place[1][t % csize]);
		//	total_num++;
		//}
		//ex_indexes.push_back(total_num);
		//for (int t = indexes[1]; t > indexes[0]; t--)
		//{
		//	extracted.push_back(four_place[3][t]);
		//	total_num++;
		//}
		//ex_indexes.push_back(total_num);
		//for (int t = indexes[2]; t > indexes[1]; t--)
		//{
		//	extracted.push_back(four_place[2][t]);
		//}
		draw_poly(countname, conf_trans(extracted), shift1, 5);
		imshow(spath, countname);
		//提取后检测是否自交
		int index_1, index_2;
		bool self_inter = self_intersect(conf_trans(extracted), index_1, index_2);
		return !self_inter;
	}

	void Tiling_opt::match_candidate(int inner_index)
	{
		inPat inner= all_inner_conts[inner_index];
		prototile_mid = protoTile(inner.in_contour);
		int cand_num = 10;
		vector<pair<int, bool>> candidate_patterns = compare_TAR(prototile_mid.contour_f, cand_num);//当前中间图案对应的候选图案
		cout << "This is the " << inner_index << "th inner contour" << endl <<
			"candidate_patterns size: " << candidate_patterns.size() << endl;

		double sc_inner = 0;
		vector<vector<double>> inner_tar = compute_TAR(prototile_mid.contour, sc_inner);

		candidate_contours.swap(vector<vector<Point_f>>());
		cand_paths.swap(vector<vector<pair<int, int>>>());

		for (int j = 0; j < cand_num; j++) //只要候选图案里的前cand_num个
		{
			//将所有的结果保存下来
			prototile_second = protoTile(contour_dataset[candidate_patterns[j].first]);
			vector<vector<double>> cand_tar;
			if (candidate_patterns[j].second)
			{
				//std::cout << "it is flip" << endl;
				cand_tar = all_con_tars_flip[candidate_patterns[j].first];
				prototile_second.Flip_contour(prototile_second.contour_f);
				prototile_second.contour = conf_trans(prototile_second.contour_f);
			}
			else
			{
				//std::cout << "it is not flip" << endl;
				cand_tar = all_con_tars[candidate_patterns[j].first];
			}
			int cmsize = prototile_mid.contour_f.size();
			vector<int> mid_index(cmsize, 0), sec_f(cmsize, 0);
			FOR(n, 0, cmsize)
			{
				mid_index[n] = prototile_mid.contour_f[n].type;
				sec_f[n] = prototile_second.contour_f[n].type;
			}
			vector<pair<int, int>> path;
			int shift = 0;
			int width = match_window_width; //
			double re = tar_mismatch_fea(inner_tar, cand_tar, mid_index, sec_f, path, shift, width);
			//double re = tar_mismatch(inner_tar, cand_tar, path, shift, width);
			vector<pair<int, int>> path_min = path;
			/*for (int n = 0; n < path.size(); n++)
			if (prototile_mid->contour_f[path[n].first].type>1)
			path_min.push_back(path[n]);*/

			vector<Point_f> contour_tem_f;
			int c2size = prototile_second.contour_f.size();
			for (int i = shift; i < shift + c2size; i++)
			{
				contour_tem_f.push_back(prototile_second.contour_f[i % c2size]);
			}
			prototile_second.contour_f = contour_tem_f;
			vector<Point2f> contour_cand = conf_trans(prototile_second.contour_f);
			//std::cout << "Mid num: "<<prototile_mid.contour.size() << "   sec num: " << c2size << "  after shifting: " << contour_cand.size() << endl;

			double scale = arcLength(prototile_mid.contour, true) / arcLength(contour_cand, true);
			//std::cout << "scale: " << scale << endl;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] * scale;
			}
			Point2f cen1 = center_p(prototile_mid.contour);
			Point2f cen2 = center_p(contour_cand);
			Point2f shift2 = cen1 - cen2;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] + shift2;
			}
			double length_min = 1000000;
			int angle_min = 0;
			int times = 0;
			double angl_al = 0;
			Point2f sh_al = Point2f(0, 0);
			Point2f shift_t = Point2f(1000, 1000);
			vector<Point2f> contour_tem;
			//Point2f v0 = contour_inner[0] - cen1;
			//Point2f v1 = contour_cand[0] - cen1;
			//double angle1 = acos(cos_two_vector(v0, v1)) / PI * 180;
			//if (sin_two_vector(v0, v1) < 0) angle1 = -angle1;

			//cout << j << "  path_min_size:  " << path_min.size() << "  cent:" << cen2 << endl;
			while (times < 3 || length_2p(shift_t, Point2f(0, 0)) > 100 || angle_min > 10)
			{
				length_min = 1000000;
				angle_min = 0;
				cen2 = center_p(contour_cand);
				//找到距离最近的角度
				for (int angle = 0; angle < 360; angle = angle + 5)
				{
					double leng = 0;
					Mat rot_mat(2, 3, CV_32FC1);
					rot_mat = getRotationMatrix2D(cen2, angle, 1);
					cv::transform(contour_cand, contour_tem, rot_mat);
					for (int m = 0; m < path_min.size(); m++)
					{
						leng += length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]);
					}
					if (leng < length_min)
					{
						angle_min = angle;
						length_min = leng;
					}
				}
				Mat rot_mat1(2, 3, CV_32FC1);
				rot_mat1 = getRotationMatrix2D(cen2, angle_min, 1);
				cv::transform(contour_cand, contour_tem, rot_mat1);
				angl_al += angle_min;
				//contour_cand = contour_mid;
				//移动重心
				shift_t = Point2f(0, 0);
				for (int m = 0; m < path_min.size(); m++)
				{
					shift_t += prototile_mid.contour[path_min[m].first] - contour_tem[path_min[m].second];
				}
				shift_t = Point2f(shift_t.x / path_min.size(), shift_t.y / path_min.size());
				for (int i = 0; i < contour_cand.size(); i++)
				{
					contour_cand[i] = contour_tem[i] + shift_t;
				}
				sh_al += shift_t;
				times++;
			}
			//cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << endl;
			prototile_second.contour = contour_cand;
			candidate_contours.push_back(set_flags(contour_cand, prototile_second.contour_f));
			cand_paths.push_back(path);
		}
	}

	void Tiling_opt::match_candidate_(int inner_index)
	{
		inPat inner = all_inner_conts[inner_index];
		prototile_mid = protoTile(inner.in_contour);
		int cand_num = 10;
		vector<pair<int, bool>> candidate_patterns = compare_TAR(prototile_mid.contour_f, cand_num);//当前中间图案对应的候选图案
		cout << "This is the " << inner_index << "th inner contour" << endl <<
			"candidate_patterns size: " << candidate_patterns.size() << endl;

		double sc_inner = 0;
		vector<vector<double>> inner_tar = compute_TAR(prototile_mid.contour, sc_inner);
		//vector<vector<Point_f>> candidate_contour;
		//vector<vector<pair<int, int>>> cand_path;
		candidate_contours.swap(vector<vector<Point_f>>());
		cand_paths.swap(vector<vector<pair<int, int>>>());
		cand_fea_paths.swap(vector<vector<pair<int, int>>>());
		path_shift.swap(vector<int>());

		for (int j = 0; j < 1; j++) //只要候选图案里的前cand_num个
		{
			//将所有的结果保存下来
			prototile_second = protoTile(contour_dataset[candidate_patterns[j].first]);
			vector<vector<double>> cand_tar;
			if (candidate_patterns[j].second)
			{
				//std::cout << "it is flip" << endl;
				cand_tar = all_con_tars_flip[candidate_patterns[j].first];
				prototile_second.Flip_contour(prototile_second.contour_f);
				prototile_second.contour = conf_trans(prototile_second.contour_f);
				prototile_second.feature_points = prototile_second.getFeatures(prototile_second.contour_f);
			}
			else
			{
				//std::cout << "it is not flip" << endl;
				cand_tar = all_con_tars[candidate_patterns[j].first];
			}
			int cmsize = prototile_mid.contour_f.size();
			vector<int> mid_index(cmsize, 0), sec_f(cmsize, 0);
			FOR(n, 0, cmsize)
			{
				mid_index[n] = prototile_mid.contour_f[n].type;
				sec_f[n] = prototile_second.contour_f[n].type;
			}
			vector<pair<int, int>> path;
			int shift = 0;
			int width = match_window_width; //
			double re = tar_mismatch_fea(inner_tar, cand_tar, mid_index, sec_f, path, shift, width); 
			vector<pair<int, int>> path_min = path;
			vector<Point_f> contour_tem_f;
			int c2size = prototile_second.contour_f.size();
			//移动点顺序
			for (int i = shift; i < shift + c2size; i++)
			{
				contour_tem_f.push_back(prototile_second.contour_f[i % c2size]);
			}
			prototile_second.contour_f = contour_tem_f;

			vector<Point2f> contour_cand = conf_trans(prototile_second.contour_f);
			//std::cout << "Mid num: "<<prototile_mid.contour.size() << "   sec num: " << c2size << "  after shifting: " << contour_cand.size() << endl;

			//移动到重心重叠的地方
			Point2f cen1 = center_p(prototile_mid.contour);
			Point2f cen2 = center_p(contour_cand);
			Point2f shift2 = cen1 - cen2;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] + shift2;
			}
			//two types of calculating the scale: i.based on the arclength  ii.based in uniform scale
			double scale = arcLength(prototile_mid.contour, true) / arcLength(contour_cand, true);
			//std::cout << "scale: " << scale << endl;
			/*double s_1 = 0, s_2 = 0;
			FOR(mn, 0, contour_cand.size())
			{
				s_1 += pow(length_2p(prototile_mid.contour[mn], cen1), 2);
				s_2 += pow(length_2p(contour_cand[mn], cen1), 2);
			}
			double scale2 = sqrt(s_1 / s_2);
			cout << "scale2: " << scale2 << endl;*/
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] * scale;
			}

			double length_min = 1000000;
			int angle_min = 0;
			int times = 0;
			double angl_al = 0;
			Point2f sh_al = Point2f(0, 0);
			Point2f shift_t = Point2f(1000, 1000);
			vector<Point2f> contour_tem;
			//Point2f v0 = contour_inner[0] - cen1;
			//Point2f v1 = contour_cand[0] - cen1;
			//double angle1 = acos(cos_two_vector(v0, v1)) / PI * 180;
			//if (sin_two_vector(v0, v1) < 0) angle1 = -angle1;

			//cout << j << "  path_min_size:  " << path_min.size() << "  cent:" << cen2 << endl;
			while (times < 3 || length_2p(shift_t, Point2f(0, 0)) > 10 || angle_min > 10)// 
			{
				length_min = 1000000;
				angle_min = 0;
				cen2 = center_p(contour_cand);
				//找到距离最近的角度
				for (int angle = 0; angle < 360; angle = angle + 2)
				{
					double leng = 0;
					Mat rot_mat(2, 3, CV_32FC1);
					rot_mat = getRotationMatrix2D(cen2, angle, 1);
					cv::transform(contour_cand, contour_tem, rot_mat);
					for (int m = 0; m < path_min.size(); m++)
					{
						leng += pow(length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]),2);
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
				cv::transform(contour_cand, contour_tem, rot_mat1);
				angl_al += angle_min;
				//contour_cand = contour_tem;
				//移动重心
				shift_t = Point2f(0, 0);
				for (int m = 0; m < path_min.size(); m++)
				{
					shift_t += prototile_mid.contour[path_min[m].first] - contour_tem[path_min[m].second];
				}
				shift_t = Point2f(shift_t.x / path_min.size(), shift_t.y / path_min.size());
				for (int i = 0; i < contour_cand.size(); i++)
				{
					contour_cand[i] = contour_tem[i] + shift_t;
				}
				sh_al += shift_t;
				times++;
				cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << "   shift_t: "<< shift_t<<endl;
			}
			//cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << endl;
			prototile_second.contour = contour_cand;
			prototile_second.contour_f = set_flags(contour_cand, prototile_second.contour_f);
			candidate_contours.push_back(prototile_second.contour_f);
			cand_paths.push_back(path);

			//filter out inappropriate matches
			vector<pair<int, int>> path_fea;
			int path_size = path.size();
			int sec_size = prototile_second.contour_f.size();
			cout << path_size << "   " << sec_size << endl;
			FOR(m, 0, path_size)
			{
				int mar = 4;
				//cout << "m: " << m << "   " << path[m].first<< "  " << path[m].second << endl;
				//cout << "m: " <<m<< "   " << prototile_mid.contour_f[path[m].first].type <<"  "<< prototile_second.contour_f[path[m].second].type << endl;
				pair<int, int> one_pair = path[m];
				int sec_index = one_pair.second;
				if (prototile_mid.contour_f[one_pair.first].type == fixed_p || prototile_mid.contour_f[one_pair.first].type == fea_p)  //must have a match
				{
					vector<int> waiting_merge;
					FOR(n, 1-mar, mar)
					{
						int tem_sec = (sec_index + n + sec_size) % sec_size;
						if (prototile_second.contour_f[tem_sec].type == fea_p)
							waiting_merge.push_back(tem_sec);
					}
					//if (m == 51) cout << waiting_merge.size()<<"    "<<waiting_merge[0] <<"   "<<shift << endl;
					if (waiting_merge.empty())
					{
						if(prototile_mid.contour_f[one_pair.first].type == fixed_p)
							path_fea.push_back(one_pair);
						else continue;
					}
					else
					{
						int min_index = 0;
						if (waiting_merge.size() != 1)
						{
							double min_dis_tar = 10000;
							FOR(t, 0, waiting_merge.size())
							{
								//calculate the pair with less tar distance
								int cand_index = (waiting_merge[t] + shift) % cand_tar.size();
								int pcsize = prototile_mid.contour.size();
								double dis_cos1 = cos_2v(prototile_mid.contour[(one_pair.first - 1 + pcsize) % pcsize] - prototile_mid.contour[one_pair.first], prototile_mid.contour[(one_pair.first + 1) % pcsize] - prototile_mid.contour[one_pair.first]);
								double dis_cos2 = cos_2v(prototile_second.contour[(waiting_merge[t] - 1 + pcsize) % pcsize] - prototile_second.contour[waiting_merge[t]], prototile_second.contour[(waiting_merge[t] + 1) % pcsize] - prototile_second.contour[waiting_merge[t]]);
								//cout << cand_index<<"   "<<prototile_second.contour[(waiting_merge[t] - 1 + pcsize) % pcsize] << "   " << prototile_second.contour[waiting_merge[t]] << "   " << prototile_second.contour[(waiting_merge[t] + 1) % pcsize] << endl;
								//cout << cand_index << "   " << prototile_second.contour[cand_index - 1 ] << "   " << prototile_second.contour[cand_index] << "   " << prototile_second.contour[cand_index + 1] << endl;
								double dis_tar = tar_length_2p(inner_tar[one_pair.first], cand_tar[cand_index]) + 0.02*abs(dis_cos1 - dis_cos2);
								//cout << "distar" << dis_tar <<"  dis_cos1:"<< dis_cos1<<"   "<< dis_cos2<< endl;
								if (dis_tar < min_dis_tar)
								{
									min_dis_tar = dis_tar;
									min_index = t;
								}
							}
						}
						cout << "m: " << m << "   " << min_index << endl;
						int repet_index = -1;
						FOR(pin, 0, path_fea.size())
						{
							if (waiting_merge[min_index] == path_fea[pin].second)
							{
								repet_index = pin;
							}
						}
						//cout<<
						if (path_fea.empty() || repet_index == -1) 
							path_fea.push_back(make_pair(one_pair.first, waiting_merge[min_index]));
						else if (prototile_mid.contour_f[one_pair.first].type>prototile_mid.contour_f[path_fea[repet_index].first].type)
						{
							path_fea.pop_back();
							path_fea.push_back(make_pair(one_pair.first, waiting_merge[min_index]));
						}
						else if (prototile_mid.contour_f[one_pair.first].type == prototile_mid.contour_f[path_fea[repet_index].first].type)
						{
							int cand_index = (waiting_merge[min_index] + shift) % cand_tar.size();
							double one_tar = tar_length_2p(inner_tar[one_pair.first], cand_tar[cand_index]);
							double back_tar = tar_length_2p(inner_tar[path_fea[repet_index].first], cand_tar[cand_index]);
							if (one_tar < back_tar)
							{
								path_fea.pop_back();
								path_fea.push_back(make_pair(one_pair.first, waiting_merge[min_index]));
							}
						}
					}				
				}			
			}
			FOR(mm, 0, path_fea.size()) cout << path_fea[mm].first<<"   "<< path_fea[mm].second<< endl;
			cand_fea_paths.push_back(path_fea);
		}
	}

	vector<pair<int, bool>> Tiling_opt::compare_TAR(vector<Point_f> contour_mid, int chosen_num, int window_width) //chosen_num  选择的最终结果的数目
	{
		//int all_types = 2;
		//std::cout << "inner_c num:" << inner_c.size() << endl;
		vector<Point2f> contour_mid_ = conf_trans(contour_mid);
		double iso_mid = isoperimetric_inequality(contour_mid_);
		cout << "Compare with " << all_types << " types candadites." << endl;

		vector<pair<int, bool>> all_final;
		vector<double> all_result;
		double shape_com_mid;
		vector<vector<double>> tar_mid = compute_TAR(contour_mid_, shape_com_mid);
		std::cout << "contour_mid: " << contour_mid.size() << "  tar_mid: " << tar_mid.size() << "  each tar: "<< tar_mid[0].size()<<endl;
		//FOR(g, 0, tar_mid[0].size()) cout << tar_mid[0][g] << endl;
		int ttrt = 0;
		for (int index = 0; index < all_types; index++)
		{
			//first filter by isoperimetric_inequality
			double iso_cand = isoperimetric_inequality(conf_trans(contour_dataset[index]));
			double diss = abs(iso_mid - iso_cand) / iso_mid;
			if (diss > SCDisThres)
			{
				ttrt++;
				continue;
			}
			vector<Point_f> sec = contour_dataset[index];
			vector<vector<double>> tar_sec = all_con_tars[index];
			vector<vector<double>> tar_sec_f = all_con_tars_flip[index];
			vector<pair<int, int>> path;
			int shift = 0;
			double dis_sc = abs(shape_com_mid - all_shape_complexity[index]);
			int cmsize = contour_mid.size();
			vector<int> mid_index(cmsize,0), sec_f(cmsize,0), sec_fflip(cmsize, 0);
			FOR(n, 0, cmsize)
			{
				mid_index[n]= contour_mid[n].type;
				sec_f[n] = sec[n].type;
				sec_fflip[n] = sec[cmsize - 1 - n].type;
				//cout << n << "   " << mid_index[n] << "    " << sec_f[n] << "     " << sec_fflip[n] << endl;
			}
			double re = tar_mismatch_fea(tar_mid, tar_sec, mid_index, sec_f, path, shift, window_width);
			//double re = tar_mismatch(tar_mid, tar_sec, path, shift, window_width);
			//re = re / (1 + shape_com_mid + all_shape_complexity[index]);  //ori
			re = re / (1 + shape_com_mid + all_shape_complexity[index] - 5 * dis_sc - 2 * diss);
			re = 1 - 0.5*re / path.size();
			//cout << "re: " << re << endl;
			double re2 = tar_mismatch_fea(tar_mid, tar_sec_f, mid_index, sec_fflip, path, shift, window_width);
			//double re2 = tar_mismatch(tar_mid, tar_sec_f, path, shift, window_width);
			//re2 = re2 / (1 + shape_com_mid + all_shape_complexity[index]);
			re2 = re2 / (1 + shape_com_mid + all_shape_complexity[index] - 5 * dis_sc - 2 * diss);
			re2 = 1 - 0.5*re2 / path.size();
			//cout << "re2: " << re2 << endl;
			if (re > re2)
			{
				all_result.push_back(re);
				all_final.push_back(make_pair(index, false));
				//std::cout << re << endl;
			}
			else
			{
				all_result.push_back(re2);
				all_final.push_back(make_pair(index, true));
				//std::cout << re2 << endl;
			}
			//cout << "cand: " << index << "  " << all_result.back() << "  diss:" << diss << endl;
			//system("pause");
		}
		//sort
		double temp;
		pair<int, bool> tempp;
		int all_size = all_result.size();
		for (int i = 0; i < all_size - 1; i++)
			for (int j = 0; j < all_size - 1 - i; j++)
				if (all_result[j] < all_result[j + 1])
				{
					temp = all_result[j];
					all_result[j] = all_result[j + 1];
					all_result[j + 1] = temp;
					tempp = all_final[j];
					all_final[j] = all_final[j + 1];
					all_final[j + 1] = tempp;
				}
		//std::cout << "the fianl order: " << endl;
		vector<pair<int, bool>> all_total_mid;
		int asize = all_final.size();
		chosen_num = min(chosen_num, asize);
		for (int t = 0; t < chosen_num; t++)
		{
			all_total_mid.push_back(all_final[t]);
			std::cout << "order: " << all_final[t].first << "  flip: " << all_final[t].second << " value: " << all_result[t] << " complxeity: " << all_shape_complexity[all_final[t].first] << "     mid sc:" << shape_com_mid << endl;
		}
		cout << "The num of ones beyond ShapeComplexity DisThres: " << ttrt << endl;
		return all_total_mid;
	}

	void Tiling_opt::feature_match(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<int> first_fea, vector<int> sec_fea, vector<pair<int, int>>& path, int &sec_shift, int width)
	{
		vector<vector<double>> first_fea_arr;
		vector<vector<double>> second_fea_arr;
		int fsize = first_fea.size();
		int ssize = sec_fea.size();
		FOR(i, 0, fsize)
		{
			first_fea_arr.push_back(first_arr[first_fea[i]]);
		}
		FOR(i, 0, ssize)
		{
			second_fea_arr.push_back(second_arr[sec_fea[i]]);
		}
		double re = tar_mismatch(first_fea_arr, second_fea_arr, path, sec_shift, width);
		vector<int> sec_fea2;
		vector<vector<double>> second_fea_arr2;
		FOR(i, 0, ssize)
		{
			sec_fea2.push_back(sec_fea[(i + sec_shift)% ssize]);
			second_fea_arr2.push_back(second_fea_arr[(i + sec_shift) % ssize]);
		}
		sec_fea = sec_fea2;
		vector<pair<int, int>> waiting_merge;
		vector<pair<int, int>> path_new;
		FOR(i, 0, path.size())
		{
			if(waiting_merge.empty()) 
				waiting_merge.push_back(path[i]);
			else
			{
				if (path[i].first == waiting_merge[0].first || path[i].second == waiting_merge[0].second)
				{
					waiting_merge.push_back(path[i]);
				}
				else
				{
					int min_index = 0;
					double min_dis_tar = 10000;
					FOR(j, 0, waiting_merge.size())
					{
						//calculate the pair with less tar distance
						double dis_tar= tar_length_2p(first_fea_arr[waiting_merge[j].first], second_fea_arr2[waiting_merge[j].second]);
						if (dis_tar < min_dis_tar)
						{
							min_dis_tar = dis_tar;
							min_index = j;
						}
					}
					path_new.push_back(waiting_merge[min_index]);
					waiting_merge.swap(vector<pair<int, int>>());
					waiting_merge.push_back(path[i]);
				}
			}
		}
		//处理waiting_merge中剩下的pair
		int min_index = 0;
		double min_dis_tar = 10000;
		FOR(j, 0, waiting_merge.size())
		{
			//calculate the pair with less tar distance
			double dis_tar = tar_length_2p(first_fea_arr[waiting_merge[j].first], second_fea_arr2[waiting_merge[j].second]);
			if (dis_tar < min_dis_tar)
			{
				min_dis_tar = dis_tar;
				min_index = j;
			}
		}
		path_new.push_back(waiting_merge[min_index]);
		path = path_new;
	}

	vector<Point_f> Tiling_opt::morphing(vector<Point_f> contour1, vector<Point_f> contour2, vector<pair<int, int>> final_pair, double ratio)
	{
		int fpsize = final_pair.size();
		int cnum1 = contour1.size();
		int cnum2 = contour2.size();
		Point2f cen1 = center_p(conf_trans(contour1));
		cout << "Matched feature point pair : " << fpsize << endl;
		//将两个轮廓分段存储
		vector<vector<Point_f>> contour1_seg;
		vector<vector<Point_f>> contour2_seg;
		vector<Point_f> final_con;
		for (int g = 0; g < fpsize; g++)
		{
			vector<Point_f> c_seg1;
			vector<Point_f> c_seg2;
			//cout << final_pair[g].first << "  :  " << final_pair[g].second << endl;
			int s_index1 = final_pair[g].first;
			int e_index1 = final_pair[(g + 1) % fpsize].first;
			e_index1 > s_index1 ? e_index1 : e_index1 += cnum1;
			//e_index1 = e_index1 + (e_index1 > s_index1 ? 0 : cnum1);
			int s_index2 = final_pair[g].second;
			int e_index2 = final_pair[(g + 1) % fpsize].second;

			e_index2> s_index2 ? e_index2 : e_index2 += cnum2;
			for (int n = s_index1; n <= e_index1; n++)
				c_seg1.push_back(contour1[n%cnum1]);
			for (int n = s_index2; n <= e_index2; n++)
				c_seg2.push_back(contour2[n%cnum2]);
			contour1_seg.push_back(c_seg1);
			contour2_seg.push_back(c_seg2);
		}

		//cout << "contour1_seg:  " << contour1_seg.size() << "    contour2_seg:  " << contour2_seg.size() << endl;
		//	<< contour1_seg.back().back().point << "   " << contour1_seg.back().back().type << endl;
		Mat draw_seg = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		Point2f shift_seg = Point2f(0.5*draw_seg.cols, 0.5*draw_seg.rows) - cen1;
		Point_f start_new = Point_f(ratio*contour1_seg[0][0].point + (1 - ratio)*contour2_seg[0][0].point, contour1_seg[0][0].type);
		double num_e = 0;
		FOR(seg_index, 0, fpsize)
		{
			if (seg_index == fpsize - 1) contour1_seg[seg_index].back().type = 4;
			vector<Point_f> morph_seg = morph_segment(contour1_seg[seg_index], contour2_seg[seg_index], start_new, ratio, num_e);
			final_con.insert(final_con.end(), morph_seg.begin(), morph_seg.end() - 1);
			start_new = morph_seg.back();

			/*if (seg_index == 10)
			{
				FOR(TT, 0, contour1_seg[seg_index].size())  cout << TT << "   " << contour1_seg[seg_index][TT].point << endl;
				FOR(TT, 0, contour2_seg[seg_index].size())  cout << TT << "   " << contour2_seg[seg_index][TT].point << endl;
			}*/
				//FOR(TT, 0, morph_seg.size())  cout << TT << "   " << morph_seg[TT].point << endl;
			line(draw_seg, morph_seg[0].point + shift_seg, morph_seg.back().point + shift_seg, Scalar(50, 225, 50), 1);
			cout << contour1_seg[seg_index].size()<<"  "<< contour2_seg[seg_index].size()<<"   "<<morph_seg.size()<<"   "<< morph_seg[0].point << "   " << morph_seg.back().point << endl;
			draw_contour_points(draw_seg, conf_trans(contour1_seg[seg_index]), shift_seg, 5, 2);
			draw_contour_points(draw_seg, conf_trans(contour2_seg[seg_index]), shift_seg+(contour1_seg[seg_index][0].point- contour2_seg[seg_index][0].point), 7, 2);
			
			//cout << morph_seg.size() << "   "<< morph_seg[0].point<< morph_seg[1].point<< morph_seg[2].point<<endl<<final_con.size() << final_con[0].point<< final_con[1].point<<endl;
		}
		cout << "contour1:  " << cnum1 << "    contour2:  " << cnum2 <<"   final: "<< final_con.size()<< endl;
		//FOR(j, 0, final_con.size()) cout << final_con[j].point << endl;
		string save = SavePath;
		imshow("frame.png", draw_seg);
		

		Mat draw_result = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		Point2f shift1 = Point2f(0.5*draw_result.cols, 0.5*draw_result.rows) - cen1;
		int dd = 0;
		FOR(tt, 0, contour1_seg.size())
		{
			FOR(i, 0, contour1_seg[tt].size()) 				
				circle(draw_result, contour1_seg[tt][i].point + shift1, 2, Scalar(255, 0, 0), -1);
			FOR(i, 0, contour2_seg[tt].size())
				circle(draw_result, contour2_seg[tt][i].point + shift1, 2, Scalar(0, 255, 0), -1);
		}
		FOR(i, 0, final_con.size())
			circle(draw_result, final_con[i].point + shift1, 2, Scalar(0, 0, 255), -1);
		imshow("morphing result show", draw_result);

		return final_con;

	}

	vector<Point_f> Tiling_opt::morph_segment(vector<Point_f> seg1, vector<Point_f> seg2, Point_f start, double ratio, double &num_error) //start 是上一个的尾端，是固定的
	{
		Point_f end = Point_f(ratio*seg1.back().point + (1- ratio)*seg2.back().point, seg1.back().type);//seg1.back();  //如果end.type==3, end 不需要变
		double length_ave = length_2p(start.point, end.point);
		double angle_ave;
		cout << "1----start: " << start.point << "     --end: " << end.point <<"  type:"<< end.type<< endl;
		//确定框架的参数 
		if (end.type < 4) // 确定新的end点 
		{
			double angle1 = cos_2v(Point2f(1, 0), seg1.back().point - seg1[0].point);
			angle1 = acos(angle1);
			if (sin_2v(Point2f(1, 0), seg1.back().point - seg1[0].point) < 0) angle1 = -angle1;
			double angle_12 = cos_2v(Point2f(1, 0), seg2.back().point - seg2[0].point);
			angle_12 = acos(angle_12);
			if (sin_2v(Point2f(1, 0), seg2.back().point - seg2[0].point) < 0) angle_12 = -angle_12;
			angle_ave = ratio*angle1 + (1 - ratio)*angle_12;
			length_ave = ratio*length_2p(seg1.back().point, seg1[0].point) + (1 - ratio)* length_2p(seg2.back().point, seg2[0].point);
			//cout << "angle1: " << angle1 << " angle_12: " << angle_12 << " length: " << length_2p(seg1.back().point, seg1[0].point) << "   " << length_2p(seg2.back().point, seg2[0].point) << "  " << length_ave << endl;
			Point2f vec_fin = Point2f(cos(angle_ave), sin(angle_ave));
			end.point = start.point + length_ave*vec_fin;
			//cout << start.point << "   " << length_ave << "    " << vec_fin << endl;

			//double angle1 = cos_2v(Point2f(1, 0), seg1.back().point - seg1[0].point);
			//angle1 = (angle1 > 0.999) ? 0.999 : (angle1 < -0.999) ? -0.999 : angle1; //防止出现acos(1)的情况，会返回错误值

		}
		cout << "2----start: " << start.point << "     --end: " << end.point << endl;
		Point2f mid_des = start.point + 0.5*length_ave * vertical_vec(end.point - start.point);
		vector<Point2f> aff_des;
		aff_des.push_back(start.point);
		aff_des.push_back(end.point);
		aff_des.push_back(mid_des);

		int number_new = 0.5*(seg1.size() + seg2.size());
		double diff = 0.5*(seg1.size() + seg2.size()) - number_new;
		if ((diff + num_error)==1.0)
		{
			number_new += 1;
			num_error = 0;
		}
		else if (diff > num_error) num_error = diff;
		//按照number_new进行重采样,这里我们在原样的基础上进行变形
		vector<Point2f> seg_sam1 = sampling_seg(conf_trans(seg1), number_new);
		vector<Point2f> seg_sam2 = sampling_seg(conf_trans(seg2), number_new);
		if (seg_sam1.size() != seg_sam2.size() || seg_sam1.size() != number_new)
		{
			cout << "seg_sam1.size() != seg_sam2.size() != number_new" << endl;
		}
		//construct Affine matrix
		vector<Point2f> aff_seg1;
		aff_seg1.push_back(seg_sam1[0]);
		aff_seg1.push_back(seg_sam1.back());
		aff_seg1.push_back(seg_sam1[0] + 0.5*length_2p(seg_sam1.back(), seg_sam1[0]) * vertical_vec(seg_sam1.back() - seg_sam1[0]));
		Mat M_seg1 = getAffineTransform(aff_seg1, aff_des);
		//cout << "M_seg1: " << M_seg1 << endl;
		cv::transform(seg_sam1, seg_sam1, M_seg1);
		vector<Point2f> aff_seg2;
		aff_seg2.push_back(seg_sam2[0]);
		aff_seg2.push_back(seg_sam2.back());
		aff_seg2.push_back(seg_sam2[0] + 0.5*length_2p(seg_sam2.back(), seg_sam2[0]) * vertical_vec(seg_sam2.back() - seg_sam2[0]));
		Mat M_seg2 = getAffineTransform(aff_seg2, aff_des);
		cv::transform(seg_sam2, seg_sam2, M_seg2);

		vector<Point2f> contour_morphed;
		for (int j = 0; j < number_new; j++)
		{
			contour_morphed.push_back(ratio*seg_sam1[j] + (1 - ratio)*seg_sam2[j]);
		}
		//reassign type of point_f
		vector<Point2f> contour_final= contour_morphed;
		vector<Point_f> contour_final_;
		FOR(m, 0, contour_final.size())  contour_final_.push_back(p2fea(contour_final[m], 0));
		contour_final_[0].type = seg1[0].type;
		contour_final_.back().type = seg1.back().type;

		//cout << "Points num of each seg: " << seg1.size() << "  " << seg2.size() << "   " << contour_final_.size() << endl;
		return contour_final_;

	}


	vector<Point_f> Tiling_opt::morphing_dir(vector<Point_f> c_mid, vector<Point_f> c_cand, vector<pair<int, int>> path, double ratio)
	{
		Mat draw = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
		draw_contour_points(draw, conf_trans(c_mid), OP, 9);
		draw_contour_points(draw, conf_trans(c_cand), OP, 8);
		cout << "mid_c.size: " << c_mid.size() << "  cand_c.size: " << c_cand.size() << endl;
		int psize = path.size();
		//vector<int> cand_ap = c_mid.anchor_points;// c_cand.anchor_points;
		vector<Point2f> final_mid;
		vector<int> pair_indexs;
		int index = 0; //按照c_mid的顺序
		FOR(i, 0, psize)
		{
			if (path[i].first == index)
			{
				pair_indexs.push_back(path[i].second);
			}
			else
			{
				Point2f new_p;
				/*if (index2 == cand_ap[0] || index2 == cand_ap[1] || index2 == cand_ap[2] || index2 == cand_ap[3])
				{
					new_p = c_cand.contour[index2];
				}
				else
				{
					FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid.contour[pair_indexs[t]] + (1 - ratio) * c_cand.contour[index2]);
					new_p = 1.0 / pair_indexs.size()*new_p;
				}*/
				FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid[index].point + (1 - ratio) * c_cand[pair_indexs[t]].point);
				new_p = 1.0 / pair_indexs.size()*new_p;
				final_mid.push_back(new_p);
				index = path[i].first;
				pair_indexs.swap(vector<int>());
				pair_indexs.push_back(path[i].second);
			}
		}
		if (!pair_indexs.empty())
		{
			Point2f new_p;
			/*if (index2 == cand_ap[0] || index2 == cand_ap[1] || index2 == cand_ap[2] || index2 == cand_ap[3])
			{
				new_p = c_cand.contour[index2];
			}
			else
			{
				FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid.contour[pair_indexs[t]] + (1 - ratio) * c_cand.contour[index2]);
				new_p = 1.0 / pair_indexs.size()*new_p;
			}*/
			FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid[index].point + (1 - ratio) * c_cand[pair_indexs[t]].point);
			new_p = 1.0 / pair_indexs.size()*new_p;
			final_mid.push_back(new_p);
		}

		draw_contour_points(draw, final_mid, OP,4);
		imshow("final", draw);
		prototile_fin = protoTile(final_mid);
		prototile_fin.contour_f = set_flags(final_mid, c_mid);

		cout << "final_mid size: " << final_mid.size() << endl;
		return prototile_fin.contour_f;
	}
}