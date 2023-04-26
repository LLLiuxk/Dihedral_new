#include "tilingOpt.h"
#include  <stdio.h>
#include  <stdlib.h>

string image_id;

using namespace Tiling_tiles;



int main()
{
	clock_t start, midtime, finish;
	start = clock();
	Tiling_opt tiling;
	
	int f = 1;

	if (f == 0)
	{
		bool control_parameter = true;
		if (control_parameter)  tiling.load_para("para.txt");
		else image_id = "192";
		tiling.tiliing_generation(image_id);

	}
	if (f == 1) //specify
	{
		// 5: 6,11,16,29   13: 1,12,27,36  14: 0,4,17,24   22: 8,19,28,36   28: 6,23,25,35  31: 0,15,18,27   37: 10,18,25,37
		//44:6,19,24,31    71:4,14,22,32    70: 0,8,19,28    192: 0,9,18,27   157:6,12,23,31
		vector<int> anc_points = { 4,14,22,32 };
		tiling.tiliing_gen_specify2("71", anc_points);


	}
	if (f == 2)
	{
		//tiling.load_dataset(true);
		string filepath = DefaultPath;
		string filepath1 = filepath + "contour/189.txt";
		string filepath2 = filepath + "contour/191.txt";
		tiling.prototile_first = protoTile(filepath1);
		tiling.prototile_second = protoTile(filepath2);
		double sc_inner = 0;
		vector<vector<double>> first_arr = compute_TAR(tiling.prototile_first.contour, sc_inner);
		vector<vector<double>> second_arr= compute_TAR(tiling.prototile_second.contour, sc_inner);
		vector<pair<int, int>> path;
		int sec_shift = 0;
		int width = WindowsWidth;
		tiling.feature_match(first_arr, second_arr, tiling.prototile_first.feature_points, tiling.prototile_second.feature_points, path, sec_shift, width);
		cout << "sec_shift: " << sec_shift << endl;

		Mat draw = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
		draw_contour_points(draw, tiling.prototile_first.contour, OP, 5, 2);
		draw_contour_points(draw, tiling.prototile_second.contour, Point2f(400,0), 7, 2);
		FOR(i, 0, tiling.prototile_first.feature_points.size())
			circle(draw, tiling.prototile_first.contour[tiling.prototile_first.feature_points[i]], 3, Scalar(0,0,255), -1);
		FOR(i, 0, tiling.prototile_second.feature_points.size())
			circle(draw, tiling.prototile_second.contour[tiling.prototile_second.feature_points[i]]+Point2f(400, 0), 3, Scalar(0, 125, 255), -1);
		circle(draw, tiling.prototile_first.contour[tiling.prototile_first.feature_points[0]], 5, Scalar(120, 0, 255), -1);
		circle(draw, tiling.prototile_second.contour[tiling.prototile_second.feature_points[0]] + Point2f(400, 0), 5, Scalar(0, 125, 255), -1);
		cout << "path size: " << path.size() << endl;
		FOR(i, 0, path.size())
		{
			int tsfsize=tiling.prototile_second.feature_points.size();
			Point2f f1 = tiling.prototile_first.contour[tiling.prototile_first.feature_points[path[i].first]];
			Point2f f2 = tiling.prototile_second.contour[tiling.prototile_second.feature_points[(path[i].second + sec_shift)% tsfsize]]+ Point2f(400, 0);
			line(draw, f1,f2, colorbar[6].second);
		}
		imshow("fea match", draw);
	}
	else if (f == 10) //save contours
	{
		//tiling.load_dataset(false);
		string filepath = DefaultPath;
		//string filepath1 = filepath + "contour/189.txt";
		string filepath1 = "D:/vs2015project/Dihedral_new/Dihedral_new/mid_shape.txt";
		tiling.prototile_first = protoTile(filepath1);
		MatrixXd V;
		MatrixXi F;
		vector<Point2f> con = tiling.prototile_first.contour;//{ Point2f(100,100),Point2f(200,100),Point2f(150,200) };// 
		//int index_ = add_points(con, 0.1);

		//FOR(m, 0, 5) con.push_back(Point2f(300 + 20 * m, 300 + 15 * m));
		triangulateContour(con,V,F);
		cout << "new V size:  "<< V.size()<<"  new F size: " << F.size() << endl;
		Mat image = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));

		for (size_t i = 0; i < F.rows(); i++)
		{
			Point pt1(V.row(F(i, 0)).x(), V.row(F(i, 0)).y());
			Point pt2(V.row(F(i, 1)).x(), V.row(F(i, 1)).y());
			Point pt3(V.row(F(i, 2)).x(), V.row(F(i, 2)).y());
			cout << F(i, 0) << "   " << F(i, 1) << "   " << F(i, 2) << endl;
			cout << pt1 << "    " << pt2 <<"   "<<pt3<< endl;
			line(image, pt1, pt2, Scalar(0, 255, 0), 1, LINE_AA);
			line(image, pt2, pt3, Scalar(0, 255, 0), 1, LINE_AA);
			line(image, pt3, pt1, Scalar(0, 255, 0), 1, LINE_AA);
		}
		//imwrite("Triangulation.png", image);
		imshow("Triangulation", image);

		write_obj("duck.obj", V,F);
		string obj_path = ParaPath;
		string para_path = ParaPath;
		string command = "D:/vs2015project/ARAP_Deform/x64/Debug/ARAP_Deform.exe "+ obj_path+"duck.obj "+ para_path+"deform_para.txt";
		cout << command << endl;
		system(command.c_str());

	}
	else if (f == 11) //save contours
	{
	}
	finish = clock();
	cout << endl << "All time consumption: " << (double)(finish - start) / CLOCKS_PER_SEC << " s " << endl;
	waitKey(0);
}