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
		bool control_parameter = false;
		if (control_parameter)  tiling.load_para("para.txt");
		else image_id = "44";
		tiling.tiliing_generation(image_id);

	}
	if (f == 1) //specify
	{
		vector<int> anc_points = {6,19,24,31};
		tiling.tiliing_gen_specify("44", anc_points);


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
		tiling.load_dataset(false);
	}

	waitKey(0);
}