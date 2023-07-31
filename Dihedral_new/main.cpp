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
	
	int f = 100;
	// 100: test   0: generate all  1: generate specify  2: show pure contours  3: obj  4: test matching  5: change color  6: dilate
	if (f == 100)  //test item
	{
		// 绘制折线图
		//std::vector<double> c = { 1.0, 2.0, 3.0, 4.0, 5.0 };
		//std::vector<double> v = { 1.0, 4.0, 9.0, 16.0, 25.0 };
		//string path = "";
		//drawLineGraph(c, v, path);

		//vector<vector<Point2f>> tex = load_texture("texture1.txt");
		//vector<Mat> images;
		//for (int g = 0; g < tex.size(); g++)
		//{
		//	Mat image = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		//	draw_contour(image, tex[g], Point2f(0, 0));
		//	images.push_back(image);
		//}
		//write_avi( images, "output.avi",1);


		Mat a = Mat(600, 800, CV_8UC3, Scalar(255, 255, 255));
		vector<Point2f> con = { Point2f(100,100),Point2f(300,100),Point2f(300,300) ,Point2f(100,300) };
		draw_contour(a, con, OP);
		// 创建图像 c
		cv::Mat c;

		// 计算交集并复制到 c
		int new_width = a.cols + 200;
		int new_height = a.rows + 200;
		cv::Mat b(new_height, new_width, a.type(), cv::Scalar(255, 255, 255));
		a.copyTo(b(cv::Rect(200, 0, a.cols, a.rows)));

		for (int y = 0; y < a.rows; ++y) {
			for (int x = 0; x < a.cols; ++x) {
				cv::Vec3b pixel_a = a.at<cv::Vec3b>(y, x);
				cv::Vec3b pixel_b = b.at<cv::Vec3b>(y, x);

				// 判断是否都是黑色
				if (pixel_a[0] == 0 || pixel_a[1] == 0 || pixel_a[2] == 0 ||
					pixel_b[0] == 0 || pixel_b[1] == 0 || pixel_b[2] == 0) {
					// 设置为黑色
					b.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 0);
				}
				else {
					// 设置为白色
					b.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 255, 255);
				}
			}
		}

		// 计算 a 和 b 的交集
		//cv::bitwise_and(a, b, c);

		// 显示原图 a 和结果图像 c
		cv::imshow("Image a", a);
		cv::imshow("Intersection Image c", b);

	}
	if (f == 0)  //generate all
	{
		//// generate from para.txt
		//bool control_parameter = true;
		//if (control_parameter)  tiling.load_para("para.txt");
		//else image_id = "29";
		//tiling.tiliing_generation(image_id);

		for (int g =101; g <= 173; g++)
		{
			image_id = to_string(g);
			cout << image_id << endl;
			Tiling_opt tiling_;
			tiling_.tiliing_generation(image_id);
		}
	}
	if (f == 1) //specify
	{
		// 5: 6,11,16,29     14: 0,4,17,24      28: 6,23,25,35  31: 0,15,18,27   37: 10,18,25,37
		//       70: 0,8,19,28    192: 0,9,18,27   157:6,12,23,31
		//Escher 22: 8,19,28,36    71: 4,14,22,32    172: 6,13,23,32    bad:171: 7,19,29,33
		//New  5: 6,11,16,29    13: 1,12,27,36    37: 10,18,25,37   44:6,19,24,31  44:6,19,25,39  53:0,12,17,36
		vector<int> anc_points = { 4,14,22,32 };
		tiling.tiliing_gen_specify("71", anc_points);

		/*Mat draw2 = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		
		vector<Point2f> old_edge = {Point2f(400,120),Point2f(350,150),Point2f(300,180),Point2f(250,200),Point2f(200,150),Point2f(150,120),Point2f(100,100)};
		vector<Point2f> new_edge = { Point2f(400,170),Point2f(350,200),Point2f(300,230),Point2f(250,300),Point2f(200,200),Point2f(150,170), Point2f(100,150)};
		draw_contour_points(draw2, old_edge, OP, 2, 3);
		draw_contour_points(draw2, new_edge, OP, 5,3);
		bound_recover(old_edge, new_edge);
		draw_contour_points(draw2, new_edge, OP, 7,3);
		imshow("show:", draw2);*/
	}
	if (f == 2) //show pure contours
	{
		string savepath = SaveSpecPath;
		savepath = savepath + "71/";
		vector<Point2f> con1 = load_point_file(savepath + "c1.txt");
		vector<Point2f> con2 = load_point_file(savepath + "c2.txt");
		Mat two_c = Mat(1000, 1600, CV_8UC3, Scalar(255, 255, 255));
		Mat two_c2 = Mat(1000, 1600, CV_8UC3, Scalar(255, 255, 255));
		draw_contour(two_c, con1, OP ,0,3);
		draw_contour(two_c2, con2, OP,0,3);
		imshow("123", two_c);
		imshow("1234", two_c2);

	}
	if (f == 3)  //generate obj
	{
		string savepath = SaveSpecPath;
		savepath = savepath + "71/";
		vector<Point2f> con1 = load_point_file(savepath + "c1.txt");
		vector<Point2f> con2 = load_point_file(savepath + "c2.txt");
		contour2obj(savepath + "c1.obj", con1, OP, 15);
		contour2obj(savepath + "c2.obj", con2, OP, 15);
	}
	else if (f == 4) //test feature matching
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
	else if (f == 5)  //change color
	{
		Mat src;
		string name = "D:\\input.png";
		src = imread(name, IMREAD_COLOR);
		if (src.empty())
		{
			return -1;
		}
		Mat dst = src;
		cvtColor(src, src, COLOR_BGR2GRAY);

		threshold(src, src, 155, 255, cv::THRESH_BINARY); //255 white
														  //imshow("zala", src);
														  //imshow("nine", dst);
		Vec3b color1, color2;
		//Scalar t = colorbar[6].second;
		//color1 = Vec3b(t.val[0], t.val[1], t.val[2]);
		//t = colorbar[7].second;
		//color2 = Vec3b(t.val[0], t.val[1], t.val[2]);

		color1 = Vec3b(255, 255, 255); //Vec3b(130, 149, 225);
		color2 = Vec3b(0, 0, 0);

		for (int i = 0; i < src.rows; i++)
			for (int j = 0; j < src.cols; j++)
			{
				//cout << (int)src.at<uchar>(i, j) << " ";
				if ((int)src.at<uchar>(i, j) == 0) dst.at<Vec3b>(i, j) = color1;
				if ((int)src.at<uchar>(i, j) == 255) dst.at<Vec3b>(i, j) = color2;
			}
		imwrite("D:\\dst.png", dst);
	}
	else if (f == 6)  //加粗，扩张白色
	{
		Mat src, dilation_dst;
		int dilation_elem = 2;
		int dilation_size = 2;
		//int const max_elem = 2;
		//int const max_kernel_size = 21;

		string name = "D:\\dst.png";//"D:\\result.png";
									//string name = "D:\\print.png";
		src = imread(name, IMREAD_COLOR);
		if (src.empty())
		{
			return -1;
		}
		int dilation_type = 0;
		namedWindow("Erosion Demo", WINDOW_AUTOSIZE);
		if (dilation_elem == 0) { dilation_type = MORPH_RECT; }
		else if (dilation_elem == 1) { dilation_type = MORPH_CROSS; }
		else if (dilation_elem == 2) { dilation_type = MORPH_ELLIPSE; }
		Mat element = getStructuringElement(dilation_type,
			Size(2 * dilation_size + 1, 2 * dilation_size + 1),
			cv::Point(dilation_size, dilation_size));
		dilate(src, dilation_dst, element);
		imwrite("D:\\Erosion Demo2.png", dilation_dst);
	}

	finish = clock();
	cout << endl << "All time consumption: " << (double)(finish - start) / CLOCKS_PER_SEC << " s " << endl;
	waitKey(0);
}