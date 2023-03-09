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
	
	int f = 0;

	if (f == 0)
	{
		bool control_parameter = true;
		if (control_parameter)  tiling.load_para("para.txt");
		else image_id = "46";
		tiling.tiliing_generation(image_id);

	}
	if (f == 1) //specify
	{
		tiling.tiliing_gen_specify("46");

	}
	else if (f == 10) //save contours
	{
		tiling.load_dataset(false);
	}

	waitKey(0);
}