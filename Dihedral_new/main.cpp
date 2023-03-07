#include "tilingOpt.h"
#include  <stdio.h>
#include  <stdlib.h>



using namespace Tiling_tiles;

int main()
{
	clock_t start, midtime, finish;
	start = clock();
	Tiling_opt tiling;
	
	int f = 0;

	if (f == 0)
	{
		tiling.tiliing_generation("70");

	}
	else if (f == 10) //save contours
	{
		tiling.load_dataset(false);
	}

	waitKey(0);
}