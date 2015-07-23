#include "ctvm.h"

BoostDoubleVector Gradient2D(BoostDoubleVector U, unsigned long pixel) // pixel = actual pixel number -1 ??
{
	unsigned long n = U.size();
	unsigned long l = sqrt(n);
	unsigned long rows = 0;
	unsigned long cols = 0;
	float q = 0;

	BoostDoubleVector Di (2);
	BoostDoubleMatrix X (l, l);

	X = VectorToMatrix(U, l, l);

	for (unsigned long i = 0; i < l; ++i)
	{
		if (!pixel)
		{
			rows = i;
			cols = 0;
			break;
		}
		else
		{
			for (unsigned long j = 0; j < l; ++j)
			{
				q = j / pixel;
				if (q)
				{
					rows = i;
					cols = j;
					break;
				}
			}
			pixel = pixel - l;
		}


		return Di;
	}
}