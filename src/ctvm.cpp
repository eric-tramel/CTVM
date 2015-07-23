#include "ctvm.h"

BoostDoubleVector 2DGradient(BoostDoubleVector U, unsigned long rank) // rank = actual rank number -1! (first rank = U(0))
{
	unsigned long n = U.size();
	unsigned long pixel = rank;
	unsigned long l = sqrt(n); // !
	unsigned long rows = 0;
	unsigned long cols = 0;
	float q = 0;

	BoostDoubleVector Di (2);
	BoostDoubleMatrix X (l, l);

	X = VectorToMatrix(U, l, l);

	for (unsigned long i = 0; i < l; ++i) /* Find pixel place in the matrix */ 
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
	}
	if (cols == l-1 || rows == l-1)
	{
		Di(0) = X(rows, cols) - X(rows, cols);
		Di(1) = X(rows, cols) - X(rows, cols);
	}
	else if (cols == l-1)
	{
		Di(0) = X(rows, cols) - X(rows, cols);
		Di(1) = X(rows + 1, cols) - X(rows, cols);
	}
	else if (rows == l-1)
	{
		Di(0) = X(rows, cols + 1) - X(rows, cols);
		Di(1) = X(rows, cols) - X(rows, cols);
	}
	else
	{
		Di(0) = X(rows, cols + 1) - X(rows, cols);
		Di(1) = X(rows + 1, cols) - X(rows, cols);
	}

return Di;
}