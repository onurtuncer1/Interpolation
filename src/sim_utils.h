#if !defined(SIM_UTILS_H)

// TODO(batuhan): std::move ? memmove ? not available ?
#define SIM_SWAP(X, Y)			do { auto temp = std::move(X); X = std::move(Y); Y = std::move(temp); } while(0,0);
#define SIM_MIN(X, Y)			(((X) < (Y)) ? (X) : (Y))
#define SIM_MAX(X, Y)			(((X) > (Y)) ? (X) : (Y))

#if defined(_DEBUG)
#define SIM_ASSERT(expr)		if (!(expr)) { *(int*)0 = 0; }
#else
#define SIM_ASSERT(expr)
#endif

#define SIM_ARRAY_COUNT(arr)    sizeof(arr) / sizeof(arr[0])

#define SIM_UTILS_H
#endif