#ifdef FORTRAN_SAME
#define FC_FUNC(name,NAME) name
#elif FORTRANUNDERSCORE
#define FC_FUNC(name,NAME) name ##_
#elif FORTRANDOUBLEUNDERSCORE
#define FC_FUNC(name,NAME)  name ##__
#endif
