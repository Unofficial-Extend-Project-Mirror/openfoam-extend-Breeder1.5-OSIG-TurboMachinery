#include <sys/types.h>
#if !defined(_WIN32) || defined(__NUTC__)
#include <sys/times.h>
#endif
#include <sys/stat.h>
#include <time.h>

double elapsed_time (void)
{
#if defined(_WIN32) && !defined(__NUTC__)
    return (double)clock() / (double)CLK_TCK;
#else
    struct tms clock;
    double usrtime;
    double systime;

    times(&clock);
    usrtime = (double) clock.tms_utime / (double)CLK_TCK;
    systime = (double) clock.tms_stime / (double)CLK_TCK;
    return (usrtime + systime);
#endif
}

double file_size (char *fname)
{
    struct stat st;

    if (stat (fname, &st)) return 0.0;
    return (double)st.st_size / 1048576.0;
}

