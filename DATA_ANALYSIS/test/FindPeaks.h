// FindPeaks.h

void FindPeaks(std::vector <std::vector<double> > &trace, double thr_to_find_peaks);
void DLED(std::vector <std::vector<double> > &trace, int dleddt);
int find_peak_fix_time(int mintp, int maxtp);
void FindPeaksRisingFalling(std::vector <std::vector<double> > &t, double thr, int max_peak_width, int rising_points, int falling_points);
