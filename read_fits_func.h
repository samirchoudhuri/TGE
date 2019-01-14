
void printerror( int status);
int read_fits_header(char* in_file);
int read_ranpar(long grp);
int read_data(long group);
//int write_fits(double *gval,char* out_file);
//int write_fitscorrlconv(char* out_file);
int close_fits();
//float winf(float u);
//void madfilter (int n, float[], float []);
//void Filterwindow(int size, int xpos, int ypos, double *Winfre, double *Winfim, double *Wufil);
int readfits(char* ,double ,double ,int ,int );
