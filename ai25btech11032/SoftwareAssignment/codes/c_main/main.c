#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#include "../c_libs/stb_image_write.h"

double *read_image(const char *filename, int *w, int *h) {
    int ch;
    unsigned char *data = stbi_load(filename, w, h, &ch, 1);
    if (!data) {
        fprintf(stderr, "Failed to load %s\n", filename);
        return NULL;
    }

    size_t N = (size_t)(*w) * (size_t)(*h);
    double *img = malloc(N * sizeof(double));
    if (!img) { stbi_image_free(data); return NULL; }

    for (size_t i = 0; i < N; i++)
        img[i] = data[i] / 255.0;

    stbi_image_free(data);
    return img;
}
void write_image_jpg(const char *filename, const double *img, int w, int h) {
    size_t N = (size_t)w * (size_t)h;
    unsigned char *buf = malloc(N);
    for (size_t i = 0; i < N; i++) {
        double x = img[i];
        if (x < 0) x = 0;
        if (x > 1) x = 1;
        buf[i] = (unsigned char)(x * 255.0 + 0.5);
    }

    stbi_write_jpg(filename, w, h, 1, buf, 90);

    free(buf);
}


double dot(int n,double *x,double *y){
	double result=0.0;
	for(int i=0;i<n;i++){
		result=result+x[i]*y[i];
	}
	return result;
}

double norm(int n,double *x){
	double result=dot(n,x,x);
	return sqrt(result);
}

void scale(double a, int n, double *x) {
    for (int i = 0; i < n; i++)
        x[i] = a * x[i];
}

void axpy(double a, int n, double *x, double *y) {
    for (int i=0;i<n;i++){
	    y[i] += a * x[i];
    }
}

void matvec(int m, int n, double *A, double *x, double *y){
	for(int i=0;i<m;i++){
		y[i]=0;
		for(int j=0;j<n;j++){
			y[i]+=A[i*n+j]*x[j];
		}
	}
}

void matvecT(int m, int n, double *A, double *x, double *y){
	for(int i=0;i<n;i++){
		y[i]=0;
		for(int j=0;j<m;j++){
			y[i]+=A[j*n+i]*x[j];
		}
	}
}

void copy(int n, double *x, double *y){
	for(int i=0;i<n;i++){
		y[i]=x[i];
	}
}

double *zeros(int n){
        double *z = malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
		z[i] = 0.0;
	}
	return z;
}

void randominitialunitvector(int n, double *x){
	for(int i=0;i<n;i++){
		double r=rand()/(double)RAND_MAX;
		x[i]=2.0*r-1;
	}
	double m=norm(n,x);
	if(m<1e-12){
		x[0]=1.0;
		for(int i=1;i<n;i++){
			x[i]=0.0;
		}
		m=1.0;
	}
	scale(1.0/m,n,x);
}

void projection_out(double *q, int k, int n, double *x){
	double c=0.0;
	for(int i=0;i<k;i++){
		double *p=q+i*n;
		c=dot(n,p,x);
		axpy(-c,n,p,x);
	}
}

double *V_basis=NULL;
int V_count=0;

double power_iteration(double *A, int m, int n, double *u, double *v, int max, double err){
	randominitialunitvector(n,v);
	double *prev=zeros(n);

	double sigma=0.0;

	for(int i=0;i<max;i++){
		copy(n,v,prev);

		//if (V_basis && V_count > 0) {
                   //projection_out(V_basis, V_count, n, v);
                   //double nv = norm(n, v);
                   //if (nv < 1e-12) {
                      //v[0] = 1.0;
                      //for(int t = 1; t < n; t++) v[t] = 0.0;
                      //nv = 1.0;
                   //}
                   //scale(1.0/nv, n, v);
                //}

		matvec(m,n,A,v,u);
		double normofu=norm(m,u);
		if(normofu<1e-12){
			u[0]=1.0;
			for(int i=1;i<m;i++){
				u[i]=0.0;
			}
			normofu=1.0;

		}
		scale(1.0/normofu,m,u);

		matvecT(m,n,A,u,v);
		double normofv=norm(n,v);
		if(normofv<1e-12){
			v[0]=1.0;
			for(int i=1;i<n;i++){
				v[i]=0.0;
			}
			normofv=1.0;
		}
		scale(1.0/normofv,n,v);

		sigma=normofu;

		double diff=fabs(dot(n,v,prev))-1.0;
		if(fabs(diff)<err){
			break;
		}
	}
	free(prev);
        return sigma;
}

void add_rank1(double *A, int m, int n, double sigma, double *u, double *v){
	if(sigma==0.0) return;
	for(int i=0;i<m;i++){
		double s=sigma*u[i];
		double *Ai=A+(size_t)i*n;
		if(s==0.0){
			continue;
		}
		for(int j=0;j<n;j++){
			Ai[j]+=s*v[j];
		}
	}
}
int top_k_svd(double *A, int m, int n, int k, double *Sigma, double *U, double *V, int max, double err){
	double *u=zeros(m);
	double *v=zeros(n);
	if (!u || !v) { free(u); free(v); return -1; }
	V_basis=V;

	for(int i=0;i<k;i++){
		V_count=i;

		double sigma=power_iteration(A,m,n,u,v,max,err);

		double nu=norm(m,u);
		if(nu<1e-12){
			u[0]=1.0;
			for(int j=1;j<m;j++){
				u[j]=0.0;
			}
			nu=1.0;
		}
		scale(1.0/nu,m,u);
		copy(m,u,U+i*m);

		double nv=norm(n,v);
		if(nv<1e-12){
			v[0]=1.0;
			for(int j=1;j<n;j++){
				v[j]=0.0;
			}
			nv=1.0;
		}
		scale(1.0/nv,n,v);
		copy(n,v,V+i*n);

		Sigma[i]=sigma;

		add_rank1(A, m, n, -sigma, u, v);

	}
	V_basis=NULL;
	V_count=0;
	free(u);
	free(v);
	return 0;
}


void reconstruct_rank_k(double *Ak, int m, int n, int k, double *Sigma, double *U, double *V){
	for(int i=0;i<m*n;i++){
		Ak[i]=0.0;
	}
	for(int i=0;i<k;i++){
		double *u=U+i*m;
		double *v=V+i*n;
		double sigma=Sigma[i];
		add_rank1(Ak,m,n,sigma,u,v);
	}
}

double frobenius_norm(int m, int n, double *A){
	double e=0.0;
	for(int i=0;i<m*n;i++){
		double d=A[i];
		e+=d*d;
	}
	return sqrt(e);
}

double frobenius_error(int m, int n, double *A, double *Ak){
	double e=0.0;
	for(int i=0;i<m*n;i++){
		double d=A[i]-Ak[i];
		e+=d*d;
	}
	return sqrt(e);
}

int main() {
    int w, h;

    char infile[256],outfile[256];
    printf("Enter input file name: ");
    scanf("%s",infile);
    printf("Enter output file name: ");
    scanf("%s",outfile);

    char input_path[512], output_path[512];
    snprintf(input_path, sizeof(input_path), "../../figs/%s", infile);
    snprintf(output_path, sizeof(output_path), "../../figs/%s", outfile);

    double *A = read_image(input_path, &w, &h);
    if (!A) return 1;

    double *A0 = malloc((size_t)h*w*sizeof(double));
    for (size_t t = 0; t < (size_t)h*w; ++t) A0[t] = A[t];

    

    int m = h;
    int n = w;
    int k;
    printf("Enter rank k: ");
    scanf("%d",&k);

    double *Sigma = malloc(k * sizeof(double));
    double *U = malloc((size_t)m * k * sizeof(double));
    double *V = malloc((size_t)n * k * sizeof(double));

    size_t N = (size_t)m * n;
    //double mu = 0.0;
    //for (size_t t = 0; t < N; ++t) mu += A[t];
    //mu /= (double)N;
    //for (size_t t = 0; t < N; ++t) A[t] -= mu;
    

    top_k_svd(A, m, n, k, Sigma, U, V, 100, 1e-10);

    double *Ak = calloc((size_t)m * n, sizeof(double));
    reconstruct_rank_k(Ak, m, n, k, Sigma, U, V);
    //for (size_t t = 0; t < (size_t)m * n; ++t) Ak[t] += mu;


    double mn = 1e300, mx = -1e300, mean = 0.0;
    for (size_t t = 0; t < N; ++t) {
        double z = Ak[t];
        if (z < mn) mn = z;
        if (z > mx) mx = z;
        mean += z;
    }
    mean /= (double)N;
    printf("Ak stats: min=%.6f max=%.6f mean=%.6f\n", mn, mx, mean);

    double err = frobenius_error(m, n, A0, Ak);
    printf("Frobenius error  ||A - Ak||_F = %.6f  (pixel-scale = %.3f)\n",err, err * 255.0);
    double abs = frobenius_norm(m,n,A0);
    printf("Relative Frobenius error  ||A - Ak||_F / ||A||_F = %.6f\n",err/abs);

    double ratio = ((double)m*n)/((double)k*(m+n+1));
    printf("Compression ratio = %.6f\n", ratio);
    

    write_image_jpg(output_path, Ak, n, m);

    free(A);
    free(Ak);
    free(Sigma);
    free(U);
    free(V);
    return 0;
}






