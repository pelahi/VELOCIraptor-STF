#include <iostream>
#include <fstream>
#include <Morphology.h>

using namespace std;
using namespace Math;
namespace NBody
{
    static void CalculateTensor(Matrix& M, const Double_t q, const Double_t s, const Int_t nbodies, Particle *p);
    static void CalculateTensorWithMass(Matrix& M, const Double_t q, const Double_t s, const Int_t nbodies, Particle *p);
    static void CalculateTensorWithMassInDenom(Matrix& M, const Double_t q, const Double_t s, const Int_t nbodies, Particle *p);
    static void RotateParticles(const Int_t nbodies, Particle *p, const Matrix& M, const Coordinate& e);

    void GetGlobalMorphology(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, int verbose)
    {
        // Calculate the axial ratios q and s.
        Double_t oldq;
        Double_t olds;

        Coordinate e;
        // Iterative procedure.  See Dubinski and Carlberg (1991).
        int i=0;
        do
        {
            Matrix M = Matrix(0.0);
            CalculateTensor(M, q, s, nbodies, p);
            e = M.Eigenvalues();
            oldq = q;
            olds = s;
            q = sqrt(e[1] / e[0]);
            s = sqrt(e[2] / e[0]);

            if (verbose)
            {
                cout<<"Iteration "<<i<<endl;
                cout << "Interia Tensor is: "<<endl;
                for (int j=0;j<3;j++){
                    for (int k=0;k<3;k++)cout<<M(j,k)<<" ";
                    cout<<endl;
                }
                cout << "Eigenvalues are :"<<e[0]<<" "<<e[1]<<" "<<e[2]<<endl;
                cout << "Axis ratios are :" << q <<" "<<s<<endl;
            }
            RotateParticles(nbodies, p, M, e);
            i++;
        } while (fabs(olds - s) > Error || fabs(oldq - q) > Error);
    }
    void GetGlobalMorphology(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, Matrix& eigenvec, int verbose)
    {
        // Calculate the axial ratios q and s.
        Double_t oldq;
        Double_t olds;

        Coordinate e;
        Matrix R=Matrix(0.);
        Matrix oldR=Matrix(0.);
        // Iterative procedure.  See Dubinski and Carlberg (1991).
        int i=0;
        do
        {
            Matrix M = Matrix(0.0);
            Matrix eigenvecp=Matrix(0.);
            CalculateTensor(M, q, s, nbodies, p);
            e = M.Eigenvalues();

            oldq = q;
            olds = s;
            q = sqrt(e[1] / e[0]);
            s = sqrt(e[2] / e[0]);

            //if i=0 get initial estimate.
            if (i==0){
                eigenvec=M.Eigenvectors(e);
                //store rotation matrix.
                oldR=eigenvec.Transpose();
            }
            //afterwards system has been rotated so must take all rotations into account and 
            //rotate final answer back to initial coordinate system.
            //calculate new matrix 
            else if (i==1)
            {
                eigenvecp=M.Eigenvectors(e);
                R=oldR;
                oldR=eigenvecp.Transpose();
                //assuming vectors stored in row
                for (int j=0;j<3;j++){
                    Coordinate ep(eigenvecp(0,j),eigenvecp(1,j),eigenvecp(2,j));
                    eigenvec(0,j)=R(0,0)*ep[0]+R(0,1)*ep[1]+R(0,2)*ep[2];
                    eigenvec(1,j)=R(1,0)*ep[0]+R(1,1)*ep[1]+R(1,2)*ep[2];
                    eigenvec(2,j)=R(2,0)*ep[0]+R(2,1)*ep[1]+R(2,2)*ep[2];
                }
            }
            else
            {
                eigenvecp=M.Eigenvectors(e);
                R=R*oldR;
                oldR=eigenvecp.Transpose();
                //assuming vectors stored in row
                for (int j=0;j<3;j++){
                    Coordinate ep(eigenvecp(0,j),eigenvecp(1,j),eigenvecp(2,j));
                    eigenvec(0,j)=R(0,0)*ep[0]+R(0,1)*ep[1]+R(0,2)*ep[2];
                    eigenvec(1,j)=R(1,0)*ep[0]+R(1,1)*ep[1]+R(1,2)*ep[2];
                    eigenvec(2,j)=R(2,0)*ep[0]+R(2,1)*ep[1]+R(2,2)*ep[2];
                }
            }

            if (verbose)
            {
                cout<<"Iteration "<<i<<endl;
                cout << "Interia Tensor is: "<<endl;
                for (int j=0;j<3;j++){
                    for (int k=0;k<3;k++)cout<<M(j,k)<<" ";
                    cout<<endl;
                }
                cout << "Eigenvectors are: "<<endl;
                for (int j=0;j<3;j++){
                    for (int k=0;k<3;k++)cout<<eigenvec(j,k)<<" ";
                    cout<<endl;
                }
                cout << "Eigenvalues are :"<<e[0]<<" "<<e[1]<<" "<<e[2]<<endl;
                cout << "Axis ratios are :" << q <<" "<<s<<endl;
            }
            RotateParticles(nbodies, p, M, e);
            i++;
        } while (fabs(olds - s) > Error || fabs(oldq - q) > Error);
    }
    void GetGlobalMorphologyWithMass(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, int mdenom, int verbose)
    {
        // Calculate the axial ratios q and s.
        Double_t oldq;
        Double_t olds;

        Coordinate e;
        // Iterative procedure.  See Dubinski and Carlberg (1991).
        int i=0;
        do
        {
            Matrix M = Matrix(0.0);
            if (mdenom==1) CalculateTensorWithMassInDenom(M, q, s, nbodies, p);
            else CalculateTensorWithMass(M, q, s, nbodies, p);
            e = M.Eigenvalues();
            
            oldq = q;
            olds = s;
            q = sqrt(e[1] / e[0]);
            s = sqrt(e[2] / e[0]);

            if (verbose)
            {
                cout<<"Iteration "<<i<<endl;
				cout << "Interia Tensor is: "<<endl;
				for (int j=0;j<3;j++){
					for (int k=0;k<3;k++)cout<<M(j,k)<<" ";
					cout<<endl;
				}
                cout << "Eigenvalues are :"<<e[0]<<" "<<e[1]<<" "<<e[2]<<endl;
                cout << "Axis ratios are :" << q <<" "<<s<<endl;
            }

            RotateParticles(nbodies, p, M, e);

            i++;
        } while (fabs(olds - s)/olds > Error || fabs(oldq - q)/oldq > Error);
    }
    void GetGlobalMorphologyWithMass(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, Matrix& eigenvec, int mdenom, int verbose)
    {
        // Calculate the axial ratios q and s.
        Double_t oldq;
        Double_t olds;

        Coordinate e;
        Matrix R=Matrix(0.);
        Matrix oldR=Matrix(0.);
        // Iterative procedure.  See Dubinski and Carlberg (1991).
        int i=0;
        do
        {
            Matrix M = Matrix(0.0);
            Matrix eigenvecp=Matrix(0.);
            if (mdenom==1) CalculateTensorWithMassInDenom(M, q, s, nbodies, p);
            else CalculateTensorWithMass(M, q, s, nbodies, p);
            e = M.Eigenvalues();
            oldq = q;
            olds = s;
            q = sqrt(e[1] / e[0]);
            s = sqrt(e[2] / e[0]);

            //if i=0 get initial estimate.
            if (i==0){
                eigenvec=M.Eigenvectors(e);
                //store rotation matrix.
                oldR=eigenvec.Transpose();
            }
            //afterwards system has been rotated so must take all rotations into account and 
            //rotate final answer back to initial coordinate system.
            //calculate new matrix 
            else if (i==1)
            {
                eigenvecp=M.Eigenvectors(e);
                R=oldR;
                oldR=eigenvecp.Transpose();
                //assuming vectors stored in row
                for (int j=0;j<3;j++){
                    Coordinate ep(eigenvecp(0,j),eigenvecp(1,j),eigenvecp(2,j));
                    eigenvec(0,j)=R(0,0)*ep[0]+R(0,1)*ep[1]+R(0,2)*ep[2];
                    eigenvec(1,j)=R(1,0)*ep[0]+R(1,1)*ep[1]+R(1,2)*ep[2];
                    eigenvec(2,j)=R(2,0)*ep[0]+R(2,1)*ep[1]+R(2,2)*ep[2];
                }
            }
            else
            {
                eigenvecp=M.Eigenvectors(e);
                R=R*oldR;
                oldR=eigenvecp.Transpose();
                //assuming vectors stored in row
                for (int j=0;j<3;j++){
                    Coordinate ep(eigenvecp(0,j),eigenvecp(1,j),eigenvecp(2,j));
                    eigenvec(0,j)=R(0,0)*ep[0]+R(0,1)*ep[1]+R(0,2)*ep[2];
                    eigenvec(1,j)=R(1,0)*ep[0]+R(1,1)*ep[1]+R(1,2)*ep[2];
                    eigenvec(2,j)=R(2,0)*ep[0]+R(2,1)*ep[1]+R(2,2)*ep[2];
                }
            }

            if (verbose)
            {
                cout<<"Iteration "<<i<<endl;
                cout << "Interia Tensor is: "<<endl;
				for (int j=0;j<3;j++){
					for (int k=0;k<3;k++)cout<<M(j,k)<<" ";
					cout<<endl;
				}
                cout << "Eigenvectors are: "<<endl;
				for (int j=0;j<3;j++){
					for (int k=0;k<3;k++)cout<<eigenvec(j,k)<<" ";
					cout<<endl;
				}
                cout << "Eigenvalues are :"<<e[0]<<" "<<e[1]<<" "<<e[2]<<endl;
                cout << "Axis ratios are :" << q <<" "<<s<<endl;
            }

            RotateParticles(nbodies, p, M, e);

            i++;
        } while (fabs(olds - s) > Error || fabs(oldq - q) > Error);
    }

	// System reference interface
    void GetGlobalMorphology(System& S, Double_t& q, Double_t& s, Double_t Error, int verbose)
    {
		GetGlobalMorphology(S.GetNumParts(), S.Parts(),q,s,Error,verbose);
    }
    void GetGlobalMorphology(System& S, Double_t& q, Double_t& s, Double_t Error, Matrix& eigenvec, int verbose)
    {
	GetGlobalMorphology(S.GetNumParts(), S.Parts(), q, s, Error, eigenvec, verbose);
    }
    void GetGlobalMorphologyWithMass(System& S, Double_t& q, Double_t& s, Double_t Error, int mdenom, int verbose)
    {
		GetGlobalMorphologyWithMass(S.GetNumParts(), S.Parts(), q, s, Error, mdenom, verbose);
    }

    void GetGlobalMorphologyWithMass(System& S, Double_t& q, Double_t& s, Double_t Error, Matrix& eigenvec, int mdenom, int verbose)
    {
		GetGlobalMorphologyWithMass(S.GetNumParts(), S.Parts(), q, s, Error, eigenvec, mdenom, verbose);
    }

/*
    void GetInertiaTensor(System& S, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec, Matrix&M)
    {
 		Double_t r2,det;
        Coordinate e;
        for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)M(j, k) = 0;
        for (Int_t i = 0; i < S.GetNumParts(); i++)
        {
			r2=S.Part(i).X()*S.Part(i).X()+S.Part(i).Y()*S.Part(i).Y()+S.Part(i).Z()*S.Part(i).Z();
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    M(j, k) += S.Part(i).GetMass()*((j==k)*r2-S.Part(i).GetPosition(j) * 
                               S.Part(i).GetPosition(k));
                }
            }
        }
		det=M.Det();
        for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) M(j, k) /= det;
		e = M.Eigenvalues();
		a=(e[0]);b=(e[1]);c=(e[2]);
        eigenvec=M.Eigenvectors(e);
        RotateParticles(S, M, e);
        for (Int_t i = 0; i < S.GetNumParts(); i++)
        {
			r2=S.Part(i).X()*S.Part(i).X()+S.Part(i).Y()*S.Part(i).Y()+S.Part(i).Z()*S.Part(i).Z();
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    M(j, k) += S.Part(i).GetMass()*((j==k)*r2-S.Part(i).GetPosition(j) * 
                               S.Part(i).GetPosition(k));
                }
            }
        }
		det=M.Det();
        for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) M(j, k) /= det;
		e = M.Eigenvalues();
		a=(e[0]);b=(e[1]);c=(e[2]);
        eigenvec=M.Eigenvectors(e);
    }
*/
    void GetInertiaTensor(const Int_t n, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec)
    {
        Double_t r2,det;
        Coordinate e;
        Matrix M(0.);
        
        for (Int_t i = 0; i < n; i++)
        {
            r2=p[i].X()*p[i].X()+p[i].Y()*p[i].Y()+p[i].Z()*p[i].Z();
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    M(j, k) += p[i].GetMass()*((j==k)*r2-p[i].GetPosition(j) * 
                            p[i].GetPosition(k));
                }
            }
        }
        det=M.Det();
        for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) M(j, k) /= det;
        e = M.Eigenvalues();
        a=e[0];b=e[1];c=e[2];
        eigenvec=M.Eigenvectors(e);
    }

    void GetInertiaTensor(const Int_t n, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec, Matrix &I)
    {
        Double_t r2,det;
        Coordinate e;
        I=Matrix(0.);
        
        for (Int_t i = 0; i < n; i++)
        {
            r2=p[i].X()*p[i].X()+p[i].Y()*p[i].Y()+p[i].Z()*p[i].Z();
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    I(j, k) += p[i].GetMass()*((j==k)*r2-p[i].GetPosition(j) * 
                            p[i].GetPosition(k));
                }
            }
        }
        det=I.Det();
        for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) I(j, k) /= det;
        e = I.Eigenvalues();
        a=e[0];b=e[1];c=e[2];
        eigenvec=I.Eigenvectors(e);
    }
    void GetInertiaTensor(System& S, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec)
    {
		GetInertiaTensor(S.GetNumParts(), S.Parts(), a, b, c, eigenvec);
    }

    void GetInertiaTensor(System& S, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec, Matrix &I)
    {
		GetInertiaTensor(S.GetNumParts(), S.Parts(), a, b, c, eigenvec, I);
    }

    static void CalculateTensor(Matrix& M, const Double_t q, const Double_t s, const Int_t n, Particle *p)
    {
        Int_t i;
        int j,k;
        Double_t a2;
#ifdef USEOPENMP
        int nthreads;
    #pragma omp parallel 
    {
        nthreads=omp_get_num_threads();
    }
        Matrix *Mt=new Matrix[nthreads];
    #pragma omp parallel default(shared) \
    private(i,j,k,a2)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i = 0; i < n; i++)
        {
            a2 = (p[i].X() * p[i].X() + 
                p[i].Y() * p[i].Y() / q / q +
                p[i].Z() * p[i].Z() / s / s);
            //if particle is at 0 then set a2 to a non-zero number
            //if(!finite(1./a2))a2=1.0;
            //if(a2==0)a2=1.0;
            if(a2!=0){
                for (j = 0; j < 3; j++)
                    for (k = 0; k < 3; k++)
#ifdef USEOPENMP
                    Mt[omp_get_thread_num()](j,k)+= p[i].GetPosition(j)*p[i].GetPosition(k) / a2;
#else
                    M(j, k) += p[i].GetPosition(j)*p[i].GetPosition(k) / a2;
#endif
            }
        }
#ifdef USEOPENMP
    }
        for (int l=0;l<nthreads;l++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    M(j,k)+=Mt[l](j,k);
#endif
    }
    static void CalculateTensorWithMass(Matrix& M, const Double_t q, const Double_t s, const Int_t n, Particle *p)
    {
        Int_t i;
        int j,k;
        //Double_t TotalMass=0;
        Double_t a2;
#ifdef USEOPENMP
        int nthreads;
    #pragma omp parallel 
    {
        nthreads=omp_get_num_threads();
    }
        Matrix *Mt=new Matrix[nthreads];
    #pragma omp parallel default(shared) \
    private(i,j,k,a2)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i = 0; i < n; i++)
        {
            a2 = (p[i].X() * p[i].X() + 
                p[i].Y() * p[i].Y() / q / q +
                p[i].Z() * p[i].Z() / s / s);
            //if particle is at 0 then set a2 to a non-zero number
            //if(!finite(1./a2))a2=1.0;
            if(a2!=0){
                for (j = 0; j < 3; j++)
                    for (k = 0; k < 3; k++)
#ifdef USEOPENMP
                    Mt[omp_get_thread_num()](j,k)+= p[i].GetMass()*p[i].GetPosition(j)*p[i].GetPosition(k) / a2;
#else
                    M(j, k) += p[i].GetMass()*p[i].GetPosition(j)*p[i].GetPosition(k) / a2;
#endif
            }
        }
#ifdef USEOPENMP
    }
        for (int l=0;l<nthreads;l++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    M(j,k)+=Mt[l](j,k);
#endif
    }
    static void CalculateTensorWithMassInDenom(Matrix& M, const Double_t q, const Double_t s, const Int_t n, Particle *p)
    {
        Int_t i;
        int j,k;
        Double_t a2;
#ifdef USEOPENMP
        int nthreads;
    #pragma omp parallel 
    {
        nthreads=omp_get_num_threads();
    }
        Matrix *Mt=new Matrix[nthreads];
    #pragma omp parallel default(shared) \
    private(i,j,k,a2)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i = 0; i < n; i++)
        {
            a2 = p[i].GetMass()*(p[i].X() * p[i].X() + 
                p[i].Y() * p[i].Y() / q / q +
                p[i].Z() * p[i].Z() / s / s);
            //if particle is at 0 then set a2 to a non-zero number
            //if(!finite(1./a2))a2=1.0;
            if(a2!=0){
                for (j = 0; j < 3; j++)
                    for (k = 0; k < 3; k++)
#ifdef USEOPENMP
                    Mt[omp_get_thread_num()](j,k)+= p[i].GetPosition(j)*p[i].GetPosition(k) / a2;
#else
                    M(j, k) += p[i].GetPosition(j)*p[i].GetPosition(k) / a2;
#endif
            }
        }
#ifdef USEOPENMP
    }
        for (int l=0;l<nthreads;l++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    M(j,k)+=Mt[l](j,k);
#endif
    }


    static void RotateParticles(const Int_t n, Particle *p, const Matrix& M, const Coordinate& e)
    {
        Matrix U= M.Eigenvectors(e);
        Int_t i;
        int j;
        Coordinate temp;
        // Now rotate the vector

#ifdef USEOPENMP
    #pragma omp parallel default(shared) \
    private(i,j,temp)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i=0; i<n; i++)
        {
            temp[0]=temp[1]=temp[2]=0.;
            for (j=0; j<3; j++)
            {
                temp[0] += U(0, j) * (p[i].GetPosition(j));
                temp[1] += U(1, j) * (p[i].GetPosition(j));
                temp[2] += U(2, j) * (p[i].GetPosition(j));
            }
            p[i].SetPosition(temp[0], temp[1], temp[2]);
        }
#ifdef USEOPENMP
    }
#endif
    }

    ///rotate coordinate system
    inline Coordinate Rotate(const Matrix R, Coordinate x)
    {
        Coordinate xp(0,0,0);
        for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            xp[i]+=R(i,j)*x[j];
        return xp;
    }

}
