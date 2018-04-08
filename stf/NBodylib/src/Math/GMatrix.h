/*! \file GMatrix.h
 *  \brief this file defines the \ref Math::GMatrix class that is part of the Math namespace
 */

#ifndef GMATRIX_H
#define GMATRIX_H

#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <Precision.h>
#include <Coordinate.h>
#include <Matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;
namespace Math
{

/*! \class Math::GMatrix
    \brief General n x m matrix
    General matrix class that eventually should replace \ref Math::Matrix class.
    will have to adjust error checks for this class and adjust inverse
    and make most calls inline to increase speed.
*/
class GMatrix
{
    private:

    int row,col;
    Double_t *matrix;

    protected:

    ///useful comparison functions which can be used with qsort
    int compare_Double(const void* a,const void* b) {
        Double_t arg1 = *((Double_t*) a);
        Double_t arg2 = *((Double_t*) b);
        if( arg1 < arg2 ) return -1;
        else if( arg1==arg2 ) return 0;
        else return 1;
    }

    public:

    /// \name Constructors & Destructors
    //@{
    /// Constructor (where data can be is initilized) otherwise set to zero
    GMatrix(int r, int c, Double_t *data=NULL)
    {
        if (r*c>0){
            row=r;col=c;
            matrix=new Double_t[row*col];
            if (data!=NULL){
                for (int i = 0; i < row; i++)
                    for (int j = 0; j < col; j++)
                        matrix[i*col+j] = data[i*col+j];
            }
            else{
                for (int i = 0; i < row; i++)
                    for (int j = 0; j < col; j++)
                        matrix[i*col+j] = 0.;
            }
        }
        else {
            row=0;col=0;matrix=NULL;
            printf("Error, net size 0 with row,col=%d,%d\n",row,col);
        }
    }

    /// Constructor based on Coordinate.
    GMatrix(const Coordinate &c)
    {
        row=3;col=1;
        matrix=new Double_t[row*col];
        for (int i = 0; i < row; i++) matrix[i] = c[i];
    }

    /// Constructor Matrix
    GMatrix(const Matrix &m)
    {
        row=3;col=3;
        matrix=new Double_t[row*col];
        for (int i=0; i<row; i++)
            for (int j=0;j<col;j++) matrix[i*col+j] = m(i,j);
    }

    /// Copy Constructor
    GMatrix(const GMatrix &m)
    {
        row=m.row;col=m.col;
        matrix=new Double_t[row*col];
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                matrix[i*col+j] = m.matrix[i*col+j];
    }

    /// Identity Matrix of dimension N
    GMatrix(int N)
    {
        row=N;col=N;
        matrix=new Double_t[row*col];
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                matrix[i*col+j] = (i==j);
    }

    ~GMatrix(){
        if (row>0&&col>0) delete[] matrix;
    }
    //@}

    /// \name Get & Set functions
    //@{
    ///get row dimensions
    int Row()const {return row;}
    ///get col dimensions
    int Col()const {return col;}
    ///get rank
    int Rank()const {return min(row,col);}
    //@}

    /// \name Overloaded Operators
    //@{
    /// How to access a matrix element: M(i, j)
    Double_t& operator () (int i, int j) {
        //if (j<col&&i<row)
        return matrix[j + col * i];
        //else {printf("Error, (i,j)=(%d, %d) exceed matrix dimesions of %d %d \n", i,j,row,col); return 0;}
    }
    const Double_t& operator() (int i, int j) const
    {
        //if (j<col&&i<row)
        return matrix[j + col * i];
        //else {printf("Error, (i,j)=(%d, %d) exceed matrix dimesions of %d %d \n", i,j,row,col); return 0;}
    }

    /// Matrix addition
    GMatrix operator + (const GMatrix& m)
    {
        if (row==m.row&&col==m.col){
            Double_t data[row*col];
            for (int i = 0; i < row; i++)
                for (int j = 0; j < col; j++)
                    data[i*col+j] = matrix[i*col+j]+m.matrix[i*col+j];
            return GMatrix(row,col,data);
        }
        else {
            printf("Error, incompatible dimensions (%d,%d) + (%d,%d)\n",row,col,m.row,m.col);
            GMatrix(0,0);
        }
    }

    /// Matrix subtraction
    GMatrix operator - (const GMatrix& m)
    {
        if (row==m.row&&col==m.col){
            Double_t data[row*col];
            for (int i = 0; i < row; i++)
                for (int j = 0; j < col; j++)
                    data[i*col+j] = matrix[i*col+j]-m.matrix[i*col+j];
            return GMatrix(row,col,data);
        }
        else {
            printf("Error, incompatible dimensions (%d,%d) + (%d,%d)\n",row,col,m.row,m.col);
            GMatrix(0,0);
        }
    }

    /// Multiply by a scalar: a * M
    friend GMatrix operator * (Double_t a, const GMatrix& m);

    /// Multiply by a scalar: M * a
    GMatrix operator * (Double_t a)
    {
        Double_t data[row*col];
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                data[i*col+j] = matrix[i*col+j]*a;
        return GMatrix(row,col,data);
    }

    /// Multiply by another Gmatrix
    GMatrix operator * (const GMatrix& m)
    {
        if (col==m.row){
            Double_t data[row*m.col];
            for (int i = 0; i < row; i++)
                for (int j = 0; j < m.col; j++){
                    data[i*m.col+j]=0.;
                    for (int k=0;k<col;k++)
                        data[i*m.col+j]+=matrix[i*col+k]*m.matrix[k*m.col+j];
                }
            return GMatrix(row,m.col,data);
        }
        else {
            printf("Error, incompatible dimensions (%d,%d) + (%d,%d)\n",row,col,m.row,m.col);
            return GMatrix(0,0);
        }

    }
    /// Asignment constructor
    GMatrix& operator = (const GMatrix m)
    {
        delete[] matrix;
        row=m.row;
        col=m.col;
        matrix=new Double_t[row*col];
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                matrix[i*col+j] = m.matrix[i*col+j];
        return *this;
    }
    //@}

    /// \name Mathematical Operations
    //@{

    /// Trace the matrix.
    Double_t Trace() const
    {
        Double_t trace=1.0;
        for (int i = 0; i < row; i++)trace*=matrix[i*col+i];
        return trace;
    }

    /// return the Transpose of the matrix.
    GMatrix Transpose() const
    {
        Double_t data[col*row];
        for (int i = 0; i < col; i++)
            for (int j = 0; j < row; j++)
                data[i*row+j]=(*this)(i,j);
        return GMatrix(col,row,data);
    }

    /// Transpose matrix in place (changes the matrix permanently!)
    void TransposeInPlace()
    {
        Double_t data[col*row];
        int temp;
        for (int i = 0; i < col; i++)
            for (int j = 0; j < row; j++)
                data[i*row+j]=(*this)(i,j);
        temp=row;row=col;col=temp;
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                matrix[i*col+j]=data[i*col+j];
    }
    ///Get the diagonals
    GMatrix Diag(){
        GMatrix D(row,col);
        for (int i = 0; i < col; i++) D(i,i)=matrix[i*col+i];
        return D;
    }
    /// Return Pivot matrix P of A, where P*A=B, such that B(i,i)!=0
    GMatrix Pivot() const
    {
        Double_t a[row*col];
        GMatrix P(row,col);
        for (int j=0;j<row;j++) P(j,j)=1;
        if (Trace()==0) {
            for (int j=0;j<row;j++) for (int k=0;k<row;k++) a[j*col+k]=matrix[j*col+k];
            for (int j=0;j<row;j++) if (a[j*col+j]*a[j*col+j]<1e-32) {
                int k=j+1,ipivot;
                if (k==row)k=0;
                ipivot=(P(k,k)==1);
                while (!(a[j*col+k]*a[j*col+k]>1e-32&&ipivot)){
                    k++;
                    if (k==row)k=0;
                    ipivot=(P(k,k)==1);
                }
                P(j,j)=0;P(k,k)=0;
                P(j,k)=1;P(k,j)=1;
                for (int l=0;l<row;l++) {
                    Double_t temp=a[l*col+j];
                    a[l*col+j]=a[l*col+k];
                    a[l*col+k]=temp;
                }
            }
        }
        return P;
    }
    /// Return Pivot matrix P of A, where P*A=B, such that B(i,i)!=0
    GMatrix PivotInPlace()
    {
        GMatrix P(row,col);
        for (int j=0;j<row;j++) P(j,j)=1;
        if (Trace()==0) {
            for (int j=0;j<row;j++) if (matrix[j*col+j]*matrix[j*col+j]<1e-32) {
                int k=j+1,ipivot;
                if (k==row)k=0;
                ipivot=(P(k,k)==1);
                while (!(matrix[j*col+k]*matrix[j*col+k]>1e-32&&ipivot)){
                    k++;
                    if (k==row)k=0;
                    ipivot=(P(k,k)==1);
                }
                P(j,j)=0;P(k,k)=0;
                P(j,k)=1;P(k,j)=1;
                for (int l=0;l<row;l++) {
                    Double_t temp=matrix[l*col+j];
                    matrix[l*col+j]=matrix[l*col+k];
                    matrix[l*col+k]=temp;
                }
            }
        }
        return P;
    }

    /// return a SubMatrix of the matrix from row ra to  row rb and col ca to cb
    GMatrix SubMatrix(int ra, int rb, int ca, int cb) const
    {
        int r=(rb-ra+1), c=(cb-ca+1);
        Double_t data[r*c];
        for (int i = ra; i <= rb; i++)
            for (int j = ca; j <=cb; j++)
                data[(i-ra)*c+(j-ca)]=(*this)(i,j);
        return GMatrix((rb-ra+1),(cb-ca+1),data);
    }

    /// return a submatrix ignore column r, c
    GMatrix ExtractMatrix(int r, int c) const
    {
        Double_t data[(row-1)*(col-1)];
        int tempi=0;
        for (int i = 0; i < row; i++)
        for (int j = 0; j <col; j++)
            if (i!=r&&j!=c) data[tempi++]=(*this)(i,j);
        return GMatrix((row-1),(col-1),data);
    }

    /// Check if matrix square symmetric
    int isSymmetric() const
    {
        if (row!=col) return 0;
        for (int i=0;i<row;i++)
        for (int j=i+1;j<col-1;j++)
            if (((*this)(i,j)-(*this)(j,i))*((*this)(i,j)-(*this)(j,i))/(*this)(i,j)/(*this)(i,j)>1e-8) return 0;
        return 1;
    }
    /// Check if matrix is zero matrix
    int isZero() const
    {
        Double_t sm=0;
        for (int i=0;i<row;i++)
            for (int j=0;j<col;j++)
                sm+=matrix[i*row+j]*matrix[i*row+j];
        if (sm==0.0) return 1;
        else return 0;
    }

    /// Determinant of matrix, for now recursively implemented.
    Double_t Det() const
    {
        Double_t total =0.;
        int tempi;
        if (row==col) {
            if (row==2) return matrix[0]*matrix[3]-matrix[1]*matrix[2];
            else {
                for (int i = 0; i < col; i++) {
                    Double_t data[(row-1)*(row-1)];
                    tempi=0;
                    for (int j=1;j<row;j++)
                        for (int k=0;k<col;k++)
                            if (k!=i) data[tempi++]=(*this)(j,k);
                    total+=matrix[i]*(GMatrix(row-1,col-1,data)).Det();
                }
                return total;
            }
        }
        else {
            printf("not a square matrix %d, %d returning 0\n",row,col);
            return 0.;
        }
    }

    /*! \brief Lower-Upper decompistion

        LU decomposition of matrix stored in such a fashion that
        matrix return is (for 3x3):
        \f[
        \left(\begin{array}{ccc}
            L11 & U12 &U13\\
            L21 & L22 & U23\\
            L31 & L32 & L33
        \end{array}\right)
        \f]
        L is lower triangular matrix, U is upper unit triangular matrix. \n
        Crout's algorithm is used to determine the L and U matrices. System of equations to solve are : \n
        \f[
        \begin{array}{lcl}
            i<j : & l_{i1}*u_{1j}+l_{i2}*u_{2j}+...+l_{ii}*u_{ij}&=a_{ij} \\
            i=j : & l_{i1}*u_{1j}+l_{i2}*u_{2j}+...+l_{ii}*u_{jj}&=a_{ij} \\
            i>j : & l_{i1}*u_{1j}+l_{i2}*u_{2j}+...+l_{ij}*u_{ij}&=a_{ij}
            \end{array}
        \f]
        with condition that \f$ l_{ii}=1 \f$. \n
        If matrix singular and no decomposition, print error and return zero matrix.
        Here matrix is singular if diagonal element is zero. \n
        Must alter this to pivot method so that subroutine reorder matrix if necessary and diagonals
        can be zero.
    */
    GMatrix LUDecomp() const
    {
        if (row!=col) {
            printf("only implemented for square, returning zeros\n");
            return GMatrix(row,col);
        }
        else {
            Double_t A[row*col];
            Int_t i, j, k, p, p_k, p_row, p_col;
            for (j=0;j<row;j++) for (k=0;k<col;k++) A[j*col+k]=matrix[j*col+k];
            for (k = 0, p_k = 0; k < row; p_k += row, k++) {
                for (i = k, p_row = p_k; i < row; p_row += row, i++) {
                    for (p = 0, p_col = 0; p < k; p_col += row, p++)
                        A[p_row + k] -= A[p_row + p] * A[p_col + k];
                }
                if ( A[p_k + k] == 0.0 ) {
                    printf("Singular matrix, returning zero matrix\n");
                    return GMatrix(row,col);
                }
                for (j = k+1; j < row; j++) {
                    for (p = 0, p_col = 0; p < k; p_col += row,  p++)
                        A[p_k + j] -= A[p_k + p] * A[p_col + j];
                    A[p_k + j] /= A[p_k + k];
                }
            }
            return GMatrix(row,col,A);
        }
    }
    /*! \brief Calculate the Inverse of the lower triangular matrix

        This routine calculates the inverse of the lower triangular matrix L.
        The superdiagonal part of the matrix is not addressed.
        The algorithm follows: \n
        Let M be the inverse of L, then \f$ L M = I \f$, \n
        \f$ M[i][i] = 1.0 / L[i][i] \f$ for i = 0, ..., row-1, and \n
        \f$ M[i][j] = -[(L[i][j] M[j][j] + ... + L[i][i-1] M[i-1][j])] / L[i][i], \f$
        for i = 1, ..., row-1, j = 0, ..., i - 1.
    */
    GMatrix Lower_Triangular_Inverse()
    {
        Double_t L[row*col];
        int i, j, k, p_i, p_j, p_k;
        Double_t sum;
        for (j=0;j<row;j++) for (k=0;k<=j;k++) L[j*col+k]=matrix[j*col+k];
        for (j=0;j<row;j++) for (k=j+1;k<col;k++) L[j*col+k]=0;
        //Invert the diagonal elements of the lower triangular matrix L.
        for (k = 0, p_k = 0; k < row; p_k += (row + 1), k++) {
            if (L[p_k] == 0.0) return GMatrix(row,row);//error, returning zero matrix
            else L[p_k] = 1.0 / L[p_k];
        }
        //Invert the remaining lower triangular matrix L row by row.
        for (i = 1, p_i = 0 + row; i < row; i++, p_i += row) {
            for (j = 0, p_j = 0; j < i; p_j += row, j++) {
                sum = 0.0;
                for (k = j, p_k = p_j; k < i; k++, p_k += row)
                    sum += L[p_i + k] * L[p_k + j];
                L[p_i + j] = - L[p_i + i] * sum;
            }
        }
        return GMatrix(row,row,L);
    }
    /*! \brief Calculate the Inverse of the lower triangular matrix

        This routine calculates the inverse of the unit upper triangular matrix U.
        The subdiagonal part of the matrix is not addressed.
        The diagonal is assumed to consist of 1's and is not addressed.
        The algorithm follows: \n
        Let M be the inverse of U, then \f$ U M = I \f$, \n
        \f$ M[i][j] = -( U[i][i+1] M[i+1][j] + ... + U[i][j] M[j][j] ), \f$
        for i = row-2, ... , 0,  j = row-1, ..., i+1.
    */
    GMatrix Unit_Upper_Triangular_Inverse()
    {
        Double_t U[row*col];
        int i, j, k, p_i, p_k;
        for (j=0;j<row;j++) for (k=j+1;k<col;k++) U[j*col+k]=matrix[j*col+k];
        for (j=0;j<row;j++) for (k=0;k<=j;k++) U[j*col+k]=(j==k);
        //         Invert the superdiagonal part of the matrix U row by row where
        //         the diagonal elements are assumed to be 1.0.
        for (i = row - 2, p_i = 0 + row * (row - 2); i >=0; p_i -= row, i-- ) {
            for (j = row - 1; j > i; j--) {
                U[p_i + j] = -U[p_i + j];
                for (k = i + 1, p_k = p_i + row; k < j; p_k += row, k++ )
                    U[p_i + j] -= U[p_i + k] * U[p_k + j];
            }
        }
        return GMatrix(row,row,U);
    }

    /// \brief Adjugate (classical adjoint) of matrix
    /// This is simply the transpose of the cofactor of the matrix
    GMatrix Adjugate() const
    {
        GMatrix cofactorT(row,col);
        for (int i=0;i<row;i++)
            for (int j=0;j<col;j++)
                cofactorT(j,i)=(ExtractMatrix(i,j)).Det()*pow(-1.,i+j);
        return cofactorT;
    }

    /// Inverse of matrix
    GMatrix Inverse() const
    {
        if (row!=col) {
            printf("only implemented for square, returning zeros\n");
            return GMatrix(row,col);
        }
        GMatrix LU(LUDecomp());
        Double_t d=1.0;
        for (int i=0;i<row;i++) d*=LU(i,i);
        if (d  == 0.0) return GMatrix(row, col);
        else {
            GMatrix Uinv=LU.Unit_Upper_Triangular_Inverse();
            GMatrix Linv=LU.Lower_Triangular_Inverse();
            return Uinv*Linv;
        }
    }

    /// Inverse of matrix using a pivot (needed in case diagonals are zero)
    GMatrix InversewithPivot() const
    {
        if (row!=col) {
            printf("only implemented for square, returning zeros\n");
            return GMatrix(row,col);
        }
        GMatrix A(*this);
        GMatrix Pivot(A.PivotInPlace());
        GMatrix LU(A.LUDecomp());
        Double_t d=1.0;
        for (int i=0;i<row;i++) d*=LU(i,i);
        if (d  == 0.0) return GMatrix(row, col);
        else {
            GMatrix Uinv=LU.Unit_Upper_Triangular_Inverse();
            GMatrix Linv=LU.Lower_Triangular_Inverse();
            return Uinv*Linv*Pivot;
        }
    }

    /// helpful rotation definition used in \ref Math::GMatrix::Jacobi (see Numerical recipes)
    inline void ROTATE(GMatrix &m, int i, int j, int k, int l, Double_t sinrot, Double_t tanrot) const {
        Double_t tempa=m(i,j),tempb=m(k,l);
        m(i,j)=tempa-sinrot*(tempb+tempa*tanrot);m(k,l)=tempb+sinrot*(tempa-tempb*tanrot);
    }
    /// Applies jacobi transformations to a copy of a matrix and returns effectively the eigenvalues of the matrix
    GMatrix Jacobi(Double_t tol=1e-2) const
    {
        GMatrix m(*this);
        Double_t eigen[row], b[row], z[row];
        GMatrix Prot(row);//identity
        int nrot=0;
        Double_t crot,srot,trot2,thres,t,theta,sm,temp1,temp2;

        //initialize vectors to diagonal elements of matrix
        for (int i=0;i<row;i++){ eigen[i]=b[i]=m(i,i);z[i]=0.;}

        do {
            //first check to see if off diagonal elements are zero
            // could check to within some precision relative to trace
            sm=0.;
            for (int ip=0;ip<row-1;ip++)
            for (int iq=ip+1;iq<col;iq++)
                sm+=sqrt(m(ip,iq)*m(ip,iq));
            if (sm==0.) return GMatrix(row,1,eigen);
            //set the threshold to be 0 until several rotations have been applied
            nrot<4? thres=tol*sm/(Double_t)(row*row): thres=0.;
            //now start applying Jacobi rotations
            for (int ip=0;ip<row-1;ip++){
                for (int iq=ip+1;iq<col;iq++){
                    //get absolute value of off diagonal element and check if it is small in comparison to
                    //threshold factor of 1/5 is from numerical recipes
                    temp1=sqrt(m(ip,iq)*m(ip,iq))*0.2;
                    //if smaller than tolerance set to 0. after for sweeps and skip the rotation.
                    if (nrot > 4 && (sqrt(eigen[ip]*eigen[ip])+temp1)==(sqrt(eigen[ip]*eigen[ip])+tol)
                        && (sqrt(eigen[iq]*eigen[iq])+temp1)==(sqrt(eigen[iq]*eigen[iq])+tol))
                        m(ip,iq)=0.0;
                    else if (thres<temp1*5.0) {
                        temp2=eigen[iq]-eigen[ip];
                        //find the solution to t^2+2*t*theta-1=0
                        //where theta=(cos^2-sin^2)/(2*s*c)=(m(iq,iq)-m(ip,ip))/(2*m(ip,iq))
                        //note that nr has ad hoc check if theta^2 would overflow the computer
                        /*if ((sqrt(temp2*temp2)+temp1) == sqrt(temp2*temp2)) t=(m(ip,iq))/temp2;
                        else {
                            theta=0.5*temp2/m(ip,iq);
                            t=1.0/(sqrt(theta*theta)+sqrt(1.0+theta*theta));
                            if (theta < 0.0) t = -t;
                        }*/
                        //instead I just determine the value of t
                        theta=0.5*temp2/m(ip,iq);
                        t=1.0/(sqrt(theta*theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                        //calculate cos and sin of rotation angle and tan(rot/2)
                        crot=1.0/sqrt(1.+t*t); srot=t*crot; trot2=srot/(1.0+crot);
                        //not alter values in matrix
                        temp2=t*m(ip,iq);
                        z[ip]-= temp2; z[iq] += temp2;
                        eigen[ip]-= temp2; eigen[iq] += temp2;
                        m(ip,iq)=0.0;
                        //now rotate the ipth values for 0<=j<ip
                        for (int j=0;j<ip;j++) {ROTATE(m,j,ip,j,iq,srot,trot2);}
                        //rotate the values for ip<j<iq
                        for (int j=ip+1;j<iq;j++) {ROTATE(m,ip,j,j,iq,srot,trot2);}
                        //rotate the values for iq<j<=row
                        for (int j=iq+1;j<row;j++) {ROTATE(m,ip,j,iq,j,srot,trot2);}
                        //store the rotations in Prot
                        for (int j=0;j<row;j++) {ROTATE(Prot,j,ip,j,iq,srot,trot2);}
                    }
                }
            }
            //now update vectors used in rotation and reinitialize b
            for (int ip=0;ip<row;ip++) {
                b[ip] += z[ip];
                eigen[ip]=b[ip];
                z[ip]=0.0;
            }
        }while (nrot++<10);
        printf("Too many iterations in routine jacobi, giving zeros\n");
        for (int i=0;i<row;i++) eigen[i]=0.;
        return GMatrix(row,1,eigen);
    }
    /// same as \ref Math::GMatrix::Jacobi above with different arguments and also returns Prot
    void Jacobi(GMatrix &eigenval, GMatrix &Prot, Double_t tol=1e-2) const
    {
        GMatrix m((*this));
        Double_t eigen[row], b[row], z[row];
        for (int i=0;i<row;i++) Prot(i,i)=1.;
        int nrot=0;
        Double_t crot,srot,trot2,thres,t,theta,sm,temp1,temp2;

        //initialize vectors to diagonal elements of matrix
        for (int i=0;i<row;i++){eigen[i]=b[i]=m(i,i);z[i]=0.;}

        do {
            //first check to see if off diagonal elements are zero
            // could check to within some precision relative to trace
            sm=0.;
            for (int ip=0;ip<row-1;ip++)
            for (int iq=ip+1;iq<col;iq++)
                sm+=sqrt(m(ip,iq)*m(ip,iq));
            if (sm==0.) break;
            //set the threshold to be 0 until several rotations have been applied
            nrot<4? thres=tol*sm/(Double_t)(row*row): thres=0.;
            //now start applying Jacobi rotations
            for (int ip=0;ip<row-1;ip++){
                for (int iq=ip+1;iq<col;iq++){
                    //get absolute value of off diagonal element and check if it is small in comparison to
                    //threshold factor of 1/5 is from numerical recipes
                    temp1=sqrt(m(ip,iq)*m(ip,iq))*0.2;
                    //if smaller than tolerance set to 0. after for sweeps and skip the rotation.
                    if (nrot > 4 && (sqrt(eigen[ip]*eigen[ip])+temp1)==(sqrt(eigen[ip]*eigen[ip])+tol)
                        && (sqrt(eigen[iq]*eigen[iq])+temp1)==(sqrt(eigen[iq]*eigen[iq])+tol))
                        m(ip,iq)=0.0;
                    else if (thres<temp1*5.0) {
                        temp2=eigen[iq]-eigen[ip];
                        //find the solution to t^2+2*t*theta-1=0
                        //where theta=(cos^2-sin^2)/(2*s*c)=(m(iq,iq)-m(ip,ip))/(2*m(ip,iq))
                        //note that nr has ad hoc check if theta^2 would overflow the computer
                        /*if ((sqrt(temp2*temp2)+temp1) == sqrt(temp2*temp2)) t=(m(ip,iq))/temp2;
                        else {
                            theta=0.5*temp2/m(ip,iq);
                            t=1.0/(sqrt(theta*theta)+sqrt(1.0+theta*theta));
                            if (theta < 0.0) t = -t;
                        }*/
                        //instead I just determine the value of t
                        theta=0.5*temp2/m(ip,iq);
                        t=1.0/(sqrt(theta*theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                        //calculate cos and sin of rotation angle and tan(rot/2)
                        crot=1.0/sqrt(1.+t*t); srot=t*crot; trot2=srot/(1.0+crot);
                        //not alter values in matrix
                        temp2=t*m(ip,iq);
                        z[ip]-= temp2; z[iq] += temp2;
                        eigen[ip]-= temp2; eigen[iq] += temp2;
                        m(ip,iq)=0.0;
                        //now rotate the ipth values for 0<=j<ip
                        for (int j=0;j<ip;j++) {ROTATE(m,j,ip,j,iq,srot,trot2);}
                        //rotate the values for ip<j<iq
                        for (int j=ip+1;j<iq;j++) {ROTATE(m,ip,j,j,iq,srot,trot2);}
                        //rotate the values for iq<j<=row
                        for (int j=iq+1;j<row;j++) {ROTATE(m,ip,j,iq,j,srot,trot2);}
                        //store the rotations in Prot
                        for (int j=0;j<row;j++) {ROTATE(Prot,j,ip,j,iq,srot,trot2);}
                    }
                }
            }
            //now update vectors used in rotation and reinitialize b
            for (int ip=0;ip<row;ip++) {
                b[ip] += z[ip];
                eigen[ip]=b[ip];
                z[ip]=0.0;
            }
        }while (nrot++<10);
        if (nrot>=10) {
            printf("Too many iterations in routine jacobi, giving zeros\n");
            for (int i=0;i<row;i++) {
                eigenval(i,0)=0.;
                for (int j=0;j<row;j++)
                Prot(i,j)=(i==j);
            }
        }
        for (int i=0;i<row;i++) eigenval(i,0)=eigen[i];
    }

    /*! \brief Calculate eigenvalues of matrix of rank N. Stored as a GMatrix, e1, e2 ...eN.

        The user must know what to expect from this function, whether real eigenvalues
        or complex ones.
        Have to implement the Francis QR algorithm or use successive Jacobi transformations
        to get eigenvalues and eigenvectors (note Jacobi slower than QR)
        but to ensure real eigenvalues, check that the matrix is symmetric
    */
    GMatrix Eigenvalues(Double_t tol=1e-2) const
    {
        if (isSymmetric()==0||isZero()==1) {
            printf("only implemented for nonzero symmetric matrices, returning null matrix\n");
            return GMatrix(0,0);
        }
        else {
            //using Jacobi transformations
            //issue is then must order eigen values
            GMatrix eigen(Jacobi(tol));
            //now order the eigen vaules e0>e1>e2..>eN
            //qsort(eigen.matrix,row,sizeof(Double_t),compare_Double);
            return eigen;
        }
    }

    // Computes the eigenvectors of the matrix, given the eigenvalues.  Sadly, this function will only work
    // properly on symmetric matrices -- i.e., all real eigenvalues.
    //GMatrix Eigenvectors(const Coordinate& e) const;

    /// Calculate eigenvalues \e and eigenvector of matrix of rank N. Stored as a GMatrix, e1, e2 ...eN.
    void Eigenvalvec(GMatrix &eigenval, GMatrix &eigenvector, Double_t tol=1e-2) const
    {
        if (isSymmetric()==0||isZero()==1) {
            printf("only implemented for nonzero symmetric matrices, returning null matrix\n");
            for (int i=0;i<row;i++) {
                eigenval(i,0)=0.;
                for (int j=0;j<row;j++)
                    eigenvector(i,j)=0.;
            }
        }
        else {
            //using Jacobi transformations
            Jacobi(eigenval,eigenvector,tol);
        }
    }

    //@}
};

}
#endif
