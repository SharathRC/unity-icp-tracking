using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// using alglib;

// using CenterSpace.NMath.Core;
// using CenterSpace.NMath.Matrix;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;

using Accord.Imaging.Filters;

public class IcpAlgorithm : MonoBehaviour
{
    private IcpUtilities icp_utilities = new IcpUtilities();
    private Utilities _utilities = new Utilities();

    internal Matrix<double> present_transform;
    internal Matrix<double> translate_mat;
    internal Matrix<double> rotation_mat;


    void best_fit_transform(Matrix4x4 A, Matrix4x4 B)
    {
        // int m = 4; //Matrix4x4 size
        // centroid_A = np.mean(A, axis=0)
        // centroid_B = np.mean(B, axis=0)
        // AA = A - centroid_A
        // BB = B - centroid_B
        // rotation matrix
        // H = np.dot(AA.T, BB)
        // U, S, Vt = np.linalg.svd(H)
        // R = np.dot(Vt.T, U.T)
        // # special reflection case
        // if np.linalg.det(R) < 0:
        //      Vt[m-1,:] *= -1
        //      R = np.dot(Vt.T, U.T)
        // # translation
        // t = centroid_B.T - np.dot(R,centroid_A.T)
        // # homogeneous transformation
        // T = np.identity(m+1)
        // T[:m, :m] = R
        // T[:m, m] = t

        Vector4 centroid_A = icp_utilities.matrix4x4_mean(A, 0);
        Vector4 centroid_B = icp_utilities.matrix4x4_mean(B, 0);

        Matrix4x4 AA = icp_utilities.matrix4x4_vector_sub(A, centroid_A);
        Matrix4x4 BB = icp_utilities.matrix4x4_vector_sub(B, centroid_B);


        Matrix4x4 H = AA.transpose * BB;

        Matrix<double> new_H = DenseMatrix.OfArray(new double[,] {
                                                                {H[0,0], H[0,1], H[0,2], H[0,3]},
                                                                {H[1,0], H[1,1], H[1,2], H[1,3]},
                                                                {H[2,0], H[2,1], H[2,2], H[2,3]},
                                                                {H[3,0], H[3,1], H[3,2], H[3,3]}   });
        
        var svd = new_H.Svd(true);

        var U = svd.U;
        var S = svd.W;
        var Vt = svd.VT;

        var R = Vt.Transpose().Multiply(U.Transpose());

        if(R.Determinant() < 0)
        {
            Matrix<double> cust_identity = DenseMatrix.OfArray(new double[,] {
                                                                            {1, 0, 0, 0},
                                                                            {0, 1, 0, 0},
                                                                            {0, 0, 1, 0},
                                                                            {0, 0, 0, -1}   });
            Vt = (Vt.Transpose().Multiply(cust_identity)).Transpose();
            R = Vt.Transpose().Multiply(U.Transpose());
        }

        Vector<double> cent_A = new DenseVector(new[] { Convert.ToDouble(centroid_A[0]), 
                                                        Convert.ToDouble(centroid_A[1]), 
                                                        Convert.ToDouble(centroid_A[2]), 
                                                        Convert.ToDouble(centroid_A[3]) });

        Vector<double> cent_B = new DenseVector(new[] { Convert.ToDouble(centroid_B[0]), 
                                                        Convert.ToDouble(centroid_B[1]), 
                                                        Convert.ToDouble(centroid_B[2]), 
                                                        Convert.ToDouble(centroid_B[3]) });
        
        var cent_A_mat = cent_A.ToColumnMatrix();
        var cent_B_mat = cent_B.ToColumnMatrix();

        var t = cent_B_mat - R.Multiply(cent_A_mat);

        var T = R.Append(t);

        Vector<double> last_col = new DenseVector(new[] { 0.0, 0.0, 0.0, 1.0 });
        var last_col_mat = last_col.ToRowMatrix();
        T = T.Append(last_col_mat);
        
        present_transform = T;
        translate_mat = t;
        rotation_mat = R;

        // return T, R, t
    }

    void nearest_neighbor( Matrix<double> src,  Matrix<double>dst)
    {
        
    }
}
