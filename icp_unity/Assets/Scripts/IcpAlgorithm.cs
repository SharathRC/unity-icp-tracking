using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// using alglib;

// using CenterSpace.NMath.Core;
// using CenterSpace.NMath.Matrix;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;

public class IcpAlgorithm : MonoBehaviour
{
    private IcpUtilities icp_utilities = new IcpUtilities();
    private Utilities _utilities = new Utilities();


    void best_fit_transform(Matrix4x4 A, Matrix4x4 B)
    {
        int m = 4; //Matrix4x4 size
        Vector4 centroid_A = icp_utilities.matrix4x4_mean(A, 0);
        Vector4 centroid_B = icp_utilities.matrix4x4_mean(B, 0);

        Matrix4x4 AA = icp_utilities.matrix4x4_vector_sub(A, centroid_A);
        Matrix4x4 BB = icp_utilities.matrix4x4_vector_sub(B, centroid_B);

        // rotation matrix
        // H = np.dot(AA.T, BB)
        // U, S, Vt = np.linalg.svd(H)
        // R = np.dot(Vt.T, U.T)

        Matrix4x4 H = AA.transpose * BB;

        Matrix<double> new_mat = DenseMatrix.OfArray(new double[,] {
                                                                {H[0,0], H[0,1], H[0,2], H[0,3]},
                                                                {H[1,0], H[1,1], H[1,2], H[1,3]},
                                                                {H[2,0], H[2,1], H[2,2], H[2,3]},
                                                                {H[3,0], H[3,1], H[3,2], H[3,3]}   });
        
        var svd = new_mat.Svd(true);

        var U = svd.U;
        var S = svd.W;
        var Vt = svd.VT;

        // # special reflection case
        // if np.linalg.det(R) < 0:
        // Vt[m-1,:] *= -1
        // R = np.dot(Vt.T, U.T)

        // # translation
        // t = centroid_B.T - np.dot(R,centroid_A.T)

        // # homogeneous transformation
        // T = np.identity(m+1)
        // T[:m, :m] = R
        // T[:m, m] = t

        // return T, R, t
    }
}
