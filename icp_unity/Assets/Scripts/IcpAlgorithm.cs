using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// using CenterSpace.NMath.Core;
// using CenterSpace.NMath.Matrix;

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
        Matrix4x4 H = AA.transpose * BB;
        // var svd = new DoubleSVDecomp(H);
        // U, S, Vt = np.linalg.svd(H)
        // R = np.dot(Vt.T, U.T)

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
