using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class IcpUtilities : MonoBehaviour
{
    internal Vector4 matrix4x4_mean(Matrix4x4 x, int axis)
    {   
        // List<float> mean_values = new List<float>();
        Vector4 mean_values = new Vector4();
        if (axis == 1)
        {
            x = x.transpose;
        }
        for(var i = 0; i < 4; i++)
        {
            Vector4 vals = x.GetColumn(i);
            float mean_val = (vals[0] + vals[1] + vals[2] + vals[3])/4;
            // mean_values.Add(mean_val);
            mean_values[i] = mean_val;
        }

        return mean_values;
    }

    internal Matrix4x4 matrix4x4_vector_sub(Matrix4x4 x, Vector4 v)
    {   
        Matrix4x4 new_mat = new Matrix4x4();
        for(var i = 0; i<4; i++)
        {
            Vector4 col = x.GetColumn(i);
            Vector4 new_col = col - v;
            new_mat.SetColumn(i, new_col);
        }

        return new_mat;
    }


}
