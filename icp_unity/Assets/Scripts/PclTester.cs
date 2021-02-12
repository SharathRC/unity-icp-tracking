using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using PclSharp;
using PclSharp.Std;
using PclSharp.Registration;
using PclSharp.IO;
using PclSharp.Struct;

public class PclTester : MonoBehaviour
{
    private IcpUtilities icp_utilities = new IcpUtilities();
    private Utilities _utilities = new Utilities();

    // internal Registration icp = new Registration();

    // internal PclSharp.PointCloud  cloud_in = new pcl.PointCloud(5,1);

    PointCloud<PointXYZ> cloud_in = new PointCloud<PointXYZ>(5,1);


}
