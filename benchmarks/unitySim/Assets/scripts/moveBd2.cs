using UnityEngine;

public class moveBd2 : MonoBehaviour
{

    public kinematics kin;  //access to kinematics
    private int bodyID; //body id
    void Start()
    {
        bodyID = 3;
    }

    void FixedUpdate()
    {
        Vector3 r = kin.r(bodyID);
        Quaternion p = kin.p(bodyID);
        transform.position = r;
        transform.rotation = p;
        //transform.localPosition = r;
        //transform.localRotation = p;
        //Transform.SetPositionAndRotation(r,p) from future import
    }



}
