using UnityEngine;

public class moveBd : MonoBehaviour {

    public kinematics kin;  //access to kinematics
    public int bodyID;      //body ID number
    void Start ()  
    {
        
    }

    void FixedUpdate()
    {
        Vector3    r = kin.r(bodyID);
        Quaternion p = kin.p(bodyID);
        transform.position = r;
        transform.rotation = p;
        //transform.localPosition = r;
        //transform.localRotation = p;
        //Transform.SetPositionAndRotation(r,p) from future import
    }



}
