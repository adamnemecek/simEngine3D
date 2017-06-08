using UnityEngine;

public class moveBd1 : MonoBehaviour {

    public kinematics kin;  //access to kinematics
    private int bodyID = 2; //body id
    void Start ()  
    {
        
    }

    void FixedUpdate()
    {
        Vector3    r = kin.r(bodyID);
        //Quaternion p = kin.p(bodyID);
        transform.position = r;
       // transform.rotation = p;
        //Transform.SetPositionAndRotation(r,p) from future import
    }



}
