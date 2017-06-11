using UnityEngine;
using System.Collections;

public class orientcamera : MonoBehaviour {

    // Use this for initialization
    public Transform ground;

	void Start () {
        transform.position = new Vector3(-10, 0, 0);
        Vector3 worldUP = new Vector3(0, 0, 1);
        transform.LookAt(ground, worldUP);
	
	}
	
	// Update is called once per frame
	void Update () {
	
	}
}
