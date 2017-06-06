using UnityEngine;

public class kinematics : MonoBehaviour
{

    public TextAsset csv; //attach the q data csv manually
    public string[,] txt_q; //output of split csvGrid
    public float[,] q;
    // Use this for initialization
    void Start()  //could use awake? 
    {
        txt_q = CSVReader.SplitCsvGrid(csv.text);
        q = str2floatArray(txt_q);
        //Debug.Log(txt_q.GetUpperBound(0));
        //Debug.Log(txt_q.GetUpperBound(1));
        //CSVReader.DebugOutputGrid(CSVReader.SplitCsvGrid(csv.text));
        Debug.Log(q[2, 0]);



    }


    //transposes txt_q and converts each element to a float
    static public float[,] str2floatArray(string[,] txt_q)
    {
        float[,] q = new float[txt_q.GetUpperBound(1), txt_q.GetUpperBound(0)]; //this could be wrong, look into further
        for (int row = 1; row < txt_q.GetUpperBound(1) - 1; row += 1) // x = 1 to get rid of row of labels!
        {
            for (int col = 0; col < txt_q.GetUpperBound(0); col += 1)
            {
                string element = txt_q[col, row];
                if (element != null && element != "")  //protect parse from null
                {
                    q[row, col] = float.Parse(txt_q[col, row]);
                }
            }
        }
        return q;


    }

    //returns the position vector (rx,ry,rz) of the specified body, at the specified frame
    public Vector3 r(int bodyID, int frameNum)
    {
        Vector3 pos = new Vector3(q[3 * (bodyID - 1) + 1, frameNum],
                                  q[3 * (bodyID - 1) + 2, frameNum],
                                  q[3 * (bodyID - 1) + 3, frameNum]);
        return pos;

    }


    //returns the position euler parameter vector of the specified body, at the specified frame
    //inputs: 
    //       nb - number of bodies in the system
    //       bodyID - IDnumber of requested body
    //       frameNum - frame number to write
    public Vector4 p(int nb, int bodyID, int frameNum)
    {
        //vector4 = [e1 e2 e3 e0] 
        int rDOFs = 3 * nb;
        Vector4 ep = new Vector4(q[rDOFs + 4 * (bodyID - 1) + 2, frameNum],
                                 q[rDOFs + 4 * (bodyID - 1) + 3, frameNum],
                                 q[rDOFs + 4 * (bodyID - 1) + 4, frameNum],
                                 q[rDOFs + 4 * (bodyID - 1) + 1, frameNum]);
        return ep;
    }

}
