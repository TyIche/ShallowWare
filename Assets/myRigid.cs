using UnityEngine;
using System.Collections;

public class myRigid : MonoBehaviour 
{
	bool launched 		= true;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	readonly Vector3 g 	= new Vector3(0,-9.8f,0);
	float mass;									// mass
	readonly float muT = 0.5f;
	readonly float muN = 0.5f;
	readonly float muNerf = 0.1f;
	readonly float EPSILON = 0.001f;
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.5f;				// for collision

    public void catched()
    {
        v = Vector3.zero;
        w = Vector3.zero;
        launched = false;
 
    }
    public void released()
    {
        v = Vector3.zero;
        w = Vector3.zero;
        launched = true;
    }
	// Use this for initialization
	void Start () 
	{		
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		Vector3 quq = Vector3.zero;
		for(int i = 0;i < vertices.Length;i++) quq+=vertices[i];
		// Debug.Log(quq);

		float m=80;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
	}
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}
	Vector3 add(Vector3 a,Vector3 b)
	{
		a[0] += b[0];a[1] += b[1];a[2] += b[2];
		return a;
	}
	Matrix4x4 add(Matrix4x4 a,Matrix4x4 b)
	{
		for(int i = 0 ;i < 4;i ++) for(int j = 0;j < 4;j++)
		a[i,j] += b[i,j];
		return a;
	}
	float max(float a,float b)
	{
		return (a > b)?a:b;
	}
	Matrix4x4 mul(float a,Matrix4x4 b)
	{
		for(int i = 0 ; i < 4;i ++) for(int j = 0; j< 4;j ++)
		b[i,j] *= a;
		return b;
	}
	// In this function, update v and w by the impulse due to the collision with
	//a plane <P, N>
	void Collision_Impulse(Vector3 P, Vector3 N)
	{
		// P = P + EPSILON * N;
		N.Normalize();
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] x = mesh.vertices;
		Vector3 c = transform.position;

		Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);
		Matrix4x4 I = R*I_ref*R.inverse;
		Vector3 avg = new Vector3(0,0,0);
		int cnt = 0;
		for(int i = 0;i < x.Length;i++)
		{
			Vector3 xnow1 = add(c , R * x[i]);
			Vector3 vnow1 = add(v , Vector3.Cross(w,R * x[i]));
			if(Vector3.Dot(add(xnow1,-1*P),N)>0 || Vector3.Dot(vnow1,N)>0) continue;
			cnt ++;
			avg = add(avg,x[i]);
		}
		
		if(cnt == 0) return;
		avg[0] /= cnt;avg[1] /= cnt;avg[2] /= cnt;
		Vector3 xnow = add(c , R * avg);
		Vector3 vnow = add(v , Vector3.Cross(w,R * avg));
		Vector3 vN = Vector3.Dot(vnow,N) * N;
		Vector3 vT = add(vnow , -1 * vN);
		float  a;Vector3 vnew;
		
		if(vN.magnitude <= 1) 
		{
			a = max(0,1-muNerf*(1+muNerf*(vN.magnitude/1))*vN.magnitude/vT.magnitude);
			vnew = add(-1*muNerf*(vN.magnitude/1)*vN , a * vT);
		}
		else
		{
			a = max(0,1-muT*(1+muN)*vN.magnitude/vT.magnitude);
			vnew = add(-1*muN*vN , a * vT);
		}

		  

		Matrix4x4 E = new Matrix4x4();
		E[0,0] = E[1,1] = E[2,2] = E[3,3] = 1;
		Matrix4x4 tmp = Get_Cross_Matrix(R * avg);
		Matrix4x4 K = add(mul(1/mass , E) , mul(-1.0f , tmp*I.inverse*tmp));
		Vector3 j = K.inverse * (add(vnew ,-1 * vnow));
		v = add(v , j / mass);
		w = add(w , I.inverse*(Vector3.Cross(R*avg,j)));
		// Debug.Log(v.magnitude);
		// Debug.Log("position");
		// Debug.Log(Vector3.Dot(xnow-P,N)+2*EPSILON);
		// Debug.Log("velocity");
		// Debug.Log(v.magnitude - 10*EPSILON);
		if(Vector3.Dot(xnow-P,N) >= -2*EPSILON && v.magnitude <= 10*EPSILON)
		{
			v = new Vector3(0,0,0);
			w = v;
		} 
		// Debug.Log(v);
		// Debug.Log(w);
	}
	Quaternion add(Quaternion a,Quaternion b)
	{
		a[0] += b[0];a[1] += b[1];a[2] += b[2];a[3] += b[3];
		return a;
	}

    public void takeImpulse(Vector3 j,Vector3 pos)
    {
        Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);
		Matrix4x4 I = R*I_ref*R.inverse;
        v = add(v , j / mass);
		w = add(w , I.inverse*(Vector3.Cross((pos - transform.position),j)));
    }
	// Update is called once per frame
	void Update () 
	{
        // Debug.Log(launched);
		//Game Control
		if(Input.GetKey("l"))
		{
			transform.position = new Vector3 (0, 2, 0);
			v = Vector3.zero;
			// w = new Vector3 (0,0,0.1f);
			launched=true;
			// Debug.Log(v);
		}
		if(!launched) return;
		// Part I: Update velocities
		v = v + g * dt;
		
		// Part II: Collision Impulse
		// Collision_Impulse(new Vector3(0, -0.5f, 0), new Vector3(0, 1, 0));
		// Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));
		// Collision_Impulse(new Vector3(-7, 0.01f, 0), new Vector3(1, 1, 0));
		v = linear_decay * v;
		w = angular_decay * w;
		// Part III: Update position & orientation
		//Update linear status
		Vector3 x    = transform.position;
		x = x + v * dt; 	
		

		//Update angular status
		Quaternion q = transform.rotation;
		Quaternion tmp;
		tmp.w = 0;
		tmp.x = w.x * dt / 2;
		tmp.y = w.y * dt / 2;
		tmp.z = w.z * dt / 2;
		q = add (q , tmp * q);
		q.Normalize();
		
		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}