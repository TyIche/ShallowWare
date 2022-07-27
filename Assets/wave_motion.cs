using UnityEngine;
using System.Collections;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.004f;
	float damping 	= 0.996f;
	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag=true;

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;


	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		Vector3[] X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}
	}

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}
	int max(int a,int b)
	{
		return a > b ?a:b;
	}
	int min(int a,int b)
	{
		return a < b ?a:b;
	}
	void Shallow_Wave(float[,] old_h, float[,] h, float [,] new_h)
	{		
		//Step 1:
		//TODO: Compute new_h based on the shallow wave model.
		
		for(int i = 0;i < size;i++) for(int j = 0;j < size;j++)
		{
			int [] dx = {0,0,-1,1};int [] dy = {1,-1,0,0};
			int num = 0;float tmp = 0.0f;
			for(int o = 0; o < 4;o++)
			{
				int xx = i + dx[o];int yy = j + dy[o];
				if(xx < 0 || xx >= size|| yy < 0 || yy >= size) continue;	
				tmp += h[xx,yy];num += 1;
			}
			new_h[i,j] = h[i,j] + (h[i,j] - old_h[i,j]) 
			* damping + (tmp - num * h[i,j])*rate;
		}
		//Step 2: Block->Water coupling
		//TODO: for block 1, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).

		//TODO: for block 2, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).
		cube_v = GameObject.Find("Cube").GetComponent<Transform>().position;
		cube_w = GameObject.Find("Block").GetComponent<Transform>().position;
		int maxx = 0,minx = size-1;int maxy = 0,miny = size-1;
		for(int i = 0;i < size;i++) for(int j = 0;j < size;j++)
		{
			RaycastHit hit;
			float xx = i*0.1f-size*0.05f;
			float zz = j*0.1f-size*0.05f;
			
			if((xx -cube_v.x)*(xx -cube_v.x)+(zz - cube_v.z) * (zz - cube_v.z) > 0.5f * 0.5f * 3&&(xx -cube_w.x)*(xx -cube_w.x)+(zz - cube_w.z) * (zz - cube_w.z) > 0.5f * 0.5f * 3) continue;
			if(Physics.Raycast(new Vector3(xx,-100,zz),
			new Vector3(0,1,0),out hit))
			{
				if(hit.point.y > new_h[i,j]) 
				{
					cg_mask[i,j] = false;b[i,j] = 0;low_h[i,j] = h[i,j];
				}
				else
				{				
					maxx = max(maxx,i);maxy = max(maxy,j);
					minx = min(minx,i);miny = min(miny,j);
					cg_mask[i,j] = true;
					b[i,j] = (new_h[i,j] - hit.point.y)/rate;
					// low_h[i,j] = hit.point.y;
					// Debug.Log("quq" + hit.point);
				}
			}
			else {cg_mask[i,j] = false;b[i,j] = 0;low_h[i,j] = h[i,j];}
		}
		Conjugate_Gradient(cg_mask,b,vh,minx,maxx,miny,maxy);
		// Conjugate_Gradient(cg_mask,b,vh,0,size,0,size);
		// Debug.Log(minx);
		// Debug.Log(maxx);
		// Debug.Log(miny);
		// Debug.Log(maxy);
		// TODO: Diminish vh.

		GameObject cube1 = GameObject.Find("Cube");
		GameObject cube2 = GameObject.Find("Block");
		for(int i = 0;i < size;i++) for(int j = 0;j < size;j++)
		{
			// if(!cg_mask[i,j]) continue;
			vh[i,j] = gamma*vh[i,j];

			float xx = i*0.1f-size*0.05f;
			float zz = j*0.1f-size*0.05f;
			if((xx -cube_v.x)*(xx -cube_v.x)+(zz - cube_v.z) * (zz - cube_v.z) <= 0.5f * 0.5f * 3)
			{
				cube1.GetComponent<myRigid>().takeImpulse(
					new Vector3(0,vh[i,j],0),new Vector3(xx,low_h[i,j],zz));
				
			}
			if((xx -cube_w.x)*(xx -cube_w.x)+(zz - cube_w.z) * (zz - cube_w.z) <= 0.5f * 0.5f * 3)
			{
				cube2.GetComponent<myRigid>().takeImpulse(
					new Vector3(0,vh[i,j],0),new Vector3(xx,low_h[i,j],zz));
				// Debug.Log("L"+vh[i,j]);	
			}
		}
		//TODO: Update new_h by vh.
		for(int i = 0;i < size;i++) for(int j = 0;j < size;j++)
		{
			// if(!cg_mask[i,j]) continue;
			int [] dx = {0,0,-1,1};int [] dy = {1,-1,0,0};
			int num = 0;float tmp = 0.0f;
			for(int o = 0; o < 4;o++)
			{
				int xx = i + dx[o];int yy = j + dy[o];
				if(xx < 0 || xx >= size|| yy < 0 || yy >= size) continue;	
				tmp += vh[xx,yy];num += 1;
			}
			new_h[i,j] = new_h[i,j] + (tmp - num * vh[i,j])*rate;
		}

		//Step 3
		//TODO: old_h <- h; h <- new_h;
		for(int i = 0;i < size;i++) for(int j = 0;j < size;j++)
		{
			old_h[i,j] = h[i,j];
			h[i,j] = new_h[i,j];
		}
		//Step 4: Water->Block coupling.
		//More TODO here.
	}
	

	// Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

		//TODO: Load X.y into h.
		for(int i = 0;i < size;i++) for(int j = 0;j < size;j++) 
		h[i,j] = X[i*size + j].y;
		if (Input.GetKeyDown ("r")) 
		{
			//TODO: Add random water.
			int tmpx = Random.Range(0,size);
			int tmpy = Random.Range(0,size);

			float tmph = Random.Range(0.0f,0.1f);
			h[tmpx,tmpy] += tmph;

			if(tmpy > 0) h[tmpx,tmpy-1] -= tmph;
			else h[tmpx,tmpy+1] -= tmph;
		}
	
		for(int l=0; l<8; l++)
		{
			Shallow_Wave(old_h, h, new_h);
		}

		
		//TODO: Store h back into X.y and recalculate normal.
		for(int i = 0;i < size;i++) for(int j = 0;j < size;j++) 
		X[i*size + j].y = h[i,j];
		
		mesh.vertices = X;
		mesh.RecalculateNormals ();
	}
}
