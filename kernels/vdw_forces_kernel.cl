__kernel void vdwforces(__global float* ax, __global float* ay, __global float* az,
						__global float* bx, __global float* by, __global float* bz,
						__global float* RAB,__global float* epsilon,__global float* energy)
{
    	unsigned int indx = get_global_id(0);
		
		float rab, rab7, RAB7;
		float abx, aby, abz;
	    float erep, erep7, eattr;

		abx = ax[indx] - bx[indx];
		aby = ay[indx] - by[indx];
		abz = az[indx] - bz[indx];
		rab = sqrt(abx * abx + aby * aby + abz * abz);
    	rab7 = pow (rab, 7);
    	RAB7 = pow(RAB[indx], 7);
    	erep = (1.07 * RAB[indx]) / (rab + 0.07 * RAB[indx]);
    	erep7 = pow (erep, 7);
    	eattr = (((1.12 * RAB7) / (rab7 + 0.12 * RAB7)) - 2.0);
    	energy[indx] = epsilon[indx] * erep7 * eattr;
}