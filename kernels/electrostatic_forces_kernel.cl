__kernel void forces(__global float* a, __global float* b, __global float* c, __global float* d)
{
    unsigned int indx = get_global_id(0);
    unsigned int i = 3*indx + 0;
    unsigned int j = 3*indx + 1;
    unsigned int k = 3*indx + 2;
    
    d[indx] = c[indx] / ( sqrt(pow(a[i] - b[i], 2) + pow(a[j] - b[j], 2) + pow(a[k] - b[k], 2)) + 0.05 ); 
}