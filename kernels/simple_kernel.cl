__kernel void forces(__global double* a, __global double* b, __global double* c)
{
    unsigned int indx = get_global_id(0);
    c[indx] = a[indx] + b[indx];
}
