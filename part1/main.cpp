/*
 * Adventures in OpenCL tutorial series
 * Part 1
 *
 * author: Ian Johnson
 * code based on advisor Gordon Erlebacher's work
 * NVIDIA's examples
 * as well as various blogs and resources on the internet
 */
#include <stdio.h>

#include "cll.h"


int main(int argc, char** argv)
{
    printf("Hello, OpenCL\n");
    //initialize our CL object, this sets up the context
    CL example;
    
    //load and build our CL program from the file
    example.loadProgram("part1.cl");

    //initialize the kernel and send data from the CPU to the GPU
    example.popCorn();
    //execute the kernel
    example.runKernel();
}
