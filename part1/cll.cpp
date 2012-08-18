#include <iostream>
#include <cstring>

#include "cll.h"
#include "util.h"

CL::CL()
{
    std::cout << "Initialize OpenCL object and context" << std::endl;
    //setup devices and context
    
    //this function is defined in util.cpp
    //it comes from the NVIDIA SDK example code
    err = oclGetPlatformID(&platform);
    //oclErrorString is also defined in util.cpp and comes from the NVIDIA SDK
    std::cout << "oclGetPlatformID: " << oclErrorString(err) << std::endl;

    // Get the number of GPU devices available to the platform
    // we should probably expose the device type to the user
    // the other common option is CL_DEVICE_TYPE_CPU
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
    std::cout << "clGetDeviceIDs (get number of devices): " << oclErrorString(err) << std::endl;


    // Create the device list
    devices = new cl_device_id [numDevices];
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
    std::cout << "clGetDeviceIDs (create device list): " << oclErrorString(err) << std::endl;
 
    //create the context
    context = clCreateContext(0, 1, devices, NULL, NULL, &err);

    //for right now we just use the first available device
    //later you may have criteria (such as support for different extensions)
    //that you want to use to select the device
    deviceUsed = 0;
    
    //create the command queue we will use to execute OpenCL commands
    command_queue = clCreateCommandQueue(context, devices[deviceUsed], 0, &err);

    cl_a = 0;
    cl_b = 0;
    cl_c = 0;
}

CL::~CL()
{
    std::cout << "Releasing OpenCL memory" << std::endl;
    if(kernel)clReleaseKernel(kernel); 
    if(program)clReleaseProgram(program);
    if(command_queue)clReleaseCommandQueue(command_queue);
   
    //need to release any other OpenCL memory objects here
    if(cl_a)clReleaseMemObject(cl_a);
    if(cl_b)clReleaseMemObject(cl_b);
    if(cl_c)clReleaseMemObject(cl_c);

    if(context)clReleaseContext(context);
    
    if(devices)delete(devices);
    std::cout << "OpenCL memory released" << std::endl;
}


void CL::loadProgram(const char* relative_path)
{
 // Program Setup
    int pl;
    size_t program_length;
    std::cout << "load the program" << std::endl;
    char full_path[500]; // probably long enough to hold a full path
    
    //CL_SOURCE_DIR is set in the CMakeLists.txt
//    char const* base_path = getenv( CL_SOURCE_DIR );

    full_path[0] = '\0'; //strcat searches for '\0' to cat after
    strcat( full_path, CL_SOURCE_DIR ); // Copy base path into full path
    strcat( full_path, "/"); // Assume a backslash is used to separate path names
    strcat( full_path, relative_path); // Copy the relative path to the full_path
    std::cout << "path: " << full_path << std::endl;

    //file_contents is defined in util.cpp
    //it loads the contents of the file at the given path
    char* cSourceCL = file_contents(full_path, &pl);
    //printf("file: %s\n", cSourceCL);
    program_length = (size_t)pl;

    // create the program
    program = clCreateProgramWithSource(context, 1,
                      (const char **) &cSourceCL, &program_length, &err);
    std::cout << "clCreateProgramWithSource: " << oclErrorString(err) << std::endl;

    buildExecutable();

    //Free buffer returned by file_contents
    free(cSourceCL);
}

//----------------------------------------------------------------------
void CL::buildExecutable()
{
    // Build the program executable
    
    std::cout << "building the program" << std::endl;
    // build the program
    //err = clBuildProgram(program, 0, NULL, "-cl-nv-verbose", NULL, NULL);
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    std::cout << "clBuildProgram: " << oclErrorString(err) << std::endl;
	//if(err != CL_SUCCESS){
		cl_build_status build_status;
		err = clGetProgramBuildInfo(program, devices[deviceUsed], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL);

		char *build_log;
		size_t ret_val_size;
		err = clGetProgramBuildInfo(program, devices[deviceUsed], CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);

		build_log = new char[ret_val_size+1];
		err = clGetProgramBuildInfo(program, devices[deviceUsed], CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);
		build_log[ret_val_size] = '\0';
		std::cout <<"BUILD LOG: " << std::endl << build_log << std::endl;
	//}
    std::cout << "program built" << std::endl;
}

