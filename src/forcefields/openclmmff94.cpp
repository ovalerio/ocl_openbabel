#include <iostream>
#include <cstring>

#include "openclmmff94.h"
#include "util.h"

OpenCLMMFF94::OpenCLMMFF94()
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

    //use the first available device
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
    full_path[0] = '\0'; //strcat searches for '\0' to cat after
    strcat( full_path, CL_SOURCE_DIR ); // Copy base path into full path
    strcat( full_path, "/"); // Assume a backslash is used to separate path names
    strcat( full_path, relative_path); // Copy the relative path to the full_path
    std::cout << "path: " << full_path << std::endl;

    //file_contents is defined in util.cpp
    //it loads the contents of the file at the given path
    char* cSourceCL = file_contents(full_path, &pl);
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


// Using the approach from the book
cl::Program createProgram(cl::Context &context, std::string source_code, const std::string params= "") {

	cl_int r;
	cl::Program::Binaries binaries;
	cl::Program::Sources sources;
	cl::Program program;
	std::vector <cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES> ();

	try {

		binaries = loadBinaries(context, params);
		program = cl::Program(context, devices, binaries);
		try {
			program.build(devices, NULL); // program build no params
		} catch (cl::Error e) {
			std::cout << "Compilation error (\"" << params < "\" log:" <<
			std::endl << program.getBuildInfo<CL_PROGRAM_BUILD_LOG> (devices [0]) <<
			std::endl;
			throw e;
		}
	} catch (cl::Error &e) {
		sources.insert(sources.end(),
				std::make_pair(source_code.c_str(), source_code.length()));

		program = cl::Program(context, sources, &r);
		try {
			program.build(devices, params.c_str()); // program build with params
		} catch (cl::Error e) {
			std::cout << "Compilation error (\"" << params < "\" log:" <<
			std::endl << program.getBuildInfo<CL_PROGRAM_BUILD_LOG> (devices [0]) <<
			std::endl;
			throw e;
		}
		saveBinaries(context, program, params);
	}
	return program;

}//createProgram()


cl::Program::Binaries loadBinaries(cl::context &context, std::string params = "") {
	cl::Program::Binaries binaries;
	std::vector <cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES> ();
	std::string bin_file_name;
	std::vector<cl::Device>::iterator d;
	char *bin_data = NULL;

	for (d = devices.begin(); d != devices.end(); d++) {
		bin_file_name = generateFileName(*d, params);
		std::ifstream bin_file(bin_file_name.c_str(), std::ios::in | std::ios::binary);
		if(bin_file.is_open()){
			size_t bin_size;
			char new_line;
			bin_file.read((char*)&bin_size, sizeof(bin_size));
			bin_file.read((char*)&new_line, sizeof(new_line));
			bin_data = new char[bin_size];
			bin_file.read(bin_data, bin_size);
			binaries.insert(binaries.begin(), std::make_pair(bin_data, bin_size));
			bin_file.close();
		} else {
			throw cl::Error(-1001, "binariesNotAvailable" );
			break;
		}
	}
	return binaries;
}//loadBinaries()


void saveBinaries(cl::Context &context, cl::Program &program, std::string params = "") {

	std::string bin_file_name;
	std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES> ();
	std::vector<size_t> bin_sizes;
	std::vector<char *> bin_binaries;
	std::vector<cl_device_id> bin_devices;

	bin_devices = program.getInfo<CL_PROGRAM_DEVICES> ();
	bin_sizes = program.getInfo<CL_PROGRAM_BINARY_SIZES> ();

	for (::size_t i=0; i < bin_sizes.size(); i++) {
		bin_binaries.insert(bin_binaries.end(), new char [bin_sizes [i]]);
	}

	program.getInfo(CL_PROGRAM_BINARIES, &bin_binaries);

	for (::size_t i = 0; i < devices.size(); i++) {
		bin_file_name = generateFileName(devices[i], params);
		std::ofstream bin_file(bin_file_name.c_str(), std::ios:out | std::ios::binary);

		if (bin_file.is_open()) {
			char new_line = '\n';
			bin_file.write((char*)&bin_sizes [i], sizeof(size_t));
			bin_file.write((char*)&new_line, sizeof(new_line));
			bin_file.write(bin_binaries[i], bin_sizes[i]);
			bin_file.close();
		} else {
			std::cout << "Error writing binary program." << std::endl;
			exit(-2);
		}
	}

	for(::size_t i=0; i < bin_sizes.size(); i++) {
		delete bin_binaries[i];
	}

}

std::string generateFileName(cl::Device &device, std::string params) {
	std::string bin_file_name;
	std::string vendor_name, device_name, app_ver = __TIME__ " " __DATE__;
	device_name = device.getInfo<CL_DEVICE_NAME> ();
	vendor_name = device.getInfo<CL_DEVICE_VENDOR> ();
	bin_file_name = "compiled_" + app_ver + " " + vendor_name + " " + device_name
			+ params + ".bin";

	for(::size_t j=0; j < bin_file_name.length(); j++) {
		if (!(((bin_file_name[j] >= 'a') && (bin_file_name[j] <= 'z')) ||
			  ((bin_file_name[j] >= 'A') && (bin_file_name[j] <= 'Z')) ||
			  ((bin_file_name[j] >= '0') && (bin_file_name[j] <= '9')) ||
			  (bin_file_name[j] == '.') || (bin_file_name[j] == '_'))
			) {
			bin_file_name[j] = '_';
		}
	}
	return bin_file_name;
}

