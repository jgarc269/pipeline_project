#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
  
    state.image_width = width;
    state.image_height = height;
    state.image_color = NULL; 
    state.image_depth = NULL;

    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
	
    for(int i = 0; i < width * height; i++)
    {
        state.image_depth[i] = std::numeric_limits<float>::max();
        state.image_color[i] = make_pixel(0, 0, 0);
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;
    
    const int triangle_edges = 3;
    const data_geometry *dg[3];
    data_geometry temp[3];
    data_vertex in[3];
    unsigned int k = 0;
const auto first_vertex = state.vertex_data;
    auto is_first_vertex = true;
    switch(type)
    {	
        case render_type::triangle:
	{ 
	    	for(int i = 0; i < state.num_vertices / triangle_edges; i++)
	        {
	            for(int j = 0; j < triangle_edges; k = k + state.floats_per_vertex)
	            {
			    in[j].data = &state.vertex_data[k];
		            temp[j].data = in[j].data;
                    	    state.vertex_shader(in[j], temp[j], state.uniform_data);
                    	    dg[j] = &temp[j]; 
	             }
			
			rasterize_triangle(state, dg);
	        }
			            
		break;
	}	

	    case render_type::indexed:
	    break;

	    case render_type::fan:
	    break;
        
	    case render_type::strip:
	    break;

    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    
    float a[3];
    float b[3];
    float c[3];
    float width_div;
    float height_div;

    width_div = state.image_width / 2.0f;
    height_div = state.image_height / 2.0f;
	
    //Recall we are using homogeneous coordinates, so you will need to divide the position in the data_geometry​ by w.
    //Calculate the pixel coordinates of the resulting ​data_geomtry ​position. In particular, the data_geomtry x​ and y positions
    //should be in Normalized Device Coordinates (NDC) with each dimension going from -1 to +1. You will need to transform x from 0 to width
    //and y from 0 to height. You will also need to account for the fact that the NDC (-1, -1) corresponds to the bottom left corner of the screen
    // but not the center of the bottom left pixel. Given (x, y) in NDC, what equation gives you (i, j) in pixel space 
    // (use w to denote width and h to denote height). 
    vec3 vertices[3];
    for(int temp = 0; temp < 3; temp++)
    {
        float i = (width_div * ((*in)[temp].gl_Position[0]/(*in)[temp].gl_Position[3]) + (width_div - 0.5f));
        float j = (height_div * ((*in)[temp].gl_Position[1]/(*in)[temp].gl_Position[3]) + (height_div - 0.5f));
        float k = (width_div * ((*in)[temp].gl_Position[2]/(*in)[temp].gl_Position[3]) + (width_div - 0.5f));
        vertices[temp][0] = i;
        vertices[temp][1] = j;
        vertices[temp][2] = k;
    } 
	
//    state.image_color[a[0], a[1], a[2]] = make_pixel(255,255,255);
//    state.image_color[b[0], b[1], b[2]] = make_pixel(255,255,255);
//    state.image_color[c[0], c[1], c[2]] = make_pixel(255,255,255);

    // Min and Max of Triangle
    float min_a = std::min(std::min(a[0], a[1]), a[2]);
    float min_b = std::min(std::min(b[0], b[1]), b[2]);
 
    float max_a = std::max(std::max(a[0], a[1]), a[2]);
    float max_b = std::max(std::max(b[0], b[1]), b[2]);

    //check if triangles goes out of bounds and +1 for extra check
    if(min_a < 0)
    {
	min_a = 0;
    }
    if(min_b < 0)
    {
	min_b = 0;
    }
    if(max_a > state.image_width + 1)
    {
	max_a = state.image_width;
    }
    if(max_b > state.image_height + 1)
    {
	max_b = state.image_height;
    }

    //You can calculate the area of the triangle using the formula:
    //AREA(abc) = 0.5 * ((bxcy − cxby ) − (axcy − cxay) + (axby − bxay))
    float area = (0.5f * ((a[1] * b[2] - a[2] * b[1]) - (a[0] * b[2] - a[2] * b[0]) + (a[0] * b[1] - a[1] * b[0])));
   

    //To rasterize the triangle, you can iterate over all pixels of the image. 
    //Say you are in the pixel with indices (i, j). 
    //You can use the barycentric coordinates of this pixel (i, j) to know if this pixel falls inside the triangle or not. 
    //Barycentric coordinates can be calculated using triangle areas. 

/*
    for(int j = min_b; j < max_b; j++)
    {
	for(int i = min_a; i < max_a; i++)
	{
	    float alpha = (0.5f * ((a[1] * b[2] - a[2] * b[1]) + (a[1] - b[2]) * i + (a[2] - a[1]) * j)) / area;
            float beta =  (0.5f * ((a[2] * b[0] - a[0] * b[2]) + (a[2] - b[0]) * i + (a[0] - a[2]) * j)) / area;
            float gamma = (0.5f * ((a[0] * b[1] - a[1] * b[0]) + (a[0] - b[1]) * i + (a[1] - a[0]) * j)) / area;

	    state.image_color[i + j * state.image_width] = make_pixel(255,255,255);
	}
    }
*/
}

