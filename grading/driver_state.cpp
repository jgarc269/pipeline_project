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
    const int triangle_points = 3;
    const data_geometry *dg[3];
    data_geometry temp[3];
    data_vertex in[3];

    switch(type) 
    {
        case render_type::triangle: 
	{
            unsigned int k = 0;
            for (int i = 0; i < state.num_vertices / triangle_points; i++)
	    {
                for (int j = 0; j < triangle_points; j++, k = k + state.floats_per_vertex) 
		{
                    in[j].data = &state.vertex_data[k];
                    temp[j].data = in[j].data;
                    state.vertex_shader(in[j], temp[j], state.uniform_data);
                    dg[j] = &temp[j];
                }
               
                    rasterize_triangle(state, dg);
            }
	}
            break;

 
   	    case render_type::indexed:
	    break;

	    case render_type::fan:
	    break;
        
	    case render_type::strip:
	    break;

	default:
	    break;
    }  
}
void clip_triangle(driver_state& state, const data_geometry* in[3], int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
 //   std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}


void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    float a[3];
    float b[3];
    float c[3];

    float width_div;
    float height_div;

    width_div = state.image_width / 2.0f;
    height_div = state.image_height / 2.0f;

    for(int temp = 0; temp < 3; temp++)
    {
        float i = (width_div * ((*in)[temp].gl_Position[0] / (*in)[temp].gl_Position[3]) + (width_div - 0.5f));
        float j = (height_div * ((*in)[temp].gl_Position[1] / (*in)[temp].gl_Position[3]) + (height_div - 0.5f));
        float k = (width_div * ((*in)[temp].gl_Position[2] / (*in)[temp].gl_Position[3]) + (width_div - 0.5f));
       
	a[temp] = i;
        b[temp] = j;
        c[temp] = k;
    } 

    float min_a = std::min(std::min(a[0], a[1]), a[2]);
    float min_b = std::min(std::min(b[0], b[1]), b[2]);
 
    float max_a = std::max(std::max(a[0], a[1]), a[2]);
    float max_b = std::max(std::max(b[0], b[1]), b[2]);

    if(min_a < 0)
    {
	min_a = 0;
    }
    if(min_b < 0)
    {
	min_b = 0;
    }
    if(max_a > state.image_width)
    {
	max_a = state.image_width;
    }
    if(max_b > state.image_height)
    {
	max_b = state.image_height;
    }

    float area = (0.5f * ((a[1] * b[2] - a[2] * b[1]) - (a[0] * b[2] - a[2] * b[0]) + (a[0] * b[1] - a[1] * b[0])));

    for(int j = min_b; j < max_b; j++)
    {
	for(int i = min_a; i < max_a; i++)
	{
	    float alpha = (0.5f * ((a[1] * b[2] - a[2] * b[1]) + (a[1] - b[2]) * i + (a[2] - a[1]) * j)) / area;
            float beta =  (0.5f * ((a[2] * b[0] - a[0] * b[2]) + (a[2] - b[0]) * i + (a[0] - a[2]) * j)) / area;
            float gamma = (0.5f * ((a[0] * b[1] - a[1] * b[0]) + (a[0] - b[1]) * i + (a[1] - a[0]) * j)) / area;
	   
	   if(alpha >= 0 && beta >= 0 && gamma >= 0) 
	   {
	   	float alpha2 = alpha;
                float beta2 = beta;
                float gamma2 = gamma;
                float z = alpha2 * c[0] + beta2 * c[1] + gamma2 * c[2];

                if(z < state.image_depth[i + j * state.image_width])
		{
                    state.image_depth[i + j * state.image_width] = z;
		    state.image_color[i + j * state.image_width] = make_pixel(255,255,255);
		}
	   }
	}
    }

}
