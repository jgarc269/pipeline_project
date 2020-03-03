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
    unsigned pixel;
    
    // dividing both height and width by 2
    float width = state.image_width / 2.0f;
    float height = state.image_height / 2.0f;
    
    // Calculate the pixel coordinates of the resulting ​data_geomtry ​position.
    // In particular, the data_geomtry x​ and y positions should be in Normalized Device Coordinates (NDC)
    // with each dimension going from -1 to +1. You will need to transform x from 0 to width and y from 0 to height. 
    // You will also need to account for the fact that the NDC (-1, -1) corresponds to the bottom left corner of the screen
    //  but not the center of the bottom left pixel. Given (x, y) in NDC, what equation gives you (i, j) in pixel space 
    //  (use w to denote width and h to denote height).
    for(int temp = 0; temp < 3; temp++)
    {
        a[temp] = (width * ((*in)[temp].gl_Position[0] / (*in)[temp].gl_Position[3]) + (width - 0.5f));
        b[temp] = (height * ((*in)[temp].gl_Position[1] / (*in)[temp].gl_Position[3]) + (height - 0.5f));
        c[temp] = (width * ((*in)[temp].gl_Position[2] / (*in)[temp].gl_Position[3]) + (width - 0.5f));
    }

    // min and max
    float min_a = std::min(std::min(a[0], a[1]), a[2]);
    float min_b = std::min(std::min(b[0], b[1]), b[2]);

    float max_a = std::max(std::max(a[0], a[1]), a[2]);
    float max_b = std::max(std::max(b[0], b[1]), b[2]);

    // checking to see if tri goes out of bounds
    if(min_a < 0)
        min_a = 0;
    if(min_b < 0)
        min_b = 0;
    if(max_a > state.image_width)
        max_a = state.image_width;
    if(max_b > state.image_height)
        max_b = state.image_height;

    // You can calculate the area of the triangle using the formula:
    // AREA(abc) = 0.5 * ((bxcy − cxby ) − (axcy − cxay) + (axby − bxay))
    float AREA_abc = (0.5f * ((a[1]*b[2] - a[2]*b[1]) - (a[0]*b[2] - a[2]*b[0]) + (a[0]*b[1] - a[1]*b[0])));
 
    // Use the fragment shader to calculate the pixel color rather than setting to white. 
    // See ​data_output in common.h a​ nd the ​fragment_shader​ function in ​driver_state.h.
    auto *data = new float[state.floats_per_vertex];
    data_fragment frag{data};
    data_output output; 

    // To rasterize the triangle, you can iterate over all pixels of the image. 
    // Say you are in the pixel with indices (i, j). You can use the barycentric coordinates of this pixel (i, j)
    // to know if this pixel falls inside the triangle or not.                                                                                                              
     for(int j = min_b; j < max_b; j++)
     {
        for(int i = min_a; i < max_a; i++) 
	{
	    // Barycentric coordinates can be calculated using triangle areas.
            float alpha = (0.5f * ((a[1] * b[2] - a[2] * b[1]) + (b[1] - b[2]) * i + (a[2] - a[1]) * j)) / AREA_abc;
            float beta =  (0.5f * ((a[2] * b[0] - a[0] * b[2]) + (b[2] - b[0]) * i + (a[0] - a[2]) * j)) / AREA_abc;
            float gamma = (0.5f * ((a[0] * b[1] - a[1] * b[0]) + (b[0] - b[1]) * i + (a[1] - a[0]) * j)) / AREA_abc;
            
	    if (alpha >= 0 && beta >= 0 && gamma >= 0)
	    {
                float w = alpha;
                float x = beta;
                float y = gamma;
                float z = w * c[0] + x * c[1] + y * c[2];
		
		pixel = i + j * state.image_width;
                if(z < state.image_depth[pixel])
	    	{
                    state.image_depth[pixel] = z;     
                  
		   // Implement color interpolation by checking ​interp_rules in ​driver_state
		   // before sending the color to the fragment_shader. 
		   // You have one interp_rule for each float in the ​data_geometry.data​. 
		   // If the rule type is noperspective (​ see interpolation types in ​common.h)​ , 
		   // then interpolate the float from the 3 vertices using the barycentric coordinates.
                   for (int k = 0; k < state.floats_per_vertex; k++)
		   {
                        switch (state.interp_rules[k])
			{

			    float temp; // Did not let me delcare this in smooth for some reason
					// so moved it outside of the everything
			    case interp_type::invalid:
			    break;

                            case interp_type::flat: 
                                frag.data[k] = (*in)[0].data[k];
                                break;

                            case interp_type::smooth:                               
                                temp = (w / (*in)[0].gl_Position[3]) + (x / (*in)[1].gl_Position[3]) + (y / (*in)[2].gl_Position[3]);

                                alpha = w / (temp * (*in)[0].gl_Position[3]);
                                beta = x / (temp * (*in)[1].gl_Position[3]);
                                gamma = y / (temp * (*in)[2].gl_Position[3]);

                            case interp_type::noperspective:
                                frag.data[k] = alpha * (*in)[0].data[k] + beta * (*in)[1].data[k] + gamma * (*in)[2].data[k];
                                break;

                            default:
                                break;
                        }
                    }

			state.fragment_shader(frag, output, state.uniform_data);

                        state.image_color[pixel] = make_pixel(output.output_color[0] * 255, output.output_color[1] * 255, output.output_color[2] * 255);
			
                   // state.image_color[pixel] = make_pixel(255,255,255);
                }
            }
        }
   }

	delete [] data;

}
