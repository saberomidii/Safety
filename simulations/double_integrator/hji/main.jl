# Worst case for Double-Integrator system 
# Saber Omidi 
using Dates 
using DelimitedFiles 
using FileIO
using Base:mkpath


using Distributions 
using StatsFuns
using Statistics 
using Random 
using NearestNeighbors 
using JuMP

using SparseArrays 
using LinearAlgebra 

println("Packages are imported.")

# Random Variable 

sigma = 1.0
mean = 0
no_samples = 100


D_list = spzeros(no_samples,1)
for index in 1:no_samples
	D_list[index]= rand(Normal(mean,sigma))
end 

println("Disturbance List is created.")


# Signed Distance Function 

function s_d(x::Float64, y::Float64)

	# points x and y coordinate
	box_p_x=[0.0,4.0]  
	box_p_y=[-3.0,3.0]

	dx = 0.0 
	dy = 0.0	

	if (minimum(box_p_x) <= x <= maximum(box_p_x)) &&
	   (minimum(box_p_y) <= y <= maximum(box_p_y))
	  	
#	   println("x and y  are  inside the box") 	
	  
	   dx=minimum([abs(x-box_p_x[1]),abs(x-box_p_x[2])])
	   dy=minimum([abs(y-box_p_y[1]),abs(y-box_p_y[2])])
	
	   m_dis= [dx,dy]
	   return -norm(m_dis)
        
        elseif (x > maximum(box_p_x) || x < minimum(box_p_x)) && 
	       (minimum(box_p_y) <= y <=  maximum(box_p_y))
	  
 #              println("x is out of the box and y is inside the box")
 
	   
	       if x > maximum(box_p_x)
                 dx = abs(x-maximum(box_p_x))
	       end

	       if x < minimum(box_p_x)
                 dx = abs(x-minimum(box_p_x))
	       end

	       dy = minimum([abs(y-box_p_y[1]),abs(y-box_p_y[2])])
           
	       return norm([dx,dy])
     
        elseif (y > maximum(box_p_y)|| y < minimum(box_p_y)) && 
	       (minimum(box_p_x) <= x <= maximum(box_p_y)) 

  #             println("y is out of the box and x is inside the box")
	    
	       if y > maximum(box_p_y)
	      dy = abs(y-maximum(box_p_y))
	    end

	    if y < minimum(box_p_y)
              dy = abs(y-minimum(box_p_y))
	    end

	    dx = minimum([abs(x-box_p_x[1]),abs(x-box_p_x[2])])

	    return norm([dx,dy])

       else 
   #       	  println("x and y are out the box")
       		   if x > maximum(box_p_x)
	   	      dx = abs(x-maximum(box_p_x))
		   end

		   if x < minimum(box_p_x)
	   	      dx = abs(x-minimum(box_p_x))
	           end

	           if y > maximum(box_p_y)
	  	      dy = abs(y-maximum(box_p_y))
		   end

		   if y < minimum(box_p_y)
	              dy = abs(y-minimum(box_p_y))
	           end
                  
		   return norm([dx,dy]) 
       end
end

#println("testing signed-distance function for four different cases",
#	"[0.0,0.0]","[-5.0,2.0]","[1.0,-4.0]","[5.0,5.0]")
#println("singed distance test:",s_d(0.0,0.0))
#println("signed distance test",s_d(-5.0,2.0))
#println("signed distance test",s_d(1.0,-4.0))
#println("signed distance test:",s_d(5.0,5.0))


# Double Integrator Dynamic 
function di_dynamic(x1::Float64,x2::Float64, u::Float64, d::Float64)
	
	 dt= 0.1
	 x_next = x1 + x2*dt
	 v_next = x2 + (u+d)*dt

	 return [x_next,v_next]
end



# Grids
x1_min = -1.0
x1_max =  5.0

x2_min = -5.0
x2_max =  5.0

u_min  = -2.0
u_max  =  2.0

d_min = -Inf 
d_max =  Inf 

num_points_actions = 11

num_points_state_1 = 61
num_points_state_2 = 101


x1 = collect(LinRange(x1_min, x1_max, num_points_state_1))
x2 = collect(LinRange(x2_min, x2_max, num_points_state_2))

state_2d = [(x,v) for x in x1 for v in x2]

nstates = length(state_2d)

actions = collect(LinRange(u_min,u_max, num_points_actions))

nactions = length(actions)

println("Number of states =$nstates")
println("Number of actions=$nactions")


# Value iteration dp
max_iteration= 1000

L=norm([5.0,5.0])


for state_index in 1:nstates

	s = state_2d[state_index]
        
	x= s[1]
	v= s[2]

	l=-Inf

	for action_index in actions 
		for dis_index in D_list
		next_states = di_dynamic(s[1],s[2],action_index,dis_index)
		signed_distance=s_d(next_states[1],next_states[2])
		
			if l > signed_distance 
				l = min(max(signed_distance,-L),L)
			end
		println("Value of l:",l)
		end
	end

end 


