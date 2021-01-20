clear all; clc;
%This script has only been tested to work with MATLAB R2019b
%Custom blocks configuration
config = [  0 1 0; 
            0 2 0; 
            0 1 0;
            0 2 0];
        % 0 = empty space
        % 1 = non-heat generating cube
        % 2 = cube with heat gen
        
%Cube face config
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Front Face is Face 1, turn right is Face 2, ect...
    % Top Plate is Face 5
    % Bottom Plate is Face 6
    %          +--------+
    %         /   5    /|
    %        /        / |
    %       +--------+  |   3 is the back face
    %       |        | 2|
    %  4->  |    1   |  +
    %       |        | /
    %       |        |/
    %       +--------+
    %          ^
    %          |
    %          6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%faces being radiated by the sun        
sun =       [1 0 0 0 0 0];
% 1 = sun is illuminating the face
% 0 = sun is not illuminating face

%faces being radiated by the earth
earth =     [0 0 1 0 0 0];
% 1 = sun is illuminating the face
% 0 = sun is not illuminating face

%to determine if it is 2D or 3D config (wont support 3D config unless manually processed(will compute))
dim = length(size(config));
dx = .1;            %1cm
dy = .1;            %1cm
dz = .03;           %0.3cm
dt = 1;             %1 second
edot = 15000;       %450W/m^2 /0.03m = 15000 W/m^3
delta = 5.67e-8;
T_space = 3;
Q_sun = 1357;   %W/m^2
Q_earth = 469;  %W/m^2
% Q_space = 0.08; %https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-in-space
%not using these
% T_sun = 5778;
% T_earth = 253.7;

%Want to know if var is a int
isaninteger = @(x)isfinite(x) & x==floor(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dim == 1
    disp('funny enough this isnt supported. Add some 0s to make it 2D and it will work');
    NumObj = 1;
elseif dim == 2
    while 1
       prompt = ['How big do you want the cubes? (1 = 1x1x1cm)', newline,'Size selection:   '];
       x = input(prompt);
       if isaninteger(x)
           if x < 3
               disp('Too small!!');
           else
               break
           end
       else
           disp(['Please enter a integer!',newline]);
       end
    end
    %The size of the cubes
    block_size =  x;
    
    %getting config dimension size since it might not be square
    NN = size(config);
    %row limit
    Nr = NN(1);
    %column limit
    Nc = NN(2);
    
    %need to label each cube
    label_pointer = 1;
    labeled_blocks = zeros(Nr,Nc);
    for row = 1:1:Nr
        for col = 1:1:Nc
            if config(row,col) == 1 || config(row,col) == 2
                labeled_blocks(row,col) = label_pointer;
                label_pointer = label_pointer + 1;
            end
        end
    end
    %Getting number of objects in the problem
    NumObj = label_pointer - 1;
    %
    k = zeros(NumObj,1);
    Cp = zeros(NumObj,1);
    rho = zeros(NumObj,1);
    Emiss = zeros(NumObj,1);
    %Different material prompt
    while 1
        prompt = ['Do you want all the cubes to be the same material?',newline,'1: Yes',newline,'2: No',newline,'1 or 2:    '];
        x = input(prompt);
        if isaninteger(x)
           if x == 1
               %Same mat
               while 1
                   prompt = ['1: Copper',newline,'2: Aluminum',newline,'3: Gold',newline,'4: Iron',newline,'5: Silver',newline,'Enter selection:    '];
                   x = input(prompt);
                   if isaninteger(x)
                       if x == 1
                           k(:) = 413;
                           Cp(:) = 376.812;
                           rho(:) = 8960;
                           Emiss(:) = 0.87;
                           break
                       elseif x == 2
                           k(:) = 237;
                           Cp(:) = 921.096;
                           rho(:) = 2710;
                           Emiss(:) = 0.11;
                           break
                       elseif x == 3
                           k(:) = 327;
                           Cp(:) = 125.604;
                           rho(:) = 19300;
                           Emiss(:) = 0.02;
                           break
                       elseif x == 4
                           k(:) = 94;
                           Cp(:) = 460.548;
                           rho(:) = 7874;
                           Emiss(:) = 0.74;
                           break
                       elseif x == 5
                           k(:) = 403;
                           Cp(:) = 238.6476;
                           rho(:) = 10497;
                           Emiss(:) = 0.01;
                           break
                       else
                           disp('Please select an option for the list.');
                       end
                   end
               end
               break
           elseif x == 2
               i = 1;
               while 1
                   clc;
                   disp(['Looking at cube ',num2str(i),', what material do you want to set it as?',newline]);
                   disp(num2str(labeled_blocks));
                   disp(newline);
                   prompt = ['1: Copper',newline,'2: Aluminum',newline,'3: Gold',newline,'4: Iron',newline,'5: Silver',newline,'Enter selection:    '];
                   x = input(prompt);
                   if isaninteger(x)
                       if x == 1
                           k(i) = 413;
                           Cp(i) = 376.812;
                           rho(i) = 8960;
                           Emiss(i) = 0.87;
                           i = i + 1;
                           if i == NumObj+1
                               break
                           end
                       elseif x == 2
                           k(i) = 237;
                           Cp(i) = 921.096;
                           rho(i) = 2710;
                           Emiss(i) = 0.11;
                           i = i + 1;
                           if i == NumObj+1
                               break
                           end
                       elseif x == 3
                           k(i) = 327;
                           Cp(i) = 125.604;
                           rho(i) = 19300;
                           Emiss(i) = 0.02;
                           i = i + 1;
                           if i == NumObj+1
                               break
                           end
                       elseif x == 4
                           k(i) = 94;
                           Cp(i) = 460.548;
                           rho(i) = 7874;
                           Emiss(i) = 0.74;
                           i = i + 1;
                           if i == NumObj+1
                               break
                           end
                       elseif x == 5
                           k(i) = 403;
                           Cp(i) = 238.6476;
                           rho(i) = 10497;
                           Emiss(i) = 0.01;
                           i = i + 1;
                           if i == NumObj+1
                               break
                           end
                       else
                           disp('Please select an option for the list.');
                       end
                   end
                       
               end
               break
           else
               disp(['Please pick from the listed numbers!',newline]);
           end
       else
           disp(['Please enter a integer!',newline]);
       end
            
    end
    %Simulation Length prompt
    while 1
        prompt = ['How long would you like to run this simulation? (1=1 second)',newline];
        x = input(prompt);
        if isaninteger(x)
            disp(['Simulating configuration for: ', num2str(x), ' seconds']);
            cycles = x;
            break
        end
    end
    %what block has heat gen
    edot_obj = zeros(label_pointer-1,1);
    pointer = 1;
    for row = 1:1:Nr
        for col = 1:1:Nc
            if config(row,col) == 1
                pointer = pointer +1;
            elseif config(row,col) == 2
                edot_obj(pointer,1) = 1;
            end
        end
    end
    %Metal Properties
    %Copper - 413 W/m*k                 Cp = 376.812        Emiss=0.87
    %Aluminum - 237 W/m*K               Cp = 921.096        Emiss=0.11
    %Gold - 327 W/m*K                   Cp = 125.604        Emiss=0.02
    %Iron - 94 W/m*k                    Cp = 460.548        Emiss=0.74
    %Silver - 403 W/m*K                 Cp = 238.6476       Emiss=0.01
    %Values from: https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
    %Emiss values from: https://www.omega.co.uk/literature/transactions/volume1/emissivitya.html
    
    %want to figure out what faces are touching
    connected_faces = zeros(label_pointer - 1,6);%keeps track of faces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Front Face is Face 1, turn right is Face 2, ect...
    % Top Plate is Face 5
    % Bottom Plate is Face 6
    %      +--------+
    %     /   5    /|
    %    /        / |
    %   +--------+  |
    %   |        | 2|
    %   |    1   |  +
    %   |        | /
    %   |        |/
    %   +--------+
    %       ^
    %       |
    %       6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %2D interp, generates connected faces
    for row = 1:1:Nr
        for col = 1:1:Nc
            if config(row,col) == 1 || config(row,col) == 2
                %Cases
                %Above: connected_faces(labeled_blocks(row,col),5) =
                %labeled_blocks(row-1,col);
                
                %Below: connected_faces(labeled_blocks(row,col),6) =
                %labeled_blocks(row+1,col);
                
                %Left: connected_faces(labeled_blocks(row,col),4) = 
                %labeled_blocks(row,col-1);
                
                %Right: connected_faces(labeled_blocks(row,col),2) =
                %labeled_blocks(row,col+1);
                
                %
                %
                %top coner
                if(row == 1 && col == 1)
                    %check to the right
                    if config(row,col+1) == 1 || config(row,col+1) == 2
                        connected_faces(labeled_blocks(row,col),2) = labeled_blocks(row,col+1);
                    end
                    %check to under
                    if config(row+1,col) == 1 || config(row+1,col) == 2
                        connected_faces(labeled_blocks(row,col),6) = labeled_blocks(row+1,col);
                    end
                elseif (row == 1 && col > 1 && col <Nc)
                    %Checking to the left
                    if config(row,col-1) == 1 || config(row,col-1) == 2
                        connected_faces(labeled_blocks(row,col),4) = labeled_blocks(row,col-1);
                    end
                    %Checking to the right
                    if config(row,col+1) == 1 || config(row,col+1) == 2
                        connected_faces(labeled_blocks(row,col),2) = labeled_blocks(row,col+1);
                    end
                    %Checking Under
                    if config(row+1,col) == 1 || config(row+1,col) == 2
                        connected_faces(labeled_blocks(row,col),6) = labeled_blocks(row+1,col);
                    end
                elseif (row == 1 && col == Nc)
                    %Code
                    %left
                    if config(row,col-1) == 1 || config(row,col-1) == 2
                        connected_faces(labeled_blocks(row,col),4) = labeled_blocks(row,col-1);
                    end
                    %under
                    if config(row+1,col) == 1 || config(row+1,col) == 2
                        connected_faces(labeled_blocks(row,col),6) = labeled_blocks(row+1,col);
                    end
                elseif (row > 1 && row < Nr && col == 1)
                    %above
                    if config(row-1,col) == 1 || config(row-1,col) == 2
                        connected_faces(labeled_blocks(row,col),5) = labeled_blocks(row-1,col);
                    end
                    %under
                    if config(row+1,col) == 1 || config(row+1,col) == 2
                        connected_faces(labeled_blocks(row,col),6) = labeled_blocks(row+1,col);
                    end
                    %right
                    if config(row,col+1) == 1 || config(row,col+1) == 2
                        connected_faces(labeled_blocks(row,col),2) = labeled_blocks(row,col+1);
                    end
                elseif (row == Nr && col == 1)
                    %Code
                    %above
                    if config(row-1,col) == 1 || config(row-1,col) == 2
                        connected_faces(labeled_blocks(row,col),5) = labeled_blocks(row-1,col);
                    end
                    %right
                    if config(row,col+1) == 1 || config(row,col+1) == 2
                        connected_faces(labeled_blocks(row,col),2) = labeled_blocks(row,col+1);
                    end
                elseif (row > 1 && row < Nr && col == Nc)
                    %up
                    if config(row-1,col) == 1 || config(row-1,col) == 2
                        connected_faces(labeled_blocks(row,col),5) = labeled_blocks(row-1,col);
                    end
                    %left
                    if config(row,col-1) == 1 || config(row,col-1) == 2
                        connected_faces(labeled_blocks(row,col),4) = labeled_blocks(row,col-1);
                    end
                    %below
                    if config(row+1,col) == 1 || config(row+1,col) == 2
                        connected_faces(labeled_blocks(row,col),6) = labeled_blocks(row+1,col);
                    end
                elseif (row == Nr && col > 1 && col < Nc)
                    %left
                    if config(row,col-1) == 1 || config(row,col-1) == 2
                        connected_faces(labeled_blocks(row,col),4) = labeled_blocks(row,col-1);
                    end
                    %right
                    if config(row,col+1) == 1 || config(row,col+1) == 2
                        connected_faces(labeled_blocks(row,col),2) = labeled_blocks(row,col+1);
                    end
                    %above
                    if config(row-1,col) == 1 || config(row-1,col) == 2
                        connected_faces(labeled_blocks(row,col),5) = labeled_blocks(row-1,col);
                    end
                elseif (row == Nr && col == Nc)
                    %Code
                    %above
                    if config(row-1,col) == 1 || config(row-1,col) == 2
                        connected_faces(labeled_blocks(row,col),5) = labeled_blocks(row-1,col);
                    end
                    %left
                    if config(row,col-1) == 1 || config(row,col-1) == 2
                        connected_faces(labeled_blocks(row,col),4) = labeled_blocks(row,col-1);
                    end
                else
                    %above
                    if config(row-1,col) == 1 || config(row-1,col) == 2
                        connected_faces(labeled_blocks(row,col),5) = labeled_blocks(row-1,col);
                    end
                    %below
                    if config(row+1,col) == 1 || config(row+1,col) == 2
                        connected_faces(labeled_blocks(row,col),6) = labeled_blocks(row+1,col);
                    end
                    %left
                    if config(row,col-1) == 1 || config(row,col-1) == 2
                        connected_faces(labeled_blocks(row,col),4) = labeled_blocks(row,col-1);
                    end
                    %right
                    if config(row,col+1) == 1 || config(row,col+1) == 2
                        connected_faces(labeled_blocks(row,col),2) = labeled_blocks(row,col+1);
                    end
                end
                %TODO
                
                
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Numbers needed throughout the program
    N = block_size;
    PointsPerCube = N^2+2*N*(N-1)+N*(N-2)+2*(N-2)^2;
    A = 2*(N+2*(N-1)+(N-2));
    %NumObj is defined from dim parsing
    %Connected_faces is defined from parsing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [f1,f2,f3,f4,f5,f6] = template_faces(N);
    bors = generate_L_neighbors(N,PointsPerCube,f1,f2,f3,f4,f5,f6);
    bors_T = generate_neighbors(bors,PointsPerCube,NumObj,connected_faces,f1,f2,f3,f4,f5,f6,N);
    %function [outputs] = function_name(inputs)
    num_bors = generate_NL_neighbors(PointsPerCube,NumObj,bors_T);
elseif dim == 3
    %Three Dim interp
    disp('not supported yet');
else
    disp('Anything larger than 3D is a bit much');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%All Non-local neighbors accounted for
errors = 0;

T = zeros(cycles+1,PointsPerCube*NumObj);
%Setting inital temperature at 0C or 273K
T(1,:) = 273;
%Calculating temperatures
%This loop could have been about 85% smaller but I didnt write functions
for j = 1:cycles
    for i = 1:PointsPerCube*NumObj
        %want to know where I am
        %determining what block and am at and the point on the standard
        %format
        if(mod(i,PointsPerCube) == 0)
            block = i/PointsPerCube;
            standard_point = PointsPerCube;
        else
            block = ceil(i/PointsPerCube);
            standard_point = mod(i,PointsPerCube);
        end
        if(standard_point <= 2*(N+2*(N-1)+(N-2))+4*(N-2))
            %Edge or corner range
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(standard_point == 1)
                %Corner 1
                %Face 1,4,5
                %All formulas have been calculated that only one unit face
                %is exposed
                corners = [1,4,5];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif(standard_point == N)
                %Coner 2: 1,2,5
                corners = [1,2,5];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point == N+(N-1))
                %Coner 3: 2,3,5
                corners = [2,3,5];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point == N+2*(N-1))
                %Coner 4: 3,4,5
                corners = [3,4,5];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point == N+3*(N-1))
                %Coner 5: 1,4,6
                corners = [1,4,6];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point == N+4*(N-1))
                %corner 6: 1,2,6
                corners = [1,2,6];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point == N+5*(N-1))
                %Coner 7: 2,3,6
                corners = [2,3,6];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point == N+6*(N-1))
                %Coner 8: 3,4,6
                corners = [3,4,6];
                sun_faces = sunCheckCorner(sun,connected_faces,block,corners);
                earth_faces = earthCheckCorner(earth,connected_faces,block,corners);
                space_faces = spaceCheckCorner(connected_faces,block,corners);
                T(j+1,i) = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= 2 && standard_point <=N-1)
                %Edge 1
                edges = [1,5];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= N+1 && standard_point <=N+(N-2))
                %Edge 2
                edges = [2,5];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= 2*N && standard_point <=(N-1)*3)
                %Edge 3
                edges = [3,5];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= 2*N+(N-1) && standard_point <=N+2*(N-1)+(N-2))
                %Edge 4
                edges = [4,5];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= 2*N+2*(N-1) && standard_point <=(N-1)*5)
                %Edge 5
                edges = [1,6];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= N+2*(N-1)+(N-2)+(N-1) && standard_point <=2*N+2*(N-1)+2*(N-2))
                %Edge 6
                edges = [2,6];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= 2*N+4*(N-1) && standard_point <=2*N+3*(N-1)+2*(N-2))
                %Edge 7
                edges = [3,6];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= A-(N-3) && standard_point <= A)
                %Edge 8
                edges = [4,6];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= A+1 && standard_point <=A+(N-2))
                %Edgy 9
                edges = [1,4];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= A+(N-2)+1 && standard_point <=A+2*(N-2))
                %Edge 10
                edges = [1,2];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= A+2*(N-2)+1 && standard_point <=A+3*(N-2))
                %Edge 11
                edges = [2,3];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            elseif(standard_point >= A+3*(N-2)+1 && standard_point <=A+4*(N-2))
                %Edge 12
                edges = [3,4];
                sun_faces = sunCheckEdge(sun,connected_faces,block,edges);
                earth_faces = earthCheckEdge(earth,connected_faces,block,edges);
                space_faces = spaceCheckEdge(connected_faces,block,edges);
                T(j+1,i) = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces);
            else
                %In the range of the corners and edges but did fit one of
                %conditions
                disp(['error: Unaccounted Node -',num2str(j)])
            end
        else
            %None edge or corner nodes
            if(standard_point >=A+4*(N-2)+1 && standard_point <=A+4*(N-2)+(N-2)^2)
                %Face 1
                E_bool_sun = sunCheckFace(sun,connected_faces,block,1);
                E_bool_earth = earthCheckFace(earth,connected_faces,block,1);
                E_bool_space = spaceCheckFace(connected_faces,block,1);
                T(j+1,i) = faceTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,E_bool_sun,Q_earth,E_bool_earth,edot,edot_obj,T_space,E_bool_space);
            elseif(standard_point >= A+4*(N-2)+(N-2)^2+1 && standard_point <=A+4*(N-2)+2*(N-2)^2 )
                %Face 2
                E_bool_sun = sunCheckFace(sun,connected_faces,block,2);
                E_bool_earth = earthCheckFace(earth,connected_faces,block,2);
                E_bool_space = spaceCheckFace(connected_faces,block,2);
                T(j+1,i) = faceTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,E_bool_sun,Q_earth,E_bool_earth,edot,edot_obj,T_space,E_bool_space);
            elseif(standard_point >= A+4*(N-2)+2*(N-2)^2+1  && standard_point <= A+4*(N-2)+3*(N-2)^2 )
                %Face 3
                E_bool_sun = sunCheckFace(sun,connected_faces,block,3);
                E_bool_earth = earthCheckFace(earth,connected_faces,block,3);
                E_bool_space = spaceCheckFace(connected_faces,block,3);
                T(j+1,i) = faceTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,E_bool_sun,Q_earth,E_bool_earth,edot,edot_obj,T_space,E_bool_space);
            elseif(standard_point >= A+4*(N-2)+3*(N-2)^2 + 1 && standard_point <= A+4*(N-2)+4*(N-2)^2 )
                %Face 4
                E_bool_sun = sunCheckFace(sun,connected_faces,block,4);
                E_bool_earth = earthCheckFace(earth,connected_faces,block,4);
                E_bool_space = spaceCheckFace(connected_faces,block,4);
                T(j+1,i) = faceTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,E_bool_sun,Q_earth,E_bool_earth,edot,edot_obj,T_space,E_bool_space);
            elseif(standard_point >= A+4*(N-2)+4*(N-2)^2 + 1 && standard_point <= A+4*(N-2)+5*(N-2)^2 )
                %Face 5
                E_bool_sun = sunCheckFace(sun,connected_faces,block,5);
                E_bool_earth = earthCheckFace(earth,connected_faces,block,5);
                E_bool_space = spaceCheckFace(connected_faces,block,5);
                T(j+1,i) = faceTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,E_bool_sun,Q_earth,E_bool_earth,edot,edot_obj,T_space,E_bool_space);
            elseif(standard_point >= A+4*(N-2)+5*(N-2)^2 + 1 && standard_point <= PointsPerCube )
                %Face 6
                E_bool_sun = sunCheckFace(sun,connected_faces,block,6);
                E_bool_earth = earthCheckFace(earth,connected_faces,block,6);
                E_bool_space = spaceCheckFace(connected_faces,block,6);
                T(j+1,i) = faceTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,E_bool_sun,Q_earth,E_bool_earth,edot,edot_obj,T_space,E_bool_space);
            else
                disp('error: Out of Bounds');
            end
        end

        if T(j,i)<0
            errors = errors + 1;
        elseif T(j,i) > 1e5
            errors = errors + 1;
        end 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If any temperaturs were below zero or above 10,000k.

disp([newline,'Errors Detected: ',num2str(errors),newline])

if errors == 0
    disp('Please refer to variable ''T'' to examine the temperatures at each node');
else
    disp('There seems to have been an error! A number in variable T appears to be too large or below zero')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [outputs] = function_name(inputs)
function [cube] = generate(p,N)
    A = zeros(N,N,N);
    
    counter = p;
    
    depth = 1;
    row = 1;
    %top row
    for col = 1:1:N
        A(row,col,depth) = counter;
        counter = counter +1;
    end
    %row = 1
    %col = 10
    for depth = 2:1:N
        A(row,col,depth) = counter;
        counter = counter +1;
    end
    %depth = 10
    for col = N-1:-1:1
        A(row,col,depth) = counter;
        counter = counter +1;
    end
    %col = 1;
    for depth = N-1:-1:2
        A(row,col,depth) = counter;
        counter = counter +1;
    end
    
    %Bottom row
    %depth = 2
    depth = 1;
    row = N;
    for col = 1:1:N
        A(row,col,depth) = counter;
        counter = counter + 1;
    end
    %col = 10
    %row = 10;
    for depth = 2:1:N
        A(row,col,depth) = counter;
        counter = counter +1;
    end
    
    for col = N-1:-1:1
        A(row,col,depth) = counter;
        counter = counter +1;
    end
    %col =1
    for depth = N-1:-1:2
        A(row,col,depth) = counter;
        counter = counter +1;
    end
    depth = 1;
    for row = 2:1:N-1
        A(row,col,depth) = counter;
        counter = counter + 1;
    end
    col = N;
    for row = 2:1:N-1
        A(row,col,depth) = counter;
        counter = counter + 1;
    end
    depth = N;
    for row = 2:1:N-1
        A(row,col,depth) = counter;
        counter = counter + 1;
    end
    col = 1;
    for row = 2:1:N-1
        A(row,col,depth) = counter;
        counter = counter + 1;
    end
    
    %midpoint face 1
    depth = 1;
    for row = 2:1:N-1
        for col = 2:1:N-1
            A(row,col,depth) = counter;
            counter = counter +1 ;
        end
    end
    
    %midpoint face 2
    col = N;
    for row = 2:1:N-1
        for depth = 2:1:N-1
            A(row,col,depth) = counter;
            counter = counter + 1;
        end
    end
    
    %midpoint face 3
    depth = N;
    for row = 2:1:N-1
        for col = N-1:-1:2
            A(row,col,depth) = counter;
            counter = counter + 1;
        end
    end
    
    %midpoints face 4
    
    col = 1;
    for row = 2:1:N-1
        for depth = N-1:-1:2
            A(row,col,depth) = counter;
            counter = counter + 1;
        end
    end
    
    %midpoints face 5
    
    row = 1;
    for depth = N-1:-1:2
        for col = 2:1:N-1
            A(row,col,depth) = counter;
            counter = counter + 1;
        end
    end
    
    row = N;
    for depth = 2:1:N-1
        for col = 2:1:N-1
            A(row,col,depth) = counter;
            counter = counter + 1;
        end
    end       
        cube = A;
end
function[f1,f2,f3,f4,f5,f6] = template_faces(N)
    %Inputs: N
    %Outputs: f1,f2,f3,f4,f5,f6
    %Generating a template cube faces
    bob = generate(1,N);
    %Setting up faces to reference later
    f1 = bob(:,:,1);
    f2 = zeros(N,N);
    for row = 1:1:N
        for depth = 1:1:N
            f2(row,depth) = bob(row,N,depth);
        end
    end
    f3 = zeros(N,N);
    for row = 1:1:N
        for col = N:-1:1
            f3(row,col) = bob(row,col,depth);
        end
    end
    f4 = zeros(N,N);
    i = 1;
    col = 1;
    for row = 1:1:N
        for depth = N:-1:1
            f4(row,i) = bob(row,col,depth);
            i = i +1 ;
        end
        i=1;
    end
    f5 = zeros(N,N);
    row = 1;
    i = 1;
    y = 1;
    for depth = N:-1:1
        for col = 1:1:N
            %Code
            f5(y,i) = bob(row,col,depth);
            i = i + 1;
        end
        i = 1;
        y = y +1;
    end
    f6 = zeros(N,N);
    for row = 1:1:N
        for col = 1:1:N
            %row = 10;
            f6(row,col) = bob(N,col,row);
        end
    end

end
function [bors] = generate_L_neighbors(N,PointsPerCube,f1,f2,f3,f4,f5,f6)
    % Inputs: N, PointsPerCube
    % Outputs: bors
    %Generater Neighbors (What nodes are touching)
    bors = zeros(PointsPerCube,7);
    %code here
    A = 2*(N+2*(N-1)+(N-2));
    %Corner 1
    bors(1,1) = 2;
    bors(1,2) = (N+2*(N-1)+(N-2));
    bors(1,3) = 2*(N+2*(N-1)+(N-2))+ 1;
    %Corner 2
    bors(N,1) = (N-1);
    bors(N,2) = (N+1);
    bors(N,3) = A + (N-2) + 1;
    %Corner 3
    bors(N+(N-1),1) = N + (N-1) + 1;
    bors(N+(N-1),2) = (N+(N-1)-1);
    bors(N+(N-1),3) = A + 2*(N-2)+1;
    %Corner 4
    bors(N+2*(N-1),1) = N + 2*(N-1) - 1; 
    bors(N+2*(N-1),2) = N + 2*(N-1) + 1;
    bors(N+2*(N-1),3) = A + 3*(N-2) + 1;
    %Corner 5
    bors(N+3*(N-1),1) = N+3*(N-1) + 1; 
    bors(N+3*(N-1),2) = A;
    bors(N+3*(N-1),3) = A + (N-2);
    %Corner 6
    bors(N+4*(N-1),1) = N+4*(N-1) - 1;
    bors(N+4*(N-1),2) = N+4*(N-1) + 1;
    bors(N+4*(N-1),3) = A + 2*(N-2);
    %Corner 7
    bors(N+5*(N-1),1) = N+5*(N-1) - 1;
    bors(N+5*(N-1),2) = N+5*(N-1) + 1;
    bors(N+5*(N-1),3) = A + 3*(N-2);
    %Corner 8
    bors(N+6*(N-1),1) = N+6*(N-1) - 1;
    bors(N+6*(N-1),2) = N+6*(N-1) + 1;
    bors(N+6*(N-1),3) = A + 4*(N-2);

    for P = 2:N-1
        bors(P,1) = P - 1; 
        bors(P,2) = P + 1;
        bors(P,3) = (P-1) + A + 4*(N-2);
        bors(P,4) = PointsPerCube - (N-2)^2 - (N-2) + (P-1);
    end
    for P = N+1:N+(N-2)
        bors(P,1) = P - 1;
        bors(P,2) = P + 1;
        bors(P,3) = A + 4*(N-2) + (N-2)^2 + (P - N);
        bors(P,4) = PointsPerCube - (N-2)^2 - (P-(N+1))*(N-2);
    end
    for P = 2*N:(N-1)*3
        bors(P,1) = P - 1;
        bors(P,2) = P + 1;
        bors(P,3) = A + 4*(N-2) + 2*(N-2)^2 + (P-(N+(N-1)));
        bors(P,4) = A + 4*(N-2) + 4*(N-2)^2 + (N-2) + (2*N-P);
    end
    for P = 2*N+(N-1):N+2*(N-1)+(N-2)
        bors(P,1) = P - 1;
        bors(P,2) = P + 1;
        bors(P,3) = A + 4*(N-2) + 3*(N-2)^2 + (P - (2*N+(N-2)));
        bors(P,4) = A + 4*(N-2) + 4*(N-2)^2 + 1 + (N-2)*(P-(2*N+(N-1)));
    end

    for P = 2*N+2*(N-1):(N-1)*5
        bors(P,1) = P - 1;
        bors(P,2) = P + 1;
        bors(P,3) = A + 4*(N-2) + (N-2)^2-(N-2)+P-(2*N+2*(N-1)-1);
        bors(P,4) = PointsPerCube - (N-2)^2 + P - (2*N+2*(N-1)-1);
    end
    for P = N+2*(N-1)+(N-2)+(N-1):2*N+2*(N-1)+2*(N-2)
        bors(P,1) = P - 1;
        bors(P,2) = P + 1;
        bors(P,3) = A + 4*(N-2) + 2*(N-2)^2 - (N-2) + P - (N+2*(N-1)+(N-2)+N);
        bors(P,4) = PointsPerCube - (N-2)^2 + (N-2) + (N-2)*(P-(N+2*(N-1)+(N-2)+(N+1)));
    end
    for P = 2*N+4*(N-1):2*N + 3*(N-1) + 2*(N-2)
        bors(P,1) = P - 1;
        bors(P,2) = P + 1;
        bors(P,3) = A + 4*(N-2)+3*(N-2)^2 - (N-2) + (P - (2*N+4*(N-1))+1);
        bors(P,4) = PointsPerCube - (P-(2*N+4*(N-1)));
    end
    for P = A-(N-3):A
        bors(P,1) = P - 1;
        bors(P,2) = P + 1;
        bors(P,3) = PointsPerCube - (N-3) - (N-2)*(P-(A-(N-3)));
        bors(P,4) = A + 4*(N-2) + 4*(N-2)^2 - (N-3) + (P-(A-(N-3)));
    end

    for P = A+1:A+(N-2)
        if(P==A+1)
            bors(P,1) = 1;
        else
            bors(P,1) = P - 1;
        end

        if(P==A+(N-2))
            bors(P,2) = A/2+1;
        else
            bors(P,2) = P + 1;
        end
        bors(P,3) = (A+4*(N-2)+1) + (N-2)*(P-(A+1));
        bors(P,4) = A + 4*(N-2) + 3*(N-2)^2 + (N-2) + (N-2)*(P-(A+1));
    end
    for P = A+(N-2)+1:A+2*(N-2)
        if (P==A+(N-2)+1)
            bors(P,1) = N;
        else
            bors(P,1) = P - 1;
        end
        if(P==A+2*(N-2))
            bors(P,2) = A/2+(N-1)+1;
        else
            bors(P,2) = P + 1;
        end
        bors(P,3) = A + 5*(N-2) + (N-2)*(P-(A+(N-2)+1));
        bors(P,4) = A + 4*(N-2) + (N-2)^2 + 1 + (N-2)*(P-(A+(N-2)+1));
    end
    for P = A+2*(N-2)+1:A+3*(N-2)
        if(P==A+2*(N-2)+1)
            bors(P,1) = N+(N-1);
        else
            bors(P,1) = P - 1;
        end
        if(P==A+3*(N-2))
            bors(P,2) = A/2 + N + (N-1);
        else
            bors(P,2) = P + 1;
        end
        bors(P,3) = A +4*(N-2)+(N-2)^2+(N-2)+(N-2)*(P-(A+2*(N-2)+1));
        bors(P,4) = A + 4*(N-2)+2*(N-2)^2 + 1 + (N-2)*(P-(A+2*(N-2)+1));
    end
    for P = A+3*(N-2)+1:A+4*(N-2)
        if(P==A+3*(N-2)+1)
            bors(P,1) = N+2*(N-1);
        else
            bors(P,1) = P - 1;
        end
        if(P==A+4*(N-2))
            bors(P,2) = A - (N-2);
        else
            bors(P,2) = P + 1;
        end
        bors(P,3) = A + 4*(N-2)+2*(N-2)^2+(N-2)+(N-2)*(P-(A+3*(N-2)+1));
        bors(P,4) = A +4*(N-2)+3*(N-2)^2+1+(N-2)*(P-(A+3*(N-2)+1));
    end

    for row = 2:N-1
        for col = 2:N-1
            bors(f1(row,col),1) = f1(row-1,col);
            bors(f1(row,col),2) = f1(row+1,col);
            bors(f1(row,col),3) = f1(row,col-1);
            bors(f1(row,col),4) = f1(row,col+1);
        end
    end
    for row = 2:N-1
        for col = 2:N-1
            bors(f2(row,col),1) = f2(row-1,col);
            bors(f2(row,col),2) = f2(row+1,col);
            bors(f2(row,col),3) = f2(row,col-1);
            bors(f2(row,col),4) = f2(row,col+1);
        end
    end
    for row = 2:N-1
        for col = 2:N-1
            bors(f3(row,col),1) = f3(row-1,col);
            bors(f3(row,col),2) = f3(row+1,col);
            bors(f3(row,col),3) = f3(row,col-1);
            bors(f3(row,col),4) = f3(row,col+1);
        end
    end
    for row = 2:N-1
        for col = 2:N-1
            bors(f4(row,col),1) = f4(row-1,col);
            bors(f4(row,col),2) = f4(row+1,col);
            bors(f4(row,col),3) = f4(row,col-1);
            bors(f4(row,col),4) = f4(row,col+1);
        end
    end
    for row = 2:N-1
        for col = 2:N-1
            bors(f5(row,col),1) = f5(row-1,col);
            bors(f5(row,col),2) = f5(row+1,col);
            bors(f5(row,col),3) = f5(row,col-1);
            bors(f5(row,col),4) = f5(row,col+1);
        end
    end
    for row = 2:N-1
        for col = 2:N-1
            bors(f6(row,col),1) = f6(row-1,col);
            bors(f6(row,col),2) = f6(row+1,col);
            bors(f6(row,col),3) = f6(row,col-1);
            bors(f6(row,col),4) = f6(row,col+1);
        end
    end
end
function [bors_T] = generate_neighbors(bors,PointsPerCube,NumObj,connected_faces,f1,f2,f3,f4,f5,f6,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs: bors, PointsPerCube, NumObj
    % Outputs: bors_T, local
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    local = ceil(bors(:,:)/(PointsPerCube+1))*PointsPerCube;
    %All local neighbors in system as a template
    bors_T = zeros(PointsPerCube*NumObj,7);
    %Applying local node neighbors throughout the system
    for i = 1:NumObj
        if(i==1)
            bors_T(1+PointsPerCube*(i-1):i*PointsPerCube,:) = bors(:,:);
        else
              bors_T(1+PointsPerCube*(i-1):i*PointsPerCube,:) = bors(:,:)+local*(i-1);
        end
    end
    
    for object = 1:NumObj
        for face = 1:6
            %Face 1
            if (face == 1 && connected_faces(object,face) ~= 0)
                %Face 3 of object
                %Knowns object face 1, touching object face 3
                obj_face_1 = f1 + (PointsPerCube*(object-1));
                obj_face_3 = f3 + (PointsPerCube*(connected_faces(object,face)-1));
                for row = 1:N
                    for col = 1:N
                        if(bors_T(obj_face_1(row,col),5) == 0)
                            bors_T(obj_face_1(row,col),5) = obj_face_3(row,col);
                        elseif(bors_T(obj_face_1(row,col),6) == 0)
                            bors_T(obj_face_1(row,col),6) = obj_face_3(row,col);
                        else
                            bors_T(obj_face_1(row,col),7) = obj_face_3(row,col);
                        end
                    end
                end
            end
            %Face 2
            if (face == 2 && connected_faces(object,face) ~= 0)
                %Face 4 of object
                obj_face_2 = f2 + (PointsPerCube*(object-1));
                obj_face_4 = f4 + (PointsPerCube*(connected_faces(object,face)-1));
                for row = 1:N
                    for col = 1:N
                        if(bors_T(obj_face_2(row,col),5) == 0)
                            bors_T(obj_face_2(row,col),5) = obj_face_4(row,col);
                        elseif(bors_T(obj_face_2(row,col),6) == 0)
                            bors_T(obj_face_2(row,col),6) = obj_face_4(row,col);
                        else
                            bors_T(obj_face_2(row,col),7) = obj_face_4(row,col);
                        end
                    end
                end
            end
            %Face 3
            if (face == 3 && connected_faces(object,face) ~= 0)
                %Face 1 of object
                obj_face_3 = f3 + (PointsPerCube*(object-1));
                obj_face_1 = f1 + (PointsPerCube*(connected_faces(object,face)-1));
                for row = 1:N
                    for col = 1:N
                        if(bors_T(obj_face_3(row,col),5) == 0)
                            bors_T(obj_face_3(row,col),5) = obj_face_1(row,col);
                        elseif(bors_T(obj_face_3(row,col),6) == 0)
                            bors_T(obj_face_3(row,col),6) = obj_face_1(row,col);
                        else
                            bors_T(obj_face_3(row,col),7) = obj_face_1(row,col);
                        end
                    end
                end
            end
            %Face 4
            if (face == 4 && connected_faces(object,face) ~= 0)
                %Face 2 of object
                obj_face_4 = f4 + (PointsPerCube*(object-1));
                obj_face_2 = f2 + (PointsPerCube*(connected_faces(object,face)-1));
                for row = 1:N
                    for col = 1:N
                        if(bors_T(obj_face_4(row,col),5) == 0)
                            bors_T(obj_face_4(row,col),5) = obj_face_2(row,col);
                        elseif(bors_T(obj_face_4(row,col),6) == 0)
                            bors_T(obj_face_4(row,col),6) = obj_face_2(row,col);
                        else
                            bors_T(obj_face_4(row,col),7) = obj_face_2(row,col);
                        end
                    end
                end
            end
            %Face 5
            if (face == 5 && connected_faces(object,face) ~= 0)
                %Face 6 of object
                obj_face_5 = f5 + (PointsPerCube*(object-1));
                obj_face_6 = f6 + (PointsPerCube*(connected_faces(object,face)-1));
                for row = 1:N
                    for col = 1:N
                        if(bors_T(obj_face_5(row,col),5) == 0)
                            bors_T(obj_face_5(row,col),5) = obj_face_6(row,col);
                        elseif(bors_T(obj_face_5(row,col),6) == 0)
                            bors_T(obj_face_5(row,col),6) = obj_face_6(row,col);
                        else
                            bors_T(obj_face_5(row,col),7) = obj_face_6(row,col);
                        end
                    end
                end
            end
            %Face 6
            if (face == 6 && connected_faces(object,face) ~= 0)
                %Face 5 of object
                obj_face_6 = f6 + (PointsPerCube*(object-1));
                obj_face_5 = f5 + (PointsPerCube*(connected_faces(object,face)-1));
                for row = 1:N
                    for col = 1:N
                        if(bors_T(obj_face_6(row,col),5) == 0)
                            bors_T(obj_face_6(row,col),5) = obj_face_5(row,col);
                        elseif(bors_T(obj_face_6(row,col),6) == 0)
                            bors_T(obj_face_6(row,col),6) = obj_face_5(row,col);
                        else
                            bors_T(obj_face_6(row,col),7) = obj_face_5(row,col);
                        end
                    end
                end
            end
        end
    end

end
function [num_bors] = generate_NL_neighbors(PointsPerCube,NumObj,bors_T)
    %Generating non-local neighbors for each cube

    num_bors = zeros(PointsPerCube*NumObj,1);
    %Generating 
    for i = 1:PointsPerCube*NumObj
        j = 0;
        if bors_T(i,5) ~= 0
            j = j + 1;
        end
        if bors_T(i,6) ~= 0
            j = j + 1;
        end
        if bors_T(i,7) ~= 0
            j = j + 1;
        end
        num_bors(i) = j;
    end
end
function [sun_faces] = sunCheckCorner(sun,connected_faces,block,corners)
    sun_faces = 0;
    if sun(corners(1)) == 1 && connected_faces(block,corners(1)) == 0  %light from sun
        sun_faces = sun_faces + 1;
    end
    if sun(corners(2)) == 1 && connected_faces(block,corners(2)) == 0
        sun_faces = sun_faces + 1;
    end
    if sun(corners(3)) == 1 && connected_faces(block,corners(3)) == 0
        sun_faces = sun_faces + 1;
    end
end
function [earth_faces] = earthCheckCorner(earth,connected_faces,block,corners)
    earth_faces = 0;
    if earth(corners(1)) == 1 && connected_faces(block,corners(1)) == 0  %light from earth
        earth_faces = earth_faces + 1;
    end
    if earth(corners(2)) == 1 && connected_faces(block,corners(2)) == 0
        earth_faces = earth_faces + 1;
    end
    if earth(corners(3)) == 1 && connected_faces(block,corners(3)) == 0
        earth_faces = earth_faces + 1;
    end
end
function [space_faces] = spaceCheckCorner(connected_faces,block,corners)
    space_faces = 0;
    if connected_faces(block,corners(1)) == 0
        space_faces = space_faces + 1;
    end
    if connected_faces(block,corners(2)) == 0
        space_faces = space_faces + 1;
    end
    if connected_faces(block,corners(3)) == 0
        space_faces = space_faces + 1;
    end
end
function [sun_faces] = sunCheckEdge(sun,connected_faces,block,edges)
    sun_faces = 0;
    if sun(edges(1)) == 1 && connected_faces(block,edges(1)) == 0 %light from sun
            sun_faces = sun_faces + 1;
    end
    if sun(edges(2)) == 1 && connected_faces(block,edges(2)) == 0
            sun_faces = sun_faces + 1;
    end
end
function [earth_faces] = earthCheckEdge(earth,connected_faces,block,edges)
    earth_faces = 0;
    if earth(edges(1)) == 1 && connected_faces(block,edges(1)) == 0 %light from earth
            earth_faces = earth_faces + 1;
    end
    if earth(edges(2)) == 1 && connected_faces(block,edges(2)) == 0
            earth_faces = earth_faces + 1;
    end
end
function [space_faces] = spaceCheckEdge(connected_faces,block,edges)
    space_faces = 0;
    if connected_faces(block,edges(1)) == 0
        space_faces = space_faces + 1;
    end
    if connected_faces(block,edges(2)) == 0
        space_faces = space_faces + 1;
    end
end
function [E_bool_sun] = sunCheckFace(sun,connected_faces,block,face)
    E_bool_sun = 0;
    if sun(face) == 1 && connected_faces(block,face) == 0
        E_bool_sun = 1;
    end
end
function [E_bool_earth] = earthCheckFace(earth,connected_faces,block,face)
    E_bool_earth = 0;
    if earth(face) == 1 && connected_faces(block,face) == 0
        E_bool_earth = 1;
    end
end
function [E_bool_space] = spaceCheckFace(connected_faces,block,face)
    E_bool_space = 0;    
    if connected_faces(block,face) == 0
        E_bool_space = 1;
    end
end
function [N_Temperature] = cornerTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces)
    if num_bors(i) == 3
        tau_x = (dt*k(block))/(rho(block)*Cp(block)*dx^2);
        tau_z = (dt*k(block))/(rho(block)*Cp(block)*dz^2);
        N_Temperature = (4/3)*tau_x*(T(j,bors_T(i,1))+ T(j,bors_T(i,2))+ T(j,bors_T(i,3))-3*T(j,i))...
                + (1/3)*tau_z*(T(j,bors_T(i,5)) + T(j,bors_T(i,6)) + T(j,bors_T(i,7)) - 3*T(j,i))...
                + dt/(rho(block)*Cp(block))*edot*edot_obj(block)...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_sun*sun_faces...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_earth*earth_faces...
                + (Emiss(block)*delta*dt)/(3*rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*space_faces...
                + T(j,i);
    elseif num_bors(i) == 2
        tau_x = (dt*k(block))/(rho(block)*Cp(block)*dx^2);
        tau_z = (dt*k(block))/(rho(block)*Cp(block)*dz^2);
        N_Temperature = (4/3)*tau_x*(T(j,bors_T(i,1))+ T(j,bors_T(i,2))+ T(j,bors_T(i,3))-3*T(j,i))...
                + (1/3)*tau_z*(T(j,bors_T(i,5)) + T(j,bors_T(i,6)) - 2*T(j,i))...
                + dt/(rho(block)*Cp(block))*edot*edot_obj(block)...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_sun*sun_faces...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_earth*earth_faces...
                + (Emiss(block)*delta*dt)/(3*rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*space_faces...
                + T(j,i);
    elseif num_bors(i) == 1
        tau_x = (dt*k(block))/(rho(block)*Cp(block)*dx^2);
        tau_z = (dt*k(block))/(rho(block)*Cp(block)*dz^2);
        N_Temperature = (4/3)*tau_x*(T(j,bors_T(i,1))+ T(j,bors_T(i,2))+ T(j,bors_T(i,3))-3*T(j,i))...
                + (1/3)*tau_z*(T(j,bors_T(i,5)) - T(j,i))...
                + dt/(rho(block)*Cp(block))*edot*edot_obj(block)...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_sun*sun_faces...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_earth*earth_faces...
                + (Emiss(block)*delta*dt)/(3*rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*space_faces...
                + T(j,i);                   
    else
        tau_x = (dt*k(block))/(rho(block)*Cp(block)*dx^2);
        N_Temperature = (4/3)*tau_x*(T(j,bors_T(i,1))+ T(j,bors_T(i,2))+ T(j,bors_T(i,3))-3*T(j,i))...
                + dt/(rho(block)*Cp(block))*edot*edot_obj(block)...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_sun*sun_faces...
                + dt/(3*rho(block)*Cp(block)*dz)*Q_earth*earth_faces...
                + (Emiss(block)*delta*dt)/(3*rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*space_faces...
                + T(j,i);
    end

end
function [N_Temperature] = edgeTemp(i,j,block,bors_T,num_bors,T,dt,dx,dy,dz,delta,Cp,k,rho,Emiss,Q_sun,sun_faces,Q_earth,earth_faces,edot,edot_obj,T_space,space_faces)
    if num_bors(i) == 2
        tau_x = (k(block)*dt)/(rho(block)*Cp(block)*dx^2);
        tau_y = (k(block)*dt)/(rho(block)*Cp(block)*dy^2);
        tau_z = (k(block)*dt)/(rho(block)*Cp(block)*dz^2);
        %point 1 and 2 are along the edge
        N_Temperature = tau_x*(T(j,bors_T(i,1))+T(j,bors_T(i,2))-2*T(j,i))...
                        + tau_y/2*(T(j,bors_T(i,3))+T(j,bors_T(i,4))-2*T(j,i))...
                        + tau_z/2*(T(j,bors_T(i,5)) + T(j,bors_T(i,6)) -2*T(j,i))...
                        + dt/(2*rho(block)*Cp(block)*dz)*Q_earth*earth_faces...
                        + dt/(2*rho(block)*Cp(block)*dz)*Q_sun*sun_faces...
                        + dt/(rho(block)*Cp(block))*edot*edot_obj(block)...
                        + (Emiss(block)*delta*dt)/(2*rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*space_faces...
                        + T(j,i);        
    elseif num_bors(i) == 1
        tau_x = (k(block)*dt)/(rho(block)*Cp(block)*dx^2);
        tau_y = (k(block)*dt)/(rho(block)*Cp(block)*dy^2);
        tau_z = (k(block)*dt)/(rho(block)*Cp(block)*dz^2);            
        N_Temperature = tau_x*(T(j,bors_T(i,1))+T(j,bors_T(i,2))-2*T(j,i))...
                        + tau_y/2*(T(j,bors_T(i,3))+T(j,bors_T(i,4))-2*T(j,i))...
                        + tau_z/2*(T(j,bors_T(i,5))-T(j,i))...
                        + dt/(2*rho(block)*Cp(block)*dz)*Q_earth*earth_faces...
                        + dt/(2*rho(block)*Cp(block)*dz)*Q_sun*sun_faces...
                        + dt/(rho(block)*Cp(block))*edot*edot_obj(block)...
                        + (Emiss(block)*delta*dt)/(2*rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*space_faces...
                        + T(j,i);
    else
        tau_x = (k(block)*dt)/(rho(block)*Cp(block)*dx^2);
        tau_y = (k(block)*dt)/(rho(block)*Cp(block)*dy^2);
        tau_z = (k(block)*dt)/(rho(block)*Cp(block)*dz^2);      
        N_Temperature = tau_x*(T(j,bors_T(i,1))+T(j,bors_T(i,2))-2*T(j,i))...
                        + tau_y/2*(T(j,bors_T(i,3))+T(j,bors_T(i,4))-2*T(j,i))...
                        + dt/(2*rho(block)*Cp(block)*dz)*Q_earth*earth_faces...
                        + dt/(2*rho(block)*Cp(block)*dz)*Q_sun*sun_faces...
                        + dt/(rho(block)*Cp(block))*edot*edot_obj(block)...
                        + (Emiss(block)*delta*dt)/(2*rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*space_faces...
                        + T(j,i);
    end
end
function [N_temperature] = faceTemp(i,j,block,bors_T,num_bors,T,dt,dx,dz,delta,Cp,k,rho,Emiss,Q_sun,E_bool_sun,Q_earth,E_bool_earth,edot,edot_obj,T_space,E_bool_space)
    if num_bors(i) == 1
        tau = (k(block)*dt)/(rho(block)*Cp(block)*dx^2);
        tau_z = (k(block)*dt)/(rho(block)*Cp(block)*dz^2);
        N_temperature = tau*(T(j,bors_T(i,1))+T(j,bors_T(i,2))+T(j,bors_T(i,3))+T(j,bors_T(i,4))-4*T(j,i))...
                        + tau_z*(T(j,bors_T(i,5))-T(j,i))...
                        +(2*edot*dt)/(rho(block)*Cp(block))*edot_obj(block)...
                        + (dt/(dz*Cp(block)*rho(block)))*(Q_sun*E_bool_sun + Q_earth*E_bool_earth)...
                        + (Emiss(block)*delta*dt)/(rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*E_bool_space...
                        + T(j,i);
    else
        tau = (k(block)*dt)/(rho(block)*Cp(block)*dx^2);
        N_temperature = tau*(T(j,bors_T(i,1))+T(j,bors_T(i,2))+T(j,bors_T(i,3))+T(j,bors_T(i,4))-4*T(j,i))...
                        +(2*edot*dt)/(rho(block)*Cp(block))*edot_obj(block)...
                        + (dt/(dz*Cp(block)*rho(block)))*(Q_sun*E_bool_sun + Q_earth*E_bool_earth)...
                        + (Emiss(block)*delta*dt)/(rho(block)*Cp(block)*dz)*(T_space^4-T(j,i)^4)*E_bool_space...
                        + T(j,i);
    end
end





















