function [ Dgrid,Sgrid,sigma,x,y,flag] = InputData(usern,userm,usera,userb,userD,userS,usersigma)

%Read and/or process the input data. Check all input values
%for correctness/sensibility. Print an error message and terminate
%if one or more errors occurs, otherwise print notification of
%successful input checking

%% Check if grid dimensions are greater 1.
if usern > 1
    n = usern;
else
    flag = 1;
    error('Grid dimension n must be greater than 1.')
end
if userm > 1
    m = userm;
else
    flag = 1;
    error('Grid dimension m must be greater than 1.')
end

x = linspace(-a,a,n);
y = linspace(-b,b,m);

%% Check if boundaries are positive.
if usera > 0
    a = usera;
else
    flag = 1;
    error('Boundary a must be positive.')
end

if userb > 0
    b = userb;
else
    flag = 1;
    error('Boundary b must be positive.')
end

%% Check diffusion coefficient to see if it is a constant, matrix, or function
if isa(userD,'function_handle')
    Dgrid = zeros(n,m);
    try
        for i = 1:n
            for j = 1:m
                Dgrid(i,j) = userD(x(i),y(j));
            end
        end
    catch
        flag = 1;
        error('If D is a function handle, D(x,y) should be a function of x and y.');
    end
    userD = Dgrid;
end
    
if isa(userD,'double')
        if all(userD > 0) & all(~isnan(userD)) & all(~isinf(userD))
            flag = 0;
            if max(size(userD)) == 1
                Dgrid = ones(n,m)*userD;
            elseif all(size(userD) == [n,m]) 
                Dgrid = userD;
            else
                flag = 1;
                error('Input diffusion coefficient has wrong dimensions.');
            end
        else
             flag = 1;
             error('The input diffusion coefficient must be a positive number.')
        end
else
	flag = 1;
    error('Diffusion coefficient input must be a double or a function handle.')  
end

%% Check the source to see if it is a constant, matrix, or function

    
if isa(userS,'double')
        if all(userS > 0) & all(~isnan(userS)) & all(~isinf(userS))
            flag = 0;
            if max(size(userS)) == 1
                Dgrid = ones(n,m)*userS;
            elseif all(size(userS) == [n,m]) 
                Sgrid = userS;
            else
                flag = 1;
                error('Input diffusion coefficient has wrong dimensions.');
            end
        else
             flag = 1;
             error('The input diffusion coefficient must be a positive number.')
        end
else
	flag = 1;
    error('Diffusion coefficient input must be a double or a function handle.')  
end

%% Check absorption coefficient.
if isa(usersigma,'function_handle')
    grid = zeros(n,m);
    try
        for i = 1:n
            for j = 1:m
                sigma(i,j) = usersigma(x(i),y(j));
            end
        end
    catch
        flag = 1;
        error('If D is a function handle, D(x,y) should be a function of x and y.');
    end
    usersigma = sigma;
end
    
if isa(usersigma,'double')
        if all(usersigma > 0) & all(~isnan(usersigma)) & all(~isinf(usersigma))
            flag = 0;
            if max(size(usersigma)) == 1
                sigma = ones(n,m)*usersigma;
            elseif all(size(usersigma) == [n,m]) 
                sigma = usersigma;
            else
                flag = 1;
                error('Input diffusion coefficient has wrong dimensions.');
            end
        else
             flag = 1;
             error('The input diffusion coefficient must be a positive number.')
        end
else
	flag = 1;
    error('Diffusion coefficient input must be a double or a function handle.')  
end
end

