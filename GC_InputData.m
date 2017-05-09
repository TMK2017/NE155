function [ Dgrid,Sgrid,sigma,x,y,flag] = GC_InputData(userx,usery,userD,userS,usersigma)
%Read and/or process the input data. Check all input values
%for correctness/sensibility. Print an error message and terminate
%if one or more errors occurs, otherwise print notification of
%successful input checking

%% Check X 
if isa(userx,'double')
    if isvector(userx)
        if issorted(userx)
            usern = length(userx)-1;
            x = transpose(userx(:));
        else
            flag=1;
            error('X must be sorted (ordered from X_min to X_max)')
        end
    else
        flag = 1;
        error('X must be a vector')
    end
else
    flag = 1;
    error('X must be of class double')
end


%% Check Y
if isa(usery,'double')
    if isvector(usery)
        if issorted(usery)
            userm = length(usery)-1;
            y = transpose(usery(:));
        else
            flag=1;
            error('Y must be sorted (ordered from Y_min to Y_max)')
        end
    else
        flag = 1;
        error('Y must be a vector')
    end
else
    flag = 1;
    error('Y must be of class double')
end


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


del = x(2:end) - x(1:end-1);
epsil = y(2:end) - y(1:end-1);


%% Check diffusion coefficient to see if it is a constant, matrix, or function
if isa(userD,'function_handle')
    Dgrid = zeros(n,m);
    try
        for i = 1:n
            for j = 1:m
                Dgrid(i,j) = userD(x(i)+del(i)/2,y(j)+epsil(j)/2);
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

% Check that D, S, sigma are same size
[nD,mD] = size(Dgrid);

%% Check the source to see if it is a constant, matrix, or function
if isa(userS,'function_handle')
    Sgrid = zeros(n,m);
    try
        for i = 1:n
            for j = 1:m
                Sgrid(i,j) = userS(x(i)+del(i)/2,y(j)+epsil(j)/2);
            end
        end
    catch
        flag = 1;
        error('If S is a function handle, S(x,y) should be a function of x and y.');
    end
    userS = Sgrid;
end
    
if isa(userS,'double')
        if all(userS > 0) & all(~isnan(userS)) & all(~isinf(userS))
            flag = 0;
            if max(size(userS)) == 1
                Sgrid = ones(n,m)*userS;
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

% Check that D, S, sigma are same size
[nS,mS] = size(Sgrid);


%% Check absorption coefficient.
if isa(usersigma,'function_handle')
    sigma = zeros(n,m);
    try
        for i = 1:n
            for j = 1:m
                sigma(i,j) = usersigma(x(i)+del(i)/2,y(j)+epsil(j)/2);
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

% Check that D, S, sigma are same size
[nsig,msig] = size(sigma);

if ~isequal(n,nD,nS,nsig) | ~isequal(m,mD,mS,msig)
    flag = 1;
    error('Dimensions of D, S, and sigma must be equal')
end
    


end

