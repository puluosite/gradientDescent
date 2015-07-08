clear all;
close all;
clc;

%% old gradient descent
err = 1.0;
lamda = 1e-4;
max_iter = 300;
iter = 1;
all_error = [];
all_y = [];

init_x = 0;
step = 0.1;% when step = 1, init_x = 0, wierd thing happens, f(x_new) - f(x) = 0 while none of them are opt
while (err > lamda && iter <= max_iter)
    [y, gradient] = func(init_x);
    all_y = [all_y, y];
    % fixed step
    new_x = init_x - gradient*step;
    [new_y, new_gradient] = func(new_x);
    new_err = abs(new_y - y)^2;
    all_error = [all_error, new_err];
    
    err = new_err;
    init_x = new_x;
end

plot(all_y,'bo-');
hold on;


%% backtracing line search
% still have to change the alpha beta for multiple times to see if x_final
% is the right one!
err = 1.0;
lamda = 1e-4;
max_iter = 300;
iter = 1;
all_error = [];
all_y = [];
% para for backtracing
% when alpha = 0.5, beta = 0.5, init_x = 0, only one iter will take to opt
alpha = 0.2; % 0 - 0.5
beta = 0.7;

init_x = 0;

while (err > lamda && iter <= max_iter)
    [y, gradient] = func(init_x);
    all_y = [all_y, y];
    % step is caled in the inner loop
    step = 1;
    % from convex optimization chap 9
    while (func(init_x - step*gradient) > func(init_x) - alpha*step*(gradient')*gradient)
        f_x = func(init_x - step*gradient)
        f_alpha_x = func(init_x) - alpha*step*(gradient')*gradient
        
        step = step * beta;
    end
    iter
    step
    new_x = init_x - gradient*step;
    [new_y, new_gradient] = func(new_x);
    new_err = abs(new_y - y)^2;
    all_error = [all_error, new_err];
  
    err = new_err;
    init_x = new_x;
end

plot(all_y,'ro-');

%% exact line search