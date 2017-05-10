%Test Cases
GC_Master(linspace(-1,1,100),linspace(-1,1,100),0.1,@(x,y) abs(cos(x))+2,2,'Test1.xlsx'); %Test 1
GC_Master( linspace(-1,1,100),linspace(-1,1,100), -0.5, 4, 2, 'Test2.xlsx'); %Test 2
GC_Master( linspace(-3,-1,500), linspace(0,1,100), 0.1, 7, 2, 'Test3.xlsx')%Test 3
GC_Master( linspace(-1,1,500),linspace(-1,1,900), 0.1, @(x,y) abs(cos(x)) + 2, 2,'Test4.xlsx')%Test 4