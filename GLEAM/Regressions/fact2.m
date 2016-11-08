function df=fact2(n)

if mod(n,2)
    df = 1;
    for i = n : -2 : 3
        df = df * i;
    end
    
else
    df = 2;
    for i = n : -2 : 4
        df = df * i;
    end
    
end

 
