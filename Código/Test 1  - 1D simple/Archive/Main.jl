using Plots



x = range(0,10,length=100)
y = sin.(x)

data = []

## Random data
for i = 0:10
    push!(data,0.1*randn())
end

println(data)

#plot(x,y)
