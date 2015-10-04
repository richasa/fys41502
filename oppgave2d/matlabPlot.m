fid=fopen('text.txt')
    s = textscan(fid,'%f %f %f')
fclose(fid)

x=s{1};
y=s{2};
z=s{3};


plot(x,y,'red')
hold on
plot(x,z,'blue')


title('the behavior as wave funciton for wr = 5')
xlabel('the relative coordinat r')
ylabel('the probability distriubtion')
