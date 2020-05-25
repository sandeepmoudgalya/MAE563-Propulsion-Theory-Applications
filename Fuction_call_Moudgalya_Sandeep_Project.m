%For calling the main function
clc
clear
for n=0.1:0.1:2.5
    n
    o=Moudgalya_Sandeep_Project(n)
    plot(n,o,'o')
    hold on
end
title('Effect of Mach no at diffuser exit vs Thrust')
xlabel('Mach no at diffuser exit') 
ylabel('Thrust (N)') 

