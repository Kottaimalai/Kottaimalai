function R=GCC(matrix,link_arrange,np)
[p,q]=size(matrix);
F_NA=arrange(matrix);
alpha_L=1;
alpha_Theta=0.2;
 beta=3;
vL=1.00;
vTheta=20;
L=zeros(p,q);
Y=zeros(p,q);
Y0=zeros(p,q);
Theta=zeros(p,q);
center_x=round(link_arrange/2);
center_y=round(link_arrange/2);
W=zeros(link_arrange,link_arrange);
for i=1:link_arrange
    for j=1:link_arrange
        if (i==center_x)&&(j==center_y)
            W(i,j)=0;
        else
            W(i,j)=1./sqrt((i-center_x).^2+(j-center_y).^2);
        end
    end
end
F=F_NA;
for n=1:np
    K=conv2(Y,W,'same');
    L=exp(-alpha_L)*L+vL*K;
    Theta=exp(-alpha_Theta)*Theta+vTheta*Y;
    U=F.*(1+beta*L);
    Y=im2double(U>Theta);
    Y0=Y0+Y;   
end
R=Y0;
