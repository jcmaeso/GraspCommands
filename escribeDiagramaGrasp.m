function [ ] = escribeDiagramaGrasp( Theta,Fi,Etheta,Ephi,fname,freqs)
%writes a radiation pattern in a file to be used as a feed in Grasp 8.
%Fernando
fid = fopen(fname,'w');
dTh=Theta(1,2)-Theta(1,1);
nTh=size(Theta,2);
for f = freqs
    for n=1:size(Theta,1)
    fprintf(fid,'Frecuencia %f GHz\n',f);
    fprintf(fid,'%18.10E%18.10E%5i%18.10E%5i%5i%5i\n',Theta(n,1)*180/pi,dTh*180/pi,nTh,Fi(n,1)*180/pi,1,1,2);
    for m=1:size(Theta,2)
    fprintf(fid,'%18.10E%18.10E%18.10E%18.10E\n',real(Etheta(n,m)),imag(Etheta(n,m)),real(Ephi(n,m)),imag(Ephi(n,m)));
    end
    end
end
fclose(fid);
end