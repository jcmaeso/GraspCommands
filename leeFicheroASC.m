
function [ THin,PHin,e_th,e_ph ] = leeFicheroASC( sFich_asc )
%LEEFICHEROASC Summary of this function goes here
%   Detailed explanation goes here
[Cmptes, Args, frec, Formato_in, arg3_in, arg4_in, ...
             datos] = LeerASC(sFich_asc);
         [THin, PHin]= meshgrid(arg3_in*pi/180, arg4_in*pi/180); 
         %Mallas de los Datos adquiridos
 temp= zeros(length(THin(:,1)), length(PHin(1,:)));
 if size(datos,2) == 6
    temp(:)= datos(:,3);
    c11= temp;
    temp(:)= datos(:,4);
    c12= temp;
    temp(:)= datos(:,5);
    c21= temp;
    temp(:)= datos(:,6);
    c22= temp;
 else
    temp(:)= datos(:,3);
    c11= temp;
    temp(:)= datos(:,4);
    c12= temp;
    c21= zeros(size(datos,1),1);
    c22= zeros(size(datos,1),1);
 end

 %Formacion de Componentes segun Formato de Entrada
 [c1, c2] = Formato(c11, c12, c21, c22, Formato_in);

 e_th= c1;
 e_ph= c2;



end

function [Cmptes, Args, frec, Formato_in, arg3_in, arg4_in, ...
             datos] = LeerASC(sFich_asc)

%
%LeerASC
%   Analiza los ficheros *.ASC con los datos numericos de la
%   medida de PROCENCA v3.0.
%
%   [Cmptes, Args, frec, Formato_in, arg3_in, arg4_in, ...
%       datos] = LeerASC(sFich_asc)
%
%       Cmptes: vector de tipo cell con las etiquetas de las
%           componentes del fichero de medida.
%       Args: vector de tipo cell con las etiquetas de las
%           coordenadas del fichero de medida.
%       frec: frecuencia en GHz.
%       Formato_in: 'REIM' o 'DBG'. Indica el formato de entrada
%           de las columnas del fichero de medida.
%       arg3_in, arg4_in: vectores con la lista de coordenadas.
%       datos: matriz con las columnas del fichero de medida.
%
%       sFichero: fichero (*.asc) de medida.
%


% Leemos todos los datos del fichero de medida
fid= fopen(sFich_asc);

buffer= fscanf(fid,'%s',9);
n3 = fscanf(fid,'%d',1); %nº de muestras en argumento 3
buffer= fscanf(fid,'%s',9);
n4 = fscanf(fid,'%d',1); %nº de muestras en argumento 4
buffer= fscanf(fid,'%s',19);
Arg3= fscanf(fid,'%s',1); %etiqueta de la 1ª columna
buffer= fscanf(fid,'%s',11);
Arg4= fscanf(fid,'%s',1); %etiqueta de la 2ª columna
Args= [cellstr(Arg3) cellstr(Arg4)];
buffer= fscanf(fid,'%s',7);
frec= fscanf(fid,'%s',1); %frecuencia
buffer= fscanf(fid,'%s',2);
Cmpte1= fscanf(fid,'%s',1); %Componente 1
ncol= 6;
switch upper(Cmpte1)
    case 'FF'
        Cmpte1= 'FF COMPONENTE A'; %Componente 1
        Cmpte2= 'FF COMPONENTE E'; %Componente 2
        buffer= fscanf(fid,'%s',7);
    case 'COMPONENTE'
        buffer= fscanf(fid,'%s',1);
        switch upper(buffer)
            case 'TETA'
                Cmpte1= 'COMPONENTE TETA'; %Componente 1
                Cmpte2= 'COMPONENTE FI';   %Componente 2
                buffer= fscanf(fid,'%s',4);
            case 'CP-X'
                Cmpte1= 'COMPONENTE CP-X'; %Componente 1
                Cmpte2= 'COMPONENTE CP-Y'; %Componente 2
                buffer= fscanf(fid,'%s',4);
            case 'RH'
                Cmpte1= 'COMPONENTE RH'; %Componente 1
                Cmpte2= 'COMPONENTE LH'; %Componente 2
                buffer= fscanf(fid,'%s',4);
        end
    case 'COMP.'
        Cmpte1= 'COMP. EJE-MAYOR'; %Componente 1
        Cmpte2= 'COMP. EJE-MENOR'; %Componente 2
        buffer= fscanf(fid,'%s',5);
    case 'CP'
        buffer= fscanf(fid,'%s',1);
        bHayPtos= find(buffer == '.');
        if bHayPtos,
            %Estamos analizando un fichero de galibo
            ncol= 5; 
            Cmpte1= 'CP_GAL'; %Componente 1
            Cmpte2= 'XP_GAL'; %Componente 2
        else
            %Estamos analizando un fichero de un medida
            Cmpte2= buffer; %Componente 2
        end
        buffer= fscanf(fid,'%s',2);
    otherwise
        Cmpte2= fscanf(fid,'%s',1); %Componente 2
        buffer= fscanf(fid,'%s',2);
end
Cmptes= [cellstr(Cmpte1) cellstr(Cmpte2)];
if length(char(Cmptes(2)))> 15   %si es mayor que 5 es porque no hay
                                %segunda componente y son muchos guiones
    ncol = 4;
    buffer= fscanf(fid,'%s',1); %formato
else
    buffer= fscanf(fid,'%s',2);
    buffer= fscanf(fid,'%s',1); %formato
end


switch upper(buffer)
    case 'P.'
        %Si la etiqueta de la columna empieza con 'P.'
        %el fichero viene en formato RE / IM
        Formato_in= 'REIM';
        if ncol == 6
            buffer= fscanf(fid,'%s',13);
        else
            buffer= fscanf(fid,'%s',7);
        end
    otherwise
        %En caso contrario consideraremos dB / º
        Formato_in= 'DBG';
        if ncol == 6,
            buffer= fscanf(fid,'%s',9);
        elseif ncol == 4,
            buffer= fscanf(fid,'%s',5);
        end
end

%Almacenamos todos los datos numericos
if ncol == 6 || ncol==4,
%     datos= fscanf(fid,'%f',[ncol,n3*n4]);
    datos = zeros(n3*n4,ncol);
    for p = 1:n3*n4
        for r = 1:ncol
            datos(p,r) = fscanf (fid,'%f',1);
            if datos(p,r) == .222507386
                buf = fscanf(fid,'%f',1);
            end
        end
    end
    
elseif ncol == 5,
    %Estamos analizando un fichero de galibo
    temp= [];
    for n=1:n3,
        for m=1:n4,
            temp= [temp; fscanf(fid,'%f',[ncol-2,1])];
            buffer= fscanf(fid,'%s',[1,1]); %Ignoramos la 4ª columna
        end
    end
    datos= zeros(ncol-2,n3*n4);
    datos(:)= temp;
    clear temp

end
% datos= datos';
 [f,c] = find(datos==.222507386);
 if f~=0
     for p=1:size(f)
         if f(p)>n4
             datos(f(p),c(p)) = datos(f(p)-n4,c(p));
         else
             datos(f(p),c(p)) = datos(f(p)+n4,c(p));
         end 
     end
 end

fclose(fid);

arg3_in= datos(:,1);
arg3_in= arg3_in(1:n4:length(arg3_in)); %lista de arg3 adquiridas

arg4_in= datos(:,2);
arg4_in= arg4_in(1:n4); %lista de arg4 adquiridas
end


function [c1out, c2out] = Formato(c11, c12, c21, c22, Format_in)

%
%Formato
%   Forma las dos componentes de la medida en funcion del
%   formato (real-imag. o dB-grados).
%
%   [c1out, c2out] = Formato(c11, c12, c21, c22, Format_in)
%
%       c1out, c2out: matrices con las componentes de la medida
%           en formato real-imag.
%
%       c11, c12: matrices con la tercera y cuarta columna del
%           fichero de medida respectivamente.
%       c21, c22: matrices con la quinta y sexta columna del
%           fichero de medida respectivamente.
%       Format_in: 'REIM' o 'DBG'. Indica el formato de entrada
%           de las matrices anteriores.
%

switch upper(Format_in)
    %Segun el formato de entrada
    %formamos las componentes
    case 'REIM'
    %Z= X+j*Y
        c1out= c11+j*c12;
        c2out= c21+j*c22;
    case 'DBG'
    %Z= R.*exp(j*fase)
        c1out= 10.^(c11/20).*exp(j*c12*pi/180);
        c2out= 10.^(c21/20).*exp(j*c22*pi/180);
end
end