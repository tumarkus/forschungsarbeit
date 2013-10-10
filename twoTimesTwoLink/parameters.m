function [ PAR ] = parameters( PROP )
%parameters Calculates the PAR-Matrix for the links and the payload
%   Detailed explanation goes here
xs=0;
for i=1:size(PROP,1)
    xs=xs+PROP(i,1)*PROP(i,2);
end
J=0;
for i=1:size(PROP,1)
   J=J+PROP(i,2)*(PROP(i,1)-xs)^2; 
end
J=J/sum(PROP(:,2))*PROP(end,4)*PROP(end,3)^2;
xs=xs/sum(PROP(:,2))*PROP(end,3);
PAR=[xs,J];
end

