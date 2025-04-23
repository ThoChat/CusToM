function [MC]=CalcAngMom(Human_model,j,c_global)
if j==0
    MC=0;
else
    c1 = Human_model(j).R*Human_model(j).c;
    c_j = c1 + Human_model(j).p;
    v = Human_model(j).v0 + cross(Human_model(j).w,c_j);
    quantite_mvt = Human_model(j).m*v;
    I=Human_model(j).R*Human_model(j).I*Human_model(j).R'; % Inertia tensor
    MC = cross((c_j-c_global'),quantite_mvt) + I*Human_model(j).w;
    MC=MC+CalcAngMom(Human_model,Human_model(j).sister,c_global)+CalcAngMom(Human_model,Human_model(j).child,c_global);
end
end