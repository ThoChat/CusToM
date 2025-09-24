function [OsteoArticularModel]= UpperTrunkCervical_KinRhythm(OsteoArticularModel,k,Mass,AttachmentPoint,varargin)
%% Implementation of a model of the thoracic spine with 12 thoracic vertebrae
% Ce modèle inclut 12 solides associés aux vertèbres lombaires et articulés
% entre eux suivant 3 degrés de liberté en rotation.
% Auteurs : A. Schuster, G. Dumont, C. Pontonnier (2023)
% Special : - Contraintes cinématiques
%           - Contraintes géométriques
% Fichier résultat : '.C3D' full_constrains
%
%% Les données anthropométriques pour ce modèle sont issues des travaux de
% Bayoglu et al. (2017a,b). Les contraintes cinématiques sont issues des
% travaux de Alemi et al. (2021).
% (Twente spine model: a complete and coherent dataset for musculo-skeletal
% modeling of the lumbar region of the human spine ;
%  Twente spine model: a complete and coherent dataset for musculo-skeletal
% modeling of the thoracic and cervical regions of the human spine
%  Alemi, Mohammad Mehdi, Katelyn A. Burkhart, Andrew C. Lynch, Brett T. 
% Allaire, Seyed Javad Mousavi, Chaofei Zhang, Mary L. Bouxsein, et Dennis
% E. Anderson. « The Influence of Kinematic Constraints on Model 
% Performance During Inverse Kinematics Analysis of the Thoracolumbar Spine
% ». Frontiers in Bioengineering and Biotechnology 9 (2021): 688041.
% https://doi.org/10.3389/fbioe.2021.688041.)
%
%% Comit 2025, authors: Minh Nguyen, Aurélien Schuster, Aurélie Tomezzoli
% You can run this model with the markersets
% - Marker_setSpine
% - Marker_setSpinePOSSUP
%
% In this case, virtual markers corresponding to the posterosuperior corners of the vertebrae are used to drive the kinematic model, as if they were skin markers. 
% % Before running the analysis in CUSTOM, you MUST compute the predicted vertebral positions using the regressions available here:
% Tomezzoli, A., A. Agouram, B. Chalamet, et al. ‘Predicting Cervico-Thoraco-Lumbar Vertebra Positions from Cutaneous Markers: Combining Local Frame and Postural Predictors Improves Robustness to Posture’. Journal of Biomechanics 164 (February 2024): 111961. https://doi.org/10.1016/j.jbiomech.2024.111961.
% 
% You MUST use simultaneously these model parts:
% UpperTrunkCervical_KinRhythm + LowerTrunk_KinRhythm + Skull    OR
% UpperTrunkCervical_Limits    + LowerTrunk_Limits   + Skull.m  OR 
% UpperTrunkCervical_NoLimits  + LowerTrunk_NoLimits + Skull.m
%
% Dof outside the sagittal plane have been disabled.



list_solid={'T12_J1' 'T12_J2' 'T12'...
            'T11_J1' 'T11_J2' 'T11'...
            'T10_J1' 'T10_J2' 'T10'...
            'T9_J1' 'T9_J2' 'T9'...
            'T8_J1' 'T8_J2' 'T8'...
            'T7_J1' 'T7_J2' 'T7'...
            'T6_J1' 'T6_J2' 'T6'...
            'T5_J1' 'T5_J2' 'T5'...
            'T4_J1' 'T4_J2' 'T4'...
            'T3_J1' 'T3_J2' 'T3'...
            'T2_J1' 'T2_J2' 'T2'...
            'T1_J1' 'T1_J2' 'T1'...
            'RClavicle_J1' 'RClavicle_J2' 'RClavicle'...
            'LClavicle_J1' 'LClavicle_J2' 'LClavicle'...
            'C7_J1' 'C7_J2' 'C7'...
            'C6_J1' 'C6_J2' 'C6'...
            'C5_J1' 'C5_J2' 'C5'...
            'C4_J1' 'C4_J2' 'C4'...
            'C3_J1' 'C3_J2' 'C3'...
            'C2_J1' 'C2_J2' 'C2'...
            'C1_J1' 'C1_J2' 'C1'...
            };


    
%% solid numbering incremation

s=size(OsteoArticularModel,2)+1;  %#ok<NASGU> % number of the first solid
for i=1:size(list_solid,2)      % each solid numbering: s_"nom du solide"
    if i==1
        eval(strcat('s_',list_solid{i},'=s;'))
    else
        eval(strcat('s_',list_solid{i},'=s_',list_solid{i-1},'+1;'))
    end
end  
    
% find the number of the mother from the attachment point: 'attachment_pt'
if ~numel(AttachmentPoint)
    s_mother=0;
    pos_attachment_pt=[0 0 0]';
else
    test=0;
    for i=1:numel(OsteoArticularModel)
        for j=1:size(OsteoArticularModel(i).anat_position,1)
            if strcmp(AttachmentPoint,OsteoArticularModel(i).anat_position{j,1})
               s_mother=i;
               pos_attachment_pt=OsteoArticularModel(i).anat_position{j,2}+OsteoArticularModel(s_mother).c;
               test=1;
               break
            end
        end
        if i==numel(OsteoArticularModel) && test==0
            error([AttachmentPoint ' is no existent'])        
        end       
    end
    if OsteoArticularModel(s_mother).child == 0      % if the mother don't have any child
        OsteoArticularModel(s_mother).child = eval(['s_' list_solid{1}]);    % the child of this mother is this solid
    else
        [OsteoArticularModel]=sister_actualize(OsteoArticularModel,OsteoArticularModel(s_mother).child,eval(['s_' list_solid{1}]));   % recherche de la derni re soeur
    end
end      


%% Definition of anatomical landmarks for the segments(CoM wrt global RF ; Node wrt local RF on CoM)

%Matrice de rotation de la base du CT-Scan vers la base suivie par les recommendations de l'ISB
M = [0 -1 0 ;...
     0 0 1  ;...
    -1 0 0  ];

%Matrice de symétrie pour le changement gauche/droite

S = [1 0 0 ;...
     0 1 0 ;...
     0 0 -1];

%Coeff de regression pour le cadavre de 1m54 :

% k_new = k * (1.80/1.54);
k_new = k ;

%Position des articulations telle que la moyenne des deux corps de %%%%%%%%
%vertèbres (sup+inf)/2 dans le repère du CT-scan: %%%%%%%%%%%%%%%%%%%%%%%%%
%Calcul des longueurs des vertèbres (entre les deux corps vertébraux): %%%%

%Calcul des centres de liaisons comme la moyenne entre les deux extrémités
%sur les corps vertébraux :

% Lombaires %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calcul des centres de liaisons comme la moyenne entre les deux extrémités
%sur les corps vertébraux :

%L5
L5_inf = M*k_new*[-2.9904 -79.4992 692.74]'*0.001 ;
L5_sup = M*k_new*[-5.3757 -96.2936 709.719]'*0.001 ;
L5_length = norm(L5_sup-L5_inf);
L4_inf = M*k_new*[-3.4144 -100.365 710.328]'*0.001 ;
L5_L4CT = (L5_sup+L4_inf)/2;        

%L4
L4_sup = M*k_new*[-3.5914 -108.569 734.236]'*0.001;
L4_length = norm(L4_sup-L4_inf);
L3_inf = M*k_new*[-1.3649 -113.081 743.652]'*0.001 ;
L4_L3CT = (L4_sup+L3_inf)/2 ;

%L3
L3_sup = M*k_new*[-1.3364 -112.546 767.612]'*0.001;
L3_length = norm(L3_sup-L3_inf);
L2_inf = M*k_new*[-2.2401 -111.985 775.05]'*0.001 ;
L3_L2CT = (L3_sup+L2_inf)/2 ;

%L2
L2_sup = M*k_new*[-2.4577 -109.706 800.086]'*0.001;
L1_inf = M*k_new*[-2.0244 -104.134 807.88]'*0.001 ;
L2_L1CT = (L2_sup+L1_inf)/2;
 
%L1
L1_sup = M*k_new*[-1.2769 -101.421 832.204]'*0.001;
T12_inf = M*k_new*[-0.4485 -97.7597 839.097]'*0.001;
 
%L1
L1_sup = M*k_new*[-1.2769 -101.421 832.204]'*0.001;
T12_inf = M*k_new*[-0.4485 -97.7597 839.097]'*0.001;
L1_T12CT = (L1_sup+T12_inf)/2;

%T12
T12_sup = M*k_new*[-1.9339 -92.7972 862.658]'*0.001; %Centre plateau sup repère CT
T12_length = norm(T12_sup-T12_inf);
T11_inf = M*k_new*[-3.6646 -89.9496 867.425]'*0.001;
T12_T11CT = (T12_sup+T11_inf)/2;%Joint nodes (global)

%T11
T11_sup = M*k_new*[-8.0478 -88.3005 890.497]'*0.001;
T11_length = norm(T11_sup-T11_inf);
T10_inf = M*k_new*[-2.2081 -86.9715 894.928]'*0.001;
T11_T10CT = (T11_sup+T10_inf)/2;

%T10
T10_sup = M*k_new*[-5.9396 -86.0169 916.308]'*0.001;
T10_length = norm(T10_sup-T10_inf);
T9_inf = M*k_new*[-5.9926 -86.5605 920.522]'*0.001;
T10_T9CT = (T10_sup+T9_inf)/2;

%T9
T9_sup = M*k_new*[-8.578 -87.4164 941.138]'*0.001;
T9_length = norm(T9_sup-T9_inf);
T8_inf = M*k_new*[-8.1813 -87.8752 944.919]'*0.001;
T9_T8CT = (T9_sup+T8_inf)/2;

%T8
T8_sup = M*k_new*[-8.4684 -88.1977 965.049]'*0.001;
T8_length = norm(T8_sup-T8_inf);
T7_inf = M*k_new*[-8.5026 -91.1361 968.422]'*0.001;
T8_T7CT = (T8_sup+T7_inf)/2;

%T7
T7_sup = M*k_new*[-9.8847 -90.9386 989.018]'*0.001;
T7_length = norm(T7_sup-T7_inf);
T6_inf = M*k_new*[-9.4087 -97.3257 991.257]'*0.001;
T7_T6CT = (T7_sup+T6_inf)/2;

%T6
T6_sup = M*k_new*[-7.9465 -101.947 1009.79]'*0.001;
T6_length = norm(T6_sup-T6_inf);
T5_inf = M*k_new*[-7.103 -103.244 1012.33]'*0.001;
T6_T5CT = (T6_sup+T5_inf)/2;

%T5
T5_sup = M*k_new*[-5.3964 -109.896 1028.45]'*0.001;
T5_length = norm(T5_sup-T5_inf);
T4_inf = M*k_new*[-4.137 -112.299 1031.48]'*0.001;
T5_T4CT = (T5_sup+T4_inf)/2;

%T4
T4_sup = M*k_new*[-1.0067 -121.573 1045.26]'*0.001;
T4_length = norm(T4_sup-T4_inf);
T3_inf = M*k_new*[-1.0567 -124.563 1048.33]'*0.001;
T4_T3CT = (T4_sup+T3_inf)/2;

%T3
T3_sup = M*k_new*[2.3591 -134.426 1061.02]'*0.001;
T3_length = norm(T3_sup-T3_inf);
T2_inf = M*k_new*[-1.099 -138.551 1062.85]'*0.001;
T3_T2CT = (T3_sup+T2_inf)/2;

%T2
T2_sup = M*k_new*[5.4555 -148.744 1075.04]'*0.001;
T2_length = norm(T2_sup-T2_inf);
T1_inf = M*k_new*[6.2343 -151.013 1078.71]'*0.001;  
T2_T1CT = (T2_sup+T1_inf)/2;

%Entre T1 et C7
T1_sup = M*k_new*[11.2313 -159.726 1089.99]'*0.001;
T1_length = norm(T1_sup-T1_inf);
C7_inf = M*k_new*[14.35 -162.76 1092.81]'*0.001;
T1_C7CT = (T1_sup+C7_inf)/2;

%Entre C7 et C6
C7_sup = M*k_new*[19.53 -169.97 1102.67]'*0.001;
C7_length = norm(C7_sup-C7_inf);
C6_inf = M*k_new*[19.01 -171.51 1104.71]'*0.001;
C7_C6CT = (C7_sup+C6_inf)/2;

%Entre C6 et C5
C6_sup = M*k_new*[22.28 -177.59 1114.87]'*0.001;
C6_length = norm(C6_sup-C6_inf);
C5_inf = M*k_new*[22.41 -178.98 1116.93]'*0.001;
C6_C5CT = (C6_sup+C5_inf)/2;

%Entre C5 et C4
C5_sup = M*k_new*[26.44 -184.74 1126.51]'*0.001;
C5_length = norm(C5_sup-C5_inf);
C4_inf = M*k_new*[26.35 -185.55 1130.56]'*0.001;
C5_C4CT = (C5_sup+C4_inf)/2;

%Entre C4 et C3
C4_sup = M*k_new*[28.44 -191.74 1138.51]'*0.001;
C4_length = norm(C4_sup-C4_inf);
C3_inf = M*k_new*[28.29 -193.451 1145]'*0.001;
C4_C3CT = (C4_sup+C3_inf)/2;

%Entre C3 et C2
C3_sup = M*k_new*[28.9698 -194.677 1158.64]'*0.001;
C3_length = norm(C3_sup-C3_inf);
C2_inf = M*k_new*[29.1789 -193.565 1162.76]'*0.001;
C3_C2CT = (C3_sup+C2_inf)/2;

%Entre C2 et C1
C2_sup = M*k_new*[30.70935 -192.2005 1178.9]'*0.001;
C1_inf = M*k_new*[29.5063 -192.8585 1180.81]'*0.001;
C1_sup = M*k_new*[30.77 -197.2 1197.52]'*0.001;
C2_C1CT = (C1_inf+C2_sup)/2;

%Position du CoM du segment enfant par rapport au repère associé au segment parent
% avec les paramètres de Pearsall et al. (1996) (Chaîne cinématique)
% Dans Pearsall et al. on a en fonction de la longueur les coefficients :
% x transverse ;
% y anteropostérieur ;
% z longitudinal mais vers le bas : on a fait le changement pour 'y' au préalable


MassT1      = 2.7/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT2      = 2.6/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT3      = 3.3/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT4      = 3.1/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT5      = 3.2/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT6      = 3.2/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT7      = 3.4/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT8      = 3.6/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT9      = 3.8/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT10     = 4.8/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT11     = 5.0/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT12     = 6.0/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);

% Thoraciques %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%L1
L1_T12JointNode = [0 norm(L1_T12CT - L2_L1CT) 0]'; %Repère segment parent
%T12
centre_geom_T12 = [0 norm(M*k_new*[-0.88880004 -95.19249946 850.85706413]'*0.001 - L1_T12CT) 0]';%Repère CT-Scan
CoM_T12 = centre_geom_T12 + T12_length*[1.640	0.000	-0.004]'; %Données de Pearsal et al. traduites en pourcentage de la longueur du segment
CoM_T12CT = CoM_T12 + L1_T12CT;
T12_T11JointNode = [0 norm(T12_T11CT-L1_T12CT) 0]';

%T11
centre_geom_T11 = [0 norm(M*k_new*[-3.836850765 -90.205750485 879.28863221]'*0.001-T12_T11CT) 0]';%Repère CT-Scan
CoM_T11 = centre_geom_T11 + T11_length*[1.760	0.000	-0.003]'; 
CoM_T11CT = CoM_T11 + T12_T11CT;
T11_T10JointNode = [0 norm(T11_T10CT-T12_T11CT) 0]';


%T10
centre_geom_T10 = [0 norm(M*k_new*[-4.74510698 -88.076947225 905.66991757]'*0.001-T11_T10CT) 0]';%Repère CT-Scan
CoM_T10 = centre_geom_T10 + T10_length*[2.300	0.000	-0.001]'; 
CoM_T10CT = CoM_T10 + T11_T10CT;
T10_T9JointNode = [0 norm(T10_T9CT-T11_T10CT) 0]';


%T9
centre_geom_T9 = [0 norm(M*k_new*[-5.831965935 -87.096696205 930.94956846]'*0.001-T10_T9CT) 0]';
CoM_T9 = centre_geom_T9 + T9_length*[2.300	0.000	0.000]'; 
CoM_T9CT = CoM_T9 + T10_T9CT;
T9_T8JointNode = [0 norm(T9_T8CT-T10_T9CT) 0]';


%T8
centre_geom_T8 = [0 norm(M*k_new*[-7.09154992 -86.621049795 955.103456975]'*0.001-T9_T8CT) 0]';
CoM_T8 = centre_geom_T8 + T8_length*[2.300	0.000	0.000]'; 
CoM_T8CT = CoM_T8 + T9_T8CT;
T8_T7JointNode = [0 norm(T8_T7CT-T9_T8CT) 0]';


%T7
centre_geom_T7 = [0 norm(M*k_new*[-7.872599875 -90.458102525 978.58589888]'*0.001-T8_T7CT) 0]';
CoM_T7 = centre_geom_T7 + T7_length*[2.150	0.000	0.002]'; 
CoM_T7CT = CoM_T7 + T8_T7CT;
T7_T6JointNode = [0 norm(T7_T6CT-T8_T7CT) 0]';


%T6
centre_geom_T6 = [0 norm(M*k_new*[-7.65406297 -98.26729827 1000.734251955]'*0.001-T7_T6CT) 0]';
CoM_T6 = centre_geom_T6 + T6_length*[1.950	0.000	0.003]'; 
CoM_T6CT = CoM_T6 + T7_T6CT;
T6_T5JointNode = [0 norm(T6_T5CT-T7_T6CT) 0]';


%T5
centre_geom_T5 = [0 norm(M*k_new*[-6.150750095 -106.459300965 1020.276397465]'*0.001-T6_T5CT) 0]';
CoM_T5 = centre_geom_T5 + T5_length*[1.700	0.000	0.002]'; 
CoM_T5CT = CoM_T5 + T6_T5CT;
T5_T4JointNode = [0 norm(T5_T4CT-T6_T5CT) 0]';


%T4
centre_geom_T4 = [0 norm(M*k_new*[-1.886976905 -117.97578271 1037.791405535]'*0.001-T5_T4CT) 0]';
CoM_T4 = centre_geom_T4 + T4_length*[1.400	0.000	0.001]';
CoM_T4CT = CoM_T4 + T5_T4CT;
T4_T3JointNode = [0 norm(T4_T3CT-T5_T4CT) 0]';


%T3
centre_geom_T3 = [0 norm(M*k_new*[2.029917365 -130.95841247 1053.39871453]'*0.001-T4_T3CT) 0]';
CoM_T3 = centre_geom_T3 + T3_length*[1.333	0.000	-0.002]'; 
CoM_T3CT = CoM_T3 + T4_T3CT;
T3_T2JointNode = [0 norm(T3_T2CT-T4_T3CT) 0]';


%T2
centre_geom_T2 = [0 norm(M*k_new*[4.85567504 -142.766151575 1069.3564713]'*0.001-T3_T2CT) 0]';
CoM_T2 = centre_geom_T2 + T2_length*[0.867	0.000	-0.005]'; 
CoM_T2CT = CoM_T2 + T3_T2CT;
T2_T1JointNode = [0 norm(T2_T1CT-T3_T2CT) 0]';

%T1
centre_geom_T1 = [0 norm(M*k_new*[9.70759988 -154.560297725 1084.445774555]'*0.001-T2_T1CT) 0]';%Repère CT-Scan
CoM_T1 = centre_geom_T1 + T1_length*[0.800	0.000	-0.006]'; 
CoM_T1CT = CoM_T1 + T2_T1CT;
T1_C7JointNode = [0 norm(T1_C7CT-T2_T1CT) 0]';

%C7
centre_geom_C7 = [0 norm(M*k_new*[17.232949845 -165.89352861 1097.79238701]'*0.001-T1_C7CT) 0]';%Repère CT-Scan
CoM_C7 = centre_geom_C7; 
CoM_C7CT = CoM_C7 + T1_C7CT;
C7_C6JointNode = [0 norm(C7_C6CT-T1_C7CT) 0]';

%C6
centre_geom_C6 = [0 norm(M*k_new*[22.199249825 -174.04229939 1109.831154345]'*0.001-C7_C6CT) 0]';%Repère CT-Scan
CoM_C6 = centre_geom_C6; 
CoM_C6CT = CoM_C6 + C7_C6CT;
C6_C5JointNode = [0 norm(C6_C5CT-C7_C6CT) 0]';

%C5
centre_geom_C5 = [0 norm(M*k_new*[23.95494934 -182.11659789 1121.572822335]'*0.001-C6_C5CT) 0]';%Repère CT-Scan
CoM_C5 = centre_geom_C5; 
CoM_C5CT = CoM_C6 + C6_C5CT;
C5_C4JointNode = [0 norm(C5_C4CT - C6_C5CT) 0]';

%C4
centre_geom_C4 = [0 norm(M*k_new*[27.027212395 -188.740900685 1136.014707235]'*0.001-C5_C4CT) 0]';%Repère CT-Scan
CoM_C4 = centre_geom_C4; 
CoM_C4CT = CoM_C4 + C5_C4CT;
C4_C3JointNode = [0 norm(C4_C3CT - C5_C4CT) 0]';

%C3
centre_geom_C3 = [0 norm(M*k_new*[28.592766545 -193.79058529 1151.91299904]'*0.001-C4_C3CT) 0]';%Repère CT-Scan
CoM_C3 = centre_geom_C3; 
CoM_C3CT = CoM_C3 + C4_C3CT;
C3_C2JointNode = [0 norm(C3_C2CT - C4_C3CT) 0]';

%C2
centre_geom_C2 = [0 norm(M*k_new*[28.463409795 -180.365311845 1171.135298765]'*0.001-C3_C2CT) 0]';%Repère CT-Scan
CoM_C2 = centre_geom_C2; 
CoM_C2CT = CoM_C2 + C3_C2CT;
C2_C1JointNode = [0 norm(C2_C1CT - C3_C2CT) 0]';

%C1
centre_geom_C1 = [0 norm(M*k_new*[31.24669939 -202.598270025 1191.16777182]'*0.001-C2_C1CT) 0]';%Repère CT-Scan
CoM_C1 = centre_geom_C1; 
CoM_C1CT = CoM_C1 + C2_C1CT;
C1_HeadJointNode = [0 norm(C1_sup - C2_C1CT) 0]';


%% Définition des théta pour la posture de la colonne : %%%%%%%%

%% Takada 
cyphose = -34.9; % cyphose de 45° supT1infT12 selon proportions de Kuntz et al. 2007
cyphose = cyphose*pi/180;

lordose_cer = 40;
lordose_cer = lordose_cer*pi/180;

theta_T12 = [0 0 0.024390244*cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T11 = [0 0 0.06097561 *cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T10 = [0 0 0.073170732*cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T9  = [0 0 0.073170732*cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T8  = [0 0 0.097560976*cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T7  = [0 0 0.12195122 *cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T6  = [0 0 0.12195122 *cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T5  = [0 0 0.12195122 *cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T4  = [0 0 0.12195122 *cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T3  = [0 0 0.085365854*cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T2  = [0 0 0.073170732*cyphose]'; %Orientation du vecteur par rapport à l'axe z
theta_T1  = [0 0 0.024390244*cyphose]'; %Orientation du vecteur par rapport à l'axe z

theta_C7 = [0 0 0.077*lordose_cer]';
theta_C6 = [0 0 0.145*lordose_cer]';
theta_C5 = [0 0 0.181*lordose_cer]';
theta_C4 = [0 0 0.179*lordose_cer]';
theta_C3 = [0 0 0.157*lordose_cer]';
theta_C2 = [0 0 0.117*lordose_cer]';
theta_C1 = [0 0 0.144*lordose_cer]';

%% Issu du modèle de bruno dernière version
% theta_T12 = [0 0 -0.034907000000000001]'; %Orientation du vecteur par rapport à l'axe z
% theta_T11 = [0 0 -0.069813]'; %Orientation du vecteur par rapport à l'axe z
% theta_T10 = [0 0 -0.087265999999999996]'; %Orientation du vecteur par rapport à l'axe z
% theta_T9  = [0 0 -0.069813]'; %Orientation du vecteur par rapport à l'axe z
% theta_T8  = [0 0 -0.087265999999999996]'; %Orientation du vecteur par rapport à l'axe z
% theta_T7  = [0 0 -0.10471999999999999]'; %Orientation du vecteur par rapport à l'axe z
% theta_T6  = [0 0 -0.10471999999999999]'; %Orientation du vecteur par rapport à l'axe z
% theta_T5  = [0 0 -0.10471999999999999]'; %Orientation du vecteur par rapport à l'axe z
% theta_T4  = [0 0 -0.10471999999999999]'; %Orientation du vecteur par rapport à l'axe z
% theta_T3  = [0 0 -0.069813]'; %Orientation du vecteur par rapport à l'axe z
% theta_T2  = [0 0 -0.052359999999999997]'; %Orientation du vecteur par rapport à l'axe z
% theta_T1  = [0 0 -0.017453]'; %Orientation du vecteur par rapport à l'axe z

% Anatomical landmarks in the case of rib cage modelling %%%%%%%%%%%%%%%%%%
% Sternum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sternum_topCT = M*k_new*[11.4366 -202.505 1050.14]'*0.001;
% T1_SternumJointNodeCT = M*k_new*[9.9997 -203.4495 1030.735]'*0.001; % Point entre les deux articulations sterno-costales au niveau de R1
% T1_SternumJointNode = T1_SternumJointNodeCT - T2_T1CT; % Centre de l'articulation entre les côtes R1 et le manubrium dans le repère T1 comme étant une liaison avec T1
% Sternum_T7JointNodeCT = M*k_new*[6.2591 -235.2 911.254]'*0.001; % Point entre les deux articulations sterno-costales au niveau de R7
% Sternum_T7JointNode = Sternum_T7JointNodeCT - T1_SternumJointNodeCT; % Centre articulaire entre les côtes R7 et le sternum dans le repère sternum

% Centre de masse du thorax avec la formule du barycentre %%%%%%%%%%%%%%%%%
% 
% CoM_ThoraxCT = (MassT1*CoM_T1CT+MassT2*CoM_T2CT+MassT3*CoM_T3CT+MassT4*CoM_T4CT+...
%                 MassT5*CoM_T5CT+MassT6*CoM_T6CT+MassT7*CoM_T7CT+MassT8*CoM_T8CT+...
%                 MassT9*CoM_T9CT+MassT10*CoM_T10CT+MassT11*CoM_T11CT+MassT12*CoM_T12CT)/MassThorax;
% 
% CoM_Thorax = CoM_ThoraxCT - L1_ThoraxCT;%Repère segment parent

%% Attachement tête %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Thorax_C7JointNodeCT = M*k_new*[13.2581 -161.816 1093.62]'*0.001; % Correspond au corps inférieur C7
%Thorax_C7JointNode = Thorax_C7JointNodeCT - T2_T1CT;
%C1_HeadJointNodeCT = M*k_new*[29.6668 -195.846 1197.5]'*0.001; %Correspond au sommet de la C2
%C1_HeadJointNode = C1_HeadJointNodeCT - T2_T1CT;
%NeckNode = [0 C1_HeadJointNode(2) 0]'; %On verticalise du fait de la position allongée et de la scoliose

%% Attachement Bras %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Modèle de Twente:
Thorax_RClavJointNodeCT = M*k_new*[-10.7812 -205.974 1051.98]'*0.001; %Point d'attache RClavicule dans le repère du CT-scan. Comme dans OpenSim
Thorax_scjJointRightNode_supine = Thorax_RClavJointNodeCT - T2_T1CT; % Point d'attache RClavicle dans le repère T1 sans prendre en compte l'orientation comme on a des contraintes cinématiques
x_component = (Thorax_scjJointRightNode_supine(1)^2+Thorax_scjJointRightNode_supine(2)^2)^(1/2); %Réajustement par rapport à la position du cadavre dans l'étude de twente (supine vers debout)
z_component = Thorax_scjJointRightNode_supine(3);
Thorax_scjJointRightNode = [x_component 0 z_component]';
Thorax_scjJointLeftNode = S*Thorax_scjJointRightNode;

% % Comme la position de l'articulation sterno-claviculaire est mal
% % représentée dans le modèle de Twente, on adapte cette partie du modèle de
% % Holzbaur :
% Thorax_RClavJointNodeCT = M*k_new*[-10.7812 -205.974 1051.98]'*0.001; %Point d'attache RClavicule dans le repère du CT-scan. Comme dans OpenSim
% % on a adapté l'articulation sterno-claviculaire du modèle de holzbaur :
% Thorax_RClavJointNode_repere_holzbaur = k*(([0.044784 0.052614 0.015413]'+[0.045522 0.051582 0.026099]')/2); %Point d'attache RClavicule dans le repère du CT-scan. Comme dans OpenSim
% T1T2_repere_holzbaur = k*(([0.013874 0.061512 0]'+[-0.004786 0.071156 0]')/2); %Point d'attache RClavicule dans le repère du CT-scan. Comme dans OpenSim
% Thorax_scjJointRightNode = Thorax_RClavJointNode_repere_holzbaur - T1T2_repere_holzbaur; % Point d'attache RClavicle dans le repère T1
% Thorax_scjJointLeftNode = S*Thorax_scjJointRightNode;

RScapula_HumerusJointNodeCT = M*k_new*[-130.311 -184.256 1118.01]'*0.001;%Tète superieure humérus
Thorax_ShoulderRightNode = RScapula_HumerusJointNodeCT - Thorax_RClavJointNodeCT; %Bras droit dans le repère RClavicle
Thorax_ShoulderLeftNode = S*Thorax_ShoulderRightNode;  %Bras gauche dans le repère LClavicle

% Ajout des marqueurs du cluster pour enlever les marqueurs de la colonne %

CLAV_CT = M*k_new*[11.2451 -206.483 1046.82]'*0.001;
Cl3  = CLAV_CT +  [0 -44.887 -22.18]'*0.001;
Cl4  = Cl3  +  [0  0  30]'*0.001;
Cl2  = Cl3  +  [0  -20  -40]'*0.001;
Cl1  = Cl2  +  [0  -70   0 ]'*0.001;
Cl6  = Cl1  +  [0   0    120]'*0.001;
Cl5  = Cl6  +  [0   80  -15]'*0.001;

% Partie adaptée du modèle de OpenSim que l'on retrouve aussi dans
% UpperTrunkClavicle pour les clavicules
Thorax_osim2antoine     = [k k Thorax_ShoulderRightNode(3)/0.17]; % scaling coef based on shoulder width
CoM_Thorax_osim         = Thorax_osim2antoine.*[-0.0591 -0.1486 0];
CoM_Clavicle_osim       = Thorax_osim2antoine.*[-0.011096 0.0063723 0.054168];
Scapula_CoM_osim        = Thorax_osim2antoine.*[-0.054694 -0.035032 -0.043734];
Clavicle_acJointNode_osim = Thorax_osim2antoine.*[-0.02924 0.02024 0.12005]-CoM_Clavicle_osim;
Scapula_acJointNode_osim     = Thorax_osim2antoine.*[-0.01357 0.00011 -0.01523]-Scapula_CoM_osim;
Clavicle2Scapula        = Scapula_acJointNode_osim - Clavicle_acJointNode_osim;


%% Thoracic region :

%% Definition of anatomical landmarks for the muscles and markers (wrt local RF on CoM)

T12_position_set = {
    % Node (repère local)
    % 'Thorax_Origin', [0 0 0]'-CoM_T12;...
    'T12_Origin',[0 0 0]'-CoM_T12;...
    'T12_T11Node', T12_T11JointNode-CoM_T12;...
    % 'BARYCENTRE', centre_geom_T12CT-CoM_T12CT;...
    % Markers (repère global)
    'T12', M*k_new*[-5.1556 -39.5187 833.379]'*0.001- CoM_T12CT;...
    'T12POSSUP', M*k_new*[-2.95 -80.73 858.31]'*0.001- CoM_T12CT;...
    };

T11_position_set = {
    % Node (repère local)
    'T11_Origin', [0 0 0]'-CoM_T11;...
    'T11_T10Node', T11_T10JointNode-CoM_T11;...
    % Markers (repère global)
    'T11', M*k_new*0.001*[-4.2833 -32.5733 869.766]' - CoM_T11CT;...
    'T11POSSUP', M*k_new*0.001*[-3.445 -76.365 887.435]' - CoM_T11CT;... 
    };

T10_position_set = {
    % Node (repère local)
    'T10_Origin', [0 0 0]'-CoM_T10;...
    'T10_T9Node', T10_T9JointNode-CoM_T10;...
    % Markers (repère global)
    'T10', M*k_new*[-7.3757 -34.1193 884.079]'*0.001 - CoM_T10CT;...
    'T10_visu', M*k_new*[-7.3757 -34.1193 884.079]'*0.001 - CoM_T10CT;...
    'T10POSSUP', M*k_new*[-5.325 -74.385 916.365]'*0.001 - CoM_T10CT;...
    };

T9_position_set = {
    % Node (repère local)
    'T9_Origin', [0 0 0]'-CoM_T9;...
    'T9_T8Node', T9_T8JointNode-CoM_T9;...
    % Markers (repère global)
    'T9', M*k_new*0.001*[-6.3262 -34.5678 899.045]' - CoM_T9CT;...
    'T9visu', M*k_new*0.001*[-6.3262 -34.5678 899.045]' - CoM_T9CT;...
    'T9POSSUP', M*k_new*0.001*[-6.08 -74.485 941.285]' - CoM_T9CT;... 
    };

T8_position_set = {
    % Node (repère local)
    'T8_Origin', [0 0 0]'-CoM_T8;...
    'T8_T7Node', T8_T7JointNode-CoM_T8;...
    % Markers (repère global)
    'T8', M*k_new*[-6.7407 -31.9621 926.954]'*0.001- CoM_T8CT;...
    'T8visu', M*k_new*[-6.7407 -31.9621 926.954]'*0.001- CoM_T8CT;...
    'T8POSSUP', M*k_new*[-7.225 -77.765 966.7]'*0.001- CoM_T8CT;...
    };

T7_position_set = {
    % Node (repère local)
    'T7_Origin', [0 0 0]'-CoM_T7;...
    'T7_T6Node', T7_T6JointNode-CoM_T7;...
    % Markers (repère global)
    'T7', M*k_new*0.001*[-5.4078 -37.3452 941.961]' - CoM_T7CT;...
    'STRN', M*k_new*[3.4527 -236.945 876.125]'*0.001 - CoM_T7CT;... %Xyphoid process
    'STRNvisu', M*k_new*[3.4527 -236.945 876.125]'*0.001 - CoM_T7CT;... %Xyphoid process
    'T7POSSUP', M*k_new*0.001*[-7.9 -82.46 989.585]' - CoM_T7CT;...
    };

T6_position_set = {
    % Node (repère local)
    'T6_Origin', [0 0 0]'-CoM_T6;...
    'T6_T5Node', T6_T5JointNode-CoM_T6;...
    % Markers (repère global)
    'T6', M*k_new*0.001*[-6.3339 -40.6978 971.015]' - CoM_T6CT;...
    'T6POSSUP', M*k_new*0.001*[-9.41 -90.695 1012.37]' - CoM_T6CT;...
    };

T5_position_set = {
    % Node (repère local)
    'T5_Origin', [0 0 0]'-CoM_T5;...
    'T5_T4Node', T5_T4JointNode-CoM_T5;...
    % Markers (repère global)
    'T5', M*k_new*0.001*[-6.8294 -48.5802 1000.22]' - CoM_T5CT;...
    'T5POSSUP', M*k_new*0.001*[-8.335 -97.075 1032.205]' - CoM_T5CT;...
    };

T4_position_set = {
    % Node (repère local)
    'T4_Origin', [0 0 0]'-CoM_T4;...
    'T4_T3Node', T4_T3JointNode-CoM_T4;...
    % Markers (repère global)
    'T4', M*k_new*0.001*[-10.1477 -60.1759 1035.06]' - CoM_T4CT;...
    'T4visu', M*k_new*0.001*[-10.1477 -60.1759 1035.06]' - CoM_T4CT;...
    'T4POSSUP', M*k_new*0.001*[-2 -105.12 1049.495]' - CoM_T4CT;...
    };

T3_position_set = {
    % Node (repère local)
    'T3_Origin', [0 0 0]'-CoM_T3;...
    'T3_T2Node', T3_T2JointNode-CoM_T3;...
    % Markers (repère global)
    'T3', M*k_new*0.001*[-1.0426 -73.0513 1057.13]' - CoM_T3CT;...
    'T3POSSUP', M*k_new*0.001*[0.96 -118.01 1064.805]' - CoM_T3CT;...
    };

T2_position_set = {
    % Node (repère local)
    'T2_Origin', [0 0 0]'-CoM_T2;...
    'T2_T1Node', T2_T1JointNode-CoM_T2;...
    % Markers (repère global)
    'T2', M*k_new*0.001*[-3.586 -87.832 1082.41]' - CoM_T2CT;...
    'T2POSSUP', M*k_new*0.001*[-3.335 -131.245 1079.735]' - CoM_T2CT;...
    };

T1_position_set = {
    % Node (repère local)
    'T1_Origin', [0 0 0]'-CoM_T1;...
    %'CervicalNode', Thorax_C7JointNode-CoM_T1;...
    'T1_C7Node', T1_C7JointNode-CoM_T1;...
    'Thorax_scjJointRightNode', Thorax_scjJointRightNode-CoM_T1; ...
    'Thorax_scjJointLeftNode', Thorax_scjJointLeftNode-CoM_T1; ...  
    'T1', M*k_new*0.001*[7.5311 -105.019 1102.14]' - CoM_T1CT;...
    % 'CLAV', M*k_new*[11.2451 -206.483 1046.82]'*0.001 - CoM_T1CT;...
    'CLAV', [norm(M*k_new*[11.2451 -206.483 1046.82]'*0.001 - CoM_T1CT) 0 0]';... réajustement par rapport à la position du cadavre dans l'étude de twente (supination vers debout)
    'CLAVvisu', [norm(M*k_new*[11.2451 -206.483 1046.82]'*0.001 - CoM_T1CT) 0 0]';...
    'T1POSSUP', M*k_new*0.001*[6.8652 -143.845 1094.21]' - CoM_T1CT;...
    };

RClavicle_position_set = {
    % Node (repère local)
    'RClavicle_Origin', [0 0 0]';...
    'Thorax_ShoulderRightNode', Thorax_ShoulderRightNode; ...
    % Markers
    'RSHO',M*k_new*[-159.888 -182.649 1104.84]'*0.001 - Thorax_RClavJointNodeCT;... %Tête humérale avant
    };

LClavicle_position_set = {
    % Node (repère local)
    'LClavicle_Origin', [0 0 0]';...
    'Thorax_ShoulderLeftNode', Thorax_ShoulderLeftNode; ...
    % Markers
    'LSHO',S*(M*k_new*[-159.888 -182.649 1104.84]'*0.001 - Thorax_RClavJointNodeCT);... %Tête humérale avant
    };

%
C7_position_set = {
    % Node (repère local)
    'C7_Origin', [0 0 0]'-CoM_C7;...
    'C7_C6Node', C7_C6JointNode-CoM_C7;...
    % Markers (repère global)
    'C7', M*k_new*0.001*[11.3472 -115.937 1117.14]'- CoM_C7CT;... %Scoliosis adjustment
    'C7POSSUP', M*k_new*0.001*[10.83 -152.8 1107]' - CoM_C7CT;...
    };
    
C6_position_set = {
    % Node (repère local)
    'C6_Origin', [0 0 0]'-CoM_C6;... %__Coordonnées de l'origine du repère global exprimé dans le repère local de C6
    'C6_C5Node', C6_C5JointNode-CoM_C6;...
    % Markers (repère global)
    'C6', M*k_new*0.001*[20.8273 -139.225 1130.94]' - CoM_C6CT;...
    'C6POSSUP', M*k_new*0.001*[20.24 -169.75 1118.66]' - CoM_C6CT;...
    };

C5_position_set = {
    % Node (repère local)
    'C5_Origin', [0 0 0]'-CoM_C5;...
    'C5_C4Node', C5_C4JointNode-CoM_C5;...
    % Markers (repère global)
    'C5', M*k_new*0.001*[21.0193 -153.035 1137.8]' - CoM_C5CT;...
    'C5POSSUP', M*k_new*0.001*[20.305 -179.02 1130.8]' - CoM_C5CT;...
    };

C4_position_set = {
    % Node (repère local)
    'C4_Origin', [0 0 0]'-CoM_C4;...
    'C4_C3Node', C4_C3JointNode-CoM_C4;...
    %'NeckNode', NeckNode - CoM_C4;...
    % Markers (repère global)
    'C4', M*k_new*0.001*[24.89 -150.15 1146.94]' - CoM_C4CT;...
    'C4POSSUP', M*k_new*0.001*[25.125 -185.295 1143.96]' - CoM_C4CT;...
    };

%
C3_position_set = {
    % Node (repère local)
    'C3_Origin', [0 0 0]'-CoM_C3;...
    'C3_C2Node',C3_C2JointNode-CoM_C3;...
    %'NeckNode', NeckNode - CoM_C3;...
    % Markers (repère global)
    'C3', M*k_new*0.001*[28.2856 -155.289 1158.88]' - CoM_C3CT;...
    'C3POSSUP', M*k_new*0.001*[28.315 -187.755 1159.215]' - CoM_C3CT;...
    };

%
C2_position_set = {
    % Node (repère local)
    'C2_Origin', [0 0 0]'-CoM_C2;...
    'C2_C1Node',C2_C1JointNode-CoM_C2;...
    %'NeckNode', NeckNode-CoM_C2;...
    % Markers (repère global)
    'C2', M*k_new*0.001*[28.8443 -154.133 1169.28]' - CoM_C2CT;...
    'C2POSSUP', M*k_new*0.001*[29.4915 -187.458 1178.79]' - CoM_C2CT;...
    };

%
C1_position_set = {
    % Node (repère local)
    'C1_Origin', [0 0 0]'-CoM_C1;...
    'NeckNode', C1_HeadJointNode - CoM_C1;...
    % Markers (repère global)
    'C1', M*k_new*0.001*[28.2751 -163.984 1185.61]' - CoM_C1CT;...
    'C1POSSUP', M*k_new*0.001*[29.50535 -190.725 1185.61]' - CoM_C1CT;...
    };
%}


%%                     Scaling inertial parameters
% From "Segmental Inertial Parameters of the Human Trunk as Determined from
% Computed Tomography", Pearsall et al. (1996)

%Each vertebrae Inertia matrix at the CoM :

IT1 = [20.2    0    0 ;...
          0 67.0    0 ;...
          0    0 87.2]/10000;
IT2 = [23.5    0    0 ;...
          0 67.3    0 ;...
          0    0 90.7]/10000;
IT3 = [31.4    0    0 ;...
          0 83.7    0 ;...
          0    0 115.1]/10000;
IT4 = [33.5    0    0 ;...
          0 83.0    0 ;...
          0    0 116.5]/10000;
IT5 = [35.1    0    0 ;...
          0 80.2    0 ;...
          0    0 115.3]/10000;
IT6 = [38.5    0    0 ;...
          0 78.0    0 ;...
          0    0 116.5]/10000;
IT7 = [40.8    0    0 ;...
          0 74.3    0 ;...
          0    0 115.2]/10000;
IT8 = [44.0    0    0 ;...
          0 72.5    0 ;...
          0    0 116.5]/10000;
IT9 = [46.6    0    0 ;...
          0 71.8    0 ;...
          0    0 118.4]/10000;
IT10 = [60.8    0    0 ;...
          0 89.0    0 ;...
          0    0 149.8]/10000;
IT11 = [61.5    0    0 ;...
          0 90.5    0 ;...
          0    0 152.0]/10000;
IT12 = [60.3    0    0 ;...
          0 110.5    0 ;...
          0    0 180.8]/10000;


Shoulder_Mass = Mass.Clavicle_Mass + Mass.Scapula_Mass;

%% ["Adjustments to McConville et al. and Young et al. body segment inertial parameters"] R. Dumas
% ------------------------- Thorax ----------------------------------------
% [I_Thorax]=rgyration2inertia([27 25 28 18 2 4*1i], UpperTrunk_Mass, [0 0 0], Length_Thorax);

% Generic Inertia extraced from (Klein Breteler et al. 1999)
Clavicle_Mass_generic=0.156;
I_Clavicle_generic=[0.00024259 0.00025526 0.00004442 -0.00001898 -0.00006994 0.00005371];
I_Clavicle=(norm(Thorax_osim2antoine)^2*Mass.Clavicle_Mass/Clavicle_Mass_generic).*I_Clavicle_generic;
Scapula_Mass_generic=0.70396;
I_Scapula_generic=[0.0012429 0.0011504 0.0013651 0.0004494 0.00040922 0.0002411];
I_Scapula=(norm(Thorax_osim2antoine)^2*Mass.Scapula_Mass/Scapula_Mass_generic)*I_Scapula_generic;
I_Shoulder=I_Clavicle+I_Scapula;

%% Equations for kinematic constrains

% Explication :
% Les 3 angles (FE/LB/AR) au niveau de L5S1 sont libres mais doivent rester
% dans les limites physiologiques. De même pour les angles au niveau de
% T12/L1. Les angles de FE et de AR sont libres au niveau de T8/T9 et en LB
% au niveau de T4/T5.
% On a dans le cas de la FE (resp. AR) :
% Dans ce cas (9ddl) :
% ThetaS1L1 = ThetaL5S1/ContribL5S1(%) donc les 
% Theta,i = (Contrib,i/ContribL5S1)*ThetaL5S1
% ThetaL1T9 = ThetaT12L1/ContribT12L1(%) donc les 
% Theta,j = (Contrib,j/ContribT12L1)*ThetaT12L1
% ThetaT9T1 = ThetaT8T9/ContribT8T9(%) donc les
% Theta,k = (Contrib,k/ContribT8T9)*ThetaT8T9
% Pour le cas de la LB
% ThetaS1L1 = ThetaL5S1/ContribL5S1(%) donc les 
% Theta,i = (Contrib,i/ContribL5S1)*ThetaL5S1
% ThetaL1T5 = ThetaT12L1/ContribT12L1(%) donc les 
% Theta,j = (Contrib,j/ContribT12L1)*ThetaT12L1
% ThetaT5T1 = ThetaT4T5/ContribT4T5(%) donc les
% Theta,k = (Contrib,k/ContribT4T5)*ThetaT4T5

L5S1_FE = 0.15; L5S1_AR = 0.04;

T12L1_FE  = 0.04/L5S1_FE;  T12L1_LB  = 0.11 ;           T12L1_AR  = 0.02/L5S1_AR;
T11T12_FE = 0.05/L5S1_FE;  T11T12_LB = 0.12/T12L1_LB ;  T11T12_AR = 0.04/L5S1_AR;
T10T11_FE = 0.05/L5S1_FE;  T10T11_LB = 0.10/T12L1_LB ;  T10T11_AR = 0.06/L5S1_AR;
T9T10_FE  = 0.03/L5S1_FE;  T9T10_LB  = 0.08/T12L1_LB ;  T9T10_AR  = 0.07/L5S1_AR;
T8T9_FE   = 0.13        ;  T8T9_LB   = 0.07/T12L1_LB ;  T8T9_AR   = 0.09/L5S1_AR;
T7T8_FE   = 0.12/T8T9_FE;  T7T8_LB   = 0.09/T12L1_LB ;  T7T8_AR   = 0.09/L5S1_AR;
T6T7_FE   = 0.10/T8T9_FE;  T6T7_LB   = 0.07/T12L1_LB ;  T6T7_AR   = 0.08/L5S1_AR;
T5T6_FE   = 0.11/T8T9_FE;  T5T6_LB   = 0.06/T12L1_LB ;  T5T6_AR   = 0.08/L5S1_AR;
T4T5_FE   = 0.06/T8T9_FE;  T4T5_LB   = 0.06/T12L1_LB ;  T4T5_AR   = 0.07/L5S1_AR;
T3T4_FE   = 0.11/T8T9_FE;  T3T4_LB   = 0.08/T12L1_LB ;  T3T4_AR   = 0.07/L5S1_AR;
T2T3_FE   = 0.17/T8T9_FE;  T2T3_LB   = 0.07/T12L1_LB ;  T2T3_AR   = 0.07/L5S1_AR;
T1T2_FE   = 0.20/T8T9_FE;  T1T2_LB   = 0.08/T12L1_LB ;  T1T2_AR   = 0.07/L5S1_AR;


C7T1_FE = 0.077;    C7T1_LB = 1;    C7T1_AR = 1;

C6C7_FE = 0.145/C7T1_FE;    C6C7_LB = 1/C7T1_LB;    C6C7_AR = 1/C7T1_AR;
C5C6_FE = 0.181/C7T1_FE;    C5C6_LB = 1/C7T1_LB;    C5C6_AR = 1/C7T1_AR;
C4C5_FE = 0.179/C7T1_FE;    C4C5_LB = 1/C7T1_LB;    C4C5_AR = 1/C7T1_AR;
C3C4_FE = 0.157/C7T1_FE;    C3C4_LB = 1/C7T1_LB;    C3C4_AR = 1/C7T1_AR;
C2C3_FE = 0.117/C7T1_FE;    C2C3_LB = 1/C7T1_LB;    C2C3_AR = 1/C7T1_AR;
C1C2_FE = 0.144/C7T1_FE;    C1C2_LB = 1/C7T1_LB;    C1C2_AR = 1/C7T1_AR;
%C0C1_FE = 1/C7T1_FE;    C0C1_LB = 1/C7T1_LB;    C0C1_AR = 1/C7T1_AR;

                    %% %% "Human_model" structure generation

num_solid = 0;
%% UpperTrunk

    %T12_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T12_J2;                   
    OsteoArticularModel(incr_solid).mother=s_mother;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= (4.8-2.9176)*(2*pi/360); % From White and Punjabi with René Louis' proportion "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= (-7.2-2.9176)*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-15,T12L1_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=pos_attachment_pt;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T12_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T12;                   
    OsteoArticularModel(incr_solid).mother=s_T12_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_inf= 1.5*(-4-0.1825)*(2*pi/360); % From White and Punjabi with 1.5 coef.
    OsteoArticularModel(incr_solid).limit_sup= 1.5*(4-0.1825)*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= inf;
    % OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T12
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T11_J1;
    OsteoArticularModel(incr_solid).mother=s_T12_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= (-5-0.1081)*(2*pi/360); % From René Louis
    % OsteoArticularModel(incr_solid).limit_sup=  (5-0.1081)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-15,T12L1_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T12;
    OsteoArticularModel(incr_solid).m=MassT12;               
    OsteoArticularModel(incr_solid).I=IT12;
    OsteoArticularModel(incr_solid).anat_position=T12_position_set;
    OsteoArticularModel(incr_solid).L={'Thorax_Origin';'T12_T11JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T12.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T12(3);

    %T11_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T11_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T12;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-18,T11T12_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T12_T11JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T11_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T11;                   
    OsteoArticularModel(incr_solid).mother=s_T11_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-3,T11T12_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T11
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T10_J1;
    OsteoArticularModel(incr_solid).mother=s_T11_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-18,T11T12_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T11;
    OsteoArticularModel(incr_solid).m=MassT11;               
    OsteoArticularModel(incr_solid).I=IT11;
    OsteoArticularModel(incr_solid).anat_position=T11_position_set;
    OsteoArticularModel(incr_solid).L={'T12_T11JointNode';'T11_T10JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T11.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T11(3);

    %T10_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T10_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T11;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-21,T10T11_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T11_T10JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T10_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T10;                   
    OsteoArticularModel(incr_solid).mother=s_T10_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-6,T10T11_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T10
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T9_J1;
    OsteoArticularModel(incr_solid).mother=s_T10_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-3*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-21,T10T11_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T10;
    OsteoArticularModel(incr_solid).m=MassT10;               
    OsteoArticularModel(incr_solid).I=IT10;
    OsteoArticularModel(incr_solid).anat_position=T10_position_set;
    OsteoArticularModel(incr_solid).L={'T11_T10JointNode';'T10_T9JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T10.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T10(3);

    %T9_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T9_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T10;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-24,T9T10_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T10_T9JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T9_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T9;                   
    OsteoArticularModel(incr_solid).mother=s_T9_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-9,T9T10_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T9
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T8_J1;
    OsteoArticularModel(incr_solid).mother=s_T9_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-24,T9T10_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T9;
    OsteoArticularModel(incr_solid).m=MassT9;               
    OsteoArticularModel(incr_solid).I=IT9;
    OsteoArticularModel(incr_solid).anat_position=T9_position_set;
    OsteoArticularModel(incr_solid).L={'T10_T9JointNode';'T9_T8JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T9.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T9(3);

    %T8_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T8_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T9;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_sup=  3*2.4*(2*pi/360); % From White and punjabi with René Louis' "Anatomy of the spine" proportion and coef 2*1.5.
    OsteoArticularModel(incr_solid).limit_inf=  3*(-3.6)*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= inf;
    % OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T9_T8JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T8_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T8;                   
    OsteoArticularModel(incr_solid).mother=s_T8_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-12, T8T9_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T8
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T7_J1;  
    OsteoArticularModel(incr_solid).mother=s_T8_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= (-3.5+0.4417)*(2*pi/360); %from white and punjabi
    % OsteoArticularModel(incr_solid).limit_sup=  (3.5+0.4417)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-27,T8T9_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T8;
    OsteoArticularModel(incr_solid).m=MassT8;               
    OsteoArticularModel(incr_solid).I=IT8;
    OsteoArticularModel(incr_solid).anat_position=T8_position_set;
    OsteoArticularModel(incr_solid).L={'T9_T8JointNode';'T8_T7JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T8.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T8(3);

    %T7_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T7_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T8;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-3, T7T8_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T8_T7JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T7_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T7;                   
    OsteoArticularModel(incr_solid).mother=s_T7_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-15, T7T8_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T7
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T6_J1; 
    OsteoArticularModel(incr_solid).mother=s_T7_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-30, T7T8_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T7;
    OsteoArticularModel(incr_solid).m=MassT7;               
    OsteoArticularModel(incr_solid).I=IT7;
    OsteoArticularModel(incr_solid).anat_position=T7_position_set;
    OsteoArticularModel(incr_solid).L={'T8_T7JointNode';'T7_T6JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T7.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T7(3);

    %T6_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T6_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T7;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-6,T6T7_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T7_T6JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T6_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T6;                   
    OsteoArticularModel(incr_solid).mother=s_T6_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-18,T6T7_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T6
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T5_J1;                   
    OsteoArticularModel(incr_solid).mother=s_T6_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-4*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 4*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-33,T6T7_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T6;
    OsteoArticularModel(incr_solid).m=MassT6;               
    OsteoArticularModel(incr_solid).I=IT6;
    OsteoArticularModel(incr_solid).anat_position=T6_position_set;
    OsteoArticularModel(incr_solid).L={'T7_T6JointNode';'T6_T5JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T6.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T6(3);

    %T5_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T5_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T6;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-9, T5T6_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T6_T5JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T5_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T5;                   
    OsteoArticularModel(incr_solid).mother=s_T5_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-21, T5T6_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T5
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T4_J1;                   
    OsteoArticularModel(incr_solid).mother=s_T5_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-3*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-36, T5T6_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T5;
    OsteoArticularModel(incr_solid).m=MassT5;               
    OsteoArticularModel(incr_solid).I=IT5;
    OsteoArticularModel(incr_solid).anat_position=T5_position_set;
    OsteoArticularModel(incr_solid).L={'T6_T5JointNode';'T5_T4JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T5.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T5(3);

    %T4_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T4_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T5;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-12, T4T5_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T5_T4JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T4_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T4;                   
    OsteoArticularModel(incr_solid).mother=s_T4_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2.5*(-3+0.0076)*(2*pi/360); %Coef. 2.5 * limites de White et Punjabi
    % OsteoArticularModel(incr_solid).limit_sup= 2.5*(3+0.0076)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-24, T4T5_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T4
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T3_J1;
    OsteoArticularModel(incr_solid).mother=s_T4_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-3*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-39, T4T5_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T4;
    OsteoArticularModel(incr_solid).m=MassT4;               
    OsteoArticularModel(incr_solid).I=IT4;
    OsteoArticularModel(incr_solid).anat_position=T4_position_set;
    OsteoArticularModel(incr_solid).L={'T5_T4JointNode';'T4_T3JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T4.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T4(3);

    %T3_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T3_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T4;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-15, T3T4_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T4_T3JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);
     
    %T3_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T3;                   
    OsteoArticularModel(incr_solid).mother=s_T3_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-27, T3T4_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);
     
    %T3
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T2_J1;                   
    OsteoArticularModel(incr_solid).mother=s_T3_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-3*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-42, T3T4_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T3;
    OsteoArticularModel(incr_solid).m=MassT3;               
    OsteoArticularModel(incr_solid).I=IT3;
    OsteoArticularModel(incr_solid).anat_position=T3_position_set;
    OsteoArticularModel(incr_solid).L={'T4_T3JointNode';'T3_T2JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T3.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T3(3);
     
    %T2_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T2_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T3;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-18, T2T3_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T3_T2JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %T2_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T2;                   
    OsteoArticularModel(incr_solid).mother=s_T2_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-30, T2T3_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);
     
    %T2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T1_J1;                   
    OsteoArticularModel(incr_solid).mother=s_T2_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-3*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-45, T2T3_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T2;
    OsteoArticularModel(incr_solid).m=MassT2;               
    OsteoArticularModel(incr_solid).I=IT2;
    OsteoArticularModel(incr_solid).anat_position=T2_position_set;
    OsteoArticularModel(incr_solid).L={'T3_T2JointNode';'T2_T1JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T2.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T2(3);
     
    %T1_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T1_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T2;           
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_sup= 2 *  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf= 2 *  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-21, T1T2_FE]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T2_T1JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);
     
    %T1_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T1;                   
    OsteoArticularModel(incr_solid).mother=s_T1_J1;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-33, T1T2_LB]';
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

     
    %T1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_RClavicle_J1;                  
    OsteoArticularModel(incr_solid).mother=s_T1_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    % OsteoArticularModel(incr_solid).limit_inf= 2 *-3*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2 * 3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_inf= -inf;
    OsteoArticularModel(incr_solid).limit_sup= inf;
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-48, T1T2_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T1;
    OsteoArticularModel(incr_solid).m=MassT1;               
    OsteoArticularModel(incr_solid).I=IT1;
    OsteoArticularModel(incr_solid).anat_position=T1_position_set;
    OsteoArticularModel(incr_solid).L={'T2_T1JointNode';'CervicalNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/T1.mat'];
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_T1(3);

    %% Rclavicle
    % RClavicle_J1
    num_solid=num_solid+1;        % number of the solid ...
    name=list_solid{num_solid}; % solid name
    eval(['incr_solid=s_' name ';'])  % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;               % solid name
    OsteoArticularModel(incr_solid).sister=s_LClavicle_J1;              
    OsteoArticularModel(incr_solid).child=s_RClavicle_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T1;           
    OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).limit_inf=-pi/6;
    OsteoArticularModel(incr_solid).limit_sup=pi/6;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).b=Thorax_scjJointRightNode;  
    OsteoArticularModel(incr_solid).I=zeros(3,3);
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).comment='Right Clavicle Protraction(+)/Retraction(-)';
    OsteoArticularModel(incr_solid).FunctionalAngle='Right Clavicle Protraction(+)/Retraction(-)';

    % RClavicle_J2
    num_solid=num_solid+1;        % number of the solid ...
    name=list_solid{num_solid}; % solid name
    eval(['incr_solid=s_' name ';'])  % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;               % solid name
    OsteoArticularModel(incr_solid).sister=0;              
    OsteoArticularModel(incr_solid).child=s_RClavicle;                   
    OsteoArticularModel(incr_solid).mother=s_RClavicle_J1;           
    OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).limit_inf=-pi/4;
    OsteoArticularModel(incr_solid).limit_sup=0.2;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).b=[0 0 0]';  
    OsteoArticularModel(incr_solid).I=zeros(3,3);
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).comment='Right Clavicle Depression(+)/Elevation(-)';
    OsteoArticularModel(incr_solid).FunctionalAngle='Right Clavicle Depression(+)/Elevation(-)';

    % RClavicle
    num_solid=num_solid+1;        % number of the solid ...
    name=list_solid{num_solid}; % solid name
    eval(['incr_solid=s_' name ';'])  % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;               % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=0;                   
    OsteoArticularModel(incr_solid).mother=s_RClavicle_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 1]';    
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).limit_inf=-pi/4;
    OsteoArticularModel(incr_solid).limit_sup=pi/4;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).m=Shoulder_Mass;                 
    OsteoArticularModel(incr_solid).b=[0 0 0]';  
    OsteoArticularModel(incr_solid).I=[I_Shoulder(1) I_Shoulder(4) I_Shoulder(5); I_Shoulder(4) I_Shoulder(2) I_Shoulder(6); I_Shoulder(5) I_Shoulder(6) I_Shoulder(3)];
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).anat_position=RClavicle_position_set;
    OsteoArticularModel(incr_solid).visual_file = ['Holzbaur/clavicle_r.mat'];
    OsteoArticularModel(incr_solid).comment='Right Clavicle Axial rotation Forward(-)/Backward(+)';
    OsteoArticularModel(incr_solid).FunctionalAngle='Right Clavicle Axial rotation Forward(-)/Backward(+)';
    
    % Wrapping 1
    OsteoArticularModel(incr_solid).wrap(1).name='WrapRThoraxGH';
    OsteoArticularModel(incr_solid).wrap(1).anat_position='WrapRThoraxGH';
    OsteoArticularModel(incr_solid).wrap(1).type='S'; % C: Cylinder or S: Sphere
    OsteoArticularModel(incr_solid).wrap(1).R=k*0.035;
    OsteoArticularModel(incr_solid).wrap(1).orientation=eye(3);
    OsteoArticularModel(incr_solid).wrap(1).location=Thorax_osim2antoine'.*([-0.0058 -0.0378 0.0096]')-Scapula_CoM_osim'-Clavicle2Scapula';
    OsteoArticularModel(incr_solid).wrap(1).h=0;
    OsteoArticularModel(incr_solid).wrap(1).num_solid=incr_solid;
    
 
    %% Lclavicle
    % LClavicle_J1
    num_solid=num_solid+1;        % number of the solid ...
    name=list_solid{num_solid}; % solid name
    eval(['incr_solid=s_' name ';'])  % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;               % solid name
    %OsteoArticularModel(incr_solid).sister=0;
    OsteoArticularModel(incr_solid).sister=s_C7_J1;
    OsteoArticularModel(incr_solid).child=s_LClavicle_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T1;           
    OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).limit_inf=-pi/6;
    OsteoArticularModel(incr_solid).limit_sup=pi/6;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).b=Thorax_scjJointLeftNode;  
    OsteoArticularModel(incr_solid).I=zeros(3,3);
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).comment='Left Clavicle Protraction(-)/Retraction(+)';
    OsteoArticularModel(incr_solid).FunctionalAngle='Left Clavicle Protraction(+)/Retraction(-)';
    
    % LClavicle_J2
    num_solid=num_solid+1;        % number of the solid ...
    name=list_solid{num_solid}; % solid name
    eval(['incr_solid=s_' name ';'])  % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;               % solid name
    OsteoArticularModel(incr_solid).sister=0;              
    OsteoArticularModel(incr_solid).child=s_LClavicle;                   
    OsteoArticularModel(incr_solid).mother=s_LClavicle_J1;           
    OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).limit_inf=-pi/4;
    OsteoArticularModel(incr_solid).limit_sup=0.2;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).b=[0 0 0]';  
    OsteoArticularModel(incr_solid).I=zeros(3,3);
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).comment='Left Clavicle Depression(-)/Elevation(+)';
    OsteoArticularModel(incr_solid).FunctionalAngle='Left Clavicle Depression(-)/Elevation(+)';

    % LClavicle
    num_solid=num_solid+1;        % number of the solid ...
    name=list_solid{num_solid}; % solid name
    eval(['incr_solid=s_' name ';'])  % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;               % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=0;                   
    OsteoArticularModel(incr_solid).mother=s_LClavicle_J2;           
    OsteoArticularModel(incr_solid).a=[0 0 1]';    
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).limit_inf=-pi/4;
    OsteoArticularModel(incr_solid).limit_sup=pi/4;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).m=Shoulder_Mass;                 
    OsteoArticularModel(incr_solid).b=[0 0 0]';  
    OsteoArticularModel(incr_solid).I=[1 1 -1; 1 1 -1; -1 -1 1].*[I_Shoulder(1) I_Shoulder(4) I_Shoulder(5); I_Shoulder(4) I_Shoulder(2) I_Shoulder(6); I_Shoulder(5) I_Shoulder(6) I_Shoulder(3)];
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).anat_position=LClavicle_position_set;
    OsteoArticularModel(incr_solid).visual_file = ['Holzbaur/clavicle_l.mat'];
    OsteoArticularModel(incr_solid).comment='Left Clavicle Axial rotation Forward(-)/Backward(+)';
    OsteoArticularModel(incr_solid).FunctionalAngle='Left Clavicle Axial rotation Forward(-)/Backward(+)';
  
    % Wrapping 2
    OsteoArticularModel(incr_solid).wrap(1).name='WrapLThoraxGH';
    OsteoArticularModel(incr_solid).wrap(1).anat_position='WrapLThoraxGH';
    OsteoArticularModel(incr_solid).wrap(1).type='S'; % C: Cylinder or S: Sphere
    OsteoArticularModel(incr_solid).wrap(1).R=k*0.035;
    OsteoArticularModel(incr_solid).wrap(1).orientation=eye(3);
    OsteoArticularModel(incr_solid).wrap(1).location=[1 0 0; 0 1 0; 0 0 -1]*(Thorax_osim2antoine'.*([-0.0058 -0.0378 0.0096]')-Scapula_CoM_osim'-Clavicle2Scapula');
    OsteoArticularModel(incr_solid).wrap(1).h=0;
    OsteoArticularModel(incr_solid).wrap(1).num_solid=incr_solid;

    %% Cervical vertebra
    %
    %C7_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C7_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';
    OsteoArticularModel(incr_solid).joint=1;


    OsteoArticularModel(incr_solid).limit_sup=  (25/4)*(2*pi/360); % From René Louis "Anatomy of the spine"
    OsteoArticularModel(incr_solid).limit_inf=  -(15/4)*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=T1_C7JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C7_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C7;                   
    OsteoArticularModel(incr_solid).mother=s_C7_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1;  


    OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C7
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C6_J1;
    OsteoArticularModel(incr_solid).mother=s_C7_J2;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_C7(3);
    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);

    
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C7;
    OsteoArticularModel(incr_solid).anat_position=C7_position_set;
    OsteoArticularModel(incr_solid).L={'T1_C7JointNode';'C7_C6JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];



    %C6_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C6_J2;                   
    OsteoArticularModel(incr_solid).mother=s_C7;
    OsteoArticularModel(incr_solid).a=[0 0 1]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-3, C6C7_FE]';

    % OsteoArticularModel(incr_solid).limit_sup=  (45/4)*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -(27/4)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=C7_C6JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C6_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C6;                   
    OsteoArticularModel(incr_solid).mother=s_C6_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 
    
    
    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-3, C6C7_LB]';
    
    
    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  0.00001;
    OsteoArticularModel(incr_solid).limit_inf=  -0.00001;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C6
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C5_J1;
    OsteoArticularModel(incr_solid).mother=s_C6_J2;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1;


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-3, C6C7_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_C6(3);


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  0.00001;
    OsteoArticularModel(incr_solid).limit_inf=  -0.00001;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C6;
    %OsteoArticularModel(incr_solid).m=MassT11;               
    %OsteoArticularModel(incr_solid).I=IT11;
    OsteoArticularModel(incr_solid).anat_position=C6_position_set;
    OsteoArticularModel(incr_solid).L={'C7_C6JointNode';'C6_C5JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];



    %C5_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C5_J2;                   
    OsteoArticularModel(incr_solid).mother=s_C6;
    OsteoArticularModel(incr_solid).a=[0 0 1]';
    OsteoArticularModel(incr_solid).joint=1;  


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-6, C5C6_FE]';


    % OsteoArticularModel(incr_solid).limit_sup=  (55/4)*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -(33/4)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=C6_C5JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C5_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C5;                   
    OsteoArticularModel(incr_solid).mother=s_C5_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-6, C5C6_LB]';


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C5
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C4_J1;
    OsteoArticularModel(incr_solid).mother=s_C5_J2;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1;


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-6, C5C6_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_C5(3);


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C5;
    %OsteoArticularModel(incr_solid).m=MassT11;               
    %OsteoArticularModel(incr_solid).I=IT11;
    OsteoArticularModel(incr_solid).anat_position=C5_position_set;
    OsteoArticularModel(incr_solid).L={'C6_C5JointNode';'C5_C4JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];



    %C4_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C4_J2;                   
    OsteoArticularModel(incr_solid).mother=s_C5;
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-9, C4C5_FE]';


    % OsteoArticularModel(incr_solid).limit_sup=  (25/2)*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -(15/2)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=C5_C4JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C4_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C4;                   
    OsteoArticularModel(incr_solid).mother=s_C4_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-9, C4C5_LB]';


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C4
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).mother=s_C4_J2;
    OsteoArticularModel(incr_solid).child=s_C3_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-9, C4C5_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_C4(3);


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C4;
    OsteoArticularModel(incr_solid).anat_position=C4_position_set;
    OsteoArticularModel(incr_solid).L={'C5_C4JointNode';'C4_C3JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];


%
    %C3_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C3_J2;                   
    OsteoArticularModel(incr_solid).mother=s_C4;
    OsteoArticularModel(incr_solid).a=[0 0 1]';
    OsteoArticularModel(incr_solid).joint=1;


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-12, C3C4_FE]';


    % OsteoArticularModel(incr_solid).limit_sup=  (75/8)*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -(45/8)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=C4_C3JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C3_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C3;                   
    OsteoArticularModel(incr_solid).mother=s_C3_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1;


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-12, C3C4_LB]';


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C3
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).mother=s_C3_J2;
    OsteoArticularModel(incr_solid).child=s_C2_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1;  


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-12, C3C4_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_C3(3);


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C3;
    OsteoArticularModel(incr_solid).anat_position=C3_position_set;
    OsteoArticularModel(incr_solid).L={'C4_C3JointNode';'C3_C2JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];
    

%
    %C2_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C2_J2;                   
    OsteoArticularModel(incr_solid).mother=s_C3;
    OsteoArticularModel(incr_solid).a=[0 0 1]';
    OsteoArticularModel(incr_solid).joint=1;  


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-15, C2C3_FE]';


    % OsteoArticularModel(incr_solid).limit_sup=  (75/8)*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -(45/8)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=C3_C2JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C2_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C2;                   
    OsteoArticularModel(incr_solid).mother=s_C2_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-15, C2C3_LB]';


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).mother=s_C2_J2;
    OsteoArticularModel(incr_solid).child=s_C1_J1; 
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-15, C2C3_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_C2(3);


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C2;
    OsteoArticularModel(incr_solid).anat_position=C2_position_set;
    OsteoArticularModel(incr_solid).L={'C3_C2JointNode';'C2_C1JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];


%
    %C1_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C1_J2;                   
    OsteoArticularModel(incr_solid).mother=s_C2;
    OsteoArticularModel(incr_solid).a=[0 0 1]';
    OsteoArticularModel(incr_solid).joint=1;  


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-18, C1C2_FE]';


    % OsteoArticularModel(incr_solid).limit_sup=  0*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  0*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=C2_C1JointNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C1_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C1;                   
    OsteoArticularModel(incr_solid).mother=s_C1_J1;
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1;


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-18, C1C2_LB]';


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %C1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).mother=s_C1_J2;
    OsteoArticularModel(incr_solid).child=0; 
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    OsteoArticularModel(incr_solid).joint=1; 


    OsteoArticularModel(incr_solid).linear_constraint=[incr_solid-18, C1C2_AR]';
    OsteoArticularModel(incr_solid).calib_k_constraint = num_solid+19; %Same homothetic coeff.
    OsteoArticularModel(incr_solid).u = [0;0;1];
    OsteoArticularModel(incr_solid).theta = theta_C1(3);


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C1;
    OsteoArticularModel(incr_solid).anat_position=C1_position_set;
    OsteoArticularModel(incr_solid).L={'C2_C1JointNode';'NeckNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];
%}
end


