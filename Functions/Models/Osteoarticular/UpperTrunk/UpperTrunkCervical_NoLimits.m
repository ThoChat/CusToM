function [OsteoArticularModel]= UpperTrunkCervical_NoLimits(OsteoArticularModel,k,Mass,AttachmentPoint,varargin)
%% Implementation of a model of the thoracic spine with 12 thoracic vertebrae
% Ce modèle inclut 12 solides associés aux vertèbres lombaires et articulés
% entre eux suivant 3 degrés de liberté en rotation.
% Auteurs : A. Schuster, G. Dumont, C. Pontonnier (2023)
% Special : - Contraintes cinématiques
%           - Contraintes géométriques
% Fichier résultat : '.C3D' full_constrains
%
%% Les données anthropométriques pour ce modèle sont issues des travaux de
% Bayoglu et al. (2017a,b). 
% (Twente spine model: a complete and coherent dataset for musculo-skeletal
% modeling of the lumbar region of the human spine ;
%  Twente spine model: a complete and coherent dataset for musculo-skeletal
% modeling of the thoracic and cervical regions of the human spine
%
%% Comit 2025, authors: Minh Nguyen, Aurélien Schuster, Aurélie Tomezzoli
% Version sans limites articulaires 
%
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

%__J1 et J2 de chaque solide sont des "solides virtuels" représentant
%respectivement la rotation axiale et la flexion latérale de ce solide.


%% Section 1 : solid numbering incremation / Préliminaires
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
        [OsteoArticularModel]=sister_actualize(OsteoArticularModel,OsteoArticularModel(s_mother).child,eval(['s_' list_solid{1}]));   % recherche de la derni�re soeur
    end
end      

%Matrice de rotation de la base du CT-Scan vers la base suivie par les recommendations de l'ISB
M = [0 -1 0 ;...
     0 0 1  ;...
    -1 0 0  ];


%Matrice de symétrie pour le changement gauche/droite
S = [1 0 0 ;...
     0 1 0 ;...
     0 0 -1];

%Coeff de regression pour le cadavre de 1m54 :
k_new = k * (1.80/1.54);



%disp(a)

%% Section 2 : Calcul de position centre articulation (centre intervertébraux)
%__Position des articulations telle que la moyenne des deux corps de vertèbres (sup+inf)/2 dans le repère du CT-scan
%__Calcul des longueurs des vertèbres (entre les deux corps vertébraux)
%__Calcul des centres de liaisons comme la moyenne entre les deux extrémités sur les corps vertébraux :
%{
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
%}

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


%
MassT1 = 2.7/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT2 = 2.6/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT3 = 3.3/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT4 = 3.1/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT5 = 3.2/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT6 = 3.2/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT7 = 3.4/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT8 = 3.6/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT9 = 3.8/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT10 = 4.8/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT11 = 5.0/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
MassT12 = 6.0/100*(Mass.Pelvis_Mass+Mass.Trunk_Mass);
%}

%% Section 3 : Définition des repères anatomiques pour les segments (CoM wrt global RF ; Node wrt local RF on CoM)
%L1
L1_T12JointNode = L1_T12CT - L2_L1CT; %Repère segment parent

%T12
centre_geom_T12CT = M*k_new*[-0.88880004 -95.19249946 850.85706413]'*0.001;%Repère CT-Scan
CoM_T12CT = centre_geom_T12CT + M*[-0.6 -4.1 0]'*T12_length;
CoM_T12 = CoM_T12CT - L1_T12CT;
T12_T11JointNode = T12_T11CT-L1_T12CT;

%T11
centre_geom_T11CT = M*k_new*[-3.836850765 -90.205750485 879.28863221]'*0.001;%Repère CT-Scan
CoM_T11CT = centre_geom_T11CT + M*[-0.5 -4.4 0]'*T11_length; 
CoM_T11 = CoM_T11CT - T12_T11CT;
T11_T10JointNode = T11_T10CT-T12_T11CT;

%T10
centre_geom_T10CT = M*k_new*[-4.74510698 -88.076947225 905.66991757]'*0.001;%Repère CT-Scan
CoM_T10CT = centre_geom_T10CT + M*[-0.2 -4.6 0]'*T10_length; 
CoM_T10 = CoM_T10CT - T11_T10CT;
T10_T9JointNode = T10_T9CT-T11_T10CT;

%T9
centre_geom_T9CT = M*k_new*[-5.831965935 -87.096696205 930.94956846]'*0.001;%Repère CT-Scan
CoM_T9CT = centre_geom_T9CT + M*[-0.1 -4.6 0]'*T9_length; 
CoM_T9 = CoM_T9CT - T10_T9CT;
T9_T8JointNode = T9_T8CT-T10_T9CT;

%T8
centre_geom_T8CT = M*k_new*[-7.09154992 -86.621049795 955.103456975]'*0.001;%Repère CT-Scan
CoM_T8CT = centre_geom_T8CT + M*[-0.1 -4.6 0]'*T8_length; 
CoM_T8 = CoM_T8CT - T9_T8CT;
T8_T7JointNode = T8_T7CT-T9_T8CT;

%T7
centre_geom_T7CT = M*k_new*[-7.872599875 -90.458102525 978.58589888]'*0.001;%Repère CT-Scan
CoM_T7CT = centre_geom_T7CT + M*[0.5 -4.3 0]'*T7_length; 
CoM_T7 = CoM_T7CT - T8_T7CT;
T7_T6JointNode = T7_T6CT-T8_T7CT;

%T6
centre_geom_T6CT = M*k_new*[-7.65406297 -98.26729827 1000.734251955]'*0.001;%Repère CT-Scan
CoM_T6CT = centre_geom_T6CT + M*[0.6 -3.9 0]'*T6_length; 
CoM_T6 = CoM_T6CT - T7_T6CT;
T6_T5JointNode = T6_T5CT-T7_T6CT;

%T5
centre_geom_T5CT = M*k_new*[-6.150750095 -106.459300965 1020.276397465]'*0.001;%Repère CT-Scan
CoM_T5CT = centre_geom_T5CT + M*[0.4 -3.4 0]'*T5_length; 
CoM_T5 = CoM_T5CT - T6_T5CT;
T5_T4JointNode = T5_T4CT-T6_T5CT;

%T4
centre_geom_T4CT = M*k_new*[-1.886976905 -117.97578271 1037.791405535]'*0.001;%Repère CT-Scan
CoM_T4CT = centre_geom_T4CT + M*[0.1 -2.8 0]'*T4_length;
CoM_T4 = CoM_T4CT - T5_T4CT;
T4_T3JointNode = T4_T3CT-T5_T4CT;

%T3
centre_geom_T3CT = M*k_new*[2.029917365 -130.95841247 1053.39871453]'*0.001;%Repère CT-Scan
CoM_T3CT = centre_geom_T3CT + M*[-0.2 -2.0 0]'*T3_length; 
CoM_T3 = CoM_T3CT - T4_T3CT;
T3_T2JointNode = T3_T2CT-T4_T3CT;

%T2
centre_geom_T2CT = M*k_new*[4.85567504 -142.766151575 1069.3564713]'*0.001;%Repère CT-Scan
CoM_T2CT = centre_geom_T2CT + M*[-0.4 -1.3 0]'*T2_length; 
CoM_T2 = CoM_T2CT - T3_T2CT;
T2_T1JointNode = T2_T1CT-T3_T2CT;

%T1
centre_geom_T1CT = M*k_new*[9.70759988 -154.560297725 1084.445774555]'*0.001;%Repère CT-Scan
CoM_T1CT = centre_geom_T1CT + M*[-0.5 -0.8 0]'*T1_length;
CoM_T1 = CoM_T1CT - T2_T1CT;
T1_C7JointNode = T1_C7CT-T2_T1CT;

%C7
centre_geom_C7CT = M*k_new*[17.232949845 -165.89352861 1097.79238701]'*0.001;
CoM_C7CT = centre_geom_C7CT;
CoM_C7 = CoM_C7CT - T1_C7CT;
C7_C6JointNode = C7_C6CT-T1_C7CT;

%C6
centre_geom_C6CT = M*k_new*[22.199249825 -174.04229939 1109.831154345]'*0.001;
CoM_C6CT = centre_geom_C6CT;
CoM_C6 = CoM_C6CT - C7_C6CT;
C6_C5JointNode = C6_C5CT-C7_C6CT;

%C5
centre_geom_C5CT = M*k_new*[23.95494934 -182.11659789 1121.572822335]'*0.001;
CoM_C5CT = centre_geom_C5CT;
CoM_C5 = CoM_C5CT - C6_C5CT;
C5_C4JointNode = C5_C4CT - C6_C5CT;

%C4
centre_geom_C4CT = M*k_new*[27.027212395 -188.740900685 1136.014707235]'*0.001;
CoM_C4CT = centre_geom_C4CT;
CoM_C4 = CoM_C4CT - C5_C4CT;
C4_C3JointNode = C4_C3CT - C5_C4CT;

%C3
centre_geom_C3CT = M*k_new*[28.592766545 -193.79058529 1151.91299904]'*0.001;
CoM_C3CT = centre_geom_C3CT;
CoM_C3 = CoM_C3CT - C4_C3CT;
C3_C2JointNode = C3_C2CT - C4_C3CT;

%C2
centre_geom_C2CT = M*k_new*[28.463409795 -180.365311845 1171.135298765]'*0.001;
CoM_C2CT = centre_geom_C2CT;
CoM_C2 = CoM_C2CT - C3_C2CT;
C2_C1JointNode = C2_C1CT - C3_C2CT;

%C1
centre_geom_C1CT = M*k_new*[31.24669939 -202.598270025 1191.16777182]'*0.001;
CoM_C1CT = centre_geom_C1CT;
CoM_C1 = CoM_C1CT - C2_C1CT;
C1_HeadJointNode = C1_sup - C2_C1CT; %__Attachement de la tête


% Attachement Bras %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Thorax_RClavJointNodeCT = M*k_new*[-10.7812 -205.974 1051.98]'*0.001; %Point d'attache RClavicule. Comme dans OpenSim
Thorax_scjJointRightNode = Thorax_RClavJointNodeCT - T2_T1CT; % Point d'attache RClavicle dans le repère T1
Thorax_scjJointLeftNode = S*Thorax_scjJointRightNode;

RScapula_HumerusJointNodeCT = M*k_new*[-130.311 -184.256 1118.01]'*0.001;%Tète superieure humérus
Thorax_ShoulderRightNode = RScapula_HumerusJointNodeCT - Thorax_RClavJointNodeCT; %Bras droit dans le repère RClavicle
Thorax_ShoulderLeftNode = S*Thorax_ShoulderRightNode;  %Bras gauche dans le repère LClavicle


%% Section 2 :Définition de la position des marqueurs attachés au modèle (wrt local RF on CoM)
T12_position_set = {
    % Node (repère local)
    'Thorax_Origin', L1_T12JointNode-CoM_T12;...
    'T12_Origin', L1_T12JointNode-CoM_T12;...
    'T12_T11Node', T12_T11JointNode-CoM_T12;...
    'BARYCENTRE', centre_geom_T12CT-CoM_T12CT;...
    % Markers (repère global)
    'T12', M*k_new*[-5.1905 -40.3098 841.586]'*0.001- CoM_T12CT;...
    'T12POSSUP', M*k_new*[-2.95 -80.73 858.31]'*0.001- CoM_T12CT;...
    };

T11_position_set = {
    % Node (repère local)
    'T11_Origin', [0 0 0]'-CoM_T11;...
    'T11_T10Node', T11_T10JointNode-CoM_T11;...
    % Markers (repère global)
    'T11', M*k_new*0.001*[-6.255 -34.7616 870.901]' - CoM_T11CT;...
    'T11POSSUP', M*k_new*0.001*[-3.445 -76.365 887.435]' - CoM_T11CT;... 
    };

T10_position_set = {
    % Node (repère local)
    'T10_Origin', [0 0 0]'-CoM_T10;...
    'T10_T9Node', T10_T9JointNode-CoM_T10;...
    % Markers (repère global)
    'T10', M*k_new*[-6.5779 -35.4064 887.06]'*0.001 - CoM_T10CT;...
    'T10POSSUP', M*k_new*[-5.325 -74.385 916.365]'*0.001 - CoM_T10CT;...
    };

T9_position_set = {
    % Node (repère local)
    'T9_Origin', [0 0 0]'-CoM_T9;...
    'T9_T8Node', T9_T8JointNode-CoM_T9;...
    % Markers (repère global)
    'T9', M*k_new*0.001*[-9.0464 -33.595 903.156]' - CoM_T9CT;...
    'T9POSSUP', M*k_new*0.001*[-6.08 -74.485 941.285]' - CoM_T9CT;... 
    };

T8_position_set = {
    % Node (repère local)
    'T8_Origin', [0 0 0]'-CoM_T8;...
    'T8_T7Node', T8_T7JointNode-CoM_T8;...
    % Markers (repère global)
    'T8', M*k_new*[-6.7407 -31.9621 926.954]'*0.001- CoM_T8CT;...
    'T8POSSUP', M*k_new*[-7.225 -77.765 966.7]'*0.001- CoM_T8CT;...
    };

T7_position_set = {
    % Node (repère local)
    'T7_Origin', [0 0 0]'-CoM_T7;...
    'T7_T6Node', T7_T6JointNode-CoM_T7;...
    % Markers (repère global)
    'T7', M*k_new*0.001*[-6.3325 -35.1669 956.617]' - CoM_T7CT;...
    'T7POSSUP', M*k_new*0.001*[-7.9 -82.46 989.585]' - CoM_T7CT;...
    };

T6_position_set = {
    % Node (repère local)
    'T6_Origin', [0 0 0]'-CoM_T6;...
    'T6_T5Node', T6_T5JointNode-CoM_T6;...
    % Markers (repère global)
    'T6', M*k_new*0.001*[-7.6284 -40.8255 983.7]' - CoM_T6CT;...
    'T6POSSUP', M*k_new*0.001*[-6.41 -90.695 1012.37]' - CoM_T6CT;...
    };

T5_position_set = {
    % Node (repère local)
    'T5_Origin', [0 0 0]'-CoM_T5;...
    'T5_T4Node', T5_T4JointNode-CoM_T5;...
    % Markers (repère global)
    'T5', M*k_new*0.001*[-7.5567 -48.8584 1009.16]' - CoM_T5CT;...
    'T5POSSUP', M*k_new*0.001*[-3.335 -102.075 1032.205]' - CoM_T5CT;...
    };

T4_position_set = {
    % Node (repère local)
    'T4_Origin', [0 0 0]'-CoM_T4;...
    'T4_T3Node', T4_T3JointNode-CoM_T4;...
    % Markers (repère global)
    'T4', M*k_new*0.001*[-10.1477 -60.1759 1035.06]' - CoM_T4CT;...
    'T4POSSUP', M*k_new*0.001*[0.28 -115.12 1049.495]' - CoM_T4CT;...
    };

T3_position_set = {
    % Node (repère local)
    'T3_Origin', [0 0 0]'-CoM_T3;...
    'T3_T2Node', T3_T2JointNode-CoM_T3;...
    % Markers (repère global)
    'T3', M*k_new*0.001*[-6.0426 -73.0513 1057.13]' - CoM_T3CT;...
    'T3POSSUP', M*k_new*0.001*[3.96 -128.01 1064.805]' - CoM_T3CT;...
    'CLAV', M*k_new*[11.2451 -206.483 1046.82]'*0.001 - CoM_T3CT;... %Protubérance manubrium entre les deux clavicules
    'STRN', M*k_new*[5.4031 -235.578 913.896]'*0.001 - CoM_T3CT;... %__STRN déplacé du repère T10 à T3 (Minh)
    };

T2_position_set = {
    % Node (repère local)
    'T2_Origin', [0 0 0]'-CoM_T2;...
    'T2_T1Node', T2_T1JointNode-CoM_T2;...
    % Markers (repère global)
    'T2', M*k_new*0.001*[-1.1896 -87.057 1081.09]' - CoM_T2CT;...
    'T2POSSUP', M*k_new*0.001*[6.335 -141.245 1079.735]' - CoM_T2CT;...
    };

T1_position_set = {
    % Node (repère local)
    'T1_Origin', [0 0 0]'-CoM_T1;...
    'T1_C7Node', T1_C7JointNode-CoM_T1;...
    %'NeckNode', NeckNode-CoM_T1;...
    'Thorax_scjJointRightNode', Thorax_scjJointRightNode-CoM_T1; ...
    'Thorax_scjJointLeftNode', Thorax_scjJointLeftNode-CoM_T1; ...
    %'C7POSSUP', M*k_new*0.001*[19.8 -162.8 1107]' - CoM_T1CT;...
    'T1', M*k_new*0.001*[1.5311 -105.019 1102.14]' - CoM_T1CT;...
    'T1POSSUP', M*k_new*0.001*[12.83 -153.845 1094.21]' - CoM_T1CT;...
    };

%
C7_position_set = {
    % Node (repère local)
    'C7_Origin', [0 0 0]'-CoM_C7;...
    'C7_C6Node', C7_C6JointNode-CoM_C7;...
    % Markers (repère global)
    'C7', M*k_new*0.001*[14.3472 -115.937 1117.14]'- CoM_C7CT;... %Scoliosis adjustment
    'C7POSSUP', M*k_new*0.001*[19.8 -162.8 1107]' - CoM_C7CT;...
    };
    
C6_position_set = {
    % Node (repère local)
    'C6_Origin', [0 0 0]'-CoM_C6;... %__Coordonnées de l'origine du repère global exprimé dans le repère local de C6
    'C6_C5Node', C6_C5JointNode-CoM_C6;...
    % Markers (repère global)
    'C6', M*k_new*0.001*[20.8273 -139.225 1130.94]' - CoM_C6CT;...
    'C6POSSUP', M*k_new*0.001*[23.24 -169.75 1118.66]' - CoM_C6CT;...
    };

C5_position_set = {
    % Node (repère local)
    'C5_Origin', [0 0 0]'-CoM_C5;...
    'C5_C4Node', C5_C4JointNode-CoM_C5;...
    % Markers (repère global)
    'C5', M*k_new*0.001*[21.0193 -153.035 1137.8]' - CoM_C5CT;...
    'C5POSSUP', M*k_new*0.001*[25.305 -179.02 1130.8]' - CoM_C5CT;...
    };

C4_position_set = {
    % Node (repère local)
    'C4_Origin', [0 0 0]'-CoM_C4;...
    'C4_C3Node', C4_C3JointNode-CoM_C4;...
    %'NeckNode', NeckNode - CoM_C4;...
    % Markers (repère global)
    'C4', M*k_new*0.001*[17.89 -150.15 1146.94]' - CoM_C4CT;...
    'C4POSSUP', M*k_new*0.001*[27.125 -185.295 1143.96]' - CoM_C4CT;...
    };

%
C3_position_set = {
    % Node (repère local)
    'C3_Origin', [0 0 0]'-CoM_C3;...
    'C3_C2Node',C3_C2JointNode-CoM_C3;...
    %'NeckNode', NeckNode - CoM_C3;...
    % Markers (repère global)
    'C3', M*k_new*0.001*[30.2856 -155.289 1158.88]' - CoM_C3CT;...
    'C3POSSUP', M*k_new*0.001*[28.315 -187.755 1159.215]' - CoM_C3CT;...
    };

%
C2_position_set = {
    % Node (repère local)
    'C2_Origin', [0 0 0]'-CoM_C2;...
    'C2_C1Node',C2_C1JointNode-CoM_C2;...
    %'NeckNode', NeckNode-CoM_C2;...
    % Markers (repère global)
    'C2', M*k_new*0.001*[26.8443 -154.133 1169.28]' - CoM_C2CT;...
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


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup=  inf;
    % OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    OsteoArticularModel(incr_solid).limit_inf=-5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 5*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T12;
    OsteoArticularModel(incr_solid).m=MassT12;               
    OsteoArticularModel(incr_solid).I=IT12;
    OsteoArticularModel(incr_solid).anat_position=T12_position_set;
    OsteoArticularModel(incr_solid).L={'Thorax_Origin';'T12_T11JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR12.mat'];

    %T11_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T11_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T12;
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';

    
    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]'; 
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe
    

    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T11;
    OsteoArticularModel(incr_solid).m=MassT11;               
    OsteoArticularModel(incr_solid).I=IT11;
    OsteoArticularModel(incr_solid).anat_position=T11_position_set;
    OsteoArticularModel(incr_solid).L={'T12_T11JointNode';'T11_T10JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];

    %T10_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T10_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T11;
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1; 


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    OsteoArticularModel(incr_solid).limit_inf=-3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 3*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T10;
    OsteoArticularModel(incr_solid).m=MassT10;               
    OsteoArticularModel(incr_solid).I=IT10;
    OsteoArticularModel(incr_solid).anat_position=T10_position_set;
    OsteoArticularModel(incr_solid).L={'T11_T10JointNode';'T10_T9JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR10.mat'];

    %T9_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T9_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T10;
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T9;
    OsteoArticularModel(incr_solid).m=MassT9;               
    OsteoArticularModel(incr_solid).I=IT9;
    OsteoArticularModel(incr_solid).anat_position=T9_position_set;
    OsteoArticularModel(incr_solid).L={'T10_T9JointNode';'T9_T8JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR9.mat'];

    %T8_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T8_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T9; 
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).a=[0 0 1]';



    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    

    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1; 


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T8;
    OsteoArticularModel(incr_solid).m=MassT8;               
    OsteoArticularModel(incr_solid).I=IT8;
    OsteoArticularModel(incr_solid).anat_position=T8_position_set;
    OsteoArticularModel(incr_solid).L={'T9_T8JointNode';'T8_T7JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR8.mat'];

    %T7_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T7_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T8; 
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;  


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1; 


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);

    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T7;
    OsteoArticularModel(incr_solid).m=MassT7;               
    OsteoArticularModel(incr_solid).I=IT7;
    OsteoArticularModel(incr_solid).anat_position=T7_position_set;
    OsteoArticularModel(incr_solid).L={'T8_T7JointNode';'T7_T6JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR7.mat'];

    %T6_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T6_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T7;
    OsteoArticularModel(incr_solid).joint=1; 
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;  


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    OsteoArticularModel(incr_solid).limit_inf=-4*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 4*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T6;
    OsteoArticularModel(incr_solid).m=MassT6;               
    OsteoArticularModel(incr_solid).I=IT6;
    OsteoArticularModel(incr_solid).anat_position=T6_position_set;
    OsteoArticularModel(incr_solid).L={'T7_T6JointNode';'T6_T5JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR6.mat'];

    %T5_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T5_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T6;
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;  


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1; 


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    OsteoArticularModel(incr_solid).limit_inf=-3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 3*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T5;
    OsteoArticularModel(incr_solid).m=MassT5;               
    OsteoArticularModel(incr_solid).I=IT5;
    OsteoArticularModel(incr_solid).anat_position=T5_position_set;
    OsteoArticularModel(incr_solid).L={'T6_T5JointNode';'T5_T4JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR5.mat'];

    %T4_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T4_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T5;
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    OsteoArticularModel(incr_solid).limit_inf=-3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 3*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T4;
    OsteoArticularModel(incr_solid).m=MassT4;               
    OsteoArticularModel(incr_solid).I=IT4;
    OsteoArticularModel(incr_solid).anat_position=T4_position_set;
    OsteoArticularModel(incr_solid).L={'T5_T4JointNode';'T4_T3JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR4.mat'];

    %T3_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T3_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T4;
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;  


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1; 


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    OsteoArticularModel(incr_solid).limit_inf=-3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 3*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T3;
    OsteoArticularModel(incr_solid).m=MassT3;               
    OsteoArticularModel(incr_solid).I=IT3;
    OsteoArticularModel(incr_solid).anat_position=T3_position_set;
    OsteoArticularModel(incr_solid).L={'T4_T3JointNode';'T3_T2JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR3.mat'];
     
    %T2_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T2_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T3;
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]'; 
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    OsteoArticularModel(incr_solid).limit_inf=-3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 3*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T2;
    OsteoArticularModel(incr_solid).m=MassT2;               
    OsteoArticularModel(incr_solid).I=IT2;
    OsteoArticularModel(incr_solid).anat_position=T2_position_set;
    OsteoArticularModel(incr_solid).L={'T3_T2JointNode';'T2_T1JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR2.mat'];
     
    %T1_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_T1_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T2;
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 


    % OsteoArticularModel(incr_solid).limit_sup=  1.67*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -2.5*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    OsteoArticularModel(incr_solid).limit_inf=-3*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 3*(2*pi/360);

    
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_T1;
    OsteoArticularModel(incr_solid).m=MassT1;               
    OsteoArticularModel(incr_solid).I=IT1;
    OsteoArticularModel(incr_solid).anat_position=T1_position_set;
    OsteoArticularModel(incr_solid).L={'T2_T1JointNode';'T1_C7JointNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR1.mat'];

    %RClavicle_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=s_LClavicle_J1;                
    OsteoArticularModel(incr_solid).child=s_RClavicle_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T1; 
    OsteoArticularModel(incr_solid).a=[0 0 1]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_inf=-pi/2; % Comme dans UpperTrunkClavicle
    OsteoArticularModel(incr_solid).limit_sup= pi/2;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=Thorax_scjJointRightNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);
     
    %RClavicle_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_RClavicle;                   
    OsteoArticularModel(incr_solid).mother=s_RClavicle_J1;
    OsteoArticularModel(incr_solid).a=[1 0 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_inf=-pi/2;
    OsteoArticularModel(incr_solid).limit_sup= pi/2;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %RClavicle
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=0;                  
    OsteoArticularModel(incr_solid).mother=s_RClavicle_J2;
    OsteoArticularModel(incr_solid).a=[0 1 0]'; 
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_inf=-pi/2;
    OsteoArticularModel(incr_solid).limit_sup= pi/2;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;  
    OsteoArticularModel(incr_solid).I=zeros(3,3);
    OsteoArticularModel(incr_solid).anat_position=RClavicle_position_set;
    OsteoArticularModel(incr_solid).L={'RClavicle_Origin';'Thorax_ShoulderRightNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/clav_scap_R.mat'];
    % OsteoArticularModel(incr_solid).visual_file = ['Twente/RClavicle.mat'];

    %LClavicle_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=s_C7_J1;
    OsteoArticularModel(incr_solid).child=s_LClavicle_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_inf=-pi/2;
    OsteoArticularModel(incr_solid).limit_sup= pi/2;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=Thorax_scjJointLeftNode;
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);
     
    %LClavicle_J2
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_LClavicle;                   
    OsteoArticularModel(incr_solid).mother=s_LClavicle_J1;
    OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_inf=-pi/2;
    OsteoArticularModel(incr_solid).limit_sup= pi/2;
    OsteoArticularModel(incr_solid).Visual=0;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;                 
    OsteoArticularModel(incr_solid).I=zeros(3,3);

    %LClavicle
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=0;                  
    OsteoArticularModel(incr_solid).mother=s_LClavicle_J2;
    OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).joint=1;  
    OsteoArticularModel(incr_solid).limit_inf=-pi/2;
    OsteoArticularModel(incr_solid).limit_sup= pi/2;
    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=[0 0 0]';
    OsteoArticularModel(incr_solid).m=0;  
    OsteoArticularModel(incr_solid).I=zeros(3,3);
    OsteoArticularModel(incr_solid).anat_position=LClavicle_position_set;
    OsteoArticularModel(incr_solid).L={'LClavicle_Origin';'Thorax_ShoulderLeftNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/clav_scap_L.mat'];
    % OsteoArticularModel(incr_solid).visual_file = ['Twente/LClavicle.mat'];


%
    %C7_J1
    num_solid=num_solid+1;                        % number of the solid ...
    name=list_solid{num_solid};                                % solid name
    eval(['incr_solid=s_' name ';'])     % number of the solid in the model
    OsteoArticularModel(incr_solid).name=name;                 % solid name
    OsteoArticularModel(incr_solid).sister=0;                
    OsteoArticularModel(incr_solid).child=s_C7_J2;                   
    OsteoArticularModel(incr_solid).mother=s_T1;
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


    % OsteoArticularModel(incr_solid).limit_sup=  (25/4)*(2*pi/360); % From René Louis "Anatomy of the spine"
    % OsteoArticularModel(incr_solid).limit_inf=  -(15/4)*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1; 


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 


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
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


    % OsteoArticularModel(incr_solid).limit_inf=-1.67*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 1.67*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  inf;
    OsteoArticularModel(incr_solid).limit_inf=  -inf;


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


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
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  0.00001;
    OsteoArticularModel(incr_solid).limit_inf=  -0.00001;


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
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);


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
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';
    


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; %__Déactiver la rotation suivant cet axe


    % OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    % OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup=  0.00001;
    OsteoArticularModel(incr_solid).limit_inf=  -0.00001;


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
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    


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
    OsteoArticularModel(incr_solid).joint=1; 


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]'; 
    

    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);


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
    OsteoArticularModel(incr_solid).joint=1;
    OsteoArticularModel(incr_solid).a=[0 0 1]';


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[1 0 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


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
    OsteoArticularModel(incr_solid).joint=1;


    %OsteoArticularModel(incr_solid).a=[0 1 0]';
    OsteoArticularModel(incr_solid).a=[0 0 0]';


    OsteoArticularModel(incr_solid).limit_inf=-2*(2*pi/360);
    OsteoArticularModel(incr_solid).limit_sup= 2*(2*pi/360);


    OsteoArticularModel(incr_solid).Visual=1;
    OsteoArticularModel(incr_solid).b=[0 0 0]';
    OsteoArticularModel(incr_solid).c=CoM_C1;
    OsteoArticularModel(incr_solid).anat_position=C1_position_set;
    OsteoArticularModel(incr_solid).L={'C2_C1JointNode';'NeckNode'};
    OsteoArticularModel(incr_solid).visual_file = ['Twente/TR11.mat'];
%}
end


