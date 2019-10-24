% Calculate homogenous joint transformation matrices for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5PPPRR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:17:12
% EndTime: 2019-10-24 10:17:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t51 = cos(qJ(4));
t50 = cos(qJ(5));
t49 = sin(qJ(4));
t48 = sin(qJ(5));
t47 = cos(pkin(7));
t46 = cos(pkin(8));
t45 = cos(pkin(9));
t44 = sin(pkin(7));
t43 = sin(pkin(8));
t42 = sin(pkin(9));
t1 = [t47, -t44, 0, 0; t44, t47, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t46, -t43, 0, pkin(1); 0, 0, -1, -qJ(2); t43, t46, 0, 0; 0, 0, 0, 1; t45, -t42, 0, pkin(2); 0, 0, -1, -qJ(3); t42, t45, 0, 0; 0, 0, 0, 1; t51, -t49, 0, pkin(3); t49, t51, 0, 0; 0, 0, 1, pkin(5); 0, 0, 0, 1; t50, -t48, 0, pkin(4); 0, 0, -1, -pkin(6); t48, t50, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
