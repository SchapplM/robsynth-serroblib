% Calculate homogenous joint transformation matrices for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:47
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RRPPR2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:47:37
% EndTime: 2019-10-24 10:47:37
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t41 = cos(qJ(1));
t40 = cos(qJ(2));
t39 = cos(qJ(5));
t38 = sin(qJ(1));
t37 = sin(qJ(2));
t36 = sin(qJ(5));
t35 = cos(pkin(8));
t34 = cos(pkin(9));
t33 = sin(pkin(8));
t32 = sin(pkin(9));
t1 = [0, 0, 1, pkin(5); -t38, -t41, 0, 0; t41, -t38, 0, 0; 0, 0, 0, 1; t40, -t37, 0, pkin(1); t37, t40, 0, 0; 0, 0, 1, pkin(6); 0, 0, 0, 1; t35, -t33, 0, pkin(2); t33, t35, 0, 0; 0, 0, 1, qJ(3); 0, 0, 0, 1; t34, -t32, 0, pkin(3); 0, 0, -1, -qJ(4); t32, t34, 0, 0; 0, 0, 0, 1; t39, -t36, 0, pkin(4); 0, 0, -1, -pkin(7); t36, t39, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
