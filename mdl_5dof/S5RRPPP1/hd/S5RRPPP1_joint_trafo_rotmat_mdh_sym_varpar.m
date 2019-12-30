% Calculate homogenous joint transformation matrices for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RRPPP1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 18:07:26
% EndTime: 2019-12-29 18:07:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (8->8), mult. (6->6), div. (0->0), fcn. (22->8), ass. (0->9)
t66 = cos(qJ(1));
t65 = cos(qJ(2));
t64 = sin(qJ(1));
t63 = sin(qJ(2));
t62 = cos(pkin(5));
t61 = cos(pkin(8));
t60 = sin(pkin(5));
t59 = sin(pkin(8));
t1 = [t66, -t64, 0, 0; t64, t66, 0, 0; 0, 0, 1, pkin(6); 0, 0, 0, 1; t65, -t63, 0, pkin(1); 0, 0, -1, -pkin(7); t63, t65, 0, 0; 0, 0, 0, 1; t61, -t59, 0, pkin(2); t62 * t59, t62 * t61, -t60, -t60 * qJ(3); t60 * t59, t60 * t61, t62, t62 * qJ(3); 0, 0, 0, 1; 0, -1, 0, pkin(3); 0, 0, -1, -qJ(4); 1, 0, 0, 0; 0, 0, 0, 1; 1, 0, 0, pkin(4); 0, 0, -1, -qJ(5); 0, 1, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
