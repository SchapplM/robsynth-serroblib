% Calculate homogenous joint transformation matrices for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_mdh [4x4x4]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S4RRRP4_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 13:55:56
% EndTime: 2019-12-29 13:55:56
% DurationCPUTime: 0.05s
% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (12->6), ass. (0->7)
t25 = cos(qJ(1));
t24 = cos(qJ(2));
t23 = cos(qJ(3));
t22 = sin(qJ(1));
t21 = sin(qJ(2));
t20 = sin(qJ(3));
t1 = [t25, -t22, 0, 0; t22, t25, 0, 0; 0, 0, 1, pkin(4); 0, 0, 0, 1; t24, -t21, 0, pkin(1); 0, 0, -1, -pkin(5); t21, t24, 0, 0; 0, 0, 0, 1; t23, -t20, 0, pkin(2); t20, t23, 0, 0; 0, 0, 1, pkin(6); 0, 0, 0, 1; 1, 0, 0, pkin(3); 0, 1, 0, 0; 0, 0, 1, qJ(4); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,4);             % numerisch
else,                         T_mdh = sym('xx', [4,4,4]); end % symbolisch

for i = 1:4
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
