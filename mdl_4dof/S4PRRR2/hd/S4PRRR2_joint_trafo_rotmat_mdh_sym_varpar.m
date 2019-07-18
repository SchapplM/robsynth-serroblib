% Calculate homogenous joint transformation matrices for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% T_mdh [4x4x4]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S4PRRR2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:19
% EndTime: 2019-07-18 13:27:19
% DurationCPUTime: 0.02s
% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (12->6), ass. (0->7)
t18 = cos(qJ(2));
t17 = cos(qJ(3));
t16 = cos(qJ(4));
t15 = sin(qJ(2));
t14 = sin(qJ(3));
t13 = sin(qJ(4));
t1 = [0, -1, 0, 0; 0, 0, -1, -qJ(1); 1, 0, 0, 0; 0, 0, 0, 1; t18, -t15, 0, 0; t15, t18, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t14, 0, pkin(1); t14, t17, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t13, 0, pkin(2); t13, t16, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,4);             % numerisch
else,                         T_mdh = sym('xx', [4,4,4]); end % symbolisch

for i = 1:4
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
