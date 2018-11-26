% Calculate homogenous joint transformation matrices for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% T_mdh [4x4x7]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 19:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S7RRRRRRR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 19:00:01
% EndTime: 2018-11-26 19:00:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (14->14), mult. (0->0), div. (0->0), fcn. (28->14), ass. (0->15)
t66 = cos(qJ(1));
t65 = cos(qJ(2));
t64 = cos(qJ(3));
t63 = cos(qJ(4));
t62 = cos(qJ(5));
t61 = cos(qJ(6));
t60 = cos(qJ(7));
t59 = sin(qJ(1));
t58 = sin(qJ(2));
t57 = sin(qJ(3));
t56 = sin(qJ(4));
t55 = sin(qJ(5));
t54 = sin(qJ(6));
t53 = sin(qJ(7));
t1 = [t66, -t59, 0, 0; t59, t66, 0, 0; 0, 0, 1, pkin(1); 0, 0, 0, 1; t65, -t58, 0, 0; 0, 0, -1, 0; t58, t65, 0, 0; 0, 0, 0, 1; t64, -t57, 0, 0; 0, 0, 1, pkin(2); -t57, -t64, 0, 0; 0, 0, 0, 1; t63, -t56, 0, 0; 0, 0, 1, 0; -t56, -t63, 0, 0; 0, 0, 0, 1; t62, -t55, 0, 0; 0, 0, -1, -pkin(3); t55, t62, 0, 0; 0, 0, 0, 1; t61, -t54, 0, 0; 0, 0, -1, 0; t54, t61, 0, 0; 0, 0, 0, 1; t60, -t53, 0, 0; 0, 0, 1, pkin(4); -t53, -t60, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,7);             % numerisch
else,                         T_mdh = sym('xx', [4,4,7]); end % symbolisch

for i = 1:7
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
