% Calculate homogenous joint transformation matrices for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6RRRPRR15_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:02:58
% EndTime: 2018-11-23 18:02:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->12), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t112 = cos(qJ(1));
t111 = cos(qJ(2));
t110 = cos(qJ(3));
t109 = cos(qJ(5));
t108 = cos(qJ(6));
t107 = sin(qJ(1));
t106 = sin(qJ(2));
t105 = sin(qJ(3));
t104 = sin(qJ(5));
t103 = sin(qJ(6));
t102 = cos(pkin(6));
t101 = cos(pkin(7));
t100 = sin(pkin(6));
t99 = sin(pkin(7));
t1 = [t112, -t107, 0, 0; t107, t112, 0, 0; 0, 0, 1, pkin(8); 0, 0, 0, 1; t111, -t106, 0, pkin(1); t102 * t106, t102 * t111, -t100, -t100 * pkin(9); t100 * t106, t100 * t111, t102, t102 * pkin(9); 0, 0, 0, 1; t110, -t105, 0, pkin(2); t101 * t105, t101 * t110, -t99, -t99 * pkin(10); t99 * t105, t99 * t110, t101, t101 * pkin(10); 0, 0, 0, 1; 0, -1, 0, pkin(3); 0, 0, -1, -qJ(4); 1, 0, 0, 0; 0, 0, 0, 1; t109, -t104, 0, pkin(4); 0, 0, -1, -pkin(11); t104, t109, 0, 0; 0, 0, 0, 1; t108, -t103, 0, pkin(5); 0, 0, -1, -pkin(12); t103, t108, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
