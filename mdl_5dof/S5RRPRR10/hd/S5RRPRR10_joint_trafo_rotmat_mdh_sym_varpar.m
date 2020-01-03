% Calculate homogenous joint transformation matrices for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RRPRR10_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:38
% EndTime: 2019-12-31 20:23:38
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t73 = cos(qJ(1));
t72 = cos(qJ(2));
t71 = cos(qJ(4));
t70 = cos(qJ(5));
t69 = sin(qJ(1));
t68 = sin(qJ(2));
t67 = sin(qJ(4));
t66 = sin(qJ(5));
t65 = cos(pkin(5));
t64 = cos(pkin(10));
t63 = sin(pkin(5));
t62 = sin(pkin(10));
t1 = [t73, -t69, 0, 0; t69, t73, 0, 0; 0, 0, 1, pkin(6); 0, 0, 0, 1; t72, -t68, 0, pkin(1); t65 * t68, t65 * t72, -t63, -t63 * pkin(7); t63 * t68, t63 * t72, t65, t65 * pkin(7); 0, 0, 0, 1; t64, -t62, 0, pkin(2); t62, t64, 0, 0; 0, 0, 1, qJ(3); 0, 0, 0, 1; t71, -t67, 0, pkin(3); 0, 0, -1, -pkin(8); t67, t71, 0, 0; 0, 0, 0, 1; t70, -t66, 0, pkin(4); 0, 0, -1, -pkin(9); t66, t70, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
