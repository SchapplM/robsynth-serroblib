% Calculate homogenous joint transformation matrices for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:31
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5PRRPR6_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:31:45
% EndTime: 2019-10-24 10:31:45
% DurationCPUTime: 0.04s
% Computational Cost: add. (9->9), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t65 = cos(qJ(2));
t64 = cos(qJ(3));
t63 = cos(qJ(5));
t62 = sin(qJ(2));
t61 = sin(qJ(3));
t60 = sin(qJ(5));
t59 = cos(pkin(5));
t58 = cos(pkin(9));
t57 = cos(pkin(10));
t56 = sin(pkin(5));
t55 = sin(pkin(9));
t54 = sin(pkin(10));
t1 = [t58, -t55, 0, 0; t55, t58, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t65, -t62, 0, pkin(1); t59 * t62, t59 * t65, -t56, -t56 * pkin(6); t56 * t62, t56 * t65, t59, t59 * pkin(6); 0, 0, 0, 1; t64, -t61, 0, pkin(2); 0, 0, -1, -pkin(7); t61, t64, 0, 0; 0, 0, 0, 1; t57, -t54, 0, pkin(3); 0, 0, -1, -qJ(4); t54, t57, 0, 0; 0, 0, 0, 1; t63, -t60, 0, pkin(4); t60, t63, 0, 0; 0, 0, 1, pkin(8); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
