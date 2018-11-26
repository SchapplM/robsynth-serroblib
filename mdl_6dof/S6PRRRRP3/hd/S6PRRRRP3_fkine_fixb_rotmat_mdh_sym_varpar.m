% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:28:41
% EndTime: 2018-11-23 15:28:41
% DurationCPUTime: 0.22s
% Computational Cost: add. (597->79), mult. (635->92), div. (0->0), fcn. (717->16), ass. (0->59)
t43 = sin(qJ(4));
t64 = t43 * pkin(4) + pkin(8);
t48 = -pkin(10) - pkin(9);
t39 = qJ(4) + qJ(5);
t32 = sin(t39);
t72 = pkin(5) * t32 + t64;
t46 = cos(qJ(4));
t31 = t46 * pkin(4) + pkin(3);
t71 = cos(qJ(3));
t40 = sin(pkin(11));
t41 = sin(pkin(6));
t70 = t40 * t41;
t42 = cos(pkin(11));
t69 = t42 * t41;
t68 = cos(pkin(6));
t67 = pkin(6) - qJ(2);
t66 = pkin(6) + qJ(2);
t65 = qJ(1) + 0;
t63 = t41 * t71;
t62 = t42 * pkin(1) + pkin(7) * t70 + 0;
t61 = t68 * pkin(7) + t65;
t60 = cos(t66);
t59 = sin(t67);
t55 = sin(t66) / 0.2e1;
t23 = t55 - t59 / 0.2e1;
t47 = cos(qJ(2));
t16 = -t40 * t23 + t42 * t47;
t58 = t16 * pkin(2) + t62;
t56 = cos(t67) / 0.2e1;
t24 = t56 - t60 / 0.2e1;
t57 = t24 * pkin(2) + t61;
t54 = t40 * pkin(1) - pkin(7) * t69 + 0;
t14 = t42 * t23 + t40 * t47;
t53 = t14 * pkin(2) + t54;
t45 = sin(qJ(2));
t49 = t56 + t60 / 0.2e1;
t15 = t40 * t49 + t42 * t45;
t52 = t15 * pkin(8) + t58;
t22 = t55 + t59 / 0.2e1;
t51 = -t22 * pkin(8) + t57;
t13 = t40 * t45 - t42 * t49;
t50 = t13 * pkin(8) + t53;
t44 = sin(qJ(3));
t38 = -qJ(6) + t48;
t33 = cos(t39);
t21 = pkin(5) * t33 + t31;
t18 = t24 * t71 + t68 * t44;
t17 = t24 * t44 - t68 * t71;
t10 = t16 * t71 + t44 * t70;
t9 = t16 * t44 - t40 * t63;
t8 = t14 * t71 - t44 * t69;
t7 = t14 * t44 + t42 * t63;
t6 = t18 * t33 - t22 * t32;
t5 = -t18 * t32 - t22 * t33;
t4 = t10 * t33 + t15 * t32;
t3 = -t10 * t32 + t15 * t33;
t2 = t13 * t32 + t8 * t33;
t1 = t13 * t33 - t8 * t32;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t40, 0, 0; t40, t42, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t16, -t15, t70, t62; t14, -t13, -t69, t54; t24, t22, t68, t61; 0, 0, 0, 1; t10, -t9, t15, t52; t8, -t7, t13, t50; t18, -t17, -t22, t51; 0, 0, 0, 1; t10 * t46 + t15 * t43, -t10 * t43 + t15 * t46, t9, t10 * pkin(3) + t9 * pkin(9) + t52; t13 * t43 + t8 * t46, t13 * t46 - t8 * t43, t7, t8 * pkin(3) + t7 * pkin(9) + t50; t18 * t46 - t22 * t43, -t18 * t43 - t22 * t46, t17, t18 * pkin(3) + t17 * pkin(9) + t51; 0, 0, 0, 1; t4, t3, t9, t10 * t31 + t64 * t15 - t9 * t48 + t58; t2, t1, t7, t64 * t13 + t8 * t31 - t7 * t48 + t53; t6, t5, t17, -t17 * t48 + t18 * t31 - t64 * t22 + t57; 0, 0, 0, 1; t4, t3, t9, t10 * t21 + t72 * t15 - t9 * t38 + t58; t2, t1, t7, t72 * t13 + t8 * t21 - t7 * t38 + t53; t6, t5, t17, -t17 * t38 + t18 * t21 - t72 * t22 + t57; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
