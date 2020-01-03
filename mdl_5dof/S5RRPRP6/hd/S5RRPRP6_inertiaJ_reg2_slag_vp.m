% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t29 = t44 * t49 + t45 * t47;
t77 = -0.2e1 * t29;
t63 = -qJ(3) - pkin(6);
t32 = t63 * t49;
t58 = t63 * t47;
t12 = -t44 * t32 - t45 * t58;
t76 = t12 ^ 2;
t27 = t44 * t47 - t45 * t49;
t25 = t27 ^ 2;
t39 = -t49 * pkin(2) - pkin(1);
t75 = 0.2e1 * t39;
t46 = sin(qJ(4));
t74 = 0.2e1 * t46;
t48 = cos(qJ(4));
t73 = -0.2e1 * t48;
t72 = 0.2e1 * t49;
t71 = t27 * pkin(4);
t70 = t44 * pkin(2);
t69 = t45 * pkin(2);
t68 = t48 * pkin(4);
t67 = t46 * t27;
t66 = t46 * t29;
t65 = t46 * t48;
t14 = -t45 * t32 + t44 * t58;
t64 = t48 * t14;
t22 = t48 * t29;
t40 = t46 ^ 2;
t42 = t48 ^ 2;
t33 = t40 + t42;
t41 = t47 ^ 2;
t43 = t49 ^ 2;
t62 = t41 + t43;
t61 = qJ(5) * t29;
t37 = pkin(7) + t70;
t60 = qJ(5) + t37;
t59 = t27 * t77;
t38 = -pkin(3) - t69;
t7 = t27 * pkin(3) - t29 * pkin(7) + t39;
t3 = -t46 * t14 + t48 * t7;
t52 = -t48 * t61 + t3;
t1 = t52 + t71;
t2 = t64 + (t7 - t61) * t46;
t57 = t1 * t48 + t2 * t46;
t4 = t46 * t7 + t64;
t56 = t3 * t48 + t4 * t46;
t55 = -t3 * t46 + t4 * t48;
t23 = t60 * t46;
t24 = t60 * t48;
t54 = -t23 * t48 + t24 * t46;
t53 = -t27 * t37 + t29 * t38;
t35 = 0.2e1 * t65;
t31 = t38 - t68;
t26 = t29 ^ 2;
t21 = t48 * t27;
t20 = t42 * t26;
t17 = t40 * t26;
t16 = t46 * t22;
t15 = -0.2e1 * t26 * t65;
t11 = 0.2e1 * t27 * t22;
t10 = t46 * t59;
t9 = t33 * t29;
t8 = (-t40 + t42) * t29;
t5 = pkin(4) * t66 + t12;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t41, t47 * t72, 0, t43, 0, 0, pkin(1) * t72, -0.2e1 * pkin(1) * t47, 0.2e1 * t62 * pkin(6), t62 * pkin(6) ^ 2 + pkin(1) ^ 2, t26, t59, 0, t25, 0, 0, t27 * t75, t29 * t75, 0.2e1 * t12 * t29 - 0.2e1 * t14 * t27, t14 ^ 2 + t39 ^ 2 + t76, t20, t15, t11, t17, t10, t25, 0.2e1 * t12 * t66 + 0.2e1 * t3 * t27, 0.2e1 * t12 * t22 - 0.2e1 * t4 * t27, t56 * t77, t3 ^ 2 + t4 ^ 2 + t76, t20, t15, t11, t17, t10, t25, 0.2e1 * t1 * t27 + 0.2e1 * t5 * t66, -0.2e1 * t2 * t27 + 0.2e1 * t5 * t22, t57 * t77, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t49, 0, -t47 * pkin(6), -t49 * pkin(6), 0, 0, 0, 0, t29, 0, -t27, 0, -t12, -t14, (-t27 * t44 - t29 * t45) * pkin(2), (-t12 * t45 + t14 * t44) * pkin(2), t16, t8, t67, -t16, t21, 0, -t12 * t48 + t53 * t46, t12 * t46 + t53 * t48, t55, t12 * t38 + t55 * t37, t16, t8, t67, -t16, t21, 0, -t23 * t27 + t31 * t66 - t5 * t48, t31 * t22 - t24 * t27 + t5 * t46, -t1 * t46 + t2 * t48 - t54 * t29, -t1 * t23 + t2 * t24 + t5 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t70, 0, (t44 ^ 2 + t45 ^ 2) * pkin(2) ^ 2, t40, t35, 0, t42, 0, 0, t38 * t73, t38 * t74, 0.2e1 * t33 * t37, t33 * t37 ^ 2 + t38 ^ 2, t40, t35, 0, t42, 0, 0, t31 * t73, t31 * t74, 0.2e1 * t23 * t46 + 0.2e1 * t24 * t48, t23 ^ 2 + t24 ^ 2 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, t39, 0, 0, 0, 0, 0, 0, t21, -t67, -t9, t56, 0, 0, 0, 0, 0, 0, t21, -t67, -t9, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t66, t27, t3, -t4, 0, 0, 0, 0, t22, 0, -t66, t27, t52 + 0.2e1 * t71, -t2, -pkin(4) * t22, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t48, 0, -t46 * t37, -t48 * t37, 0, 0, 0, 0, t46, 0, t48, 0, -t23, -t24, -t46 * pkin(4), -t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t46, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t46, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t22, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t46, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
