% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR13_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t87 = -t56 * t53 + t59 * t54;
t46 = -t54 * pkin(3) - pkin(2);
t27 = -pkin(4) * t87 + t46;
t86 = 0.2e1 * t27;
t85 = 0.2e1 * t46;
t60 = cos(qJ(2));
t84 = -0.2e1 * t60;
t83 = 0.2e1 * t60;
t57 = sin(qJ(2));
t51 = t57 ^ 2;
t82 = t51 * pkin(6);
t55 = sin(qJ(5));
t81 = t55 * pkin(4);
t48 = t57 * pkin(6);
t58 = cos(qJ(5));
t80 = t58 * pkin(4);
t37 = t59 * t53 + t56 * t54;
t28 = t37 * t57;
t39 = -t60 * pkin(2) - t57 * qJ(3) - pkin(1);
t33 = t54 * t39;
t74 = t54 * t57;
t19 = -pkin(7) * t74 + t33 + (-pkin(6) * t53 - pkin(3)) * t60;
t73 = t54 * t60;
t25 = pkin(6) * t73 + t53 * t39;
t76 = t53 * t57;
t21 = -pkin(7) * t76 + t25;
t9 = t56 * t19 + t59 * t21;
t7 = -t28 * pkin(8) + t9;
t79 = t58 * t7;
t78 = t60 * pkin(4);
t77 = t53 * t54;
t75 = t53 * t60;
t71 = t57 * t60;
t69 = pkin(7) + qJ(3);
t38 = pkin(3) * t76 + t48;
t49 = t53 ^ 2;
t50 = t54 ^ 2;
t68 = t49 + t50;
t67 = 0.2e1 * t71;
t66 = t53 * t74;
t30 = t87 * t57;
t8 = t59 * t19 - t56 * t21;
t4 = -t30 * pkin(8) - t78 + t8;
t1 = t58 * t4 - t55 * t7;
t40 = t69 * t53;
t41 = t69 * t54;
t22 = -t59 * t40 - t56 * t41;
t65 = -pkin(2) * t57 + qJ(3) * t60;
t24 = -pkin(6) * t75 + t33;
t64 = -t24 * t53 + t25 * t54;
t23 = -t56 * t40 + t59 * t41;
t62 = pkin(6) ^ 2;
t52 = t60 ^ 2;
t47 = t51 * t62;
t20 = t28 * pkin(4) + t38;
t18 = t58 * t37 + t55 * t87;
t16 = t55 * t37 - t58 * t87;
t14 = pkin(8) * t87 + t23;
t13 = -t37 * pkin(8) + t22;
t12 = -t55 * t28 + t58 * t30;
t10 = t58 * t28 + t55 * t30;
t6 = t55 * t13 + t58 * t14;
t5 = t58 * t13 - t55 * t14;
t2 = t55 * t4 + t79;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t51, t67, 0, t52, 0, 0, pkin(1) * t83, -0.2e1 * pkin(1) * t57, 0.2e1 * (t51 + t52) * pkin(6), pkin(1) ^ 2 + t52 * t62 + t47, t50 * t51, -0.2e1 * t51 * t77, -0.2e1 * t54 * t71, t49 * t51, t53 * t67, t52, -0.2e1 * t24 * t60 + 0.2e1 * t53 * t82, 0.2e1 * t25 * t60 + 0.2e1 * t54 * t82, 0.2e1 * (-t24 * t54 - t25 * t53) * t57, t24 ^ 2 + t25 ^ 2 + t47, t30 ^ 2, -0.2e1 * t30 * t28, t30 * t84, t28 ^ 2, -t28 * t84, t52, 0.2e1 * t38 * t28 - 0.2e1 * t8 * t60, 0.2e1 * t38 * t30 + 0.2e1 * t9 * t60, -0.2e1 * t9 * t28 - 0.2e1 * t8 * t30, t38 ^ 2 + t8 ^ 2 + t9 ^ 2, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t84, t10 ^ 2, t10 * t83, t52, -0.2e1 * t1 * t60 + 0.2e1 * t20 * t10, 0.2e1 * t20 * t12 + 0.2e1 * t2 * t60, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t60, 0, -t48, -t60 * pkin(6), 0, 0, t66, (-t49 + t50) * t57, -t75, -t66, -t73, 0, -pkin(6) * t74 + t53 * t65, pkin(6) * t76 + t54 * t65, t64, -pkin(2) * t48 + qJ(3) * t64, t30 * t37, -t37 * t28 + t30 * t87, -t37 * t60, -t28 * t87, -t87 * t60, 0, -t22 * t60 + t46 * t28 - t38 * t87, t23 * t60 + t46 * t30 + t38 * t37, -t22 * t30 - t23 * t28 - t8 * t37 + t87 * t9, t8 * t22 + t9 * t23 + t38 * t46, t12 * t18, -t18 * t10 - t12 * t16, -t18 * t60, t10 * t16, t16 * t60, 0, t27 * t10 + t20 * t16 - t5 * t60, t27 * t12 + t20 * t18 + t6 * t60, -t1 * t18 - t6 * t10 - t5 * t12 - t2 * t16, t1 * t5 + t2 * t6 + t20 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t49, 0.2e1 * t77, 0, t50, 0, 0, 0.2e1 * pkin(2) * t54, -0.2e1 * pkin(2) * t53, 0.2e1 * t68 * qJ(3), qJ(3) ^ 2 * t68 + pkin(2) ^ 2, t37 ^ 2, 0.2e1 * t37 * t87, 0, t87 ^ 2, 0, 0, -t87 * t85, t37 * t85, -0.2e1 * t22 * t37 + 0.2e1 * t23 * t87, t22 ^ 2 + t23 ^ 2 + t46 ^ 2, t18 ^ 2, -0.2e1 * t18 * t16, 0, t16 ^ 2, 0, 0, t16 * t86, t18 * t86, -0.2e1 * t6 * t16 - 0.2e1 * t5 * t18, t27 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t74, 0, t48, 0, 0, 0, 0, 0, 0, t28, t30, 0, t38, 0, 0, 0, 0, 0, 0, t10, t12, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t53, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t87, t37, 0, t46, 0, 0, 0, 0, 0, 0, t16, t18, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t28, -t60, t8, -t9, 0, 0, 0, 0, t12, 0, -t10, -t60, -t58 * t78 + t1, -t79 + (-t4 + t78) * t55, (-t10 * t55 - t12 * t58) * pkin(4), (t1 * t58 + t2 * t55) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t87, 0, t22, -t23, 0, 0, 0, 0, t18, 0, -t16, 0, t5, -t6, (-t16 * t55 - t18 * t58) * pkin(4), (t5 * t58 + t55 * t6) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t80, -0.2e1 * t81, 0, (t55 ^ 2 + t58 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, -t60, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t80, -t81, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
