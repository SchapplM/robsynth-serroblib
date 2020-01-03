% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t53 = cos(pkin(5));
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t51 = sin(pkin(5));
t56 = sin(qJ(2));
t77 = t51 * t56;
t32 = -t53 * t58 + t55 * t77;
t33 = t53 * t55 + t58 * t77;
t50 = sin(pkin(10));
t52 = cos(pkin(10));
t20 = -t50 * t32 + t52 * t33;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t59 = cos(qJ(2));
t76 = t51 * t59;
t14 = t54 * t20 + t57 * t76;
t84 = -0.2e1 * t14;
t83 = -0.2e1 * t33;
t82 = 0.2e1 * t58;
t81 = pkin(1) * t56;
t80 = pkin(1) * t59;
t65 = pkin(7) * t76;
t28 = t65 + (pkin(8) + t81) * t53;
t29 = (-pkin(2) * t59 - pkin(8) * t56 - pkin(1)) * t51;
t17 = t58 * t28 + t55 * t29;
t13 = -t32 * qJ(4) + t17;
t16 = -t55 * t28 + t58 * t29;
t9 = -pkin(3) * t76 - t33 * qJ(4) + t16;
t6 = t52 * t13 + t50 * t9;
t15 = t57 * t20 - t54 * t76;
t79 = t15 * t54;
t47 = t51 ^ 2;
t78 = t47 * t59;
t75 = t53 * t56;
t19 = t52 * t32 + t50 * t33;
t74 = t54 * t19;
t37 = t50 * t55 - t52 * t58;
t73 = t54 * t37;
t38 = t50 * t58 + t52 * t55;
t72 = t54 * t38;
t44 = t50 * pkin(3) + pkin(9);
t71 = t54 * t44;
t70 = t54 * t57;
t69 = t57 * t38;
t68 = t57 * t44;
t67 = -qJ(4) - pkin(8);
t66 = 0.2e1 * t76;
t64 = t55 * t76;
t63 = t58 * t76;
t46 = -t58 * pkin(3) - pkin(2);
t62 = t67 * t55;
t5 = -t50 * t13 + t52 * t9;
t45 = -t52 * pkin(3) - pkin(4);
t61 = -t37 * t44 + t38 * t45;
t41 = pkin(7) * t77;
t27 = t41 + (-pkin(2) - t80) * t53;
t21 = t32 * pkin(3) + t27;
t49 = t57 ^ 2;
t48 = t54 ^ 2;
t40 = t67 * t58;
t36 = t38 ^ 2;
t35 = pkin(1) * t75 + t65;
t34 = t53 * t80 - t41;
t31 = t57 * t37;
t25 = -t52 * t40 + t50 * t62;
t23 = -t50 * t40 - t52 * t62;
t22 = t37 * pkin(4) - t38 * pkin(9) + t46;
t18 = t57 * t19;
t11 = t54 * t22 + t57 * t25;
t10 = t57 * t22 - t54 * t25;
t7 = t19 * pkin(4) - t20 * pkin(9) + t21;
t4 = -pkin(9) * t76 + t6;
t3 = pkin(4) * t76 - t5;
t2 = t57 * t4 + t54 * t7;
t1 = -t54 * t4 + t57 * t7;
t8 = [1, 0, 0, t47 * t56 ^ 2, 0.2e1 * t56 * t78, 0.2e1 * t51 * t75, t53 * t66, t53 ^ 2, 0.2e1 * pkin(1) * t78 + 0.2e1 * t34 * t53, -0.2e1 * t35 * t53 - 0.2e1 * t47 * t81, t33 ^ 2, t32 * t83, t76 * t83, t32 * t66, t47 * t59 ^ 2, -0.2e1 * t16 * t76 + 0.2e1 * t27 * t32, 0.2e1 * t17 * t76 + 0.2e1 * t27 * t33, -0.2e1 * t6 * t19 - 0.2e1 * t5 * t20, t21 ^ 2 + t5 ^ 2 + t6 ^ 2, t15 ^ 2, t15 * t84, 0.2e1 * t15 * t19, t19 * t84, t19 ^ 2, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t14, 0.2e1 * t3 * t15 - 0.2e1 * t2 * t19; 0, 0, 0, 0, 0, t77, t76, t53, t34, -t35, t33 * t55, -t55 * t32 + t33 * t58, -t64, -t63, 0, -pkin(2) * t32 + pkin(8) * t64 - t27 * t58, -pkin(2) * t33 + pkin(8) * t63 + t27 * t55, -t25 * t19 + t23 * t20 - t6 * t37 - t5 * t38, t21 * t46 - t5 * t23 + t6 * t25, t15 * t69, (-t14 * t57 - t79) * t38, t15 * t37 + t19 * t69, -t14 * t37 - t19 * t72, t19 * t37, t1 * t37 + t10 * t19 + t23 * t14 + t3 * t72, -t11 * t19 + t23 * t15 - t2 * t37 + t3 * t69; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t55 ^ 2, t55 * t82, 0, 0, 0, pkin(2) * t82, -0.2e1 * pkin(2) * t55, 0.2e1 * t23 * t38 - 0.2e1 * t25 * t37, t23 ^ 2 + t25 ^ 2 + t46 ^ 2, t49 * t36, -0.2e1 * t36 * t70, 0.2e1 * t37 * t69, -0.2e1 * t37 * t72, t37 ^ 2, 0.2e1 * t10 * t37 + 0.2e1 * t23 * t72, -0.2e1 * t11 * t37 + 0.2e1 * t23 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, -t76, t16, -t17, (-t19 * t50 - t20 * t52) * pkin(3), (t5 * t52 + t50 * t6) * pkin(3), t79, -t54 * t14 + t15 * t57, t74, t18, 0, t45 * t14 - t19 * t71 - t3 * t57, t45 * t15 - t19 * t68 + t3 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t55 * pkin(8), -t58 * pkin(8), (-t37 * t50 - t38 * t52) * pkin(3), (-t23 * t52 + t25 * t50) * pkin(3), t54 * t69, (-t48 + t49) * t38, t73, t31, 0, -t23 * t57 + t54 * t61, t23 * t54 + t57 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t50 ^ 2 + t52 ^ 2) * pkin(3) ^ 2, t48, 0.2e1 * t70, 0, 0, 0, -0.2e1 * t45 * t57, 0.2e1 * t45 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, t18, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, t31, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t19, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t72, t37, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t57, 0, -t71, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
