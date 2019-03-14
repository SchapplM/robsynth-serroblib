% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t50 = cos(pkin(6));
t53 = sin(qJ(4));
t57 = cos(qJ(4));
t49 = sin(pkin(6));
t58 = cos(qJ(2));
t81 = t49 * t58;
t22 = t50 * t53 + t57 * t81;
t21 = t22 ^ 2;
t56 = cos(qJ(5));
t41 = -t56 * pkin(5) - pkin(4);
t93 = 0.2e1 * t41;
t92 = 0.2e1 * t53;
t91 = 0.2e1 * t56;
t90 = 2 * qJ(3);
t89 = -pkin(10) - pkin(9);
t51 = sin(qJ(6));
t88 = t51 * pkin(5);
t55 = cos(qJ(6));
t87 = t55 * pkin(5);
t31 = t53 * pkin(4) - t57 * pkin(9) + qJ(3);
t52 = sin(qJ(5));
t59 = -pkin(2) - pkin(8);
t74 = t56 * t59;
t68 = t53 * t74;
t8 = t68 + (-pkin(10) * t57 + t31) * t52;
t86 = t55 * t8;
t85 = t57 * pkin(4);
t84 = t22 * t57;
t29 = t51 * t56 + t55 * t52;
t83 = t29 * t53;
t54 = sin(qJ(2));
t82 = t49 * t54;
t80 = t52 * t53;
t79 = t52 * t56;
t78 = t52 * t57;
t77 = t52 * t59;
t76 = t53 * t59;
t75 = t56 * t53;
t40 = t56 * t57;
t73 = t57 * t29;
t72 = t57 * t53;
t71 = t57 * t59;
t44 = t52 ^ 2;
t46 = t56 ^ 2;
t70 = t44 + t46;
t45 = t53 ^ 2;
t47 = t57 ^ 2;
t37 = t45 + t47;
t69 = -0.2e1 * t72;
t67 = t52 * t40;
t25 = t56 * t31;
t7 = -pkin(10) * t40 + t25 + (pkin(5) - t77) * t53;
t1 = -t51 * t8 + t55 * t7;
t66 = t70 * t53;
t65 = -pkin(9) * t53 - t85;
t24 = t50 * t57 - t53 * t81;
t10 = t24 * t56 + t52 * t82;
t9 = -t24 * t52 + t56 * t82;
t64 = t10 * t56 - t9 * t52;
t13 = -t52 * t76 + t25;
t14 = t52 * t31 + t68;
t63 = -t13 * t52 + t14 * t56;
t5 = t24 * t53 - t84;
t60 = qJ(3) ^ 2;
t48 = t59 ^ 2;
t43 = t49 ^ 2;
t42 = t47 * t48;
t38 = t43 * t54 ^ 2;
t36 = qJ(3) * t82;
t33 = t89 * t56;
t32 = t89 * t52;
t30 = t37 * t59;
t27 = t51 * t52 - t55 * t56;
t26 = (pkin(5) * t52 - t59) * t57;
t20 = t43 * t58 ^ 2 + t50 ^ 2 + t38;
t19 = t55 * t40 - t51 * t78;
t18 = -t51 * t80 + t55 * t75;
t12 = t51 * t32 - t55 * t33;
t11 = t55 * t32 + t51 * t33;
t4 = t55 * t10 + t51 * t9;
t3 = -t51 * t10 + t55 * t9;
t2 = t51 * t7 + t86;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 ^ 2 + t21 + t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t9 ^ 2 + t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t82, pkin(2) * t81 + t36, 0, 0, 0, 0, 0, 0, t53 * t82, t57 * t82, -t5, t5 * t59 + t36, 0, 0, 0, 0, 0, 0, t22 * t78 + t9 * t53, -t10 * t53 + t22 * t40 (-t10 * t52 - t56 * t9) * t57, t10 * t14 + t9 * t13 - t22 * t71, 0, 0, 0, 0, 0, 0, t22 * t73 + t3 * t53, t22 * t19 - t4 * t53, -t3 * t19 - t4 * t73, t3 * t1 + t4 * t2 + t22 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t90, pkin(2) ^ 2 + t60, t47, t69, 0, t45, 0, 0, t53 * t90, t57 * t90, -0.2e1 * t30, t45 * t48 + t42 + t60, t46 * t47, -0.2e1 * t47 * t79, t72 * t91, t44 * t47, t52 * t69, t45, 0.2e1 * t13 * t53 - 0.2e1 * t47 * t77, -0.2e1 * t14 * t53 - 0.2e1 * t47 * t74, 0.2e1 * (-t13 * t56 - t14 * t52) * t57, t13 ^ 2 + t14 ^ 2 + t42, t19 ^ 2, -0.2e1 * t19 * t73, t19 * t92, t73 ^ 2, -t73 * t92, t45, 0.2e1 * t1 * t53 + 0.2e1 * t26 * t73, 0.2e1 * t26 * t19 - 0.2e1 * t2 * t53, -0.2e1 * t1 * t19 - 0.2e1 * t2 * t73, t1 ^ 2 + t2 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t53 - t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t18 - t3 * t83 - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t37, t30, 0, 0, 0, 0, 0, 0, -t37 * t52, -t37 * t56, 0, t47 * t59 + t63 * t53, 0, 0, 0, 0, 0, 0, -t53 * t83 - t57 * t73, -t18 * t53 - t57 * t19, -t18 * t73 + t19 * t83, -t1 * t83 + t2 * t18 - t26 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t45 + t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 + t83 ^ 2 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t24, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t56, t22 * t52, t64, -t22 * pkin(4) + t64 * pkin(9), 0, 0, 0, 0, 0, 0, t22 * t27, t22 * t29, -t4 * t27 - t3 * t29, t3 * t11 + t4 * t12 + t22 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t53, 0, t71, -t76, 0, 0, t67 (-t44 + t46) * t57, t80, -t67, t75, 0, t65 * t52 + t56 * t71, -t52 * t71 + t65 * t56, t63, pkin(4) * t71 + t63 * pkin(9), t19 * t29, -t19 * t27 - t29 * t73, t83, t73 * t27, -t27 * t53, 0, t11 * t53 + t26 * t27 + t41 * t73, -t12 * t53 + t41 * t19 + t26 * t29, -t1 * t29 - t11 * t19 - t12 * t73 - t2 * t27, t1 * t11 + t2 * t12 + t26 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t53, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t78, t66, pkin(9) * t66 + t85, 0, 0, 0, 0, 0, 0, -t57 * t27, -t73, -t18 * t27 + t29 * t83, -t11 * t83 + t18 * t12 - t57 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t44, 0.2e1 * t79, 0, t46, 0, 0, pkin(4) * t91, -0.2e1 * pkin(4) * t52, 0.2e1 * t70 * pkin(9), t70 * pkin(9) ^ 2 + pkin(4) ^ 2, t29 ^ 2, -0.2e1 * t29 * t27, 0, t27 ^ 2, 0, 0, t27 * t93, t29 * t93, -0.2e1 * t11 * t29 - 0.2e1 * t12 * t27, t11 ^ 2 + t12 ^ 2 + t41 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0 (t3 * t55 + t4 * t51) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t78, t53, t13, -t14, 0, 0, 0, 0, t19, 0, -t73, t53, t53 * t87 + t1, -t86 + (-t53 * pkin(5) - t7) * t51 (-t19 * t55 - t51 * t73) * pkin(5) (t1 * t55 + t2 * t51) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t75, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t18, 0 (t18 * t51 - t55 * t83) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t56, 0, -t52 * pkin(9), -t56 * pkin(9), 0, 0, 0, 0, t29, 0, -t27, 0, t11, -t12 (-t27 * t51 - t29 * t55) * pkin(5) (t11 * t55 + t12 * t51) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t87, -0.2e1 * t88, 0 (t51 ^ 2 + t55 ^ 2) * pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t73, t53, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t27, 0, t11, -t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t87, -t88, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t6;