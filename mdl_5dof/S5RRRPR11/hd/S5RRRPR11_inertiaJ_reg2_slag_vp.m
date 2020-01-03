% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t53 = sin(qJ(3));
t48 = t53 ^ 2;
t56 = cos(qJ(3));
t50 = t56 ^ 2;
t90 = t48 + t50;
t54 = sin(qJ(2));
t41 = t56 * t54;
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t73 = t55 * t53;
t14 = t52 * t41 - t54 * t73;
t89 = -0.2e1 * t14;
t69 = t53 * qJ(4);
t83 = pkin(3) + pkin(4);
t17 = t83 * t56 + pkin(2) + t69;
t88 = 0.2e1 * t17;
t62 = -t56 * pkin(3) - t69;
t29 = -pkin(2) + t62;
t87 = -0.2e1 * t29;
t86 = -0.2e1 * t54;
t85 = 0.2e1 * t54;
t57 = cos(qJ(2));
t84 = 0.2e1 * t57;
t82 = pkin(2) * t53;
t81 = pkin(2) * t56;
t49 = t54 ^ 2;
t80 = t49 * pkin(6);
t79 = t53 * pkin(7);
t78 = t54 * pkin(6);
t77 = t53 * t54;
t76 = t53 * t56;
t75 = t53 * t57;
t74 = t54 * t57;
t72 = t56 * t57;
t30 = -t57 * pkin(2) - t54 * pkin(7) - pkin(1);
t71 = pkin(6) * t75 - t56 * t30;
t12 = pkin(6) * t72 + t53 * t30;
t70 = t90 * pkin(7) ^ 2;
t68 = t56 * qJ(4);
t67 = t57 * qJ(4);
t66 = t53 * t74;
t65 = t49 * t76;
t47 = t57 * pkin(3);
t10 = t47 + t71;
t64 = (pkin(7) - pkin(8)) * t53;
t9 = -t67 + t12;
t3 = t57 * pkin(4) - pkin(8) * t41 + t10;
t4 = pkin(8) * t77 + t9;
t1 = t55 * t3 - t52 * t4;
t2 = t52 * t3 + t55 * t4;
t63 = t10 * t53 + t9 * t56;
t61 = -pkin(3) * t53 + t68;
t60 = t12 * t56 + t53 * t71;
t20 = t52 * t53 + t55 * t56;
t59 = pkin(6) ^ 2;
t51 = t57 ^ 2;
t46 = t56 * pkin(7);
t44 = t49 * t59;
t40 = t50 * t49;
t39 = t48 * t49;
t36 = pkin(7) * t75;
t34 = t53 * t41;
t32 = -t56 * pkin(8) + t46;
t31 = t72 * t86;
t28 = t55 * qJ(4) - t52 * t83;
t26 = t52 * qJ(4) + t55 * t83;
t25 = 0.2e1 * t90 * pkin(7);
t23 = -t52 * t56 + t73;
t22 = (t48 - t50) * t54;
t16 = t20 * t54;
t13 = (pkin(6) - t61) * t54;
t8 = t55 * t32 + t52 * t64;
t6 = t52 * t32 - t55 * t64;
t5 = (-t83 * t53 - pkin(6) + t68) * t54;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t49, 0.2e1 * t74, 0, t51, 0, 0, pkin(1) * t84, pkin(1) * t86, 0.2e1 * (t49 + t51) * pkin(6), pkin(1) ^ 2 + t51 * t59 + t44, t40, -0.2e1 * t65, t31, t39, 0.2e1 * t66, t51, 0.2e1 * t53 * t80 + 0.2e1 * t57 * t71, 0.2e1 * t12 * t57 + 0.2e1 * t56 * t80, (-t12 * t53 + t56 * t71) * t85, t12 ^ 2 + t71 ^ 2 + t44, t40, t31, 0.2e1 * t65, t51, -0.2e1 * t66, t39, 0.2e1 * t10 * t57 + 0.2e1 * t13 * t77, (t10 * t56 - t53 * t9) * t85, -0.2e1 * t13 * t41 - 0.2e1 * t9 * t57, t10 ^ 2 + t13 ^ 2 + t9 ^ 2, t16 ^ 2, t16 * t89, t16 * t84, t14 ^ 2, t57 * t89, t51, 0.2e1 * t1 * t57 + 0.2e1 * t5 * t14, 0.2e1 * t5 * t16 - 0.2e1 * t2 * t57, -0.2e1 * t1 * t16 - 0.2e1 * t2 * t14, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t57, 0, -t78, -t57 * pkin(6), 0, 0, t34, -t22, -t75, -t34, -t72, 0, t36 + (-pkin(6) * t56 - t82) * t54, pkin(7) * t72 + (pkin(6) * t53 - t81) * t54, t60, -pkin(2) * t78 + pkin(7) * t60, t34, -t75, t22, 0, t72, -t34, -t13 * t56 + t29 * t77 + t36, t63, -t13 * t53 + (-pkin(7) * t57 - t29 * t54) * t56, pkin(7) * t63 + t13 * t29, t16 * t23, -t23 * t14 - t16 * t20, t23 * t57, t14 * t20, -t20 * t57, 0, t17 * t14 + t5 * t20 - t6 * t57, t17 * t16 + t5 * t23 - t8 * t57, -t1 * t23 - t8 * t14 + t6 * t16 - t2 * t20, -t1 * t6 + t5 * t17 + t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t48, 0.2e1 * t76, 0, t50, 0, 0, 0.2e1 * t81, -0.2e1 * t82, t25, pkin(2) ^ 2 + t70, t48, 0, -0.2e1 * t76, 0, 0, t50, t56 * t87, t25, t53 * t87, t29 ^ 2 + t70, t23 ^ 2, -0.2e1 * t23 * t20, 0, t20 ^ 2, 0, 0, t20 * t88, t23 * t88, -0.2e1 * t8 * t20 + 0.2e1 * t6 * t23, t17 ^ 2 + t6 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t77, -t57, -t71, -t12, 0, 0, 0, t41, 0, -t57, t77, 0, -0.2e1 * t47 - t71, t62 * t54, -0.2e1 * t67 + t12, -t10 * pkin(3) + t9 * qJ(4), 0, 0, -t16, 0, t14, -t57, -t26 * t57 - t1, -t28 * t57 + t2, -t28 * t14 + t26 * t16, -t1 * t26 + t2 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t56, 0, -t79, -t46, 0, 0, 0, t53, 0, 0, -t56, 0, -t79, t61, t46, t61 * pkin(7), 0, 0, -t23, 0, t20, 0, t6, t8, -t28 * t20 + t26 * t23, t6 * t26 + t8 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, 0.2e1 * t28, 0, t26 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t41, 0, t10, 0, 0, 0, 0, 0, 0, t55 * t57, -t52 * t57, -t52 * t14 - t55 * t16, t1 * t55 + t2 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t20 - t55 * t23, t8 * t52 - t6 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t55, t52, 0, -t26 * t55 + t28 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 ^ 2 + t55 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, t57, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t20, 0, -t6, -t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t26, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;
