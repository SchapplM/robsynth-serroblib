% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:38
% EndTime: 2019-12-31 21:24:42
% DurationCPUTime: 0.92s
% Computational Cost: add. (1187->131), mult. (2470->270), div. (0->0), fcn. (2700->8), ass. (0->74)
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t59 = sin(qJ(3));
t61 = cos(qJ(3));
t36 = t56 * t59 - t57 * t61;
t49 = -t61 * pkin(3) - pkin(2);
t29 = t36 * pkin(4) + t49;
t89 = 0.2e1 * t29;
t88 = 0.2e1 * t49;
t60 = sin(qJ(2));
t87 = -0.2e1 * t60;
t62 = cos(qJ(2));
t86 = -0.2e1 * t62;
t85 = 0.2e1 * t62;
t84 = pkin(2) * t61;
t83 = pkin(6) * t59;
t53 = t60 ^ 2;
t82 = t53 * pkin(6);
t81 = t56 * pkin(3);
t80 = t57 * pkin(3);
t51 = t60 * pkin(6);
t79 = t62 * pkin(3);
t78 = cos(qJ(5));
t77 = t59 * t60;
t76 = t59 * t61;
t75 = t59 * t62;
t74 = t61 * t60;
t73 = t61 * t62;
t72 = -qJ(4) - pkin(7);
t42 = -t62 * pkin(2) - t60 * pkin(7) - pkin(1);
t39 = t61 * t42;
t70 = qJ(4) * t60;
t20 = -t61 * t70 + t39 + (-pkin(3) - t83) * t62;
t68 = pkin(6) * t73;
t23 = t68 + (t42 - t70) * t59;
t9 = t56 * t20 + t57 * t23;
t41 = pkin(3) * t77 + t51;
t52 = t59 ^ 2;
t54 = t61 ^ 2;
t71 = t52 + t54;
t69 = t60 * t85;
t67 = t59 * t74;
t32 = -t56 * t77 + t57 * t74;
t8 = t57 * t20 - t56 * t23;
t4 = -t62 * pkin(4) - t32 * pkin(8) + t8;
t58 = sin(qJ(5));
t38 = t56 * t61 + t57 * t59;
t30 = t38 * t60;
t7 = -t30 * pkin(8) + t9;
t1 = t78 * t4 - t58 * t7;
t43 = t72 * t59;
t44 = t72 * t61;
t24 = t57 * t43 + t56 * t44;
t27 = -pkin(6) * t75 + t39;
t28 = t59 * t42 + t68;
t66 = -t27 * t59 + t28 * t61;
t25 = t56 * t43 - t57 * t44;
t2 = t58 * t4 + t78 * t7;
t64 = pkin(6) ^ 2;
t55 = t62 ^ 2;
t50 = t53 * t64;
t48 = pkin(4) + t80;
t34 = t58 * t48 + t78 * t81;
t33 = t78 * t48 - t58 * t81;
t21 = t30 * pkin(4) + t41;
t19 = -t58 * t36 + t78 * t38;
t17 = t78 * t36 + t58 * t38;
t14 = -t36 * pkin(8) + t25;
t13 = -t38 * pkin(8) + t24;
t12 = -t58 * t30 + t78 * t32;
t10 = t78 * t30 + t58 * t32;
t6 = t58 * t13 + t78 * t14;
t5 = t78 * t13 - t58 * t14;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t53, t69, 0, t55, 0, 0, pkin(1) * t85, pkin(1) * t87, 0.2e1 * (t53 + t55) * pkin(6), pkin(1) ^ 2 + t55 * t64 + t50, t54 * t53, -0.2e1 * t53 * t76, t73 * t87, t52 * t53, t59 * t69, t55, -0.2e1 * t27 * t62 + 0.2e1 * t59 * t82, 0.2e1 * t28 * t62 + 0.2e1 * t61 * t82, 0.2e1 * (-t27 * t61 - t28 * t59) * t60, t27 ^ 2 + t28 ^ 2 + t50, t32 ^ 2, -0.2e1 * t32 * t30, t32 * t86, t30 ^ 2, -t30 * t86, t55, 0.2e1 * t41 * t30 - 0.2e1 * t8 * t62, 0.2e1 * t41 * t32 + 0.2e1 * t9 * t62, -0.2e1 * t9 * t30 - 0.2e1 * t8 * t32, t41 ^ 2 + t8 ^ 2 + t9 ^ 2, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t86, t10 ^ 2, t10 * t85, t55, -0.2e1 * t1 * t62 + 0.2e1 * t21 * t10, 0.2e1 * t21 * t12 + 0.2e1 * t2 * t62, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t62, 0, -t51, -t62 * pkin(6), 0, 0, t67, (-t52 + t54) * t60, -t75, -t67, -t73, 0, -pkin(6) * t74 + (-pkin(2) * t60 + pkin(7) * t62) * t59, pkin(7) * t73 + (t83 - t84) * t60, t66, -pkin(2) * t51 + pkin(7) * t66, t32 * t38, -t38 * t30 - t32 * t36, -t62 * t38, t30 * t36, t62 * t36, 0, -t24 * t62 + t49 * t30 + t41 * t36, t25 * t62 + t49 * t32 + t41 * t38, -t24 * t32 - t25 * t30 - t9 * t36 - t8 * t38, t8 * t24 + t9 * t25 + t41 * t49, t12 * t19, -t19 * t10 - t12 * t17, -t19 * t62, t10 * t17, t17 * t62, 0, t29 * t10 + t21 * t17 - t5 * t62, t29 * t12 + t21 * t19 + t6 * t62, -t1 * t19 - t6 * t10 - t5 * t12 - t2 * t17, t1 * t5 + t2 * t6 + t21 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t52, 0.2e1 * t76, 0, t54, 0, 0, 0.2e1 * t84, -0.2e1 * pkin(2) * t59, 0.2e1 * t71 * pkin(7), t71 * pkin(7) ^ 2 + pkin(2) ^ 2, t38 ^ 2, -0.2e1 * t38 * t36, 0, t36 ^ 2, 0, 0, t36 * t88, t38 * t88, -0.2e1 * t24 * t38 - 0.2e1 * t25 * t36, t24 ^ 2 + t25 ^ 2 + t49 ^ 2, t19 ^ 2, -0.2e1 * t19 * t17, 0, t17 ^ 2, 0, 0, t17 * t89, t19 * t89, -0.2e1 * t6 * t17 - 0.2e1 * t5 * t19, t29 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t77, -t62, t27, -t28, 0, 0, 0, 0, t32, 0, -t30, -t62, -t57 * t79 + t8, t56 * t79 - t9, (-t30 * t56 - t32 * t57) * pkin(3), (t56 * t9 + t57 * t8) * pkin(3), 0, 0, t12, 0, -t10, -t62, -t33 * t62 + t1, t34 * t62 - t2, -t34 * t10 - t33 * t12, t1 * t33 + t2 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t61, 0, -t59 * pkin(7), -t61 * pkin(7), 0, 0, 0, 0, t38, 0, -t36, 0, t24, -t25, (-t36 * t56 - t38 * t57) * pkin(3), (t24 * t57 + t25 * t56) * pkin(3), 0, 0, t19, 0, -t17, 0, t5, -t6, -t34 * t17 - t33 * t19, t5 * t33 + t6 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t80, -0.2e1 * t81, 0, (t56 ^ 2 + t57 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t34, 0, t33 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, t41, 0, 0, 0, 0, 0, 0, t10, t12, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t38, 0, t49, 0, 0, 0, 0, 0, 0, t17, t19, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, -t62, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t33, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
