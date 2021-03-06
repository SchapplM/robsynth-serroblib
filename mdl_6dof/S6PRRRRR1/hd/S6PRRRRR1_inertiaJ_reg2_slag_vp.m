% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:25:57
% EndTime: 2019-05-05 10:26:04
% DurationCPUTime: 2.07s
% Computational Cost: add. (1846->164), mult. (3846->284), div. (0->0), fcn. (4714->12), ass. (0->97)
t76 = cos(qJ(4));
t60 = t76 * pkin(3);
t55 = t60 + pkin(4);
t70 = sin(qJ(5));
t71 = sin(qJ(4));
t109 = t71 * pkin(3);
t75 = cos(qJ(5));
t92 = t75 * t109;
t39 = t55 * t70 + t92;
t37 = pkin(11) + t39;
t69 = sin(qJ(6));
t63 = t69 ^ 2;
t74 = cos(qJ(6));
t65 = t74 ^ 2;
t95 = t63 + t65;
t98 = t95 * t37;
t110 = t70 * pkin(4);
t53 = pkin(11) + t110;
t119 = t95 * t53;
t112 = -pkin(9) - pkin(8);
t72 = sin(qJ(3));
t90 = t112 * t72;
t77 = cos(qJ(3));
t91 = t112 * t77;
t30 = t71 * t90 - t76 * t91;
t43 = -t71 * t72 + t76 * t77;
t19 = t43 * pkin(10) + t30;
t29 = t71 * t91 + t76 * t90;
t44 = t71 * t77 + t72 * t76;
t83 = -t44 * pkin(10) + t29;
t10 = t19 * t70 - t75 * t83;
t118 = t10 ^ 2;
t67 = sin(pkin(6));
t73 = sin(qJ(2));
t103 = t67 * t73;
t68 = cos(pkin(6));
t40 = -t72 * t103 + t68 * t77;
t41 = t77 * t103 + t68 * t72;
t20 = t40 * t76 - t41 * t71;
t21 = t40 * t71 + t41 * t76;
t13 = -t75 * t20 + t21 * t70;
t117 = t13 ^ 2;
t26 = -t75 * t43 + t44 * t70;
t116 = t26 ^ 2;
t56 = -pkin(3) * t77 - pkin(2);
t34 = -pkin(4) * t43 + t56;
t115 = 0.2e1 * t34;
t114 = 0.2e1 * t44;
t113 = 0.2e1 * t77;
t111 = pkin(5) * t69;
t108 = t10 * t13;
t107 = t10 * t74;
t106 = t13 * t74;
t88 = t70 * t109 - t75 * t55;
t36 = -pkin(5) + t88;
t105 = t36 * t74;
t59 = t75 * pkin(4);
t54 = -t59 - pkin(5);
t104 = t54 * t74;
t78 = cos(qJ(2));
t102 = t67 * t78;
t28 = t43 * t70 + t44 * t75;
t101 = t69 * t28;
t100 = t69 * t74;
t99 = t74 * t28;
t96 = t95 * pkin(11);
t64 = t72 ^ 2;
t66 = t77 ^ 2;
t94 = t64 + t66;
t93 = -0.2e1 * t28 * t26;
t87 = -pkin(5) * t28 - pkin(11) * t26;
t12 = t75 * t19 + t70 * t83;
t9 = pkin(5) * t26 - pkin(11) * t28 + t34;
t3 = -t12 * t69 + t74 * t9;
t4 = t12 * t74 + t69 * t9;
t1 = -t3 * t69 + t4 * t74;
t15 = t20 * t70 + t21 * t75;
t5 = -t74 * t102 - t15 * t69;
t6 = -t69 * t102 + t15 * t74;
t2 = -t5 * t69 + t6 * t74;
t86 = -t26 * t37 + t28 * t36;
t85 = -t26 * t53 + t28 * t54;
t84 = -t40 * t72 + t41 * t77;
t62 = t67 ^ 2;
t61 = pkin(5) * t74;
t51 = t62 * t78 ^ 2;
t49 = 0.2e1 * t100;
t48 = t54 * t69;
t33 = t36 * t69;
t25 = t28 ^ 2;
t24 = t74 * t26;
t23 = t69 * t26;
t22 = t69 * t99;
t16 = (-t63 + t65) * t28;
t8 = t13 * t69;
t7 = t10 * t69;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t73 ^ 2 + t68 ^ 2 + t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 ^ 2 + t41 ^ 2 + t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 ^ 2 + t21 ^ 2 + t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t117 + t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t103, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t102, -t72 * t102, t84, pkin(2) * t102 + t84 * pkin(8), 0, 0, 0, 0, 0, 0, t43 * t102, -t44 * t102, -t20 * t44 + t21 * t43, -t56 * t102 + t20 * t29 + t21 * t30, 0, 0, 0, 0, 0, 0, -t26 * t102, -t28 * t102, t13 * t28 - t15 * t26, -t34 * t102 + t12 * t15 + t108, 0, 0, 0, 0, 0, 0, t13 * t101 + t26 * t5, t13 * t99 - t26 * t6 (-t5 * t74 - t6 * t69) * t28, t3 * t5 + t4 * t6 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t64, t72 * t113, 0, t66, 0, 0, pkin(2) * t113, -0.2e1 * pkin(2) * t72, 0.2e1 * t94 * pkin(8), t94 * pkin(8) ^ 2 + pkin(2) ^ 2, t44 ^ 2, t43 * t114, 0, t43 ^ 2, 0, 0, -0.2e1 * t56 * t43, t56 * t114, -0.2e1 * t29 * t44 + 0.2e1 * t30 * t43, t29 ^ 2 + t30 ^ 2 + t56 ^ 2, t25, t93, 0, t116, 0, 0, t26 * t115, t28 * t115, 0.2e1 * t10 * t28 - 0.2e1 * t12 * t26, t12 ^ 2 + t34 ^ 2 + t118, t65 * t25, -0.2e1 * t25 * t100, 0.2e1 * t26 * t99, t63 * t25, t69 * t93, t116, 0.2e1 * t10 * t101 + 0.2e1 * t26 * t3, 0.2e1 * t10 * t99 - 0.2e1 * t26 * t4, 0.2e1 * (-t3 * t74 - t4 * t69) * t28, t3 ^ 2 + t4 ^ 2 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t21, 0 (t20 * t76 + t21 * t71) * pkin(3), 0, 0, 0, 0, 0, 0, -t13, -t15, 0, t13 * t88 + t15 * t39, 0, 0, 0, 0, 0, 0, -t106, t8, t2, t13 * t36 + t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t77, 0, -t72 * pkin(8), -t77 * pkin(8), 0, 0, 0, 0, t44, 0, t43, 0, t29, -t30 (t43 * t71 - t44 * t76) * pkin(3) (t29 * t76 + t30 * t71) * pkin(3), 0, 0, t28, 0, -t26, 0, -t10, -t12, -t26 * t39 + t28 * t88, t10 * t88 + t12 * t39, t22, t16, t23, -t22, t24, 0, t86 * t69 - t107, t86 * t74 + t7, t1, t1 * t37 + t10 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t109, 0 (t71 ^ 2 + t76 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t88, -0.2e1 * t39, 0, t39 ^ 2 + t88 ^ 2, t63, t49, 0, t65, 0, 0, -0.2e1 * t105, 0.2e1 * t33, 0.2e1 * t98, t95 * t37 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t21, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0 (-t13 * t75 + t15 * t70) * pkin(4), 0, 0, 0, 0, 0, 0, -t106, t8, t2, t13 * t54 + t2 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t43, 0, t29, -t30, 0, 0, 0, 0, t28, 0, -t26, 0, -t10, -t12 (-t26 * t70 - t28 * t75) * pkin(4) (-t10 * t75 + t12 * t70) * pkin(4), t22, t16, t23, -t22, t24, 0, t85 * t69 - t107, t85 * t74 + t7, t1, t1 * t53 + t10 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t60, -t109, 0, 0, 0, 0, 0, 0, 0, 1, t59 - t88, -t92 + (-pkin(4) - t55) * t70, 0 (t39 * t70 - t75 * t88) * pkin(4), t63, t49, 0, t65, 0, 0 (-t36 - t54) * t74, t48 + t33, t119 + t98, t119 * t37 + t36 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t110, 0 (t70 ^ 2 + t75 ^ 2) * pkin(4) ^ 2, t63, t49, 0, t65, 0, 0, -0.2e1 * t104, 0.2e1 * t48, 0.2e1 * t119, t95 * t53 ^ 2 + t54 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t106, t8, t2, -pkin(5) * t13 + t2 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, 0, -t10, -t12, 0, 0, t22, t16, t23, -t22, t24, 0, t87 * t69 - t107, t87 * t74 + t7, t1, -pkin(5) * t10 + t1 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t88, -t39, 0, 0, t63, t49, 0, t65, 0, 0, t61 - t105, t33 - t111, t96 + t98, -pkin(5) * t36 + pkin(11) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t110, 0, 0, t63, t49, 0, t65, 0, 0, t61 - t104, t48 - t111, t96 + t119, -pkin(5) * t54 + pkin(11) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t63, t49, 0, t65, 0, 0, 0.2e1 * t61, -0.2e1 * t111, 0.2e1 * t96, t95 * pkin(11) ^ 2 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, -t101, t26, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t74, 0, -t69 * t37, -t74 * t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t74, 0, -t69 * t53, -t74 * t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t74, 0, -t69 * pkin(11), -t74 * pkin(11), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t11;
