% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:14:53
% EndTime: 2019-05-06 03:14:59
% DurationCPUTime: 1.79s
% Computational Cost: add. (2994->132), mult. (5643->238), div. (0->0), fcn. (7022->10), ass. (0->87)
t68 = cos(qJ(4));
t54 = t68 * pkin(3);
t50 = t54 + pkin(4);
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t64 = sin(qJ(4));
t99 = t64 * pkin(3);
t85 = t67 * t99;
t35 = t63 * t50 + t85;
t33 = pkin(10) + t35;
t62 = sin(qJ(6));
t58 = t62 ^ 2;
t66 = cos(qJ(6));
t59 = t66 ^ 2;
t87 = t58 + t59;
t91 = t87 * t33;
t100 = t63 * pkin(4);
t48 = pkin(10) + t100;
t114 = t87 * t48;
t92 = pkin(7) + qJ(2);
t98 = cos(qJ(3));
t109 = t92 * t98;
t65 = sin(qJ(3));
t110 = t65 * t92;
t113 = t109 * t68 - t110 * t64;
t112 = -t109 * t64 - t110 * t68;
t60 = sin(pkin(11));
t61 = cos(pkin(11));
t37 = t65 * t60 - t61 * t98;
t39 = t60 * t98 + t65 * t61;
t23 = -t68 * t37 - t64 * t39;
t74 = t64 * t37 - t68 * t39;
t20 = t63 * t23 - t67 * t74;
t111 = -0.2e1 * t20;
t25 = -t109 * t60 - t110 * t61;
t26 = t109 * t61 - t110 * t60;
t11 = t68 * (-t37 * pkin(8) + t26) + t64 * (-t39 * pkin(8) + t25) + t23 * pkin(9);
t12 = t74 * pkin(8) + t112 * t61 - t113 * t60;
t73 = pkin(9) * t74 + t12;
t5 = t63 * t11 - t67 * t73;
t108 = t5 ^ 2;
t18 = -t67 * t23 - t63 * t74;
t107 = t18 ^ 2;
t47 = -t61 * pkin(2) - pkin(1);
t27 = t37 * pkin(3) + t47;
t21 = -t23 * pkin(4) + t27;
t106 = 0.2e1 * t21;
t105 = -0.2e1 * t74;
t104 = 0.2e1 * t47;
t103 = 0.2e1 * t61;
t102 = pkin(5) * t62;
t101 = t5 * t66;
t81 = -t67 * t50 + t63 * t99;
t32 = -pkin(5) + t81;
t97 = t32 * t66;
t53 = t67 * pkin(4);
t49 = -t53 - pkin(5);
t96 = t49 * t66;
t15 = t62 * t18;
t95 = t62 * t20;
t94 = t62 * t66;
t93 = t66 * t20;
t89 = pkin(10) * t87;
t56 = t60 ^ 2;
t57 = t61 ^ 2;
t88 = t56 + t57;
t86 = t18 * t111;
t80 = -pkin(5) * t20 - pkin(10) * t18;
t7 = t67 * t11 + t63 * t73;
t8 = t18 * pkin(5) - t20 * pkin(10) + t21;
t2 = -t62 * t7 + t66 * t8;
t3 = t62 * t8 + t66 * t7;
t79 = t2 * t66 + t3 * t62;
t1 = -t2 * t62 + t3 * t66;
t76 = -t18 * t33 + t20 * t32;
t75 = -t18 * t48 + t20 * t49;
t55 = pkin(5) * t66;
t44 = 0.2e1 * t94;
t43 = t49 * t62;
t30 = t32 * t62;
t17 = t20 ^ 2;
t16 = t66 * t18;
t14 = t62 * t93;
t13 = t23 * pkin(8) + t112 * t60 + t113 * t61;
t9 = (-t58 + t59) * t20;
t4 = t5 * t62;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t56, t60 * t103, 0, t57, 0, 0, pkin(1) * t103, -0.2e1 * pkin(1) * t60, 0.2e1 * t88 * qJ(2), qJ(2) ^ 2 * t88 + pkin(1) ^ 2, t39 ^ 2, -0.2e1 * t39 * t37, 0, t37 ^ 2, 0, 0, t37 * t104, t39 * t104, -0.2e1 * t25 * t39 - 0.2e1 * t26 * t37, t25 ^ 2 + t26 ^ 2 + t47 ^ 2, t74 ^ 2, t23 * t105, 0, t23 ^ 2, 0, 0, -0.2e1 * t27 * t23, t27 * t105, 0.2e1 * t12 * t74 + 0.2e1 * t13 * t23, t12 ^ 2 + t13 ^ 2 + t27 ^ 2, t17, t86, 0, t107, 0, 0, t18 * t106, t20 * t106, -0.2e1 * t7 * t18 + 0.2e1 * t5 * t20, t21 ^ 2 + t7 ^ 2 + t108, t59 * t17, -0.2e1 * t17 * t94, 0.2e1 * t18 * t93, t58 * t17, t62 * t86, t107, 0.2e1 * t2 * t18 + 0.2e1 * t5 * t95, -0.2e1 * t3 * t18 + 0.2e1 * t5 * t93, t79 * t111, t2 ^ 2 + t3 ^ 2 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t60, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t37, t39, 0, t47, 0, 0, 0, 0, 0, 0, -t23, -t74, 0, t27, 0, 0, 0, 0, 0, 0, t18, t20, 0, t21, 0, 0, 0, 0, 0, 0, t16, -t15, -t87 * t20, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t37, 0, t25, -t26, 0, 0, 0, 0, -t74, 0, t23, 0, t12, -t13 (t23 * t64 + t68 * t74) * pkin(3) (t12 * t68 + t13 * t64) * pkin(3), 0, 0, t20, 0, -t18, 0, -t5, -t7, -t35 * t18 + t20 * t81, t7 * t35 + t5 * t81, t14, t9, t15, -t14, t16, 0, t62 * t76 - t101, t66 * t76 + t4, t1, t1 * t33 + t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t54, -0.2e1 * t99, 0 (t64 ^ 2 + t68 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t81, -0.2e1 * t35, 0, t35 ^ 2 + t81 ^ 2, t58, t44, 0, t59, 0, 0, -0.2e1 * t97, 0.2e1 * t30, 0.2e1 * t91, t33 ^ 2 * t87 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, t23, 0, t12, -t13, 0, 0, 0, 0, t20, 0, -t18, 0, -t5, -t7 (-t18 * t63 - t20 * t67) * pkin(4) (-t5 * t67 + t63 * t7) * pkin(4), t14, t9, t15, -t14, t16, 0, t62 * t75 - t101, t66 * t75 + t4, t1, t1 * t48 + t5 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t54, -t99, 0, 0, 0, 0, 0, 0, 0, 1, t53 - t81, -t85 + (-pkin(4) - t50) * t63, 0 (t35 * t63 - t67 * t81) * pkin(4), t58, t44, 0, t59, 0, 0 (-t32 - t49) * t66, t43 + t30, t114 + t91, t114 * t33 + t32 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t53, -0.2e1 * t100, 0 (t63 ^ 2 + t67 ^ 2) * pkin(4) ^ 2, t58, t44, 0, t59, 0, 0, -0.2e1 * t96, 0.2e1 * t43, 0.2e1 * t114, t48 ^ 2 * t87 + t49 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, -t5, -t7, 0, 0, t14, t9, t15, -t14, t16, 0, t62 * t80 - t101, t66 * t80 + t4, t1, -t5 * pkin(5) + pkin(10) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t81, -t35, 0, 0, t58, t44, 0, t59, 0, 0, t55 - t97, t30 - t102, t89 + t91, -t32 * pkin(5) + pkin(10) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t53, -t100, 0, 0, t58, t44, 0, t59, 0, 0, t55 - t96, t43 - t102, t89 + t114, -t49 * pkin(5) + pkin(10) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t58, t44, 0, t59, 0, 0, 0.2e1 * t55, -0.2e1 * t102, 0.2e1 * t89, pkin(10) ^ 2 * t87 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, -t95, t18, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t62, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t66, 0, -t62 * t33, -t66 * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t66, 0, -t62 * t48, -t66 * t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t66, 0, -t62 * pkin(10), -t66 * pkin(10), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t6;
