% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:33
% DurationCPUTime: 1.42s
% Computational Cost: add. (907->180), mult. (2710->342), div. (0->0), fcn. (2313->8), ass. (0->125)
t69 = sin(qJ(3));
t140 = -0.4e1 * t69;
t68 = sin(qJ(4));
t62 = t68 ^ 2;
t71 = cos(qJ(4));
t64 = t71 ^ 2;
t121 = t62 - t64;
t72 = cos(qJ(3));
t139 = t69 * qJD(5) + (pkin(7) * qJD(4) + qJ(5) * qJD(3)) * t72;
t138 = 0.2e1 * qJD(4);
t137 = pkin(7) * t68;
t136 = t71 * pkin(4);
t66 = sin(pkin(5));
t70 = sin(qJ(2));
t129 = t66 * t70;
t67 = cos(pkin(5));
t34 = t72 * t129 + t67 * t69;
t73 = cos(qJ(2));
t117 = qJD(2) * t73;
t96 = t66 * t117;
t17 = t34 * qJD(3) + t69 * t96;
t135 = t17 * t69;
t110 = t72 * qJD(3);
t103 = pkin(7) * t110;
t61 = qJD(4) * t71;
t37 = t68 * t110 + t69 * t61;
t27 = t37 * pkin(4) + t103;
t134 = t27 * t68;
t133 = t27 * t71;
t33 = t69 * t129 - t67 * t72;
t132 = t33 * t17;
t124 = -qJ(5) - pkin(8);
t48 = t124 * t68;
t131 = t48 * t69;
t49 = t124 * t71;
t130 = t49 * t69;
t128 = t66 * t73;
t127 = t68 * t72;
t126 = t69 * t71;
t125 = t71 * t72;
t86 = -t72 * pkin(3) - t69 * pkin(8);
t81 = -pkin(2) + t86;
t77 = qJD(4) * t81;
t85 = pkin(3) * t69 - pkin(8) * t72;
t78 = t85 * qJD(3);
t123 = -t68 * t78 - t71 * t77;
t59 = t69 * qJD(3);
t94 = t68 * t59;
t122 = pkin(7) * t94 + t71 * t78;
t56 = pkin(7) * t125;
t29 = t68 * t81 + t56;
t63 = t69 ^ 2;
t120 = -t72 ^ 2 + t63;
t119 = qJ(5) * t69;
t118 = qJD(2) * t70;
t116 = qJD(3) * t68;
t115 = qJD(3) * t71;
t114 = qJD(4) * t68;
t113 = qJD(4) * t69;
t112 = qJD(4) * t72;
t109 = pkin(7) * t127;
t108 = -0.2e1 * pkin(2) * qJD(3);
t107 = t71 * t128;
t106 = -0.2e1 * t114;
t105 = pkin(4) * t114;
t104 = pkin(4) * t59;
t102 = t68 * t113;
t101 = t68 * t112;
t100 = t71 * t112;
t99 = t33 * t114;
t45 = (pkin(4) * t68 + pkin(7)) * t69;
t98 = t45 * t114;
t97 = t66 * t118;
t95 = t68 * t61;
t93 = t69 * t110;
t92 = t71 * t110;
t58 = -pkin(3) - t136;
t91 = -t58 + t136;
t90 = t120 * qJD(3);
t89 = 0.2e1 * t93;
t88 = t68 * t92;
t87 = t63 * t95;
t84 = pkin(4) * t62 + t58 * t71;
t18 = -t34 * t68 - t107;
t79 = t68 * t128 - t34 * t71;
t83 = -t18 * t71 + t68 * t79;
t43 = t71 * t81;
t28 = t43 - t109;
t82 = -t28 * t71 - t29 * t68;
t10 = t17 * t68 + t33 * t61;
t11 = -t17 * t71 + t99;
t35 = t92 - t102;
t36 = t71 * t59 + t101;
t16 = -qJD(3) * t33 + t72 * t96;
t7 = t79 * qJD(4) - t16 * t68 + t71 * t97;
t8 = -qJD(4) * t107 - t34 * t114 + t16 * t71 + t68 * t97;
t3 = t83 * qJD(4) - t7 * t68 + t8 * t71;
t12 = t36 * pkin(7) + t123;
t13 = -t29 * qJD(4) + t122;
t76 = t82 * qJD(4) - t12 * t71 - t13 * t68;
t75 = t16 * t72 + t135 + (t33 * t72 - t34 * t69) * qJD(3);
t74 = qJ(5) * t102 - t139 * t71 - t68 * t77 + t122;
t54 = -0.2e1 * t93;
t53 = -0.2e1 * t95;
t52 = 0.2e1 * t95;
t44 = -0.2e1 * t121 * qJD(4);
t38 = t94 - t100;
t32 = -t68 * qJD(5) + t124 * t61;
t31 = -t71 * qJD(5) - t124 * t114;
t26 = 0.2e1 * t64 * t93 - 0.2e1 * t87;
t25 = 0.2e1 * t62 * t93 + 0.2e1 * t87;
t24 = t121 * t113 - t88;
t23 = -t121 * t110 + t95 * t140;
t22 = 0.2e1 * t69 * t100 - 0.2e1 * t68 * t90;
t21 = 0.2e1 * t69 * t101 + 0.2e1 * t120 * t115;
t20 = -t68 * t119 + t29;
t15 = t121 * t63 * t138 + t88 * t140;
t14 = -t71 * t119 + t43 + (-pkin(4) - t137) * t72;
t9 = (pkin(7) * qJD(3) + qJ(5) * qJD(4)) * t126 + t139 * t68 + t123;
t6 = t74 + t104;
t5 = (t33 * t115 + t8) * t72 + (qJD(3) * t79 - t11) * t69;
t4 = (t33 * t116 - t7) * t72 + (qJD(3) * t18 + t10) * t69;
t2 = 0.2e1 * t18 * t7 - 0.2e1 * t79 * t8 + 0.2e1 * t132;
t1 = t83 * t110 + (-t68 * t8 - t7 * t71 + (t18 * t68 + t71 * t79) * qJD(4)) * t69;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t66 ^ 2 * t70 * t117 + 0.2e1 * t34 * t16 + 0.2e1 * t132, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, 0, 0, 0, 0, 0, 0, 0, 0, (-t72 * t118 - t73 * t59) * t66, (-t73 * t110 + t69 * t118) * t66, t75, -pkin(2) * t97 + t75 * pkin(7), 0, 0, 0, 0, 0, 0, t4, t5, t1, t79 * t12 + t18 * t13 + t7 * t28 + t8 * t29 + (t33 * t110 + t135) * pkin(7), 0, 0, 0, 0, 0, 0, t4, t5, t1, t7 * t14 + t17 * t45 + t18 * t6 + t8 * t20 + t33 * t27 + t79 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -0.2e1 * t90, 0, t54, 0, 0, t69 * t108, t72 * t108, 0, 0, t26, t15, t21, t25, t22, t54, 0.2e1 * t28 * t59 - 0.2e1 * t13 * t72 + 0.2e1 * (t63 * t61 + t68 * t89) * pkin(7), -0.2e1 * t29 * t59 - 0.2e1 * t12 * t72 + 0.2e1 * (-t63 * t114 + t71 * t89) * pkin(7), 0.2e1 * t82 * t110 + 0.2e1 * (t12 * t68 - t13 * t71 + (t28 * t68 - t29 * t71) * qJD(4)) * t69, 0.2e1 * pkin(7) ^ 2 * t93 - 0.2e1 * t29 * t12 + 0.2e1 * t28 * t13, t26, t15, t21, t25, t22, t54, 0.2e1 * (t45 * t116 - t6) * t72 + 0.2e1 * (qJD(3) * t14 + t45 * t61 + t134) * t69, 0.2e1 * (t45 * t115 - t9) * t72 + 0.2e1 * (-qJD(3) * t20 + t133 - t98) * t69, 0.2e1 * (-t14 * t71 - t20 * t68) * t110 + 0.2e1 * (-t6 * t71 + t68 * t9 + (t14 * t68 - t20 * t71) * qJD(4)) * t69, 0.2e1 * t14 * t6 - 0.2e1 * t20 * t9 + 0.2e1 * t45 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t16, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, t3, -t17 * pkin(3) + t3 * pkin(8), 0, 0, 0, 0, 0, 0, t11, t10, t3, pkin(4) * t99 + t17 * t58 + t18 * t32 + t31 * t79 + t7 * t48 - t8 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, -t59, 0, -t103, pkin(7) * t59, 0, 0, -t24, t23, t38, t24, t36, 0, (pkin(8) * t125 + (-pkin(3) * t71 + t137) * t69) * qJD(4) + (t86 * t68 - t56) * qJD(3), (pkin(7) * t126 + t85 * t68) * qJD(4) + (t86 * t71 + t109) * qJD(3), t76, -pkin(3) * t103 + t76 * pkin(8), -t24, t23, t38, t24, t36, 0, -t133 - t32 * t72 + (t58 * t127 + t131) * qJD(3) + (t45 * t68 + t84 * t69) * qJD(4), t134 - t31 * t72 + (t58 * t125 + t130) * qJD(3) + (t91 * t69 * t68 + t45 * t71) * qJD(4), (-t48 * t110 - t32 * t69 - t9 + (-t14 + t130) * qJD(4)) * t71 + (t49 * t110 + t31 * t69 - t6 + (-t20 + t131) * qJD(4)) * t68, pkin(4) * t98 + t14 * t32 - t20 * t31 + t27 * t58 + t6 * t48 + t9 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t44, 0, t53, 0, 0, pkin(3) * t106, -0.2e1 * pkin(3) * t61, 0, 0, t52, t44, 0, t53, 0, 0, t91 * t106, t84 * t138, -0.2e1 * t31 * t71 - 0.2e1 * t32 * t68 + 0.2e1 * (-t48 * t71 + t49 * t68) * qJD(4), 0.2e1 * t105 * t58 + 0.2e1 * t49 * t31 + 0.2e1 * t48 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, t7 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t37, t59, t13, t12, 0, 0, 0, 0, t35, 0, -t37, t59, t74 + 0.2e1 * t104, t9, -t35 * pkin(4), t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t114, 0, -pkin(8) * t61, pkin(8) * t114, 0, 0, 0, 0, t61, 0, -t114, 0, t32, t31, -pkin(4) * t61, t32 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t35, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t61, 0, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t19;
