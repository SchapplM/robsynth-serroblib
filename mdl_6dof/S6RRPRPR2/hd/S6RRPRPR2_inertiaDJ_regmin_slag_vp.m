% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:55
% EndTime: 2019-03-09 10:13:58
% DurationCPUTime: 1.09s
% Computational Cost: add. (2349->162), mult. (5192->280), div. (0->0), fcn. (5155->8), ass. (0->108)
t79 = sin(qJ(6));
t77 = t79 ^ 2;
t82 = cos(qJ(6));
t127 = -t82 ^ 2 + t77;
t107 = t127 * qJD(6);
t85 = 2 * qJD(5);
t135 = pkin(4) + pkin(9);
t133 = cos(qJ(4));
t80 = sin(qJ(4));
t125 = sin(pkin(10));
t126 = cos(pkin(10));
t81 = sin(qJ(2));
t83 = cos(qJ(2));
t92 = -t125 * t83 - t126 * t81;
t109 = t125 * t81;
t110 = t126 * t83;
t93 = t109 - t110;
t37 = -t133 * t92 - t80 * t93;
t89 = qJD(2) * t92;
t70 = qJD(2) * t109;
t98 = qJD(2) * t110 - t70;
t24 = qJD(4) * t37 - t133 * t89 + t80 * t98;
t111 = qJD(4) * t133;
t124 = qJD(4) * t80;
t130 = -qJ(3) - pkin(7);
t108 = qJD(2) * t130;
t58 = t83 * qJD(3) + t108 * t81;
t59 = -t81 * qJD(3) + t108 * t83;
t32 = t125 * t59 + t126 * t58;
t28 = pkin(8) * t89 + t32;
t67 = t130 * t81;
t68 = t130 * t83;
t38 = t125 * t68 + t126 * t67;
t33 = pkin(8) * t92 + t38;
t39 = t125 * t67 - t126 * t68;
t34 = -pkin(8) * t93 + t39;
t31 = -t125 * t58 + t126 * t59;
t87 = -pkin(8) * t98 + t31;
t8 = -t33 * t111 + t124 * t34 - t133 * t28 - t80 * t87;
t5 = -pkin(5) * t24 - t8;
t3 = t5 * t79;
t4 = t5 * t82;
t121 = qJD(6) * t82;
t88 = t133 * t93;
t36 = -t80 * t92 + t88;
t99 = t133 * t34 + t80 * t33;
t16 = -t36 * pkin(5) + t99;
t134 = t16 * t121 + t3;
t132 = t24 * t79;
t131 = t24 * t82;
t112 = pkin(2) * t125;
t71 = t80 * t112;
t72 = pkin(2) * t126 + pkin(3);
t49 = qJD(4) * t71 - t72 * t111;
t44 = -qJD(5) + t49;
t91 = t112 * t133 + t80 * t72;
t55 = qJ(5) + t91;
t129 = t55 * t121 - t44 * t79;
t117 = qJ(5) * qJD(6);
t128 = qJD(5) * t79 + t82 * t117;
t123 = qJD(6) * t16;
t122 = qJD(6) * t79;
t120 = qJD(6) * t135;
t119 = t81 * qJD(2);
t118 = t83 * qJD(2);
t116 = -0.2e1 * pkin(1) * qJD(2);
t115 = t79 * t131;
t74 = pkin(2) * t119;
t114 = t79 * t121;
t113 = -t83 * pkin(2) - pkin(1);
t21 = -t133 * t33 + t80 * t34;
t45 = pkin(3) * t93 + t113;
t86 = -t37 * qJ(5) + t45;
t14 = t135 * t36 + t86;
t15 = t37 * pkin(5) + t21;
t106 = t14 * t82 + t15 * t79;
t105 = t14 * t79 - t15 * t82;
t104 = -qJ(5) * t24 - qJD(5) * t36;
t56 = -t133 * t72 - pkin(4) + t71;
t23 = qJD(4) * t88 - t124 * t92 - t133 * t98 - t80 * t89;
t103 = t121 * t37 - t23 * t79;
t102 = t122 * t37 + t23 * t82;
t101 = t121 * t36 + t132;
t100 = t122 * t36 - t131;
t54 = -pkin(9) + t56;
t97 = qJD(6) * (t36 * t55 - t37 * t54);
t96 = qJD(6) * (qJ(5) * t36 + t135 * t37);
t50 = t91 * qJD(4);
t95 = -t24 * t55 + t36 * t44 + t37 * t50;
t94 = t135 * t23 + t104;
t90 = -t23 * t54 + t95;
t9 = qJD(4) * t99 - t133 * t87 + t80 * t28;
t42 = -pkin(3) * t89 + t74;
t10 = t24 * pkin(4) + t23 * qJ(5) - t37 * qJD(5) + t42;
t76 = qJD(5) * t82;
t69 = -0.2e1 * t114;
t65 = 0.2e1 * t107;
t41 = t44 * t82;
t35 = t36 ^ 2;
t20 = t36 * pkin(4) + t86;
t19 = -0.2e1 * t37 * t23;
t12 = -t107 * t36 + t115;
t11 = -0.4e1 * t114 * t36 - t127 * t24;
t7 = t24 * pkin(9) + t10;
t6 = -t23 * pkin(5) + t9;
t2 = -qJD(6) * t106 + t82 * t6 - t79 * t7;
t1 = qJD(6) * t105 - t79 * t6 - t82 * t7;
t13 = [0, 0, 0, 0.2e1 * t81 * t118, 0.2e1 * (-t81 ^ 2 + t83 ^ 2) * qJD(2), 0, 0, 0, t81 * t116, t83 * t116, 0.2e1 * t31 * t92 - 0.2e1 * t32 * t93 - 0.2e1 * t38 * t98 + 0.2e1 * t39 * t89, 0.2e1 * t113 * t74 + 0.2e1 * t38 * t31 + 0.2e1 * t39 * t32, t19, 0.2e1 * t23 * t36 - 0.2e1 * t24 * t37, 0, 0, 0, 0.2e1 * t24 * t45 + 0.2e1 * t36 * t42, -0.2e1 * t23 * t45 + 0.2e1 * t37 * t42, -0.2e1 * t21 * t23 - 0.2e1 * t24 * t99 + 0.2e1 * t36 * t8 + 0.2e1 * t37 * t9, -0.2e1 * t10 * t36 - 0.2e1 * t20 * t24, -0.2e1 * t10 * t37 + 0.2e1 * t20 * t23, 0.2e1 * t10 * t20 + 0.2e1 * t21 * t9 - 0.2e1 * t8 * t99, 0.2e1 * t24 * t36 * t77 + 0.2e1 * t114 * t35, -0.2e1 * t107 * t35 + 0.4e1 * t115 * t36, 0.2e1 * t103 * t36 + 0.2e1 * t132 * t37, -0.2e1 * t102 * t36 + 0.2e1 * t131 * t37, t19, 0.2e1 * t100 * t16 + 0.2e1 * t105 * t23 + 0.2e1 * t2 * t37 - 0.2e1 * t36 * t4, 0.2e1 * t1 * t37 + 0.2e1 * t101 * t16 + 0.2e1 * t106 * t23 + 0.2e1 * t3 * t36; 0, 0, 0, 0, 0, t118, -t119, 0, -pkin(7) * t118, pkin(7) * t119 (t126 * t70 + (-t126 ^ 2 * t83 + t125 * t92) * qJD(2)) * pkin(2) (t125 * t32 + t126 * t31) * pkin(2), 0, 0, -t23, -t24, 0, -t9, t8, -t23 * t56 + t95, t9, -t8, t21 * t50 - t44 * t99 - t55 * t8 + t56 * t9, t12, t11, -t102, -t103, 0, t79 * t97 + t82 * t90 + t134, t4 + t82 * t97 + (-t90 - t123) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t50, 0.2e1 * t49, 0, 0.2e1 * t50, -0.2e1 * t44, -0.2e1 * t44 * t55 + 0.2e1 * t50 * t56, t69, t65, 0, 0, 0, 0.2e1 * t129, -0.2e1 * t122 * t55 - 0.2e1 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, t24, -t23, 0, -t24, t23, t10, 0, 0, 0, 0, 0, -t103, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, -t9, t8, pkin(4) * t23 + t104, t9, -t8, -pkin(4) * t9 - qJ(5) * t8 + qJD(5) * t99, t12, t11, -t102, -t103, 0, t79 * t96 + t82 * t94 + t134, t4 + t82 * t96 + (-t94 - t123) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49, 0, t50, t85 - t49, -pkin(4) * t50 - qJ(5) * t44 + qJD(5) * t55, t69, t65, 0, 0, 0, t128 + t129, -t41 + t76 + (-qJ(5) - t55) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, qJ(5) * t85, t69, t65, 0, 0, 0, 0.2e1 * t128, -0.2e1 * t117 * t79 + 0.2e1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, t9, 0, 0, 0, 0, 0, -t102, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t100, -t23, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t121, 0, -t122 * t54 + t50 * t82, -t121 * t54 - t50 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t121, 0, t79 * t120, t82 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
