% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:01
% EndTime: 2019-03-08 23:14:05
% DurationCPUTime: 1.29s
% Computational Cost: add. (1311->182), mult. (3359->313), div. (0->0), fcn. (3182->10), ass. (0->118)
t73 = sin(qJ(2));
t123 = qJD(2) * t73;
t141 = qJD(3) + qJD(4);
t71 = sin(qJ(4));
t75 = cos(qJ(3));
t130 = t71 * t75;
t137 = cos(qJ(4));
t72 = sin(qJ(3));
t45 = t137 * t72 + t130;
t33 = t141 * t45;
t108 = t137 * t75;
t131 = t71 * t72;
t44 = -t108 + t131;
t69 = sin(pkin(6));
t76 = cos(qJ(2));
t143 = (t44 * t123 - t33 * t76) * t69;
t102 = qJD(3) * t108;
t106 = t137 * qJD(4);
t32 = -t75 * t106 + t141 * t131 - t102;
t142 = (t45 * t123 + t32 * t76) * t69;
t70 = sin(qJ(6));
t67 = t70 ^ 2;
t74 = cos(qJ(6));
t125 = -t74 ^ 2 + t67;
t105 = t125 * qJD(6);
t122 = qJD(2) * t76;
t110 = t69 * t122;
t124 = cos(pkin(6));
t134 = t69 * t73;
t86 = -t124 * t75 + t72 * t134;
t140 = t86 * qJD(3) - t75 * t110;
t78 = 2 * qJD(5);
t139 = pkin(4) + pkin(10);
t138 = -pkin(9) - pkin(8);
t136 = t44 * t70;
t135 = t44 * t74;
t133 = t69 * t76;
t132 = t70 * t33;
t129 = t74 * t33;
t118 = qJD(6) * t74;
t107 = qJD(3) * t138;
t121 = qJD(4) * t71;
t47 = t72 * t107;
t51 = t138 * t72;
t52 = t138 * t75;
t17 = -t51 * t106 - t107 * t130 - t52 * t121 - t137 * t47;
t12 = -t33 * pkin(5) - t17;
t90 = -t137 * t52 + t71 * t51;
t25 = -t44 * pkin(5) + t90;
t128 = t25 * t118 + t12 * t70;
t103 = pkin(3) * t106;
t54 = t103 + qJD(5);
t58 = t71 * pkin(3) + qJ(5);
t127 = t58 * t118 + t54 * t70;
t114 = qJ(5) * qJD(6);
t126 = qJD(5) * t70 + t74 * t114;
t120 = qJD(6) * t25;
t119 = qJD(6) * t70;
t117 = qJD(6) * t139;
t116 = t72 * qJD(3);
t115 = t75 * qJD(3);
t113 = -0.2e1 * pkin(2) * qJD(3);
t112 = t70 * t129;
t64 = pkin(3) * t116;
t63 = pkin(3) * t121;
t111 = t69 * t123;
t109 = t70 * t118;
t62 = -t75 * pkin(3) - pkin(2);
t61 = -t137 * pkin(3) - pkin(4);
t34 = -t137 * t51 - t71 * t52;
t98 = -t45 * qJ(5) + t62;
t21 = t139 * t44 + t98;
t24 = t45 * pkin(5) + t34;
t101 = t74 * t21 + t70 * t24;
t100 = t70 * t21 - t74 * t24;
t99 = -qJ(5) * t33 - qJD(5) * t44;
t37 = t124 * t72 + t75 * t134;
t82 = t137 * t86;
t26 = t71 * t37 + t82;
t97 = t74 * t133 - t70 * t26;
t96 = t70 * t133 + t74 * t26;
t93 = t45 * t118 - t70 * t32;
t19 = -t45 * t119 - t74 * t32;
t92 = t44 * t118 + t132;
t91 = t44 * t119 - t129;
t89 = t32 * qJ(5) - t45 * qJD(5) + t64;
t57 = -pkin(10) + t61;
t88 = qJD(6) * (t44 * t58 - t45 * t57);
t87 = qJD(6) * (qJ(5) * t44 + t139 * t45);
t85 = t139 * t32 + t99;
t84 = t71 * t86;
t83 = -t58 * t33 - t54 * t44 + t45 * t63;
t81 = -t32 * t57 + t83;
t18 = t90 * qJD(4) - t138 * t102 + t71 * t47;
t80 = t37 * qJD(3) + t72 * t110;
t66 = qJD(5) * t74;
t53 = -0.2e1 * t109;
t49 = t54 * t74;
t43 = 0.2e1 * t105;
t42 = t44 ^ 2;
t28 = t44 * pkin(4) + t98;
t27 = t137 * t37 - t84;
t23 = -0.2e1 * t45 * t32;
t16 = -t44 * t105 + t112;
t15 = t33 * pkin(4) + t89;
t14 = -0.4e1 * t44 * t109 - t125 * t33;
t13 = -t32 * pkin(5) + t18;
t11 = t12 * t74;
t9 = -qJD(4) * t84 + t37 * t106 + t137 * t80 - t140 * t71;
t8 = qJD(4) * t82 + t37 * t121 + t137 * t140 + t71 * t80;
t7 = t139 * t33 + t89;
t6 = -t27 * t119 - t8 * t74;
t5 = t27 * t118 - t8 * t70;
t4 = t96 * qJD(6) + t74 * t111 + t70 * t9;
t3 = t97 * qJD(6) - t70 * t111 + t74 * t9;
t2 = -t101 * qJD(6) + t74 * t13 - t70 * t7;
t1 = t100 * qJD(6) - t70 * t13 - t74 * t7;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t69 ^ 2 * t73 * t122 + 0.2e1 * t26 * t9 - 0.2e1 * t27 * t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t111, -t110, 0, 0, 0, 0, 0 (-t76 * t116 - t75 * t123) * t69 (-t76 * t115 + t72 * t123) * t69, 0, 0, 0, 0, 0, t143, t142, -t26 * t32 - t27 * t33 + t8 * t44 + t9 * t45, -t143, -t142, -t27 * t17 + t26 * t18 + t9 * t34 - t8 * t90 + (t28 * t123 - t15 * t76) * t69, 0, 0, 0, 0, 0, t8 * t135 + t91 * t27 + t3 * t45 - t96 * t32, -t8 * t136 + t27 * t92 - t32 * t97 - t4 * t45; 0, 0, 0, 0, 0.2e1 * t72 * t115, 0.2e1 * (-t72 ^ 2 + t75 ^ 2) * qJD(3), 0, 0, 0, t72 * t113, t75 * t113, t23, 0.2e1 * t32 * t44 - 0.2e1 * t45 * t33, 0, 0, 0, 0.2e1 * t62 * t33 + 0.2e1 * t44 * t64, -0.2e1 * t62 * t32 + 0.2e1 * t45 * t64, 0.2e1 * t17 * t44 + 0.2e1 * t18 * t45 - 0.2e1 * t34 * t32 - 0.2e1 * t33 * t90, -0.2e1 * t15 * t44 - 0.2e1 * t28 * t33, -0.2e1 * t15 * t45 + 0.2e1 * t28 * t32, 0.2e1 * t28 * t15 - 0.2e1 * t17 * t90 + 0.2e1 * t34 * t18, 0.2e1 * t67 * t44 * t33 + 0.2e1 * t42 * t109, -0.2e1 * t42 * t105 + 0.4e1 * t44 * t112, 0.2e1 * t45 * t132 + 0.2e1 * t93 * t44, 0.2e1 * t45 * t129 + 0.2e1 * t19 * t44, t23, 0.2e1 * t100 * t32 - 0.2e1 * t12 * t135 + 0.2e1 * t2 * t45 + 0.2e1 * t91 * t25, 0.2e1 * t1 * t45 + 0.2e1 * t101 * t32 + 0.2e1 * t12 * t136 + 0.2e1 * t25 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t140, 0, 0, 0, 0, 0, -t9, t8, 0, t9, -t8, t26 * t63 + t27 * t54 - t8 * t58 + t9 * t61, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, t115, -t116, 0, -pkin(8) * t115, pkin(8) * t116, 0, 0, -t32, -t33, 0, -t18, t17, -t61 * t32 + t83, t18, -t17, -t17 * t58 + t18 * t61 + t34 * t63 + t54 * t90, t16, t14, t19, -t93, 0, t70 * t88 + t74 * t81 + t128, t11 + t74 * t88 + (-t81 - t120) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t63, -0.2e1 * t103, 0, 0.2e1 * t63, 0.2e1 * t54, 0.2e1 * t58 * t54 + 0.2e1 * t61 * t63, t53, t43, 0, 0, 0, 0.2e1 * t127, -0.2e1 * t58 * t119 + 0.2e1 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, 0, t9, -t8, -t9 * pkin(4) - t8 * qJ(5) + t27 * qJD(5), 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, -t18, t17, pkin(4) * t32 + t99, t18, -t17, -t18 * pkin(4) - t17 * qJ(5) + qJD(5) * t90, t16, t14, t19, -t93, 0, t70 * t87 + t74 * t85 + t128, t11 + t74 * t87 + (-t85 - t120) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t103, 0, t63, t78 + t103, -pkin(4) * t63 + t54 * qJ(5) + t58 * qJD(5), t53, t43, 0, 0, 0, t126 + t127, t49 + t66 + (-qJ(5) - t58) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, qJ(5) * t78, t53, t43, 0, 0, 0, 0.2e1 * t126, -0.2e1 * t114 * t70 + 0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, t18, 0, 0, 0, 0, 0, t19, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t91, -t32, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t118, 0, -t57 * t119 + t63 * t74, -t118 * t57 - t63 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t118, 0, t70 * t117, t74 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
