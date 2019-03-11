% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:06
% EndTime: 2019-03-09 05:57:11
% DurationCPUTime: 1.59s
% Computational Cost: add. (2585->198), mult. (5445->323), div. (0->0), fcn. (4916->8), ass. (0->125)
t82 = cos(qJ(5));
t79 = sin(qJ(5));
t96 = t82 * pkin(5) + t79 * qJ(6);
t151 = t96 * qJD(5) - t82 * qJD(6);
t145 = cos(qJ(4));
t102 = qJD(4) * t145;
t103 = qJD(3) * t145;
t80 = sin(qJ(4));
t81 = sin(qJ(3));
t137 = t80 * t81;
t152 = qJD(3) + qJD(4);
t83 = cos(qJ(3));
t39 = t152 * t137 + (-t103 - t102) * t83;
t60 = t145 * t81 + t80 * t83;
t95 = t79 * pkin(5) - t82 * qJ(6);
t156 = -t151 * t60 + t95 * t39;
t59 = -t145 * t83 + t137;
t70 = -cos(pkin(10)) * pkin(1) - pkin(2);
t63 = -pkin(3) * t83 + t70;
t33 = t59 * pkin(4) - t60 * pkin(9) + t63;
t69 = sin(pkin(10)) * pkin(1) + pkin(7);
t146 = pkin(8) + t69;
t55 = t146 * t81;
t56 = t146 * t83;
t35 = t145 * t56 - t80 * t55;
t154 = t79 * t33 + t82 * t35;
t76 = t79 ^ 2;
t77 = t82 ^ 2;
t153 = t76 + t77;
t101 = (t76 - t77) * qJD(5);
t126 = qJD(4) * t80;
t89 = t146 * t103;
t99 = qJD(3) * t80 * t146;
t12 = t55 * t102 + t56 * t126 + t81 * t89 + t83 * t99;
t121 = t81 * qJD(3);
t116 = pkin(3) * t121;
t40 = t152 * t60;
t19 = pkin(4) * t40 + pkin(9) * t39 + t116;
t5 = -qJD(5) * t154 + t12 * t79 + t82 * t19;
t150 = 0.2e1 * qJD(6);
t149 = pkin(5) * t40;
t148 = pkin(9) * t40;
t147 = pkin(9) * t59;
t144 = t39 * t76;
t72 = pkin(3) * t80 + pkin(9);
t143 = t40 * t72;
t142 = t59 * t40;
t141 = t59 * t72;
t140 = t60 * t39;
t139 = t60 * t79;
t138 = t60 * t82;
t36 = t77 * t39;
t136 = t82 * t39;
t135 = t82 * t40;
t13 = t35 * qJD(4) - t81 * t99 + t83 * t89;
t34 = t145 * t55 + t80 * t56;
t75 = qJD(5) * t82;
t134 = t13 * t79 + t34 * t75;
t133 = t60 * t135 - t59 * t136;
t113 = pkin(3) * t126;
t122 = qJD(6) * t79;
t124 = qJD(5) * t79;
t52 = pkin(5) * t124 - qJ(6) * t75 - t122;
t41 = t52 + t113;
t131 = -t41 - t52;
t110 = t145 * pkin(3);
t73 = -t110 - pkin(4);
t130 = t79 * t113 + t73 * t75;
t128 = pkin(3) * qJD(4);
t127 = qJ(6) * t40;
t18 = t95 * t60 + t34;
t125 = qJD(5) * t18;
t123 = qJD(6) * t59;
t119 = t83 * qJD(3);
t118 = t79 * t136;
t117 = 0.2e1 * t119;
t115 = pkin(4) * t124;
t114 = pkin(4) * t75;
t112 = pkin(9) * t124;
t111 = pkin(9) * t75;
t109 = t60 * t124;
t108 = t79 * t75;
t107 = t79 * t145;
t106 = t82 * t145;
t105 = t153 * t39;
t100 = pkin(3) * t102;
t7 = qJ(6) * t59 + t154;
t94 = t33 * t82 - t35 * t79;
t8 = -pkin(5) * t59 - t94;
t98 = -t7 * t82 - t79 * t8;
t97 = t7 * t79 - t8 * t82;
t92 = t39 * t59 - t40 * t60;
t91 = -t60 * t73 + t141;
t90 = -t82 * t113 + t73 * t124;
t64 = -pkin(4) - t96;
t24 = -t39 * t79 + t60 * t75;
t22 = t109 + t136;
t21 = t59 * t124 - t135;
t4 = t82 * t12 + t35 * t124 - t79 * t19 - t33 * t75;
t88 = -t39 * t64 + t52 * t60 - t148;
t6 = t13 - t156;
t87 = -t6 + (t60 * t64 - t147) * qJD(5);
t54 = -t110 + t64;
t86 = -t6 + (t54 * t60 - t141) * qJD(5);
t50 = t153 * t145 * t128;
t2 = t123 - t4 + t127;
t3 = -t149 - t5;
t1 = -t97 * qJD(5) + t2 * t82 + t3 * t79;
t85 = -t59 * t100 - t39 * t54 + t41 * t60 - t143;
t84 = -t39 * t73 - t143 + (-t145 * t59 + t60 * t80) * t128;
t66 = 0.2e1 * t108;
t58 = -0.2e1 * t101;
t57 = t60 ^ 2;
t53 = t64 * t124;
t47 = t54 * t124;
t46 = t79 * t100 + t72 * t75;
t45 = -t82 * t100 + t72 * t124;
t30 = t34 * t124;
t25 = t60 * t36;
t23 = t40 * t79 + t59 * t75;
t20 = -t36 - t144;
t17 = -t60 * t101 - t118;
t14 = t18 * t124;
t9 = -0.4e1 * t60 * t108 + t144 - t36;
t10 = [0, 0, 0, 0, t81 * t117, 0.2e1 * (-t81 ^ 2 + t83 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t70 * t121, t70 * t117, -0.2e1 * t140, 0.2e1 * t92, 0, 0, 0, 0.2e1 * t59 * t116 + 0.2e1 * t40 * t63, 0.2e1 * t60 * t116 - 0.2e1 * t39 * t63, -0.2e1 * t57 * t108 - 0.2e1 * t25, 0.2e1 * t57 * t101 + 0.4e1 * t60 * t118, -0.2e1 * t59 * t109 + 0.2e1 * t133, -0.2e1 * t40 * t139 - 0.2e1 * t24 * t59, 0.2e1 * t142, 0.2e1 * t13 * t139 + 0.2e1 * t24 * t34 + 0.2e1 * t40 * t94 + 0.2e1 * t5 * t59, 0.2e1 * t13 * t138 - 0.2e1 * t154 * t40 - 0.2e1 * t22 * t34 + 0.2e1 * t4 * t59, 0.2e1 * t6 * t139 + 0.2e1 * t18 * t24 - 0.2e1 * t3 * t59 - 0.2e1 * t40 * t8, 0.2e1 * t97 * t39 + 0.2e1 * (t98 * qJD(5) - t2 * t79 + t3 * t82) * t60, -0.2e1 * t6 * t138 + 0.2e1 * t18 * t22 + 0.2e1 * t2 * t59 + 0.2e1 * t40 * t7, 0.2e1 * t18 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t3 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t92 + t133, 0, 0, 0, t1 * t60 + t18 * t40 + t98 * t39 + t59 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t76 * t140 + 0.2e1 * t142 - 0.2e1 * t25; 0, 0, 0, 0, 0, 0, t119, -t121, 0, -t69 * t119, t69 * t121, 0, 0, -t39, -t40, 0, -t13, t12, t17, t9, t23, -t21, 0, t30 + (-qJD(5) * t91 - t13) * t82 + t84 * t79, t91 * t124 + t82 * t84 + t134, t79 * t85 + t82 * t86 + t14, t1, t86 * t79 + (-t85 - t125) * t82, t18 * t41 + t6 * t54 + (t7 * t106 + t8 * t107) * t128 + t1 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t119, 0, 0, 0, 0, 0, -t40, t39, 0, 0, 0, 0, 0, t21, t23, t21, t20, -t23, -t72 * t105 + t40 * t54 + t59 * t41 + t60 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t113, -0.2e1 * t100, t66, t58, 0, 0, 0, 0.2e1 * t90, 0.2e1 * t130, -0.2e1 * t41 * t82 + 0.2e1 * t47, 0.2e1 * t50, -0.2e1 * t41 * t79 - 0.2e1 * t54 * t75, 0.2e1 * t54 * t41 + 0.2e1 * t50 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t40, 0, -t13, t12, t17, t9, t23, -t21, 0, t30 + (pkin(4) * t39 - t148) * t79 + (-t13 + (-pkin(4) * t60 - t147) * qJD(5)) * t82, pkin(4) * t22 + pkin(9) * t21 + t134, t79 * t88 + t82 * t87 + t14, t1, t87 * t79 + (-t88 - t125) * t82, pkin(9) * t1 + t18 * t52 + t6 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t39, 0, 0, 0, 0, 0, t21, t23, t21, t20, -t23, -pkin(9) * t105 + t40 * t64 + t52 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t100, t66, t58, 0, 0, 0, t90 - t115, -t114 + t130, t131 * t82 + t47 + t53, t50, t131 * t79 + (-t54 - t64) * t75, pkin(9) * t50 + t41 * t64 + t54 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t58, 0, 0, 0, -0.2e1 * t115, -0.2e1 * t114, -0.2e1 * t52 * t82 + 0.2e1 * t53, 0, -0.2e1 * t52 * t79 - 0.2e1 * t64 * t75, 0.2e1 * t64 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t24, t40, t5, t4, t5 + 0.2e1 * t149, t96 * t39 + (qJD(5) * t95 - t122) * t60, 0.2e1 * t123 - t4 + 0.2e1 * t127, -pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t22, -t24, 0, -t22, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t124, 0, -t46, t45, -t46, -t151, -t45 (-pkin(5) * t107 + qJ(6) * t106) * t128 - t151 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t124, 0, -t111, t112, -t111, -t151, -t112, -t151 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, qJ(6) * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t22, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
