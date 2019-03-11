% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:57
% EndTime: 2019-03-09 10:10:02
% DurationCPUTime: 1.66s
% Computational Cost: add. (4122->181), mult. (9194->338), div. (0->0), fcn. (9361->10), ass. (0->126)
t108 = sin(pkin(11));
t110 = cos(pkin(11));
t144 = t108 ^ 2 + t110 ^ 2;
t166 = t144 * qJD(5);
t167 = 0.2e1 * t166;
t111 = cos(pkin(10));
t100 = t111 * pkin(2) + pkin(3);
t113 = sin(qJ(4));
t109 = sin(pkin(10));
t163 = t109 * pkin(2);
t164 = cos(qJ(4));
t122 = -t164 * t100 + t113 * t163;
t75 = t122 * qJD(4);
t71 = qJD(5) - t75;
t165 = t144 * t71;
t112 = sin(qJ(6));
t115 = cos(qJ(6));
t146 = t115 * t108;
t88 = t112 * t110 + t146;
t85 = t88 * qJD(6);
t162 = t110 * pkin(5);
t105 = t110 * pkin(9);
t161 = -qJ(3) - pkin(7);
t114 = sin(qJ(2));
t116 = cos(qJ(2));
t87 = t109 * t116 + t111 * t114;
t123 = t87 * qJD(2);
t124 = t109 * t114 - t111 * t116;
t60 = t113 * t87 + t164 * t124;
t140 = t116 * qJD(2);
t141 = t114 * qJD(2);
t83 = -t109 * t141 + t111 * t140;
t41 = -t60 * qJD(4) - t113 * t123 + t164 * t83;
t153 = t108 * t41;
t134 = qJD(2) * t161;
t80 = t116 * qJD(3) + t114 * t134;
t81 = -t114 * qJD(3) + t116 * t134;
t52 = t109 * t81 + t111 * t80;
t117 = -pkin(8) * t123 + t52;
t51 = -t109 * t80 + t111 * t81;
t120 = -t83 * pkin(8) + t51;
t135 = qJD(4) * t164;
t143 = qJD(4) * t113;
t94 = t161 * t114;
t95 = t161 * t116;
t65 = t109 * t95 + t111 * t94;
t53 = -t87 * pkin(8) + t65;
t66 = t109 * t94 - t111 * t95;
t54 = -t124 * pkin(8) + t66;
t20 = t113 * t117 - t164 * t120 + t54 * t135 + t53 * t143;
t12 = pkin(5) * t153 + t20;
t61 = -t113 * t124 + t164 * t87;
t152 = t108 * t61;
t37 = t113 * t54 - t164 * t53;
t30 = pkin(5) * t152 + t37;
t147 = t112 * t108;
t86 = -t115 * t110 + t147;
t160 = t12 * t86 + t30 * t85;
t136 = qJD(6) * t147;
t142 = qJD(6) * t115;
t84 = -t110 * t142 + t136;
t159 = t12 * t88 - t30 * t84;
t19 = -t113 * t120 - t164 * t117 - t53 * t135 + t54 * t143;
t42 = t61 * qJD(4) + t113 * t83 + t164 * t123;
t102 = pkin(2) * t141;
t69 = pkin(3) * t123 + t102;
t23 = t42 * pkin(4) - t41 * qJ(5) - t61 * qJD(5) + t69;
t8 = t108 * t23 - t110 * t19;
t138 = -t116 * pkin(2) - pkin(1);
t72 = t124 * pkin(3) + t138;
t36 = t60 * pkin(4) - t61 * qJ(5) + t72;
t38 = t113 * t53 + t164 * t54;
t25 = t108 * t36 + t110 * t38;
t79 = -pkin(4) + t122;
t70 = t79 - t162;
t121 = t113 * t100 + t164 * t163;
t76 = t121 * qJD(4);
t158 = t70 * t85 + t76 * t86;
t157 = -t70 * t84 + t76 * t88;
t101 = -pkin(4) - t162;
t155 = t101 * t84;
t154 = t101 * t85;
t151 = t110 * t41;
t150 = t20 * t110;
t149 = t76 * t108;
t148 = t76 * t110;
t139 = -0.2e1 * pkin(1) * qJD(2);
t7 = t108 * t19 + t110 * t23;
t24 = -t108 * t38 + t110 * t36;
t3 = -t7 * t108 + t8 * t110;
t131 = t20 * t61 + t37 * t41;
t28 = t88 * t42 - t84 * t60;
t13 = t60 * pkin(5) - t61 * t105 + t24;
t16 = -pkin(9) * t152 + t25;
t130 = t112 * t16 - t115 * t13;
t129 = t112 * t13 + t115 * t16;
t78 = qJ(5) + t121;
t63 = (-pkin(9) - t78) * t108;
t64 = t110 * t78 + t105;
t128 = t112 * t64 - t115 * t63;
t127 = t112 * t63 + t115 * t64;
t92 = (-pkin(9) - qJ(5)) * t108;
t93 = t110 * qJ(5) + t105;
t126 = t112 * t93 - t115 * t92;
t125 = t112 * t92 + t115 * t93;
t119 = -pkin(4) * t41 - qJ(5) * t42 - qJD(5) * t60;
t118 = t41 * t79 - t42 * t78 - t60 * t71 + t61 * t76;
t62 = -0.2e1 * t88 * t84;
t49 = -t88 * qJD(5) - t125 * qJD(6);
t48 = t86 * qJD(5) + t126 * qJD(6);
t43 = 0.2e1 * t84 * t86 - 0.2e1 * t88 * t85;
t40 = t86 * t61;
t39 = t88 * t61;
t32 = -t127 * qJD(6) - t88 * t71;
t31 = t128 * qJD(6) + t86 * t71;
t29 = -t86 * t42 - t85 * t60;
t18 = t20 * t108;
t15 = t41 * t146 - t61 * t136 + (t112 * t41 + t61 * t142) * t110;
t14 = -t86 * t41 - t61 * t85;
t9 = t14 * t88 + t40 * t84;
t6 = -pkin(9) * t153 + t8;
t5 = t42 * pkin(5) - pkin(9) * t151 + t7;
t4 = -t14 * t86 - t88 * t15 + t84 * t39 + t40 * t85;
t2 = -t129 * qJD(6) - t112 * t6 + t115 * t5;
t1 = t130 * qJD(6) - t112 * t5 - t115 * t6;
t10 = [0, 0, 0, 0.2e1 * t114 * t140, 0.2e1 * (-t114 ^ 2 + t116 ^ 2) * qJD(2), 0, 0, 0, t114 * t139, t116 * t139, -0.2e1 * t66 * t123 - 0.2e1 * t52 * t124 - 0.2e1 * t51 * t87 - 0.2e1 * t65 * t83, 0.2e1 * t138 * t102 + 0.2e1 * t65 * t51 + 0.2e1 * t66 * t52, 0.2e1 * t61 * t41, -0.2e1 * t41 * t60 - 0.2e1 * t61 * t42, 0, 0, 0, 0.2e1 * t72 * t42 + 0.2e1 * t69 * t60, 0.2e1 * t72 * t41 + 0.2e1 * t69 * t61, 0.2e1 * t131 * t108 + 0.2e1 * t24 * t42 + 0.2e1 * t7 * t60, 0.2e1 * t131 * t110 - 0.2e1 * t25 * t42 - 0.2e1 * t8 * t60, 0.2e1 * (-t24 * t41 - t61 * t7) * t110 + 0.2e1 * (-t25 * t41 - t61 * t8) * t108, 0.2e1 * t37 * t20 + 0.2e1 * t24 * t7 + 0.2e1 * t25 * t8, -0.2e1 * t40 * t14, -0.2e1 * t14 * t39 + 0.2e1 * t40 * t15, 0.2e1 * t14 * t60 - 0.2e1 * t40 * t42, -0.2e1 * t15 * t60 - 0.2e1 * t39 * t42, 0.2e1 * t60 * t42, 0.2e1 * t12 * t39 - 0.2e1 * t130 * t42 + 0.2e1 * t30 * t15 + 0.2e1 * t2 * t60, 0.2e1 * t1 * t60 - 0.2e1 * t12 * t40 - 0.2e1 * t129 * t42 + 0.2e1 * t30 * t14; 0, 0, 0, 0, 0, t140, -t141, 0, -pkin(7) * t140, pkin(7) * t141 (-t109 * t123 - t111 * t83) * pkin(2) (t52 * t109 + t111 * t51) * pkin(2), 0, 0, t41, -t42, 0, -t20, t19, t108 * t118 - t150, t110 * t118 + t18, t3, t20 * t79 + t37 * t76 + (t25 * t71 + t78 * t8) * t110 + (-t24 * t71 - t7 * t78) * t108, t9, t4, t28, t29, 0, -t128 * t42 + t70 * t15 + t32 * t60 + t76 * t39 + t160, -t127 * t42 + t70 * t14 + t31 * t60 - t76 * t40 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t76, 0.2e1 * t75, -0.2e1 * t148, 0.2e1 * t149, 0.2e1 * t165, 0.2e1 * t165 * t78 + 0.2e1 * t79 * t76, t62, t43, 0, 0, 0, 0.2e1 * t158, 0.2e1 * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, t42, t41, t110 * t42, -t108 * t42, -t144 * t41, t8 * t108 + t7 * t110, 0, 0, 0, 0, 0, t29, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t42, 0, -t20, t19, t119 * t108 - t150, t119 * t110 + t18, t3, -t20 * pkin(4) + (-t108 * t24 + t110 * t25) * qJD(5) + t3 * qJ(5), t9, t4, t28, t29, 0, t101 * t15 - t126 * t42 + t49 * t60 + t160, t101 * t14 - t125 * t42 + t48 * t60 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75, -t148, t149, t166 + t165, -t76 * pkin(4) + qJ(5) * t165 + t166 * t78, t62, t43, 0, 0, 0, t154 + t158, -t155 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, qJ(5) * t167, t62, t43, 0, 0, 0, 0.2e1 * t154, -0.2e1 * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t151, 0, t20, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, t85, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t42, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t85, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t85, 0, t49, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
