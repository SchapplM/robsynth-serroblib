% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:36
% EndTime: 2019-03-09 08:52:40
% DurationCPUTime: 1.38s
% Computational Cost: add. (3393->200), mult. (7836->368), div. (0->0), fcn. (7998->10), ass. (0->122)
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t111 = sin(qJ(2));
t113 = cos(qJ(2));
t88 = t106 * t113 + t108 * t111;
t107 = cos(pkin(11));
t110 = sin(qJ(5));
t105 = sin(pkin(11));
t112 = cos(qJ(5));
t148 = t112 * t105;
t89 = t107 * t110 + t148;
t55 = t89 * t88;
t144 = qJD(5) * t112;
t147 = t112 * t107;
t99 = pkin(2) * t106 + qJ(4);
t159 = pkin(8) + t99;
t85 = t159 * t105;
t86 = t159 * t107;
t46 = (qJD(4) * t105 + qJD(5) * t86) * t110 - qJD(4) * t147 + t85 * t144;
t145 = qJD(5) * t110;
t132 = t105 * t145;
t84 = -t107 * t144 + t132;
t162 = -0.2e1 * t84;
t82 = t88 * qJD(2);
t161 = t82 * pkin(5);
t87 = t106 * t111 - t108 * t113;
t160 = t87 * pkin(5);
t158 = cos(qJ(6));
t156 = -qJ(3) - pkin(7);
t129 = qJD(2) * t156;
t79 = t113 * qJD(3) + t111 * t129;
t80 = -t111 * qJD(3) + t113 * t129;
t53 = t106 * t79 - t108 * t80;
t92 = t156 * t111;
t93 = t156 * t113;
t64 = -t106 * t93 - t108 * t92;
t157 = t64 * t53;
t134 = -pkin(2) * t113 - pkin(1);
t59 = pkin(3) * t87 - qJ(4) * t88 + t134;
t65 = t106 * t92 - t108 * t93;
t39 = -t105 * t65 + t107 * t59;
t26 = -pkin(8) * t107 * t88 + pkin(4) * t87 + t39;
t150 = t105 * t88;
t40 = t105 * t59 + t107 * t65;
t29 = -pkin(8) * t150 + t40;
t155 = t110 * t26 + t112 * t29;
t142 = t111 * qJD(2);
t102 = pkin(2) * t142;
t141 = t113 * qJD(2);
t83 = -t106 * t142 + t108 * t141;
t45 = pkin(3) * t82 - qJ(4) * t83 - qJD(4) * t88 + t102;
t54 = t106 * t80 + t108 * t79;
t23 = t105 * t45 + t107 * t54;
t153 = -t110 * t85 + t112 * t86;
t152 = t105 * t144 + t107 * t145;
t151 = t105 * t83;
t149 = t107 * t83;
t146 = t105 ^ 2 + t107 ^ 2;
t109 = sin(qJ(6));
t143 = qJD(6) * t109;
t140 = -0.2e1 * pkin(1) * qJD(2);
t139 = t158 * pkin(5);
t31 = t83 * t148 - t88 * t132 + (t110 * t83 + t144 * t88) * t107;
t22 = -t105 * t54 + t107 * t45;
t17 = pkin(4) * t82 - pkin(8) * t149 + t22;
t20 = -pkin(8) * t151 + t23;
t6 = -t110 * t17 - t112 * t20 - t144 * t26 + t145 * t29;
t5 = -pkin(9) * t31 - t6;
t138 = t158 * t5;
t137 = pkin(5) * t143;
t136 = pkin(5) * t152;
t13 = -pkin(9) * t55 + t155;
t135 = t158 * t13;
t101 = -t108 * pkin(2) - pkin(3);
t122 = t105 * t110 - t147;
t30 = -qJD(5) * t55 - t122 * t83;
t7 = -qJD(5) * t155 - t110 * t20 + t112 * t17;
t4 = -t30 * pkin(9) + t161 + t7;
t133 = -t109 * t5 + t158 * t4;
t131 = -t110 * t29 + t112 * t26;
t130 = -t110 * t86 - t112 * t85;
t128 = qJD(6) * t139;
t41 = pkin(4) * t151 + t53;
t51 = pkin(4) * t150 + t64;
t127 = 0.2e1 * t146 * qJD(4);
t114 = t158 * t122;
t35 = qJD(6) * t114 + t109 * t152 + t143 * t89 + t158 * t84;
t61 = -t109 * t122 + t158 * t89;
t126 = t35 * t87 - t61 * t82;
t125 = t53 * t88 + t64 * t83;
t124 = t82 * t89 - t84 * t87;
t123 = -t105 * t22 + t107 * t23;
t91 = -t107 * pkin(4) + t101;
t56 = t122 * t88;
t12 = pkin(9) * t56 + t131 + t160;
t120 = t109 * t13 - t12 * t158;
t119 = t109 * t12 + t135;
t48 = -pkin(9) * t89 + t130;
t49 = -pkin(9) * t122 + t153;
t118 = t109 * t49 - t158 * t48;
t117 = t109 * t48 + t158 * t49;
t116 = t109 * t56 - t158 * t55;
t33 = -t109 * t55 - t158 * t56;
t115 = -qJD(4) * t87 + t101 * t83 - t82 * t99;
t47 = -qJD(4) * t89 - qJD(5) * t153;
t66 = pkin(5) * t122 + t91;
t62 = 0.2e1 * t87 * t82;
t60 = t109 * t89 + t114;
t44 = -t122 * t82 - t152 * t87;
t38 = t84 * pkin(9) + t47;
t37 = -pkin(9) * t152 - t46;
t36 = qJD(6) * t61 - t109 * t84 + t152 * t158;
t34 = pkin(5) * t55 + t51;
t18 = -t36 * t87 - t60 * t82;
t14 = pkin(5) * t31 + t41;
t11 = -qJD(6) * t117 - t109 * t37 + t158 * t38;
t10 = qJD(6) * t118 - t109 * t38 - t158 * t37;
t9 = qJD(6) * t33 + t109 * t30 + t158 * t31;
t8 = qJD(6) * t116 - t109 * t31 + t158 * t30;
t2 = -qJD(6) * t119 + t133;
t1 = qJD(6) * t120 - t109 * t4 - t138;
t3 = [0, 0, 0, 0.2e1 * t111 * t141, 0.2e1 * (-t111 ^ 2 + t113 ^ 2) * qJD(2), 0, 0, 0, t111 * t140, t113 * t140, -0.2e1 * t54 * t87 - 0.2e1 * t65 * t82 + 0.2e1 * t125, 0.2e1 * t102 * t134 + 0.2e1 * t65 * t54 + 0.2e1 * t157, 0.2e1 * t105 * t125 + 0.2e1 * t22 * t87 + 0.2e1 * t39 * t82, 0.2e1 * t107 * t125 - 0.2e1 * t23 * t87 - 0.2e1 * t40 * t82, 0.2e1 * (-t22 * t88 - t39 * t83) * t107 + 0.2e1 * (-t23 * t88 - t40 * t83) * t105, 0.2e1 * t22 * t39 + 0.2e1 * t23 * t40 + 0.2e1 * t157, -0.2e1 * t56 * t30, -0.2e1 * t30 * t55 + 0.2e1 * t31 * t56, 0.2e1 * t30 * t87 - 0.2e1 * t56 * t82, -0.2e1 * t31 * t87 - 0.2e1 * t55 * t82, t62, 0.2e1 * t131 * t82 + 0.2e1 * t31 * t51 + 0.2e1 * t41 * t55 + 0.2e1 * t7 * t87, -0.2e1 * t155 * t82 + 0.2e1 * t30 * t51 - 0.2e1 * t41 * t56 + 0.2e1 * t6 * t87, 0.2e1 * t33 * t8, 0.2e1 * t116 * t8 - 0.2e1 * t33 * t9, 0.2e1 * t33 * t82 + 0.2e1 * t8 * t87, 0.2e1 * t116 * t82 - 0.2e1 * t87 * t9, t62, -0.2e1 * t116 * t14 - 0.2e1 * t120 * t82 + 0.2e1 * t2 * t87 + 0.2e1 * t34 * t9, 0.2e1 * t1 * t87 - 0.2e1 * t119 * t82 + 0.2e1 * t14 * t33 + 0.2e1 * t34 * t8; 0, 0, 0, 0, 0, t141, -t142, 0, -pkin(7) * t141, pkin(7) * t142 (-t106 * t82 - t108 * t83) * pkin(2) (t106 * t54 - t108 * t53) * pkin(2), t105 * t115 - t53 * t107, t53 * t105 + t107 * t115, t123, t53 * t101 + t123 * t99 + (-t105 * t39 + t107 * t40) * qJD(4), t30 * t89 + t56 * t84, -t122 * t30 + t152 * t56 - t31 * t89 + t55 * t84, t124, t44, 0, t122 * t41 + t130 * t82 + t152 * t51 + t91 * t31 + t47 * t87, -t153 * t82 + t30 * t91 + t41 * t89 + t46 * t87 - t51 * t84, -t33 * t35 + t61 * t8, -t116 * t35 - t33 * t36 - t60 * t8 - t61 * t9, -t126, t18, 0, t11 * t87 - t116 * t136 - t118 * t82 + t14 * t60 + t34 * t36 + t66 * t9, t10 * t87 - t117 * t82 + t136 * t33 + t14 * t61 - t34 * t35 + t66 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t99 * t127, t89 * t162, 0.2e1 * t122 * t84 - 0.2e1 * t152 * t89, 0, 0, 0, 0.2e1 * t91 * t152, t91 * t162, -0.2e1 * t61 * t35, 0.2e1 * t35 * t60 - 0.2e1 * t36 * t61, 0, 0, 0, 0.2e1 * t136 * t60 + 0.2e1 * t36 * t66, 0.2e1 * t136 * t61 - 0.2e1 * t35 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t107 * t82, -t105 * t82, -t146 * t83, t105 * t23 + t107 * t22, 0, 0, 0, 0, 0, t44, -t124, 0, 0, 0, 0, 0, t18, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t149, 0, t53, 0, 0, 0, 0, 0, t31, t30, 0, 0, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t84, 0, 0, 0, 0, 0, t36, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31, t82, t7, t6, 0, 0, t8, -t9, t82, t82 * t139 + (-t135 + (-t12 - t160) * t109) * qJD(6) + t133, -t138 + (-t4 - t161) * t109 + (-t139 * t87 + t120) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t152, 0, t47, t46, 0, 0, -t35, -t36, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t84, 0, 0, 0, 0, 0, -t36, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t137, -0.2e1 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t82, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, -t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
