% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:17
% EndTime: 2019-03-08 23:59:22
% DurationCPUTime: 1.62s
% Computational Cost: add. (2161->230), mult. (5326->393), div. (0->0), fcn. (5104->10), ass. (0->133)
t84 = sin(pkin(6));
t91 = cos(qJ(2));
t150 = t84 * t91;
t123 = qJD(2) * t150;
t87 = sin(qJ(3));
t154 = sin(qJ(2));
t121 = t84 * t154;
t141 = cos(pkin(6));
t90 = cos(qJ(3));
t98 = t90 * t121 + t141 * t87;
t94 = t98 * qJD(3) + t87 * t123;
t160 = qJD(4) * t98 + t94;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t59 = t86 * t87 - t89 * t90;
t60 = t86 * t90 + t87 * t89;
t77 = -pkin(3) * t90 - pkin(2);
t35 = pkin(4) * t59 - pkin(10) * t60 + t77;
t156 = -pkin(9) - pkin(8);
t69 = t156 * t87;
t70 = t156 * t90;
t46 = t69 * t86 - t70 * t89;
t88 = cos(qJ(5));
t44 = t88 * t46;
t85 = sin(qJ(5));
t146 = t85 * t35 + t44;
t99 = t87 * t121 - t141 * t90;
t157 = t99 * qJD(3) - t90 * t123;
t159 = -qJD(4) * t99 - t157;
t83 = t88 ^ 2;
t144 = t85 ^ 2 - t83;
t114 = t144 * qJD(5);
t158 = qJD(3) + qJD(4);
t155 = pkin(3) * t89;
t42 = t158 * t59;
t153 = t60 * t42;
t152 = t60 * t85;
t151 = t60 * t88;
t149 = t85 * t88;
t148 = -qJ(6) - pkin(10);
t120 = qJD(3) * t156;
t109 = t90 * t120;
t110 = t87 * t120;
t24 = t46 * qJD(4) - t89 * t109 + t86 * t110;
t45 = -t69 * t89 - t86 * t70;
t80 = qJD(5) * t88;
t147 = t24 * t85 + t45 * t80;
t137 = qJD(4) * t86;
t128 = pkin(3) * t137;
t75 = -pkin(4) - t155;
t145 = t85 * t128 + t75 * t80;
t143 = qJ(6) * t60;
t81 = t88 * qJ(6);
t74 = pkin(3) * t86 + pkin(10);
t142 = -qJ(6) - t74;
t140 = qJD(3) * t87;
t139 = qJD(3) * t90;
t138 = qJD(3) * t91;
t136 = qJD(4) * t89;
t135 = qJD(5) * t85;
t134 = -0.2e1 * pkin(2) * qJD(3);
t133 = t85 * t150;
t131 = pkin(3) * t140;
t43 = t158 * t60;
t20 = pkin(4) * t43 + pkin(10) * t42 + t131;
t23 = -t86 * t109 - t89 * t110 - t69 * t136 - t70 * t137;
t132 = t85 * t20 - t88 * t23 + t35 * t80;
t130 = pkin(4) * t135;
t129 = pkin(4) * t80;
t127 = pkin(3) * t136;
t78 = pkin(5) * t135;
t126 = pkin(5) * t80;
t30 = t86 * t98 + t89 * t99;
t125 = t30 * t135;
t124 = t60 * t80;
t122 = t85 * t80;
t76 = -pkin(5) * t88 - pkin(4);
t119 = -0.4e1 * t60 * t149;
t118 = t88 * t20 + t23 * t85;
t117 = t88 * t35 - t46 * t85;
t116 = qJD(2) * t154;
t115 = qJD(5) * t148;
t113 = qJD(5) * t142;
t112 = t88 * t127;
t108 = t84 * t116;
t31 = -t86 * t99 + t89 * t98;
t103 = t88 * t150 + t85 * t31;
t26 = t88 * t31 - t133;
t107 = t103 * t88 - t26 * t85;
t106 = t59 * t74 - t60 * t75;
t105 = qJ(6) * t42 - qJD(6) * t60;
t104 = -t88 * t128 + t75 * t135;
t102 = -t42 * t85 + t124;
t101 = t60 * t135 + t42 * t88;
t100 = t59 * t135 - t43 * t88;
t95 = -t42 * t75 - t43 * t74 + (-t59 * t89 + t60 * t86) * qJD(4) * pkin(3);
t92 = -t159 * t89 + t160 * t86;
t79 = t88 * qJD(6);
t72 = 0.2e1 * t122;
t68 = pkin(10) * t88 + t81;
t67 = t148 * t85;
t66 = t76 - t155;
t62 = t78 + t128;
t58 = -0.2e1 * t114;
t57 = t60 ^ 2;
t55 = t74 * t88 + t81;
t54 = t142 * t85;
t50 = -qJD(6) * t85 + t88 * t115;
t49 = t85 * t115 + t79;
t48 = t49 * t88;
t41 = (-qJD(6) - t127) * t85 + t88 * t113;
t40 = t85 * t113 + t112 + t79;
t38 = t45 * t135;
t36 = t40 * t88;
t29 = pkin(5) * t152 + t45;
t28 = t43 * t85 + t59 * t80;
t19 = -t60 * t114 - t42 * t149;
t16 = -t85 * t143 + t146;
t15 = qJD(5) * t119 + t144 * t42;
t14 = pkin(5) * t59 - t60 * t81 + t117;
t13 = t159 * t86 + t160 * t89;
t11 = t102 * pkin(5) + t24;
t10 = -t13 * t88 + t125;
t9 = t13 * t85 + t30 * t80;
t8 = qJD(5) * t133 + t88 * t108 - t31 * t80 + t85 * t92;
t7 = t103 * qJD(5) - t85 * t108 + t88 * t92;
t6 = -t146 * qJD(5) + t118;
t5 = t46 * t135 - t132;
t4 = -qJ(6) * t124 + (-qJD(5) * t46 + t105) * t85 + t132;
t3 = t4 * t88;
t2 = pkin(5) * t43 + t105 * t88 + (-t44 + (-t35 + t143) * t85) * qJD(5) + t118;
t1 = t107 * qJD(5) - t7 * t88 - t8 * t85;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t103 * t8 + 0.2e1 * t13 * t30 - 0.2e1 * t26 * t7; 0, 0, -t108, -t123, 0, 0, 0, 0, 0 (-t90 * t116 - t87 * t138) * t84 (t87 * t116 - t90 * t138) * t84, 0, 0, 0, 0, 0 (t59 * t116 - t43 * t91) * t84 (t60 * t116 + t42 * t91) * t84, 0, 0, 0, 0, 0, t102 * t30 - t103 * t43 + t13 * t152 + t59 * t8, -t101 * t30 + t13 * t151 - t26 * t43 + t59 * t7, -t107 * t42 + (t7 * t85 - t8 * t88 + (-t103 * t85 - t26 * t88) * qJD(5)) * t60, -t103 * t2 + t11 * t30 + t13 * t29 + t14 * t8 - t16 * t7 + t26 * t4; 0, 0, 0, 0, 0.2e1 * t87 * t139, 0.2e1 * (-t87 ^ 2 + t90 ^ 2) * qJD(3), 0, 0, 0, t87 * t134, t90 * t134, -0.2e1 * t153, 0.2e1 * t42 * t59 - 0.2e1 * t43 * t60, 0, 0, 0, 0.2e1 * t59 * t131 + 0.2e1 * t43 * t77, 0.2e1 * t60 * t131 - 0.2e1 * t42 * t77, -0.2e1 * t57 * t122 - 0.2e1 * t83 * t153, 0.2e1 * t57 * t114 - t42 * t119, -0.2e1 * t101 * t59 + 0.2e1 * t43 * t151, -0.2e1 * t102 * t59 - 0.2e1 * t43 * t152, 0.2e1 * t59 * t43, 0.2e1 * t102 * t45 + 0.2e1 * t117 * t43 + 0.2e1 * t24 * t152 + 0.2e1 * t6 * t59, -0.2e1 * t101 * t45 - 0.2e1 * t146 * t43 + 0.2e1 * t24 * t151 + 0.2e1 * t5 * t59, -0.2e1 * (-t14 * t88 - t16 * t85) * t42 + 0.2e1 * (-t2 * t88 - t4 * t85 + (t14 * t85 - t16 * t88) * qJD(5)) * t60, 0.2e1 * t11 * t29 + 0.2e1 * t14 * t2 + 0.2e1 * t16 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t157, 0, 0, 0, 0, 0, -t13, t92, 0, 0, 0, 0, 0, t10, t9, t1, -t103 * t41 + t13 * t66 + t26 * t40 + t30 * t62 + t54 * t8 - t55 * t7; 0, 0, 0, 0, 0, 0, t139, -t140, 0, -pkin(8) * t139, pkin(8) * t140, 0, 0, -t42, -t43, 0, -t24, t23, t19, t15, t28, -t100, 0, t38 + (-t106 * qJD(5) - t24) * t88 + t95 * t85, t106 * t135 + t95 * t88 + t147, t3 + (-t41 * t60 + t42 * t54 + (-t55 * t60 - t14) * qJD(5)) * t88 + (-t40 * t60 + t42 * t55 - t2 + (t54 * t60 - t16) * qJD(5)) * t85, t11 * t66 + t14 * t41 + t16 * t40 + t2 * t54 + t29 * t62 + t4 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t128, -0.2e1 * t127, t72, t58, 0, 0, 0, 0.2e1 * t104, 0.2e1 * t145, -0.2e1 * t41 * t85 + 0.2e1 * t36 + 0.2e1 * (-t54 * t88 - t55 * t85) * qJD(5), 0.2e1 * t40 * t55 + 0.2e1 * t41 * t54 + 0.2e1 * t62 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t92, 0, 0, 0, 0, 0, t10, t9, t1, pkin(5) * t125 - t103 * t50 + t13 * t76 + t26 * t49 + t67 * t8 - t68 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t43, 0, -t24, t23, t19, t15, t28, -t100, 0, t38 + (pkin(4) * t42 - pkin(10) * t43) * t85 + (-t24 + (-pkin(4) * t60 - pkin(10) * t59) * qJD(5)) * t88, t101 * pkin(4) + pkin(10) * t100 + t147, t3 + (t42 * t67 - t50 * t60 + (-t60 * t68 - t14) * qJD(5)) * t88 + (t42 * t68 - t49 * t60 - t2 + (t60 * t67 - t16) * qJD(5)) * t85, t11 * t76 + t14 * t50 + t16 * t49 + t2 * t67 + t29 * t78 + t4 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t127, t72, t58, 0, 0, 0, t104 - t130, -t129 + t145, t36 + t48 + (-t41 - t50) * t85 + ((-t54 - t67) * t88 + (-t55 - t68) * t85) * qJD(5), t40 * t68 + t41 * t67 + t49 * t55 + t50 * t54 + t62 * t76 + t66 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t58, 0, 0, 0, -0.2e1 * t130, -0.2e1 * t129, -0.2e1 * t50 * t85 + 0.2e1 * t48 + 0.2e1 * (-t67 * t88 - t68 * t85) * qJD(5), 0.2e1 * t49 * t68 + 0.2e1 * t50 * t67 + 0.2e1 * t76 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t102, t43, t6, t5, t101 * pkin(5), t2 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t135, 0, -t85 * t127 - t74 * t80, t135 * t74 - t112, -t126, t41 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t135, 0, -pkin(10) * t80, pkin(10) * t135, -t126, t50 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
