% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRP2
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
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:25
% EndTime: 2019-03-09 00:04:30
% DurationCPUTime: 1.87s
% Computational Cost: add. (2795->231), mult. (6867->384), div. (0->0), fcn. (6635->10), ass. (0->141)
t150 = sin(pkin(6));
t124 = t150 * cos(qJ(2));
t110 = qJD(2) * t124;
t87 = sin(qJ(3));
t123 = t150 * sin(qJ(2));
t151 = cos(pkin(6));
t90 = cos(qJ(3));
t98 = t90 * t123 + t151 * t87;
t93 = t98 * qJD(3) + t87 * t110;
t179 = qJD(4) * t98 + t93;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t64 = t86 * t87 - t89 * t90;
t65 = t86 * t90 + t89 * t87;
t80 = -t90 * pkin(3) - pkin(2);
t39 = t64 * pkin(4) - t65 * pkin(10) + t80;
t171 = -pkin(9) - pkin(8);
t73 = t171 * t87;
t74 = t171 * t90;
t47 = t86 * t73 - t89 * t74;
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t178 = t85 * t39 + t88 * t47;
t99 = t87 * t123 - t151 * t90;
t174 = t99 * qJD(3) - t90 * t110;
t177 = -qJD(4) * t99 - t174;
t83 = t85 ^ 2;
t84 = t88 ^ 2;
t153 = t83 - t84;
t127 = t153 * qJD(5);
t176 = qJD(3) + qJD(4);
t144 = qJD(5) * t85;
t43 = t176 * t64;
t159 = t88 * t43;
t107 = t65 * t144 + t159;
t14 = t177 * t86 + t179 * t89;
t162 = t65 * t88;
t112 = t85 * t124;
t36 = -t86 * t99 + t89 * t98;
t31 = t88 * t36 - t112;
t35 = t86 * t98 + t89 * t99;
t44 = t176 * t65;
t109 = qJD(2) * t123;
t30 = t88 * t124 + t85 * t36;
t91 = -t177 * t89 + t179 * t86;
t8 = t30 * qJD(5) - t85 * t109 + t88 * t91;
t175 = t107 * t35 - t14 * t162 + t31 * t44 - t8 * t64;
t149 = qJD(3) * t87;
t137 = pkin(3) * t149;
t22 = t44 * pkin(4) + t43 * pkin(10) + t137;
t130 = qJD(3) * t171;
t125 = t90 * t130;
t126 = t87 * t130;
t146 = qJD(4) * t89;
t147 = qJD(4) * t86;
t27 = -t86 * t125 - t89 * t126 - t73 * t146 - t74 * t147;
t6 = -qJD(5) * t178 + t88 * t22 + t85 * t27;
t122 = t88 * pkin(5) + t85 * qJ(6);
t173 = t122 * qJD(5) - t88 * qJD(6);
t172 = 0.2e1 * qJD(6);
t170 = pkin(10) * t44;
t169 = pkin(10) * t64;
t168 = t44 * pkin(5);
t167 = t89 * pkin(3);
t78 = t86 * pkin(3) + pkin(10);
t166 = t44 * t78;
t165 = t64 * t78;
t164 = t65 * t43;
t163 = t65 * t85;
t161 = t85 * t43;
t160 = t85 * t44;
t158 = t88 * t44;
t28 = t47 * qJD(4) - t89 * t125 + t86 * t126;
t46 = -t89 * t73 - t86 * t74;
t82 = qJD(5) * t88;
t157 = t28 * t85 + t46 * t82;
t136 = pkin(3) * t147;
t142 = t85 * qJD(6);
t55 = pkin(5) * t144 - qJ(6) * t82 - t142;
t48 = t55 + t136;
t155 = -t48 - t55;
t79 = -pkin(4) - t167;
t154 = t85 * t136 + t79 * t82;
t152 = t44 * qJ(6);
t148 = qJD(3) * t90;
t121 = t85 * pkin(5) - t88 * qJ(6);
t26 = t121 * t65 + t46;
t145 = qJD(5) * t26;
t143 = t64 * qJD(6);
t140 = -0.2e1 * pkin(2) * qJD(3);
t139 = pkin(4) * t144;
t138 = pkin(4) * t82;
t135 = pkin(3) * t146;
t134 = pkin(10) * t144;
t133 = pkin(10) * t82;
t132 = t35 * t82;
t131 = t85 * t82;
t129 = -0.4e1 * t85 * t162;
t16 = t64 * qJ(6) + t178;
t116 = t88 * t39 - t85 * t47;
t17 = -t64 * pkin(5) - t116;
t120 = t16 * t88 + t17 * t85;
t119 = -t16 * t85 + t17 * t88;
t118 = t30 * t88 - t31 * t85;
t117 = t30 * t85 + t31 * t88;
t114 = -t65 * t79 + t165;
t53 = (t83 + t84) * t135;
t113 = -t88 * t136 + t79 * t144;
t69 = -pkin(4) - t122;
t111 = qJD(3) * t124;
t108 = t65 * t82 - t161;
t106 = t64 * t144 - t158;
t5 = t47 * t144 - t85 * t22 + t88 * t27 - t39 * t82;
t104 = -t43 * t69 + t55 * t65 - t170;
t9 = -qJD(5) * t112 - t88 * t109 + t36 * t82 - t85 * t91;
t103 = t65 * t132 + t14 * t163 - t35 * t161 - t30 * t44 - t9 * t64;
t7 = -t121 * t43 + t173 * t65 + t28;
t102 = -t7 + (t65 * t69 - t169) * qJD(5);
t60 = t69 - t167;
t101 = -t7 + (t60 * t65 - t165) * qJD(5);
t97 = -t64 * t135 - t43 * t60 + t48 * t65 - t166;
t3 = t143 - t5 + t152;
t4 = -t168 - t6;
t1 = t119 * qJD(5) + t3 * t88 + t4 * t85;
t2 = t118 * qJD(5) - t8 * t88 + t9 * t85;
t94 = -t43 * t79 - t166 + (-t64 * t89 + t65 * t86) * qJD(4) * pkin(3);
t76 = 0.2e1 * t131;
t63 = -0.2e1 * t127;
t62 = t65 ^ 2;
t59 = t69 * t144;
t52 = t60 * t144;
t51 = t85 * t135 + t78 * t82;
t50 = -t88 * t135 + t78 * t144;
t41 = t46 * t144;
t34 = t64 * t82 + t160;
t23 = t26 * t144;
t21 = -t65 * t127 - t85 * t159;
t15 = qJD(5) * t129 + t153 * t43;
t11 = -t14 * t88 + t35 * t144;
t10 = t14 * t85 + t132;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t35 * t14 + 0.2e1 * t30 * t9 - 0.2e1 * t31 * t8; 0, 0, -t109, -t110, 0, 0, 0, 0, 0, -t90 * t109 - t87 * t111, t87 * t109 - t90 * t111, 0, 0, 0, 0, 0, t109 * t64 - t124 * t44, t109 * t65 + t124 * t43, 0, 0, 0, 0, 0, t103, -t175, t103, -t118 * t43 + (-qJD(5) * t117 + t8 * t85 + t88 * t9) * t65, t175, t14 * t26 - t8 * t16 + t9 * t17 + t31 * t3 + t30 * t4 + t35 * t7; 0, 0, 0, 0, 0.2e1 * t87 * t148, 0.2e1 * (-t87 ^ 2 + t90 ^ 2) * qJD(3), 0, 0, 0, t87 * t140, t90 * t140, -0.2e1 * t164, 0.2e1 * t43 * t64 - 0.2e1 * t65 * t44, 0, 0, 0, 0.2e1 * t137 * t64 + 0.2e1 * t80 * t44, 0.2e1 * t137 * t65 - 0.2e1 * t80 * t43, -0.2e1 * t131 * t62 - 0.2e1 * t84 * t164, 0.2e1 * t62 * t127 - t43 * t129, -0.2e1 * t107 * t64 + 0.2e1 * t65 * t158, -0.2e1 * t108 * t64 - 0.2e1 * t65 * t160, 0.2e1 * t64 * t44, 0.2e1 * t108 * t46 + 0.2e1 * t116 * t44 + 0.2e1 * t28 * t163 + 0.2e1 * t6 * t64, -0.2e1 * t107 * t46 + 0.2e1 * t28 * t162 - 0.2e1 * t178 * t44 + 0.2e1 * t5 * t64, 0.2e1 * t108 * t26 + 0.2e1 * t7 * t163 - 0.2e1 * t17 * t44 - 0.2e1 * t4 * t64, -0.2e1 * t119 * t43 + 0.2e1 * (-qJD(5) * t120 - t3 * t85 + t4 * t88) * t65, 0.2e1 * t107 * t26 + 0.2e1 * t16 * t44 - 0.2e1 * t7 * t162 + 0.2e1 * t3 * t64, 0.2e1 * t16 * t3 + 0.2e1 * t17 * t4 + 0.2e1 * t26 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t174, 0, 0, 0, 0, 0, -t14, t91, 0, 0, 0, 0, 0, t11, t10, t11, t2, -t10, t117 * t135 + t14 * t60 + t2 * t78 + t35 * t48; 0, 0, 0, 0, 0, 0, t148, -t149, 0, -pkin(8) * t148, pkin(8) * t149, 0, 0, -t43, -t44, 0, -t28, t27, t21, t15, t34, -t106, 0, t41 + (-qJD(5) * t114 - t28) * t88 + t94 * t85, t114 * t144 + t88 * t94 + t157, t101 * t88 + t85 * t97 + t23, t1, t101 * t85 + (-t97 - t145) * t88, t1 * t78 + t120 * t135 + t26 * t48 + t7 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t136, -0.2e1 * t135, t76, t63, 0, 0, 0, 0.2e1 * t113, 0.2e1 * t154, -0.2e1 * t48 * t88 + 0.2e1 * t52, 0.2e1 * t53, -0.2e1 * t48 * t85 - 0.2e1 * t60 * t82, 0.2e1 * t60 * t48 + 0.2e1 * t53 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t91, 0, 0, 0, 0, 0, t11, t10, t11, t2, -t10, pkin(10) * t2 + t14 * t69 + t35 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, -t28, t27, t21, t15, t34, -t106, 0, t41 + (pkin(4) * t43 - t170) * t85 + (-t28 + (-pkin(4) * t65 - t169) * qJD(5)) * t88, pkin(4) * t107 + pkin(10) * t106 + t157, t102 * t88 + t104 * t85 + t23, t1, t102 * t85 + (-t104 - t145) * t88, pkin(10) * t1 + t26 * t55 + t7 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t135, t76, t63, 0, 0, 0, t113 - t139, -t138 + t154, t155 * t88 + t52 + t59, t53, t155 * t85 + (-t60 - t69) * t82, pkin(10) * t53 + t48 * t69 + t60 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t63, 0, 0, 0, -0.2e1 * t139, -0.2e1 * t138, -0.2e1 * t55 * t88 + 0.2e1 * t59, 0, -0.2e1 * t55 * t85 - 0.2e1 * t69 * t82, 0.2e1 * t69 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, -t9, 0, -t8, -t9 * pkin(5) - t8 * qJ(6) + t31 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t108, t44, t6, t5, t6 + 0.2e1 * t168, t122 * t43 + (qJD(5) * t121 - t142) * t65, 0.2e1 * t143 - t5 + 0.2e1 * t152, -t4 * pkin(5) + t3 * qJ(6) + t16 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t144, 0, -t51, t50, -t51, -t173, -t50, -t121 * t135 - t173 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t144, 0, -t133, t134, -t133, -t173, -t134, -t173 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, qJ(6) * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t107, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
