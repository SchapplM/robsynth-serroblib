% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 22:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:01:37
% EndTime: 2019-05-04 22:01:42
% DurationCPUTime: 1.46s
% Computational Cost: add. (5531->230), mult. (9844->324), div. (0->0), fcn. (6359->12), ass. (0->155)
t137 = sin(qJ(6));
t140 = cos(qJ(6));
t138 = sin(qJ(5));
t167 = qJD(2) * t138;
t92 = t140 * qJD(5) + t137 * t167;
t93 = -t137 * qJD(5) + t140 * t167;
t73 = t92 * t93;
t164 = qJD(2) * qJD(5);
t115 = t138 * t164;
t141 = cos(qJ(5));
t162 = t141 * qJDD(2);
t98 = t115 - t162;
t90 = qJDD(6) - t98;
t179 = -t73 + t90;
t182 = t137 * t179;
t181 = t140 * t179;
t126 = -g(3) + qJDD(1);
t131 = sin(pkin(6));
t130 = sin(pkin(10));
t133 = cos(pkin(10));
t103 = t130 * g(1) - t133 * g(2);
t134 = cos(pkin(6));
t170 = t134 * t103;
t180 = t126 * t131 + t170;
t113 = t141 * qJD(2) + qJD(6);
t160 = t141 * t164;
t163 = t138 * qJDD(2);
t96 = -t160 - t163;
t158 = -t140 * qJDD(5) + t137 * t96;
t47 = (qJD(6) - t113) * t93 - t158;
t88 = t92 ^ 2;
t89 = t93 ^ 2;
t111 = t113 ^ 2;
t144 = qJD(2) ^ 2;
t178 = pkin(2) + pkin(3);
t129 = sin(pkin(11));
t132 = cos(pkin(11));
t128 = qJDD(2) * pkin(2);
t104 = -t133 * g(1) - t130 * g(2);
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t64 = -t139 * t104 + t180 * t142;
t147 = -qJDD(3) + t64;
t55 = -t144 * qJ(3) - t128 - t147;
t145 = -qJDD(2) * pkin(3) + t55;
t122 = qJDD(2) * qJ(3);
t65 = t142 * t104 + t180 * t139;
t156 = 0.2e1 * qJD(3) * qJD(2) + t65;
t154 = t122 + t156;
t53 = -t178 * t144 + t154;
t39 = t129 * t145 + t132 * t53;
t143 = qJD(5) ^ 2;
t155 = pkin(5) * t141 + pkin(9) * t138;
t146 = t144 * t155;
t114 = t134 * t126;
t157 = -t131 * t103 + t114;
t150 = qJDD(4) - t157;
t37 = -t144 * pkin(4) - qJDD(2) * pkin(8) + t39;
t28 = t138 * t37 - t141 * t150;
t22 = -qJDD(5) * pkin(5) - t143 * pkin(9) - t138 * t146 + t28;
t177 = t137 * t22;
t62 = t73 + t90;
t176 = t137 * t62;
t175 = t140 * t22;
t174 = t140 * t62;
t173 = t113 * t137;
t172 = t113 * t140;
t112 = t138 * t144 * t141;
t105 = qJDD(5) + t112;
t169 = t138 * t105;
t106 = qJDD(5) - t112;
t168 = t141 * t106;
t165 = qJD(6) + t113;
t161 = t141 * t73;
t29 = t138 * t150 + t141 * t37;
t23 = -t143 * pkin(5) + qJDD(5) * pkin(9) - t141 * t146 + t29;
t152 = -t96 + t160;
t153 = -t98 - t115;
t159 = t129 * t53 - t132 * t145;
t36 = qJDD(2) * pkin(4) - t144 * pkin(8) + t159;
t27 = pkin(5) * t153 + pkin(9) * t152 + t36;
t10 = t137 * t23 - t140 * t27;
t11 = t137 * t27 + t140 * t23;
t6 = t137 * t10 + t140 * t11;
t151 = t140 * t10 - t137 * t11;
t13 = t138 * t28 + t141 * t29;
t149 = -t137 * qJDD(5) - t140 * t96;
t148 = pkin(4) + t155;
t67 = t92 * qJD(6) - t149;
t125 = t141 ^ 2;
t124 = t138 ^ 2;
t120 = t125 * t144;
t119 = t124 * t144;
t110 = -t120 - t143;
t109 = -t119 - t143;
t102 = t119 + t120;
t101 = (-t124 - t125) * qJDD(2);
t100 = t132 * qJDD(2) + t129 * t144;
t99 = -t129 * qJDD(2) + t132 * t144;
t97 = 0.2e1 * t115 - t162;
t95 = 0.2e1 * t160 + t163;
t85 = (qJDD(2) * t142 - t139 * t144) * t131;
t84 = (qJDD(2) * t139 + t142 * t144) * t131;
t82 = t113 * t92;
t81 = -t89 + t111;
t80 = t88 - t111;
t79 = -t138 * t109 - t168;
t78 = t141 * t110 - t169;
t77 = -t138 * t106 + t141 * t109;
t76 = t141 * t105 + t138 * t110;
t75 = t134 * t157;
t72 = t89 - t88;
t71 = -t89 - t111;
t70 = t132 * t101 - t129 * t102;
t69 = t129 * t101 + t132 * t102;
t68 = -t111 - t88;
t66 = t93 * qJD(6) - t158;
t60 = t88 + t89;
t59 = -t129 * t95 + t132 * t79;
t58 = -t129 * t97 + t132 * t78;
t57 = t129 * t79 + t132 * t95;
t56 = t129 * t78 + t132 * t97;
t54 = -t144 * pkin(2) + t154;
t50 = -t165 * t92 + t149;
t49 = t67 - t82;
t48 = t67 + t82;
t46 = t165 * t93 - t158;
t43 = -t137 * t71 - t174;
t42 = t140 * t71 - t176;
t41 = t140 * t68 - t182;
t40 = t137 * t68 + t181;
t35 = t137 * t49 + t140 * t47;
t34 = t137 * t47 - t140 * t49;
t33 = -t138 * t50 + t141 * t43;
t32 = t138 * t43 + t141 * t50;
t31 = -t138 * t46 + t141 * t41;
t30 = t138 * t41 + t141 * t46;
t25 = -t138 * t60 + t141 * t35;
t24 = t138 * t35 + t141 * t60;
t21 = t129 * t42 + t132 * t33;
t20 = t129 * t33 - t132 * t42;
t19 = t129 * t159 + t132 * t39;
t18 = t129 * t39 - t132 * t159;
t17 = t129 * t40 + t132 * t31;
t16 = t129 * t31 - t132 * t40;
t15 = t129 * t34 + t132 * t25;
t14 = t129 * t25 - t132 * t34;
t12 = t138 * t29 - t141 * t28;
t8 = t129 * t36 + t132 * t13;
t7 = t129 * t13 - t132 * t36;
t4 = t138 * t22 + t141 * t6;
t3 = t138 * t6 - t141 * t22;
t2 = -t129 * t151 + t132 * t4;
t1 = t129 * t4 + t132 * t151;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, 0, t85, -t84, 0, t75 + (t139 * t65 + t142 * t64) * t131, 0, 0, 0, 0, 0, 0, t85, 0, t84, t75 + (t139 * t54 - t142 * t55) * t131, 0, 0, 0, 0, 0, 0, (t100 * t142 - t139 * t99) * t131, (t100 * t139 + t142 * t99) * t131, 0, t134 * (-qJDD(4) + t114) + (t139 * t19 - t142 * t18 - t170) * t131, 0, 0, 0, 0, 0, 0, -t134 * t76 + (t139 * t58 - t142 * t56) * t131, -t134 * t77 + (t139 * t59 - t142 * t57) * t131, (t139 * t70 - t142 * t69) * t131, -t134 * t12 + (t139 * t8 - t142 * t7) * t131, 0, 0, 0, 0, 0, 0, -t134 * t30 + (t139 * t17 - t142 * t16) * t131, -t134 * t32 + (t139 * t21 - t142 * t20) * t131, -t134 * t24 + (t139 * t15 - t14 * t142) * t131, -t134 * t3 + (-t1 * t142 + t139 * t2) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t64, -t65, 0, 0, 0, 0, 0, qJDD(2), 0, 0, 0.2e1 * t128 + t147, 0, 0.2e1 * t122 + t156, -pkin(2) * t55 + qJ(3) * t54, 0, 0, 0, 0, 0, qJDD(2), -qJ(3) * t99 + t178 * t100 + t159, qJ(3) * t100 + t178 * t99 + t39, 0, qJ(3) * t19 - t178 * t18, t152 * t138, -t138 * t97 + t141 * t95, -t169 - t141 * (-t119 + t143), t153 * t141, -t138 * (t120 - t143) - t168, 0, -pkin(4) * t97 - pkin(8) * t78 + qJ(3) * t58 + t141 * t36 - t178 * t56, -pkin(4) * t95 - pkin(8) * t79 + qJ(3) * t59 - t138 * t36 - t178 * t57, -pkin(4) * t102 - pkin(8) * t101 + qJ(3) * t70 - t178 * t69 - t13, pkin(4) * t36 - pkin(8) * t13 + qJ(3) * t8 - t178 * t7, -t138 * (t140 * t67 + t93 * t173) + t161, -t138 * (-t137 * t48 + t140 * t46) + t141 * t72, -t138 * (-t137 * t81 + t181) + t141 * t49, -t138 * (-t137 * t66 - t92 * t172) - t161, -t138 * (t140 * t80 - t176) + t141 * t47, t141 * t90 - t138 * (-t137 * t93 + t140 * t92) * t113, qJ(3) * t17 - pkin(8) * t31 - t138 * (-pkin(9) * t40 + t177) - t141 * (-pkin(5) * t40 + t10) + pkin(4) * t40 - t178 * t16, qJ(3) * t21 - pkin(8) * t33 - t138 * (-pkin(9) * t42 + t175) - t141 * (-pkin(5) * t42 + t11) + pkin(4) * t42 - t178 * t20, -pkin(8) * t25 + qJ(3) * t15 - t138 * t151 - t178 * t14 + t148 * t34, -pkin(8) * t4 + qJ(3) * t2 - t178 * t1 - t148 * t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t144, t55, 0, 0, 0, 0, 0, 0, -t100, -t99, 0, t18, 0, 0, 0, 0, 0, 0, t56, t57, t69, t7, 0, 0, 0, 0, 0, 0, t16, t20, t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, 0, 0, 0, 0, 0, 0, t76, t77, 0, t12, 0, 0, 0, 0, 0, 0, t30, t32, t24, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t119 - t120, -t163, t112, -t162, qJDD(5), -t28, -t29, 0, 0, t137 * t67 - t93 * t172, t137 * t46 + t140 * t48, t140 * t81 + t182, t140 * t66 - t92 * t173, t137 * t80 + t174, (t137 * t92 + t140 * t93) * t113, pkin(5) * t46 + pkin(9) * t41 - t175, pkin(5) * t50 + pkin(9) * t43 + t177, pkin(5) * t60 + pkin(9) * t35 + t6, -pkin(5) * t22 + pkin(9) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, t49, -t73, t47, t90, -t10, -t11, 0, 0;];
tauJ_reg  = t5;
