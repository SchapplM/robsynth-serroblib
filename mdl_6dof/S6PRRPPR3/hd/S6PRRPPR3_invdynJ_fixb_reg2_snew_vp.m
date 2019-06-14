% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:00:33
% EndTime: 2019-05-05 03:00:40
% DurationCPUTime: 2.42s
% Computational Cost: add. (4528->281), mult. (9333->345), div. (0->0), fcn. (5467->10), ass. (0->167)
t139 = sin(qJ(3));
t142 = cos(qJ(3));
t145 = qJD(2) ^ 2;
t221 = t142 * t145;
t115 = t139 * t221;
t106 = qJDD(3) - t115;
t144 = qJD(3) ^ 2;
t128 = t139 ^ 2;
t195 = t128 * t145;
t110 = t144 + t195;
t133 = sin(pkin(6));
t134 = cos(pkin(6));
t140 = sin(qJ(2));
t143 = cos(qJ(2));
t189 = t142 * t106;
t66 = -t139 * t110 + t189;
t182 = qJD(2) * qJD(3);
t174 = t142 * t182;
t180 = t139 * qJDD(2);
t94 = 0.2e1 * t174 + t180;
t28 = t134 * (t139 * t106 + t142 * t110) + (t140 * t66 + t143 * t94) * t133;
t214 = pkin(2) * t94 + pkin(8) * t66;
t105 = qJDD(3) + t115;
t129 = t142 ^ 2;
t194 = t129 * t145;
t111 = t144 + t194;
t192 = t139 * t105;
t64 = t142 * t111 + t192;
t122 = t139 * t182;
t179 = t142 * qJDD(2);
t97 = -0.2e1 * t122 + t179;
t233 = (t140 * t64 - t143 * t97) * t133 - t134 * (t142 * t105 - t139 * t111);
t232 = pkin(8) * t64;
t95 = t174 + t180;
t231 = t95 + t174;
t138 = sin(qJ(6));
t141 = cos(qJ(6));
t185 = qJD(2) * t142;
t89 = -t141 * qJD(3) + t138 * t185;
t90 = t138 * qJD(3) + t141 * t185;
t57 = t89 * t90;
t83 = qJDD(6) + t95;
t225 = -t57 + t83;
t230 = t138 * t225;
t229 = t141 * t225;
t201 = sin(pkin(10));
t202 = cos(pkin(10));
t157 = t201 * g(1) - t202 * g(2);
t188 = -g(3) + qJDD(1);
t228 = t133 * t188 + t134 * t157;
t226 = pkin(8) - qJ(5);
t224 = pkin(3) * t110 + qJ(4) * t106;
t186 = t128 + t129;
t100 = t186 * qJDD(2);
t101 = t186 * t145;
t53 = (t100 * t140 + t101 * t143) * t133;
t62 = t189 + t139 * (-t144 + t194);
t184 = t139 * qJD(2);
t104 = -qJD(3) * pkin(4) - qJ(5) * t184;
t220 = t104 * t184 + qJDD(5);
t72 = -t133 * t157 + t134 * t188;
t69 = t142 * t72;
t171 = qJDD(3) * pkin(3) + t144 * qJ(4) - qJDD(4) + t69;
t103 = -t202 * g(1) - t201 * g(2);
t49 = t143 * t103 + t228 * t140;
t43 = -t145 * pkin(2) + qJDD(2) * pkin(8) + t49;
t193 = t139 * qJ(4);
t166 = -t142 * pkin(3) - t193;
t92 = t166 * qJD(2);
t173 = qJD(2) * t92 + t43;
t23 = t173 * t139 - t171;
t81 = t89 ^ 2;
t82 = t90 ^ 2;
t116 = qJD(6) + t184;
t113 = t116 ^ 2;
t219 = pkin(3) + pkin(4);
t218 = pkin(4) + pkin(9);
t96 = -t122 + t179;
t217 = t96 * pkin(4);
t93 = pkin(8) * t100;
t216 = pkin(5) + qJ(4);
t33 = t139 * t72 + t142 * t43;
t215 = pkin(2) * t97 - t232;
t213 = pkin(2) * t101 + t93;
t212 = t116 * t89;
t167 = pkin(5) * t139 + pkin(9) * t142;
t181 = qJD(4) * qJD(3);
t124 = 0.2e1 * t181;
t168 = -t144 * pkin(3) + qJDD(3) * qJ(4) + t92 * t185 + t33;
t176 = 0.2e1 * qJD(5) * qJD(2);
t149 = pkin(4) * t194 + t96 * qJ(5) - qJD(3) * t104 + t142 * t176 - t168;
t17 = t124 - t149;
t12 = qJDD(3) * pkin(5) - t144 * pkin(9) - t167 * t221 + t17;
t211 = t138 * t12;
t47 = t57 + t83;
t210 = t138 * t47;
t161 = t141 * qJDD(3) - t138 * t96;
t153 = (qJD(6) - t116) * t90 - t161;
t51 = t89 * qJD(6) - t138 * qJDD(3) - t141 * t96;
t40 = t51 - t212;
t21 = t138 * t40 + t141 * t153;
t209 = t139 * t21;
t208 = t139 * t43;
t207 = t141 * t12;
t206 = t141 * t47;
t205 = t142 * t94;
t200 = qJ(4) * t101;
t199 = qJ(4) * t111;
t198 = qJ(4) * t142;
t197 = t116 * t138;
t196 = t116 * t141;
t187 = -0.2e1 * qJD(5) + t92;
t183 = pkin(3) + t218;
t178 = -t82 - t113;
t177 = t139 * t57;
t175 = qJ(5) * qJD(3) * t142;
t121 = 0.2e1 * qJD(4) * t184;
t169 = t140 * t103 - t228 * t143;
t42 = -qJDD(2) * pkin(2) - t145 * pkin(8) + t169;
t154 = -t96 * pkin(3) - t231 * qJ(4) + t42;
t148 = qJ(5) * t194 + t154 - t220;
t10 = t121 + t218 * t96 + (pkin(5) * t142 + (-pkin(3) - pkin(9)) * t139) * t182 - t148 + t95 * pkin(5);
t159 = -pkin(4) * t105 - t95 * qJ(5) - t171;
t152 = t159 + t208;
t13 = -t144 * pkin(5) - qJDD(3) * pkin(9) + (t175 + (-t167 * qJD(2) + t187) * t139) * qJD(2) + t152;
t5 = -t141 * t10 + t138 * t13;
t32 = -t69 + t208;
t14 = t139 * t32 + t142 * t33;
t170 = -t142 * t219 - pkin(2);
t6 = t138 * t10 + t141 * t13;
t2 = t138 * t6 - t141 * t5;
t3 = t138 * t5 + t141 * t6;
t1 = t142 * t12 + t139 * t3;
t18 = (t187 * t139 + t175) * qJD(2) + t152;
t7 = t139 * t18 + t142 * t17;
t55 = t139 * t97 + t205;
t162 = qJD(2) * (pkin(3) * qJD(3) - 0.2e1 * qJD(4));
t22 = t124 + t168;
t160 = t139 * t162;
t158 = t51 + t212;
t156 = t139 * t216 + t142 * t183 + pkin(2);
t151 = t121 - t154;
t147 = t148 - t217;
t146 = -pkin(3) * t122 + qJ(4) * t94 + t151;
t102 = (t128 - t129) * t145;
t71 = -t82 + t113;
t70 = t81 - t113;
t63 = t192 + t142 * (t144 - t195);
t61 = t231 * t139;
t60 = (t96 - t122) * t142;
t56 = t82 - t81;
t52 = -t113 - t81;
t50 = t90 * qJD(6) - t161;
t45 = -t81 - t82;
t35 = (-qJD(6) - t116) * t90 + t161;
t31 = -t138 * t178 - t206;
t30 = t141 * t178 - t210;
t26 = t141 * t52 - t230;
t25 = t138 * t52 + t229;
t24 = t160 + t154;
t20 = t138 * t153 - t141 * t40;
t19 = t160 + t147;
t16 = t139 * t31 + t142 * t158;
t15 = t139 * t26 + t142 * t35;
t11 = t142 * t45 + t209;
t8 = t139 * t23 + t142 * t22;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, 0, 0, 0, 0, 0, (qJDD(2) * t143 - t140 * t145) * t133, (-qJDD(2) * t140 - t143 * t145) * t133, 0, t134 * t72 + (t140 * t49 - t143 * t169) * t133, 0, 0, 0, 0, 0, 0, -t233, -t28, t53, t134 * (t139 * t33 - t142 * t32) + (t140 * t14 - t143 * t42) * t133, 0, 0, 0, 0, 0, 0, -t233, t53, t28, t134 * (t139 * t22 - t142 * t23) + (t140 * t8 - t143 * t24) * t133, 0, 0, 0, 0, 0, 0, t28, t233, -t53, t134 * (t139 * t17 - t142 * t18) + (t140 * t7 - t143 * t19) * t133, 0, 0, 0, 0, 0, 0, t134 * (t139 * t35 - t142 * t26) + (t140 * t15 + t143 * t25) * t133, t134 * (t139 * t158 - t142 * t31) + (t140 * t16 + t143 * t30) * t133, t134 * (t139 * t45 - t142 * t21) + (t140 * t11 + t143 * t20) * t133, t134 * (t139 * t12 - t142 * t3) + (t140 * t1 + t143 * t2) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t169, -t49, 0, 0, t61, t55, t63, t60, t62, 0, -t142 * t42 + t215, t139 * t42 - t214, t14 + t213, -pkin(2) * t42 + pkin(8) * t14, t61, t63, -t55, 0, -t62, t60, t97 * t193 + t142 * ((t97 - t122) * pkin(3) + t151) + t215, t142 * (pkin(3) * t101 + t22) + (t23 + t200) * t139 + t213, pkin(3) * t205 + t139 * t146 + t214, pkin(8) * t8 + (-pkin(2) + t166) * t24, t60, t55, t62, t61, t63, 0, t139 * ((t110 - t194) * qJ(5) + t146 + t217 + t220) + t142 * (-qJ(5) * t106 + t219 * t94) + t214, t142 * (-qJ(5) * t111 + t147) + t232 + t170 * t97 + (-qJ(4) * t97 - qJ(5) * t105 + t142 * t162) * t139, t142 * (qJ(5) * t179 + t149 - 0.2e1 * t181) - t93 + t170 * t101 + (-qJ(5) * t174 - t159 - t200 + (qJ(5) * qJDD(2) - t173 + t176) * t139) * t139, (t170 - t193) * t19 + t226 * t7, t177 + t142 * (-t141 * t51 - t197 * t90), t139 * t56 + t142 * (t138 * t158 + t141 * t35), t139 * t40 + t142 * (t138 * t71 - t229), -t177 + t142 * (t138 * t50 + t196 * t89), t139 * t153 + t142 * (-t141 * t70 + t210), t139 * t83 + t142 * (t138 * t90 - t141 * t89) * t116, t139 * (-qJ(5) * t26 - t5) + t142 * (-qJ(5) * t35 - t211) + pkin(8) * t15 + t156 * t25, t139 * (-qJ(5) * t31 - t6) + t142 * (-qJ(5) * t158 - t207) + pkin(8) * t16 + t156 * t30, -qJ(5) * t209 + t142 * (-qJ(5) * t45 + t2) + pkin(8) * t11 + t156 * t20, t226 * t1 + t156 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, t102, t180, t115, t179, qJDD(3), -t32, -t33, 0, 0, -t115, t180, -t102, qJDD(3), -t179, t115, pkin(3) * t105 - t199 - t23, (-pkin(3) * t139 + t198) * qJDD(2), t22 + t224, -pkin(3) * t23 + qJ(4) * t22, t115, t102, t179, -t115, t180, qJDD(3), pkin(4) * t110 + t17 + t224, -t219 * t105 + t18 + t199, (t219 * t139 - t198) * qJDD(2), qJ(4) * t17 - t219 * t18, -t138 * t51 + t196 * t90, t138 * t35 - t141 * t158, -t141 * t71 - t230, -t141 * t50 + t197 * t89, -t138 * t70 - t206, (-t138 * t89 - t141 * t90) * t116, -t183 * t26 + t216 * t35 + t207, t158 * t216 - t183 * t31 - t211, -t183 * t21 + t216 * t45 - t3, t216 * t12 - t183 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t180, -t110, t23, 0, 0, 0, 0, 0, 0, -t110, t105, -t180, t18, 0, 0, 0, 0, 0, 0, t26, t31, t21, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t97, -t101, -t19, 0, 0, 0, 0, 0, 0, t25, t30, t20, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t56, t40, -t57, t153, t83, -t5, -t6, 0, 0;];
tauJ_reg  = t4;
