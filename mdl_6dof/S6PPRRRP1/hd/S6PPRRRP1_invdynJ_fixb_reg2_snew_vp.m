% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPRRRP1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPRRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:29:04
% EndTime: 2019-05-04 20:29:12
% DurationCPUTime: 2.59s
% Computational Cost: add. (10073->280), mult. (18584->406), div. (0->0), fcn. (14888->14), ass. (0->202)
t166 = sin(qJ(4));
t169 = cos(qJ(4));
t239 = pkin(4) * t169;
t199 = -pkin(10) * t166 - t239;
t132 = t199 * qJD(3);
t247 = qJD(4) ^ 2;
t171 = qJD(3) ^ 2;
t167 = sin(qJ(3));
t170 = cos(qJ(3));
t158 = sin(pkin(7));
t161 = cos(pkin(7));
t157 = sin(pkin(12));
t160 = cos(pkin(12));
t159 = sin(pkin(6));
t162 = cos(pkin(6));
t225 = sin(pkin(11));
t226 = cos(pkin(11));
t180 = g(1) * t225 - g(2) * t226;
t215 = -g(3) + qJDD(1);
t177 = t159 * t215 + t162 * t180;
t181 = -g(1) * t226 - g(2) * t225;
t174 = -t157 * t181 + t160 * t177;
t176 = -t159 * t180 + t162 * t215 + qJDD(2);
t253 = t158 * t176 + t161 * t174;
t92 = t157 * t177 + t160 * t181;
t58 = t253 * t167 + t170 * t92;
t54 = -t171 * pkin(3) + qJDD(3) * pkin(9) + t58;
t172 = -t158 * t174 + t161 * t176;
t75 = t169 * t172;
t33 = (qJD(3) * t132 + t54) * t166 - qJDD(4) * pkin(4) - t247 * pkin(10) - t75;
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t213 = qJD(3) * t166;
t129 = -t168 * qJD(4) + t165 * t213;
t131 = qJD(4) * t165 + t168 * t213;
t108 = t131 * t129;
t208 = qJD(3) * qJD(4);
t148 = t166 * t208;
t206 = t169 * qJDD(3);
t135 = -t148 + t206;
t128 = -qJDD(5) + t135;
t249 = -t108 - t128;
t257 = pkin(5) * t249;
t203 = t169 * t208;
t207 = t166 * qJDD(3);
t134 = t203 + t207;
t101 = -qJD(5) * t129 + qJDD(4) * t165 + t134 * t168;
t212 = qJD(3) * t169;
t145 = -qJD(5) + t212;
t118 = t129 * t145;
t88 = t101 - t118;
t256 = qJ(6) * t88;
t255 = t165 * t249;
t254 = t168 * t249;
t127 = t131 ^ 2;
t143 = t145 ^ 2;
t105 = -t127 - t143;
t126 = t129 ^ 2;
t201 = -t168 * qJDD(4) + t165 * t134;
t100 = -qJD(5) * t131 - t201;
t113 = -pkin(5) * t145 - qJ(6) * t131;
t37 = t166 * t172 + t169 * t54;
t34 = -pkin(4) * t247 + qJDD(4) * pkin(10) + t132 * t212 + t37;
t195 = -t135 + t148;
t196 = t134 + t203;
t200 = t167 * t92 - t253 * t170;
t53 = -qJDD(3) * pkin(3) - t171 * pkin(9) + t200;
t41 = pkin(4) * t195 - pkin(10) * t196 + t53;
t26 = t165 * t41 + t168 * t34;
t183 = t100 * qJ(6) - 0.2e1 * qJD(6) * t129 + t145 * t113 + t26;
t252 = -t183 + (t105 + t126) * pkin(5);
t31 = -t100 * pkin(5) - t126 * qJ(6) + t131 * t113 + qJDD(6) + t33;
t250 = t101 + t118;
t248 = t167 * t58 - t170 * t200;
t84 = (qJD(5) + t145) * t131 + t201;
t102 = -t143 - t126;
t67 = t102 * t165 + t254;
t246 = pkin(4) * t67;
t94 = -t108 + t128;
t231 = t165 * t94;
t73 = t105 * t168 + t231;
t245 = pkin(4) * t73;
t210 = qJD(6) * t131;
t122 = -0.2e1 * t210;
t25 = t165 * t34 - t168 * t41;
t182 = -t25 - t256 + t257;
t18 = t122 + t182;
t244 = pkin(5) * t18;
t243 = pkin(5) * t88;
t60 = -t165 * t84 - t168 * t88;
t242 = pkin(10) * t60;
t241 = pkin(10) * t67;
t240 = pkin(10) * t73;
t61 = t165 * t88 - t168 * t84;
t93 = -t126 - t127;
t43 = t166 * t93 + t169 * t61;
t238 = -pkin(3) * t60 + pkin(9) * t43;
t68 = t102 * t168 - t255;
t83 = (qJD(5) - t145) * t131 + t201;
t47 = t166 * t83 + t169 * t68;
t237 = -pkin(3) * t67 + pkin(9) * t47;
t228 = t168 * t94;
t74 = -t105 * t165 + t228;
t50 = t166 * t250 + t169 * t74;
t236 = -pkin(3) * t73 + pkin(9) * t50;
t235 = -pkin(4) * t93 + pkin(10) * t61;
t234 = -pkin(4) * t83 + pkin(10) * t68;
t233 = -pkin(4) * t250 + pkin(10) * t74;
t232 = t165 * t33;
t229 = t168 * t33;
t224 = qJ(6) * t165;
t223 = qJ(6) * t168;
t222 = t145 * t165;
t221 = t145 * t168;
t220 = t157 * t159;
t219 = t159 * t160;
t218 = t160 * t161;
t144 = t166 * t171 * t169;
t139 = qJDD(4) + t144;
t217 = t166 * t139;
t140 = qJDD(4) - t144;
t216 = t169 * t140;
t204 = t169 * t108;
t15 = t165 * t25 + t168 * t26;
t36 = t166 * t54 - t75;
t22 = t166 * t36 + t169 * t37;
t19 = -pkin(5) * t126 + t183;
t10 = -t165 * t18 + t168 * t19;
t4 = t10 * t169 + t166 * t31;
t9 = t165 * t19 + t168 * t18;
t198 = t167 * t4 - t170 * t9;
t14 = t165 * t26 - t168 * t25;
t7 = t15 * t169 + t166 * t33;
t197 = -t14 * t170 + t167 * t7;
t194 = t167 * t22 - t170 * t53;
t193 = t167 * t43 - t170 * t60;
t192 = t167 * t47 - t170 * t67;
t191 = t167 * t50 - t170 * t73;
t154 = t169 ^ 2;
t152 = t154 * t171;
t142 = -t152 - t247;
t111 = t142 * t169 - t217;
t136 = -0.2e1 * t148 + t206;
t189 = t111 * t167 + t136 * t170;
t153 = t166 ^ 2;
t150 = t153 * t171;
t141 = -t150 - t247;
t112 = -t141 * t166 - t216;
t133 = 0.2e1 * t203 + t207;
t188 = t112 * t167 - t133 * t170;
t137 = (t153 + t154) * qJDD(3);
t138 = t150 + t152;
t187 = t137 * t167 + t138 * t170;
t185 = qJDD(3) * t170 - t167 * t171;
t184 = -qJDD(3) * t167 - t170 * t171;
t178 = t182 + t257;
t123 = 0.2e1 * t210;
t120 = t185 * t158;
t119 = t184 * t158;
t115 = -t127 + t143;
t114 = t126 - t143;
t110 = -t140 * t166 + t141 * t169;
t109 = t139 * t169 + t142 * t166;
t106 = t127 - t126;
t103 = t187 * t158;
t90 = (t129 * t165 + t131 * t168) * t145;
t80 = t101 * t165 - t131 * t221;
t79 = t100 * t168 - t129 * t222;
t78 = t169 * t128 + t166 * (t129 * t168 - t131 * t165) * t145;
t77 = t114 * t165 - t228;
t76 = t115 * t168 + t255;
t70 = t161 * t110 + t158 * t188;
t69 = t161 * t109 + t158 * t189;
t64 = t166 * (t101 * t168 + t131 * t222) - t204;
t63 = t166 * (-t100 * t165 - t129 * t221) + t204;
t62 = -pkin(5) * t250 + qJ(6) * t94;
t59 = -t165 * t83 + t168 * t250;
t52 = t166 * (t114 * t168 + t231) + t169 * t84;
t51 = t166 * (-t115 * t165 + t254) - t169 * t88;
t49 = t166 * t74 - t169 * t250;
t46 = t166 * t68 - t169 * t83;
t44 = t166 * (-t165 * t250 - t168 * t83) - t169 * t106;
t42 = t166 * t61 - t169 * t93;
t30 = t248 * t158 + t161 * t172;
t29 = -qJ(6) * t105 + t31;
t28 = t158 * t191 + t161 * t49;
t27 = t158 * t192 + t161 * t46;
t23 = -pkin(5) * t83 + qJ(6) * t102 - t31;
t21 = t166 * t37 - t169 * t36;
t20 = t158 * t193 + t161 * t42;
t17 = t123 - t182 + t256;
t16 = -qJ(6) * t84 + (-t126 - t93) * pkin(5) + t183;
t13 = -pkin(5) * t31 + qJ(6) * t19;
t12 = t158 * t194 + t161 * t21;
t11 = t162 * t28 + (t157 * (t167 * t73 + t170 * t50) + t160 * (-t158 * t49 + t161 * t191)) * t159;
t8 = t162 * t27 + (t157 * (t167 * t67 + t170 * t47) + t160 * (-t158 * t46 + t161 * t192)) * t159;
t6 = t15 * t166 - t169 * t33;
t5 = t162 * t20 + (t157 * (t167 * t60 + t170 * t43) + t160 * (-t158 * t42 + t161 * t193)) * t159;
t3 = t10 * t166 - t169 * t31;
t2 = t158 * t197 + t161 * t6;
t1 = t158 * t198 + t161 * t3;
t24 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t215, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162 * t176 + t174 * t219 + t220 * t92, 0, 0, 0, 0, 0, 0, t162 * t120 + (t157 * t184 + t185 * t218) * t159, t162 * t119 + (-t157 * t185 + t184 * t218) * t159, 0, (t167 * t200 + t170 * t58) * t220 + (-t158 * t172 + t248 * t161) * t219 + t162 * t30, 0, 0, 0, 0, 0, 0, t162 * t69 + (t157 * (t111 * t170 - t136 * t167) + t160 * (-t158 * t109 + t161 * t189)) * t159, t162 * t70 + (t157 * (t112 * t170 + t133 * t167) + t160 * (-t158 * t110 + t161 * t188)) * t159, t162 * t103 + (t157 * (t137 * t170 - t138 * t167) + t187 * t218) * t159, t162 * t12 + (t157 * (t167 * t53 + t170 * t22) + t160 * (-t158 * t21 + t161 * t194)) * t159, 0, 0, 0, 0, 0, 0, t8, t11, t5, t162 * t2 + (t157 * (t14 * t167 + t170 * t7) + t160 * (-t158 * t6 + t161 * t197)) * t159, 0, 0, 0, 0, 0, 0, t8, t11, t5, t162 * t1 + (t157 * (t167 * t9 + t170 * t4) + t160 * (-t158 * t3 + t161 * t198)) * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, 0, 0, 0, 0, 0, 0, t120, t119, 0, t30, 0, 0, 0, 0, 0, 0, t69, t70, t103, t12, 0, 0, 0, 0, 0, 0, t27, t28, t20, t2, 0, 0, 0, 0, 0, 0, t27, t28, t20, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t200, -t58, 0, 0, t196 * t166, t133 * t169 + t136 * t166, t217 + t169 * (-t150 + t247), -t195 * t169, t166 * (t152 - t247) + t216, 0, pkin(3) * t136 + pkin(9) * t111 - t169 * t53, -pkin(3) * t133 + pkin(9) * t112 + t166 * t53, pkin(3) * t138 + pkin(9) * t137 + t22, -pkin(3) * t53 + pkin(9) * t22, t64, t44, t51, t63, t52, t78, t166 * (t232 - t241) + t169 * (t25 - t246) + t237, t166 * (t229 - t240) + t169 * (t26 - t245) + t236, t166 * (-t14 - t242) - t60 * t239 + t238, pkin(9) * t7 + (-pkin(3) + t199) * t14, t64, t44, t51, t63, t52, t78, t166 * (-t165 * t23 - t223 * t249 - t241) + t169 * (t123 - t178 - t246) + t237, t166 * (-t165 * t62 + t168 * t29 - t240) + t169 * (-t245 - t252) + t236, t166 * (-t16 * t165 + t168 * t17 - t242) + t169 * (-pkin(4) * t60 + t243) + t238, t166 * (-pkin(10) * t9 - t13 * t165 - t18 * t223) + t169 * (-pkin(4) * t9 - t244) - pkin(3) * t9 + pkin(9) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t150 - t152, t207, t144, t206, qJDD(4), -t36, -t37, 0, 0, t80, t59, t76, t79, t77, t90, -t229 + t234, t232 + t233, t15 + t235, -pkin(4) * t33 + pkin(10) * t15, t80, t59, t76, t79, t77, t90, t168 * t23 - t224 * t249 + t234, t165 * t29 + t168 * t62 + t233, t16 * t168 + t165 * t17 + t235, -pkin(4) * t31 + pkin(10) * t10 + t13 * t168 - t18 * t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t106, t88, -t108, -t84, -t128, -t25, -t26, 0, 0, t108, t106, t88, -t108, -t84, -t128, t122 + t178, t252, -t243, t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t250, t93, t31;];
tauJ_reg  = t24;
