% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:56
% EndTime: 2019-03-09 02:24:01
% DurationCPUTime: 2.92s
% Computational Cost: add. (2646->370), mult. (5035->502), div. (0->0), fcn. (3419->14), ass. (0->197)
t145 = sin(qJ(4));
t224 = qJD(1) * t145;
t116 = qJD(5) + t224;
t149 = cos(qJ(4));
t135 = qJ(1) + pkin(10);
t125 = sin(t135);
t126 = cos(t135);
t178 = g(1) * t125 - g(2) * t126;
t254 = g(3) * t145;
t158 = t178 * t149 - t254;
t211 = qJD(2) * qJD(4);
t222 = qJD(4) * t145;
t142 = cos(pkin(10));
t119 = -pkin(1) * t142 - pkin(2);
t115 = -pkin(7) + t119;
t88 = t115 * qJD(1) + qJD(3);
t199 = -t145 * qJDD(2) - t149 * t211 - t88 * t222;
t85 = t115 * qJDD(1) + qJDD(3);
t24 = -qJDD(4) * pkin(4) - t149 * t85 - t199;
t274 = qJD(5) * pkin(8) * t116 + t158 + t24;
t141 = sin(pkin(10));
t117 = pkin(1) * t141 + qJ(3);
t227 = qJDD(1) * t117;
t143 = sin(qJ(6));
t147 = cos(qJ(6));
t144 = sin(qJ(5));
t148 = cos(qJ(5));
t213 = t148 * qJD(4);
t223 = qJD(1) * t149;
t93 = t144 * t223 - t213;
t95 = qJD(4) * t144 + t148 * t223;
t174 = t143 * t93 - t147 * t95;
t41 = t143 * t95 + t147 * t93;
t273 = t174 * t41;
t272 = t174 ^ 2 - t41 ^ 2;
t113 = qJD(6) + t116;
t215 = qJD(6) * t147;
t216 = qJD(6) * t143;
t217 = qJD(5) * t149;
t194 = t144 * t217;
t163 = t145 * t213 + t194;
t205 = t149 * qJDD(1);
t37 = -qJD(1) * t163 + qJD(5) * t213 + t144 * qJDD(4) + t148 * t205;
t192 = t144 * t222;
t38 = -qJD(1) * t192 + qJD(5) * t95 - t148 * qJDD(4) + t144 * t205;
t7 = -t143 * t38 + t147 * t37 - t93 * t215 - t95 * t216;
t271 = t113 * t41 + t7;
t59 = t149 * qJD(2) + t145 * t88;
t52 = qJD(4) * pkin(8) + t59;
t179 = pkin(4) * t145 - pkin(8) * t149;
t84 = t117 + t179;
t60 = t84 * qJD(1);
t22 = t144 * t60 + t148 * t52;
t16 = -pkin(9) * t93 + t22;
t13 = t16 * t216;
t140 = qJ(5) + qJ(6);
t134 = cos(t140);
t253 = g(3) * t149;
t58 = -t145 * qJD(2) + t149 * t88;
t51 = -qJD(4) * pkin(4) - t58;
t29 = pkin(5) * t93 + t51;
t133 = sin(t140);
t231 = t134 * t145;
t55 = t125 * t231 + t126 * t133;
t57 = -t125 * t133 + t126 * t231;
t270 = g(1) * t55 - g(2) * t57 + t134 * t253 + t29 * t41 + t13;
t266 = qJD(4) * t58;
t23 = qJDD(4) * pkin(8) + qJDD(2) * t149 + t145 * t85 + t266;
t180 = pkin(4) * t149 + pkin(8) * t145;
t91 = t180 * qJD(4) + qJD(3);
t34 = t91 * qJD(1) + t179 * qJDD(1) + t227;
t28 = t148 * t34;
t212 = qJD(1) * qJD(4);
t190 = t149 * t212;
t206 = t145 * qJDD(1);
t90 = qJDD(5) + t190 + t206;
t2 = pkin(5) * t90 - pkin(9) * t37 - qJD(5) * t22 - t144 * t23 + t28;
t218 = qJD(5) * t148;
t204 = -t144 * t34 - t148 * t23 - t60 * t218;
t220 = qJD(5) * t144;
t168 = -t52 * t220 - t204;
t3 = -pkin(9) * t38 + t168;
t198 = -t143 * t3 + t147 * t2;
t21 = -t144 * t52 + t148 * t60;
t15 = -pkin(9) * t95 + t21;
t12 = pkin(5) * t116 + t15;
t239 = t147 * t16;
t5 = t143 * t12 + t239;
t232 = t133 * t145;
t54 = -t125 * t232 + t126 * t134;
t56 = t125 * t134 + t126 * t232;
t269 = -g(1) * t54 - g(2) * t56 - qJD(6) * t5 + t133 * t253 + t29 * t174 + t198;
t156 = qJD(6) * t174 - t143 * t37 - t147 * t38;
t268 = -t113 * t174 + t156;
t97 = t143 * t148 + t144 * t147;
t71 = t97 * t149;
t103 = qJD(1) * t117;
t265 = -qJD(1) * t103 - t178;
t264 = -qJD(6) * t148 - t218;
t263 = t119 * qJDD(1);
t262 = qJD(5) + qJD(6);
t228 = t148 * t149;
t261 = -t116 * t163 + t90 * t228;
t260 = 0.2e1 * t103 * qJD(4) + qJDD(4) * t115;
t259 = pkin(8) + pkin(9);
t221 = qJD(4) * t149;
t258 = t7 * t145 - t174 * t221;
t256 = pkin(9) * t149;
t83 = qJDD(6) + t90;
t96 = t143 * t144 - t147 * t148;
t252 = t83 * t96;
t251 = t83 * t97;
t193 = t144 * t216;
t18 = -t149 * t193 + (t262 * t228 - t192) * t147 - t163 * t143;
t250 = -t18 * t113 - t71 * t83;
t249 = t37 * t145 + t95 * t221;
t166 = t96 * t145;
t248 = -qJD(1) * t166 - t262 * t96;
t164 = qJD(1) * t97;
t247 = t145 * t164 + t262 * t97;
t98 = t180 * qJD(1);
t246 = t144 * t98 + t148 * t58;
t229 = t145 * t148;
t92 = t115 * t229;
t245 = t144 * t84 + t92;
t244 = t116 * t93;
t243 = t116 * t95;
t242 = t144 * t37;
t241 = t144 * t90;
t240 = t145 * t38;
t238 = t148 * t90;
t237 = t148 * t95;
t236 = t149 * t37;
t234 = qJD(4) * t93;
t233 = t115 * t144;
t230 = t144 * t145;
t139 = qJDD(2) - g(3);
t138 = t149 ^ 2;
t226 = t145 ^ 2 - t138;
t151 = qJD(4) ^ 2;
t152 = qJD(1) ^ 2;
t225 = -t151 - t152;
t219 = qJD(5) * t145;
t210 = qJD(3) * qJD(1);
t208 = qJDD(4) * t145;
t207 = qJDD(4) * t149;
t195 = t149 * t213;
t203 = t115 * t195 + t144 * t91 + t84 * t218;
t202 = pkin(9) * t229;
t197 = qJD(5) * t259;
t196 = t144 * t224;
t189 = pkin(5) - t233;
t188 = qJD(6) * t12 + t3;
t186 = -qJD(5) * t60 - t23;
t185 = t115 * t116 + t52;
t183 = qJD(1) + t219;
t182 = -t59 + (t196 + t220) * pkin(5);
t108 = t259 * t144;
t181 = pkin(9) * t196 + qJD(6) * t108 + t144 * t197 + t246;
t146 = sin(qJ(1));
t150 = cos(qJ(1));
t177 = g(1) * t146 - g(2) * t150;
t109 = t259 * t148;
t82 = t148 * t98;
t176 = qJD(6) * t109 + t148 * t197 - t144 * t58 + t82 + (pkin(5) * t149 + t202) * qJD(1);
t17 = qJD(4) * t166 - t262 * t71;
t72 = t96 * t149;
t175 = -t113 * t17 + t72 * t83;
t172 = t145 * t156 - t41 * t221;
t170 = t116 * t218 + t241;
t167 = t113 * t96;
t165 = qJDD(3) + t263;
t162 = -t148 * t217 + t192;
t161 = -g(1) * t126 - g(2) * t125 + t227;
t160 = -pkin(8) * t90 + t116 * t51;
t159 = t85 + t265;
t89 = t210 + t227;
t155 = -t115 * t151 + t161 + t210 + t89;
t124 = -pkin(5) * t148 - pkin(4);
t102 = -t145 * t151 + t207;
t101 = -t149 * t151 - t208;
t87 = t116 * t192;
t78 = t148 * t91;
t73 = (pkin(5) * t144 - t115) * t149;
t70 = t148 * t84;
t66 = -t125 * t144 + t126 * t229;
t65 = t125 * t148 + t126 * t230;
t64 = t125 * t229 + t126 * t144;
t63 = -t125 * t230 + t126 * t148;
t48 = -pkin(5) * t162 + t115 * t222;
t35 = -t144 * t256 + t245;
t26 = -pkin(9) * t228 + t189 * t145 + t70;
t11 = pkin(9) * t162 - t219 * t233 + t203;
t10 = t78 + (-t92 + (-t84 + t256) * t144) * qJD(5) + (t189 * t149 + t202) * qJD(4);
t9 = pkin(5) * t38 + t24;
t4 = t147 * t12 - t143 * t16;
t1 = [qJDD(1), t177, g(1) * t150 + g(2) * t146 (t177 + (t141 ^ 2 + t142 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) + 0.2e1 * t263 - t178, t161 + 0.2e1 * t210 + t227, t89 * t117 + t103 * qJD(3) + t165 * t119 - g(1) * (-pkin(1) * t146 - pkin(2) * t125 + qJ(3) * t126) - g(2) * (pkin(1) * t150 + pkin(2) * t126 + qJ(3) * t125) qJDD(1) * t138 - 0.2e1 * t145 * t190, -0.2e1 * t145 * t205 + 0.2e1 * t226 * t212, t102, t101, 0, t155 * t145 + t260 * t149, -t260 * t145 + t155 * t149, -t95 * t194 + (-t95 * t222 + t236) * t148 (t144 * t95 + t148 * t93) * t222 + (-t242 - t148 * t38 + (t144 * t93 - t237) * qJD(5)) * t149, t249 + t261, -t240 + t87 + (-t170 - t234) * t149, t116 * t221 + t145 * t90, -g(1) * t66 - g(2) * t64 + t78 * t116 + t70 * t90 + (t115 * t234 - t185 * t218 + t28) * t145 + (t21 * qJD(4) - t115 * t38 + t51 * t218) * t149 + ((-qJD(5) * t84 - t115 * t221) * t116 + t149 * t24 + (-qJD(4) * t51 - t115 * t90 + t186) * t145) * t144, -t203 * t116 - t245 * t90 + g(1) * t65 - g(2) * t63 + (t185 * t220 + (t115 * t95 - t148 * t51) * qJD(4) + t204) * t145 + (-t22 * qJD(4) - t115 * t37 + t24 * t148 - t51 * t220) * t149, -t17 * t174 - t7 * t72, -t156 * t72 - t17 * t41 + t174 * t18 - t7 * t71, -t175 + t258, t172 + t250, t113 * t221 + t145 * t83 (t10 * t147 - t11 * t143) * t113 + (-t143 * t35 + t147 * t26) * t83 + t198 * t145 + t4 * t221 + t48 * t41 - t73 * t156 + t9 * t71 + t29 * t18 - g(1) * t57 - g(2) * t55 + ((-t143 * t26 - t147 * t35) * t113 - t5 * t145) * qJD(6), -t5 * t221 + g(1) * t56 - g(2) * t54 + t13 * t145 + t29 * t17 - t48 * t174 + t73 * t7 - t9 * t72 + (-(-qJD(6) * t35 + t10) * t113 - t26 * t83 - t2 * t145) * t143 + (-(qJD(6) * t26 + t11) * t113 - t35 * t83 - t188 * t145) * t147; 0, 0, 0, t139, 0, 0, t139, 0, 0, 0, 0, 0, t101, -t102, 0, 0, 0, 0, 0, t240 + t87 + (-t170 + t234) * t149, t249 - t261, 0, 0, 0, 0, 0, -t172 + t250, t175 + t258; 0, 0, 0, 0, qJDD(1), -t152, t165 + t265, 0, 0, 0, 0, 0, t225 * t145 + t207, t225 * t149 - t208, 0, 0, 0, 0, 0, -t149 * t38 + (t234 - t241) * t145 + (-t144 * t221 - t183 * t148) * t116, -t236 + (qJD(4) * t95 - t238) * t145 + (t183 * t144 - t195) * t116, 0, 0, 0, 0, 0, qJD(1) * t167 + (-qJD(4) * t113 * t97 + t156) * t149 + ((t143 * t220 + t264 * t147 + t193) * t113 - t251 + qJD(4) * t41) * t145, t113 * t164 + (qJD(4) * t167 - t7) * t149 + (-(t264 * t143 - t144 * t215 - t147 * t220) * t113 + t252 - qJD(4) * t174) * t145; 0, 0, 0, 0, 0, 0, 0, t149 * t152 * t145, -t226 * t152, t205, -t206, qJDD(4), qJD(4) * t59 + t149 * t159 + t199 + t254, t266 + (-qJD(4) * t88 - t139) * t149 + (-t159 + t211) * t145, t116 * t237 + t242 (t37 - t244) * t148 + (-t38 - t243) * t144 (t116 * t229 - t149 * t95) * qJD(1) + t170, -t116 * t220 + t238 + (-t116 * t230 + t149 * t93) * qJD(1), -t116 * t223, -t21 * t223 - pkin(4) * t38 - t82 * t116 - t59 * t93 + (t58 * t116 + t160) * t144 - t274 * t148, -pkin(4) * t37 + t246 * t116 + t274 * t144 + t160 * t148 + t22 * t223 - t59 * t95, -t174 * t248 + t7 * t97, t156 * t97 + t174 * t247 - t248 * t41 - t7 * t96, t248 * t113 + t174 * t223 + t251, -t247 * t113 + t41 * t223 - t252, -t113 * t223 (-t108 * t147 - t109 * t143) * t83 - t124 * t156 + t9 * t96 - t4 * t223 + t182 * t41 + t247 * t29 + (t143 * t181 - t147 * t176) * t113 - t158 * t134 -(-t108 * t143 + t109 * t147) * t83 + t124 * t7 + t9 * t97 + t5 * t223 - t182 * t174 + t248 * t29 + (t143 * t176 + t147 * t181) * t113 + t158 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t93, -t93 ^ 2 + t95 ^ 2, t37 + t244, t243 - t38, t90, -t52 * t218 - g(1) * t63 - g(2) * t65 + t116 * t22 - t51 * t95 + t28 + (t186 + t253) * t144, g(1) * t64 - g(2) * t66 + g(3) * t228 + t116 * t21 + t51 * t93 - t168, -t273, t272, t271, t268, t83 -(-t143 * t15 - t239) * t113 + (-t113 * t216 + t147 * t83 - t95 * t41) * pkin(5) + t269 (-t16 * t113 - t2) * t143 + (t15 * t113 - t188) * t147 + (-t113 * t215 - t143 * t83 + t174 * t95) * pkin(5) + t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, t272, t271, t268, t83, t113 * t5 + t269, t113 * t4 - t143 * t2 - t147 * t188 + t270;];
tau_reg  = t1;
