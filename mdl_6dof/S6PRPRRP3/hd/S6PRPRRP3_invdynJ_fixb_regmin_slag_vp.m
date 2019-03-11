% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:48
% EndTime: 2019-03-08 20:07:57
% DurationCPUTime: 3.29s
% Computational Cost: add. (3851->386), mult. (9068->520), div. (0->0), fcn. (7430->14), ass. (0->208)
t150 = cos(pkin(11));
t157 = cos(qJ(4));
t148 = sin(pkin(11));
t154 = sin(qJ(4));
t239 = t148 * t154;
t112 = -t157 * t150 + t239;
t264 = pkin(8) + qJ(3);
t118 = t264 * t148;
t119 = t264 * t150;
t188 = -t118 * t157 - t119 * t154;
t41 = -t112 * qJD(3) + t188 * qJD(4);
t149 = sin(pkin(6));
t158 = cos(qJ(2));
t237 = t149 * t158;
t171 = t112 * t237;
t84 = qJD(1) * t171;
t289 = t41 + t84;
t109 = t112 * qJD(4);
t113 = t148 * t157 + t150 * t154;
t110 = t113 * qJD(4);
t155 = sin(qJ(2));
t231 = qJD(1) * t149;
t215 = t155 * t231;
t288 = pkin(4) * t110 + pkin(9) * t109 - t215;
t115 = qJD(2) * qJ(3) + t215;
t250 = cos(pkin(6));
t204 = qJD(1) * t250;
t130 = t150 * t204;
t76 = t130 + (-pkin(8) * qJD(2) - t115) * t148;
t230 = qJD(2) * t150;
t92 = t150 * t115 + t148 * t204;
t77 = pkin(8) * t230 + t92;
t33 = t154 * t76 + t157 * t77;
t287 = qJD(4) * t33;
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t286 = -t153 * t84 + t288 * t156;
t228 = qJD(5) * t153;
t222 = t150 * qJDD(2);
t223 = t148 * qJDD(2);
t195 = t154 * t223 - t157 * t222;
t67 = qJD(2) * t110 + t195;
t61 = qJDD(5) + t67;
t134 = t157 * t230;
t213 = qJD(2) * t239;
t106 = -t134 + t213;
t97 = qJD(5) + t106;
t285 = t156 * t61 - t97 * t228;
t227 = qJD(5) * t156;
t137 = pkin(3) * t150 + pkin(2);
t64 = pkin(4) * t112 - pkin(9) * t113 - t137;
t284 = t288 * t153 + t289 * t156 + t64 * t227;
t283 = -t154 * t77 + t157 * t76;
t147 = pkin(11) + qJ(4);
t140 = sin(t147);
t249 = cos(pkin(10));
t198 = t250 * t249;
t248 = sin(pkin(10));
t102 = t248 * t155 - t158 * t198;
t197 = t250 * t248;
t104 = t249 * t155 + t158 * t197;
t201 = g(1) * t104 + g(2) * t102;
t281 = g(3) * t237 - t201;
t282 = t281 * t140;
t159 = qJD(2) ^ 2;
t173 = (qJDD(2) * t158 - t155 * t159) * t149;
t103 = t155 * t198 + t248 * t158;
t105 = -t155 * t197 + t249 * t158;
t200 = g(1) * t105 + g(2) * t103;
t108 = t113 * qJD(2);
t214 = t158 * t231;
t196 = qJD(3) - t214;
t229 = qJD(2) * t155;
t210 = qJD(1) * t229;
t280 = t149 * t210 + qJDD(3);
t216 = qJD(4) * t134 + t154 * t222 + t157 * t223;
t169 = qJD(4) * t213 - t216;
t162 = -t156 * qJDD(4) - t153 * t169;
t88 = qJD(4) * t153 + t108 * t156;
t27 = t88 * qJD(5) + t162;
t279 = pkin(5) * t27 + qJDD(6);
t29 = qJD(4) * pkin(9) + t33;
t96 = -t137 * qJD(2) + t196;
t36 = pkin(4) * t106 - pkin(9) * t108 + t96;
t16 = t153 * t36 + t156 * t29;
t225 = qJDD(1) * t149;
t209 = t158 * t225;
t178 = -t209 + t280;
t79 = -t137 * qJDD(2) + t178;
t21 = t67 * pkin(4) + t169 * pkin(9) + t79;
t19 = t156 * t21;
t203 = qJDD(1) * t250;
t128 = t150 * t203;
t224 = qJDD(2) * qJ(3);
t90 = t155 * t225 + t224 + (qJD(3) + t214) * qJD(2);
t43 = t128 + (-pkin(8) * qJDD(2) - t90) * t148;
t59 = t148 * t203 + t150 * t90;
t44 = pkin(8) * t222 + t59;
t192 = t154 * t43 + t157 * t44;
t8 = qJDD(4) * pkin(9) + qJD(4) * t283 + t192;
t164 = -t16 * qJD(5) - t153 * t8 + t19;
t226 = t156 * qJD(4);
t26 = -qJD(5) * t226 - t153 * qJDD(4) + t108 * t228 + t156 * t169;
t1 = pkin(5) * t61 + qJ(6) * t26 - qJD(6) * t88 + t164;
t86 = t108 * t153 - t226;
t11 = -qJ(6) * t86 + t16;
t278 = t11 * t97 + t1;
t246 = qJDD(2) * pkin(2);
t93 = t178 - t246;
t184 = t201 - t93;
t277 = (-g(3) * t158 + t210) * t149 + t184 + t246;
t194 = t148 * (-t115 * t148 + t130) - t150 * t92;
t276 = t194 * t158 - (-qJD(2) * pkin(2) + t196) * t155;
t275 = t88 ^ 2;
t15 = -t153 * t29 + t156 * t36;
t10 = -qJ(6) * t88 + t15;
t6 = pkin(5) * t97 + t10;
t273 = t10 - t6;
t187 = qJ(6) * t109 - qJD(6) * t113;
t82 = -t118 * t154 + t119 * t157;
t75 = t156 * t82;
t272 = pkin(5) * t110 - t153 * t41 + t187 * t156 + (-t75 + (qJ(6) * t113 - t64) * t153) * qJD(5) + t286;
t211 = t113 * t227;
t271 = -qJ(6) * t211 + (-qJD(5) * t82 + t187) * t153 + t284;
t270 = pkin(5) * t153;
t265 = g(3) * t149;
t151 = -qJ(6) - pkin(9);
t263 = -t153 * t27 - t86 * t227;
t62 = pkin(4) * t108 + pkin(9) * t106;
t262 = t153 * t62 + t156 * t283;
t172 = t113 * t237;
t261 = -qJD(1) * t172 + t113 * qJD(3) + t82 * qJD(4);
t260 = t153 * t64 + t75;
t207 = qJD(5) * t151;
t236 = t153 * t106;
t259 = -qJ(6) * t236 + qJD(6) * t156 + t153 * t207 - t262;
t247 = qJ(6) * t156;
t53 = t156 * t62;
t258 = -pkin(5) * t108 - t106 * t247 + t156 * t207 - t53 + (-qJD(6) + t283) * t153;
t257 = t108 * t86;
t256 = t108 * t88;
t255 = t153 * t26;
t254 = t153 * t61;
t253 = t153 * t88;
t252 = t156 * t86;
t245 = qJDD(4) * pkin(4);
t244 = t109 * t153;
t243 = t109 * t156;
t242 = t113 * t153;
t141 = cos(t147);
t240 = t141 * t153;
t238 = t149 * t155;
t235 = t153 * t158;
t234 = t158 * t159;
t233 = qJDD(1) - g(3);
t232 = t148 ^ 2 + t150 ^ 2;
t220 = g(3) * t238;
t218 = t149 * t235;
t217 = t156 * t237;
t212 = t149 * t229;
t206 = t149 * t249;
t205 = t149 * t248;
t202 = t156 * t97;
t180 = -t153 * t21 - t156 * t8 - t36 * t227 + t29 * t228;
t2 = -qJ(6) * t27 - qJD(6) * t86 - t180;
t199 = -t97 * t6 + t2;
t193 = -t252 - t253;
t190 = t154 * t44 - t157 * t43 + t287;
t100 = -t148 * t238 + t250 * t150;
t101 = t250 * t148 + t150 * t238;
t189 = t100 * t157 - t101 * t154;
t48 = t100 * t154 + t101 * t157;
t186 = -t97 * t236 + t285;
t28 = -qJD(4) * pkin(4) - t283;
t182 = -t156 * t26 - t88 * t228;
t37 = -t153 * t48 - t217;
t181 = -t156 * t48 + t218;
t176 = -pkin(9) * t61 + t97 * t28;
t68 = t103 * t140 + t141 * t206;
t70 = t105 * t140 - t141 * t205;
t94 = t140 * t238 - t250 * t141;
t175 = g(1) * t70 + g(2) * t68 + g(3) * t94;
t69 = t103 * t141 - t140 * t206;
t71 = t105 * t141 + t140 * t205;
t95 = t250 * t140 + t141 * t238;
t174 = g(1) * t71 + g(2) * t69 + g(3) * t95;
t9 = t190 - t245;
t167 = t281 * t141;
t166 = -t281 + t209;
t58 = -t148 * t90 + t128;
t165 = -t148 * t58 + t150 * t59 - t200;
t163 = qJD(5) * pkin(9) * t97 - t175 + t9;
t161 = t175 - t190;
t160 = -g(1) * (t104 * t156 - t153 * t71) - g(2) * (t102 * t156 - t153 * t69) - g(3) * (-t153 * t95 - t217);
t139 = pkin(5) * t156 + pkin(4);
t121 = t151 * t156;
t120 = t151 * t153;
t85 = t86 ^ 2;
t56 = t156 * t64;
t31 = qJD(2) * t172 + t48 * qJD(4);
t30 = -qJD(2) * t171 + t189 * qJD(4);
t23 = -qJ(6) * t242 + t260;
t22 = pkin(5) * t86 + qJD(6) + t28;
t20 = pkin(5) * t112 - t113 * t247 - t153 * t82 + t56;
t14 = t181 * qJD(5) - t153 * t30 + t156 * t212;
t13 = t37 * qJD(5) + t153 * t212 + t156 * t30;
t3 = t9 + t279;
t4 = [t233, 0, t173 (-qJDD(2) * t155 - t234) * t149, t150 * t173, -t148 * t173, t232 * t149 * t234 + (-t100 * t148 + t101 * t150) * qJDD(2), t100 * t58 + t101 * t59 - g(3) + (-qJD(2) * t276 - t158 * t93) * t149, 0, 0, 0, 0, 0, -qJD(4) * t31 + qJDD(4) * t189 + (t106 * t229 - t158 * t67) * t149, -t30 * qJD(4) - t48 * qJDD(4) + (t108 * t229 + t158 * t169) * t149, 0, 0, 0, 0, 0, t14 * t97 - t189 * t27 + t31 * t86 + t37 * t61, -t13 * t97 + t181 * t61 + t189 * t26 + t31 * t88, -t13 * t86 - t14 * t88 + t181 * t27 + t26 * t37, t1 * t37 + t11 * t13 + t14 * t6 - t181 * t2 - t189 * t3 + t22 * t31 - g(3); 0, qJDD(2), t166, -t233 * t238 + t200, t277 * t150, -t277 * t148, -t220 + t165 + (qJD(2) * t196 + t224) * t232, -t194 * qJD(3) + t184 * pkin(2) + t165 * qJ(3) + (-g(3) * (pkin(2) * t158 + qJ(3) * t155) + t276 * qJD(1)) * t149, -t108 * t109 - t169 * t113, t109 * t106 - t108 * t110 + t169 * t112 - t113 * t67, -qJD(4) * t109 + qJDD(4) * t113, -qJD(4) * t110 - qJDD(4) * t112, 0, -t261 * qJD(4) + qJDD(4) * t188 - t106 * t215 + t110 * t96 + t112 * t79 - t137 * t67 - t167, -t82 * qJDD(4) - t137 * t216 + t79 * t113 - t96 * t109 - t108 * t215 + (t137 * t213 - t289) * qJD(4) + t282, t182 * t113 - t88 * t243, -t193 * t109 + (t255 - t156 * t27 + (t153 * t86 - t156 * t88) * qJD(5)) * t113, t110 * t88 - t112 * t26 + t113 * t285 - t97 * t243, t97 * t244 - t110 * t86 - t112 * t27 + (-t227 * t97 - t254) * t113, t110 * t97 + t112 * t61, t15 * t110 + t19 * t112 - t188 * t27 + t56 * t61 + t286 * t97 + t261 * t86 + (-t167 + (-t112 * t29 + t113 * t28 - t82 * t97) * qJD(5)) * t156 + ((-qJD(5) * t64 - t41) * t97 - t82 * t61 + (-qJD(5) * t36 - t8) * t112 + t9 * t113 - t28 * t109 - t220 - t200) * t153, -t260 * t61 + t180 * t112 - t16 * t110 + t188 * t26 - t28 * t243 - g(1) * (t104 * t240 + t105 * t156) - g(2) * (t102 * t240 + t103 * t156) + (t228 * t82 - t284) * t97 + t261 * t88 - (-t141 * t235 + t155 * t156) * t265 + (t9 * t156 - t228 * t28) * t113, t20 * t26 - t23 * t27 - t272 * t88 - t271 * t86 - (-t11 * t153 - t156 * t6) * t109 - t282 + (-t1 * t156 - t153 * t2 + (-t11 * t156 + t153 * t6) * qJD(5)) * t113, t2 * t23 + t1 * t20 + t3 * (pkin(5) * t242 - t188) + t272 * t6 + ((t211 - t244) * pkin(5) + t261) * t22 + t271 * t11 + (-t155 * t265 - t200) * (t264 + t270) + (-t158 * t265 + t201) * (t139 * t141 - t140 * t151 + t137); 0, 0, 0, 0, -t222, t223, -t232 * t159, t194 * qJD(2) - t166 - t246 + t280, 0, 0, 0, 0, 0, 0.2e1 * qJD(4) * t108 + t195 (-t106 - t213) * qJD(4) + t216, 0, 0, 0, 0, 0, t186 - t257, -t156 * t97 ^ 2 - t254 - t256 (-t252 + t253) * t106 - t182 + t263, -t108 * t22 + t199 * t153 + t156 * t278 + t281; 0, 0, 0, 0, 0, 0, 0, 0, t108 * t106, -t106 ^ 2 + t108 ^ 2 (t106 - t213) * qJD(4) + t216, -t195, qJDD(4), -t108 * t96 + t161 + t287, t106 * t96 + t174 - t192, t202 * t88 - t255, t106 * t193 + t182 + t263, t202 * t97 + t254 - t256, t186 + t257, -t97 * t108, -pkin(4) * t27 - t15 * t108 - t33 * t86 - t53 * t97 + (t283 * t97 + t176) * t153 - t163 * t156, pkin(4) * t26 + t16 * t108 + t163 * t153 + t176 * t156 + t262 * t97 - t33 * t88, t120 * t26 + t121 * t27 - t153 * t278 + t199 * t156 - t258 * t88 - t259 * t86 - t174, -t2 * t121 + t1 * t120 - t3 * t139 - g(1) * (-t139 * t70 - t151 * t71) - g(2) * (-t139 * t68 - t151 * t69) - g(3) * (-t139 * t94 - t151 * t95) + t258 * t6 + (t97 * t270 - t33) * t22 + t259 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t86, -t85 + t275, t86 * t97 - t26, -t162 + (-qJD(5) + t97) * t88, t61, t16 * t97 - t28 * t88 + t160 + t164, t15 * t97 + t28 * t86 - g(1) * (-t104 * t153 - t156 * t71) - g(2) * (-t102 * t153 - t156 * t69) - g(3) * (-t156 * t95 + t218) + t180, pkin(5) * t26 + t273 * t86, -t273 * t11 + (-t22 * t88 + t1 + t160) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 - t275, t11 * t86 + t6 * t88 - t161 - t245 + t279;];
tau_reg  = t4;
