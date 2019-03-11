% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:49
% EndTime: 2019-03-08 20:48:59
% DurationCPUTime: 3.82s
% Computational Cost: add. (2829->436), mult. (5984->623), div. (0->0), fcn. (4657->14), ass. (0->229)
t150 = sin(qJ(5));
t154 = cos(qJ(5));
t236 = t154 * qJD(4);
t155 = cos(qJ(4));
t251 = qJD(2) * t155;
t104 = t150 * t251 - t236;
t248 = qJD(4) * t150;
t106 = t154 * t251 + t248;
t149 = sin(qJ(6));
t153 = cos(qJ(6));
t185 = t104 * t149 - t106 * t153;
t42 = t104 * t153 + t106 * t149;
t312 = t185 * t42;
t151 = sin(qJ(4));
t228 = t155 * qJDD(2);
t235 = qJD(2) * qJD(4);
t311 = -t151 * t235 + t228;
t310 = t185 ^ 2 - t42 ^ 2;
t252 = qJD(2) * t151;
t133 = qJD(5) + t252;
t127 = qJD(6) + t133;
t238 = qJD(6) * t153;
t239 = qJD(6) * t149;
t216 = t151 * t236;
t240 = qJD(5) * t155;
t172 = -t150 * t240 - t216;
t37 = qJD(2) * t172 + qJD(5) * t236 + qJDD(4) * t150 + t154 * t228;
t38 = qJD(5) * t106 - t154 * qJDD(4) + t150 * t311;
t6 = -t104 * t238 - t106 * t239 - t149 * t38 + t153 * t37;
t309 = t127 * t42 + t6;
t148 = cos(pkin(6));
t256 = qJD(1) * t155;
t130 = t148 * t256;
t157 = -pkin(2) - pkin(8);
t156 = cos(qJ(2));
t147 = sin(pkin(6));
t257 = qJD(1) * t147;
t212 = t156 * t257;
t190 = qJD(3) - t212;
t93 = qJD(2) * t157 + t190;
t64 = t151 * t93 + t130;
t50 = qJD(4) * pkin(9) + t64;
t113 = pkin(4) * t151 - pkin(9) * t155 + qJ(3);
t152 = sin(qJ(2));
t213 = t152 * t257;
t71 = qJD(2) * t113 + t213;
t21 = t150 * t71 + t154 * t50;
t12 = -pkin(10) * t104 + t21;
t10 = t12 * t239;
t145 = qJ(5) + qJ(6);
t141 = sin(t145);
t142 = cos(t145);
t271 = t147 * t152;
t269 = t148 * t151;
t303 = -qJD(1) * t269 + t155 * t93;
t49 = -qJD(4) * pkin(4) - t303;
t33 = pkin(5) * t104 + t49;
t146 = sin(pkin(11));
t274 = t146 * t147;
t281 = cos(pkin(11));
t201 = t281 * t152;
t272 = t146 * t156;
t87 = t148 * t272 + t201;
t55 = t151 * t87 + t155 * t274;
t202 = t147 * t281;
t200 = t281 * t156;
t273 = t146 * t152;
t85 = -t148 * t200 + t273;
t57 = -t151 * t85 + t155 * t202;
t86 = t148 * t201 + t272;
t88 = -t148 * t273 + t200;
t270 = t147 * t156;
t92 = t148 * t155 - t151 * t270;
t308 = t33 * t42 - g(1) * (-t141 * t88 - t142 * t55) - g(2) * (-t141 * t86 + t142 * t57) - g(3) * (-t141 * t271 - t142 * t92) + t10;
t208 = t155 * t235;
t229 = t151 * qJDD(2);
t100 = qJDD(5) + t208 + t229;
t233 = qJDD(1) * t148;
t205 = t155 * t233;
t123 = qJD(2) * t213;
t234 = qJDD(1) * t147;
t206 = t156 * t234;
t181 = qJDD(3) + t123 - t206;
t66 = qJDD(2) * t157 + t181;
t17 = qJDD(4) * pkin(9) + qJD(4) * t303 + t151 * t66 + t205;
t193 = pkin(4) * t155 + pkin(9) * t151;
t101 = qJD(4) * t193 + qJD(3);
t207 = t152 * t234;
t32 = t207 + t113 * qJDD(2) + (t101 + t212) * qJD(2);
t31 = t154 * t32;
t165 = -qJD(5) * t21 - t150 * t17 + t31;
t2 = t100 * pkin(5) - t37 * pkin(10) + t165;
t241 = qJD(5) * t154;
t227 = -t150 * t32 - t154 * t17 - t241 * t71;
t243 = qJD(5) * t150;
t180 = t243 * t50 + t227;
t3 = -pkin(10) * t38 - t180;
t223 = -t149 * t3 + t153 * t2;
t285 = t12 * t153;
t20 = -t150 * t50 + t154 * t71;
t11 = -pkin(10) * t106 + t20;
t9 = pkin(5) * t133 + t11;
t5 = t149 * t9 + t285;
t307 = t33 * t185 - g(1) * (-t141 * t55 + t142 * t88) - g(2) * (t141 * t57 + t142 * t86) - g(3) * (-t141 * t92 + t142 * t271) - t5 * qJD(6) + t223;
t163 = qJD(6) * t185 - t149 * t37 - t153 * t38;
t306 = -t127 * t185 + t163;
t266 = t151 * t152;
t305 = -(-t150 * t266 + t154 * t156) * t257 + t154 * t101;
t215 = t155 * t236;
t264 = t152 * t154;
t304 = (t150 * t156 + t151 * t264) * t257 - t150 * t101 - t113 * t241 - t157 * t215;
t108 = t149 * t154 + t150 * t153;
t78 = t108 * t155;
t168 = g(1) * t87 + g(2) * t85 - g(3) * t270;
t302 = -t149 * t243 - t150 * t239;
t301 = -qJD(6) * t154 - t241;
t300 = qJD(5) + qJD(6);
t210 = qJD(6) * t9 + t3;
t299 = t149 * t2 + t153 * t210;
t194 = g(1) * t88 + g(2) * t86;
t260 = qJDD(1) - g(3);
t298 = -t260 * t271 + t194;
t255 = qJD(2) * qJ(3);
t112 = t213 + t255;
t297 = qJD(4) * (-t112 + t213 - t255) - qJDD(4) * t157;
t169 = g(3) * t271 + t194;
t296 = (t133 * t157 + t50) * qJD(5) + t169;
t295 = pkin(9) + pkin(10);
t294 = pkin(10) * t155;
t265 = t151 * t154;
t126 = t157 * t265;
t267 = t150 * t157;
t204 = pkin(5) - t267;
t226 = pkin(10) * t265;
t292 = (-t126 + (-t113 + t294) * t150) * qJD(5) + (t155 * t204 + t226) * qJD(4) + t305;
t247 = qJD(4) * t151;
t217 = t150 * t247;
t171 = -t154 * t240 + t217;
t242 = qJD(5) * t151;
t291 = -pkin(10) * t171 + t242 * t267 + t304;
t109 = t193 * qJD(2);
t290 = t109 * t150 + t154 * t303;
t107 = t149 * t150 - t153 * t154;
t175 = t107 * t151;
t289 = -qJD(2) * t175 - t107 * t300;
t174 = t108 * qJD(2);
t288 = t108 * t300 + t151 * t174;
t97 = qJDD(6) + t100;
t287 = t107 * t97;
t286 = t108 * t97;
t283 = t37 * t150;
t282 = t113 * t150 + t126;
t280 = qJDD(2) * pkin(2);
t279 = t104 * t133;
t278 = t106 * t133;
t277 = t106 * t154;
t276 = t141 * t151;
t275 = t142 * t151;
t268 = t150 * t100;
t263 = t154 * t100;
t262 = t154 * t155;
t159 = qJD(2) ^ 2;
t261 = t156 * t159;
t144 = t155 ^ 2;
t259 = t151 ^ 2 - t144;
t158 = qJD(4) ^ 2;
t258 = -t158 - t159;
t254 = qJD(2) * t112;
t253 = qJD(2) * t147;
t250 = qJD(4) * t104;
t249 = qJD(4) * t106;
t246 = qJD(4) * t155;
t245 = qJD(4) * t157;
t244 = qJD(5) * t133;
t232 = qJDD(2) * qJ(3);
t231 = qJDD(4) * t151;
t224 = -qJD(4) * t130 - t151 * t233 - t247 * t93;
t222 = qJD(5) * t295;
t221 = t152 * t256;
t220 = t152 * t253;
t219 = t156 * t253;
t218 = t150 * t252;
t198 = -t66 + t254;
t197 = qJD(2) + t242;
t196 = -t64 + (t218 + t243) * pkin(5);
t121 = t295 * t150;
t192 = pkin(10) * t218 + qJD(6) * t121 + t150 * t222 + t290;
t122 = t295 * t154;
t95 = t154 * t109;
t191 = qJD(6) * t122 + t154 * t222 - t150 * t303 + t95 + (pkin(5) * t155 + t226) * qJD(2);
t99 = t154 * t113;
t40 = -pkin(10) * t262 + t151 * t204 + t99;
t53 = -t150 * t294 + t282;
t189 = t149 * t40 + t153 * t53;
t60 = t147 * t264 - t150 * t92;
t61 = t150 * t271 + t154 * t92;
t188 = -t149 * t61 + t153 * t60;
t187 = t149 * t60 + t153 * t61;
t186 = (-qJD(2) * pkin(2) + t190) * t152 + t112 * t156;
t184 = qJDD(2) * t152 + t261;
t91 = t155 * t270 + t269;
t179 = t133 * t241 + t268;
t178 = -t133 * t243 + t263;
t177 = g(1) * (-t151 * t274 + t155 * t87) + g(2) * (t151 * t202 + t155 * t85) - g(3) * t91;
t176 = t127 * t107;
t18 = -qJDD(4) * pkin(4) - t155 * t66 - t224;
t167 = -pkin(9) * t100 + t133 * t49;
t166 = t206 + t168;
t164 = qJDD(3) - t166;
t162 = pkin(9) * t244 + t177 + t18;
t67 = t207 + t232 + (qJD(3) + t212) * qJD(2);
t160 = qJD(2) * t190 - t157 * t158 - t169 + t232 + t67;
t139 = qJDD(4) * t155;
t136 = -pkin(5) * t154 - pkin(4);
t102 = (pkin(5) * t150 - t157) * t155;
t84 = t184 * t147;
t83 = (-qJDD(2) * t156 + t152 * t159) * t147;
t79 = t107 * t155;
t72 = t181 - t280;
t68 = -pkin(5) * t171 + t151 * t245;
t59 = qJD(4) * t92 - t155 * t220;
t58 = -qJD(4) * t91 + t151 * t220;
t23 = -t149 * t216 + t302 * t155 + (t262 * t300 - t217) * t153;
t22 = qJD(4) * t175 - t300 * t78;
t14 = qJD(5) * t60 + t150 * t219 + t58 * t154;
t13 = -qJD(5) * t61 - t58 * t150 + t154 * t219;
t8 = pkin(5) * t38 + t18;
t4 = -t12 * t149 + t153 * t9;
t1 = [t260, 0, -t83, -t84, t83, t84, t148 ^ 2 * qJDD(1) - g(3) + (qJD(2) * t186 + t152 * t67 - t156 * t72) * t147, 0, 0, 0, 0, 0, -t59 * qJD(4) - t91 * qJDD(4) + (t151 * t184 + t152 * t208) * t147, -t58 * qJD(4) - t92 * qJDD(4) + (t152 * t311 + t155 * t261) * t147, 0, 0, 0, 0, 0, t100 * t60 + t104 * t59 + t13 * t133 + t38 * t91, -t100 * t61 + t106 * t59 - t133 * t14 + t37 * t91, 0, 0, 0, 0, 0 (-qJD(6) * t187 + t13 * t153 - t14 * t149) * t127 + t188 * t97 + t59 * t42 - t91 * t163 -(qJD(6) * t188 + t149 * t13 + t153 * t14) * t127 - t187 * t97 - t59 * t185 + t91 * t6; 0, qJDD(2), t166, t298, t164 - 0.2e1 * t280, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t232 - t298, t67 * qJ(3) + t112 * qJD(3) - t72 * pkin(2) - g(1) * (-pkin(2) * t87 + qJ(3) * t88) - g(2) * (-pkin(2) * t85 + qJ(3) * t86) + (-g(3) * (pkin(2) * t156 + qJ(3) * t152) - t186 * qJD(1)) * t147, qJDD(2) * t144 - 0.2e1 * t151 * t208, -0.2e1 * t151 * t228 + 0.2e1 * t235 * t259, -t151 * t158 + t139, -t155 * t158 - t231, 0, t151 * t160 - t155 * t297, t151 * t297 + t155 * t160, t106 * t172 + t262 * t37 (t104 * t154 + t106 * t150) * t247 + (-t283 - t154 * t38 + (t104 * t150 - t277) * qJD(5)) * t155 (-t133 * t236 + t37) * t151 + (t178 + t249) * t155 (t133 * t248 - t38) * t151 + (-t179 - t250) * t155, t100 * t151 + t133 * t246, t99 * t100 + t305 * t133 + (-t113 * t244 + t168) * t150 + (t104 * t245 + t31 + (-qJD(4) * t49 - qJD(5) * t71 - t100 * t157 - t17) * t150 - t296 * t154) * t151 + (t104 * t213 + t49 * t241 + t18 * t150 - t157 * t38 + (-t133 * t267 + t20) * qJD(4)) * t155, -t282 * t100 + t304 * t133 + t168 * t154 + ((t106 * t157 - t154 * t49) * qJD(4) + t296 * t150 + t227) * t151 + (-qJD(4) * t21 + t106 * t213 + t18 * t154 - t157 * t37 - t243 * t49) * t155, -t185 * t22 - t6 * t79, -t163 * t79 + t185 * t23 - t22 * t42 - t6 * t78, t127 * t22 + t151 * t6 - t185 * t246 - t79 * t97, -t127 * t23 + t151 * t163 - t246 * t42 - t78 * t97, t127 * t246 + t151 * t97 (-t149 * t53 + t153 * t40) * t97 + t223 * t151 + t4 * t246 + t68 * t42 - t102 * t163 + t8 * t78 + t33 * t23 - g(1) * (-t141 * t87 + t275 * t88) - g(2) * (-t141 * t85 + t275 * t86) + (t149 * t291 + t153 * t292) * t127 + (-t127 * t189 - t151 * t5) * qJD(6) + (t42 * t221 - g(3) * (t141 * t156 + t142 * t266)) * t147, -t189 * t97 - (-t10 + t299) * t151 - t5 * t246 - t68 * t185 + t102 * t6 - t8 * t79 + t33 * t22 - g(1) * (-t142 * t87 - t276 * t88) - g(2) * (-t142 * t85 - t276 * t86) + ((-qJD(6) * t40 + t291) * t153 + (qJD(6) * t53 - t292) * t149) * t127 + (-t185 * t221 - g(3) * (-t141 * t266 + t142 * t156)) * t147; 0, 0, 0, 0, qJDD(2), -t159, t123 + t164 - t254 - t280, 0, 0, 0, 0, 0, t151 * t258 + t139, t155 * t258 - t231, 0, 0, 0, 0, 0, -t155 * t38 + (t250 - t268) * t151 + (-t150 * t246 - t154 * t197) * t133, -t155 * t37 + (t249 - t263) * t151 + (t150 * t197 - t215) * t133, 0, 0, 0, 0, 0, qJD(2) * t176 + (-qJD(4) * t108 * t127 + t163) * t155 + ((t153 * t301 - t302) * t127 - t286 + qJD(4) * t42) * t151, t127 * t174 + (qJD(4) * t176 - t6) * t155 + (-(t149 * t301 - t150 * t238 - t153 * t243) * t127 + t287 - qJD(4) * t185) * t151; 0, 0, 0, 0, 0, 0, 0, t155 * t159 * t151, -t259 * t159, t228, -t229, qJDD(4), t64 * qJD(4) - t155 * t198 - t177 + t224, g(1) * t55 - g(2) * t57 + g(3) * t92 + t151 * t198 - t205, t133 * t277 + t283 (t37 - t279) * t154 + (-t38 - t278) * t150 (-t106 * t155 + t133 * t265) * qJD(2) + t179 (-t133 * t150 * t151 + t104 * t155) * qJD(2) + t178, -t133 * t251, -t20 * t251 - pkin(4) * t38 - t64 * t104 - t95 * t133 + (t133 * t303 + t167) * t150 - t162 * t154, -pkin(4) * t37 - t64 * t106 + t133 * t290 + t150 * t162 + t154 * t167 + t21 * t251, t108 * t6 - t185 * t289, -t107 * t6 + t108 * t163 + t185 * t288 - t289 * t42, t127 * t289 + t185 * t251 + t286, -t127 * t288 + t251 * t42 - t287, -t127 * t251 (-t121 * t153 - t122 * t149) * t97 - t136 * t163 + t8 * t107 - t4 * t251 + t196 * t42 + t288 * t33 + (t149 * t192 - t153 * t191) * t127 - t177 * t142 -(-t121 * t149 + t122 * t153) * t97 + t136 * t6 + t8 * t108 + t5 * t251 - t196 * t185 + t289 * t33 + (t149 * t191 + t153 * t192) * t127 + t177 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t104, -t104 ^ 2 + t106 ^ 2, t37 + t279, t278 - t38, t100, t21 * t133 - t49 * t106 - g(1) * (-t150 * t55 + t154 * t88) - g(2) * (t150 * t57 + t154 * t86) - g(3) * t60 + t165, t20 * t133 + t49 * t104 - g(1) * (-t150 * t88 - t154 * t55) - g(2) * (-t150 * t86 + t154 * t57) + g(3) * t61 + t180, -t312, t310, t309, t306, t97 -(-t11 * t149 - t285) * t127 + (-t106 * t42 - t127 * t239 + t153 * t97) * pkin(5) + t307 (-t12 * t127 - t2) * t149 + (t11 * t127 - t210) * t153 + (t106 * t185 - t127 * t238 - t149 * t97) * pkin(5) + t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t312, t310, t309, t306, t97, t5 * t127 + t307, t4 * t127 - t299 + t308;];
tau_reg  = t1;
