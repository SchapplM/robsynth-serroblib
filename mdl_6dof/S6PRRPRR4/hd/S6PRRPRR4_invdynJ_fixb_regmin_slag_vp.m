% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:31
% EndTime: 2019-03-08 22:14:43
% DurationCPUTime: 5.11s
% Computational Cost: add. (3191->417), mult. (7225->560), div. (0->0), fcn. (5782->12), ass. (0->235)
t171 = cos(pkin(6));
t176 = sin(qJ(2));
t291 = sin(pkin(11));
t239 = t291 * t176;
t170 = cos(pkin(11));
t180 = cos(qJ(2));
t282 = t170 * t180;
t91 = t171 * t282 - t239;
t238 = t291 * t180;
t283 = t170 * t176;
t93 = -t171 * t238 - t283;
t224 = g(1) * t93 + g(2) * t91;
t169 = sin(pkin(6));
t284 = t169 * t180;
t258 = g(3) * t284;
t196 = -t224 - t258;
t175 = sin(qJ(3));
t156 = t175 * qJD(4);
t179 = cos(qJ(3));
t259 = t179 * qJDD(2);
t155 = t175 * qJDD(2);
t261 = qJD(2) * qJD(3);
t243 = t179 * t261;
t334 = -t243 - t155;
t262 = qJD(1) * qJD(2);
t247 = t176 * t262;
t122 = t169 * t247;
t131 = qJDD(1) * t284;
t168 = qJDD(2) * pkin(2);
t83 = t122 - t131 - t168;
t201 = pkin(3) * t259 - t334 * qJ(4) + qJD(2) * t156 - t83;
t244 = t175 * t261;
t37 = pkin(3) * t244 - t201;
t340 = -t37 + t196;
t161 = qJDD(3) - qJDD(5);
t174 = sin(qJ(5));
t178 = cos(qJ(5));
t181 = -pkin(3) - pkin(4);
t286 = t169 * t176;
t251 = qJD(1) * t286;
t113 = qJD(2) * pkin(8) + t251;
t268 = qJD(3) * t175;
t245 = qJD(1) * t268;
t260 = qJDD(1) * t171;
t267 = qJD(3) * t179;
t246 = t180 * t262;
t289 = qJDD(2) * pkin(8);
t84 = t289 + (qJDD(1) * t176 + t246) * t169;
t228 = t113 * t267 + t171 * t245 + t175 * t84 - t179 * t260;
t217 = -qJDD(4) - t228;
t18 = t334 * pkin(9) + t181 * qJDD(3) - t217;
t163 = qJDD(3) * qJ(4);
t164 = qJD(3) * qJD(4);
t273 = qJD(1) * t171;
t135 = t179 * t273;
t257 = qJD(3) * t135 + t175 * t260 + t179 * t84;
t29 = -t113 * t268 + t163 + t164 + t257;
t333 = t244 - t259;
t19 = t333 * pkin(9) + t29;
t265 = qJD(5) * t178;
t266 = qJD(5) * t174;
t252 = t181 * qJD(3);
t271 = qJD(2) * t175;
t101 = t175 * t113;
t79 = -t101 + t135;
t335 = qJD(4) - t79;
t336 = -pkin(9) * t271 + t335;
t53 = t252 + t336;
t165 = qJD(3) * qJ(4);
t269 = qJD(2) * t179;
t80 = t179 * t113 + t175 * t273;
t65 = -pkin(9) * t269 + t80;
t55 = t165 + t65;
t216 = t174 * t19 - t178 * t18 + t55 * t265 + t53 * t266;
t2 = pkin(5) * t161 + t216;
t96 = -t171 * t179 + t175 * t286;
t285 = t169 * t179;
t255 = t176 * t285;
t97 = t171 * t175 + t255;
t218 = t174 * t97 - t96 * t178;
t92 = t171 * t283 + t238;
t58 = t170 * t285 + t175 * t92;
t59 = -t170 * t169 * t175 + t179 * t92;
t240 = t169 * t291;
t94 = -t171 * t239 + t282;
t60 = t175 * t94 - t179 * t240;
t61 = t175 * t240 + t94 * t179;
t207 = g(1) * (t174 * t61 - t178 * t60) + g(2) * (t174 * t59 - t178 * t58) + g(3) * t218;
t339 = t2 - t207;
t104 = t174 * t175 + t178 * t179;
t323 = t104 * qJD(2);
t332 = qJD(6) + t323;
t248 = t174 * t269;
t100 = t178 * t271 - t248;
t162 = qJD(3) - qJD(5);
t173 = sin(qJ(6));
t177 = cos(qJ(6));
t71 = t100 * t173 + t177 * t162;
t338 = t332 * t71;
t325 = t177 * t332;
t134 = qJD(1) * t284;
t172 = qJD(2) * pkin(2);
t114 = -t134 - t172;
t81 = -pkin(3) * t269 - qJ(4) * t271 + t114;
t66 = pkin(4) * t269 - t81;
t337 = -t66 * t100 + t207 - t216;
t324 = t332 - qJD(6);
t198 = t174 * t267 + t175 * t265;
t43 = qJD(2) * t198 - qJD(5) * t248 + qJDD(2) * t104 - t178 * t244;
t38 = qJDD(6) + t43;
t297 = t173 * t38;
t73 = t100 * t177 - t162 * t173;
t331 = -t73 * t100 + t325 * t332 + t297;
t263 = qJD(6) * t177;
t264 = qJD(6) * t173;
t202 = t333 * t174 - t334 * t178;
t204 = t104 * qJD(5);
t42 = -qJD(2) * t204 + t202;
t9 = -t100 * t264 - t173 * t161 - t162 * t263 + t177 * t42;
t304 = t9 * t173;
t330 = t325 * t73 + t304;
t235 = t177 * t161 + t173 * t42;
t10 = qJD(6) * t73 + t235;
t329 = (-t9 + t338) * t177 + (t332 * t73 + t10) * t173;
t20 = -t174 * t55 + t178 * t53;
t16 = pkin(5) * t162 - t20;
t328 = t16 * t332;
t70 = t165 + t80;
t67 = -qJD(3) * pkin(3) + t335;
t327 = t100 * t162 + t43;
t226 = t175 * t252;
t277 = qJ(4) * t267 + t156;
t326 = t226 + t277 + t251;
t276 = t178 * qJ(4) + t174 * t181;
t211 = t179 * pkin(3) + t175 * qJ(4) + pkin(2);
t54 = t100 * pkin(5) + pkin(10) * t323;
t320 = (pkin(10) * qJD(6) + t54) * t332 + t339;
t108 = -pkin(10) + t276;
t148 = qJ(4) * t269;
t88 = t181 * t271 + t148;
t319 = (qJD(6) * t108 - t54 + t88) * t332 - t339;
t208 = t174 * t18 + t178 * t19 + t53 * t265 - t266 * t55;
t1 = -pkin(10) * t161 + t208;
t213 = t174 * t179 - t175 * t178;
t223 = g(1) * t94 + g(2) * t92;
t314 = pkin(8) - pkin(9);
t109 = t314 * t268;
t121 = t314 * t179;
t110 = qJD(3) * t121;
t120 = t314 * t175;
t214 = t120 * t178 - t121 * t174;
t78 = t104 * t284;
t303 = -qJD(1) * t78 + qJD(5) * t214 - t178 * t109 + t174 * t110;
t33 = pkin(5) * t323 - pkin(10) * t100 + t66;
t103 = t179 * pkin(4) + t211;
t46 = pkin(5) * t104 + pkin(10) * t213 + t103;
t56 = qJD(3) * t104 - t204;
t75 = t120 * t174 + t121 * t178;
t318 = -(qJD(6) * t33 + t1) * t104 + t16 * t56 - t2 * t213 + (-qJD(6) * t46 - t303) * t332 + g(3) * t286 - t75 * t38 + t223;
t295 = t177 * t38;
t296 = t173 * t332;
t317 = -t71 * t100 + t296 * t332 - t295;
t316 = -g(1) * t61 - g(2) * t59 - g(3) * t97 - (t79 + t101) * qJD(3) + t257;
t24 = t174 * t58 + t178 * t59;
t27 = t174 * t60 + t178 * t61;
t51 = t174 * t96 + t178 * t97;
t315 = g(1) * t27 + g(2) * t24 + g(3) * t51 + t323 * t66 - t208;
t182 = qJD(3) ^ 2;
t308 = pkin(8) * t182;
t307 = g(3) * t180;
t21 = t174 * t53 + t178 * t55;
t17 = -pkin(10) * t162 + t21;
t219 = t17 * t173 - t177 * t33;
t306 = t219 * t100;
t6 = t17 * t177 + t173 * t33;
t305 = t6 * t100;
t302 = qJD(5) * t75 - t174 * t109 - t178 * t110 - t213 * t134;
t215 = -qJ(4) * t174 + t178 * t181;
t301 = qJD(5) * t215 - t174 * t65 + t336 * t178;
t300 = t276 * qJD(5) + t336 * t174 + t178 * t65;
t299 = t100 * t323;
t298 = t213 * t16;
t294 = t177 * t73;
t293 = t332 * t100;
t292 = t323 * t162;
t290 = pkin(8) * qJDD(3);
t288 = qJDD(3) * pkin(3);
t280 = qJDD(1) - g(3);
t279 = t179 * t122 + t245 * t284;
t166 = t175 ^ 2;
t167 = t179 ^ 2;
t275 = t166 - t167;
t274 = t166 + t167;
t270 = qJD(2) * t176;
t256 = t332 * t271;
t183 = qJD(2) ^ 2;
t254 = t175 * t183 * t179;
t253 = t100 ^ 2 - t323 ^ 2;
t250 = t169 * t270;
t249 = qJD(2) * t284;
t237 = t114 - t172;
t231 = -qJD(2) * t211 + t81;
t229 = t332 * t162;
t227 = t162 ^ 2;
t225 = t175 * t249;
t57 = -t178 * t268 - t179 * t266 + t198;
t222 = t57 * pkin(5) - t56 * pkin(10) + t326;
t212 = qJDD(2) * t180 - t176 * t183;
t210 = -t263 * t332 - t297;
t209 = -t264 * t332 + t295;
t39 = -t173 * t51 + t177 * t284;
t40 = t173 * t284 + t177 * t51;
t205 = g(3) * t78 + t224 * t104;
t195 = -pkin(10) * t38 + t20 * t332 + t328;
t193 = t46 * t38 - t205;
t191 = -t168 + t224 + t83 + t308;
t190 = -t108 * t38 - t301 * t332 - t328;
t189 = g(1) * t60 + g(2) * t58 + g(3) * t96 - t228;
t188 = qJD(3) * t80 + t189;
t187 = qJDD(2) * t211 - t308 + t340;
t28 = pkin(4) * t259 + qJD(2) * t226 + t201;
t63 = qJD(3) * t97 + t225;
t186 = -t183 * t255 - t96 * qJDD(3) + t259 * t284 + (-t63 - t225) * qJD(3);
t32 = -t217 - t288;
t185 = t32 * t175 + t29 * t179 + (-t175 * t70 + t179 * t67) * qJD(3) - t223;
t107 = pkin(5) - t215;
t106 = pkin(3) * t271 - t148;
t89 = pkin(3) * t268 - t277;
t62 = -qJD(3) * t96 + t179 * t249;
t12 = t62 * qJD(3) + t97 * qJDD(3) + (t175 * t212 + t180 * t243) * t169;
t8 = -qJD(5) * t218 + t63 * t174 + t62 * t178;
t7 = qJD(5) * t51 + t62 * t174 - t63 * t178;
t4 = pkin(5) * t43 - pkin(10) * t42 + t28;
t3 = t177 * t4;
t5 = [t280, 0, t212 * t169 (-qJDD(2) * t176 - t180 * t183) * t169, 0, 0, 0, 0, 0, t186, -t12, t186 (t175 * t96 + t179 * t97) * qJDD(2) + (t175 * t63 + t179 * t62 + (-t175 * t97 + t179 * t96) * qJD(3)) * qJD(2), t12, t29 * t97 + t32 * t96 + t70 * t62 + t67 * t63 - g(3) + (-t180 * t37 + t270 * t81) * t169, 0, 0, 0, 0, 0, t218 * t161 + t7 * t162 + (t180 * t43 - t270 * t323) * t169, t51 * t161 + t8 * t162 + (-t100 * t270 + t180 * t42) * t169, 0, 0, 0, 0, 0 (-qJD(6) * t40 - t8 * t173 - t177 * t250) * t332 + t39 * t38 + t7 * t71 + t218 * t10 -(qJD(6) * t39 - t173 * t250 + t8 * t177) * t332 - t40 * t38 + t7 * t73 + t218 * t9; 0, qJDD(2), t131 + t196, -t280 * t286 + t223, qJDD(2) * t166 + 0.2e1 * t175 * t243, 0.2e1 * t175 * t259 - 0.2e1 * t261 * t275, qJDD(3) * t175 + t179 * t182, qJDD(3) * t179 - t175 * t182, 0 (qJD(3) * t237 - t290) * t175 + (-t191 - t258) * t179 + t279 (-t290 + (t237 + t134) * qJD(3)) * t179 + ((-t247 + t307) * t169 + t191) * t175 (qJD(3) * t231 - t290) * t175 + (-qJD(2) * t89 + t187) * t179 + t279, t274 * t289 + (-g(3) * t176 - t246 * t274) * t169 + t185 (t290 + (-t231 - t134) * qJD(3)) * t179 + ((-t89 + t251) * qJD(2) + t187) * t175, t81 * t89 + t185 * pkin(8) + ((-pkin(8) * g(3) - qJD(1) * t81) * t176 + (-t175 * t67 - t179 * t70) * qJD(1) * t180) * t169 + t340 * t211, t100 * t56 - t213 * t42, -t100 * t57 - t104 * t42 + t213 * t43 - t323 * t56, t161 * t213 - t162 * t56, t104 * t161 + t162 * t57, 0, t103 * t43 + t28 * t104 - t161 * t214 + t302 * t162 + t323 * t326 + t66 * t57 - t205, t326 * t100 + t103 * t42 + t75 * t161 + t303 * t162 + t66 * t56 + (t307 * t169 + t224 - t28) * t213, t56 * t294 - (t9 * t177 - t264 * t73) * t213 (-t173 * t73 - t177 * t71) * t56 - (-t10 * t177 - t304 + (t173 * t71 - t294) * qJD(6)) * t213, t9 * t104 - t209 * t213 + t325 * t56 + t73 * t57, -t10 * t104 - t210 * t213 - t296 * t56 - t71 * t57, t104 * t38 + t332 * t57, -t214 * t10 + t3 * t104 - t219 * t57 + t302 * t71 + (t222 * t332 + (-t104 * t17 - t332 * t75 - t298) * qJD(6) + t193) * t177 + t318 * t173, -t6 * t57 - t214 * t9 + t302 * t73 + (-(-qJD(6) * t17 + t4) * t104 + qJD(6) * t298 + (qJD(6) * t75 - t222) * t332 - t193) * t173 + t318 * t177; 0, 0, 0, 0, -t254, t275 * t183, t155, t259, qJDD(3), -t114 * t271 + t188, -t114 * t269 - t316, 0.2e1 * t288 - qJDD(4) + (t106 * t179 - t175 * t81) * qJD(2) + t188 (-pkin(3) * t175 + qJ(4) * t179) * qJDD(2), 0.2e1 * t163 + 0.2e1 * t164 + (t106 * t175 + t179 * t81) * qJD(2) + t316, t29 * qJ(4) - t32 * pkin(3) - t81 * t106 - t67 * t80 - g(1) * (-pkin(3) * t60 + qJ(4) * t61) - g(2) * (-pkin(3) * t58 + qJ(4) * t59) - g(3) * (-pkin(3) * t96 + qJ(4) * t97) + t335 * t70, -t299, -t253, qJD(5) * t323 - t202 + t292, t327, t161, -t161 * t215 + t162 * t300 - t323 * t88 - t337, -t88 * t100 + t161 * t276 + t162 * t301 - t315, -t330, t329, -t331, t317, t293, t107 * t10 + t190 * t173 - t319 * t177 + t300 * t71 - t306, t107 * t9 + t319 * t173 + t190 * t177 + t300 * t73 - t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t254, t155, -t166 * t183 - t182, -qJD(3) * t70 + t271 * t81 + qJDD(4) - t189 - t288, 0, 0, 0, 0, 0, -t178 * t161 - t174 * t227 - t271 * t323, -t100 * t271 + t174 * t161 - t178 * t227, 0, 0, 0, 0, 0, -t177 * t256 + (t173 * t229 - t10) * t178 + (-t162 * t71 + t210) * t174, t173 * t256 + (t177 * t229 - t9) * t178 + (-t162 * t73 - t209) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t253, t42 - t292, -t327, -t161, -t162 * t21 + t337, -t162 * t20 + t315, t330, -t329, t331, -t317, -t293, -pkin(5) * t10 + t195 * t173 - t320 * t177 - t21 * t71 + t306, -pkin(5) * t9 + t320 * t173 + t195 * t177 - t21 * t73 + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t9 + t338, t324 * t73 - t235, t38, -t173 * t1 + t3 - t16 * t73 - g(1) * (-t173 * t27 + t177 * t93) - g(2) * (-t173 * t24 + t177 * t91) - g(3) * t39 + t324 * t6, -t177 * t1 - t173 * t4 + t16 * t71 - g(1) * (-t173 * t93 - t177 * t27) - g(2) * (-t173 * t91 - t177 * t24) + g(3) * t40 - t324 * t219;];
tau_reg  = t5;
