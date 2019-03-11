% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:24:55
% EndTime: 2019-03-09 05:25:07
% DurationCPUTime: 5.20s
% Computational Cost: add. (5751->475), mult. (11769->644), div. (0->0), fcn. (8231->12), ass. (0->244)
t213 = sin(qJ(3));
t298 = qJD(1) * t213;
t190 = qJD(4) + t298;
t185 = qJD(6) + t190;
t212 = sin(qJ(4));
t216 = cos(qJ(4));
t285 = t216 * qJD(3);
t217 = cos(qJ(3));
t297 = qJD(1) * t217;
t168 = t212 * t297 - t285;
t294 = qJD(3) * t212;
t170 = t216 * t297 + t294;
t208 = sin(pkin(10));
t209 = cos(pkin(10));
t104 = t168 * t208 - t170 * t209;
t211 = sin(qJ(6));
t215 = cos(qJ(6));
t246 = -t168 * t209 - t170 * t208;
t306 = t215 * t246;
t52 = t104 * t211 + t306;
t328 = t185 * t52;
t286 = qJD(6) * t211;
t271 = t213 * t285;
t288 = qJD(4) * t217;
t231 = -t212 * t288 - t271;
t280 = t217 * qJDD(1);
t97 = qJD(1) * t231 + qJD(4) * t285 + t212 * qJDD(3) + t216 * t280;
t293 = qJD(3) * t213;
t272 = t212 * t293;
t98 = -qJD(1) * t272 + qJD(4) * t170 - t216 * qJDD(3) + t212 * t280;
t43 = -t208 * t97 - t209 * t98;
t44 = -t208 * t98 + t209 * t97;
t8 = qJD(6) * t306 + t104 * t286 + t211 * t43 + t215 * t44;
t359 = t8 - t328;
t346 = -t215 * t104 + t211 * t246;
t358 = t346 * t52;
t357 = t346 ^ 2 - t52 ^ 2;
t200 = qJ(4) + pkin(10) + qJ(6);
t192 = sin(t200);
t193 = cos(t200);
t218 = cos(qJ(1));
t214 = sin(qJ(1));
t310 = t213 * t214;
t121 = t192 * t218 + t193 * t310;
t308 = t213 * t218;
t123 = -t192 * t214 + t193 * t308;
t284 = qJD(1) * qJD(3);
t267 = t217 * t284;
t281 = t213 * qJDD(1);
t161 = qJDD(4) + t267 + t281;
t254 = pkin(3) * t217 + pkin(8) * t213;
t162 = qJD(3) * t254 + qJD(2);
t174 = pkin(3) * t213 - pkin(8) * t217 + qJ(2);
t110 = qJD(1) * t162 + qJDD(1) * t174;
t100 = t216 * t110;
t219 = -pkin(1) - pkin(7);
t183 = t219 * qJDD(1) + qJDD(2);
t187 = t219 * qJD(1) + qJD(2);
t292 = qJD(3) * t217;
t115 = qJDD(3) * pkin(8) + t183 * t213 + t187 * t292;
t147 = t174 * qJD(1);
t173 = t213 * t187;
t154 = qJD(3) * pkin(8) + t173;
t94 = t147 * t212 + t154 * t216;
t14 = pkin(4) * t161 - qJ(5) * t97 - qJD(4) * t94 - qJD(5) * t170 - t212 * t115 + t100;
t289 = qJD(4) * t216;
t277 = -t212 * t110 - t216 * t115 - t147 * t289;
t290 = qJD(4) * t212;
t235 = -t154 * t290 - t277;
t19 = -qJ(5) * t98 - qJD(5) * t168 + t235;
t4 = t209 * t14 - t19 * t208;
t2 = pkin(5) * t161 - pkin(9) * t44 + t4;
t69 = -qJ(5) * t168 + t94;
t327 = t209 * t69;
t93 = t216 * t147 - t154 * t212;
t68 = -qJ(5) * t170 + t93;
t61 = pkin(4) * t190 + t68;
t28 = t208 * t61 + t327;
t345 = pkin(9) * t246;
t21 = t28 + t345;
t20 = t21 * t286;
t334 = g(3) * t217;
t315 = t187 * t217;
t155 = -qJD(3) * pkin(3) - t315;
t112 = pkin(4) * t168 + qJD(5) + t155;
t58 = -pkin(5) * t246 + t112;
t356 = g(1) * t121 - g(2) * t123 + t193 * t334 - t211 * t2 - t58 * t52 + t20;
t329 = t185 * t346;
t9 = qJD(6) * t346 + t211 * t44 - t215 * t43;
t354 = -t9 + t329;
t210 = -qJ(5) - pkin(8);
t262 = qJD(4) * t210;
t273 = t212 * t298;
t287 = qJD(5) * t216;
t172 = t254 * qJD(1);
t304 = t216 * t217;
t302 = t212 * t172 + t187 * t304;
t353 = -qJ(5) * t273 + t212 * t262 + t287 - t302;
t153 = t216 * t172;
t309 = t213 * t216;
t311 = t212 * t217;
t352 = -qJD(5) * t212 + t216 * t262 + t187 * t311 - t153 - (pkin(4) * t217 + qJ(5) * t309) * qJD(1);
t261 = -qJDD(3) * pkin(3) + t187 * t293;
t316 = t183 * t217;
t114 = t261 - t316;
t335 = g(3) * t213;
t342 = -g(1) * t214 + g(2) * t218;
t226 = -t217 * t342 - t335;
t351 = qJD(4) * pkin(8) * t190 + t114 + t226;
t120 = -t192 * t310 + t193 * t218;
t122 = t192 * t308 + t193 * t214;
t5 = t208 * t14 + t209 * t19;
t3 = pkin(9) * t43 + t5;
t274 = t215 * t2 - t211 * t3;
t350 = -g(1) * t120 - g(2) * t122 + t192 * t334 - t58 * t346 + t274;
t349 = pkin(9) * t104;
t164 = t208 * t216 + t209 * t212;
t144 = t164 * qJD(4);
t145 = t164 * qJD(1);
t348 = t213 * t145 + t144;
t313 = t209 * t216;
t243 = t208 * t212 - t313;
t340 = qJD(4) * t243;
t347 = -t208 * t273 + t298 * t313 - t340;
t331 = -t208 * t353 + t209 * t352;
t330 = t208 * t352 + t209 * t353;
t301 = t212 * t174 + t219 * t309;
t303 = t216 * t218;
t148 = -t212 * t310 + t303;
t307 = t214 * t216;
t150 = t212 * t308 + t307;
t341 = -g(1) * t148 - g(2) * t150;
t339 = pkin(4) * t208;
t142 = t216 * t162;
t266 = -t212 * t219 + pkin(4);
t42 = qJ(5) * t271 + t142 - t301 * qJD(4) + (qJ(5) * t290 + qJD(3) * t266 - t287) * t217;
t269 = t216 * t288;
t291 = qJD(3) * t219;
t270 = t217 * t291;
t276 = t212 * t162 + t174 * t289 + t216 * t270;
t46 = -qJ(5) * t269 + (-qJD(5) * t217 + (qJ(5) * qJD(3) - qJD(4) * t219) * t213) * t212 + t276;
t16 = t208 * t42 + t209 * t46;
t247 = -t164 * t211 - t215 * t243;
t333 = qJD(6) * t247 - t211 * t348 + t215 * t347;
t103 = t164 * t215 - t211 * t243;
t332 = qJD(6) * t103 + t211 * t347 + t215 * t348;
t64 = t208 * t69;
t32 = t209 * t68 - t64;
t326 = t212 * t97;
t27 = t209 * t61 - t64;
t18 = pkin(5) * t190 + t27 + t349;
t325 = t215 * t18;
t160 = t216 * t174;
t101 = -qJ(5) * t304 + t213 * t266 + t160;
t113 = -qJ(5) * t311 + t301;
t55 = t208 * t101 + t209 * t113;
t324 = -t173 + t348 * pkin(5) + (t273 + t290) * pkin(4);
t134 = t164 * t217;
t323 = -t243 * qJD(1) + qJD(3) * t134 - t213 * t340;
t136 = t243 * t217;
t322 = -qJD(3) * t136 - t144 * t213 - t145;
t321 = pkin(1) * qJDD(1);
t221 = qJD(1) ^ 2;
t320 = qJ(2) * t221;
t319 = t168 * t190;
t318 = t170 * t190;
t317 = t170 * t216;
t314 = t190 * t212;
t312 = t212 * t161;
t305 = t216 * t161;
t180 = t210 * t212;
t181 = t210 * t216;
t117 = t208 * t180 - t209 * t181;
t207 = t217 ^ 2;
t300 = t213 ^ 2 - t207;
t220 = qJD(3) ^ 2;
t299 = -t220 - t221;
t296 = qJD(3) * t168;
t295 = qJD(3) * t170;
t283 = qJDD(1) * qJ(2);
t282 = qJDD(3) * t213;
t278 = 0.2e1 * qJD(1) * qJD(2);
t196 = pkin(4) * t216 + pkin(3);
t275 = pkin(4) * t212 + pkin(7);
t268 = g(2) * (t218 * pkin(1) + t214 * qJ(2));
t264 = qJD(6) * t18 + t3;
t15 = -t208 * t46 + t209 * t42;
t31 = -t208 * t68 - t327;
t54 = t209 * t101 - t113 * t208;
t260 = t190 * t219 + t154;
t116 = t209 * t180 + t181 * t208;
t259 = pkin(4) * t311 - t217 * t219;
t258 = -qJD(4) * t147 - t115;
t257 = qJD(4) * t213 + qJD(1);
t88 = -pkin(9) * t243 + t117;
t256 = pkin(5) * t297 + pkin(9) * t347 + qJD(6) * t88 - t331;
t87 = -pkin(9) * t164 + t116;
t255 = pkin(9) * t348 - qJD(6) * t87 - t330;
t253 = g(1) * t218 + g(2) * t214;
t251 = qJDD(2) + t342;
t133 = t164 * t213;
t250 = qJD(6) * t133 - t322;
t135 = t243 * t213;
t249 = -qJD(6) * t135 + t323;
t7 = t211 * t18 + t215 * t21;
t248 = -t215 * t134 + t136 * t211;
t77 = -t134 * t211 - t136 * t215;
t244 = t196 * t213 + t210 * t217;
t242 = -t183 - t342;
t194 = pkin(4) * t209 + pkin(5);
t241 = t194 * t211 + t215 * t339;
t240 = t194 * t215 - t211 * t339;
t238 = t190 * t289 + t312;
t237 = -t190 * t290 + t305;
t236 = t213 * t291 + (t269 - t272) * pkin(4);
t234 = 0.2e1 * qJ(2) * t284 + qJDD(3) * t219;
t232 = pkin(4) * t98 + qJDD(5) + t261;
t230 = t242 + t320;
t229 = -pkin(8) * t161 + t155 * t190;
t56 = t232 - t316;
t224 = -t253 + t278 + 0.2e1 * t283;
t222 = -t219 * t220 + t224;
t203 = t218 * qJ(2);
t199 = qJDD(3) * t217;
t157 = qJDD(6) + t161;
t151 = -t212 * t214 + t213 * t303;
t149 = t212 * t218 + t213 * t307;
t127 = pkin(5) * t243 - t196;
t109 = pkin(5) * t134 + t259;
t84 = t164 * t288 - t208 * t272 + t209 * t271;
t82 = t164 * t293 + t217 * t340;
t70 = pkin(4) * t170 - pkin(5) * t104;
t57 = -pkin(5) * t82 + t236;
t34 = -pkin(9) * t134 + t55;
t33 = pkin(5) * t213 + pkin(9) * t136 + t54;
t26 = qJD(6) * t77 - t211 * t84 - t215 * t82;
t25 = qJD(6) * t248 + t211 * t82 - t215 * t84;
t24 = t32 + t349;
t23 = t31 - t345;
t22 = -pkin(5) * t43 + t56;
t11 = pkin(9) * t82 + t16;
t10 = pkin(5) * t292 + pkin(9) * t84 + t15;
t6 = -t21 * t211 + t325;
t1 = [qJDD(1), -t342, t253, t251 - 0.2e1 * t321, t224 -(qJDD(2) - t321) * pkin(1) - g(1) * (-pkin(1) * t214 + t203) - t268 + (t278 + t283) * qJ(2), qJDD(1) * t207 - 0.2e1 * t213 * t267, -0.2e1 * t213 * t280 + 0.2e1 * t284 * t300, -t213 * t220 + t199, -t217 * t220 - t282, 0, t213 * t222 + t217 * t234, -t213 * t234 + t217 * t222, t170 * t231 + t304 * t97 (t168 * t216 + t170 * t212) * t293 + (-t326 - t216 * t98 + (t168 * t212 - t317) * qJD(4)) * t217 (-t190 * t285 + t97) * t213 + (t237 + t295) * t217 (t190 * t294 - t98) * t213 + (-t238 - t296) * t217, t213 * t161 + t190 * t292, -g(1) * t151 - g(2) * t149 + t142 * t190 + t160 * t161 + (t168 * t291 - t260 * t289 + t100) * t213 + (t93 * qJD(3) + t155 * t289 - t219 * t98) * t217 + ((-qJD(4) * t174 - t270) * t190 + t114 * t217 + (-t155 * qJD(3) - t219 * t161 + t258) * t213) * t212, -t276 * t190 - t301 * t161 + g(1) * t150 - g(2) * t148 + (t260 * t290 + (-t155 * t216 + t170 * t219) * qJD(3) + t277) * t213 + (-qJD(3) * t94 + t114 * t216 - t155 * t290 - t219 * t97) * t217, t104 * t15 - t134 * t5 + t136 * t4 + t16 * t246 + t217 * t253 + t27 * t84 + t28 * t82 + t43 * t55 - t44 * t54, t5 * t55 + t28 * t16 + t4 * t54 + t27 * t15 + t56 * t259 + t112 * t236 - g(1) * t203 - t268 + (-g(1) * t244 - g(2) * t275) * t218 + (-g(1) * (-pkin(1) - t275) - g(2) * t244) * t214, t25 * t346 + t77 * t8, t248 * t8 + t25 * t52 - t26 * t346 - t77 * t9, t157 * t77 + t185 * t25 + t213 * t8 + t292 * t346, t157 * t248 - t185 * t26 - t213 * t9 + t292 * t52, t157 * t213 + t185 * t292 (t10 * t215 - t11 * t211) * t185 + (-t211 * t34 + t215 * t33) * t157 + t274 * t213 + t6 * t292 - t57 * t52 + t109 * t9 - t22 * t248 + t58 * t26 - g(1) * t123 - g(2) * t121 + ((-t211 * t33 - t215 * t34) * t185 - t7 * t213) * qJD(6), -t7 * t292 + g(1) * t122 - g(2) * t120 + t109 * t8 + t20 * t213 + t22 * t77 + t58 * t25 + t57 * t346 + (-(-qJD(6) * t34 + t10) * t185 - t33 * t157 - t2 * t213) * t211 + (-(qJD(6) * t33 + t11) * t185 - t34 * t157 - t264 * t213) * t215; 0, 0, 0, qJDD(1), -t221, t251 - t320 - t321, 0, 0, 0, 0, 0, t213 * t299 + t199, t217 * t299 - t282, 0, 0, 0, 0, 0, -t217 * t98 + (t296 - t312) * t213 + (-t212 * t292 - t216 * t257) * t190, -t217 * t97 + (t295 - t305) * t213 + (t212 * t257 - t217 * t285) * t190, -t104 * t323 + t133 * t44 - t135 * t43 + t246 * t322, t112 * t293 - t133 * t4 - t135 * t5 - t217 * t56 - t27 * t323 + t28 * t322 + t342, 0, 0, 0, 0, 0 (-t133 * t215 + t135 * t211) * t157 - t52 * t293 - t217 * t9 + (t211 * t250 - t215 * t249) * t185 -(-t133 * t211 - t135 * t215) * t157 + t346 * t293 - t217 * t8 + (t211 * t249 + t215 * t250) * t185; 0, 0, 0, 0, 0, 0, t217 * t221 * t213, -t300 * t221, t280, -t281, qJDD(3), -t217 * t230 + t335, t213 * t230 + t334, t190 * t317 + t326 (t97 - t319) * t216 + (-t98 - t318) * t212 (-t170 * t217 + t190 * t309) * qJD(1) + t238 (t168 * t217 - t213 * t314) * qJD(1) + t237, -t190 * t297, -t93 * t297 - t168 * t173 - pkin(3) * t98 - t153 * t190 + (t190 * t315 + t229) * t212 - t351 * t216, -pkin(3) * t97 - t170 * t173 + t302 * t190 + t212 * t351 + t229 * t216 + t94 * t297, t331 * t104 - t116 * t44 + t117 * t43 - t164 * t4 + t213 * t342 - t243 * t5 + t330 * t246 - t27 * t347 - t28 * t348 - t334, t5 * t117 + t4 * t116 - t56 * t196 + g(3) * t244 + t330 * t28 + t331 * t27 + (pkin(4) * t314 - t173) * t112 + t342 * (t196 * t217 - t210 * t213) t103 * t8 + t333 * t346, -t103 * t9 + t247 * t8 - t332 * t346 + t333 * t52, t103 * t157 + t185 * t333 - t297 * t346, t157 * t247 - t185 * t332 - t297 * t52, -t185 * t297 (-t211 * t88 + t215 * t87) * t157 + t127 * t9 - t22 * t247 - t6 * t297 + t332 * t58 - t324 * t52 + (t211 * t255 - t215 * t256) * t185 - t226 * t193 -(t211 * t87 + t215 * t88) * t157 + t127 * t8 + t22 * t103 + t7 * t297 + t333 * t58 + t324 * t346 + (t211 * t256 + t215 * t255) * t185 + t226 * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170 * t168, -t168 ^ 2 + t170 ^ 2, t97 + t319, t318 - t98, t161, -t154 * t289 - t155 * t170 + t190 * t94 + t100 + (t258 + t334) * t212 + t341, g(1) * t149 - g(2) * t151 + g(3) * t304 + t155 * t168 + t190 * t93 - t235 (t208 * t43 - t209 * t44) * pkin(4) + (-t32 + t27) * t246 + (-t28 - t31) * t104, -t27 * t31 - t28 * t32 + (g(3) * t311 - t112 * t170 + t5 * t208 + t4 * t209 + t341) * pkin(4), -t358, t357, t359, t354, t157, t240 * t157 - (-t211 * t24 + t215 * t23) * t185 + t70 * t52 + (-t185 * t241 - t7) * qJD(6) + t350, -t241 * t157 - t215 * t3 + (t211 * t23 + t215 * t24) * t185 - t70 * t346 + (-t185 * t240 - t325) * qJD(6) + t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 ^ 2 - t246 ^ 2, -t104 * t27 + t217 * t242 - t246 * t28 + t232 - t335, 0, 0, 0, 0, 0, t9 + t329, t8 + t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t358, t357, t359, t354, t157 (-qJD(6) + t185) * t7 + t350, t185 * t6 - t215 * t264 + t356;];
tau_reg  = t1;
