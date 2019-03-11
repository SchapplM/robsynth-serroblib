% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:05
% EndTime: 2019-03-09 09:37:18
% DurationCPUTime: 5.85s
% Computational Cost: add. (5595->529), mult. (11988->684), div. (0->0), fcn. (8530->14), ass. (0->268)
t239 = sin(qJ(2));
t335 = qJD(1) * t239;
t200 = qJD(5) + t335;
t189 = qJD(6) + t200;
t237 = sin(qJ(6));
t241 = cos(qJ(6));
t235 = cos(pkin(10));
t243 = cos(qJ(2));
t334 = qJD(1) * t243;
t311 = t235 * t334;
t234 = sin(pkin(10));
t333 = qJD(2) * t234;
t165 = t311 + t333;
t309 = t234 * t334;
t167 = qJD(2) * t235 - t309;
t238 = sin(qJ(5));
t242 = cos(qJ(5));
t94 = t165 * t238 - t167 * t242;
t95 = t242 * t165 + t167 * t238;
t375 = t237 * t94 - t241 * t95;
t386 = t189 * t375;
t170 = t234 * t242 + t235 * t238;
t153 = t170 * qJD(5);
t263 = t170 * t239;
t340 = -qJD(1) * t263 - t153;
t370 = pkin(3) + pkin(7);
t385 = t200 * t95;
t279 = t237 * t95 + t241 * t94;
t384 = t279 * t375;
t383 = t94 * t200;
t382 = t189 * t279;
t312 = t235 * t335;
t329 = qJD(5) * t242;
t330 = qJD(5) * t238;
t348 = t234 * t238;
t339 = -t234 * t330 + t235 * t329 + t242 * t312 - t335 * t348;
t274 = -t235 * t242 + t348;
t211 = pkin(7) * t335;
t381 = qJD(3) + t211;
t380 = t279 ^ 2 - t375 ^ 2;
t327 = qJD(6) * t241;
t328 = qJD(6) * t237;
t321 = qJD(1) * qJD(2);
t308 = t239 * t321;
t319 = t243 * qJDD(1);
t116 = qJDD(2) * t234 + (-t308 + t319) * t235;
t293 = t235 * qJDD(2) - t234 * t319;
t117 = t234 * t308 + t293;
t33 = -t238 * t116 + t242 * t117 - t165 * t329 - t167 * t330;
t34 = -qJD(5) * t94 + t242 * t116 + t238 * t117;
t7 = -t237 * t34 + t241 * t33 - t95 * t327 + t328 * t94;
t379 = t7 - t386;
t220 = t239 * qJ(3);
t303 = -pkin(1) - t220;
t364 = pkin(2) + qJ(4);
t259 = -t243 * t364 + t303;
t127 = t259 * qJD(1);
t323 = pkin(3) * t335 + t381;
t137 = -qJD(2) * t364 + t323;
t64 = -t127 * t234 + t235 * t137;
t44 = pkin(4) * t335 - pkin(8) * t167 + t64;
t65 = t235 * t127 + t234 * t137;
t51 = -pkin(8) * t165 + t65;
t18 = t238 * t44 + t242 * t51;
t13 = -pkin(9) * t95 + t18;
t11 = t13 * t328;
t228 = pkin(10) + qJ(5);
t218 = qJ(6) + t228;
t206 = sin(t218);
t207 = cos(t218);
t240 = sin(qJ(1));
t244 = cos(qJ(1));
t344 = t239 * t244;
t119 = t206 * t344 + t207 * t240;
t345 = t239 * t240;
t121 = -t206 * t345 + t207 * t244;
t227 = g(3) * t243;
t212 = pkin(7) * t334;
t177 = pkin(3) * t334 + t212;
t231 = qJD(2) * qJ(3);
t148 = qJD(4) + t231 + t177;
t103 = pkin(4) * t165 + t148;
t45 = pkin(5) * t95 + t103;
t378 = g(1) * t119 - g(2) * t121 - t206 * t227 - t375 * t45 + t11;
t118 = -t206 * t240 + t207 * t344;
t120 = t206 * t244 + t207 * t345;
t307 = t243 * t321;
t320 = t239 * qJDD(1);
t261 = t307 + t320;
t169 = qJDD(5) + t261;
t199 = pkin(2) * t308;
t354 = qJ(3) * t243;
t277 = qJ(4) * t239 - t354;
t325 = t239 * qJD(3);
t248 = qJD(2) * t277 - t243 * qJD(4) - t325;
t57 = qJD(1) * t248 + qJDD(1) * t259 + t199;
t198 = pkin(7) * t307;
t208 = pkin(7) * t320;
t304 = qJDD(3) + t198 + t208;
t88 = pkin(3) * t261 - qJD(2) * qJD(4) - qJDD(2) * t364 + t304;
t28 = -t234 * t57 + t235 * t88;
t19 = pkin(4) * t261 - t117 * pkin(8) + t28;
t29 = t234 * t88 + t235 * t57;
t25 = -pkin(8) * t116 + t29;
t301 = t242 * t19 - t238 * t25;
t253 = -qJD(5) * t18 + t301;
t2 = t169 * pkin(5) - t33 * pkin(9) + t253;
t268 = t238 * t19 + t242 * t25 + t44 * t329 - t330 * t51;
t3 = -pkin(9) * t34 + t268;
t314 = t241 * t2 - t237 * t3;
t17 = -t238 * t51 + t242 * t44;
t12 = pkin(9) * t94 + t17;
t10 = pkin(5) * t200 + t12;
t358 = t13 * t241;
t5 = t10 * t237 + t358;
t377 = -g(1) * t118 - g(2) * t120 - qJD(6) * t5 + t207 * t227 + t45 * t279 + t314;
t252 = qJD(6) * t279 - t237 * t33 - t241 * t34;
t376 = t252 - t382;
t374 = t340 * t241;
t100 = -t170 * t237 - t241 * t274;
t363 = -pkin(8) - t364;
t179 = t363 * t234;
t180 = t363 * t235;
t338 = t242 * t179 + t238 * t180;
t272 = -pkin(8) * t234 * t239 + pkin(4) * t243;
t215 = pkin(2) * t335;
t144 = qJD(1) * t277 + t215;
t86 = -t234 * t144 + t235 * t177;
t60 = qJD(1) * t272 + t86;
t87 = t235 * t144 + t234 * t177;
t69 = pkin(8) * t312 + t87;
t373 = qJD(4) * t274 - t338 * qJD(5) + t238 * t69 - t242 * t60;
t286 = g(1) * t244 + g(2) * t240;
t315 = -pkin(4) * t235 - pkin(3);
t324 = -t315 * t335 + t381;
t331 = qJD(2) * t243;
t372 = t364 * t331 + (qJD(4) - t148) * t239;
t371 = qJD(4) * t170 + t179 * t330 - t180 * t329 + t238 * t60 + t242 * t69;
t369 = g(1) * t240;
t366 = g(2) * t244;
t365 = g(3) * t239;
t224 = t243 * pkin(2);
t99 = t241 * t170 - t237 * t274;
t362 = -qJD(6) * t99 - t237 * t339 + t374;
t361 = t100 * qJD(6) + t237 * t340 + t241 * t339;
t337 = t224 + t220;
t353 = qJ(4) * t243;
t284 = t337 + t353;
t164 = -pkin(1) - t284;
t186 = t370 * t239;
t173 = t235 * t186;
t78 = t239 * pkin(4) + t173 + (pkin(8) * t243 - t164) * t234;
t102 = t235 * t164 + t234 * t186;
t346 = t235 * t243;
t84 = -pkin(8) * t346 + t102;
t359 = t238 * t78 + t242 * t84;
t163 = qJDD(6) + t169;
t357 = t99 * t163;
t356 = pkin(5) * t339 + t324;
t355 = pkin(7) * qJDD(2);
t352 = qJDD(2) * pkin(2);
t351 = t100 * t163;
t232 = t239 ^ 2;
t246 = qJD(1) ^ 2;
t349 = t232 * t246;
t343 = t239 * t246;
t342 = t240 * t243;
t341 = t243 * t244;
t205 = t234 * pkin(4) + qJ(3);
t332 = qJD(2) * t239;
t214 = pkin(2) * t332;
t105 = t214 + t248;
t178 = t370 * t331;
t68 = t235 * t105 + t234 * t178;
t187 = t370 * t243;
t233 = t243 ^ 2;
t336 = t232 - t233;
t326 = t148 * qJD(2);
t318 = t243 * t343;
t209 = pkin(7) * t319;
t229 = qJDD(2) * qJ(3);
t230 = qJD(2) * qJD(3);
t317 = t209 + t229 + t230;
t152 = pkin(4) * t346 + t187;
t316 = -g(1) * t344 - g(2) * t345 + t227;
t313 = t370 * qJD(2);
t306 = t234 * t320;
t305 = t235 * t320;
t302 = qJD(6) * t10 + t3;
t67 = -t234 * t105 + t235 * t178;
t52 = qJD(2) * t272 + t67;
t56 = pkin(8) * t235 * t332 + t68;
t299 = -t238 * t56 + t242 * t52;
t297 = -t238 * t84 + t242 * t78;
t296 = -qJD(2) * pkin(2) + qJD(3);
t294 = -t179 * t238 + t242 * t180;
t292 = t244 * pkin(1) + pkin(2) * t341 + t240 * pkin(7) + qJ(3) * t344;
t291 = -t208 - t316;
t290 = qJD(1) * t313;
t73 = -pkin(9) * t170 + t338;
t289 = pkin(5) * t334 + pkin(9) * t340 + qJD(6) * t73 - t373;
t72 = pkin(9) * t274 + t294;
t288 = pkin(9) * t339 - qJD(6) * t72 + t371;
t245 = qJD(2) ^ 2;
t287 = pkin(7) * t245 + t366;
t285 = -t366 + t369;
t283 = -t274 * t169 + t200 * t340;
t281 = qJD(6) * t274 - t339;
t280 = t29 * t234 + t28 * t235;
t140 = t274 * t243;
t141 = t170 * t243;
t276 = t241 * t140 + t141 * t237;
t76 = t140 * t237 - t141 * t241;
t181 = t211 + t296;
t185 = -t212 - t231;
t275 = t181 * t243 + t185 * t239;
t273 = pkin(3) * t319 + qJDD(4) + t317;
t271 = t303 - t224;
t270 = t286 * t243;
t269 = -0.2e1 * pkin(1) * t321 - t355;
t267 = t238 * t52 + t242 * t56 + t78 * t329 - t330 * t84;
t159 = t271 * qJD(1);
t266 = t159 * t335 + qJDD(3) - t291;
t265 = (-t234 * t64 + t235 * t65) * t239;
t130 = (-pkin(7) + t315) * t332;
t264 = -qJ(3) * t331 - t325;
t262 = -t170 * t169 - t200 * t339;
t93 = -t239 * t290 + t273;
t260 = -t243 * t93 + t286;
t258 = 0.2e1 * qJDD(1) * pkin(1) - t287;
t182 = -pkin(1) - t337;
t257 = t355 + (-qJD(1) * t182 - t159) * qJD(2);
t255 = -t270 - t365;
t254 = t93 + t255;
t146 = t214 + t264;
t92 = qJD(1) * t264 + qJDD(1) * t271 + t199;
t251 = qJD(1) * t146 + qJDD(1) * t182 + t287 + t92;
t128 = pkin(7) * t308 - t317;
t143 = t304 - t352;
t249 = qJD(2) * t275 - t128 * t243 + t143 * t239;
t55 = t116 * pkin(4) + t93;
t225 = t244 * pkin(7);
t217 = cos(t228);
t216 = sin(t228);
t203 = g(1) * t342;
t197 = qJ(3) * t341;
t194 = qJ(3) * t342;
t176 = t239 * t313;
t174 = -qJ(3) * t334 + t215;
t135 = -t216 * t345 + t217 * t244;
t134 = t216 * t244 + t217 * t345;
t133 = t216 * t344 + t217 * t240;
t132 = -t216 * t240 + t217 * t344;
t124 = pkin(5) * t170 + t205;
t101 = -t164 * t234 + t173;
t91 = -pkin(5) * t140 + t152;
t81 = t243 * t153 - t274 * t332;
t80 = qJD(2) * t263 + qJD(5) * t140;
t50 = -t81 * pkin(5) + t130;
t27 = pkin(9) * t140 + t359;
t26 = pkin(5) * t239 + pkin(9) * t141 + t297;
t21 = qJD(6) * t76 + t237 * t80 - t241 * t81;
t20 = qJD(6) * t276 + t237 * t81 + t241 * t80;
t14 = t34 * pkin(5) + t55;
t9 = pkin(9) * t81 + t267;
t6 = pkin(5) * t331 - t80 * pkin(9) - qJD(5) * t359 + t299;
t4 = t10 * t241 - t13 * t237;
t1 = [qJDD(1), t285, t286, qJDD(1) * t232 + 0.2e1 * t239 * t307, 0.2e1 * t239 * t319 - 0.2e1 * t321 * t336, qJDD(2) * t239 + t243 * t245, qJDD(2) * t243 - t239 * t245, 0, t239 * t269 + t243 * t258 + t203, t269 * t243 + (-t258 - t369) * t239 (t232 + t233) * qJDD(1) * pkin(7) + t249 - t286, t239 * t257 + t243 * t251 - t203, t257 * t243 + (-t251 + t369) * t239, pkin(7) * t249 - g(1) * t225 - g(2) * t292 + t159 * t146 + t92 * t182 - t271 * t369, t187 * t116 - t176 * t165 + (qJD(1) * t101 + t64) * t331 - t260 * t235 + (t67 * qJD(1) + t101 * qJDD(1) + t234 * t285 - t235 * t326 + t28) * t239, t187 * t117 - t176 * t167 + (-qJD(1) * t102 - t65) * t331 + t260 * t234 + (-t68 * qJD(1) - t102 * qJDD(1) + t234 * t326 + t235 * t285 - t29) * t239, -t101 * t117 - t102 * t116 - t68 * t165 - t67 * t167 + t203 + qJD(2) * t265 + (t234 * t28 - t235 * t29 - t366) * t243, t29 * t102 + t65 * t68 + t28 * t101 + t64 * t67 + t93 * t187 - t148 * t176 - g(1) * (t244 * pkin(3) + t225) - g(2) * (qJ(4) * t341 + t292) + (-g(1) * (t271 - t353) - g(2) * pkin(3)) * t240, -t141 * t33 - t80 * t94, t140 * t33 + t141 * t34 - t80 * t95 - t81 * t94, -t141 * t169 + t200 * t80 + t239 * t33 - t331 * t94, t140 * t169 + t200 * t81 - t239 * t34 - t331 * t95, t169 * t239 + t200 * t331, t299 * t200 + t297 * t169 + t301 * t239 + t17 * t331 + t130 * t95 + t152 * t34 - t55 * t140 - t103 * t81 - g(1) * t135 - g(2) * t133 + (-t18 * t239 - t200 * t359) * qJD(5), g(1) * t134 - g(2) * t132 + t103 * t80 - t130 * t94 - t55 * t141 + t152 * t33 - t169 * t359 - t18 * t331 - t200 * t267 - t239 * t268, -t20 * t279 + t7 * t76, t20 * t375 + t21 * t279 + t252 * t76 + t276 * t7, t163 * t76 + t189 * t20 + t239 * t7 - t279 * t331, t163 * t276 - t189 * t21 + t239 * t252 + t331 * t375, t163 * t239 + t189 * t331 (-t237 * t9 + t241 * t6) * t189 + (-t237 * t27 + t241 * t26) * t163 + t314 * t239 + t4 * t331 - t50 * t375 - t91 * t252 - t14 * t276 + t45 * t21 - g(1) * t121 - g(2) * t119 + ((-t237 * t26 - t241 * t27) * t189 - t5 * t239) * qJD(6), -t5 * t331 + g(1) * t120 - g(2) * t118 + t11 * t239 + t14 * t76 + t45 * t20 - t50 * t279 + t91 * t7 + (-(-qJD(6) * t27 + t6) * t189 - t26 * t163 - t2 * t239) * t237 + (-(qJD(6) * t26 + t9) * t189 - t27 * t163 - t302 * t239) * t241; 0, 0, 0, -t318, t336 * t246, t320, t319, qJDD(2), pkin(1) * t343 + t291, t365 - t209 + (pkin(1) * t246 + t286) * t243 (-pkin(2) * t239 + t354) * qJDD(1) + ((-t185 - t231) * t239 + (-t181 + t296) * t243) * qJD(1), -t174 * t334 + t266 - 0.2e1 * t352, t209 + 0.2e1 * t229 + 0.2e1 * t230 + (qJD(1) * t174 - g(3)) * t239 + (qJD(1) * t159 - t286) * t243, -t128 * qJ(3) - t185 * qJD(3) - t143 * pkin(2) - t159 * t174 - g(1) * (-pkin(2) * t344 + t197) - g(2) * (-pkin(2) * t345 + t194) - g(3) * t337 - t275 * qJD(1) * pkin(7), -t364 * t305 + qJ(3) * t116 + t323 * t165 + t254 * t234 + (-t372 * t235 - t239 * t86 - t243 * t64) * qJD(1), t364 * t306 + qJ(3) * t117 + t323 * t167 + t254 * t235 + (t372 * t234 + t239 * t87 + t243 * t65) * qJD(1), t165 * t87 + t167 * t86 + (qJD(4) * t167 + t117 * t364 - t335 * t65 - t28) * t235 + (qJD(4) * t165 + t116 * t364 + t335 * t64 - t29) * t234 - t316, t93 * qJ(3) - t65 * t87 - t64 * t86 - g(1) * t197 - g(2) * t194 - g(3) * t284 + t323 * t148 + (-t234 * t65 - t235 * t64) * qJD(4) + (t239 * t286 - t280) * t364, -t274 * t33 - t340 * t94, -t33 * t170 + t274 * t34 + t339 * t94 - t340 * t95, t334 * t94 + t283, t334 * t95 + t262, -t200 * t334, t339 * t103 + t294 * t169 - t17 * t334 + t55 * t170 + t373 * t200 + t205 * t34 + t255 * t216 + t324 * t95, t340 * t103 - t338 * t169 + t18 * t334 + t371 * t200 + t205 * t33 + t255 * t217 - t274 * t55 - t324 * t94, t7 * t100 - t279 * t362, t100 * t252 + t279 * t361 + t362 * t375 - t7 * t99, t189 * t362 + t279 * t334 + t351, -t189 * t361 - t334 * t375 - t357, -t189 * t334 (-t237 * t73 + t241 * t72) * t163 - t124 * t252 + t14 * t99 - t4 * t334 + t361 * t45 - t356 * t375 + (t237 * t288 - t241 * t289) * t189 + t255 * t206 -(t237 * t72 + t241 * t73) * t163 + t124 * t7 + t14 * t100 + t5 * t334 + t362 * t45 - t356 * t279 + (t237 * t289 + t241 * t288) * t189 + t255 * t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, qJDD(2) + t318, -t245 - t349, qJD(2) * t185 + t198 + t266 - t352, t305 - t234 * t349 + (-t165 + t311) * qJD(2), -t306 - t235 * t349 + (-t167 - t309) * qJD(2), -t234 * t116 - t235 * t117 + (-t165 * t235 + t167 * t234) * t335, qJD(1) * t265 + t280 + t316 - t326, 0, 0, 0, 0, 0, -qJD(2) * t95 + t283, qJD(2) * t94 + t262, 0, 0, 0, 0, 0, t351 + qJD(2) * t375 + (-t170 * t327 + t237 * t281 + t374) * t189, -t357 + qJD(2) * t279 + (t281 * t241 + (qJD(6) * t170 - t340) * t237) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167 * t335 + t116 (-t165 + t333) * t335 + t293, -t165 ^ 2 - t167 ^ 2, t65 * t165 + t64 * t167 - t270 + (-g(3) - t290) * t239 + t273, 0, 0, 0, 0, 0, t34 - t383, t33 - t385, 0, 0, 0, 0, 0, -t252 - t382, t7 + t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 * t95, t94 ^ 2 - t95 ^ 2, t33 + t385, -t34 - t383, t169, -g(1) * t132 - g(2) * t134 + t103 * t94 + t18 * t200 + t217 * t227 + t253, g(1) * t133 - g(2) * t135 + t103 * t95 + t17 * t200 - t216 * t227 - t268, t384, t380, t379, t376, t163 -(-t12 * t237 - t358) * t189 + (t163 * t241 - t189 * t328 - t375 * t94) * pkin(5) + t377 (-t13 * t189 - t2) * t237 + (t12 * t189 - t302) * t241 + (-t163 * t237 - t189 * t327 - t279 * t94) * pkin(5) + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t384, t380, t379, t376, t163, t5 * t189 + t377, t4 * t189 - t237 * t2 - t241 * t302 + t378;];
tau_reg  = t1;
