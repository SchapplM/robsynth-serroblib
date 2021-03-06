% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:41
% EndTime: 2019-03-08 22:40:59
% DurationCPUTime: 7.92s
% Computational Cost: add. (5092->579), mult. (12990->827), div. (0->0), fcn. (11076->14), ass. (0->292)
t207 = sin(pkin(7));
t214 = sin(qJ(3));
t330 = qJD(3) * t214;
t303 = t207 * t330;
t189 = pkin(3) * t303;
t218 = cos(qJ(3));
t371 = qJ(4) * t218;
t269 = pkin(10) * t214 - t371;
t328 = qJD(4) * t214;
t230 = qJD(3) * t269 - t328;
t208 = sin(pkin(6));
t215 = sin(qJ(2));
t336 = qJD(1) * t215;
t308 = t208 * t336;
t401 = -t189 + (-t230 + t308) * t207;
t210 = cos(pkin(7));
t349 = t215 * t218;
t219 = cos(qJ(2));
t350 = t214 * t219;
t253 = t210 * t349 + t350;
t125 = t253 * t208;
t110 = qJD(1) * t125;
t356 = t210 * t214;
t201 = pkin(2) * t356;
t362 = t207 * t218;
t387 = pkin(4) + pkin(9);
t400 = (t362 * t387 + t201) * qJD(3) - t110;
t321 = qJD(2) * qJD(3);
t300 = t214 * t321;
t279 = t207 * t300;
t317 = qJDD(2) * t218;
t397 = -t207 * t317 + t279;
t335 = qJD(2) * t207;
t153 = pkin(9) * t335 + t308;
t359 = t208 * t219;
t307 = qJD(1) * t359;
t166 = qJD(2) * pkin(2) + t307;
t211 = cos(pkin(6));
t337 = qJD(1) * t211;
t309 = t207 * t337;
t355 = t210 * t218;
t69 = t214 * t153 - t166 * t355 - t218 * t309;
t344 = qJD(4) + t69;
t352 = t214 * t215;
t311 = t210 * t352;
t285 = t208 * t311;
t329 = qJD(3) * t218;
t301 = t210 * t329;
t399 = -pkin(2) * t301 - qJD(1) * t285 + t218 * t307;
t331 = qJD(2) * t219;
t386 = pkin(9) * t207;
t398 = qJDD(2) * t386 + (qJD(1) * t331 + qJDD(1) * t215) * t208 + qJD(3) * t309;
t334 = qJD(2) * t210;
t196 = qJD(3) + t334;
t213 = sin(qJ(5));
t217 = cos(qJ(5));
t332 = qJD(2) * t218;
t305 = t207 * t332;
t127 = t196 * t213 + t217 * t305;
t117 = qJD(6) + t127;
t333 = qJD(2) * t214;
t306 = t207 * t333;
t180 = qJD(5) + t306;
t282 = t213 * t305;
t129 = t196 * t217 - t282;
t212 = sin(qJ(6));
t216 = cos(qJ(6));
t78 = t129 * t212 - t216 * t180;
t396 = t180 * t78;
t220 = -pkin(3) - pkin(10);
t372 = qJ(4) * t214;
t257 = t218 * t220 - t372;
t102 = (-pkin(2) + t257) * t207;
t326 = qJD(5) * t217;
t327 = qJD(5) * t213;
t364 = t207 * t214;
t198 = pkin(9) * t364;
t310 = -pkin(2) * t218 - pkin(3);
t89 = pkin(4) * t364 + t198 + (-pkin(10) + t310) * t210;
t395 = t102 * t327 - t400 * t213 + t217 * t401 - t89 * t326;
t202 = t210 * qJD(4);
t374 = -t303 * t387 + t202 - t399;
t394 = pkin(9) * t303 + t399;
t345 = pkin(4) * t306 + t344;
t193 = t210 * t337;
t270 = -pkin(3) * t218 - t372;
t77 = t193 + (qJD(2) * t270 - t166) * t207;
t393 = t77 * t306 + qJDD(4);
t319 = qJDD(2) * t210;
t195 = qJDD(3) + t319;
t299 = t218 * t321;
t318 = qJDD(2) * t214;
t240 = t299 + t318;
t233 = t240 * t207;
t194 = qJDD(1) * t359;
t360 = t208 * t215;
t304 = qJD(2) * t360;
t281 = qJD(1) * t304;
t130 = qJDD(2) * pkin(2) + t194 - t281;
t141 = t166 * t356;
t320 = qJDD(1) * t211;
t298 = t207 * t320;
t259 = qJD(3) * t141 - t130 * t355 + t153 * t329 + t214 * t398 - t218 * t298;
t237 = qJDD(4) + t259;
t11 = pkin(4) * t233 + t195 * t220 + t237;
t42 = t196 * t220 + t345;
t66 = t193 + (qJD(2) * t257 - t166) * t207;
t21 = t213 * t42 + t217 * t66;
t192 = t210 * t320;
t340 = pkin(3) * t279 + t192;
t27 = (qJD(2) * t230 + qJDD(2) * t257 - t130) * t207 + t340;
t228 = -qJD(5) * t21 + t217 * t11 - t213 * t27;
t134 = qJDD(5) + t233;
t54 = -qJD(5) * t127 + t217 * t195 + t213 * t397;
t80 = t129 * t216 + t180 * t212;
t24 = qJD(6) * t80 - t216 * t134 + t212 * t54;
t373 = sin(pkin(12));
t293 = t373 * t215;
t209 = cos(pkin(12));
t357 = t209 * t219;
t146 = t211 * t357 - t293;
t292 = t373 * t219;
t358 = t209 * t215;
t147 = t211 * t358 + t292;
t361 = t208 * t209;
t315 = t207 * t361;
t60 = -t146 * t355 + t147 * t214 + t218 * t315;
t148 = -t211 * t292 - t358;
t149 = -t211 * t293 + t357;
t294 = t208 * t373;
t277 = t207 * t294;
t62 = -t148 * t355 + t149 * t214 - t218 * t277;
t145 = -t207 * t359 + t210 * t211;
t312 = t208 * t352;
t313 = t210 * t359;
t90 = -t211 * t362 - t218 * t313 + t312;
t64 = t145 * t213 - t90 * t217;
t93 = -t146 * t207 - t210 * t361;
t94 = -t148 * t207 + t210 * t294;
t246 = g(1) * (-t213 * t94 + t217 * t62) + g(2) * (-t213 * t93 + t217 * t60) - g(3) * t64;
t4 = -t134 * pkin(5) - t228;
t392 = t117 * (pkin(5) * t129 + pkin(11) * t117) + t246 + t4;
t172 = pkin(5) * t213 - pkin(11) * t217 + qJ(4);
t275 = pkin(5) * t217 + pkin(11) * t213;
t55 = -qJD(5) * t282 + t213 * t195 + (qJD(5) * t196 - t397) * t217;
t50 = qJDD(6) + t55;
t391 = t117 * (t275 * qJD(5) - (-pkin(4) - t275) * t306 + t344) + t172 * t50;
t126 = (t218 * t219 - t311) * t208;
t339 = pkin(9) * t362 + t201;
t341 = t339 * qJD(3) - t110;
t72 = t146 * t218 - t147 * t356;
t74 = t148 * t218 - t149 * t356;
t390 = g(1) * t74 + g(2) * t72 + g(3) * t126 + t341 * t196;
t239 = t300 - t317;
t204 = t207 ^ 2;
t221 = qJD(2) ^ 2;
t348 = t215 * t221;
t287 = t204 * t208 * t348;
t368 = t145 * t207;
t252 = t210 * t350 + t349;
t51 = t211 * t303 + (qJD(2) * t253 + qJD(3) * t252) * t208;
t389 = -t90 * t195 - t51 * t196 - t218 * t287 + t239 * t368;
t16 = pkin(11) * t180 + t21;
t185 = t196 * qJ(4);
t70 = t218 * t153 + t214 * t309 + t141;
t59 = pkin(4) * t305 + t70;
t47 = t185 + t59;
t25 = pkin(5) * t127 - pkin(11) * t129 + t47;
t268 = t16 * t212 - t216 * t25;
t316 = t213 * t11 + t217 * t27 + t42 * t326;
t258 = -t327 * t66 + t316;
t3 = pkin(11) * t134 + t258;
t181 = t195 * qJ(4);
t184 = t196 * qJD(4);
t260 = -t130 * t356 + t153 * t330 - t166 * t301 - t214 * t298 - t218 * t398;
t17 = -t181 - t184 + t260;
t12 = -pkin(4) * t397 - t17;
t6 = pkin(5) * t55 - pkin(11) * t54 + t12;
t1 = -t268 * qJD(6) + t212 * t6 + t216 * t3;
t267 = t217 * t102 + t213 * t89;
t388 = -qJD(5) * t267 + t213 * t401 + t400 * t217;
t385 = t195 * pkin(3);
t302 = t207 * t329;
t384 = -pkin(5) * t302 - t388;
t191 = pkin(3) * t306;
t107 = t269 * t335 + t191;
t383 = t217 * t107 + t213 * t59;
t382 = t117 * t78;
t380 = t212 * t50;
t379 = t216 * t50;
t323 = qJD(6) * t216;
t324 = qJD(6) * t212;
t23 = -t129 * t324 + t212 * t134 + t180 * t323 + t216 * t54;
t378 = t217 * t23;
t377 = t23 * t212;
t376 = t80 * t117;
t370 = t127 * t180;
t369 = t129 * t180;
t367 = t180 * t213;
t366 = t180 * t217;
t365 = t204 * t221;
t363 = t207 * t217;
t354 = t212 * t220;
t353 = t213 * t134;
t351 = t214 * t216;
t347 = t216 * t220;
t346 = t220 * t134;
t343 = qJDD(1) - g(3);
t342 = -t202 + t394;
t205 = t214 ^ 2;
t338 = -t218 ^ 2 + t205;
t325 = qJD(5) * t220;
t322 = qJD(3) - t196;
t314 = t213 * t362;
t131 = -t210 * qJ(4) - t339;
t290 = t117 * t216;
t289 = t196 + t334;
t288 = t195 + t319;
t286 = t214 * t218 * t365;
t101 = pkin(4) * t362 - t131;
t283 = t207 * t304;
t150 = t210 * t213 + t217 * t362;
t151 = t210 * t217 - t314;
t48 = pkin(5) * t150 - pkin(11) * t151 + t101;
t278 = -pkin(11) * t302 - qJD(6) * t48 + t395;
t44 = pkin(11) * t364 + t267;
t97 = -qJD(5) * t150 + t213 * t303;
t98 = -qJD(5) * t314 + t210 * t326 - t217 * t303;
t276 = -t98 * pkin(5) + t97 * pkin(11) + qJD(6) * t44 - t374;
t274 = g(1) * t149 + g(2) * t147;
t113 = (t212 * t218 + t213 * t351) * t335;
t271 = -t216 * t327 - t113;
t8 = t16 * t216 + t212 * t25;
t65 = t145 * t217 + t213 * t90;
t91 = t208 * t252 + t211 * t364;
t31 = t212 * t91 + t216 * t65;
t30 = -t212 * t65 + t216 * t91;
t20 = -t213 * t66 + t217 * t42;
t266 = -t102 * t213 + t217 * t89;
t132 = (-pkin(2) + t270) * t207;
t264 = qJD(3) * (-qJD(2) * t132 - t77);
t263 = t214 * t281;
t262 = t218 * t281;
t256 = -t117 * t323 - t380;
t255 = -t117 * t324 + t379;
t254 = -t151 * t212 + t207 * t351;
t100 = t151 * t216 + t212 * t364;
t251 = t180 * t80;
t109 = t189 + (-qJ(4) * t329 - t328) * t207;
t37 = (-pkin(3) * t317 - qJ(4) * t240 - qJD(2) * t328 - t130) * t207 + t340;
t247 = qJD(2) * t109 + qJDD(2) * t132 + t37;
t245 = g(1) * t62 + g(2) * t60 + g(3) * t90;
t61 = t146 * t356 + t147 * t218 - t214 * t315;
t63 = t149 * t218 + (t148 * t210 + t277) * t214;
t244 = g(1) * t63 + g(2) * t61 + g(3) * t91;
t71 = t146 * t214 + t147 * t355;
t73 = t148 * t214 + t149 * t355;
t242 = g(1) * t73 + g(2) * t71 + g(3) * t125;
t236 = t12 - t244;
t234 = -g(3) * t360 - t274;
t15 = -pkin(5) * t180 - t20;
t232 = -pkin(11) * t50 + (t15 + t20) * t117;
t231 = (t322 * t332 + t318) * t207;
t2 = -qJD(6) * t8 - t212 * t3 + t216 * t6;
t229 = qJD(6) * t117 * t220 + t244;
t52 = -qJD(2) * t285 - qJD(3) * t312 + (t208 * t331 + (t207 * t211 + t313) * qJD(3)) * t218;
t227 = t91 * t195 + t52 * t196 - t214 * t287;
t225 = (pkin(11) * t305 - qJD(6) * t172 + t383) * t117 + t245;
t224 = -t244 - t260;
t223 = t245 - t259;
t222 = t196 * t70 + t223;
t139 = -qJ(4) * t305 + t191;
t133 = t210 * t310 + t198;
t124 = t217 * t134;
t115 = -t207 * t166 + t193;
t112 = t212 * t213 * t306 - t216 * t305;
t87 = -t207 * t130 + t192;
t82 = t125 * t213 + t360 * t363;
t57 = -t185 - t70;
t56 = -pkin(3) * t196 + t344;
t46 = t149 * t363 + t213 * t73;
t45 = t147 * t363 + t213 * t71;
t43 = -pkin(5) * t364 - t266;
t41 = qJD(6) * t100 + t212 * t97 - t216 * t302;
t40 = qJD(6) * t254 + t212 * t302 + t216 * t97;
t36 = t213 * t62 + t217 * t94;
t34 = t213 * t60 + t217 * t93;
t28 = -pkin(5) * t305 + t107 * t213 - t217 * t59;
t22 = t237 - t385;
t19 = -qJD(5) * t64 + t51 * t213 + t217 * t283;
t18 = qJD(5) * t65 + t213 * t283 - t51 * t217;
t5 = [t343, 0 (qJDD(2) * t219 - t348) * t208 (-qJDD(2) * t215 - t219 * t221) * t208, 0, 0, 0, 0, 0, t389, t145 * t233 - t227 ((t214 * t90 + t218 * t91) * qJDD(2) + (t214 * t51 + t218 * t52 + (-t214 * t91 + t218 * t90) * qJD(3)) * qJD(2)) * t207, -t389, -t240 * t368 + t227, t145 * t37 - t17 * t91 + t22 * t90 + t283 * t77 + t51 * t56 - t52 * t57 - g(3), 0, 0, 0, 0, 0, t127 * t52 - t134 * t64 - t18 * t180 + t55 * t91, t129 * t52 - t134 * t65 - t180 * t19 + t54 * t91, 0, 0, 0, 0, 0 (-qJD(6) * t31 - t19 * t212 + t52 * t216) * t117 + t30 * t50 + t18 * t78 + t64 * t24 -(qJD(6) * t30 + t19 * t216 + t52 * t212) * t117 - t31 * t50 + t18 * t80 + t64 * t23; 0, qJDD(2), -g(1) * t148 - g(2) * t146 - g(3) * t359 + t194, -t343 * t360 + t274 (qJDD(2) * t205 + 0.2e1 * t214 * t299) * t204, 0.2e1 * (t214 * t317 - t321 * t338) * t204 (t214 * t288 + t289 * t329) * t207 (t218 * t288 - t289 * t330) * t207, t195 * t210 (pkin(2) * t355 - t198) * t195 - t259 * t210 + (t115 * t330 - t218 * t87) * t207 + (-pkin(2) * t239 + t262) * t204 - t390, -t339 * t195 + t260 * t210 + (t115 * t329 + t214 * t87) * t207 + t394 * t196 + (-pkin(2) * t240 - t263) * t204 + t242 ((qJD(3) * t56 - qJDD(2) * t131 - t17) * t218 + (qJD(3) * t57 + qJDD(2) * t133 + t22) * t214 + ((qJD(3) * t133 - t342) * t218 + (qJD(3) * t131 + t341) * t214) * qJD(2) + t234) * t207, -t204 * t262 + t133 * t195 + t210 * t22 + (t214 * t264 + t218 * t247) * t207 + t390, t204 * t263 - t131 * t195 - t17 * t210 - t342 * t196 + (-t214 * t247 + t218 * t264) * t207 - t242, t37 * t132 + t77 * t109 + t17 * t131 + t22 * t133 - g(1) * (pkin(2) * t148 + pkin(3) * t74 + qJ(4) * t73 + t149 * t386) - g(2) * (pkin(2) * t146 + pkin(3) * t72 + qJ(4) * t71 + t147 * t386) - g(3) * (pkin(3) * t126 + qJ(4) * t125) + t342 * t57 + t341 * t56 + (-t77 * t207 * t336 - g(3) * (pkin(2) * t219 + t215 * t386)) * t208, t129 * t97 + t151 * t54, -t127 * t97 - t129 * t98 - t150 * t54 - t151 * t55, t151 * t134 + t97 * t180 + (t129 * t329 + t214 * t54) * t207, -t150 * t134 - t98 * t180 + (-t127 * t329 - t214 * t55) * t207 (t134 * t214 + t180 * t329) * t207, t266 * t134 + t101 * t55 + t12 * t150 + t47 * t98 - g(1) * t46 - g(2) * t45 - g(3) * t82 + (t20 * t329 + t214 * t228) * t207 + t388 * t180 + t374 * t127, -t267 * t134 + t101 * t54 + t12 * t151 + t47 * t97 - t242 * t217 + (-t316 * t214 - t21 * t329 + (qJD(5) * t214 * t66 - t234) * t213) * t207 + t395 * t180 + t374 * t129, t100 * t23 + t40 * t80, -t100 * t24 + t23 * t254 - t40 * t78 - t41 * t80, t100 * t50 + t117 * t40 + t150 * t23 + t80 * t98, -t117 * t41 - t150 * t24 + t254 * t50 - t78 * t98, t117 * t98 + t150 * t50 (-t212 * t44 + t216 * t48) * t50 + t2 * t150 - t268 * t98 + t43 * t24 - t4 * t254 + t15 * t41 - g(1) * (t212 * t74 + t216 * t46) - g(2) * (t212 * t72 + t216 * t45) - g(3) * (t126 * t212 + t216 * t82) + t384 * t78 + (t212 * t278 - t216 * t276) * t117 -(t212 * t48 + t216 * t44) * t50 - t1 * t150 - t8 * t98 + t43 * t23 + t4 * t100 + t15 * t40 - g(1) * (-t212 * t46 + t216 * t74) - g(2) * (-t212 * t45 + t216 * t72) - g(3) * (t126 * t216 - t212 * t82) + t384 * t80 + (t212 * t276 + t216 * t278) * t117; 0, 0, 0, 0, -t286, t338 * t365, t231 (-t322 * t333 + t317) * t207, t195, -t115 * t306 + t222, -t115 * t305 - t196 * t69 - t224 ((-pkin(3) * t214 + t371) * qJDD(2) + ((-qJ(4) * qJD(3) - t57 - t70) * t214 + (-pkin(3) * qJD(3) + t344 - t56) * t218) * qJD(2)) * t207, -t139 * t305 - t222 - 0.2e1 * t385 + t393, 0.2e1 * t181 + t184 + t344 * t196 + (t139 * t214 + t218 * t77) * t335 + t224, -t17 * qJ(4) - t22 * pkin(3) - t77 * t139 - t56 * t70 - g(1) * (-pkin(3) * t62 + qJ(4) * t63) - g(2) * (-pkin(3) * t60 + qJ(4) * t61) - g(3) * (-pkin(3) * t90 + qJ(4) * t91) - t344 * t57, -t129 * t367 + t54 * t217 (-t55 - t369) * t217 + (-t54 + t370) * t213, -t180 * t327 + t124 + (-t129 * t218 - t214 * t367) * t335, -t180 * t326 - t353 + (t127 * t218 - t214 * t366) * t335, -t180 * t305, -t20 * t305 + qJ(4) * t55 + t345 * t127 + (t346 + (t47 - t59) * t180) * t217 + ((t107 - t325) * t180 + t236) * t213, qJ(4) * t54 + t383 * t180 + t21 * t305 + t345 * t129 + (-t180 * t47 - t346) * t213 + (-t180 * t325 + t236) * t217, t216 * t378 + (-t217 * t324 + t271) * t80, t80 * t112 + t113 * t78 + (t212 * t80 + t216 * t78) * t327 + (-t377 - t216 * t24 + (t212 * t78 - t216 * t80) * qJD(6)) * t217, t23 * t213 + t271 * t117 + (t251 + t255) * t217, -t24 * t213 + (t212 * t327 + t112) * t117 + (t256 - t396) * t217, t117 * t366 + t50 * t213, -t15 * t112 - t28 * t78 + t391 * t216 + t225 * t212 + (-t50 * t354 + t2 + (-t15 * t212 + t220 * t78) * qJD(5) - t229 * t216) * t213 + (-t268 * t306 + t15 * t323 + t4 * t212 - t220 * t24 + (-t117 * t354 - t268) * qJD(5)) * t217, -t15 * t113 - t28 * t80 - t391 * t212 + t225 * t216 + (-t50 * t347 - t1 + (-t15 * t216 + t220 * t80) * qJD(5) + t229 * t212) * t213 + (-t8 * t306 - t15 * t324 + t4 * t216 - t220 * t23 + (-t117 * t347 - t8) * qJD(5)) * t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, t195 + t286, -t196 ^ 2 - t205 * t365, t196 * t57 - t223 - t385 + t393, 0, 0, 0, 0, 0, -t196 * t127 - t180 * t367 + t124, -t196 * t129 - t180 * t366 - t353, 0, 0, 0, 0, 0, -t217 * t24 + (-t216 * t196 - t212 * t366) * t117 + (t256 + t396) * t213, -t378 + (t212 * t196 - t216 * t366) * t117 + (t251 - t255) * t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t127, -t127 ^ 2 + t129 ^ 2, t54 + t370, t369 - t55, t134, -t47 * t129 + t21 * t180 + t228 - t246, g(1) * t36 + g(2) * t34 + g(3) * t65 + t127 * t47 + t180 * t20 - t258, t290 * t80 + t377 (t23 - t382) * t216 + (-t24 - t376) * t212, t117 * t290 - t80 * t129 + t380, -t117 ^ 2 * t212 + t78 * t129 + t379, -t117 * t129, -pkin(5) * t24 + t129 * t268 - t21 * t78 + t232 * t212 - t216 * t392, -pkin(5) * t23 + t8 * t129 - t21 * t80 + t212 * t392 + t232 * t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t78, -t78 ^ 2 + t80 ^ 2, t23 + t382, -t24 + t376, t50, t8 * t117 - t15 * t80 - g(1) * (-t212 * t36 + t216 * t63) - g(2) * (-t212 * t34 + t216 * t61) - g(3) * t30 + t2, -t268 * t117 + t15 * t78 - g(1) * (-t212 * t63 - t216 * t36) - g(2) * (-t212 * t61 - t216 * t34) + g(3) * t31 - t1;];
tau_reg  = t5;
