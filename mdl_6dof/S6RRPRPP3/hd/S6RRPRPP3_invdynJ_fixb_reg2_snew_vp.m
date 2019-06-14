% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRPP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:32:14
% EndTime: 2019-05-06 12:32:38
% DurationCPUTime: 9.85s
% Computational Cost: add. (22719->468), mult. (49436->556), div. (0->0), fcn. (34864->8), ass. (0->280)
t297 = cos(qJ(2));
t341 = t297 * qJD(1);
t281 = -qJD(4) + t341;
t343 = qJD(5) * t281;
t296 = cos(qJ(4));
t279 = t281 ^ 2;
t291 = sin(pkin(9));
t292 = cos(pkin(9));
t294 = sin(qJ(2));
t345 = qJD(1) * t294;
t264 = t292 * qJD(2) - t291 * t345;
t265 = t291 * qJD(2) + t292 * t345;
t293 = sin(qJ(4));
t237 = t293 * t264 + t296 * t265;
t392 = t237 ^ 2;
t403 = -t392 - t279;
t340 = qJD(1) * qJD(2);
t283 = t294 * t340;
t338 = t297 * qJDD(1);
t269 = -t283 + t338;
t266 = -qJDD(4) + t269;
t235 = -t296 * t264 + t293 * t265;
t358 = t237 * t235;
t311 = t266 - t358;
t411 = t293 * t311;
t104 = t296 * t403 + t411;
t410 = t296 * t311;
t433 = t293 * t403 - t410;
t67 = t104 * t292 - t291 * t433;
t484 = -pkin(2) * t67 - pkin(3) * t104;
t497 = 0.2e1 * t343 - t484;
t334 = t297 * t340;
t339 = t294 * qJDD(1);
t312 = t334 + t339;
t245 = t291 * qJDD(2) + t292 * t312;
t307 = -t292 * qJDD(2) + t291 * t312;
t168 = qJD(4) * t237 + t245 * t293 + t296 * t307;
t357 = t237 * t281;
t409 = t168 - t357;
t310 = t266 + t358;
t172 = t296 * t310;
t393 = t235 ^ 2;
t405 = -t279 - t393;
t100 = -t405 * t293 + t172;
t365 = t310 * t293;
t435 = t405 * t296 + t365;
t469 = t100 * t291 + t292 * t435;
t62 = t100 * t292 - t291 * t435;
t496 = pkin(7) * (t294 * t409 + t297 * t469) + pkin(1) * t62;
t303 = t296 * t245 - t293 * t307;
t169 = -t235 * qJD(4) + t303;
t359 = t235 * t281;
t407 = t169 + t359;
t468 = t104 * t291 + t292 * t433;
t486 = t297 * t468;
t490 = pkin(1) * t67;
t495 = t490 + pkin(7) * (-t294 * t407 + t486);
t489 = qJ(3) * t62;
t488 = qJ(3) * t67;
t483 = -pkin(2) * t62 - pkin(3) * t100;
t209 = -t392 + t279;
t124 = t209 * t293 + t172;
t404 = -t359 + t169;
t436 = t296 * t209 - t365;
t493 = t294 * (t124 * t292 + t291 * t436) + t297 * t404;
t186 = t393 - t392;
t434 = -t409 * t293 + t407 * t296;
t419 = t409 * t296;
t84 = t407 * t293 + t419;
t492 = -t297 * t186 + t294 * (t291 * t434 + t292 * t84);
t487 = qJ(3) * t468;
t471 = -t291 * t84 + t292 * t434;
t482 = -pkin(2) * t409 + qJ(3) * t469;
t481 = pkin(2) * t407 + t487;
t208 = t393 - t279;
t121 = t208 * t293 - t410;
t127 = t208 * t296 + t411;
t408 = t168 + t357;
t431 = t294 * (t121 * t291 - t127 * t292) - t297 * t408;
t476 = t124 * t291 - t292 * t436;
t460 = pkin(8) * t100;
t459 = pkin(8) * t104;
t438 = t404 * t293 - t296 * t408;
t439 = -t293 * t408 - t296 * t404;
t444 = t291 * t438 + t292 * t439;
t473 = qJ(3) * t444;
t472 = -pkin(2) * t444 - pkin(3) * t439;
t402 = -t392 - t393;
t445 = -t291 * t439 + t292 * t438;
t470 = -pkin(2) * t402 + qJ(3) * t445;
t467 = t121 * t292 + t127 * t291;
t464 = pkin(7) * (t294 * t402 + t297 * t445) - pkin(1) * t444;
t458 = pkin(8) * t433;
t457 = pkin(8) * t435;
t456 = pkin(8) * t439;
t449 = -pkin(3) * t409 + t457;
t448 = -pkin(3) * t407 - t458;
t447 = -pkin(3) * t402 + pkin(8) * t438;
t442 = (t409 - t357) * pkin(4);
t441 = qJ(5) * t407;
t299 = qJD(1) ^ 2;
t295 = sin(qJ(1));
t298 = cos(qJ(1));
t324 = t298 * g(1) + t295 * g(2);
t377 = qJDD(1) * pkin(7);
t257 = -t299 * pkin(1) - t324 + t377;
t319 = -t297 * pkin(2) - t294 * qJ(3);
t328 = t299 * t319 + t257;
t383 = t297 * g(3);
t390 = qJD(2) ^ 2;
t201 = -qJDD(2) * pkin(2) - t390 * qJ(3) + t294 * t328 + qJDD(3) + t383;
t246 = -pkin(3) * t341 - t265 * pkin(8);
t391 = t264 ^ 2;
t148 = t307 * pkin(3) - t391 * pkin(8) + t265 * t246 + t201;
t302 = t168 * pkin(4) + t148 - t441;
t301 = 0.2e1 * qJD(5) * t237 - t302;
t424 = qJ(5) * t402;
t423 = qJ(5) * t405;
t313 = (t235 * t293 + t237 * t296) * t281;
t415 = t291 * t313;
t356 = t264 * t265;
t314 = -t269 + t356;
t414 = t291 * t314;
t413 = t292 * t313;
t412 = t292 * t314;
t406 = t393 * pkin(5) - 0.2e1 * qJD(6) * t235;
t250 = t264 * t341;
t223 = -t250 - t245;
t251 = t265 * t341;
t221 = -t307 - t251;
t185 = t235 * pkin(4) - t237 * qJ(5);
t331 = t295 * g(1) - t298 * g(2);
t256 = qJDD(1) * pkin(1) + t299 * pkin(7) + t331;
t268 = 0.2e1 * t334 + t339;
t195 = (-t269 + t283) * pkin(2) - t268 * qJ(3) - t256;
t384 = t294 * g(3);
t202 = -pkin(2) * t390 + qJDD(2) * qJ(3) + t297 * t328 - t384;
t389 = 2 * qJD(3);
t146 = -t292 * t195 + t291 * t202 + t265 * t389;
t95 = t314 * pkin(3) + pkin(8) * t223 - t146;
t147 = t291 * t195 + t292 * t202 + t264 * t389;
t97 = -pkin(3) * t391 - pkin(8) * t307 + t246 * t341 + t147;
t56 = t293 * t97 - t296 * t95;
t46 = t266 * pkin(4) - t279 * qJ(5) + t237 * t185 + qJDD(5) + t56;
t306 = t169 * pkin(5) + qJ(6) * t310 + t46;
t388 = 0.2e1 * qJD(6);
t38 = (-pkin(5) * t235 + t388) * t281 + t306;
t382 = pkin(4) + qJ(6);
t275 = -0.2e1 * t343;
t206 = t237 * pkin(5) + t281 * qJ(6);
t57 = t293 * t95 + t296 * t97;
t320 = -t279 * pkin(4) - t266 * qJ(5) - t235 * t185 + t57;
t308 = -t168 * pkin(5) - qJ(6) * t393 - t281 * t206 + qJDD(6) + t320;
t39 = t275 + t308;
t399 = qJ(5) * t39 - t38 * t382;
t134 = qJ(5) * t408;
t398 = -t382 * t404 - t134;
t173 = qJ(5) * t311;
t397 = -t382 * t403 - t173;
t396 = -t310 * t382 + t423;
t354 = t281 * t293;
t207 = t237 * t354;
t353 = t281 * t296;
t337 = t235 * t353;
t321 = -t207 + t337;
t395 = t291 * t321 + t413;
t315 = t168 * t293 - t337;
t322 = -t296 * t168 - t235 * t354;
t75 = t291 * t315 + t292 * t322;
t323 = t293 * t169 - t237 * t353;
t346 = t296 * t169 + t207;
t76 = t291 * t346 + t292 * t323;
t253 = t297 * t266;
t394 = t294 * (t292 * t321 - t415) + t253;
t336 = t297 * t358;
t327 = t294 * (-t291 * t322 + t292 * t315) + t336;
t326 = t294 * (-t291 * t323 + t292 * t346) - t336;
t263 = t265 ^ 2;
t387 = pkin(4) * t281;
t386 = pkin(4) * t293;
t385 = pkin(4) * t296;
t36 = t293 * t57 - t296 * t56;
t380 = t291 * t36;
t379 = t292 * t36;
t370 = t148 * t293;
t369 = t148 * t296;
t368 = t168 * qJ(6);
t362 = t201 * t291;
t361 = t201 * t292;
t224 = t269 + t356;
t352 = t291 * t224;
t351 = t292 * t224;
t280 = t297 * t299 * t294;
t350 = t294 * (qJDD(2) + t280);
t347 = t297 * (-t280 + qJDD(2));
t335 = t297 * t356;
t333 = qJ(5) * t293 + pkin(3);
t332 = -pkin(4) * t404 - t134;
t37 = t293 * t56 + t296 * t57;
t330 = -0.2e1 * qJD(5) - t387;
t91 = t146 * t291 + t292 * t147;
t242 = t294 * t257 + t383;
t243 = t297 * t257 - t384;
t329 = t294 * t242 + t297 * t243;
t45 = t275 + t320;
t325 = -pkin(4) * t46 + qJ(5) * t45;
t318 = t146 * t292 - t147 * t291;
t317 = -pkin(1) + t319;
t309 = -pkin(4) * t403 - t173 + t320;
t305 = pkin(4) * t310 - t423 + t46;
t304 = t281 * t388 + t306;
t300 = t301 + t406;
t289 = t297 ^ 2;
t288 = t294 ^ 2;
t285 = t289 * t299;
t284 = t288 * t299;
t270 = -0.2e1 * t283 + t338;
t260 = t297 * t269;
t249 = -t263 - t285;
t248 = -t263 + t285;
t247 = -t285 + t391;
t239 = -t285 - t391;
t222 = -t250 + t245;
t220 = -t251 + t307;
t213 = -t263 - t391;
t192 = -t291 * t249 + t351;
t191 = t292 * t249 + t352;
t182 = t292 * t239 - t414;
t181 = t291 * t239 + t412;
t167 = t221 * t292 - t223 * t291;
t160 = (t235 * t296 - t237 * t293) * t281;
t140 = (-qJD(4) + t281) * t235 + t303;
t92 = pkin(5) * t310 - qJ(5) * t409;
t77 = t369 - t459;
t66 = t370 + t460;
t59 = -pkin(5) * t311 + t382 * t407;
t58 = t370 + t448;
t54 = t237 * t330 + t302;
t53 = -t369 + t449;
t44 = -t301 + t442;
t43 = pkin(4) * t357 + qJ(5) * t140 + t301;
t42 = t46 - t424;
t41 = -pkin(4) * t402 + t45;
t40 = (-t206 + t330) * t237 + t302 + t368 - t406;
t35 = -t140 * t386 + t296 * t43 + t459;
t34 = qJ(5) * t419 - t293 * t44 - t460;
t33 = (t206 + t387) * t237 + t300 + t441 + pkin(5) * t403 - t368;
t32 = -pkin(3) * t148 + pkin(8) * t37;
t31 = t458 + t293 * t43 + (pkin(3) + t385) * t140;
t30 = (-t409 - t168) * qJ(6) - t442 + t300 + pkin(5) * t405 + t237 * t206;
t29 = t296 * t44 + t333 * t409 - t457;
t28 = -t424 + (t404 - t359) * pkin(5) + t304;
t27 = -pkin(5) * t408 - t382 * t402 + t39;
t26 = -t36 - t456;
t25 = t293 * t46 + t296 * t45;
t24 = t293 * t45 - t296 * t46;
t23 = t37 + t447;
t22 = -t293 * t30 + t296 * t92 + t460;
t21 = -t293 * t59 + t296 * t33 + t459;
t20 = pkin(5) * t38 - qJ(5) * t40;
t19 = t293 * t92 + t296 * t30 + t449;
t18 = t293 * t38 + t296 * t39;
t17 = t293 * t39 - t296 * t38;
t16 = t293 * t33 + t296 * t59 - t448;
t15 = -t293 * t41 + t296 * t42 - t456;
t14 = t292 * t37 - t380;
t13 = t291 * t37 + t379;
t12 = t293 * t42 + t296 * t41 + t447;
t11 = -pkin(8) * t24 + (-qJ(5) * t296 + t386) * t54;
t10 = pkin(5) * t39 - t382 * t40;
t9 = -t24 * t291 + t25 * t292;
t8 = t24 * t292 + t25 * t291;
t7 = -t27 * t293 + t28 * t296 - t456;
t6 = pkin(8) * t25 + (-t333 - t385) * t54;
t5 = t27 * t296 + t28 * t293 + t447;
t4 = -t17 * t291 + t18 * t292;
t3 = t17 * t292 + t18 * t291;
t2 = -pkin(8) * t17 - t10 * t293 + t20 * t296;
t1 = -pkin(3) * t40 + pkin(8) * t18 + t10 * t296 + t20 * t293;
t47 = [0, 0, 0, 0, 0, qJDD(1), t331, t324, 0, 0, t268 * t294, t268 * t297 + t270 * t294, t350 + t297 * (-t284 + t390), -t294 * t334 + t260, t294 * (t285 - t390) + t347, 0, t297 * t256 + pkin(1) * t270 + pkin(7) * (t297 * (-t285 - t390) - t350), -t294 * t256 - pkin(1) * t268 + pkin(7) * (-t347 - t294 * (-t284 - t390)), pkin(1) * (t284 + t285) + (t288 + t289) * t377 + t329, pkin(1) * t256 + pkin(7) * t329, t294 * (t245 * t292 + t251 * t291) + t335, t294 * (-t220 * t292 - t222 * t291) + t297 * (-t263 + t391), t294 * (-t248 * t291 + t412) + t297 * t223, t294 * (t250 * t292 + t291 * t307) - t335, t294 * (t247 * t292 + t352) - t297 * t221, t260 + t294 * (-t264 * t292 - t265 * t291) * t341, t294 * (-qJ(3) * t181 + t362) + t297 * (-pkin(2) * t181 + t146) - pkin(1) * t181 + pkin(7) * (t182 * t297 + t220 * t294), t294 * (-qJ(3) * t191 + t361) + t297 * (-pkin(2) * t191 + t147) - pkin(1) * t191 + pkin(7) * (t192 * t297 + t222 * t294), t294 * t318 + pkin(7) * (t167 * t297 + t213 * t294) + t317 * (t221 * t291 + t223 * t292), pkin(7) * (t201 * t294 + t297 * t91) - t317 * t318, t326, -t492, -t493, t327, -t431, t394, t294 * (-t291 * t53 + t292 * t66 + t489) + t297 * (-t483 + t56) + t496, t294 * (-t291 * t58 + t292 * t77 - t488) + t297 * (t484 + t57) - t495, t294 * (-t23 * t291 + t26 * t292 - t473) + t297 * t472 + t464, t294 * (-pkin(8) * t379 - qJ(3) * t13 - t291 * t32) + t297 * (-pkin(2) * t13 - pkin(3) * t36) - pkin(1) * t13 + pkin(7) * (t14 * t297 + t148 * t294), t294 * (t160 * t292 - t415) + t253, t493, t431, t326, -t492, t327, t294 * (-t12 * t291 + t15 * t292 - t473) + t297 * (-t332 + t472) + t464, t294 * (-t29 * t291 + t292 * t34 - t489) + t297 * (-t305 + t483) - t496, t294 * (-t291 * t31 + t292 * t35 + t488) + t297 * (-t309 + t497) + t490 + pkin(7) * (-t140 * t294 + t486), t294 * (-qJ(3) * t8 + t11 * t292 - t291 * t6) + t297 * (-pkin(2) * t8 - pkin(3) * t24 - t325) - pkin(1) * t8 + pkin(7) * (t294 * t54 + t297 * t9), t394, t431, -t493, t327, t492, t326, t294 * (-t291 * t5 + t292 * t7 - t473) + t297 * (-t398 + t472) + t464, t294 * (-t16 * t291 + t21 * t292 + t488) + t297 * (-t308 - t397 + t497) + t495, t294 * (-t19 * t291 + t22 * t292 + t489) + t297 * (t38 - t396 - t483) + t496, t294 * (-qJ(3) * t3 - t1 * t291 + t2 * t292) + t297 * (-pkin(2) * t3 - pkin(3) * t17 - t399) - pkin(1) * t3 + pkin(7) * (t294 * t40 + t297 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, t284 - t285, t339, t280, t338, qJDD(2), -t242, -t243, 0, 0, t245 * t291 - t251 * t292, -t220 * t291 + t222 * t292, t248 * t292 + t414, t250 * t291 - t292 * t307, t247 * t291 - t351, (-t264 * t291 + t265 * t292) * t341, -pkin(2) * t220 + qJ(3) * t182 - t361, -pkin(2) * t222 + qJ(3) * t192 + t362, -pkin(2) * t213 + qJ(3) * t167 + t91, -pkin(2) * t201 + qJ(3) * t91, t76, t471, -t476, t75, t467, t395, t291 * t66 + t292 * t53 + t482, t291 * t77 + t292 * t58 - t481, t23 * t292 + t26 * t291 + t470, -pkin(2) * t148 - pkin(8) * t380 + qJ(3) * t14 + t292 * t32, t160 * t291 + t413, t476, -t467, t76, t471, t75, t12 * t292 + t15 * t291 + t470, t29 * t292 + t291 * t34 - t482, pkin(2) * t140 + t291 * t35 + t292 * t31 + t487, -pkin(2) * t54 + qJ(3) * t9 + t11 * t291 + t292 * t6, t395, -t467, -t476, t75, -t471, t76, t291 * t7 + t292 * t5 + t470, t16 * t292 + t21 * t291 + t481, t19 * t292 + t22 * t291 + t482, -pkin(2) * t40 + qJ(3) * t4 + t1 * t292 + t2 * t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t222, t213, t201, 0, 0, 0, 0, 0, 0, t409, t407, t402, t148, 0, 0, 0, 0, 0, 0, t402, -t409, -t140, t54, 0, 0, 0, 0, 0, 0, t402, -t407, t409, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t358, -t186, t404, -t358, -t408, -t266, -t56, -t57, 0, 0, -t266, -t404, t408, t358, -t186, -t358, t332, t305, t275 + t309, t325, -t266, t408, t404, -t358, t186, t358, t398, t39 + t397, pkin(5) * t359 - t304 + t396, t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, -t310, t403, t46, 0, 0, 0, 0, 0, 0, t404, t403, t310, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t408, -t311, t405, t39;];
tauJ_reg  = t47;
