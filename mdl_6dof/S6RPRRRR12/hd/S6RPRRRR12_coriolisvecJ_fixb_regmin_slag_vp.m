% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:00:42
% EndTime: 2019-03-09 08:01:15
% DurationCPUTime: 12.14s
% Computational Cost: add. (33521->633), mult. (115882->951), div. (0->0), fcn. (102809->16), ass. (0->305)
t294 = cos(qJ(3));
t285 = cos(pkin(7));
t282 = sin(pkin(6));
t283 = cos(pkin(14));
t430 = t282 * t283;
t281 = sin(pkin(7));
t286 = cos(pkin(6));
t432 = t281 * t286;
t323 = t285 * t430 + t432;
t313 = t323 * qJD(1);
t290 = sin(qJ(3));
t279 = sin(pkin(14));
t405 = qJD(1) * t282;
t381 = t279 * t405;
t359 = t290 * t381;
t211 = t294 * t313 - t359;
t424 = t285 * t290;
t219 = t282 * (t279 * t294 + t283 * t424) + t290 * t432;
t214 = t219 * qJD(1);
t293 = cos(qJ(4));
t284 = cos(pkin(8));
t289 = sin(qJ(4));
t427 = t284 * t289;
t170 = t211 * t293 - t214 * t427;
t280 = sin(pkin(8));
t399 = qJD(4) * t293;
t375 = t280 * t399;
t479 = -t170 + t375;
t380 = t283 * t405;
t404 = qJD(1) * t286;
t390 = pkin(1) * t404;
t246 = qJ(2) * t380 + t279 * t390;
t205 = pkin(10) * t313 + t246;
t267 = t283 * t390;
t439 = t279 * t282;
t310 = pkin(2) * t286 + (-pkin(10) * t285 - qJ(2)) * t439;
t209 = qJD(1) * t310 + t267;
t235 = (-pkin(10) * t279 * t281 - pkin(2) * t283 - pkin(1)) * t282;
t226 = qJD(1) * t235 + qJD(2);
t423 = t285 * t294;
t431 = t281 * t294;
t412 = t209 * t423 + t226 * t431;
t353 = -t290 * t205 + t412;
t460 = pkin(11) * t284;
t126 = -t214 * t460 + t353;
t305 = -t290 * (t209 * t285 + t226 * t281) - t294 * t205;
t428 = t284 * t211;
t127 = -pkin(11) * t428 + t305;
t206 = t280 * t211;
t172 = pkin(3) * t214 - pkin(11) * t206;
t435 = t280 * t289;
t271 = pkin(11) * t435;
t426 = t284 * t293;
t470 = pkin(3) * t426 - t271;
t478 = t470 * qJD(4) - t293 * t126 - t127 * t427 - t172 * t435;
t419 = t293 * t214;
t169 = t211 * t289 + t284 * t419;
t78 = -t127 * t280 + t284 * t172;
t477 = pkin(4) * t169 - pkin(12) * t170 - (pkin(4) * t289 - pkin(12) * t293) * t280 * qJD(4) + t78;
t441 = t214 * t280;
t476 = -pkin(12) * t441 + t478;
t400 = qJD(4) * t289;
t376 = t280 * t400;
t475 = t169 - t376;
t429 = t283 * t294;
t233 = (-t279 * t424 + t429) * t405;
t315 = t282 * (-t279 * t423 - t283 * t290);
t232 = qJD(1) * t315;
t362 = t281 * t381;
t312 = t232 * t284 + t280 * t362;
t418 = t293 * t294;
t422 = t289 * t290;
t321 = t284 * t418 - t422;
t474 = -t233 * t293 - t289 * t312 + t285 * t375 + (t321 * qJD(4) + (-t284 * t422 + t418) * qJD(3)) * t281;
t402 = qJD(3) * t290;
t473 = t284 * t362 + (-t281 * t402 - t232) * t280;
t288 = sin(qJ(5));
t292 = cos(qJ(5));
t251 = -t292 * t284 + t288 * t435;
t436 = t280 * t288;
t411 = -qJD(5) * t251 - t214 * t436 + t292 * t479;
t434 = t280 * t292;
t252 = t284 * t288 + t289 * t434;
t410 = qJD(5) * t252 + t214 * t434 + t288 * t479;
t433 = t280 * t293;
t407 = pkin(3) * t427 + pkin(11) * t433;
t472 = t407 * qJD(4) - t289 * t126 + t293 * (t127 * t284 + t172 * t280);
t258 = t281 * t380;
t352 = t285 * t404 - t258;
t328 = qJD(3) + t352;
t316 = t328 * t280;
t230 = t293 * t316;
t151 = -t211 * t426 + t214 * t289 - t230;
t150 = qJD(5) + t151;
t302 = t316 + t428;
t153 = t289 * t302 + t419;
t306 = t284 * t328 - t206;
t303 = -qJD(4) - t306;
t179 = t292 * t303;
t106 = t153 * t288 + t179;
t105 = qJD(6) + t106;
t449 = pkin(4) * t441 + t472;
t471 = t282 ^ 2 * (t279 ^ 2 + t283 ^ 2);
t420 = t290 * t293;
t421 = t289 * t294;
t320 = t284 * t421 + t420;
t413 = -t233 * t289 + t293 * t312 + t285 * t376 + (t320 * qJD(4) + (t284 * t420 + t421) * qJD(3)) * t281;
t243 = pkin(12) * t284 + t407;
t244 = (-pkin(4) * t293 - pkin(12) * t289 - pkin(3)) * t280;
t409 = t292 * t243 + t288 * t244;
t395 = qJD(5) * t292;
t397 = qJD(5) * t288;
t469 = t243 * t397 - t244 * t395 + t288 * t477 - t476 * t292;
t213 = t219 * qJD(3);
t108 = t292 * t153 - t288 * t303;
t102 = pkin(11) * t302 - t305;
t104 = pkin(3) * t328 + t126;
t177 = -t209 * t281 + t285 * t226;
t135 = -pkin(3) * t211 - pkin(11) * t441 + t177;
t51 = t293 * t102 + t104 * t427 + t135 * t435;
t46 = -pkin(12) * t303 + t51;
t69 = -t104 * t280 + t284 * t135;
t48 = pkin(4) * t151 - pkin(12) * t153 + t69;
t17 = t288 * t48 + t292 * t46;
t387 = t283 * t423;
t364 = t282 * t387;
t401 = qJD(3) * t294;
t378 = t281 * t401;
t203 = t378 * t404 + (qJD(1) * t364 - t359) * qJD(3);
t204 = qJD(1) * t213;
t392 = qJD(1) * qJD(2);
t372 = t282 * t392;
t358 = t279 * t372;
t332 = t281 * t358;
t461 = pkin(11) * t280;
t160 = pkin(3) * t204 - t203 * t461 + t332;
t373 = t284 * t399;
t331 = t285 * t358;
t403 = qJD(2) * t282;
t260 = t403 * t429;
t377 = t285 * t401;
t382 = qJD(1) * t260 + t209 * t377 + t226 * t378;
t298 = (-qJD(3) * t205 - t331) * t290 + t382;
t443 = t204 * t284;
t94 = -pkin(11) * t443 + t298;
t307 = qJD(2) * t315;
t304 = qJD(1) * t307;
t296 = qJD(3) * t305 + t304;
t95 = -t203 * t460 + t296;
t311 = t102 * t400 - t104 * t373 - t135 * t375 - t160 * t435 - t293 * t94 - t95 * t427;
t437 = t280 * t204;
t21 = pkin(12) * t437 - t311;
t70 = t284 * t160 - t280 * t95;
t90 = qJD(4) * t230 + t293 * t203 - t204 * t427 + t211 * t373 - t214 * t400;
t367 = t289 * t203 + t204 * t426;
t91 = qJD(4) * t153 + t367;
t42 = pkin(4) * t91 - pkin(12) * t90 + t70;
t6 = -qJD(5) * t17 - t21 * t288 + t292 * t42;
t4 = -pkin(5) * t91 - t6;
t468 = t105 * (pkin(5) * t108 + pkin(13) * t105) + t4;
t287 = sin(qJ(6));
t291 = cos(qJ(6));
t62 = -qJD(5) * t179 + t292 * t90 + (-qJD(5) * t153 + t437) * t288;
t73 = t108 * t291 + t150 * t287;
t32 = qJD(6) * t73 + t287 * t62 - t291 * t91;
t50 = -t289 * t102 + t293 * (t104 * t284 + t135 * t280);
t385 = t286 * t431;
t438 = t279 * t290;
t218 = t282 * t438 - t364 - t385;
t425 = t285 * t286;
t322 = t281 * t430 - t425;
t314 = t322 * t280;
t301 = -t218 * t284 - t314;
t462 = pkin(1) * t286;
t408 = qJ(2) * t430 + t279 * t462;
t216 = pkin(10) * t323 + t408;
t270 = t283 * t462;
t220 = t270 + t310;
t335 = t220 * t285 + t235 * t281;
t465 = -t294 * t216 - t290 * t335;
t125 = pkin(11) * t301 - t465;
t130 = -pkin(3) * t322 - t290 * t216 - t219 * t460 + t294 * t335;
t180 = -t220 * t281 + t285 * t235;
t142 = pkin(3) * t218 - t219 * t461 + t180;
t467 = -t289 * t125 + (t130 * t284 + t142 * t280) * t293;
t379 = t279 * t403;
t299 = t220 * t377 + t235 * t378 + t260 + (-qJD(3) * t216 - t285 * t379) * t290;
t114 = -t213 * t460 + t299;
t212 = (t385 + (t387 - t438) * t282) * qJD(3);
t297 = qJD(3) * t465 + t307;
t115 = -t212 * t460 + t297;
t361 = t281 * t379;
t165 = pkin(3) * t213 - t212 * t461 + t361;
t374 = t284 * t400;
t464 = (t115 * t284 + t165 * t280) * t293 - t289 * t114 - t125 * t399 - t130 * t374 - t142 * t376;
t63 = qJD(5) * t108 - t204 * t434 + t288 * t90;
t308 = t293 * t114 + t115 * t427 - t125 * t400 + t130 * t373 + t142 * t375 + t165 * t435;
t442 = t213 * t280;
t27 = pkin(12) * t442 + t308;
t184 = -t218 * t280 + t284 * t322;
t384 = t293 * t125 + t130 * t427 + t142 * t435;
t54 = -pkin(12) * t184 + t384;
t440 = t219 * t289;
t161 = t218 * t426 + t293 * t314 + t440;
t162 = t219 * t293 + t289 * t301;
t74 = -t130 * t280 + t284 * t142;
t57 = pkin(4) * t161 - pkin(12) * t162 + t74;
t342 = t288 * t57 + t292 * t54;
t109 = qJD(4) * t162 + t212 * t289 + t213 * t426;
t110 = -t213 * t427 + t212 * t293 + (t293 * t301 - t440) * qJD(4);
t75 = -t115 * t280 + t284 * t165;
t44 = pkin(4) * t109 - pkin(12) * t110 + t75;
t463 = -qJD(5) * t342 - t27 * t288 + t292 * t44;
t330 = t102 * t399 + t104 * t374 + t135 * t376 - t160 * t433 + t289 * t94 - t95 * t426;
t22 = -pkin(4) * t437 + t330;
t12 = pkin(5) * t63 - pkin(13) * t62 + t22;
t5 = t292 * t21 + t288 * t42 + t48 * t395 - t397 * t46;
t3 = pkin(13) * t91 + t5;
t15 = pkin(13) * t150 + t17;
t45 = pkin(4) * t303 - t50;
t25 = t106 * pkin(5) - t108 * pkin(13) + t45;
t346 = t15 * t287 - t25 * t291;
t1 = -qJD(6) * t346 + t12 * t287 + t291 * t3;
t82 = pkin(4) * t153 + pkin(12) * t151;
t459 = t288 * t82 + t292 * t50;
t71 = t108 * t287 - t291 * t150;
t457 = t105 * t71;
t456 = t105 * t73;
t393 = qJD(6) * t291;
t394 = qJD(6) * t287;
t31 = -t108 * t394 + t150 * t393 + t287 * t91 + t291 * t62;
t455 = t287 * t31;
t454 = t287 * t63;
t453 = t291 * t63;
t452 = t475 * pkin(5) + t409 * qJD(5) + t476 * t288 + t477 * t292;
t227 = t252 * t287 + t291 * t433;
t451 = -qJD(6) * t227 - t287 * t475 + t291 * t411;
t389 = t287 * t433;
t450 = -qJD(6) * t389 + t252 * t393 + t287 * t411 + t291 * t475;
t448 = -t51 + t150 * (pkin(5) * t288 - pkin(13) * t292);
t447 = t106 * t150;
t446 = t108 * t150;
t445 = t151 * t292;
t444 = t204 * t280 ^ 2;
t223 = t281 * t320 + t285 * t435;
t250 = -t280 * t431 + t284 * t285;
t185 = t223 * t288 - t250 * t292;
t415 = qJD(5) * t185 + t288 * t473 - t292 * t474;
t186 = t223 * t292 + t250 * t288;
t414 = qJD(5) * t186 + t288 * t474 + t292 * t473;
t398 = qJD(5) * t287;
t396 = qJD(5) * t291;
t295 = qJD(1) ^ 2;
t388 = t282 * t286 * t295;
t263 = -pkin(5) * t292 - pkin(13) * t288 - pkin(4);
t368 = pkin(13) * t153 - qJD(6) * t263 + t459;
t366 = t292 * t150;
t365 = t105 * t291;
t80 = t153 * t287 - t291 * t445;
t356 = t291 * t395 - t80;
t242 = t271 + (-pkin(3) * t293 - pkin(4)) * t284;
t187 = pkin(5) * t251 - pkin(13) * t252 + t242;
t355 = pkin(13) * t475 - qJD(6) * t187 + t469;
t189 = -pkin(13) * t433 + t409;
t354 = -pkin(5) * t410 + pkin(13) * t411 + qJD(6) * t189 - t449;
t222 = -t281 * t321 - t285 * t433;
t349 = -qJD(6) * t222 + t415;
t348 = qJD(6) * t186 - t413;
t347 = -0.2e1 * t286 * t372;
t10 = t15 * t291 + t25 * t287;
t19 = pkin(13) * t161 + t342;
t128 = t162 * t288 + t184 * t292;
t129 = t162 * t292 - t184 * t288;
t53 = pkin(4) * t184 - t467;
t35 = pkin(5) * t128 - pkin(13) * t129 + t53;
t345 = t19 * t291 + t287 * t35;
t344 = -t19 * t287 + t291 * t35;
t16 = -t288 * t46 + t292 * t48;
t341 = -t288 * t54 + t292 * t57;
t77 = t129 * t291 + t161 * t287;
t76 = t129 * t287 - t161 * t291;
t334 = -t243 * t288 + t244 * t292;
t333 = (-qJ(2) * t381 + t267) * t279 - t246 * t283;
t326 = -t105 * t393 - t454;
t325 = -t105 * t394 + t453;
t319 = t292 * t27 + t288 * t44 + t57 * t395 - t397 * t54;
t318 = -pkin(12) * t91 + t150 * t45;
t14 = -pkin(5) * t150 - t16;
t309 = -pkin(13) * t63 + (t14 + t16) * t105;
t2 = -qJD(6) * t10 + t291 * t12 - t287 * t3;
t300 = qJD(4) * t303;
t28 = -pkin(4) * t442 - t464;
t228 = t252 * t291 - t389;
t188 = pkin(5) * t433 - t334;
t79 = -t291 * t153 - t287 * t445;
t67 = -qJD(5) * t128 + t110 * t292 + t213 * t436;
t66 = qJD(5) * t129 + t110 * t288 - t213 * t434;
t37 = -qJD(6) * t76 + t109 * t287 + t291 * t67;
t36 = qJD(6) * t77 - t109 * t291 + t287 * t67;
t33 = -pkin(5) * t153 + t288 * t50 - t292 * t82;
t18 = -pkin(5) * t161 - t341;
t13 = pkin(5) * t66 - pkin(13) * t67 + t28;
t8 = -pkin(5) * t109 - t463;
t7 = pkin(13) * t109 + t319;
t9 = [0, 0, 0, t279 * t347, t283 * t347, 0.2e1 * t392 * t471 ((t283 * t408 + (qJ(2) * t439 - t270) * t279) * qJD(1) - t333) * t403, t203 * t219 + t212 * t214, -t203 * t218 - t204 * t219 + t211 * t212 - t213 * t214, -t203 * t322 + t212 * t328, t204 * t322 - t213 * t328, 0, t177 * t213 + t180 * t204 - t211 * t361 + t218 * t332 - t296 * t322 + t297 * t328, t177 * t212 + t180 * t203 + t214 * t361 + t219 * t332 + t298 * t322 - t299 * t328, t110 * t153 + t162 * t90, -t109 * t153 - t110 * t151 - t161 * t90 - t162 * t91, -t184 * t90 - t303 * t110 + (t153 * t213 + t162 * t204) * t280, t184 * t91 + t303 * t109 + (-t151 * t213 - t161 * t204) * t280 (-t204 * t184 - t213 * t303) * t280, t69 * t109 + t75 * t151 + t70 * t161 + t184 * t330 - t303 * t464 + t437 * t467 + t50 * t442 + t74 * t91, t75 * t153 + t74 * t90 + t70 * t162 + t69 * t110 + t308 * t303 - t311 * t184 + (-t204 * t384 - t51 * t213) * t280, t108 * t67 + t129 * t62, -t106 * t67 - t108 * t66 - t128 * t62 - t129 * t63, t108 * t109 + t129 * t91 + t150 * t67 + t161 * t62, -t106 * t109 - t128 * t91 - t150 * t66 - t161 * t63, t109 * t150 + t161 * t91, t28 * t106 + t16 * t109 + t22 * t128 + t150 * t463 + t6 * t161 + t341 * t91 + t45 * t66 + t53 * t63, t28 * t108 - t17 * t109 + t22 * t129 - t150 * t319 - t5 * t161 - t342 * t91 + t45 * t67 + t53 * t62, t31 * t77 + t37 * t73, -t31 * t76 - t32 * t77 - t36 * t73 - t37 * t71, t105 * t37 + t128 * t31 + t63 * t77 + t66 * t73, -t105 * t36 - t128 * t32 - t63 * t76 - t66 * t71, t105 * t66 + t128 * t63 (-qJD(6) * t345 + t13 * t291 - t287 * t7) * t105 + t344 * t63 + t2 * t128 - t346 * t66 + t8 * t71 + t18 * t32 + t4 * t76 + t14 * t36 -(qJD(6) * t344 + t13 * t287 + t291 * t7) * t105 - t345 * t63 - t1 * t128 - t10 * t66 + t8 * t73 + t18 * t31 + t4 * t77 + t14 * t37; 0, 0, 0, t279 * t388, t283 * t388, -t295 * t471, t333 * t405, 0, 0, 0, 0, 0, t285 * t204 - t232 * t328 + (t211 * t381 - t328 * t402) * t281, t285 * t203 + t233 * t328 + (-t214 * t381 - t328 * t401) * t281, 0, 0, 0, 0, 0, -t151 * t473 - t222 * t437 + t250 * t91 + t303 * t413, -t153 * t473 - t223 * t437 + t250 * t90 + t303 * t474, 0, 0, 0, 0, 0, t106 * t413 - t150 * t414 - t185 * t91 + t222 * t63, t108 * t413 + t150 * t415 - t186 * t91 + t222 * t62, 0, 0, 0, 0, 0 (-t186 * t287 + t222 * t291) * t63 + t185 * t32 + t414 * t71 + (t287 * t349 - t291 * t348) * t105 -(t186 * t291 + t222 * t287) * t63 + t185 * t31 + t414 * t73 + (t287 * t348 + t291 * t349) * t105; 0, 0, 0, 0, 0, 0, 0, -t214 * t211, -t211 ^ 2 + t214 ^ 2, -t211 * t328 + t203 (qJD(3) - t258) * t214 + (t214 * t425 - t213) * qJD(1), 0, -t177 * t214 - t305 * t352 + t304, qJD(3) * t412 - t177 * t211 + t290 * t331 + t352 * t353 - t382, t153 * t479 + t90 * t435, t151 * t170 + t153 * t169 + (-t289 * t91 + t293 * t90 + (-t151 * t293 - t153 * t289) * qJD(4)) * t280, t284 * t90 + t289 * t444 + t303 * t170 + (-t214 * t153 - t293 * t300) * t280, -t284 * t91 + t293 * t444 - t303 * t169 + (t214 * t151 + t289 * t300) * t280 (t214 * t303 + t443) * t280, -t280 * pkin(3) * t91 - t78 * t151 - t330 * t284 + t303 * t472 - t70 * t433 + t470 * t437 - t50 * t441 - t475 * t69, t311 * t284 - t78 * t153 - t69 * t170 + (-pkin(3) * t90 - t204 * t407 + t51 * t214 + t70 * t289 + t399 * t69) * t280 + t478 * t303, t108 * t411 + t252 * t62, -t106 * t411 - t108 * t410 - t251 * t62 - t252 * t63, -t108 * t475 + t150 * t411 + t252 * t91 - t433 * t62, t106 * t475 - t150 * t410 - t251 * t91 + t433 * t63, -t150 * t475 - t433 * t91, t334 * t91 - t6 * t433 + t242 * t63 + t22 * t251 + t410 * t45 - t475 * t16 + ((-qJD(5) * t243 - t477) * t292 + (-qJD(5) * t244 - t476) * t288) * t150 + t449 * t106, t449 * t108 + t150 * t469 + t17 * t475 + t22 * t252 + t242 * t62 - t409 * t91 + t411 * t45 + t5 * t433, t228 * t31 + t451 * t73, -t227 * t31 - t228 * t32 - t450 * t73 - t451 * t71, t105 * t451 + t228 * t63 + t251 * t31 + t410 * t73, -t105 * t450 - t227 * t63 - t251 * t32 - t410 * t71, t105 * t410 + t251 * t63 (t187 * t291 - t189 * t287) * t63 + t2 * t251 + t188 * t32 + t4 * t227 - t410 * t346 + t452 * t71 + t450 * t14 + (t287 * t355 - t291 * t354) * t105 -(t187 * t287 + t189 * t291) * t63 - t1 * t251 + t188 * t31 + t4 * t228 + t452 * t73 + t451 * t14 + (t287 * t354 + t291 * t355) * t105 - t410 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151 * t153, -t151 ^ 2 + t153 ^ 2, -t151 * t303 + t90, t306 * t153 - t367, t437, -t69 * t153 - t303 * t51 - t330, t69 * t151 - t303 * t50 + t311, t108 * t366 + t288 * t62 (t62 - t447) * t292 + (-t63 - t446) * t288, -t108 * t153 + t150 * t366 + t288 * t91, -t150 ^ 2 * t288 + t106 * t153 + t292 * t91, -t150 * t153, -pkin(4) * t63 - t51 * t106 - t16 * t153 + (-t22 + (-pkin(12) * qJD(5) - t82) * t150) * t292 + (t50 * t150 + t318) * t288, -pkin(4) * t62 - t51 * t108 + t17 * t153 + t22 * t288 + (pkin(12) * t397 + t459) * t150 + t318 * t292, t288 * t291 * t31 + (-t288 * t394 + t356) * t73, t71 * t80 + t73 * t79 + (-t287 * t73 - t291 * t71) * t395 + (-t455 - t291 * t32 + (t287 * t71 - t291 * t73) * qJD(6)) * t288, -t292 * t31 + t356 * t105 + (t150 * t73 + t325) * t288, t292 * t32 + (-t287 * t395 + t79) * t105 + (-t150 * t71 + t326) * t288, t105 * t150 * t288 - t292 * t63, t263 * t453 - t14 * t79 - t33 * t71 + (t287 * t368 + t291 * t448) * t105 + (t14 * t398 - t2 + (qJD(5) * t71 + t326) * pkin(12)) * t292 + (t14 * t393 + t4 * t287 - t150 * t346 + (t105 * t398 + t32) * pkin(12)) * t288, -t263 * t454 - t14 * t80 - t33 * t73 + (-t287 * t448 + t291 * t368) * t105 + (t14 * t396 + t1 + (qJD(5) * t73 - t325) * pkin(12)) * t292 + (-t14 * t394 + t4 * t291 - t150 * t10 + (t105 * t396 + t31) * pkin(12)) * t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t106, -t106 ^ 2 + t108 ^ 2, t62 + t447, t446 - t63, t91, -t108 * t45 + t150 * t17 + t6, t106 * t45 + t150 * t16 - t5, t365 * t73 + t455 (t31 - t457) * t291 + (-t32 - t456) * t287, t105 * t365 - t108 * t73 + t454, -t105 ^ 2 * t287 + t108 * t71 + t453, -t105 * t108, -pkin(5) * t32 + t108 * t346 - t17 * t71 + t309 * t287 - t291 * t468, -pkin(5) * t31 + t10 * t108 - t17 * t73 + t287 * t468 + t309 * t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t31 + t457, -t32 + t456, t63, t10 * t105 - t14 * t73 + t2, -t105 * t346 + t14 * t71 - t1;];
tauc_reg  = t9;
