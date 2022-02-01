% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRP1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:20:03
% EndTime: 2022-01-20 10:20:11
% DurationCPUTime: 5.66s
% Computational Cost: add. (18067->384), mult. (25982->523), div. (0->0), fcn. (14576->8), ass. (0->279)
t424 = qJD(1) + qJD(2);
t422 = t424 ^ 2;
t423 = qJDD(1) + qJDD(2);
t428 = sin(pkin(8));
t429 = cos(pkin(8));
t385 = t429 * t422 + t428 * t423;
t388 = t428 * t422 - t429 * t423;
t431 = sin(qJ(2));
t434 = cos(qJ(2));
t326 = t434 * t385 - t431 * t388;
t427 = g(3) - qJDD(3);
t364 = qJ(3) * t385 - t429 * t427;
t508 = qJ(3) * t388 - t428 * t427;
t274 = pkin(6) * t326 + t434 * t364 - t431 * t508;
t330 = t431 * t385 + t434 * t388;
t432 = sin(qJ(1));
t435 = cos(qJ(1));
t288 = t432 * t326 + t435 * t330;
t516 = pkin(6) * t330 + t431 * t364 + t434 * t508;
t525 = pkin(5) * t288 + t432 * t274 + t435 * t516;
t507 = t435 * t326 - t432 * t330;
t524 = pkin(5) * t507 + t435 * t274 - t432 * t516;
t410 = t435 * g(1) + t432 * g(2);
t437 = qJD(1) ^ 2;
t440 = -t437 * pkin(1) - t410;
t409 = t432 * g(1) - t435 * g(2);
t442 = qJDD(1) * pkin(1) + t409;
t342 = t431 * t442 + t434 * t440;
t334 = -t422 * pkin(2) + t342;
t439 = -t431 * t440 + t434 * t442;
t438 = t423 * pkin(2) + t439;
t292 = t428 * t334 - t429 * t438;
t293 = t429 * t334 + t428 * t438;
t453 = t428 * t292 + t429 * t293;
t237 = t429 * t292 - t428 * t293;
t469 = t434 * t237;
t196 = -t431 * t453 + t469;
t475 = t431 * t237;
t510 = t434 * t453 + t475;
t178 = t432 * t196 + t435 * t510;
t521 = t435 * t196 - t432 * t510;
t392 = t434 * t422 + t431 * t423;
t395 = t431 * t422 - t434 * t423;
t337 = t432 * t392 + t435 * t395;
t368 = pkin(6) * t392 - t434 * g(3);
t509 = pkin(6) * t395 - t431 * g(3);
t518 = pkin(5) * t337 + t432 * t368 + t435 * t509;
t443 = t435 * t392 - t432 * t395;
t517 = pkin(5) * t443 + t435 * t368 - t432 * t509;
t452 = t434 * t342 - t431 * t439;
t297 = -t431 * t342 - t434 * t439;
t468 = t435 * t297;
t511 = -t432 * t452 + t468;
t474 = t432 * t297;
t246 = t435 * t452 + t474;
t430 = sin(qJ(4));
t433 = cos(qJ(4));
t408 = t433 * t422 * t430;
t464 = qJDD(4) + t408;
t506 = t464 * pkin(4);
t282 = -t422 * pkin(3) + t423 * pkin(7) + t293;
t265 = t430 * t282 + t433 * t427;
t465 = qJD(4) * t433;
t458 = t424 * t465;
t476 = t430 * t423;
t380 = t458 + t476;
t371 = t380 * qJ(5);
t499 = -t371 - t265 + t506;
t436 = qJD(4) ^ 2;
t425 = t430 ^ 2;
t483 = t425 * t422;
t405 = -t436 - t483;
t401 = qJDD(4) - t408;
t477 = t430 * t401;
t350 = t433 * t405 - t477;
t498 = pkin(3) * t350;
t426 = t433 ^ 2;
t415 = t426 * t422;
t407 = -t415 - t436;
t478 = t430 * t464;
t352 = t433 * t407 - t478;
t466 = qJD(4) * t424;
t459 = t430 * t466;
t470 = t433 * t423;
t382 = -0.2e1 * t459 + t470;
t307 = t428 * t352 + t429 * t382;
t309 = t429 * t352 - t428 * t382;
t260 = t434 * t307 + t431 * t309;
t262 = -t431 * t307 + t434 * t309;
t212 = t435 * t260 + t432 * t262;
t497 = pkin(5) * t212;
t471 = t433 * t401;
t354 = -t430 * t405 - t471;
t379 = 0.2e1 * t458 + t476;
t308 = t428 * t354 - t429 * t379;
t310 = t429 * t354 + t428 * t379;
t261 = t434 * t308 + t431 * t310;
t263 = -t431 * t308 + t434 * t310;
t213 = t435 * t261 + t432 * t263;
t496 = pkin(5) * t213;
t467 = t425 + t426;
t390 = t467 * t423;
t396 = t415 + t483;
t332 = t428 * t390 + t429 * t396;
t333 = t429 * t390 - t428 * t396;
t290 = t434 * t332 + t431 * t333;
t291 = -t431 * t332 + t434 * t333;
t236 = t435 * t290 + t432 * t291;
t495 = pkin(5) * t236;
t494 = pkin(6) * t260;
t493 = pkin(6) * t261;
t492 = pkin(6) * t290;
t391 = t433 * t464;
t348 = t430 * t407 + t391;
t491 = pkin(7) * t348;
t490 = pkin(7) * t350;
t487 = qJ(3) * t307;
t486 = qJ(3) * t308;
t485 = qJ(3) * t332;
t484 = t424 * t430;
t242 = (qJ(5) * t465 - 0.2e1 * qJD(5) * t430) * t424 + t499;
t480 = t430 * t242;
t281 = -t423 * pkin(3) - t422 * pkin(7) + t292;
t479 = t430 * t281;
t473 = t433 * t242;
t472 = t433 * t281;
t266 = t433 * t282 - t430 * t427;
t463 = 0.2e1 * qJD(5) * t424;
t462 = t428 * t476;
t461 = t429 * t476;
t457 = -pkin(1) * t348 + pkin(6) * t262;
t456 = -pkin(1) * t350 + pkin(6) * t263;
t455 = -pkin(2) * t348 + qJ(3) * t309;
t454 = -pkin(2) * t350 + qJ(3) * t310;
t220 = t430 * t265 + t433 * t266;
t360 = -t432 * t409 - t435 * t410;
t450 = t428 * t408;
t449 = t429 * t408;
t448 = pkin(1) * t260 + pkin(2) * t307 + pkin(3) * t382 + pkin(7) * t352;
t447 = pkin(1) * t261 + pkin(2) * t308 - pkin(3) * t379 + pkin(7) * t354;
t446 = pkin(1) * t290 + pkin(2) * t332 + pkin(3) * t396 + pkin(7) * t390;
t248 = -pkin(3) * t348 + t265;
t403 = t435 * qJDD(1) - t432 * t437;
t445 = -pkin(5) * t403 - t432 * g(3);
t219 = t433 * t265 - t430 * t266;
t359 = t435 * t409 - t432 * t410;
t381 = -t459 + t470;
t399 = qJD(4) * pkin(4) - qJ(5) * t484;
t441 = t381 * qJ(5) - qJD(4) * t399 + t433 * t463 + t266;
t247 = -t381 * pkin(4) - qJ(5) * t415 + t399 * t484 + qJDD(5) + t281;
t412 = t430 * t463;
t406 = t415 - t436;
t404 = t436 - t483;
t402 = t432 * qJDD(1) + t435 * t437;
t397 = t415 - t483;
t377 = -pkin(5) * t402 + t435 * g(3);
t376 = t467 * t466;
t358 = t428 * qJDD(4) + t429 * t376;
t357 = -t429 * qJDD(4) + t428 * t376;
t356 = t433 * t380 - t425 * t466;
t355 = -t430 * t381 - t426 * t466;
t353 = -t430 * t404 + t391;
t351 = t433 * t406 - t477;
t349 = t433 * t404 + t478;
t347 = t430 * t406 + t471;
t344 = (t380 + t458) * t430;
t343 = (t381 - t459) * t433;
t340 = -pkin(4) * t379 - qJ(5) * t401;
t324 = qJ(3) * t333;
t323 = -t430 * t379 + t433 * t382;
t322 = t433 * t379 + t430 * t382;
t318 = t429 * t353 + t462;
t317 = t429 * t351 + t428 * t470;
t316 = t428 * t353 - t461;
t315 = t428 * t351 - t429 * t470;
t314 = t429 * t356 - t450;
t313 = t429 * t355 + t450;
t312 = t428 * t356 + t449;
t311 = t428 * t355 - t449;
t302 = -t431 * t357 + t434 * t358;
t301 = t434 * t357 + t431 * t358;
t300 = t429 * t323 - t428 * t397;
t299 = t428 * t323 + t429 * t397;
t294 = pkin(1) * g(3) + pkin(6) * t452;
t284 = pkin(6) * t291;
t280 = -t431 * t316 + t434 * t318;
t279 = -t431 * t315 + t434 * t317;
t278 = t434 * t316 + t431 * t318;
t277 = t434 * t315 + t431 * t317;
t270 = -t431 * t312 + t434 * t314;
t269 = -t431 * t311 + t434 * t313;
t268 = t434 * t312 + t431 * t314;
t267 = t434 * t311 + t431 * t313;
t255 = t472 - t490;
t254 = t479 - t491;
t253 = -t432 * t301 + t435 * t302;
t252 = t435 * t301 + t432 * t302;
t251 = -t431 * t299 + t434 * t300;
t250 = t434 * t299 + t431 * t300;
t249 = t266 - t498;
t244 = -qJ(5) * t405 + t247;
t243 = -pkin(4) * t415 + t441;
t241 = t412 + (-t458 + t476) * qJ(5) - t499;
t240 = pkin(4) * t382 + qJ(5) * t407 - t247;
t239 = -t432 * t290 + t435 * t291;
t233 = pkin(5) * t239;
t232 = qJ(5) * t470 + (t396 - t415) * pkin(4) + t441;
t231 = pkin(2) * t427 + qJ(3) * t453;
t230 = -t498 + (-t405 - t415) * pkin(4) + t441;
t229 = -t432 * t278 + t435 * t280;
t228 = -t432 * t277 + t435 * t279;
t227 = t435 * t278 + t432 * t280;
t226 = t435 * t277 + t432 * t279;
t225 = -qJ(5) * t458 + t248 + t371 + t412 - 0.2e1 * t506;
t224 = -t432 * t268 + t435 * t270;
t223 = -t432 * t267 + t435 * t269;
t222 = t435 * t268 + t432 * t270;
t221 = t435 * t267 + t432 * t269;
t217 = -qJ(5) * t391 - t430 * t240 - t491;
t216 = t433 * t244 - t430 * t340 - t490;
t215 = -t432 * t261 + t435 * t263;
t214 = -t432 * t260 + t435 * t262;
t211 = pkin(5) * t215;
t210 = pkin(5) * t214;
t209 = t429 * t219 - t485;
t208 = t428 * t219 + t324;
t207 = -t432 * t250 + t435 * t251;
t206 = t435 * t250 + t432 * t251;
t205 = -pkin(4) * t247 + qJ(5) * t243;
t204 = t433 * t243 - t480;
t203 = t430 * t243 + t473;
t202 = t429 * t220 + t428 * t281;
t201 = t428 * t220 - t429 * t281;
t200 = -t428 * t249 + t429 * t255 - t486;
t199 = -t428 * t248 + t429 * t254 - t487;
t198 = -t430 * t232 + t433 * t241;
t193 = t429 * t249 + t428 * t255 + t454;
t192 = t429 * t248 + t428 * t254 + t455;
t191 = -pkin(4) * t462 + t429 * t198 - t485;
t190 = pkin(4) * t461 + t428 * t198 + t324;
t189 = t429 * t204 + t428 * t247;
t188 = t428 * t204 - t429 * t247;
t187 = -pkin(3) * t203 - pkin(4) * t242;
t186 = t429 * t216 - t428 * t230 - t486;
t185 = t429 * t217 - t428 * t225 - t487;
t184 = t428 * t216 + t429 * t230 + t454;
t183 = t428 * t217 + t429 * t225 + t455;
t182 = -t431 * t208 + t434 * t209 - t492;
t181 = t434 * t208 + t431 * t209 + t284;
t180 = -t431 * t201 + t434 * t202;
t179 = t434 * t201 + t431 * t202;
t176 = pkin(6) * t196 + qJ(3) * t469 - t431 * t231;
t175 = pkin(1) * t427 + pkin(6) * t510 + qJ(3) * t475 + t434 * t231;
t174 = -qJ(3) * t201 - (pkin(3) * t428 - pkin(7) * t429) * t219;
t173 = -pkin(7) * t203 - qJ(5) * t473 - t430 * t205;
t172 = -t431 * t193 + t434 * t200 - t493;
t171 = -t431 * t192 + t434 * t199 - t494;
t170 = t434 * t193 + t431 * t200 + t456;
t169 = t434 * t192 + t431 * t199 + t457;
t168 = -t431 * t190 + t434 * t191 - t492;
t167 = t434 * t190 + t431 * t191 + t284;
t166 = -t431 * t188 + t434 * t189;
t165 = t434 * t188 + t431 * t189;
t164 = qJ(3) * t202 - (-pkin(3) * t429 - pkin(7) * t428 - pkin(2)) * t219;
t163 = -t431 * t184 + t434 * t186 - t493;
t162 = -t431 * t183 + t434 * t185 - t494;
t161 = t434 * t184 + t431 * t186 + t456;
t160 = t434 * t183 + t431 * t185 + t457;
t159 = -t432 * t179 + t435 * t180;
t158 = t435 * t179 + t432 * t180;
t157 = -qJ(3) * t188 + t429 * t173 - t428 * t187;
t156 = -t432 * t165 + t435 * t166;
t155 = t435 * t165 + t432 * t166;
t154 = -pkin(2) * t203 + qJ(3) * t189 + t428 * t173 + t429 * t187;
t153 = -pkin(6) * t179 - t431 * t164 + t434 * t174;
t152 = pkin(1) * t219 + pkin(6) * t180 + t434 * t164 + t431 * t174;
t151 = -pkin(6) * t165 - t431 * t154 + t434 * t157;
t150 = -pkin(1) * t203 + pkin(6) * t166 + t434 * t154 + t431 * t157;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t402, -t403, 0, t360, 0, 0, 0, 0, 0, 0, -t443, t337, 0, t246, 0, 0, 0, 0, 0, 0, -t507, t288, 0, t178, 0, 0, 0, 0, 0, 0, t214, t215, t239, t159, 0, 0, 0, 0, 0, 0, t214, t215, t239, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t403, -t402, 0, t359, 0, 0, 0, 0, 0, 0, -t337, -t443, 0, -t511, 0, 0, 0, 0, 0, 0, -t288, -t507, 0, -t521, 0, 0, 0, 0, 0, 0, t212, t213, t236, t158, 0, 0, 0, 0, 0, 0, t212, t213, t236, t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t427, 0, 0, 0, 0, 0, 0, t348, t350, 0, -t219, 0, 0, 0, 0, 0, 0, t348, t350, 0, t203; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t403, 0, -t402, 0, t445, -t377, -t359, -pkin(5) * t359, 0, 0, -t337, 0, -t443, 0, t518, t517, t511, pkin(5) * t511 + pkin(6) * t468 - t432 * t294, 0, 0, -t288, 0, -t507, 0, t525, t524, t521, pkin(5) * t521 - t432 * t175 + t435 * t176, t224, t207, t229, t223, t228, t253, -t432 * t169 + t435 * t171 - t497, -t432 * t170 + t435 * t172 - t496, -t432 * t181 + t435 * t182 - t495, -pkin(5) * t158 - t432 * t152 + t435 * t153, t224, t207, t229, t223, t228, t253, -t432 * t160 + t435 * t162 - t497, -t432 * t161 + t435 * t163 - t496, -t432 * t167 + t435 * t168 - t495, -pkin(5) * t155 - t432 * t150 + t435 * t151; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t402, 0, t403, 0, t377, t445, t360, pkin(5) * t360, 0, 0, t443, 0, -t337, 0, -t517, t518, t246, pkin(5) * t246 + pkin(6) * t474 + t435 * t294, 0, 0, t507, 0, -t288, 0, -t524, t525, t178, pkin(5) * t178 + t435 * t175 + t432 * t176, t222, t206, t227, t221, t226, t252, t435 * t169 + t432 * t171 + t210, t435 * t170 + t432 * t172 + t211, t435 * t181 + t432 * t182 + t233, pkin(5) * t159 + t435 * t152 + t432 * t153, t222, t206, t227, t221, t226, t252, t435 * t160 + t432 * t162 + t210, t435 * t161 + t432 * t163 + t211, t435 * t167 + t432 * t168 + t233, pkin(5) * t156 + t435 * t150 + t432 * t151; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t409, t410, 0, 0, 0, 0, 0, 0, 0, t423, -pkin(1) * t395 + t439, -pkin(1) * t392 - t342, 0, -pkin(1) * t297, 0, 0, 0, 0, 0, t423, -pkin(1) * t330 - pkin(2) * t388 - t292, -pkin(1) * t326 - pkin(2) * t385 - t293, 0, -pkin(1) * t196 - pkin(2) * t237, t344, t322, t349, t343, t347, 0, t448 - t472, t447 + t479, t220 + t446, pkin(1) * t179 + pkin(2) * t201 - pkin(3) * t281 + pkin(7) * t220, t344, t322, t349, t343, t347, 0, -qJ(5) * t478 + t433 * t240 + t448, t430 * t244 + t433 * t340 + t447, t433 * t232 + t430 * t241 + t446, pkin(1) * t165 + pkin(2) * t188 - pkin(3) * t247 + pkin(7) * t204 - qJ(5) * t480 + t433 * t205;];
tauB_reg = t1;
