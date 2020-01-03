% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPP3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:45
% EndTime: 2019-12-31 20:53:53
% DurationCPUTime: 5.07s
% Computational Cost: add. (10325->391), mult. (14758->470), div. (0->0), fcn. (7402->6), ass. (0->260)
t456 = qJD(3) ^ 2;
t442 = qJD(1) + qJD(2);
t440 = t442 ^ 2;
t450 = sin(qJ(3));
t446 = t450 ^ 2;
t518 = t446 * t440;
t417 = -t456 - t518;
t453 = cos(qJ(3));
t485 = t453 * t440 * t450;
t412 = qJDD(3) - t485;
t499 = t453 * t412;
t362 = t417 * t450 + t499;
t490 = qJD(3) * t442;
t478 = t453 * t490;
t441 = qJDD(1) + qJDD(2);
t509 = t450 * t441;
t390 = 0.2e1 * t478 + t509;
t451 = sin(qJ(2));
t454 = cos(qJ(2));
t308 = t362 * t451 + t390 * t454;
t312 = t362 * t454 - t390 * t451;
t452 = sin(qJ(1));
t455 = cos(qJ(1));
t261 = t308 * t455 + t312 * t452;
t537 = pkin(5) * t261;
t266 = t308 * t452 - t312 * t455;
t258 = pkin(5) * t266;
t447 = t453 ^ 2;
t517 = t447 * t440;
t419 = -t456 - t517;
t411 = qJDD(3) + t485;
t515 = t450 * t411;
t359 = -t419 * t453 + t515;
t428 = t450 * t490;
t498 = t453 * t441;
t392 = -0.2e1 * t428 + t498;
t307 = t359 * t451 - t392 * t454;
t311 = t359 * t454 + t392 * t451;
t260 = t307 * t455 + t311 * t452;
t538 = pkin(5) * t260;
t263 = t307 * t452 - t311 * t455;
t257 = pkin(5) * t263;
t534 = pkin(6) * t308;
t481 = pkin(1) * t308 + pkin(2) * t390 + pkin(7) * t362;
t514 = t450 * t412;
t355 = -t417 * t453 + t514;
t476 = -pkin(1) * t355 + pkin(6) * t312;
t416 = t456 - t518;
t500 = t453 * t411;
t361 = -t416 * t450 + t500;
t494 = t454 * t441;
t324 = t361 * t451 - t450 * t494;
t505 = t451 * t441;
t328 = t361 * t454 + t450 * t505;
t278 = t324 * t455 + t328 * t452;
t562 = t324 * t452 - t328 * t455;
t535 = pkin(6) * t307;
t351 = t419 * t450 + t500;
t540 = pkin(2) * t351;
t561 = qJ(4) * t419 + t540;
t560 = pkin(1) * t307 + pkin(7) * t359;
t477 = -pkin(1) * t351 - pkin(6) * t311;
t348 = t450 * t390;
t501 = t453 * t392;
t334 = -t501 + t348;
t407 = (t446 - t447) * t440;
t296 = t334 * t451 + t407 * t454;
t298 = t334 * t454 - t407 * t451;
t559 = t296 * t455 + t298 * t452;
t252 = t296 * t452 - t298 * t455;
t418 = -t456 + t517;
t360 = -t418 * t453 + t514;
t323 = t360 * t451 + t453 * t494;
t327 = t360 * t454 - t451 * t498;
t277 = t323 * t455 + t327 * t452;
t558 = t323 * t452 - t327 * t455;
t402 = t440 * t454 + t505;
t405 = t440 * t451 - t494;
t341 = t402 * t452 + t405 * t455;
t376 = pkin(6) * t402 - g(3) * t454;
t547 = pkin(6) * t405 - g(3) * t451;
t557 = pkin(5) * t341 + t376 * t452 + t455 * t547;
t465 = t402 * t455 - t405 * t452;
t556 = pkin(5) * t465 + t376 * t455 - t452 * t547;
t532 = pkin(7) * t351;
t422 = g(1) * t455 + g(2) * t452;
t457 = qJD(1) ^ 2;
t409 = -pkin(1) * t457 - t422;
t421 = g(1) * t452 - g(2) * t455;
t463 = qJDD(1) * pkin(1) + t421;
t345 = t451 * t409 - t454 * t463;
t346 = t409 * t454 + t451 * t463;
t474 = t345 * t451 + t346 * t454;
t293 = t345 * t454 - t346 * t451;
t493 = t455 * t293;
t552 = -t452 * t474 + t493;
t504 = t452 * t293;
t248 = t455 * t474 + t504;
t549 = pkin(2) * t355;
t531 = pkin(7) * t355;
t529 = t453 * g(3);
t464 = -qJDD(3) * pkin(3) - qJ(4) * t456 + qJDD(4) + t529;
t331 = -pkin(2) * t440 + pkin(7) * t441 + t346;
t523 = qJ(4) * t450;
t469 = -pkin(3) * t453 - t523;
t521 = t469 * t440;
t475 = t331 + t521;
t466 = t475 * t450;
t283 = -t464 - t466;
t546 = (-t392 + t428) * pkin(3);
t519 = t442 * t450;
t410 = pkin(4) * t519 - qJD(3) * qJ(5);
t543 = -pkin(4) * t517 - t410 * t519;
t542 = t418 * t450 + t499;
t488 = qJD(4) * qJD(3);
t541 = -qJ(4) * t412 - 0.2e1 * t488 - t549;
t491 = t446 + t447;
t397 = t491 * t441;
t406 = t491 * t440;
t337 = t397 * t451 + t406 * t454;
t340 = t397 * t454 - t406 * t451;
t284 = t337 * t455 + t340 * t452;
t536 = pkin(5) * t284;
t533 = pkin(6) * t337;
t527 = pkin(3) + qJ(5);
t526 = qJ(4) * t406;
t522 = qJ(4) * t453;
t391 = -t428 + t498;
t520 = t391 * qJ(5);
t330 = -pkin(2) * t441 - pkin(7) * t440 + t345;
t516 = t450 * t330;
t503 = t453 * t330;
t502 = t453 * t390;
t492 = g(3) * t450 - t331 * t453;
t489 = qJD(5) * t453;
t482 = pkin(2) * t392 - t560;
t480 = pkin(1) * t337 + pkin(2) * t406 + pkin(7) * t397;
t479 = qJD(4) * t519;
t314 = t450 * t331 + t529;
t269 = t314 * t450 - t453 * t492;
t371 = -t421 * t452 - t422 * t455;
t473 = t451 * t485;
t472 = t454 * t485;
t415 = qJDD(1) * t455 - t452 * t457;
t471 = -pkin(5) * t415 - g(3) * t452;
t468 = pkin(3) * t456 - qJDD(3) * qJ(4) - t453 * t521 + t492;
t467 = (pkin(3) * qJD(3) - 0.2e1 * qJD(4)) * t450;
t268 = t314 * t453 + t450 * t492;
t332 = t392 * t450 + t502;
t354 = t416 * t453 + t515;
t370 = t421 * t455 - t422 * t452;
t436 = 0.2e1 * t488;
t281 = t436 - t468;
t462 = t478 + t509;
t461 = -t391 * pkin(3) + t330 + (-t462 - t478) * qJ(4);
t460 = t461 + t543;
t459 = pkin(4) * t391 - qJ(5) * t517 + qJD(3) * t410 + qJDD(5) - t468;
t254 = t436 + t459;
t458 = -0.2e1 * qJD(5) * qJD(3) + t464 - qJ(5) * t411 + (t462 - t478) * pkin(4);
t427 = 0.2e1 * t479;
t255 = -pkin(3) * t428 + qJ(4) * t390 + t427 - t461;
t253 = t466 + t458;
t430 = pkin(3) * t509;
t426 = 0.2e1 * t442 * t489;
t414 = qJDD(1) * t452 + t455 * t457;
t387 = -pkin(5) * t414 + g(3) * t455;
t386 = -qJ(4) * t498 + t430;
t385 = t491 * t490;
t369 = t430 + (qJ(5) * t450 - t522) * t441;
t368 = qJDD(3) * t451 + t385 * t454;
t367 = -qJDD(3) * t454 + t385 * t451;
t366 = -t446 * t490 + t453 * t462;
t365 = -t391 * t450 - t447 * t490;
t347 = (t391 - t428) * t453;
t344 = -pkin(4) * t411 + qJ(4) * t392;
t335 = pkin(6) * t340;
t321 = t366 * t454 - t473;
t320 = t365 * t454 + t473;
t319 = t366 * t451 + t472;
t318 = t365 * t451 - t472;
t316 = pkin(4) * t412 + t390 * t527;
t300 = -t367 * t452 + t368 * t455;
t299 = t367 * t455 + t368 * t452;
t290 = pkin(1) * g(3) + pkin(6) * t474;
t289 = t503 + t531;
t288 = t516 - t532;
t287 = -t492 + t549;
t286 = t314 - t540;
t285 = -t337 * t452 + t340 * t455;
t282 = pkin(5) * t285;
t276 = -t319 * t452 + t321 * t455;
t275 = -t318 * t452 + t320 * t455;
t274 = t319 * t455 + t321 * t452;
t273 = t318 * t455 + t320 * t452;
t272 = t526 - t283;
t271 = pkin(3) * t406 + t281;
t270 = t442 * t467 + t461;
t256 = t461 - 0.2e1 * t479 + t546;
t250 = pkin(3) * t411 + t283 + t561;
t249 = pkin(3) * t417 + t468 + t541;
t246 = t268 * t454 - t533;
t245 = t268 * t451 + t335;
t244 = t526 + (pkin(4) * t441 + t475) * t450 + t458;
t243 = -t520 + (t467 - 0.2e1 * t489) * t442 + t460;
t242 = t269 * t454 + t330 * t451;
t241 = t269 * t451 - t330 * t454;
t240 = pkin(4) * t498 + t406 * t527 + t254;
t239 = t281 * t453 - t283 * t450;
t238 = t281 * t450 + t283 * t453;
t237 = -pkin(3) * t348 + t255 * t453 - t531;
t236 = -qJ(4) * t501 - t256 * t450 + t532;
t235 = pkin(4) * t417 + t255 + t426 + t520 - t543;
t234 = t417 * t527 - t459 + t541;
t233 = -t411 * t527 + t253 - t561;
t232 = pkin(4) * t419 + t426 + t427 + (t391 + t392) * qJ(5) - t546 - t460;
t231 = -t271 * t450 + t272 * t453;
t230 = -t287 * t451 + t289 * t454 + t534;
t229 = -t286 * t451 + t288 * t454 + t535;
t228 = t287 * t454 + t289 * t451 - t476;
t227 = t286 * t454 + t288 * t451 + t477;
t226 = t253 * t450 + t254 * t453;
t225 = -t253 * t453 + t254 * t450;
t224 = t231 * t454 - t386 * t451 - t533;
t223 = t231 * t451 + t386 * t454 + t335;
t222 = -t232 * t450 + t344 * t453 - t532;
t221 = t235 * t453 - t316 * t450 - t531;
t220 = t239 * t454 + t270 * t451;
t219 = t239 * t451 - t270 * t454;
t218 = pkin(4) * t253 - qJ(4) * t243;
t217 = -t241 * t452 + t242 * t455;
t216 = t241 * t455 + t242 * t452;
t215 = -t240 * t450 + t244 * t453;
t214 = -pkin(2) * t238 - pkin(3) * t283 - qJ(4) * t281;
t213 = t236 * t454 - t250 * t451 - t535;
t212 = t237 * t454 - t249 * t451 - t534;
t211 = -pkin(6) * t241 - (pkin(2) * t451 - pkin(7) * t454) * t268;
t210 = -pkin(7) * t238 + (pkin(3) * t450 - t522) * t270;
t209 = t215 * t454 - t369 * t451 - t533;
t208 = t215 * t451 + t369 * t454 + t335;
t207 = t236 * t451 + t250 * t454 - t477;
t206 = t237 * t451 + t249 * t454 + t476;
t205 = t226 * t454 + t243 * t451;
t204 = t226 * t451 - t243 * t454;
t203 = pkin(4) * t254 - t243 * t527;
t202 = pkin(6) * t242 - (-pkin(2) * t454 - pkin(7) * t451 - pkin(1)) * t268;
t201 = t221 * t454 - t234 * t451 - t534;
t200 = t222 * t454 - t233 * t451 + t535;
t199 = t221 * t451 + t234 * t454 + t476;
t198 = t222 * t451 + t233 * t454 + t477;
t197 = -t219 * t452 + t220 * t455;
t196 = t219 * t455 + t220 * t452;
t195 = -pkin(2) * t225 - qJ(4) * t254 + t253 * t527;
t194 = -t204 * t452 + t205 * t455;
t193 = t204 * t455 + t205 * t452;
t192 = -pkin(7) * t225 - t203 * t450 + t218 * t453;
t191 = -pkin(6) * t219 + t210 * t454 - t214 * t451;
t190 = -pkin(1) * t238 + pkin(6) * t220 + t210 * t451 + t214 * t454;
t189 = -pkin(6) * t204 + t192 * t454 - t195 * t451;
t188 = -pkin(1) * t225 + pkin(6) * t205 + t192 * t451 + t195 * t454;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t414, -t415, 0, t371, 0, 0, 0, 0, 0, 0, -t465, t341, 0, t248, 0, 0, 0, 0, 0, 0, t263, t266, t285, t217, 0, 0, 0, 0, 0, 0, t285, -t263, -t266, t197, 0, 0, 0, 0, 0, 0, t285, -t266, t263, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t415, -t414, 0, t370, 0, 0, 0, 0, 0, 0, -t341, -t465, 0, -t552, 0, 0, 0, 0, 0, 0, -t260, -t261, t284, t216, 0, 0, 0, 0, 0, 0, t284, t260, t261, t196, 0, 0, 0, 0, 0, 0, t284, t261, -t260, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t351, -t355, 0, -t268, 0, 0, 0, 0, 0, 0, 0, -t351, t355, t238, 0, 0, 0, 0, 0, 0, 0, t355, t351, t225; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t415, 0, -t414, 0, t471, -t387, -t370, -pkin(5) * t370, 0, 0, -t341, 0, -t465, 0, t557, t556, t552, pkin(5) * t552 + pkin(6) * t493 - t290 * t452, t276, t252, -t562, t275, t558, t300, -t227 * t452 + t229 * t455 + t538, -t228 * t452 + t230 * t455 + t537, -t245 * t452 + t246 * t455 - t536, -pkin(5) * t216 - t202 * t452 + t211 * t455, t300, t562, -t558, t276, t252, t275, -t223 * t452 + t224 * t455 - t536, -t207 * t452 + t213 * t455 - t538, -t206 * t452 + t212 * t455 - t537, -pkin(5) * t196 - t190 * t452 + t191 * t455, t300, -t558, -t562, t275, -t252, t276, -t208 * t452 + t209 * t455 - t536, -t199 * t452 + t201 * t455 - t537, -t198 * t452 + t200 * t455 + t538, -pkin(5) * t193 - t188 * t452 + t189 * t455; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t414, 0, t415, 0, t387, t471, t371, pkin(5) * t371, 0, 0, t465, 0, -t341, 0, -t556, t557, t248, pkin(5) * t248 + pkin(6) * t504 + t290 * t455, t274, -t559, t278, t273, -t277, t299, t227 * t455 + t229 * t452 + t257, t228 * t455 + t230 * t452 + t258, t245 * t455 + t246 * t452 + t282, pkin(5) * t217 + t202 * t455 + t211 * t452, t299, -t278, t277, t274, -t559, t273, t223 * t455 + t224 * t452 + t282, t207 * t455 + t213 * t452 - t257, t206 * t455 + t212 * t452 - t258, pkin(5) * t197 + t190 * t455 + t191 * t452, t299, t277, t278, t273, t559, t274, t208 * t455 + t209 * t452 + t282, t199 * t455 + t201 * t452 - t258, t198 * t455 + t200 * t452 + t257, pkin(5) * t194 + t188 * t455 + t189 * t452; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t421, t422, 0, 0, 0, 0, 0, 0, 0, t441, -pkin(1) * t405 - t345, -pkin(1) * t402 - t346, 0, -pkin(1) * t293, t348, t332, t354, t347, t542, 0, t482 - t503, -t481 + t516, t269 + t480, pkin(1) * t241 - pkin(2) * t330 + pkin(7) * t269, 0, -t354, -t542, t348, t332, t347, t271 * t453 + t272 * t450 + t480, t453 * t256 + (-pkin(2) - t523) * t392 + t560, pkin(3) * t502 + t255 * t450 + t481, pkin(1) * t219 + pkin(7) * t239 + (-pkin(2) + t469) * t270, 0, -t542, t354, t347, -t332, t348, t240 * t453 + t244 * t450 + t480, t235 * t450 + t316 * t453 + t481, t232 * t453 + t344 * t450 + t482, pkin(1) * t204 - pkin(2) * t243 + pkin(7) * t226 + t203 * t453 + t218 * t450;];
tauB_reg = t1;
