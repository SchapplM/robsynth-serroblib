% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:52:02
% EndTime: 2019-12-31 20:52:10
% DurationCPUTime: 5.20s
% Computational Cost: add. (10376->391), mult. (14817->474), div. (0->0), fcn. (7453->6), ass. (0->259)
t457 = qJD(3) ^ 2;
t443 = qJD(1) + qJD(2);
t441 = t443 ^ 2;
t451 = sin(qJ(3));
t446 = t451 ^ 2;
t519 = t446 * t441;
t420 = t457 + t519;
t454 = cos(qJ(3));
t488 = t454 * t441 * t451;
t416 = qJDD(3) - t488;
t502 = t454 * t416;
t368 = -t451 * t420 + t502;
t493 = qJD(3) * t454;
t481 = t443 * t493;
t442 = qJDD(1) + qJDD(2);
t510 = t451 * t442;
t395 = 0.2e1 * t481 + t510;
t452 = sin(qJ(2));
t455 = cos(qJ(2));
t315 = t452 * t368 + t455 * t395;
t318 = t455 * t368 - t452 * t395;
t453 = sin(qJ(1));
t456 = cos(qJ(1));
t268 = t456 * t315 + t453 * t318;
t536 = pkin(5) * t268;
t272 = t453 * t315 - t456 * t318;
t266 = pkin(5) * t272;
t533 = pkin(6) * t315;
t484 = pkin(1) * t315 + pkin(2) * t395 + pkin(7) * t368;
t514 = t451 * t416;
t362 = t454 * t420 + t514;
t478 = -pkin(1) * t362 + pkin(6) * t318;
t415 = qJDD(3) + t488;
t403 = t454 * t415;
t421 = -t457 + t519;
t367 = t451 * t421 + t403;
t498 = t455 * t442;
t330 = t452 * t367 - t451 * t498;
t507 = t452 * t442;
t334 = t455 * t367 + t451 * t507;
t284 = t456 * t330 + t453 * t334;
t559 = t453 * t330 - t456 * t334;
t447 = t454 ^ 2;
t495 = t446 + t447;
t402 = t495 * t442;
t410 = t495 * t441;
t343 = t452 * t402 + t455 * t410;
t347 = t455 * t402 - t452 * t410;
t290 = t456 * t343 + t453 * t347;
t535 = pkin(5) * t290;
t293 = t453 * t343 - t456 * t347;
t288 = pkin(5) * t293;
t356 = t451 * t395;
t494 = qJD(3) * t443;
t432 = t451 * t494;
t501 = t454 * t442;
t397 = -0.2e1 * t432 + t501;
t503 = t454 * t397;
t340 = -t503 + t356;
t411 = (t446 - t447) * t441;
t304 = t452 * t340 + t455 * t411;
t306 = t455 * t340 - t452 * t411;
t261 = t456 * t304 + t453 * t306;
t558 = t453 * t304 - t456 * t306;
t518 = t447 * t441;
t422 = -t457 + t518;
t366 = -t454 * t422 + t514;
t329 = t452 * t366 + t454 * t498;
t333 = t455 * t366 - t452 * t501;
t557 = t456 * t329 + t453 * t333;
t285 = t453 * t329 - t456 * t333;
t406 = t455 * t441 + t507;
t409 = t452 * t441 - t498;
t349 = t453 * t406 + t456 * t409;
t382 = pkin(6) * t406 - t455 * g(3);
t548 = pkin(6) * t409 - t452 * g(3);
t556 = pkin(5) * t349 + t453 * t382 + t456 * t548;
t468 = t456 * t406 - t453 * t409;
t555 = pkin(5) * t468 + t456 * t382 - t453 * t548;
t426 = t456 * g(1) + t453 * g(2);
t458 = qJD(1) ^ 2;
t413 = -t458 * pkin(1) - t426;
t425 = t453 * g(1) - t456 * g(2);
t466 = qJDD(1) * pkin(1) + t425;
t353 = t452 * t413 - t455 * t466;
t354 = t455 * t413 + t452 * t466;
t476 = t452 * t353 + t455 * t354;
t301 = t455 * t353 - t452 * t354;
t497 = t456 * t301;
t552 = -t453 * t476 + t497;
t506 = t453 * t301;
t256 = t456 * t476 + t506;
t549 = pkin(2) * t362;
t532 = pkin(6) * t343;
t341 = pkin(6) * t347;
t530 = pkin(7) * t362;
t483 = pkin(1) * t343 + pkin(2) * t410 + pkin(7) * t402;
t359 = t451 * t422 + t502;
t423 = -t457 - t518;
t360 = t451 * t423 + t403;
t539 = pkin(2) * t360;
t544 = -qJ(4) * t423 - t539;
t521 = t443 * t451;
t414 = -qJD(3) * pkin(4) - qJ(5) * t521;
t543 = t414 * t521 + qJDD(5);
t492 = qJD(4) * qJD(3);
t438 = -0.2e1 * t492;
t542 = -qJ(4) * t416 + t438 - t549;
t396 = -t432 + t501;
t541 = t396 * pkin(4) + t543;
t337 = -t441 * pkin(2) + t442 * pkin(7) + t354;
t472 = -pkin(3) * t454 - qJ(4) * t451;
t394 = t472 * t443;
t477 = t394 * t443 + t337;
t480 = qJDD(3) * pkin(3) + t457 * qJ(4) - qJDD(4);
t527 = t454 * g(3);
t460 = t451 * t477 - t480 + t527;
t540 = pkin(3) + pkin(4);
t515 = t451 * t415;
t365 = t454 * t423 - t515;
t314 = t452 * t365 + t455 * t397;
t317 = t455 * t365 - t452 * t397;
t267 = t456 * t314 + t453 * t317;
t537 = pkin(5) * t267;
t534 = pkin(6) * t314;
t531 = pkin(7) * t360;
t525 = qJ(4) * t410;
t522 = qJ(4) * t454;
t520 = t443 * t454;
t336 = -t442 * pkin(2) - t441 * pkin(7) + t353;
t517 = t451 * t336;
t516 = t451 * t397;
t505 = t454 * t336;
t504 = t454 * t395;
t496 = t451 * g(3) - t454 * t337;
t491 = 0.2e1 * t521;
t485 = pkin(1) * t314 + pkin(2) * t397 + pkin(7) * t365;
t482 = qJD(5) * t520;
t479 = -pkin(1) * t360 + pkin(6) * t317;
t320 = t451 * t337 + t527;
t275 = t451 * t320 - t454 * t496;
t376 = -t453 * t425 - t456 * t426;
t475 = t452 * t488;
t474 = t455 * t488;
t418 = t456 * qJDD(1) - t453 * t458;
t473 = -pkin(5) * t418 - t453 * g(3);
t471 = pkin(3) * t451 - t522;
t470 = t457 * pkin(3) - qJDD(3) * qJ(4) - t394 * t520 + t496;
t274 = t454 * t320 + t451 * t496;
t469 = t504 + t516;
t361 = -t454 * t421 + t515;
t375 = t456 * t425 - t453 * t426;
t437 = 0.2e1 * t492;
t287 = t437 - t470;
t464 = t481 + t510;
t465 = pkin(4) * t415 + t464 * qJ(5) + t480;
t463 = -t396 * pkin(3) + t336 + (-t464 - t481) * qJ(4);
t462 = pkin(4) * t518 - qJD(3) * t414 + t470;
t461 = qJD(4) * t491 - t463;
t459 = t396 * qJ(5) + t462;
t276 = (pkin(3) * qJD(3) - 0.2e1 * qJD(4)) * t521 + t463;
t264 = (t397 - t432) * pkin(3) + t461;
t263 = -pkin(3) * t432 + qJ(4) * t395 + t461;
t260 = (qJ(5) * t493 + (-0.2e1 * qJD(5) + t394) * t451) * t443 - t465 + t320;
t430 = 0.2e1 * t482;
t417 = t453 * qJDD(1) + t456 * t458;
t392 = -pkin(5) * t417 + t456 * g(3);
t391 = t471 * t442;
t390 = t495 * t494;
t377 = (-t540 * t451 + t522) * t442;
t374 = t452 * qJDD(3) + t455 * t390;
t373 = -t455 * qJDD(3) + t452 * t390;
t372 = -t446 * t494 + t454 * t464;
t371 = -t451 * t396 - t447 * t494;
t355 = (t396 - t432) * t454;
t352 = qJ(4) * t397 + qJ(5) * t415;
t327 = t455 * t372 - t475;
t326 = t455 * t371 + t475;
t325 = t452 * t372 + t474;
t324 = t452 * t371 - t474;
t322 = -qJ(5) * t416 + t540 * t395;
t308 = -t453 * t373 + t456 * t374;
t307 = t456 * t373 + t453 * t374;
t298 = pkin(1) * g(3) + pkin(6) * t476;
t297 = t505 + t530;
t296 = t517 - t531;
t295 = -t496 + t549;
t294 = t320 - t539;
t282 = -t453 * t325 + t456 * t327;
t281 = -t453 * t324 + t456 * t326;
t280 = t456 * t325 + t453 * t327;
t279 = t456 * t324 + t453 * t326;
t278 = t460 + t525;
t277 = pkin(3) * t410 + t287;
t270 = -t453 * t314 + t456 * t317;
t265 = pkin(5) * t270;
t259 = t437 - t459 - 0.2e1 * t482;
t258 = -pkin(3) * t415 + t460 + t544;
t257 = -pkin(3) * t420 + t470 + t542;
t254 = t455 * t274 - t532;
t253 = t452 * t274 + t341;
t252 = qJ(5) * t518 + t276 - t541;
t251 = -t525 + qJD(5) * t491 + (-qJ(5) * t494 - g(3)) * t454 + (qJ(5) * t442 - t477) * t451 + t465;
t250 = t455 * t275 + t452 * t336;
t249 = t452 * t275 - t455 * t336;
t248 = t430 + t438 - t540 * t410 + (t396 + t501) * qJ(5) + t462;
t247 = t454 * t287 + t451 * t460;
t246 = t451 * t287 - t454 * t460;
t245 = -pkin(3) * t356 + t454 * t263 - t530;
t244 = qJ(4) * t503 - t451 * t264 - t531;
t243 = t263 + (t420 - t518) * qJ(5) + t541;
t242 = (-t423 - t518) * qJ(5) + (t396 + t397) * pkin(4) + t264 + t543;
t241 = -t540 * t420 + t430 + t459 + t542;
t240 = -t540 * t415 + t260 + t544;
t239 = -t451 * t277 + t454 * t278;
t238 = -t452 * t295 + t455 * t297 + t533;
t237 = -t452 * t294 + t455 * t296 - t534;
t236 = t455 * t295 + t452 * t297 - t478;
t235 = t455 * t294 + t452 * t296 + t479;
t234 = -t451 * t242 + t454 * t352 - t531;
t233 = t454 * t243 - t451 * t322 - t530;
t232 = t455 * t239 - t452 * t391 - t532;
t231 = t452 * t239 + t455 * t391 + t341;
t230 = t454 * t259 + t451 * t260;
t229 = t451 * t259 - t454 * t260;
t228 = t455 * t247 + t452 * t276;
t227 = t452 * t247 - t455 * t276;
t226 = -qJ(4) * t252 - qJ(5) * t260;
t225 = -t453 * t249 + t456 * t250;
t224 = t456 * t249 + t453 * t250;
t223 = -t451 * t248 + t454 * t251;
t222 = -pkin(2) * t246 + pkin(3) * t460 - qJ(4) * t287;
t221 = t455 * t244 - t452 * t258 - t534;
t220 = t455 * t245 - t452 * t257 - t533;
t219 = -pkin(6) * t249 - (pkin(2) * t452 - pkin(7) * t455) * t274;
t218 = -pkin(7) * t246 + t276 * t471;
t217 = t455 * t223 - t452 * t377 + t532;
t216 = t452 * t223 + t455 * t377 - t341;
t215 = t452 * t244 + t455 * t258 + t479;
t214 = t452 * t245 + t455 * t257 + t478;
t213 = t455 * t230 + t452 * t252;
t212 = t452 * t230 - t455 * t252;
t211 = -qJ(5) * t259 - t540 * t252;
t210 = pkin(6) * t250 - (-pkin(2) * t455 - pkin(7) * t452 - pkin(1)) * t274;
t209 = t455 * t233 - t452 * t241 - t533;
t208 = t455 * t234 - t452 * t240 - t534;
t207 = t452 * t233 + t455 * t241 + t478;
t206 = t452 * t234 + t455 * t240 + t479;
t205 = -t453 * t227 + t456 * t228;
t204 = t456 * t227 + t453 * t228;
t203 = -pkin(2) * t229 - qJ(4) * t259 + t540 * t260;
t202 = -t453 * t212 + t456 * t213;
t201 = t456 * t212 + t453 * t213;
t200 = -pkin(7) * t229 - t451 * t211 + t454 * t226;
t199 = -pkin(6) * t227 + t455 * t218 - t452 * t222;
t198 = -pkin(1) * t246 + pkin(6) * t228 + t452 * t218 + t455 * t222;
t197 = -pkin(6) * t212 + t455 * t200 - t452 * t203;
t196 = -pkin(1) * t229 + pkin(6) * t213 + t452 * t200 + t455 * t203;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t417, -t418, 0, t376, 0, 0, 0, 0, 0, 0, -t468, t349, 0, t256, 0, 0, 0, 0, 0, 0, t270, t272, -t293, t225, 0, 0, 0, 0, 0, 0, t270, -t293, -t272, t205, 0, 0, 0, 0, 0, 0, t270, -t272, t293, t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t418, -t417, 0, t375, 0, 0, 0, 0, 0, 0, -t349, -t468, 0, -t552, 0, 0, 0, 0, 0, 0, t267, -t268, t290, t224, 0, 0, 0, 0, 0, 0, t267, t290, t268, t204, 0, 0, 0, 0, 0, 0, t267, t268, -t290, t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t360, -t362, 0, -t274, 0, 0, 0, 0, 0, 0, t360, 0, t362, t246, 0, 0, 0, 0, 0, 0, t360, t362, 0, t229; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t418, 0, -t417, 0, t473, -t392, -t375, -pkin(5) * t375, 0, 0, -t349, 0, -t468, 0, t556, t555, t552, pkin(5) * t552 + pkin(6) * t497 - t453 * t298, t282, t558, -t559, t281, t285, t308, -t453 * t235 + t456 * t237 - t537, -t453 * t236 + t456 * t238 + t536, -t453 * t253 + t456 * t254 - t535, -pkin(5) * t224 - t453 * t210 + t456 * t219, t282, -t559, -t558, t308, -t285, t281, -t453 * t215 + t456 * t221 - t537, -t453 * t231 + t456 * t232 - t535, -t453 * t214 + t456 * t220 - t536, -pkin(5) * t204 - t453 * t198 + t456 * t199, t282, -t558, t559, t281, t285, t308, -t453 * t206 + t456 * t208 - t537, -t453 * t207 + t456 * t209 - t536, -t453 * t216 + t456 * t217 + t535, -pkin(5) * t201 - t453 * t196 + t456 * t197; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t417, 0, t418, 0, t392, t473, t376, pkin(5) * t376, 0, 0, t468, 0, -t349, 0, -t555, t556, t256, pkin(5) * t256 + pkin(6) * t506 + t456 * t298, t280, -t261, t284, t279, -t557, t307, t456 * t235 + t453 * t237 + t265, t456 * t236 + t453 * t238 + t266, t456 * t253 + t453 * t254 - t288, pkin(5) * t225 + t456 * t210 + t453 * t219, t280, t284, t261, t307, t557, t279, t456 * t215 + t453 * t221 + t265, t456 * t231 + t453 * t232 - t288, t456 * t214 + t453 * t220 - t266, pkin(5) * t205 + t456 * t198 + t453 * t199, t280, t261, -t284, t279, -t557, t307, t456 * t206 + t453 * t208 + t265, t456 * t207 + t453 * t209 - t266, t456 * t216 + t453 * t217 + t288, pkin(5) * t202 + t456 * t196 + t453 * t197; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t425, t426, 0, 0, 0, 0, 0, 0, 0, t442, -pkin(1) * t409 - t353, -pkin(1) * t406 - t354, 0, -pkin(1) * t301, t356, t469, t361, t355, t359, 0, t485 - t505, -t484 + t517, t275 + t483, pkin(1) * t249 - pkin(2) * t336 + pkin(7) * t275, t356, t361, -t469, 0, -t359, t355, qJ(4) * t516 + t454 * t264 + t485, t454 * t277 + t451 * t278 + t483, pkin(3) * t504 + t451 * t263 + t484, pkin(1) * t227 + pkin(7) * t247 + (-pkin(2) + t472) * t276, t356, -t469, -t361, t355, t359, 0, t454 * t242 + t451 * t352 + t485, t451 * t243 + t454 * t322 + t484, t454 * t248 + t451 * t251 - t483, pkin(1) * t212 - pkin(2) * t252 + pkin(7) * t230 + t454 * t211 + t451 * t226;];
tauB_reg = t1;
