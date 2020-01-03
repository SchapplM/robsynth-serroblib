% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPPR6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:56
% EndTime: 2019-12-31 17:48:03
% DurationCPUTime: 6.37s
% Computational Cost: add. (13982->481), mult. (37326->710), div. (0->0), fcn. (22759->8), ass. (0->343)
t495 = sin(pkin(7));
t555 = qJD(1) * t495;
t536 = qJD(3) * t555;
t477 = -0.2e1 * t536;
t497 = cos(pkin(7));
t502 = qJD(1) ^ 2;
t499 = sin(qJ(1));
t501 = cos(qJ(1));
t471 = t499 * g(1) - t501 * g(2);
t519 = qJDD(2) - t471;
t584 = qJ(3) * t495;
t530 = -pkin(1) - t584;
t554 = t497 * qJD(1);
t490 = t495 ^ 2;
t493 = t497 ^ 2;
t557 = t490 + t493;
t588 = pkin(2) + qJ(4);
t605 = t477 + (-t557 * pkin(3) - qJ(2)) * t502 + (-t588 * t497 + t530) * qJDD(1) + t519 - 0.2e1 * qJD(4) * t554;
t472 = t501 * g(1) + t499 * g(2);
t444 = -t502 * pkin(1) + qJDD(1) * qJ(2) - t472;
t552 = qJD(1) * qJD(2);
t527 = t444 + 0.2e1 * t552;
t459 = t557 * t502;
t446 = t495 * t459;
t547 = t501 * qJDD(1);
t407 = -t499 * t446 + t495 * t547;
t602 = pkin(5) * t407;
t492 = t497 * t493;
t449 = (t490 * t497 + t492) * t502;
t531 = t497 * t547;
t409 = -t499 * t449 + t531;
t601 = pkin(5) * t409;
t548 = t499 * qJDD(1);
t410 = t501 * t446 + t495 * t548;
t600 = pkin(5) * t410;
t532 = t497 * t548;
t412 = t501 * t449 + t532;
t599 = pkin(5) * t412;
t498 = sin(qJ(5));
t494 = sin(pkin(8));
t500 = cos(qJ(5));
t433 = (t494 * t497 * t498 + t495 * t500) * qJD(1);
t434 = t500 * t494 * t554 - t498 * t555;
t388 = t433 * t434;
t496 = cos(pkin(8));
t550 = qJDD(1) * t497;
t533 = t496 * t550;
t513 = qJDD(5) + t533;
t595 = -t388 + t513;
t598 = t498 * t595;
t597 = t500 * t595;
t589 = t497 * g(3);
t593 = pkin(2) * t497;
t517 = -t584 - t593;
t452 = t517 * qJD(1);
t594 = t452 * t555 + qJDD(3);
t518 = t589 + t594;
t568 = t497 * t502;
t587 = pkin(3) * qJDD(1);
t353 = (-qJ(4) * t568 + t527 + t587) * t495 + t518;
t297 = -t496 * t353 + t605 * t494;
t577 = t493 * t502;
t298 = t494 * t353 + t605 * t496;
t398 = -t495 * g(3) + t527 * t497;
t468 = t496 * t554 + qJD(5);
t400 = t468 * t433;
t534 = t494 * t550;
t549 = t495 * qJDD(1);
t540 = -t433 * qJD(5) - t498 * t549 + t500 * t534;
t508 = t540 - t400;
t559 = t498 * t534 + t500 * t549;
t337 = (qJD(5) - t468) * t434 + t559;
t430 = t433 ^ 2;
t431 = t434 ^ 2;
t463 = t468 ^ 2;
t592 = pkin(4) * t494;
t591 = pkin(4) * t496;
t484 = t490 * qJDD(1);
t485 = t493 * qJDD(1);
t455 = t485 + t484;
t401 = t499 * t455 + t501 * t459;
t590 = pkin(5) * t401;
t586 = qJ(2) * t446;
t585 = qJ(2) * t449;
t583 = qJDD(1) * pkin(1);
t582 = t468 * t498;
t581 = t468 * t500;
t489 = t494 ^ 2;
t580 = t489 * t493;
t579 = t490 * t502;
t491 = t496 ^ 2;
t578 = t491 * t493;
t512 = -qJ(4) * t577 + qJDD(4) + t398;
t556 = qJD(1) * t452;
t359 = (t556 + t587) * t497 + t512;
t576 = t494 * t359;
t569 = t496 * t502;
t544 = t494 * t569;
t458 = t493 * t544;
t440 = -t458 + t549;
t575 = t494 * t440;
t574 = t495 * t497;
t573 = t495 * t502;
t572 = t496 * t359;
t439 = t458 + t549;
t571 = t496 * t439;
t570 = t496 * t440;
t522 = pkin(6) * t494 + t591;
t506 = t522 * t577;
t279 = -pkin(4) * t549 - pkin(6) * t579 - t494 * t506 + t297;
t567 = t498 * t279;
t367 = t388 + t513;
t566 = t498 * t367;
t507 = t502 * qJ(2) - t519;
t438 = t507 + t583;
t565 = t499 * t438;
t564 = t500 * t279;
t563 = t500 * t367;
t562 = t501 * t438;
t280 = -pkin(4) * t579 + pkin(6) * t549 - t496 * t506 + t298;
t511 = pkin(3) + t522;
t521 = pkin(6) * t496 - t592;
t316 = (qJDD(1) * t511 + t521 * t573 + t556) * t497 + t512;
t255 = t500 * t280 + t498 * t316;
t560 = pkin(1) * t459 + qJ(2) * t455;
t551 = qJDD(1) * t494;
t546 = t430 + t431;
t545 = t494 * t388;
t543 = t495 * t569;
t542 = t495 * t568;
t541 = t496 * t388;
t539 = pkin(1) + t593;
t535 = t496 * t551;
t529 = (t489 + t491) * t502;
t528 = t438 + t583;
t254 = t498 * t280 - t500 * t316;
t232 = t498 * t254 + t500 * t255;
t397 = t495 * t527 + t589;
t344 = t495 * t397 + t497 * t398;
t415 = -t499 * t471 - t501 * t472;
t526 = t492 * t544;
t525 = t494 * t543;
t523 = t495 * t531;
t467 = -t499 * t502 + t547;
t520 = -pkin(5) * t467 - t499 * g(3);
t516 = pkin(2) * t495 - qJ(3) * t497;
t515 = t495 * t458;
t514 = (-t490 - t580) * t502;
t231 = -t500 * t254 + t498 * t255;
t256 = -t496 * t297 + t494 * t298;
t257 = t494 * t297 + t496 * t298;
t343 = t497 * t397 - t495 * t398;
t414 = t501 * t471 - t499 * t472;
t466 = t501 * t502 + t548;
t510 = pkin(1) - t517;
t509 = -0.2e1 * t495 * t552 - t518;
t369 = t452 * t554 + t398;
t505 = t507 + 0.2e1 * t536;
t473 = 0.2e1 * t497 * t549;
t460 = (-t490 + t493) * t502;
t457 = t494 * t542;
t456 = t485 - t484;
t450 = t516 * qJDD(1);
t448 = (-t490 - t578) * t502;
t447 = (-t490 + t578) * t502;
t445 = (t490 - t580) * t502;
t443 = -pkin(5) * t466 + t501 * g(3);
t442 = (t489 - t491) * t577;
t441 = t493 * t529;
t429 = -t457 - t533;
t428 = -t457 + t533;
t427 = (-t543 + t551) * t497;
t426 = (t543 + t551) * t497;
t425 = t494 * t439;
t423 = -t499 * t542 + t523;
t422 = t466 * t574;
t421 = t529 * t574;
t420 = (qJDD(1) * t489 + t525) * t497;
t419 = (-t491 * t573 - t535) * t497;
t418 = (-t489 * t573 + t535) * t497;
t417 = (qJDD(1) * t491 - t525) * t497;
t404 = t501 * t456 - t499 * t460;
t403 = t501 * t455 - t499 * t459;
t402 = t499 * t456 + t501 * t460;
t399 = pkin(5) * t403;
t396 = -t431 + t463;
t395 = t430 - t463;
t393 = t477 + (t530 - 0.2e1 * t593) * qJDD(1) - t507;
t392 = qJDD(1) * t510 + t505;
t391 = (t539 + 0.2e1 * t584) * qJDD(1) + t505;
t390 = -t495 * t420 - t526;
t389 = -t495 * t417 + t526;
t387 = -t494 * t514 - t570;
t386 = t496 * t448 - t425;
t385 = t496 * t514 - t575;
t384 = -t496 * t447 + t575;
t383 = -t496 * t445 - t425;
t382 = t494 * t448 + t571;
t381 = -t494 * t447 - t570;
t380 = t494 * t445 - t571;
t379 = -t431 + t430;
t377 = t434 * qJD(5) + t559;
t376 = -t431 - t463;
t375 = -pkin(2) * t484 + t497 * t391;
t374 = -qJ(3) * t485 - t495 * t393;
t373 = -t494 * t427 + t496 * t429;
t372 = t496 * t427 + t494 * t429;
t371 = t496 * t426 + t494 * t428;
t370 = -t494 * t426 + t496 * t428;
t365 = -t463 - t430;
t364 = -t495 * t444 + t509;
t361 = pkin(2) * t459 + t369;
t358 = qJ(3) * t459 + t397 + t594;
t357 = -pkin(2) * t579 + (qJ(3) * t573 - t556) * t497 - t398;
t356 = (t433 * t500 - t434 * t498) * t468;
t355 = (t433 * t498 + t434 * t500) * t468;
t354 = -qJ(3) * t577 + (pkin(2) * t568 - t444) * t495 + t509;
t352 = -t495 * t381 + t497 * t429;
t351 = t495 * t382 + t497 * t428;
t350 = t495 * t385 - t497 * t426;
t349 = -t495 * t383 - t497 * t427;
t348 = -t497 * t382 + t495 * t428;
t347 = -t497 * t385 - t495 * t426;
t340 = t400 + t540;
t338 = (qJD(5) + t468) * t434 + t559;
t336 = t495 * t372 - t497 * t441;
t335 = -t495 * t371 + t497 * t442;
t334 = -t497 * t372 - t495 * t441;
t331 = t434 * t582 - t500 * t540;
t330 = -t434 * t581 - t498 * t540;
t329 = -t498 * t377 - t433 * t581;
t328 = t500 * t377 - t433 * t582;
t327 = -t494 * t356 + t496 * t513;
t326 = -t496 * t356 - t494 * t513;
t325 = t500 * t395 - t566;
t324 = -t498 * t396 + t597;
t323 = t498 * t395 + t563;
t322 = t500 * t396 + t598;
t321 = t501 * t344 - t565;
t320 = t499 * t344 + t562;
t319 = -t498 * t376 - t563;
t318 = t500 * t376 - t566;
t317 = pkin(3) * t372 - qJ(3) * t373;
t315 = t500 * t365 - t598;
t314 = t498 * t365 + t597;
t313 = -t495 * t364 + t497 * t369;
t312 = t497 * t364 + t495 * t369;
t309 = t501 * t351 + t499 * t386;
t308 = t501 * t350 + t499 * t387;
t307 = t499 * t351 - t501 * t386;
t306 = t499 * t350 - t501 * t387;
t305 = -t494 * t331 + t541;
t304 = -t494 * t329 - t541;
t303 = -t496 * t329 + t545;
t302 = -t496 * t331 - t545;
t301 = t501 * t336 + t499 * t373;
t300 = t499 * t336 - t501 * t373;
t299 = t497 * t358 - t495 * t361;
t295 = t500 * t338 + t498 * t508;
t294 = t500 * t337 - t498 * t340;
t293 = t498 * t338 - t500 * t508;
t292 = t498 * t337 + t500 * t340;
t291 = t501 * t313 - t499 * t392;
t290 = t499 * t313 + t501 * t392;
t289 = -t495 * t327 + t497 * t355;
t288 = -t494 * t325 + t337 * t496;
t287 = -t494 * t324 - t496 * t340;
t286 = -t496 * t324 + t494 * t340;
t285 = -t496 * t325 - t337 * t494;
t284 = pkin(3) * t428 - t588 * t386 + t572;
t283 = -pkin(3) * t426 - t588 * t387 - t576;
t282 = t496 * t319 - t494 * t508;
t281 = t494 * t319 + t496 * t508;
t276 = t496 * t315 - t494 * t338;
t275 = t494 * t315 + t496 * t338;
t274 = -qJ(2) * t312 - t392 * t516;
t273 = -t496 * t295 + t494 * t379;
t272 = -t494 * t295 - t496 * t379;
t271 = -t495 * t305 + t497 * t330;
t270 = -t495 * t304 + t497 * t328;
t269 = -pkin(1) * t312 - pkin(2) * t364 - qJ(3) * t369;
t268 = t496 * t294 - t494 * t546;
t267 = t494 * t294 + t496 * t546;
t266 = pkin(3) * t385 - qJ(3) * t387 - t298;
t265 = pkin(3) * t382 - qJ(3) * t386 - t297;
t264 = -pkin(1) * t348 - qJ(3) * t428 + t588 * t382 - t576;
t263 = -pkin(1) * t347 + qJ(3) * t426 + t588 * t385 - t572;
t262 = -t495 * t288 + t497 * t323;
t261 = -t495 * t287 + t497 * t322;
t260 = -pkin(6) * t318 + t564;
t259 = t495 * t281 + t497 * t318;
t258 = -t497 * t281 + t495 * t318;
t253 = -pkin(6) * t314 + t567;
t251 = t495 * t275 + t497 * t314;
t250 = -t497 * t275 + t495 * t314;
t249 = -t495 * t272 + t497 * t293;
t248 = t495 * t256 + t497 * t359;
t247 = -t497 * t256 + t495 * t359;
t246 = t495 * t267 + t497 * t292;
t245 = -t497 * t267 + t495 * t292;
t244 = -pkin(4) * t318 + t255;
t243 = -pkin(4) * t314 + t254;
t242 = -pkin(3) * t441 - t588 * t373 - t257;
t241 = -qJ(2) * t347 + t497 * t266 - t495 * t283;
t240 = -qJ(2) * t348 + t497 * t265 - t495 * t284;
t239 = t501 * t259 + t499 * t282;
t238 = t499 * t259 - t501 * t282;
t237 = -pkin(1) * t334 + qJ(3) * t441 + t588 * t372 + t256;
t236 = t501 * t251 + t499 * t276;
t235 = t499 * t251 - t501 * t276;
t234 = pkin(3) * t256 - qJ(3) * t257;
t233 = -qJ(2) * t334 - t495 * t242 + t497 * t317;
t230 = t501 * t246 + t499 * t268;
t229 = t499 * t246 - t501 * t268;
t228 = pkin(3) * t359 - t588 * t257;
t227 = t501 * t248 + t499 * t257;
t226 = t499 * t248 - t501 * t257;
t225 = pkin(3) * t281 + pkin(4) * t508 + pkin(6) * t319 - qJ(3) * t282 + t567;
t224 = -pkin(6) * t292 - t231;
t223 = pkin(3) * t275 + pkin(4) * t338 + pkin(6) * t315 - qJ(3) * t276 - t564;
t222 = t496 * t232 + t494 * t279;
t221 = t494 * t232 - t496 * t279;
t220 = -pkin(1) * t247 - qJ(3) * t359 + t588 * t256;
t219 = pkin(3) * t318 - t496 * t244 - t494 * t260 - t588 * t282;
t218 = pkin(3) * t314 - t496 * t243 - t494 * t253 - t588 * t276;
t217 = pkin(3) * t267 + pkin(4) * t546 + pkin(6) * t294 - qJ(3) * t268 + t232;
t216 = t495 * t221 + t497 * t231;
t215 = -t497 * t221 + t495 * t231;
t214 = -pkin(1) * t258 - qJ(3) * t318 + t494 * t244 - t496 * t260 + t588 * t281;
t213 = -t494 * t224 + (pkin(3) + t591) * t292 - t588 * t268;
t212 = -pkin(1) * t250 - qJ(3) * t314 + t494 * t243 - t496 * t253 + t588 * t275;
t211 = -qJ(2) * t247 - t495 * t228 + t497 * t234;
t210 = -pkin(1) * t245 - t496 * t224 + (-qJ(3) - t592) * t292 + t588 * t267;
t209 = -qJ(2) * t258 - t495 * t219 + t497 * t225;
t208 = t501 * t216 + t499 * t222;
t207 = t499 * t216 - t501 * t222;
t206 = -qJ(2) * t250 - t495 * t218 + t497 * t223;
t205 = pkin(3) * t221 - pkin(4) * t279 + pkin(6) * t232 - qJ(3) * t222;
t204 = -qJ(2) * t245 - t495 * t213 + t497 * t217;
t203 = -t588 * t222 + t511 * t231;
t202 = -pkin(1) * t215 + t588 * t221 + (-qJ(3) + t521) * t231;
t201 = -qJ(2) * t215 - t495 * t203 + t497 * t205;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t466, -t467, 0, t415, 0, 0, 0, 0, 0, 0, -t412, t410, t403, t321, 0, 0, 0, 0, 0, 0, t403, t412, -t410, t291, 0, 0, 0, 0, 0, 0, t309, t308, t301, t227, 0, 0, 0, 0, 0, 0, t236, t239, t230, t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t467, -t466, 0, t414, 0, 0, 0, 0, 0, 0, t409, -t407, t401, t320, 0, 0, 0, 0, 0, 0, t401, -t409, t407, t290, 0, 0, 0, 0, 0, 0, t307, t306, t300, t226, 0, 0, 0, 0, 0, 0, t235, t238, t229, t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t343, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, 0, 0, 0, 0, 0, 0, t348, t347, t334, t247, 0, 0, 0, 0, 0, 0, t250, t258, t245, t215; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t467, 0, -t466, 0, t520, -t443, -t414, -pkin(5) * t414, t423, t404, t410, -t423, t412, 0, -t499 * t397 - t495 * t562 - t601, -t499 * t398 - t497 * t562 + t602, t501 * t343 - t590, -pkin(5) * t320 - (pkin(1) * t499 - qJ(2) * t501) * t343, 0, -t410, -t412, t423, t404, -t423, t501 * t299 - t499 * t450 - t590, -t499 * t354 + t501 * t374 + t601, -t499 * t357 + t501 * t375 - t602, -pkin(5) * t290 - t499 * t269 + t501 * t274, t501 * t390 - t499 * t418, t501 * t335 - t499 * t370, t501 * t349 - t499 * t380, t501 * t389 - t499 * t419, t501 * t352 - t499 * t384, -t499 * t421 + t523, -pkin(5) * t307 + t501 * t240 - t499 * t264, -pkin(5) * t306 + t501 * t241 - t499 * t263, -pkin(5) * t300 + t501 * t233 - t499 * t237, -pkin(5) * t226 + t501 * t211 - t499 * t220, t501 * t271 - t499 * t302, t501 * t249 - t499 * t273, t501 * t261 - t499 * t286, t501 * t270 - t499 * t303, t501 * t262 - t499 * t285, t501 * t289 - t499 * t326, -pkin(5) * t235 + t501 * t206 - t499 * t212, -pkin(5) * t238 + t501 * t209 - t499 * t214, -pkin(5) * t229 + t501 * t204 - t499 * t210, -pkin(5) * t207 + t501 * t201 - t499 * t202; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t466, 0, t467, 0, t443, t520, t415, pkin(5) * t415, t422, t402, -t407, -t422, -t409, 0, t501 * t397 - t495 * t565 - t599, t501 * t398 - t497 * t565 + t600, t499 * t343 + t399, pkin(5) * t321 - (-pkin(1) * t501 - qJ(2) * t499) * t343, 0, t407, t409, t422, t402, -t422, t499 * t299 + t501 * t450 + t399, t501 * t354 + t499 * t374 + t599, t501 * t357 + t499 * t375 - t600, pkin(5) * t291 + t501 * t269 + t499 * t274, t499 * t390 + t501 * t418, t499 * t335 + t501 * t370, t499 * t349 + t501 * t380, t499 * t389 + t501 * t419, t499 * t352 + t501 * t384, t501 * t421 + t495 * t532, pkin(5) * t309 + t499 * t240 + t501 * t264, pkin(5) * t308 + t499 * t241 + t501 * t263, pkin(5) * t301 + t499 * t233 + t501 * t237, pkin(5) * t227 + t499 * t211 + t501 * t220, t499 * t271 + t501 * t302, t499 * t249 + t501 * t273, t499 * t261 + t501 * t286, t499 * t270 + t501 * t303, t499 * t262 + t501 * t285, t499 * t289 + t501 * t326, pkin(5) * t236 + t499 * t206 + t501 * t212, pkin(5) * t239 + t499 * t209 + t501 * t214, pkin(5) * t230 + t499 * t204 + t501 * t210, pkin(5) * t208 + t499 * t201 + t501 * t202; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t471, t472, 0, 0, t484, t473, 0, t485, 0, 0, t497 * t528 - t585, -t495 * t528 + t586, t344 + t560, pkin(1) * t438 + qJ(2) * t344, 0, 0, 0, t484, t473, t485, t495 * t358 + t497 * t361 + t560, t585 + (qJDD(1) * t530 + t393) * t497, -t586 + (qJDD(1) * t539 + t391) * t495, qJ(2) * t313 + t392 * t510, t497 * t420 - t515, t497 * t371 + t495 * t442, t497 * t383 - t495 * t427, t497 * t417 + t515, t497 * t381 + t495 * t429, t484, -pkin(1) * t386 + qJ(2) * t351 + t495 * t265 + t497 * t284, -pkin(1) * t387 + qJ(2) * t350 + t495 * t266 + t497 * t283, -pkin(1) * t373 + qJ(2) * t336 + t497 * t242 + t495 * t317, -pkin(1) * t257 + qJ(2) * t248 + t497 * t228 + t495 * t234, t497 * t305 + t495 * t330, t497 * t272 + t495 * t293, t497 * t287 + t495 * t322, t497 * t304 + t495 * t328, t497 * t288 + t495 * t323, t497 * t327 + t495 * t355, -pkin(1) * t276 + qJ(2) * t251 + t497 * t218 + t495 * t223, -pkin(1) * t282 + qJ(2) * t259 + t497 * t219 + t495 * t225, -pkin(1) * t268 + qJ(2) * t246 + t497 * t213 + t495 * t217, -pkin(1) * t222 + qJ(2) * t216 + t497 * t203 + t495 * t205;];
tauB_reg = t1;
