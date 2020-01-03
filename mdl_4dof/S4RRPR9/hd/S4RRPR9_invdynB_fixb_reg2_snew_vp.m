% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPR9_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:09
% EndTime: 2019-12-31 17:10:15
% DurationCPUTime: 5.21s
% Computational Cost: add. (18014->435), mult. (39630->657), div. (0->0), fcn. (26459->8), ass. (0->307)
t476 = sin(pkin(7));
t479 = sin(qJ(2));
t504 = qJD(1) * qJD(2);
t466 = t479 * t504;
t482 = cos(qJ(2));
t501 = t482 * qJDD(1);
t446 = -t466 + t501;
t477 = cos(pkin(7));
t507 = qJD(1) * t479;
t436 = -t477 * qJD(2) + t476 * t507;
t438 = t476 * qJD(2) + t477 * t507;
t528 = t438 * t436;
t487 = -t446 - t528;
t539 = t476 * t487;
t538 = t477 * t487;
t478 = sin(qJ(4));
t439 = -qJDD(4) + t446;
t481 = cos(qJ(4));
t395 = t481 * t436 + t478 * t438;
t397 = -t478 * t436 + t481 * t438;
t529 = t397 * t395;
t486 = -t439 - t529;
t537 = t478 * t486;
t536 = t481 * t486;
t506 = t482 * qJD(1);
t463 = -qJD(4) + t506;
t384 = t395 * t463;
t495 = t482 * t504;
t503 = t479 * qJDD(1);
t445 = t495 + t503;
t419 = t476 * qJDD(2) + t477 * t445;
t493 = -t477 * qJDD(2) + t476 * t445;
t485 = t395 * qJD(4) - t481 * t419 + t478 * t493;
t535 = t384 - t485;
t425 = t436 * t506;
t389 = -t419 + t425;
t426 = t438 * t506;
t387 = -t493 - t426;
t494 = t478 * t419 + t481 * t493;
t313 = (qJD(4) + t463) * t397 + t494;
t393 = t395 ^ 2;
t394 = t397 ^ 2;
t534 = t436 ^ 2;
t435 = t438 ^ 2;
t461 = t463 ^ 2;
t533 = qJD(2) ^ 2;
t532 = pkin(2) * t479;
t531 = pkin(2) * t482;
t530 = t482 * g(3);
t527 = t463 * t478;
t526 = t463 * t481;
t473 = t479 ^ 2;
t484 = qJD(1) ^ 2;
t525 = t473 * t484;
t480 = sin(qJ(1));
t483 = cos(qJ(1));
t455 = t480 * g(1) - t483 * g(2);
t431 = qJDD(1) * pkin(1) + t484 * pkin(5) + t455;
t488 = t445 + t495;
t358 = -t488 * qJ(3) + (-t446 + t466) * pkin(2) - t431;
t456 = t483 * g(1) + t480 * g(2);
t432 = -t484 * pkin(1) + qJDD(1) * pkin(5) - t456;
t414 = -t479 * g(3) + t482 * t432;
t489 = -qJ(3) * t479 - t531;
t443 = t489 * qJD(1);
t367 = -t533 * pkin(2) + qJDD(2) * qJ(3) + t443 * t506 + t414;
t320 = 0.2e1 * qJD(3) * t438 - t477 * t358 + t476 * t367;
t287 = t487 * pkin(3) + t389 * pkin(6) - t320;
t321 = -0.2e1 * qJD(3) * t436 + t476 * t358 + t477 * t367;
t420 = -pkin(3) * t506 - t438 * pkin(6);
t292 = -t534 * pkin(3) - pkin(6) * t493 + t420 * t506 + t321;
t251 = -t481 * t287 + t478 * t292;
t252 = t478 * t287 + t481 * t292;
t224 = -t481 * t251 + t478 * t252;
t524 = t476 * t224;
t366 = t530 + qJDD(3) - t533 * qJ(3) - qJDD(2) * pkin(2) + (qJD(1) * t443 + t432) * t479;
t523 = t476 * t366;
t390 = t446 - t528;
t522 = t476 * t390;
t521 = t477 * t224;
t520 = t477 * t366;
t519 = t477 * t390;
t323 = pkin(3) * t493 - t534 * pkin(6) + t438 * t420 + t366;
t518 = t478 * t323;
t338 = t439 - t529;
t517 = t478 * t338;
t516 = t479 * t431;
t515 = t479 * t446;
t462 = t482 * t484 * t479;
t453 = qJDD(2) + t462;
t514 = t479 * t453;
t454 = qJDD(2) - t462;
t513 = t479 * t454;
t512 = t481 * t323;
t511 = t481 * t338;
t510 = t482 * t431;
t509 = t482 * t454;
t474 = t482 ^ 2;
t508 = t473 + t474;
t502 = t480 * qJDD(1);
t500 = t483 * qJDD(1);
t499 = t479 * t529;
t498 = t479 * t528;
t497 = t482 * t529;
t496 = t482 * t528;
t225 = t478 * t251 + t481 * t252;
t413 = t479 * t432 + t530;
t361 = t479 * t413 + t482 * t414;
t405 = -t480 * t455 - t483 * t456;
t492 = t480 * t462;
t491 = t483 * t462;
t450 = -t480 * t484 + t500;
t490 = -pkin(4) * t450 - t480 * g(3);
t279 = -t477 * t320 + t476 * t321;
t280 = t476 * t320 + t477 * t321;
t360 = t482 * t413 - t479 * t414;
t404 = t483 * t455 - t480 * t456;
t470 = t474 * t484;
t460 = -t470 - t533;
t459 = t470 - t533;
t458 = -t525 - t533;
t457 = -t525 + t533;
t452 = t470 - t525;
t451 = t470 + t525;
t449 = t483 * t484 + t502;
t448 = t508 * qJDD(1);
t447 = -0.2e1 * t466 + t501;
t444 = 0.2e1 * t495 + t503;
t441 = t482 * t453;
t440 = t508 * t504;
t433 = t482 * t446;
t427 = -pkin(4) * t449 + t483 * g(3);
t423 = -t435 - t470;
t422 = -t435 + t470;
t421 = -t470 + t534;
t417 = t482 * t445 - t473 * t504;
t416 = -t474 * t504 - t515;
t412 = -t479 * t458 - t509;
t411 = -t479 * t457 + t441;
t410 = t482 * t460 - t514;
t409 = t482 * t459 - t513;
t408 = t482 * t458 - t513;
t407 = t479 * t460 + t441;
t402 = t483 * t448 - t480 * t451;
t401 = t480 * t448 + t483 * t451;
t400 = -t435 + t534;
t399 = -t479 * t444 + t482 * t447;
t398 = -t470 - t534;
t388 = t419 + t425;
t386 = t426 - t493;
t382 = t435 + t534;
t381 = (t436 * t477 - t438 * t476) * t506;
t380 = (-t436 * t476 - t438 * t477) * t506;
t379 = -t394 + t461;
t378 = t393 - t461;
t377 = t483 * t412 + t480 * t444;
t376 = t483 * t410 - t480 * t447;
t375 = t480 * t412 - t483 * t444;
t374 = t480 * t410 + t483 * t447;
t373 = -pkin(5) * t408 - t510;
t372 = -pkin(5) * t407 - t516;
t371 = t477 * t419 + t476 * t426;
t370 = -t476 * t419 + t477 * t426;
t369 = -t477 * t425 + t476 * t493;
t368 = t476 * t425 + t477 * t493;
t365 = -pkin(1) * t408 + t414;
t364 = -pkin(1) * t407 + t413;
t363 = -t394 - t461;
t355 = t482 * t381 - t515;
t354 = t477 * t421 + t522;
t353 = -t476 * t422 + t538;
t352 = -t476 * t423 + t519;
t351 = -t476 * t421 + t519;
t350 = -t477 * t422 - t539;
t349 = t477 * t423 + t522;
t348 = -t394 + t393;
t347 = t483 * t361 - t480 * t431;
t346 = t480 * t361 + t483 * t431;
t345 = -t461 - t393;
t344 = t477 * t398 - t539;
t343 = t476 * t398 + t538;
t342 = t482 * t371 + t498;
t341 = t482 * t369 - t498;
t336 = -t397 * qJD(4) - t494;
t335 = t477 * t387 - t476 * t389;
t334 = t477 * t386 - t476 * t388;
t333 = t476 * t387 + t477 * t389;
t332 = -t476 * t386 - t477 * t388;
t331 = (t395 * t481 - t397 * t478) * t463;
t330 = (t395 * t478 + t397 * t481) * t463;
t329 = -t393 - t394;
t328 = t482 * t354 + t387 * t479;
t327 = t482 * t353 - t479 * t389;
t326 = t482 * t352 + t479 * t388;
t325 = t479 * t352 - t482 * t388;
t324 = -qJ(3) * t349 + t520;
t322 = t482 * t334 - t479 * t400;
t319 = t482 * t344 - t479 * t386;
t318 = t479 * t344 + t482 * t386;
t317 = t384 + t485;
t312 = (qJD(4) - t463) * t397 + t494;
t311 = t481 * t378 + t517;
t310 = -t478 * t379 + t536;
t309 = t478 * t378 - t511;
t308 = t481 * t379 + t537;
t307 = t397 * t527 - t481 * t485;
t306 = -t397 * t526 - t478 * t485;
t305 = -t478 * t336 - t395 * t526;
t304 = t481 * t336 - t395 * t527;
t303 = -qJ(3) * t343 + t523;
t302 = t482 * t335 - t479 * t382;
t301 = t479 * t335 + t482 * t382;
t300 = -t478 * t363 + t511;
t299 = t481 * t363 + t517;
t298 = t481 * t345 - t537;
t297 = t478 * t345 + t536;
t296 = -t476 * t330 + t477 * t331;
t295 = -t477 * t330 - t476 * t331;
t294 = t483 * t326 + t480 * t349;
t293 = t480 * t326 - t483 * t349;
t290 = -pkin(2) * t349 + t321;
t289 = t482 * t296 - t479 * t439;
t288 = -pkin(2) * t343 + t320;
t284 = t483 * t319 + t480 * t343;
t283 = t480 * t319 - t483 * t343;
t282 = t483 * t302 + t480 * t333;
t281 = t480 * t302 - t483 * t333;
t278 = -t313 * t481 - t478 * t317;
t277 = -t481 * t312 - t478 * t535;
t276 = -t313 * t478 + t481 * t317;
t275 = -t478 * t312 + t481 * t535;
t274 = -pkin(6) * t299 + t512;
t273 = -t476 * t309 + t477 * t311;
t272 = -t476 * t308 + t477 * t310;
t271 = -t477 * t309 - t476 * t311;
t270 = -t477 * t308 - t476 * t310;
t269 = -t476 * t306 + t477 * t307;
t268 = -t476 * t304 + t477 * t305;
t267 = -t477 * t306 - t476 * t307;
t266 = -t477 * t304 - t476 * t305;
t265 = -pkin(1) * t325 + pkin(2) * t388 - qJ(3) * t352 - t523;
t264 = -t476 * t299 + t477 * t300;
t263 = t477 * t299 + t476 * t300;
t262 = -pkin(6) * t297 + t518;
t261 = -pkin(1) * t318 - pkin(2) * t386 - qJ(3) * t344 + t520;
t260 = t482 * t280 + t479 * t366;
t259 = t479 * t280 - t482 * t366;
t258 = -t476 * t297 + t477 * t298;
t257 = t477 * t297 + t476 * t298;
t256 = -qJ(3) * t333 - t279;
t255 = t482 * t269 + t499;
t254 = t482 * t268 - t499;
t253 = -pkin(3) * t535 + pkin(6) * t300 + t518;
t250 = t482 * t273 - t479 * t313;
t249 = t482 * t272 - t479 * t317;
t247 = t482 * t264 + t479 * t535;
t246 = t479 * t264 - t482 * t535;
t245 = -pkin(3) * t312 + pkin(6) * t298 - t512;
t244 = -pkin(5) * t325 - t479 * t290 + t482 * t324;
t243 = t482 * t258 + t479 * t312;
t242 = t479 * t258 - t482 * t312;
t241 = -pkin(5) * t318 - t479 * t288 + t482 * t303;
t240 = -pkin(1) * t301 - pkin(2) * t382 - qJ(3) * t335 - t280;
t239 = -t476 * t276 + t477 * t278;
t238 = -t476 * t275 + t477 * t277;
t237 = t477 * t276 + t476 * t278;
t236 = -t477 * t275 - t476 * t277;
t235 = -pkin(5) * t301 + t482 * t256 + t333 * t532;
t234 = t483 * t260 + t480 * t279;
t233 = t480 * t260 - t483 * t279;
t232 = t482 * t238 - t479 * t348;
t231 = -pkin(1) * t259 + pkin(2) * t366 - qJ(3) * t280;
t230 = t482 * t239 + t479 * t329;
t229 = t479 * t239 - t482 * t329;
t228 = t483 * t247 + t480 * t263;
t227 = t480 * t247 - t483 * t263;
t226 = -pkin(2) * t237 - pkin(3) * t276;
t223 = t483 * t243 + t480 * t257;
t222 = t480 * t243 - t483 * t257;
t221 = -pkin(2) * t263 - pkin(3) * t299 + t252;
t220 = -pkin(5) * t259 + (-qJ(3) * t482 + t532) * t279;
t219 = -pkin(2) * t257 - pkin(3) * t297 + t251;
t218 = -pkin(3) * t323 + pkin(6) * t225;
t217 = -qJ(3) * t263 - t476 * t253 + t477 * t274;
t216 = -qJ(3) * t257 - t476 * t245 + t477 * t262;
t215 = -pkin(6) * t276 - t224;
t214 = t483 * t230 + t480 * t237;
t213 = t480 * t230 - t483 * t237;
t212 = -pkin(3) * t329 + pkin(6) * t278 + t225;
t211 = -pkin(1) * t246 + pkin(2) * t535 - qJ(3) * t264 - t477 * t253 - t476 * t274;
t210 = t477 * t225 - t524;
t209 = t476 * t225 + t521;
t208 = -pkin(1) * t242 + pkin(2) * t312 - qJ(3) * t258 - t477 * t245 - t476 * t262;
t207 = t482 * t210 + t479 * t323;
t206 = t479 * t210 - t482 * t323;
t205 = -pkin(5) * t246 + t482 * t217 - t479 * t221;
t204 = -pkin(5) * t242 + t482 * t216 - t479 * t219;
t203 = -pkin(2) * t209 - pkin(3) * t224;
t202 = -qJ(3) * t237 - t476 * t212 + t477 * t215;
t201 = -pkin(6) * t521 - qJ(3) * t209 - t476 * t218;
t200 = t483 * t207 + t480 * t209;
t199 = t480 * t207 - t483 * t209;
t198 = -pkin(1) * t229 + pkin(2) * t329 - qJ(3) * t239 - t477 * t212 - t476 * t215;
t197 = -pkin(5) * t229 + t482 * t202 - t479 * t226;
t196 = -pkin(1) * t206 + pkin(2) * t323 + pkin(6) * t524 - qJ(3) * t210 - t477 * t218;
t195 = -pkin(5) * t206 + t482 * t201 - t479 * t203;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t449, -t450, 0, t405, 0, 0, 0, 0, 0, 0, t376, t377, t402, t347, 0, 0, 0, 0, 0, 0, t284, t294, t282, t234, 0, 0, 0, 0, 0, 0, t223, t228, t214, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t450, -t449, 0, t404, 0, 0, 0, 0, 0, 0, t374, t375, t401, t346, 0, 0, 0, 0, 0, 0, t283, t293, t281, t233, 0, 0, 0, 0, 0, 0, t222, t227, t213, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t407, t408, 0, -t360, 0, 0, 0, 0, 0, 0, t318, t325, t301, t259, 0, 0, 0, 0, 0, 0, t242, t246, t229, t206; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t450, 0, -t449, 0, t490, -t427, -t404, -pkin(4) * t404, t483 * t417 - t492, t483 * t399 - t480 * t452, t483 * t411 + t479 * t502, t483 * t416 + t492, t483 * t409 + t480 * t501, t480 * qJDD(2) + t483 * t440, -pkin(4) * t374 - t480 * t364 + t483 * t372, -pkin(4) * t375 - t480 * t365 + t483 * t373, -pkin(4) * t401 + t483 * t360, -pkin(4) * t346 - (pkin(1) * t480 - pkin(5) * t483) * t360, t483 * t342 - t480 * t370, t483 * t322 - t480 * t332, t483 * t327 - t480 * t350, t483 * t341 - t480 * t368, t483 * t328 - t480 * t351, t483 * t355 - t480 * t380, -pkin(4) * t283 + t483 * t241 - t480 * t261, -pkin(4) * t293 + t483 * t244 - t480 * t265, -pkin(4) * t281 + t483 * t235 - t480 * t240, -pkin(4) * t233 + t483 * t220 - t480 * t231, t483 * t255 - t480 * t267, t483 * t232 - t480 * t236, t483 * t249 - t480 * t270, t483 * t254 - t480 * t266, t483 * t250 - t480 * t271, t483 * t289 - t480 * t295, -pkin(4) * t222 + t483 * t204 - t480 * t208, -pkin(4) * t227 + t483 * t205 - t480 * t211, -pkin(4) * t213 + t483 * t197 - t480 * t198, -pkin(4) * t199 + t483 * t195 - t480 * t196; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t449, 0, t450, 0, t427, t490, t405, pkin(4) * t405, t480 * t417 + t491, t480 * t399 + t483 * t452, t480 * t411 - t479 * t500, t480 * t416 - t491, t480 * t409 - t482 * t500, -t483 * qJDD(2) + t480 * t440, pkin(4) * t376 + t483 * t364 + t480 * t372, pkin(4) * t377 + t483 * t365 + t480 * t373, pkin(4) * t402 + t480 * t360, pkin(4) * t347 - (-pkin(1) * t483 - pkin(5) * t480) * t360, t480 * t342 + t483 * t370, t480 * t322 + t483 * t332, t480 * t327 + t483 * t350, t480 * t341 + t483 * t368, t480 * t328 + t483 * t351, t480 * t355 + t483 * t380, pkin(4) * t284 + t480 * t241 + t483 * t261, pkin(4) * t294 + t480 * t244 + t483 * t265, pkin(4) * t282 + t480 * t235 + t483 * t240, pkin(4) * t234 + t480 * t220 + t483 * t231, t480 * t255 + t483 * t267, t480 * t232 + t483 * t236, t480 * t249 + t483 * t270, t480 * t254 + t483 * t266, t480 * t250 + t483 * t271, t480 * t289 + t483 * t295, pkin(4) * t223 + t480 * t204 + t483 * t208, pkin(4) * t228 + t480 * t205 + t483 * t211, pkin(4) * t214 + t480 * t197 + t483 * t198, pkin(4) * t200 + t480 * t195 + t483 * t196; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t455, t456, 0, 0, t488 * t479, t482 * t444 + t479 * t447, t482 * t457 + t514, -t479 * t495 + t433, t479 * t459 + t509, 0, pkin(1) * t447 + pkin(5) * t410 + t510, -pkin(1) * t444 + pkin(5) * t412 - t516, pkin(1) * t451 + pkin(5) * t448 + t361, pkin(1) * t431 + pkin(5) * t361, t479 * t371 - t496, t479 * t334 + t482 * t400, t479 * t353 + t482 * t389, t479 * t369 + t496, t479 * t354 - t387 * t482, t479 * t381 + t433, -pkin(1) * t343 + pkin(5) * t319 + t482 * t288 + t479 * t303, -pkin(1) * t349 + pkin(5) * t326 + t482 * t290 + t479 * t324, pkin(5) * t302 + t479 * t256 + (-pkin(1) - t531) * t333, pkin(5) * t260 + (-pkin(1) + t489) * t279, t479 * t269 - t497, t479 * t238 + t482 * t348, t479 * t272 + t482 * t317, t479 * t268 + t497, t479 * t273 + t482 * t313, t479 * t296 + t482 * t439, -pkin(1) * t257 + pkin(5) * t243 + t479 * t216 + t482 * t219, -pkin(1) * t263 + pkin(5) * t247 + t479 * t217 + t482 * t221, -pkin(1) * t237 + pkin(5) * t230 + t479 * t202 + t482 * t226, -pkin(1) * t209 + pkin(5) * t207 + t479 * t201 + t482 * t203;];
tauB_reg = t1;
