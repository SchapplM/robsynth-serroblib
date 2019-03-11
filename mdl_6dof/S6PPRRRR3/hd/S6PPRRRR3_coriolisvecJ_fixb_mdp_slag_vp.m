% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:55
% EndTime: 2019-03-08 19:12:05
% DurationCPUTime: 6.41s
% Computational Cost: add. (4876->452), mult. (13527->686), div. (0->0), fcn. (12131->16), ass. (0->202)
t410 = cos(pkin(6));
t397 = qJD(1) * t410 + qJD(2);
t405 = sin(pkin(7));
t418 = cos(qJ(3));
t512 = t405 * t418;
t384 = t397 * t512;
t409 = cos(pkin(7));
t407 = cos(pkin(14));
t406 = sin(pkin(6));
t492 = qJD(1) * t406;
t473 = t407 * t492;
t403 = sin(pkin(14));
t414 = sin(qJ(3));
t517 = t403 * t414;
t342 = t409 * t418 * t473 - t492 * t517 + t384;
t513 = t405 * t414;
t511 = t407 * t409;
t534 = (t403 * t418 + t414 * t511) * t406;
t343 = -qJD(1) * t534 - t397 * t513;
t404 = sin(pkin(8));
t413 = sin(qJ(4));
t516 = t404 * t413;
t398 = pkin(10) * t516;
t417 = cos(qJ(4));
t408 = cos(pkin(8));
t509 = t408 * t417;
t510 = t408 * t413;
t537 = t342 * t417 + t343 * t510 - (pkin(3) * t509 - t398) * qJD(4);
t491 = qJD(3) * t404;
t334 = pkin(10) * t491 - t343;
t337 = qJD(3) * pkin(3) + t342;
t364 = t397 * t409 - t405 * t473;
t440 = t337 * t408 + t364 * t404;
t302 = -t413 * t334 + t417 * t440;
t412 = sin(qJ(5));
t416 = cos(qJ(5));
t480 = t408 * qJD(3);
t460 = qJD(4) + t480;
t472 = t413 * t491;
t536 = -t412 * t472 + t416 * t460;
t400 = t404 ^ 2;
t518 = t400 * t417;
t535 = t413 * MDP(6) * t518 - (t413 ^ 2 - t417 ^ 2) * MDP(7) * t400;
t349 = t410 * t513 + t534;
t490 = qJD(3) * t417;
t471 = t404 * t490;
t396 = -qJD(5) + t471;
t361 = qJD(6) - t536;
t452 = pkin(4) * t413 - pkin(11) * t417;
t432 = t452 * qJD(4);
t375 = t404 * t432;
t514 = t404 * t417;
t478 = pkin(10) * t514;
t372 = t478 + (pkin(3) * t413 + pkin(11)) * t408;
t453 = -pkin(4) * t417 - pkin(11) * t413;
t373 = (-pkin(3) + t453) * t404;
t496 = t372 * t416 + t373 * t412;
t515 = t404 * t416;
t532 = qJD(5) * t496 - t343 * t515 - t375 * t416 - t412 * t537;
t483 = qJD(5) * t416;
t485 = qJD(5) * t412;
t507 = t412 * t404;
t531 = t343 * t507 - t372 * t485 + t373 * t483 + t412 * t375 - t416 * t537;
t495 = -t342 * t413 + t343 * t509 + (pkin(3) * t510 + t478) * qJD(4);
t332 = t417 * t334;
t303 = t337 * t510 + t364 * t516 + t332;
t298 = pkin(11) * t460 + t303;
t357 = t408 * t364;
t311 = t357 + (qJD(3) * t453 - t337) * t404;
t283 = t298 * t416 + t311 * t412;
t436 = t418 * t511 - t517;
t426 = t436 * t406;
t335 = (qJD(1) * t426 + t384) * qJD(3);
t336 = qJD(3) * t343;
t506 = t413 * t336;
t285 = qJD(4) * t302 + t417 * t335 + t408 * t506;
t320 = (qJD(3) * t432 - t336) * t404;
t422 = -qJD(5) * t283 - t285 * t412 + t320 * t416;
t479 = qJD(3) * qJD(4);
t467 = t404 * t479;
t455 = t413 * t467;
t273 = -pkin(5) * t455 - t422;
t367 = t412 * t460 + t416 * t472;
t530 = t361 * (pkin(5) * t367 + pkin(12) * t361) + t273;
t348 = t410 * t512 + t426;
t378 = -t405 * t406 * t407 + t409 * t410;
t439 = t348 * t408 + t378 * t404;
t529 = -t349 * t413 + t417 * t439;
t419 = qJD(3) ^ 2;
t454 = t417 * t467;
t344 = qJD(5) * t536 + t416 * t454;
t411 = sin(qJ(6));
t415 = cos(qJ(6));
t481 = qJD(6) * t415;
t474 = t344 * t415 - t396 * t481 + t411 * t455;
t482 = qJD(6) * t411;
t309 = -t367 * t482 + t474;
t527 = t309 * t411;
t521 = t367 * t411;
t339 = t415 * t396 + t521;
t526 = t339 * t361;
t341 = t367 * t415 - t396 * t411;
t525 = t341 * t361;
t523 = t536 * t396;
t522 = t367 * t396;
t520 = t396 * t412;
t519 = t396 * t416;
t437 = t412 * t454;
t345 = qJD(5) * t367 + t437;
t508 = t411 * t345;
t505 = t413 * t414;
t504 = t413 * t418;
t503 = t414 * t417;
t502 = t415 * t345;
t501 = t416 * t417;
t500 = t417 * t418;
t499 = t303 + t396 * (pkin(5) * t412 - pkin(12) * t416);
t374 = t452 * t491;
t498 = t416 * t302 + t374 * t412;
t489 = qJD(4) * t413;
t470 = t404 * t489;
t497 = -pkin(5) * t470 + t532;
t493 = MDP(5) * qJD(3);
t488 = qJD(4) * t416;
t487 = qJD(4) * t417;
t486 = qJD(5) * t411;
t484 = qJD(5) * t415;
t476 = t411 * t514;
t469 = t404 * t487;
t428 = t285 * t416 - t298 * t485 + t311 * t483 + t320 * t412;
t272 = pkin(12) * t455 + t428;
t329 = t336 * t509;
t286 = t413 * t335 - t329 + (t413 * t440 + t332) * qJD(4);
t279 = pkin(5) * t345 - pkin(12) * t344 + t286;
t465 = -t272 * t411 + t279 * t415;
t463 = t344 * t411 - t415 * t455;
t462 = t361 * t415;
t395 = -pkin(5) * t416 - pkin(12) * t412 - pkin(4);
t461 = pkin(12) * t472 - qJD(6) * t395 + t498;
t457 = t491 * t513;
t360 = (t411 * t413 + t415 * t501) * t491;
t450 = t415 * t483 - t360;
t371 = t398 + (-pkin(3) * t417 - pkin(4)) * t408;
t380 = -t408 * t416 + t413 * t507;
t381 = t408 * t412 + t413 * t515;
t326 = pkin(5) * t380 - pkin(12) * t381 + t371;
t449 = -pkin(12) * t470 - qJD(6) * t326 - t531;
t328 = -pkin(12) * t514 + t496;
t352 = -qJD(5) * t380 + t416 * t469;
t353 = qJD(5) * t381 + t412 * t469;
t448 = -pkin(5) * t353 + pkin(12) * t352 + qJD(6) * t328 - t495;
t447 = t272 * t415 + t279 * t411;
t281 = -pkin(12) * t396 + t283;
t297 = -pkin(4) * t460 - t302;
t289 = -pkin(5) * t536 - t367 * pkin(12) + t297;
t275 = t281 * t415 + t289 * t411;
t446 = t281 * t411 - t289 * t415;
t307 = t349 * t417 + t413 * t439;
t323 = -t348 * t404 + t378 * t408;
t291 = t307 * t416 + t323 * t412;
t445 = t291 * t415 - t411 * t529;
t444 = -t291 * t411 - t415 * t529;
t282 = -t298 * t412 + t311 * t416;
t443 = -t302 * t412 + t374 * t416;
t290 = t307 * t412 - t323 * t416;
t433 = t408 * t504 + t503;
t351 = t405 * t433 + t409 * t516;
t379 = -t404 * t512 + t408 * t409;
t325 = t351 * t416 + t379 * t412;
t434 = t408 * t500 - t505;
t350 = -t405 * t434 - t409 * t514;
t442 = t325 * t415 + t350 * t411;
t441 = -t325 * t411 + t350 * t415;
t324 = t351 * t412 - t379 * t416;
t438 = -t372 * t412 + t373 * t416;
t354 = t381 * t411 + t415 * t514;
t431 = -t361 * t481 - t508;
t430 = -t361 * t482 + t502;
t280 = pkin(5) * t396 - t282;
t423 = -pkin(12) * t345 + (t280 + t282) * t361;
t317 = -t337 * t404 + t357;
t421 = -qJD(4) * t440 - t317 * t491 - t335;
t359 = t411 * t416 * t471 - t415 * t472;
t355 = t381 * t415 - t476;
t347 = t349 * qJD(3);
t346 = t348 * qJD(3);
t327 = pkin(5) * t514 - t438;
t322 = t409 * t470 + (t433 * qJD(4) + (t408 * t503 + t504) * qJD(3)) * t405;
t321 = t409 * t469 + (t434 * qJD(4) + (-t408 * t505 + t500) * qJD(3)) * t405;
t316 = -qJD(6) * t476 + t352 * t411 + t381 * t481 - t415 * t470;
t315 = -qJD(6) * t354 + t352 * t415 + t411 * t470;
t310 = qJD(6) * t341 + t463;
t300 = qJD(5) * t325 + t321 * t412 - t416 * t457;
t299 = -qJD(5) * t324 + t321 * t416 + t412 * t457;
t292 = -pkin(5) * t472 - t443;
t288 = qJD(4) * t529 + t346 * t417 - t347 * t510;
t287 = qJD(4) * t307 + t346 * t413 + t347 * t509;
t277 = -qJD(5) * t290 + t288 * t416 + t347 * t507;
t276 = qJD(5) * t291 + t288 * t412 - t347 * t515;
t271 = -qJD(6) * t275 + t465;
t270 = -qJD(6) * t446 + t447;
t1 = [-t347 * qJD(3) * MDP(4) - t346 * t493 + (-t287 * qJD(4) + (-t287 * t408 + t323 * t470 - t347 * t518) * qJD(3)) * MDP(11) + (-t288 * qJD(4) + (t347 * t400 * t413 - t288 * t408 + t323 * t469) * qJD(3)) * MDP(12) + (t276 * t396 - t287 * t536 - t290 * t455 - t345 * t529) * MDP(18) + (t277 * t396 + t287 * t367 - t291 * t455 - t344 * t529) * MDP(19) + ((-qJD(6) * t445 - t277 * t411 + t287 * t415) * t361 + t444 * t345 + t276 * t339 + t290 * t310) * MDP(25) + (-(qJD(6) * t444 + t277 * t415 + t287 * t411) * t361 - t445 * t345 + t276 * t341 + t290 * t309) * MDP(26); (-t322 * t460 + t379 * t455) * MDP(11) + (-t321 * t460 + t379 * t454) * MDP(12) + (t300 * t396 - t322 * t536 - t324 * t455 + t345 * t350) * MDP(18) + (t299 * t396 + t322 * t367 - t325 * t455 + t344 * t350) * MDP(19) + ((-qJD(6) * t442 - t299 * t411 + t322 * t415) * t361 + t441 * t345 + t300 * t339 + t324 * t310) * MDP(25) + (-(qJD(6) * t441 + t299 * t415 + t322 * t411) * t361 - t442 * t345 + t300 * t341 + t324 * t309) * MDP(26) + (-MDP(5) * t418 + (-MDP(4) + (-MDP(11) * t417 + MDP(12) * t413) * t400) * t414) * t419 * t405; (-t436 * t492 + t342 - t384) * t493 + (t317 * t470 - t286 * t408 + (t336 * t417 + (-pkin(3) * t489 - t343 * t417) * qJD(3)) * t400 - t495 * t460) * MDP(11) + (t317 * t469 - t285 * t408 + (-t506 + (-pkin(3) * t487 + t343 * t413) * qJD(3)) * t400 + t537 * t460) * MDP(12) + (t344 * t381 + t352 * t367) * MDP(13) + (-t344 * t380 - t345 * t381 + t352 * t536 - t353 * t367) * MDP(14) + (-t352 * t396 + (-t344 * t417 + (qJD(3) * t381 + t367) * t489) * t404) * MDP(15) + (t353 * t396 + (t345 * t417 + (-qJD(3) * t380 + t536) * t489) * t404) * MDP(16) + (-t396 * t404 - t400 * t490) * MDP(17) * t489 + (t286 * t380 + t297 * t353 + t371 * t345 + t532 * t396 - t495 * t536 + (-t422 * t417 + (qJD(3) * t438 + t282) * t489) * t404) * MDP(18) + (t286 * t381 + t297 * t352 + t371 * t344 + t531 * t396 + t495 * t367 + (t428 * t417 + (-qJD(3) * t496 - t283) * t489) * t404) * MDP(19) + (t309 * t355 + t315 * t341) * MDP(20) + (-t309 * t354 - t310 * t355 - t315 * t339 - t316 * t341) * MDP(21) + (t309 * t380 + t315 * t361 + t341 * t353 + t345 * t355) * MDP(22) + (-t310 * t380 - t316 * t361 - t339 * t353 - t345 * t354) * MDP(23) + (t345 * t380 + t353 * t361) * MDP(24) + ((t326 * t415 - t328 * t411) * t345 + t271 * t380 - t446 * t353 + t327 * t310 + t273 * t354 + t280 * t316 + (t411 * t449 - t415 * t448) * t361 + t497 * t339) * MDP(25) + (-(t326 * t411 + t328 * t415) * t345 - t270 * t380 - t275 * t353 + t327 * t309 + t273 * t355 + t280 * t315 + (t411 * t448 + t415 * t449) * t361 + t497 * t341) * MDP(26) + 0.2e1 * t535 * t479 + (MDP(8) * t469 - MDP(9) * t470) * (qJD(4) + 0.2e1 * t480); (t303 * t460 - t334 * t487 + t413 * t421 + t329) * MDP(11) + (t302 * t460 + (qJD(4) * t334 - t336 * t408) * t413 + t421 * t417) * MDP(12) + (t344 * t412 - t367 * t519) * MDP(13) + ((t344 - t523) * t416 + (-t345 + t522) * t412) * MDP(14) + (-t396 * t483 + (t396 * t501 + (qJD(4) * t412 - t367) * t413) * t491) * MDP(15) + (t396 * t485 + (-t417 * t520 + (-t536 + t488) * t413) * t491) * MDP(16) + t396 * MDP(17) * t472 + (-pkin(4) * t345 - t286 * t416 + t443 * t396 + t303 * t536 + (pkin(11) * t519 + t297 * t412) * qJD(5) + (-t282 * t413 + (-pkin(11) * t489 - t297 * t417) * t412) * t491) * MDP(18) + (-pkin(4) * t344 + t286 * t412 - t498 * t396 - t303 * t367 + (-pkin(11) * t520 + t297 * t416) * qJD(5) + (-t297 * t501 + (-pkin(11) * t488 + t283) * t413) * t491) * MDP(19) + (t309 * t412 * t415 + (-t412 * t482 + t450) * t341) * MDP(20) + (t339 * t360 + t341 * t359 + (-t339 * t415 - t341 * t411) * t483 + (-t527 - t310 * t415 + (t339 * t411 - t341 * t415) * qJD(6)) * t412) * MDP(21) + (-t309 * t416 + t450 * t361 + (-t341 * t396 + t430) * t412) * MDP(22) + (t310 * t416 + (-t411 * t483 + t359) * t361 + (t339 * t396 + t431) * t412) * MDP(23) + (-t345 * t416 - t361 * t520) * MDP(24) + (t395 * t502 - t280 * t359 - t292 * t339 + (t411 * t461 - t415 * t499) * t361 + (t280 * t486 - t271 + (qJD(5) * t339 + t431) * pkin(11)) * t416 + (t280 * t481 + t273 * t411 + t396 * t446 + (t361 * t486 + t310) * pkin(11)) * t412) * MDP(25) + (-t395 * t508 - t280 * t360 - t292 * t341 + (t411 * t499 + t415 * t461) * t361 + (t280 * t484 + t270 + (qJD(5) * t341 - t430) * pkin(11)) * t416 + (-t280 * t482 + t273 * t415 + t396 * t275 + (t361 * t484 + t309) * pkin(11)) * t412) * MDP(26) + ((-MDP(8) * t417 + MDP(9) * t413) * t404 * t408 - t535) * t419; -t536 ^ 2 * MDP(14) + (t344 + t523) * MDP(15) + (-t437 - t522) * MDP(16) + MDP(17) * t455 + (-t283 * t396 + t422) * MDP(18) + (-t282 * t396 - t297 * t536 - t428) * MDP(19) + (t341 * t462 + t527) * MDP(20) + ((t309 - t526) * t415 + (-t310 - t525) * t411) * MDP(21) + (t361 * t462 + t508) * MDP(22) + (-t361 ^ 2 * t411 + t502) * MDP(23) + (-pkin(5) * t310 - t283 * t339 + t423 * t411 - t415 * t530) * MDP(25) + (-pkin(5) * t309 - t283 * t341 + t411 * t530 + t423 * t415) * MDP(26) + (-MDP(13) * t536 + t367 * MDP(14) - MDP(16) * qJD(5) - t297 * MDP(18) - t341 * MDP(22) + t339 * MDP(23) - t361 * MDP(24) + MDP(25) * t446 + t275 * MDP(26)) * t367; t341 * t339 * MDP(20) + (-t339 ^ 2 + t341 ^ 2) * MDP(21) + (t474 + t526) * MDP(22) + (-t463 + t525) * MDP(23) + t345 * MDP(24) + (t275 * t361 - t280 * t341 + t465) * MDP(25) + (t280 * t339 - t361 * t446 - t447) * MDP(26) + (-MDP(22) * t521 - MDP(23) * t341 - MDP(25) * t275 + MDP(26) * t446) * qJD(6);];
tauc  = t1;
