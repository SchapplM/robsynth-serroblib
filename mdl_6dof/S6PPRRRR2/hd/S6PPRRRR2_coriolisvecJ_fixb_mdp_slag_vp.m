% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:36
% EndTime: 2019-03-08 19:05:44
% DurationCPUTime: 4.68s
% Computational Cost: add. (2983->406), mult. (7871->595), div. (0->0), fcn. (6694->14), ass. (0->195)
t403 = cos(pkin(6));
t389 = qJD(1) * t403 + qJD(2);
t398 = sin(pkin(13));
t400 = sin(pkin(6));
t407 = sin(qJ(3));
t411 = cos(qJ(3));
t401 = cos(pkin(13));
t402 = cos(pkin(7));
t499 = t401 * t402;
t416 = (t398 * t411 + t407 * t499) * t400;
t399 = sin(pkin(7));
t501 = t399 * t407;
t322 = qJD(1) * t416 + t389 * t501;
t406 = sin(qJ(4));
t410 = cos(qJ(4));
t439 = pkin(4) * t406 - pkin(10) * t410;
t377 = t439 * qJD(4);
t533 = t322 - t377;
t320 = qJD(3) * pkin(9) + t322;
t484 = qJD(1) * t400;
t461 = t401 * t484;
t348 = t389 * t402 - t399 * t461;
t520 = -t406 * t320 + t348 * t410;
t297 = -qJD(4) * pkin(4) - t520;
t405 = sin(qJ(5));
t409 = cos(qJ(5));
t469 = t409 * qJD(4);
t482 = qJD(3) * t406;
t368 = t405 * t482 - t469;
t290 = pkin(5) * t368 + t297;
t408 = cos(qJ(6));
t478 = qJD(4) * t405;
t370 = t409 * t482 + t478;
t404 = sin(qJ(6));
t503 = t370 * t404;
t326 = t408 * t368 + t503;
t532 = t290 * t326;
t480 = qJD(3) * t410;
t390 = -qJD(5) + t480;
t386 = -qJD(6) + t390;
t531 = t326 * t386;
t430 = t368 * t404 - t408 * t370;
t530 = t386 * t430;
t500 = t399 * t411;
t523 = (-t398 * t407 + t411 * t499) * t400;
t529 = t403 * t500 + t523;
t485 = MDP(12) * t410;
t468 = qJD(3) * qJD(4);
t451 = t410 * t468;
t475 = qJD(5) * t405;
t455 = t406 * t475;
t467 = qJD(4) * qJD(5);
t345 = -qJD(3) * t455 + (t451 + t467) * t409;
t300 = t410 * t320 + t406 * t348;
t298 = qJD(4) * pkin(10) + t300;
t441 = t402 * t461;
t462 = t398 * t484;
t321 = -t407 * t462 + t411 * (t389 * t399 + t441);
t380 = -pkin(4) * t410 - pkin(10) * t406 - pkin(3);
t312 = qJD(3) * t380 - t321;
t508 = t312 * t405;
t278 = t298 * t409 + t508;
t527 = qJD(1) * t523 + t389 * t500;
t317 = t527 * qJD(3);
t280 = qJD(4) * t520 + t317 * t410;
t481 = qJD(3) * t407;
t460 = t399 * t481;
t479 = qJD(3) * t411;
t318 = t389 * t460 + t441 * t481 + t462 * t479;
t308 = qJD(3) * t377 + t318;
t445 = t280 * t405 - t409 * t308;
t414 = -qJD(5) * t278 - t445;
t452 = t406 * t468;
t266 = pkin(5) * t452 - pkin(11) * t345 + t414;
t474 = qJD(5) * t409;
t454 = t406 * t474;
t476 = qJD(4) * t410;
t458 = t405 * t476;
t419 = t454 + t458;
t346 = qJD(3) * t419 + t405 * t467;
t421 = t409 * t280 - t298 * t475 + t405 * t308 + t312 * t474;
t267 = -pkin(11) * t346 + t421;
t447 = t408 * t266 - t404 * t267;
t528 = t290 * t430 + t447;
t449 = MDP(24) * t482;
t526 = qJD(4) * t449 + (-t326 ^ 2 + t430 ^ 2) * MDP(21) - t326 * t430 * MDP(20);
t525 = MDP(6) * t406;
t396 = t406 ^ 2;
t524 = MDP(7) * (-t410 ^ 2 + t396);
t372 = t404 * t409 + t405 * t408;
t351 = t372 * t406;
t477 = qJD(4) * t406;
t497 = t405 * t410;
t515 = pkin(9) * t405;
t522 = -t321 * t497 + t409 * t533 - t477 * t515;
t495 = t409 * t410;
t521 = -t321 * t495 + t380 * t474 - t405 * t533;
t519 = qJD(3) * t322 - t318;
t518 = t410 * t469 - t455;
t517 = qJD(5) + qJD(6);
t443 = t345 * t404 + t408 * t346;
t289 = -qJD(6) * t430 + t443;
t516 = pkin(10) + pkin(11);
t514 = qJD(3) * pkin(3);
t277 = -t298 * t405 + t409 * t312;
t274 = -pkin(11) * t370 + t277;
t272 = -pkin(5) * t390 + t274;
t513 = t272 * t408;
t275 = -pkin(11) * t368 + t278;
t512 = t275 * t408;
t281 = t406 * t317 + t320 * t476 + t348 * t477;
t511 = t281 * t405;
t510 = t281 * t409;
t509 = t297 * t405;
t507 = t345 * t405;
t505 = t368 * t390;
t504 = t370 * t390;
t502 = t390 * t409;
t498 = t405 * t406;
t496 = t406 * t409;
t391 = pkin(9) * t495;
t428 = pkin(5) * t406 - pkin(11) * t495;
t494 = -t428 * qJD(4) - (-t391 + (pkin(11) * t406 - t380) * t405) * qJD(5) + t522;
t493 = -t419 * pkin(11) + (-t406 * t469 - t410 * t475) * pkin(9) + t521;
t374 = t439 * qJD(3);
t492 = t405 * t374 + t409 * t520;
t371 = t404 * t405 - t408 * t409;
t422 = t371 * t410;
t491 = qJD(3) * t422 - t371 * t517;
t490 = (-t480 + t517) * t372;
t487 = t405 * t380 + t391;
t473 = qJD(6) * t404;
t472 = qJD(6) * t408;
t471 = t297 * qJD(5);
t464 = t408 * t345 - t404 * t346 - t368 * t472;
t463 = qJD(5) * t516;
t459 = t399 * t479;
t456 = t390 * t475;
t453 = t405 * t480;
t448 = MDP(17) * t477;
t273 = t275 * t473;
t446 = t404 * t266 - t273;
t444 = t409 * t374 - t405 * t520;
t442 = qJD(6) * t272 + t267;
t440 = -t300 + (-t453 + t475) * pkin(5);
t384 = t516 * t409;
t438 = qJD(3) * t428 + qJD(6) * t384 + t409 * t463 + t444;
t383 = t516 * t405;
t437 = pkin(11) * t453 - qJD(6) * t383 - t405 * t463 - t492;
t269 = t272 * t404 + t512;
t332 = t403 * t501 + t416;
t354 = -t399 * t400 * t401 + t402 * t403;
t311 = t332 * t410 + t354 * t406;
t286 = -t311 * t405 - t409 * t529;
t287 = t311 * t409 - t405 * t529;
t436 = t286 * t408 - t287 * t404;
t435 = t286 * t404 + t287 * t408;
t367 = t409 * t380;
t330 = -pkin(11) * t496 + t367 + (-pkin(5) - t515) * t410;
t339 = -pkin(11) * t498 + t487;
t433 = t330 * t404 + t339 * t408;
t310 = t332 * t406 - t354 * t410;
t358 = t402 * t406 + t410 * t501;
t337 = -t358 * t405 - t409 * t500;
t425 = -t358 * t409 + t405 * t500;
t432 = t337 * t408 + t404 * t425;
t431 = t337 * t404 - t408 * t425;
t429 = qJD(3) * t396 - t390 * t410;
t412 = qJD(4) ^ 2;
t427 = pkin(9) * t412 - t519;
t319 = -t321 - t514;
t426 = qJD(4) * (t319 + t321 - t514);
t357 = -t402 * t410 + t406 * t501;
t423 = -MDP(11) * t410 + MDP(12) * t406 - MDP(4);
t288 = -t370 * t473 + t464;
t413 = qJD(3) ^ 2;
t394 = -pkin(5) * t409 - pkin(4);
t378 = (pkin(5) * t405 + pkin(9)) * t406;
t352 = t371 * t406;
t347 = pkin(5) * t419 + pkin(9) * t476;
t336 = qJD(4) * t358 + t406 * t459;
t335 = -qJD(4) * t357 + t410 * t459;
t324 = t332 * qJD(3);
t323 = t529 * qJD(3);
t304 = -t473 * t498 + (t496 * t517 + t458) * t408 + t518 * t404;
t303 = -qJD(4) * t422 - t351 * t517;
t296 = qJD(5) * t425 - t335 * t405 + t409 * t460;
t295 = qJD(5) * t337 + t335 * t409 + t405 * t460;
t285 = -qJD(4) * t310 + t323 * t410;
t284 = qJD(4) * t311 + t323 * t406;
t276 = pkin(5) * t346 + t281;
t271 = qJD(5) * t286 + t285 * t409 + t324 * t405;
t270 = -qJD(5) * t287 - t285 * t405 + t324 * t409;
t268 = -t275 * t404 + t513;
t1 = [(-t270 * t390 + t284 * t368 + t310 * t346) * MDP(18) + (t271 * t390 + t284 * t370 + t310 * t345) * MDP(19) + (-(-qJD(6) * t435 + t270 * t408 - t271 * t404) * t386 + t284 * t326 + t310 * t289) * MDP(25) + ((qJD(6) * t436 + t270 * t404 + t271 * t408) * t386 - t284 * t430 + t310 * t288) * MDP(26) + (-MDP(11) * t284 - MDP(12) * t285) * qJD(4) + (-t323 * MDP(5) + t423 * t324 + (-t529 * t485 + (-MDP(11) * t529 + t286 * MDP(18) - t287 * MDP(19) + MDP(25) * t436 - MDP(26) * t435) * t406) * qJD(4)) * qJD(3); (-t296 * t390 + t336 * t368 + t346 * t357) * MDP(18) + (t295 * t390 + t336 * t370 + t345 * t357) * MDP(19) + (-(-qJD(6) * t431 - t295 * t404 + t296 * t408) * t386 + t336 * t326 + t357 * t289) * MDP(25) + ((qJD(6) * t432 + t295 * t408 + t296 * t404) * t386 - t336 * t430 + t357 * t288) * MDP(26) + (-t336 * MDP(11) - t335 * MDP(12) + (t337 * MDP(18) + MDP(19) * t425 + MDP(25) * t432 - MDP(26) * t431) * t482) * qJD(4) + ((-MDP(11) * t406 - t485) * t411 * t468 + (-t411 * MDP(5) + t407 * t423) * t413) * t399; t519 * MDP(4) + (t321 - t527) * qJD(3) * MDP(5) + 0.2e1 * t451 * t525 - 0.2e1 * t468 * t524 + (t406 * t426 - t410 * t427) * MDP(11) + (t406 * t427 + t410 * t426) * MDP(12) + (t345 * t496 + t370 * t518) * MDP(13) + ((-t368 * t409 - t370 * t405) * t476 + (-t507 - t346 * t409 + (t368 * t405 - t370 * t409) * qJD(5)) * t406) * MDP(14) + (t390 * t455 - t345 * t410 + (t370 * t406 + t409 * t429) * qJD(4)) * MDP(15) + (t390 * t454 + t346 * t410 + (-t368 * t406 - t405 * t429) * qJD(4)) * MDP(16) + (-t390 - t480) * t448 + ((t380 * t475 + t522) * t390 + ((pkin(9) * t368 + t509) * qJD(4) + (t508 + (pkin(9) * t390 + t298) * t409) * qJD(5) + t445) * t410 + (t409 * t471 + pkin(9) * t346 + t511 - t321 * t368 + ((-pkin(9) * t497 + t367) * qJD(3) + t277) * qJD(4)) * t406) * MDP(18) + (t521 * t390 + (t297 * t469 + (qJD(4) * t370 - t456) * pkin(9) + t421) * t410 + (-t405 * t471 + pkin(9) * t345 + t510 - t321 * t370 + (-pkin(9) * t502 - qJD(3) * t487 - t278) * qJD(4)) * t406) * MDP(19) + (-t288 * t352 - t303 * t430) * MDP(20) + (-t288 * t351 + t289 * t352 - t303 * t326 + t304 * t430) * MDP(21) + (-t288 * t410 - t303 * t386 + (-qJD(3) * t352 - t430) * t477) * MDP(22) + (t289 * t410 + t304 * t386 + (-qJD(3) * t351 - t326) * t477) * MDP(23) + (-t386 - t480) * MDP(24) * t477 + (-t447 * t410 + t347 * t326 + t378 * t289 + t276 * t351 + t290 * t304 + (t404 * t493 + t408 * t494) * t386 + (t269 * t410 + t386 * t433) * qJD(6) + (-t321 * t326 + ((t330 * t408 - t339 * t404) * qJD(3) + t268) * qJD(4)) * t406) * MDP(25) + ((t442 * t408 + t446) * t410 - t347 * t430 + t378 * t288 - t276 * t352 + t290 * t303 + ((qJD(6) * t330 + t493) * t408 + (-qJD(6) * t339 - t494) * t404) * t386 + (t321 * t430 + (-qJD(3) * t433 - t269) * qJD(4)) * t406) * MDP(26) + (MDP(8) * t410 - MDP(9) * t406) * t412; (qJD(4) * t300 - t319 * t482 - t281) * MDP(11) + (-qJD(3) * t319 - t317) * t485 + (-t370 * t502 + t507) * MDP(13) + ((t345 + t505) * t409 + (-t346 + t504) * t405) * MDP(14) + (-t390 * t474 + (t390 * t495 + (-t370 + t478) * t406) * qJD(3)) * MDP(15) + (t456 + (-t390 * t497 + (t368 + t469) * t406) * qJD(3)) * MDP(16) + t390 * MDP(17) * t482 + (-pkin(4) * t346 - t510 + t444 * t390 - t300 * t368 + (pkin(10) * t502 + t509) * qJD(5) + (-t277 * t406 + (-pkin(10) * t477 - t297 * t410) * t405) * qJD(3)) * MDP(18) + (-pkin(4) * t345 + t511 - t492 * t390 - t300 * t370 + (-pkin(10) * t390 * t405 + t297 * t409) * qJD(5) + (-t297 * t495 + (-pkin(10) * t469 + t278) * t406) * qJD(3)) * MDP(19) + (t288 * t372 - t430 * t491) * MDP(20) + (-t288 * t371 - t289 * t372 - t326 * t491 + t430 * t490) * MDP(21) + (-t491 * t386 + (qJD(4) * t372 + t430) * t482) * MDP(22) + (t490 * t386 + (-qJD(4) * t371 + t326) * t482) * MDP(23) + t386 * t449 + (t276 * t371 + t394 * t289 + (t404 * t437 + t408 * t438) * t386 + t440 * t326 + t490 * t290 + ((-t383 * t408 - t384 * t404) * qJD(4) - t268) * t482) * MDP(25) + (t276 * t372 + t394 * t288 + (-t404 * t438 + t408 * t437) * t386 - t440 * t430 + t491 * t290 + (-(-t383 * t404 + t384 * t408) * qJD(4) + t269) * t482) * MDP(26) + (-t410 * t525 + t524) * t413; t370 * t368 * MDP(13) + (-t368 ^ 2 + t370 ^ 2) * MDP(14) + (t345 - t505) * MDP(15) + (-t346 - t504) * MDP(16) + qJD(3) * t448 + (-t278 * t390 - t297 * t370 + t414) * MDP(18) + (-t277 * t390 + t297 * t368 - t421) * MDP(19) + (t288 - t531) * MDP(22) + (-t289 + t530) * MDP(23) + ((-t274 * t404 - t512) * t386 - t269 * qJD(6) + (-t326 * t370 + t386 * t473 + t408 * t452) * pkin(5) + t528) * MDP(25) + (t532 + t273 + (t275 * t386 - t266) * t404 + (-t274 * t386 - t442) * t408 + (t370 * t430 + t386 * t472 - t404 * t452) * pkin(5)) * MDP(26) + t526; (t464 - t531) * MDP(22) + (-t443 + t530) * MDP(23) + (-t269 * t386 + t528) * MDP(25) + (-t408 * t267 - t268 * t386 - t446 + t532) * MDP(26) + (-MDP(22) * t503 + MDP(23) * t430 - MDP(25) * t269 - MDP(26) * t513) * qJD(6) + t526;];
tauc  = t1;
