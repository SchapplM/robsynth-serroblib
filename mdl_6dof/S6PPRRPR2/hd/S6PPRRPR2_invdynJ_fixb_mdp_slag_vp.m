% Calculate vector of inverse dynamics joint torques for
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:51:04
% EndTime: 2019-03-08 18:51:11
% DurationCPUTime: 4.84s
% Computational Cost: add. (2557->423), mult. (6511->587), div. (0->0), fcn. (6063->14), ass. (0->200)
t539 = cos(pkin(6));
t387 = t539 * qJDD(1) + qJDD(2);
t390 = qJD(1) * t539 + qJD(2);
t402 = sin(pkin(7));
t405 = cos(pkin(7));
t411 = cos(qJ(3));
t404 = cos(pkin(12));
t403 = sin(pkin(6));
t515 = qJD(1) * t403;
t489 = t404 * t515;
t473 = t405 * t489;
t498 = qJDD(1) * t403;
t481 = t404 * t498;
t408 = sin(qJ(3));
t513 = qJD(3) * t408;
t488 = t402 * t513;
t401 = sin(pkin(12));
t490 = t401 * t515;
t524 = t401 * t408;
t494 = t403 * t524;
t511 = qJD(3) * t411;
t421 = -t411 * (t387 * t402 + t405 * t481) + qJDD(1) * t494 + t390 * t488 + t473 * t513 + t490 * t511;
t407 = sin(qJ(4));
t499 = qJD(3) * qJD(4);
t484 = t407 * t499;
t419 = pkin(4) * t484 + t421;
t503 = qJD(5) * t407;
t410 = cos(qJ(4));
t504 = qJD(4) * t410;
t448 = -qJ(5) * t504 - t503;
t480 = -qJ(5) * t407 - pkin(3);
t457 = pkin(4) * t410 - t480;
t563 = qJDD(3) * t457;
t294 = qJD(3) * t448 + t419 - t563;
t538 = cos(pkin(11));
t467 = t539 * t538;
t537 = sin(pkin(11));
t357 = t401 * t467 + t404 * t537;
t427 = t401 * t537 - t404 * t467;
t478 = t403 * t538;
t558 = t402 * t478 + t405 * t427;
t317 = t357 * t408 + t411 * t558;
t466 = t539 * t537;
t358 = -t401 * t466 + t404 * t538;
t428 = t401 * t538 + t404 * t466;
t477 = t403 * t537;
t559 = -t402 * t477 + t428 * t405;
t319 = t358 * t408 + t411 * t559;
t479 = t402 * t539;
t469 = t411 * t479;
t521 = t404 * t405;
t492 = t411 * t521;
t339 = -t403 * t492 - t469 + t494;
t445 = g(1) * t319 + g(2) * t317 + g(3) * t339;
t568 = -t294 + t445;
t452 = t401 * t411 + t408 * t521;
t440 = t452 * t403;
t523 = t402 * t408;
t331 = qJD(1) * t440 + t390 * t523;
t329 = qJD(3) * pkin(9) + t331;
t352 = t390 * t405 - t402 * t489;
t519 = -t407 * t329 + t410 * t352;
t564 = -qJD(5) + t519;
t567 = qJD(3) * t331 - t421 + t445;
t311 = -qJD(4) * pkin(4) - t564;
t496 = qJDD(3) * t410;
t566 = t484 - t496;
t512 = qJD(3) * t410;
t565 = qJD(3) * t457;
t562 = MDP(11) - MDP(14);
t561 = MDP(12) - MDP(15);
t506 = qJD(4) * t407;
t560 = -pkin(4) * t506 + t331;
t483 = t410 * t499;
t497 = qJDD(3) * t407;
t447 = t483 + t497;
t365 = qJDD(6) + t447;
t409 = cos(qJ(6));
t359 = t409 * t365;
t514 = qJD(3) * t407;
t392 = qJD(6) + t514;
t406 = sin(qJ(6));
t502 = qJD(6) * t406;
t557 = -t392 * t502 + t359;
t314 = t410 * t329 + t407 * t352;
t312 = -qJD(4) * qJ(5) - t314;
t500 = pkin(5) * t514 - t564;
t453 = t492 - t524;
t533 = qJDD(3) * pkin(9);
t305 = t533 + (t387 * t408 + t390 * t511) * t402 + (qJD(1) * qJD(3) * t453 + qJDD(1) * t452) * t403;
t556 = t407 * t305 + t329 * t504 + t352 * t506;
t318 = t357 * t411 - t408 * t558;
t416 = t402 * t427 - t405 * t478;
t297 = t318 * t407 - t410 * t416;
t320 = t358 * t411 - t408 * t559;
t417 = t402 * t428 + t405 * t477;
t299 = t320 * t407 - t410 * t417;
t340 = t408 * t479 + t440;
t356 = -t403 * t404 * t402 + t405 * t539;
t321 = t340 * t407 - t356 * t410;
t555 = -g(1) * t299 - g(2) * t297 - g(3) * t321;
t393 = pkin(5) * t512;
t308 = -t312 + t393;
t412 = -pkin(4) - pkin(10);
t554 = t412 * t365 + (t308 - t393 - t314) * t392;
t330 = -t408 * t490 + t411 * (t390 * t402 + t473);
t518 = t448 - t560;
t413 = qJD(4) ^ 2;
t547 = pkin(9) * t413;
t551 = qJD(3) * t518 + t547 - t563 - t568;
t548 = pkin(5) + pkin(9);
t540 = qJD(3) * pkin(3);
t536 = pkin(9) * qJDD(4);
t535 = qJ(5) * t410;
t532 = qJDD(4) * pkin(4);
t507 = qJD(4) * t406;
t366 = t409 * t512 + t507;
t456 = t409 * qJDD(4) + t566 * t406;
t332 = -qJD(6) * t366 + t456;
t531 = t332 * t409;
t350 = t387 * t405 - t402 * t481;
t530 = t350 * t410;
t529 = t365 * t406;
t528 = t366 * t392;
t485 = t406 * t512;
t505 = qJD(4) * t409;
t368 = -t485 + t505;
t527 = t368 * t392;
t526 = t392 * t406;
t525 = t392 * t409;
t522 = t402 * t411;
t465 = pkin(10) * t407 - t535;
t430 = qJD(4) * t465 - t503;
t520 = -t430 + t560;
t399 = t407 ^ 2;
t400 = t410 ^ 2;
t517 = t399 - t400;
t509 = qJD(4) * t366;
t508 = qJD(4) * t368;
t501 = qJD(6) * t410;
t495 = qJDD(4) * qJ(5);
t382 = t548 * t410;
t493 = t407 * t523;
t414 = qJD(3) ^ 2;
t491 = t407 * t410 * t414;
t487 = t402 * t511;
t482 = t411 * t499;
t475 = -t410 * t305 + t329 * t506 - t407 * t350 - t352 * t504;
t464 = -qJD(6) * t485 + qJDD(4) * t406 - t409 * t484;
t307 = qJD(4) * t412 + t500;
t364 = t410 * t412 + t480;
t315 = qJD(3) * t364 - t330;
t291 = t307 * t409 - t315 * t406;
t292 = t307 * t406 + t315 * t409;
t301 = t321 * t409 - t339 * t406;
t302 = t321 * t406 + t339 * t409;
t322 = t340 * t410 + t356 * t407;
t461 = qJDD(3) * t411 - t408 * t414;
t360 = -t405 * t410 + t493;
t455 = -t360 * t406 + t409 * t522;
t454 = t360 * t409 + t406 * t522;
t361 = t405 * t407 + t410 * t523;
t449 = -qJD(6) * t525 - t529;
t298 = t318 * t410 + t407 * t416;
t300 = t320 * t410 + t407 * t417;
t446 = -g(1) * t300 - g(2) * t298 - g(3) * t322;
t444 = g(1) * t320 + g(2) * t318 + g(3) * t340;
t443 = qJDD(5) - t530 + t556;
t293 = qJD(3) * t430 + qJDD(3) * t364 + t419;
t439 = -t293 + t445;
t328 = -t330 - t540;
t438 = -t536 + (t328 + t330 - t540) * qJD(4);
t437 = -g(1) * t477 + g(2) * t478 - g(3) * t539;
t381 = t548 * t407;
t434 = -t381 * t365 + t444;
t290 = t443 - t532;
t323 = -t330 - t565;
t433 = t536 + (-t323 - t330 + t565) * qJD(4);
t431 = qJD(4) * qJD(5) - t475 + t495;
t426 = t446 - t475;
t423 = -qJD(4) * t314 + t555 + t556;
t288 = -t566 * pkin(5) + t431;
t395 = pkin(4) * t514;
t420 = t288 + (qJD(3) * t465 - qJD(6) * t412 + t395) * t392 + t446;
t418 = 0.2e1 * qJDD(3) * pkin(3) - t547 + t567;
t415 = t431 * t410 + t290 * t407 + (t311 * t410 + t312 * t407) * qJD(4) - t444;
t375 = qJD(4) * t382;
t374 = t548 * t506;
t369 = -qJ(5) * t512 + t395;
t342 = qJD(4) * t361 + t407 * t487;
t341 = qJD(4) * t493 - t405 * t504 - t410 * t487;
t336 = t340 * qJD(3);
t335 = (t403 * t453 + t469) * qJD(3);
t333 = (qJD(4) * qJD(6) + t496) * t409 + t464;
t316 = t323 * t514;
t296 = -t340 * t506 + (qJD(4) * t356 + t335) * t410;
t295 = qJD(4) * t322 + t335 * t407;
t287 = t447 * pkin(5) + qJDD(4) * t412 + t443;
t286 = t409 * t287;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t387 * t539 - g(3) + (t401 ^ 2 + t404 ^ 2) * t403 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(3) * t336 - qJDD(3) * t339) * MDP(4) + (-qJD(3) * t335 - qJDD(3) * t340) * MDP(5) + ((t321 * t407 + t322 * t410) * qJDD(3) + (t295 * t407 + t296 * t410 + (t321 * t410 - t322 * t407) * qJD(4)) * qJD(3)) * MDP(13) + (t290 * t321 + t294 * t339 + t295 * t311 - t296 * t312 + t322 * t431 + t323 * t336 - g(3)) * MDP(16) + ((-qJD(6) * t302 + t295 * t409 - t336 * t406) * t392 + t301 * t365 + t296 * t366 + t322 * t333) * MDP(22) + (-(qJD(6) * t301 + t295 * t406 + t336 * t409) * t392 - t302 * t365 + t296 * t368 + t322 * t332) * MDP(23) + t561 * (qJD(3) * (t336 * t407 + t339 * t504) - qJD(4) * t296 - qJDD(4) * t322 + t339 * t497) + t562 * (qJD(3) * (-t336 * t410 + t339 * t506) - qJD(4) * t295 - qJDD(4) * t321 - t339 * t496); (t437 + t387) * MDP(2) + ((t360 * t407 + t361 * t410) * qJDD(3) + (-t341 * t410 + t342 * t407 + (t360 * t410 - t361 * t407) * qJD(4)) * qJD(3)) * MDP(13) + (t290 * t360 + t311 * t342 + t312 * t341 + t361 * t431 + t437) * MDP(16) + ((qJD(6) * t455 + t342 * t409 - t406 * t488) * t392 + t454 * t365 - t341 * t366 + t361 * t333) * MDP(22) + (-(qJD(6) * t454 + t342 * t406 + t409 * t488) * t392 + t455 * t365 - t341 * t368 + t361 * t332) * MDP(23) - t561 * (-qJD(4) * t341 + qJDD(4) * t361) + t562 * (-qJD(4) * t342 - qJDD(4) * t360) + (t461 * MDP(4) + (-qJDD(3) * t408 - t411 * t414) * MDP(5) + (-t294 * t411 + t323 * t513) * MDP(16) - t561 * (t407 * t461 + t410 * t482) + t562 * (-t407 * t482 + t410 * t461)) * t402; qJDD(3) * MDP(3) + t567 * MDP(4) + (-t387 * t523 - t452 * t498 + (-t390 * t522 - t453 * t515 + t330) * qJD(3) + t444) * MDP(5) + (qJDD(3) * t399 + 0.2e1 * t407 * t483) * MDP(6) + 0.2e1 * (t407 * t496 - t499 * t517) * MDP(7) + (qJDD(4) * t407 + t410 * t413) * MDP(8) + (qJDD(4) * t410 - t407 * t413) * MDP(9) + (t407 * t438 + t410 * t418) * MDP(11) + (-t407 * t418 + t410 * t438) * MDP(12) + (t415 + (-qJD(3) * t330 + t533) * (t399 + t400)) * MDP(13) + (t433 * t407 + t410 * t551) * MDP(14) + (-t407 * t551 + t433 * t410) * MDP(15) + ((-t311 * t407 + t312 * t410) * t330 + t518 * t323 + t415 * pkin(9) + t568 * t457) * MDP(16) + (-t332 * t406 * t410 + (t406 * t506 - t409 * t501) * t368) * MDP(17) + ((-t366 * t406 + t368 * t409) * t506 + (-t531 + t333 * t406 + (t366 * t409 + t368 * t406) * qJD(6)) * t410) * MDP(18) + ((t392 * t507 + t332) * t407 + (t449 + t508) * t410) * MDP(19) + ((t392 * t505 - t333) * t407 + (-t509 - t557) * t410) * MDP(20) + (t365 * t407 + t392 * t504) * MDP(21) + (-t364 * t529 + t382 * t333 - t374 * t366 - t434 * t409 + (-t308 * t505 + t406 * t439 + t286) * t407 + ((-t330 * t407 + t375) * t409 + t520 * t406) * t392 + ((-t364 * t409 - t381 * t406) * t392 - t292 * t407) * qJD(6) + (qJD(4) * t291 + t288 * t409 - t308 * t502 - t330 * t366) * t410) * MDP(22) + (-t292 * t504 + t382 * t332 + (-t330 * t410 - t374) * t368 + (-t308 * t501 - t364 * t365 + (-qJD(6) * t381 + t520) * t392 + (-qJD(6) * t307 + t439) * t407) * t409 + (-(-qJD(6) * t364 + t375) * t392 - t288 * t410 + (qJD(4) * t308 + qJD(6) * t315 + t330 * t392 - t287) * t407 + t434) * t406) * MDP(23); -MDP(6) * t491 + t517 * t414 * MDP(7) + MDP(8) * t497 + MDP(9) * t496 + qJDD(4) * MDP(10) + (-t328 * t514 - t423 + t530) * MDP(11) + (qJD(4) * t519 - t328 * t512 - t426) * MDP(12) + (-pkin(4) * t407 + t535) * qJDD(3) * MDP(13) + (-0.2e1 * t532 + qJDD(5) + t316 + (-qJD(3) * t369 - t350) * t410 + t423) * MDP(14) + (0.2e1 * t495 + (0.2e1 * qJD(5) - t519) * qJD(4) + (t323 * t410 + t369 * t407) * qJD(3) + t426) * MDP(15) + (t431 * qJ(5) - t290 * pkin(4) - t323 * t369 - t311 * t314 - g(1) * (-pkin(4) * t299 + qJ(5) * t300) - g(2) * (-pkin(4) * t297 + qJ(5) * t298) - g(3) * (-pkin(4) * t321 + qJ(5) * t322) + t564 * t312) * MDP(16) + (-t368 * t526 + t531) * MDP(17) + ((-t333 - t527) * t409 + (-t332 + t528) * t406) * MDP(18) + ((-t368 * t410 - t407 * t526) * qJD(3) + t557) * MDP(19) + ((t366 * t410 - t407 * t525) * qJD(3) + t449) * MDP(20) - t392 * MDP(21) * t512 + (qJ(5) * t333 - t291 * t512 + t500 * t366 + t420 * t406 + t409 * t554) * MDP(22) + (qJ(5) * t332 + t292 * t512 + t500 * t368 - t406 * t554 + t420 * t409) * MDP(23); MDP(13) * t497 + (qJDD(4) + t491) * MDP(14) + (-t399 * t414 - t413) * MDP(15) + (qJD(4) * t312 + t290 + t316 + t555) * MDP(16) + (t359 - t509) * MDP(22) + (-t508 - t529) * MDP(23) + (-MDP(22) * t526 - MDP(23) * t525) * t392; t368 * t366 * MDP(17) + (-t366 ^ 2 + t368 ^ 2) * MDP(18) + (t456 + t528) * MDP(19) + (-t409 * t496 - t464 + t527) * MDP(20) + t365 * MDP(21) + (-t406 * t293 + t286 + t292 * t392 - t308 * t368 - g(1) * (t299 * t409 - t319 * t406) - g(2) * (t297 * t409 - t317 * t406) - g(3) * t301) * MDP(22) + (-t409 * t293 - t406 * t287 + t291 * t392 + t308 * t366 - g(1) * (-t299 * t406 - t319 * t409) - g(2) * (-t297 * t406 - t317 * t409) + g(3) * t302) * MDP(23) + (-MDP(19) * t366 - MDP(20) * t505 - MDP(22) * t292 - MDP(23) * t291) * qJD(6);];
tau  = t1;
