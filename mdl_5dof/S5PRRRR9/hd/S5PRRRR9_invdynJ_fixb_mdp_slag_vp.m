% Calculate vector of inverse dynamics joint torques for
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:24
% EndTime: 2019-12-05 17:21:33
% DurationCPUTime: 5.48s
% Computational Cost: add. (2397->448), mult. (5522->636), div. (0->0), fcn. (4346->14), ass. (0->197)
t422 = sin(qJ(4));
t426 = cos(qJ(4));
t493 = t426 * qJD(3);
t423 = sin(qJ(3));
t505 = qJD(2) * t423;
t378 = t422 * t505 - t493;
t425 = cos(qJ(5));
t502 = qJD(3) * t422;
t380 = t426 * t505 + t502;
t421 = sin(qJ(5));
t534 = t380 * t421;
t319 = t378 * t425 + t534;
t427 = cos(qJ(3));
t504 = qJD(2) * t427;
t403 = -qJD(4) + t504;
t399 = -qJD(5) + t403;
t560 = t319 * t399;
t452 = t378 * t421 - t380 * t425;
t559 = t399 * t452;
t491 = qJD(2) * qJD(3);
t473 = t427 * t491;
t489 = qJDD(2) * t423;
t558 = t473 + t489;
t497 = qJD(4) * t423;
t557 = -qJD(2) * t497 + qJDD(3);
t424 = sin(qJ(2));
t419 = sin(pkin(5));
t509 = qJD(1) * t419;
t388 = qJD(2) * pkin(7) + t424 * t509;
t420 = cos(pkin(5));
t508 = qJD(1) * t420;
t401 = t423 * t508;
t347 = t388 * t427 + t401;
t341 = qJD(3) * pkin(8) + t347;
t390 = -pkin(3) * t427 - pkin(8) * t423 - pkin(2);
t428 = cos(qJ(2));
t507 = qJD(1) * t428;
t485 = t419 * t507;
t349 = qJD(2) * t390 - t485;
t307 = t341 * t426 + t349 * t422;
t293 = -pkin(9) * t378 + t307;
t495 = qJD(5) * t421;
t289 = t293 * t495;
t551 = -t388 * t423 + t427 * t508;
t340 = -qJD(3) * pkin(3) - t551;
t312 = pkin(4) * t378 + t340;
t541 = cos(pkin(10));
t468 = t541 * t424;
t418 = sin(pkin(10));
t529 = t418 * t428;
t359 = t420 * t468 + t529;
t469 = t419 * t541;
t330 = t359 * t427 - t423 * t469;
t467 = t541 * t428;
t530 = t418 * t424;
t361 = -t420 * t530 + t467;
t332 = t418 * t419 * t423 + t361 * t427;
t358 = -t420 * t467 + t530;
t360 = t420 * t529 + t468;
t527 = t419 * t427;
t366 = t420 * t423 + t424 * t527;
t417 = qJ(4) + qJ(5);
t413 = sin(t417);
t414 = cos(t417);
t526 = t419 * t428;
t556 = t312 * t319 - g(1) * (-t332 * t414 - t360 * t413) - g(2) * (-t330 * t414 - t358 * t413) - g(3) * (-t366 * t414 + t413 * t526) + t289;
t316 = qJD(4) * t493 + t422 * t557 + t426 * t558;
t411 = t427 * qJDD(2);
t550 = -t423 * t491 + t411;
t375 = qJDD(4) - t550;
t492 = qJD(1) * qJD(2);
t351 = qJDD(2) * pkin(7) + (qJDD(1) * t424 + t428 * t492) * t419;
t490 = qJDD(1) * t420;
t471 = t423 * t490;
t297 = qJDD(3) * pkin(8) + qJD(3) * t551 + t351 * t427 + t471;
t461 = pkin(3) * t423 - pkin(8) * t427;
t386 = t461 * qJD(3);
t475 = t424 * t492;
t456 = -qJDD(1) * t526 + t419 * t475;
t311 = qJD(2) * t386 + qJDD(2) * t390 + t456;
t310 = t426 * t311;
t433 = -qJD(4) * t307 - t422 * t297 + t310;
t281 = pkin(4) * t375 - pkin(9) * t316 + t433;
t317 = (qJD(3) * (qJD(4) + t504) + t489) * t422 - t557 * t426;
t496 = qJD(4) * t426;
t488 = t297 * t426 + t311 * t422 + t349 * t496;
t498 = qJD(4) * t422;
t443 = t341 * t498 - t488;
t282 = -pkin(9) * t317 - t443;
t466 = t281 * t425 - t421 * t282;
t555 = t312 * t452 - g(1) * (-t332 * t413 + t360 * t414) - g(2) * (-t330 * t413 + t358 * t414) - g(3) * (-t366 * t413 - t414 * t526) + t466;
t371 = qJDD(5) + t375;
t554 = t371 * MDP(23) + (-t319 ^ 2 + t452 ^ 2) * MDP(20) - t319 * MDP(19) * t452;
t382 = t421 * t426 + t422 * t425;
t354 = t382 * t423;
t501 = qJD(3) * t423;
t521 = t427 * t428;
t543 = pkin(7) * t422;
t553 = (-t422 * t521 + t424 * t426) * t509 - t426 * t386 - t501 * t543;
t552 = -(t422 * t424 + t426 * t521) * t509 + t422 * t386 + t390 * t496;
t549 = -t422 * t497 + t427 * t493;
t548 = qJD(4) + qJD(5);
t429 = qJD(3) ^ 2;
t460 = g(1) * t360 + g(2) * t358;
t547 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t429 + t419 * (-g(3) * t428 + t475) - t456 + t460;
t545 = -g(3) * t526 + qJD(4) * (pkin(7) * t403 + t341) + t460;
t465 = t316 * t421 + t317 * t425;
t286 = -qJD(5) * t452 + t465;
t544 = pkin(8) + pkin(9);
t542 = qJD(2) * pkin(2);
t306 = -t341 * t422 + t349 * t426;
t292 = -pkin(9) * t380 + t306;
t288 = -pkin(4) * t403 + t292;
t539 = t288 * t425;
t538 = t293 * t425;
t537 = t316 * t422;
t536 = t378 * t403;
t535 = t380 * t403;
t533 = t380 * t426;
t532 = t413 * t427;
t531 = t414 * t427;
t528 = t419 * t424;
t525 = t421 * t281;
t524 = t422 * t423;
t523 = t423 * t426;
t522 = t426 * t427;
t520 = qJDD(1) - g(3);
t404 = pkin(7) * t522;
t451 = pkin(4) * t423 - pkin(9) * t522;
t519 = -t451 * qJD(3) - (-t404 + (pkin(9) * t423 - t390) * t422) * qJD(4) + t553;
t500 = qJD(3) * t427;
t478 = t422 * t500;
t440 = t423 * t496 + t478;
t518 = -t440 * pkin(9) + (-t423 * t493 - t427 * t498) * pkin(7) + t552;
t381 = t421 * t422 - t425 * t426;
t444 = t381 * t427;
t517 = qJD(2) * t444 - t381 * t548;
t516 = (-t504 + t548) * t382;
t383 = t461 * qJD(2);
t515 = t383 * t422 + t426 * t551;
t512 = t390 * t422 + t404;
t415 = t423 ^ 2;
t511 = -t427 ^ 2 + t415;
t506 = qJD(2) * t419;
t503 = qJD(3) * t378;
t499 = qJD(4) * t403;
t494 = qJD(5) * t425;
t487 = t316 * t425 - t317 * t421 - t378 * t494;
t486 = qJD(4) * t544;
t483 = t423 * t507;
t482 = t424 * t506;
t481 = t428 * t506;
t480 = t422 * t504;
t479 = t403 * t493;
t464 = qJD(5) * t288 + t282;
t462 = -t347 + (-t480 + t498) * pkin(4);
t459 = g(1) * t361 + g(2) * t359;
t369 = t426 * t383;
t395 = t544 * t426;
t458 = qJD(2) * t451 + qJD(5) * t395 - t422 * t551 + t426 * t486 + t369;
t394 = t544 * t422;
t457 = pkin(9) * t480 - qJD(5) * t394 - t422 * t486 - t515;
t284 = t288 * t421 + t538;
t377 = t426 * t390;
t325 = -pkin(9) * t523 + t377 + (-pkin(4) - t543) * t427;
t337 = -pkin(9) * t524 + t512;
t455 = t325 * t421 + t337 * t425;
t335 = -t366 * t422 - t426 * t526;
t448 = -t366 * t426 + t422 * t526;
t454 = t335 * t425 + t421 * t448;
t453 = t335 * t421 - t425 * t448;
t365 = -t420 * t427 + t423 * t528;
t446 = t375 * t422 - t403 * t496;
t445 = t375 * t426 + t403 * t498;
t285 = -t380 * t495 + t487;
t442 = g(1) * (-t361 * t423 + t418 * t527) + g(2) * (-t359 * t423 - t427 * t469) - g(3) * t365;
t441 = qJD(3) * t401 + t351 * t423 + t388 * t500 - t427 * t490;
t438 = -g(3) * t528 - t459;
t436 = -pkin(8) * t375 - t340 * t403;
t298 = -qJDD(3) * pkin(3) + t441;
t432 = pkin(8) * t499 - t298 - t442;
t389 = -t485 - t542;
t431 = -pkin(7) * qJDD(3) + (t389 + t485 - t542) * qJD(3);
t430 = qJD(2) ^ 2;
t408 = -pkin(4) * t426 - pkin(3);
t387 = (pkin(4) * t422 + pkin(7)) * t423;
t355 = t381 * t423;
t348 = pkin(4) * t440 + pkin(7) * t500;
t334 = qJD(3) * t366 + t423 * t481;
t333 = -qJD(3) * t365 + t427 * t481;
t300 = -t495 * t524 + (t523 * t548 + t478) * t425 + t549 * t421;
t299 = -qJD(3) * t444 - t354 * t548;
t291 = qJD(4) * t335 + t333 * t426 + t422 * t482;
t290 = qJD(4) * t448 - t333 * t422 + t426 * t482;
t287 = pkin(4) * t317 + t298;
t283 = -t293 * t421 + t539;
t1 = [t520 * MDP(1) + (-qJD(3) * t334 - qJDD(3) * t365) * MDP(10) + (-qJD(3) * t333 - qJDD(3) * t366) * MDP(11) + (-t290 * t403 + t317 * t365 + t334 * t378 + t335 * t375) * MDP(17) + (t291 * t403 + t316 * t365 + t334 * t380 + t375 * t448) * MDP(18) + (t334 * t319 + t365 * t286 - (-qJD(5) * t453 + t290 * t425 - t291 * t421) * t399 + t454 * t371) * MDP(24) + (-t334 * t452 + t365 * t285 + (qJD(5) * t454 + t290 * t421 + t291 * t425) * t399 - t453 * t371) * MDP(25) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t427 + MDP(11) * t423 - MDP(3)) * t430) * t424 + (t550 * MDP(10) - MDP(11) * t558 + qJDD(2) * MDP(3) - t430 * MDP(4)) * t428) * t419; qJDD(2) * MDP(2) + (t520 * t526 + t460) * MDP(3) + (-t520 * t528 + t459) * MDP(4) + (qJDD(2) * t415 + 0.2e1 * t423 * t473) * MDP(5) + 0.2e1 * (t411 * t423 - t491 * t511) * MDP(6) + (qJDD(3) * t423 + t427 * t429) * MDP(7) + (qJDD(3) * t427 - t423 * t429) * MDP(8) + (t431 * t423 + t427 * t547) * MDP(10) + (-t423 * t547 + t431 * t427) * MDP(11) + (t316 * t523 + t380 * t549) * MDP(12) + ((-t378 * t426 - t380 * t422) * t500 + (-t537 - t317 * t426 + (t378 * t422 - t533) * qJD(4)) * t423) * MDP(13) + ((-t316 - t479) * t427 + (qJD(3) * t380 + t445) * t423) * MDP(14) + ((t403 * t502 + t317) * t427 + (-t446 - t503) * t423) * MDP(15) + (-t375 * t427 - t403 * t501) * MDP(16) + (t377 * t375 + t553 * t403 + (t390 * t499 + t438) * t422 + (pkin(7) * t503 - t310 + (-pkin(7) * t375 + qJD(3) * t340 + qJD(4) * t349 + t297) * t422 + t545 * t426) * t427 + (pkin(7) * t317 + qJD(3) * t306 + t298 * t422 + t340 * t496 - t378 * t485) * t423) * MDP(17) + (-t512 * t375 + t552 * t403 + t438 * t426 + ((pkin(7) * t380 + t340 * t426) * qJD(3) - t545 * t422 + t488) * t427 + (-t380 * t485 - t340 * t498 - t307 * qJD(3) + t298 * t426 + (t316 - t479) * pkin(7)) * t423) * MDP(18) + (-t285 * t355 - t299 * t452) * MDP(19) + (-t285 * t354 + t286 * t355 - t299 * t319 + t300 * t452) * MDP(20) + (-t285 * t427 - t299 * t399 - t355 * t371 - t452 * t501) * MDP(21) + (t286 * t427 + t300 * t399 - t319 * t501 - t354 * t371) * MDP(22) + (-t371 * t427 - t399 * t501) * MDP(23) + (t348 * t319 + t387 * t286 + t287 * t354 + t312 * t300 + (t325 * t425 - t337 * t421) * t371 - t466 * t427 + t283 * t501 - g(1) * (-t360 * t531 + t361 * t413) - g(2) * (-t358 * t531 + t359 * t413) + (t421 * t518 + t425 * t519) * t399 + (t284 * t427 + t399 * t455) * qJD(5) + (-t319 * t483 - g(3) * (t413 * t424 + t414 * t521)) * t419) * MDP(24) + (-t348 * t452 + t387 * t285 - t287 * t355 + t312 * t299 - t455 * t371 + (t464 * t425 - t289 + t525) * t427 - t284 * t501 - g(1) * (t360 * t532 + t361 * t414) - g(2) * (t358 * t532 + t359 * t414) + ((qJD(5) * t325 + t518) * t425 + (-qJD(5) * t337 - t519) * t421) * t399 + (t452 * t483 - g(3) * (-t413 * t521 + t414 * t424)) * t419) * MDP(25); MDP(7) * t489 + MDP(8) * t411 + qJDD(3) * MDP(9) + (qJD(3) * t347 - t441 - t442) * MDP(10) + (-t471 + g(1) * t332 + g(2) * t330 + g(3) * t366 + (-qJD(2) * t389 - t351) * t427) * MDP(11) + (-t403 * t533 + t537) * MDP(12) + ((t316 + t536) * t426 + (-t317 + t535) * t422) * MDP(13) + ((-t380 * t423 + t403 * t522) * qJD(2) + t446) * MDP(14) + ((-t403 * t422 * t427 + t378 * t423) * qJD(2) + t445) * MDP(15) + (-pkin(3) * t317 - t347 * t378 + t369 * t403 + (-t403 * t551 + t436) * t422 + t432 * t426) * MDP(17) + (-pkin(3) * t316 - t347 * t380 - t403 * t515 - t422 * t432 + t426 * t436) * MDP(18) + (t285 * t382 - t452 * t517) * MDP(19) + (-t285 * t381 - t286 * t382 - t319 * t517 + t452 * t516) * MDP(20) + (t371 * t382 - t399 * t517) * MDP(21) + (-t371 * t381 + t399 * t516) * MDP(22) + (t408 * t286 + t287 * t381 + (-t394 * t425 - t395 * t421) * t371 + (t421 * t457 + t425 * t458) * t399 + t462 * t319 + t516 * t312 - t442 * t414) * MDP(24) + (t408 * t285 + t287 * t382 - (-t394 * t421 + t395 * t425) * t371 + (-t421 * t458 + t425 * t457) * t399 - t462 * t452 + t517 * t312 + t442 * t413) * MDP(25) + (-MDP(10) * t389 + t403 * MDP(16) - t306 * MDP(17) + MDP(18) * t307 + MDP(21) * t452 + t319 * MDP(22) + t399 * MDP(23) - t283 * MDP(24) + t284 * MDP(25)) * t505 + (-MDP(5) * t423 * t427 + MDP(6) * t511) * t430; t380 * t378 * MDP(12) + (-t378 ^ 2 + t380 ^ 2) * MDP(13) + (t316 - t536) * MDP(14) + (-t317 - t535) * MDP(15) + t375 * MDP(16) + (-t307 * t403 - t340 * t380 - g(1) * (-t332 * t422 + t360 * t426) - g(2) * (-t330 * t422 + t358 * t426) - g(3) * t335 + t433) * MDP(17) + (-t306 * t403 + t340 * t378 - g(1) * (-t332 * t426 - t360 * t422) - g(2) * (-t330 * t426 - t358 * t422) - g(3) * t448 + t443) * MDP(18) + (t285 - t560) * MDP(21) + (-t286 + t559) * MDP(22) + ((-t292 * t421 - t538) * t399 - t284 * qJD(5) + (-t319 * t380 + t371 * t425 + t399 * t495) * pkin(4) + t555) * MDP(24) + ((t293 * t399 - t281) * t421 + (-t292 * t399 - t464) * t425 + (-t371 * t421 + t380 * t452 + t399 * t494) * pkin(4) + t556) * MDP(25) + t554; (t487 - t560) * MDP(21) + (-t465 + t559) * MDP(22) + (-t284 * t399 + t555) * MDP(24) + (-t425 * t282 - t283 * t399 - t525 + t556) * MDP(25) + (-MDP(21) * t534 + MDP(22) * t452 - MDP(24) * t284 - MDP(25) * t539) * qJD(5) + t554;];
tau = t1;
