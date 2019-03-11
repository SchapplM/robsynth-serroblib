% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:21:07
% EndTime: 2019-03-08 20:21:13
% DurationCPUTime: 4.76s
% Computational Cost: add. (3102->474), mult. (6302->607), div. (0->0), fcn. (4617->10), ass. (0->201)
t417 = sin(qJ(4));
t420 = cos(qJ(4));
t486 = qJDD(2) * t420;
t491 = qJD(2) * qJD(4);
t564 = -t417 * t491 + t486;
t404 = qJD(2) * t417 + qJD(5);
t464 = pkin(4) * t420 + pkin(9) * t417;
t377 = qJD(4) * t464 + qJD(3);
t419 = cos(qJ(5));
t413 = sin(pkin(6));
t507 = qJD(1) * t413;
t421 = cos(qJ(2));
t416 = sin(qJ(5));
t418 = sin(qJ(2));
t521 = t416 * t418;
t560 = t417 * t521 - t419 * t421;
t563 = -t377 * t419 - t560 * t507;
t518 = t418 * t419;
t350 = (t416 * t421 + t417 * t518) * t413;
t463 = pkin(4) * t417 - pkin(9) * t420;
t385 = qJ(3) + t463;
t422 = -pkin(2) - pkin(8);
t493 = t419 * qJD(4);
t494 = qJD(5) * t419;
t562 = -t420 * t422 * t493 + qJD(1) * t350 - t416 * t377 - t385 * t494;
t519 = t417 * t419;
t511 = t416 * t385 + t422 * t519;
t476 = t421 * t507;
t459 = qJD(3) - t476;
t372 = qJD(2) * t422 + t459;
t415 = cos(pkin(6));
t525 = t415 * t417;
t561 = -qJD(1) * t525 + t372 * t420;
t499 = qJD(4) * t416;
t501 = qJD(2) * t420;
t381 = t419 * t501 + t499;
t496 = qJD(5) * t381;
t322 = -t419 * qJDD(4) + t564 * t416 + t496;
t412 = sin(pkin(10));
t414 = cos(pkin(10));
t524 = t415 * t418;
t367 = t412 * t421 + t414 * t524;
t369 = -t412 * t524 + t414 * t421;
t559 = -g(1) * t369 - g(2) * t367;
t484 = MDP(20) + MDP(22);
t483 = MDP(21) - MDP(24);
t523 = t415 * t421;
t366 = t412 * t418 - t414 * t523;
t368 = t412 * t523 + t414 * t418;
t526 = t413 * t421;
t370 = t420 * t526 + t525;
t529 = t413 * t417;
t443 = g(3) * t370 - g(2) * (t366 * t420 + t414 * t529) - g(1) * (t368 * t420 - t412 * t529);
t325 = -qJD(4) * pkin(4) - t561;
t379 = t416 * t501 - t493;
t301 = pkin(5) * t379 - qJ(6) * t381 + t325;
t472 = t420 * t491;
t487 = qJDD(2) * t417;
t446 = t472 + t487;
t376 = qJDD(5) + t446;
t553 = pkin(9) * t376;
t558 = t301 * t404 - t553;
t516 = qJDD(1) - g(3);
t528 = t413 * t418;
t557 = -t516 * t528 - t559;
t495 = qJD(5) * t416;
t440 = t417 * t493 + t420 * t495;
t321 = qJD(2) * t440 - qJD(5) * t493 - t416 * qJDD(4) - t419 * t486;
t477 = t418 * t507;
t505 = qJD(2) * qJ(3);
t384 = t477 + t505;
t556 = qJD(4) * (-t384 + t477 - t505) - qJDD(4) * t422;
t555 = t381 ^ 2;
t554 = pkin(5) * t376;
t547 = pkin(9) * qJD(5);
t546 = qJ(6) * t376;
t545 = qJDD(2) * pkin(2);
t506 = qJD(1) * t420;
t399 = t415 * t506;
t339 = t417 * t372 + t399;
t326 = qJD(4) * pkin(9) + t339;
t347 = qJD(2) * t385 + t477;
t300 = t326 * t419 + t347 * t416;
t297 = qJ(6) * t404 + t300;
t544 = t297 * t404;
t543 = t300 * t404;
t542 = t321 * t416;
t541 = t322 * t419;
t538 = t379 * t381;
t537 = t379 * t404;
t536 = t379 * t416;
t535 = t379 * t419;
t534 = t379 * t420;
t533 = t381 * t404;
t532 = t381 * t416;
t531 = t381 * t419;
t530 = t385 * t419;
t527 = t413 * t420;
t522 = t416 * t417;
t520 = t416 * t422;
t497 = qJD(4) * t420;
t515 = qJ(6) * t497 + (-t422 * t495 + qJD(6)) * t417 - t562;
t468 = -pkin(5) + t520;
t514 = t511 * qJD(5) + t468 * t497 + t563;
t383 = t464 * qJD(2);
t513 = t416 * t383 + t419 * t561;
t460 = pkin(5) * t416 - qJ(6) * t419;
t512 = -qJD(6) * t416 + t404 * t460 - t339;
t510 = pkin(2) * t526 + qJ(3) * t528;
t411 = t420 ^ 2;
t509 = t417 ^ 2 - t411;
t423 = qJD(4) ^ 2;
t424 = qJD(2) ^ 2;
t508 = -t423 - t424;
t504 = qJD(2) * t384;
t503 = qJD(2) * t413;
t500 = qJD(4) * t379;
t498 = qJD(4) * t417;
t299 = -t326 * t416 + t347 * t419;
t492 = qJD(6) - t299;
t490 = qJDD(1) * t413;
t489 = qJDD(1) * t415;
t488 = qJDD(2) * qJ(3);
t475 = t418 * t503;
t390 = qJD(1) * t475;
t470 = t421 * t490;
t450 = qJDD(3) + t390 - t470;
t342 = qJDD(2) * t422 + t450;
t469 = t420 * t489;
t295 = qJDD(4) * pkin(9) + qJD(4) * t561 + t342 * t417 + t469;
t471 = t418 * t490;
t313 = t471 + t385 * qJDD(2) + (t377 + t476) * qJD(2);
t481 = -t419 * t295 - t416 * t313 - t347 * t494;
t480 = -qJD(4) * t399 - t372 * t498 - t417 * t489;
t474 = t421 * t503;
t467 = -t342 + t504;
t466 = -t416 * t295 + t419 * t313 - t326 * t494 - t347 * t495;
t465 = t381 * t477;
t461 = pkin(5) * t419 + qJ(6) * t416;
t445 = -t326 * t495 - t481;
t287 = qJD(6) * t404 + t445 + t546;
t288 = qJDD(6) - t466 - t554;
t458 = t287 * t419 + t288 * t416;
t294 = -pkin(5) * t404 + t492;
t457 = t294 * t419 - t297 * t416;
t456 = t294 * t416 + t297 * t419;
t455 = pkin(4) + t461;
t454 = -g(1) * t368 - g(2) * t366 + g(3) * t526;
t453 = -t422 + t460;
t371 = t415 * t420 - t417 * t526;
t451 = -t371 * t416 + t413 * t518;
t335 = t371 * t419 + t413 * t521;
t448 = t376 * t416 + t404 * t494;
t447 = t376 * t419 - t404 * t495;
t316 = t366 * t419 + t367 * t522;
t318 = t368 * t419 + t369 * t522;
t349 = t560 * t413;
t444 = g(1) * t318 + g(2) * t316 + g(3) * t349;
t329 = t368 * t417 + t412 * t527;
t331 = -t366 * t417 + t414 * t527;
t442 = g(1) * t329 - g(2) * t331 + g(3) * t371;
t438 = g(3) * t528 - t559;
t317 = -t366 * t416 + t367 * t519;
t319 = -t368 * t416 + t369 * t519;
t437 = -g(1) * t319 - g(2) * t317 - g(3) * t350 + t477 * t534;
t296 = -qJDD(4) * pkin(4) - t342 * t420 - t480;
t436 = t325 * t404 - t553;
t435 = -t416 * t484 - t419 * t483;
t434 = -t454 + t470;
t432 = -t404 * t547 + t443;
t431 = qJDD(3) - t434;
t307 = t329 * t416 - t369 * t419;
t309 = -t331 * t416 - t367 * t419;
t430 = g(1) * t307 + g(2) * t309 - g(3) * t451 + t466;
t289 = pkin(5) * t322 + qJ(6) * t321 - qJD(6) * t381 + t296;
t429 = -t289 + t432;
t308 = t329 * t419 + t369 * t416;
t310 = -t331 * t419 + t367 * t416;
t428 = -g(1) * t308 - g(2) * t310 - g(3) * t335 + t445;
t427 = t301 * t381 + qJDD(6) - t430;
t343 = t471 + t488 + (qJD(3) + t476) * qJD(2);
t426 = qJD(2) * t459 - t422 * t423 + t343 - t438 + t488;
t425 = (t531 + t536) * MDP(23) + t457 * MDP(25) + (t416 * t483 - t419 * t484) * t404;
t408 = qJDD(4) * t420;
t378 = -qJD(2) * pkin(2) + t459;
t360 = t368 * pkin(2);
t359 = t366 * pkin(2);
t351 = t453 * t420;
t348 = t450 - t545;
t341 = t417 * t468 - t530;
t340 = qJ(6) * t417 + t511;
t333 = qJD(4) * t371 - t420 * t475;
t332 = -qJD(4) * t370 + t417 * t475;
t327 = pkin(5) * t381 + qJ(6) * t379;
t314 = (qJD(5) * t461 - qJD(6) * t419) * t420 - t453 * t498;
t312 = -pkin(5) * t501 - t383 * t419 + t416 * t561;
t311 = qJ(6) * t501 + t513;
t303 = -t321 + t537;
t291 = qJD(5) * t451 + t332 * t419 + t416 * t474;
t290 = qJD(5) * t335 + t332 * t416 - t419 * t474;
t1 = [t516 * MDP(1) + (qJDD(1) * t415 ^ 2 - g(3)) * MDP(7) + (-qJD(4) * t333 - qJDD(4) * t370) * MDP(13) + (-qJD(4) * t332 - qJDD(4) * t371) * MDP(14) + (t290 * t381 - t291 * t379 + t321 * t451 - t322 * t335) * MDP(23) + (t287 * t335 - t288 * t451 + t289 * t370 + t290 * t294 + t291 * t297 + t301 * t333 - g(3)) * MDP(25) + t484 * (-t290 * t404 + t370 * t322 + t333 * t379 + t376 * t451) - t483 * (t291 * t404 + t321 * t370 - t333 * t381 + t335 * t376) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t418 + t421 * t424) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t421 + t418 * t424) + ((-t348 + t504) * MDP(7) + (MDP(13) * t417 + MDP(14) * t420) * t424) * t421 + ((qJD(2) * t378 + t343) * MDP(7) + t446 * MDP(13) + t564 * MDP(14)) * t418) * t413; qJDD(2) * MDP(2) + t434 * MDP(3) + t557 * MDP(4) + (t431 - 0.2e1 * t545) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t488 - t557) * MDP(6) + (t343 * qJ(3) + t384 * qJD(3) - t348 * pkin(2) - g(1) * (qJ(3) * t369 - t360) - g(2) * (qJ(3) * t367 - t359) - g(3) * t510 + (-t378 * t418 - t384 * t421) * t507) * MDP(7) + (qJDD(2) * t411 - 0.2e1 * t417 * t472) * MDP(8) + 0.2e1 * (-t417 * t486 + t491 * t509) * MDP(9) + (-t417 * t423 + t408) * MDP(10) + (-qJDD(4) * t417 - t420 * t423) * MDP(11) + (t426 * t417 - t556 * t420) * MDP(13) + (t556 * t417 + t426 * t420) * MDP(14) + (-t321 * t419 * t420 - t381 * t440) * MDP(15) + ((t532 + t535) * t498 + (t542 - t541 + (-t531 + t536) * qJD(5)) * t420) * MDP(16) + ((-t404 * t493 - t321) * t417 + (qJD(4) * t381 + t447) * t420) * MDP(17) + ((t404 * t499 - t322) * t417 + (-t448 - t500) * t420) * MDP(18) + (t376 * t417 + t404 * t497) * MDP(19) + (t376 * t530 + (-t385 * t495 - t563) * t404 + (-t325 * t499 + (-t448 + t500) * t422 + t466) * t417 + (t325 * t494 + t296 * t416 - t422 * t322 + (-t404 * t520 + t299) * qJD(4)) * t420 + t437) * MDP(20) + (-t511 * t376 + t562 * t404 + ((t404 * t422 + t326) * t495 + (-t325 * t419 + t381 * t422) * qJD(4) + t481) * t417 + (-qJD(4) * t300 + t296 * t419 + t321 * t422 - t325 * t495 + t465) * t420 + t444) * MDP(21) + (t314 * t379 + t322 * t351 - t341 * t376 + (-t301 * t499 - t288) * t417 - t514 * t404 + (-qJD(4) * t294 + t289 * t416 + t301 * t494) * t420 + t437) * MDP(22) + (-t321 * t341 - t322 * t340 + t514 * t381 - t515 * t379 - t457 * t498 + (-qJD(5) * t456 - t287 * t416 + t288 * t419 + t438) * t420) * MDP(23) + (-t314 * t381 + t321 * t351 + t340 * t376 + (t301 * t493 + t287) * t417 + t515 * t404 + (qJD(4) * t297 - t289 * t419 + t301 * t495 - t465) * t420 - t444) * MDP(24) + (t287 * t340 + t289 * t351 + t301 * t314 + t288 * t341 - g(1) * (pkin(5) * t319 - pkin(8) * t368 + qJ(6) * t318 - t360) - g(2) * (pkin(5) * t317 - pkin(8) * t366 + qJ(6) * t316 - t359) - g(3) * (pkin(5) * t350 + qJ(6) * t349 + t510) + t515 * t297 + t514 * t294 + (-g(3) * pkin(8) * t421 + (-g(3) * t463 + t301 * t506) * t418) * t413 + t559 * t385) * MDP(25); qJDD(2) * MDP(5) - t424 * MDP(6) + (t390 + t431 - t545) * MDP(7) + t408 * MDP(13) + t454 * MDP(25) + t484 * t379 * t498 + (-t384 * MDP(7) + t425) * qJD(2) + (t508 * MDP(14) - t289 * MDP(25) - t484 * t322 + t483 * t321 + ((t532 - t535) * MDP(23) + t456 * MDP(25) + t435 * t404) * qJD(4)) * t420 + (t508 * MDP(13) - qJDD(4) * MDP(14) + (-t541 - t542) * MDP(23) + t458 * MDP(25) + (t301 * MDP(25) + t381 * t483) * qJD(4) + t435 * t376 + t425 * qJD(5)) * t417; MDP(10) * t486 - MDP(11) * t487 + qJDD(4) * MDP(12) + (qJD(4) * t339 - t420 * t467 + t443 + t480) * MDP(13) + (t467 * t417 + t442 - t469) * MDP(14) + (t404 * t531 - t542) * MDP(15) + ((-t321 - t537) * t419 + (-t322 - t533) * t416) * MDP(16) + ((-t381 * t420 + t404 * t519) * qJD(2) + t448) * MDP(17) + ((-t404 * t522 + t534) * qJD(2) + t447) * MDP(18) - t404 * MDP(19) * t501 + (-t299 * t501 - pkin(4) * t322 - t339 * t379 + (t404 * t561 + t436) * t416 + (-t296 + (-t383 - t547) * t404 + t443) * t419) * MDP(20) + (pkin(4) * t321 + t513 * t404 + t300 * t501 - t339 * t381 + t436 * t419 + (t296 - t432) * t416) * MDP(21) + (t294 * t501 + t312 * t404 - t322 * t455 + t512 * t379 + t558 * t416 + t429 * t419) * MDP(22) + (t311 * t379 - t312 * t381 + (t287 + t404 * t294 + (-t322 + t496) * pkin(9)) * t419 + (t288 - t544 + (qJD(5) * t379 - t321) * pkin(9)) * t416 - t442) * MDP(23) + (-t297 * t501 - t311 * t404 - t321 * t455 - t512 * t381 + t429 * t416 - t558 * t419) * MDP(24) + (-t294 * t312 - t297 * t311 + t512 * t301 + (qJD(5) * t457 - t442 + t458) * pkin(9) + (-t289 + t443) * t455) * MDP(25) + (t420 * t417 * MDP(8) - t509 * MDP(9)) * t424; MDP(15) * t538 + (-t379 ^ 2 + t555) * MDP(16) + t303 * MDP(17) + (t533 - t322) * MDP(18) + t376 * MDP(19) + (-t325 * t381 + t430 + t543) * MDP(20) + (t299 * t404 + t325 * t379 - t428) * MDP(21) + (-t327 * t379 - t427 + t543 + 0.2e1 * t554) * MDP(22) + (pkin(5) * t321 - qJ(6) * t322 + (t297 - t300) * t381 + (t294 - t492) * t379) * MDP(23) + (0.2e1 * t546 - t301 * t379 + t327 * t381 + (0.2e1 * qJD(6) - t299) * t404 + t428) * MDP(24) + (t287 * qJ(6) - t288 * pkin(5) - t301 * t327 - t294 * t300 - g(1) * (-pkin(5) * t307 + qJ(6) * t308) - g(2) * (-pkin(5) * t309 + qJ(6) * t310) - g(3) * (pkin(5) * t451 + qJ(6) * t335) + t492 * t297) * MDP(25); (-t376 + t538) * MDP(22) + t303 * MDP(23) + (-t404 ^ 2 - t555) * MDP(24) + (t427 - t544 - t554) * MDP(25);];
tau  = t1;
