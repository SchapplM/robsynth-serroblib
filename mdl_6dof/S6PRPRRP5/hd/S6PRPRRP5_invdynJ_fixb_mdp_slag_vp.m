% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP5
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
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:52:05
% EndTime: 2021-01-16 01:52:18
% DurationCPUTime: 6.73s
% Computational Cost: add. (3037->511), mult. (6240->670), div. (0->0), fcn. (4531->10), ass. (0->233)
t429 = sin(qJ(5));
t432 = cos(qJ(5));
t430 = sin(qJ(4));
t510 = t432 * qJD(4);
t491 = t430 * t510;
t433 = cos(qJ(4));
t512 = qJD(5) * t433;
t450 = t429 * t512 + t491;
t504 = qJDD(2) * t433;
t320 = qJD(2) * t450 - qJD(5) * t510 - t429 * qJDD(4) - t432 * t504;
t521 = qJD(2) * t433;
t389 = t429 * t521 - t510;
t522 = qJD(2) * t430;
t413 = qJD(5) + t522;
t559 = t389 * t413;
t591 = -t320 - t559;
t425 = sin(pkin(6));
t434 = cos(qJ(2));
t431 = sin(qJ(2));
t544 = t430 * t431;
t356 = (-t429 * t544 + t432 * t434) * t425;
t471 = pkin(4) * t433 + pkin(9) * t430;
t386 = qJD(4) * t471 + qJD(3);
t590 = -qJD(1) * t356 + t432 * t386;
t541 = t431 * t432;
t357 = (t429 * t434 + t430 * t541) * t425;
t395 = pkin(4) * t430 - pkin(9) * t433 + qJ(3);
t435 = pkin(2) + pkin(8);
t515 = qJD(4) * t435;
t490 = t433 * t515;
t513 = qJD(5) * t432;
t589 = qJD(1) * t357 - t429 * t386 - t395 * t513 + t432 * t490;
t509 = qJD(2) * qJD(4);
t588 = -t430 * t509 + t504;
t481 = t429 * t435 + pkin(5);
t511 = qJD(6) * t432;
t514 = qJD(5) * t429;
t543 = t430 * t432;
t405 = t435 * t543;
t585 = -t429 * t395 + t405;
t536 = qJ(6) * t491 + t585 * qJD(5) + (qJ(6) * t514 + qJD(4) * t481 - t511) * t433 + t590;
t487 = t432 * t512;
t587 = -qJ(6) * t487 + (-qJD(6) * t433 + (qJ(6) * qJD(4) + qJD(5) * t435) * t430) * t429 - t589;
t518 = qJD(4) * t429;
t391 = t432 * t521 + t518;
t583 = qJD(5) * t391;
t321 = -t432 * qJDD(4) + t429 * t588 + t583;
t556 = t391 * t413;
t586 = t321 + t556;
t527 = qJD(1) * t425;
t495 = t434 * t527;
t468 = qJD(3) - t495;
t378 = -qJD(2) * t435 + t468;
t427 = cos(pkin(6));
t526 = qJD(1) * t427;
t584 = t378 * t433 - t430 * t526;
t424 = sin(pkin(10));
t426 = cos(pkin(10));
t547 = t427 * t431;
t373 = t424 * t434 + t426 * t547;
t462 = t424 * t547 - t426 * t434;
t469 = -g(1) * t462 + g(2) * t373;
t537 = qJDD(1) - g(3);
t551 = t425 * t431;
t582 = -t537 * t551 + t469;
t496 = t431 * t527;
t525 = qJD(2) * qJ(3);
t394 = t496 + t525;
t581 = qJD(4) * (-t394 + t496 - t525) + qJDD(4) * t435;
t546 = t427 * t434;
t374 = t424 * t546 + t426 * t431;
t550 = t425 * t433;
t334 = t374 * t430 + t424 * t550;
t372 = t424 * t431 - t426 * t546;
t336 = -t372 * t430 + t426 * t550;
t363 = t373 * t432;
t563 = t462 * t432;
t542 = t430 * t434;
t377 = -t425 * t542 + t427 * t433;
t339 = -t377 * t429 + t425 * t541;
t576 = g(3) * t339;
t580 = -g(1) * (-t334 * t429 - t563) - t576 - g(2) * (t336 * t429 + t363);
t579 = t391 ^ 2;
t578 = pkin(5) * t389;
t577 = g(1) * t374;
t340 = t377 * t432 + t429 * t551;
t575 = g(3) * t340;
t538 = t433 * t434;
t376 = t425 * t538 + t427 * t430;
t574 = g(3) * t376;
t485 = t433 * t509;
t505 = qJDD(2) * t430;
t453 = t485 + t505;
t385 = qJDD(5) + t453;
t573 = t385 * pkin(5);
t428 = -qJ(6) - pkin(9);
t572 = qJDD(2) * pkin(2);
t409 = t433 * t526;
t343 = t430 * t378 + t409;
t331 = qJD(4) * pkin(9) + t343;
t352 = qJD(2) * t395 + t496;
t314 = t331 * t432 + t352 * t429;
t306 = -qJ(6) * t389 + t314;
t571 = t306 * t413;
t570 = t320 * t429;
t569 = t321 * t432;
t361 = t372 * t429;
t567 = t373 * t429;
t566 = t374 * t429;
t565 = t374 * t432;
t564 = t462 * t429;
t561 = t385 * t429;
t560 = t385 * t432;
t558 = t389 * t429;
t557 = t389 * t432;
t555 = t391 * t429;
t554 = t391 * t432;
t553 = t424 * t427;
t552 = t425 * t430;
t549 = t425 * t434;
t548 = t426 * t427;
t545 = t427 * t435;
t540 = t431 * t433;
t539 = t432 * t433;
t313 = -t331 * t429 + t432 * t352;
t305 = -qJ(6) * t391 + t313;
t302 = pkin(5) * t413 + t305;
t534 = -t305 + t302;
t393 = t471 * qJD(2);
t533 = t429 * t393 + t432 * t584;
t477 = qJD(5) * t428;
t492 = t429 * t522;
t532 = -qJ(6) * t492 + t429 * t477 + t511 - t533;
t380 = t432 * t393;
t531 = t432 * t477 - t380 - (pkin(5) * t433 + qJ(6) * t543) * qJD(2) + (-qJD(6) + t584) * t429;
t423 = t433 ^ 2;
t529 = t430 ^ 2 - t423;
t436 = qJD(4) ^ 2;
t437 = qJD(2) ^ 2;
t528 = -t436 - t437;
t524 = qJD(2) * t394;
t523 = qJD(2) * t425;
t520 = qJD(4) * t389;
t519 = qJD(4) * t391;
t517 = qJD(4) * t430;
t516 = qJD(4) * t433;
t508 = qJDD(1) * t425;
t507 = qJDD(1) * t427;
t506 = qJDD(2) * qJ(3);
t502 = MDP(20) + MDP(22);
t501 = MDP(21) + MDP(23);
t500 = qJ(3) * t548;
t499 = -qJD(4) * t409 - t378 * t517 - t430 * t507;
t494 = t431 * t523;
t493 = t434 * t523;
t330 = -qJD(4) * pkin(4) - t584;
t476 = -qJD(6) - t578;
t319 = t330 - t476;
t489 = t319 * t513;
t488 = t413 * t514;
t484 = t431 * t508;
t483 = t434 * t508;
t482 = t433 * t507;
t480 = pkin(5) * t429 + t435;
t472 = t433 * t496;
t479 = -g(3) * t357 + t389 * t472;
t478 = -g(3) * t356 + t391 * t472;
t401 = qJD(1) * t494;
t459 = qJDD(3) + t401 - t483;
t345 = -qJDD(2) * t435 + t459;
t475 = -t345 + t524;
t474 = qJD(5) * t430 + qJD(2);
t310 = qJDD(4) * pkin(9) + qJD(4) * t584 + t345 * t430 + t482;
t318 = t484 + t395 * qJDD(2) + (t386 + t495) * qJD(2);
t473 = -t432 * t310 - t429 * t318 + t331 * t514 - t352 * t513;
t333 = -t374 * t433 + t424 * t552;
t335 = t372 * t433 + t426 * t552;
t470 = g(1) * t333 - g(2) * t335;
t467 = t302 * t432 + t306 * t429;
t466 = t302 * t429 - t306 * t432;
t415 = pkin(5) * t432 + pkin(4);
t465 = t415 * t430 + t428 * t433;
t464 = -g(2) * t372 + g(3) * t549 - t577;
t461 = t427 * t538 - t552;
t460 = t427 * t542 + t550;
t457 = t413 * t513 + t561;
t456 = -t488 + t560;
t454 = t321 * qJ(6) + t473;
t452 = t470 + t574;
t451 = g(1) * t334 - g(2) * t336 + g(3) * t377;
t448 = g(3) * t551 + t469;
t311 = -qJDD(4) * pkin(4) - t345 * t433 - t499;
t301 = pkin(5) * t321 + qJDD(6) + t311;
t447 = -t301 + t452;
t446 = -pkin(9) * t385 + t330 * t413;
t445 = pkin(9) * qJD(5) * t413 + g(1) * (t424 * t461 + t426 * t540) - g(2) * (-t424 * t540 + t426 * t461) + t311;
t444 = -t464 + t483;
t443 = qJDD(3) - t444;
t317 = t432 * t318;
t442 = -qJD(5) * t314 - t429 * t310 + t317;
t440 = t320 * qJ(6) + t442;
t439 = (t554 + t558) * MDP(24) - t467 * MDP(25);
t346 = t484 + t506 + (qJD(3) + t495) * qJD(2);
t438 = qJD(2) * t468 + t435 * t436 + t346 - t448 + t506;
t420 = qJDD(4) * t433;
t418 = t426 * qJ(3);
t416 = t424 * qJ(3);
t410 = qJ(3) * t551;
t407 = qJ(3) * t553;
t398 = t428 * t432;
t397 = t428 * t429;
t388 = -qJD(2) * pkin(2) + t468;
t387 = t480 * t433;
t384 = t389 ^ 2;
t383 = t432 * t395;
t362 = t372 * t432;
t359 = t429 * t574;
t355 = t373 * t430;
t354 = t462 * t430;
t353 = t459 - t572;
t347 = -t430 * t515 + (-t429 * t517 + t487) * pkin(5);
t338 = qJD(4) * t377 - t433 * t494;
t337 = -qJD(4) * t376 + t430 * t494;
t332 = -qJ(6) * t429 * t433 - t585;
t327 = -t424 * t544 + t426 * t460;
t325 = t424 * t460 + t426 * t544;
t324 = -qJ(6) * t539 + t430 * t481 + t383;
t322 = -pkin(5) * t492 + t343;
t308 = qJD(5) * t339 + t337 * t432 + t429 * t493;
t307 = -qJD(5) * t340 - t337 * t429 + t432 * t493;
t296 = -t389 * qJD(6) - t454;
t295 = -t391 * qJD(6) + t440 + t573;
t1 = [t537 * MDP(1) + (qJDD(1) * t427 ^ 2 - g(3)) * MDP(7) + (-qJD(4) * t338 - qJDD(4) * t376) * MDP(13) + (-qJD(4) * t337 - qJDD(4) * t377) * MDP(14) + (-t307 * t391 - t308 * t389 + t320 * t339 - t321 * t340) * MDP(24) + (t295 * t339 + t296 * t340 + t301 * t376 + t302 * t307 + t306 * t308 + t319 * t338 - g(3)) * MDP(25) + t501 * (-t308 * t413 - t320 * t376 + t338 * t391 - t340 * t385) + t502 * (t307 * t413 + t321 * t376 + t338 * t389 + t339 * t385) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t431 + t434 * t437) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t434 + t431 * t437) + ((-t353 + t524) * MDP(7) + (MDP(13) * t430 + MDP(14) * t433) * t437) * t434 + ((qJD(2) * t388 + t346) * MDP(7) + t453 * MDP(13) + t588 * MDP(14)) * t431) * t425; qJDD(2) * MDP(2) + t444 * MDP(3) + t582 * MDP(4) + (t443 - 0.2e1 * t572) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t506 - t582) * MDP(6) + (t346 * qJ(3) + t394 * qJD(3) - t353 * pkin(2) - g(1) * (-(pkin(2) * t426 + t407) * t431 + (-pkin(2) * t553 + t418) * t434) - g(2) * (-(pkin(2) * t424 - t500) * t431 + (pkin(2) * t548 + t416) * t434) - g(3) * (pkin(2) * t549 + t410) + (-t388 * t431 - t394 * t434) * t527) * MDP(7) + (qJDD(2) * t423 - 0.2e1 * t430 * t485) * MDP(8) + 0.2e1 * (-t430 * t504 + t529 * t509) * MDP(9) + (-t430 * t436 + t420) * MDP(10) + (-qJDD(4) * t430 - t433 * t436) * MDP(11) + (t438 * t430 - t581 * t433) * MDP(13) + (t581 * t430 + t438 * t433) * MDP(14) + (-t320 * t539 - t391 * t450) * MDP(15) + ((t555 + t557) * t517 + (t570 - t569 + (-t554 + t558) * qJD(5)) * t433) * MDP(16) + ((-t413 * t510 - t320) * t430 + (t456 + t519) * t433) * MDP(17) + ((t413 * t518 - t321) * t430 + (-t457 - t520) * t433) * MDP(18) + (t385 * t430 + t413 * t516) * MDP(19) + (t383 * t385 + g(1) * t354 * t432 - g(2) * (t355 * t432 - t361) + (-t331 * t513 - t389 * t515 + t317) * t430 + (qJD(5) * t405 + t590) * t413 + (qJD(4) * t313 + t321 * t435 + t330 * t513) * t433 + ((-qJD(5) * t395 + t490) * t413 + t311 * t433 + t577 + (-qJD(4) * t330 - qJD(5) * t352 + t385 * t435 - t310) * t430) * t429 + t479) * MDP(20) + (t585 * t385 - g(1) * (t354 * t429 - t565) - g(2) * (-t355 * t429 - t362) + t589 * t413 + (-t435 * t488 + (-t330 * t432 - t391 * t435) * qJD(4) + t473) * t430 + (-qJD(4) * t314 + t311 * t432 - t320 * t435 - t330 * t514) * t433 + t478) * MDP(21) + (g(1) * t566 + g(2) * t361 + t387 * t321 + t324 * t385 + t347 * t389 + (-t319 * t518 - t432 * t469 + t295) * t430 + t536 * t413 + (qJD(4) * t302 + t301 * t429 + t489) * t433 + t479) * MDP(22) + (g(1) * t565 + g(2) * t362 - t387 * t320 - t332 * t385 + t347 * t391 + (-t319 * t510 + t429 * t469 - t296) * t430 - t587 * t413 + (-qJD(4) * t306 + t301 * t432 - t319 * t514) * t433 + t478) * MDP(23) + (t320 * t324 - t321 * t332 - t536 * t391 - t587 * t389 + t467 * t517 + (qJD(5) * t466 - t295 * t432 - t296 * t429 + t448) * t433) * MDP(24) + (t296 * t332 + t295 * t324 + t301 * t387 - g(1) * (-pkin(5) * t566 - (t426 * t435 + t407) * t431 + (-t424 * t545 + t418) * t434 - t465 * t462) - g(2) * (-pkin(5) * t361 - (t424 * t435 - t500) * t431 + (t426 * t545 + t416) * t434 + t465 * t373) - g(3) * (t410 + (t431 * t465 + t434 * t480) * t425) + (t347 + t472) * t319 + t587 * t306 + t536 * t302) * MDP(25); qJDD(2) * MDP(5) - t437 * MDP(6) + (t401 + t443 - t572) * MDP(7) + t420 * MDP(13) + t464 * MDP(25) + t502 * (-t321 * t433 + (t520 - t561) * t430 + (-t429 * t516 - t432 * t474) * t413) + t501 * (t320 * t433 + (t519 - t560) * t430 + (t429 * t474 - t433 * t510) * t413) + (-t394 * MDP(7) + t439) * qJD(2) + (t528 * MDP(14) - t301 * MDP(25) + ((t555 - t557) * MDP(24) - t466 * MDP(25)) * qJD(4)) * t433 + (t528 * MDP(13) - qJDD(4) * MDP(14) + (-t569 - t570) * MDP(24) + (qJD(4) * t319 - t295 * t429 + t296 * t432) * MDP(25) + t439 * qJD(5)) * t430; MDP(10) * t504 - MDP(11) * t505 + qJDD(4) * MDP(12) + (qJD(4) * t343 - t433 * t475 + t452 + t499) * MDP(13) + (t475 * t430 + t451 - t482) * MDP(14) + (t413 * t554 - t570) * MDP(15) + (-t586 * t429 + t591 * t432) * MDP(16) + ((-t391 * t433 + t413 * t543) * qJD(2) + t457) * MDP(17) + ((-t413 * t429 * t430 + t389 * t433) * qJD(2) + t456) * MDP(18) - t413 * MDP(19) * t521 + (-t313 * t521 - pkin(4) * t321 - t343 * t389 - t380 * t413 + (t413 * t584 + t446) * t429 + (-t445 + t574) * t432) * MDP(20) + (pkin(4) * t320 + t314 * t521 - t343 * t391 + t413 * t533 + t429 * t445 + t432 * t446 - t359) * MDP(21) + (-t302 * t521 - t321 * t415 - t322 * t389 + t385 * t397 + t531 * t413 + (t319 * t522 + (t319 + t578) * qJD(5)) * t429 + t447 * t432) * MDP(22) + (t489 + t320 * t415 - t322 * t391 + t385 * t398 - t359 - t532 * t413 + (t306 * t433 + t319 * t543) * qJD(2) + (pkin(5) * t583 + t301 - t470) * t429) * MDP(23) + (t320 * t397 + t321 * t398 - t531 * t391 - t532 * t389 + (-t302 * t413 + t296) * t432 + (-t295 - t571) * t429 - t451) * MDP(24) + (-t296 * t398 + t295 * t397 - t301 * t415 - g(1) * (-t333 * t415 - t334 * t428) - g(2) * (t335 * t415 + t336 * t428) - g(3) * (-t376 * t415 - t377 * t428) + (pkin(5) * t514 - t322) * t319 + t532 * t306 + t531 * t302) * MDP(25) + (MDP(8) * t430 * t433 - MDP(9) * t529) * t437; t391 * t389 * MDP(15) + (-t384 + t579) * MDP(16) + (-t320 + t559) * MDP(17) + (-t321 + t556) * MDP(18) + t385 * MDP(19) + (t314 * t413 - t330 * t391 - g(1) * (-t325 * t429 - t563) - g(2) * (t327 * t429 + t363) - t576 + t442) * MDP(20) + (t313 * t413 + t330 * t389 - g(1) * (-t325 * t432 + t564) - g(2) * (t327 * t432 - t567) + t575 + t473) * MDP(21) + (0.2e1 * t573 + t571 + (-t319 + t476) * t391 + t440 + t580) * MDP(22) + (t305 * t413 - t579 * pkin(5) - g(1) * (-t334 * t432 + t564) - g(2) * (t336 * t432 - t567) + t575 + (qJD(6) + t319) * t389 + t454) * MDP(23) + (pkin(5) * t320 - t389 * t534) * MDP(24) + (t534 * t306 + (-t319 * t391 + t295 + t580) * pkin(5)) * MDP(25); t586 * MDP(22) + t591 * MDP(23) + (-t384 - t579) * MDP(24) + (t302 * t391 + t306 * t389 - t447) * MDP(25);];
tau = t1;
