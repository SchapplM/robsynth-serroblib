% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:34
% EndTime: 2019-03-09 02:21:42
% DurationCPUTime: 6.26s
% Computational Cost: add. (4503->464), mult. (9847->592), div. (0->0), fcn. (7655->18), ass. (0->216)
t475 = cos(pkin(11));
t483 = cos(qJ(4));
t540 = qJD(1) * t483;
t448 = t475 * t540;
t473 = sin(pkin(11));
t479 = sin(qJ(4));
t541 = qJD(1) * t479;
t524 = t473 * t541;
t423 = t448 - t524;
t587 = qJD(5) + qJD(6);
t600 = t423 - t587;
t474 = sin(pkin(10));
t449 = pkin(1) * t474 + qJ(3);
t429 = qJD(1) * qJD(3) + qJDD(1) * t449;
t456 = t475 * qJDD(2);
t581 = pkin(7) * qJDD(1);
t391 = t456 + (-t429 - t581) * t473;
t409 = t473 * qJDD(2) + t475 * t429;
t392 = t475 * t581 + t409;
t441 = t449 * qJD(1);
t458 = t475 * qJD(2);
t403 = t458 + (-pkin(7) * qJD(1) - t441) * t473;
t417 = t473 * qJD(2) + t475 * t441;
t542 = qJD(1) * t475;
t404 = pkin(7) * t542 + t417;
t344 = t479 * t403 + t483 * t404;
t598 = qJD(4) * t344;
t498 = t391 * t483 - t479 * t392 - t598;
t317 = -qJDD(4) * pkin(4) - t498;
t420 = qJD(5) - t423;
t470 = pkin(11) + qJ(4);
t459 = sin(t470);
t461 = cos(t470);
t471 = qJ(1) + pkin(10);
t460 = sin(t471);
t462 = cos(t471);
t514 = g(1) * t462 + g(2) * t460;
t489 = -g(3) * t461 + t459 * t514;
t599 = -qJD(5) * pkin(8) * t420 - t317 + t489;
t432 = t473 * t483 + t475 * t479;
t424 = t432 * qJD(1);
t478 = sin(qJ(5));
t482 = cos(qJ(5));
t535 = t482 * qJD(4);
t405 = t424 * t478 - t535;
t481 = cos(qJ(6));
t407 = qJD(4) * t478 + t424 * t482;
t477 = sin(qJ(6));
t568 = t407 * t477;
t346 = t481 * t405 + t568;
t418 = qJD(6) + t420;
t597 = t346 * t418;
t501 = t405 * t477 - t481 * t407;
t596 = t418 * t501;
t552 = t477 * t482;
t434 = t478 * t481 + t552;
t545 = t600 * t434;
t426 = t432 * qJD(4);
t531 = qJDD(1) * t483;
t447 = t475 * t531;
t532 = qJDD(1) * t479;
t595 = -qJD(1) * t426 - t473 * t532 + t447;
t539 = qJD(5) * t478;
t566 = t423 * t478;
t594 = t539 - t566;
t341 = qJD(4) * pkin(8) + t344;
t476 = cos(pkin(10));
t451 = -pkin(1) * t476 - pkin(2);
t440 = -pkin(3) * t475 + t451;
t422 = qJD(1) * t440 + qJD(3);
t355 = -pkin(4) * t423 - pkin(8) * t424 + t422;
t322 = t341 * t482 + t355 * t478;
t312 = -pkin(9) * t405 + t322;
t537 = qJD(6) * t477;
t310 = t312 * t537;
t589 = t403 * t483 - t479 * t404;
t340 = -qJD(4) * pkin(4) - t589;
t329 = pkin(5) * t405 + t340;
t472 = qJ(5) + qJ(6);
t466 = sin(t472);
t557 = t462 * t466;
t467 = cos(t472);
t560 = t460 * t467;
t400 = -t461 * t560 + t557;
t556 = t462 * t467;
t561 = t460 * t466;
t402 = t461 * t556 + t561;
t584 = g(3) * t459;
t593 = g(1) * t402 - g(2) * t400 + t329 * t346 + t467 * t584 + t310;
t399 = t461 * t561 + t556;
t401 = -t461 * t557 + t560;
t503 = t391 * t479 + t392 * t483;
t316 = qJDD(4) * pkin(8) + qJD(4) * t589 + t503;
t526 = qJD(4) * t448 + t473 * t531 + t475 * t532;
t382 = -qJD(4) * t524 + t526;
t419 = qJDD(1) * t440 + qJDD(3);
t328 = -pkin(4) * t595 - pkin(8) * t382 + t419;
t327 = t482 * t328;
t337 = qJD(5) * t535 + t478 * qJDD(4) + t482 * t382 - t424 * t539;
t379 = qJDD(5) - t595;
t299 = pkin(5) * t379 - pkin(9) * t337 - qJD(5) * t322 - t316 * t478 + t327;
t338 = t407 * qJD(5) - t482 * qJDD(4) + t382 * t478;
t538 = qJD(5) * t482;
t492 = t482 * t316 + t478 * t328 - t341 * t539 + t355 * t538;
t300 = -pkin(9) * t338 + t492;
t521 = t481 * t299 - t477 * t300;
t592 = -g(1) * t401 + g(2) * t399 + t329 * t501 + t466 * t584 + t521;
t377 = qJDD(6) + t379;
t591 = t377 * MDP(27) + (-t346 ^ 2 + t501 ^ 2) * MDP(24) - t346 * MDP(23) * t501;
t372 = t434 * t432;
t433 = t477 * t478 - t481 * t482;
t546 = t600 * t433;
t586 = -t377 * t434 - t418 * t546;
t520 = t337 * t477 + t481 * t338;
t307 = -qJD(6) * t501 + t520;
t585 = pkin(8) + pkin(9);
t582 = pkin(7) + t449;
t321 = -t341 * t478 + t482 * t355;
t311 = -pkin(9) * t407 + t321;
t309 = pkin(5) * t420 + t311;
t580 = t309 * t481;
t579 = t312 * t481;
t578 = t337 * t478;
t577 = t346 * t424;
t576 = t501 * t424;
t574 = t379 * t478;
t572 = t405 * t420;
t571 = t405 * t424;
t570 = t407 * t420;
t569 = t407 * t424;
t567 = (-t441 * t473 + t458) * t473;
t431 = t473 * t479 - t483 * t475;
t425 = t431 * qJD(4);
t565 = t425 * t478;
t564 = t425 * t482;
t563 = t432 * t478;
t562 = t432 * t482;
t559 = t460 * t478;
t558 = t460 * t482;
t555 = t462 * t478;
t554 = t462 * t482;
t553 = t475 * MDP(5);
t427 = t582 * t473;
t428 = t582 * t475;
t375 = -t427 * t479 + t428 * t483;
t364 = t482 * t375;
t367 = t482 * t379;
t536 = qJD(6) * t481;
t527 = t481 * t337 - t477 * t338 - t405 * t536;
t306 = -t407 * t537 + t527;
t551 = t306 * t431 - t426 * t501;
t523 = t432 * t539;
t319 = -t425 * t552 - t477 * t523 - t537 * t563 + (t562 * t587 - t565) * t481;
t550 = -t319 * t418 - t372 * t377;
t549 = t337 * t431 + t407 * t426;
t380 = pkin(4) * t424 - pkin(8) * t423;
t548 = t478 * t380 + t482 * t589;
t365 = pkin(4) * t431 - pkin(8) * t432 + t440;
t547 = t478 * t365 + t364;
t544 = t473 ^ 2 + t475 ^ 2;
t533 = qJDD(1) * t451;
t529 = t379 * t563;
t528 = t432 * t367;
t525 = qJD(5) * t585;
t522 = t432 * t538;
t518 = t420 * t482;
t517 = -qJD(5) * t355 - t316;
t516 = qJD(6) * t309 + t300;
t515 = t594 * pkin(5) - t344;
t513 = g(1) * t460 - g(2) * t462;
t480 = sin(qJ(1));
t484 = cos(qJ(1));
t512 = g(1) * t480 - g(2) * t484;
t511 = -t433 * t377 + t545 * t418;
t510 = -t341 * t538 + t327;
t370 = t482 * t380;
t443 = t585 * t482;
t509 = pkin(5) * t424 + qJD(6) * t443 - t589 * t478 + t370 + (-pkin(9) * t423 + t525) * t482;
t442 = t585 * t478;
t508 = -pkin(9) * t566 + qJD(6) * t442 + t478 * t525 + t548;
t506 = -t307 * t431 - t346 * t426;
t302 = t309 * t477 + t579;
t318 = -t372 * t587 + t433 * t425;
t373 = t433 * t432;
t505 = -t318 * t418 + t373 * t377;
t504 = -t338 * t431 - t405 * t426;
t408 = -t429 * t473 + t456;
t500 = -t408 * t473 + t409 * t475;
t499 = -t427 * t483 - t428 * t479;
t497 = -t594 * t420 + t367;
t496 = t522 - t565;
t495 = t523 + t564;
t494 = -pkin(8) * t379 + t340 * t420;
t351 = -qJD(3) * t431 + qJD(4) * t499;
t381 = pkin(4) * t426 + pkin(8) * t425;
t491 = t482 * t351 + t365 * t538 - t375 * t539 + t478 * t381;
t352 = qJD(3) * t432 + qJD(4) * t375;
t454 = -pkin(5) * t482 - pkin(4);
t437 = qJDD(3) + t533;
t413 = t461 * t554 + t559;
t412 = -t461 * t555 + t558;
t411 = -t461 * t558 + t555;
t410 = t461 * t559 + t554;
t385 = -qJD(4) * t425 + qJDD(4) * t432;
t384 = -qJD(4) * t426 - qJDD(4) * t431;
t371 = t482 * t381;
t362 = t482 * t365;
t353 = pkin(5) * t563 - t499;
t325 = pkin(5) * t496 + t352;
t324 = -pkin(9) * t563 + t547;
t323 = pkin(5) * t431 - pkin(9) * t562 - t375 * t478 + t362;
t308 = pkin(5) * t338 + t317;
t304 = -pkin(9) * t496 + t491;
t303 = pkin(9) * t564 + pkin(5) * t426 - t351 * t478 + t371 + (-t364 + (pkin(9) * t432 - t365) * t478) * qJD(5);
t301 = -t312 * t477 + t580;
t1 = [(-qJD(4) * t351 - qJDD(4) * t375 + t382 * t440 + t419 * t432 - t422 * t425 - t459 * t513) * MDP(15) + (t382 * t432 - t424 * t425) * MDP(9) + (-(-t405 * t482 - t407 * t478) * t425 + (-t578 - t338 * t482 + (t405 * t478 - t407 * t482) * qJD(5)) * t432) * MDP(17) + (-t382 * t431 - t423 * t425 - t424 * t426 + t432 * t595) * MDP(10) + (-qJD(4) * t352 + qJDD(4) * t499 + t419 * t431 + t422 * t426 - t440 * t595 + t461 * t513) * MDP(14) + qJDD(1) * MDP(1) + (-MDP(6) * t473 + t553) * (-t437 + t513 - t533) + (-g(1) * t399 - g(2) * t401 - t302 * t426 + t353 * t306 - t308 * t373 + t310 * t431 + t329 * t318 - t325 * t501 + (-(-qJD(6) * t324 + t303) * t418 - t323 * t377 - t299 * t431) * t477 + (-(qJD(6) * t323 + t304) * t418 - t324 * t377 - t516 * t431) * t481) * MDP(29) + (-t306 * t373 - t318 * t501) * MDP(23) + (-t306 * t372 + t307 * t373 - t318 * t346 + t319 * t501) * MDP(24) + (-g(1) * t410 - g(2) * t412 + t317 * t562 - t322 * t426 - t337 * t499 - t340 * t495 + t352 * t407 - t379 * t547 - t420 * t491 - t431 * t492) * MDP(22) + ((-t375 * t538 + t371) * t420 + t362 * t379 + t510 * t431 + t321 * t426 + t352 * t405 - t499 * t338 + t340 * t522 - g(1) * t411 - g(2) * t413 + ((-qJD(5) * t365 - t351) * t420 - t375 * t379 + t517 * t431 + t317 * t432 - t340 * t425) * t478) * MDP(21) + (t512 + (t474 ^ 2 + t476 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t429 * t544 + t500 - t514) * MDP(7) + t512 * MDP(2) + ((t303 * t481 - t304 * t477) * t418 + (t323 * t481 - t324 * t477) * t377 + t521 * t431 + t301 * t426 + t325 * t346 + t353 * t307 + t308 * t372 + t329 * t319 - g(1) * t400 - g(2) * t402 + ((-t323 * t477 - t324 * t481) * t418 - t302 * t431) * qJD(6)) * MDP(28) + t384 * MDP(12) + t385 * MDP(11) + (-t420 * t496 + t504 - t529) * MDP(19) + (t377 * t431 + t418 * t426) * MDP(27) + (t379 * t431 + t420 * t426) * MDP(20) + (-t420 * t495 + t528 + t549) * MDP(18) + (t506 + t550) * MDP(26) + (-t505 + t551) * MDP(25) + (t337 * t562 - t407 * t495) * MDP(16) + (t437 * t451 - g(1) * (-pkin(1) * t480 - pkin(2) * t460 + qJ(3) * t462) - g(2) * (pkin(1) * t484 + pkin(2) * t462 + qJ(3) * t460) + t500 * t449 + (t417 * t475 - t567) * qJD(3)) * MDP(8) + (g(1) * t484 + g(2) * t480) * MDP(3); (qJDD(2) - g(3)) * MDP(4) + (t408 * t475 + t409 * t473 - g(3)) * MDP(8) + t384 * MDP(14) - t385 * MDP(15) + (-t504 - t529) * MDP(21) + (-t528 + t549) * MDP(22) + (-t506 + t550) * MDP(28) + (t505 + t551) * MDP(29) + (-MDP(21) * t496 + t495 * MDP(22)) * t420; (qJD(1) * t567 - t417 * t542 + qJDD(3) - t513) * MDP(8) - t447 * MDP(14) + t526 * MDP(15) + (t497 - t571) * MDP(21) + (-t420 ^ 2 * t482 - t569 - t574) * MDP(22) + (t511 - t577) * MDP(28) + (t576 + t586) * MDP(29) - t544 * MDP(7) * qJD(1) ^ 2 + (-t553 + t451 * MDP(8) + (MDP(14) * t479 + MDP(6)) * t473) * qJDD(1) + ((t473 * t540 + t475 * t541 + t424) * MDP(14) + (t423 - t524) * MDP(15)) * qJD(4); -t423 ^ 2 * MDP(10) + ((-t423 - t524) * qJD(4) + t526) * MDP(11) + t595 * MDP(12) + qJDD(4) * MDP(13) + (t489 + t498 + t598) * MDP(14) + (-t422 * t423 + t514 * t461 - t503 + t584) * MDP(15) + (t407 * t518 + t578) * MDP(16) + ((t337 - t572) * t482 + (-t338 - t570) * t478) * MDP(17) + (t420 * t518 - t569 + t574) * MDP(18) + (t497 + t571) * MDP(19) + (-pkin(4) * t338 - t344 * t405 - t370 * t420 + (t420 * t589 + t494) * t478 + t599 * t482) * MDP(21) + (-pkin(4) * t337 - t344 * t407 + t548 * t420 - t478 * t599 + t494 * t482) * MDP(22) + (t306 * t434 - t501 * t546) * MDP(23) + (-t306 * t433 - t307 * t434 - t346 * t546 - t501 * t545) * MDP(24) + (t576 - t586) * MDP(25) + (t511 + t577) * MDP(26) + ((-t442 * t481 - t443 * t477) * t377 + t454 * t307 + t308 * t433 + (t477 * t508 - t481 * t509) * t418 + t515 * t346 - t545 * t329 + t489 * t467) * MDP(28) + (-(-t442 * t477 + t443 * t481) * t377 + t454 * t306 + t308 * t434 + (t477 * t509 + t481 * t508) * t418 - t515 * t501 + t546 * t329 - t489 * t466) * MDP(29) + (t424 * MDP(10) + MDP(12) * qJD(4) - t422 * MDP(14) - t420 * MDP(20) - t321 * MDP(21) + t322 * MDP(22) - t418 * MDP(27) - t301 * MDP(28) + t302 * MDP(29) - MDP(9) * t423) * t424; t407 * t405 * MDP(16) + (-t405 ^ 2 + t407 ^ 2) * MDP(17) + (t337 + t572) * MDP(18) + (-t338 + t570) * MDP(19) + t379 * MDP(20) + (-g(1) * t412 + g(2) * t410 + t322 * t420 - t340 * t407 + (t517 + t584) * t478 + t510) * MDP(21) + (g(1) * t413 - g(2) * t411 + t321 * t420 + t340 * t405 + t482 * t584 - t492) * MDP(22) + (t306 + t597) * MDP(25) + (-t307 - t596) * MDP(26) + (-(-t311 * t477 - t579) * t418 - t302 * qJD(6) + (-t346 * t407 + t377 * t481 - t418 * t537) * pkin(5) + t592) * MDP(28) + ((-t312 * t418 - t299) * t477 + (t311 * t418 - t516) * t481 + (-t377 * t477 + t407 * t501 - t418 * t536) * pkin(5) + t593) * MDP(29) + t591; (t527 + t597) * MDP(25) + (-t520 - t596) * MDP(26) + (t302 * t418 + t592) * MDP(28) + (-t477 * t299 - t481 * t300 + t301 * t418 + t593) * MDP(29) + (-MDP(25) * t568 + MDP(26) * t501 - MDP(28) * t302 - MDP(29) * t580) * qJD(6) + t591;];
tau  = t1;
