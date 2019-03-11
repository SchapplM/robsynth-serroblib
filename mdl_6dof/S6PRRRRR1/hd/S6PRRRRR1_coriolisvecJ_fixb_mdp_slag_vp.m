% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:40:09
% EndTime: 2019-03-09 00:40:17
% DurationCPUTime: 5.82s
% Computational Cost: add. (6076->403), mult. (15379->564), div. (0->0), fcn. (12258->12), ass. (0->198)
t474 = cos(qJ(6));
t545 = qJD(6) * t474;
t471 = sin(qJ(4));
t472 = sin(qJ(3));
t552 = qJD(2) * t472;
t531 = t471 * t552;
t476 = cos(qJ(4));
t477 = cos(qJ(3));
t550 = qJD(2) * t477;
t532 = t476 * t550;
t432 = -t531 + t532;
t433 = -t471 * t550 - t476 * t552;
t470 = sin(qJ(5));
t475 = cos(qJ(5));
t394 = t475 * t432 + t433 * t470;
t605 = t394 * t474;
t610 = t545 - t605;
t439 = t471 * t477 + t472 * t476;
t464 = qJD(3) + qJD(4);
t412 = t464 * t439;
t400 = t412 * qJD(2);
t473 = sin(qJ(2));
t467 = sin(pkin(6));
t555 = qJD(1) * t467;
t536 = t473 * t555;
t589 = pkin(8) + pkin(9);
t526 = qJD(2) * t589 + t536;
t468 = cos(pkin(6));
t554 = qJD(1) * t468;
t416 = t472 * t554 + t477 * t526;
t478 = cos(qJ(2));
t553 = qJD(2) * t467;
t530 = qJD(1) * t553;
t510 = t478 * t530;
t381 = -qJD(3) * t416 - t472 * t510;
t549 = qJD(4) * t471;
t521 = t471 * t381 - t416 * t549;
t415 = -t526 * t472 + t477 * t554;
t380 = t415 * qJD(3) + t477 * t510;
t585 = qJD(3) * pkin(3);
t410 = t415 + t585;
t595 = t476 * (qJD(4) * t410 + t380);
t322 = -pkin(10) * t400 + t521 + t595;
t427 = t433 * pkin(10);
t407 = t471 * t416;
t519 = t476 * t410 - t407;
t363 = t427 + t519;
t356 = pkin(4) * t464 + t363;
t539 = qJD(2) * qJD(3);
t529 = t477 * t539;
t399 = qJD(4) * t532 - t464 * t531 + t476 * t529;
t409 = t476 * t416;
t497 = -t410 * t471 - t409;
t522 = -t471 * t380 + t476 * t381;
t485 = qJD(4) * t497 + t522;
t323 = -pkin(10) * t399 + t485;
t587 = pkin(10) * t432;
t364 = -t497 + t587;
t548 = qJD(5) * t470;
t523 = t323 * t470 - t364 * t548;
t307 = t475 * (qJD(5) * t356 + t322) + t523;
t460 = -pkin(3) * t477 - pkin(2);
t535 = t478 * t555;
t426 = qJD(2) * t460 - t535;
t402 = -pkin(4) * t432 + t426;
t576 = t394 * t402;
t609 = -t307 - t576;
t537 = qJD(3) * t589;
t440 = t472 * t537;
t441 = t477 * t537;
t445 = t589 * t472;
t446 = t589 * t477;
t495 = t445 * t471 - t446 * t476;
t608 = qJD(4) * t495 + t439 * t535 + t471 * t440 - t476 * t441;
t438 = t471 * t472 - t476 * t477;
t572 = t445 * t476;
t607 = qJD(4) * t572 - t438 * t535 + t476 * t440 + t471 * t441 + t446 * t549;
t496 = t432 * t470 - t475 * t433;
t469 = sin(qJ(6));
t546 = qJD(6) * t469;
t547 = qJD(5) * t475;
t345 = t475 * t399 - t470 * t400 + t432 * t547 + t433 * t548;
t463 = qJD(5) + t464;
t560 = t474 * t345 + t463 * t545;
t335 = -t496 * t546 + t560;
t331 = t335 * t469;
t332 = t335 * t474;
t377 = t463 * t469 + t474 * t496;
t582 = t345 * t469;
t336 = t377 * qJD(6) + t582;
t346 = qJD(5) * t496 + t399 * t470 + t475 * t400;
t344 = t474 * t346;
t574 = t496 * t469;
t375 = -t474 * t463 + t574;
t540 = -qJD(6) + t394;
t342 = t469 * t346;
t561 = -t540 * t545 + t342;
t604 = t540 * t469;
t606 = -t346 * MDP(22) - t394 ^ 2 * MDP(20) + (-t394 * t463 + t345) * MDP(21) + (-MDP(19) * t394 + MDP(20) * t496 + MDP(22) * t463 + MDP(30) * t540) * t496 + (t610 * t377 + t331) * MDP(26) + (-t377 * t496 + t540 * t605 + t561) * MDP(28) + (t375 * t496 - t540 * t604 + t344) * MDP(29) + (-t469 * t336 - t610 * t375 + t377 * t604 + t332) * MDP(27);
t580 = t364 * t470;
t328 = t356 * t475 - t580;
t326 = -pkin(5) * t463 - t328;
t584 = t326 * t394;
t603 = -pkin(10) * t412 - t607;
t411 = t464 * t438;
t602 = -pkin(10) * t411 - t608;
t579 = t364 * t475;
t329 = t356 * t470 + t579;
t524 = t322 * t470 - t475 * t323;
t308 = qJD(5) * t329 + t524;
t578 = t496 * t402;
t601 = -t308 - t578;
t489 = -t472 * t585 + t536;
t367 = pkin(5) * t496 - pkin(11) * t394;
t597 = MDP(5) * t472;
t596 = MDP(6) * (t472 ^ 2 - t477 ^ 2);
t505 = pkin(4) * t412 - t489;
t594 = t472 * MDP(10) + t477 * MDP(11);
t327 = pkin(11) * t463 + t329;
t347 = -pkin(5) * t394 - pkin(11) * t496 + t402;
t500 = t327 * t469 - t347 * t474;
t494 = -t308 * t474 + t326 * t546 + t496 * t500;
t313 = t327 * t474 + t347 * t469;
t507 = t308 * t469 + t313 * t496 + t326 * t545;
t405 = t475 * t438 + t439 * t470;
t351 = -qJD(5) * t405 - t411 * t475 - t412 * t470;
t406 = -t438 * t470 + t439 * t475;
t421 = pkin(4) * t438 + t460;
t358 = pkin(5) * t405 - pkin(11) * t406 + t421;
t383 = -pkin(10) * t439 - t446 * t471 - t572;
t384 = -pkin(10) * t438 - t495;
t360 = t383 * t470 + t384 * t475;
t499 = t383 * t475 - t384 * t470;
t565 = -qJD(5) * t499 + t470 * t602 - t475 * t603;
t591 = -(qJD(6) * t347 + t307) * t405 + t308 * t406 + t326 * t351 - (-qJD(6) * t358 + t565) * t540 - t360 * t346;
t588 = pkin(4) * t433;
t586 = qJD(2) * pkin(2);
t583 = t326 * t406;
t581 = t358 * t346;
t573 = t426 * t433;
t571 = t467 * t473;
t570 = t470 * t471;
t569 = t471 * t475;
t479 = qJD(3) ^ 2;
t568 = t472 * t479;
t567 = t477 * t479;
t564 = qJD(5) * t360 + t470 * t603 + t475 * t602;
t518 = -t415 * t471 - t409;
t365 = t518 - t587;
t558 = t476 * t415 - t407;
t366 = t427 + t558;
t459 = pkin(3) * t476 + pkin(4);
t562 = t365 * t470 + t366 * t475 - t459 * t547 - (-t471 * t548 + (t475 * t476 - t570) * qJD(4)) * pkin(3);
t559 = t365 * t475 - t366 * t470 + t459 * t548 + (t471 * t547 + (t470 * t476 + t569) * qJD(4)) * pkin(3);
t461 = pkin(3) * t552;
t428 = qJD(3) * t461 + t473 * t530;
t551 = qJD(2) * t473;
t544 = qJD(6) * t478;
t533 = t478 * t553;
t528 = -pkin(3) * t464 - t410;
t527 = -pkin(4) * t463 - t356;
t353 = t367 - t588;
t425 = pkin(3) * t569 + t459 * t470 + pkin(11);
t512 = qJD(6) * t425 + t353 + t461;
t457 = pkin(4) * t470 + pkin(11);
t511 = qJD(6) * t457 + t353;
t374 = pkin(4) * t400 + t428;
t333 = t363 * t470 + t579;
t509 = pkin(4) * t548 - t333;
t334 = t363 * t475 - t580;
t508 = -pkin(4) * t547 + t334;
t352 = qJD(5) * t406 - t411 * t470 + t475 * t412;
t506 = pkin(5) * t352 - pkin(11) * t351 + t505;
t502 = -t346 * t425 - t584;
t501 = -t346 * t457 - t584;
t429 = t468 * t477 - t472 * t571;
t430 = t468 * t472 + t477 * t571;
t386 = t429 * t476 - t430 * t471;
t387 = t429 * t471 + t430 * t476;
t498 = t386 * t475 - t387 * t470;
t362 = t386 * t470 + t387 * t475;
t493 = -t426 * t432 - t521;
t492 = t351 * t474 - t406 * t546;
t488 = -0.2e1 * qJD(3) * t586;
t481 = t433 * t432 * MDP(12) + t399 * MDP(14) + (-t432 ^ 2 + t433 ^ 2) * MDP(13) + (-t432 * MDP(14) + (-qJD(2) * t439 - t433) * MDP(15)) * t464 + t606;
t480 = qJD(2) ^ 2;
t458 = -pkin(4) * t475 - pkin(5);
t424 = pkin(3) * t570 - t459 * t475 - pkin(5);
t419 = t461 - t588;
t414 = -qJD(3) * t430 - t472 * t533;
t413 = qJD(3) * t429 + t477 * t533;
t349 = -qJD(4) * t387 - t413 * t471 + t414 * t476;
t348 = qJD(4) * t386 + t413 * t476 + t414 * t471;
t317 = pkin(5) * t346 - pkin(11) * t345 + t374;
t316 = t474 * t317;
t315 = qJD(5) * t362 + t348 * t470 - t349 * t475;
t314 = qJD(5) * t498 + t348 * t475 + t349 * t470;
t1 = [(-(-t314 * t469 - t362 * t545) * t540 - t362 * t342 + t315 * t375 - t498 * t336) * MDP(31) + ((t314 * t474 - t362 * t546) * t540 - t362 * t344 + t315 * t377 - t498 * t335) * MDP(32) + (MDP(17) * t349 - MDP(18) * t348) * t464 + (-MDP(24) * t315 - MDP(25) * t314) * t463 + (MDP(10) * t414 - MDP(11) * t413) * qJD(3) + ((-t400 * t478 - t432 * t551) * MDP(17) + (-t399 * t478 - t433 * t551) * MDP(18) + (-t346 * t478 - t394 * t551) * MDP(24) + (-t345 * t478 + t496 * t551) * MDP(25) + (-(t469 * t544 + t474 * t551) * t540 - t478 * t344) * MDP(31) + ((t469 * t551 - t474 * t544) * t540 + t478 * t342) * MDP(32) - t594 * t478 * t539 + (-t478 * MDP(4) + (-MDP(10) * t477 + MDP(11) * t472 - MDP(3)) * t473) * t480) * t467; (t332 * t406 + t377 * t492) * MDP(26) + ((-t375 * t474 - t377 * t469) * t351 + (-t331 - t336 * t474 + (t375 * t469 - t377 * t474) * qJD(6)) * t406) * MDP(27) + (-t313 * t352 - t499 * t335 + t564 * t377 + (-t581 - (-qJD(6) * t327 + t317) * t405 - qJD(6) * t583 - (qJD(6) * t360 - t506) * t540) * t469 + t591 * t474) * MDP(32) + (-t500 * t352 + t316 * t405 - t499 * t336 + t564 * t375 + (t581 - t506 * t540 + (-t327 * t405 + t360 * t540 + t583) * qJD(6)) * t474 + t591 * t469) * MDP(31) + (-t406 * t342 - t336 * t405 - t352 * t375 - (-t351 * t469 - t406 * t545) * t540) * MDP(29) + (t460 * t399 - t426 * t411 + t428 * t439 + t489 * t433) * MDP(18) + (t335 * t405 + t344 * t406 + t352 * t377 - t492 * t540) * MDP(28) - MDP(8) * t568 + (pkin(8) * t568 + t477 * t488) * MDP(11) + (t346 * t421 + t352 * t402 + t374 * t405 - t394 * t505) * MDP(24) + (t345 * t421 + t351 * t402 + t374 * t406 + t496 * t505) * MDP(25) + (-pkin(8) * t567 + t472 * t488) * MDP(10) + MDP(7) * t567 - 0.2e1 * t539 * t596 + 0.2e1 * t529 * t597 + (t346 * t405 - t352 * t540) * MDP(30) + (t345 * t406 + t351 * t496) * MDP(19) + (-t345 * t405 - t346 * t406 + t351 * t394 - t352 * t496) * MDP(20) + (t399 * t439 + t411 * t433) * MDP(12) + (-t399 * t438 - t400 * t439 - t411 * t432 + t412 * t433) * MDP(13) + (t460 * t400 + t426 * t412 + t428 * t438 + t489 * t432) * MDP(17) + (-t411 * MDP(14) - t412 * MDP(15) + MDP(17) * t608 + MDP(18) * t607) * t464 + (MDP(21) * t351 - MDP(22) * t352 - MDP(24) * t564 + MDP(25) * t565) * t463; t481 + (t424 * t336 + t502 * t469 + t559 * t375 - (t469 * t562 - t474 * t512) * t540 + t494) * MDP(31) + (t424 * t335 + t502 * t474 + t559 * t377 - (t469 * t512 + t474 * t562) * t540 + t507) * MDP(32) + (-t419 * t496 + t463 * t562 + t609) * MDP(25) + (t394 * t419 - t463 * t559 + t601) * MDP(24) + (t558 * t464 + t433 * t461 + (qJD(4) * t528 - t380) * t476 + t493) * MDP(18) + (-t518 * t464 + t432 * t461 + t573 + (t471 * t528 - t409) * qJD(4) + t522) * MDP(17) + t594 * qJD(2) * t586 + (-t477 * t597 + t596) * t480; (t496 * t588 + t334 * t463 - t576 + (qJD(5) * t527 - t322) * t475 - t523) * MDP(25) + t481 + (t458 * t336 + t501 * t469 + t509 * t375 - (t469 * t508 - t474 * t511) * t540 + t494) * MDP(31) + (-t394 * t588 + t333 * t463 - t578 + (t470 * t527 - t579) * qJD(5) - t524) * MDP(24) + (t458 * t335 + t501 * t474 + t509 * t377 - (t469 * t511 + t474 * t508) * t540 + t507) * MDP(32) + (-t464 * t497 + t485 + t573) * MDP(17) + (t464 * t519 + t493 - t595) * MDP(18); (t329 * t463 + t601) * MDP(24) + (t328 * t463 + t609) * MDP(25) + (-pkin(5) * t336 + (-t328 * t469 + t367 * t474) * t540 - t329 * t375 - t469 * t584 - t561 * pkin(11) + t494) * MDP(31) + (-pkin(5) * t335 - (t328 * t474 + t367 * t469) * t540 - t329 * t377 - t326 * t605 + (-t540 * t546 - t344) * pkin(11) + t507) * MDP(32) + t606; t377 * t375 * MDP(26) + (-t375 ^ 2 + t377 ^ 2) * MDP(27) + (-t375 * t540 + t560) * MDP(28) + (-t377 * t540 - t582) * MDP(29) + t346 * MDP(30) + (-t307 * t469 - t313 * t540 - t326 * t377 + t316) * MDP(31) + (-t307 * t474 - t317 * t469 + t326 * t375 + t500 * t540) * MDP(32) + (-MDP(28) * t574 - MDP(29) * t377 - MDP(31) * t313 + MDP(32) * t500) * qJD(6);];
tauc  = t1;
