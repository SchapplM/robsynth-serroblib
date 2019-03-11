% Calculate vector of inverse dynamics joint torques for
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:43
% EndTime: 2019-03-08 19:01:51
% DurationCPUTime: 6.41s
% Computational Cost: add. (4182->442), mult. (10251->644), div. (0->0), fcn. (9436->18), ass. (0->212)
t445 = qJD(4) + qJD(5);
t455 = sin(qJ(5));
t459 = cos(qJ(4));
t585 = cos(qJ(5));
t531 = qJD(3) * t585;
t456 = sin(qJ(4));
t553 = qJD(3) * t456;
t594 = -t455 * t553 + t459 * t531;
t595 = t445 * t594;
t586 = pkin(9) + pkin(10);
t443 = qJDD(4) + qJDD(5);
t580 = cos(pkin(6));
t430 = qJD(1) * t580 + qJD(2);
t450 = sin(pkin(7));
t453 = cos(pkin(7));
t452 = cos(pkin(13));
t451 = sin(pkin(6));
t554 = qJD(1) * t451;
t535 = t452 * t554;
t396 = t430 * t453 - t450 * t535;
t449 = sin(pkin(13));
t457 = sin(qJ(3));
t460 = cos(qJ(3));
t563 = t452 * t453;
t492 = t449 * t460 + t457 * t563;
t480 = t492 * t451;
t565 = t450 * t457;
t372 = qJD(1) * t480 + t430 * t565;
t522 = t586 * qJD(3) + t372;
t345 = t396 * t456 + t459 * t522;
t423 = t580 * qJDD(1) + qJDD(2);
t545 = qJDD(1) * t451;
t527 = t452 * t545;
t395 = t423 * t453 - t450 * t527;
t386 = t459 * t395;
t540 = t460 * t563;
t566 = t449 * t457;
t493 = t540 - t566;
t551 = qJD(3) * t460;
t339 = qJDD(3) * pkin(9) + (t423 * t457 + t430 * t551) * t450 + (qJD(1) * qJD(3) * t493 + qJDD(1) * t492) * t451;
t521 = pkin(10) * qJDD(3) + t339;
t308 = qJDD(4) * pkin(4) - t345 * qJD(4) - t521 * t456 + t386;
t344 = t459 * t396 - t522 * t456;
t310 = t344 * qJD(4) + t456 * t395 + t521 * t459;
t343 = qJD(4) * pkin(4) + t344;
t537 = t585 * t345;
t317 = t455 * t343 + t537;
t587 = qJD(5) * t317 - t585 * t308 + t455 * t310;
t300 = -t443 * pkin(5) + t587;
t579 = cos(pkin(12));
t505 = t580 * t579;
t578 = sin(pkin(12));
t398 = t449 * t505 + t452 * t578;
t472 = t449 * t578 - t452 * t505;
t524 = t451 * t579;
t590 = t450 * t524 + t472 * t453;
t355 = t398 * t460 - t590 * t457;
t504 = t580 * t578;
t399 = -t449 * t504 + t452 * t579;
t473 = t449 * t579 + t452 * t504;
t523 = t451 * t578;
t589 = -t450 * t523 + t473 * t453;
t357 = t399 * t460 - t589 * t457;
t525 = t450 * t580;
t379 = t457 * t525 + t480;
t380 = t450 * t472 - t453 * t524;
t381 = t450 * t473 + t453 * t523;
t397 = -t451 * t452 * t450 + t453 * t580;
t448 = qJ(4) + qJ(5);
t441 = sin(t448);
t442 = cos(t448);
t486 = -g(3) * (-t379 * t441 + t397 * t442) - g(2) * (-t355 * t441 + t380 * t442) - g(1) * (-t357 * t441 + t381 * t442);
t479 = -t300 + t486;
t513 = t453 * t535;
t552 = qJD(3) * t457;
t534 = t450 * t552;
t536 = t449 * t554;
t541 = t451 * t566;
t468 = -t460 * (t423 * t450 + t453 * t527) + qJDD(1) * t541 + t430 * t534 + t513 * t552 + t536 * t551;
t354 = t398 * t457 + t590 * t460;
t356 = t399 * t457 + t589 * t460;
t512 = t460 * t525;
t378 = -t451 * t540 - t512 + t541;
t485 = g(1) * t356 + g(2) * t354 + g(3) * t378;
t593 = qJD(3) * t372 - t468 + t485;
t560 = t455 * t459;
t409 = t456 * t585 + t560;
t383 = t445 * t409;
t526 = qJDD(3) * t585;
t544 = qJDD(3) * t456;
t503 = t455 * t544 - t459 * t526;
t363 = qJD(3) * t383 + t503;
t361 = qJDD(6) + t363;
t439 = -pkin(4) * t459 - pkin(3);
t488 = -t455 * t456 + t459 * t585;
t377 = -pkin(5) * t488 - pkin(11) * t409 + t439;
t477 = t485 * t442;
t592 = t377 * t361 + t477;
t550 = qJD(4) * t456;
t509 = pkin(4) * t550 - t372;
t358 = -t379 * t456 + t397 * t459;
t359 = t379 * t459 + t397 * t456;
t323 = t455 * t358 + t359 * t585;
t458 = cos(qJ(6));
t375 = t378 * t458;
t454 = sin(qJ(6));
t591 = -t323 * t454 + t375;
t530 = qJD(5) * t585;
t549 = qJD(5) * t455;
t469 = t455 * t308 + t310 * t585 + t343 * t530 - t345 * t549;
t299 = t443 * pkin(11) + t469;
t561 = t455 * t345;
t316 = t343 * t585 - t561;
t314 = -t445 * pkin(5) - t316;
t371 = -t457 * t536 + t460 * (t430 * t450 + t513);
t365 = qJD(3) * t439 - t371;
t405 = -qJD(3) * t560 - t456 * t531;
t342 = -pkin(5) * t594 + pkin(11) * t405 + t365;
t382 = t445 * t488;
t418 = t586 * t456;
t419 = t586 * t459;
t392 = -t455 * t418 + t419 * t585;
t400 = qJD(6) - t594;
t507 = g(1) * t357 + g(2) * t355;
t484 = g(3) * t379 + t507;
t538 = qJD(4) * t586;
t413 = t456 * t538;
t414 = t459 * t538;
t489 = -t418 * t585 - t455 * t419;
t558 = -qJD(5) * t489 + t488 * t371 + t413 * t585 + t455 * t414;
t588 = (qJD(6) * t342 + t299) * t488 + t300 * t409 + t314 * t382 + (-qJD(6) * t377 + t558) * t400 - t392 * t361 - t484;
t581 = qJD(3) * pkin(3);
t576 = t314 * t594;
t575 = t314 * t409;
t543 = qJDD(3) * t459;
t362 = t455 * t543 + t456 * t526 + t595;
t547 = qJD(6) * t458;
t539 = t458 * t362 + t454 * t443 + t445 * t547;
t548 = qJD(6) * t454;
t335 = t405 * t548 + t539;
t573 = t335 * t454;
t571 = t378 * t454;
t567 = t405 * t454;
t387 = -t458 * t445 - t567;
t570 = t387 * t400;
t498 = t405 * t458 - t445 * t454;
t569 = t498 * t400;
t564 = t450 * t460;
t562 = t454 * t361;
t559 = t458 * t361;
t557 = pkin(5) * t383 - pkin(11) * t382 + t509;
t556 = qJD(5) * t392 - t409 * t371 - t455 * t413 + t414 * t585;
t446 = t456 ^ 2;
t555 = -t459 ^ 2 + t446;
t546 = qJD(3) * qJD(4);
t533 = t450 * t551;
t529 = t456 * t546;
t528 = t460 * t546;
t315 = t445 * pkin(11) + t317;
t500 = t315 * t454 - t342 * t458;
t520 = t314 * t548 - t405 * t500;
t518 = t362 * t454 - t458 * t443;
t517 = t400 * t458;
t376 = -pkin(5) * t405 - pkin(11) * t594;
t437 = pkin(4) * t455 + pkin(11);
t515 = pkin(4) * t553 + qJD(6) * t437 + t376;
t318 = t455 * t344 + t537;
t508 = t549 * pkin(4) - t318;
t506 = -g(1) * t381 - g(2) * t380;
t501 = -t437 * t361 - t576;
t305 = t315 * t458 + t342 * t454;
t499 = t323 * t458 + t571;
t462 = qJD(3) ^ 2;
t497 = qJDD(3) * t460 - t457 * t462;
t319 = t344 * t585 - t561;
t496 = -t530 * pkin(4) + t319;
t401 = t453 * t459 - t456 * t565;
t402 = t453 * t456 + t459 * t565;
t370 = t455 * t401 + t402 * t585;
t495 = -t370 * t454 - t458 * t564;
t494 = -t370 * t458 + t454 * t564;
t491 = t358 * t585 - t455 * t359;
t490 = t401 * t585 - t455 * t402;
t487 = t382 * t458 - t409 * t548;
t481 = -t305 * t405 + t314 * t547 - t479 * t454;
t366 = -t371 - t581;
t478 = -qJD(3) * t366 - t339 + t507;
t476 = -pkin(9) * qJDD(4) + (t366 + t371 - t581) * qJD(4);
t461 = qJD(4) ^ 2;
t467 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t461 + t593;
t336 = -qJD(6) * t498 + t518;
t466 = ((t335 - t570) * t458 + (-t336 + t569) * t454) * MDP(21) + (-t498 * t517 + t573) * MDP(20) + (-t400 ^ 2 * t454 - t387 * t405 + t559) * MDP(23) + (t400 * t517 - t405 * t498 + t562) * MDP(22) + (t362 - t595) * MDP(15) + (-t503 + (-qJD(3) * t409 - t405) * t445) * MDP(16) + (t405 ^ 2 - t594 ^ 2) * MDP(14) + t443 * MDP(17) + (MDP(13) * t594 + t400 * MDP(24)) * t405;
t334 = pkin(4) * t529 + qJDD(3) * t439 + t468;
t464 = t365 * t405 + t486 - t587;
t329 = t355 * t442 + t380 * t441;
t331 = t357 * t442 + t381 * t441;
t351 = t379 * t442 + t397 * t441;
t463 = g(1) * t331 + g(2) * t329 + g(3) * t351 - t365 * t594 - t469;
t438 = -t585 * pkin(4) - pkin(5);
t385 = -qJD(4) * t402 - t456 * t533;
t384 = qJD(4) * t401 + t459 * t533;
t374 = t379 * qJD(3);
t373 = (t451 * t493 + t512) * qJD(3);
t333 = qJD(4) * t358 + t373 * t459;
t332 = -qJD(4) * t359 - t373 * t456;
t325 = qJD(5) * t370 + t455 * t384 - t385 * t585;
t324 = qJD(5) * t490 + t384 * t585 + t455 * t385;
t309 = pkin(5) * t363 - pkin(11) * t362 + t334;
t307 = t458 * t309;
t303 = qJD(5) * t323 - t332 * t585 + t455 * t333;
t302 = qJD(5) * t491 + t455 * t332 + t333 * t585;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t423 * t580 - g(3) + (t449 ^ 2 + t452 ^ 2) * t451 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(3) * t374 - qJDD(3) * t378) * MDP(4) + (-qJD(3) * t373 - qJDD(3) * t379) * MDP(5) + (-t378 * t543 + qJD(4) * t332 + qJDD(4) * t358 + (-t374 * t459 + t378 * t550) * qJD(3)) * MDP(11) + (t378 * t544 - qJD(4) * t333 - qJDD(4) * t359 + (qJD(4) * t378 * t459 + t374 * t456) * qJD(3)) * MDP(12) + (-t303 * t445 + t363 * t378 - t374 * t594 + t443 * t491) * MDP(18) + (-t302 * t445 - t323 * t443 + t362 * t378 - t374 * t405) * MDP(19) + ((-qJD(6) * t499 - t302 * t454 + t374 * t458) * t400 + t591 * t361 + t303 * t387 - t491 * t336) * MDP(25) + (-(t591 * qJD(6) + t302 * t458 + t374 * t454) * t400 - t499 * t361 - t303 * t498 - t491 * t335) * MDP(26); (-g(3) * t580 + (-g(1) * t578 + g(2) * t579) * t451 + t423) * MDP(2) + (qJD(4) * t385 + qJDD(4) * t401) * MDP(11) + (-qJD(4) * t384 - qJDD(4) * t402) * MDP(12) + (-t325 * t445 + t443 * t490) * MDP(18) + (-t324 * t445 - t370 * t443) * MDP(19) + ((qJD(6) * t494 - t324 * t454 + t458 * t534) * t400 + t495 * t361 + t325 * t387 - t490 * t336) * MDP(25) + (-(qJD(6) * t495 + t324 * t458 + t454 * t534) * t400 + t494 * t361 - t325 * t498 - t490 * t335) * MDP(26) + (t497 * MDP(4) + (-qJDD(3) * t457 - t460 * t462) * MDP(5) + (-t456 * t528 + t459 * t497) * MDP(11) + (-t456 * t497 - t459 * t528) * MDP(12) + (-t363 * t460 - t552 * t594) * MDP(18) + (-t362 * t460 - t405 * t552) * MDP(19)) * t450; qJDD(3) * MDP(3) + t593 * MDP(4) + (-t423 * t565 - t492 * t545 + (-t430 * t564 - t493 * t554 + t371) * qJD(3) + t484) * MDP(5) + (qJDD(3) * t446 + 0.2e1 * t459 * t529) * MDP(6) + 0.2e1 * (t456 * t543 - t546 * t555) * MDP(7) + (qJDD(4) * t456 + t459 * t461) * MDP(8) + (qJDD(4) * t459 - t456 * t461) * MDP(9) + (t456 * t476 + t459 * t467) * MDP(11) + (-t456 * t467 + t459 * t476) * MDP(12) + (t362 * t409 - t382 * t405) * MDP(13) + (t362 * t488 - t363 * t409 + t382 * t594 + t383 * t405) * MDP(14) + (t382 * t445 + t409 * t443) * MDP(15) + (-t383 * t445 + t443 * t488) * MDP(16) + (-t334 * t488 + t363 * t439 + t365 * t383 + t443 * t489 - t445 * t556 - t509 * t594 + t477) * MDP(18) + (t334 * t409 + t362 * t439 + t365 * t382 - t392 * t443 - t405 * t509 - t441 * t485 + t445 * t558) * MDP(19) + (t335 * t409 * t458 - t487 * t498) * MDP(20) + ((-t387 * t458 + t454 * t498) * t382 + (-t573 - t336 * t458 + (t387 * t454 + t458 * t498) * qJD(6)) * t409) * MDP(21) + (-t335 * t488 - t383 * t498 + t400 * t487 + t409 * t559) * MDP(22) + (-t409 * t562 + t336 * t488 - t383 * t387 + (-t382 * t454 - t409 * t547) * t400) * MDP(23) + (-t361 * t488 + t383 * t400) * MDP(24) + (-t500 * t383 - t307 * t488 - t489 * t336 + t556 * t387 + (t557 * t400 + (t315 * t488 - t392 * t400 + t575) * qJD(6) + t592) * t458 + t588 * t454) * MDP(25) + (-t305 * t383 - t489 * t335 - t556 * t498 + ((-qJD(6) * t315 + t309) * t488 - qJD(6) * t575 + (qJD(6) * t392 - t557) * t400 - t592) * t454 + t588 * t458) * MDP(26); (t319 * t445 + (t405 * t553 - t443 * t455 - t445 * t530) * pkin(4) + t463) * MDP(19) + (t318 * t445 + (t443 * t585 - t445 * t549 + t553 * t594) * pkin(4) + t464) * MDP(18) + t466 + MDP(9) * t543 + MDP(8) * t544 + (-g(3) * t358 + t456 * t478 + t459 * t506 + t386) * MDP(11) + (g(3) * t359 + (-t395 - t506) * t456 + t478 * t459) * MDP(12) + qJDD(4) * MDP(10) + (t438 * t336 + t508 * t387 + (t400 * t496 + t501) * t454 + (-t400 * t515 + t479) * t458 + t520) * MDP(25) + (t438 * t335 + t501 * t458 - t508 * t498 + (t454 * t515 + t496 * t458) * t400 + t481) * MDP(26) + (-MDP(6) * t456 * t459 + MDP(7) * t555) * t462; (t316 * t445 + t463) * MDP(19) + (t317 * t445 + t464) * MDP(18) + t466 + (-pkin(5) * t335 + (t316 * t458 + t376 * t454) * t400 + t317 * t498 - t458 * t576 + (t400 * t548 - t559) * pkin(11) + t481) * MDP(26) + (-pkin(5) * t336 - t317 * t387 + (-pkin(11) * t361 + t316 * t400 - t576) * t454 + ((-pkin(11) * qJD(6) - t376) * t400 + t479) * t458 + t520) * MDP(25); -t498 * t387 * MDP(20) + (-t387 ^ 2 + t498 ^ 2) * MDP(21) + (t539 + t570) * MDP(22) + (-t518 - t569) * MDP(23) + t361 * MDP(24) + (-t454 * t299 + t307 + t305 * t400 + t314 * t498 - g(1) * (-t331 * t454 + t356 * t458) - g(2) * (-t329 * t454 + t354 * t458) - g(3) * (-t351 * t454 + t375)) * MDP(25) + (-t458 * t299 - t454 * t309 - t500 * t400 + t314 * t387 - g(1) * (-t331 * t458 - t356 * t454) - g(2) * (-t329 * t458 - t354 * t454) - g(3) * (-t351 * t458 - t571)) * MDP(26) + (MDP(22) * t567 + MDP(23) * t498 - MDP(25) * t305 + MDP(26) * t500) * qJD(6);];
tau  = t1;
