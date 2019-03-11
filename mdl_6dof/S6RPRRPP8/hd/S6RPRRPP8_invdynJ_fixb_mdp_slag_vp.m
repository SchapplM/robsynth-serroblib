% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:56:05
% EndTime: 2019-03-09 04:56:14
% DurationCPUTime: 6.60s
% Computational Cost: add. (3984->572), mult. (7334->672), div. (0->0), fcn. (4315->6), ass. (0->230)
t479 = sin(qJ(4));
t652 = t479 * qJ(5) + pkin(3);
t483 = cos(qJ(3));
t573 = qJD(1) * qJD(3);
t550 = t483 * t573;
t480 = sin(qJ(3));
t571 = qJDD(1) * t480;
t413 = qJDD(4) + t550 + t571;
t634 = pkin(4) + qJ(6);
t552 = t634 * t413;
t579 = qJD(4) * t483;
t555 = t479 * t579;
t482 = cos(qJ(4));
t577 = t482 * qJD(3);
t506 = t480 * t577 + t555;
t472 = t483 * pkin(8);
t649 = pkin(3) * t480 + qJ(2) - t472;
t397 = t649 * qJD(1);
t485 = -pkin(1) - pkin(7);
t444 = qJD(1) * t485 + qJD(2);
t423 = t480 * t444;
t404 = qJD(3) * pkin(8) + t423;
t357 = -t482 * t397 + t404 * t479;
t586 = qJD(3) * t479;
t587 = qJD(1) * t483;
t417 = t482 * t587 + t586;
t515 = pkin(5) * t417 + t357;
t576 = qJD(5) + t515;
t567 = MDP(23) + MDP(26);
t403 = t413 * qJ(5);
t589 = qJD(1) * t480;
t449 = qJD(4) + t589;
t435 = t449 * qJD(5);
t651 = -t403 - t435;
t585 = qJD(3) * t480;
t560 = t479 * t585;
t570 = qJDD(1) * t483;
t582 = qJD(4) * t417;
t361 = -qJD(1) * t560 - t482 * qJDD(3) + t479 * t570 + t582;
t581 = qJD(4) * t479;
t588 = qJD(1) * t482;
t629 = qJ(5) * t480;
t648 = pkin(4) * t581 - qJD(5) * t479 - t588 * t629 - t423;
t568 = MDP(21) + MDP(25);
t647 = -pkin(5) * t361 + qJDD(6);
t360 = qJD(1) * t506 - qJD(4) * t577 - t479 * qJDD(3) - t482 * t570;
t415 = t479 * t587 - t577;
t646 = t415 ^ 2;
t412 = t417 ^ 2;
t446 = t449 ^ 2;
t645 = 0.2e1 * t403;
t644 = pkin(5) + pkin(8);
t643 = pkin(4) * t413;
t642 = pkin(4) * t482;
t641 = pkin(5) * t360;
t639 = pkin(5) * t415;
t638 = pkin(5) * t483;
t637 = pkin(8) * t413;
t481 = sin(qJ(1));
t636 = g(1) * t481;
t484 = cos(qJ(1));
t475 = g(2) * t484;
t635 = g(3) * t480;
t633 = pkin(8) * qJD(4);
t632 = pkin(1) * qJDD(1);
t631 = qJ(5) * t361;
t630 = qJ(5) * t415;
t628 = qJ(6) * t479;
t340 = -t449 * t634 + t576;
t627 = t340 * t449;
t358 = t479 * t397 + t482 * t404;
t351 = -qJ(5) * t449 - t358;
t626 = t351 * t449;
t625 = t358 * t449;
t624 = t360 * t479;
t623 = t415 * t417;
t622 = t415 * t449;
t621 = t415 * t479;
t620 = t417 * t449;
t619 = t417 * t479;
t618 = t417 * t482;
t616 = t479 * t480;
t615 = t479 * t483;
t614 = t480 * t481;
t613 = t480 * t482;
t612 = t480 * t484;
t611 = t480 * t485;
t610 = t481 * t479;
t609 = t481 * t482;
t608 = t481 * t483;
t607 = t482 * t483;
t606 = t483 * t484;
t605 = t484 * t482;
t486 = qJD(3) ^ 2;
t604 = t485 * t486;
t551 = t634 * t480;
t603 = t479 * qJD(1) * t551 - qJD(6) * t482 + (-qJ(5) * t482 + t628) * qJD(4) + t648;
t580 = qJD(4) * t482;
t602 = pkin(4) * t479 * t589 - qJ(5) * t580 + t648;
t530 = pkin(3) * t483 + pkin(8) * t480;
t420 = t530 * qJD(1);
t601 = t479 * t420 + t444 * t607;
t600 = -t644 * t581 - (pkin(5) * t616 + qJ(5) * t483) * qJD(1) - t601;
t434 = t644 * t482;
t507 = -pkin(5) * t613 - t483 * t634;
t409 = t444 * t615;
t541 = -t420 * t482 + t409;
t599 = -qJD(1) * t507 + qJD(4) * t434 - t541;
t565 = g(2) * t606;
t597 = g(3) * t616 + t479 * t565;
t596 = g(3) * t613 + t482 * t565;
t595 = t479 * t649 + t482 * t611;
t594 = -pkin(4) * t615 + qJ(5) * t607;
t593 = g(1) * t606 + g(2) * t608;
t592 = t484 * pkin(1) + t481 * qJ(2);
t477 = t483 ^ 2;
t591 = t480 ^ 2 - t477;
t487 = qJD(1) ^ 2;
t590 = -t486 - t487;
t584 = qJD(3) * t483;
t583 = qJD(3) * t485;
t578 = qJD(4) * t485;
t575 = qJD(5) + t357;
t347 = t358 - t639;
t574 = -qJD(6) - t347;
t572 = qJDD(1) * qJ(2);
t569 = MDP(19) + MDP(27);
t566 = g(1) * t608;
t460 = g(2) * t612;
t564 = 0.2e1 * qJD(1) * qJD(2);
t563 = t481 * t607;
t414 = qJD(3) * t530 + qJD(2);
t558 = t483 * t583;
t562 = t479 * t414 + t482 * t558 + t580 * t649;
t405 = -qJD(3) * pkin(3) - t444 * t483;
t505 = -qJ(5) * t417 + t405;
t343 = t415 * t634 + t505;
t557 = t343 * t581;
t556 = t343 * t580;
t554 = t479 * t578;
t553 = t480 * t580;
t547 = -t475 + t636;
t546 = g(3) * t483 - t460;
t545 = MDP(22) - t569;
t544 = MDP(20) - t567;
t398 = t480 * t610 - t605;
t399 = t479 * t484 + t480 * t609;
t543 = -t398 * pkin(4) + qJ(5) * t399;
t400 = t479 * t612 + t609;
t401 = t480 * t605 - t610;
t542 = t400 * pkin(4) - qJ(5) * t401;
t441 = t479 * t611;
t540 = t482 * t649 - t441;
t538 = t482 * t568;
t537 = qJDD(2) - t632;
t367 = qJD(1) * t414 + qJDD(1) * t649;
t437 = qJDD(1) * t485 + qJDD(2);
t373 = qJDD(3) * pkin(8) + t437 * t480 + t444 * t584;
t535 = t479 * t367 + t482 * t373 + t397 * t580 - t404 * t581;
t534 = t482 * t367 - t479 * t373 - t397 * t581 - t404 * t580;
t533 = pkin(4) * t563 + pkin(8) * t614 + t652 * t608;
t532 = t506 * qJ(5) + t480 * t583 + t579 * t642;
t372 = -qJDD(3) * pkin(3) - t437 * t483 + t444 * t585;
t494 = qJ(5) * t360 - qJD(5) * t417 + t372;
t332 = qJD(6) * t415 + t361 * t634 + t494;
t531 = -t332 - t566;
t529 = g(1) * t400 + g(2) * t398;
t528 = -g(1) * t401 - g(2) * t399;
t527 = g(1) * t484 + g(2) * t481;
t376 = -t595 - t629;
t526 = qJ(2) * t487 + t636;
t334 = -t535 + t651;
t520 = -qJDD(5) + t534;
t335 = -t520 - t643;
t525 = -t334 * t482 + t335 * t479;
t342 = qJD(6) - t351 - t639;
t524 = t340 * t482 - t342 * t479;
t523 = t340 * t479 + t342 * t482;
t350 = -pkin(4) * t449 + t575;
t522 = t350 * t482 + t351 * t479;
t521 = t350 * t479 - t351 * t482;
t519 = t652 + t642;
t517 = t414 * t482 - t479 * t558 - t485 * t553 - t581 * t649;
t516 = t564 + 0.2e1 * t572;
t514 = -t437 + t526;
t512 = t413 * t479 + t449 * t580;
t511 = t413 * t482 - t449 * t581;
t510 = -g(1) * t614 - t546;
t509 = 0.2e1 * qJ(2) * t573 + qJDD(3) * t485;
t508 = t449 * t633 + t566;
t408 = -t482 * t634 - t652;
t504 = pkin(3) * t614 + t399 * pkin(4) + t484 * pkin(7) + qJ(5) * t398 + t592;
t336 = pkin(4) * t361 + t494;
t503 = t336 + t508;
t502 = t405 * t449 - t637;
t354 = pkin(4) * t415 + t505;
t501 = -t354 * t449 + t637;
t471 = t484 * qJ(2);
t500 = pkin(3) * t612 + t401 * pkin(4) - pkin(8) * t606 + qJ(5) * t400 + t471;
t499 = g(1) * t398 - g(2) * t400 + g(3) * t615 + t534;
t498 = t516 - t527;
t496 = -qJDD(5) + t499;
t495 = t479 * t545 - t482 * t544;
t493 = g(1) * t399 - g(2) * t401 + g(3) * t607 - t535;
t492 = t354 * t417 - t496;
t345 = -t360 + t622;
t491 = t343 * t417 - t496 - t641;
t490 = t620 - t361;
t489 = -t343 * t415 - t493 + t647;
t488 = t522 * MDP(24) + t524 * MDP(28) + t568 * t621 + (t479 * t544 + t482 * t545) * t449;
t468 = qJDD(3) * t483;
t433 = t644 * t479;
t378 = -t483 * t485 - t594;
t377 = -pkin(4) * t480 - t540;
t375 = (-t485 + t628) * t483 - t594;
t370 = pkin(4) * t417 + t630;
t366 = -pkin(5) * t615 - t376;
t365 = -pkin(4) * t587 + t541;
t364 = -qJ(5) * t587 - t601;
t356 = t441 + (-t649 + t638) * t482 - t551;
t352 = t417 * t634 + t630;
t349 = -pkin(4) * t560 - qJD(5) * t607 + t532;
t344 = -pkin(4) * t584 - t517;
t341 = -qJ(5) * t584 + (-qJD(5) + t554) * t480 - t562;
t339 = (qJ(6) * qJD(4) - qJD(5)) * t607 + (-qJD(3) * t551 + qJD(6) * t483) * t479 + t532;
t338 = (-pkin(5) * t580 + qJ(5) * qJD(3)) * t483 + (qJD(5) + (pkin(5) * qJD(3) - t578) * t479) * t480 + t562;
t337 = -pkin(5) * t555 + qJD(3) * t507 - qJD(6) * t480 - t517;
t333 = -t334 + t647;
t331 = -qJD(6) * t449 - t520 - t552 - t641;
t1 = [(-qJDD(3) * t480 - t483 * t486) * MDP(10) + (-t480 * t486 + t468) * MDP(9) + qJDD(1) * MDP(1) + 0.2e1 * (-t480 * t570 + t573 * t591) * MDP(8) + (-t537 * pkin(1) - g(1) * (-pkin(1) * t481 + t471) - g(2) * t592 + (t564 + t572) * qJ(2)) * MDP(6) + (t341 * t415 + t344 * t417 - t360 * t377 + t361 * t376 - t522 * t585 + (-qJD(4) * t521 + t334 * t479 + t335 * t482) * t483 + t593) * MDP(21) + (t337 * t417 - t338 * t415 - t356 * t360 - t361 * t366 - t524 * t585 + (-qJD(4) * t523 + t331 * t482 - t333 * t479) * t483 + t593) * MDP(25) + (-t562 * t449 - t595 * t413 + (t449 * t554 + (-t405 * t482 + t417 * t485) * qJD(3) - t535) * t480 + (-qJD(3) * t358 + t360 * t485 + t372 * t482 - t405 * t581) * t483 + t529) * MDP(20) + (t332 * t375 + t343 * t339 + t331 * t356 + t340 * t337 + t333 * t366 + t342 * t338 - g(1) * (-pkin(5) * t606 + qJ(6) * t401 + t500) - g(2) * (qJ(6) * t399 + t504) + (g(2) * t483 * t644 - g(1) * t485) * t481) * MDP(28) + (t336 * t378 + t354 * t349 + t334 * t376 + t351 * t341 + t335 * t377 + t350 * t344 - g(1) * (t481 * t485 + t500) - g(2) * (-pkin(8) * t608 + t504)) * MDP(24) + (-t337 * t449 + t339 * t415 - t356 * t413 + t361 * t375 + (-t343 * t586 - t331) * t480 + (-qJD(3) * t340 + t332 * t479 + t556) * t483 + t528) * MDP(27) + ((t449 * t586 - t361) * t480 + (-qJD(3) * t415 - t512) * t483) * MDP(17) + (t344 * t449 - t349 * t415 - t361 * t378 + t377 * t413 + (t354 * t586 + t335) * t480 + (qJD(3) * t350 - t336 * t479 - t354 * t580) * t483 - t528) * MDP(22) + (-t341 * t449 - t349 * t417 + t360 * t378 - t376 * t413 + (t354 * t577 - t334) * t480 + (-qJD(3) * t351 - t336 * t482 + t354 * t581) * t483 - t529) * MDP(23) + (t413 * t480 + t449 * t584) * MDP(18) + (t517 * t449 + t540 * t413 + ((-t405 * t479 + t415 * t485) * qJD(3) + t534) * t480 + (-qJD(3) * t357 - t361 * t485 + t372 * t479 + t405 * t580) * t483 + t528) * MDP(19) + ((-t449 * t577 - t360) * t480 + (qJD(3) * t417 + t511) * t483) * MDP(16) + (t338 * t449 - t339 * t417 + t360 * t375 + t366 * t413 + (t343 * t577 + t333) * t480 + (qJD(3) * t342 - t332 * t482 + t557) * t483 - t529) * MDP(26) + t498 * MDP(5) + ((t415 * t482 + t619) * t585 + (t624 - t361 * t482 + (-t618 + t621) * qJD(4)) * t483) * MDP(15) + (qJDD(2) - t547 - 0.2e1 * t632) * MDP(4) + (qJDD(1) * t477 - 0.2e1 * t480 * t550) * MDP(7) + t547 * MDP(2) + t527 * MDP(3) + (-t509 * t480 + (t516 - t604) * t483 - t593) * MDP(13) + (t509 * t483 + (t498 - t604) * t480) * MDP(12) + (-t360 * t607 - t417 * t506) * MDP(14); qJDD(1) * MDP(4) - t487 * MDP(5) + (t475 - t526 + t537) * MDP(6) + t468 * MDP(12) + t488 * qJD(1) + (t590 * MDP(13) - t336 * MDP(24) - t332 * MDP(28) + t545 * t361 + t544 * t360 + (MDP(24) * t521 + MDP(28) * t523 - t415 * t538 + t449 * t495) * qJD(3)) * t483 + (t590 * MDP(12) - qJDD(3) * MDP(13) + t525 * MDP(24) + (t331 * t479 + t333 * t482) * MDP(28) - t361 * t538 + (MDP(20) * t417 - MDP(22) * t415 + MDP(24) * t354 + MDP(28) * t343) * qJD(3) + t495 * t413 + t488 * qJD(4)) * t480 + (t415 * t569 - t417 * t567) * t585 + (-MDP(24) - MDP(28)) * t547 + t568 * (-t360 * t616 + t584 * t619 + (t553 + t588) * t417); MDP(9) * t570 - MDP(10) * t571 + qJDD(3) * MDP(11) + (t635 + (-t514 + t475) * t483) * MDP(12) + (t480 * t514 + t546) * MDP(13) + (t449 * t618 - t624) * MDP(14) + ((-t360 - t622) * t482 + (-t361 - t620) * t479) * MDP(15) + ((-t417 * t483 + t449 * t613) * qJD(1) + t512) * MDP(16) + ((t415 * t483 - t449 * t616) * qJD(1) + t511) * MDP(17) - t449 * MDP(18) * t587 + (t357 * t587 - t415 * t423 - pkin(3) * t361 + t409 * t449 + (-t566 - t372 + (-t420 - t633) * t449) * t482 + t502 * t479 + t596) * MDP(19) + (pkin(3) * t360 - t417 * t423 + t601 * t449 + t358 * t587 + t502 * t482 + (t372 + t508) * t479 - t597) * MDP(20) + (-t364 * t415 - t365 * t417 + (-t334 + t449 * t350 + (-t361 + t582) * pkin(8)) * t482 + (t335 + t626 + (qJD(4) * t415 - t360) * pkin(8)) * t479 + t510) * MDP(21) + (-t350 * t587 + t361 * t519 - t365 * t449 - t415 * t602 + t479 * t501 + t482 * t503 - t596) * MDP(22) + (t351 * t587 - t360 * t519 + t364 * t449 - t417 * t602 - t479 * t503 + t482 * t501 + t597) * MDP(23) + (-t351 * t364 - t350 * t365 - g(1) * t533 - g(3) * t472 + t602 * t354 + (qJD(4) * t522 + t460 + t525) * pkin(8) + (-t336 + t635 + t565) * t519) * MDP(24) + (-t360 * t433 - t361 * t434 + t599 * t417 - t600 * t415 + (t333 + t627) * t482 + (-t342 * t449 + t331) * t479 + t510) * MDP(25) + (-t556 + t360 * t408 + t413 * t434 + t531 * t479 + t600 * t449 - t603 * t417 + (-t342 * t483 - t343 * t613) * qJD(1) + t597) * MDP(26) + (t557 + t361 * t408 - t413 * t433 + t531 * t482 - t599 * t449 + t603 * t415 + (t340 * t483 + t343 * t616) * qJD(1) + t596) * MDP(27) + (t332 * t408 + t331 * t433 + t333 * t434 - g(1) * (qJ(6) * t563 + t533) - g(3) * (t472 + t638) + t603 * t343 + t600 * t342 + t599 * t340 + (-pkin(5) * t636 - g(3) * t408) * t480 - (t408 * t483 - t480 * t644) * t475) * MDP(28) + (MDP(7) * t480 * t483 - MDP(8) * t591) * t487; MDP(14) * t623 + (t412 - t646) * MDP(15) + t345 * MDP(16) + t490 * MDP(17) + t413 * MDP(18) + (-t405 * t417 + t499 + t625) * MDP(19) + (-t357 * t449 + t405 * t415 + t493) * MDP(20) + (pkin(4) * t360 - t631 + (-t351 - t358) * t417 + (t350 - t575) * t415) * MDP(21) + (t370 * t415 + t492 - t625 - 0.2e1 * t643) * MDP(22) + (-t354 * t415 + t370 * t417 + t449 * t575 + t435 - t493 + t645) * MDP(23) + (-t335 * pkin(4) - g(1) * t543 - g(2) * t542 - g(3) * t594 - t334 * qJ(5) - t350 * t358 - t351 * t575 - t354 * t370) * MDP(24) + (-t631 + t360 * t634 + (t342 + t574) * t417 + (t340 - t576) * t415) * MDP(25) + (t352 * t417 + t449 * t515 + 0.2e1 * t435 + t489 + t645) * MDP(26) + (-t352 * t415 + (0.2e1 * qJD(6) + t347) * t449 + 0.2e1 * t552 - t491) * MDP(27) + (-t331 * t634 + t333 * qJ(5) - t343 * t352 - g(1) * (-qJ(6) * t398 + t543) - g(2) * (qJ(6) * t400 + t542) - g(3) * (-qJ(6) * t615 + t594) + t576 * t342 + t574 * t340) * MDP(28); (t492 + t626 - t643) * MDP(24) + ((-qJD(6) - t342) * t449 - t552 + t491) * MDP(28) + (MDP(22) - MDP(27)) * (t413 - t623) + t568 * t345 + t567 * (-t412 - t446); t490 * MDP(25) + (t413 + t623) * MDP(26) + (-t446 - t646) * MDP(27) + (t489 + t627 - t651) * MDP(28);];
tau  = t1;
