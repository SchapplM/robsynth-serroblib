% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRPR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:24
% EndTime: 2019-03-09 05:37:38
% DurationCPUTime: 9.34s
% Computational Cost: add. (3870->578), mult. (7406->734), div. (0->0), fcn. (4759->8), ass. (0->237)
t498 = sin(qJ(3));
t644 = g(3) * t498;
t502 = cos(qJ(3));
t503 = cos(qJ(1));
t499 = sin(qJ(1));
t645 = g(1) * t499;
t665 = g(2) * t503 - t645;
t662 = t665 * t502;
t675 = t644 + t662;
t606 = qJD(1) * t498;
t476 = qJD(4) + t606;
t497 = sin(qJ(4));
t501 = cos(qJ(4));
t590 = t501 * qJD(3);
t605 = qJD(1) * t502;
t443 = t497 * t605 - t590;
t603 = qJD(3) * t497;
t445 = t501 * t605 + t603;
t496 = sin(qJ(6));
t500 = cos(qJ(6));
t392 = t443 * t496 + t445 * t500;
t586 = qJD(1) * qJD(3);
t567 = t502 * t586;
t584 = qJDD(1) * t498;
t441 = qJDD(4) + t567 + t584;
t436 = -qJDD(6) + t441;
t538 = -t500 * t443 + t445 * t496;
t674 = MDP(25) * t392 * t538 + (t392 ^ 2 - t538 ^ 2) * MDP(26) - t436 * MDP(29);
t587 = qJD(6) - t476;
t673 = t392 * t587;
t653 = pkin(4) + pkin(5);
t578 = t653 * t497;
t638 = qJ(5) * t501;
t528 = t578 - t638;
t663 = qJD(4) - qJD(6);
t505 = -pkin(1) - pkin(7);
t471 = qJD(1) * t505 + qJD(2);
t433 = -qJD(3) * pkin(3) - t471 * t502;
t518 = qJ(5) * t445 - t433;
t359 = -t443 * t653 + t518;
t618 = t503 * t501;
t621 = t499 * t497;
t426 = t498 * t621 - t618;
t620 = t499 * t501;
t427 = t497 * t503 + t498 * t620;
t375 = t426 * t500 - t427 * t496;
t619 = t501 * t502;
t625 = t497 * t502;
t418 = t496 * t619 - t500 * t625;
t623 = t498 * t503;
t428 = t497 * t623 + t620;
t429 = t498 * t618 - t621;
t563 = -t428 * t500 + t429 * t496;
t596 = qJD(4) * t502;
t569 = t497 * t596;
t519 = t590 * t498 + t569;
t583 = qJDD(1) * t502;
t380 = qJD(1) * t519 - qJD(4) * t590 - t497 * qJDD(3) - t501 * t583;
t557 = pkin(3) * t502 + pkin(8) * t498;
t442 = qJD(3) * t557 + qJD(2);
t648 = pkin(8) * t502;
t556 = pkin(3) * t498 - t648;
t455 = qJ(2) + t556;
t387 = qJD(1) * t442 + qJDD(1) * t455;
t466 = qJDD(1) * t505 + qJDD(2);
t601 = qJD(3) * t502;
t400 = qJDD(3) * pkin(8) + t466 * t498 + t471 * t601;
t425 = t455 * qJD(1);
t454 = t498 * t471;
t432 = qJD(3) * pkin(8) + t454;
t597 = qJD(4) * t501;
t599 = qJD(4) * t497;
t558 = t501 * t387 - t497 * t400 - t425 * t599 - t432 * t597;
t536 = qJDD(5) - t558;
t339 = pkin(9) * t380 - t441 * t653 + t536;
t431 = t441 * qJ(5);
t463 = t476 * qJD(5);
t575 = -t497 * t387 - t501 * t400 - t425 * t597;
t523 = -t432 * t599 - t575;
t345 = t431 + t463 + t523;
t602 = qJD(3) * t498;
t572 = t497 * t602;
t600 = qJD(4) * t445;
t381 = -qJD(1) * t572 - t501 * qJDD(3) + t497 * t583 + t600;
t341 = pkin(9) * t381 + t345;
t565 = t500 * t339 - t496 * t341;
t670 = g(1) * t375 + g(2) * t563 - g(3) * t418 + t359 * t392 - t565;
t667 = t587 * t538;
t666 = qJD(5) * t497 + t454;
t507 = qJD(1) ^ 2;
t520 = -qJ(2) * t507 + t665;
t377 = t501 * t425 - t497 * t432;
t588 = qJD(5) - t377;
t369 = pkin(4) * t443 - t518;
t649 = pkin(8) * t441;
t661 = t369 * t476 - t649;
t639 = qJ(5) * t497;
t658 = -t501 * t653 - t639;
t642 = pkin(8) * qJD(4);
t580 = t476 * t642;
t657 = -t580 + t675;
t376 = t426 * t496 + t427 * t500;
t448 = t496 * t497 + t500 * t501;
t419 = t448 * t502;
t539 = t428 * t496 + t429 * t500;
t589 = -pkin(9) * t445 + t588;
t353 = -t476 * t653 + t589;
t592 = qJD(6) * t500;
t577 = t496 * t339 + t500 * t341 + t353 * t592;
t655 = -g(1) * t376 + g(2) * t539 - g(3) * t419 - t359 * t538 + t577;
t654 = t445 ^ 2;
t652 = pkin(8) - pkin(9);
t651 = pkin(4) * t441;
t650 = pkin(4) * t497;
t643 = g(3) * t502;
t641 = pkin(1) * qJDD(1);
t640 = qJ(5) * t443;
t378 = t497 * t425 + t501 * t432;
t465 = t476 * qJ(5);
t367 = t465 + t378;
t637 = t367 * t476;
t636 = t378 * t476;
t635 = t380 * t497;
t632 = t441 * t501;
t631 = t443 * t445;
t630 = t443 * t476;
t629 = t445 * t476;
t628 = t445 * t501;
t363 = pkin(9) * t443 + t378;
t356 = t363 + t465;
t627 = t496 * t356;
t626 = t497 * t498;
t624 = t498 * t501;
t622 = t498 * t505;
t395 = t663 * t448;
t522 = t448 * qJD(1);
t617 = t498 * t522 + t395;
t593 = qJD(6) * t496;
t396 = t496 * t597 + t497 * t592 - t500 * t599 - t501 * t593;
t537 = t496 * t501 - t497 * t500;
t416 = t537 * t498;
t616 = qJD(1) * t416 + t396;
t615 = -t476 * t528 + t666;
t547 = -t638 + t650;
t614 = t476 * t547 - t666;
t450 = t557 * qJD(1);
t613 = t497 * t450 + t471 * t619;
t612 = g(2) * t502 * t618 + g(3) * t624;
t611 = t497 * t455 + t501 * t622;
t610 = t503 * pkin(1) + t499 * qJ(2);
t495 = t502 ^ 2;
t609 = t498 ^ 2 - t495;
t506 = qJD(3) ^ 2;
t608 = -t506 - t507;
t604 = qJD(3) * t445;
t598 = qJD(4) * t498;
t594 = qJD(5) * t501;
t585 = qJDD(1) * qJ(2);
t582 = qJDD(3) * t498;
t581 = t502 * t645;
t579 = 0.2e1 * qJD(1) * qJD(2);
t462 = t652 * t501;
t576 = -t500 * t380 + t496 * t381 + t443 * t592;
t570 = t505 * t601;
t574 = t497 * t442 + t455 * t597 + t501 * t570;
t384 = qJ(5) * t605 + t613;
t402 = t498 * qJ(5) + t611;
t568 = t505 * t598;
t564 = -t380 * t496 - t500 * t381;
t438 = t471 * t625;
t562 = -t450 * t501 + t438;
t468 = t497 * t622;
t561 = t455 * t501 - t468;
t560 = qJDD(2) - t641;
t559 = qJD(1) + t598;
t555 = g(1) * t428 + g(2) * t426;
t554 = -g(1) * t429 - g(2) * t427;
t553 = g(1) * t503 + g(2) * t499;
t453 = t502 * t466;
t551 = qJDD(3) * pkin(3) - t471 * t602 + t453;
t521 = pkin(9) * t624 - t502 * t653;
t550 = qJD(1) * t521 - t462 * t663 + t562;
t461 = t652 * t497;
t549 = pkin(9) * t497 * t606 + qJD(6) * t461 - t652 * t599 - t384;
t548 = pkin(4) * t501 + t639;
t545 = qJ(5) * t500 - t496 * t653;
t544 = -qJ(5) * t496 - t500 * t653;
t344 = t496 * t353 + t500 * t356;
t366 = -pkin(4) * t476 + t588;
t543 = t366 * t501 - t367 * t497;
t542 = t366 * t497 + t367 * t501;
t374 = t468 + (-pkin(9) * t502 - t455) * t501 - t653 * t498;
t386 = pkin(9) * t625 + t402;
t541 = t374 * t500 - t386 * t496;
t540 = t374 * t496 + t386 * t500;
t533 = pkin(3) + t548;
t532 = -t455 * t599 - t497 * t570 + (t442 - t568) * t501;
t527 = t441 * t497 + t476 * t597;
t526 = -t476 * t599 + t632;
t525 = t587 * t537;
t524 = 0.2e1 * qJ(2) * t586 + qJDD(3) * t505;
t348 = -t445 * t593 + t576;
t517 = t433 * t476 - t649;
t516 = t498 * t665 - t643;
t515 = g(1) * t426 - g(2) * t428 + g(3) * t625 + t558;
t355 = qJ(5) * t601 + t498 * qJD(5) - t497 * t568 + t574;
t514 = -qJ(5) * t380 + qJD(5) * t445 + t551;
t513 = -t553 + t579 + 0.2e1 * t585;
t349 = qJD(6) * t392 + t564;
t346 = t536 - t651;
t512 = qJD(4) * t543 + t345 * t501 + t346 * t497;
t511 = -t505 * t506 + t513;
t510 = t369 * t445 + qJDD(5) - t515;
t509 = g(1) * t427 - g(2) * t429 + g(3) * t619 + t377 * t476 - t523;
t491 = t503 * qJ(2);
t486 = qJDD(3) * t502;
t475 = qJ(5) * t619;
t437 = pkin(3) - t658;
t417 = t448 * t498;
t407 = -t475 + (-t505 + t650) * t502;
t403 = -pkin(4) * t498 - t561;
t401 = t475 + (t505 - t578) * t502;
t397 = pkin(4) * t445 + t640;
t385 = -pkin(4) * t605 + t562;
t370 = -t445 * t653 - t640;
t365 = (qJD(4) * t548 - t594) * t502 + (t505 - t547) * t602;
t361 = -t380 + t630;
t360 = -pkin(4) * t601 - t532;
t358 = t502 * t537 * t663 - t448 * t602;
t357 = qJD(3) * t416 + t395 * t502;
t354 = (qJD(4) * t658 + t594) * t502 + (-t505 + t528) * t602;
t351 = (t501 * t596 - t572) * pkin(9) + t355;
t350 = pkin(9) * t569 + qJD(3) * t521 - t532;
t347 = pkin(4) * t381 - t514;
t343 = t353 * t500 - t627;
t342 = -t381 * t653 + t514;
t1 = [(t348 * t419 + t358 * t392) * MDP(25) + (-t498 * t506 + t486) * MDP(9) + (qJDD(2) + t665 - 0.2e1 * t641) * MDP(4) - t665 * MDP(2) + (t436 * t498 - t587 * t601) * MDP(29) + ((t350 * t500 - t351 * t496) * t587 - t541 * t436 - t565 * t498 - t343 * t601 + t354 * t538 + t401 * t349 + t342 * t418 - t359 * t357 - g(1) * t539 - g(2) * t376 + (t344 * t498 - t540 * t587) * qJD(6)) * MDP(30) + (t349 * t498 + t357 * t587 + t418 * t436 + t538 * t601) * MDP(28) + (-t348 * t498 + t358 * t587 - t392 * t601 - t419 * t436) * MDP(27) + (-(qJD(6) * t541 + t350 * t496 + t351 * t500) * t587 + t540 * t436 + (-t356 * t593 + t577) * t498 + t344 * t601 + t354 * t392 + t401 * t348 + t342 * t419 + t359 * t358 + g(1) * t563 - g(2) * t375) * MDP(31) + (t345 * t402 + t367 * t355 + t347 * t407 + t369 * t365 + t346 * t403 + t366 * t360 - g(1) * (pkin(3) * t623 + pkin(4) * t429 + qJ(5) * t428 - t503 * t648 + t491) - g(2) * (pkin(4) * t427 + pkin(7) * t503 + qJ(5) * t426 + t610) + (-g(1) * t505 - g(2) * t556) * t499) * MDP(24) + (-t380 * t619 - t445 * t519) * MDP(14) + (-t502 * t506 - t582) * MDP(10) + qJDD(1) * MDP(1) + ((t443 * t501 + t445 * t497) * t602 + (t635 - t381 * t501 + (t443 * t497 - t628) * qJD(4)) * t502) * MDP(15) + (-t574 * t476 - t611 * t441 + ((t476 * t505 + t432) * t599 + (-t433 * t501 + t445 * t505) * qJD(3) + t575) * t498 + (-qJD(3) * t378 + t380 * t505 - t433 * t599 - t501 * t551) * t502 + t555) * MDP(20) + (t532 * t476 + t561 * t441 + ((-t433 * t497 + t443 * t505) * qJD(3) + t558) * t498 + (qJD(3) * t377 - t381 * t505 + t433 * t597 - t497 * t551) * t502 + t554) * MDP(19) + 0.2e1 * (-t498 * t583 + t586 * t609) * MDP(8) + (-t560 * pkin(1) - g(1) * (-pkin(1) * t499 + t491) - g(2) * t610 + (t579 + t585) * qJ(2)) * MDP(6) + t513 * MDP(5) + (-t360 * t476 + t365 * t443 + t381 * t407 - t403 * t441 + (-t369 * t603 - t346) * t498 + (-qJD(3) * t366 + t347 * t497 + t369 * t597) * t502 + t554) * MDP(21) + ((t476 * t603 - t381) * t498 + (-qJD(3) * t443 - t527) * t502) * MDP(17) + ((-t476 * t590 - t380) * t498 + (t526 + t604) * t502) * MDP(16) + (t498 * t511 + t502 * t524) * MDP(12) + (-t498 * t524 + t502 * t511) * MDP(13) + (-t348 * t418 - t349 * t419 + t357 * t392 - t358 * t538) * MDP(26) + t553 * MDP(3) + (t355 * t476 - t365 * t445 + t380 * t407 + t402 * t441 + (t369 * t590 + t345) * t498 + (qJD(3) * t367 - t347 * t501 + t369 * t599) * t502 - t555) * MDP(23) + (qJDD(1) * t495 - 0.2e1 * t498 * t567) * MDP(7) + (t441 * t498 + t476 * t601) * MDP(18) + (-t355 * t443 + t360 * t445 - t380 * t403 - t381 * t402 - t543 * t602 + (-qJD(4) * t542 - t345 * t497 + t346 * t501 + t553) * t502) * MDP(22); qJDD(1) * MDP(4) - t507 * MDP(5) + (t560 + t520) * MDP(6) + (t498 * t608 + t486) * MDP(12) + (t502 * t608 - t582) * MDP(13) + ((-t381 * t498 - t443 * t601 + t445 * t559) * t501 + (-t380 * t498 + t443 * t559 + t445 * t601) * t497) * MDP(22) + (t543 * qJD(1) + (qJD(3) * t542 - t347) * t502 + (qJD(3) * t369 + t512) * t498 + t665) * MDP(24) + (t416 * t436 + t587 * t522 + (-qJD(3) * t525 + t349) * t502 + (-qJD(3) * t538 + t395 * t587) * t498) * MDP(30) + (t417 * t436 - qJD(1) * t525 + (-qJD(3) * t448 * t587 + t348) * t502 + (-qJD(3) * t392 - t396 * t587) * t498) * MDP(31) + (MDP(19) + MDP(21)) * (-t441 * t626 - t381 * t502 + t443 * t602 + (-t497 * t601 - t501 * t559) * t476) + (MDP(20) - MDP(23)) * (t476 * (t497 * t559 - t502 * t590) + (t604 - t632) * t498 + t380 * t502); (pkin(3) * t380 + t613 * t476 - t445 * t454 + t517 * t501 + (-t551 - t657) * t497) * MDP(20) + MDP(9) * t583 + ((t461 * t496 + t462 * t500) * t436 + t437 * t348 - t342 * t537 - g(3) * t416 - (-t496 * t550 + t500 * t549) * t587 + t615 * t392 + t617 * t359 + (-t344 * qJD(1) - t537 * t665) * t502) * MDP(31) + (-t366 * t385 - t367 * t384 + t614 * t369 + (t512 + t516) * pkin(8) + (-t347 + t675) * t533) * MDP(24) + (-(t461 * t500 - t462 * t496) * t436 + t437 * t349 + g(3) * t417 - (t496 * t549 + t500 * t550) * t587 + t615 * t538 + t616 * t359 + (t342 + t662) * t448) * MDP(30) + (-t381 * t533 + t385 * t476 + t614 * t443 + (-t347 - t580 - t581) * t501 + t661 * t497 + t612) * MDP(21) + (-t380 * t533 - t384 * t476 - t614 * t445 - t661 * t501 + (-t347 + t657) * t497) * MDP(23) - MDP(10) * t584 + qJDD(3) * MDP(11) + (t476 * t628 - t635) * MDP(14) + (t643 + (-t466 - t520) * t498) * MDP(13) + (t502 * t520 + t453 + t644) * MDP(12) + (t384 * t443 - t385 * t445 + (t345 + t476 * t366 + (-t381 + t600) * pkin(8)) * t501 + (t346 - t637 + (qJD(4) * t443 - t380) * pkin(8)) * t497 + t516) * MDP(22) + ((-t380 - t630) * t501 + (-t381 - t629) * t497) * MDP(15) + ((-t445 * t502 + t476 * t624) * qJD(1) + t527) * MDP(16) + ((t443 * t502 - t476 * t626) * qJD(1) + t526) * MDP(17) + (t436 * t448 - t587 * t616) * MDP(28) + (t436 * t537 + t587 * t617) * MDP(27) + (-t443 * t454 - pkin(3) * t381 + t438 * t476 + (-t581 + t551 + (-t450 - t642) * t476) * t501 + t517 * t497 + t612) * MDP(19) + (-t348 * t537 + t392 * t617) * MDP(25) + (-t348 * t448 + t349 * t537 - t392 * t616 - t538 * t617) * MDP(26) + (-MDP(18) * t476 - MDP(19) * t377 + MDP(20) * t378 + MDP(21) * t366 - MDP(23) * t367 + MDP(27) * t392 - MDP(28) * t538 + MDP(29) * t587 + MDP(30) * t343) * t605 + (MDP(7) * t498 * t502 - MDP(8) * t609) * t507; MDP(14) * t631 + (-t443 ^ 2 + t654) * MDP(15) + t361 * MDP(16) + (t629 - t381) * MDP(17) + t441 * MDP(18) + (-t433 * t445 + t515 + t636) * MDP(19) + (t433 * t443 + t509) * MDP(20) + (-t397 * t443 - t510 + t636 + 0.2e1 * t651) * MDP(21) + (pkin(4) * t380 - qJ(5) * t381 + (t367 - t378) * t445 + (t366 - t588) * t443) * MDP(22) + (-t369 * t443 + t397 * t445 + 0.2e1 * t431 + 0.2e1 * t463 - t509) * MDP(23) + (t345 * qJ(5) - t346 * pkin(4) - t369 * t397 - t366 * t378 - g(1) * (-pkin(4) * t426 + qJ(5) * t427) - g(2) * (pkin(4) * t428 - qJ(5) * t429) - g(3) * (-pkin(4) * t625 + t475) + t588 * t367) * MDP(24) + (-t348 - t667) * MDP(27) + (t349 - t673) * MDP(28) + (-t544 * t436 - t370 * t538 - (t500 * t363 + t496 * t589) * t587 + (-t545 * t587 + t344) * qJD(6) + t670) * MDP(30) + (t545 * t436 - t370 * t392 - (-t496 * t363 + t500 * t589) * t587 + (-t544 * t587 - t627) * qJD(6) + t655) * MDP(31) - t674; (-t441 + t631) * MDP(21) + t361 * MDP(22) + (-t476 ^ 2 - t654) * MDP(23) + (t510 - t637 - t651) * MDP(24) + (-t436 * t500 - t445 * t538) * MDP(30) + (-t392 * t445 + t436 * t496) * MDP(31) - (MDP(30) * t496 + MDP(31) * t500) * t587 ^ 2; (t576 + t667) * MDP(27) + (-t564 + t673) * MDP(28) + (t344 * t587 - t670) * MDP(30) + (t343 * t587 - t655) * MDP(31) + ((-MDP(28) * t445 - MDP(30) * t356) * t500 + (-MDP(27) * t445 - MDP(28) * t443 - MDP(30) * t353 + MDP(31) * t356) * t496) * qJD(6) + t674;];
tau  = t1;
