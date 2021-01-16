% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:24
% EndTime: 2021-01-16 02:06:50
% DurationCPUTime: 12.89s
% Computational Cost: add. (5156->567), mult. (12292->777), div. (0->0), fcn. (9902->18), ass. (0->246)
t548 = cos(qJ(3));
t667 = cos(pkin(11));
t597 = t667 * t548;
t516 = qJD(2) * t597;
t538 = sin(pkin(11));
t545 = sin(qJ(3));
t629 = qJD(2) * t545;
t484 = t538 * t629 - t516;
t477 = qJD(6) + t484;
t598 = t667 * t545;
t498 = t538 * t548 + t598;
t487 = t498 * qJD(2);
t537 = sin(pkin(12));
t541 = cos(pkin(12));
t465 = t541 * qJD(3) - t487 * t537;
t547 = cos(qJ(6));
t464 = qJD(3) * t537 + t487 * t541;
t544 = sin(qJ(6));
t662 = t464 * t544;
t685 = t465 * t547 - t662;
t688 = t685 * t477;
t546 = sin(qJ(2));
t540 = sin(pkin(6));
t631 = qJD(1) * t540;
t611 = t546 * t631;
t669 = qJD(3) * pkin(3);
t687 = t545 * t669 - t611;
t486 = t498 * qJD(3);
t652 = t538 * t545;
t568 = t597 - t652;
t489 = t568 * qJD(3);
t686 = pkin(4) * t486 - qJ(5) * t489 - qJD(5) * t498 + t687;
t543 = qJ(4) + pkin(8);
t603 = qJD(3) * t543;
t478 = qJD(4) * t548 - t545 * t603;
t564 = -qJD(4) * t545 - t548 * t603;
t549 = cos(qJ(2));
t610 = t549 * t631;
t637 = t478 * t667 + t538 * t564 - t568 * t610;
t577 = -t464 * t547 - t465 * t544;
t683 = t477 * t577;
t499 = t537 * t547 + t541 * t544;
t491 = t499 * qJD(6);
t634 = t499 * t484 + t491;
t542 = cos(pkin(6));
t648 = t540 * t546;
t492 = t542 * t548 - t545 * t648;
t534 = qJ(3) + pkin(11);
t530 = sin(t534);
t668 = cos(pkin(10));
t600 = t668 * t549;
t539 = sin(pkin(10));
t650 = t539 * t546;
t481 = -t542 * t600 + t650;
t601 = t668 * t546;
t649 = t539 * t549;
t483 = t542 * t649 + t601;
t589 = g(1) * t483 + g(2) * t481;
t646 = t540 * t549;
t681 = g(3) * t646 - t589;
t682 = t681 * t530;
t618 = t542 * qJDD(1);
t515 = t548 * t618;
t622 = qJD(1) * qJD(2);
t470 = qJDD(2) * pkin(8) + (qJDD(1) * t546 + t549 * t622) * t540;
t630 = qJD(1) * t542;
t554 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t630 + t470;
t593 = qJD(2) * t543 + t611;
t575 = t593 * qJD(3);
t386 = qJDD(3) * pkin(3) - t545 * t554 - t548 * t575 + t515;
t387 = (-t575 + t618) * t545 + t554 * t548;
t356 = t667 * t386 - t538 * t387;
t355 = -qJDD(3) * pkin(4) + qJDD(5) - t356;
t482 = t542 * t601 + t649;
t532 = cos(t534);
t602 = t540 * t668;
t443 = t482 * t530 + t532 * t602;
t480 = t542 * t650 - t600;
t651 = t539 * t540;
t445 = t480 * t530 + t532 * t651;
t472 = t530 * t648 - t542 * t532;
t566 = g(1) * t445 - g(2) * t443 - g(3) * t472;
t561 = t355 + t566;
t641 = -t637 * t537 + t686 * t541;
t640 = t686 * t537 + t637 * t541;
t638 = t478 * t538 - t498 * t610 - t667 * t564;
t456 = t545 * t630 + t548 * t593;
t447 = t538 * t456;
t455 = -t545 * t593 + t548 * t630;
t399 = t455 * t667 - t447;
t616 = pkin(3) * t629;
t426 = pkin(4) * t487 + qJ(5) * t484 + t616;
t371 = t541 * t399 + t537 * t426;
t680 = qJD(5) * t541 - t371;
t370 = -t399 * t537 + t541 * t426;
t679 = -qJD(5) * t537 - t370;
t620 = qJDD(2) * t545;
t584 = -qJDD(2) * t597 + t538 * t620;
t439 = qJD(2) * t486 + t584;
t435 = qJDD(6) + t439;
t497 = t537 * t544 - t547 * t541;
t635 = t477 * t497;
t678 = -t435 * t499 + t477 * t635;
t550 = qJD(3) ^ 2;
t604 = qJDD(1) * t646;
t607 = t546 * t622;
t585 = t540 * t607 - t604;
t677 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t550 + t540 * (-g(3) * t549 + t607) - t585 + t589;
t479 = t484 ^ 2;
t676 = pkin(3) * t538;
t675 = pkin(3) * t545;
t674 = pkin(9) * t541;
t672 = g(3) * t540;
t521 = qJ(5) + t676;
t671 = pkin(9) + t521;
t670 = qJD(2) * pkin(2);
t665 = t685 * t487;
t664 = t577 * t487;
t661 = t484 * t537;
t660 = t489 * t537;
t659 = t498 * t537;
t658 = t498 * t541;
t657 = t521 * t537;
t656 = t521 * t541;
t533 = pkin(12) + qJ(6);
t529 = sin(t533);
t655 = t529 * t532;
t531 = cos(t533);
t654 = t531 * t532;
t653 = t532 * t549;
t647 = t540 * t548;
t644 = qJDD(1) - g(3);
t357 = t538 * t386 + t667 * t387;
t354 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t357;
t526 = pkin(3) * t548 + pkin(2);
t621 = qJD(2) * qJD(3);
t606 = t545 * t621;
t438 = pkin(3) * t606 - qJDD(2) * t526 + qJDD(4) + t585;
t556 = qJDD(2) * t498 - t538 * t606;
t440 = qJD(3) * t516 + t556;
t369 = pkin(4) * t439 - qJ(5) * t440 - qJD(5) * t487 + t438;
t350 = t541 * t354 + t537 * t369;
t643 = pkin(5) * t486 - t489 * t674 + t641;
t642 = pkin(9) * t660 - t640;
t451 = t455 + t669;
t599 = t667 * t456;
t395 = t538 * t451 + t599;
t390 = qJD(3) * qJ(5) + t395;
t476 = -qJD(2) * t526 + qJD(4) - t610;
t409 = pkin(4) * t484 - qJ(5) * t487 + t476;
t365 = t541 * t390 + t537 * t409;
t639 = pkin(5) * t660 + t638;
t436 = -pkin(4) * t568 - qJ(5) * t498 - t526;
t508 = t543 * t548;
t461 = t508 * t667 - t543 * t652;
t392 = t537 * t436 + t541 * t461;
t636 = -t480 * t543 - t483 * t526;
t633 = t526 * t646 + t543 * t648;
t535 = t545 ^ 2;
t632 = -t548 ^ 2 + t535;
t628 = qJD(2) * t546;
t625 = qJD(6) * t544;
t624 = qJD(6) * t547;
t394 = t451 * t667 - t447;
t389 = -qJD(3) * pkin(4) + qJD(5) - t394;
t623 = -qJD(5) + t389;
t619 = qJDD(2) * t548;
t617 = g(3) * t648;
t422 = -t541 * qJDD(3) + t440 * t537;
t423 = qJDD(3) * t537 + t440 * t541;
t613 = -t544 * t422 + t547 * t423 + t465 * t624;
t612 = t667 * pkin(3);
t609 = t540 * t628;
t608 = qJD(2) * t646;
t605 = t548 * t621;
t349 = -t354 * t537 + t541 * t369;
t345 = pkin(5) * t439 - pkin(9) * t423 + t349;
t346 = -pkin(9) * t422 + t350;
t596 = t547 * t345 - t544 * t346;
t364 = -t390 * t537 + t541 * t409;
t595 = t547 * t422 + t544 * t423;
t391 = t541 * t436 - t461 * t537;
t397 = t455 * t538 + t599;
t594 = -t481 * t526 + t482 * t543;
t442 = t480 * t532 - t530 * t651;
t460 = t508 * t538 + t543 * t598;
t592 = t545 * t608;
t591 = t539 * pkin(3) * t647 + t480 * t675;
t525 = -t612 - pkin(4);
t590 = g(1) * t480 - g(2) * t482;
t587 = -t497 * t435 - t477 * t634;
t586 = pkin(4) * t532 + qJ(5) * t530;
t583 = t544 * t345 + t547 * t346;
t582 = -t349 * t541 - t350 * t537;
t352 = pkin(5) * t484 - pkin(9) * t464 + t364;
t358 = pkin(9) * t465 + t365;
t347 = t352 * t547 - t358 * t544;
t348 = t352 * t544 + t358 * t547;
t374 = -pkin(5) * t568 - pkin(9) * t658 + t391;
t376 = -pkin(9) * t659 + t392;
t581 = t374 * t547 - t376 * t544;
t580 = t374 * t544 + t376 * t547;
t493 = t542 * t545 + t546 * t647;
t432 = t538 * t492 + t493 * t667;
t410 = -t432 * t537 - t541 * t646;
t411 = t432 * t541 - t537 * t646;
t579 = t410 * t547 - t411 * t544;
t578 = t410 * t544 + t411 * t547;
t576 = t492 * pkin(3);
t551 = qJD(2) ^ 2;
t574 = qJDD(2) * t549 - t546 * t551;
t572 = -g(1) * t539 + g(2) * t668;
t495 = t671 * t541;
t570 = pkin(5) * t487 + qJD(6) * t495 + t484 * t674 - t679;
t494 = t671 * t537;
t569 = pkin(9) * t661 + qJD(6) * t494 - t680;
t359 = -t464 * t625 + t613;
t444 = t482 * t532 - t530 * t602;
t473 = t530 * t542 + t532 * t648;
t567 = g(1) * t442 - g(2) * t444 - g(3) * t473;
t565 = t574 * t540;
t563 = t493 * qJD(3);
t505 = -t610 - t670;
t560 = -qJD(2) * t505 - t470 - t590;
t559 = (-t482 * t545 - t548 * t602) * pkin(3);
t558 = t681 * t532;
t557 = t355 * t498 + t389 * t489 + t590;
t360 = -qJD(6) * t577 + t595;
t553 = -pkin(8) * qJDD(3) + (t505 + t610 - t670) * qJD(3);
t552 = -t563 - t592;
t506 = -t541 * pkin(5) + t525;
t454 = qJD(3) * t492 + t548 * t608;
t431 = -t492 * t667 + t493 * t538;
t430 = t497 * t498;
t429 = t499 * t498;
t421 = pkin(5) * t659 + t460;
t398 = t454 * t667 + t538 * t552;
t396 = t454 * t538 - t552 * t667;
t383 = t398 * t541 + t537 * t609;
t382 = -t398 * t537 + t541 * t609;
t379 = -pkin(5) * t661 + t397;
t378 = t489 * t499 + t624 * t658 - t625 * t659;
t377 = -t489 * t497 - t491 * t498;
t375 = -pkin(5) * t465 + t389;
t351 = pkin(5) * t422 + t355;
t1 = [t644 * MDP(1) + MDP(3) * t565 + (-qJDD(2) * t546 - t549 * t551) * t540 * MDP(4) + (t492 * qJDD(3) + t548 * t565 + (-t563 - 0.2e1 * t592) * qJD(3)) * MDP(10) + (-qJD(3) * t454 - qJDD(3) * t493 + (-t545 * t574 - t549 * t605) * t540) * MDP(11) + (-qJD(3) * t396 - qJDD(3) * t431 + (-t439 * t549 + t484 * t628) * t540) * MDP(12) + (-qJD(3) * t398 - qJDD(3) * t432 + (-t440 * t549 + t487 * t628) * t540) * MDP(13) + (t396 * t487 - t398 * t484 + t431 * t440 - t432 * t439) * MDP(14) + (-t356 * t431 + t357 * t432 - t394 * t396 + t395 * t398 - g(3) + (-t438 * t549 + t476 * t628) * t540) * MDP(15) + (t382 * t484 - t396 * t465 + t410 * t439 + t422 * t431) * MDP(16) + (-t383 * t484 + t396 * t464 - t411 * t439 + t423 * t431) * MDP(17) + (-t382 * t464 + t383 * t465 - t410 * t423 - t411 * t422) * MDP(18) + (t349 * t410 + t350 * t411 + t355 * t431 + t364 * t382 + t365 * t383 + t389 * t396 - g(3)) * MDP(19) + ((-qJD(6) * t578 + t382 * t547 - t383 * t544) * t477 + t579 * t435 - t396 * t685 + t431 * t360) * MDP(25) + (-(qJD(6) * t579 + t382 * t544 + t383 * t547) * t477 - t578 * t435 - t396 * t577 + t431 * t359) * MDP(26); qJDD(2) * MDP(2) + (-t681 + t604) * MDP(3) + (-t644 * t648 - t590) * MDP(4) + (qJDD(2) * t535 + 0.2e1 * t545 * t605) * MDP(5) + 0.2e1 * (t545 * t619 - t621 * t632) * MDP(6) + (qJDD(3) * t545 + t548 * t550) * MDP(7) + (qJDD(3) * t548 - t545 * t550) * MDP(8) + (t553 * t545 + t548 * t677) * MDP(10) + (-t545 * t677 + t553 * t548) * MDP(11) + (-t484 * t611 - qJDD(3) * t460 - t438 * t568 - t439 * t526 + t476 * t486 - t558 + (t484 * t675 - t638) * qJD(3)) * MDP(12) + (-t487 * t611 - qJDD(3) * t461 + t438 * t498 - t440 * t526 + t476 * t489 + t682 + (t487 * t675 - t637) * qJD(3)) * MDP(13) + (-t356 * t498 + t357 * t568 - t394 * t489 - t395 * t486 - t439 * t461 + t440 * t460 - t484 * t637 + t487 * t638 + t590 - t617) * MDP(14) + (-g(1) * t636 - g(2) * t594 - g(3) * t633 - t356 * t460 + t357 * t461 - t638 * t394 + t637 * t395 - t438 * t526 + t476 * t687) * MDP(15) + (-t349 * t568 + t364 * t486 + t391 * t439 + t460 * t422 - t541 * t558 + (t557 - t617) * t537 + t641 * t484 - t638 * t465) * MDP(16) + (t350 * t568 - t365 * t486 - t392 * t439 + t460 * t423 - t589 * t537 * t532 + t557 * t541 - (-t537 * t653 + t541 * t546) * t672 - t640 * t484 + t638 * t464) * MDP(17) + (-t391 * t423 - t392 * t422 + t582 * t498 + (-t364 * t541 - t365 * t537) * t489 - t641 * t464 + t640 * t465 - t682) * MDP(18) + (t350 * t392 + t349 * t391 + t355 * t460 - g(1) * (-t483 * t586 + t636) - g(2) * (-t481 * t586 + t594) - g(3) * (t586 * t646 + t633) + t638 * t389 + t640 * t365 + t641 * t364) * MDP(19) + (-t359 * t430 - t377 * t577) * MDP(20) + (-t359 * t429 + t360 * t430 + t377 * t685 + t378 * t577) * MDP(21) + (-t359 * t568 + t377 * t477 - t430 * t435 - t486 * t577) * MDP(22) + (t360 * t568 - t378 * t477 - t429 * t435 + t486 * t685) * MDP(23) + (-t435 * t568 + t477 * t486) * MDP(24) + (t581 * t435 - t596 * t568 + t347 * t486 + t421 * t360 + t351 * t429 + t375 * t378 - g(1) * (-t480 * t529 - t483 * t654) - g(2) * (-t481 * t654 + t482 * t529) - (t529 * t546 + t531 * t653) * t672 + (t544 * t642 + t547 * t643) * t477 - t639 * t685 + (t348 * t568 - t477 * t580) * qJD(6)) * MDP(25) + (-t580 * t435 + t583 * t568 - t348 * t486 + t421 * t359 - t351 * t430 + t375 * t377 - g(1) * (-t480 * t531 + t483 * t655) - g(2) * (t481 * t655 + t482 * t531) - (-t529 * t653 + t531 * t546) * t672 + (-t544 * t643 + t547 * t642) * t477 - t639 * t577 + (t347 * t568 - t477 * t581) * qJD(6)) * MDP(26); MDP(7) * t620 + MDP(8) * t619 + qJDD(3) * MDP(9) + (-g(3) * t492 + t545 * t560 + t572 * t647 + t515) * MDP(10) + (g(3) * t493 + (-t540 * t572 - t618) * t545 + t560 * t548) * MDP(11) + (t397 * qJD(3) - t476 * t487 + (qJDD(3) * t667 - t484 * t629) * pkin(3) - t566 + t356) * MDP(12) + (qJD(3) * t399 + t476 * t484 + (-qJDD(3) * t538 - t487 * t629) * pkin(3) - t567 - t357) * MDP(13) + ((t395 - t397) * t487 + (-t394 + t399) * t484 + (-t439 * t538 - t440 * t667) * pkin(3)) * MDP(14) + (-g(1) * t591 - g(2) * t559 - g(3) * t576 + t356 * t612 + t357 * t676 + t394 * t397 - t395 * t399 - t476 * t616) * MDP(15) + (-t439 * t657 - t364 * t487 + t397 * t465 + t422 * t525 + (t537 * t623 - t370) * t484 - t561 * t541) * MDP(16) + (-t439 * t656 + t365 * t487 - t397 * t464 + t423 * t525 + (t541 * t623 + t371) * t484 + t561 * t537) * MDP(17) + (t370 * t464 - t371 * t465 + (qJD(5) * t465 - t364 * t484 - t422 * t521 + t350) * t541 + (qJD(5) * t464 - t365 * t484 + t423 * t521 - t349) * t537 + t567) * MDP(18) + (t350 * t656 - t349 * t657 + t355 * t525 - t389 * t397 - g(1) * (pkin(4) * t445 - qJ(5) * t442 + t591) - g(2) * (-t443 * pkin(4) + t444 * qJ(5) + t559) - g(3) * (-pkin(4) * t472 + qJ(5) * t473 + t576) + t680 * t365 + t679 * t364) * MDP(19) + (t359 * t499 + t577 * t635) * MDP(20) + (-t359 * t497 - t360 * t499 + t577 * t634 - t635 * t685) * MDP(21) + (t664 - t678) * MDP(22) + (t587 - t665) * MDP(23) - t477 * t487 * MDP(24) + ((-t494 * t547 - t495 * t544) * t435 + t506 * t360 + t351 * t497 - t347 * t487 + t379 * t685 + (t544 * t569 - t547 * t570) * t477 + t634 * t375 - t566 * t531) * MDP(25) + (-(-t494 * t544 + t495 * t547) * t435 + t506 * t359 + t351 * t499 + t348 * t487 + t379 * t577 + (t544 * t570 + t547 * t569) * t477 - t635 * t375 + t566 * t529) * MDP(26) + (-MDP(5) * t545 * t548 + MDP(6) * t632) * t551; (0.2e1 * t487 * qJD(3) + t584) * MDP(12) + ((t516 - t484) * qJD(3) + t556) * MDP(13) + (-t487 ^ 2 - t479) * MDP(14) + (t394 * t487 + t395 * t484 + t438 + t681) * MDP(15) + (t439 * t541 + t465 * t487 - t479 * t537) * MDP(16) + (-t439 * t537 - t464 * t487 - t479 * t541) * MDP(17) + (-t422 * t537 - t423 * t541 + (t464 * t537 + t465 * t541) * t484) * MDP(18) + (-t389 * t487 + (-t364 * t537 + t365 * t541) * t484 + t681 - t582) * MDP(19) + (t587 + t665) * MDP(25) + (t664 + t678) * MDP(26); (t464 * t484 + t422) * MDP(16) + (t465 * t484 + t423) * MDP(17) + (-t464 ^ 2 - t465 ^ 2) * MDP(18) + (t364 * t464 - t365 * t465 + t561) * MDP(19) + (t360 - t683) * MDP(25) + (t359 + t688) * MDP(26); t577 * t685 * MDP(20) + (t577 ^ 2 - t685 ^ 2) * MDP(21) + (t613 - t688) * MDP(22) + (-t595 - t683) * MDP(23) + t435 * MDP(24) + (t348 * t477 + t375 * t577 - g(1) * (t442 * t529 + t483 * t531) - g(2) * (-t444 * t529 + t481 * t531) - g(3) * (-t473 * t529 - t531 * t646) + t596) * MDP(25) + (t347 * t477 - t375 * t685 - g(1) * (t442 * t531 - t483 * t529) - g(2) * (-t444 * t531 - t481 * t529) - g(3) * (-t473 * t531 + t529 * t646) - t583) * MDP(26) + (-MDP(22) * t662 + MDP(23) * t577 - MDP(25) * t348 - MDP(26) * t347) * qJD(6);];
tau = t1;
