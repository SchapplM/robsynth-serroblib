% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:46
% EndTime: 2019-03-09 00:35:04
% DurationCPUTime: 12.55s
% Computational Cost: add. (8498->589), mult. (22674->804), div. (0->0), fcn. (18316->12), ass. (0->236)
t522 = sin(pkin(7));
t528 = sin(qJ(3));
t652 = t522 * t528;
t515 = pkin(9) * t652;
t524 = cos(pkin(7));
t532 = cos(qJ(3));
t533 = cos(qJ(2));
t643 = t532 * t533;
t529 = sin(qJ(2));
t647 = t528 * t529;
t551 = -t524 * t647 + t643;
t523 = sin(pkin(6));
t624 = qJD(1) * t523;
t648 = t524 * t532;
t630 = t551 * t624 - (pkin(2) * t648 - t515) * qJD(3);
t571 = pkin(3) * t528 - pkin(10) * t532;
t550 = t571 * qJD(3);
t601 = t529 * t624;
t691 = (t550 - t601) * t522;
t619 = qJD(2) * t532;
t514 = t522 * t619;
t565 = t514 - qJD(4);
t651 = t522 * t532;
t606 = pkin(9) * t651;
t481 = t606 + (pkin(2) * t528 + pkin(10)) * t524;
t572 = -pkin(3) * t532 - pkin(10) * t528;
t482 = (-pkin(2) + t572) * t522;
t527 = sin(qJ(4));
t531 = cos(qJ(4));
t613 = qJD(4) * t531;
t615 = qJD(4) * t527;
t690 = -t481 * t615 + t482 * t613 + t691 * t527 - t630 * t531;
t645 = t529 * t532;
t646 = t528 * t533;
t553 = t524 * t645 + t646;
t649 = t524 * t528;
t628 = -t553 * t624 + (pkin(2) * t649 + t606) * qJD(3);
t618 = qJD(3) * t528;
t599 = t522 * t618;
t689 = -pkin(11) * t599 - t690;
t493 = -t531 * t524 + t527 * t652;
t617 = qJD(3) * t532;
t598 = t522 * t617;
t455 = -qJD(4) * t493 + t531 * t598;
t494 = t524 * t527 + t531 * t652;
t456 = qJD(4) * t494 + t527 * t598;
t688 = pkin(4) * t456 - pkin(11) * t455 + t628;
t525 = cos(pkin(6));
t623 = qJD(1) * t525;
t602 = t522 * t623;
t505 = t528 * t602;
t622 = qJD(2) * t522;
t498 = pkin(9) * t622 + t601;
t668 = qJD(2) * pkin(2);
t504 = t533 * t624 + t668;
t677 = t532 * t498 + t504 * t649;
t436 = t505 + t677;
t687 = t436 + t565 * (pkin(4) * t527 - pkin(11) * t531);
t686 = pkin(10) * t615;
t620 = qJD(2) * t524;
t587 = qJD(3) + t620;
t600 = t528 * t622;
t537 = -t527 * t600 + t531 * t587;
t607 = qJD(2) * qJD(3);
t593 = t522 * t607;
t574 = t532 * t593;
t558 = t531 * t574;
t535 = t537 * qJD(4) + t558;
t684 = -qJD(5) * t565 + t535;
t683 = t481 * t613 + t482 * t615 - t630 * t527 - t691 * t531;
t682 = t527 * t514 - t615;
t519 = t522 ^ 2;
t681 = (t528 * t532 * MDP(5) - (t528 ^ 2 - t532 ^ 2) * MDP(6)) * t519;
t469 = qJD(5) - t537;
t634 = -pkin(4) * t599 + t683;
t480 = t515 + (-pkin(2) * t532 - pkin(3)) * t524;
t426 = pkin(4) * t493 - pkin(11) * t494 + t480;
t629 = t531 * t481 + t527 * t482;
t428 = -pkin(11) * t651 + t629;
t526 = sin(qJ(5));
t530 = cos(qJ(5));
t611 = qJD(5) * t530;
t612 = qJD(5) * t526;
t679 = -t426 * t611 + t428 * t612 - t688 * t526 + t530 * t689;
t678 = t526 * t426 + t530 * t428;
t489 = t528 * t498;
t654 = t504 * t524;
t435 = t532 * (t602 + t654) - t489;
t484 = t571 * t622;
t631 = t531 * t435 + t527 * t484;
t403 = pkin(11) * t600 + t631;
t512 = -pkin(4) * t531 - pkin(11) * t527 - pkin(3);
t676 = t530 * t403 - t512 * t611 + t526 * t687;
t423 = pkin(10) * t587 + t436;
t513 = t524 * t623;
t442 = t513 + (qJD(2) * t572 - t504) * t522;
t387 = -t527 * t423 + t531 * t442;
t381 = pkin(4) * t565 - t387;
t475 = t527 * t587 + t531 * t600;
t444 = t475 * t526 + t530 * t565;
t446 = t530 * t475 - t526 * t565;
t363 = t444 * pkin(5) - t446 * qJ(6) + t381;
t559 = t527 * t574;
t448 = qJD(4) * t475 + t559;
t669 = pkin(11) * t448;
t675 = t363 * t469 - t669;
t674 = -qJD(5) * t678 + t526 * t689 + t688 * t530;
t673 = t446 ^ 2;
t672 = t469 ^ 2;
t534 = qJD(2) ^ 2;
t671 = pkin(5) * t448;
t670 = pkin(10) * t531;
t666 = qJ(6) * t448;
t575 = t528 * t593;
t540 = t551 * qJD(2);
t577 = t525 * t598;
t405 = (t504 * t648 - t489) * qJD(3) + (t523 * t540 + t577) * qJD(1);
t451 = (t550 + t601) * t622;
t583 = -t527 * t405 - t423 * t613 - t442 * t615 + t531 * t451;
t360 = -pkin(4) * t575 - t583;
t395 = t475 * t612 - t526 * t575 - t530 * t684;
t396 = t475 * t611 + t526 * t684 - t530 * t575;
t349 = pkin(5) * t396 + qJ(6) * t395 - qJD(6) * t446 + t360;
t665 = t349 * t526;
t664 = t349 * t530;
t388 = t531 * t423 + t527 * t442;
t382 = -pkin(11) * t565 + t388;
t422 = -pkin(3) * t587 - t435;
t390 = -pkin(4) * t537 - t475 * pkin(11) + t422;
t362 = t382 * t530 + t390 * t526;
t356 = qJ(6) * t469 + t362;
t663 = t356 * t469;
t662 = t360 * t526;
t661 = t362 * t469;
t660 = t395 * t526;
t659 = t444 * t469;
t658 = t446 * t444;
t657 = t446 * t469;
t656 = t448 * t526;
t655 = t448 * t530;
t653 = t512 * t530;
t650 = t523 * t534;
t644 = t530 * t532;
t642 = qJ(6) * t456 + qJD(6) * t493 - t679;
t641 = -pkin(5) * t456 - t674;
t457 = t494 * t526 + t522 * t644;
t407 = qJD(5) * t457 - t530 * t455 - t526 * t599;
t604 = t526 * t651;
t408 = -qJD(5) * t604 + t455 * t526 + t494 * t611 - t530 * t599;
t458 = t494 * t530 - t604;
t640 = pkin(5) * t408 + qJ(6) * t407 - qJD(6) * t458 + t634;
t610 = qJD(5) * t531;
t614 = qJD(4) * t530;
t639 = qJD(6) * t531 - (-t526 * t610 - t527 * t614) * pkin(10) + t676 + t682 * qJ(6);
t638 = -t512 * t612 + (t403 + t686) * t526 + (-pkin(10) * t610 - t687) * t530 - t682 * pkin(5);
t566 = pkin(5) * t526 - qJ(6) * t530;
t637 = qJD(6) * t526 - t469 * t566 + t388;
t431 = t527 * t435;
t402 = -pkin(4) * t600 - t484 * t531 + t431;
t579 = t531 * t514;
t466 = t526 * t579 - t530 * t600;
t467 = (t526 * t528 + t531 * t644) * t622;
t557 = pkin(10) + t566;
t567 = pkin(5) * t530 + qJ(6) * t526;
t636 = pkin(5) * t466 - qJ(6) * t467 + t402 - (qJD(5) * t567 - qJD(6) * t530) * t527 - t557 * t613;
t434 = pkin(4) * t475 - pkin(11) * t537;
t635 = t530 * t387 + t526 * t434;
t626 = t526 * t512 + t530 * t670;
t621 = qJD(2) * t523;
t616 = qJD(4) * t526;
t609 = t422 * qJD(4);
t361 = -t382 * t526 + t390 * t530;
t608 = qJD(6) - t361;
t605 = pkin(11) * t612;
t603 = t529 * t650;
t597 = t469 * t612;
t595 = qJD(3) * t654;
t594 = qJD(1) * t621;
t592 = -t527 * t481 + t482 * t531;
t591 = t532 * t565;
t590 = t469 * t530;
t589 = qJD(4) * t565;
t585 = t519 * t603;
t545 = -t531 * t405 + t423 * t615 - t442 * t613 - t527 * t451;
t359 = pkin(11) * t575 - t545;
t560 = t524 * t529 * t594;
t373 = t448 * pkin(4) - pkin(11) * t535 + qJD(3) * t505 + t498 * t617 + t528 * t595 + t532 * t560 + t594 * t646;
t348 = -t526 * t359 + t530 * t373 - t382 * t611 - t390 * t612;
t581 = t522 * t529 * t621;
t578 = t525 * t599;
t573 = MDP(16) * t600;
t427 = pkin(4) * t651 - t592;
t569 = t526 * t613 - t466;
t568 = t530 * t613 - t467;
t355 = -pkin(5) * t469 + t608;
t564 = t355 * t530 - t356 * t526;
t552 = t524 * t646 + t645;
t454 = t523 * t552 + t525 * t652;
t492 = -t522 * t523 * t533 + t524 * t525;
t425 = t454 * t531 + t492 * t527;
t554 = t524 * t643 - t647;
t453 = -t523 * t554 - t525 * t651;
t394 = t425 * t530 + t453 * t526;
t393 = t425 * t526 - t453 * t530;
t562 = t426 * t530 - t428 * t526;
t424 = t454 * t527 - t492 * t531;
t468 = -t504 * t522 + t513;
t555 = t468 * t522 - t519 * t668;
t548 = -t469 * t611 - t656;
t547 = t381 * t469 - t669;
t546 = t363 * t446 - t348;
t347 = t530 * t359 + t526 * t373 - t382 * t612 + t390 * t611;
t541 = t553 * qJD(2);
t539 = qJD(3) * t498 + t560;
t536 = -t468 * t622 - t595 + (-qJD(3) * t522 * t525 - t533 * t621) * qJD(1);
t507 = -pkin(4) - t567;
t483 = t557 * t527;
t460 = -t653 + (pkin(10) * t526 + pkin(5)) * t531;
t459 = -qJ(6) * t531 + t626;
t418 = t577 + (qJD(3) * t554 + t540) * t523;
t417 = t578 + (qJD(3) * t552 + t541) * t523;
t406 = t677 * qJD(3) + (t523 * t541 + t578) * qJD(1);
t399 = pkin(5) * t446 + qJ(6) * t444;
t383 = pkin(5) * t457 - qJ(6) * t458 + t427;
t378 = -qJD(4) * t424 + t418 * t531 + t527 * t581;
t377 = qJD(4) * t425 + t418 * t527 - t531 * t581;
t376 = -pkin(5) * t493 - t562;
t375 = qJ(6) * t493 + t678;
t372 = -t395 + t659;
t365 = -pkin(5) * t475 + t387 * t526 - t434 * t530;
t364 = qJ(6) * t475 + t635;
t353 = -qJD(5) * t393 + t378 * t530 + t417 * t526;
t352 = qJD(5) * t394 + t378 * t526 - t417 * t530;
t346 = -t348 - t671;
t345 = qJD(6) * t469 + t347 + t666;
t1 = [-MDP(3) * t603 - t533 * MDP(4) * t650 + (-t417 * t587 + t492 * t575 - t532 * t585) * MDP(10) + (-t418 * t587 + t492 * t574 + t528 * t585) * MDP(11) + (t377 * t565 - t417 * t537 - t424 * t575 + t453 * t448) * MDP(17) + (t378 * t565 + t417 * t475 - t425 * t575 + t453 * t535) * MDP(18) + (t352 * t446 - t353 * t444 - t393 * t395 - t394 * t396) * MDP(27) + (t345 * t394 + t346 * t393 + t349 * t424 + t352 * t355 + t353 * t356 + t363 * t377) * MDP(29) + (MDP(24) + MDP(26)) * (-t352 * t469 + t377 * t444 - t393 * t448 + t424 * t396) + (-MDP(25) + MDP(28)) * (t353 * t469 - t377 * t446 + t394 * t448 + t395 * t424); (t395 * t457 - t396 * t458 + t407 * t444 - t408 * t446) * MDP(20) + (-t395 * t493 - t407 * t469 + t446 * t456 + t448 * t458) * MDP(21) + (-t396 * t493 - t408 * t469 - t444 * t456 - t448 * t457) * MDP(22) + (t448 * t493 + t456 * t469) * MDP(23) + (t348 * t493 + t360 * t457 + t361 * t456 + t381 * t408 + t427 * t396 + t634 * t444 + t562 * t448 + t469 * t674) * MDP(24) + (-t347 * t493 + t360 * t458 - t362 * t456 - t381 * t407 - t427 * t395 + t634 * t446 - t448 * t678 + t469 * t679) * MDP(25) + (-t346 * t493 + t349 * t457 - t355 * t456 + t363 * t408 - t376 * t448 + t383 * t396 + t444 * t640 - t469 * t641) * MDP(26) + (-t345 * t457 + t346 * t458 - t355 * t407 - t356 * t408 - t375 * t396 - t376 * t395 - t444 * t642 + t446 * t641) * MDP(27) + (t345 * t493 - t349 * t458 + t356 * t456 + t363 * t407 + t375 * t448 + t383 * t395 - t446 * t640 + t469 * t642) * MDP(28) + (t345 * t375 + t346 * t376 + t349 * t383 + t355 * t641 + t356 * t642 + t363 * t640) * MDP(29) + ((-qJD(2) * t628 - t406) * t524 + (t528 * t555 - t628) * qJD(3)) * MDP(10) + ((qJD(2) * t630 - t405) * t524 + (t532 * t555 + t630) * qJD(3)) * MDP(11) + (t475 * t455 + t494 * t535) * MDP(12) + (-t494 * t448 + t455 * t537 - t475 * t456 - t493 * t535) * MDP(13) + (-t455 * t565 + t475 * t599 + t494 * t575 - t535 * t651) * MDP(14) + (t456 * t565 + (t448 * t532 + (-qJD(2) * t493 + t537) * t618) * t522) * MDP(15) + (-t519 * t619 - t522 * t565) * MDP(16) * t618 + (t480 * t448 + t406 * t493 + t422 * t456 - t628 * t537 + (-t583 * t532 + (qJD(2) * t592 + t387) * t618) * t522 + t683 * t565) * MDP(17) + (-t388 * t599 + t406 * t494 + t422 * t455 + t628 * t475 + t480 * t535 - t545 * t651 + t565 * t690 - t575 * t629) * MDP(18) + (-t395 * t458 - t407 * t446) * MDP(19) + 0.2e1 * t681 * t607 + (MDP(7) * t598 - MDP(8) * t599) * (qJD(3) + 0.2e1 * t620); t565 * t573 + (t436 * t587 + t528 * t536 - t532 * t539) * MDP(10) + (t435 * t587 + t528 * t539 + t532 * t536) * MDP(11) + (-qJD(4) * t527 ^ 2 * t600 + ((qJD(4) * t587 + t574) * t527 - t565 * t475) * t531) * MDP(12) + (-t527 * t448 + t531 * t535 + t682 * t475 - (t579 - t613) * t537) * MDP(13) + (-t531 * t589 + (t531 * t591 + (qJD(3) * t527 - t475) * t528) * t622) * MDP(14) + (t527 * t589 + (-t527 * t591 + (qJD(3) * t531 - t537) * t528) * t622) * MDP(15) + (-pkin(3) * t448 + t527 * t609 - t431 * t565 + t436 * t537 + (pkin(10) * t589 + t484 * t565 - t406) * t531 + (-t387 * t528 + (-pkin(10) * t618 - t422 * t532) * t527) * t622) * MDP(17) + (-pkin(3) * t535 + t388 * t600 + t406 * t527 - t422 * t579 - t436 * t475 + t531 * t609 - t575 * t670 + (-t631 - t686) * t565) * MDP(18) + (-t395 * t527 * t530 + (-t527 * t612 + t568) * t446) * MDP(19) + (t444 * t467 + t446 * t466 + (-t444 * t530 - t446 * t526) * t613 + (t660 - t396 * t530 + (t444 * t526 - t446 * t530) * qJD(5)) * t527) * MDP(20) + (t395 * t531 + t568 * t469 + (-t446 * t565 - t597 + t655) * t527) * MDP(21) + (t396 * t531 - t569 * t469 + (t444 * t565 + t548) * t527) * MDP(22) + (-t469 * t527 * t565 - t448 * t531) * MDP(23) + (t448 * t653 - t381 * t466 - t402 * t444 + (-t687 * t530 + (-qJD(5) * t512 + t403) * t526) * t469 + (t381 * t616 - t348 + (qJD(4) * t444 + t548) * pkin(10)) * t531 + (t381 * t611 + t662 - t565 * t361 + (t469 * t616 + t396) * pkin(10)) * t527) * MDP(24) + (-t626 * t448 - t402 * t446 - t381 * t467 + t676 * t469 + (t381 * t614 + t347 + (qJD(4) * t446 + t597) * pkin(10)) * t531 + (-t381 * t612 + t360 * t530 + t565 * t362 + (t469 * t614 - t395) * pkin(10)) * t527) * MDP(25) + (t346 * t531 + t396 * t483 - t448 * t460 + t638 * t469 - t636 * t444 + t569 * t363 + (t355 * t565 + t363 * t611 + t665) * t527) * MDP(26) + (-t355 * t467 + t356 * t466 - t395 * t460 - t396 * t459 - t638 * t446 + t639 * t444 + t564 * t613 + (-t345 * t526 + t346 * t530 + (-t355 * t526 - t356 * t530) * qJD(5)) * t527) * MDP(27) + (-t345 * t531 + t395 * t483 + t448 * t459 - t639 * t469 + t636 * t446 - t568 * t363 + (-t356 * t565 + t363 * t612 - t664) * t527) * MDP(28) + (t345 * t459 + t346 * t460 + t349 * t483 - t355 * t638 - t356 * t639 - t363 * t636) * MDP(29) + ((-MDP(7) * t532 + MDP(8) * t528) * t522 * t524 - t681) * t534; -t537 ^ 2 * MDP(13) + (t537 * t514 + t558) * MDP(14) - t559 * MDP(15) + qJD(3) * t573 + (-t388 * t565 + t583) * MDP(17) + (-t387 * t565 - t422 * t537 + t545) * MDP(18) + (t446 * t590 - t660) * MDP(19) + ((-t395 - t659) * t530 + (-t396 - t657) * t526) * MDP(20) + (t469 * t590 + t656) * MDP(21) + (-t526 * t672 + t655) * MDP(22) + (-pkin(4) * t396 - t388 * t444 + (-t360 + (-pkin(11) * qJD(5) - t434) * t469) * t530 + (t387 * t469 + t547) * t526) * MDP(24) + (pkin(4) * t395 + t662 - t388 * t446 + (t605 + t635) * t469 + t547 * t530) * MDP(25) + (-t664 + t396 * t507 + (-pkin(11) * t611 + t365) * t469 - t637 * t444 + t675 * t526) * MDP(26) + (t364 * t444 - t365 * t446 + (t345 + t469 * t355 + (qJD(5) * t446 - t396) * pkin(11)) * t530 + (t346 - t663 + (qJD(5) * t444 - t395) * pkin(11)) * t526) * MDP(27) + (-t665 + t395 * t507 + (-t364 - t605) * t469 + t637 * t446 - t675 * t530) * MDP(28) + (t349 * t507 - t355 * t365 - t356 * t364 - t637 * t363 + (qJD(5) * t564 + t345 * t530 + t346 * t526) * pkin(11)) * MDP(29) + (-MDP(12) * t537 + MDP(13) * t475 - t514 * MDP(15) - t422 * MDP(17) - t446 * MDP(21) + t444 * MDP(22) - t469 * MDP(23) - t361 * MDP(24) + t362 * MDP(25) + t355 * MDP(26) - t356 * MDP(28)) * t475; MDP(19) * t658 + (-t444 ^ 2 + t673) * MDP(20) + t372 * MDP(21) + (-t396 + t657) * MDP(22) + t448 * MDP(23) + (-t381 * t446 + t348 + t661) * MDP(24) + (t361 * t469 + t381 * t444 - t347) * MDP(25) + (-t399 * t444 - t546 + t661 + 0.2e1 * t671) * MDP(26) + (pkin(5) * t395 - qJ(6) * t396 + (t356 - t362) * t446 + (t355 - t608) * t444) * MDP(27) + (0.2e1 * t666 - t363 * t444 + t399 * t446 + (0.2e1 * qJD(6) - t361) * t469 + t347) * MDP(28) + (-pkin(5) * t346 + qJ(6) * t345 - t355 * t362 + t356 * t608 - t363 * t399) * MDP(29); (-t448 + t658) * MDP(26) + t372 * MDP(27) + (-t672 - t673) * MDP(28) + (t546 - t663 - t671) * MDP(29);];
tauc  = t1;
