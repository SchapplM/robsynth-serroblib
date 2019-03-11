% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:21
% EndTime: 2019-03-09 05:57:30
% DurationCPUTime: 7.37s
% Computational Cost: add. (7230->519), mult. (14872->646), div. (0->0), fcn. (10279->14), ass. (0->232)
t535 = sin(qJ(3));
t534 = sin(qJ(4));
t538 = cos(qJ(3));
t645 = t534 * t538;
t681 = cos(qJ(4));
t479 = t535 * t681 + t645;
t620 = qJD(3) + qJD(4);
t443 = t620 * t479;
t611 = t681 * t538;
t623 = qJDD(1) * t535;
t426 = qJD(1) * t443 - qJDD(1) * t611 + t534 * t623;
t423 = qJDD(5) + t426;
t527 = qJ(1) + pkin(10);
t520 = sin(t527);
t521 = cos(t527);
t592 = g(1) * t521 + g(2) * t520;
t537 = cos(qJ(5));
t626 = qJD(5) * t537;
t609 = qJD(1) * t681;
t630 = qJD(1) * t535;
t467 = t534 * t630 - t538 * t609;
t655 = t467 * t537;
t699 = t626 + t655;
t533 = sin(qJ(5));
t627 = qJD(5) * t533;
t698 = t467 * t533 + t627;
t531 = sin(pkin(10));
t509 = pkin(1) * t531 + pkin(7);
t673 = pkin(8) + t509;
t530 = qJ(3) + qJ(4);
t524 = sin(t530);
t651 = t521 * t524;
t653 = t520 * t524;
t695 = g(1) * t651 + g(2) * t653;
t532 = cos(pkin(10));
t510 = -t532 * pkin(1) - pkin(2);
t526 = t538 * pkin(3);
t687 = t510 - t526;
t694 = -pkin(9) * t479 + t687;
t525 = cos(t530);
t693 = t592 * t525;
t466 = qJD(5) + t467;
t610 = t466 * t627;
t644 = t537 * t423;
t692 = pkin(9) * (t610 - t644);
t603 = t673 * qJD(1);
t454 = qJD(2) * t535 + t538 * t603;
t522 = t538 * qJDD(2);
t490 = t509 * qJDD(1);
t602 = pkin(8) * qJDD(1) + t490;
t417 = qJDD(3) * pkin(3) - qJD(3) * t454 - t602 * t535 + t522;
t453 = t538 * qJD(2) - t603 * t535;
t424 = qJD(3) * t453 + t535 * qJDD(2) + t602 * t538;
t671 = qJD(3) * pkin(3);
t448 = t453 + t671;
t608 = t681 * qJD(4);
t628 = qJD(4) * t534;
t556 = t534 * t417 + t424 * t681 + t448 * t608 - t454 * t628;
t619 = qJDD(3) + qJDD(4);
t362 = pkin(9) * t619 + t556;
t624 = qJD(1) * qJD(3);
t607 = t535 * t624;
t507 = pkin(3) * t607;
t573 = -t534 * t535 + t611;
t442 = t620 * t573;
t547 = t442 * qJD(1);
t375 = t426 * pkin(4) - pkin(9) * t547 + qJDD(1) * t694 + t507;
t447 = t681 * t454;
t412 = t534 * t448 + t447;
t406 = pkin(9) * t620 + t412;
t469 = -qJD(1) * t645 - t535 * t609;
t470 = t687 * qJD(1);
t429 = pkin(4) * t467 + pkin(9) * t469 + t470;
t569 = t537 * t362 + t533 * t375 - t406 * t627 + t429 * t626;
t670 = qJ(6) * t423;
t351 = qJD(6) * t466 + t569 + t670;
t598 = t533 * t362 - t537 * t375 + t406 * t626 + t429 * t627;
t680 = pkin(5) * t423;
t353 = qJDD(6) + t598 - t680;
t584 = t351 * t537 + t353 * t533;
t545 = t479 * qJDD(1) + t547;
t600 = t537 * t620;
t391 = -qJD(5) * t600 - t469 * t627 - t533 * t619 - t537 * t545;
t567 = t537 * t469 - t533 * t620;
t641 = t391 * t573 - t443 * t567;
t392 = -qJD(5) * t567 + t533 * t545 - t537 * t619;
t449 = -t469 * t533 - t600;
t689 = -t392 * t573 + t443 * t449;
t418 = t534 * t453 + t447;
t595 = pkin(3) * t628 - t418;
t435 = -pkin(4) * t573 + t694;
t472 = t673 * t535;
t473 = t673 * t538;
t438 = -t534 * t472 + t473 * t681;
t635 = t533 * t435 + t537 * t438;
t688 = t698 * pkin(5) - qJ(6) * t699 - qJD(6) * t533;
t686 = t525 * pkin(4) + t524 * pkin(9);
t684 = MDP(24) + MDP(26);
t683 = t567 ^ 2;
t682 = t466 ^ 2;
t679 = pkin(5) * t469;
t678 = pkin(9) * t423;
t513 = g(3) * t524;
t674 = g(3) * t525;
t672 = pkin(9) * qJD(5);
t597 = -t681 * t417 + t534 * t424 + t448 * t628 + t454 * t608;
t363 = -pkin(4) * t619 + t597;
t354 = t392 * pkin(5) + t391 * qJ(6) + qJD(6) * t567 + t363;
t669 = t354 * t533;
t380 = t406 * t537 + t429 * t533;
t668 = t380 * t466;
t446 = t534 * t454;
t411 = t448 * t681 - t446;
t405 = -pkin(4) * t620 - t411;
t381 = t449 * pkin(5) + qJ(6) * t567 + t405;
t667 = t381 * t467;
t666 = t391 * t533;
t665 = t405 * t467;
t516 = pkin(3) * t534 + pkin(9);
t664 = t423 * t516;
t663 = t442 * t533;
t662 = t442 * t537;
t661 = t449 * t466;
t660 = t449 * t533;
t659 = t567 * t449;
t658 = t567 * t466;
t657 = t567 * t537;
t601 = t466 * t537;
t654 = t479 * t537;
t649 = t525 * t533;
t648 = t525 * t537;
t647 = t533 * qJ(6);
t646 = t533 * t423;
t643 = qJDD(2) - g(3);
t642 = -t392 * t654 - t449 * t662;
t639 = -t688 - t595;
t638 = t442 * t601 + t479 * t644;
t439 = -pkin(4) * t469 + pkin(9) * t467;
t637 = t537 * t411 + t533 * t439;
t419 = t453 * t681 - t446;
t433 = pkin(3) * t630 + t439;
t636 = t537 * t419 + t533 * t433;
t634 = -t412 + t688;
t633 = t695 * t537;
t528 = t535 ^ 2;
t632 = -t538 ^ 2 + t528;
t631 = MDP(27) * t533;
t493 = qJD(1) * t510;
t379 = -t406 * t533 + t429 * t537;
t625 = qJD(6) - t379;
t622 = qJDD(1) * t538;
t618 = t681 * pkin(3);
t617 = t535 * t671;
t615 = t513 + t693;
t614 = g(3) * t649 - t533 * t695;
t613 = t533 * t681;
t612 = t537 * t681;
t606 = -t354 - t674;
t605 = -t363 - t674;
t604 = qJD(3) * t673;
t599 = pkin(3) * t608;
t596 = pkin(5) * t648 + t525 * t647 + t686;
t457 = t520 * t649 + t521 * t537;
t459 = -t520 * t537 + t521 * t649;
t594 = -g(1) * t457 + g(2) * t459;
t458 = t520 * t648 - t521 * t533;
t460 = t520 * t533 + t521 * t648;
t593 = g(1) * t458 - g(2) * t460;
t591 = g(1) * t520 - g(2) * t521;
t536 = sin(qJ(1));
t539 = cos(qJ(1));
t590 = g(1) * t536 - g(2) * t539;
t588 = t537 * pkin(5) + t647;
t587 = pkin(5) * t533 - qJ(6) * t537;
t369 = -pkin(5) * t466 + t625;
t370 = qJ(6) * t466 + t380;
t583 = t369 * t537 - t370 * t533;
t582 = t369 * t533 + t370 * t537;
t581 = -t664 + t665;
t580 = -t657 + t660;
t578 = pkin(4) + t588;
t577 = t526 + pkin(2) + t686;
t576 = -t369 * t469 + t381 * t627 + t633;
t575 = t379 * t469 + t405 * t627 + t633;
t574 = -t472 * t681 - t534 * t473;
t571 = t479 * t626 + t663;
t570 = t479 * t627 - t662;
t463 = t535 * t604;
t464 = t538 * t604;
t394 = qJD(4) * t574 - t463 * t681 - t534 * t464;
t400 = pkin(4) * t443 - pkin(9) * t442 + t617;
t568 = t537 * t394 + t533 * t400 + t435 * t626 - t438 * t627;
t566 = -t537 * MDP(25) - t533 * t684;
t565 = t363 * t533 - t380 * t469 + t405 * t626 + t614;
t563 = -qJD(1) * t493 - t490 + t592;
t562 = 0.2e1 * qJD(3) * t493 - qJDD(3) * t509;
t561 = t370 * t469 - t381 * t655 - t614 - t669;
t560 = -t516 * t627 + t537 * t599;
t559 = g(1) * t459 + g(2) * t457 + t533 * t513 - t598;
t558 = t470 * t469 - t597 - t674 + t695;
t541 = qJD(3) ^ 2;
t557 = -0.2e1 * qJDD(1) * t510 - t509 * t541 + t591;
t555 = qJD(5) * t583 + t584;
t554 = qJD(5) * t580 - t392 * t537 - t666;
t553 = -t381 * t567 + qJDD(6) - t559;
t552 = t369 * t699 - t698 * t370 + t584 - t615;
t551 = -g(1) * t460 - g(2) * t458 - t513 * t537 + t569;
t395 = qJD(4) * t438 - t534 * t463 + t464 * t681;
t550 = ((-t391 - t661) * t537 + (-t392 + t658) * t533) * MDP(20) + (-t567 * t601 - t666) * MDP(19) + (-t449 * t469 - t533 * t682 + t644) * MDP(22) + (t466 * t601 - t469 * t567 + t646) * MDP(21) + (t467 * t620 + t545) * MDP(14) + (-t469 * t620 - t426) * MDP(15) + (-t467 ^ 2 + t469 ^ 2) * MDP(13) + t619 * MDP(16) + (-MDP(12) * t467 + MDP(23) * t466) * t469;
t549 = t470 * t467 - t556 + t615;
t546 = t592 * t524 * t578 - pkin(9) * t693;
t540 = -pkin(8) - pkin(7);
t517 = -t618 - pkin(4);
t489 = qJDD(3) * t538 - t535 * t541;
t488 = qJDD(3) * t535 + t538 * t541;
t471 = -t618 - t578;
t462 = t469 * qJ(6);
t455 = qJDD(1) * t687 + t507;
t425 = -t443 * t620 + t573 * t619;
t415 = -pkin(5) * t567 + qJ(6) * t449;
t399 = t479 * t587 - t574;
t389 = pkin(5) * t573 - t435 * t537 + t438 * t533;
t388 = -qJ(6) * t573 + t635;
t383 = t411 * t533 - t439 * t537 + t679;
t382 = -t462 + t637;
t378 = t419 * t533 - t433 * t537 + t679;
t377 = -t462 + t636;
t372 = -t391 + t661;
t358 = t587 * t442 + (qJD(5) * t588 - qJD(6) * t537) * t479 + t395;
t357 = -pkin(5) * t443 + qJD(5) * t635 + t394 * t533 - t400 * t537;
t356 = qJ(6) * t443 - qJD(6) * t573 + t568;
t1 = [(-t356 * t449 - t357 * t567 - t388 * t392 - t389 * t391 + t591 * t524 + t583 * t442 + (-qJD(5) * t582 - t351 * t533 + t353 * t537) * t479) * MDP(27) + (-t391 * t654 + t567 * t570) * MDP(19) + (t567 * t663 + (t666 + (t657 + t660) * qJD(5)) * t479 + t642) * MDP(20) + (t598 * t573 + t379 * t443 + t395 * t449 - t574 * t392 + ((-qJD(5) * t438 + t400) * t466 + t435 * t423 + t405 * qJD(5) * t479) * t537 + ((-qJD(5) * t435 - t394) * t466 - t438 * t423 + t363 * t479 + t405 * t442) * t533 + t593) * MDP(24) + (t363 * t654 - t380 * t443 + t391 * t574 - t395 * t567 - t405 * t570 - t423 * t635 - t466 * t568 + t569 * t573 + t594) * MDP(25) + (-t351 * t573 - t354 * t654 + t356 * t466 + t358 * t567 + t370 * t443 + t381 * t570 + t388 * t423 + t391 * t399 - t594) * MDP(28) + (-t479 * t426 - t442 * t467 + t469 * t443 + t545 * t573) * MDP(13) + (t353 * t573 - t357 * t466 + t358 * t449 - t369 * t443 + t381 * t571 - t389 * t423 + t392 * t399 + t479 * t669 + t593) * MDP(26) + (-t423 * t573 + t443 * t466) * MDP(23) + (-t395 * t620 + t426 * t687 + t470 * t443 - t455 * t573 + t467 * t617 + t525 * t591 + t574 * t619) * MDP(17) + (-g(1) * t653 + g(2) * t651 - t394 * t620 - t438 * t619 + t470 * t442 + t455 * t479 - t469 * t617 + t545 * t687) * MDP(18) + (t535 * t562 + t538 * t557) * MDP(10) + (-t535 * t557 + t538 * t562) * MDP(11) + (-t466 * t571 - t479 * t646 - t689) * MDP(22) + t590 * MDP(2) + (-t469 * t442 + t479 * t545) * MDP(12) + qJDD(1) * MDP(1) + (t351 * t388 + t370 * t356 + t354 * t399 + t381 * t358 + t353 * t389 + t369 * t357 - g(1) * (-pkin(1) * t536 - pkin(5) * t458 - qJ(6) * t457) - g(2) * (pkin(1) * t539 + pkin(5) * t460 + qJ(6) * t459) + (g(1) * t540 - g(2) * t577) * t521 + (g(1) * t577 + g(2) * t540) * t520) * MDP(29) + (qJDD(1) * t528 + 0.2e1 * t538 * t607) * MDP(5) + (t590 + (t531 ^ 2 + t532 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t442 * t620 + t479 * t619) * MDP(14) + 0.2e1 * (t535 * t622 - t624 * t632) * MDP(6) + (-t479 * t610 + t638 + t641) * MDP(21) + t425 * MDP(15) + t488 * MDP(7) + t489 * MDP(8) + (g(1) * t539 + g(2) * t536) * MDP(3); t643 * MDP(4) + t489 * MDP(10) - t488 * MDP(11) + t425 * MDP(17) + t641 * MDP(25) + t642 * MDP(27) + (t638 - t641) * MDP(28) + (-t354 * t573 + t381 * t443 - g(3)) * MDP(29) + (-MDP(18) * t620 + MDP(29) * t582 + t466 * t566 - t567 * t631) * t442 + (-t619 * MDP(18) - t391 * t631 + t584 * MDP(29) + t566 * t423 + (t580 * MDP(27) + t583 * MDP(29) + (-t684 * t537 + (MDP(25) - MDP(28)) * t533) * t466) * qJD(5)) * t479 + t684 * t689; MDP(8) * t622 + (t471 * t392 + t606 * t537 + (-t664 + t667) * t533 - t639 * t449 + (-t516 * t626 - t533 * t599 + t378) * t466 + t576) * MDP(26) + (t419 * t620 + (t469 * t630 - t534 * t619 - t608 * t620) * pkin(3) + t549) * MDP(18) + (t377 * t449 + t378 * t567 + (-t449 * t612 - t567 * t613) * qJD(4) * pkin(3) + t554 * t516 + t552) * MDP(27) + (t471 * t391 + (-qJD(5) * t381 + t664) * t537 - t639 * t567 + (-t377 + t560) * t466 + t561) * MDP(28) + t550 + (t517 * t392 + t605 * t537 + t581 * t533 + t595 * t449 + ((-qJD(5) * t516 - t433) * t537 + (-t599 + t419) * t533) * t466 + t575) * MDP(24) + MDP(7) * t623 + (t354 * t471 - t370 * t377 - t369 * t378 - g(3) * (t526 + t596) - t639 * t381 + (t592 * t535 + (t369 * t613 + t370 * t612) * qJD(4)) * pkin(3) + t555 * t516 + t546) * MDP(29) + (-t517 * t391 + t581 * t537 - t595 * t567 + (-t560 + t636) * t466 + t565) * MDP(25) + (t418 * t620 + (-t467 * t630 + t619 * t681 - t620 * t628) * pkin(3) + t558) * MDP(17) + (-t535 * t643 + t538 * t563) * MDP(11) + qJDD(3) * MDP(9) + (-g(3) * t538 + t535 * t563 + t522) * MDP(10) + (-MDP(5) * t535 * t538 + MDP(6) * t632) * qJD(1) ^ 2; (t411 * t620 + t549) * MDP(18) + (pkin(9) * t554 + t382 * t449 + t383 * t567 + t552) * MDP(27) + (-t381 * t626 - t382 * t466 - t391 * t578 + t567 * t634 + t561 - t692) * MDP(28) + t550 + (pkin(9) * t555 - g(3) * t596 - t354 * t578 - t369 * t383 - t370 * t382 + t381 * t634 + t546) * MDP(29) + (t383 * t466 - t392 * t578 + (t667 - t678) * t533 + t634 * t449 + (-t466 * t672 + t606) * t537 + t576) * MDP(26) + (-pkin(4) * t392 - t412 * t449 + (t411 * t466 + t665 - t678) * t533 + ((-t439 - t672) * t466 + t605) * t537 + t575) * MDP(24) + (t412 * t620 + t558) * MDP(17) + (pkin(4) * t391 + t405 * t655 + t412 * t567 + t466 * t637 + t565 + t692) * MDP(25); -MDP(19) * t659 + (-t449 ^ 2 + t683) * MDP(20) + t372 * MDP(21) + (-t392 - t658) * MDP(22) + t423 * MDP(23) + (t405 * t567 + t559 + t668) * MDP(24) + (t379 * t466 + t405 * t449 - t551) * MDP(25) + (-t415 * t449 - t553 + t668 + 0.2e1 * t680) * MDP(26) + (pkin(5) * t391 - qJ(6) * t392 - (t370 - t380) * t567 + (t369 - t625) * t449) * MDP(27) + (0.2e1 * t670 - t381 * t449 - t415 * t567 + (0.2e1 * qJD(6) - t379) * t466 + t551) * MDP(28) + (t351 * qJ(6) - t353 * pkin(5) - t381 * t415 - t369 * t380 - g(1) * (-pkin(5) * t459 + qJ(6) * t460) - g(2) * (-pkin(5) * t457 + qJ(6) * t458) + t587 * t513 + t625 * t370) * MDP(29); t372 * MDP(27) + (-t682 - t683) * MDP(28) + (-t370 * t466 + t553 - t680) * MDP(29) + (-t659 - t423) * MDP(26);];
tau  = t1;
