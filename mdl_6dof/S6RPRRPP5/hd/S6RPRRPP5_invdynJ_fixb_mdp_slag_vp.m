% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:55
% EndTime: 2019-03-09 04:45:07
% DurationCPUTime: 8.57s
% Computational Cost: add. (6498->576), mult. (14957->652), div. (0->0), fcn. (11118->10), ass. (0->239)
t550 = sin(pkin(9));
t690 = pkin(7) + qJ(2);
t511 = t690 * t550;
t502 = qJD(1) * t511;
t551 = cos(pkin(9));
t512 = t690 * t551;
t503 = qJD(1) * t512;
t554 = sin(qJ(3));
t698 = cos(qJ(3));
t447 = -t554 * t502 + t698 * t503;
t728 = qJD(3) * t447;
t628 = t698 * t551;
t587 = -t554 * t550 + t628;
t727 = t587 * qJD(1);
t556 = cos(qJ(4));
t553 = sin(qJ(4));
t644 = qJD(3) * t553;
t501 = t550 * t698 + t554 * t551;
t708 = t501 * qJD(1);
t459 = t556 * t708 + t644;
t480 = qJD(4) - t727;
t672 = t480 * t553;
t726 = t459 * t672;
t495 = t501 * qJD(3);
t634 = qJDD(1) * t550;
t604 = -qJDD(1) * t628 + t554 * t634;
t445 = qJD(1) * t495 + t604;
t439 = qJDD(4) + t445;
t725 = t439 * qJ(5) + t480 * qJD(5);
t548 = pkin(9) + qJ(3);
t540 = sin(t548);
t557 = cos(qJ(1));
t665 = t540 * t557;
t555 = sin(qJ(1));
t667 = t540 * t555;
t713 = g(1) * t665 + g(2) * t667;
t683 = qJDD(1) * pkin(1);
t539 = qJDD(2) - t683;
t711 = g(1) * t555 - g(2) * t557;
t597 = -t539 + t711;
t494 = t587 * qJD(3);
t575 = t501 * qJDD(1);
t563 = qJD(1) * t494 + t575;
t724 = -qJD(3) * qJD(4) - t563;
t723 = qJD(3) * t708;
t722 = MDP(24) + MDP(27);
t661 = t553 * qJ(5);
t699 = pkin(4) + pkin(5);
t591 = -t556 * t699 - t661;
t498 = pkin(3) - t591;
t636 = qJD(1) * qJD(2);
t701 = qJDD(1) * t690 + t636;
t467 = t701 * t550;
t468 = t701 * t551;
t589 = -t554 * t467 + t468 * t698;
t715 = -t698 * t502 - t554 * t503;
t399 = qJDD(3) * pkin(8) + qJD(3) * t715 + t589;
t633 = qJDD(1) * t551;
t401 = -pkin(2) * t633 + t445 * pkin(3) - pkin(8) * t563 + t539;
t535 = pkin(2) * t551 + pkin(1);
t509 = -qJD(1) * t535 + qJD(2);
t419 = -pkin(3) * t727 - pkin(8) * t708 + t509;
t441 = qJD(3) * pkin(8) + t447;
t642 = qJD(4) * t556;
t643 = qJD(4) * t553;
t617 = t553 * t399 - t556 * t401 + t419 * t643 + t441 * t642;
t435 = t439 * pkin(4);
t705 = t435 - qJDD(5);
t367 = t617 - t705;
t391 = t553 * t419 + t556 * t441;
t470 = t480 * qJ(5);
t385 = t470 + t391;
t681 = t385 * t480;
t721 = -t367 + t681;
t607 = -t556 * qJDD(3) + t642 * t708;
t410 = t553 * t575 + (qJD(4) + t727) * t644 + t607;
t457 = -t556 * qJD(3) + t553 * t708;
t720 = t410 * qJ(6) + t457 * qJD(6);
t409 = -t553 * qJDD(3) + t556 * t724 + t708 * t643;
t678 = t457 * t480;
t572 = t409 - t678;
t719 = 0.2e1 * t725;
t718 = qJD(5) * t553 + t447;
t716 = t495 * qJ(5) - qJD(5) * t587;
t714 = -t698 * t511 - t554 * t512;
t541 = cos(t548);
t712 = -t541 * pkin(3) - t540 * pkin(8);
t609 = g(1) * t557 + g(2) * t555;
t453 = -t554 * t511 + t512 * t698;
t709 = t453 * qJD(3);
t390 = t556 * t419 - t553 * t441;
t638 = qJD(5) - t390;
t707 = qJ(2) * qJDD(1);
t706 = MDP(22) + MDP(26);
t427 = t556 * t439;
t671 = t727 * t553;
t593 = t427 + (-t643 + t671) * t480;
t581 = t556 * t399 + t553 * t401 + t419 * t642 - t441 * t643;
t366 = t581 + t725;
t364 = t366 + t720;
t380 = qJ(6) * t459 + t390;
t639 = qJD(5) - t380;
t374 = -t480 * t699 + t639;
t704 = t374 * t480 + t364;
t612 = qJD(3) * pkin(3) + t715;
t578 = qJ(5) * t459 + t612;
t392 = pkin(4) * t457 - t578;
t696 = pkin(8) * t439;
t703 = t392 * t480 - t696;
t379 = -t457 * t699 + qJD(6) + t578;
t656 = t556 * t557;
t660 = t553 * t555;
t481 = t541 * t660 + t656;
t655 = t557 * t553;
t657 = t555 * t556;
t483 = t541 * t655 - t657;
t668 = t540 * t553;
t571 = g(1) * t483 + g(2) * t481 + g(3) * t668 - t617;
t568 = t571 + t705;
t684 = qJ(6) * t409;
t702 = (qJD(6) + t379) * t459 + t568 - t684;
t700 = t457 ^ 2;
t456 = t459 ^ 2;
t475 = t480 ^ 2;
t697 = pkin(5) * t439;
t693 = g(2) * t690;
t691 = g(3) * t541;
t689 = pkin(8) - qJ(6);
t688 = pkin(8) * qJD(4);
t687 = qJ(5) * t410;
t686 = qJ(5) * t457;
t685 = qJ(5) * t556;
t381 = qJ(6) * t457 + t391;
t376 = t381 + t470;
t682 = t376 * t480;
t680 = t391 * t480;
t679 = t409 * t553;
t677 = t457 * t727;
t676 = t457 * t708;
t675 = t459 * t457;
t674 = t459 * t480;
t673 = t459 * t708;
t618 = t480 * t556;
t670 = t501 * t553;
t669 = t501 * t556;
t666 = t540 * t556;
t664 = t541 * t555;
t663 = t541 * t556;
t662 = t541 * t557;
t425 = t553 * t439;
t654 = -t553 * t410 - t457 * t642;
t442 = pkin(3) * t708 - pkin(8) * t727;
t653 = t553 * t442 + t556 * t715;
t444 = -pkin(3) * t587 - pkin(8) * t501 - t535;
t652 = t553 * t444 + t556 * t453;
t632 = t699 * t553;
t592 = -t632 + t685;
t651 = t480 * t592 + t718;
t387 = qJ(5) * t708 + t653;
t650 = -qJ(6) * t671 - qJD(6) * t556 - t643 * t689 - t387;
t605 = pkin(4) * t553 - t685;
t649 = t480 * t605 - t718;
t437 = t553 * t715;
t514 = t689 * t556;
t648 = qJD(4) * t514 - qJD(6) * t553 - t437 - (-qJ(6) * t727 - t442) * t556 + t699 * t708;
t647 = t713 * t553;
t646 = t713 * t556;
t645 = t550 ^ 2 + t551 ^ 2;
t640 = qJD(5) * t556;
t420 = t587 * qJD(2) + qJD(3) * t714;
t631 = t553 * t420 + t444 * t643 + t453 * t642;
t443 = pkin(3) * t495 - pkin(8) * t494;
t630 = t556 * t420 + t553 * t443 + t444 * t642;
t397 = -qJ(5) * t587 + t652;
t629 = g(1) * t662 + g(2) * t664 + g(3) * t540;
t627 = t501 * t643;
t626 = t501 * t642;
t616 = -t698 * t467 - t554 * t468 - t728;
t400 = -qJDD(3) * pkin(3) - t616;
t368 = t410 * pkin(4) + t409 * qJ(5) - t459 * qJD(5) + t400;
t365 = -pkin(5) * t410 + qJDD(6) - t368;
t624 = t365 - t691;
t623 = t645 * qJD(1) ^ 2;
t482 = t541 * t657 - t655;
t622 = -t481 * pkin(4) + qJ(5) * t482;
t484 = t541 * t656 + t660;
t621 = -t483 * pkin(4) + qJ(5) * t484;
t449 = t553 * t453;
t620 = t444 * t556 - t449;
t615 = pkin(4) * t663 + t541 * t661 - t712;
t614 = 0.2e1 * t645;
t613 = -g(1) * t667 + g(2) * t665;
t611 = g(1) * t481 - g(2) * t483;
t610 = g(1) * t482 - g(2) * t484;
t606 = pkin(4) * t556 + t661;
t384 = -pkin(4) * t480 + t638;
t603 = t384 * t556 - t385 * t553;
t363 = -qJD(6) * t459 + t367 + t684 - t697;
t602 = -t363 + t682;
t601 = t384 * t480 + t366;
t600 = -qJ(6) * t494 - qJD(6) * t501;
t599 = pkin(3) + t606;
t598 = -t535 + t712;
t596 = t480 * t642 - t618 * t727 + t425;
t595 = t443 * t556 - t631;
t594 = -t480 * t688 - t691;
t590 = -t482 * pkin(4) - qJ(5) * t481 + t557 * t690;
t586 = t494 * t553 + t626;
t585 = -t494 * t556 + t627;
t584 = -t480 * t612 - t696;
t583 = pkin(3) * t662 + t484 * pkin(4) + pkin(8) * t665 + qJ(5) * t483 + t557 * t535;
t582 = -t368 + t594;
t580 = -t453 * t643 + t630;
t569 = t614 * t636 - t609;
t566 = t392 * t459 - t568;
t565 = g(1) * t484 + g(2) * t482 + g(3) * t666 - t581;
t421 = t501 * qJD(2) + t709;
t564 = t390 * t480 + t565;
t561 = t553 * t724 - t607;
t521 = pkin(8) * t662;
t518 = pkin(8) * t664;
t515 = qJ(5) * t666;
t513 = t689 * t553;
t508 = -qJDD(1) * t535 + qJDD(2);
t412 = pkin(4) * t459 + t686;
t411 = t501 * t605 - t714;
t403 = -t459 * t699 - t686;
t402 = t501 * t592 + t714;
t398 = pkin(4) * t587 - t620;
t389 = -pkin(4) * t708 - t442 * t556 + t437;
t386 = qJ(6) * t670 + t397;
t378 = t449 + (-qJ(6) * t501 - t444) * t556 + t699 * t587;
t375 = t605 * t494 + (qJD(4) * t606 - t640) * t501 + t421;
t373 = -pkin(4) * t495 - t595;
t372 = -t709 + t592 * t494 + (qJD(4) * t591 - qJD(2) + t640) * t501;
t371 = t580 + t716;
t370 = qJ(6) * t626 + (-qJD(4) * t453 - t600) * t553 + t630 + t716;
t369 = qJ(6) * t627 - t699 * t495 + (-t443 + t600) * t556 + t631;
t1 = [(-t391 * t495 + t400 * t669 + t409 * t714 + t421 * t459 - t439 * t652 - t480 * t580 + t581 * t587 + t585 * t612 - t611) * MDP(21) + (t390 * t495 + t400 * t670 - t410 * t714 + t421 * t457 + t439 * t620 + t480 * t595 - t586 * t612 + t587 * t617 + t610) * MDP(20) + (-qJD(3) * t421 + qJDD(3) * t714 - t445 * t535 + t495 * t509 - t508 * t587 + t541 * t711) * MDP(13) + (-t501 * t445 + t494 * t727 - t495 * t708 + t563 * t587) * MDP(9) + (t366 * t397 + t385 * t371 + t368 * t411 + t392 * t375 + t367 * t398 + t384 * t373 - g(1) * t590 - g(2) * t583 + (-g(1) * t598 - t693) * t555) * MDP(25) + (t364 * t386 + t376 * t370 + t363 * t378 + t374 * t369 + t365 * t402 + t379 * t372 - g(1) * (-pkin(5) * t482 + t590) - g(2) * (pkin(5) * t484 - qJ(6) * t665 + t583) + (-g(1) * (qJ(6) * t540 + t598) - t693) * t555) * MDP(29) + (qJD(3) * t494 + qJDD(3) * t501) * MDP(10) + t609 * MDP(3) + (-t409 * t669 - t459 * t585) * MDP(15) + ((-t457 * t556 - t459 * t553) * t494 + (t679 - t410 * t556 + (t457 * t553 - t459 * t556) * qJD(4)) * t501) * MDP(16) + (t551 * MDP(4) - t550 * MDP(5)) * (t597 + t683) + qJDD(1) * MDP(1) + (pkin(1) * t597 + (t645 * t707 + t569) * qJ(2)) * MDP(7) + (t614 * t707 + t569) * MDP(6) + (t494 * t708 + t501 * t563) * MDP(8) + (t363 * t587 - t365 * t670 - t369 * t480 - t372 * t457 - t374 * t495 - t378 * t439 - t379 * t586 - t402 * t410 + t610) * MDP(26) + (t367 * t587 + t368 * t670 - t373 * t480 + t375 * t457 - t384 * t495 + t392 * t586 - t398 * t439 + t410 * t411 + t610) * MDP(22) + (t410 * t587 - t425 * t501 - t457 * t495 - t480 * t586) * MDP(18) + (-t439 * t587 + t480 * t495) * MDP(19) + (-qJD(3) * t495 + qJDD(3) * t587) * MDP(11) + (t409 * t587 + t427 * t501 + t459 * t495 - t480 * t585) * MDP(17) + (-t366 * t587 - t368 * t669 + t371 * t480 - t375 * t459 + t385 * t495 + t392 * t585 + t397 * t439 + t409 * t411 + t611) * MDP(24) + (-t364 * t587 + t365 * t669 + t370 * t480 + t372 * t459 + t376 * t495 - t379 * t585 + t386 * t439 - t402 * t409 + t611) * MDP(27) + t711 * MDP(2) + (-t369 * t459 + t370 * t457 + t378 * t409 + t386 * t410 + (-t374 * t556 + t376 * t553) * t494 + (-t363 * t556 + t364 * t553 + (t374 * t553 + t376 * t556) * qJD(4)) * t501 + t613) * MDP(28) + (-t371 * t457 + t373 * t459 - t397 * t410 - t398 * t409 + t603 * t494 + (-t366 * t553 + t367 * t556 + (-t384 * t553 - t385 * t556) * qJD(4)) * t501 - t613) * MDP(23) + (-t420 * qJD(3) - t453 * qJDD(3) + t509 * t494 + t508 * t501 - t535 * t563 + t613) * MDP(14); -MDP(4) * t633 + MDP(5) * t634 - MDP(6) * t623 + (-qJ(2) * t623 - t597) * MDP(7) + (t604 + 0.2e1 * t723) * MDP(13) + (0.2e1 * t727 * qJD(3) + t575) * MDP(14) + (-t475 * t556 - t425 - t673) * MDP(21) + ((t409 + t677) * t556 + t726 + t654) * MDP(23) + (-t392 * t708 + t601 * t553 + t556 * t721 - t711) * MDP(25) + (-t572 * t556 + (t410 - t674) * t553) * MDP(28) + (t379 * t708 + t553 * t704 + t602 * t556 - t711) * MDP(29) + t722 * (t596 + t673) + (MDP(20) + t706) * (t593 - t676); -t727 ^ 2 * MDP(9) + t575 * MDP(10) - t604 * MDP(11) + qJDD(3) * MDP(12) + (t616 - t691 + t713 + t728) * MDP(13) + (-t509 * t727 - t589 + t629) * MDP(14) + (t459 * t618 - t679) * MDP(15) + ((-t409 + t677) * t556 - t726 + t654) * MDP(16) + (t596 - t673) * MDP(17) + (t593 + t676) * MDP(18) + (-pkin(3) * t410 + t437 * t480 - t447 * t457 + (-t691 - t400 + (-t442 - t688) * t480) * t556 + t584 * t553 + t646) * MDP(20) + (pkin(3) * t409 + t653 * t480 - t447 * t459 + t584 * t556 + (t400 - t594) * t553 - t647) * MDP(21) + (t389 * t480 - t410 * t599 + t649 * t457 + t553 * t703 + t582 * t556 + t646) * MDP(22) + (t387 * t457 - t389 * t459 + ((qJD(4) * t459 - t410) * pkin(8) + t601) * t556 + ((qJD(4) * t457 - t409) * pkin(8) - t721) * t553 - t629) * MDP(23) + (-t387 * t480 - t409 * t599 - t649 * t459 + t582 * t553 - t556 * t703 + t647) * MDP(24) + (-t385 * t387 - t384 * t389 - g(1) * t521 - g(2) * t518 - g(3) * t615 + t649 * t392 + (qJD(4) * t603 + t366 * t556 + t367 * t553) * pkin(8) + (t540 * t609 - t368) * t599) * MDP(25) + (-t379 * t672 - t410 * t498 - t439 * t513 - t457 * t651 - t480 * t648 + t556 * t624 + t646) * MDP(26) + (t379 * t618 - t409 * t498 + t439 * t514 + t459 * t651 + t480 * t650 + t553 * t624 + t647) * MDP(27) + (t409 * t513 + t410 * t514 + t650 * t457 - t648 * t459 + t602 * t553 - t556 * t704 + t629) * MDP(28) + (t364 * t514 + t363 * t513 + t365 * t498 - g(1) * (-qJ(6) * t662 + t521) - g(2) * (-qJ(6) * t664 + t518) - g(3) * (pkin(5) * t663 + t615) + t651 * t379 + t650 * t376 + t648 * t374 + (g(3) * qJ(6) + t498 * t609) * t540) * MDP(29) + (-t509 * MDP(13) - t480 * MDP(19) - t390 * MDP(20) + t391 * MDP(21) + t384 * MDP(22) - t385 * MDP(24) + t374 * MDP(26) - t376 * MDP(27) - MDP(8) * t727 + MDP(9) * t708) * t708; MDP(15) * t675 + (t456 - t700) * MDP(16) - t572 * MDP(17) + (t561 + t674) * MDP(18) + t439 * MDP(19) + (t459 * t612 + t571 + t680) * MDP(20) + (-t457 * t612 + t564) * MDP(21) + (-t412 * t457 + t435 - t566 + t680) * MDP(22) + (pkin(4) * t409 - t687 + (t385 - t391) * t459 + (t384 - t638) * t457) * MDP(23) + (-t392 * t457 + t412 * t459 - t564 + t719) * MDP(24) + (t366 * qJ(5) - t367 * pkin(4) - t392 * t412 - t384 * t391 - g(1) * t621 - g(2) * t622 - g(3) * (-pkin(4) * t668 + t515) + t638 * t385) * MDP(25) + ((pkin(5) + t699) * t439 + t381 * t480 + t403 * t457 + t702) * MDP(26) + (t379 * t457 - t380 * t480 - t403 * t459 - t565 + t719 + t720) * MDP(27) + (t687 - t409 * t699 + (-t376 + t381) * t459 + (-t374 + t639) * t457) * MDP(28) + (t364 * qJ(5) - t363 * t699 - t374 * t381 - t379 * t403 - g(1) * (-pkin(5) * t483 + t621) - g(2) * (-pkin(5) * t481 + t622) - g(3) * (-t540 * t632 + t515) + t639 * t376) * MDP(29); (t566 - t681) * MDP(25) + (-t682 - t697 - t702) * MDP(29) + t706 * (-qJDD(4) - t604 + t675 - t723) + t722 * (-t475 - t456) + (-MDP(23) + MDP(28)) * t572; (t561 - t674) * MDP(26) + (-t409 - t678) * MDP(27) + (-t456 - t700) * MDP(28) + (t374 * t459 - t376 * t457 + t624 + t713) * MDP(29);];
tau  = t1;
