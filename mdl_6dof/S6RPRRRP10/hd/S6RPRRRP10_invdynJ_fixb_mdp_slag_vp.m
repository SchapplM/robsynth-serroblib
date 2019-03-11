% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP10_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:50
% EndTime: 2019-03-09 06:33:03
% DurationCPUTime: 9.28s
% Computational Cost: add. (6696->622), mult. (12851->763), div. (0->0), fcn. (8393->10), ass. (0->259)
t554 = sin(qJ(4));
t559 = cos(qJ(3));
t665 = qJD(1) * t559;
t624 = t554 * t665;
t558 = cos(qJ(4));
t662 = qJD(3) * t558;
t497 = -t624 + t662;
t664 = qJD(3) * t554;
t498 = t558 * t665 + t664;
t553 = sin(qJ(5));
t557 = cos(qJ(5));
t439 = -t557 * t497 + t498 * t553;
t596 = t497 * t553 + t557 * t498;
t748 = t439 * t596;
t555 = sin(qJ(3));
t608 = pkin(3) * t559 + pkin(8) * t555;
t502 = t608 * qJD(1);
t486 = t558 * t502;
t562 = -pkin(1) - pkin(7);
t527 = qJD(1) * t562 + qJD(2);
t719 = pkin(9) + pkin(8);
t631 = qJD(4) * t719;
t693 = t555 * t558;
t639 = pkin(9) * t693;
t697 = t554 * t559;
t747 = -t527 * t697 + t486 + (pkin(4) * t559 + t639) * qJD(1) + t558 * t631;
t666 = qJD(1) * t555;
t630 = t554 * t666;
t687 = t558 * t559;
t673 = t554 * t502 + t527 * t687;
t746 = pkin(9) * t630 + t554 * t631 + t673;
t614 = qJD(4) - t666;
t643 = qJDD(1) * t559;
t581 = qJD(3) * t614 + t643;
t647 = qJD(1) * qJD(4);
t621 = t559 * t647;
t601 = -qJDD(3) + t621;
t588 = t601 * t554;
t567 = t581 * t558 - t588;
t593 = qJD(3) * qJD(4) + t643;
t633 = t554 * t593 + t558 * t621;
t591 = qJDD(3) * t558 - t633;
t648 = qJD(1) * qJD(3);
t623 = t555 * t648;
t573 = t554 * t623 + t591;
t653 = qJD(5) * t557;
t654 = qJD(5) * t553;
t389 = -t497 * t653 + t498 * t654 - t553 * t573 - t557 * t567;
t536 = qJD(4) + t666;
t525 = qJD(5) + t536;
t381 = t439 * t525 - t389;
t390 = qJD(5) * t596 + t553 * t567 - t557 * t573;
t622 = t559 * t648;
t644 = qJDD(1) * t555;
t494 = qJDD(4) + t622 + t644;
t490 = qJDD(5) + t494;
t720 = t596 ^ 2;
t745 = t490 * MDP(25) + (t525 * t596 - t390) * MDP(24) + MDP(21) * t748 + (-t439 ^ 2 + t720) * MDP(22) + t381 * MDP(23);
t516 = t719 * t554;
t517 = t719 * t558;
t453 = -t516 * t553 + t517 * t557;
t552 = qJ(4) + qJ(5);
t544 = sin(t552);
t556 = sin(qJ(1));
t560 = cos(qJ(1));
t741 = -g(1) * t556 + g(2) * t560;
t590 = t741 * t559;
t713 = g(3) * t555;
t744 = t453 * t490 - (-t590 - t713) * t544;
t703 = t527 * t559;
t711 = qJD(3) * pkin(3);
t488 = -t703 - t711;
t449 = -pkin(4) * t497 + t488;
t393 = pkin(5) * t439 - qJ(6) * t596 + t449;
t743 = t393 * t439;
t742 = t439 * t449;
t640 = qJD(4) + qJD(5);
t656 = qJD(4) * t558;
t688 = t557 * t558;
t699 = t553 * t554;
t675 = t553 * t630 - t557 * t656 - t558 * t653 + t640 * t699 - t666 * t688;
t501 = t553 * t558 + t554 * t557;
t447 = t640 * t501;
t479 = t501 * qJD(1);
t674 = t555 * t479 + t447;
t663 = qJD(3) * t555;
t629 = t554 * t663;
t655 = qJD(4) * t559;
t739 = t558 * t655 - t629;
t403 = pkin(5) * t596 + qJ(6) * t439;
t595 = -t516 * t557 - t517 * t553;
t735 = -qJD(5) * t595 + t553 * t747 + t557 * t746;
t734 = -qJD(5) * t453 + t553 * t746 - t557 * t747;
t507 = pkin(3) * t555 - pkin(8) * t559 + qJ(2);
t480 = t507 * qJD(1);
t506 = t555 * t527;
t487 = qJD(3) * pkin(8) + t506;
t430 = t554 * t480 + t558 * t487;
t420 = pkin(9) * t497 + t430;
t493 = t558 * t507;
t696 = t554 * t562;
t620 = pkin(4) - t696;
t436 = -pkin(9) * t687 + t555 * t620 + t493;
t685 = t558 * t562;
t524 = t555 * t685;
t672 = t554 * t507 + t524;
t448 = -pkin(9) * t697 + t672;
t733 = t553 * t436 + t557 * t448;
t658 = qJD(4) * t554;
t610 = -t506 + (t630 + t658) * pkin(4);
t472 = t490 * qJ(6);
t508 = t525 * qJD(6);
t732 = t472 + t508;
t564 = qJD(1) ^ 2;
t730 = qJ(2) * t564 - t741;
t686 = t558 * t560;
t694 = t555 * t556;
t481 = -t554 * t694 + t686;
t689 = t556 * t558;
t692 = t555 * t560;
t483 = t554 * t692 + t689;
t729 = -g(1) * t481 - g(2) * t483;
t551 = t559 ^ 2;
t728 = qJDD(1) * t551;
t475 = t490 * pkin(5);
t727 = t475 - qJDD(6);
t545 = cos(t552);
t683 = t560 * t545;
t690 = t556 * t544;
t461 = t555 * t690 - t683;
t701 = t545 * t556;
t463 = t544 * t692 + t701;
t495 = qJD(3) * t608 + qJD(2);
t437 = qJD(1) * t495 + qJDD(1) * t507;
t434 = t558 * t437;
t523 = qJDD(1) * t562 + qJDD(2);
t661 = qJD(3) * t559;
t451 = qJDD(3) * pkin(8) + t523 * t555 + t527 * t661;
t642 = qJDD(3) * t554;
t380 = -t554 * t451 + t434 - (t642 + (-t623 + t643) * t558) * pkin(9) + t494 * pkin(4) - t420 * qJD(4);
t635 = -t554 * t437 - t558 * t451 - t480 * t656;
t585 = -t487 * t658 - t635;
t384 = pkin(9) * t573 + t585;
t429 = t558 * t480 - t487 * t554;
t419 = -pkin(9) * t498 + t429;
t413 = pkin(4) * t536 + t419;
t612 = -t557 * t380 + t553 * t384 + t413 * t654 + t420 * t653;
t702 = t544 * t559;
t577 = g(1) * t461 - g(2) * t463 + g(3) * t702 - t612;
t568 = t393 * t596 - t577 - t727;
t726 = -t596 * t449 + t577;
t474 = t558 * t495;
t401 = t474 + (-t524 + (pkin(9) * t559 - t507) * t554) * qJD(4) + (t559 * t620 + t639) * qJD(3);
t660 = qJD(3) * t562;
t634 = t554 * t495 + t507 * t656 + t660 * t687;
t657 = qJD(4) * t555;
t404 = -pkin(9) * t739 - t657 * t696 + t634;
t723 = -qJD(5) * t733 + t401 * t557 - t404 * t553;
t722 = t555 * t640;
t718 = pkin(8) * t494;
t712 = g(3) * t559;
t710 = pkin(1) * qJDD(1);
t709 = qJDD(3) * pkin(3);
t706 = t420 * t557;
t388 = t413 * t553 + t706;
t708 = t388 * t525;
t707 = t420 * t553;
t704 = t498 * t536;
t700 = t545 * t559;
t698 = t554 * t494;
t695 = t555 * t494;
t684 = t559 * t719;
t682 = pkin(5) * t674 + qJ(6) * t675 - qJD(6) * t501 + t610;
t500 = -t688 + t699;
t469 = t500 * t559;
t681 = -qJD(3) * t469 - t501 * t722 - t479;
t467 = t501 * t559;
t680 = qJD(3) * t467 + (-qJD(1) - t722) * t500;
t679 = -qJ(6) * t665 - t735;
t678 = pkin(5) * t665 - t734;
t392 = t419 * t557 - t707;
t671 = pkin(4) * t653 + qJD(6) - t392;
t670 = t560 * pkin(1) + t556 * qJ(2);
t669 = t555 ^ 2 - t551;
t563 = qJD(3) ^ 2;
t668 = -t563 - t564;
t659 = qJD(4) * t497;
t652 = t488 * qJD(4);
t651 = t498 * qJD(3);
t650 = t498 * qJD(4);
t387 = t413 * t557 - t707;
t649 = qJD(6) - t387;
t645 = qJDD(1) * qJ(2);
t641 = qJDD(3) * t555;
t638 = g(1) * t701;
t636 = 0.2e1 * qJD(1) * qJD(2);
t542 = pkin(4) * t558 + pkin(3);
t632 = pkin(4) * t554 + pkin(7);
t628 = t555 * t651;
t627 = t555 * t662;
t626 = t554 * t655;
t618 = t536 * t562 + t487;
t537 = pkin(4) * t697;
t496 = -t559 * t562 + t537;
t617 = -qJD(4) * t480 - t451;
t616 = qJDD(2) - t710;
t615 = qJD(1) + t657;
t613 = t553 * t380 + t557 * t384 + t413 * t653 - t420 * t654;
t611 = -pkin(5) * t702 + qJ(6) * t700;
t391 = t419 * t553 + t706;
t609 = pkin(4) * t654 - t391;
t607 = g(1) * t463 + g(2) * t461;
t462 = t544 * t560 + t545 * t694;
t464 = t555 * t683 - t690;
t606 = -g(1) * t464 - g(2) * t462;
t605 = g(1) * t560 + g(2) * t556;
t603 = g(2) * t559 * t683 + t490 * t595 + t545 * t713;
t598 = t436 * t557 - t448 * t553;
t594 = t542 * t555 - t684;
t592 = pkin(5) * t545 + qJ(6) * t544 + t542;
t450 = -t559 * t523 + t527 * t663 - t709;
t589 = t536 * t656 + t698;
t454 = pkin(4) * t739 + t555 * t660;
t587 = t601 * t559;
t586 = 0.2e1 * qJ(2) * t648 + qJDD(3) * t562;
t584 = t553 * t401 + t557 * t404 + t436 * t653 - t448 * t654;
t583 = -t523 + t730;
t580 = -qJD(4) * pkin(8) * t536 - t450 + t709 + t713;
t579 = t593 - t623;
t578 = g(1) * t462 - g(2) * t464 + g(3) * t700 - t613;
t575 = -t605 + t636 + 0.2e1 * t645;
t572 = -g(1) * (-t461 * pkin(5) + qJ(6) * t462) - g(2) * (t463 * pkin(5) - qJ(6) * t464);
t570 = t387 * t525 + t578;
t569 = -t562 * t563 + t575;
t405 = -pkin(4) * t573 + t450;
t547 = t560 * qJ(2);
t543 = qJDD(3) * t559;
t541 = -pkin(4) * t557 - pkin(5);
t538 = pkin(4) * t553 + qJ(6);
t484 = -t554 * t556 + t555 * t686;
t482 = t554 * t560 + t555 * t689;
t468 = t500 * t555;
t466 = t501 * t555;
t435 = pkin(5) * t500 - qJ(6) * t501 - t542;
t418 = pkin(5) * t467 + qJ(6) * t469 + t496;
t412 = -t654 * t697 + (t640 * t687 - t629) * t557 + (-t626 - t627) * t553;
t410 = t447 * t559 - t553 * t629 + t557 * t627;
t400 = -pkin(5) * t555 - t598;
t399 = qJ(6) * t555 + t733;
t397 = pkin(4) * t498 + t403;
t386 = qJ(6) * t525 + t388;
t385 = -pkin(5) * t525 + t649;
t377 = pkin(5) * t412 + qJ(6) * t410 + qJD(6) * t469 + t454;
t376 = -pkin(5) * t661 - t723;
t375 = qJ(6) * t661 + qJD(6) * t555 + t584;
t374 = t390 * pkin(5) + t389 * qJ(6) - qJD(6) * t596 + t405;
t373 = t612 - t727;
t372 = t613 + t732;
t1 = [(t536 * t661 + t695) * MDP(18) + ((-t507 * t658 + t474) * t536 + t493 * t494 - g(1) * t484 - g(2) * t482 + (-t497 * t660 + t434 - t618 * t656 + (-qJD(3) * t488 - t494 * t562 + t617) * t554) * t555 + (t562 * t591 + t450 * t554 + t558 * t652 + (t429 + (-t536 + t666) * t696) * qJD(3)) * t559) * MDP(19) + (-t372 * t467 - t373 * t469 - t375 * t439 + t376 * t596 - t385 * t410 - t386 * t412 - t389 * t400 - t390 * t399 + t559 * t605) * MDP(29) + (-t498 * t626 + (-t554 * t587 + t579 * t687 - t628) * t558) * MDP(14) + (-t555 * t563 + t543) * MDP(9) + (-t389 * t555 - t410 * t525 - t469 * t490 + t596 * t661) * MDP(23) + (t372 * t555 + t374 * t469 + t375 * t525 - t377 * t596 + t386 * t661 + t389 * t418 + t393 * t410 + t399 * t490 - t607) * MDP(30) + (t389 * t469 - t410 * t596) * MDP(21) + (t389 * t467 + t390 * t469 + t410 * t439 - t412 * t596) * MDP(22) + t575 * MDP(5) + ((-t497 * t558 + t498 * t554) * t663 + ((t591 - t650) * t558 + (-t659 + t588 + (-t643 + (-qJD(4) + 0.2e1 * t666) * qJD(3)) * t558) * t554) * t559) * MDP(15) + (t387 * t661 + t496 * t390 + t405 * t467 + t449 * t412 + t454 * t439 + t598 * t490 + t525 * t723 - t612 * t555 + t606) * MDP(26) + (-0.2e1 * t555 * t622 + t728) * MDP(7) + (-t390 * t555 - t412 * t525 - t439 * t661 - t467 * t490) * MDP(24) + (-t373 * t555 + t374 * t467 - t376 * t525 + t377 * t439 - t385 * t661 + t390 * t418 + t393 * t412 - t400 * t490 + t606) * MDP(28) - t741 * MDP(2) + (qJDD(2) + t741 - 0.2e1 * t710) * MDP(4) + (-t388 * t661 - t496 * t389 - t405 * t469 - t449 * t410 + t454 * t596 - t490 * t733 - t525 * t584 - t555 * t613 + t607) * MDP(27) + (((t536 + t666) * t664 + t591) * t555 + (qJD(3) * t497 - t589) * t559) * MDP(17) + 0.2e1 * (-t555 * t643 + t648 * t669) * MDP(8) + (t372 * t399 + t386 * t375 + t374 * t418 + t393 * t377 + t373 * t400 + t385 * t376 - g(1) * (pkin(5) * t464 + qJ(6) * t463 + t547) - g(2) * (pkin(5) * t462 + qJ(6) * t461 + t670) + (-g(1) * t594 - g(2) * t632) * t560 + (-g(1) * (-pkin(1) - t632) - g(2) * t594) * t556) * MDP(31) + (-t616 * pkin(1) - g(1) * (-pkin(1) * t556 + t547) - g(2) * t670 + (t636 + t645) * qJ(2)) * MDP(6) + (-t634 * t536 - t672 * t494 + g(1) * t483 - g(2) * t481 + (t618 * t658 + (-t488 * t558 + t498 * t562) * qJD(3) + t635) * t555 + (-t554 * t652 + t450 * t558 + (-t642 - (qJDD(1) * t558 - t554 * t647) * t559) * t562 + (-t614 * t685 - t430) * qJD(3)) * t559) * MDP(20) + (-t559 * t563 - t641) * MDP(10) + qJDD(1) * MDP(1) + (t490 * t555 + t525 * t661) * MDP(25) + (t559 * t651 + (-t536 * t655 - t555 * t601) * t554 + ((t494 + t644) * t559 + (-t536 + t614) * t663) * t558) * MDP(16) + t605 * MDP(3) + (t555 * t569 + t559 * t586) * MDP(12) + (-t555 * t586 + t559 * t569) * MDP(13); qJDD(1) * MDP(4) - t564 * MDP(5) + (t616 - t730) * MDP(6) + (t555 * t668 + t543) * MDP(12) + (t559 * t668 - t641) * MDP(13) + (t559 * t591 + (-t698 + (-t497 + t624) * qJD(3)) * t555 + (-t554 * t661 - t558 * t615) * t536) * MDP(19) + (t628 + (t536 * t615 + t587) * t554 + (-t728 - t695 + (-t536 - t614) * t661) * t558) * MDP(20) + (-t389 * t466 + t390 * t468 - t439 * t681 + t596 * t680) * MDP(29) + (-t372 * t468 + t373 * t466 - t374 * t559 + t385 * t680 + t386 * t681 + t393 * t663 + t741) * MDP(31) + (MDP(26) + MDP(28)) * (-t390 * t559 + t439 * t663 - t466 * t490 - t525 * t680) + (-MDP(27) + MDP(30)) * (-t389 * t559 - t468 * t490 + t525 * t681 - t596 * t663); (-pkin(3) * t633 - t486 * t536 + t497 * t506 + (t536 * t703 - t718 + t652 + (t488 + t711) * t666) * t554 + (t590 + t580) * t558) * MDP(19) + (t673 * t536 - t498 * t506 + (-pkin(3) * t579 + t488 * t536 - t718) * t558 + ((pkin(3) * t647 - t741) * t559 - t580) * t554) * MDP(20) + (-t372 * t500 + t373 * t501 - t385 * t675 - t386 * t674 + t389 * t595 - t390 * t453 - t439 * t679 + t555 * t741 + t596 * t678 - t712) * MDP(29) + (t555 * t583 + t712) * MDP(13) + (t542 * t389 + t405 * t501 - t675 * t449 + t735 * t525 + t610 * t596 - t744) * MDP(27) + (-t374 * t501 + t389 * t435 + t675 * t393 + t679 * t525 - t682 * t596 + t744) * MDP(30) + (-t559 * t583 + t713) * MDP(12) + (-g(3) * t684 + t372 * t453 - t373 * t595 + t374 * t435 + t678 * t385 + t679 * t386 + t682 * t393 + t592 * t713 + t741 * (t555 * t719 + t559 * t592)) * MDP(31) + (-t601 * t554 ^ 2 + (t554 * t581 + t704) * t558) * MDP(14) + ((-t498 * t559 + t536 * t693) * qJD(1) + t589) * MDP(16) + (-t536 * t658 + t558 * t494 + (-t536 * t554 * t555 - t497 * t559) * qJD(1)) * MDP(17) - MDP(10) * t644 + qJDD(3) * MDP(11) + ((-t650 + (-t498 + t664) * t666 - t633) * t554 + (t659 + 0.2e1 * t642 + t593 * t558 + (-t626 + (t497 - t662) * t555) * qJD(1)) * t558) * MDP(15) + (-t542 * t390 + t405 * t500 + (-qJD(1) * t387 - t638) * t559 + t734 * t525 + t674 * t449 + t610 * t439 + t603) * MDP(26) + (-t490 * t500 - t525 * t674) * MDP(24) + (t389 * t500 - t390 * t501 + t439 * t675 - t596 * t674) * MDP(22) + (t490 * t501 - t525 * t675) * MDP(23) + (-t389 * t501 - t596 * t675) * MDP(21) + (t374 * t500 + t390 * t435 + (qJD(1) * t385 - t638) * t559 - t678 * t525 + t682 * t439 + t674 * t393 + t603) * MDP(28) + MDP(9) * t643 + (-t536 * MDP(18) - t429 * MDP(19) + t430 * MDP(20) - MDP(23) * t596 + t439 * MDP(24) - t525 * MDP(25) + t388 * MDP(27) - t386 * MDP(30)) * t665 + (MDP(7) * t555 * t559 - MDP(8) * t669) * t564; -t498 * t497 * MDP(14) + (-t497 ^ 2 + t498 ^ 2) * MDP(15) + (-t497 * t536 + t567) * MDP(16) + (t573 + t704) * MDP(17) + t494 * MDP(18) + (-t487 * t656 + t430 * t536 - t488 * t498 + t434 + (t617 + t712) * t554 + t729) * MDP(19) + (g(1) * t482 - g(2) * t484 + g(3) * t687 + t429 * t536 - t488 * t497 - t585) * MDP(20) + (t391 * t525 + (-t439 * t498 + t557 * t490 - t525 * t654) * pkin(4) + t726) * MDP(26) + (t392 * t525 + t742 + (-t553 * t490 - t498 * t596 - t525 * t653) * pkin(4) + t578) * MDP(27) + (-t397 * t439 - t490 * t541 - t525 * t609 - t568) * MDP(28) + (-t389 * t541 - t390 * t538 + (t386 + t609) * t596 + (t385 - t671) * t439) * MDP(29) + (t397 * t596 + t490 * t538 + t525 * t671 - t578 + t732 - t743) * MDP(30) + (t372 * t538 + t373 * t541 - t393 * t397 - t385 * t391 - g(3) * (-t537 + t611) + t671 * t386 + (t385 * t654 + t729) * pkin(4) + t572) * MDP(31) + t745; (t708 + t726) * MDP(26) + (t570 + t742) * MDP(27) + (-t403 * t439 + t475 - t568 + t708) * MDP(28) + (pkin(5) * t389 - qJ(6) * t390 + (t386 - t388) * t596 + (t385 - t649) * t439) * MDP(29) + (t403 * t596 + 0.2e1 * t472 + 0.2e1 * t508 - t570 - t743) * MDP(30) + (-t373 * pkin(5) - g(3) * t611 + t372 * qJ(6) - t385 * t388 + t386 * t649 - t393 * t403 + t572) * MDP(31) + t745; (-t490 + t748) * MDP(28) + t381 * MDP(29) + (-t525 ^ 2 - t720) * MDP(30) + (-t386 * t525 + t568) * MDP(31);];
tau  = t1;
