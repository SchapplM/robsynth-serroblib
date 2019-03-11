% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:37:03
% EndTime: 2019-03-08 23:37:18
% DurationCPUTime: 11.30s
% Computational Cost: add. (4410->631), mult. (9987->832), div. (0->0), fcn. (7828->12), ass. (0->259)
t588 = sin(qJ(4));
t592 = cos(qJ(4));
t589 = sin(qJ(3));
t682 = qJDD(2) * t589;
t593 = cos(qJ(3));
t701 = qJD(2) * t593;
t695 = qJD(4) * t589;
t779 = qJD(2) * t695 - qJDD(3);
t457 = t588 * (qJD(3) * (qJD(4) + t701) + t682) + t779 * t592;
t688 = t592 * qJD(3);
t683 = qJD(2) * qJD(3);
t664 = t593 * t683;
t782 = -t664 - t682;
t456 = -qJD(4) * t688 + t588 * t779 + t592 * t782;
t699 = qJD(3) * t588;
t702 = qJD(2) * t589;
t533 = t592 * t702 + t699;
t586 = sin(pkin(6));
t590 = sin(qJ(2));
t594 = cos(qJ(2));
t684 = qJD(1) * qJD(2);
t500 = qJDD(2) * pkin(8) + (qJDD(1) * t590 + t594 * t684) * t586;
t704 = qJD(1) * t586;
t541 = qJD(2) * pkin(8) + t590 * t704;
t743 = cos(pkin(6));
t660 = qJD(1) * t743;
t565 = t589 * t660;
t655 = qJDD(1) * t743;
t697 = qJD(3) * t593;
t652 = -qJD(3) * t565 - t589 * t500 - t541 * t697 + t593 * t655;
t618 = qJDD(3) * pkin(3) + t652;
t601 = -qJ(5) * t456 + qJD(5) * t533 + t618;
t757 = pkin(4) + pkin(5);
t398 = -t457 * t757 + t601;
t585 = sin(pkin(11));
t742 = cos(pkin(11));
t642 = t743 * t742;
t513 = t585 * t594 + t590 * t642;
t662 = t585 * t743;
t515 = -t590 * t662 + t594 * t742;
t730 = t586 * t590;
t518 = t589 * t730 - t593 * t743;
t661 = t586 * t742;
t729 = t586 * t593;
t766 = g(3) * t518 - g(2) * (-t513 * t589 - t593 * t661) - g(1) * (-t515 * t589 + t585 * t729);
t786 = t398 + t766;
t531 = t588 * t702 - t688;
t587 = sin(qJ(6));
t591 = cos(qJ(6));
t463 = t531 * t587 + t533 * t591;
t580 = t593 * qJDD(2);
t771 = -t589 * t683 + t580;
t529 = qJDD(4) - t771;
t523 = -qJDD(6) + t529;
t630 = -t591 * t531 + t533 * t587;
t785 = t630 * MDP(23) * t463 + (t463 ^ 2 - t630 ^ 2) * MDP(24) - t523 * MDP(27);
t781 = qJD(4) - t701;
t685 = -qJD(6) + t781;
t784 = t463 * t685;
t719 = t593 * t594;
t498 = (t588 * t590 + t592 * t719) * t586;
t648 = pkin(3) * t589 - pkin(9) * t593;
t540 = t648 * qJD(3);
t627 = pkin(3) * t593 + pkin(9) * t589 + pkin(2);
t694 = qJD(4) * t592;
t783 = -qJD(1) * t498 + t588 * t540 - t627 * t694;
t780 = qJD(6) - qJD(4);
t772 = -t589 * t541 + t593 * t660;
t649 = qJD(3) * pkin(3) + t772;
t610 = qJ(5) * t533 + t649;
t414 = -t531 * t757 + t610;
t519 = t589 * t743 + t590 * t729;
t728 = t586 * t594;
t477 = t519 * t588 + t592 * t728;
t478 = t519 * t592 - t588 * t728;
t416 = t477 * t591 - t478 * t587;
t470 = t513 * t593 - t589 * t661;
t512 = t585 * t590 - t594 * t642;
t434 = t470 * t588 - t512 * t592;
t435 = t470 * t592 + t512 * t588;
t472 = t585 * t586 * t589 + t515 * t593;
t514 = t590 * t742 + t594 * t662;
t436 = t472 * t588 - t514 * t592;
t437 = t472 * t592 + t514 * t588;
t645 = t589 * t655;
t423 = qJDD(3) * pkin(9) + qJD(3) * t772 + t593 * t500 + t645;
t666 = t590 * t684;
t639 = -qJDD(1) * t728 + t586 * t666;
t445 = qJD(2) * t540 - qJDD(2) * t627 + t639;
t494 = t593 * t541 + t565;
t484 = qJD(3) * pkin(9) + t494;
t675 = t594 * t704;
t495 = -qJD(2) * t627 - t675;
t696 = qJD(4) * t588;
t653 = t588 * t423 - t592 * t445 + t484 * t694 + t495 * t696;
t628 = qJDD(5) + t653;
t394 = pkin(10) * t456 - t529 * t757 + t628;
t522 = t529 * qJ(5);
t558 = t781 * qJD(5);
t616 = t592 * t423 + t588 * t445 - t484 * t696 + t495 * t694;
t397 = t522 + t558 + t616;
t396 = pkin(10) * t457 + t397;
t659 = t591 * t394 - t587 * t396;
t777 = t414 * t463 + g(1) * (t436 * t591 - t437 * t587) + g(2) * (t434 * t591 - t435 * t587) - t659 + g(3) * t416;
t775 = t685 * t630;
t720 = t592 * t593;
t573 = pkin(8) * t720;
t768 = t588 * t719 - t590 * t592;
t774 = qJD(4) * t573 - t540 * t592 - t627 * t696 - t768 * t704;
t773 = qJD(5) * t588 + t494;
t667 = t588 * t695;
t770 = t593 * t688 - t667;
t698 = qJD(3) * t589;
t769 = qJ(5) * t698 + t783;
t430 = -t588 * t484 + t592 * t495;
t686 = qJD(5) - t430;
t432 = pkin(4) * t531 - t610;
t754 = pkin(9) * t529;
t765 = -t432 * t781 + t754;
t596 = qJD(3) ^ 2;
t749 = g(2) * t512;
t751 = g(1) * t514;
t647 = t749 + t751;
t764 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t596 + t586 * (-g(3) * t594 + t666) - t639 + t647;
t740 = qJ(5) * t588;
t763 = -t592 * t757 - t740;
t417 = t477 * t587 + t478 * t591;
t687 = pkin(10) * t533 - t686;
t406 = -t757 * t781 - t687;
t690 = qJD(6) * t591;
t679 = t587 * t394 + t591 * t396 + t406 * t690;
t759 = -g(1) * (t436 * t587 + t437 * t591) - g(2) * (t434 * t587 + t435 * t591) - t414 * t630 + t679 - g(3) * t417;
t758 = t533 ^ 2;
t756 = pkin(9) - pkin(10);
t755 = pkin(4) * t529;
t753 = pkin(10) * t589;
t745 = pkin(9) * qJD(4);
t744 = qJD(2) * pkin(2);
t741 = qJ(5) * t531;
t739 = qJ(5) * t592;
t431 = t592 * t484 + t588 * t495;
t560 = t781 * qJ(5);
t418 = t560 + t431;
t737 = t418 * t781;
t736 = t431 * t781;
t735 = t456 * t588;
t734 = t531 * t533;
t733 = t531 * t781;
t732 = t533 * t781;
t731 = t533 * t592;
t413 = pkin(10) * t531 + t431;
t407 = t413 + t560;
t727 = t587 * t407;
t726 = t587 * t592;
t725 = t588 * t591;
t724 = t588 * t593;
t723 = t589 * t591;
t722 = t589 * t592;
t718 = qJDD(1) - g(3);
t676 = -pkin(8) * t588 - pkin(4);
t681 = pkin(10) * t720;
t717 = pkin(10) * t667 + (-t681 + (-pkin(5) + t676) * t589) * qJD(3) + t774;
t716 = -(-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t722 - (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t588) * t593 - t769;
t715 = -qJD(5) * t593 + (-t589 * t688 - t593 * t696) * pkin(8) + t769;
t714 = t676 * t698 + t774;
t535 = t587 * t588 + t591 * t592;
t617 = t593 * t535;
t713 = -qJD(2) * t617 - t535 * t780;
t672 = t588 * t701;
t691 = qJD(6) * t587;
t712 = t587 * t694 + t588 * t690 - t592 * t691 - t701 * t726 + (t672 - t696) * t591;
t537 = t648 * qJD(2);
t711 = t588 * t537 + t592 * t772;
t622 = -t588 * t757 + t739;
t710 = t622 * t781 + t773;
t640 = pkin(4) * t588 - t739;
t709 = t640 * t781 - t773;
t707 = -t588 * t627 + t573;
t583 = t589 ^ 2;
t706 = -t593 ^ 2 + t583;
t703 = qJD(2) * t586;
t700 = qJD(3) * t533;
t692 = qJD(5) * t592;
t553 = t756 * t592;
t678 = -t591 * t456 + t587 * t457 + t531 * t690;
t441 = qJ(5) * t702 + t711;
t674 = t590 * t703;
t673 = t594 * t703;
t671 = t781 * t699;
t670 = t781 * t688;
t668 = t781 * t696;
t658 = -t456 * t587 - t591 * t457;
t481 = t588 * t772;
t657 = -t537 * t592 + t481;
t572 = pkin(8) * t724;
t656 = -t592 * t627 - t572;
t651 = t531 * t675;
t650 = t533 * t675;
t488 = -qJ(5) * t593 + t707;
t644 = (-t589 * t757 - t681) * qJD(2) + t657 + t780 * t553;
t552 = t756 * t588;
t643 = pkin(10) * t672 - qJD(6) * t552 + t756 * t696 + t441;
t641 = pkin(4) * t592 + t740;
t638 = qJ(5) * t591 - t587 * t757;
t637 = qJ(5) * t587 + t591 * t757;
t401 = t587 * t406 + t591 * t407;
t415 = -pkin(4) * t781 + t686;
t636 = t415 * t592 - t418 * t588;
t582 = t593 * pkin(4);
t458 = pkin(5) * t593 + t572 + t582 + (t627 - t753) * t592;
t464 = t588 * t753 + t488;
t631 = t458 * t587 + t464 * t591;
t629 = -t725 + t726;
t626 = pkin(3) + t641;
t625 = pkin(8) + t640;
t615 = -pkin(8) + t622;
t624 = t615 * t697 + (qJD(4) * t763 + t675 + t692) * t589;
t620 = t529 * t588 + t694 * t781;
t619 = t529 * t592 - t668;
t509 = t535 * t589;
t403 = -t533 * t691 + t678;
t447 = -t512 * t724 - t513 * t592;
t449 = -t514 * t724 - t515 * t592;
t497 = t768 * t586;
t614 = g(1) * t449 + g(2) * t447 + g(3) * t497;
t448 = -t512 * t720 + t513 * t588;
t450 = -t514 * t720 + t515 * t588;
t613 = -g(1) * t450 - g(2) * t448 - g(3) * t498;
t611 = g(1) * t472 + g(2) * t470 + g(3) * t519;
t607 = -t649 * t781 - t754;
t605 = -t745 * t781 + t766;
t603 = g(1) * t436 + g(2) * t434 + g(3) * t477 - t653;
t404 = qJD(6) * t463 + t658;
t402 = pkin(4) * t457 - t601;
t602 = -t402 + t605;
t542 = -t675 - t744;
t600 = -pkin(8) * qJDD(3) + (t542 + t675 - t744) * qJD(3);
t599 = t432 * t533 + qJDD(5) - t603;
t598 = g(1) * t437 + g(2) * t435 + g(3) * t478 + t430 * t781 - t616;
t597 = qJD(2) ^ 2;
t528 = pkin(3) - t763;
t508 = t587 * t722 - t588 * t723;
t507 = t625 * t589;
t489 = t582 - t656;
t487 = t615 * t589;
t476 = qJD(3) * t519 + t589 * t673;
t475 = -qJD(3) * t518 + t593 * t673;
t468 = pkin(4) * t533 + t741;
t446 = -t533 * t757 - t741;
t444 = (qJD(4) * t641 - t692) * t589 + t625 * t697;
t443 = -pkin(4) * t702 + t657;
t428 = -t456 + t733;
t426 = qJD(6) * t509 + t587 * t770 - t694 * t723 - t697 * t725;
t425 = -t589 * t629 * t780 + qJD(3) * t617;
t410 = -qJD(4) * t477 + t475 * t592 + t588 * t674;
t409 = qJD(4) * t478 + t475 * t588 - t592 * t674;
t400 = t406 * t591 - t727;
t399 = t628 - t755;
t1 = [t718 * MDP(1) + (-qJD(3) * t476 - qJDD(3) * t518) * MDP(10) + (-qJD(3) * t475 - qJDD(3) * t519) * MDP(11) + (t409 * t533 - t410 * t531 - t456 * t477 - t457 * t478) * MDP(20) + (t397 * t478 + t399 * t477 + t402 * t518 + t409 * t415 + t410 * t418 + t432 * t476 - g(3)) * MDP(22) + (-(-qJD(6) * t417 + t409 * t591 - t410 * t587) * t685 - t416 * t523 - t476 * t630 - t518 * t404) * MDP(28) + ((qJD(6) * t416 + t409 * t587 + t410 * t591) * t685 + t417 * t523 - t476 * t463 - t518 * t403) * MDP(29) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t593 + MDP(11) * t589 - MDP(3)) * t597) * t590 + (t771 * MDP(10) + MDP(11) * t782 + qJDD(2) * MDP(3) - t597 * MDP(4)) * t594) * t586 + (MDP(17) + MDP(19)) * (-t409 * t781 + t518 * t457 + t476 * t531 - t477 * t529) + (MDP(18) - MDP(21)) * (-t410 * t781 - t456 * t518 + t476 * t533 - t478 * t529); (t403 * t509 + t425 * t463) * MDP(23) + 0.2e1 * (t580 * t589 - t683 * t706) * MDP(6) + (t631 * t523 - (-t407 * t691 + t679) * t593 + t401 * t698 + t487 * t403 + t398 * t509 + t414 * t425 - g(1) * (t449 * t591 - t450 * t587) - g(2) * (t447 * t591 - t448 * t587) - g(3) * (t497 * t591 - t498 * t587) - ((-qJD(6) * t458 + t716) * t591 + (qJD(6) * t464 - t717) * t587) * t685 + t624 * t463) * MDP(29) + (-(t458 * t591 - t464 * t587) * t523 + t659 * t593 - t400 * t698 + t487 * t404 + t398 * t508 + t414 * t426 - g(1) * (t449 * t587 + t450 * t591) - g(2) * (t447 * t587 + t448 * t591) - g(3) * (t497 * t587 + t498 * t591) - (t587 * t716 + t591 * t717) * t685 + t624 * t630 + (-t401 * t593 + t631 * t685) * qJD(6)) * MDP(28) + (-t404 * t593 + t426 * t685 + t508 * t523 + t630 * t698) * MDP(26) + (-t523 * t593 + t685 * t698) * MDP(27) + ((-t531 * t592 - t533 * t588) * t697 + (t735 - t457 * t592 + (t531 * t588 - t731) * qJD(4)) * t589) * MDP(13) + (-t456 * t489 - t457 * t488 + t714 * t533 - t715 * t531 + t636 * t697 + (-g(3) * t728 - t397 * t588 + t399 * t592 + (-t415 * t588 - t418 * t592) * qJD(4) + t647) * t589) * MDP(20) + (t718 * t728 + t647) * MDP(3) + (g(1) * t515 + g(2) * t513 - t718 * t730) * MDP(4) + ((t457 - t671) * t593 + (-qJD(3) * t531 - t620) * t589) * MDP(15) + ((t456 + t670) * t593 + (t619 + t700) * t589) * MDP(14) + (-t403 * t508 - t404 * t509 - t425 * t630 - t426 * t463) * MDP(24) + (qJDD(3) * t589 + t593 * t596) * MDP(7) + (qJDD(3) * t593 - t589 * t596) * MDP(8) + (t656 * t529 - t774 * t781 + ((pkin(8) * t531 - t588 * t649) * qJD(3) + t653) * t593 + (-t651 - t649 * t694 + t430 * qJD(3) - t618 * t588 + (t457 + t671) * pkin(8)) * t589 + t613) * MDP(17) + (t444 * t531 + t457 * t507 - t489 * t529 + (t432 * t699 + t399) * t593 - t714 * t781 + (-qJD(3) * t415 + t402 * t588 + t432 * t694 - t651) * t589 + t613) * MDP(19) + (-t444 * t533 + t456 * t507 + t488 * t529 + (-t432 * t688 - t397) * t593 + t715 * t781 + (qJD(3) * t418 - t402 * t592 + t432 * t696 + t650) * t589 - t614) * MDP(21) + (-t529 * t593 + t698 * t781) * MDP(16) + (-t707 * t529 - t783 * t781 + (-t649 * t688 + (t668 + t700) * pkin(8) + t616) * t593 + (-t650 + t649 * t696 - t431 * qJD(3) - t618 * t592 + (-t456 + t670) * pkin(8)) * t589 + t614) * MDP(18) + (t397 * t488 + t402 * t507 + t432 * t444 + t399 * t489 - g(1) * (pkin(4) * t450 + pkin(8) * t515 + qJ(5) * t449) - g(2) * (pkin(4) * t448 + pkin(8) * t513 + qJ(5) * t447) - g(3) * (pkin(4) * t498 + qJ(5) * t497) + t627 * t751 + t627 * t749 + t715 * t418 + t714 * t415 + (-g(3) * pkin(8) * t590 + (-t432 * t589 * qJD(1) - g(3) * t627) * t594) * t586) * MDP(22) + qJDD(2) * MDP(2) + (qJDD(2) * t583 + 0.2e1 * t589 * t664) * MDP(5) + (t600 * t589 + t593 * t764) * MDP(10) + (-t589 * t764 + t600 * t593) * MDP(11) + (-t456 * t722 + t533 * t770) * MDP(12) + (t403 * t593 - t425 * t685 - t463 * t698 - t509 * t523) * MDP(25); (-t403 * t535 + t404 * t629 - t463 * t712 - t630 * t713) * MDP(24) + (t523 * t629 - t685 * t713) * MDP(25) + (t523 * t535 + t685 * t712) * MDP(26) + (qJD(3) * t494 + t652 + t766) * MDP(10) + (t441 * t531 - t443 * t533 + (t397 + t781 * t415 + (qJD(4) * t533 - t457) * pkin(9)) * t592 + (t399 - t737 + (qJD(4) * t531 - t456) * pkin(9)) * t588 - t611) * MDP(20) + (t731 * t781 - t735) * MDP(12) + ((-t456 - t733) * t592 + (-t457 - t732) * t588) * MDP(13) + ((t531 * t589 + t724 * t781) * qJD(2) + t619) * MDP(15) + ((-t533 * t589 - t720 * t781) * qJD(2) + t620) * MDP(14) + (pkin(3) * t456 + t711 * t781 - t494 * t533 + t607 * t592 + (-t618 - t605) * t588) * MDP(18) + (-pkin(3) * t457 + t481 * t781 - t494 * t531 + t607 * t588 + (t618 - (t537 + t745) * t781 + t766) * t592) * MDP(17) + (-t403 * t629 + t463 * t713) * MDP(23) + qJDD(3) * MDP(9) + (-t645 + (-qJD(2) * t542 - t500) * t593 + t611) * MDP(11) + (-t441 * t781 - t456 * t626 - t709 * t533 + t602 * t588 + t592 * t765) * MDP(21) + (t443 * t781 - t457 * t626 + t709 * t531 - t588 * t765 + t602 * t592) * MDP(19) + (-(t552 * t591 - t553 * t587) * t523 + t528 * t404 - (t587 * t643 - t591 * t644) * t685 + t710 * t630 + t712 * t414 + t786 * t535) * MDP(28) + (-t415 * t443 - t418 * t441 + t709 * t432 + (qJD(4) * t636 + t397 * t592 + t399 * t588 - t611) * pkin(9) + (-t402 + t766) * t626) * MDP(22) + ((t552 * t587 + t553 * t591) * t523 + t528 * t403 - (t587 * t644 + t591 * t643) * t685 + t710 * t463 + t713 * t414 - t786 * t629) * MDP(29) + MDP(8) * t580 + MDP(7) * t682 + (-MDP(10) * t542 - MDP(16) * t781 - MDP(17) * t430 + MDP(18) * t431 + t415 * MDP(19) - t418 * MDP(21) + t463 * MDP(25) - MDP(26) * t630 - MDP(27) * t685 + t400 * MDP(28) - t401 * MDP(29)) * t702 + (-MDP(5) * t589 * t593 + MDP(6) * t706) * t597; MDP(12) * t734 + (-t531 ^ 2 + t758) * MDP(13) + t428 * MDP(14) + (-t457 + t732) * MDP(15) + t529 * MDP(16) + (t533 * t649 + t603 + t736) * MDP(17) + (-t531 * t649 + t598) * MDP(18) + (-t468 * t531 - t599 + t736 + 0.2e1 * t755) * MDP(19) + (pkin(4) * t456 - qJ(5) * t457 + (t418 - t431) * t533 + (t415 - t686) * t531) * MDP(20) + (-t432 * t531 + t468 * t533 + 0.2e1 * t522 + 0.2e1 * t558 - t598) * MDP(21) + (t397 * qJ(5) - t399 * pkin(4) - t432 * t468 - t415 * t431 - g(1) * (-pkin(4) * t436 + qJ(5) * t437) - g(2) * (-pkin(4) * t434 + qJ(5) * t435) - g(3) * (-pkin(4) * t477 + qJ(5) * t478) + t686 * t418) * MDP(22) + (-t403 + t775) * MDP(25) + (t404 + t784) * MDP(26) + (t637 * t523 - t446 * t630 - (-t413 * t591 + t587 * t687) * t685 + (t638 * t685 + t401) * qJD(6) + t777) * MDP(28) + (t638 * t523 - t446 * t463 - (t413 * t587 + t591 * t687) * t685 + (-t637 * t685 - t727) * qJD(6) + t759) * MDP(29) - t785; (-t529 + t734) * MDP(19) + t428 * MDP(20) + (-t781 ^ 2 - t758) * MDP(21) + (t599 - t737 - t755) * MDP(22) + (-t523 * t591 - t533 * t630) * MDP(28) + (-t463 * t533 + t523 * t587) * MDP(29) - (MDP(28) * t587 + MDP(29) * t591) * t685 ^ 2; (t678 - t775) * MDP(25) + (-t658 - t784) * MDP(26) + (-t401 * t685 - t777) * MDP(28) + (-t400 * t685 - t759) * MDP(29) + ((-MDP(26) * t533 - MDP(28) * t407) * t591 + (-MDP(25) * t533 - MDP(26) * t531 - MDP(28) * t406 + MDP(29) * t407) * t587) * qJD(6) + t785;];
tau  = t1;
