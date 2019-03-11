% Calculate vector of inverse dynamics joint torques for
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:57
% EndTime: 2019-03-08 22:53:08
% DurationCPUTime: 8.96s
% Computational Cost: add. (4550->612), mult. (10231->744), div. (0->0), fcn. (7781->10), ass. (0->254)
t777 = MDP(19) + MDP(23);
t776 = MDP(20) - MDP(25);
t591 = cos(qJ(3));
t775 = -pkin(3) * t591 - pkin(2);
t587 = sin(qJ(4));
t739 = qJ(5) * t587;
t774 = pkin(3) + t739;
t586 = sin(pkin(6));
t589 = sin(qJ(2));
t590 = cos(qJ(4));
t592 = cos(qJ(2));
t714 = t591 * t592;
t491 = (t587 * t589 + t590 * t714) * t586;
t588 = sin(qJ(3));
t639 = pkin(3) * t588 - pkin(9) * t591;
t540 = t639 * qJD(3);
t544 = -pkin(9) * t588 + t775;
t687 = qJD(4) * t590;
t773 = -qJD(1) * t491 + t587 * t540 + t544 * t687;
t579 = t591 * qJDD(2);
t678 = qJD(2) * qJD(3);
t759 = -t588 * t678 + t579;
t533 = qJDD(4) - t759;
t694 = qJD(2) * t591;
t768 = qJD(4) - t694;
t772 = t533 * qJ(5) + qJD(5) * t768;
t658 = t591 * t678;
t676 = qJDD(2) * t588;
t622 = -t658 - t676;
t771 = -qJD(3) * qJD(4) + t622;
t715 = t590 * t591;
t572 = pkin(8) * t715;
t689 = qJD(4) * t587;
t697 = qJD(1) * t586;
t758 = t587 * t714 - t589 * t590;
t770 = qJD(4) * t572 + t544 * t689 - t758 * t697;
t767 = MDP(21) + MDP(24);
t688 = qJD(4) * t588;
t657 = qJD(2) * t688;
t627 = (-qJDD(3) + t657) * t590;
t449 = t587 * (qJD(3) * (qJD(4) + t694) + t676) + t627;
t683 = t590 * qJD(3);
t695 = qJD(2) * t588;
t535 = t587 * t695 - t683;
t766 = t449 * qJ(6) + t535 * qJD(6);
t765 = -t540 * t590 + t770;
t541 = qJD(2) * pkin(8) + t589 * t697;
t743 = cos(pkin(6));
t654 = qJD(1) * t743;
t563 = t588 * t654;
t486 = t591 * t541 + t563;
t764 = qJD(5) * t587 + t486;
t763 = 0.2e1 * t772;
t762 = -t588 * t541 + t591 * t654;
t691 = qJD(3) * t588;
t760 = qJ(5) * t691 - qJD(5) * t591 + t773;
t476 = qJD(3) * pkin(9) + t486;
t668 = t592 * t697;
t488 = qJD(2) * t544 - t668;
t416 = -t587 * t476 + t590 * t488;
t681 = qJD(5) - t416;
t524 = t533 * pkin(4);
t757 = t524 - qJDD(5);
t692 = qJD(3) * t587;
t537 = t590 * t695 + t692;
t640 = qJD(3) * pkin(3) + t762;
t615 = qJ(5) * t537 + t640;
t418 = pkin(4) * t535 - t615;
t747 = pkin(9) * t533;
t755 = -t418 * t768 + t747;
t750 = pkin(4) + pkin(5);
t402 = -t535 * t750 + qJD(6) + t615;
t585 = sin(pkin(10));
t742 = cos(pkin(10));
t635 = t743 * t742;
t514 = t585 * t592 + t589 * t635;
t655 = t586 * t742;
t462 = t514 * t591 - t588 * t655;
t513 = t585 * t589 - t592 * t635;
t422 = t462 * t587 - t513 * t590;
t656 = t585 * t743;
t516 = -t589 * t656 + t592 * t742;
t464 = t585 * t586 * t588 + t516 * t591;
t515 = t589 * t742 + t592 * t656;
t424 = t464 * t587 - t515 * t590;
t720 = t586 * t591;
t519 = t588 * t743 + t589 * t720;
t719 = t586 * t592;
t469 = t519 * t587 + t590 * t719;
t679 = qJD(1) * qJD(2);
t493 = qJDD(2) * pkin(8) + (qJDD(1) * t589 + t592 * t679) * t586;
t648 = qJDD(1) * t743;
t636 = t588 * t648;
t411 = qJDD(3) * pkin(9) + qJD(3) * t762 + t591 * t493 + t636;
t660 = t589 * t679;
t554 = t586 * t660;
t632 = -qJDD(1) * t719 + t554;
t433 = qJD(2) * t540 + qJDD(2) * t544 + t632;
t646 = t587 * t411 - t590 * t433 + t476 * t687 + t488 * t689;
t604 = g(1) * t424 + g(2) * t422 + g(3) * t469 - t646;
t600 = t604 + t757;
t645 = t587 * qJDD(3) - t590 * t771;
t448 = t587 * t657 - t645;
t737 = qJ(6) * t448;
t754 = (qJD(6) + t402) * t537 + t600 - t737;
t753 = -t590 * t750 - t739;
t751 = t535 ^ 2;
t532 = t537 ^ 2;
t748 = pkin(5) * t533;
t746 = pkin(9) - qJ(6);
t745 = pkin(9) * qJD(4);
t744 = qJD(2) * pkin(2);
t741 = qJ(5) * t449;
t740 = qJ(5) * t535;
t738 = qJ(5) * t590;
t736 = qJ(6) * t588;
t417 = t590 * t476 + t587 * t488;
t404 = qJ(6) * t535 + t417;
t560 = t768 * qJ(5);
t399 = t404 + t560;
t734 = t399 * t768;
t406 = t560 + t417;
t733 = t406 * t768;
t732 = t417 * t768;
t731 = t448 * t587;
t461 = -t514 * t588 - t591 * t655;
t730 = t461 * t590;
t463 = -t516 * t588 + t585 * t720;
t729 = t463 * t590;
t728 = t513 * t588;
t727 = t515 * t588;
t721 = t586 * t589;
t518 = t588 * t721 - t591 * t743;
t726 = t518 * t590;
t725 = t535 * t537;
t724 = t535 * t768;
t723 = t537 * t768;
t722 = t537 * t590;
t718 = t587 * t591;
t717 = t588 * t590;
t713 = qJDD(1) - g(3);
t669 = -pkin(8) * t587 - pkin(4);
t684 = qJD(6) * t590;
t690 = qJD(3) * t591;
t712 = (-qJ(6) * t690 - t540) * t590 + (qJ(6) * t689 - t684 + (-pkin(5) + t669) * qJD(3)) * t588 + t770;
t711 = (-pkin(8) * qJD(3) + qJ(6) * qJD(4)) * t717 + (qJD(6) * t588 + (-pkin(8) * qJD(4) + qJ(6) * qJD(3)) * t591) * t587 + t760;
t709 = (-t588 * t683 - t591 * t689) * pkin(8) + t760;
t707 = t669 * t691 + t765;
t539 = t639 * qJD(2);
t706 = t587 * t539 + t590 * t762;
t625 = -t587 * t750 + t738;
t705 = t625 * t768 + t764;
t429 = qJ(5) * t695 + t706;
t704 = -qJ(6) * t587 * t694 - t689 * t746 - t429 - t684;
t633 = pkin(4) * t587 - t738;
t703 = t633 * t768 - t764;
t550 = t746 * t590;
t473 = t587 * t762;
t650 = -t539 * t590 + t473;
t702 = qJD(4) * t550 - qJD(6) * t587 - (-qJ(6) * t715 - t588 * t750) * qJD(2) - t650;
t699 = t587 * t544 + t572;
t582 = t588 ^ 2;
t698 = -t591 ^ 2 + t582;
t696 = qJD(2) * t586;
t693 = qJD(3) * t537;
t685 = qJD(5) * t590;
t403 = qJ(6) * t537 + t416;
t682 = qJD(5) - t403;
t675 = t588 * t719;
t673 = pkin(4) * t730 + t461 * t774;
t672 = pkin(4) * t729 + t463 * t774;
t671 = -pkin(4) * t726 - t518 * t774;
t670 = -g(1) * t727 - g(2) * t728 + g(3) * t675;
t667 = t589 * t696;
t666 = t592 * t696;
t665 = t768 * t692;
t664 = t768 * t683;
t663 = t402 * t689;
t662 = t402 * t687;
t661 = t768 * t689;
t423 = t462 * t590 + t513 * t587;
t653 = -t422 * pkin(4) + qJ(5) * t423;
t425 = t464 * t590 + t515 * t587;
t652 = -t424 * pkin(4) + qJ(5) * t425;
t470 = t519 * t590 - t587 * t719;
t651 = -t469 * pkin(4) + qJ(5) * t470;
t571 = pkin(8) * t718;
t649 = t544 * t590 - t571;
t644 = -qJD(3) * t563 - t588 * t493 - t541 * t690 + t591 * t648;
t643 = t537 * t668;
t642 = t588 * t668;
t641 = t535 * t668;
t638 = g(1) * t515 + g(2) * t513;
t480 = -qJ(5) * t591 + t699;
t634 = pkin(4) * t590 + t739;
t412 = -qJDD(3) * pkin(3) - t644;
t405 = -pkin(4) * t768 + t681;
t631 = t405 * t590 - t406 * t587;
t629 = pkin(8) + t633;
t594 = qJD(3) ^ 2;
t626 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t594 - t632;
t392 = t646 - t757;
t624 = t533 * t587 + t687 * t768;
t623 = t533 * t590 - t661;
t621 = t590 * t411 + t587 * t433 - t476 * t689 + t488 * t687;
t620 = -pkin(8) + t625;
t439 = -t513 * t718 - t514 * t590;
t441 = -t515 * t718 - t516 * t590;
t490 = t758 * t586;
t619 = g(1) * t441 + g(2) * t439 + g(3) * t490;
t440 = -t513 * t715 + t514 * t587;
t442 = -t515 * t715 + t516 * t587;
t618 = -g(1) * t442 - g(2) * t440 - g(3) * t491;
t617 = g(1) * t463 + g(2) * t461 - g(3) * t518;
t616 = g(1) * t464 + g(2) * t462 + g(3) * t519;
t613 = t586 * pkin(3) * t714 + pkin(2) * t719 + t491 * pkin(4) + pkin(8) * t721 + pkin(9) * t675 + qJ(5) * t490;
t393 = t449 * pkin(4) + t448 * qJ(5) - t537 * qJD(5) + t412;
t391 = -pkin(5) * t449 + qJDD(6) - t393;
t611 = t391 - t617;
t610 = -t640 * t768 - t747;
t390 = t621 + t772;
t608 = t440 * pkin(4) + pkin(8) * t514 - pkin(9) * t728 + qJ(5) * t439 + t513 * t775;
t607 = t442 * pkin(4) + pkin(8) * t516 - pkin(9) * t727 + qJ(5) * t441 + t515 * t775;
t605 = -t745 * t768 - t617;
t603 = -t393 + t605;
t602 = t448 - t724;
t542 = -t668 - t744;
t601 = -pkin(8) * qJDD(3) + (t542 + t668 - t744) * qJD(3);
t599 = g(1) * t425 + g(2) * t423 + g(3) * t470 - t621;
t597 = t418 * t537 - t600;
t596 = t416 * t768 + t599;
t595 = qJD(2) ^ 2;
t581 = t591 * pkin(4);
t549 = t746 * t587;
t543 = -pkin(3) - t634;
t531 = pkin(3) - t753;
t506 = t629 * t588;
t481 = t581 - t649;
t479 = t620 * t588;
t468 = qJD(3) * t519 + t588 * t666;
t467 = -qJD(3) * t518 + t591 * t666;
t460 = pkin(4) * t537 + t740;
t455 = t587 * t736 + t480;
t450 = pkin(5) * t591 + t571 + t581 + (-t544 - t736) * t590;
t434 = -t537 * t750 - t740;
t432 = (qJD(4) * t634 - t685) * t588 + t629 * t690;
t431 = -pkin(4) * t695 + t650;
t409 = (qJD(4) * t753 + t685) * t588 + t620 * t690;
t401 = -qJD(4) * t469 + t467 * t590 + t587 * t667;
t400 = qJD(4) * t470 + t467 * t587 - t590 * t667;
t394 = -t750 * t768 + t682;
t389 = t390 + t766;
t388 = -qJD(6) * t537 + t392 + t737 - t748;
t1 = [t713 * MDP(1) + (-qJD(3) * t468 - qJDD(3) * t518) * MDP(10) + (-qJD(3) * t467 - qJDD(3) * t519) * MDP(11) + (t390 * t470 + t392 * t469 + t393 * t518 + t400 * t405 + t401 * t406 + t418 * t468 - g(3)) * MDP(22) + (t388 * t469 + t389 * t470 - t391 * t518 + t394 * t400 + t399 * t401 - t402 * t468 - g(3)) * MDP(26) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t591 + MDP(11) * t588 - MDP(3)) * t595) * t589 + (MDP(10) * t759 + t622 * MDP(11) + qJDD(2) * MDP(3) - t595 * MDP(4)) * t592) * t586 + t776 * (t400 * t537 - t401 * t535 - t448 * t469 - t449 * t470) + (MDP(18) - t767) * (-t401 * t768 - t448 * t518 + t468 * t537 - t470 * t533) + (MDP(17) + t777) * (-t400 * t768 + t518 * t449 + t468 * t535 - t469 * t533); qJDD(2) * MDP(2) + (t713 * t719 + t638) * MDP(3) + (g(1) * t516 + g(2) * t514 - t713 * t721) * MDP(4) + (qJDD(2) * t582 + 0.2e1 * t588 * t658) * MDP(5) + 0.2e1 * (t579 * t588 - t678 * t698) * MDP(6) + (qJDD(3) * t588 + t591 * t594) * MDP(7) + (qJDD(3) * t591 - t588 * t594) * MDP(8) + (t601 * t588 + ((-g(3) * t592 + t660) * t586 + t626 + t638) * t591) * MDP(10) + (t601 * t591 + (-t626 - t554) * t588 + t670) * MDP(11) + (-t448 * t717 + (-t587 * t688 + t591 * t683) * t537) * MDP(12) + ((-t535 * t590 - t537 * t587) * t690 + (t731 - t449 * t590 + (t535 * t587 - t722) * qJD(4)) * t588) * MDP(13) + ((t448 + t664) * t591 + (t623 + t693) * t588) * MDP(14) + ((t449 - t665) * t591 + (-qJD(3) * t535 - t624) * t588) * MDP(15) + (-t533 * t591 + t691 * t768) * MDP(16) + (t649 * t533 - t765 * t768 + ((pkin(8) * t535 - t587 * t640) * qJD(3) + t646) * t591 + (-t641 - t640 * t687 + t416 * qJD(3) + t412 * t587 + (t449 + t665) * pkin(8)) * t588 + t618) * MDP(17) + (-t699 * t533 - t773 * t768 + (-t640 * t683 + (t661 + t693) * pkin(8) + t621) * t591 + (-t643 + t640 * t689 - qJD(3) * t417 + t412 * t590 + (-t448 + t664) * pkin(8)) * t588 + t619) * MDP(18) + (t432 * t535 + t449 * t506 - t481 * t533 + (t418 * t692 + t392) * t591 - t707 * t768 + (-qJD(3) * t405 + t393 * t587 + t418 * t687 - t641) * t588 + t618) * MDP(19) + (-t448 * t481 - t449 * t480 + t707 * t537 - t709 * t535 + t631 * t690 + (-t390 * t587 + t392 * t590 + (-t405 * t587 - t406 * t590) * qJD(4)) * t588 - t670) * MDP(20) + (-t432 * t537 + t448 * t506 + t480 * t533 + (-t418 * t683 - t390) * t591 + t709 * t768 + (qJD(3) * t406 - t393 * t590 + t418 * t689 + t643) * t588 - t619) * MDP(21) + (t390 * t480 + t393 * t506 + t392 * t481 - g(1) * t607 - g(2) * t608 - g(3) * t613 + (t432 - t642) * t418 + t709 * t406 + t707 * t405) * MDP(22) + (-t409 * t535 - t449 * t479 - t450 * t533 + (-t402 * t692 + t388) * t591 - t712 * t768 + (-qJD(3) * t394 - t391 * t587 - t641 - t662) * t588 + t618) * MDP(23) + (t409 * t537 - t448 * t479 + t455 * t533 + (t402 * t683 - t389) * t591 + t711 * t768 + (qJD(3) * t399 + t391 * t590 + t643 - t663) * t588 - t619) * MDP(24) + (t448 * t450 + t449 * t455 - t712 * t537 + t711 * t535 + (-t394 * t590 + t399 * t587) * t690 + (-t388 * t590 + t389 * t587 + (t394 * t587 + t399 * t590) * qJD(4)) * t588 + t670) * MDP(25) + (t389 * t455 + t388 * t450 + t391 * t479 - g(1) * (pkin(5) * t442 + qJ(6) * t727 + t607) - g(2) * (pkin(5) * t440 + qJ(6) * t728 + t608) - g(3) * (pkin(5) * t491 - qJ(6) * t675 + t613) + (t409 + t642) * t402 + t711 * t399 + t712 * t394) * MDP(26); MDP(7) * t676 + MDP(8) * t579 + qJDD(3) * MDP(9) + (qJD(3) * t486 - t542 * t695 - t617 + t644) * MDP(10) + (-t636 + (-qJD(2) * t542 - t493) * t591 + t616) * MDP(11) + (t722 * t768 - t731) * MDP(12) + ((-t448 - t724) * t590 + (-t449 - t723) * t587) * MDP(13) + ((-t537 * t588 - t715 * t768) * qJD(2) + t624) * MDP(14) + ((t535 * t588 + t718 * t768) * qJD(2) + t623) * MDP(15) - t768 * MDP(16) * t695 + (-t416 * t695 - pkin(3) * t449 + t473 * t768 - t486 * t535 + t610 * t587 + (-t412 - (t539 + t745) * t768 - t617) * t590) * MDP(17) + (pkin(3) * t448 + t706 * t768 + t417 * t695 - t486 * t537 + t610 * t590 + (t412 - t605) * t587) * MDP(18) + (t405 * t695 + t431 * t768 + t449 * t543 + t703 * t535 - t587 * t755 + t603 * t590) * MDP(19) + (t429 * t535 - t431 * t537 + (t390 + t768 * t405 + (qJD(4) * t537 - t449) * pkin(9)) * t590 + (t392 - t733 + (qJD(4) * t535 - t448) * pkin(9)) * t587 - t616) * MDP(20) + (-t406 * t695 - t429 * t768 + t448 * t543 - t703 * t537 + t603 * t587 + t590 * t755) * MDP(21) + (t393 * t543 - t406 * t429 - t405 * t431 - g(1) * t672 - g(2) * t673 - g(3) * t671 + t703 * t418 + (qJD(4) * t631 + t390 * t590 + t392 * t587 - t616) * pkin(9)) * MDP(22) + (-t663 - t449 * t531 - t533 * t549 - t702 * t768 - t705 * t535 + (t394 * t588 + t402 * t718) * qJD(2) + t611 * t590) * MDP(23) + (t662 - t448 * t531 + t533 * t550 + t704 * t768 + t705 * t537 + (-t399 * t588 - t402 * t715) * qJD(2) + t611 * t587) * MDP(24) + (t448 * t549 + t449 * t550 - t702 * t537 + t704 * t535 + (-t394 * t768 - t389) * t590 + (-t388 + t734) * t587 + t616) * MDP(25) + (t389 * t550 + t388 * t549 + t391 * t531 - g(1) * (pkin(5) * t729 + t464 * t746 + t672) - g(2) * (pkin(5) * t730 + t462 * t746 + t673) - g(3) * (-pkin(5) * t726 + t519 * t746 + t671) + t705 * t402 + t704 * t399 + t702 * t394) * MDP(26) + (-MDP(5) * t588 * t591 + MDP(6) * t698) * t595; MDP(12) * t725 + (t532 - t751) * MDP(13) - t602 * MDP(14) + (-t449 + t723) * MDP(15) + t533 * MDP(16) + (t537 * t640 + t604 + t732) * MDP(17) + (-t535 * t640 + t596) * MDP(18) + (-t460 * t535 + t524 - t597 + t732) * MDP(19) + (pkin(4) * t448 - t741 + (t406 - t417) * t537 + (t405 - t681) * t535) * MDP(20) + (-t418 * t535 + t460 * t537 - t596 + t763) * MDP(21) + (-t392 * pkin(4) - g(1) * t652 - g(2) * t653 - g(3) * t651 + t390 * qJ(5) - t405 * t417 + t406 * t681 - t418 * t460) * MDP(22) + (t404 * t768 + t434 * t535 + (pkin(5) + t750) * t533 + t754) * MDP(23) + (t402 * t535 - t403 * t768 - t434 * t537 - t599 + t763 + t766) * MDP(24) + (t741 - t448 * t750 + (-t399 + t404) * t537 + (-t394 + t682) * t535) * MDP(25) + (t389 * qJ(5) - t388 * t750 - t394 * t404 - t402 * t434 - g(1) * (-pkin(5) * t424 + t652) - g(2) * (-pkin(5) * t422 + t653) - g(3) * (-pkin(5) * t469 + t651) + t682 * t399) * MDP(26); (t597 - t733) * MDP(22) + (-t734 - t748 - t754) * MDP(26) + t777 * (-t533 + t725) + t767 * (-t768 ^ 2 - t532) - t776 * t602; (-t627 - t723) * MDP(23) + (t645 - t724) * MDP(24) + (-t532 - t751) * MDP(25) + (t394 * t537 - t399 * t535 + t611) * MDP(26) + (MDP(23) * t771 - MDP(24) * t657) * t587;];
tau  = t1;
