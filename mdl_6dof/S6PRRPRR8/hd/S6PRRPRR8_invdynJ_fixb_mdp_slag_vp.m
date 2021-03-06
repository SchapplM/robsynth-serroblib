% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:59
% EndTime: 2019-03-08 22:41:14
% DurationCPUTime: 11.12s
% Computational Cost: add. (5176->667), mult. (13085->917), div. (0->0), fcn. (11076->14), ass. (0->299)
t608 = sin(pkin(7));
t615 = sin(qJ(3));
t735 = qJD(3) * t615;
t705 = t608 * t735;
t590 = pkin(3) * t705;
t619 = cos(qJ(3));
t788 = qJ(4) * t619;
t672 = pkin(10) * t615 - t788;
t733 = qJD(4) * t615;
t629 = qJD(3) * t672 - t733;
t616 = sin(qJ(2));
t609 = sin(pkin(6));
t742 = qJD(1) * t609;
t710 = t616 * t742;
t686 = t608 * t710;
t810 = -t608 * t629 - t590 + t686;
t611 = cos(pkin(7));
t760 = t616 * t619;
t620 = cos(qJ(2));
t761 = t615 * t620;
t655 = t611 * t760 + t761;
t526 = t655 * t609;
t511 = qJD(1) * t526;
t767 = t611 * t615;
t602 = pkin(2) * t767;
t773 = t608 * t619;
t794 = pkin(4) + pkin(9);
t809 = (t773 * t794 + t602) * qJD(3) - t511;
t724 = qJD(2) * qJD(3);
t702 = t615 * t724;
t681 = t608 * t702;
t720 = qJDD(2) * t619;
t806 = -t608 * t720 + t681;
t763 = t615 * t616;
t715 = t611 * t763;
t687 = t609 * t715;
t734 = qJD(3) * t619;
t703 = t611 * t734;
t709 = t620 * t742;
t808 = -pkin(2) * t703 - qJD(1) * t687 + t619 * t709;
t612 = cos(pkin(6));
t741 = qJD(1) * t612;
t711 = t608 * t741;
t736 = qJD(2) * t620;
t792 = pkin(9) * t608;
t807 = qJDD(2) * t792 + (qJD(1) * t736 + qJDD(1) * t616) * t609 + qJD(3) * t711;
t740 = qJD(2) * t608;
t554 = pkin(9) * t740 + t710;
t567 = qJD(2) * pkin(2) + t709;
t766 = t611 * t619;
t470 = t554 * t615 - t567 * t766 - t619 * t711;
t725 = qJD(4) + t470;
t739 = qJD(2) * t611;
t597 = qJD(3) + t739;
t614 = sin(qJ(5));
t618 = cos(qJ(5));
t737 = qJD(2) * t619;
t707 = t608 * t737;
t528 = t597 * t614 + t618 * t707;
t518 = qJD(6) + t528;
t738 = qJD(2) * t615;
t708 = t608 * t738;
t581 = qJD(5) + t708;
t617 = cos(qJ(6));
t684 = t614 * t707;
t530 = t597 * t618 - t684;
t613 = sin(qJ(6));
t782 = t530 * t613;
t479 = -t581 * t617 + t782;
t805 = t479 * t581;
t775 = t608 * t615;
t599 = pkin(9) * t775;
t712 = -pkin(2) * t619 - pkin(3);
t490 = pkin(4) * t775 + t599 + (-pkin(10) + t712) * t611;
t621 = -pkin(3) - pkin(10);
t789 = qJ(4) * t615;
t657 = t619 * t621 - t789;
t503 = (-pkin(2) + t657) * t608;
t731 = qJD(5) * t618;
t732 = qJD(5) * t614;
t804 = -t490 * t731 + t503 * t732 - t614 * t809 + t618 * t810;
t603 = t611 * qJD(4);
t751 = -t705 * t794 + t603 - t808;
t803 = t490 * t614 + t503 * t618;
t802 = pkin(9) * t705 + t808;
t726 = pkin(4) * t708 + t725;
t801 = MDP(10) - MDP(13);
t800 = MDP(11) - MDP(14);
t594 = t611 * t741;
t673 = -pkin(3) * t619 - t789;
t478 = t594 + (qJD(2) * t673 - t567) * t608;
t799 = t478 * t708 + qJDD(4);
t701 = t619 * t724;
t721 = qJDD(2) * t615;
t646 = t701 + t721;
t632 = t646 * t608;
t535 = qJDD(5) + t632;
t722 = qJDD(2) * t611;
t596 = qJDD(3) + t722;
t770 = t609 * t620;
t595 = qJDD(1) * t770;
t771 = t609 * t616;
t706 = qJD(2) * t771;
t683 = qJD(1) * t706;
t531 = qJDD(2) * pkin(2) + t595 - t683;
t542 = t567 * t767;
t723 = qJDD(1) * t612;
t700 = t608 * t723;
t658 = qJD(3) * t542 - t531 * t766 + t554 * t734 + t615 * t807 - t619 * t700;
t637 = qJDD(4) + t658;
t412 = pkin(4) * t632 + t596 * t621 + t637;
t443 = t597 * t621 + t726;
t467 = t594 + (qJD(2) * t657 - t567) * t608;
t422 = t443 * t614 + t467 * t618;
t593 = t611 * t723;
t746 = pkin(3) * t681 + t593;
t428 = (qJD(2) * t629 + qJDD(2) * t657 - t531) * t608 + t746;
t627 = -qJD(5) * t422 + t412 * t618 - t428 * t614;
t405 = -pkin(5) * t535 - t627;
t791 = sin(pkin(12));
t697 = t791 * t616;
t610 = cos(pkin(12));
t768 = t610 * t620;
t547 = t612 * t768 - t697;
t696 = t791 * t620;
t769 = t610 * t616;
t548 = t612 * t769 + t696;
t772 = t609 * t610;
t719 = t608 * t772;
t461 = -t547 * t766 + t548 * t615 + t619 * t719;
t549 = -t612 * t696 - t769;
t550 = -t612 * t697 + t768;
t698 = t609 * t791;
t680 = t608 * t698;
t463 = -t549 * t766 + t550 * t615 - t619 * t680;
t494 = -t547 * t608 - t611 * t772;
t495 = -t549 * t608 + t611 * t698;
t716 = t609 * t763;
t717 = t611 * t770;
t491 = -t612 * t773 - t619 * t717 + t716;
t546 = -t608 * t770 + t611 * t612;
t664 = t491 * t618 - t546 * t614;
t642 = g(1) * (t463 * t618 - t495 * t614) + g(2) * (t461 * t618 - t494 * t614) + g(3) * t664;
t798 = t518 * (pkin(5) * t530 + pkin(11) * t518) + t405 + t642;
t456 = -qJD(5) * t684 + t596 * t614 + (qJD(5) * t597 - t806) * t618;
t451 = qJDD(6) + t456;
t573 = pkin(5) * t614 - pkin(11) * t618 + qJ(4);
t679 = pkin(5) * t618 + pkin(11) * t614;
t797 = t518 * ((-pkin(4) - t679) * t708 - qJD(5) * t679 - t725) - t573 * t451;
t473 = t547 * t619 - t548 * t767;
t475 = t549 * t619 - t550 * t767;
t527 = (t619 * t620 - t715) * t609;
t745 = pkin(9) * t773 + t602;
t747 = qJD(3) * t745 - t511;
t796 = -g(1) * t475 - g(2) * t473 - g(3) * t527 - t747 * t597;
t795 = -qJD(5) * t803 + t614 * t810 + t618 * t809;
t793 = pkin(3) * t596;
t790 = MDP(8) * t608;
t455 = -qJD(5) * t528 + t618 * t596 + t614 * t806;
t728 = qJD(6) * t617;
t713 = t455 * t617 + t535 * t613 + t581 * t728;
t729 = qJD(6) * t613;
t424 = -t530 * t729 + t713;
t787 = t424 * t613;
t786 = t479 * t518;
t481 = t530 * t617 + t581 * t613;
t785 = t481 * t518;
t784 = t528 * t581;
t783 = t530 * t581;
t781 = t535 * t614;
t779 = t581 * t614;
t778 = t581 * t618;
t777 = t596 * MDP(9);
t582 = t596 * qJ(4);
t605 = t608 ^ 2;
t622 = qJD(2) ^ 2;
t776 = t605 * t622;
t774 = t608 * t618;
t765 = t613 * t451;
t764 = t613 * t621;
t762 = t615 * t617;
t759 = t616 * t622;
t758 = t617 * t451;
t757 = t617 * t621;
t756 = t618 * t424;
t755 = t621 * t535;
t754 = qJDD(1) - g(3);
t704 = t608 * t734;
t753 = -pkin(5) * t704 - t795;
t471 = t554 * t619 + t615 * t711 + t542;
t460 = pkin(4) * t707 + t471;
t592 = pkin(3) * t708;
t508 = t672 * t740 + t592;
t749 = t460 * t614 + t508 * t618;
t748 = -t603 + t802;
t606 = t615 ^ 2;
t744 = -t619 ^ 2 + t606;
t743 = MDP(12) * t608;
t730 = qJD(5) * t621;
t727 = qJD(3) - t597;
t718 = t614 * t773;
t714 = t412 * t614 + t428 * t618 + t443 * t731;
t532 = -qJ(4) * t611 - t745;
t644 = -t467 * t732 + t714;
t404 = pkin(11) * t535 + t644;
t585 = t597 * qJD(4);
t659 = -t531 * t767 + t554 * t735 - t567 * t703 - t615 * t700 - t619 * t807;
t418 = -t582 - t585 + t659;
t413 = -pkin(4) * t806 - t418;
t407 = pkin(5) * t456 - pkin(11) * t455 + t413;
t694 = -t404 * t613 + t407 * t617;
t692 = t455 * t613 - t535 * t617;
t691 = t518 * t617;
t690 = t597 + t739;
t689 = t596 + t722;
t688 = t615 * t619 * t776;
t502 = pkin(4) * t773 - t532;
t685 = t608 * t706;
t678 = g(1) * t550 + g(2) * t548;
t514 = (t613 * t619 + t614 * t762) * t740;
t676 = -t617 * t732 - t514;
t551 = t611 * t614 + t618 * t773;
t552 = t611 * t618 - t718;
t449 = pkin(5) * t551 - pkin(11) * t552 + t502;
t675 = -pkin(11) * t704 - qJD(6) * t449 + t804;
t445 = pkin(11) * t775 + t803;
t498 = -qJD(5) * t551 + t614 * t705;
t499 = -qJD(5) * t718 + t611 * t731 - t618 * t705;
t674 = -pkin(5) * t499 + pkin(11) * t498 + qJD(6) * t445 - t751;
t586 = t597 * qJ(4);
t448 = t586 + t460;
t671 = t404 * t617 + t407 * t613;
t417 = pkin(11) * t581 + t422;
t426 = pkin(5) * t528 - pkin(11) * t530 + t448;
t409 = t417 * t617 + t426 * t613;
t670 = t417 * t613 - t426 * t617;
t421 = t443 * t618 - t467 * t614;
t466 = t491 * t614 + t546 * t618;
t654 = t611 * t761 + t760;
t492 = t609 * t654 + t612 * t775;
t432 = t466 * t617 + t492 * t613;
t431 = -t466 * t613 + t492 * t617;
t666 = t490 * t618 - t503 * t614;
t663 = t615 * t683;
t662 = t619 * t683;
t533 = (-pkin(2) + t673) * t608;
t661 = qJD(3) * (-qJD(2) * t533 - t478);
t656 = -t552 * t613 + t608 * t762;
t501 = t552 * t617 + t613 * t775;
t653 = -t518 * t728 - t765;
t652 = -t518 * t729 + t758;
t649 = t581 * t481;
t438 = (-pkin(3) * t720 - qJ(4) * t646 - qJD(2) * t733 - t531) * t608 + t746;
t510 = t590 + (-qJ(4) * t734 - t733) * t608;
t647 = qJD(2) * t510 + qJDD(2) * t533 + t438;
t645 = t702 - t720;
t641 = g(1) * t463 + g(2) * t461 + g(3) * t491;
t462 = t547 * t767 + t548 * t619 - t615 * t719;
t464 = t550 * t619 + (t549 * t611 + t680) * t615;
t640 = g(1) * t464 + g(2) * t462 + g(3) * t492;
t472 = t547 * t615 + t548 * t766;
t474 = t549 * t615 + t550 * t766;
t639 = g(1) * t474 + g(2) * t472 + g(3) * t526;
t635 = -g(3) * t771 - t678;
t633 = t413 - t640;
t416 = -pkin(5) * t581 - t421;
t631 = -pkin(11) * t451 + (t416 + t421) * t518;
t630 = t608 * (t727 * t737 + t721);
t628 = qJD(6) * t518 * t621 + t640;
t626 = (pkin(11) * t707 - qJD(6) * t573 + t749) * t518 + t641;
t625 = t641 - t658;
t624 = -t640 - t659;
t623 = t471 * t597 + t625;
t540 = -qJ(4) * t707 + t592;
t534 = t611 * t712 + t599;
t525 = t618 * t535;
t516 = -t567 * t608 + t594;
t513 = t613 * t614 * t708 - t617 * t707;
t488 = -t531 * t608 + t593;
t483 = t526 * t614 + t771 * t774;
t458 = -t586 - t471;
t457 = -pkin(3) * t597 + t725;
t453 = -qJD(2) * t687 - qJD(3) * t716 + (t609 * t736 + (t608 * t612 + t717) * qJD(3)) * t619;
t452 = t612 * t705 + (qJD(2) * t655 + qJD(3) * t654) * t609;
t447 = t474 * t614 + t550 * t774;
t446 = t472 * t614 + t548 * t774;
t444 = -pkin(5) * t775 - t666;
t442 = qJD(6) * t501 + t498 * t613 - t617 * t704;
t441 = qJD(6) * t656 + t498 * t617 + t613 * t704;
t437 = t463 * t614 + t495 * t618;
t435 = t461 * t614 + t494 * t618;
t429 = -pkin(5) * t707 - t460 * t618 + t508 * t614;
t425 = qJD(6) * t481 + t692;
t423 = t637 - t793;
t420 = qJD(5) * t664 + t452 * t614 + t618 * t685;
t419 = qJD(5) * t466 - t452 * t618 + t614 * t685;
t403 = -qJD(6) * t409 + t694;
t402 = -t670 * qJD(6) + t671;
t1 = [t754 * MDP(1) + (-t418 * t492 + t423 * t491 + t438 * t546 + t452 * t457 - t453 * t458 - g(3)) * MDP(15) + (-t419 * t581 + t453 * t528 + t456 * t492 + t535 * t664) * MDP(21) + (-t420 * t581 + t453 * t530 + t455 * t492 - t466 * t535) * MDP(22) + ((-qJD(6) * t432 - t420 * t613 + t453 * t617) * t518 + t431 * t451 + t419 * t479 - t664 * t425) * MDP(28) + (-(qJD(6) * t431 + t420 * t617 + t453 * t613) * t518 - t432 * t451 + t419 * t481 - t664 * t424) * MDP(29) + ((qJDD(2) * t620 - t759) * MDP(3) + (-qJDD(2) * t616 - t620 * t622) * MDP(4) + (t615 * t800 - t619 * t801) * t605 * t759) * t609 + ((t491 * t615 + t492 * t619) * MDP(12) * qJDD(2) + ((t452 * t615 + t453 * t619 + t491 * t734 - t492 * t735) * MDP(12) + t478 * MDP(15) * t771) * qJD(2) + (t645 * t801 + t646 * t800) * t546) * t608 - t801 * (t452 * t597 + t491 * t596) - t800 * (t453 * t597 + t492 * t596); (t438 * t533 + t418 * t532 + t423 * t534 - g(1) * (pkin(2) * t549 + pkin(3) * t475 + qJ(4) * t474 + t550 * t792) - g(2) * (pkin(2) * t547 + pkin(3) * t473 + qJ(4) * t472 + t548 * t792) - g(3) * (pkin(3) * t527 + qJ(4) * t526 + (pkin(2) * t620 + t616 * t792) * t609) + (t510 - t686) * t478 + t748 * t458 + t747 * t457) * MDP(15) + (-t754 * t771 + t678) * MDP(4) + (-g(1) * t549 - g(2) * t547 - g(3) * t770 + t595) * MDP(3) + (-(t445 * t617 + t449 * t613) * t451 - t402 * t551 - t409 * t499 + t444 * t424 + t405 * t501 + t416 * t441 - g(1) * (-t447 * t613 + t475 * t617) - g(2) * (-t446 * t613 + t473 * t617) - g(3) * (-t483 * t613 + t527 * t617) + (t613 * t674 + t617 * t675) * t518 + t753 * t481) * MDP(29) + (t424 * t656 - t425 * t501 - t441 * t479 - t442 * t481) * MDP(24) + (-t425 * t551 - t442 * t518 + t451 * t656 - t479 * t499) * MDP(26) + ((-t445 * t613 + t449 * t617) * t451 + t403 * t551 - t670 * t499 + t444 * t425 - t405 * t656 + t416 * t442 - g(1) * (t447 * t617 + t475 * t613) - g(2) * (t446 * t617 + t473 * t613) - g(3) * (t483 * t617 + t527 * t613) + (t613 * t675 - t617 * t674) * t518 + t753 * t479) * MDP(28) + (t413 * t552 + t448 * t498 + t502 * t455 + t751 * t530 - t535 * t803 + t804 * t581 - t639 * t618) * MDP(22) + ((qJD(3) * t457 - qJDD(2) * t532 - t418) * t619 + (qJD(3) * t458 + qJDD(2) * t534 + t423) * t615 + ((qJD(3) * t534 - t748) * t619 + (qJD(3) * t532 + t747) * t615) * qJD(2) + t635) * t743 + ((-pkin(2) * t645 + t662) * MDP(10) + (-pkin(2) * t646 - t663) * MDP(11) - t662 * MDP(13) + t663 * MDP(14) + 0.2e1 * (t615 * t720 - t724 * t744) * MDP(6) + (qJDD(2) * t606 + 0.2e1 * t615 * t701) * MDP(5)) * t605 + ((-t488 * t619 + t516 * t735) * MDP(10) + (t488 * t615 + t516 * t734) * MDP(11) + (t615 * t661 + t619 * t647) * MDP(13) + (-t615 * t647 + t619 * t661) * MDP(14) + (-t714 * t615 - t422 * t734 + (qJD(5) * t467 * t615 - t635) * t614) * MDP(22) + (t421 * t734 + t615 * t627) * MDP(21) + (-t456 * t615 - t528 * t734) * MDP(19) + (t455 * t615 + t530 * t734) * MDP(18) + (t615 * t689 + t690 * t734) * MDP(7) + (t535 * t615 + t581 * t734) * MDP(20)) * t608 + (-t499 * t581 - t535 * t551) * MDP(19) + (t498 * t581 + t535 * t552) * MDP(18) + qJDD(2) * MDP(2) + t611 * t777 + (t619 * t689 - t690 * t735) * t790 + (-g(1) * t447 - g(2) * t446 - g(3) * t483 + t413 * t551 + t448 * t499 + t502 * t456 + t751 * t528 + t666 * t535 + t795 * t581) * MDP(21) + ((pkin(2) * t766 - t599) * t596 - t658 * t611 + t796) * MDP(10) + (t423 * t611 + t534 * t596 - t796) * MDP(13) + (-t745 * t596 + t802 * t597 + t659 * t611 + t639) * MDP(11) + (t424 * t501 + t441 * t481) * MDP(23) + (t451 * t551 + t499 * t518) * MDP(27) + (t424 * t551 + t441 * t518 + t451 * t501 + t481 * t499) * MDP(25) + (t455 * t552 + t498 * t530) * MDP(16) + (-t455 * t551 - t456 * t552 - t498 * t528 - t499 * t530) * MDP(17) + (-t418 * t611 - t532 * t596 - t748 * t597 - t639) * MDP(14); (-t421 * t707 + qJ(4) * t456 + t726 * t528 + (t755 + (t448 - t460) * t581) * t618 + ((t508 - t730) * t581 + t633) * t614) * MDP(21) + ((-pkin(3) * t615 + t788) * qJDD(2) + ((-qJ(4) * qJD(3) - t458 - t471) * t615 + (-pkin(3) * qJD(3) - t457 + t725) * t619) * qJD(2)) * t743 + (t479 * t514 + t481 * t513 + (t479 * t617 + t481 * t613) * t732 + (-t787 - t425 * t617 + (t479 * t613 - t481 * t617) * qJD(6)) * t618) * MDP(24) + (-t581 * t731 - t781 + (t528 * t619 - t615 * t778) * t740) * MDP(19) + ((-t456 - t783) * t618 + (-t455 + t784) * t614) * MDP(17) + (0.2e1 * t582 + t585 + t725 * t597 + (t478 * t619 + t540 * t615) * t740 + t624) * MDP(14) + (t451 * t614 + t518 * t778) * MDP(27) + (-t581 * t732 + t525 + (-t530 * t619 - t615 * t779) * t740) * MDP(18) + (qJ(4) * t455 + t749 * t581 + t422 * t707 + t726 * t530 + (-t448 * t581 - t755) * t614 + (-t581 * t730 + t633) * t618) * MDP(22) + (t617 * t756 + (-t618 * t729 + t676) * t481) * MDP(23) + t777 + (-t416 * t513 - t429 * t479 - t797 * t617 + t626 * t613 + (-t451 * t764 + t403 + (-t416 * t613 + t479 * t621) * qJD(5) - t628 * t617) * t614 + (-t670 * t708 + t416 * t728 + t405 * t613 - t621 * t425 + (-t518 * t764 - t670) * qJD(5)) * t618) * MDP(28) + (t455 * t618 - t530 * t779) * MDP(16) + MDP(7) * t630 + (-t418 * qJ(4) - t423 * pkin(3) - t478 * t540 - t457 * t471 - g(1) * (-pkin(3) * t463 + qJ(4) * t464) - g(2) * (-pkin(3) * t461 + qJ(4) * t462) - g(3) * (-pkin(3) * t491 + qJ(4) * t492) - t725 * t458) * MDP(15) + (-t516 * t708 + t623) * MDP(10) + (-t470 * t597 - t516 * t707 - t624) * MDP(11) - MDP(5) * t688 - t581 * MDP(20) * t707 + (t424 * t614 + t676 * t518 + (t649 + t652) * t618) * MDP(25) + (-t727 * t738 + t720) * t790 + (-t416 * t514 - t429 * t481 + t797 * t613 + t626 * t617 + (-t451 * t757 - t402 + (-t416 * t617 + t481 * t621) * qJD(5) + t628 * t613) * t614 + (-t409 * t708 - t416 * t729 + t405 * t617 - t621 * t424 + (-t518 * t757 - t409) * qJD(5)) * t618) * MDP(29) + (-t540 * t707 - t623 - 0.2e1 * t793 + t799) * MDP(13) + t744 * MDP(6) * t776 + (-t425 * t614 + (t613 * t732 + t513) * t518 + (t653 - t805) * t618) * MDP(26); MDP(12) * t630 + (t596 + t688) * MDP(13) + (-t597 ^ 2 - t606 * t776) * MDP(14) + (t458 * t597 - t625 - t793 + t799) * MDP(15) + (-t528 * t597 - t581 * t779 + t525) * MDP(21) + (-t530 * t597 - t581 * t778 - t781) * MDP(22) + (-t618 * t425 + (-t597 * t617 - t613 * t778) * t518 + (t653 + t805) * t614) * MDP(28) + (-t756 + (t597 * t613 - t617 * t778) * t518 + (t649 - t652) * t614) * MDP(29); -t528 ^ 2 * MDP(17) + (t455 + t784) * MDP(18) + (-t456 + t783) * MDP(19) + t535 * MDP(20) + (t422 * t581 + t627 - t642) * MDP(21) + (g(1) * t437 + g(2) * t435 + g(3) * t466 + t421 * t581 + t448 * t528 - t644) * MDP(22) + (t481 * t691 + t787) * MDP(23) + ((t424 - t786) * t617 + (-t425 - t785) * t613) * MDP(24) + (t518 * t691 + t765) * MDP(25) + (-t518 ^ 2 * t613 + t758) * MDP(26) + (-pkin(5) * t425 - t422 * t479 + t631 * t613 - t617 * t798) * MDP(28) + (-pkin(5) * t424 - t422 * t481 + t613 * t798 + t631 * t617) * MDP(29) + (MDP(16) * t528 + t530 * MDP(17) - t448 * MDP(21) - t481 * MDP(25) + t479 * MDP(26) - t518 * MDP(27) + MDP(28) * t670 + MDP(29) * t409) * t530; t481 * t479 * MDP(23) + (-t479 ^ 2 + t481 ^ 2) * MDP(24) + (t713 + t786) * MDP(25) + (-t692 + t785) * MDP(26) + t451 * MDP(27) + (t409 * t518 - t416 * t481 - g(1) * (-t437 * t613 + t464 * t617) - g(2) * (-t435 * t613 + t462 * t617) - g(3) * t431 + t694) * MDP(28) + (-t670 * t518 + t416 * t479 - g(1) * (-t437 * t617 - t464 * t613) - g(2) * (-t435 * t617 - t462 * t613) + g(3) * t432 - t671) * MDP(29) + (-MDP(25) * t782 - MDP(26) * t481 - MDP(28) * t409 + MDP(29) * t670) * qJD(6);];
tau  = t1;
