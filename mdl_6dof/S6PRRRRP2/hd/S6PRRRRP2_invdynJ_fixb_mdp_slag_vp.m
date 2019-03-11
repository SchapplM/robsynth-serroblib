% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:39
% EndTime: 2019-03-09 00:04:52
% DurationCPUTime: 10.73s
% Computational Cost: add. (7715->596), mult. (17191->762), div. (0->0), fcn. (13446->14), ass. (0->273)
t605 = sin(qJ(3));
t604 = sin(qJ(4));
t608 = cos(qJ(3));
t731 = t604 * t608;
t777 = cos(qJ(4));
t564 = t605 * t777 + t731;
t703 = qJD(3) + qJD(4);
t519 = t703 * t564;
t691 = t777 * t608;
t706 = qJDD(2) * t605;
t479 = qJD(2) * t519 - qJDD(2) * t691 + t604 * t706;
t475 = qJDD(5) + t479;
t791 = pkin(9) + pkin(8);
t692 = qJD(3) * t791;
t565 = t605 * t692;
t566 = t608 * t692;
t573 = t791 * t605;
t574 = t791 * t608;
t650 = -t573 * t777 - t604 * t574;
t465 = qJD(4) * t650 - t565 * t777 - t604 * t566;
t649 = -t604 * t605 + t691;
t609 = cos(qJ(2));
t601 = sin(pkin(6));
t718 = qJD(1) * t601;
t689 = t609 * t718;
t526 = t649 * t689;
t792 = -t465 + t526;
t603 = sin(qJ(5));
t607 = cos(qJ(5));
t663 = t607 * pkin(5) + t603 * qJ(6);
t570 = -pkin(4) - t663;
t606 = sin(qJ(2));
t690 = t606 * t718;
t676 = qJD(2) * t791 + t690;
t602 = cos(pkin(6));
t717 = qJD(1) * t602;
t524 = t605 * t717 + t608 * t676;
t516 = t777 * t524;
t523 = -t676 * t605 + t608 * t717;
t461 = t604 * t523 + t516;
t713 = qJD(4) * t604;
t667 = pkin(3) * t713 - t461;
t710 = qJD(5) * t607;
t711 = qJD(5) * t603;
t684 = qJD(2) * t777;
t715 = qJD(2) * t605;
t559 = t604 * t715 - t608 * t684;
t743 = t559 * t607;
t744 = t559 * t603;
t790 = -qJD(6) * t603 + (-t710 - t743) * qJ(6) + (t711 + t744) * pkin(5);
t764 = cos(pkin(11));
t678 = t764 * t606;
t600 = sin(pkin(11));
t737 = t600 * t609;
t554 = t602 * t678 + t737;
t599 = qJ(3) + qJ(4);
t595 = sin(t599);
t596 = cos(t599);
t679 = t601 * t764;
t507 = -t554 * t595 - t596 * t679;
t677 = t764 * t609;
t738 = t600 * t606;
t556 = -t602 * t738 + t677;
t739 = t600 * t601;
t509 = -t556 * t595 + t596 * t739;
t736 = t601 * t606;
t542 = -t595 * t736 + t596 * t602;
t788 = g(1) * t509 + g(2) * t507 + g(3) * t542;
t683 = qJD(4) * t777;
t670 = pkin(3) * t683;
t552 = qJD(5) + t559;
t753 = t475 * t607;
t787 = pkin(10) * (t552 * t711 - t753);
t553 = -t602 * t677 + t738;
t555 = t602 * t737 + t678;
t666 = g(1) * t555 + g(2) * t553;
t734 = t601 * t609;
t636 = -g(3) * t734 + t666;
t786 = t636 * t595;
t704 = t602 * qJDD(1);
t579 = t608 * t704;
t708 = qJD(1) * qJD(2);
t541 = qJDD(2) * pkin(8) + (qJDD(1) * t606 + t609 * t708) * t601;
t675 = pkin(9) * qJDD(2) + t541;
t449 = qJDD(3) * pkin(3) - qJD(3) * t524 - t675 * t605 + t579;
t456 = qJD(3) * t523 + t605 * t704 + t675 * t608;
t765 = qJD(3) * pkin(3);
t517 = t523 + t765;
t624 = t604 * t449 + t456 * t777 + t517 * t683 - t524 * t713;
t702 = qJDD(3) + qJDD(4);
t409 = pkin(10) * t702 + t624;
t682 = t606 * t708;
t661 = -qJDD(1) * t734 + t601 * t682;
t762 = qJDD(2) * pkin(2);
t540 = t661 - t762;
t707 = qJD(2) * qJD(3);
t681 = t605 * t707;
t587 = pkin(3) * t681;
t518 = t703 * t649;
t614 = t518 * qJD(2) + t564 * qJDD(2);
t705 = qJDD(2) * t608;
t424 = -pkin(3) * t705 + t479 * pkin(4) - pkin(10) * t614 + t540 + t587;
t460 = t604 * t517 + t516;
t452 = pkin(10) * t703 + t460;
t593 = t608 * pkin(3) + pkin(2);
t550 = -qJD(2) * t593 - t689;
t561 = -qJD(2) * t731 - t605 * t684;
t478 = pkin(4) * t559 + pkin(10) * t561 + t550;
t645 = t607 * t409 + t603 * t424 - t452 * t711 + t478 * t710;
t763 = qJ(6) * t475;
t399 = qJD(6) * t552 + t645 + t763;
t669 = t603 * t409 - t607 * t424 + t452 * t710 + t478 * t711;
t776 = pkin(5) * t475;
t401 = qJDD(6) + t669 - t776;
t785 = t399 * t607 + t401 * t603;
t725 = -t790 - t667;
t532 = -t604 * t573 + t574 * t777;
t631 = t564 * t734;
t722 = -qJD(1) * t631 + qJD(4) * t532 - t604 * t565 + t566 * t777;
t699 = t605 * t765;
t457 = pkin(4) * t519 - pkin(10) * t518 + t699;
t505 = -pkin(4) * t649 - pkin(10) * t564 - t593;
t784 = -t505 * t710 + t532 * t711 + t792 * t607 + (-t457 + t690) * t603;
t721 = t603 * t505 + t607 * t532;
t557 = t602 * t608 - t605 * t736;
t735 = t601 * t608;
t558 = t602 * t605 + t606 * t735;
t490 = t604 * t557 + t558 * t777;
t577 = t603 * t734;
t477 = t490 * t607 - t577;
t782 = t690 - t699;
t611 = qJD(3) ^ 2;
t781 = -pkin(8) * t611 + t601 * (-g(3) * t609 + t682) - t540 + t666 + t762;
t639 = t607 * t561 - t603 * t703;
t441 = -qJD(5) * t639 + t603 * t614 - t607 * t702;
t428 = -t452 * t603 + t478 * t607;
t709 = qJD(6) - t428;
t419 = -pkin(5) * t552 + t709;
t429 = t452 * t607 + t478 * t603;
t420 = qJ(6) * t552 + t429;
t780 = t419 * t710 - t420 * t711 + t785;
t779 = t639 ^ 2;
t778 = t552 ^ 2;
t775 = pkin(5) * t561;
t774 = pkin(10) * t475;
t767 = pkin(10) * qJD(5);
t766 = qJD(2) * pkin(2);
t668 = -t777 * t449 + t604 * t456 + t517 * t713 + t524 * t683;
t410 = -pkin(4) * t702 + t668;
t671 = t607 * t703;
t440 = -qJD(5) * t671 - t561 * t711 - t603 * t702 - t607 * t614;
t402 = t441 * pkin(5) + t440 * qJ(6) + qJD(6) * t639 + t410;
t761 = t402 * t603;
t760 = t419 * t603;
t759 = t429 * t552;
t515 = t604 * t524;
t459 = t517 * t777 - t515;
t451 = -pkin(4) * t703 - t459;
t527 = -t561 * t603 - t671;
t433 = t527 * pkin(5) + qJ(6) * t639 + t451;
t758 = t433 * t559;
t757 = t440 * t603;
t756 = t441 * t607;
t755 = t451 * t559;
t591 = pkin(3) * t604 + pkin(10);
t754 = t475 * t591;
t751 = t527 * t552;
t750 = t527 * t603;
t749 = t527 * t607;
t748 = t639 * t527;
t747 = t639 * t552;
t746 = t639 * t603;
t745 = t639 * t607;
t742 = t564 * t607;
t741 = t596 * t603;
t740 = t596 * t607;
t732 = t603 * t475;
t730 = t607 * t609;
t729 = qJDD(1) - g(3);
t728 = qJ(6) * t519 - qJD(6) * t649 - t784;
t486 = t526 * t603 - t607 * t690;
t727 = -pkin(5) * t519 + qJD(5) * t721 - t457 * t607 + t465 * t603 - t486;
t662 = pkin(5) * t603 - qJ(6) * t607;
t726 = t662 * t518 + (qJD(5) * t663 - qJD(6) * t607) * t564 + t722;
t504 = -pkin(4) * t561 + pkin(10) * t559;
t724 = t607 * t459 + t603 * t504;
t462 = t523 * t777 - t515;
t488 = pkin(3) * t715 + t504;
t723 = t607 * t462 + t603 * t488;
t720 = -t460 + t790;
t597 = t605 ^ 2;
t719 = -t608 ^ 2 + t597;
t716 = qJD(2) * t601;
t714 = qJD(2) * t606;
t712 = qJD(5) * t591;
t700 = t777 * pkin(3);
t697 = t601 * t730;
t696 = t570 * t507;
t695 = t570 * t509;
t694 = t788 * t603;
t693 = t570 * t542;
t687 = t601 * t714;
t686 = t609 * t716;
t680 = t609 * t707;
t674 = -t419 * t561 + t433 * t711;
t673 = t428 * t561 + t451 * t711;
t672 = t552 * t607;
t665 = g(1) * t556 + g(2) * t554;
t658 = t607 * t670;
t657 = t419 * t607 - t420 * t603;
t656 = -t754 + t755;
t612 = qJD(2) ^ 2;
t655 = qJDD(2) * t609 - t606 * t612;
t654 = pkin(4) * t596 + pkin(10) * t595 + t593;
t653 = -g(1) * t600 + g(2) * t764;
t476 = t490 * t603 + t697;
t651 = t557 * t777 - t604 * t558;
t647 = t518 * t603 + t564 * t710;
t646 = -t518 * t607 + t564 * t711;
t481 = -t553 * t741 - t554 * t607;
t483 = -t555 * t741 - t556 * t607;
t536 = t577 * t596 - t607 * t736;
t643 = g(1) * t483 + g(2) * t481 + g(3) * t536;
t482 = -t553 * t740 + t554 * t603;
t484 = -t555 * t740 + t556 * t603;
t537 = (t596 * t730 + t603 * t606) * t601;
t642 = -g(1) * t484 - g(2) * t482 - g(3) * t537;
t508 = t554 * t596 - t595 * t679;
t510 = t556 * t596 + t595 * t739;
t543 = t595 * t602 + t596 * t736;
t640 = g(1) * t510 + g(2) * t508 + g(3) * t543;
t637 = t410 * t603 - t429 * t561 + t451 * t710 + t694;
t635 = -t402 - t788;
t634 = -t410 - t788;
t568 = -t689 - t766;
t632 = -qJD(2) * t568 - t541 + t665;
t630 = t420 * t561 - t433 * t743 - t694 - t761;
t629 = -t591 * t711 + t658;
t467 = t508 * t603 - t553 * t607;
t469 = t510 * t603 - t555 * t607;
t513 = t543 * t603 + t697;
t626 = g(1) * t469 + g(2) * t467 + g(3) * t513 - t669;
t625 = -pkin(8) * qJDD(3) + (t568 + t689 - t766) * qJD(3);
t623 = -t757 - t756 + (-t745 + t750) * qJD(5);
t622 = t550 * t561 - t668 - t788;
t468 = t508 * t607 + t553 * t603;
t470 = t510 * t607 + t555 * t603;
t514 = t543 * t607 - t577;
t621 = -g(1) * t470 - g(2) * t468 - g(3) * t514 + t645;
t620 = -t433 * t639 + qJDD(6) - t626;
t619 = ((-t440 - t751) * t607 + (-t441 + t747) * t603) * MDP(20) + (-t639 * t672 - t757) * MDP(19) + (-t527 * t561 - t603 * t778 + t753) * MDP(22) + (t552 * t672 - t561 * t639 + t732) * MDP(21) + (t559 * t703 + t614) * MDP(14) + (-t561 * t703 - t479) * MDP(15) + (-t559 ^ 2 + t561 ^ 2) * MDP(13) + t702 * MDP(16) + (-MDP(12) * t559 + MDP(23) * t552) * t561;
t618 = t419 * t743 - t420 * t744 - t640 + t780;
t615 = t550 * t559 - t624 + t640;
t592 = -t700 - pkin(4);
t562 = -t700 + t570;
t549 = t561 * qJ(6);
t521 = -qJD(3) * t558 - t605 * t686;
t520 = qJD(3) * t557 + t608 * t686;
t506 = -qJDD(2) * t593 + t587 + t661;
t471 = -pkin(5) * t639 + qJ(6) * t527;
t464 = t564 * t662 - t650;
t443 = pkin(5) * t649 - t505 * t607 + t532 * t603;
t442 = -qJ(6) * t649 + t721;
t438 = qJD(4) * t490 + t604 * t520 - t521 * t777;
t437 = qJD(4) * t651 + t520 * t777 + t604 * t521;
t435 = t459 * t603 - t504 * t607 + t775;
t434 = -t549 + t724;
t432 = t462 * t603 - t488 * t607 + t775;
t431 = -t549 + t723;
t425 = -t440 + t751;
t413 = qJD(5) * t477 + t437 * t603 - t607 * t687;
t412 = -qJD(5) * t476 + t437 * t607 + t603 * t687;
t1 = [t729 * MDP(1) + (qJD(3) * t521 + qJDD(3) * t557) * MDP(10) + (-qJD(3) * t520 - qJDD(3) * t558) * MDP(11) + (-t438 * t703 + t651 * t702) * MDP(17) + (-t437 * t703 - t490 * t702 - qJDD(2) * t631 + (-t518 * t609 - t606 * t561) * t716) * MDP(18) + (-t412 * t527 - t413 * t639 - t440 * t476 - t441 * t477) * MDP(27) + (t399 * t477 + t401 * t476 - t402 * t651 + t412 * t420 + t413 * t419 + t433 * t438 - g(3)) * MDP(29) + (MDP(24) + MDP(26)) * (-t413 * t552 + t438 * t527 - t441 * t651 - t475 * t476) + (-MDP(25) + MDP(28)) * (t412 * t552 + t438 * t639 - t440 * t651 + t475 * t477) + (t655 * MDP(3) + (-qJDD(2) * t606 - t609 * t612) * MDP(4) + (-t605 * t680 + t608 * t655) * MDP(10) + (-t605 * t655 - t608 * t680) * MDP(11) + (-t479 * t609 + t559 * t714) * MDP(17)) * t601; (-t440 * t742 + t639 * t646) * MDP(19) + (-t440 * t443 - t441 * t442 - t727 * t639 - t728 * t527 + t657 * t518 + t786 + (-t399 * t603 + t401 * t607 + (-t420 * t607 - t760) * qJD(5)) * t564) * MDP(27) + (-t399 * t649 - t402 * t742 + t420 * t519 + t433 * t646 + t440 * t464 + t442 * t475 + t552 * t728 + t639 * t726 - t643) * MDP(28) + (t440 * t649 + t475 * t742 - t519 * t639 - t552 * t646) * MDP(21) + (t410 * t742 - t429 * t519 + t440 * t650 - t646 * t451 - t721 * t475 + t552 * t784 - t639 * t722 + t645 * t649 + t643) * MDP(25) + (t669 * t649 + t428 * t519 - t650 * t441 + t486 * t552 + t722 * t527 + ((-qJD(5) * t532 + t457) * t552 + t505 * t475 + t451 * qJD(5) * t564) * t607 + ((-qJD(5) * t505 - t465) * t552 - t532 * t475 + t410 * t564 + t451 * t518) * t603 + t642) * MDP(24) + (-t593 * t479 - t506 * t649 + t550 * t519 - t559 * t782 + t636 * t596 + t650 * t702 - t703 * t722) * MDP(17) + (t441 * t649 - t519 * t527 - t552 * t647 - t564 * t732) * MDP(22) + (t401 * t649 - t419 * t519 + t433 * t647 + t441 * t464 - t443 * t475 + t527 * t726 - t552 * t727 + t564 * t761 + t642) * MDP(26) + (-t519 * t703 + t649 * t702) * MDP(15) + (-t475 * t649 + t519 * t552) * MDP(23) + (-t564 * t479 - t518 * t559 + t561 * t519 + t614 * t649) * MDP(13) + (t399 * t442 + t402 * t464 + t401 * t443 - g(1) * (pkin(5) * t484 + qJ(6) * t483 + t556 * t791) - g(2) * (pkin(5) * t482 + qJ(6) * t481 + t554 * t791) + t726 * t433 + t728 * t420 + t727 * t419 + t666 * t654 + (-pkin(5) * t537 - qJ(6) * t536 - (t606 * t791 + t609 * t654) * t601) * g(3)) * MDP(29) + (qJDD(2) * t597 + 0.2e1 * t608 * t681) * MDP(5) + (-t729 * t736 + t665) * MDP(4) + (t729 * t734 + t666) * MDP(3) + qJDD(2) * MDP(2) + 0.2e1 * (t605 * t705 - t707 * t719) * MDP(6) + ((t746 - t749) * t518 + (t757 - t756 + (t745 + t750) * qJD(5)) * t564) * MDP(20) + (qJDD(3) * t605 + t608 * t611) * MDP(7) + (qJDD(3) * t608 - t605 * t611) * MDP(8) + (t506 * t564 + t550 * t518 - t532 * t702 + t782 * t561 - t593 * t614 + t792 * t703 - t786) * MDP(18) + (t518 * t703 + t564 * t702) * MDP(14) + (t625 * t605 + t608 * t781) * MDP(10) + (-t605 * t781 + t625 * t608) * MDP(11) + (-t561 * t518 + t564 * t614) * MDP(12); (t461 * t703 + (-t559 * t715 + t702 * t777 - t703 * t713) * pkin(3) + t622) * MDP(17) + (t402 * t562 + t670 * t760 - t419 * t432 - g(1) * (pkin(10) * t510 + (-t556 * t605 + t600 * t735) * pkin(3) - t695) - g(2) * (t508 * pkin(10) + (-t554 * t605 - t608 * t679) * pkin(3) - t696) - g(3) * (pkin(3) * t557 + pkin(10) * t543 - t693) - t725 * t433 + (t658 - t431) * t420 + t780 * t591) * MDP(29) + qJDD(3) * MDP(9) + (-t592 * t440 + t656 * t607 - t667 * t639 + (-t629 + t723) * t552 + t637) * MDP(25) + (g(3) * t558 + (-t601 * t653 - t704) * t605 + t632 * t608) * MDP(11) + (t431 * t527 + t432 * t639 + t623 * t591 + t618 + (-t746 - t749) * t670) * MDP(27) + (t462 * t703 + (t561 * t715 - t604 * t702 - t683 * t703) * pkin(3) + t615) * MDP(18) + (-g(3) * t557 + t605 * t632 + t653 * t735 + t579) * MDP(10) + (t592 * t441 + t667 * t527 + ((-t670 + t462) * t552 + t656) * t603 + ((-t488 - t712) * t552 + t634) * t607 + t673) * MDP(24) + MDP(7) * t706 + MDP(8) * t705 + t619 + (t562 * t440 + (-qJD(5) * t433 + t754) * t607 - t725 * t639 + (-t431 + t629) * t552 + t630) * MDP(28) + (t432 * t552 + t562 * t441 - t725 * t527 + (-t552 * t670 - t754 + t758) * t603 + (-t552 * t712 + t635) * t607 + t674) * MDP(26) + (-MDP(5) * t605 * t608 + MDP(6) * t719) * t612; (-pkin(4) * t441 - t460 * t527 + (t459 * t552 + t755 - t774) * t603 + ((-t504 - t767) * t552 + t634) * t607 + t673) * MDP(24) + (pkin(4) * t440 + t451 * t743 + t460 * t639 + t552 * t724 + t637 + t787) * MDP(25) + (t402 * t570 - t420 * t434 - t419 * t435 + g(1) * t695 + g(2) * t696 + g(3) * t693 + t720 * t433 + (qJD(5) * t657 - t640 + t785) * pkin(10)) * MDP(29) + (t435 * t552 + t441 * t570 + (t758 - t774) * t603 + t720 * t527 + (-t552 * t767 + t635) * t607 + t674) * MDP(26) + (pkin(10) * t623 + t434 * t527 + t435 * t639 + t618) * MDP(27) + (t459 * t703 + t615) * MDP(18) + t619 + (-t433 * t710 - t434 * t552 + t440 * t570 + t639 * t720 + t630 - t787) * MDP(28) + (t460 * t703 + t622) * MDP(17); -MDP(19) * t748 + (-t527 ^ 2 + t779) * MDP(20) + t425 * MDP(21) + (-t441 - t747) * MDP(22) + t475 * MDP(23) + (t451 * t639 + t626 + t759) * MDP(24) + (t428 * t552 + t451 * t527 - t621) * MDP(25) + (-t471 * t527 - t620 + t759 + 0.2e1 * t776) * MDP(26) + (pkin(5) * t440 - qJ(6) * t441 - (t420 - t429) * t639 + (t419 - t709) * t527) * MDP(27) + (0.2e1 * t763 - t433 * t527 - t471 * t639 + (0.2e1 * qJD(6) - t428) * t552 + t621) * MDP(28) + (t399 * qJ(6) - t401 * pkin(5) - t433 * t471 - t419 * t429 - g(1) * (-pkin(5) * t469 + qJ(6) * t470) - g(2) * (-pkin(5) * t467 + qJ(6) * t468) - g(3) * (-pkin(5) * t513 + qJ(6) * t514) + t709 * t420) * MDP(29); t425 * MDP(27) + (-t778 - t779) * MDP(28) + (-t420 * t552 + t620 - t776) * MDP(29) + (-t748 - t475) * MDP(26);];
tau  = t1;
