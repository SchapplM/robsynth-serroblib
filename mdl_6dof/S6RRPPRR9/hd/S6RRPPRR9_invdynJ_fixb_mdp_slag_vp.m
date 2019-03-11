% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:31
% EndTime: 2019-03-09 09:32:45
% DurationCPUTime: 13.46s
% Computational Cost: add. (4605->665), mult. (11189->821), div. (0->0), fcn. (8317->10), ass. (0->291)
t596 = sin(pkin(6));
t606 = cos(qJ(2));
t705 = qJD(1) * qJD(2);
t688 = t606 * t705;
t602 = sin(qJ(2));
t703 = qJDD(1) * t602;
t789 = t688 + t703;
t787 = t789 * t596;
t597 = cos(pkin(6));
t719 = qJD(1) * t597;
t576 = qJD(2) + t719;
t601 = sin(qJ(5));
t605 = cos(qJ(5));
t720 = qJD(1) * t596;
t694 = t602 * t720;
t499 = t576 * t605 + t601 * t694;
t704 = qJDD(1) * t597;
t574 = qJDD(2) + t704;
t433 = qJD(5) * t499 + t574 * t601 - t605 * t787;
t693 = t606 * t720;
t539 = qJD(5) + t693;
t758 = t499 * t539;
t788 = -t433 - t758;
t702 = qJDD(1) * t606;
t572 = t596 * t702;
t689 = t602 * t705;
t668 = t596 * t689;
t786 = -t572 + t668;
t700 = pkin(1) * t719;
t727 = -pkin(8) * t694 + t606 * t700;
t708 = qJD(3) - t727;
t532 = t605 * t694;
t497 = t576 * t601 - t532;
t493 = qJD(6) + t497;
t604 = cos(qJ(6));
t600 = sin(qJ(6));
t757 = t499 * t600;
t456 = -t604 * t539 + t757;
t785 = t456 * t539;
t714 = qJD(5) * t601;
t432 = qJD(5) * t532 + t605 * t574 - t576 * t714 + t601 * t787;
t759 = t497 * t539;
t784 = -t432 + t759;
t768 = pkin(1) * t602;
t584 = t597 * t768;
t746 = t596 * t606;
t723 = pkin(8) * t746 + t584;
t500 = -t597 * qJ(3) - t723;
t483 = pkin(3) * t746 - t500;
t455 = pkin(4) * t746 - pkin(9) * t597 + t483;
t687 = -qJ(4) * t606 - pkin(1);
t749 = t596 * t602;
t724 = pkin(2) * t746 + qJ(3) * t749;
t465 = (pkin(9) * t602 + t687) * t596 - t724;
t783 = t601 * t455 + t605 * t465;
t541 = t574 * qJ(3);
t543 = t576 * qJD(3);
t782 = -t541 - t543;
t562 = t576 * pkin(2);
t781 = -t576 * qJ(4) - t562;
t607 = cos(qJ(1));
t734 = t606 * t607;
t603 = sin(qJ(1));
t739 = t602 * t603;
t516 = -t597 * t734 + t739;
t737 = t603 * t606;
t738 = t602 * t607;
t518 = t597 * t737 + t738;
t627 = -g(1) * t518 - g(2) * t516 + g(3) * t746;
t508 = pkin(8) * t693 + t602 * t700;
t487 = pkin(3) * t693 + t508;
t706 = qJD(4) + t487;
t779 = pkin(3) * t572 + qJDD(4);
t547 = t574 * pkin(2);
t778 = t547 - qJDD(3);
t503 = -qJDD(5) + t786;
t769 = pkin(1) * t597;
t699 = qJD(2) * t769;
t673 = qJD(1) * t699;
t698 = pkin(1) * t704;
t671 = pkin(8) * t786 - t602 * t698 - t606 * t673;
t434 = t671 + t782;
t624 = -t434 + t779;
t635 = -t689 + t702;
t410 = -pkin(9) * t574 + (-pkin(3) * t689 + pkin(4) * t635) * t596 + t624;
t549 = t576 * qJ(3);
t707 = pkin(4) * t693 + t706;
t437 = -pkin(9) * t576 + t549 + t707;
t598 = qJ(3) - pkin(9);
t599 = pkin(2) + qJ(4);
t667 = -t599 * t606 - pkin(1);
t623 = -t598 * t602 + t667;
t449 = t623 * t720;
t414 = t437 * t601 + t449 * t605;
t716 = qJD(3) * t602;
t617 = -t716 + (-qJD(2) * t598 - qJD(4)) * t606;
t533 = pkin(2) * t668;
t730 = qJ(4) * t668 + t533;
t415 = (qJD(1) * t617 + qJDD(1) * t623) * t596 + t730;
t616 = -qJD(5) * t414 + t605 * t410 - t415 * t601;
t397 = pkin(5) * t503 - t616;
t519 = -t597 * t739 + t734;
t748 = t596 * t603;
t476 = t519 * t605 - t601 * t748;
t747 = t596 * t605;
t645 = -t597 * t601 + t602 * t747;
t740 = t601 * t607;
t517 = t597 * t738 + t737;
t754 = t517 * t605;
t632 = g(1) * t476 + g(2) * (t596 * t740 + t754) + g(3) * t645;
t777 = t493 * (pkin(5) * t499 + pkin(10) * t493) + t397 + t632;
t429 = qJDD(6) + t433;
t526 = pkin(5) * t601 - pkin(10) * t605 + t599;
t666 = pkin(5) * t605 + pkin(10) * t601;
t776 = t493 * ((-pkin(4) - t666) * t693 - qJD(5) * t666 - t706) - t526 * t429;
t745 = t596 * t607;
t646 = -t517 * t601 + t605 * t745;
t775 = -t516 * t604 + t600 * t646;
t774 = t516 * t600 + t604 * t646;
t718 = qJD(2) * t602;
t692 = t596 * t718;
t565 = pkin(2) * t692;
t728 = qJ(4) * t692 + t565;
t441 = t596 * t617 + t728;
t770 = pkin(3) + pkin(8);
t680 = t596 * (-pkin(4) - t770);
t569 = t606 * t699;
t588 = t597 * qJD(3);
t726 = t569 + t588;
t451 = t680 * t718 + t726;
t773 = -qJD(5) * t783 - t441 * t601 + t451 * t605;
t772 = -t787 * pkin(3) + t574 * qJ(4) + t576 * qJD(4);
t573 = t576 ^ 2;
t771 = 0.2e1 * t541;
t767 = qJ(3) * t602;
t711 = qJD(6) * t604;
t695 = t604 * t432 - t600 * t503 + t539 * t711;
t712 = qJD(6) * t600;
t404 = -t499 * t712 + t695;
t766 = t404 * t600;
t765 = t429 * t604;
t764 = t456 * t493;
t763 = t456 * t499;
t458 = t499 * t604 + t539 * t600;
t762 = t458 * t493;
t761 = t458 * t499;
t760 = t493 * t600;
t679 = t493 * t604;
t752 = t539 * t605;
t751 = t574 * MDP(8);
t593 = t596 ^ 2;
t750 = t593 * qJD(1) ^ 2;
t744 = t598 * t503;
t743 = t600 * t429;
t742 = t601 * t503;
t741 = t601 * t606;
t736 = t604 * t606;
t735 = t605 * t404;
t550 = qJ(4) * t694;
t570 = pkin(2) * t694;
t462 = -t598 * t693 + t550 + t570;
t674 = (-pkin(3) - pkin(4)) * t749;
t463 = qJD(1) * t674 + t727;
t732 = t605 * t462 + t601 * t463;
t729 = qJ(3) * t572 + qJD(3) * t693;
t725 = -pkin(8) * t749 + t606 * t769;
t594 = t602 ^ 2;
t595 = t606 ^ 2;
t722 = t594 - t595;
t721 = qJ(3) * qJD(2);
t717 = qJD(2) * t606;
t715 = qJD(5) * t598;
t713 = qJD(5) * t605;
t710 = qJD(3) - t463;
t486 = -pkin(3) * t694 + t727;
t709 = qJD(3) - t486;
t701 = g(3) * t749;
t696 = t606 * t750;
t502 = -t597 * pkin(2) - t725;
t691 = t596 * t717;
t690 = 0.2e1 * pkin(1) * t593;
t685 = -t516 * pkin(2) + qJ(3) * t517;
t684 = -t518 * pkin(2) + qJ(3) * t519;
t634 = t601 * t410 + t605 * t415 + t437 * t713 - t449 * t714;
t396 = -pkin(10) * t503 + t634;
t670 = pkin(8) * t787 + t602 * t673 - t606 * t698;
t442 = t670 - t778;
t422 = t442 - t772;
t411 = -pkin(4) * t787 - t422;
t399 = pkin(5) * t433 - pkin(10) * t432 + t411;
t683 = -t600 * t396 + t604 * t399;
t681 = t432 * t600 + t604 * t503;
t630 = t667 - t767;
t466 = t630 * t720;
t482 = t596 * t687 - t724;
t678 = qJD(1) * t482 + t466;
t677 = t576 + t719;
t676 = -qJD(4) - t721;
t675 = t574 + t704;
t672 = t602 * t696;
t669 = t597 * qJ(4) - t502;
t664 = g(1) * t516 - g(2) * t518;
t663 = g(1) * t519 + g(2) * t517;
t662 = g(1) * t517 - g(2) * t519;
t661 = g(1) * t607 + g(2) * t603;
t491 = (-t600 * t602 + t601 * t736) * t720;
t659 = -t604 * t714 - t491;
t658 = t604 * t396 + t600 * t399;
t407 = pkin(10) * t539 + t414;
t436 = -t710 - t781;
t417 = pkin(5) * t497 - pkin(10) * t499 + t436;
t401 = t407 * t604 + t417 * t600;
t657 = t407 * t600 - t417 * t604;
t419 = pkin(10) * t746 + t783;
t454 = t674 + t669;
t515 = t597 * t605 + t601 * t749;
t427 = -pkin(5) * t645 - pkin(10) * t515 + t454;
t656 = t419 * t604 + t427 * t600;
t655 = -t419 * t600 + t427 * t604;
t413 = t437 * t605 - t449 * t601;
t653 = t455 * t605 - t465 * t601;
t649 = -pkin(2) * t606 - pkin(1) - t767;
t489 = t649 * t720;
t501 = -pkin(1) * t596 - t724;
t651 = qJD(2) * (-qJD(1) * t501 - t489);
t650 = -pkin(8) * t692 + t569;
t648 = t607 * pkin(1) + t519 * pkin(2) + pkin(8) * t748 + qJ(3) * t518;
t505 = -qJ(3) * t693 + t570;
t647 = -t515 * t600 + t596 * t736;
t475 = t515 * t604 + t600 * t746;
t644 = -t493 * t711 - t743;
t643 = -t493 * t712 + t765;
t641 = t539 * t601;
t640 = t539 * t458;
t639 = -qJ(3) * t717 - t716;
t621 = t606 * t676 - t716;
t423 = (qJD(1) * t621 + qJDD(1) * t630) * t596 + t730;
t447 = t596 * t621 + t728;
t638 = -qJD(1) * t447 - qJDD(1) * t482 - t423;
t438 = t533 + (qJD(1) * t639 + qJDD(1) * t649) * t596;
t488 = t596 * t639 + t565;
t637 = qJD(1) * t488 + qJDD(1) * t501 + t438;
t467 = -t576 * t693 + t787;
t633 = t605 * t441 + t601 * t451 + t455 * t713 - t465 * t714;
t631 = -pkin(1) * t603 - t517 * pkin(2) + pkin(8) * t745 - qJ(3) * t516;
t509 = t723 * qJD(2);
t629 = t509 * t576 - t662;
t628 = t572 + (-qJD(2) + t576) * t694;
t406 = -pkin(5) * t539 - t413;
t625 = -qJD(3) * t493 - qJD(5) * t406 - t429 * t598;
t622 = -t663 - t671;
t620 = t411 - t627;
t619 = -pkin(10) * t429 + (t406 + t413) * t493;
t618 = -t670 - t627;
t615 = qJD(6) * t493 * t598 + t627;
t587 = t597 * qJD(4);
t452 = t587 + (t606 * t680 - t584) * qJD(2);
t614 = t618 + t778;
t613 = t466 * t693 + t622 + t779;
t612 = t508 * t576 + t618;
t611 = t701 + (-pkin(10) * t694 - qJD(6) * t526 + t732) * t493 + t663;
t609 = t614 + t772;
t496 = t605 * t503;
t492 = -t588 - t650;
t490 = (t600 * t741 + t602 * t604) * t720;
t485 = -t549 - t508;
t484 = t505 + t550;
t481 = -t562 + t708;
t477 = t519 * t601 + t603 * t747;
t473 = qJD(5) * t515 - t605 * t691;
t472 = qJD(5) * t645 + t601 * t691;
t471 = -t587 + (t746 * t770 + t584) * qJD(2);
t470 = -t692 * t770 + t726;
t469 = pkin(3) * t749 - t669;
t468 = t489 * t694;
t461 = t549 + t706;
t444 = t709 + t781;
t440 = t477 * t604 - t518 * t600;
t439 = -t477 * t600 - t518 * t604;
t426 = qJD(6) * t647 + t472 * t604 - t600 * t692;
t425 = qJD(6) * t475 + t472 * t600 + t604 * t692;
t424 = -pkin(3) * t668 + t624;
t420 = pkin(5) * t694 + t462 * t601 - t463 * t605;
t418 = -pkin(5) * t746 - t653;
t416 = pkin(5) * t473 - pkin(10) * t472 + t452;
t405 = qJD(6) * t458 + t681;
t403 = pkin(5) * t692 - t773;
t402 = -pkin(10) * t692 + t633;
t395 = -t401 * qJD(6) + t683;
t394 = -t657 * qJD(6) + t658;
t1 = [(t423 * t482 + t466 * t447 + t422 * t469 + t444 * t471 + t424 * t483 + t461 * t470 - g(1) * (pkin(3) * t745 - qJ(4) * t517 + t631) - g(2) * (pkin(3) * t748 + qJ(4) * t519 + t648)) * MDP(18) + (-t574 * t723 - t576 * t650 + t597 * t671 - t690 * t789 - t664) * MDP(10) + (t404 * t647 - t405 * t475 - t425 * t458 - t426 * t456) * MDP(27) + (-(qJD(6) * t655 + t402 * t604 + t416 * t600) * t493 - t656 * t429 + t394 * t645 - t401 * t473 + t403 * t458 + t418 * t404 + t397 * t475 + t406 * t426 + g(1) * t775 - g(2) * t439) * MDP(32) + (-g(1) * t646 - g(2) * t477 - t411 * t645 + t454 * t433 + t436 * t473 + t452 * t497 - t653 * t503 + t773 * t539) * MDP(24) + ((-qJD(6) * t656 - t402 * t600 + t416 * t604) * t493 + t655 * t429 - t395 * t645 - t657 * t473 + t403 * t456 + t418 * t405 - t397 * t647 + t406 * t425 - g(1) * t774 - g(2) * t440) * MDP(31) + (t405 * t645 - t425 * t493 + t429 * t647 - t456 * t473) * MDP(29) + (-t473 * t539 - t503 * t645) * MDP(22) + (-t429 * t645 + t473 * t493) * MDP(30) + (-t404 * t645 + t426 * t493 + t429 * t475 + t458 * t473) * MDP(28) + (t432 * t645 - t433 * t515 - t472 * t497 - t473 * t499) * MDP(20) + (g(1) * t754 - g(2) * t476 + t411 * t515 + t454 * t432 + t436 * t472 + t452 * t499 + t503 * t783 - t633 * t539) * MDP(25) + (t472 * t539 - t503 * t515) * MDP(21) + (g(1) * t603 - g(2) * t607) * MDP(2) + (-g(1) * t631 - g(2) * t648 + t434 * t500 + t438 * t501 + t442 * t502 + t481 * t509 + t485 * t492 + t489 * t488) * MDP(14) + t661 * MDP(3) + qJDD(1) * MDP(1) + t597 * t751 + (t574 * t725 - t597 * t670 + t635 * t690 - t629) * MDP(9) + (t432 * t515 + t472 * t499) * MDP(19) + (0.2e1 * (t602 * t702 - t705 * t722) * MDP(5) + (qJDD(1) * t594 + 0.2e1 * t602 * t688) * MDP(4)) * t593 + (-t434 * t597 - t492 * t576 - t500 * t574 + t664) * MDP(13) + (t442 * t597 + t502 * t574 + t629) * MDP(12) + (t424 * t597 + t470 * t576 + t483 * t574 + t664) * MDP(16) + (-t422 * t597 - t469 * t574 - t471 * t576 + t662) * MDP(17) + ((g(1) * t740 + t414 * t718 - t606 * t634) * MDP(25) + (-t602 * t637 + t606 * t651) * MDP(13) + (t602 * t651 + t606 * t637) * MDP(12) + (t602 * t675 + t677 * t717) * MDP(6) + (t606 * t675 - t677 * t718) * MDP(7) + (-t503 * t606 - t539 * t718) * MDP(23) + (t602 * t638 - t678 * t717) * MDP(16) + (-t413 * t718 + t606 * t616) * MDP(24) + ((qJD(2) * t481 - qJDD(1) * t500 - t434 + (qJD(2) * t502 - t492) * qJD(1)) * t606 + (qJD(2) * t485 + qJDD(1) * t502 + t442 + (qJD(2) * t500 + t509) * qJD(1)) * t602 - t661) * MDP(11) + ((qJD(2) * t444 + qJDD(1) * t483 + t424 + (qJD(2) * t469 + t470) * qJD(1)) * t606 + (-qJD(2) * t461 + qJDD(1) * t469 + t422 + (-qJD(2) * t483 + t471) * qJD(1)) * t602 - t661) * MDP(15) + (t432 * t606 - t499 * t718) * MDP(21) + (t606 * t638 + t678 * t718) * MDP(17) + (-t433 * t606 + t497 * t718) * MDP(22)) * t596 + (t404 * t475 + t426 * t458) * MDP(26); (-t414 * t694 + t599 * t432 + t707 * t499 + t744 * t601 + t620 * t605 + (t732 + (-qJD(3) - t436) * t601 - t715 * t605) * t539) * MDP(25) + (-t405 * t601 + (t600 * t714 + t490) * t493 + (t644 - t785) * t605) * MDP(29) + (t413 * t694 + t599 * t433 + t707 * t497 + (-t744 + (t436 + t710) * t539) * t605 + ((t462 - t715) * t539 + t620) * t601) * MDP(24) + (t784 * t601 + t788 * t605) * MDP(20) + (-t406 * t491 - t420 * t458 + t776 * t600 + t611 * t604 + (t458 * t715 + t600 * t615 + t604 * t625 - t394) * t601 + (-t401 * t693 - t406 * t712 - qJD(3) * t458 + t397 * t604 - t598 * t404 + (-t598 * t679 - t401) * qJD(5)) * t605) * MDP(32) + (-t486 * t576 + t771 + 0.2e1 * t543 + (-g(3) + (-pkin(3) * qJD(2) + t484) * qJD(1)) * t749 + t613) * MDP(16) + (t771 + t543 + t708 * t576 + (-g(3) * t602 + (t489 * t606 + t505 * t602) * qJD(1)) * t596 + t622) * MDP(13) + (t750 * t768 + t612) * MDP(9) + (t456 * t491 + t458 * t490 + (t456 * t604 + t458 * t600) * t714 + (-t766 - t405 * t604 + (t456 * t600 - t458 * t604) * qJD(6)) * t605) * MDP(27) + (-t422 * t599 + t424 * qJ(3) - t466 * t484 - g(1) * (-qJ(4) * t518 + t684) - g(2) * (-qJ(4) * t516 + t685) - g(3) * (qJ(4) * t746 + t724) + t709 * t461 - t706 * t444) * MDP(18) + (-t539 * t713 + t742 + (-t497 * t602 - t606 * t752) * t720) * MDP(22) + (t429 * t601 + t493 * t752) * MDP(30) + (-t539 * t714 - t496 + (t499 * t602 - t539 * t741) * t720) * MDP(21) + (-t406 * t490 - t420 * t456 - t776 * t604 + t611 * t600 + (t456 * t715 + t600 * t625 - t604 * t615 + t395) * t601 + (-t657 * t693 + t406 * t711 - qJD(3) * t456 + t397 * t600 - t598 * t405 + (-t598 * t760 - t657) * qJD(5)) * t605) * MDP(31) + (pkin(1) * t696 + t576 * t727 - t622 + t701) * MDP(10) + ((-pkin(2) * t703 + ((-pkin(2) * qJD(2) - t481 - t727) * t606 + (-t485 - t508 - t721) * t602) * qJD(1)) * t596 + t729) * MDP(11) + (-t505 * t693 + qJDD(3) + t468 - 0.2e1 * t547 - t612) * MDP(12) - MDP(4) * t672 + (t604 * t735 + (-t605 * t712 + t659) * t458) * MDP(26) + t751 + (t432 * t605 - t499 * t641) * MDP(19) + (t404 * t601 + t659 * t493 + (t640 + t643) * t605) * MDP(28) + t467 * MDP(6) + ((-t599 * t703 + ((-qJD(2) * t599 - t444 - t486) * t606 + (t461 - t487 + t676) * t602) * qJD(1)) * t596 + t729) * MDP(15) + (t609 + t706 * t576 + t574 * t599 + (-t466 * t602 + t484 * t606) * t720) * MDP(17) + (-t442 * pkin(2) - g(1) * t684 - g(2) * t685 - g(3) * t724 - t434 * qJ(3) - t481 * t508 - t485 * t708 - t489 * t505) * MDP(14) + t539 * MDP(23) * t694 + t722 * MDP(5) * t750 + t628 * MDP(7); (t485 * t576 + t468 - t614) * MDP(14) + (-t461 * t576 + t466 * t694 - t609) * MDP(18) + t788 * MDP(24) + t784 * MDP(25) + (t763 - t765) * MDP(31) + (t743 + t761) * MDP(32) + (MDP(12) - MDP(17)) * (t574 + t672) + (MDP(31) * t760 + MDP(32) * t679) * t493 + (MDP(11) + MDP(15)) * t467 + (MDP(16) + MDP(13)) * (-t594 * t750 - t573); t628 * MDP(15) + (t574 - t672) * MDP(16) + (-t595 * t750 - t573) * MDP(17) + (t444 * t576 + (-pkin(3) * t705 - g(3)) * t749 + t613 - t782) * MDP(18) + (-t497 * t576 - t539 * t641 - t496) * MDP(24) + (-t499 * t576 - t539 * t752 + t742) * MDP(25) + (-t605 * t405 + (-t576 * t604 - t600 * t752) * t493 + (t644 + t785) * t601) * MDP(31) + (-t735 + (t576 * t600 - t604 * t752) * t493 + (t640 - t643) * t601) * MDP(32); -t497 ^ 2 * MDP(20) + (t432 + t759) * MDP(21) + (-t433 + t758) * MDP(22) - t503 * MDP(23) + (t414 * t539 + t616 - t632) * MDP(24) + (g(1) * t477 - g(2) * t646 + g(3) * t515 + t413 * t539 + t436 * t497 - t634) * MDP(25) + (t458 * t679 + t766) * MDP(26) + ((t404 - t764) * t604 + (-t405 - t762) * t600) * MDP(27) + (t493 * t679 + t743 - t761) * MDP(28) + (-t493 ^ 2 * t600 + t763 + t765) * MDP(29) + (-pkin(5) * t405 - t414 * t456 + t619 * t600 - t604 * t777) * MDP(31) + (-pkin(5) * t404 - t414 * t458 + t600 * t777 + t619 * t604) * MDP(32) + (MDP(19) * t497 + t499 * MDP(20) - t436 * MDP(24) - t493 * MDP(30) + MDP(31) * t657 + t401 * MDP(32)) * t499; t458 * t456 * MDP(26) + (-t456 ^ 2 + t458 ^ 2) * MDP(27) + (t695 + t764) * MDP(28) + (-t681 + t762) * MDP(29) + t429 * MDP(30) + (-g(1) * t439 - g(2) * t775 - g(3) * t647 + t401 * t493 - t406 * t458 + t683) * MDP(31) + (g(1) * t440 - g(2) * t774 + g(3) * t475 + t406 * t456 - t493 * t657 - t658) * MDP(32) + (-MDP(28) * t757 - MDP(29) * t458 - MDP(31) * t401 + MDP(32) * t657) * qJD(6);];
tau  = t1;
