% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:31
% EndTime: 2019-03-09 11:45:46
% DurationCPUTime: 11.26s
% Computational Cost: add. (12332->603), mult. (29223->726), div. (0->0), fcn. (22460->14), ass. (0->264)
t626 = sin(qJ(5));
t630 = cos(qJ(5));
t623 = sin(pkin(10));
t624 = cos(pkin(10));
t628 = sin(qJ(2));
t631 = cos(qJ(2));
t566 = -t623 * t628 + t624 * t631;
t556 = t566 * qJD(1);
t567 = t623 * t631 + t624 * t628;
t558 = t567 * qJD(1);
t627 = sin(qJ(4));
t772 = cos(qJ(4));
t657 = -t627 * t556 - t558 * t772;
t702 = qJDD(2) + qJDD(4);
t714 = qJD(5) * t626;
t557 = t567 * qJD(2);
t648 = qJD(1) * t557;
t651 = t566 * qJDD(1);
t637 = t651 - t648;
t707 = qJD(1) * qJD(2);
t690 = t631 * t707;
t691 = t628 * t707;
t523 = qJDD(1) * t567 - t623 * t691 + t624 * t690;
t696 = t772 * t523;
t649 = t627 * t558;
t778 = t772 * t556 - t649;
t790 = t778 * qJD(4);
t635 = t627 * t637 + t696 + t790;
t704 = qJD(2) + qJD(4);
t793 = qJD(5) * t704 + t635;
t441 = -t626 * t702 - t630 * t793 - t657 * t714;
t439 = t441 * t626;
t500 = t626 * t704 - t630 * t657;
t509 = qJD(5) - t778;
t642 = t627 * t651 + t696;
t684 = t627 * t523 - t772 * t637;
t791 = t657 * qJD(4);
t457 = t684 - t791;
t456 = qJDD(5) + t457;
t449 = t626 * t456;
t713 = qJD(5) * t630;
t795 = t778 * t630;
t805 = t713 - t795;
t664 = t509 * t805 + t449;
t442 = t626 * t793 - t630 * t702 - t657 * t713;
t498 = -t626 * t657 - t630 * t704;
t665 = -t626 * t442 - t498 * t805;
t711 = t657 * qJD(2);
t745 = t500 * t657;
t803 = t509 * t626;
t785 = -t441 * t630 - t500 * t803;
t807 = (t665 + t785) * MDP(21) + (t664 + t745) * MDP(22) + t702 * MDP(17) + ((-t778 - t649) * qJD(2) + t642) * MDP(15) - t778 ^ 2 * MDP(14) + (MDP(13) * t778 + MDP(14) * t657 + MDP(24) * t509) * t657 + (-t711 - t684) * MDP(16) + (t500 * t805 - t439) * MDP(20);
t450 = t630 * t456;
t804 = -t509 * t803 + t450;
t617 = t631 * pkin(2);
t612 = t617 + pkin(1);
t535 = -pkin(3) * t566 - t612;
t799 = pkin(5) * t657;
t798 = qJ(6) * t657;
t625 = -qJ(3) - pkin(7);
t592 = t625 * t631;
t577 = qJD(1) * t592;
t561 = t623 * t577;
t591 = t625 * t628;
t576 = qJD(1) * t591;
t756 = qJD(2) * pkin(2);
t565 = t576 + t756;
t521 = t624 * t565 + t561;
t767 = pkin(8) * t558;
t493 = qJD(2) * pkin(3) + t521 - t767;
t734 = t624 * t577;
t522 = t623 * t565 - t734;
t768 = pkin(8) * t556;
t497 = t522 + t768;
t460 = t493 * t772 - t627 * t497;
t451 = -pkin(4) * t704 - t460;
t430 = t498 * pkin(5) - t500 * qJ(6) + t451;
t797 = t430 * t778;
t796 = t451 * t778;
t748 = t498 * t657;
t620 = qJ(2) + pkin(10);
t616 = qJ(4) + t620;
t606 = sin(t616);
t632 = cos(qJ(1));
t737 = t606 * t632;
t629 = sin(qJ(1));
t738 = t606 * t629;
t794 = g(1) * t737 + g(2) * t738;
t672 = t630 * pkin(5) + t626 * qJ(6);
t461 = t627 * t493 + t772 * t497;
t452 = pkin(9) * t704 + t461;
t579 = -qJD(1) * t612 + qJD(3);
t529 = -pkin(3) * t556 + t579;
t466 = -pkin(4) * t778 + pkin(9) * t657 + t529;
t427 = -t452 * t626 + t466 * t630;
t708 = qJD(6) - t427;
t420 = -pkin(5) * t509 + t708;
t719 = t794 * t630;
t789 = -t420 * t657 + t430 * t714 + t719;
t428 = t452 * t630 + t466 * t626;
t421 = qJ(6) * t509 + t428;
t607 = cos(t616);
t761 = g(3) * t626;
t699 = -t607 * t761 + t626 * t794;
t686 = qJD(2) * t625;
t553 = -qJD(3) * t628 + t631 * t686;
t520 = qJDD(2) * pkin(2) + qJD(1) * t553 + qJDD(1) * t591;
t552 = qJD(3) * t631 + t628 * t686;
t526 = qJD(1) * t552 - qJDD(1) * t592;
t481 = t624 * t520 - t526 * t623;
t459 = qJDD(2) * pkin(3) - pkin(8) * t523 + t481;
t482 = t623 * t520 + t624 * t526;
t464 = pkin(8) * t637 + t482;
t692 = qJD(4) * t772;
t715 = qJD(4) * t627;
t679 = -t772 * t459 + t627 * t464 + t493 * t715 + t497 * t692;
t415 = -pkin(4) * t702 + t679;
t407 = t442 * pkin(5) + t441 * qJ(6) - t500 * qJD(6) + t415;
t753 = t407 * t626;
t788 = t421 * t657 + t699 - t753;
t787 = t427 * t657 + t451 * t714 + t719;
t786 = t415 * t626 - t428 * t657 + t451 * t713 - t699;
t643 = t627 * t459 + t464 * t772 + t493 * t692 - t497 * t715;
t598 = g(3) * t606;
t735 = t607 * t632;
t736 = t607 * t629;
t698 = g(1) * t735 + g(2) * t736 + t598;
t784 = -t529 * t778 - t643 + t698;
t762 = g(3) * t607;
t783 = t529 * t657 - t679 - t762 + t794;
t671 = pkin(5) * t626 - qJ(6) * t630;
t782 = pkin(5) * t714 - qJ(6) * t713 - t626 * qJD(6) - t671 * t778;
t480 = -pkin(4) * t657 - pkin(9) * t778;
t414 = pkin(9) * t702 + t643;
t703 = pkin(2) * t691 + qJDD(3);
t705 = qJDD(1) * t631;
t754 = qJDD(1) * pkin(1);
t418 = -pkin(2) * t705 - pkin(3) * t637 + t457 * pkin(4) - pkin(9) * t635 + t703 - t754;
t653 = t630 * t414 + t626 * t418 - t452 * t714 + t466 * t713;
t755 = qJ(6) * t456;
t404 = qJD(6) * t509 + t653 + t755;
t680 = t626 * t414 - t630 * t418 + t452 * t713 + t466 * t714;
t769 = pkin(5) * t456;
t406 = qJDD(6) + t680 - t769;
t780 = t404 * t630 + t406 * t626;
t525 = t627 * t566 + t567 * t772;
t656 = t566 * t772 - t627 * t567;
t475 = -pkin(4) * t656 - pkin(9) * t525 + t535;
t531 = t624 * t591 + t592 * t623;
t507 = -pkin(8) * t567 + t531;
t532 = t623 * t591 - t624 * t592;
t508 = pkin(8) * t566 + t532;
t477 = t627 * t507 + t508 * t772;
t724 = t626 * t475 + t630 * t477;
t527 = -t576 * t623 + t734;
t502 = t527 - t768;
t528 = t624 * t576 + t561;
t503 = t528 - t767;
t469 = t627 * t502 + t503 * t772;
t608 = pkin(2) * t624 + pkin(3);
t771 = pkin(2) * t623;
t673 = t772 * t608 - t627 * t771;
t539 = t673 * qJD(4);
t779 = -t539 + t469;
t718 = t627 * t608 + t772 * t771;
t721 = t718 * qJD(4) + t502 * t772 - t627 * t503;
t777 = t607 * pkin(4) + t606 * pkin(9);
t776 = g(1) * t629 - g(2) * t632;
t775 = t509 * t714 - t450;
t675 = g(1) * t632 + g(2) * t629;
t774 = t500 ^ 2;
t773 = t509 ^ 2;
t766 = pkin(9) * t456;
t760 = g(3) * t631;
t759 = t628 * pkin(2);
t757 = pkin(9) * qJD(5);
t752 = t421 * t509;
t751 = t428 * t509;
t750 = t442 * t630;
t551 = pkin(9) + t718;
t749 = t456 * t551;
t747 = t498 * t626;
t746 = t500 * t498;
t744 = t500 * t630;
t740 = t525 * t630;
t739 = t539 * t630;
t732 = t626 * t629;
t731 = t629 * t630;
t730 = t630 * t632;
t729 = t632 * t626;
t726 = t630 * t460 + t626 * t480;
t533 = pkin(3) * t558 + qJD(1) * t759;
t470 = t480 + t533;
t725 = t630 * t469 + t626 * t470;
t722 = t721 + t782;
t506 = t624 * t552 + t623 * t553;
t720 = -t461 + t782;
t717 = pkin(3) * cos(t620) + t617;
t621 = t628 ^ 2;
t716 = -t631 ^ 2 + t621;
t706 = qJDD(1) * t628;
t615 = t628 * t756;
t693 = t551 * t714;
t534 = pkin(3) * t557 + t615;
t688 = -t407 - t762;
t687 = -t415 - t762;
t424 = t469 * t626 - t470 * t630 + t799;
t685 = t539 * t626 - t424;
t505 = -t552 * t623 + t624 * t553;
t678 = t607 * t672 + t777;
t546 = t607 * t732 + t730;
t548 = t607 * t729 - t731;
t677 = -g(1) * t546 + g(2) * t548;
t547 = t607 * t731 - t729;
t549 = t607 * t730 + t732;
t676 = g(1) * t547 - g(2) * t549;
t669 = t420 * t630 - t421 * t626;
t668 = -t749 - t796;
t667 = pkin(4) + t672;
t666 = pkin(1) + t717 + t777;
t662 = t778 * t803 - t775;
t661 = t420 * t713 - t698 + t780;
t660 = -0.2e1 * pkin(1) * t707 - pkin(7) * qJDD(2);
t659 = t507 * t772 - t627 * t508;
t560 = t566 * qJD(2);
t483 = qJD(4) * t656 - t627 * t557 + t560 * t772;
t655 = t483 * t626 + t525 * t713;
t654 = -t483 * t630 + t525 * t714;
t486 = -pkin(8) * t560 + t505;
t487 = -pkin(8) * t557 + t506;
t434 = qJD(4) * t659 + t627 * t486 + t487 * t772;
t484 = qJD(4) * t525 + t557 * t772 + t627 * t560;
t444 = pkin(4) * t484 - pkin(9) * t483 + t534;
t652 = t630 * t434 + t626 * t444 + t475 * t713 - t477 * t714;
t650 = -qJDD(1) * t612 + t703;
t633 = qJD(2) ^ 2;
t646 = -pkin(7) * t633 + 0.2e1 * t754 + t776;
t634 = qJD(1) ^ 2;
t645 = pkin(1) * t634 - pkin(7) * qJDD(1) + t675;
t644 = g(1) * t548 + g(2) * t546 + t606 * t761 - t680;
t641 = qJD(5) * t669 + t780;
t640 = t430 * t500 + qJDD(6) - t644;
t639 = -g(1) * t549 - g(2) * t547 - t598 * t630 + t653;
t435 = qJD(4) * t477 - t486 * t772 + t627 * t487;
t636 = t675 * t667 * t606;
t619 = -pkin(8) + t625;
t584 = pkin(9) * t735;
t582 = pkin(9) * t736;
t578 = -pkin(3) * sin(t620) - t759;
t550 = -pkin(4) - t673;
t530 = -t667 - t673;
t496 = pkin(3) * t648 + qJDD(1) * t535 + t703;
t471 = pkin(5) * t500 + qJ(6) * t498;
t445 = t525 * t671 - t659;
t432 = pkin(5) * t656 - t475 * t630 + t477 * t626;
t431 = -qJ(6) * t656 + t724;
t426 = t460 * t626 - t480 * t630 + t799;
t425 = t726 - t798;
t423 = t725 - t798;
t422 = t498 * t509 - t441;
t410 = t671 * t483 + (qJD(5) * t672 - qJD(6) * t630) * t525 + t435;
t409 = -pkin(5) * t484 + qJD(5) * t724 + t434 * t626 - t444 * t630;
t408 = qJ(6) * t484 - qJD(6) * t656 + t652;
t1 = [t675 * MDP(3) + (-t481 * t567 + t482 * t566 - t505 * t558 + t506 * t556 - t521 * t560 - t522 * t557 - t531 * t523 + t532 * t637 - t675) * MDP(11) + (t404 * t431 + t421 * t408 + t407 * t445 + t430 * t410 + t406 * t432 + t420 * t409 - g(1) * (-pkin(5) * t547 - qJ(6) * t546) - g(2) * (pkin(5) * t549 + qJ(6) * t548) + (g(1) * t619 - g(2) * t666) * t632 + (g(1) * t666 + g(2) * t619) * t629) * MDP(30) + (t628 * t660 + t631 * t646) * MDP(9) + (-t628 * t646 + t631 * t660) * MDP(10) + (-t483 * t657 + t525 * t635) * MDP(13) + (-g(1) * t738 + g(2) * t737 - t434 * t704 - t477 * t702 + t529 * t483 + t496 * t525 - t534 * t657 + t535 * t635) * MDP(19) + (t442 * t656 - t449 * t525 - t484 * t498 - t509 * t655) * MDP(23) + (-t404 * t656 - t407 * t740 + t408 * t509 - t410 * t500 + t421 * t484 + t430 * t654 + t431 * t456 + t441 * t445 - t677) * MDP(29) + (t406 * t656 - t409 * t509 + t410 * t498 - t420 * t484 + t430 * t655 - t432 * t456 + t442 * t445 + t525 * t753 + t676) * MDP(27) + (-t456 * t656 + t484 * t509) * MDP(24) + (-t525 * t457 + t483 * t778 + t484 * t657 + t635 * t656) * MDP(14) + (t680 * t656 + t427 * t484 + t435 * t498 - t659 * t442 + ((-qJD(5) * t477 + t444) * t509 + t475 * t456 + t451 * qJD(5) * t525) * t630 + ((-qJD(5) * t475 - t434) * t509 - t477 * t456 + t415 * t525 + t451 * t483) * t626 + t676) * MDP(25) + (t415 * t740 - t428 * t484 + t435 * t500 + t441 * t659 - t451 * t654 - t456 * t724 - t509 * t652 + t653 * t656 + t677) * MDP(26) + (-t484 * t704 + t656 * t702) * MDP(16) + (t441 * t656 + t450 * t525 + t484 * t500 - t509 * t654) * MDP(22) + (-t408 * t498 + t409 * t500 - t431 * t442 - t432 * t441 + t776 * t606 + t669 * t483 + (-t404 * t626 + t406 * t630 + (-t420 * t626 - t421 * t630) * qJD(5)) * t525) * MDP(28) + (-t435 * t704 + t535 * t457 + t529 * t484 - t496 * t656 - t534 * t778 + t607 * t776 + t659 * t702) * MDP(18) + t776 * MDP(2) + qJDD(1) * MDP(1) + (qJDD(2) * t628 + t631 * t633) * MDP(6) + (qJDD(2) * t631 - t628 * t633) * MDP(7) + (qJDD(1) * t621 + 0.2e1 * t628 * t690) * MDP(4) + (t482 * t532 + t522 * t506 + t481 * t531 + t521 * t505 - t650 * t612 + t579 * t615 - g(1) * (-t612 * t629 - t625 * t632) - g(2) * (t612 * t632 - t625 * t629)) * MDP(12) + (t483 * t704 + t525 * t702) * MDP(15) + 0.2e1 * (t628 * t705 - t707 * t716) * MDP(5) + (-t441 * t740 - t500 * t654) * MDP(20) + ((-t498 * t630 - t500 * t626) * t483 + (t439 - t750 + (-t744 + t747) * qJD(5)) * t525) * MDP(21); ((t522 + t527) * t558 + (-t528 + t521) * t556 + (-t624 * t523 + ((-t690 - t706) * t623 + (-t691 + t705) * t624) * t623) * pkin(2)) * MDP(11) + (t662 - t748) * MDP(23) + (g(3) * t628 + t631 * t645) * MDP(10) + (t423 * t498 - t424 * t500 + (-t420 * t778 - t498 * t539 + (qJD(5) * t500 - t442) * t551) * t630 + (t421 * t778 - t441 * t551 + t500 * t539 + (t498 * t551 - t421) * qJD(5)) * t626 + t661) * MDP(28) + qJDD(2) * MDP(8) + (t533 * t778 + t673 * t702 - t704 * t721 + t783) * MDP(18) + (t533 * t657 - t702 * t718 + t704 * t779 + t784) * MDP(19) + (-t550 * t441 + t668 * t630 + t721 * t500 + (t693 + t725 - t739) * t509 + t786) * MDP(26) + (t550 * t442 + t687 * t630 + t668 * t626 + t721 * t498 + ((-qJD(5) * t551 - t470) * t630 + t779 * t626) * t509 + t787) * MDP(25) + (t441 * t530 - t722 * t500 + t749 * t630 + (-t423 - t693 + (-t430 + t539) * t630) * t509 + t788) * MDP(29) + (-MDP(4) * t628 * t631 + MDP(5) * t716) * t634 + (t442 * t530 + t688 * t630 + (-t749 - t797) * t626 + t722 * t498 + (-t551 * t713 - t685) * t509 + t789) * MDP(27) + MDP(7) * t705 + MDP(6) * t706 + (t407 * t530 - g(1) * (t578 * t632 + t584) - g(2) * (t578 * t629 + t582) - g(3) * (t678 + t717) + t722 * t430 + (-t423 + t739) * t421 + t685 * t420 + t636 + t641 * t551) * MDP(30) + (-t521 * t527 - t522 * t528 + (-t760 + t481 * t624 + t482 * t623 + (-qJD(1) * t579 + t675) * t628) * pkin(2)) * MDP(12) + (t628 * t645 - t760) * MDP(9) + t807; (-t556 ^ 2 - t558 ^ 2) * MDP(11) + (t521 * t558 - t522 * t556 + t650 - t776) * MDP(12) + (-t711 + t684 - 0.2e1 * t791) * MDP(18) + ((t778 - t649) * qJD(2) + 0.2e1 * t790 + t642) * MDP(19) + (t662 + t748) * MDP(25) + (-t630 * t773 - t449 + t745) * MDP(26) + (t748 + t804) * MDP(27) + (t665 - t785) * MDP(28) + (t664 - t745) * MDP(29) + (t430 * t657 + (-t406 + t752) * t630 + (t420 * t509 + t404) * t626 - t776) * MDP(30); (t461 * t704 + t783) * MDP(18) + (t460 * t704 + t784) * MDP(19) + (-t748 + t804) * MDP(23) + (-pkin(4) * t442 - t461 * t498 + (t460 * t509 - t766 - t796) * t626 + ((-t480 - t757) * t509 + t687) * t630 + t787) * MDP(25) + (pkin(4) * t441 + pkin(9) * t775 - t451 * t795 - t461 * t500 + t726 * t509 + t786) * MDP(26) + (t426 * t509 - t442 * t667 + (-t766 - t797) * t626 + t720 * t498 + (-t509 * t757 + t688) * t630 + t789) * MDP(27) + (-t420 * t795 + t425 * t498 - t426 * t500 - t421 * t803 + (-t439 - t750 + (t744 + t747) * qJD(5)) * pkin(9) + t661) * MDP(28) + (-t441 * t667 + (-pkin(9) * t714 - t425) * t509 - t720 * t500 + (-t430 * t509 + t766) * t630 + t788) * MDP(29) + (pkin(9) * t641 - g(1) * t584 - g(2) * t582 - g(3) * t678 - t407 * t667 - t420 * t426 - t421 * t425 + t430 * t720 + t636) * MDP(30) + t807; MDP(20) * t746 + (-t498 ^ 2 + t774) * MDP(21) + t422 * MDP(22) + (t500 * t509 - t442) * MDP(23) + t456 * MDP(24) + (-t451 * t500 + t644 + t751) * MDP(25) + (t427 * t509 + t451 * t498 - t639) * MDP(26) + (-t471 * t498 - t640 + t751 + 0.2e1 * t769) * MDP(27) + (pkin(5) * t441 - qJ(6) * t442 + (t421 - t428) * t500 + (t420 - t708) * t498) * MDP(28) + (0.2e1 * t755 - t430 * t498 + t471 * t500 + (0.2e1 * qJD(6) - t427) * t509 + t639) * MDP(29) + (t404 * qJ(6) - t406 * pkin(5) - t430 * t471 - t420 * t428 - g(1) * (-pkin(5) * t548 + qJ(6) * t549) - g(2) * (-pkin(5) * t546 + qJ(6) * t547) + t671 * t598 + t708 * t421) * MDP(30); (t746 - t456) * MDP(27) + t422 * MDP(28) + (-t773 - t774) * MDP(29) + (t640 - t752 - t769) * MDP(30);];
tau  = t1;
