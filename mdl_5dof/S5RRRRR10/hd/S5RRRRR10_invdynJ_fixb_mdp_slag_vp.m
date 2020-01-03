% Calculate vector of inverse dynamics joint torques for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:36:21
% EndTime: 2019-12-31 22:36:47
% DurationCPUTime: 14.76s
% Computational Cost: add. (8383->595), mult. (20415->830), div. (0->0), fcn. (16565->14), ass. (0->262)
t629 = sin(pkin(5));
t639 = cos(qJ(2));
t736 = qJD(1) * t639;
t708 = t629 * t736;
t598 = -qJD(3) + t708;
t817 = qJD(4) - t598;
t636 = cos(qJ(5));
t728 = qJD(5) * t636;
t630 = cos(pkin(5));
t737 = qJD(1) * t630;
t614 = qJD(2) + t737;
t638 = cos(qJ(3));
t633 = sin(qJ(3));
t634 = sin(qJ(2));
t738 = qJD(1) * t629;
t709 = t634 * t738;
t683 = t633 * t709;
t549 = t614 * t638 - t683;
t550 = t614 * t633 + t638 * t709;
t632 = sin(qJ(4));
t637 = cos(qJ(4));
t497 = -t637 * t549 + t550 * t632;
t812 = t497 * t636;
t819 = t728 + t812;
t725 = qJD(1) * qJD(2);
t703 = t639 * t725;
t723 = qJDD(1) * t634;
t818 = t703 + t723;
t677 = pkin(2) * t634 - pkin(8) * t639;
t566 = t677 * t738;
t554 = t638 * t566;
t720 = pkin(1) * t737;
t565 = -pkin(7) * t709 + t639 * t720;
t789 = pkin(8) + pkin(9);
t710 = qJD(3) * t789;
t754 = t638 * t639;
t816 = -t565 * t633 + t554 + (pkin(3) * t634 - pkin(9) * t754) * t738 + t638 * t710;
t684 = t633 * t708;
t743 = t638 * t565 + t633 * t566;
t815 = pkin(9) * t684 - t633 * t710 - t743;
t663 = t549 * t632 + t637 * t550;
t724 = qJDD(1) * t630;
t613 = qJDD(2) + t724;
t732 = qJD(3) * t638;
t809 = t818 * t629;
t486 = -qJD(3) * t683 + t633 * t613 + t614 * t732 + t638 * t809;
t734 = qJD(2) * t639;
t705 = t633 * t734;
t733 = qJD(3) * t633;
t487 = -t638 * t613 + t614 * t733 + t629 * (qJD(1) * (t634 * t732 + t705) + t633 * t723);
t730 = qJD(4) * t637;
t731 = qJD(4) * t632;
t440 = t637 * t486 - t632 * t487 + t549 * t730 - t550 * t731;
t722 = qJDD(1) * t639;
t612 = t629 * t722;
t704 = t634 * t725;
t682 = t629 * t704;
t563 = qJDD(3) - t612 + t682;
t557 = qJDD(4) + t563;
t631 = sin(qJ(5));
t712 = t636 * t440 + t631 * t557 + t728 * t817;
t729 = qJD(5) * t631;
t420 = -t663 * t729 + t712;
t418 = t420 * t631;
t473 = t631 * t817 + t636 * t663;
t698 = t440 * t631 - t636 * t557;
t421 = qJD(5) * t473 + t698;
t441 = qJD(4) * t663 + t486 * t632 + t637 * t487;
t437 = qJDD(5) + t441;
t434 = t631 * t437;
t435 = t636 * t437;
t777 = t663 * t631;
t471 = -t636 * t817 + t777;
t726 = -qJD(5) - t497;
t811 = t726 * t631;
t814 = t557 * MDP(22) - t441 * MDP(21) - t497 ^ 2 * MDP(19) + (t497 * t817 + t440) * MDP(20) + (MDP(18) * t497 + t663 * MDP(19) + MDP(21) * t817 + t726 * MDP(29)) * t663 + (t473 * t819 + t418) * MDP(25) + (-t473 * t663 - t726 * t819 + t434) * MDP(27) + (t471 * t663 - t726 * t811 + t435) * MDP(28) + (t420 * t636 - t631 * t421 - t471 * t819 + t473 * t811) * MDP(26);
t568 = pkin(7) * t708 + t634 * t720;
t532 = pkin(8) * t614 + t568;
t661 = -pkin(2) * t639 - pkin(8) * t634 - pkin(1);
t562 = t661 * t629;
t544 = qJD(1) * t562;
t483 = -t532 * t633 + t638 * t544;
t463 = -pkin(9) * t550 + t483;
t460 = -pkin(3) * t598 + t463;
t484 = t532 * t638 + t544 * t633;
t464 = pkin(9) * t549 + t484;
t781 = t464 * t632;
t432 = t460 * t637 - t781;
t430 = -pkin(4) * t817 - t432;
t813 = t430 * t497;
t579 = t632 * t633 - t637 * t638;
t810 = t817 * t579;
t580 = t632 * t638 + t633 * t637;
t745 = t817 * t580;
t640 = cos(qJ(1));
t752 = t639 * t640;
t635 = sin(qJ(1));
t758 = t634 * t635;
t577 = -t630 * t758 + t752;
t628 = qJ(3) + qJ(4);
t623 = sin(t628);
t624 = cos(t628);
t764 = t629 * t635;
t519 = -t577 * t623 + t624 * t764;
t761 = t629 * t640;
t765 = t629 * t634;
t756 = t635 * t639;
t757 = t634 * t640;
t575 = t630 * t757 + t756;
t772 = t575 * t623;
t808 = g(3) * (-t623 * t765 + t624 * t630) + g(2) * (-t624 * t761 - t772) + g(1) * t519;
t458 = pkin(4) * t663 + pkin(10) * t497;
t531 = -pkin(2) * t614 - t565;
t501 = -pkin(3) * t549 + t531;
t517 = t575 * t624 - t623 * t761;
t520 = t577 * t624 + t623 * t764;
t559 = t623 * t630 + t624 * t765;
t719 = pkin(1) * qJD(2) * t630;
t687 = qJD(1) * t719;
t717 = pkin(1) * t724;
t711 = -pkin(7) * t612 - t634 * t717 - t639 * t687;
t650 = -pkin(7) * t682 - t711;
t505 = pkin(8) * t613 + t650;
t660 = t677 * qJD(2);
t507 = (qJD(1) * t660 + qJDD(1) * t661) * t629;
t645 = -qJD(3) * t484 - t633 * t505 + t638 * t507;
t425 = pkin(3) * t563 - pkin(9) * t486 + t645;
t659 = -t638 * t505 - t633 * t507 + t532 * t733 - t544 * t732;
t427 = -pkin(9) * t487 - t659;
t686 = -t632 * t425 - t637 * t427 - t460 * t730 + t464 * t731;
t807 = g(1) * t520 + g(2) * t517 + g(3) * t559 + t497 * t501 + t686;
t574 = -t630 * t752 + t758;
t806 = t517 * t631 - t574 * t636;
t805 = t517 * t636 + t574 * t631;
t763 = t629 * t638;
t572 = t630 * t633 + t634 * t763;
t762 = t629 * t639;
t788 = pkin(1) * t634;
t742 = pkin(7) * t762 + t630 * t788;
t561 = pkin(8) * t630 + t742;
t694 = -t561 * t633 + t638 * t562;
t469 = -pkin(3) * t762 - pkin(9) * t572 + t694;
t571 = -t630 * t638 + t633 * t765;
t744 = t638 * t561 + t633 * t562;
t476 = -pkin(9) * t571 + t744;
t801 = t632 * t469 + t637 * t476;
t599 = t789 * t633;
t600 = t789 * t638;
t662 = -t599 * t637 - t600 * t632;
t800 = qJD(4) * t662 - t632 * t816 + t815 * t637;
t537 = -t599 * t632 + t600 * t637;
t799 = qJD(4) * t537 + t815 * t632 + t637 * t816;
t680 = -t568 + (-t684 + t733) * pkin(3);
t780 = t464 * t637;
t433 = t460 * t632 + t780;
t431 = pkin(10) * t817 + t433;
t444 = pkin(4) * t497 - pkin(10) * t663 + t501;
t670 = t431 * t631 - t444 * t636;
t797 = t430 * t729 + t663 * t670;
t647 = -qJD(4) * t433 + t637 * t425 - t427 * t632;
t410 = -pkin(4) * t557 - t647;
t408 = t410 * t631;
t416 = t431 * t636 + t444 * t631;
t796 = t416 * t663 + t430 * t728 + t631 * t808 + t408;
t793 = -t663 * t501 + t647 - t808;
t706 = t629 * t734;
t525 = -qJD(3) * t571 + t638 * t706;
t567 = t629 * t660;
t615 = pkin(7) * t765;
t760 = t630 * t639;
t569 = (t760 * pkin(1) - t615) * qJD(2);
t646 = -qJD(3) * t744 + t638 * t567 - t569 * t633;
t735 = qJD(2) * t634;
t707 = t629 * t735;
t450 = pkin(3) * t707 - pkin(9) * t525 + t646;
t524 = qJD(3) * t572 + t629 * t705;
t658 = -t561 * t733 + t562 * t732 + t633 * t567 + t638 * t569;
t453 = -pkin(9) * t524 + t658;
t791 = -qJD(4) * t801 + t450 * t637 - t453 * t632;
t786 = g(1) * t640;
t783 = g(3) * t629;
t782 = MDP(6) * t629;
t776 = t549 * t598;
t775 = t550 * t598;
t771 = t575 * t633;
t770 = t580 * t636;
t769 = t613 * MDP(8);
t768 = t624 * t631;
t767 = t624 * t636;
t625 = t629 ^ 2;
t766 = t625 * qJD(1) ^ 2;
t759 = t631 * t639;
t755 = t636 * t639;
t753 = t638 * t640;
t747 = pkin(4) * t709 + t799;
t570 = pkin(7) * t706 + t634 * t719;
t626 = t634 ^ 2;
t741 = -t639 ^ 2 + t626;
t727 = qJD(2) - t614;
t721 = 0.2e1 * t625;
t716 = t639 * t766;
t715 = t629 * t759;
t714 = t629 * t755;
t622 = -pkin(3) * t638 - pkin(2);
t409 = pkin(10) * t557 - t686;
t685 = pkin(7) * t809 + t634 * t687 - t639 * t717;
t506 = -pkin(2) * t613 + t685;
t457 = pkin(3) * t487 + t506;
t412 = pkin(4) * t441 - pkin(10) * t440 + t457;
t700 = -t631 * t409 + t636 * t412;
t696 = t631 * t810 - t636 * t709;
t695 = t631 * t709 + t636 * t810;
t693 = t575 * t638 - t633 * t761;
t690 = t614 + t737;
t620 = pkin(3) * t632 + pkin(10);
t689 = pkin(3) * t550 + qJD(5) * t620 + t458;
t688 = t613 + t724;
t502 = pkin(3) * t524 + t570;
t442 = t463 * t632 + t780;
t679 = pkin(3) * t731 - t442;
t443 = t463 * t637 - t781;
t678 = -pkin(3) * t730 + t443;
t515 = pkin(4) * t579 - pkin(10) * t580 + t622;
t675 = pkin(10) * t709 - qJD(5) * t515 - t800;
t674 = -pkin(4) * t745 - pkin(10) * t810 + qJD(5) * t537 - t680;
t672 = t636 * t409 + t631 * t412;
t671 = -t620 * t437 + t813;
t446 = -pkin(10) * t762 + t801;
t510 = t637 * t571 + t572 * t632;
t511 = -t571 * t632 + t572 * t637;
t560 = t615 + (-pkin(1) * t639 - pkin(2)) * t630;
t512 = pkin(3) * t571 + t560;
t454 = pkin(4) * t510 - pkin(10) * t511 + t512;
t669 = t446 * t636 + t454 * t631;
t668 = -t446 * t631 + t454 * t636;
t666 = t469 * t637 - t476 * t632;
t492 = t511 * t631 + t714;
t657 = t632 * t450 + t637 * t453 + t469 * t730 - t476 * t731;
t655 = t580 * t728 - t696;
t654 = -t580 * t729 - t695;
t576 = t630 * t756 + t757;
t652 = -g(1) * t576 - g(2) * t574 + g(3) * t762;
t651 = -t410 - t808;
t649 = -pkin(8) * t563 - t531 * t598;
t643 = pkin(8) * qJD(3) * t598 - t506 - t652;
t621 = -pkin(3) * t637 - pkin(4);
t528 = t577 * t638 + t633 * t764;
t527 = -t577 * t633 + t635 * t763;
t493 = t511 * t636 - t715;
t490 = t520 * t636 + t576 * t631;
t489 = -t520 * t631 + t576 * t636;
t456 = qJD(4) * t511 + t637 * t524 + t525 * t632;
t455 = -qJD(4) * t510 - t524 * t632 + t525 * t637;
t445 = pkin(4) * t762 - t666;
t439 = -qJD(5) * t715 + t455 * t631 + t511 * t728 - t636 * t707;
t438 = -qJD(5) * t492 + t455 * t636 + t631 * t707;
t422 = pkin(4) * t456 - pkin(10) * t455 + t502;
t414 = -pkin(4) * t707 - t791;
t413 = pkin(10) * t707 + t657;
t407 = -t416 * qJD(5) + t700;
t406 = -t670 * qJD(5) + t672;
t1 = [(t634 * t688 + t690 * t734) * t782 + t630 * t769 + qJDD(1) * MDP(1) + (g(2) * t635 + t786) * MDP(3) + (-t570 * t614 - t615 * t613 - t685 * t630 + g(1) * t575 - g(2) * t577 + (t613 * t760 + (-t704 + t722) * t721) * pkin(1)) * MDP(9) + (-t525 * t598 + t563 * t572) * MDP(13) + (t524 * t598 - t563 * t571) * MDP(14) + (-t456 * t817 - t510 * t557) * MDP(21) + (t455 * t817 + t511 * t557) * MDP(20) + ((-t433 * t735 - t624 * t786 - t639 * t686) * MDP(24) + (-g(1) * t753 - t484 * t735 - t639 * t659) * MDP(17) + (t441 * t639 - t497 * t735) * MDP(21) + (t483 * t735 - t639 * t645) * MDP(16) + (t432 * t735 - t639 * t647) * MDP(23) + (-t486 * t639 + t550 * t735) * MDP(13) + (t487 * t639 + t549 * t735) * MDP(14) + (-t440 * t639 + t663 * t735) * MDP(20) + (t639 * t688 - t690 * t735) * MDP(7) + (-t563 * t639 - t598 * t735) * MDP(15) + (-t557 * t639 + t735 * t817) * MDP(22)) * t629 + (g(1) * t517 - g(2) * t520 + t512 * t441 + t501 * t456 + t457 * t510 + t502 * t497 + t666 * t557 + t791 * t817) * MDP(23) + (-g(1) * t772 - g(2) * t519 + t512 * t440 + t501 * t455 + t457 * t511 + t502 * t663 - t557 * t801 - t657 * t817) * MDP(24) + (-pkin(1) * t721 * t818 - g(1) * t574 + g(2) * t576 - t569 * t614 - t742 * t613 - t650 * t630) * MDP(10) + (t440 * t511 + t455 * t663) * MDP(18) + (-t440 * t510 - t441 * t511 - t455 * t497 - t456 * t663) * MDP(19) + (g(1) * t635 - g(2) * t640) * MDP(2) + (-t486 * t571 - t487 * t572 - t524 * t550 + t525 * t549) * MDP(12) + (t486 * t572 + t525 * t550) * MDP(11) + (-(-qJD(5) * t669 - t413 * t631 + t422 * t636) * t726 + t668 * t437 + t407 * t510 - t670 * t456 + t414 * t471 + t445 * t421 + t410 * t492 + t430 * t439 + g(1) * t805 - g(2) * t490) * MDP(30) + ((qJD(5) * t668 + t413 * t636 + t422 * t631) * t726 - t669 * t437 - t406 * t510 - t416 * t456 + t414 * t473 + t445 * t420 + t410 * t493 + t430 * t438 - g(1) * t806 - g(2) * t489) * MDP(31) + (t437 * t510 - t456 * t726) * MDP(29) + (-t421 * t510 - t437 * t492 + t439 * t726 - t456 * t471) * MDP(28) + (t420 * t510 + t437 * t493 - t438 * t726 + t456 * t473) * MDP(27) + ((qJDD(1) * t626 + 0.2e1 * t634 * t703) * MDP(4) + 0.2e1 * (t634 * t722 - t725 * t741) * MDP(5)) * t625 + (-g(1) * t771 - g(2) * t527 + t560 * t486 + t506 * t572 + t531 * t525 + t570 * t550 - t744 * t563 + t658 * t598) * MDP(17) + (g(1) * t693 - g(2) * t528 + t560 * t487 + t506 * t571 + t531 * t524 - t570 * t549 + t694 * t563 - t646 * t598) * MDP(16) + (t420 * t493 + t438 * t473) * MDP(25) + (-t420 * t492 - t421 * t493 - t438 * t471 - t439 * t473) * MDP(26); (t727 * t736 + t723) * t782 + (-pkin(2) * t486 - t568 * t550 - t598 * t743 - t633 * t643 + t638 * t649) * MDP(17) + (-pkin(2) * t487 + t568 * t549 + t554 * t598 + (-t565 * t598 + t649) * t633 + t643 * t638) * MDP(16) + (-t557 * t579 - t745 * t817) * MDP(21) + (t598 * MDP(15) - t483 * MDP(16) + t484 * MDP(17) - MDP(20) * t663 + t497 * MDP(21) - MDP(22) * t817 - t432 * MDP(23) + t433 * MDP(24) - MDP(7) * t727) * t709 + (t557 * t580 - t810 * t817) * MDP(20) + (t622 * t440 + t457 * t580 - t501 * t810 - t537 * t557 + t652 * t623 + t680 * t663 - t800 * t817) * MDP(24) + (t622 * t441 + t457 * t579 + t680 * t497 + t745 * t501 + t557 * t662 - t652 * t624 - t799 * t817) * MDP(23) + t769 + t741 * MDP(5) * t766 + (t568 * t614 + t766 * t788 - t652 - t685) * MDP(9) + (t696 * t473 + t695 * t471 + (-t418 - t421 * t636 + (t471 * t631 - t473 * t636) * qJD(5)) * t580) * MDP(26) + (t486 * t633 - t638 * t775) * MDP(11) + ((t486 - t776) * t638 + (-t487 + t775) * t633) * MDP(12) + (t420 * t770 + t473 * t654) * MDP(25) + (pkin(1) * t716 + g(1) * t577 + g(2) * t575 + t565 * t614 + (pkin(7) * t725 + g(3)) * t765 + t711) * MDP(10) + t612 * MDP(7) + (t440 * t580 - t663 * t810) * MDP(18) + (-t440 * t579 - t441 * t580 + t497 * t810 - t663 * t745) * MDP(19) + (-t598 * t732 + t563 * t633 + (-t550 * t634 + t598 * t754) * t738) * MDP(13) + (t598 * t733 + t563 * t638 + (-t598 * t633 * t639 - t549 * t634) * t738) * MDP(14) - t634 * MDP(4) * t716 + (-t421 * t579 - t434 * t580 - t471 * t745 + t655 * t726) * MDP(28) + (t420 * t579 + t435 * t580 + t473 * t745 - t654 * t726) * MDP(27) + (t437 * t579 - t726 * t745) * MDP(29) + (-(t515 * t631 + t537 * t636) * t437 - t406 * t579 - t662 * t420 + t410 * t770 - g(1) * (t576 * t768 + t577 * t636) - g(2) * (t574 * t768 + t575 * t636) - (-t624 * t759 + t634 * t636) * t783 - (t631 * t674 + t636 * t675) * t726 + t747 * t473 - t745 * t416 + t654 * t430) * MDP(31) + ((t515 * t636 - t537 * t631) * t437 + t407 * t579 - t662 * t421 + t580 * t408 - g(1) * (-t576 * t767 + t577 * t631) - g(2) * (-t574 * t767 + t575 * t631) - (t624 * t755 + t631 * t634) * t783 - (t631 * t675 - t636 * t674) * t726 + t747 * t471 - t745 * t670 + t655 * t430) * MDP(30); (-t487 - t775) * MDP(14) + (t486 + t776) * MDP(13) + (-t484 * t598 - t531 * t550 - g(1) * t527 - g(2) * (-t629 * t753 - t771) + g(3) * t571 + t645) * MDP(16) + (t442 * t817 + (-t497 * t550 + t557 * t637 - t731 * t817) * pkin(3) + t793) * MDP(23) + (g(1) * t528 + g(2) * t693 + g(3) * t572 - t483 * t598 - t531 * t549 + t659) * MDP(17) + (t621 * t421 + t679 * t471 + (-t678 * t726 + t671) * t631 + (t689 * t726 + t651) * t636 + t797) * MDP(30) + (t621 * t420 + t671 * t636 + t679 * t473 - (t631 * t689 + t636 * t678) * t726 + t796) * MDP(31) + (t443 * t817 + (-t550 * t663 - t557 * t632 - t730 * t817) * pkin(3) + t807) * MDP(24) - t550 * t549 * MDP(11) + (-t549 ^ 2 + t550 ^ 2) * MDP(12) + t563 * MDP(15) + t814; (t433 * t817 + t793) * MDP(23) + (t432 * t817 + t807) * MDP(24) + (-pkin(4) * t421 - t433 * t471 + (-pkin(10) * t437 - t432 * t726 + t813) * t631 + (-(-pkin(10) * qJD(5) - t458) * t726 + t651) * t636 + t797) * MDP(30) + (-pkin(4) * t420 - (t432 * t636 + t458 * t631) * t726 - t433 * t473 + t430 * t812 + (-t726 * t729 - t435) * pkin(10) + t796) * MDP(31) + t814; t473 * t471 * MDP(25) + (-t471 ^ 2 + t473 ^ 2) * MDP(26) + (-t471 * t726 + t712) * MDP(27) + (-t473 * t726 - t698) * MDP(28) + t437 * MDP(29) + (-t416 * t726 - t430 * t473 - g(1) * t489 + g(2) * t806 - g(3) * (-t559 * t631 - t714) + t700) * MDP(30) + (t670 * t726 + t430 * t471 + g(1) * t490 + g(2) * t805 - g(3) * (-t559 * t636 + t715) - t672) * MDP(31) + (-MDP(27) * t777 - MDP(28) * t473 - MDP(30) * t416 + MDP(31) * t670) * qJD(5);];
tau = t1;
