% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:17:00
% EndTime: 2019-03-09 12:17:21
% DurationCPUTime: 13.45s
% Computational Cost: add. (18085->711), mult. (34818->892), div. (0->0), fcn. (34300->6), ass. (0->391)
t518 = cos(qJ(2));
t853 = pkin(7) - pkin(8);
t480 = t853 * t518;
t514 = sin(qJ(4));
t517 = cos(qJ(4));
t515 = sin(qJ(2));
t602 = t853 * t515;
t340 = t517 * t480 + t514 * t602;
t516 = cos(qJ(5));
t712 = t516 * mrSges(6,1);
t513 = sin(qJ(5));
t715 = t513 * mrSges(6,2);
t591 = t712 - t715;
t692 = t340 * t591;
t720 = t340 * mrSges(5,1);
t449 = t514 * t515 + t517 * t518;
t587 = t513 * pkin(5) - qJ(6) * t516;
t133 = -t449 * t587 + t340;
t711 = t516 * mrSges(7,1);
t714 = t513 * mrSges(7,3);
t461 = t711 + t714;
t873 = t133 * t461;
t876 = -t692 - t720 - t873;
t450 = -t518 * t514 + t515 * t517;
t803 = t514 * t480 - t517 * t602;
t132 = t450 * t587 + t803;
t875 = t132 * t133;
t519 = -pkin(2) - pkin(3);
t456 = -t514 * qJ(3) + t517 * t519;
t679 = t513 * qJ(6);
t588 = t516 * pkin(5) + t679;
t458 = -pkin(4) - t588;
t363 = -t456 - t458;
t874 = t133 * t363;
t700 = t133 * t458;
t510 = t513 ^ 2;
t511 = t516 ^ 2;
t657 = t510 + t511;
t866 = t657 * t456;
t872 = t363 + t866;
t454 = pkin(4) - t456;
t871 = t454 + t866;
t814 = -t518 * pkin(2) - t515 * qJ(3);
t459 = -pkin(1) + t814;
t408 = t518 * pkin(3) - t459;
t202 = pkin(4) * t449 - pkin(9) * t450 + t408;
t107 = t202 * t516 - t340 * t513;
t662 = t516 * t340;
t676 = t513 * t202;
t108 = t662 + t676;
t708 = t516 * mrSges(7,3);
t716 = t513 * mrSges(7,1);
t590 = -t708 + t716;
t268 = t590 * t450;
t710 = t516 * mrSges(6,2);
t717 = t513 * mrSges(6,1);
t466 = t710 + t717;
t269 = t466 * t450;
t690 = t449 * t513;
t281 = -mrSges(6,2) * t450 + mrSges(6,3) * t690;
t688 = t449 * t516;
t650 = mrSges(7,2) * t688;
t718 = t450 * mrSges(7,1);
t284 = t650 + t718;
t285 = mrSges(6,1) * t450 + mrSges(6,3) * t688;
t290 = mrSges(7,2) * t690 + mrSges(7,3) * t450;
t691 = t449 * qJ(6);
t80 = t108 + t691;
t81 = -pkin(5) * t449 - t107;
t820 = t590 * t449;
t822 = t466 * t449;
t870 = t107 * t285 + t108 * t281 - t132 * t820 + t133 * t268 + t340 * t269 - t81 * t284 + t80 * t290 - t803 * t822;
t728 = Ifges(7,6) * t513;
t473 = Ifges(7,4) * t516 + t728;
t501 = Ifges(6,6) * t513;
t732 = Ifges(7,5) * t513;
t477 = Ifges(7,1) * t516 + t732;
t172 = -Ifges(7,4) * t450 + t449 * t477;
t737 = Ifges(6,4) * t513;
t479 = Ifges(6,1) * t516 - t737;
t175 = -Ifges(6,5) * t450 + t449 * t479;
t623 = t172 / 0.2e1 + t175 / 0.2e1;
t503 = Ifges(7,5) * t516;
t813 = Ifges(7,3) * t513 + t503;
t164 = -Ifges(7,6) * t450 + t449 * t813;
t736 = Ifges(6,4) * t516;
t475 = -Ifges(6,2) * t513 + t736;
t169 = -Ifges(6,6) * t450 + t449 * t475;
t624 = t169 / 0.2e1 - t164 / 0.2e1;
t733 = Ifges(6,5) * t516;
t862 = Ifges(6,3) + Ifges(7,2);
t869 = t624 * t513 - t623 * t516 + t408 * mrSges(5,1) - Ifges(5,4) * t450 + (-t501 + t733 + t473) * t450 / 0.2e1 + (Ifges(5,2) + t862) * t449 / 0.2e1;
t868 = m(5) * t408;
t614 = t290 / 0.2e1 + t281 / 0.2e1;
t867 = t614 * t516;
t693 = t803 * t514;
t694 = t340 * t517;
t864 = t694 + t693;
t863 = t269 / 0.2e1 + t268 / 0.2e1;
t744 = pkin(4) * t340;
t825 = Ifges(7,4) + Ifges(6,5);
t861 = t340 * t454;
t860 = t340 * t803;
t827 = mrSges(7,2) + mrSges(6,3);
t548 = t827 * t657;
t806 = -t164 / 0.4e1 + t169 / 0.4e1;
t859 = t172 / 0.4e1 + t175 / 0.4e1;
t751 = t516 / 0.2e1;
t754 = -t513 / 0.2e1;
t854 = -Ifges(6,6) / 0.2e1;
t802 = Ifges(7,6) * t751 + t516 * t854 + t754 * t825 + Ifges(5,6);
t847 = t803 * mrSges(5,2);
t857 = -t802 * t450 + t847;
t174 = t449 * Ifges(7,4) + t450 * t477;
t177 = t449 * Ifges(6,5) + t450 * t479;
t649 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t546 = t449 * t649 + t174 / 0.2e1 + t177 / 0.2e1;
t648 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t686 = t450 * t516;
t387 = Ifges(7,5) * t686;
t687 = t450 * t513;
t166 = t449 * Ifges(7,6) + Ifges(7,3) * t687 + t387;
t171 = t449 * Ifges(6,6) + t450 * t475;
t647 = Ifges(7,6) / 0.2e1 + t854;
t804 = t449 * t647 + t166 / 0.2e1 - t171 / 0.2e1;
t830 = -t450 / 0.2e1;
t856 = -t804 * t513 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + t648) * t450 - t408 * mrSges(5,2) + Ifges(5,1) * t830 + Ifges(5,4) * t449 - t546 * t516;
t786 = m(7) / 0.2e1;
t855 = 0.2e1 * t786;
t852 = -t461 / 0.2e1;
t851 = -t513 / 0.4e1;
t848 = m(7) * t587;
t828 = mrSges(6,1) + mrSges(7,1);
t846 = t107 + t81;
t845 = -qJ(6) * t290 / 0.2e1;
t457 = t517 * qJ(3) + t514 * t519;
t843 = t457 * t803;
t842 = t513 * t803;
t841 = t516 * t803;
t609 = -t501 / 0.4e1 + t733 / 0.4e1 + t473 / 0.4e1;
t840 = (t513 * t647 + t516 * t649 + t609) * t449;
t838 = t285 + t284;
t651 = mrSges(6,3) * t687;
t282 = -mrSges(6,2) * t449 - t651;
t719 = t449 * mrSges(7,3);
t289 = -mrSges(7,2) * t687 + t719;
t818 = t289 + t282;
t817 = t290 + t281;
t476 = Ifges(7,1) * t513 - t503;
t478 = Ifges(6,1) * t513 + t736;
t606 = -t478 / 0.2e1 - t476 / 0.2e1;
t468 = -Ifges(7,3) * t516 + t732;
t474 = Ifges(6,2) * t516 + t737;
t607 = t474 / 0.2e1 - t468 / 0.2e1;
t532 = (-t479 / 0.2e1 - t477 / 0.2e1 + t607) * t513 + (-t475 / 0.2e1 + t813 / 0.2e1 + t606) * t516;
t837 = -t587 * t461 - t532;
t455 = -pkin(9) + t457;
t601 = t657 * t455;
t553 = (-t457 + t601) * t517;
t834 = -t268 / 0.2e1;
t749 = t517 / 0.2e1;
t562 = m(7) * t588;
t826 = Ifges(3,4) - Ifges(4,5);
t824 = t591 + mrSges(5,1);
t767 = -t284 / 0.2e1;
t617 = t767 - t285 / 0.2e1;
t819 = t617 * t513;
t446 = t514 * t461;
t447 = t514 * t591;
t816 = -t447 - t446;
t815 = t590 + t848;
t812 = (-mrSges(5,2) + t548) * t456;
t661 = t516 * t461;
t628 = t661 / 0.2e1;
t784 = mrSges(7,1) / 0.2e1;
t811 = (t784 + t628) * t450 + t650;
t810 = (-t363 + t458) * t786 - t461;
t703 = t108 * t513;
t640 = t703 / 0.2e1;
t641 = -t703 / 0.2e1;
t809 = (t640 + t641) * mrSges(6,3);
t608 = -t474 / 0.4e1 + t468 / 0.4e1;
t759 = -t466 / 0.2e1;
t760 = -t590 / 0.2e1;
t564 = -t848 / 0.2e1 + t759 + t760;
t627 = t516 * t749;
t671 = t513 * t517;
t630 = -t671 / 0.2e1;
t644 = -t710 / 0.2e1;
t747 = t517 * t848;
t583 = -t747 / 0.2e1 + mrSges(7,3) * t627 + t517 * t644 + t828 * t630;
t535 = t517 * t564 + t583;
t808 = t535 * qJD(5);
t493 = m(7) * qJ(6) + mrSges(7,3);
t805 = m(7) * t363 + t461;
t605 = t478 / 0.4e1 + t476 / 0.4e1;
t800 = -t513 * (t475 / 0.4e1 - t813 / 0.4e1 + t605) + t516 * (t479 / 0.4e1 + t477 / 0.4e1 + t608);
t599 = mrSges(7,2) * pkin(5) - t825;
t655 = qJD(5) * t516;
t799 = -qJD(5) * (mrSges(7,2) * t679 + t501 - t728) - t599 * t655;
t797 = -t609 * t449 + t825 * t688 / 0.2e1 + t862 * t830;
t261 = t588 * t450;
t270 = t450 * t468;
t271 = t450 * t474;
t272 = -Ifges(7,1) * t687 + t387;
t273 = t450 * t478;
t796 = -t132 * t760 + t261 * t852 - t803 * t759 - t587 * t834 - t516 * t271 / 0.4e1 + (t273 + t171) * t851 + (t272 + t166) * t513 / 0.4e1 + (t270 + t177 + t174) * t516 / 0.4e1;
t795 = 0.2e1 * pkin(9);
t794 = 2 * qJD(4);
t793 = m(5) / 0.2e1;
t792 = -m(6) / 0.2e1;
t791 = -m(6) / 0.4e1;
t790 = m(6) / 0.2e1;
t789 = m(6) / 0.4e1;
t788 = -m(7) / 0.2e1;
t787 = -m(7) / 0.4e1;
t785 = m(7) / 0.4e1;
t783 = -mrSges(7,3) / 0.2e1;
t742 = pkin(9) * t449;
t300 = pkin(4) * t450 + t742;
t496 = t518 * qJ(3);
t445 = t515 * t519 + t496;
t204 = -t300 + t445;
t109 = t204 * t516 - t842;
t743 = pkin(5) * t450;
t83 = -t109 + t743;
t780 = t83 / 0.2e1;
t779 = m(7) * t83;
t128 = t300 * t516 + t842;
t106 = -t128 - t743;
t777 = t106 / 0.2e1;
t776 = -t109 / 0.2e1;
t110 = t513 * t204 + t841;
t775 = t110 / 0.2e1;
t774 = -t128 / 0.2e1;
t129 = t513 * t300 - t841;
t773 = t129 / 0.2e1;
t263 = t591 * t450;
t772 = -t263 / 0.2e1;
t771 = -t822 / 0.2e1;
t770 = t822 / 0.2e1;
t769 = t284 / 0.2e1;
t766 = -t363 / 0.2e1;
t765 = t454 / 0.2e1;
t764 = t456 / 0.2e1;
t762 = t458 / 0.2e1;
t761 = -t591 / 0.2e1;
t748 = m(7) * t106;
t745 = m(7) * t516;
t706 = qJ(6) * t450;
t82 = t110 - t706;
t740 = t82 * mrSges(7,2);
t739 = t83 * mrSges(7,2);
t738 = m(7) * qJD(6);
t725 = t109 * mrSges(6,3);
t724 = t110 * mrSges(6,3);
t287 = mrSges(6,1) * t449 - mrSges(6,3) * t686;
t288 = -mrSges(7,1) * t449 + mrSges(7,2) * t686;
t463 = -t518 * mrSges(4,1) - t515 * mrSges(4,3);
t603 = m(4) * t459 + t463;
t3 = (-mrSges(3,2) * pkin(1) - mrSges(4,3) * t459 + t826 * t518) * t518 + (-pkin(1) * mrSges(3,1) + t459 * mrSges(4,1) + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t518 - t826 * t515) * t515 + t603 * (pkin(2) * t515 - t496) + m(7) * (t80 * t82 + t81 * t83 - t875) + (t445 * mrSges(5,2) - t869) * t450 + t445 * t868 + m(6) * (t107 * t109 + t108 * t110 - t860) + t110 * t282 + t109 * t287 + t83 * t288 + t82 * t289 + (t445 * mrSges(5,1) - t856) * t449 - t870;
t723 = t3 * qJD(1);
t713 = t513 * t80;
t709 = t516 * mrSges(7,2);
t262 = t461 * t450;
t9 = -t108 * t288 - m(7) * (t108 * t81 + t132 * t261) - t261 * t268 - t132 * t262 - t803 * t263 + t108 * t287 + ((t270 / 0.2e1 - t271 / 0.2e1 + t81 * mrSges(7,2) + t546) * t513 + (-t272 / 0.2e1 + t273 / 0.2e1 + t80 * mrSges(7,2) + t108 * mrSges(6,3) - t804) * t516) * t450 + (-m(7) * t80 - t651 - t818) * t107;
t707 = t9 * qJD(1);
t575 = t107 * t516 + t703;
t585 = -t516 * t81 + t713;
t663 = t516 * t288;
t664 = t516 * t287;
t672 = t513 * t515;
t16 = -t818 * t672 + (-m(6) * t575 - m(7) * t585 - mrSges(5,1) * t449 - mrSges(5,2) * t450 + t603 + t663 - t664 - t868) * t515;
t705 = qJD(1) * t16;
t35 = m(7) * (-t132 * t686 + t449 * t80) - t268 * t686 + t449 * t289;
t704 = qJD(1) * t35;
t702 = t128 * t513;
t701 = t129 * t516;
t698 = t132 * t587;
t697 = t132 * t513;
t689 = t449 * t514;
t685 = t450 * t517;
t681 = t456 * t513;
t675 = t513 * t268;
t669 = t515 * t516;
t660 = t866 * pkin(9);
t659 = t657 * pkin(9) * t517;
t658 = -t514 * mrSges(5,1) - t517 * mrSges(5,2);
t160 = (-t515 / 0.2e1 - t689 / 0.2e1 - t685 / 0.2e1) * t745;
t656 = qJD(1) * t160;
t654 = m(7) * t671;
t653 = t779 / 0.2e1;
t652 = t513 * t738;
t645 = mrSges(7,2) * t754;
t643 = -t709 / 0.2e1;
t639 = -t690 / 0.2e1;
t638 = t690 / 0.2e1;
t633 = -t681 / 0.2e1;
t632 = t456 * t751;
t631 = t513 * t288 / 0.2e1;
t620 = -t820 / 0.2e1 + t771;
t619 = -t269 / 0.2e1 + t834;
t618 = t285 / 0.2e1 + t769;
t616 = -t288 / 0.2e1 + t287 / 0.2e1;
t615 = -t289 / 0.2e1 - t282 / 0.2e1;
t611 = t446 / 0.2e1 + t447 / 0.2e1;
t610 = t762 + t766;
t592 = t517 * t601;
t586 = t513 * t81 + t516 * t80;
t584 = t83 * t513 + t82 * t516;
t105 = t129 + t706;
t544 = -pkin(9) * t614 + t806;
t545 = t618 * pkin(9) + t859;
t558 = t454 * t770 - t766 * t820;
t565 = pkin(4) * t771 + t762 * t820;
t570 = t843 + t861;
t571 = t132 * t457 + t874;
t1 = t619 * t457 - t700 * t786 + t571 * t788 + (t615 * t456 - t614 * t455 + (t773 + t775) * mrSges(6,3) + (t105 / 0.2e1 + t82 / 0.2e1) * mrSges(7,2) + (t105 * t455 + t456 * t80) * t788 + (t110 * t789 + t785 * t82) * t795 + t544 - t806) * t516 + (t616 * t456 - t617 * t455 + (t774 + t776) * mrSges(6,3) + (t777 + t780) * mrSges(7,2) + (t106 * t455 + t456 * t81) * t788 + (t779 / 0.4e1 + t109 * t791) * t795 + t545 - t859) * t513 + t558 + t565 + (t570 - t744 + (t108 * t456 + t129 * t455) * t516 + (-t107 * t456 - t128 * t455) * t513) * t792;
t29 = (-m(7) - m(6)) * t455 * t866 + (-m(6) * t454 - t805 - t824) * t457 + t812;
t582 = -t1 * qJD(1) - t29 * qJD(2);
t40 = t363 * t815 + t454 * t466 - t837;
t529 = t663 / 0.2e1 + (t575 - t585) * t786 - t664 / 0.2e1 + t818 * t754;
t520 = (t261 * t363 - t698) * t786 + t263 * t765 + t529 * t455 - t796 + t363 * t262 / 0.2e1 + t809 + t846 * t643 + (t713 / 0.2e1 + t641) * mrSges(7,2);
t541 = t827 * (-t510 / 0.2e1 - t511 / 0.2e1);
t526 = t455 * t541 - t800;
t530 = (-pkin(5) * t83 + qJ(6) * t82) * t788 + pkin(5) * t769 - t845 + mrSges(6,1) * t776 + mrSges(6,2) * t775 + t82 * t783 + mrSges(7,1) * t780;
t5 = t520 - t840 + (t526 + t648) * t450 + t530;
t581 = t5 * qJD(1) - t40 * qJD(2);
t563 = t287 * t630 + t514 * t863 + t517 * t631 + t818 * t627;
t572 = t701 - t702;
t574 = -t107 * t513 + t108 * t516;
t576 = t105 * t516 + t106 * t513;
t13 = (t820 / 0.2e1 + t770 + (-t133 + t586) * t786 + (-t340 + t574) * t790) * t517 + (t867 + t819 + (t132 + t576) * t786 + (t803 + t572) * t790) * t514 + t563;
t4 = t105 * t289 + t106 * t288 + t128 * t287 + t129 * t282 + m(6) * (t107 * t128 + t108 * t129 + t860) + m(7) * (t105 * t80 + t106 * t81 + t875) + t869 * t450 + t856 * t449 + t870;
t580 = t4 * qJD(1) + t13 * qJD(3);
t524 = -m(5) * t864 / 0.2e1 + (t517 * t574 + t693) * t792 + (t132 * t514 + t517 * t586) * t788;
t573 = -t109 * t513 + t110 * t516;
t525 = t864 * t793 + (t514 * t573 + t694) * t790 + (t133 * t517 + t514 * t584) * t786;
t540 = t513 * t618 - t867;
t11 = (t513 * t616 + t516 * t615 + t620) * t517 + (t540 + t619) * t514 + t524 + t525;
t537 = t517 * t548 + t658;
t39 = -mrSges(4,3) - m(4) * qJ(3) - m(7) * (t363 * t514 + t592) - m(6) * (t454 * t514 + t592) - m(5) * (-t456 * t514 + t457 * t517) + t537 + t816;
t578 = -t11 * qJD(1) - t39 * qJD(2);
t522 = (t261 * t788 - t262 / 0.2e1 + t772) * t517 + (t450 * t541 + t529) * t514;
t15 = (-t562 / 0.2e1 - t714 / 0.2e1 - t712 / 0.2e1 - t711 / 0.2e1 + t715 / 0.2e1) * t515 + t522;
t577 = t15 * qJD(1) - qJD(2) * t535;
t238 = t805 * t513;
t533 = (t697 + (-t363 * t450 + t449 * t455) * t516) * t786 + t675 / 0.2e1;
t30 = t653 - t533 + t811;
t569 = -qJD(1) * t30 + qJD(2) * t238;
t46 = t719 + 0.2e1 * (t691 / 0.2e1 + t676 / 0.4e1 + t662 / 0.4e1 - t108 / 0.4e1) * m(7);
t568 = qJD(1) * t46 + qJD(5) * t493;
t567 = (t458 * t514 + t659) * t786 + (-pkin(4) * t514 + t659) * t790;
t349 = (-0.1e1 + t657) * t517 * t514;
t137 = 0.4e1 * (t785 + t789) * t349;
t539 = (t852 + t761) * t514 + t567;
t26 = 0.2e1 * (t787 * t872 + t791 * t871) * t514 + t537 + t539 - t611 + 0.2e1 * (t787 + t791) * t553;
t550 = t13 * qJD(1) - t26 * qJD(2) + t137 * qJD(3);
t308 = (m(7) * t458 - t461) * t513;
t534 = (-t697 + (-t450 * t458 + t742) * t516) * t786 - t675 / 0.2e1;
t33 = -t748 / 0.2e1 + t534 + t811;
t84 = (t461 + (t764 - t610) * m(7)) * t513;
t549 = -qJD(1) * t33 + qJD(2) * t84 + qJD(4) * t308;
t523 = (t765 + pkin(4) / 0.2e1) * t466 - t610 * t590 - t810 * t587 + t532;
t19 = (t848 / 0.2e1 + t710 / 0.2e1 - t708 / 0.2e1 + t717 / 0.2e1 + t716 / 0.2e1) * t456 + t523;
t49 = -pkin(4) * t466 + t458 * t815 + t837;
t521 = (t261 * t458 + t698) * t786 + t262 * t762 + pkin(4) * t772 + t529 * pkin(9) + t796 + mrSges(7,2) * t640 + t80 * t645 + t809 + t846 * t709 / 0.2e1;
t527 = pkin(9) * t541 + t800;
t528 = (-pkin(5) * t106 + qJ(6) * t105) * t788 + pkin(5) * t767 + t845 + t105 * t783 + mrSges(7,1) * t777 + mrSges(6,1) * t774 + mrSges(6,2) * t773;
t7 = t521 + (t527 - t648) * t450 + t840 + t528;
t76 = t747 / 0.2e1 + ((t783 + mrSges(6,2) / 0.2e1) * t516 + (t784 + mrSges(6,1) / 0.2e1) * t513 + t564) * t517;
t542 = t7 * qJD(1) + t19 * qJD(2) + t76 * qJD(3) + t49 * qJD(4);
t538 = t513 * t607 + t516 * t606 - Ifges(5,5);
t531 = (-t461 - t591 - t562) * qJD(5);
t460 = (m(7) * pkin(9) + mrSges(7,2)) * t516;
t354 = (m(7) * t455 - mrSges(7,2)) * t516;
t161 = (t685 + t689) * t745 / 0.2e1 + t669 * t788;
t85 = t513 * t810 + t681 * t786;
t78 = t583 + (t466 + t815) * t749;
t41 = t80 * t855 + t289;
t34 = t661 * t830 + t718 / 0.2e1 + t653 + t533;
t32 = t450 * t628 - t718 / 0.2e1 + t748 / 0.2e1 + t534;
t28 = (t514 * t872 + t553) * t786 + (t514 * t871 + t553) * t790 + t539 + t611;
t20 = mrSges(7,3) * t632 + t456 * t644 + t633 * t828 - t764 * t848 + t523;
t14 = t515 * t562 / 0.2e1 + t522 + (mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1) * t672 + t828 * t669 / 0.2e1;
t12 = t13 * qJD(4);
t10 = (m(4) * pkin(7) + mrSges(4,2)) * t518 + (t450 * mrSges(5,3) + t540) * t514 - t524 + t525 + t563 + (-t449 * mrSges(5,3) + t620) * t517;
t8 = Ifges(6,6) * t638 + Ifges(7,6) * t639 + t527 * t450 + t521 - t528 - t797;
t6 = Ifges(6,6) * t639 + Ifges(7,6) * t638 + t526 * t450 + t520 - t530 + t797;
t2 = (t478 + t476) * t688 / 0.4e1 + (t513 * t608 + Ifges(5,5)) * t449 + t608 * t690 + (t605 * t449 + t724 / 0.2e1 + t740 / 0.2e1 + t544 + t806) * t516 + t565 - t558 + (-t725 / 0.2e1 + t739 / 0.2e1 + t545) * t513 - t857 + t692 / 0.2e1 + t457 * t863 + t873 + (-t701 / 0.2e1 + t702 / 0.2e1) * mrSges(6,3) + (pkin(9) * t573 + t456 * t574 + t570 + t744) * t790 + (-t175 - t172) * t851 - (t761 - mrSges(5,1) / 0.2e1) * t340 + (pkin(9) * t584 + t456 * t586 + t571 - t700) * t786 + (t572 * t790 + t576 * t786 + t751 * t817 + t819) * t455 + t720 / 0.2e1 + t818 * t632 + t456 * t631 + t287 * t633 + t105 * t643 + t106 * t645;
t17 = [qJD(2) * t3 - qJD(3) * t16 + qJD(4) * t4 - qJD(5) * t9 + qJD(6) * t35, t10 * qJD(3) + t2 * qJD(4) + t6 * qJD(5) + t34 * qJD(6) + t723 + (t363 * t820 + t454 * t822 + t847 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t518 + (-mrSges(4,2) * qJ(3) - Ifges(3,6) + Ifges(4,6)) * t515 + (-t455 * t817 - t624 - t724 - t740) * t516 + (t455 * t838 - t623 + t725 - t739) * t513 + 0.2e1 * (t455 * t573 - t861) * t790 + (t455 * t584 - t874) * t855 + 0.2e1 * (t340 * t456 + t843) * t793 + (t457 * mrSges(5,3) - t802) * t450 + (-t456 * mrSges(5,3) + t538) * t449 + (m(4) * t814 - t518 * mrSges(3,1) + t515 * mrSges(3,2) + t463) * pkin(7) + t876) * qJD(2), qJD(2) * t10 + qJD(5) * t14 + qJD(6) * t161 + t12 - t705, t2 * qJD(2) + t8 * qJD(5) + t32 * qJD(6) + ((pkin(9) * t572 - t744) * t790 + (pkin(9) * t576 + t700) * t786) * t794 + t580 + (pkin(4) * t822 - t458 * t820 + (t105 * mrSges(7,2) + t129 * mrSges(6,3) + pkin(9) * t817 - t624) * t516 + (t106 * mrSges(7,2) - t128 * mrSges(6,3) - pkin(9) * t838 - t623) * t513 + t538 * t449 + t857 + t876) * qJD(4), t6 * qJD(2) + t14 * qJD(3) + t8 * qJD(4) + t41 * qJD(6) - t707 + ((-m(7) * pkin(5) - t828) * t108 + (-mrSges(6,2) + t493) * t107 + ((-mrSges(7,2) * qJ(6) - Ifges(6,6) + Ifges(7,6)) * t516 + t599 * t513) * t450) * qJD(5), qJD(2) * t34 + qJD(3) * t161 + qJD(4) * t32 + qJD(5) * t41 + t704; -qJD(3) * t11 - qJD(4) * t1 + qJD(5) * t5 - qJD(6) * t30 - t723, -qJD(3) * t39 - qJD(4) * t29 - qJD(5) * t40 + qJD(6) * t238, 0.2e1 * (t786 + t790) * t349 * qJD(3) + t28 * qJD(4) + t78 * qJD(5) + t578, t28 * qJD(3) + t20 * qJD(5) + t85 * qJD(6) + ((t457 * t458 + t660) * t786 + (-pkin(4) * t457 + t660) * t790) * t794 + t582 + ((-t461 - t824) * t457 + t812) * qJD(4), t78 * qJD(3) + t20 * qJD(4) + t354 * qJD(6) + t455 * t531 + t581 - t799, qJD(4) * t85 + qJD(5) * t354 + t569; qJD(2) * t11 + qJD(5) * t15 - qJD(6) * t160 + t12 + t705, -qJD(4) * t26 - t517 * t652 - t578 - t808, t137 * qJD(4) (t658 + t816) * qJD(4) + t808 + t567 * t794 + (qJD(4) * t548 + t652) * t517 + t550, t535 * qJD(4) + (t516 * t738 + t531) * t514 + t577, -t656 + (t514 * t655 + (-qJD(2) + qJD(4)) * t671) * m(7); qJD(2) * t1 + qJD(5) * t7 + qJD(6) * t33 - t580, qJD(3) * t26 + qJD(5) * t19 - qJD(6) * t84 - t582, qJD(5) * t76 - t550, qJD(5) * t49 - qJD(6) * t308, pkin(9) * t531 + t460 * qJD(6) + t542 + t799, qJD(5) * t460 - t549; -qJD(2) * t5 - qJD(3) * t15 - qJD(4) * t7 + qJD(6) * t46 + t707, qJD(3) * t535 - qJD(4) * t19 - t581, -qJD(4) * t76 - t577, -t542, t493 * qJD(6), t568; qJD(2) * t30 + qJD(3) * t160 - qJD(4) * t33 - qJD(5) * t46 - t704, qJD(3) * t654 + qJD(4) * t84 - t569, qJD(2) * t654 + t656, t549, -t568, 0;];
Cq  = t17;
