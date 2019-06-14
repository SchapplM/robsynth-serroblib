% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:07:13
% EndTime: 2019-05-06 19:07:35
% DurationCPUTime: 9.17s
% Computational Cost: add. (130494->361), mult. (290960->436), div. (0->0), fcn. (209824->10), ass. (0->152)
t841 = -2 * qJD(3);
t840 = Ifges(3,1) + Ifges(4,2);
t839 = Ifges(6,1) + Ifges(7,1);
t828 = Ifges(3,4) + Ifges(4,6);
t827 = Ifges(6,4) + Ifges(7,4);
t826 = Ifges(3,5) - Ifges(4,4);
t825 = Ifges(6,5) + Ifges(7,5);
t838 = Ifges(3,2) + Ifges(4,3);
t837 = Ifges(6,2) + Ifges(7,2);
t824 = Ifges(3,6) - Ifges(4,5);
t836 = Ifges(6,6) + Ifges(7,6);
t835 = Ifges(3,3) + Ifges(4,1);
t834 = Ifges(6,3) + Ifges(7,3);
t777 = cos(pkin(6));
t771 = t777 * qJD(1) + qJD(2);
t780 = sin(qJ(2));
t776 = sin(pkin(6));
t808 = qJD(1) * t776;
t799 = t780 * t808;
t833 = (pkin(2) * t771 + t841) * t799;
t781 = sin(qJ(1));
t785 = cos(qJ(1));
t766 = t781 * g(1) - t785 * g(2);
t786 = qJD(1) ^ 2;
t751 = t786 * t776 * pkin(8) + qJDD(1) * pkin(1) + t766;
t767 = -t785 * g(1) - t781 * g(2);
t805 = qJDD(1) * t776;
t752 = -t786 * pkin(1) + pkin(8) * t805 + t767;
t784 = cos(qJ(2));
t819 = t777 * t780;
t821 = t776 * t780;
t717 = -g(3) * t821 + t751 * t819 + t784 * t752;
t753 = (-pkin(2) * t784 - qJ(3) * t780) * t808;
t769 = t771 ^ 2;
t770 = t777 * qJDD(1) + qJDD(2);
t807 = qJD(1) * t784;
t800 = t776 * t807;
t687 = t769 * pkin(2) - t770 * qJ(3) - t753 * t800 + t771 * t841 - t717;
t832 = -pkin(2) - pkin(9);
t831 = t777 * g(3);
t830 = mrSges(3,1) - mrSges(4,2);
t829 = -mrSges(6,2) - mrSges(7,2);
t822 = t776 ^ 2 * t786;
t820 = t776 * t784;
t818 = t777 * t784;
t731 = -t776 * t751 - t831;
t747 = t771 * mrSges(3,1) - mrSges(3,3) * t799;
t748 = -t771 * mrSges(3,2) + mrSges(3,3) * t800;
t750 = mrSges(4,1) * t799 + t771 * mrSges(4,2);
t757 = (qJD(2) * t807 + qJDD(1) * t780) * t776;
t758 = -qJD(2) * t799 + t784 * t805;
t688 = -t758 * pkin(2) + (-t771 * t800 - t757) * qJ(3) + t731 + t833;
t749 = -mrSges(4,1) * t800 - t771 * mrSges(4,3);
t756 = pkin(3) * t799 - t771 * pkin(9);
t804 = t784 ^ 2 * t822;
t677 = -pkin(3) * t804 - t831 - t757 * qJ(3) + t832 * t758 + (-t751 + (-qJ(3) * t771 * t784 - t756 * t780) * qJD(1)) * t776 + t833;
t809 = g(3) * t820 + t780 * t752;
t793 = -t769 * qJ(3) + t753 * t799 + qJDD(3) + t809;
t679 = t757 * pkin(3) + t832 * t770 + (-pkin(3) * t771 * t808 - pkin(9) * t780 * t822 - t751 * t777) * t784 + t793;
t779 = sin(qJ(4));
t783 = cos(qJ(4));
t672 = t783 * t677 + t779 * t679;
t741 = t783 * t771 - t779 * t800;
t714 = -t741 * qJD(4) - t783 * t758 - t779 * t770;
t740 = -t779 * t771 - t783 * t800;
t718 = -t740 * mrSges(5,1) + t741 * mrSges(5,2);
t762 = qJD(4) + t799;
t724 = t762 * mrSges(5,1) - t741 * mrSges(5,3);
t746 = qJDD(4) + t757;
t719 = -t740 * pkin(4) - t741 * pkin(10);
t760 = t762 ^ 2;
t667 = -t760 * pkin(4) + t746 * pkin(10) + t740 * t719 + t672;
t676 = t758 * pkin(3) - pkin(9) * t804 + t771 * t756 - t687;
t715 = t740 * qJD(4) - t779 * t758 + t783 * t770;
t670 = (-t740 * t762 - t715) * pkin(10) + (t741 * t762 - t714) * pkin(4) + t676;
t778 = sin(qJ(5));
t782 = cos(qJ(5));
t662 = -t778 * t667 + t782 * t670;
t721 = -t778 * t741 + t782 * t762;
t684 = t721 * qJD(5) + t782 * t715 + t778 * t746;
t722 = t782 * t741 + t778 * t762;
t699 = -t721 * mrSges(7,1) + t722 * mrSges(7,2);
t700 = -t721 * mrSges(6,1) + t722 * mrSges(6,2);
t739 = qJD(5) - t740;
t703 = -t739 * mrSges(6,2) + t721 * mrSges(6,3);
t712 = qJDD(5) - t714;
t659 = -0.2e1 * qJD(6) * t722 + (t721 * t739 - t684) * qJ(6) + (t721 * t722 + t712) * pkin(5) + t662;
t702 = -t739 * mrSges(7,2) + t721 * mrSges(7,3);
t802 = m(7) * t659 + t712 * mrSges(7,1) + t739 * t702;
t652 = m(6) * t662 + t712 * mrSges(6,1) + t739 * t703 + (-t699 - t700) * t722 + (-mrSges(6,3) - mrSges(7,3)) * t684 + t802;
t663 = t782 * t667 + t778 * t670;
t683 = -t722 * qJD(5) - t778 * t715 + t782 * t746;
t704 = t739 * pkin(5) - t722 * qJ(6);
t720 = t721 ^ 2;
t661 = -t720 * pkin(5) + t683 * qJ(6) + 0.2e1 * qJD(6) * t721 - t739 * t704 + t663;
t801 = m(7) * t661 + t683 * mrSges(7,3) + t721 * t699;
t705 = t739 * mrSges(7,1) - t722 * mrSges(7,3);
t813 = -t739 * mrSges(6,1) + t722 * mrSges(6,3) - t705;
t655 = m(6) * t663 + t683 * mrSges(6,3) + t721 * t700 + t829 * t712 + t813 * t739 + t801;
t796 = -t778 * t652 + t782 * t655;
t648 = m(5) * t672 - t746 * mrSges(5,2) + t714 * mrSges(5,3) + t740 * t718 - t762 * t724 + t796;
t671 = -t779 * t677 + t783 * t679;
t723 = -t762 * mrSges(5,2) + t740 * mrSges(5,3);
t666 = -t746 * pkin(4) - t760 * pkin(10) + t741 * t719 - t671;
t664 = -t683 * pkin(5) - t720 * qJ(6) + t722 * t704 + qJDD(6) + t666;
t795 = m(7) * t664 - t683 * mrSges(7,1) - t721 * t702;
t787 = -m(6) * t666 + t683 * mrSges(6,1) + t829 * t684 + t721 * t703 + t813 * t722 - t795;
t656 = m(5) * t671 + t746 * mrSges(5,1) - t715 * mrSges(5,3) - t741 * t718 + t762 * t723 + t787;
t797 = t783 * t648 - t779 * t656;
t794 = m(4) * t688 - t757 * mrSges(4,3) + t749 * t800 + t797;
t638 = m(3) * t731 + t757 * mrSges(3,2) - t830 * t758 + (-t748 * t784 + (t747 - t750) * t780) * t808 + t794;
t803 = t751 * t818;
t716 = t803 - t809;
t754 = (mrSges(4,2) * t784 - mrSges(4,3) * t780) * t808;
t755 = (-mrSges(3,1) * t784 + mrSges(3,2) * t780) * t808;
t641 = t779 * t648 + t783 * t656;
t697 = -t770 * pkin(2) + t793 - t803;
t791 = -m(4) * t697 - t757 * mrSges(4,1) - t641;
t639 = m(3) * t716 - t757 * mrSges(3,3) + (t748 - t749) * t771 + t830 * t770 + (-t754 - t755) * t799 + t791;
t650 = t782 * t652 + t778 * t655;
t789 = -m(5) * t676 + t714 * mrSges(5,1) - t715 * mrSges(5,2) + t740 * t723 - t741 * t724 - t650;
t788 = -m(4) * t687 + t770 * mrSges(4,3) + t771 * t750 + t754 * t800 - t789;
t646 = (mrSges(3,3) + mrSges(4,1)) * t758 + t788 + t755 * t800 + m(3) * t717 - t770 * mrSges(3,2) - t771 * t747;
t628 = -t776 * t638 + t639 * t818 + t646 * t819;
t626 = m(2) * t766 + qJDD(1) * mrSges(2,1) - t786 * mrSges(2,2) + t628;
t633 = -t780 * t639 + t784 * t646;
t632 = m(2) * t767 - t786 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t633;
t817 = t785 * t626 + t781 * t632;
t816 = t836 * t721 + t825 * t722 + t834 * t739;
t815 = -t837 * t721 - t827 * t722 - t836 * t739;
t814 = t827 * t721 + t839 * t722 + t825 * t739;
t812 = (t826 * t780 + t824 * t784) * t808 + t835 * t771;
t811 = (-t828 * t780 - t838 * t784) * t808 - t824 * t771;
t810 = (t840 * t780 + t828 * t784) * t808 + t826 * t771;
t627 = t777 * t638 + t639 * t820 + t646 * t821;
t798 = -t781 * t626 + t785 * t632;
t642 = -mrSges(6,1) * t666 + mrSges(6,3) * t663 - mrSges(7,1) * t664 + mrSges(7,3) * t661 - pkin(5) * t795 + qJ(6) * t801 + (-qJ(6) * t705 + t814) * t739 + (-pkin(5) * t705 - t816) * t722 + (-qJ(6) * mrSges(7,2) + t836) * t712 + (-pkin(5) * mrSges(7,2) + t827) * t684 + t837 * t683;
t657 = -t684 * mrSges(7,3) - t722 * t699 + t802;
t649 = mrSges(6,2) * t666 + mrSges(7,2) * t664 - mrSges(6,3) * t662 - mrSges(7,3) * t659 - qJ(6) * t657 + t827 * t683 + t839 * t684 + t825 * t712 + t816 * t721 + t815 * t739;
t708 = Ifges(5,5) * t741 + Ifges(5,6) * t740 + Ifges(5,3) * t762;
t709 = Ifges(5,4) * t741 + Ifges(5,2) * t740 + Ifges(5,6) * t762;
t629 = mrSges(5,2) * t676 - mrSges(5,3) * t671 + Ifges(5,1) * t715 + Ifges(5,4) * t714 + Ifges(5,5) * t746 - pkin(10) * t650 - t778 * t642 + t782 * t649 + t740 * t708 - t762 * t709;
t710 = Ifges(5,1) * t741 + Ifges(5,4) * t740 + Ifges(5,5) * t762;
t634 = -mrSges(5,1) * t676 - mrSges(6,1) * t662 - mrSges(7,1) * t659 + mrSges(6,2) * t663 + mrSges(7,2) * t661 + mrSges(5,3) * t672 + Ifges(5,4) * t715 + Ifges(5,2) * t714 + Ifges(5,6) * t746 - pkin(4) * t650 - pkin(5) * t657 - t741 * t708 + t762 * t710 + t815 * t722 + t814 * t721 - t834 * t712 - t825 * t684 - t836 * t683;
t640 = t758 * mrSges(4,2) - t750 * t799 + t794;
t623 = -mrSges(3,1) * t731 - mrSges(4,1) * t687 + mrSges(4,2) * t688 + mrSges(3,3) * t717 - pkin(2) * t640 - pkin(3) * t789 - pkin(9) * t797 - t779 * t629 - t783 * t634 + t828 * t757 + t838 * t758 + t824 * t770 + t810 * t771 - t812 * t799;
t624 = pkin(3) * t641 - qJ(3) * t640 + mrSges(5,1) * t671 - mrSges(5,2) * t672 - mrSges(4,3) * t688 + mrSges(4,1) * t697 + Ifges(5,6) * t714 + Ifges(5,5) * t715 - mrSges(3,3) * t716 + pkin(4) * t787 + mrSges(3,2) * t731 - t740 * t710 + t741 * t709 + Ifges(5,3) * t746 + t778 * t649 + pkin(10) * t796 + t782 * t642 + t811 * t771 + t826 * t770 + t828 * t758 + t840 * t757 + t812 * t800;
t792 = pkin(8) * t633 + t623 * t784 + t624 * t780;
t622 = mrSges(3,1) * t716 - mrSges(3,2) * t717 + mrSges(4,2) * t697 - mrSges(4,3) * t687 + t783 * t629 - t779 * t634 - pkin(9) * t641 + pkin(2) * (-t771 * t749 + t791) + qJ(3) * t788 + (-pkin(2) * mrSges(4,2) + t835) * t770 + (qJ(3) * mrSges(4,1) + t824) * t758 + t826 * t757 + (-t810 * t784 + (-pkin(2) * t754 - t811) * t780) * t808;
t621 = -mrSges(2,2) * g(3) - mrSges(2,3) * t766 + Ifges(2,5) * qJDD(1) - t786 * Ifges(2,6) - t780 * t623 + t784 * t624 + (-t627 * t776 - t628 * t777) * pkin(8);
t620 = mrSges(2,1) * g(3) + mrSges(2,3) * t767 + t786 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t627 - t776 * t622 + t792 * t777;
t1 = [-m(1) * g(1) + t798; -m(1) * g(2) + t817; (-m(1) - m(2)) * g(3) + t627; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t817 - t781 * t620 + t785 * t621; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t798 + t785 * t620 + t781 * t621; -mrSges(1,1) * g(2) + mrSges(2,1) * t766 + mrSges(1,2) * g(1) - mrSges(2,2) * t767 + Ifges(2,3) * qJDD(1) + pkin(1) * t628 + t777 * t622 + t792 * t776;];
tauB  = t1;
