% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 03:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR15_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 03:22:57
% EndTime: 2019-05-08 03:24:05
% DurationCPUTime: 50.86s
% Computational Cost: add. (838460->399), mult. (2072939->510), div. (0->0), fcn. (1735395->14), ass. (0->175)
t862 = Ifges(5,1) + Ifges(6,2);
t853 = Ifges(5,4) + Ifges(6,6);
t852 = Ifges(5,5) - Ifges(6,4);
t861 = Ifges(5,2) + Ifges(6,3);
t851 = Ifges(5,6) - Ifges(6,5);
t860 = -Ifges(5,3) - Ifges(6,1);
t804 = cos(pkin(6));
t799 = t804 * qJD(1) + qJD(2);
t801 = sin(pkin(7));
t803 = cos(pkin(7));
t802 = sin(pkin(6));
t812 = cos(qJ(2));
t834 = qJD(1) * t812;
t831 = t802 * t834;
t827 = t803 * t831;
t782 = (t799 * t801 + t827) * pkin(10);
t808 = sin(qJ(2));
t836 = qJD(1) * t802;
t856 = pkin(10) * t801;
t786 = (-pkin(2) * t812 - t808 * t856) * t836;
t833 = qJD(1) * qJD(2);
t792 = (qJDD(1) * t808 + t812 * t833) * t802;
t798 = t804 * qJDD(1) + qJDD(2);
t809 = sin(qJ(1));
t813 = cos(qJ(1));
t796 = t809 * g(1) - t813 * g(2);
t814 = qJD(1) ^ 2;
t857 = pkin(9) * t802;
t789 = qJDD(1) * pkin(1) + t814 * t857 + t796;
t797 = -t813 * g(1) - t809 * g(2);
t790 = -t814 * pkin(1) + qJDD(1) * t857 + t797;
t843 = t804 * t812;
t828 = t789 * t843 - t808 * t790;
t835 = qJD(1) * t808;
t855 = pkin(10) * t803;
t733 = -t792 * t855 + t798 * pkin(2) + t799 * t782 + (-g(3) * t812 - t786 * t835) * t802 + t828;
t832 = t802 * t835;
t785 = t799 * pkin(2) - t832 * t855;
t793 = (qJDD(1) * t812 - t808 * t833) * t802;
t825 = t793 * t803 + t798 * t801;
t844 = t804 * t808;
t837 = t789 * t844 + t812 * t790;
t734 = -t799 * t785 + (-g(3) * t808 + t786 * t834) * t802 + t825 * pkin(10) + t837;
t854 = t804 * g(3);
t740 = -t792 * t856 - t793 * pkin(2) - t854 + (-t789 + (-t782 * t812 + t785 * t808) * qJD(1)) * t802;
t807 = sin(qJ(3));
t811 = cos(qJ(3));
t703 = -t807 * t734 + (t733 * t803 + t740 * t801) * t811;
t845 = t803 * t807;
t849 = t801 * t807;
t772 = t799 * t849 + (t808 * t811 + t812 * t845) * t836;
t755 = -t772 * qJD(3) - t807 * t792 + t825 * t811;
t859 = -2 * qJD(5);
t858 = cos(qJ(4));
t783 = t803 * t799 - t801 * t831 + qJD(3);
t806 = sin(qJ(4));
t762 = t806 * t772 - t858 * t783;
t848 = t801 * t811;
t771 = t799 * t848 - t807 * t832 + t811 * t827;
t769 = qJD(4) - t771;
t850 = t762 * t769;
t847 = t802 * t808;
t846 = t802 * t812;
t704 = t733 * t845 + t811 * t734 + t740 * t849;
t757 = -t771 * mrSges(4,1) + t772 * mrSges(4,2);
t765 = t783 * mrSges(4,1) - t772 * mrSges(4,3);
t773 = -t801 * t793 + t803 * t798 + qJDD(3);
t758 = -t771 * pkin(3) - t772 * pkin(11);
t781 = t783 ^ 2;
t697 = -t781 * pkin(3) + t773 * pkin(11) + t771 * t758 + t704;
t709 = -t801 * t733 + t803 * t740;
t756 = t771 * qJD(3) + t811 * t792 + t825 * t807;
t699 = (-t771 * t783 - t756) * pkin(11) + (t772 * t783 - t755) * pkin(3) + t709;
t692 = -t806 * t697 + t858 * t699;
t716 = -t762 * qJD(4) + t858 * t756 + t806 * t773;
t763 = t858 * t772 + t806 * t783;
t736 = t762 * mrSges(5,1) + t763 * mrSges(5,2);
t744 = t762 * mrSges(6,1) - t769 * mrSges(6,3);
t746 = -t769 * mrSges(5,2) - t762 * mrSges(5,3);
t754 = qJDD(4) - t755;
t735 = t762 * pkin(4) - t763 * qJ(5);
t768 = t769 ^ 2;
t690 = -t754 * pkin(4) - t768 * qJ(5) + t763 * t735 + qJDD(5) - t692;
t685 = (t762 * t763 - t754) * pkin(12) + (t716 + t850) * pkin(5) + t690;
t715 = t763 * qJD(4) + t806 * t756 - t858 * t773;
t748 = t763 * pkin(5) - t769 * pkin(12);
t761 = t762 ^ 2;
t696 = -t773 * pkin(3) - t781 * pkin(11) + t772 * t758 - t703;
t815 = (-t716 + t850) * qJ(5) + t696 + (t769 * pkin(4) + t859) * t763;
t688 = -t761 * pkin(5) - t763 * t748 + (pkin(4) + pkin(12)) * t715 + t815;
t805 = sin(qJ(6));
t810 = cos(qJ(6));
t683 = t810 * t685 - t805 * t688;
t742 = t810 * t762 - t805 * t769;
t702 = t742 * qJD(6) + t805 * t715 + t810 * t754;
t743 = t805 * t762 + t810 * t769;
t710 = -mrSges(7,1) * t742 + mrSges(7,2) * t743;
t714 = qJDD(6) + t716;
t760 = qJD(6) + t763;
t719 = -t760 * mrSges(7,2) + t742 * mrSges(7,3);
t681 = m(7) * t683 + t714 * mrSges(7,1) - t702 * mrSges(7,3) - t743 * t710 + t760 * t719;
t684 = t805 * t685 + t810 * t688;
t701 = -t743 * qJD(6) + t810 * t715 - t805 * t754;
t720 = t760 * mrSges(7,1) - t743 * mrSges(7,3);
t682 = m(7) * t684 - t714 * mrSges(7,2) + t701 * mrSges(7,3) + t742 * t710 - t760 * t720;
t673 = t810 * t681 + t805 * t682;
t737 = -t762 * mrSges(6,2) - t763 * mrSges(6,3);
t819 = -m(6) * t690 - t716 * mrSges(6,1) - t763 * t737 - t673;
t671 = m(5) * t692 - t716 * mrSges(5,3) - t763 * t736 + (-t744 + t746) * t769 + (mrSges(5,1) - mrSges(6,2)) * t754 + t819;
t693 = t858 * t697 + t806 * t699;
t747 = t769 * mrSges(5,1) - t763 * mrSges(5,3);
t818 = -t768 * pkin(4) + t754 * qJ(5) - t762 * t735 + t693;
t689 = t769 * t859 - t818;
t745 = t763 * mrSges(6,1) + t769 * mrSges(6,2);
t687 = -t715 * pkin(5) - t761 * pkin(12) + ((2 * qJD(5)) + t748) * t769 + t818;
t820 = -m(7) * t687 + t701 * mrSges(7,1) - t702 * mrSges(7,2) + t742 * t719 - t743 * t720;
t817 = -m(6) * t689 + t754 * mrSges(6,3) + t769 * t745 - t820;
t678 = m(5) * t693 - t754 * mrSges(5,2) - t769 * t747 + (-t736 - t737) * t762 + (-mrSges(5,3) - mrSges(6,1)) * t715 + t817;
t829 = -t806 * t671 + t858 * t678;
t663 = m(4) * t704 - t773 * mrSges(4,2) + t755 * mrSges(4,3) + t771 * t757 - t783 * t765 + t829;
t666 = t858 * t671 + t806 * t678;
t764 = -t783 * mrSges(4,2) + t771 * mrSges(4,3);
t665 = m(4) * t709 - t755 * mrSges(4,1) + t756 * mrSges(4,2) - t771 * t764 + t772 * t765 + t666;
t691 = t715 * pkin(4) + t815;
t841 = -t805 * t681 + t810 * t682;
t823 = -m(6) * t691 + t715 * mrSges(6,2) + t762 * t744 - t841;
t816 = -m(5) * t696 - t715 * mrSges(5,1) - t762 * t746 + (t745 - t747) * t763 + (-mrSges(5,2) + mrSges(6,3)) * t716 + t823;
t669 = m(4) * t703 + t773 * mrSges(4,1) - t756 * mrSges(4,3) - t772 * t757 + t783 * t764 + t816;
t652 = t803 * t811 * t669 + t663 * t845 - t801 * t665;
t766 = -g(3) * t846 + t828;
t788 = -t799 * mrSges(3,2) + mrSges(3,3) * t831;
t791 = (-mrSges(3,1) * t812 + mrSges(3,2) * t808) * t836;
t648 = m(3) * t766 + t798 * mrSges(3,1) - t792 * mrSges(3,3) + t799 * t788 - t791 * t832 + t652;
t651 = t663 * t849 + t803 * t665 + t669 * t848;
t777 = -t802 * t789 - t854;
t787 = t799 * mrSges(3,1) - mrSges(3,3) * t832;
t650 = m(3) * t777 - t793 * mrSges(3,1) + t792 * mrSges(3,2) + (t787 * t808 - t788 * t812) * t836 + t651;
t659 = t811 * t663 - t807 * t669;
t767 = -g(3) * t847 + t837;
t658 = m(3) * t767 - t798 * mrSges(3,2) + t793 * mrSges(3,3) - t799 * t787 + t791 * t831 + t659;
t638 = t648 * t843 - t802 * t650 + t658 * t844;
t636 = m(2) * t796 + qJDD(1) * mrSges(2,1) - t814 * mrSges(2,2) + t638;
t644 = -t808 * t648 + t812 * t658;
t643 = m(2) * t797 - t814 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t644;
t842 = t813 * t636 + t809 * t643;
t840 = t851 * t762 - t852 * t763 + t860 * t769;
t839 = t861 * t762 - t853 * t763 - t851 * t769;
t838 = -t853 * t762 + t862 * t763 + t852 * t769;
t637 = t648 * t846 + t804 * t650 + t658 * t847;
t830 = -t809 * t636 + t813 * t643;
t672 = -t716 * mrSges(6,3) - t763 * t745 - t823;
t705 = Ifges(7,5) * t743 + Ifges(7,6) * t742 + Ifges(7,3) * t760;
t707 = Ifges(7,1) * t743 + Ifges(7,4) * t742 + Ifges(7,5) * t760;
t674 = -mrSges(7,1) * t687 + mrSges(7,3) * t684 + Ifges(7,4) * t702 + Ifges(7,2) * t701 + Ifges(7,6) * t714 - t743 * t705 + t760 * t707;
t706 = Ifges(7,4) * t743 + Ifges(7,2) * t742 + Ifges(7,6) * t760;
t675 = mrSges(7,2) * t687 - mrSges(7,3) * t683 + Ifges(7,1) * t702 + Ifges(7,4) * t701 + Ifges(7,5) * t714 + t742 * t705 - t760 * t706;
t653 = -mrSges(5,1) * t696 - mrSges(6,1) * t689 + mrSges(6,2) * t691 + mrSges(5,3) * t693 - pkin(4) * t672 - pkin(5) * t820 - pkin(12) * t841 - t810 * t674 - t805 * t675 - t861 * t715 + t853 * t716 + t851 * t754 + t840 * t763 + t838 * t769;
t654 = mrSges(6,1) * t690 + mrSges(7,1) * t683 + mrSges(5,2) * t696 - mrSges(7,2) * t684 - mrSges(5,3) * t692 - mrSges(6,3) * t691 + Ifges(7,5) * t702 + Ifges(7,6) * t701 + Ifges(7,3) * t714 + pkin(5) * t673 - qJ(5) * t672 + t743 * t706 - t742 * t707 + t839 * t769 + t840 * t762 + t852 * t754 + t862 * t716 - t853 * t715;
t751 = Ifges(4,4) * t772 + Ifges(4,2) * t771 + Ifges(4,6) * t783;
t752 = Ifges(4,1) * t772 + Ifges(4,4) * t771 + Ifges(4,5) * t783;
t639 = mrSges(4,1) * t703 - mrSges(4,2) * t704 + Ifges(4,5) * t756 + Ifges(4,6) * t755 + Ifges(4,3) * t773 + pkin(3) * t816 + pkin(11) * t829 + t858 * t653 + t806 * t654 + t772 * t751 - t771 * t752;
t774 = Ifges(3,3) * t799 + (Ifges(3,5) * t808 + Ifges(3,6) * t812) * t836;
t776 = Ifges(3,5) * t799 + (Ifges(3,1) * t808 + Ifges(3,4) * t812) * t836;
t750 = Ifges(4,5) * t772 + Ifges(4,6) * t771 + Ifges(4,3) * t783;
t640 = mrSges(4,2) * t709 - mrSges(4,3) * t703 + Ifges(4,1) * t756 + Ifges(4,4) * t755 + Ifges(4,5) * t773 - pkin(11) * t666 - t806 * t653 + t858 * t654 + t771 * t750 - t783 * t751;
t645 = -mrSges(4,1) * t709 + mrSges(4,3) * t704 + mrSges(5,2) * t693 - mrSges(5,1) * t692 + mrSges(6,3) * t689 - mrSges(6,2) * t690 + t783 * t752 - t772 * t750 + Ifges(4,6) * t773 + Ifges(4,2) * t755 + Ifges(4,4) * t756 + pkin(12) * t673 - pkin(3) * t666 - qJ(5) * t817 + t805 * t674 - t810 * t675 - pkin(4) * (-t769 * t744 + t819) + t839 * t763 + (qJ(5) * t737 - t838) * t762 + (pkin(4) * mrSges(6,2) + t860) * t754 - t852 * t716 + (qJ(5) * mrSges(6,1) + t851) * t715;
t821 = pkin(10) * t659 + t640 * t807 + t645 * t811;
t633 = -mrSges(3,1) * t777 + mrSges(3,3) * t767 + Ifges(3,4) * t792 + Ifges(3,2) * t793 + Ifges(3,6) * t798 - pkin(2) * t651 - t801 * t639 - t774 * t832 + t799 * t776 + t821 * t803;
t775 = Ifges(3,6) * t799 + (Ifges(3,4) * t808 + Ifges(3,2) * t812) * t836;
t634 = t774 * t831 + mrSges(3,2) * t777 - mrSges(3,3) * t766 + Ifges(3,1) * t792 + Ifges(3,4) * t793 + Ifges(3,5) * t798 + t811 * t640 - t807 * t645 - t799 * t775 + (-t651 * t801 - t652 * t803) * pkin(10);
t822 = pkin(9) * t644 + t633 * t812 + t634 * t808;
t632 = mrSges(3,1) * t766 - mrSges(3,2) * t767 + Ifges(3,5) * t792 + Ifges(3,6) * t793 + Ifges(3,3) * t798 + pkin(2) * t652 + t803 * t639 + (t775 * t808 - t776 * t812) * t836 + t821 * t801;
t631 = -mrSges(2,2) * g(3) - mrSges(2,3) * t796 + Ifges(2,5) * qJDD(1) - t814 * Ifges(2,6) - t808 * t633 + t812 * t634 + (-t637 * t802 - t638 * t804) * pkin(9);
t630 = mrSges(2,1) * g(3) + mrSges(2,3) * t797 + t814 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t637 - t802 * t632 + t822 * t804;
t1 = [-m(1) * g(1) + t830; -m(1) * g(2) + t842; (-m(1) - m(2)) * g(3) + t637; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t842 - t809 * t630 + t813 * t631; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t830 + t813 * t630 + t809 * t631; -mrSges(1,1) * g(2) + mrSges(2,1) * t796 + mrSges(1,2) * g(1) - mrSges(2,2) * t797 + Ifges(2,3) * qJDD(1) + pkin(1) * t638 + t804 * t632 + t822 * t802;];
tauB  = t1;
