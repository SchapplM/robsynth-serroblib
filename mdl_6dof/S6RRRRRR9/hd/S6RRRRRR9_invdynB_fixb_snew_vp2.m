% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 16:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 15:50:50
% EndTime: 2019-05-08 15:55:10
% DurationCPUTime: 153.41s
% Computational Cost: add. (2446005->417), mult. (6016711->555), div. (0->0), fcn. (5122475->16), ass. (0->181)
t802 = cos(pkin(6));
t797 = qJD(1) * t802 + qJD(2);
t799 = sin(pkin(7));
t801 = cos(pkin(7));
t800 = sin(pkin(6));
t813 = cos(qJ(2));
t834 = qJD(1) * t813;
t830 = t800 * t834;
t782 = (t797 * t799 + t801 * t830) * pkin(10);
t807 = sin(qJ(2));
t836 = qJD(1) * t800;
t849 = pkin(10) * t799;
t786 = (-pkin(2) * t813 - t807 * t849) * t836;
t833 = qJD(1) * qJD(2);
t792 = (qJDD(1) * t807 + t813 * t833) * t800;
t796 = qJDD(1) * t802 + qJDD(2);
t808 = sin(qJ(1));
t814 = cos(qJ(1));
t794 = t808 * g(1) - g(2) * t814;
t815 = qJD(1) ^ 2;
t850 = pkin(9) * t800;
t789 = qJDD(1) * pkin(1) + t815 * t850 + t794;
t795 = -g(1) * t814 - g(2) * t808;
t790 = -pkin(1) * t815 + qJDD(1) * t850 + t795;
t839 = t802 * t813;
t825 = t789 * t839 - t807 * t790;
t835 = qJD(1) * t807;
t848 = pkin(10) * t801;
t738 = -t792 * t848 + t796 * pkin(2) + t797 * t782 + (-g(3) * t813 - t786 * t835) * t800 + t825;
t831 = t800 * t835;
t785 = pkin(2) * t797 - t831 * t848;
t793 = (qJDD(1) * t813 - t807 * t833) * t800;
t823 = t793 * t801 + t796 * t799;
t840 = t802 * t807;
t837 = t789 * t840 + t813 * t790;
t739 = -t797 * t785 + (-g(3) * t807 + t786 * t834) * t800 + t823 * pkin(10) + t837;
t847 = t802 * g(3);
t744 = -t792 * t849 - t793 * pkin(2) - t847 + (-t789 + (-t782 * t813 + t785 * t807) * qJD(1)) * t800;
t806 = sin(qJ(3));
t812 = cos(qJ(3));
t710 = -t806 * t739 + (t738 * t801 + t744 * t799) * t812;
t841 = t801 * t813;
t846 = t799 * t806;
t773 = t797 * t846 + (t806 * t841 + t807 * t812) * t836;
t756 = -t773 * qJD(3) - t806 * t792 + t812 * t823;
t845 = t799 * t812;
t772 = (-t806 * t807 + t812 * t841) * t836 + t797 * t845;
t844 = t800 * t807;
t843 = t800 * t813;
t842 = t801 * t806;
t711 = t738 * t842 + t812 * t739 + t744 * t846;
t758 = -mrSges(4,1) * t772 + mrSges(4,2) * t773;
t783 = t797 * t801 - t799 * t830 + qJD(3);
t766 = mrSges(4,1) * t783 - mrSges(4,3) * t773;
t774 = -t793 * t799 + t796 * t801 + qJDD(3);
t759 = -pkin(3) * t772 - pkin(11) * t773;
t781 = t783 ^ 2;
t699 = -pkin(3) * t781 + pkin(11) * t774 + t759 * t772 + t711;
t718 = -t799 * t738 + t801 * t744;
t757 = t772 * qJD(3) + t812 * t792 + t806 * t823;
t701 = (-t772 * t783 - t757) * pkin(11) + (t773 * t783 - t756) * pkin(3) + t718;
t805 = sin(qJ(4));
t811 = cos(qJ(4));
t692 = t811 * t699 + t805 * t701;
t764 = t773 * t811 + t783 * t805;
t725 = -t764 * qJD(4) - t805 * t757 + t774 * t811;
t763 = -t805 * t773 + t783 * t811;
t740 = -mrSges(5,1) * t763 + mrSges(5,2) * t764;
t771 = qJD(4) - t772;
t750 = mrSges(5,1) * t771 - mrSges(5,3) * t764;
t755 = qJDD(4) - t756;
t741 = -pkin(4) * t763 - pkin(12) * t764;
t770 = t771 ^ 2;
t687 = -pkin(4) * t770 + pkin(12) * t755 + t741 * t763 + t692;
t698 = -t774 * pkin(3) - t781 * pkin(11) + t773 * t759 - t710;
t726 = qJD(4) * t763 + t757 * t811 + t774 * t805;
t690 = (-t763 * t771 - t726) * pkin(12) + (t764 * t771 - t725) * pkin(4) + t698;
t804 = sin(qJ(5));
t810 = cos(qJ(5));
t682 = -t804 * t687 + t810 * t690;
t747 = -t764 * t804 + t771 * t810;
t709 = qJD(5) * t747 + t726 * t810 + t755 * t804;
t724 = qJDD(5) - t725;
t748 = t764 * t810 + t771 * t804;
t762 = qJD(5) - t763;
t680 = (t747 * t762 - t709) * pkin(13) + (t747 * t748 + t724) * pkin(5) + t682;
t683 = t810 * t687 + t804 * t690;
t708 = -qJD(5) * t748 - t726 * t804 + t755 * t810;
t730 = pkin(5) * t762 - pkin(13) * t748;
t746 = t747 ^ 2;
t681 = -pkin(5) * t746 + pkin(13) * t708 - t730 * t762 + t683;
t803 = sin(qJ(6));
t809 = cos(qJ(6));
t678 = t680 * t809 - t681 * t803;
t719 = t747 * t809 - t748 * t803;
t695 = qJD(6) * t719 + t708 * t803 + t709 * t809;
t720 = t747 * t803 + t748 * t809;
t706 = -mrSges(7,1) * t719 + mrSges(7,2) * t720;
t760 = qJD(6) + t762;
t712 = -mrSges(7,2) * t760 + mrSges(7,3) * t719;
t722 = qJDD(6) + t724;
t676 = m(7) * t678 + mrSges(7,1) * t722 - mrSges(7,3) * t695 - t706 * t720 + t712 * t760;
t679 = t680 * t803 + t681 * t809;
t694 = -qJD(6) * t720 + t708 * t809 - t709 * t803;
t713 = mrSges(7,1) * t760 - mrSges(7,3) * t720;
t677 = m(7) * t679 - mrSges(7,2) * t722 + mrSges(7,3) * t694 + t706 * t719 - t713 * t760;
t668 = t809 * t676 + t803 * t677;
t721 = -mrSges(6,1) * t747 + mrSges(6,2) * t748;
t728 = -mrSges(6,2) * t762 + mrSges(6,3) * t747;
t666 = m(6) * t682 + mrSges(6,1) * t724 - mrSges(6,3) * t709 - t721 * t748 + t728 * t762 + t668;
t729 = mrSges(6,1) * t762 - mrSges(6,3) * t748;
t826 = -t676 * t803 + t809 * t677;
t667 = m(6) * t683 - mrSges(6,2) * t724 + mrSges(6,3) * t708 + t721 * t747 - t729 * t762 + t826;
t827 = -t666 * t804 + t810 * t667;
t663 = m(5) * t692 - mrSges(5,2) * t755 + mrSges(5,3) * t725 + t740 * t763 - t750 * t771 + t827;
t691 = -t805 * t699 + t701 * t811;
t749 = -mrSges(5,2) * t771 + mrSges(5,3) * t763;
t686 = -pkin(4) * t755 - pkin(12) * t770 + t764 * t741 - t691;
t684 = -pkin(5) * t708 - pkin(13) * t746 + t730 * t748 + t686;
t818 = m(7) * t684 - t694 * mrSges(7,1) + mrSges(7,2) * t695 - t719 * t712 + t713 * t720;
t816 = -m(6) * t686 + t708 * mrSges(6,1) - mrSges(6,2) * t709 + t747 * t728 - t729 * t748 - t818;
t672 = m(5) * t691 + mrSges(5,1) * t755 - mrSges(5,3) * t726 - t740 * t764 + t749 * t771 + t816;
t828 = t811 * t663 - t672 * t805;
t652 = m(4) * t711 - mrSges(4,2) * t774 + mrSges(4,3) * t756 + t758 * t772 - t766 * t783 + t828;
t655 = t805 * t663 + t811 * t672;
t765 = -mrSges(4,2) * t783 + mrSges(4,3) * t772;
t654 = m(4) * t718 - mrSges(4,1) * t756 + mrSges(4,2) * t757 - t765 * t772 + t766 * t773 + t655;
t664 = t666 * t810 + t667 * t804;
t817 = -m(5) * t698 + t725 * mrSges(5,1) - mrSges(5,2) * t726 + t763 * t749 - t750 * t764 - t664;
t660 = m(4) * t710 + mrSges(4,1) * t774 - mrSges(4,3) * t757 - t758 * t773 + t765 * t783 + t817;
t641 = t801 * t812 * t660 + t652 * t842 - t654 * t799;
t767 = -g(3) * t843 + t825;
t788 = -mrSges(3,2) * t797 + mrSges(3,3) * t830;
t791 = (-mrSges(3,1) * t813 + mrSges(3,2) * t807) * t836;
t637 = m(3) * t767 + mrSges(3,1) * t796 - mrSges(3,3) * t792 + t788 * t797 - t791 * t831 + t641;
t640 = t652 * t846 + t801 * t654 + t660 * t845;
t778 = -t800 * t789 - t847;
t787 = mrSges(3,1) * t797 - mrSges(3,3) * t831;
t639 = m(3) * t778 - t793 * mrSges(3,1) + t792 * mrSges(3,2) + (t787 * t807 - t788 * t813) * t836 + t640;
t647 = t812 * t652 - t660 * t806;
t768 = -g(3) * t844 + t837;
t646 = m(3) * t768 - mrSges(3,2) * t796 + mrSges(3,3) * t793 - t787 * t797 + t791 * t830 + t647;
t627 = t637 * t839 - t639 * t800 + t646 * t840;
t625 = m(2) * t794 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t815 + t627;
t633 = -t637 * t807 + t813 * t646;
t632 = m(2) * t795 - mrSges(2,1) * t815 - qJDD(1) * mrSges(2,2) + t633;
t838 = t814 * t625 + t808 * t632;
t626 = t637 * t843 + t802 * t639 + t646 * t844;
t829 = -t625 * t808 + t814 * t632;
t702 = Ifges(7,5) * t720 + Ifges(7,6) * t719 + Ifges(7,3) * t760;
t704 = Ifges(7,1) * t720 + Ifges(7,4) * t719 + Ifges(7,5) * t760;
t669 = -mrSges(7,1) * t684 + mrSges(7,3) * t679 + Ifges(7,4) * t695 + Ifges(7,2) * t694 + Ifges(7,6) * t722 - t702 * t720 + t704 * t760;
t703 = Ifges(7,4) * t720 + Ifges(7,2) * t719 + Ifges(7,6) * t760;
t670 = mrSges(7,2) * t684 - mrSges(7,3) * t678 + Ifges(7,1) * t695 + Ifges(7,4) * t694 + Ifges(7,5) * t722 + t702 * t719 - t703 * t760;
t714 = Ifges(6,5) * t748 + Ifges(6,6) * t747 + Ifges(6,3) * t762;
t716 = Ifges(6,1) * t748 + Ifges(6,4) * t747 + Ifges(6,5) * t762;
t656 = -mrSges(6,1) * t686 + mrSges(6,3) * t683 + Ifges(6,4) * t709 + Ifges(6,2) * t708 + Ifges(6,6) * t724 - pkin(5) * t818 + pkin(13) * t826 + t809 * t669 + t803 * t670 - t748 * t714 + t762 * t716;
t715 = Ifges(6,4) * t748 + Ifges(6,2) * t747 + Ifges(6,6) * t762;
t657 = mrSges(6,2) * t686 - mrSges(6,3) * t682 + Ifges(6,1) * t709 + Ifges(6,4) * t708 + Ifges(6,5) * t724 - pkin(13) * t668 - t669 * t803 + t670 * t809 + t714 * t747 - t715 * t762;
t731 = Ifges(5,5) * t764 + Ifges(5,6) * t763 + Ifges(5,3) * t771;
t732 = Ifges(5,4) * t764 + Ifges(5,2) * t763 + Ifges(5,6) * t771;
t642 = mrSges(5,2) * t698 - mrSges(5,3) * t691 + Ifges(5,1) * t726 + Ifges(5,4) * t725 + Ifges(5,5) * t755 - pkin(12) * t664 - t656 * t804 + t657 * t810 + t731 * t763 - t732 * t771;
t733 = Ifges(5,1) * t764 + Ifges(5,4) * t763 + Ifges(5,5) * t771;
t648 = Ifges(5,4) * t726 + Ifges(5,2) * t725 + Ifges(5,6) * t755 - t764 * t731 + t771 * t733 - mrSges(5,1) * t698 + mrSges(5,3) * t692 - Ifges(6,5) * t709 - Ifges(6,6) * t708 - Ifges(6,3) * t724 - t748 * t715 + t747 * t716 - mrSges(6,1) * t682 + mrSges(6,2) * t683 - Ifges(7,5) * t695 - Ifges(7,6) * t694 - Ifges(7,3) * t722 - t720 * t703 + t719 * t704 - mrSges(7,1) * t678 + mrSges(7,2) * t679 - pkin(5) * t668 - pkin(4) * t664;
t752 = Ifges(4,4) * t773 + Ifges(4,2) * t772 + Ifges(4,6) * t783;
t753 = Ifges(4,1) * t773 + Ifges(4,4) * t772 + Ifges(4,5) * t783;
t628 = mrSges(4,1) * t710 - mrSges(4,2) * t711 + Ifges(4,5) * t757 + Ifges(4,6) * t756 + Ifges(4,3) * t774 + pkin(3) * t817 + pkin(11) * t828 + t805 * t642 + t811 * t648 + t773 * t752 - t772 * t753;
t775 = Ifges(3,3) * t797 + (Ifges(3,5) * t807 + Ifges(3,6) * t813) * t836;
t777 = Ifges(3,5) * t797 + (Ifges(3,1) * t807 + Ifges(3,4) * t813) * t836;
t751 = Ifges(4,5) * t773 + Ifges(4,6) * t772 + Ifges(4,3) * t783;
t629 = mrSges(4,2) * t718 - mrSges(4,3) * t710 + Ifges(4,1) * t757 + Ifges(4,4) * t756 + Ifges(4,5) * t774 - pkin(11) * t655 + t642 * t811 - t648 * t805 + t751 * t772 - t752 * t783;
t634 = Ifges(4,4) * t757 + Ifges(4,2) * t756 + Ifges(4,6) * t774 - t773 * t751 + t783 * t753 - mrSges(4,1) * t718 + mrSges(4,3) * t711 - Ifges(5,5) * t726 - Ifges(5,6) * t725 - Ifges(5,3) * t755 - t764 * t732 + t763 * t733 - mrSges(5,1) * t691 + mrSges(5,2) * t692 - t804 * t657 - t810 * t656 - pkin(4) * t816 - pkin(12) * t827 - pkin(3) * t655;
t819 = pkin(10) * t647 + t629 * t806 + t634 * t812;
t622 = -mrSges(3,1) * t778 + mrSges(3,3) * t768 + Ifges(3,4) * t792 + Ifges(3,2) * t793 + Ifges(3,6) * t796 - pkin(2) * t640 - t799 * t628 - t775 * t831 + t797 * t777 + t801 * t819;
t776 = Ifges(3,6) * t797 + (Ifges(3,4) * t807 + Ifges(3,2) * t813) * t836;
t623 = t775 * t830 + mrSges(3,2) * t778 - mrSges(3,3) * t767 + Ifges(3,1) * t792 + Ifges(3,4) * t793 + Ifges(3,5) * t796 + t812 * t629 - t806 * t634 - t797 * t776 + (-t640 * t799 - t641 * t801) * pkin(10);
t820 = pkin(9) * t633 + t622 * t813 + t623 * t807;
t621 = mrSges(3,1) * t767 - mrSges(3,2) * t768 + Ifges(3,5) * t792 + Ifges(3,6) * t793 + Ifges(3,3) * t796 + pkin(2) * t641 + t801 * t628 + (t776 * t807 - t777 * t813) * t836 + t819 * t799;
t620 = -mrSges(2,2) * g(3) - mrSges(2,3) * t794 + Ifges(2,5) * qJDD(1) - t815 * Ifges(2,6) - t807 * t622 + t813 * t623 + (-t626 * t800 - t627 * t802) * pkin(9);
t619 = mrSges(2,1) * g(3) + mrSges(2,3) * t795 + t815 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t626 - t800 * t621 + t802 * t820;
t1 = [-m(1) * g(1) + t829; -m(1) * g(2) + t838; (-m(1) - m(2)) * g(3) + t626; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t838 - t808 * t619 + t814 * t620; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t829 + t814 * t619 + t808 * t620; -mrSges(1,1) * g(2) + mrSges(2,1) * t794 + mrSges(1,2) * g(1) - mrSges(2,2) * t795 + Ifges(2,3) * qJDD(1) + pkin(1) * t627 + t802 * t621 + t800 * t820;];
tauB  = t1;
