% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:38:15
% EndTime: 2019-05-07 15:40:45
% DurationCPUTime: 136.97s
% Computational Cost: add. (2205010->417), mult. (5518599->554), div. (0->0), fcn. (4676567->16), ass. (0->180)
t848 = cos(qJ(3));
t802 = sin(pkin(6));
t847 = pkin(9) * t802;
t801 = sin(pkin(7));
t846 = pkin(10) * t801;
t804 = cos(pkin(7));
t845 = pkin(10) * t804;
t805 = cos(pkin(6));
t844 = t805 * g(3);
t808 = sin(qJ(3));
t843 = t801 * t808;
t809 = sin(qJ(2));
t842 = t802 * t809;
t813 = cos(qJ(2));
t841 = t802 * t813;
t840 = t804 * t808;
t839 = t805 * t809;
t838 = t805 * t813;
t798 = qJD(1) * t805 + qJD(2);
t833 = qJD(1) * t813;
t828 = t802 * t833;
t822 = t804 * t828;
t779 = (t798 * t801 + t822) * pkin(10);
t835 = qJD(1) * t802;
t784 = (-pkin(2) * t813 - t809 * t846) * t835;
t832 = qJD(1) * qJD(2);
t790 = (qJDD(1) * t809 + t813 * t832) * t802;
t797 = qJDD(1) * t805 + qJDD(2);
t810 = sin(qJ(1));
t814 = cos(qJ(1));
t795 = t810 * g(1) - g(2) * t814;
t815 = qJD(1) ^ 2;
t787 = qJDD(1) * pkin(1) + t815 * t847 + t795;
t796 = -g(1) * t814 - g(2) * t810;
t788 = -pkin(1) * t815 + qJDD(1) * t847 + t796;
t823 = t787 * t838 - t809 * t788;
t834 = qJD(1) * t809;
t737 = -t790 * t845 + t797 * pkin(2) + t798 * t779 + (-g(3) * t813 - t784 * t834) * t802 + t823;
t829 = t802 * t834;
t783 = pkin(2) * t798 - t829 * t845;
t791 = (qJDD(1) * t813 - t809 * t832) * t802;
t821 = t791 * t804 + t797 * t801;
t836 = t787 * t839 + t813 * t788;
t738 = -t798 * t783 + (-g(3) * t809 + t784 * t833) * t802 + t821 * pkin(10) + t836;
t745 = -t790 * t846 - t791 * pkin(2) - t844 + (-t787 + (-t779 * t813 + t783 * t809) * qJD(1)) * t802;
t711 = t737 * t840 + t848 * t738 + t745 * t843;
t770 = t798 * t843 + (t848 * t809 + t813 * t840) * t835;
t830 = t804 * t848;
t831 = t801 * t848;
t754 = qJD(3) * t770 + t790 * t808 - t791 * t830 - t797 * t831;
t769 = -t798 * t831 + t808 * t829 - t848 * t822;
t757 = mrSges(4,1) * t769 + mrSges(4,2) * t770;
t781 = t798 * t804 - t801 * t828 + qJD(3);
t764 = mrSges(4,1) * t781 - mrSges(4,3) * t770;
t771 = -t791 * t801 + t797 * t804 + qJDD(3);
t756 = pkin(3) * t769 - qJ(4) * t770;
t778 = t781 ^ 2;
t698 = -pkin(3) * t778 + qJ(4) * t771 - t756 * t769 + t711;
t721 = -t801 * t737 + t804 * t745;
t755 = -t769 * qJD(3) + t848 * t790 + t821 * t808;
t701 = (t769 * t781 - t755) * qJ(4) + (t770 * t781 + t754) * pkin(3) + t721;
t800 = sin(pkin(13));
t803 = cos(pkin(13));
t762 = t770 * t803 + t781 * t800;
t690 = -0.2e1 * qJD(4) * t762 - t800 * t698 + t803 * t701;
t742 = t755 * t803 + t771 * t800;
t761 = -t770 * t800 + t781 * t803;
t687 = (t761 * t769 - t742) * pkin(11) + (t761 * t762 + t754) * pkin(4) + t690;
t691 = 0.2e1 * qJD(4) * t761 + t803 * t698 + t800 * t701;
t741 = -t755 * t800 + t771 * t803;
t749 = pkin(4) * t769 - pkin(11) * t762;
t760 = t761 ^ 2;
t689 = -pkin(4) * t760 + pkin(11) * t741 - t749 * t769 + t691;
t807 = sin(qJ(5));
t812 = cos(qJ(5));
t684 = t807 * t687 + t812 * t689;
t736 = t761 * t807 + t762 * t812;
t708 = -qJD(5) * t736 + t741 * t812 - t742 * t807;
t735 = t761 * t812 - t762 * t807;
t719 = -mrSges(6,1) * t735 + mrSges(6,2) * t736;
t768 = qJD(5) + t769;
t725 = mrSges(6,1) * t768 - mrSges(6,3) * t736;
t753 = qJDD(5) + t754;
t720 = -pkin(5) * t735 - pkin(12) * t736;
t767 = t768 ^ 2;
t682 = -pkin(5) * t767 + pkin(12) * t753 + t720 * t735 + t684;
t710 = t737 * t830 - t808 * t738 + t745 * t831;
t697 = -t771 * pkin(3) - t778 * qJ(4) + t770 * t756 + qJDD(4) - t710;
t692 = -t741 * pkin(4) - t760 * pkin(11) + t762 * t749 + t697;
t709 = qJD(5) * t735 + t741 * t807 + t742 * t812;
t685 = (-t735 * t768 - t709) * pkin(12) + (t736 * t768 - t708) * pkin(5) + t692;
t806 = sin(qJ(6));
t811 = cos(qJ(6));
t679 = -t682 * t806 + t685 * t811;
t722 = -t736 * t806 + t768 * t811;
t695 = qJD(6) * t722 + t709 * t811 + t753 * t806;
t707 = qJDD(6) - t708;
t723 = t736 * t811 + t768 * t806;
t712 = -mrSges(7,1) * t722 + mrSges(7,2) * t723;
t734 = qJD(6) - t735;
t713 = -mrSges(7,2) * t734 + mrSges(7,3) * t722;
t677 = m(7) * t679 + mrSges(7,1) * t707 - mrSges(7,3) * t695 - t712 * t723 + t713 * t734;
t680 = t682 * t811 + t685 * t806;
t694 = -qJD(6) * t723 - t709 * t806 + t753 * t811;
t714 = mrSges(7,1) * t734 - mrSges(7,3) * t723;
t678 = m(7) * t680 - mrSges(7,2) * t707 + mrSges(7,3) * t694 + t712 * t722 - t714 * t734;
t824 = -t677 * t806 + t811 * t678;
t665 = m(6) * t684 - mrSges(6,2) * t753 + mrSges(6,3) * t708 + t719 * t735 - t725 * t768 + t824;
t683 = t687 * t812 - t689 * t807;
t724 = -mrSges(6,2) * t768 + mrSges(6,3) * t735;
t681 = -pkin(5) * t753 - pkin(12) * t767 + t720 * t736 - t683;
t818 = -m(7) * t681 + t694 * mrSges(7,1) - mrSges(7,2) * t695 + t722 * t713 - t714 * t723;
t673 = m(6) * t683 + mrSges(6,1) * t753 - mrSges(6,3) * t709 - t719 * t736 + t724 * t768 + t818;
t662 = t807 * t665 + t812 * t673;
t740 = -mrSges(5,1) * t761 + mrSges(5,2) * t762;
t747 = -mrSges(5,2) * t769 + mrSges(5,3) * t761;
t660 = m(5) * t690 + mrSges(5,1) * t754 - mrSges(5,3) * t742 - t740 * t762 + t747 * t769 + t662;
t748 = mrSges(5,1) * t769 - mrSges(5,3) * t762;
t825 = t812 * t665 - t673 * t807;
t661 = m(5) * t691 - mrSges(5,2) * t754 + mrSges(5,3) * t741 + t740 * t761 - t748 * t769 + t825;
t826 = -t660 * t800 + t803 * t661;
t651 = m(4) * t711 - mrSges(4,2) * t771 - mrSges(4,3) * t754 - t757 * t769 - t764 * t781 + t826;
t654 = t803 * t660 + t800 * t661;
t763 = -mrSges(4,2) * t781 - mrSges(4,3) * t769;
t653 = m(4) * t721 + mrSges(4,1) * t754 + mrSges(4,2) * t755 + t763 * t769 + t764 * t770 + t654;
t669 = t811 * t677 + t806 * t678;
t817 = m(6) * t692 - t708 * mrSges(6,1) + mrSges(6,2) * t709 - t735 * t724 + t725 * t736 + t669;
t816 = -m(5) * t697 + t741 * mrSges(5,1) - mrSges(5,2) * t742 + t761 * t747 - t748 * t762 - t817;
t668 = m(4) * t710 + mrSges(4,1) * t771 - mrSges(4,3) * t755 - t757 * t770 + t763 * t781 + t816;
t640 = t651 * t840 - t653 * t801 + t668 * t830;
t765 = -g(3) * t841 + t823;
t786 = -mrSges(3,2) * t798 + mrSges(3,3) * t828;
t789 = (-mrSges(3,1) * t813 + mrSges(3,2) * t809) * t835;
t636 = m(3) * t765 + mrSges(3,1) * t797 - mrSges(3,3) * t790 + t786 * t798 - t789 * t829 + t640;
t639 = t651 * t843 + t804 * t653 + t668 * t831;
t775 = -t802 * t787 - t844;
t785 = mrSges(3,1) * t798 - mrSges(3,3) * t829;
t638 = m(3) * t775 - t791 * mrSges(3,1) + t790 * mrSges(3,2) + (t785 * t809 - t786 * t813) * t835 + t639;
t647 = t848 * t651 - t668 * t808;
t766 = -g(3) * t842 + t836;
t646 = m(3) * t766 - mrSges(3,2) * t797 + mrSges(3,3) * t791 - t785 * t798 + t789 * t828 + t647;
t626 = t636 * t838 - t638 * t802 + t646 * t839;
t624 = m(2) * t795 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t815 + t626;
t632 = -t636 * t809 + t813 * t646;
t631 = m(2) * t796 - mrSges(2,1) * t815 - qJDD(1) * mrSges(2,2) + t632;
t837 = t814 * t624 + t810 * t631;
t625 = t636 * t841 + t805 * t638 + t646 * t842;
t827 = -t624 * t810 + t814 * t631;
t702 = Ifges(7,5) * t723 + Ifges(7,6) * t722 + Ifges(7,3) * t734;
t704 = Ifges(7,1) * t723 + Ifges(7,4) * t722 + Ifges(7,5) * t734;
t670 = -mrSges(7,1) * t681 + mrSges(7,3) * t680 + Ifges(7,4) * t695 + Ifges(7,2) * t694 + Ifges(7,6) * t707 - t702 * t723 + t704 * t734;
t703 = Ifges(7,4) * t723 + Ifges(7,2) * t722 + Ifges(7,6) * t734;
t671 = mrSges(7,2) * t681 - mrSges(7,3) * t679 + Ifges(7,1) * t695 + Ifges(7,4) * t694 + Ifges(7,5) * t707 + t702 * t722 - t703 * t734;
t715 = Ifges(6,5) * t736 + Ifges(6,6) * t735 + Ifges(6,3) * t768;
t716 = Ifges(6,4) * t736 + Ifges(6,2) * t735 + Ifges(6,6) * t768;
t655 = mrSges(6,2) * t692 - mrSges(6,3) * t683 + Ifges(6,1) * t709 + Ifges(6,4) * t708 + Ifges(6,5) * t753 - pkin(12) * t669 - t670 * t806 + t671 * t811 + t715 * t735 - t716 * t768;
t717 = Ifges(6,1) * t736 + Ifges(6,4) * t735 + Ifges(6,5) * t768;
t656 = -mrSges(6,1) * t692 - mrSges(7,1) * t679 + mrSges(7,2) * t680 + mrSges(6,3) * t684 + Ifges(6,4) * t709 - Ifges(7,5) * t695 + Ifges(6,2) * t708 + Ifges(6,6) * t753 - Ifges(7,6) * t694 - Ifges(7,3) * t707 - pkin(5) * t669 - t703 * t723 + t704 * t722 - t715 * t736 + t717 * t768;
t726 = Ifges(5,5) * t762 + Ifges(5,6) * t761 + Ifges(5,3) * t769;
t728 = Ifges(5,1) * t762 + Ifges(5,4) * t761 + Ifges(5,5) * t769;
t641 = -mrSges(5,1) * t697 + mrSges(5,3) * t691 + Ifges(5,4) * t742 + Ifges(5,2) * t741 + Ifges(5,6) * t754 - pkin(4) * t817 + pkin(11) * t825 + t807 * t655 + t812 * t656 - t762 * t726 + t769 * t728;
t727 = Ifges(5,4) * t762 + Ifges(5,2) * t761 + Ifges(5,6) * t769;
t642 = mrSges(5,2) * t697 - mrSges(5,3) * t690 + Ifges(5,1) * t742 + Ifges(5,4) * t741 + Ifges(5,5) * t754 - pkin(11) * t662 + t655 * t812 - t656 * t807 + t726 * t761 - t727 * t769;
t751 = Ifges(4,4) * t770 - Ifges(4,2) * t769 + Ifges(4,6) * t781;
t752 = Ifges(4,1) * t770 - Ifges(4,4) * t769 + Ifges(4,5) * t781;
t627 = mrSges(4,1) * t710 - mrSges(4,2) * t711 + Ifges(4,5) * t755 - Ifges(4,6) * t754 + Ifges(4,3) * t771 + pkin(3) * t816 + qJ(4) * t826 + t803 * t641 + t800 * t642 + t770 * t751 + t769 * t752;
t772 = Ifges(3,3) * t798 + (Ifges(3,5) * t809 + Ifges(3,6) * t813) * t835;
t774 = Ifges(3,5) * t798 + (Ifges(3,1) * t809 + Ifges(3,4) * t813) * t835;
t750 = Ifges(4,5) * t770 - Ifges(4,6) * t769 + Ifges(4,3) * t781;
t628 = mrSges(4,2) * t721 - mrSges(4,3) * t710 + Ifges(4,1) * t755 - Ifges(4,4) * t754 + Ifges(4,5) * t771 - qJ(4) * t654 - t641 * t800 + t642 * t803 - t750 * t769 - t751 * t781;
t633 = (-Ifges(5,3) - Ifges(4,2)) * t754 - t811 * t670 - t806 * t671 + Ifges(4,6) * t771 + t781 * t752 - t770 * t750 + t761 * t728 - t762 * t727 - Ifges(6,3) * t753 + Ifges(4,4) * t755 - Ifges(5,5) * t742 + t735 * t717 - t736 * t716 - Ifges(5,6) * t741 - mrSges(4,1) * t721 - Ifges(6,6) * t708 - Ifges(6,5) * t709 + mrSges(4,3) * t711 - mrSges(5,1) * t690 + mrSges(5,2) * t691 + mrSges(6,2) * t684 - mrSges(6,1) * t683 - pkin(4) * t662 - pkin(5) * t818 - pkin(12) * t824 - pkin(3) * t654;
t819 = pkin(10) * t647 + t628 * t808 + t848 * t633;
t621 = -mrSges(3,1) * t775 + mrSges(3,3) * t766 + Ifges(3,4) * t790 + Ifges(3,2) * t791 + Ifges(3,6) * t797 - pkin(2) * t639 - t801 * t627 - t772 * t829 + t798 * t774 + t804 * t819;
t773 = Ifges(3,6) * t798 + (Ifges(3,4) * t809 + Ifges(3,2) * t813) * t835;
t622 = t772 * t828 + mrSges(3,2) * t775 - mrSges(3,3) * t765 + t848 * t628 + Ifges(3,1) * t790 + Ifges(3,4) * t791 + Ifges(3,5) * t797 - t808 * t633 - t798 * t773 + (-t639 * t801 - t640 * t804) * pkin(10);
t820 = pkin(9) * t632 + t621 * t813 + t622 * t809;
t620 = mrSges(3,1) * t765 - mrSges(3,2) * t766 + Ifges(3,5) * t790 + Ifges(3,6) * t791 + Ifges(3,3) * t797 + pkin(2) * t640 + t804 * t627 + (t773 * t809 - t774 * t813) * t835 + t819 * t801;
t619 = -mrSges(2,2) * g(3) - mrSges(2,3) * t795 + Ifges(2,5) * qJDD(1) - t815 * Ifges(2,6) - t809 * t621 + t813 * t622 + (-t625 * t802 - t626 * t805) * pkin(9);
t618 = mrSges(2,1) * g(3) + mrSges(2,3) * t796 + t815 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t625 - t802 * t620 + t805 * t820;
t1 = [-m(1) * g(1) + t827; -m(1) * g(2) + t837; (-m(1) - m(2)) * g(3) + t625; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t837 - t810 * t618 + t814 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t827 + t814 * t618 + t810 * t619; -mrSges(1,1) * g(2) + mrSges(2,1) * t795 + mrSges(1,2) * g(1) - mrSges(2,2) * t796 + Ifges(2,3) * qJDD(1) + pkin(1) * t626 + t805 * t620 + t802 * t820;];
tauB  = t1;
