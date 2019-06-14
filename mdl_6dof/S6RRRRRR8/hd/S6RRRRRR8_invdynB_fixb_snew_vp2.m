% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR8
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
% Datum: 2019-05-08 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 13:40:47
% EndTime: 2019-05-08 13:44:05
% DurationCPUTime: 143.64s
% Computational Cost: add. (2369386->416), mult. (5809371->551), div. (0->0), fcn. (4956642->16), ass. (0->181)
t812 = cos(pkin(6));
t807 = t812 * qJD(1) + qJD(2);
t809 = sin(pkin(7));
t811 = cos(pkin(7));
t810 = sin(pkin(6));
t823 = cos(qJ(2));
t842 = qJD(1) * t823;
t839 = t810 * t842;
t831 = t807 * t809 + t811 * t839;
t790 = t831 * pkin(10);
t817 = sin(qJ(2));
t844 = qJD(1) * t810;
t857 = pkin(10) * t809;
t795 = (-pkin(2) * t823 - t817 * t857) * t844;
t841 = qJD(1) * qJD(2);
t801 = (qJDD(1) * t817 + t823 * t841) * t810;
t806 = t812 * qJDD(1) + qJDD(2);
t818 = sin(qJ(1));
t824 = cos(qJ(1));
t804 = t818 * g(1) - t824 * g(2);
t825 = qJD(1) ^ 2;
t858 = pkin(9) * t810;
t798 = qJDD(1) * pkin(1) + t825 * t858 + t804;
t805 = -t824 * g(1) - t818 * g(2);
t799 = -t825 * pkin(1) + qJDD(1) * t858 + t805;
t848 = t812 * t823;
t834 = t798 * t848 - t817 * t799;
t843 = qJD(1) * t817;
t856 = pkin(10) * t811;
t748 = -t801 * t856 + t806 * pkin(2) + t807 * t790 + (-g(3) * t823 - t795 * t843) * t810 + t834;
t840 = t810 * t843;
t793 = t807 * pkin(2) - t840 * t856;
t802 = (qJDD(1) * t823 - t817 * t841) * t810;
t832 = t802 * t811 + t806 * t809;
t849 = t812 * t817;
t845 = t798 * t849 + t823 * t799;
t749 = -t807 * t793 + (-g(3) * t817 + t795 * t842) * t810 + t832 * pkin(10) + t845;
t855 = t812 * g(3);
t755 = -t801 * t857 - t802 * pkin(2) - t855 + (-t798 + (-t790 * t823 + t793 * t817) * qJD(1)) * t810;
t816 = sin(qJ(3));
t822 = cos(qJ(3));
t720 = -t816 * t749 + (t748 * t811 + t755 * t809) * t822;
t780 = -t816 * t840 + t822 * t831;
t850 = t811 * t816;
t853 = t809 * t816;
t781 = t807 * t853 + (t817 * t822 + t823 * t850) * t844;
t765 = -t781 * qJD(3) - t816 * t801 + t822 * t832;
t766 = t780 * qJD(3) + t822 * t801 + t816 * t832;
t767 = -t780 * mrSges(4,1) + t781 * mrSges(4,2);
t791 = t811 * t807 - t809 * t839 + qJD(3);
t772 = -t791 * mrSges(4,2) + t780 * mrSges(4,3);
t782 = -t809 * t802 + t811 * t806 + qJDD(3);
t768 = -t780 * pkin(3) - t781 * pkin(11);
t789 = t791 ^ 2;
t707 = -t782 * pkin(3) - t789 * pkin(11) + t781 * t768 - t720;
t815 = sin(qJ(4));
t821 = cos(qJ(4));
t771 = t821 * t781 + t815 * t791;
t733 = -t771 * qJD(4) - t815 * t766 + t821 * t782;
t770 = -t815 * t781 + t821 * t791;
t734 = t770 * qJD(4) + t821 * t766 + t815 * t782;
t779 = qJD(4) - t780;
t757 = -t779 * mrSges(5,2) + t770 * mrSges(5,3);
t758 = t779 * mrSges(5,1) - t771 * mrSges(5,3);
t721 = t748 * t850 + t822 * t749 + t755 * t853;
t708 = -t789 * pkin(3) + t782 * pkin(11) + t780 * t768 + t721;
t731 = -t809 * t748 + t811 * t755;
t711 = (-t780 * t791 - t766) * pkin(11) + (t781 * t791 - t765) * pkin(3) + t731;
t700 = -t815 * t708 + t821 * t711;
t764 = qJDD(4) - t765;
t697 = (t770 * t779 - t734) * pkin(12) + (t770 * t771 + t764) * pkin(4) + t700;
t701 = t821 * t708 + t815 * t711;
t759 = t779 * pkin(4) - t771 * pkin(12);
t769 = t770 ^ 2;
t699 = -t769 * pkin(4) + t733 * pkin(12) - t779 * t759 + t701;
t814 = sin(qJ(5));
t820 = cos(qJ(5));
t694 = t814 * t697 + t820 * t699;
t750 = t820 * t770 - t814 * t771;
t751 = t814 * t770 + t820 * t771;
t730 = -t750 * pkin(5) - t751 * pkin(13);
t763 = qJDD(5) + t764;
t777 = qJD(5) + t779;
t776 = t777 ^ 2;
t692 = -t776 * pkin(5) + t763 * pkin(13) + t750 * t730 + t694;
t702 = -t733 * pkin(4) - t769 * pkin(12) + t771 * t759 + t707;
t714 = -t751 * qJD(5) + t820 * t733 - t814 * t734;
t715 = t750 * qJD(5) + t814 * t733 + t820 * t734;
t695 = (-t750 * t777 - t715) * pkin(13) + (t751 * t777 - t714) * pkin(5) + t702;
t813 = sin(qJ(6));
t819 = cos(qJ(6));
t689 = -t813 * t692 + t819 * t695;
t735 = -t813 * t751 + t819 * t777;
t705 = t735 * qJD(6) + t819 * t715 + t813 * t763;
t713 = qJDD(6) - t714;
t736 = t819 * t751 + t813 * t777;
t722 = -t735 * mrSges(7,1) + t736 * mrSges(7,2);
t747 = qJD(6) - t750;
t723 = -t747 * mrSges(7,2) + t735 * mrSges(7,3);
t687 = m(7) * t689 + t713 * mrSges(7,1) - t705 * mrSges(7,3) - t736 * t722 + t747 * t723;
t690 = t819 * t692 + t813 * t695;
t704 = -t736 * qJD(6) - t813 * t715 + t819 * t763;
t724 = t747 * mrSges(7,1) - t736 * mrSges(7,3);
t688 = m(7) * t690 - t713 * mrSges(7,2) + t704 * mrSges(7,3) + t735 * t722 - t747 * t724;
t679 = t819 * t687 + t813 * t688;
t737 = -t777 * mrSges(6,2) + t750 * mrSges(6,3);
t738 = t777 * mrSges(6,1) - t751 * mrSges(6,3);
t827 = m(6) * t702 - t714 * mrSges(6,1) + t715 * mrSges(6,2) - t750 * t737 + t751 * t738 + t679;
t826 = -m(5) * t707 + t733 * mrSges(5,1) - t734 * mrSges(5,2) + t770 * t757 - t771 * t758 - t827;
t675 = m(4) * t720 + t782 * mrSges(4,1) - t766 * mrSges(4,3) - t781 * t767 + t791 * t772 + t826;
t854 = t675 * t822;
t852 = t810 * t817;
t851 = t810 * t823;
t773 = t791 * mrSges(4,1) - t781 * mrSges(4,3);
t729 = -t750 * mrSges(6,1) + t751 * mrSges(6,2);
t835 = -t813 * t687 + t819 * t688;
t678 = m(6) * t694 - t763 * mrSges(6,2) + t714 * mrSges(6,3) + t750 * t729 - t777 * t738 + t835;
t693 = t820 * t697 - t814 * t699;
t691 = -t763 * pkin(5) - t776 * pkin(13) + t751 * t730 - t693;
t828 = -m(7) * t691 + t704 * mrSges(7,1) - t705 * mrSges(7,2) + t735 * t723 - t736 * t724;
t683 = m(6) * t693 + t763 * mrSges(6,1) - t715 * mrSges(6,3) - t751 * t729 + t777 * t737 + t828;
t672 = t814 * t678 + t820 * t683;
t752 = -t770 * mrSges(5,1) + t771 * mrSges(5,2);
t670 = m(5) * t700 + t764 * mrSges(5,1) - t734 * mrSges(5,3) - t771 * t752 + t779 * t757 + t672;
t836 = t820 * t678 - t814 * t683;
t671 = m(5) * t701 - t764 * mrSges(5,2) + t733 * mrSges(5,3) + t770 * t752 - t779 * t758 + t836;
t837 = -t815 * t670 + t821 * t671;
t661 = m(4) * t721 - t782 * mrSges(4,2) + t765 * mrSges(4,3) + t780 * t767 - t791 * t773 + t837;
t664 = t821 * t670 + t815 * t671;
t663 = m(4) * t731 - t765 * mrSges(4,1) + t766 * mrSges(4,2) - t780 * t772 + t781 * t773 + t664;
t650 = t661 * t850 - t809 * t663 + t811 * t854;
t774 = -g(3) * t851 + t834;
t797 = -t807 * mrSges(3,2) + mrSges(3,3) * t839;
t800 = (-mrSges(3,1) * t823 + mrSges(3,2) * t817) * t844;
t646 = m(3) * t774 + t806 * mrSges(3,1) - t801 * mrSges(3,3) + t807 * t797 - t800 * t840 + t650;
t649 = t661 * t853 + t811 * t663 + t809 * t854;
t786 = -t810 * t798 - t855;
t796 = t807 * mrSges(3,1) - mrSges(3,3) * t840;
t648 = m(3) * t786 - t802 * mrSges(3,1) + t801 * mrSges(3,2) + (t796 * t817 - t797 * t823) * t844 + t649;
t657 = t822 * t661 - t816 * t675;
t775 = -g(3) * t852 + t845;
t656 = m(3) * t775 - t806 * mrSges(3,2) + t802 * mrSges(3,3) - t807 * t796 + t800 * t839 + t657;
t636 = t646 * t848 - t810 * t648 + t656 * t849;
t634 = m(2) * t804 + qJDD(1) * mrSges(2,1) - t825 * mrSges(2,2) + t636;
t642 = -t817 * t646 + t823 * t656;
t641 = m(2) * t805 - t825 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t642;
t847 = t824 * t634 + t818 * t641;
t635 = t646 * t851 + t812 * t648 + t656 * t852;
t838 = -t818 * t634 + t824 * t641;
t716 = Ifges(7,5) * t736 + Ifges(7,6) * t735 + Ifges(7,3) * t747;
t718 = Ifges(7,1) * t736 + Ifges(7,4) * t735 + Ifges(7,5) * t747;
t680 = -mrSges(7,1) * t691 + mrSges(7,3) * t690 + Ifges(7,4) * t705 + Ifges(7,2) * t704 + Ifges(7,6) * t713 - t736 * t716 + t747 * t718;
t717 = Ifges(7,4) * t736 + Ifges(7,2) * t735 + Ifges(7,6) * t747;
t681 = mrSges(7,2) * t691 - mrSges(7,3) * t689 + Ifges(7,1) * t705 + Ifges(7,4) * t704 + Ifges(7,5) * t713 + t735 * t716 - t747 * t717;
t725 = Ifges(6,5) * t751 + Ifges(6,6) * t750 + Ifges(6,3) * t777;
t726 = Ifges(6,4) * t751 + Ifges(6,2) * t750 + Ifges(6,6) * t777;
t665 = mrSges(6,2) * t702 - mrSges(6,3) * t693 + Ifges(6,1) * t715 + Ifges(6,4) * t714 + Ifges(6,5) * t763 - pkin(13) * t679 - t813 * t680 + t819 * t681 + t750 * t725 - t777 * t726;
t727 = Ifges(6,1) * t751 + Ifges(6,4) * t750 + Ifges(6,5) * t777;
t666 = -mrSges(6,1) * t702 - mrSges(7,1) * t689 + mrSges(7,2) * t690 + mrSges(6,3) * t694 + Ifges(6,4) * t715 - Ifges(7,5) * t705 + Ifges(6,2) * t714 + Ifges(6,6) * t763 - Ifges(7,6) * t704 - Ifges(7,3) * t713 - pkin(5) * t679 - t736 * t717 + t735 * t718 - t751 * t725 + t777 * t727;
t739 = Ifges(5,5) * t771 + Ifges(5,6) * t770 + Ifges(5,3) * t779;
t741 = Ifges(5,1) * t771 + Ifges(5,4) * t770 + Ifges(5,5) * t779;
t651 = -mrSges(5,1) * t707 + mrSges(5,3) * t701 + Ifges(5,4) * t734 + Ifges(5,2) * t733 + Ifges(5,6) * t764 - pkin(4) * t827 + pkin(12) * t836 + t814 * t665 + t820 * t666 - t771 * t739 + t779 * t741;
t740 = Ifges(5,4) * t771 + Ifges(5,2) * t770 + Ifges(5,6) * t779;
t652 = mrSges(5,2) * t707 - mrSges(5,3) * t700 + Ifges(5,1) * t734 + Ifges(5,4) * t733 + Ifges(5,5) * t764 - pkin(12) * t672 + t820 * t665 - t814 * t666 + t770 * t739 - t779 * t740;
t761 = Ifges(4,4) * t781 + Ifges(4,2) * t780 + Ifges(4,6) * t791;
t762 = Ifges(4,1) * t781 + Ifges(4,4) * t780 + Ifges(4,5) * t791;
t637 = mrSges(4,1) * t720 - mrSges(4,2) * t721 + Ifges(4,5) * t766 + Ifges(4,6) * t765 + Ifges(4,3) * t782 + pkin(3) * t826 + pkin(11) * t837 + t821 * t651 + t815 * t652 + t781 * t761 - t780 * t762;
t783 = Ifges(3,3) * t807 + (Ifges(3,5) * t817 + Ifges(3,6) * t823) * t844;
t785 = Ifges(3,5) * t807 + (Ifges(3,1) * t817 + Ifges(3,4) * t823) * t844;
t760 = Ifges(4,5) * t781 + Ifges(4,6) * t780 + Ifges(4,3) * t791;
t638 = mrSges(4,2) * t731 - mrSges(4,3) * t720 + Ifges(4,1) * t766 + Ifges(4,4) * t765 + Ifges(4,5) * t782 - pkin(11) * t664 - t815 * t651 + t821 * t652 + t780 * t760 - t791 * t761;
t643 = -pkin(5) * t828 - pkin(13) * t835 - pkin(3) * t664 - t819 * t680 - t813 * t681 + t791 * t762 - t781 * t760 + Ifges(4,6) * t782 + t770 * t741 - t771 * t740 - Ifges(6,3) * t763 - Ifges(5,3) * t764 + Ifges(4,2) * t765 + Ifges(4,4) * t766 - t751 * t726 + t750 * t727 - mrSges(4,1) * t731 - Ifges(5,6) * t733 - Ifges(5,5) * t734 + mrSges(4,3) * t721 - Ifges(6,5) * t715 - Ifges(6,6) * t714 - mrSges(5,1) * t700 + mrSges(5,2) * t701 - mrSges(6,1) * t693 + mrSges(6,2) * t694 - pkin(4) * t672;
t829 = pkin(10) * t657 + t638 * t816 + t643 * t822;
t631 = -mrSges(3,1) * t786 + mrSges(3,3) * t775 + Ifges(3,4) * t801 + Ifges(3,2) * t802 + Ifges(3,6) * t806 - pkin(2) * t649 - t809 * t637 - t783 * t840 + t807 * t785 + t811 * t829;
t784 = Ifges(3,6) * t807 + (Ifges(3,4) * t817 + Ifges(3,2) * t823) * t844;
t632 = t783 * t839 + mrSges(3,2) * t786 - mrSges(3,3) * t774 + Ifges(3,1) * t801 + Ifges(3,4) * t802 + Ifges(3,5) * t806 + t822 * t638 - t816 * t643 - t807 * t784 + (-t649 * t809 - t650 * t811) * pkin(10);
t830 = pkin(9) * t642 + t631 * t823 + t632 * t817;
t630 = mrSges(3,1) * t774 - mrSges(3,2) * t775 + Ifges(3,5) * t801 + Ifges(3,6) * t802 + Ifges(3,3) * t806 + pkin(2) * t650 + t811 * t637 + (t784 * t817 - t785 * t823) * t844 + t829 * t809;
t629 = -mrSges(2,2) * g(3) - mrSges(2,3) * t804 + Ifges(2,5) * qJDD(1) - t825 * Ifges(2,6) - t817 * t631 + t823 * t632 + (-t635 * t810 - t636 * t812) * pkin(9);
t628 = mrSges(2,1) * g(3) + mrSges(2,3) * t805 + t825 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t635 - t810 * t630 + t812 * t830;
t1 = [-m(1) * g(1) + t838; -m(1) * g(2) + t847; (-m(1) - m(2)) * g(3) + t635; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t847 - t818 * t628 + t824 * t629; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t838 + t824 * t628 + t818 * t629; -mrSges(1,1) * g(2) + mrSges(2,1) * t804 + mrSges(1,2) * g(1) - mrSges(2,2) * t805 + Ifges(2,3) * qJDD(1) + pkin(1) * t636 + t812 * t630 + t810 * t830;];
tauB  = t1;
