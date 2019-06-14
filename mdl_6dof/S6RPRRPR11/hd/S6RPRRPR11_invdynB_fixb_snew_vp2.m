% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-06 00:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:09:29
% EndTime: 2019-05-06 00:11:18
% DurationCPUTime: 112.25s
% Computational Cost: add. (1744086->402), mult. (5459090->551), div. (0->0), fcn. (4688557->16), ass. (0->182)
t800 = sin(pkin(12));
t802 = sin(pkin(6));
t804 = cos(pkin(12));
t806 = cos(pkin(6));
t809 = sin(qJ(3));
t805 = cos(pkin(7));
t812 = cos(qJ(3));
t846 = t805 * t812;
t801 = sin(pkin(7));
t851 = t801 * t812;
t817 = t802 * (-t800 * t809 + t804 * t846) + t806 * t851;
t772 = t817 * qJD(1);
t847 = t805 * t809;
t852 = t801 * t809;
t819 = t806 * t852 + (t800 * t812 + t804 * t847) * t802;
t773 = t819 * qJD(1);
t761 = -t773 * qJD(3) + qJDD(1) * t817;
t849 = t802 * t805;
t786 = (t801 * t806 + t804 * t849) * qJD(1) * pkin(9);
t810 = sin(qJ(1));
t813 = cos(qJ(1));
t797 = -t813 * g(1) - t810 * g(2);
t814 = qJD(1) ^ 2;
t855 = qJ(2) * t802;
t790 = -t814 * pkin(1) + qJDD(1) * t855 + t797;
t858 = pkin(9) * t800;
t830 = -pkin(2) * t804 - t801 * t858;
t844 = qJD(1) * t802;
t856 = pkin(9) * qJDD(1);
t825 = qJD(1) * t830 * t844 + t805 * t856;
t796 = t810 * g(1) - t813 * g(2);
t789 = qJDD(1) * pkin(1) + t814 * t855 + t796;
t840 = qJD(2) * t844;
t848 = t804 * t806;
t850 = t802 * t804;
t831 = -g(3) * t850 + t789 * t848 - 0.2e1 * t800 * t840;
t740 = (pkin(2) * qJDD(1) + qJD(1) * t786) * t806 + (-t802 * t825 - t790) * t800 + t831;
t791 = (pkin(2) * t806 - t849 * t858) * qJD(1);
t853 = t800 * t806;
t841 = t789 * t853 + (t790 + 0.2e1 * t840) * t804;
t741 = (-qJD(1) * t791 + t801 * t856) * t806 + (-g(3) * t800 + t804 * t825) * t802 + t841;
t839 = -t806 * g(3) + qJDD(2);
t752 = (-t789 + t830 * qJDD(1) + (-t786 * t804 + t791 * t800) * qJD(1)) * t802 + t839;
t710 = -t809 * t741 + (t740 * t805 + t752 * t801) * t812;
t859 = cos(qJ(4));
t857 = Ifges(3,3) * t806;
t854 = t800 * t802;
t711 = t740 * t847 + t812 * t741 + t752 * t852;
t759 = -t772 * mrSges(4,1) + t773 * mrSges(4,2);
t826 = -t801 * t850 + t805 * t806;
t787 = qJD(1) * t826 + qJD(3);
t769 = t787 * mrSges(4,1) - t773 * mrSges(4,3);
t784 = qJDD(1) * t826 + qJDD(3);
t760 = -t772 * pkin(3) - t773 * pkin(10);
t783 = t787 ^ 2;
t702 = -t783 * pkin(3) + t784 * pkin(10) + t772 * t760 + t711;
t719 = -t801 * t740 + t805 * t752;
t762 = t772 * qJD(3) + qJDD(1) * t819;
t708 = (-t772 * t787 - t762) * pkin(10) + (t773 * t787 - t761) * pkin(3) + t719;
t808 = sin(qJ(4));
t695 = t702 * t859 + t808 * t708;
t767 = t773 * t859 + t808 * t787;
t735 = t767 * qJD(4) + t808 * t762 - t784 * t859;
t766 = t808 * t773 - t787 * t859;
t743 = t766 * mrSges(5,1) + t767 * mrSges(5,2);
t771 = qJD(4) - t772;
t754 = t771 * mrSges(5,1) - t767 * mrSges(5,3);
t758 = qJDD(4) - t761;
t742 = t766 * pkin(4) - t767 * qJ(5);
t770 = t771 ^ 2;
t690 = -t770 * pkin(4) + t758 * qJ(5) - t766 * t742 + t695;
t701 = -t784 * pkin(3) - t783 * pkin(10) + t773 * t760 - t710;
t736 = -t766 * qJD(4) + t762 * t859 + t808 * t784;
t693 = (t766 * t771 - t736) * qJ(5) + (t767 * t771 + t735) * pkin(4) + t701;
t799 = sin(pkin(13));
t803 = cos(pkin(13));
t751 = t803 * t767 + t799 * t771;
t685 = -0.2e1 * qJD(5) * t751 - t799 * t690 + t803 * t693;
t721 = t803 * t736 + t799 * t758;
t750 = -t799 * t767 + t803 * t771;
t683 = (t750 * t766 - t721) * pkin(11) + (t750 * t751 + t735) * pkin(5) + t685;
t686 = 0.2e1 * qJD(5) * t750 + t803 * t690 + t799 * t693;
t720 = -t799 * t736 + t803 * t758;
t728 = t766 * pkin(5) - t751 * pkin(11);
t749 = t750 ^ 2;
t684 = -t749 * pkin(5) + t720 * pkin(11) - t766 * t728 + t686;
t807 = sin(qJ(6));
t811 = cos(qJ(6));
t681 = t811 * t683 - t807 * t684;
t722 = t811 * t750 - t807 * t751;
t698 = t722 * qJD(6) + t807 * t720 + t811 * t721;
t723 = t807 * t750 + t811 * t751;
t709 = -t722 * mrSges(7,1) + t723 * mrSges(7,2);
t763 = qJD(6) + t766;
t712 = -t763 * mrSges(7,2) + t722 * mrSges(7,3);
t733 = qJDD(6) + t735;
t679 = m(7) * t681 + t733 * mrSges(7,1) - t698 * mrSges(7,3) - t723 * t709 + t763 * t712;
t682 = t807 * t683 + t811 * t684;
t697 = -t723 * qJD(6) + t811 * t720 - t807 * t721;
t713 = t763 * mrSges(7,1) - t723 * mrSges(7,3);
t680 = m(7) * t682 - t733 * mrSges(7,2) + t697 * mrSges(7,3) + t722 * t709 - t763 * t713;
t671 = t811 * t679 + t807 * t680;
t724 = -t750 * mrSges(6,1) + t751 * mrSges(6,2);
t726 = -t766 * mrSges(6,2) + t750 * mrSges(6,3);
t669 = m(6) * t685 + t735 * mrSges(6,1) - t721 * mrSges(6,3) - t751 * t724 + t766 * t726 + t671;
t727 = t766 * mrSges(6,1) - t751 * mrSges(6,3);
t835 = -t807 * t679 + t811 * t680;
t670 = m(6) * t686 - t735 * mrSges(6,2) + t720 * mrSges(6,3) + t750 * t724 - t766 * t727 + t835;
t836 = -t799 * t669 + t803 * t670;
t666 = m(5) * t695 - t758 * mrSges(5,2) - t735 * mrSges(5,3) - t766 * t743 - t771 * t754 + t836;
t694 = -t808 * t702 + t708 * t859;
t753 = -t771 * mrSges(5,2) - t766 * mrSges(5,3);
t689 = -t758 * pkin(4) - t770 * qJ(5) + t767 * t742 + qJDD(5) - t694;
t687 = -t720 * pkin(5) - t749 * pkin(11) + t751 * t728 + t689;
t822 = m(7) * t687 - t697 * mrSges(7,1) + t698 * mrSges(7,2) - t722 * t712 + t723 * t713;
t815 = -m(6) * t689 + t720 * mrSges(6,1) - t721 * mrSges(6,2) + t750 * t726 - t751 * t727 - t822;
t675 = m(5) * t694 + t758 * mrSges(5,1) - t736 * mrSges(5,3) - t767 * t743 + t771 * t753 + t815;
t837 = t666 * t859 - t808 * t675;
t655 = m(4) * t711 - t784 * mrSges(4,2) + t761 * mrSges(4,3) + t772 * t759 - t787 * t769 + t837;
t658 = t808 * t666 + t675 * t859;
t768 = -t787 * mrSges(4,2) + t772 * mrSges(4,3);
t657 = m(4) * t719 - t761 * mrSges(4,1) + t762 * mrSges(4,2) - t772 * t768 + t773 * t769 + t658;
t667 = t803 * t669 + t799 * t670;
t816 = -m(5) * t701 - t735 * mrSges(5,1) - t736 * mrSges(5,2) - t766 * t753 - t767 * t754 - t667;
t663 = m(4) * t710 + t784 * mrSges(4,1) - t762 * mrSges(4,3) - t773 * t759 + t787 * t768 + t816;
t644 = t655 * t847 - t801 * t657 + t663 * t846;
t764 = -t800 * t790 + t831;
t834 = -mrSges(3,1) * t804 + mrSges(3,2) * t800;
t788 = t834 * t844;
t828 = -mrSges(3,2) * t806 + mrSges(3,3) * t850;
t793 = t828 * qJD(1);
t829 = mrSges(3,1) * t806 - mrSges(3,3) * t854;
t640 = m(3) * t764 + t829 * qJDD(1) + (-t788 * t854 + t793 * t806) * qJD(1) + t644;
t643 = t655 * t852 + t805 * t657 + t663 * t851;
t774 = -t802 * t789 + t839;
t792 = t829 * qJD(1);
t642 = m(3) * t774 + (t834 * qJDD(1) + (t792 * t800 - t793 * t804) * qJD(1)) * t802 + t643;
t650 = t812 * t655 - t809 * t663;
t765 = -g(3) * t854 + t841;
t649 = m(3) * t765 + t828 * qJDD(1) + (t788 * t850 - t792 * t806) * qJD(1) + t650;
t630 = t640 * t848 - t802 * t642 + t649 * t853;
t628 = m(2) * t796 + qJDD(1) * mrSges(2,1) - t814 * mrSges(2,2) + t630;
t636 = -t800 * t640 + t804 * t649;
t635 = m(2) * t797 - t814 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t636;
t845 = t813 * t628 + t810 * t635;
t629 = t640 * t850 + t806 * t642 + t649 * t854;
t838 = -t810 * t628 + t813 * t635;
t833 = Ifges(3,5) * t800 + Ifges(3,6) * t804;
t704 = Ifges(7,5) * t723 + Ifges(7,6) * t722 + Ifges(7,3) * t763;
t706 = Ifges(7,1) * t723 + Ifges(7,4) * t722 + Ifges(7,5) * t763;
t672 = -mrSges(7,1) * t687 + mrSges(7,3) * t682 + Ifges(7,4) * t698 + Ifges(7,2) * t697 + Ifges(7,6) * t733 - t723 * t704 + t763 * t706;
t705 = Ifges(7,4) * t723 + Ifges(7,2) * t722 + Ifges(7,6) * t763;
t673 = mrSges(7,2) * t687 - mrSges(7,3) * t681 + Ifges(7,1) * t698 + Ifges(7,4) * t697 + Ifges(7,5) * t733 + t722 * t704 - t763 * t705;
t714 = Ifges(6,5) * t751 + Ifges(6,6) * t750 + Ifges(6,3) * t766;
t716 = Ifges(6,1) * t751 + Ifges(6,4) * t750 + Ifges(6,5) * t766;
t659 = -mrSges(6,1) * t689 + mrSges(6,3) * t686 + Ifges(6,4) * t721 + Ifges(6,2) * t720 + Ifges(6,6) * t735 - pkin(5) * t822 + pkin(11) * t835 + t811 * t672 + t807 * t673 - t751 * t714 + t766 * t716;
t715 = Ifges(6,4) * t751 + Ifges(6,2) * t750 + Ifges(6,6) * t766;
t660 = mrSges(6,2) * t689 - mrSges(6,3) * t685 + Ifges(6,1) * t721 + Ifges(6,4) * t720 + Ifges(6,5) * t735 - pkin(11) * t671 - t807 * t672 + t811 * t673 + t750 * t714 - t766 * t715;
t729 = Ifges(5,5) * t767 - Ifges(5,6) * t766 + Ifges(5,3) * t771;
t730 = Ifges(5,4) * t767 - Ifges(5,2) * t766 + Ifges(5,6) * t771;
t645 = mrSges(5,2) * t701 - mrSges(5,3) * t694 + Ifges(5,1) * t736 - Ifges(5,4) * t735 + Ifges(5,5) * t758 - qJ(5) * t667 - t799 * t659 + t803 * t660 - t766 * t729 - t771 * t730;
t731 = Ifges(5,1) * t767 - Ifges(5,4) * t766 + Ifges(5,5) * t771;
t651 = Ifges(5,4) * t736 + Ifges(5,6) * t758 - t767 * t729 + t771 * t731 - mrSges(5,1) * t701 + mrSges(5,3) * t695 - Ifges(6,5) * t721 - Ifges(6,6) * t720 - t751 * t715 + t750 * t716 - mrSges(6,1) * t685 + mrSges(6,2) * t686 - Ifges(7,5) * t698 - Ifges(7,6) * t697 - Ifges(7,3) * t733 - t723 * t705 + t722 * t706 - mrSges(7,1) * t681 + mrSges(7,2) * t682 - pkin(5) * t671 - pkin(4) * t667 + (-Ifges(5,2) - Ifges(6,3)) * t735;
t755 = Ifges(4,5) * t773 + Ifges(4,6) * t772 + Ifges(4,3) * t787;
t756 = Ifges(4,4) * t773 + Ifges(4,2) * t772 + Ifges(4,6) * t787;
t632 = mrSges(4,2) * t719 - mrSges(4,3) * t710 + Ifges(4,1) * t762 + Ifges(4,4) * t761 + Ifges(4,5) * t784 - pkin(10) * t658 + t645 * t859 - t808 * t651 + t772 * t755 - t787 * t756;
t757 = Ifges(4,1) * t773 + Ifges(4,4) * t772 + Ifges(4,5) * t787;
t637 = Ifges(4,4) * t762 + Ifges(4,2) * t761 + Ifges(4,6) * t784 - t773 * t755 + t787 * t757 - mrSges(4,1) * t719 + mrSges(4,3) * t711 - Ifges(5,5) * t736 + Ifges(5,6) * t735 - Ifges(5,3) * t758 - t767 * t730 - t766 * t731 - mrSges(5,1) * t694 + mrSges(5,2) * t695 - t799 * t660 - t803 * t659 - pkin(4) * t815 - qJ(5) * t836 - pkin(3) * t658;
t824 = pkin(9) * t650 + t632 * t809 + t637 * t812;
t631 = mrSges(4,1) * t710 - mrSges(4,2) * t711 + Ifges(4,5) * t762 + Ifges(4,6) * t761 + Ifges(4,3) * t784 + pkin(3) * t816 + pkin(10) * t837 + t808 * t645 + t651 * t859 + t773 * t756 - t772 * t757;
t777 = (t802 * t833 + t857) * qJD(1);
t821 = Ifges(3,5) * t806 + (Ifges(3,1) * t800 + Ifges(3,4) * t804) * t802;
t779 = t821 * qJD(1);
t820 = Ifges(3,6) * t806 + (Ifges(3,4) * t800 + Ifges(3,2) * t804) * t802;
t625 = -mrSges(3,1) * t774 + mrSges(3,3) * t765 - pkin(2) * t643 - t801 * t631 + (-t777 * t854 + t779 * t806) * qJD(1) + t824 * t805 + t820 * qJDD(1);
t778 = t820 * qJD(1);
t626 = mrSges(3,2) * t774 - mrSges(3,3) * t764 + t812 * t632 - t809 * t637 + (t777 * t850 - t778 * t806) * qJD(1) + (-t643 * t801 - t644 * t805) * pkin(9) + t821 * qJDD(1);
t823 = qJ(2) * t636 + t625 * t804 + t626 * t800;
t624 = qJDD(1) * t857 + mrSges(3,1) * t764 - mrSges(3,2) * t765 + pkin(2) * t644 + t805 * t631 + t824 * t801 + (t833 * qJDD(1) + (t778 * t800 - t779 * t804) * qJD(1)) * t802;
t623 = -mrSges(2,2) * g(3) - mrSges(2,3) * t796 + Ifges(2,5) * qJDD(1) - t814 * Ifges(2,6) - t800 * t625 + t804 * t626 + (-t629 * t802 - t630 * t806) * qJ(2);
t622 = mrSges(2,1) * g(3) + mrSges(2,3) * t797 + t814 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t629 - t802 * t624 + t806 * t823;
t1 = [-m(1) * g(1) + t838; -m(1) * g(2) + t845; (-m(1) - m(2)) * g(3) + t629; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t845 - t810 * t622 + t813 * t623; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t838 + t813 * t622 + t810 * t623; -mrSges(1,1) * g(2) + mrSges(2,1) * t796 + mrSges(1,2) * g(1) - mrSges(2,2) * t797 + Ifges(2,3) * qJDD(1) + pkin(1) * t630 + t806 * t624 + t802 * t823;];
tauB  = t1;
