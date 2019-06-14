% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:05:08
% EndTime: 2019-05-05 22:05:26
% DurationCPUTime: 17.99s
% Computational Cost: add. (306602->343), mult. (617599->429), div. (0->0), fcn. (419852->12), ass. (0->138)
t837 = sin(qJ(1));
t841 = cos(qJ(1));
t819 = t837 * g(1) - g(2) * t841;
t810 = qJDD(1) * pkin(1) + t819;
t820 = -g(1) * t841 - g(2) * t837;
t843 = qJD(1) ^ 2;
t812 = -pkin(1) * t843 + t820;
t831 = sin(pkin(10));
t833 = cos(pkin(10));
t789 = t833 * t810 - t831 * t812;
t773 = -qJDD(1) * pkin(2) - t843 * pkin(7) - t789;
t836 = sin(qJ(3));
t840 = cos(qJ(3));
t857 = qJD(1) * qJD(3);
t856 = t840 * t857;
t814 = qJDD(1) * t836 + t856;
t824 = t836 * t857;
t815 = qJDD(1) * t840 - t824;
t756 = (-t814 - t856) * pkin(8) + (-t815 + t824) * pkin(3) + t773;
t790 = t831 * t810 + t833 * t812;
t774 = -pkin(2) * t843 + qJDD(1) * pkin(7) + t790;
t829 = -g(3) + qJDD(2);
t767 = t840 * t774 + t836 * t829;
t813 = (-pkin(3) * t840 - pkin(8) * t836) * qJD(1);
t842 = qJD(3) ^ 2;
t858 = qJD(1) * t840;
t764 = -pkin(3) * t842 + qJDD(3) * pkin(8) + t813 * t858 + t767;
t835 = sin(qJ(4));
t839 = cos(qJ(4));
t745 = t839 * t756 - t835 * t764;
t859 = qJD(1) * t836;
t808 = qJD(3) * t839 - t835 * t859;
t784 = qJD(4) * t808 + qJDD(3) * t835 + t814 * t839;
t807 = qJDD(4) - t815;
t809 = qJD(3) * t835 + t839 * t859;
t822 = qJD(4) - t858;
t729 = (t808 * t822 - t784) * qJ(5) + (t808 * t809 + t807) * pkin(4) + t745;
t746 = t835 * t756 + t839 * t764;
t783 = -qJD(4) * t809 + qJDD(3) * t839 - t814 * t835;
t793 = pkin(4) * t822 - qJ(5) * t809;
t806 = t808 ^ 2;
t737 = -pkin(4) * t806 + qJ(5) * t783 - t793 * t822 + t746;
t830 = sin(pkin(11));
t832 = cos(pkin(11));
t788 = t808 * t830 + t809 * t832;
t723 = -0.2e1 * qJD(5) * t788 + t832 * t729 - t830 * t737;
t758 = t783 * t830 + t784 * t832;
t787 = t808 * t832 - t809 * t830;
t721 = (t787 * t822 - t758) * pkin(9) + (t787 * t788 + t807) * pkin(5) + t723;
t724 = 0.2e1 * qJD(5) * t787 + t830 * t729 + t832 * t737;
t757 = t783 * t832 - t784 * t830;
t770 = pkin(5) * t822 - pkin(9) * t788;
t785 = t787 ^ 2;
t722 = -pkin(5) * t785 + pkin(9) * t757 - t770 * t822 + t724;
t834 = sin(qJ(6));
t838 = cos(qJ(6));
t719 = t721 * t838 - t722 * t834;
t761 = t787 * t838 - t788 * t834;
t736 = qJD(6) * t761 + t757 * t834 + t758 * t838;
t762 = t787 * t834 + t788 * t838;
t747 = -mrSges(7,1) * t761 + mrSges(7,2) * t762;
t821 = qJD(6) + t822;
t748 = -mrSges(7,2) * t821 + mrSges(7,3) * t761;
t803 = qJDD(6) + t807;
t712 = m(7) * t719 + mrSges(7,1) * t803 - mrSges(7,3) * t736 - t747 * t762 + t748 * t821;
t720 = t721 * t834 + t722 * t838;
t735 = -qJD(6) * t762 + t757 * t838 - t758 * t834;
t749 = mrSges(7,1) * t821 - mrSges(7,3) * t762;
t713 = m(7) * t720 - mrSges(7,2) * t803 + mrSges(7,3) * t735 + t747 * t761 - t749 * t821;
t706 = t838 * t712 + t834 * t713;
t765 = -mrSges(6,1) * t787 + mrSges(6,2) * t788;
t768 = -mrSges(6,2) * t822 + mrSges(6,3) * t787;
t704 = m(6) * t723 + mrSges(6,1) * t807 - mrSges(6,3) * t758 - t765 * t788 + t768 * t822 + t706;
t769 = mrSges(6,1) * t822 - mrSges(6,3) * t788;
t851 = -t712 * t834 + t838 * t713;
t705 = m(6) * t724 - mrSges(6,2) * t807 + mrSges(6,3) * t757 + t765 * t787 - t769 * t822 + t851;
t700 = t832 * t704 + t830 * t705;
t754 = Ifges(6,4) * t788 + Ifges(6,2) * t787 + Ifges(6,6) * t822;
t755 = Ifges(6,1) * t788 + Ifges(6,4) * t787 + Ifges(6,5) * t822;
t776 = Ifges(5,4) * t809 + Ifges(5,2) * t808 + Ifges(5,6) * t822;
t777 = Ifges(5,1) * t809 + Ifges(5,4) * t808 + Ifges(5,5) * t822;
t740 = Ifges(7,4) * t762 + Ifges(7,2) * t761 + Ifges(7,6) * t821;
t741 = Ifges(7,1) * t762 + Ifges(7,4) * t761 + Ifges(7,5) * t821;
t848 = -mrSges(7,1) * t719 + mrSges(7,2) * t720 - Ifges(7,5) * t736 - Ifges(7,6) * t735 - Ifges(7,3) * t803 - t762 * t740 + t761 * t741;
t863 = mrSges(5,1) * t745 + mrSges(6,1) * t723 - mrSges(5,2) * t746 - mrSges(6,2) * t724 + Ifges(5,5) * t784 + Ifges(6,5) * t758 + Ifges(5,6) * t783 + Ifges(6,6) * t757 + pkin(4) * t700 + pkin(5) * t706 + t788 * t754 - t787 * t755 + t809 * t776 - t808 * t777 + (Ifges(5,3) + Ifges(6,3)) * t807 - t848;
t766 = -t836 * t774 + t829 * t840;
t763 = -qJDD(3) * pkin(3) - pkin(8) * t842 + t813 * t859 - t766;
t742 = -pkin(4) * t783 - qJ(5) * t806 + t809 * t793 + qJDD(5) + t763;
t726 = -pkin(5) * t757 - pkin(9) * t785 + t770 * t788 + t742;
t739 = Ifges(7,5) * t762 + Ifges(7,6) * t761 + Ifges(7,3) * t821;
t707 = -mrSges(7,1) * t726 + mrSges(7,3) * t720 + Ifges(7,4) * t736 + Ifges(7,2) * t735 + Ifges(7,6) * t803 - t739 * t762 + t741 * t821;
t708 = mrSges(7,2) * t726 - mrSges(7,3) * t719 + Ifges(7,1) * t736 + Ifges(7,4) * t735 + Ifges(7,5) * t803 + t739 * t761 - t740 * t821;
t753 = Ifges(6,5) * t788 + Ifges(6,6) * t787 + Ifges(6,3) * t822;
t850 = m(7) * t726 - t735 * mrSges(7,1) + t736 * mrSges(7,2) - t761 * t748 + t762 * t749;
t695 = -mrSges(6,1) * t742 + mrSges(6,3) * t724 + Ifges(6,4) * t758 + Ifges(6,2) * t757 + Ifges(6,6) * t807 - pkin(5) * t850 + pkin(9) * t851 + t838 * t707 + t834 * t708 - t788 * t753 + t822 * t755;
t696 = mrSges(6,2) * t742 - mrSges(6,3) * t723 + Ifges(6,1) * t758 + Ifges(6,4) * t757 + Ifges(6,5) * t807 - pkin(9) * t706 - t707 * t834 + t708 * t838 + t753 * t787 - t754 * t822;
t717 = m(6) * t742 - t757 * mrSges(6,1) + t758 * mrSges(6,2) - t787 * t768 + t788 * t769 + t850;
t775 = Ifges(5,5) * t809 + Ifges(5,6) * t808 + Ifges(5,3) * t822;
t852 = -t704 * t830 + t832 * t705;
t679 = -mrSges(5,1) * t763 + mrSges(5,3) * t746 + Ifges(5,4) * t784 + Ifges(5,2) * t783 + Ifges(5,6) * t807 - pkin(4) * t717 + qJ(5) * t852 + t832 * t695 + t830 * t696 - t809 * t775 + t822 * t777;
t680 = mrSges(5,2) * t763 - mrSges(5,3) * t745 + Ifges(5,1) * t784 + Ifges(5,4) * t783 + Ifges(5,5) * t807 - qJ(5) * t700 - t695 * t830 + t696 * t832 + t775 * t808 - t776 * t822;
t791 = -mrSges(5,1) * t808 + mrSges(5,2) * t809;
t792 = -mrSges(5,2) * t822 + mrSges(5,3) * t808;
t698 = m(5) * t745 + mrSges(5,1) * t807 - mrSges(5,3) * t784 - t791 * t809 + t792 * t822 + t700;
t794 = mrSges(5,1) * t822 - mrSges(5,3) * t809;
t699 = m(5) * t746 - mrSges(5,2) * t807 + mrSges(5,3) * t783 + t791 * t808 - t794 * t822 + t852;
t694 = -t698 * t835 + t839 * t699;
t716 = -m(5) * t763 + t783 * mrSges(5,1) - t784 * mrSges(5,2) + t808 * t792 - t809 * t794 - t717;
t801 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t836 + Ifges(4,2) * t840) * qJD(1);
t802 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t836 + Ifges(4,4) * t840) * qJD(1);
t862 = mrSges(4,1) * t766 - mrSges(4,2) * t767 + Ifges(4,5) * t814 + Ifges(4,6) * t815 + Ifges(4,3) * qJDD(3) + pkin(3) * t716 + pkin(8) * t694 + t839 * t679 + t835 * t680 + (t801 * t836 - t802 * t840) * qJD(1);
t811 = (-mrSges(4,1) * t840 + mrSges(4,2) * t836) * qJD(1);
t817 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t859;
t692 = m(4) * t767 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t815 - qJD(3) * t817 + t811 * t858 + t694;
t818 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t858;
t715 = m(4) * t766 + qJDD(3) * mrSges(4,1) - t814 * mrSges(4,3) + qJD(3) * t818 - t811 * t859 + t716;
t853 = t840 * t692 - t715 * t836;
t683 = m(3) * t790 - mrSges(3,1) * t843 - qJDD(1) * mrSges(3,2) + t853;
t693 = t698 * t839 + t699 * t835;
t846 = -m(4) * t773 + t815 * mrSges(4,1) - mrSges(4,2) * t814 - t817 * t859 + t818 * t858 - t693;
t688 = m(3) * t789 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t843 + t846;
t676 = t831 * t683 + t833 * t688;
t673 = m(2) * t819 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t843 + t676;
t854 = t833 * t683 - t688 * t831;
t674 = m(2) * t820 - mrSges(2,1) * t843 - qJDD(1) * mrSges(2,2) + t854;
t860 = t841 * t673 + t837 * t674;
t686 = t836 * t692 + t840 * t715;
t684 = m(3) * t829 + t686;
t855 = -t673 * t837 + t841 * t674;
t800 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t836 + Ifges(4,6) * t840) * qJD(1);
t669 = mrSges(4,2) * t773 - mrSges(4,3) * t766 + Ifges(4,1) * t814 + Ifges(4,4) * t815 + Ifges(4,5) * qJDD(3) - pkin(8) * t693 - qJD(3) * t801 - t679 * t835 + t680 * t839 + t800 * t858;
t678 = -mrSges(4,1) * t773 + mrSges(4,3) * t767 + Ifges(4,4) * t814 + Ifges(4,2) * t815 + Ifges(4,6) * qJDD(3) - pkin(3) * t693 + qJD(3) * t802 - t800 * t859 - t863;
t847 = mrSges(2,1) * t819 + mrSges(3,1) * t789 - mrSges(2,2) * t820 - mrSges(3,2) * t790 + pkin(1) * t676 + pkin(2) * t846 + pkin(7) * t853 + t836 * t669 + t840 * t678 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t667 = -mrSges(3,1) * t829 + mrSges(3,3) * t790 + t843 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t686 - t862;
t666 = mrSges(3,2) * t829 - mrSges(3,3) * t789 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t843 - pkin(7) * t686 + t669 * t840 - t678 * t836;
t665 = -mrSges(2,2) * g(3) - mrSges(2,3) * t819 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t843 - qJ(2) * t676 + t666 * t833 - t667 * t831;
t664 = mrSges(2,1) * g(3) + mrSges(2,3) * t820 + t843 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t684 + qJ(2) * t854 + t831 * t666 + t833 * t667;
t1 = [-m(1) * g(1) + t855; -m(1) * g(2) + t860; (-m(1) - m(2)) * g(3) + t684; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t860 - t837 * t664 + t841 * t665; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t855 + t841 * t664 + t837 * t665; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t847; t847; t684; t862; t863; t717; -t848;];
tauJB  = t1;
