% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR1
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
% Datum: 2019-05-05 21:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:56:13
% EndTime: 2019-05-05 21:56:28
% DurationCPUTime: 15.77s
% Computational Cost: add. (255569->345), mult. (547314->434), div. (0->0), fcn. (381767->12), ass. (0->140)
t869 = sin(qJ(1));
t873 = cos(qJ(1));
t847 = t869 * g(1) - g(2) * t873;
t838 = qJDD(1) * pkin(1) + t847;
t848 = -g(1) * t873 - g(2) * t869;
t874 = qJD(1) ^ 2;
t840 = -pkin(1) * t874 + t848;
t863 = sin(pkin(10));
t865 = cos(pkin(10));
t822 = t863 * t838 + t865 * t840;
t818 = -pkin(2) * t874 + qJDD(1) * pkin(7) + t822;
t861 = -g(3) + qJDD(2);
t868 = sin(qJ(3));
t872 = cos(qJ(3));
t801 = -t868 * t818 + t872 * t861;
t892 = qJD(1) * qJD(3);
t891 = t872 * t892;
t841 = qJDD(1) * t868 + t891;
t792 = (-t841 + t891) * pkin(8) + (t868 * t872 * t874 + qJDD(3)) * pkin(3) + t801;
t802 = t872 * t818 + t868 * t861;
t842 = qJDD(1) * t872 - t868 * t892;
t894 = qJD(1) * t868;
t846 = qJD(3) * pkin(3) - pkin(8) * t894;
t860 = t872 ^ 2;
t793 = -pkin(3) * t860 * t874 + pkin(8) * t842 - qJD(3) * t846 + t802;
t867 = sin(qJ(4));
t871 = cos(qJ(4));
t765 = t871 * t792 - t867 * t793;
t833 = (-t867 * t868 + t871 * t872) * qJD(1);
t804 = qJD(4) * t833 + t841 * t871 + t842 * t867;
t834 = (t867 * t872 + t868 * t871) * qJD(1);
t856 = qJDD(3) + qJDD(4);
t857 = qJD(3) + qJD(4);
t756 = (t833 * t857 - t804) * qJ(5) + (t833 * t834 + t856) * pkin(4) + t765;
t766 = t867 * t792 + t871 * t793;
t803 = -qJD(4) * t834 - t841 * t867 + t842 * t871;
t824 = pkin(4) * t857 - qJ(5) * t834;
t826 = t833 ^ 2;
t758 = -pkin(4) * t826 + qJ(5) * t803 - t824 * t857 + t766;
t862 = sin(pkin(11));
t864 = cos(pkin(11));
t815 = t833 * t864 - t834 * t862;
t896 = 2 * qJD(5);
t753 = t862 * t756 + t864 * t758 + t815 * t896;
t777 = t803 * t864 - t804 * t862;
t816 = t833 * t862 + t834 * t864;
t790 = -mrSges(6,1) * t815 + mrSges(6,2) * t816;
t806 = mrSges(6,1) * t857 - mrSges(6,3) * t816;
t791 = -pkin(5) * t815 - pkin(9) * t816;
t855 = t857 ^ 2;
t750 = -pkin(5) * t855 + pkin(9) * t856 + t791 * t815 + t753;
t821 = t865 * t838 - t863 * t840;
t882 = -qJDD(1) * pkin(2) - t821;
t794 = -t842 * pkin(3) + t846 * t894 + (-pkin(8) * t860 - pkin(7)) * t874 + t882;
t760 = -t803 * pkin(4) - t826 * qJ(5) + t834 * t824 + qJDD(5) + t794;
t778 = t803 * t862 + t804 * t864;
t754 = (-t815 * t857 - t778) * pkin(9) + (t816 * t857 - t777) * pkin(5) + t760;
t866 = sin(qJ(6));
t870 = cos(qJ(6));
t747 = -t750 * t866 + t754 * t870;
t799 = -t816 * t866 + t857 * t870;
t764 = qJD(6) * t799 + t778 * t870 + t856 * t866;
t776 = qJDD(6) - t777;
t800 = t816 * t870 + t857 * t866;
t779 = -mrSges(7,1) * t799 + mrSges(7,2) * t800;
t809 = qJD(6) - t815;
t780 = -mrSges(7,2) * t809 + mrSges(7,3) * t799;
t743 = m(7) * t747 + mrSges(7,1) * t776 - mrSges(7,3) * t764 - t779 * t800 + t780 * t809;
t748 = t750 * t870 + t754 * t866;
t763 = -qJD(6) * t800 - t778 * t866 + t856 * t870;
t781 = mrSges(7,1) * t809 - mrSges(7,3) * t800;
t744 = m(7) * t748 - mrSges(7,2) * t776 + mrSges(7,3) * t763 + t779 * t799 - t781 * t809;
t885 = -t743 * t866 + t870 * t744;
t729 = m(6) * t753 - mrSges(6,2) * t856 + mrSges(6,3) * t777 + t790 * t815 - t806 * t857 + t885;
t884 = -t864 * t756 + t862 * t758;
t752 = -0.2e1 * qJD(5) * t816 - t884;
t805 = -mrSges(6,2) * t857 + mrSges(6,3) * t815;
t749 = -t856 * pkin(5) - t855 * pkin(9) + (t896 + t791) * t816 + t884;
t881 = -m(7) * t749 + t763 * mrSges(7,1) - mrSges(7,2) * t764 + t799 * t780 - t781 * t800;
t739 = m(6) * t752 + mrSges(6,1) * t856 - mrSges(6,3) * t778 - t790 * t816 + t805 * t857 + t881;
t723 = t862 * t729 + t864 * t739;
t819 = -mrSges(5,1) * t833 + mrSges(5,2) * t834;
t823 = -mrSges(5,2) * t857 + mrSges(5,3) * t833;
t720 = m(5) * t765 + mrSges(5,1) * t856 - mrSges(5,3) * t804 - t819 * t834 + t823 * t857 + t723;
t825 = mrSges(5,1) * t857 - mrSges(5,3) * t834;
t886 = t864 * t729 - t739 * t862;
t721 = m(5) * t766 - mrSges(5,2) * t856 + mrSges(5,3) * t803 + t819 * t833 - t825 * t857 + t886;
t714 = t871 * t720 + t867 * t721;
t831 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t868 + Ifges(4,2) * t872) * qJD(1);
t832 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t868 + Ifges(4,4) * t872) * qJD(1);
t767 = Ifges(7,5) * t800 + Ifges(7,6) * t799 + Ifges(7,3) * t809;
t769 = Ifges(7,1) * t800 + Ifges(7,4) * t799 + Ifges(7,5) * t809;
t736 = -mrSges(7,1) * t749 + mrSges(7,3) * t748 + Ifges(7,4) * t764 + Ifges(7,2) * t763 + Ifges(7,6) * t776 - t767 * t800 + t769 * t809;
t768 = Ifges(7,4) * t800 + Ifges(7,2) * t799 + Ifges(7,6) * t809;
t737 = mrSges(7,2) * t749 - mrSges(7,3) * t747 + Ifges(7,1) * t764 + Ifges(7,4) * t763 + Ifges(7,5) * t776 + t767 * t799 - t768 * t809;
t783 = Ifges(6,4) * t816 + Ifges(6,2) * t815 + Ifges(6,6) * t857;
t784 = Ifges(6,1) * t816 + Ifges(6,4) * t815 + Ifges(6,5) * t857;
t811 = Ifges(5,4) * t834 + Ifges(5,2) * t833 + Ifges(5,6) * t857;
t812 = Ifges(5,1) * t834 + Ifges(5,4) * t833 + Ifges(5,5) * t857;
t877 = -mrSges(5,1) * t765 - mrSges(6,1) * t752 + mrSges(5,2) * t766 + mrSges(6,2) * t753 - pkin(4) * t723 - pkin(5) * t881 - pkin(9) * t885 - t870 * t736 - t866 * t737 - t816 * t783 + t815 * t784 + t833 * t812 - Ifges(6,6) * t777 - Ifges(6,5) * t778 - t834 * t811 - Ifges(5,6) * t803 - Ifges(5,5) * t804 + (-Ifges(5,3) - Ifges(6,3)) * t856;
t897 = mrSges(4,1) * t801 - mrSges(4,2) * t802 + Ifges(4,5) * t841 + Ifges(4,6) * t842 + Ifges(4,3) * qJDD(3) + pkin(3) * t714 + (t831 * t868 - t832 * t872) * qJD(1) - t877;
t839 = (-mrSges(4,1) * t872 + mrSges(4,2) * t868) * qJD(1);
t893 = qJD(1) * t872;
t845 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t893;
t712 = m(4) * t801 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t841 + qJD(3) * t845 - t839 * t894 + t714;
t844 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t894;
t887 = -t720 * t867 + t871 * t721;
t713 = m(4) * t802 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t842 - qJD(3) * t844 + t839 * t893 + t887;
t888 = -t712 * t868 + t872 * t713;
t704 = m(3) * t822 - mrSges(3,1) * t874 - qJDD(1) * mrSges(3,2) + t888;
t817 = -t874 * pkin(7) + t882;
t732 = t870 * t743 + t866 * t744;
t730 = m(6) * t760 - t777 * mrSges(6,1) + t778 * mrSges(6,2) - t815 * t805 + t816 * t806 + t732;
t879 = m(5) * t794 - t803 * mrSges(5,1) + mrSges(5,2) * t804 - t833 * t823 + t825 * t834 + t730;
t876 = -m(4) * t817 + t842 * mrSges(4,1) - mrSges(4,2) * t841 - t844 * t894 + t845 * t893 - t879;
t725 = m(3) * t821 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t874 + t876;
t700 = t863 * t704 + t865 * t725;
t697 = m(2) * t847 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t874 + t700;
t889 = t865 * t704 - t725 * t863;
t698 = m(2) * t848 - mrSges(2,1) * t874 - qJDD(1) * mrSges(2,2) + t889;
t895 = t873 * t697 + t869 * t698;
t707 = t872 * t712 + t868 * t713;
t705 = m(3) * t861 + t707;
t890 = -t697 * t869 + t873 * t698;
t782 = Ifges(6,5) * t816 + Ifges(6,6) * t815 + Ifges(6,3) * t857;
t715 = mrSges(6,2) * t760 - mrSges(6,3) * t752 + Ifges(6,1) * t778 + Ifges(6,4) * t777 + Ifges(6,5) * t856 - pkin(9) * t732 - t736 * t866 + t737 * t870 + t782 * t815 - t783 * t857;
t878 = mrSges(7,1) * t747 - mrSges(7,2) * t748 + Ifges(7,5) * t764 + Ifges(7,6) * t763 + Ifges(7,3) * t776 + t768 * t800 - t769 * t799;
t716 = -mrSges(6,1) * t760 + mrSges(6,3) * t753 + Ifges(6,4) * t778 + Ifges(6,2) * t777 + Ifges(6,6) * t856 - pkin(5) * t732 - t782 * t816 + t784 * t857 - t878;
t810 = Ifges(5,5) * t834 + Ifges(5,6) * t833 + Ifges(5,3) * t857;
t701 = -mrSges(5,1) * t794 + mrSges(5,3) * t766 + Ifges(5,4) * t804 + Ifges(5,2) * t803 + Ifges(5,6) * t856 - pkin(4) * t730 + qJ(5) * t886 + t862 * t715 + t864 * t716 - t834 * t810 + t857 * t812;
t708 = mrSges(5,2) * t794 - mrSges(5,3) * t765 + Ifges(5,1) * t804 + Ifges(5,4) * t803 + Ifges(5,5) * t856 - qJ(5) * t723 + t715 * t864 - t716 * t862 + t810 * t833 - t811 * t857;
t830 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t868 + Ifges(4,6) * t872) * qJD(1);
t690 = -mrSges(4,1) * t817 + mrSges(4,3) * t802 + Ifges(4,4) * t841 + Ifges(4,2) * t842 + Ifges(4,6) * qJDD(3) - pkin(3) * t879 + pkin(8) * t887 + qJD(3) * t832 + t871 * t701 + t867 * t708 - t830 * t894;
t692 = mrSges(4,2) * t817 - mrSges(4,3) * t801 + Ifges(4,1) * t841 + Ifges(4,4) * t842 + Ifges(4,5) * qJDD(3) - pkin(8) * t714 - qJD(3) * t831 - t701 * t867 + t708 * t871 + t830 * t893;
t880 = mrSges(2,1) * t847 + mrSges(3,1) * t821 - mrSges(2,2) * t848 - mrSges(3,2) * t822 + pkin(1) * t700 + pkin(2) * t876 + pkin(7) * t888 + t872 * t690 + t868 * t692 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t693 = -mrSges(3,1) * t861 + mrSges(3,3) * t822 + t874 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t707 - t897;
t688 = mrSges(3,2) * t861 - mrSges(3,3) * t821 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t874 - pkin(7) * t707 - t690 * t868 + t692 * t872;
t687 = -mrSges(2,2) * g(3) - mrSges(2,3) * t847 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t874 - qJ(2) * t700 + t688 * t865 - t693 * t863;
t686 = mrSges(2,1) * g(3) + mrSges(2,3) * t848 + t874 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t705 + qJ(2) * t889 + t863 * t688 + t865 * t693;
t1 = [-m(1) * g(1) + t890; -m(1) * g(2) + t895; (-m(1) - m(2)) * g(3) + t705; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t895 - t869 * t686 + t873 * t687; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t890 + t873 * t686 + t869 * t687; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t880; t880; t705; t897; -t877; t730; t878;];
tauJB  = t1;
