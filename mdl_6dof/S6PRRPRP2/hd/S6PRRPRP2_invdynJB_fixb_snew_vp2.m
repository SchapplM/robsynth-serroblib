% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:44:14
% EndTime: 2019-05-05 03:44:24
% DurationCPUTime: 9.98s
% Computational Cost: add. (158622->317), mult. (338525->400), div. (0->0), fcn. (237870->12), ass. (0->136)
t914 = -2 * qJD(4);
t913 = Ifges(6,1) + Ifges(7,1);
t906 = Ifges(6,4) - Ifges(7,5);
t905 = -Ifges(6,5) - Ifges(7,4);
t912 = Ifges(6,2) + Ifges(7,3);
t904 = Ifges(6,6) - Ifges(7,6);
t911 = -Ifges(6,3) - Ifges(7,2);
t862 = sin(pkin(10));
t865 = cos(pkin(10));
t850 = g(1) * t862 - g(2) * t865;
t851 = -g(1) * t865 - g(2) * t862;
t860 = -g(3) + qJDD(1);
t871 = cos(qJ(2));
t866 = cos(pkin(6));
t869 = sin(qJ(2));
t899 = t866 * t869;
t863 = sin(pkin(6));
t900 = t863 * t869;
t815 = t850 * t899 + t871 * t851 + t860 * t900;
t873 = qJD(2) ^ 2;
t810 = -pkin(2) * t873 + qJDD(2) * pkin(8) + t815;
t831 = -t850 * t863 + t860 * t866;
t868 = sin(qJ(3));
t870 = cos(qJ(3));
t784 = -t868 * t810 + t870 * t831;
t890 = qJD(2) * qJD(3);
t888 = t870 * t890;
t848 = qJDD(2) * t868 + t888;
t781 = (-t848 + t888) * qJ(4) + (t868 * t870 * t873 + qJDD(3)) * pkin(3) + t784;
t785 = t870 * t810 + t868 * t831;
t849 = qJDD(2) * t870 - t868 * t890;
t893 = qJD(2) * t868;
t852 = qJD(3) * pkin(3) - qJ(4) * t893;
t859 = t870 ^ 2;
t782 = -pkin(3) * t859 * t873 + qJ(4) * t849 - qJD(3) * t852 + t785;
t861 = sin(pkin(11));
t864 = cos(pkin(11));
t837 = (t861 * t870 + t864 * t868) * qJD(2);
t774 = t864 * t781 - t861 * t782 + t837 * t914;
t836 = (t861 * t868 - t864 * t870) * qJD(2);
t814 = -t869 * t851 + (t850 * t866 + t860 * t863) * t871;
t775 = t861 * t781 + t864 * t782 + t836 * t914;
t818 = pkin(4) * t836 - pkin(9) * t837;
t872 = qJD(3) ^ 2;
t773 = -pkin(4) * t872 + qJDD(3) * pkin(9) - t818 * t836 + t775;
t877 = -qJDD(2) * pkin(2) - t814;
t783 = -t849 * pkin(3) + qJDD(4) + t852 * t893 + (-qJ(4) * t859 - pkin(8)) * t873 + t877;
t823 = -t848 * t861 + t849 * t864;
t824 = t848 * t864 + t849 * t861;
t777 = (qJD(3) * t836 - t824) * pkin(9) + (qJD(3) * t837 - t823) * pkin(4) + t783;
t867 = sin(qJ(5));
t908 = cos(qJ(5));
t770 = t773 * t908 + t867 * t777;
t826 = t867 * qJD(3) + t837 * t908;
t796 = qJD(5) * t826 - qJDD(3) * t908 + t824 * t867;
t835 = qJD(5) + t836;
t807 = mrSges(6,1) * t835 - mrSges(6,3) * t826;
t822 = qJDD(5) - t823;
t825 = -qJD(3) * t908 + t837 * t867;
t800 = pkin(5) * t825 - qJ(6) * t826;
t834 = t835 ^ 2;
t766 = -pkin(5) * t834 + qJ(6) * t822 + 0.2e1 * qJD(6) * t835 - t800 * t825 + t770;
t808 = -mrSges(7,1) * t835 + mrSges(7,2) * t826;
t889 = m(7) * t766 + t822 * mrSges(7,3) + t835 * t808;
t801 = mrSges(7,1) * t825 - mrSges(7,3) * t826;
t894 = -mrSges(6,1) * t825 - mrSges(6,2) * t826 - t801;
t907 = -mrSges(6,3) - mrSges(7,2);
t758 = m(6) * t770 - t822 * mrSges(6,2) + t796 * t907 - t835 * t807 + t825 * t894 + t889;
t769 = -t867 * t773 + t777 * t908;
t797 = -t825 * qJD(5) + t867 * qJDD(3) + t824 * t908;
t806 = -mrSges(6,2) * t835 - mrSges(6,3) * t825;
t767 = -t822 * pkin(5) - t834 * qJ(6) + t826 * t800 + qJDD(6) - t769;
t805 = -mrSges(7,2) * t825 + mrSges(7,3) * t835;
t882 = -m(7) * t767 + t822 * mrSges(7,1) + t835 * t805;
t760 = m(6) * t769 + t822 * mrSges(6,1) + t797 * t907 + t835 * t806 + t826 * t894 + t882;
t753 = t758 * t908 - t760 * t867;
t817 = mrSges(5,1) * t836 + mrSges(5,2) * t837;
t830 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t837;
t748 = m(5) * t775 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t823 - qJD(3) * t830 - t817 * t836 + t753;
t772 = -qJDD(3) * pkin(4) - t872 * pkin(9) + t837 * t818 - t774;
t768 = -0.2e1 * qJD(6) * t826 + (t825 * t835 - t797) * qJ(6) + (t826 * t835 + t796) * pkin(5) + t772;
t764 = m(7) * t768 + t796 * mrSges(7,1) - t797 * mrSges(7,3) + t805 * t825 - t826 * t808;
t761 = -m(6) * t772 - t796 * mrSges(6,1) - t797 * mrSges(6,2) - t825 * t806 - t807 * t826 - t764;
t829 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t836;
t755 = m(5) * t774 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t824 + qJD(3) * t829 - t817 * t837 + t761;
t742 = t861 * t748 + t864 * t755;
t895 = t906 * t825 - t913 * t826 + t905 * t835;
t897 = t904 * t825 + t905 * t826 + t911 * t835;
t749 = -mrSges(6,1) * t772 - mrSges(7,1) * t768 + mrSges(7,2) * t766 + mrSges(6,3) * t770 - pkin(5) * t764 - t912 * t796 + t906 * t797 + t904 * t822 + t897 * t826 - t895 * t835;
t896 = t912 * t825 - t906 * t826 - t904 * t835;
t750 = mrSges(6,2) * t772 + mrSges(7,2) * t767 - mrSges(6,3) * t769 - mrSges(7,3) * t768 - qJ(6) * t764 - t906 * t796 + t913 * t797 - t905 * t822 + t897 * t825 + t896 * t835;
t812 = Ifges(5,4) * t837 - Ifges(5,2) * t836 + Ifges(5,6) * qJD(3);
t813 = Ifges(5,1) * t837 - Ifges(5,4) * t836 + Ifges(5,5) * qJD(3);
t840 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t868 + Ifges(4,2) * t870) * qJD(2);
t841 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t868 + Ifges(4,4) * t870) * qJD(2);
t910 = (t840 * t868 - t841 * t870) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t784 + mrSges(5,1) * t774 - mrSges(4,2) * t785 - mrSges(5,2) * t775 + Ifges(4,5) * t848 + Ifges(5,5) * t824 + Ifges(4,6) * t849 + Ifges(5,6) * t823 + pkin(3) * t742 + pkin(4) * t761 + pkin(9) * t753 + t749 * t908 + t867 * t750 + t837 * t812 + t836 * t813;
t763 = t797 * mrSges(7,2) + t826 * t801 - t882;
t909 = -t796 * t904 - t797 * t905 - t911 * t822 - t825 * t895 - t826 * t896 + mrSges(6,1) * t769 - mrSges(7,1) * t767 - mrSges(6,2) * t770 + mrSges(7,3) * t766 - pkin(5) * t763 + qJ(6) * (-t796 * mrSges(7,2) - t825 * t801 + t889);
t752 = t867 * t758 + t760 * t908;
t751 = m(5) * t783 - t823 * mrSges(5,1) + mrSges(5,2) * t824 + t836 * t829 + t830 * t837 + t752;
t809 = -t873 * pkin(8) + t877;
t853 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t893;
t892 = qJD(2) * t870;
t854 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t892;
t875 = -m(4) * t809 + t849 * mrSges(4,1) - mrSges(4,2) * t848 - t853 * t893 + t854 * t892 - t751;
t745 = m(3) * t814 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t873 + t875;
t901 = t745 * t871;
t847 = (-mrSges(4,1) * t870 + mrSges(4,2) * t868) * qJD(2);
t740 = m(4) * t784 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t848 + qJD(3) * t854 - t847 * t893 + t742;
t885 = t864 * t748 - t755 * t861;
t741 = m(4) * t785 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t849 - qJD(3) * t853 + t847 * t892 + t885;
t886 = -t740 * t868 + t870 * t741;
t731 = m(3) * t815 - mrSges(3,1) * t873 - qJDD(2) * mrSges(3,2) + t886;
t734 = t870 * t740 + t868 * t741;
t733 = m(3) * t831 + t734;
t722 = t731 * t899 - t733 * t863 + t866 * t901;
t720 = m(2) * t850 + t722;
t727 = t871 * t731 - t745 * t869;
t726 = m(2) * t851 + t727;
t898 = t865 * t720 + t862 * t726;
t721 = t731 * t900 + t866 * t733 + t863 * t901;
t887 = -t720 * t862 + t865 * t726;
t883 = m(2) * t860 + t721;
t811 = Ifges(5,5) * t837 - Ifges(5,6) * t836 + Ifges(5,3) * qJD(3);
t735 = mrSges(5,2) * t783 - mrSges(5,3) * t774 + Ifges(5,1) * t824 + Ifges(5,4) * t823 + Ifges(5,5) * qJDD(3) - pkin(9) * t752 - qJD(3) * t812 - t867 * t749 + t750 * t908 - t836 * t811;
t736 = -mrSges(5,1) * t783 + mrSges(5,3) * t775 + Ifges(5,4) * t824 + Ifges(5,2) * t823 + Ifges(5,6) * qJDD(3) - pkin(4) * t752 + qJD(3) * t813 - t837 * t811 - t909;
t839 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t868 + Ifges(4,6) * t870) * qJD(2);
t718 = -mrSges(4,1) * t809 + mrSges(4,3) * t785 + Ifges(4,4) * t848 + Ifges(4,2) * t849 + Ifges(4,6) * qJDD(3) - pkin(3) * t751 + qJ(4) * t885 + qJD(3) * t841 + t861 * t735 + t864 * t736 - t839 * t893;
t723 = mrSges(4,2) * t809 - mrSges(4,3) * t784 + Ifges(4,1) * t848 + Ifges(4,4) * t849 + Ifges(4,5) * qJDD(3) - qJ(4) * t742 - qJD(3) * t840 + t735 * t864 - t736 * t861 + t839 * t892;
t716 = mrSges(3,2) * t831 - mrSges(3,3) * t814 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t873 - pkin(8) * t734 - t718 * t868 + t723 * t870;
t717 = -mrSges(3,1) * t831 + mrSges(3,3) * t815 + t873 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t734 - t910;
t878 = pkin(7) * t727 + t716 * t869 + t717 * t871;
t715 = mrSges(3,1) * t814 - mrSges(3,2) * t815 + Ifges(3,3) * qJDD(2) + pkin(2) * t875 + pkin(8) * t886 + t870 * t718 + t868 * t723;
t714 = mrSges(2,2) * t860 - mrSges(2,3) * t850 + t871 * t716 - t869 * t717 + (-t721 * t863 - t722 * t866) * pkin(7);
t713 = -mrSges(2,1) * t860 + mrSges(2,3) * t851 - pkin(1) * t721 - t863 * t715 + t866 * t878;
t1 = [-m(1) * g(1) + t887; -m(1) * g(2) + t898; -m(1) * g(3) + t883; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t898 - t862 * t713 + t865 * t714; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t887 + t865 * t713 + t862 * t714; -mrSges(1,1) * g(2) + mrSges(2,1) * t850 + mrSges(1,2) * g(1) - mrSges(2,2) * t851 + pkin(1) * t722 + t866 * t715 + t863 * t878; t883; t715; t910; t751; t909; t763;];
tauJB  = t1;
