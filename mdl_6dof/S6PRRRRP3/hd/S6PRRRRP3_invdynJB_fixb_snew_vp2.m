% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:43:34
% EndTime: 2019-05-05 09:43:47
% DurationCPUTime: 12.14s
% Computational Cost: add. (203909->319), mult. (396161->395), div. (0->0), fcn. (282193->12), ass. (0->137)
t906 = Ifges(6,4) + Ifges(7,4);
t915 = Ifges(6,2) + Ifges(7,2);
t911 = Ifges(6,6) + Ifges(7,6);
t872 = sin(qJ(4));
t876 = cos(qJ(4));
t873 = sin(qJ(3));
t899 = qJD(2) * t873;
t852 = qJD(3) * t872 + t876 * t899;
t877 = cos(qJ(3));
t897 = qJD(2) * qJD(3);
t894 = t877 * t897;
t855 = qJDD(2) * t873 + t894;
t825 = -qJD(4) * t852 + qJDD(3) * t876 - t855 * t872;
t851 = qJD(3) * t876 - t872 * t899;
t826 = qJD(4) * t851 + qJDD(3) * t872 + t855 * t876;
t871 = sin(qJ(5));
t875 = cos(qJ(5));
t828 = t851 * t875 - t852 * t871;
t786 = qJD(5) * t828 + t825 * t871 + t826 * t875;
t829 = t851 * t871 + t852 * t875;
t805 = -mrSges(7,1) * t828 + mrSges(7,2) * t829;
t867 = sin(pkin(11));
t869 = cos(pkin(11));
t857 = g(1) * t867 - g(2) * t869;
t858 = -g(1) * t869 - g(2) * t867;
t866 = -g(3) + qJDD(1);
t878 = cos(qJ(2));
t870 = cos(pkin(6));
t874 = sin(qJ(2));
t903 = t870 * t874;
t868 = sin(pkin(6));
t904 = t868 * t874;
t818 = t857 * t903 + t858 * t878 + t866 * t904;
t880 = qJD(2) ^ 2;
t811 = -pkin(2) * t880 + qJDD(2) * pkin(8) + t818;
t835 = -t857 * t868 + t866 * t870;
t804 = t811 * t877 + t835 * t873;
t854 = (-pkin(3) * t877 - pkin(9) * t873) * qJD(2);
t879 = qJD(3) ^ 2;
t898 = t877 * qJD(2);
t790 = -pkin(3) * t879 + qJDD(3) * pkin(9) + t854 * t898 + t804;
t817 = -t874 * t858 + (t857 * t870 + t866 * t868) * t878;
t810 = -qJDD(2) * pkin(2) - t880 * pkin(8) - t817;
t864 = t873 * t897;
t856 = qJDD(2) * t877 - t864;
t793 = (-t855 - t894) * pkin(9) + (-t856 + t864) * pkin(3) + t810;
t772 = -t872 * t790 + t793 * t876;
t848 = qJDD(4) - t856;
t863 = qJD(4) - t898;
t769 = (t851 * t863 - t826) * pkin(10) + (t851 * t852 + t848) * pkin(4) + t772;
t773 = t790 * t876 + t793 * t872;
t834 = pkin(4) * t863 - pkin(10) * t852;
t847 = t851 ^ 2;
t771 = -pkin(4) * t847 + pkin(10) * t825 - t834 * t863 + t773;
t763 = t769 * t875 - t871 * t771;
t844 = qJDD(5) + t848;
t862 = qJD(5) + t863;
t758 = -0.2e1 * qJD(6) * t829 + (t828 * t862 - t786) * qJ(6) + (t828 * t829 + t844) * pkin(5) + t763;
t812 = -mrSges(7,2) * t862 + mrSges(7,3) * t828;
t896 = m(7) * t758 + mrSges(7,1) * t844 + t812 * t862;
t755 = -t786 * mrSges(7,3) - t829 * t805 + t896;
t764 = t769 * t871 + t771 * t875;
t785 = -qJD(5) * t829 + t825 * t875 - t826 * t871;
t814 = pkin(5) * t862 - qJ(6) * t829;
t827 = t828 ^ 2;
t760 = -pkin(5) * t827 + qJ(6) * t785 + 0.2e1 * qJD(6) * t828 - t814 * t862 + t764;
t912 = Ifges(6,5) + Ifges(7,5);
t913 = Ifges(6,1) + Ifges(7,1);
t900 = -t828 * t906 - t829 * t913 - t862 * t912;
t909 = t828 * t915 + t829 * t906 + t862 * t911;
t910 = Ifges(6,3) + Ifges(7,3);
t914 = mrSges(6,1) * t763 + mrSges(7,1) * t758 - mrSges(6,2) * t764 - mrSges(7,2) * t760 + pkin(5) * t755 + t785 * t911 + t786 * t912 + t828 * t900 + t829 * t909 + t844 * t910;
t806 = -mrSges(6,1) * t828 + mrSges(6,2) * t829;
t813 = -mrSges(6,2) * t862 + mrSges(6,3) * t828;
t747 = m(6) * t763 + t844 * mrSges(6,1) + t862 * t813 + (-t805 - t806) * t829 + (-mrSges(6,3) - mrSges(7,3)) * t786 + t896;
t815 = mrSges(7,1) * t862 - mrSges(7,3) * t829;
t816 = mrSges(6,1) * t862 - mrSges(6,3) * t829;
t895 = m(7) * t760 + mrSges(7,3) * t785 + t805 * t828;
t750 = m(6) * t764 + t785 * mrSges(6,3) + t828 * t806 + (-t815 - t816) * t862 + (-mrSges(6,2) - mrSges(7,2)) * t844 + t895;
t745 = t747 * t875 + t750 * t871;
t820 = Ifges(5,4) * t852 + Ifges(5,2) * t851 + Ifges(5,6) * t863;
t821 = Ifges(5,1) * t852 + Ifges(5,4) * t851 + Ifges(5,5) * t863;
t908 = mrSges(5,1) * t772 - mrSges(5,2) * t773 + Ifges(5,5) * t826 + Ifges(5,6) * t825 + Ifges(5,3) * t848 + pkin(4) * t745 + t852 * t820 - t851 * t821 + t914;
t803 = -t811 * t873 + t835 * t877;
t789 = -qJDD(3) * pkin(3) - pkin(9) * t879 + t854 * t899 - t803;
t774 = -pkin(4) * t825 - pkin(10) * t847 + t834 * t852 + t789;
t766 = -pkin(5) * t785 - qJ(6) * t827 + t814 * t829 + qJDD(6) + t774;
t761 = m(7) * t766 - t785 * mrSges(7,1) + t786 * mrSges(7,2) - t828 * t812 + t829 * t815;
t901 = -t828 * t911 - t829 * t912 - t862 * t910;
t740 = -mrSges(6,1) * t774 + mrSges(6,3) * t764 - mrSges(7,1) * t766 + mrSges(7,3) * t760 - pkin(5) * t761 + qJ(6) * t895 + (-qJ(6) * t815 - t900) * t862 + (-mrSges(7,2) * qJ(6) + t911) * t844 + t901 * t829 + t906 * t786 + t915 * t785;
t744 = mrSges(6,2) * t774 + mrSges(7,2) * t766 - mrSges(6,3) * t763 - mrSges(7,3) * t758 - qJ(6) * t755 + t785 * t906 + t786 * t913 - t828 * t901 + t844 * t912 - t862 * t909;
t819 = Ifges(5,5) * t852 + Ifges(5,6) * t851 + Ifges(5,3) * t863;
t885 = m(6) * t774 - mrSges(6,1) * t785 + mrSges(6,2) * t786 - t813 * t828 + t816 * t829 + t761;
t890 = -t747 * t871 + t750 * t875;
t723 = -mrSges(5,1) * t789 + mrSges(5,3) * t773 + Ifges(5,4) * t826 + Ifges(5,2) * t825 + Ifges(5,6) * t848 - pkin(4) * t885 + pkin(10) * t890 + t875 * t740 + t871 * t744 - t852 * t819 + t863 * t821;
t724 = mrSges(5,2) * t789 - mrSges(5,3) * t772 + Ifges(5,1) * t826 + Ifges(5,4) * t825 + Ifges(5,5) * t848 - pkin(10) * t745 - t740 * t871 + t744 * t875 + t819 * t851 - t820 * t863;
t830 = -mrSges(5,1) * t851 + mrSges(5,2) * t852;
t832 = -mrSges(5,2) * t863 + mrSges(5,3) * t851;
t742 = m(5) * t772 + mrSges(5,1) * t848 - mrSges(5,3) * t826 - t830 * t852 + t832 * t863 + t745;
t833 = mrSges(5,1) * t863 - mrSges(5,3) * t852;
t743 = m(5) * t773 - mrSges(5,2) * t848 + mrSges(5,3) * t825 + t830 * t851 - t833 * t863 + t890;
t739 = -t742 * t872 + t743 * t876;
t753 = -m(5) * t789 + mrSges(5,1) * t825 - mrSges(5,2) * t826 + t832 * t851 - t833 * t852 - t885;
t842 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t873 + Ifges(4,2) * t877) * qJD(2);
t843 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t873 + Ifges(4,4) * t877) * qJD(2);
t907 = mrSges(4,1) * t803 - mrSges(4,2) * t804 + Ifges(4,5) * t855 + Ifges(4,6) * t856 + Ifges(4,3) * qJDD(3) + pkin(3) * t753 + pkin(9) * t739 + t876 * t723 + t872 * t724 + (t842 * t873 - t843 * t877) * qJD(2);
t738 = t742 * t876 + t743 * t872;
t859 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t899;
t860 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t898;
t884 = -m(4) * t810 + mrSges(4,1) * t856 - mrSges(4,2) * t855 - t859 * t899 + t860 * t898 - t738;
t734 = m(3) * t817 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t880 + t884;
t905 = t734 * t878;
t853 = (-mrSges(4,1) * t877 + mrSges(4,2) * t873) * qJD(2);
t737 = m(4) * t804 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t856 - qJD(3) * t859 + t853 * t898 + t739;
t752 = m(4) * t803 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t855 + qJD(3) * t860 - t853 * t899 + t753;
t891 = t737 * t877 - t752 * t873;
t728 = m(3) * t818 - mrSges(3,1) * t880 - qJDD(2) * mrSges(3,2) + t891;
t731 = t737 * t873 + t752 * t877;
t730 = m(3) * t835 + t731;
t717 = t728 * t903 - t730 * t868 + t870 * t905;
t715 = m(2) * t857 + t717;
t721 = t728 * t878 - t734 * t874;
t720 = m(2) * t858 + t721;
t902 = t715 * t869 + t720 * t867;
t716 = t728 * t904 + t730 * t870 + t868 * t905;
t892 = -t715 * t867 + t720 * t869;
t889 = m(2) * t866 + t716;
t841 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t873 + Ifges(4,6) * t877) * qJD(2);
t713 = mrSges(4,2) * t810 - mrSges(4,3) * t803 + Ifges(4,1) * t855 + Ifges(4,4) * t856 + Ifges(4,5) * qJDD(3) - pkin(9) * t738 - qJD(3) * t842 - t723 * t872 + t724 * t876 + t841 * t898;
t722 = -mrSges(4,1) * t810 + mrSges(4,3) * t804 + Ifges(4,4) * t855 + Ifges(4,2) * t856 + Ifges(4,6) * qJDD(3) - pkin(3) * t738 + qJD(3) * t843 - t841 * t899 - t908;
t711 = mrSges(3,2) * t835 - mrSges(3,3) * t817 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t880 - pkin(8) * t731 + t713 * t877 - t722 * t873;
t712 = -mrSges(3,1) * t835 + mrSges(3,3) * t818 + t880 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t731 - t907;
t886 = pkin(7) * t721 + t711 * t874 + t712 * t878;
t710 = mrSges(3,1) * t817 - mrSges(3,2) * t818 + Ifges(3,3) * qJDD(2) + pkin(2) * t884 + pkin(8) * t891 + t873 * t713 + t877 * t722;
t709 = mrSges(2,2) * t866 - mrSges(2,3) * t857 + t878 * t711 - t874 * t712 + (-t716 * t868 - t717 * t870) * pkin(7);
t708 = -mrSges(2,1) * t866 + mrSges(2,3) * t858 - pkin(1) * t716 - t868 * t710 + t870 * t886;
t1 = [-m(1) * g(1) + t892; -m(1) * g(2) + t902; -m(1) * g(3) + t889; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t902 - t708 * t867 + t709 * t869; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t892 + t869 * t708 + t867 * t709; -mrSges(1,1) * g(2) + mrSges(2,1) * t857 + mrSges(1,2) * g(1) - mrSges(2,2) * t858 + pkin(1) * t717 + t870 * t710 + t868 * t886; t889; t710; t907; t908; t914; t761;];
tauJB  = t1;
