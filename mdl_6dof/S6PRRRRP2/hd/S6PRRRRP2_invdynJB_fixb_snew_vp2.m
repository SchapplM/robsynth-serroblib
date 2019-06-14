% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRP2
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
% Datum: 2019-05-05 09:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:35:28
% EndTime: 2019-05-05 09:35:39
% DurationCPUTime: 10.78s
% Computational Cost: add. (184723->319), mult. (360853->400), div. (0->0), fcn. (256866->12), ass. (0->138)
t920 = Ifges(6,1) + Ifges(7,1);
t913 = Ifges(6,4) - Ifges(7,5);
t912 = -Ifges(6,5) - Ifges(7,4);
t919 = Ifges(6,2) + Ifges(7,3);
t911 = Ifges(6,6) - Ifges(7,6);
t918 = -Ifges(6,3) - Ifges(7,2);
t873 = sin(qJ(4));
t874 = sin(qJ(3));
t876 = cos(qJ(4));
t877 = cos(qJ(3));
t843 = (t873 * t874 - t876 * t877) * qJD(2);
t868 = sin(pkin(11));
t870 = cos(pkin(11));
t853 = g(1) * t868 - g(2) * t870;
t854 = -g(1) * t870 - g(2) * t868;
t867 = -g(3) + qJDD(1);
t869 = sin(pkin(6));
t871 = cos(pkin(6));
t875 = sin(qJ(2));
t878 = cos(qJ(2));
t825 = -t875 * t854 + (t853 * t871 + t867 * t869) * t878;
t907 = t871 * t875;
t908 = t869 * t875;
t826 = t853 * t907 + t878 * t854 + t867 * t908;
t879 = qJD(2) ^ 2;
t821 = -pkin(2) * t879 + qJDD(2) * pkin(8) + t826;
t836 = -t853 * t869 + t867 * t871;
t798 = -t874 * t821 + t877 * t836;
t899 = qJD(2) * qJD(3);
t897 = t877 * t899;
t851 = qJDD(2) * t874 + t897;
t783 = (-t851 + t897) * pkin(9) + (t874 * t877 * t879 + qJDD(3)) * pkin(3) + t798;
t799 = t877 * t821 + t874 * t836;
t852 = qJDD(2) * t877 - t874 * t899;
t901 = qJD(2) * t874;
t858 = qJD(3) * pkin(3) - pkin(9) * t901;
t866 = t877 ^ 2;
t784 = -pkin(3) * t866 * t879 + pkin(9) * t852 - qJD(3) * t858 + t799;
t779 = t873 * t783 + t876 * t784;
t844 = (t873 * t877 + t874 * t876) * qJD(2);
t814 = -qJD(4) * t844 - t851 * t873 + t852 * t876;
t828 = mrSges(5,1) * t843 + mrSges(5,2) * t844;
t865 = qJD(3) + qJD(4);
t835 = mrSges(5,1) * t865 - mrSges(5,3) * t844;
t864 = qJDD(3) + qJDD(4);
t829 = pkin(4) * t843 - pkin(10) * t844;
t863 = t865 ^ 2;
t774 = -pkin(4) * t863 + pkin(10) * t864 - t829 * t843 + t779;
t886 = -qJDD(2) * pkin(2) - t825;
t789 = -t852 * pkin(3) + t858 * t901 + (-pkin(9) * t866 - pkin(8)) * t879 + t886;
t815 = -qJD(4) * t843 + t851 * t876 + t852 * t873;
t776 = (t843 * t865 - t815) * pkin(10) + (t844 * t865 - t814) * pkin(4) + t789;
t872 = sin(qJ(5));
t915 = cos(qJ(5));
t771 = t915 * t774 + t872 * t776;
t831 = t844 * t915 + t872 * t865;
t787 = qJD(5) * t831 + t815 * t872 - t864 * t915;
t812 = qJDD(5) - t814;
t838 = qJD(5) + t843;
t818 = mrSges(6,1) * t838 - mrSges(6,3) * t831;
t830 = t844 * t872 - t865 * t915;
t802 = pkin(5) * t830 - qJ(6) * t831;
t837 = t838 ^ 2;
t767 = -pkin(5) * t837 + qJ(6) * t812 + 0.2e1 * qJD(6) * t838 - t802 * t830 + t771;
t819 = -mrSges(7,1) * t838 + mrSges(7,2) * t831;
t898 = m(7) * t767 + t812 * mrSges(7,3) + t838 * t819;
t803 = mrSges(7,1) * t830 - mrSges(7,3) * t831;
t902 = -mrSges(6,1) * t830 - mrSges(6,2) * t831 - t803;
t914 = -mrSges(6,3) - mrSges(7,2);
t758 = m(6) * t771 - t812 * mrSges(6,2) + t787 * t914 - t838 * t818 + t830 * t902 + t898;
t770 = -t872 * t774 + t776 * t915;
t788 = -t830 * qJD(5) + t815 * t915 + t872 * t864;
t817 = -mrSges(6,2) * t838 - mrSges(6,3) * t830;
t768 = -t812 * pkin(5) - t837 * qJ(6) + t831 * t802 + qJDD(6) - t770;
t816 = -mrSges(7,2) * t830 + mrSges(7,3) * t838;
t891 = -m(7) * t768 + t812 * mrSges(7,1) + t838 * t816;
t760 = m(6) * t770 + t812 * mrSges(6,1) + t788 * t914 + t838 * t817 + t831 * t902 + t891;
t893 = t915 * t758 - t760 * t872;
t746 = m(5) * t779 - mrSges(5,2) * t864 + mrSges(5,3) * t814 - t828 * t843 - t835 * t865 + t893;
t778 = t876 * t783 - t873 * t784;
t834 = -mrSges(5,2) * t865 - mrSges(5,3) * t843;
t773 = -t864 * pkin(4) - t863 * pkin(10) + t844 * t829 - t778;
t769 = -0.2e1 * qJD(6) * t831 + (t830 * t838 - t788) * qJ(6) + (t831 * t838 + t787) * pkin(5) + t773;
t765 = m(7) * t769 + mrSges(7,1) * t787 - t788 * mrSges(7,3) + t816 * t830 - t831 * t819;
t882 = -m(6) * t773 - t787 * mrSges(6,1) - mrSges(6,2) * t788 - t830 * t817 - t818 * t831 - t765;
t755 = m(5) * t778 + mrSges(5,1) * t864 - mrSges(5,3) * t815 - t828 * t844 + t834 * t865 + t882;
t740 = t873 * t746 + t876 * t755;
t841 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t874 + Ifges(4,2) * t877) * qJD(2);
t842 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t874 + Ifges(4,4) * t877) * qJD(2);
t903 = t913 * t830 - t831 * t920 + t912 * t838;
t905 = t830 * t911 + t831 * t912 + t838 * t918;
t748 = -mrSges(6,1) * t773 - mrSges(7,1) * t769 + mrSges(7,2) * t767 + mrSges(6,3) * t771 - pkin(5) * t765 - t787 * t919 + t913 * t788 + t911 * t812 + t905 * t831 - t903 * t838;
t904 = t830 * t919 - t831 * t913 - t838 * t911;
t750 = mrSges(6,2) * t773 + mrSges(7,2) * t768 - mrSges(6,3) * t770 - mrSges(7,3) * t769 - qJ(6) * t765 - t913 * t787 + t788 * t920 - t912 * t812 + t905 * t830 + t904 * t838;
t823 = Ifges(5,4) * t844 - Ifges(5,2) * t843 + Ifges(5,6) * t865;
t824 = Ifges(5,1) * t844 - Ifges(5,4) * t843 + Ifges(5,5) * t865;
t884 = -mrSges(5,1) * t778 + mrSges(5,2) * t779 - Ifges(5,5) * t815 - Ifges(5,6) * t814 - Ifges(5,3) * t864 - pkin(4) * t882 - pkin(10) * t893 - t915 * t748 - t872 * t750 - t844 * t823 - t843 * t824;
t917 = mrSges(4,1) * t798 - mrSges(4,2) * t799 + Ifges(4,5) * t851 + Ifges(4,6) * t852 + Ifges(4,3) * qJDD(3) + pkin(3) * t740 + (t841 * t874 - t842 * t877) * qJD(2) - t884;
t764 = t788 * mrSges(7,2) + t831 * t803 - t891;
t916 = -t787 * t911 - t788 * t912 - t918 * t812 - t830 * t903 - t831 * t904 + mrSges(6,1) * t770 - mrSges(7,1) * t768 - mrSges(6,2) * t771 + mrSges(7,3) * t767 - pkin(5) * t764 + qJ(6) * (-t787 * mrSges(7,2) - t830 * t803 + t898);
t820 = -t879 * pkin(8) + t886;
t855 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t901;
t900 = qJD(2) * t877;
t856 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t900;
t752 = t872 * t758 + t915 * t760;
t885 = m(5) * t789 - t814 * mrSges(5,1) + mrSges(5,2) * t815 + t843 * t834 + t835 * t844 + t752;
t881 = -m(4) * t820 + t852 * mrSges(4,1) - mrSges(4,2) * t851 - t855 * t901 + t856 * t900 - t885;
t743 = m(3) * t825 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t879 + t881;
t909 = t743 * t878;
t850 = (-mrSges(4,1) * t877 + mrSges(4,2) * t874) * qJD(2);
t738 = m(4) * t798 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t851 + qJD(3) * t856 - t850 * t901 + t740;
t894 = t876 * t746 - t755 * t873;
t739 = m(4) * t799 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t852 - qJD(3) * t855 + t850 * t900 + t894;
t895 = -t738 * t874 + t877 * t739;
t729 = m(3) * t826 - mrSges(3,1) * t879 - qJDD(2) * mrSges(3,2) + t895;
t732 = t877 * t738 + t874 * t739;
t731 = m(3) * t836 + t732;
t720 = t729 * t907 - t731 * t869 + t871 * t909;
t718 = m(2) * t853 + t720;
t725 = t878 * t729 - t743 * t875;
t724 = m(2) * t854 + t725;
t906 = t870 * t718 + t868 * t724;
t719 = t729 * t908 + t871 * t731 + t869 * t909;
t896 = -t718 * t868 + t870 * t724;
t892 = m(2) * t867 + t719;
t822 = Ifges(5,5) * t844 - Ifges(5,6) * t843 + Ifges(5,3) * t865;
t733 = mrSges(5,2) * t789 - mrSges(5,3) * t778 + Ifges(5,1) * t815 + Ifges(5,4) * t814 + Ifges(5,5) * t864 - pkin(10) * t752 - t872 * t748 + t750 * t915 - t843 * t822 - t865 * t823;
t734 = -mrSges(5,1) * t789 + mrSges(5,3) * t779 + Ifges(5,4) * t815 + Ifges(5,2) * t814 + Ifges(5,6) * t864 - pkin(4) * t752 - t844 * t822 + t865 * t824 - t916;
t840 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t874 + Ifges(4,6) * t877) * qJD(2);
t716 = -mrSges(4,1) * t820 + mrSges(4,3) * t799 + Ifges(4,4) * t851 + Ifges(4,2) * t852 + Ifges(4,6) * qJDD(3) - pkin(3) * t885 + pkin(9) * t894 + qJD(3) * t842 + t873 * t733 + t876 * t734 - t840 * t901;
t721 = mrSges(4,2) * t820 - mrSges(4,3) * t798 + Ifges(4,1) * t851 + Ifges(4,4) * t852 + Ifges(4,5) * qJDD(3) - pkin(9) * t740 - qJD(3) * t841 + t733 * t876 - t734 * t873 + t840 * t900;
t714 = mrSges(3,2) * t836 - mrSges(3,3) * t825 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t879 - pkin(8) * t732 - t716 * t874 + t721 * t877;
t715 = -mrSges(3,1) * t836 + mrSges(3,3) * t826 + t879 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t732 - t917;
t887 = pkin(7) * t725 + t714 * t875 + t715 * t878;
t713 = mrSges(3,1) * t825 - mrSges(3,2) * t826 + Ifges(3,3) * qJDD(2) + pkin(2) * t881 + pkin(8) * t895 + t877 * t716 + t874 * t721;
t712 = mrSges(2,2) * t867 - mrSges(2,3) * t853 + t878 * t714 - t875 * t715 + (-t719 * t869 - t720 * t871) * pkin(7);
t711 = -mrSges(2,1) * t867 + mrSges(2,3) * t854 - pkin(1) * t719 - t869 * t713 + t871 * t887;
t1 = [-m(1) * g(1) + t896; -m(1) * g(2) + t906; -m(1) * g(3) + t892; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t906 - t868 * t711 + t870 * t712; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t896 + t870 * t711 + t868 * t712; -mrSges(1,1) * g(2) + mrSges(2,1) * t853 + mrSges(1,2) * g(1) - mrSges(2,2) * t854 + pkin(1) * t720 + t871 * t713 + t869 * t887; t892; t713; t917; -t884; t916; t764;];
tauJB  = t1;
