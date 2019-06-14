% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:44:53
% EndTime: 2019-05-05 08:45:45
% DurationCPUTime: 50.05s
% Computational Cost: add. (877234->352), mult. (1840550->468), div. (0->0), fcn. (1482770->16), ass. (0->157)
t895 = sin(pkin(7));
t903 = sin(qJ(3));
t906 = cos(qJ(3));
t924 = qJD(2) * qJD(3);
t876 = (-qJDD(2) * t906 + t903 * t924) * t895;
t894 = sin(pkin(12));
t898 = cos(pkin(12));
t884 = g(1) * t894 - g(2) * t898;
t885 = -g(1) * t898 - g(2) * t894;
t892 = -g(3) + qJDD(1);
t904 = sin(qJ(2));
t900 = cos(pkin(6));
t907 = cos(qJ(2));
t927 = t900 * t907;
t896 = sin(pkin(6));
t930 = t896 * t907;
t854 = t884 * t927 - t885 * t904 + t892 * t930;
t908 = qJD(2) ^ 2;
t934 = pkin(9) * t895;
t847 = qJDD(2) * pkin(2) + t908 * t934 + t854;
t928 = t900 * t904;
t931 = t896 * t904;
t855 = t884 * t928 + t907 * t885 + t892 * t931;
t848 = -pkin(2) * t908 + qJDD(2) * t934 + t855;
t869 = -t884 * t896 + t892 * t900;
t899 = cos(pkin(7));
t818 = -t903 * t848 + (t847 * t899 + t869 * t895) * t906;
t935 = cos(qJ(4));
t890 = qJD(2) * t899 + qJD(3);
t925 = qJD(2) * t895;
t922 = t906 * t925;
t872 = -mrSges(4,2) * t890 + mrSges(4,3) * t922;
t873 = (-mrSges(4,1) * t906 + mrSges(4,2) * t903) * t925;
t875 = (qJDD(2) * t903 + t906 * t924) * t895;
t889 = qJDD(2) * t899 + qJDD(3);
t929 = t899 * t903;
t932 = t895 * t903;
t819 = t847 * t929 + t906 * t848 + t869 * t932;
t874 = (-pkin(3) * t906 - pkin(10) * t903) * t925;
t888 = t890 ^ 2;
t814 = -pkin(3) * t888 + pkin(10) * t889 + t874 * t922 + t819;
t865 = t899 * t869;
t817 = t876 * pkin(3) - t875 * pkin(10) + t865 + (-t847 + (pkin(3) * t903 - pkin(10) * t906) * t890 * qJD(2)) * t895;
t902 = sin(qJ(4));
t803 = t935 * t814 + t902 * t817;
t923 = t903 * t925;
t867 = -t935 * t890 + t902 * t923;
t868 = t902 * t890 + t923 * t935;
t849 = pkin(4) * t867 - qJ(5) * t868;
t870 = qJDD(4) + t876;
t883 = qJD(4) - t922;
t882 = t883 ^ 2;
t798 = -pkin(4) * t882 + qJ(5) * t870 - t849 * t867 + t803;
t813 = -t889 * pkin(3) - t888 * pkin(10) + t874 * t923 - t818;
t842 = qJD(4) * t868 + t875 * t902 - t935 * t889;
t843 = -t867 * qJD(4) + t875 * t935 + t902 * t889;
t801 = (t867 * t883 - t843) * qJ(5) + (t868 * t883 + t842) * pkin(4) + t813;
t893 = sin(pkin(13));
t897 = cos(pkin(13));
t857 = t868 * t897 + t883 * t893;
t793 = -0.2e1 * qJD(5) * t857 - t893 * t798 + t897 * t801;
t830 = t843 * t897 + t870 * t893;
t856 = -t868 * t893 + t883 * t897;
t791 = (t856 * t867 - t830) * pkin(11) + (t856 * t857 + t842) * pkin(5) + t793;
t794 = 0.2e1 * qJD(5) * t856 + t897 * t798 + t893 * t801;
t829 = -t843 * t893 + t870 * t897;
t835 = pkin(5) * t867 - pkin(11) * t857;
t853 = t856 ^ 2;
t792 = -pkin(5) * t853 + pkin(11) * t829 - t835 * t867 + t794;
t901 = sin(qJ(6));
t905 = cos(qJ(6));
t789 = t791 * t905 - t792 * t901;
t826 = t856 * t905 - t857 * t901;
t806 = qJD(6) * t826 + t829 * t901 + t830 * t905;
t827 = t856 * t901 + t857 * t905;
t815 = -mrSges(7,1) * t826 + mrSges(7,2) * t827;
t866 = qJD(6) + t867;
t820 = -mrSges(7,2) * t866 + mrSges(7,3) * t826;
t840 = qJDD(6) + t842;
t786 = m(7) * t789 + mrSges(7,1) * t840 - mrSges(7,3) * t806 - t815 * t827 + t820 * t866;
t790 = t791 * t901 + t792 * t905;
t805 = -qJD(6) * t827 + t829 * t905 - t830 * t901;
t821 = mrSges(7,1) * t866 - mrSges(7,3) * t827;
t787 = m(7) * t790 - mrSges(7,2) * t840 + mrSges(7,3) * t805 + t815 * t826 - t821 * t866;
t778 = t905 * t786 + t901 * t787;
t831 = -mrSges(6,1) * t856 + mrSges(6,2) * t857;
t917 = -mrSges(6,2) * t867 + mrSges(6,3) * t856;
t776 = m(6) * t793 + t842 * mrSges(6,1) - t830 * mrSges(6,3) - t857 * t831 + t867 * t917 + t778;
t834 = mrSges(6,1) * t867 - mrSges(6,3) * t857;
t919 = -t786 * t901 + t905 * t787;
t777 = m(6) * t794 - mrSges(6,2) * t842 + mrSges(6,3) * t829 + t831 * t856 - t834 * t867 + t919;
t773 = t776 * t897 + t777 * t893;
t858 = -mrSges(5,2) * t883 - mrSges(5,3) * t867;
t859 = mrSges(5,1) * t883 - mrSges(5,3) * t868;
t911 = -m(5) * t813 - t842 * mrSges(5,1) - mrSges(5,2) * t843 - t867 * t858 - t859 * t868 - t773;
t769 = m(4) * t818 + mrSges(4,1) * t889 - mrSges(4,3) * t875 + t872 * t890 - t873 * t923 + t911;
t933 = t769 * t906;
t871 = mrSges(4,1) * t890 - mrSges(4,3) * t923;
t774 = -t776 * t893 + t897 * t777;
t850 = mrSges(5,1) * t867 + mrSges(5,2) * t868;
t772 = m(5) * t803 - mrSges(5,2) * t870 - mrSges(5,3) * t842 - t850 * t867 - t859 * t883 + t774;
t802 = -t902 * t814 + t817 * t935;
t797 = -t870 * pkin(4) - t882 * qJ(5) + t868 * t849 + qJDD(5) - t802;
t795 = -t829 * pkin(5) - t853 * pkin(11) + t857 * t835 + t797;
t912 = m(7) * t795 - t805 * mrSges(7,1) + mrSges(7,2) * t806 - t826 * t820 + t821 * t827;
t788 = m(6) * t797 - t829 * mrSges(6,1) + mrSges(6,2) * t830 + t834 * t857 - t856 * t917 + t912;
t782 = m(5) * t802 + mrSges(5,1) * t870 - mrSges(5,3) * t843 - t850 * t868 + t858 * t883 - t788;
t920 = t935 * t772 - t782 * t902;
t761 = m(4) * t819 - mrSges(4,2) * t889 - mrSges(4,3) * t876 - t871 * t890 + t873 * t922 + t920;
t764 = t902 * t772 + t935 * t782;
t832 = -t895 * t847 + t865;
t763 = m(4) * t832 + t876 * mrSges(4,1) + t875 * mrSges(4,2) + (t871 * t903 - t872 * t906) * t925 + t764;
t750 = t761 * t929 - t763 * t895 + t899 * t933;
t746 = m(3) * t854 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t908 + t750;
t749 = t761 * t932 + t899 * t763 + t895 * t933;
t748 = m(3) * t869 + t749;
t756 = t906 * t761 - t769 * t903;
t755 = m(3) * t855 - mrSges(3,1) * t908 - qJDD(2) * mrSges(3,2) + t756;
t736 = t746 * t927 - t748 * t896 + t755 * t928;
t734 = m(2) * t884 + t736;
t742 = -t746 * t904 + t907 * t755;
t741 = m(2) * t885 + t742;
t926 = t898 * t734 + t894 * t741;
t735 = t746 * t930 + t900 * t748 + t755 * t931;
t921 = -t734 * t894 + t898 * t741;
t918 = m(2) * t892 + t735;
t807 = Ifges(7,5) * t827 + Ifges(7,6) * t826 + Ifges(7,3) * t866;
t809 = Ifges(7,1) * t827 + Ifges(7,4) * t826 + Ifges(7,5) * t866;
t779 = -mrSges(7,1) * t795 + mrSges(7,3) * t790 + Ifges(7,4) * t806 + Ifges(7,2) * t805 + Ifges(7,6) * t840 - t807 * t827 + t809 * t866;
t808 = Ifges(7,4) * t827 + Ifges(7,2) * t826 + Ifges(7,6) * t866;
t780 = mrSges(7,2) * t795 - mrSges(7,3) * t789 + Ifges(7,1) * t806 + Ifges(7,4) * t805 + Ifges(7,5) * t840 + t807 * t826 - t808 * t866;
t822 = Ifges(6,5) * t857 + Ifges(6,6) * t856 + Ifges(6,3) * t867;
t824 = Ifges(6,1) * t857 + Ifges(6,4) * t856 + Ifges(6,5) * t867;
t765 = -mrSges(6,1) * t797 + mrSges(6,3) * t794 + Ifges(6,4) * t830 + Ifges(6,2) * t829 + Ifges(6,6) * t842 - pkin(5) * t912 + pkin(11) * t919 + t905 * t779 + t901 * t780 - t857 * t822 + t867 * t824;
t823 = Ifges(6,4) * t857 + Ifges(6,2) * t856 + Ifges(6,6) * t867;
t766 = mrSges(6,2) * t797 - mrSges(6,3) * t793 + Ifges(6,1) * t830 + Ifges(6,4) * t829 + Ifges(6,5) * t842 - pkin(11) * t778 - t779 * t901 + t780 * t905 + t822 * t856 - t823 * t867;
t836 = Ifges(5,5) * t868 - Ifges(5,6) * t867 + Ifges(5,3) * t883;
t837 = Ifges(5,4) * t868 - Ifges(5,2) * t867 + Ifges(5,6) * t883;
t751 = mrSges(5,2) * t813 - mrSges(5,3) * t802 + Ifges(5,1) * t843 - Ifges(5,4) * t842 + Ifges(5,5) * t870 - qJ(5) * t773 - t765 * t893 + t766 * t897 - t836 * t867 - t837 * t883;
t838 = Ifges(5,1) * t868 - Ifges(5,4) * t867 + Ifges(5,5) * t883;
t910 = mrSges(7,1) * t789 - mrSges(7,2) * t790 + Ifges(7,5) * t806 + Ifges(7,6) * t805 + Ifges(7,3) * t840 + t827 * t808 - t826 * t809;
t757 = -t910 + (-Ifges(5,2) - Ifges(6,3)) * t842 + t883 * t838 - t868 * t836 + Ifges(5,6) * t870 + t856 * t824 - t857 * t823 + Ifges(5,4) * t843 - Ifges(6,6) * t829 - Ifges(6,5) * t830 - mrSges(5,1) * t813 + mrSges(5,3) * t803 - mrSges(6,1) * t793 + mrSges(6,2) * t794 - pkin(4) * t773 - pkin(5) * t778;
t862 = Ifges(4,6) * t890 + (Ifges(4,4) * t903 + Ifges(4,2) * t906) * t925;
t863 = Ifges(4,5) * t890 + (Ifges(4,1) * t903 + Ifges(4,4) * t906) * t925;
t737 = Ifges(4,5) * t875 - Ifges(4,6) * t876 + Ifges(4,3) * t889 + mrSges(4,1) * t818 - mrSges(4,2) * t819 + t902 * t751 + t935 * t757 + pkin(3) * t911 + pkin(10) * t920 + (t862 * t903 - t863 * t906) * t925;
t861 = Ifges(4,3) * t890 + (Ifges(4,5) * t903 + Ifges(4,6) * t906) * t925;
t738 = mrSges(4,2) * t832 - mrSges(4,3) * t818 + Ifges(4,1) * t875 - Ifges(4,4) * t876 + Ifges(4,5) * t889 - pkin(10) * t764 + t751 * t935 - t902 * t757 + t861 * t922 - t890 * t862;
t909 = mrSges(5,1) * t802 - mrSges(5,2) * t803 + Ifges(5,5) * t843 - Ifges(5,6) * t842 + Ifges(5,3) * t870 - pkin(4) * t788 + qJ(5) * t774 + t897 * t765 + t893 * t766 + t868 * t837 + t867 * t838;
t743 = -mrSges(4,1) * t832 + mrSges(4,3) * t819 + Ifges(4,4) * t875 - Ifges(4,2) * t876 + Ifges(4,6) * t889 - pkin(3) * t764 - t861 * t923 + t890 * t863 - t909;
t913 = pkin(9) * t756 + t738 * t903 + t743 * t906;
t731 = -mrSges(3,1) * t869 + mrSges(3,3) * t855 + t908 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t749 - t895 * t737 + t899 * t913;
t732 = mrSges(3,2) * t869 - mrSges(3,3) * t854 + Ifges(3,5) * qJDD(2) - t908 * Ifges(3,6) + t906 * t738 - t903 * t743 + (-t749 * t895 - t750 * t899) * pkin(9);
t914 = pkin(8) * t742 + t731 * t907 + t732 * t904;
t730 = mrSges(3,1) * t854 - mrSges(3,2) * t855 + Ifges(3,3) * qJDD(2) + pkin(2) * t750 + t899 * t737 + t895 * t913;
t729 = mrSges(2,2) * t892 - mrSges(2,3) * t884 - t904 * t731 + t907 * t732 + (-t735 * t896 - t736 * t900) * pkin(8);
t728 = -mrSges(2,1) * t892 + mrSges(2,3) * t885 - pkin(1) * t735 - t896 * t730 + t900 * t914;
t1 = [-m(1) * g(1) + t921; -m(1) * g(2) + t926; -m(1) * g(3) + t918; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t926 - t894 * t728 + t898 * t729; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t921 + t898 * t728 + t894 * t729; -mrSges(1,1) * g(2) + mrSges(2,1) * t884 + mrSges(1,2) * g(1) - mrSges(2,2) * t885 + pkin(1) * t736 + t900 * t730 + t896 * t914; t918; t730; t737; t909; t788; t910;];
tauJB  = t1;
