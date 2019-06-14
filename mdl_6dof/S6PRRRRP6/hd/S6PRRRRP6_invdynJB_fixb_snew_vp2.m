% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRP6
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:12:32
% EndTime: 2019-05-05 10:12:56
% DurationCPUTime: 22.60s
% Computational Cost: add. (399671->329), mult. (818250->425), div. (0->0), fcn. (648352->14), ass. (0->150)
t937 = Ifges(6,1) + Ifges(7,1);
t930 = Ifges(6,4) - Ifges(7,5);
t929 = -Ifges(6,5) - Ifges(7,4);
t936 = Ifges(6,2) + Ifges(7,3);
t928 = Ifges(6,6) - Ifges(7,6);
t935 = -Ifges(6,3) - Ifges(7,2);
t886 = sin(pkin(7));
t893 = sin(qJ(3));
t896 = cos(qJ(3));
t913 = qJD(2) * qJD(3);
t869 = (-qJDD(2) * t896 + t893 * t913) * t886;
t885 = sin(pkin(12));
t888 = cos(pkin(12));
t876 = g(1) * t885 - g(2) * t888;
t877 = -g(1) * t888 - g(2) * t885;
t884 = -g(3) + qJDD(1);
t894 = sin(qJ(2));
t890 = cos(pkin(6));
t897 = cos(qJ(2));
t920 = t890 * t897;
t887 = sin(pkin(6));
t923 = t887 * t897;
t845 = t876 * t920 - t877 * t894 + t884 * t923;
t898 = qJD(2) ^ 2;
t932 = pkin(9) * t886;
t841 = qJDD(2) * pkin(2) + t898 * t932 + t845;
t921 = t890 * t894;
t924 = t887 * t894;
t846 = t876 * t921 + t897 * t877 + t884 * t924;
t842 = -pkin(2) * t898 + qJDD(2) * t932 + t846;
t862 = -t876 * t887 + t884 * t890;
t889 = cos(pkin(7));
t803 = -t893 * t842 + (t841 * t889 + t862 * t886) * t896;
t882 = qJD(2) * t889 + qJD(3);
t892 = sin(qJ(4));
t895 = cos(qJ(4));
t914 = qJD(2) * t886;
t911 = t893 * t914;
t860 = t882 * t895 - t892 * t911;
t868 = (qJDD(2) * t893 + t896 * t913) * t886;
t881 = qJDD(2) * t889 + qJDD(3);
t837 = qJD(4) * t860 + t868 * t895 + t881 * t892;
t861 = t882 * t892 + t895 * t911;
t910 = t896 * t914;
t875 = qJD(4) - t910;
t891 = sin(qJ(5));
t933 = cos(qJ(5));
t847 = t891 * t861 - t875 * t933;
t863 = qJDD(4) + t869;
t808 = -t847 * qJD(5) + t837 * t933 + t891 * t863;
t848 = t861 * t933 + t891 * t875;
t820 = mrSges(7,1) * t847 - mrSges(7,3) * t848;
t922 = t889 * t893;
t925 = t886 * t893;
t804 = t841 * t922 + t896 * t842 + t862 * t925;
t867 = (-pkin(3) * t896 - pkin(10) * t893) * t914;
t880 = t882 ^ 2;
t800 = -pkin(3) * t880 + pkin(10) * t881 + t867 * t910 + t804;
t856 = t889 * t862;
t802 = t869 * pkin(3) - t868 * pkin(10) + t856 + (-t841 + (pkin(3) * t893 - pkin(10) * t896) * t882 * qJD(2)) * t886;
t796 = t895 * t800 + t892 * t802;
t844 = -pkin(4) * t860 - pkin(11) * t861;
t874 = t875 ^ 2;
t792 = -pkin(4) * t874 + pkin(11) * t863 + t844 * t860 + t796;
t799 = -t881 * pkin(3) - t880 * pkin(10) + t867 * t911 - t803;
t836 = -qJD(4) * t861 - t868 * t892 + t881 * t895;
t794 = (-t860 * t875 - t837) * pkin(11) + (t861 * t875 - t836) * pkin(4) + t799;
t788 = -t891 * t792 + t794 * t933;
t819 = pkin(5) * t847 - qJ(6) * t848;
t834 = qJDD(5) - t836;
t858 = qJD(5) - t860;
t857 = t858 ^ 2;
t786 = -t834 * pkin(5) - t857 * qJ(6) + t848 * t819 + qJDD(6) - t788;
t824 = -mrSges(7,2) * t847 + mrSges(7,3) * t858;
t906 = -m(7) * t786 + t834 * mrSges(7,1) + t858 * t824;
t782 = t808 * mrSges(7,2) + t848 * t820 - t906;
t789 = t792 * t933 + t891 * t794;
t785 = -pkin(5) * t857 + qJ(6) * t834 + 0.2e1 * qJD(6) * t858 - t819 * t847 + t789;
t807 = t848 * qJD(5) + t891 * t837 - t863 * t933;
t827 = -mrSges(7,1) * t858 + mrSges(7,2) * t848;
t912 = m(7) * t785 + t834 * mrSges(7,3) + t858 * t827;
t916 = t930 * t847 - t937 * t848 + t929 * t858;
t917 = t936 * t847 - t930 * t848 - t928 * t858;
t934 = -t807 * t928 - t808 * t929 - t935 * t834 - t847 * t916 - t848 * t917 + mrSges(6,1) * t788 - mrSges(7,1) * t786 - mrSges(6,2) * t789 + mrSges(7,3) * t785 - pkin(5) * t782 + qJ(6) * (-t807 * mrSges(7,2) - t847 * t820 + t912);
t931 = -mrSges(6,3) - mrSges(7,2);
t865 = -mrSges(4,2) * t882 + mrSges(4,3) * t910;
t866 = (-mrSges(4,1) * t896 + mrSges(4,2) * t893) * t914;
t826 = mrSges(6,1) * t858 - mrSges(6,3) * t848;
t915 = -mrSges(6,1) * t847 - mrSges(6,2) * t848 - t820;
t778 = m(6) * t789 - t834 * mrSges(6,2) + t807 * t931 - t858 * t826 + t847 * t915 + t912;
t825 = -mrSges(6,2) * t858 - mrSges(6,3) * t847;
t779 = m(6) * t788 + t834 * mrSges(6,1) + t808 * t931 + t858 * t825 + t848 * t915 + t906;
t773 = t891 * t778 + t779 * t933;
t849 = -mrSges(5,2) * t875 + mrSges(5,3) * t860;
t850 = mrSges(5,1) * t875 - mrSges(5,3) * t861;
t900 = -m(5) * t799 + t836 * mrSges(5,1) - t837 * mrSges(5,2) + t860 * t849 - t861 * t850 - t773;
t767 = m(4) * t803 + t881 * mrSges(4,1) - t868 * mrSges(4,3) + t882 * t865 - t866 * t911 + t900;
t926 = t767 * t896;
t864 = mrSges(4,1) * t882 - mrSges(4,3) * t911;
t774 = t778 * t933 - t779 * t891;
t843 = -mrSges(5,1) * t860 + mrSges(5,2) * t861;
t770 = m(5) * t796 - mrSges(5,2) * t863 + mrSges(5,3) * t836 + t843 * t860 - t850 * t875 + t774;
t795 = -t892 * t800 + t895 * t802;
t791 = -t863 * pkin(4) - t874 * pkin(11) + t861 * t844 - t795;
t787 = -0.2e1 * qJD(6) * t848 + (t847 * t858 - t808) * qJ(6) + (t848 * t858 + t807) * pkin(5) + t791;
t783 = m(7) * t787 + mrSges(7,1) * t807 - t808 * mrSges(7,3) + t824 * t847 - t848 * t827;
t780 = -m(6) * t791 - t807 * mrSges(6,1) - mrSges(6,2) * t808 - t847 * t825 - t826 * t848 - t783;
t776 = m(5) * t795 + mrSges(5,1) * t863 - mrSges(5,3) * t837 - t843 * t861 + t849 * t875 + t780;
t908 = t895 * t770 - t776 * t892;
t761 = m(4) * t804 - mrSges(4,2) * t881 - mrSges(4,3) * t869 - t864 * t882 + t866 * t910 + t908;
t764 = t892 * t770 + t895 * t776;
t822 = -t886 * t841 + t856;
t763 = m(4) * t822 + t869 * mrSges(4,1) + t868 * mrSges(4,2) + (t864 * t893 - t865 * t896) * t914 + t764;
t750 = t761 * t922 - t763 * t886 + t889 * t926;
t746 = m(3) * t845 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t898 + t750;
t749 = t761 * t925 + t889 * t763 + t886 * t926;
t748 = m(3) * t862 + t749;
t756 = t896 * t761 - t767 * t893;
t755 = m(3) * t846 - mrSges(3,1) * t898 - qJDD(2) * mrSges(3,2) + t756;
t736 = t746 * t920 - t748 * t887 + t755 * t921;
t734 = m(2) * t876 + t736;
t742 = -t746 * t894 + t897 * t755;
t741 = m(2) * t877 + t742;
t919 = t888 * t734 + t885 * t741;
t918 = t928 * t847 + t929 * t848 + t935 * t858;
t735 = t746 * t923 + t890 * t748 + t755 * t924;
t909 = -t734 * t885 + t888 * t741;
t907 = m(2) * t884 + t735;
t771 = -mrSges(6,1) * t791 - mrSges(7,1) * t787 + mrSges(7,2) * t785 + mrSges(6,3) * t789 - pkin(5) * t783 - t936 * t807 + t930 * t808 + t928 * t834 + t918 * t848 - t916 * t858;
t772 = mrSges(6,2) * t791 + mrSges(7,2) * t786 - mrSges(6,3) * t788 - mrSges(7,3) * t787 - qJ(6) * t783 - t930 * t807 + t937 * t808 - t929 * t834 + t918 * t847 + t917 * t858;
t830 = Ifges(5,5) * t861 + Ifges(5,6) * t860 + Ifges(5,3) * t875;
t831 = Ifges(5,4) * t861 + Ifges(5,2) * t860 + Ifges(5,6) * t875;
t751 = mrSges(5,2) * t799 - mrSges(5,3) * t795 + Ifges(5,1) * t837 + Ifges(5,4) * t836 + Ifges(5,5) * t863 - pkin(11) * t773 - t891 * t771 + t772 * t933 + t860 * t830 - t875 * t831;
t832 = Ifges(5,1) * t861 + Ifges(5,4) * t860 + Ifges(5,5) * t875;
t757 = -mrSges(5,1) * t799 + mrSges(5,3) * t796 + Ifges(5,4) * t837 + Ifges(5,2) * t836 + Ifges(5,6) * t863 - pkin(4) * t773 - t861 * t830 + t875 * t832 - t934;
t853 = Ifges(4,6) * t882 + (Ifges(4,4) * t893 + Ifges(4,2) * t896) * t914;
t854 = Ifges(4,5) * t882 + (Ifges(4,1) * t893 + Ifges(4,4) * t896) * t914;
t737 = Ifges(4,5) * t868 - Ifges(4,6) * t869 + Ifges(4,3) * t881 + mrSges(4,1) * t803 - mrSges(4,2) * t804 + t892 * t751 + t895 * t757 + pkin(3) * t900 + pkin(10) * t908 + (t853 * t893 - t854 * t896) * t914;
t852 = Ifges(4,3) * t882 + (Ifges(4,5) * t893 + Ifges(4,6) * t896) * t914;
t738 = mrSges(4,2) * t822 - mrSges(4,3) * t803 + Ifges(4,1) * t868 - Ifges(4,4) * t869 + Ifges(4,5) * t881 - pkin(10) * t764 + t751 * t895 - t757 * t892 + t852 * t910 - t853 * t882;
t899 = mrSges(5,1) * t795 - mrSges(5,2) * t796 + Ifges(5,5) * t837 + Ifges(5,6) * t836 + Ifges(5,3) * t863 + pkin(4) * t780 + pkin(11) * t774 + t771 * t933 + t891 * t772 + t861 * t831 - t860 * t832;
t743 = -mrSges(4,1) * t822 + mrSges(4,3) * t804 + Ifges(4,4) * t868 - Ifges(4,2) * t869 + Ifges(4,6) * t881 - pkin(3) * t764 - t852 * t911 + t882 * t854 - t899;
t902 = pkin(9) * t756 + t738 * t893 + t743 * t896;
t731 = -mrSges(3,1) * t862 + mrSges(3,3) * t846 + t898 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t749 - t886 * t737 + t889 * t902;
t732 = mrSges(3,2) * t862 - mrSges(3,3) * t845 + Ifges(3,5) * qJDD(2) - t898 * Ifges(3,6) + t896 * t738 - t893 * t743 + (-t749 * t886 - t750 * t889) * pkin(9);
t903 = pkin(8) * t742 + t731 * t897 + t732 * t894;
t730 = mrSges(3,1) * t845 - mrSges(3,2) * t846 + Ifges(3,3) * qJDD(2) + pkin(2) * t750 + t889 * t737 + t886 * t902;
t729 = mrSges(2,2) * t884 - mrSges(2,3) * t876 - t894 * t731 + t897 * t732 + (-t735 * t887 - t736 * t890) * pkin(8);
t728 = -mrSges(2,1) * t884 + mrSges(2,3) * t877 - pkin(1) * t735 - t887 * t730 + t890 * t903;
t1 = [-m(1) * g(1) + t909; -m(1) * g(2) + t919; -m(1) * g(3) + t907; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t919 - t885 * t728 + t888 * t729; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t909 + t888 * t728 + t885 * t729; -mrSges(1,1) * g(2) + mrSges(2,1) * t876 + mrSges(1,2) * g(1) - mrSges(2,2) * t877 + pkin(1) * t736 + t890 * t730 + t887 * t903; t907; t730; t737; t899; t934; t782;];
tauJB  = t1;
