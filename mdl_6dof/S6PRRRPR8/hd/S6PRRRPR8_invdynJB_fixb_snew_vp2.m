% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 09:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:06:37
% EndTime: 2019-05-05 09:06:55
% DurationCPUTime: 17.84s
% Computational Cost: add. (309147->335), mult. (639330->424), div. (0->0), fcn. (497462->14), ass. (0->154)
t974 = Ifges(5,1) + Ifges(6,2);
t967 = Ifges(5,4) + Ifges(6,6);
t966 = Ifges(5,5) - Ifges(6,4);
t973 = -Ifges(5,2) - Ifges(6,3);
t965 = Ifges(5,6) - Ifges(6,5);
t972 = Ifges(5,3) + Ifges(6,1);
t918 = sin(pkin(7));
t925 = sin(qJ(3));
t928 = cos(qJ(3));
t949 = qJD(2) * qJD(3);
t899 = (-qJDD(2) * t928 + t925 * t949) * t918;
t917 = sin(pkin(12));
t920 = cos(pkin(12));
t908 = g(1) * t917 - g(2) * t920;
t909 = -g(1) * t920 - g(2) * t917;
t916 = -g(3) + qJDD(1);
t926 = sin(qJ(2));
t922 = cos(pkin(6));
t929 = cos(qJ(2));
t956 = t922 * t929;
t919 = sin(pkin(6));
t959 = t919 * t929;
t870 = t908 * t956 - t909 * t926 + t916 * t959;
t930 = qJD(2) ^ 2;
t968 = pkin(9) * t918;
t864 = qJDD(2) * pkin(2) + t930 * t968 + t870;
t957 = t922 * t926;
t960 = t919 * t926;
t871 = t908 * t957 + t929 * t909 + t916 * t960;
t865 = -pkin(2) * t930 + qJDD(2) * t968 + t871;
t891 = -t908 * t919 + t916 * t922;
t921 = cos(pkin(7));
t829 = -t925 * t865 + (t864 * t921 + t891 * t918) * t928;
t914 = qJD(2) * t921 + qJD(3);
t924 = sin(qJ(4));
t950 = qJD(2) * t918;
t948 = t925 * t950;
t969 = cos(qJ(4));
t889 = -t969 * t914 + t924 * t948;
t947 = t928 * t950;
t906 = -qJD(4) + t947;
t876 = mrSges(6,1) * t889 + mrSges(6,3) * t906;
t893 = qJDD(4) + t899;
t958 = t921 * t925;
t961 = t918 * t925;
t830 = t864 * t958 + t928 * t865 + t891 * t961;
t897 = (-pkin(3) * t928 - pkin(10) * t925) * t950;
t912 = t914 ^ 2;
t913 = qJDD(2) * t921 + qJDD(3);
t825 = -pkin(3) * t912 + pkin(10) * t913 + t897 * t947 + t830;
t885 = t921 * t891;
t898 = (qJDD(2) * t925 + t928 * t949) * t918;
t827 = t899 * pkin(3) - t898 * pkin(10) + t885 + (-t864 + (pkin(3) * t925 - pkin(10) * t928) * t914 * qJD(2)) * t918;
t820 = -t924 * t825 + t969 * t827;
t890 = t924 * t914 + t969 * t948;
t866 = pkin(4) * t889 - qJ(5) * t890;
t905 = t906 ^ 2;
t818 = -t893 * pkin(4) - t905 * qJ(5) + t890 * t866 + qJDD(5) - t820;
t860 = -t889 * qJD(4) + t969 * t898 + t924 * t913;
t962 = t889 * t906;
t813 = (t889 * t890 - t893) * pkin(11) + (t860 - t962) * pkin(5) + t818;
t859 = qJD(4) * t890 + t898 * t924 - t969 * t913;
t878 = pkin(5) * t890 + pkin(11) * t906;
t888 = t889 ^ 2;
t824 = -t913 * pkin(3) - t912 * pkin(10) + t897 * t948 - t829;
t970 = -2 * qJD(5);
t932 = (-t860 - t962) * qJ(5) + t824 + (-t906 * pkin(4) + t970) * t890;
t816 = -t888 * pkin(5) - t890 * t878 + (pkin(4) + pkin(11)) * t859 + t932;
t923 = sin(qJ(6));
t927 = cos(qJ(6));
t811 = t813 * t927 - t816 * t923;
t872 = t889 * t927 + t906 * t923;
t835 = qJD(6) * t872 + t859 * t923 + t893 * t927;
t873 = t889 * t923 - t906 * t927;
t840 = -mrSges(7,1) * t872 + mrSges(7,2) * t873;
t887 = qJD(6) + t890;
t844 = -mrSges(7,2) * t887 + mrSges(7,3) * t872;
t856 = qJDD(6) + t860;
t808 = m(7) * t811 + mrSges(7,1) * t856 - mrSges(7,3) * t835 - t840 * t873 + t844 * t887;
t812 = t813 * t923 + t816 * t927;
t834 = -qJD(6) * t873 + t859 * t927 - t893 * t923;
t845 = mrSges(7,1) * t887 - mrSges(7,3) * t873;
t809 = m(7) * t812 - mrSges(7,2) * t856 + mrSges(7,3) * t834 + t840 * t872 - t845 * t887;
t799 = t927 * t808 + t923 * t809;
t868 = -mrSges(6,2) * t889 - mrSges(6,3) * t890;
t937 = -m(6) * t818 - t860 * mrSges(6,1) - t890 * t868 - t799;
t797 = t893 * mrSges(6,2) - t906 * t876 - t937;
t821 = t969 * t825 + t924 * t827;
t936 = -t905 * pkin(4) + t893 * qJ(5) - t889 * t866 + t821;
t815 = -t859 * pkin(5) - t888 * pkin(11) + (t970 - t878) * t906 + t936;
t836 = Ifges(7,5) * t873 + Ifges(7,6) * t872 + Ifges(7,3) * t887;
t838 = Ifges(7,1) * t873 + Ifges(7,4) * t872 + Ifges(7,5) * t887;
t800 = -mrSges(7,1) * t815 + mrSges(7,3) * t812 + Ifges(7,4) * t835 + Ifges(7,2) * t834 + Ifges(7,6) * t856 - t836 * t873 + t838 * t887;
t837 = Ifges(7,4) * t873 + Ifges(7,2) * t872 + Ifges(7,6) * t887;
t801 = mrSges(7,2) * t815 - mrSges(7,3) * t811 + Ifges(7,1) * t835 + Ifges(7,4) * t834 + Ifges(7,5) * t856 + t836 * t872 - t837 * t887;
t817 = 0.2e1 * qJD(5) * t906 - t936;
t877 = mrSges(6,1) * t890 - mrSges(6,2) * t906;
t938 = -m(7) * t815 + t834 * mrSges(7,1) - t835 * mrSges(7,2) + t872 * t844 - t873 * t845;
t934 = -m(6) * t817 + t893 * mrSges(6,3) - t906 * t877 - t938;
t951 = -t967 * t889 + t974 * t890 - t966 * t906;
t952 = t973 * t889 + t967 * t890 - t965 * t906;
t971 = -t965 * t859 + t966 * t860 + t951 * t889 + t952 * t890 + t972 * t893 + mrSges(5,1) * t820 - mrSges(5,2) * t821 + mrSges(6,2) * t818 - mrSges(6,3) * t817 - pkin(4) * t797 - pkin(11) * t799 + qJ(5) * (-t859 * mrSges(6,1) - t889 * t868 + t934) - t923 * t800 + t927 * t801;
t895 = -mrSges(4,2) * t914 + mrSges(4,3) * t947;
t896 = (-mrSges(4,1) * t928 + mrSges(4,2) * t925) * t950;
t874 = mrSges(5,2) * t906 - mrSges(5,3) * t889;
t875 = -mrSges(5,1) * t906 - mrSges(5,3) * t890;
t819 = t859 * pkin(4) + t932;
t954 = -t923 * t808 + t927 * t809;
t942 = -m(6) * t819 + t859 * mrSges(6,2) + t889 * t876 - t954;
t933 = -m(5) * t824 - t859 * mrSges(5,1) - t889 * t874 + (-t875 + t877) * t890 + (-mrSges(5,2) + mrSges(6,3)) * t860 + t942;
t794 = m(4) * t829 + t913 * mrSges(4,1) - t898 * mrSges(4,3) + t914 * t895 - t896 * t948 + t933;
t963 = t794 * t928;
t894 = mrSges(4,1) * t914 - mrSges(4,3) * t948;
t867 = mrSges(5,1) * t889 + mrSges(5,2) * t890;
t796 = m(5) * t820 - t860 * mrSges(5,3) - t890 * t867 + (-t874 + t876) * t906 + (mrSges(5,1) - mrSges(6,2)) * t893 + t937;
t804 = m(5) * t821 - t893 * mrSges(5,2) + t906 * t875 + (-t867 - t868) * t889 + (-mrSges(5,3) - mrSges(6,1)) * t859 + t934;
t945 = -t796 * t924 + t969 * t804;
t788 = m(4) * t830 - mrSges(4,2) * t913 - mrSges(4,3) * t899 - t894 * t914 + t896 * t947 + t945;
t791 = t969 * t796 + t924 * t804;
t841 = -t918 * t864 + t885;
t790 = m(4) * t841 + t899 * mrSges(4,1) + t898 * mrSges(4,2) + (t894 * t925 - t895 * t928) * t950 + t791;
t777 = t788 * t958 - t790 * t918 + t921 * t963;
t773 = m(3) * t870 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t930 + t777;
t776 = t788 * t961 + t921 * t790 + t918 * t963;
t775 = m(3) * t891 + t776;
t784 = t928 * t788 - t794 * t925;
t783 = m(3) * t871 - mrSges(3,1) * t930 - qJDD(2) * mrSges(3,2) + t784;
t763 = t773 * t956 - t775 * t919 + t783 * t957;
t761 = m(2) * t908 + t763;
t769 = -t773 * t926 + t929 * t783;
t768 = m(2) * t909 + t769;
t955 = t920 * t761 + t917 * t768;
t953 = t965 * t889 - t966 * t890 + t972 * t906;
t762 = t773 * t959 + t922 * t775 + t783 * t960;
t946 = -t761 * t917 + t920 * t768;
t944 = m(2) * t916 + t762;
t798 = -t860 * mrSges(6,3) - t890 * t877 - t942;
t778 = -mrSges(5,1) * t824 - mrSges(6,1) * t817 + mrSges(6,2) * t819 + mrSges(5,3) * t821 - pkin(4) * t798 - pkin(5) * t938 - pkin(11) * t954 - t927 * t800 - t923 * t801 + t973 * t859 + t967 * t860 + t953 * t890 + t965 * t893 - t951 * t906;
t935 = mrSges(7,1) * t811 - mrSges(7,2) * t812 + Ifges(7,5) * t835 + Ifges(7,6) * t834 + Ifges(7,3) * t856 + t873 * t837 - t872 * t838;
t779 = mrSges(6,1) * t818 + mrSges(5,2) * t824 - mrSges(5,3) * t820 - mrSges(6,3) * t819 + pkin(5) * t799 - qJ(5) * t798 - t967 * t859 + t974 * t860 + t953 * t889 + t966 * t893 + t952 * t906 + t935;
t882 = Ifges(4,6) * t914 + (Ifges(4,4) * t925 + Ifges(4,2) * t928) * t950;
t883 = Ifges(4,5) * t914 + (Ifges(4,1) * t925 + Ifges(4,4) * t928) * t950;
t764 = Ifges(4,5) * t898 - Ifges(4,6) * t899 + Ifges(4,3) * t913 + mrSges(4,1) * t829 - mrSges(4,2) * t830 + t924 * t779 + t969 * t778 + pkin(3) * t933 + pkin(10) * t945 + (t882 * t925 - t883 * t928) * t950;
t881 = Ifges(4,3) * t914 + (Ifges(4,5) * t925 + Ifges(4,6) * t928) * t950;
t765 = mrSges(4,2) * t841 - mrSges(4,3) * t829 + Ifges(4,1) * t898 - Ifges(4,4) * t899 + Ifges(4,5) * t913 - pkin(10) * t791 - t924 * t778 + t969 * t779 + t881 * t947 - t914 * t882;
t770 = -mrSges(4,1) * t841 + mrSges(4,3) * t830 + Ifges(4,4) * t898 - Ifges(4,2) * t899 + Ifges(4,6) * t913 - pkin(3) * t791 - t881 * t948 + t914 * t883 - t971;
t939 = pkin(9) * t784 + t765 * t925 + t770 * t928;
t758 = -mrSges(3,1) * t891 + mrSges(3,3) * t871 + t930 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t776 - t918 * t764 + t939 * t921;
t759 = mrSges(3,2) * t891 - mrSges(3,3) * t870 + Ifges(3,5) * qJDD(2) - t930 * Ifges(3,6) + t928 * t765 - t925 * t770 + (-t776 * t918 - t777 * t921) * pkin(9);
t940 = pkin(8) * t769 + t758 * t929 + t759 * t926;
t757 = mrSges(3,1) * t870 - mrSges(3,2) * t871 + Ifges(3,3) * qJDD(2) + pkin(2) * t777 + t921 * t764 + t939 * t918;
t756 = mrSges(2,2) * t916 - mrSges(2,3) * t908 - t926 * t758 + t929 * t759 + (-t762 * t919 - t763 * t922) * pkin(8);
t755 = -mrSges(2,1) * t916 + mrSges(2,3) * t909 - pkin(1) * t762 - t919 * t757 + t940 * t922;
t1 = [-m(1) * g(1) + t946; -m(1) * g(2) + t955; -m(1) * g(3) + t944; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t955 - t917 * t755 + t920 * t756; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t946 + t920 * t755 + t917 * t756; -mrSges(1,1) * g(2) + mrSges(2,1) * t908 + mrSges(1,2) * g(1) - mrSges(2,2) * t909 + pkin(1) * t763 + t922 * t757 + t940 * t919; t944; t757; t764; t971; t797; t935;];
tauJB  = t1;
