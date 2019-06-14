% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:17:52
% EndTime: 2019-05-06 17:18:12
% DurationCPUTime: 17.22s
% Computational Cost: add. (243240->366), mult. (563066->448), div. (0->0), fcn. (414827->10), ass. (0->145)
t974 = Ifges(6,1) + Ifges(7,1);
t967 = Ifges(6,4) + Ifges(7,4);
t966 = Ifges(6,5) + Ifges(7,5);
t973 = Ifges(6,2) + Ifges(7,2);
t965 = Ifges(6,6) + Ifges(7,6);
t972 = Ifges(6,3) + Ifges(7,3);
t930 = sin(qJ(2));
t934 = cos(qJ(2));
t954 = qJD(1) * qJD(2);
t911 = qJDD(1) * t930 + t934 * t954;
t931 = sin(qJ(1));
t935 = cos(qJ(1));
t918 = -g(1) * t935 - g(2) * t931;
t936 = qJD(1) ^ 2;
t906 = -pkin(1) * t936 + qJDD(1) * pkin(7) + t918;
t962 = t930 * t906;
t969 = pkin(2) * t936;
t871 = qJDD(2) * pkin(2) - t911 * qJ(3) - t962 + (qJ(3) * t954 + t930 * t969 - g(3)) * t934;
t892 = -g(3) * t930 + t934 * t906;
t912 = qJDD(1) * t934 - t930 * t954;
t956 = qJD(1) * t930;
t914 = qJD(2) * pkin(2) - qJ(3) * t956;
t925 = t934 ^ 2;
t872 = qJ(3) * t912 - qJD(2) * t914 - t925 * t969 + t892;
t926 = sin(pkin(10));
t927 = cos(pkin(10));
t901 = (t926 * t934 + t927 * t930) * qJD(1);
t842 = -0.2e1 * qJD(3) * t901 + t927 * t871 - t926 * t872;
t890 = t911 * t927 + t912 * t926;
t900 = (-t926 * t930 + t927 * t934) * qJD(1);
t821 = (qJD(2) * t900 - t890) * pkin(8) + (t900 * t901 + qJDD(2)) * pkin(3) + t842;
t843 = 0.2e1 * qJD(3) * t900 + t926 * t871 + t927 * t872;
t889 = -t911 * t926 + t912 * t927;
t895 = qJD(2) * pkin(3) - pkin(8) * t901;
t899 = t900 ^ 2;
t824 = -pkin(3) * t899 + pkin(8) * t889 - qJD(2) * t895 + t843;
t929 = sin(qJ(4));
t933 = cos(qJ(4));
t819 = t929 * t821 + t933 * t824;
t884 = t900 * t929 + t901 * t933;
t852 = -qJD(4) * t884 + t889 * t933 - t890 * t929;
t883 = t900 * t933 - t901 * t929;
t866 = -mrSges(5,1) * t883 + mrSges(5,2) * t884;
t923 = qJD(2) + qJD(4);
t878 = mrSges(5,1) * t923 - mrSges(5,3) * t884;
t922 = qJDD(2) + qJDD(4);
t867 = -pkin(4) * t883 - pkin(9) * t884;
t921 = t923 ^ 2;
t813 = -pkin(4) * t921 + pkin(9) * t922 + t867 * t883 + t819;
t917 = t931 * g(1) - t935 * g(2);
t944 = -qJDD(1) * pkin(1) - t917;
t874 = -t912 * pkin(2) + qJDD(3) + t914 * t956 + (-qJ(3) * t925 - pkin(7)) * t936 + t944;
t840 = -t889 * pkin(3) - t899 * pkin(8) + t901 * t895 + t874;
t853 = qJD(4) * t883 + t889 * t929 + t890 * t933;
t816 = (-t883 * t923 - t853) * pkin(9) + (t884 * t923 - t852) * pkin(4) + t840;
t928 = sin(qJ(5));
t932 = cos(qJ(5));
t808 = -t928 * t813 + t932 * t816;
t875 = -t884 * t928 + t923 * t932;
t829 = qJD(5) * t875 + t853 * t932 + t922 * t928;
t851 = qJDD(5) - t852;
t876 = t884 * t932 + t923 * t928;
t854 = -mrSges(7,1) * t875 + mrSges(7,2) * t876;
t855 = -mrSges(6,1) * t875 + mrSges(6,2) * t876;
t879 = qJD(5) - t883;
t857 = -mrSges(6,2) * t879 + mrSges(6,3) * t875;
t805 = -0.2e1 * qJD(6) * t876 + (t875 * t879 - t829) * qJ(6) + (t875 * t876 + t851) * pkin(5) + t808;
t856 = -mrSges(7,2) * t879 + mrSges(7,3) * t875;
t953 = m(7) * t805 + t851 * mrSges(7,1) + t879 * t856;
t794 = m(6) * t808 + t851 * mrSges(6,1) + t879 * t857 + (-t854 - t855) * t876 + (-mrSges(6,3) - mrSges(7,3)) * t829 + t953;
t809 = t932 * t813 + t928 * t816;
t828 = -qJD(5) * t876 - t853 * t928 + t922 * t932;
t858 = pkin(5) * t879 - qJ(6) * t876;
t873 = t875 ^ 2;
t807 = -pkin(5) * t873 + t828 * qJ(6) + 0.2e1 * qJD(6) * t875 - t858 * t879 + t809;
t952 = m(7) * t807 + t828 * mrSges(7,3) + t875 * t854;
t859 = mrSges(7,1) * t879 - mrSges(7,3) * t876;
t957 = -mrSges(6,1) * t879 + mrSges(6,3) * t876 - t859;
t968 = -mrSges(6,2) - mrSges(7,2);
t797 = m(6) * t809 + t828 * mrSges(6,3) + t851 * t968 + t875 * t855 + t879 * t957 + t952;
t947 = -t794 * t928 + t932 * t797;
t786 = m(5) * t819 - mrSges(5,2) * t922 + mrSges(5,3) * t852 + t866 * t883 - t878 * t923 + t947;
t818 = t821 * t933 - t929 * t824;
t877 = -mrSges(5,2) * t923 + mrSges(5,3) * t883;
t812 = -pkin(4) * t922 - pkin(9) * t921 + t884 * t867 - t818;
t810 = -t828 * pkin(5) - qJ(6) * t873 + t858 * t876 + qJDD(6) + t812;
t946 = -m(7) * t810 + t828 * mrSges(7,1) + t875 * t856;
t939 = -m(6) * t812 + t828 * mrSges(6,1) + t829 * t968 + t875 * t857 + t876 * t957 + t946;
t799 = m(5) * t818 + t922 * mrSges(5,1) - t853 * mrSges(5,3) - t884 * t866 + t923 * t877 + t939;
t778 = t929 * t786 + t933 * t799;
t887 = -mrSges(4,1) * t900 + mrSges(4,2) * t901;
t893 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t900;
t776 = m(4) * t842 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t890 + qJD(2) * t893 - t887 * t901 + t778;
t894 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t901;
t948 = t933 * t786 - t799 * t929;
t777 = m(4) * t843 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t889 - qJD(2) * t894 + t887 * t900 + t948;
t771 = t927 * t776 + t926 * t777;
t881 = Ifges(4,4) * t901 + Ifges(4,2) * t900 + Ifges(4,6) * qJD(2);
t882 = Ifges(4,1) * t901 + Ifges(4,4) * t900 + Ifges(4,5) * qJD(2);
t891 = -t934 * g(3) - t962;
t903 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t930 + Ifges(3,2) * t934) * qJD(1);
t904 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t930 + Ifges(3,4) * t934) * qJD(1);
t803 = t829 * mrSges(7,2) + t876 * t859 - t946;
t958 = t967 * t875 + t974 * t876 + t966 * t879;
t960 = -t965 * t875 - t966 * t876 - t972 * t879;
t780 = -mrSges(6,1) * t812 + mrSges(6,3) * t809 - mrSges(7,1) * t810 + mrSges(7,3) * t807 - pkin(5) * t803 + qJ(6) * t952 + (-qJ(6) * t859 + t958) * t879 + t960 * t876 + (-mrSges(7,2) * qJ(6) + t965) * t851 + t967 * t829 + t973 * t828;
t802 = -t829 * mrSges(7,3) - t876 * t854 + t953;
t959 = -t973 * t875 - t967 * t876 - t965 * t879;
t789 = mrSges(6,2) * t812 + mrSges(7,2) * t810 - mrSges(6,3) * t808 - mrSges(7,3) * t805 - qJ(6) * t802 + t967 * t828 + t974 * t829 + t966 * t851 - t960 * t875 + t959 * t879;
t862 = Ifges(5,4) * t884 + Ifges(5,2) * t883 + Ifges(5,6) * t923;
t863 = Ifges(5,1) * t884 + Ifges(5,4) * t883 + Ifges(5,5) * t923;
t940 = -mrSges(5,1) * t818 + mrSges(5,2) * t819 - Ifges(5,5) * t853 - Ifges(5,6) * t852 - Ifges(5,3) * t922 - pkin(4) * t939 - pkin(9) * t947 - t932 * t780 - t928 * t789 - t884 * t862 + t883 * t863;
t971 = mrSges(3,1) * t891 + mrSges(4,1) * t842 - mrSges(3,2) * t892 - mrSges(4,2) * t843 + Ifges(3,5) * t911 + Ifges(4,5) * t890 + Ifges(3,6) * t912 + Ifges(4,6) * t889 + pkin(2) * t771 + pkin(3) * t778 + (t903 * t930 - t904 * t934) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t901 * t881 - t900 * t882 - t940;
t970 = mrSges(6,1) * t808 + mrSges(7,1) * t805 - mrSges(6,2) * t809 - mrSges(7,2) * t807 + pkin(5) * t802 + t828 * t965 + t829 * t966 + t972 * t851 - t875 * t958 - t876 * t959;
t910 = (-mrSges(3,1) * t934 + mrSges(3,2) * t930) * qJD(1);
t955 = qJD(1) * t934;
t916 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t955;
t769 = m(3) * t891 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t911 + qJD(2) * t916 - t910 * t956 + t771;
t915 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t956;
t949 = -t776 * t926 + t927 * t777;
t770 = m(3) * t892 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t912 - qJD(2) * t915 + t910 * t955 + t949;
t950 = -t769 * t930 + t934 * t770;
t762 = m(2) * t918 - mrSges(2,1) * t936 - qJDD(1) * mrSges(2,2) + t950;
t791 = t932 * t794 + t928 * t797;
t943 = m(5) * t840 - t852 * mrSges(5,1) + t853 * mrSges(5,2) - t883 * t877 + t884 * t878 + t791;
t787 = m(4) * t874 - t889 * mrSges(4,1) + mrSges(4,2) * t890 - t900 * t893 + t894 * t901 + t943;
t905 = -t936 * pkin(7) + t944;
t938 = -m(3) * t905 + t912 * mrSges(3,1) - mrSges(3,2) * t911 - t915 * t956 + t916 * t955 - t787;
t782 = m(2) * t917 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t936 + t938;
t961 = t931 * t762 + t935 * t782;
t764 = t934 * t769 + t930 * t770;
t951 = t935 * t762 - t782 * t931;
t861 = Ifges(5,5) * t884 + Ifges(5,6) * t883 + Ifges(5,3) * t923;
t765 = mrSges(5,2) * t840 - mrSges(5,3) * t818 + Ifges(5,1) * t853 + Ifges(5,4) * t852 + Ifges(5,5) * t922 - pkin(9) * t791 - t780 * t928 + t789 * t932 + t861 * t883 - t862 * t923;
t772 = -mrSges(5,1) * t840 + mrSges(5,3) * t819 + Ifges(5,4) * t853 + Ifges(5,2) * t852 + Ifges(5,6) * t922 - pkin(4) * t791 - t884 * t861 + t923 * t863 - t970;
t880 = Ifges(4,5) * t901 + Ifges(4,6) * t900 + Ifges(4,3) * qJD(2);
t758 = -mrSges(4,1) * t874 + mrSges(4,3) * t843 + Ifges(4,4) * t890 + Ifges(4,2) * t889 + Ifges(4,6) * qJDD(2) - pkin(3) * t943 + pkin(8) * t948 + qJD(2) * t882 + t929 * t765 + t933 * t772 - t901 * t880;
t759 = mrSges(4,2) * t874 - mrSges(4,3) * t842 + Ifges(4,1) * t890 + Ifges(4,4) * t889 + Ifges(4,5) * qJDD(2) - pkin(8) * t778 - qJD(2) * t881 + t765 * t933 - t772 * t929 + t880 * t900;
t902 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t930 + Ifges(3,6) * t934) * qJD(1);
t754 = -mrSges(3,1) * t905 + mrSges(3,3) * t892 + Ifges(3,4) * t911 + Ifges(3,2) * t912 + Ifges(3,6) * qJDD(2) - pkin(2) * t787 + qJ(3) * t949 + qJD(2) * t904 + t927 * t758 + t926 * t759 - t902 * t956;
t756 = mrSges(3,2) * t905 - mrSges(3,3) * t891 + Ifges(3,1) * t911 + Ifges(3,4) * t912 + Ifges(3,5) * qJDD(2) - qJ(3) * t771 - qJD(2) * t903 - t758 * t926 + t759 * t927 + t902 * t955;
t942 = mrSges(2,1) * t917 - mrSges(2,2) * t918 + Ifges(2,3) * qJDD(1) + pkin(1) * t938 + pkin(7) * t950 + t934 * t754 + t930 * t756;
t757 = mrSges(2,1) * g(3) + mrSges(2,3) * t918 + t936 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t764 - t971;
t752 = -mrSges(2,2) * g(3) - mrSges(2,3) * t917 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t936 - pkin(7) * t764 - t754 * t930 + t756 * t934;
t1 = [-m(1) * g(1) + t951; -m(1) * g(2) + t961; (-m(1) - m(2)) * g(3) + t764; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t961 + t935 * t752 - t931 * t757; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t951 + t931 * t752 + t935 * t757; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t942; t942; t971; t787; -t940; t970; t803;];
tauJB  = t1;
