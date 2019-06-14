% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:47:50
% EndTime: 2019-05-05 17:48:03
% DurationCPUTime: 12.76s
% Computational Cost: add. (185033->342), mult. (454225->418), div. (0->0), fcn. (339341->10), ass. (0->146)
t972 = Ifges(6,1) + Ifges(7,1);
t964 = Ifges(6,4) - Ifges(7,5);
t963 = Ifges(7,4) + Ifges(6,5);
t971 = Ifges(6,2) + Ifges(7,3);
t970 = Ifges(7,2) + Ifges(6,3);
t961 = Ifges(6,6) - Ifges(7,6);
t926 = qJD(1) ^ 2;
t918 = sin(pkin(9));
t920 = cos(pkin(9));
t922 = sin(qJ(3));
t968 = cos(qJ(3));
t933 = t918 * t968 + t920 * t922;
t899 = t933 * qJD(1);
t917 = sin(pkin(10));
t919 = cos(pkin(10));
t890 = qJD(3) * t919 - t899 * t917;
t891 = qJD(3) * t917 + t899 * t919;
t921 = sin(qJ(5));
t967 = cos(qJ(5));
t856 = -t890 * t967 + t921 * t891;
t946 = t920 * t968;
t953 = qJD(1) * t918;
t898 = -qJD(1) * t946 + t922 * t953;
t951 = t898 * qJD(3);
t883 = qJDD(1) * t933 - t951;
t868 = qJDD(3) * t919 - t883 * t917;
t869 = qJDD(3) * t917 + t883 * t919;
t827 = -t856 * qJD(5) + t921 * t868 + t869 * t967;
t857 = t921 * t890 + t891 * t967;
t841 = mrSges(7,1) * t856 - mrSges(7,3) * t857;
t923 = sin(qJ(1));
t924 = cos(qJ(1));
t905 = -g(1) * t924 - g(2) * t923;
t900 = -pkin(1) * t926 + qJDD(1) * qJ(2) + t905;
t950 = qJD(1) * qJD(2);
t945 = -t920 * g(3) - 0.2e1 * t918 * t950;
t966 = pkin(2) * t920;
t864 = (-pkin(7) * qJDD(1) + t926 * t966 - t900) * t918 + t945;
t885 = -g(3) * t918 + (t900 + 0.2e1 * t950) * t920;
t948 = qJDD(1) * t920;
t915 = t920 ^ 2;
t959 = t915 * t926;
t870 = -pkin(2) * t959 + pkin(7) * t948 + t885;
t846 = t922 * t864 + t968 * t870;
t877 = pkin(3) * t898 - qJ(4) * t899;
t925 = qJD(3) ^ 2;
t828 = -pkin(3) * t925 + qJDD(3) * qJ(4) - t877 * t898 + t846;
t914 = t918 ^ 2;
t904 = t923 * g(1) - t924 * g(2);
t939 = qJDD(2) - t904;
t881 = (-pkin(1) - t966) * qJDD(1) + (-qJ(2) + (-t914 - t915) * pkin(7)) * t926 + t939;
t949 = qJDD(1) * t918;
t952 = qJD(3) * t899;
t882 = -qJDD(1) * t946 + t922 * t949 + t952;
t831 = (-t883 + t951) * qJ(4) + (t882 + t952) * pkin(3) + t881;
t819 = -0.2e1 * qJD(4) * t891 - t917 * t828 + t919 * t831;
t816 = (t890 * t898 - t869) * pkin(8) + (t890 * t891 + t882) * pkin(4) + t819;
t820 = 0.2e1 * qJD(4) * t890 + t919 * t828 + t917 * t831;
t867 = pkin(4) * t898 - pkin(8) * t891;
t889 = t890 ^ 2;
t818 = -pkin(4) * t889 + pkin(8) * t868 - t867 * t898 + t820;
t812 = t816 * t967 - t921 * t818;
t840 = pkin(5) * t856 - qJ(6) * t857;
t880 = qJDD(5) + t882;
t896 = qJD(5) + t898;
t895 = t896 ^ 2;
t811 = -t880 * pkin(5) - t895 * qJ(6) + t857 * t840 + qJDD(6) - t812;
t847 = -mrSges(7,2) * t856 + mrSges(7,3) * t896;
t940 = -m(7) * t811 + t880 * mrSges(7,1) + t896 * t847;
t807 = t827 * mrSges(7,2) + t857 * t841 - t940;
t813 = t921 * t816 + t967 * t818;
t810 = -pkin(5) * t895 + qJ(6) * t880 + 0.2e1 * qJD(6) * t896 - t840 * t856 + t813;
t826 = t857 * qJD(5) - t868 * t967 + t921 * t869;
t850 = -mrSges(7,1) * t896 + mrSges(7,2) * t857;
t947 = m(7) * t810 + t880 * mrSges(7,3) + t896 * t850;
t955 = -t964 * t856 + t972 * t857 + t963 * t896;
t957 = t971 * t856 - t964 * t857 - t961 * t896;
t969 = t961 * t826 - t963 * t827 - t955 * t856 + t957 * t857 - t970 * t880 - mrSges(6,1) * t812 + mrSges(7,1) * t811 + mrSges(6,2) * t813 - mrSges(7,3) * t810 + pkin(5) * t807 - qJ(6) * (-t826 * mrSges(7,2) - t856 * t841 + t947);
t965 = -mrSges(6,3) - mrSges(7,2);
t960 = mrSges(3,2) * t918;
t849 = mrSges(6,1) * t896 - mrSges(6,3) * t857;
t954 = -mrSges(6,1) * t856 - mrSges(6,2) * t857 - t841;
t800 = m(6) * t813 - t880 * mrSges(6,2) + t826 * t965 - t896 * t849 + t856 * t954 + t947;
t848 = -mrSges(6,2) * t896 - mrSges(6,3) * t856;
t802 = m(6) * t812 + t880 * mrSges(6,1) + t827 * t965 + t896 * t848 + t857 * t954 + t940;
t795 = t921 * t800 + t967 * t802;
t858 = -mrSges(5,1) * t890 + mrSges(5,2) * t891;
t938 = -mrSges(5,2) * t898 + mrSges(5,3) * t890;
t793 = m(5) * t819 + t882 * mrSges(5,1) - t869 * mrSges(5,3) - t891 * t858 + t898 * t938 + t795;
t866 = mrSges(5,1) * t898 - mrSges(5,3) * t891;
t941 = t967 * t800 - t802 * t921;
t794 = m(5) * t820 - mrSges(5,2) * t882 + mrSges(5,3) * t868 + t858 * t890 - t866 * t898 + t941;
t789 = -t793 * t917 + t919 * t794;
t878 = mrSges(4,1) * t898 + mrSges(4,2) * t899;
t893 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t899;
t787 = m(4) * t846 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t882 - qJD(3) * t893 - t878 * t898 + t789;
t845 = t864 * t968 - t922 * t870;
t825 = -qJDD(3) * pkin(3) - t925 * qJ(4) + t899 * t877 + qJDD(4) - t845;
t821 = -t868 * pkin(4) - t889 * pkin(8) + t891 * t867 + t825;
t814 = -0.2e1 * qJD(6) * t857 + (t856 * t896 - t827) * qJ(6) + (t857 * t896 + t826) * pkin(5) + t821;
t808 = m(7) * t814 + t826 * mrSges(7,1) - t827 * mrSges(7,3) + t856 * t847 - t857 * t850;
t928 = m(6) * t821 + t826 * mrSges(6,1) + t827 * mrSges(6,2) + t856 * t848 + t857 * t849 + t808;
t805 = m(5) * t825 - t868 * mrSges(5,1) + t869 * mrSges(5,2) + t891 * t866 - t890 * t938 + t928;
t892 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t898;
t804 = m(4) * t845 + qJDD(3) * mrSges(4,1) - t883 * mrSges(4,3) + qJD(3) * t892 - t899 * t878 - t805;
t780 = t922 * t787 + t968 * t804;
t884 = -t918 * t900 + t945;
t934 = mrSges(3,3) * qJDD(1) + t926 * (-mrSges(3,1) * t920 + t960);
t778 = m(3) * t884 - t918 * t934 + t780;
t942 = t968 * t787 - t922 * t804;
t779 = m(3) * t885 + t920 * t934 + t942;
t943 = -t778 * t918 + t920 * t779;
t769 = m(2) * t905 - mrSges(2,1) * t926 - qJDD(1) * mrSges(2,2) + t943;
t897 = -qJDD(1) * pkin(1) - t926 * qJ(2) + t939;
t788 = t919 * t793 + t917 * t794;
t931 = m(4) * t881 + t882 * mrSges(4,1) + t883 * mrSges(4,2) + t898 * t892 + t899 * t893 + t788;
t929 = -m(3) * t897 + mrSges(3,1) * t948 - t931 + (t914 * t926 + t959) * mrSges(3,3);
t782 = t929 + (mrSges(2,1) - t960) * qJDD(1) - t926 * mrSges(2,2) + m(2) * t904;
t958 = t923 * t769 + t924 * t782;
t771 = t920 * t778 + t918 * t779;
t956 = t961 * t856 - t963 * t857 - t970 * t896;
t944 = t924 * t769 - t782 * t923;
t937 = Ifges(3,1) * t918 + Ifges(3,4) * t920;
t936 = Ifges(3,4) * t918 + Ifges(3,2) * t920;
t935 = Ifges(3,5) * t918 + Ifges(3,6) * t920;
t796 = -mrSges(6,1) * t821 - mrSges(7,1) * t814 + mrSges(7,2) * t810 + mrSges(6,3) * t813 - pkin(5) * t808 - t971 * t826 + t964 * t827 + t956 * t857 + t961 * t880 + t955 * t896;
t797 = mrSges(6,2) * t821 + mrSges(7,2) * t811 - mrSges(6,3) * t812 - mrSges(7,3) * t814 - qJ(6) * t808 - t964 * t826 + t972 * t827 + t956 * t856 + t963 * t880 + t957 * t896;
t851 = Ifges(5,5) * t891 + Ifges(5,6) * t890 + Ifges(5,3) * t898;
t853 = Ifges(5,1) * t891 + Ifges(5,4) * t890 + Ifges(5,5) * t898;
t773 = -mrSges(5,1) * t825 + mrSges(5,3) * t820 + Ifges(5,4) * t869 + Ifges(5,2) * t868 + Ifges(5,6) * t882 - pkin(4) * t928 + pkin(8) * t941 + t796 * t967 + t921 * t797 - t891 * t851 + t898 * t853;
t852 = Ifges(5,4) * t891 + Ifges(5,2) * t890 + Ifges(5,6) * t898;
t774 = mrSges(5,2) * t825 - mrSges(5,3) * t819 + Ifges(5,1) * t869 + Ifges(5,4) * t868 + Ifges(5,5) * t882 - pkin(8) * t795 - t921 * t796 + t797 * t967 + t890 * t851 - t898 * t852;
t871 = Ifges(4,5) * t899 - Ifges(4,6) * t898 + Ifges(4,3) * qJD(3);
t872 = Ifges(4,4) * t899 - Ifges(4,2) * t898 + Ifges(4,6) * qJD(3);
t766 = mrSges(4,2) * t881 - mrSges(4,3) * t845 + Ifges(4,1) * t883 - Ifges(4,4) * t882 + Ifges(4,5) * qJDD(3) - qJ(4) * t788 - qJD(3) * t872 - t773 * t917 + t774 * t919 - t871 * t898;
t873 = Ifges(4,1) * t899 - Ifges(4,4) * t898 + Ifges(4,5) * qJD(3);
t772 = -pkin(3) * t788 + Ifges(4,6) * qJDD(3) + t969 + (-Ifges(5,3) - Ifges(4,2)) * t882 - t899 * t871 + t890 * t853 - t891 * t852 + Ifges(4,4) * t883 - mrSges(4,1) * t881 + qJD(3) * t873 - Ifges(5,6) * t868 - Ifges(5,5) * t869 + mrSges(4,3) * t846 - mrSges(5,1) * t819 + mrSges(5,2) * t820 - pkin(4) * t795;
t902 = t935 * qJD(1);
t762 = -mrSges(3,1) * t897 + mrSges(3,3) * t885 - pkin(2) * t931 + pkin(7) * t942 + qJDD(1) * t936 + t922 * t766 + t772 * t968 - t902 * t953;
t765 = t920 * qJD(1) * t902 + mrSges(3,2) * t897 - mrSges(3,3) * t884 - pkin(7) * t780 + qJDD(1) * t937 + t766 * t968 - t922 * t772;
t784 = mrSges(3,2) * t949 - t929;
t932 = mrSges(2,1) * t904 - mrSges(2,2) * t905 + Ifges(2,3) * qJDD(1) - pkin(1) * t784 + qJ(2) * t943 + t920 * t762 + t918 * t765;
t927 = mrSges(4,1) * t845 - mrSges(4,2) * t846 + Ifges(4,5) * t883 - Ifges(4,6) * t882 + Ifges(4,3) * qJDD(3) - pkin(3) * t805 + qJ(4) * t789 + t919 * t773 + t917 * t774 + t899 * t872 + t898 * t873;
t763 = -pkin(1) * t771 + mrSges(2,1) * g(3) + (Ifges(2,6) - t935) * qJDD(1) - pkin(2) * t780 + mrSges(2,3) * t905 - mrSges(3,1) * t884 + mrSges(3,2) * t885 - t927 + (-t918 * t936 + t920 * t937 + Ifges(2,5)) * t926;
t760 = -mrSges(2,2) * g(3) - mrSges(2,3) * t904 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t926 - qJ(2) * t771 - t762 * t918 + t765 * t920;
t1 = [-m(1) * g(1) + t944; -m(1) * g(2) + t958; (-m(1) - m(2)) * g(3) + t771; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t958 + t924 * t760 - t923 * t763; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t944 + t923 * t760 + t924 * t763; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t932; t932; t784; t927; t805; -t969; t807;];
tauJB  = t1;
