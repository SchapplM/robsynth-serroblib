% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:24:35
% EndTime: 2019-05-06 01:24:51
% DurationCPUTime: 15.16s
% Computational Cost: add. (223130->344), mult. (535249->414), div. (0->0), fcn. (417579->10), ass. (0->147)
t982 = Ifges(6,1) + Ifges(7,1);
t976 = Ifges(6,4) + Ifges(7,4);
t975 = Ifges(6,5) + Ifges(7,5);
t981 = Ifges(6,2) + Ifges(7,2);
t974 = Ifges(6,6) + Ifges(7,6);
t980 = Ifges(6,3) + Ifges(7,3);
t937 = qJD(1) ^ 2;
t927 = sin(pkin(10));
t928 = cos(pkin(10));
t931 = sin(qJ(3));
t935 = cos(qJ(3));
t947 = -t927 * t931 + t928 * t935;
t906 = t947 * qJD(1);
t948 = t927 * t935 + t928 * t931;
t907 = t948 * qJD(1);
t930 = sin(qJ(4));
t934 = cos(qJ(4));
t889 = t934 * t906 - t930 * t907;
t896 = -t907 * qJD(3) + t947 * qJDD(1);
t964 = t906 * qJD(3);
t897 = t948 * qJDD(1) + t964;
t858 = t889 * qJD(4) + t930 * t896 + t934 * t897;
t890 = t930 * t906 + t934 * t907;
t925 = qJD(3) + qJD(4);
t929 = sin(qJ(5));
t933 = cos(qJ(5));
t876 = -t929 * t890 + t933 * t925;
t922 = qJDD(3) + qJDD(4);
t833 = t876 * qJD(5) + t933 * t858 + t929 * t922;
t877 = t933 * t890 + t929 * t925;
t859 = -t876 * mrSges(7,1) + t877 * mrSges(7,2);
t932 = sin(qJ(1));
t936 = cos(qJ(1));
t913 = -t936 * g(1) - t932 * g(2);
t908 = -t937 * pkin(1) + qJDD(1) * qJ(2) + t913;
t963 = qJD(1) * qJD(2);
t959 = -t928 * g(3) - 0.2e1 * t927 * t963;
t978 = pkin(2) * t928;
t883 = (-pkin(7) * qJDD(1) + t937 * t978 - t908) * t927 + t959;
t899 = -t927 * g(3) + (t908 + 0.2e1 * t963) * t928;
t962 = qJDD(1) * t928;
t924 = t928 ^ 2;
t971 = t924 * t937;
t884 = -pkin(2) * t971 + pkin(7) * t962 + t899;
t861 = t935 * t883 - t931 * t884;
t828 = (-t897 + t964) * pkin(8) + (t906 * t907 + qJDD(3)) * pkin(3) + t861;
t862 = t931 * t883 + t935 * t884;
t902 = qJD(3) * pkin(3) - t907 * pkin(8);
t905 = t906 ^ 2;
t836 = -t905 * pkin(3) + t896 * pkin(8) - qJD(3) * t902 + t862;
t826 = t930 * t828 + t934 * t836;
t874 = -t889 * pkin(4) - t890 * pkin(9);
t921 = t925 ^ 2;
t820 = -t921 * pkin(4) + t922 * pkin(9) + t889 * t874 + t826;
t923 = t927 ^ 2;
t912 = t932 * g(1) - t936 * g(2);
t952 = qJDD(2) - t912;
t895 = (-pkin(1) - t978) * qJDD(1) + (-qJ(2) + (-t923 - t924) * pkin(7)) * t937 + t952;
t849 = -t896 * pkin(3) - t905 * pkin(8) + t907 * t902 + t895;
t857 = -t890 * qJD(4) + t934 * t896 - t930 * t897;
t823 = (-t889 * t925 - t858) * pkin(9) + (t890 * t925 - t857) * pkin(4) + t849;
t815 = -t929 * t820 + t933 * t823;
t856 = qJDD(5) - t857;
t885 = qJD(5) - t889;
t812 = -0.2e1 * qJD(6) * t877 + (t876 * t885 - t833) * qJ(6) + (t876 * t877 + t856) * pkin(5) + t815;
t863 = -t885 * mrSges(7,2) + t876 * mrSges(7,3);
t961 = m(7) * t812 + t856 * mrSges(7,1) + t885 * t863;
t809 = -t833 * mrSges(7,3) - t877 * t859 + t961;
t816 = t933 * t820 + t929 * t823;
t832 = -t877 * qJD(5) - t929 * t858 + t933 * t922;
t865 = t885 * pkin(5) - t877 * qJ(6);
t875 = t876 ^ 2;
t814 = -t875 * pkin(5) + t832 * qJ(6) + 0.2e1 * qJD(6) * t876 - t885 * t865 + t816;
t967 = t976 * t876 + t982 * t877 + t975 * t885;
t968 = -t981 * t876 - t976 * t877 - t974 * t885;
t979 = mrSges(6,1) * t815 + mrSges(7,1) * t812 - mrSges(6,2) * t816 - mrSges(7,2) * t814 + pkin(5) * t809 + t974 * t832 + t975 * t833 + t980 * t856 - t967 * t876 - t968 * t877;
t977 = -mrSges(6,2) - mrSges(7,2);
t972 = mrSges(3,2) * t927;
t873 = -t889 * mrSges(5,1) + t890 * mrSges(5,2);
t881 = t925 * mrSges(5,1) - t890 * mrSges(5,3);
t860 = -t876 * mrSges(6,1) + t877 * mrSges(6,2);
t864 = -t885 * mrSges(6,2) + t876 * mrSges(6,3);
t801 = m(6) * t815 + t856 * mrSges(6,1) + t885 * t864 + (-t859 - t860) * t877 + (-mrSges(6,3) - mrSges(7,3)) * t833 + t961;
t960 = m(7) * t814 + t832 * mrSges(7,3) + t876 * t859;
t866 = t885 * mrSges(7,1) - t877 * mrSges(7,3);
t966 = -t885 * mrSges(6,1) + t877 * mrSges(6,3) - t866;
t804 = m(6) * t816 + t832 * mrSges(6,3) + t977 * t856 + t876 * t860 + t966 * t885 + t960;
t954 = -t929 * t801 + t933 * t804;
t794 = m(5) * t826 - t922 * mrSges(5,2) + t857 * mrSges(5,3) + t889 * t873 - t925 * t881 + t954;
t825 = t934 * t828 - t930 * t836;
t880 = -t925 * mrSges(5,2) + t889 * mrSges(5,3);
t819 = -t922 * pkin(4) - t921 * pkin(9) + t890 * t874 - t825;
t817 = -t832 * pkin(5) - t875 * qJ(6) + t877 * t865 + qJDD(6) + t819;
t953 = -m(7) * t817 + t832 * mrSges(7,1) + t876 * t863;
t940 = -m(6) * t819 + t832 * mrSges(6,1) + t977 * t833 + t876 * t864 + t966 * t877 + t953;
t806 = m(5) * t825 + t922 * mrSges(5,1) - t858 * mrSges(5,3) - t890 * t873 + t925 * t880 + t940;
t785 = t930 * t794 + t934 * t806;
t893 = -t906 * mrSges(4,1) + t907 * mrSges(4,2);
t900 = -qJD(3) * mrSges(4,2) + t906 * mrSges(4,3);
t783 = m(4) * t861 + qJDD(3) * mrSges(4,1) - t897 * mrSges(4,3) + qJD(3) * t900 - t907 * t893 + t785;
t901 = qJD(3) * mrSges(4,1) - t907 * mrSges(4,3);
t955 = t934 * t794 - t930 * t806;
t784 = m(4) * t862 - qJDD(3) * mrSges(4,2) + t896 * mrSges(4,3) - qJD(3) * t901 + t906 * t893 + t955;
t778 = t935 * t783 + t931 * t784;
t898 = -t927 * t908 + t959;
t946 = mrSges(3,3) * qJDD(1) + t937 * (-mrSges(3,1) * t928 + t972);
t776 = m(3) * t898 - t946 * t927 + t778;
t956 = -t931 * t783 + t935 * t784;
t777 = m(3) * t899 + t946 * t928 + t956;
t957 = -t927 * t776 + t928 * t777;
t769 = m(2) * t913 - t937 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t957;
t904 = -qJDD(1) * pkin(1) - t937 * qJ(2) + t952;
t798 = t933 * t801 + t929 * t804;
t945 = m(5) * t849 - t857 * mrSges(5,1) + t858 * mrSges(5,2) - t889 * t880 + t890 * t881 + t798;
t941 = m(4) * t895 - t896 * mrSges(4,1) + t897 * mrSges(4,2) - t906 * t900 + t907 * t901 + t945;
t939 = -m(3) * t904 + mrSges(3,1) * t962 - t941 + (t923 * t937 + t971) * mrSges(3,3);
t789 = -t937 * mrSges(2,2) + m(2) * t912 + t939 + (mrSges(2,1) - t972) * qJDD(1);
t970 = t932 * t769 + t936 * t789;
t771 = t928 * t776 + t927 * t777;
t969 = -t974 * t876 - t975 * t877 - t980 * t885;
t949 = Ifges(3,5) * t927 + Ifges(3,6) * t928;
t965 = t937 * t949;
t958 = t936 * t769 - t932 * t789;
t951 = Ifges(3,1) * t927 + Ifges(3,4) * t928;
t950 = Ifges(3,4) * t927 + Ifges(3,2) * t928;
t810 = t833 * mrSges(7,2) + t877 * t866 - t953;
t787 = -mrSges(6,1) * t819 + mrSges(6,3) * t816 - mrSges(7,1) * t817 + mrSges(7,3) * t814 - pkin(5) * t810 + qJ(6) * t960 + (-qJ(6) * t866 + t967) * t885 + t969 * t877 + (-qJ(6) * mrSges(7,2) + t974) * t856 + t976 * t833 + t981 * t832;
t796 = mrSges(6,2) * t819 + mrSges(7,2) * t817 - mrSges(6,3) * t815 - mrSges(7,3) * t812 - qJ(6) * t809 + t976 * t832 + t982 * t833 + t975 * t856 - t969 * t876 + t968 * t885;
t868 = Ifges(5,5) * t890 + Ifges(5,6) * t889 + Ifges(5,3) * t925;
t869 = Ifges(5,4) * t890 + Ifges(5,2) * t889 + Ifges(5,6) * t925;
t772 = mrSges(5,2) * t849 - mrSges(5,3) * t825 + Ifges(5,1) * t858 + Ifges(5,4) * t857 + Ifges(5,5) * t922 - pkin(9) * t798 - t929 * t787 + t933 * t796 + t889 * t868 - t925 * t869;
t870 = Ifges(5,1) * t890 + Ifges(5,4) * t889 + Ifges(5,5) * t925;
t779 = -mrSges(5,1) * t849 + mrSges(5,3) * t826 + Ifges(5,4) * t858 + Ifges(5,2) * t857 + Ifges(5,6) * t922 - pkin(4) * t798 - t890 * t868 + t925 * t870 - t979;
t886 = Ifges(4,5) * t907 + Ifges(4,6) * t906 + Ifges(4,3) * qJD(3);
t888 = Ifges(4,1) * t907 + Ifges(4,4) * t906 + Ifges(4,5) * qJD(3);
t765 = -mrSges(4,1) * t895 + mrSges(4,3) * t862 + Ifges(4,4) * t897 + Ifges(4,2) * t896 + Ifges(4,6) * qJDD(3) - pkin(3) * t945 + pkin(8) * t955 + qJD(3) * t888 + t930 * t772 + t934 * t779 - t907 * t886;
t887 = Ifges(4,4) * t907 + Ifges(4,2) * t906 + Ifges(4,6) * qJD(3);
t766 = mrSges(4,2) * t895 - mrSges(4,3) * t861 + Ifges(4,1) * t897 + Ifges(4,4) * t896 + Ifges(4,5) * qJDD(3) - pkin(8) * t785 - qJD(3) * t887 + t934 * t772 - t930 * t779 + t906 * t886;
t761 = -mrSges(3,1) * t904 + mrSges(3,3) * t899 - pkin(2) * t941 + pkin(7) * t956 + t950 * qJDD(1) + t935 * t765 + t931 * t766 - t927 * t965;
t763 = mrSges(3,2) * t904 - mrSges(3,3) * t898 - pkin(7) * t778 + t951 * qJDD(1) - t931 * t765 + t935 * t766 + t928 * t965;
t791 = qJDD(1) * t972 - t939;
t944 = mrSges(2,1) * t912 - mrSges(2,2) * t913 + Ifges(2,3) * qJDD(1) - pkin(1) * t791 + qJ(2) * t957 + t928 * t761 + t927 * t763;
t942 = -mrSges(5,1) * t825 + mrSges(5,2) * t826 - Ifges(5,5) * t858 - Ifges(5,6) * t857 - Ifges(5,3) * t922 - pkin(4) * t940 - pkin(9) * t954 - t933 * t787 - t929 * t796 - t890 * t869 + t889 * t870;
t938 = mrSges(4,1) * t861 - mrSges(4,2) * t862 + Ifges(4,5) * t897 + Ifges(4,6) * t896 + Ifges(4,3) * qJDD(3) + pkin(3) * t785 + t907 * t887 - t906 * t888 - t942;
t764 = -t938 + mrSges(2,3) * t913 - mrSges(3,1) * t898 + mrSges(3,2) * t899 + mrSges(2,1) * g(3) - pkin(2) * t778 - pkin(1) * t771 + (Ifges(2,6) - t949) * qJDD(1) + (-t927 * t950 + t928 * t951 + Ifges(2,5)) * t937;
t759 = -mrSges(2,2) * g(3) - mrSges(2,3) * t912 + Ifges(2,5) * qJDD(1) - t937 * Ifges(2,6) - qJ(2) * t771 - t927 * t761 + t928 * t763;
t1 = [-m(1) * g(1) + t958; -m(1) * g(2) + t970; (-m(1) - m(2)) * g(3) + t771; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t970 + t936 * t759 - t932 * t764; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t958 + t932 * t759 + t936 * t764; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t944; t944; t791; t938; -t942; t979; t810;];
tauJB  = t1;
