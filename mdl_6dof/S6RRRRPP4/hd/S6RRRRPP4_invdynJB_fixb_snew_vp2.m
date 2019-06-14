% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:16:18
% EndTime: 2019-05-07 18:16:39
% DurationCPUTime: 19.12s
% Computational Cost: add. (312626->363), mult. (644243->442), div. (0->0), fcn. (464854->10), ass. (0->145)
t988 = Ifges(6,1) + Ifges(7,1);
t978 = Ifges(6,4) - Ifges(7,5);
t986 = Ifges(7,4) + Ifges(6,5);
t987 = Ifges(6,2) + Ifges(7,3);
t984 = Ifges(6,6) - Ifges(7,6);
t985 = -Ifges(7,2) - Ifges(6,3);
t947 = sin(qJ(3));
t951 = cos(qJ(3));
t948 = sin(qJ(2));
t973 = qJD(1) * t948;
t927 = qJD(2) * t951 - t947 * t973;
t928 = qJD(2) * t947 + t951 * t973;
t946 = sin(qJ(4));
t950 = cos(qJ(4));
t901 = t927 * t950 - t928 * t946;
t902 = t927 * t946 + t928 * t950;
t945 = sin(pkin(10));
t977 = cos(pkin(10));
t879 = -t901 * t977 + t945 * t902;
t880 = t945 * t901 + t902 * t977;
t952 = cos(qJ(2));
t972 = qJD(1) * t952;
t940 = qJD(3) - t972;
t939 = qJD(4) + t940;
t983 = t879 * t987 - t880 * t978 - t939 * t984;
t982 = -t978 * t879 + t880 * t988 + t986 * t939;
t949 = sin(qJ(1));
t953 = cos(qJ(1));
t937 = -g(1) * t953 - g(2) * t949;
t955 = qJD(1) ^ 2;
t921 = -pkin(1) * t955 + qJDD(1) * pkin(7) + t937;
t906 = -t952 * g(3) - t948 * t921;
t930 = (-pkin(2) * t952 - pkin(8) * t948) * qJD(1);
t954 = qJD(2) ^ 2;
t886 = -qJDD(2) * pkin(2) - t954 * pkin(8) + t930 * t973 - t906;
t971 = qJD(1) * qJD(2);
t969 = t952 * t971;
t931 = qJDD(1) * t948 + t969;
t898 = -qJD(3) * t928 + qJDD(2) * t951 - t931 * t947;
t908 = pkin(3) * t940 - pkin(9) * t928;
t925 = t927 ^ 2;
t857 = -t898 * pkin(3) - t925 * pkin(9) + t928 * t908 + t886;
t899 = qJD(3) * t927 + qJDD(2) * t947 + t931 * t951;
t863 = -qJD(4) * t902 + t898 * t950 - t899 * t946;
t889 = pkin(4) * t939 - qJ(5) * t902;
t900 = t901 ^ 2;
t826 = -t863 * pkin(4) - t900 * qJ(5) + t902 * t889 + qJDD(5) + t857;
t864 = qJD(4) * t901 + t898 * t946 + t899 * t950;
t837 = -t863 * t977 + t945 * t864;
t838 = t945 * t863 + t864 * t977;
t817 = t826 - 0.2e1 * qJD(6) * t880 + (t879 * t939 - t838) * qJ(6) + (t880 * t939 + t837) * pkin(5);
t869 = -mrSges(7,2) * t879 + mrSges(7,3) * t939;
t872 = -mrSges(7,1) * t939 + mrSges(7,2) * t880;
t808 = m(7) * t817 + t837 * mrSges(7,1) - t838 * mrSges(7,3) + t879 * t869 - t880 * t872;
t936 = t949 * g(1) - t953 * g(2);
t920 = -qJDD(1) * pkin(1) - t955 * pkin(7) - t936;
t941 = t948 * t971;
t932 = qJDD(1) * t952 - t941;
t884 = (-t931 - t969) * pkin(8) + (-t932 + t941) * pkin(2) + t920;
t907 = -g(3) * t948 + t952 * t921;
t887 = -pkin(2) * t954 + qJDD(2) * pkin(8) + t930 * t972 + t907;
t865 = t951 * t884 - t947 * t887;
t926 = qJDD(3) - t932;
t841 = (t927 * t940 - t899) * pkin(9) + (t927 * t928 + t926) * pkin(3) + t865;
t866 = t947 * t884 + t951 * t887;
t843 = -pkin(3) * t925 + pkin(9) * t898 - t908 * t940 + t866;
t823 = t950 * t841 - t946 * t843;
t922 = qJDD(4) + t926;
t819 = (t901 * t939 - t864) * qJ(5) + (t901 * t902 + t922) * pkin(4) + t823;
t824 = t946 * t841 + t950 * t843;
t821 = -pkin(4) * t900 + qJ(5) * t863 - t889 * t939 + t824;
t980 = -2 * qJD(5);
t815 = t945 * t819 + t977 * t821 + t879 * t980;
t854 = pkin(5) * t879 - qJ(6) * t880;
t938 = t939 ^ 2;
t811 = -pkin(5) * t938 + qJ(6) * t922 + 0.2e1 * qJD(6) * t939 - t854 * t879 + t815;
t975 = t984 * t879 - t986 * t880 + t985 * t939;
t793 = -mrSges(6,1) * t826 - mrSges(7,1) * t817 + mrSges(7,2) * t811 + mrSges(6,3) * t815 - pkin(5) * t808 - t837 * t987 + t978 * t838 + t975 * t880 + t984 * t922 + t982 * t939;
t962 = t819 * t977 - t945 * t821;
t812 = -t922 * pkin(5) - t938 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t854) * t880 - t962;
t814 = t880 * t980 + t962;
t794 = mrSges(6,2) * t826 + mrSges(7,2) * t812 - mrSges(6,3) * t814 - mrSges(7,3) * t817 - qJ(6) * t808 - t978 * t837 + t838 * t988 + t975 * t879 + t986 * t922 + t983 * t939;
t870 = -mrSges(6,2) * t939 - mrSges(6,3) * t879;
t871 = mrSges(6,1) * t939 - mrSges(6,3) * t880;
t803 = m(6) * t826 + t837 * mrSges(6,1) + t838 * mrSges(6,2) + t879 * t870 + t880 * t871 + t808;
t873 = Ifges(5,5) * t902 + Ifges(5,6) * t901 + Ifges(5,3) * t939;
t875 = Ifges(5,1) * t902 + Ifges(5,4) * t901 + Ifges(5,5) * t939;
t970 = m(7) * t811 + t922 * mrSges(7,3) + t939 * t872;
t855 = mrSges(7,1) * t879 - mrSges(7,3) * t880;
t974 = -mrSges(6,1) * t879 - mrSges(6,2) * t880 - t855;
t979 = -mrSges(6,3) - mrSges(7,2);
t797 = m(6) * t815 - t922 * mrSges(6,2) + t837 * t979 - t939 * t871 + t879 * t974 + t970;
t964 = -m(7) * t812 + t922 * mrSges(7,1) + t939 * t869;
t799 = m(6) * t814 + t922 * mrSges(6,1) + t838 * t979 + t939 * t870 + t880 * t974 + t964;
t965 = t977 * t797 - t799 * t945;
t780 = -mrSges(5,1) * t857 + mrSges(5,3) * t824 + Ifges(5,4) * t864 + Ifges(5,2) * t863 + Ifges(5,6) * t922 - pkin(4) * t803 + qJ(5) * t965 + t793 * t977 + t945 * t794 - t902 * t873 + t939 * t875;
t792 = t945 * t797 + t977 * t799;
t874 = Ifges(5,4) * t902 + Ifges(5,2) * t901 + Ifges(5,6) * t939;
t781 = mrSges(5,2) * t857 - mrSges(5,3) * t823 + Ifges(5,1) * t864 + Ifges(5,4) * t863 + Ifges(5,5) * t922 - qJ(5) * t792 - t945 * t793 + t794 * t977 + t901 * t873 - t939 * t874;
t891 = Ifges(4,5) * t928 + Ifges(4,6) * t927 + Ifges(4,3) * t940;
t893 = Ifges(4,1) * t928 + Ifges(4,4) * t927 + Ifges(4,5) * t940;
t888 = -mrSges(5,2) * t939 + mrSges(5,3) * t901;
t890 = mrSges(5,1) * t939 - mrSges(5,3) * t902;
t959 = m(5) * t857 - t863 * mrSges(5,1) + t864 * mrSges(5,2) - t901 * t888 + t902 * t890 + t803;
t881 = -mrSges(5,1) * t901 + mrSges(5,2) * t902;
t789 = m(5) * t823 + mrSges(5,1) * t922 - mrSges(5,3) * t864 - t881 * t902 + t888 * t939 + t792;
t790 = m(5) * t824 - mrSges(5,2) * t922 + mrSges(5,3) * t863 + t881 * t901 - t890 * t939 + t965;
t966 = -t789 * t946 + t950 * t790;
t765 = -mrSges(4,1) * t886 + mrSges(4,3) * t866 + Ifges(4,4) * t899 + Ifges(4,2) * t898 + Ifges(4,6) * t926 - pkin(3) * t959 + pkin(9) * t966 + t950 * t780 + t946 * t781 - t928 * t891 + t940 * t893;
t785 = t950 * t789 + t946 * t790;
t892 = Ifges(4,4) * t928 + Ifges(4,2) * t927 + Ifges(4,6) * t940;
t766 = mrSges(4,2) * t886 - mrSges(4,3) * t865 + Ifges(4,1) * t899 + Ifges(4,4) * t898 + Ifges(4,5) * t926 - pkin(9) * t785 - t780 * t946 + t781 * t950 + t891 * t927 - t892 * t940;
t903 = -mrSges(4,1) * t927 + mrSges(4,2) * t928;
t904 = -mrSges(4,2) * t940 + mrSges(4,3) * t927;
t783 = m(4) * t865 + mrSges(4,1) * t926 - mrSges(4,3) * t899 - t903 * t928 + t904 * t940 + t785;
t905 = mrSges(4,1) * t940 - mrSges(4,3) * t928;
t784 = m(4) * t866 - mrSges(4,2) * t926 + mrSges(4,3) * t898 + t903 * t927 - t905 * t940 + t966;
t779 = -t783 * t947 + t951 * t784;
t802 = -m(4) * t886 + t898 * mrSges(4,1) - t899 * mrSges(4,2) + t927 * t904 - t928 * t905 - t959;
t918 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t948 + Ifges(3,2) * t952) * qJD(1);
t919 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t948 + Ifges(3,4) * t952) * qJD(1);
t981 = mrSges(3,1) * t906 - mrSges(3,2) * t907 + Ifges(3,5) * t931 + Ifges(3,6) * t932 + Ifges(3,3) * qJDD(2) + pkin(2) * t802 + pkin(8) * t779 + t951 * t765 + t947 * t766 + (t918 * t948 - t919 * t952) * qJD(1);
t929 = (-mrSges(3,1) * t952 + mrSges(3,2) * t948) * qJD(1);
t934 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t973;
t777 = m(3) * t907 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t932 - qJD(2) * t934 + t929 * t972 + t779;
t935 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t972;
t801 = m(3) * t906 + qJDD(2) * mrSges(3,1) - t931 * mrSges(3,3) + qJD(2) * t935 - t929 * t973 + t802;
t967 = t952 * t777 - t801 * t948;
t769 = m(2) * t937 - mrSges(2,1) * t955 - qJDD(1) * mrSges(2,2) + t967;
t778 = t783 * t951 + t784 * t947;
t960 = -m(3) * t920 + t932 * mrSges(3,1) - mrSges(3,2) * t931 - t934 * t973 + t935 * t972 - t778;
t773 = m(2) * t936 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t955 + t960;
t976 = t949 * t769 + t953 * t773;
t771 = t948 * t777 + t952 * t801;
t968 = t953 * t769 - t773 * t949;
t917 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t948 + Ifges(3,6) * t952) * qJD(1);
t762 = mrSges(3,2) * t920 - mrSges(3,3) * t906 + Ifges(3,1) * t931 + Ifges(3,4) * t932 + Ifges(3,5) * qJDD(2) - pkin(8) * t778 - qJD(2) * t918 - t765 * t947 + t766 * t951 + t917 * t972;
t807 = t838 * mrSges(7,2) + t880 * t855 - t964;
t957 = -mrSges(5,1) * t823 - mrSges(6,1) * t814 + mrSges(7,1) * t812 + mrSges(5,2) * t824 + mrSges(6,2) * t815 - mrSges(7,3) * t811 - Ifges(5,5) * t864 - Ifges(5,6) * t863 - pkin(4) * t792 + pkin(5) * t807 - qJ(6) * t970 - t902 * t874 + t901 * t875 + t983 * t880 + (qJ(6) * t855 - t982) * t879 - t986 * t838 + (qJ(6) * mrSges(7,2) + t984) * t837 + (-Ifges(5,3) + t985) * t922;
t956 = mrSges(4,1) * t865 - mrSges(4,2) * t866 + Ifges(4,5) * t899 + Ifges(4,6) * t898 + Ifges(4,3) * t926 + pkin(3) * t785 + t928 * t892 - t927 * t893 - t957;
t764 = -mrSges(3,1) * t920 + mrSges(3,3) * t907 + Ifges(3,4) * t931 + Ifges(3,2) * t932 + Ifges(3,6) * qJDD(2) - pkin(2) * t778 + qJD(2) * t919 - t917 * t973 - t956;
t961 = mrSges(2,1) * t936 - mrSges(2,2) * t937 + Ifges(2,3) * qJDD(1) + pkin(1) * t960 + pkin(7) * t967 + t948 * t762 + t952 * t764;
t760 = mrSges(2,1) * g(3) + mrSges(2,3) * t937 + t955 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t771 - t981;
t759 = -mrSges(2,2) * g(3) - mrSges(2,3) * t936 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t955 - pkin(7) * t771 + t762 * t952 - t764 * t948;
t1 = [-m(1) * g(1) + t968; -m(1) * g(2) + t976; (-m(1) - m(2)) * g(3) + t771; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t976 + t953 * t759 - t949 * t760; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t968 + t949 * t759 + t953 * t760; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t961; t961; t981; t956; -t957; t803; t807;];
tauJB  = t1;
