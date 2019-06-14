% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:07:23
% EndTime: 2019-05-05 04:07:32
% DurationCPUTime: 5.59s
% Computational Cost: add. (58077->302), mult. (113087->359), div. (0->0), fcn. (66869->10), ass. (0->137)
t966 = Ifges(6,1) + Ifges(7,1);
t947 = Ifges(6,4) - Ifges(7,5);
t962 = Ifges(7,4) + Ifges(6,5);
t965 = Ifges(6,2) + Ifges(7,3);
t960 = Ifges(6,6) - Ifges(7,6);
t964 = -2 * qJD(4);
t963 = Ifges(4,1) + Ifges(5,2);
t948 = Ifges(4,4) + Ifges(5,6);
t946 = Ifges(4,5) - Ifges(5,4);
t961 = Ifges(4,2) + Ifges(5,3);
t945 = Ifges(4,6) - Ifges(5,5);
t959 = Ifges(4,3) + Ifges(5,1);
t958 = Ifges(6,3) + Ifges(7,2);
t901 = sin(qJ(5));
t904 = cos(qJ(5));
t905 = cos(qJ(3));
t931 = qJD(2) * t905;
t869 = qJD(3) * t901 + t904 * t931;
t870 = qJD(3) * t904 - t901 * t931;
t902 = sin(qJ(3));
t932 = qJD(2) * t902;
t888 = qJD(5) + t932;
t957 = t869 * t965 - t870 * t947 - t888 * t960;
t956 = -t947 * t869 + t870 * t966 + t962 * t888;
t872 = (mrSges(5,2) * t905 - mrSges(5,3) * t902) * qJD(2);
t882 = -mrSges(5,1) * t931 - qJD(3) * mrSges(5,3);
t897 = sin(pkin(10));
t899 = cos(pkin(10));
t878 = g(1) * t897 - g(2) * t899;
t896 = -g(3) + qJDD(1);
t898 = sin(pkin(6));
t900 = cos(pkin(6));
t845 = -t878 * t898 + t896 * t900;
t930 = qJD(2) * qJD(3);
t927 = t905 * t930;
t874 = qJDD(2) * t902 + t927;
t879 = -g(1) * t899 - g(2) * t897;
t906 = cos(qJ(2));
t903 = sin(qJ(2));
t941 = t900 * t903;
t942 = t898 * t903;
t818 = t878 * t941 + t906 * t879 + t896 * t942;
t908 = qJD(2) ^ 2;
t816 = -pkin(2) * t908 + qJDD(2) * pkin(8) + t818;
t813 = t902 * t816;
t871 = (-pkin(3) * t905 - qJ(4) * t902) * qJD(2);
t907 = qJD(3) ^ 2;
t920 = -t907 * qJ(4) + t871 * t932 + qJDD(4) + t813;
t951 = pkin(9) * t908;
t952 = -pkin(3) - pkin(9);
t803 = t874 * pkin(4) + t952 * qJDD(3) + (-pkin(4) * t930 - t902 * t951 - t845) * t905 + t920;
t926 = t902 * t930;
t875 = qJDD(2) * t905 - t926;
t885 = pkin(4) * t932 - qJD(3) * pkin(9);
t895 = t905 ^ 2;
t817 = -t903 * t879 + (t878 * t900 + t896 * t898) * t906;
t915 = -qJDD(2) * pkin(2) - t817;
t913 = pkin(3) * t926 + t932 * t964 + (-t874 - t927) * qJ(4) + t915;
t805 = -t885 * t932 + (-pkin(4) * t895 - pkin(8)) * t908 + t952 * t875 + t913;
t798 = t901 * t803 + t904 * t805;
t833 = qJD(5) * t870 + qJDD(3) * t901 + t904 * t875;
t843 = mrSges(6,1) * t888 - mrSges(6,3) * t870;
t866 = qJDD(5) + t874;
t837 = pkin(5) * t869 - qJ(6) * t870;
t886 = t888 ^ 2;
t793 = -pkin(5) * t886 + qJ(6) * t866 + 0.2e1 * qJD(6) * t888 - t837 * t869 + t798;
t844 = -mrSges(7,1) * t888 + mrSges(7,2) * t870;
t928 = m(7) * t793 + t866 * mrSges(7,3) + t888 * t844;
t838 = mrSges(7,1) * t869 - mrSges(7,3) * t870;
t936 = -mrSges(6,1) * t869 - mrSges(6,2) * t870 - t838;
t949 = -mrSges(6,3) - mrSges(7,2);
t783 = m(6) * t798 - t866 * mrSges(6,2) + t949 * t833 - t888 * t843 + t936 * t869 + t928;
t797 = t803 * t904 - t805 * t901;
t834 = -qJD(5) * t869 + qJDD(3) * t904 - t875 * t901;
t841 = -mrSges(6,2) * t888 - mrSges(6,3) * t869;
t794 = -pkin(5) * t866 - qJ(6) * t886 + t837 * t870 + qJDD(6) - t797;
t842 = -mrSges(7,2) * t869 + mrSges(7,3) * t888;
t922 = -m(7) * t794 + t866 * mrSges(7,1) + t888 * t842;
t785 = m(6) * t797 + t866 * mrSges(6,1) + t949 * t834 + t888 * t841 + t936 * t870 + t922;
t778 = t901 * t783 + t904 * t785;
t940 = t905 * t845;
t807 = -qJDD(3) * pkin(3) + t920 - t940;
t916 = -m(5) * t807 - t874 * mrSges(5,1) - t778;
t775 = qJDD(3) * mrSges(5,2) + qJD(3) * t882 + t872 * t932 - t916;
t810 = t905 * t816 + t902 * t845;
t806 = t907 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t964 - t871 * t931 - t810;
t802 = t875 * pkin(4) + qJD(3) * t885 - t895 * t951 - t806;
t799 = -0.2e1 * qJD(6) * t870 + (t869 * t888 - t834) * qJ(6) + (t870 * t888 + t833) * pkin(5) + t802;
t790 = m(7) * t799 + t833 * mrSges(7,1) - t834 * mrSges(7,3) + t869 * t842 - t870 * t844;
t937 = t960 * t869 - t962 * t870 - t958 * t888;
t776 = -mrSges(6,1) * t802 - mrSges(7,1) * t799 + mrSges(7,2) * t793 + mrSges(6,3) * t798 - pkin(5) * t790 - t833 * t965 + t947 * t834 + t960 * t866 + t937 * t870 + t956 * t888;
t777 = mrSges(6,2) * t802 + mrSges(7,2) * t794 - mrSges(6,3) * t797 - mrSges(7,3) * t799 - qJ(6) * t790 - t947 * t833 + t834 * t966 + t962 * t866 + t937 * t869 + t957 * t888;
t809 = -t813 + t940;
t883 = mrSges(5,1) * t932 + qJD(3) * mrSges(5,2);
t914 = m(6) * t802 + t833 * mrSges(6,1) + t834 * mrSges(6,2) + t869 * t841 + t870 * t843 + t790;
t912 = -m(5) * t806 + qJDD(3) * mrSges(5,3) + qJD(3) * t883 + t872 * t931 + t914;
t933 = t946 * qJD(3) + (t963 * t902 + t948 * t905) * qJD(2);
t934 = t945 * qJD(3) + (t948 * t902 + t961 * t905) * qJD(2);
t955 = (t934 * t902 - t933 * t905) * qJD(2) + t959 * qJDD(3) + t946 * t874 + t945 * t875 + mrSges(4,1) * t809 - mrSges(4,2) * t810 + mrSges(5,2) * t807 - mrSges(5,3) * t806 - pkin(3) * t775 - pkin(9) * t778 + qJ(4) * (mrSges(5,1) * t875 + t912) - t901 * t776 + t904 * t777;
t950 = t908 * pkin(8);
t815 = t915 - t950;
t880 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t932;
t881 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t931;
t808 = -t875 * pkin(3) + t913 - t950;
t938 = t904 * t783 - t901 * t785;
t919 = -m(5) * t808 - t875 * mrSges(5,2) + t883 * t932 - t938;
t911 = -m(4) * t815 + t881 * t931 + t875 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t874 + (-t880 * t902 - t882 * t905) * qJD(2) + t919;
t771 = m(3) * t817 + qJDD(2) * mrSges(3,1) - t908 * mrSges(3,2) + t911;
t943 = t771 * t906;
t873 = (-mrSges(4,1) * t905 + mrSges(4,2) * t902) * qJD(2);
t773 = m(4) * t809 - t874 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t881 - t882) * qJD(3) + (-t872 - t873) * t932 + t916;
t781 = t912 + (mrSges(4,3) + mrSges(5,1)) * t875 - qJDD(3) * mrSges(4,2) - qJD(3) * t880 + t873 * t931 + m(4) * t810;
t924 = -t773 * t902 + t905 * t781;
t765 = m(3) * t818 - mrSges(3,1) * t908 - qJDD(2) * mrSges(3,2) + t924;
t768 = t905 * t773 + t902 * t781;
t767 = m(3) * t845 + t768;
t756 = t765 * t941 - t767 * t898 + t900 * t943;
t754 = m(2) * t878 + t756;
t761 = t906 * t765 - t771 * t903;
t760 = m(2) * t879 + t761;
t939 = t899 * t754 + t897 * t760;
t935 = t959 * qJD(3) + (t946 * t902 + t945 * t905) * qJD(2);
t755 = t765 * t942 + t900 * t767 + t898 * t943;
t925 = -t754 * t897 + t899 * t760;
t923 = m(2) * t896 + t755;
t774 = -t874 * mrSges(5,3) + t882 * t931 - t919;
t752 = -mrSges(4,1) * t815 - mrSges(5,1) * t806 + mrSges(5,2) * t808 + mrSges(4,3) * t810 - pkin(3) * t774 + pkin(4) * t914 - pkin(9) * t938 + t933 * qJD(3) + t945 * qJDD(3) - t904 * t776 - t901 * t777 + t948 * t874 + t961 * t875 - t935 * t932;
t789 = t834 * mrSges(7,2) + t870 * t838 - t922;
t909 = mrSges(6,1) * t797 - mrSges(7,1) * t794 - mrSges(6,2) * t798 + mrSges(7,3) * t793 - pkin(5) * t789 + qJ(6) * t928 - t957 * t870 + (-qJ(6) * t838 + t956) * t869 + t958 * t866 + t962 * t834 + (-mrSges(7,2) * qJ(6) - t960) * t833;
t757 = mrSges(5,1) * t807 + mrSges(4,2) * t815 - mrSges(4,3) * t809 - mrSges(5,3) * t808 + pkin(4) * t778 - qJ(4) * t774 - t934 * qJD(3) + t946 * qJDD(3) + t963 * t874 + t948 * t875 + t935 * t931 + t909;
t750 = mrSges(3,2) * t845 - mrSges(3,3) * t817 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t908 - pkin(8) * t768 - t752 * t902 + t757 * t905;
t751 = -mrSges(3,1) * t845 + mrSges(3,3) * t818 + t908 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t768 - t955;
t917 = pkin(7) * t761 + t750 * t903 + t751 * t906;
t749 = mrSges(3,1) * t817 - mrSges(3,2) * t818 + Ifges(3,3) * qJDD(2) + pkin(2) * t911 + pkin(8) * t924 + t905 * t752 + t902 * t757;
t748 = mrSges(2,2) * t896 - mrSges(2,3) * t878 + t906 * t750 - t903 * t751 + (-t755 * t898 - t756 * t900) * pkin(7);
t747 = -mrSges(2,1) * t896 + mrSges(2,3) * t879 - pkin(1) * t755 - t898 * t749 + t900 * t917;
t1 = [-m(1) * g(1) + t925; -m(1) * g(2) + t939; -m(1) * g(3) + t923; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t939 - t897 * t747 + t899 * t748; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t925 + t899 * t747 + t897 * t748; -mrSges(1,1) * g(2) + mrSges(2,1) * t878 + mrSges(1,2) * g(1) - mrSges(2,2) * t879 + pkin(1) * t756 + t900 * t749 + t898 * t917; t923; t749; t955; t775; t909; t789;];
tauJB  = t1;
