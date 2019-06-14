% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-05-05 03:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:25:11
% EndTime: 2019-05-05 03:25:19
% DurationCPUTime: 7.97s
% Computational Cost: add. (123276->325), mult. (257235->402), div. (0->0), fcn. (161019->12), ass. (0->144)
t955 = -2 * qJD(4);
t954 = Ifges(4,1) + Ifges(5,2);
t947 = Ifges(4,4) + Ifges(5,6);
t946 = Ifges(4,5) - Ifges(5,4);
t953 = Ifges(4,2) + Ifges(5,3);
t945 = Ifges(4,6) - Ifges(5,5);
t952 = Ifges(4,3) + Ifges(5,1);
t897 = sin(pkin(10));
t900 = cos(pkin(10));
t878 = g(1) * t897 - g(2) * t900;
t895 = -g(3) + qJDD(1);
t898 = sin(pkin(6));
t901 = cos(pkin(6));
t846 = -t878 * t898 + t895 * t901;
t903 = sin(qJ(3));
t906 = cos(qJ(3));
t930 = qJD(2) * qJD(3);
t928 = t906 * t930;
t874 = qJDD(2) * t903 + t928;
t879 = -g(1) * t900 - g(2) * t897;
t907 = cos(qJ(2));
t904 = sin(qJ(2));
t939 = t901 * t904;
t940 = t898 * t904;
t831 = t878 * t939 + t907 * t879 + t895 * t940;
t909 = qJD(2) ^ 2;
t827 = -pkin(2) * t909 + qJDD(2) * pkin(8) + t831;
t824 = t903 * t827;
t871 = (-pkin(3) * t906 - qJ(4) * t903) * qJD(2);
t908 = qJD(3) ^ 2;
t932 = qJD(2) * t903;
t921 = -t908 * qJ(4) + t871 * t932 + qJDD(4) + t824;
t942 = qJ(5) * t909;
t943 = -pkin(3) - qJ(5);
t805 = t874 * pkin(4) + t943 * qJDD(3) + (-pkin(4) * t930 - t903 * t942 - t846) * t906 + t921;
t927 = t903 * t930;
t875 = qJDD(2) * t906 - t927;
t882 = pkin(4) * t932 - qJD(3) * qJ(5);
t894 = t906 ^ 2;
t830 = -t904 * t879 + (t878 * t901 + t895 * t898) * t907;
t914 = -qJDD(2) * pkin(2) - t830;
t913 = pkin(3) * t927 + t932 * t955 + (-t874 - t928) * qJ(4) + t914;
t807 = -t882 * t932 + (-pkin(4) * t894 - pkin(8)) * t909 + t943 * t875 + t913;
t896 = sin(pkin(11));
t899 = cos(pkin(11));
t931 = qJD(2) * t906;
t864 = qJD(3) * t899 - t896 * t931;
t797 = -0.2e1 * qJD(5) * t864 + t899 * t805 - t896 * t807;
t844 = qJDD(3) * t899 - t875 * t896;
t863 = -qJD(3) * t896 - t899 * t931;
t795 = (t863 * t932 - t844) * pkin(9) + (t863 * t864 + t874) * pkin(5) + t797;
t798 = 0.2e1 * qJD(5) * t863 + t896 * t805 + t899 * t807;
t843 = -qJDD(3) * t896 - t875 * t899;
t845 = pkin(5) * t932 - pkin(9) * t864;
t862 = t863 ^ 2;
t796 = -pkin(5) * t862 + pkin(9) * t843 - t845 * t932 + t798;
t902 = sin(qJ(6));
t905 = cos(qJ(6));
t794 = t795 * t902 + t796 * t905;
t822 = t906 * t827 + t903 * t846;
t809 = pkin(3) * t908 - qJDD(3) * qJ(4) + qJD(3) * t955 - t871 * t931 - t822;
t804 = pkin(4) * t875 + qJD(3) * t882 - t894 * t942 + qJDD(5) - t809;
t800 = -pkin(5) * t843 - pkin(9) * t862 + t845 * t864 + t804;
t837 = t863 * t902 + t864 * t905;
t815 = -qJD(6) * t837 + t843 * t905 - t844 * t902;
t836 = t863 * t905 - t864 * t902;
t816 = qJD(6) * t836 + t843 * t902 + t844 * t905;
t887 = qJD(6) + t932;
t817 = Ifges(7,5) * t837 + Ifges(7,6) * t836 + Ifges(7,3) * t887;
t819 = Ifges(7,1) * t837 + Ifges(7,4) * t836 + Ifges(7,5) * t887;
t868 = qJDD(6) + t874;
t780 = -mrSges(7,1) * t800 + mrSges(7,3) * t794 + Ifges(7,4) * t816 + Ifges(7,2) * t815 + Ifges(7,6) * t868 - t817 * t837 + t819 * t887;
t793 = t795 * t905 - t796 * t902;
t818 = Ifges(7,4) * t837 + Ifges(7,2) * t836 + Ifges(7,6) * t887;
t781 = mrSges(7,2) * t800 - mrSges(7,3) * t793 + Ifges(7,1) * t816 + Ifges(7,4) * t815 + Ifges(7,5) * t868 + t817 * t836 - t818 * t887;
t832 = Ifges(6,5) * t864 + Ifges(6,6) * t863 + Ifges(6,3) * t932;
t834 = Ifges(6,1) * t864 + Ifges(6,4) * t863 + Ifges(6,5) * t932;
t828 = -mrSges(7,2) * t887 + mrSges(7,3) * t836;
t829 = mrSges(7,1) * t887 - mrSges(7,3) * t837;
t917 = m(7) * t800 - t815 * mrSges(7,1) + t816 * mrSges(7,2) - t836 * t828 + t837 * t829;
t823 = -mrSges(7,1) * t836 + mrSges(7,2) * t837;
t788 = m(7) * t793 + mrSges(7,1) * t868 - t816 * mrSges(7,3) - t823 * t837 + t828 * t887;
t789 = m(7) * t794 - mrSges(7,2) * t868 + t815 * mrSges(7,3) + t823 * t836 - t829 * t887;
t924 = -t788 * t902 + t905 * t789;
t765 = -mrSges(6,1) * t804 + mrSges(6,3) * t798 + Ifges(6,4) * t844 + Ifges(6,2) * t843 + Ifges(6,6) * t874 - pkin(5) * t917 + pkin(9) * t924 + t905 * t780 + t902 * t781 - t864 * t832 + t834 * t932;
t779 = t905 * t788 + t902 * t789;
t833 = Ifges(6,4) * t864 + Ifges(6,2) * t863 + Ifges(6,6) * t932;
t766 = mrSges(6,2) * t804 - mrSges(6,3) * t797 + Ifges(6,1) * t844 + Ifges(6,4) * t843 + Ifges(6,5) * t874 - pkin(9) * t779 - t780 * t902 + t781 * t905 + t832 * t863 - t833 * t932;
t872 = (mrSges(5,2) * t906 - mrSges(5,3) * t903) * qJD(2);
t883 = -mrSges(5,1) * t931 - qJD(3) * mrSges(5,3);
t838 = -mrSges(6,1) * t863 + mrSges(6,2) * t864;
t841 = -mrSges(6,2) * t932 + mrSges(6,3) * t863;
t777 = m(6) * t797 + mrSges(6,1) * t874 - mrSges(6,3) * t844 - t838 * t864 + t841 * t932 + t779;
t842 = mrSges(6,1) * t932 - mrSges(6,3) * t864;
t778 = m(6) * t798 - mrSges(6,2) * t874 + mrSges(6,3) * t843 + t838 * t863 - t842 * t932 + t924;
t774 = t899 * t777 + t896 * t778;
t938 = t906 * t846;
t810 = -qJDD(3) * pkin(3) + t921 - t938;
t916 = -m(5) * t810 - t874 * mrSges(5,1) - t774;
t773 = qJDD(3) * mrSges(5,2) + qJD(3) * t883 + t872 * t932 - t916;
t821 = -t824 + t938;
t791 = m(6) * t804 - t843 * mrSges(6,1) + t844 * mrSges(6,2) - t863 * t841 + t864 * t842 + t917;
t884 = mrSges(5,1) * t932 + qJD(3) * mrSges(5,2);
t912 = -m(5) * t809 + qJDD(3) * mrSges(5,3) + qJD(3) * t884 + t872 * t931 + t791;
t933 = t946 * qJD(3) + (t954 * t903 + t947 * t906) * qJD(2);
t934 = t945 * qJD(3) + (t947 * t903 + t953 * t906) * qJD(2);
t951 = (t903 * t934 - t906 * t933) * qJD(2) + t952 * qJDD(3) + t946 * t874 + t945 * t875 + mrSges(4,1) * t821 - mrSges(4,2) * t822 + mrSges(5,2) * t810 - mrSges(5,3) * t809 - pkin(3) * t773 + qJ(4) * (mrSges(5,1) * t875 + t912) - qJ(5) * t774 - t896 * t765 + t899 * t766;
t948 = t909 * pkin(8);
t826 = t914 - t948;
t880 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t932;
t881 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t931;
t811 = -t875 * pkin(3) + t913 - t948;
t936 = -t896 * t777 + t899 * t778;
t920 = -m(5) * t811 - t875 * mrSges(5,2) + t884 * t932 - t936;
t911 = -m(4) * t826 + t881 * t931 + t875 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t874 + (-t880 * t903 - t883 * t906) * qJD(2) + t920;
t769 = m(3) * t830 + qJDD(2) * mrSges(3,1) - t909 * mrSges(3,2) + t911;
t941 = t769 * t907;
t873 = (-mrSges(4,1) * t906 + mrSges(4,2) * t903) * qJD(2);
t771 = m(4) * t821 - t874 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t881 - t883) * qJD(3) + (-t872 - t873) * t932 + t916;
t784 = t912 - qJD(3) * t880 + m(4) * t822 - qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t875 + t873 * t931;
t925 = -t771 * t903 + t906 * t784;
t761 = m(3) * t831 - mrSges(3,1) * t909 - qJDD(2) * mrSges(3,2) + t925;
t764 = t906 * t771 + t903 * t784;
t763 = m(3) * t846 + t764;
t752 = t761 * t939 - t763 * t898 + t901 * t941;
t750 = m(2) * t878 + t752;
t757 = t907 * t761 - t769 * t904;
t756 = m(2) * t879 + t757;
t937 = t900 * t750 + t897 * t756;
t935 = t952 * qJD(3) + (t946 * t903 + t945 * t906) * qJD(2);
t751 = t761 * t940 + t901 * t763 + t898 * t941;
t926 = -t750 * t897 + t900 * t756;
t923 = m(2) * t895 + t751;
t772 = -t874 * mrSges(5,3) + t883 * t931 - t920;
t748 = -mrSges(4,1) * t826 - mrSges(5,1) * t809 + mrSges(5,2) * t811 + mrSges(4,3) * t822 - pkin(3) * t772 + pkin(4) * t791 - qJ(5) * t936 + t933 * qJD(3) + t945 * qJDD(3) - t899 * t765 - t896 * t766 + t947 * t874 + t953 * t875 - t935 * t932;
t915 = mrSges(7,1) * t793 - mrSges(7,2) * t794 + Ifges(7,5) * t816 + Ifges(7,6) * t815 + Ifges(7,3) * t868 + t837 * t818 - t836 * t819;
t753 = t915 + t935 * t931 - t863 * t834 + t864 * t833 + Ifges(6,5) * t844 + Ifges(6,6) * t843 + mrSges(4,2) * t826 - mrSges(5,3) * t811 - mrSges(4,3) * t821 + mrSges(5,1) * t810 + mrSges(6,1) * t797 - mrSges(6,2) * t798 + pkin(5) * t779 + pkin(4) * t774 - qJ(4) * t772 + t946 * qJDD(3) + t947 * t875 + (Ifges(6,3) + t954) * t874 - t934 * qJD(3);
t746 = mrSges(3,2) * t846 - mrSges(3,3) * t830 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t909 - pkin(8) * t764 - t748 * t903 + t753 * t906;
t747 = -mrSges(3,1) * t846 + mrSges(3,3) * t831 + t909 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t764 - t951;
t918 = pkin(7) * t757 + t746 * t904 + t747 * t907;
t745 = mrSges(3,1) * t830 - mrSges(3,2) * t831 + Ifges(3,3) * qJDD(2) + pkin(2) * t911 + pkin(8) * t925 + t906 * t748 + t903 * t753;
t744 = mrSges(2,2) * t895 - mrSges(2,3) * t878 + t907 * t746 - t904 * t747 + (-t751 * t898 - t752 * t901) * pkin(7);
t743 = -mrSges(2,1) * t895 + mrSges(2,3) * t879 - pkin(1) * t751 - t898 * t745 + t901 * t918;
t1 = [-m(1) * g(1) + t926; -m(1) * g(2) + t937; -m(1) * g(3) + t923; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t937 - t897 * t743 + t900 * t744; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t926 + t900 * t743 + t897 * t744; -mrSges(1,1) * g(2) + mrSges(2,1) * t878 + mrSges(1,2) * g(1) - mrSges(2,2) * t879 + pkin(1) * t752 + t901 * t745 + t898 * t918; t923; t745; t951; t773; t791; t915;];
tauJB  = t1;
