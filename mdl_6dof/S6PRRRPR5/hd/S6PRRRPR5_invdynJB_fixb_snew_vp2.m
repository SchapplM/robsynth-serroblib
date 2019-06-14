% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR5
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
% Datum: 2019-05-05 08:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:02:20
% EndTime: 2019-05-05 08:03:09
% DurationCPUTime: 47.65s
% Computational Cost: add. (837242->353), mult. (1772097->468), div. (0->0), fcn. (1433736->16), ass. (0->159)
t903 = sin(pkin(7));
t911 = sin(qJ(3));
t915 = cos(qJ(3));
t932 = qJD(2) * qJD(3);
t886 = (-qJDD(2) * t915 + t911 * t932) * t903;
t902 = sin(pkin(12));
t906 = cos(pkin(12));
t892 = g(1) * t902 - g(2) * t906;
t893 = -g(1) * t906 - g(2) * t902;
t900 = -g(3) + qJDD(1);
t912 = sin(qJ(2));
t908 = cos(pkin(6));
t916 = cos(qJ(2));
t935 = t908 * t916;
t904 = sin(pkin(6));
t938 = t904 * t916;
t865 = t892 * t935 - t893 * t912 + t900 * t938;
t917 = qJD(2) ^ 2;
t943 = pkin(9) * t903;
t862 = qJDD(2) * pkin(2) + t917 * t943 + t865;
t936 = t908 * t912;
t939 = t904 * t912;
t866 = t892 * t936 + t916 * t893 + t900 * t939;
t863 = -pkin(2) * t917 + qJDD(2) * t943 + t866;
t879 = -t892 * t904 + t900 * t908;
t907 = cos(pkin(7));
t828 = -t911 * t863 + (t862 * t907 + t879 * t903) * t915;
t937 = t907 * t911;
t940 = t903 * t911;
t829 = t862 * t937 + t915 * t863 + t879 * t940;
t933 = qJD(2) * t903;
t884 = (-pkin(3) * t915 - pkin(10) * t911) * t933;
t898 = qJD(2) * t907 + qJD(3);
t896 = t898 ^ 2;
t897 = qJDD(2) * t907 + qJDD(3);
t930 = t915 * t933;
t823 = -pkin(3) * t896 + pkin(10) * t897 + t884 * t930 + t829;
t875 = t907 * t879;
t885 = (qJDD(2) * t911 + t915 * t932) * t903;
t826 = pkin(3) * t886 - pkin(10) * t885 + t875 + (-t862 + (pkin(3) * t911 - pkin(10) * t915) * t898 * qJD(2)) * t903;
t910 = sin(qJ(4));
t914 = cos(qJ(4));
t811 = -t823 * t910 + t914 * t826;
t931 = t911 * t933;
t877 = t898 * t914 - t910 * t931;
t854 = qJD(4) * t877 + t885 * t914 + t897 * t910;
t878 = t898 * t910 + t914 * t931;
t880 = qJDD(4) + t886;
t891 = qJD(4) - t930;
t808 = (t877 * t891 - t854) * qJ(5) + (t877 * t878 + t880) * pkin(4) + t811;
t812 = t914 * t823 + t910 * t826;
t853 = -qJD(4) * t878 - t885 * t910 + t897 * t914;
t868 = pkin(4) * t891 - qJ(5) * t878;
t876 = t877 ^ 2;
t810 = -pkin(4) * t876 + qJ(5) * t853 - t868 * t891 + t812;
t901 = sin(pkin(13));
t905 = cos(pkin(13));
t860 = t877 * t905 - t878 * t901;
t944 = 2 * qJD(5);
t805 = t901 * t808 + t905 * t810 + t860 * t944;
t861 = t877 * t901 + t878 * t905;
t841 = -pkin(5) * t860 - pkin(11) * t861;
t890 = t891 ^ 2;
t803 = -pkin(5) * t890 + pkin(11) * t880 + t841 * t860 + t805;
t822 = -pkin(3) * t897 - pkin(10) * t896 + t884 * t931 - t828;
t813 = -pkin(4) * t853 - qJ(5) * t876 + t878 * t868 + qJDD(5) + t822;
t834 = t853 * t905 - t854 * t901;
t835 = t853 * t901 + t854 * t905;
t806 = (-t860 * t891 - t835) * pkin(11) + (t861 * t891 - t834) * pkin(5) + t813;
t909 = sin(qJ(6));
t913 = cos(qJ(6));
t800 = -t803 * t909 + t806 * t913;
t843 = -t861 * t909 + t891 * t913;
t816 = qJD(6) * t843 + t835 * t913 + t880 * t909;
t844 = t861 * t913 + t891 * t909;
t827 = -mrSges(7,1) * t843 + mrSges(7,2) * t844;
t857 = qJD(6) - t860;
t830 = -mrSges(7,2) * t857 + mrSges(7,3) * t843;
t833 = qJDD(6) - t834;
t797 = m(7) * t800 + mrSges(7,1) * t833 - mrSges(7,3) * t816 - t827 * t844 + t830 * t857;
t801 = t803 * t913 + t806 * t909;
t815 = -qJD(6) * t844 - t835 * t909 + t880 * t913;
t831 = mrSges(7,1) * t857 - mrSges(7,3) * t844;
t798 = m(7) * t801 - mrSges(7,2) * t833 + mrSges(7,3) * t815 + t827 * t843 - t831 * t857;
t789 = -t797 * t909 + t913 * t798;
t840 = -mrSges(6,1) * t860 + mrSges(6,2) * t861;
t846 = mrSges(6,1) * t891 - mrSges(6,3) * t861;
t786 = m(6) * t805 - mrSges(6,2) * t880 + mrSges(6,3) * t834 + t840 * t860 - t846 * t891 + t789;
t925 = -t808 * t905 + t810 * t901;
t802 = -pkin(5) * t880 - pkin(11) * t890 + (t944 + t841) * t861 + t925;
t799 = -m(7) * t802 + t815 * mrSges(7,1) - mrSges(7,2) * t816 + t843 * t830 - t831 * t844;
t804 = -0.2e1 * qJD(5) * t861 - t925;
t845 = -mrSges(6,2) * t891 + mrSges(6,3) * t860;
t793 = m(6) * t804 + mrSges(6,1) * t880 - mrSges(6,3) * t835 - t840 * t861 + t845 * t891 + t799;
t780 = t901 * t786 + t905 * t793;
t817 = Ifges(7,5) * t844 + Ifges(7,6) * t843 + Ifges(7,3) * t857;
t819 = Ifges(7,1) * t844 + Ifges(7,4) * t843 + Ifges(7,5) * t857;
t790 = -mrSges(7,1) * t802 + mrSges(7,3) * t801 + Ifges(7,4) * t816 + Ifges(7,2) * t815 + Ifges(7,6) * t833 - t817 * t844 + t819 * t857;
t818 = Ifges(7,4) * t844 + Ifges(7,2) * t843 + Ifges(7,6) * t857;
t791 = mrSges(7,2) * t802 - mrSges(7,3) * t800 + Ifges(7,1) * t816 + Ifges(7,4) * t815 + Ifges(7,5) * t833 + t817 * t843 - t818 * t857;
t837 = Ifges(6,4) * t861 + Ifges(6,2) * t860 + Ifges(6,6) * t891;
t838 = Ifges(6,1) * t861 + Ifges(6,4) * t860 + Ifges(6,5) * t891;
t848 = Ifges(5,4) * t878 + Ifges(5,2) * t877 + Ifges(5,6) * t891;
t849 = Ifges(5,1) * t878 + Ifges(5,4) * t877 + Ifges(5,5) * t891;
t945 = Ifges(5,5) * t854 + Ifges(5,6) * t853 + t878 * t848 - t877 * t849 + mrSges(5,1) * t811 - mrSges(5,2) * t812 + Ifges(6,5) * t835 + Ifges(6,6) * t834 + t861 * t837 - t860 * t838 + mrSges(6,1) * t804 - mrSges(6,2) * t805 + t909 * t791 + t913 * t790 + pkin(5) * t799 + pkin(11) * t789 + pkin(4) * t780 + (Ifges(5,3) + Ifges(6,3)) * t880;
t882 = -mrSges(4,2) * t898 + mrSges(4,3) * t930;
t883 = (-mrSges(4,1) * t915 + mrSges(4,2) * t911) * t933;
t788 = t913 * t797 + t909 * t798;
t787 = m(6) * t813 - t834 * mrSges(6,1) + mrSges(6,2) * t835 - t860 * t845 + t846 * t861 + t788;
t867 = -mrSges(5,2) * t891 + mrSges(5,3) * t877;
t869 = mrSges(5,1) * t891 - mrSges(5,3) * t878;
t919 = -m(5) * t822 + t853 * mrSges(5,1) - mrSges(5,2) * t854 + t877 * t867 - t869 * t878 - t787;
t783 = m(4) * t828 + mrSges(4,1) * t897 - mrSges(4,3) * t885 + t882 * t898 - t883 * t931 + t919;
t941 = t783 * t915;
t881 = mrSges(4,1) * t898 - mrSges(4,3) * t931;
t864 = -mrSges(5,1) * t877 + mrSges(5,2) * t878;
t778 = m(5) * t811 + mrSges(5,1) * t880 - mrSges(5,3) * t854 - t864 * t878 + t867 * t891 + t780;
t927 = t905 * t786 - t793 * t901;
t779 = m(5) * t812 - mrSges(5,2) * t880 + mrSges(5,3) * t853 + t864 * t877 - t869 * t891 + t927;
t928 = -t778 * t910 + t914 * t779;
t769 = m(4) * t829 - mrSges(4,2) * t897 - mrSges(4,3) * t886 - t881 * t898 + t883 * t930 + t928;
t772 = t914 * t778 + t910 * t779;
t842 = -t903 * t862 + t875;
t771 = m(4) * t842 + mrSges(4,1) * t886 + mrSges(4,2) * t885 + (t881 * t911 - t882 * t915) * t933 + t772;
t758 = t769 * t937 - t771 * t903 + t907 * t941;
t754 = m(3) * t865 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t917 + t758;
t757 = t769 * t940 + t907 * t771 + t903 * t941;
t756 = m(3) * t879 + t757;
t765 = t915 * t769 - t783 * t911;
t764 = m(3) * t866 - mrSges(3,1) * t917 - qJDD(2) * mrSges(3,2) + t765;
t744 = t754 * t935 - t756 * t904 + t764 * t936;
t742 = m(2) * t892 + t744;
t750 = -t754 * t912 + t916 * t764;
t749 = m(2) * t893 + t750;
t934 = t906 * t742 + t902 * t749;
t743 = t754 * t938 + t908 * t756 + t764 * t939;
t929 = -t742 * t902 + t906 * t749;
t926 = m(2) * t900 + t743;
t836 = Ifges(6,5) * t861 + Ifges(6,6) * t860 + Ifges(6,3) * t891;
t773 = mrSges(6,2) * t813 - mrSges(6,3) * t804 + Ifges(6,1) * t835 + Ifges(6,4) * t834 + Ifges(6,5) * t880 - pkin(11) * t788 - t790 * t909 + t791 * t913 + t836 * t860 - t837 * t891;
t920 = mrSges(7,1) * t800 - mrSges(7,2) * t801 + Ifges(7,5) * t816 + Ifges(7,6) * t815 + Ifges(7,3) * t833 + t818 * t844 - t819 * t843;
t774 = -mrSges(6,1) * t813 + mrSges(6,3) * t805 + Ifges(6,4) * t835 + Ifges(6,2) * t834 + Ifges(6,6) * t880 - pkin(5) * t788 - t836 * t861 + t838 * t891 - t920;
t847 = Ifges(5,5) * t878 + Ifges(5,6) * t877 + Ifges(5,3) * t891;
t759 = -mrSges(5,1) * t822 + mrSges(5,3) * t812 + Ifges(5,4) * t854 + Ifges(5,2) * t853 + Ifges(5,6) * t880 - pkin(4) * t787 + qJ(5) * t927 + t901 * t773 + t905 * t774 - t878 * t847 + t891 * t849;
t760 = mrSges(5,2) * t822 - mrSges(5,3) * t811 + Ifges(5,1) * t854 + Ifges(5,4) * t853 + Ifges(5,5) * t880 - qJ(5) * t780 + t773 * t905 - t774 * t901 + t847 * t877 - t848 * t891;
t872 = Ifges(4,6) * t898 + (Ifges(4,4) * t911 + Ifges(4,2) * t915) * t933;
t873 = Ifges(4,5) * t898 + (Ifges(4,1) * t911 + Ifges(4,4) * t915) * t933;
t745 = Ifges(4,5) * t885 - Ifges(4,6) * t886 + Ifges(4,3) * t897 + mrSges(4,1) * t828 - mrSges(4,2) * t829 + t910 * t760 + t914 * t759 + pkin(3) * t919 + pkin(10) * t928 + (t872 * t911 - t873 * t915) * t933;
t871 = Ifges(4,3) * t898 + (Ifges(4,5) * t911 + Ifges(4,6) * t915) * t933;
t746 = mrSges(4,2) * t842 - mrSges(4,3) * t828 + Ifges(4,1) * t885 - Ifges(4,4) * t886 + Ifges(4,5) * t897 - pkin(10) * t772 - t759 * t910 + t760 * t914 + t871 * t930 - t872 * t898;
t751 = -mrSges(4,1) * t842 + mrSges(4,3) * t829 + Ifges(4,4) * t885 - Ifges(4,2) * t886 + Ifges(4,6) * t897 - pkin(3) * t772 - t871 * t931 + t898 * t873 - t945;
t921 = pkin(9) * t765 + t746 * t911 + t751 * t915;
t739 = -mrSges(3,1) * t879 + mrSges(3,3) * t866 + Ifges(3,5) * t917 + Ifges(3,6) * qJDD(2) - pkin(2) * t757 - t745 * t903 + t921 * t907;
t740 = mrSges(3,2) * t879 - mrSges(3,3) * t865 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t917 + t746 * t915 - t751 * t911 + (-t757 * t903 - t758 * t907) * pkin(9);
t922 = pkin(8) * t750 + t739 * t916 + t740 * t912;
t738 = mrSges(3,1) * t865 - mrSges(3,2) * t866 + Ifges(3,3) * qJDD(2) + pkin(2) * t758 + t745 * t907 + t921 * t903;
t737 = mrSges(2,2) * t900 - mrSges(2,3) * t892 - t739 * t912 + t740 * t916 + (-t743 * t904 - t744 * t908) * pkin(8);
t736 = -mrSges(2,1) * t900 + mrSges(2,3) * t893 - pkin(1) * t743 - t738 * t904 + t922 * t908;
t1 = [-m(1) * g(1) + t929; -m(1) * g(2) + t934; -m(1) * g(3) + t926; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t934 - t902 * t736 + t906 * t737; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t929 + t906 * t736 + t902 * t737; -mrSges(1,1) * g(2) + mrSges(2,1) * t892 + mrSges(1,2) * g(1) - mrSges(2,2) * t893 + pkin(1) * t744 + t738 * t908 + t904 * t922; t926; t738; t745; t945; t787; t920;];
tauJB  = t1;
