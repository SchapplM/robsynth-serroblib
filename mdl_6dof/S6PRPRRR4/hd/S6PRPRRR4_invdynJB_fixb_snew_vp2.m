% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 01:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:58:21
% EndTime: 2019-05-05 00:58:41
% DurationCPUTime: 19.86s
% Computational Cost: add. (328302->322), mult. (731489->412), div. (0->0), fcn. (559996->14), ass. (0->151)
t892 = qJD(2) ^ 2;
t878 = sin(pkin(11));
t881 = cos(pkin(11));
t863 = g(1) * t878 - g(2) * t881;
t864 = -g(1) * t881 - g(2) * t878;
t876 = -g(3) + qJDD(1);
t879 = sin(pkin(6));
t882 = cos(pkin(6));
t886 = sin(qJ(2));
t890 = cos(qJ(2));
t835 = -t886 * t864 + (t863 * t882 + t876 * t879) * t890;
t880 = cos(pkin(12));
t925 = pkin(3) * t880;
t877 = sin(pkin(12));
t924 = mrSges(4,2) * t877;
t898 = qJDD(3) - t835;
t826 = -qJDD(2) * pkin(2) - t892 * qJ(3) + t898;
t874 = t877 ^ 2;
t920 = t882 * t886;
t921 = t879 * t886;
t836 = t863 * t920 + t890 * t864 + t876 * t921;
t830 = -pkin(2) * t892 + qJDD(2) * qJ(3) + t836;
t852 = -t863 * t879 + t876 * t882;
t914 = qJD(2) * qJD(3);
t918 = t880 * t852 - 0.2e1 * t877 * t914;
t807 = (-pkin(8) * qJDD(2) + t892 * t925 - t830) * t877 + t918;
t810 = t877 * t852 + (t830 + 0.2e1 * t914) * t880;
t912 = qJDD(2) * t880;
t875 = t880 ^ 2;
t922 = t875 * t892;
t808 = -pkin(3) * t922 + pkin(8) * t912 + t810;
t885 = sin(qJ(4));
t889 = cos(qJ(4));
t794 = t885 * t807 + t889 * t808;
t916 = qJD(2) * t880;
t917 = qJD(2) * t877;
t856 = -t885 * t917 + t889 * t916;
t902 = t877 * t889 + t880 * t885;
t857 = t902 * qJD(2);
t842 = -pkin(4) * t856 - pkin(9) * t857;
t891 = qJD(4) ^ 2;
t787 = -pkin(4) * t891 + qJDD(4) * pkin(9) + t842 * t856 + t794;
t820 = (-pkin(2) - t925) * qJDD(2) + (-qJ(3) + (-t874 - t875) * pkin(8)) * t892 + t898;
t854 = t857 * qJD(4);
t913 = qJDD(2) * t877;
t843 = -t885 * t913 + t889 * t912 - t854;
t915 = t856 * qJD(4);
t844 = qJDD(2) * t902 + t915;
t798 = (-t844 - t915) * pkin(9) + (-t843 + t854) * pkin(4) + t820;
t884 = sin(qJ(5));
t888 = cos(qJ(5));
t782 = -t884 * t787 + t888 * t798;
t846 = qJD(4) * t888 - t857 * t884;
t819 = qJD(5) * t846 + qJDD(4) * t884 + t844 * t888;
t841 = qJDD(5) - t843;
t847 = qJD(4) * t884 + t857 * t888;
t855 = qJD(5) - t856;
t780 = (t846 * t855 - t819) * pkin(10) + (t846 * t847 + t841) * pkin(5) + t782;
t783 = t888 * t787 + t884 * t798;
t818 = -qJD(5) * t847 + qJDD(4) * t888 - t844 * t884;
t829 = pkin(5) * t855 - pkin(10) * t847;
t845 = t846 ^ 2;
t781 = -pkin(5) * t845 + pkin(10) * t818 - t829 * t855 + t783;
t883 = sin(qJ(6));
t887 = cos(qJ(6));
t778 = t780 * t887 - t781 * t883;
t821 = t846 * t887 - t847 * t883;
t792 = qJD(6) * t821 + t818 * t883 + t819 * t887;
t822 = t846 * t883 + t847 * t887;
t803 = -mrSges(7,1) * t821 + mrSges(7,2) * t822;
t853 = qJD(6) + t855;
t811 = -mrSges(7,2) * t853 + mrSges(7,3) * t821;
t838 = qJDD(6) + t841;
t774 = m(7) * t778 + mrSges(7,1) * t838 - t792 * mrSges(7,3) - t803 * t822 + t811 * t853;
t779 = t780 * t883 + t781 * t887;
t791 = -qJD(6) * t822 + t818 * t887 - t819 * t883;
t812 = mrSges(7,1) * t853 - mrSges(7,3) * t822;
t775 = m(7) * t779 - mrSges(7,2) * t838 + t791 * mrSges(7,3) + t803 * t821 - t812 * t853;
t766 = t887 * t774 + t883 * t775;
t823 = -mrSges(6,1) * t846 + mrSges(6,2) * t847;
t827 = -mrSges(6,2) * t855 + mrSges(6,3) * t846;
t764 = m(6) * t782 + mrSges(6,1) * t841 - mrSges(6,3) * t819 - t823 * t847 + t827 * t855 + t766;
t828 = mrSges(6,1) * t855 - mrSges(6,3) * t847;
t908 = -t774 * t883 + t887 * t775;
t765 = m(6) * t783 - mrSges(6,2) * t841 + mrSges(6,3) * t818 + t823 * t846 - t828 * t855 + t908;
t759 = t888 * t764 + t884 * t765;
t850 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t856;
t851 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t857;
t896 = m(5) * t820 - t843 * mrSges(5,1) + mrSges(5,2) * t844 - t856 * t850 + t851 * t857 + t759;
t895 = -m(4) * t826 + mrSges(4,1) * t912 - t896 + (t874 * t892 + t922) * mrSges(4,3);
t754 = t895 + m(3) * t835 - mrSges(3,2) * t892 + (mrSges(3,1) - t924) * qJDD(2);
t923 = t754 * t890;
t760 = -t764 * t884 + t888 * t765;
t839 = -mrSges(5,1) * t856 + mrSges(5,2) * t857;
t757 = m(5) * t794 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t843 - qJD(4) * t851 + t839 * t856 + t760;
t793 = t807 * t889 - t885 * t808;
t786 = -qJDD(4) * pkin(4) - pkin(9) * t891 + t857 * t842 - t793;
t784 = -pkin(5) * t818 - pkin(10) * t845 + t829 * t847 + t786;
t899 = m(7) * t784 - t791 * mrSges(7,1) + t792 * mrSges(7,2) - t821 * t811 + t812 * t822;
t776 = -m(6) * t786 + t818 * mrSges(6,1) - mrSges(6,2) * t819 + t846 * t827 - t828 * t847 - t899;
t770 = m(5) * t793 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t844 + qJD(4) * t850 - t839 * t857 + t776;
t749 = t885 * t757 + t889 * t770;
t809 = -t830 * t877 + t918;
t901 = mrSges(4,3) * qJDD(2) + t892 * (-mrSges(4,1) * t880 + t924);
t747 = m(4) * t809 - t877 * t901 + t749;
t909 = t889 * t757 - t885 * t770;
t748 = m(4) * t810 + t880 * t901 + t909;
t910 = -t747 * t877 + t880 * t748;
t739 = m(3) * t836 - mrSges(3,1) * t892 - qJDD(2) * mrSges(3,2) + t910;
t742 = t880 * t747 + t877 * t748;
t741 = m(3) * t852 + t742;
t730 = t739 * t920 - t741 * t879 + t882 * t923;
t728 = m(2) * t863 + t730;
t734 = t890 * t739 - t754 * t886;
t733 = m(2) * t864 + t734;
t919 = t881 * t728 + t878 * t733;
t729 = t739 * t921 + t882 * t741 + t879 * t923;
t911 = -t728 * t878 + t881 * t733;
t907 = m(2) * t876 + t729;
t906 = Ifges(4,1) * t877 + Ifges(4,4) * t880;
t905 = Ifges(4,4) * t877 + Ifges(4,2) * t880;
t904 = Ifges(4,5) * t877 + Ifges(4,6) * t880;
t799 = Ifges(7,5) * t822 + Ifges(7,6) * t821 + Ifges(7,3) * t853;
t801 = Ifges(7,1) * t822 + Ifges(7,4) * t821 + Ifges(7,5) * t853;
t767 = -mrSges(7,1) * t784 + mrSges(7,3) * t779 + Ifges(7,4) * t792 + Ifges(7,2) * t791 + Ifges(7,6) * t838 - t799 * t822 + t801 * t853;
t800 = Ifges(7,4) * t822 + Ifges(7,2) * t821 + Ifges(7,6) * t853;
t768 = mrSges(7,2) * t784 - mrSges(7,3) * t778 + Ifges(7,1) * t792 + Ifges(7,4) * t791 + Ifges(7,5) * t838 + t799 * t821 - t800 * t853;
t813 = Ifges(6,5) * t847 + Ifges(6,6) * t846 + Ifges(6,3) * t855;
t815 = Ifges(6,1) * t847 + Ifges(6,4) * t846 + Ifges(6,5) * t855;
t750 = -mrSges(6,1) * t786 + mrSges(6,3) * t783 + Ifges(6,4) * t819 + Ifges(6,2) * t818 + Ifges(6,6) * t841 - pkin(5) * t899 + pkin(10) * t908 + t887 * t767 + t883 * t768 - t847 * t813 + t855 * t815;
t814 = Ifges(6,4) * t847 + Ifges(6,2) * t846 + Ifges(6,6) * t855;
t751 = mrSges(6,2) * t786 - mrSges(6,3) * t782 + Ifges(6,1) * t819 + Ifges(6,4) * t818 + Ifges(6,5) * t841 - pkin(10) * t766 - t767 * t883 + t768 * t887 + t813 * t846 - t814 * t855;
t831 = Ifges(5,5) * t857 + Ifges(5,6) * t856 + Ifges(5,3) * qJD(4);
t832 = Ifges(5,4) * t857 + Ifges(5,2) * t856 + Ifges(5,6) * qJD(4);
t735 = mrSges(5,2) * t820 - mrSges(5,3) * t793 + Ifges(5,1) * t844 + Ifges(5,4) * t843 + Ifges(5,5) * qJDD(4) - pkin(9) * t759 - qJD(4) * t832 - t750 * t884 + t751 * t888 + t831 * t856;
t833 = Ifges(5,1) * t857 + Ifges(5,4) * t856 + Ifges(5,5) * qJD(4);
t897 = -mrSges(7,1) * t778 + mrSges(7,2) * t779 - Ifges(7,5) * t792 - Ifges(7,6) * t791 - Ifges(7,3) * t838 - t822 * t800 + t821 * t801;
t893 = mrSges(6,1) * t782 - mrSges(6,2) * t783 + Ifges(6,5) * t819 + Ifges(6,6) * t818 + Ifges(6,3) * t841 + pkin(5) * t766 + t847 * t814 - t846 * t815 - t897;
t743 = -mrSges(5,1) * t820 + mrSges(5,3) * t794 + Ifges(5,4) * t844 + Ifges(5,2) * t843 + Ifges(5,6) * qJDD(4) - pkin(4) * t759 + qJD(4) * t833 - t857 * t831 - t893;
t862 = t904 * qJD(2);
t724 = -mrSges(4,1) * t826 + mrSges(4,3) * t810 - pkin(3) * t896 + pkin(8) * t909 + qJDD(2) * t905 + t885 * t735 + t889 * t743 - t862 * t917;
t726 = mrSges(4,2) * t826 - mrSges(4,3) * t809 - pkin(8) * t749 + qJDD(2) * t906 + t735 * t889 - t743 * t885 + t862 * t916;
t723 = mrSges(3,2) * t852 - mrSges(3,3) * t835 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t892 - qJ(3) * t742 - t724 * t877 + t726 * t880;
t894 = mrSges(5,1) * t793 - mrSges(5,2) * t794 + Ifges(5,5) * t844 + Ifges(5,6) * t843 + Ifges(5,3) * qJDD(4) + pkin(4) * t776 + pkin(9) * t760 + t888 * t750 + t884 * t751 + t857 * t832 - t856 * t833;
t725 = (Ifges(3,6) - t904) * qJDD(2) - mrSges(3,1) * t852 + mrSges(3,3) * t836 - mrSges(4,1) * t809 + mrSges(4,2) * t810 - t894 - pkin(3) * t749 - pkin(2) * t742 + (-t877 * t905 + t880 * t906 + Ifges(3,5)) * t892;
t900 = pkin(7) * t734 + t723 * t886 + t725 * t890;
t758 = mrSges(4,2) * t913 - t895;
t722 = mrSges(3,1) * t835 - mrSges(3,2) * t836 + Ifges(3,3) * qJDD(2) - pkin(2) * t758 + qJ(3) * t910 + t880 * t724 + t877 * t726;
t721 = mrSges(2,2) * t876 - mrSges(2,3) * t863 + t890 * t723 - t886 * t725 + (-t729 * t879 - t730 * t882) * pkin(7);
t720 = -mrSges(2,1) * t876 + mrSges(2,3) * t864 - pkin(1) * t729 - t722 * t879 + t882 * t900;
t1 = [-m(1) * g(1) + t911; -m(1) * g(2) + t919; -m(1) * g(3) + t907; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t919 - t878 * t720 + t881 * t721; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t911 + t881 * t720 + t878 * t721; -mrSges(1,1) * g(2) + mrSges(2,1) * t863 + mrSges(1,2) * g(1) - mrSges(2,2) * t864 + pkin(1) * t730 + t722 * t882 + t879 * t900; t907; t722; t758; t894; t893; -t897;];
tauJB  = t1;
