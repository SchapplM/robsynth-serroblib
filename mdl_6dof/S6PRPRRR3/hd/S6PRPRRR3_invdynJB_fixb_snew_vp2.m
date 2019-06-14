% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRR3
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
% Datum: 2019-05-05 00:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:41:07
% EndTime: 2019-05-05 00:41:28
% DurationCPUTime: 21.19s
% Computational Cost: add. (359629->321), mult. (812234->410), div. (0->0), fcn. (638122->14), ass. (0->148)
t892 = qJD(2) ^ 2;
t879 = sin(pkin(11));
t882 = cos(pkin(11));
t862 = g(1) * t879 - g(2) * t882;
t863 = -g(1) * t882 - g(2) * t879;
t877 = -g(3) + qJDD(1);
t880 = sin(pkin(6));
t883 = cos(pkin(6));
t887 = sin(qJ(2));
t891 = cos(qJ(2));
t840 = -t887 * t863 + (t862 * t883 + t877 * t880) * t891;
t881 = cos(pkin(12));
t926 = pkin(3) * t881;
t878 = sin(pkin(12));
t925 = mrSges(4,2) * t878;
t898 = qJDD(3) - t840;
t830 = -qJDD(2) * pkin(2) - t892 * qJ(3) + t898;
t874 = t878 ^ 2;
t875 = t881 ^ 2;
t826 = (-pkin(2) - t926) * qJDD(2) + (-qJ(3) + (-t874 - t875) * pkin(8)) * t892 + t898;
t886 = sin(qJ(4));
t890 = cos(qJ(4));
t904 = t878 * t890 + t881 * t886;
t856 = t904 * qJD(2);
t903 = -t878 * t886 + t881 * t890;
t846 = -t856 * qJD(4) + t903 * qJDD(2);
t855 = t903 * qJD(2);
t917 = t855 * qJD(4);
t847 = t904 * qJDD(2) + t917;
t850 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t855;
t851 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t856;
t921 = t883 * t887;
t922 = t880 * t887;
t841 = t862 * t921 + t891 * t863 + t877 * t922;
t833 = -pkin(2) * t892 + qJDD(2) * qJ(3) + t841;
t853 = -t862 * t880 + t877 * t883;
t916 = qJD(2) * qJD(3);
t919 = t881 * t853 - 0.2e1 * t878 * t916;
t813 = (-pkin(8) * qJDD(2) + t892 * t926 - t833) * t878 + t919;
t825 = t878 * t853 + (t833 + 0.2e1 * t916) * t881;
t915 = qJDD(2) * t881;
t923 = t875 * t892;
t816 = -pkin(3) * t923 + pkin(8) * t915 + t825;
t790 = t890 * t813 - t886 * t816;
t787 = (-t847 + t917) * pkin(9) + (t855 * t856 + qJDD(4)) * pkin(4) + t790;
t791 = t886 * t813 + t890 * t816;
t852 = qJD(4) * pkin(4) - pkin(9) * t856;
t854 = t855 ^ 2;
t789 = -pkin(4) * t854 + pkin(9) * t846 - qJD(4) * t852 + t791;
t885 = sin(qJ(5));
t889 = cos(qJ(5));
t784 = t885 * t787 + t889 * t789;
t838 = t855 * t889 - t856 * t885;
t839 = t855 * t885 + t856 * t889;
t823 = -pkin(5) * t838 - pkin(10) * t839;
t876 = qJD(4) + qJD(5);
t872 = t876 ^ 2;
t873 = qJDD(4) + qJDD(5);
t781 = -pkin(5) * t872 + pkin(10) * t873 + t823 * t838 + t784;
t796 = -t846 * pkin(4) - t854 * pkin(9) + t856 * t852 + t826;
t807 = -qJD(5) * t839 + t846 * t889 - t847 * t885;
t808 = qJD(5) * t838 + t846 * t885 + t847 * t889;
t785 = (-t838 * t876 - t808) * pkin(10) + (t839 * t876 - t807) * pkin(5) + t796;
t884 = sin(qJ(6));
t888 = cos(qJ(6));
t778 = -t781 * t884 + t785 * t888;
t827 = -t839 * t884 + t876 * t888;
t794 = qJD(6) * t827 + t808 * t888 + t873 * t884;
t806 = qJDD(6) - t807;
t828 = t839 * t888 + t876 * t884;
t809 = -mrSges(7,1) * t827 + mrSges(7,2) * t828;
t834 = qJD(6) - t838;
t814 = -mrSges(7,2) * t834 + mrSges(7,3) * t827;
t774 = m(7) * t778 + mrSges(7,1) * t806 - t794 * mrSges(7,3) - t809 * t828 + t814 * t834;
t779 = t781 * t888 + t785 * t884;
t793 = -qJD(6) * t828 - t808 * t884 + t873 * t888;
t815 = mrSges(7,1) * t834 - mrSges(7,3) * t828;
t775 = m(7) * t779 - mrSges(7,2) * t806 + t793 * mrSges(7,3) + t809 * t827 - t815 * t834;
t763 = t888 * t774 + t884 * t775;
t831 = -mrSges(6,2) * t876 + mrSges(6,3) * t838;
t832 = mrSges(6,1) * t876 - mrSges(6,3) * t839;
t901 = m(6) * t796 - t807 * mrSges(6,1) + t808 * mrSges(6,2) - t838 * t831 + t839 * t832 + t763;
t896 = m(5) * t826 - t846 * mrSges(5,1) + t847 * mrSges(5,2) - t855 * t850 + t856 * t851 + t901;
t894 = -m(4) * t830 + mrSges(4,1) * t915 - t896 + (t874 * t892 + t923) * mrSges(4,3);
t757 = t894 + (mrSges(3,1) - t925) * qJDD(2) - t892 * mrSges(3,2) + m(3) * t840;
t924 = t757 * t891;
t822 = -mrSges(6,1) * t838 + mrSges(6,2) * t839;
t910 = -t774 * t884 + t888 * t775;
t760 = m(6) * t784 - mrSges(6,2) * t873 + mrSges(6,3) * t807 + t822 * t838 - t832 * t876 + t910;
t783 = t787 * t889 - t789 * t885;
t780 = -pkin(5) * t873 - pkin(10) * t872 + t823 * t839 - t783;
t899 = -m(7) * t780 + t793 * mrSges(7,1) - t794 * mrSges(7,2) + t827 * t814 - t815 * t828;
t770 = m(6) * t783 + mrSges(6,1) * t873 - mrSges(6,3) * t808 - t822 * t839 + t831 * t876 + t899;
t754 = t885 * t760 + t889 * t770;
t844 = -mrSges(5,1) * t855 + mrSges(5,2) * t856;
t752 = m(5) * t790 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t847 + qJD(4) * t850 - t844 * t856 + t754;
t911 = t889 * t760 - t770 * t885;
t753 = m(5) * t791 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t846 - qJD(4) * t851 + t844 * t855 + t911;
t746 = t890 * t752 + t886 * t753;
t824 = -t833 * t878 + t919;
t902 = mrSges(4,3) * qJDD(2) + t892 * (-mrSges(4,1) * t881 + t925);
t744 = m(4) * t824 - t902 * t878 + t746;
t912 = -t886 * t752 + t890 * t753;
t745 = m(4) * t825 + t902 * t881 + t912;
t913 = -t744 * t878 + t881 * t745;
t736 = m(3) * t841 - mrSges(3,1) * t892 - qJDD(2) * mrSges(3,2) + t913;
t739 = t881 * t744 + t878 * t745;
t738 = m(3) * t853 + t739;
t727 = t736 * t921 - t738 * t880 + t883 * t924;
t725 = m(2) * t862 + t727;
t731 = t891 * t736 - t757 * t887;
t730 = m(2) * t863 + t731;
t920 = t882 * t725 + t879 * t730;
t906 = Ifges(4,5) * t878 + Ifges(4,6) * t881;
t918 = t892 * t906;
t726 = t736 * t922 + t883 * t738 + t880 * t924;
t914 = -t725 * t879 + t882 * t730;
t909 = m(2) * t877 + t726;
t908 = Ifges(4,1) * t878 + Ifges(4,4) * t881;
t907 = Ifges(4,4) * t878 + Ifges(4,2) * t881;
t797 = Ifges(7,5) * t828 + Ifges(7,6) * t827 + Ifges(7,3) * t834;
t799 = Ifges(7,1) * t828 + Ifges(7,4) * t827 + Ifges(7,5) * t834;
t767 = -mrSges(7,1) * t780 + mrSges(7,3) * t779 + Ifges(7,4) * t794 + Ifges(7,2) * t793 + Ifges(7,6) * t806 - t797 * t828 + t799 * t834;
t798 = Ifges(7,4) * t828 + Ifges(7,2) * t827 + Ifges(7,6) * t834;
t768 = mrSges(7,2) * t780 - mrSges(7,3) * t778 + Ifges(7,1) * t794 + Ifges(7,4) * t793 + Ifges(7,5) * t806 + t797 * t827 - t798 * t834;
t817 = Ifges(6,5) * t839 + Ifges(6,6) * t838 + Ifges(6,3) * t876;
t818 = Ifges(6,4) * t839 + Ifges(6,2) * t838 + Ifges(6,6) * t876;
t747 = mrSges(6,2) * t796 - mrSges(6,3) * t783 + Ifges(6,1) * t808 + Ifges(6,4) * t807 + Ifges(6,5) * t873 - pkin(10) * t763 - t767 * t884 + t768 * t888 + t817 * t838 - t818 * t876;
t819 = Ifges(6,1) * t839 + Ifges(6,4) * t838 + Ifges(6,5) * t876;
t895 = mrSges(7,1) * t778 - mrSges(7,2) * t779 + Ifges(7,5) * t794 + Ifges(7,6) * t793 + Ifges(7,3) * t806 + t798 * t828 - t799 * t827;
t748 = -mrSges(6,1) * t796 + mrSges(6,3) * t784 + Ifges(6,4) * t808 + Ifges(6,2) * t807 + Ifges(6,6) * t873 - pkin(5) * t763 - t817 * t839 + t819 * t876 - t895;
t835 = Ifges(5,5) * t856 + Ifges(5,6) * t855 + Ifges(5,3) * qJD(4);
t837 = Ifges(5,1) * t856 + Ifges(5,4) * t855 + Ifges(5,5) * qJD(4);
t732 = -mrSges(5,1) * t826 + mrSges(5,3) * t791 + Ifges(5,4) * t847 + Ifges(5,2) * t846 + Ifges(5,6) * qJDD(4) - pkin(4) * t901 + pkin(9) * t911 + qJD(4) * t837 + t885 * t747 + t889 * t748 - t856 * t835;
t836 = Ifges(5,4) * t856 + Ifges(5,2) * t855 + Ifges(5,6) * qJD(4);
t740 = mrSges(5,2) * t826 - mrSges(5,3) * t790 + Ifges(5,1) * t847 + Ifges(5,4) * t846 + Ifges(5,5) * qJDD(4) - pkin(9) * t754 - qJD(4) * t836 + t747 * t889 - t748 * t885 + t835 * t855;
t721 = -mrSges(4,1) * t830 + mrSges(4,3) * t825 - pkin(3) * t896 + pkin(8) * t912 + t907 * qJDD(2) + t890 * t732 + t886 * t740 - t878 * t918;
t722 = mrSges(4,2) * t830 - mrSges(4,3) * t824 - pkin(8) * t746 + t908 * qJDD(2) - t886 * t732 + t890 * t740 + t881 * t918;
t720 = mrSges(3,2) * t853 - mrSges(3,3) * t840 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t892 - qJ(3) * t739 - t721 * t878 + t722 * t881;
t897 = -mrSges(6,1) * t783 + mrSges(6,2) * t784 - Ifges(6,5) * t808 - Ifges(6,6) * t807 - Ifges(6,3) * t873 - pkin(5) * t899 - pkin(10) * t910 - t888 * t767 - t884 * t768 - t839 * t818 + t838 * t819;
t893 = mrSges(5,1) * t790 - mrSges(5,2) * t791 + Ifges(5,5) * t847 + Ifges(5,6) * t846 + Ifges(5,3) * qJDD(4) + pkin(4) * t754 + t856 * t836 - t855 * t837 - t897;
t723 = -t893 + (Ifges(3,6) - t906) * qJDD(2) - mrSges(3,1) * t853 + mrSges(3,3) * t841 - mrSges(4,1) * t824 + mrSges(4,2) * t825 - pkin(3) * t746 - pkin(2) * t739 + (-t878 * t907 + t881 * t908 + Ifges(3,5)) * t892;
t900 = pkin(7) * t731 + t720 * t887 + t723 * t891;
t761 = qJDD(2) * t925 - t894;
t719 = mrSges(3,1) * t840 - mrSges(3,2) * t841 + Ifges(3,3) * qJDD(2) - pkin(2) * t761 + qJ(3) * t913 + t881 * t721 + t878 * t722;
t718 = mrSges(2,2) * t877 - mrSges(2,3) * t862 + t891 * t720 - t887 * t723 + (-t726 * t880 - t727 * t883) * pkin(7);
t717 = -mrSges(2,1) * t877 + mrSges(2,3) * t863 - pkin(1) * t726 - t880 * t719 + t900 * t883;
t1 = [-m(1) * g(1) + t914; -m(1) * g(2) + t920; -m(1) * g(3) + t909; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t920 - t879 * t717 + t882 * t718; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t914 + t882 * t717 + t879 * t718; -mrSges(1,1) * g(2) + mrSges(2,1) * t862 + mrSges(1,2) * g(1) - mrSges(2,2) * t863 + pkin(1) * t727 + t883 * t719 + t900 * t880; t909; t719; t761; t893; -t897; t895;];
tauJB  = t1;
