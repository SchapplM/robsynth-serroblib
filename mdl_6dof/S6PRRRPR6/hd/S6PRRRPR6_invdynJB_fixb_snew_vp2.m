% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 08:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:22:56
% EndTime: 2019-05-05 08:23:05
% DurationCPUTime: 8.34s
% Computational Cost: add. (137722->321), mult. (262678->395), div. (0->0), fcn. (178630->12), ass. (0->141)
t937 = Ifges(5,1) + Ifges(6,1);
t928 = Ifges(5,4) - Ifges(6,5);
t927 = Ifges(5,5) + Ifges(6,4);
t936 = -Ifges(5,2) - Ifges(6,3);
t926 = Ifges(5,6) - Ifges(6,6);
t935 = Ifges(5,3) + Ifges(6,2);
t887 = sin(qJ(4));
t888 = sin(qJ(3));
t915 = qJD(2) * t888;
t930 = cos(qJ(4));
t862 = -qJD(3) * t930 + t887 * t915;
t891 = cos(qJ(3));
t913 = qJD(2) * qJD(3);
t911 = t891 * t913;
t866 = qJDD(2) * t888 + t911;
t830 = -t862 * qJD(4) + t887 * qJDD(3) + t866 * t930;
t882 = sin(pkin(11));
t884 = cos(pkin(11));
t868 = g(1) * t882 - g(2) * t884;
t869 = -g(1) * t884 - g(2) * t882;
t881 = -g(3) + qJDD(1);
t892 = cos(qJ(2));
t885 = cos(pkin(6));
t889 = sin(qJ(2));
t921 = t885 * t889;
t883 = sin(pkin(6));
t922 = t883 * t889;
t818 = t868 * t921 + t892 * t869 + t881 * t922;
t894 = qJD(2) ^ 2;
t812 = -pkin(2) * t894 + qJDD(2) * pkin(8) + t818;
t845 = -t868 * t883 + t881 * t885;
t806 = -t888 * t812 + t891 * t845;
t865 = (-pkin(3) * t891 - pkin(9) * t888) * qJD(2);
t893 = qJD(3) ^ 2;
t900 = qJDD(3) * pkin(3) + t893 * pkin(9) - t865 * t915 + t806;
t914 = qJD(2) * t891;
t876 = qJD(4) - t914;
t923 = t862 * t876;
t934 = (-t830 + t923) * qJ(5) - t900;
t817 = -t889 * t869 + (t868 * t885 + t881 * t883) * t892;
t807 = t891 * t812 + t888 * t845;
t798 = -pkin(3) * t893 + qJDD(3) * pkin(9) + t865 * t914 + t807;
t811 = -qJDD(2) * pkin(2) - t894 * pkin(8) - t817;
t912 = t888 * t913;
t867 = t891 * qJDD(2) - t912;
t800 = (-t866 - t911) * pkin(9) + (-t867 + t912) * pkin(3) + t811;
t786 = -t887 * t798 + t800 * t930;
t863 = t887 * qJD(3) + t915 * t930;
t835 = pkin(4) * t862 - qJ(5) * t863;
t859 = qJDD(4) - t867;
t875 = t876 ^ 2;
t784 = -t859 * pkin(4) - t875 * qJ(5) + t863 * t835 + qJDD(5) - t786;
t778 = (-t830 - t923) * pkin(10) + (t862 * t863 - t859) * pkin(5) + t784;
t787 = t930 * t798 + t887 * t800;
t931 = 2 * qJD(5);
t783 = -pkin(4) * t875 + t859 * qJ(5) - t862 * t835 + t876 * t931 + t787;
t829 = t863 * qJD(4) - qJDD(3) * t930 + t887 * t866;
t844 = -pkin(5) * t876 - pkin(10) * t863;
t858 = t862 ^ 2;
t779 = -pkin(5) * t858 + pkin(10) * t829 + t844 * t876 + t783;
t886 = sin(qJ(6));
t890 = cos(qJ(6));
t777 = t778 * t886 + t779 * t890;
t781 = -t858 * pkin(10) + (-pkin(4) - pkin(5)) * t829 + (-pkin(4) * t876 + t844 + t931) * t863 - t934;
t832 = t862 * t886 + t863 * t890;
t793 = -qJD(6) * t832 + t829 * t890 - t830 * t886;
t831 = t862 * t890 - t863 * t886;
t794 = qJD(6) * t831 + t829 * t886 + t830 * t890;
t874 = qJD(6) - t876;
t801 = Ifges(7,5) * t832 + Ifges(7,6) * t831 + Ifges(7,3) * t874;
t803 = Ifges(7,1) * t832 + Ifges(7,4) * t831 + Ifges(7,5) * t874;
t855 = qJDD(6) - t859;
t766 = -mrSges(7,1) * t781 + mrSges(7,3) * t777 + Ifges(7,4) * t794 + Ifges(7,2) * t793 + Ifges(7,6) * t855 - t801 * t832 + t803 * t874;
t776 = t778 * t890 - t779 * t886;
t802 = Ifges(7,4) * t832 + Ifges(7,2) * t831 + Ifges(7,6) * t874;
t767 = mrSges(7,2) * t781 - mrSges(7,3) * t776 + Ifges(7,1) * t794 + Ifges(7,4) * t793 + Ifges(7,5) * t855 + t801 * t831 - t802 * t874;
t785 = -0.2e1 * qJD(5) * t863 + (t863 * t876 + t829) * pkin(4) + t934;
t842 = -mrSges(6,1) * t876 + mrSges(6,2) * t863;
t843 = -mrSges(6,2) * t862 + mrSges(6,3) * t876;
t813 = -mrSges(7,2) * t874 + mrSges(7,3) * t831;
t814 = mrSges(7,1) * t874 - mrSges(7,3) * t832;
t905 = -m(7) * t781 + t793 * mrSges(7,1) - t794 * mrSges(7,2) + t831 * t813 - t832 * t814;
t771 = m(6) * t785 + t829 * mrSges(6,1) - t830 * mrSges(6,3) - t863 * t842 + t862 * t843 + t905;
t808 = -mrSges(7,1) * t831 + mrSges(7,2) * t832;
t773 = m(7) * t776 + mrSges(7,1) * t855 - mrSges(7,3) * t794 - t808 * t832 + t813 * t874;
t774 = m(7) * t777 - mrSges(7,2) * t855 + mrSges(7,3) * t793 + t808 * t831 - t814 * t874;
t908 = -t886 * t773 + t890 * t774;
t917 = -t928 * t862 + t863 * t937 + t927 * t876;
t919 = t862 * t926 - t863 * t927 - t876 * t935;
t743 = mrSges(5,1) * t900 - mrSges(6,1) * t785 + mrSges(6,2) * t783 + mrSges(5,3) * t787 - pkin(4) * t771 - pkin(5) * t905 - pkin(10) * t908 - t890 * t766 - t886 * t767 + t829 * t936 + t928 * t830 + t926 * t859 + t919 * t863 + t917 * t876;
t765 = t890 * t773 + t886 * t774;
t918 = t862 * t936 + t863 * t928 + t876 * t926;
t744 = -mrSges(5,2) * t900 + mrSges(6,2) * t784 - mrSges(5,3) * t786 - mrSges(6,3) * t785 - pkin(10) * t765 - qJ(5) * t771 - t886 * t766 + t890 * t767 - t928 * t829 + t830 * t937 + t927 * t859 + t919 * t862 - t918 * t876;
t841 = mrSges(5,1) * t876 - mrSges(5,3) * t863;
t902 = m(6) * t783 + t859 * mrSges(6,3) + t876 * t842 + t908;
t836 = mrSges(6,1) * t862 - mrSges(6,3) * t863;
t916 = -mrSges(5,1) * t862 - mrSges(5,2) * t863 - t836;
t929 = -mrSges(5,3) - mrSges(6,2);
t761 = m(5) * t787 - t859 * mrSges(5,2) + t829 * t929 - t876 * t841 + t862 * t916 + t902;
t840 = -mrSges(5,2) * t876 - mrSges(5,3) * t862;
t899 = -m(6) * t784 + t859 * mrSges(6,1) + t876 * t843 - t765;
t762 = m(5) * t786 + t859 * mrSges(5,1) + t830 * t929 + t876 * t840 + t863 * t916 + t899;
t759 = t930 * t761 - t762 * t887;
t770 = m(5) * t900 - t829 * mrSges(5,1) - t830 * mrSges(5,2) - t862 * t840 - t863 * t841 - t771;
t850 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t888 + Ifges(4,2) * t891) * qJD(2);
t851 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t888 + Ifges(4,4) * t891) * qJD(2);
t933 = mrSges(4,1) * t806 - mrSges(4,2) * t807 + Ifges(4,5) * t866 + Ifges(4,6) * t867 + Ifges(4,3) * qJDD(3) + pkin(3) * t770 + pkin(9) * t759 + (t850 * t888 - t851 * t891) * qJD(2) + t743 * t930 + t887 * t744;
t764 = t830 * mrSges(6,2) + t863 * t836 - t899;
t898 = -mrSges(7,1) * t776 + mrSges(7,2) * t777 - Ifges(7,5) * t794 - Ifges(7,6) * t793 - Ifges(7,3) * t855 - t832 * t802 + t831 * t803;
t932 = -t829 * t926 + t830 * t927 + t935 * t859 + t862 * t917 + t863 * t918 + mrSges(5,1) * t786 - mrSges(6,1) * t784 - mrSges(5,2) * t787 + mrSges(6,3) * t783 - pkin(4) * t764 - pkin(5) * t765 + qJ(5) * (-t829 * mrSges(6,2) - t862 * t836 + t902) + t898;
t758 = t887 * t761 + t762 * t930;
t870 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t915;
t871 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t914;
t897 = -m(4) * t811 + t867 * mrSges(4,1) - t866 * mrSges(4,2) - t870 * t915 + t871 * t914 - t758;
t754 = m(3) * t817 + qJDD(2) * mrSges(3,1) - t894 * mrSges(3,2) + t897;
t924 = t754 * t892;
t864 = (-mrSges(4,1) * t891 + mrSges(4,2) * t888) * qJD(2);
t757 = m(4) * t807 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t867 - qJD(3) * t870 + t864 * t914 + t759;
t769 = m(4) * t806 + qJDD(3) * mrSges(4,1) - t866 * mrSges(4,3) + qJD(3) * t871 - t864 * t915 + t770;
t909 = t891 * t757 - t769 * t888;
t748 = m(3) * t818 - mrSges(3,1) * t894 - qJDD(2) * mrSges(3,2) + t909;
t751 = t888 * t757 + t891 * t769;
t750 = m(3) * t845 + t751;
t737 = t748 * t921 - t750 * t883 + t885 * t924;
t735 = m(2) * t868 + t737;
t742 = t892 * t748 - t754 * t889;
t741 = m(2) * t869 + t742;
t920 = t884 * t735 + t882 * t741;
t736 = t748 * t922 + t885 * t750 + t883 * t924;
t910 = -t735 * t882 + t884 * t741;
t906 = m(2) * t881 + t736;
t849 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t888 + Ifges(4,6) * t891) * qJD(2);
t733 = mrSges(4,2) * t811 - mrSges(4,3) * t806 + Ifges(4,1) * t866 + Ifges(4,4) * t867 + Ifges(4,5) * qJDD(3) - pkin(9) * t758 - qJD(3) * t850 - t887 * t743 + t744 * t930 + t849 * t914;
t738 = -mrSges(4,1) * t811 + mrSges(4,3) * t807 + Ifges(4,4) * t866 + Ifges(4,2) * t867 + Ifges(4,6) * qJDD(3) - pkin(3) * t758 + qJD(3) * t851 - t849 * t915 - t932;
t731 = mrSges(3,2) * t845 - mrSges(3,3) * t817 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t894 - pkin(8) * t751 + t733 * t891 - t738 * t888;
t732 = -mrSges(3,1) * t845 + mrSges(3,3) * t818 + t894 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t751 - t933;
t901 = pkin(7) * t742 + t731 * t889 + t732 * t892;
t730 = mrSges(3,1) * t817 - mrSges(3,2) * t818 + Ifges(3,3) * qJDD(2) + pkin(2) * t897 + pkin(8) * t909 + t888 * t733 + t891 * t738;
t729 = mrSges(2,2) * t881 - mrSges(2,3) * t868 + t892 * t731 - t889 * t732 + (-t736 * t883 - t737 * t885) * pkin(7);
t728 = -mrSges(2,1) * t881 + mrSges(2,3) * t869 - pkin(1) * t736 - t883 * t730 + t885 * t901;
t1 = [-m(1) * g(1) + t910; -m(1) * g(2) + t920; -m(1) * g(3) + t906; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t920 - t882 * t728 + t884 * t729; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t910 + t884 * t728 + t882 * t729; -mrSges(1,1) * g(2) + mrSges(2,1) * t868 + mrSges(1,2) * g(1) - mrSges(2,2) * t869 + pkin(1) * t737 + t885 * t730 + t883 * t901; t906; t730; t933; t932; t764; -t898;];
tauJB  = t1;
