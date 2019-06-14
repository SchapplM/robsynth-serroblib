% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:43:13
% EndTime: 2019-05-05 18:43:21
% DurationCPUTime: 7.28s
% Computational Cost: add. (87011->328), mult. (172672->392), div. (0->0), fcn. (99620->10), ass. (0->140)
t939 = -2 * qJD(4);
t938 = Ifges(4,1) + Ifges(5,2);
t929 = Ifges(4,4) + Ifges(5,6);
t928 = Ifges(4,5) - Ifges(5,4);
t937 = Ifges(4,2) + Ifges(5,3);
t927 = Ifges(4,6) - Ifges(5,5);
t936 = Ifges(4,3) + Ifges(5,1);
t893 = cos(qJ(3));
t889 = sin(qJ(3));
t918 = qJD(1) * qJD(3);
t916 = t889 * t918;
t858 = t893 * qJDD(1) - t916;
t875 = t889 * qJD(1);
t866 = pkin(4) * t875 - qJD(3) * pkin(8);
t883 = t893 ^ 2;
t896 = qJD(1) ^ 2;
t915 = t893 * t918;
t857 = t889 * qJDD(1) + t915;
t890 = sin(qJ(1));
t894 = cos(qJ(1));
t867 = t890 * g(1) - t894 * g(2);
t852 = qJDD(1) * pkin(1) + t867;
t868 = -t894 * g(1) - t890 * g(2);
t856 = -t896 * pkin(1) + t868;
t885 = sin(pkin(10));
t886 = cos(pkin(10));
t822 = t886 * t852 - t885 * t856;
t909 = -qJDD(1) * pkin(2) - t822;
t902 = pkin(3) * t916 + t875 * t939 + (-t857 - t915) * qJ(4) + t909;
t932 = -pkin(3) - pkin(8);
t781 = -t866 * t875 + (-pkin(4) * t883 - pkin(7)) * t896 + t932 * t858 + t902;
t884 = -g(3) + qJDD(2);
t823 = t885 * t852 + t886 * t856;
t809 = -t896 * pkin(2) + qJDD(1) * pkin(7) + t823;
t806 = t889 * t809;
t853 = (-pkin(3) * t893 - qJ(4) * t889) * qJD(1);
t895 = qJD(3) ^ 2;
t910 = -t895 * qJ(4) + t853 * t875 + qJDD(4) + t806;
t931 = pkin(8) * t896;
t792 = t857 * pkin(4) + t932 * qJDD(3) + (-pkin(4) * t918 - t889 * t931 - t884) * t893 + t910;
t888 = sin(qJ(5));
t892 = cos(qJ(5));
t776 = -t888 * t781 + t892 * t792;
t919 = qJD(1) * t893;
t850 = -t888 * qJD(3) - t892 * t919;
t818 = t850 * qJD(5) + t892 * qJDD(3) - t888 * t858;
t849 = qJDD(5) + t857;
t851 = t892 * qJD(3) - t888 * t919;
t871 = t875 + qJD(5);
t773 = (t850 * t871 - t818) * pkin(9) + (t850 * t851 + t849) * pkin(5) + t776;
t777 = t892 * t781 + t888 * t792;
t817 = -t851 * qJD(5) - t888 * qJDD(3) - t892 * t858;
t827 = t871 * pkin(5) - t851 * pkin(9);
t848 = t850 ^ 2;
t774 = -t848 * pkin(5) + t817 * pkin(9) - t871 * t827 + t777;
t887 = sin(qJ(6));
t891 = cos(qJ(6));
t772 = t887 * t773 + t891 * t774;
t803 = t893 * t809 + t889 * t884;
t798 = t895 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t939 - t853 * t919 - t803;
t791 = t858 * pkin(4) + qJD(3) * t866 - t883 * t931 - t798;
t779 = -t817 * pkin(5) - t848 * pkin(9) + t851 * t827 + t791;
t821 = t887 * t850 + t891 * t851;
t786 = -t821 * qJD(6) + t891 * t817 - t887 * t818;
t820 = t891 * t850 - t887 * t851;
t787 = t820 * qJD(6) + t887 * t817 + t891 * t818;
t869 = qJD(6) + t871;
t794 = Ifges(7,5) * t821 + Ifges(7,6) * t820 + Ifges(7,3) * t869;
t796 = Ifges(7,1) * t821 + Ifges(7,4) * t820 + Ifges(7,5) * t869;
t842 = qJDD(6) + t849;
t759 = -mrSges(7,1) * t779 + mrSges(7,3) * t772 + Ifges(7,4) * t787 + Ifges(7,2) * t786 + Ifges(7,6) * t842 - t821 * t794 + t869 * t796;
t771 = t891 * t773 - t887 * t774;
t795 = Ifges(7,4) * t821 + Ifges(7,2) * t820 + Ifges(7,6) * t869;
t760 = mrSges(7,2) * t779 - mrSges(7,3) * t771 + Ifges(7,1) * t787 + Ifges(7,4) * t786 + Ifges(7,5) * t842 + t820 * t794 - t869 * t795;
t810 = Ifges(6,5) * t851 + Ifges(6,6) * t850 + Ifges(6,3) * t871;
t812 = Ifges(6,1) * t851 + Ifges(6,4) * t850 + Ifges(6,5) * t871;
t804 = -t869 * mrSges(7,2) + t820 * mrSges(7,3);
t805 = t869 * mrSges(7,1) - t821 * mrSges(7,3);
t906 = m(7) * t779 - t786 * mrSges(7,1) + t787 * mrSges(7,2) - t820 * t804 + t821 * t805;
t800 = -t820 * mrSges(7,1) + t821 * mrSges(7,2);
t768 = m(7) * t771 + t842 * mrSges(7,1) - t787 * mrSges(7,3) - t821 * t800 + t869 * t804;
t769 = m(7) * t772 - t842 * mrSges(7,2) + t786 * mrSges(7,3) + t820 * t800 - t869 * t805;
t911 = -t887 * t768 + t891 * t769;
t743 = -mrSges(6,1) * t791 + mrSges(6,3) * t777 + Ifges(6,4) * t818 + Ifges(6,2) * t817 + Ifges(6,6) * t849 - pkin(5) * t906 + pkin(9) * t911 + t891 * t759 + t887 * t760 - t851 * t810 + t871 * t812;
t758 = t891 * t768 + t887 * t769;
t811 = Ifges(6,4) * t851 + Ifges(6,2) * t850 + Ifges(6,6) * t871;
t744 = mrSges(6,2) * t791 - mrSges(6,3) * t776 + Ifges(6,1) * t818 + Ifges(6,4) * t817 + Ifges(6,5) * t849 - pkin(9) * t758 - t887 * t759 + t891 * t760 + t850 * t810 - t871 * t811;
t854 = (mrSges(5,2) * t893 - mrSges(5,3) * t889) * qJD(1);
t864 = -mrSges(5,1) * t919 - qJD(3) * mrSges(5,3);
t824 = -t850 * mrSges(6,1) + t851 * mrSges(6,2);
t825 = -t871 * mrSges(6,2) + t850 * mrSges(6,3);
t755 = m(6) * t776 + t849 * mrSges(6,1) - t818 * mrSges(6,3) - t851 * t824 + t871 * t825 + t758;
t826 = t871 * mrSges(6,1) - t851 * mrSges(6,3);
t756 = m(6) * t777 - t849 * mrSges(6,2) + t817 * mrSges(6,3) + t850 * t824 - t871 * t826 + t911;
t752 = t892 * t755 + t888 * t756;
t925 = t893 * t884;
t799 = -qJDD(3) * pkin(3) + t910 - t925;
t905 = -m(5) * t799 - t857 * mrSges(5,1) - t752;
t751 = qJDD(3) * mrSges(5,2) + qJD(3) * t864 + t854 * t875 - t905;
t802 = -t806 + t925;
t865 = mrSges(5,1) * t875 + qJD(3) * mrSges(5,2);
t901 = -m(6) * t791 + t817 * mrSges(6,1) - t818 * mrSges(6,2) + t850 * t825 - t851 * t826 - t906;
t899 = -m(5) * t798 + qJDD(3) * mrSges(5,3) + qJD(3) * t865 + t854 * t919 - t901;
t920 = t928 * qJD(3) + (t938 * t889 + t929 * t893) * qJD(1);
t921 = t927 * qJD(3) + (t929 * t889 + t937 * t893) * qJD(1);
t935 = (t889 * t921 - t893 * t920) * qJD(1) + t936 * qJDD(3) + t928 * t857 + t927 * t858 + mrSges(4,1) * t802 - mrSges(4,2) * t803 + mrSges(5,2) * t799 - mrSges(5,3) * t798 - pkin(3) * t751 - pkin(8) * t752 + qJ(4) * (t858 * mrSges(5,1) + t899) - t888 * t743 + t892 * t744;
t930 = t896 * pkin(7);
t855 = (-mrSges(4,1) * t893 + mrSges(4,2) * t889) * qJD(1);
t863 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t919;
t749 = m(4) * t802 - t857 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t863 - t864) * qJD(3) + (-t854 - t855) * t875 + t905;
t862 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t875;
t763 = -qJDD(3) * mrSges(4,2) + t855 * t919 + t899 + (mrSges(4,3) + mrSges(5,1)) * t858 + m(4) * t803 - qJD(3) * t862;
t912 = -t889 * t749 + t893 * t763;
t739 = m(3) * t823 - t896 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t912;
t808 = t909 - t930;
t793 = -t858 * pkin(3) + t902 - t930;
t923 = -t888 * t755 + t892 * t756;
t908 = -m(5) * t793 - t858 * mrSges(5,2) + t865 * t875 - t923;
t898 = -m(4) * t808 + t863 * t919 + t858 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t857 + (-t862 * t889 - t864 * t893) * qJD(1) + t908;
t746 = m(3) * t822 + qJDD(1) * mrSges(3,1) - t896 * mrSges(3,2) + t898;
t736 = t885 * t739 + t886 * t746;
t733 = m(2) * t867 + qJDD(1) * mrSges(2,1) - t896 * mrSges(2,2) + t736;
t913 = t886 * t739 - t885 * t746;
t734 = m(2) * t868 - t896 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t913;
t924 = t894 * t733 + t890 * t734;
t742 = t893 * t749 + t889 * t763;
t922 = t936 * qJD(3) + (t928 * t889 + t927 * t893) * qJD(1);
t740 = m(3) * t884 + t742;
t914 = -t890 * t733 + t894 * t734;
t904 = mrSges(7,1) * t771 - mrSges(7,2) * t772 + Ifges(7,5) * t787 + Ifges(7,6) * t786 + Ifges(7,3) * t842 + t821 * t795 - t820 * t796;
t750 = -t857 * mrSges(5,3) + t864 * t919 - t908;
t727 = -mrSges(4,1) * t808 - mrSges(5,1) * t798 + mrSges(5,2) * t793 + mrSges(4,3) * t803 - pkin(3) * t750 - pkin(4) * t901 - pkin(8) * t923 + t920 * qJD(3) + t927 * qJDD(3) - t892 * t743 - t888 * t744 + t929 * t857 + t937 * t858 - t922 * t875;
t900 = mrSges(6,1) * t776 - mrSges(6,2) * t777 + Ifges(6,5) * t818 + Ifges(6,6) * t817 + Ifges(6,3) * t849 + pkin(5) * t758 + t851 * t811 - t850 * t812 + t904;
t729 = mrSges(5,1) * t799 + mrSges(4,2) * t808 - mrSges(4,3) * t802 - mrSges(5,3) * t793 + pkin(4) * t752 - qJ(4) * t750 - t921 * qJD(3) + t928 * qJDD(3) + t938 * t857 + t929 * t858 + t922 * t919 + t900;
t903 = mrSges(2,1) * t867 + mrSges(3,1) * t822 - mrSges(2,2) * t868 - mrSges(3,2) * t823 + pkin(1) * t736 + pkin(2) * t898 + pkin(7) * t912 + t893 * t727 + t889 * t729 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t725 = -mrSges(3,1) * t884 + mrSges(3,3) * t823 + t896 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t742 - t935;
t724 = mrSges(3,2) * t884 - mrSges(3,3) * t822 + Ifges(3,5) * qJDD(1) - t896 * Ifges(3,6) - pkin(7) * t742 - t889 * t727 + t893 * t729;
t723 = -mrSges(2,2) * g(3) - mrSges(2,3) * t867 + Ifges(2,5) * qJDD(1) - t896 * Ifges(2,6) - qJ(2) * t736 + t886 * t724 - t885 * t725;
t722 = mrSges(2,1) * g(3) + mrSges(2,3) * t868 + t896 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t740 + qJ(2) * t913 + t885 * t724 + t886 * t725;
t1 = [-m(1) * g(1) + t914; -m(1) * g(2) + t924; (-m(1) - m(2)) * g(3) + t740; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t924 - t890 * t722 + t894 * t723; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t914 + t894 * t722 + t890 * t723; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t903; t903; t740; t935; t751; t900; t904;];
tauJB  = t1;
