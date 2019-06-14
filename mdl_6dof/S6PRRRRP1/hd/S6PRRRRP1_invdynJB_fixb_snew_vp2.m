% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:28:35
% EndTime: 2019-05-05 09:28:47
% DurationCPUTime: 10.77s
% Computational Cost: add. (188387->321), mult. (368182->400), div. (0->0), fcn. (262618->12), ass. (0->139)
t923 = Ifges(6,1) + Ifges(7,1);
t917 = Ifges(6,4) + Ifges(7,4);
t916 = Ifges(6,5) + Ifges(7,5);
t922 = Ifges(6,2) + Ifges(7,2);
t915 = Ifges(6,6) + Ifges(7,6);
t921 = Ifges(6,3) + Ifges(7,3);
t875 = sin(qJ(4));
t876 = sin(qJ(3));
t879 = cos(qJ(4));
t880 = cos(qJ(3));
t847 = (t875 * t876 - t879 * t880) * qJD(2);
t870 = sin(pkin(11));
t872 = cos(pkin(11));
t857 = t870 * g(1) - t872 * g(2);
t858 = -t872 * g(1) - t870 * g(2);
t869 = -g(3) + qJDD(1);
t871 = sin(pkin(6));
t873 = cos(pkin(6));
t877 = sin(qJ(2));
t881 = cos(qJ(2));
t829 = -t877 * t858 + (t857 * t873 + t869 * t871) * t881;
t911 = t873 * t877;
t912 = t871 * t877;
t830 = t857 * t911 + t881 * t858 + t869 * t912;
t882 = qJD(2) ^ 2;
t825 = -t882 * pkin(2) + qJDD(2) * pkin(8) + t830;
t841 = -t871 * t857 + t873 * t869;
t804 = -t876 * t825 + t880 * t841;
t903 = qJD(2) * qJD(3);
t900 = t880 * t903;
t855 = t876 * qJDD(2) + t900;
t787 = (-t855 + t900) * pkin(9) + (t876 * t880 * t882 + qJDD(3)) * pkin(3) + t804;
t805 = t880 * t825 + t876 * t841;
t856 = t880 * qJDD(2) - t876 * t903;
t905 = qJD(2) * t876;
t862 = qJD(3) * pkin(3) - pkin(9) * t905;
t868 = t880 ^ 2;
t788 = -t868 * t882 * pkin(3) + t856 * pkin(9) - qJD(3) * t862 + t805;
t783 = t875 * t787 + t879 * t788;
t848 = (t875 * t880 + t876 * t879) * qJD(2);
t817 = -t848 * qJD(4) - t875 * t855 + t879 * t856;
t832 = t847 * mrSges(5,1) + t848 * mrSges(5,2);
t867 = qJD(3) + qJD(4);
t840 = t867 * mrSges(5,1) - t848 * mrSges(5,3);
t866 = qJDD(3) + qJDD(4);
t833 = t847 * pkin(4) - t848 * pkin(10);
t865 = t867 ^ 2;
t777 = -t865 * pkin(4) + t866 * pkin(10) - t847 * t833 + t783;
t889 = -qJDD(2) * pkin(2) - t829;
t795 = -t856 * pkin(3) + t862 * t905 + (-pkin(9) * t868 - pkin(8)) * t882 + t889;
t818 = -t847 * qJD(4) + t879 * t855 + t875 * t856;
t780 = (t847 * t867 - t818) * pkin(10) + (t848 * t867 - t817) * pkin(4) + t795;
t874 = sin(qJ(5));
t878 = cos(qJ(5));
t772 = -t874 * t777 + t878 * t780;
t835 = -t874 * t848 + t878 * t867;
t793 = t835 * qJD(5) + t878 * t818 + t874 * t866;
t836 = t878 * t848 + t874 * t867;
t807 = -t835 * mrSges(7,1) + t836 * mrSges(7,2);
t808 = -t835 * mrSges(6,1) + t836 * mrSges(6,2);
t815 = qJDD(5) - t817;
t842 = qJD(5) + t847;
t820 = -t842 * mrSges(6,2) + t835 * mrSges(6,3);
t769 = -0.2e1 * qJD(6) * t836 + (t835 * t842 - t793) * qJ(6) + (t835 * t836 + t815) * pkin(5) + t772;
t819 = -t842 * mrSges(7,2) + t835 * mrSges(7,3);
t902 = m(7) * t769 + t815 * mrSges(7,1) + t842 * t819;
t758 = m(6) * t772 + t815 * mrSges(6,1) + t842 * t820 + (-t807 - t808) * t836 + (-mrSges(6,3) - mrSges(7,3)) * t793 + t902;
t773 = t878 * t777 + t874 * t780;
t792 = -t836 * qJD(5) - t874 * t818 + t878 * t866;
t821 = t842 * pkin(5) - t836 * qJ(6);
t834 = t835 ^ 2;
t771 = -t834 * pkin(5) + t792 * qJ(6) + 0.2e1 * qJD(6) * t835 - t842 * t821 + t773;
t901 = m(7) * t771 + t792 * mrSges(7,3) + t835 * t807;
t822 = t842 * mrSges(7,1) - t836 * mrSges(7,3);
t906 = -t842 * mrSges(6,1) + t836 * mrSges(6,3) - t822;
t918 = -mrSges(6,2) - mrSges(7,2);
t761 = m(6) * t773 + t792 * mrSges(6,3) + t835 * t808 + t918 * t815 + t906 * t842 + t901;
t896 = -t874 * t758 + t878 * t761;
t751 = m(5) * t783 - t866 * mrSges(5,2) + t817 * mrSges(5,3) - t847 * t832 - t867 * t840 + t896;
t782 = t879 * t787 - t875 * t788;
t839 = -t867 * mrSges(5,2) - t847 * mrSges(5,3);
t776 = -t866 * pkin(4) - t865 * pkin(10) + t848 * t833 - t782;
t774 = -t792 * pkin(5) - t834 * qJ(6) + t836 * t821 + qJDD(6) + t776;
t894 = -m(7) * t774 + t792 * mrSges(7,1) + t835 * t819;
t885 = -m(6) * t776 + t792 * mrSges(6,1) + t918 * t793 + t835 * t820 + t906 * t836 + t894;
t763 = m(5) * t782 + t866 * mrSges(5,1) - t818 * mrSges(5,3) - t848 * t832 + t867 * t839 + t885;
t743 = t875 * t751 + t879 * t763;
t845 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t876 + Ifges(4,2) * t880) * qJD(2);
t846 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t876 + Ifges(4,4) * t880) * qJD(2);
t767 = t793 * mrSges(7,2) + t836 * t822 - t894;
t907 = t917 * t835 + t923 * t836 + t916 * t842;
t909 = -t915 * t835 - t916 * t836 - t921 * t842;
t745 = -mrSges(6,1) * t776 + mrSges(6,3) * t773 - mrSges(7,1) * t774 + mrSges(7,3) * t771 - pkin(5) * t767 + qJ(6) * t901 + (-qJ(6) * t822 + t907) * t842 + t909 * t836 + (-qJ(6) * mrSges(7,2) + t915) * t815 + t917 * t793 + t922 * t792;
t766 = -t793 * mrSges(7,3) - t836 * t807 + t902;
t908 = -t922 * t835 - t917 * t836 - t915 * t842;
t753 = mrSges(6,2) * t776 + mrSges(7,2) * t774 - mrSges(6,3) * t772 - mrSges(7,3) * t769 - qJ(6) * t766 + t917 * t792 + t923 * t793 + t916 * t815 - t909 * t835 + t908 * t842;
t827 = Ifges(5,4) * t848 - Ifges(5,2) * t847 + Ifges(5,6) * t867;
t828 = Ifges(5,1) * t848 - Ifges(5,4) * t847 + Ifges(5,5) * t867;
t886 = -mrSges(5,1) * t782 + mrSges(5,2) * t783 - Ifges(5,5) * t818 - Ifges(5,6) * t817 - Ifges(5,3) * t866 - pkin(4) * t885 - pkin(10) * t896 - t878 * t745 - t874 * t753 - t848 * t827 - t847 * t828;
t920 = mrSges(4,1) * t804 - mrSges(4,2) * t805 + Ifges(4,5) * t855 + Ifges(4,6) * t856 + Ifges(4,3) * qJDD(3) + pkin(3) * t743 + (t876 * t845 - t880 * t846) * qJD(2) - t886;
t919 = mrSges(6,1) * t772 + mrSges(7,1) * t769 - mrSges(6,2) * t773 - mrSges(7,2) * t771 + pkin(5) * t766 + t915 * t792 + t916 * t793 + t921 * t815 - t907 * t835 - t908 * t836;
t824 = -t882 * pkin(8) + t889;
t859 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t905;
t904 = qJD(2) * t880;
t860 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t904;
t755 = t878 * t758 + t874 * t761;
t888 = m(5) * t795 - t817 * mrSges(5,1) + t818 * mrSges(5,2) + t847 * t839 + t848 * t840 + t755;
t884 = -m(4) * t824 + t856 * mrSges(4,1) - t855 * mrSges(4,2) - t859 * t905 + t860 * t904 - t888;
t748 = m(3) * t829 + qJDD(2) * mrSges(3,1) - t882 * mrSges(3,2) + t884;
t913 = t748 * t881;
t854 = (-mrSges(4,1) * t880 + mrSges(4,2) * t876) * qJD(2);
t741 = m(4) * t804 + qJDD(3) * mrSges(4,1) - t855 * mrSges(4,3) + qJD(3) * t860 - t854 * t905 + t743;
t897 = t879 * t751 - t875 * t763;
t742 = m(4) * t805 - qJDD(3) * mrSges(4,2) + t856 * mrSges(4,3) - qJD(3) * t859 + t854 * t904 + t897;
t898 = -t876 * t741 + t880 * t742;
t733 = m(3) * t830 - t882 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t898;
t736 = t880 * t741 + t876 * t742;
t735 = m(3) * t841 + t736;
t723 = t733 * t911 - t871 * t735 + t873 * t913;
t721 = m(2) * t857 + t723;
t728 = t881 * t733 - t877 * t748;
t727 = m(2) * t858 + t728;
t910 = t872 * t721 + t870 * t727;
t722 = t733 * t912 + t873 * t735 + t871 * t913;
t899 = -t870 * t721 + t872 * t727;
t895 = m(2) * t869 + t722;
t826 = Ifges(5,5) * t848 - Ifges(5,6) * t847 + Ifges(5,3) * t867;
t729 = mrSges(5,2) * t795 - mrSges(5,3) * t782 + Ifges(5,1) * t818 + Ifges(5,4) * t817 + Ifges(5,5) * t866 - pkin(10) * t755 - t874 * t745 + t878 * t753 - t847 * t826 - t867 * t827;
t737 = -mrSges(5,1) * t795 + mrSges(5,3) * t783 + Ifges(5,4) * t818 + Ifges(5,2) * t817 + Ifges(5,6) * t866 - pkin(4) * t755 - t848 * t826 + t867 * t828 - t919;
t844 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t876 + Ifges(4,6) * t880) * qJD(2);
t719 = -mrSges(4,1) * t824 + mrSges(4,3) * t805 + Ifges(4,4) * t855 + Ifges(4,2) * t856 + Ifges(4,6) * qJDD(3) - pkin(3) * t888 + pkin(9) * t897 + qJD(3) * t846 + t875 * t729 + t879 * t737 - t844 * t905;
t724 = mrSges(4,2) * t824 - mrSges(4,3) * t804 + Ifges(4,1) * t855 + Ifges(4,4) * t856 + Ifges(4,5) * qJDD(3) - pkin(9) * t743 - qJD(3) * t845 + t879 * t729 - t875 * t737 + t844 * t904;
t717 = mrSges(3,2) * t841 - mrSges(3,3) * t829 + Ifges(3,5) * qJDD(2) - t882 * Ifges(3,6) - pkin(8) * t736 - t876 * t719 + t880 * t724;
t718 = -mrSges(3,1) * t841 + mrSges(3,3) * t830 + t882 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t736 - t920;
t890 = pkin(7) * t728 + t717 * t877 + t718 * t881;
t716 = mrSges(3,1) * t829 - mrSges(3,2) * t830 + Ifges(3,3) * qJDD(2) + pkin(2) * t884 + pkin(8) * t898 + t880 * t719 + t876 * t724;
t715 = mrSges(2,2) * t869 - mrSges(2,3) * t857 + t881 * t717 - t877 * t718 + (-t722 * t871 - t723 * t873) * pkin(7);
t714 = -mrSges(2,1) * t869 + mrSges(2,3) * t858 - pkin(1) * t722 - t871 * t716 + t890 * t873;
t1 = [-m(1) * g(1) + t899; -m(1) * g(2) + t910; -m(1) * g(3) + t895; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t910 - t870 * t714 + t872 * t715; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t899 + t872 * t714 + t870 * t715; -mrSges(1,1) * g(2) + mrSges(2,1) * t857 + mrSges(1,2) * g(1) - mrSges(2,2) * t858 + pkin(1) * t723 + t873 * t716 + t890 * t871; t895; t716; t920; -t886; t919; t767;];
tauJB  = t1;
