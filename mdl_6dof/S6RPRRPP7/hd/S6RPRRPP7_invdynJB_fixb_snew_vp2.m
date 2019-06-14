% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:46:36
% EndTime: 2019-05-05 21:46:40
% DurationCPUTime: 2.75s
% Computational Cost: add. (25899->292), mult. (48912->328), div. (0->0), fcn. (27453->6), ass. (0->116)
t885 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t862 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t861 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t884 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t860 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t883 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t832 = sin(qJ(1));
t834 = cos(qJ(1));
t815 = -t834 * g(1) - t832 * g(2);
t882 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t815;
t836 = qJD(1) ^ 2;
t876 = (-pkin(1) - pkin(7));
t783 = (t876 * t836) - t882;
t831 = sin(qJ(3));
t833 = cos(qJ(3));
t864 = qJD(1) * qJD(3);
t853 = t833 * t864;
t808 = -t831 * qJDD(1) - t853;
t854 = t831 * t864;
t809 = t833 * qJDD(1) - t854;
t734 = (-t809 + t854) * pkin(8) + (-t808 + t853) * pkin(3) + t783;
t814 = t832 * g(1) - t834 * g(2);
t845 = -t836 * qJ(2) + qJDD(2) - t814;
t784 = t876 * qJDD(1) + t845;
t772 = -t833 * g(3) + t831 * t784;
t807 = (pkin(3) * t831 - pkin(8) * t833) * qJD(1);
t835 = qJD(3) ^ 2;
t865 = t831 * qJD(1);
t738 = -t835 * pkin(3) + qJDD(3) * pkin(8) - t807 * t865 + t772;
t830 = sin(qJ(4));
t875 = cos(qJ(4));
t731 = t875 * t734 - t830 * t738;
t866 = qJD(1) * t833;
t804 = -t875 * qJD(3) + t830 * t866;
t805 = t830 * qJD(3) + t875 * t866;
t767 = t804 * pkin(4) - t805 * qJ(5);
t803 = qJDD(4) - t808;
t817 = qJD(4) + t865;
t816 = t817 ^ 2;
t729 = -t803 * pkin(4) - t816 * qJ(5) + t805 * t767 + qJDD(5) - t731;
t779 = -t804 * mrSges(6,2) + t817 * mrSges(6,3);
t881 = -m(6) * t729 + t803 * mrSges(6,1) + t817 * t779;
t763 = -t804 * qJD(4) + t830 * qJDD(3) + t875 * t809;
t771 = t831 * g(3) + t833 * t784;
t844 = qJDD(3) * pkin(3) + t835 * pkin(8) - t807 * t866 + t771;
t869 = t804 * t817;
t880 = (-t763 + t869) * qJ(5) - t844;
t773 = t817 * mrSges(7,2) + t804 * mrSges(7,3);
t878 = -0.2e1 * t805;
t722 = qJD(6) * t878 + (-t763 - t869) * qJ(6) + (t804 * t805 - t803) * pkin(5) + t729;
t769 = -t804 * mrSges(7,1) + t805 * mrSges(7,2);
t848 = -m(7) * t722 + t763 * mrSges(7,3) + t805 * t769;
t720 = -t803 * mrSges(7,1) - t817 * t773 - t848;
t768 = t804 * mrSges(6,1) - t805 * mrSges(6,3);
t718 = t763 * mrSges(6,2) + t805 * t768 + t720 - t881;
t732 = t830 * t734 + t875 * t738;
t877 = 2 * qJD(5);
t728 = -t816 * pkin(4) + t803 * qJ(5) - t804 * t767 + t817 * t877 + t732;
t762 = t805 * qJD(4) - t875 * qJDD(3) + t830 * t809;
t775 = -t817 * pkin(5) - t805 * qJ(6);
t802 = t804 ^ 2;
t724 = -t802 * pkin(5) + t762 * qJ(6) + 0.2e1 * qJD(6) * t804 + t817 * t775 + t728;
t776 = -t817 * mrSges(7,1) - t805 * mrSges(7,3);
t778 = -t817 * mrSges(6,1) + t805 * mrSges(6,2);
t858 = m(7) * t724 + t762 * mrSges(7,3) + t804 * t769;
t846 = m(6) * t728 + t803 * mrSges(6,3) + t817 * t778 + t858;
t855 = -t862 * t804 + t885 * t805 + t861 * t817;
t856 = t884 * t804 + t862 * t805 + t860 * t817;
t879 = -t860 * t762 + t861 * t763 + t883 * t803 + t855 * t804 + t856 * t805 + mrSges(5,1) * t731 - mrSges(6,1) * t729 - mrSges(7,1) * t722 - mrSges(5,2) * t732 + mrSges(7,2) * t724 + mrSges(6,3) * t728 - pkin(4) * t718 - pkin(5) * t720 + qJ(5) * (-t762 * mrSges(6,2) + t803 * mrSges(7,2) - t804 * t768 + t817 * t776 + t846);
t873 = mrSges(2,1) - mrSges(3,2);
t872 = -mrSges(5,3) - mrSges(6,2);
t871 = -Ifges(3,4) + Ifges(2,5);
t870 = (Ifges(3,5) - Ifges(2,6));
t806 = (mrSges(4,1) * t831 + mrSges(4,2) * t833) * qJD(1);
t812 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t866;
t774 = -t817 * mrSges(5,2) - t804 * mrSges(5,3);
t867 = -t804 * mrSges(5,1) - t805 * mrSges(5,2) - t768;
t714 = m(5) * t731 + (t773 + t774) * t817 + t867 * t805 + (mrSges(5,1) + mrSges(7,1)) * t803 + t872 * t763 + t848 + t881;
t777 = t817 * mrSges(5,1) - t805 * mrSges(5,3);
t715 = m(5) * t732 + (t776 - t777) * t817 + t867 * t804 + (-mrSges(5,2) + mrSges(7,2)) * t803 + t872 * t762 + t846;
t850 = -t830 * t714 + t875 * t715;
t707 = m(4) * t772 - qJDD(3) * mrSges(4,2) + t808 * mrSges(4,3) - qJD(3) * t812 - t806 * t865 + t850;
t811 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t865;
t726 = -t802 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t762 + (-pkin(4) * t817 + t775 + t877) * t805 - t880;
t721 = m(7) * t726 - t762 * mrSges(7,1) + t763 * mrSges(7,2) - t804 * t773 + t805 * t776;
t730 = qJD(5) * t878 + (t805 * t817 + t762) * pkin(4) + t880;
t717 = m(6) * t730 + t762 * mrSges(6,1) - t763 * mrSges(6,3) - t805 * t778 + t804 * t779 - t721;
t837 = m(5) * t844 - t762 * mrSges(5,1) - t763 * mrSges(5,2) - t804 * t774 - t805 * t777 - t717;
t710 = m(4) * t771 + qJDD(3) * mrSges(4,1) - t809 * mrSges(4,3) + qJD(3) * t811 - t806 * t866 + t837;
t699 = t831 * t707 + t833 * t710;
t789 = -qJDD(1) * pkin(1) + t845;
t843 = -m(3) * t789 + (t836 * mrSges(3,3)) - t699;
t695 = m(2) * t814 - (t836 * mrSges(2,2)) + t873 * qJDD(1) + t843;
t787 = t836 * pkin(1) + t882;
t709 = t875 * t714 + t830 * t715;
t842 = -m(4) * t783 + t808 * mrSges(4,1) - t809 * mrSges(4,2) - t811 * t865 - t812 * t866 - t709;
t840 = -m(3) * t787 + (t836 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t842;
t704 = m(2) * t815 - (t836 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t840;
t868 = t834 * t695 + t832 * t704;
t857 = t860 * t804 - t861 * t805 - t883 * t817;
t852 = -t832 * t695 + t834 * t704;
t851 = t833 * t707 - t831 * t710;
t693 = mrSges(5,1) * t844 + mrSges(5,3) * t732 - mrSges(6,1) * t730 + mrSges(6,2) * t728 + mrSges(7,1) * t726 - mrSges(7,3) * t724 + pkin(5) * t721 - qJ(6) * t858 - pkin(4) * t717 + (-qJ(6) * t776 + t855) * t817 + t857 * t805 + (-qJ(6) * mrSges(7,2) + t860) * t803 + t862 * t763 + t884 * t762;
t701 = -mrSges(5,2) * t844 + mrSges(6,2) * t729 + mrSges(7,2) * t726 - mrSges(5,3) * t731 - mrSges(6,3) * t730 - mrSges(7,3) * t722 - qJ(5) * t717 - qJ(6) * t720 - t862 * t762 + t885 * t763 + t861 * t803 + t857 * t804 - t856 * t817;
t792 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t833 - Ifges(4,2) * t831) * qJD(1);
t793 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t833 - Ifges(4,4) * t831) * qJD(1);
t841 = mrSges(4,1) * t771 - mrSges(4,2) * t772 + Ifges(4,5) * t809 + Ifges(4,6) * t808 + Ifges(4,3) * qJDD(3) + pkin(3) * t837 + pkin(8) * t850 + t875 * t693 + t830 * t701 + t792 * t866 + t793 * t865;
t791 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t833 - Ifges(4,6) * t831) * qJD(1);
t690 = mrSges(4,2) * t783 - mrSges(4,3) * t771 + Ifges(4,1) * t809 + Ifges(4,4) * t808 + Ifges(4,5) * qJDD(3) - pkin(8) * t709 - qJD(3) * t792 - t830 * t693 + t875 * t701 - t791 * t865;
t691 = -mrSges(4,1) * t783 + mrSges(4,3) * t772 + Ifges(4,4) * t809 + Ifges(4,2) * t808 + Ifges(4,6) * qJDD(3) - pkin(3) * t709 + qJD(3) * t793 - t791 * t866 - t879;
t697 = qJDD(1) * mrSges(3,2) - t843;
t839 = mrSges(2,1) * t814 - mrSges(2,2) * t815 + mrSges(3,2) * t789 - mrSges(3,3) * t787 - pkin(1) * t697 - pkin(7) * t699 + qJ(2) * t840 + t833 * t690 - t831 * t691 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t698 = -m(3) * g(3) + t851;
t688 = (t870 * t836) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t814 + mrSges(3,1) * t789 + pkin(2) * t699 - qJ(2) * t698 + t871 * qJDD(1) + t841;
t687 = -mrSges(3,1) * t787 + mrSges(2,3) * t815 - pkin(1) * t698 - pkin(2) * t842 - pkin(7) * t851 + t873 * g(3) - t870 * qJDD(1) - t831 * t690 - t833 * t691 + t871 * t836;
t1 = [-m(1) * g(1) + t852; -m(1) * g(2) + t868; (-m(1) - m(2) - m(3)) * g(3) + t851; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t868 - t832 * t687 + t834 * t688; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t852 + t834 * t687 + t832 * t688; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t839; t839; t697; t841; t879; t718; t721;];
tauJB  = t1;
