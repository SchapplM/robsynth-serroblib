% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:56:07
% EndTime: 2019-05-06 01:56:18
% DurationCPUTime: 6.47s
% Computational Cost: add. (73858->315), mult. (142658->369), div. (0->0), fcn. (90997->8), ass. (0->129)
t890 = Ifges(6,1) + Ifges(7,1);
t877 = Ifges(6,4) - Ifges(7,5);
t888 = Ifges(7,4) + Ifges(6,5);
t889 = Ifges(6,2) + Ifges(7,3);
t887 = Ifges(6,6) - Ifges(7,6);
t886 = -Ifges(6,3) - Ifges(7,2);
t844 = sin(qJ(4));
t847 = cos(qJ(4));
t848 = cos(qJ(3));
t872 = qJD(1) * t848;
t817 = t847 * qJD(3) - t844 * t872;
t818 = t844 * qJD(3) + t847 * t872;
t843 = sin(qJ(5));
t881 = cos(qJ(5));
t782 = -t881 * t817 + t843 * t818;
t783 = t843 * t817 + t881 * t818;
t845 = sin(qJ(3));
t833 = t845 * qJD(1);
t830 = t833 + qJD(4);
t829 = qJD(5) + t830;
t885 = t782 * t889 - t783 * t877 - t829 * t887;
t884 = -t877 * t782 + t783 * t890 + t888 * t829;
t846 = sin(qJ(1));
t849 = cos(qJ(1));
t827 = -t849 * g(1) - t846 * g(2);
t883 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t827;
t882 = -pkin(1) - pkin(7);
t880 = mrSges(2,1) - mrSges(3,2);
t879 = -mrSges(6,3) - mrSges(7,2);
t878 = -Ifges(3,4) + Ifges(2,5);
t876 = Ifges(3,5) - Ifges(2,6);
t826 = t846 * g(1) - t849 * g(2);
t851 = qJD(1) ^ 2;
t861 = -t851 * qJ(2) + qJDD(2) - t826;
t794 = t882 * qJDD(1) + t861;
t787 = -t848 * g(3) + t845 * t794;
t819 = (mrSges(4,1) * t845 + mrSges(4,2) * t848) * qJD(1);
t871 = qJD(1) * qJD(3);
t831 = t848 * t871;
t821 = -t845 * qJDD(1) - t831;
t825 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t872;
t793 = t882 * t851 - t883;
t868 = t845 * t871;
t822 = t848 * qJDD(1) - t868;
t763 = (-t822 + t868) * pkin(8) + (-t821 + t831) * pkin(3) + t793;
t820 = (pkin(3) * t845 - pkin(8) * t848) * qJD(1);
t850 = qJD(3) ^ 2;
t768 = -t850 * pkin(3) + qJDD(3) * pkin(8) - t820 * t833 + t787;
t736 = t847 * t763 - t844 * t768;
t781 = t817 * qJD(4) + t844 * qJDD(3) + t847 * t822;
t816 = qJDD(4) - t821;
t732 = (t817 * t830 - t781) * pkin(9) + (t817 * t818 + t816) * pkin(4) + t736;
t737 = t844 * t763 + t847 * t768;
t780 = -t818 * qJD(4) + t847 * qJDD(3) - t844 * t822;
t792 = t830 * pkin(4) - t818 * pkin(9);
t815 = t817 ^ 2;
t734 = -t815 * pkin(4) + t780 * pkin(9) - t830 * t792 + t737;
t730 = t843 * t732 + t881 * t734;
t747 = t783 * qJD(5) - t881 * t780 + t843 * t781;
t771 = t829 * mrSges(6,1) - t783 * mrSges(6,3);
t809 = qJDD(5) + t816;
t758 = t782 * pkin(5) - t783 * qJ(6);
t828 = t829 ^ 2;
t724 = -t828 * pkin(5) + t809 * qJ(6) + 0.2e1 * qJD(6) * t829 - t782 * t758 + t730;
t772 = -t829 * mrSges(7,1) + t783 * mrSges(7,2);
t869 = m(7) * t724 + t809 * mrSges(7,3) + t829 * t772;
t759 = t782 * mrSges(7,1) - t783 * mrSges(7,3);
t873 = -t782 * mrSges(6,1) - t783 * mrSges(6,2) - t759;
t714 = m(6) * t730 - t809 * mrSges(6,2) + t879 * t747 - t829 * t771 + t873 * t782 + t869;
t729 = t881 * t732 - t843 * t734;
t748 = -t782 * qJD(5) + t843 * t780 + t881 * t781;
t770 = -t829 * mrSges(6,2) - t782 * mrSges(6,3);
t725 = -t809 * pkin(5) - t828 * qJ(6) + t783 * t758 + qJDD(6) - t729;
t769 = -t782 * mrSges(7,2) + t829 * mrSges(7,3);
t863 = -m(7) * t725 + t809 * mrSges(7,1) + t829 * t769;
t716 = m(6) * t729 + t809 * mrSges(6,1) + t879 * t748 + t829 * t770 + t873 * t783 + t863;
t709 = t843 * t714 + t881 * t716;
t785 = -t817 * mrSges(5,1) + t818 * mrSges(5,2);
t788 = -t830 * mrSges(5,2) + t817 * mrSges(5,3);
t705 = m(5) * t736 + t816 * mrSges(5,1) - t781 * mrSges(5,3) - t818 * t785 + t830 * t788 + t709;
t789 = t830 * mrSges(5,1) - t818 * mrSges(5,3);
t864 = t881 * t714 - t843 * t716;
t706 = m(5) * t737 - t816 * mrSges(5,2) + t780 * mrSges(5,3) + t817 * t785 - t830 * t789 + t864;
t865 = -t844 * t705 + t847 * t706;
t699 = m(4) * t787 - qJDD(3) * mrSges(4,2) + t821 * mrSges(4,3) - qJD(3) * t825 - t819 * t833 + t865;
t786 = t845 * g(3) + t848 * t794;
t824 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t833;
t767 = -qJDD(3) * pkin(3) - t850 * pkin(8) + t820 * t872 - t786;
t735 = -t780 * pkin(4) - t815 * pkin(9) + t818 * t792 + t767;
t727 = -0.2e1 * qJD(6) * t783 + (t782 * t829 - t748) * qJ(6) + (t783 * t829 + t747) * pkin(5) + t735;
t717 = m(7) * t727 + t747 * mrSges(7,1) - t748 * mrSges(7,3) + t782 * t769 - t783 * t772;
t856 = m(6) * t735 + t747 * mrSges(6,1) + t748 * mrSges(6,2) + t782 * t770 + t783 * t771 + t717;
t853 = -m(5) * t767 + t780 * mrSges(5,1) - t781 * mrSges(5,2) + t817 * t788 - t818 * t789 - t856;
t710 = m(4) * t786 + qJDD(3) * mrSges(4,1) - t822 * mrSges(4,3) + qJD(3) * t824 - t819 * t872 + t853;
t693 = t845 * t699 + t848 * t710;
t799 = -qJDD(1) * pkin(1) + t861;
t860 = -m(3) * t799 + t851 * mrSges(3,3) - t693;
t689 = m(2) * t826 - t851 * mrSges(2,2) + t880 * qJDD(1) + t860;
t797 = t851 * pkin(1) + t883;
t701 = t847 * t705 + t844 * t706;
t859 = -m(4) * t793 + t821 * mrSges(4,1) - t822 * mrSges(4,2) - t824 * t833 - t825 * t872 - t701;
t857 = -m(3) * t797 + t851 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t859;
t696 = m(2) * t827 - t851 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t857;
t875 = t849 * t689 + t846 * t696;
t874 = t887 * t782 - t888 * t783 + t886 * t829;
t867 = -t846 * t689 + t849 * t696;
t866 = t848 * t699 - t845 * t710;
t707 = -mrSges(6,1) * t735 - mrSges(7,1) * t727 + mrSges(7,2) * t724 + mrSges(6,3) * t730 - pkin(5) * t717 - t747 * t889 + t877 * t748 + t874 * t783 + t887 * t809 + t884 * t829;
t708 = mrSges(6,2) * t735 + mrSges(7,2) * t725 - mrSges(6,3) * t729 - mrSges(7,3) * t727 - qJ(6) * t717 - t877 * t747 + t748 * t890 + t874 * t782 + t888 * t809 + t885 * t829;
t774 = Ifges(5,5) * t818 + Ifges(5,6) * t817 + Ifges(5,3) * t830;
t776 = Ifges(5,1) * t818 + Ifges(5,4) * t817 + Ifges(5,5) * t830;
t685 = -mrSges(5,1) * t767 + mrSges(5,3) * t737 + Ifges(5,4) * t781 + Ifges(5,2) * t780 + Ifges(5,6) * t816 - pkin(4) * t856 + pkin(9) * t864 + t881 * t707 + t843 * t708 - t818 * t774 + t830 * t776;
t775 = Ifges(5,4) * t818 + Ifges(5,2) * t817 + Ifges(5,6) * t830;
t687 = mrSges(5,2) * t767 - mrSges(5,3) * t736 + Ifges(5,1) * t781 + Ifges(5,4) * t780 + Ifges(5,5) * t816 - pkin(9) * t709 - t843 * t707 + t881 * t708 + t817 * t774 - t830 * t775;
t807 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t848 - Ifges(4,2) * t845) * qJD(1);
t808 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t848 - Ifges(4,4) * t845) * qJD(1);
t858 = mrSges(4,1) * t786 - mrSges(4,2) * t787 + Ifges(4,5) * t822 + Ifges(4,6) * t821 + Ifges(4,3) * qJDD(3) + pkin(3) * t853 + pkin(8) * t865 + t847 * t685 + t844 * t687 + t807 * t872 + t808 * t833;
t806 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t848 - Ifges(4,6) * t845) * qJD(1);
t682 = mrSges(4,2) * t793 - mrSges(4,3) * t786 + Ifges(4,1) * t822 + Ifges(4,4) * t821 + Ifges(4,5) * qJDD(3) - pkin(8) * t701 - qJD(3) * t807 - t844 * t685 + t847 * t687 - t806 * t833;
t721 = t748 * mrSges(7,2) + t783 * t759 - t863;
t854 = -mrSges(6,1) * t729 + mrSges(7,1) * t725 + mrSges(6,2) * t730 - mrSges(7,3) * t724 + pkin(5) * t721 - qJ(6) * t869 + t886 * t809 + t885 * t783 + (qJ(6) * t759 - t884) * t782 - t888 * t748 + (mrSges(7,2) * qJ(6) + t887) * t747;
t852 = mrSges(5,1) * t736 - mrSges(5,2) * t737 + Ifges(5,5) * t781 + Ifges(5,6) * t780 + Ifges(5,3) * t816 + pkin(4) * t709 + t818 * t775 - t817 * t776 - t854;
t683 = -mrSges(4,1) * t793 + mrSges(4,3) * t787 + Ifges(4,4) * t822 + Ifges(4,2) * t821 + Ifges(4,6) * qJDD(3) - pkin(3) * t701 + qJD(3) * t808 - t806 * t872 - t852;
t691 = qJDD(1) * mrSges(3,2) - t860;
t855 = mrSges(2,1) * t826 - mrSges(2,2) * t827 + mrSges(3,2) * t799 - mrSges(3,3) * t797 - pkin(1) * t691 - pkin(7) * t693 + qJ(2) * t857 + t848 * t682 - t845 * t683 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t692 = -m(3) * g(3) + t866;
t680 = -mrSges(2,3) * t826 + mrSges(3,1) * t799 + pkin(2) * t693 - qJ(2) * t692 + t878 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t876 * t851 + t858;
t679 = -mrSges(3,1) * t797 + mrSges(2,3) * t827 - pkin(1) * t692 - pkin(2) * t859 - pkin(7) * t866 + t880 * g(3) - t876 * qJDD(1) - t845 * t682 - t848 * t683 + t878 * t851;
t1 = [-m(1) * g(1) + t867; -m(1) * g(2) + t875; (-m(1) - m(2) - m(3)) * g(3) + t866; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t875 - t846 * t679 + t849 * t680; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t867 + t849 * t679 + t846 * t680; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t855; t855; t691; t858; t852; -t854; t721;];
tauJB  = t1;
