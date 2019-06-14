% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:58:01
% EndTime: 2019-05-06 02:58:20
% DurationCPUTime: 18.66s
% Computational Cost: add. (324149->345), mult. (639061->428), div. (0->0), fcn. (440157->12), ass. (0->141)
t846 = sin(qJ(1));
t851 = cos(qJ(1));
t829 = t846 * g(1) - g(2) * t851;
t819 = qJDD(1) * pkin(1) + t829;
t830 = -g(1) * t851 - g(2) * t846;
t853 = qJD(1) ^ 2;
t821 = -pkin(1) * t853 + t830;
t840 = sin(pkin(11));
t841 = cos(pkin(11));
t796 = t841 * t819 - t840 * t821;
t782 = -qJDD(1) * pkin(2) - t853 * pkin(7) - t796;
t845 = sin(qJ(3));
t850 = cos(qJ(3));
t869 = qJD(1) * qJD(3);
t868 = t850 * t869;
t823 = qJDD(1) * t845 + t868;
t834 = t845 * t869;
t824 = qJDD(1) * t850 - t834;
t767 = (-t823 - t868) * pkin(8) + (-t824 + t834) * pkin(3) + t782;
t797 = t840 * t819 + t841 * t821;
t783 = -pkin(2) * t853 + qJDD(1) * pkin(7) + t797;
t839 = -g(3) + qJDD(2);
t776 = t850 * t783 + t845 * t839;
t822 = (-pkin(3) * t850 - pkin(8) * t845) * qJD(1);
t852 = qJD(3) ^ 2;
t870 = qJD(1) * t850;
t773 = -pkin(3) * t852 + qJDD(3) * pkin(8) + t822 * t870 + t776;
t844 = sin(qJ(4));
t849 = cos(qJ(4));
t751 = t849 * t767 - t844 * t773;
t871 = qJD(1) * t845;
t817 = qJD(3) * t849 - t844 * t871;
t791 = qJD(4) * t817 + qJDD(3) * t844 + t823 * t849;
t816 = qJDD(4) - t824;
t818 = qJD(3) * t844 + t849 * t871;
t832 = qJD(4) - t870;
t741 = (t817 * t832 - t791) * pkin(9) + (t817 * t818 + t816) * pkin(4) + t751;
t752 = t844 * t767 + t849 * t773;
t790 = -qJD(4) * t818 + qJDD(3) * t849 - t823 * t844;
t801 = pkin(4) * t832 - pkin(9) * t818;
t815 = t817 ^ 2;
t743 = -pkin(4) * t815 + pkin(9) * t790 - t801 * t832 + t752;
t843 = sin(qJ(5));
t848 = cos(qJ(5));
t729 = t848 * t741 - t843 * t743;
t794 = t817 * t848 - t818 * t843;
t759 = qJD(5) * t794 + t790 * t843 + t791 * t848;
t795 = t817 * t843 + t818 * t848;
t812 = qJDD(5) + t816;
t831 = qJD(5) + t832;
t726 = (t794 * t831 - t759) * pkin(10) + (t794 * t795 + t812) * pkin(5) + t729;
t730 = t843 * t741 + t848 * t743;
t758 = -qJD(5) * t795 + t790 * t848 - t791 * t843;
t779 = pkin(5) * t831 - pkin(10) * t795;
t793 = t794 ^ 2;
t727 = -pkin(5) * t793 + pkin(10) * t758 - t779 * t831 + t730;
t842 = sin(qJ(6));
t847 = cos(qJ(6));
t725 = t726 * t842 + t727 * t847;
t775 = -t845 * t783 + t839 * t850;
t772 = -qJDD(3) * pkin(3) - pkin(8) * t852 + t822 * t871 - t775;
t750 = -pkin(4) * t790 - pkin(9) * t815 + t818 * t801 + t772;
t732 = -pkin(5) * t758 - pkin(10) * t793 + t779 * t795 + t750;
t771 = t794 * t842 + t795 * t847;
t737 = -qJD(6) * t771 + t758 * t847 - t759 * t842;
t770 = t794 * t847 - t795 * t842;
t738 = qJD(6) * t770 + t758 * t842 + t759 * t847;
t826 = qJD(6) + t831;
t745 = Ifges(7,5) * t771 + Ifges(7,6) * t770 + Ifges(7,3) * t826;
t747 = Ifges(7,1) * t771 + Ifges(7,4) * t770 + Ifges(7,5) * t826;
t805 = qJDD(6) + t812;
t713 = -mrSges(7,1) * t732 + mrSges(7,3) * t725 + Ifges(7,4) * t738 + Ifges(7,2) * t737 + Ifges(7,6) * t805 - t745 * t771 + t747 * t826;
t724 = t726 * t847 - t727 * t842;
t746 = Ifges(7,4) * t771 + Ifges(7,2) * t770 + Ifges(7,6) * t826;
t714 = mrSges(7,2) * t732 - mrSges(7,3) * t724 + Ifges(7,1) * t738 + Ifges(7,4) * t737 + Ifges(7,5) * t805 + t745 * t770 - t746 * t826;
t764 = Ifges(6,5) * t795 + Ifges(6,6) * t794 + Ifges(6,3) * t831;
t766 = Ifges(6,1) * t795 + Ifges(6,4) * t794 + Ifges(6,5) * t831;
t760 = -mrSges(7,2) * t826 + mrSges(7,3) * t770;
t761 = mrSges(7,1) * t826 - mrSges(7,3) * t771;
t862 = m(7) * t732 - t737 * mrSges(7,1) + t738 * mrSges(7,2) - t770 * t760 + t771 * t761;
t753 = -mrSges(7,1) * t770 + mrSges(7,2) * t771;
t720 = m(7) * t724 + mrSges(7,1) * t805 - mrSges(7,3) * t738 - t753 * t771 + t760 * t826;
t721 = m(7) * t725 - mrSges(7,2) * t805 + mrSges(7,3) * t737 + t753 * t770 - t761 * t826;
t863 = -t720 * t842 + t847 * t721;
t700 = -mrSges(6,1) * t750 + mrSges(6,3) * t730 + Ifges(6,4) * t759 + Ifges(6,2) * t758 + Ifges(6,6) * t812 - pkin(5) * t862 + pkin(10) * t863 + t847 * t713 + t842 * t714 - t795 * t764 + t831 * t766;
t712 = t847 * t720 + t842 * t721;
t765 = Ifges(6,4) * t795 + Ifges(6,2) * t794 + Ifges(6,6) * t831;
t701 = mrSges(6,2) * t750 - mrSges(6,3) * t729 + Ifges(6,1) * t759 + Ifges(6,4) * t758 + Ifges(6,5) * t812 - pkin(10) * t712 - t713 * t842 + t714 * t847 + t764 * t794 - t765 * t831;
t784 = Ifges(5,5) * t818 + Ifges(5,6) * t817 + Ifges(5,3) * t832;
t786 = Ifges(5,1) * t818 + Ifges(5,4) * t817 + Ifges(5,5) * t832;
t777 = -mrSges(6,2) * t831 + mrSges(6,3) * t794;
t778 = mrSges(6,1) * t831 - mrSges(6,3) * t795;
t858 = m(6) * t750 - t758 * mrSges(6,1) + t759 * mrSges(6,2) - t794 * t777 + t795 * t778 + t862;
t774 = -mrSges(6,1) * t794 + mrSges(6,2) * t795;
t709 = m(6) * t729 + mrSges(6,1) * t812 - mrSges(6,3) * t759 - t774 * t795 + t777 * t831 + t712;
t710 = m(6) * t730 - mrSges(6,2) * t812 + mrSges(6,3) * t758 + t774 * t794 - t778 * t831 + t863;
t864 = -t709 * t843 + t848 * t710;
t684 = -mrSges(5,1) * t772 + mrSges(5,3) * t752 + Ifges(5,4) * t791 + Ifges(5,2) * t790 + Ifges(5,6) * t816 - pkin(4) * t858 + pkin(9) * t864 + t848 * t700 + t843 * t701 - t818 * t784 + t832 * t786;
t705 = t848 * t709 + t843 * t710;
t785 = Ifges(5,4) * t818 + Ifges(5,2) * t817 + Ifges(5,6) * t832;
t685 = mrSges(5,2) * t772 - mrSges(5,3) * t751 + Ifges(5,1) * t791 + Ifges(5,4) * t790 + Ifges(5,5) * t816 - pkin(9) * t705 - t700 * t843 + t701 * t848 + t784 * t817 - t785 * t832;
t798 = -mrSges(5,1) * t817 + mrSges(5,2) * t818;
t799 = -mrSges(5,2) * t832 + mrSges(5,3) * t817;
t703 = m(5) * t751 + mrSges(5,1) * t816 - mrSges(5,3) * t791 - t798 * t818 + t799 * t832 + t705;
t800 = mrSges(5,1) * t832 - mrSges(5,3) * t818;
t704 = m(5) * t752 - mrSges(5,2) * t816 + mrSges(5,3) * t790 + t798 * t817 - t800 * t832 + t864;
t699 = -t703 * t844 + t849 * t704;
t722 = -m(5) * t772 + t790 * mrSges(5,1) - t791 * mrSges(5,2) + t817 * t799 - t818 * t800 - t858;
t810 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t845 + Ifges(4,2) * t850) * qJD(1);
t811 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t845 + Ifges(4,4) * t850) * qJD(1);
t873 = mrSges(4,1) * t775 - mrSges(4,2) * t776 + Ifges(4,5) * t823 + Ifges(4,6) * t824 + Ifges(4,3) * qJDD(3) + pkin(3) * t722 + pkin(8) * t699 + t849 * t684 + t844 * t685 + (t810 * t845 - t811 * t850) * qJD(1);
t820 = (-mrSges(4,1) * t850 + mrSges(4,2) * t845) * qJD(1);
t827 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t871;
t697 = m(4) * t776 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t824 - qJD(3) * t827 + t820 * t870 + t699;
t828 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t870;
t716 = m(4) * t775 + qJDD(3) * mrSges(4,1) - t823 * mrSges(4,3) + qJD(3) * t828 - t820 * t871 + t722;
t865 = t850 * t697 - t716 * t845;
t688 = m(3) * t797 - mrSges(3,1) * t853 - qJDD(1) * mrSges(3,2) + t865;
t698 = t703 * t849 + t704 * t844;
t857 = -m(4) * t782 + t824 * mrSges(4,1) - mrSges(4,2) * t823 - t827 * t871 + t828 * t870 - t698;
t693 = m(3) * t796 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t853 + t857;
t681 = t840 * t688 + t841 * t693;
t678 = m(2) * t829 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t853 + t681;
t866 = t841 * t688 - t693 * t840;
t679 = m(2) * t830 - mrSges(2,1) * t853 - qJDD(1) * mrSges(2,2) + t866;
t872 = t851 * t678 + t846 * t679;
t691 = t845 * t697 + t850 * t716;
t689 = m(3) * t839 + t691;
t867 = -t678 * t846 + t851 * t679;
t860 = -mrSges(7,1) * t724 + mrSges(7,2) * t725 - Ifges(7,5) * t738 - Ifges(7,6) * t737 - Ifges(7,3) * t805 - t771 * t746 + t770 * t747;
t809 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t845 + Ifges(4,6) * t850) * qJD(1);
t674 = mrSges(4,2) * t782 - mrSges(4,3) * t775 + Ifges(4,1) * t823 + Ifges(4,4) * t824 + Ifges(4,5) * qJDD(3) - pkin(8) * t698 - qJD(3) * t810 - t684 * t844 + t685 * t849 + t809 * t870;
t856 = -mrSges(6,1) * t729 + mrSges(6,2) * t730 - Ifges(6,5) * t759 - Ifges(6,6) * t758 - Ifges(6,3) * t812 - pkin(5) * t712 - t795 * t765 + t794 * t766 + t860;
t854 = mrSges(5,1) * t751 - mrSges(5,2) * t752 + Ifges(5,5) * t791 + Ifges(5,6) * t790 + Ifges(5,3) * t816 + pkin(4) * t705 + t818 * t785 - t817 * t786 - t856;
t683 = -mrSges(4,1) * t782 + mrSges(4,3) * t776 + Ifges(4,4) * t823 + Ifges(4,2) * t824 + Ifges(4,6) * qJDD(3) - pkin(3) * t698 + qJD(3) * t811 - t809 * t871 - t854;
t859 = mrSges(2,1) * t829 + mrSges(3,1) * t796 - mrSges(2,2) * t830 - mrSges(3,2) * t797 + pkin(1) * t681 + pkin(2) * t857 + pkin(7) * t865 + t845 * t674 + t850 * t683 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t672 = -mrSges(3,1) * t839 + mrSges(3,3) * t797 + t853 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t691 - t873;
t671 = mrSges(3,2) * t839 - mrSges(3,3) * t796 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t853 - pkin(7) * t691 + t674 * t850 - t683 * t845;
t670 = -mrSges(2,2) * g(3) - mrSges(2,3) * t829 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t853 - qJ(2) * t681 + t671 * t841 - t672 * t840;
t669 = mrSges(2,1) * g(3) + mrSges(2,3) * t830 + t853 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t689 + qJ(2) * t866 + t840 * t671 + t841 * t672;
t1 = [-m(1) * g(1) + t867; -m(1) * g(2) + t872; (-m(1) - m(2)) * g(3) + t689; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t872 - t846 * t669 + t851 * t670; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t867 + t851 * t669 + t846 * t670; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t859; t859; t689; t873; t854; -t856; -t860;];
tauJB  = t1;
