% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:13:20
% EndTime: 2019-05-06 01:13:31
% DurationCPUTime: 8.45s
% Computational Cost: add. (133053->321), mult. (258270->386), div. (0->0), fcn. (170169->10), ass. (0->130)
t875 = Ifges(6,4) + Ifges(7,4);
t884 = Ifges(6,2) + Ifges(7,2);
t880 = Ifges(6,6) + Ifges(7,6);
t845 = sin(qJ(4));
t849 = cos(qJ(4));
t846 = sin(qJ(3));
t871 = qJD(1) * t846;
t821 = qJD(3) * t845 + t849 * t871;
t850 = cos(qJ(3));
t869 = qJD(1) * qJD(3);
t866 = t850 * t869;
t826 = qJDD(1) * t846 + t866;
t793 = -qJD(4) * t821 + qJDD(3) * t849 - t826 * t845;
t820 = qJD(3) * t849 - t845 * t871;
t794 = qJD(4) * t820 + qJDD(3) * t845 + t826 * t849;
t844 = sin(qJ(5));
t848 = cos(qJ(5));
t797 = t820 * t848 - t821 * t844;
t756 = qJD(5) * t797 + t793 * t844 + t794 * t848;
t798 = t820 * t844 + t821 * t848;
t773 = -mrSges(7,1) * t797 + mrSges(7,2) * t798;
t847 = sin(qJ(1));
t851 = cos(qJ(1));
t831 = t847 * g(1) - g(2) * t851;
t822 = qJDD(1) * pkin(1) + t831;
t832 = -g(1) * t851 - g(2) * t847;
t853 = qJD(1) ^ 2;
t824 = -pkin(1) * t853 + t832;
t842 = sin(pkin(10));
t843 = cos(pkin(10));
t799 = t843 * t822 - t842 * t824;
t785 = -qJDD(1) * pkin(2) - t853 * pkin(7) - t799;
t836 = t846 * t869;
t827 = qJDD(1) * t850 - t836;
t766 = (-t826 - t866) * pkin(8) + (-t827 + t836) * pkin(3) + t785;
t800 = t842 * t822 + t843 * t824;
t786 = -pkin(2) * t853 + qJDD(1) * pkin(7) + t800;
t841 = -g(3) + qJDD(2);
t777 = t850 * t786 + t846 * t841;
t825 = (-pkin(3) * t850 - pkin(8) * t846) * qJD(1);
t852 = qJD(3) ^ 2;
t870 = qJD(1) * t850;
t772 = -pkin(3) * t852 + qJDD(3) * pkin(8) + t825 * t870 + t777;
t743 = t849 * t766 - t845 * t772;
t819 = qJDD(4) - t827;
t834 = qJD(4) - t870;
t739 = (t820 * t834 - t794) * pkin(9) + (t820 * t821 + t819) * pkin(4) + t743;
t744 = t845 * t766 + t849 * t772;
t804 = pkin(4) * t834 - pkin(9) * t821;
t818 = t820 ^ 2;
t741 = -pkin(4) * t818 + pkin(9) * t793 - t804 * t834 + t744;
t733 = t848 * t739 - t844 * t741;
t815 = qJDD(5) + t819;
t833 = qJD(5) + t834;
t728 = -0.2e1 * qJD(6) * t798 + (t797 * t833 - t756) * qJ(6) + (t797 * t798 + t815) * pkin(5) + t733;
t778 = -mrSges(7,2) * t833 + mrSges(7,3) * t797;
t868 = m(7) * t728 + t815 * mrSges(7,1) + t833 * t778;
t725 = -t756 * mrSges(7,3) - t798 * t773 + t868;
t734 = t844 * t739 + t848 * t741;
t755 = -qJD(5) * t798 + t793 * t848 - t794 * t844;
t780 = pkin(5) * t833 - qJ(6) * t798;
t796 = t797 ^ 2;
t731 = -pkin(5) * t796 + qJ(6) * t755 + 0.2e1 * qJD(6) * t797 - t780 * t833 + t734;
t881 = Ifges(6,5) + Ifges(7,5);
t882 = Ifges(6,1) + Ifges(7,1);
t872 = -t875 * t797 - t882 * t798 - t881 * t833;
t878 = t884 * t797 + t875 * t798 + t880 * t833;
t879 = Ifges(6,3) + Ifges(7,3);
t883 = mrSges(6,1) * t733 + mrSges(7,1) * t728 - mrSges(6,2) * t734 - mrSges(7,2) * t731 + pkin(5) * t725 + t880 * t755 + t881 * t756 + t872 * t797 + t878 * t798 + t879 * t815;
t774 = -mrSges(6,1) * t797 + mrSges(6,2) * t798;
t779 = -mrSges(6,2) * t833 + mrSges(6,3) * t797;
t717 = m(6) * t733 + t815 * mrSges(6,1) + t833 * t779 + (-t773 - t774) * t798 + (-mrSges(6,3) - mrSges(7,3)) * t756 + t868;
t781 = mrSges(7,1) * t833 - mrSges(7,3) * t798;
t782 = mrSges(6,1) * t833 - mrSges(6,3) * t798;
t867 = m(7) * t731 + t755 * mrSges(7,3) + t797 * t773;
t720 = m(6) * t734 + t755 * mrSges(6,3) + t797 * t774 + (-t781 - t782) * t833 + (-mrSges(6,2) - mrSges(7,2)) * t815 + t867;
t715 = t848 * t717 + t844 * t720;
t788 = Ifges(5,4) * t821 + Ifges(5,2) * t820 + Ifges(5,6) * t834;
t789 = Ifges(5,1) * t821 + Ifges(5,4) * t820 + Ifges(5,5) * t834;
t877 = mrSges(5,1) * t743 - mrSges(5,2) * t744 + Ifges(5,5) * t794 + Ifges(5,6) * t793 + Ifges(5,3) * t819 + pkin(4) * t715 + t821 * t788 - t820 * t789 + t883;
t776 = -t846 * t786 + t841 * t850;
t771 = -qJDD(3) * pkin(3) - pkin(8) * t852 + t825 * t871 - t776;
t742 = -pkin(4) * t793 - pkin(9) * t818 + t821 * t804 + t771;
t736 = -pkin(5) * t755 - qJ(6) * t796 + t780 * t798 + qJDD(6) + t742;
t729 = m(7) * t736 - t755 * mrSges(7,1) + t756 * mrSges(7,2) - t797 * t778 + t798 * t781;
t873 = -t880 * t797 - t881 * t798 - t879 * t833;
t710 = -mrSges(6,1) * t742 + mrSges(6,3) * t734 - mrSges(7,1) * t736 + mrSges(7,3) * t731 - pkin(5) * t729 + qJ(6) * t867 + (-qJ(6) * t781 - t872) * t833 + (-mrSges(7,2) * qJ(6) + t880) * t815 + t873 * t798 + t875 * t756 + t884 * t755;
t714 = mrSges(6,2) * t742 + mrSges(7,2) * t736 - mrSges(6,3) * t733 - mrSges(7,3) * t728 - qJ(6) * t725 + t875 * t755 + t882 * t756 - t873 * t797 + t881 * t815 - t878 * t833;
t787 = Ifges(5,5) * t821 + Ifges(5,6) * t820 + Ifges(5,3) * t834;
t858 = m(6) * t742 - t755 * mrSges(6,1) + t756 * mrSges(6,2) - t797 * t779 + t798 * t782 + t729;
t861 = -t717 * t844 + t848 * t720;
t694 = -mrSges(5,1) * t771 + mrSges(5,3) * t744 + Ifges(5,4) * t794 + Ifges(5,2) * t793 + Ifges(5,6) * t819 - pkin(4) * t858 + pkin(9) * t861 + t848 * t710 + t844 * t714 - t821 * t787 + t834 * t789;
t695 = mrSges(5,2) * t771 - mrSges(5,3) * t743 + Ifges(5,1) * t794 + Ifges(5,4) * t793 + Ifges(5,5) * t819 - pkin(9) * t715 - t710 * t844 + t714 * t848 + t787 * t820 - t788 * t834;
t801 = -mrSges(5,1) * t820 + mrSges(5,2) * t821;
t802 = -mrSges(5,2) * t834 + mrSges(5,3) * t820;
t712 = m(5) * t743 + mrSges(5,1) * t819 - mrSges(5,3) * t794 - t801 * t821 + t802 * t834 + t715;
t803 = mrSges(5,1) * t834 - mrSges(5,3) * t821;
t713 = m(5) * t744 - mrSges(5,2) * t819 + mrSges(5,3) * t793 + t801 * t820 - t803 * t834 + t861;
t709 = -t712 * t845 + t849 * t713;
t723 = -m(5) * t771 + t793 * mrSges(5,1) - t794 * mrSges(5,2) + t820 * t802 - t821 * t803 - t858;
t813 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t846 + Ifges(4,2) * t850) * qJD(1);
t814 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t846 + Ifges(4,4) * t850) * qJD(1);
t876 = mrSges(4,1) * t776 - mrSges(4,2) * t777 + Ifges(4,5) * t826 + Ifges(4,6) * t827 + Ifges(4,3) * qJDD(3) + pkin(3) * t723 + pkin(8) * t709 + t849 * t694 + t845 * t695 + (t813 * t846 - t814 * t850) * qJD(1);
t823 = (-mrSges(4,1) * t850 + mrSges(4,2) * t846) * qJD(1);
t829 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t871;
t707 = m(4) * t777 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t827 - qJD(3) * t829 + t823 * t870 + t709;
t830 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t870;
t722 = m(4) * t776 + qJDD(3) * mrSges(4,1) - t826 * mrSges(4,3) + qJD(3) * t830 - t823 * t871 + t723;
t862 = t850 * t707 - t722 * t846;
t698 = m(3) * t800 - mrSges(3,1) * t853 - qJDD(1) * mrSges(3,2) + t862;
t708 = t712 * t849 + t713 * t845;
t857 = -m(4) * t785 + t827 * mrSges(4,1) - mrSges(4,2) * t826 - t829 * t871 + t830 * t870 - t708;
t703 = m(3) * t799 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t853 + t857;
t691 = t842 * t698 + t843 * t703;
t688 = m(2) * t831 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t853 + t691;
t863 = t843 * t698 - t703 * t842;
t689 = m(2) * t832 - mrSges(2,1) * t853 - qJDD(1) * mrSges(2,2) + t863;
t874 = t851 * t688 + t847 * t689;
t701 = t846 * t707 + t850 * t722;
t699 = m(3) * t841 + t701;
t864 = -t688 * t847 + t851 * t689;
t812 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t846 + Ifges(4,6) * t850) * qJD(1);
t684 = mrSges(4,2) * t785 - mrSges(4,3) * t776 + Ifges(4,1) * t826 + Ifges(4,4) * t827 + Ifges(4,5) * qJDD(3) - pkin(8) * t708 - qJD(3) * t813 - t694 * t845 + t695 * t849 + t812 * t870;
t693 = -mrSges(4,1) * t785 + mrSges(4,3) * t777 + Ifges(4,4) * t826 + Ifges(4,2) * t827 + Ifges(4,6) * qJDD(3) - pkin(3) * t708 + qJD(3) * t814 - t812 * t871 - t877;
t859 = mrSges(2,1) * t831 + mrSges(3,1) * t799 - mrSges(2,2) * t832 - mrSges(3,2) * t800 + pkin(1) * t691 + pkin(2) * t857 + pkin(7) * t862 + t846 * t684 + t850 * t693 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t682 = -mrSges(3,1) * t841 + mrSges(3,3) * t800 + t853 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t701 - t876;
t681 = mrSges(3,2) * t841 - mrSges(3,3) * t799 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t853 - pkin(7) * t701 + t684 * t850 - t693 * t846;
t680 = -mrSges(2,2) * g(3) - mrSges(2,3) * t831 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t853 - qJ(2) * t691 + t681 * t843 - t682 * t842;
t679 = mrSges(2,1) * g(3) + mrSges(2,3) * t832 + t853 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t699 + qJ(2) * t863 + t842 * t681 + t843 * t682;
t1 = [-m(1) * g(1) + t864; -m(1) * g(2) + t874; (-m(1) - m(2)) * g(3) + t699; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t874 - t847 * t679 + t851 * t680; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t864 + t851 * t679 + t847 * t680; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t859; t859; t699; t876; t877; t883; t729;];
tauJB  = t1;
