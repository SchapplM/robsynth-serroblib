% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:45:19
% EndTime: 2019-12-05 18:45:26
% DurationCPUTime: 6.65s
% Computational Cost: add. (71308->292), mult. (157129->357), div. (0->0), fcn. (108667->8), ass. (0->117)
t810 = Ifges(5,4) + Ifges(6,4);
t821 = Ifges(5,2) + Ifges(6,2);
t816 = Ifges(5,6) + Ifges(6,6);
t779 = sin(qJ(3));
t780 = sin(qJ(2));
t783 = cos(qJ(3));
t784 = cos(qJ(2));
t752 = (t779 * t784 + t780 * t783) * qJD(1);
t803 = qJD(1) * qJD(2);
t759 = t780 * qJDD(1) + t784 * t803;
t760 = t784 * qJDD(1) - t780 * t803;
t724 = -t752 * qJD(3) - t779 * t759 + t783 * t760;
t751 = (-t779 * t780 + t783 * t784) * qJD(1);
t725 = t751 * qJD(3) + t783 * t759 + t779 * t760;
t778 = sin(qJ(4));
t782 = cos(qJ(4));
t738 = t782 * t751 - t778 * t752;
t695 = t738 * qJD(4) + t778 * t724 + t782 * t725;
t739 = t778 * t751 + t782 * t752;
t712 = -t738 * mrSges(6,1) + t739 * mrSges(6,2);
t781 = sin(qJ(1));
t785 = cos(qJ(1));
t766 = -t785 * g(1) - t781 * g(2);
t786 = qJD(1) ^ 2;
t754 = -t786 * pkin(1) + qJDD(1) * pkin(6) + t766;
t809 = t780 * t754;
t811 = pkin(2) * t786;
t718 = qJDD(2) * pkin(2) - t759 * pkin(7) - t809 + (pkin(7) * t803 + t780 * t811 - g(3)) * t784;
t742 = -t780 * g(3) + t784 * t754;
t805 = qJD(1) * t780;
t764 = qJD(2) * pkin(2) - pkin(7) * t805;
t777 = t784 ^ 2;
t719 = t760 * pkin(7) - qJD(2) * t764 - t777 * t811 + t742;
t700 = t783 * t718 - t779 * t719;
t774 = qJDD(2) + qJDD(3);
t775 = qJD(2) + qJD(3);
t680 = (t751 * t775 - t725) * pkin(8) + (t751 * t752 + t774) * pkin(3) + t700;
t701 = t779 * t718 + t783 * t719;
t745 = t775 * pkin(3) - t752 * pkin(8);
t747 = t751 ^ 2;
t682 = -t747 * pkin(3) + t724 * pkin(8) - t775 * t745 + t701;
t674 = t782 * t680 - t778 * t682;
t771 = qJDD(4) + t774;
t772 = qJD(4) + t775;
t670 = -0.2e1 * qJD(5) * t739 + (t738 * t772 - t695) * qJ(5) + (t738 * t739 + t771) * pkin(4) + t674;
t727 = -t772 * mrSges(6,2) + t738 * mrSges(6,3);
t802 = m(6) * t670 + t771 * mrSges(6,1) + t772 * t727;
t666 = -t695 * mrSges(6,3) - t739 * t712 + t802;
t675 = t778 * t680 + t782 * t682;
t694 = -t739 * qJD(4) + t782 * t724 - t778 * t725;
t729 = t772 * pkin(4) - t739 * qJ(5);
t737 = t738 ^ 2;
t672 = -t737 * pkin(4) + t694 * qJ(5) + 0.2e1 * qJD(5) * t738 - t772 * t729 + t675;
t817 = Ifges(5,5) + Ifges(6,5);
t818 = Ifges(5,1) + Ifges(6,1);
t806 = t810 * t738 + t818 * t739 + t817 * t772;
t814 = t821 * t738 + t810 * t739 + t816 * t772;
t815 = Ifges(5,3) + Ifges(6,3);
t820 = mrSges(5,1) * t674 + mrSges(6,1) * t670 - mrSges(5,2) * t675 - mrSges(6,2) * t672 + pkin(4) * t666 + t816 * t694 + t817 * t695 - t806 * t738 + t814 * t739 + t815 * t771;
t713 = -t738 * mrSges(5,1) + t739 * mrSges(5,2);
t728 = -t772 * mrSges(5,2) + t738 * mrSges(5,3);
t658 = m(5) * t674 + t771 * mrSges(5,1) + t772 * t728 + (-t712 - t713) * t739 + (-mrSges(5,3) - mrSges(6,3)) * t695 + t802;
t730 = t772 * mrSges(6,1) - t739 * mrSges(6,3);
t731 = t772 * mrSges(5,1) - t739 * mrSges(5,3);
t801 = m(6) * t672 + t694 * mrSges(6,3) + t738 * t712;
t663 = m(5) * t675 + t694 * mrSges(5,3) + t738 * t713 + (-t730 - t731) * t772 + (-mrSges(5,2) - mrSges(6,2)) * t771 + t801;
t656 = t782 * t658 + t778 * t663;
t733 = Ifges(4,4) * t752 + Ifges(4,2) * t751 + Ifges(4,6) * t775;
t734 = Ifges(4,1) * t752 + Ifges(4,4) * t751 + Ifges(4,5) * t775;
t819 = mrSges(4,1) * t700 - mrSges(4,2) * t701 + Ifges(4,5) * t725 + Ifges(4,6) * t724 + Ifges(4,3) * t774 + pkin(3) * t656 + t752 * t733 - t751 * t734 + t820;
t740 = -t751 * mrSges(4,1) + t752 * mrSges(4,2);
t743 = -t775 * mrSges(4,2) + t751 * mrSges(4,3);
t652 = m(4) * t700 + t774 * mrSges(4,1) - t725 * mrSges(4,3) - t752 * t740 + t775 * t743 + t656;
t744 = t775 * mrSges(4,1) - t752 * mrSges(4,3);
t796 = -t778 * t658 + t782 * t663;
t653 = m(4) * t701 - t774 * mrSges(4,2) + t724 * mrSges(4,3) + t751 * t740 - t775 * t744 + t796;
t647 = t783 * t652 + t779 * t653;
t741 = -t784 * g(3) - t809;
t749 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t780 + Ifges(3,2) * t784) * qJD(1);
t750 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t780 + Ifges(3,4) * t784) * qJD(1);
t813 = mrSges(3,1) * t741 - mrSges(3,2) * t742 + Ifges(3,5) * t759 + Ifges(3,6) * t760 + Ifges(3,3) * qJDD(2) + pkin(2) * t647 + (t780 * t749 - t784 * t750) * qJD(1) + t819;
t758 = (-mrSges(3,1) * t784 + mrSges(3,2) * t780) * qJD(1);
t804 = qJD(1) * t784;
t763 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t804;
t645 = m(3) * t741 + qJDD(2) * mrSges(3,1) - t759 * mrSges(3,3) + qJD(2) * t763 - t758 * t805 + t647;
t762 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t805;
t797 = -t779 * t652 + t783 * t653;
t646 = m(3) * t742 - qJDD(2) * mrSges(3,2) + t760 * mrSges(3,3) - qJD(2) * t762 + t758 * t804 + t797;
t798 = -t780 * t645 + t784 * t646;
t637 = m(2) * t766 - t786 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t798;
t765 = t781 * g(1) - t785 * g(2);
t794 = -qJDD(1) * pkin(1) - t765;
t753 = -t786 * pkin(6) + t794;
t726 = -t760 * pkin(2) + t764 * t805 + (-pkin(7) * t777 - pkin(6)) * t786 + t794;
t684 = -t724 * pkin(3) - t747 * pkin(8) + t752 * t745 + t726;
t677 = -t694 * pkin(4) - t737 * qJ(5) + t739 * t729 + qJDD(5) + t684;
t667 = m(6) * t677 - t694 * mrSges(6,1) + t695 * mrSges(6,2) - t738 * t727 + t739 * t730;
t792 = m(5) * t684 - t694 * mrSges(5,1) + t695 * mrSges(5,2) - t738 * t728 + t739 * t731 + t667;
t790 = m(4) * t726 - t724 * mrSges(4,1) + t725 * mrSges(4,2) - t751 * t743 + t752 * t744 + t792;
t788 = -m(3) * t753 + t760 * mrSges(3,1) - t759 * mrSges(3,2) - t762 * t805 + t763 * t804 - t790;
t660 = m(2) * t765 + qJDD(1) * mrSges(2,1) - t786 * mrSges(2,2) + t788;
t808 = t781 * t637 + t785 * t660;
t639 = t784 * t645 + t780 * t646;
t807 = -t816 * t738 - t817 * t739 - t815 * t772;
t799 = t785 * t637 - t781 * t660;
t648 = -mrSges(5,1) * t684 + mrSges(5,3) * t675 - mrSges(6,1) * t677 + mrSges(6,3) * t672 - pkin(4) * t667 + qJ(5) * t801 + (-qJ(5) * t730 + t806) * t772 + (-qJ(5) * mrSges(6,2) + t816) * t771 + t807 * t739 + t810 * t695 + t821 * t694;
t654 = mrSges(5,2) * t684 + mrSges(6,2) * t677 - mrSges(5,3) * t674 - mrSges(6,3) * t670 - qJ(5) * t666 + t810 * t694 + t818 * t695 - t807 * t738 + t817 * t771 - t814 * t772;
t732 = Ifges(4,5) * t752 + Ifges(4,6) * t751 + Ifges(4,3) * t775;
t640 = -mrSges(4,1) * t726 + mrSges(4,3) * t701 + Ifges(4,4) * t725 + Ifges(4,2) * t724 + Ifges(4,6) * t774 - pkin(3) * t792 + pkin(8) * t796 + t782 * t648 + t778 * t654 - t752 * t732 + t775 * t734;
t641 = mrSges(4,2) * t726 - mrSges(4,3) * t700 + Ifges(4,1) * t725 + Ifges(4,4) * t724 + Ifges(4,5) * t774 - pkin(8) * t656 - t778 * t648 + t782 * t654 + t751 * t732 - t775 * t733;
t748 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t780 + Ifges(3,6) * t784) * qJD(1);
t631 = -mrSges(3,1) * t753 + mrSges(3,3) * t742 + Ifges(3,4) * t759 + Ifges(3,2) * t760 + Ifges(3,6) * qJDD(2) - pkin(2) * t790 + pkin(7) * t797 + qJD(2) * t750 + t783 * t640 + t779 * t641 - t748 * t805;
t633 = mrSges(3,2) * t753 - mrSges(3,3) * t741 + Ifges(3,1) * t759 + Ifges(3,4) * t760 + Ifges(3,5) * qJDD(2) - pkin(7) * t647 - qJD(2) * t749 - t779 * t640 + t783 * t641 + t748 * t804;
t793 = mrSges(2,1) * t765 - mrSges(2,2) * t766 + Ifges(2,3) * qJDD(1) + pkin(1) * t788 + pkin(6) * t798 + t784 * t631 + t780 * t633;
t634 = mrSges(2,1) * g(3) + mrSges(2,3) * t766 + t786 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t639 - t813;
t629 = -mrSges(2,2) * g(3) - mrSges(2,3) * t765 + Ifges(2,5) * qJDD(1) - t786 * Ifges(2,6) - pkin(6) * t639 - t780 * t631 + t784 * t633;
t1 = [-m(1) * g(1) + t799; -m(1) * g(2) + t808; (-m(1) - m(2)) * g(3) + t639; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t808 + t785 * t629 - t781 * t634; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t799 + t781 * t629 + t785 * t634; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t793; t793; t813; t819; t820; t667;];
tauJB = t1;
