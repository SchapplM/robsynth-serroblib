% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:36:52
% EndTime: 2019-12-31 20:37:07
% DurationCPUTime: 14.82s
% Computational Cost: add. (234849->326), mult. (536655->428), div. (0->0), fcn. (420103->12), ass. (0->137)
t779 = sin(pkin(5));
t810 = pkin(7) * t779;
t781 = cos(pkin(5));
t809 = t781 * g(3);
t784 = sin(qJ(2));
t808 = t779 * t784;
t788 = cos(qJ(2));
t807 = t779 * t788;
t806 = t781 * t784;
t805 = t781 * t788;
t785 = sin(qJ(1));
t789 = cos(qJ(1));
t769 = t785 * g(1) - t789 * g(2);
t790 = qJD(1) ^ 2;
t760 = qJDD(1) * pkin(1) + t790 * t810 + t769;
t770 = -t789 * g(1) - t785 * g(2);
t800 = qJDD(1) * t779;
t761 = -t790 * pkin(1) + pkin(7) * t800 + t770;
t803 = t760 * t806 + t788 * t761;
t733 = -g(3) * t808 + t803;
t774 = t781 * qJD(1) + qJD(2);
t802 = qJD(1) * t779;
t799 = t784 * t802;
t758 = t774 * mrSges(3,1) - mrSges(3,3) * t799;
t763 = (-mrSges(3,1) * t788 + mrSges(3,2) * t784) * t802;
t765 = -qJD(2) * t799 + t788 * t800;
t773 = t781 * qJDD(1) + qJDD(2);
t762 = (-pkin(2) * t788 - qJ(3) * t784) * t802;
t772 = t774 ^ 2;
t801 = qJD(1) * t788;
t721 = -t772 * pkin(2) + t773 * qJ(3) + (-g(3) * t784 + t762 * t801) * t779 + t803;
t764 = (qJD(2) * t801 + qJDD(1) * t784) * t779;
t722 = -t765 * pkin(2) - t809 - t764 * qJ(3) + (-t760 + (pkin(2) * t784 - qJ(3) * t788) * t774 * qJD(1)) * t779;
t778 = sin(pkin(10));
t780 = cos(pkin(10));
t754 = t778 * t774 + t780 * t799;
t694 = -0.2e1 * qJD(3) * t754 - t778 * t721 + t780 * t722;
t742 = t780 * t764 + t778 * t773;
t753 = t780 * t774 - t778 * t799;
t798 = t779 * t801;
t691 = (-t753 * t798 - t742) * pkin(8) + (t753 * t754 - t765) * pkin(3) + t694;
t695 = 0.2e1 * qJD(3) * t753 + t780 * t721 + t778 * t722;
t741 = -t778 * t764 + t780 * t773;
t743 = -pkin(3) * t798 - t754 * pkin(8);
t752 = t753 ^ 2;
t693 = -t752 * pkin(3) + t741 * pkin(8) + t743 * t798 + t695;
t783 = sin(qJ(4));
t787 = cos(qJ(4));
t688 = t783 * t691 + t787 * t693;
t735 = t787 * t753 - t783 * t754;
t736 = t783 * t753 + t787 * t754;
t716 = -t735 * pkin(4) - t736 * pkin(9);
t757 = qJDD(4) - t765;
t768 = qJD(4) - t798;
t767 = t768 ^ 2;
t686 = -t767 * pkin(4) + t757 * pkin(9) + t735 * t716 + t688;
t732 = -g(3) * t807 + t760 * t805 - t784 * t761;
t720 = -t773 * pkin(2) - t772 * qJ(3) + t762 * t799 + qJDD(3) - t732;
t699 = -t741 * pkin(3) - t752 * pkin(8) + t754 * t743 + t720;
t707 = -t736 * qJD(4) + t787 * t741 - t783 * t742;
t708 = t735 * qJD(4) + t783 * t741 + t787 * t742;
t689 = (-t735 * t768 - t708) * pkin(9) + (t736 * t768 - t707) * pkin(4) + t699;
t782 = sin(qJ(5));
t786 = cos(qJ(5));
t683 = -t782 * t686 + t786 * t689;
t723 = -t782 * t736 + t786 * t768;
t698 = t723 * qJD(5) + t786 * t708 + t782 * t757;
t724 = t786 * t736 + t782 * t768;
t704 = -t723 * mrSges(6,1) + t724 * mrSges(6,2);
t706 = qJDD(5) - t707;
t734 = qJD(5) - t735;
t709 = -t734 * mrSges(6,2) + t723 * mrSges(6,3);
t680 = m(6) * t683 + t706 * mrSges(6,1) - t698 * mrSges(6,3) - t724 * t704 + t734 * t709;
t684 = t786 * t686 + t782 * t689;
t697 = -t724 * qJD(5) - t782 * t708 + t786 * t757;
t710 = t734 * mrSges(6,1) - t724 * mrSges(6,3);
t681 = m(6) * t684 - t706 * mrSges(6,2) + t697 * mrSges(6,3) + t723 * t704 - t734 * t710;
t672 = -t782 * t680 + t786 * t681;
t715 = -t735 * mrSges(5,1) + t736 * mrSges(5,2);
t726 = t768 * mrSges(5,1) - t736 * mrSges(5,3);
t669 = m(5) * t688 - t757 * mrSges(5,2) + t707 * mrSges(5,3) + t735 * t715 - t768 * t726 + t672;
t687 = t787 * t691 - t783 * t693;
t685 = -t757 * pkin(4) - t767 * pkin(9) + t736 * t716 - t687;
t682 = -m(6) * t685 + t697 * mrSges(6,1) - t698 * mrSges(6,2) + t723 * t709 - t724 * t710;
t725 = -t768 * mrSges(5,2) + t735 * mrSges(5,3);
t676 = m(5) * t687 + t757 * mrSges(5,1) - t708 * mrSges(5,3) - t736 * t715 + t768 * t725 + t682;
t663 = t783 * t669 + t787 * t676;
t737 = -t753 * mrSges(4,1) + t754 * mrSges(4,2);
t739 = mrSges(4,2) * t798 + t753 * mrSges(4,3);
t661 = m(4) * t694 - t765 * mrSges(4,1) - t742 * mrSges(4,3) - t754 * t737 - t739 * t798 + t663;
t740 = -mrSges(4,1) * t798 - t754 * mrSges(4,3);
t795 = t787 * t669 - t783 * t676;
t662 = m(4) * t695 + t765 * mrSges(4,2) + t741 * mrSges(4,3) + t753 * t737 + t740 * t798 + t795;
t796 = -t778 * t661 + t780 * t662;
t652 = m(3) * t733 - t773 * mrSges(3,2) + t765 * mrSges(3,3) - t774 * t758 + t763 * t798 + t796;
t655 = t780 * t661 + t778 * t662;
t747 = -t779 * t760 - t809;
t759 = -t774 * mrSges(3,2) + mrSges(3,3) * t798;
t654 = m(3) * t747 - t765 * mrSges(3,1) + t764 * mrSges(3,2) + (t758 * t784 - t759 * t788) * t802 + t655;
t671 = t786 * t680 + t782 * t681;
t793 = m(5) * t699 - t707 * mrSges(5,1) + t708 * mrSges(5,2) - t735 * t725 + t736 * t726 + t671;
t670 = m(4) * t720 - t741 * mrSges(4,1) + t742 * mrSges(4,2) - t753 * t739 + t754 * t740 + t793;
t666 = m(3) * t732 + t773 * mrSges(3,1) - t764 * mrSges(3,3) + t774 * t759 - t763 * t799 - t670;
t641 = t652 * t806 - t779 * t654 + t666 * t805;
t638 = m(2) * t769 + qJDD(1) * mrSges(2,1) - t790 * mrSges(2,2) + t641;
t648 = t788 * t652 - t784 * t666;
t646 = m(2) * t770 - t790 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t648;
t804 = t789 * t638 + t785 * t646;
t640 = t652 * t808 + t781 * t654 + t666 * t807;
t797 = -t785 * t638 + t789 * t646;
t700 = Ifges(6,5) * t724 + Ifges(6,6) * t723 + Ifges(6,3) * t734;
t702 = Ifges(6,1) * t724 + Ifges(6,4) * t723 + Ifges(6,5) * t734;
t673 = -mrSges(6,1) * t685 + mrSges(6,3) * t684 + Ifges(6,4) * t698 + Ifges(6,2) * t697 + Ifges(6,6) * t706 - t724 * t700 + t734 * t702;
t701 = Ifges(6,4) * t724 + Ifges(6,2) * t723 + Ifges(6,6) * t734;
t674 = mrSges(6,2) * t685 - mrSges(6,3) * t683 + Ifges(6,1) * t698 + Ifges(6,4) * t697 + Ifges(6,5) * t706 + t723 * t700 - t734 * t701;
t711 = Ifges(5,5) * t736 + Ifges(5,6) * t735 + Ifges(5,3) * t768;
t712 = Ifges(5,4) * t736 + Ifges(5,2) * t735 + Ifges(5,6) * t768;
t656 = mrSges(5,2) * t699 - mrSges(5,3) * t687 + Ifges(5,1) * t708 + Ifges(5,4) * t707 + Ifges(5,5) * t757 - pkin(9) * t671 - t782 * t673 + t786 * t674 + t735 * t711 - t768 * t712;
t713 = Ifges(5,1) * t736 + Ifges(5,4) * t735 + Ifges(5,5) * t768;
t792 = mrSges(6,1) * t683 - mrSges(6,2) * t684 + Ifges(6,5) * t698 + Ifges(6,6) * t697 + Ifges(6,3) * t706 + t724 * t701 - t723 * t702;
t657 = -mrSges(5,1) * t699 + mrSges(5,3) * t688 + Ifges(5,4) * t708 + Ifges(5,2) * t707 + Ifges(5,6) * t757 - pkin(4) * t671 - t736 * t711 + t768 * t713 - t792;
t727 = Ifges(4,5) * t754 + Ifges(4,6) * t753 - Ifges(4,3) * t798;
t729 = Ifges(4,1) * t754 + Ifges(4,4) * t753 - Ifges(4,5) * t798;
t642 = -mrSges(4,1) * t720 + mrSges(4,3) * t695 + Ifges(4,4) * t742 + Ifges(4,2) * t741 - Ifges(4,6) * t765 - pkin(3) * t793 + pkin(8) * t795 + t783 * t656 + t787 * t657 - t754 * t727 - t729 * t798;
t728 = Ifges(4,4) * t754 + Ifges(4,2) * t753 - Ifges(4,6) * t798;
t643 = mrSges(4,2) * t720 - mrSges(4,3) * t694 + Ifges(4,1) * t742 + Ifges(4,4) * t741 - Ifges(4,5) * t765 - pkin(8) * t663 + t787 * t656 - t783 * t657 + t753 * t727 + t728 * t798;
t745 = Ifges(3,6) * t774 + (Ifges(3,4) * t784 + Ifges(3,2) * t788) * t802;
t746 = Ifges(3,5) * t774 + (Ifges(3,1) * t784 + Ifges(3,4) * t788) * t802;
t632 = Ifges(3,5) * t764 + Ifges(3,6) * t765 + Ifges(3,3) * t773 + mrSges(3,1) * t732 - mrSges(3,2) * t733 + t778 * t643 + t780 * t642 - pkin(2) * t670 + qJ(3) * t796 + (t745 * t784 - t746 * t788) * t802;
t744 = Ifges(3,3) * t774 + (Ifges(3,5) * t784 + Ifges(3,6) * t788) * t802;
t634 = mrSges(3,2) * t747 - mrSges(3,3) * t732 + Ifges(3,1) * t764 + Ifges(3,4) * t765 + Ifges(3,5) * t773 - qJ(3) * t655 - t778 * t642 + t780 * t643 + t744 * t798 - t774 * t745;
t791 = mrSges(5,1) * t687 - mrSges(5,2) * t688 + Ifges(5,5) * t708 + Ifges(5,6) * t707 + Ifges(5,3) * t757 + pkin(4) * t682 + pkin(9) * t672 + t786 * t673 + t782 * t674 + t736 * t712 - t735 * t713;
t636 = -t744 * t799 - t791 + (Ifges(3,2) + Ifges(4,3)) * t765 + t774 * t746 + Ifges(3,6) * t773 + Ifges(3,4) * t764 + t753 * t729 - t754 * t728 - Ifges(4,6) * t741 - Ifges(4,5) * t742 - mrSges(3,1) * t747 + mrSges(3,3) * t733 - mrSges(4,1) * t694 + mrSges(4,2) * t695 - pkin(3) * t663 - pkin(2) * t655;
t794 = mrSges(2,1) * t769 - mrSges(2,2) * t770 + Ifges(2,3) * qJDD(1) + pkin(1) * t641 + t781 * t632 + t634 * t808 + t636 * t807 + t648 * t810;
t630 = -mrSges(2,2) * g(3) - mrSges(2,3) * t769 + Ifges(2,5) * qJDD(1) - t790 * Ifges(2,6) + t788 * t634 - t784 * t636 + (-t640 * t779 - t641 * t781) * pkin(7);
t629 = mrSges(2,1) * g(3) + mrSges(2,3) * t770 + t790 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t640 - t779 * t632 + (pkin(7) * t648 + t634 * t784 + t636 * t788) * t781;
t1 = [-m(1) * g(1) + t797; -m(1) * g(2) + t804; (-m(1) - m(2)) * g(3) + t640; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t804 - t785 * t629 + t789 * t630; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t797 + t789 * t629 + t785 * t630; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t794; t794; t632; t670; t791; t792;];
tauJB = t1;
