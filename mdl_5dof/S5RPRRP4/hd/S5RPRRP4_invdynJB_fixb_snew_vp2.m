% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:01
% EndTime: 2019-12-05 18:06:09
% DurationCPUTime: 4.25s
% Computational Cost: add. (35957->268), mult. (86495->338), div. (0->0), fcn. (56474->8), ass. (0->118)
t802 = Ifges(5,4) + Ifges(6,4);
t812 = Ifges(5,2) + Ifges(6,2);
t807 = Ifges(5,6) + Ifges(6,6);
t761 = sin(qJ(1));
t764 = cos(qJ(1));
t742 = t761 * g(2) - t764 * g(3);
t765 = qJD(1) ^ 2;
t811 = -t765 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t742;
t759 = sin(qJ(4));
t760 = sin(qJ(3));
t762 = cos(qJ(4));
t763 = cos(qJ(3));
t757 = sin(pkin(8));
t793 = t757 * qJD(1);
t720 = (-t759 * t763 - t760 * t762) * t793;
t790 = qJD(1) * qJD(3);
t729 = (-qJDD(1) * t760 - t763 * t790) * t757;
t730 = (qJDD(1) * t763 - t760 * t790) * t757;
t689 = t720 * qJD(4) + t759 * t729 + t762 * t730;
t721 = (-t759 * t760 + t762 * t763) * t793;
t702 = -t720 * mrSges(6,1) + t721 * mrSges(6,2);
t758 = cos(pkin(8));
t713 = -t757 * g(1) + t811 * t758;
t777 = -pkin(2) * t758 - pkin(6) * t757;
t736 = t777 * qJD(1);
t792 = t758 * qJD(1);
t701 = t736 * t792 + t713;
t743 = t764 * g(2) + t761 * g(3);
t772 = -t765 * qJ(2) + qJDD(2) - t743;
t714 = (-pkin(1) + t777) * qJDD(1) + t772;
t711 = t763 * t714;
t789 = t758 * qJDD(1);
t745 = qJDD(3) - t789;
t746 = qJD(3) - t792;
t798 = t757 ^ 2 * t765;
t672 = t745 * pkin(3) - t730 * pkin(7) + t711 + (-pkin(3) * t763 * t798 - pkin(7) * t746 * t793 - t701) * t760;
t676 = t763 * t701 + t760 * t714;
t784 = t763 * t793;
t728 = t746 * pkin(3) - pkin(7) * t784;
t788 = t760 ^ 2 * t798;
t673 = -pkin(3) * t788 + t729 * pkin(7) - t746 * t728 + t676;
t665 = t762 * t672 - t759 * t673;
t741 = qJDD(4) + t745;
t744 = qJD(4) + t746;
t661 = -0.2e1 * qJD(5) * t721 + (t720 * t744 - t689) * qJ(5) + (t720 * t721 + t741) * pkin(4) + t665;
t705 = -t744 * mrSges(6,2) + t720 * mrSges(6,3);
t787 = m(6) * t661 + t741 * mrSges(6,1) + t744 * t705;
t657 = -t689 * mrSges(6,3) - t721 * t702 + t787;
t666 = t759 * t672 + t762 * t673;
t688 = -t721 * qJD(4) + t762 * t729 - t759 * t730;
t707 = t744 * pkin(4) - t721 * qJ(5);
t719 = t720 ^ 2;
t663 = -t719 * pkin(4) + t688 * qJ(5) + 0.2e1 * qJD(5) * t720 - t744 * t707 + t666;
t808 = Ifges(5,5) + Ifges(6,5);
t809 = Ifges(5,1) + Ifges(6,1);
t796 = -t802 * t720 - t809 * t721 - t808 * t744;
t805 = t812 * t720 + t802 * t721 + t807 * t744;
t806 = Ifges(5,3) + Ifges(6,3);
t810 = mrSges(5,1) * t665 + mrSges(6,1) * t661 - mrSges(5,2) * t666 - mrSges(6,2) * t663 + pkin(4) * t657 + t807 * t688 + t808 * t689 + t796 * t720 + t805 * t721 + t806 * t741;
t703 = -t720 * mrSges(5,1) + t721 * mrSges(5,2);
t706 = -t744 * mrSges(5,2) + t720 * mrSges(5,3);
t650 = m(5) * t665 + t741 * mrSges(5,1) + t744 * t706 + (-t702 - t703) * t721 + (-mrSges(5,3) - mrSges(6,3)) * t689 + t787;
t708 = t744 * mrSges(6,1) - t721 * mrSges(6,3);
t709 = t744 * mrSges(5,1) - t721 * mrSges(5,3);
t786 = m(6) * t663 + t688 * mrSges(6,3) + t720 * t702;
t653 = m(5) * t666 + t688 * mrSges(5,3) + t720 * t703 + (-t708 - t709) * t744 + (-mrSges(5,2) - mrSges(6,2)) * t741 + t786;
t648 = t762 * t650 + t759 * t653;
t675 = -t760 * t701 + t711;
t804 = -mrSges(4,1) * t675 + mrSges(4,2) * t676 - Ifges(4,5) * t730 - Ifges(4,6) * t729 - Ifges(4,3) * t745 - pkin(3) * t648 - t810;
t712 = -t758 * g(1) - t811 * t757;
t785 = t760 * t793;
t725 = -t746 * mrSges(4,2) - mrSges(4,3) * t785;
t726 = t746 * mrSges(4,1) - mrSges(4,3) * t784;
t803 = -t725 * t760 - t726 * t763;
t801 = mrSges(3,2) * t757;
t734 = (-mrSges(3,1) * t758 + t801) * qJD(1);
t727 = (mrSges(4,1) * t760 + mrSges(4,2) * t763) * t793;
t645 = m(4) * t675 + t745 * mrSges(4,1) - t730 * mrSges(4,3) + t746 * t725 - t727 * t784 + t648;
t778 = -t759 * t650 + t762 * t653;
t646 = m(4) * t676 - t745 * mrSges(4,2) + t729 * mrSges(4,3) - t746 * t726 - t727 * t785 + t778;
t779 = -t760 * t645 + t763 * t646;
t794 = qJDD(1) * mrSges(3,3);
t639 = m(3) * t713 + (qJD(1) * t734 + t794) * t758 + t779;
t700 = t736 * t793 - t712;
t674 = -t729 * pkin(3) - pkin(7) * t788 + t728 * t784 + t700;
t668 = -t688 * pkin(4) - t719 * qJ(5) + t721 * t707 + qJDD(5) + t674;
t658 = m(6) * t668 - t688 * mrSges(6,1) + t689 * mrSges(6,2) - t720 * t705 + t721 * t708;
t769 = m(5) * t674 - t688 * mrSges(5,1) + t689 * mrSges(5,2) - t720 * t706 + t721 * t709 + t658;
t767 = -m(4) * t700 + t729 * mrSges(4,1) - t730 * mrSges(4,2) - t769;
t655 = t767 + (-t794 + (-t734 + t803) * qJD(1)) * t757 + m(3) * t712;
t635 = t757 * t639 + t758 * t655;
t797 = -t807 * t720 - t808 * t721 - t806 * t744;
t780 = t758 * t639 - t757 * t655;
t633 = m(2) * t742 - t765 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t780;
t642 = t763 * t645 + t760 * t646;
t732 = -qJDD(1) * pkin(1) + t772;
t770 = -m(3) * t732 + mrSges(3,1) * t789 - t642 + (t758 ^ 2 * t765 + t798) * mrSges(3,3);
t636 = m(2) * t743 - t765 * mrSges(2,2) + (mrSges(2,1) - t801) * qJDD(1) + t770;
t781 = t764 * t633 - t761 * t636;
t776 = Ifges(3,1) * t757 + Ifges(3,4) * t758;
t775 = Ifges(3,5) * t757 + Ifges(3,6) * t758;
t774 = -t761 * t633 - t764 * t636;
t716 = Ifges(4,6) * t746 + (Ifges(4,4) * t763 - Ifges(4,2) * t760) * t793;
t717 = Ifges(4,5) * t746 + (Ifges(4,1) * t763 - Ifges(4,4) * t760) * t793;
t773 = t716 * t763 + t717 * t760;
t643 = -mrSges(5,1) * t674 + mrSges(5,3) * t666 - mrSges(6,1) * t668 + mrSges(6,3) * t663 - pkin(4) * t658 + qJ(5) * t786 + (-qJ(5) * t708 - t796) * t744 + (-qJ(5) * mrSges(6,2) + t807) * t741 + t797 * t721 + t802 * t689 + t812 * t688;
t647 = mrSges(5,2) * t674 + mrSges(6,2) * t668 - mrSges(5,3) * t665 - mrSges(6,3) * t661 - qJ(5) * t657 + t802 * t688 + t809 * t689 - t797 * t720 + t808 * t741 - t805 * t744;
t715 = Ifges(4,3) * t746 + (Ifges(4,5) * t763 - Ifges(4,6) * t760) * t793;
t630 = -mrSges(4,1) * t700 + mrSges(4,3) * t676 + Ifges(4,4) * t730 + Ifges(4,2) * t729 + Ifges(4,6) * t745 - pkin(3) * t769 + pkin(7) * t778 + t762 * t643 + t759 * t647 - t715 * t784 + t746 * t717;
t631 = mrSges(4,2) * t700 - mrSges(4,3) * t675 + Ifges(4,1) * t730 + Ifges(4,4) * t729 + Ifges(4,5) * t745 - pkin(7) * t648 - t759 * t643 + t762 * t647 - t715 * t785 - t746 * t716;
t735 = t775 * qJD(1);
t627 = mrSges(3,2) * t732 - mrSges(3,3) * t712 - pkin(6) * t642 + t776 * qJDD(1) - t760 * t630 + t763 * t631 + t735 * t792;
t629 = Ifges(3,2) * t789 + (Ifges(3,4) * qJDD(1) + (-t735 - t773) * qJD(1)) * t757 - mrSges(3,1) * t732 + mrSges(3,3) * t713 - pkin(2) * t642 + t804;
t641 = qJDD(1) * t801 - t770;
t771 = mrSges(2,1) * t743 - mrSges(2,2) * t742 + Ifges(2,3) * qJDD(1) - pkin(1) * t641 + qJ(2) * t780 + t757 * t627 + t758 * t629;
t625 = t765 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t742 - mrSges(3,1) * t712 + mrSges(3,2) * t713 - t760 * t631 - t763 * t630 - pkin(2) * t767 - pkin(6) * t779 - pkin(1) * t635 + (Ifges(2,6) - t775) * qJDD(1) + (-pkin(2) * t803 * t757 + (-t757 * (Ifges(3,4) * t757 + Ifges(3,2) * t758) + t758 * t776) * qJD(1)) * qJD(1);
t624 = -mrSges(2,2) * g(1) - mrSges(2,3) * t743 + Ifges(2,5) * qJDD(1) - t765 * Ifges(2,6) - qJ(2) * t635 + t758 * t627 - t757 * t629;
t1 = [(-m(1) - m(2)) * g(1) + t635; -m(1) * g(2) + t774; -m(1) * g(3) + t781; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t771; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t781 - t761 * t624 - t764 * t625; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t774 + t764 * t624 - t761 * t625; t771; t641; t773 * t793 - t804; t810; t658;];
tauJB = t1;
