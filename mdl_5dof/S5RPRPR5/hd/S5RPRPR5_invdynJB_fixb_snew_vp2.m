% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:42:27
% EndTime: 2020-01-03 11:42:38
% DurationCPUTime: 8.30s
% Computational Cost: add. (80930->290), mult. (207955->381), div. (0->0), fcn. (141748->10), ass. (0->126)
t765 = sin(qJ(1));
t768 = cos(qJ(1));
t742 = -t765 * g(2) + t768 * g(3);
t769 = qJD(1) ^ 2;
t805 = -t769 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t742;
t760 = sin(pkin(8));
t762 = cos(pkin(8));
t712 = -t762 * g(1) - t805 * t760;
t793 = t762 * qJD(1);
t746 = qJD(3) - t793;
t764 = sin(qJ(3));
t794 = t760 * qJD(1);
t788 = t764 * t794;
t727 = -t746 * mrSges(4,2) - mrSges(4,3) * t788;
t767 = cos(qJ(3));
t787 = t767 * t794;
t729 = t746 * mrSges(4,1) - mrSges(4,3) * t787;
t804 = -t727 * t764 - t729 * t767;
t713 = -t760 * g(1) + t805 * t762;
t780 = -pkin(2) * t762 - pkin(6) * t760;
t738 = t780 * qJD(1);
t702 = t738 * t793 + t713;
t743 = -t768 * g(2) - t765 * g(3);
t775 = -t769 * qJ(2) + qJDD(2) - t743;
t714 = (-pkin(1) + t780) * qJDD(1) + t775;
t711 = t767 * t714;
t791 = qJD(1) * qJD(3);
t732 = (qJDD(1) * t767 - t764 * t791) * t760;
t790 = t762 * qJDD(1);
t745 = qJDD(3) - t790;
t798 = t760 ^ 2 * t769;
t681 = t745 * pkin(3) - t732 * qJ(4) + t711 + (-pkin(3) * t767 * t798 - qJ(4) * t746 * t794 - t702) * t764;
t691 = t767 * t702 + t764 * t714;
t728 = t746 * pkin(3) - qJ(4) * t787;
t731 = (-qJDD(1) * t764 - t767 * t791) * t760;
t789 = t764 ^ 2 * t798;
t682 = -pkin(3) * t789 + t731 * qJ(4) - t746 * t728 + t691;
t759 = sin(pkin(9));
t761 = cos(pkin(9));
t723 = (-t759 * t764 + t761 * t767) * t794;
t667 = -0.2e1 * qJD(4) * t723 + t761 * t681 - t759 * t682;
t706 = t759 * t731 + t761 * t732;
t722 = (-t759 * t767 - t761 * t764) * t794;
t665 = (t722 * t746 - t706) * pkin(7) + (t722 * t723 + t745) * pkin(4) + t667;
t668 = 0.2e1 * qJD(4) * t722 + t759 * t681 + t761 * t682;
t705 = t761 * t731 - t759 * t732;
t709 = t746 * pkin(4) - t723 * pkin(7);
t721 = t722 ^ 2;
t666 = -t721 * pkin(4) + t705 * pkin(7) - t746 * t709 + t668;
t763 = sin(qJ(5));
t766 = cos(qJ(5));
t663 = t766 * t665 - t763 * t666;
t699 = t766 * t722 - t763 * t723;
t677 = t699 * qJD(5) + t763 * t705 + t766 * t706;
t700 = t763 * t722 + t766 * t723;
t688 = -t699 * mrSges(6,1) + t700 * mrSges(6,2);
t744 = qJD(5) + t746;
t692 = -t744 * mrSges(6,2) + t699 * mrSges(6,3);
t741 = qJDD(5) + t745;
t659 = m(6) * t663 + t741 * mrSges(6,1) - t677 * mrSges(6,3) - t700 * t688 + t744 * t692;
t664 = t763 * t665 + t766 * t666;
t676 = -t700 * qJD(5) + t766 * t705 - t763 * t706;
t693 = t744 * mrSges(6,1) - t700 * mrSges(6,3);
t660 = m(6) * t664 - t741 * mrSges(6,2) + t676 * mrSges(6,3) + t699 * t688 - t744 * t693;
t651 = t766 * t659 + t763 * t660;
t703 = -t722 * mrSges(5,1) + t723 * mrSges(5,2);
t707 = -t746 * mrSges(5,2) + t722 * mrSges(5,3);
t649 = m(5) * t667 + t745 * mrSges(5,1) - t706 * mrSges(5,3) - t723 * t703 + t746 * t707 + t651;
t708 = t746 * mrSges(5,1) - t723 * mrSges(5,3);
t781 = -t763 * t659 + t766 * t660;
t650 = m(5) * t668 - t745 * mrSges(5,2) + t705 * mrSges(5,3) + t722 * t703 - t746 * t708 + t781;
t645 = t761 * t649 + t759 * t650;
t690 = -t764 * t702 + t711;
t695 = Ifges(5,4) * t723 + Ifges(5,2) * t722 + Ifges(5,6) * t746;
t696 = Ifges(5,1) * t723 + Ifges(5,4) * t722 + Ifges(5,5) * t746;
t684 = Ifges(6,4) * t700 + Ifges(6,2) * t699 + Ifges(6,6) * t744;
t685 = Ifges(6,1) * t700 + Ifges(6,4) * t699 + Ifges(6,5) * t744;
t773 = -mrSges(6,1) * t663 + mrSges(6,2) * t664 - Ifges(6,5) * t677 - Ifges(6,6) * t676 - Ifges(6,3) * t741 - t700 * t684 + t699 * t685;
t803 = -mrSges(4,1) * t690 - mrSges(5,1) * t667 + mrSges(4,2) * t691 + mrSges(5,2) * t668 - Ifges(4,5) * t732 - Ifges(5,5) * t706 - Ifges(4,6) * t731 - Ifges(5,6) * t705 - pkin(3) * t645 - pkin(4) * t651 - t723 * t695 + t722 * t696 - (Ifges(4,3) + Ifges(5,3)) * t745 + t773;
t801 = mrSges(3,2) * t760;
t736 = (-mrSges(3,1) * t762 + t801) * qJD(1);
t730 = (mrSges(4,1) * t764 + mrSges(4,2) * t767) * t794;
t643 = m(4) * t690 + t745 * mrSges(4,1) - t732 * mrSges(4,3) + t746 * t727 - t730 * t787 + t645;
t782 = -t759 * t649 + t761 * t650;
t644 = m(4) * t691 - t745 * mrSges(4,2) + t731 * mrSges(4,3) - t746 * t729 - t730 * t788 + t782;
t783 = -t764 * t643 + t767 * t644;
t795 = qJDD(1) * mrSges(3,3);
t636 = m(3) * t713 + (qJD(1) * t736 + t795) * t762 + t783;
t701 = t738 * t794 - t712;
t689 = -t731 * pkin(3) - qJ(4) * t789 + t728 * t787 + qJDD(4) + t701;
t670 = -t705 * pkin(4) - t721 * pkin(7) + t723 * t709 + t689;
t777 = m(6) * t670 - t676 * mrSges(6,1) + t677 * mrSges(6,2) - t699 * t692 + t700 * t693;
t661 = m(5) * t689 - t705 * mrSges(5,1) + t706 * mrSges(5,2) - t722 * t707 + t723 * t708 + t777;
t771 = -m(4) * t701 + t731 * mrSges(4,1) - t732 * mrSges(4,2) - t661;
t658 = t771 + (-t795 + (-t736 + t804) * qJD(1)) * t760 + m(3) * t712;
t784 = t762 * t636 - t760 * t658;
t628 = m(2) * t742 - t769 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t784;
t639 = t767 * t643 + t764 * t644;
t734 = -qJDD(1) * pkin(1) + t775;
t772 = -m(3) * t734 + mrSges(3,1) * t790 - t639 + (t762 ^ 2 * t769 + t798) * mrSges(3,3);
t633 = m(2) * t743 - t769 * mrSges(2,2) + (mrSges(2,1) - t801) * qJDD(1) + t772;
t797 = t765 * t628 + t768 * t633;
t630 = t760 * t636 + t762 * t658;
t785 = -t768 * t628 + t765 * t633;
t779 = Ifges(3,1) * t760 + Ifges(3,4) * t762;
t778 = Ifges(3,5) * t760 + Ifges(3,6) * t762;
t716 = Ifges(4,6) * t746 + (Ifges(4,4) * t767 - Ifges(4,2) * t764) * t794;
t717 = Ifges(4,5) * t746 + (Ifges(4,1) * t767 - Ifges(4,4) * t764) * t794;
t776 = t716 * t767 + t717 * t764;
t683 = Ifges(6,5) * t700 + Ifges(6,6) * t699 + Ifges(6,3) * t744;
t652 = -mrSges(6,1) * t670 + mrSges(6,3) * t664 + Ifges(6,4) * t677 + Ifges(6,2) * t676 + Ifges(6,6) * t741 - t700 * t683 + t744 * t685;
t653 = mrSges(6,2) * t670 - mrSges(6,3) * t663 + Ifges(6,1) * t677 + Ifges(6,4) * t676 + Ifges(6,5) * t741 + t699 * t683 - t744 * t684;
t694 = Ifges(5,5) * t723 + Ifges(5,6) * t722 + Ifges(5,3) * t746;
t640 = -mrSges(5,1) * t689 + mrSges(5,3) * t668 + Ifges(5,4) * t706 + Ifges(5,2) * t705 + Ifges(5,6) * t745 - pkin(4) * t777 + pkin(7) * t781 + t766 * t652 + t763 * t653 - t723 * t694 + t746 * t696;
t641 = mrSges(5,2) * t689 - mrSges(5,3) * t667 + Ifges(5,1) * t706 + Ifges(5,4) * t705 + Ifges(5,5) * t745 - pkin(7) * t651 - t763 * t652 + t766 * t653 + t722 * t694 - t746 * t695;
t715 = Ifges(4,3) * t746 + (Ifges(4,5) * t767 - Ifges(4,6) * t764) * t794;
t625 = -mrSges(4,1) * t701 + mrSges(4,3) * t691 + Ifges(4,4) * t732 + Ifges(4,2) * t731 + Ifges(4,6) * t745 - pkin(3) * t661 + qJ(4) * t782 + t761 * t640 + t759 * t641 - t715 * t787 + t746 * t717;
t626 = mrSges(4,2) * t701 - mrSges(4,3) * t690 + Ifges(4,1) * t732 + Ifges(4,4) * t731 + Ifges(4,5) * t745 - qJ(4) * t645 - t759 * t640 + t761 * t641 - t715 * t788 - t746 * t716;
t737 = t778 * qJD(1);
t622 = mrSges(3,2) * t734 - mrSges(3,3) * t712 - pkin(6) * t639 + t779 * qJDD(1) - t764 * t625 + t767 * t626 + t737 * t793;
t624 = (Ifges(3,4) * qJDD(1) + (-t737 - t776) * qJD(1)) * t760 - mrSges(3,1) * t734 + mrSges(3,3) * t713 - pkin(2) * t639 + Ifges(3,2) * t790 + t803;
t638 = qJDD(1) * t801 - t772;
t774 = mrSges(2,1) * t743 - mrSges(2,2) * t742 + Ifges(2,3) * qJDD(1) - pkin(1) * t638 + qJ(2) * t784 + t760 * t622 + t762 * t624;
t620 = t769 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t742 - mrSges(3,1) * t712 + mrSges(3,2) * t713 - t764 * t626 - t767 * t625 - pkin(2) * t771 - pkin(6) * t783 - pkin(1) * t630 + (Ifges(2,6) - t778) * qJDD(1) + (-pkin(2) * t804 * t760 + (-t760 * (Ifges(3,4) * t760 + Ifges(3,2) * t762) + t762 * t779) * qJD(1)) * qJD(1);
t619 = -mrSges(2,2) * g(1) - mrSges(2,3) * t743 + Ifges(2,5) * qJDD(1) - t769 * Ifges(2,6) - qJ(2) * t630 + t762 * t622 - t760 * t624;
t1 = [(-m(1) - m(2)) * g(1) + t630; -m(1) * g(2) + t797; -m(1) * g(3) + t785; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t774; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t785 + t765 * t619 + t768 * t620; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t797 - t768 * t619 + t765 * t620; t774; t638; t776 * t794 - t803; t661; -t773;];
tauJB = t1;
