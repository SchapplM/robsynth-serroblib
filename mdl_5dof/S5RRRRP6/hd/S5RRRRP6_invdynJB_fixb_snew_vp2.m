% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:53:18
% EndTime: 2019-12-31 21:53:25
% DurationCPUTime: 5.08s
% Computational Cost: add. (56552->293), mult. (114191->357), div. (0->0), fcn. (76745->8), ass. (0->119)
t784 = Ifges(5,1) + Ifges(6,1);
t777 = Ifges(5,4) + Ifges(6,4);
t776 = Ifges(5,5) + Ifges(6,5);
t783 = Ifges(5,2) + Ifges(6,2);
t775 = Ifges(5,6) + Ifges(6,6);
t782 = Ifges(5,3) + Ifges(6,3);
t740 = sin(qJ(3));
t741 = sin(qJ(2));
t744 = cos(qJ(3));
t745 = cos(qJ(2));
t716 = (t740 * t741 - t744 * t745) * qJD(1);
t765 = qJD(1) * qJD(2);
t724 = t741 * qJDD(1) + t745 * t765;
t742 = sin(qJ(1));
t746 = cos(qJ(1));
t731 = -t746 * g(1) - t742 * g(2);
t747 = qJD(1) ^ 2;
t719 = -t747 * pkin(1) + qJDD(1) * pkin(6) + t731;
t773 = t741 * t719;
t779 = pkin(2) * t747;
t682 = qJDD(2) * pkin(2) - t724 * pkin(7) - t773 + (pkin(7) * t765 + t741 * t779 - g(3)) * t745;
t708 = -t741 * g(3) + t745 * t719;
t725 = t745 * qJDD(1) - t741 * t765;
t767 = qJD(1) * t741;
t729 = qJD(2) * pkin(2) - pkin(7) * t767;
t738 = t745 ^ 2;
t683 = t725 * pkin(7) - qJD(2) * t729 - t738 * t779 + t708;
t660 = t740 * t682 + t744 * t683;
t717 = (t740 * t745 + t741 * t744) * qJD(1);
t690 = -t717 * qJD(3) - t740 * t724 + t744 * t725;
t702 = t716 * mrSges(4,1) + t717 * mrSges(4,2);
t736 = qJD(2) + qJD(3);
t710 = t736 * mrSges(4,1) - t717 * mrSges(4,3);
t735 = qJDD(2) + qJDD(3);
t691 = -t716 * qJD(3) + t744 * t724 + t740 * t725;
t730 = t742 * g(1) - t746 * g(2);
t755 = -qJDD(1) * pkin(1) - t730;
t692 = -t725 * pkin(2) + t729 * t767 + (-pkin(7) * t738 - pkin(6)) * t747 + t755;
t654 = (t716 * t736 - t691) * pkin(8) + (t717 * t736 - t690) * pkin(3) + t692;
t703 = t716 * pkin(3) - t717 * pkin(8);
t734 = t736 ^ 2;
t657 = -t734 * pkin(3) + t735 * pkin(8) - t716 * t703 + t660;
t739 = sin(qJ(4));
t743 = cos(qJ(4));
t649 = t743 * t654 - t739 * t657;
t705 = -t739 * t717 + t743 * t736;
t665 = t705 * qJD(4) + t743 * t691 + t739 * t735;
t706 = t743 * t717 + t739 * t736;
t679 = -t705 * mrSges(6,1) + t706 * mrSges(6,2);
t680 = -t705 * mrSges(5,1) + t706 * mrSges(5,2);
t689 = qJDD(4) - t690;
t712 = qJD(4) + t716;
t694 = -t712 * mrSges(5,2) + t705 * mrSges(5,3);
t646 = -0.2e1 * qJD(5) * t706 + (t705 * t712 - t665) * qJ(5) + (t705 * t706 + t689) * pkin(4) + t649;
t693 = -t712 * mrSges(6,2) + t705 * mrSges(6,3);
t764 = m(6) * t646 + t689 * mrSges(6,1) + t712 * t693;
t635 = m(5) * t649 + t689 * mrSges(5,1) + t712 * t694 + (-t679 - t680) * t706 + (-mrSges(5,3) - mrSges(6,3)) * t665 + t764;
t650 = t739 * t654 + t743 * t657;
t664 = -t706 * qJD(4) - t739 * t691 + t743 * t735;
t695 = t712 * pkin(4) - t706 * qJ(5);
t704 = t705 ^ 2;
t648 = -t704 * pkin(4) + t664 * qJ(5) + 0.2e1 * qJD(5) * t705 - t712 * t695 + t650;
t763 = m(6) * t648 + t664 * mrSges(6,3) + t705 * t679;
t696 = t712 * mrSges(6,1) - t706 * mrSges(6,3);
t768 = -t712 * mrSges(5,1) + t706 * mrSges(5,3) - t696;
t778 = -mrSges(5,2) - mrSges(6,2);
t638 = m(5) * t650 + t664 * mrSges(5,3) + t705 * t680 + t778 * t689 + t768 * t712 + t763;
t759 = -t739 * t635 + t743 * t638;
t628 = m(4) * t660 - t735 * mrSges(4,2) + t690 * mrSges(4,3) - t716 * t702 - t736 * t710 + t759;
t659 = t744 * t682 - t740 * t683;
t709 = -t736 * mrSges(4,2) - t716 * mrSges(4,3);
t656 = -t735 * pkin(3) - t734 * pkin(8) + t717 * t703 - t659;
t651 = -t664 * pkin(4) - t704 * qJ(5) + t706 * t695 + qJDD(5) + t656;
t758 = -m(6) * t651 + t664 * mrSges(6,1) + t705 * t693;
t750 = -m(5) * t656 + t664 * mrSges(5,1) + t778 * t665 + t705 * t694 + t768 * t706 + t758;
t640 = m(4) * t659 + t735 * mrSges(4,1) - t691 * mrSges(4,3) - t717 * t702 + t736 * t709 + t750;
t620 = t740 * t628 + t744 * t640;
t707 = -t745 * g(3) - t773;
t714 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t741 + Ifges(3,2) * t745) * qJD(1);
t715 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t741 + Ifges(3,4) * t745) * qJD(1);
t644 = t665 * mrSges(6,2) + t706 * t696 - t758;
t769 = t777 * t705 + t784 * t706 + t776 * t712;
t771 = -t775 * t705 - t776 * t706 - t782 * t712;
t622 = -mrSges(5,1) * t656 + mrSges(5,3) * t650 - mrSges(6,1) * t651 + mrSges(6,3) * t648 - pkin(4) * t644 + qJ(5) * t763 + (-qJ(5) * t696 + t769) * t712 + t771 * t706 + (-qJ(5) * mrSges(6,2) + t775) * t689 + t777 * t665 + t783 * t664;
t643 = -t665 * mrSges(6,3) - t706 * t679 + t764;
t770 = -t783 * t705 - t777 * t706 - t775 * t712;
t630 = mrSges(5,2) * t656 + mrSges(6,2) * t651 - mrSges(5,3) * t649 - mrSges(6,3) * t646 - qJ(5) * t643 + t777 * t664 + t784 * t665 + t776 * t689 - t771 * t705 + t770 * t712;
t699 = Ifges(4,4) * t717 - Ifges(4,2) * t716 + Ifges(4,6) * t736;
t700 = Ifges(4,1) * t717 - Ifges(4,4) * t716 + Ifges(4,5) * t736;
t751 = -mrSges(4,1) * t659 + mrSges(4,2) * t660 - Ifges(4,5) * t691 - Ifges(4,6) * t690 - Ifges(4,3) * t735 - pkin(3) * t750 - pkin(8) * t759 - t743 * t622 - t739 * t630 - t717 * t699 - t716 * t700;
t781 = mrSges(3,1) * t707 - mrSges(3,2) * t708 + Ifges(3,5) * t724 + Ifges(3,6) * t725 + Ifges(3,3) * qJDD(2) + pkin(2) * t620 + (t741 * t714 - t745 * t715) * qJD(1) - t751;
t780 = mrSges(5,1) * t649 + mrSges(6,1) * t646 - mrSges(5,2) * t650 - mrSges(6,2) * t648 + pkin(4) * t643 + t775 * t664 + t776 * t665 + t782 * t689 - t769 * t705 - t770 * t706;
t723 = (-mrSges(3,1) * t745 + mrSges(3,2) * t741) * qJD(1);
t766 = qJD(1) * t745;
t728 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t766;
t618 = m(3) * t707 + qJDD(2) * mrSges(3,1) - t724 * mrSges(3,3) + qJD(2) * t728 - t723 * t767 + t620;
t727 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t767;
t760 = t744 * t628 - t740 * t640;
t619 = m(3) * t708 - qJDD(2) * mrSges(3,2) + t725 * mrSges(3,3) - qJD(2) * t727 + t723 * t766 + t760;
t761 = -t741 * t618 + t745 * t619;
t611 = m(2) * t731 - t747 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t761;
t718 = -t747 * pkin(6) + t755;
t632 = t743 * t635 + t739 * t638;
t753 = m(4) * t692 - t690 * mrSges(4,1) + t691 * mrSges(4,2) + t716 * t709 + t717 * t710 + t632;
t749 = -m(3) * t718 + t725 * mrSges(3,1) - t724 * mrSges(3,2) - t727 * t767 + t728 * t766 - t753;
t624 = m(2) * t730 + qJDD(1) * mrSges(2,1) - t747 * mrSges(2,2) + t749;
t772 = t742 * t611 + t746 * t624;
t613 = t745 * t618 + t741 * t619;
t762 = t746 * t611 - t742 * t624;
t698 = Ifges(4,5) * t717 - Ifges(4,6) * t716 + Ifges(4,3) * t736;
t608 = mrSges(4,2) * t692 - mrSges(4,3) * t659 + Ifges(4,1) * t691 + Ifges(4,4) * t690 + Ifges(4,5) * t735 - pkin(8) * t632 - t739 * t622 + t743 * t630 - t716 * t698 - t736 * t699;
t614 = -mrSges(4,1) * t692 + mrSges(4,3) * t660 + Ifges(4,4) * t691 + Ifges(4,2) * t690 + Ifges(4,6) * t735 - pkin(3) * t632 - t717 * t698 + t736 * t700 - t780;
t713 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t741 + Ifges(3,6) * t745) * qJD(1);
t605 = -mrSges(3,1) * t718 + mrSges(3,3) * t708 + Ifges(3,4) * t724 + Ifges(3,2) * t725 + Ifges(3,6) * qJDD(2) - pkin(2) * t753 + pkin(7) * t760 + qJD(2) * t715 + t740 * t608 + t744 * t614 - t713 * t767;
t607 = mrSges(3,2) * t718 - mrSges(3,3) * t707 + Ifges(3,1) * t724 + Ifges(3,4) * t725 + Ifges(3,5) * qJDD(2) - pkin(7) * t620 - qJD(2) * t714 + t744 * t608 - t740 * t614 + t713 * t766;
t754 = mrSges(2,1) * t730 - mrSges(2,2) * t731 + Ifges(2,3) * qJDD(1) + pkin(1) * t749 + pkin(6) * t761 + t745 * t605 + t741 * t607;
t603 = mrSges(2,1) * g(3) + mrSges(2,3) * t731 + t747 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t613 - t781;
t602 = -mrSges(2,2) * g(3) - mrSges(2,3) * t730 + Ifges(2,5) * qJDD(1) - t747 * Ifges(2,6) - pkin(6) * t613 - t741 * t605 + t745 * t607;
t1 = [-m(1) * g(1) + t762; -m(1) * g(2) + t772; (-m(1) - m(2)) * g(3) + t613; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t772 + t746 * t602 - t742 * t603; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t762 + t742 * t602 + t746 * t603; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t754; t754; t781; -t751; t780; t644;];
tauJB = t1;
