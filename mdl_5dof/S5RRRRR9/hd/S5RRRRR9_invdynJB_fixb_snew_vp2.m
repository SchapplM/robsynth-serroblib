% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:28:37
% EndTime: 2019-12-31 22:28:47
% DurationCPUTime: 9.69s
% Computational Cost: add. (152980->314), mult. (310659->393), div. (0->0), fcn. (218081->10), ass. (0->127)
t747 = sin(qJ(1));
t752 = cos(qJ(1));
t735 = t747 * g(1) - t752 * g(2);
t754 = qJD(1) ^ 2;
t718 = -qJDD(1) * pkin(1) - t754 * pkin(6) - t735;
t746 = sin(qJ(2));
t751 = cos(qJ(2));
t769 = qJD(1) * qJD(2);
t768 = t751 * t769;
t729 = t746 * qJDD(1) + t768;
t739 = t746 * t769;
t730 = t751 * qJDD(1) - t739;
t685 = (-t729 - t768) * pkin(7) + (-t730 + t739) * pkin(2) + t718;
t736 = -t752 * g(1) - t747 * g(2);
t719 = -t754 * pkin(1) + qJDD(1) * pkin(6) + t736;
t707 = -t746 * g(3) + t751 * t719;
t728 = (-pkin(2) * t751 - pkin(7) * t746) * qJD(1);
t753 = qJD(2) ^ 2;
t770 = t751 * qJD(1);
t688 = -t753 * pkin(2) + qJDD(2) * pkin(7) + t728 * t770 + t707;
t745 = sin(qJ(3));
t750 = cos(qJ(3));
t672 = t750 * t685 - t745 * t688;
t771 = qJD(1) * t746;
t725 = t750 * qJD(2) - t745 * t771;
t699 = t725 * qJD(3) + t745 * qJDD(2) + t750 * t729;
t724 = qJDD(3) - t730;
t726 = t745 * qJD(2) + t750 * t771;
t738 = qJD(3) - t770;
t656 = (t725 * t738 - t699) * pkin(8) + (t725 * t726 + t724) * pkin(3) + t672;
t673 = t745 * t685 + t750 * t688;
t698 = -t726 * qJD(3) + t750 * qJDD(2) - t745 * t729;
t708 = t738 * pkin(3) - t726 * pkin(8);
t723 = t725 ^ 2;
t658 = -t723 * pkin(3) + t698 * pkin(8) - t738 * t708 + t673;
t744 = sin(qJ(4));
t749 = cos(qJ(4));
t643 = t749 * t656 - t744 * t658;
t701 = t749 * t725 - t744 * t726;
t671 = t701 * qJD(4) + t744 * t698 + t749 * t699;
t702 = t744 * t725 + t749 * t726;
t720 = qJDD(4) + t724;
t737 = qJD(4) + t738;
t640 = (t701 * t737 - t671) * pkin(9) + (t701 * t702 + t720) * pkin(4) + t643;
t644 = t744 * t656 + t749 * t658;
t670 = -t702 * qJD(4) + t749 * t698 - t744 * t699;
t691 = t737 * pkin(4) - t702 * pkin(9);
t700 = t701 ^ 2;
t641 = -t700 * pkin(4) + t670 * pkin(9) - t737 * t691 + t644;
t743 = sin(qJ(5));
t748 = cos(qJ(5));
t639 = t743 * t640 + t748 * t641;
t706 = -t751 * g(3) - t746 * t719;
t687 = -qJDD(2) * pkin(2) - t753 * pkin(7) + t728 * t771 - t706;
t665 = -t698 * pkin(3) - t723 * pkin(8) + t726 * t708 + t687;
t646 = -t670 * pkin(4) - t700 * pkin(9) + t702 * t691 + t665;
t681 = t743 * t701 + t748 * t702;
t651 = -t681 * qJD(5) + t748 * t670 - t743 * t671;
t680 = t748 * t701 - t743 * t702;
t652 = t680 * qJD(5) + t743 * t670 + t748 * t671;
t732 = qJD(5) + t737;
t659 = Ifges(6,5) * t681 + Ifges(6,6) * t680 + Ifges(6,3) * t732;
t661 = Ifges(6,1) * t681 + Ifges(6,4) * t680 + Ifges(6,5) * t732;
t714 = qJDD(5) + t720;
t627 = -mrSges(6,1) * t646 + mrSges(6,3) * t639 + Ifges(6,4) * t652 + Ifges(6,2) * t651 + Ifges(6,6) * t714 - t681 * t659 + t732 * t661;
t638 = t748 * t640 - t743 * t641;
t660 = Ifges(6,4) * t681 + Ifges(6,2) * t680 + Ifges(6,6) * t732;
t628 = mrSges(6,2) * t646 - mrSges(6,3) * t638 + Ifges(6,1) * t652 + Ifges(6,4) * t651 + Ifges(6,5) * t714 + t680 * t659 - t732 * t660;
t676 = Ifges(5,5) * t702 + Ifges(5,6) * t701 + Ifges(5,3) * t737;
t678 = Ifges(5,1) * t702 + Ifges(5,4) * t701 + Ifges(5,5) * t737;
t674 = -t732 * mrSges(6,2) + t680 * mrSges(6,3);
t675 = t732 * mrSges(6,1) - t681 * mrSges(6,3);
t763 = m(6) * t646 - t651 * mrSges(6,1) + t652 * mrSges(6,2) - t680 * t674 + t681 * t675;
t664 = -mrSges(6,1) * t680 + mrSges(6,2) * t681;
t635 = m(6) * t638 + t714 * mrSges(6,1) - t652 * mrSges(6,3) - t681 * t664 + t732 * t674;
t636 = m(6) * t639 - t714 * mrSges(6,2) + t651 * mrSges(6,3) + t680 * t664 - t732 * t675;
t764 = -t743 * t635 + t748 * t636;
t614 = -mrSges(5,1) * t665 + mrSges(5,3) * t644 + Ifges(5,4) * t671 + Ifges(5,2) * t670 + Ifges(5,6) * t720 - pkin(4) * t763 + pkin(9) * t764 + t748 * t627 + t743 * t628 - t702 * t676 + t737 * t678;
t626 = t748 * t635 + t743 * t636;
t677 = Ifges(5,4) * t702 + Ifges(5,2) * t701 + Ifges(5,6) * t737;
t615 = mrSges(5,2) * t665 - mrSges(5,3) * t643 + Ifges(5,1) * t671 + Ifges(5,4) * t670 + Ifges(5,5) * t720 - pkin(9) * t626 - t743 * t627 + t748 * t628 + t701 * t676 - t737 * t677;
t692 = Ifges(4,5) * t726 + Ifges(4,6) * t725 + Ifges(4,3) * t738;
t694 = Ifges(4,1) * t726 + Ifges(4,4) * t725 + Ifges(4,5) * t738;
t689 = -t737 * mrSges(5,2) + t701 * mrSges(5,3);
t690 = t737 * mrSges(5,1) - t702 * mrSges(5,3);
t759 = m(5) * t665 - t670 * mrSges(5,1) + t671 * mrSges(5,2) - t701 * t689 + t702 * t690 + t763;
t682 = -mrSges(5,1) * t701 + mrSges(5,2) * t702;
t623 = m(5) * t643 + t720 * mrSges(5,1) - t671 * mrSges(5,3) - t702 * t682 + t737 * t689 + t626;
t624 = m(5) * t644 - t720 * mrSges(5,2) + t670 * mrSges(5,3) + t701 * t682 - t737 * t690 + t764;
t765 = -t744 * t623 + t749 * t624;
t599 = -mrSges(4,1) * t687 + mrSges(4,3) * t673 + Ifges(4,4) * t699 + Ifges(4,2) * t698 + Ifges(4,6) * t724 - pkin(3) * t759 + pkin(8) * t765 + t749 * t614 + t744 * t615 - t726 * t692 + t738 * t694;
t619 = t749 * t623 + t744 * t624;
t693 = Ifges(4,4) * t726 + Ifges(4,2) * t725 + Ifges(4,6) * t738;
t600 = mrSges(4,2) * t687 - mrSges(4,3) * t672 + Ifges(4,1) * t699 + Ifges(4,4) * t698 + Ifges(4,5) * t724 - pkin(8) * t619 - t744 * t614 + t749 * t615 + t725 * t692 - t738 * t693;
t703 = -t725 * mrSges(4,1) + t726 * mrSges(4,2);
t704 = -t738 * mrSges(4,2) + t725 * mrSges(4,3);
t617 = m(4) * t672 + t724 * mrSges(4,1) - t699 * mrSges(4,3) - t726 * t703 + t738 * t704 + t619;
t705 = t738 * mrSges(4,1) - t726 * mrSges(4,3);
t618 = m(4) * t673 - t724 * mrSges(4,2) + t698 * mrSges(4,3) + t725 * t703 - t738 * t705 + t765;
t613 = -t745 * t617 + t750 * t618;
t631 = -m(4) * t687 + t698 * mrSges(4,1) - t699 * mrSges(4,2) + t725 * t704 - t726 * t705 - t759;
t716 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t746 + Ifges(3,2) * t751) * qJD(1);
t717 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t746 + Ifges(3,4) * t751) * qJD(1);
t773 = mrSges(3,1) * t706 - mrSges(3,2) * t707 + Ifges(3,5) * t729 + Ifges(3,6) * t730 + Ifges(3,3) * qJDD(2) + pkin(2) * t631 + pkin(7) * t613 + t750 * t599 + t745 * t600 + (t746 * t716 - t751 * t717) * qJD(1);
t727 = (-mrSges(3,1) * t751 + mrSges(3,2) * t746) * qJD(1);
t733 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t771;
t611 = m(3) * t707 - qJDD(2) * mrSges(3,2) + t730 * mrSges(3,3) - qJD(2) * t733 + t727 * t770 + t613;
t734 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t770;
t630 = m(3) * t706 + qJDD(2) * mrSges(3,1) - t729 * mrSges(3,3) + qJD(2) * t734 - t727 * t771 + t631;
t766 = t751 * t611 - t746 * t630;
t603 = m(2) * t736 - t754 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t766;
t612 = t750 * t617 + t745 * t618;
t758 = -m(3) * t718 + t730 * mrSges(3,1) - t729 * mrSges(3,2) - t733 * t771 + t734 * t770 - t612;
t607 = m(2) * t735 + qJDD(1) * mrSges(2,1) - t754 * mrSges(2,2) + t758;
t772 = t747 * t603 + t752 * t607;
t605 = t746 * t611 + t751 * t630;
t767 = t752 * t603 - t747 * t607;
t715 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t746 + Ifges(3,6) * t751) * qJD(1);
t596 = mrSges(3,2) * t718 - mrSges(3,3) * t706 + Ifges(3,1) * t729 + Ifges(3,4) * t730 + Ifges(3,5) * qJDD(2) - pkin(7) * t612 - qJD(2) * t716 - t745 * t599 + t750 * t600 + t715 * t770;
t760 = -mrSges(6,1) * t638 + mrSges(6,2) * t639 - Ifges(6,5) * t652 - Ifges(6,6) * t651 - Ifges(6,3) * t714 - t681 * t660 + t680 * t661;
t757 = -mrSges(5,1) * t643 + mrSges(5,2) * t644 - Ifges(5,5) * t671 - Ifges(5,6) * t670 - Ifges(5,3) * t720 - pkin(4) * t626 - t702 * t677 + t701 * t678 + t760;
t755 = mrSges(4,1) * t672 - mrSges(4,2) * t673 + Ifges(4,5) * t699 + Ifges(4,6) * t698 + Ifges(4,3) * t724 + pkin(3) * t619 + t726 * t693 - t725 * t694 - t757;
t598 = -mrSges(3,1) * t718 + mrSges(3,3) * t707 + Ifges(3,4) * t729 + Ifges(3,2) * t730 + Ifges(3,6) * qJDD(2) - pkin(2) * t612 + qJD(2) * t717 - t715 * t771 - t755;
t761 = mrSges(2,1) * t735 - mrSges(2,2) * t736 + Ifges(2,3) * qJDD(1) + pkin(1) * t758 + pkin(6) * t766 + t746 * t596 + t751 * t598;
t594 = mrSges(2,1) * g(3) + mrSges(2,3) * t736 + t754 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t605 - t773;
t593 = -mrSges(2,2) * g(3) - mrSges(2,3) * t735 + Ifges(2,5) * qJDD(1) - t754 * Ifges(2,6) - pkin(6) * t605 + t751 * t596 - t746 * t598;
t1 = [-m(1) * g(1) + t767; -m(1) * g(2) + t772; (-m(1) - m(2)) * g(3) + t605; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t772 + t752 * t593 - t747 * t594; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t767 + t747 * t593 + t752 * t594; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t761; t761; t773; t755; -t757; -t760;];
tauJB = t1;
