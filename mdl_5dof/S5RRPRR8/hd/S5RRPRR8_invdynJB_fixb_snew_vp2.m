% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:17:28
% EndTime: 2019-12-31 20:17:37
% DurationCPUTime: 8.12s
% Computational Cost: add. (113426->314), mult. (263005->400), div. (0->0), fcn. (186760->10), ass. (0->125)
t744 = sin(qJ(2));
t748 = cos(qJ(2));
t765 = qJD(1) * qJD(2);
t725 = t744 * qJDD(1) + t748 * t765;
t745 = sin(qJ(1));
t749 = cos(qJ(1));
t732 = -t749 * g(1) - t745 * g(2);
t750 = qJD(1) ^ 2;
t720 = -t750 * pkin(1) + qJDD(1) * pkin(6) + t732;
t769 = t744 * t720;
t771 = pkin(2) * t750;
t686 = qJDD(2) * pkin(2) - t725 * qJ(3) - t769 + (qJ(3) * t765 + t744 * t771 - g(3)) * t748;
t706 = -t744 * g(3) + t748 * t720;
t726 = t748 * qJDD(1) - t744 * t765;
t767 = qJD(1) * t744;
t728 = qJD(2) * pkin(2) - qJ(3) * t767;
t739 = t748 ^ 2;
t687 = t726 * qJ(3) - qJD(2) * t728 - t739 * t771 + t706;
t740 = sin(pkin(9));
t741 = cos(pkin(9));
t715 = (t740 * t748 + t741 * t744) * qJD(1);
t663 = -0.2e1 * qJD(3) * t715 + t741 * t686 - t740 * t687;
t704 = t741 * t725 + t740 * t726;
t714 = (-t740 * t744 + t741 * t748) * qJD(1);
t651 = (qJD(2) * t714 - t704) * pkin(7) + (t714 * t715 + qJDD(2)) * pkin(3) + t663;
t664 = 0.2e1 * qJD(3) * t714 + t740 * t686 + t741 * t687;
t703 = -t740 * t725 + t741 * t726;
t709 = qJD(2) * pkin(3) - t715 * pkin(7);
t713 = t714 ^ 2;
t653 = -t713 * pkin(3) + t703 * pkin(7) - qJD(2) * t709 + t664;
t743 = sin(qJ(4));
t747 = cos(qJ(4));
t649 = t743 * t651 + t747 * t653;
t698 = t743 * t714 + t747 * t715;
t671 = -t698 * qJD(4) + t747 * t703 - t743 * t704;
t697 = t747 * t714 - t743 * t715;
t681 = -t697 * mrSges(5,1) + t698 * mrSges(5,2);
t737 = qJD(2) + qJD(4);
t692 = t737 * mrSges(5,1) - t698 * mrSges(5,3);
t736 = qJDD(2) + qJDD(4);
t682 = -t697 * pkin(4) - t698 * pkin(8);
t735 = t737 ^ 2;
t645 = -t735 * pkin(4) + t736 * pkin(8) + t697 * t682 + t649;
t731 = t745 * g(1) - t749 * g(2);
t758 = -qJDD(1) * pkin(1) - t731;
t688 = -t726 * pkin(2) + qJDD(3) + t728 * t767 + (-qJ(3) * t739 - pkin(6)) * t750 + t758;
t662 = -t703 * pkin(3) - t713 * pkin(7) + t715 * t709 + t688;
t672 = t697 * qJD(4) + t743 * t703 + t747 * t704;
t646 = (-t697 * t737 - t672) * pkin(8) + (t698 * t737 - t671) * pkin(4) + t662;
t742 = sin(qJ(5));
t746 = cos(qJ(5));
t642 = -t742 * t645 + t746 * t646;
t689 = -t742 * t698 + t746 * t737;
t656 = t689 * qJD(5) + t746 * t672 + t742 * t736;
t670 = qJDD(5) - t671;
t690 = t746 * t698 + t742 * t737;
t673 = -t689 * mrSges(6,1) + t690 * mrSges(6,2);
t693 = qJD(5) - t697;
t674 = -t693 * mrSges(6,2) + t689 * mrSges(6,3);
t638 = m(6) * t642 + t670 * mrSges(6,1) - t656 * mrSges(6,3) - t690 * t673 + t693 * t674;
t643 = t746 * t645 + t742 * t646;
t655 = -t690 * qJD(5) - t742 * t672 + t746 * t736;
t675 = t693 * mrSges(6,1) - t690 * mrSges(6,3);
t639 = m(6) * t643 - t670 * mrSges(6,2) + t655 * mrSges(6,3) + t689 * t673 - t693 * t675;
t760 = -t742 * t638 + t746 * t639;
t624 = m(5) * t649 - t736 * mrSges(5,2) + t671 * mrSges(5,3) + t697 * t681 - t737 * t692 + t760;
t648 = t747 * t651 - t743 * t653;
t691 = -t737 * mrSges(5,2) + t697 * mrSges(5,3);
t644 = -t736 * pkin(4) - t735 * pkin(8) + t698 * t682 - t648;
t755 = -m(6) * t644 + t655 * mrSges(6,1) - t656 * mrSges(6,2) + t689 * t674 - t690 * t675;
t634 = m(5) * t648 + t736 * mrSges(5,1) - t672 * mrSges(5,3) - t698 * t681 + t737 * t691 + t755;
t618 = t743 * t624 + t747 * t634;
t701 = -t714 * mrSges(4,1) + t715 * mrSges(4,2);
t707 = -qJD(2) * mrSges(4,2) + t714 * mrSges(4,3);
t616 = m(4) * t663 + qJDD(2) * mrSges(4,1) - t704 * mrSges(4,3) + qJD(2) * t707 - t715 * t701 + t618;
t708 = qJD(2) * mrSges(4,1) - t715 * mrSges(4,3);
t761 = t747 * t624 - t743 * t634;
t617 = m(4) * t664 - qJDD(2) * mrSges(4,2) + t703 * mrSges(4,3) - qJD(2) * t708 + t714 * t701 + t761;
t610 = t741 * t616 + t740 * t617;
t695 = Ifges(4,4) * t715 + Ifges(4,2) * t714 + Ifges(4,6) * qJD(2);
t696 = Ifges(4,1) * t715 + Ifges(4,4) * t714 + Ifges(4,5) * qJD(2);
t705 = -t748 * g(3) - t769;
t717 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t744 + Ifges(3,2) * t748) * qJD(1);
t718 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t744 + Ifges(3,4) * t748) * qJD(1);
t657 = Ifges(6,5) * t690 + Ifges(6,6) * t689 + Ifges(6,3) * t693;
t659 = Ifges(6,1) * t690 + Ifges(6,4) * t689 + Ifges(6,5) * t693;
t631 = -mrSges(6,1) * t644 + mrSges(6,3) * t643 + Ifges(6,4) * t656 + Ifges(6,2) * t655 + Ifges(6,6) * t670 - t690 * t657 + t693 * t659;
t658 = Ifges(6,4) * t690 + Ifges(6,2) * t689 + Ifges(6,6) * t693;
t632 = mrSges(6,2) * t644 - mrSges(6,3) * t642 + Ifges(6,1) * t656 + Ifges(6,4) * t655 + Ifges(6,5) * t670 + t689 * t657 - t693 * t658;
t677 = Ifges(5,4) * t698 + Ifges(5,2) * t697 + Ifges(5,6) * t737;
t678 = Ifges(5,1) * t698 + Ifges(5,4) * t697 + Ifges(5,5) * t737;
t754 = -mrSges(5,1) * t648 + mrSges(5,2) * t649 - Ifges(5,5) * t672 - Ifges(5,6) * t671 - Ifges(5,3) * t736 - pkin(4) * t755 - pkin(8) * t760 - t746 * t631 - t742 * t632 - t698 * t677 + t697 * t678;
t772 = mrSges(3,1) * t705 + mrSges(4,1) * t663 - mrSges(3,2) * t706 - mrSges(4,2) * t664 + Ifges(3,5) * t725 + Ifges(4,5) * t704 + Ifges(3,6) * t726 + Ifges(4,6) * t703 + pkin(2) * t610 + pkin(3) * t618 + (t744 * t717 - t748 * t718) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t715 * t695 - t714 * t696 - t754;
t724 = (-mrSges(3,1) * t748 + mrSges(3,2) * t744) * qJD(1);
t766 = qJD(1) * t748;
t730 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t766;
t608 = m(3) * t705 + qJDD(2) * mrSges(3,1) - t725 * mrSges(3,3) + qJD(2) * t730 - t724 * t767 + t610;
t729 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t767;
t762 = -t740 * t616 + t741 * t617;
t609 = m(3) * t706 - qJDD(2) * mrSges(3,2) + t726 * mrSges(3,3) - qJD(2) * t729 + t724 * t766 + t762;
t763 = -t744 * t608 + t748 * t609;
t601 = m(2) * t732 - t750 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t763;
t627 = t746 * t638 + t742 * t639;
t757 = m(5) * t662 - t671 * mrSges(5,1) + t672 * mrSges(5,2) - t697 * t691 + t698 * t692 + t627;
t625 = m(4) * t688 - t703 * mrSges(4,1) + t704 * mrSges(4,2) - t714 * t707 + t715 * t708 + t757;
t719 = -t750 * pkin(6) + t758;
t752 = -m(3) * t719 + t726 * mrSges(3,1) - t725 * mrSges(3,2) - t729 * t767 + t730 * t766 - t625;
t620 = m(2) * t731 + qJDD(1) * mrSges(2,1) - t750 * mrSges(2,2) + t752;
t768 = t745 * t601 + t749 * t620;
t603 = t748 * t608 + t744 * t609;
t764 = t749 * t601 - t745 * t620;
t676 = Ifges(5,5) * t698 + Ifges(5,6) * t697 + Ifges(5,3) * t737;
t611 = mrSges(5,2) * t662 - mrSges(5,3) * t648 + Ifges(5,1) * t672 + Ifges(5,4) * t671 + Ifges(5,5) * t736 - pkin(8) * t627 - t742 * t631 + t746 * t632 + t697 * t676 - t737 * t677;
t753 = mrSges(6,1) * t642 - mrSges(6,2) * t643 + Ifges(6,5) * t656 + Ifges(6,6) * t655 + Ifges(6,3) * t670 + t690 * t658 - t689 * t659;
t612 = -mrSges(5,1) * t662 + mrSges(5,3) * t649 + Ifges(5,4) * t672 + Ifges(5,2) * t671 + Ifges(5,6) * t736 - pkin(4) * t627 - t698 * t676 + t737 * t678 - t753;
t694 = Ifges(4,5) * t715 + Ifges(4,6) * t714 + Ifges(4,3) * qJD(2);
t598 = -mrSges(4,1) * t688 + mrSges(4,3) * t664 + Ifges(4,4) * t704 + Ifges(4,2) * t703 + Ifges(4,6) * qJDD(2) - pkin(3) * t757 + pkin(7) * t761 + qJD(2) * t696 + t743 * t611 + t747 * t612 - t715 * t694;
t604 = mrSges(4,2) * t688 - mrSges(4,3) * t663 + Ifges(4,1) * t704 + Ifges(4,4) * t703 + Ifges(4,5) * qJDD(2) - pkin(7) * t618 - qJD(2) * t695 + t747 * t611 - t743 * t612 + t714 * t694;
t716 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t744 + Ifges(3,6) * t748) * qJD(1);
t594 = -mrSges(3,1) * t719 + mrSges(3,3) * t706 + Ifges(3,4) * t725 + Ifges(3,2) * t726 + Ifges(3,6) * qJDD(2) - pkin(2) * t625 + qJ(3) * t762 + qJD(2) * t718 + t741 * t598 + t740 * t604 - t716 * t767;
t596 = mrSges(3,2) * t719 - mrSges(3,3) * t705 + Ifges(3,1) * t725 + Ifges(3,4) * t726 + Ifges(3,5) * qJDD(2) - qJ(3) * t610 - qJD(2) * t717 - t740 * t598 + t741 * t604 + t716 * t766;
t756 = mrSges(2,1) * t731 - mrSges(2,2) * t732 + Ifges(2,3) * qJDD(1) + pkin(1) * t752 + pkin(6) * t763 + t748 * t594 + t744 * t596;
t597 = mrSges(2,1) * g(3) + mrSges(2,3) * t732 + t750 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t603 - t772;
t592 = -mrSges(2,2) * g(3) - mrSges(2,3) * t731 + Ifges(2,5) * qJDD(1) - t750 * Ifges(2,6) - pkin(6) * t603 - t744 * t594 + t748 * t596;
t1 = [-m(1) * g(1) + t764; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t603; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t768 + t749 * t592 - t745 * t597; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t764 + t745 * t592 + t749 * t597; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t756; t756; t772; t625; -t754; t753;];
tauJB = t1;
