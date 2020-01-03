% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:44
% EndTime: 2019-12-31 21:16:52
% DurationCPUTime: 7.99s
% Computational Cost: add. (122513->314), mult. (255186->399), div. (0->0), fcn. (177314->10), ass. (0->125)
t745 = sin(qJ(2));
t748 = cos(qJ(2));
t766 = qJD(1) * qJD(2);
t725 = t745 * qJDD(1) + t748 * t766;
t746 = sin(qJ(1));
t749 = cos(qJ(1));
t732 = -t749 * g(1) - t746 * g(2);
t750 = qJD(1) ^ 2;
t719 = -t750 * pkin(1) + qJDD(1) * pkin(6) + t732;
t770 = t745 * t719;
t771 = pkin(2) * t750;
t684 = qJDD(2) * pkin(2) - t725 * pkin(7) - t770 + (pkin(7) * t766 + t745 * t771 - g(3)) * t748;
t708 = -t745 * g(3) + t748 * t719;
t726 = t748 * qJDD(1) - t745 * t766;
t768 = qJD(1) * t745;
t730 = qJD(2) * pkin(2) - pkin(7) * t768;
t740 = t748 ^ 2;
t685 = t726 * pkin(7) - qJD(2) * t730 - t740 * t771 + t708;
t744 = sin(qJ(3));
t772 = cos(qJ(3));
t667 = t744 * t684 + t772 * t685;
t717 = (t744 * t748 + t772 * t745) * qJD(1);
t691 = t717 * qJD(3) + t744 * t725 - t772 * t726;
t767 = qJD(1) * t748;
t716 = t744 * t768 - t772 * t767;
t701 = t716 * mrSges(4,1) + t717 * mrSges(4,2);
t738 = qJD(2) + qJD(3);
t710 = t738 * mrSges(4,1) - t717 * mrSges(4,3);
t737 = qJDD(2) + qJDD(3);
t692 = -t716 * qJD(3) + t772 * t725 + t744 * t726;
t731 = t746 * g(1) - t749 * g(2);
t758 = -qJDD(1) * pkin(1) - t731;
t693 = -t726 * pkin(2) + t730 * t768 + (-pkin(7) * t740 - pkin(6)) * t750 + t758;
t656 = (t716 * t738 - t692) * qJ(4) + (t717 * t738 + t691) * pkin(3) + t693;
t700 = t716 * pkin(3) - t717 * qJ(4);
t736 = t738 ^ 2;
t659 = -t736 * pkin(3) + t737 * qJ(4) - t716 * t700 + t667;
t741 = sin(pkin(9));
t742 = cos(pkin(9));
t706 = t742 * t717 + t741 * t738;
t648 = -0.2e1 * qJD(4) * t706 + t742 * t656 - t741 * t659;
t678 = t742 * t692 + t741 * t737;
t705 = -t741 * t717 + t742 * t738;
t646 = (t716 * t705 - t678) * pkin(8) + (t705 * t706 + t691) * pkin(4) + t648;
t649 = 0.2e1 * qJD(4) * t705 + t741 * t656 + t742 * t659;
t677 = -t741 * t692 + t742 * t737;
t695 = t716 * pkin(4) - t706 * pkin(8);
t704 = t705 ^ 2;
t647 = -t704 * pkin(4) + t677 * pkin(8) - t716 * t695 + t649;
t743 = sin(qJ(5));
t747 = cos(qJ(5));
t644 = t747 * t646 - t743 * t647;
t675 = t747 * t705 - t743 * t706;
t655 = t675 * qJD(5) + t743 * t677 + t747 * t678;
t676 = t743 * t705 + t747 * t706;
t664 = -t675 * mrSges(6,1) + t676 * mrSges(6,2);
t712 = qJD(5) + t716;
t668 = -t712 * mrSges(6,2) + t675 * mrSges(6,3);
t690 = qJDD(5) + t691;
t640 = m(6) * t644 + t690 * mrSges(6,1) - t655 * mrSges(6,3) - t676 * t664 + t712 * t668;
t645 = t743 * t646 + t747 * t647;
t654 = -t676 * qJD(5) + t747 * t677 - t743 * t678;
t669 = t712 * mrSges(6,1) - t676 * mrSges(6,3);
t641 = m(6) * t645 - t690 * mrSges(6,2) + t654 * mrSges(6,3) + t675 * t664 - t712 * t669;
t632 = t747 * t640 + t743 * t641;
t680 = -t705 * mrSges(5,1) + t706 * mrSges(5,2);
t760 = -t716 * mrSges(5,2) + t705 * mrSges(5,3);
t630 = m(5) * t648 + t691 * mrSges(5,1) - t678 * mrSges(5,3) - t706 * t680 + t716 * t760 + t632;
t694 = t716 * mrSges(5,1) - t706 * mrSges(5,3);
t761 = -t743 * t640 + t747 * t641;
t631 = m(5) * t649 - t691 * mrSges(5,2) + t677 * mrSges(5,3) + t705 * t680 - t716 * t694 + t761;
t762 = -t741 * t630 + t742 * t631;
t623 = m(4) * t667 - t737 * mrSges(4,2) - t691 * mrSges(4,3) - t716 * t701 - t738 * t710 + t762;
t666 = t772 * t684 - t744 * t685;
t658 = -t737 * pkin(3) - t736 * qJ(4) + t717 * t700 + qJDD(4) - t666;
t650 = -t677 * pkin(4) - t704 * pkin(8) + t706 * t695 + t658;
t756 = m(6) * t650 - t654 * mrSges(6,1) + t655 * mrSges(6,2) - t675 * t668 + t676 * t669;
t643 = m(5) * t658 - t677 * mrSges(5,1) + t678 * mrSges(5,2) + t706 * t694 - t705 * t760 + t756;
t709 = -t738 * mrSges(4,2) - t716 * mrSges(4,3);
t636 = m(4) * t666 + t737 * mrSges(4,1) - t692 * mrSges(4,3) - t717 * t701 + t738 * t709 - t643;
t613 = t744 * t623 + t772 * t636;
t707 = -t748 * g(3) - t770;
t714 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t745 + Ifges(3,2) * t748) * qJD(1);
t715 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t745 + Ifges(3,4) * t748) * qJD(1);
t660 = Ifges(6,5) * t676 + Ifges(6,6) * t675 + Ifges(6,3) * t712;
t662 = Ifges(6,1) * t676 + Ifges(6,4) * t675 + Ifges(6,5) * t712;
t633 = -mrSges(6,1) * t650 + mrSges(6,3) * t645 + Ifges(6,4) * t655 + Ifges(6,2) * t654 + Ifges(6,6) * t690 - t676 * t660 + t712 * t662;
t661 = Ifges(6,4) * t676 + Ifges(6,2) * t675 + Ifges(6,6) * t712;
t634 = mrSges(6,2) * t650 - mrSges(6,3) * t644 + Ifges(6,1) * t655 + Ifges(6,4) * t654 + Ifges(6,5) * t690 + t675 * t660 - t712 * t661;
t670 = Ifges(5,5) * t706 + Ifges(5,6) * t705 + Ifges(5,3) * t716;
t672 = Ifges(5,1) * t706 + Ifges(5,4) * t705 + Ifges(5,5) * t716;
t615 = -mrSges(5,1) * t658 + mrSges(5,3) * t649 + Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t691 - pkin(4) * t756 + pkin(8) * t761 + t747 * t633 + t743 * t634 - t706 * t670 + t716 * t672;
t671 = Ifges(5,4) * t706 + Ifges(5,2) * t705 + Ifges(5,6) * t716;
t617 = mrSges(5,2) * t658 - mrSges(5,3) * t648 + Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t691 - pkin(8) * t632 - t743 * t633 + t747 * t634 + t705 * t670 - t716 * t671;
t697 = Ifges(4,4) * t717 - Ifges(4,2) * t716 + Ifges(4,6) * t738;
t698 = Ifges(4,1) * t717 - Ifges(4,4) * t716 + Ifges(4,5) * t738;
t754 = -mrSges(4,1) * t666 + mrSges(4,2) * t667 - Ifges(4,5) * t692 + Ifges(4,6) * t691 - Ifges(4,3) * t737 + pkin(3) * t643 - qJ(4) * t762 - t742 * t615 - t741 * t617 - t717 * t697 - t716 * t698;
t773 = mrSges(3,1) * t707 - mrSges(3,2) * t708 + Ifges(3,5) * t725 + Ifges(3,6) * t726 + Ifges(3,3) * qJDD(2) + pkin(2) * t613 + (t745 * t714 - t748 * t715) * qJD(1) - t754;
t724 = (-mrSges(3,1) * t748 + mrSges(3,2) * t745) * qJD(1);
t729 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t767;
t611 = m(3) * t707 + qJDD(2) * mrSges(3,1) - t725 * mrSges(3,3) + qJD(2) * t729 - t724 * t768 + t613;
t728 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t768;
t763 = t772 * t623 - t744 * t636;
t612 = m(3) * t708 - qJDD(2) * mrSges(3,2) + t726 * mrSges(3,3) - qJD(2) * t728 + t724 * t767 + t763;
t764 = -t745 * t611 + t748 * t612;
t604 = m(2) * t732 - t750 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t764;
t718 = -t750 * pkin(6) + t758;
t625 = t742 * t630 + t741 * t631;
t755 = m(4) * t693 + t691 * mrSges(4,1) + t692 * mrSges(4,2) + t716 * t709 + t717 * t710 + t625;
t752 = -m(3) * t718 + t726 * mrSges(3,1) - t725 * mrSges(3,2) - t728 * t768 + t729 * t767 - t755;
t619 = m(2) * t731 + qJDD(1) * mrSges(2,1) - t750 * mrSges(2,2) + t752;
t769 = t746 * t604 + t749 * t619;
t606 = t748 * t611 + t745 * t612;
t765 = t749 * t604 - t746 * t619;
t696 = Ifges(4,5) * t717 - Ifges(4,6) * t716 + Ifges(4,3) * t738;
t601 = mrSges(4,2) * t693 - mrSges(4,3) * t666 + Ifges(4,1) * t692 - Ifges(4,4) * t691 + Ifges(4,5) * t737 - qJ(4) * t625 - t741 * t615 + t742 * t617 - t716 * t696 - t738 * t697;
t753 = mrSges(6,1) * t644 - mrSges(6,2) * t645 + Ifges(6,5) * t655 + Ifges(6,6) * t654 + Ifges(6,3) * t690 + t676 * t661 - t675 * t662;
t607 = -t753 + (-Ifges(4,2) - Ifges(5,3)) * t691 + Ifges(4,6) * t737 + t738 * t698 - t717 * t696 + t705 * t672 - t706 * t671 - mrSges(4,1) * t693 + Ifges(4,4) * t692 - Ifges(5,6) * t677 - Ifges(5,5) * t678 + mrSges(4,3) * t667 - mrSges(5,1) * t648 + mrSges(5,2) * t649 - pkin(4) * t632 - pkin(3) * t625;
t713 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t745 + Ifges(3,6) * t748) * qJD(1);
t597 = -mrSges(3,1) * t718 + mrSges(3,3) * t708 + Ifges(3,4) * t725 + Ifges(3,2) * t726 + Ifges(3,6) * qJDD(2) - pkin(2) * t755 + pkin(7) * t763 + qJD(2) * t715 + t744 * t601 + t772 * t607 - t713 * t768;
t600 = mrSges(3,2) * t718 - mrSges(3,3) * t707 + Ifges(3,1) * t725 + Ifges(3,4) * t726 + Ifges(3,5) * qJDD(2) - pkin(7) * t613 - qJD(2) * t714 + t772 * t601 - t744 * t607 + t713 * t767;
t757 = mrSges(2,1) * t731 - mrSges(2,2) * t732 + Ifges(2,3) * qJDD(1) + pkin(1) * t752 + pkin(6) * t764 + t748 * t597 + t745 * t600;
t598 = mrSges(2,1) * g(3) + mrSges(2,3) * t732 + t750 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t606 - t773;
t595 = -mrSges(2,2) * g(3) - mrSges(2,3) * t731 + Ifges(2,5) * qJDD(1) - t750 * Ifges(2,6) - pkin(6) * t606 - t745 * t597 + t748 * t600;
t1 = [-m(1) * g(1) + t765; -m(1) * g(2) + t769; (-m(1) - m(2)) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t769 + t749 * t595 - t746 * t598; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t765 + t746 * t595 + t749 * t598; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t757; t757; t773; -t754; t643; t753;];
tauJB = t1;
