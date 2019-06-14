% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 15:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:14:51
% EndTime: 2019-05-06 15:15:49
% DurationCPUTime: 52.08s
% Computational Cost: add. (827012->397), mult. (1893377->517), div. (0->0), fcn. (1539848->14), ass. (0->160)
t785 = cos(qJ(4));
t754 = cos(pkin(6));
t784 = g(3) * t754;
t751 = sin(pkin(6));
t757 = sin(qJ(2));
t783 = t751 * t757;
t760 = cos(qJ(2));
t782 = t751 * t760;
t781 = t754 * t757;
t780 = t754 * t760;
t758 = sin(qJ(1));
t761 = cos(qJ(1));
t741 = t758 * g(1) - g(2) * t761;
t762 = qJD(1) ^ 2;
t733 = pkin(8) * t751 * t762 + qJDD(1) * pkin(1) + t741;
t742 = -g(1) * t761 - g(2) * t758;
t775 = qJDD(1) * t751;
t734 = -pkin(1) * t762 + pkin(8) * t775 + t742;
t778 = t733 * t781 + t760 * t734;
t704 = -g(3) * t783 + t778;
t746 = qJD(1) * t754 + qJD(2);
t777 = qJD(1) * t751;
t774 = t757 * t777;
t731 = mrSges(3,1) * t746 - mrSges(3,3) * t774;
t736 = (-mrSges(3,1) * t760 + mrSges(3,2) * t757) * t777;
t738 = -qJD(2) * t774 + t760 * t775;
t745 = qJDD(1) * t754 + qJDD(2);
t735 = (-pkin(2) * t760 - qJ(3) * t757) * t777;
t744 = t746 ^ 2;
t776 = qJD(1) * t760;
t689 = -pkin(2) * t744 + qJ(3) * t745 + (-g(3) * t757 + t735 * t776) * t751 + t778;
t737 = (qJD(2) * t776 + qJDD(1) * t757) * t751;
t690 = -pkin(2) * t738 - t784 - qJ(3) * t737 + (-t733 + (pkin(2) * t757 - qJ(3) * t760) * t746 * qJD(1)) * t751;
t750 = sin(pkin(11));
t753 = cos(pkin(11));
t727 = t746 * t750 + t753 * t774;
t656 = -0.2e1 * qJD(3) * t727 - t689 * t750 + t753 * t690;
t714 = t737 * t753 + t745 * t750;
t726 = t746 * t753 - t750 * t774;
t773 = t751 * t776;
t647 = (-t726 * t773 - t714) * pkin(9) + (t726 * t727 - t738) * pkin(3) + t656;
t657 = 0.2e1 * qJD(3) * t726 + t753 * t689 + t750 * t690;
t713 = -t737 * t750 + t745 * t753;
t715 = -pkin(3) * t773 - pkin(9) * t727;
t725 = t726 ^ 2;
t654 = -pkin(3) * t725 + pkin(9) * t713 + t715 * t773 + t657;
t756 = sin(qJ(4));
t639 = t756 * t647 + t654 * t785;
t707 = t756 * t726 + t727 * t785;
t674 = qJD(4) * t707 - t713 * t785 + t714 * t756;
t706 = -t726 * t785 + t727 * t756;
t684 = mrSges(5,1) * t706 + mrSges(5,2) * t707;
t740 = qJD(4) - t773;
t697 = mrSges(5,1) * t740 - mrSges(5,3) * t707;
t730 = qJDD(4) - t738;
t683 = pkin(4) * t706 - qJ(5) * t707;
t739 = t740 ^ 2;
t637 = -pkin(4) * t739 + qJ(5) * t730 - t683 * t706 + t639;
t703 = -g(3) * t782 + t733 * t780 - t757 * t734;
t688 = -pkin(2) * t745 - qJ(3) * t744 + t735 * t774 + qJDD(3) - t703;
t658 = -pkin(3) * t713 - pkin(9) * t725 + t727 * t715 + t688;
t675 = -t706 * qJD(4) + t756 * t713 + t714 * t785;
t642 = (t706 * t740 - t675) * qJ(5) + (t707 * t740 + t674) * pkin(4) + t658;
t749 = sin(pkin(12));
t752 = cos(pkin(12));
t695 = t707 * t752 + t740 * t749;
t632 = -0.2e1 * qJD(5) * t695 - t637 * t749 + t752 * t642;
t668 = t675 * t752 + t730 * t749;
t694 = -t707 * t749 + t740 * t752;
t630 = (t694 * t706 - t668) * pkin(10) + (t694 * t695 + t674) * pkin(5) + t632;
t633 = 0.2e1 * qJD(5) * t694 + t752 * t637 + t749 * t642;
t667 = -t675 * t749 + t730 * t752;
t678 = pkin(5) * t706 - pkin(10) * t695;
t693 = t694 ^ 2;
t631 = -pkin(5) * t693 + pkin(10) * t667 - t678 * t706 + t633;
t755 = sin(qJ(6));
t759 = cos(qJ(6));
t628 = t630 * t759 - t631 * t755;
t669 = t694 * t759 - t695 * t755;
t645 = qJD(6) * t669 + t667 * t755 + t668 * t759;
t670 = t694 * t755 + t695 * t759;
t655 = -mrSges(7,1) * t669 + mrSges(7,2) * t670;
t705 = qJD(6) + t706;
t659 = -mrSges(7,2) * t705 + mrSges(7,3) * t669;
t673 = qJDD(6) + t674;
t626 = m(7) * t628 + mrSges(7,1) * t673 - mrSges(7,3) * t645 - t655 * t670 + t659 * t705;
t629 = t630 * t755 + t631 * t759;
t644 = -qJD(6) * t670 + t667 * t759 - t668 * t755;
t660 = mrSges(7,1) * t705 - mrSges(7,3) * t670;
t627 = m(7) * t629 - mrSges(7,2) * t673 + mrSges(7,3) * t644 + t655 * t669 - t660 * t705;
t618 = t759 * t626 + t755 * t627;
t671 = -mrSges(6,1) * t694 + mrSges(6,2) * t695;
t676 = -mrSges(6,2) * t706 + mrSges(6,3) * t694;
t616 = m(6) * t632 + mrSges(6,1) * t674 - mrSges(6,3) * t668 - t671 * t695 + t676 * t706 + t618;
t677 = mrSges(6,1) * t706 - mrSges(6,3) * t695;
t768 = -t626 * t755 + t759 * t627;
t617 = m(6) * t633 - mrSges(6,2) * t674 + mrSges(6,3) * t667 + t671 * t694 - t677 * t706 + t768;
t769 = -t616 * t749 + t752 * t617;
t611 = m(5) * t639 - mrSges(5,2) * t730 - mrSges(5,3) * t674 - t684 * t706 - t697 * t740 + t769;
t638 = t647 * t785 - t756 * t654;
t696 = -mrSges(5,2) * t740 - mrSges(5,3) * t706;
t636 = -t730 * pkin(4) - t739 * qJ(5) + t707 * t683 + qJDD(5) - t638;
t634 = -t667 * pkin(5) - t693 * pkin(10) + t695 * t678 + t636;
t766 = m(7) * t634 - t644 * mrSges(7,1) + mrSges(7,2) * t645 - t669 * t659 + t660 * t670;
t764 = -m(6) * t636 + t667 * mrSges(6,1) - mrSges(6,2) * t668 + t694 * t676 - t677 * t695 - t766;
t622 = m(5) * t638 + mrSges(5,1) * t730 - mrSges(5,3) * t675 - t684 * t707 + t696 * t740 + t764;
t603 = t756 * t611 + t622 * t785;
t708 = -mrSges(4,1) * t726 + mrSges(4,2) * t727;
t711 = mrSges(4,2) * t773 + mrSges(4,3) * t726;
t601 = m(4) * t656 - mrSges(4,1) * t738 - mrSges(4,3) * t714 - t708 * t727 - t711 * t773 + t603;
t712 = -mrSges(4,1) * t773 - mrSges(4,3) * t727;
t770 = t611 * t785 - t622 * t756;
t602 = m(4) * t657 + mrSges(4,2) * t738 + mrSges(4,3) * t713 + t708 * t726 + t712 * t773 + t770;
t771 = -t601 * t750 + t753 * t602;
t593 = m(3) * t704 - mrSges(3,2) * t745 + mrSges(3,3) * t738 - t731 * t746 + t736 * t773 + t771;
t596 = t753 * t601 + t750 * t602;
t719 = -t733 * t751 - t784;
t732 = -mrSges(3,2) * t746 + mrSges(3,3) * t773;
t595 = m(3) * t719 - mrSges(3,1) * t738 + mrSges(3,2) * t737 + (t731 * t757 - t732 * t760) * t777 + t596;
t612 = t752 * t616 + t749 * t617;
t765 = m(5) * t658 + t674 * mrSges(5,1) + mrSges(5,2) * t675 + t706 * t696 + t697 * t707 + t612;
t763 = -m(4) * t688 + t713 * mrSges(4,1) - mrSges(4,2) * t714 + t726 * t711 - t712 * t727 - t765;
t608 = m(3) * t703 + mrSges(3,1) * t745 - mrSges(3,3) * t737 + t732 * t746 - t736 * t774 + t763;
t584 = t593 * t781 - t595 * t751 + t608 * t780;
t582 = m(2) * t741 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t762 + t584;
t588 = t760 * t593 - t608 * t757;
t587 = m(2) * t742 - mrSges(2,1) * t762 - qJDD(1) * mrSges(2,2) + t588;
t779 = t761 * t582 + t758 * t587;
t583 = t593 * t783 + t754 * t595 + t608 * t782;
t772 = -t582 * t758 + t761 * t587;
t648 = Ifges(7,5) * t670 + Ifges(7,6) * t669 + Ifges(7,3) * t705;
t650 = Ifges(7,1) * t670 + Ifges(7,4) * t669 + Ifges(7,5) * t705;
t619 = -mrSges(7,1) * t634 + mrSges(7,3) * t629 + Ifges(7,4) * t645 + Ifges(7,2) * t644 + Ifges(7,6) * t673 - t648 * t670 + t650 * t705;
t649 = Ifges(7,4) * t670 + Ifges(7,2) * t669 + Ifges(7,6) * t705;
t620 = mrSges(7,2) * t634 - mrSges(7,3) * t628 + Ifges(7,1) * t645 + Ifges(7,4) * t644 + Ifges(7,5) * t673 + t648 * t669 - t649 * t705;
t661 = Ifges(6,5) * t695 + Ifges(6,6) * t694 + Ifges(6,3) * t706;
t663 = Ifges(6,1) * t695 + Ifges(6,4) * t694 + Ifges(6,5) * t706;
t604 = -mrSges(6,1) * t636 + mrSges(6,3) * t633 + Ifges(6,4) * t668 + Ifges(6,2) * t667 + Ifges(6,6) * t674 - pkin(5) * t766 + pkin(10) * t768 + t759 * t619 + t755 * t620 - t695 * t661 + t706 * t663;
t662 = Ifges(6,4) * t695 + Ifges(6,2) * t694 + Ifges(6,6) * t706;
t605 = mrSges(6,2) * t636 - mrSges(6,3) * t632 + Ifges(6,1) * t668 + Ifges(6,4) * t667 + Ifges(6,5) * t674 - pkin(10) * t618 - t619 * t755 + t620 * t759 + t661 * t694 - t662 * t706;
t679 = Ifges(5,5) * t707 - Ifges(5,6) * t706 + Ifges(5,3) * t740;
t680 = Ifges(5,4) * t707 - Ifges(5,2) * t706 + Ifges(5,6) * t740;
t589 = mrSges(5,2) * t658 - mrSges(5,3) * t638 + Ifges(5,1) * t675 - Ifges(5,4) * t674 + Ifges(5,5) * t730 - qJ(5) * t612 - t604 * t749 + t605 * t752 - t679 * t706 - t680 * t740;
t681 = Ifges(5,1) * t707 - Ifges(5,4) * t706 + Ifges(5,5) * t740;
t597 = Ifges(5,4) * t675 + Ifges(5,6) * t730 - t707 * t679 + t740 * t681 - mrSges(5,1) * t658 + mrSges(5,3) * t639 - Ifges(6,5) * t668 - Ifges(6,6) * t667 - t695 * t662 + t694 * t663 - mrSges(6,1) * t632 + mrSges(6,2) * t633 - Ifges(7,5) * t645 - Ifges(7,6) * t644 - Ifges(7,3) * t673 - t670 * t649 + t669 * t650 - mrSges(7,1) * t628 + mrSges(7,2) * t629 - pkin(5) * t618 - pkin(4) * t612 + (-Ifges(5,2) - Ifges(6,3)) * t674;
t698 = Ifges(4,5) * t727 + Ifges(4,6) * t726 - Ifges(4,3) * t773;
t700 = Ifges(4,1) * t727 + Ifges(4,4) * t726 - Ifges(4,5) * t773;
t578 = -mrSges(4,1) * t688 + mrSges(4,3) * t657 + Ifges(4,4) * t714 + Ifges(4,2) * t713 - Ifges(4,6) * t738 - pkin(3) * t765 + pkin(9) * t770 + t756 * t589 + t597 * t785 - t727 * t698 - t700 * t773;
t699 = Ifges(4,4) * t727 + Ifges(4,2) * t726 - Ifges(4,6) * t773;
t580 = mrSges(4,2) * t688 - mrSges(4,3) * t656 + Ifges(4,1) * t714 + Ifges(4,4) * t713 - Ifges(4,5) * t738 - pkin(9) * t603 + t589 * t785 - t756 * t597 + t726 * t698 + t699 * t773;
t716 = Ifges(3,3) * t746 + (Ifges(3,5) * t757 + Ifges(3,6) * t760) * t777;
t717 = Ifges(3,6) * t746 + (Ifges(3,4) * t757 + Ifges(3,2) * t760) * t777;
t577 = mrSges(3,2) * t719 - mrSges(3,3) * t703 + Ifges(3,1) * t737 + Ifges(3,4) * t738 + Ifges(3,5) * t745 - qJ(3) * t596 - t578 * t750 + t580 * t753 + t716 * t773 - t717 * t746;
t718 = Ifges(3,5) * t746 + (Ifges(3,1) * t757 + Ifges(3,4) * t760) * t777;
t579 = (Ifges(3,2) + Ifges(4,3)) * t738 - t749 * t605 - t752 * t604 + Ifges(3,6) * t745 + t746 * t718 - Ifges(5,3) * t730 + Ifges(3,4) * t737 + t726 * t700 - t727 * t699 - Ifges(4,6) * t713 - Ifges(4,5) * t714 - mrSges(3,1) * t719 - t706 * t681 - t707 * t680 + mrSges(3,3) * t704 + Ifges(5,6) * t674 - Ifges(5,5) * t675 - mrSges(4,1) * t656 + mrSges(4,2) * t657 - pkin(2) * t596 - mrSges(5,1) * t638 + mrSges(5,2) * t639 - t716 * t774 - qJ(5) * t769 - pkin(3) * t603 - pkin(4) * t764;
t767 = pkin(8) * t588 + t577 * t757 + t579 * t760;
t576 = Ifges(3,5) * t737 + Ifges(3,6) * t738 + Ifges(3,3) * t745 + mrSges(3,1) * t703 - mrSges(3,2) * t704 + t750 * t580 + t753 * t578 + pkin(2) * t763 + qJ(3) * t771 + (t717 * t757 - t718 * t760) * t777;
t575 = -mrSges(2,2) * g(3) - mrSges(2,3) * t741 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t762 + t577 * t760 - t579 * t757 + (-t583 * t751 - t584 * t754) * pkin(8);
t574 = mrSges(2,1) * g(3) + mrSges(2,3) * t742 + Ifges(2,5) * t762 + Ifges(2,6) * qJDD(1) - pkin(1) * t583 - t576 * t751 + t754 * t767;
t1 = [-m(1) * g(1) + t772; -m(1) * g(2) + t779; (-m(1) - m(2)) * g(3) + t583; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t779 - t758 * t574 + t761 * t575; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t772 + t761 * t574 + t758 * t575; -mrSges(1,1) * g(2) + mrSges(2,1) * t741 + mrSges(1,2) * g(1) - mrSges(2,2) * t742 + Ifges(2,3) * qJDD(1) + pkin(1) * t584 + t576 * t754 + t751 * t767;];
tauB  = t1;
