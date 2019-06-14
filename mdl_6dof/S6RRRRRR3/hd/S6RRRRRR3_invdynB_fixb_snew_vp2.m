% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 09:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:54:17
% EndTime: 2019-05-08 08:55:09
% DurationCPUTime: 31.46s
% Computational Cost: add. (527763->388), mult. (1051182->488), div. (0->0), fcn. (783476->12), ass. (0->151)
t749 = qJD(1) ^ 2;
t767 = pkin(2) * t749;
t742 = sin(qJ(1));
t748 = cos(qJ(1));
t730 = -g(1) * t748 - g(2) * t742;
t718 = -pkin(1) * t749 + qJDD(1) * pkin(7) + t730;
t741 = sin(qJ(2));
t766 = t741 * t718;
t747 = cos(qJ(2));
t762 = qJD(1) * qJD(2);
t724 = qJDD(1) * t741 + t747 * t762;
t681 = qJDD(2) * pkin(2) - t724 * pkin(8) - t766 + (pkin(8) * t762 + t741 * t767 - g(3)) * t747;
t704 = -g(3) * t741 + t747 * t718;
t725 = qJDD(1) * t747 - t741 * t762;
t764 = qJD(1) * t741;
t728 = qJD(2) * pkin(2) - pkin(8) * t764;
t736 = t747 ^ 2;
t682 = pkin(8) * t725 - qJD(2) * t728 - t736 * t767 + t704;
t740 = sin(qJ(3));
t746 = cos(qJ(3));
t662 = t740 * t681 + t746 * t682;
t716 = (t740 * t747 + t741 * t746) * qJD(1);
t688 = -t716 * qJD(3) - t740 * t724 + t725 * t746;
t763 = qJD(1) * t747;
t715 = -t740 * t764 + t746 * t763;
t698 = -mrSges(4,1) * t715 + mrSges(4,2) * t716;
t735 = qJD(2) + qJD(3);
t706 = mrSges(4,1) * t735 - mrSges(4,3) * t716;
t734 = qJDD(2) + qJDD(3);
t689 = qJD(3) * t715 + t724 * t746 + t725 * t740;
t729 = t742 * g(1) - t748 * g(2);
t754 = -qJDD(1) * pkin(1) - t729;
t690 = -t725 * pkin(2) + t728 * t764 + (-pkin(8) * t736 - pkin(7)) * t749 + t754;
t648 = (-t715 * t735 - t689) * pkin(9) + (t716 * t735 - t688) * pkin(3) + t690;
t699 = -pkin(3) * t715 - pkin(9) * t716;
t733 = t735 ^ 2;
t651 = -pkin(3) * t733 + pkin(9) * t734 + t699 * t715 + t662;
t739 = sin(qJ(4));
t745 = cos(qJ(4));
t634 = t745 * t648 - t739 * t651;
t701 = -t716 * t739 + t735 * t745;
t665 = qJD(4) * t701 + t689 * t745 + t734 * t739;
t687 = qJDD(4) - t688;
t702 = t716 * t745 + t735 * t739;
t711 = qJD(4) - t715;
t631 = (t701 * t711 - t665) * pkin(10) + (t701 * t702 + t687) * pkin(4) + t634;
t635 = t739 * t648 + t745 * t651;
t664 = -qJD(4) * t702 - t689 * t739 + t734 * t745;
t693 = pkin(4) * t711 - pkin(10) * t702;
t700 = t701 ^ 2;
t633 = -pkin(4) * t700 + pkin(10) * t664 - t693 * t711 + t635;
t738 = sin(qJ(5));
t744 = cos(qJ(5));
t621 = t744 * t631 - t738 * t633;
t675 = t701 * t744 - t702 * t738;
t645 = qJD(5) * t675 + t664 * t738 + t665 * t744;
t676 = t701 * t738 + t702 * t744;
t684 = qJDD(5) + t687;
t709 = qJD(5) + t711;
t619 = (t675 * t709 - t645) * pkin(11) + (t675 * t676 + t684) * pkin(5) + t621;
t622 = t738 * t631 + t744 * t633;
t644 = -qJD(5) * t676 + t664 * t744 - t665 * t738;
t668 = pkin(5) * t709 - pkin(11) * t676;
t674 = t675 ^ 2;
t620 = -pkin(5) * t674 + pkin(11) * t644 - t668 * t709 + t622;
t737 = sin(qJ(6));
t743 = cos(qJ(6));
t617 = t619 * t743 - t620 * t737;
t658 = t675 * t743 - t676 * t737;
t628 = qJD(6) * t658 + t644 * t737 + t645 * t743;
t659 = t675 * t737 + t676 * t743;
t642 = -mrSges(7,1) * t658 + mrSges(7,2) * t659;
t707 = qJD(6) + t709;
t652 = -mrSges(7,2) * t707 + mrSges(7,3) * t658;
t683 = qJDD(6) + t684;
t615 = m(7) * t617 + mrSges(7,1) * t683 - mrSges(7,3) * t628 - t642 * t659 + t652 * t707;
t618 = t619 * t737 + t620 * t743;
t627 = -qJD(6) * t659 + t644 * t743 - t645 * t737;
t653 = mrSges(7,1) * t707 - mrSges(7,3) * t659;
t616 = m(7) * t618 - mrSges(7,2) * t683 + mrSges(7,3) * t627 + t642 * t658 - t653 * t707;
t607 = t743 * t615 + t737 * t616;
t660 = -mrSges(6,1) * t675 + mrSges(6,2) * t676;
t666 = -mrSges(6,2) * t709 + mrSges(6,3) * t675;
t605 = m(6) * t621 + mrSges(6,1) * t684 - mrSges(6,3) * t645 - t660 * t676 + t666 * t709 + t607;
t667 = mrSges(6,1) * t709 - mrSges(6,3) * t676;
t756 = -t615 * t737 + t743 * t616;
t606 = m(6) * t622 - mrSges(6,2) * t684 + mrSges(6,3) * t644 + t660 * t675 - t667 * t709 + t756;
t601 = t744 * t605 + t738 * t606;
t680 = -mrSges(5,1) * t701 + mrSges(5,2) * t702;
t691 = -mrSges(5,2) * t711 + mrSges(5,3) * t701;
t599 = m(5) * t634 + mrSges(5,1) * t687 - mrSges(5,3) * t665 - t680 * t702 + t691 * t711 + t601;
t692 = mrSges(5,1) * t711 - mrSges(5,3) * t702;
t757 = -t605 * t738 + t744 * t606;
t600 = m(5) * t635 - mrSges(5,2) * t687 + mrSges(5,3) * t664 + t680 * t701 - t692 * t711 + t757;
t758 = -t599 * t739 + t745 * t600;
t592 = m(4) * t662 - mrSges(4,2) * t734 + mrSges(4,3) * t688 + t698 * t715 - t706 * t735 + t758;
t661 = t681 * t746 - t740 * t682;
t705 = -mrSges(4,2) * t735 + mrSges(4,3) * t715;
t650 = -pkin(3) * t734 - pkin(9) * t733 + t716 * t699 - t661;
t636 = -pkin(4) * t664 - pkin(10) * t700 + t702 * t693 + t650;
t624 = -pkin(5) * t644 - pkin(11) * t674 + t668 * t676 + t636;
t755 = m(7) * t624 - t627 * mrSges(7,1) + t628 * mrSges(7,2) - t658 * t652 + t659 * t653;
t752 = m(6) * t636 - t644 * mrSges(6,1) + mrSges(6,2) * t645 - t675 * t666 + t667 * t676 + t755;
t750 = -m(5) * t650 + t664 * mrSges(5,1) - mrSges(5,2) * t665 + t701 * t691 - t692 * t702 - t752;
t611 = m(4) * t661 + mrSges(4,1) * t734 - mrSges(4,3) * t689 - t698 * t716 + t705 * t735 + t750;
t587 = t740 * t592 + t746 * t611;
t703 = -t747 * g(3) - t766;
t723 = (-mrSges(3,1) * t747 + mrSges(3,2) * t741) * qJD(1);
t727 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t763;
t585 = m(3) * t703 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t724 + qJD(2) * t727 - t723 * t764 + t587;
t726 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t764;
t759 = t746 * t592 - t611 * t740;
t586 = m(3) * t704 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t725 - qJD(2) * t726 + t723 * t763 + t759;
t760 = -t585 * t741 + t747 * t586;
t577 = m(2) * t730 - mrSges(2,1) * t749 - qJDD(1) * mrSges(2,2) + t760;
t717 = -t749 * pkin(7) + t754;
t593 = t745 * t599 + t739 * t600;
t753 = m(4) * t690 - t688 * mrSges(4,1) + mrSges(4,2) * t689 - t715 * t705 + t706 * t716 + t593;
t751 = -m(3) * t717 + t725 * mrSges(3,1) - mrSges(3,2) * t724 - t726 * t764 + t727 * t763 - t753;
t589 = m(2) * t729 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t749 + t751;
t765 = t742 * t577 + t748 * t589;
t578 = t747 * t585 + t741 * t586;
t761 = t748 * t577 - t589 * t742;
t714 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t741 + Ifges(3,4) * t747) * qJD(1);
t713 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t741 + Ifges(3,2) * t747) * qJD(1);
t712 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t741 + Ifges(3,6) * t747) * qJD(1);
t696 = Ifges(4,1) * t716 + Ifges(4,4) * t715 + Ifges(4,5) * t735;
t695 = Ifges(4,4) * t716 + Ifges(4,2) * t715 + Ifges(4,6) * t735;
t694 = Ifges(4,5) * t716 + Ifges(4,6) * t715 + Ifges(4,3) * t735;
t671 = Ifges(5,1) * t702 + Ifges(5,4) * t701 + Ifges(5,5) * t711;
t670 = Ifges(5,4) * t702 + Ifges(5,2) * t701 + Ifges(5,6) * t711;
t669 = Ifges(5,5) * t702 + Ifges(5,6) * t701 + Ifges(5,3) * t711;
t656 = Ifges(6,1) * t676 + Ifges(6,4) * t675 + Ifges(6,5) * t709;
t655 = Ifges(6,4) * t676 + Ifges(6,2) * t675 + Ifges(6,6) * t709;
t654 = Ifges(6,5) * t676 + Ifges(6,6) * t675 + Ifges(6,3) * t709;
t639 = Ifges(7,1) * t659 + Ifges(7,4) * t658 + Ifges(7,5) * t707;
t638 = Ifges(7,4) * t659 + Ifges(7,2) * t658 + Ifges(7,6) * t707;
t637 = Ifges(7,5) * t659 + Ifges(7,6) * t658 + Ifges(7,3) * t707;
t609 = mrSges(7,2) * t624 - mrSges(7,3) * t617 + Ifges(7,1) * t628 + Ifges(7,4) * t627 + Ifges(7,5) * t683 + t637 * t658 - t638 * t707;
t608 = -mrSges(7,1) * t624 + mrSges(7,3) * t618 + Ifges(7,4) * t628 + Ifges(7,2) * t627 + Ifges(7,6) * t683 - t637 * t659 + t639 * t707;
t595 = mrSges(6,2) * t636 - mrSges(6,3) * t621 + Ifges(6,1) * t645 + Ifges(6,4) * t644 + Ifges(6,5) * t684 - pkin(11) * t607 - t608 * t737 + t609 * t743 + t654 * t675 - t655 * t709;
t594 = -mrSges(6,1) * t636 + mrSges(6,3) * t622 + Ifges(6,4) * t645 + Ifges(6,2) * t644 + Ifges(6,6) * t684 - pkin(5) * t755 + pkin(11) * t756 + t743 * t608 + t737 * t609 - t676 * t654 + t709 * t656;
t581 = mrSges(5,2) * t650 - mrSges(5,3) * t634 + Ifges(5,1) * t665 + Ifges(5,4) * t664 + Ifges(5,5) * t687 - pkin(10) * t601 - t594 * t738 + t595 * t744 + t669 * t701 - t670 * t711;
t580 = -mrSges(5,1) * t650 + mrSges(5,3) * t635 + Ifges(5,4) * t665 + Ifges(5,2) * t664 + Ifges(5,6) * t687 - pkin(4) * t752 + pkin(10) * t757 + t744 * t594 + t738 * t595 - t702 * t669 + t711 * t671;
t579 = Ifges(4,6) * t734 + t735 * t696 - t716 * t694 + t701 * t671 - t702 * t670 - mrSges(4,1) * t690 - Ifges(7,3) * t683 - Ifges(6,3) * t684 - Ifges(5,3) * t687 + Ifges(4,2) * t688 + Ifges(4,4) * t689 + t675 * t656 - t676 * t655 + mrSges(4,3) * t662 - Ifges(5,6) * t664 - Ifges(5,5) * t665 + t658 * t639 - t659 * t638 - Ifges(6,6) * t644 - Ifges(6,5) * t645 + mrSges(5,2) * t635 - mrSges(5,1) * t634 - Ifges(7,6) * t627 - Ifges(7,5) * t628 - mrSges(6,1) * t621 + mrSges(6,2) * t622 + mrSges(7,2) * t618 - mrSges(7,1) * t617 - pkin(5) * t607 - pkin(4) * t601 - pkin(3) * t593;
t574 = mrSges(4,2) * t690 - mrSges(4,3) * t661 + Ifges(4,1) * t689 + Ifges(4,4) * t688 + Ifges(4,5) * t734 - pkin(9) * t593 - t580 * t739 + t581 * t745 + t694 * t715 - t695 * t735;
t573 = mrSges(3,2) * t717 - mrSges(3,3) * t703 + Ifges(3,1) * t724 + Ifges(3,4) * t725 + Ifges(3,5) * qJDD(2) - pkin(8) * t587 - qJD(2) * t713 + t574 * t746 - t579 * t740 + t712 * t763;
t572 = -pkin(1) * t578 + mrSges(2,3) * t730 - pkin(2) * t587 - Ifges(3,6) * t725 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t703 + mrSges(3,2) * t704 - Ifges(3,5) * t724 - t745 * t580 - pkin(3) * t750 - pkin(9) * t758 - t739 * t581 - Ifges(4,5) * t689 - Ifges(4,6) * t688 - Ifges(4,3) * t734 - mrSges(4,1) * t661 + mrSges(4,2) * t662 + t749 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - t716 * t695 + t715 * t696 + (-t713 * t741 + t714 * t747) * qJD(1);
t571 = -mrSges(3,1) * t717 + mrSges(3,3) * t704 + Ifges(3,4) * t724 + Ifges(3,2) * t725 + Ifges(3,6) * qJDD(2) - pkin(2) * t753 + pkin(8) * t759 + qJD(2) * t714 + t740 * t574 + t746 * t579 - t712 * t764;
t570 = -mrSges(2,2) * g(3) - mrSges(2,3) * t729 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t749 - pkin(7) * t578 - t571 * t741 + t573 * t747;
t1 = [-m(1) * g(1) + t761; -m(1) * g(2) + t765; (-m(1) - m(2)) * g(3) + t578; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t765 + t748 * t570 - t742 * t572; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t761 + t742 * t570 + t748 * t572; -mrSges(1,1) * g(2) + mrSges(2,1) * t729 + mrSges(1,2) * g(1) - mrSges(2,2) * t730 + Ifges(2,3) * qJDD(1) + pkin(1) * t751 + pkin(7) * t760 + t747 * t571 + t741 * t573;];
tauB  = t1;
