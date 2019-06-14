% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR2
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
% Datum: 2019-05-08 08:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:26:23
% EndTime: 2019-05-08 08:27:06
% DurationCPUTime: 30.13s
% Computational Cost: add. (481175->388), mult. (1037100->489), div. (0->0), fcn. (791480->12), ass. (0->151)
t747 = qJD(1) ^ 2;
t765 = pkin(2) * t747;
t740 = sin(qJ(1));
t746 = cos(qJ(1));
t727 = -g(1) * t746 - g(2) * t740;
t716 = -pkin(1) * t747 + qJDD(1) * pkin(7) + t727;
t739 = sin(qJ(2));
t764 = t739 * t716;
t745 = cos(qJ(2));
t760 = qJD(1) * qJD(2);
t721 = qJDD(1) * t739 + t745 * t760;
t679 = qJDD(2) * pkin(2) - t721 * pkin(8) - t764 + (pkin(8) * t760 + t739 * t765 - g(3)) * t745;
t703 = -g(3) * t739 + t745 * t716;
t722 = qJDD(1) * t745 - t739 * t760;
t762 = qJD(1) * t739;
t725 = qJD(2) * pkin(2) - pkin(8) * t762;
t734 = t745 ^ 2;
t680 = pkin(8) * t722 - qJD(2) * t725 - t734 * t765 + t703;
t738 = sin(qJ(3));
t744 = cos(qJ(3));
t661 = t744 * t679 - t738 * t680;
t713 = (-t738 * t739 + t744 * t745) * qJD(1);
t687 = qJD(3) * t713 + t721 * t744 + t722 * t738;
t714 = (t738 * t745 + t739 * t744) * qJD(1);
t732 = qJDD(2) + qJDD(3);
t733 = qJD(2) + qJD(3);
t638 = (t713 * t733 - t687) * pkin(9) + (t713 * t714 + t732) * pkin(3) + t661;
t662 = t738 * t679 + t744 * t680;
t686 = -qJD(3) * t714 - t721 * t738 + t722 * t744;
t706 = pkin(3) * t733 - pkin(9) * t714;
t709 = t713 ^ 2;
t643 = -pkin(3) * t709 + pkin(9) * t686 - t706 * t733 + t662;
t737 = sin(qJ(4));
t743 = cos(qJ(4));
t632 = t737 * t638 + t743 * t643;
t700 = t713 * t737 + t714 * t743;
t658 = -t700 * qJD(4) + t686 * t743 - t737 * t687;
t699 = t713 * t743 - t737 * t714;
t674 = -mrSges(5,1) * t699 + mrSges(5,2) * t700;
t730 = qJD(4) + t733;
t690 = mrSges(5,1) * t730 - mrSges(5,3) * t700;
t729 = qJDD(4) + t732;
t675 = -pkin(4) * t699 - pkin(10) * t700;
t728 = t730 ^ 2;
t624 = -pkin(4) * t728 + pkin(10) * t729 + t675 * t699 + t632;
t726 = t740 * g(1) - t746 * g(2);
t753 = -qJDD(1) * pkin(1) - t726;
t688 = -t722 * pkin(2) + t725 * t762 + (-pkin(8) * t734 - pkin(7)) * t747 + t753;
t647 = -t686 * pkin(3) - t709 * pkin(9) + t714 * t706 + t688;
t659 = qJD(4) * t699 + t686 * t737 + t687 * t743;
t630 = (-t699 * t730 - t659) * pkin(10) + (t700 * t730 - t658) * pkin(4) + t647;
t736 = sin(qJ(5));
t742 = cos(qJ(5));
t619 = -t736 * t624 + t742 * t630;
t683 = -t700 * t736 + t730 * t742;
t645 = qJD(5) * t683 + t659 * t742 + t729 * t736;
t656 = qJDD(5) - t658;
t684 = t700 * t742 + t730 * t736;
t696 = qJD(5) - t699;
t617 = (t683 * t696 - t645) * pkin(11) + (t683 * t684 + t656) * pkin(5) + t619;
t620 = t742 * t624 + t736 * t630;
t644 = -qJD(5) * t684 - t659 * t736 + t729 * t742;
t668 = pkin(5) * t696 - pkin(11) * t684;
t682 = t683 ^ 2;
t618 = -pkin(5) * t682 + pkin(11) * t644 - t668 * t696 + t620;
t735 = sin(qJ(6));
t741 = cos(qJ(6));
t615 = t617 * t741 - t618 * t735;
t663 = t683 * t741 - t684 * t735;
t629 = qJD(6) * t663 + t644 * t735 + t645 * t741;
t664 = t683 * t735 + t684 * t741;
t639 = -mrSges(7,1) * t663 + mrSges(7,2) * t664;
t691 = qJD(6) + t696;
t648 = -mrSges(7,2) * t691 + mrSges(7,3) * t663;
t650 = qJDD(6) + t656;
t613 = m(7) * t615 + mrSges(7,1) * t650 - mrSges(7,3) * t629 - t639 * t664 + t648 * t691;
t616 = t617 * t735 + t618 * t741;
t628 = -qJD(6) * t664 + t644 * t741 - t645 * t735;
t649 = mrSges(7,1) * t691 - mrSges(7,3) * t664;
t614 = m(7) * t616 - mrSges(7,2) * t650 + mrSges(7,3) * t628 + t639 * t663 - t649 * t691;
t605 = t741 * t613 + t735 * t614;
t665 = -mrSges(6,1) * t683 + mrSges(6,2) * t684;
t666 = -mrSges(6,2) * t696 + mrSges(6,3) * t683;
t603 = m(6) * t619 + mrSges(6,1) * t656 - mrSges(6,3) * t645 - t665 * t684 + t666 * t696 + t605;
t667 = mrSges(6,1) * t696 - mrSges(6,3) * t684;
t754 = -t613 * t735 + t741 * t614;
t604 = m(6) * t620 - mrSges(6,2) * t656 + mrSges(6,3) * t644 + t665 * t683 - t667 * t696 + t754;
t755 = -t603 * t736 + t742 * t604;
t598 = m(5) * t632 - mrSges(5,2) * t729 + mrSges(5,3) * t658 + t674 * t699 - t690 * t730 + t755;
t631 = t638 * t743 - t737 * t643;
t689 = -mrSges(5,2) * t730 + mrSges(5,3) * t699;
t623 = -pkin(4) * t729 - pkin(10) * t728 + t700 * t675 - t631;
t621 = -pkin(5) * t644 - pkin(11) * t682 + t668 * t684 + t623;
t751 = m(7) * t621 - t628 * mrSges(7,1) + mrSges(7,2) * t629 - t663 * t648 + t649 * t664;
t749 = -m(6) * t623 + t644 * mrSges(6,1) - mrSges(6,2) * t645 + t683 * t666 - t667 * t684 - t751;
t609 = m(5) * t631 + mrSges(5,1) * t729 - mrSges(5,3) * t659 - t674 * t700 + t689 * t730 + t749;
t591 = t737 * t598 + t743 * t609;
t701 = -mrSges(4,1) * t713 + mrSges(4,2) * t714;
t704 = -mrSges(4,2) * t733 + mrSges(4,3) * t713;
t589 = m(4) * t661 + mrSges(4,1) * t732 - mrSges(4,3) * t687 - t701 * t714 + t704 * t733 + t591;
t705 = mrSges(4,1) * t733 - mrSges(4,3) * t714;
t756 = t743 * t598 - t609 * t737;
t590 = m(4) * t662 - mrSges(4,2) * t732 + mrSges(4,3) * t686 + t701 * t713 - t705 * t733 + t756;
t584 = t744 * t589 + t738 * t590;
t702 = -t745 * g(3) - t764;
t720 = (-mrSges(3,1) * t745 + mrSges(3,2) * t739) * qJD(1);
t761 = qJD(1) * t745;
t724 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t761;
t582 = m(3) * t702 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t721 + qJD(2) * t724 - t720 * t762 + t584;
t723 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t762;
t757 = -t589 * t738 + t744 * t590;
t583 = m(3) * t703 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t722 - qJD(2) * t723 + t720 * t761 + t757;
t758 = -t582 * t739 + t745 * t583;
t576 = m(2) * t727 - mrSges(2,1) * t747 - qJDD(1) * mrSges(2,2) + t758;
t715 = -t747 * pkin(7) + t753;
t599 = t742 * t603 + t736 * t604;
t752 = m(5) * t647 - t658 * mrSges(5,1) + t659 * mrSges(5,2) - t699 * t689 + t700 * t690 + t599;
t750 = m(4) * t688 - t686 * mrSges(4,1) + mrSges(4,2) * t687 - t713 * t704 + t705 * t714 + t752;
t748 = -m(3) * t715 + t722 * mrSges(3,1) - mrSges(3,2) * t721 - t723 * t762 + t724 * t761 - t750;
t595 = m(2) * t726 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t747 + t748;
t763 = t740 * t576 + t746 * t595;
t577 = t745 * t582 + t739 * t583;
t759 = t746 * t576 - t595 * t740;
t712 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t739 + Ifges(3,4) * t745) * qJD(1);
t711 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t739 + Ifges(3,2) * t745) * qJD(1);
t710 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t739 + Ifges(3,6) * t745) * qJD(1);
t694 = Ifges(4,1) * t714 + Ifges(4,4) * t713 + Ifges(4,5) * t733;
t693 = Ifges(4,4) * t714 + Ifges(4,2) * t713 + Ifges(4,6) * t733;
t692 = Ifges(4,5) * t714 + Ifges(4,6) * t713 + Ifges(4,3) * t733;
t671 = Ifges(5,1) * t700 + Ifges(5,4) * t699 + Ifges(5,5) * t730;
t670 = Ifges(5,4) * t700 + Ifges(5,2) * t699 + Ifges(5,6) * t730;
t669 = Ifges(5,5) * t700 + Ifges(5,6) * t699 + Ifges(5,3) * t730;
t653 = Ifges(6,1) * t684 + Ifges(6,4) * t683 + Ifges(6,5) * t696;
t652 = Ifges(6,4) * t684 + Ifges(6,2) * t683 + Ifges(6,6) * t696;
t651 = Ifges(6,5) * t684 + Ifges(6,6) * t683 + Ifges(6,3) * t696;
t635 = Ifges(7,1) * t664 + Ifges(7,4) * t663 + Ifges(7,5) * t691;
t634 = Ifges(7,4) * t664 + Ifges(7,2) * t663 + Ifges(7,6) * t691;
t633 = Ifges(7,5) * t664 + Ifges(7,6) * t663 + Ifges(7,3) * t691;
t607 = mrSges(7,2) * t621 - mrSges(7,3) * t615 + Ifges(7,1) * t629 + Ifges(7,4) * t628 + Ifges(7,5) * t650 + t633 * t663 - t634 * t691;
t606 = -mrSges(7,1) * t621 + mrSges(7,3) * t616 + Ifges(7,4) * t629 + Ifges(7,2) * t628 + Ifges(7,6) * t650 - t633 * t664 + t635 * t691;
t593 = mrSges(6,2) * t623 - mrSges(6,3) * t619 + Ifges(6,1) * t645 + Ifges(6,4) * t644 + Ifges(6,5) * t656 - pkin(11) * t605 - t606 * t735 + t607 * t741 + t651 * t683 - t652 * t696;
t592 = -mrSges(6,1) * t623 + mrSges(6,3) * t620 + Ifges(6,4) * t645 + Ifges(6,2) * t644 + Ifges(6,6) * t656 - pkin(5) * t751 + pkin(11) * t754 + t741 * t606 + t735 * t607 - t684 * t651 + t696 * t653;
t585 = Ifges(5,4) * t659 + Ifges(5,2) * t658 + Ifges(5,6) * t729 - t700 * t669 + t730 * t671 - mrSges(5,1) * t647 + mrSges(5,3) * t632 - Ifges(6,5) * t645 - Ifges(6,6) * t644 - Ifges(6,3) * t656 - t684 * t652 + t683 * t653 - mrSges(6,1) * t619 + mrSges(6,2) * t620 - Ifges(7,5) * t629 - Ifges(7,6) * t628 - Ifges(7,3) * t650 - t664 * t634 + t663 * t635 - mrSges(7,1) * t615 + mrSges(7,2) * t616 - pkin(5) * t605 - pkin(4) * t599;
t578 = mrSges(5,2) * t647 - mrSges(5,3) * t631 + Ifges(5,1) * t659 + Ifges(5,4) * t658 + Ifges(5,5) * t729 - pkin(10) * t599 - t592 * t736 + t593 * t742 + t669 * t699 - t670 * t730;
t573 = mrSges(4,2) * t688 - mrSges(4,3) * t661 + Ifges(4,1) * t687 + Ifges(4,4) * t686 + Ifges(4,5) * t732 - pkin(9) * t591 + t578 * t743 - t585 * t737 + t692 * t713 - t693 * t733;
t572 = -mrSges(4,1) * t688 + mrSges(4,3) * t662 + Ifges(4,4) * t687 + Ifges(4,2) * t686 + Ifges(4,6) * t732 - pkin(3) * t752 + pkin(9) * t756 + t737 * t578 + t743 * t585 - t714 * t692 + t733 * t694;
t571 = -pkin(2) * t584 - Ifges(3,3) * qJDD(2) + (-t711 * t739 + t712 * t745) * qJD(1) + t747 * Ifges(2,5) - t742 * t592 - Ifges(4,3) * t732 - t736 * t593 + mrSges(2,3) * t727 - Ifges(5,3) * t729 - t714 * t693 - Ifges(3,5) * t721 - Ifges(3,6) * t722 - mrSges(3,1) * t702 + mrSges(3,2) * t703 + t713 * t694 + t699 * t671 - t700 * t670 - Ifges(4,6) * t686 - Ifges(4,5) * t687 - Ifges(5,6) * t658 - Ifges(5,5) * t659 - mrSges(4,1) * t661 + mrSges(4,2) * t662 - mrSges(5,1) * t631 + mrSges(5,2) * t632 + mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) - pkin(3) * t591 - pkin(1) * t577 - pkin(4) * t749 - pkin(10) * t755;
t570 = mrSges(3,2) * t715 - mrSges(3,3) * t702 + Ifges(3,1) * t721 + Ifges(3,4) * t722 + Ifges(3,5) * qJDD(2) - pkin(8) * t584 - qJD(2) * t711 - t572 * t738 + t573 * t744 + t710 * t761;
t569 = -mrSges(3,1) * t715 + mrSges(3,3) * t703 + Ifges(3,4) * t721 + Ifges(3,2) * t722 + Ifges(3,6) * qJDD(2) - pkin(2) * t750 + pkin(8) * t757 + qJD(2) * t712 + t744 * t572 + t738 * t573 - t710 * t762;
t568 = -mrSges(2,2) * g(3) - mrSges(2,3) * t726 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t747 - pkin(7) * t577 - t569 * t739 + t570 * t745;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t763; (-m(1) - m(2)) * g(3) + t577; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t763 + t746 * t568 - t740 * t571; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t759 + t740 * t568 + t746 * t571; -mrSges(1,1) * g(2) + mrSges(2,1) * t726 + mrSges(1,2) * g(1) - mrSges(2,2) * t727 + Ifges(2,3) * qJDD(1) + pkin(1) * t748 + pkin(7) * t758 + t745 * t569 + t739 * t570;];
tauB  = t1;
