% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 10:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:57:46
% EndTime: 2019-05-07 09:58:15
% DurationCPUTime: 27.77s
% Computational Cost: add. (432177->387), mult. (971090->490), div. (0->0), fcn. (729760->12), ass. (0->149)
t742 = sin(qJ(2));
t747 = cos(qJ(2));
t763 = qJD(1) * qJD(2);
t725 = qJDD(1) * t742 + t747 * t763;
t743 = sin(qJ(1));
t748 = cos(qJ(1));
t731 = -g(1) * t748 - g(2) * t743;
t749 = qJD(1) ^ 2;
t720 = -pkin(1) * t749 + qJDD(1) * pkin(7) + t731;
t768 = t742 * t720;
t769 = pkin(2) * t749;
t682 = qJDD(2) * pkin(2) - t725 * pkin(8) - t768 + (pkin(8) * t763 + t742 * t769 - g(3)) * t747;
t707 = -g(3) * t742 + t747 * t720;
t726 = qJDD(1) * t747 - t742 * t763;
t766 = qJD(1) * t742;
t729 = qJD(2) * pkin(2) - pkin(8) * t766;
t736 = t747 ^ 2;
t683 = pkin(8) * t726 - qJD(2) * t729 - t736 * t769 + t707;
t741 = sin(qJ(3));
t746 = cos(qJ(3));
t658 = t746 * t682 - t741 * t683;
t717 = (-t741 * t742 + t746 * t747) * qJD(1);
t690 = qJD(3) * t717 + t725 * t746 + t726 * t741;
t718 = (t741 * t747 + t742 * t746) * qJD(1);
t734 = qJDD(2) + qJDD(3);
t735 = qJD(2) + qJD(3);
t641 = (t717 * t735 - t690) * qJ(4) + (t717 * t718 + t734) * pkin(3) + t658;
t659 = t741 * t682 + t746 * t683;
t689 = -qJD(3) * t718 - t725 * t741 + t726 * t746;
t709 = pkin(3) * t735 - qJ(4) * t718;
t713 = t717 ^ 2;
t645 = -pkin(3) * t713 + qJ(4) * t689 - t709 * t735 + t659;
t737 = sin(pkin(11));
t738 = cos(pkin(11));
t704 = t717 * t737 + t718 * t738;
t628 = -0.2e1 * qJD(4) * t704 + t641 * t738 - t737 * t645;
t703 = t717 * t738 - t737 * t718;
t629 = 0.2e1 * qJD(4) * t703 + t737 * t641 + t738 * t645;
t666 = t689 * t738 - t737 * t690;
t677 = -mrSges(5,1) * t703 + mrSges(5,2) * t704;
t693 = mrSges(5,1) * t735 - mrSges(5,3) * t704;
t678 = -pkin(4) * t703 - pkin(9) * t704;
t733 = t735 ^ 2;
t627 = -pkin(4) * t733 + pkin(9) * t734 + t678 * t703 + t629;
t730 = t743 * g(1) - t748 * g(2);
t755 = -qJDD(1) * pkin(1) - t730;
t691 = -t726 * pkin(2) + t729 * t766 + (-pkin(8) * t736 - pkin(7)) * t749 + t755;
t650 = -t689 * pkin(3) - t713 * qJ(4) + t718 * t709 + qJDD(4) + t691;
t667 = t689 * t737 + t690 * t738;
t632 = (-t703 * t735 - t667) * pkin(9) + (t704 * t735 - t666) * pkin(4) + t650;
t740 = sin(qJ(5));
t745 = cos(qJ(5));
t622 = -t740 * t627 + t745 * t632;
t687 = -t704 * t740 + t735 * t745;
t648 = qJD(5) * t687 + t667 * t745 + t734 * t740;
t665 = qJDD(5) - t666;
t688 = t704 * t745 + t735 * t740;
t697 = qJD(5) - t703;
t620 = (t687 * t697 - t648) * pkin(10) + (t687 * t688 + t665) * pkin(5) + t622;
t623 = t745 * t627 + t740 * t632;
t647 = -qJD(5) * t688 - t667 * t740 + t734 * t745;
t671 = pkin(5) * t697 - pkin(10) * t688;
t685 = t687 ^ 2;
t621 = -pkin(5) * t685 + pkin(10) * t647 - t671 * t697 + t623;
t739 = sin(qJ(6));
t744 = cos(qJ(6));
t618 = t620 * t744 - t621 * t739;
t661 = t687 * t744 - t688 * t739;
t635 = qJD(6) * t661 + t647 * t739 + t648 * t744;
t662 = t687 * t739 + t688 * t744;
t642 = -mrSges(7,1) * t661 + mrSges(7,2) * t662;
t694 = qJD(6) + t697;
t651 = -mrSges(7,2) * t694 + mrSges(7,3) * t661;
t660 = qJDD(6) + t665;
t616 = m(7) * t618 + mrSges(7,1) * t660 - mrSges(7,3) * t635 - t642 * t662 + t651 * t694;
t619 = t620 * t739 + t621 * t744;
t634 = -qJD(6) * t662 + t647 * t744 - t648 * t739;
t652 = mrSges(7,1) * t694 - mrSges(7,3) * t662;
t617 = m(7) * t619 - mrSges(7,2) * t660 + mrSges(7,3) * t634 + t642 * t661 - t652 * t694;
t608 = t744 * t616 + t739 * t617;
t668 = -mrSges(6,1) * t687 + mrSges(6,2) * t688;
t669 = -mrSges(6,2) * t697 + mrSges(6,3) * t687;
t606 = m(6) * t622 + mrSges(6,1) * t665 - mrSges(6,3) * t648 - t668 * t688 + t669 * t697 + t608;
t670 = mrSges(6,1) * t697 - mrSges(6,3) * t688;
t757 = -t616 * t739 + t744 * t617;
t607 = m(6) * t623 - mrSges(6,2) * t665 + mrSges(6,3) * t647 + t668 * t687 - t670 * t697 + t757;
t758 = -t606 * t740 + t745 * t607;
t601 = m(5) * t629 - mrSges(5,2) * t734 + mrSges(5,3) * t666 + t677 * t703 - t693 * t735 + t758;
t692 = -mrSges(5,2) * t735 + mrSges(5,3) * t703;
t626 = -pkin(4) * t734 - pkin(9) * t733 + t704 * t678 - t628;
t624 = -pkin(5) * t647 - pkin(10) * t685 + t671 * t688 + t626;
t753 = m(7) * t624 - t634 * mrSges(7,1) + mrSges(7,2) * t635 - t661 * t651 + t652 * t662;
t751 = -m(6) * t626 + t647 * mrSges(6,1) - mrSges(6,2) * t648 + t687 * t669 - t670 * t688 - t753;
t612 = m(5) * t628 + mrSges(5,1) * t734 - mrSges(5,3) * t667 - t677 * t704 + t692 * t735 + t751;
t594 = t737 * t601 + t738 * t612;
t705 = -mrSges(4,1) * t717 + mrSges(4,2) * t718;
t708 = -mrSges(4,2) * t735 + mrSges(4,3) * t717;
t592 = m(4) * t658 + mrSges(4,1) * t734 - mrSges(4,3) * t690 - t705 * t718 + t708 * t735 + t594;
t710 = mrSges(4,1) * t735 - mrSges(4,3) * t718;
t759 = t738 * t601 - t612 * t737;
t593 = m(4) * t659 - mrSges(4,2) * t734 + mrSges(4,3) * t689 + t705 * t717 - t710 * t735 + t759;
t587 = t746 * t592 + t741 * t593;
t706 = -t747 * g(3) - t768;
t724 = (-mrSges(3,1) * t747 + mrSges(3,2) * t742) * qJD(1);
t765 = qJD(1) * t747;
t728 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t765;
t585 = m(3) * t706 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t725 + qJD(2) * t728 - t724 * t766 + t587;
t727 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t766;
t760 = -t592 * t741 + t746 * t593;
t586 = m(3) * t707 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t726 - qJD(2) * t727 + t724 * t765 + t760;
t761 = -t585 * t742 + t747 * t586;
t579 = m(2) * t731 - mrSges(2,1) * t749 - qJDD(1) * mrSges(2,2) + t761;
t719 = -t749 * pkin(7) + t755;
t602 = t745 * t606 + t740 * t607;
t754 = m(5) * t650 - t666 * mrSges(5,1) + t667 * mrSges(5,2) - t703 * t692 + t704 * t693 + t602;
t752 = m(4) * t691 - t689 * mrSges(4,1) + mrSges(4,2) * t690 - t717 * t708 + t710 * t718 + t754;
t750 = -m(3) * t719 + t726 * mrSges(3,1) - mrSges(3,2) * t725 - t727 * t766 + t728 * t765 - t752;
t598 = m(2) * t730 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t749 + t750;
t767 = t743 * t579 + t748 * t598;
t580 = t747 * t585 + t742 * t586;
t762 = t748 * t579 - t598 * t743;
t716 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t742 + Ifges(3,4) * t747) * qJD(1);
t715 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t742 + Ifges(3,2) * t747) * qJD(1);
t714 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t742 + Ifges(3,6) * t747) * qJD(1);
t700 = Ifges(4,1) * t718 + Ifges(4,4) * t717 + Ifges(4,5) * t735;
t699 = Ifges(4,4) * t718 + Ifges(4,2) * t717 + Ifges(4,6) * t735;
t698 = Ifges(4,5) * t718 + Ifges(4,6) * t717 + Ifges(4,3) * t735;
t674 = Ifges(5,1) * t704 + Ifges(5,4) * t703 + Ifges(5,5) * t735;
t673 = Ifges(5,4) * t704 + Ifges(5,2) * t703 + Ifges(5,6) * t735;
t672 = Ifges(5,5) * t704 + Ifges(5,6) * t703 + Ifges(5,3) * t735;
t655 = Ifges(6,1) * t688 + Ifges(6,4) * t687 + Ifges(6,5) * t697;
t654 = Ifges(6,4) * t688 + Ifges(6,2) * t687 + Ifges(6,6) * t697;
t653 = Ifges(6,5) * t688 + Ifges(6,6) * t687 + Ifges(6,3) * t697;
t638 = Ifges(7,1) * t662 + Ifges(7,4) * t661 + Ifges(7,5) * t694;
t637 = Ifges(7,4) * t662 + Ifges(7,2) * t661 + Ifges(7,6) * t694;
t636 = Ifges(7,5) * t662 + Ifges(7,6) * t661 + Ifges(7,3) * t694;
t610 = mrSges(7,2) * t624 - mrSges(7,3) * t618 + Ifges(7,1) * t635 + Ifges(7,4) * t634 + Ifges(7,5) * t660 + t636 * t661 - t637 * t694;
t609 = -mrSges(7,1) * t624 + mrSges(7,3) * t619 + Ifges(7,4) * t635 + Ifges(7,2) * t634 + Ifges(7,6) * t660 - t636 * t662 + t638 * t694;
t596 = mrSges(6,2) * t626 - mrSges(6,3) * t622 + Ifges(6,1) * t648 + Ifges(6,4) * t647 + Ifges(6,5) * t665 - pkin(10) * t608 - t609 * t739 + t610 * t744 + t653 * t687 - t654 * t697;
t595 = -mrSges(6,1) * t626 + mrSges(6,3) * t623 + Ifges(6,4) * t648 + Ifges(6,2) * t647 + Ifges(6,6) * t665 - pkin(5) * t753 + pkin(10) * t757 + t744 * t609 + t739 * t610 - t688 * t653 + t697 * t655;
t588 = Ifges(5,4) * t667 + Ifges(5,2) * t666 + Ifges(5,6) * t734 - t704 * t672 + t735 * t674 - mrSges(5,1) * t650 + mrSges(5,3) * t629 - Ifges(6,5) * t648 - Ifges(6,6) * t647 - Ifges(6,3) * t665 - t688 * t654 + t687 * t655 - mrSges(6,1) * t622 + mrSges(6,2) * t623 - Ifges(7,5) * t635 - Ifges(7,6) * t634 - Ifges(7,3) * t660 - t662 * t637 + t661 * t638 - mrSges(7,1) * t618 + mrSges(7,2) * t619 - pkin(5) * t608 - pkin(4) * t602;
t581 = mrSges(5,2) * t650 - mrSges(5,3) * t628 + Ifges(5,1) * t667 + Ifges(5,4) * t666 + Ifges(5,5) * t734 - pkin(9) * t602 - t595 * t740 + t596 * t745 + t672 * t703 - t673 * t735;
t576 = mrSges(4,2) * t691 - mrSges(4,3) * t658 + Ifges(4,1) * t690 + Ifges(4,4) * t689 + Ifges(4,5) * t734 - qJ(4) * t594 + t581 * t738 - t588 * t737 + t698 * t717 - t699 * t735;
t575 = -mrSges(4,1) * t691 + mrSges(4,3) * t659 + Ifges(4,4) * t690 + Ifges(4,2) * t689 + Ifges(4,6) * t734 - pkin(3) * t754 + qJ(4) * t759 + t737 * t581 + t738 * t588 - t718 * t698 + t735 * t700;
t574 = mrSges(2,1) * g(3) - pkin(4) * t751 + Ifges(2,6) * qJDD(1) - pkin(3) * t594 - pkin(9) * t758 - pkin(2) * t587 + (-Ifges(4,3) - Ifges(5,3)) * t734 + (-t715 * t742 + t716 * t747) * qJD(1) - Ifges(3,3) * qJDD(2) + t749 * Ifges(2,5) - t745 * t595 - t740 * t596 + mrSges(2,3) * t731 - t718 * t699 - Ifges(3,5) * t725 - Ifges(3,6) * t726 - mrSges(3,1) * t706 + mrSges(3,2) * t707 + t717 * t700 + t703 * t674 - t704 * t673 - Ifges(4,6) * t689 - Ifges(4,5) * t690 - Ifges(5,6) * t666 - Ifges(5,5) * t667 - mrSges(4,1) * t658 + mrSges(4,2) * t659 - mrSges(5,1) * t628 + mrSges(5,2) * t629 - pkin(1) * t580;
t573 = mrSges(3,2) * t719 - mrSges(3,3) * t706 + Ifges(3,1) * t725 + Ifges(3,4) * t726 + Ifges(3,5) * qJDD(2) - pkin(8) * t587 - qJD(2) * t715 - t575 * t741 + t576 * t746 + t714 * t765;
t572 = -mrSges(3,1) * t719 + mrSges(3,3) * t707 + Ifges(3,4) * t725 + Ifges(3,2) * t726 + Ifges(3,6) * qJDD(2) - pkin(2) * t752 + pkin(8) * t760 + qJD(2) * t716 + t746 * t575 + t741 * t576 - t714 * t766;
t571 = -mrSges(2,2) * g(3) - mrSges(2,3) * t730 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t749 - pkin(7) * t580 - t572 * t742 + t573 * t747;
t1 = [-m(1) * g(1) + t762; -m(1) * g(2) + t767; (-m(1) - m(2)) * g(3) + t580; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t767 + t748 * t571 - t743 * t574; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t762 + t743 * t571 + t748 * t574; -mrSges(1,1) * g(2) + mrSges(2,1) * t730 + mrSges(1,2) * g(1) - mrSges(2,2) * t731 + Ifges(2,3) * qJDD(1) + pkin(1) * t750 + pkin(7) * t761 + t747 * t572 + t742 * t573;];
tauB  = t1;
