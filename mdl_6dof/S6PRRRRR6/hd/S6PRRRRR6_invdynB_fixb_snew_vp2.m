% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:35:31
% EndTime: 2019-05-05 12:37:23
% DurationCPUTime: 111.72s
% Computational Cost: add. (1959013->371), mult. (4608945->504), div. (0->0), fcn. (3892652->18), ass. (0->172)
t733 = sin(pkin(14));
t737 = cos(pkin(14));
t726 = g(1) * t733 - g(2) * t737;
t727 = -g(1) * t737 - g(2) * t733;
t732 = -g(3) + qJDD(1);
t750 = cos(qJ(2));
t740 = cos(pkin(6));
t745 = sin(qJ(2));
t772 = t740 * t745;
t736 = sin(pkin(6));
t778 = t736 * t745;
t701 = t726 * t772 + t750 * t727 + t732 * t778;
t751 = qJD(2) ^ 2;
t735 = sin(pkin(7));
t785 = pkin(10) * t735;
t698 = -pkin(2) * t751 + qJDD(2) * t785 + t701;
t739 = cos(pkin(7));
t731 = qJD(2) * t739 + qJD(3);
t734 = sin(pkin(8));
t738 = cos(pkin(8));
t749 = cos(qJ(3));
t768 = qJD(2) * t735;
t764 = t749 * t768;
t713 = (t731 * t734 + t738 * t764) * pkin(11);
t744 = sin(qJ(3));
t784 = pkin(11) * t734;
t717 = (-pkin(3) * t749 - t744 * t784) * t768;
t767 = qJD(2) * qJD(3);
t721 = (qJDD(2) * t744 + t749 * t767) * t735;
t730 = qJDD(2) * t739 + qJDD(3);
t771 = t740 * t750;
t777 = t736 * t750;
t700 = t726 * t771 - t727 * t745 + t732 * t777;
t697 = qJDD(2) * pkin(2) + t751 * t785 + t700;
t716 = -t726 * t736 + t732 * t740;
t773 = t739 * t749;
t779 = t735 * t749;
t769 = t697 * t773 + t716 * t779;
t783 = pkin(11) * t738;
t656 = -t721 * t783 + t730 * pkin(3) + t731 * t713 + (-t717 * t768 - t698) * t744 + t769;
t774 = t739 * t744;
t780 = t735 * t744;
t674 = t697 * t774 + t749 * t698 + t716 * t780;
t765 = t744 * t768;
t715 = pkin(3) * t731 - t765 * t783;
t722 = (qJDD(2) * t749 - t744 * t767) * t735;
t759 = t722 * t738 + t730 * t734;
t657 = pkin(11) * t759 - t731 * t715 + t717 * t764 + t674;
t711 = t739 * t716;
t670 = -t721 * t784 - t722 * pkin(3) + t711 + (-t697 + (-t713 * t749 + t715 * t744) * qJD(2)) * t735;
t743 = sin(qJ(4));
t748 = cos(qJ(4));
t642 = -t743 * t657 + (t656 * t738 + t670 * t734) * t748;
t775 = t738 * t749;
t782 = t734 * t743;
t704 = t731 * t782 + (t743 * t775 + t744 * t748) * t768;
t684 = -t704 * qJD(4) - t743 * t721 + t748 * t759;
t781 = t734 * t748;
t703 = (-t743 * t744 + t748 * t775) * t768 + t731 * t781;
t776 = t738 * t743;
t643 = t656 * t776 + t748 * t657 + t670 * t782;
t686 = -mrSges(5,1) * t703 + mrSges(5,2) * t704;
t714 = t731 * t738 - t734 * t764 + qJD(4);
t693 = mrSges(5,1) * t714 - mrSges(5,3) * t704;
t705 = -t722 * t734 + t730 * t738 + qJDD(4);
t687 = -pkin(4) * t703 - pkin(12) * t704;
t712 = t714 ^ 2;
t639 = -pkin(4) * t712 + pkin(12) * t705 + t687 * t703 + t643;
t647 = -t734 * t656 + t738 * t670;
t685 = t703 * qJD(4) + t748 * t721 + t743 * t759;
t641 = (-t703 * t714 - t685) * pkin(12) + (t704 * t714 - t684) * pkin(4) + t647;
t742 = sin(qJ(5));
t747 = cos(qJ(5));
t635 = t747 * t639 + t742 * t641;
t691 = t704 * t747 + t714 * t742;
t660 = -qJD(5) * t691 - t685 * t742 + t705 * t747;
t690 = -t704 * t742 + t714 * t747;
t671 = -mrSges(6,1) * t690 + mrSges(6,2) * t691;
t702 = qJD(5) - t703;
t679 = mrSges(6,1) * t702 - mrSges(6,3) * t691;
t683 = qJDD(5) - t684;
t672 = -pkin(5) * t690 - pkin(13) * t691;
t699 = t702 ^ 2;
t633 = -pkin(5) * t699 + pkin(13) * t683 + t672 * t690 + t635;
t638 = -t705 * pkin(4) - t712 * pkin(12) + t704 * t687 - t642;
t661 = qJD(5) * t690 + t685 * t747 + t705 * t742;
t636 = (-t690 * t702 - t661) * pkin(13) + (t691 * t702 - t660) * pkin(5) + t638;
t741 = sin(qJ(6));
t746 = cos(qJ(6));
t630 = -t633 * t741 + t636 * t746;
t676 = -t691 * t741 + t702 * t746;
t646 = qJD(6) * t676 + t661 * t746 + t683 * t741;
t677 = t691 * t746 + t702 * t741;
t652 = -mrSges(7,1) * t676 + mrSges(7,2) * t677;
t659 = qJDD(6) - t660;
t689 = qJD(6) - t690;
t662 = -mrSges(7,2) * t689 + mrSges(7,3) * t676;
t628 = m(7) * t630 + mrSges(7,1) * t659 - mrSges(7,3) * t646 - t652 * t677 + t662 * t689;
t631 = t633 * t746 + t636 * t741;
t645 = -qJD(6) * t677 - t661 * t741 + t683 * t746;
t663 = mrSges(7,1) * t689 - mrSges(7,3) * t677;
t629 = m(7) * t631 - mrSges(7,2) * t659 + mrSges(7,3) * t645 + t652 * t676 - t663 * t689;
t761 = -t628 * t741 + t746 * t629;
t621 = m(6) * t635 - mrSges(6,2) * t683 + mrSges(6,3) * t660 + t671 * t690 - t679 * t702 + t761;
t634 = -t639 * t742 + t641 * t747;
t678 = -mrSges(6,2) * t702 + mrSges(6,3) * t690;
t632 = -pkin(5) * t683 - pkin(13) * t699 + t672 * t691 - t634;
t753 = -m(7) * t632 + t645 * mrSges(7,1) - mrSges(7,2) * t646 + t676 * t662 - t663 * t677;
t626 = m(6) * t634 + mrSges(6,1) * t683 - mrSges(6,3) * t661 - t671 * t691 + t678 * t702 + t753;
t762 = t747 * t621 - t626 * t742;
t612 = m(5) * t643 - mrSges(5,2) * t705 + mrSges(5,3) * t684 + t686 * t703 - t693 * t714 + t762;
t615 = t742 * t621 + t747 * t626;
t692 = -mrSges(5,2) * t714 + mrSges(5,3) * t703;
t614 = m(5) * t647 - mrSges(5,1) * t684 + mrSges(5,2) * t685 - t692 * t703 + t693 * t704 + t615;
t622 = t628 * t746 + t629 * t741;
t752 = -m(6) * t638 + t660 * mrSges(6,1) - mrSges(6,2) * t661 + t690 * t678 - t679 * t691 - t622;
t618 = m(5) * t642 + mrSges(5,1) * t705 - mrSges(5,3) * t685 - t686 * t704 + t692 * t714 + t752;
t601 = t738 * t748 * t618 + t612 * t776 - t614 * t734;
t673 = -t698 * t744 + t769;
t719 = -mrSges(4,2) * t731 + mrSges(4,3) * t764;
t720 = (-mrSges(4,1) * t749 + mrSges(4,2) * t744) * t768;
t597 = m(4) * t673 + mrSges(4,1) * t730 - mrSges(4,3) * t721 + t719 * t731 - t720 * t765 + t601;
t600 = t612 * t782 + t738 * t614 + t618 * t781;
t688 = -t735 * t697 + t711;
t718 = mrSges(4,1) * t731 - mrSges(4,3) * t765;
t599 = m(4) * t688 - t722 * mrSges(4,1) + t721 * mrSges(4,2) + (t718 * t744 - t719 * t749) * t768 + t600;
t606 = t748 * t612 - t618 * t743;
t605 = m(4) * t674 - mrSges(4,2) * t730 + mrSges(4,3) * t722 - t718 * t731 + t720 * t764 + t606;
t586 = t597 * t773 - t599 * t735 + t605 * t774;
t582 = m(3) * t700 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t751 + t586;
t585 = t597 * t779 + t739 * t599 + t605 * t780;
t584 = m(3) * t716 + t585;
t591 = -t597 * t744 + t749 * t605;
t590 = m(3) * t701 - mrSges(3,1) * t751 - qJDD(2) * mrSges(3,2) + t591;
t572 = t582 * t771 - t584 * t736 + t590 * t772;
t570 = m(2) * t726 + t572;
t579 = -t582 * t745 + t750 * t590;
t578 = m(2) * t727 + t579;
t770 = t737 * t570 + t733 * t578;
t571 = t582 * t777 + t740 * t584 + t590 * t778;
t763 = -t570 * t733 + t737 * t578;
t648 = Ifges(7,5) * t677 + Ifges(7,6) * t676 + Ifges(7,3) * t689;
t650 = Ifges(7,1) * t677 + Ifges(7,4) * t676 + Ifges(7,5) * t689;
t623 = -mrSges(7,1) * t632 + mrSges(7,3) * t631 + Ifges(7,4) * t646 + Ifges(7,2) * t645 + Ifges(7,6) * t659 - t648 * t677 + t650 * t689;
t649 = Ifges(7,4) * t677 + Ifges(7,2) * t676 + Ifges(7,6) * t689;
t624 = mrSges(7,2) * t632 - mrSges(7,3) * t630 + Ifges(7,1) * t646 + Ifges(7,4) * t645 + Ifges(7,5) * t659 + t648 * t676 - t649 * t689;
t666 = Ifges(6,5) * t691 + Ifges(6,6) * t690 + Ifges(6,3) * t702;
t667 = Ifges(6,4) * t691 + Ifges(6,2) * t690 + Ifges(6,6) * t702;
t607 = mrSges(6,2) * t638 - mrSges(6,3) * t634 + Ifges(6,1) * t661 + Ifges(6,4) * t660 + Ifges(6,5) * t683 - pkin(13) * t622 - t623 * t741 + t624 * t746 + t666 * t690 - t667 * t702;
t668 = Ifges(6,1) * t691 + Ifges(6,4) * t690 + Ifges(6,5) * t702;
t608 = -mrSges(6,1) * t638 - mrSges(7,1) * t630 + mrSges(7,2) * t631 + mrSges(6,3) * t635 + Ifges(6,4) * t661 - Ifges(7,5) * t646 + Ifges(6,2) * t660 + Ifges(6,6) * t683 - Ifges(7,6) * t645 - Ifges(7,3) * t659 - pkin(5) * t622 - t649 * t677 + t650 * t676 - t666 * t691 + t668 * t702;
t681 = Ifges(5,4) * t704 + Ifges(5,2) * t703 + Ifges(5,6) * t714;
t682 = Ifges(5,1) * t704 + Ifges(5,4) * t703 + Ifges(5,5) * t714;
t592 = mrSges(5,1) * t642 - mrSges(5,2) * t643 + Ifges(5,5) * t685 + Ifges(5,6) * t684 + Ifges(5,3) * t705 + pkin(4) * t752 + pkin(12) * t762 + t742 * t607 + t747 * t608 + t704 * t681 - t703 * t682;
t709 = Ifges(4,6) * t731 + (Ifges(4,4) * t744 + Ifges(4,2) * t749) * t768;
t710 = Ifges(4,5) * t731 + (Ifges(4,1) * t744 + Ifges(4,4) * t749) * t768;
t680 = Ifges(5,5) * t704 + Ifges(5,6) * t703 + Ifges(5,3) * t714;
t593 = mrSges(5,2) * t647 - mrSges(5,3) * t642 + Ifges(5,1) * t685 + Ifges(5,4) * t684 + Ifges(5,5) * t705 - pkin(12) * t615 + t607 * t747 - t608 * t742 + t680 * t703 - t681 * t714;
t594 = Ifges(5,4) * t685 + Ifges(5,2) * t684 + Ifges(5,6) * t705 - t704 * t680 + t714 * t682 - mrSges(5,1) * t647 + mrSges(5,3) * t643 - Ifges(6,5) * t661 - Ifges(6,6) * t660 - Ifges(6,3) * t683 - t691 * t667 + t690 * t668 - mrSges(6,1) * t634 + mrSges(6,2) * t635 - t741 * t624 - t746 * t623 - pkin(5) * t753 - pkin(13) * t761 - pkin(4) * t615;
t754 = pkin(11) * t606 + t593 * t743 + t594 * t748;
t573 = mrSges(4,1) * t673 - mrSges(4,2) * t674 + Ifges(4,5) * t721 + Ifges(4,6) * t722 + Ifges(4,3) * t730 + pkin(3) * t601 + t738 * t592 + (t709 * t744 - t710 * t749) * t768 + t754 * t734;
t708 = Ifges(4,3) * t731 + (Ifges(4,5) * t744 + Ifges(4,6) * t749) * t768;
t574 = -mrSges(4,1) * t688 + mrSges(4,3) * t674 + Ifges(4,4) * t721 + Ifges(4,2) * t722 + Ifges(4,6) * t730 - pkin(3) * t600 - t734 * t592 - t708 * t765 + t731 * t710 + t738 * t754;
t575 = t708 * t764 + mrSges(4,2) * t688 - mrSges(4,3) * t673 + Ifges(4,1) * t721 + Ifges(4,4) * t722 + Ifges(4,5) * t730 + t748 * t593 - t743 * t594 - t731 * t709 + (-t600 * t734 - t601 * t738) * pkin(11);
t755 = pkin(10) * t591 + t574 * t749 + t575 * t744;
t567 = -mrSges(3,1) * t716 + mrSges(3,3) * t701 + t751 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t585 - t735 * t573 + t739 * t755;
t568 = mrSges(3,2) * t716 - mrSges(3,3) * t700 + Ifges(3,5) * qJDD(2) - t751 * Ifges(3,6) - t744 * t574 + t749 * t575 + (-t585 * t735 - t586 * t739) * pkin(10);
t756 = pkin(9) * t579 + t567 * t750 + t568 * t745;
t566 = mrSges(3,1) * t700 - mrSges(3,2) * t701 + Ifges(3,3) * qJDD(2) + pkin(2) * t586 + t739 * t573 + t735 * t755;
t565 = mrSges(2,2) * t732 - mrSges(2,3) * t726 - t745 * t567 + t750 * t568 + (-t571 * t736 - t572 * t740) * pkin(9);
t564 = -mrSges(2,1) * t732 + mrSges(2,3) * t727 - pkin(1) * t571 - t736 * t566 + t740 * t756;
t1 = [-m(1) * g(1) + t763; -m(1) * g(2) + t770; -m(1) * g(3) + m(2) * t732 + t571; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t770 - t733 * t564 + t737 * t565; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t763 + t737 * t564 + t733 * t565; -mrSges(1,1) * g(2) + mrSges(2,1) * t726 + mrSges(1,2) * g(1) - mrSges(2,2) * t727 + pkin(1) * t572 + t740 * t566 + t736 * t756;];
tauB  = t1;
