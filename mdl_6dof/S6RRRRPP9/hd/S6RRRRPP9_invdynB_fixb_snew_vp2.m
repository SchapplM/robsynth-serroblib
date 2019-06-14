% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:12:49
% EndTime: 2019-05-07 19:13:10
% DurationCPUTime: 10.81s
% Computational Cost: add. (164099->355), mult. (347985->433), div. (0->0), fcn. (268794->10), ass. (0->145)
t779 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t778 = -Ifges(6,1) - Ifges(5,3) - Ifges(7,1);
t760 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t759 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t777 = Ifges(5,2) + Ifges(7,2) + Ifges(6,3);
t758 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t728 = sin(pkin(6));
t732 = sin(qJ(2));
t735 = cos(qJ(2));
t761 = qJD(1) * qJD(2);
t716 = (-qJDD(1) * t735 + t732 * t761) * t728;
t776 = -2 * qJD(5);
t775 = cos(qJ(4));
t774 = pkin(8) * t728;
t729 = cos(pkin(6));
t773 = t729 * g(3);
t772 = -mrSges(7,1) - mrSges(5,3);
t725 = qJD(1) * t729 + qJD(2);
t731 = sin(qJ(3));
t734 = cos(qJ(3));
t763 = qJD(1) * t728;
t753 = t732 * t763;
t705 = t725 * t731 + t734 * t753;
t762 = qJD(1) * t735;
t752 = t728 * t762;
t720 = qJD(3) - t752;
t730 = sin(qJ(4));
t691 = t730 * t705 - t720 * t775;
t704 = t734 * t725 - t731 * t753;
t702 = qJD(4) - t704;
t771 = t691 * t702;
t770 = t728 * t732;
t769 = t728 * t735;
t768 = t729 * t732;
t767 = t729 * t735;
t733 = sin(qJ(1));
t736 = cos(qJ(1));
t721 = t733 * g(1) - g(2) * t736;
t737 = qJD(1) ^ 2;
t711 = qJDD(1) * pkin(1) + t737 * t774 + t721;
t722 = -g(1) * t736 - g(2) * t733;
t712 = -pkin(1) * t737 + qJDD(1) * t774 + t722;
t764 = t711 * t768 + t735 * t712;
t686 = -g(3) * t770 + t764;
t709 = mrSges(3,1) * t725 - mrSges(3,3) * t753;
t713 = (-mrSges(3,1) * t735 + mrSges(3,2) * t732) * t763;
t724 = qJDD(1) * t729 + qJDD(2);
t714 = (-pkin(2) * t735 - pkin(9) * t732) * t763;
t723 = t725 ^ 2;
t656 = -t723 * pkin(2) + t724 * pkin(9) + (-g(3) * t732 + t714 * t762) * t728 + t764;
t715 = (qJDD(1) * t732 + t735 * t761) * t728;
t657 = t716 * pkin(2) - t715 * pkin(9) - t773 + (-t711 + (pkin(2) * t732 - pkin(9) * t735) * t725 * qJD(1)) * t728;
t629 = t734 * t656 + t731 * t657;
t683 = -qJD(3) * t705 - t715 * t731 + t724 * t734;
t687 = -mrSges(4,1) * t704 + mrSges(4,2) * t705;
t694 = mrSges(4,1) * t720 - mrSges(4,3) * t705;
t708 = qJDD(3) + t716;
t688 = -pkin(3) * t704 - pkin(10) * t705;
t719 = t720 ^ 2;
t625 = -pkin(3) * t719 + pkin(10) * t708 + t688 * t704 + t629;
t685 = -g(3) * t769 + t711 * t767 - t732 * t712;
t655 = -t724 * pkin(2) - t723 * pkin(9) + t714 * t753 - t685;
t684 = qJD(3) * t704 + t715 * t734 + t724 * t731;
t627 = (-t704 * t720 - t684) * pkin(10) + (t705 * t720 - t683) * pkin(3) + t655;
t620 = -t730 * t625 + t627 * t775;
t636 = -t691 * qJD(4) + t684 * t775 + t730 * t708;
t692 = t705 * t775 + t730 * t720;
t661 = -mrSges(7,2) * t692 + mrSges(7,3) * t691;
t663 = mrSges(5,1) * t691 + mrSges(5,2) * t692;
t668 = mrSges(6,1) * t691 - mrSges(6,3) * t702;
t671 = -mrSges(5,2) * t702 - mrSges(5,3) * t691;
t681 = qJDD(4) - t683;
t662 = pkin(4) * t691 - qJ(5) * t692;
t701 = t702 ^ 2;
t618 = -t681 * pkin(4) - t701 * qJ(5) + t692 * t662 + qJDD(5) - t620;
t664 = -mrSges(6,2) * t691 - mrSges(6,3) * t692;
t612 = -0.2e1 * qJD(6) * t702 + (t691 * t692 - t681) * qJ(6) + (t636 + t771) * pkin(5) + t618;
t669 = -mrSges(7,1) * t691 + mrSges(7,2) * t702;
t748 = m(7) * t612 - t681 * mrSges(7,3) - t702 * t669;
t743 = -m(6) * t618 - t636 * mrSges(6,1) - t692 * t664 - t748;
t607 = m(5) * t620 + (-t668 + t671) * t702 + (-t661 - t663) * t692 + (mrSges(5,1) - mrSges(6,2)) * t681 + t772 * t636 + t743;
t621 = t775 * t625 + t730 * t627;
t635 = t692 * qJD(4) + t730 * t684 - t708 * t775;
t672 = mrSges(5,1) * t702 - mrSges(5,3) * t692;
t741 = -t701 * pkin(4) + t681 * qJ(5) - t691 * t662 + t621;
t617 = t702 * t776 - t741;
t670 = mrSges(6,1) * t692 + mrSges(6,2) * t702;
t666 = pkin(5) * t692 - qJ(6) * t702;
t690 = t691 ^ 2;
t614 = -t635 * pkin(5) - t690 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t666) * t702 + t741;
t667 = mrSges(7,1) * t692 - mrSges(7,3) * t702;
t757 = m(7) * t614 + t681 * mrSges(7,2) + t702 * t667;
t745 = -m(6) * t617 + t681 * mrSges(6,3) + t702 * t670 + t757;
t765 = -t661 - t664;
t609 = m(5) * t621 - t681 * mrSges(5,2) - t702 * t672 + (-t663 + t765) * t691 + (-mrSges(6,1) + t772) * t635 + t745;
t749 = -t607 * t730 + t775 * t609;
t603 = m(4) * t629 - mrSges(4,2) * t708 + mrSges(4,3) * t683 + t687 * t704 - t694 * t720 + t749;
t628 = -t731 * t656 + t734 * t657;
t693 = -mrSges(4,2) * t720 + mrSges(4,3) * t704;
t624 = -t708 * pkin(3) - t719 * pkin(10) + t705 * t688 - t628;
t739 = (-t636 + t771) * qJ(5) + t624 + (pkin(4) * t702 + t776) * t692;
t619 = t635 * pkin(4) + t739;
t616 = -t690 * pkin(5) + 0.2e1 * qJD(6) * t691 - t692 * t666 + (pkin(4) + qJ(6)) * t635 + t739;
t747 = m(7) * t616 - t636 * mrSges(7,2) + t635 * mrSges(7,3) - t692 * t667 + t691 * t669;
t742 = -m(6) * t619 + t635 * mrSges(6,2) + t691 * t668 - t747;
t738 = -m(5) * t624 - t635 * mrSges(5,1) - t691 * t671 + (t670 - t672) * t692 + (-mrSges(5,2) + mrSges(6,3)) * t636 + t742;
t606 = m(4) * t628 + t708 * mrSges(4,1) - t684 * mrSges(4,3) - t705 * t687 + t720 * t693 + t738;
t750 = t734 * t603 - t606 * t731;
t593 = m(3) * t686 - mrSges(3,2) * t724 - mrSges(3,3) * t716 - t709 * t725 + t713 * t752 + t750;
t596 = t731 * t603 + t734 * t606;
t698 = -t728 * t711 - t773;
t710 = -mrSges(3,2) * t725 + mrSges(3,3) * t752;
t595 = m(3) * t698 + t716 * mrSges(3,1) + t715 * mrSges(3,2) + (t709 * t732 - t710 * t735) * t763 + t596;
t604 = t607 * t775 + t730 * t609;
t740 = -m(4) * t655 + t683 * mrSges(4,1) - t684 * mrSges(4,2) + t704 * t693 - t705 * t694 - t604;
t600 = m(3) * t685 + t724 * mrSges(3,1) - t715 * mrSges(3,3) + t725 * t710 - t713 * t753 + t740;
t582 = t593 * t768 - t595 * t728 + t600 * t767;
t580 = m(2) * t721 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t737 + t582;
t588 = t735 * t593 - t600 * t732;
t587 = m(2) * t722 - mrSges(2,1) * t737 - qJDD(1) * mrSges(2,2) + t588;
t766 = t736 * t580 + t733 * t587;
t581 = t593 * t770 + t729 * t595 + t600 * t769;
t756 = t758 * t691 - t759 * t692 + t778 * t702;
t755 = t777 * t691 - t760 * t692 - t758 * t702;
t754 = -t760 * t691 + t779 * t692 + t759 * t702;
t751 = -t580 * t733 + t736 * t587;
t610 = -t636 * mrSges(6,3) - t692 * t670 - t742;
t589 = -mrSges(5,1) * t624 + mrSges(5,3) * t621 - mrSges(6,1) * t617 + mrSges(6,2) * t619 + mrSges(7,1) * t614 - mrSges(7,3) * t616 - pkin(5) * (t691 * t661 - t757) - qJ(6) * t747 - pkin(4) * t610 + t754 * t702 + t756 * t692 + t758 * t681 + t760 * t636 + (-mrSges(7,1) * pkin(5) - t777) * t635;
t611 = t636 * mrSges(7,1) + t692 * t661 + t748;
t597 = mrSges(6,1) * t618 + mrSges(7,1) * t612 + mrSges(5,2) * t624 - mrSges(7,2) * t616 - mrSges(5,3) * t620 - mrSges(6,3) * t619 + pkin(5) * t611 - qJ(5) * t610 - t760 * t635 + t779 * t636 + t759 * t681 + t756 * t691 + t755 * t702;
t677 = Ifges(4,5) * t705 + Ifges(4,6) * t704 + Ifges(4,3) * t720;
t678 = Ifges(4,4) * t705 + Ifges(4,2) * t704 + Ifges(4,6) * t720;
t583 = mrSges(4,2) * t655 - mrSges(4,3) * t628 + Ifges(4,1) * t684 + Ifges(4,4) * t683 + Ifges(4,5) * t708 - pkin(10) * t604 - t730 * t589 + t597 * t775 + t704 * t677 - t720 * t678;
t679 = Ifges(4,1) * t705 + Ifges(4,4) * t704 + Ifges(4,5) * t720;
t584 = -qJ(5) * t745 + t720 * t679 + Ifges(4,6) * t708 - pkin(4) * (-t702 * t668 + t743) - t705 * t677 + Ifges(4,4) * t684 + Ifges(4,2) * t683 - mrSges(4,1) * t655 + mrSges(4,3) * t629 + mrSges(5,2) * t621 - mrSges(5,1) * t620 + mrSges(6,3) * t617 - mrSges(6,2) * t618 - mrSges(7,2) * t614 + qJ(6) * t611 + mrSges(7,3) * t612 - pkin(3) * t604 + (pkin(4) * t661 + t755) * t692 + (mrSges(6,2) * pkin(4) + t778) * t681 + (mrSges(7,1) * pkin(4) - t759) * t636 + (-qJ(5) * t765 - t754) * t691 + (-qJ(5) * (-mrSges(6,1) - mrSges(7,1)) + t758) * t635;
t695 = Ifges(3,3) * t725 + (Ifges(3,5) * t732 + Ifges(3,6) * t735) * t763;
t696 = Ifges(3,6) * t725 + (Ifges(3,4) * t732 + Ifges(3,2) * t735) * t763;
t577 = mrSges(3,2) * t698 - mrSges(3,3) * t685 + Ifges(3,1) * t715 - Ifges(3,4) * t716 + Ifges(3,5) * t724 - pkin(9) * t596 + t583 * t734 - t584 * t731 + t695 * t752 - t696 * t725;
t697 = Ifges(3,5) * t725 + (Ifges(3,1) * t732 + Ifges(3,4) * t735) * t763;
t578 = Ifges(3,4) * t715 - Ifges(3,2) * t716 + Ifges(3,6) * t724 - t695 * t753 + t725 * t697 - mrSges(3,1) * t698 + mrSges(3,3) * t686 - Ifges(4,5) * t684 - Ifges(4,6) * t683 - Ifges(4,3) * t708 - t705 * t678 + t704 * t679 - mrSges(4,1) * t628 + mrSges(4,2) * t629 - t730 * t597 - t775 * t589 - pkin(3) * t738 - pkin(10) * t749 - pkin(2) * t596;
t744 = pkin(8) * t588 + t577 * t732 + t578 * t735;
t576 = Ifges(3,5) * t715 - Ifges(3,6) * t716 + Ifges(3,3) * t724 + mrSges(3,1) * t685 - mrSges(3,2) * t686 + t731 * t583 + t734 * t584 + pkin(2) * t740 + pkin(9) * t750 + (t696 * t732 - t697 * t735) * t763;
t575 = -mrSges(2,2) * g(3) - mrSges(2,3) * t721 + Ifges(2,5) * qJDD(1) - t737 * Ifges(2,6) + t735 * t577 - t732 * t578 + (-t581 * t728 - t582 * t729) * pkin(8);
t574 = mrSges(2,1) * g(3) + mrSges(2,3) * t722 + t737 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t581 - t728 * t576 + t729 * t744;
t1 = [-m(1) * g(1) + t751; -m(1) * g(2) + t766; (-m(1) - m(2)) * g(3) + t581; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t766 - t733 * t574 + t736 * t575; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t751 + t736 * t574 + t733 * t575; -mrSges(1,1) * g(2) + mrSges(2,1) * t721 + mrSges(1,2) * g(1) - mrSges(2,2) * t722 + Ifges(2,3) * qJDD(1) + pkin(1) * t582 + t729 * t576 + t728 * t744;];
tauB  = t1;
