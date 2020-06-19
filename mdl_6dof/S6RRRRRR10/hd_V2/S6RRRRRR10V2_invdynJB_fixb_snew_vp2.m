% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% m [7x1]
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 10:42
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR10V2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 10:18:07
% EndTime: 2020-06-19 10:18:37
% DurationCPUTime: 29.69s
% Computational Cost: add. (208247->355), mult. (410941->449), div. (0->0), fcn. (313839->12), ass. (0->140)
t743 = sin(qJ(3));
t744 = sin(qJ(2));
t749 = cos(qJ(3));
t750 = cos(qJ(2));
t719 = (t743 * t744 - t749 * t750) * qJD(1);
t745 = sin(qJ(1));
t751 = cos(qJ(1));
t733 = -t751 * g(1) - t745 * g(2);
t752 = qJD(1) ^ 2;
t726 = -t752 * pkin(1) + t733;
t713 = -t750 * g(3) - t744 * t726;
t706 = (t744 * t750 * t752 + qJDD(2)) * pkin(2) + t713;
t714 = -t744 * g(3) + t750 * t726;
t707 = (-t750 ^ 2 * t752 - qJD(2) ^ 2) * pkin(2) + t714;
t683 = t743 * t706 + t749 * t707;
t720 = (t743 * t750 + t744 * t749) * qJD(1);
t769 = qJD(1) * qJD(2);
t727 = t744 * qJDD(1) + t750 * t769;
t768 = t744 * t769;
t728 = t750 * qJDD(1) - t768;
t693 = -t720 * qJD(3) - t743 * t727 + t749 * t728;
t701 = t719 * mrSges(4,1) + t720 * mrSges(4,2);
t738 = qJD(2) + qJD(3);
t712 = t738 * mrSges(4,1) - t720 * mrSges(4,3);
t737 = qJDD(2) + qJDD(3);
t703 = t719 * pkin(3) - t720 * pkin(5);
t736 = t738 ^ 2;
t667 = -t736 * pkin(3) + t737 * pkin(5) - t719 * t703 + t683;
t742 = sin(qJ(4));
t748 = cos(qJ(4));
t694 = -t719 * qJD(3) + t749 * t727 + t743 * t728;
t732 = t745 * g(1) - t751 * g(2);
t724 = -qJDD(1) * pkin(1) - t732;
t704 = (-t728 + t768) * pkin(2) + t724;
t755 = (t719 * t738 - t694) * pkin(5) + (t720 * t738 - t693) * pkin(3) + t704;
t651 = t748 * t667 + t742 * t755;
t682 = t749 * t706 - t743 * t707;
t666 = -t737 * pkin(3) - t736 * pkin(5) + t720 * t703 - t682;
t741 = sin(qJ(5));
t747 = cos(qJ(5));
t646 = t747 * t651 + t741 * t666;
t710 = t748 * t720 + t742 * t738;
t671 = -t710 * qJD(4) - t742 * t694 + t748 * t737;
t670 = qJDD(5) - t671;
t715 = qJD(4) + t719;
t687 = -t741 * t710 + t747 * t715;
t688 = t747 * t710 + t741 * t715;
t642 = (-t687 * t688 + t670) * pkin(6) + t646;
t650 = t742 * t667 - t748 * t755;
t709 = -t742 * t720 + t748 * t738;
t672 = t709 * qJD(4) + t748 * t694 + t742 * t737;
t692 = qJDD(4) - t693;
t657 = t687 * qJD(5) + t747 * t672 + t741 * t692;
t708 = qJD(5) - t709;
t643 = (-t687 * t708 - t657) * pkin(6) + t650;
t740 = sin(qJ(6));
t746 = cos(qJ(6));
t640 = -t740 * t642 + t746 * t643;
t673 = -t740 * t688 + t746 * t708;
t648 = t673 * qJD(6) + t746 * t657 + t740 * t670;
t656 = -t688 * qJD(5) - t741 * t672 + t747 * t692;
t655 = qJDD(6) - t656;
t674 = t746 * t688 + t740 * t708;
t659 = -t673 * mrSges(7,1) + t674 * mrSges(7,2);
t686 = qJD(6) - t687;
t660 = -t686 * mrSges(7,2) + t673 * mrSges(7,3);
t638 = m(7) * t640 + t655 * mrSges(7,1) - t648 * mrSges(7,3) - t674 * t659 + t686 * t660;
t641 = t746 * t642 + t740 * t643;
t647 = -t674 * qJD(6) - t740 * t657 + t746 * t670;
t661 = t686 * mrSges(7,1) - t674 * mrSges(7,3);
t639 = m(7) * t641 - t655 * mrSges(7,2) + t647 * mrSges(7,3) + t673 * t659 - t686 * t661;
t633 = -t740 * t638 + t746 * t639;
t668 = -t687 * mrSges(6,1) + t688 * mrSges(6,2);
t676 = t708 * mrSges(6,1) - t688 * mrSges(6,3);
t632 = m(6) * t646 - t670 * mrSges(6,2) + t656 * mrSges(6,3) + t687 * t668 - t708 * t676 + t633;
t645 = -t741 * t651 + t747 * t666;
t644 = (-t688 ^ 2 - t708 ^ 2) * pkin(6) - t645;
t675 = -t708 * mrSges(6,2) + t687 * mrSges(6,3);
t636 = m(6) * t645 - m(7) * t644 + t670 * mrSges(6,1) + t647 * mrSges(7,1) - t648 * mrSges(7,2) - t657 * mrSges(6,3) + t673 * t660 - t674 * t661 - t688 * t668 + t708 * t675;
t684 = -t709 * mrSges(5,1) + t710 * mrSges(5,2);
t696 = t715 * mrSges(5,1) - t710 * mrSges(5,3);
t626 = m(5) * t651 - t692 * mrSges(5,2) + t671 * mrSges(5,3) + t747 * t632 - t741 * t636 + t709 * t684 - t715 * t696;
t695 = -t715 * mrSges(5,2) + t709 * mrSges(5,3);
t765 = t746 * t638 + t740 * t639;
t630 = t692 * mrSges(5,1) + t656 * mrSges(6,1) - t657 * mrSges(6,2) - t672 * mrSges(5,3) + t687 * t675 - t688 * t676 - t710 * t684 + t715 * t695 + (-m(5) - m(6)) * t650 - t765;
t766 = t748 * t626 - t742 * t630;
t615 = m(4) * t683 - t737 * mrSges(4,2) + t693 * mrSges(4,3) - t719 * t701 - t738 * t712 + t766;
t711 = -t738 * mrSges(4,2) - t719 * mrSges(4,3);
t759 = -m(5) * t666 + t671 * mrSges(5,1) - t672 * mrSges(5,2) - t741 * t632 - t747 * t636 + t709 * t695 - t710 * t696;
t623 = m(4) * t682 + t737 * mrSges(4,1) - t694 * mrSges(4,3) - t720 * t701 + t738 * t711 + t759;
t608 = t743 * t615 + t749 * t623;
t717 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t744 + Ifges(3,2) * t750) * qJD(1);
t718 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t744 + Ifges(3,4) * t750) * qJD(1);
t652 = Ifges(7,5) * t674 + Ifges(7,6) * t673 + Ifges(7,3) * t686;
t654 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t686;
t634 = -mrSges(7,1) * t644 + mrSges(7,3) * t641 + Ifges(7,4) * t648 + Ifges(7,2) * t647 + Ifges(7,6) * t655 - t674 * t652 + t686 * t654;
t653 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t686;
t635 = mrSges(7,2) * t644 - mrSges(7,3) * t640 + Ifges(7,1) * t648 + Ifges(7,4) * t647 + Ifges(7,5) * t655 + t673 * t652 - t686 * t653;
t662 = Ifges(6,5) * t688 + Ifges(6,6) * t687 + Ifges(6,3) * t708;
t663 = Ifges(6,4) * t688 + Ifges(6,2) * t687 + Ifges(6,6) * t708;
t621 = mrSges(6,2) * t650 - mrSges(6,3) * t645 + Ifges(6,1) * t657 + Ifges(6,4) * t656 + Ifges(6,5) * t670 - pkin(6) * t765 - t740 * t634 + t746 * t635 + t687 * t662 - t708 * t663;
t664 = Ifges(6,1) * t688 + Ifges(6,4) * t687 + Ifges(6,5) * t708;
t758 = mrSges(7,1) * t640 - mrSges(7,2) * t641 + Ifges(7,5) * t648 + Ifges(7,6) * t647 + Ifges(7,3) * t655 + t674 * t653 - t673 * t654;
t631 = -mrSges(6,1) * t650 + mrSges(6,3) * t646 + Ifges(6,4) * t657 + Ifges(6,2) * t656 + Ifges(6,6) * t670 - t688 * t662 + t708 * t664 - t758;
t677 = Ifges(5,5) * t710 + Ifges(5,6) * t709 + Ifges(5,3) * t715;
t678 = Ifges(5,4) * t710 + Ifges(5,2) * t709 + Ifges(5,6) * t715;
t610 = mrSges(5,2) * t666 + mrSges(5,3) * t650 + Ifges(5,1) * t672 + Ifges(5,4) * t671 + Ifges(5,5) * t692 + t747 * t621 - t741 * t631 + t709 * t677 - t715 * t678;
t679 = Ifges(5,1) * t710 + Ifges(5,4) * t709 + Ifges(5,5) * t715;
t754 = mrSges(6,1) * t645 - mrSges(6,2) * t646 + Ifges(6,5) * t657 + Ifges(6,6) * t656 + Ifges(6,3) * t670 + pkin(6) * t633 + t746 * t634 + t740 * t635 + t688 * t663 - t687 * t664;
t620 = -mrSges(5,1) * t666 + mrSges(5,3) * t651 + Ifges(5,4) * t672 + Ifges(5,2) * t671 + Ifges(5,6) * t692 - t710 * t677 + t715 * t679 - t754;
t698 = Ifges(4,4) * t720 - Ifges(4,2) * t719 + Ifges(4,6) * t738;
t699 = Ifges(4,1) * t720 - Ifges(4,4) * t719 + Ifges(4,5) * t738;
t760 = -mrSges(4,1) * t682 + mrSges(4,2) * t683 - Ifges(4,5) * t694 - Ifges(4,6) * t693 - Ifges(4,3) * t737 - pkin(3) * t759 - pkin(5) * t766 - t742 * t610 - t748 * t620 - t720 * t698 - t719 * t699;
t774 = mrSges(3,1) * t713 - mrSges(3,2) * t714 + Ifges(3,5) * t727 + Ifges(3,6) * t728 + Ifges(3,3) * qJDD(2) + pkin(2) * t608 + (t744 * t717 - t750 * t718) * qJD(1) - t760;
t725 = (-mrSges(3,1) * t750 + mrSges(3,2) * t744) * qJD(1);
t770 = qJD(1) * t750;
t731 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t770;
t771 = qJD(1) * t744;
t606 = m(3) * t713 + qJDD(2) * mrSges(3,1) - t727 * mrSges(3,3) + qJD(2) * t731 - t725 * t771 + t608;
t730 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t771;
t607 = m(3) * t714 - qJDD(2) * mrSges(3,2) + t728 * mrSges(3,3) - qJD(2) * t730 + t749 * t615 - t743 * t623 + t725 * t770;
t601 = m(2) * t733 - t752 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t744 * t606 + t750 * t607;
t617 = t742 * t626 + t748 * t630;
t761 = m(4) * t704 - t693 * mrSges(4,1) + t694 * mrSges(4,2) + t719 * t711 + t720 * t712 + t617;
t757 = -m(3) * t724 + t728 * mrSges(3,1) - t727 * mrSges(3,2) - t730 * t771 + t731 * t770 - t761;
t612 = m(2) * t732 + qJDD(1) * mrSges(2,1) - t752 * mrSges(2,2) + t757;
t773 = t745 * t601 + t751 * t612;
t772 = t750 * t606 + t744 * t607;
t767 = t751 * t601 - t745 * t612;
t697 = Ifges(4,5) * t720 - Ifges(4,6) * t719 + Ifges(4,3) * t738;
t602 = mrSges(4,2) * t704 - mrSges(4,3) * t682 + Ifges(4,1) * t694 + Ifges(4,4) * t693 + Ifges(4,5) * t737 - pkin(5) * t617 + t748 * t610 - t742 * t620 - t719 * t697 - t738 * t698;
t756 = mrSges(5,1) * t650 + mrSges(5,2) * t651 - Ifges(5,5) * t672 - Ifges(5,6) * t671 - Ifges(5,3) * t692 - t741 * t621 - t747 * t631 - t710 * t678 + t709 * t679;
t603 = -mrSges(4,1) * t704 + mrSges(4,3) * t683 + Ifges(4,4) * t694 + Ifges(4,2) * t693 + Ifges(4,6) * t737 - pkin(3) * t617 - t720 * t697 + t738 * t699 + t756;
t716 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t744 + Ifges(3,6) * t750) * qJD(1);
t596 = -mrSges(3,1) * t724 + mrSges(3,3) * t714 + Ifges(3,4) * t727 + Ifges(3,2) * t728 + Ifges(3,6) * qJDD(2) - pkin(2) * t761 + qJD(2) * t718 + t743 * t602 + t749 * t603 - t716 * t771;
t598 = mrSges(3,2) * t724 - mrSges(3,3) * t713 + Ifges(3,1) * t727 + Ifges(3,4) * t728 + Ifges(3,5) * qJDD(2) - qJD(2) * t717 + t749 * t602 - t743 * t603 + t716 * t770;
t762 = mrSges(2,1) * t732 - mrSges(2,2) * t733 + Ifges(2,3) * qJDD(1) + pkin(1) * t757 + t750 * t596 + t744 * t598;
t594 = mrSges(2,1) * g(3) + mrSges(2,3) * t733 + t752 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t772 - t774;
t593 = -mrSges(2,2) * g(3) - mrSges(2,3) * t732 + Ifges(2,5) * qJDD(1) - t752 * Ifges(2,6) - t744 * t596 + t750 * t598;
t1 = [-m(1) * g(1) + t767; -m(1) * g(2) + t773; (-m(1) - m(2)) * g(3) + t772; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t773 + t751 * t593 - t745 * t594; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t767 + t745 * t593 + t751 * t594; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t762; t762; t774; -t760; -t756; t754; t758;];
tauJB = t1;
