% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:24:22
% EndTime: 2019-05-07 18:24:36
% DurationCPUTime: 7.05s
% Computational Cost: add. (79868->341), mult. (158815->402), div. (0->0), fcn. (109124->8), ass. (0->128)
t774 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t759 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t773 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t772 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t771 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t757 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t730 = sin(qJ(3));
t733 = cos(qJ(3));
t731 = sin(qJ(2));
t762 = qJD(1) * t731;
t712 = qJD(2) * t730 + t733 * t762;
t734 = cos(qJ(2));
t760 = qJD(1) * qJD(2);
t752 = t734 * t760;
t715 = qJDD(1) * t731 + t752;
t683 = -qJD(3) * t712 + qJDD(2) * t733 - t715 * t730;
t711 = qJD(2) * t733 - t730 * t762;
t684 = qJD(3) * t711 + qJDD(2) * t730 + t715 * t733;
t729 = sin(qJ(4));
t768 = cos(qJ(4));
t686 = -t768 * t711 + t729 * t712;
t636 = -t686 * qJD(4) + t729 * t683 + t768 * t684;
t761 = qJD(1) * t734;
t725 = qJD(3) - t761;
t693 = pkin(3) * t725 - pkin(9) * t712;
t709 = t711 ^ 2;
t732 = sin(qJ(1));
t735 = cos(qJ(1));
t721 = -g(1) * t735 - g(2) * t732;
t737 = qJD(1) ^ 2;
t705 = -pkin(1) * t737 + qJDD(1) * pkin(7) + t721;
t691 = -t734 * g(3) - t731 * t705;
t714 = (-pkin(2) * t734 - pkin(8) * t731) * qJD(1);
t736 = qJD(2) ^ 2;
t743 = qJDD(2) * pkin(2) + t736 * pkin(8) - t714 * t762 + t691;
t741 = t683 * pkin(3) + t709 * pkin(9) - t712 * t693 + t743;
t723 = -qJD(4) - t725;
t766 = t686 * t723;
t770 = -(t636 + t766) * qJ(5) - t741;
t769 = -2 * qJD(5);
t767 = -mrSges(5,3) - mrSges(6,2);
t670 = -mrSges(7,2) * t723 + mrSges(7,3) * t686;
t765 = t723 * t670;
t692 = -g(3) * t731 + t734 * t705;
t713 = (-mrSges(3,1) * t734 + mrSges(3,2) * t731) * qJD(1);
t726 = t731 * t760;
t716 = qJDD(1) * t734 - t726;
t718 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t762;
t720 = t732 * g(1) - t735 * g(2);
t704 = -qJDD(1) * pkin(1) - t737 * pkin(7) - t720;
t663 = (-t715 - t752) * pkin(8) + (-t716 + t726) * pkin(2) + t704;
t668 = -pkin(2) * t736 + qJDD(2) * pkin(8) + t714 * t761 + t692;
t637 = t733 * t663 - t730 * t668;
t710 = qJDD(3) - t716;
t623 = (t711 * t725 - t684) * pkin(9) + (t711 * t712 + t710) * pkin(3) + t637;
t638 = t730 * t663 + t733 * t668;
t626 = -pkin(3) * t709 + pkin(9) * t683 - t693 * t725 + t638;
t620 = t768 * t623 - t729 * t626;
t671 = mrSges(5,2) * t723 - mrSges(5,3) * t686;
t687 = t729 * t711 + t768 * t712;
t706 = -qJDD(4) - t710;
t657 = pkin(4) * t686 - qJ(5) * t687;
t722 = t723 ^ 2;
t617 = t706 * pkin(4) - t722 * qJ(5) + t687 * t657 + qJDD(5) - t620;
t669 = -mrSges(6,2) * t686 - mrSges(6,3) * t723;
t610 = -0.2e1 * qJD(6) * t687 + (-t636 + t766) * qJ(6) + (t686 * t687 + t706) * pkin(5) + t617;
t659 = -mrSges(7,1) * t686 + mrSges(7,2) * t687;
t746 = -m(7) * t610 + t636 * mrSges(7,3) + t687 * t659;
t742 = -m(6) * t617 - t706 * mrSges(6,1) - t723 * t669 + t746;
t658 = mrSges(6,1) * t686 - mrSges(6,3) * t687;
t763 = -mrSges(5,1) * t686 - mrSges(5,2) * t687 - t658;
t604 = m(5) * t620 + (-t670 - t671) * t723 + (-mrSges(5,1) - mrSges(7,1)) * t706 + t763 * t687 + t767 * t636 + t742;
t621 = t729 * t623 + t768 * t626;
t635 = t687 * qJD(4) - t768 * t683 + t729 * t684;
t673 = mrSges(7,1) * t723 - mrSges(7,3) * t687;
t674 = -mrSges(5,1) * t723 - mrSges(5,3) * t687;
t616 = -pkin(4) * t722 - t706 * qJ(5) - t686 * t657 + t723 * t769 + t621;
t675 = mrSges(6,1) * t723 + mrSges(6,2) * t687;
t672 = pkin(5) * t723 - qJ(6) * t687;
t685 = t686 ^ 2;
t612 = -pkin(5) * t685 + qJ(6) * t635 + 0.2e1 * qJD(6) * t686 - t672 * t723 + t616;
t756 = m(7) * t612 + t635 * mrSges(7,3) + t686 * t659;
t744 = m(6) * t616 - t706 * mrSges(6,3) - t723 * t675 + t756;
t607 = m(5) * t621 + (-t673 + t674) * t723 + (mrSges(5,2) - mrSges(7,2)) * t706 + t763 * t686 + t767 * t635 + t744;
t600 = t768 * t604 + t729 * t607;
t688 = -mrSges(4,1) * t711 + mrSges(4,2) * t712;
t689 = -mrSges(4,2) * t725 + mrSges(4,3) * t711;
t598 = m(4) * t637 + mrSges(4,1) * t710 - mrSges(4,3) * t684 - t688 * t712 + t689 * t725 + t600;
t690 = mrSges(4,1) * t725 - mrSges(4,3) * t712;
t748 = -t604 * t729 + t768 * t607;
t599 = m(4) * t638 - mrSges(4,2) * t710 + mrSges(4,3) * t683 + t688 * t711 - t690 * t725 + t748;
t749 = -t598 * t730 + t733 * t599;
t593 = m(3) * t692 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t716 - qJD(2) * t718 + t713 * t761 + t749;
t719 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t761;
t619 = t687 * t769 + (-t687 * t723 + t635) * pkin(4) + t770;
t614 = -t685 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t635 + (pkin(4) * t723 + (2 * qJD(5)) + t672) * t687 - t770;
t745 = -m(7) * t614 + t635 * mrSges(7,1) - t636 * mrSges(7,2) + t686 * t670 - t687 * t673;
t608 = m(6) * t619 + t635 * mrSges(6,1) - t636 * mrSges(6,3) + t686 * t669 - t687 * t675 + t745;
t739 = -m(5) * t741 + t635 * mrSges(5,1) + t636 * mrSges(5,2) + t686 * t671 + t687 * t674 + t608;
t738 = m(4) * t743 + t683 * mrSges(4,1) - t684 * mrSges(4,2) + t711 * t689 - t712 * t690 - t739;
t602 = m(3) * t691 + qJDD(2) * mrSges(3,1) - t715 * mrSges(3,3) + qJD(2) * t719 - t713 * t762 + t738;
t750 = t734 * t593 - t602 * t731;
t587 = m(2) * t721 - mrSges(2,1) * t737 - qJDD(1) * mrSges(2,2) + t750;
t594 = t598 * t733 + t599 * t730;
t740 = -m(3) * t704 + t716 * mrSges(3,1) - mrSges(3,2) * t715 - t718 * t762 + t719 * t761 - t594;
t590 = m(2) * t720 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t737 + t740;
t764 = t732 * t587 + t735 * t590;
t588 = t731 * t593 + t734 * t602;
t755 = t686 * t757 - t687 * t773 + t723 * t771;
t754 = t686 * t772 + t687 * t759 - t723 * t757;
t753 = t759 * t686 - t687 * t774 + t773 * t723;
t751 = t735 * t587 - t590 * t732;
t703 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t731 + Ifges(3,4) * t734) * qJD(1);
t702 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t731 + Ifges(3,2) * t734) * qJD(1);
t701 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t731 + Ifges(3,6) * t734) * qJD(1);
t679 = Ifges(4,1) * t712 + Ifges(4,4) * t711 + Ifges(4,5) * t725;
t678 = Ifges(4,4) * t712 + Ifges(4,2) * t711 + Ifges(4,6) * t725;
t677 = Ifges(4,5) * t712 + Ifges(4,6) * t711 + Ifges(4,3) * t725;
t609 = t706 * mrSges(7,1) - t746 + t765;
t596 = -mrSges(5,2) * t741 + mrSges(6,2) * t617 + mrSges(7,2) * t614 - mrSges(5,3) * t620 - mrSges(6,3) * t619 - mrSges(7,3) * t610 - qJ(5) * t608 - qJ(6) * t609 - t759 * t635 + t636 * t774 + t755 * t686 - t706 * t773 + t754 * t723;
t595 = mrSges(5,1) * t741 + mrSges(5,3) * t621 - mrSges(6,1) * t619 + mrSges(6,2) * t616 + mrSges(7,1) * t614 - mrSges(7,3) * t612 - pkin(5) * t745 - qJ(6) * t756 - pkin(4) * t608 + (qJ(6) * t673 + t753) * t723 + (mrSges(7,2) * qJ(6) - t757) * t706 + t755 * t687 + t759 * t636 + t772 * t635;
t584 = -mrSges(4,2) * t743 - mrSges(4,3) * t637 + Ifges(4,1) * t684 + Ifges(4,4) * t683 + Ifges(4,5) * t710 - pkin(9) * t600 - t729 * t595 + t768 * t596 + t711 * t677 - t725 * t678;
t583 = mrSges(4,1) * t743 + mrSges(4,3) * t638 + Ifges(4,4) * t684 + Ifges(4,2) * t683 + Ifges(4,6) * t710 - pkin(3) * t739 + pkin(9) * t748 + t768 * t595 + t729 * t596 - t712 * t677 + t725 * t679;
t582 = -pkin(4) * (t742 - t765) + (mrSges(6,2) * pkin(4) - t773) * t636 + Ifges(3,6) * qJDD(2) + (mrSges(7,1) * pkin(4) + mrSges(7,2) * qJ(5) + t771) * t706 - t701 * t762 + (qJ(5) * t658 + t753) * t686 + (pkin(4) * t658 - t754) * t687 + (mrSges(6,2) * qJ(5) + t757) * t635 + Ifges(3,2) * t716 - Ifges(4,3) * t710 + t711 * t679 - t712 * t678 + Ifges(3,4) * t715 + qJD(2) * t703 - mrSges(3,1) * t704 + mrSges(3,3) * t692 - Ifges(4,6) * t683 - Ifges(4,5) * t684 - mrSges(4,1) * t637 + mrSges(4,2) * t638 - mrSges(5,1) * t620 + mrSges(5,2) * t621 - mrSges(6,3) * t616 + mrSges(6,1) * t617 - mrSges(7,2) * t612 + pkin(5) * t609 + mrSges(7,1) * t610 - pkin(3) * t600 - pkin(2) * t594 - qJ(5) * (-t723 * t673 + t744);
t581 = mrSges(3,2) * t704 - mrSges(3,3) * t691 + Ifges(3,1) * t715 + Ifges(3,4) * t716 + Ifges(3,5) * qJDD(2) - pkin(8) * t594 - qJD(2) * t702 - t583 * t730 + t584 * t733 + t701 * t761;
t580 = Ifges(2,6) * qJDD(1) + t737 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t721 - Ifges(3,5) * t715 - Ifges(3,6) * t716 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t691 + mrSges(3,2) * t692 - t730 * t584 - t733 * t583 - pkin(2) * t738 - pkin(8) * t749 - pkin(1) * t588 + (-t702 * t731 + t703 * t734) * qJD(1);
t579 = -mrSges(2,2) * g(3) - mrSges(2,3) * t720 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t737 - pkin(7) * t588 + t581 * t734 - t582 * t731;
t1 = [-m(1) * g(1) + t751; -m(1) * g(2) + t764; (-m(1) - m(2)) * g(3) + t588; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t764 + t735 * t579 - t732 * t580; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t751 + t732 * t579 + t735 * t580; -mrSges(1,1) * g(2) + mrSges(2,1) * t720 + mrSges(1,2) * g(1) - mrSges(2,2) * t721 + Ifges(2,3) * qJDD(1) + pkin(1) * t740 + pkin(7) * t750 + t731 * t581 + t734 * t582;];
tauB  = t1;
