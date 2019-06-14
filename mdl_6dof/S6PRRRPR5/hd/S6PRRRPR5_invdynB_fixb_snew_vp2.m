% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:01:35
% EndTime: 2019-05-05 08:02:20
% DurationCPUTime: 44.16s
% Computational Cost: add. (759680->353), mult. (1607939->468), div. (0->0), fcn. (1300853->16), ass. (0->156)
t695 = sin(pkin(7));
t703 = sin(qJ(3));
t707 = cos(qJ(3));
t724 = qJD(2) * qJD(3);
t679 = (-qJDD(2) * t707 + t703 * t724) * t695;
t694 = sin(pkin(12));
t698 = cos(pkin(12));
t685 = g(1) * t694 - g(2) * t698;
t686 = -g(1) * t698 - g(2) * t694;
t692 = -g(3) + qJDD(1);
t704 = sin(qJ(2));
t700 = cos(pkin(6));
t708 = cos(qJ(2));
t727 = t700 * t708;
t696 = sin(pkin(6));
t730 = t696 * t708;
t658 = t685 * t727 - t686 * t704 + t692 * t730;
t709 = qJD(2) ^ 2;
t734 = pkin(9) * t695;
t655 = qJDD(2) * pkin(2) + t709 * t734 + t658;
t728 = t700 * t704;
t731 = t696 * t704;
t659 = t685 * t728 + t708 * t686 + t692 * t731;
t656 = -pkin(2) * t709 + qJDD(2) * t734 + t659;
t672 = -t685 * t696 + t692 * t700;
t699 = cos(pkin(7));
t621 = -t703 * t656 + (t655 * t699 + t672 * t695) * t707;
t735 = 2 * qJD(5);
t691 = qJD(2) * t699 + qJD(3);
t725 = qJD(2) * t695;
t722 = t707 * t725;
t675 = -mrSges(4,2) * t691 + mrSges(4,3) * t722;
t676 = (-mrSges(4,1) * t707 + mrSges(4,2) * t703) * t725;
t678 = (qJDD(2) * t703 + t707 * t724) * t695;
t690 = qJDD(2) * t699 + qJDD(3);
t677 = (-pkin(3) * t707 - pkin(10) * t703) * t725;
t689 = t691 ^ 2;
t723 = t703 * t725;
t615 = -t690 * pkin(3) - t689 * pkin(10) + t677 * t723 - t621;
t702 = sin(qJ(4));
t706 = cos(qJ(4));
t671 = t691 * t702 + t706 * t723;
t646 = -qJD(4) * t671 - t678 * t702 + t690 * t706;
t670 = t691 * t706 - t702 * t723;
t647 = qJD(4) * t670 + t678 * t706 + t690 * t702;
t684 = qJD(4) - t722;
t660 = -mrSges(5,2) * t684 + mrSges(5,3) * t670;
t662 = mrSges(5,1) * t684 - mrSges(5,3) * t671;
t729 = t699 * t703;
t732 = t695 * t703;
t622 = t655 * t729 + t707 * t656 + t672 * t732;
t616 = -pkin(3) * t689 + pkin(10) * t690 + t677 * t722 + t622;
t668 = t699 * t672;
t619 = t679 * pkin(3) - t678 * pkin(10) + t668 + (-t655 + (pkin(3) * t703 - pkin(10) * t707) * t691 * qJD(2)) * t695;
t604 = -t702 * t616 + t706 * t619;
t673 = qJDD(4) + t679;
t601 = (t670 * t684 - t647) * qJ(5) + (t670 * t671 + t673) * pkin(4) + t604;
t605 = t706 * t616 + t702 * t619;
t661 = pkin(4) * t684 - qJ(5) * t671;
t669 = t670 ^ 2;
t603 = -pkin(4) * t669 + qJ(5) * t646 - t661 * t684 + t605;
t693 = sin(pkin(13));
t697 = cos(pkin(13));
t653 = t670 * t697 - t671 * t693;
t598 = t693 * t601 + t697 * t603 + t653 * t735;
t654 = t670 * t693 + t671 * t697;
t634 = -pkin(5) * t653 - pkin(11) * t654;
t683 = t684 ^ 2;
t596 = -pkin(5) * t683 + pkin(11) * t673 + t634 * t653 + t598;
t606 = -t646 * pkin(4) - t669 * qJ(5) + t671 * t661 + qJDD(5) + t615;
t627 = t646 * t697 - t647 * t693;
t628 = t646 * t693 + t647 * t697;
t599 = (-t653 * t684 - t628) * pkin(11) + (t654 * t684 - t627) * pkin(5) + t606;
t701 = sin(qJ(6));
t705 = cos(qJ(6));
t593 = -t596 * t701 + t599 * t705;
t636 = -t654 * t701 + t684 * t705;
t609 = qJD(6) * t636 + t628 * t705 + t673 * t701;
t637 = t654 * t705 + t684 * t701;
t620 = -mrSges(7,1) * t636 + mrSges(7,2) * t637;
t650 = qJD(6) - t653;
t623 = -mrSges(7,2) * t650 + mrSges(7,3) * t636;
t626 = qJDD(6) - t627;
t591 = m(7) * t593 + mrSges(7,1) * t626 - mrSges(7,3) * t609 - t620 * t637 + t623 * t650;
t594 = t596 * t705 + t599 * t701;
t608 = -qJD(6) * t637 - t628 * t701 + t673 * t705;
t624 = mrSges(7,1) * t650 - mrSges(7,3) * t637;
t592 = m(7) * t594 - mrSges(7,2) * t626 + mrSges(7,3) * t608 + t620 * t636 - t624 * t650;
t583 = t705 * t591 + t701 * t592;
t638 = -mrSges(6,2) * t684 + mrSges(6,3) * t653;
t639 = mrSges(6,1) * t684 - mrSges(6,3) * t654;
t711 = m(6) * t606 - t627 * mrSges(6,1) + mrSges(6,2) * t628 - t653 * t638 + t639 * t654 + t583;
t710 = -m(5) * t615 + t646 * mrSges(5,1) - mrSges(5,2) * t647 + t670 * t660 - t662 * t671 - t711;
t579 = m(4) * t621 + mrSges(4,1) * t690 - mrSges(4,3) * t678 + t675 * t691 - t676 * t723 + t710;
t733 = t579 * t707;
t674 = mrSges(4,1) * t691 - mrSges(4,3) * t723;
t633 = -mrSges(6,1) * t653 + mrSges(6,2) * t654;
t718 = -t591 * t701 + t705 * t592;
t582 = m(6) * t598 - mrSges(6,2) * t673 + mrSges(6,3) * t627 + t633 * t653 - t639 * t684 + t718;
t717 = -t697 * t601 + t693 * t603;
t597 = -0.2e1 * qJD(5) * t654 - t717;
t595 = -t673 * pkin(5) - t683 * pkin(11) + (t735 + t634) * t654 + t717;
t712 = -m(7) * t595 + t608 * mrSges(7,1) - mrSges(7,2) * t609 + t636 * t623 - t624 * t637;
t587 = m(6) * t597 + mrSges(6,1) * t673 - mrSges(6,3) * t628 - t633 * t654 + t638 * t684 + t712;
t576 = t693 * t582 + t697 * t587;
t657 = -mrSges(5,1) * t670 + mrSges(5,2) * t671;
t574 = m(5) * t604 + mrSges(5,1) * t673 - mrSges(5,3) * t647 - t657 * t671 + t660 * t684 + t576;
t719 = t697 * t582 - t587 * t693;
t575 = m(5) * t605 - mrSges(5,2) * t673 + mrSges(5,3) * t646 + t657 * t670 - t662 * t684 + t719;
t720 = -t574 * t702 + t706 * t575;
t565 = m(4) * t622 - mrSges(4,2) * t690 - mrSges(4,3) * t679 - t674 * t691 + t676 * t722 + t720;
t568 = t706 * t574 + t702 * t575;
t635 = -t695 * t655 + t668;
t567 = m(4) * t635 + t679 * mrSges(4,1) + t678 * mrSges(4,2) + (t674 * t703 - t675 * t707) * t725 + t568;
t554 = t565 * t729 - t567 * t695 + t699 * t733;
t550 = m(3) * t658 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t709 + t554;
t553 = t565 * t732 + t699 * t567 + t695 * t733;
t552 = m(3) * t672 + t553;
t561 = t707 * t565 - t579 * t703;
t560 = m(3) * t659 - mrSges(3,1) * t709 - qJDD(2) * mrSges(3,2) + t561;
t540 = t550 * t727 - t552 * t696 + t560 * t728;
t538 = m(2) * t685 + t540;
t546 = -t550 * t704 + t708 * t560;
t545 = m(2) * t686 + t546;
t726 = t698 * t538 + t694 * t545;
t539 = t550 * t730 + t700 * t552 + t560 * t731;
t721 = -t538 * t694 + t698 * t545;
t610 = Ifges(7,5) * t637 + Ifges(7,6) * t636 + Ifges(7,3) * t650;
t612 = Ifges(7,1) * t637 + Ifges(7,4) * t636 + Ifges(7,5) * t650;
t584 = -mrSges(7,1) * t595 + mrSges(7,3) * t594 + Ifges(7,4) * t609 + Ifges(7,2) * t608 + Ifges(7,6) * t626 - t610 * t637 + t612 * t650;
t611 = Ifges(7,4) * t637 + Ifges(7,2) * t636 + Ifges(7,6) * t650;
t585 = mrSges(7,2) * t595 - mrSges(7,3) * t593 + Ifges(7,1) * t609 + Ifges(7,4) * t608 + Ifges(7,5) * t626 + t610 * t636 - t611 * t650;
t629 = Ifges(6,5) * t654 + Ifges(6,6) * t653 + Ifges(6,3) * t684;
t630 = Ifges(6,4) * t654 + Ifges(6,2) * t653 + Ifges(6,6) * t684;
t569 = mrSges(6,2) * t606 - mrSges(6,3) * t597 + Ifges(6,1) * t628 + Ifges(6,4) * t627 + Ifges(6,5) * t673 - pkin(11) * t583 - t584 * t701 + t585 * t705 + t629 * t653 - t630 * t684;
t631 = Ifges(6,1) * t654 + Ifges(6,4) * t653 + Ifges(6,5) * t684;
t570 = -mrSges(6,1) * t606 - mrSges(7,1) * t593 + mrSges(7,2) * t594 + mrSges(6,3) * t598 + Ifges(6,4) * t628 - Ifges(7,5) * t609 + Ifges(6,2) * t627 + Ifges(6,6) * t673 - Ifges(7,6) * t608 - Ifges(7,3) * t626 - pkin(5) * t583 - t611 * t637 + t612 * t636 - t629 * t654 + t631 * t684;
t640 = Ifges(5,5) * t671 + Ifges(5,6) * t670 + Ifges(5,3) * t684;
t642 = Ifges(5,1) * t671 + Ifges(5,4) * t670 + Ifges(5,5) * t684;
t555 = -mrSges(5,1) * t615 + mrSges(5,3) * t605 + Ifges(5,4) * t647 + Ifges(5,2) * t646 + Ifges(5,6) * t673 - pkin(4) * t711 + qJ(5) * t719 + t693 * t569 + t697 * t570 - t671 * t640 + t684 * t642;
t641 = Ifges(5,4) * t671 + Ifges(5,2) * t670 + Ifges(5,6) * t684;
t556 = mrSges(5,2) * t615 - mrSges(5,3) * t604 + Ifges(5,1) * t647 + Ifges(5,4) * t646 + Ifges(5,5) * t673 - qJ(5) * t576 + t569 * t697 - t570 * t693 + t640 * t670 - t641 * t684;
t665 = Ifges(4,6) * t691 + (Ifges(4,4) * t703 + Ifges(4,2) * t707) * t725;
t666 = Ifges(4,5) * t691 + (Ifges(4,1) * t703 + Ifges(4,4) * t707) * t725;
t541 = Ifges(4,5) * t678 - Ifges(4,6) * t679 + Ifges(4,3) * t690 + mrSges(4,1) * t621 - mrSges(4,2) * t622 + t702 * t556 + t706 * t555 + pkin(3) * t710 + pkin(10) * t720 + (t665 * t703 - t666 * t707) * t725;
t664 = Ifges(4,3) * t691 + (Ifges(4,5) * t703 + Ifges(4,6) * t707) * t725;
t542 = mrSges(4,2) * t635 - mrSges(4,3) * t621 + Ifges(4,1) * t678 - Ifges(4,4) * t679 + Ifges(4,5) * t690 - pkin(10) * t568 - t555 * t702 + t556 * t706 + t664 * t722 - t665 * t691;
t547 = (-Ifges(5,3) - Ifges(6,3)) * t673 - t705 * t584 - t701 * t585 + Ifges(4,6) * t690 + t691 * t666 + Ifges(4,4) * t678 - Ifges(4,2) * t679 + t670 * t642 - t671 * t641 - Ifges(5,5) * t647 + t653 * t631 - t654 * t630 - Ifges(5,6) * t646 - Ifges(6,6) * t627 - Ifges(6,5) * t628 - mrSges(4,1) * t635 + mrSges(4,3) * t622 - mrSges(5,1) * t604 + mrSges(5,2) * t605 - mrSges(6,1) * t597 + mrSges(6,2) * t598 - pkin(4) * t576 - pkin(3) * t568 - t664 * t723 - pkin(11) * t718 - pkin(5) * t712;
t713 = pkin(9) * t561 + t542 * t703 + t547 * t707;
t535 = -mrSges(3,1) * t672 + mrSges(3,3) * t659 + t709 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t553 - t695 * t541 + t699 * t713;
t536 = mrSges(3,2) * t672 - mrSges(3,3) * t658 + Ifges(3,5) * qJDD(2) - t709 * Ifges(3,6) + t707 * t542 - t703 * t547 + (-t553 * t695 - t554 * t699) * pkin(9);
t714 = pkin(8) * t546 + t535 * t708 + t536 * t704;
t534 = mrSges(3,1) * t658 - mrSges(3,2) * t659 + Ifges(3,3) * qJDD(2) + pkin(2) * t554 + t699 * t541 + t695 * t713;
t533 = mrSges(2,2) * t692 - mrSges(2,3) * t685 - t704 * t535 + t708 * t536 + (-t539 * t696 - t540 * t700) * pkin(8);
t532 = -mrSges(2,1) * t692 + mrSges(2,3) * t686 - pkin(1) * t539 - t696 * t534 + t700 * t714;
t1 = [-m(1) * g(1) + t721; -m(1) * g(2) + t726; -m(1) * g(3) + m(2) * t692 + t539; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t726 - t694 * t532 + t698 * t533; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t721 + t698 * t532 + t694 * t533; -mrSges(1,1) * g(2) + mrSges(2,1) * t685 + mrSges(1,2) * g(1) - mrSges(2,2) * t686 + pkin(1) * t540 + t700 * t534 + t696 * t714;];
tauB  = t1;
