% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 23:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:57:23
% EndTime: 2019-05-04 22:57:33
% DurationCPUTime: 7.35s
% Computational Cost: add. (101922->301), mult. (227735->371), div. (0->0), fcn. (164016->12), ass. (0->138)
t729 = Ifges(5,1) + Ifges(6,2);
t723 = Ifges(5,4) + Ifges(6,6);
t722 = Ifges(5,5) - Ifges(6,4);
t728 = Ifges(5,2) + Ifges(6,3);
t721 = Ifges(5,6) - Ifges(6,5);
t727 = -Ifges(5,3) - Ifges(6,1);
t684 = qJD(2) ^ 2;
t673 = sin(pkin(10));
t676 = cos(pkin(10));
t658 = g(1) * t673 - g(2) * t676;
t659 = -g(1) * t676 - g(2) * t673;
t671 = -g(3) + qJDD(1);
t674 = sin(pkin(6));
t677 = cos(pkin(6));
t680 = sin(qJ(2));
t682 = cos(qJ(2));
t621 = -t680 * t659 + (t658 * t677 + t671 * t674) * t682;
t726 = -2 * qJD(5);
t725 = cos(qJ(4));
t675 = cos(pkin(11));
t724 = pkin(3) * t675;
t672 = sin(pkin(11));
t720 = mrSges(4,2) * t672;
t690 = qJDD(3) - t621;
t611 = -qJDD(2) * pkin(2) - t684 * qJ(3) + t690;
t668 = t672 ^ 2;
t679 = sin(qJ(4));
t703 = t675 * t725;
t705 = qJDD(2) * t672;
t694 = t725 * t672 + t675 * t679;
t652 = t694 * qJD(2);
t707 = t652 * qJD(4);
t634 = -qJDD(2) * t703 + t679 * t705 + t707;
t669 = t675 ^ 2;
t606 = (-pkin(2) - t724) * qJDD(2) + (-qJ(3) + (-t668 - t669) * pkin(8)) * t684 + t690;
t709 = qJD(2) * t672;
t651 = -qJD(2) * t703 + t679 * t709;
t708 = t651 * qJD(4);
t635 = t694 * qJDD(2) - t708;
t685 = pkin(4) * t707 + t652 * t726 + (-t635 + t708) * qJ(5) + t606;
t591 = t634 * pkin(4) + t685;
t643 = mrSges(6,1) * t651 - qJD(4) * mrSges(6,3);
t644 = mrSges(6,1) * t652 + qJD(4) * mrSges(6,2);
t716 = t677 * t680;
t717 = t674 * t680;
t622 = t658 * t716 + t682 * t659 + t671 * t717;
t614 = -pkin(2) * t684 + qJDD(2) * qJ(3) + t622;
t646 = -t658 * t674 + t671 * t677;
t706 = qJD(2) * qJD(3);
t710 = t675 * t646 - 0.2e1 * t672 * t706;
t595 = (-pkin(8) * qJDD(2) + t684 * t724 - t614) * t672 + t710;
t598 = t672 * t646 + (t614 + 0.2e1 * t706) * t675;
t704 = qJDD(2) * t675;
t718 = t669 * t684;
t596 = -pkin(3) * t718 + pkin(8) * t704 + t598;
t588 = t725 * t595 - t679 * t596;
t626 = pkin(4) * t651 - qJ(5) * t652;
t683 = qJD(4) ^ 2;
t587 = -qJDD(4) * pkin(4) - t683 * qJ(5) + t652 * t626 + qJDD(5) - t588;
t582 = (t651 * t652 - qJDD(4)) * pkin(9) + (t635 + t708) * pkin(5) + t587;
t645 = pkin(5) * t652 - qJD(4) * pkin(9);
t650 = t651 ^ 2;
t585 = -t650 * pkin(5) - t652 * t645 + (pkin(4) + pkin(9)) * t634 + t685;
t678 = sin(qJ(6));
t681 = cos(qJ(6));
t580 = t582 * t681 - t585 * t678;
t636 = -qJD(4) * t678 + t651 * t681;
t605 = qJD(6) * t636 + qJDD(4) * t681 + t634 * t678;
t637 = qJD(4) * t681 + t651 * t678;
t607 = -mrSges(7,1) * t636 + mrSges(7,2) * t637;
t649 = qJD(6) + t652;
t612 = -mrSges(7,2) * t649 + mrSges(7,3) * t636;
t633 = qJDD(6) + t635;
t578 = m(7) * t580 + mrSges(7,1) * t633 - mrSges(7,3) * t605 - t607 * t637 + t612 * t649;
t581 = t582 * t678 + t585 * t681;
t604 = -qJD(6) * t637 - qJDD(4) * t678 + t634 * t681;
t613 = mrSges(7,1) * t649 - mrSges(7,3) * t637;
t579 = m(7) * t581 - mrSges(7,2) * t633 + mrSges(7,3) * t604 + t607 * t636 - t613 * t649;
t714 = -t678 * t578 + t681 * t579;
t569 = m(6) * t591 - t634 * mrSges(6,2) - t635 * mrSges(6,3) - t651 * t643 - t652 * t644 + t714;
t641 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t651;
t642 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t652;
t687 = m(5) * t606 + t634 * mrSges(5,1) + t635 * mrSges(5,2) + t651 * t641 + t652 * t642 + t569;
t686 = -m(4) * t611 + mrSges(4,1) * t704 - t687 + (t668 * t684 + t718) * mrSges(4,3);
t568 = t686 + (mrSges(3,1) - t720) * qJDD(2) - t684 * mrSges(3,2) + m(3) * t621;
t719 = t568 * t682;
t627 = mrSges(5,1) * t651 + mrSges(5,2) * t652;
t570 = t681 * t578 + t678 * t579;
t628 = -mrSges(6,2) * t651 - mrSges(6,3) * t652;
t691 = -m(6) * t587 - t635 * mrSges(6,1) - t652 * t628 - t570;
t565 = m(5) * t588 - t635 * mrSges(5,3) - t652 * t627 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t641 - t643) * qJD(4) + t691;
t589 = t679 * t595 + t725 * t596;
t689 = -t683 * pkin(4) + qJDD(4) * qJ(5) - t651 * t626 + t589;
t586 = qJD(4) * t726 - t689;
t584 = -t634 * pkin(5) - t650 * pkin(9) + ((2 * qJD(5)) + t645) * qJD(4) + t689;
t692 = -m(7) * t584 + t604 * mrSges(7,1) - t605 * mrSges(7,2) + t636 * t612 - t637 * t613;
t688 = -m(6) * t586 + qJDD(4) * mrSges(6,3) + qJD(4) * t644 - t692;
t575 = m(5) * t589 - qJDD(4) * mrSges(5,2) - qJD(4) * t642 + (-t627 - t628) * t651 + (-mrSges(5,3) - mrSges(6,1)) * t634 + t688;
t563 = t725 * t565 + t679 * t575;
t597 = -t614 * t672 + t710;
t695 = mrSges(4,3) * qJDD(2) + t684 * (-mrSges(4,1) * t675 + t720);
t561 = m(4) * t597 - t695 * t672 + t563;
t700 = -t679 * t565 + t725 * t575;
t562 = m(4) * t598 + t695 * t675 + t700;
t701 = -t561 * t672 + t675 * t562;
t553 = m(3) * t622 - mrSges(3,1) * t684 - qJDD(2) * mrSges(3,2) + t701;
t556 = t675 * t561 + t672 * t562;
t555 = m(3) * t646 + t556;
t544 = t553 * t716 - t555 * t674 + t677 * t719;
t542 = m(2) * t658 + t544;
t548 = t682 * t553 - t568 * t680;
t547 = m(2) * t659 + t548;
t715 = t676 * t542 + t673 * t547;
t713 = qJD(4) * t727 + t651 * t721 - t652 * t722;
t712 = -qJD(4) * t721 + t651 * t728 - t652 * t723;
t711 = t722 * qJD(4) - t723 * t651 + t652 * t729;
t543 = t553 * t717 + t677 * t555 + t674 * t719;
t702 = -t542 * t673 + t676 * t547;
t699 = Ifges(4,1) * t672 + Ifges(4,4) * t675;
t698 = Ifges(4,4) * t672 + Ifges(4,2) * t675;
t697 = Ifges(4,5) * t672 + Ifges(4,6) * t675;
t599 = Ifges(7,5) * t637 + Ifges(7,6) * t636 + Ifges(7,3) * t649;
t601 = Ifges(7,1) * t637 + Ifges(7,4) * t636 + Ifges(7,5) * t649;
t571 = -mrSges(7,1) * t584 + mrSges(7,3) * t581 + Ifges(7,4) * t605 + Ifges(7,2) * t604 + Ifges(7,6) * t633 - t599 * t637 + t601 * t649;
t600 = Ifges(7,4) * t637 + Ifges(7,2) * t636 + Ifges(7,6) * t649;
t572 = mrSges(7,2) * t584 - mrSges(7,3) * t580 + Ifges(7,1) * t605 + Ifges(7,4) * t604 + Ifges(7,5) * t633 + t599 * t636 - t600 * t649;
t549 = -mrSges(5,1) * t606 - mrSges(6,1) * t586 + mrSges(6,2) * t591 + mrSges(5,3) * t589 - pkin(4) * t569 - pkin(5) * t692 - pkin(9) * t714 + t711 * qJD(4) + t721 * qJDD(4) - t681 * t571 - t678 * t572 - t634 * t728 + t723 * t635 + t713 * t652;
t557 = mrSges(6,1) * t587 + mrSges(7,1) * t580 + mrSges(5,2) * t606 - mrSges(7,2) * t581 - mrSges(5,3) * t588 - mrSges(6,3) * t591 + Ifges(7,5) * t605 + Ifges(7,6) * t604 + Ifges(7,3) * t633 + pkin(5) * t570 - qJ(5) * t569 + t637 * t600 - t636 * t601 + t713 * t651 + t729 * t635 - t723 * t634 + t722 * qJDD(4) + t712 * qJD(4);
t657 = t697 * qJD(2);
t538 = -mrSges(4,1) * t611 + mrSges(4,3) * t598 - pkin(3) * t687 + pkin(8) * t700 + t698 * qJDD(2) + t725 * t549 + t679 * t557 - t657 * t709;
t540 = t675 * qJD(2) * t657 + mrSges(4,2) * t611 - mrSges(4,3) * t597 - pkin(8) * t563 + t699 * qJDD(2) - t679 * t549 + t725 * t557;
t537 = mrSges(3,2) * t646 - mrSges(3,3) * t621 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t684 - qJ(3) * t556 - t538 * t672 + t540 * t675;
t539 = t678 * t571 - pkin(4) * (-qJD(4) * t643 + t691) - t681 * t572 - qJ(5) * t688 - mrSges(3,1) * t646 + mrSges(3,3) * t622 - mrSges(4,1) * t597 + mrSges(4,2) * t598 + mrSges(6,3) * t586 - mrSges(6,2) * t587 - mrSges(5,1) * t588 + mrSges(5,2) * t589 + pkin(9) * t570 - pkin(3) * t563 - pkin(2) * t556 + t712 * t652 + (qJ(5) * t628 - t711) * t651 - t722 * t635 + (mrSges(6,1) * qJ(5) + t721) * t634 + (mrSges(6,2) * pkin(4) + t727) * qJDD(4) + (Ifges(3,6) - t697) * qJDD(2) + (-t672 * t698 + t675 * t699 + Ifges(3,5)) * t684;
t693 = pkin(7) * t548 + t537 * t680 + t539 * t682;
t536 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t621 - mrSges(3,2) * t622 + t672 * t540 + t675 * t538 + pkin(2) * (-mrSges(4,2) * t705 + t686) + qJ(3) * t701;
t535 = mrSges(2,2) * t671 - mrSges(2,3) * t658 + t682 * t537 - t680 * t539 + (-t543 * t674 - t544 * t677) * pkin(7);
t534 = -mrSges(2,1) * t671 + mrSges(2,3) * t659 - pkin(1) * t543 - t674 * t536 + t693 * t677;
t1 = [-m(1) * g(1) + t702; -m(1) * g(2) + t715; -m(1) * g(3) + m(2) * t671 + t543; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t715 - t673 * t534 + t676 * t535; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t702 + t676 * t534 + t673 * t535; -mrSges(1,1) * g(2) + mrSges(2,1) * t658 + mrSges(1,2) * g(1) - mrSges(2,2) * t659 + pkin(1) * t544 + t677 * t536 + t693 * t674;];
tauB  = t1;
