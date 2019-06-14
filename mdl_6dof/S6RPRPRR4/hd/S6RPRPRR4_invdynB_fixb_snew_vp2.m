% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:43:05
% EndTime: 2019-05-05 18:43:13
% DurationCPUTime: 6.07s
% Computational Cost: add. (72599->328), mult. (144038->392), div. (0->0), fcn. (83090->10), ass. (0->135)
t721 = -2 * qJD(4);
t720 = Ifges(4,1) + Ifges(5,2);
t714 = Ifges(4,4) + Ifges(5,6);
t713 = Ifges(4,5) - Ifges(5,4);
t719 = Ifges(4,2) + Ifges(5,3);
t712 = Ifges(4,6) - Ifges(5,5);
t718 = (Ifges(4,3) + Ifges(5,1));
t679 = sin(qJ(1));
t683 = cos(qJ(1));
t658 = t679 * g(1) - t683 * g(2);
t644 = qJDD(1) * pkin(1) + t658;
t659 = -t683 * g(1) - t679 * g(2);
t685 = qJD(1) ^ 2;
t648 = -t685 * pkin(1) + t659;
t674 = sin(pkin(10));
t675 = cos(pkin(10));
t617 = t674 * t644 + t675 * t648;
t606 = -t685 * pkin(2) + qJDD(1) * pkin(7) + t617;
t673 = -g(3) + qJDD(2);
t678 = sin(qJ(3));
t682 = cos(qJ(3));
t600 = t682 * t606 + t678 * t673;
t645 = (-pkin(3) * t682 - qJ(4) * t678) * qJD(1);
t684 = qJD(3) ^ 2;
t705 = qJD(1) * t682;
t596 = t684 * pkin(3) - qJDD(3) * qJ(4) + (qJD(3) * t721) - t645 * t705 - t600;
t717 = -pkin(3) - pkin(8);
t716 = pkin(8) * t685;
t715 = t685 * pkin(7);
t711 = t682 * t673;
t603 = t678 * t606;
t599 = -t603 + t711;
t646 = (mrSges(5,2) * t682 - mrSges(5,3) * t678) * qJD(1);
t647 = (-mrSges(4,1) * t682 + mrSges(4,2) * t678) * qJD(1);
t704 = qJD(1) * qJD(3);
t700 = t682 * t704;
t649 = t678 * qJDD(1) + t700;
t654 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t705;
t655 = -mrSges(5,1) * t705 - qJD(3) * mrSges(5,3);
t666 = t678 * qJD(1);
t701 = t678 * t704;
t650 = t682 * qJDD(1) - t701;
t657 = pkin(4) * t666 - qJD(3) * pkin(8);
t672 = t682 ^ 2;
t616 = t675 * t644 - t674 * t648;
t694 = -qJDD(1) * pkin(2) - t616;
t689 = pkin(3) * t701 + t666 * t721 + (-t649 - t700) * qJ(4) + t694;
t582 = -t657 * t666 + (-pkin(4) * t672 - pkin(7)) * t685 + t717 * t650 + t689;
t695 = -t684 * qJ(4) + t645 * t666 + qJDD(4) + t603;
t590 = t649 * pkin(4) + t717 * qJDD(3) + (-pkin(4) * t704 - t678 * t716 - t673) * t682 + t695;
t677 = sin(qJ(5));
t681 = cos(qJ(5));
t577 = -t677 * t582 + t681 * t590;
t642 = -t677 * qJD(3) - t681 * t705;
t613 = t642 * qJD(5) + t681 * qJDD(3) - t677 * t650;
t641 = qJDD(5) + t649;
t643 = t681 * qJD(3) - t677 * t705;
t662 = t666 + qJD(5);
t575 = (t642 * t662 - t613) * pkin(9) + (t642 * t643 + t641) * pkin(5) + t577;
t578 = t681 * t582 + t677 * t590;
t612 = -t643 * qJD(5) - t677 * qJDD(3) - t681 * t650;
t621 = t662 * pkin(5) - t643 * pkin(9);
t640 = t642 ^ 2;
t576 = -t640 * pkin(5) + t612 * pkin(9) - t662 * t621 + t578;
t676 = sin(qJ(6));
t680 = cos(qJ(6));
t573 = t680 * t575 - t676 * t576;
t614 = t680 * t642 - t676 * t643;
t585 = t614 * qJD(6) + t676 * t612 + t680 * t613;
t615 = t676 * t642 + t680 * t643;
t598 = -t614 * mrSges(7,1) + t615 * mrSges(7,2);
t660 = qJD(6) + t662;
t601 = -t660 * mrSges(7,2) + t614 * mrSges(7,3);
t634 = qJDD(6) + t641;
t571 = m(7) * t573 + t634 * mrSges(7,1) - t585 * mrSges(7,3) - t615 * t598 + t660 * t601;
t574 = t676 * t575 + t680 * t576;
t584 = -t615 * qJD(6) + t680 * t612 - t676 * t613;
t602 = t660 * mrSges(7,1) - t615 * mrSges(7,3);
t572 = m(7) * t574 - t634 * mrSges(7,2) + t584 * mrSges(7,3) + t614 * t598 - t660 * t602;
t562 = t680 * t571 + t676 * t572;
t618 = -t642 * mrSges(6,1) + t643 * mrSges(6,2);
t619 = -t662 * mrSges(6,2) + t642 * mrSges(6,3);
t560 = m(6) * t577 + t641 * mrSges(6,1) - t613 * mrSges(6,3) - t643 * t618 + t662 * t619 + t562;
t620 = t662 * mrSges(6,1) - t643 * mrSges(6,3);
t696 = -t676 * t571 + t680 * t572;
t561 = m(6) * t578 - t641 * mrSges(6,2) + t612 * mrSges(6,3) + t642 * t618 - t662 * t620 + t696;
t557 = t681 * t560 + t677 * t561;
t597 = -qJDD(3) * pkin(3) + t695 - t711;
t690 = -m(5) * t597 - t649 * mrSges(5,1) - t557;
t555 = m(4) * t599 - t649 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t654 - t655) * qJD(3) + (-t646 - t647) * t666 + t690;
t653 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t666;
t656 = mrSges(5,1) * t666 + qJD(3) * mrSges(5,2);
t589 = t650 * pkin(4) + qJD(3) * t657 - t672 * t716 - t596;
t580 = -t612 * pkin(5) - t640 * pkin(9) + t643 * t621 + t589;
t691 = m(7) * t580 - t584 * mrSges(7,1) + t585 * mrSges(7,2) - t614 * t601 + t615 * t602;
t688 = -m(6) * t589 + t612 * mrSges(6,1) - t613 * mrSges(6,2) + t642 * t619 - t643 * t620 - t691;
t687 = -m(5) * t596 + qJDD(3) * mrSges(5,3) + qJD(3) * t656 + t646 * t705 - t688;
t567 = t687 - qJDD(3) * mrSges(4,2) - qJD(3) * t653 + m(4) * t600 + (mrSges(4,3) + mrSges(5,1)) * t650 + t647 * t705;
t697 = -t678 * t555 + t682 * t567;
t548 = m(3) * t617 - t685 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t697;
t605 = t694 - t715;
t591 = -t650 * pkin(3) + t689 - t715;
t709 = -t677 * t560 + t681 * t561;
t693 = -m(5) * t591 - t650 * mrSges(5,2) + t656 * t666 - t709;
t686 = -m(4) * t605 + t654 * t705 + t650 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t649 + (-t653 * t678 - t655 * t682) * qJD(1) + t693;
t553 = m(3) * t616 + qJDD(1) * mrSges(3,1) - t685 * mrSges(3,2) + t686;
t545 = t674 * t548 + t675 * t553;
t543 = m(2) * t658 + qJDD(1) * mrSges(2,1) - t685 * mrSges(2,2) + t545;
t698 = t675 * t548 - t674 * t553;
t544 = m(2) * t659 - t685 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t698;
t710 = t683 * t543 + t679 * t544;
t549 = t682 * t555 + t678 * t567;
t708 = (t718 * qJD(3)) + (t713 * t678 + t712 * t682) * qJD(1);
t707 = -t712 * qJD(3) + (-t714 * t678 - t719 * t682) * qJD(1);
t706 = t713 * qJD(3) + (t720 * t678 + t714 * t682) * qJD(1);
t702 = m(3) * t673 + t549;
t699 = -t679 * t543 + t683 * t544;
t609 = Ifges(6,1) * t643 + Ifges(6,4) * t642 + Ifges(6,5) * t662;
t608 = Ifges(6,4) * t643 + Ifges(6,2) * t642 + Ifges(6,6) * t662;
t607 = Ifges(6,5) * t643 + Ifges(6,6) * t642 + Ifges(6,3) * t662;
t594 = Ifges(7,1) * t615 + Ifges(7,4) * t614 + Ifges(7,5) * t660;
t593 = Ifges(7,4) * t615 + Ifges(7,2) * t614 + Ifges(7,6) * t660;
t592 = Ifges(7,5) * t615 + Ifges(7,6) * t614 + Ifges(7,3) * t660;
t564 = mrSges(7,2) * t580 - mrSges(7,3) * t573 + Ifges(7,1) * t585 + Ifges(7,4) * t584 + Ifges(7,5) * t634 + t614 * t592 - t660 * t593;
t563 = -mrSges(7,1) * t580 + mrSges(7,3) * t574 + Ifges(7,4) * t585 + Ifges(7,2) * t584 + Ifges(7,6) * t634 - t615 * t592 + t660 * t594;
t556 = -t649 * mrSges(5,3) + t655 * t705 - t693;
t551 = mrSges(6,2) * t589 - mrSges(6,3) * t577 + Ifges(6,1) * t613 + Ifges(6,4) * t612 + Ifges(6,5) * t641 - pkin(9) * t562 - t676 * t563 + t680 * t564 + t642 * t607 - t662 * t608;
t550 = -mrSges(6,1) * t589 + mrSges(6,3) * t578 + Ifges(6,4) * t613 + Ifges(6,2) * t612 + Ifges(6,6) * t641 - pkin(5) * t691 + pkin(9) * t696 + t680 * t563 + t676 * t564 - t643 * t607 + t662 * t609;
t539 = t643 * t608 + Ifges(6,3) * t641 - t642 * t609 + Ifges(7,3) * t634 + Ifges(6,6) * t612 + Ifges(6,5) * t613 - t614 * t594 + t615 * t593 - mrSges(4,3) * t599 + mrSges(4,2) * t605 - mrSges(5,3) * t591 + mrSges(5,1) * t597 + Ifges(7,6) * t584 + Ifges(7,5) * t585 + mrSges(6,1) * t577 - mrSges(6,2) * t578 + mrSges(7,1) * t573 - mrSges(7,2) * t574 + pkin(5) * t562 + pkin(4) * t557 - qJ(4) * t556 + t720 * t649 + t707 * qJD(3) + t708 * t705 + t713 * qJDD(3) + t714 * t650;
t538 = -mrSges(4,1) * t605 - mrSges(5,1) * t596 + mrSges(5,2) * t591 + mrSges(4,3) * t600 - pkin(3) * t556 - pkin(4) * t688 - pkin(8) * t709 + t706 * qJD(3) + t712 * qJDD(3) - t681 * t550 - t677 * t551 + t714 * t649 + t719 * t650 - t708 * t666;
t537 = -pkin(2) * t549 - mrSges(3,1) * t673 + mrSges(3,3) * t617 - pkin(3) * (-qJD(3) * t655 + t690) - qJ(4) * t687 - t681 * t551 + t677 * t550 + pkin(8) * t557 - mrSges(4,1) * t599 + mrSges(4,2) * t600 - mrSges(5,2) * t597 + mrSges(5,3) * t596 + t685 * Ifges(3,5) + Ifges(3,6) * qJDD(1) + (-qJ(4) * mrSges(5,1) - t712) * t650 - t713 * t649 + (pkin(3) * mrSges(5,2) - t718) * qJDD(3) + (t706 * t682 + (pkin(3) * t646 + t707) * t678) * qJD(1);
t536 = mrSges(3,2) * t673 - mrSges(3,3) * t616 + Ifges(3,5) * qJDD(1) - t685 * Ifges(3,6) - pkin(7) * t549 - t678 * t538 + t682 * t539;
t535 = -mrSges(2,2) * g(3) - mrSges(2,3) * t658 + Ifges(2,5) * qJDD(1) - t685 * Ifges(2,6) - qJ(2) * t545 + t675 * t536 - t674 * t537;
t534 = mrSges(2,1) * g(3) + mrSges(2,3) * t659 + t685 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t702 + qJ(2) * t698 + t674 * t536 + t675 * t537;
t1 = [-m(1) * g(1) + t699; -m(1) * g(2) + t710; (-m(1) - m(2)) * g(3) + t702; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t710 - t679 * t534 + t683 * t535; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t699 + t683 * t534 + t679 * t535; pkin(1) * t545 + mrSges(2,1) * t658 - mrSges(2,2) * t659 + t678 * t539 + t682 * t538 + pkin(2) * t686 + pkin(7) * t697 + mrSges(3,1) * t616 - mrSges(3,2) * t617 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
