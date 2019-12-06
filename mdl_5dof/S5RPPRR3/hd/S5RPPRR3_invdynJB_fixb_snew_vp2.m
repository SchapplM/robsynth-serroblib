% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:59
% EndTime: 2019-12-05 17:42:04
% DurationCPUTime: 4.70s
% Computational Cost: add. (55727->249), mult. (123458->311), div. (0->0), fcn. (83128->10), ass. (0->112)
t684 = qJD(1) ^ 2;
t676 = cos(pkin(9));
t711 = pkin(3) * t676;
t674 = sin(pkin(9));
t710 = mrSges(4,2) * t674;
t669 = t676 ^ 2;
t709 = t669 * t684;
t680 = sin(qJ(1));
t683 = cos(qJ(1));
t654 = t683 * g(2) + t680 * g(3);
t650 = qJDD(1) * pkin(1) + t654;
t653 = t680 * g(2) - t683 * g(3);
t651 = -t684 * pkin(1) + t653;
t675 = sin(pkin(8));
t677 = cos(pkin(8));
t638 = t675 * t650 + t677 * t651;
t630 = -t684 * pkin(2) + qJDD(1) * qJ(3) + t638;
t673 = -g(1) + qJDD(2);
t705 = qJD(1) * qJD(3);
t708 = t676 * t673 - 0.2e1 * t674 * t705;
t616 = (-pkin(6) * qJDD(1) + t684 * t711 - t630) * t674 + t708;
t620 = t674 * t673 + (t630 + 0.2e1 * t705) * t676;
t704 = qJDD(1) * t676;
t617 = -pkin(3) * t709 + pkin(6) * t704 + t620;
t679 = sin(qJ(4));
t682 = cos(qJ(4));
t598 = t682 * t616 - t679 * t617;
t692 = t674 * t682 + t676 * t679;
t691 = -t674 * t679 + t676 * t682;
t643 = t691 * qJD(1);
t706 = t643 * qJD(4);
t635 = t692 * qJDD(1) + t706;
t644 = t692 * qJD(1);
t594 = (-t635 + t706) * pkin(7) + (t643 * t644 + qJDD(4)) * pkin(4) + t598;
t599 = t679 * t616 + t682 * t617;
t634 = -t644 * qJD(4) + t691 * qJDD(1);
t641 = qJD(4) * pkin(4) - t644 * pkin(7);
t642 = t643 ^ 2;
t595 = -t642 * pkin(4) + t634 * pkin(7) - qJD(4) * t641 + t599;
t678 = sin(qJ(5));
t681 = cos(qJ(5));
t592 = t681 * t594 - t678 * t595;
t628 = t681 * t643 - t678 * t644;
t606 = t628 * qJD(5) + t678 * t634 + t681 * t635;
t629 = t678 * t643 + t681 * t644;
t612 = -t628 * mrSges(6,1) + t629 * mrSges(6,2);
t670 = qJD(4) + qJD(5);
t621 = -t670 * mrSges(6,2) + t628 * mrSges(6,3);
t667 = qJDD(4) + qJDD(5);
t589 = m(6) * t592 + t667 * mrSges(6,1) - t606 * mrSges(6,3) - t629 * t612 + t670 * t621;
t593 = t678 * t594 + t681 * t595;
t605 = -t629 * qJD(5) + t681 * t634 - t678 * t635;
t622 = t670 * mrSges(6,1) - t629 * mrSges(6,3);
t590 = m(6) * t593 - t667 * mrSges(6,2) + t605 * mrSges(6,3) + t628 * t612 - t670 * t622;
t579 = t681 * t589 + t678 * t590;
t632 = -t643 * mrSges(5,1) + t644 * mrSges(5,2);
t639 = -qJD(4) * mrSges(5,2) + t643 * mrSges(5,3);
t577 = m(5) * t598 + qJDD(4) * mrSges(5,1) - t635 * mrSges(5,3) + qJD(4) * t639 - t644 * t632 + t579;
t640 = qJD(4) * mrSges(5,1) - t644 * mrSges(5,3);
t699 = -t678 * t589 + t681 * t590;
t578 = m(5) * t599 - qJDD(4) * mrSges(5,2) + t634 * mrSges(5,3) - qJD(4) * t640 + t643 * t632 + t699;
t573 = t682 * t577 + t679 * t578;
t619 = -t674 * t630 + t708;
t690 = mrSges(4,3) * qJDD(1) + t684 * (-mrSges(4,1) * t676 + t710);
t571 = m(4) * t619 - t690 * t674 + t573;
t700 = -t679 * t577 + t682 * t578;
t572 = m(4) * t620 + t690 * t676 + t700;
t701 = -t674 * t571 + t676 * t572;
t562 = m(3) * t638 - t684 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t701;
t637 = t677 * t650 - t675 * t651;
t695 = qJDD(3) - t637;
t624 = -qJDD(1) * pkin(2) - t684 * qJ(3) + t695;
t668 = t674 ^ 2;
t618 = (-pkin(2) - t711) * qJDD(1) + (-qJ(3) + (-t668 - t669) * pkin(6)) * t684 + t695;
t597 = -t634 * pkin(4) - t642 * pkin(7) + t644 * t641 + t618;
t694 = m(6) * t597 - t605 * mrSges(6,1) + t606 * mrSges(6,2) - t628 * t621 + t629 * t622;
t687 = m(5) * t618 - t634 * mrSges(5,1) + t635 * mrSges(5,2) - t643 * t639 + t644 * t640 + t694;
t686 = -m(4) * t624 + mrSges(4,1) * t704 - t687 + (t668 * t684 + t709) * mrSges(4,3);
t583 = t686 + (mrSges(3,1) - t710) * qJDD(1) - t684 * mrSges(3,2) + m(3) * t637;
t559 = t675 * t562 + t677 * t583;
t565 = t676 * t571 + t674 * t572;
t696 = Ifges(4,5) * t674 + Ifges(4,6) * t676;
t707 = t684 * t696;
t563 = m(3) * t673 + t565;
t702 = t677 * t562 - t675 * t583;
t556 = m(2) * t653 - t684 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t702;
t557 = m(2) * t654 + qJDD(1) * mrSges(2,1) - t684 * mrSges(2,2) + t559;
t703 = t683 * t556 - t680 * t557;
t698 = Ifges(4,1) * t674 + Ifges(4,4) * t676;
t697 = Ifges(4,4) * t674 + Ifges(4,2) * t676;
t693 = -t680 * t556 - t683 * t557;
t608 = Ifges(6,4) * t629 + Ifges(6,2) * t628 + Ifges(6,6) * t670;
t609 = Ifges(6,1) * t629 + Ifges(6,4) * t628 + Ifges(6,5) * t670;
t689 = -mrSges(6,1) * t592 + mrSges(6,2) * t593 - Ifges(6,5) * t606 - Ifges(6,6) * t605 - Ifges(6,3) * t667 - t629 * t608 + t628 * t609;
t607 = Ifges(6,5) * t629 + Ifges(6,6) * t628 + Ifges(6,3) * t670;
t580 = -mrSges(6,1) * t597 + mrSges(6,3) * t593 + Ifges(6,4) * t606 + Ifges(6,2) * t605 + Ifges(6,6) * t667 - t629 * t607 + t670 * t609;
t581 = mrSges(6,2) * t597 - mrSges(6,3) * t592 + Ifges(6,1) * t606 + Ifges(6,4) * t605 + Ifges(6,5) * t667 + t628 * t607 - t670 * t608;
t625 = Ifges(5,5) * t644 + Ifges(5,6) * t643 + Ifges(5,3) * qJD(4);
t627 = Ifges(5,1) * t644 + Ifges(5,4) * t643 + Ifges(5,5) * qJD(4);
t566 = -mrSges(5,1) * t618 + mrSges(5,3) * t599 + Ifges(5,4) * t635 + Ifges(5,2) * t634 + Ifges(5,6) * qJDD(4) - pkin(4) * t694 + pkin(7) * t699 + qJD(4) * t627 + t681 * t580 + t678 * t581 - t644 * t625;
t626 = Ifges(5,4) * t644 + Ifges(5,2) * t643 + Ifges(5,6) * qJD(4);
t567 = mrSges(5,2) * t618 - mrSges(5,3) * t598 + Ifges(5,1) * t635 + Ifges(5,4) * t634 + Ifges(5,5) * qJDD(4) - pkin(7) * t579 - qJD(4) * t626 - t678 * t580 + t681 * t581 + t643 * t625;
t552 = -mrSges(4,1) * t624 + mrSges(4,3) * t620 - pkin(3) * t687 + pkin(6) * t700 + t697 * qJDD(1) + t682 * t566 + t679 * t567 - t674 * t707;
t554 = mrSges(4,2) * t624 - mrSges(4,3) * t619 - pkin(6) * t573 + t698 * qJDD(1) - t679 * t566 + t682 * t567 + t676 * t707;
t585 = qJDD(1) * t710 - t686;
t688 = mrSges(2,1) * t654 + mrSges(3,1) * t637 - mrSges(2,2) * t653 - mrSges(3,2) * t638 + pkin(1) * t559 - pkin(2) * t585 + qJ(3) * t701 + t676 * t552 + t674 * t554 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t685 = mrSges(5,1) * t598 - mrSges(5,2) * t599 + Ifges(5,5) * t635 + Ifges(5,6) * t634 + Ifges(5,3) * qJDD(4) + pkin(4) * t579 + t644 * t626 - t643 * t627 - t689;
t550 = -t685 + (Ifges(3,6) - t696) * qJDD(1) - mrSges(3,1) * t673 + mrSges(3,3) * t638 - mrSges(4,1) * t619 + mrSges(4,2) * t620 - pkin(2) * t565 - pkin(3) * t573 + (-t674 * t697 + t676 * t698 + Ifges(3,5)) * t684;
t549 = mrSges(3,2) * t673 - mrSges(3,3) * t637 + Ifges(3,5) * qJDD(1) - t684 * Ifges(3,6) - qJ(3) * t565 - t674 * t552 + t676 * t554;
t548 = -mrSges(2,2) * g(1) - mrSges(2,3) * t654 + Ifges(2,5) * qJDD(1) - t684 * Ifges(2,6) - qJ(2) * t559 + t677 * t549 - t675 * t550;
t547 = mrSges(2,1) * g(1) + mrSges(2,3) * t653 + t684 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t563 + qJ(2) * t702 + t675 * t549 + t677 * t550;
t1 = [(-m(1) - m(2)) * g(1) + t563; -m(1) * g(2) + t693; -m(1) * g(3) + t703; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t688; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t703 - t683 * t547 - t680 * t548; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t693 - t680 * t547 + t683 * t548; t688; t563; t585; t685; -t689;];
tauJB = t1;
