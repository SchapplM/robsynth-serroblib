% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:13
% EndTime: 2020-01-03 12:09:19
% DurationCPUTime: 6.17s
% Computational Cost: add. (107420->273), mult. (143754->347), div. (0->0), fcn. (91374->10), ass. (0->112)
t662 = qJDD(1) + qJDD(2);
t670 = sin(qJ(3));
t674 = cos(qJ(3));
t664 = qJD(1) + qJD(2);
t691 = qJD(3) * t664;
t644 = t670 * t662 + t674 * t691;
t672 = sin(qJ(1));
t676 = cos(qJ(1));
t656 = -t676 * g(2) - t672 * g(3);
t649 = qJDD(1) * pkin(1) + t656;
t655 = -t672 * g(2) + t676 * g(3);
t677 = qJD(1) ^ 2;
t650 = -t677 * pkin(1) + t655;
t671 = sin(qJ(2));
t675 = cos(qJ(2));
t628 = t671 * t649 + t675 * t650;
t660 = t664 ^ 2;
t622 = -t660 * pkin(2) + t662 * pkin(7) + t628;
t693 = t670 * t622;
t697 = pkin(3) * t660;
t606 = qJDD(3) * pkin(3) - t644 * qJ(4) - t693 + (qJ(4) * t691 + t670 * t697 - g(1)) * t674;
t612 = -t670 * g(1) + t674 * t622;
t645 = t674 * t662 - t670 * t691;
t695 = t664 * t670;
t651 = qJD(3) * pkin(3) - qJ(4) * t695;
t666 = t674 ^ 2;
t607 = t645 * qJ(4) - qJD(3) * t651 - t666 * t697 + t612;
t667 = sin(pkin(9));
t668 = cos(pkin(9));
t636 = (t667 * t674 + t668 * t670) * t664;
t586 = -0.2e1 * qJD(4) * t636 + t668 * t606 - t667 * t607;
t625 = t668 * t644 + t667 * t645;
t635 = (-t667 * t670 + t668 * t674) * t664;
t584 = (qJD(3) * t635 - t625) * pkin(8) + (t635 * t636 + qJDD(3)) * pkin(4) + t586;
t587 = 0.2e1 * qJD(4) * t635 + t667 * t606 + t668 * t607;
t624 = -t667 * t644 + t668 * t645;
t631 = qJD(3) * pkin(4) - t636 * pkin(8);
t634 = t635 ^ 2;
t585 = -t634 * pkin(4) + t624 * pkin(8) - qJD(3) * t631 + t587;
t669 = sin(qJ(5));
t673 = cos(qJ(5));
t582 = t673 * t584 - t669 * t585;
t616 = t673 * t635 - t669 * t636;
t596 = t616 * qJD(5) + t669 * t624 + t673 * t625;
t617 = t669 * t635 + t673 * t636;
t602 = -t616 * mrSges(6,1) + t617 * mrSges(6,2);
t663 = qJD(3) + qJD(5);
t609 = -t663 * mrSges(6,2) + t616 * mrSges(6,3);
t661 = qJDD(3) + qJDD(5);
t578 = m(6) * t582 + t661 * mrSges(6,1) - t596 * mrSges(6,3) - t617 * t602 + t663 * t609;
t583 = t669 * t584 + t673 * t585;
t595 = -t617 * qJD(5) + t673 * t624 - t669 * t625;
t610 = t663 * mrSges(6,1) - t617 * mrSges(6,3);
t579 = m(6) * t583 - t661 * mrSges(6,2) + t595 * mrSges(6,3) + t616 * t602 - t663 * t610;
t569 = t673 * t578 + t669 * t579;
t620 = -t635 * mrSges(5,1) + t636 * mrSges(5,2);
t629 = -qJD(3) * mrSges(5,2) + t635 * mrSges(5,3);
t567 = m(5) * t586 + qJDD(3) * mrSges(5,1) - t625 * mrSges(5,3) + qJD(3) * t629 - t636 * t620 + t569;
t630 = qJD(3) * mrSges(5,1) - t636 * mrSges(5,3);
t686 = -t669 * t578 + t673 * t579;
t568 = m(5) * t587 - qJDD(3) * mrSges(5,2) + t624 * mrSges(5,3) - qJD(3) * t630 + t635 * t620 + t686;
t563 = t668 * t567 + t667 * t568;
t611 = -t674 * g(1) - t693;
t614 = Ifges(5,4) * t636 + Ifges(5,2) * t635 + Ifges(5,6) * qJD(3);
t615 = Ifges(5,1) * t636 + Ifges(5,4) * t635 + Ifges(5,5) * qJD(3);
t638 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t670 + Ifges(4,2) * t674) * t664;
t639 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t670 + Ifges(4,4) * t674) * t664;
t598 = Ifges(6,4) * t617 + Ifges(6,2) * t616 + Ifges(6,6) * t663;
t599 = Ifges(6,1) * t617 + Ifges(6,4) * t616 + Ifges(6,5) * t663;
t681 = -mrSges(6,1) * t582 + mrSges(6,2) * t583 - Ifges(6,5) * t596 - Ifges(6,6) * t595 - Ifges(6,3) * t661 - t617 * t598 + t616 * t599;
t698 = mrSges(4,1) * t611 + mrSges(5,1) * t586 - mrSges(4,2) * t612 - mrSges(5,2) * t587 + Ifges(4,5) * t644 + Ifges(5,5) * t625 + Ifges(4,6) * t645 + Ifges(5,6) * t624 + pkin(3) * t563 + pkin(4) * t569 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t636 * t614 - t635 * t615 + (t670 * t638 - t674 * t639) * t664 - t681;
t694 = t664 * t674;
t643 = (-mrSges(4,1) * t674 + mrSges(4,2) * t670) * t664;
t653 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t694;
t561 = m(4) * t611 + qJDD(3) * mrSges(4,1) - t644 * mrSges(4,3) + qJD(3) * t653 - t643 * t695 + t563;
t652 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t695;
t687 = -t667 * t567 + t668 * t568;
t562 = m(4) * t612 - qJDD(3) * mrSges(4,2) + t645 * mrSges(4,3) - qJD(3) * t652 + t643 * t694 + t687;
t688 = -t670 * t561 + t674 * t562;
t553 = m(3) * t628 - t660 * mrSges(3,1) - t662 * mrSges(3,2) + t688;
t627 = t675 * t649 - t671 * t650;
t683 = -t662 * pkin(2) - t627;
t608 = -t645 * pkin(3) + qJDD(4) + t651 * t695 + (-qJ(4) * t666 - pkin(7)) * t660 + t683;
t589 = -t624 * pkin(4) - t634 * pkin(8) + t636 * t631 + t608;
t685 = m(6) * t589 - t595 * mrSges(6,1) + t596 * mrSges(6,2) - t616 * t609 + t617 * t610;
t580 = m(5) * t608 - t624 * mrSges(5,1) + t625 * mrSges(5,2) - t635 * t629 + t636 * t630 + t685;
t621 = -t660 * pkin(7) + t683;
t679 = -m(4) * t621 + t645 * mrSges(4,1) - t644 * mrSges(4,2) - t652 * t695 + t653 * t694 - t580;
t573 = m(3) * t627 + t662 * mrSges(3,1) - t660 * mrSges(3,2) + t679;
t689 = t675 * t553 - t671 * t573;
t547 = m(2) * t655 - t677 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t689;
t550 = t671 * t553 + t675 * t573;
t548 = m(2) * t656 + qJDD(1) * mrSges(2,1) - t677 * mrSges(2,2) + t550;
t692 = t672 * t547 + t676 * t548;
t555 = t674 * t561 + t670 * t562;
t690 = -t676 * t547 + t672 * t548;
t597 = Ifges(6,5) * t617 + Ifges(6,6) * t616 + Ifges(6,3) * t663;
t570 = -mrSges(6,1) * t589 + mrSges(6,3) * t583 + Ifges(6,4) * t596 + Ifges(6,2) * t595 + Ifges(6,6) * t661 - t617 * t597 + t663 * t599;
t571 = mrSges(6,2) * t589 - mrSges(6,3) * t582 + Ifges(6,1) * t596 + Ifges(6,4) * t595 + Ifges(6,5) * t661 + t616 * t597 - t663 * t598;
t613 = Ifges(5,5) * t636 + Ifges(5,6) * t635 + Ifges(5,3) * qJD(3);
t556 = -mrSges(5,1) * t608 + mrSges(5,3) * t587 + Ifges(5,4) * t625 + Ifges(5,2) * t624 + Ifges(5,6) * qJDD(3) - pkin(4) * t685 + pkin(8) * t686 + qJD(3) * t615 + t673 * t570 + t669 * t571 - t636 * t613;
t557 = mrSges(5,2) * t608 - mrSges(5,3) * t586 + Ifges(5,1) * t625 + Ifges(5,4) * t624 + Ifges(5,5) * qJDD(3) - pkin(8) * t569 - qJD(3) * t614 - t669 * t570 + t673 * t571 + t635 * t613;
t637 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t670 + Ifges(4,6) * t674) * t664;
t541 = -mrSges(4,1) * t621 + mrSges(4,3) * t612 + Ifges(4,4) * t644 + Ifges(4,2) * t645 + Ifges(4,6) * qJDD(3) - pkin(3) * t580 + qJ(4) * t687 + qJD(3) * t639 + t668 * t556 + t667 * t557 - t637 * t695;
t543 = mrSges(4,2) * t621 - mrSges(4,3) * t611 + Ifges(4,1) * t644 + Ifges(4,4) * t645 + Ifges(4,5) * qJDD(3) - qJ(4) * t563 - qJD(3) * t638 - t667 * t556 + t668 * t557 + t637 * t694;
t682 = mrSges(3,1) * t627 - mrSges(3,2) * t628 + Ifges(3,3) * t662 + pkin(2) * t679 + pkin(7) * t688 + t674 * t541 + t670 * t543;
t680 = mrSges(2,1) * t656 - mrSges(2,2) * t655 + Ifges(2,3) * qJDD(1) + pkin(1) * t550 + t682;
t539 = mrSges(3,1) * g(1) + mrSges(3,3) * t628 + t660 * Ifges(3,5) + Ifges(3,6) * t662 - pkin(2) * t555 - t698;
t538 = -mrSges(3,2) * g(1) - mrSges(3,3) * t627 + Ifges(3,5) * t662 - t660 * Ifges(3,6) - pkin(7) * t555 - t670 * t541 + t674 * t543;
t537 = -mrSges(2,2) * g(1) - mrSges(2,3) * t656 + Ifges(2,5) * qJDD(1) - t677 * Ifges(2,6) - pkin(6) * t550 + t675 * t538 - t671 * t539;
t536 = Ifges(2,6) * qJDD(1) + t677 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t655 + t671 * t538 + t675 * t539 - pkin(1) * (-m(3) * g(1) + t555) + pkin(6) * t689;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t555; -m(1) * g(2) + t692; -m(1) * g(3) + t690; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t680; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t690 + t676 * t536 + t672 * t537; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t692 + t672 * t536 - t676 * t537; t680; t682; t698; t580; -t681;];
tauJB = t1;
