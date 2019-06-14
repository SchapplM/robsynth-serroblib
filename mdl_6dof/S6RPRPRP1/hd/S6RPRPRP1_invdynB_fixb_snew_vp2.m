% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:28:48
% EndTime: 2019-05-05 17:28:58
% DurationCPUTime: 7.56s
% Computational Cost: add. (85412->321), mult. (183513->391), div. (0->0), fcn. (119185->10), ass. (0->126)
t686 = -2 * qJD(4);
t685 = Ifges(6,1) + Ifges(7,1);
t680 = Ifges(6,4) + Ifges(7,4);
t679 = Ifges(6,5) + Ifges(7,5);
t684 = Ifges(6,2) + Ifges(7,2);
t683 = Ifges(6,6) + Ifges(7,6);
t682 = Ifges(6,3) + Ifges(7,3);
t647 = sin(qJ(1));
t650 = cos(qJ(1));
t633 = t647 * g(1) - g(2) * t650;
t625 = qJDD(1) * pkin(1) + t633;
t634 = -g(1) * t650 - g(2) * t647;
t652 = qJD(1) ^ 2;
t627 = -pkin(1) * t652 + t634;
t642 = sin(pkin(9));
t644 = cos(pkin(9));
t603 = t642 * t625 + t644 * t627;
t595 = -pkin(2) * t652 + qJDD(1) * pkin(7) + t603;
t640 = -g(3) + qJDD(2);
t646 = sin(qJ(3));
t649 = cos(qJ(3));
t583 = -t646 * t595 + t649 * t640;
t669 = qJD(1) * qJD(3);
t665 = t649 * t669;
t628 = qJDD(1) * t646 + t665;
t562 = (-t628 + t665) * qJ(4) + (t646 * t649 * t652 + qJDD(3)) * pkin(3) + t583;
t584 = t649 * t595 + t646 * t640;
t629 = qJDD(1) * t649 - t646 * t669;
t672 = qJD(1) * t646;
t630 = qJD(3) * pkin(3) - qJ(4) * t672;
t639 = t649 ^ 2;
t563 = -pkin(3) * t639 * t652 + qJ(4) * t629 - qJD(3) * t630 + t584;
t641 = sin(pkin(10));
t643 = cos(pkin(10));
t615 = (t641 * t649 + t643 * t646) * qJD(1);
t554 = t562 * t643 - t641 * t563 + t615 * t686;
t614 = (t641 * t646 - t643 * t649) * qJD(1);
t681 = -mrSges(6,2) - mrSges(7,2);
t555 = t641 * t562 + t643 * t563 + t614 * t686;
t597 = mrSges(5,1) * t614 + mrSges(5,2) * t615;
t604 = -t628 * t641 + t629 * t643;
t610 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t615;
t598 = pkin(4) * t614 - pkin(8) * t615;
t651 = qJD(3) ^ 2;
t553 = -pkin(4) * t651 + qJDD(3) * pkin(8) - t598 * t614 + t555;
t602 = t644 * t625 - t642 * t627;
t656 = -qJDD(1) * pkin(2) - t602;
t565 = -t629 * pkin(3) + qJDD(4) + t630 * t672 + (-qJ(4) * t639 - pkin(7)) * t652 + t656;
t605 = t628 * t643 + t629 * t641;
t558 = (qJD(3) * t614 - t605) * pkin(8) + (qJD(3) * t615 - t604) * pkin(4) + t565;
t645 = sin(qJ(5));
t648 = cos(qJ(5));
t548 = -t645 * t553 + t648 * t558;
t607 = qJD(3) * t648 - t615 * t645;
t578 = qJD(5) * t607 + qJDD(3) * t645 + t605 * t648;
t608 = qJD(3) * t645 + t615 * t648;
t580 = -mrSges(7,1) * t607 + mrSges(7,2) * t608;
t581 = -mrSges(6,1) * t607 + mrSges(6,2) * t608;
t613 = qJD(5) + t614;
t586 = -mrSges(6,2) * t613 + mrSges(6,3) * t607;
t601 = qJDD(5) - t604;
t545 = -0.2e1 * qJD(6) * t608 + (t607 * t613 - t578) * qJ(6) + (t607 * t608 + t601) * pkin(5) + t548;
t585 = -mrSges(7,2) * t613 + mrSges(7,3) * t607;
t668 = m(7) * t545 + t601 * mrSges(7,1) + t613 * t585;
t537 = m(6) * t548 + t601 * mrSges(6,1) + t613 * t586 + (-t580 - t581) * t608 + (-mrSges(6,3) - mrSges(7,3)) * t578 + t668;
t549 = t648 * t553 + t645 * t558;
t577 = -qJD(5) * t608 + qJDD(3) * t648 - t605 * t645;
t587 = pkin(5) * t613 - qJ(6) * t608;
t606 = t607 ^ 2;
t547 = -pkin(5) * t606 + qJ(6) * t577 + 0.2e1 * qJD(6) * t607 - t587 * t613 + t549;
t667 = m(7) * t547 + t577 * mrSges(7,3) + t607 * t580;
t588 = mrSges(7,1) * t613 - mrSges(7,3) * t608;
t673 = -mrSges(6,1) * t613 + mrSges(6,3) * t608 - t588;
t540 = m(6) * t549 + t577 * mrSges(6,3) + t607 * t581 + t601 * t681 + t673 * t613 + t667;
t660 = -t537 * t645 + t648 * t540;
t533 = m(5) * t555 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t604 - qJD(3) * t610 - t597 * t614 + t660;
t609 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t614;
t552 = -qJDD(3) * pkin(4) - pkin(8) * t651 + t615 * t598 - t554;
t550 = -pkin(5) * t577 - qJ(6) * t606 + t587 * t608 + qJDD(6) + t552;
t658 = m(7) * t550 - t577 * mrSges(7,1) - t607 * t585;
t654 = -m(6) * t552 + t577 * mrSges(6,1) + t578 * t681 + t607 * t586 + t673 * t608 - t658;
t542 = m(5) * t554 + qJDD(3) * mrSges(5,1) - t605 * mrSges(5,3) + qJD(3) * t609 - t615 * t597 + t654;
t527 = t641 * t533 + t643 * t542;
t626 = (-mrSges(4,1) * t649 + mrSges(4,2) * t646) * qJD(1);
t671 = qJD(1) * t649;
t632 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t671;
t525 = m(4) * t583 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t628 + qJD(3) * t632 - t626 * t672 + t527;
t631 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t672;
t661 = t643 * t533 - t542 * t641;
t526 = m(4) * t584 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t629 - qJD(3) * t631 + t626 * t671 + t661;
t662 = -t525 * t646 + t649 * t526;
t519 = m(3) * t603 - mrSges(3,1) * t652 - qJDD(1) * mrSges(3,2) + t662;
t594 = -t652 * pkin(7) + t656;
t535 = t648 * t537 + t645 * t540;
t655 = m(5) * t565 - t604 * mrSges(5,1) + mrSges(5,2) * t605 + t614 * t609 + t610 * t615 + t535;
t653 = -m(4) * t594 + t629 * mrSges(4,1) - mrSges(4,2) * t628 - t631 * t672 + t632 * t671 - t655;
t530 = m(3) * t602 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t652 + t653;
t515 = t642 * t519 + t644 * t530;
t513 = m(2) * t633 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t652 + t515;
t663 = t644 * t519 - t530 * t642;
t514 = m(2) * t634 - mrSges(2,1) * t652 - qJDD(1) * mrSges(2,2) + t663;
t677 = t650 * t513 + t647 * t514;
t520 = t649 * t525 + t646 * t526;
t676 = t683 * t607 + t679 * t608 + t682 * t613;
t675 = -t684 * t607 - t680 * t608 - t683 * t613;
t674 = t680 * t607 + t685 * t608 + t679 * t613;
t666 = m(3) * t640 + t520;
t664 = -t513 * t647 + t650 * t514;
t621 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t646 + Ifges(4,4) * t649) * qJD(1);
t620 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t646 + Ifges(4,2) * t649) * qJD(1);
t619 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t646 + Ifges(4,6) * t649) * qJD(1);
t593 = Ifges(5,1) * t615 - Ifges(5,4) * t614 + Ifges(5,5) * qJD(3);
t592 = Ifges(5,4) * t615 - Ifges(5,2) * t614 + Ifges(5,6) * qJD(3);
t591 = Ifges(5,5) * t615 - Ifges(5,6) * t614 + Ifges(5,3) * qJD(3);
t543 = -t578 * mrSges(7,3) - t608 * t580 + t668;
t534 = mrSges(6,2) * t552 + mrSges(7,2) * t550 - mrSges(6,3) * t548 - mrSges(7,3) * t545 - qJ(6) * t543 + t680 * t577 + t685 * t578 + t679 * t601 + t676 * t607 + t675 * t613;
t528 = -mrSges(6,1) * t552 + mrSges(6,3) * t549 - mrSges(7,1) * t550 + mrSges(7,3) * t547 - pkin(5) * t658 + qJ(6) * t667 + (-qJ(6) * t588 + t674) * t613 + (-pkin(5) * t588 - t676) * t608 + (-mrSges(7,2) * qJ(6) + t683) * t601 + (-mrSges(7,2) * pkin(5) + t680) * t578 + t684 * t577;
t521 = -mrSges(5,1) * t565 - mrSges(6,1) * t548 - mrSges(7,1) * t545 + mrSges(6,2) * t549 + mrSges(7,2) * t547 + mrSges(5,3) * t555 + Ifges(5,4) * t605 + Ifges(5,2) * t604 + Ifges(5,6) * qJDD(3) - pkin(4) * t535 - pkin(5) * t543 + qJD(3) * t593 - t615 * t591 + t675 * t608 + t674 * t607 - t682 * t601 - t679 * t578 - t683 * t577;
t516 = mrSges(5,2) * t565 - mrSges(5,3) * t554 + Ifges(5,1) * t605 + Ifges(5,4) * t604 + Ifges(5,5) * qJDD(3) - pkin(8) * t535 - qJD(3) * t592 - t528 * t645 + t534 * t648 - t591 * t614;
t509 = mrSges(4,2) * t594 - mrSges(4,3) * t583 + Ifges(4,1) * t628 + Ifges(4,4) * t629 + Ifges(4,5) * qJDD(3) - qJ(4) * t527 - qJD(3) * t620 + t516 * t643 - t521 * t641 + t619 * t671;
t508 = -mrSges(4,1) * t594 + mrSges(4,3) * t584 + Ifges(4,4) * t628 + Ifges(4,2) * t629 + Ifges(4,6) * qJDD(3) - pkin(3) * t655 + qJ(4) * t661 + qJD(3) * t621 + t641 * t516 + t643 * t521 - t619 * t672;
t507 = -pkin(2) * t520 - mrSges(3,1) * t640 + mrSges(3,3) * t603 - pkin(3) * t527 - mrSges(4,1) * t583 + mrSges(4,2) * t584 - Ifges(4,5) * t628 - Ifges(4,6) * t629 - t648 * t528 - pkin(4) * t654 - pkin(8) * t660 - Ifges(5,5) * t605 - Ifges(5,6) * t604 - mrSges(5,1) * t554 + mrSges(5,2) * t555 - t645 * t534 + t652 * Ifges(3,5) - t615 * t592 - t614 * t593 + Ifges(3,6) * qJDD(1) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t620 * t646 + t621 * t649) * qJD(1);
t506 = mrSges(3,2) * t640 - mrSges(3,3) * t602 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t652 - pkin(7) * t520 - t508 * t646 + t509 * t649;
t505 = -mrSges(2,2) * g(3) - mrSges(2,3) * t633 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t652 - qJ(2) * t515 + t506 * t644 - t507 * t642;
t504 = mrSges(2,1) * g(3) + mrSges(2,3) * t634 + t652 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t666 + qJ(2) * t663 + t642 * t506 + t644 * t507;
t1 = [-m(1) * g(1) + t664; -m(1) * g(2) + t677; (-m(1) - m(2)) * g(3) + t666; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t677 - t647 * t504 + t650 * t505; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t664 + t650 * t504 + t647 * t505; pkin(1) * t515 + mrSges(2,1) * t633 - mrSges(2,2) * t634 + pkin(7) * t662 + t646 * t509 + t649 * t508 + pkin(2) * t653 + mrSges(3,1) * t602 - mrSges(3,2) * t603 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
