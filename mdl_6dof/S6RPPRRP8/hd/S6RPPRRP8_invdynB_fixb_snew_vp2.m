% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:06:33
% EndTime: 2019-05-05 15:06:37
% DurationCPUTime: 4.29s
% Computational Cost: add. (40718->290), mult. (90915->340), div. (0->0), fcn. (61518->8), ass. (0->122)
t692 = Ifges(6,1) + Ifges(7,1);
t681 = Ifges(6,4) - Ifges(7,5);
t691 = -Ifges(6,5) - Ifges(7,4);
t690 = Ifges(6,2) + Ifges(7,3);
t678 = Ifges(6,6) - Ifges(7,6);
t689 = -Ifges(6,3) - Ifges(7,2);
t638 = sin(qJ(1));
t640 = cos(qJ(1));
t616 = t638 * g(1) - t640 * g(2);
t642 = qJD(1) ^ 2;
t649 = -t642 * qJ(2) + qJDD(2) - t616;
t677 = -pkin(1) - qJ(3);
t688 = -(2 * qJD(1) * qJD(3)) + t677 * qJDD(1) + t649;
t634 = sin(pkin(9));
t628 = t634 ^ 2;
t635 = cos(pkin(9));
t670 = t635 ^ 2 + t628;
t662 = t670 * mrSges(4,3);
t617 = -t640 * g(1) - t638 * g(2);
t687 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t617;
t637 = sin(qJ(4));
t639 = cos(qJ(4));
t653 = t634 * t639 + t635 * t637;
t612 = t653 * qJD(1);
t652 = -t634 * t637 + t635 * t639;
t613 = t652 * qJD(1);
t667 = t613 * qJD(4);
t597 = -t653 * qJDD(1) - t667;
t686 = cos(qJ(5));
t685 = pkin(3) * t642;
t684 = mrSges(2,1) - mrSges(3,2);
t683 = -mrSges(6,3) - mrSges(7,2);
t682 = -Ifges(3,4) + Ifges(2,5);
t679 = -Ifges(2,6) + Ifges(3,5);
t676 = mrSges(4,2) * t635;
t594 = t634 * g(3) + t635 * t688;
t577 = (-pkin(7) * qJDD(1) - t634 * t685) * t635 + t594;
t595 = -g(3) * t635 + t634 * t688;
t665 = qJDD(1) * t634;
t578 = -pkin(7) * t665 - t628 * t685 + t595;
t555 = t637 * t577 + t639 * t578;
t591 = mrSges(5,1) * t612 + mrSges(5,2) * t613;
t606 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t613;
t596 = pkin(4) * t612 - pkin(8) * t613;
t641 = qJD(4) ^ 2;
t551 = -pkin(4) * t641 + qJDD(4) * pkin(8) - t596 * t612 + t555;
t648 = qJDD(3) + t687;
t584 = pkin(3) * t665 + (-t670 * pkin(7) + t677) * t642 + t648;
t668 = t612 * qJD(4);
t598 = t652 * qJDD(1) - t668;
t553 = (-t598 + t668) * pkin(8) + (-t597 + t667) * pkin(4) + t584;
t636 = sin(qJ(5));
t548 = t686 * t551 + t636 * t553;
t601 = t636 * qJD(4) + t686 * t613;
t566 = t601 * qJD(5) - t686 * qJDD(4) + t636 * t598;
t610 = qJD(5) + t612;
t581 = mrSges(6,1) * t610 - mrSges(6,3) * t601;
t593 = qJDD(5) - t597;
t600 = -t686 * qJD(4) + t636 * t613;
t571 = pkin(5) * t600 - qJ(6) * t601;
t609 = t610 ^ 2;
t544 = -pkin(5) * t609 + qJ(6) * t593 + 0.2e1 * qJD(6) * t610 - t571 * t600 + t548;
t582 = -mrSges(7,1) * t610 + mrSges(7,2) * t601;
t663 = m(7) * t544 + t593 * mrSges(7,3) + t610 * t582;
t572 = mrSges(7,1) * t600 - mrSges(7,3) * t601;
t671 = -mrSges(6,1) * t600 - mrSges(6,2) * t601 - t572;
t539 = m(6) * t548 - t593 * mrSges(6,2) + t683 * t566 - t610 * t581 + t671 * t600 + t663;
t547 = -t636 * t551 + t686 * t553;
t567 = -t600 * qJD(5) + t636 * qJDD(4) + t686 * t598;
t580 = -mrSges(6,2) * t610 - mrSges(6,3) * t600;
t545 = -t593 * pkin(5) - t609 * qJ(6) + t601 * t571 + qJDD(6) - t547;
t579 = -mrSges(7,2) * t600 + mrSges(7,3) * t610;
t657 = -m(7) * t545 + t593 * mrSges(7,1) + t610 * t579;
t541 = m(6) * t547 + t593 * mrSges(6,1) + t683 * t567 + t610 * t580 + t671 * t601 + t657;
t658 = t686 * t539 - t541 * t636;
t532 = m(5) * t555 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t597 - qJD(4) * t606 - t591 * t612 + t658;
t554 = t639 * t577 - t637 * t578;
t605 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t612;
t550 = -qJDD(4) * pkin(4) - t641 * pkin(8) + t613 * t596 - t554;
t546 = -0.2e1 * qJD(6) * t601 + (t600 * t610 - t567) * qJ(6) + (t601 * t610 + t566) * pkin(5) + t550;
t542 = m(7) * t546 + mrSges(7,1) * t566 - t567 * mrSges(7,3) + t579 * t600 - t601 * t582;
t643 = -m(6) * t550 - t566 * mrSges(6,1) - mrSges(6,2) * t567 - t600 * t580 - t581 * t601 - t542;
t536 = m(5) * t554 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t598 + qJD(4) * t605 - t591 * t613 + t643;
t525 = t637 * t532 + t639 * t536;
t651 = -qJDD(1) * mrSges(4,3) - t642 * (mrSges(4,1) * t634 + t676);
t523 = m(4) * t594 + t651 * t635 + t525;
t659 = t639 * t532 - t637 * t536;
t524 = m(4) * t595 + t651 * t634 + t659;
t520 = t635 * t523 + t634 * t524;
t611 = -qJDD(1) * pkin(1) + t649;
t647 = -m(3) * t611 + t642 * mrSges(3,3) - t520;
t517 = m(2) * t616 - t642 * mrSges(2,2) + t684 * qJDD(1) + t647;
t608 = t642 * pkin(1) - t687;
t604 = t677 * t642 + t648;
t534 = t636 * t539 + t686 * t541;
t646 = m(5) * t584 - t597 * mrSges(5,1) + t598 * mrSges(5,2) + t612 * t605 + t613 * t606 + t534;
t645 = -m(4) * t604 - mrSges(4,1) * t665 - qJDD(1) * t676 - t646;
t644 = -m(3) * t608 + t642 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t645;
t528 = (-mrSges(2,1) - t662) * t642 + t644 - qJDD(1) * mrSges(2,2) + m(2) * t617;
t675 = t640 * t517 + t638 * t528;
t674 = t600 * t690 - t601 * t681 - t610 * t678;
t673 = t600 * t678 + t601 * t691 + t610 * t689;
t672 = -t681 * t600 + t601 * t692 - t691 * t610;
t654 = Ifges(4,5) * t635 - Ifges(4,6) * t634;
t669 = t642 * t654;
t661 = -t517 * t638 + t640 * t528;
t660 = -t634 * t523 + t635 * t524;
t656 = Ifges(4,1) * t635 - Ifges(4,4) * t634;
t655 = Ifges(4,4) * t635 - Ifges(4,2) * t634;
t587 = Ifges(5,1) * t613 - Ifges(5,4) * t612 + Ifges(5,5) * qJD(4);
t586 = Ifges(5,4) * t613 - Ifges(5,2) * t612 + Ifges(5,6) * qJD(4);
t585 = Ifges(5,5) * t613 - Ifges(5,6) * t612 + Ifges(5,3) * qJD(4);
t533 = mrSges(6,2) * t550 + mrSges(7,2) * t545 - mrSges(6,3) * t547 - mrSges(7,3) * t546 - qJ(6) * t542 - t681 * t566 + t567 * t692 - t593 * t691 + t673 * t600 + t674 * t610;
t529 = -mrSges(6,1) * t550 - mrSges(7,1) * t546 + mrSges(7,2) * t544 + mrSges(6,3) * t548 - pkin(5) * t542 - t566 * t690 + t681 * t567 + t678 * t593 + t673 * t601 + t672 * t610;
t521 = Ifges(5,4) * t598 + Ifges(5,2) * t597 + Ifges(5,6) * qJDD(4) - t613 * t585 + qJD(4) * t587 - mrSges(5,1) * t584 + mrSges(5,3) * t555 - mrSges(6,1) * t547 + mrSges(6,2) * t548 + mrSges(7,1) * t545 - mrSges(7,3) * t544 - pkin(5) * t657 - qJ(6) * t663 - pkin(4) * t534 + (pkin(5) * t572 + t674) * t601 + (qJ(6) * t572 - t672) * t600 + t689 * t593 + (mrSges(7,2) * pkin(5) + t691) * t567 + (mrSges(7,2) * qJ(6) + t678) * t566;
t519 = -m(3) * g(3) + t660;
t518 = mrSges(5,2) * t584 - mrSges(5,3) * t554 + Ifges(5,1) * t598 + Ifges(5,4) * t597 + Ifges(5,5) * qJDD(4) - pkin(8) * t534 - qJD(4) * t586 - t636 * t529 + t686 * t533 - t612 * t585;
t515 = mrSges(4,2) * t604 - mrSges(4,3) * t594 - pkin(7) * t525 + t656 * qJDD(1) + t639 * t518 - t637 * t521 - t634 * t669;
t514 = -mrSges(4,1) * t604 + mrSges(4,3) * t595 - pkin(3) * t646 + pkin(7) * t659 + t655 * qJDD(1) + t637 * t518 + t639 * t521 - t635 * t669;
t513 = t686 * t529 + Ifges(5,3) * qJDD(4) + pkin(8) * t658 + t636 * t533 - mrSges(2,3) * t616 + mrSges(3,1) * t611 + t612 * t587 + t613 * t586 + pkin(4) * t643 + Ifges(5,5) * t598 + mrSges(4,1) * t594 - mrSges(4,2) * t595 + Ifges(5,6) * t597 - mrSges(5,2) * t555 + mrSges(5,1) * t554 + pkin(3) * t525 + pkin(2) * t520 - qJ(2) * t519 + (mrSges(3,3) - mrSges(2,2)) * g(3) + (t654 + t682) * qJDD(1) + (t634 * t656 + t635 * t655 + t679) * t642;
t512 = mrSges(2,3) * t617 - mrSges(3,1) * t608 - t634 * t515 - t635 * t514 - pkin(2) * t645 - qJ(3) * t660 - pkin(1) * t519 - t679 * qJDD(1) + t684 * g(3) + (-pkin(2) * t662 + t682) * t642;
t1 = [-m(1) * g(1) + t661; -m(1) * g(2) + t675; (-m(1) - m(2) - m(3)) * g(3) + t660; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t675 - t638 * t512 + t640 * t513; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t661 + t640 * t512 + t638 * t513; pkin(1) * t647 + qJ(2) * (-t642 * t662 + t644) + t635 * t515 - t634 * t514 - qJ(3) * t520 + mrSges(2,1) * t616 - mrSges(2,2) * t617 + mrSges(3,2) * t611 - mrSges(3,3) * t608 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
