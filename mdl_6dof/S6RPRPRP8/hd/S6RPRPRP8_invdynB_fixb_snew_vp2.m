% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:00:48
% EndTime: 2019-05-05 18:00:54
% DurationCPUTime: 4.52s
% Computational Cost: add. (44757->314), mult. (96671->374), div. (0->0), fcn. (61670->8), ass. (0->121)
t685 = -2 * qJD(4);
t684 = Ifges(6,1) + Ifges(7,1);
t676 = Ifges(6,4) - Ifges(7,5);
t683 = -Ifges(6,5) - Ifges(7,4);
t682 = Ifges(6,2) + Ifges(7,3);
t672 = Ifges(6,6) - Ifges(7,6);
t681 = -Ifges(6,3) - Ifges(7,2);
t640 = sin(qJ(1));
t642 = cos(qJ(1));
t623 = t640 * g(1) - t642 * g(2);
t644 = qJD(1) ^ 2;
t650 = -t644 * qJ(2) + qJDD(2) - t623;
t680 = -pkin(1) - pkin(7);
t601 = qJDD(1) * t680 + t650;
t639 = sin(qJ(3));
t641 = cos(qJ(3));
t591 = t639 * g(3) + t641 * t601;
t663 = qJD(1) * qJD(3);
t660 = t639 * t663;
t619 = t641 * qJDD(1) - t660;
t567 = (-t619 - t660) * qJ(4) + (-t639 * t641 * t644 + qJDD(3)) * pkin(3) + t591;
t592 = -t641 * g(3) + t639 * t601;
t618 = -t639 * qJDD(1) - t641 * t663;
t665 = qJD(1) * t641;
t621 = qJD(3) * pkin(3) - qJ(4) * t665;
t633 = t639 ^ 2;
t568 = -t633 * t644 * pkin(3) + t618 * qJ(4) - qJD(3) * t621 + t592;
t636 = sin(pkin(9));
t637 = cos(pkin(9));
t666 = qJD(1) * t639;
t609 = -t636 * t666 + t637 * t665;
t548 = t637 * t567 - t636 * t568 + t609 * t685;
t624 = -t642 * g(1) - t640 * g(2);
t651 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t624;
t608 = (t636 * t641 + t637 * t639) * qJD(1);
t679 = cos(qJ(5));
t678 = mrSges(2,1) - mrSges(3,2);
t677 = -mrSges(6,3) - mrSges(7,2);
t674 = Ifges(2,5) - Ifges(3,4);
t673 = -Ifges(2,6) + Ifges(3,5);
t549 = t636 * t567 + t637 * t568 + t608 * t685;
t583 = t608 * mrSges(5,1) + t609 * mrSges(5,2);
t589 = t637 * t618 - t636 * t619;
t600 = qJD(3) * mrSges(5,1) - t609 * mrSges(5,3);
t584 = t608 * pkin(4) - t609 * pkin(8);
t643 = qJD(3) ^ 2;
t545 = -t643 * pkin(4) + qJDD(3) * pkin(8) - t608 * t584 + t549;
t570 = -t618 * pkin(3) + qJDD(4) + t621 * t665 + (-qJ(4) * t633 + t680) * t644 + t651;
t590 = t636 * t618 + t637 * t619;
t547 = (qJD(3) * t608 - t590) * pkin(8) + (qJD(3) * t609 - t589) * pkin(4) + t570;
t638 = sin(qJ(5));
t542 = t679 * t545 + t638 * t547;
t594 = t638 * qJD(3) + t609 * t679;
t560 = t594 * qJD(5) - qJDD(3) * t679 + t638 * t590;
t606 = qJD(5) + t608;
t577 = t606 * mrSges(6,1) - t594 * mrSges(6,3);
t588 = qJDD(5) - t589;
t593 = -qJD(3) * t679 + t638 * t609;
t571 = t593 * pkin(5) - t594 * qJ(6);
t605 = t606 ^ 2;
t538 = -t605 * pkin(5) + t588 * qJ(6) + 0.2e1 * qJD(6) * t606 - t593 * t571 + t542;
t578 = -t606 * mrSges(7,1) + t594 * mrSges(7,2);
t661 = m(7) * t538 + t588 * mrSges(7,3) + t606 * t578;
t572 = t593 * mrSges(7,1) - t594 * mrSges(7,3);
t667 = -t593 * mrSges(6,1) - t594 * mrSges(6,2) - t572;
t533 = m(6) * t542 - t588 * mrSges(6,2) + t560 * t677 - t606 * t577 + t593 * t667 + t661;
t541 = -t638 * t545 + t547 * t679;
t561 = -t593 * qJD(5) + t638 * qJDD(3) + t590 * t679;
t576 = -t606 * mrSges(6,2) - t593 * mrSges(6,3);
t539 = -t588 * pkin(5) - t605 * qJ(6) + t594 * t571 + qJDD(6) - t541;
t575 = -t593 * mrSges(7,2) + t606 * mrSges(7,3);
t654 = -m(7) * t539 + t588 * mrSges(7,1) + t606 * t575;
t535 = m(6) * t541 + t588 * mrSges(6,1) + t561 * t677 + t606 * t576 + t594 * t667 + t654;
t656 = t679 * t533 - t638 * t535;
t525 = m(5) * t549 - qJDD(3) * mrSges(5,2) + t589 * mrSges(5,3) - qJD(3) * t600 - t608 * t583 + t656;
t599 = -qJD(3) * mrSges(5,2) - t608 * mrSges(5,3);
t544 = -qJDD(3) * pkin(4) - t643 * pkin(8) + t609 * t584 - t548;
t540 = -0.2e1 * qJD(6) * t594 + (t593 * t606 - t561) * qJ(6) + (t594 * t606 + t560) * pkin(5) + t544;
t536 = m(7) * t540 + t560 * mrSges(7,1) - t561 * mrSges(7,3) + t593 * t575 - t594 * t578;
t646 = -m(6) * t544 - t560 * mrSges(6,1) - t561 * mrSges(6,2) - t593 * t576 - t594 * t577 - t536;
t530 = m(5) * t548 + qJDD(3) * mrSges(5,1) - t590 * mrSges(5,3) + qJD(3) * t599 - t609 * t583 + t646;
t519 = t636 * t525 + t637 * t530;
t617 = (mrSges(4,1) * t639 + mrSges(4,2) * t641) * qJD(1);
t620 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t666;
t517 = m(4) * t591 + qJDD(3) * mrSges(4,1) - t619 * mrSges(4,3) + qJD(3) * t620 - t617 * t665 + t519;
t622 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t665;
t657 = t637 * t525 - t636 * t530;
t518 = m(4) * t592 - qJDD(3) * mrSges(4,2) + t618 * mrSges(4,3) - qJD(3) * t622 - t617 * t666 + t657;
t513 = t641 * t517 + t639 * t518;
t607 = -qJDD(1) * pkin(1) + t650;
t649 = -m(3) * t607 + t644 * mrSges(3,3) - t513;
t511 = m(2) * t623 - t644 * mrSges(2,2) + qJDD(1) * t678 + t649;
t602 = t644 * pkin(1) - t651;
t598 = t644 * t680 + t651;
t528 = t638 * t533 + t679 * t535;
t648 = m(5) * t570 - t589 * mrSges(5,1) + t590 * mrSges(5,2) + t608 * t599 + t609 * t600 + t528;
t647 = -m(4) * t598 + t618 * mrSges(4,1) - t619 * mrSges(4,2) - t620 * t666 - t622 * t665 - t648;
t645 = -m(3) * t602 + t644 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t647;
t522 = m(2) * t624 - t644 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t645;
t671 = t642 * t511 + t640 * t522;
t670 = t682 * t593 - t676 * t594 - t672 * t606;
t669 = t672 * t593 + t683 * t594 + t681 * t606;
t668 = -t676 * t593 + t684 * t594 - t683 * t606;
t659 = -t640 * t511 + t642 * t522;
t658 = -t639 * t517 + t641 * t518;
t612 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t641 - Ifges(4,4) * t639) * qJD(1);
t611 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t641 - Ifges(4,2) * t639) * qJD(1);
t610 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t641 - Ifges(4,6) * t639) * qJD(1);
t581 = Ifges(5,1) * t609 - Ifges(5,4) * t608 + Ifges(5,5) * qJD(3);
t580 = Ifges(5,4) * t609 - Ifges(5,2) * t608 + Ifges(5,6) * qJD(3);
t579 = Ifges(5,5) * t609 - Ifges(5,6) * t608 + Ifges(5,3) * qJD(3);
t527 = mrSges(6,2) * t544 + mrSges(7,2) * t539 - mrSges(6,3) * t541 - mrSges(7,3) * t540 - qJ(6) * t536 - t676 * t560 + t684 * t561 - t588 * t683 + t669 * t593 + t670 * t606;
t526 = -mrSges(6,1) * t544 - mrSges(7,1) * t540 + mrSges(7,2) * t538 + mrSges(6,3) * t542 - pkin(5) * t536 - t682 * t560 + t676 * t561 + t672 * t588 + t669 * t594 + t668 * t606;
t515 = Ifges(5,4) * t590 + Ifges(5,2) * t589 + Ifges(5,6) * qJDD(3) - t609 * t579 + qJD(3) * t581 - mrSges(5,1) * t570 + mrSges(5,3) * t549 - mrSges(6,1) * t541 + mrSges(6,2) * t542 + mrSges(7,1) * t539 - mrSges(7,3) * t538 - pkin(5) * t654 - qJ(6) * t661 - pkin(4) * t528 + (pkin(5) * t572 + t670) * t594 + (qJ(6) * t572 - t668) * t593 + t681 * t588 + (pkin(5) * mrSges(7,2) + t683) * t561 + (qJ(6) * mrSges(7,2) + t672) * t560;
t514 = mrSges(5,2) * t570 - mrSges(5,3) * t548 + Ifges(5,1) * t590 + Ifges(5,4) * t589 + Ifges(5,5) * qJDD(3) - pkin(8) * t528 - qJD(3) * t580 - t638 * t526 + t527 * t679 - t608 * t579;
t512 = -m(3) * g(3) + t658;
t509 = mrSges(4,2) * t598 - mrSges(4,3) * t591 + Ifges(4,1) * t619 + Ifges(4,4) * t618 + Ifges(4,5) * qJDD(3) - qJ(4) * t519 - qJD(3) * t611 + t637 * t514 - t636 * t515 - t610 * t666;
t508 = -mrSges(4,1) * t598 + mrSges(4,3) * t592 + Ifges(4,4) * t619 + Ifges(4,2) * t618 + Ifges(4,6) * qJDD(3) - pkin(3) * t648 + qJ(4) * t657 + qJD(3) * t612 + t636 * t514 + t637 * t515 - t610 * t665;
t507 = -qJ(2) * t512 + pkin(2) * t513 + pkin(8) * t656 + t638 * t527 + mrSges(5,1) * t548 - mrSges(5,2) * t549 + Ifges(5,5) * t590 + mrSges(4,1) * t591 - mrSges(4,2) * t592 + Ifges(4,6) * t618 + Ifges(4,5) * t619 - mrSges(2,3) * t623 + mrSges(3,1) * t607 + t608 * t581 + t609 * t580 + pkin(4) * t646 + pkin(3) * t519 + Ifges(5,6) * t589 + t679 * t526 + t673 * t644 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t674 * qJDD(1) + (t641 * t611 + t639 * t612) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t506 = -mrSges(3,1) * t602 + mrSges(2,3) * t624 - pkin(1) * t512 - pkin(2) * t647 - pkin(7) * t658 + g(3) * t678 - qJDD(1) * t673 - t641 * t508 - t639 * t509 + t644 * t674;
t1 = [-m(1) * g(1) + t659; -m(1) * g(2) + t671; (-m(1) - m(2) - m(3)) * g(3) + t658; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t671 - t640 * t506 + t642 * t507; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t659 + t642 * t506 + t640 * t507; pkin(1) * t649 + qJ(2) * t645 - pkin(7) * t513 + mrSges(2,1) * t623 - mrSges(2,2) * t624 + t641 * t509 - t639 * t508 + mrSges(3,2) * t607 - mrSges(3,3) * t602 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
