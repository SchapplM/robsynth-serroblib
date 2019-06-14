% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:12:49
% EndTime: 2019-05-05 19:12:57
% DurationCPUTime: 7.29s
% Computational Cost: add. (111189->338), mult. (247225->418), div. (0->0), fcn. (172174->10), ass. (0->130)
t639 = sin(qJ(1));
t643 = cos(qJ(1));
t622 = -t643 * g(1) - t639 * g(2);
t652 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t622;
t668 = -pkin(1) - pkin(7);
t667 = mrSges(2,1) - mrSges(3,2);
t666 = Ifges(2,5) - Ifges(3,4);
t665 = (-Ifges(2,6) + Ifges(3,5));
t621 = t639 * g(1) - t643 * g(2);
t644 = qJD(1) ^ 2;
t651 = -t644 * qJ(2) + qJDD(2) - t621;
t600 = qJDD(1) * t668 + t651;
t638 = sin(qJ(3));
t642 = cos(qJ(3));
t591 = t638 * g(3) + t642 * t600;
t661 = qJD(1) * qJD(3);
t659 = t638 * t661;
t617 = t642 * qJDD(1) - t659;
t571 = (-t617 - t659) * qJ(4) + (-t638 * t642 * t644 + qJDD(3)) * pkin(3) + t591;
t592 = -t642 * g(3) + t638 * t600;
t616 = -t638 * qJDD(1) - t642 * t661;
t662 = qJD(1) * t642;
t619 = qJD(3) * pkin(3) - qJ(4) * t662;
t631 = t638 ^ 2;
t572 = -t631 * t644 * pkin(3) + t616 * qJ(4) - qJD(3) * t619 + t592;
t634 = sin(pkin(10));
t635 = cos(pkin(10));
t607 = (-t634 * t638 + t635 * t642) * qJD(1);
t551 = -0.2e1 * qJD(4) * t607 + t635 * t571 - t634 * t572;
t590 = t634 * t616 + t635 * t617;
t606 = (-t634 * t642 - t635 * t638) * qJD(1);
t541 = (qJD(3) * t606 - t590) * pkin(8) + (t606 * t607 + qJDD(3)) * pkin(4) + t551;
t552 = 0.2e1 * qJD(4) * t606 + t634 * t571 + t635 * t572;
t589 = t635 * t616 - t634 * t617;
t599 = qJD(3) * pkin(4) - t607 * pkin(8);
t605 = t606 ^ 2;
t543 = -t605 * pkin(4) + t589 * pkin(8) - qJD(3) * t599 + t552;
t637 = sin(qJ(5));
t641 = cos(qJ(5));
t538 = t637 * t541 + t641 * t543;
t584 = t637 * t606 + t641 * t607;
t557 = -t584 * qJD(5) + t641 * t589 - t637 * t590;
t583 = t641 * t606 - t637 * t607;
t566 = -t583 * mrSges(6,1) + t584 * mrSges(6,2);
t628 = qJD(3) + qJD(5);
t578 = t628 * mrSges(6,1) - t584 * mrSges(6,3);
t627 = qJDD(3) + qJDD(5);
t567 = -t583 * pkin(5) - t584 * pkin(9);
t626 = t628 ^ 2;
t536 = -t626 * pkin(5) + t627 * pkin(9) + t583 * t567 + t538;
t574 = -t616 * pkin(3) + qJDD(4) + t619 * t662 + (-qJ(4) * t631 + t668) * t644 + t652;
t554 = -t589 * pkin(4) - t605 * pkin(8) + t607 * t599 + t574;
t558 = t583 * qJD(5) + t637 * t589 + t641 * t590;
t539 = t554 + (-t583 * t628 - t558) * pkin(9) + (t584 * t628 - t557) * pkin(5);
t636 = sin(qJ(6));
t640 = cos(qJ(6));
t533 = -t636 * t536 + t640 * t539;
t575 = -t636 * t584 + t640 * t628;
t546 = t575 * qJD(6) + t640 * t558 + t636 * t627;
t556 = qJDD(6) - t557;
t576 = t640 * t584 + t636 * t628;
t559 = -t575 * mrSges(7,1) + t576 * mrSges(7,2);
t579 = qJD(6) - t583;
t560 = -t579 * mrSges(7,2) + t575 * mrSges(7,3);
t531 = m(7) * t533 + t556 * mrSges(7,1) - t546 * mrSges(7,3) - t576 * t559 + t579 * t560;
t534 = t640 * t536 + t636 * t539;
t545 = -t576 * qJD(6) - t636 * t558 + t640 * t627;
t561 = t579 * mrSges(7,1) - t576 * mrSges(7,3);
t532 = m(7) * t534 - t556 * mrSges(7,2) + t545 * mrSges(7,3) + t575 * t559 - t579 * t561;
t654 = -t636 * t531 + t640 * t532;
t522 = m(6) * t538 - t627 * mrSges(6,2) + t557 * mrSges(6,3) + t583 * t566 - t628 * t578 + t654;
t537 = t641 * t541 - t637 * t543;
t577 = -t628 * mrSges(6,2) + t583 * mrSges(6,3);
t535 = -t627 * pkin(5) - t626 * pkin(9) + t584 * t567 - t537;
t649 = -m(7) * t535 + t545 * mrSges(7,1) - t546 * mrSges(7,2) + t575 * t560 - t576 * t561;
t527 = m(6) * t537 + t627 * mrSges(6,1) - t558 * mrSges(6,3) - t584 * t566 + t628 * t577 + t649;
t516 = t637 * t522 + t641 * t527;
t587 = -t606 * mrSges(5,1) + t607 * mrSges(5,2);
t597 = -qJD(3) * mrSges(5,2) + t606 * mrSges(5,3);
t514 = m(5) * t551 + qJDD(3) * mrSges(5,1) - t590 * mrSges(5,3) + qJD(3) * t597 - t607 * t587 + t516;
t598 = qJD(3) * mrSges(5,1) - t607 * mrSges(5,3);
t655 = t641 * t522 - t637 * t527;
t515 = m(5) * t552 - qJDD(3) * mrSges(5,2) + t589 * mrSges(5,3) - qJD(3) * t598 + t606 * t587 + t655;
t508 = t635 * t514 + t634 * t515;
t615 = (mrSges(4,1) * t638 + mrSges(4,2) * t642) * qJD(1);
t663 = qJD(1) * t638;
t618 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t663;
t506 = m(4) * t591 + qJDD(3) * mrSges(4,1) - t617 * mrSges(4,3) + qJD(3) * t618 - t615 * t662 + t508;
t620 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t662;
t656 = -t634 * t514 + t635 * t515;
t507 = m(4) * t592 - qJDD(3) * mrSges(4,2) + t616 * mrSges(4,3) - qJD(3) * t620 - t615 * t663 + t656;
t503 = t642 * t506 + t638 * t507;
t604 = -qJDD(1) * pkin(1) + t651;
t650 = -m(3) * t604 + (t644 * mrSges(3,3)) - t503;
t501 = m(2) * t621 - (t644 * mrSges(2,2)) + qJDD(1) * t667 + t650;
t601 = t644 * pkin(1) - t652;
t596 = t644 * t668 + t652;
t523 = t640 * t531 + t636 * t532;
t648 = m(6) * t554 - t557 * mrSges(6,1) + t558 * mrSges(6,2) - t583 * t577 + t584 * t578 + t523;
t647 = m(5) * t574 - t589 * mrSges(5,1) + t590 * mrSges(5,2) - t606 * t597 + t607 * t598 + t648;
t646 = -m(4) * t596 + t616 * mrSges(4,1) - t617 * mrSges(4,2) - t618 * t663 - t620 * t662 - t647;
t645 = -m(3) * t601 + (t644 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t646;
t519 = m(2) * t622 - (t644 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t645;
t664 = t643 * t501 + t639 * t519;
t658 = -t639 * t501 + t643 * t519;
t657 = -t638 * t506 + t642 * t507;
t610 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t642 - Ifges(4,4) * t638) * qJD(1);
t609 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t642 - Ifges(4,2) * t638) * qJD(1);
t608 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t642 - Ifges(4,6) * t638) * qJD(1);
t582 = Ifges(5,1) * t607 + Ifges(5,4) * t606 + (Ifges(5,5) * qJD(3));
t581 = Ifges(5,4) * t607 + Ifges(5,2) * t606 + (Ifges(5,6) * qJD(3));
t580 = Ifges(5,5) * t607 + Ifges(5,6) * t606 + (Ifges(5,3) * qJD(3));
t564 = Ifges(6,1) * t584 + Ifges(6,4) * t583 + Ifges(6,5) * t628;
t563 = Ifges(6,4) * t584 + Ifges(6,2) * t583 + Ifges(6,6) * t628;
t562 = Ifges(6,5) * t584 + Ifges(6,6) * t583 + Ifges(6,3) * t628;
t549 = Ifges(7,1) * t576 + Ifges(7,4) * t575 + Ifges(7,5) * t579;
t548 = Ifges(7,4) * t576 + Ifges(7,2) * t575 + Ifges(7,6) * t579;
t547 = Ifges(7,5) * t576 + Ifges(7,6) * t575 + Ifges(7,3) * t579;
t525 = mrSges(7,2) * t535 - mrSges(7,3) * t533 + Ifges(7,1) * t546 + Ifges(7,4) * t545 + Ifges(7,5) * t556 + t575 * t547 - t579 * t548;
t524 = -mrSges(7,1) * t535 + mrSges(7,3) * t534 + Ifges(7,4) * t546 + Ifges(7,2) * t545 + Ifges(7,6) * t556 - t576 * t547 + t579 * t549;
t510 = -mrSges(6,1) * t554 - mrSges(7,1) * t533 + mrSges(7,2) * t534 + mrSges(6,3) * t538 + Ifges(6,4) * t558 - Ifges(7,5) * t546 + Ifges(6,2) * t557 + Ifges(6,6) * t627 - Ifges(7,6) * t545 - Ifges(7,3) * t556 - pkin(5) * t523 - t576 * t548 + t575 * t549 - t584 * t562 + t628 * t564;
t509 = mrSges(6,2) * t554 - mrSges(6,3) * t537 + Ifges(6,1) * t558 + Ifges(6,4) * t557 + Ifges(6,5) * t627 - pkin(9) * t523 - t636 * t524 + t640 * t525 + t583 * t562 - t628 * t563;
t504 = mrSges(5,2) * t574 - mrSges(5,3) * t551 + Ifges(5,1) * t590 + Ifges(5,4) * t589 + Ifges(5,5) * qJDD(3) - pkin(8) * t516 - qJD(3) * t581 + t641 * t509 - t637 * t510 + t606 * t580;
t502 = -m(3) * g(3) + t657;
t499 = -mrSges(5,1) * t574 + mrSges(5,3) * t552 + Ifges(5,4) * t590 + Ifges(5,2) * t589 + Ifges(5,6) * qJDD(3) - pkin(4) * t648 + pkin(8) * t655 + qJD(3) * t582 + t637 * t509 + t641 * t510 - t607 * t580;
t498 = mrSges(4,2) * t596 - mrSges(4,3) * t591 + Ifges(4,1) * t617 + Ifges(4,4) * t616 + Ifges(4,5) * qJDD(3) - qJ(4) * t508 - qJD(3) * t609 - t634 * t499 + t635 * t504 - t608 * t663;
t497 = -mrSges(4,1) * t596 + mrSges(4,3) * t592 + Ifges(4,4) * t617 + Ifges(4,2) * t616 + Ifges(4,6) * qJDD(3) - pkin(3) * t647 + qJ(4) * t656 + qJD(3) * t610 + t635 * t499 + t634 * t504 - t608 * t662;
t496 = pkin(4) * t516 + pkin(5) * t649 + (t665 * t644) + t666 * qJDD(1) + (t642 * t609 + t638 * t610) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t640 * t524 + t636 * t525 + Ifges(6,3) * t627 + Ifges(4,6) * t616 + Ifges(4,5) * t617 - mrSges(2,3) * t621 + mrSges(3,1) * t604 - t606 * t582 + t607 * t581 + Ifges(5,6) * t589 + Ifges(5,5) * t590 + mrSges(4,1) * t591 - mrSges(4,2) * t592 - t583 * t564 + t584 * t563 + Ifges(6,6) * t557 + Ifges(6,5) * t558 + mrSges(5,1) * t551 - mrSges(5,2) * t552 + pkin(9) * t654 + pkin(3) * t508 + pkin(2) * t503 - qJ(2) * t502 - mrSges(6,2) * t538 + mrSges(6,1) * t537;
t495 = -mrSges(3,1) * t601 + mrSges(2,3) * t622 - pkin(1) * t502 - pkin(2) * t646 - pkin(7) * t657 + g(3) * t667 - qJDD(1) * t665 - t642 * t497 - t638 * t498 + t644 * t666;
t1 = [-m(1) * g(1) + t658; -m(1) * g(2) + t664; (-m(1) - m(2) - m(3)) * g(3) + t657; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t664 - t639 * t495 + t643 * t496; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t658 + t643 * t495 + t639 * t496; pkin(1) * t650 + qJ(2) * t645 + t642 * t498 - t638 * t497 - pkin(7) * t503 + mrSges(2,1) * t621 - mrSges(2,2) * t622 + mrSges(3,2) * t604 - mrSges(3,3) * t601 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
