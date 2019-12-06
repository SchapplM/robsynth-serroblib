% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:24
% EndTime: 2019-12-05 18:16:27
% DurationCPUTime: 3.20s
% Computational Cost: add. (49986->229), mult. (67543->290), div. (0->0), fcn. (38166->10), ass. (0->100)
t593 = sin(qJ(1));
t597 = cos(qJ(1));
t569 = t597 * g(2) + t593 * g(3);
t562 = qJDD(1) * pkin(1) + t569;
t568 = t593 * g(2) - t597 * g(3);
t598 = qJD(1) ^ 2;
t563 = -t598 * pkin(1) + t568;
t588 = sin(pkin(9));
t589 = cos(pkin(9));
t545 = t589 * t562 - t588 * t563;
t542 = qJDD(1) * pkin(2) + t545;
t546 = t588 * t562 + t589 * t563;
t543 = -t598 * pkin(2) + t546;
t592 = sin(qJ(3));
t596 = cos(qJ(3));
t527 = t592 * t542 + t596 * t543;
t583 = qJD(1) + qJD(3);
t579 = t583 ^ 2;
t581 = qJDD(1) + qJDD(3);
t524 = -t579 * pkin(3) + t581 * pkin(7) + t527;
t587 = -g(1) + qJDD(2);
t591 = sin(qJ(4));
t595 = cos(qJ(4));
t520 = -t591 * t524 + t595 * t587;
t615 = qJD(4) * t583;
t613 = t595 * t615;
t557 = t591 * t581 + t613;
t517 = (-t557 + t613) * pkin(8) + (t579 * t591 * t595 + qJDD(4)) * pkin(4) + t520;
t521 = t595 * t524 + t591 * t587;
t558 = t595 * t581 - t591 * t615;
t617 = t583 * t591;
t566 = qJD(4) * pkin(4) - pkin(8) * t617;
t586 = t595 ^ 2;
t518 = -t586 * t579 * pkin(4) + t558 * pkin(8) - qJD(4) * t566 + t521;
t590 = sin(qJ(5));
t594 = cos(qJ(5));
t515 = t594 * t517 - t590 * t518;
t552 = (-t590 * t591 + t594 * t595) * t583;
t533 = t552 * qJD(5) + t594 * t557 + t590 * t558;
t553 = (t590 * t595 + t591 * t594) * t583;
t538 = -t552 * mrSges(6,1) + t553 * mrSges(6,2);
t582 = qJD(4) + qJD(5);
t547 = -t582 * mrSges(6,2) + t552 * mrSges(6,3);
t580 = qJDD(4) + qJDD(5);
t512 = m(6) * t515 + t580 * mrSges(6,1) - t533 * mrSges(6,3) - t553 * t538 + t582 * t547;
t516 = t590 * t517 + t594 * t518;
t532 = -t553 * qJD(5) - t590 * t557 + t594 * t558;
t548 = t582 * mrSges(6,1) - t553 * mrSges(6,3);
t513 = m(6) * t516 - t580 * mrSges(6,2) + t532 * mrSges(6,3) + t552 * t538 - t582 * t548;
t503 = t594 * t512 + t590 * t513;
t550 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t591 + Ifges(5,2) * t595) * t583;
t551 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t591 + Ifges(5,4) * t595) * t583;
t535 = Ifges(6,4) * t553 + Ifges(6,2) * t552 + Ifges(6,6) * t582;
t536 = Ifges(6,1) * t553 + Ifges(6,4) * t552 + Ifges(6,5) * t582;
t602 = -mrSges(6,1) * t515 + mrSges(6,2) * t516 - Ifges(6,5) * t533 - Ifges(6,6) * t532 - Ifges(6,3) * t580 - t553 * t535 + t552 * t536;
t618 = mrSges(5,1) * t520 - mrSges(5,2) * t521 + Ifges(5,5) * t557 + Ifges(5,6) * t558 + Ifges(5,3) * qJDD(4) + pkin(4) * t503 + (t591 * t550 - t595 * t551) * t583 - t602;
t616 = t583 * t595;
t556 = (-mrSges(5,1) * t595 + mrSges(5,2) * t591) * t583;
t565 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t616;
t501 = m(5) * t520 + qJDD(4) * mrSges(5,1) - t557 * mrSges(5,3) + qJD(4) * t565 - t556 * t617 + t503;
t564 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t617;
t608 = -t590 * t512 + t594 * t513;
t502 = m(5) * t521 - qJDD(4) * mrSges(5,2) + t558 * mrSges(5,3) - qJD(4) * t564 + t556 * t616 + t608;
t609 = -t591 * t501 + t595 * t502;
t494 = m(4) * t527 - t579 * mrSges(4,1) - t581 * mrSges(4,2) + t609;
t526 = t596 * t542 - t592 * t543;
t605 = -t581 * pkin(3) - t526;
t523 = -t579 * pkin(7) + t605;
t519 = t566 * t617 - t558 * pkin(4) + (-pkin(8) * t586 - pkin(7)) * t579 + t605;
t603 = m(6) * t519 - t532 * mrSges(6,1) + t533 * mrSges(6,2) - t552 * t547 + t553 * t548;
t600 = -m(5) * t523 + t558 * mrSges(5,1) - t557 * mrSges(5,2) - t564 * t617 + t565 * t616 - t603;
t507 = m(4) * t526 + t581 * mrSges(4,1) - t579 * mrSges(4,2) + t600;
t489 = t592 * t494 + t596 * t507;
t486 = m(3) * t545 + qJDD(1) * mrSges(3,1) - t598 * mrSges(3,2) + t489;
t610 = t596 * t494 - t592 * t507;
t487 = m(3) * t546 - t598 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t610;
t479 = t589 * t486 + t588 * t487;
t497 = t595 * t501 + t591 * t502;
t614 = m(4) * t587 + t497;
t611 = -t588 * t486 + t589 * t487;
t476 = m(2) * t568 - t598 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t611;
t477 = m(2) * t569 + qJDD(1) * mrSges(2,1) - t598 * mrSges(2,2) + t479;
t612 = t597 * t476 - t593 * t477;
t495 = m(3) * t587 + t614;
t607 = -t593 * t476 - t597 * t477;
t534 = Ifges(6,5) * t553 + Ifges(6,6) * t552 + Ifges(6,3) * t582;
t504 = -mrSges(6,1) * t519 + mrSges(6,3) * t516 + Ifges(6,4) * t533 + Ifges(6,2) * t532 + Ifges(6,6) * t580 - t553 * t534 + t582 * t536;
t505 = mrSges(6,2) * t519 - mrSges(6,3) * t515 + Ifges(6,1) * t533 + Ifges(6,4) * t532 + Ifges(6,5) * t580 + t552 * t534 - t582 * t535;
t549 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t591 + Ifges(5,6) * t595) * t583;
t482 = -mrSges(5,1) * t523 + mrSges(5,3) * t521 + Ifges(5,4) * t557 + Ifges(5,2) * t558 + Ifges(5,6) * qJDD(4) - pkin(4) * t603 + pkin(8) * t608 + qJD(4) * t551 + t594 * t504 + t590 * t505 - t549 * t617;
t491 = mrSges(5,2) * t523 - mrSges(5,3) * t520 + Ifges(5,1) * t557 + Ifges(5,4) * t558 + Ifges(5,5) * qJDD(4) - pkin(8) * t503 - qJD(4) * t550 - t590 * t504 + t594 * t505 + t549 * t616;
t604 = mrSges(4,1) * t526 - mrSges(4,2) * t527 + Ifges(4,3) * t581 + pkin(3) * t600 + pkin(7) * t609 + t595 * t482 + t591 * t491;
t601 = mrSges(2,1) * t569 + mrSges(3,1) * t545 - mrSges(2,2) * t568 - mrSges(3,2) * t546 + pkin(1) * t479 + pkin(2) * t489 + t604 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t480 = -mrSges(4,1) * t587 + mrSges(4,3) * t527 + t579 * Ifges(4,5) + Ifges(4,6) * t581 - pkin(3) * t497 - t618;
t474 = mrSges(4,2) * t587 - mrSges(4,3) * t526 + Ifges(4,5) * t581 - t579 * Ifges(4,6) - pkin(7) * t497 - t591 * t482 + t595 * t491;
t473 = mrSges(3,2) * t587 - mrSges(3,3) * t545 + Ifges(3,5) * qJDD(1) - t598 * Ifges(3,6) - pkin(6) * t489 + t596 * t474 - t592 * t480;
t472 = -mrSges(3,1) * t587 + mrSges(3,3) * t546 + t598 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t614 + pkin(6) * t610 + t592 * t474 + t596 * t480;
t471 = -mrSges(2,2) * g(1) - mrSges(2,3) * t569 + Ifges(2,5) * qJDD(1) - t598 * Ifges(2,6) - qJ(2) * t479 - t588 * t472 + t589 * t473;
t470 = mrSges(2,1) * g(1) + mrSges(2,3) * t568 + t598 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t495 + qJ(2) * t611 + t589 * t472 + t588 * t473;
t1 = [(-m(1) - m(2)) * g(1) + t495; -m(1) * g(2) + t607; -m(1) * g(3) + t612; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t601; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t612 - t597 * t470 - t593 * t471; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t607 - t593 * t470 + t597 * t471; t601; t495; t604; t618; -t602;];
tauJB = t1;
