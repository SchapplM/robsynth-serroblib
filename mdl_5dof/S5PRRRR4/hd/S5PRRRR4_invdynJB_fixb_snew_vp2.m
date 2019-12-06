% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:54
% EndTime: 2019-12-05 17:07:57
% DurationCPUTime: 2.89s
% Computational Cost: add. (44055->218), mult. (58801->281), div. (0->0), fcn. (37282->10), ass. (0->99)
t559 = sin(pkin(9));
t560 = cos(pkin(9));
t542 = t559 * g(1) - t560 * g(2);
t543 = -t560 * g(1) - t559 * g(2);
t564 = sin(qJ(2));
t568 = cos(qJ(2));
t524 = t568 * t542 - t564 * t543;
t519 = qJDD(2) * pkin(2) + t524;
t525 = t564 * t542 + t568 * t543;
t569 = qJD(2) ^ 2;
t520 = -t569 * pkin(2) + t525;
t563 = sin(qJ(3));
t567 = cos(qJ(3));
t505 = t563 * t519 + t567 * t520;
t555 = qJD(2) + qJD(3);
t551 = t555 ^ 2;
t553 = qJDD(2) + qJDD(3);
t501 = -t551 * pkin(3) + t553 * pkin(7) + t505;
t558 = -g(3) + qJDD(1);
t562 = sin(qJ(4));
t566 = cos(qJ(4));
t497 = -t562 * t501 + t566 * t558;
t587 = qJD(4) * t555;
t585 = t566 * t587;
t534 = t562 * t553 + t585;
t494 = (-t534 + t585) * pkin(8) + (t551 * t562 * t566 + qJDD(4)) * pkin(4) + t497;
t498 = t566 * t501 + t562 * t558;
t535 = t566 * t553 - t562 * t587;
t590 = t555 * t562;
t541 = qJD(4) * pkin(4) - pkin(8) * t590;
t557 = t566 ^ 2;
t495 = -t557 * t551 * pkin(4) + t535 * pkin(8) - qJD(4) * t541 + t498;
t561 = sin(qJ(5));
t565 = cos(qJ(5));
t492 = t565 * t494 - t561 * t495;
t529 = (-t561 * t562 + t565 * t566) * t555;
t510 = t529 * qJD(5) + t565 * t534 + t561 * t535;
t530 = (t561 * t566 + t562 * t565) * t555;
t515 = -t529 * mrSges(6,1) + t530 * mrSges(6,2);
t554 = qJD(4) + qJD(5);
t522 = -t554 * mrSges(6,2) + t529 * mrSges(6,3);
t552 = qJDD(4) + qJDD(5);
t489 = m(6) * t492 + t552 * mrSges(6,1) - t510 * mrSges(6,3) - t530 * t515 + t554 * t522;
t493 = t561 * t494 + t565 * t495;
t509 = -t530 * qJD(5) - t561 * t534 + t565 * t535;
t523 = t554 * mrSges(6,1) - t530 * mrSges(6,3);
t490 = m(6) * t493 - t552 * mrSges(6,2) + t509 * mrSges(6,3) + t529 * t515 - t554 * t523;
t480 = t565 * t489 + t561 * t490;
t527 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t562 + Ifges(5,2) * t566) * t555;
t528 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t562 + Ifges(5,4) * t566) * t555;
t512 = Ifges(6,4) * t530 + Ifges(6,2) * t529 + Ifges(6,6) * t554;
t513 = Ifges(6,1) * t530 + Ifges(6,4) * t529 + Ifges(6,5) * t554;
t573 = -mrSges(6,1) * t492 + mrSges(6,2) * t493 - Ifges(6,5) * t510 - Ifges(6,6) * t509 - Ifges(6,3) * t552 - t530 * t512 + t529 * t513;
t591 = mrSges(5,1) * t497 - mrSges(5,2) * t498 + Ifges(5,5) * t534 + Ifges(5,6) * t535 + Ifges(5,3) * qJDD(4) + pkin(4) * t480 + (t562 * t527 - t566 * t528) * t555 - t573;
t589 = t555 * t566;
t533 = (-mrSges(5,1) * t566 + mrSges(5,2) * t562) * t555;
t540 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t589;
t478 = m(5) * t497 + qJDD(4) * mrSges(5,1) - t534 * mrSges(5,3) + qJD(4) * t540 - t533 * t590 + t480;
t539 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t590;
t580 = -t561 * t489 + t565 * t490;
t479 = m(5) * t498 - qJDD(4) * mrSges(5,2) + t535 * mrSges(5,3) - qJD(4) * t539 + t533 * t589 + t580;
t581 = -t562 * t478 + t566 * t479;
t472 = m(4) * t505 - t551 * mrSges(4,1) - t553 * mrSges(4,2) + t581;
t504 = t567 * t519 - t563 * t520;
t576 = -t553 * pkin(3) - t504;
t500 = -t551 * pkin(7) + t576;
t496 = t541 * t590 - t535 * pkin(4) + (-pkin(8) * t557 - pkin(7)) * t551 + t576;
t574 = m(6) * t496 - t509 * mrSges(6,1) + t510 * mrSges(6,2) - t529 * t522 + t530 * t523;
t571 = -m(5) * t500 + t535 * mrSges(5,1) - t534 * mrSges(5,2) - t539 * t590 + t540 * t589 - t574;
t484 = m(4) * t504 + t553 * mrSges(4,1) - t551 * mrSges(4,2) + t571;
t467 = t563 * t472 + t567 * t484;
t464 = m(3) * t524 + qJDD(2) * mrSges(3,1) - t569 * mrSges(3,2) + t467;
t582 = t567 * t472 - t563 * t484;
t465 = m(3) * t525 - t569 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t582;
t457 = t568 * t464 + t564 * t465;
t455 = m(2) * t542 + t457;
t583 = -t564 * t464 + t568 * t465;
t456 = m(2) * t543 + t583;
t588 = t560 * t455 + t559 * t456;
t474 = t566 * t478 + t562 * t479;
t586 = m(4) * t558 + t474;
t584 = -t559 * t455 + t560 * t456;
t579 = m(3) * t558 + t586;
t578 = m(2) * t558 + t579;
t511 = Ifges(6,5) * t530 + Ifges(6,6) * t529 + Ifges(6,3) * t554;
t481 = -mrSges(6,1) * t496 + mrSges(6,3) * t493 + Ifges(6,4) * t510 + Ifges(6,2) * t509 + Ifges(6,6) * t552 - t530 * t511 + t554 * t513;
t482 = mrSges(6,2) * t496 - mrSges(6,3) * t492 + Ifges(6,1) * t510 + Ifges(6,4) * t509 + Ifges(6,5) * t552 + t529 * t511 - t554 * t512;
t526 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t562 + Ifges(5,6) * t566) * t555;
t460 = -mrSges(5,1) * t500 + mrSges(5,3) * t498 + Ifges(5,4) * t534 + Ifges(5,2) * t535 + Ifges(5,6) * qJDD(4) - pkin(4) * t574 + pkin(8) * t580 + qJD(4) * t528 + t565 * t481 + t561 * t482 - t526 * t590;
t469 = mrSges(5,2) * t500 - mrSges(5,3) * t497 + Ifges(5,1) * t534 + Ifges(5,4) * t535 + Ifges(5,5) * qJDD(4) - pkin(8) * t480 - qJD(4) * t527 - t561 * t481 + t565 * t482 + t526 * t589;
t575 = mrSges(4,1) * t504 - mrSges(4,2) * t505 + Ifges(4,3) * t553 + pkin(3) * t571 + pkin(7) * t581 + t566 * t460 + t562 * t469;
t572 = mrSges(3,1) * t524 - mrSges(3,2) * t525 + Ifges(3,3) * qJDD(2) + pkin(2) * t467 + t575;
t458 = -mrSges(4,1) * t558 + mrSges(4,3) * t505 + t551 * Ifges(4,5) + Ifges(4,6) * t553 - pkin(3) * t474 - t591;
t451 = mrSges(4,2) * t558 - mrSges(4,3) * t504 + Ifges(4,5) * t553 - t551 * Ifges(4,6) - pkin(7) * t474 - t562 * t460 + t566 * t469;
t450 = mrSges(3,2) * t558 - mrSges(3,3) * t524 + Ifges(3,5) * qJDD(2) - t569 * Ifges(3,6) - pkin(6) * t467 + t567 * t451 - t563 * t458;
t449 = -mrSges(3,1) * t558 + mrSges(3,3) * t525 + t569 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t586 + pkin(6) * t582 + t563 * t451 + t567 * t458;
t448 = mrSges(2,2) * t558 - mrSges(2,3) * t542 - pkin(5) * t457 - t564 * t449 + t568 * t450;
t447 = -mrSges(2,1) * t558 + mrSges(2,3) * t543 - pkin(1) * t579 + pkin(5) * t583 + t568 * t449 + t564 * t450;
t1 = [-m(1) * g(1) + t584; -m(1) * g(2) + t588; -m(1) * g(3) + t578; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t588 - t559 * t447 + t560 * t448; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t584 + t560 * t447 + t559 * t448; -mrSges(1,1) * g(2) + mrSges(2,1) * t542 + mrSges(1,2) * g(1) - mrSges(2,2) * t543 + pkin(1) * t457 + t572; t578; t572; t575; t591; -t573;];
tauJB = t1;
