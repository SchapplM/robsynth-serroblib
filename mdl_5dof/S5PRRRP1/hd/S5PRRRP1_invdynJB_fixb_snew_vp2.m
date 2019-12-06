% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:01
% EndTime: 2019-12-05 16:40:03
% DurationCPUTime: 1.70s
% Computational Cost: add. (18661->196), mult. (25105->238), div. (0->0), fcn. (13985->8), ass. (0->93)
t603 = Ifges(5,1) + Ifges(6,1);
t595 = Ifges(5,4) + Ifges(6,4);
t594 = Ifges(5,5) + Ifges(6,5);
t602 = Ifges(5,2) + Ifges(6,2);
t593 = Ifges(5,6) + Ifges(6,6);
t601 = Ifges(5,3) + Ifges(6,3);
t553 = qJD(2) + qJD(3);
t560 = sin(qJ(4));
t563 = cos(qJ(4));
t530 = (-mrSges(6,1) * t563 + mrSges(6,2) * t560) * t553;
t552 = qJDD(2) + qJDD(3);
t584 = qJD(4) * t553;
t579 = t563 * t584;
t532 = t552 * t560 + t579;
t558 = sin(pkin(8));
t559 = cos(pkin(8));
t543 = g(1) * t558 - g(2) * t559;
t544 = -g(1) * t559 - g(2) * t558;
t562 = sin(qJ(2));
t565 = cos(qJ(2));
t516 = t565 * t543 - t544 * t562;
t513 = qJDD(2) * pkin(2) + t516;
t517 = t562 * t543 + t565 * t544;
t566 = qJD(2) ^ 2;
t514 = -pkin(2) * t566 + t517;
t561 = sin(qJ(3));
t564 = cos(qJ(3));
t509 = t561 * t513 + t564 * t514;
t551 = t553 ^ 2;
t506 = -pkin(3) * t551 + pkin(7) * t552 + t509;
t557 = -g(3) + qJDD(1);
t546 = t563 * t557;
t583 = qJD(5) * t553;
t597 = pkin(4) * t551;
t499 = qJDD(4) * pkin(4) + t546 + (-t532 + t579) * qJ(5) + (t563 * t597 - t506 - 0.2e1 * t583) * t560;
t590 = t553 * t563;
t541 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t590;
t582 = m(6) * t499 + qJDD(4) * mrSges(6,1) + qJD(4) * t541;
t591 = t553 * t560;
t496 = -t532 * mrSges(6,3) - t530 * t591 + t582;
t503 = t563 * t506 + t560 * t557;
t533 = t552 * t563 - t560 * t584;
t538 = qJD(4) * pkin(4) - qJ(5) * t591;
t556 = t563 ^ 2;
t500 = qJ(5) * t533 - qJD(4) * t538 - t556 * t597 + 0.2e1 * t563 * t583 + t503;
t502 = -t560 * t506 + t546;
t586 = (t560 * t603 + t563 * t595) * t553 + t594 * qJD(4);
t587 = (t560 * t595 + t563 * t602) * t553 + t593 * qJD(4);
t600 = mrSges(5,1) * t502 + mrSges(6,1) * t499 - mrSges(5,2) * t503 - mrSges(6,2) * t500 + pkin(4) * t496 + t601 * qJDD(4) + t594 * t532 + t593 * t533 + (t560 * t587 - t563 * t586) * t553;
t596 = -mrSges(5,2) - mrSges(6,2);
t531 = (-mrSges(5,1) * t563 + mrSges(5,2) * t560) * t553;
t542 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t590;
t493 = m(5) * t502 + qJDD(4) * mrSges(5,1) + qJD(4) * t542 + (-t530 - t531) * t591 + (-mrSges(5,3) - mrSges(6,3)) * t532 + t582;
t581 = m(6) * t500 + t533 * mrSges(6,3) + t530 * t590;
t539 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t591;
t585 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t591 - t539;
t494 = m(5) * t503 + t533 * mrSges(5,3) + qJD(4) * t585 + qJDD(4) * t596 + t531 * t590 + t581;
t575 = -t493 * t560 + t494 * t563;
t484 = m(4) * t509 - mrSges(4,1) * t551 - mrSges(4,2) * t552 + t575;
t508 = t564 * t513 - t561 * t514;
t571 = -t552 * pkin(3) - t508;
t505 = -pkin(7) * t551 + t571;
t501 = t538 * t591 - t533 * pkin(4) + qJDD(5) + (-qJ(5) * t556 - pkin(7)) * t551 + t571;
t573 = -m(6) * t501 + t533 * mrSges(6,1) + t541 * t590;
t567 = -m(5) * t505 + t533 * mrSges(5,1) + t532 * t596 + t542 * t590 + t585 * t591 + t573;
t488 = m(4) * t508 + t552 * mrSges(4,1) - t551 * mrSges(4,2) + t567;
t477 = t484 * t561 + t488 * t564;
t474 = m(3) * t516 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t566 + t477;
t576 = t484 * t564 - t488 * t561;
t475 = m(3) * t517 - mrSges(3,1) * t566 - qJDD(2) * mrSges(3,2) + t576;
t469 = t474 * t565 + t475 * t562;
t467 = m(2) * t543 + t469;
t577 = -t474 * t562 + t475 * t565;
t468 = m(2) * t544 + t577;
t589 = t467 * t559 + t468 * t558;
t486 = t493 * t563 + t494 * t560;
t588 = (t560 * t594 + t563 * t593) * t553 + t601 * qJD(4);
t580 = m(4) * t557 + t486;
t578 = -t467 * t558 + t468 * t559;
t574 = m(3) * t557 + t580;
t572 = m(2) * t557 + t574;
t495 = t532 * mrSges(6,2) + t539 * t591 - t573;
t479 = -mrSges(5,1) * t505 + mrSges(5,3) * t503 - mrSges(6,1) * t501 + mrSges(6,3) * t500 - pkin(4) * t495 + qJ(5) * t581 - t588 * t591 + t602 * t533 + t595 * t532 + (-mrSges(6,2) * qJ(5) + t593) * qJDD(4) + (-qJ(5) * t539 + t586) * qJD(4);
t481 = mrSges(5,2) * t505 + mrSges(6,2) * t501 - mrSges(5,3) * t502 - mrSges(6,3) * t499 - qJ(5) * t496 - t587 * qJD(4) + t594 * qJDD(4) + t532 * t603 + t595 * t533 + t588 * t590;
t570 = mrSges(4,1) * t508 - mrSges(4,2) * t509 + Ifges(4,3) * t552 + pkin(3) * t567 + pkin(7) * t575 + t479 * t563 + t481 * t560;
t568 = mrSges(3,1) * t516 - mrSges(3,2) * t517 + Ifges(3,3) * qJDD(2) + pkin(2) * t477 + t570;
t470 = -mrSges(4,1) * t557 + mrSges(4,3) * t509 + t551 * Ifges(4,5) + Ifges(4,6) * t552 - pkin(3) * t486 - t600;
t463 = mrSges(4,2) * t557 - mrSges(4,3) * t508 + Ifges(4,5) * t552 - Ifges(4,6) * t551 - pkin(7) * t486 - t479 * t560 + t481 * t563;
t462 = mrSges(3,2) * t557 - mrSges(3,3) * t516 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t566 - pkin(6) * t477 + t463 * t564 - t470 * t561;
t461 = -mrSges(3,1) * t557 + mrSges(3,3) * t517 + t566 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t580 + pkin(6) * t576 + t561 * t463 + t564 * t470;
t460 = mrSges(2,2) * t557 - mrSges(2,3) * t543 - pkin(5) * t469 - t461 * t562 + t462 * t565;
t459 = -mrSges(2,1) * t557 + mrSges(2,3) * t544 - pkin(1) * t574 + pkin(5) * t577 + t565 * t461 + t562 * t462;
t1 = [-m(1) * g(1) + t578; -m(1) * g(2) + t589; -m(1) * g(3) + t572; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t589 - t459 * t558 + t460 * t559; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t578 + t559 * t459 + t558 * t460; -mrSges(1,1) * g(2) + mrSges(2,1) * t543 + mrSges(1,2) * g(1) - mrSges(2,2) * t544 + pkin(1) * t469 + t568; t572; t568; t570; t600; t495;];
tauJB = t1;
