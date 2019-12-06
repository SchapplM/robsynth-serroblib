% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR6
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:51
% EndTime: 2019-12-05 17:09:54
% DurationCPUTime: 2.67s
% Computational Cost: add. (38644->220), mult. (49072->281), div. (0->0), fcn. (30201->10), ass. (0->98)
t561 = sin(pkin(9));
t590 = cos(pkin(9));
t547 = -t590 * g(1) - t561 * g(2);
t560 = -g(3) + qJDD(1);
t565 = sin(qJ(2));
t569 = cos(qJ(2));
t528 = -t565 * t547 + t569 * t560;
t526 = qJDD(2) * pkin(2) + t528;
t529 = t569 * t547 + t565 * t560;
t570 = qJD(2) ^ 2;
t527 = -t570 * pkin(2) + t529;
t564 = sin(qJ(3));
t568 = cos(qJ(3));
t515 = t564 * t526 + t568 * t527;
t558 = qJD(2) + qJD(3);
t554 = t558 ^ 2;
t556 = qJDD(2) + qJDD(3);
t507 = -t554 * pkin(3) + t556 * pkin(7) + t515;
t546 = t561 * g(1) - t590 * g(2);
t563 = sin(qJ(4));
t567 = cos(qJ(4));
t502 = -t563 * t507 - t567 * t546;
t586 = qJD(4) * t558;
t584 = t567 * t586;
t538 = t563 * t556 + t584;
t499 = (-t538 + t584) * pkin(8) + (t554 * t563 * t567 + qJDD(4)) * pkin(4) + t502;
t503 = t567 * t507 - t563 * t546;
t539 = t567 * t556 - t563 * t586;
t589 = t558 * t563;
t545 = qJD(4) * pkin(4) - pkin(8) * t589;
t559 = t567 ^ 2;
t500 = -t559 * t554 * pkin(4) + t539 * pkin(8) - qJD(4) * t545 + t503;
t562 = sin(qJ(5));
t566 = cos(qJ(5));
t497 = t566 * t499 - t562 * t500;
t533 = (-t562 * t563 + t566 * t567) * t558;
t512 = t533 * qJD(5) + t566 * t538 + t562 * t539;
t534 = (t562 * t567 + t563 * t566) * t558;
t520 = -t533 * mrSges(6,1) + t534 * mrSges(6,2);
t557 = qJD(4) + qJD(5);
t521 = -t557 * mrSges(6,2) + t533 * mrSges(6,3);
t555 = qJDD(4) + qJDD(5);
t494 = m(6) * t497 + t555 * mrSges(6,1) - t512 * mrSges(6,3) - t534 * t520 + t557 * t521;
t498 = t562 * t499 + t566 * t500;
t511 = -t534 * qJD(5) - t562 * t538 + t566 * t539;
t522 = t557 * mrSges(6,1) - t534 * mrSges(6,3);
t495 = m(6) * t498 - t555 * mrSges(6,2) + t511 * mrSges(6,3) + t533 * t520 - t557 * t522;
t484 = t566 * t494 + t562 * t495;
t531 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t563 + Ifges(5,2) * t567) * t558;
t532 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t563 + Ifges(5,4) * t567) * t558;
t517 = Ifges(6,4) * t534 + Ifges(6,2) * t533 + Ifges(6,6) * t557;
t518 = Ifges(6,1) * t534 + Ifges(6,4) * t533 + Ifges(6,5) * t557;
t574 = -mrSges(6,1) * t497 + mrSges(6,2) * t498 - Ifges(6,5) * t512 - Ifges(6,6) * t511 - Ifges(6,3) * t555 - t534 * t517 + t533 * t518;
t592 = mrSges(5,1) * t502 - mrSges(5,2) * t503 + Ifges(5,5) * t538 + Ifges(5,6) * t539 + Ifges(5,3) * qJDD(4) + pkin(4) * t484 + (t563 * t531 - t567 * t532) * t558 - t574;
t591 = m(3) + m(4);
t588 = t558 * t567;
t537 = (-mrSges(5,1) * t567 + mrSges(5,2) * t563) * t558;
t544 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t588;
t482 = m(5) * t502 + qJDD(4) * mrSges(5,1) - t538 * mrSges(5,3) + qJD(4) * t544 - t537 * t589 + t484;
t543 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t589;
t579 = -t562 * t494 + t566 * t495;
t483 = m(5) * t503 - qJDD(4) * mrSges(5,2) + t539 * mrSges(5,3) - qJD(4) * t543 + t537 * t588 + t579;
t580 = -t563 * t482 + t567 * t483;
t473 = m(4) * t515 - t554 * mrSges(4,1) - t556 * mrSges(4,2) + t580;
t514 = t568 * t526 - t564 * t527;
t577 = -t556 * pkin(3) - t514;
t506 = -t554 * pkin(7) + t577;
t501 = t545 * t589 - t539 * pkin(4) + (-pkin(8) * t559 - pkin(7)) * t554 + t577;
t575 = m(6) * t501 - t511 * mrSges(6,1) + t512 * mrSges(6,2) - t533 * t521 + t534 * t522;
t572 = -m(5) * t506 + t539 * mrSges(5,1) - t538 * mrSges(5,2) - t543 * t589 + t544 * t588 - t575;
t488 = m(4) * t514 + t556 * mrSges(4,1) - t554 * mrSges(4,2) + t572;
t468 = t564 * t473 + t568 * t488;
t466 = m(3) * t528 + qJDD(2) * mrSges(3,1) - t570 * mrSges(3,2) + t468;
t581 = t568 * t473 - t564 * t488;
t467 = m(3) * t529 - t570 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t581;
t582 = -t565 * t466 + t569 * t467;
t458 = m(2) * t547 + t582;
t477 = t567 * t482 + t563 * t483;
t475 = (m(2) + t591) * t546 - t477;
t587 = t561 * t458 + t590 * t475;
t459 = t569 * t466 + t565 * t467;
t585 = m(2) * t560 + t459;
t583 = t590 * t458 - t561 * t475;
t516 = Ifges(6,5) * t534 + Ifges(6,6) * t533 + Ifges(6,3) * t557;
t485 = -mrSges(6,1) * t501 + mrSges(6,3) * t498 + Ifges(6,4) * t512 + Ifges(6,2) * t511 + Ifges(6,6) * t555 - t534 * t516 + t557 * t518;
t486 = mrSges(6,2) * t501 - mrSges(6,3) * t497 + Ifges(6,1) * t512 + Ifges(6,4) * t511 + Ifges(6,5) * t555 + t533 * t516 - t557 * t517;
t530 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t563 + Ifges(5,6) * t567) * t558;
t462 = -mrSges(5,1) * t506 + mrSges(5,3) * t503 + Ifges(5,4) * t538 + Ifges(5,2) * t539 + Ifges(5,6) * qJDD(4) - pkin(4) * t575 + pkin(8) * t579 + qJD(4) * t532 + t566 * t485 + t562 * t486 - t530 * t589;
t470 = mrSges(5,2) * t506 - mrSges(5,3) * t502 + Ifges(5,1) * t538 + Ifges(5,4) * t539 + Ifges(5,5) * qJDD(4) - pkin(8) * t484 - qJD(4) * t531 - t562 * t485 + t566 * t486 + t530 * t588;
t576 = mrSges(4,1) * t514 - mrSges(4,2) * t515 + Ifges(4,3) * t556 + pkin(3) * t572 + pkin(7) * t580 + t567 * t462 + t563 * t470;
t573 = mrSges(3,1) * t528 - mrSges(3,2) * t529 + Ifges(3,3) * qJDD(2) + pkin(2) * t468 + t576;
t460 = mrSges(4,1) * t546 + mrSges(4,3) * t515 + t554 * Ifges(4,5) + Ifges(4,6) * t556 - pkin(3) * t477 - t592;
t455 = -mrSges(4,2) * t546 - mrSges(4,3) * t514 + Ifges(4,5) * t556 - t554 * Ifges(4,6) - pkin(7) * t477 - t563 * t462 + t567 * t470;
t454 = -mrSges(3,2) * t546 - mrSges(3,3) * t528 + Ifges(3,5) * qJDD(2) - t570 * Ifges(3,6) - pkin(6) * t468 + t568 * t455 - t564 * t460;
t453 = -mrSges(2,1) * t560 + mrSges(2,3) * t547 - pkin(1) * t459 - t573;
t452 = Ifges(3,6) * qJDD(2) + t570 * Ifges(3,5) + mrSges(3,1) * t546 + mrSges(3,3) * t529 + t564 * t455 + t568 * t460 - pkin(2) * (-m(4) * t546 + t477) + pkin(6) * t581;
t451 = mrSges(2,2) * t560 - mrSges(2,3) * t546 - pkin(5) * t459 - t565 * t452 + t569 * t454;
t1 = [-m(1) * g(1) + t583; -m(1) * g(2) + t587; -m(1) * g(3) + t585; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t587 + t590 * t451 - t561 * t453; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t583 + t561 * t451 + t590 * t453; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t547 + t565 * t454 + t569 * t452 - pkin(1) * t477 + pkin(5) * t582 + (pkin(1) * t591 + mrSges(2,1)) * t546; t585; t573; t576; t592; -t574;];
tauJB = t1;
