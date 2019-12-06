% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:18
% EndTime: 2019-12-05 18:20:20
% DurationCPUTime: 2.69s
% Computational Cost: add. (39285->208), mult. (51725->271), div. (0->0), fcn. (28744->10), ass. (0->103)
t604 = 2 * qJD(4);
t569 = sin(qJ(1));
t572 = cos(qJ(1));
t547 = t572 * g(2) + t569 * g(3);
t540 = qJDD(1) * pkin(1) + t547;
t546 = t569 * g(2) - t572 * g(3);
t573 = qJD(1) ^ 2;
t541 = -t573 * pkin(1) + t546;
t568 = sin(qJ(2));
t571 = cos(qJ(2));
t526 = t571 * t540 - t568 * t541;
t559 = qJDD(1) + qJDD(2);
t520 = t559 * pkin(2) + t526;
t527 = t568 * t540 + t571 * t541;
t560 = (qJD(1) + qJD(2));
t558 = t560 ^ 2;
t521 = -t558 * pkin(2) + t527;
t564 = sin(pkin(8));
t566 = cos(pkin(8));
t516 = t564 * t520 + t566 * t521;
t513 = -t558 * pkin(3) + t559 * qJ(4) + t516;
t603 = (t560 * t604) + t513;
t563 = sin(pkin(9));
t602 = mrSges(5,2) * t563;
t600 = mrSges(5,3) * t559;
t599 = t563 * t560;
t567 = sin(qJ(5));
t598 = t563 * t567;
t570 = cos(qJ(5));
t597 = t563 * t570;
t565 = cos(pkin(9));
t596 = t565 * t559;
t595 = t565 * t560;
t562 = -g(1) + qJDD(3);
t594 = t565 * t562;
t509 = t563 * t562 + t603 * t565;
t534 = (-mrSges(5,1) * t565 + t602) * t560;
t584 = -pkin(4) * t565 - pkin(7) * t563;
t536 = t584 * t560;
t507 = t536 * t595 + t509;
t515 = t566 * t520 - t564 * t521;
t578 = -t558 * qJ(4) + qJDD(4) - t515;
t510 = (-pkin(3) + t584) * t559 + t578;
t504 = -t567 * t507 + t570 * t510;
t544 = qJD(5) - t595;
t591 = t560 * t598;
t529 = -t544 * mrSges(6,2) - mrSges(6,3) * t591;
t531 = (mrSges(6,1) * t567 + mrSges(6,2) * t570) * t599;
t592 = qJD(5) * t560;
t533 = (t559 * t570 - t567 * t592) * t563;
t543 = qJDD(5) - t596;
t590 = t560 * t597;
t502 = m(6) * t504 + t543 * mrSges(6,1) - t533 * mrSges(6,3) + t544 * t529 - t531 * t590;
t505 = t570 * t507 + t567 * t510;
t530 = t544 * mrSges(6,1) - mrSges(6,3) * t590;
t532 = (-t559 * t567 - t570 * t592) * t563;
t503 = m(6) * t505 - t543 * mrSges(6,2) + t532 * mrSges(6,3) - t544 * t530 - t531 * t591;
t585 = -t567 * t502 + t570 * t503;
t493 = m(5) * t509 + (t534 * t560 + t600) * t565 + t585;
t508 = -t603 * t563 + t594;
t506 = -t594 + (t513 + (t604 + t536) * t560) * t563;
t579 = -m(6) * t506 + t532 * mrSges(6,1) - t533 * mrSges(6,2);
t500 = m(5) * t508 + (-t600 + (-t529 * t567 - t530 * t570 - t534) * t560) * t563 + t579;
t586 = t565 * t493 - t563 * t500;
t485 = m(4) * t516 - t558 * mrSges(4,1) - t559 * mrSges(4,2) + t586;
t496 = t570 * t502 + t567 * t503;
t512 = -t559 * pkin(3) + t578;
t577 = -m(5) * t512 + mrSges(5,1) * t596 - t496 + (t563 ^ 2 + t565 ^ 2) * mrSges(5,3) * t558;
t490 = m(4) * t515 - t558 * mrSges(4,2) + (mrSges(4,1) - t602) * t559 + t577;
t478 = t564 * t485 + t566 * t490;
t475 = m(3) * t526 + t559 * mrSges(3,1) - t558 * mrSges(3,2) + t478;
t587 = t566 * t485 - t564 * t490;
t476 = m(3) * t527 - t558 * mrSges(3,1) - t559 * mrSges(3,2) + t587;
t470 = t571 * t475 + t568 * t476;
t488 = t563 * t493 + t565 * t500;
t486 = m(4) * t562 + t488;
t588 = -t568 * t475 + t571 * t476;
t467 = m(2) * t546 - t573 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t588;
t468 = m(2) * t547 + qJDD(1) * mrSges(2,1) - t573 * mrSges(2,2) + t470;
t589 = t572 * t467 - t569 * t468;
t583 = Ifges(5,1) * t563 + Ifges(5,4) * t565;
t582 = Ifges(5,5) * t563 + Ifges(5,6) * t565;
t581 = -t569 * t467 - t572 * t468;
t523 = Ifges(6,6) * t544 + (Ifges(6,4) * t570 - Ifges(6,2) * t567) * t599;
t524 = Ifges(6,5) * t544 + (Ifges(6,1) * t570 - Ifges(6,4) * t567) * t599;
t580 = t523 * t570 + t524 * t567;
t576 = mrSges(6,1) * t504 - mrSges(6,2) * t505 + Ifges(6,5) * t533 + Ifges(6,6) * t532 + Ifges(6,3) * t543;
t522 = Ifges(6,3) * t544 + (Ifges(6,5) * t570 - Ifges(6,6) * t567) * t599;
t497 = -mrSges(6,1) * t506 + mrSges(6,3) * t505 + Ifges(6,4) * t533 + Ifges(6,2) * t532 + Ifges(6,6) * t543 - t522 * t590 + t544 * t524;
t498 = mrSges(6,2) * t506 - mrSges(6,3) * t504 + Ifges(6,1) * t533 + Ifges(6,4) * t532 + Ifges(6,5) * t543 - t522 * t591 - t544 * t523;
t535 = t582 * t560;
t480 = mrSges(5,2) * t512 - mrSges(5,3) * t508 - pkin(7) * t496 - t567 * t497 + t570 * t498 + t535 * t595 + t583 * t559;
t482 = Ifges(5,2) * t596 - mrSges(5,1) * t512 + mrSges(5,3) * t509 - pkin(4) * t496 + (Ifges(5,4) * t559 + (-t535 - t580) * t560) * t563 - t576;
t495 = t559 * t602 - t577;
t575 = mrSges(3,1) * t526 + mrSges(4,1) * t515 - mrSges(3,2) * t527 - mrSges(4,2) * t516 + pkin(2) * t478 - pkin(3) * t495 + qJ(4) * t586 + t563 * t480 + t565 * t482 + (Ifges(3,3) + Ifges(4,3)) * t559;
t574 = mrSges(2,1) * t547 - mrSges(2,2) * t546 + Ifges(2,3) * qJDD(1) + pkin(1) * t470 + t575;
t471 = t558 * Ifges(4,5) - mrSges(4,1) * t562 + mrSges(4,3) * t516 - mrSges(5,1) * t508 + mrSges(5,2) * t509 - t567 * t498 - t570 * t497 - pkin(4) * t579 - pkin(7) * t585 - pkin(3) * t488 + (Ifges(4,6) - t582) * t559 + (-pkin(4) * (-t529 * t598 - t530 * t597) + (-t563 * (Ifges(5,4) * t563 + Ifges(5,2) * t565) + t565 * t583) * t560) * t560;
t465 = mrSges(4,2) * t562 - mrSges(4,3) * t515 + Ifges(4,5) * t559 - t558 * Ifges(4,6) - qJ(4) * t488 + t565 * t480 - t563 * t482;
t464 = -mrSges(3,2) * g(1) - mrSges(3,3) * t526 + Ifges(3,5) * t559 - t558 * Ifges(3,6) - qJ(3) * t478 + t566 * t465 - t564 * t471;
t463 = mrSges(3,1) * g(1) + mrSges(3,3) * t527 + t558 * Ifges(3,5) + Ifges(3,6) * t559 - pkin(2) * t486 + qJ(3) * t587 + t564 * t465 + t566 * t471;
t462 = -mrSges(2,2) * g(1) - mrSges(2,3) * t547 + Ifges(2,5) * qJDD(1) - t573 * Ifges(2,6) - pkin(6) * t470 - t568 * t463 + t571 * t464;
t461 = Ifges(2,6) * qJDD(1) + t573 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t546 + t568 * t464 + t571 * t463 - pkin(1) * (-m(3) * g(1) + t486) + pkin(6) * t588;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t486; -m(1) * g(2) + t581; -m(1) * g(3) + t589; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t574; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t589 - t572 * t461 - t569 * t462; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t581 - t569 * t461 + t572 * t462; t574; t575; t486; t495; t580 * t599 + t576;];
tauJB = t1;
