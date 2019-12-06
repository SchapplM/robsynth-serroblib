% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:50
% EndTime: 2019-12-05 16:41:52
% DurationCPUTime: 1.63s
% Computational Cost: add. (18285->193), mult. (24393->240), div. (0->0), fcn. (13624->8), ass. (0->88)
t579 = Ifges(5,1) + Ifges(6,1);
t572 = Ifges(5,4) - Ifges(6,5);
t571 = Ifges(5,5) + Ifges(6,4);
t578 = Ifges(5,2) + Ifges(6,3);
t570 = Ifges(5,6) - Ifges(6,6);
t577 = Ifges(5,3) + Ifges(6,2);
t533 = qJD(2) + qJD(3);
t540 = sin(qJ(4));
t543 = cos(qJ(4));
t510 = (-mrSges(6,1) * t543 - mrSges(6,3) * t540) * t533;
t532 = qJDD(2) + qJDD(3);
t561 = qJD(4) * t533;
t512 = t540 * t532 + t543 * t561;
t538 = sin(pkin(8));
t539 = cos(pkin(8));
t523 = t538 * g(1) - t539 * g(2);
t524 = -t539 * g(1) - t538 * g(2);
t542 = sin(qJ(2));
t545 = cos(qJ(2));
t495 = t545 * t523 - t542 * t524;
t492 = qJDD(2) * pkin(2) + t495;
t496 = t542 * t523 + t545 * t524;
t547 = qJD(2) ^ 2;
t493 = -t547 * pkin(2) + t496;
t541 = sin(qJ(3));
t544 = cos(qJ(3));
t488 = t541 * t492 + t544 * t493;
t531 = t533 ^ 2;
t485 = -t531 * pkin(3) + t532 * pkin(7) + t488;
t509 = (-pkin(4) * t543 - qJ(5) * t540) * t533;
t546 = qJD(4) ^ 2;
t537 = -g(3) + qJDD(1);
t566 = t543 * t537;
t480 = -qJDD(4) * pkin(4) - t546 * qJ(5) - t566 + qJDD(5) + (t509 * t533 + t485) * t540;
t567 = t533 * t543;
t522 = mrSges(6,2) * t567 + qJD(4) * mrSges(6,3);
t553 = -m(6) * t480 + qJDD(4) * mrSges(6,1) + qJD(4) * t522;
t568 = t533 * t540;
t476 = t512 * mrSges(6,2) + t510 * t568 - t553;
t482 = t543 * t485 + t540 * t537;
t479 = -t546 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t509 * t567 + t482;
t481 = -t540 * t485 + t566;
t513 = t543 * t532 - t540 * t561;
t520 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t568;
t555 = m(6) * t479 + qJDD(4) * mrSges(6,3) + qJD(4) * t520 + t510 * t567;
t562 = (t579 * t540 + t572 * t543) * t533 + t571 * qJD(4);
t564 = (-t572 * t540 - t578 * t543) * t533 - t570 * qJD(4);
t576 = -(t540 * t564 + t543 * t562) * t533 + t577 * qJDD(4) + t571 * t512 + t570 * t513 + mrSges(5,1) * t481 - mrSges(6,1) * t480 - mrSges(5,2) * t482 + mrSges(6,3) * t479 - pkin(4) * t476 + qJ(5) * (t513 * mrSges(6,2) + t555);
t573 = mrSges(5,3) + mrSges(6,2);
t511 = (-mrSges(5,1) * t543 + mrSges(5,2) * t540) * t533;
t519 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t568;
t472 = m(5) * t482 - qJDD(4) * mrSges(5,2) - qJD(4) * t519 + t511 * t567 + t573 * t513 + t555;
t521 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t567;
t473 = m(5) * t481 + qJDD(4) * mrSges(5,1) + qJD(4) * t521 + (-t510 - t511) * t568 - t573 * t512 + t553;
t556 = t543 * t472 - t540 * t473;
t463 = m(4) * t488 - t531 * mrSges(4,1) - t532 * mrSges(4,2) + t556;
t487 = t544 * t492 - t541 * t493;
t484 = -t532 * pkin(3) - t531 * pkin(7) - t487;
t477 = -t513 * pkin(4) - t512 * qJ(5) + (-0.2e1 * qJD(5) * t540 + (pkin(4) * t540 - qJ(5) * t543) * qJD(4)) * t533 + t484;
t474 = m(6) * t477 - t513 * mrSges(6,1) - t512 * mrSges(6,3) - t520 * t568 - t522 * t567;
t548 = -m(5) * t484 + t513 * mrSges(5,1) - t512 * mrSges(5,2) - t519 * t568 + t521 * t567 - t474;
t467 = m(4) * t487 + t532 * mrSges(4,1) - t531 * mrSges(4,2) + t548;
t456 = t541 * t463 + t544 * t467;
t453 = m(3) * t495 + qJDD(2) * mrSges(3,1) - t547 * mrSges(3,2) + t456;
t557 = t544 * t463 - t541 * t467;
t454 = m(3) * t496 - t547 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t557;
t447 = t545 * t453 + t542 * t454;
t445 = m(2) * t523 + t447;
t558 = -t542 * t453 + t545 * t454;
t446 = m(2) * t524 + t558;
t565 = t539 * t445 + t538 * t446;
t465 = t540 * t472 + t543 * t473;
t563 = (t571 * t540 + t570 * t543) * t533 + t577 * qJD(4);
t560 = m(4) * t537 + t465;
t559 = -t538 * t445 + t539 * t446;
t554 = m(3) * t537 + t560;
t552 = m(2) * t537 + t554;
t459 = -mrSges(5,1) * t484 - mrSges(6,1) * t477 + mrSges(6,2) * t479 + mrSges(5,3) * t482 - pkin(4) * t474 + t562 * qJD(4) + t570 * qJDD(4) + t572 * t512 + t578 * t513 - t563 * t568;
t460 = mrSges(5,2) * t484 + mrSges(6,2) * t480 - mrSges(5,3) * t481 - mrSges(6,3) * t477 - qJ(5) * t474 + t564 * qJD(4) + t571 * qJDD(4) + t579 * t512 + t572 * t513 + t563 * t567;
t551 = mrSges(4,1) * t487 - mrSges(4,2) * t488 + Ifges(4,3) * t532 + pkin(3) * t548 + pkin(7) * t556 + t543 * t459 + t540 * t460;
t550 = mrSges(3,1) * t495 - mrSges(3,2) * t496 + Ifges(3,3) * qJDD(2) + pkin(2) * t456 + t551;
t449 = -mrSges(4,1) * t537 + mrSges(4,3) * t488 + t531 * Ifges(4,5) + Ifges(4,6) * t532 - pkin(3) * t465 - t576;
t448 = mrSges(4,2) * t537 - mrSges(4,3) * t487 + Ifges(4,5) * t532 - t531 * Ifges(4,6) - pkin(7) * t465 - t540 * t459 + t543 * t460;
t441 = mrSges(3,2) * t537 - mrSges(3,3) * t495 + Ifges(3,5) * qJDD(2) - t547 * Ifges(3,6) - pkin(6) * t456 + t544 * t448 - t541 * t449;
t440 = -mrSges(3,1) * t537 + mrSges(3,3) * t496 + t547 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t560 + pkin(6) * t557 + t541 * t448 + t544 * t449;
t439 = mrSges(2,2) * t537 - mrSges(2,3) * t523 - pkin(5) * t447 - t542 * t440 + t545 * t441;
t438 = -mrSges(2,1) * t537 + mrSges(2,3) * t524 - pkin(1) * t554 + pkin(5) * t558 + t545 * t440 + t542 * t441;
t1 = [-m(1) * g(1) + t559; -m(1) * g(2) + t565; -m(1) * g(3) + t552; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t565 - t538 * t438 + t539 * t439; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t559 + t539 * t438 + t538 * t439; -mrSges(1,1) * g(2) + mrSges(2,1) * t523 + mrSges(1,2) * g(1) - mrSges(2,2) * t524 + pkin(1) * t447 + t550; t552; t550; t551; t576; t476;];
tauJB = t1;
