% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:28
% EndTime: 2019-12-31 17:46:29
% DurationCPUTime: 1.34s
% Computational Cost: add. (15868->201), mult. (30041->243), div. (0->0), fcn. (14346->8), ass. (0->93)
t536 = qJD(1) ^ 2;
t528 = sin(pkin(8));
t566 = t528 ^ 2;
t565 = -pkin(1) - pkin(2);
t564 = mrSges(2,1) + mrSges(3,1);
t563 = Ifges(3,4) + Ifges(2,5);
t562 = Ifges(2,6) - Ifges(3,6);
t530 = cos(pkin(8));
t561 = mrSges(5,1) * t530;
t560 = mrSges(5,2) * t528;
t518 = t530 ^ 2;
t559 = t518 * t536;
t533 = sin(qJ(1));
t535 = cos(qJ(1));
t508 = -t535 * g(1) - t533 * g(2);
t542 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t508;
t500 = -t536 * pkin(1) + t542;
t494 = t565 * t536 + t542;
t507 = t533 * g(1) - t535 * g(2);
t541 = -t536 * qJ(2) + qJDD(2) - t507;
t499 = t565 * qJDD(1) + t541;
t529 = sin(pkin(7));
t531 = cos(pkin(7));
t482 = t531 * t494 + t529 * t499;
t479 = -t536 * pkin(3) - qJDD(1) * qJ(4) + t482;
t525 = g(3) + qJDD(3);
t555 = qJD(1) * qJD(4);
t557 = t530 * t525 + 0.2e1 * t528 * t555;
t472 = (pkin(4) * t530 * t536 + pkin(6) * qJDD(1) - t479) * t528 + t557;
t476 = t528 * t525 + (t479 - (2 * t555)) * t530;
t554 = qJDD(1) * t530;
t473 = -pkin(4) * t559 - pkin(6) * t554 + t476;
t532 = sin(qJ(5));
t534 = cos(qJ(5));
t470 = t534 * t472 - t532 * t473;
t546 = t528 * t532 - t530 * t534;
t502 = t546 * qJD(1);
t547 = -t528 * t534 - t530 * t532;
t503 = t547 * qJD(1);
t487 = -t502 * mrSges(6,1) + t503 * mrSges(6,2);
t490 = t502 * qJD(5) + t547 * qJDD(1);
t495 = -qJD(5) * mrSges(6,2) + t502 * mrSges(6,3);
t468 = m(6) * t470 + qJDD(5) * mrSges(6,1) - t490 * mrSges(6,3) + qJD(5) * t495 - t503 * t487;
t471 = t532 * t472 + t534 * t473;
t489 = -t503 * qJD(5) + t546 * qJDD(1);
t496 = qJD(5) * mrSges(6,1) - t503 * mrSges(6,3);
t469 = m(6) * t471 - qJDD(5) * mrSges(6,2) + t489 * mrSges(6,3) - qJD(5) * t496 + t502 * t487;
t460 = t534 * t468 + t532 * t469;
t475 = -t528 * t479 + t557;
t545 = mrSges(5,3) * qJDD(1) + t536 * (-t560 + t561);
t458 = m(5) * t475 + t545 * t528 + t460;
t551 = -t532 * t468 + t534 * t469;
t459 = m(5) * t476 - t545 * t530 + t551;
t456 = -t528 * t458 + t530 * t459;
t452 = m(4) * t482 - t536 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t456;
t481 = -t529 * t494 + t531 * t499;
t543 = qJDD(1) * pkin(3) + qJDD(4) - t481;
t478 = -t536 * qJ(4) + t543;
t474 = pkin(4) * t554 + (-qJ(4) + (-t518 - t566) * pkin(6)) * t536 + t543;
t540 = m(6) * t474 - t489 * mrSges(6,1) + t490 * mrSges(6,2) - t502 * t495 + t503 * t496;
t539 = -m(5) * t478 + qJDD(1) * t560 - t540 + (t536 * t566 + t559) * mrSges(5,3);
t463 = m(4) * t481 - t536 * mrSges(4,2) + (-mrSges(4,1) - t561) * qJDD(1) + t539;
t552 = t531 * t452 - t529 * t463;
t544 = m(3) * t500 + qJDD(1) * mrSges(3,3) + t552;
t443 = m(2) * t508 - qJDD(1) * mrSges(2,2) - t564 * t536 + t544;
t449 = t529 * t452 + t531 * t463;
t501 = -qJDD(1) * pkin(1) + t541;
t447 = m(3) * t501 - qJDD(1) * mrSges(3,1) - t536 * mrSges(3,3) + t449;
t444 = m(2) * t507 + qJDD(1) * mrSges(2,1) - t536 * mrSges(2,2) - t447;
t558 = t533 * t443 + t535 * t444;
t548 = -Ifges(5,5) * t528 - Ifges(5,6) * t530;
t556 = t536 * t548;
t553 = t535 * t443 - t533 * t444;
t550 = -Ifges(5,1) * t528 - Ifges(5,4) * t530;
t549 = -Ifges(5,4) * t528 - Ifges(5,2) * t530;
t455 = t530 * t458 + t528 * t459;
t454 = m(4) * t525 + t455;
t484 = Ifges(6,4) * t503 + Ifges(6,2) * t502 + Ifges(6,6) * qJD(5);
t485 = Ifges(6,1) * t503 + Ifges(6,4) * t502 + Ifges(6,5) * qJD(5);
t538 = mrSges(6,1) * t470 - mrSges(6,2) * t471 + Ifges(6,5) * t490 + Ifges(6,6) * t489 + Ifges(6,3) * qJDD(5) + t503 * t484 - t502 * t485;
t483 = Ifges(6,5) * t503 + Ifges(6,6) * t502 + Ifges(6,3) * qJD(5);
t461 = -mrSges(6,1) * t474 + mrSges(6,3) * t471 + Ifges(6,4) * t490 + Ifges(6,2) * t489 + Ifges(6,6) * qJDD(5) + qJD(5) * t485 - t503 * t483;
t462 = mrSges(6,2) * t474 - mrSges(6,3) * t470 + Ifges(6,1) * t490 + Ifges(6,4) * t489 + Ifges(6,5) * qJDD(5) - qJD(5) * t484 + t502 * t483;
t448 = -mrSges(5,1) * t478 + mrSges(5,3) * t476 - pkin(4) * t540 + pkin(6) * t551 + t549 * qJDD(1) + t534 * t461 + t532 * t462 + t528 * t556;
t450 = mrSges(5,2) * t478 - mrSges(5,3) * t475 - pkin(6) * t460 + t550 * qJDD(1) - t532 * t461 + t534 * t462 - t530 * t556;
t464 = mrSges(5,1) * t554 - t539;
t537 = -mrSges(3,1) * t501 - mrSges(4,1) * t481 - mrSges(2,2) * t508 - pkin(2) * t449 + pkin(3) * t464 - qJ(4) * t456 - t530 * t448 - t528 * t450 + qJ(2) * (-t536 * mrSges(3,1) + t544) - pkin(1) * t447 + mrSges(4,2) * t482 + mrSges(3,3) * t500 + mrSges(2,1) * t507 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t453 = -m(3) * g(3) - t454;
t439 = -mrSges(4,1) * t525 - mrSges(5,1) * t475 + mrSges(5,2) * t476 + mrSges(4,3) * t482 - pkin(3) * t455 - pkin(4) * t460 + (-Ifges(4,6) - t548) * qJDD(1) - t538 + (t528 * t549 - t530 * t550 + Ifges(4,5)) * t536;
t438 = mrSges(4,2) * t525 - mrSges(4,3) * t481 - Ifges(4,5) * qJDD(1) - t536 * Ifges(4,6) - qJ(4) * t455 - t528 * t448 + t530 * t450;
t437 = mrSges(3,2) * t501 - mrSges(2,3) * t507 - qJ(2) * t453 - qJ(3) * t449 + t531 * t438 - t529 * t439 - t562 * t536 + t563 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t436 = mrSges(3,2) * t500 + mrSges(2,3) * t508 - pkin(1) * t453 + pkin(2) * t454 + t564 * g(3) - qJ(3) * t552 + t562 * qJDD(1) - t529 * t438 - t531 * t439 + t563 * t536;
t1 = [-m(1) * g(1) + t553; -m(1) * g(2) + t558; (-m(1) - m(2) - m(3)) * g(3) - t454; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t558 - t533 * t436 + t535 * t437; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t553 + t535 * t436 + t533 * t437; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t537; t537; t447; t454; t464; t538;];
tauJB = t1;
