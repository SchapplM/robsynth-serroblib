% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:32:04
% EndTime: 2019-12-05 18:32:07
% DurationCPUTime: 3.15s
% Computational Cost: add. (53053->230), mult. (68329->290), div. (0->0), fcn. (38611->10), ass. (0->99)
t591 = sin(qJ(1));
t595 = cos(qJ(1));
t568 = t595 * g(2) + t591 * g(3);
t561 = qJDD(1) * pkin(1) + t568;
t567 = t591 * g(2) - g(3) * t595;
t596 = qJD(1) ^ 2;
t562 = -pkin(1) * t596 + t567;
t590 = sin(qJ(2));
t594 = cos(qJ(2));
t544 = t594 * t561 - t562 * t590;
t580 = qJDD(1) + qJDD(2);
t541 = pkin(2) * t580 + t544;
t545 = t590 * t561 + t594 * t562;
t582 = qJD(1) + qJD(2);
t578 = t582 ^ 2;
t542 = -pkin(2) * t578 + t545;
t586 = sin(pkin(9));
t587 = cos(pkin(9));
t526 = t586 * t541 + t587 * t542;
t523 = -pkin(3) * t578 + pkin(7) * t580 + t526;
t585 = -g(1) + qJDD(3);
t589 = sin(qJ(4));
t593 = cos(qJ(4));
t519 = -t589 * t523 + t593 * t585;
t612 = qJD(4) * t582;
t611 = t593 * t612;
t556 = t580 * t589 + t611;
t516 = (-t556 + t611) * pkin(8) + (t578 * t589 * t593 + qJDD(4)) * pkin(4) + t519;
t520 = t593 * t523 + t589 * t585;
t557 = t580 * t593 - t589 * t612;
t614 = t582 * t589;
t565 = qJD(4) * pkin(4) - pkin(8) * t614;
t584 = t593 ^ 2;
t517 = -pkin(4) * t578 * t584 + pkin(8) * t557 - qJD(4) * t565 + t520;
t588 = sin(qJ(5));
t592 = cos(qJ(5));
t514 = t516 * t592 - t517 * t588;
t551 = (-t588 * t589 + t592 * t593) * t582;
t532 = qJD(5) * t551 + t556 * t592 + t557 * t588;
t552 = (t588 * t593 + t589 * t592) * t582;
t537 = -mrSges(6,1) * t551 + mrSges(6,2) * t552;
t581 = qJD(4) + qJD(5);
t546 = -mrSges(6,2) * t581 + mrSges(6,3) * t551;
t579 = qJDD(4) + qJDD(5);
t511 = m(6) * t514 + mrSges(6,1) * t579 - mrSges(6,3) * t532 - t537 * t552 + t546 * t581;
t515 = t516 * t588 + t517 * t592;
t531 = -qJD(5) * t552 - t556 * t588 + t557 * t592;
t547 = mrSges(6,1) * t581 - mrSges(6,3) * t552;
t512 = m(6) * t515 - mrSges(6,2) * t579 + mrSges(6,3) * t531 + t537 * t551 - t547 * t581;
t502 = t592 * t511 + t588 * t512;
t549 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t589 + Ifges(5,2) * t593) * t582;
t550 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t589 + Ifges(5,4) * t593) * t582;
t534 = Ifges(6,4) * t552 + Ifges(6,2) * t551 + Ifges(6,6) * t581;
t535 = Ifges(6,1) * t552 + Ifges(6,4) * t551 + Ifges(6,5) * t581;
t601 = -mrSges(6,1) * t514 + mrSges(6,2) * t515 - Ifges(6,5) * t532 - Ifges(6,6) * t531 - Ifges(6,3) * t579 - t552 * t534 + t551 * t535;
t615 = mrSges(5,1) * t519 - mrSges(5,2) * t520 + Ifges(5,5) * t556 + Ifges(5,6) * t557 + Ifges(5,3) * qJDD(4) + pkin(4) * t502 + (t549 * t589 - t550 * t593) * t582 - t601;
t613 = t582 * t593;
t555 = (-mrSges(5,1) * t593 + mrSges(5,2) * t589) * t582;
t564 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t613;
t500 = m(5) * t519 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t556 + qJD(4) * t564 - t555 * t614 + t502;
t563 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t614;
t606 = -t511 * t588 + t592 * t512;
t501 = m(5) * t520 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t557 - qJD(4) * t563 + t555 * t613 + t606;
t607 = -t500 * t589 + t593 * t501;
t493 = m(4) * t526 - mrSges(4,1) * t578 - mrSges(4,2) * t580 + t607;
t525 = t587 * t541 - t586 * t542;
t603 = -t580 * pkin(3) - t525;
t522 = -pkin(7) * t578 + t603;
t518 = t565 * t614 - t557 * pkin(4) + (-pkin(8) * t584 - pkin(7)) * t578 + t603;
t602 = m(6) * t518 - t531 * mrSges(6,1) + mrSges(6,2) * t532 - t551 * t546 + t547 * t552;
t598 = -m(5) * t522 + t557 * mrSges(5,1) - mrSges(5,2) * t556 - t563 * t614 + t564 * t613 - t602;
t506 = m(4) * t525 + mrSges(4,1) * t580 - mrSges(4,2) * t578 + t598;
t488 = t586 * t493 + t587 * t506;
t485 = m(3) * t544 + mrSges(3,1) * t580 - mrSges(3,2) * t578 + t488;
t608 = t587 * t493 - t506 * t586;
t486 = m(3) * t545 - mrSges(3,1) * t578 - mrSges(3,2) * t580 + t608;
t478 = t594 * t485 + t590 * t486;
t496 = t593 * t500 + t589 * t501;
t494 = m(4) * t585 + t496;
t609 = -t485 * t590 + t594 * t486;
t475 = m(2) * t567 - mrSges(2,1) * t596 - qJDD(1) * mrSges(2,2) + t609;
t476 = m(2) * t568 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t596 + t478;
t610 = t595 * t475 - t476 * t591;
t605 = -t475 * t591 - t476 * t595;
t533 = Ifges(6,5) * t552 + Ifges(6,6) * t551 + Ifges(6,3) * t581;
t503 = -mrSges(6,1) * t518 + mrSges(6,3) * t515 + Ifges(6,4) * t532 + Ifges(6,2) * t531 + Ifges(6,6) * t579 - t533 * t552 + t535 * t581;
t504 = mrSges(6,2) * t518 - mrSges(6,3) * t514 + Ifges(6,1) * t532 + Ifges(6,4) * t531 + Ifges(6,5) * t579 + t533 * t551 - t534 * t581;
t548 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t589 + Ifges(5,6) * t593) * t582;
t481 = -mrSges(5,1) * t522 + mrSges(5,3) * t520 + Ifges(5,4) * t556 + Ifges(5,2) * t557 + Ifges(5,6) * qJDD(4) - pkin(4) * t602 + pkin(8) * t606 + qJD(4) * t550 + t592 * t503 + t588 * t504 - t548 * t614;
t490 = mrSges(5,2) * t522 - mrSges(5,3) * t519 + Ifges(5,1) * t556 + Ifges(5,4) * t557 + Ifges(5,5) * qJDD(4) - pkin(8) * t502 - qJD(4) * t549 - t503 * t588 + t504 * t592 + t548 * t613;
t600 = mrSges(3,1) * t544 + mrSges(4,1) * t525 - mrSges(3,2) * t545 - mrSges(4,2) * t526 + pkin(2) * t488 + pkin(3) * t598 + pkin(7) * t607 + t593 * t481 + t589 * t490 + (Ifges(4,3) + Ifges(3,3)) * t580;
t599 = mrSges(2,1) * t568 - mrSges(2,2) * t567 + Ifges(2,3) * qJDD(1) + pkin(1) * t478 + t600;
t479 = -mrSges(4,1) * t585 + mrSges(4,3) * t526 + t578 * Ifges(4,5) + Ifges(4,6) * t580 - pkin(3) * t496 - t615;
t473 = mrSges(4,2) * t585 - mrSges(4,3) * t525 + Ifges(4,5) * t580 - Ifges(4,6) * t578 - pkin(7) * t496 - t481 * t589 + t490 * t593;
t472 = -mrSges(3,2) * g(1) - mrSges(3,3) * t544 + Ifges(3,5) * t580 - Ifges(3,6) * t578 - qJ(3) * t488 + t473 * t587 - t479 * t586;
t471 = mrSges(3,1) * g(1) + mrSges(3,3) * t545 + t578 * Ifges(3,5) + Ifges(3,6) * t580 - pkin(2) * t494 + qJ(3) * t608 + t586 * t473 + t587 * t479;
t470 = -mrSges(2,2) * g(1) - mrSges(2,3) * t568 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t596 - pkin(6) * t478 - t471 * t590 + t472 * t594;
t469 = Ifges(2,6) * qJDD(1) + t596 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t567 + t590 * t472 + t594 * t471 - pkin(1) * (-m(3) * g(1) + t494) + pkin(6) * t609;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t494; -m(1) * g(2) + t605; -m(1) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t599; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t610 - t595 * t469 - t591 * t470; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t605 - t591 * t469 + t595 * t470; t599; t600; t494; t615; -t601;];
tauJB = t1;
