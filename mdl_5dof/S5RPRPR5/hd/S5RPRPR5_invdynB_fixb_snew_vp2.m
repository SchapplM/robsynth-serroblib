% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:56:28
% EndTime: 2019-12-05 17:56:38
% DurationCPUTime: 6.68s
% Computational Cost: add. (65060->288), mult. (167253->380), div. (0->0), fcn. (114048->10), ass. (0->121)
t573 = sin(qJ(1));
t576 = cos(qJ(1));
t553 = t573 * g(2) - g(3) * t576;
t577 = qJD(1) ^ 2;
t608 = -pkin(1) * t577 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t553;
t568 = sin(pkin(8));
t570 = cos(pkin(8));
t525 = -t570 * g(1) - t608 * t568;
t600 = qJD(1) * t570;
t557 = qJD(3) - t600;
t572 = sin(qJ(3));
t601 = qJD(1) * t568;
t594 = t572 * t601;
t540 = -mrSges(4,2) * t557 - mrSges(4,3) * t594;
t575 = cos(qJ(3));
t593 = t575 * t601;
t542 = mrSges(4,1) * t557 - mrSges(4,3) * t593;
t607 = -t540 * t572 - t542 * t575;
t606 = mrSges(3,2) * t568;
t603 = t568 ^ 2 * t577;
t526 = -g(1) * t568 + t608 * t570;
t549 = (-mrSges(3,1) * t570 + t606) * qJD(1);
t586 = -pkin(2) * t570 - pkin(6) * t568;
t551 = t586 * qJD(1);
t515 = t551 * t600 + t526;
t554 = t576 * g(2) + t573 * g(3);
t581 = -t577 * qJ(2) + qJDD(2) - t554;
t527 = (-pkin(1) + t586) * qJDD(1) + t581;
t524 = t575 * t527;
t597 = qJD(1) * qJD(3);
t545 = (qJDD(1) * t575 - t572 * t597) * t568;
t596 = qJDD(1) * t570;
t556 = qJDD(3) - t596;
t494 = t556 * pkin(3) - t545 * qJ(4) + t524 + (-pkin(3) * t575 * t603 - qJ(4) * t557 * t601 - t515) * t572;
t504 = t575 * t515 + t572 * t527;
t541 = pkin(3) * t557 - qJ(4) * t593;
t544 = (-qJDD(1) * t572 - t575 * t597) * t568;
t595 = t572 ^ 2 * t603;
t495 = -pkin(3) * t595 + qJ(4) * t544 - t541 * t557 + t504;
t567 = sin(pkin(9));
t569 = cos(pkin(9));
t536 = (-t567 * t572 + t569 * t575) * t601;
t483 = -0.2e1 * qJD(4) * t536 + t569 * t494 - t567 * t495;
t519 = t544 * t567 + t545 * t569;
t535 = (-t567 * t575 - t569 * t572) * t601;
t481 = (t535 * t557 - t519) * pkin(7) + (t535 * t536 + t556) * pkin(4) + t483;
t484 = 0.2e1 * qJD(4) * t535 + t567 * t494 + t569 * t495;
t518 = t544 * t569 - t545 * t567;
t522 = pkin(4) * t557 - pkin(7) * t536;
t534 = t535 ^ 2;
t482 = -pkin(4) * t534 + pkin(7) * t518 - t522 * t557 + t484;
t571 = sin(qJ(5));
t574 = cos(qJ(5));
t479 = t481 * t574 - t482 * t571;
t512 = t535 * t574 - t536 * t571;
t490 = qJD(5) * t512 + t518 * t571 + t519 * t574;
t513 = t535 * t571 + t536 * t574;
t501 = -mrSges(6,1) * t512 + mrSges(6,2) * t513;
t555 = qJD(5) + t557;
t505 = -mrSges(6,2) * t555 + mrSges(6,3) * t512;
t552 = qJDD(5) + t556;
t477 = m(6) * t479 + mrSges(6,1) * t552 - mrSges(6,3) * t490 - t501 * t513 + t505 * t555;
t480 = t481 * t571 + t482 * t574;
t489 = -qJD(5) * t513 + t518 * t574 - t519 * t571;
t506 = mrSges(6,1) * t555 - mrSges(6,3) * t513;
t478 = m(6) * t480 - mrSges(6,2) * t552 + mrSges(6,3) * t489 + t501 * t512 - t506 * t555;
t469 = t574 * t477 + t571 * t478;
t516 = -mrSges(5,1) * t535 + mrSges(5,2) * t536;
t520 = -mrSges(5,2) * t557 + mrSges(5,3) * t535;
t467 = m(5) * t483 + mrSges(5,1) * t556 - mrSges(5,3) * t519 - t516 * t536 + t520 * t557 + t469;
t521 = mrSges(5,1) * t557 - mrSges(5,3) * t536;
t587 = -t477 * t571 + t574 * t478;
t468 = m(5) * t484 - mrSges(5,2) * t556 + mrSges(5,3) * t518 + t516 * t535 - t521 * t557 + t587;
t463 = t569 * t467 + t567 * t468;
t503 = -t572 * t515 + t524;
t543 = (mrSges(4,1) * t572 + mrSges(4,2) * t575) * t601;
t461 = m(4) * t503 + mrSges(4,1) * t556 - mrSges(4,3) * t545 + t540 * t557 - t543 * t593 + t463;
t588 = -t467 * t567 + t569 * t468;
t462 = m(4) * t504 - mrSges(4,2) * t556 + mrSges(4,3) * t544 - t542 * t557 - t543 * t594 + t588;
t589 = -t572 * t461 + t575 * t462;
t599 = qJDD(1) * mrSges(3,3);
t456 = m(3) * t526 + (qJD(1) * t549 + t599) * t570 + t589;
t514 = t551 * t601 - t525;
t502 = -pkin(3) * t544 - qJ(4) * t595 + t541 * t593 + qJDD(4) + t514;
t486 = -pkin(4) * t518 - pkin(7) * t534 + t522 * t536 + t502;
t583 = m(6) * t486 - t489 * mrSges(6,1) + t490 * mrSges(6,2) - t512 * t505 + t513 * t506;
t579 = m(5) * t502 - t518 * mrSges(5,1) + t519 * mrSges(5,2) - t535 * t520 + t536 * t521 + t583;
t578 = -m(4) * t514 + t544 * mrSges(4,1) - t545 * mrSges(4,2) - t579;
t476 = t578 + (-t599 + (-t549 + t607) * qJD(1)) * t568 + m(3) * t525;
t452 = t568 * t456 + t570 * t476;
t590 = t570 * t456 - t476 * t568;
t451 = m(2) * t553 - mrSges(2,1) * t577 - qJDD(1) * mrSges(2,2) + t590;
t457 = t575 * t461 + t572 * t462;
t547 = -qJDD(1) * pkin(1) + t581;
t580 = -m(3) * t547 + mrSges(3,1) * t596 - t457 + (t570 ^ 2 * t577 + t603) * mrSges(3,3);
t453 = m(2) * t554 - t577 * mrSges(2,2) + (mrSges(2,1) - t606) * qJDD(1) + t580;
t591 = t576 * t451 - t453 * t573;
t585 = Ifges(3,1) * t568 + Ifges(3,4) * t570;
t584 = Ifges(3,5) * t568 + Ifges(3,6) * t570;
t582 = -t451 * t573 - t453 * t576;
t550 = t584 * qJD(1);
t530 = Ifges(4,5) * t557 + (Ifges(4,1) * t575 - Ifges(4,4) * t572) * t601;
t529 = Ifges(4,6) * t557 + (Ifges(4,4) * t575 - Ifges(4,2) * t572) * t601;
t528 = Ifges(4,3) * t557 + (Ifges(4,5) * t575 - Ifges(4,6) * t572) * t601;
t509 = Ifges(5,1) * t536 + Ifges(5,4) * t535 + Ifges(5,5) * t557;
t508 = Ifges(5,4) * t536 + Ifges(5,2) * t535 + Ifges(5,6) * t557;
t507 = Ifges(5,5) * t536 + Ifges(5,6) * t535 + Ifges(5,3) * t557;
t498 = Ifges(6,1) * t513 + Ifges(6,4) * t512 + Ifges(6,5) * t555;
t497 = Ifges(6,4) * t513 + Ifges(6,2) * t512 + Ifges(6,6) * t555;
t496 = Ifges(6,5) * t513 + Ifges(6,6) * t512 + Ifges(6,3) * t555;
t471 = mrSges(6,2) * t486 - mrSges(6,3) * t479 + Ifges(6,1) * t490 + Ifges(6,4) * t489 + Ifges(6,5) * t552 + t496 * t512 - t497 * t555;
t470 = -mrSges(6,1) * t486 + mrSges(6,3) * t480 + Ifges(6,4) * t490 + Ifges(6,2) * t489 + Ifges(6,6) * t552 - t496 * t513 + t498 * t555;
t459 = mrSges(5,2) * t502 - mrSges(5,3) * t483 + Ifges(5,1) * t519 + Ifges(5,4) * t518 + Ifges(5,5) * t556 - pkin(7) * t469 - t470 * t571 + t471 * t574 + t507 * t535 - t508 * t557;
t458 = -mrSges(5,1) * t502 + mrSges(5,3) * t484 + Ifges(5,4) * t519 + Ifges(5,2) * t518 + Ifges(5,6) * t556 - pkin(4) * t583 + pkin(7) * t587 + t574 * t470 + t571 * t471 - t536 * t507 + t557 * t509;
t449 = mrSges(4,2) * t514 - mrSges(4,3) * t503 + Ifges(4,1) * t545 + Ifges(4,4) * t544 + Ifges(4,5) * t556 - qJ(4) * t463 - t458 * t567 + t459 * t569 - t528 * t594 - t529 * t557;
t448 = -mrSges(4,1) * t514 + mrSges(4,3) * t504 + Ifges(4,4) * t545 + Ifges(4,2) * t544 + Ifges(4,6) * t556 - pkin(3) * t579 + qJ(4) * t588 + t569 * t458 + t567 * t459 - t528 * t593 + t557 * t530;
t447 = -pkin(3) * t463 - pkin(2) * t457 + (Ifges(3,4) * qJDD(1) + (-t529 * t575 - t530 * t572 - t550) * qJD(1)) * t568 + (-Ifges(4,3) - Ifges(5,3)) * t556 - pkin(4) * t469 - mrSges(6,1) * t479 + mrSges(6,2) * t480 - mrSges(5,1) * t483 + mrSges(5,2) * t484 - Ifges(6,6) * t489 - Ifges(6,5) * t490 - mrSges(4,1) * t503 + mrSges(4,2) * t504 + t512 * t498 - t513 * t497 - Ifges(5,6) * t518 - Ifges(5,5) * t519 + mrSges(3,3) * t526 + t535 * t509 - t536 * t508 - Ifges(4,6) * t544 - Ifges(4,5) * t545 - mrSges(3,1) * t547 - Ifges(6,3) * t552 + Ifges(3,2) * t596;
t446 = mrSges(3,2) * t547 - mrSges(3,3) * t525 - pkin(6) * t457 + t585 * qJDD(1) - t572 * t448 + t575 * t449 + t550 * t600;
t445 = t577 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t553 - mrSges(3,1) * t525 + mrSges(3,2) * t526 - t572 * t449 - t575 * t448 - pkin(2) * t578 - pkin(6) * t589 - pkin(1) * t452 + (Ifges(2,6) - t584) * qJDD(1) + (-pkin(2) * t607 * t568 + (-t568 * (Ifges(3,4) * t568 + Ifges(3,2) * t570) + t570 * t585) * qJD(1)) * qJD(1);
t444 = -mrSges(2,2) * g(1) - mrSges(2,3) * t554 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t577 - qJ(2) * t452 + t446 * t570 - t447 * t568;
t1 = [(-m(1) - m(2)) * g(1) + t452; -m(1) * g(2) + t582; -m(1) * g(3) + t591; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t554 - mrSges(2,2) * t553 + t568 * t446 + t570 * t447 + pkin(1) * (-qJDD(1) * t606 + t580) + qJ(2) * t590; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t591 - t573 * t444 - t576 * t445; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t582 + t576 * t444 - t573 * t445;];
tauB = t1;
