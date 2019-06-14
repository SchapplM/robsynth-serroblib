% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-05-05 14:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:21:14
% EndTime: 2019-05-05 14:21:18
% DurationCPUTime: 2.72s
% Computational Cost: add. (29647->288), mult. (59477->343), div. (0->0), fcn. (33437->8), ass. (0->106)
t559 = sin(qJ(1));
t562 = cos(qJ(1));
t538 = t559 * g(1) - t562 * g(2);
t564 = qJD(1) ^ 2;
t517 = -qJDD(1) * pkin(1) - t564 * qJ(2) + qJDD(2) - t538;
t508 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t517;
t539 = -t562 * g(1) - t559 * g(2);
t589 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t539;
t588 = -m(3) - m(4);
t587 = mrSges(2,1) - mrSges(3,2);
t516 = pkin(1) * t564 - t589;
t509 = qJDD(3) + (-pkin(1) - qJ(3)) * t564 + t589;
t506 = -qJDD(1) * pkin(7) + t509;
t558 = sin(qJ(4));
t561 = cos(qJ(4));
t499 = -g(3) * t561 + t558 * t506;
t533 = (mrSges(5,1) * t558 + mrSges(5,2) * t561) * qJD(1);
t583 = qJD(1) * qJD(4);
t577 = t561 * t583;
t534 = qJDD(1) * t558 + t577;
t584 = qJD(1) * t561;
t537 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t584;
t505 = -t564 * pkin(7) - t508;
t578 = t558 * t583;
t535 = qJDD(1) * t561 - t578;
t488 = (-t535 + t578) * qJ(5) + (t534 + t577) * pkin(4) + t505;
t532 = (pkin(4) * t558 - qJ(5) * t561) * qJD(1);
t563 = qJD(4) ^ 2;
t585 = qJD(1) * t558;
t491 = -pkin(4) * t563 + qJDD(4) * qJ(5) - t532 * t585 + t499;
t555 = sin(pkin(9));
t556 = cos(pkin(9));
t528 = qJD(4) * t555 + t556 * t584;
t475 = -0.2e1 * qJD(5) * t528 + t556 * t488 - t555 * t491;
t514 = qJDD(4) * t555 + t535 * t556;
t527 = qJD(4) * t556 - t555 * t584;
t473 = (t527 * t585 - t514) * pkin(8) + (t527 * t528 + t534) * pkin(5) + t475;
t476 = 0.2e1 * qJD(5) * t527 + t555 * t488 + t556 * t491;
t513 = qJDD(4) * t556 - t535 * t555;
t515 = pkin(5) * t585 - pkin(8) * t528;
t526 = t527 ^ 2;
t474 = -pkin(5) * t526 + pkin(8) * t513 - t515 * t585 + t476;
t557 = sin(qJ(6));
t560 = cos(qJ(6));
t471 = t473 * t560 - t474 * t557;
t500 = t527 * t560 - t528 * t557;
t480 = qJD(6) * t500 + t513 * t557 + t514 * t560;
t501 = t527 * t557 + t528 * t560;
t485 = -mrSges(7,1) * t500 + mrSges(7,2) * t501;
t540 = qJD(6) + t585;
t492 = -mrSges(7,2) * t540 + mrSges(7,3) * t500;
t531 = qJDD(6) + t534;
t469 = m(7) * t471 + mrSges(7,1) * t531 - mrSges(7,3) * t480 - t485 * t501 + t492 * t540;
t472 = t473 * t557 + t474 * t560;
t479 = -qJD(6) * t501 + t513 * t560 - t514 * t557;
t493 = mrSges(7,1) * t540 - mrSges(7,3) * t501;
t470 = m(7) * t472 - mrSges(7,2) * t531 + mrSges(7,3) * t479 + t485 * t500 - t493 * t540;
t461 = t560 * t469 + t557 * t470;
t502 = -mrSges(6,1) * t527 + mrSges(6,2) * t528;
t511 = -mrSges(6,2) * t585 + mrSges(6,3) * t527;
t459 = m(6) * t475 + mrSges(6,1) * t534 - mrSges(6,3) * t514 - t502 * t528 + t511 * t585 + t461;
t512 = mrSges(6,1) * t585 - mrSges(6,3) * t528;
t573 = -t469 * t557 + t560 * t470;
t460 = m(6) * t476 - mrSges(6,2) * t534 + mrSges(6,3) * t513 + t502 * t527 - t512 * t585 + t573;
t574 = -t459 * t555 + t556 * t460;
t454 = m(5) * t499 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t534 - qJD(4) * t537 - t533 * t585 + t574;
t498 = g(3) * t558 + t506 * t561;
t536 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t585;
t490 = -qJDD(4) * pkin(4) - qJ(5) * t563 + t532 * t584 + qJDD(5) - t498;
t477 = -pkin(5) * t513 - pkin(8) * t526 + t515 * t528 + t490;
t567 = m(7) * t477 - t479 * mrSges(7,1) + mrSges(7,2) * t480 - t500 * t492 + t493 * t501;
t565 = -m(6) * t490 + t513 * mrSges(6,1) - mrSges(6,2) * t514 + t527 * t511 - t512 * t528 - t567;
t465 = m(5) * t498 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t535 + qJD(4) * t536 - t533 * t584 + t565;
t447 = t558 * t454 + t561 * t465;
t572 = -m(4) * t509 - qJDD(1) * mrSges(4,2) - t447;
t568 = -m(3) * t516 + t564 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t572;
t444 = m(2) * t539 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t564 + t568;
t455 = t556 * t459 + t555 * t460;
t569 = -m(5) * t505 - t534 * mrSges(5,1) - t535 * mrSges(5,2) - t536 * t585 - t537 * t584 - t455;
t451 = m(4) * t508 - t564 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t569;
t566 = -m(3) * t517 + t564 * mrSges(3,3) - t451;
t450 = m(2) * t538 - t564 * mrSges(2,2) + t587 * qJDD(1) + t566;
t586 = t559 * t444 + t562 * t450;
t580 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t579 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t576 = t562 * t444 - t450 * t559;
t575 = t561 * t454 - t558 * t465;
t523 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t561 - Ifges(5,4) * t558) * qJD(1);
t522 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t561 - Ifges(5,2) * t558) * qJD(1);
t521 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t561 - Ifges(5,6) * t558) * qJD(1);
t496 = Ifges(6,1) * t528 + Ifges(6,4) * t527 + Ifges(6,5) * t585;
t495 = Ifges(6,4) * t528 + Ifges(6,2) * t527 + Ifges(6,6) * t585;
t494 = Ifges(6,5) * t528 + Ifges(6,6) * t527 + Ifges(6,3) * t585;
t483 = Ifges(7,1) * t501 + Ifges(7,4) * t500 + Ifges(7,5) * t540;
t482 = Ifges(7,4) * t501 + Ifges(7,2) * t500 + Ifges(7,6) * t540;
t481 = Ifges(7,5) * t501 + Ifges(7,6) * t500 + Ifges(7,3) * t540;
t463 = mrSges(7,2) * t477 - mrSges(7,3) * t471 + Ifges(7,1) * t480 + Ifges(7,4) * t479 + Ifges(7,5) * t531 + t481 * t500 - t482 * t540;
t462 = -mrSges(7,1) * t477 + mrSges(7,3) * t472 + Ifges(7,4) * t480 + Ifges(7,2) * t479 + Ifges(7,6) * t531 - t481 * t501 + t483 * t540;
t448 = mrSges(6,2) * t490 - mrSges(6,3) * t475 + Ifges(6,1) * t514 + Ifges(6,4) * t513 + Ifges(6,5) * t534 - pkin(8) * t461 - t462 * t557 + t463 * t560 + t494 * t527 - t495 * t585;
t446 = t588 * g(3) + t575;
t445 = -mrSges(6,1) * t490 + mrSges(6,3) * t476 + Ifges(6,4) * t514 + Ifges(6,2) * t513 + Ifges(6,6) * t534 - pkin(5) * t567 + pkin(8) * t573 + t560 * t462 + t557 * t463 - t528 * t494 + t496 * t585;
t441 = Ifges(5,4) * t535 + Ifges(5,6) * qJDD(4) - t521 * t584 + qJD(4) * t523 - mrSges(5,1) * t505 + mrSges(5,3) * t499 - Ifges(6,5) * t514 - Ifges(6,6) * t513 - t528 * t495 + t527 * t496 - mrSges(6,1) * t475 + mrSges(6,2) * t476 - Ifges(7,5) * t480 - Ifges(7,6) * t479 - Ifges(7,3) * t531 - t501 * t482 + t500 * t483 - mrSges(7,1) * t471 + mrSges(7,2) * t472 - pkin(5) * t461 - pkin(4) * t455 + (-Ifges(5,2) - Ifges(6,3)) * t534;
t440 = mrSges(5,2) * t505 - mrSges(5,3) * t498 + Ifges(5,1) * t535 - Ifges(5,4) * t534 + Ifges(5,5) * qJDD(4) - qJ(5) * t455 - qJD(4) * t522 - t445 * t555 + t448 * t556 - t521 * t585;
t439 = -pkin(2) * t572 + Ifges(5,3) * qJDD(4) + qJ(5) * t574 + t555 * t448 + t556 * t445 - qJ(3) * t575 - Ifges(5,6) * t534 + Ifges(5,5) * t535 + mrSges(2,3) * t539 + pkin(4) * t565 - mrSges(3,1) * t516 + mrSges(4,1) * t509 + mrSges(5,1) * t498 - mrSges(5,2) * t499 + pkin(3) * t447 - pkin(1) * t446 + (t522 * t561 + t523 * t558) * qJD(1) + (-mrSges(4,3) * pkin(2) + t580) * t564 + t579 * qJDD(1) + (m(4) * qJ(3) + mrSges(4,3) + t587) * g(3);
t438 = -qJ(2) * t446 - mrSges(2,3) * t538 + pkin(2) * t451 + mrSges(3,1) * t517 + pkin(7) * t575 + t558 * t440 + t561 * t441 + pkin(3) * t569 + mrSges(4,1) * t508 - t579 * t564 + t580 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t576; -m(1) * g(2) + t586; (-m(1) - m(2) + t588) * g(3) + t575; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t586 + t562 * t438 - t559 * t439; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t576 + t559 * t438 + t562 * t439; qJ(2) * (-mrSges(4,3) * t564 + t568) + pkin(1) * t566 + mrSges(2,1) * t538 - mrSges(2,2) * t539 - qJ(3) * t451 + mrSges(3,2) * t517 - mrSges(3,3) * t516 - t558 * t441 - pkin(7) * t447 + t561 * t440 + mrSges(4,2) * t509 - mrSges(4,3) * t508 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
