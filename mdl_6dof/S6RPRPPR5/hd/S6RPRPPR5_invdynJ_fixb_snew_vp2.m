% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:58:10
% EndTime: 2019-05-05 16:58:15
% DurationCPUTime: 3.74s
% Computational Cost: add. (31778->301), mult. (79011->370), div. (0->0), fcn. (57142->10), ass. (0->129)
t574 = -2 * qJD(4);
t573 = Ifges(4,1) + Ifges(5,2);
t567 = Ifges(4,5) - Ifges(5,4);
t572 = -Ifges(4,2) - Ifges(5,3);
t566 = Ifges(4,6) - Ifges(5,5);
t565 = -Ifges(5,6) - Ifges(4,4);
t571 = (Ifges(4,3) + Ifges(5,1));
t533 = qJD(1) ^ 2;
t529 = sin(qJ(1));
t531 = cos(qJ(1));
t546 = -g(1) * t531 - g(2) * t529;
t510 = -pkin(1) * t533 + qJDD(1) * qJ(2) + t546;
t524 = sin(pkin(9));
t526 = cos(pkin(9));
t554 = qJD(1) * qJD(2);
t550 = -g(3) * t526 - 0.2e1 * t524 * t554;
t564 = pkin(7) * qJDD(1);
t569 = pkin(2) * t533;
t463 = (t526 * t569 - t510 - t564) * t524 + t550;
t489 = -g(3) * t524 + (t510 + 0.2e1 * t554) * t526;
t520 = t526 ^ 2;
t470 = -t520 * t569 + t526 * t564 + t489;
t528 = sin(qJ(3));
t570 = cos(qJ(3));
t447 = t528 * t463 + t570 * t470;
t552 = t526 * t570;
t557 = qJD(1) * t524;
t508 = -qJD(1) * t552 + t528 * t557;
t542 = t524 * t570 + t526 * t528;
t509 = t542 * qJD(1);
t479 = pkin(3) * t508 - qJ(4) * t509;
t532 = qJD(3) ^ 2;
t434 = pkin(3) * t532 - qJDD(3) * qJ(4) + (qJD(3) * t574) + t508 * t479 - t447;
t568 = mrSges(4,1) - mrSges(5,2);
t446 = t463 * t570 - t528 * t470;
t480 = mrSges(4,1) * t508 + mrSges(4,2) * t509;
t555 = t508 * qJD(3);
t487 = qJDD(1) * t542 - t555;
t556 = qJD(3) * t509;
t486 = t556 + (t524 * t528 - t552) * qJDD(1);
t499 = pkin(4) * t509 - (qJD(3) * qJ(5));
t507 = t508 ^ 2;
t551 = g(1) * t529 - t531 * g(2);
t545 = qJDD(2) - t551;
t558 = -t524 ^ 2 - t520;
t485 = (-pkin(2) * t526 - pkin(1)) * qJDD(1) + (pkin(7) * t558 - qJ(2)) * t533 + t545;
t534 = pkin(3) * t556 + t509 * t574 + (-t487 + t555) * qJ(4) + t485;
t425 = -pkin(4) * t507 - t499 * t509 + (pkin(3) + qJ(5)) * t486 + t534;
t438 = -qJDD(3) * pkin(3) - t532 * qJ(4) + t509 * t479 + qJDD(4) - t446;
t428 = (t508 * t509 - qJDD(3)) * qJ(5) + (t487 + t555) * pkin(4) + t438;
t523 = sin(pkin(10));
t525 = cos(pkin(10));
t494 = qJD(3) * t525 + t508 * t523;
t420 = -0.2e1 * qJD(5) * t494 - t425 * t523 + t525 * t428;
t469 = qJDD(3) * t525 + t486 * t523;
t493 = -qJD(3) * t523 + t508 * t525;
t418 = (t493 * t509 - t469) * pkin(8) + (t493 * t494 + t487) * pkin(5) + t420;
t421 = 0.2e1 * qJD(5) * t493 + t525 * t425 + t523 * t428;
t467 = pkin(5) * t509 - pkin(8) * t494;
t468 = -qJDD(3) * t523 + t486 * t525;
t492 = t493 ^ 2;
t419 = -pkin(5) * t492 + pkin(8) * t468 - t467 * t509 + t421;
t527 = sin(qJ(6));
t530 = cos(qJ(6));
t416 = t418 * t530 - t419 * t527;
t454 = t493 * t530 - t494 * t527;
t440 = qJD(6) * t454 + t468 * t527 + t469 * t530;
t455 = t493 * t527 + t494 * t530;
t445 = -mrSges(7,1) * t454 + mrSges(7,2) * t455;
t505 = qJD(6) + t509;
t448 = -mrSges(7,2) * t505 + mrSges(7,3) * t454;
t484 = qJDD(6) + t487;
t412 = m(7) * t416 + mrSges(7,1) * t484 - mrSges(7,3) * t440 - t445 * t455 + t448 * t505;
t417 = t418 * t527 + t419 * t530;
t439 = -qJD(6) * t455 + t468 * t530 - t469 * t527;
t449 = mrSges(7,1) * t505 - mrSges(7,3) * t455;
t413 = m(7) * t417 - mrSges(7,2) * t484 + mrSges(7,3) * t439 + t445 * t454 - t449 * t505;
t403 = t530 * t412 + t527 * t413;
t456 = -mrSges(6,1) * t493 + mrSges(6,2) * t494;
t465 = -mrSges(6,2) * t509 + mrSges(6,3) * t493;
t401 = m(6) * t420 + mrSges(6,1) * t487 - mrSges(6,3) * t469 - t456 * t494 + t465 * t509 + t403;
t466 = mrSges(6,1) * t509 - mrSges(6,3) * t494;
t547 = -t412 * t527 + t530 * t413;
t402 = m(6) * t421 - mrSges(6,2) * t487 + mrSges(6,3) * t468 + t456 * t493 - t466 * t509 + t547;
t399 = t401 * t525 + t402 * t523;
t481 = -mrSges(5,2) * t508 - mrSges(5,3) * t509;
t538 = -m(5) * t438 - t487 * mrSges(5,1) - t509 * t481 - t399;
t500 = mrSges(5,1) * t508 - qJD(3) * mrSges(5,3);
t559 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t508 - t500;
t396 = m(4) * t446 - mrSges(4,3) * t487 + qJD(3) * t559 + qJDD(3) * t568 - t480 * t509 + t538;
t498 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t509;
t430 = -pkin(4) * t486 - qJ(5) * t507 + qJD(3) * t499 + qJDD(5) - t434;
t423 = -pkin(5) * t468 - pkin(8) * t492 + t467 * t494 + t430;
t539 = m(7) * t423 - mrSges(7,1) * t439 + t440 * mrSges(7,2) - t448 * t454 + t455 * t449;
t414 = m(6) * t430 - mrSges(6,1) * t468 + t469 * mrSges(6,2) - t465 * t493 + t494 * t466 + t539;
t501 = mrSges(5,1) * t509 + qJD(3) * mrSges(5,2);
t536 = -m(5) * t434 + qJDD(3) * mrSges(5,3) + qJD(3) * t501 + t414;
t408 = (-t480 - t481) * t508 - qJDD(3) * mrSges(4,2) + m(4) * t447 - qJD(3) * t498 + t536 + (-mrSges(4,3) - mrSges(5,1)) * t486;
t563 = t570 * t396 + t528 * t408;
t562 = -(qJD(3) * t571) + t508 * t566 - t509 * t567;
t561 = qJD(3) * t566 + t508 * t572 - t509 * t565;
t560 = qJD(3) * t567 + t508 * t565 + t509 * t573;
t549 = -t528 * t396 + t570 * t408;
t548 = -t401 * t523 + t525 * t402;
t544 = -mrSges(3,1) * t526 + mrSges(3,2) * t524;
t543 = qJDD(1) * mrSges(3,3) + t533 * t544;
t433 = pkin(3) * t486 + t534;
t540 = m(5) * t433 - t487 * mrSges(5,3) - t509 * t501 + t548;
t442 = Ifges(7,4) * t455 + Ifges(7,2) * t454 + Ifges(7,6) * t505;
t443 = Ifges(7,1) * t455 + Ifges(7,4) * t454 + Ifges(7,5) * t505;
t537 = mrSges(7,1) * t416 - mrSges(7,2) * t417 + Ifges(7,5) * t440 + Ifges(7,6) * t439 + Ifges(7,3) * t484 + t455 * t442 - t454 * t443;
t535 = m(4) * t485 + mrSges(4,2) * t487 + t486 * t568 + t498 * t509 + t508 * t559 + t540;
t512 = (Ifges(3,5) * t524 + Ifges(3,6) * t526) * qJD(1);
t506 = -qJDD(1) * pkin(1) - qJ(2) * t533 + t545;
t488 = -t510 * t524 + t550;
t452 = Ifges(6,1) * t494 + Ifges(6,4) * t493 + Ifges(6,5) * t509;
t451 = Ifges(6,4) * t494 + Ifges(6,2) * t493 + Ifges(6,6) * t509;
t450 = Ifges(6,5) * t494 + Ifges(6,6) * t493 + Ifges(6,3) * t509;
t441 = Ifges(7,5) * t455 + Ifges(7,6) * t454 + Ifges(7,3) * t505;
t405 = mrSges(7,2) * t423 - mrSges(7,3) * t416 + Ifges(7,1) * t440 + Ifges(7,4) * t439 + Ifges(7,5) * t484 + t441 * t454 - t442 * t505;
t404 = -mrSges(7,1) * t423 + mrSges(7,3) * t417 + Ifges(7,4) * t440 + Ifges(7,2) * t439 + Ifges(7,6) * t484 - t441 * t455 + t443 * t505;
t398 = qJDD(3) * mrSges(5,2) + qJD(3) * t500 - t538;
t397 = -mrSges(5,2) * t486 - t500 * t508 + t540;
t394 = mrSges(3,3) * t533 * t558 + m(3) * t506 + qJDD(1) * t544 + t535;
t393 = mrSges(6,2) * t430 - mrSges(6,3) * t420 + Ifges(6,1) * t469 + Ifges(6,4) * t468 + Ifges(6,5) * t487 - pkin(8) * t403 - t404 * t527 + t405 * t530 + t450 * t493 - t451 * t509;
t392 = -mrSges(6,1) * t430 + mrSges(6,3) * t421 + Ifges(6,4) * t469 + Ifges(6,2) * t468 + Ifges(6,6) * t487 - pkin(5) * t539 + pkin(8) * t547 + t530 * t404 + t527 * t405 - t494 * t450 + t509 * t452;
t391 = t537 + t562 * t508 + (Ifges(6,3) + t573) * t487 + t565 * t486 - t561 * qJD(3) - t493 * t452 + t494 * t451 + mrSges(4,2) * t485 + Ifges(6,6) * t468 + Ifges(6,5) * t469 - mrSges(4,3) * t446 + mrSges(5,1) * t438 - mrSges(5,3) * t433 - mrSges(6,2) * t421 + mrSges(6,1) * t420 + t567 * qJDD(3) + pkin(5) * t403 - qJ(4) * t397 + pkin(4) * t399;
t390 = -mrSges(4,1) * t485 - mrSges(5,1) * t434 + mrSges(5,2) * t433 + mrSges(4,3) * t447 - pkin(3) * t397 + pkin(4) * t414 - qJ(5) * t548 + t560 * qJD(3) + t566 * qJDD(3) - t525 * t392 - t523 * t393 + t486 * t572 - t565 * t487 + t562 * t509;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t551 - mrSges(2,2) * t546 + t524 * (t526 * qJD(1) * t512 + mrSges(3,2) * t506 - mrSges(3,3) * t488 + t570 * t391 - t528 * t390 - pkin(7) * t563 + (Ifges(3,1) * t524 + Ifges(3,4) * t526) * qJDD(1)) + t526 * (-t512 * t557 - mrSges(3,1) * t506 + mrSges(3,3) * t489 + t528 * t391 + t570 * t390 - pkin(2) * t535 + pkin(7) * t549 + (Ifges(3,4) * t524 + Ifges(3,2) * t526) * qJDD(1)) - pkin(1) * t394 + qJ(2) * ((m(3) * t489 + t526 * t543 + t549) * t526 + (-m(3) * t488 + t524 * t543 - t563) * t524); t394; mrSges(4,1) * t446 - mrSges(4,2) * t447 + mrSges(5,2) * t438 - mrSges(5,3) * t434 + t525 * t393 - t523 * t392 - qJ(5) * t399 - pkin(3) * t398 + qJ(4) * t536 + t561 * t509 + (-qJ(4) * t481 + t560) * t508 + t567 * t487 + (-mrSges(5,1) * qJ(4) - t566) * t486 + t571 * qJDD(3); t398; t414; t537;];
tauJ  = t1;
