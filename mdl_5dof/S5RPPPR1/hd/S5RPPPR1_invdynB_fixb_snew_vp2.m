% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:57
% EndTime: 2019-12-05 17:29:01
% DurationCPUTime: 3.33s
% Computational Cost: add. (31387->236), mult. (72890->326), div. (0->0), fcn. (44132->10), ass. (0->114)
t520 = sin(qJ(1));
t522 = cos(qJ(1));
t497 = t522 * g(2) + t520 * g(3);
t492 = qJDD(1) * pkin(1) + t497;
t496 = t520 * g(2) - t522 * g(3);
t523 = qJD(1) ^ 2;
t493 = -t523 * pkin(1) + t496;
t515 = sin(pkin(7));
t518 = cos(pkin(7));
t475 = t515 * t492 + t518 * t493;
t563 = -t523 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t475;
t474 = t518 * t492 - t515 * t493;
t529 = -t523 * qJ(3) + qJDD(3) - t474;
t514 = sin(pkin(8));
t517 = cos(pkin(8));
t536 = -pkin(3) * t517 - qJ(4) * t514;
t554 = t514 * qJD(1);
t562 = (-pkin(2) + t536) * qJDD(1) + t529 - 0.2e1 * qJD(4) * t554;
t512 = -g(1) + qJDD(2);
t460 = t517 * t512 - t563 * t514;
t561 = pkin(6) * t514;
t560 = mrSges(4,2) * t514;
t559 = Ifges(4,6) * t517;
t511 = t514 ^ 2;
t558 = t511 * t523;
t513 = sin(pkin(9));
t557 = t513 * t514;
t516 = cos(pkin(9));
t556 = t514 * t516;
t461 = t514 * t512 + t563 * t517;
t490 = (-mrSges(4,1) * t517 + t560) * qJD(1);
t489 = t536 * qJD(1);
t553 = t517 * qJD(1);
t454 = t489 * t553 + t461;
t534 = -pkin(4) * t517 - pkin(6) * t556;
t555 = t562 * t516;
t447 = t534 * qJDD(1) + (-t454 + (-pkin(4) * t511 * t516 + t517 * t561) * t523) * t513 + t555;
t450 = t516 * t454 + t562 * t513;
t488 = t534 * qJD(1);
t549 = t513 ^ 2 * t558;
t551 = qJDD(1) * t513;
t448 = -pkin(4) * t549 + t488 * t553 - t551 * t561 + t450;
t519 = sin(qJ(5));
t521 = cos(qJ(5));
t445 = t521 * t447 - t519 * t448;
t531 = (-t513 * t521 - t516 * t519) * t514;
t479 = qJD(1) * t531;
t530 = (-t513 * t519 + t516 * t521) * t514;
t480 = qJD(1) * t530;
t464 = -t479 * mrSges(6,1) + t480 * mrSges(6,2);
t467 = t479 * qJD(5) + qJDD(1) * t530;
t499 = qJD(5) - t553;
t472 = -t499 * mrSges(6,2) + t479 * mrSges(6,3);
t550 = t517 * qJDD(1);
t498 = qJDD(5) - t550;
t443 = m(6) * t445 + t498 * mrSges(6,1) - t467 * mrSges(6,3) - t480 * t464 + t499 * t472;
t446 = t519 * t447 + t521 * t448;
t466 = -t480 * qJD(5) + qJDD(1) * t531;
t473 = t499 * mrSges(6,1) - t480 * mrSges(6,3);
t444 = m(6) * t446 - t498 * mrSges(6,2) + t466 * mrSges(6,3) + t479 * t464 - t499 * t473;
t435 = t521 * t443 + t519 * t444;
t449 = -t513 * t454 + t555;
t539 = mrSges(5,1) * t513 + mrSges(5,2) * t516;
t481 = t539 * t554;
t532 = mrSges(5,2) * t517 - mrSges(5,3) * t557;
t483 = t532 * qJD(1);
t533 = -mrSges(5,1) * t517 - mrSges(5,3) * t556;
t433 = m(5) * t449 + t533 * qJDD(1) + (-t481 * t556 - t483 * t517) * qJD(1) + t435;
t484 = t533 * qJD(1);
t542 = -t519 * t443 + t521 * t444;
t434 = m(5) * t450 + t532 * qJDD(1) + (-t481 * t557 + t484 * t517) * qJD(1) + t542;
t543 = -t513 * t433 + t516 * t434;
t430 = m(4) * t461 + (qJDD(1) * mrSges(4,3) + qJD(1) * t490) * t517 + t543;
t453 = t489 * t554 + qJDD(4) - t460;
t451 = -pkin(6) * t549 + (qJD(1) * t488 * t516 + pkin(4) * t551) * t514 + t453;
t528 = m(6) * t451 - t466 * mrSges(6,1) + t467 * mrSges(6,2) - t479 * t472 + t480 * t473;
t524 = -m(5) * t453 - t528;
t439 = m(4) * t460 + ((-mrSges(4,3) - t539) * qJDD(1) + (-t483 * t513 - t484 * t516 - t490) * qJD(1)) * t514 + t524;
t544 = t517 * t430 - t514 * t439;
t422 = m(3) * t475 - t523 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t544;
t431 = t516 * t433 + t513 * t434;
t470 = -qJDD(1) * pkin(2) + t529;
t525 = -m(4) * t470 + mrSges(4,1) * t550 - t431 + (t517 ^ 2 * t523 + t558) * mrSges(4,3);
t427 = m(3) * t474 - t523 * mrSges(3,2) + (mrSges(3,1) - t560) * qJDD(1) + t525;
t418 = t515 * t422 + t518 * t427;
t423 = t514 * t430 + t517 * t439;
t548 = m(3) * t512 + t423;
t545 = t518 * t422 - t515 * t427;
t416 = m(2) * t496 - t523 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t545;
t417 = m(2) * t497 + qJDD(1) * mrSges(2,1) - t523 * mrSges(2,2) + t418;
t546 = t522 * t416 - t520 * t417;
t538 = Ifges(4,1) * t514 + Ifges(4,4) * t517;
t537 = Ifges(5,5) * t516 - Ifges(5,6) * t513;
t535 = -t520 * t416 - t522 * t417;
t527 = -Ifges(5,5) * t517 + (Ifges(5,1) * t516 - Ifges(5,4) * t513) * t514;
t526 = -Ifges(5,6) * t517 + (Ifges(5,4) * t516 - Ifges(5,2) * t513) * t514;
t491 = (Ifges(4,5) * t514 + t559) * qJD(1);
t478 = t527 * qJD(1);
t477 = t526 * qJD(1);
t476 = (-Ifges(5,3) * t517 + t537 * t514) * qJD(1);
t457 = Ifges(6,1) * t480 + Ifges(6,4) * t479 + Ifges(6,5) * t499;
t456 = Ifges(6,4) * t480 + Ifges(6,2) * t479 + Ifges(6,6) * t499;
t455 = Ifges(6,5) * t480 + Ifges(6,6) * t479 + Ifges(6,3) * t499;
t437 = mrSges(6,2) * t451 - mrSges(6,3) * t445 + Ifges(6,1) * t467 + Ifges(6,4) * t466 + Ifges(6,5) * t498 + t479 * t455 - t499 * t456;
t436 = -mrSges(6,1) * t451 + mrSges(6,3) * t446 + Ifges(6,4) * t467 + Ifges(6,2) * t466 + Ifges(6,6) * t498 - t480 * t455 + t499 * t457;
t425 = mrSges(5,2) * t453 - mrSges(5,3) * t449 - pkin(6) * t435 - t519 * t436 + t521 * t437 + (-t476 * t557 + t477 * t517) * qJD(1) + t527 * qJDD(1);
t424 = -mrSges(5,1) * t453 + mrSges(5,3) * t450 + t519 * t437 + t521 * t436 - pkin(4) * t528 + pkin(6) * t542 + (-t476 * t556 - t517 * t478) * qJD(1) + t526 * qJDD(1);
t419 = -mrSges(4,1) * t470 - mrSges(5,1) * t449 - mrSges(6,1) * t445 + mrSges(5,2) * t450 + mrSges(6,2) * t446 + mrSges(4,3) * t461 - Ifges(6,5) * t467 - Ifges(6,6) * t466 - Ifges(6,3) * t498 - pkin(3) * t431 - pkin(4) * t435 - t480 * t456 + t479 * t457 + (Ifges(4,2) + Ifges(5,3)) * t550 + ((Ifges(4,4) - t537) * qJDD(1) + (-t477 * t516 - t478 * t513 - t491) * qJD(1)) * t514;
t414 = mrSges(4,2) * t470 - mrSges(4,3) * t460 - qJ(4) * t431 + t538 * qJDD(1) - t513 * t424 + t516 * t425 + t491 * t553;
t413 = t523 * Ifges(3,5) - mrSges(3,1) * t512 + mrSges(3,3) * t475 - mrSges(4,1) * t460 + mrSges(4,2) * t461 - t513 * t425 - t516 * t424 - pkin(3) * t524 - qJ(4) * t543 - pkin(2) * t423 + (-t559 + Ifges(3,6) + (pkin(3) * t539 - Ifges(4,5)) * t514) * qJDD(1) + (-pkin(3) * (-t483 * t557 - t484 * t556) + (-t514 * (Ifges(4,4) * t514 + Ifges(4,2) * t517) + t517 * t538) * qJD(1)) * qJD(1);
t412 = mrSges(3,2) * t512 - mrSges(3,3) * t474 + Ifges(3,5) * qJDD(1) - t523 * Ifges(3,6) - qJ(3) * t423 + t517 * t414 - t514 * t419;
t411 = -mrSges(2,2) * g(1) - mrSges(2,3) * t497 + Ifges(2,5) * qJDD(1) - t523 * Ifges(2,6) - qJ(2) * t418 + t518 * t412 - t515 * t413;
t410 = mrSges(2,1) * g(1) + mrSges(2,3) * t496 + t523 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t548 + qJ(2) * t545 + t515 * t412 + t518 * t413;
t1 = [(-m(1) - m(2)) * g(1) + t548; -m(1) * g(2) + t535; -m(1) * g(3) + t546; pkin(1) * t418 + qJ(3) * t544 + t514 * t414 + t517 * t419 + pkin(2) * t525 + mrSges(3,1) * t474 - mrSges(3,2) * t475 + mrSges(2,1) * t497 - mrSges(2,2) * t496 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(2) * t560 + Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t546 - t522 * t410 - t520 * t411; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t535 - t520 * t410 + t522 * t411;];
tauB = t1;
