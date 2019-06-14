% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 15:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:40:33
% EndTime: 2019-05-05 15:40:40
% DurationCPUTime: 4.77s
% Computational Cost: add. (67435->297), mult. (122628->359), div. (0->0), fcn. (67217->10), ass. (0->119)
t561 = -pkin(1) - pkin(2);
t560 = mrSges(2,1) + mrSges(3,1);
t559 = Ifges(3,4) + Ifges(2,5);
t558 = Ifges(2,6) - Ifges(3,6);
t533 = sin(qJ(1));
t537 = cos(qJ(1));
t514 = -t537 * g(1) - t533 * g(2);
t539 = qJD(1) ^ 2;
t546 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t514;
t495 = -pkin(1) * t539 + t546;
t492 = t561 * t539 + t546;
t513 = t533 * g(1) - t537 * g(2);
t545 = -t539 * qJ(2) + qJDD(2) - t513;
t494 = t561 * qJDD(1) + t545;
t528 = sin(pkin(10));
t529 = cos(pkin(10));
t473 = t529 * t492 + t528 * t494;
t471 = -pkin(3) * t539 - qJDD(1) * pkin(7) + t473;
t525 = g(3) + qJDD(3);
t532 = sin(qJ(4));
t536 = cos(qJ(4));
t467 = t536 * t471 + t532 * t525;
t507 = (mrSges(5,1) * t536 - mrSges(5,2) * t532) * qJD(1);
t555 = qJD(1) * qJD(4);
t554 = t532 * t555;
t510 = -t536 * qJDD(1) + t554;
t556 = qJD(1) * t532;
t511 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t556;
t519 = t536 * qJD(1);
t472 = -t528 * t492 + t529 * t494;
t470 = qJDD(1) * pkin(3) - t539 * pkin(7) - t472;
t553 = t536 * t555;
t509 = -qJDD(1) * t532 - t553;
t457 = (-t509 + t553) * pkin(8) + (-t510 - t554) * pkin(4) + t470;
t508 = (pkin(4) * t536 + pkin(8) * t532) * qJD(1);
t538 = qJD(4) ^ 2;
t460 = -pkin(4) * t538 + qJDD(4) * pkin(8) - t508 * t519 + t467;
t531 = sin(qJ(5));
t535 = cos(qJ(5));
t449 = t535 * t457 - t531 * t460;
t505 = qJD(4) * t535 + t531 * t556;
t482 = qJD(5) * t505 + qJDD(4) * t531 + t509 * t535;
t504 = qJDD(5) - t510;
t506 = qJD(4) * t531 - t535 * t556;
t516 = t519 + qJD(5);
t447 = (t505 * t516 - t482) * pkin(9) + (t505 * t506 + t504) * pkin(5) + t449;
t450 = t531 * t457 + t535 * t460;
t481 = -qJD(5) * t506 + qJDD(4) * t535 - t509 * t531;
t491 = pkin(5) * t516 - pkin(9) * t506;
t503 = t505 ^ 2;
t448 = -pkin(5) * t503 + pkin(9) * t481 - t491 * t516 + t450;
t530 = sin(qJ(6));
t534 = cos(qJ(6));
t445 = t447 * t534 - t448 * t530;
t483 = t505 * t534 - t506 * t530;
t454 = qJD(6) * t483 + t481 * t530 + t482 * t534;
t484 = t505 * t530 + t506 * t534;
t465 = -mrSges(7,1) * t483 + mrSges(7,2) * t484;
t515 = qJD(6) + t516;
t474 = -mrSges(7,2) * t515 + mrSges(7,3) * t483;
t500 = qJDD(6) + t504;
t443 = m(7) * t445 + mrSges(7,1) * t500 - mrSges(7,3) * t454 - t465 * t484 + t474 * t515;
t446 = t447 * t530 + t448 * t534;
t453 = -qJD(6) * t484 + t481 * t534 - t482 * t530;
t475 = mrSges(7,1) * t515 - mrSges(7,3) * t484;
t444 = m(7) * t446 - mrSges(7,2) * t500 + mrSges(7,3) * t453 + t465 * t483 - t475 * t515;
t436 = t534 * t443 + t530 * t444;
t485 = -mrSges(6,1) * t505 + mrSges(6,2) * t506;
t489 = -mrSges(6,2) * t516 + mrSges(6,3) * t505;
t434 = m(6) * t449 + mrSges(6,1) * t504 - mrSges(6,3) * t482 - t485 * t506 + t489 * t516 + t436;
t490 = mrSges(6,1) * t516 - mrSges(6,3) * t506;
t548 = -t443 * t530 + t534 * t444;
t435 = m(6) * t450 - mrSges(6,2) * t504 + mrSges(6,3) * t481 + t485 * t505 - t490 * t516 + t548;
t549 = -t434 * t531 + t535 * t435;
t431 = m(5) * t467 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t510 - qJD(4) * t511 - t507 * t519 + t549;
t466 = -t532 * t471 + t525 * t536;
t512 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t519;
t459 = -qJDD(4) * pkin(4) - pkin(8) * t538 - t508 * t556 - t466;
t451 = -pkin(5) * t481 - pkin(9) * t503 + t491 * t506 + t459;
t543 = m(7) * t451 - t453 * mrSges(7,1) + mrSges(7,2) * t454 - t483 * t474 + t475 * t484;
t540 = -m(6) * t459 + t481 * mrSges(6,1) - mrSges(6,2) * t482 + t505 * t489 - t490 * t506 - t543;
t439 = m(5) * t466 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t509 + qJD(4) * t512 + t507 * t556 + t540;
t550 = t536 * t431 - t439 * t532;
t424 = m(4) * t473 - mrSges(4,1) * t539 + qJDD(1) * mrSges(4,2) + t550;
t432 = t434 * t535 + t435 * t531;
t541 = -m(5) * t470 + t510 * mrSges(5,1) - mrSges(5,2) * t509 + t511 * t556 - t512 * t519 - t432;
t429 = m(4) * t472 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t539 + t541;
t551 = t529 * t424 - t528 * t429;
t547 = m(3) * t495 + qJDD(1) * mrSges(3,3) + t551;
t419 = m(2) * t514 - qJDD(1) * mrSges(2,2) - t560 * t539 + t547;
t421 = t424 * t528 + t429 * t529;
t496 = -qJDD(1) * pkin(1) + t545;
t542 = -m(3) * t496 + qJDD(1) * mrSges(3,1) + t539 * mrSges(3,3) - t421;
t420 = m(2) * t513 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t539 + t542;
t557 = t533 * t419 + t537 * t420;
t552 = t537 * t419 - t420 * t533;
t427 = t532 * t431 + t536 * t439;
t544 = -m(4) * t525 - t427;
t499 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t532 - Ifges(5,4) * t536) * qJD(1);
t498 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t532 - Ifges(5,2) * t536) * qJD(1);
t497 = (Ifges(5,3) * qJD(4)) + (-Ifges(5,5) * t532 - Ifges(5,6) * t536) * qJD(1);
t478 = Ifges(6,1) * t506 + Ifges(6,4) * t505 + Ifges(6,5) * t516;
t477 = Ifges(6,4) * t506 + Ifges(6,2) * t505 + Ifges(6,6) * t516;
t476 = Ifges(6,5) * t506 + Ifges(6,6) * t505 + Ifges(6,3) * t516;
t463 = Ifges(7,1) * t484 + Ifges(7,4) * t483 + Ifges(7,5) * t515;
t462 = Ifges(7,4) * t484 + Ifges(7,2) * t483 + Ifges(7,6) * t515;
t461 = Ifges(7,5) * t484 + Ifges(7,6) * t483 + Ifges(7,3) * t515;
t438 = mrSges(7,2) * t451 - mrSges(7,3) * t445 + Ifges(7,1) * t454 + Ifges(7,4) * t453 + Ifges(7,5) * t500 + t461 * t483 - t462 * t515;
t437 = -mrSges(7,1) * t451 + mrSges(7,3) * t446 + Ifges(7,4) * t454 + Ifges(7,2) * t453 + Ifges(7,6) * t500 - t461 * t484 + t463 * t515;
t428 = mrSges(6,2) * t459 - mrSges(6,3) * t449 + Ifges(6,1) * t482 + Ifges(6,4) * t481 + Ifges(6,5) * t504 - pkin(9) * t436 - t437 * t530 + t438 * t534 + t476 * t505 - t477 * t516;
t426 = -m(3) * g(3) + t544;
t425 = -mrSges(6,1) * t459 + mrSges(6,3) * t450 + Ifges(6,4) * t482 + Ifges(6,2) * t481 + Ifges(6,6) * t504 - pkin(5) * t543 + pkin(9) * t548 + t534 * t437 + t530 * t438 - t506 * t476 + t516 * t478;
t422 = Ifges(5,4) * t509 + Ifges(5,2) * t510 + Ifges(5,6) * qJDD(4) + t497 * t556 + qJD(4) * t499 - mrSges(5,1) * t470 + mrSges(5,3) * t467 - Ifges(6,5) * t482 - Ifges(6,6) * t481 - Ifges(6,3) * t504 - t506 * t477 + t505 * t478 - mrSges(6,1) * t449 + mrSges(6,2) * t450 - Ifges(7,5) * t454 - Ifges(7,6) * t453 - Ifges(7,3) * t500 - t484 * t462 + t483 * t463 - mrSges(7,1) * t445 + mrSges(7,2) * t446 - pkin(5) * t436 - pkin(4) * t432;
t415 = mrSges(5,2) * t470 - mrSges(5,3) * t466 + Ifges(5,1) * t509 + Ifges(5,4) * t510 + Ifges(5,5) * qJDD(4) - pkin(8) * t432 - qJD(4) * t498 - t425 * t531 + t428 * t535 - t497 * t519;
t414 = -Ifges(4,6) * qJDD(1) + t539 * Ifges(4,5) - mrSges(4,1) * t525 + mrSges(4,3) * t473 - Ifges(5,5) * t509 - Ifges(5,6) * t510 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t466 + mrSges(5,2) * t467 - t531 * t428 - t535 * t425 - pkin(4) * t540 - pkin(8) * t549 - pkin(3) * t427 + (t498 * t532 - t499 * t536) * qJD(1);
t413 = mrSges(4,2) * t525 - mrSges(4,3) * t472 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t539 - pkin(7) * t427 + t415 * t536 - t422 * t532;
t412 = mrSges(3,2) * t496 - mrSges(2,3) * t513 - qJ(2) * t426 - qJ(3) * t421 + t529 * t413 - t528 * t414 - t558 * t539 + t559 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t411 = mrSges(3,2) * t495 + mrSges(2,3) * t514 - pkin(1) * t426 - pkin(2) * t544 + t560 * g(3) - qJ(3) * t551 + t558 * qJDD(1) - t528 * t413 - t529 * t414 + t559 * t539;
t1 = [-m(1) * g(1) + t552; -m(1) * g(2) + t557; (-m(1) - m(2) - m(3)) * g(3) + t544; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t557 - t533 * t411 + t537 * t412; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t552 + t537 * t411 + t533 * t412; pkin(1) * t542 + qJ(2) * (-mrSges(3,1) * t539 + t547) + mrSges(2,1) * t513 - mrSges(2,2) * t514 - pkin(2) * t421 - mrSges(3,1) * t496 + mrSges(3,3) * t495 - t536 * t422 - pkin(3) * t541 - pkin(7) * t550 - t532 * t415 - mrSges(4,1) * t472 + mrSges(4,2) * t473 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB  = t1;
