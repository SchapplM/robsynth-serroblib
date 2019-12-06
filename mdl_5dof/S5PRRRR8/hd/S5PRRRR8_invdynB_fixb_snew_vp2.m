% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:15:33
% EndTime: 2019-12-05 17:15:42
% DurationCPUTime: 6.42s
% Computational Cost: add. (78869->268), mult. (154414->352), div. (0->0), fcn. (107786->12), ass. (0->115)
t527 = sin(qJ(4));
t528 = sin(qJ(3));
t531 = cos(qJ(4));
t532 = cos(qJ(3));
t501 = (t527 * t528 - t531 * t532) * qJD(2);
t522 = sin(pkin(10));
t524 = cos(pkin(10));
t511 = t522 * g(1) - t524 * g(2);
t512 = -t524 * g(1) - t522 * g(2);
t521 = -g(3) + qJDD(1);
t523 = sin(pkin(5));
t525 = cos(pkin(5));
t529 = sin(qJ(2));
t533 = cos(qJ(2));
t484 = -t529 * t512 + (t511 * t525 + t521 * t523) * t533;
t534 = qJD(2) ^ 2;
t537 = -qJDD(2) * pkin(2) - t484;
t479 = -t534 * pkin(7) + t537;
t547 = qJD(2) * qJD(3);
t546 = t532 * t547;
t509 = t528 * qJDD(2) + t546;
t510 = t532 * qJDD(2) - t528 * t547;
t549 = qJD(2) * t528;
t513 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t549;
t548 = qJD(2) * t532;
t514 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t548;
t551 = t525 * t529;
t552 = t523 * t529;
t485 = t511 * t551 + t533 * t512 + t521 * t552;
t480 = -t534 * pkin(2) + qJDD(2) * pkin(7) + t485;
t495 = -t523 * t511 + t525 * t521;
t469 = -t528 * t480 + t532 * t495;
t459 = (-t509 + t546) * pkin(8) + (t528 * t532 * t534 + qJDD(3)) * pkin(3) + t469;
t470 = t532 * t480 + t528 * t495;
t516 = qJD(3) * pkin(3) - pkin(8) * t549;
t520 = t532 ^ 2;
t460 = -t520 * t534 * pkin(3) + t510 * pkin(8) - qJD(3) * t516 + t470;
t456 = t527 * t459 + t531 * t460;
t502 = (t527 * t532 + t528 * t531) * qJD(2);
t488 = t501 * pkin(4) - t502 * pkin(9);
t519 = qJD(3) + qJD(4);
t517 = t519 ^ 2;
t518 = qJDD(3) + qJDD(4);
t453 = -t517 * pkin(4) + t518 * pkin(9) - t501 * t488 + t456;
t464 = -t510 * pkin(3) + t516 * t549 + (-pkin(8) * t520 - pkin(7)) * t534 + t537;
t475 = -t502 * qJD(4) - t527 * t509 + t531 * t510;
t476 = -t501 * qJD(4) + t531 * t509 + t527 * t510;
t454 = (t501 * t519 - t476) * pkin(9) + (t502 * t519 - t475) * pkin(4) + t464;
t526 = sin(qJ(5));
t530 = cos(qJ(5));
t450 = -t526 * t453 + t530 * t454;
t489 = -t526 * t502 + t530 * t519;
t463 = t489 * qJD(5) + t530 * t476 + t526 * t518;
t490 = t530 * t502 + t526 * t519;
t471 = -t489 * mrSges(6,1) + t490 * mrSges(6,2);
t473 = qJDD(5) - t475;
t496 = qJD(5) + t501;
t477 = -t496 * mrSges(6,2) + t489 * mrSges(6,3);
t448 = m(6) * t450 + t473 * mrSges(6,1) - t463 * mrSges(6,3) - t490 * t471 + t496 * t477;
t451 = t530 * t453 + t526 * t454;
t462 = -t490 * qJD(5) - t526 * t476 + t530 * t518;
t478 = t496 * mrSges(6,1) - t490 * mrSges(6,3);
t449 = m(6) * t451 - t473 * mrSges(6,2) + t462 * mrSges(6,3) + t489 * t471 - t496 * t478;
t440 = t530 * t448 + t526 * t449;
t493 = -t519 * mrSges(5,2) - t501 * mrSges(5,3);
t494 = t519 * mrSges(5,1) - t502 * mrSges(5,3);
t536 = m(5) * t464 - t475 * mrSges(5,1) + t476 * mrSges(5,2) + t501 * t493 + t502 * t494 + t440;
t535 = -m(4) * t479 + t510 * mrSges(4,1) - t509 * mrSges(4,2) - t513 * t549 + t514 * t548 - t536;
t436 = m(3) * t484 + qJDD(2) * mrSges(3,1) - t534 * mrSges(3,2) + t535;
t553 = t436 * t533;
t487 = t501 * mrSges(5,1) + t502 * mrSges(5,2);
t542 = -t526 * t448 + t530 * t449;
t439 = m(5) * t456 - t518 * mrSges(5,2) + t475 * mrSges(5,3) - t501 * t487 - t519 * t494 + t542;
t455 = t531 * t459 - t527 * t460;
t452 = -t518 * pkin(4) - t517 * pkin(9) + t502 * t488 - t455;
t538 = -m(6) * t452 + t462 * mrSges(6,1) - t463 * mrSges(6,2) + t489 * t477 - t490 * t478;
t444 = m(5) * t455 + t518 * mrSges(5,1) - t476 * mrSges(5,3) - t502 * t487 + t519 * t493 + t538;
t433 = t527 * t439 + t531 * t444;
t508 = (-mrSges(4,1) * t532 + mrSges(4,2) * t528) * qJD(2);
t431 = m(4) * t469 + qJDD(3) * mrSges(4,1) - t509 * mrSges(4,3) + qJD(3) * t514 - t508 * t549 + t433;
t543 = t531 * t439 - t527 * t444;
t432 = m(4) * t470 - qJDD(3) * mrSges(4,2) + t510 * mrSges(4,3) - qJD(3) * t513 + t508 * t548 + t543;
t544 = -t528 * t431 + t532 * t432;
t422 = m(3) * t485 - t534 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t544;
t425 = t532 * t431 + t528 * t432;
t424 = m(3) * t495 + t425;
t412 = t422 * t551 - t523 * t424 + t525 * t553;
t410 = m(2) * t511 + t412;
t418 = t533 * t422 - t529 * t436;
t417 = m(2) * t512 + t418;
t550 = t524 * t410 + t522 * t417;
t411 = t422 * t552 + t525 * t424 + t523 * t553;
t545 = -t522 * t410 + t524 * t417;
t465 = Ifges(6,5) * t490 + Ifges(6,6) * t489 + Ifges(6,3) * t496;
t467 = Ifges(6,1) * t490 + Ifges(6,4) * t489 + Ifges(6,5) * t496;
t441 = -mrSges(6,1) * t452 + mrSges(6,3) * t451 + Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t473 - t490 * t465 + t496 * t467;
t466 = Ifges(6,4) * t490 + Ifges(6,2) * t489 + Ifges(6,6) * t496;
t442 = mrSges(6,2) * t452 - mrSges(6,3) * t450 + Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t473 + t489 * t465 - t496 * t466;
t481 = Ifges(5,5) * t502 - Ifges(5,6) * t501 + Ifges(5,3) * t519;
t482 = Ifges(5,4) * t502 - Ifges(5,2) * t501 + Ifges(5,6) * t519;
t426 = mrSges(5,2) * t464 - mrSges(5,3) * t455 + Ifges(5,1) * t476 + Ifges(5,4) * t475 + Ifges(5,5) * t518 - pkin(9) * t440 - t526 * t441 + t530 * t442 - t501 * t481 - t519 * t482;
t483 = Ifges(5,1) * t502 - Ifges(5,4) * t501 + Ifges(5,5) * t519;
t427 = -mrSges(5,1) * t464 - mrSges(6,1) * t450 + mrSges(6,2) * t451 + mrSges(5,3) * t456 + Ifges(5,4) * t476 - Ifges(6,5) * t463 + Ifges(5,2) * t475 + Ifges(5,6) * t518 - Ifges(6,6) * t462 - Ifges(6,3) * t473 - pkin(4) * t440 - t490 * t466 + t489 * t467 - t502 * t481 + t519 * t483;
t498 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t528 + Ifges(4,6) * t532) * qJD(2);
t500 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t528 + Ifges(4,4) * t532) * qJD(2);
t413 = -mrSges(4,1) * t479 + mrSges(4,3) * t470 + Ifges(4,4) * t509 + Ifges(4,2) * t510 + Ifges(4,6) * qJDD(3) - pkin(3) * t536 + pkin(8) * t543 + qJD(3) * t500 + t527 * t426 + t531 * t427 - t498 * t549;
t499 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t528 + Ifges(4,2) * t532) * qJD(2);
t414 = mrSges(4,2) * t479 - mrSges(4,3) * t469 + Ifges(4,1) * t509 + Ifges(4,4) * t510 + Ifges(4,5) * qJDD(3) - pkin(8) * t433 - qJD(3) * t499 + t531 * t426 - t527 * t427 + t498 * t548;
t407 = mrSges(3,2) * t495 - mrSges(3,3) * t484 + Ifges(3,5) * qJDD(2) - t534 * Ifges(3,6) - pkin(7) * t425 - t528 * t413 + t532 * t414;
t408 = -pkin(2) * t425 + Ifges(3,6) * qJDD(2) + mrSges(3,3) * t485 - mrSges(3,1) * t495 - pkin(3) * t433 - Ifges(4,5) * t509 - Ifges(4,6) * t510 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t469 + mrSges(4,2) * t470 - pkin(9) * t542 - Ifges(5,5) * t476 - Ifges(5,6) * t475 - Ifges(5,3) * t518 - mrSges(5,1) * t455 + mrSges(5,2) * t456 - t526 * t442 - t530 * t441 - pkin(4) * t538 + t534 * Ifges(3,5) - t502 * t482 - t501 * t483 + (-t528 * t499 + t532 * t500) * qJD(2);
t539 = pkin(6) * t418 + t407 * t529 + t408 * t533;
t406 = mrSges(3,1) * t484 - mrSges(3,2) * t485 + Ifges(3,3) * qJDD(2) + pkin(2) * t535 + pkin(7) * t544 + t532 * t413 + t528 * t414;
t405 = mrSges(2,2) * t521 - mrSges(2,3) * t511 + t533 * t407 - t529 * t408 + (-t411 * t523 - t412 * t525) * pkin(6);
t404 = -mrSges(2,1) * t521 + mrSges(2,3) * t512 - pkin(1) * t411 - t523 * t406 + t539 * t525;
t1 = [-m(1) * g(1) + t545; -m(1) * g(2) + t550; -m(1) * g(3) + m(2) * t521 + t411; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t550 - t522 * t404 + t524 * t405; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t545 + t524 * t404 + t522 * t405; -mrSges(1,1) * g(2) + mrSges(2,1) * t511 + mrSges(1,2) * g(1) - mrSges(2,2) * t512 + pkin(1) * t412 + t525 * t406 + t539 * t523;];
tauB = t1;
