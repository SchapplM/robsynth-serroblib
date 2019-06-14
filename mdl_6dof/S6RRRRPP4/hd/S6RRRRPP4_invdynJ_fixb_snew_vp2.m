% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:15:46
% EndTime: 2019-05-07 18:15:58
% DurationCPUTime: 7.46s
% Computational Cost: add. (79211->325), mult. (163352->398), div. (0->0), fcn. (118111->10), ass. (0->129)
t560 = Ifges(6,1) + Ifges(7,1);
t551 = Ifges(6,4) - Ifges(7,5);
t558 = Ifges(7,4) + Ifges(6,5);
t559 = Ifges(6,2) + Ifges(7,3);
t556 = Ifges(6,6) - Ifges(7,6);
t557 = -Ifges(7,2) - Ifges(6,3);
t523 = sin(qJ(3));
t527 = cos(qJ(3));
t524 = sin(qJ(2));
t547 = qJD(1) * t524;
t507 = qJD(2) * t527 - t523 * t547;
t508 = qJD(2) * t523 + t527 * t547;
t522 = sin(qJ(4));
t526 = cos(qJ(4));
t483 = t507 * t526 - t508 * t522;
t484 = t507 * t522 + t508 * t526;
t521 = sin(pkin(10));
t550 = cos(pkin(10));
t461 = -t483 * t550 + t521 * t484;
t462 = t521 * t483 + t484 * t550;
t528 = cos(qJ(2));
t546 = qJD(1) * t528;
t517 = qJD(3) - t546;
t516 = qJD(4) + t517;
t555 = t559 * t461 - t551 * t462 - t556 * t516;
t554 = -t551 * t461 + t560 * t462 + t558 * t516;
t553 = -2 * qJD(5);
t552 = -mrSges(6,3) - mrSges(7,2);
t531 = qJD(1) ^ 2;
t525 = sin(qJ(1));
t529 = cos(qJ(1));
t542 = t525 * g(1) - t529 * g(2);
t502 = -qJDD(1) * pkin(1) - t531 * pkin(7) - t542;
t545 = qJD(1) * qJD(2);
t543 = t528 * t545;
t511 = qJDD(1) * t524 + t543;
t518 = t524 * t545;
t512 = qJDD(1) * t528 - t518;
t466 = (-t511 - t543) * pkin(8) + (-t512 + t518) * pkin(2) + t502;
t537 = -g(1) * t529 - g(2) * t525;
t503 = -pkin(1) * t531 + qJDD(1) * pkin(7) + t537;
t489 = -g(3) * t524 + t528 * t503;
t510 = (-pkin(2) * t528 - pkin(8) * t524) * qJD(1);
t530 = qJD(2) ^ 2;
t469 = -pkin(2) * t530 + qJDD(2) * pkin(8) + t510 * t546 + t489;
t447 = t527 * t466 - t523 * t469;
t481 = qJD(3) * t507 + qJDD(2) * t523 + t511 * t527;
t506 = qJDD(3) - t512;
t423 = (t507 * t517 - t481) * pkin(9) + (t507 * t508 + t506) * pkin(3) + t447;
t448 = t523 * t466 + t527 * t469;
t480 = -qJD(3) * t508 + qJDD(2) * t527 - t511 * t523;
t490 = pkin(3) * t517 - pkin(9) * t508;
t505 = t507 ^ 2;
t425 = -pkin(3) * t505 + pkin(9) * t480 - t490 * t517 + t448;
t405 = t526 * t423 - t522 * t425;
t446 = qJD(4) * t483 + t480 * t522 + t481 * t526;
t504 = qJDD(4) + t506;
t401 = (t483 * t516 - t446) * qJ(5) + (t483 * t484 + t504) * pkin(4) + t405;
t406 = t522 * t423 + t526 * t425;
t445 = -qJD(4) * t484 + t480 * t526 - t481 * t522;
t471 = pkin(4) * t516 - qJ(5) * t484;
t482 = t483 ^ 2;
t403 = -pkin(4) * t482 + qJ(5) * t445 - t471 * t516 + t406;
t397 = t521 * t401 + t550 * t403 + t461 * t553;
t419 = -t445 * t550 + t521 * t446;
t453 = mrSges(6,1) * t516 - mrSges(6,3) * t462;
t436 = pkin(5) * t461 - qJ(6) * t462;
t515 = t516 ^ 2;
t393 = -pkin(5) * t515 + qJ(6) * t504 + 0.2e1 * qJD(6) * t516 - t436 * t461 + t397;
t454 = -mrSges(7,1) * t516 + mrSges(7,2) * t462;
t544 = m(7) * t393 + t504 * mrSges(7,3) + t516 * t454;
t437 = mrSges(7,1) * t461 - mrSges(7,3) * t462;
t548 = -mrSges(6,1) * t461 - mrSges(6,2) * t462 - t437;
t382 = m(6) * t397 - t504 * mrSges(6,2) + t419 * t552 - t516 * t453 + t461 * t548 + t544;
t536 = t401 * t550 - t521 * t403;
t396 = t462 * t553 + t536;
t420 = t521 * t445 + t446 * t550;
t452 = -mrSges(6,2) * t516 - mrSges(6,3) * t461;
t394 = -t504 * pkin(5) - t515 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t436) * t462 - t536;
t451 = -mrSges(7,2) * t461 + mrSges(7,3) * t516;
t538 = -m(7) * t394 + t504 * mrSges(7,1) + t516 * t451;
t384 = m(6) * t396 + t504 * mrSges(6,1) + t420 * t552 + t516 * t452 + t462 * t548 + t538;
t377 = t521 * t382 + t550 * t384;
t463 = -mrSges(5,1) * t483 + mrSges(5,2) * t484;
t470 = -mrSges(5,2) * t516 + mrSges(5,3) * t483;
t374 = m(5) * t405 + mrSges(5,1) * t504 - mrSges(5,3) * t446 - t463 * t484 + t470 * t516 + t377;
t472 = mrSges(5,1) * t516 - mrSges(5,3) * t484;
t539 = t550 * t382 - t384 * t521;
t375 = m(5) * t406 - mrSges(5,2) * t504 + mrSges(5,3) * t445 + t463 * t483 - t472 * t516 + t539;
t370 = t526 * t374 + t522 * t375;
t549 = t556 * t461 - t558 * t462 + t557 * t516;
t488 = -t528 * g(3) - t524 * t503;
t485 = -mrSges(4,1) * t507 + mrSges(4,2) * t508;
t486 = -mrSges(4,2) * t517 + mrSges(4,3) * t507;
t368 = m(4) * t447 + mrSges(4,1) * t506 - mrSges(4,3) * t481 - t485 * t508 + t486 * t517 + t370;
t487 = mrSges(4,1) * t517 - mrSges(4,3) * t508;
t540 = -t374 * t522 + t526 * t375;
t369 = m(4) * t448 - mrSges(4,2) * t506 + mrSges(4,3) * t480 + t485 * t507 - t487 * t517 + t540;
t541 = -t368 * t523 + t527 * t369;
t468 = -qJDD(2) * pkin(2) - t530 * pkin(8) + t510 * t547 - t488;
t439 = -t480 * pkin(3) - t505 * pkin(9) + t508 * t490 + t468;
t408 = -t445 * pkin(4) - t482 * qJ(5) + t484 * t471 + qJDD(5) + t439;
t399 = t408 + (t461 * t516 - t420) * qJ(6) - 0.2e1 * qJD(6) * t462 + (t462 * t516 + t419) * pkin(5);
t390 = m(7) * t399 + t419 * mrSges(7,1) - t420 * mrSges(7,3) + t461 * t451 - t462 * t454;
t364 = t527 * t368 + t523 * t369;
t385 = m(6) * t408 + t419 * mrSges(6,1) + t420 * mrSges(6,2) + t461 * t452 + t462 * t453 + t390;
t535 = m(5) * t439 - t445 * mrSges(5,1) + t446 * mrSges(5,2) - t483 * t470 + t484 * t472 + t385;
t534 = -m(4) * t468 + t480 * mrSges(4,1) - t481 * mrSges(4,2) + t507 * t486 - t508 * t487 - t535;
t389 = t420 * mrSges(7,2) + t462 * t437 - t538;
t456 = Ifges(5,4) * t484 + Ifges(5,2) * t483 + Ifges(5,6) * t516;
t457 = Ifges(5,1) * t484 + Ifges(5,4) * t483 + Ifges(5,5) * t516;
t533 = -mrSges(5,1) * t405 - mrSges(6,1) * t396 + mrSges(7,1) * t394 + mrSges(5,2) * t406 + mrSges(6,2) * t397 - mrSges(7,3) * t393 - Ifges(5,5) * t446 - Ifges(5,6) * t445 - pkin(4) * t377 + pkin(5) * t389 - qJ(6) * t544 - t484 * t456 + t483 * t457 + t555 * t462 + (qJ(6) * t437 - t554) * t461 - t558 * t420 + (qJ(6) * mrSges(7,2) + t556) * t419 + (-Ifges(5,3) + t557) * t504;
t474 = Ifges(4,4) * t508 + Ifges(4,2) * t507 + Ifges(4,6) * t517;
t475 = Ifges(4,1) * t508 + Ifges(4,4) * t507 + Ifges(4,5) * t517;
t532 = mrSges(4,1) * t447 - mrSges(4,2) * t448 + Ifges(4,5) * t481 + Ifges(4,6) * t480 + Ifges(4,3) * t506 + pkin(3) * t370 + t508 * t474 - t507 * t475 - t533;
t514 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t546;
t513 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t547;
t509 = (-mrSges(3,1) * t528 + mrSges(3,2) * t524) * qJD(1);
t501 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t524 + Ifges(3,4) * t528) * qJD(1);
t500 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t524 + Ifges(3,2) * t528) * qJD(1);
t473 = Ifges(4,5) * t508 + Ifges(4,6) * t507 + Ifges(4,3) * t517;
t455 = Ifges(5,5) * t484 + Ifges(5,6) * t483 + Ifges(5,3) * t516;
t379 = mrSges(6,2) * t408 + mrSges(7,2) * t394 - mrSges(6,3) * t396 - mrSges(7,3) * t399 - qJ(6) * t390 - t551 * t419 + t560 * t420 + t549 * t461 + t558 * t504 + t555 * t516;
t378 = -mrSges(6,1) * t408 - mrSges(7,1) * t399 + mrSges(7,2) * t393 + mrSges(6,3) * t397 - pkin(5) * t390 - t559 * t419 + t551 * t420 + t549 * t462 + t556 * t504 + t554 * t516;
t366 = mrSges(5,2) * t439 - mrSges(5,3) * t405 + Ifges(5,1) * t446 + Ifges(5,4) * t445 + Ifges(5,5) * t504 - qJ(5) * t377 - t521 * t378 + t379 * t550 + t483 * t455 - t516 * t456;
t365 = -mrSges(5,1) * t439 + mrSges(5,3) * t406 + Ifges(5,4) * t446 + Ifges(5,2) * t445 + Ifges(5,6) * t504 - pkin(4) * t385 + qJ(5) * t539 + t378 * t550 + t521 * t379 - t484 * t455 + t516 * t457;
t363 = mrSges(4,2) * t468 - mrSges(4,3) * t447 + Ifges(4,1) * t481 + Ifges(4,4) * t480 + Ifges(4,5) * t506 - pkin(9) * t370 - t365 * t522 + t366 * t526 + t473 * t507 - t474 * t517;
t362 = -mrSges(4,1) * t468 + mrSges(4,3) * t448 + Ifges(4,4) * t481 + Ifges(4,2) * t480 + Ifges(4,6) * t506 - pkin(3) * t535 + pkin(9) * t540 + t526 * t365 + t522 * t366 - t508 * t473 + t517 * t475;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t542 - mrSges(2,2) * t537 + t524 * (mrSges(3,2) * t502 - mrSges(3,3) * t488 + Ifges(3,1) * t511 + Ifges(3,4) * t512 + Ifges(3,5) * qJDD(2) - pkin(8) * t364 - qJD(2) * t500 - t523 * t362 + t527 * t363) + t528 * (-mrSges(3,1) * t502 + mrSges(3,3) * t489 + Ifges(3,4) * t511 + Ifges(3,2) * t512 + Ifges(3,6) * qJDD(2) - pkin(2) * t364 + qJD(2) * t501 - t532) + pkin(1) * (-m(3) * t502 + t512 * mrSges(3,1) - t511 * mrSges(3,2) + (-t513 * t524 + t514 * t528) * qJD(1) - t364) + pkin(7) * (t528 * (m(3) * t489 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t512 - qJD(2) * t513 + t509 * t546 + t541) - t524 * (m(3) * t488 + qJDD(2) * mrSges(3,1) - t511 * mrSges(3,3) + qJD(2) * t514 - t509 * t547 + t534)); Ifges(3,5) * t511 + Ifges(3,6) * t512 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t488 - mrSges(3,2) * t489 + t523 * t363 + t527 * t362 + pkin(2) * t534 + pkin(8) * t541 + (t524 * t500 - t528 * t501) * qJD(1); t532; -t533; t385; t389;];
tauJ  = t1;
