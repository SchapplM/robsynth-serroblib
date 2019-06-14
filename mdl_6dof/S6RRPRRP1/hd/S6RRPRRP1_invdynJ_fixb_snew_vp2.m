% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:17:26
% EndTime: 2019-05-06 17:17:35
% DurationCPUTime: 5.85s
% Computational Cost: add. (56347->327), mult. (130462->404), div. (0->0), fcn. (96271->10), ass. (0->129)
t554 = Ifges(6,1) + Ifges(7,1);
t548 = Ifges(6,4) + Ifges(7,4);
t547 = Ifges(6,5) + Ifges(7,5);
t553 = Ifges(6,2) + Ifges(7,2);
t546 = Ifges(6,6) + Ifges(7,6);
t552 = Ifges(6,3) + Ifges(7,3);
t513 = sin(pkin(10));
t514 = cos(pkin(10));
t517 = sin(qJ(2));
t521 = cos(qJ(2));
t493 = (-t513 * t517 + t514 * t521) * qJD(1);
t494 = (t513 * t521 + t514 * t517) * qJD(1);
t516 = sin(qJ(4));
t520 = cos(qJ(4));
t478 = t493 * t520 - t494 * t516;
t537 = qJD(1) * qJD(2);
t502 = qJDD(1) * t517 + t521 * t537;
t503 = qJDD(1) * t521 - t517 * t537;
t482 = -t502 * t513 + t503 * t514;
t483 = t502 * t514 + t503 * t513;
t449 = qJD(4) * t478 + t482 * t516 + t483 * t520;
t479 = t493 * t516 + t494 * t520;
t511 = qJD(2) + qJD(4);
t515 = sin(qJ(5));
t519 = cos(qJ(5));
t470 = -t479 * t515 + t511 * t519;
t510 = qJDD(2) + qJDD(4);
t426 = qJD(5) * t470 + t449 * t519 + t510 * t515;
t471 = t479 * t519 + t511 * t515;
t450 = -mrSges(7,1) * t470 + mrSges(7,2) * t471;
t523 = qJD(1) ^ 2;
t518 = sin(qJ(1));
t522 = cos(qJ(1));
t529 = -g(1) * t522 - g(2) * t518;
t499 = -pkin(1) * t523 + qJDD(1) * pkin(7) + t529;
t544 = t499 * t517;
t550 = pkin(2) * t523;
t466 = qJDD(2) * pkin(2) - qJ(3) * t502 - t544 + (qJ(3) * t537 + t517 * t550 - g(3)) * t521;
t485 = -g(3) * t517 + t521 * t499;
t539 = qJD(1) * t517;
t504 = qJD(2) * pkin(2) - qJ(3) * t539;
t512 = t521 ^ 2;
t467 = qJ(3) * t503 - qJD(2) * t504 - t512 * t550 + t485;
t439 = -0.2e1 * qJD(3) * t494 + t514 * t466 - t467 * t513;
t418 = (qJD(2) * t493 - t483) * pkin(8) + (t493 * t494 + qJDD(2)) * pkin(3) + t439;
t440 = 0.2e1 * qJD(3) * t493 + t513 * t466 + t514 * t467;
t488 = qJD(2) * pkin(3) - pkin(8) * t494;
t492 = t493 ^ 2;
t421 = -pkin(3) * t492 + pkin(8) * t482 - qJD(2) * t488 + t440;
t416 = t516 * t418 + t520 * t421;
t462 = -pkin(4) * t478 - pkin(9) * t479;
t509 = t511 ^ 2;
t410 = -pkin(4) * t509 + pkin(9) * t510 + t462 * t478 + t416;
t534 = g(1) * t518 - t522 * g(2);
t528 = -qJDD(1) * pkin(1) - t534;
t469 = -pkin(2) * t503 + qJDD(3) + t504 * t539 + (-qJ(3) * t512 - pkin(7)) * t523 + t528;
t437 = -pkin(3) * t482 - pkin(8) * t492 + t494 * t488 + t469;
t448 = -qJD(4) * t479 + t482 * t520 - t483 * t516;
t413 = (-t478 * t511 - t449) * pkin(9) + (t479 * t511 - t448) * pkin(4) + t437;
t405 = -t410 * t515 + t519 * t413;
t447 = qJDD(5) - t448;
t474 = qJD(5) - t478;
t402 = -0.2e1 * qJD(6) * t471 + (t470 * t474 - t426) * qJ(6) + (t470 * t471 + t447) * pkin(5) + t405;
t452 = -mrSges(7,2) * t474 + mrSges(7,3) * t470;
t536 = m(7) * t402 + t447 * mrSges(7,1) + t474 * t452;
t399 = -mrSges(7,3) * t426 - t450 * t471 + t536;
t406 = t519 * t410 + t515 * t413;
t425 = -qJD(5) * t471 - t449 * t515 + t510 * t519;
t454 = pkin(5) * t474 - qJ(6) * t471;
t468 = t470 ^ 2;
t404 = -pkin(5) * t468 + qJ(6) * t425 + 0.2e1 * qJD(6) * t470 - t454 * t474 + t406;
t541 = t548 * t470 + t471 * t554 + t547 * t474;
t542 = -t470 * t553 - t471 * t548 - t474 * t546;
t551 = mrSges(6,1) * t405 + mrSges(7,1) * t402 - mrSges(6,2) * t406 - mrSges(7,2) * t404 + pkin(5) * t399 + t425 * t546 + t426 * t547 + t447 * t552 - t470 * t541 - t471 * t542;
t549 = -mrSges(6,2) - mrSges(7,2);
t461 = -mrSges(5,1) * t478 + mrSges(5,2) * t479;
t473 = mrSges(5,1) * t511 - mrSges(5,3) * t479;
t451 = -mrSges(6,1) * t470 + mrSges(6,2) * t471;
t453 = -mrSges(6,2) * t474 + mrSges(6,3) * t470;
t392 = m(6) * t405 + mrSges(6,1) * t447 + t453 * t474 + (-t450 - t451) * t471 + (-mrSges(6,3) - mrSges(7,3)) * t426 + t536;
t535 = m(7) * t404 + t425 * mrSges(7,3) + t470 * t450;
t455 = mrSges(7,1) * t474 - mrSges(7,3) * t471;
t540 = -mrSges(6,1) * t474 + mrSges(6,3) * t471 - t455;
t395 = m(6) * t406 + mrSges(6,3) * t425 + t447 * t549 + t451 * t470 + t474 * t540 + t535;
t531 = -t392 * t515 + t519 * t395;
t385 = m(5) * t416 - mrSges(5,2) * t510 + mrSges(5,3) * t448 + t461 * t478 - t473 * t511 + t531;
t415 = t418 * t520 - t516 * t421;
t472 = -mrSges(5,2) * t511 + mrSges(5,3) * t478;
t409 = -pkin(4) * t510 - pkin(9) * t509 + t479 * t462 - t415;
t407 = -pkin(5) * t425 - qJ(6) * t468 + t454 * t471 + qJDD(6) + t409;
t530 = -m(7) * t407 + t425 * mrSges(7,1) + t470 * t452;
t524 = -m(6) * t409 + t425 * mrSges(6,1) + t426 * t549 + t470 * t453 + t471 * t540 + t530;
t397 = m(5) * t415 + mrSges(5,1) * t510 - mrSges(5,3) * t449 - t461 * t479 + t472 * t511 + t524;
t380 = t516 * t385 + t520 * t397;
t481 = -mrSges(4,1) * t493 + mrSges(4,2) * t494;
t486 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t493;
t378 = m(4) * t439 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t483 + qJD(2) * t486 - t481 * t494 + t380;
t487 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t494;
t532 = t520 * t385 - t397 * t516;
t379 = m(4) * t440 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t482 - qJD(2) * t487 + t481 * t493 + t532;
t373 = t514 * t378 + t513 * t379;
t390 = t519 * t392 + t515 * t395;
t543 = -t470 * t546 - t471 * t547 - t474 * t552;
t538 = qJD(1) * t521;
t533 = -t378 * t513 + t514 * t379;
t527 = -m(5) * t437 + mrSges(5,1) * t448 - t449 * mrSges(5,2) + t472 * t478 - t479 * t473 - t390;
t400 = mrSges(7,2) * t426 + t455 * t471 - t530;
t382 = -mrSges(6,1) * t409 + mrSges(6,3) * t406 - mrSges(7,1) * t407 + mrSges(7,3) * t404 - pkin(5) * t400 + qJ(6) * t535 + (-qJ(6) * t455 + t541) * t474 + t543 * t471 + (-mrSges(7,2) * qJ(6) + t546) * t447 + t548 * t426 + t553 * t425;
t388 = mrSges(6,2) * t409 + mrSges(7,2) * t407 - mrSges(6,3) * t405 - mrSges(7,3) * t402 - qJ(6) * t399 + t548 * t425 + t426 * t554 + t547 * t447 - t543 * t470 + t542 * t474;
t458 = Ifges(5,4) * t479 + Ifges(5,2) * t478 + Ifges(5,6) * t511;
t459 = Ifges(5,1) * t479 + Ifges(5,4) * t478 + Ifges(5,5) * t511;
t525 = mrSges(5,1) * t415 - mrSges(5,2) * t416 + Ifges(5,5) * t449 + Ifges(5,6) * t448 + Ifges(5,3) * t510 + pkin(4) * t524 + pkin(9) * t531 + t519 * t382 + t515 * t388 + t479 * t458 - t478 * t459;
t386 = m(4) * t469 - mrSges(4,1) * t482 + mrSges(4,2) * t483 - t486 * t493 + t487 * t494 - t527;
t506 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t538;
t505 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t539;
t501 = (-mrSges(3,1) * t521 + mrSges(3,2) * t517) * qJD(1);
t498 = -pkin(7) * t523 + t528;
t497 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t517 + Ifges(3,4) * t521) * qJD(1);
t496 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t517 + Ifges(3,2) * t521) * qJD(1);
t484 = -g(3) * t521 - t544;
t477 = Ifges(4,1) * t494 + Ifges(4,4) * t493 + Ifges(4,5) * qJD(2);
t476 = Ifges(4,4) * t494 + Ifges(4,2) * t493 + Ifges(4,6) * qJD(2);
t475 = Ifges(4,5) * t494 + Ifges(4,6) * t493 + Ifges(4,3) * qJD(2);
t457 = Ifges(5,5) * t479 + Ifges(5,6) * t478 + Ifges(5,3) * t511;
t374 = -mrSges(5,1) * t437 + mrSges(5,3) * t416 + Ifges(5,4) * t449 + Ifges(5,2) * t448 + Ifges(5,6) * t510 - pkin(4) * t390 - t479 * t457 + t511 * t459 - t551;
t372 = mrSges(5,2) * t437 - mrSges(5,3) * t415 + Ifges(5,1) * t449 + Ifges(5,4) * t448 + Ifges(5,5) * t510 - pkin(9) * t390 - t382 * t515 + t388 * t519 + t457 * t478 - t458 * t511;
t371 = mrSges(4,2) * t469 - mrSges(4,3) * t439 + Ifges(4,1) * t483 + Ifges(4,4) * t482 + Ifges(4,5) * qJDD(2) - pkin(8) * t380 - qJD(2) * t476 + t372 * t520 - t374 * t516 + t475 * t493;
t370 = -mrSges(4,1) * t469 + mrSges(4,3) * t440 + Ifges(4,4) * t483 + Ifges(4,2) * t482 + Ifges(4,6) * qJDD(2) + pkin(3) * t527 + pkin(8) * t532 + qJD(2) * t477 + t516 * t372 + t520 * t374 - t494 * t475;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t534 - mrSges(2,2) * t529 + t517 * (mrSges(3,2) * t498 - mrSges(3,3) * t484 + Ifges(3,1) * t502 + Ifges(3,4) * t503 + Ifges(3,5) * qJDD(2) - qJ(3) * t373 - qJD(2) * t496 - t513 * t370 + t514 * t371) + t521 * (-mrSges(3,1) * t498 + mrSges(3,3) * t485 + Ifges(3,4) * t502 + Ifges(3,2) * t503 + Ifges(3,6) * qJDD(2) - pkin(2) * t386 + qJ(3) * t533 + qJD(2) * t497 + t514 * t370 + t513 * t371) + pkin(1) * ((-t505 * t517 + t506 * t521) * qJD(1) - m(3) * t498 + mrSges(3,1) * t503 - mrSges(3,2) * t502 - t386) + pkin(7) * (t521 * (m(3) * t485 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t503 - qJD(2) * t505 + t501 * t538 + t533) - t517 * (m(3) * t484 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t502 + qJD(2) * t506 - t501 * t539 + t373)); (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + Ifges(3,5) * t502 + Ifges(3,6) * t503 - t493 * t477 + t494 * t476 + (t517 * t496 - t521 * t497) * qJD(1) + Ifges(4,6) * t482 + Ifges(4,5) * t483 + mrSges(3,1) * t484 - mrSges(3,2) * t485 + mrSges(4,1) * t439 - mrSges(4,2) * t440 + t525 + pkin(3) * t380 + pkin(2) * t373; t386; t525; t551; t400;];
tauJ  = t1;
