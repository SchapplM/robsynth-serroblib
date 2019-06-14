% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:24:16
% EndTime: 2019-05-06 01:24:23
% DurationCPUTime: 4.94s
% Computational Cost: add. (49656->297), mult. (119147->363), div. (0->0), fcn. (93113->10), ass. (0->127)
t558 = Ifges(6,1) + Ifges(7,1);
t552 = Ifges(6,4) + Ifges(7,4);
t551 = Ifges(6,5) + Ifges(7,5);
t557 = Ifges(6,2) + Ifges(7,2);
t550 = Ifges(6,6) + Ifges(7,6);
t556 = Ifges(6,3) + Ifges(7,3);
t519 = qJD(1) ^ 2;
t509 = sin(pkin(10));
t510 = cos(pkin(10));
t513 = sin(qJ(3));
t517 = cos(qJ(3));
t526 = -t509 * t513 + t510 * t517;
t495 = t526 * qJD(1);
t527 = t509 * t517 + t510 * t513;
t496 = t527 * qJD(1);
t512 = sin(qJ(4));
t516 = cos(qJ(4));
t480 = t495 * t516 - t496 * t512;
t485 = -t496 * qJD(3) + t526 * qJDD(1);
t540 = t495 * qJD(3);
t486 = t527 * qJDD(1) + t540;
t450 = qJD(4) * t480 + t485 * t512 + t486 * t516;
t481 = t495 * t512 + t496 * t516;
t508 = qJD(3) + qJD(4);
t511 = sin(qJ(5));
t515 = cos(qJ(5));
t467 = -t481 * t511 + t508 * t515;
t505 = qJDD(3) + qJDD(4);
t426 = qJD(5) * t467 + t450 * t515 + t505 * t511;
t468 = t481 * t515 + t508 * t511;
t451 = -mrSges(7,1) * t467 + mrSges(7,2) * t468;
t514 = sin(qJ(1));
t518 = cos(qJ(1));
t530 = -g(1) * t518 - g(2) * t514;
t497 = -pkin(1) * t519 + qJDD(1) * qJ(2) + t530;
t539 = qJD(1) * qJD(2);
t535 = -t510 * g(3) - 0.2e1 * t509 * t539;
t548 = pkin(7) * qJDD(1);
t554 = pkin(2) * t519;
t474 = (t510 * t554 - t497 - t548) * t509 + t535;
t488 = -g(3) * t509 + (t497 + 0.2e1 * t539) * t510;
t507 = t510 ^ 2;
t475 = -t507 * t554 + t510 * t548 + t488;
t453 = t517 * t474 - t513 * t475;
t421 = (-t486 + t540) * pkin(8) + (t495 * t496 + qJDD(3)) * pkin(3) + t453;
t454 = t513 * t474 + t517 * t475;
t491 = qJD(3) * pkin(3) - pkin(8) * t496;
t494 = t495 ^ 2;
t429 = -pkin(3) * t494 + pkin(8) * t485 - qJD(3) * t491 + t454;
t419 = t512 * t421 + t516 * t429;
t465 = -pkin(4) * t480 - pkin(9) * t481;
t504 = t508 ^ 2;
t413 = -pkin(4) * t504 + pkin(9) * t505 + t465 * t480 + t419;
t536 = t514 * g(1) - t518 * g(2);
t529 = qJDD(2) - t536;
t542 = -t509 ^ 2 - t507;
t484 = (-pkin(2) * t510 - pkin(1)) * qJDD(1) + (t542 * pkin(7) - qJ(2)) * t519 + t529;
t442 = -t485 * pkin(3) - t494 * pkin(8) + t496 * t491 + t484;
t449 = -qJD(4) * t481 + t485 * t516 - t486 * t512;
t416 = (-t480 * t508 - t450) * pkin(9) + (t481 * t508 - t449) * pkin(4) + t442;
t408 = -t511 * t413 + t515 * t416;
t448 = qJDD(5) - t449;
t476 = qJD(5) - t480;
t405 = -0.2e1 * qJD(6) * t468 + (t467 * t476 - t426) * qJ(6) + (t467 * t468 + t448) * pkin(5) + t408;
t455 = -mrSges(7,2) * t476 + mrSges(7,3) * t467;
t538 = m(7) * t405 + t448 * mrSges(7,1) + t476 * t455;
t402 = -t426 * mrSges(7,3) - t468 * t451 + t538;
t409 = t515 * t413 + t511 * t416;
t425 = -qJD(5) * t468 - t450 * t511 + t505 * t515;
t457 = pkin(5) * t476 - qJ(6) * t468;
t466 = t467 ^ 2;
t407 = -pkin(5) * t466 + qJ(6) * t425 + 0.2e1 * qJD(6) * t467 - t457 * t476 + t409;
t544 = t552 * t467 + t558 * t468 + t551 * t476;
t545 = -t557 * t467 - t552 * t468 - t550 * t476;
t555 = mrSges(6,1) * t408 + mrSges(7,1) * t405 - mrSges(6,2) * t409 - mrSges(7,2) * t407 + pkin(5) * t402 + t550 * t425 + t551 * t426 + t556 * t448 - t544 * t467 - t545 * t468;
t553 = -mrSges(6,2) - mrSges(7,2);
t464 = -mrSges(5,1) * t480 + mrSges(5,2) * t481;
t472 = mrSges(5,1) * t508 - mrSges(5,3) * t481;
t452 = -mrSges(6,1) * t467 + mrSges(6,2) * t468;
t456 = -mrSges(6,2) * t476 + mrSges(6,3) * t467;
t395 = m(6) * t408 + t448 * mrSges(6,1) + t476 * t456 + (-t451 - t452) * t468 + (-mrSges(6,3) - mrSges(7,3)) * t426 + t538;
t537 = m(7) * t407 + t425 * mrSges(7,3) + t467 * t451;
t458 = mrSges(7,1) * t476 - mrSges(7,3) * t468;
t543 = -mrSges(6,1) * t476 + mrSges(6,3) * t468 - t458;
t398 = m(6) * t409 + t425 * mrSges(6,3) + t553 * t448 + t467 * t452 + t543 * t476 + t537;
t532 = -t395 * t511 + t515 * t398;
t389 = m(5) * t419 - mrSges(5,2) * t505 + mrSges(5,3) * t449 + t464 * t480 - t472 * t508 + t532;
t418 = t421 * t516 - t512 * t429;
t471 = -mrSges(5,2) * t508 + mrSges(5,3) * t480;
t412 = -pkin(4) * t505 - pkin(9) * t504 + t481 * t465 - t418;
t410 = -pkin(5) * t425 - qJ(6) * t466 + t457 * t468 + qJDD(6) + t412;
t531 = -m(7) * t410 + t425 * mrSges(7,1) + t467 * t455;
t521 = -m(6) * t412 + t425 * mrSges(6,1) + t553 * t426 + t467 * t456 + t543 * t468 + t531;
t400 = m(5) * t418 + t505 * mrSges(5,1) - t450 * mrSges(5,3) - t481 * t464 + t508 * t471 + t521;
t383 = t512 * t389 + t516 * t400;
t483 = -mrSges(4,1) * t495 + mrSges(4,2) * t496;
t489 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t495;
t381 = m(4) * t453 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t486 + qJD(3) * t489 - t483 * t496 + t383;
t490 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t496;
t533 = t516 * t389 - t400 * t512;
t382 = m(4) * t454 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t485 - qJD(3) * t490 + t483 * t495 + t533;
t547 = t517 * t381 + t513 * t382;
t393 = t515 * t395 + t511 * t398;
t546 = -t550 * t467 - t551 * t468 - t556 * t476;
t534 = -t513 * t381 + t517 * t382;
t528 = -mrSges(3,1) * t510 + mrSges(3,2) * t509;
t525 = mrSges(3,3) * qJDD(1) + t519 * t528;
t524 = m(5) * t442 - t449 * mrSges(5,1) + t450 * mrSges(5,2) - t480 * t471 + t481 * t472 + t393;
t403 = t426 * mrSges(7,2) + t468 * t458 - t531;
t385 = -mrSges(6,1) * t412 + mrSges(6,3) * t409 - mrSges(7,1) * t410 + mrSges(7,3) * t407 - pkin(5) * t403 + qJ(6) * t537 + (-qJ(6) * t458 + t544) * t476 + t546 * t468 + (-mrSges(7,2) * qJ(6) + t550) * t448 + t552 * t426 + t557 * t425;
t391 = mrSges(6,2) * t412 + mrSges(7,2) * t410 - mrSges(6,3) * t408 - mrSges(7,3) * t405 - qJ(6) * t402 + t552 * t425 + t558 * t426 + t551 * t448 - t546 * t467 + t545 * t476;
t461 = Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * t508;
t462 = Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * t508;
t522 = mrSges(5,1) * t418 - mrSges(5,2) * t419 + Ifges(5,5) * t450 + Ifges(5,6) * t449 + Ifges(5,3) * t505 + pkin(4) * t521 + pkin(9) * t532 + t515 * t385 + t511 * t391 + t481 * t461 - t462 * t480;
t520 = m(4) * t484 - t485 * mrSges(4,1) + t486 * mrSges(4,2) - t495 * t489 + t496 * t490 + t524;
t493 = -qJDD(1) * pkin(1) - t519 * qJ(2) + t529;
t487 = -t509 * t497 + t535;
t479 = Ifges(4,1) * t496 + Ifges(4,4) * t495 + Ifges(4,5) * qJD(3);
t478 = Ifges(4,4) * t496 + Ifges(4,2) * t495 + Ifges(4,6) * qJD(3);
t477 = Ifges(4,5) * t496 + Ifges(4,6) * t495 + Ifges(4,3) * qJD(3);
t460 = Ifges(5,5) * t481 + Ifges(5,6) * t480 + Ifges(5,3) * t508;
t386 = t542 * t519 * mrSges(3,3) + m(3) * t493 + t528 * qJDD(1) + t520;
t377 = -mrSges(5,1) * t442 + mrSges(5,3) * t419 + Ifges(5,4) * t450 + Ifges(5,2) * t449 + Ifges(5,6) * t505 - pkin(4) * t393 - t481 * t460 + t508 * t462 - t555;
t376 = mrSges(5,2) * t442 - mrSges(5,3) * t418 + Ifges(5,1) * t450 + Ifges(5,4) * t449 + Ifges(5,5) * t505 - pkin(9) * t393 - t385 * t511 + t391 * t515 + t460 * t480 - t461 * t508;
t375 = mrSges(4,2) * t484 - mrSges(4,3) * t453 + Ifges(4,1) * t486 + Ifges(4,4) * t485 + Ifges(4,5) * qJDD(3) - pkin(8) * t383 - qJD(3) * t478 + t376 * t516 - t377 * t512 + t477 * t495;
t374 = -mrSges(4,1) * t484 + mrSges(4,3) * t454 + Ifges(4,4) * t486 + Ifges(4,2) * t485 + Ifges(4,6) * qJDD(3) - pkin(3) * t524 + pkin(8) * t533 + qJD(3) * t479 + t512 * t376 + t516 * t377 - t496 * t477;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t536 - mrSges(2,2) * t530 + t509 * (mrSges(3,2) * t493 - mrSges(3,3) * t487 + t517 * t375 - t513 * t374 - pkin(7) * t547 + (Ifges(3,1) * t509 + Ifges(3,4) * t510) * qJDD(1)) + t510 * (-mrSges(3,1) * t493 + mrSges(3,3) * t488 + t513 * t375 + t517 * t374 - pkin(2) * t520 + pkin(7) * t534 + (Ifges(3,4) * t509 + Ifges(3,2) * t510) * qJDD(1)) - pkin(1) * t386 + qJ(2) * ((m(3) * t488 + t525 * t510 + t534) * t510 + (-m(3) * t487 + t525 * t509 - t547) * t509); t386; mrSges(4,1) * t453 - mrSges(4,2) * t454 + Ifges(4,5) * t486 + Ifges(4,6) * t485 + Ifges(4,3) * qJDD(3) + pkin(3) * t383 + t478 * t496 - t479 * t495 + t522; t522; t555; t403;];
tauJ  = t1;
