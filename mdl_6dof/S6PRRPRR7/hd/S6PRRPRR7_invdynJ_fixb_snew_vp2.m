% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR7
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:01:01
% EndTime: 2019-05-05 06:01:06
% DurationCPUTime: 2.25s
% Computational Cost: add. (14865->272), mult. (29542->332), div. (0->0), fcn. (18704->12), ass. (0->121)
t502 = sin(qJ(3));
t506 = cos(qJ(3));
t539 = Ifges(4,4) + Ifges(5,6);
t551 = t502 * t539 + t506 * (Ifges(4,2) + Ifges(5,3));
t550 = t502 * (Ifges(5,2) + Ifges(4,1)) + t506 * t539;
t549 = -2 * qJD(4);
t538 = (Ifges(4,5) - Ifges(5,4));
t537 = (Ifges(4,6) - Ifges(5,5));
t496 = sin(pkin(11));
t498 = cos(pkin(11));
t477 = g(1) * t496 - g(2) * t498;
t495 = -g(3) + qJDD(1);
t497 = sin(pkin(6));
t499 = cos(pkin(6));
t548 = t477 * t499 + t495 * t497;
t478 = -g(1) * t498 - g(2) * t496;
t503 = sin(qJ(2));
t507 = cos(qJ(2));
t428 = t507 * t478 + t548 * t503;
t509 = qJD(2) ^ 2;
t424 = -pkin(2) * t509 + qJDD(2) * pkin(8) + t428;
t445 = -t477 * t497 + t495 * t499;
t418 = t506 * t424 + t502 * t445;
t470 = (-pkin(3) * t506 - qJ(4) * t502) * qJD(2);
t508 = qJD(3) ^ 2;
t530 = qJD(2) * t506;
t410 = pkin(3) * t508 - qJDD(3) * qJ(4) + (qJD(3) * t549) - t470 * t530 - t418;
t427 = -t503 * t478 + t548 * t507;
t514 = -qJDD(2) * pkin(2) - t427;
t542 = pkin(8) * t509;
t423 = t514 - t542;
t529 = qJD(2) * qJD(3);
t527 = t506 * t529;
t473 = qJDD(2) * t502 + t527;
t526 = t502 * t529;
t474 = qJDD(2) * t506 - t526;
t490 = t502 * qJD(2);
t479 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t490;
t480 = -(qJD(3) * mrSges(4,2)) + mrSges(4,3) * t530;
t481 = -mrSges(5,1) * t530 - (qJD(3) * mrSges(5,3));
t512 = pkin(3) * t526 + t490 * t549 + (-t473 - t527) * qJ(4) + t514;
t412 = -pkin(3) * t474 + t512 - t542;
t482 = mrSges(5,1) * t490 + (qJD(3) * mrSges(5,2));
t421 = t502 * t424;
t522 = -qJ(4) * t508 + t470 * t490 + qJDD(4) + t421;
t541 = pkin(9) * t509;
t543 = -pkin(3) - pkin(9);
t401 = pkin(4) * t473 + t543 * qJDD(3) + (-pkin(4) * t529 - t502 * t541 - t445) * t506 + t522;
t484 = pkin(4) * t490 - qJD(3) * pkin(9);
t494 = t506 ^ 2;
t403 = -t484 * t490 + (-pkin(4) * t494 - pkin(8)) * t509 + t543 * t474 + t512;
t501 = sin(qJ(5));
t505 = cos(qJ(5));
t393 = t505 * t401 - t403 * t501;
t468 = -qJD(3) * t501 - t505 * t530;
t437 = qJD(5) * t468 + qJDD(3) * t505 - t474 * t501;
t465 = qJDD(5) + t473;
t469 = qJD(3) * t505 - t501 * t530;
t487 = t490 + qJD(5);
t390 = (t468 * t487 - t437) * pkin(10) + (t468 * t469 + t465) * pkin(5) + t393;
t394 = t501 * t401 + t505 * t403;
t436 = -qJD(5) * t469 - qJDD(3) * t501 - t474 * t505;
t444 = pkin(5) * t487 - pkin(10) * t469;
t464 = t468 ^ 2;
t391 = -pkin(5) * t464 + pkin(10) * t436 - t444 * t487 + t394;
t500 = sin(qJ(6));
t504 = cos(qJ(6));
t388 = t390 * t504 - t391 * t500;
t438 = t468 * t504 - t469 * t500;
t409 = qJD(6) * t438 + t436 * t500 + t437 * t504;
t439 = t468 * t500 + t469 * t504;
t419 = -mrSges(7,1) * t438 + mrSges(7,2) * t439;
t485 = qJD(6) + t487;
t425 = -mrSges(7,2) * t485 + mrSges(7,3) * t438;
t458 = qJDD(6) + t465;
t385 = m(7) * t388 + mrSges(7,1) * t458 - mrSges(7,3) * t409 - t419 * t439 + t425 * t485;
t389 = t390 * t500 + t391 * t504;
t408 = -qJD(6) * t439 + t436 * t504 - t437 * t500;
t426 = mrSges(7,1) * t485 - mrSges(7,3) * t439;
t386 = m(7) * t389 - mrSges(7,2) * t458 + mrSges(7,3) * t408 + t419 * t438 - t426 * t485;
t377 = t504 * t385 + t500 * t386;
t440 = -mrSges(6,1) * t468 + mrSges(6,2) * t469;
t442 = -mrSges(6,2) * t487 + mrSges(6,3) * t468;
t374 = m(6) * t393 + mrSges(6,1) * t465 - mrSges(6,3) * t437 - t440 * t469 + t442 * t487 + t377;
t443 = mrSges(6,1) * t487 - mrSges(6,3) * t469;
t524 = -t385 * t500 + t504 * t386;
t375 = m(6) * t394 - mrSges(6,2) * t465 + mrSges(6,3) * t436 + t440 * t468 - t443 * t487 + t524;
t533 = -t501 * t374 + t505 * t375;
t521 = m(5) * t412 + t474 * mrSges(5,2) - t482 * t490 + t533;
t545 = -m(4) * t423 + t474 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t473 + t480 * t530 + (-t479 * t502 - t481 * t506) * qJD(2) - t521;
t544 = (t551 * qJD(2) + (t537 * qJD(3))) * t502 - t506 * (t550 * qJD(2) + (t538 * qJD(3)));
t536 = t445 * t506;
t417 = -t421 + t536;
t471 = (mrSges(5,2) * t506 - mrSges(5,3) * t502) * qJD(2);
t472 = (-mrSges(4,1) * t506 + mrSges(4,2) * t502) * qJD(2);
t371 = t374 * t505 + t375 * t501;
t411 = -qJDD(3) * pkin(3) + t522 - t536;
t517 = -m(5) * t411 - t473 * mrSges(5,1) - t371;
t368 = m(4) * t417 - mrSges(4,3) * t473 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t480 - t481) * qJD(3) + (-t471 - t472) * t490 + t517;
t400 = pkin(4) * t474 + qJD(3) * t484 - t494 * t541 - t410;
t396 = -pkin(5) * t436 - pkin(10) * t464 + t444 * t469 + t400;
t518 = m(7) * t396 - mrSges(7,1) * t408 + t409 * mrSges(7,2) - t425 * t438 + t439 * t426;
t513 = -m(6) * t400 + mrSges(6,1) * t436 - t437 * mrSges(6,2) + t442 * t468 - t469 * t443 - t518;
t510 = -m(5) * t410 + qJDD(3) * mrSges(5,3) + qJD(3) * t482 + t471 * t530 - t513;
t381 = t510 - qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t474 + m(4) * t418 - qJD(3) * t479 + t472 * t530;
t525 = -t368 * t502 + t506 * t381;
t414 = Ifges(7,4) * t439 + Ifges(7,2) * t438 + Ifges(7,6) * t485;
t415 = Ifges(7,1) * t439 + Ifges(7,4) * t438 + Ifges(7,5) * t485;
t516 = mrSges(7,1) * t388 - mrSges(7,2) * t389 + Ifges(7,5) * t409 + Ifges(7,6) * t408 + Ifges(7,3) * t458 + t439 * t414 - t438 * t415;
t430 = Ifges(6,4) * t469 + Ifges(6,2) * t468 + Ifges(6,6) * t487;
t431 = Ifges(6,1) * t469 + Ifges(6,4) * t468 + Ifges(6,5) * t487;
t511 = mrSges(6,1) * t393 - mrSges(6,2) * t394 + Ifges(6,5) * t437 + Ifges(6,6) * t436 + Ifges(6,3) * t465 + pkin(5) * t377 + t469 * t430 - t468 * t431 + t516;
t429 = Ifges(6,5) * t469 + Ifges(6,6) * t468 + Ifges(6,3) * t487;
t413 = Ifges(7,5) * t439 + Ifges(7,6) * t438 + Ifges(7,3) * t485;
t379 = mrSges(7,2) * t396 - mrSges(7,3) * t388 + Ifges(7,1) * t409 + Ifges(7,4) * t408 + Ifges(7,5) * t458 + t413 * t438 - t414 * t485;
t378 = -mrSges(7,1) * t396 + mrSges(7,3) * t389 + Ifges(7,4) * t409 + Ifges(7,2) * t408 + Ifges(7,6) * t458 - t413 * t439 + t415 * t485;
t370 = qJDD(3) * mrSges(5,2) + qJD(3) * t481 + t471 * t490 - t517;
t369 = -mrSges(5,3) * t473 + t481 * t530 + t521;
t367 = mrSges(6,2) * t400 - mrSges(6,3) * t393 + Ifges(6,1) * t437 + Ifges(6,4) * t436 + Ifges(6,5) * t465 - pkin(10) * t377 - t378 * t500 + t379 * t504 + t429 * t468 - t430 * t487;
t366 = -mrSges(6,1) * t400 + mrSges(6,3) * t394 + Ifges(6,4) * t437 + Ifges(6,2) * t436 + Ifges(6,6) * t465 - pkin(5) * t518 + pkin(10) * t524 + t504 * t378 + t500 * t379 - t469 * t429 + t487 * t431;
t1 = [m(2) * t495 + t499 * (m(3) * t445 + t368 * t506 + t381 * t502) + (t503 * (m(3) * t428 - mrSges(3,1) * t509 - qJDD(2) * mrSges(3,2) + t525) + t507 * (m(3) * t427 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t509 + t545)) * t497; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t427 - mrSges(3,2) * t428 + t502 * (mrSges(5,1) * t411 + mrSges(4,2) * t423 - mrSges(4,3) * t417 - mrSges(5,3) * t412 + pkin(4) * t371 - qJ(4) * t369 + t511) + t506 * (-mrSges(4,1) * t423 - mrSges(5,1) * t410 + mrSges(5,2) * t412 + mrSges(4,3) * t418 - pkin(3) * t369 - pkin(4) * t513 - pkin(9) * t533 - t505 * t366 - t501 * t367) + pkin(8) * t525 + t551 * t474 + (t502 * t538 + t506 * t537) * qJDD(3) - t544 * qJD(3) + t550 * t473 + t545 * pkin(2); mrSges(4,1) * t417 - mrSges(4,2) * t418 + mrSges(5,2) * t411 - mrSges(5,3) * t410 + t505 * t367 - t501 * t366 - pkin(9) * t371 - pkin(3) * t370 + qJ(4) * t510 + (mrSges(5,1) * qJ(4) + t537) * t474 + t538 * t473 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + t544 * qJD(2); t370; t511; t516;];
tauJ  = t1;
