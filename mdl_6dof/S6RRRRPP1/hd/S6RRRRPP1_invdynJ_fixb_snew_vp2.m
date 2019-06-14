% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP1
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
% Datum: 2019-05-07 18:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:54:44
% EndTime: 2019-05-07 17:54:53
% DurationCPUTime: 6.25s
% Computational Cost: add. (69394->326), mult. (140184->404), div. (0->0), fcn. (101175->10), ass. (0->130)
t555 = Ifges(6,1) + Ifges(7,1);
t548 = Ifges(6,4) - Ifges(7,5);
t547 = Ifges(6,5) + Ifges(7,4);
t554 = -Ifges(6,2) - Ifges(7,3);
t553 = -Ifges(7,2) - Ifges(6,3);
t546 = Ifges(6,6) - Ifges(7,6);
t514 = sin(qJ(3));
t515 = sin(qJ(2));
t518 = cos(qJ(3));
t519 = cos(qJ(2));
t495 = (t514 * t515 - t518 * t519) * qJD(1);
t496 = (t514 * t519 + t515 * t518) * qJD(1);
t537 = qJD(1) * qJD(2);
t501 = qJDD(1) * t515 + t519 * t537;
t502 = qJDD(1) * t519 - t515 * t537;
t470 = -qJD(3) * t496 - t501 * t514 + t502 * t518;
t471 = -qJD(3) * t495 + t501 * t518 + t502 * t514;
t539 = qJD(1) * t515;
t505 = qJD(2) * pkin(2) - pkin(8) * t539;
t511 = t519 ^ 2;
t521 = qJD(1) ^ 2;
t516 = sin(qJ(1));
t520 = cos(qJ(1));
t534 = t516 * g(1) - t520 * g(2);
t527 = -qJDD(1) * pkin(1) - t534;
t472 = -t502 * pkin(2) + t505 * t539 + (-pkin(8) * t511 - pkin(7)) * t521 + t527;
t510 = qJD(2) + qJD(3);
t415 = (t495 * t510 - t471) * pkin(9) + (t496 * t510 - t470) * pkin(3) + t472;
t529 = -g(1) * t520 - g(2) * t516;
t498 = -pkin(1) * t521 + qJDD(1) * pkin(7) + t529;
t544 = t515 * t498;
t550 = pkin(2) * t521;
t462 = qJDD(2) * pkin(2) - t501 * pkin(8) - t544 + (pkin(8) * t537 + t515 * t550 - g(3)) * t519;
t486 = -g(3) * t515 + t519 * t498;
t463 = pkin(8) * t502 - qJD(2) * t505 - t511 * t550 + t486;
t438 = t514 * t462 + t518 * t463;
t480 = pkin(3) * t495 - pkin(9) * t496;
t508 = t510 ^ 2;
t509 = qJDD(2) + qJDD(3);
t423 = -pkin(3) * t508 + pkin(9) * t509 - t480 * t495 + t438;
t513 = sin(qJ(4));
t517 = cos(qJ(4));
t410 = t517 * t415 - t513 * t423;
t483 = -t496 * t513 + t510 * t517;
t444 = qJD(4) * t483 + t471 * t517 + t509 * t513;
t469 = qJDD(4) - t470;
t484 = t496 * t517 + t510 * t513;
t491 = qJD(4) + t495;
t407 = (t483 * t491 - t444) * qJ(5) + (t483 * t484 + t469) * pkin(4) + t410;
t411 = t513 * t415 + t517 * t423;
t443 = -qJD(4) * t484 - t471 * t513 + t509 * t517;
t474 = pkin(4) * t491 - qJ(5) * t484;
t482 = t483 ^ 2;
t409 = -pkin(4) * t482 + qJ(5) * t443 - t474 * t491 + t411;
t512 = sin(pkin(10));
t545 = cos(pkin(10));
t455 = -t545 * t483 + t484 * t512;
t551 = -2 * qJD(5);
t403 = t512 * t407 + t545 * t409 + t455 * t551;
t419 = -t545 * t443 + t444 * t512;
t456 = t512 * t483 + t545 * t484;
t447 = mrSges(6,1) * t491 - mrSges(6,3) * t456;
t433 = pkin(5) * t455 - qJ(6) * t456;
t490 = t491 ^ 2;
t400 = -pkin(5) * t490 + qJ(6) * t469 + 0.2e1 * qJD(6) * t491 - t433 * t455 + t403;
t448 = -mrSges(7,1) * t491 + mrSges(7,2) * t456;
t535 = m(7) * t400 + t469 * mrSges(7,3) + t491 * t448;
t434 = mrSges(7,1) * t455 - mrSges(7,3) * t456;
t540 = -mrSges(6,1) * t455 - mrSges(6,2) * t456 - t434;
t549 = -mrSges(6,3) - mrSges(7,2);
t389 = m(6) * t403 - t469 * mrSges(6,2) + t549 * t419 - t491 * t447 + t540 * t455 + t535;
t526 = t545 * t407 - t512 * t409;
t402 = t456 * t551 + t526;
t420 = t512 * t443 + t545 * t444;
t446 = -mrSges(6,2) * t491 - mrSges(6,3) * t455;
t401 = -t469 * pkin(5) - t490 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t433) * t456 - t526;
t445 = -mrSges(7,2) * t455 + mrSges(7,3) * t491;
t530 = -m(7) * t401 + t469 * mrSges(7,1) + t491 * t445;
t391 = m(6) * t402 + t469 * mrSges(6,1) + t549 * t420 + t491 * t446 + t540 * t456 + t530;
t384 = t512 * t389 + t545 * t391;
t397 = t420 * mrSges(7,2) + t456 * t434 - t530;
t450 = Ifges(5,4) * t484 + Ifges(5,2) * t483 + Ifges(5,6) * t491;
t451 = Ifges(5,1) * t484 + Ifges(5,4) * t483 + Ifges(5,5) * t491;
t541 = -t548 * t455 + t555 * t456 + t547 * t491;
t542 = t554 * t455 + t548 * t456 + t546 * t491;
t552 = -t546 * t419 + t547 * t420 + t541 * t455 + t542 * t456 + (Ifges(5,3) - t553) * t469 + mrSges(5,1) * t410 + mrSges(6,1) * t402 - mrSges(7,1) * t401 - mrSges(5,2) * t411 - mrSges(6,2) * t403 + mrSges(7,3) * t400 + Ifges(5,5) * t444 + Ifges(5,6) * t443 + pkin(4) * t384 - pkin(5) * t397 + qJ(6) * (-t419 * mrSges(7,2) - t455 * t434 + t535) + t484 * t450 - t483 * t451;
t479 = mrSges(4,1) * t495 + mrSges(4,2) * t496;
t488 = mrSges(4,1) * t510 - mrSges(4,3) * t496;
t460 = -mrSges(5,1) * t483 + mrSges(5,2) * t484;
t473 = -mrSges(5,2) * t491 + mrSges(5,3) * t483;
t382 = m(5) * t410 + mrSges(5,1) * t469 - mrSges(5,3) * t444 - t460 * t484 + t473 * t491 + t384;
t475 = mrSges(5,1) * t491 - mrSges(5,3) * t484;
t531 = t545 * t389 - t391 * t512;
t383 = m(5) * t411 - mrSges(5,2) * t469 + mrSges(5,3) * t443 + t460 * t483 - t475 * t491 + t531;
t532 = -t382 * t513 + t517 * t383;
t376 = m(4) * t438 - mrSges(4,2) * t509 + mrSges(4,3) * t470 - t479 * t495 - t488 * t510 + t532;
t437 = t518 * t462 - t514 * t463;
t487 = -mrSges(4,2) * t510 - mrSges(4,3) * t495;
t422 = -t509 * pkin(3) - t508 * pkin(9) + t496 * t480 - t437;
t412 = -t443 * pkin(4) - t482 * qJ(5) + t484 * t474 + qJDD(5) + t422;
t405 = -0.2e1 * qJD(6) * t456 + (t455 * t491 - t420) * qJ(6) + (t456 * t491 + t419) * pkin(5) + t412;
t398 = m(7) * t405 + t419 * mrSges(7,1) - t420 * mrSges(7,3) + t455 * t445 - t456 * t448;
t395 = m(6) * t412 + t419 * mrSges(6,1) + mrSges(6,2) * t420 + t455 * t446 + t447 * t456 + t398;
t523 = -m(5) * t422 + t443 * mrSges(5,1) - mrSges(5,2) * t444 + t483 * t473 - t475 * t484 - t395;
t393 = m(4) * t437 + mrSges(4,1) * t509 - mrSges(4,3) * t471 - t479 * t496 + t487 * t510 + t523;
t373 = t514 * t376 + t518 * t393;
t378 = t517 * t382 + t513 * t383;
t543 = t546 * t455 - t547 * t456 + t553 * t491;
t538 = qJD(1) * t519;
t533 = t518 * t376 - t393 * t514;
t385 = -mrSges(6,1) * t412 - mrSges(7,1) * t405 + mrSges(7,2) * t400 + mrSges(6,3) * t403 - pkin(5) * t398 + t554 * t419 + t548 * t420 + t543 * t456 + t546 * t469 + t541 * t491;
t386 = mrSges(6,2) * t412 + mrSges(7,2) * t401 - mrSges(6,3) * t402 - mrSges(7,3) * t405 - qJ(6) * t398 - t548 * t419 + t555 * t420 + t543 * t455 + t547 * t469 - t542 * t491;
t449 = Ifges(5,5) * t484 + Ifges(5,6) * t483 + Ifges(5,3) * t491;
t370 = -mrSges(5,1) * t422 + mrSges(5,3) * t411 + Ifges(5,4) * t444 + Ifges(5,2) * t443 + Ifges(5,6) * t469 - pkin(4) * t395 + qJ(5) * t531 + t545 * t385 + t512 * t386 - t484 * t449 + t491 * t451;
t372 = mrSges(5,2) * t422 - mrSges(5,3) * t410 + Ifges(5,1) * t444 + Ifges(5,4) * t443 + Ifges(5,5) * t469 - qJ(5) * t384 - t512 * t385 + t545 * t386 + t483 * t449 - t491 * t450;
t477 = Ifges(4,4) * t496 - Ifges(4,2) * t495 + Ifges(4,6) * t510;
t478 = Ifges(4,1) * t496 - Ifges(4,4) * t495 + Ifges(4,5) * t510;
t525 = mrSges(4,1) * t437 - mrSges(4,2) * t438 + Ifges(4,5) * t471 + Ifges(4,6) * t470 + Ifges(4,3) * t509 + pkin(3) * t523 + pkin(9) * t532 + t517 * t370 + t513 * t372 + t496 * t477 + t495 * t478;
t524 = m(4) * t472 - t470 * mrSges(4,1) + t471 * mrSges(4,2) + t495 * t487 + t496 * t488 + t378;
t504 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t538;
t503 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t539;
t500 = (-mrSges(3,1) * t519 + mrSges(3,2) * t515) * qJD(1);
t497 = -t521 * pkin(7) + t527;
t494 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t515 + Ifges(3,4) * t519) * qJD(1);
t493 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t515 + Ifges(3,2) * t519) * qJD(1);
t485 = -t519 * g(3) - t544;
t476 = Ifges(4,5) * t496 - Ifges(4,6) * t495 + Ifges(4,3) * t510;
t368 = -mrSges(4,1) * t472 + mrSges(4,3) * t438 + Ifges(4,4) * t471 + Ifges(4,2) * t470 + Ifges(4,6) * t509 - pkin(3) * t378 - t496 * t476 + t510 * t478 - t552;
t367 = mrSges(4,2) * t472 - mrSges(4,3) * t437 + Ifges(4,1) * t471 + Ifges(4,4) * t470 + Ifges(4,5) * t509 - pkin(9) * t378 - t370 * t513 + t372 * t517 - t476 * t495 - t477 * t510;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t534 - mrSges(2,2) * t529 + t515 * (mrSges(3,2) * t497 - mrSges(3,3) * t485 + Ifges(3,1) * t501 + Ifges(3,4) * t502 + Ifges(3,5) * qJDD(2) - pkin(8) * t373 - qJD(2) * t493 + t518 * t367 - t514 * t368) + t519 * (-mrSges(3,1) * t497 + mrSges(3,3) * t486 + Ifges(3,4) * t501 + Ifges(3,2) * t502 + Ifges(3,6) * qJDD(2) - pkin(2) * t524 + pkin(8) * t533 + qJD(2) * t494 + t514 * t367 + t518 * t368) + pkin(1) * (-m(3) * t497 + t502 * mrSges(3,1) - t501 * mrSges(3,2) + (-t503 * t515 + t504 * t519) * qJD(1) - t524) + pkin(7) * (t519 * (m(3) * t486 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t502 - qJD(2) * t503 + t500 * t538 + t533) - t515 * (m(3) * t485 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t501 + qJD(2) * t504 - t500 * t539 + t373)); t525 + Ifges(3,3) * qJDD(2) + (t515 * t493 - t519 * t494) * qJD(1) + Ifges(3,6) * t502 + Ifges(3,5) * t501 + mrSges(3,1) * t485 - mrSges(3,2) * t486 + pkin(2) * t373; t525; t552; t395; t397;];
tauJ  = t1;
