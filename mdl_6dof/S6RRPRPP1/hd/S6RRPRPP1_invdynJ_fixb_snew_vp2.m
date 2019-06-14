% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-05-06 12:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:16:47
% EndTime: 2019-05-06 12:16:54
% DurationCPUTime: 5.15s
% Computational Cost: add. (50371->325), mult. (116312->404), div. (0->0), fcn. (82625->10), ass. (0->128)
t554 = -2 * qJD(3);
t553 = Ifges(6,1) + Ifges(7,1);
t546 = Ifges(6,4) - Ifges(7,5);
t545 = Ifges(6,5) + Ifges(7,4);
t552 = -Ifges(6,2) - Ifges(7,3);
t551 = -Ifges(7,2) - Ifges(6,3);
t544 = Ifges(6,6) - Ifges(7,6);
t513 = sin(qJ(2));
t516 = cos(qJ(2));
t534 = qJD(1) * qJD(2);
t502 = qJDD(1) * t513 + t516 * t534;
t519 = qJD(1) ^ 2;
t514 = sin(qJ(1));
t517 = cos(qJ(1));
t525 = -g(1) * t517 - g(2) * t514;
t499 = -pkin(1) * t519 + qJDD(1) * pkin(7) + t525;
t542 = t513 * t499;
t548 = pkin(2) * t519;
t462 = qJDD(2) * pkin(2) - t502 * qJ(3) - t542 + (qJ(3) * t534 + t513 * t548 - g(3)) * t516;
t485 = -g(3) * t513 + t516 * t499;
t503 = qJDD(1) * t516 - t513 * t534;
t537 = qJD(1) * t513;
t504 = qJD(2) * pkin(2) - qJ(3) * t537;
t508 = t516 ^ 2;
t463 = qJ(3) * t503 - qJD(2) * t504 - t508 * t548 + t485;
t510 = sin(pkin(9));
t511 = cos(pkin(9));
t494 = (t510 * t516 + t511 * t513) * qJD(1);
t436 = t511 * t462 - t510 * t463 + t494 * t554;
t493 = (t510 * t513 - t511 * t516) * qJD(1);
t437 = t510 * t462 + t511 * t463 + t493 * t554;
t474 = pkin(3) * t493 - pkin(8) * t494;
t518 = qJD(2) ^ 2;
t418 = -pkin(3) * t518 + qJDD(2) * pkin(8) - t474 * t493 + t437;
t531 = t514 * g(1) - t517 * g(2);
t523 = -qJDD(1) * pkin(1) - t531;
t466 = -t503 * pkin(2) + qJDD(3) + t504 * t537 + (-qJ(3) * t508 - pkin(7)) * t519 + t523;
t478 = -t502 * t510 + t503 * t511;
t479 = t502 * t511 + t503 * t510;
t421 = (qJD(2) * t493 - t479) * pkin(8) + (qJD(2) * t494 - t478) * pkin(3) + t466;
t512 = sin(qJ(4));
t515 = cos(qJ(4));
t414 = -t512 * t418 + t515 * t421;
t482 = qJD(2) * t515 - t494 * t512;
t455 = qJD(4) * t482 + qJDD(2) * t512 + t479 * t515;
t477 = qJDD(4) - t478;
t483 = qJD(2) * t512 + t494 * t515;
t492 = qJD(4) + t493;
t410 = (t482 * t492 - t455) * qJ(5) + (t482 * t483 + t477) * pkin(4) + t414;
t415 = t515 * t418 + t512 * t421;
t454 = -qJD(4) * t483 + qJDD(2) * t515 - t479 * t512;
t468 = pkin(4) * t492 - qJ(5) * t483;
t481 = t482 ^ 2;
t412 = -pkin(4) * t481 + qJ(5) * t454 - t468 * t492 + t415;
t509 = sin(pkin(10));
t543 = cos(pkin(10));
t460 = -t543 * t482 + t509 * t483;
t549 = -2 * qJD(5);
t406 = t509 * t410 + t543 * t412 + t460 * t549;
t425 = -t543 * t454 + t509 * t455;
t461 = t509 * t482 + t543 * t483;
t445 = mrSges(6,1) * t492 - mrSges(6,3) * t461;
t438 = pkin(5) * t460 - qJ(6) * t461;
t491 = t492 ^ 2;
t403 = -pkin(5) * t491 + qJ(6) * t477 + 0.2e1 * qJD(6) * t492 - t438 * t460 + t406;
t446 = -mrSges(7,1) * t492 + mrSges(7,2) * t461;
t532 = m(7) * t403 + t477 * mrSges(7,3) + t492 * t446;
t439 = mrSges(7,1) * t460 - mrSges(7,3) * t461;
t538 = -mrSges(6,1) * t460 - mrSges(6,2) * t461 - t439;
t547 = -mrSges(6,3) - mrSges(7,2);
t393 = m(6) * t406 - t477 * mrSges(6,2) + t547 * t425 - t492 * t445 + t538 * t460 + t532;
t522 = t543 * t410 - t509 * t412;
t405 = t461 * t549 + t522;
t426 = t509 * t454 + t543 * t455;
t444 = -mrSges(6,2) * t492 - mrSges(6,3) * t460;
t404 = -t477 * pkin(5) - t491 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t438) * t461 - t522;
t443 = -mrSges(7,2) * t460 + mrSges(7,3) * t492;
t526 = -m(7) * t404 + t477 * mrSges(7,1) + t492 * t443;
t395 = m(6) * t405 + t477 * mrSges(6,1) + t547 * t426 + t492 * t444 + t538 * t461 + t526;
t388 = t509 * t393 + t543 * t395;
t400 = t426 * mrSges(7,2) + t461 * t439 - t526;
t448 = Ifges(5,4) * t483 + Ifges(5,2) * t482 + Ifges(5,6) * t492;
t449 = Ifges(5,1) * t483 + Ifges(5,4) * t482 + Ifges(5,5) * t492;
t539 = -t546 * t460 + t553 * t461 + t545 * t492;
t540 = t552 * t460 + t546 * t461 + t544 * t492;
t550 = -t544 * t425 + t545 * t426 + t539 * t460 + t540 * t461 + (Ifges(5,3) - t551) * t477 + mrSges(5,1) * t414 + mrSges(6,1) * t405 - mrSges(7,1) * t404 - mrSges(5,2) * t415 - mrSges(6,2) * t406 + mrSges(7,3) * t403 + Ifges(5,5) * t455 + Ifges(5,6) * t454 + pkin(4) * t388 - pkin(5) * t400 + qJ(6) * (-t425 * mrSges(7,2) - t460 * t439 + t532) + t483 * t448 - t482 * t449;
t473 = mrSges(4,1) * t493 + mrSges(4,2) * t494;
t487 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t494;
t464 = -mrSges(5,1) * t482 + mrSges(5,2) * t483;
t467 = -mrSges(5,2) * t492 + mrSges(5,3) * t482;
t386 = m(5) * t414 + mrSges(5,1) * t477 - mrSges(5,3) * t455 - t464 * t483 + t467 * t492 + t388;
t469 = mrSges(5,1) * t492 - mrSges(5,3) * t483;
t528 = t543 * t393 - t395 * t509;
t387 = m(5) * t415 - mrSges(5,2) * t477 + mrSges(5,3) * t454 + t464 * t482 - t469 * t492 + t528;
t529 = -t386 * t512 + t515 * t387;
t380 = m(4) * t437 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t478 - qJD(2) * t487 - t473 * t493 + t529;
t486 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t493;
t417 = -qJDD(2) * pkin(3) - t518 * pkin(8) + t494 * t474 - t436;
t413 = -t454 * pkin(4) - t481 * qJ(5) + t483 * t468 + qJDD(5) + t417;
t408 = -0.2e1 * qJD(6) * t461 + (t460 * t492 - t426) * qJ(6) + (t461 * t492 + t425) * pkin(5) + t413;
t401 = m(7) * t408 + t425 * mrSges(7,1) - t426 * mrSges(7,3) + t460 * t443 - t461 * t446;
t398 = m(6) * t413 + t425 * mrSges(6,1) + mrSges(6,2) * t426 + t460 * t444 + t445 * t461 + t401;
t521 = -m(5) * t417 + t454 * mrSges(5,1) - mrSges(5,2) * t455 + t482 * t467 - t469 * t483 - t398;
t397 = m(4) * t436 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t479 + qJD(2) * t486 - t473 * t494 + t521;
t377 = t510 * t380 + t511 * t397;
t382 = t515 * t386 + t512 * t387;
t541 = t544 * t460 - t545 * t461 + t551 * t492;
t536 = qJD(1) * t516;
t530 = t511 * t380 - t397 * t510;
t381 = m(4) * t466 - t478 * mrSges(4,1) + t479 * mrSges(4,2) + t493 * t486 + t494 * t487 + t382;
t506 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t536;
t505 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t537;
t501 = (-mrSges(3,1) * t516 + mrSges(3,2) * t513) * qJD(1);
t498 = -t519 * pkin(7) + t523;
t497 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t513 + Ifges(3,4) * t516) * qJD(1);
t496 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t513 + Ifges(3,2) * t516) * qJD(1);
t484 = -t516 * g(3) - t542;
t472 = Ifges(4,1) * t494 - Ifges(4,4) * t493 + Ifges(4,5) * qJD(2);
t471 = Ifges(4,4) * t494 - Ifges(4,2) * t493 + Ifges(4,6) * qJD(2);
t470 = Ifges(4,5) * t494 - Ifges(4,6) * t493 + Ifges(4,3) * qJD(2);
t447 = Ifges(5,5) * t483 + Ifges(5,6) * t482 + Ifges(5,3) * t492;
t390 = mrSges(6,2) * t413 + mrSges(7,2) * t404 - mrSges(6,3) * t405 - mrSges(7,3) * t408 - qJ(6) * t401 - t546 * t425 + t553 * t426 + t541 * t460 + t545 * t477 - t540 * t492;
t389 = -mrSges(6,1) * t413 - mrSges(7,1) * t408 + mrSges(7,2) * t403 + mrSges(6,3) * t406 - pkin(5) * t401 + t552 * t425 + t546 * t426 + t541 * t461 + t544 * t477 + t539 * t492;
t376 = mrSges(5,2) * t417 - mrSges(5,3) * t414 + Ifges(5,1) * t455 + Ifges(5,4) * t454 + Ifges(5,5) * t477 - qJ(5) * t388 - t509 * t389 + t543 * t390 + t482 * t447 - t492 * t448;
t375 = -mrSges(5,1) * t417 + mrSges(5,3) * t415 + Ifges(5,4) * t455 + Ifges(5,2) * t454 + Ifges(5,6) * t477 - pkin(4) * t398 + qJ(5) * t528 + t543 * t389 + t509 * t390 - t483 * t447 + t492 * t449;
t374 = -mrSges(4,1) * t466 + mrSges(4,3) * t437 + Ifges(4,4) * t479 + Ifges(4,2) * t478 + Ifges(4,6) * qJDD(2) - pkin(3) * t382 + qJD(2) * t472 - t494 * t470 - t550;
t373 = mrSges(4,2) * t466 - mrSges(4,3) * t436 + Ifges(4,1) * t479 + Ifges(4,4) * t478 + Ifges(4,5) * qJDD(2) - pkin(8) * t382 - qJD(2) * t471 - t375 * t512 + t376 * t515 - t470 * t493;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t531 - mrSges(2,2) * t525 + t513 * (mrSges(3,2) * t498 - mrSges(3,3) * t484 + Ifges(3,1) * t502 + Ifges(3,4) * t503 + Ifges(3,5) * qJDD(2) - qJ(3) * t377 - qJD(2) * t496 + t511 * t373 - t510 * t374) + t516 * (-mrSges(3,1) * t498 + mrSges(3,3) * t485 + Ifges(3,4) * t502 + Ifges(3,2) * t503 + Ifges(3,6) * qJDD(2) - pkin(2) * t381 + qJ(3) * t530 + qJD(2) * t497 + t510 * t373 + t511 * t374) + pkin(1) * (-m(3) * t498 + t503 * mrSges(3,1) - t502 * mrSges(3,2) + (-t505 * t513 + t506 * t516) * qJD(1) - t381) + pkin(7) * (t516 * (m(3) * t485 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t503 - qJD(2) * t505 + t501 * t536 + t530) - t513 * (m(3) * t484 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t502 + qJD(2) * t506 - t501 * t537 + t377)); Ifges(3,5) * t502 + Ifges(3,6) * t503 + mrSges(3,1) * t484 - mrSges(3,2) * t485 + Ifges(4,5) * t479 + Ifges(4,6) * t478 + t494 * t471 + t493 * t472 + mrSges(4,1) * t436 - mrSges(4,2) * t437 + t512 * t376 + t515 * t375 + pkin(3) * t521 + pkin(8) * t529 + pkin(2) * t377 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t513 * t496 - t516 * t497) * qJD(1); t381; t550; t398; t400;];
tauJ  = t1;
