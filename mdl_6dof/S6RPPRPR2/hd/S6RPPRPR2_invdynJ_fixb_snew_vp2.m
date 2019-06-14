% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-05-05 14:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:02:03
% EndTime: 2019-05-05 14:02:05
% DurationCPUTime: 1.69s
% Computational Cost: add. (10324->250), mult. (23303->304), div. (0->0), fcn. (15494->10), ass. (0->115)
t514 = Ifges(5,1) + Ifges(6,2);
t508 = Ifges(5,4) + Ifges(6,6);
t507 = Ifges(5,5) - Ifges(6,4);
t513 = -Ifges(5,2) - Ifges(6,3);
t506 = Ifges(5,6) - Ifges(6,5);
t512 = Ifges(5,3) + Ifges(6,1);
t475 = qJD(1) ^ 2;
t511 = -2 * qJD(5);
t510 = cos(qJ(4));
t467 = cos(pkin(10));
t509 = t467 * pkin(3);
t465 = sin(pkin(10));
t505 = t465 * mrSges(4,2);
t462 = t467 ^ 2;
t504 = t462 * t475;
t471 = sin(qJ(1));
t473 = cos(qJ(1));
t490 = t471 * g(1) - g(2) * t473;
t449 = qJDD(1) * pkin(1) + t490;
t487 = -g(1) * t473 - g(2) * t471;
t450 = -pkin(1) * t475 + t487;
t466 = sin(pkin(9));
t468 = cos(pkin(9));
t429 = t466 * t449 + t468 * t450;
t413 = -pkin(2) * t475 + qJDD(1) * qJ(3) + t429;
t464 = -g(3) + qJDD(2);
t494 = qJD(1) * qJD(3);
t498 = t467 * t464 - 0.2e1 * t465 * t494;
t393 = (-pkin(7) * qJDD(1) + t475 * t509 - t413) * t465 + t498;
t399 = t465 * t464 + (t413 + 0.2e1 * t494) * t467;
t492 = qJDD(1) * t467;
t396 = -pkin(3) * t504 + pkin(7) * t492 + t399;
t470 = sin(qJ(4));
t380 = t393 * t510 - t470 * t396;
t491 = t467 * t510;
t497 = qJD(1) * t465;
t442 = -qJD(1) * t491 + t470 * t497;
t484 = t465 * t510 + t467 * t470;
t443 = t484 * qJD(1);
t419 = mrSges(5,1) * t442 + mrSges(5,2) * t443;
t495 = t442 * qJD(4);
t427 = qJDD(1) * t484 - t495;
t433 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t442;
t435 = mrSges(6,1) * t442 - qJD(4) * mrSges(6,3);
t418 = pkin(4) * t442 - qJ(5) * t443;
t474 = qJD(4) ^ 2;
t377 = -qJDD(4) * pkin(4) - t474 * qJ(5) + t443 * t418 + qJDD(5) - t380;
t372 = (t442 * t443 - qJDD(4)) * pkin(8) + (t427 + t495) * pkin(5) + t377;
t493 = qJDD(1) * t465;
t496 = qJD(4) * t443;
t426 = -qJDD(1) * t491 + t470 * t493 + t496;
t437 = pkin(5) * t443 - qJD(4) * pkin(8);
t441 = t442 ^ 2;
t461 = t465 ^ 2;
t428 = t449 * t468 - t466 * t450;
t486 = qJDD(3) - t428;
t397 = (-pkin(2) - t509) * qJDD(1) + (-qJ(3) + (-t461 - t462) * pkin(7)) * t475 + t486;
t476 = pkin(4) * t496 + t443 * t511 + (-t427 + t495) * qJ(5) + t397;
t375 = -pkin(5) * t441 - t437 * t443 + (pkin(4) + pkin(8)) * t426 + t476;
t469 = sin(qJ(6));
t472 = cos(qJ(6));
t370 = t372 * t472 - t375 * t469;
t430 = -qJD(4) * t469 + t442 * t472;
t395 = qJD(6) * t430 + qJDD(4) * t472 + t426 * t469;
t431 = qJD(4) * t472 + t442 * t469;
t400 = -mrSges(7,1) * t430 + mrSges(7,2) * t431;
t440 = qJD(6) + t443;
t403 = -mrSges(7,2) * t440 + mrSges(7,3) * t430;
t425 = qJDD(6) + t427;
t367 = m(7) * t370 + mrSges(7,1) * t425 - mrSges(7,3) * t395 - t400 * t431 + t403 * t440;
t371 = t372 * t469 + t375 * t472;
t394 = -qJD(6) * t431 - qJDD(4) * t469 + t426 * t472;
t404 = mrSges(7,1) * t440 - mrSges(7,3) * t431;
t368 = m(7) * t371 - mrSges(7,2) * t425 + mrSges(7,3) * t394 + t400 * t430 - t404 * t440;
t359 = t367 * t472 + t368 * t469;
t420 = -mrSges(6,2) * t442 - mrSges(6,3) * t443;
t482 = -m(6) * t377 - t427 * mrSges(6,1) - t443 * t420 - t359;
t356 = m(5) * t380 - mrSges(5,3) * t427 - t419 * t443 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t433 - t435) * qJD(4) + t482;
t381 = t470 * t393 + t510 * t396;
t434 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t443;
t481 = -pkin(4) * t474 + qJDD(4) * qJ(5) - t418 * t442 + t381;
t376 = qJD(4) * t511 - t481;
t436 = mrSges(6,1) * t443 + qJD(4) * mrSges(6,2);
t374 = -pkin(5) * t426 - pkin(8) * t441 + ((2 * qJD(5)) + t437) * qJD(4) + t481;
t483 = -m(7) * t374 + mrSges(7,1) * t394 - t395 * mrSges(7,2) + t403 * t430 - t431 * t404;
t479 = -m(6) * t376 + qJDD(4) * mrSges(6,3) + qJD(4) * t436 - t483;
t364 = m(5) * t381 - qJDD(4) * mrSges(5,2) - qJD(4) * t434 + (-t419 - t420) * t442 + (-mrSges(5,3) - mrSges(6,1)) * t426 + t479;
t503 = t510 * t356 + t470 * t364;
t502 = -t469 * t367 + t472 * t368;
t501 = -t512 * qJD(4) + t506 * t442 - t507 * t443;
t500 = t506 * qJD(4) + t513 * t442 + t508 * t443;
t499 = t507 * qJD(4) - t508 * t442 + t514 * t443;
t398 = -t413 * t465 + t498;
t485 = mrSges(4,3) * qJDD(1) + t475 * (-t467 * mrSges(4,1) + t505);
t352 = m(4) * t398 - t465 * t485 + t503;
t488 = -t356 * t470 + t510 * t364;
t353 = m(4) * t399 + t467 * t485 + t488;
t489 = -t352 * t465 + t467 * t353;
t379 = pkin(4) * t426 + t476;
t357 = m(6) * t379 - t426 * mrSges(6,2) - t427 * mrSges(6,3) - t442 * t435 - t443 * t436 + t502;
t385 = Ifges(7,4) * t431 + Ifges(7,2) * t430 + Ifges(7,6) * t440;
t386 = Ifges(7,1) * t431 + Ifges(7,4) * t430 + Ifges(7,5) * t440;
t480 = mrSges(7,1) * t370 - mrSges(7,2) * t371 + Ifges(7,5) * t395 + Ifges(7,6) * t394 + Ifges(7,3) * t425 + t431 * t385 - t430 * t386;
t478 = m(5) * t397 + t426 * mrSges(5,1) + mrSges(5,2) * t427 + t442 * t433 + t434 * t443 + t357;
t406 = -qJDD(1) * pkin(2) - qJ(3) * t475 + t486;
t477 = -m(4) * t406 + mrSges(4,1) * t492 - t478 + (t461 * t475 + t504) * mrSges(4,3);
t448 = (Ifges(4,5) * t465 + Ifges(4,6) * t467) * qJD(1);
t384 = Ifges(7,5) * t431 + Ifges(7,6) * t430 + Ifges(7,3) * t440;
t361 = mrSges(7,2) * t374 - mrSges(7,3) * t370 + Ifges(7,1) * t395 + Ifges(7,4) * t394 + Ifges(7,5) * t425 + t384 * t430 - t385 * t440;
t360 = -mrSges(7,1) * t374 + mrSges(7,3) * t371 + Ifges(7,4) * t395 + Ifges(7,2) * t394 + Ifges(7,6) * t425 - t384 * t431 + t386 * t440;
t358 = qJDD(4) * mrSges(6,2) + qJD(4) * t435 - t482;
t354 = mrSges(4,2) * t493 - t477;
t350 = mrSges(6,1) * t377 + mrSges(5,2) * t397 - mrSges(5,3) * t380 - mrSges(6,3) * t379 + pkin(5) * t359 - qJ(5) * t357 - t500 * qJD(4) + t507 * qJDD(4) - t508 * t426 + t514 * t427 + t501 * t442 + t480;
t349 = -mrSges(5,1) * t397 - mrSges(6,1) * t376 + mrSges(6,2) * t379 + mrSges(5,3) * t381 - pkin(4) * t357 - pkin(5) * t483 - pkin(8) * t502 + t499 * qJD(4) + t506 * qJDD(4) - t472 * t360 - t469 * t361 + t513 * t426 + t508 * t427 + t501 * t443;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t490 - mrSges(2,2) * t487 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t428 - mrSges(3,2) * t429 + t465 * (t467 * qJD(1) * t448 + mrSges(4,2) * t406 - mrSges(4,3) * t398 + t510 * t350 - t470 * t349 - pkin(7) * t503 + (Ifges(4,1) * t465 + Ifges(4,4) * t467) * qJDD(1)) + t467 * (-t448 * t497 - mrSges(4,1) * t406 + mrSges(4,3) * t399 + t470 * t350 + t510 * t349 - pkin(3) * t478 + pkin(7) * t488 + (Ifges(4,4) * t465 + Ifges(4,2) * t467) * qJDD(1)) - pkin(2) * t354 + qJ(3) * t489 + pkin(1) * (t466 * (m(3) * t429 - mrSges(3,1) * t475 - qJDD(1) * mrSges(3,2) + t489) + t468 * (t477 + (mrSges(3,1) - t505) * qJDD(1) + m(3) * t428 - mrSges(3,2) * t475)); m(3) * t464 + t352 * t467 + t353 * t465; t354; mrSges(5,1) * t380 - mrSges(5,2) * t381 + mrSges(6,2) * t377 - mrSges(6,3) * t376 + t472 * t361 - t469 * t360 - pkin(8) * t359 - pkin(4) * t358 + qJ(5) * t479 + t500 * t443 + (-qJ(5) * t420 + t499) * t442 + t507 * t427 + (-qJ(5) * mrSges(6,1) - t506) * t426 + t512 * qJDD(4); t358; t480;];
tauJ  = t1;
