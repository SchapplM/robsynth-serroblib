% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:13:38
% EndTime: 2019-05-05 07:13:42
% DurationCPUTime: 3.93s
% Computational Cost: add. (42908->290), mult. (86387->375), div. (0->0), fcn. (63115->14), ass. (0->121)
t486 = sin(pkin(11));
t489 = cos(pkin(11));
t472 = g(1) * t486 - g(2) * t489;
t484 = -g(3) + qJDD(1);
t487 = sin(pkin(6));
t490 = cos(pkin(6));
t517 = t472 * t490 + t484 * t487;
t473 = -g(1) * t489 - g(2) * t486;
t494 = sin(qJ(2));
t497 = cos(qJ(2));
t441 = -t494 * t473 + t497 * t517;
t516 = cos(qJ(4));
t442 = t497 * t473 + t494 * t517;
t498 = qJD(2) ^ 2;
t437 = -pkin(2) * t498 + qJDD(2) * pkin(8) + t442;
t455 = -t472 * t487 + t484 * t490;
t493 = sin(qJ(3));
t496 = cos(qJ(3));
t416 = -t493 * t437 + t496 * t455;
t511 = qJD(2) * qJD(3);
t510 = t496 * t511;
t470 = qJDD(2) * t493 + t510;
t407 = (-t470 + t510) * pkin(9) + (t493 * t496 * t498 + qJDD(3)) * pkin(3) + t416;
t417 = t496 * t437 + t493 * t455;
t471 = qJDD(2) * t496 - t493 * t511;
t513 = qJD(2) * t493;
t477 = qJD(3) * pkin(3) - pkin(9) * t513;
t483 = t496 ^ 2;
t408 = -pkin(3) * t483 * t498 + pkin(9) * t471 - qJD(3) * t477 + t417;
t492 = sin(qJ(4));
t395 = t492 * t407 + t516 * t408;
t462 = (t492 * t496 + t493 * t516) * qJD(2);
t431 = qJD(4) * t462 + t470 * t492 - t516 * t471;
t512 = qJD(2) * t496;
t461 = t492 * t513 - t516 * t512;
t445 = mrSges(5,1) * t461 + mrSges(5,2) * t462;
t482 = qJD(3) + qJD(4);
t454 = mrSges(5,1) * t482 - mrSges(5,3) * t462;
t481 = qJDD(3) + qJDD(4);
t444 = pkin(4) * t461 - qJ(5) * t462;
t480 = t482 ^ 2;
t389 = -pkin(4) * t480 + qJ(5) * t481 - t444 * t461 + t395;
t503 = -qJDD(2) * pkin(2) - t441;
t411 = -t471 * pkin(3) + t477 * t513 + (-pkin(9) * t483 - pkin(8)) * t498 + t503;
t432 = -t461 * qJD(4) + t470 * t516 + t492 * t471;
t392 = (t461 * t482 - t432) * qJ(5) + (t462 * t482 + t431) * pkin(4) + t411;
t485 = sin(pkin(12));
t488 = cos(pkin(12));
t450 = t462 * t488 + t482 * t485;
t384 = -0.2e1 * qJD(5) * t450 - t485 * t389 + t488 * t392;
t422 = t432 * t488 + t481 * t485;
t449 = -t462 * t485 + t482 * t488;
t382 = (t449 * t461 - t422) * pkin(10) + (t449 * t450 + t431) * pkin(5) + t384;
t385 = 0.2e1 * qJD(5) * t449 + t488 * t389 + t485 * t392;
t421 = -t432 * t485 + t481 * t488;
t435 = pkin(5) * t461 - pkin(10) * t450;
t448 = t449 ^ 2;
t383 = -pkin(5) * t448 + pkin(10) * t421 - t435 * t461 + t385;
t491 = sin(qJ(6));
t495 = cos(qJ(6));
t380 = t382 * t495 - t383 * t491;
t419 = t449 * t495 - t450 * t491;
t398 = qJD(6) * t419 + t421 * t491 + t422 * t495;
t420 = t449 * t491 + t450 * t495;
t403 = -mrSges(7,1) * t419 + mrSges(7,2) * t420;
t456 = qJD(6) + t461;
t409 = -mrSges(7,2) * t456 + mrSges(7,3) * t419;
t429 = qJDD(6) + t431;
t376 = m(7) * t380 + mrSges(7,1) * t429 - mrSges(7,3) * t398 - t403 * t420 + t409 * t456;
t381 = t382 * t491 + t383 * t495;
t397 = -qJD(6) * t420 + t421 * t495 - t422 * t491;
t410 = mrSges(7,1) * t456 - mrSges(7,3) * t420;
t377 = m(7) * t381 - mrSges(7,2) * t429 + mrSges(7,3) * t397 + t403 * t419 - t410 * t456;
t368 = t495 * t376 + t491 * t377;
t423 = -mrSges(6,1) * t449 + mrSges(6,2) * t450;
t433 = -mrSges(6,2) * t461 + mrSges(6,3) * t449;
t366 = m(6) * t384 + mrSges(6,1) * t431 - mrSges(6,3) * t422 - t423 * t450 + t433 * t461 + t368;
t434 = mrSges(6,1) * t461 - mrSges(6,3) * t450;
t506 = -t376 * t491 + t495 * t377;
t367 = m(6) * t385 - mrSges(6,2) * t431 + mrSges(6,3) * t421 + t423 * t449 - t434 * t461 + t506;
t507 = -t366 * t485 + t488 * t367;
t360 = m(5) * t395 - mrSges(5,2) * t481 - mrSges(5,3) * t431 - t445 * t461 - t454 * t482 + t507;
t394 = t516 * t407 - t492 * t408;
t388 = -t481 * pkin(4) - t480 * qJ(5) + t462 * t444 + qJDD(5) - t394;
t386 = -t421 * pkin(5) - t448 * pkin(10) + t450 * t435 + t388;
t504 = m(7) * t386 - t397 * mrSges(7,1) + mrSges(7,2) * t398 - t419 * t409 + t410 * t420;
t379 = m(6) * t388 - t421 * mrSges(6,1) + mrSges(6,2) * t422 - t449 * t433 + t434 * t450 + t504;
t453 = -mrSges(5,2) * t482 - mrSges(5,3) * t461;
t372 = m(5) * t394 + mrSges(5,1) * t481 - mrSges(5,3) * t432 - t445 * t462 + t453 * t482 - t379;
t353 = t492 * t360 + t516 * t372;
t362 = t488 * t366 + t485 * t367;
t469 = (-mrSges(4,1) * t496 + mrSges(4,2) * t493) * qJD(2);
t475 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t512;
t351 = m(4) * t416 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t470 + qJD(3) * t475 - t469 * t513 + t353;
t474 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t513;
t508 = t516 * t360 - t372 * t492;
t352 = m(4) * t417 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t471 - qJD(3) * t474 + t469 * t512 + t508;
t509 = -t351 * t493 + t496 * t352;
t502 = m(5) * t411 + t431 * mrSges(5,1) + mrSges(5,2) * t432 + t461 * t453 + t454 * t462 + t362;
t399 = Ifges(7,5) * t420 + Ifges(7,6) * t419 + Ifges(7,3) * t456;
t401 = Ifges(7,1) * t420 + Ifges(7,4) * t419 + Ifges(7,5) * t456;
t369 = -mrSges(7,1) * t386 + mrSges(7,3) * t381 + Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t429 - t399 * t420 + t401 * t456;
t400 = Ifges(7,4) * t420 + Ifges(7,2) * t419 + Ifges(7,6) * t456;
t370 = mrSges(7,2) * t386 - mrSges(7,3) * t380 + Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t429 + t399 * t419 - t400 * t456;
t412 = Ifges(6,5) * t450 + Ifges(6,6) * t449 + Ifges(6,3) * t461;
t414 = Ifges(6,1) * t450 + Ifges(6,4) * t449 + Ifges(6,5) * t461;
t355 = -mrSges(6,1) * t388 + mrSges(6,3) * t385 + Ifges(6,4) * t422 + Ifges(6,2) * t421 + Ifges(6,6) * t431 - pkin(5) * t504 + pkin(10) * t506 + t495 * t369 + t491 * t370 - t450 * t412 + t461 * t414;
t413 = Ifges(6,4) * t450 + Ifges(6,2) * t449 + Ifges(6,6) * t461;
t357 = mrSges(6,2) * t388 - mrSges(6,3) * t384 + Ifges(6,1) * t422 + Ifges(6,4) * t421 + Ifges(6,5) * t431 - pkin(10) * t368 - t369 * t491 + t370 * t495 + t412 * t449 - t413 * t461;
t439 = Ifges(5,4) * t462 - Ifges(5,2) * t461 + Ifges(5,6) * t482;
t440 = Ifges(5,1) * t462 - Ifges(5,4) * t461 + Ifges(5,5) * t482;
t501 = mrSges(5,1) * t394 - mrSges(5,2) * t395 + Ifges(5,5) * t432 - Ifges(5,6) * t431 + Ifges(5,3) * t481 - pkin(4) * t379 + qJ(5) * t507 + t488 * t355 + t485 * t357 + t462 * t439 + t461 * t440;
t500 = mrSges(7,1) * t380 - mrSges(7,2) * t381 + Ifges(7,5) * t398 + Ifges(7,6) * t397 + Ifges(7,3) * t429 + t420 * t400 - t419 * t401;
t436 = -t498 * pkin(8) + t503;
t499 = -m(4) * t436 + t471 * mrSges(4,1) - mrSges(4,2) * t470 - t474 * t513 + t475 * t512 - t502;
t460 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t493 + Ifges(4,4) * t496) * qJD(2);
t459 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t493 + Ifges(4,2) * t496) * qJD(2);
t438 = Ifges(5,5) * t462 - Ifges(5,6) * t461 + Ifges(5,3) * t482;
t349 = Ifges(5,6) * t481 + t482 * t440 - t462 * t438 + t449 * t414 - t450 * t413 + Ifges(5,4) * t432 - Ifges(6,6) * t421 - Ifges(6,5) * t422 - mrSges(5,1) * t411 + mrSges(5,3) * t395 - mrSges(6,1) * t384 + mrSges(6,2) * t385 - pkin(5) * t368 - pkin(4) * t362 + (-Ifges(5,2) - Ifges(6,3)) * t431 - t500;
t348 = mrSges(5,2) * t411 - mrSges(5,3) * t394 + Ifges(5,1) * t432 - Ifges(5,4) * t431 + Ifges(5,5) * t481 - qJ(5) * t362 - t355 * t485 + t357 * t488 - t438 * t461 - t439 * t482;
t1 = [m(2) * t484 + t490 * (m(3) * t455 + t351 * t496 + t352 * t493) + (t494 * (m(3) * t442 - mrSges(3,1) * t498 - qJDD(2) * mrSges(3,2) + t509) + t497 * (m(3) * t441 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t498 + t499)) * t487; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t441 - mrSges(3,2) * t442 + t493 * (mrSges(4,2) * t436 - mrSges(4,3) * t416 + Ifges(4,1) * t470 + Ifges(4,4) * t471 + Ifges(4,5) * qJDD(3) - pkin(9) * t353 - qJD(3) * t459 + t348 * t516 - t492 * t349) + t496 * (-mrSges(4,1) * t436 + mrSges(4,3) * t417 + Ifges(4,4) * t470 + Ifges(4,2) * t471 + Ifges(4,6) * qJDD(3) - pkin(3) * t502 + pkin(9) * t508 + qJD(3) * t460 + t492 * t348 + t516 * t349) + pkin(2) * t499 + pkin(8) * t509; Ifges(4,5) * t470 + Ifges(4,6) * t471 + mrSges(4,1) * t416 - mrSges(4,2) * t417 + t501 + (t459 * t493 - t460 * t496) * qJD(2) + pkin(3) * t353 + Ifges(4,3) * qJDD(3); t501; t379; t500;];
tauJ  = t1;
