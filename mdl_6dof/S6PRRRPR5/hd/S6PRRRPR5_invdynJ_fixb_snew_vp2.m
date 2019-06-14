% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:01:28
% EndTime: 2019-05-05 08:01:35
% DurationCPUTime: 6.55s
% Computational Cost: add. (77562->301), mult. (164158->397), div. (0->0), fcn. (132883->16), ass. (0->134)
t492 = sin(pkin(12));
t496 = cos(pkin(12));
t483 = g(1) * t492 - g(2) * t496;
t490 = -g(3) + qJDD(1);
t494 = sin(pkin(6));
t498 = cos(pkin(6));
t531 = t483 * t498 + t490 * t494;
t484 = -g(1) * t496 - g(2) * t492;
t502 = sin(qJ(2));
t506 = cos(qJ(2));
t456 = -t484 * t502 + t506 * t531;
t507 = qJD(2) ^ 2;
t493 = sin(pkin(7));
t527 = pkin(9) * t493;
t453 = qJDD(2) * pkin(2) + t507 * t527 + t456;
t457 = t506 * t484 + t502 * t531;
t454 = -pkin(2) * t507 + qJDD(2) * t527 + t457;
t501 = sin(qJ(3));
t505 = cos(qJ(3));
t497 = cos(pkin(7));
t522 = t497 * t501;
t470 = -t483 * t494 + t490 * t498;
t525 = t470 * t493;
t420 = t453 * t522 + t454 * t505 + t501 * t525;
t520 = qJD(2) * t493;
t475 = (-pkin(3) * t505 - pkin(10) * t501) * t520;
t489 = qJD(2) * t497 + qJD(3);
t487 = t489 ^ 2;
t488 = qJDD(2) * t497 + qJDD(3);
t517 = t505 * t520;
t414 = -pkin(3) * t487 + pkin(10) * t488 + t475 * t517 + t420;
t466 = t497 * t470;
t519 = qJD(2) * qJD(3);
t476 = (qJDD(2) * t501 + t505 * t519) * t493;
t477 = (-qJDD(2) * t505 + t501 * t519) * t493;
t417 = t477 * pkin(3) - t476 * pkin(10) + t466 + (-t453 + (pkin(3) * t501 - pkin(10) * t505) * t489 * qJD(2)) * t493;
t500 = sin(qJ(4));
t504 = cos(qJ(4));
t402 = -t500 * t414 + t417 * t504;
t518 = t501 * t520;
t468 = t489 * t504 - t500 * t518;
t445 = qJD(4) * t468 + t476 * t504 + t488 * t500;
t469 = t489 * t500 + t504 * t518;
t471 = qJDD(4) + t477;
t482 = qJD(4) - t517;
t399 = (t468 * t482 - t445) * qJ(5) + (t468 * t469 + t471) * pkin(4) + t402;
t403 = t414 * t504 + t417 * t500;
t444 = -qJD(4) * t469 - t476 * t500 + t488 * t504;
t459 = pkin(4) * t482 - qJ(5) * t469;
t467 = t468 ^ 2;
t401 = -pkin(4) * t467 + qJ(5) * t444 - t459 * t482 + t403;
t491 = sin(pkin(13));
t495 = cos(pkin(13));
t451 = t468 * t495 - t469 * t491;
t528 = 2 * qJD(5);
t396 = t399 * t491 + t401 * t495 + t451 * t528;
t452 = t468 * t491 + t469 * t495;
t432 = -pkin(5) * t451 - pkin(11) * t452;
t481 = t482 ^ 2;
t394 = -pkin(5) * t481 + pkin(11) * t471 + t432 * t451 + t396;
t419 = -t501 * t454 + t505 * (t453 * t497 + t525);
t413 = -t488 * pkin(3) - t487 * pkin(10) + t475 * t518 - t419;
t404 = -t444 * pkin(4) - t467 * qJ(5) + t459 * t469 + qJDD(5) + t413;
t425 = t444 * t495 - t445 * t491;
t426 = t444 * t491 + t445 * t495;
t397 = (-t451 * t482 - t426) * pkin(11) + (t452 * t482 - t425) * pkin(5) + t404;
t499 = sin(qJ(6));
t503 = cos(qJ(6));
t391 = -t394 * t499 + t397 * t503;
t434 = -t452 * t499 + t482 * t503;
t407 = qJD(6) * t434 + t426 * t503 + t471 * t499;
t435 = t452 * t503 + t482 * t499;
t418 = -mrSges(7,1) * t434 + mrSges(7,2) * t435;
t448 = qJD(6) - t451;
t421 = -mrSges(7,2) * t448 + mrSges(7,3) * t434;
t424 = qJDD(6) - t425;
t388 = m(7) * t391 + mrSges(7,1) * t424 - mrSges(7,3) * t407 - t418 * t435 + t421 * t448;
t392 = t394 * t503 + t397 * t499;
t406 = -qJD(6) * t435 - t426 * t499 + t471 * t503;
t422 = mrSges(7,1) * t448 - mrSges(7,3) * t435;
t389 = m(7) * t392 - mrSges(7,2) * t424 + mrSges(7,3) * t406 + t418 * t434 - t422 * t448;
t380 = -t388 * t499 + t389 * t503;
t431 = -mrSges(6,1) * t451 + mrSges(6,2) * t452;
t437 = mrSges(6,1) * t482 - mrSges(6,3) * t452;
t377 = m(6) * t396 - mrSges(6,2) * t471 + mrSges(6,3) * t425 + t431 * t451 - t437 * t482 + t380;
t513 = -t495 * t399 + t491 * t401;
t393 = -t471 * pkin(5) - t481 * pkin(11) + (t528 + t432) * t452 + t513;
t390 = -m(7) * t393 + mrSges(7,1) * t406 - mrSges(7,2) * t407 + t421 * t434 - t422 * t435;
t395 = -0.2e1 * qJD(5) * t452 - t513;
t436 = -mrSges(6,2) * t482 + mrSges(6,3) * t451;
t384 = m(6) * t395 + mrSges(6,1) * t471 - mrSges(6,3) * t426 - t431 * t452 + t436 * t482 + t390;
t372 = t377 * t491 + t384 * t495;
t408 = Ifges(7,5) * t435 + Ifges(7,6) * t434 + Ifges(7,3) * t448;
t410 = Ifges(7,1) * t435 + Ifges(7,4) * t434 + Ifges(7,5) * t448;
t381 = -mrSges(7,1) * t393 + mrSges(7,3) * t392 + Ifges(7,4) * t407 + Ifges(7,2) * t406 + Ifges(7,6) * t424 - t408 * t435 + t410 * t448;
t409 = Ifges(7,4) * t435 + Ifges(7,2) * t434 + Ifges(7,6) * t448;
t382 = mrSges(7,2) * t393 - mrSges(7,3) * t391 + Ifges(7,1) * t407 + Ifges(7,4) * t406 + Ifges(7,5) * t424 + t408 * t434 - t409 * t448;
t428 = Ifges(6,4) * t452 + Ifges(6,2) * t451 + Ifges(6,6) * t482;
t429 = Ifges(6,1) * t452 + Ifges(6,4) * t451 + Ifges(6,5) * t482;
t439 = Ifges(5,4) * t469 + Ifges(5,2) * t468 + Ifges(5,6) * t482;
t440 = Ifges(5,1) * t469 + Ifges(5,4) * t468 + Ifges(5,5) * t482;
t530 = Ifges(5,5) * t445 + Ifges(5,6) * t444 + t469 * t439 - t468 * t440 + mrSges(5,1) * t402 - mrSges(5,2) * t403 + Ifges(6,5) * t426 + Ifges(6,6) * t425 + t452 * t428 - t451 * t429 + mrSges(6,1) * t395 - mrSges(6,2) * t396 + t499 * t382 + t503 * t381 + pkin(5) * t390 + pkin(11) * t380 + pkin(4) * t372 + (Ifges(5,3) + Ifges(6,3)) * t471;
t473 = -mrSges(4,2) * t489 + mrSges(4,3) * t517;
t474 = (-mrSges(4,1) * t505 + mrSges(4,2) * t501) * t520;
t379 = t388 * t503 + t389 * t499;
t378 = m(6) * t404 - mrSges(6,1) * t425 + mrSges(6,2) * t426 - t436 * t451 + t437 * t452 + t379;
t458 = -mrSges(5,2) * t482 + mrSges(5,3) * t468;
t460 = mrSges(5,1) * t482 - mrSges(5,3) * t469;
t509 = -m(5) * t413 + mrSges(5,1) * t444 - mrSges(5,2) * t445 + t458 * t468 - t460 * t469 - t378;
t374 = m(4) * t419 + mrSges(4,1) * t488 - mrSges(4,3) * t476 + t473 * t489 - t474 * t518 + t509;
t526 = t374 * t505;
t472 = mrSges(4,1) * t489 - mrSges(4,3) * t518;
t455 = -mrSges(5,1) * t468 + mrSges(5,2) * t469;
t370 = m(5) * t402 + mrSges(5,1) * t471 - mrSges(5,3) * t445 - t455 * t469 + t458 * t482 + t372;
t514 = t377 * t495 - t384 * t491;
t371 = m(5) * t403 - mrSges(5,2) * t471 + mrSges(5,3) * t444 + t455 * t468 - t460 * t482 + t514;
t515 = -t370 * t500 + t371 * t504;
t362 = m(4) * t420 - mrSges(4,2) * t488 - mrSges(4,3) * t477 - t472 * t489 + t474 * t517 + t515;
t521 = t362 * t522 + t497 * t526;
t364 = t370 * t504 + t371 * t500;
t516 = t362 * t505 - t374 * t501;
t510 = mrSges(7,1) * t391 - mrSges(7,2) * t392 + Ifges(7,5) * t407 + Ifges(7,6) * t406 + Ifges(7,3) * t424 + t409 * t435 - t410 * t434;
t464 = Ifges(4,5) * t489 + (Ifges(4,1) * t501 + Ifges(4,4) * t505) * t520;
t463 = Ifges(4,6) * t489 + (Ifges(4,4) * t501 + Ifges(4,2) * t505) * t520;
t438 = Ifges(5,5) * t469 + Ifges(5,6) * t468 + Ifges(5,3) * t482;
t433 = -t493 * t453 + t466;
t427 = Ifges(6,5) * t452 + Ifges(6,6) * t451 + Ifges(6,3) * t482;
t366 = -mrSges(6,1) * t404 + mrSges(6,3) * t396 + Ifges(6,4) * t426 + Ifges(6,2) * t425 + Ifges(6,6) * t471 - pkin(5) * t379 - t427 * t452 + t429 * t482 - t510;
t365 = mrSges(6,2) * t404 - mrSges(6,3) * t395 + Ifges(6,1) * t426 + Ifges(6,4) * t425 + Ifges(6,5) * t471 - pkin(11) * t379 - t381 * t499 + t382 * t503 + t427 * t451 - t428 * t482;
t363 = m(4) * t433 + t477 * mrSges(4,1) + t476 * mrSges(4,2) + (t472 * t501 - t473 * t505) * t520 + t364;
t359 = mrSges(5,2) * t413 - mrSges(5,3) * t402 + Ifges(5,1) * t445 + Ifges(5,4) * t444 + Ifges(5,5) * t471 - qJ(5) * t372 + t365 * t495 - t366 * t491 + t438 * t468 - t439 * t482;
t358 = -mrSges(5,1) * t413 + mrSges(5,3) * t403 + Ifges(5,4) * t445 + Ifges(5,2) * t444 + Ifges(5,6) * t471 - pkin(4) * t378 + qJ(5) * t514 + t491 * t365 + t495 * t366 - t469 * t438 + t482 * t440;
t357 = Ifges(4,5) * t476 - Ifges(4,6) * t477 + Ifges(4,3) * t488 + mrSges(4,1) * t419 - mrSges(4,2) * t420 + t500 * t359 + t504 * t358 + pkin(3) * t509 + pkin(10) * t515 + (t463 * t501 - t464 * t505) * t520;
t1 = [m(2) * t490 + t498 * (m(3) * t470 + t497 * t363 + (t362 * t501 + t526) * t493) + (t502 * (m(3) * t457 - mrSges(3,1) * t507 - qJDD(2) * mrSges(3,2) + t516) + t506 * (m(3) * t456 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t507 - t363 * t493 + t521)) * t494; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t456 - mrSges(3,2) * t457 + t497 * t357 + pkin(2) * t521 + (t501 * (mrSges(4,2) * t433 - mrSges(4,3) * t419 + Ifges(4,1) * t476 - Ifges(4,4) * t477 + Ifges(4,5) * t488 - pkin(10) * t364 - t358 * t500 + t359 * t504 - t463 * t489) + t505 * (-mrSges(4,1) * t433 + mrSges(4,3) * t420 + Ifges(4,4) * t476 - Ifges(4,2) * t477 + Ifges(4,6) * t488 - pkin(3) * t364 + t489 * t464 - t530) - pkin(2) * t363 + pkin(9) * t516) * t493; t357; t530; t378; t510;];
tauJ  = t1;
