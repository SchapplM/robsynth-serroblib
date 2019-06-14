% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:46:20
% EndTime: 2019-05-04 22:46:23
% DurationCPUTime: 3.19s
% Computational Cost: add. (29717->263), mult. (68362->340), div. (0->0), fcn. (51957->14), ass. (0->121)
t482 = sin(pkin(10));
t486 = cos(pkin(10));
t467 = g(1) * t482 - g(2) * t486;
t479 = -g(3) + qJDD(1);
t483 = sin(pkin(6));
t487 = cos(pkin(6));
t522 = t467 * t487 + t479 * t483;
t494 = qJD(2) ^ 2;
t468 = -g(1) * t486 - g(2) * t482;
t490 = sin(qJ(2));
t492 = cos(qJ(2));
t440 = -t468 * t490 + t492 * t522;
t485 = cos(pkin(11));
t478 = t485 ^ 2;
t521 = 0.2e1 * t485;
t520 = cos(qJ(4));
t519 = pkin(3) * t485;
t481 = sin(pkin(11));
t518 = mrSges(4,2) * t481;
t516 = t478 * t494;
t441 = t492 * t468 + t490 * t522;
t431 = -pkin(2) * t494 + qJDD(2) * qJ(3) + t441;
t458 = -t467 * t483 + t479 * t487;
t510 = qJD(2) * qJD(3);
t513 = t458 * t485 - 0.2e1 * t481 * t510;
t413 = (-pkin(8) * qJDD(2) + t494 * t519 - t431) * t481 + t513;
t416 = t431 * t485 + t458 * t481 + t510 * t521;
t508 = qJDD(2) * t485;
t414 = -pkin(3) * t516 + pkin(8) * t508 + t416;
t489 = sin(qJ(4));
t398 = t413 * t489 + t414 * t520;
t507 = t485 * t520;
t460 = (t481 * t489 - t507) * qJD(2);
t500 = t481 * t520 + t485 * t489;
t461 = t500 * qJD(2);
t444 = mrSges(5,1) * t460 + mrSges(5,2) * t461;
t509 = qJDD(2) * t481;
t512 = qJD(4) * t461;
t447 = -qJDD(2) * t507 + t489 * t509 + t512;
t457 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t461;
t443 = pkin(4) * t460 - qJ(5) * t461;
t493 = qJD(4) ^ 2;
t396 = -pkin(4) * t493 + qJDD(4) * qJ(5) - t443 * t460 + t398;
t477 = t481 ^ 2;
t498 = qJDD(3) - t440;
t423 = (-pkin(2) - t519) * qJDD(2) + (-qJ(3) + (-t477 - t478) * pkin(8)) * t494 + t498;
t511 = t460 * qJD(4);
t448 = qJDD(2) * t500 - t511;
t401 = (-t448 + t511) * qJ(5) + (t447 + t512) * pkin(4) + t423;
t480 = sin(pkin(12));
t484 = cos(pkin(12));
t453 = qJD(4) * t480 + t461 * t484;
t391 = -0.2e1 * qJD(5) * t453 - t480 * t396 + t401 * t484;
t436 = qJDD(4) * t480 + t448 * t484;
t452 = qJD(4) * t484 - t461 * t480;
t389 = (t452 * t460 - t436) * pkin(9) + (t452 * t453 + t447) * pkin(5) + t391;
t392 = 0.2e1 * qJD(5) * t452 + t396 * t484 + t401 * t480;
t434 = pkin(5) * t460 - pkin(9) * t453;
t435 = qJDD(4) * t484 - t448 * t480;
t451 = t452 ^ 2;
t390 = -pkin(5) * t451 + pkin(9) * t435 - t434 * t460 + t392;
t488 = sin(qJ(6));
t491 = cos(qJ(6));
t387 = t389 * t491 - t390 * t488;
t424 = t452 * t491 - t453 * t488;
t404 = qJD(6) * t424 + t435 * t488 + t436 * t491;
t425 = t452 * t488 + t453 * t491;
t409 = -mrSges(7,1) * t424 + mrSges(7,2) * t425;
t459 = qJD(6) + t460;
t417 = -mrSges(7,2) * t459 + mrSges(7,3) * t424;
t446 = qJDD(6) + t447;
t384 = m(7) * t387 + mrSges(7,1) * t446 - mrSges(7,3) * t404 - t409 * t425 + t417 * t459;
t388 = t389 * t488 + t390 * t491;
t403 = -qJD(6) * t425 + t435 * t491 - t436 * t488;
t418 = mrSges(7,1) * t459 - mrSges(7,3) * t425;
t385 = m(7) * t388 - mrSges(7,2) * t446 + mrSges(7,3) * t403 + t409 * t424 - t418 * t459;
t376 = t384 * t491 + t385 * t488;
t426 = -mrSges(6,1) * t452 + mrSges(6,2) * t453;
t432 = -mrSges(6,2) * t460 + mrSges(6,3) * t452;
t374 = m(6) * t391 + mrSges(6,1) * t447 - mrSges(6,3) * t436 - t426 * t453 + t432 * t460 + t376;
t433 = mrSges(6,1) * t460 - mrSges(6,3) * t453;
t503 = -t384 * t488 + t385 * t491;
t375 = m(6) * t392 - mrSges(6,2) * t447 + mrSges(6,3) * t435 + t426 * t452 - t433 * t460 + t503;
t504 = -t374 * t480 + t375 * t484;
t368 = m(5) * t398 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t447 - qJD(4) * t457 - t444 * t460 + t504;
t397 = t413 * t520 - t414 * t489;
t395 = -qJDD(4) * pkin(4) - qJ(5) * t493 + t443 * t461 + qJDD(5) - t397;
t393 = -pkin(5) * t435 - pkin(9) * t451 + t434 * t453 + t395;
t499 = m(7) * t393 - mrSges(7,1) * t403 + mrSges(7,2) * t404 - t417 * t424 + t418 * t425;
t386 = m(6) * t395 - mrSges(6,1) * t435 + mrSges(6,2) * t436 - t432 * t452 + t433 * t453 + t499;
t456 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t460;
t380 = m(5) * t397 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t448 + qJD(4) * t456 - t444 * t461 - t386;
t514 = t368 * t489 + t380 * t520;
t370 = t374 * t484 + t375 * t480;
t415 = -t431 * t481 + t513;
t501 = mrSges(4,3) * qJDD(2) + t494 * (-mrSges(4,1) * t485 + t518);
t362 = m(4) * t415 - t481 * t501 + t514;
t505 = t368 * t520 - t489 * t380;
t363 = m(4) * t416 + t485 * t501 + t505;
t506 = -t362 * t481 + t363 * t485;
t497 = m(5) * t423 + mrSges(5,1) * t447 + t448 * mrSges(5,2) + t456 * t460 + t461 * t457 + t370;
t429 = -qJDD(2) * pkin(2) - t494 * qJ(3) + t498;
t496 = -m(4) * t429 + mrSges(4,1) * t508 - t497 + (t477 * t494 + t516) * mrSges(4,3);
t406 = Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t459;
t407 = Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t459;
t495 = mrSges(7,1) * t387 - mrSges(7,2) * t388 + Ifges(7,5) * t404 + Ifges(7,6) * t403 + Ifges(7,3) * t446 + t425 * t406 - t424 * t407;
t439 = Ifges(5,1) * t461 - Ifges(5,4) * t460 + Ifges(5,5) * qJD(4);
t438 = Ifges(5,4) * t461 - Ifges(5,2) * t460 + Ifges(5,6) * qJD(4);
t437 = Ifges(5,5) * t461 - Ifges(5,6) * t460 + Ifges(5,3) * qJD(4);
t421 = Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t460;
t420 = Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t460;
t419 = Ifges(6,5) * t453 + Ifges(6,6) * t452 + Ifges(6,3) * t460;
t405 = Ifges(7,5) * t425 + Ifges(7,6) * t424 + Ifges(7,3) * t459;
t378 = mrSges(7,2) * t393 - mrSges(7,3) * t387 + Ifges(7,1) * t404 + Ifges(7,4) * t403 + Ifges(7,5) * t446 + t405 * t424 - t406 * t459;
t377 = -mrSges(7,1) * t393 + mrSges(7,3) * t388 + Ifges(7,4) * t404 + Ifges(7,2) * t403 + Ifges(7,6) * t446 - t405 * t425 + t407 * t459;
t369 = mrSges(4,2) * t509 - t496;
t365 = mrSges(6,2) * t395 - mrSges(6,3) * t391 + Ifges(6,1) * t436 + Ifges(6,4) * t435 + Ifges(6,5) * t447 - pkin(9) * t376 - t377 * t488 + t378 * t491 + t419 * t452 - t420 * t460;
t364 = -mrSges(6,1) * t395 + mrSges(6,3) * t392 + Ifges(6,4) * t436 + Ifges(6,2) * t435 + Ifges(6,6) * t447 - pkin(5) * t499 + pkin(9) * t503 + t491 * t377 + t488 * t378 - t453 * t419 + t460 * t421;
t360 = Ifges(5,6) * qJDD(4) - t461 * t437 + t452 * t421 - t453 * t420 + qJD(4) * t439 + Ifges(5,4) * t448 - Ifges(6,6) * t435 - Ifges(6,5) * t436 - mrSges(5,1) * t423 + mrSges(5,3) * t398 - mrSges(6,1) * t391 + mrSges(6,2) * t392 - pkin(5) * t376 + (-Ifges(5,2) - Ifges(6,3)) * t447 - pkin(4) * t370 - t495;
t359 = mrSges(5,2) * t423 - mrSges(5,3) * t397 + Ifges(5,1) * t448 - Ifges(5,4) * t447 + Ifges(5,5) * qJDD(4) - qJ(5) * t370 - qJD(4) * t438 - t364 * t480 + t365 * t484 - t437 * t460;
t1 = [m(2) * t479 + t487 * (m(3) * t458 + t362 * t485 + t363 * t481) + (t490 * (m(3) * t441 - mrSges(3,1) * t494 - qJDD(2) * mrSges(3,2) + t506) + t492 * (-t494 * mrSges(3,2) + m(3) * t440 + (mrSges(3,1) - t518) * qJDD(2) + t496)) * t483; mrSges(3,1) * t440 - mrSges(3,2) * t441 + t481 * (mrSges(4,2) * t429 - mrSges(4,3) * t415 - pkin(8) * t514 + t359 * t520 - t489 * t360) + t485 * (-mrSges(4,1) * t429 + mrSges(4,3) * t416 - pkin(3) * t497 + pkin(8) * t505 + t489 * t359 + t360 * t520) - pkin(2) * t369 + qJ(3) * t506 + (Ifges(4,2) * t478 + Ifges(3,3) + (Ifges(4,1) * t481 + Ifges(4,4) * t521) * t481) * qJDD(2); t369; mrSges(5,1) * t397 - mrSges(5,2) * t398 + Ifges(5,5) * t448 - Ifges(5,6) * t447 + Ifges(5,3) * qJDD(4) - pkin(4) * t386 + qJ(5) * t504 + t484 * t364 + t480 * t365 + t461 * t438 + t460 * t439; t386; t495;];
tauJ  = t1;
