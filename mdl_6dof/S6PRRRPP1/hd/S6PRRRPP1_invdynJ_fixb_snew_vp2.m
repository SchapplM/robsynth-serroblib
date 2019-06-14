% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-05-05 06:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:34:00
% EndTime: 2019-05-05 06:34:04
% DurationCPUTime: 2.43s
% Computational Cost: add. (19798->266), mult. (39230->329), div. (0->0), fcn. (27598->12), ass. (0->113)
t504 = Ifges(6,1) + Ifges(7,1);
t494 = Ifges(6,4) - Ifges(7,5);
t493 = Ifges(6,5) + Ifges(7,4);
t503 = -Ifges(6,2) - Ifges(7,3);
t492 = Ifges(6,6) - Ifges(7,6);
t502 = -Ifges(7,2) - Ifges(6,3);
t457 = sin(pkin(10));
t459 = cos(pkin(10));
t448 = g(1) * t457 - g(2) * t459;
t455 = -g(3) + qJDD(1);
t458 = sin(pkin(6));
t460 = cos(pkin(6));
t501 = t448 * t460 + t455 * t458;
t449 = -g(1) * t459 - g(2) * t457;
t463 = sin(qJ(2));
t466 = cos(qJ(2));
t409 = t466 * t449 + t501 * t463;
t468 = qJD(2) ^ 2;
t403 = -pkin(2) * t468 + qJDD(2) * pkin(8) + t409;
t427 = -t448 * t458 + t455 * t460;
t462 = sin(qJ(3));
t465 = cos(qJ(3));
t394 = t465 * t403 + t462 * t427;
t445 = (-pkin(3) * t465 - pkin(9) * t462) * qJD(2);
t467 = qJD(3) ^ 2;
t483 = qJD(2) * t465;
t375 = -pkin(3) * t467 + qJDD(3) * pkin(9) + t445 * t483 + t394;
t408 = -t463 * t449 + t501 * t466;
t402 = -qJDD(2) * pkin(2) - t468 * pkin(8) - t408;
t482 = qJD(2) * qJD(3);
t478 = t465 * t482;
t446 = qJDD(2) * t462 + t478;
t479 = t462 * t482;
t447 = qJDD(2) * t465 - t479;
t378 = (-t446 - t478) * pkin(9) + (-t447 + t479) * pkin(3) + t402;
t461 = sin(qJ(4));
t464 = cos(qJ(4));
t370 = -t461 * t375 + t464 * t378;
t484 = qJD(2) * t462;
t442 = qJD(3) * t464 - t461 * t484;
t419 = qJD(4) * t442 + qJDD(3) * t461 + t446 * t464;
t439 = qJDD(4) - t447;
t443 = qJD(3) * t461 + t464 * t484;
t454 = qJD(4) - t483;
t367 = (t442 * t454 - t419) * qJ(5) + (t442 * t443 + t439) * pkin(4) + t370;
t371 = t464 * t375 + t461 * t378;
t418 = -qJD(4) * t443 + qJDD(3) * t464 - t446 * t461;
t425 = pkin(4) * t454 - qJ(5) * t443;
t438 = t442 ^ 2;
t369 = -pkin(4) * t438 + qJ(5) * t418 - t425 * t454 + t371;
t456 = sin(pkin(11));
t491 = cos(pkin(11));
t420 = -t491 * t442 + t443 * t456;
t496 = -2 * qJD(5);
t363 = t456 * t367 + t491 * t369 + t420 * t496;
t388 = -t491 * t418 + t419 * t456;
t421 = t456 * t442 + t491 * t443;
t405 = mrSges(6,1) * t454 - mrSges(6,3) * t421;
t395 = pkin(5) * t420 - qJ(6) * t421;
t453 = t454 ^ 2;
t360 = -pkin(5) * t453 + qJ(6) * t439 + 0.2e1 * qJD(6) * t454 - t395 * t420 + t363;
t406 = -mrSges(7,1) * t454 + mrSges(7,2) * t421;
t480 = m(7) * t360 + t439 * mrSges(7,3) + t454 * t406;
t396 = mrSges(7,1) * t420 - mrSges(7,3) * t421;
t485 = -mrSges(6,1) * t420 - mrSges(6,2) * t421 - t396;
t495 = -mrSges(6,3) - mrSges(7,2);
t351 = m(6) * t363 - t439 * mrSges(6,2) + t495 * t388 - t454 * t405 + t485 * t420 + t480;
t472 = t491 * t367 - t456 * t369;
t362 = t421 * t496 + t472;
t389 = t456 * t418 + t491 * t419;
t404 = -mrSges(6,2) * t454 - mrSges(6,3) * t420;
t361 = -t439 * pkin(5) - t453 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t395) * t421 - t472;
t407 = -mrSges(7,2) * t420 + mrSges(7,3) * t454;
t474 = -m(7) * t361 + t439 * mrSges(7,1) + t454 * t407;
t353 = m(6) * t362 + mrSges(6,1) * t439 + t495 * t389 + t404 * t454 + t485 * t421 + t474;
t346 = t456 * t351 + t491 * t353;
t357 = mrSges(7,2) * t389 + t396 * t421 - t474;
t412 = Ifges(5,4) * t443 + Ifges(5,2) * t442 + Ifges(5,6) * t454;
t413 = Ifges(5,1) * t443 + Ifges(5,4) * t442 + Ifges(5,5) * t454;
t486 = -t494 * t420 + t421 * t504 + t493 * t454;
t487 = t420 * t503 + t421 * t494 + t454 * t492;
t500 = -t492 * t388 + t493 * t389 + (Ifges(5,3) - t502) * t439 + mrSges(5,1) * t370 + mrSges(6,1) * t362 - mrSges(7,1) * t361 - mrSges(5,2) * t371 - mrSges(6,2) * t363 + mrSges(7,3) * t360 + Ifges(5,5) * t419 + Ifges(5,6) * t418 + pkin(4) * t346 - pkin(5) * t357 + qJ(6) * (-t388 * mrSges(7,2) - t420 * t396 + t480) + t443 * t412 - t442 * t413 + t487 * t421 + t486 * t420;
t488 = t492 * t420 - t493 * t421 + t502 * t454;
t444 = (-mrSges(4,1) * t465 + mrSges(4,2) * t462) * qJD(2);
t450 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t484;
t422 = -mrSges(5,1) * t442 + mrSges(5,2) * t443;
t424 = -mrSges(5,2) * t454 + mrSges(5,3) * t442;
t344 = m(5) * t370 + mrSges(5,1) * t439 - mrSges(5,3) * t419 - t422 * t443 + t424 * t454 + t346;
t426 = mrSges(5,1) * t454 - mrSges(5,3) * t443;
t475 = t491 * t351 - t353 * t456;
t345 = m(5) * t371 - mrSges(5,2) * t439 + mrSges(5,3) * t418 + t422 * t442 - t426 * t454 + t475;
t476 = -t344 * t461 + t464 * t345;
t341 = m(4) * t394 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t447 - qJD(3) * t450 + t444 * t483 + t476;
t393 = -t462 * t403 + t427 * t465;
t451 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t483;
t374 = -qJDD(3) * pkin(3) - pkin(9) * t467 + t445 * t484 - t393;
t372 = -pkin(4) * t418 - qJ(5) * t438 + t443 * t425 + qJDD(5) + t374;
t365 = -0.2e1 * qJD(6) * t421 + (t420 * t454 - t389) * qJ(6) + (t421 * t454 + t388) * pkin(5) + t372;
t358 = m(7) * t365 + t388 * mrSges(7,1) - t389 * mrSges(7,3) - t421 * t406 + t420 * t407;
t355 = m(6) * t372 + t388 * mrSges(6,1) + t389 * mrSges(6,2) + t420 * t404 + t421 * t405 + t358;
t470 = -m(5) * t374 + t418 * mrSges(5,1) - t419 * mrSges(5,2) + t442 * t424 - t443 * t426 - t355;
t354 = m(4) * t393 + qJDD(3) * mrSges(4,1) - t446 * mrSges(4,3) + qJD(3) * t451 - t444 * t484 + t470;
t477 = t465 * t341 - t354 * t462;
t342 = t344 * t464 + t345 * t461;
t471 = -m(4) * t402 + t447 * mrSges(4,1) - mrSges(4,2) * t446 - t450 * t484 + t451 * t483 - t342;
t433 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t462 + Ifges(4,4) * t465) * qJD(2);
t432 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t462 + Ifges(4,2) * t465) * qJD(2);
t411 = Ifges(5,5) * t443 + Ifges(5,6) * t442 + Ifges(5,3) * t454;
t348 = mrSges(6,2) * t372 + mrSges(7,2) * t361 - mrSges(6,3) * t362 - mrSges(7,3) * t365 - qJ(6) * t358 - t494 * t388 + t389 * t504 + t488 * t420 + t493 * t439 - t487 * t454;
t347 = -mrSges(6,1) * t372 - mrSges(7,1) * t365 + mrSges(7,2) * t360 + mrSges(6,3) * t363 - pkin(5) * t358 + t388 * t503 + t494 * t389 + t488 * t421 + t492 * t439 + t486 * t454;
t339 = mrSges(5,2) * t374 - mrSges(5,3) * t370 + Ifges(5,1) * t419 + Ifges(5,4) * t418 + Ifges(5,5) * t439 - qJ(5) * t346 - t456 * t347 + t491 * t348 + t442 * t411 - t454 * t412;
t338 = -mrSges(5,1) * t374 + mrSges(5,3) * t371 + Ifges(5,4) * t419 + Ifges(5,2) * t418 + Ifges(5,6) * t439 - pkin(4) * t355 + qJ(5) * t475 + t491 * t347 + t456 * t348 - t443 * t411 + t454 * t413;
t1 = [m(2) * t455 + t460 * (m(3) * t427 + t341 * t462 + t354 * t465) + (t463 * (m(3) * t409 - mrSges(3,1) * t468 - qJDD(2) * mrSges(3,2) + t477) + t466 * (m(3) * t408 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t468 + t471)) * t458; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t408 - mrSges(3,2) * t409 + t462 * (mrSges(4,2) * t402 - mrSges(4,3) * t393 + Ifges(4,1) * t446 + Ifges(4,4) * t447 + Ifges(4,5) * qJDD(3) - pkin(9) * t342 - qJD(3) * t432 - t338 * t461 + t339 * t464) + pkin(2) * t471 + pkin(8) * t477 + (-mrSges(4,1) * t402 + mrSges(4,3) * t394 + Ifges(4,4) * t446 + Ifges(4,2) * t447 + Ifges(4,6) * qJDD(3) - pkin(3) * t342 + qJD(3) * t433 - t500) * t465; Ifges(4,5) * t446 + Ifges(4,6) * t447 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t393 - mrSges(4,2) * t394 + t461 * t339 + t464 * t338 + pkin(3) * t470 + pkin(9) * t476 + (t432 * t462 - t433 * t465) * qJD(2); t500; t355; t357;];
tauJ  = t1;
