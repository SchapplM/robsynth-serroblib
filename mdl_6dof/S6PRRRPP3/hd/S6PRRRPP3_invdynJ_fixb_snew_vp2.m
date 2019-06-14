% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:50:04
% EndTime: 2019-05-05 06:50:07
% DurationCPUTime: 1.45s
% Computational Cost: add. (7517->245), mult. (14171->284), div. (0->0), fcn. (9177->10), ass. (0->105)
t460 = sin(pkin(10));
t462 = cos(pkin(10));
t449 = g(1) * t460 - g(2) * t462;
t459 = -g(3) + qJDD(1);
t461 = sin(pkin(6));
t463 = cos(pkin(6));
t510 = t449 * t463 + t459 * t461;
t509 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t492 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t491 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t508 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t490 = -Ifges(5,6) + Ifges(6,5) + Ifges(7,4);
t507 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t450 = -g(1) * t462 - g(2) * t460;
t466 = sin(qJ(2));
t468 = cos(qJ(2));
t390 = t468 * t450 + t510 * t466;
t470 = qJD(2) ^ 2;
t386 = -pkin(2) * t470 + qJDD(2) * pkin(8) + t390;
t427 = -t449 * t461 + t459 * t463;
t465 = sin(qJ(3));
t467 = cos(qJ(3));
t382 = t467 * t386 + t465 * t427;
t446 = (-pkin(3) * t467 - pkin(9) * t465) * qJD(2);
t469 = qJD(3) ^ 2;
t494 = qJD(2) * t467;
t378 = -pkin(3) * t469 + qJDD(3) * pkin(9) + t446 * t494 + t382;
t389 = -t466 * t450 + t510 * t468;
t385 = -qJDD(2) * pkin(2) - t470 * pkin(8) - t389;
t493 = qJD(2) * qJD(3);
t483 = t467 * t493;
t447 = qJDD(2) * t465 + t483;
t484 = t465 * t493;
t448 = qJDD(2) * t467 - t484;
t380 = (-t447 - t483) * pkin(9) + (-t448 + t484) * pkin(3) + t385;
t464 = sin(qJ(4));
t502 = cos(qJ(4));
t373 = -t464 * t378 + t502 * t380;
t495 = qJD(2) * t465;
t443 = -t502 * qJD(3) + t464 * t495;
t444 = t464 * qJD(3) + t502 * t495;
t415 = pkin(4) * t443 - qJ(5) * t444;
t440 = qJDD(4) - t448;
t455 = -qJD(4) + t494;
t454 = t455 ^ 2;
t371 = -t440 * pkin(4) - t454 * qJ(5) + t444 * t415 + qJDD(5) - t373;
t410 = -t443 * qJD(4) + t464 * qJDD(3) + t502 * t447;
t417 = -mrSges(6,2) * t443 - mrSges(6,3) * t444;
t506 = -m(6) * t371 - t410 * mrSges(6,1) - t444 * t417;
t414 = -mrSges(7,2) * t444 + mrSges(7,3) * t443;
t499 = t443 * t455;
t503 = 2 * qJD(6);
t365 = t455 * t503 + (t443 * t444 - t440) * qJ(6) + (t410 - t499) * pkin(5) + t371;
t424 = -mrSges(7,1) * t443 - mrSges(7,2) * t455;
t480 = -m(7) * t365 + t440 * mrSges(7,3) - t455 * t424;
t363 = t410 * mrSges(7,1) + t444 * t414 - t480;
t423 = mrSges(6,1) * t443 + mrSges(6,3) * t455;
t360 = t440 * mrSges(6,2) - t455 * t423 + t363 - t506;
t409 = t444 * qJD(4) - t502 * qJDD(3) + t464 * t447;
t421 = pkin(5) * t444 + qJ(6) * t455;
t439 = t443 ^ 2;
t374 = t502 * t378 + t464 * t380;
t475 = -t454 * pkin(4) + t440 * qJ(5) - t443 * t415 + t374;
t504 = -2 * qJD(5);
t367 = -t409 * pkin(5) - t439 * qJ(6) + qJDD(6) + (t504 - t421) * t455 + t475;
t370 = 0.2e1 * qJD(5) * t455 - t475;
t425 = mrSges(6,1) * t444 - mrSges(6,2) * t455;
t422 = mrSges(7,1) * t444 + mrSges(7,3) * t455;
t488 = m(7) * t367 + t440 * mrSges(7,2) - t455 * t422;
t477 = -m(6) * t370 + t440 * mrSges(6,3) - t455 * t425 + t488;
t485 = -t492 * t443 + t509 * t444 - t491 * t455;
t486 = t508 * t443 + t492 * t444 + t490 * t455;
t496 = -t414 - t417;
t505 = t490 * t409 + t491 * t410 + t507 * t440 + t485 * t443 + t486 * t444 + mrSges(5,1) * t373 - mrSges(5,2) * t374 + mrSges(6,2) * t371 + mrSges(7,2) * t367 - mrSges(6,3) * t370 - mrSges(7,3) * t365 - pkin(4) * t360 + qJ(5) * (t496 * t443 + (-mrSges(6,1) - mrSges(7,1)) * t409 + t477) - qJ(6) * t363;
t500 = -mrSges(7,1) - mrSges(5,3);
t487 = -t490 * t443 - t491 * t444 + t507 * t455;
t445 = (-mrSges(4,1) * t467 + mrSges(4,2) * t465) * qJD(2);
t451 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t495;
t416 = mrSges(5,1) * t443 + mrSges(5,2) * t444;
t419 = mrSges(5,2) * t455 - mrSges(5,3) * t443;
t357 = m(5) * t373 + (-t419 + t423) * t455 + (-t414 - t416) * t444 + (mrSges(5,1) - mrSges(6,2)) * t440 + t500 * t410 + t480 + t506;
t420 = -mrSges(5,1) * t455 - mrSges(5,3) * t444;
t359 = m(5) * t374 - t440 * mrSges(5,2) + t455 * t420 + (-t416 + t496) * t443 + (-mrSges(6,1) + t500) * t409 + t477;
t481 = -t357 * t464 + t502 * t359;
t354 = m(4) * t382 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t448 - qJD(3) * t451 + t445 * t494 + t481;
t381 = -t465 * t386 + t467 * t427;
t452 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t494;
t377 = -qJDD(3) * pkin(3) - t469 * pkin(9) + t446 * t495 - t381;
t473 = (-t410 - t499) * qJ(5) + t377 + (-pkin(4) * t455 + t504) * t444;
t372 = t409 * pkin(4) + t473;
t369 = -t439 * pkin(5) + t443 * t503 - t444 * t421 + (pkin(4) + qJ(6)) * t409 + t473;
t479 = m(7) * t369 - t410 * mrSges(7,2) + t409 * mrSges(7,3) - t444 * t422 + t443 * t424;
t476 = -m(6) * t372 + t409 * mrSges(6,2) + t443 * t423 - t479;
t472 = -m(5) * t377 - t409 * mrSges(5,1) - t443 * t419 + (-t420 + t425) * t444 + (-mrSges(5,2) + mrSges(6,3)) * t410 + t476;
t356 = m(4) * t381 + qJDD(3) * mrSges(4,1) - t447 * mrSges(4,3) + qJD(3) * t452 - t445 * t495 + t472;
t482 = t467 * t354 - t356 * t465;
t355 = t502 * t357 + t464 * t359;
t474 = -m(4) * t385 + t448 * mrSges(4,1) - t447 * mrSges(4,2) - t451 * t495 + t452 * t494 - t355;
t433 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t465 + Ifges(4,4) * t467) * qJD(2);
t432 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t465 + Ifges(4,2) * t467) * qJD(2);
t364 = -t409 * mrSges(7,1) - t443 * t414 + t488;
t362 = -t410 * mrSges(6,3) - t444 * t425 - t476;
t352 = mrSges(6,1) * t371 + mrSges(7,1) * t365 + mrSges(5,2) * t377 - mrSges(7,2) * t369 - mrSges(5,3) * t373 - mrSges(6,3) * t372 + pkin(5) * t363 - qJ(5) * t362 - t492 * t409 + t509 * t410 + t491 * t440 + t487 * t443 + t486 * t455;
t351 = -mrSges(5,1) * t377 - mrSges(6,1) * t370 + mrSges(7,1) * t367 + mrSges(6,2) * t372 + mrSges(5,3) * t374 - mrSges(7,3) * t369 - pkin(4) * t362 + pkin(5) * t364 - qJ(6) * t479 + t508 * t409 + t492 * t410 - t490 * t440 + t487 * t444 - t485 * t455;
t1 = [m(2) * t459 + t463 * (m(3) * t427 + t354 * t465 + t356 * t467) + (t466 * (m(3) * t390 - mrSges(3,1) * t470 - qJDD(2) * mrSges(3,2) + t482) + t468 * (m(3) * t389 + qJDD(2) * mrSges(3,1) - t470 * mrSges(3,2) + t474)) * t461; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t389 - mrSges(3,2) * t390 + t465 * (mrSges(4,2) * t385 - mrSges(4,3) * t381 + Ifges(4,1) * t447 + Ifges(4,4) * t448 + Ifges(4,5) * qJDD(3) - pkin(9) * t355 - qJD(3) * t432 - t464 * t351 + t502 * t352) + pkin(2) * t474 + pkin(8) * t482 + (-mrSges(4,1) * t385 + mrSges(4,3) * t382 + Ifges(4,4) * t447 + Ifges(4,2) * t448 + Ifges(4,6) * qJDD(3) - pkin(3) * t355 + qJD(3) * t433 - t505) * t467; Ifges(4,5) * t447 + Ifges(4,6) * t448 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t381 - mrSges(4,2) * t382 + t464 * t352 + t502 * t351 + pkin(3) * t472 + pkin(9) * t481 + (t432 * t465 - t433 * t467) * qJD(2); t505; t360; t364;];
tauJ  = t1;
