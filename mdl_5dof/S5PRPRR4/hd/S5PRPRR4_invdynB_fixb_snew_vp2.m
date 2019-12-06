% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:50:09
% EndTime: 2019-12-05 15:50:14
% DurationCPUTime: 3.56s
% Computational Cost: add. (45148->225), mult. (81573->293), div. (0->0), fcn. (55733->12), ass. (0->103)
t440 = sin(pkin(5));
t446 = sin(qJ(2));
t470 = t440 * t446;
t449 = cos(qJ(2));
t469 = t440 * t449;
t443 = cos(pkin(5));
t468 = t443 * t446;
t467 = t443 * t449;
t439 = sin(pkin(9));
t442 = cos(pkin(9));
t429 = g(1) * t439 - g(2) * t442;
t437 = -g(3) + qJDD(1);
t413 = -t429 * t440 + t443 * t437;
t412 = qJDD(3) + t413;
t448 = cos(qJ(4));
t466 = t448 * t412;
t430 = -g(1) * t442 - g(2) * t439;
t398 = t429 * t467 - t430 * t446 + t437 * t469;
t396 = qJDD(2) * pkin(2) + t398;
t399 = t429 * t468 + t449 * t430 + t437 * t470;
t451 = qJD(2) ^ 2;
t397 = -pkin(2) * t451 + t399;
t438 = sin(pkin(10));
t441 = cos(pkin(10));
t392 = t438 * t396 + t441 * t397;
t390 = -pkin(3) * t451 + qJDD(2) * pkin(7) + t392;
t445 = sin(qJ(4));
t387 = t448 * t390 + t445 * t412;
t425 = (-mrSges(5,1) * t448 + mrSges(5,2) * t445) * qJD(2);
t462 = qJD(2) * qJD(4);
t460 = t445 * t462;
t428 = qJDD(2) * t448 - t460;
t464 = qJD(2) * t445;
t431 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t464;
t426 = (-pkin(4) * t448 - pkin(8) * t445) * qJD(2);
t450 = qJD(4) ^ 2;
t463 = qJD(2) * t448;
t384 = -pkin(4) * t450 + qJDD(4) * pkin(8) + t426 * t463 + t387;
t391 = t441 * t396 - t438 * t397;
t389 = -qJDD(2) * pkin(3) - t451 * pkin(7) - t391;
t459 = t448 * t462;
t427 = qJDD(2) * t445 + t459;
t385 = (-t427 - t459) * pkin(8) + (-t428 + t460) * pkin(4) + t389;
t444 = sin(qJ(5));
t447 = cos(qJ(5));
t381 = -t384 * t444 + t385 * t447;
t423 = qJD(4) * t447 - t444 * t464;
t406 = qJD(5) * t423 + qJDD(4) * t444 + t427 * t447;
t424 = qJD(4) * t444 + t447 * t464;
t407 = -mrSges(6,1) * t423 + mrSges(6,2) * t424;
t435 = qJD(5) - t463;
t410 = -mrSges(6,2) * t435 + mrSges(6,3) * t423;
t421 = qJDD(5) - t428;
t379 = m(6) * t381 + mrSges(6,1) * t421 - mrSges(6,3) * t406 - t407 * t424 + t410 * t435;
t382 = t384 * t447 + t385 * t444;
t405 = -qJD(5) * t424 + qJDD(4) * t447 - t427 * t444;
t411 = mrSges(6,1) * t435 - mrSges(6,3) * t424;
t380 = m(6) * t382 - mrSges(6,2) * t421 + mrSges(6,3) * t405 + t407 * t423 - t411 * t435;
t455 = -t379 * t444 + t447 * t380;
t372 = m(5) * t387 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t428 - qJD(4) * t431 + t425 * t463 + t455;
t386 = -t390 * t445 + t466;
t432 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t463;
t383 = -qJDD(4) * pkin(4) - t450 * pkin(8) - t466 + (qJD(2) * t426 + t390) * t445;
t453 = -m(6) * t383 + t405 * mrSges(6,1) - mrSges(6,2) * t406 + t423 * t410 - t411 * t424;
t377 = m(5) * t386 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t427 + qJD(4) * t432 - t425 * t464 + t453;
t456 = t448 * t372 - t377 * t445;
t364 = m(4) * t392 - mrSges(4,1) * t451 - qJDD(2) * mrSges(4,2) + t456;
t373 = t379 * t447 + t380 * t444;
t452 = -m(5) * t389 + t428 * mrSges(5,1) - mrSges(5,2) * t427 - t431 * t464 + t432 * t463 - t373;
t369 = m(4) * t391 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t451 + t452;
t359 = t438 * t364 + t441 * t369;
t357 = m(3) * t398 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t451 + t359;
t457 = t441 * t364 - t369 * t438;
t358 = m(3) * t399 - mrSges(3,1) * t451 - qJDD(2) * mrSges(3,2) + t457;
t367 = t445 * t372 + t448 * t377;
t461 = m(4) * t412 + t367;
t366 = m(3) * t413 + t461;
t345 = t357 * t467 + t358 * t468 - t366 * t440;
t343 = m(2) * t429 + t345;
t350 = -t357 * t446 + t449 * t358;
t349 = m(2) * t430 + t350;
t465 = t442 * t343 + t439 * t349;
t344 = t357 * t469 + t358 * t470 + t443 * t366;
t458 = -t343 * t439 + t442 * t349;
t400 = Ifges(6,5) * t424 + Ifges(6,6) * t423 + Ifges(6,3) * t435;
t402 = Ifges(6,1) * t424 + Ifges(6,4) * t423 + Ifges(6,5) * t435;
t374 = -mrSges(6,1) * t383 + mrSges(6,3) * t382 + Ifges(6,4) * t406 + Ifges(6,2) * t405 + Ifges(6,6) * t421 - t400 * t424 + t402 * t435;
t401 = Ifges(6,4) * t424 + Ifges(6,2) * t423 + Ifges(6,6) * t435;
t375 = mrSges(6,2) * t383 - mrSges(6,3) * t381 + Ifges(6,1) * t406 + Ifges(6,4) * t405 + Ifges(6,5) * t421 + t400 * t423 - t401 * t435;
t416 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t445 + Ifges(5,6) * t448) * qJD(2);
t417 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t445 + Ifges(5,2) * t448) * qJD(2);
t360 = mrSges(5,2) * t389 - mrSges(5,3) * t386 + Ifges(5,1) * t427 + Ifges(5,4) * t428 + Ifges(5,5) * qJDD(4) - pkin(8) * t373 - qJD(4) * t417 - t374 * t444 + t375 * t447 + t416 * t463;
t418 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t445 + Ifges(5,4) * t448) * qJD(2);
t361 = -mrSges(5,1) * t389 - mrSges(6,1) * t381 + mrSges(6,2) * t382 + mrSges(5,3) * t387 + Ifges(5,4) * t427 - Ifges(6,5) * t406 + Ifges(5,2) * t428 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t405 - Ifges(6,3) * t421 - pkin(4) * t373 + qJD(4) * t418 - t401 * t424 + t402 * t423 - t416 * t464;
t346 = mrSges(4,2) * t412 - mrSges(4,3) * t391 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t451 - pkin(7) * t367 + t360 * t448 - t361 * t445;
t351 = Ifges(4,6) * qJDD(2) + t451 * Ifges(4,5) - mrSges(4,1) * t412 + mrSges(4,3) * t392 - Ifges(5,5) * t427 - Ifges(5,6) * t428 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t386 + mrSges(5,2) * t387 - t444 * t375 - t447 * t374 - pkin(4) * t453 - pkin(8) * t455 - pkin(3) * t367 + (-t417 * t445 + t418 * t448) * qJD(2);
t339 = -mrSges(3,1) * t413 + mrSges(3,3) * t399 + t451 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t461 + qJ(3) * t457 + t438 * t346 + t441 * t351;
t340 = mrSges(3,2) * t413 - mrSges(3,3) * t398 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t451 - qJ(3) * t359 + t346 * t441 - t351 * t438;
t454 = pkin(6) * t350 + t339 * t449 + t340 * t446;
t341 = mrSges(3,1) * t398 - mrSges(3,2) * t399 + mrSges(4,1) * t391 - mrSges(4,2) * t392 + t445 * t360 + t448 * t361 + pkin(3) * t452 + pkin(7) * t456 + pkin(2) * t359 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t338 = mrSges(2,2) * t437 - mrSges(2,3) * t429 - t446 * t339 + t449 * t340 + (-t344 * t440 - t345 * t443) * pkin(6);
t337 = -mrSges(2,1) * t437 + mrSges(2,3) * t430 - pkin(1) * t344 - t440 * t341 + t454 * t443;
t1 = [-m(1) * g(1) + t458; -m(1) * g(2) + t465; -m(1) * g(3) + m(2) * t437 + t344; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t465 - t439 * t337 + t442 * t338; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t458 + t442 * t337 + t439 * t338; -mrSges(1,1) * g(2) + mrSges(2,1) * t429 + mrSges(1,2) * g(1) - mrSges(2,2) * t430 + pkin(1) * t345 + t443 * t341 + t454 * t440;];
tauB = t1;
