% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:18:11
% EndTime: 2019-12-31 19:18:18
% DurationCPUTime: 6.94s
% Computational Cost: add. (82794->282), mult. (258005->401), div. (0->0), fcn. (215445->14), ass. (0->139)
t452 = sin(pkin(11));
t454 = sin(pkin(5));
t455 = cos(pkin(11));
t457 = cos(pkin(5));
t460 = sin(qJ(3));
t456 = cos(pkin(6));
t464 = cos(qJ(3));
t492 = t456 * t464;
t453 = sin(pkin(6));
t497 = t453 * t464;
t470 = (-t452 * t460 + t455 * t492) * t454 + t457 * t497;
t429 = t470 * qJD(1);
t493 = t456 * t460;
t498 = t453 * t460;
t472 = t457 * t498 + (t452 * t464 + t455 * t493) * t454;
t430 = t472 * qJD(1);
t418 = -t430 * qJD(3) + t470 * qJDD(1);
t494 = t454 * t456;
t441 = (t453 * t457 + t455 * t494) * qJD(1) * pkin(8);
t466 = qJD(1) ^ 2;
t461 = sin(qJ(1));
t465 = cos(qJ(1));
t483 = -t465 * g(1) - t461 * g(2);
t496 = t454 * qJ(2);
t445 = -t466 * pkin(1) + qJDD(1) * t496 + t483;
t502 = pkin(8) * t453;
t479 = -t455 * pkin(2) - t452 * t502;
t491 = qJD(1) * t454;
t501 = pkin(8) * qJDD(1);
t474 = qJD(1) * t479 * t491 + t456 * t501;
t487 = qJD(2) * t491;
t495 = t454 * t455;
t486 = t461 * g(1) - t465 * g(2);
t444 = qJDD(1) * pkin(1) + t466 * t496 + t486;
t500 = t444 * t457;
t480 = -g(3) * t495 - 0.2e1 * t452 * t487 + t455 * t500;
t400 = (pkin(2) * qJDD(1) + qJD(1) * t441) * t457 + (-t474 * t454 - t445) * t452 + t480;
t503 = pkin(8) * t452;
t446 = (t457 * pkin(2) - t494 * t503) * qJD(1);
t504 = 0.2e1 * t455;
t488 = t455 * t445 + t452 * t500 + t487 * t504;
t401 = (-qJD(1) * t446 + t453 * t501) * t457 + (-g(3) * t452 + t474 * t455) * t454 + t488;
t485 = -t457 * g(3) + qJDD(2);
t409 = (-t444 + t479 * qJDD(1) + (-t441 * t455 + t446 * t452) * qJD(1)) * t454 + t485;
t376 = -t460 * t401 + (t400 * t456 + t409 * t453) * t464;
t499 = t452 * t454;
t377 = t400 * t493 + t464 * t401 + t409 * t498;
t417 = -t429 * pkin(3) - t430 * pkin(9);
t475 = -t453 * t495 + t457 * t456;
t442 = t475 * qJD(1) + qJD(3);
t438 = t442 ^ 2;
t439 = t475 * qJDD(1) + qJDD(3);
t373 = -t438 * pkin(3) + t439 * pkin(9) + t429 * t417 + t377;
t385 = -t453 * t400 + t456 * t409;
t419 = t429 * qJD(3) + t472 * qJDD(1);
t375 = (-t429 * t442 - t419) * pkin(9) + (t430 * t442 - t418) * pkin(3) + t385;
t459 = sin(qJ(4));
t463 = cos(qJ(4));
t370 = t463 * t373 + t459 * t375;
t423 = -t459 * t430 + t463 * t442;
t424 = t463 * t430 + t459 * t442;
t403 = -t423 * pkin(4) - t424 * pkin(10);
t415 = qJDD(4) - t418;
t428 = qJD(4) - t429;
t427 = t428 ^ 2;
t367 = -t427 * pkin(4) + t415 * pkin(10) + t423 * t403 + t370;
t372 = -t439 * pkin(3) - t438 * pkin(9) + t430 * t417 - t376;
t395 = -t424 * qJD(4) - t459 * t419 + t463 * t439;
t396 = t423 * qJD(4) + t463 * t419 + t459 * t439;
t368 = (-t423 * t428 - t396) * pkin(10) + (t424 * t428 - t395) * pkin(4) + t372;
t458 = sin(qJ(5));
t462 = cos(qJ(5));
t364 = -t458 * t367 + t462 * t368;
t407 = -t458 * t424 + t462 * t428;
t380 = t407 * qJD(5) + t462 * t396 + t458 * t415;
t408 = t462 * t424 + t458 * t428;
t386 = -t407 * mrSges(6,1) + t408 * mrSges(6,2);
t420 = qJD(5) - t423;
t387 = -t420 * mrSges(6,2) + t407 * mrSges(6,3);
t393 = qJDD(5) - t395;
t361 = m(6) * t364 + t393 * mrSges(6,1) - t380 * mrSges(6,3) - t408 * t386 + t420 * t387;
t365 = t462 * t367 + t458 * t368;
t379 = -t408 * qJD(5) - t458 * t396 + t462 * t415;
t388 = t420 * mrSges(6,1) - t408 * mrSges(6,3);
t362 = m(6) * t365 - t393 * mrSges(6,2) + t379 * mrSges(6,3) + t407 * t386 - t420 * t388;
t355 = -t458 * t361 + t462 * t362;
t402 = -t423 * mrSges(5,1) + t424 * mrSges(5,2);
t411 = t428 * mrSges(5,1) - t424 * mrSges(5,3);
t353 = m(5) * t370 - t415 * mrSges(5,2) + t395 * mrSges(5,3) + t423 * t402 - t428 * t411 + t355;
t369 = -t459 * t373 + t463 * t375;
t366 = -t415 * pkin(4) - t427 * pkin(10) + t424 * t403 - t369;
t363 = -m(6) * t366 + t379 * mrSges(6,1) - t380 * mrSges(6,2) + t407 * t387 - t408 * t388;
t410 = -t428 * mrSges(5,2) + t423 * mrSges(5,3);
t359 = m(5) * t369 + t415 * mrSges(5,1) - t396 * mrSges(5,3) - t424 * t402 + t428 * t410 + t363;
t347 = t459 * t353 + t463 * t359;
t416 = -t429 * mrSges(4,1) + t430 * mrSges(4,2);
t426 = t442 * mrSges(4,1) - t430 * mrSges(4,3);
t484 = t463 * t353 - t459 * t359;
t344 = m(4) * t377 - t439 * mrSges(4,2) + t418 * mrSges(4,3) + t429 * t416 - t442 * t426 + t484;
t425 = -t442 * mrSges(4,2) + t429 * mrSges(4,3);
t346 = m(4) * t385 - t418 * mrSges(4,1) + t419 * mrSges(4,2) - t429 * t425 + t430 * t426 + t347;
t354 = t462 * t361 + t458 * t362;
t469 = -m(5) * t372 + t395 * mrSges(5,1) - t396 * mrSges(5,2) + t423 * t410 - t424 * t411 - t354;
t350 = m(4) * t376 + t439 * mrSges(4,1) - t419 * mrSges(4,3) - t430 * t416 + t442 * t425 + t469;
t335 = t344 * t498 + t456 * t346 + t350 * t497;
t338 = t464 * t344 - t460 * t350;
t336 = t344 * t493 - t453 * t346 + t350 * t492;
t482 = -t455 * mrSges(3,1) + t452 * mrSges(3,2);
t478 = t457 * mrSges(3,1) - mrSges(3,3) * t499;
t477 = -t457 * mrSges(3,2) + mrSges(3,3) * t495;
t381 = Ifges(6,5) * t408 + Ifges(6,6) * t407 + Ifges(6,3) * t420;
t383 = Ifges(6,1) * t408 + Ifges(6,4) * t407 + Ifges(6,5) * t420;
t356 = -mrSges(6,1) * t366 + mrSges(6,3) * t365 + Ifges(6,4) * t380 + Ifges(6,2) * t379 + Ifges(6,6) * t393 - t408 * t381 + t420 * t383;
t382 = Ifges(6,4) * t408 + Ifges(6,2) * t407 + Ifges(6,6) * t420;
t357 = mrSges(6,2) * t366 - mrSges(6,3) * t364 + Ifges(6,1) * t380 + Ifges(6,4) * t379 + Ifges(6,5) * t393 + t407 * t381 - t420 * t382;
t389 = Ifges(5,5) * t424 + Ifges(5,6) * t423 + Ifges(5,3) * t428;
t390 = Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * t428;
t339 = mrSges(5,2) * t372 - mrSges(5,3) * t369 + Ifges(5,1) * t396 + Ifges(5,4) * t395 + Ifges(5,5) * t415 - pkin(10) * t354 - t458 * t356 + t462 * t357 + t423 * t389 - t428 * t390;
t391 = Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * t428;
t468 = mrSges(6,1) * t364 - mrSges(6,2) * t365 + Ifges(6,5) * t380 + Ifges(6,6) * t379 + Ifges(6,3) * t393 + t408 * t382 - t407 * t383;
t340 = -mrSges(5,1) * t372 + mrSges(5,3) * t370 + Ifges(5,4) * t396 + Ifges(5,2) * t395 + Ifges(5,6) * t415 - pkin(4) * t354 - t424 * t389 + t428 * t391 - t468;
t412 = Ifges(4,5) * t430 + Ifges(4,6) * t429 + Ifges(4,3) * t442;
t413 = Ifges(4,4) * t430 + Ifges(4,2) * t429 + Ifges(4,6) * t442;
t331 = mrSges(4,2) * t385 - mrSges(4,3) * t376 + Ifges(4,1) * t419 + Ifges(4,4) * t418 + Ifges(4,5) * t439 - pkin(9) * t347 + t463 * t339 - t459 * t340 + t429 * t412 - t442 * t413;
t414 = Ifges(4,1) * t430 + Ifges(4,4) * t429 + Ifges(4,5) * t442;
t467 = mrSges(5,1) * t369 - mrSges(5,2) * t370 + Ifges(5,5) * t396 + Ifges(5,6) * t395 + Ifges(5,3) * t415 + pkin(4) * t363 + pkin(10) * t355 + t462 * t356 + t458 * t357 + t424 * t390 - t423 * t391;
t332 = -mrSges(4,1) * t385 + mrSges(4,3) * t377 + Ifges(4,4) * t419 + Ifges(4,2) * t418 + Ifges(4,6) * t439 - pkin(3) * t347 - t430 * t412 + t442 * t414 - t467;
t473 = pkin(8) * t338 + t331 * t460 + t332 * t464;
t448 = t477 * qJD(1);
t447 = t478 * qJD(1);
t443 = t482 * t491;
t431 = -t454 * t444 + t485;
t422 = -g(3) * t499 + t488;
t421 = -t452 * t445 + t480;
t337 = m(3) * t422 + t477 * qJDD(1) + (t443 * t495 - t447 * t457) * qJD(1) + t338;
t334 = m(3) * t431 + (t482 * qJDD(1) + (t447 * t452 - t448 * t455) * qJD(1)) * t454 + t335;
t333 = m(3) * t421 + t478 * qJDD(1) + (-t443 * t499 + t448 * t457) * qJD(1) + t336;
t330 = mrSges(4,1) * t376 - mrSges(4,2) * t377 + Ifges(4,5) * t419 + Ifges(4,6) * t418 + Ifges(4,3) * t439 + pkin(3) * t469 + pkin(9) * t484 + t459 * t339 + t463 * t340 + t430 * t413 - t429 * t414;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t486 - mrSges(2,2) * t483 + (mrSges(3,1) * t421 - mrSges(3,2) * t422 + pkin(2) * t336 + t456 * t330 + pkin(1) * (t333 * t455 + t337 * t452) + Ifges(3,3) * t457 * qJDD(1) + t473 * t453) * t457 + (t452 * (mrSges(3,2) * t431 - mrSges(3,3) * t421 + t464 * t331 - t460 * t332 - t335 * t502) + t455 * (-mrSges(3,1) * t431 + mrSges(3,3) * t422 - pkin(2) * t335 - t453 * t330) - pkin(1) * t334 + qJ(2) * (-t452 * t333 + t455 * t337) + (-t336 * t503 + t455 * t473) * t456 + ((Ifges(3,2) * t455 ^ 2 + (Ifges(3,1) * t452 + Ifges(3,4) * t504) * t452) * t454 + 0.2e1 * t457 * (Ifges(3,5) * t452 + Ifges(3,6) * t455)) * qJDD(1)) * t454; t334; t330; t467; t468;];
tauJ = t1;
