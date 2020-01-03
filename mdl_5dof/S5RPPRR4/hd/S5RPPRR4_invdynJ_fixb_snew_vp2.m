% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:01
% EndTime: 2020-01-03 11:31:04
% DurationCPUTime: 2.39s
% Computational Cost: add. (13893->231), mult. (38009->320), div. (0->0), fcn. (26039->10), ass. (0->112)
t439 = qJD(1) ^ 2;
t435 = sin(qJ(1));
t438 = cos(qJ(1));
t459 = -g(2) * t435 + t438 * g(3);
t476 = -pkin(1) * t439 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t459;
t468 = -t438 * g(2) - t435 * g(3);
t447 = -qJ(2) * t439 + qJDD(2) - t468;
t430 = sin(pkin(8));
t432 = cos(pkin(8));
t453 = -pkin(2) * t432 - qJ(3) * t430;
t467 = qJD(1) * t430;
t475 = (-pkin(1) + t453) * qJDD(1) + t447 - 0.2e1 * qJD(3) * t467;
t391 = -t432 * g(1) - t476 * t430;
t431 = cos(pkin(9));
t474 = Ifges(4,4) * t431;
t429 = sin(pkin(9));
t473 = Ifges(4,2) * t429;
t472 = t429 * t430;
t471 = t430 * t431;
t392 = -g(1) * t430 + t476 * t432;
t410 = t453 * qJD(1);
t466 = qJD(1) * t432;
t381 = t410 * t466 + t392;
t428 = t430 ^ 2;
t451 = -pkin(3) * t432 - pkin(6) * t471;
t470 = t475 * t431;
t360 = t451 * qJDD(1) + (-t381 + (-pkin(3) * t428 * t431 + pkin(6) * t430 * t432) * t439) * t429 + t470;
t368 = t431 * t381 + t475 * t429;
t409 = t451 * qJD(1);
t464 = qJDD(1) * t430;
t460 = t429 * t464;
t462 = t429 ^ 2 * t428 * t439;
t361 = -pkin(3) * t462 - pkin(6) * t460 + t409 * t466 + t368;
t434 = sin(qJ(4));
t437 = cos(qJ(4));
t347 = t437 * t360 - t361 * t434;
t446 = t430 * (-t429 * t437 - t431 * t434);
t399 = qJD(1) * t446;
t445 = t430 * (-t429 * t434 + t431 * t437);
t385 = qJD(4) * t399 + qJDD(1) * t445;
t400 = qJD(1) * t445;
t463 = qJDD(1) * t432;
t419 = qJDD(4) - t463;
t420 = qJD(4) - t466;
t345 = (t399 * t420 - t385) * pkin(7) + (t399 * t400 + t419) * pkin(4) + t347;
t348 = t434 * t360 + t437 * t361;
t384 = -qJD(4) * t400 + qJDD(1) * t446;
t390 = pkin(4) * t420 - pkin(7) * t400;
t398 = t399 ^ 2;
t346 = -pkin(4) * t398 + pkin(7) * t384 - t390 * t420 + t348;
t433 = sin(qJ(5));
t436 = cos(qJ(5));
t343 = t345 * t436 - t346 * t433;
t378 = t399 * t436 - t400 * t433;
t356 = qJD(5) * t378 + t384 * t433 + t385 * t436;
t379 = t399 * t433 + t400 * t436;
t366 = -mrSges(6,1) * t378 + mrSges(6,2) * t379;
t418 = qJD(5) + t420;
t371 = -mrSges(6,2) * t418 + mrSges(6,3) * t378;
t416 = qJDD(5) + t419;
t340 = m(6) * t343 + mrSges(6,1) * t416 - mrSges(6,3) * t356 - t366 * t379 + t371 * t418;
t344 = t345 * t433 + t346 * t436;
t355 = -qJD(5) * t379 + t384 * t436 - t385 * t433;
t372 = mrSges(6,1) * t418 - mrSges(6,3) * t379;
t341 = m(6) * t344 - mrSges(6,2) * t416 + mrSges(6,3) * t355 + t366 * t378 - t372 * t418;
t334 = t436 * t340 + t433 * t341;
t382 = -mrSges(5,1) * t399 + mrSges(5,2) * t400;
t386 = -mrSges(5,2) * t420 + mrSges(5,3) * t399;
t332 = m(5) * t347 + mrSges(5,1) * t419 - mrSges(5,3) * t385 - t382 * t400 + t386 * t420 + t334;
t387 = mrSges(5,1) * t420 - mrSges(5,3) * t400;
t457 = -t340 * t433 + t436 * t341;
t333 = m(5) * t348 - mrSges(5,2) * t419 + mrSges(5,3) * t384 + t382 * t399 - t387 * t420 + t457;
t328 = t437 * t332 + t434 * t333;
t458 = -t332 * t434 + t437 * t333;
t455 = -mrSges(3,1) * t432 + mrSges(3,2) * t430;
t454 = mrSges(4,1) * t429 + mrSges(4,2) * t431;
t367 = -t381 * t429 + t470;
t403 = t454 * t467;
t448 = mrSges(4,2) * t432 - mrSges(4,3) * t472;
t406 = t448 * qJD(1);
t449 = -mrSges(4,1) * t432 - mrSges(4,3) * t471;
t326 = m(4) * t367 + t449 * qJDD(1) + (-t403 * t471 - t406 * t432) * qJD(1) + t328;
t407 = t449 * qJD(1);
t327 = m(4) * t368 + t448 * qJDD(1) + (-t403 * t472 + t407 * t432) * qJD(1) + t458;
t323 = t326 * t431 + t327 * t429;
t452 = t406 * t429 + t407 * t431;
t380 = t410 * t467 + qJDD(3) - t391;
t450 = -Ifges(4,5) * t431 + Ifges(4,6) * t429 + Ifges(3,4);
t370 = t431 * t409 * t467 + pkin(3) * t460 - pkin(6) * t462 + t380;
t350 = -pkin(4) * t384 - pkin(7) * t398 + t390 * t400 + t370;
t444 = m(6) * t350 - t355 * mrSges(6,1) + t356 * mrSges(6,2) - t378 * t371 + t379 * t372;
t363 = Ifges(6,4) * t379 + Ifges(6,2) * t378 + Ifges(6,6) * t418;
t364 = Ifges(6,1) * t379 + Ifges(6,4) * t378 + Ifges(6,5) * t418;
t443 = -mrSges(6,1) * t343 + mrSges(6,2) * t344 - Ifges(6,5) * t356 - Ifges(6,6) * t355 - Ifges(6,3) * t416 - t379 * t363 + t378 * t364;
t442 = m(5) * t370 - t384 * mrSges(5,1) + t385 * mrSges(5,2) - t399 * t386 + t400 * t387 + t444;
t441 = m(4) * t380 + t442;
t374 = Ifges(5,4) * t400 + Ifges(5,2) * t399 + Ifges(5,6) * t420;
t375 = Ifges(5,1) * t400 + Ifges(5,4) * t399 + Ifges(5,5) * t420;
t440 = mrSges(5,1) * t347 - mrSges(5,2) * t348 + Ifges(5,5) * t385 + Ifges(5,6) * t384 + Ifges(5,3) * t419 + pkin(4) * t334 + t400 * t374 - t399 * t375 - t443;
t412 = (Ifges(3,5) * t430 + Ifges(3,6) * t432) * qJD(1);
t411 = t455 * qJD(1);
t405 = -qJDD(1) * pkin(1) + t447;
t396 = (-Ifges(4,5) * t432 + (Ifges(4,1) * t431 - Ifges(4,4) * t429) * t430) * qJD(1);
t395 = (-Ifges(4,6) * t432 + (-t473 + t474) * t430) * qJD(1);
t373 = Ifges(5,5) * t400 + Ifges(5,6) * t399 + Ifges(5,3) * t420;
t362 = Ifges(6,5) * t379 + Ifges(6,6) * t378 + Ifges(6,3) * t418;
t336 = mrSges(6,2) * t350 - mrSges(6,3) * t343 + Ifges(6,1) * t356 + Ifges(6,4) * t355 + Ifges(6,5) * t416 + t362 * t378 - t363 * t418;
t335 = -mrSges(6,1) * t350 + mrSges(6,3) * t344 + Ifges(6,4) * t356 + Ifges(6,2) * t355 + Ifges(6,6) * t416 - t362 * t379 + t364 * t418;
t325 = mrSges(5,2) * t370 - mrSges(5,3) * t347 + Ifges(5,1) * t385 + Ifges(5,4) * t384 + Ifges(5,5) * t419 - pkin(7) * t334 - t335 * t433 + t336 * t436 + t373 * t399 - t374 * t420;
t324 = -mrSges(5,1) * t370 + mrSges(5,3) * t348 + Ifges(5,4) * t385 + Ifges(5,2) * t384 + Ifges(5,6) * t419 - pkin(4) * t444 + pkin(7) * t457 + t436 * t335 + t433 * t336 - t400 * t373 + t420 * t375;
t322 = m(3) * t405 + t455 * qJDD(1) + (-t432 ^ 2 - t428) * t439 * mrSges(3,3) + t323;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t468 - mrSges(2,2) * t459 + t430 * (t412 * t466 + mrSges(3,2) * t405 - mrSges(3,3) * t391 + t431 * (mrSges(4,2) * t380 - mrSges(4,3) * t367 - pkin(6) * t328 - t324 * t434 + t437 * t325 + t395 * t466) - t429 * (-mrSges(4,1) * t380 + mrSges(4,3) * t368 - pkin(3) * t442 + pkin(6) * t458 + t437 * t324 + t434 * t325 - t396 * t466) - qJ(3) * t323 + (t432 * t450 + (Ifges(4,1) * t431 ^ 2 + Ifges(3,1) + (t473 - 0.2e1 * t474) * t429) * t430) * qJDD(1)) + t432 * (-t440 + (t450 * qJDD(1) + (-t395 * t431 - t396 * t429 - t412) * qJD(1)) * t430 + (Ifges(3,2) + Ifges(4,3)) * t463 - mrSges(3,1) * t405 + mrSges(3,3) * t392 - mrSges(4,1) * t367 + mrSges(4,2) * t368 - pkin(3) * t328 - pkin(2) * t323) - pkin(1) * t322 + qJ(2) * ((m(3) * t392 - t429 * t326 + t431 * t327 + (qJDD(1) * mrSges(3,3) + qJD(1) * t411) * t432) * t432 + (t441 - m(3) * t391 + (t411 + t452) * t467 + (mrSges(3,3) + t454) * t464) * t430); t322; (t452 * qJD(1) + t454 * qJDD(1)) * t430 + t441; t440; -t443;];
tauJ = t1;
