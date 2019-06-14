% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:31:05
% EndTime: 2019-05-05 01:31:07
% DurationCPUTime: 1.42s
% Computational Cost: add. (11932->234), mult. (22392->293), div. (0->0), fcn. (14930->12), ass. (0->102)
t430 = sin(pkin(11));
t432 = cos(pkin(11));
t414 = g(1) * t430 - g(2) * t432;
t427 = -g(3) + qJDD(1);
t431 = sin(pkin(6));
t433 = cos(pkin(6));
t464 = t414 * t433 + t427 * t431;
t415 = -g(1) * t432 - g(2) * t430;
t437 = sin(qJ(2));
t441 = cos(qJ(2));
t377 = t441 * t415 + t464 * t437;
t463 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t377;
t376 = -t437 * t415 + t441 * t464;
t462 = -pkin(2) - pkin(8);
t443 = qJD(2) ^ 2;
t447 = -t443 * qJ(3) + qJDD(3) - t376;
t371 = qJDD(2) * t462 + t447;
t393 = -t414 * t431 + t427 * t433;
t436 = sin(qJ(4));
t440 = cos(qJ(4));
t365 = t436 * t371 + t440 * t393;
t411 = (pkin(4) * t436 - pkin(9) * t440) * qJD(2);
t423 = t436 * qJD(2);
t442 = qJD(4) ^ 2;
t357 = -pkin(4) * t442 + qJDD(4) * pkin(9) - t411 * t423 + t365;
t370 = t443 * t462 - t463;
t458 = qJD(2) * qJD(4);
t421 = t440 * t458;
t412 = -t436 * qJDD(2) - t421;
t456 = t436 * t458;
t413 = qJDD(2) * t440 - t456;
t360 = (-t413 + t456) * pkin(9) + (-t412 + t421) * pkin(4) + t370;
t435 = sin(qJ(5));
t439 = cos(qJ(5));
t346 = -t435 * t357 + t439 * t360;
t459 = qJD(2) * t440;
t408 = qJD(4) * t439 - t435 * t459;
t384 = qJD(5) * t408 + qJDD(4) * t435 + t413 * t439;
t405 = qJDD(5) - t412;
t409 = qJD(4) * t435 + t439 * t459;
t420 = t423 + qJD(5);
t344 = (t408 * t420 - t384) * pkin(10) + (t408 * t409 + t405) * pkin(5) + t346;
t347 = t439 * t357 + t435 * t360;
t383 = -qJD(5) * t409 + qJDD(4) * t439 - t413 * t435;
t392 = pkin(5) * t420 - pkin(10) * t409;
t404 = t408 ^ 2;
t345 = -pkin(5) * t404 + pkin(10) * t383 - t392 * t420 + t347;
t434 = sin(qJ(6));
t438 = cos(qJ(6));
t342 = t344 * t438 - t345 * t434;
t385 = t408 * t438 - t409 * t434;
t354 = qJD(6) * t385 + t383 * t434 + t384 * t438;
t386 = t408 * t434 + t409 * t438;
t367 = -mrSges(7,1) * t385 + mrSges(7,2) * t386;
t419 = qJD(6) + t420;
t374 = -mrSges(7,2) * t419 + mrSges(7,3) * t385;
t400 = qJDD(6) + t405;
t339 = m(7) * t342 + mrSges(7,1) * t400 - mrSges(7,3) * t354 - t367 * t386 + t374 * t419;
t343 = t344 * t434 + t345 * t438;
t353 = -qJD(6) * t386 + t383 * t438 - t384 * t434;
t375 = mrSges(7,1) * t419 - mrSges(7,3) * t386;
t340 = m(7) * t343 - mrSges(7,2) * t400 + mrSges(7,3) * t353 + t367 * t385 - t375 * t419;
t332 = t438 * t339 + t434 * t340;
t387 = -mrSges(6,1) * t408 + mrSges(6,2) * t409;
t390 = -mrSges(6,2) * t420 + mrSges(6,3) * t408;
t330 = m(6) * t346 + mrSges(6,1) * t405 - mrSges(6,3) * t384 - t387 * t409 + t390 * t420 + t332;
t391 = mrSges(6,1) * t420 - mrSges(6,3) * t409;
t454 = -t339 * t434 + t438 * t340;
t331 = m(6) * t347 - mrSges(6,2) * t405 + mrSges(6,3) * t383 + t387 * t408 - t391 * t420 + t454;
t326 = t439 * t330 + t435 * t331;
t455 = -t330 * t435 + t439 * t331;
t364 = t371 * t440 - t436 * t393;
t410 = (mrSges(5,1) * t436 + mrSges(5,2) * t440) * qJD(2);
t417 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t459;
t325 = m(5) * t365 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t412 - qJD(4) * t417 - t410 * t423 + t455;
t416 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t423;
t356 = -qJDD(4) * pkin(4) - pkin(9) * t442 + t411 * t459 - t364;
t348 = -pkin(5) * t383 - pkin(10) * t404 + t392 * t409 + t356;
t449 = m(7) * t348 - t353 * mrSges(7,1) + mrSges(7,2) * t354 - t385 * t374 + t375 * t386;
t445 = -m(6) * t356 + t383 * mrSges(6,1) - mrSges(6,2) * t384 + t408 * t390 - t391 * t409 - t449;
t335 = m(5) * t364 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t413 + qJD(4) * t416 - t410 * t459 + t445;
t452 = t436 * t325 + t440 * t335;
t373 = -qJDD(2) * pkin(2) + t447;
t450 = -m(4) * t373 + (t443 * mrSges(4,3)) - t452;
t362 = Ifges(7,4) * t386 + Ifges(7,2) * t385 + Ifges(7,6) * t419;
t363 = Ifges(7,1) * t386 + Ifges(7,4) * t385 + Ifges(7,5) * t419;
t448 = -mrSges(7,1) * t342 + mrSges(7,2) * t343 - Ifges(7,5) * t354 - Ifges(7,6) * t353 - Ifges(7,3) * t400 - t386 * t362 + t385 * t363;
t372 = t443 * pkin(2) + t463;
t446 = -m(4) * t372 + m(5) * t370 - mrSges(5,1) * t412 + (t443 * mrSges(4,2)) + t413 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t416 * t423 + t417 * t459 + t326;
t379 = Ifges(6,4) * t409 + Ifges(6,2) * t408 + Ifges(6,6) * t420;
t380 = Ifges(6,1) * t409 + Ifges(6,4) * t408 + Ifges(6,5) * t420;
t444 = mrSges(6,1) * t346 - mrSges(6,2) * t347 + Ifges(6,5) * t384 + Ifges(6,6) * t383 + Ifges(6,3) * t405 + pkin(5) * t332 + t409 * t379 - t408 * t380 - t448;
t399 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t440 - Ifges(5,4) * t436) * qJD(2);
t398 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t440 - Ifges(5,2) * t436) * qJD(2);
t378 = Ifges(6,5) * t409 + Ifges(6,6) * t408 + Ifges(6,3) * t420;
t361 = Ifges(7,5) * t386 + Ifges(7,6) * t385 + Ifges(7,3) * t419;
t334 = mrSges(7,2) * t348 - mrSges(7,3) * t342 + Ifges(7,1) * t354 + Ifges(7,4) * t353 + Ifges(7,5) * t400 + t361 * t385 - t362 * t419;
t333 = -mrSges(7,1) * t348 + mrSges(7,3) * t343 + Ifges(7,4) * t354 + Ifges(7,2) * t353 + Ifges(7,6) * t400 - t361 * t386 + t363 * t419;
t324 = mrSges(6,2) * t356 - mrSges(6,3) * t346 + Ifges(6,1) * t384 + Ifges(6,4) * t383 + Ifges(6,5) * t405 - pkin(10) * t332 - t333 * t434 + t334 * t438 + t378 * t408 - t379 * t420;
t323 = -mrSges(6,1) * t356 + mrSges(6,3) * t347 + Ifges(6,4) * t384 + Ifges(6,2) * t383 + Ifges(6,6) * t405 - pkin(5) * t449 + pkin(10) * t454 + t438 * t333 + t434 * t334 - t409 * t378 + t420 * t380;
t322 = qJDD(2) * mrSges(4,2) - t450;
t1 = [m(2) * t427 + t433 * (t440 * t325 - t436 * t335 + (m(3) + m(4)) * t393) + (t437 * (m(3) * t377 - (mrSges(3,1) * t443) - qJDD(2) * mrSges(3,2) + t446) + t441 * (m(3) * t376 - t443 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t450)) * t431; mrSges(3,1) * t376 - mrSges(3,2) * t377 + mrSges(4,2) * t373 - mrSges(4,3) * t372 + t440 * (mrSges(5,2) * t370 - mrSges(5,3) * t364 + Ifges(5,1) * t413 + Ifges(5,4) * t412 + Ifges(5,5) * qJDD(4) - pkin(9) * t326 - qJD(4) * t398 - t323 * t435 + t324 * t439) - t436 * (-mrSges(5,1) * t370 + mrSges(5,3) * t365 + Ifges(5,4) * t413 + Ifges(5,2) * t412 + Ifges(5,6) * qJDD(4) - pkin(4) * t326 + qJD(4) * t399 - t444) - pkin(8) * t452 - pkin(2) * t322 + qJ(3) * t446 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t322; Ifges(5,5) * t413 + Ifges(5,6) * t412 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t364 - mrSges(5,2) * t365 + t435 * t324 + t439 * t323 + pkin(4) * t445 + pkin(9) * t455 + (t398 * t440 + t399 * t436) * qJD(2); t444; -t448;];
tauJ  = t1;
