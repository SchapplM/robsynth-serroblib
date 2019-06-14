% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-04 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:08:09
% EndTime: 2019-05-04 23:08:11
% DurationCPUTime: 1.40s
% Computational Cost: add. (10483->232), mult. (20778->294), div. (0->0), fcn. (13708->12), ass. (0->99)
t422 = sin(pkin(10));
t425 = cos(pkin(10));
t409 = t422 * g(1) - t425 * g(2);
t418 = -g(3) + qJDD(1);
t423 = sin(pkin(6));
t426 = cos(pkin(6));
t455 = t409 * t426 + t418 * t423;
t410 = -t425 * g(1) - t422 * g(2);
t429 = sin(qJ(2));
t432 = cos(qJ(2));
t372 = t432 * t410 + t455 * t429;
t454 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t372;
t371 = -t429 * t410 + t455 * t432;
t453 = -pkin(2) - pkin(8);
t434 = qJD(2) ^ 2;
t437 = -t434 * qJ(3) + qJDD(3) - t371;
t366 = t453 * qJDD(2) + t437;
t388 = -t423 * t409 + t426 * t418;
t428 = sin(qJ(4));
t431 = cos(qJ(4));
t361 = t428 * t366 + t431 * t388;
t405 = (t428 * pkin(4) - t431 * qJ(5)) * qJD(2);
t433 = qJD(4) ^ 2;
t450 = t428 * qJD(2);
t349 = -t433 * pkin(4) + qJDD(4) * qJ(5) - t405 * t450 + t361;
t365 = t453 * t434 - t454;
t448 = qJD(2) * qJD(4);
t445 = t431 * t448;
t407 = t428 * qJDD(2) + t445;
t446 = t428 * t448;
t408 = t431 * qJDD(2) - t446;
t352 = (-t408 + t446) * qJ(5) + (t407 + t445) * pkin(4) + t365;
t421 = sin(pkin(11));
t424 = cos(pkin(11));
t449 = t431 * qJD(2);
t400 = t421 * qJD(4) + t424 * t449;
t344 = -0.2e1 * qJD(5) * t400 - t421 * t349 + t424 * t352;
t386 = t421 * qJDD(4) + t424 * t408;
t399 = t424 * qJD(4) - t421 * t449;
t342 = (t399 * t450 - t386) * pkin(9) + (t399 * t400 + t407) * pkin(5) + t344;
t345 = 0.2e1 * qJD(5) * t399 + t424 * t349 + t421 * t352;
t385 = t424 * qJDD(4) - t421 * t408;
t387 = pkin(5) * t450 - t400 * pkin(9);
t398 = t399 ^ 2;
t343 = -t398 * pkin(5) + t385 * pkin(9) - t387 * t450 + t345;
t427 = sin(qJ(6));
t430 = cos(qJ(6));
t340 = t430 * t342 - t427 * t343;
t377 = t430 * t399 - t427 * t400;
t355 = t377 * qJD(6) + t427 * t385 + t430 * t386;
t378 = t427 * t399 + t430 * t400;
t362 = -t377 * mrSges(7,1) + t378 * mrSges(7,2);
t414 = qJD(6) + t450;
t369 = -t414 * mrSges(7,2) + t377 * mrSges(7,3);
t402 = qJDD(6) + t407;
t337 = m(7) * t340 + t402 * mrSges(7,1) - t355 * mrSges(7,3) - t378 * t362 + t414 * t369;
t341 = t427 * t342 + t430 * t343;
t354 = -t378 * qJD(6) + t430 * t385 - t427 * t386;
t370 = t414 * mrSges(7,1) - t378 * mrSges(7,3);
t338 = m(7) * t341 - t402 * mrSges(7,2) + t354 * mrSges(7,3) + t377 * t362 - t414 * t370;
t330 = t430 * t337 + t427 * t338;
t379 = -t399 * mrSges(6,1) + t400 * mrSges(6,2);
t383 = -mrSges(6,2) * t450 + t399 * mrSges(6,3);
t328 = m(6) * t344 + t407 * mrSges(6,1) - t386 * mrSges(6,3) - t400 * t379 + t383 * t450 + t330;
t384 = mrSges(6,1) * t450 - t400 * mrSges(6,3);
t443 = -t427 * t337 + t430 * t338;
t329 = m(6) * t345 - t407 * mrSges(6,2) + t385 * mrSges(6,3) + t399 * t379 - t384 * t450 + t443;
t324 = t424 * t328 + t421 * t329;
t444 = -t421 * t328 + t424 * t329;
t360 = t431 * t366 - t428 * t388;
t406 = (t428 * mrSges(5,1) + t431 * mrSges(5,2)) * qJD(2);
t412 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t449;
t323 = m(5) * t361 - qJDD(4) * mrSges(5,2) - t407 * mrSges(5,3) - qJD(4) * t412 - t406 * t450 + t444;
t348 = -qJDD(4) * pkin(4) - t433 * qJ(5) + t405 * t449 + qJDD(5) - t360;
t346 = -t385 * pkin(5) - t398 * pkin(9) + t400 * t387 + t348;
t438 = m(7) * t346 - t354 * mrSges(7,1) + t355 * mrSges(7,2) - t377 * t369 + t378 * t370;
t339 = m(6) * t348 - t385 * mrSges(6,1) + t386 * mrSges(6,2) - t399 * t383 + t400 * t384 + t438;
t411 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t450;
t333 = m(5) * t360 + qJDD(4) * mrSges(5,1) - t408 * mrSges(5,3) + qJD(4) * t411 - t406 * t449 - t339;
t441 = t428 * t323 + t431 * t333;
t368 = -qJDD(2) * pkin(2) + t437;
t439 = -m(4) * t368 + (t434 * mrSges(4,3)) - t441;
t367 = t434 * pkin(2) + t454;
t436 = -m(4) * t367 + m(5) * t365 + t407 * mrSges(5,1) + (t434 * mrSges(4,2)) + t408 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t411 * t450 + t412 * t449 + t324;
t357 = Ifges(7,4) * t378 + Ifges(7,2) * t377 + Ifges(7,6) * t414;
t358 = Ifges(7,1) * t378 + Ifges(7,4) * t377 + Ifges(7,5) * t414;
t435 = mrSges(7,1) * t340 - mrSges(7,2) * t341 + Ifges(7,5) * t355 + Ifges(7,6) * t354 + Ifges(7,3) * t402 + t378 * t357 - t377 * t358;
t395 = (Ifges(5,5) * qJD(4)) + (t431 * Ifges(5,1) - t428 * Ifges(5,4)) * qJD(2);
t394 = (Ifges(5,6) * qJD(4)) + (t431 * Ifges(5,4) - Ifges(5,2) * t428) * qJD(2);
t375 = Ifges(6,1) * t400 + Ifges(6,4) * t399 + Ifges(6,5) * t450;
t374 = Ifges(6,4) * t400 + Ifges(6,2) * t399 + Ifges(6,6) * t450;
t373 = Ifges(6,5) * t400 + Ifges(6,6) * t399 + Ifges(6,3) * t450;
t356 = Ifges(7,5) * t378 + Ifges(7,6) * t377 + Ifges(7,3) * t414;
t332 = mrSges(7,2) * t346 - mrSges(7,3) * t340 + Ifges(7,1) * t355 + Ifges(7,4) * t354 + Ifges(7,5) * t402 + t377 * t356 - t414 * t357;
t331 = -mrSges(7,1) * t346 + mrSges(7,3) * t341 + Ifges(7,4) * t355 + Ifges(7,2) * t354 + Ifges(7,6) * t402 - t378 * t356 + t414 * t358;
t322 = mrSges(6,2) * t348 - mrSges(6,3) * t344 + Ifges(6,1) * t386 + Ifges(6,4) * t385 + Ifges(6,5) * t407 - pkin(9) * t330 - t427 * t331 + t430 * t332 + t399 * t373 - t374 * t450;
t321 = -mrSges(6,1) * t348 + mrSges(6,3) * t345 + Ifges(6,4) * t386 + Ifges(6,2) * t385 + Ifges(6,6) * t407 - pkin(5) * t438 + pkin(9) * t443 + t430 * t331 + t427 * t332 - t400 * t373 + t375 * t450;
t320 = qJDD(2) * mrSges(4,2) - t439;
t1 = [m(2) * t418 + t426 * (t431 * t323 - t428 * t333 + (m(3) + m(4)) * t388) + (t429 * (m(3) * t372 - (t434 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t436) + t432 * (m(3) * t371 - t434 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t439)) * t423; mrSges(3,1) * t371 - mrSges(3,2) * t372 + mrSges(4,2) * t368 - mrSges(4,3) * t367 + t431 * (mrSges(5,2) * t365 - mrSges(5,3) * t360 + Ifges(5,1) * t408 - Ifges(5,4) * t407 + Ifges(5,5) * qJDD(4) - qJ(5) * t324 - qJD(4) * t394 - t421 * t321 + t424 * t322) - t428 * (-mrSges(5,1) * t365 - mrSges(6,1) * t344 + mrSges(6,2) * t345 + mrSges(5,3) * t361 + Ifges(5,4) * t408 - Ifges(6,5) * t386 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t385 - pkin(4) * t324 - pkin(5) * t330 + qJD(4) * t395 - t400 * t374 + t399 * t375 - t435 + (-Ifges(5,2) - Ifges(6,3)) * t407) - pkin(8) * t441 - pkin(2) * t320 + qJ(3) * t436 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t320; Ifges(5,5) * t408 - Ifges(5,6) * t407 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t360 - mrSges(5,2) * t361 + t421 * t322 + t424 * t321 - pkin(4) * t339 + qJ(5) * t444 + (t431 * t394 + t428 * t395) * qJD(2); t339; t435;];
tauJ  = t1;
