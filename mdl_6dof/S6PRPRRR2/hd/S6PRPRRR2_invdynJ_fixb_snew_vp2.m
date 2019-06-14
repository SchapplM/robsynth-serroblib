% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:22:11
% EndTime: 2019-05-05 00:22:13
% DurationCPUTime: 1.55s
% Computational Cost: add. (16893->236), mult. (31432->302), div. (0->0), fcn. (22366->14), ass. (0->108)
t419 = sin(pkin(11));
t422 = cos(pkin(11));
t407 = g(1) * t419 - g(2) * t422;
t423 = cos(pkin(6));
t451 = t407 * t423;
t420 = sin(pkin(6));
t427 = sin(qJ(2));
t450 = t420 * t427;
t431 = cos(qJ(2));
t449 = t420 * t431;
t408 = -g(1) * t422 - g(2) * t419;
t417 = -g(3) + qJDD(1);
t370 = -t408 * t427 + t417 * t449 + t431 * t451;
t368 = qJDD(2) * pkin(2) + t370;
t371 = t431 * t408 + t417 * t450 + t427 * t451;
t433 = qJD(2) ^ 2;
t369 = -pkin(2) * t433 + t371;
t418 = sin(pkin(12));
t421 = cos(pkin(12));
t357 = t418 * t368 + t421 * t369;
t355 = -pkin(3) * t433 + qJDD(2) * pkin(8) + t357;
t439 = -t407 * t420 + t423 * t417;
t387 = qJDD(3) + t439;
t426 = sin(qJ(4));
t430 = cos(qJ(4));
t345 = t430 * t355 + t426 * t387;
t403 = (-mrSges(5,1) * t430 + mrSges(5,2) * t426) * qJD(2);
t445 = qJD(2) * qJD(4);
t416 = t426 * t445;
t406 = qJDD(2) * t430 - t416;
t447 = qJD(2) * t426;
t409 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t447;
t404 = (-pkin(4) * t430 - pkin(9) * t426) * qJD(2);
t432 = qJD(4) ^ 2;
t446 = qJD(2) * t430;
t340 = -pkin(4) * t432 + qJDD(4) * pkin(9) + t404 * t446 + t345;
t356 = t421 * t368 - t418 * t369;
t354 = -qJDD(2) * pkin(3) - t433 * pkin(8) - t356;
t443 = t430 * t445;
t405 = qJDD(2) * t426 + t443;
t343 = (-t405 - t443) * pkin(9) + (-t406 + t416) * pkin(4) + t354;
t425 = sin(qJ(5));
t429 = cos(qJ(5));
t335 = -t425 * t340 + t429 * t343;
t401 = qJD(4) * t429 - t425 * t447;
t378 = qJD(5) * t401 + qJDD(4) * t425 + t405 * t429;
t399 = qJDD(5) - t406;
t402 = qJD(4) * t425 + t429 * t447;
t414 = qJD(5) - t446;
t333 = (t401 * t414 - t378) * pkin(10) + (t401 * t402 + t399) * pkin(5) + t335;
t336 = t429 * t340 + t425 * t343;
t377 = -qJD(5) * t402 + qJDD(4) * t429 - t405 * t425;
t386 = pkin(5) * t414 - pkin(10) * t402;
t398 = t401 ^ 2;
t334 = -pkin(5) * t398 + pkin(10) * t377 - t386 * t414 + t336;
t424 = sin(qJ(6));
t428 = cos(qJ(6));
t331 = t333 * t428 - t334 * t424;
t379 = t401 * t428 - t402 * t424;
t351 = qJD(6) * t379 + t377 * t424 + t378 * t428;
t380 = t401 * t424 + t402 * t428;
t362 = -mrSges(7,1) * t379 + mrSges(7,2) * t380;
t413 = qJD(6) + t414;
t366 = -mrSges(7,2) * t413 + mrSges(7,3) * t379;
t395 = qJDD(6) + t399;
t328 = m(7) * t331 + mrSges(7,1) * t395 - mrSges(7,3) * t351 - t362 * t380 + t366 * t413;
t332 = t333 * t424 + t334 * t428;
t350 = -qJD(6) * t380 + t377 * t428 - t378 * t424;
t367 = mrSges(7,1) * t413 - mrSges(7,3) * t380;
t329 = m(7) * t332 - mrSges(7,2) * t395 + mrSges(7,3) * t350 + t362 * t379 - t367 * t413;
t320 = t428 * t328 + t424 * t329;
t381 = -mrSges(6,1) * t401 + mrSges(6,2) * t402;
t384 = -mrSges(6,2) * t414 + mrSges(6,3) * t401;
t318 = m(6) * t335 + mrSges(6,1) * t399 - mrSges(6,3) * t378 - t381 * t402 + t384 * t414 + t320;
t385 = mrSges(6,1) * t414 - mrSges(6,3) * t402;
t440 = -t328 * t424 + t428 * t329;
t319 = m(6) * t336 - mrSges(6,2) * t399 + mrSges(6,3) * t377 + t381 * t401 - t385 * t414 + t440;
t441 = -t318 * t425 + t429 * t319;
t315 = m(5) * t345 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t406 - qJD(4) * t409 + t403 * t446 + t441;
t344 = -t426 * t355 + t387 * t430;
t410 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t446;
t339 = -qJDD(4) * pkin(4) - pkin(9) * t432 + t404 * t447 - t344;
t337 = -pkin(5) * t377 - pkin(10) * t398 + t386 * t402 + t339;
t438 = m(7) * t337 - t350 * mrSges(7,1) + mrSges(7,2) * t351 - t379 * t366 + t367 * t380;
t435 = -m(6) * t339 + t377 * mrSges(6,1) - mrSges(6,2) * t378 + t401 * t384 - t385 * t402 - t438;
t324 = m(5) * t344 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t405 + qJD(4) * t410 - t403 * t447 + t435;
t442 = t430 * t315 - t324 * t426;
t308 = m(4) * t357 - mrSges(4,1) * t433 - qJDD(2) * mrSges(4,2) + t442;
t316 = t318 * t429 + t319 * t425;
t436 = -m(5) * t354 + t406 * mrSges(5,1) - mrSges(5,2) * t405 - t409 * t447 + t410 * t446 - t316;
t312 = m(4) * t356 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t433 + t436;
t448 = t418 * t308 + t421 * t312;
t444 = m(4) * t387 + t426 * t315 + t430 * t324;
t359 = Ifges(7,4) * t380 + Ifges(7,2) * t379 + Ifges(7,6) * t413;
t360 = Ifges(7,1) * t380 + Ifges(7,4) * t379 + Ifges(7,5) * t413;
t437 = -mrSges(7,1) * t331 + mrSges(7,2) * t332 - Ifges(7,5) * t351 - Ifges(7,6) * t350 - Ifges(7,3) * t395 - t380 * t359 + t379 * t360;
t373 = Ifges(6,4) * t402 + Ifges(6,2) * t401 + Ifges(6,6) * t414;
t374 = Ifges(6,1) * t402 + Ifges(6,4) * t401 + Ifges(6,5) * t414;
t434 = mrSges(6,1) * t335 - mrSges(6,2) * t336 + Ifges(6,5) * t378 + Ifges(6,6) * t377 + Ifges(6,3) * t399 + pkin(5) * t320 + t402 * t373 - t401 * t374 - t437;
t394 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t426 + Ifges(5,4) * t430) * qJD(2);
t393 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t426 + Ifges(5,2) * t430) * qJD(2);
t372 = Ifges(6,5) * t402 + Ifges(6,6) * t401 + Ifges(6,3) * t414;
t358 = Ifges(7,5) * t380 + Ifges(7,6) * t379 + Ifges(7,3) * t413;
t322 = mrSges(7,2) * t337 - mrSges(7,3) * t331 + Ifges(7,1) * t351 + Ifges(7,4) * t350 + Ifges(7,5) * t395 + t358 * t379 - t359 * t413;
t321 = -mrSges(7,1) * t337 + mrSges(7,3) * t332 + Ifges(7,4) * t351 + Ifges(7,2) * t350 + Ifges(7,6) * t395 - t358 * t380 + t360 * t413;
t310 = mrSges(6,2) * t339 - mrSges(6,3) * t335 + Ifges(6,1) * t378 + Ifges(6,4) * t377 + Ifges(6,5) * t399 - pkin(10) * t320 - t321 * t424 + t322 * t428 + t372 * t401 - t373 * t414;
t309 = -mrSges(6,1) * t339 + mrSges(6,3) * t336 + Ifges(6,4) * t378 + Ifges(6,2) * t377 + Ifges(6,6) * t399 - pkin(5) * t438 + pkin(10) * t440 + t428 * t321 + t424 * t322 - t402 * t372 + t414 * t374;
t1 = [m(2) * t417 + (m(3) * t371 - mrSges(3,1) * t433 - qJDD(2) * mrSges(3,2) + t308 * t421 - t312 * t418) * t450 + (m(3) * t370 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t433 + t448) * t449 + t423 * (m(3) * t439 + t444); mrSges(3,1) * t370 - mrSges(3,2) * t371 + mrSges(4,1) * t356 - mrSges(4,2) * t357 + t426 * (mrSges(5,2) * t354 - mrSges(5,3) * t344 + Ifges(5,1) * t405 + Ifges(5,4) * t406 + Ifges(5,5) * qJDD(4) - pkin(9) * t316 - qJD(4) * t393 - t309 * t425 + t310 * t429) + t430 * (-mrSges(5,1) * t354 + mrSges(5,3) * t345 + Ifges(5,4) * t405 + Ifges(5,2) * t406 + Ifges(5,6) * qJDD(4) - pkin(4) * t316 + qJD(4) * t394 - t434) + pkin(3) * t436 + pkin(8) * t442 + pkin(2) * t448 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t444; Ifges(5,5) * t405 + Ifges(5,6) * t406 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t344 - mrSges(5,2) * t345 + t425 * t310 + t429 * t309 + pkin(4) * t435 + pkin(9) * t441 + (t393 * t426 - t394 * t430) * qJD(2); t434; -t437;];
tauJ  = t1;
