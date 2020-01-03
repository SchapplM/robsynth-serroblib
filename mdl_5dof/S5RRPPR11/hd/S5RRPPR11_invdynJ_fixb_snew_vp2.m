% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:39
% EndTime: 2019-12-31 19:46:42
% DurationCPUTime: 1.75s
% Computational Cost: add. (8990->260), mult. (19649->318), div. (0->0), fcn. (11188->8), ass. (0->104)
t454 = -2 * qJD(3);
t453 = Ifges(3,1) + Ifges(4,2);
t448 = Ifges(3,4) + Ifges(4,6);
t447 = Ifges(3,5) - Ifges(4,4);
t452 = Ifges(3,2) + Ifges(4,3);
t446 = Ifges(3,6) - Ifges(4,5);
t451 = (Ifges(3,3) + Ifges(4,1));
t423 = qJD(1) ^ 2;
t418 = sin(qJ(1));
t421 = cos(qJ(1));
t432 = -g(1) * t421 - g(2) * t418;
t385 = -pkin(1) * t423 + qJDD(1) * pkin(6) + t432;
t417 = sin(qJ(2));
t420 = cos(qJ(2));
t364 = -g(3) * t417 + t420 * t385;
t392 = (-pkin(2) * t420 - qJ(3) * t417) * qJD(1);
t422 = qJD(2) ^ 2;
t440 = qJD(1) * t420;
t352 = pkin(2) * t422 - qJDD(2) * qJ(3) + (qJD(2) * t454) - t392 * t440 - t364;
t450 = pkin(6) * t423;
t449 = mrSges(3,1) - mrSges(4,2);
t439 = qJD(1) * qJD(2);
t436 = t417 * t439;
t396 = qJDD(1) * t420 - t436;
t441 = qJD(1) * t417;
t401 = pkin(3) * t441 - (qJD(2) * qJ(4));
t413 = t420 ^ 2;
t437 = t420 * t439;
t395 = qJDD(1) * t417 + t437;
t435 = g(1) * t418 - t421 * g(2);
t430 = -qJDD(1) * pkin(1) - t435;
t425 = pkin(2) * t436 + t441 * t454 + (-t395 - t437) * qJ(3) + t430;
t335 = -t401 * t441 + (-pkin(3) * t413 - pkin(6)) * t423 + (-pkin(2) - qJ(4)) * t396 + t425;
t363 = -t420 * g(3) - t417 * t385;
t353 = -qJDD(2) * pkin(2) - qJ(3) * t422 + t392 * t441 + qJDD(3) - t363;
t348 = (-t417 * t420 * t423 - qJDD(2)) * qJ(4) + (t395 - t437) * pkin(3) + t353;
t414 = sin(pkin(8));
t415 = cos(pkin(8));
t389 = qJD(2) * t415 - t414 * t440;
t329 = -0.2e1 * qJD(4) * t389 - t335 * t414 + t348 * t415;
t369 = qJDD(2) * t415 - t396 * t414;
t388 = -qJD(2) * t414 - t415 * t440;
t327 = (t388 * t441 - t369) * pkin(7) + (t388 * t389 + t395) * pkin(4) + t329;
t330 = 0.2e1 * qJD(4) * t388 + t335 * t415 + t348 * t414;
t368 = -qJDD(2) * t414 - t396 * t415;
t370 = pkin(4) * t441 - pkin(7) * t389;
t387 = t388 ^ 2;
t328 = -pkin(4) * t387 + pkin(7) * t368 - t370 * t441 + t330;
t416 = sin(qJ(5));
t419 = cos(qJ(5));
t325 = t327 * t419 - t328 * t416;
t360 = t388 * t419 - t389 * t416;
t340 = qJD(5) * t360 + t368 * t416 + t369 * t419;
t361 = t388 * t416 + t389 * t419;
t350 = -mrSges(6,1) * t360 + mrSges(6,2) * t361;
t405 = qJD(5) + t441;
t354 = -mrSges(6,2) * t405 + mrSges(6,3) * t360;
t391 = qJDD(5) + t395;
t321 = m(6) * t325 + mrSges(6,1) * t391 - mrSges(6,3) * t340 - t350 * t361 + t354 * t405;
t326 = t327 * t416 + t328 * t419;
t339 = -qJD(5) * t361 + t368 * t419 - t369 * t416;
t355 = mrSges(6,1) * t405 - mrSges(6,3) * t361;
t322 = m(6) * t326 - mrSges(6,2) * t391 + mrSges(6,3) * t339 + t350 * t360 - t355 * t405;
t315 = t321 * t419 + t322 * t416;
t445 = (t451 * qJD(2)) + (t417 * t447 + t420 * t446) * qJD(1);
t444 = t446 * qJD(2) + (t417 * t448 + t420 * t452) * qJD(1);
t443 = t447 * qJD(2) + (t417 * t453 + t420 * t448) * qJD(1);
t402 = -mrSges(4,1) * t440 - (qJD(2) * mrSges(4,3));
t442 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t440 - t402;
t362 = -mrSges(5,1) * t388 + mrSges(5,2) * t389;
t366 = -mrSges(5,2) * t441 + mrSges(5,3) * t388;
t313 = m(5) * t329 + mrSges(5,1) * t395 - mrSges(5,3) * t369 - t362 * t389 + t366 * t441 + t315;
t367 = mrSges(5,1) * t441 - mrSges(5,3) * t389;
t433 = -t321 * t416 + t322 * t419;
t314 = m(5) * t330 - mrSges(5,2) * t395 + mrSges(5,3) * t368 + t362 * t388 - t367 * t441 + t433;
t434 = -t313 * t414 + t314 * t415;
t311 = t415 * t313 + t414 * t314;
t351 = -pkin(2) * t396 + t425 - t450;
t431 = m(4) * t351 + t434;
t344 = -qJ(4) * t413 * t423 + pkin(3) * t396 + qJD(2) * t401 + qJDD(4) - t352;
t332 = -pkin(4) * t368 - pkin(7) * t387 + t370 * t389 + t344;
t428 = m(6) * t332 - t339 * mrSges(6,1) + mrSges(6,2) * t340 - t360 * t354 + t361 * t355;
t427 = m(4) * t353 + t395 * mrSges(4,1) + t311;
t346 = Ifges(6,4) * t361 + Ifges(6,2) * t360 + Ifges(6,6) * t405;
t347 = Ifges(6,1) * t361 + Ifges(6,4) * t360 + Ifges(6,5) * t405;
t426 = mrSges(6,1) * t325 - mrSges(6,2) * t326 + Ifges(6,5) * t340 + Ifges(6,6) * t339 + Ifges(6,3) * t391 + t346 * t361 - t360 * t347;
t323 = m(5) * t344 - t368 * mrSges(5,1) + t369 * mrSges(5,2) - t388 * t366 + t389 * t367 + t428;
t393 = (mrSges(4,2) * t420 - mrSges(4,3) * t417) * qJD(1);
t403 = mrSges(4,1) * t441 + qJD(2) * mrSges(4,2);
t424 = -m(4) * t352 + qJDD(2) * mrSges(4,3) + qJD(2) * t403 + t393 * t440 + t323;
t399 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t441;
t394 = (-mrSges(3,1) * t420 + mrSges(3,2) * t417) * qJD(1);
t384 = t430 - t450;
t358 = Ifges(5,1) * t389 + Ifges(5,4) * t388 + Ifges(5,5) * t441;
t357 = Ifges(5,4) * t389 + Ifges(5,2) * t388 + Ifges(5,6) * t441;
t356 = Ifges(5,5) * t389 + Ifges(5,6) * t388 + Ifges(5,3) * t441;
t345 = Ifges(6,5) * t361 + Ifges(6,6) * t360 + Ifges(6,3) * t405;
t317 = mrSges(6,2) * t332 - mrSges(6,3) * t325 + Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t391 + t345 * t360 - t346 * t405;
t316 = -mrSges(6,1) * t332 + mrSges(6,3) * t326 + Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t391 - t345 * t361 + t347 * t405;
t310 = qJDD(2) * mrSges(4,2) + qJD(2) * t402 + t393 * t441 + t427;
t309 = mrSges(4,2) * t396 - mrSges(4,3) * t395 + (t402 * t420 - t403 * t417) * qJD(1) + t431;
t308 = mrSges(5,2) * t344 - mrSges(5,3) * t329 + Ifges(5,1) * t369 + Ifges(5,4) * t368 + Ifges(5,5) * t395 - pkin(7) * t315 - t316 * t416 + t317 * t419 + t356 * t388 - t357 * t441;
t307 = -mrSges(5,1) * t344 + mrSges(5,3) * t330 + Ifges(5,4) * t369 + Ifges(5,2) * t368 + Ifges(5,6) * t395 - pkin(4) * t428 + pkin(7) * t433 + t419 * t316 + t416 * t317 - t389 * t356 + t358 * t441;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t435 - mrSges(2,2) * t432 + t417 * (t426 + t445 * t440 + t448 * t396 + (Ifges(5,3) + t453) * t395 - t444 * qJD(2) - t388 * t358 + t389 * t357 + mrSges(3,2) * t384 + Ifges(5,6) * t368 + Ifges(5,5) * t369 - mrSges(3,3) * t363 - mrSges(4,3) * t351 + mrSges(4,1) * t353 + mrSges(5,1) * t329 - mrSges(5,2) * t330 + pkin(4) * t315 + pkin(3) * t311 - qJ(3) * t309 + t447 * qJDD(2)) + t420 * (-mrSges(3,1) * t384 - mrSges(4,1) * t352 + mrSges(4,2) * t351 + mrSges(3,3) * t364 - pkin(2) * t309 + pkin(3) * t323 - qJ(4) * t434 + t443 * qJD(2) + t446 * qJDD(2) - t415 * t307 - t414 * t308 + t448 * t395 + t452 * t396 - t445 * t441) + pkin(1) * (-m(3) * t384 + t449 * t396 + (-mrSges(3,2) + mrSges(4,3)) * t395 + (t442 * t420 + (-t399 + t403) * t417) * qJD(1) - t431) + pkin(6) * (t420 * (t424 - qJDD(2) * mrSges(3,2) + t394 * t440 + (mrSges(3,3) + mrSges(4,1)) * t396 - qJD(2) * t399 + m(3) * t364) + (-m(3) * t363 + t395 * mrSges(3,3) - t449 * qJDD(2) - t442 * qJD(2) + (t393 + t394) * t441 + t427) * t417); mrSges(3,1) * t363 - mrSges(3,2) * t364 + mrSges(4,2) * t353 - mrSges(4,3) * t352 + t415 * t308 - t414 * t307 - qJ(4) * t311 - pkin(2) * t310 + qJ(3) * t424 + (mrSges(4,1) * qJ(3) + t446) * t396 + t447 * t395 + t451 * qJDD(2) + (t444 * t417 - t443 * t420) * qJD(1); t310; t323; t426;];
tauJ = t1;
