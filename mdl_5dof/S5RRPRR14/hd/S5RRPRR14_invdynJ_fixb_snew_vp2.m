% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:36:35
% EndTime: 2019-12-31 20:36:40
% DurationCPUTime: 4.52s
% Computational Cost: add. (47765->287), mult. (109063->381), div. (0->0), fcn. (85494->12), ass. (0->123)
t427 = sin(pkin(5));
t456 = pkin(7) * t427;
t429 = cos(pkin(5));
t455 = t429 * g(3);
t438 = qJD(1) ^ 2;
t433 = sin(qJ(1));
t437 = cos(qJ(1));
t445 = t433 * g(1) - t437 * g(2);
t412 = qJDD(1) * pkin(1) + t438 * t456 + t445;
t454 = t412 * t429;
t432 = sin(qJ(2));
t453 = t427 * t432;
t436 = cos(qJ(2));
t452 = t427 * t436;
t450 = qJD(1) * t427;
t414 = (-t436 * pkin(2) - t432 * qJ(3)) * t450;
t423 = t429 * qJD(1) + qJD(2);
t421 = t423 ^ 2;
t422 = t429 * qJDD(1) + qJDD(2);
t449 = qJD(1) * t436;
t442 = -t437 * g(1) - t433 * g(2);
t448 = qJDD(1) * t427;
t413 = -t438 * pkin(1) + pkin(7) * t448 + t442;
t451 = t436 * t413 + t432 * t454;
t373 = -t421 * pkin(2) + t422 * qJ(3) + (-g(3) * t432 + t414 * t449) * t427 + t451;
t416 = (qJD(2) * t449 + qJDD(1) * t432) * t427;
t447 = t432 * t450;
t417 = -qJD(2) * t447 + t436 * t448;
t374 = -t417 * pkin(2) - t455 - t416 * qJ(3) + (-t412 + (pkin(2) * t432 - qJ(3) * t436) * t423 * qJD(1)) * t427;
t426 = sin(pkin(10));
t428 = cos(pkin(10));
t406 = t426 * t423 + t428 * t447;
t346 = -0.2e1 * qJD(3) * t406 - t426 * t373 + t428 * t374;
t394 = t428 * t416 + t426 * t422;
t405 = t428 * t423 - t426 * t447;
t446 = t427 * t449;
t343 = (-t405 * t446 - t394) * pkin(8) + (t405 * t406 - t417) * pkin(3) + t346;
t347 = 0.2e1 * qJD(3) * t405 + t428 * t373 + t426 * t374;
t393 = -t426 * t416 + t428 * t422;
t395 = -pkin(3) * t446 - t406 * pkin(8);
t404 = t405 ^ 2;
t345 = -t404 * pkin(3) + t393 * pkin(8) + t395 * t446 + t347;
t431 = sin(qJ(4));
t435 = cos(qJ(4));
t340 = t431 * t343 + t435 * t345;
t387 = t435 * t405 - t431 * t406;
t388 = t431 * t405 + t435 * t406;
t368 = -t387 * pkin(4) - t388 * pkin(9);
t409 = qJDD(4) - t417;
t419 = qJD(4) - t446;
t418 = t419 ^ 2;
t338 = -t418 * pkin(4) + t409 * pkin(9) + t387 * t368 + t340;
t384 = -g(3) * t452 - t432 * t413 + t436 * t454;
t372 = -t422 * pkin(2) - t421 * qJ(3) + t414 * t447 + qJDD(3) - t384;
t351 = -t393 * pkin(3) - t404 * pkin(8) + t406 * t395 + t372;
t359 = -t388 * qJD(4) + t435 * t393 - t431 * t394;
t360 = t387 * qJD(4) + t431 * t393 + t435 * t394;
t341 = (-t387 * t419 - t360) * pkin(9) + (t388 * t419 - t359) * pkin(4) + t351;
t430 = sin(qJ(5));
t434 = cos(qJ(5));
t335 = -t430 * t338 + t434 * t341;
t375 = -t430 * t388 + t434 * t419;
t350 = t375 * qJD(5) + t434 * t360 + t430 * t409;
t376 = t434 * t388 + t430 * t419;
t356 = -t375 * mrSges(6,1) + t376 * mrSges(6,2);
t358 = qJDD(5) - t359;
t386 = qJD(5) - t387;
t361 = -t386 * mrSges(6,2) + t375 * mrSges(6,3);
t332 = m(6) * t335 + t358 * mrSges(6,1) - t350 * mrSges(6,3) - t376 * t356 + t386 * t361;
t336 = t434 * t338 + t430 * t341;
t349 = -t376 * qJD(5) - t430 * t360 + t434 * t409;
t362 = t386 * mrSges(6,1) - t376 * mrSges(6,3);
t333 = m(6) * t336 - t358 * mrSges(6,2) + t349 * mrSges(6,3) + t375 * t356 - t386 * t362;
t324 = -t430 * t332 + t434 * t333;
t367 = -t387 * mrSges(5,1) + t388 * mrSges(5,2);
t378 = t419 * mrSges(5,1) - t388 * mrSges(5,3);
t321 = m(5) * t340 - t409 * mrSges(5,2) + t359 * mrSges(5,3) + t387 * t367 - t419 * t378 + t324;
t339 = t435 * t343 - t431 * t345;
t337 = -t409 * pkin(4) - t418 * pkin(9) + t388 * t368 - t339;
t334 = -m(6) * t337 + t349 * mrSges(6,1) - t350 * mrSges(6,2) + t375 * t361 - t376 * t362;
t377 = -t419 * mrSges(5,2) + t387 * mrSges(5,3);
t328 = m(5) * t339 + t409 * mrSges(5,1) - t360 * mrSges(5,3) - t388 * t367 + t419 * t377 + t334;
t317 = t431 * t321 + t435 * t328;
t389 = -t405 * mrSges(4,1) + t406 * mrSges(4,2);
t391 = mrSges(4,2) * t446 + t405 * mrSges(4,3);
t315 = m(4) * t346 - t417 * mrSges(4,1) - t394 * mrSges(4,3) - t406 * t389 - t391 * t446 + t317;
t392 = -mrSges(4,1) * t446 - t406 * mrSges(4,3);
t443 = t435 * t321 - t431 * t328;
t316 = m(4) * t347 + t417 * mrSges(4,2) + t393 * mrSges(4,3) + t405 * t389 + t392 * t446 + t443;
t309 = t428 * t315 + t426 * t316;
t323 = t434 * t332 + t430 * t333;
t444 = -t426 * t315 + t428 * t316;
t441 = m(5) * t351 - t359 * mrSges(5,1) + t360 * mrSges(5,2) - t387 * t377 + t388 * t378 + t323;
t353 = Ifges(6,4) * t376 + Ifges(6,2) * t375 + Ifges(6,6) * t386;
t354 = Ifges(6,1) * t376 + Ifges(6,4) * t375 + Ifges(6,5) * t386;
t440 = mrSges(6,1) * t335 - mrSges(6,2) * t336 + Ifges(6,5) * t350 + Ifges(6,6) * t349 + Ifges(6,3) * t358 + t376 * t353 - t375 * t354;
t322 = m(4) * t372 - t393 * mrSges(4,1) + t394 * mrSges(4,2) - t405 * t391 + t406 * t392 + t441;
t352 = Ifges(6,5) * t376 + Ifges(6,6) * t375 + Ifges(6,3) * t386;
t325 = -mrSges(6,1) * t337 + mrSges(6,3) * t336 + Ifges(6,4) * t350 + Ifges(6,2) * t349 + Ifges(6,6) * t358 - t376 * t352 + t386 * t354;
t326 = mrSges(6,2) * t337 - mrSges(6,3) * t335 + Ifges(6,1) * t350 + Ifges(6,4) * t349 + Ifges(6,5) * t358 + t375 * t352 - t386 * t353;
t364 = Ifges(5,4) * t388 + Ifges(5,2) * t387 + Ifges(5,6) * t419;
t365 = Ifges(5,1) * t388 + Ifges(5,4) * t387 + Ifges(5,5) * t419;
t439 = mrSges(5,1) * t339 - mrSges(5,2) * t340 + Ifges(5,5) * t360 + Ifges(5,6) * t359 + Ifges(5,3) * t409 + pkin(4) * t334 + pkin(9) * t324 + t434 * t325 + t430 * t326 + t388 * t364 - t387 * t365;
t415 = (-t436 * mrSges(3,1) + t432 * mrSges(3,2)) * t450;
t411 = -t423 * mrSges(3,2) + mrSges(3,3) * t446;
t410 = t423 * mrSges(3,1) - mrSges(3,3) * t447;
t399 = -t427 * t412 - t455;
t398 = Ifges(3,5) * t423 + (t432 * Ifges(3,1) + t436 * Ifges(3,4)) * t450;
t397 = Ifges(3,6) * t423 + (t432 * Ifges(3,4) + Ifges(3,2) * t436) * t450;
t396 = Ifges(3,3) * t423 + (t432 * Ifges(3,5) + t436 * Ifges(3,6)) * t450;
t385 = -g(3) * t453 + t451;
t381 = Ifges(4,1) * t406 + Ifges(4,4) * t405 - Ifges(4,5) * t446;
t380 = Ifges(4,4) * t406 + Ifges(4,2) * t405 - Ifges(4,6) * t446;
t379 = Ifges(4,5) * t406 + Ifges(4,6) * t405 - Ifges(4,3) * t446;
t363 = Ifges(5,5) * t388 + Ifges(5,6) * t387 + Ifges(5,3) * t419;
t318 = m(3) * t384 + t422 * mrSges(3,1) - t416 * mrSges(3,3) + t423 * t411 - t415 * t447 - t322;
t311 = -mrSges(5,1) * t351 + mrSges(5,3) * t340 + Ifges(5,4) * t360 + Ifges(5,2) * t359 + Ifges(5,6) * t409 - pkin(4) * t323 - t388 * t363 + t419 * t365 - t440;
t310 = mrSges(5,2) * t351 - mrSges(5,3) * t339 + Ifges(5,1) * t360 + Ifges(5,4) * t359 + Ifges(5,5) * t409 - pkin(9) * t323 - t430 * t325 + t434 * t326 + t387 * t363 - t419 * t364;
t308 = m(3) * t385 - t422 * mrSges(3,2) + t417 * mrSges(3,3) - t423 * t410 + t415 * t446 + t444;
t307 = mrSges(4,2) * t372 - mrSges(4,3) * t346 + Ifges(4,1) * t394 + Ifges(4,4) * t393 - Ifges(4,5) * t417 - pkin(8) * t317 + t435 * t310 - t431 * t311 + t405 * t379 + t380 * t446;
t306 = -mrSges(4,1) * t372 + mrSges(4,3) * t347 + Ifges(4,4) * t394 + Ifges(4,2) * t393 - Ifges(4,6) * t417 - pkin(3) * t441 + pkin(8) * t443 + t431 * t310 + t435 * t311 - t406 * t379 - t381 * t446;
t305 = Ifges(3,5) * t416 + Ifges(3,6) * t417 + Ifges(3,3) * t422 + mrSges(3,1) * t384 - mrSges(3,2) * t385 + t426 * t307 + t428 * t306 - pkin(2) * t322 + qJ(3) * t444 + (t432 * t397 - t436 * t398) * t450;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t445 - mrSges(2,2) * t442 + (mrSges(3,2) * t399 - mrSges(3,3) * t384 + Ifges(3,1) * t416 + Ifges(3,4) * t417 + Ifges(3,5) * t422 - qJ(3) * t309 - t426 * t306 + t428 * t307 + t396 * t446 - t423 * t397) * t453 + (-t439 + (Ifges(3,2) + Ifges(4,3)) * t417 - t396 * t447 + Ifges(3,6) * t422 + t423 * t398 + Ifges(3,4) * t416 - t406 * t380 - Ifges(4,6) * t393 - Ifges(4,5) * t394 - mrSges(3,1) * t399 + t405 * t381 + mrSges(3,3) * t385 - mrSges(4,1) * t346 + mrSges(4,2) * t347 - pkin(3) * t317 - pkin(2) * t309) * t452 + t429 * t305 + pkin(1) * ((t432 * t308 + t436 * t318) * t429 + (-m(3) * t399 + t417 * mrSges(3,1) - t416 * mrSges(3,2) + (-t410 * t432 + t411 * t436) * t450 - t309) * t427) + (t436 * t308 - t432 * t318) * t456; t305; t322; t439; t440;];
tauJ = t1;
