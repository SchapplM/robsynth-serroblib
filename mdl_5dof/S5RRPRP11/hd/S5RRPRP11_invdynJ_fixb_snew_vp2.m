% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:32
% EndTime: 2019-12-31 20:12:34
% DurationCPUTime: 1.57s
% Computational Cost: add. (4187->237), mult. (8467->275), div. (0->0), fcn. (4393->6), ass. (0->97)
t459 = Ifges(5,1) + Ifges(6,1);
t444 = Ifges(5,4) - Ifges(6,5);
t455 = Ifges(6,4) + Ifges(5,5);
t458 = Ifges(5,2) + Ifges(6,3);
t453 = Ifges(5,6) - Ifges(6,6);
t457 = -2 * qJD(3);
t456 = Ifges(3,1) + Ifges(4,2);
t445 = Ifges(3,4) + Ifges(4,6);
t443 = Ifges(3,5) - Ifges(4,4);
t454 = Ifges(3,2) + Ifges(4,3);
t442 = Ifges(3,6) - Ifges(4,5);
t452 = Ifges(3,3) + Ifges(4,1);
t451 = Ifges(5,3) + Ifges(6,2);
t409 = sin(qJ(4));
t412 = cos(qJ(4));
t413 = cos(qJ(2));
t434 = qJD(1) * t413;
t384 = qJD(2) * t409 + t412 * t434;
t385 = qJD(2) * t412 - t409 * t434;
t410 = sin(qJ(2));
t435 = qJD(1) * t410;
t400 = qJD(4) + t435;
t450 = t458 * t384 - t444 * t385 - t453 * t400;
t449 = -t444 * t384 + t459 * t385 + t455 * t400;
t416 = qJD(1) ^ 2;
t411 = sin(qJ(1));
t414 = cos(qJ(1));
t425 = -g(1) * t414 - g(2) * t411;
t376 = -pkin(1) * t416 + qJDD(1) * pkin(6) + t425;
t363 = -g(3) * t410 + t413 * t376;
t386 = (-t413 * pkin(2) - t410 * qJ(3)) * qJD(1);
t415 = qJD(2) ^ 2;
t335 = pkin(2) * t415 - qJDD(2) * qJ(3) + qJD(2) * t457 - t386 * t434 - t363;
t448 = pkin(6) * t416;
t447 = mrSges(3,1) - mrSges(4,2);
t446 = -mrSges(5,3) - mrSges(6,2);
t433 = qJD(1) * qJD(2);
t429 = t410 * t433;
t390 = qJDD(1) * t413 - t429;
t397 = pkin(3) * t435 - qJD(2) * pkin(7);
t408 = t413 ^ 2;
t430 = t413 * t433;
t389 = qJDD(1) * t410 + t430;
t428 = g(1) * t411 - t414 * g(2);
t423 = -qJDD(1) * pkin(1) - t428;
t420 = pkin(2) * t429 + t435 * t457 + (-t389 - t430) * qJ(3) + t423;
t327 = -t397 * t435 + (-pkin(3) * t408 - pkin(6)) * t416 + (-pkin(2) - pkin(7)) * t390 + t420;
t362 = -t413 * g(3) - t410 * t376;
t336 = -qJDD(2) * pkin(2) - qJ(3) * t415 + t386 * t435 + qJDD(3) - t362;
t331 = (-t410 * t413 * t416 - qJDD(2)) * pkin(7) + (t389 - t430) * pkin(3) + t336;
t325 = t412 * t327 + t409 * t331;
t441 = t453 * t384 - t455 * t385 - t451 * t400;
t356 = mrSges(6,1) * t384 - mrSges(6,3) * t385;
t440 = -mrSges(5,1) * t384 - mrSges(5,2) * t385 - t356;
t439 = t452 * qJD(2) + (t410 * t443 + t413 * t442) * qJD(1);
t438 = t442 * qJD(2) + (t410 * t445 + t413 * t454) * qJD(1);
t437 = t443 * qJD(2) + (t410 * t456 + t413 * t445) * qJD(1);
t395 = -mrSges(4,1) * t434 - qJD(2) * mrSges(4,3);
t436 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t434 - t395;
t355 = pkin(4) * t384 - qJ(5) * t385;
t383 = qJDD(4) + t389;
t398 = t400 ^ 2;
t319 = -pkin(4) * t398 + qJ(5) * t383 + 0.2e1 * qJD(5) * t400 - t355 * t384 + t325;
t361 = -mrSges(6,1) * t400 + mrSges(6,2) * t385;
t431 = m(6) * t319 + t383 * mrSges(6,3) + t400 * t361;
t351 = qJD(4) * t385 + qJDD(2) * t409 + t412 * t390;
t360 = mrSges(5,1) * t400 - mrSges(5,3) * t385;
t311 = m(5) * t325 - mrSges(5,2) * t383 + t446 * t351 - t360 * t400 + t440 * t384 + t431;
t324 = -t327 * t409 + t331 * t412;
t352 = -qJD(4) * t384 + qJDD(2) * t412 - t390 * t409;
t358 = -mrSges(5,2) * t400 - mrSges(5,3) * t384;
t320 = -pkin(4) * t383 - qJ(5) * t398 + t355 * t385 + qJDD(5) - t324;
t359 = -mrSges(6,2) * t384 + mrSges(6,3) * t400;
t426 = -m(6) * t320 + t383 * mrSges(6,1) + t400 * t359;
t312 = m(5) * t324 + mrSges(5,1) * t383 + t446 * t352 + t358 * t400 + t440 * t385 + t426;
t427 = t412 * t311 - t312 * t409;
t309 = t409 * t311 + t412 * t312;
t332 = -pkin(2) * t390 + t420 - t448;
t424 = m(4) * t332 + t427;
t330 = -pkin(7) * t408 * t416 + pkin(3) * t390 + qJD(2) * t397 - t335;
t322 = -0.2e1 * qJD(5) * t385 + (t384 * t400 - t352) * qJ(5) + (t385 * t400 + t351) * pkin(4) + t330;
t316 = m(6) * t322 + t351 * mrSges(6,1) - t352 * mrSges(6,3) + t384 * t359 - t385 * t361;
t421 = m(4) * t336 + t389 * mrSges(4,1) + t309;
t419 = m(5) * t330 + t351 * mrSges(5,1) + t352 * mrSges(5,2) + t384 * t358 + t385 * t360 + t316;
t387 = (t413 * mrSges(4,2) - t410 * mrSges(4,3)) * qJD(1);
t396 = mrSges(4,1) * t435 + qJD(2) * mrSges(4,2);
t418 = -m(4) * t335 + qJDD(2) * mrSges(4,3) + qJD(2) * t396 + t387 * t434 + t419;
t315 = mrSges(6,2) * t352 + t356 * t385 - t426;
t417 = mrSges(5,1) * t324 - mrSges(6,1) * t320 - mrSges(5,2) * t325 + mrSges(6,3) * t319 - pkin(4) * t315 + qJ(5) * t431 - t450 * t385 + (-qJ(5) * t356 + t449) * t384 + t451 * t383 + t455 * t352 + (-qJ(5) * mrSges(6,2) - t453) * t351;
t393 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t435;
t388 = (-t413 * mrSges(3,1) + t410 * mrSges(3,2)) * qJD(1);
t375 = t423 - t448;
t308 = mrSges(5,2) * t330 + mrSges(6,2) * t320 - mrSges(5,3) * t324 - mrSges(6,3) * t322 - qJ(5) * t316 - t444 * t351 + t459 * t352 + t455 * t383 + t441 * t384 + t450 * t400;
t307 = -mrSges(5,1) * t330 - mrSges(6,1) * t322 + mrSges(6,2) * t319 + mrSges(5,3) * t325 - pkin(4) * t316 - t458 * t351 + t444 * t352 + t453 * t383 + t441 * t385 + t449 * t400;
t306 = qJDD(2) * mrSges(4,2) + qJD(2) * t395 + t387 * t435 + t421;
t305 = mrSges(4,2) * t390 - mrSges(4,3) * t389 + (t395 * t413 - t396 * t410) * qJD(1) + t424;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t428 - mrSges(2,2) * t425 + t410 * (mrSges(4,1) * t336 + mrSges(3,2) * t375 - mrSges(3,3) * t362 - mrSges(4,3) * t332 + pkin(3) * t309 - qJ(3) * t305 - t438 * qJD(2) + t443 * qJDD(2) + t456 * t389 + t445 * t390 + t439 * t434 + t417) + t413 * (-mrSges(3,1) * t375 - mrSges(4,1) * t335 + mrSges(4,2) * t332 + mrSges(3,3) * t363 - pkin(2) * t305 + pkin(3) * t419 - pkin(7) * t427 + t437 * qJD(2) + t442 * qJDD(2) - t412 * t307 - t409 * t308 + t445 * t389 + t454 * t390 - t439 * t435) + pkin(1) * (-m(3) * t375 + t447 * t390 + (-mrSges(3,2) + mrSges(4,3)) * t389 + (t436 * t413 + (-t393 + t396) * t410) * qJD(1) - t424) + pkin(6) * (t413 * (t418 + t388 * t434 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t390 - qJD(2) * t393 + m(3) * t363) + (-m(3) * t362 + t389 * mrSges(3,3) - t447 * qJDD(2) - t436 * qJD(2) + (t387 + t388) * t435 + t421) * t410); mrSges(3,1) * t362 - mrSges(3,2) * t363 + mrSges(4,2) * t336 - mrSges(4,3) * t335 + t412 * t308 - t409 * t307 - pkin(7) * t309 - pkin(2) * t306 + qJ(3) * t418 + (qJ(3) * mrSges(4,1) + t442) * t390 + t443 * t389 + t452 * qJDD(2) + (t438 * t410 - t437 * t413) * qJD(1); t306; t417; t315;];
tauJ = t1;
