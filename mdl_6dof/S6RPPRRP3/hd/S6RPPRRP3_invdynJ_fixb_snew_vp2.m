% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-05-05 14:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:50:53
% EndTime: 2019-05-05 14:50:55
% DurationCPUTime: 1.24s
% Computational Cost: add. (4854->211), mult. (8816->250), div. (0->0), fcn. (4798->8), ass. (0->90)
t442 = Ifges(6,1) + Ifges(7,1);
t430 = Ifges(6,4) - Ifges(7,5);
t438 = -Ifges(6,5) - Ifges(7,4);
t441 = Ifges(6,2) + Ifges(7,3);
t428 = Ifges(6,6) - Ifges(7,6);
t399 = sin(qJ(5));
t402 = cos(qJ(4));
t422 = qJD(1) * t402;
t432 = cos(qJ(5));
t374 = -qJD(4) * t432 + t399 * t422;
t400 = sin(qJ(4));
t421 = qJD(1) * qJD(4);
t417 = t400 * t421;
t381 = qJDD(1) * t402 - t417;
t349 = -t374 * qJD(5) + t399 * qJDD(4) + t381 * t432;
t375 = t399 * qJD(4) + t422 * t432;
t355 = mrSges(7,1) * t374 - mrSges(7,3) * t375;
t405 = qJD(1) ^ 2;
t433 = -pkin(2) - pkin(7);
t401 = sin(qJ(1));
t403 = cos(qJ(1));
t415 = t401 * g(1) - g(2) * t403;
t376 = qJDD(1) * pkin(1) + t415;
t412 = -g(1) * t403 - g(2) * t401;
t378 = -pkin(1) * t405 + t412;
t397 = sin(pkin(9));
t398 = cos(pkin(9));
t351 = t397 * t376 + t398 * t378;
t434 = -qJDD(1) * qJ(3) - 0.2e1 * qJD(3) * qJD(1) - t351;
t334 = t405 * t433 - t434;
t416 = t402 * t421;
t380 = -qJDD(1) * t400 - t416;
t326 = (-t381 + t417) * pkin(8) + (-t380 + t416) * pkin(4) + t334;
t350 = t398 * t376 - t397 * t378;
t410 = -t405 * qJ(3) + qJDD(3) - t350;
t335 = qJDD(1) * t433 + t410;
t394 = -g(3) + qJDD(2);
t331 = t400 * t335 + t402 * t394;
t379 = (pkin(4) * t400 - pkin(8) * t402) * qJD(1);
t404 = qJD(4) ^ 2;
t423 = qJD(1) * t400;
t329 = -pkin(4) * t404 + qJDD(4) * pkin(8) - t379 * t423 + t331;
t323 = t326 * t432 - t399 * t329;
t354 = pkin(5) * t374 - qJ(6) * t375;
t373 = qJDD(5) - t380;
t385 = qJD(5) + t423;
t384 = t385 ^ 2;
t321 = -t373 * pkin(5) - t384 * qJ(6) + t375 * t354 + qJDD(6) - t323;
t360 = -mrSges(7,2) * t374 + mrSges(7,3) * t385;
t413 = -m(7) * t321 + t373 * mrSges(7,1) + t385 * t360;
t317 = t349 * mrSges(7,2) + t375 * t355 - t413;
t324 = t399 * t326 + t432 * t329;
t320 = -pkin(5) * t384 + qJ(6) * t373 + 0.2e1 * qJD(6) * t385 - t354 * t374 + t324;
t348 = t375 * qJD(5) - qJDD(4) * t432 + t399 * t381;
t359 = -mrSges(7,1) * t385 + mrSges(7,2) * t375;
t419 = m(7) * t320 + t373 * mrSges(7,3) + t385 * t359;
t426 = t441 * t374 - t430 * t375 - t428 * t385;
t435 = t430 * t374 - t442 * t375 + t438 * t385;
t437 = -Ifges(6,3) - Ifges(7,2);
t440 = -t438 * t349 - t435 * t374 - t428 * t348 - t437 * t373 + mrSges(6,1) * t323 - mrSges(7,1) * t321 - mrSges(6,2) * t324 + mrSges(7,3) * t320 - pkin(5) * t317 + qJ(6) * (-t348 * mrSges(7,2) - t374 * t355 + t419) - t426 * t375;
t431 = -mrSges(6,3) - mrSges(7,2);
t358 = mrSges(6,1) * t385 - mrSges(6,3) * t375;
t424 = -mrSges(6,1) * t374 - mrSges(6,2) * t375 - t355;
t313 = m(6) * t324 - t373 * mrSges(6,2) + t348 * t431 - t385 * t358 + t374 * t424 + t419;
t357 = -mrSges(6,2) * t385 - mrSges(6,3) * t374;
t315 = m(6) * t323 + t373 * mrSges(6,1) + t349 * t431 + t385 * t357 + t375 * t424 + t413;
t309 = t399 * t313 + t432 * t315;
t427 = t428 * t374 + t438 * t375 + t437 * t385;
t414 = t432 * t313 - t315 * t399;
t330 = t402 * t335 - t400 * t394;
t377 = (mrSges(5,1) * t400 + mrSges(5,2) * t402) * qJD(1);
t383 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t422;
t306 = m(5) * t331 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t380 - qJD(4) * t383 - t377 * t423 + t414;
t382 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t423;
t328 = -qJDD(4) * pkin(4) - t404 * pkin(8) + t379 * t422 - t330;
t322 = -0.2e1 * qJD(6) * t375 + (t374 * t385 - t349) * qJ(6) + (t375 * t385 + t348) * pkin(5) + t328;
t318 = m(7) * t322 + mrSges(7,1) * t348 - t349 * mrSges(7,3) - t375 * t359 + t360 * t374;
t406 = -m(6) * t328 - t348 * mrSges(6,1) - mrSges(6,2) * t349 - t374 * t357 - t358 * t375 - t318;
t310 = m(5) * t330 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t381 + qJD(4) * t382 - t377 * t422 + t406;
t411 = t400 * t306 + t402 * t310;
t337 = -qJDD(1) * pkin(2) + t410;
t409 = m(4) * t337 - t405 * mrSges(4,3) + t411;
t336 = t405 * pkin(2) + t434;
t408 = -m(4) * t336 + m(5) * t334 - mrSges(5,1) * t380 + t405 * mrSges(4,2) + t381 * mrSges(5,2) + qJDD(1) * mrSges(4,3) + t382 * t423 + t383 * t422 + t309;
t367 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t402 - Ifges(5,4) * t400) * qJD(1);
t366 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t402 - Ifges(5,2) * t400) * qJD(1);
t308 = mrSges(6,2) * t328 + mrSges(7,2) * t321 - mrSges(6,3) * t323 - mrSges(7,3) * t322 - qJ(6) * t318 - t430 * t348 + t442 * t349 - t438 * t373 + t427 * t374 + t426 * t385;
t307 = -mrSges(6,1) * t328 - mrSges(7,1) * t322 + mrSges(7,2) * t320 + mrSges(6,3) * t324 - pkin(5) * t318 - t441 * t348 + t430 * t349 + t428 * t373 + t427 * t375 - t435 * t385;
t305 = qJDD(1) * mrSges(4,2) + t409;
t1 = [pkin(1) * (t397 * (m(3) * t351 - mrSges(3,1) * t405 + t408) + t398 * (m(3) * t350 - t405 * mrSges(3,2) - t409)) + mrSges(2,1) * t415 - mrSges(2,2) * t412 - pkin(2) * t305 + qJ(3) * t408 + t402 * (mrSges(5,2) * t334 - mrSges(5,3) * t330 + Ifges(5,1) * t381 + Ifges(5,4) * t380 + Ifges(5,5) * qJDD(4) - pkin(8) * t309 - qJD(4) * t366 - t399 * t307 + t308 * t432) - t400 * (-mrSges(5,1) * t334 + mrSges(5,3) * t331 + Ifges(5,4) * t381 + Ifges(5,2) * t380 + Ifges(5,6) * qJDD(4) - pkin(4) * t309 + qJD(4) * t367 - t440) - pkin(7) * t411 + mrSges(3,1) * t350 - mrSges(3,2) * t351 + mrSges(4,2) * t337 - mrSges(4,3) * t336 + (pkin(1) * (-t397 * mrSges(3,2) + t398 * (mrSges(3,1) - mrSges(4,2))) + Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1); t402 * t306 - t400 * t310 + (m(3) + m(4)) * t394; t305; Ifges(5,5) * t381 + Ifges(5,6) * t380 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t330 - mrSges(5,2) * t331 + t399 * t308 + t432 * t307 + pkin(4) * t406 + pkin(8) * t414 + (t366 * t402 + t367 * t400) * qJD(1); t440; t317;];
tauJ  = t1;
