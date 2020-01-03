% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:45
% EndTime: 2019-12-31 20:57:48
% DurationCPUTime: 1.69s
% Computational Cost: add. (5101->235), mult. (10717->272), div. (0->0), fcn. (6533->6), ass. (0->92)
t448 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t480 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t471 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t479 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t473 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t427 = sin(qJ(3));
t430 = cos(qJ(2));
t451 = qJD(1) * t430;
t428 = sin(qJ(2));
t452 = qJD(1) * t428;
t464 = cos(qJ(3));
t402 = t427 * t452 - t451 * t464;
t450 = qJD(1) * qJD(2);
t409 = qJDD(1) * t428 + t430 * t450;
t410 = qJDD(1) * t430 - t428 * t450;
t369 = -t402 * qJD(3) + t409 * t464 + t427 * t410;
t413 = qJD(2) * pkin(2) - pkin(7) * t452;
t426 = t430 ^ 2;
t432 = qJD(1) ^ 2;
t429 = sin(qJ(1));
t431 = cos(qJ(1));
t453 = t429 * g(1) - t431 * g(2);
t438 = qJDD(1) * pkin(1) + t453;
t370 = -pkin(2) * t410 + t413 * t452 - (pkin(7) * t426 + pkin(6)) * t432 - t438;
t425 = qJD(2) + qJD(3);
t458 = t402 * t425;
t475 = t370 + (-t369 + t458) * qJ(4);
t403 = (t427 * t430 + t428 * t464) * qJD(1);
t474 = t480 * t402 + t448 * t403 + t471 * t425;
t472 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t470 = -t448 * t402 + t479 * t403 + t473 * t425;
t439 = -g(1) * t431 - g(2) * t429;
t405 = -pkin(1) * t432 + qJDD(1) * pkin(6) + t439;
t457 = t405 * t428;
t462 = pkin(2) * t432;
t354 = qJDD(2) * pkin(2) - pkin(7) * t409 - t457 + (pkin(7) * t450 + t428 * t462 - g(3)) * t430;
t389 = -g(3) * t428 + t430 * t405;
t355 = pkin(7) * t410 - qJD(2) * t413 - t426 * t462 + t389;
t346 = t354 * t464 - t427 * t355;
t382 = pkin(3) * t402 - qJ(4) * t403;
t423 = t425 ^ 2;
t424 = qJDD(2) + qJDD(3);
t344 = -t424 * pkin(3) - t423 * qJ(4) + t403 * t382 + qJDD(4) - t346;
t396 = -mrSges(5,2) * t402 + mrSges(5,3) * t425;
t467 = -m(5) * t344 + t424 * mrSges(5,1) + t425 * t396;
t466 = -0.2e1 * t403;
t465 = 2 * qJD(4);
t460 = -mrSges(6,2) - mrSges(5,3);
t459 = -mrSges(4,3) - mrSges(5,2);
t390 = mrSges(6,2) * t425 + mrSges(6,3) * t402;
t391 = -mrSges(4,2) * t425 - mrSges(4,3) * t402;
t335 = qJD(5) * t466 + (-t369 - t458) * qJ(5) + (t402 * t403 - t424) * pkin(4) + t344;
t384 = -mrSges(6,1) * t402 + mrSges(6,2) * t403;
t440 = -m(6) * t335 + t369 * mrSges(6,3) + t403 * t384;
t383 = mrSges(5,1) * t402 - mrSges(5,3) * t403;
t455 = -mrSges(4,1) * t402 - mrSges(4,2) * t403 - t383;
t324 = m(4) * t346 + (t390 + t391) * t425 + (mrSges(4,1) + mrSges(6,1)) * t424 + t455 * t403 + t459 * t369 + t440 + t467;
t347 = t427 * t354 + t464 * t355;
t368 = qJD(3) * t403 + t409 * t427 - t410 * t464;
t343 = -pkin(3) * t423 + t424 * qJ(4) - t402 * t382 + t425 * t465 + t347;
t395 = -mrSges(5,1) * t425 + mrSges(5,2) * t403;
t392 = -pkin(4) * t425 - qJ(5) * t403;
t398 = t402 ^ 2;
t338 = -pkin(4) * t398 + qJ(5) * t368 + 0.2e1 * qJD(5) * t402 + t392 * t425 + t343;
t447 = m(6) * t338 + t368 * mrSges(6,3) + t402 * t384;
t437 = m(5) * t343 + t424 * mrSges(5,3) + t425 * t395 + t447;
t393 = -mrSges(6,1) * t425 - mrSges(6,3) * t403;
t454 = -mrSges(4,1) * t425 + mrSges(4,3) * t403 + t393;
t327 = m(4) * t347 + t454 * t425 + (-mrSges(4,2) + mrSges(6,2)) * t424 + t455 * t402 + t459 * t368 + t437;
t322 = t464 * t324 + t427 * t327;
t446 = t471 * t402 - t473 * t403 + t472 * t425;
t443 = -t324 * t427 + t464 * t327;
t334 = -qJ(5) * t398 + qJDD(5) + (-pkin(3) - pkin(4)) * t368 + (-pkin(3) * t425 + t392 + t465) * t403 - t475;
t441 = m(6) * t334 - t368 * mrSges(6,1) - t402 * t390;
t340 = qJD(4) * t466 + (t403 * t425 + t368) * pkin(3) + t475;
t436 = m(5) * t340 + t368 * mrSges(5,1) + t402 * t396 - t441;
t333 = -mrSges(6,1) * t424 - t390 * t425 - t440;
t434 = m(4) * t370 + t368 * mrSges(4,1) + t402 * t391 + t436;
t330 = mrSges(5,2) * t369 + t383 * t403 + t333 - t467;
t433 = -mrSges(5,1) * t344 - mrSges(6,1) * t335 - mrSges(4,2) * t347 - pkin(4) * t333 - pkin(3) * t330 + qJ(4) * (t393 * t425 + t437) + mrSges(6,2) * t338 + mrSges(5,3) * t343 + mrSges(4,1) * t346 + t474 * t403 + (qJ(4) * mrSges(6,2) - t472) * t424 + (-qJ(4) * t383 + t470) * t402 + t473 * t369 + (-qJ(4) * mrSges(5,2) - t471) * t368;
t412 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t451;
t411 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t452;
t408 = (-mrSges(3,1) * t430 + mrSges(3,2) * t428) * qJD(1);
t404 = -pkin(6) * t432 - t438;
t401 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t428 + Ifges(3,4) * t430) * qJD(1);
t400 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t428 + Ifges(3,2) * t430) * qJD(1);
t388 = -g(3) * t430 - t457;
t332 = t369 * mrSges(6,2) + t403 * t393 + t441;
t328 = (-t393 - t395) * t403 + t460 * t369 + t436;
t321 = mrSges(4,2) * t370 + mrSges(5,2) * t344 + mrSges(6,2) * t334 - mrSges(4,3) * t346 - mrSges(5,3) * t340 - mrSges(6,3) * t335 - qJ(4) * t328 - qJ(5) * t333 - t448 * t368 + t479 * t369 + t446 * t402 + t473 * t424 - t474 * t425;
t320 = -mrSges(4,1) * t370 + mrSges(4,3) * t347 - mrSges(5,1) * t340 + mrSges(5,2) * t343 + mrSges(6,1) * t334 - mrSges(6,3) * t338 + pkin(4) * t332 - qJ(5) * t447 - pkin(3) * t328 + (-qJ(5) * t393 + t470) * t425 + (-qJ(5) * mrSges(6,2) + t471) * t424 + t446 * t403 + t448 * t369 + t480 * t368;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t453 - mrSges(2,2) * t439 + t428 * (mrSges(3,2) * t404 - mrSges(3,3) * t388 + Ifges(3,1) * t409 + Ifges(3,4) * t410 + Ifges(3,5) * qJDD(2) - pkin(7) * t322 - qJD(2) * t400 - t427 * t320 + t464 * t321) + t430 * (-mrSges(3,1) * t404 + mrSges(3,3) * t389 + Ifges(3,4) * t409 + Ifges(3,2) * t410 + Ifges(3,6) * qJDD(2) - pkin(2) * t434 + pkin(7) * t443 + qJD(2) * t401 + t464 * t320 + t427 * t321) + pkin(1) * (-m(3) * t404 + t410 * mrSges(3,1) - t409 * mrSges(3,2) - t411 * t452 + t412 * t451 - t434) + pkin(6) * (t430 * (m(3) * t389 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t410 - qJD(2) * t411 + t408 * t451 + t443) - t428 * (m(3) * t388 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t409 + qJD(2) * t412 - t408 * t452 + t322)) + (-t369 * (mrSges(4,2) + t460) + t403 * (t395 + t454)) * (t430 * pkin(2) + pkin(1)); t433 + pkin(2) * t322 + Ifges(3,3) * qJDD(2) + (t428 * t400 - t430 * t401) * qJD(1) + mrSges(3,1) * t388 - mrSges(3,2) * t389 + Ifges(3,5) * t409 + Ifges(3,6) * t410; t433; t330; t332;];
tauJ = t1;
