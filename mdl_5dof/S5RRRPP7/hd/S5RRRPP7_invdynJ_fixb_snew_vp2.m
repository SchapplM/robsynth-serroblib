% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP7
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:04:06
% EndTime: 2019-12-31 21:04:09
% DurationCPUTime: 1.17s
% Computational Cost: add. (4959->232), mult. (9665->273), div. (0->0), fcn. (5681->6), ass. (0->88)
t443 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t428 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t427 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t442 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t426 = -Ifges(5,6) + Ifges(6,6) + Ifges(4,6);
t441 = Ifges(6,3) + Ifges(4,3) + Ifges(5,2);
t410 = qJD(1) ^ 2;
t406 = sin(qJ(1));
t408 = cos(qJ(1));
t419 = g(1) * t406 - t408 * g(2);
t380 = -qJDD(1) * pkin(1) - pkin(6) * t410 - t419;
t405 = sin(qJ(2));
t407 = cos(qJ(2));
t429 = qJD(1) * qJD(2);
t420 = t407 * t429;
t391 = qJDD(1) * t405 + t420;
t399 = t405 * t429;
t392 = qJDD(1) * t407 - t399;
t329 = (-t391 - t420) * pkin(7) + (-t392 + t399) * pkin(2) + t380;
t415 = -g(1) * t408 - g(2) * t406;
t381 = -pkin(1) * t410 + qJDD(1) * pkin(6) + t415;
t373 = -g(3) * t405 + t407 * t381;
t390 = (-t407 * pkin(2) - t405 * pkin(7)) * qJD(1);
t409 = qJD(2) ^ 2;
t430 = qJD(1) * t407;
t333 = -pkin(2) * t409 + qJDD(2) * pkin(7) + t390 * t430 + t373;
t404 = sin(qJ(3));
t436 = cos(qJ(3));
t326 = t436 * t329 - t404 * t333;
t431 = qJD(1) * t405;
t387 = -t436 * qJD(2) + t404 * t431;
t388 = t404 * qJD(2) + t436 * t431;
t361 = pkin(3) * t387 - qJ(4) * t388;
t386 = -qJDD(3) + t392;
t397 = -qJD(3) + t430;
t396 = t397 ^ 2;
t325 = t386 * pkin(3) - t396 * qJ(4) + t388 * t361 + qJDD(4) - t326;
t371 = -mrSges(5,2) * t387 - mrSges(5,3) * t397;
t440 = -m(5) * t325 - t386 * mrSges(5,1) - t397 * t371;
t358 = -t387 * qJD(3) + t404 * qJDD(2) + t436 * t391;
t372 = -t407 * g(3) - t405 * t381;
t413 = qJDD(2) * pkin(2) + pkin(7) * t409 - t390 * t431 + t372;
t433 = t387 * t397;
t439 = -(t358 + t433) * qJ(4) - t413;
t365 = -mrSges(6,2) * t397 + mrSges(6,3) * t387;
t317 = -0.2e1 * qJD(5) * t388 + (-t358 + t433) * qJ(5) + (t387 * t388 + t386) * pkin(4) + t325;
t363 = -mrSges(6,1) * t387 + mrSges(6,2) * t388;
t416 = -m(6) * t317 + t358 * mrSges(6,3) + t388 * t363;
t315 = mrSges(6,1) * t386 + t365 * t397 - t416;
t362 = mrSges(5,1) * t387 - mrSges(5,3) * t388;
t313 = mrSges(5,2) * t358 + t362 * t388 + t315 - t440;
t327 = t404 * t329 + t436 * t333;
t437 = -2 * qJD(4);
t323 = -pkin(3) * t396 - t386 * qJ(4) - t387 * t361 + t397 * t437 + t327;
t357 = qJD(3) * t388 - t436 * qJDD(2) + t391 * t404;
t367 = pkin(4) * t397 - qJ(5) * t388;
t385 = t387 ^ 2;
t319 = -pkin(4) * t385 + qJ(5) * t357 + 0.2e1 * qJD(5) * t387 - t367 * t397 + t323;
t368 = mrSges(6,1) * t397 - mrSges(6,3) * t388;
t370 = mrSges(5,1) * t397 + mrSges(5,2) * t388;
t424 = m(6) * t319 + t357 * mrSges(6,3) + t387 * t363;
t414 = m(5) * t323 - t386 * mrSges(5,3) - t397 * t370 + t424;
t421 = -t428 * t387 + t443 * t388 - t427 * t397;
t422 = t442 * t387 + t428 * t388 - t426 * t397;
t438 = -t357 * t426 + t358 * t427 - t441 * t386 + t421 * t387 + t388 * t422 + mrSges(4,1) * t326 - mrSges(5,1) * t325 - mrSges(6,1) * t317 - mrSges(4,2) * t327 + mrSges(6,2) * t319 + mrSges(5,3) * t323 - pkin(3) * t313 - pkin(4) * t315 + qJ(4) * (-mrSges(5,2) * t357 - mrSges(6,2) * t386 - t362 * t387 - t368 * t397 + t414);
t434 = -mrSges(4,3) - mrSges(5,2);
t432 = -mrSges(4,1) * t387 - mrSges(4,2) * t388 - t362;
t423 = t426 * t387 - t427 * t388 + t441 * t397;
t366 = mrSges(4,2) * t397 - mrSges(4,3) * t387;
t310 = m(4) * t326 + (-t365 - t366) * t397 + t432 * t388 + (-mrSges(4,1) - mrSges(6,1)) * t386 + t434 * t358 + t416 + t440;
t369 = -mrSges(4,1) * t397 - mrSges(4,3) * t388;
t311 = m(4) * t327 + (-t368 + t369) * t397 + t432 * t387 + (mrSges(4,2) - mrSges(6,2)) * t386 + t434 * t357 + t414;
t418 = -t310 * t404 + t436 * t311;
t321 = -qJ(5) * t385 + qJDD(5) + (-pkin(3) - pkin(4)) * t357 + (pkin(3) * t397 + (2 * qJD(4)) + t367) * t388 - t439;
t316 = m(6) * t321 - t357 * mrSges(6,1) + t358 * mrSges(6,2) - t387 * t365 + t388 * t368;
t308 = t436 * t310 + t404 * t311;
t324 = t388 * t437 + (-t388 * t397 + t357) * pkin(3) + t439;
t312 = m(5) * t324 + t357 * mrSges(5,1) - t358 * mrSges(5,3) - t388 * t370 + t387 * t371 - t316;
t411 = m(4) * t413 - t357 * mrSges(4,1) - t358 * mrSges(4,2) - t387 * t366 - t388 * t369 - t312;
t394 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t430;
t393 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t431;
t389 = (-t407 * mrSges(3,1) + t405 * mrSges(3,2)) * qJD(1);
t379 = Ifges(3,5) * qJD(2) + (t405 * Ifges(3,1) + t407 * Ifges(3,4)) * qJD(1);
t378 = Ifges(3,6) * qJD(2) + (t405 * Ifges(3,4) + t407 * Ifges(3,2)) * qJD(1);
t377 = Ifges(3,3) * qJD(2) + (t405 * Ifges(3,5) + t407 * Ifges(3,6)) * qJD(1);
t307 = -mrSges(4,2) * t413 + mrSges(5,2) * t325 + mrSges(6,2) * t321 - mrSges(4,3) * t326 - mrSges(5,3) * t324 - mrSges(6,3) * t317 - qJ(4) * t312 - qJ(5) * t315 - t428 * t357 + t443 * t358 - t427 * t386 + t423 * t387 + t422 * t397;
t306 = mrSges(4,1) * t413 + mrSges(4,3) * t327 - mrSges(5,1) * t324 + mrSges(5,2) * t323 + mrSges(6,1) * t321 - mrSges(6,3) * t319 + pkin(4) * t316 - qJ(5) * t424 - pkin(3) * t312 + (qJ(5) * t368 - t421) * t397 + t423 * t388 + (qJ(5) * mrSges(6,2) - t426) * t386 + t428 * t358 + t442 * t357;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t419 - mrSges(2,2) * t415 + t405 * (mrSges(3,2) * t380 - mrSges(3,3) * t372 + Ifges(3,1) * t391 + Ifges(3,4) * t392 + Ifges(3,5) * qJDD(2) - pkin(7) * t308 - qJD(2) * t378 - t404 * t306 + t436 * t307 + t377 * t430) + t407 * (-mrSges(3,1) * t380 + mrSges(3,3) * t373 + Ifges(3,4) * t391 + Ifges(3,2) * t392 + Ifges(3,6) * qJDD(2) - pkin(2) * t308 + qJD(2) * t379 - t377 * t431 - t438) + pkin(1) * (-m(3) * t380 + t392 * mrSges(3,1) - t391 * mrSges(3,2) + (-t393 * t405 + t394 * t407) * qJD(1) - t308) + pkin(6) * (t407 * (m(3) * t373 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t392 - qJD(2) * t393 + t389 * t430 + t418) - t405 * (m(3) * t372 + qJDD(2) * mrSges(3,1) - t391 * mrSges(3,3) + qJD(2) * t394 - t389 * t431 + t411)); Ifges(3,5) * t391 + Ifges(3,6) * t392 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t372 - mrSges(3,2) * t373 + t404 * t307 + t436 * t306 + pkin(2) * t411 + pkin(7) * t418 + (t405 * t378 - t407 * t379) * qJD(1); t438; t313; t316;];
tauJ = t1;
