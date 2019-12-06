% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:45
% EndTime: 2019-12-05 16:45:47
% DurationCPUTime: 1.35s
% Computational Cost: add. (13925->195), mult. (17324->240), div. (0->0), fcn. (9126->8), ass. (0->82)
t441 = Ifges(5,1) + Ifges(6,1);
t436 = Ifges(5,4) - Ifges(6,5);
t435 = Ifges(5,5) + Ifges(6,4);
t440 = Ifges(5,2) + Ifges(6,3);
t434 = Ifges(5,6) - Ifges(6,6);
t439 = Ifges(5,3) + Ifges(6,2);
t438 = m(3) + m(4);
t437 = mrSges(5,3) + mrSges(6,2);
t433 = cos(pkin(8));
t409 = sin(pkin(8));
t397 = g(1) * t409 - t433 * g(2);
t413 = cos(qJ(4));
t432 = t397 * t413;
t405 = qJD(2) + qJD(3);
t410 = sin(qJ(4));
t431 = t405 * t410;
t430 = t405 * t413;
t398 = -t433 * g(1) - t409 * g(2);
t408 = -g(3) + qJDD(1);
t412 = sin(qJ(2));
t415 = cos(qJ(2));
t371 = -t398 * t412 + t415 * t408;
t368 = qJDD(2) * pkin(2) + t371;
t372 = t415 * t398 + t412 * t408;
t417 = qJD(2) ^ 2;
t369 = -pkin(2) * t417 + t372;
t411 = sin(qJ(3));
t414 = cos(qJ(3));
t364 = t411 * t368 + t414 * t369;
t403 = t405 ^ 2;
t404 = qJDD(2) + qJDD(3);
t362 = -pkin(3) * t403 + pkin(7) * t404 + t364;
t359 = t413 * t362 - t410 * t397;
t386 = (-mrSges(5,1) * t413 + mrSges(5,2) * t410) * t405;
t425 = qJD(4) * t405;
t388 = t404 * t413 - t410 * t425;
t393 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t431;
t384 = (-pkin(4) * t413 - qJ(5) * t410) * t405;
t416 = qJD(4) ^ 2;
t356 = -pkin(4) * t416 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t384 * t430 + t359;
t385 = (-mrSges(6,1) * t413 - mrSges(6,3) * t410) * t405;
t394 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t431;
t420 = m(6) * t356 + qJDD(4) * mrSges(6,3) + qJD(4) * t394 + t385 * t430;
t351 = m(5) * t359 - qJDD(4) * mrSges(5,2) - qJD(4) * t393 + t386 * t430 + t437 * t388 + t420;
t358 = -t362 * t410 - t432;
t387 = t404 * t410 + t413 * t425;
t395 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t430;
t357 = -qJDD(4) * pkin(4) - qJ(5) * t416 + t432 + qJDD(5) + (t384 * t405 + t362) * t410;
t396 = mrSges(6,2) * t430 + qJD(4) * mrSges(6,3);
t419 = -m(6) * t357 + qJDD(4) * mrSges(6,1) + qJD(4) * t396;
t352 = m(5) * t358 + qJDD(4) * mrSges(5,1) + qJD(4) * t395 + (-t385 - t386) * t431 - t437 * t387 + t419;
t421 = t413 * t351 - t410 * t352;
t342 = m(4) * t364 - mrSges(4,1) * t403 - mrSges(4,2) * t404 + t421;
t363 = t368 * t414 - t411 * t369;
t361 = -pkin(3) * t404 - pkin(7) * t403 - t363;
t354 = -pkin(4) * t388 - qJ(5) * t387 + (-0.2e1 * qJD(5) * t410 + (pkin(4) * t410 - qJ(5) * t413) * qJD(4)) * t405 + t361;
t353 = m(6) * t354 - mrSges(6,1) * t388 - t387 * mrSges(6,3) - t394 * t431 - t396 * t430;
t418 = -m(5) * t361 + t388 * mrSges(5,1) - mrSges(5,2) * t387 - t393 * t431 + t395 * t430 - t353;
t347 = m(4) * t363 + mrSges(4,1) * t404 - mrSges(4,2) * t403 + t418;
t337 = t411 * t342 + t414 * t347;
t335 = m(3) * t371 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t417 + t337;
t422 = t414 * t342 - t347 * t411;
t336 = m(3) * t372 - mrSges(3,1) * t417 - qJDD(2) * mrSges(3,2) + t422;
t423 = -t335 * t412 + t415 * t336;
t328 = m(2) * t398 + t423;
t345 = t410 * t351 + t413 * t352;
t344 = (m(2) + t438) * t397 - t345;
t429 = t409 * t328 + t433 * t344;
t329 = t415 * t335 + t412 * t336;
t428 = (-t436 * t410 - t440 * t413) * t405 - t434 * qJD(4);
t427 = (t435 * t410 + t434 * t413) * t405 + t439 * qJD(4);
t426 = (t441 * t410 + t436 * t413) * t405 + t435 * qJD(4);
t424 = t433 * t328 - t344 * t409;
t339 = mrSges(5,2) * t361 + mrSges(6,2) * t357 - mrSges(5,3) * t358 - mrSges(6,3) * t354 - qJ(5) * t353 + t428 * qJD(4) + t435 * qJDD(4) + t441 * t387 + t436 * t388 + t427 * t430;
t338 = -mrSges(5,1) * t361 - mrSges(6,1) * t354 + mrSges(6,2) * t356 + mrSges(5,3) * t359 - pkin(4) * t353 + t426 * qJD(4) + t434 * qJDD(4) + t436 * t387 + t440 * t388 - t427 * t431;
t331 = Ifges(4,6) * t404 + t403 * Ifges(4,5) + mrSges(4,1) * t397 + mrSges(4,3) * t364 - mrSges(5,1) * t358 + mrSges(5,2) * t359 + mrSges(6,1) * t357 - mrSges(6,3) * t356 - pkin(4) * t419 - qJ(5) * t420 - pkin(3) * t345 + (-mrSges(6,2) * qJ(5) - t434) * t388 + (mrSges(6,2) * pkin(4) - t435) * t387 - t439 * qJDD(4) + (t426 * t413 + (pkin(4) * t385 + t428) * t410) * t405;
t330 = -mrSges(4,2) * t397 - mrSges(4,3) * t363 + Ifges(4,5) * t404 - Ifges(4,6) * t403 - pkin(7) * t345 - t338 * t410 + t339 * t413;
t325 = -mrSges(3,2) * t397 - mrSges(3,3) * t371 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t417 - pkin(6) * t337 + t330 * t414 - t331 * t411;
t324 = Ifges(3,6) * qJDD(2) + t417 * Ifges(3,5) + mrSges(3,1) * t397 + mrSges(3,3) * t372 + t411 * t330 + t414 * t331 - pkin(2) * (-m(4) * t397 + t345) + pkin(6) * t422;
t323 = -mrSges(2,1) * t408 - mrSges(3,1) * t371 - mrSges(4,1) * t363 + mrSges(3,2) * t372 + mrSges(4,2) * t364 + mrSges(2,3) * t398 - Ifges(3,3) * qJDD(2) - Ifges(4,3) * t404 - pkin(1) * t329 - pkin(2) * t337 - pkin(3) * t418 - pkin(7) * t421 - t413 * t338 - t410 * t339;
t322 = mrSges(2,2) * t408 - mrSges(2,3) * t397 - pkin(5) * t329 - t324 * t412 + t325 * t415;
t1 = [-m(1) * g(1) + t424; -m(1) * g(2) + t429; -m(1) * g(3) + m(2) * t408 + t329; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t429 + t433 * t322 - t409 * t323; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t424 + t409 * t322 + t433 * t323; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t398 + t412 * t325 + t415 * t324 - pkin(1) * t345 + pkin(5) * t423 + (pkin(1) * t438 + mrSges(2,1)) * t397;];
tauB = t1;
