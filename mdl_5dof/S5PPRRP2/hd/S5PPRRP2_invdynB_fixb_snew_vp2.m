% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:41
% EndTime: 2019-12-05 15:08:42
% DurationCPUTime: 1.13s
% Computational Cost: add. (9242->182), mult. (15900->228), div. (0->0), fcn. (9126->8), ass. (0->78)
t404 = Ifges(5,1) + Ifges(6,1);
t400 = Ifges(5,4) - Ifges(6,5);
t399 = Ifges(6,4) + Ifges(5,5);
t403 = Ifges(5,2) + Ifges(6,3);
t398 = Ifges(5,6) - Ifges(6,6);
t402 = Ifges(5,3) + Ifges(6,2);
t401 = mrSges(5,3) + mrSges(6,2);
t373 = sin(pkin(7));
t375 = cos(pkin(7));
t362 = t373 * g(1) - t375 * g(2);
t361 = qJDD(2) - t362;
t378 = cos(qJ(4));
t397 = t378 * t361;
t363 = -t375 * g(1) - t373 * g(2);
t371 = -g(3) + qJDD(1);
t372 = sin(pkin(8));
t374 = cos(pkin(8));
t339 = -t372 * t363 + t374 * t371;
t340 = t374 * t363 + t372 * t371;
t377 = sin(qJ(3));
t379 = cos(qJ(3));
t335 = t377 * t339 + t379 * t340;
t381 = qJD(3) ^ 2;
t333 = -t381 * pkin(3) + qJDD(3) * pkin(6) + t335;
t376 = sin(qJ(4));
t330 = t378 * t333 + t376 * t361;
t356 = (-mrSges(5,1) * t378 + mrSges(5,2) * t376) * qJD(3);
t390 = qJD(3) * qJD(4);
t358 = t378 * qJDD(3) - t376 * t390;
t392 = qJD(3) * t376;
t364 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t392;
t354 = (-pkin(4) * t378 - qJ(5) * t376) * qJD(3);
t380 = qJD(4) ^ 2;
t391 = qJD(3) * t378;
t327 = -t380 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t354 * t391 + t330;
t355 = (-mrSges(6,1) * t378 - mrSges(6,3) * t376) * qJD(3);
t365 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t392;
t385 = m(6) * t327 + qJDD(4) * mrSges(6,3) + qJD(4) * t365 + t355 * t391;
t322 = m(5) * t330 - qJDD(4) * mrSges(5,2) - qJD(4) * t364 + t356 * t391 + t401 * t358 + t385;
t329 = -t376 * t333 + t397;
t357 = t376 * qJDD(3) + t378 * t390;
t366 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t391;
t328 = -qJDD(4) * pkin(4) - t380 * qJ(5) - t397 + qJDD(5) + (qJD(3) * t354 + t333) * t376;
t367 = mrSges(6,2) * t391 + qJD(4) * mrSges(6,3);
t384 = -m(6) * t328 + qJDD(4) * mrSges(6,1) + qJD(4) * t367;
t323 = m(5) * t329 + qJDD(4) * mrSges(5,1) + qJD(4) * t366 - t401 * t357 + (-t355 - t356) * t392 + t384;
t386 = t378 * t322 - t376 * t323;
t313 = m(4) * t335 - t381 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t386;
t334 = t379 * t339 - t377 * t340;
t332 = -qJDD(3) * pkin(3) - t381 * pkin(6) - t334;
t325 = -t358 * pkin(4) - t357 * qJ(5) + (-0.2e1 * qJD(5) * t376 + (pkin(4) * t376 - qJ(5) * t378) * qJD(4)) * qJD(3) + t332;
t324 = m(6) * t325 - t358 * mrSges(6,1) - t357 * mrSges(6,3) - t365 * t392 - t367 * t391;
t382 = -m(5) * t332 + t358 * mrSges(5,1) - t357 * mrSges(5,2) - t364 * t392 + t366 * t391 - t324;
t318 = m(4) * t334 + qJDD(3) * mrSges(4,1) - t381 * mrSges(4,2) + t382;
t308 = t377 * t313 + t379 * t318;
t306 = m(3) * t339 + t308;
t387 = t379 * t313 - t377 * t318;
t307 = m(3) * t340 + t387;
t388 = -t372 * t306 + t374 * t307;
t299 = m(2) * t363 + t388;
t316 = t376 * t322 + t378 * t323;
t383 = (-m(3) - m(4)) * t361 - t316;
t315 = m(2) * t362 + t383;
t396 = t373 * t299 + t375 * t315;
t300 = t374 * t306 + t372 * t307;
t395 = -t398 * qJD(4) + (-t400 * t376 - t403 * t378) * qJD(3);
t394 = t402 * qJD(4) + (t399 * t376 + t398 * t378) * qJD(3);
t393 = t399 * qJD(4) + (t404 * t376 + t400 * t378) * qJD(3);
t389 = t375 * t299 - t373 * t315;
t310 = mrSges(5,2) * t332 + mrSges(6,2) * t328 - mrSges(5,3) * t329 - mrSges(6,3) * t325 - qJ(5) * t324 + t395 * qJD(4) + t399 * qJDD(4) + t404 * t357 + t400 * t358 + t394 * t391;
t309 = -mrSges(5,1) * t332 - mrSges(6,1) * t325 + mrSges(6,2) * t327 + mrSges(5,3) * t330 - pkin(4) * t324 + t393 * qJD(4) + t398 * qJDD(4) + t400 * t357 + t403 * t358 - t394 * t392;
t302 = Ifges(4,6) * qJDD(3) + t381 * Ifges(4,5) - mrSges(4,1) * t361 + mrSges(4,3) * t335 - mrSges(5,1) * t329 + mrSges(5,2) * t330 + mrSges(6,1) * t328 - mrSges(6,3) * t327 - pkin(4) * t384 - qJ(5) * t385 - pkin(3) * t316 + (-qJ(5) * mrSges(6,2) - t398) * t358 + (pkin(4) * mrSges(6,2) - t399) * t357 - t402 * qJDD(4) + (t393 * t378 + (pkin(4) * t355 + t395) * t376) * qJD(3);
t301 = mrSges(4,2) * t361 - mrSges(4,3) * t334 + Ifges(4,5) * qJDD(3) - t381 * Ifges(4,6) - pkin(6) * t316 - t376 * t309 + t378 * t310;
t296 = mrSges(3,2) * t361 - mrSges(3,3) * t339 - pkin(5) * t308 + t379 * t301 - t377 * t302;
t295 = -mrSges(3,1) * t361 + mrSges(3,3) * t340 + t377 * t301 + t379 * t302 - pkin(2) * (m(4) * t361 + t316) + pkin(5) * t387;
t294 = -mrSges(2,1) * t371 - mrSges(3,1) * t339 - mrSges(4,1) * t334 + mrSges(3,2) * t340 + mrSges(4,2) * t335 + mrSges(2,3) * t363 - Ifges(4,3) * qJDD(3) - pkin(1) * t300 - pkin(2) * t308 - pkin(3) * t382 - pkin(6) * t386 - t378 * t309 - t376 * t310;
t293 = mrSges(2,2) * t371 - mrSges(2,3) * t362 - qJ(2) * t300 - t372 * t295 + t374 * t296;
t1 = [-m(1) * g(1) + t389; -m(1) * g(2) + t396; -m(1) * g(3) + m(2) * t371 + t300; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t396 + t375 * t293 - t373 * t294; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t389 + t373 * t293 + t375 * t294; -mrSges(1,1) * g(2) + mrSges(2,1) * t362 + mrSges(1,2) * g(1) - mrSges(2,2) * t363 + pkin(1) * t383 + qJ(2) * t388 + t374 * t295 + t372 * t296;];
tauB = t1;
