% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:10
% EndTime: 2019-12-31 21:11:11
% DurationCPUTime: 0.90s
% Computational Cost: add. (5952->203), mult. (7538->253), div. (0->0), fcn. (3882->8), ass. (0->87)
t408 = Ifges(4,1) + Ifges(5,1);
t402 = Ifges(4,4) - Ifges(5,5);
t401 = Ifges(4,5) + Ifges(5,4);
t407 = Ifges(4,2) + Ifges(5,3);
t400 = Ifges(4,6) - Ifges(5,6);
t406 = Ifges(4,3) + Ifges(5,2);
t405 = 2 * qJD(4);
t367 = qJD(1) + qJD(2);
t363 = t367 ^ 2;
t404 = t363 * pkin(7);
t403 = mrSges(4,3) + mrSges(5,2);
t373 = sin(qJ(3));
t399 = t367 * t373;
t377 = cos(qJ(3));
t398 = t367 * t377;
t375 = sin(qJ(1));
t379 = cos(qJ(1));
t392 = t375 * g(1) - t379 * g(2);
t351 = qJDD(1) * pkin(1) + t392;
t388 = -t379 * g(1) - t375 * g(2);
t352 = -qJD(1) ^ 2 * pkin(1) + t388;
t374 = sin(qJ(2));
t378 = cos(qJ(2));
t321 = t374 * t351 + t378 * t352;
t365 = qJDD(1) + qJDD(2);
t318 = -t363 * pkin(2) + t365 * pkin(7) + t321;
t307 = -t377 * g(3) - t373 * t318;
t397 = (t373 * t401 + t377 * t400) * t367 + t406 * qJD(3);
t396 = (-t373 * t402 - t407 * t377) * t367 - t400 * qJD(3);
t395 = (t408 * t373 + t377 * t402) * t367 + t401 * qJD(3);
t320 = t378 * t351 - t374 * t352;
t394 = qJD(3) * t377;
t393 = t367 * t394;
t308 = -t373 * g(3) + t377 * t318;
t341 = (-mrSges(5,1) * t377 - mrSges(5,3) * t373) * t367;
t342 = (-mrSges(4,1) * t377 + mrSges(4,2) * t373) * t367;
t343 = t373 * t365 + t393;
t344 = -qJD(3) * t399 + t377 * t365;
t353 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t399;
t355 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t398;
t340 = (-pkin(3) * t377 - qJ(4) * t373) * t367;
t380 = qJD(3) ^ 2;
t301 = -t380 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t405 + t340 * t398 + t308;
t357 = -qJD(3) * pkin(4) - pkin(8) * t399;
t371 = t377 ^ 2;
t297 = -t371 * t363 * pkin(4) - t344 * pkin(8) + qJD(3) * t357 + t301;
t302 = -qJDD(3) * pkin(3) - t380 * qJ(4) + t340 * t399 + qJDD(4) - t307;
t298 = (-t343 + t393) * pkin(8) + (-t363 * t373 * t377 - qJDD(3)) * pkin(4) + t302;
t372 = sin(qJ(5));
t376 = cos(qJ(5));
t293 = -t372 * t297 + t376 * t298;
t333 = (-t372 * t373 - t376 * t377) * t367;
t306 = t333 * qJD(5) + t376 * t343 - t372 * t344;
t334 = (-t372 * t377 + t373 * t376) * t367;
t316 = -t333 * mrSges(6,1) + t334 * mrSges(6,2);
t366 = -qJD(3) + qJD(5);
t322 = -t366 * mrSges(6,2) + t333 * mrSges(6,3);
t364 = -qJDD(3) + qJDD(5);
t291 = m(6) * t293 + t364 * mrSges(6,1) - t306 * mrSges(6,3) - t334 * t316 + t366 * t322;
t294 = t376 * t297 + t372 * t298;
t305 = -t334 * qJD(5) - t372 * t343 - t376 * t344;
t323 = t366 * mrSges(6,1) - t334 * mrSges(6,3);
t292 = m(6) * t294 - t364 * mrSges(6,2) + t305 * mrSges(6,3) + t333 * t316 - t366 * t323;
t285 = t376 * t291 + t372 * t292;
t356 = mrSges(5,2) * t398 + qJD(3) * mrSges(5,3);
t383 = -m(5) * t302 + qJDD(3) * mrSges(5,1) + qJD(3) * t356 - t285;
t354 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t399;
t390 = -t372 * t291 + t376 * t292;
t385 = m(5) * t301 + qJDD(3) * mrSges(5,3) + qJD(3) * t354 + t341 * t398 + t390;
t391 = -t373 * (m(4) * t307 + qJDD(3) * mrSges(4,1) + qJD(3) * t355 + (-t341 - t342) * t399 - t403 * t343 + t383) + t377 * (m(4) * t308 - qJDD(3) * mrSges(4,2) - qJD(3) * t353 + t342 * t398 + t403 * t344 + t385);
t389 = t365 * pkin(2) + t320;
t386 = -t343 * qJ(4) - t389;
t296 = (-pkin(8) * t371 + pkin(7)) * t363 + (pkin(3) + pkin(4)) * t344 + (qJ(4) * t394 + (-pkin(3) * qJD(3) + t357 + t405) * t373) * t367 - t386;
t387 = -m(6) * t296 + t305 * mrSges(6,1) - t306 * mrSges(6,2) + t333 * t322 - t334 * t323;
t309 = Ifges(6,5) * t334 + Ifges(6,6) * t333 + Ifges(6,3) * t366;
t311 = Ifges(6,1) * t334 + Ifges(6,4) * t333 + Ifges(6,5) * t366;
t287 = -mrSges(6,1) * t296 + mrSges(6,3) * t294 + Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t364 - t334 * t309 + t366 * t311;
t310 = Ifges(6,4) * t334 + Ifges(6,2) * t333 + Ifges(6,6) * t366;
t288 = mrSges(6,2) * t296 - mrSges(6,3) * t293 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t364 + t333 * t309 - t366 * t310;
t299 = -t344 * pkin(3) - t404 + (-0.2e1 * qJD(4) * t373 + (pkin(3) * t373 - qJ(4) * t377) * qJD(3)) * t367 + t386;
t289 = m(5) * t299 - t344 * mrSges(5,1) - t343 * mrSges(5,3) - t354 * t399 - t356 * t398 + t387;
t317 = -t389 - t404;
t381 = -m(4) * t317 + t344 * mrSges(4,1) - t343 * mrSges(4,2) - t353 * t399 + t355 * t398 - t289;
t384 = -mrSges(3,2) * t321 + t377 * (-mrSges(4,1) * t317 - mrSges(5,1) * t299 + mrSges(5,2) * t301 + mrSges(4,3) * t308 - pkin(3) * t289 - pkin(4) * t387 - pkin(8) * t390 + t395 * qJD(3) + t400 * qJDD(3) - t376 * t287 - t372 * t288 + t402 * t343 + t407 * t344 - t397 * t399) + t373 * (mrSges(4,2) * t317 + mrSges(5,2) * t302 - mrSges(4,3) * t307 - mrSges(5,3) * t299 - pkin(8) * t285 - qJ(4) * t289 + t396 * qJD(3) + t401 * qJDD(3) - t372 * t287 + t376 * t288 + t408 * t343 + t402 * t344 + t397 * t398) + pkin(7) * t391 + pkin(2) * t381 + mrSges(3,1) * t320 + Ifges(3,3) * t365;
t382 = mrSges(6,1) * t293 - mrSges(6,2) * t294 + Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t364 + t334 * t310 - t333 * t311;
t284 = t343 * mrSges(5,2) + t341 * t399 - t383;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t392 - mrSges(2,2) * t388 + pkin(1) * (t374 * (m(3) * t321 - t363 * mrSges(3,1) - t365 * mrSges(3,2) + t391) + t378 * (m(3) * t320 + t365 * mrSges(3,1) - t363 * mrSges(3,2) + t381)) + t384; t384; mrSges(4,1) * t307 - mrSges(4,2) * t308 + mrSges(5,3) * t301 - mrSges(5,1) * t302 - pkin(4) * t285 - pkin(3) * t284 + (-t396 * t373 - t395 * t377) * t367 - t382 + t401 * t343 + t406 * qJDD(3) + qJ(4) * t385 + (qJ(4) * mrSges(5,2) + t400) * t344; t284; t382;];
tauJ = t1;
