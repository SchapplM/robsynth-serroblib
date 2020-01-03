% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP3
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:14
% EndTime: 2019-12-31 20:53:15
% DurationCPUTime: 0.64s
% Computational Cost: add. (2873->179), mult. (3601->208), div. (0->0), fcn. (1541->6), ass. (0->76)
t401 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t389 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t388 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t400 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t387 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t399 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t362 = sin(qJ(3));
t357 = qJD(1) + qJD(2);
t390 = qJD(3) * t357;
t382 = t362 * t390;
t393 = t357 * t362;
t397 = -2 * qJD(4);
t398 = pkin(3) * t382 + t393 * t397;
t396 = -2 * qJD(5);
t395 = -mrSges(6,1) - mrSges(4,3);
t356 = qJDD(1) + qJDD(2);
t365 = cos(qJ(3));
t383 = t365 * t390;
t332 = t356 * t362 + t383;
t394 = mrSges(6,1) * t332;
t392 = t357 * t365;
t364 = sin(qJ(1));
t367 = cos(qJ(1));
t381 = t364 * g(1) - g(2) * t367;
t340 = qJDD(1) * pkin(1) + t381;
t377 = -g(1) * t367 - g(2) * t364;
t341 = -qJD(1) ^ 2 * pkin(1) + t377;
t363 = sin(qJ(2));
t366 = cos(qJ(2));
t306 = t363 * t340 + t366 * t341;
t355 = t357 ^ 2;
t303 = -pkin(2) * t355 + pkin(7) * t356 + t306;
t298 = -t365 * g(3) - t362 * t303;
t329 = (mrSges(5,2) * t365 - mrSges(5,3) * t362) * t357;
t331 = (-mrSges(6,2) * t362 - mrSges(6,3) * t365) * t357;
t391 = t329 + t331;
t386 = (t362 * t388 + t365 * t387) * t357 + t399 * qJD(3);
t385 = (t362 * t389 + t400 * t365) * t357 + t387 * qJD(3);
t384 = (t401 * t362 + t365 * t389) * t357 + t388 * qJD(3);
t299 = -g(3) * t362 + t365 * t303;
t330 = (-mrSges(4,1) * t365 + mrSges(4,2) * t362) * t357;
t333 = t356 * t365 - t382;
t342 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t393;
t343 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t392;
t346 = -mrSges(5,1) * t392 - qJD(3) * mrSges(5,3);
t328 = (-pkin(3) * t365 - qJ(4) * t362) * t357;
t368 = qJD(3) ^ 2;
t370 = -pkin(3) * t368 + qJDD(3) * qJ(4) + t328 * t392 + t299;
t296 = qJD(3) * t397 - t370;
t348 = mrSges(5,1) * t393 + qJD(3) * mrSges(5,2);
t344 = pkin(4) * t393 - qJD(3) * qJ(5);
t361 = t365 ^ 2;
t294 = -qJ(5) * t355 * t361 + pkin(4) * t333 + qJDD(5) + ((2 * qJD(4)) + t344) * qJD(3) + t370;
t345 = mrSges(6,1) * t393 - qJD(3) * mrSges(6,3);
t379 = m(6) * t294 + qJDD(3) * mrSges(6,2) + qJD(3) * t345 + t331 * t392;
t371 = -m(5) * t296 + qJDD(3) * mrSges(5,3) + qJD(3) * t348 + t329 * t392 + t379;
t297 = -qJDD(3) * pkin(3) - qJ(4) * t368 + t328 * t393 + qJDD(4) - t298;
t293 = qJD(3) * t396 + (-t355 * t362 * t365 - qJDD(3)) * qJ(5) + (t332 - t383) * pkin(4) + t297;
t347 = mrSges(6,1) * t392 + qJD(3) * mrSges(6,2);
t378 = -m(6) * t293 + qJDD(3) * mrSges(6,3) + qJD(3) * t347;
t374 = m(5) * t297 + t332 * mrSges(5,1) - t378;
t380 = -(m(4) * t298 + t395 * t332 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t343 - t346) * qJD(3) + (-t330 - t391) * t393 - t374) * t362 + t365 * (t330 * t392 + m(4) * t299 - qJDD(3) * mrSges(4,2) - qJD(3) * t342 + (mrSges(5,1) - t395) * t333 + t371);
t305 = t340 * t366 - t363 * t341;
t375 = -pkin(2) * t356 - t305;
t291 = -qJ(4) * t332 + (-pkin(4) * t361 - pkin(7)) * t355 + (-pkin(3) - qJ(5)) * t333 + (-t344 * t362 + (-qJ(4) * qJD(3) + t396) * t365) * t357 + t375 + t398;
t376 = m(6) * t291 - t332 * mrSges(6,2) - t333 * mrSges(6,3) - t345 * t393 - t347 * t392;
t302 = -pkin(7) * t355 + t375;
t295 = -pkin(3) * t333 + (-t332 - t383) * qJ(4) + t302 + t398;
t372 = -m(5) * t295 - t333 * mrSges(5,2) + t348 * t393 - t376;
t286 = -mrSges(5,3) * t332 + t346 * t392 - t372;
t288 = t331 * t393 - t378 + t394;
t289 = mrSges(6,1) * t333 + t379;
t369 = -m(4) * t302 + t343 * t392 + t333 * mrSges(4,1) + (-t342 * t362 - t346 * t365) * t357 + (-mrSges(4,2) + mrSges(5,3)) * t332 + t372;
t373 = -mrSges(3,2) * t306 + t365 * (-mrSges(4,1) * t302 - mrSges(5,1) * t296 + mrSges(6,1) * t294 + mrSges(5,2) * t295 + mrSges(4,3) * t299 - mrSges(6,3) * t291 - pkin(3) * t286 + pkin(4) * t289 - qJ(5) * t376 + t384 * qJD(3) + t387 * qJDD(3) + t389 * t332 + t400 * t333 - t386 * t393) + t362 * (mrSges(5,1) * t297 + mrSges(6,1) * t293 + mrSges(4,2) * t302 - mrSges(6,2) * t291 - mrSges(4,3) * t298 - mrSges(5,3) * t295 + pkin(4) * t288 - qJ(4) * t286 - t385 * qJD(3) + t388 * qJDD(3) + t401 * t332 + t389 * t333 + t386 * t392) + pkin(7) * t380 + pkin(2) * t369 + mrSges(3,1) * t305 + Ifges(3,3) * t356;
t287 = qJDD(3) * mrSges(5,2) + qJD(3) * t346 + t391 * t393 + t374 + t394;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t381 - mrSges(2,2) * t377 + pkin(1) * (t363 * (m(3) * t306 - mrSges(3,1) * t355 - mrSges(3,2) * t356 + t380) + t366 * (m(3) * t305 + mrSges(3,1) * t356 - mrSges(3,2) * t355 + t369)) + t373; t373; mrSges(4,1) * t298 - mrSges(4,2) * t299 + mrSges(5,2) * t297 - mrSges(5,3) * t296 + mrSges(6,2) * t294 - mrSges(6,3) * t293 - qJ(5) * t288 - pkin(3) * t287 + qJ(4) * t371 + t388 * t332 + t399 * qJDD(3) + (qJ(4) * (mrSges(5,1) + mrSges(6,1)) + t387) * t333 + (t385 * t362 - t384 * t365) * t357; t287; t289;];
tauJ = t1;
