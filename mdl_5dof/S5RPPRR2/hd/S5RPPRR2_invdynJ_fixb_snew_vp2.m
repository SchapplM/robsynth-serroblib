% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:47
% EndTime: 2019-12-05 17:39:48
% DurationCPUTime: 1.10s
% Computational Cost: add. (6529->184), mult. (14706->236), div. (0->0), fcn. (9884->8), ass. (0->85)
t372 = qJD(1) ^ 2;
t368 = sin(qJ(1));
t371 = cos(qJ(1));
t388 = t368 * g(1) - t371 * g(2);
t378 = -t372 * qJ(2) + qJDD(2) - t388;
t397 = -pkin(1) - qJ(3);
t401 = -(2 * qJD(1) * qJD(3)) + t397 * qJDD(1) + t378;
t365 = cos(pkin(8));
t400 = t365 ^ 2;
t384 = -t371 * g(1) - t368 * g(2);
t399 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t384;
t398 = pkin(3) * t372;
t396 = t365 * mrSges(4,2);
t364 = sin(pkin(8));
t331 = t364 * g(3) + t401 * t365;
t317 = (-pkin(6) * qJDD(1) - t364 * t398) * t365 + t331;
t332 = -t365 * g(3) + t401 * t364;
t360 = t364 ^ 2;
t390 = t364 * qJDD(1);
t318 = -pkin(6) * t390 - t360 * t398 + t332;
t367 = sin(qJ(4));
t370 = cos(qJ(4));
t307 = t370 * t317 - t367 * t318;
t381 = -t364 * t367 + t365 * t370;
t382 = -t364 * t370 - t365 * t367;
t346 = t382 * qJD(1);
t392 = t346 * qJD(4);
t334 = t381 * qJDD(1) + t392;
t347 = t381 * qJD(1);
t297 = (-t334 + t392) * pkin(7) + (t346 * t347 + qJDD(4)) * pkin(4) + t307;
t308 = t367 * t317 + t370 * t318;
t333 = -t347 * qJD(4) + t382 * qJDD(1);
t341 = qJD(4) * pkin(4) - t347 * pkin(7);
t345 = t346 ^ 2;
t298 = -t345 * pkin(4) + t333 * pkin(7) - qJD(4) * t341 + t308;
t366 = sin(qJ(5));
t369 = cos(qJ(5));
t295 = t369 * t297 - t366 * t298;
t326 = t369 * t346 - t366 * t347;
t306 = t326 * qJD(5) + t366 * t333 + t369 * t334;
t327 = t366 * t346 + t369 * t347;
t313 = -t326 * mrSges(6,1) + t327 * mrSges(6,2);
t362 = qJD(4) + qJD(5);
t319 = -t362 * mrSges(6,2) + t326 * mrSges(6,3);
t359 = qJDD(4) + qJDD(5);
t292 = m(6) * t295 + t359 * mrSges(6,1) - t306 * mrSges(6,3) - t327 * t313 + t362 * t319;
t296 = t366 * t297 + t369 * t298;
t305 = -t327 * qJD(5) + t369 * t333 - t366 * t334;
t320 = t362 * mrSges(6,1) - t327 * mrSges(6,3);
t293 = m(6) * t296 - t359 * mrSges(6,2) + t305 * mrSges(6,3) + t326 * t313 - t362 * t320;
t286 = t369 * t292 + t366 * t293;
t329 = -t346 * mrSges(5,1) + t347 * mrSges(5,2);
t339 = -qJD(4) * mrSges(5,2) + t346 * mrSges(5,3);
t284 = m(5) * t307 + qJDD(4) * mrSges(5,1) - t334 * mrSges(5,3) + qJD(4) * t339 - t347 * t329 + t286;
t340 = qJD(4) * mrSges(5,1) - t347 * mrSges(5,3);
t385 = -t366 * t292 + t369 * t293;
t285 = m(5) * t308 - qJDD(4) * mrSges(5,2) + t333 * mrSges(5,3) - qJD(4) * t340 + t346 * t329 + t385;
t395 = t370 * t284 + t367 * t285;
t394 = -t360 - t400;
t387 = t394 * mrSges(4,3);
t386 = -t367 * t284 + t370 * t285;
t380 = -qJDD(1) * mrSges(4,3) - t372 * (t364 * mrSges(4,1) + t396);
t383 = t365 * (m(4) * t331 + t380 * t365 + t395) + t364 * (m(4) * t332 + t380 * t364 + t386);
t377 = qJDD(3) + t399;
t322 = pkin(3) * t390 + (t394 * pkin(6) + t397) * t372 + t377;
t300 = -t333 * pkin(4) - t345 * pkin(7) + t347 * t341 + t322;
t376 = m(6) * t300 - t305 * mrSges(6,1) + t306 * mrSges(6,2) - t326 * t319 + t327 * t320;
t310 = Ifges(6,4) * t327 + Ifges(6,2) * t326 + Ifges(6,6) * t362;
t311 = Ifges(6,1) * t327 + Ifges(6,4) * t326 + Ifges(6,5) * t362;
t375 = mrSges(6,1) * t295 - mrSges(6,2) * t296 + Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t359 + t327 * t310 - t326 * t311;
t374 = m(5) * t322 - t333 * mrSges(5,1) + t334 * mrSges(5,2) - t346 * t339 + t347 * t340 + t376;
t338 = t397 * t372 + t377;
t373 = m(4) * t338 + mrSges(4,1) * t390 + qJDD(1) * t396 + t374;
t344 = -qJDD(1) * pkin(1) + t378;
t343 = t372 * pkin(1) - t399;
t325 = Ifges(5,1) * t347 + Ifges(5,4) * t346 + Ifges(5,5) * qJD(4);
t324 = Ifges(5,4) * t347 + Ifges(5,2) * t346 + Ifges(5,6) * qJD(4);
t323 = Ifges(5,5) * t347 + Ifges(5,6) * t346 + Ifges(5,3) * qJD(4);
t309 = Ifges(6,5) * t327 + Ifges(6,6) * t326 + Ifges(6,3) * t362;
t288 = mrSges(6,2) * t300 - mrSges(6,3) * t295 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t359 + t326 * t309 - t362 * t310;
t287 = -mrSges(6,1) * t300 + mrSges(6,3) * t296 + Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t359 - t327 * t309 + t362 * t311;
t278 = mrSges(5,2) * t322 - mrSges(5,3) * t307 + Ifges(5,1) * t334 + Ifges(5,4) * t333 + Ifges(5,5) * qJDD(4) - pkin(7) * t286 - qJD(4) * t324 - t366 * t287 + t369 * t288 + t346 * t323;
t277 = -mrSges(5,1) * t322 + mrSges(5,3) * t308 + Ifges(5,4) * t334 + Ifges(5,2) * t333 + Ifges(5,6) * qJDD(4) - pkin(4) * t376 + pkin(7) * t385 + qJD(4) * t325 + t369 * t287 + t366 * t288 - t347 * t323;
t276 = m(3) * t344 + qJDD(1) * mrSges(3,2) - t372 * mrSges(3,3) + t383;
t1 = [mrSges(2,1) * t388 - mrSges(2,2) * t384 + mrSges(3,2) * t344 - mrSges(3,3) * t343 + t365 * (mrSges(4,2) * t338 - mrSges(4,3) * t331 - pkin(6) * t395 - t367 * t277 + t370 * t278) - t364 * (-mrSges(4,1) * t338 + mrSges(4,3) * t332 - pkin(3) * t374 + pkin(6) * t386 + t370 * t277 + t367 * t278) - qJ(3) * t383 - pkin(1) * t276 + qJ(2) * (t373 + (mrSges(3,2) + t387) * t372 - m(3) * t343) + (Ifges(4,1) * t400 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t365 + Ifges(4,2) * t364) * t364) * qJDD(1); t276; t372 * t387 + t373; mrSges(5,1) * t307 - mrSges(5,2) * t308 + Ifges(5,5) * t334 + Ifges(5,6) * t333 + Ifges(5,3) * qJDD(4) + pkin(4) * t286 + t347 * t324 - t346 * t325 + t375; t375;];
tauJ = t1;
