% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:27
% EndTime: 2019-12-05 18:30:29
% DurationCPUTime: 2.16s
% Computational Cost: add. (39679->186), mult. (47778->233), div. (0->0), fcn. (25560->10), ass. (0->81)
t372 = qJD(1) + qJD(2);
t364 = qJD(4) + t372;
t376 = sin(qJ(5));
t396 = t364 * t376;
t380 = cos(qJ(5));
t395 = t364 * t380;
t379 = sin(qJ(1));
t383 = cos(qJ(1));
t361 = t383 * g(2) + t379 * g(3);
t358 = qJDD(1) * pkin(1) + t361;
t360 = t379 * g(2) - t383 * g(3);
t384 = qJD(1) ^ 2;
t359 = -t384 * pkin(1) + t360;
t378 = sin(qJ(2));
t382 = cos(qJ(2));
t343 = t382 * t358 - t378 * t359;
t371 = qJDD(1) + qJDD(2);
t341 = t371 * pkin(2) + t343;
t344 = t378 * t358 + t382 * t359;
t370 = t372 ^ 2;
t342 = -t370 * pkin(2) + t344;
t374 = sin(pkin(9));
t375 = cos(pkin(9));
t336 = t375 * t341 - t374 * t342;
t334 = t371 * pkin(3) + t336;
t337 = t374 * t341 + t375 * t342;
t335 = -t370 * pkin(3) + t337;
t377 = sin(qJ(4));
t381 = cos(qJ(4));
t331 = t377 * t334 + t381 * t335;
t362 = t364 ^ 2;
t363 = qJDD(4) + t371;
t329 = -t362 * pkin(4) + t363 * pkin(8) + t331;
t373 = -g(1) + qJDD(3);
t326 = -t376 * t329 + t380 * t373;
t350 = (-mrSges(6,1) * t380 + mrSges(6,2) * t376) * t364;
t394 = qJD(5) * t364;
t351 = t376 * t363 + t380 * t394;
t357 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t395;
t324 = m(6) * t326 + qJDD(5) * mrSges(6,1) - t351 * mrSges(6,3) + qJD(5) * t357 - t350 * t396;
t327 = t380 * t329 + t376 * t373;
t352 = t380 * t363 - t376 * t394;
t356 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t396;
t325 = m(6) * t327 - qJDD(5) * mrSges(6,2) + t352 * mrSges(6,3) - qJD(5) * t356 + t350 * t395;
t388 = -t376 * t324 + t380 * t325;
t315 = m(5) * t331 - t362 * mrSges(5,1) - t363 * mrSges(5,2) + t388;
t330 = t381 * t334 - t377 * t335;
t328 = -t363 * pkin(4) - t362 * pkin(8) - t330;
t385 = -m(6) * t328 + t352 * mrSges(6,1) - t351 * mrSges(6,2) - t356 * t396 + t357 * t395;
t320 = m(5) * t330 + t363 * mrSges(5,1) - t362 * mrSges(5,2) + t385;
t312 = t377 * t315 + t381 * t320;
t309 = m(4) * t336 + t371 * mrSges(4,1) - t370 * mrSges(4,2) + t312;
t389 = t381 * t315 - t377 * t320;
t310 = m(4) * t337 - t370 * mrSges(4,1) - t371 * mrSges(4,2) + t389;
t304 = t375 * t309 + t374 * t310;
t302 = m(3) * t343 + t371 * mrSges(3,1) - t370 * mrSges(3,2) + t304;
t390 = -t374 * t309 + t375 * t310;
t303 = m(3) * t344 - t370 * mrSges(3,1) - t371 * mrSges(3,2) + t390;
t296 = t382 * t302 + t378 * t303;
t316 = t380 * t324 + t376 * t325;
t393 = m(5) * t373 + t316;
t391 = -t378 * t302 + t382 * t303;
t294 = m(2) * t360 - t384 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t391;
t295 = m(2) * t361 + qJDD(1) * mrSges(2,1) - t384 * mrSges(2,2) + t296;
t392 = t383 * t294 - t379 * t295;
t387 = m(4) * t373 + t393;
t386 = -t379 * t294 - t383 * t295;
t347 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t376 + Ifges(6,4) * t380) * t364;
t346 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t376 + Ifges(6,2) * t380) * t364;
t345 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t376 + Ifges(6,6) * t380) * t364;
t318 = mrSges(6,2) * t328 - mrSges(6,3) * t326 + Ifges(6,1) * t351 + Ifges(6,4) * t352 + Ifges(6,5) * qJDD(5) - qJD(5) * t346 + t345 * t395;
t317 = -mrSges(6,1) * t328 + mrSges(6,3) * t327 + Ifges(6,4) * t351 + Ifges(6,2) * t352 + Ifges(6,6) * qJDD(5) + qJD(5) * t347 - t345 * t396;
t311 = -mrSges(5,1) * t373 - mrSges(6,1) * t326 + mrSges(6,2) * t327 + mrSges(5,3) * t331 + t362 * Ifges(5,5) - Ifges(6,5) * t351 + Ifges(5,6) * t363 - Ifges(6,6) * t352 - Ifges(6,3) * qJDD(5) - pkin(4) * t316 + (-t346 * t376 + t347 * t380) * t364;
t305 = mrSges(5,2) * t373 - mrSges(5,3) * t330 + Ifges(5,5) * t363 - t362 * Ifges(5,6) - pkin(8) * t316 - t376 * t317 + t380 * t318;
t298 = mrSges(4,2) * t373 - mrSges(4,3) * t336 + Ifges(4,5) * t371 - t370 * Ifges(4,6) - pkin(7) * t312 + t381 * t305 - t377 * t311;
t297 = -mrSges(4,1) * t373 + mrSges(4,3) * t337 + t370 * Ifges(4,5) + Ifges(4,6) * t371 - pkin(3) * t393 + pkin(7) * t389 + t377 * t305 + t381 * t311;
t292 = -mrSges(3,2) * g(1) - mrSges(3,3) * t343 + Ifges(3,5) * t371 - t370 * Ifges(3,6) - qJ(3) * t304 - t374 * t297 + t375 * t298;
t291 = mrSges(3,1) * g(1) + mrSges(3,3) * t344 + t370 * Ifges(3,5) + Ifges(3,6) * t371 - pkin(2) * t387 + qJ(3) * t390 + t375 * t297 + t374 * t298;
t290 = -mrSges(2,2) * g(1) - mrSges(2,3) * t361 + Ifges(2,5) * qJDD(1) - t384 * Ifges(2,6) - pkin(6) * t296 - t378 * t291 + t382 * t292;
t289 = Ifges(2,6) * qJDD(1) + t384 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t360 + t378 * t292 + t382 * t291 - pkin(1) * (-m(3) * g(1) + t387) + pkin(6) * t391;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t387; -m(1) * g(2) + t386; -m(1) * g(3) + t392; pkin(1) * t296 + pkin(2) * t304 - mrSges(3,2) * t344 + mrSges(3,1) * t343 + pkin(3) * t312 + mrSges(4,1) * t336 - mrSges(4,2) * t337 + pkin(8) * t388 + mrSges(5,1) * t330 - mrSges(5,2) * t331 + t376 * t318 + t380 * t317 + pkin(4) * t385 + mrSges(2,1) * t361 - mrSges(2,2) * t360 + Ifges(2,3) * qJDD(1) + Ifges(5,3) * t363 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(3,3) + Ifges(4,3)) * t371; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t392 - t383 * t289 - t379 * t290; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t386 - t379 * t289 + t383 * t290;];
tauB = t1;
