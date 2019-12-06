% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:34
% EndTime: 2019-12-05 15:16:36
% DurationCPUTime: 1.84s
% Computational Cost: add. (20762->206), mult. (37661->269), div. (0->0), fcn. (24678->10), ass. (0->89)
t371 = sin(pkin(8));
t373 = cos(pkin(8));
t361 = -t373 * g(1) - t371 * g(2);
t369 = -g(3) + qJDD(1);
t370 = sin(pkin(9));
t372 = cos(pkin(9));
t347 = t372 * t361 + t370 * t369;
t360 = t371 * g(1) - t373 * g(2);
t359 = qJDD(2) - t360;
t376 = sin(qJ(3));
t379 = cos(qJ(3));
t337 = t379 * t347 + t376 * t359;
t380 = qJD(3) ^ 2;
t332 = -t380 * pkin(3) + qJDD(3) * pkin(6) + t337;
t346 = t370 * t361 - t372 * t369;
t375 = sin(qJ(4));
t378 = cos(qJ(4));
t325 = -t375 * t332 + t378 * t346;
t391 = qJD(3) * qJD(4);
t390 = t378 * t391;
t357 = t375 * qJDD(3) + t390;
t322 = (-t357 + t390) * pkin(7) + (t375 * t378 * t380 + qJDD(4)) * pkin(4) + t325;
t326 = t378 * t332 + t375 * t346;
t358 = t378 * qJDD(3) - t375 * t391;
t393 = qJD(3) * t375;
t364 = qJD(4) * pkin(4) - pkin(7) * t393;
t368 = t378 ^ 2;
t323 = -t368 * t380 * pkin(4) + t358 * pkin(7) - qJD(4) * t364 + t326;
t374 = sin(qJ(5));
t377 = cos(qJ(5));
t320 = t377 * t322 - t374 * t323;
t351 = (-t374 * t375 + t377 * t378) * qJD(3);
t330 = t351 * qJD(5) + t377 * t357 + t374 * t358;
t352 = (t374 * t378 + t375 * t377) * qJD(3);
t339 = -t351 * mrSges(6,1) + t352 * mrSges(6,2);
t367 = qJD(4) + qJD(5);
t344 = -t367 * mrSges(6,2) + t351 * mrSges(6,3);
t366 = qJDD(4) + qJDD(5);
t318 = m(6) * t320 + t366 * mrSges(6,1) - t330 * mrSges(6,3) - t352 * t339 + t367 * t344;
t321 = t374 * t322 + t377 * t323;
t329 = -t352 * qJD(5) - t374 * t357 + t377 * t358;
t345 = t367 * mrSges(6,1) - t352 * mrSges(6,3);
t319 = m(6) * t321 - t366 * mrSges(6,2) + t329 * mrSges(6,3) + t351 * t339 - t367 * t345;
t311 = t377 * t318 + t374 * t319;
t356 = (-mrSges(5,1) * t378 + mrSges(5,2) * t375) * qJD(3);
t392 = qJD(3) * t378;
t363 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t392;
t309 = m(5) * t325 + qJDD(4) * mrSges(5,1) - t357 * mrSges(5,3) + qJD(4) * t363 - t356 * t393 + t311;
t362 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t393;
t385 = -t374 * t318 + t377 * t319;
t310 = m(5) * t326 - qJDD(4) * mrSges(5,2) + t358 * mrSges(5,3) - qJD(4) * t362 + t356 * t392 + t385;
t386 = -t375 * t309 + t378 * t310;
t305 = m(4) * t337 - t380 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t386;
t336 = -t376 * t347 + t379 * t359;
t384 = -qJDD(3) * pkin(3) - t336;
t331 = -t380 * pkin(6) + t384;
t324 = t364 * t393 - t358 * pkin(4) + (-pkin(7) * t368 - pkin(6)) * t380 + t384;
t382 = m(6) * t324 - t329 * mrSges(6,1) + t330 * mrSges(6,2) - t351 * t344 + t352 * t345;
t381 = -m(5) * t331 + t358 * mrSges(5,1) - t357 * mrSges(5,2) - t362 * t393 + t363 * t392 - t382;
t314 = m(4) * t336 + qJDD(3) * mrSges(4,1) - t380 * mrSges(4,2) + t381;
t387 = t379 * t305 - t376 * t314;
t300 = m(3) * t347 + t387;
t307 = t378 * t309 + t375 * t310;
t306 = (-m(3) - m(4)) * t346 - t307;
t388 = t372 * t300 - t370 * t306;
t292 = m(2) * t361 + t388;
t301 = t376 * t305 + t379 * t314;
t383 = -m(3) * t359 - t301;
t299 = m(2) * t360 + t383;
t394 = t371 * t292 + t373 * t299;
t293 = t370 * t300 + t372 * t306;
t389 = t373 * t292 - t371 * t299;
t350 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t375 + Ifges(5,4) * t378) * qJD(3);
t349 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t375 + Ifges(5,2) * t378) * qJD(3);
t348 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t375 + Ifges(5,6) * t378) * qJD(3);
t335 = Ifges(6,1) * t352 + Ifges(6,4) * t351 + Ifges(6,5) * t367;
t334 = Ifges(6,4) * t352 + Ifges(6,2) * t351 + Ifges(6,6) * t367;
t333 = Ifges(6,5) * t352 + Ifges(6,6) * t351 + Ifges(6,3) * t367;
t313 = mrSges(6,2) * t324 - mrSges(6,3) * t320 + Ifges(6,1) * t330 + Ifges(6,4) * t329 + Ifges(6,5) * t366 + t351 * t333 - t367 * t334;
t312 = -mrSges(6,1) * t324 + mrSges(6,3) * t321 + Ifges(6,4) * t330 + Ifges(6,2) * t329 + Ifges(6,6) * t366 - t352 * t333 + t367 * t335;
t302 = mrSges(5,2) * t331 - mrSges(5,3) * t325 + Ifges(5,1) * t357 + Ifges(5,4) * t358 + Ifges(5,5) * qJDD(4) - pkin(7) * t311 - qJD(4) * t349 - t374 * t312 + t377 * t313 + t348 * t392;
t295 = -mrSges(5,1) * t331 + mrSges(5,3) * t326 + Ifges(5,4) * t357 + Ifges(5,2) * t358 + Ifges(5,6) * qJDD(4) - pkin(4) * t382 + pkin(7) * t385 + qJD(4) * t350 + t377 * t312 + t374 * t313 - t348 * t393;
t294 = Ifges(4,6) * qJDD(3) + t380 * Ifges(4,5) - mrSges(4,1) * t346 + mrSges(4,3) * t337 - Ifges(5,5) * t357 - Ifges(5,6) * t358 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t325 + mrSges(5,2) * t326 - Ifges(6,5) * t330 - Ifges(6,6) * t329 - Ifges(6,3) * t366 - t352 * t334 + t351 * t335 - mrSges(6,1) * t320 + mrSges(6,2) * t321 - pkin(4) * t311 - pkin(3) * t307 + (-t375 * t349 + t378 * t350) * qJD(3);
t289 = mrSges(4,2) * t346 - mrSges(4,3) * t336 + Ifges(4,5) * qJDD(3) - t380 * Ifges(4,6) - pkin(6) * t307 - t375 * t295 + t378 * t302;
t288 = -mrSges(3,1) * t359 - mrSges(4,1) * t336 + mrSges(4,2) * t337 + mrSges(3,3) * t347 - Ifges(4,3) * qJDD(3) - pkin(2) * t301 - pkin(3) * t381 - pkin(6) * t386 - t378 * t295 - t375 * t302;
t287 = mrSges(3,2) * t359 + mrSges(3,3) * t346 - pkin(5) * t301 + t379 * t289 - t376 * t294;
t286 = -mrSges(2,1) * t369 + mrSges(2,3) * t361 + mrSges(3,1) * t346 + mrSges(3,2) * t347 - t376 * t289 - t379 * t294 - pkin(2) * (-m(4) * t346 - t307) - pkin(5) * t387 - pkin(1) * t293;
t285 = mrSges(2,2) * t369 - mrSges(2,3) * t360 - qJ(2) * t293 + t372 * t287 - t370 * t288;
t1 = [-m(1) * g(1) + t389; -m(1) * g(2) + t394; -m(1) * g(3) + m(2) * t369 + t293; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t394 + t373 * t285 - t371 * t286; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t389 + t371 * t285 + t373 * t286; -mrSges(1,1) * g(2) + mrSges(2,1) * t360 + mrSges(1,2) * g(1) - mrSges(2,2) * t361 + pkin(1) * t383 + qJ(2) * t388 + t370 * t287 + t372 * t288;];
tauB = t1;
