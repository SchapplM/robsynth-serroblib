% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:53
% EndTime: 2019-12-05 15:42:54
% DurationCPUTime: 0.83s
% Computational Cost: add. (5849->176), mult. (13383->231), div. (0->0), fcn. (9566->10), ass. (0->86)
t355 = qJD(2) ^ 2;
t347 = cos(pkin(9));
t342 = t347 ^ 2;
t375 = 0.2e1 * t347;
t374 = pkin(3) * t355;
t373 = pkin(6) * qJDD(2);
t346 = sin(pkin(8));
t348 = cos(pkin(8));
t333 = t346 * g(1) - t348 * g(2);
t334 = -t348 * g(1) - t346 * g(2);
t351 = sin(qJ(2));
t354 = cos(qJ(2));
t371 = t351 * t333 + t354 * t334;
t321 = -t355 * pkin(2) + qJDD(2) * qJ(3) + t371;
t345 = sin(pkin(9));
t344 = -g(3) + qJDD(1);
t367 = qJD(2) * qJD(3);
t370 = t347 * t344 - 0.2e1 * t345 * t367;
t304 = (t347 * t374 - t321 - t373) * t345 + t370;
t308 = t347 * t321 + t345 * t344 + t367 * t375;
t305 = -t342 * t374 + t347 * t373 + t308;
t350 = sin(qJ(4));
t353 = cos(qJ(4));
t288 = t353 * t304 - t350 * t305;
t361 = t345 * t353 + t347 * t350;
t360 = -t345 * t350 + t347 * t353;
t326 = t360 * qJD(2);
t368 = t326 * qJD(4);
t320 = t361 * qJDD(2) + t368;
t327 = t361 * qJD(2);
t284 = (-t320 + t368) * pkin(7) + (t326 * t327 + qJDD(4)) * pkin(4) + t288;
t289 = t350 * t304 + t353 * t305;
t319 = -t327 * qJD(4) + t360 * qJDD(2);
t324 = qJD(4) * pkin(4) - t327 * pkin(7);
t325 = t326 ^ 2;
t285 = -t325 * pkin(4) + t319 * pkin(7) - qJD(4) * t324 + t289;
t349 = sin(qJ(5));
t352 = cos(qJ(5));
t282 = t352 * t284 - t349 * t285;
t314 = t352 * t326 - t349 * t327;
t295 = t314 * qJD(5) + t349 * t319 + t352 * t320;
t315 = t349 * t326 + t352 * t327;
t300 = -t314 * mrSges(6,1) + t315 * mrSges(6,2);
t343 = qJD(4) + qJD(5);
t309 = -t343 * mrSges(6,2) + t314 * mrSges(6,3);
t340 = qJDD(4) + qJDD(5);
t279 = m(6) * t282 + t340 * mrSges(6,1) - t295 * mrSges(6,3) - t315 * t300 + t343 * t309;
t283 = t349 * t284 + t352 * t285;
t294 = -t315 * qJD(5) + t352 * t319 - t349 * t320;
t310 = t343 * mrSges(6,1) - t315 * mrSges(6,3);
t280 = m(6) * t283 - t340 * mrSges(6,2) + t294 * mrSges(6,3) + t314 * t300 - t343 * t310;
t272 = t352 * t279 + t349 * t280;
t317 = -t326 * mrSges(5,1) + t327 * mrSges(5,2);
t322 = -qJD(4) * mrSges(5,2) + t326 * mrSges(5,3);
t270 = m(5) * t288 + qJDD(4) * mrSges(5,1) - t320 * mrSges(5,3) + qJD(4) * t322 - t327 * t317 + t272;
t323 = qJD(4) * mrSges(5,1) - t327 * mrSges(5,3);
t365 = -t349 * t279 + t352 * t280;
t271 = m(5) * t289 - qJDD(4) * mrSges(5,2) + t319 * mrSges(5,3) - qJD(4) * t323 + t326 * t317 + t365;
t372 = t353 * t270 + t350 * t271;
t369 = -t345 ^ 2 - t342;
t366 = -t350 * t270 + t353 * t271;
t364 = t354 * t333 - t351 * t334;
t363 = -t347 * mrSges(4,1) + t345 * mrSges(4,2);
t362 = qJDD(3) - t364;
t359 = mrSges(4,3) * qJDD(2) + t355 * t363;
t306 = (-pkin(3) * t347 - pkin(2)) * qJDD(2) + (t369 * pkin(6) - qJ(3)) * t355 + t362;
t287 = -t319 * pkin(4) - t325 * pkin(7) + t327 * t324 + t306;
t358 = m(6) * t287 - t294 * mrSges(6,1) + t295 * mrSges(6,2) - t314 * t309 + t315 * t310;
t297 = Ifges(6,4) * t315 + Ifges(6,2) * t314 + Ifges(6,6) * t343;
t298 = Ifges(6,1) * t315 + Ifges(6,4) * t314 + Ifges(6,5) * t343;
t357 = mrSges(6,1) * t282 - mrSges(6,2) * t283 + Ifges(6,5) * t295 + Ifges(6,6) * t294 + Ifges(6,3) * t340 + t315 * t297 - t314 * t298;
t356 = m(5) * t306 - t319 * mrSges(5,1) + t320 * mrSges(5,2) - t326 * t322 + t327 * t323 + t358;
t318 = -qJDD(2) * pkin(2) - t355 * qJ(3) + t362;
t313 = Ifges(5,1) * t327 + Ifges(5,4) * t326 + Ifges(5,5) * qJD(4);
t312 = Ifges(5,4) * t327 + Ifges(5,2) * t326 + Ifges(5,6) * qJD(4);
t311 = Ifges(5,5) * t327 + Ifges(5,6) * t326 + Ifges(5,3) * qJD(4);
t307 = -t345 * t321 + t370;
t296 = Ifges(6,5) * t315 + Ifges(6,6) * t314 + Ifges(6,3) * t343;
t275 = t369 * t355 * mrSges(4,3) + m(4) * t318 + t363 * qJDD(2) + t356;
t274 = mrSges(6,2) * t287 - mrSges(6,3) * t282 + Ifges(6,1) * t295 + Ifges(6,4) * t294 + Ifges(6,5) * t340 + t314 * t296 - t343 * t297;
t273 = -mrSges(6,1) * t287 + mrSges(6,3) * t283 + Ifges(6,4) * t295 + Ifges(6,2) * t294 + Ifges(6,6) * t340 - t315 * t296 + t343 * t298;
t266 = m(4) * t308 + t359 * t347 + t366;
t265 = m(4) * t307 - t359 * t345 + t372;
t264 = mrSges(5,2) * t306 - mrSges(5,3) * t288 + Ifges(5,1) * t320 + Ifges(5,4) * t319 + Ifges(5,5) * qJDD(4) - pkin(7) * t272 - qJD(4) * t312 - t349 * t273 + t352 * t274 + t326 * t311;
t263 = -mrSges(5,1) * t306 + mrSges(5,3) * t289 + Ifges(5,4) * t320 + Ifges(5,2) * t319 + Ifges(5,6) * qJDD(4) - pkin(4) * t358 + pkin(7) * t365 + qJD(4) * t313 + t352 * t273 + t349 * t274 - t327 * t311;
t1 = [t347 * t265 + t345 * t266 + (m(2) + m(3)) * t344; mrSges(3,1) * t364 - mrSges(3,2) * t371 + t345 * (mrSges(4,2) * t318 - mrSges(4,3) * t307 - pkin(6) * t372 - t350 * t263 + t353 * t264) + t347 * (-mrSges(4,1) * t318 + mrSges(4,3) * t308 - pkin(3) * t356 + pkin(6) * t366 + t353 * t263 + t350 * t264) - pkin(2) * t275 + qJ(3) * (-t345 * t265 + t347 * t266) + (Ifges(4,2) * t342 + Ifges(3,3) + (Ifges(4,1) * t345 + Ifges(4,4) * t375) * t345) * qJDD(2); t275; mrSges(5,1) * t288 - mrSges(5,2) * t289 + Ifges(5,5) * t320 + Ifges(5,6) * t319 + Ifges(5,3) * qJDD(4) + pkin(4) * t272 + t327 * t312 - t326 * t313 + t357; t357;];
tauJ = t1;
