% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:00
% EndTime: 2019-12-05 15:22:01
% DurationCPUTime: 0.89s
% Computational Cost: add. (3355->151), mult. (8150->220), div. (0->0), fcn. (5317->10), ass. (0->86)
t355 = qJD(2) ^ 2;
t347 = sin(pkin(7));
t350 = cos(pkin(7));
t335 = g(1) * t347 - g(2) * t350;
t336 = -g(1) * t350 - g(2) * t347;
t352 = sin(qJ(2));
t354 = cos(qJ(2));
t370 = t335 * t354 - t352 * t336;
t359 = -qJ(3) * t355 + qJDD(3) - t370;
t346 = sin(pkin(8));
t349 = cos(pkin(8));
t366 = -pkin(3) * t349 - qJ(4) * t346;
t378 = qJD(2) * t346;
t386 = (-pkin(2) + t366) * qJDD(2) + t359 - 0.2e1 * qJD(4) * t378;
t385 = 0.2e1 * t349;
t379 = t352 * t335 + t354 * t336;
t314 = -pkin(2) * t355 + qJDD(2) * qJ(3) + t379;
t344 = -g(3) + qJDD(1);
t376 = qJD(2) * qJD(3);
t304 = t344 * t349 + (-t314 - 0.2e1 * t376) * t346;
t384 = t349 ^ 2;
t383 = pkin(6) * t346;
t345 = sin(pkin(9));
t382 = t345 * t346;
t348 = cos(pkin(9));
t381 = t346 * t348;
t305 = t349 * t314 + t346 * t344 + t376 * t385;
t330 = t366 * qJD(2);
t377 = qJD(2) * t349;
t297 = t330 * t377 + t305;
t343 = t346 ^ 2;
t364 = -pkin(4) * t349 - pkin(6) * t381;
t380 = t386 * t348;
t289 = t364 * qJDD(2) + (-t297 + (-pkin(4) * t343 * t348 + t349 * t383) * t355) * t345 + t380;
t292 = t348 * t297 + t386 * t345;
t326 = t364 * qJD(2);
t374 = t345 ^ 2 * t343 * t355;
t375 = qJDD(2) * t345;
t290 = -pkin(4) * t374 + t326 * t377 - t375 * t383 + t292;
t351 = sin(qJ(5));
t353 = cos(qJ(5));
t287 = t289 * t353 - t290 * t351;
t361 = (-t345 * t353 - t348 * t351) * t346;
t320 = qJD(2) * t361;
t360 = (-t345 * t351 + t348 * t353) * t346;
t321 = qJD(2) * t360;
t306 = -mrSges(6,1) * t320 + mrSges(6,2) * t321;
t310 = qJD(5) * t320 + qJDD(2) * t360;
t338 = qJD(5) - t377;
t315 = -mrSges(6,2) * t338 + mrSges(6,3) * t320;
t337 = -qJDD(2) * t349 + qJDD(5);
t285 = m(6) * t287 + mrSges(6,1) * t337 - mrSges(6,3) * t310 - t306 * t321 + t315 * t338;
t288 = t289 * t351 + t290 * t353;
t309 = -qJD(5) * t321 + qJDD(2) * t361;
t316 = mrSges(6,1) * t338 - mrSges(6,3) * t321;
t286 = m(6) * t288 - mrSges(6,2) * t337 + mrSges(6,3) * t309 + t306 * t320 - t316 * t338;
t278 = t353 * t285 + t351 * t286;
t372 = -t285 * t351 + t353 * t286;
t368 = -mrSges(4,1) * t349 + mrSges(4,2) * t346;
t367 = t345 * mrSges(5,1) + t348 * mrSges(5,2);
t291 = -t297 * t345 + t380;
t322 = t367 * t378;
t362 = mrSges(5,2) * t349 - mrSges(5,3) * t382;
t324 = t362 * qJD(2);
t363 = -mrSges(5,1) * t349 - mrSges(5,3) * t381;
t276 = m(5) * t291 + t363 * qJDD(2) + (-t322 * t381 - t324 * t349) * qJD(2) + t278;
t325 = t363 * qJD(2);
t277 = m(5) * t292 + t362 * qJDD(2) + (-t322 * t382 + t325 * t349) * qJD(2) + t372;
t275 = t276 * t348 + t277 * t345;
t365 = t324 * t345 + t325 * t348;
t296 = t330 * t378 + qJDD(4) - t304;
t294 = -pkin(6) * t374 + (qJD(2) * t326 * t348 + pkin(4) * t375) * t346 + t296;
t358 = -m(6) * t294 + mrSges(6,1) * t309 - t310 * mrSges(6,2) + t315 * t320 - t321 * t316;
t357 = m(5) * t296 - t358;
t299 = Ifges(6,4) * t321 + Ifges(6,2) * t320 + Ifges(6,6) * t338;
t300 = Ifges(6,1) * t321 + Ifges(6,4) * t320 + Ifges(6,5) * t338;
t356 = mrSges(6,1) * t287 - mrSges(6,2) * t288 + Ifges(6,5) * t310 + Ifges(6,6) * t309 + Ifges(6,3) * t337 + t321 * t299 - t320 * t300;
t331 = t368 * qJD(2);
t313 = -qJDD(2) * pkin(2) + t359;
t298 = Ifges(6,5) * t321 + Ifges(6,6) * t320 + Ifges(6,3) * t338;
t281 = m(4) * t304 + ((-mrSges(4,3) - t367) * qJDD(2) + (-t331 - t365) * qJD(2)) * t346 - t357;
t280 = mrSges(6,2) * t294 - mrSges(6,3) * t287 + Ifges(6,1) * t310 + Ifges(6,4) * t309 + Ifges(6,5) * t337 + t298 * t320 - t299 * t338;
t279 = -mrSges(6,1) * t294 + mrSges(6,3) * t288 + Ifges(6,4) * t310 + Ifges(6,2) * t309 + Ifges(6,6) * t337 - t298 * t321 + t300 * t338;
t274 = m(4) * t313 + t368 * qJDD(2) + (-t343 - t384) * t355 * mrSges(4,3) + t275;
t273 = m(4) * t305 - t276 * t345 + t277 * t348 + (qJDD(2) * mrSges(4,3) + qJD(2) * t331) * t349;
t1 = [t273 * t346 + t281 * t349 + (m(2) + m(3)) * t344; mrSges(3,1) * t370 - mrSges(3,2) * t379 - pkin(2) * t274 + (-mrSges(4,1) * t313 - mrSges(5,1) * t291 + mrSges(5,2) * t292 + mrSges(4,3) * t305 - pkin(3) * t275 - pkin(4) * t278 + qJ(3) * t273 - t356) * t349 + (mrSges(4,2) * t313 - mrSges(4,3) * t304 + t348 * (mrSges(5,2) * t296 - mrSges(5,3) * t291 - pkin(6) * t278 - t351 * t279 + t353 * t280) - t345 * (-mrSges(5,1) * t296 + mrSges(5,3) * t292 + pkin(4) * t358 + pkin(6) * t372 + t353 * t279 + t351 * t280) - qJ(4) * t275 - qJ(3) * t281) * t346 + (Ifges(3,3) + (Ifges(4,2) + Ifges(5,3)) * t384 + ((Ifges(5,1) * t348 ^ 2 + Ifges(4,1) + (-0.2e1 * Ifges(5,4) * t348 + Ifges(5,2) * t345) * t345) * t346 + (-Ifges(5,5) * t348 + Ifges(5,6) * t345 + Ifges(4,4)) * t385) * t346) * qJDD(2); t274; (t365 * qJD(2) + t367 * qJDD(2)) * t346 + t357; t356;];
tauJ = t1;
