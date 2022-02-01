% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:40
% EndTime: 2022-01-23 09:22:41
% DurationCPUTime: 1.12s
% Computational Cost: add. (8911->218), mult. (19231->285), div. (0->0), fcn. (12248->10), ass. (0->90)
t369 = sin(qJ(1));
t372 = cos(qJ(1));
t382 = t369 * g(1) - t372 * g(2);
t347 = qJDD(1) * pkin(1) + t382;
t373 = qJD(1) ^ 2;
t378 = -t372 * g(1) - t369 * g(2);
t349 = -t373 * pkin(1) + t378;
t364 = sin(pkin(8));
t366 = cos(pkin(8));
t327 = t364 * t347 + t366 * t349;
t322 = -t373 * pkin(2) + qJDD(1) * pkin(6) + t327;
t362 = -g(3) + qJDD(2);
t368 = sin(qJ(3));
t371 = cos(qJ(3));
t311 = -t368 * t322 + t371 * t362;
t384 = qJD(1) * qJD(3);
t383 = t371 * t384;
t350 = t368 * qJDD(1) + t383;
t308 = (-t350 + t383) * qJ(4) + (t368 * t371 * t373 + qJDD(3)) * pkin(3) + t311;
t312 = t371 * t322 + t368 * t362;
t351 = t371 * qJDD(1) - t368 * t384;
t385 = t368 * qJD(1);
t352 = qJD(3) * pkin(3) - qJ(4) * t385;
t361 = t371 ^ 2;
t309 = -t361 * t373 * pkin(3) + t351 * qJ(4) - qJD(3) * t352 + t312;
t363 = sin(pkin(9));
t365 = cos(pkin(9));
t337 = (t371 * t363 + t368 * t365) * qJD(1);
t288 = -0.2e1 * qJD(4) * t337 + t365 * t308 - t363 * t309;
t329 = t365 * t350 + t363 * t351;
t336 = (-t368 * t363 + t371 * t365) * qJD(1);
t286 = (qJD(3) * t336 - t329) * pkin(7) + (t336 * t337 + qJDD(3)) * pkin(4) + t288;
t289 = 0.2e1 * qJD(4) * t336 + t363 * t308 + t365 * t309;
t328 = -t363 * t350 + t365 * t351;
t332 = qJD(3) * pkin(4) - t337 * pkin(7);
t335 = t336 ^ 2;
t287 = -t335 * pkin(4) + t328 * pkin(7) - qJD(3) * t332 + t289;
t367 = sin(qJ(5));
t370 = cos(qJ(5));
t284 = t370 * t286 - t367 * t287;
t319 = t370 * t336 - t367 * t337;
t298 = t319 * qJD(5) + t367 * t328 + t370 * t329;
t320 = t367 * t336 + t370 * t337;
t307 = -t319 * mrSges(6,1) + t320 * mrSges(6,2);
t360 = qJD(3) + qJD(5);
t313 = -t360 * mrSges(6,2) + t319 * mrSges(6,3);
t359 = qJDD(3) + qJDD(5);
t280 = m(6) * t284 + t359 * mrSges(6,1) - t298 * mrSges(6,3) - t320 * t307 + t360 * t313;
t285 = t367 * t286 + t370 * t287;
t297 = -t320 * qJD(5) + t370 * t328 - t367 * t329;
t314 = t360 * mrSges(6,1) - t320 * mrSges(6,3);
t281 = m(6) * t285 - t359 * mrSges(6,2) + t297 * mrSges(6,3) + t319 * t307 - t360 * t314;
t274 = t370 * t280 + t367 * t281;
t324 = -t336 * mrSges(5,1) + t337 * mrSges(5,2);
t330 = -qJD(3) * mrSges(5,2) + t336 * mrSges(5,3);
t272 = m(5) * t288 + qJDD(3) * mrSges(5,1) - t329 * mrSges(5,3) + qJD(3) * t330 - t337 * t324 + t274;
t331 = qJD(3) * mrSges(5,1) - t337 * mrSges(5,3);
t379 = -t367 * t280 + t370 * t281;
t273 = m(5) * t289 - qJDD(3) * mrSges(5,2) + t328 * mrSges(5,3) - qJD(3) * t331 + t336 * t324 + t379;
t268 = t365 * t272 + t363 * t273;
t386 = qJD(1) * t371;
t348 = (-t371 * mrSges(4,1) + t368 * mrSges(4,2)) * qJD(1);
t354 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t386;
t266 = m(4) * t311 + qJDD(3) * mrSges(4,1) - t350 * mrSges(4,3) + qJD(3) * t354 - t348 * t385 + t268;
t353 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t385;
t380 = -t363 * t272 + t365 * t273;
t267 = m(4) * t312 - qJDD(3) * mrSges(4,2) + t351 * mrSges(4,3) - qJD(3) * t353 + t348 * t386 + t380;
t381 = -t368 * t266 + t371 * t267;
t326 = t366 * t347 - t364 * t349;
t376 = -qJDD(1) * pkin(2) - t326;
t310 = -t351 * pkin(3) + qJDD(4) + t352 * t385 + (-qJ(4) * t361 - pkin(6)) * t373 + t376;
t291 = -t328 * pkin(4) - t335 * pkin(7) + t337 * t332 + t310;
t377 = m(6) * t291 - t297 * mrSges(6,1) + t298 * mrSges(6,2) - t319 * t313 + t320 * t314;
t300 = Ifges(6,4) * t320 + Ifges(6,2) * t319 + Ifges(6,6) * t360;
t301 = Ifges(6,1) * t320 + Ifges(6,4) * t319 + Ifges(6,5) * t360;
t375 = mrSges(6,1) * t284 - mrSges(6,2) * t285 + Ifges(6,5) * t298 + Ifges(6,6) * t297 + Ifges(6,3) * t359 + t320 * t300 - t319 * t301;
t282 = m(5) * t310 - t328 * mrSges(5,1) + t329 * mrSges(5,2) - t336 * t330 + t337 * t331 + t377;
t321 = -t373 * pkin(6) + t376;
t374 = -m(4) * t321 + t351 * mrSges(4,1) - t350 * mrSges(4,2) - t353 * t385 + t354 * t386 - t282;
t343 = Ifges(4,5) * qJD(3) + (t368 * Ifges(4,1) + t371 * Ifges(4,4)) * qJD(1);
t342 = Ifges(4,6) * qJD(3) + (t368 * Ifges(4,4) + t371 * Ifges(4,2)) * qJD(1);
t318 = Ifges(5,1) * t337 + Ifges(5,4) * t336 + Ifges(5,5) * qJD(3);
t317 = Ifges(5,4) * t337 + Ifges(5,2) * t336 + Ifges(5,6) * qJD(3);
t316 = Ifges(5,5) * t337 + Ifges(5,6) * t336 + Ifges(5,3) * qJD(3);
t299 = Ifges(6,5) * t320 + Ifges(6,6) * t319 + Ifges(6,3) * t360;
t276 = mrSges(6,2) * t291 - mrSges(6,3) * t284 + Ifges(6,1) * t298 + Ifges(6,4) * t297 + Ifges(6,5) * t359 + t319 * t299 - t360 * t300;
t275 = -mrSges(6,1) * t291 + mrSges(6,3) * t285 + Ifges(6,4) * t298 + Ifges(6,2) * t297 + Ifges(6,6) * t359 - t320 * t299 + t360 * t301;
t264 = mrSges(5,2) * t310 - mrSges(5,3) * t288 + Ifges(5,1) * t329 + Ifges(5,4) * t328 + Ifges(5,5) * qJDD(3) - pkin(7) * t274 - qJD(3) * t317 - t367 * t275 + t370 * t276 + t336 * t316;
t263 = -mrSges(5,1) * t310 + mrSges(5,3) * t289 + Ifges(5,4) * t329 + Ifges(5,2) * t328 + Ifges(5,6) * qJDD(3) - pkin(4) * t377 + pkin(7) * t379 + qJD(3) * t318 + t370 * t275 + t367 * t276 - t337 * t316;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t382 - mrSges(2,2) * t378 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t326 - mrSges(3,2) * t327 + t368 * (mrSges(4,2) * t321 - mrSges(4,3) * t311 + Ifges(4,1) * t350 + Ifges(4,4) * t351 + Ifges(4,5) * qJDD(3) - qJ(4) * t268 - qJD(3) * t342 - t363 * t263 + t365 * t264) + t371 * (-mrSges(4,1) * t321 + mrSges(4,3) * t312 + Ifges(4,4) * t350 + Ifges(4,2) * t351 + Ifges(4,6) * qJDD(3) - pkin(3) * t282 + qJ(4) * t380 + qJD(3) * t343 + t365 * t263 + t363 * t264) + pkin(2) * t374 + pkin(6) * t381 + pkin(1) * (t364 * (m(3) * t327 - t373 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t381) + t366 * (m(3) * t326 + qJDD(1) * mrSges(3,1) - t373 * mrSges(3,2) + t374)); m(3) * t362 + t371 * t266 + t368 * t267; (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t368 * t342 - t371 * t343) * qJD(1) + t375 + Ifges(4,5) * t350 + Ifges(4,6) * t351 - t336 * t318 + t337 * t317 + Ifges(5,5) * t329 + Ifges(5,6) * t328 + mrSges(4,1) * t311 - mrSges(4,2) * t312 + mrSges(5,1) * t288 - mrSges(5,2) * t289 + pkin(3) * t268 + pkin(4) * t274; t282; t375;];
tauJ = t1;
