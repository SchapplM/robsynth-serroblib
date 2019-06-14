% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPPRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:32:52
% EndTime: 2019-05-04 19:32:54
% DurationCPUTime: 1.79s
% Computational Cost: add. (21482->167), mult. (36395->229), div. (0->0), fcn. (32318->18), ass. (0->92)
t342 = sin(pkin(12));
t348 = cos(pkin(12));
t335 = -g(1) * t348 - g(2) * t342;
t341 = sin(pkin(13));
t347 = cos(pkin(13));
t334 = g(1) * t342 - g(2) * t348;
t339 = -g(3) + qJDD(1);
t345 = sin(pkin(6));
t351 = cos(pkin(6));
t363 = t334 * t351 + t339 * t345;
t310 = t347 * t335 + t363 * t341;
t340 = sin(pkin(14));
t346 = cos(pkin(14));
t309 = -t341 * t335 + t363 * t347;
t321 = -t334 * t345 + t339 * t351 + qJDD(2);
t344 = sin(pkin(7));
t350 = cos(pkin(7));
t364 = t309 * t350 + t321 * t344;
t305 = -t340 * t310 + t364 * t346;
t308 = -t309 * t344 + t321 * t350 + qJDD(3);
t343 = sin(pkin(8));
t349 = cos(pkin(8));
t378 = t305 * t349 + t308 * t343;
t306 = t346 * t310 + t364 * t340;
t354 = sin(qJ(4));
t357 = cos(qJ(4));
t298 = -t354 * t306 + t378 * t357;
t301 = -t305 * t343 + t308 * t349;
t356 = cos(qJ(5));
t377 = t301 * t356;
t299 = t357 * t306 + t378 * t354;
t359 = qJD(4) ^ 2;
t297 = -pkin(4) * t359 + qJDD(4) * pkin(10) + t299;
t353 = sin(qJ(5));
t293 = t356 * t297 + t353 * t301;
t374 = qJD(4) * t353;
t373 = qJD(4) * t356;
t372 = qJD(4) * qJD(5);
t371 = t353 * t372;
t370 = t356 * t372;
t330 = (-mrSges(6,1) * t356 + mrSges(6,2) * t353) * qJD(4);
t333 = qJDD(4) * t356 - t371;
t336 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t374;
t331 = (-pkin(5) * t356 - pkin(11) * t353) * qJD(4);
t358 = qJD(5) ^ 2;
t291 = -pkin(5) * t358 + qJDD(5) * pkin(11) + t331 * t373 + t293;
t296 = -qJDD(4) * pkin(4) - pkin(10) * t359 - t298;
t332 = qJDD(4) * t353 + t370;
t294 = (-t332 - t370) * pkin(11) + (-t333 + t371) * pkin(5) + t296;
t352 = sin(qJ(6));
t355 = cos(qJ(6));
t288 = -t291 * t352 + t294 * t355;
t328 = qJD(5) * t355 - t352 * t374;
t317 = qJD(6) * t328 + qJDD(5) * t352 + t332 * t355;
t329 = qJD(5) * t352 + t355 * t374;
t318 = -mrSges(7,1) * t328 + mrSges(7,2) * t329;
t338 = qJD(6) - t373;
t319 = -mrSges(7,2) * t338 + mrSges(7,3) * t328;
t327 = qJDD(6) - t333;
t286 = m(7) * t288 + mrSges(7,1) * t327 - mrSges(7,3) * t317 - t318 * t329 + t319 * t338;
t289 = t291 * t355 + t294 * t352;
t316 = -qJD(6) * t329 + qJDD(5) * t355 - t332 * t352;
t320 = mrSges(7,1) * t338 - mrSges(7,3) * t329;
t287 = m(7) * t289 - mrSges(7,2) * t327 + mrSges(7,3) * t316 + t318 * t328 - t320 * t338;
t368 = -t286 * t352 + t355 * t287;
t280 = m(6) * t293 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t333 - qJD(5) * t336 + t330 * t373 + t368;
t292 = -t297 * t353 + t377;
t337 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t373;
t290 = -qJDD(5) * pkin(5) - pkin(11) * t358 - t377 + (qJD(4) * t331 + t297) * t353;
t362 = -m(7) * t290 + t316 * mrSges(7,1) - mrSges(7,2) * t317 + t328 * t319 - t320 * t329;
t284 = m(6) * t292 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t332 + qJD(5) * t337 - t330 * t374 + t362;
t369 = t356 * t280 - t284 * t353;
t277 = m(5) * t301 + t280 * t353 + t284 * t356;
t276 = m(5) * t299 - mrSges(5,1) * t359 - qJDD(4) * mrSges(5,2) + t369;
t281 = t286 * t355 + t287 * t352;
t361 = -m(6) * t296 + t333 * mrSges(6,1) - mrSges(6,2) * t332 - t336 * t374 + t337 * t373 - t281;
t278 = m(5) * t298 + qJDD(4) * mrSges(5,1) - mrSges(5,2) * t359 + t361;
t366 = t276 * t354 + t278 * t357;
t273 = m(4) * t305 - t343 * t277 + t366 * t349;
t275 = m(4) * t306 + t276 * t357 - t278 * t354;
t367 = t273 * t346 + t275 * t340;
t312 = Ifges(7,4) * t329 + Ifges(7,2) * t328 + Ifges(7,6) * t338;
t313 = Ifges(7,1) * t329 + Ifges(7,4) * t328 + Ifges(7,5) * t338;
t360 = mrSges(7,1) * t288 - mrSges(7,2) * t289 + Ifges(7,5) * t317 + Ifges(7,6) * t316 + Ifges(7,3) * t327 + t312 * t329 - t328 * t313;
t324 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t353 + Ifges(6,4) * t356) * qJD(4);
t323 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t353 + Ifges(6,2) * t356) * qJD(4);
t311 = Ifges(7,5) * t329 + Ifges(7,6) * t328 + Ifges(7,3) * t338;
t283 = mrSges(7,2) * t290 - mrSges(7,3) * t288 + Ifges(7,1) * t317 + Ifges(7,4) * t316 + Ifges(7,5) * t327 + t311 * t328 - t312 * t338;
t282 = -mrSges(7,1) * t290 + mrSges(7,3) * t289 + Ifges(7,4) * t317 + Ifges(7,2) * t316 + Ifges(7,6) * t327 - t311 * t329 + t313 * t338;
t274 = m(4) * t308 + t349 * t277 + t366 * t343;
t272 = m(3) * t321 + t350 * t274 + t367 * t344;
t1 = [m(2) * t339 + t351 * t272 + (t341 * (m(3) * t310 - t273 * t340 + t275 * t346) + t347 * (m(3) * t309 - t344 * t274 + t367 * t350)) * t345; t272; t274; Ifges(5,3) * qJDD(4) + mrSges(5,1) * t298 - mrSges(5,2) * t299 + t353 * (mrSges(6,2) * t296 - mrSges(6,3) * t292 + Ifges(6,1) * t332 + Ifges(6,4) * t333 + Ifges(6,5) * qJDD(5) - pkin(11) * t281 - qJD(5) * t323 - t282 * t352 + t283 * t355) + t356 * (-mrSges(6,1) * t296 + mrSges(6,3) * t293 + Ifges(6,4) * t332 + Ifges(6,2) * t333 + Ifges(6,6) * qJDD(5) - pkin(5) * t281 + qJD(5) * t324 - t360) + pkin(4) * t361 + pkin(10) * t369; Ifges(6,5) * t332 + Ifges(6,6) * t333 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t292 - mrSges(6,2) * t293 + t352 * t283 + t355 * t282 + pkin(5) * t362 + pkin(11) * t368 + (t323 * t353 - t324 * t356) * qJD(4); t360;];
tauJ  = t1;
