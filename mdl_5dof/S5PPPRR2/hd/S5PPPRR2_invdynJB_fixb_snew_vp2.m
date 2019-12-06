% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPPRR2
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:33
% EndTime: 2019-12-05 14:59:34
% DurationCPUTime: 1.14s
% Computational Cost: add. (11675->151), mult. (17698->200), div. (0->0), fcn. (12690->10), ass. (0->74)
t347 = sin(pkin(7));
t350 = cos(pkin(7));
t340 = -g(1) * t350 - g(2) * t347;
t344 = -g(3) + qJDD(1);
t346 = sin(pkin(8));
t349 = cos(pkin(8));
t328 = t340 * t349 + t344 * t346;
t339 = g(1) * t347 - g(2) * t350;
t338 = qJDD(2) - t339;
t345 = sin(pkin(9));
t348 = cos(pkin(9));
t324 = t328 * t348 + t338 * t345;
t327 = -t340 * t346 + t344 * t349;
t326 = qJDD(3) - t327;
t352 = sin(qJ(4));
t354 = cos(qJ(4));
t321 = t354 * t324 + t352 * t326;
t355 = qJD(4) ^ 2;
t319 = -pkin(4) * t355 + qJDD(4) * pkin(6) + t321;
t323 = t328 * t345 - t348 * t338;
t351 = sin(qJ(5));
t353 = cos(qJ(5));
t316 = -t319 * t351 + t323 * t353;
t317 = t319 * t353 + t323 * t351;
t330 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t351 + Ifges(6,2) * t353) * qJD(4);
t331 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t351 + Ifges(6,4) * t353) * qJD(4);
t364 = qJD(4) * qJD(5);
t336 = qJDD(4) * t351 + t353 * t364;
t337 = qJDD(4) * t353 - t351 * t364;
t368 = mrSges(6,1) * t316 - mrSges(6,2) * t317 + Ifges(6,5) * t336 + Ifges(6,6) * t337 + Ifges(6,3) * qJDD(5) + (t330 * t351 - t331 * t353) * qJD(4);
t335 = (-mrSges(6,1) * t353 + mrSges(6,2) * t351) * qJD(4);
t365 = qJD(4) * t353;
t342 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t365;
t366 = qJD(4) * t351;
t313 = m(6) * t316 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t336 + qJD(5) * t342 - t335 * t366;
t341 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t366;
t314 = m(6) * t317 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t337 - qJD(5) * t341 + t335 * t365;
t308 = -t313 * t351 + t353 * t314;
t306 = m(5) * t321 - mrSges(5,1) * t355 - qJDD(4) * mrSges(5,2) + t308;
t320 = -t324 * t352 + t326 * t354;
t318 = -qJDD(4) * pkin(4) - pkin(6) * t355 - t320;
t315 = -m(6) * t318 + t337 * mrSges(6,1) - mrSges(6,2) * t336 - t341 * t366 + t342 * t365;
t311 = m(5) * t320 + qJDD(4) * mrSges(5,1) - mrSges(5,2) * t355 + t315;
t359 = t354 * t306 - t311 * t352;
t301 = m(4) * t324 + t359;
t307 = t353 * t313 + t351 * t314;
t304 = (-m(4) - m(5)) * t323 - t307;
t360 = t348 * t301 - t304 * t345;
t293 = m(3) * t328 + t360;
t303 = t306 * t352 + t311 * t354;
t302 = m(4) * t326 + t303;
t299 = m(3) * t327 - t302;
t361 = t349 * t293 - t299 * t346;
t287 = m(2) * t340 + t361;
t295 = t301 * t345 + t304 * t348;
t294 = m(3) * t338 + t295;
t292 = m(2) * t339 - t294;
t367 = t347 * t287 + t350 * t292;
t288 = t346 * t293 + t349 * t299;
t363 = m(2) * t344 + t288;
t362 = t350 * t287 - t292 * t347;
t329 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t351 + Ifges(6,6) * t353) * qJD(4);
t309 = -mrSges(6,1) * t318 + mrSges(6,3) * t317 + Ifges(6,4) * t336 + Ifges(6,2) * t337 + Ifges(6,6) * qJDD(5) + qJD(5) * t331 - t329 * t366;
t310 = mrSges(6,2) * t318 - mrSges(6,3) * t316 + Ifges(6,1) * t336 + Ifges(6,4) * t337 + Ifges(6,5) * qJDD(5) - qJD(5) * t330 + t329 * t365;
t356 = mrSges(5,1) * t320 - mrSges(5,2) * t321 + Ifges(5,3) * qJDD(4) + pkin(4) * t315 + pkin(6) * t308 + t309 * t353 + t310 * t351;
t297 = -mrSges(5,1) * t323 + mrSges(5,3) * t321 + t355 * Ifges(5,5) + Ifges(5,6) * qJDD(4) - pkin(4) * t307 - t368;
t296 = mrSges(5,2) * t323 - mrSges(5,3) * t320 + Ifges(5,5) * qJDD(4) - Ifges(5,6) * t355 - pkin(6) * t307 - t309 * t351 + t310 * t353;
t284 = -mrSges(4,1) * t326 + mrSges(4,3) * t324 - pkin(3) * t303 - t356;
t283 = mrSges(4,2) * t326 + mrSges(4,3) * t323 - pkin(5) * t303 + t296 * t354 - t297 * t352;
t282 = -mrSges(3,1) * t338 + mrSges(3,3) * t328 + mrSges(4,1) * t323 + mrSges(4,2) * t324 - t352 * t296 - t354 * t297 - pkin(3) * (-m(5) * t323 - t307) - pkin(5) * t359 - pkin(2) * t295;
t281 = mrSges(3,2) * t338 - mrSges(3,3) * t327 - qJ(3) * t295 + t283 * t348 - t284 * t345;
t280 = -mrSges(2,1) * t344 - mrSges(3,1) * t327 + mrSges(3,2) * t328 + mrSges(2,3) * t340 - pkin(1) * t288 + pkin(2) * t302 - qJ(3) * t360 - t345 * t283 - t348 * t284;
t279 = mrSges(2,2) * t344 - mrSges(2,3) * t339 - qJ(2) * t288 + t281 * t349 - t282 * t346;
t1 = [-m(1) * g(1) + t362; -m(1) * g(2) + t367; -m(1) * g(3) + t363; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t367 + t350 * t279 - t347 * t280; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t362 + t347 * t279 + t350 * t280; -mrSges(1,1) * g(2) + mrSges(2,1) * t339 + mrSges(1,2) * g(1) - mrSges(2,2) * t340 - pkin(1) * t294 + qJ(2) * t361 + t346 * t281 + t349 * t282; t363; t294; t302; t356; t368;];
tauJB = t1;
