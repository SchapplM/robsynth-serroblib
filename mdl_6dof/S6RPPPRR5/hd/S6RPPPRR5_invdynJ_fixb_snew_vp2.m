% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-05-05 13:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:50:15
% EndTime: 2019-05-05 13:50:16
% DurationCPUTime: 0.76s
% Computational Cost: add. (4519->182), mult. (7601->219), div. (0->0), fcn. (3107->8), ass. (0->81)
t369 = -2 * qJD(1);
t340 = sin(qJ(1));
t343 = cos(qJ(1));
t354 = -t343 * g(1) - t340 * g(2);
t368 = -qJDD(1) * qJ(2) + (qJD(2) * t369) - t354;
t367 = -pkin(1) - qJ(3);
t345 = qJD(1) ^ 2;
t366 = qJ(2) * t345;
t335 = -g(3) + qJDD(4);
t342 = cos(qJ(5));
t365 = t335 * t342;
t357 = t340 * g(1) - t343 * g(2);
t353 = qJDD(2) - t357;
t348 = (qJD(3) * t369) + t367 * qJDD(1) + t353;
t305 = (-pkin(3) - qJ(2)) * t345 + t348;
t309 = t367 * t345 + qJDD(3) - t368;
t306 = qJDD(1) * pkin(3) + t309;
t336 = sin(pkin(9));
t337 = cos(pkin(9));
t293 = t337 * t305 + t336 * t306;
t291 = -(pkin(4) * t345) + qJDD(1) * pkin(7) + t293;
t339 = sin(qJ(5));
t288 = t342 * t291 + t339 * t335;
t322 = (-t342 * mrSges(6,1) + t339 * mrSges(6,2)) * qJD(1);
t361 = qJD(1) * qJD(5);
t359 = t339 * t361;
t325 = qJDD(1) * t342 - t359;
t363 = qJD(1) * t339;
t326 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t363;
t292 = -t336 * t305 + t306 * t337;
t290 = -qJDD(1) * pkin(4) - pkin(7) * t345 - t292;
t358 = t342 * t361;
t324 = qJDD(1) * t339 + t358;
t284 = (-t324 - t358) * pkin(8) + (-t325 + t359) * pkin(5) + t290;
t323 = (-t342 * pkin(5) - t339 * pkin(8)) * qJD(1);
t344 = qJD(5) ^ 2;
t362 = qJD(1) * t342;
t286 = -pkin(5) * t344 + qJDD(5) * pkin(8) + t323 * t362 + t288;
t338 = sin(qJ(6));
t341 = cos(qJ(6));
t282 = t284 * t341 - t286 * t338;
t320 = qJD(5) * t341 - t338 * t363;
t300 = qJD(6) * t320 + qJDD(5) * t338 + t324 * t341;
t321 = qJD(5) * t338 + t341 * t363;
t304 = -mrSges(7,1) * t320 + mrSges(7,2) * t321;
t328 = qJD(6) - t362;
t310 = -mrSges(7,2) * t328 + mrSges(7,3) * t320;
t319 = qJDD(6) - t325;
t280 = m(7) * t282 + mrSges(7,1) * t319 - mrSges(7,3) * t300 - t304 * t321 + t310 * t328;
t283 = t284 * t338 + t286 * t341;
t299 = -qJD(6) * t321 + qJDD(5) * t341 - t324 * t338;
t311 = mrSges(7,1) * t328 - mrSges(7,3) * t321;
t281 = m(7) * t283 - mrSges(7,2) * t319 + mrSges(7,3) * t299 + t304 * t320 - t311 * t328;
t355 = -t280 * t338 + t341 * t281;
t274 = m(6) * t288 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t325 - qJD(5) * t326 + t322 * t362 + t355;
t287 = -t291 * t339 + t365;
t327 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t362;
t285 = -qJDD(5) * pkin(5) - pkin(8) * t344 - t365 + (qJD(1) * t323 + t291) * t339;
t349 = -m(7) * t285 + t299 * mrSges(7,1) - mrSges(7,2) * t300 + t320 * t310 - t311 * t321;
t278 = m(6) * t287 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t324 + qJD(5) * t327 - t322 * t363 + t349;
t356 = t342 * t274 - t339 * t278;
t270 = m(5) * t293 - (mrSges(5,1) * t345) - qJDD(1) * mrSges(5,2) + t356;
t275 = t341 * t280 + t338 * t281;
t347 = -m(6) * t290 + t325 * mrSges(6,1) - t324 * mrSges(6,2) - t326 * t363 + t327 * t362 - t275;
t272 = m(5) * t292 + qJDD(1) * mrSges(5,1) - mrSges(5,2) * t345 + t347;
t364 = t336 * t270 + t337 * t272;
t351 = m(4) * t309 + qJDD(1) * mrSges(4,1) - (t345 * mrSges(4,3)) + t364;
t308 = t348 - t366;
t350 = m(4) * t308 + t337 * t270 - t336 * t272;
t295 = Ifges(7,4) * t321 + Ifges(7,2) * t320 + Ifges(7,6) * t328;
t296 = Ifges(7,1) * t321 + Ifges(7,4) * t320 + Ifges(7,5) * t328;
t346 = mrSges(7,1) * t282 - mrSges(7,2) * t283 + Ifges(7,5) * t300 + Ifges(7,6) * t299 + Ifges(7,3) * t319 + t321 * t295 - t320 * t296;
t316 = (Ifges(6,5) * qJD(5)) + (t339 * Ifges(6,1) + t342 * Ifges(6,4)) * qJD(1);
t315 = (Ifges(6,6) * qJD(5)) + (t339 * Ifges(6,4) + t342 * Ifges(6,2)) * qJD(1);
t313 = -qJDD(1) * pkin(1) + t353 - t366;
t312 = pkin(1) * t345 + t368;
t294 = Ifges(7,5) * t321 + Ifges(7,6) * t320 + Ifges(7,3) * t328;
t277 = mrSges(7,2) * t285 - mrSges(7,3) * t282 + Ifges(7,1) * t300 + Ifges(7,4) * t299 + Ifges(7,5) * t319 + t294 * t320 - t295 * t328;
t276 = -mrSges(7,1) * t285 + mrSges(7,3) * t283 + Ifges(7,4) * t300 + Ifges(7,2) * t299 + Ifges(7,6) * t319 - t294 * t321 + t296 * t328;
t268 = m(3) * t313 + (-mrSges(4,1) - mrSges(3,3)) * t345 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t350;
t1 = [-pkin(1) * t268 - mrSges(2,2) * t354 + mrSges(2,1) * t357 + qJ(2) * (-m(3) * t312 + (t345 * mrSges(3,2)) + t351) - qJ(3) * (-t345 * mrSges(4,1) + t350) - mrSges(3,3) * t312 + mrSges(3,2) * t313 + pkin(3) * t364 + mrSges(4,1) * t309 - mrSges(4,3) * t308 + t339 * (mrSges(6,2) * t290 - mrSges(6,3) * t287 + Ifges(6,1) * t324 + Ifges(6,4) * t325 + Ifges(6,5) * qJDD(5) - pkin(8) * t275 - qJD(5) * t315 - t338 * t276 + t341 * t277) + t342 * (-mrSges(6,1) * t290 + mrSges(6,3) * t288 + Ifges(6,4) * t324 + Ifges(6,2) * t325 + Ifges(6,6) * qJDD(5) - pkin(5) * t275 + qJD(5) * t316 - t346) + pkin(4) * t347 + pkin(7) * t356 + mrSges(5,1) * t292 - mrSges(5,2) * t293 + (mrSges(3,3) * qJ(2) + mrSges(4,3) * qJ(3) + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + Ifges(5,3)) * qJDD(1); t268; t351; m(5) * t335 + t274 * t339 + t278 * t342; Ifges(6,5) * t324 + Ifges(6,6) * t325 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t287 - mrSges(6,2) * t288 + t338 * t277 + t341 * t276 + pkin(5) * t349 + pkin(8) * t355 + (t339 * t315 - t342 * t316) * qJD(1); t346;];
tauJ  = t1;
