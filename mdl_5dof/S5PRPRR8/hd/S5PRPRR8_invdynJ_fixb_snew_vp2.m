% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:39
% EndTime: 2019-12-05 16:02:40
% DurationCPUTime: 0.54s
% Computational Cost: add. (2237->160), mult. (4042->204), div. (0->0), fcn. (2494->10), ass. (0->74)
t314 = sin(pkin(9));
t316 = cos(pkin(9));
t302 = g(1) * t314 - g(2) * t316;
t311 = -g(3) + qJDD(1);
t315 = sin(pkin(5));
t317 = cos(pkin(5));
t346 = t302 * t317 + t311 * t315;
t303 = -g(1) * t316 - g(2) * t314;
t320 = sin(qJ(2));
t323 = cos(qJ(2));
t273 = t323 * t303 + t346 * t320;
t345 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t273;
t272 = -t320 * t303 + t346 * t323;
t344 = -pkin(2) - pkin(7);
t285 = -t302 * t315 + t311 * t317;
t319 = sin(qJ(4));
t343 = t285 * t319;
t325 = qJD(2) ^ 2;
t328 = -qJ(3) * t325 + qJDD(3) - t272;
t269 = t344 * qJDD(2) + t328;
t322 = cos(qJ(4));
t265 = t319 * t269 + t322 * t285;
t299 = (t319 * pkin(4) - t322 * pkin(8)) * qJD(2);
t324 = qJD(4) ^ 2;
t340 = qJD(2) * t319;
t262 = -pkin(4) * t324 + qJDD(4) * pkin(8) - t299 * t340 + t265;
t268 = t344 * t325 - t345;
t338 = qJD(2) * qJD(4);
t335 = t322 * t338;
t300 = -qJDD(2) * t319 - t335;
t336 = t319 * t338;
t301 = qJDD(2) * t322 - t336;
t263 = (-t301 + t336) * pkin(8) + (-t300 + t335) * pkin(4) + t268;
t318 = sin(qJ(5));
t321 = cos(qJ(5));
t259 = -t262 * t318 + t263 * t321;
t339 = qJD(2) * t322;
t296 = qJD(4) * t321 - t318 * t339;
t280 = qJD(5) * t296 + qJDD(4) * t318 + t301 * t321;
t297 = qJD(4) * t318 + t321 * t339;
t281 = -mrSges(6,1) * t296 + mrSges(6,2) * t297;
t307 = qJD(5) + t340;
t283 = -mrSges(6,2) * t307 + mrSges(6,3) * t296;
t293 = qJDD(5) - t300;
t257 = m(6) * t259 + mrSges(6,1) * t293 - mrSges(6,3) * t280 - t281 * t297 + t283 * t307;
t260 = t262 * t321 + t263 * t318;
t279 = -qJD(5) * t297 + qJDD(4) * t321 - t301 * t318;
t284 = mrSges(6,1) * t307 - mrSges(6,3) * t297;
t258 = m(6) * t260 - mrSges(6,2) * t293 + mrSges(6,3) * t279 + t281 * t296 - t284 * t307;
t250 = t321 * t257 + t318 * t258;
t334 = -t257 * t318 + t321 * t258;
t298 = (t319 * mrSges(5,1) + t322 * mrSges(5,2)) * qJD(2);
t305 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t339;
t249 = m(5) * t265 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t300 - qJD(4) * t305 - t298 * t340 + t334;
t264 = t269 * t322 - t343;
t304 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t340;
t261 = -qJDD(4) * pkin(4) - pkin(8) * t324 + t343 + (qJD(2) * t299 - t269) * t322;
t329 = -m(6) * t261 + t279 * mrSges(6,1) - mrSges(6,2) * t280 + t296 * t283 - t284 * t297;
t253 = m(5) * t264 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t301 + qJD(4) * t304 - t298 * t339 + t329;
t332 = t249 * t319 + t253 * t322;
t271 = -qJDD(2) * pkin(2) + t328;
t330 = -m(4) * t271 + (t325 * mrSges(4,3)) - t332;
t270 = pkin(2) * t325 + t345;
t327 = -m(4) * t270 + m(5) * t268 - mrSges(5,1) * t300 + (t325 * mrSges(4,2)) + t301 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t304 * t340 + t305 * t339 + t250;
t275 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t307;
t276 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t307;
t326 = mrSges(6,1) * t259 - mrSges(6,2) * t260 + Ifges(6,5) * t280 + Ifges(6,6) * t279 + Ifges(6,3) * t293 + t275 * t297 - t296 * t276;
t289 = (Ifges(5,5) * qJD(4)) + (t322 * Ifges(5,1) - t319 * Ifges(5,4)) * qJD(2);
t288 = (Ifges(5,6) * qJD(4)) + (t322 * Ifges(5,4) - t319 * Ifges(5,2)) * qJD(2);
t274 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t307;
t252 = mrSges(6,2) * t261 - mrSges(6,3) * t259 + Ifges(6,1) * t280 + Ifges(6,4) * t279 + Ifges(6,5) * t293 + t274 * t296 - t275 * t307;
t251 = -mrSges(6,1) * t261 + mrSges(6,3) * t260 + Ifges(6,4) * t280 + Ifges(6,2) * t279 + Ifges(6,6) * t293 - t274 * t297 + t276 * t307;
t248 = qJDD(2) * mrSges(4,2) - t330;
t1 = [m(2) * t311 + t317 * (t249 * t322 - t253 * t319 + (m(3) + m(4)) * t285) + (t320 * (m(3) * t273 - (mrSges(3,1) * t325) - qJDD(2) * mrSges(3,2) + t327) + t323 * (m(3) * t272 - mrSges(3,2) * t325 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t330)) * t315; mrSges(3,1) * t272 - mrSges(3,2) * t273 + mrSges(4,2) * t271 - mrSges(4,3) * t270 + t322 * (mrSges(5,2) * t268 - mrSges(5,3) * t264 + Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * qJDD(4) - pkin(8) * t250 - qJD(4) * t288 - t251 * t318 + t321 * t252) - t319 * (-mrSges(5,1) * t268 + mrSges(5,3) * t265 + Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * qJDD(4) - pkin(4) * t250 + qJD(4) * t289 - t326) - pkin(7) * t332 - pkin(2) * t248 + qJ(3) * t327 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t248; Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t264 - mrSges(5,2) * t265 + t318 * t252 + t321 * t251 + pkin(4) * t329 + pkin(8) * t334 + (t322 * t288 + t319 * t289) * qJD(2); t326;];
tauJ = t1;
