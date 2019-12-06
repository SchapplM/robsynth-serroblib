% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:58
% EndTime: 2019-12-05 16:39:59
% DurationCPUTime: 0.48s
% Computational Cost: add. (1786->134), mult. (2390->162), div. (0->0), fcn. (1321->8), ass. (0->66)
t330 = Ifges(5,1) + Ifges(6,1);
t325 = Ifges(5,4) + Ifges(6,4);
t324 = Ifges(5,5) + Ifges(6,5);
t329 = Ifges(5,2) + Ifges(6,2);
t323 = Ifges(5,6) + Ifges(6,6);
t328 = Ifges(5,3) + Ifges(6,3);
t293 = qJD(2) + qJD(3);
t291 = t293 ^ 2;
t327 = pkin(4) * t291;
t326 = -mrSges(5,2) - mrSges(6,2);
t299 = sin(qJ(4));
t322 = t293 * t299;
t302 = cos(qJ(4));
t321 = t293 * t302;
t297 = sin(pkin(8));
t298 = cos(pkin(8));
t286 = g(1) * t297 - g(2) * t298;
t287 = -g(1) * t298 - g(2) * t297;
t301 = sin(qJ(2));
t304 = cos(qJ(2));
t309 = t304 * t286 - t287 * t301;
t259 = qJDD(2) * pkin(2) + t309;
t317 = t301 * t286 + t304 * t287;
t260 = -qJD(2) ^ 2 * pkin(2) + t317;
t300 = sin(qJ(3));
t303 = cos(qJ(3));
t255 = t300 * t259 + t303 * t260;
t292 = qJDD(2) + qJDD(3);
t252 = -pkin(3) * t291 + pkin(7) * t292 + t255;
t296 = -g(3) + qJDD(1);
t249 = t302 * t252 + t299 * t296;
t320 = (t299 * t324 + t323 * t302) * t293 + t328 * qJD(4);
t319 = (t299 * t325 + t329 * t302) * t293 + t323 * qJD(4);
t318 = (-t330 * t299 - t302 * t325) * t293 - t324 * qJD(4);
t282 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t322;
t316 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t322 - t282;
t315 = qJD(4) * t293;
t314 = qJD(5) * t293;
t311 = t302 * t315;
t275 = t292 * t299 + t311;
t289 = t302 * t296;
t245 = qJDD(4) * pkin(4) + t289 + (-t275 + t311) * qJ(5) + (t302 * t327 - t252 - 0.2e1 * t314) * t299;
t284 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t321;
t313 = m(6) * t245 + qJDD(4) * mrSges(6,1) + qJD(4) * t284;
t276 = t292 * t302 - t299 * t315;
t281 = qJD(4) * pkin(4) - qJ(5) * t322;
t295 = t302 ^ 2;
t246 = qJ(5) * t276 - qJD(4) * t281 - t295 * t327 + 0.2e1 * t302 * t314 + t249;
t273 = (-mrSges(6,1) * t302 + mrSges(6,2) * t299) * t293;
t312 = m(6) * t246 + t276 * mrSges(6,3) + t273 * t321;
t248 = -t299 * t252 + t289;
t274 = (-mrSges(5,1) * t302 + mrSges(5,2) * t299) * t293;
t285 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t321;
t239 = m(5) * t248 + qJDD(4) * mrSges(5,1) + qJD(4) * t285 + (-t273 - t274) * t322 + (-mrSges(5,3) - mrSges(6,3)) * t275 + t313;
t240 = m(5) * t249 + t276 * mrSges(5,3) + t316 * qJD(4) + t326 * qJDD(4) + t274 * t321 + t312;
t310 = -t239 * t299 + t302 * t240;
t254 = t303 * t259 - t300 * t260;
t307 = -t292 * pkin(3) - t254;
t247 = t281 * t322 - t276 * pkin(4) + qJDD(5) + (-qJ(5) * t295 - pkin(7)) * t291 + t307;
t308 = -m(6) * t247 + t276 * mrSges(6,1) + t284 * t321;
t241 = t275 * mrSges(6,2) + t282 * t322 - t308;
t242 = -t275 * mrSges(6,3) - t273 * t322 + t313;
t251 = -pkin(7) * t291 + t307;
t305 = -m(5) * t251 + t276 * mrSges(5,1) + t326 * t275 + t285 * t321 + t316 * t322 + t308;
t306 = -mrSges(4,2) * t255 + t302 * (-mrSges(5,1) * t251 + mrSges(5,3) * t249 - mrSges(6,1) * t247 + mrSges(6,3) * t246 - pkin(4) * t241 + qJ(5) * t312 - t320 * t322 + t329 * t276 + t325 * t275 + (-mrSges(6,2) * qJ(5) + t323) * qJDD(4) + (-qJ(5) * t282 - t318) * qJD(4)) + t299 * (mrSges(5,2) * t251 + mrSges(6,2) * t247 - mrSges(5,3) * t248 - mrSges(6,3) * t245 - qJ(5) * t242 - t319 * qJD(4) + t324 * qJDD(4) + t330 * t275 + t325 * t276 + t320 * t321) + pkin(7) * t310 + pkin(3) * t305 + mrSges(4,1) * t254 + Ifges(4,3) * t292;
t1 = [t302 * t239 + t299 * t240 + (m(2) + m(3) + m(4)) * t296; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t309 - mrSges(3,2) * t317 + pkin(2) * (t300 * (m(4) * t255 - mrSges(4,1) * t291 - mrSges(4,2) * t292 + t310) + t303 * (m(4) * t254 + t292 * mrSges(4,1) - t291 * mrSges(4,2) + t305)) + t306; t306; mrSges(5,1) * t248 + mrSges(6,1) * t245 - mrSges(5,2) * t249 - mrSges(6,2) * t246 + pkin(4) * t242 + t323 * t276 + t324 * t275 + t328 * qJDD(4) + (t319 * t299 + t318 * t302) * t293; t241;];
tauJ = t1;
