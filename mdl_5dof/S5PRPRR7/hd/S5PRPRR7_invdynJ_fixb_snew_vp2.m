% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:51
% EndTime: 2019-12-05 15:59:51
% DurationCPUTime: 0.48s
% Computational Cost: add. (1956->155), mult. (3597->200), div. (0->0), fcn. (2087->8), ass. (0->67)
t297 = sin(pkin(8));
t298 = cos(pkin(8));
t282 = -g(1) * t298 - g(2) * t297;
t294 = -g(3) + qJDD(1);
t301 = sin(qJ(2));
t304 = cos(qJ(2));
t266 = t304 * t282 + t301 * t294;
t312 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t266;
t320 = -pkin(2) - pkin(6);
t305 = qJD(2) ^ 2;
t319 = pkin(4) * t305;
t265 = -t301 * t282 + t294 * t304;
t310 = -qJ(3) * t305 + qJDD(3) - t265;
t260 = qJDD(2) * t320 + t310;
t303 = cos(qJ(4));
t257 = t303 * t260;
t300 = sin(qJ(4));
t316 = qJD(2) * qJD(4);
t280 = qJDD(2) * t303 - t300 * t316;
t281 = -g(1) * t297 + g(2) * t298;
t239 = (qJDD(4) * pkin(4)) - pkin(7) * t280 + t257 + (-pkin(7) * t316 - t303 * t319 - t281) * t300;
t248 = t300 * t260 + t303 * t281;
t279 = -qJDD(2) * t300 - t303 * t316;
t317 = t303 * qJD(2);
t285 = (qJD(4) * pkin(4)) - pkin(7) * t317;
t293 = t300 ^ 2;
t240 = pkin(7) * t279 - qJD(4) * t285 - t293 * t319 + t248;
t299 = sin(qJ(5));
t302 = cos(qJ(5));
t237 = t239 * t302 - t240 * t299;
t270 = (-t299 * t303 - t300 * t302) * qJD(2);
t250 = qJD(5) * t270 + t279 * t299 + t280 * t302;
t271 = (-t299 * t300 + t302 * t303) * qJD(2);
t255 = -mrSges(6,1) * t270 + mrSges(6,2) * t271;
t290 = qJD(4) + qJD(5);
t263 = -mrSges(6,2) * t290 + mrSges(6,3) * t270;
t289 = qJDD(4) + qJDD(5);
t234 = m(6) * t237 + mrSges(6,1) * t289 - mrSges(6,3) * t250 - t255 * t271 + t263 * t290;
t238 = t239 * t299 + t240 * t302;
t249 = -qJD(5) * t271 + t279 * t302 - t280 * t299;
t264 = mrSges(6,1) * t290 - mrSges(6,3) * t271;
t235 = m(6) * t238 - mrSges(6,2) * t289 + mrSges(6,3) * t249 + t255 * t270 - t264 * t290;
t228 = t302 * t234 + t299 * t235;
t318 = qJD(2) * t300;
t313 = -t234 * t299 + t302 * t235;
t247 = -t281 * t300 + t257;
t278 = (mrSges(5,1) * t300 + mrSges(5,2) * t303) * qJD(2);
t283 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t318;
t284 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t317;
t311 = (m(5) * t247 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t280 + qJD(4) * t283 - t278 * t317 + t228) * t303 + (m(5) * t248 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t279 - qJD(4) * t284 - t278 * t318 + t313) * t300;
t242 = t285 * t317 - pkin(4) * t279 + (-pkin(7) * t293 + t320) * t305 + t312;
t309 = m(6) * t242 - mrSges(6,1) * t249 + t250 * mrSges(6,2) - t263 * t270 + t271 * t264;
t262 = -qJDD(2) * pkin(2) + t310;
t308 = -m(4) * t262 + (t305 * mrSges(4,3)) - t311;
t252 = Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t290;
t253 = Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t290;
t307 = mrSges(6,1) * t237 - mrSges(6,2) * t238 + Ifges(6,5) * t250 + Ifges(6,6) * t249 + (Ifges(6,3) * t289) + t271 * t252 - t253 * t270;
t259 = t305 * t320 + t312;
t261 = pkin(2) * t305 - t312;
t306 = -m(4) * t261 + m(5) * t259 - mrSges(5,1) * t279 + (t305 * mrSges(4,2)) + t280 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t283 * t318 + t284 * t317 + t309;
t269 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t303 - Ifges(5,4) * t300) * qJD(2);
t268 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t303 - Ifges(5,2) * t300) * qJD(2);
t251 = Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t290;
t230 = mrSges(6,2) * t242 - mrSges(6,3) * t237 + Ifges(6,1) * t250 + Ifges(6,4) * t249 + Ifges(6,5) * t289 + t251 * t270 - t252 * t290;
t229 = -mrSges(6,1) * t242 + mrSges(6,3) * t238 + Ifges(6,4) * t250 + Ifges(6,2) * t249 + Ifges(6,6) * t289 - t251 * t271 + t253 * t290;
t225 = qJDD(2) * mrSges(4,2) - t308;
t1 = [m(2) * t294 + t301 * (m(3) * t266 - (mrSges(3,1) * t305) - qJDD(2) * mrSges(3,2) + t306) + t304 * (m(3) * t265 - mrSges(3,2) * t305 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t308); mrSges(3,1) * t265 - mrSges(3,2) * t266 + mrSges(4,2) * t262 - mrSges(4,3) * t261 + t303 * (mrSges(5,2) * t259 - mrSges(5,3) * t247 + Ifges(5,1) * t280 + Ifges(5,4) * t279 + (Ifges(5,5) * qJDD(4)) - pkin(7) * t228 - qJD(4) * t268 - t229 * t299 + t230 * t302) - t300 * (-mrSges(5,1) * t259 + mrSges(5,3) * t248 + Ifges(5,4) * t280 + Ifges(5,2) * t279 + (Ifges(5,6) * qJDD(4)) - pkin(4) * t309 + pkin(7) * t313 + qJD(4) * t269 + t302 * t229 + t299 * t230) - pkin(6) * t311 - pkin(2) * t225 + qJ(3) * t306 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t225; mrSges(5,1) * t247 - mrSges(5,2) * t248 + Ifges(5,5) * t280 + Ifges(5,6) * t279 + (Ifges(5,3) * qJDD(4)) + pkin(4) * t228 + (t268 * t303 + t269 * t300) * qJD(2) + t307; t307;];
tauJ = t1;
