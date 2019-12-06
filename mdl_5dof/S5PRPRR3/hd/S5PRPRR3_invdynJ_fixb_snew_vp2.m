% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR3
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:50
% EndTime: 2019-12-05 15:46:51
% DurationCPUTime: 0.49s
% Computational Cost: add. (2790->157), mult. (5014->208), div. (0->0), fcn. (3089->10), ass. (0->71)
t296 = sin(pkin(8));
t298 = cos(pkin(8));
t284 = -g(1) * t298 - g(2) * t296;
t294 = -g(3) + qJDD(1);
t301 = sin(qJ(2));
t304 = cos(qJ(2));
t268 = -t284 * t301 + t304 * t294;
t264 = qJDD(2) * pkin(2) + t268;
t269 = t304 * t284 + t301 * t294;
t305 = qJD(2) ^ 2;
t265 = -pkin(2) * t305 + t269;
t295 = sin(pkin(9));
t297 = cos(pkin(9));
t250 = t295 * t264 + t297 * t265;
t247 = -pkin(3) * t305 + qJDD(2) * pkin(6) + t250;
t283 = -g(1) * t296 + g(2) * t298 + qJDD(3);
t300 = sin(qJ(4));
t303 = cos(qJ(4));
t243 = -t247 * t300 + t303 * t283;
t313 = qJD(2) * qJD(4);
t312 = t303 * t313;
t281 = qJDD(2) * t300 + t312;
t240 = (-t281 + t312) * pkin(7) + (t300 * t303 * t305 + qJDD(4)) * pkin(4) + t243;
t244 = t303 * t247 + t300 * t283;
t282 = qJDD(2) * t303 - t300 * t313;
t315 = qJD(2) * t300;
t287 = qJD(4) * pkin(4) - pkin(7) * t315;
t293 = t303 ^ 2;
t241 = -pkin(4) * t293 * t305 + pkin(7) * t282 - qJD(4) * t287 + t244;
t299 = sin(qJ(5));
t302 = cos(qJ(5));
t238 = t240 * t302 - t241 * t299;
t273 = (-t300 * t299 + t303 * t302) * qJD(2);
t255 = qJD(5) * t273 + t281 * t302 + t282 * t299;
t274 = (t303 * t299 + t300 * t302) * qJD(2);
t260 = -mrSges(6,1) * t273 + mrSges(6,2) * t274;
t292 = qJD(4) + qJD(5);
t266 = -mrSges(6,2) * t292 + mrSges(6,3) * t273;
t291 = qJDD(4) + qJDD(5);
t235 = m(6) * t238 + mrSges(6,1) * t291 - mrSges(6,3) * t255 - t260 * t274 + t266 * t292;
t239 = t240 * t299 + t241 * t302;
t254 = -qJD(5) * t274 - t281 * t299 + t282 * t302;
t267 = mrSges(6,1) * t292 - mrSges(6,3) * t274;
t236 = m(6) * t239 - mrSges(6,2) * t291 + mrSges(6,3) * t254 + t260 * t273 - t267 * t292;
t227 = t302 * t235 + t299 * t236;
t280 = (-t303 * mrSges(5,1) + t300 * mrSges(5,2)) * qJD(2);
t314 = qJD(2) * t303;
t286 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t314;
t225 = m(5) * t243 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t281 + qJD(4) * t286 - t280 * t315 + t227;
t285 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t315;
t310 = -t235 * t299 + t302 * t236;
t226 = m(5) * t244 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t282 - qJD(4) * t285 + t280 * t314 + t310;
t311 = -t225 * t300 + t303 * t226;
t223 = m(4) * t250 - mrSges(4,1) * t305 - qJDD(2) * mrSges(4,2) + t311;
t249 = t264 * t297 - t295 * t265;
t309 = -qJDD(2) * pkin(3) - t249;
t246 = -pkin(6) * t305 + t309;
t242 = t287 * t315 - pkin(4) * t282 + (-pkin(7) * t293 - pkin(6)) * t305 + t309;
t308 = m(6) * t242 - t254 * mrSges(6,1) + mrSges(6,2) * t255 - t273 * t266 + t267 * t274;
t306 = -m(5) * t246 + t282 * mrSges(5,1) - mrSges(5,2) * t281 - t285 * t315 + t286 * t314 - t308;
t231 = m(4) * t249 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t305 + t306;
t316 = t295 * t223 + t297 * t231;
t257 = Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * t292;
t258 = Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * t292;
t307 = mrSges(6,1) * t238 - mrSges(6,2) * t239 + Ifges(6,5) * t255 + Ifges(6,6) * t254 + Ifges(6,3) * t291 + t274 * t257 - t258 * t273;
t272 = Ifges(5,5) * qJD(4) + (t300 * Ifges(5,1) + t303 * Ifges(5,4)) * qJD(2);
t271 = Ifges(5,6) * qJD(4) + (t300 * Ifges(5,4) + t303 * Ifges(5,2)) * qJD(2);
t256 = Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * t292;
t229 = mrSges(6,2) * t242 - mrSges(6,3) * t238 + Ifges(6,1) * t255 + Ifges(6,4) * t254 + Ifges(6,5) * t291 + t256 * t273 - t257 * t292;
t228 = -mrSges(6,1) * t242 + mrSges(6,3) * t239 + Ifges(6,4) * t255 + Ifges(6,2) * t254 + Ifges(6,6) * t291 - t256 * t274 + t258 * t292;
t1 = [m(2) * t294 + t301 * (m(3) * t269 - mrSges(3,1) * t305 - qJDD(2) * mrSges(3,2) + t223 * t297 - t231 * t295) + t304 * (m(3) * t268 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t305 + t316); mrSges(3,1) * t268 - mrSges(3,2) * t269 + mrSges(4,1) * t249 - mrSges(4,2) * t250 + t300 * (mrSges(5,2) * t246 - mrSges(5,3) * t243 + Ifges(5,1) * t281 + Ifges(5,4) * t282 + Ifges(5,5) * qJDD(4) - pkin(7) * t227 - qJD(4) * t271 - t228 * t299 + t229 * t302) + t303 * (-mrSges(5,1) * t246 + mrSges(5,3) * t244 + Ifges(5,4) * t281 + Ifges(5,2) * t282 + Ifges(5,6) * qJDD(4) - pkin(4) * t308 + pkin(7) * t310 + qJD(4) * t272 + t302 * t228 + t299 * t229) + pkin(3) * t306 + pkin(6) * t311 + pkin(2) * t316 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); m(4) * t283 + t225 * t303 + t226 * t300; mrSges(5,1) * t243 - mrSges(5,2) * t244 + Ifges(5,5) * t281 + Ifges(5,6) * t282 + Ifges(5,3) * qJDD(4) + pkin(4) * t227 + (t300 * t271 - t303 * t272) * qJD(2) + t307; t307;];
tauJ = t1;
