% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:29
% EndTime: 2020-01-03 11:57:30
% DurationCPUTime: 0.63s
% Computational Cost: add. (4951->131), mult. (6524->179), div. (0->0), fcn. (3631->10), ass. (0->73)
t327 = 2 * qJD(4);
t301 = sin(qJ(1));
t304 = cos(qJ(1));
t310 = -t304 * g(2) - t301 * g(3);
t279 = qJDD(1) * pkin(1) + t310;
t313 = -t301 * g(2) + t304 * g(3);
t280 = -qJD(1) ^ 2 * pkin(1) + t313;
t300 = sin(qJ(2));
t303 = cos(qJ(2));
t266 = t303 * t279 - t300 * t280;
t292 = qJDD(1) + qJDD(2);
t260 = t292 * pkin(2) + t266;
t267 = t300 * t279 + t303 * t280;
t293 = qJD(1) + qJD(2);
t291 = t293 ^ 2;
t261 = -t291 * pkin(2) + t267;
t296 = sin(pkin(8));
t298 = cos(pkin(8));
t256 = t296 * t260 + t298 * t261;
t253 = -t291 * pkin(3) + t292 * qJ(4) + t256;
t326 = t293 * t327 + t253;
t295 = sin(pkin(9));
t325 = mrSges(5,2) * t295;
t323 = mrSges(5,3) * t292;
t322 = t293 * t295;
t297 = cos(pkin(9));
t321 = t297 * t292;
t320 = t297 * t293;
t294 = -g(1) + qJDD(3);
t319 = t297 * t294;
t249 = t295 * t294 + t326 * t297;
t311 = -pkin(4) * t297 - pkin(7) * t295;
t275 = t311 * t293;
t247 = t275 * t320 + t249;
t255 = t298 * t260 - t296 * t261;
t308 = -t291 * qJ(4) + qJDD(4) - t255;
t250 = (-pkin(3) + t311) * t292 + t308;
t299 = sin(qJ(5));
t302 = cos(qJ(5));
t244 = -t299 * t247 + t302 * t250;
t282 = qJD(5) - t320;
t315 = t299 * t322;
t268 = -t282 * mrSges(6,2) - mrSges(6,3) * t315;
t270 = (t299 * mrSges(6,1) + t302 * mrSges(6,2)) * t322;
t316 = qJD(5) * t293;
t272 = (t292 * t302 - t299 * t316) * t295;
t281 = qJDD(5) - t321;
t314 = t302 * t322;
t242 = m(6) * t244 + t281 * mrSges(6,1) - t272 * mrSges(6,3) + t282 * t268 - t270 * t314;
t245 = t302 * t247 + t299 * t250;
t269 = t282 * mrSges(6,1) - mrSges(6,3) * t314;
t271 = (-t292 * t299 - t302 * t316) * t295;
t243 = m(6) * t245 - t281 * mrSges(6,2) + t271 * mrSges(6,3) - t282 * t269 - t270 * t315;
t273 = (-mrSges(5,1) * t297 + t325) * t293;
t237 = m(5) * t249 - t299 * t242 + t302 * t243 + (t273 * t293 + t323) * t297;
t246 = -t319 + (t253 + (t327 + t275) * t293) * t295;
t248 = -t326 * t295 + t319;
t241 = m(5) * t248 - m(6) * t246 + t271 * mrSges(6,1) - t272 * mrSges(6,2) + (-t323 + (-t268 * t299 - t269 * t302 - t273) * t293) * t295;
t312 = t297 * t237 - t295 * t241;
t232 = m(4) * t256 - t291 * mrSges(4,1) - t292 * mrSges(4,2) + t312;
t240 = t302 * t242 + t299 * t243;
t252 = -t292 * pkin(3) + t308;
t307 = -m(5) * t252 + mrSges(5,1) * t321 - t240 + (t295 ^ 2 + t297 ^ 2) * mrSges(5,3) * t291;
t235 = m(4) * t255 - t291 * mrSges(4,2) + (mrSges(4,1) - t325) * t292 + t307;
t318 = t296 * t232 + t298 * t235;
t263 = Ifges(6,6) * t282 + (t302 * Ifges(6,4) - t299 * Ifges(6,2)) * t322;
t264 = Ifges(6,5) * t282 + (t302 * Ifges(6,1) - t299 * Ifges(6,4)) * t322;
t309 = t302 * t263 + t299 * t264;
t306 = mrSges(6,1) * t244 - mrSges(6,2) * t245 + Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t281;
t239 = t292 * t325 - t307;
t274 = (Ifges(5,5) * t295 + Ifges(5,6) * t297) * t293;
t305 = -mrSges(3,2) * t267 - mrSges(4,2) * t256 + pkin(2) * t318 + t295 * (t274 * t320 + mrSges(5,2) * t252 - mrSges(5,3) * t248 + t302 * (mrSges(6,2) * t246 - mrSges(6,3) * t244 + Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t281 - t282 * t263) - t299 * (-mrSges(6,1) * t246 + mrSges(6,3) * t245 + Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t281 + t282 * t264) - pkin(7) * t240 + (Ifges(5,1) * t295 + Ifges(5,4) * t297) * t292) + t297 * (Ifges(5,2) * t321 - mrSges(5,1) * t252 + mrSges(5,3) * t249 - pkin(4) * t240 + (Ifges(5,4) * t292 + (-t274 - t309) * t293) * t295 - t306) + qJ(4) * t312 - pkin(3) * t239 + mrSges(4,1) * t255 + mrSges(3,1) * t266 + Ifges(4,3) * t292 + Ifges(3,3) * t292;
t1 = [-mrSges(2,2) * t313 + pkin(1) * (t300 * (m(3) * t267 - t291 * mrSges(3,1) - t292 * mrSges(3,2) + t298 * t232 - t296 * t235) + t303 * (m(3) * t266 + t292 * mrSges(3,1) - t291 * mrSges(3,2) + t318)) + mrSges(2,1) * t310 + t305 + Ifges(2,3) * qJDD(1); t305; m(4) * t294 + t295 * t237 + t297 * t241; t239; t309 * t322 + t306;];
tauJ = t1;
