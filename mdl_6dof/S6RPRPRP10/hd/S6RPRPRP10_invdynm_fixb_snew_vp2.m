% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:09:35
% EndTime: 2019-05-05 18:09:45
% DurationCPUTime: 3.36s
% Computational Cost: add. (35926->343), mult. (70241->390), div. (0->0), fcn. (34991->6), ass. (0->127)
t360 = -2 * qJD(4);
t295 = sin(qJ(1));
t297 = cos(qJ(1));
t270 = t295 * g(1) - t297 * g(2);
t299 = qJD(1) ^ 2;
t320 = -t299 * qJ(2) + qJDD(2) - t270;
t354 = -pkin(1) - pkin(7);
t226 = t354 * qJDD(1) + t320;
t296 = cos(qJ(3));
t341 = t296 * t226;
t294 = sin(qJ(3));
t349 = t294 * g(3);
t217 = t341 + t349;
t336 = qJD(1) * qJD(3);
t278 = t294 * t336;
t261 = qJDD(1) * t296 - t278;
t338 = qJD(1) * t294;
t265 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t338;
t267 = mrSges(5,1) * t338 - qJD(3) * mrSges(5,3);
t331 = t296 * t336;
t260 = qJDD(1) * t294 + t331;
t337 = qJD(1) * t296;
t269 = pkin(4) * t337 - qJD(3) * pkin(8);
t290 = t294 ^ 2;
t271 = -t297 * g(1) - t295 * g(2);
t322 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t271;
t311 = pkin(3) * t331 + t337 * t360 + t322 + (-t261 + t278) * qJ(4);
t353 = pkin(3) + pkin(8);
t180 = -t269 * t337 + t353 * t260 + (-pkin(4) * t290 + t354) * t299 + t311;
t257 = (pkin(3) * t294 - qJ(4) * t296) * qJD(1);
t298 = qJD(3) ^ 2;
t319 = -t298 * qJ(4) + t257 * t337 + qJDD(4) - t341;
t350 = pkin(8) * t299;
t184 = t261 * pkin(4) - t353 * qJDD(3) + (pkin(4) * t336 + t296 * t350 - g(3)) * t294 + t319;
t293 = sin(qJ(5));
t352 = cos(qJ(5));
t178 = t352 * t180 + t293 * t184;
t256 = t352 * qJD(3) + t293 * t338;
t209 = qJD(5) * t256 + qJDD(3) * t293 - t352 * t260;
t276 = qJD(5) + t337;
t220 = mrSges(6,1) * t276 - mrSges(6,3) * t256;
t254 = qJDD(5) + t261;
t255 = qJD(3) * t293 - t352 * t338;
t214 = pkin(5) * t255 - qJ(6) * t256;
t272 = t276 ^ 2;
t171 = -pkin(5) * t272 + qJ(6) * t254 + 0.2e1 * qJD(6) * t276 - t214 * t255 + t178;
t221 = -mrSges(7,1) * t276 + mrSges(7,2) * t256;
t332 = m(7) * t171 + t254 * mrSges(7,3) + t276 * t221;
t215 = mrSges(7,1) * t255 - mrSges(7,3) * t256;
t339 = -mrSges(6,1) * t255 - mrSges(6,2) * t256 - t215;
t346 = -mrSges(6,3) - mrSges(7,2);
t159 = m(6) * t178 - t254 * mrSges(6,2) + t346 * t209 - t276 * t220 + t339 * t255 + t332;
t177 = -t293 * t180 + t352 * t184;
t210 = -t255 * qJD(5) + t352 * qJDD(3) + t293 * t260;
t219 = -mrSges(6,2) * t276 - mrSges(6,3) * t255;
t173 = -t254 * pkin(5) - t272 * qJ(6) + t256 * t214 + qJDD(6) - t177;
t222 = -mrSges(7,2) * t255 + mrSges(7,3) * t276;
t328 = -m(7) * t173 + t254 * mrSges(7,1) + t276 * t222;
t160 = m(6) * t177 + t254 * mrSges(6,1) + t346 * t210 + t276 * t219 + t339 * t256 + t328;
t154 = t293 * t159 + t352 * t160;
t189 = -qJDD(3) * pkin(3) + t319 - t349;
t316 = -m(5) * t189 - t261 * mrSges(5,1) - t154;
t258 = (-mrSges(5,2) * t294 - mrSges(5,3) * t296) * qJD(1);
t329 = qJD(1) * (-t258 - (mrSges(4,1) * t294 + mrSges(4,2) * t296) * qJD(1));
t347 = mrSges(4,1) - mrSges(5,2);
t147 = m(4) * t217 - t261 * mrSges(4,3) + t347 * qJDD(3) + (t265 - t267) * qJD(3) + t296 * t329 + t316;
t218 = -t296 * g(3) + t294 * t226;
t266 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t337;
t187 = t298 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t360 + t257 * t338 - t218;
t268 = mrSges(5,1) * t337 + qJD(3) * mrSges(5,2);
t183 = -t260 * pkin(4) + qJD(3) * t269 - t290 * t350 - t187;
t175 = -0.2e1 * qJD(6) * t256 + (t255 * t276 - t210) * qJ(6) + (t256 * t276 + t209) * pkin(5) + t183;
t168 = m(7) * t175 + t209 * mrSges(7,1) - mrSges(7,3) * t210 - t221 * t256 + t255 * t222;
t307 = m(6) * t183 + mrSges(6,1) * t209 + t210 * mrSges(6,2) + t219 * t255 + t256 * t220 + t168;
t305 = -m(5) * t187 + qJDD(3) * mrSges(5,3) + qJD(3) * t268 + t307;
t157 = -qJDD(3) * mrSges(4,2) + t305 + (-mrSges(4,3) - mrSges(5,1)) * t260 + m(4) * t218 - qJD(3) * t266 + t294 * t329;
t142 = t296 * t147 + t294 * t157;
t231 = -qJDD(1) * pkin(1) + t320;
t234 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t296 - Ifges(4,2) * t294) * qJD(1);
t235 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t296 - Ifges(4,4) * t294) * qJD(1);
t197 = Ifges(7,1) * t256 + Ifges(7,4) * t276 + Ifges(7,5) * t255;
t198 = Ifges(6,1) * t256 - Ifges(6,4) * t255 + Ifges(6,5) * t276;
t327 = -mrSges(7,1) * t175 + mrSges(7,2) * t171;
t195 = Ifges(7,4) * t256 + Ifges(7,2) * t276 + Ifges(7,6) * t255;
t340 = -Ifges(6,5) * t256 + Ifges(6,6) * t255 - Ifges(6,3) * t276 - t195;
t151 = -mrSges(6,1) * t183 + mrSges(6,3) * t178 - pkin(5) * t168 + (t197 + t198) * t276 + t340 * t256 + (Ifges(6,6) - Ifges(7,6)) * t254 + (Ifges(6,4) - Ifges(7,5)) * t210 + (-Ifges(6,2) - Ifges(7,3)) * t209 + t327;
t196 = Ifges(6,4) * t256 - Ifges(6,2) * t255 + Ifges(6,6) * t276;
t193 = Ifges(7,5) * t256 + Ifges(7,6) * t276 + Ifges(7,3) * t255;
t318 = mrSges(7,2) * t173 - mrSges(7,3) * t175 + Ifges(7,1) * t210 + Ifges(7,4) * t254 + Ifges(7,5) * t209 + t276 * t193;
t153 = mrSges(6,2) * t183 - mrSges(6,3) * t177 + Ifges(6,1) * t210 - Ifges(6,4) * t209 + Ifges(6,5) * t254 - qJ(6) * t168 - t276 * t196 + t340 * t255 + t318;
t236 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t296 + Ifges(5,3) * t294) * qJD(1);
t237 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t296 + Ifges(5,6) * t294) * qJD(1);
t357 = mrSges(5,2) * t189 - mrSges(5,3) * t187 + Ifges(5,1) * qJDD(3) - Ifges(5,4) * t261 + Ifges(5,5) * t260 - pkin(8) * t154 - t293 * t151 + t352 * t153 - (t236 * t296 + t237 * t294) * qJD(1);
t300 = -mrSges(4,2) * t218 + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t267 - t258 * t337 + t316) + qJ(4) * (-mrSges(5,1) * t260 - t258 * t338 + t305) + mrSges(4,1) * t217 + t235 * t338 + t234 * t337 - Ifges(4,6) * t260 + Ifges(4,5) * t261 + Ifges(4,3) * qJDD(3) + t357;
t359 = mrSges(3,1) * t231 + pkin(2) * t142 + t300;
t333 = t354 * t299;
t225 = t333 + t322;
t155 = t352 * t159 - t293 * t160;
t186 = t260 * pkin(3) + t311 + t333;
t309 = m(5) * t186 - t261 * mrSges(5,3) - (t267 * t294 + t268 * t296) * qJD(1) + t155;
t358 = -m(4) * t225 - t261 * mrSges(4,2) - t347 * t260 - t265 * t338 - t266 * t337 - t309;
t348 = mrSges(2,1) - mrSges(3,2);
t345 = Ifges(2,5) - Ifges(3,4);
t344 = -Ifges(2,6) + Ifges(3,5);
t343 = -Ifges(5,6) - Ifges(4,4);
t143 = -t147 * t294 + t296 * t157;
t238 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t296 + Ifges(5,5) * t294) * qJD(1);
t330 = qJD(1) * (-Ifges(4,3) * qJD(3) - (Ifges(4,5) * t296 - Ifges(4,6) * t294) * qJD(1) - t238);
t317 = -m(3) * t231 + t299 * mrSges(3,3) - t142;
t315 = mrSges(7,1) * t173 - mrSges(7,3) * t171 - Ifges(7,4) * t210 - Ifges(7,2) * t254 - Ifges(7,6) * t209 + t256 * t193 - t255 * t197;
t148 = -t260 * mrSges(5,2) + t309;
t308 = mrSges(5,1) * t187 - mrSges(5,2) * t186 - pkin(4) * t307 + pkin(8) * t155 + t352 * t151 + t293 * t153;
t136 = -mrSges(4,1) * t225 + mrSges(4,3) * t218 - pkin(3) * t148 - t343 * t261 + (-Ifges(4,2) - Ifges(5,3)) * t260 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t235 - t237) * qJD(3) + t296 * t330 - t308;
t304 = mrSges(6,2) * t178 - t255 * t198 - qJ(6) * (-t209 * mrSges(7,2) - t255 * t215 + t332) - pkin(5) * (-t210 * mrSges(7,2) - t256 * t215 + t328) - mrSges(6,1) * t177 - t256 * t196 + Ifges(6,6) * t209 - Ifges(6,5) * t210 - Ifges(6,3) * t254 + t315;
t301 = -mrSges(5,1) * t189 + mrSges(5,3) * t186 - pkin(4) * t154 + t304;
t138 = (Ifges(5,2) + Ifges(4,1)) * t261 + t343 * t260 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) + (-t234 + t236) * qJD(3) + mrSges(4,2) * t225 - mrSges(4,3) * t217 - qJ(4) * t148 - t301 + t294 * t330;
t229 = t299 * pkin(1) - t322;
t313 = mrSges(3,2) * t231 - mrSges(3,3) * t229 + Ifges(3,1) * qJDD(1) - pkin(7) * t142 - t136 * t294 + t296 * t138;
t312 = -mrSges(3,1) * t229 - pkin(2) * t358 - pkin(7) * t143 - t296 * t136 - t294 * t138;
t302 = -m(3) * t229 + t299 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t358;
t306 = -mrSges(2,2) * t271 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t317) + qJ(2) * t302 + mrSges(2,1) * t270 + Ifges(2,3) * qJDD(1) + t313;
t144 = m(2) * t271 - mrSges(2,1) * t299 - qJDD(1) * mrSges(2,2) + t302;
t141 = -m(3) * g(3) + t143;
t139 = m(2) * t270 - t299 * mrSges(2,2) + t348 * qJDD(1) + t317;
t135 = t344 * t299 + t345 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t270 - qJ(2) * t141 + t359;
t134 = mrSges(2,3) * t271 - pkin(1) * t141 + t348 * g(3) - t344 * qJDD(1) + t345 * t299 + t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t297 * t135 - t295 * t134 - pkin(6) * (t139 * t297 + t144 * t295), t135, t313, t138, t357, t153, -t195 * t255 + t318; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t295 * t135 + t297 * t134 + pkin(6) * (-t139 * t295 + t144 * t297), t134, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t299 * Ifges(3,5) - t359, t136, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t261 + Ifges(5,6) * t260 - qJD(3) * t236 + t238 * t338 + t301, t151, -t315; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t306, t306, mrSges(3,2) * g(3) + t299 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t312, t300, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t261 + Ifges(5,3) * t260 + qJD(3) * t237 + t238 * t337 + t308, -t304, Ifges(7,5) * t210 + Ifges(7,6) * t254 + Ifges(7,3) * t209 + t256 * t195 - t276 * t197 - t327;];
m_new  = t1;
