% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:50:49
% EndTime: 2019-05-05 21:50:58
% DurationCPUTime: 3.50s
% Computational Cost: add. (41322->344), mult. (77880->378), div. (0->0), fcn. (43792->6), ass. (0->126)
t292 = sin(qJ(1));
t294 = cos(qJ(1));
t276 = -g(1) * t294 - g(2) * t292;
t344 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t276;
t296 = qJD(1) ^ 2;
t340 = (-pkin(1) - pkin(7));
t239 = (t340 * t296) - t344;
t291 = sin(qJ(3));
t293 = cos(qJ(3));
t324 = qJD(1) * qJD(3);
t320 = t293 * t324;
t270 = -qJDD(1) * t291 - t320;
t321 = t291 * t324;
t271 = qJDD(1) * t293 - t321;
t181 = (-t271 + t321) * pkin(8) + (-t270 + t320) * pkin(3) + t239;
t275 = g(1) * t292 - t294 * g(2);
t316 = -qJ(2) * t296 + qJDD(2) - t275;
t241 = t340 * qJDD(1) + t316;
t229 = -g(3) * t293 + t291 * t241;
t269 = (pkin(3) * t291 - pkin(8) * t293) * qJD(1);
t295 = qJD(3) ^ 2;
t326 = qJD(1) * t291;
t188 = -pkin(3) * t295 + qJDD(3) * pkin(8) - t269 * t326 + t229;
t290 = sin(qJ(4));
t339 = cos(qJ(4));
t178 = t339 * t181 - t290 * t188;
t325 = qJD(1) * t293;
t266 = -t339 * qJD(3) + t290 * t325;
t267 = t290 * qJD(3) + t339 * t325;
t225 = pkin(4) * t266 - qJ(5) * t267;
t265 = qJDD(4) - t270;
t278 = qJD(4) + t326;
t277 = t278 ^ 2;
t176 = -t265 * pkin(4) - t277 * qJ(5) + t267 * t225 + qJDD(5) - t178;
t219 = -t266 * qJD(4) + t290 * qJDD(3) + t339 * t271;
t227 = -mrSges(6,2) * t266 - mrSges(6,3) * t267;
t343 = -m(6) * t176 - t219 * mrSges(6,1) - t267 * t227;
t224 = -mrSges(7,2) * t267 + mrSges(7,3) * t266;
t332 = t266 * t278;
t167 = -0.2e1 * qJD(6) * t278 + (t266 * t267 - t265) * qJ(6) + (t219 + t332) * pkin(5) + t176;
t235 = -mrSges(7,1) * t266 + mrSges(7,2) * t278;
t318 = -m(7) * t167 + t265 * mrSges(7,3) + t278 * t235;
t164 = mrSges(7,1) * t219 + t224 * t267 - t318;
t179 = t290 * t181 + t339 * t188;
t195 = Ifges(7,4) * t278 + Ifges(7,2) * t266 + Ifges(7,6) * t267;
t197 = Ifges(5,4) * t267 - Ifges(5,2) * t266 + Ifges(5,6) * t278;
t218 = qJD(4) * t267 - t339 * qJDD(3) + t271 * t290;
t234 = mrSges(6,1) * t266 - mrSges(6,3) * t278;
t309 = -pkin(4) * t277 + qJ(5) * t265 - t225 * t266 + t179;
t341 = -2 * qJD(5);
t174 = t278 * t341 - t309;
t193 = Ifges(6,5) * t278 - Ifges(6,6) * t267 + Ifges(6,3) * t266;
t232 = pkin(5) * t267 - qJ(6) * t278;
t264 = t266 ^ 2;
t173 = -pkin(5) * t218 - qJ(6) * t264 + qJDD(6) + ((2 * qJD(5)) + t232) * t278 + t309;
t192 = Ifges(7,5) * t278 + Ifges(7,6) * t266 + Ifges(7,3) * t267;
t312 = -mrSges(7,2) * t173 + mrSges(7,3) * t167 - Ifges(7,1) * t265 - Ifges(7,4) * t218 - Ifges(7,5) * t219 - t266 * t192;
t300 = -mrSges(6,2) * t176 + mrSges(6,3) * t174 - Ifges(6,1) * t265 + Ifges(6,4) * t219 - Ifges(6,5) * t218 + qJ(6) * t164 + t267 * t193 + t312;
t236 = mrSges(6,1) * t267 + mrSges(6,2) * t278;
t233 = mrSges(7,1) * t267 - mrSges(7,3) * t278;
t322 = m(7) * t173 + t265 * mrSges(7,2) + t278 * t233;
t315 = -m(6) * t174 + t265 * mrSges(6,3) + t278 * t236 + t322;
t327 = -t224 - t227;
t196 = Ifges(6,4) * t278 - Ifges(6,2) * t267 + Ifges(6,6) * t266;
t328 = Ifges(5,1) * t267 - Ifges(5,4) * t266 + Ifges(5,5) * t278 - t196;
t342 = t328 * t266 + (t197 - t195) * t267 + mrSges(5,1) * t178 - mrSges(5,2) * t179 + Ifges(5,5) * t219 - Ifges(5,6) * t218 + Ifges(5,3) * t265 + pkin(4) * (-mrSges(6,2) * t265 - t234 * t278 - t164 + t343) + qJ(5) * (t327 * t266 + (-mrSges(6,1) - mrSges(7,1)) * t218 + t315) - t300;
t337 = mrSges(2,1) - mrSges(3,2);
t336 = -mrSges(7,1) - mrSges(5,3);
t335 = -Ifges(3,4) + Ifges(2,5);
t334 = Ifges(5,4) + Ifges(6,6);
t333 = (Ifges(3,5) - Ifges(2,6));
t331 = t267 * t195;
t226 = mrSges(5,1) * t266 + mrSges(5,2) * t267;
t230 = -mrSges(5,2) * t278 - mrSges(5,3) * t266;
t154 = m(5) * t178 + (t230 - t234) * t278 + (-t224 - t226) * t267 + (mrSges(5,1) - mrSges(6,2)) * t265 + t336 * t219 + t318 + t343;
t231 = mrSges(5,1) * t278 - mrSges(5,3) * t267;
t157 = m(5) * t179 - mrSges(5,2) * t265 - t231 * t278 + (-t226 + t327) * t266 + (-mrSges(6,1) + t336) * t218 + t315;
t150 = t339 * t154 + t290 * t157;
t198 = Ifges(7,1) * t278 + Ifges(7,4) * t266 + Ifges(7,5) * t267;
t199 = Ifges(6,1) * t278 - Ifges(6,4) * t267 + Ifges(6,5) * t266;
t329 = t198 + t199;
t268 = (mrSges(4,1) * t291 + mrSges(4,2) * t293) * qJD(1);
t274 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t325;
t319 = -t154 * t290 + t339 * t157;
t148 = m(4) * t229 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t270 - qJD(3) * t274 - t268 * t326 + t319;
t228 = g(3) * t291 + t241 * t293;
t273 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t326;
t187 = -qJDD(3) * pkin(3) - pkin(8) * t295 + t269 * t325 - t228;
t299 = (-t219 + t332) * qJ(5) + t187 + (t278 * pkin(4) + t341) * t267;
t170 = -pkin(5) * t264 + 0.2e1 * qJD(6) * t266 - t232 * t267 + (pkin(4) + qJ(6)) * t218 + t299;
t163 = m(7) * t170 - t219 * mrSges(7,2) + t218 * mrSges(7,3) - t267 * t233 + t266 * t235;
t177 = pkin(4) * t218 + t299;
t310 = -m(6) * t177 + t218 * mrSges(6,2) + t266 * t234 - t163;
t298 = -m(5) * t187 - t218 * mrSges(5,1) - t266 * t230 + (-t231 + t236) * t267 + (-mrSges(5,2) + mrSges(6,3)) * t219 + t310;
t152 = m(4) * t228 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t271 + qJD(3) * t273 - t268 * t325 + t298;
t141 = t293 * t148 - t152 * t291;
t140 = t148 * t291 + t152 * t293;
t314 = mrSges(7,1) * t173 - mrSges(7,3) * t170 - Ifges(7,4) * t265 - Ifges(7,2) * t218 - Ifges(7,6) * t219 - t267 * t198;
t313 = -mrSges(7,1) * t167 + mrSges(7,2) * t170 - Ifges(7,5) * t265 - Ifges(7,6) * t218 - Ifges(7,3) * t219 - t278 * t195;
t246 = -qJDD(1) * pkin(1) + t316;
t311 = -m(3) * t246 + (t296 * mrSges(3,3)) - t140;
t146 = -m(4) * t239 + mrSges(4,1) * t270 - t271 * mrSges(4,2) - t273 * t326 - t274 * t325 - t150;
t160 = -mrSges(6,3) * t219 - t236 * t267 - t310;
t194 = Ifges(5,5) * t267 - Ifges(5,6) * t266 + Ifges(5,3) * t278;
t302 = mrSges(6,1) * t174 - mrSges(6,2) * t177 + pkin(5) * (mrSges(7,1) * t218 + t224 * t266 - t322) + qJ(6) * t163 - t314;
t136 = (t192 + t328) * t278 - t302 + (-t194 - t199) * t267 + (Ifges(5,6) - Ifges(6,5)) * t265 + t334 * t219 + (-Ifges(5,2) - Ifges(6,3)) * t218 + mrSges(5,3) * t179 - mrSges(5,1) * t187 - pkin(4) * t160;
t304 = -mrSges(6,1) * t176 + mrSges(6,3) * t177 - pkin(5) * t164 + t313;
t143 = -t304 + (-t197 + t193) * t278 + (-t194 - t329) * t266 + (Ifges(5,5) - Ifges(6,4)) * t265 + (Ifges(5,1) + Ifges(6,2)) * t219 - t334 * t218 - mrSges(5,3) * t178 + mrSges(5,2) * t187 - qJ(5) * t160;
t249 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t293 - Ifges(4,6) * t291) * qJD(1);
t250 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t293 - Ifges(4,2) * t291) * qJD(1);
t133 = mrSges(4,2) * t239 - mrSges(4,3) * t228 + Ifges(4,1) * t271 + Ifges(4,4) * t270 + Ifges(4,5) * qJDD(3) - pkin(8) * t150 - qJD(3) * t250 - t290 * t136 + t339 * t143 - t249 * t326;
t251 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t293 - Ifges(4,4) * t291) * qJD(1);
t134 = -mrSges(4,1) * t239 + mrSges(4,3) * t229 + Ifges(4,4) * t271 + Ifges(4,2) * t270 + Ifges(4,6) * qJDD(3) - pkin(3) * t150 + qJD(3) * t251 - t249 * t325 - t342;
t244 = pkin(1) * t296 + t344;
t308 = mrSges(3,2) * t246 - mrSges(3,3) * t244 + Ifges(3,1) * qJDD(1) - pkin(7) * t140 + t293 * t133 - t134 * t291;
t307 = -mrSges(3,1) * t244 - pkin(2) * t146 - pkin(7) * t141 - t133 * t291 - t134 * t293;
t306 = mrSges(4,1) * t228 - mrSges(4,2) * t229 + Ifges(4,5) * t271 + Ifges(4,6) * t270 + Ifges(4,3) * qJDD(3) + pkin(3) * t298 + pkin(8) * t319 + t339 * t136 + t290 * t143 + t250 * t325 + t251 * t326;
t305 = -m(3) * t244 + t296 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t146;
t303 = -mrSges(2,2) * t276 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t311) + qJ(2) * t305 + mrSges(2,1) * t275 + Ifges(2,3) * qJDD(1) + t308;
t301 = mrSges(3,1) * t246 + pkin(2) * t140 + t306;
t144 = m(2) * t276 - mrSges(2,1) * t296 - qJDD(1) * mrSges(2,2) + t305;
t139 = -m(3) * g(3) + t141;
t137 = m(2) * t275 - mrSges(2,2) * t296 + t337 * qJDD(1) + t311;
t131 = t301 + (t333 * t296) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t335 * qJDD(1) - mrSges(2,3) * t275 - qJ(2) * t139;
t130 = mrSges(2,3) * t276 - pkin(1) * t139 + t337 * g(3) - t333 * qJDD(1) + t335 * t296 + t307;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t294 * t131 - t292 * t130 - pkin(6) * (t137 * t294 + t144 * t292), t131, t308, t133, t143, -t266 * t196 - t300 - t331, -t312 - t331; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t292 * t131 + t294 * t130 + pkin(6) * (-t137 * t292 + t144 * t294), t130, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (Ifges(3,5) * t296) - t301, t134, t136, Ifges(6,4) * t265 - Ifges(6,2) * t219 + Ifges(6,6) * t218 - t278 * t193 + t329 * t266 + t304, -t278 * t192 - t314; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t303, t303, mrSges(3,2) * g(3) + Ifges(3,4) * t296 + Ifges(3,5) * qJDD(1) - t307, t306, t342, t302 + t267 * t199 + Ifges(6,5) * t265 + Ifges(6,3) * t218 - Ifges(6,6) * t219 + (-t192 + t196) * t278, -t266 * t198 - t313;];
m_new  = t1;
