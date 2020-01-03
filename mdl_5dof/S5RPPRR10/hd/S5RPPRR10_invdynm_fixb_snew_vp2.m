% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:03
% EndTime: 2019-12-31 18:04:09
% DurationCPUTime: 4.17s
% Computational Cost: add. (43159->262), mult. (103977->337), div. (0->0), fcn. (69380->8), ass. (0->112)
t272 = sin(qJ(1));
t275 = cos(qJ(1));
t239 = -g(1) * t275 - g(2) * t272;
t276 = qJD(1) ^ 2;
t315 = -pkin(1) * t276 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t239;
t268 = sin(pkin(8));
t261 = t268 ^ 2;
t269 = cos(pkin(8));
t262 = t269 ^ 2;
t304 = t262 * t276;
t314 = t261 * t276 + t304;
t212 = -t269 * g(3) - t268 * t315;
t228 = (-pkin(2) * t269 - qJ(3) * t268) * qJD(1);
t300 = qJD(1) * t268;
t192 = t228 * t300 + qJDD(3) - t212;
t183 = (-pkin(3) * t269 * t276 - pkin(6) * qJDD(1)) * t268 + t192;
t213 = -g(3) * t268 + t269 * t315;
t299 = qJD(1) * t269;
t194 = t228 * t299 + t213;
t295 = qJDD(1) * t269;
t185 = -pkin(3) * t304 - pkin(6) * t295 + t194;
t271 = sin(qJ(4));
t274 = cos(qJ(4));
t165 = t183 * t274 - t185 * t271;
t288 = t268 * t274 - t269 * t271;
t287 = -t268 * t271 - t269 * t274;
t225 = t287 * qJD(1);
t298 = qJD(4) * t225;
t211 = qJDD(1) * t288 + t298;
t226 = t288 * qJD(1);
t159 = (-t211 + t298) * pkin(7) + (t225 * t226 + qJDD(4)) * pkin(4) + t165;
t166 = t183 * t271 + t185 * t274;
t210 = -qJD(4) * t226 + qJDD(1) * t287;
t216 = qJD(4) * pkin(4) - pkin(7) * t226;
t224 = t225 ^ 2;
t160 = -pkin(4) * t224 + pkin(7) * t210 - qJD(4) * t216 + t166;
t270 = sin(qJ(5));
t273 = cos(qJ(5));
t157 = t159 * t273 - t160 * t270;
t200 = t225 * t273 - t226 * t270;
t174 = qJD(5) * t200 + t210 * t270 + t211 * t273;
t201 = t225 * t270 + t226 * t273;
t180 = -mrSges(6,1) * t200 + mrSges(6,2) * t201;
t263 = qJD(4) + qJD(5);
t195 = -mrSges(6,2) * t263 + mrSges(6,3) * t200;
t260 = qJDD(4) + qJDD(5);
t153 = m(6) * t157 + mrSges(6,1) * t260 - mrSges(6,3) * t174 - t180 * t201 + t195 * t263;
t158 = t159 * t270 + t160 * t273;
t173 = -qJD(5) * t201 + t210 * t273 - t211 * t270;
t196 = mrSges(6,1) * t263 - mrSges(6,3) * t201;
t154 = m(6) * t158 - mrSges(6,2) * t260 + mrSges(6,3) * t173 + t180 * t200 - t196 * t263;
t143 = t153 * t273 + t154 * t270;
t204 = -mrSges(5,1) * t225 + mrSges(5,2) * t226;
t214 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t225;
t140 = m(5) * t165 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t211 + qJD(4) * t214 - t204 * t226 + t143;
t215 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t226;
t291 = -t153 * t270 + t154 * t273;
t141 = m(5) * t166 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t210 - qJD(4) * t215 + t204 * t225 + t291;
t138 = -t140 * t271 + t141 * t274;
t229 = (-mrSges(4,1) * t269 - mrSges(4,3) * t268) * qJD(1);
t136 = m(4) * t194 + mrSges(4,2) * t295 + t229 * t299 + t138;
t290 = Ifges(4,5) * t268 - Ifges(4,3) * t269;
t231 = t290 * qJD(1);
t235 = (Ifges(4,1) * t268 - Ifges(4,5) * t269) * qJD(1);
t137 = t140 * t274 + t141 * t271;
t198 = Ifges(5,4) * t226 + Ifges(5,2) * t225 + Ifges(5,6) * qJD(4);
t199 = Ifges(5,1) * t226 + Ifges(5,4) * t225 + Ifges(5,5) * qJD(4);
t176 = Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t263;
t177 = Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t263;
t286 = mrSges(6,1) * t157 - mrSges(6,2) * t158 + Ifges(6,5) * t174 + Ifges(6,6) * t173 + Ifges(6,3) * t260 + t176 * t201 - t177 * t200;
t280 = mrSges(5,1) * t165 - mrSges(5,2) * t166 + Ifges(5,5) * t211 + Ifges(5,6) * t210 + Ifges(5,3) * qJDD(4) + pkin(4) * t143 + t198 * t226 - t199 * t225 + t286;
t296 = qJDD(1) * t268;
t278 = -mrSges(4,1) * t192 + mrSges(4,3) * t194 + Ifges(4,4) * t296 - pkin(3) * t137 - t280;
t284 = -m(4) * t192 - t137;
t307 = Ifges(3,1) * t268;
t313 = ((t231 - (Ifges(3,4) * t268 + Ifges(3,2) * t269) * qJD(1)) * t268 + (t235 + (Ifges(3,4) * t269 + t307) * qJD(1)) * t269) * qJD(1) - mrSges(3,1) * t212 + mrSges(3,2) * t213 - pkin(2) * ((-qJDD(1) * mrSges(4,2) - qJD(1) * t229) * t268 + t284) - qJ(3) * t136 - t278;
t238 = t272 * g(1) - t275 * g(2);
t223 = -qJDD(1) * pkin(1) - t276 * qJ(2) + qJDD(2) - t238;
t209 = -pkin(2) * t295 - qJ(3) * t296 - 0.2e1 * qJD(3) * t300 + t223;
t306 = Ifges(3,5) * t268;
t312 = (Ifges(3,6) - Ifges(4,6)) * t269 + t306;
t310 = Ifges(3,4) - Ifges(4,5);
t308 = mrSges(3,2) * t268;
t230 = (-mrSges(3,1) * t269 + t308) * qJD(1);
t133 = m(3) * t212 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t229 - t230) * qJD(1)) * t268 + t284;
t134 = m(3) * t213 + (qJDD(1) * mrSges(3,3) + qJD(1) * t230) * t269 + t136;
t292 = -t133 * t268 + t134 * t269;
t189 = pkin(3) * t295 + (-t261 - t262) * t276 * pkin(6) - t209;
t162 = -pkin(4) * t210 - pkin(7) * t224 + t216 * t226 + t189;
t289 = m(6) * t162 - mrSges(6,1) * t173 + mrSges(6,2) * t174 - t195 * t200 + t196 * t201;
t149 = -m(5) * t189 + mrSges(5,1) * t210 - mrSges(5,2) * t211 + t214 * t225 - t215 * t226 - t289;
t148 = m(4) * t209 - mrSges(4,1) * t295 - mrSges(4,2) * t314 - mrSges(4,3) * t296 + t149;
t232 = (Ifges(3,6) * t269 + t306) * qJD(1);
t233 = (Ifges(4,4) * t268 - Ifges(4,6) * t269) * qJD(1);
t175 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t263;
t144 = -mrSges(6,1) * t162 + mrSges(6,3) * t158 + Ifges(6,4) * t174 + Ifges(6,2) * t173 + Ifges(6,6) * t260 - t175 * t201 + t177 * t263;
t145 = mrSges(6,2) * t162 - mrSges(6,3) * t157 + Ifges(6,1) * t174 + Ifges(6,4) * t173 + Ifges(6,5) * t260 + t175 * t200 - t176 * t263;
t197 = Ifges(5,5) * t226 + Ifges(5,6) * t225 + Ifges(5,3) * qJD(4);
t129 = -mrSges(5,1) * t189 + mrSges(5,3) * t166 + Ifges(5,4) * t211 + Ifges(5,2) * t210 + Ifges(5,6) * qJDD(4) - pkin(4) * t289 + pkin(7) * t291 + qJD(4) * t199 + t273 * t144 + t270 * t145 - t226 * t197;
t131 = mrSges(5,2) * t189 - mrSges(5,3) * t165 + Ifges(5,1) * t211 + Ifges(5,4) * t210 + Ifges(5,5) * qJDD(4) - pkin(7) * t143 - qJD(4) * t198 - t144 * t270 + t145 * t273 + t197 * t225;
t281 = mrSges(4,1) * t209 - mrSges(4,2) * t194 + pkin(3) * t149 + pkin(6) * t138 + t129 * t274 + t131 * t271;
t123 = -mrSges(3,1) * t223 + mrSges(3,3) * t213 - pkin(2) * t148 + (-t232 - t233) * t300 + ((Ifges(3,2) + Ifges(4,3)) * t269 + t310 * t268) * qJDD(1) - t281;
t282 = mrSges(4,2) * t192 - mrSges(4,3) * t209 + Ifges(4,1) * t296 - pkin(6) * t137 - t129 * t271 + t131 * t274 + t233 * t299;
t125 = t232 * t299 + mrSges(3,2) * t223 - mrSges(3,3) * t212 - qJ(3) * t148 + (t269 * t310 + t307) * qJDD(1) + t282;
t279 = -m(3) * t223 + mrSges(3,1) * t295 + mrSges(3,3) * t314 - t148;
t283 = -mrSges(2,2) * t239 + qJ(2) * t292 + t269 * t123 + t268 * t125 + pkin(1) * (-mrSges(3,2) * t296 + t279) + mrSges(2,1) * t238 + Ifges(2,3) * qJDD(1);
t146 = t279 + m(2) * t238 - t276 * mrSges(2,2) + (mrSges(2,1) - t308) * qJDD(1);
t128 = t133 * t269 + t134 * t268;
t126 = m(2) * t239 - mrSges(2,1) * t276 - qJDD(1) * mrSges(2,2) + t292;
t121 = (Ifges(2,6) - t312) * qJDD(1) + mrSges(2,1) * g(3) + t276 * Ifges(2,5) + mrSges(2,3) * t239 - pkin(1) * t128 + t313;
t120 = -mrSges(2,2) * g(3) - mrSges(2,3) * t238 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t276 - qJ(2) * t128 - t123 * t268 + t125 * t269;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t275 * t120 - t272 * t121 - pkin(5) * (t126 * t272 + t146 * t275), t120, t125, -Ifges(4,5) * t295 + t282, t131, t145; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t272 * t120 + t275 * t121 + pkin(5) * (t126 * t275 - t146 * t272), t121, t123, -Ifges(4,6) * t295 + t278 + (-t268 * t231 - t269 * t235) * qJD(1), t129, t144; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t283, t283, qJDD(1) * t312 - t313, qJDD(1) * t290 + t233 * t300 + t281, t280, t286;];
m_new = t1;
