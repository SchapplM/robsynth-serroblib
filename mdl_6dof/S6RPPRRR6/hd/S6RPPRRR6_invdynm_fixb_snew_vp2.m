% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:57:06
% EndTime: 2019-05-05 15:57:17
% DurationCPUTime: 5.19s
% Computational Cost: add. (77321->299), mult. (145793->351), div. (0->0), fcn. (84050->8), ass. (0->116)
t278 = sin(qJ(1));
t282 = cos(qJ(1));
t247 = t278 * g(1) - t282 * g(2);
t285 = qJD(1) ^ 2;
t223 = -qJDD(1) * pkin(1) - t285 * qJ(2) + qJDD(2) - t247;
t214 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t223;
t248 = -t282 * g(1) - t278 * g(2);
t316 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t248;
t315 = mrSges(3,2) - mrSges(4,3);
t314 = Ifges(3,4) - Ifges(4,5);
t313 = Ifges(2,6) - Ifges(3,5);
t312 = mrSges(4,3) * t285;
t215 = qJDD(3) + (-pkin(1) - qJ(3)) * t285 + t316;
t210 = -qJDD(1) * pkin(7) + t215;
t277 = sin(qJ(4));
t281 = cos(qJ(4));
t203 = -g(3) * t281 + t277 * t210;
t240 = (mrSges(5,1) * t277 + mrSges(5,2) * t281) * qJD(1);
t310 = qJD(1) * qJD(4);
t251 = t281 * t310;
t242 = -t277 * qJDD(1) - t251;
t311 = qJD(1) * t281;
t246 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t311;
t253 = t277 * qJD(1);
t209 = -t285 * pkin(7) - t214;
t307 = t277 * t310;
t243 = qJDD(1) * t281 - t307;
t187 = (-t243 + t307) * pkin(8) + (-t242 + t251) * pkin(4) + t209;
t241 = (pkin(4) * t277 - pkin(8) * t281) * qJD(1);
t284 = qJD(4) ^ 2;
t190 = -pkin(4) * t284 + qJDD(4) * pkin(8) - t241 * t253 + t203;
t276 = sin(qJ(5));
t280 = cos(qJ(5));
t172 = t280 * t187 - t276 * t190;
t238 = qJD(4) * t280 - t276 * t311;
t201 = qJD(5) * t238 + qJDD(4) * t276 + t243 * t280;
t237 = qJDD(5) - t242;
t239 = qJD(4) * t276 + t280 * t311;
t250 = t253 + qJD(5);
t169 = (t238 * t250 - t201) * pkin(9) + (t238 * t239 + t237) * pkin(5) + t172;
t173 = t276 * t187 + t280 * t190;
t200 = -qJD(5) * t239 + qJDD(4) * t280 - t243 * t276;
t218 = pkin(5) * t250 - pkin(9) * t239;
t236 = t238 ^ 2;
t170 = -pkin(5) * t236 + pkin(9) * t200 - t218 * t250 + t173;
t275 = sin(qJ(6));
t279 = cos(qJ(6));
t167 = t169 * t279 - t170 * t275;
t204 = t238 * t279 - t239 * t275;
t179 = qJD(6) * t204 + t200 * t275 + t201 * t279;
t205 = t238 * t275 + t239 * t279;
t186 = -mrSges(7,1) * t204 + mrSges(7,2) * t205;
t249 = qJD(6) + t250;
t191 = -mrSges(7,2) * t249 + mrSges(7,3) * t204;
t229 = qJDD(6) + t237;
t162 = m(7) * t167 + mrSges(7,1) * t229 - mrSges(7,3) * t179 - t186 * t205 + t191 * t249;
t168 = t169 * t275 + t170 * t279;
t178 = -qJD(6) * t205 + t200 * t279 - t201 * t275;
t192 = mrSges(7,1) * t249 - mrSges(7,3) * t205;
t163 = m(7) * t168 - mrSges(7,2) * t229 + mrSges(7,3) * t178 + t186 * t204 - t192 * t249;
t154 = t279 * t162 + t275 * t163;
t207 = -mrSges(6,1) * t238 + mrSges(6,2) * t239;
t216 = -mrSges(6,2) * t250 + mrSges(6,3) * t238;
t152 = m(6) * t172 + mrSges(6,1) * t237 - mrSges(6,3) * t201 - t207 * t239 + t216 * t250 + t154;
t217 = mrSges(6,1) * t250 - mrSges(6,3) * t239;
t304 = -t162 * t275 + t279 * t163;
t153 = m(6) * t173 - mrSges(6,2) * t237 + mrSges(6,3) * t200 + t207 * t238 - t217 * t250 + t304;
t305 = -t152 * t276 + t280 * t153;
t145 = m(5) * t203 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t242 - qJD(4) * t246 - t240 * t253 + t305;
t202 = g(3) * t277 + t210 * t281;
t245 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t253;
t189 = -qJDD(4) * pkin(4) - pkin(8) * t284 + t241 * t311 - t202;
t171 = -pkin(5) * t200 - pkin(9) * t236 + t218 * t239 + t189;
t297 = m(7) * t171 - t178 * mrSges(7,1) + mrSges(7,2) * t179 - t204 * t191 + t192 * t205;
t289 = -m(6) * t189 + t200 * mrSges(6,1) - mrSges(6,2) * t201 + t238 * t216 - t217 * t239 - t297;
t158 = m(5) * t202 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t243 + qJD(4) * t245 - t240 * t311 + t289;
t134 = t277 * t145 + t281 * t158;
t147 = t280 * t152 + t276 * t153;
t306 = t281 * t145 - t277 * t158;
t303 = m(4) * t215 + qJDD(1) * mrSges(4,2) + t134;
t300 = -m(5) * t209 + t242 * mrSges(5,1) - t243 * mrSges(5,2) - t245 * t253 - t246 * t311 - t147;
t180 = Ifges(7,5) * t205 + Ifges(7,6) * t204 + Ifges(7,3) * t249;
t182 = Ifges(7,1) * t205 + Ifges(7,4) * t204 + Ifges(7,5) * t249;
t155 = -mrSges(7,1) * t171 + mrSges(7,3) * t168 + Ifges(7,4) * t179 + Ifges(7,2) * t178 + Ifges(7,6) * t229 - t180 * t205 + t182 * t249;
t181 = Ifges(7,4) * t205 + Ifges(7,2) * t204 + Ifges(7,6) * t249;
t156 = mrSges(7,2) * t171 - mrSges(7,3) * t167 + Ifges(7,1) * t179 + Ifges(7,4) * t178 + Ifges(7,5) * t229 + t180 * t204 - t181 * t249;
t193 = Ifges(6,5) * t239 + Ifges(6,6) * t238 + Ifges(6,3) * t250;
t195 = Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t250;
t127 = -mrSges(6,1) * t189 + mrSges(6,3) * t173 + Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t237 - pkin(5) * t297 + pkin(9) * t304 + t279 * t155 + t275 * t156 - t239 * t193 + t250 * t195;
t194 = Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t250;
t137 = mrSges(6,2) * t189 - mrSges(6,3) * t172 + Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t237 - pkin(9) * t154 - t155 * t275 + t156 * t279 + t193 * t238 - t194 * t250;
t226 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t281 - Ifges(5,6) * t277) * qJD(1);
t227 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t281 - Ifges(5,2) * t277) * qJD(1);
t123 = mrSges(5,2) * t209 - mrSges(5,3) * t202 + Ifges(5,1) * t243 + Ifges(5,4) * t242 + Ifges(5,5) * qJDD(4) - pkin(8) * t147 - qJD(4) * t227 - t127 * t276 + t137 * t280 - t226 * t253;
t228 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t281 - Ifges(5,4) * t277) * qJD(1);
t296 = -mrSges(7,1) * t167 + mrSges(7,2) * t168 - Ifges(7,5) * t179 - Ifges(7,6) * t178 - Ifges(7,3) * t229 - t205 * t181 + t204 * t182;
t286 = mrSges(6,1) * t172 - mrSges(6,2) * t173 + Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t237 + pkin(5) * t154 + t239 * t194 - t238 * t195 - t296;
t125 = -mrSges(5,1) * t209 + mrSges(5,3) * t203 + Ifges(5,4) * t243 + Ifges(5,2) * t242 + Ifges(5,6) * qJDD(4) - pkin(4) * t147 + qJD(4) * t228 - t226 * t311 - t286;
t299 = mrSges(4,1) * t214 + mrSges(4,2) * g(3) + t285 * Ifges(4,4) + Ifges(4,5) * qJDD(1) + pkin(3) * t300 + pkin(7) * t306 + t277 * t123 + t281 * t125;
t221 = pkin(1) * t285 - t316;
t298 = -m(3) * t221 + t285 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t303;
t295 = mrSges(4,2) * t215 - mrSges(4,3) * t214 + Ifges(4,1) * qJDD(1) - pkin(7) * t134 + t281 * t123 - t125 * t277;
t294 = mrSges(5,1) * t202 - mrSges(5,2) * t203 + Ifges(5,5) * t243 + Ifges(5,6) * t242 + Ifges(5,3) * qJDD(4) + pkin(4) * t289 + pkin(8) * t305 + t280 * t127 + t276 * t137 + t227 * t311 + t228 * t253;
t140 = m(4) * t214 - t285 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t300;
t293 = mrSges(3,1) * t223 + pkin(2) * t140 + t299;
t292 = -m(3) * t223 + t285 * mrSges(3,3) - t140;
t291 = mrSges(4,1) * t215 - Ifges(4,4) * qJDD(1) + pkin(3) * t134 + t294;
t290 = mrSges(3,2) * t223 - mrSges(3,3) * t221 + Ifges(3,1) * qJDD(1) - qJ(3) * t140 + t295;
t288 = -mrSges(2,2) * t248 + qJ(2) * (t298 - t312) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t292) + mrSges(2,1) * t247 + Ifges(2,3) * qJDD(1) + t290;
t287 = mrSges(3,1) * t221 + pkin(2) * (-t303 + t312) + qJ(3) * (-m(4) * g(3) + t306) - t291;
t138 = t292 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) - t285 * mrSges(2,2) + m(2) * t247;
t131 = (-m(3) - m(4)) * g(3) + t306;
t128 = m(2) * t248 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t285 + t298;
t120 = -t287 - pkin(1) * t131 + (Ifges(2,5) - t314) * t285 + t313 * qJDD(1) + (mrSges(2,1) - t315) * g(3) + mrSges(2,3) * t248;
t119 = t293 - qJ(2) * t131 - t313 * t285 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) - mrSges(2,3) * t247;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t282 * t119 - t278 * t120 - pkin(6) * (t128 * t278 + t138 * t282), t119, t290, t295, t123, t137, t156; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t278 * t119 + t282 * t120 + pkin(6) * (t128 * t282 - t138 * t278), t120, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t285 * Ifges(3,5) - t293, -mrSges(4,3) * g(3) - t285 * Ifges(4,5) - t291, t125, t127, t155; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t288, t288, Ifges(3,5) * qJDD(1) + t315 * g(3) + t314 * t285 + t287, t299, t294, t286, -t296;];
m_new  = t1;
