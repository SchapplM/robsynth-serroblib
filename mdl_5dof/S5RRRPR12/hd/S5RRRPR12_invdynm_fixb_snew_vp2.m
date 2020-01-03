% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:24
% EndTime: 2019-12-31 21:38:00
% DurationCPUTime: 19.16s
% Computational Cost: add. (333746->324), mult. (731779->428), div. (0->0), fcn. (569003->12), ass. (0->136)
t268 = sin(pkin(5));
t273 = sin(qJ(2));
t276 = cos(qJ(2));
t290 = qJD(1) * qJD(2);
t253 = (-qJDD(1) * t276 + t273 * t290) * t268;
t300 = cos(qJ(3));
t299 = pkin(7) * t268;
t270 = cos(pkin(5));
t298 = t270 * g(3);
t297 = t268 * t273;
t296 = t268 * t276;
t295 = t270 * t273;
t294 = t270 * t276;
t292 = qJD(1) * t268;
t251 = (-pkin(2) * t276 - pkin(8) * t273) * t292;
t263 = t270 * qJD(1) + qJD(2);
t261 = t263 ^ 2;
t262 = t270 * qJDD(1) + qJDD(2);
t291 = qJD(1) * t276;
t274 = sin(qJ(1));
t277 = cos(qJ(1));
t259 = t274 * g(1) - t277 * g(2);
t278 = qJD(1) ^ 2;
t248 = qJDD(1) * pkin(1) + t278 * t299 + t259;
t260 = -t277 * g(1) - t274 * g(2);
t249 = -t278 * pkin(1) + qJDD(1) * t299 + t260;
t293 = t248 * t295 + t276 * t249;
t202 = -t261 * pkin(2) + t262 * pkin(8) + (-g(3) * t273 + t251 * t291) * t268 + t293;
t252 = (qJDD(1) * t273 + t276 * t290) * t268;
t203 = t253 * pkin(2) - t252 * pkin(8) - t298 + (-t248 + (pkin(2) * t273 - pkin(8) * t276) * t263 * qJD(1)) * t268;
t272 = sin(qJ(3));
t186 = t300 * t202 + t272 * t203;
t289 = t273 * t292;
t241 = -t300 * t263 + t272 * t289;
t242 = t272 * t263 + t300 * t289;
t225 = t241 * pkin(3) - t242 * qJ(4);
t245 = qJDD(3) + t253;
t288 = t268 * t291;
t258 = qJD(3) - t288;
t257 = t258 ^ 2;
t175 = -t257 * pkin(3) + t245 * qJ(4) - t241 * t225 + t186;
t223 = -g(3) * t296 + t248 * t294 - t273 * t249;
t201 = -t262 * pkin(2) - t261 * pkin(8) + t251 * t289 - t223;
t221 = t242 * qJD(3) + t272 * t252 - t300 * t262;
t222 = -t241 * qJD(3) + t300 * t252 + t272 * t262;
t179 = (t241 * t258 - t222) * qJ(4) + (t242 * t258 + t221) * pkin(3) + t201;
t267 = sin(pkin(10));
t269 = cos(pkin(10));
t231 = t269 * t242 + t267 * t258;
t170 = -0.2e1 * qJD(4) * t231 - t267 * t175 + t269 * t179;
t208 = t269 * t222 + t267 * t245;
t230 = -t267 * t242 + t269 * t258;
t168 = (t230 * t241 - t208) * pkin(9) + (t230 * t231 + t221) * pkin(4) + t170;
t171 = 0.2e1 * qJD(4) * t230 + t269 * t175 + t267 * t179;
t207 = -t267 * t222 + t269 * t245;
t213 = t241 * pkin(4) - t231 * pkin(9);
t229 = t230 ^ 2;
t169 = -t229 * pkin(4) + t207 * pkin(9) - t241 * t213 + t171;
t271 = sin(qJ(5));
t275 = cos(qJ(5));
t166 = t275 * t168 - t271 * t169;
t204 = t275 * t230 - t271 * t231;
t184 = t204 * qJD(5) + t271 * t207 + t275 * t208;
t205 = t271 * t230 + t275 * t231;
t191 = -t204 * mrSges(6,1) + t205 * mrSges(6,2);
t240 = qJD(5) + t241;
t192 = -t240 * mrSges(6,2) + t204 * mrSges(6,3);
t219 = qJDD(5) + t221;
t162 = m(6) * t166 + t219 * mrSges(6,1) - t184 * mrSges(6,3) - t205 * t191 + t240 * t192;
t167 = t271 * t168 + t275 * t169;
t183 = -t205 * qJD(5) + t275 * t207 - t271 * t208;
t193 = t240 * mrSges(6,1) - t205 * mrSges(6,3);
t163 = m(6) * t167 - t219 * mrSges(6,2) + t183 * mrSges(6,3) + t204 * t191 - t240 * t193;
t154 = t275 * t162 + t271 * t163;
t209 = -t230 * mrSges(5,1) + t231 * mrSges(5,2);
t211 = -t241 * mrSges(5,2) + t230 * mrSges(5,3);
t152 = m(5) * t170 + t221 * mrSges(5,1) - t208 * mrSges(5,3) - t231 * t209 + t241 * t211 + t154;
t212 = t241 * mrSges(5,1) - t231 * mrSges(5,3);
t286 = -t271 * t162 + t275 * t163;
t153 = m(5) * t171 - t221 * mrSges(5,2) + t207 * mrSges(5,3) + t230 * t209 - t241 * t212 + t286;
t150 = -t267 * t152 + t269 * t153;
t226 = t241 * mrSges(4,1) + t242 * mrSges(4,2);
t233 = t258 * mrSges(4,1) - t242 * mrSges(4,3);
t148 = m(4) * t186 - t245 * mrSges(4,2) - t221 * mrSges(4,3) - t241 * t226 - t258 * t233 + t150;
t185 = -t272 * t202 + t300 * t203;
t174 = -t245 * pkin(3) - t257 * qJ(4) + t242 * t225 + qJDD(4) - t185;
t172 = -t207 * pkin(4) - t229 * pkin(9) + t231 * t213 + t174;
t284 = m(6) * t172 - t183 * mrSges(6,1) + t184 * mrSges(6,2) - t204 * t192 + t205 * t193;
t164 = -m(5) * t174 + t207 * mrSges(5,1) - t208 * mrSges(5,2) + t230 * t211 - t231 * t212 - t284;
t232 = -t258 * mrSges(4,2) - t241 * mrSges(4,3);
t158 = m(4) * t185 + t245 * mrSges(4,1) - t222 * mrSges(4,3) - t242 * t226 + t258 * t232 + t164;
t141 = t272 * t148 + t300 * t158;
t224 = -g(3) * t297 + t293;
t246 = t263 * mrSges(3,1) - mrSges(3,3) * t289;
t250 = (-mrSges(3,1) * t276 + mrSges(3,2) * t273) * t292;
t287 = t300 * t148 - t272 * t158;
t139 = m(3) * t224 - t262 * mrSges(3,2) - t253 * mrSges(3,3) - t263 * t246 + t250 * t288 + t287;
t247 = -t263 * mrSges(3,2) + mrSges(3,3) * t288;
t149 = t269 * t152 + t267 * t153;
t281 = -m(4) * t201 - t221 * mrSges(4,1) - t222 * mrSges(4,2) - t241 * t232 - t242 * t233 - t149;
t145 = m(3) * t223 + t262 * mrSges(3,1) - t252 * mrSges(3,3) + t263 * t247 - t250 * t289 + t281;
t135 = t276 * t139 - t273 * t145;
t237 = -t268 * t248 - t298;
t140 = m(3) * t237 + t253 * mrSges(3,1) + t252 * mrSges(3,2) + (t246 * t273 - t247 * t276) * t292 + t141;
t131 = t139 * t295 - t268 * t140 + t145 * t294;
t187 = Ifges(6,5) * t205 + Ifges(6,6) * t204 + Ifges(6,3) * t240;
t189 = Ifges(6,1) * t205 + Ifges(6,4) * t204 + Ifges(6,5) * t240;
t155 = -mrSges(6,1) * t172 + mrSges(6,3) * t167 + Ifges(6,4) * t184 + Ifges(6,2) * t183 + Ifges(6,6) * t219 - t205 * t187 + t240 * t189;
t188 = Ifges(6,4) * t205 + Ifges(6,2) * t204 + Ifges(6,6) * t240;
t156 = mrSges(6,2) * t172 - mrSges(6,3) * t166 + Ifges(6,1) * t184 + Ifges(6,4) * t183 + Ifges(6,5) * t219 + t204 * t187 - t240 * t188;
t194 = Ifges(5,5) * t231 + Ifges(5,6) * t230 + Ifges(5,3) * t241;
t196 = Ifges(5,1) * t231 + Ifges(5,4) * t230 + Ifges(5,5) * t241;
t142 = -mrSges(5,1) * t174 + mrSges(5,3) * t171 + Ifges(5,4) * t208 + Ifges(5,2) * t207 + Ifges(5,6) * t221 - pkin(4) * t284 + pkin(9) * t286 + t275 * t155 + t271 * t156 - t231 * t194 + t241 * t196;
t195 = Ifges(5,4) * t231 + Ifges(5,2) * t230 + Ifges(5,6) * t241;
t143 = mrSges(5,2) * t174 - mrSges(5,3) * t170 + Ifges(5,1) * t208 + Ifges(5,4) * t207 + Ifges(5,5) * t221 - pkin(9) * t154 - t271 * t155 + t275 * t156 + t230 * t194 - t241 * t195;
t215 = Ifges(4,5) * t242 - Ifges(4,6) * t241 + Ifges(4,3) * t258;
t216 = Ifges(4,4) * t242 - Ifges(4,2) * t241 + Ifges(4,6) * t258;
t132 = mrSges(4,2) * t201 - mrSges(4,3) * t185 + Ifges(4,1) * t222 - Ifges(4,4) * t221 + Ifges(4,5) * t245 - qJ(4) * t149 - t267 * t142 + t269 * t143 - t241 * t215 - t258 * t216;
t217 = Ifges(4,1) * t242 - Ifges(4,4) * t241 + Ifges(4,5) * t258;
t282 = -mrSges(6,1) * t166 + mrSges(6,2) * t167 - Ifges(6,5) * t184 - Ifges(6,6) * t183 - Ifges(6,3) * t219 - t205 * t188 + t204 * t189;
t280 = -mrSges(5,1) * t170 + mrSges(5,2) * t171 - Ifges(5,5) * t208 - Ifges(5,6) * t207 - pkin(4) * t154 - t231 * t195 + t230 * t196 + t282;
t136 = t280 + (-Ifges(4,2) - Ifges(5,3)) * t221 + t258 * t217 + Ifges(4,6) * t245 - t242 * t215 + Ifges(4,4) * t222 - mrSges(4,1) * t201 + mrSges(4,3) * t186 - pkin(3) * t149;
t235 = Ifges(3,6) * t263 + (Ifges(3,4) * t273 + Ifges(3,2) * t276) * t292;
t236 = Ifges(3,5) * t263 + (Ifges(3,1) * t273 + Ifges(3,4) * t276) * t292;
t123 = Ifges(3,5) * t252 - Ifges(3,6) * t253 + Ifges(3,3) * t262 + mrSges(3,1) * t223 - mrSges(3,2) * t224 + t272 * t132 + t300 * t136 + pkin(2) * t281 + pkin(8) * t287 + (t235 * t273 - t236 * t276) * t292;
t234 = Ifges(3,3) * t263 + (Ifges(3,5) * t273 + Ifges(3,6) * t276) * t292;
t125 = mrSges(3,2) * t237 - mrSges(3,3) * t223 + Ifges(3,1) * t252 - Ifges(3,4) * t253 + Ifges(3,5) * t262 - pkin(8) * t141 + t300 * t132 - t272 * t136 + t234 * t288 - t263 * t235;
t279 = mrSges(4,1) * t185 - mrSges(4,2) * t186 + Ifges(4,5) * t222 - Ifges(4,6) * t221 + Ifges(4,3) * t245 + pkin(3) * t164 + qJ(4) * t150 + t269 * t142 + t267 * t143 + t242 * t216 + t241 * t217;
t127 = -mrSges(3,1) * t237 + mrSges(3,3) * t224 + Ifges(3,4) * t252 - Ifges(3,2) * t253 + Ifges(3,6) * t262 - pkin(2) * t141 - t234 * t289 + t263 * t236 - t279;
t283 = mrSges(2,1) * t259 - mrSges(2,2) * t260 + Ifges(2,3) * qJDD(1) + pkin(1) * t131 + t270 * t123 + t125 * t297 + t127 * t296 + t135 * t299;
t133 = m(2) * t260 - t278 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t135;
t130 = t270 * t140 + (t139 * t273 + t145 * t276) * t268;
t128 = m(2) * t259 + qJDD(1) * mrSges(2,1) - t278 * mrSges(2,2) + t131;
t121 = -mrSges(2,2) * g(3) - mrSges(2,3) * t259 + Ifges(2,5) * qJDD(1) - t278 * Ifges(2,6) + t276 * t125 - t273 * t127 + (-t130 * t268 - t131 * t270) * pkin(7);
t120 = mrSges(2,1) * g(3) + mrSges(2,3) * t260 + t278 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t268 * t123 + (pkin(7) * t135 + t125 * t273 + t127 * t276) * t270;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t277 * t121 - t274 * t120 - pkin(6) * (t277 * t128 + t274 * t133), t121, t125, t132, t143, t156; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t274 * t121 + t277 * t120 + pkin(6) * (-t274 * t128 + t277 * t133), t120, t127, t136, t142, t155; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t283, t283, t123, t279, Ifges(5,3) * t221 - t280, -t282;];
m_new = t1;
