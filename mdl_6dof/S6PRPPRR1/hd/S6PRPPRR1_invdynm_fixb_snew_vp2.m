% Calculate vector of cutting torques with Newton-Euler for
% S6PRPPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-05-04 21:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:39:58
% EndTime: 2019-05-04 21:40:14
% DurationCPUTime: 14.56s
% Computational Cost: add. (274651->277), mult. (557444->356), div. (0->0), fcn. (411352->14), ass. (0->134)
t286 = qJD(2) ^ 2;
t271 = sin(pkin(12));
t275 = cos(pkin(12));
t280 = sin(qJ(5));
t283 = cos(qJ(5));
t296 = t271 * t280 - t275 * t283;
t246 = t296 * qJD(2);
t273 = sin(pkin(10));
t277 = cos(pkin(10));
t256 = g(1) * t273 - g(2) * t277;
t257 = -g(1) * t277 - g(2) * t273;
t270 = -g(3) + qJDD(1);
t281 = sin(qJ(2));
t278 = cos(pkin(6));
t284 = cos(qJ(2));
t314 = t278 * t284;
t274 = sin(pkin(6));
t316 = t274 * t284;
t225 = t256 * t314 - t257 * t281 + t270 * t316;
t220 = qJDD(2) * pkin(2) + t225;
t315 = t278 * t281;
t317 = t274 * t281;
t226 = t256 * t315 + t284 * t257 + t270 * t317;
t221 = -pkin(2) * t286 + t226;
t272 = sin(pkin(11));
t276 = cos(pkin(11));
t205 = t272 * t220 + t276 * t221;
t202 = -pkin(3) * t286 + qJDD(2) * qJ(4) + t205;
t244 = -t256 * t274 + t278 * t270;
t241 = qJDD(3) + t244;
t309 = qJD(2) * qJD(4);
t313 = t275 * t241 - 0.2e1 * t271 * t309;
t321 = pkin(4) * t275;
t195 = (-pkin(8) * qJDD(2) + t286 * t321 - t202) * t271 + t313;
t198 = t271 * t241 + (t202 + 0.2e1 * t309) * t275;
t308 = qJDD(2) * t275;
t267 = t275 ^ 2;
t318 = t267 * t286;
t196 = -pkin(4) * t318 + pkin(8) * t308 + t198;
t191 = t280 * t195 + t283 * t196;
t297 = t271 * t283 + t275 * t280;
t247 = t297 * qJD(2);
t228 = mrSges(6,1) * t246 + mrSges(6,2) * t247;
t310 = t247 * qJD(5);
t234 = -t296 * qJDD(2) - t310;
t243 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t247;
t233 = pkin(5) * t246 - pkin(9) * t247;
t285 = qJD(5) ^ 2;
t188 = -pkin(5) * t285 + qJDD(5) * pkin(9) - t233 * t246 + t191;
t266 = t271 ^ 2;
t204 = t276 * t220 - t272 * t221;
t299 = qJDD(4) - t204;
t199 = (-pkin(3) - t321) * qJDD(2) + (-qJ(4) + (-t266 - t267) * pkin(8)) * t286 + t299;
t311 = t246 * qJD(5);
t235 = t297 * qJDD(2) - t311;
t192 = (-t235 + t311) * pkin(9) + (-t234 + t310) * pkin(5) + t199;
t279 = sin(qJ(6));
t282 = cos(qJ(6));
t185 = -t188 * t279 + t192 * t282;
t238 = qJD(5) * t282 - t247 * t279;
t212 = qJD(6) * t238 + qJDD(5) * t279 + t235 * t282;
t239 = qJD(5) * t279 + t247 * t282;
t214 = -mrSges(7,1) * t238 + mrSges(7,2) * t239;
t245 = qJD(6) + t246;
t215 = -mrSges(7,2) * t245 + mrSges(7,3) * t238;
t232 = qJDD(6) - t234;
t181 = m(7) * t185 + mrSges(7,1) * t232 - mrSges(7,3) * t212 - t214 * t239 + t215 * t245;
t186 = t188 * t282 + t192 * t279;
t211 = -qJD(6) * t239 + qJDD(5) * t282 - t235 * t279;
t216 = mrSges(7,1) * t245 - mrSges(7,3) * t239;
t182 = m(7) * t186 - mrSges(7,2) * t232 + mrSges(7,3) * t211 + t214 * t238 - t216 * t245;
t303 = -t181 * t279 + t282 * t182;
t167 = m(6) * t191 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t234 - qJD(5) * t243 - t228 * t246 + t303;
t190 = t195 * t283 - t196 * t280;
t242 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t246;
t187 = -qJDD(5) * pkin(5) - pkin(9) * t285 + t233 * t247 - t190;
t292 = -m(7) * t187 + t211 * mrSges(7,1) - mrSges(7,2) * t212 + t238 * t215 - t216 * t239;
t177 = m(6) * t190 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t235 + qJD(5) * t242 - t228 * t247 + t292;
t162 = t280 * t167 + t283 * t177;
t197 = -t202 * t271 + t313;
t206 = Ifges(7,5) * t239 + Ifges(7,6) * t238 + Ifges(7,3) * t245;
t208 = Ifges(7,1) * t239 + Ifges(7,4) * t238 + Ifges(7,5) * t245;
t174 = -mrSges(7,1) * t187 + mrSges(7,3) * t186 + Ifges(7,4) * t212 + Ifges(7,2) * t211 + Ifges(7,6) * t232 - t206 * t239 + t208 * t245;
t207 = Ifges(7,4) * t239 + Ifges(7,2) * t238 + Ifges(7,6) * t245;
t175 = mrSges(7,2) * t187 - mrSges(7,3) * t185 + Ifges(7,1) * t212 + Ifges(7,4) * t211 + Ifges(7,5) * t232 + t206 * t238 - t207 * t245;
t223 = Ifges(6,4) * t247 - Ifges(6,2) * t246 + Ifges(6,6) * qJD(5);
t224 = Ifges(6,1) * t247 - Ifges(6,4) * t246 + Ifges(6,5) * qJD(5);
t290 = -mrSges(6,1) * t190 + mrSges(6,2) * t191 - Ifges(6,5) * t235 - Ifges(6,6) * t234 - Ifges(6,3) * qJDD(5) - pkin(5) * t292 - pkin(9) * t303 - t282 * t174 - t279 * t175 - t247 * t223 - t246 * t224;
t301 = Ifges(5,4) * t271 + Ifges(5,2) * t275;
t302 = Ifges(5,1) * t271 + Ifges(5,4) * t275;
t322 = -mrSges(5,1) * t197 + mrSges(5,2) * t198 - pkin(4) * t162 - (t271 * t301 - t275 * t302) * t286 + t290;
t319 = mrSges(5,2) * t271;
t295 = mrSges(5,3) * qJDD(2) + t286 * (-mrSges(5,1) * t275 + t319);
t160 = m(5) * t197 - t295 * t271 + t162;
t304 = t283 * t167 - t280 * t177;
t161 = m(5) * t198 + t295 * t275 + t304;
t305 = -t160 * t271 + t275 * t161;
t151 = m(4) * t205 - mrSges(4,1) * t286 - qJDD(2) * mrSges(4,2) + t305;
t201 = -qJDD(2) * pkin(3) - t286 * qJ(4) + t299;
t170 = t282 * t181 + t279 * t182;
t291 = m(6) * t199 - t234 * mrSges(6,1) + t235 * mrSges(6,2) + t246 * t242 + t247 * t243 + t170;
t289 = -m(5) * t201 + mrSges(5,1) * t308 - t291 + (t266 * t286 + t318) * mrSges(5,3);
t164 = t289 - t286 * mrSges(4,2) + m(4) * t204 + (mrSges(4,1) - t319) * qJDD(2);
t148 = t272 * t151 + t276 * t164;
t146 = m(3) * t225 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t286 + t148;
t306 = t276 * t151 - t164 * t272;
t147 = m(3) * t226 - mrSges(3,1) * t286 - qJDD(2) * mrSges(3,2) + t306;
t137 = -t146 * t281 + t284 * t147;
t320 = pkin(7) * t137;
t154 = t275 * t160 + t271 * t161;
t300 = Ifges(5,5) * t271 + Ifges(5,6) * t275;
t312 = t286 * t300;
t307 = m(4) * t241 + t154;
t152 = m(3) * t244 + t307;
t134 = t146 * t314 + t147 * t315 - t152 * t274;
t222 = Ifges(6,5) * t247 - Ifges(6,6) * t246 + Ifges(6,3) * qJD(5);
t155 = mrSges(6,2) * t199 - mrSges(6,3) * t190 + Ifges(6,1) * t235 + Ifges(6,4) * t234 + Ifges(6,5) * qJDD(5) - pkin(9) * t170 - qJD(5) * t223 - t174 * t279 + t175 * t282 - t222 * t246;
t288 = mrSges(7,1) * t185 - mrSges(7,2) * t186 + Ifges(7,5) * t212 + Ifges(7,6) * t211 + Ifges(7,3) * t232 + t207 * t239 - t208 * t238;
t156 = -mrSges(6,1) * t199 + mrSges(6,3) * t191 + Ifges(6,4) * t235 + Ifges(6,2) * t234 + Ifges(6,6) * qJDD(5) - pkin(5) * t170 + qJD(5) * t224 - t222 * t247 - t288;
t140 = -mrSges(5,1) * t201 + mrSges(5,3) * t198 - pkin(4) * t291 + pkin(8) * t304 + t301 * qJDD(2) + t280 * t155 + t283 * t156 - t271 * t312;
t142 = mrSges(5,2) * t201 - mrSges(5,3) * t197 - pkin(8) * t162 + t302 * qJDD(2) + t283 * t155 - t280 * t156 + t275 * t312;
t130 = mrSges(4,2) * t241 - mrSges(4,3) * t204 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t286 - qJ(4) * t154 - t140 * t271 + t142 * t275;
t138 = (Ifges(4,6) - t300) * qJDD(2) + t286 * Ifges(4,5) - mrSges(4,1) * t241 + mrSges(4,3) * t205 - pkin(3) * t154 + t322;
t125 = -mrSges(3,1) * t244 + mrSges(3,3) * t226 + t286 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t307 + qJ(3) * t306 + t272 * t130 + t276 * t138;
t127 = mrSges(3,2) * t244 - mrSges(3,3) * t225 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t286 - qJ(3) * t148 + t130 * t276 - t138 * t272;
t293 = -mrSges(4,2) * t205 + qJ(4) * t305 + t275 * t140 + t271 * t142 + pkin(3) * (-qJDD(2) * t319 + t289) + mrSges(4,1) * t204 + Ifges(4,3) * qJDD(2);
t129 = mrSges(3,1) * t225 - mrSges(3,2) * t226 + Ifges(3,3) * qJDD(2) + pkin(2) * t148 + t293;
t294 = mrSges(2,1) * t256 - mrSges(2,2) * t257 + pkin(1) * t134 + t125 * t316 + t127 * t317 + t278 * t129 + t274 * t320;
t135 = m(2) * t257 + t137;
t133 = t278 * t152 + (t146 * t284 + t147 * t281) * t274;
t131 = m(2) * t256 + t134;
t123 = mrSges(2,2) * t270 - mrSges(2,3) * t256 - t281 * t125 + t284 * t127 + (-t133 * t274 - t134 * t278) * pkin(7);
t122 = -mrSges(2,1) * t270 + mrSges(2,3) * t257 - pkin(1) * t133 - t274 * t129 + (t125 * t284 + t127 * t281 + t320) * t278;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t277 * t123 - t273 * t122 - qJ(1) * (t131 * t277 + t135 * t273), t123, t127, t130, t142, t155, t175; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t273 * t123 + t277 * t122 + qJ(1) * (-t131 * t273 + t135 * t277), t122, t125, t138, t140, t156, t174; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t294, t294, t129, t293, t300 * qJDD(2) - t322, -t290, t288;];
m_new  = t1;
