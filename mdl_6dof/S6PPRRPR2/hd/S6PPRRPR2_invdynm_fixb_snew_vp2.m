% Calculate vector of cutting torques with Newton-Euler for
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:15:49
% EndTime: 2019-05-04 20:16:05
% DurationCPUTime: 11.54s
% Computational Cost: add. (218626->299), mult. (399140->374), div. (0->0), fcn. (287624->14), ass. (0->136)
t267 = sin(pkin(11));
t271 = cos(pkin(11));
t252 = -g(1) * t271 - g(2) * t267;
t266 = sin(pkin(12));
t270 = cos(pkin(12));
t251 = g(1) * t267 - g(2) * t271;
t265 = -g(3) + qJDD(1);
t269 = sin(pkin(6));
t273 = cos(pkin(6));
t297 = t251 * t273 + t265 * t269;
t208 = -t266 * t252 + t297 * t270;
t209 = t270 * t252 + t297 * t266;
t223 = -t251 * t269 + t265 * t273 + qJDD(2);
t279 = cos(qJ(3));
t272 = cos(pkin(7));
t276 = sin(qJ(3));
t309 = t272 * t276;
t268 = sin(pkin(7));
t310 = t268 * t276;
t201 = t208 * t309 + t279 * t209 + t223 * t310;
t281 = qJD(3) ^ 2;
t199 = -pkin(3) * t281 + qJDD(3) * pkin(9) + t201;
t275 = sin(qJ(4));
t196 = t275 * t199;
t203 = -t208 * t268 + t223 * t272;
t278 = cos(qJ(4));
t308 = t278 * t203;
t193 = -t196 + t308;
t245 = (mrSges(6,2) * t278 - mrSges(6,3) * t275) * qJD(3);
t246 = (-mrSges(5,1) * t278 + mrSges(5,2) * t275) * qJD(3);
t302 = qJD(3) * qJD(4);
t301 = t278 * t302;
t247 = qJDD(3) * t275 + t301;
t303 = qJD(3) * t278;
t254 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t303;
t255 = -mrSges(6,1) * t303 - qJD(4) * mrSges(6,3);
t244 = (-pkin(4) * t278 - qJ(5) * t275) * qJD(3);
t280 = qJD(4) ^ 2;
t304 = qJD(3) * t275;
t296 = -t280 * qJ(5) + t244 * t304 + qJDD(5) + t196;
t319 = pkin(10) * t281;
t320 = -pkin(4) - pkin(10);
t188 = t247 * pkin(5) + t320 * qJDD(4) + (-pkin(5) * t302 - t275 * t319 - t203) * t278 + t296;
t300 = t275 * t302;
t248 = qJDD(3) * t278 - t300;
t257 = pkin(5) * t304 - qJD(4) * pkin(10);
t264 = t278 ^ 2;
t200 = -t276 * t209 + (t208 * t272 + t223 * t268) * t279;
t290 = -qJDD(3) * pkin(3) - t200;
t321 = -2 * qJD(5);
t284 = pkin(4) * t300 + t304 * t321 + (-t247 - t301) * qJ(5) + t290;
t192 = -t257 * t304 + (-pkin(5) * t264 - pkin(9)) * t281 + t320 * t248 + t284;
t274 = sin(qJ(6));
t277 = cos(qJ(6));
t183 = t188 * t277 - t192 * t274;
t242 = -qJD(4) * t274 - t277 * t303;
t218 = qJD(6) * t242 + qJDD(4) * t277 - t248 * t274;
t243 = qJD(4) * t277 - t274 * t303;
t219 = -mrSges(7,1) * t242 + mrSges(7,2) * t243;
t259 = qJD(6) + t304;
t221 = -mrSges(7,2) * t259 + mrSges(7,3) * t242;
t241 = qJDD(6) + t247;
t180 = m(7) * t183 + mrSges(7,1) * t241 - mrSges(7,3) * t218 - t219 * t243 + t221 * t259;
t184 = t188 * t274 + t192 * t277;
t217 = -qJD(6) * t243 - qJDD(4) * t274 - t248 * t277;
t222 = mrSges(7,1) * t259 - mrSges(7,3) * t243;
t181 = m(7) * t184 - mrSges(7,2) * t241 + mrSges(7,3) * t217 + t219 * t242 - t222 * t259;
t169 = t277 * t180 + t274 * t181;
t191 = -qJDD(4) * pkin(4) + t296 - t308;
t293 = -m(6) * t191 - t247 * mrSges(6,1) - t169;
t166 = m(5) * t193 - t247 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t254 - t255) * qJD(4) + (-t245 - t246) * t304 + t293;
t194 = t278 * t199 + t275 * t203;
t253 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t304;
t289 = -t280 * pkin(4) + qJDD(4) * qJ(5) + t244 * t303 + t194;
t187 = -t264 * t319 + t248 * pkin(5) + ((2 * qJD(5)) + t257) * qJD(4) + t289;
t185 = -m(7) * t187 + t217 * mrSges(7,1) - t218 * mrSges(7,2) + t242 * t221 - t243 * t222;
t189 = qJD(4) * t321 - t289;
t256 = mrSges(6,1) * t304 + qJD(4) * mrSges(6,2);
t288 = -m(6) * t189 + qJDD(4) * mrSges(6,3) + qJD(4) * t256 + t245 * t303 - t185;
t173 = t246 * t303 + m(5) * t194 - qJDD(4) * mrSges(5,2) - qJD(4) * t253 + (mrSges(5,3) + mrSges(6,1)) * t248 + t288;
t299 = -t166 * t275 + t278 * t173;
t158 = m(4) * t201 - mrSges(4,1) * t281 - qJDD(3) * mrSges(4,2) + t299;
t161 = t278 * t166 + t275 * t173;
t160 = m(4) * t203 + t161;
t318 = t281 * pkin(9);
t198 = t290 - t318;
t170 = -t274 * t180 + t277 * t181;
t195 = -t248 * pkin(4) + t284 - t318;
t295 = -m(6) * t195 - t248 * mrSges(6,2) + t256 * t304 - t170;
t283 = -m(5) * t198 + t254 * t303 + t248 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t247 + (-t253 * t275 - t255 * t278) * qJD(3) + t295;
t164 = m(4) * t200 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,2) + t283;
t311 = t164 * t279;
t148 = t158 * t309 - t160 * t268 + t272 * t311;
t145 = m(3) * t208 + t148;
t154 = t279 * t158 - t164 * t276;
t153 = m(3) * t209 + t154;
t325 = t145 * t270 + t153 * t266;
t168 = -t247 * mrSges(6,3) + t255 * t303 - t295;
t229 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t275 + Ifges(5,4) * t278) * qJD(3);
t231 = Ifges(6,4) * qJD(4) + (-Ifges(6,2) * t275 - Ifges(6,6) * t278) * qJD(3);
t210 = Ifges(7,5) * t243 + Ifges(7,6) * t242 + Ifges(7,3) * t259;
t212 = Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t259;
t175 = -mrSges(7,1) * t187 + mrSges(7,3) * t184 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t241 - t210 * t243 + t212 * t259;
t211 = Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t259;
t176 = mrSges(7,2) * t187 - mrSges(7,3) * t183 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t241 + t210 * t242 - t211 * t259;
t286 = -mrSges(6,1) * t189 + mrSges(6,2) * t195 - pkin(5) * t185 - pkin(10) * t170 - t277 * t175 - t274 * t176;
t232 = Ifges(6,1) * qJD(4) + (-Ifges(6,4) * t275 - Ifges(6,5) * t278) * qJD(3);
t306 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t275 + Ifges(5,6) * t278) * qJD(3) + t232;
t317 = Ifges(5,4) + Ifges(6,6);
t149 = -mrSges(5,1) * t198 + mrSges(5,3) * t194 - pkin(4) * t168 + (Ifges(5,2) + Ifges(6,3)) * t248 + t317 * t247 + (Ifges(5,6) - Ifges(6,5)) * qJDD(4) + (t229 - t231) * qJD(4) - t306 * t304 + t286;
t291 = mrSges(7,1) * t183 - mrSges(7,2) * t184 + Ifges(7,5) * t218 + Ifges(7,6) * t217 + Ifges(7,3) * t241 + t243 * t211 - t242 * t212;
t285 = mrSges(6,1) * t191 - mrSges(6,3) * t195 + pkin(5) * t169 + t291;
t230 = Ifges(6,5) * qJD(4) + (-Ifges(6,6) * t275 - Ifges(6,3) * t278) * qJD(3);
t305 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t275 + Ifges(5,2) * t278) * qJD(3) - t230;
t150 = t317 * t248 + (Ifges(5,1) + Ifges(6,2)) * t247 + (Ifges(5,5) - Ifges(6,4)) * qJDD(4) - t305 * qJD(4) + t306 * t303 + mrSges(5,2) * t198 - mrSges(5,3) * t193 - qJ(5) * t168 + t285;
t138 = mrSges(4,1) * t200 - mrSges(4,2) * t201 + Ifges(4,3) * qJDD(3) + pkin(3) * t283 + pkin(9) * t299 + t278 * t149 + t275 * t150;
t147 = t158 * t310 + t272 * t160 + t268 * t311;
t139 = mrSges(4,2) * t203 - mrSges(4,3) * t200 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t281 - pkin(9) * t161 - t149 * t275 + t150 * t278;
t287 = -mrSges(6,2) * t191 + mrSges(6,3) * t189 - Ifges(6,1) * qJDD(4) + Ifges(6,4) * t247 + Ifges(6,5) * t248 + pkin(10) * t169 + t274 * t175 - t277 * t176 - t231 * t303;
t323 = (-t278 * t229 + t305 * t275) * qJD(3) + mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t247 + Ifges(5,6) * t248 + Ifges(5,3) * qJDD(4) + pkin(4) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t255 - t245 * t304 + t293) + qJ(5) * (mrSges(6,1) * t248 + t288) - t287;
t143 = -mrSges(4,1) * t203 + mrSges(4,3) * t201 + t281 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t161 - t323;
t294 = pkin(8) * t154 + t139 * t276 + t143 * t279;
t131 = -mrSges(3,1) * t223 + mrSges(3,3) * t209 - pkin(2) * t147 - t268 * t138 + t294 * t272;
t133 = mrSges(3,2) * t223 - mrSges(3,3) * t208 + t279 * t139 - t276 * t143 + (-t147 * t268 - t148 * t272) * pkin(8);
t142 = -t145 * t266 + t270 * t153;
t324 = qJ(2) * t142 + t131 * t270 + t133 * t266;
t146 = m(3) * t223 + t147;
t137 = -t146 * t269 + t325 * t273;
t129 = mrSges(3,1) * t208 - mrSges(3,2) * t209 + pkin(2) * t148 + t272 * t138 + t294 * t268;
t292 = mrSges(2,1) * t251 - mrSges(2,2) * t252 + pkin(1) * t137 + t273 * t129 + t324 * t269;
t140 = m(2) * t252 + t142;
t136 = t273 * t146 + t325 * t269;
t134 = m(2) * t251 + t137;
t127 = mrSges(2,2) * t265 - mrSges(2,3) * t251 - t266 * t131 + t270 * t133 + (-t136 * t269 - t137 * t273) * qJ(2);
t126 = -mrSges(2,1) * t265 + mrSges(2,3) * t252 - pkin(1) * t136 - t269 * t129 + t324 * t273;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t271 * t127 - t267 * t126 - qJ(1) * (t134 * t271 + t140 * t267), t127, t133, t139, t150, -t230 * t304 - t287, t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t267 * t127 + t271 * t126 + qJ(1) * (-t134 * t267 + t140 * t271), t126, t131, t143, t149, Ifges(6,4) * qJDD(4) - Ifges(6,2) * t247 - Ifges(6,6) * t248 - qJD(4) * t230 - t232 * t303 - t285, t175; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t292, t292, t129, t138, t323, Ifges(6,5) * qJDD(4) - Ifges(6,6) * t247 - Ifges(6,3) * t248 + qJD(4) * t231 + t232 * t304 - t286, t291;];
m_new  = t1;
