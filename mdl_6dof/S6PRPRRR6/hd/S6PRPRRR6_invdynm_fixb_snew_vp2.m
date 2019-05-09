% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:30:48
% EndTime: 2019-05-05 01:31:05
% DurationCPUTime: 9.40s
% Computational Cost: add. (181276->296), mult. (339984->369), div. (0->0), fcn. (226702->12), ass. (0->131)
t277 = sin(pkin(11));
t279 = cos(pkin(11));
t259 = g(1) * t277 - g(2) * t279;
t260 = -g(1) * t279 - g(2) * t277;
t274 = -g(3) + qJDD(1);
t278 = sin(pkin(6));
t280 = cos(pkin(6));
t284 = sin(qJ(2));
t288 = cos(qJ(2));
t215 = -t284 * t260 + (t259 * t280 + t274 * t278) * t288;
t290 = qJD(2) ^ 2;
t297 = -qJ(3) * t290 + qJDD(3) - t215;
t319 = -pkin(2) - pkin(8);
t209 = t319 * qJDD(2) + t297;
t233 = -t259 * t278 + t274 * t280;
t283 = sin(qJ(4));
t287 = cos(qJ(4));
t203 = t283 * t209 + t287 * t233;
t254 = (mrSges(5,1) * t283 + mrSges(5,2) * t287) * qJD(2);
t309 = qJD(2) * qJD(4);
t266 = t287 * t309;
t256 = -t283 * qJDD(2) - t266;
t310 = qJD(2) * t287;
t262 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t310;
t268 = t283 * qJD(2);
t255 = (pkin(4) * t283 - pkin(9) * t287) * qJD(2);
t289 = qJD(4) ^ 2;
t194 = -pkin(4) * t289 + qJDD(4) * pkin(9) - t255 * t268 + t203;
t311 = t280 * t284;
t312 = t278 * t284;
t216 = t259 * t311 + t288 * t260 + t274 * t312;
t320 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t216;
t208 = t319 * t290 - t320;
t307 = t283 * t309;
t257 = qJDD(2) * t287 - t307;
t197 = (-t257 + t307) * pkin(9) + (-t256 + t266) * pkin(4) + t208;
t282 = sin(qJ(5));
t286 = cos(qJ(5));
t183 = -t194 * t282 + t286 * t197;
t252 = qJD(4) * t286 - t282 * t310;
t223 = qJD(5) * t252 + qJDD(4) * t282 + t257 * t286;
t249 = qJDD(5) - t256;
t253 = qJD(4) * t282 + t286 * t310;
t265 = t268 + qJD(5);
t181 = (t252 * t265 - t223) * pkin(10) + (t252 * t253 + t249) * pkin(5) + t183;
t184 = t286 * t194 + t282 * t197;
t222 = -qJD(5) * t253 + qJDD(4) * t286 - t257 * t282;
t231 = pkin(5) * t265 - pkin(10) * t253;
t248 = t252 ^ 2;
t182 = -pkin(5) * t248 + pkin(10) * t222 - t231 * t265 + t184;
t281 = sin(qJ(6));
t285 = cos(qJ(6));
t179 = t181 * t285 - t182 * t281;
t224 = t252 * t285 - t253 * t281;
t191 = qJD(6) * t224 + t222 * t281 + t223 * t285;
t225 = t252 * t281 + t253 * t285;
t205 = -mrSges(7,1) * t224 + mrSges(7,2) * t225;
t264 = qJD(6) + t265;
t213 = -mrSges(7,2) * t264 + mrSges(7,3) * t224;
t242 = qJDD(6) + t249;
t174 = m(7) * t179 + mrSges(7,1) * t242 - t191 * mrSges(7,3) - t205 * t225 + t213 * t264;
t180 = t181 * t281 + t182 * t285;
t190 = -qJD(6) * t225 + t222 * t285 - t223 * t281;
t214 = mrSges(7,1) * t264 - mrSges(7,3) * t225;
t175 = m(7) * t180 - mrSges(7,2) * t242 + t190 * mrSges(7,3) + t205 * t224 - t214 * t264;
t167 = t285 * t174 + t281 * t175;
t226 = -mrSges(6,1) * t252 + mrSges(6,2) * t253;
t229 = -mrSges(6,2) * t265 + mrSges(6,3) * t252;
t165 = m(6) * t183 + mrSges(6,1) * t249 - mrSges(6,3) * t223 - t226 * t253 + t229 * t265 + t167;
t230 = mrSges(6,1) * t265 - mrSges(6,3) * t253;
t305 = -t174 * t281 + t285 * t175;
t166 = m(6) * t184 - mrSges(6,2) * t249 + mrSges(6,3) * t222 + t226 * t252 - t230 * t265 + t305;
t306 = -t165 * t282 + t286 * t166;
t157 = m(5) * t203 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t256 - qJD(4) * t262 - t254 * t268 + t306;
t202 = t209 * t287 - t283 * t233;
t261 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t268;
t193 = -qJDD(4) * pkin(4) - pkin(9) * t289 + t255 * t310 - t202;
t185 = -pkin(5) * t222 - pkin(10) * t248 + t231 * t253 + t193;
t300 = m(7) * t185 - t190 * mrSges(7,1) + t191 * mrSges(7,2) - t224 * t213 + t214 * t225;
t292 = -m(6) * t193 + t222 * mrSges(6,1) - mrSges(6,2) * t223 + t252 * t229 - t230 * t253 - t300;
t170 = m(5) * t202 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t257 + qJD(4) * t261 - t254 * t310 + t292;
t150 = t287 * t157 - t170 * t283;
t147 = m(4) * t233 + t150;
t199 = Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t264;
t201 = Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t264;
t168 = -mrSges(7,1) * t185 + mrSges(7,3) * t180 + Ifges(7,4) * t191 + Ifges(7,2) * t190 + Ifges(7,6) * t242 - t199 * t225 + t201 * t264;
t200 = Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t264;
t169 = mrSges(7,2) * t185 - mrSges(7,3) * t179 + Ifges(7,1) * t191 + Ifges(7,4) * t190 + Ifges(7,5) * t242 + t199 * t224 - t200 * t264;
t217 = Ifges(6,5) * t253 + Ifges(6,6) * t252 + Ifges(6,3) * t265;
t219 = Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t265;
t148 = -mrSges(6,1) * t193 + mrSges(6,3) * t184 + Ifges(6,4) * t223 + Ifges(6,2) * t222 + Ifges(6,6) * t249 - pkin(5) * t300 + pkin(10) * t305 + t285 * t168 + t281 * t169 - t253 * t217 + t265 * t219;
t218 = Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t265;
t152 = mrSges(6,2) * t193 - mrSges(6,3) * t183 + Ifges(6,1) * t223 + Ifges(6,4) * t222 + Ifges(6,5) * t249 - pkin(10) * t167 - t168 * t281 + t169 * t285 + t217 * t252 - t218 * t265;
t160 = t286 * t165 + t282 * t166;
t239 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t287 - Ifges(5,6) * t283) * qJD(2);
t240 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t287 - Ifges(5,2) * t283) * qJD(2);
t138 = mrSges(5,2) * t208 - mrSges(5,3) * t202 + Ifges(5,1) * t257 + Ifges(5,4) * t256 + Ifges(5,5) * qJDD(4) - pkin(9) * t160 - qJD(4) * t240 - t148 * t282 + t152 * t286 - t239 * t268;
t241 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t287 - Ifges(5,4) * t283) * qJD(2);
t299 = -mrSges(7,1) * t179 + mrSges(7,2) * t180 - Ifges(7,5) * t191 - Ifges(7,6) * t190 - Ifges(7,3) * t242 - t225 * t200 + t224 * t201;
t291 = mrSges(6,1) * t183 - mrSges(6,2) * t184 + Ifges(6,5) * t223 + Ifges(6,6) * t222 + Ifges(6,3) * t249 + pkin(5) * t167 + t253 * t218 - t252 * t219 - t299;
t142 = -mrSges(5,1) * t208 + mrSges(5,3) * t203 + Ifges(5,4) * t257 + Ifges(5,2) * t256 + Ifges(5,6) * qJDD(4) - pkin(4) * t160 + qJD(4) * t241 - t239 * t310 - t291;
t158 = -m(5) * t208 + t256 * mrSges(5,1) - t257 * mrSges(5,2) - t261 * t268 - t262 * t310 - t160;
t210 = pkin(2) * t290 + t320;
t296 = -mrSges(4,1) * t210 - pkin(3) * t158 - pkin(8) * t150 - t138 * t283 - t142 * t287;
t315 = Ifges(4,5) - Ifges(3,6);
t316 = -Ifges(4,4) + Ifges(3,5);
t317 = mrSges(3,1) - mrSges(4,2);
t130 = mrSges(3,3) * t216 - pkin(2) * t147 - t315 * qJDD(2) - t317 * t233 + t316 * t290 + t296;
t149 = t157 * t283 + t170 * t287;
t212 = -qJDD(2) * pkin(2) + t297;
t302 = -m(4) * t212 + t290 * mrSges(4,3) - t149;
t144 = m(3) * t215 - mrSges(3,2) * t290 + t317 * qJDD(2) + t302;
t294 = -m(4) * t210 + t290 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t158;
t155 = m(3) * t216 - mrSges(3,1) * t290 - qJDD(2) * mrSges(3,2) + t294;
t141 = -t144 * t284 + t288 * t155;
t321 = pkin(7) * t141 + t130 * t288;
t313 = t144 * t288;
t145 = m(3) * t233 + t147;
t136 = -t145 * t278 + t155 * t311 + t280 * t313;
t298 = mrSges(4,2) * t212 - mrSges(4,3) * t210 + Ifges(4,1) * qJDD(2) - pkin(8) * t149 + t287 * t138 - t283 * t142;
t128 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t215 - mrSges(3,2) * t216 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t302) + qJ(3) * t294 + t298;
t295 = mrSges(5,1) * t202 - mrSges(5,2) * t203 + Ifges(5,5) * t257 + Ifges(5,6) * t256 + Ifges(5,3) * qJDD(4) + pkin(4) * t292 + pkin(9) * t306 + t286 * t148 + t282 * t152 + t240 * t310 + t241 * t268;
t293 = mrSges(4,1) * t212 + pkin(3) * t149 + t295;
t132 = t315 * t290 + t316 * qJDD(2) + (mrSges(3,2) - mrSges(4,3)) * t233 - mrSges(3,3) * t215 + t293 - qJ(3) * t147;
t301 = mrSges(2,1) * t259 - mrSges(2,2) * t260 + pkin(1) * t136 + t280 * t128 + t132 * t312 + t278 * t321;
t139 = m(2) * t260 + t141;
t135 = t145 * t280 + (t155 * t284 + t313) * t278;
t133 = m(2) * t259 + t136;
t126 = mrSges(2,2) * t274 - mrSges(2,3) * t259 - t130 * t284 + t132 * t288 + (-t135 * t278 - t136 * t280) * pkin(7);
t125 = -mrSges(2,1) * t274 + mrSges(2,3) * t260 - pkin(1) * t135 - t278 * t128 + (t132 * t284 + t321) * t280;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t279 * t126 - t277 * t125 - qJ(1) * (t133 * t279 + t139 * t277), t126, t132, t298, t138, t152, t169; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t277 * t126 + t279 * t125 + qJ(1) * (-t277 * t133 + t139 * t279), t125, t130, mrSges(4,3) * t233 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t290 - t293, t142, t148, t168; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t301, t301, t128, -mrSges(4,2) * t233 + Ifges(4,4) * t290 + Ifges(4,5) * qJDD(2) - t296, t295, t291, -t299;];
m_new  = t1;
