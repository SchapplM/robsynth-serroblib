% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:21:48
% EndTime: 2019-05-05 00:22:11
% DurationCPUTime: 17.44s
% Computational Cost: add. (353244->297), mult. (656940->382), div. (0->0), fcn. (467546->14), ass. (0->133)
t270 = sin(pkin(11));
t273 = cos(pkin(11));
t257 = t270 * g(1) - t273 * g(2);
t258 = -t273 * g(1) - t270 * g(2);
t268 = -g(3) + qJDD(1);
t278 = sin(qJ(2));
t274 = cos(pkin(6));
t282 = cos(qJ(2));
t301 = t274 * t282;
t271 = sin(pkin(6));
t303 = t271 * t282;
t218 = t257 * t301 - t278 * t258 + t268 * t303;
t216 = qJDD(2) * pkin(2) + t218;
t302 = t274 * t278;
t304 = t271 * t278;
t219 = t257 * t302 + t282 * t258 + t268 * t304;
t284 = qJD(2) ^ 2;
t217 = -t284 * pkin(2) + t219;
t269 = sin(pkin(12));
t272 = cos(pkin(12));
t205 = t269 * t216 + t272 * t217;
t202 = -t284 * pkin(3) + qJDD(2) * pkin(8) + t205;
t236 = -t271 * t257 + t274 * t268;
t235 = qJDD(3) + t236;
t277 = sin(qJ(4));
t281 = cos(qJ(4));
t192 = t281 * t202 + t277 * t235;
t253 = (-pkin(4) * t281 - pkin(9) * t277) * qJD(2);
t283 = qJD(4) ^ 2;
t299 = t281 * qJD(2);
t187 = -t283 * pkin(4) + qJDD(4) * pkin(9) + t253 * t299 + t192;
t204 = t272 * t216 - t269 * t217;
t201 = -qJDD(2) * pkin(3) - t284 * pkin(8) - t204;
t298 = qJD(2) * qJD(4);
t296 = t281 * t298;
t254 = t277 * qJDD(2) + t296;
t266 = t277 * t298;
t255 = t281 * qJDD(2) - t266;
t190 = (-t254 - t296) * pkin(9) + (-t255 + t266) * pkin(4) + t201;
t276 = sin(qJ(5));
t280 = cos(qJ(5));
t182 = -t276 * t187 + t280 * t190;
t300 = qJD(2) * t277;
t250 = t280 * qJD(4) - t276 * t300;
t226 = t250 * qJD(5) + t276 * qJDD(4) + t280 * t254;
t248 = qJDD(5) - t255;
t251 = t276 * qJD(4) + t280 * t300;
t264 = qJD(5) - t299;
t180 = (t250 * t264 - t226) * pkin(10) + (t250 * t251 + t248) * pkin(5) + t182;
t183 = t280 * t187 + t276 * t190;
t225 = -t251 * qJD(5) + t280 * qJDD(4) - t276 * t254;
t234 = t264 * pkin(5) - t251 * pkin(10);
t247 = t250 ^ 2;
t181 = -t247 * pkin(5) + t225 * pkin(10) - t264 * t234 + t183;
t275 = sin(qJ(6));
t279 = cos(qJ(6));
t179 = t275 * t180 + t279 * t181;
t191 = -t277 * t202 + t281 * t235;
t186 = -qJDD(4) * pkin(4) - t283 * pkin(9) + t253 * t300 - t191;
t184 = -t225 * pkin(5) - t247 * pkin(10) + t251 * t234 + t186;
t228 = t275 * t250 + t279 * t251;
t197 = -t228 * qJD(6) + t279 * t225 - t275 * t226;
t227 = t279 * t250 - t275 * t251;
t198 = t227 * qJD(6) + t275 * t225 + t279 * t226;
t263 = qJD(6) + t264;
t206 = Ifges(7,5) * t228 + Ifges(7,6) * t227 + Ifges(7,3) * t263;
t208 = Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t263;
t244 = qJDD(6) + t248;
t167 = -mrSges(7,1) * t184 + mrSges(7,3) * t179 + Ifges(7,4) * t198 + Ifges(7,2) * t197 + Ifges(7,6) * t244 - t228 * t206 + t263 * t208;
t178 = t279 * t180 - t275 * t181;
t207 = Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t263;
t168 = mrSges(7,2) * t184 - mrSges(7,3) * t178 + Ifges(7,1) * t198 + Ifges(7,4) * t197 + Ifges(7,5) * t244 + t227 * t206 - t263 * t207;
t220 = Ifges(6,5) * t251 + Ifges(6,6) * t250 + Ifges(6,3) * t264;
t222 = Ifges(6,1) * t251 + Ifges(6,4) * t250 + Ifges(6,5) * t264;
t214 = -t263 * mrSges(7,2) + t227 * mrSges(7,3);
t215 = t263 * mrSges(7,1) - t228 * mrSges(7,3);
t289 = m(7) * t184 - t197 * mrSges(7,1) + t198 * mrSges(7,2) - t227 * t214 + t228 * t215;
t210 = -t227 * mrSges(7,1) + t228 * mrSges(7,2);
t174 = m(7) * t178 + t244 * mrSges(7,1) - t198 * mrSges(7,3) - t228 * t210 + t263 * t214;
t175 = m(7) * t179 - t244 * mrSges(7,2) + t197 * mrSges(7,3) + t227 * t210 - t263 * t215;
t293 = -t275 * t174 + t279 * t175;
t153 = -mrSges(6,1) * t186 + mrSges(6,3) * t183 + Ifges(6,4) * t226 + Ifges(6,2) * t225 + Ifges(6,6) * t248 - pkin(5) * t289 + pkin(10) * t293 + t279 * t167 + t275 * t168 - t251 * t220 + t264 * t222;
t166 = t279 * t174 + t275 * t175;
t221 = Ifges(6,4) * t251 + Ifges(6,2) * t250 + Ifges(6,6) * t264;
t154 = mrSges(6,2) * t186 - mrSges(6,3) * t182 + Ifges(6,1) * t226 + Ifges(6,4) * t225 + Ifges(6,5) * t248 - pkin(10) * t166 - t275 * t167 + t279 * t168 + t250 * t220 - t264 * t221;
t229 = -t250 * mrSges(6,1) + t251 * mrSges(6,2);
t232 = -t264 * mrSges(6,2) + t250 * mrSges(6,3);
t164 = m(6) * t182 + t248 * mrSges(6,1) - t226 * mrSges(6,3) - t251 * t229 + t264 * t232 + t166;
t233 = t264 * mrSges(6,1) - t251 * mrSges(6,3);
t165 = m(6) * t183 - t248 * mrSges(6,2) + t225 * mrSges(6,3) + t250 * t229 - t264 * t233 + t293;
t162 = -t276 * t164 + t280 * t165;
t176 = -m(6) * t186 + t225 * mrSges(6,1) - t226 * mrSges(6,2) + t250 * t232 - t251 * t233 - t289;
t242 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t277 + Ifges(5,2) * t281) * qJD(2);
t243 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t277 + Ifges(5,4) * t281) * qJD(2);
t306 = mrSges(5,1) * t191 - mrSges(5,2) * t192 + Ifges(5,5) * t254 + Ifges(5,6) * t255 + Ifges(5,3) * qJDD(4) + pkin(4) * t176 + pkin(9) * t162 + t280 * t153 + t276 * t154 + (t277 * t242 - t281 * t243) * qJD(2);
t252 = (-mrSges(5,1) * t281 + mrSges(5,2) * t277) * qJD(2);
t259 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t300;
t159 = m(5) * t192 - qJDD(4) * mrSges(5,2) + t255 * mrSges(5,3) - qJD(4) * t259 + t252 * t299 + t162;
t260 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t299;
t170 = m(5) * t191 + qJDD(4) * mrSges(5,1) - t254 * mrSges(5,3) + qJD(4) * t260 - t252 * t300 + t176;
t294 = t281 * t159 - t277 * t170;
t149 = m(4) * t205 - t284 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t294;
t161 = t280 * t164 + t276 * t165;
t287 = -m(5) * t201 + t255 * mrSges(5,1) - t254 * mrSges(5,2) - t259 * t300 + t260 * t299 - t161;
t156 = m(4) * t204 + qJDD(2) * mrSges(4,1) - t284 * mrSges(4,2) + t287;
t144 = t269 * t149 + t272 * t156;
t142 = m(3) * t218 + qJDD(2) * mrSges(3,1) - t284 * mrSges(3,2) + t144;
t295 = t272 * t149 - t269 * t156;
t143 = m(3) * t219 - t284 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t295;
t136 = -t278 * t142 + t282 * t143;
t305 = pkin(7) * t136;
t152 = t277 * t159 + t281 * t170;
t297 = m(4) * t235 + t152;
t150 = m(3) * t236 + t297;
t132 = t142 * t301 + t143 * t302 - t271 * t150;
t241 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t277 + Ifges(5,6) * t281) * qJD(2);
t138 = mrSges(5,2) * t201 - mrSges(5,3) * t191 + Ifges(5,1) * t254 + Ifges(5,4) * t255 + Ifges(5,5) * qJDD(4) - pkin(9) * t161 - qJD(4) * t242 - t276 * t153 + t280 * t154 + t241 * t299;
t288 = -mrSges(7,1) * t178 + mrSges(7,2) * t179 - Ifges(7,5) * t198 - Ifges(7,6) * t197 - Ifges(7,3) * t244 - t228 * t207 + t227 * t208;
t285 = mrSges(6,1) * t182 - mrSges(6,2) * t183 + Ifges(6,5) * t226 + Ifges(6,6) * t225 + Ifges(6,3) * t248 + pkin(5) * t166 + t251 * t221 - t250 * t222 - t288;
t146 = -mrSges(5,1) * t201 + mrSges(5,3) * t192 + Ifges(5,4) * t254 + Ifges(5,2) * t255 + Ifges(5,6) * qJDD(4) - pkin(4) * t161 + qJD(4) * t243 - t241 * t300 - t285;
t128 = mrSges(4,2) * t235 - mrSges(4,3) * t204 + Ifges(4,5) * qJDD(2) - t284 * Ifges(4,6) - pkin(8) * t152 + t281 * t138 - t277 * t146;
t133 = -mrSges(4,1) * t235 + mrSges(4,3) * t205 + t284 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t152 - t306;
t123 = -mrSges(3,1) * t236 + mrSges(3,3) * t219 + t284 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t297 + qJ(3) * t295 + t269 * t128 + t272 * t133;
t125 = mrSges(3,2) * t236 - mrSges(3,3) * t218 + Ifges(3,5) * qJDD(2) - t284 * Ifges(3,6) - qJ(3) * t144 + t272 * t128 - t269 * t133;
t290 = mrSges(4,1) * t204 - mrSges(4,2) * t205 + Ifges(4,3) * qJDD(2) + pkin(3) * t287 + pkin(8) * t294 + t277 * t138 + t281 * t146;
t127 = mrSges(3,1) * t218 - mrSges(3,2) * t219 + Ifges(3,3) * qJDD(2) + pkin(2) * t144 + t290;
t291 = mrSges(2,1) * t257 - mrSges(2,2) * t258 + pkin(1) * t132 + t123 * t303 + t125 * t304 + t274 * t127 + t271 * t305;
t134 = m(2) * t258 + t136;
t131 = t274 * t150 + (t142 * t282 + t143 * t278) * t271;
t129 = m(2) * t257 + t132;
t121 = mrSges(2,2) * t268 - mrSges(2,3) * t257 - t278 * t123 + t282 * t125 + (-t131 * t271 - t132 * t274) * pkin(7);
t120 = -mrSges(2,1) * t268 + mrSges(2,3) * t258 - pkin(1) * t131 - t271 * t127 + (t123 * t282 + t125 * t278 + t305) * t274;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t273 * t121 - t270 * t120 - qJ(1) * (t273 * t129 + t270 * t134), t121, t125, t128, t138, t154, t168; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t270 * t121 + t273 * t120 + qJ(1) * (-t270 * t129 + t273 * t134), t120, t123, t133, t146, t153, t167; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t291, t291, t127, t290, t306, t285, -t288;];
m_new  = t1;
