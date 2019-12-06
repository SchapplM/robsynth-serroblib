% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:32
% EndTime: 2019-12-05 18:11:51
% DurationCPUTime: 12.56s
% Computational Cost: add. (192158->291), mult. (477246->365), div. (0->0), fcn. (362937->10), ass. (0->125)
t271 = qJD(1) ^ 2;
t266 = sin(qJ(1));
t270 = cos(qJ(1));
t244 = -t270 * g(1) - t266 * g(2);
t237 = -t271 * pkin(1) + qJDD(1) * qJ(2) + t244;
t261 = sin(pkin(9));
t262 = cos(pkin(9));
t295 = qJD(1) * qJD(2);
t293 = -t262 * g(3) - 0.2e1 * t261 * t295;
t300 = pkin(2) * t262;
t210 = (-pkin(6) * qJDD(1) + t271 * t300 - t237) * t261 + t293;
t228 = -t261 * g(3) + (t237 + 0.2e1 * t295) * t262;
t294 = qJDD(1) * t262;
t257 = t262 ^ 2;
t298 = t257 * t271;
t211 = -pkin(2) * t298 + pkin(6) * t294 + t228;
t265 = sin(qJ(3));
t269 = cos(qJ(3));
t191 = t269 * t210 - t265 * t211;
t282 = t261 * t269 + t262 * t265;
t281 = -t261 * t265 + t262 * t269;
t235 = t281 * qJD(1);
t296 = t235 * qJD(3);
t226 = qJDD(1) * t282 + t296;
t236 = t282 * qJD(1);
t172 = (-t226 + t296) * pkin(7) + (t235 * t236 + qJDD(3)) * pkin(3) + t191;
t192 = t265 * t210 + t269 * t211;
t225 = -t236 * qJD(3) + qJDD(1) * t281;
t231 = qJD(3) * pkin(3) - t236 * pkin(7);
t234 = t235 ^ 2;
t180 = -t234 * pkin(3) + t225 * pkin(7) - qJD(3) * t231 + t192;
t264 = sin(qJ(4));
t268 = cos(qJ(4));
t161 = t268 * t172 - t264 * t180;
t216 = t268 * t235 - t264 * t236;
t189 = t216 * qJD(4) + t264 * t225 + t268 * t226;
t217 = t264 * t235 + t268 * t236;
t255 = qJDD(3) + qJDD(4);
t258 = qJD(3) + qJD(4);
t156 = (t216 * t258 - t189) * pkin(8) + (t216 * t217 + t255) * pkin(4) + t161;
t162 = t264 * t172 + t268 * t180;
t188 = -t217 * qJD(4) + t268 * t225 - t264 * t226;
t208 = t258 * pkin(4) - t217 * pkin(8);
t212 = t216 ^ 2;
t157 = -t212 * pkin(4) + t188 * pkin(8) - t258 * t208 + t162;
t263 = sin(qJ(5));
t267 = cos(qJ(5));
t154 = t267 * t156 - t263 * t157;
t200 = t267 * t216 - t263 * t217;
t168 = t200 * qJD(5) + t263 * t188 + t267 * t189;
t201 = t263 * t216 + t267 * t217;
t178 = -t200 * mrSges(6,1) + t201 * mrSges(6,2);
t253 = qJD(5) + t258;
t193 = -t253 * mrSges(6,2) + t200 * mrSges(6,3);
t252 = qJDD(5) + t255;
t151 = m(6) * t154 + t252 * mrSges(6,1) - t168 * mrSges(6,3) - t201 * t178 + t253 * t193;
t155 = t263 * t156 + t267 * t157;
t167 = -t201 * qJD(5) + t267 * t188 - t263 * t189;
t194 = t253 * mrSges(6,1) - t201 * mrSges(6,3);
t152 = m(6) * t155 - t252 * mrSges(6,2) + t167 * mrSges(6,3) + t200 * t178 - t253 * t194;
t143 = t267 * t151 + t263 * t152;
t202 = -t216 * mrSges(5,1) + t217 * mrSges(5,2);
t206 = -t258 * mrSges(5,2) + t216 * mrSges(5,3);
t140 = m(5) * t161 + t255 * mrSges(5,1) - t189 * mrSges(5,3) - t217 * t202 + t258 * t206 + t143;
t207 = t258 * mrSges(5,1) - t217 * mrSges(5,3);
t289 = -t263 * t151 + t267 * t152;
t141 = m(5) * t162 - t255 * mrSges(5,2) + t188 * mrSges(5,3) + t216 * t202 - t258 * t207 + t289;
t136 = t268 * t140 + t264 * t141;
t220 = -t235 * mrSges(4,1) + t236 * mrSges(4,2);
t229 = -qJD(3) * mrSges(4,2) + t235 * mrSges(4,3);
t133 = m(4) * t191 + qJDD(3) * mrSges(4,1) - t226 * mrSges(4,3) + qJD(3) * t229 - t236 * t220 + t136;
t230 = qJD(3) * mrSges(4,1) - t236 * mrSges(4,3);
t290 = -t264 * t140 + t268 * t141;
t134 = m(4) * t192 - qJDD(3) * mrSges(4,2) + t225 * mrSges(4,3) - qJD(3) * t230 + t235 * t220 + t290;
t127 = t269 * t133 + t265 * t134;
t227 = -t261 * t237 + t293;
t214 = Ifges(4,4) * t236 + Ifges(4,2) * t235 + Ifges(4,6) * qJD(3);
t215 = Ifges(4,1) * t236 + Ifges(4,4) * t235 + Ifges(4,5) * qJD(3);
t196 = Ifges(5,4) * t217 + Ifges(5,2) * t216 + Ifges(5,6) * t258;
t197 = Ifges(5,1) * t217 + Ifges(5,4) * t216 + Ifges(5,5) * t258;
t174 = Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t253;
t175 = Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t253;
t277 = -mrSges(6,1) * t154 + mrSges(6,2) * t155 - Ifges(6,5) * t168 - Ifges(6,6) * t167 - Ifges(6,3) * t252 - t201 * t174 + t200 * t175;
t275 = -mrSges(5,1) * t161 + mrSges(5,2) * t162 - Ifges(5,5) * t189 - Ifges(5,6) * t188 - Ifges(5,3) * t255 - pkin(4) * t143 - t217 * t196 + t216 * t197 + t277;
t273 = -mrSges(4,1) * t191 + mrSges(4,2) * t192 - Ifges(4,5) * t226 - Ifges(4,6) * t225 - Ifges(4,3) * qJDD(3) - pkin(3) * t136 - t236 * t214 + t235 * t215 + t275;
t286 = Ifges(3,4) * t261 + Ifges(3,2) * t262;
t287 = Ifges(3,1) * t261 + Ifges(3,4) * t262;
t301 = -mrSges(3,1) * t227 + mrSges(3,2) * t228 - pkin(2) * t127 - (t261 * t286 - t262 * t287) * t271 + t273;
t299 = mrSges(3,2) * t261;
t285 = Ifges(3,5) * t261 + Ifges(3,6) * t262;
t297 = t271 * t285;
t243 = t266 * g(1) - t270 * g(2);
t280 = mrSges(3,3) * qJDD(1) + t271 * (-mrSges(3,1) * t262 + t299);
t125 = m(3) * t227 - t261 * t280 + t127;
t291 = -t265 * t133 + t269 * t134;
t126 = m(3) * t228 + t262 * t280 + t291;
t292 = -t261 * t125 + t262 * t126;
t288 = qJDD(2) - t243;
t256 = t261 ^ 2;
t224 = (-pkin(1) - t300) * qJDD(1) + (-qJ(2) + (-t256 - t257) * pkin(6)) * t271 + t288;
t183 = -t225 * pkin(3) - t234 * pkin(7) + t236 * t231 + t224;
t159 = -t188 * pkin(4) - t212 * pkin(8) + t217 * t208 + t183;
t284 = m(6) * t159 - t167 * mrSges(6,1) + t168 * mrSges(6,2) - t200 * t193 + t201 * t194;
t173 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t253;
t144 = -mrSges(6,1) * t159 + mrSges(6,3) * t155 + Ifges(6,4) * t168 + Ifges(6,2) * t167 + Ifges(6,6) * t252 - t201 * t173 + t253 * t175;
t145 = mrSges(6,2) * t159 - mrSges(6,3) * t154 + Ifges(6,1) * t168 + Ifges(6,4) * t167 + Ifges(6,5) * t252 + t200 * t173 - t253 * t174;
t195 = Ifges(5,5) * t217 + Ifges(5,6) * t216 + Ifges(5,3) * t258;
t128 = -mrSges(5,1) * t183 + mrSges(5,3) * t162 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * t255 - pkin(4) * t284 + pkin(8) * t289 + t267 * t144 + t263 * t145 - t217 * t195 + t258 * t197;
t129 = mrSges(5,2) * t183 - mrSges(5,3) * t161 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * t255 - pkin(8) * t143 - t263 * t144 + t267 * t145 + t216 * t195 - t258 * t196;
t213 = Ifges(4,5) * t236 + Ifges(4,6) * t235 + Ifges(4,3) * qJD(3);
t278 = m(5) * t183 - t188 * mrSges(5,1) + t189 * mrSges(5,2) - t216 * t206 + t217 * t207 + t284;
t122 = -mrSges(4,1) * t224 + mrSges(4,3) * t192 + Ifges(4,4) * t226 + Ifges(4,2) * t225 + Ifges(4,6) * qJDD(3) - pkin(3) * t278 + pkin(7) * t290 + qJD(3) * t215 + t268 * t128 + t264 * t129 - t236 * t213;
t123 = mrSges(4,2) * t224 - mrSges(4,3) * t191 + Ifges(4,1) * t226 + Ifges(4,4) * t225 + Ifges(4,5) * qJDD(3) - pkin(7) * t136 - qJD(3) * t214 - t264 * t128 + t268 * t129 + t235 * t213;
t233 = -qJDD(1) * pkin(1) - t271 * qJ(2) + t288;
t276 = m(4) * t224 - t225 * mrSges(4,1) + t226 * mrSges(4,2) - t235 * t229 + t236 * t230 + t278;
t115 = -mrSges(3,1) * t233 + mrSges(3,3) * t228 - pkin(2) * t276 + pkin(6) * t291 + qJDD(1) * t286 + t269 * t122 + t265 * t123 - t261 * t297;
t117 = mrSges(3,2) * t233 - mrSges(3,3) * t227 - pkin(6) * t127 + qJDD(1) * t287 - t265 * t122 + t269 * t123 + t262 * t297;
t274 = -m(3) * t233 + mrSges(3,1) * t294 - t276 + (t256 * t271 + t298) * mrSges(3,3);
t279 = -mrSges(2,2) * t244 + qJ(2) * t292 + t262 * t115 + t261 * t117 + pkin(1) * (-qJDD(1) * t299 + t274) + mrSges(2,1) * t243 + Ifges(2,3) * qJDD(1);
t146 = t274 + (mrSges(2,1) - t299) * qJDD(1) - t271 * mrSges(2,2) + m(2) * t243;
t121 = t262 * t125 + t261 * t126;
t119 = m(2) * t244 - t271 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t292;
t118 = mrSges(2,1) * g(3) - pkin(1) * t121 + (Ifges(2,6) - t285) * qJDD(1) + t271 * Ifges(2,5) + mrSges(2,3) * t244 + t301;
t113 = -mrSges(2,2) * g(3) - mrSges(2,3) * t243 + Ifges(2,5) * qJDD(1) - t271 * Ifges(2,6) - qJ(2) * t121 - t261 * t115 + t262 * t117;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t270 * t113 - t266 * t118 - pkin(5) * (t266 * t119 + t270 * t146), t113, t117, t123, t129, t145; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t266 * t113 + t270 * t118 + pkin(5) * (t270 * t119 - t266 * t146), t118, t115, t122, t128, t144; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t279, t279, qJDD(1) * t285 - t301, -t273, -t275, -t277;];
m_new = t1;
