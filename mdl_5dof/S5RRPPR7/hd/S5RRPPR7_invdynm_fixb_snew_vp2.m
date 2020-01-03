% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:35:02
% EndTime: 2019-12-31 19:35:12
% DurationCPUTime: 4.76s
% Computational Cost: add. (48825->315), mult. (112383->379), div. (0->0), fcn. (71241->8), ass. (0->120)
t265 = sin(qJ(2));
t268 = cos(qJ(2));
t291 = qJD(1) * qJD(2);
t248 = t265 * qJDD(1) + t268 * t291;
t266 = sin(qJ(1));
t269 = cos(qJ(1));
t255 = -t269 * g(1) - t266 * g(2);
t271 = qJD(1) ^ 2;
t243 = -t271 * pkin(1) + qJDD(1) * pkin(6) + t255;
t299 = t265 * t243;
t302 = pkin(2) * t271;
t187 = qJDD(2) * pkin(2) - t248 * qJ(3) - t299 + (qJ(3) * t291 + t265 * t302 - g(3)) * t268;
t223 = -g(3) * t265 + t268 * t243;
t249 = t268 * qJDD(1) - t265 * t291;
t295 = qJD(1) * t265;
t251 = qJD(2) * pkin(2) - qJ(3) * t295;
t262 = t268 ^ 2;
t188 = t249 * qJ(3) - qJD(2) * t251 - t262 * t302 + t223;
t263 = sin(pkin(8));
t300 = cos(pkin(8));
t237 = (t263 * t268 + t300 * t265) * qJD(1);
t306 = -2 * qJD(3);
t170 = t300 * t187 - t263 * t188 + t237 * t306;
t294 = qJD(1) * t268;
t236 = t263 * t295 - t300 * t294;
t231 = t236 * t306;
t298 = t263 * t187 + t300 * t188;
t171 = t231 + t298;
t198 = Ifges(4,4) * t237 - Ifges(4,2) * t236 + Ifges(4,6) * qJD(2);
t206 = -mrSges(5,2) * t236 - mrSges(5,3) * t237;
t217 = t248 * t263 - t300 * t249;
t218 = t300 * t248 + t263 * t249;
t226 = mrSges(5,1) * t236 - qJD(2) * mrSges(5,3);
t204 = pkin(3) * t236 - qJ(4) * t237;
t270 = qJD(2) ^ 2;
t166 = -qJDD(2) * pkin(3) - t270 * qJ(4) + t237 * t204 + qJDD(4) - t170;
t293 = qJD(2) * t236;
t160 = (t236 * t237 - qJDD(2)) * pkin(7) + (t218 + t293) * pkin(4) + t166;
t228 = pkin(4) * t237 - qJD(2) * pkin(7);
t235 = t236 ^ 2;
t254 = t266 * g(1) - t269 * g(2);
t285 = -qJDD(1) * pkin(1) - t254;
t192 = -t249 * pkin(2) + qJDD(3) + t251 * t295 + (-qJ(3) * t262 - pkin(6)) * t271 + t285;
t303 = -2 * qJD(4);
t274 = (-t218 + t293) * qJ(4) + t192 + (qJD(2) * pkin(3) + t303) * t237;
t163 = -t235 * pkin(4) - t237 * t228 + (pkin(3) + pkin(7)) * t217 + t274;
t264 = sin(qJ(5));
t267 = cos(qJ(5));
t158 = t267 * t160 - t264 * t163;
t219 = -t264 * qJD(2) + t267 * t236;
t181 = t219 * qJD(5) + t267 * qJDD(2) + t264 * t217;
t220 = t267 * qJD(2) + t264 * t236;
t189 = -mrSges(6,1) * t219 + mrSges(6,2) * t220;
t234 = qJD(5) + t237;
t193 = -mrSges(6,2) * t234 + mrSges(6,3) * t219;
t216 = qJDD(5) + t218;
t154 = m(6) * t158 + mrSges(6,1) * t216 - mrSges(6,3) * t181 - t189 * t220 + t193 * t234;
t159 = t264 * t160 + t267 * t163;
t180 = -t220 * qJD(5) - t264 * qJDD(2) + t267 * t217;
t194 = mrSges(6,1) * t234 - mrSges(6,3) * t220;
t155 = m(6) * t159 - mrSges(6,2) * t216 + mrSges(6,3) * t180 + t189 * t219 - t194 * t234;
t142 = t267 * t154 + t264 * t155;
t284 = t270 * pkin(3) - qJDD(2) * qJ(4) - t298;
t162 = -pkin(4) * t217 - pkin(7) * t235 - t204 * t236 + t231 + ((2 * qJD(4)) + t228) * qJD(2) - t284;
t173 = Ifges(6,5) * t220 + Ifges(6,6) * t219 + Ifges(6,3) * t234;
t175 = Ifges(6,1) * t220 + Ifges(6,4) * t219 + Ifges(6,5) * t234;
t145 = -mrSges(6,1) * t162 + mrSges(6,3) * t159 + Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t216 - t173 * t220 + t175 * t234;
t174 = Ifges(6,4) * t220 + Ifges(6,2) * t219 + Ifges(6,6) * t234;
t146 = mrSges(6,2) * t162 - mrSges(6,3) * t158 + Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t216 + t173 * t219 - t174 * t234;
t164 = qJD(2) * t303 + ((2 * qJD(3)) + t204) * t236 + t284;
t195 = Ifges(5,5) * qJD(2) - Ifges(5,6) * t237 + Ifges(5,3) * t236;
t279 = -mrSges(5,2) * t166 + mrSges(5,3) * t164 - Ifges(5,1) * qJDD(2) + Ifges(5,4) * t218 - Ifges(5,5) * t217 + pkin(7) * t142 + t264 * t145 - t267 * t146 + t237 * t195;
t156 = -m(6) * t162 + mrSges(6,1) * t180 - t181 * mrSges(6,2) + t193 * t219 - t220 * t194;
t227 = mrSges(5,1) * t237 + qJD(2) * mrSges(5,2);
t280 = -m(5) * t164 + qJDD(2) * mrSges(5,3) + qJD(2) * t227 - t156;
t282 = -m(5) * t166 - t218 * mrSges(5,1) - t237 * t206 - t142;
t197 = Ifges(5,4) * qJD(2) - Ifges(5,2) * t237 + Ifges(5,6) * t236;
t296 = Ifges(4,1) * t237 - Ifges(4,4) * t236 + Ifges(4,5) * qJD(2) - t197;
t307 = t296 * t236 - mrSges(4,2) * t171 + pkin(3) * (-qJDD(2) * mrSges(5,2) - qJD(2) * t226 + t282) + qJ(4) * (-mrSges(5,1) * t217 - t206 * t236 + t280) + mrSges(4,1) * t170 + t237 * t198 - Ifges(4,6) * t217 + Ifges(4,5) * t218 + Ifges(4,3) * qJDD(2) - t279;
t205 = mrSges(4,1) * t236 + mrSges(4,2) * t237;
t224 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t236;
t138 = m(4) * t170 - t218 * mrSges(4,3) - t237 * t205 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t224 - t226) * qJD(2) + t282;
t225 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t237;
t149 = m(4) * t171 - qJDD(2) * mrSges(4,2) - qJD(2) * t225 + (-t205 - t206) * t236 + (-mrSges(4,3) - mrSges(5,1)) * t217 + t280;
t134 = t300 * t138 + t263 * t149;
t222 = -t268 * g(3) - t299;
t239 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t265 + Ifges(3,2) * t268) * qJD(1);
t240 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t265 + Ifges(3,4) * t268) * qJD(1);
t305 = mrSges(3,1) * t222 - mrSges(3,2) * t223 + Ifges(3,5) * t248 + Ifges(3,6) * t249 + Ifges(3,3) * qJDD(2) + pkin(2) * t134 + (t265 * t239 - t268 * t240) * qJD(1) + t307;
t301 = Ifges(4,4) + Ifges(5,6);
t143 = -t264 * t154 + t267 * t155;
t199 = Ifges(5,1) * qJD(2) - Ifges(5,4) * t237 + Ifges(5,5) * t236;
t297 = -Ifges(4,5) * t237 + Ifges(4,6) * t236 - Ifges(4,3) * qJD(2) - t199;
t247 = (-mrSges(3,1) * t268 + mrSges(3,2) * t265) * qJD(1);
t253 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t294;
t132 = m(3) * t222 + qJDD(2) * mrSges(3,1) - t248 * mrSges(3,3) + qJD(2) * t253 - t247 * t295 + t134;
t252 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t295;
t288 = -t263 * t138 + t300 * t149;
t133 = m(3) * t223 - qJDD(2) * mrSges(3,2) + t249 * mrSges(3,3) - qJD(2) * t252 + t247 * t294 + t288;
t289 = -t132 * t265 + t268 * t133;
t168 = t217 * pkin(3) + t274;
t139 = m(5) * t168 - t217 * mrSges(5,2) - t218 * mrSges(5,3) - t236 * t226 - t237 * t227 + t143;
t278 = -mrSges(5,1) * t164 + mrSges(5,2) * t168 - pkin(4) * t156 - pkin(7) * t143 - t267 * t145 - t264 * t146;
t126 = -mrSges(4,1) * t192 + mrSges(4,3) * t171 - pkin(3) * t139 + t297 * t237 + t301 * t218 + (-Ifges(4,2) - Ifges(5,3)) * t217 + (Ifges(4,6) - Ifges(5,5)) * qJDD(2) + t296 * qJD(2) + t278;
t281 = mrSges(6,1) * t158 - mrSges(6,2) * t159 + Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t216 + t220 * t174 - t219 * t175;
t276 = mrSges(5,1) * t166 - mrSges(5,3) * t168 + pkin(4) * t142 + t281;
t130 = t297 * t236 + (Ifges(4,1) + Ifges(5,2)) * t218 - t301 * t217 + (Ifges(4,5) - Ifges(5,4)) * qJDD(2) + (-t198 + t195) * qJD(2) + mrSges(4,2) * t192 - mrSges(4,3) * t170 - qJ(4) * t139 + t276;
t238 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t265 + Ifges(3,6) * t268) * qJD(1);
t242 = -t271 * pkin(6) + t285;
t277 = m(4) * t192 + t217 * mrSges(4,1) + t218 * mrSges(4,2) + t236 * t224 + t237 * t225 + t139;
t122 = -mrSges(3,1) * t242 + mrSges(3,3) * t223 + Ifges(3,4) * t248 + Ifges(3,2) * t249 + Ifges(3,6) * qJDD(2) - pkin(2) * t277 + qJ(3) * t288 + qJD(2) * t240 + t300 * t126 + t263 * t130 - t238 * t295;
t125 = mrSges(3,2) * t242 - mrSges(3,3) * t222 + Ifges(3,1) * t248 + Ifges(3,4) * t249 + Ifges(3,5) * qJDD(2) - qJ(3) * t134 - qJD(2) * t239 - t263 * t126 + t300 * t130 + t238 * t294;
t273 = -m(3) * t242 + t249 * mrSges(3,1) - t248 * mrSges(3,2) - t252 * t295 + t253 * t294 - t277;
t283 = mrSges(2,1) * t254 - mrSges(2,2) * t255 + Ifges(2,3) * qJDD(1) + pkin(1) * t273 + pkin(6) * t289 + t268 * t122 + t265 * t125;
t135 = m(2) * t254 + qJDD(1) * mrSges(2,1) - t271 * mrSges(2,2) + t273;
t129 = t268 * t132 + t265 * t133;
t127 = m(2) * t255 - t271 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t289;
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t255 + t271 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t129 - t305;
t120 = -mrSges(2,2) * g(3) - mrSges(2,3) * t254 + Ifges(2,5) * qJDD(1) - t271 * Ifges(2,6) - pkin(6) * t129 - t265 * t122 + t268 * t125;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t120 - t266 * t123 - pkin(5) * (t266 * t127 + t269 * t135), t120, t125, t130, -t236 * t197 - t279, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t266 * t120 + t269 * t123 + pkin(5) * (t269 * t127 - t266 * t135), t123, t122, t126, Ifges(5,4) * qJDD(2) - Ifges(5,2) * t218 + Ifges(5,6) * t217 - qJD(2) * t195 + t236 * t199 - t276, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t283, t283, t305, t307, Ifges(5,5) * qJDD(2) - Ifges(5,6) * t218 + Ifges(5,3) * t217 + qJD(2) * t197 + t237 * t199 - t278, t281;];
m_new = t1;
