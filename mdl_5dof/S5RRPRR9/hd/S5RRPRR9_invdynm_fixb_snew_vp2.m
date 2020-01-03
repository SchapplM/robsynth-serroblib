% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:12
% EndTime: 2019-12-31 20:20:28
% DurationCPUTime: 8.96s
% Computational Cost: add. (144878->311), mult. (329232->397), div. (0->0), fcn. (227761->10), ass. (0->124)
t266 = sin(qJ(2));
t270 = cos(qJ(2));
t290 = qJD(1) * qJD(2);
t249 = qJDD(1) * t266 + t270 * t290;
t267 = sin(qJ(1));
t271 = cos(qJ(1));
t256 = -g(1) * t271 - g(2) * t267;
t273 = qJD(1) ^ 2;
t244 = -pkin(1) * t273 + qJDD(1) * pkin(6) + t256;
t294 = t266 * t244;
t295 = pkin(2) * t273;
t203 = qJDD(2) * pkin(2) - t249 * qJ(3) - t294 + (qJ(3) * t290 + t266 * t295 - g(3)) * t270;
t229 = -g(3) * t266 + t270 * t244;
t250 = qJDD(1) * t270 - t266 * t290;
t293 = qJD(1) * t266;
t252 = qJD(2) * pkin(2) - qJ(3) * t293;
t261 = t270 ^ 2;
t204 = qJ(3) * t250 - qJD(2) * t252 - t261 * t295 + t229;
t262 = sin(pkin(9));
t263 = cos(pkin(9));
t238 = (t262 * t270 + t263 * t266) * qJD(1);
t185 = -0.2e1 * qJD(3) * t238 + t203 * t263 - t262 * t204;
t292 = qJD(1) * t270;
t237 = -t262 * t293 + t263 * t292;
t186 = 0.2e1 * qJD(3) * t237 + t262 * t203 + t263 * t204;
t215 = -mrSges(4,1) * t237 + mrSges(4,2) * t238;
t223 = -t262 * t249 + t250 * t263;
t231 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t238;
t217 = -pkin(3) * t237 - pkin(7) * t238;
t272 = qJD(2) ^ 2;
t176 = -pkin(3) * t272 + qJDD(2) * pkin(7) + t217 * t237 + t186;
t255 = t267 * g(1) - t271 * g(2);
t283 = -qJDD(1) * pkin(1) - t255;
t207 = -t250 * pkin(2) + qJDD(3) + t252 * t293 + (-qJ(3) * t261 - pkin(6)) * t273 + t283;
t224 = t249 * t263 + t250 * t262;
t179 = (-qJD(2) * t237 - t224) * pkin(7) + (qJD(2) * t238 - t223) * pkin(3) + t207;
t265 = sin(qJ(4));
t269 = cos(qJ(4));
t166 = -t265 * t176 + t179 * t269;
t226 = qJD(2) * t269 - t238 * t265;
t196 = qJD(4) * t226 + qJDD(2) * t265 + t224 * t269;
t222 = qJDD(4) - t223;
t227 = qJD(2) * t265 + t238 * t269;
t236 = qJD(4) - t237;
t163 = (t226 * t236 - t196) * pkin(8) + (t226 * t227 + t222) * pkin(4) + t166;
t167 = t176 * t269 + t179 * t265;
t195 = -qJD(4) * t227 + qJDD(2) * t269 - t224 * t265;
t210 = pkin(4) * t236 - pkin(8) * t227;
t225 = t226 ^ 2;
t164 = -pkin(4) * t225 + pkin(8) * t195 - t210 * t236 + t167;
t264 = sin(qJ(5));
t268 = cos(qJ(5));
t161 = t163 * t268 - t164 * t264;
t200 = t226 * t268 - t227 * t264;
t172 = qJD(5) * t200 + t195 * t264 + t196 * t268;
t201 = t226 * t264 + t227 * t268;
t187 = -mrSges(6,1) * t200 + mrSges(6,2) * t201;
t232 = qJD(5) + t236;
t188 = -mrSges(6,2) * t232 + mrSges(6,3) * t200;
t218 = qJDD(5) + t222;
t156 = m(6) * t161 + mrSges(6,1) * t218 - mrSges(6,3) * t172 - t187 * t201 + t188 * t232;
t162 = t163 * t264 + t164 * t268;
t171 = -qJD(5) * t201 + t195 * t268 - t196 * t264;
t189 = mrSges(6,1) * t232 - mrSges(6,3) * t201;
t157 = m(6) * t162 - mrSges(6,2) * t218 + mrSges(6,3) * t171 + t187 * t200 - t189 * t232;
t148 = t156 * t268 + t157 * t264;
t205 = -mrSges(5,1) * t226 + mrSges(5,2) * t227;
t208 = -mrSges(5,2) * t236 + mrSges(5,3) * t226;
t146 = m(5) * t166 + mrSges(5,1) * t222 - mrSges(5,3) * t196 - t205 * t227 + t208 * t236 + t148;
t209 = mrSges(5,1) * t236 - mrSges(5,3) * t227;
t286 = -t156 * t264 + t157 * t268;
t147 = m(5) * t167 - mrSges(5,2) * t222 + mrSges(5,3) * t195 + t205 * t226 - t209 * t236 + t286;
t287 = -t146 * t265 + t147 * t269;
t139 = m(4) * t186 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t223 - qJD(2) * t231 + t215 * t237 + t287;
t230 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t237;
t175 = -qJDD(2) * pkin(3) - pkin(7) * t272 + t238 * t217 - t185;
t165 = -pkin(4) * t195 - pkin(8) * t225 + t210 * t227 + t175;
t281 = m(6) * t165 - mrSges(6,1) * t171 + mrSges(6,2) * t172 - t188 * t200 + t189 * t201;
t277 = -m(5) * t175 + mrSges(5,1) * t195 - mrSges(5,2) * t196 + t208 * t226 - t209 * t227 - t281;
t152 = m(4) * t185 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t224 + qJD(2) * t230 - t215 * t238 + t277;
t132 = t139 * t262 + t152 * t263;
t228 = -t270 * g(3) - t294;
t240 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t266 + Ifges(3,2) * t270) * qJD(1);
t241 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t266 + Ifges(3,4) * t270) * qJD(1);
t180 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t232;
t182 = Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t232;
t149 = -mrSges(6,1) * t165 + mrSges(6,3) * t162 + Ifges(6,4) * t172 + Ifges(6,2) * t171 + Ifges(6,6) * t218 - t180 * t201 + t182 * t232;
t181 = Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t232;
t150 = mrSges(6,2) * t165 - mrSges(6,3) * t161 + Ifges(6,1) * t172 + Ifges(6,4) * t171 + Ifges(6,5) * t218 + t180 * t200 - t181 * t232;
t190 = Ifges(5,5) * t227 + Ifges(5,6) * t226 + Ifges(5,3) * t236;
t192 = Ifges(5,1) * t227 + Ifges(5,4) * t226 + Ifges(5,5) * t236;
t131 = -mrSges(5,1) * t175 + mrSges(5,3) * t167 + Ifges(5,4) * t196 + Ifges(5,2) * t195 + Ifges(5,6) * t222 - pkin(4) * t281 + pkin(8) * t286 + t268 * t149 + t264 * t150 - t227 * t190 + t236 * t192;
t191 = Ifges(5,4) * t227 + Ifges(5,2) * t226 + Ifges(5,6) * t236;
t134 = mrSges(5,2) * t175 - mrSges(5,3) * t166 + Ifges(5,1) * t196 + Ifges(5,4) * t195 + Ifges(5,5) * t222 - pkin(8) * t148 - t149 * t264 + t150 * t268 + t190 * t226 - t191 * t236;
t212 = Ifges(4,4) * t238 + Ifges(4,2) * t237 + Ifges(4,6) * qJD(2);
t213 = Ifges(4,1) * t238 + Ifges(4,4) * t237 + Ifges(4,5) * qJD(2);
t278 = -mrSges(4,1) * t185 + mrSges(4,2) * t186 - Ifges(4,5) * t224 - Ifges(4,6) * t223 - Ifges(4,3) * qJDD(2) - pkin(3) * t277 - pkin(7) * t287 - t131 * t269 - t134 * t265 - t238 * t212 + t237 * t213;
t296 = mrSges(3,1) * t228 - mrSges(3,2) * t229 + Ifges(3,5) * t249 + Ifges(3,6) * t250 + Ifges(3,3) * qJDD(2) + pkin(2) * t132 + (t240 * t266 - t241 * t270) * qJD(1) - t278;
t141 = t146 * t269 + t147 * t265;
t248 = (-mrSges(3,1) * t270 + mrSges(3,2) * t266) * qJD(1);
t254 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t292;
t128 = m(3) * t228 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t249 + qJD(2) * t254 - t248 * t293 + t132;
t253 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t293;
t288 = t139 * t263 - t152 * t262;
t129 = m(3) * t229 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t250 - qJD(2) * t253 + t248 * t292 + t288;
t289 = -t128 * t266 + t129 * t270;
t211 = Ifges(4,5) * t238 + Ifges(4,6) * t237 + Ifges(4,3) * qJD(2);
t122 = mrSges(4,2) * t207 - mrSges(4,3) * t185 + Ifges(4,1) * t224 + Ifges(4,4) * t223 + Ifges(4,5) * qJDD(2) - pkin(7) * t141 - qJD(2) * t212 - t131 * t265 + t134 * t269 + t211 * t237;
t280 = -mrSges(6,1) * t161 + mrSges(6,2) * t162 - Ifges(6,5) * t172 - Ifges(6,6) * t171 - Ifges(6,3) * t218 - t181 * t201 + t182 * t200;
t274 = mrSges(5,1) * t166 - mrSges(5,2) * t167 + Ifges(5,5) * t196 + Ifges(5,6) * t195 + Ifges(5,3) * t222 + pkin(4) * t148 + t191 * t227 - t192 * t226 - t280;
t126 = -mrSges(4,1) * t207 + mrSges(4,3) * t186 + Ifges(4,4) * t224 + Ifges(4,2) * t223 + Ifges(4,6) * qJDD(2) - pkin(3) * t141 + qJD(2) * t213 - t211 * t238 - t274;
t239 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t266 + Ifges(3,6) * t270) * qJD(1);
t243 = -t273 * pkin(6) + t283;
t279 = m(4) * t207 - t223 * mrSges(4,1) + mrSges(4,2) * t224 - t237 * t230 + t231 * t238 + t141;
t118 = -mrSges(3,1) * t243 + mrSges(3,3) * t229 + Ifges(3,4) * t249 + Ifges(3,2) * t250 + Ifges(3,6) * qJDD(2) - pkin(2) * t279 + qJ(3) * t288 + qJD(2) * t241 + t262 * t122 + t263 * t126 - t239 * t293;
t121 = mrSges(3,2) * t243 - mrSges(3,3) * t228 + Ifges(3,1) * t249 + Ifges(3,4) * t250 + Ifges(3,5) * qJDD(2) - qJ(3) * t132 - qJD(2) * t240 + t122 * t263 - t126 * t262 + t239 * t292;
t276 = -m(3) * t243 + t250 * mrSges(3,1) - mrSges(3,2) * t249 - t253 * t293 + t254 * t292 - t279;
t282 = mrSges(2,1) * t255 - mrSges(2,2) * t256 + Ifges(2,3) * qJDD(1) + pkin(1) * t276 + pkin(6) * t289 + t118 * t270 + t121 * t266;
t135 = m(2) * t255 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t273 + t276;
t125 = t128 * t270 + t129 * t266;
t123 = m(2) * t256 - mrSges(2,1) * t273 - qJDD(1) * mrSges(2,2) + t289;
t119 = mrSges(2,1) * g(3) + mrSges(2,3) * t256 + t273 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t125 - t296;
t116 = -mrSges(2,2) * g(3) - mrSges(2,3) * t255 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t273 - pkin(6) * t125 - t118 * t266 + t121 * t270;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t271 * t116 - t267 * t119 - pkin(5) * (t123 * t267 + t135 * t271), t116, t121, t122, t134, t150; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t267 * t116 + t271 * t119 + pkin(5) * (t123 * t271 - t135 * t267), t119, t118, t126, t131, t149; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t282, t282, t296, -t278, t274, -t280;];
m_new = t1;
