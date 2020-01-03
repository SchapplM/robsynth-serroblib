% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:20
% EndTime: 2019-12-31 22:08:39
% DurationCPUTime: 8.99s
% Computational Cost: add. (151622->319), mult. (324725->408), div. (0->0), fcn. (247372->10), ass. (0->128)
t263 = sin(pkin(5));
t267 = sin(qJ(2));
t271 = cos(qJ(2));
t288 = qJD(1) * qJD(2);
t251 = (-qJDD(1) * t271 + t267 * t288) * t263;
t264 = cos(pkin(5));
t259 = t264 * qJD(1) + qJD(2);
t266 = sin(qJ(3));
t270 = cos(qJ(3));
t290 = qJD(1) * t263;
t285 = t267 * t290;
t239 = t270 * t259 - t266 * t285;
t250 = (qJDD(1) * t267 + t271 * t288) * t263;
t258 = t264 * qJDD(1) + qJDD(2);
t222 = t239 * qJD(3) + t270 * t250 + t266 * t258;
t240 = t266 * t259 + t270 * t285;
t289 = qJD(1) * t271;
t284 = t263 * t289;
t254 = qJD(3) - t284;
t265 = sin(qJ(4));
t269 = cos(qJ(4));
t228 = -t265 * t240 + t269 * t254;
t243 = qJDD(3) + t251;
t186 = t228 * qJD(4) + t269 * t222 + t265 * t243;
t229 = t269 * t240 + t265 * t254;
t203 = -t228 * mrSges(6,1) + t229 * mrSges(6,2);
t249 = (-pkin(2) * t271 - pkin(8) * t267) * t290;
t257 = t259 ^ 2;
t268 = sin(qJ(1));
t272 = cos(qJ(1));
t255 = t268 * g(1) - t272 * g(2);
t273 = qJD(1) ^ 2;
t300 = pkin(7) * t263;
t246 = qJDD(1) * pkin(1) + t273 * t300 + t255;
t256 = -t272 * g(1) - t268 * g(2);
t247 = -t273 * pkin(1) + qJDD(1) * t300 + t256;
t295 = t264 * t267;
t291 = t246 * t295 + t271 * t247;
t200 = -t257 * pkin(2) + t258 * pkin(8) + (-g(3) * t267 + t249 * t289) * t263 + t291;
t299 = t264 * g(3);
t201 = t251 * pkin(2) - t250 * pkin(8) - t299 + (-t246 + (pkin(2) * t267 - pkin(8) * t271) * t259 * qJD(1)) * t263;
t172 = t270 * t200 + t266 * t201;
t226 = -t239 * pkin(3) - t240 * pkin(9);
t253 = t254 ^ 2;
t167 = -t253 * pkin(3) + t243 * pkin(9) + t239 * t226 + t172;
t294 = t264 * t271;
t296 = t263 * t271;
t223 = -g(3) * t296 + t246 * t294 - t267 * t247;
t199 = -t258 * pkin(2) - t257 * pkin(8) + t249 * t285 - t223;
t221 = -t240 * qJD(3) - t266 * t250 + t270 * t258;
t170 = (-t239 * t254 - t222) * pkin(9) + (t240 * t254 - t221) * pkin(3) + t199;
t161 = -t265 * t167 + t269 * t170;
t219 = qJDD(4) - t221;
t238 = qJD(4) - t239;
t157 = -0.2e1 * qJD(5) * t229 + (t228 * t238 - t186) * qJ(5) + (t228 * t229 + t219) * pkin(4) + t161;
t206 = -t238 * mrSges(6,2) + t228 * mrSges(6,3);
t287 = m(6) * t157 + t219 * mrSges(6,1) + t238 * t206;
t154 = -t186 * mrSges(6,3) - t229 * t203 + t287;
t162 = t269 * t167 + t265 * t170;
t185 = -t229 * qJD(4) - t265 * t222 + t269 * t243;
t191 = Ifges(5,4) * t229 + Ifges(5,2) * t228 + Ifges(5,6) * t238;
t192 = Ifges(6,1) * t229 + Ifges(6,4) * t228 + Ifges(6,5) * t238;
t193 = Ifges(5,1) * t229 + Ifges(5,4) * t228 + Ifges(5,5) * t238;
t208 = t238 * pkin(4) - t229 * qJ(5);
t227 = t228 ^ 2;
t160 = -t227 * pkin(4) + t185 * qJ(5) + 0.2e1 * qJD(5) * t228 - t238 * t208 + t162;
t190 = Ifges(6,4) * t229 + Ifges(6,2) * t228 + Ifges(6,6) * t238;
t279 = -mrSges(6,1) * t157 + mrSges(6,2) * t160 - Ifges(6,5) * t186 - Ifges(6,6) * t185 - Ifges(6,3) * t219 - t229 * t190;
t301 = mrSges(5,1) * t161 - mrSges(5,2) * t162 + Ifges(5,5) * t186 + Ifges(5,6) * t185 + Ifges(5,3) * t219 + pkin(4) * t154 + t229 * t191 - (t193 + t192) * t228 - t279;
t298 = -mrSges(5,2) - mrSges(6,2);
t297 = t263 * t267;
t204 = -t228 * mrSges(5,1) + t229 * mrSges(5,2);
t207 = -t238 * mrSges(5,2) + t228 * mrSges(5,3);
t148 = m(5) * t161 + t219 * mrSges(5,1) + t238 * t207 + (-t203 - t204) * t229 + (-mrSges(5,3) - mrSges(6,3)) * t186 + t287;
t286 = m(6) * t160 + t185 * mrSges(6,3) + t228 * t203;
t209 = t238 * mrSges(6,1) - t229 * mrSges(6,3);
t292 = -t238 * mrSges(5,1) + t229 * mrSges(5,3) - t209;
t150 = m(5) * t162 + t185 * mrSges(5,3) + t228 * t204 + t298 * t219 + t292 * t238 + t286;
t147 = -t265 * t148 + t269 * t150;
t225 = -t239 * mrSges(4,1) + t240 * mrSges(4,2);
t231 = t254 * mrSges(4,1) - t240 * mrSges(4,3);
t144 = m(4) * t172 - t243 * mrSges(4,2) + t221 * mrSges(4,3) + t239 * t225 - t254 * t231 + t147;
t171 = -t266 * t200 + t270 * t201;
t166 = -t243 * pkin(3) - t253 * pkin(9) + t240 * t226 - t171;
t164 = -t185 * pkin(4) - t227 * qJ(5) + t229 * t208 + qJDD(5) + t166;
t282 = -m(6) * t164 + t185 * mrSges(6,1) + t228 * t206;
t153 = -m(5) * t166 + t185 * mrSges(5,1) + t298 * t186 + t228 * t207 + t292 * t229 + t282;
t230 = -t254 * mrSges(4,2) + t239 * mrSges(4,3);
t152 = m(4) * t171 + t243 * mrSges(4,1) - t222 * mrSges(4,3) - t240 * t225 + t254 * t230 + t153;
t138 = t266 * t144 + t270 * t152;
t224 = -g(3) * t297 + t291;
t244 = t259 * mrSges(3,1) - mrSges(3,3) * t285;
t248 = (-mrSges(3,1) * t271 + mrSges(3,2) * t267) * t290;
t283 = t270 * t144 - t266 * t152;
t136 = m(3) * t224 - t258 * mrSges(3,2) - t251 * mrSges(3,3) - t259 * t244 + t248 * t284 + t283;
t245 = -t259 * mrSges(3,2) + mrSges(3,3) * t284;
t146 = t269 * t148 + t265 * t150;
t276 = -m(4) * t199 + t221 * mrSges(4,1) - t222 * mrSges(4,2) + t239 * t230 - t240 * t231 - t146;
t141 = m(3) * t223 + t258 * mrSges(3,1) - t250 * mrSges(3,3) + t259 * t245 - t248 * t285 + t276;
t132 = t271 * t136 - t267 * t141;
t235 = -t263 * t246 - t299;
t137 = m(3) * t235 + t251 * mrSges(3,1) + t250 * mrSges(3,2) + (t244 * t267 - t245 * t271) * t290 + t138;
t128 = t136 * t295 - t263 * t137 + t141 * t294;
t280 = -mrSges(6,1) * t164 + mrSges(6,3) * t160 + Ifges(6,4) * t186 + Ifges(6,2) * t185 + Ifges(6,6) * t219 + t238 * t192;
t188 = Ifges(6,5) * t229 + Ifges(6,6) * t228 + Ifges(6,3) * t238;
t278 = mrSges(6,2) * t164 - mrSges(6,3) * t157 + Ifges(6,1) * t186 + Ifges(6,4) * t185 + Ifges(6,5) * t219 + t228 * t188;
t189 = Ifges(5,5) * t229 + Ifges(5,6) * t228 + Ifges(5,3) * t238;
t139 = Ifges(5,4) * t186 + Ifges(5,2) * t185 + Ifges(5,6) * t219 + t238 * t193 - mrSges(5,1) * t166 + mrSges(5,3) * t162 - pkin(4) * (t186 * mrSges(6,2) - t282) + qJ(5) * (-t219 * mrSges(6,2) - t238 * t209 + t286) + (-pkin(4) * t209 - t188 - t189) * t229 + t280;
t145 = mrSges(5,2) * t166 - mrSges(5,3) * t161 + Ifges(5,1) * t186 + Ifges(5,4) * t185 + Ifges(5,5) * t219 - qJ(5) * t154 + t228 * t189 + (-t190 - t191) * t238 + t278;
t215 = Ifges(4,5) * t240 + Ifges(4,6) * t239 + Ifges(4,3) * t254;
t216 = Ifges(4,4) * t240 + Ifges(4,2) * t239 + Ifges(4,6) * t254;
t129 = mrSges(4,2) * t199 - mrSges(4,3) * t171 + Ifges(4,1) * t222 + Ifges(4,4) * t221 + Ifges(4,5) * t243 - pkin(9) * t146 - t265 * t139 + t269 * t145 + t239 * t215 - t254 * t216;
t217 = Ifges(4,1) * t240 + Ifges(4,4) * t239 + Ifges(4,5) * t254;
t133 = -mrSges(4,1) * t199 + mrSges(4,3) * t172 + Ifges(4,4) * t222 + Ifges(4,2) * t221 + Ifges(4,6) * t243 - pkin(3) * t146 - t240 * t215 + t254 * t217 - t301;
t233 = Ifges(3,6) * t259 + (Ifges(3,4) * t267 + Ifges(3,2) * t271) * t290;
t234 = Ifges(3,5) * t259 + (Ifges(3,1) * t267 + Ifges(3,4) * t271) * t290;
t120 = Ifges(3,5) * t250 - Ifges(3,6) * t251 + Ifges(3,3) * t258 + mrSges(3,1) * t223 - mrSges(3,2) * t224 + t266 * t129 + t270 * t133 + pkin(2) * t276 + pkin(8) * t283 + (t233 * t267 - t234 * t271) * t290;
t232 = Ifges(3,3) * t259 + (Ifges(3,5) * t267 + Ifges(3,6) * t271) * t290;
t122 = mrSges(3,2) * t235 - mrSges(3,3) * t223 + Ifges(3,1) * t250 - Ifges(3,4) * t251 + Ifges(3,5) * t258 - pkin(8) * t138 + t270 * t129 - t266 * t133 + t232 * t284 - t259 * t233;
t274 = mrSges(4,1) * t171 - mrSges(4,2) * t172 + Ifges(4,5) * t222 + Ifges(4,6) * t221 + Ifges(4,3) * t243 + pkin(3) * t153 + pkin(9) * t147 + t269 * t139 + t265 * t145 + t240 * t216 - t239 * t217;
t124 = -mrSges(3,1) * t235 + mrSges(3,3) * t224 + Ifges(3,4) * t250 - Ifges(3,2) * t251 + Ifges(3,6) * t258 - pkin(2) * t138 - t232 * t285 + t259 * t234 - t274;
t277 = mrSges(2,1) * t255 - mrSges(2,2) * t256 + Ifges(2,3) * qJDD(1) + pkin(1) * t128 + t264 * t120 + t122 * t297 + t124 * t296 + t132 * t300;
t130 = m(2) * t256 - t273 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t132;
t127 = t264 * t137 + (t136 * t267 + t141 * t271) * t263;
t125 = m(2) * t255 + qJDD(1) * mrSges(2,1) - t273 * mrSges(2,2) + t128;
t118 = -mrSges(2,2) * g(3) - mrSges(2,3) * t255 + Ifges(2,5) * qJDD(1) - t273 * Ifges(2,6) + t271 * t122 - t267 * t124 + (-t127 * t263 - t128 * t264) * pkin(7);
t117 = mrSges(2,1) * g(3) + mrSges(2,3) * t256 + t273 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t127 - t263 * t120 + (pkin(7) * t132 + t122 * t267 + t124 * t271) * t264;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t272 * t118 - t268 * t117 - pkin(6) * (t272 * t125 + t268 * t130), t118, t122, t129, t145, -t238 * t190 + t278; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t268 * t118 + t272 * t117 + pkin(6) * (-t268 * t125 + t272 * t130), t117, t124, t133, t139, -t229 * t188 + t280; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t277, t277, t120, t274, t301, -t228 * t192 - t279;];
m_new = t1;
