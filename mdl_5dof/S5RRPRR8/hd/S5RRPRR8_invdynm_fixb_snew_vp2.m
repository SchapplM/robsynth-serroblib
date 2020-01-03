% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:17:03
% EndTime: 2019-12-31 20:17:19
% DurationCPUTime: 9.80s
% Computational Cost: add. (154909->311), mult. (359241->398), div. (0->0), fcn. (255297->10), ass. (0->124)
t264 = sin(qJ(2));
t268 = cos(qJ(2));
t286 = qJD(1) * qJD(2);
t244 = qJDD(1) * t264 + t268 * t286;
t265 = sin(qJ(1));
t269 = cos(qJ(1));
t251 = -g(1) * t269 - g(2) * t265;
t270 = qJD(1) ^ 2;
t239 = -pkin(1) * t270 + qJDD(1) * pkin(6) + t251;
t289 = t264 * t239;
t290 = pkin(2) * t270;
t203 = qJDD(2) * pkin(2) - t244 * qJ(3) - t289 + (qJ(3) * t286 + t264 * t290 - g(3)) * t268;
t225 = -g(3) * t264 + t268 * t239;
t245 = qJDD(1) * t268 - t264 * t286;
t288 = qJD(1) * t264;
t247 = qJD(2) * pkin(2) - qJ(3) * t288;
t259 = t268 ^ 2;
t204 = qJ(3) * t245 - qJD(2) * t247 - t259 * t290 + t225;
t260 = sin(pkin(9));
t261 = cos(pkin(9));
t234 = (t260 * t268 + t261 * t264) * qJD(1);
t179 = -0.2e1 * qJD(3) * t234 + t261 * t203 - t260 * t204;
t223 = t244 * t261 + t245 * t260;
t233 = (-t260 * t264 + t261 * t268) * qJD(1);
t166 = (qJD(2) * t233 - t223) * pkin(7) + (t233 * t234 + qJDD(2)) * pkin(3) + t179;
t180 = 0.2e1 * qJD(3) * t233 + t260 * t203 + t261 * t204;
t222 = -t244 * t260 + t245 * t261;
t228 = qJD(2) * pkin(3) - pkin(7) * t234;
t232 = t233 ^ 2;
t168 = -pkin(3) * t232 + pkin(7) * t222 - qJD(2) * t228 + t180;
t263 = sin(qJ(4));
t267 = cos(qJ(4));
t164 = t263 * t166 + t267 * t168;
t215 = t233 * t263 + t234 * t267;
t187 = -qJD(4) * t215 + t222 * t267 - t223 * t263;
t214 = t233 * t267 - t234 * t263;
t197 = -mrSges(5,1) * t214 + mrSges(5,2) * t215;
t256 = qJD(2) + qJD(4);
t209 = mrSges(5,1) * t256 - mrSges(5,3) * t215;
t255 = qJDD(2) + qJDD(4);
t198 = -pkin(4) * t214 - pkin(8) * t215;
t254 = t256 ^ 2;
t160 = -pkin(4) * t254 + pkin(8) * t255 + t198 * t214 + t164;
t250 = t265 * g(1) - t269 * g(2);
t280 = -qJDD(1) * pkin(1) - t250;
t205 = -t245 * pkin(2) + qJDD(3) + t247 * t288 + (-qJ(3) * t259 - pkin(6)) * t270 + t280;
t177 = -t222 * pkin(3) - t232 * pkin(7) + t234 * t228 + t205;
t188 = qJD(4) * t214 + t222 * t263 + t223 * t267;
t161 = (-t214 * t256 - t188) * pkin(8) + (t215 * t256 - t187) * pkin(4) + t177;
t262 = sin(qJ(5));
t266 = cos(qJ(5));
t157 = -t160 * t262 + t161 * t266;
t206 = -t215 * t262 + t256 * t266;
t171 = qJD(5) * t206 + t188 * t266 + t255 * t262;
t186 = qJDD(5) - t187;
t207 = t215 * t266 + t256 * t262;
t189 = -mrSges(6,1) * t206 + mrSges(6,2) * t207;
t210 = qJD(5) - t214;
t190 = -mrSges(6,2) * t210 + mrSges(6,3) * t206;
t153 = m(6) * t157 + mrSges(6,1) * t186 - t171 * mrSges(6,3) - t189 * t207 + t190 * t210;
t158 = t160 * t266 + t161 * t262;
t170 = -qJD(5) * t207 - t188 * t262 + t255 * t266;
t191 = mrSges(6,1) * t210 - mrSges(6,3) * t207;
t154 = m(6) * t158 - mrSges(6,2) * t186 + t170 * mrSges(6,3) + t189 * t206 - t191 * t210;
t282 = -t153 * t262 + t266 * t154;
t140 = m(5) * t164 - mrSges(5,2) * t255 + mrSges(5,3) * t187 + t197 * t214 - t209 * t256 + t282;
t163 = t166 * t267 - t168 * t263;
t208 = -mrSges(5,2) * t256 + mrSges(5,3) * t214;
t159 = -pkin(4) * t255 - pkin(8) * t254 + t198 * t215 - t163;
t277 = -m(6) * t159 + t170 * mrSges(6,1) - t171 * mrSges(6,2) + t206 * t190 - t191 * t207;
t149 = m(5) * t163 + mrSges(5,1) * t255 - mrSges(5,3) * t188 - t197 * t215 + t208 * t256 + t277;
t135 = t263 * t140 + t267 * t149;
t218 = -mrSges(4,1) * t233 + mrSges(4,2) * t234;
t226 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t233;
t132 = m(4) * t179 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t223 + qJD(2) * t226 - t218 * t234 + t135;
t227 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t234;
t283 = t267 * t140 - t149 * t263;
t133 = m(4) * t180 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t222 - qJD(2) * t227 + t218 * t233 + t283;
t126 = t261 * t132 + t260 * t133;
t224 = -t268 * g(3) - t289;
t236 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t264 + Ifges(3,2) * t268) * qJD(1);
t237 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t264 + Ifges(3,4) * t268) * qJD(1);
t212 = Ifges(4,4) * t234 + Ifges(4,2) * t233 + Ifges(4,6) * qJD(2);
t213 = Ifges(4,1) * t234 + Ifges(4,4) * t233 + Ifges(4,5) * qJD(2);
t172 = Ifges(6,5) * t207 + Ifges(6,6) * t206 + Ifges(6,3) * t210;
t174 = Ifges(6,1) * t207 + Ifges(6,4) * t206 + Ifges(6,5) * t210;
t146 = -mrSges(6,1) * t159 + mrSges(6,3) * t158 + Ifges(6,4) * t171 + Ifges(6,2) * t170 + Ifges(6,6) * t186 - t172 * t207 + t174 * t210;
t173 = Ifges(6,4) * t207 + Ifges(6,2) * t206 + Ifges(6,6) * t210;
t147 = mrSges(6,2) * t159 - mrSges(6,3) * t157 + Ifges(6,1) * t171 + Ifges(6,4) * t170 + Ifges(6,5) * t186 + t172 * t206 - t173 * t210;
t193 = Ifges(5,4) * t215 + Ifges(5,2) * t214 + Ifges(5,6) * t256;
t194 = Ifges(5,1) * t215 + Ifges(5,4) * t214 + Ifges(5,5) * t256;
t276 = -mrSges(5,1) * t163 + mrSges(5,2) * t164 - Ifges(5,5) * t188 - Ifges(5,6) * t187 - Ifges(5,3) * t255 - pkin(4) * t277 - pkin(8) * t282 - t266 * t146 - t262 * t147 - t215 * t193 + t214 * t194;
t273 = -mrSges(4,1) * t179 + mrSges(4,2) * t180 - Ifges(4,5) * t223 - Ifges(4,6) * t222 - Ifges(4,3) * qJDD(2) - pkin(3) * t135 - t234 * t212 + t233 * t213 + t276;
t291 = mrSges(3,1) * t224 - mrSges(3,2) * t225 + Ifges(3,5) * t244 + Ifges(3,6) * t245 + Ifges(3,3) * qJDD(2) + pkin(2) * t126 + (t236 * t264 - t237 * t268) * qJD(1) - t273;
t142 = t266 * t153 + t262 * t154;
t287 = qJD(1) * t268;
t243 = (-mrSges(3,1) * t268 + mrSges(3,2) * t264) * qJD(1);
t249 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t287;
t124 = m(3) * t224 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t244 + qJD(2) * t249 - t243 * t288 + t126;
t248 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t288;
t284 = -t132 * t260 + t261 * t133;
t125 = m(3) * t225 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t245 - qJD(2) * t248 + t243 * t287 + t284;
t285 = -t124 * t264 + t268 * t125;
t279 = m(5) * t177 - t187 * mrSges(5,1) + t188 * mrSges(5,2) - t214 * t208 + t215 * t209 + t142;
t192 = Ifges(5,5) * t215 + Ifges(5,6) * t214 + Ifges(5,3) * t256;
t127 = mrSges(5,2) * t177 - mrSges(5,3) * t163 + Ifges(5,1) * t188 + Ifges(5,4) * t187 + Ifges(5,5) * t255 - pkin(8) * t142 - t146 * t262 + t147 * t266 + t192 * t214 - t193 * t256;
t274 = mrSges(6,1) * t157 - mrSges(6,2) * t158 + Ifges(6,5) * t171 + Ifges(6,6) * t170 + Ifges(6,3) * t186 + t173 * t207 - t174 * t206;
t128 = -mrSges(5,1) * t177 + mrSges(5,3) * t164 + Ifges(5,4) * t188 + Ifges(5,2) * t187 + Ifges(5,6) * t255 - pkin(4) * t142 - t192 * t215 + t194 * t256 - t274;
t211 = Ifges(4,5) * t234 + Ifges(4,6) * t233 + Ifges(4,3) * qJD(2);
t118 = -mrSges(4,1) * t205 + mrSges(4,3) * t180 + Ifges(4,4) * t223 + Ifges(4,2) * t222 + Ifges(4,6) * qJDD(2) - pkin(3) * t279 + pkin(7) * t283 + qJD(2) * t213 + t263 * t127 + t267 * t128 - t234 * t211;
t122 = mrSges(4,2) * t205 - mrSges(4,3) * t179 + Ifges(4,1) * t223 + Ifges(4,4) * t222 + Ifges(4,5) * qJDD(2) - pkin(7) * t135 - qJD(2) * t212 + t127 * t267 - t128 * t263 + t211 * t233;
t235 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t264 + Ifges(3,6) * t268) * qJD(1);
t238 = -t270 * pkin(6) + t280;
t275 = m(4) * t205 - t222 * mrSges(4,1) + mrSges(4,2) * t223 - t233 * t226 + t227 * t234 + t279;
t114 = -mrSges(3,1) * t238 + mrSges(3,3) * t225 + Ifges(3,4) * t244 + Ifges(3,2) * t245 + Ifges(3,6) * qJDD(2) - pkin(2) * t275 + qJ(3) * t284 + qJD(2) * t237 + t261 * t118 + t260 * t122 - t235 * t288;
t116 = mrSges(3,2) * t238 - mrSges(3,3) * t224 + Ifges(3,1) * t244 + Ifges(3,4) * t245 + Ifges(3,5) * qJDD(2) - qJ(3) * t126 - qJD(2) * t236 - t118 * t260 + t122 * t261 + t235 * t287;
t272 = -m(3) * t238 + t245 * mrSges(3,1) - mrSges(3,2) * t244 - t248 * t288 + t249 * t287 - t275;
t278 = mrSges(2,1) * t250 - mrSges(2,2) * t251 + Ifges(2,3) * qJDD(1) + pkin(1) * t272 + pkin(6) * t285 + t268 * t114 + t264 * t116;
t136 = m(2) * t250 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t270 + t272;
t121 = t124 * t268 + t125 * t264;
t119 = m(2) * t251 - mrSges(2,1) * t270 - qJDD(1) * mrSges(2,2) + t285;
t117 = mrSges(2,1) * g(3) + mrSges(2,3) * t251 + t270 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t121 - t291;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t250 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t270 - pkin(6) * t121 - t114 * t264 + t116 * t268;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t112 - t265 * t117 - pkin(5) * (t119 * t265 + t136 * t269), t112, t116, t122, t127, t147; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t265 * t112 + t269 * t117 + pkin(5) * (t119 * t269 - t136 * t265), t117, t114, t118, t128, t146; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t278, t278, t291, -t273, -t276, t274;];
m_new = t1;
