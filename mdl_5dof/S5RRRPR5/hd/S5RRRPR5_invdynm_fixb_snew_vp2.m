% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:11
% EndTime: 2019-12-31 21:13:27
% DurationCPUTime: 9.84s
% Computational Cost: add. (164599->312), mult. (369064->398), div. (0->0), fcn. (261995->10), ass. (0->126)
t263 = sin(qJ(2));
t267 = cos(qJ(2));
t286 = qJD(1) * qJD(2);
t243 = t263 * qJDD(1) + t267 * t286;
t264 = sin(qJ(1));
t268 = cos(qJ(1));
t250 = -t268 * g(1) - t264 * g(2);
t269 = qJD(1) ^ 2;
t238 = -t269 * pkin(1) + qJDD(1) * pkin(6) + t250;
t289 = t263 * t238;
t290 = pkin(2) * t269;
t203 = qJDD(2) * pkin(2) - t243 * pkin(7) - t289 + (pkin(7) * t286 + t263 * t290 - g(3)) * t267;
t226 = -t263 * g(3) + t267 * t238;
t244 = t267 * qJDD(1) - t263 * t286;
t288 = qJD(1) * t263;
t248 = qJD(2) * pkin(2) - pkin(7) * t288;
t258 = t267 ^ 2;
t204 = t244 * pkin(7) - qJD(2) * t248 - t258 * t290 + t226;
t262 = sin(qJ(3));
t266 = cos(qJ(3));
t180 = t266 * t203 - t262 * t204;
t235 = (-t262 * t263 + t266 * t267) * qJD(1);
t211 = t235 * qJD(3) + t266 * t243 + t262 * t244;
t236 = (t262 * t267 + t263 * t266) * qJD(1);
t255 = qJDD(2) + qJDD(3);
t256 = qJD(2) + qJD(3);
t166 = (t235 * t256 - t211) * qJ(4) + (t235 * t236 + t255) * pkin(3) + t180;
t181 = t262 * t203 + t266 * t204;
t210 = -t236 * qJD(3) - t262 * t243 + t266 * t244;
t228 = t256 * pkin(3) - t236 * qJ(4);
t231 = t235 ^ 2;
t168 = -t231 * pkin(3) + t210 * qJ(4) - t256 * t228 + t181;
t259 = sin(pkin(9));
t260 = cos(pkin(9));
t222 = t260 * t235 - t259 * t236;
t291 = 2 * qJD(4);
t163 = t259 * t166 + t260 * t168 + t222 * t291;
t187 = t260 * t210 - t259 * t211;
t223 = t259 * t235 + t260 * t236;
t197 = -t222 * mrSges(5,1) + t223 * mrSges(5,2);
t214 = t256 * mrSges(5,1) - t223 * mrSges(5,3);
t198 = -t222 * pkin(4) - t223 * pkin(8);
t254 = t256 ^ 2;
t160 = -t254 * pkin(4) + t255 * pkin(8) + t222 * t198 + t163;
t249 = t264 * g(1) - t268 * g(2);
t279 = -qJDD(1) * pkin(1) - t249;
t212 = -t244 * pkin(2) + t248 * t288 + (-pkin(7) * t258 - pkin(6)) * t269 + t279;
t173 = -t210 * pkin(3) - t231 * qJ(4) + t236 * t228 + qJDD(4) + t212;
t188 = t259 * t210 + t260 * t211;
t164 = (-t222 * t256 - t188) * pkin(8) + (t223 * t256 - t187) * pkin(4) + t173;
t261 = sin(qJ(5));
t265 = cos(qJ(5));
t157 = -t261 * t160 + t265 * t164;
t208 = -t261 * t223 + t265 * t256;
t171 = t208 * qJD(5) + t265 * t188 + t261 * t255;
t186 = qJDD(5) - t187;
t209 = t265 * t223 + t261 * t256;
t189 = -t208 * mrSges(6,1) + t209 * mrSges(6,2);
t216 = qJD(5) - t222;
t190 = -t216 * mrSges(6,2) + t208 * mrSges(6,3);
t153 = m(6) * t157 + t186 * mrSges(6,1) - t171 * mrSges(6,3) - t209 * t189 + t216 * t190;
t158 = t265 * t160 + t261 * t164;
t170 = -t209 * qJD(5) - t261 * t188 + t265 * t255;
t191 = t216 * mrSges(6,1) - t209 * mrSges(6,3);
t154 = m(6) * t158 - t186 * mrSges(6,2) + t170 * mrSges(6,3) + t208 * t189 - t216 * t191;
t282 = -t261 * t153 + t265 * t154;
t140 = m(5) * t163 - t255 * mrSges(5,2) + t187 * mrSges(5,3) + t222 * t197 - t256 * t214 + t282;
t281 = -t260 * t166 + t259 * t168;
t162 = -0.2e1 * qJD(4) * t223 - t281;
t213 = -t256 * mrSges(5,2) + t222 * mrSges(5,3);
t159 = -t255 * pkin(4) - t254 * pkin(8) + (t291 + t198) * t223 + t281;
t276 = -m(6) * t159 + t170 * mrSges(6,1) - t171 * mrSges(6,2) + t208 * t190 - t209 * t191;
t149 = m(5) * t162 + t255 * mrSges(5,1) - t188 * mrSges(5,3) - t223 * t197 + t256 * t213 + t276;
t135 = t259 * t140 + t260 * t149;
t224 = -t235 * mrSges(4,1) + t236 * mrSges(4,2);
t227 = -t256 * mrSges(4,2) + t235 * mrSges(4,3);
t132 = m(4) * t180 + t255 * mrSges(4,1) - t211 * mrSges(4,3) - t236 * t224 + t256 * t227 + t135;
t229 = t256 * mrSges(4,1) - t236 * mrSges(4,3);
t283 = t260 * t140 - t259 * t149;
t133 = m(4) * t181 - t255 * mrSges(4,2) + t210 * mrSges(4,3) + t235 * t224 - t256 * t229 + t283;
t126 = t266 * t132 + t262 * t133;
t225 = -t267 * g(3) - t289;
t233 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t263 + Ifges(3,2) * t267) * qJD(1);
t234 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t263 + Ifges(3,4) * t267) * qJD(1);
t218 = Ifges(4,4) * t236 + Ifges(4,2) * t235 + Ifges(4,6) * t256;
t219 = Ifges(4,1) * t236 + Ifges(4,4) * t235 + Ifges(4,5) * t256;
t174 = Ifges(6,5) * t209 + Ifges(6,6) * t208 + Ifges(6,3) * t216;
t176 = Ifges(6,1) * t209 + Ifges(6,4) * t208 + Ifges(6,5) * t216;
t146 = -mrSges(6,1) * t159 + mrSges(6,3) * t158 + Ifges(6,4) * t171 + Ifges(6,2) * t170 + Ifges(6,6) * t186 - t209 * t174 + t216 * t176;
t175 = Ifges(6,4) * t209 + Ifges(6,2) * t208 + Ifges(6,6) * t216;
t147 = mrSges(6,2) * t159 - mrSges(6,3) * t157 + Ifges(6,1) * t171 + Ifges(6,4) * t170 + Ifges(6,5) * t186 + t208 * t174 - t216 * t175;
t193 = Ifges(5,4) * t223 + Ifges(5,2) * t222 + Ifges(5,6) * t256;
t194 = Ifges(5,1) * t223 + Ifges(5,4) * t222 + Ifges(5,5) * t256;
t275 = -mrSges(5,1) * t162 + mrSges(5,2) * t163 - Ifges(5,5) * t188 - Ifges(5,6) * t187 - Ifges(5,3) * t255 - pkin(4) * t276 - pkin(8) * t282 - t265 * t146 - t261 * t147 - t223 * t193 + t222 * t194;
t272 = -mrSges(4,1) * t180 + mrSges(4,2) * t181 - Ifges(4,5) * t211 - Ifges(4,6) * t210 - Ifges(4,3) * t255 - pkin(3) * t135 - t236 * t218 + t235 * t219 + t275;
t292 = mrSges(3,1) * t225 - mrSges(3,2) * t226 + Ifges(3,5) * t243 + Ifges(3,6) * t244 + Ifges(3,3) * qJDD(2) + pkin(2) * t126 + (t263 * t233 - t267 * t234) * qJD(1) - t272;
t142 = t265 * t153 + t261 * t154;
t287 = qJD(1) * t267;
t242 = (-mrSges(3,1) * t267 + mrSges(3,2) * t263) * qJD(1);
t247 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t287;
t124 = m(3) * t225 + qJDD(2) * mrSges(3,1) - t243 * mrSges(3,3) + qJD(2) * t247 - t242 * t288 + t126;
t246 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t288;
t284 = -t262 * t132 + t266 * t133;
t125 = m(3) * t226 - qJDD(2) * mrSges(3,2) + t244 * mrSges(3,3) - qJD(2) * t246 + t242 * t287 + t284;
t285 = -t263 * t124 + t267 * t125;
t278 = m(5) * t173 - t187 * mrSges(5,1) + t188 * mrSges(5,2) - t222 * t213 + t223 * t214 + t142;
t192 = Ifges(5,5) * t223 + Ifges(5,6) * t222 + Ifges(5,3) * t256;
t127 = mrSges(5,2) * t173 - mrSges(5,3) * t162 + Ifges(5,1) * t188 + Ifges(5,4) * t187 + Ifges(5,5) * t255 - pkin(8) * t142 - t261 * t146 + t265 * t147 + t222 * t192 - t256 * t193;
t273 = mrSges(6,1) * t157 - mrSges(6,2) * t158 + Ifges(6,5) * t171 + Ifges(6,6) * t170 + Ifges(6,3) * t186 + t209 * t175 - t208 * t176;
t128 = -mrSges(5,1) * t173 + mrSges(5,3) * t163 + Ifges(5,4) * t188 + Ifges(5,2) * t187 + Ifges(5,6) * t255 - pkin(4) * t142 - t223 * t192 + t256 * t194 - t273;
t217 = Ifges(4,5) * t236 + Ifges(4,6) * t235 + Ifges(4,3) * t256;
t118 = -mrSges(4,1) * t212 + mrSges(4,3) * t181 + Ifges(4,4) * t211 + Ifges(4,2) * t210 + Ifges(4,6) * t255 - pkin(3) * t278 + qJ(4) * t283 + t259 * t127 + t260 * t128 - t236 * t217 + t256 * t219;
t122 = mrSges(4,2) * t212 - mrSges(4,3) * t180 + Ifges(4,1) * t211 + Ifges(4,4) * t210 + Ifges(4,5) * t255 - qJ(4) * t135 + t260 * t127 - t259 * t128 + t235 * t217 - t256 * t218;
t232 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t263 + Ifges(3,6) * t267) * qJD(1);
t237 = -t269 * pkin(6) + t279;
t274 = m(4) * t212 - t210 * mrSges(4,1) + t211 * mrSges(4,2) - t235 * t227 + t236 * t229 + t278;
t114 = -mrSges(3,1) * t237 + mrSges(3,3) * t226 + Ifges(3,4) * t243 + Ifges(3,2) * t244 + Ifges(3,6) * qJDD(2) - pkin(2) * t274 + pkin(7) * t284 + qJD(2) * t234 + t266 * t118 + t262 * t122 - t232 * t288;
t116 = mrSges(3,2) * t237 - mrSges(3,3) * t225 + Ifges(3,1) * t243 + Ifges(3,4) * t244 + Ifges(3,5) * qJDD(2) - pkin(7) * t126 - qJD(2) * t233 - t262 * t118 + t266 * t122 + t232 * t287;
t271 = -m(3) * t237 + t244 * mrSges(3,1) - t243 * mrSges(3,2) - t246 * t288 + t247 * t287 - t274;
t277 = mrSges(2,1) * t249 - mrSges(2,2) * t250 + Ifges(2,3) * qJDD(1) + pkin(1) * t271 + pkin(6) * t285 + t267 * t114 + t263 * t116;
t136 = m(2) * t249 + qJDD(1) * mrSges(2,1) - t269 * mrSges(2,2) + t271;
t121 = t267 * t124 + t263 * t125;
t119 = m(2) * t250 - t269 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t285;
t117 = mrSges(2,1) * g(3) + mrSges(2,3) * t250 + t269 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t121 - t292;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t249 + Ifges(2,5) * qJDD(1) - t269 * Ifges(2,6) - pkin(6) * t121 - t263 * t114 + t267 * t116;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t268 * t112 - t264 * t117 - pkin(5) * (t264 * t119 + t268 * t136), t112, t116, t122, t127, t147; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t112 + t268 * t117 + pkin(5) * (t268 * t119 - t264 * t136), t117, t114, t118, t128, t146; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t277, t277, t292, -t272, -t275, t273;];
m_new = t1;
