% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:33
% EndTime: 2019-12-31 22:21:51
% DurationCPUTime: 10.74s
% Computational Cost: add. (177934->312), mult. (386478->396), div. (0->0), fcn. (278685->10), ass. (0->126)
t263 = sin(qJ(2));
t268 = cos(qJ(2));
t286 = qJD(1) * qJD(2);
t242 = t263 * qJDD(1) + t268 * t286;
t264 = sin(qJ(1));
t269 = cos(qJ(1));
t249 = -t269 * g(1) - t264 * g(2);
t270 = qJD(1) ^ 2;
t237 = -t270 * pkin(1) + qJDD(1) * pkin(6) + t249;
t289 = t263 * t237;
t290 = pkin(2) * t270;
t203 = qJDD(2) * pkin(2) - t242 * pkin(7) - t289 + (pkin(7) * t286 + t263 * t290 - g(3)) * t268;
t225 = -t263 * g(3) + t268 * t237;
t243 = t268 * qJDD(1) - t263 * t286;
t288 = qJD(1) * t263;
t247 = qJD(2) * pkin(2) - pkin(7) * t288;
t259 = t268 ^ 2;
t204 = t243 * pkin(7) - qJD(2) * t247 - t259 * t290 + t225;
t262 = sin(qJ(3));
t267 = cos(qJ(3));
t187 = t267 * t203 - t262 * t204;
t234 = (-t262 * t263 + t267 * t268) * qJD(1);
t211 = t234 * qJD(3) + t267 * t242 + t262 * t243;
t235 = (t262 * t268 + t263 * t267) * qJD(1);
t256 = qJDD(2) + qJDD(3);
t257 = qJD(2) + qJD(3);
t166 = (t234 * t257 - t211) * pkin(8) + (t234 * t235 + t256) * pkin(3) + t187;
t188 = t262 * t203 + t267 * t204;
t210 = -t235 * qJD(3) - t262 * t242 + t267 * t243;
t228 = t257 * pkin(3) - t235 * pkin(8);
t230 = t234 ^ 2;
t169 = -t230 * pkin(3) + t210 * pkin(8) - t257 * t228 + t188;
t261 = sin(qJ(4));
t266 = cos(qJ(4));
t164 = t261 * t166 + t266 * t169;
t222 = t261 * t234 + t266 * t235;
t183 = -t222 * qJD(4) + t266 * t210 - t261 * t211;
t221 = t266 * t234 - t261 * t235;
t197 = -t221 * mrSges(5,1) + t222 * mrSges(5,2);
t254 = qJD(4) + t257;
t214 = t254 * mrSges(5,1) - t222 * mrSges(5,3);
t253 = qJDD(4) + t256;
t198 = -t221 * pkin(4) - t222 * pkin(9);
t252 = t254 ^ 2;
t160 = -t252 * pkin(4) + t253 * pkin(9) + t221 * t198 + t164;
t248 = t264 * g(1) - t269 * g(2);
t280 = -qJDD(1) * pkin(1) - t248;
t212 = -t243 * pkin(2) + t247 * t288 + (-pkin(7) * t259 - pkin(6)) * t270 + t280;
t173 = -t210 * pkin(3) - t230 * pkin(8) + t235 * t228 + t212;
t184 = t221 * qJD(4) + t261 * t210 + t266 * t211;
t161 = (-t221 * t254 - t184) * pkin(9) + (t222 * t254 - t183) * pkin(4) + t173;
t260 = sin(qJ(5));
t265 = cos(qJ(5));
t157 = -t260 * t160 + t265 * t161;
t205 = -t260 * t222 + t265 * t254;
t171 = t205 * qJD(5) + t265 * t184 + t260 * t253;
t181 = qJDD(5) - t183;
t206 = t265 * t222 + t260 * t254;
t189 = -t205 * mrSges(6,1) + t206 * mrSges(6,2);
t218 = qJD(5) - t221;
t190 = -t218 * mrSges(6,2) + t205 * mrSges(6,3);
t153 = m(6) * t157 + t181 * mrSges(6,1) - t171 * mrSges(6,3) - t206 * t189 + t218 * t190;
t158 = t265 * t160 + t260 * t161;
t170 = -t206 * qJD(5) - t260 * t184 + t265 * t253;
t191 = t218 * mrSges(6,1) - t206 * mrSges(6,3);
t154 = m(6) * t158 - t181 * mrSges(6,2) + t170 * mrSges(6,3) + t205 * t189 - t218 * t191;
t282 = -t260 * t153 + t265 * t154;
t140 = m(5) * t164 - t253 * mrSges(5,2) + t183 * mrSges(5,3) + t221 * t197 - t254 * t214 + t282;
t163 = t266 * t166 - t261 * t169;
t213 = -t254 * mrSges(5,2) + t221 * mrSges(5,3);
t159 = -t253 * pkin(4) - t252 * pkin(9) + t222 * t198 - t163;
t277 = -m(6) * t159 + t170 * mrSges(6,1) - t171 * mrSges(6,2) + t205 * t190 - t206 * t191;
t149 = m(5) * t163 + t253 * mrSges(5,1) - t184 * mrSges(5,3) - t222 * t197 + t254 * t213 + t277;
t135 = t261 * t140 + t266 * t149;
t223 = -t234 * mrSges(4,1) + t235 * mrSges(4,2);
t226 = -t257 * mrSges(4,2) + t234 * mrSges(4,3);
t132 = m(4) * t187 + t256 * mrSges(4,1) - t211 * mrSges(4,3) - t235 * t223 + t257 * t226 + t135;
t227 = t257 * mrSges(4,1) - t235 * mrSges(4,3);
t283 = t266 * t140 - t261 * t149;
t133 = m(4) * t188 - t256 * mrSges(4,2) + t210 * mrSges(4,3) + t234 * t223 - t257 * t227 + t283;
t126 = t267 * t132 + t262 * t133;
t224 = -t268 * g(3) - t289;
t232 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t263 + Ifges(3,2) * t268) * qJD(1);
t233 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t263 + Ifges(3,4) * t268) * qJD(1);
t216 = Ifges(4,4) * t235 + Ifges(4,2) * t234 + Ifges(4,6) * t257;
t217 = Ifges(4,1) * t235 + Ifges(4,4) * t234 + Ifges(4,5) * t257;
t174 = Ifges(6,5) * t206 + Ifges(6,6) * t205 + Ifges(6,3) * t218;
t176 = Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t218;
t146 = -mrSges(6,1) * t159 + mrSges(6,3) * t158 + Ifges(6,4) * t171 + Ifges(6,2) * t170 + Ifges(6,6) * t181 - t206 * t174 + t218 * t176;
t175 = Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t218;
t147 = mrSges(6,2) * t159 - mrSges(6,3) * t157 + Ifges(6,1) * t171 + Ifges(6,4) * t170 + Ifges(6,5) * t181 + t205 * t174 - t218 * t175;
t193 = Ifges(5,4) * t222 + Ifges(5,2) * t221 + Ifges(5,6) * t254;
t194 = Ifges(5,1) * t222 + Ifges(5,4) * t221 + Ifges(5,5) * t254;
t276 = -mrSges(5,1) * t163 + mrSges(5,2) * t164 - Ifges(5,5) * t184 - Ifges(5,6) * t183 - Ifges(5,3) * t253 - pkin(4) * t277 - pkin(9) * t282 - t265 * t146 - t260 * t147 - t222 * t193 + t221 * t194;
t273 = -mrSges(4,1) * t187 + mrSges(4,2) * t188 - Ifges(4,5) * t211 - Ifges(4,6) * t210 - Ifges(4,3) * t256 - pkin(3) * t135 - t235 * t216 + t234 * t217 + t276;
t291 = mrSges(3,1) * t224 - mrSges(3,2) * t225 + Ifges(3,5) * t242 + Ifges(3,6) * t243 + Ifges(3,3) * qJDD(2) + pkin(2) * t126 + (t263 * t232 - t268 * t233) * qJD(1) - t273;
t142 = t265 * t153 + t260 * t154;
t287 = qJD(1) * t268;
t241 = (-mrSges(3,1) * t268 + mrSges(3,2) * t263) * qJD(1);
t246 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t287;
t124 = m(3) * t224 + qJDD(2) * mrSges(3,1) - t242 * mrSges(3,3) + qJD(2) * t246 - t241 * t288 + t126;
t245 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t288;
t284 = -t262 * t132 + t267 * t133;
t125 = m(3) * t225 - qJDD(2) * mrSges(3,2) + t243 * mrSges(3,3) - qJD(2) * t245 + t241 * t287 + t284;
t285 = -t263 * t124 + t268 * t125;
t279 = m(5) * t173 - t183 * mrSges(5,1) + t184 * mrSges(5,2) - t221 * t213 + t222 * t214 + t142;
t192 = Ifges(5,5) * t222 + Ifges(5,6) * t221 + Ifges(5,3) * t254;
t127 = mrSges(5,2) * t173 - mrSges(5,3) * t163 + Ifges(5,1) * t184 + Ifges(5,4) * t183 + Ifges(5,5) * t253 - pkin(9) * t142 - t260 * t146 + t265 * t147 + t221 * t192 - t254 * t193;
t274 = mrSges(6,1) * t157 - mrSges(6,2) * t158 + Ifges(6,5) * t171 + Ifges(6,6) * t170 + Ifges(6,3) * t181 + t206 * t175 - t205 * t176;
t128 = -mrSges(5,1) * t173 + mrSges(5,3) * t164 + Ifges(5,4) * t184 + Ifges(5,2) * t183 + Ifges(5,6) * t253 - pkin(4) * t142 - t222 * t192 + t254 * t194 - t274;
t215 = Ifges(4,5) * t235 + Ifges(4,6) * t234 + Ifges(4,3) * t257;
t118 = -mrSges(4,1) * t212 + mrSges(4,3) * t188 + Ifges(4,4) * t211 + Ifges(4,2) * t210 + Ifges(4,6) * t256 - pkin(3) * t279 + pkin(8) * t283 + t261 * t127 + t266 * t128 - t235 * t215 + t257 * t217;
t122 = mrSges(4,2) * t212 - mrSges(4,3) * t187 + Ifges(4,1) * t211 + Ifges(4,4) * t210 + Ifges(4,5) * t256 - pkin(8) * t135 + t266 * t127 - t261 * t128 + t234 * t215 - t257 * t216;
t231 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t263 + Ifges(3,6) * t268) * qJD(1);
t236 = -t270 * pkin(6) + t280;
t275 = m(4) * t212 - t210 * mrSges(4,1) + t211 * mrSges(4,2) - t234 * t226 + t235 * t227 + t279;
t114 = -mrSges(3,1) * t236 + mrSges(3,3) * t225 + Ifges(3,4) * t242 + Ifges(3,2) * t243 + Ifges(3,6) * qJDD(2) - pkin(2) * t275 + pkin(7) * t284 + qJD(2) * t233 + t267 * t118 + t262 * t122 - t231 * t288;
t116 = mrSges(3,2) * t236 - mrSges(3,3) * t224 + Ifges(3,1) * t242 + Ifges(3,4) * t243 + Ifges(3,5) * qJDD(2) - pkin(7) * t126 - qJD(2) * t232 - t262 * t118 + t267 * t122 + t231 * t287;
t272 = -m(3) * t236 + t243 * mrSges(3,1) - t242 * mrSges(3,2) - t245 * t288 + t246 * t287 - t275;
t278 = mrSges(2,1) * t248 - mrSges(2,2) * t249 + Ifges(2,3) * qJDD(1) + pkin(1) * t272 + pkin(6) * t285 + t268 * t114 + t263 * t116;
t136 = m(2) * t248 + qJDD(1) * mrSges(2,1) - t270 * mrSges(2,2) + t272;
t121 = t268 * t124 + t263 * t125;
t119 = m(2) * t249 - t270 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t285;
t117 = mrSges(2,1) * g(3) + mrSges(2,3) * t249 + t270 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t121 - t291;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t248 + Ifges(2,5) * qJDD(1) - t270 * Ifges(2,6) - pkin(6) * t121 - t263 * t114 + t268 * t116;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t112 - t264 * t117 - pkin(5) * (t264 * t119 + t269 * t136), t112, t116, t122, t127, t147; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t112 + t269 * t117 + pkin(5) * (t269 * t119 - t264 * t136), t117, t114, t118, t128, t146; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t278, t278, t291, -t273, -t276, t274;];
m_new = t1;
