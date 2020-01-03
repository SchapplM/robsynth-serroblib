% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:51
% EndTime: 2019-12-31 19:57:02
% DurationCPUTime: 5.48s
% Computational Cost: add. (62525->306), mult. (141434->379), div. (0->0), fcn. (93235->8), ass. (0->115)
t301 = -2 * qJD(3);
t263 = sin(qJ(2));
t266 = cos(qJ(2));
t290 = qJD(1) * qJD(2);
t248 = t263 * qJDD(1) + t266 * t290;
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t255 = -t267 * g(1) - t264 * g(2);
t269 = qJD(1) ^ 2;
t243 = -t269 * pkin(1) + qJDD(1) * pkin(6) + t255;
t296 = t263 * t243;
t298 = pkin(2) * t269;
t199 = qJDD(2) * pkin(2) - t248 * qJ(3) - t296 + (qJ(3) * t290 + t263 * t298 - g(3)) * t266;
t230 = -t263 * g(3) + t266 * t243;
t249 = t266 * qJDD(1) - t263 * t290;
t293 = qJD(1) * t263;
t251 = qJD(2) * pkin(2) - qJ(3) * t293;
t259 = t266 ^ 2;
t200 = t249 * qJ(3) - qJD(2) * t251 - t259 * t298 + t230;
t260 = sin(pkin(8));
t261 = cos(pkin(8));
t238 = (t260 * t266 + t261 * t263) * qJD(1);
t169 = t261 * t199 - t260 * t200 + t238 * t301;
t237 = (t260 * t263 - t261 * t266) * qJD(1);
t225 = t261 * t248 + t260 * t249;
t262 = sin(qJ(4));
t265 = cos(qJ(4));
t227 = t265 * qJD(2) - t262 * t238;
t193 = t227 * qJD(4) + t262 * qJDD(2) + t265 * t225;
t228 = t262 * qJD(2) + t265 * t238;
t201 = -t227 * mrSges(6,1) + t228 * mrSges(6,2);
t170 = t260 * t199 + t261 * t200 + t237 * t301;
t215 = t237 * pkin(3) - t238 * pkin(7);
t268 = qJD(2) ^ 2;
t164 = -t268 * pkin(3) + qJDD(2) * pkin(7) - t237 * t215 + t170;
t254 = t264 * g(1) - t267 * g(2);
t280 = -qJDD(1) * pkin(1) - t254;
t204 = -t249 * pkin(2) + qJDD(3) + t251 * t293 + (-qJ(3) * t259 - pkin(6)) * t269 + t280;
t224 = -t260 * t248 + t261 * t249;
t167 = (qJD(2) * t237 - t225) * pkin(7) + (qJD(2) * t238 - t224) * pkin(3) + t204;
t160 = -t262 * t164 + t265 * t167;
t223 = qJDD(4) - t224;
t236 = qJD(4) + t237;
t154 = -0.2e1 * qJD(5) * t228 + (t227 * t236 - t193) * qJ(5) + (t227 * t228 + t223) * pkin(4) + t160;
t205 = -t236 * mrSges(6,2) + t227 * mrSges(6,3);
t289 = m(6) * t154 + t223 * mrSges(6,1) + t236 * t205;
t151 = -t193 * mrSges(6,3) - t228 * t201 + t289;
t161 = t265 * t164 + t262 * t167;
t178 = Ifges(5,4) * t228 + Ifges(5,2) * t227 + Ifges(5,6) * t236;
t179 = Ifges(6,1) * t228 + Ifges(6,4) * t227 + Ifges(6,5) * t236;
t180 = Ifges(5,1) * t228 + Ifges(5,4) * t227 + Ifges(5,5) * t236;
t192 = -t228 * qJD(4) + t265 * qJDD(2) - t262 * t225;
t207 = t236 * pkin(4) - t228 * qJ(5);
t226 = t227 ^ 2;
t157 = -t226 * pkin(4) + t192 * qJ(5) + 0.2e1 * qJD(5) * t227 - t236 * t207 + t161;
t177 = Ifges(6,4) * t228 + Ifges(6,2) * t227 + Ifges(6,6) * t236;
t278 = -mrSges(6,1) * t154 + mrSges(6,2) * t157 - Ifges(6,5) * t193 - Ifges(6,6) * t192 - Ifges(6,3) * t223 - t228 * t177;
t300 = mrSges(5,1) * t160 - mrSges(5,2) * t161 + Ifges(5,5) * t193 + Ifges(5,6) * t192 + Ifges(5,3) * t223 + pkin(4) * t151 + t228 * t178 - (t180 + t179) * t227 - t278;
t214 = t237 * mrSges(4,1) + t238 * mrSges(4,2);
t232 = qJD(2) * mrSges(4,1) - t238 * mrSges(4,3);
t202 = -t227 * mrSges(5,1) + t228 * mrSges(5,2);
t206 = -t236 * mrSges(5,2) + t227 * mrSges(5,3);
t143 = m(5) * t160 + t223 * mrSges(5,1) + t236 * t206 + (-t201 - t202) * t228 + (-mrSges(5,3) - mrSges(6,3)) * t193 + t289;
t288 = m(6) * t157 + t192 * mrSges(6,3) + t227 * t201;
t208 = t236 * mrSges(6,1) - t228 * mrSges(6,3);
t294 = -t236 * mrSges(5,1) + t228 * mrSges(5,3) - t208;
t297 = -mrSges(5,2) - mrSges(6,2);
t146 = m(5) * t161 + t192 * mrSges(5,3) + t227 * t202 + t297 * t223 + t294 * t236 + t288;
t285 = -t262 * t143 + t265 * t146;
t136 = m(4) * t170 - qJDD(2) * mrSges(4,2) + t224 * mrSges(4,3) - qJD(2) * t232 - t237 * t214 + t285;
t231 = -qJD(2) * mrSges(4,2) - t237 * mrSges(4,3);
t163 = -qJDD(2) * pkin(3) - t268 * pkin(7) + t238 * t215 - t169;
t159 = -t192 * pkin(4) - t226 * qJ(5) + t228 * t207 + qJDD(5) + t163;
t283 = -m(6) * t159 + t192 * mrSges(6,1) + t227 * t205;
t273 = -m(5) * t163 + t192 * mrSges(5,1) + t297 * t193 + t227 * t206 + t294 * t228 + t283;
t148 = m(4) * t169 + qJDD(2) * mrSges(4,1) - t225 * mrSges(4,3) + qJD(2) * t231 - t238 * t214 + t273;
t129 = t260 * t136 + t261 * t148;
t229 = -t266 * g(3) - t296;
t240 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t263 + Ifges(3,2) * t266) * qJD(1);
t241 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t263 + Ifges(3,4) * t266) * qJD(1);
t175 = Ifges(6,5) * t228 + Ifges(6,6) * t227 + Ifges(6,3) * t236;
t176 = Ifges(5,5) * t228 + Ifges(5,6) * t227 + Ifges(5,3) * t236;
t279 = -mrSges(6,1) * t159 + mrSges(6,3) * t157 + Ifges(6,4) * t193 + Ifges(6,2) * t192 + Ifges(6,6) * t223 + t236 * t179;
t131 = Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * t223 + t236 * t180 - mrSges(5,1) * t163 + mrSges(5,3) * t161 - pkin(4) * (t193 * mrSges(6,2) - t283) + qJ(5) * (-t223 * mrSges(6,2) - t236 * t208 + t288) + (-pkin(4) * t208 - t175 - t176) * t228 + t279;
t277 = mrSges(6,2) * t159 - mrSges(6,3) * t154 + Ifges(6,1) * t193 + Ifges(6,4) * t192 + Ifges(6,5) * t223 + t227 * t175;
t138 = mrSges(5,2) * t163 - mrSges(5,3) * t160 + Ifges(5,1) * t193 + Ifges(5,4) * t192 + Ifges(5,5) * t223 - qJ(5) * t151 + t227 * t176 + (-t177 - t178) * t236 + t277;
t211 = Ifges(4,4) * t238 - Ifges(4,2) * t237 + Ifges(4,6) * qJD(2);
t212 = Ifges(4,1) * t238 - Ifges(4,4) * t237 + Ifges(4,5) * qJD(2);
t274 = -mrSges(4,1) * t169 + mrSges(4,2) * t170 - Ifges(4,5) * t225 - Ifges(4,6) * t224 - Ifges(4,3) * qJDD(2) - pkin(3) * t273 - pkin(7) * t285 - t265 * t131 - t262 * t138 - t238 * t211 - t237 * t212;
t299 = mrSges(3,1) * t229 - mrSges(3,2) * t230 + Ifges(3,5) * t248 + Ifges(3,6) * t249 + Ifges(3,3) * qJDD(2) + pkin(2) * t129 + (t263 * t240 - t266 * t241) * qJD(1) - t274;
t140 = t265 * t143 + t262 * t146;
t292 = qJD(1) * t266;
t247 = (-mrSges(3,1) * t266 + mrSges(3,2) * t263) * qJD(1);
t253 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t292;
t127 = m(3) * t229 + qJDD(2) * mrSges(3,1) - t248 * mrSges(3,3) + qJD(2) * t253 - t247 * t293 + t129;
t252 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t293;
t286 = t261 * t136 - t260 * t148;
t128 = m(3) * t230 - qJDD(2) * mrSges(3,2) + t249 * mrSges(3,3) - qJD(2) * t252 + t247 * t292 + t286;
t287 = -t263 * t127 + t266 * t128;
t210 = Ifges(4,5) * t238 - Ifges(4,6) * t237 + Ifges(4,3) * qJD(2);
t121 = mrSges(4,2) * t204 - mrSges(4,3) * t169 + Ifges(4,1) * t225 + Ifges(4,4) * t224 + Ifges(4,5) * qJDD(2) - pkin(7) * t140 - qJD(2) * t211 - t262 * t131 + t265 * t138 - t237 * t210;
t125 = -mrSges(4,1) * t204 + mrSges(4,3) * t170 + Ifges(4,4) * t225 + Ifges(4,2) * t224 + Ifges(4,6) * qJDD(2) - pkin(3) * t140 + qJD(2) * t212 - t238 * t210 - t300;
t239 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t263 + Ifges(3,6) * t266) * qJD(1);
t242 = -t269 * pkin(6) + t280;
t275 = m(4) * t204 - t224 * mrSges(4,1) + t225 * mrSges(4,2) + t237 * t231 + t238 * t232 + t140;
t118 = -mrSges(3,1) * t242 + mrSges(3,3) * t230 + Ifges(3,4) * t248 + Ifges(3,2) * t249 + Ifges(3,6) * qJDD(2) - pkin(2) * t275 + qJ(3) * t286 + qJD(2) * t241 + t260 * t121 + t261 * t125 - t239 * t293;
t120 = mrSges(3,2) * t242 - mrSges(3,3) * t229 + Ifges(3,1) * t248 + Ifges(3,4) * t249 + Ifges(3,5) * qJDD(2) - qJ(3) * t129 - qJD(2) * t240 + t261 * t121 - t260 * t125 + t239 * t292;
t272 = -m(3) * t242 + t249 * mrSges(3,1) - t248 * mrSges(3,2) - t252 * t293 + t253 * t292 - t275;
t276 = mrSges(2,1) * t254 - mrSges(2,2) * t255 + Ifges(2,3) * qJDD(1) + pkin(1) * t272 + pkin(6) * t287 + t266 * t118 + t263 * t120;
t132 = m(2) * t254 + qJDD(1) * mrSges(2,1) - t269 * mrSges(2,2) + t272;
t124 = t266 * t127 + t263 * t128;
t122 = m(2) * t255 - t269 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t287;
t116 = mrSges(2,1) * g(3) + mrSges(2,3) * t255 + t269 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t124 - t299;
t115 = -mrSges(2,2) * g(3) - mrSges(2,3) * t254 + Ifges(2,5) * qJDD(1) - t269 * Ifges(2,6) - pkin(6) * t124 - t263 * t118 + t266 * t120;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t267 * t115 - t264 * t116 - pkin(5) * (t264 * t122 + t267 * t132), t115, t120, t121, t138, -t236 * t177 + t277; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t115 + t267 * t116 + pkin(5) * (t267 * t122 - t264 * t132), t116, t118, t125, t131, -t228 * t175 + t279; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t276, t276, t299, -t274, t300, -t227 * t179 - t278;];
m_new = t1;
