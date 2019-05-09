% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-05-04 19:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:25:52
% EndTime: 2019-05-04 19:26:03
% DurationCPUTime: 5.77s
% Computational Cost: add. (93275->290), mult. (204517->373), div. (0->0), fcn. (155565->10), ass. (0->116)
t251 = sin(qJ(1));
t256 = cos(qJ(1));
t233 = -g(1) * t256 - g(2) * t251;
t257 = qJD(1) ^ 2;
t226 = -pkin(1) * t257 + t233;
t250 = sin(qJ(2));
t255 = cos(qJ(2));
t210 = t255 * g(3) - t250 * t226;
t204 = (t250 * t255 * t257 + qJDD(2)) * pkin(2) + t210;
t211 = t250 * g(3) + t255 * t226;
t205 = (-t255 ^ 2 * t257 - qJD(2) ^ 2) * pkin(2) + t211;
t249 = sin(qJ(3));
t254 = cos(qJ(3));
t179 = t254 * t204 - t249 * t205;
t216 = (t249 * t250 - t254 * t255) * qJD(1);
t269 = qJD(1) * qJD(2);
t227 = -qJDD(1) * t250 - t255 * t269;
t268 = t250 * t269;
t228 = -qJDD(1) * t255 + t268;
t188 = qJD(3) * t216 + t227 * t254 + t228 * t249;
t217 = (-t249 * t255 - t250 * t254) * qJD(1);
t198 = -mrSges(4,1) * t216 + mrSges(4,2) * t217;
t243 = qJD(2) + qJD(3);
t206 = -mrSges(4,2) * t243 + mrSges(4,3) * t216;
t242 = qJDD(2) + qJDD(3);
t164 = (t216 * t217 + t242) * pkin(3) + t179;
t180 = t249 * t204 + t254 * t205;
t171 = (-t216 ^ 2 - t243 ^ 2) * pkin(3) + t180;
t248 = sin(qJ(4));
t253 = cos(qJ(4));
t151 = t248 * t164 + t253 * t171;
t187 = -qJD(3) * t217 - t227 * t249 + t228 * t254;
t197 = t216 * t248 + t217 * t253;
t160 = -qJD(4) * t197 + t187 * t253 - t188 * t248;
t196 = t216 * t253 - t217 * t248;
t176 = -mrSges(5,1) * t196 + mrSges(5,2) * t197;
t238 = qJD(4) + t243;
t190 = mrSges(5,1) * t238 - mrSges(5,3) * t197;
t237 = qJDD(4) + t242;
t161 = qJD(4) * t196 + t187 * t248 + t188 * t253;
t232 = t251 * g(1) - t256 * g(2);
t224 = qJDD(1) * pkin(1) + t232;
t202 = (-t228 - t268) * pkin(2) + t224;
t167 = t202 + (t217 * t243 - t187) * pkin(3);
t143 = (-t196 * t238 - t161) * pkin(6) + (t197 * t238 - t160) * pkin(4) + t167;
t177 = -pkin(4) * t196 - pkin(6) * t197;
t236 = t238 ^ 2;
t145 = -pkin(4) * t236 + pkin(6) * t237 + t177 * t196 + t151;
t247 = sin(qJ(5));
t252 = cos(qJ(5));
t141 = t143 * t252 - t145 * t247;
t182 = -t197 * t247 + t238 * t252;
t148 = qJD(5) * t182 + t161 * t252 + t237 * t247;
t158 = qJDD(5) - t160;
t183 = t197 * t252 + t238 * t247;
t165 = -mrSges(6,1) * t182 + mrSges(6,2) * t183;
t194 = qJD(5) - t196;
t168 = -mrSges(6,2) * t194 + mrSges(6,3) * t182;
t138 = m(6) * t141 + mrSges(6,1) * t158 - mrSges(6,3) * t148 - t165 * t183 + t168 * t194;
t142 = t143 * t247 + t145 * t252;
t147 = -qJD(5) * t183 - t161 * t247 + t237 * t252;
t169 = mrSges(6,1) * t194 - mrSges(6,3) * t183;
t139 = m(6) * t142 - mrSges(6,2) * t158 + mrSges(6,3) * t147 + t165 * t182 - t169 * t194;
t267 = -t138 * t247 + t252 * t139;
t126 = m(5) * t151 - mrSges(5,2) * t237 + mrSges(5,3) * t160 + t176 * t196 - t190 * t238 + t267;
t150 = t164 * t253 - t171 * t248;
t189 = -mrSges(5,2) * t238 + mrSges(5,3) * t196;
t144 = -pkin(4) * t237 - pkin(6) * t236 + t177 * t197 - t150;
t266 = -m(6) * t144 + t147 * mrSges(6,1) - mrSges(6,2) * t148 + t182 * t168 - t169 * t183;
t134 = m(5) * t150 + mrSges(5,1) * t237 - mrSges(5,3) * t161 - t176 * t197 + t189 * t238 + t266;
t272 = t248 * t126 + t253 * t134;
t120 = m(4) * t179 + mrSges(4,1) * t242 - mrSges(4,3) * t188 - t198 * t217 + t206 * t243 + t272;
t207 = mrSges(4,1) * t243 - mrSges(4,3) * t217;
t121 = m(4) * t180 - mrSges(4,2) * t242 + mrSges(4,3) * t187 + t126 * t253 - t134 * t248 + t198 * t216 - t207 * t243;
t273 = t254 * t120 + t249 * t121;
t128 = t252 * t138 + t247 * t139;
t271 = qJD(1) * t250;
t270 = qJD(1) * t255;
t265 = m(5) * t167 - mrSges(5,1) * t160 + t161 * mrSges(5,2) - t189 * t196 + t197 * t190 + t128;
t152 = Ifges(6,5) * t183 + Ifges(6,6) * t182 + Ifges(6,3) * t194;
t154 = Ifges(6,1) * t183 + Ifges(6,4) * t182 + Ifges(6,5) * t194;
t131 = -mrSges(6,1) * t144 + mrSges(6,3) * t142 + Ifges(6,4) * t148 + Ifges(6,2) * t147 + Ifges(6,6) * t158 - t152 * t183 + t154 * t194;
t153 = Ifges(6,4) * t183 + Ifges(6,2) * t182 + Ifges(6,6) * t194;
t132 = mrSges(6,2) * t144 - mrSges(6,3) * t141 + Ifges(6,1) * t148 + Ifges(6,4) * t147 + Ifges(6,5) * t158 + t152 * t182 - t153 * t194;
t172 = Ifges(5,5) * t197 + Ifges(5,6) * t196 + Ifges(5,3) * t238;
t173 = Ifges(5,4) * t197 + Ifges(5,2) * t196 + Ifges(5,6) * t238;
t116 = mrSges(5,2) * t167 - mrSges(5,3) * t150 + Ifges(5,1) * t161 + Ifges(5,4) * t160 + Ifges(5,5) * t237 - pkin(6) * t128 - t131 * t247 + t132 * t252 + t172 * t196 - t173 * t238;
t174 = Ifges(5,1) * t197 + Ifges(5,4) * t196 + Ifges(5,5) * t238;
t262 = mrSges(6,1) * t141 - mrSges(6,2) * t142 + Ifges(6,5) * t148 + Ifges(6,6) * t147 + Ifges(6,3) * t158 + t153 * t183 - t154 * t182;
t117 = -mrSges(5,1) * t167 + mrSges(5,3) * t151 + Ifges(5,4) * t161 + Ifges(5,2) * t160 + Ifges(5,6) * t237 - pkin(4) * t128 - t172 * t197 + t174 * t238 - t262;
t191 = Ifges(4,5) * t217 + Ifges(4,6) * t216 + Ifges(4,3) * t243;
t193 = Ifges(4,1) * t217 + Ifges(4,4) * t216 + Ifges(4,5) * t243;
t111 = -mrSges(4,1) * t202 + mrSges(4,3) * t180 + Ifges(4,4) * t188 + Ifges(4,2) * t187 + Ifges(4,6) * t242 - pkin(3) * t265 + t248 * t116 + t253 * t117 - t217 * t191 + t243 * t193;
t192 = Ifges(4,4) * t217 + Ifges(4,2) * t216 + Ifges(4,6) * t243;
t112 = mrSges(4,2) * t202 - mrSges(4,3) * t179 + Ifges(4,1) * t188 + Ifges(4,4) * t187 + Ifges(4,5) * t242 + t116 * t253 - t117 * t248 + t191 * t216 - t192 * t243;
t213 = Ifges(3,3) * qJD(2) + (-Ifges(3,5) * t250 - Ifges(3,6) * t255) * qJD(1);
t215 = Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t250 - Ifges(3,4) * t255) * qJD(1);
t261 = m(4) * t202 - mrSges(4,1) * t187 + t188 * mrSges(4,2) - t206 * t216 + t217 * t207 + t265;
t108 = -mrSges(3,1) * t224 + mrSges(3,3) * t211 + Ifges(3,4) * t227 + Ifges(3,2) * t228 + Ifges(3,6) * qJDD(2) - pkin(2) * t261 + qJD(2) * t215 + t254 * t111 + t249 * t112 + t213 * t271;
t214 = Ifges(3,6) * qJD(2) + (-Ifges(3,4) * t250 - Ifges(3,2) * t255) * qJD(1);
t109 = mrSges(3,2) * t224 - mrSges(3,3) * t210 + Ifges(3,1) * t227 + Ifges(3,4) * t228 + Ifges(3,5) * qJDD(2) - qJD(2) * t214 - t111 * t249 + t112 * t254 - t213 * t270;
t230 = qJD(2) * mrSges(3,1) + mrSges(3,3) * t271;
t231 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t270;
t259 = m(3) * t224 - mrSges(3,1) * t228 + t227 * mrSges(3,2) - t230 * t271 + t231 * t270 + t261;
t264 = mrSges(2,1) * t232 - mrSges(2,2) * t233 + Ifges(2,3) * qJDD(1) + pkin(1) * t259 - t108 * t255 - t109 * t250;
t263 = mrSges(5,1) * t150 - mrSges(5,2) * t151 + Ifges(5,5) * t161 + Ifges(5,6) * t160 + Ifges(5,3) * t237 + pkin(4) * t266 + pkin(6) * t267 + t252 * t131 + t247 * t132 + t197 * t173 - t196 * t174;
t260 = mrSges(4,1) * t179 - mrSges(4,2) * t180 + Ifges(4,5) * t188 + Ifges(4,6) * t187 + Ifges(4,3) * t242 + pkin(3) * t272 + t217 * t192 - t216 * t193 + t263;
t258 = mrSges(3,1) * t210 - mrSges(3,2) * t211 + Ifges(3,5) * t227 + Ifges(3,6) * t228 + Ifges(3,3) * qJDD(2) + pkin(2) * t273 - t214 * t271 + t215 * t270 + t260;
t225 = (mrSges(3,1) * t255 - mrSges(3,2) * t250) * qJD(1);
t123 = m(2) * t232 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t257 + t259;
t114 = m(3) * t211 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t228 - qJD(2) * t230 - t120 * t249 + t121 * t254 - t225 * t270;
t113 = m(3) * t210 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t227 + qJD(2) * t231 + t225 * t271 + t273;
t110 = m(2) * t233 - mrSges(2,1) * t257 - qJDD(1) * mrSges(2,2) - t113 * t250 + t255 * t114;
t107 = t258 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - pkin(1) * (-t113 * t255 - t114 * t250) + t257 * Ifges(2,5) + mrSges(2,3) * t233;
t106 = -mrSges(2,2) * g(3) - mrSges(2,3) * t232 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t257 - t108 * t250 + t109 * t255;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t256 * t106 - t251 * t107 - pkin(5) * (t110 * t251 + t123 * t256), t106, t109, t112, t116, t132; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t251 * t106 + t256 * t107 + pkin(5) * (t110 * t256 - t123 * t251), t107, t108, t111, t117, t131; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t264, t264, t258, t260, t263, t262;];
m_new  = t1;
