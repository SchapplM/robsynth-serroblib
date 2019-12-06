% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:44
% EndTime: 2019-12-05 17:24:09
% DurationCPUTime: 15.76s
% Computational Cost: add. (281455->277), mult. (577968->376), div. (0->0), fcn. (452474->14), ass. (0->129)
t239 = sin(pkin(6));
t246 = sin(qJ(3));
t250 = cos(qJ(3));
t263 = qJD(2) * qJD(3);
t223 = (-qJDD(2) * t250 + t246 * t263) * t239;
t238 = sin(pkin(11));
t241 = cos(pkin(11));
t230 = t238 * g(1) - t241 * g(2);
t231 = -t241 * g(1) - t238 * g(2);
t237 = -g(3) + qJDD(1);
t247 = sin(qJ(2));
t243 = cos(pkin(5));
t251 = cos(qJ(2));
t265 = t243 * t251;
t240 = sin(pkin(5));
t268 = t240 * t251;
t201 = t230 * t265 - t247 * t231 + t237 * t268;
t252 = qJD(2) ^ 2;
t272 = pkin(8) * t239;
t197 = qJDD(2) * pkin(2) + t252 * t272 + t201;
t266 = t243 * t247;
t269 = t240 * t247;
t202 = t230 * t266 + t251 * t231 + t237 * t269;
t198 = -t252 * pkin(2) + qJDD(2) * t272 + t202;
t216 = -t240 * t230 + t243 * t237;
t242 = cos(pkin(6));
t173 = -t246 * t198 + (t197 * t242 + t216 * t239) * t250;
t267 = t242 * t246;
t270 = t239 * t246;
t174 = t197 * t267 + t250 * t198 + t216 * t270;
t236 = t242 * qJD(2) + qJD(3);
t264 = qJD(2) * t239;
t262 = t246 * t264;
t218 = t236 * mrSges(4,1) - mrSges(4,3) * t262;
t220 = (-mrSges(4,1) * t250 + mrSges(4,2) * t246) * t264;
t235 = t242 * qJDD(2) + qJDD(3);
t221 = (-pkin(3) * t250 - pkin(9) * t246) * t264;
t234 = t236 ^ 2;
t261 = t250 * t264;
t170 = -t234 * pkin(3) + t235 * pkin(9) + t221 * t261 + t174;
t212 = t242 * t216;
t222 = (qJDD(2) * t246 + t250 * t263) * t239;
t172 = t223 * pkin(3) - t222 * pkin(9) + t212 + (-t197 + (pkin(3) * t246 - pkin(9) * t250) * t236 * qJD(2)) * t239;
t245 = sin(qJ(4));
t249 = cos(qJ(4));
t167 = t249 * t170 + t245 * t172;
t214 = t249 * t236 - t245 * t262;
t215 = t245 * t236 + t249 * t262;
t200 = -t214 * pkin(4) - t215 * pkin(10);
t217 = qJDD(4) + t223;
t229 = qJD(4) - t261;
t228 = t229 ^ 2;
t164 = -t228 * pkin(4) + t217 * pkin(10) + t214 * t200 + t167;
t169 = -t235 * pkin(3) - t234 * pkin(9) + t221 * t262 - t173;
t192 = -t215 * qJD(4) - t245 * t222 + t249 * t235;
t193 = t214 * qJD(4) + t249 * t222 + t245 * t235;
t165 = (-t214 * t229 - t193) * pkin(10) + (t215 * t229 - t192) * pkin(4) + t169;
t244 = sin(qJ(5));
t248 = cos(qJ(5));
t161 = -t244 * t164 + t248 * t165;
t203 = -t244 * t215 + t248 * t229;
t177 = t203 * qJD(5) + t248 * t193 + t244 * t217;
t204 = t248 * t215 + t244 * t229;
t182 = -t203 * mrSges(6,1) + t204 * mrSges(6,2);
t213 = qJD(5) - t214;
t184 = -t213 * mrSges(6,2) + t203 * mrSges(6,3);
t190 = qJDD(5) - t192;
t158 = m(6) * t161 + t190 * mrSges(6,1) - t177 * mrSges(6,3) - t204 * t182 + t213 * t184;
t162 = t248 * t164 + t244 * t165;
t176 = -t204 * qJD(5) - t244 * t193 + t248 * t217;
t185 = t213 * mrSges(6,1) - t204 * mrSges(6,3);
t159 = m(6) * t162 - t190 * mrSges(6,2) + t176 * mrSges(6,3) + t203 * t182 - t213 * t185;
t152 = -t244 * t158 + t248 * t159;
t199 = -t214 * mrSges(5,1) + t215 * mrSges(5,2);
t206 = t229 * mrSges(5,1) - t215 * mrSges(5,3);
t150 = m(5) * t167 - t217 * mrSges(5,2) + t192 * mrSges(5,3) + t214 * t199 - t229 * t206 + t152;
t166 = -t245 * t170 + t249 * t172;
t163 = -t217 * pkin(4) - t228 * pkin(10) + t215 * t200 - t166;
t160 = -m(6) * t163 + t176 * mrSges(6,1) - t177 * mrSges(6,2) + t203 * t184 - t204 * t185;
t205 = -t229 * mrSges(5,2) + t214 * mrSges(5,3);
t156 = m(5) * t166 + t217 * mrSges(5,1) - t193 * mrSges(5,3) - t215 * t199 + t229 * t205 + t160;
t260 = t249 * t150 - t245 * t156;
t141 = m(4) * t174 - t235 * mrSges(4,2) - t223 * mrSges(4,3) - t236 * t218 + t220 * t261 + t260;
t144 = t245 * t150 + t249 * t156;
t183 = -t239 * t197 + t212;
t219 = -t236 * mrSges(4,2) + mrSges(4,3) * t261;
t143 = m(4) * t183 + t223 * mrSges(4,1) + t222 * mrSges(4,2) + (t218 * t246 - t219 * t250) * t264 + t144;
t151 = t248 * t158 + t244 * t159;
t255 = -m(5) * t169 + t192 * mrSges(5,1) - t193 * mrSges(5,2) + t214 * t205 - t215 * t206 - t151;
t147 = m(4) * t173 + t235 * mrSges(4,1) - t222 * mrSges(4,3) + t236 * t219 - t220 * t262 + t255;
t271 = t147 * t250;
t131 = t141 * t267 - t239 * t143 + t242 * t271;
t128 = m(3) * t201 + qJDD(2) * mrSges(3,1) - t252 * mrSges(3,2) + t131;
t135 = t250 * t141 - t246 * t147;
t134 = m(3) * t202 - t252 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t135;
t123 = -t247 * t128 + t251 * t134;
t273 = pkin(7) * t123;
t130 = t141 * t270 + t242 * t143 + t239 * t271;
t129 = m(3) * t216 + t130;
t120 = t128 * t265 - t240 * t129 + t134 * t266;
t178 = Ifges(6,5) * t204 + Ifges(6,6) * t203 + Ifges(6,3) * t213;
t180 = Ifges(6,1) * t204 + Ifges(6,4) * t203 + Ifges(6,5) * t213;
t153 = -mrSges(6,1) * t163 + mrSges(6,3) * t162 + Ifges(6,4) * t177 + Ifges(6,2) * t176 + Ifges(6,6) * t190 - t204 * t178 + t213 * t180;
t179 = Ifges(6,4) * t204 + Ifges(6,2) * t203 + Ifges(6,6) * t213;
t154 = mrSges(6,2) * t163 - mrSges(6,3) * t161 + Ifges(6,1) * t177 + Ifges(6,4) * t176 + Ifges(6,5) * t190 + t203 * t178 - t213 * t179;
t186 = Ifges(5,5) * t215 + Ifges(5,6) * t214 + Ifges(5,3) * t229;
t187 = Ifges(5,4) * t215 + Ifges(5,2) * t214 + Ifges(5,6) * t229;
t136 = mrSges(5,2) * t169 - mrSges(5,3) * t166 + Ifges(5,1) * t193 + Ifges(5,4) * t192 + Ifges(5,5) * t217 - pkin(10) * t151 - t244 * t153 + t248 * t154 + t214 * t186 - t229 * t187;
t188 = Ifges(5,1) * t215 + Ifges(5,4) * t214 + Ifges(5,5) * t229;
t254 = mrSges(6,1) * t161 - mrSges(6,2) * t162 + Ifges(6,5) * t177 + Ifges(6,6) * t176 + Ifges(6,3) * t190 + t204 * t179 - t203 * t180;
t137 = -mrSges(5,1) * t169 + mrSges(5,3) * t167 + Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * t217 - pkin(4) * t151 - t215 * t186 + t229 * t188 - t254;
t208 = Ifges(4,3) * t236 + (Ifges(4,5) * t246 + Ifges(4,6) * t250) * t264;
t209 = Ifges(4,6) * t236 + (Ifges(4,4) * t246 + Ifges(4,2) * t250) * t264;
t125 = mrSges(4,2) * t183 - mrSges(4,3) * t173 + Ifges(4,1) * t222 - Ifges(4,4) * t223 + Ifges(4,5) * t235 - pkin(9) * t144 + t249 * t136 - t245 * t137 + t208 * t261 - t236 * t209;
t210 = Ifges(4,5) * t236 + (Ifges(4,1) * t246 + Ifges(4,4) * t250) * t264;
t253 = mrSges(5,1) * t166 - mrSges(5,2) * t167 + Ifges(5,5) * t193 + Ifges(5,6) * t192 + Ifges(5,3) * t217 + pkin(4) * t160 + pkin(10) * t152 + t248 * t153 + t244 * t154 + t215 * t187 - t214 * t188;
t126 = -mrSges(4,1) * t183 + mrSges(4,3) * t174 + Ifges(4,4) * t222 - Ifges(4,2) * t223 + Ifges(4,6) * t235 - pkin(3) * t144 - t208 * t262 + t236 * t210 - t253;
t257 = pkin(8) * t135 + t125 * t246 + t126 * t250;
t124 = Ifges(4,5) * t222 - Ifges(4,6) * t223 + Ifges(4,3) * t235 + mrSges(4,1) * t173 - mrSges(4,2) * t174 + t245 * t136 + t249 * t137 + pkin(3) * t255 + pkin(9) * t260 + (t209 * t246 - t210 * t250) * t264;
t112 = mrSges(3,1) * t201 - mrSges(3,2) * t202 + Ifges(3,3) * qJDD(2) + pkin(2) * t131 + t242 * t124 + t257 * t239;
t114 = -mrSges(3,1) * t216 + mrSges(3,3) * t202 + t252 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t130 - t239 * t124 + t257 * t242;
t116 = mrSges(3,2) * t216 - mrSges(3,3) * t201 + Ifges(3,5) * qJDD(2) - t252 * Ifges(3,6) + t250 * t125 - t246 * t126 + (-t130 * t239 - t131 * t242) * pkin(8);
t256 = mrSges(2,1) * t230 - mrSges(2,2) * t231 + pkin(1) * t120 + t243 * t112 + t114 * t268 + t116 * t269 + t240 * t273;
t121 = m(2) * t231 + t123;
t119 = t243 * t129 + (t128 * t251 + t134 * t247) * t240;
t117 = m(2) * t230 + t120;
t110 = mrSges(2,2) * t237 - mrSges(2,3) * t230 - t247 * t114 + t251 * t116 + (-t119 * t240 - t120 * t243) * pkin(7);
t109 = -mrSges(2,1) * t237 + mrSges(2,3) * t231 - pkin(1) * t119 - t240 * t112 + (t114 * t251 + t116 * t247 + t273) * t243;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t241 * t110 - t238 * t109 - qJ(1) * (t241 * t117 + t238 * t121), t110, t116, t125, t136, t154; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t238 * t110 + t241 * t109 + qJ(1) * (-t238 * t117 + t241 * t121), t109, t114, t126, t137, t153; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t256, t256, t112, t124, t253, t254;];
m_new = t1;
