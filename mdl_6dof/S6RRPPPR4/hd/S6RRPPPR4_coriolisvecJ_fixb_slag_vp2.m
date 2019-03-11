% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:18
% EndTime: 2019-03-09 08:17:31
% DurationCPUTime: 6.85s
% Computational Cost: add. (3867->523), mult. (9259->704), div. (0->0), fcn. (5333->6), ass. (0->246)
t324 = Ifges(5,1) + Ifges(6,1);
t320 = Ifges(6,4) + Ifges(5,5);
t195 = cos(qJ(2));
t296 = pkin(3) + pkin(7);
t323 = t296 * t195;
t322 = -mrSges(3,1) + mrSges(4,2);
t321 = -Ifges(5,4) + Ifges(6,5);
t190 = cos(pkin(9));
t192 = sin(qJ(6));
t189 = sin(pkin(9));
t194 = cos(qJ(6));
t267 = t189 * t194;
t139 = -t190 * t192 + t267;
t193 = sin(qJ(2));
t268 = t189 * t193;
t295 = pkin(4) + pkin(5);
t207 = -pkin(8) * t268 - t195 * t295;
t256 = qJD(1) * t193;
t178 = pkin(2) * t256;
t214 = -qJ(3) * t195 + qJ(4) * t193;
t118 = qJD(1) * t214 + t178;
t255 = qJD(1) * t195;
t179 = pkin(7) * t255;
t144 = pkin(3) * t255 + t179;
t61 = -t189 * t118 + t190 * t144;
t37 = qJD(1) * t207 - t61;
t240 = t190 * t256;
t62 = t190 * t118 + t189 * t144;
t54 = qJ(5) * t255 + t62;
t42 = -pkin(8) * t240 + t54;
t191 = -pkin(2) - qJ(4);
t283 = pkin(8) + t191;
t147 = t283 * t189;
t148 = t283 * t190;
t86 = -t147 * t192 - t148 * t194;
t319 = -qJD(4) * t139 + qJD(6) * t86 - t192 * t37 - t194 * t42;
t212 = t192 * t189 + t194 * t190;
t87 = t147 * t194 - t148 * t192;
t318 = qJD(4) * t212 - qJD(6) * t87 + t192 * t42 - t194 * t37;
t275 = Ifges(6,5) * t190;
t279 = Ifges(5,4) * t190;
t317 = t189 * t324 - t275 + t279;
t244 = Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t204 = t207 * qJD(2);
t117 = (qJD(1) * t323 - qJD(4)) * qJD(2);
t249 = qJD(1) * qJD(2);
t238 = t193 * t249;
t167 = pkin(2) * t238;
t251 = t193 * qJD(3);
t198 = qJD(2) * t214 - qJD(4) * t195 - t251;
t74 = qJD(1) * t198 + t167;
t32 = t190 * t117 - t189 * t74;
t14 = qJD(1) * t204 - t32;
t183 = t193 * qJD(5);
t237 = t195 * t249;
t33 = t189 * t117 + t190 * t74;
t20 = qJ(5) * t237 + qJD(1) * t183 + t33;
t232 = t190 * t238;
t15 = -pkin(8) * t232 + t20;
t252 = t190 * qJD(2);
t137 = -t189 * t255 + t252;
t236 = -t193 * qJ(3) - pkin(1);
t135 = t191 * t195 + t236;
t107 = t135 * qJD(1);
t177 = pkin(3) * t256;
t176 = pkin(7) * t256;
t250 = qJD(3) + t176;
t111 = qJD(2) * t191 + t177 + t250;
t45 = -t189 * t107 + t111 * t190;
t226 = qJD(5) - t45;
t241 = t295 * t193;
t19 = -pkin(8) * t137 - qJD(1) * t241 + t226;
t136 = t189 * qJD(2) + t190 * t255;
t46 = t190 * t107 + t189 * t111;
t39 = qJ(5) * t256 + t46;
t24 = pkin(8) * t136 + t39;
t5 = t19 * t194 - t192 * t24;
t1 = qJD(6) * t5 + t14 * t192 + t15 * t194;
t6 = t19 * t192 + t194 * t24;
t2 = -qJD(6) * t6 + t14 * t194 - t15 * t192;
t209 = t212 * t193;
t106 = qJD(1) * t209;
t126 = t212 * qJD(6);
t258 = -t126 + t106;
t105 = t192 * t240 - t256 * t267;
t125 = t139 * qJD(6);
t259 = -t125 - t105;
t316 = t1 * t139 - t2 * t212 + t258 * t6 + t259 * t5;
t303 = -Ifges(7,6) / 0.2e1;
t254 = qJD(2) * t193;
t206 = t139 * t254;
t76 = t136 * t194 - t137 * t192;
t34 = qJD(1) * t206 + qJD(6) * t76;
t302 = t34 / 0.2e1;
t205 = qJD(2) * t209;
t213 = t136 * t192 + t137 * t194;
t35 = -qJD(1) * t205 - qJD(6) * t213;
t301 = t35 / 0.2e1;
t315 = t137 / 0.2e1;
t314 = -t249 / 0.2e1;
t310 = -qJD(1) / 0.2e1;
t309 = -qJD(2) / 0.2e1;
t308 = qJD(2) / 0.2e1;
t113 = t139 * t195;
t168 = qJD(6) - t256;
t307 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t34 - Ifges(7,6) * t35;
t269 = qJ(5) * t189;
t211 = t190 * t295 + t269;
t152 = -pkin(2) * t195 + t236;
t134 = t152 * qJD(1);
t151 = -qJD(2) * pkin(2) + t250;
t175 = Ifges(3,4) * t255;
t246 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t274 = Ifges(4,6) * t195;
t38 = -pkin(4) * t256 + t226;
t306 = -0.2e1 * t244 * t136 + t246 * t137 + t151 * mrSges(4,1) + t39 * mrSges(6,3) + t45 * mrSges(5,1) + t6 * mrSges(7,2) + Ifges(3,5) * t308 + t175 / 0.2e1 + Ifges(4,4) * t309 + (-Ifges(4,2) * t193 - t274) * t310 - t168 * Ifges(7,3) - t213 * Ifges(7,5) - t76 * Ifges(7,6) - t134 * mrSges(4,3) - t38 * mrSges(6,1) - t46 * mrSges(5,2) - t5 * mrSges(7,1) + t320 * t315 + (Ifges(3,1) + Ifges(5,3) + Ifges(6,2)) * t256 / 0.2e1;
t188 = qJD(2) * qJ(3);
t124 = qJD(4) + t188 + t144;
t157 = -t179 - t188;
t159 = -mrSges(4,1) * t255 - qJD(2) * mrSges(4,3);
t216 = t189 * t45 - t190 * t46;
t217 = t189 * t38 + t190 * t39;
t276 = Ifges(6,5) * t189;
t221 = -Ifges(6,3) * t190 + t276;
t280 = Ifges(5,4) * t189;
t222 = Ifges(5,2) * t190 + t280;
t281 = mrSges(6,3) * t189;
t225 = -mrSges(6,1) * t190 - t281;
t287 = t190 / 0.2e1;
t288 = -t190 / 0.2e1;
t289 = t189 / 0.2e1;
t208 = qJ(5) * t137 - t124;
t50 = pkin(4) * t136 - t208;
t305 = t217 * mrSges(6,2) - t216 * mrSges(5,3) + (m(4) * t157 + qJD(2) * mrSges(3,2) - mrSges(3,3) * t255 + t159) * pkin(7) + t124 * (-mrSges(5,1) * t190 + mrSges(5,2) * t189) + t157 * mrSges(4,1) + t50 * t225 + Ifges(3,6) * t309 + (Ifges(3,4) * t193 + t195 * Ifges(3,2)) * t310 + Ifges(4,5) * t308 + (-Ifges(4,6) * t193 - t195 * Ifges(4,3)) * qJD(1) / 0.2e1 - t134 * mrSges(4,2) + (t137 * Ifges(6,5) + Ifges(6,6) * t256) * t288 + (t137 * Ifges(5,4) + Ifges(5,6) * t256) * t287 + t317 * t315 + (t137 * t324 + t320 * t256) * t289 + (Ifges(6,3) * t288 - Ifges(5,2) * t287 + t321 * t289 - t222 / 0.2e1 + t221 / 0.2e1) * t136;
t304 = Ifges(7,4) * t302 + Ifges(7,2) * t301 + t237 * t303;
t300 = -t76 / 0.2e1;
t299 = t76 / 0.2e1;
t298 = -t213 / 0.2e1;
t297 = t213 / 0.2e1;
t294 = pkin(1) * mrSges(3,1);
t293 = pkin(1) * mrSges(3,2);
t291 = -t168 / 0.2e1;
t290 = t168 / 0.2e1;
t286 = Ifges(7,4) * t213;
t21 = -mrSges(7,1) * t76 + mrSges(7,2) * t213;
t83 = mrSges(6,1) * t136 - mrSges(6,3) * t137;
t282 = -t21 + t83;
t278 = Ifges(6,4) * t189;
t277 = Ifges(5,5) * t189;
t273 = Ifges(5,6) * t190;
t272 = Ifges(6,6) * t190;
t84 = mrSges(5,1) * t136 + mrSges(5,2) * t137;
t270 = -t159 + t84;
t146 = qJD(2) * t323;
t181 = pkin(2) * t254;
t93 = t181 + t198;
t49 = t189 * t146 + t190 * t93;
t266 = t189 * t195;
t264 = t190 * t193;
t101 = -mrSges(6,2) * t136 + mrSges(6,3) * t256;
t102 = -mrSges(5,2) * t256 - mrSges(5,3) * t136;
t263 = -t102 - t101;
t103 = mrSges(5,1) * t256 - mrSges(5,3) * t137;
t104 = -mrSges(6,1) * t256 + mrSges(6,2) * t137;
t262 = t103 - t104;
t161 = t296 * t193;
t80 = t190 * t135 + t189 * t161;
t119 = (mrSges(5,1) * t195 - mrSges(5,3) * t268) * t249;
t233 = t189 * t238;
t154 = mrSges(6,2) * t233;
t120 = -mrSges(6,1) * t237 + t154;
t261 = t119 - t120;
t121 = (-mrSges(5,2) * t195 + mrSges(5,3) * t264) * t249;
t122 = (mrSges(6,2) * t264 + mrSges(6,3) * t195) * t249;
t260 = t121 + t122;
t257 = t176 + t177;
t253 = qJD(2) * t195;
t248 = 0.3e1 / 0.2e1 * Ifges(3,4) + 0.3e1 / 0.2e1 * Ifges(4,6);
t247 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t245 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t243 = m(4) * pkin(7) + mrSges(4,1);
t71 = t193 * qJ(5) + t80;
t239 = qJD(5) * t266;
t48 = t190 * t146 - t189 * t93;
t235 = qJ(5) * t190 - qJ(3);
t79 = -t189 * t135 + t190 * t161;
t187 = qJD(2) * qJD(3);
t234 = -qJD(5) * t137 + t187;
t36 = qJ(5) * t253 + t183 + t49;
t145 = t296 * t254;
t230 = m(4) * t151 + (mrSges(4,1) + mrSges(3,3)) * t256 + t322 * qJD(2);
t9 = -t35 * mrSges(7,1) + t34 * mrSges(7,2);
t220 = pkin(4) * t190 + t269;
t25 = -pkin(4) * t237 - t32;
t219 = t189 * t20 - t190 * t25;
t218 = t189 * t33 + t190 * t32;
t47 = pkin(8) * t266 - t241 - t79;
t55 = pkin(8) * t190 * t195 + t71;
t12 = -t192 * t55 + t194 * t47;
t13 = t192 * t47 + t194 * t55;
t52 = -mrSges(7,2) * t168 + mrSges(7,3) * t76;
t53 = mrSges(7,1) * t168 - mrSges(7,3) * t213;
t215 = -t192 * t53 + t194 * t52;
t210 = -qJ(3) * t253 - t251;
t112 = t212 * t195;
t203 = (-t220 - t296) * t254;
t202 = (t211 + t296) * t254;
t201 = -t278 / 0.2e1 - t277 / 0.2e1 - t273 / 0.2e1 + t272 / 0.2e1;
t116 = -qJD(1) * t145 + t187;
t40 = qJD(1) * t203 + t234;
t200 = t116 * mrSges(5,1) + t40 * mrSges(6,1) + (Ifges(6,6) * t195 + t193 * t221) * t249 / 0.2e1 + (Ifges(5,6) * t195 + t193 * t222) * t314 - t33 * mrSges(5,3) - t20 * mrSges(6,2);
t199 = -t116 * mrSges(5,2) - t25 * mrSges(6,2) + t32 * mrSges(5,3) + t40 * mrSges(6,3) + (t193 * t317 + t195 * t320) * t314;
t166 = -t190 * qJD(5) + qJD(3);
t155 = mrSges(5,2) * t233;
t150 = pkin(7) * t238 - t187;
t149 = pkin(4) * t189 - t235;
t141 = (mrSges(4,2) * t195 - t193 * mrSges(4,3)) * qJD(1);
t127 = -t189 * t295 + t235;
t123 = t181 + t210;
t110 = -mrSges(5,1) * t232 + t155;
t109 = t225 * t238;
t108 = qJD(1) * t210 + t167;
t96 = t195 * t220 + t323;
t88 = -t220 * t256 - t257;
t85 = -t195 * t211 - t323;
t73 = Ifges(7,4) * t76;
t72 = -t193 * pkin(4) - t79;
t64 = t203 + t239;
t63 = t211 * t256 + t257;
t58 = qJD(6) * t112 + t206;
t57 = qJD(6) * t113 - t205;
t56 = -pkin(4) * t255 - t61;
t51 = t202 - t239;
t41 = -pkin(4) * t253 - t48;
t31 = -t136 * t295 + t208;
t28 = qJD(1) * t202 - t234;
t27 = -pkin(8) * t193 * t252 + t36;
t26 = t204 - t48;
t23 = mrSges(7,2) * t237 + t35 * mrSges(7,3);
t22 = -mrSges(7,1) * t237 - t34 * mrSges(7,3);
t18 = Ifges(7,1) * t213 + Ifges(7,5) * t168 + t73;
t17 = Ifges(7,2) * t76 + Ifges(7,6) * t168 + t286;
t8 = Ifges(7,1) * t34 + Ifges(7,4) * t35 - Ifges(7,5) * t237;
t4 = -qJD(6) * t13 - t192 * t27 + t194 * t26;
t3 = qJD(6) * t12 + t192 * t26 + t194 * t27;
t7 = [(t108 * mrSges(4,2) - t243 * t150 + t200 * t190 + t199 * t189 + (t247 * qJD(2) + t230 * pkin(7) + (-t152 * mrSges(4,3) - 0.2e1 * t293 + Ifges(7,5) * t113 / 0.2e1 + t112 * t303 + (t201 + t248) * t195) * qJD(1) + t306) * qJD(2)) * t195 + (-t108 * mrSges(4,3) + t32 * mrSges(5,1) - t25 * mrSges(6,1) - t33 * mrSges(5,2) + t20 * mrSges(6,3) + (t245 * qJD(2) + ((0.3e1 / 0.2e1 * t278 - 0.3e1 / 0.2e1 * t272 + 0.3e1 / 0.2e1 * t277 + 0.3e1 / 0.2e1 * t273 - t248) * t193 - t152 * mrSges(4,2) - 0.2e1 * t294 + (0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + Ifges(7,3) + t243 * pkin(7) + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1) * t190 ^ 2 + ((-Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1) * t189 + t321 * t190) * t189) * t195) * qJD(1) + t305) * qJD(2) + t307) * t193 + t28 * (-mrSges(7,1) * t112 - mrSges(7,2) * t113) + (t1 * t112 + t113 * t2 - t5 * t58 + t57 * t6) * mrSges(7,3) + (-Ifges(7,4) * t113 + Ifges(7,2) * t112) * t301 + (-Ifges(7,1) * t113 + Ifges(7,4) * t112) * t302 + t123 * t141 - t145 * t84 + t80 * t121 + t71 * t122 - t113 * t8 / 0.2e1 + t79 * t119 + t72 * t120 + t48 * t103 + t41 * t104 + t96 * t109 + t36 * t101 + t49 * t102 + t64 * t83 + t85 * t9 + t57 * t17 / 0.2e1 + t31 * (-mrSges(7,1) * t57 + mrSges(7,2) * t58) + t58 * t18 / 0.2e1 + t51 * t21 + t3 * t52 + t4 * t53 + t12 * t22 + t13 * t23 + m(7) * (t1 * t13 + t12 * t2 + t28 * t85 + t3 * t6 + t31 * t51 + t4 * t5) + m(6) * (t20 * t71 + t25 * t72 + t36 * t39 + t38 * t41 + t40 * t96 + t50 * t64) + m(5) * (t116 * t323 - t124 * t145 + t32 * t79 + t33 * t80 + t45 * t48 + t46 * t49) + t323 * t110 + (Ifges(7,4) * t58 + Ifges(7,2) * t57) * t299 + t112 * t304 + (Ifges(7,5) * t58 + Ifges(7,6) * t57) * t290 + (Ifges(7,1) * t58 + Ifges(7,4) * t57) * t297 + m(4) * (t108 * t152 + t134 * t123); (-m(4) * t134 - t141) * (-qJ(3) * t255 + t178) + t28 * (-mrSges(7,1) * t139 + mrSges(7,2) * t212) + (Ifges(7,4) * t212 + Ifges(7,2) * t139) * t301 + (Ifges(7,1) * t212 + Ifges(7,4) * t139) * t302 + t212 * t8 / 0.2e1 + t257 * t84 + (-t45 * t61 - t46 * t62 + qJ(3) * t116 + t218 * t191 + (-t189 * t46 - t190 * t45) * qJD(4) + (t257 + qJD(3)) * t124) * m(5) + (-t126 / 0.2e1 + t106 / 0.2e1) * t17 + (-Ifges(7,4) * t105 - Ifges(7,2) * t106) * t300 + (-Ifges(7,5) * t105 - Ifges(7,6) * t106) * t291 + (-Ifges(7,1) * t105 - Ifges(7,4) * t106) * t298 + ((-t175 / 0.2e1 + ((-m(4) * pkin(2) + t322) * qJD(2) - t230) * pkin(7) + (-Ifges(7,5) * t212 / 0.2e1 + t139 * t303 - pkin(2) * mrSges(4,1) + t246 * t190 - t244 * t189 + t247) * qJD(2) + (t293 - t274 / 0.2e1) * qJD(1) - t306) * t195 + (((Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1 + t201) * t193 + t294) * qJD(1) + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1) * t255 + ((Ifges(6,3) * t189 + t275) * t288 + (-Ifges(5,2) * t189 + t279) * t287 - qJ(3) * mrSges(4,1) + pkin(7) * mrSges(3,2) + t245 + (t190 * t324 + t276 - t280) * t289) * qJD(2) - t305) * t193) * qJD(1) + t318 * t53 + (t1 * t87 + t127 * t28 + t2 * t86 + t319 * t6 + t318 * t5 + (-t166 - t63) * t31) * m(7) + t319 * t52 + (Ifges(7,4) * t125 - Ifges(7,2) * t126) * t299 + (Ifges(7,5) * t125 - Ifges(7,6) * t126) * t290 + (Ifges(7,1) * t125 - Ifges(7,4) * t126) * t297 + t149 * t109 - t150 * mrSges(4,3) + t127 * t9 - t61 * t103 - t56 * t104 + qJ(3) * t110 - t88 * t83 - t54 * t101 - t62 * t102 + t86 * t22 + t87 * t23 - t63 * t21 + t316 * mrSges(7,3) + (-t38 * t56 - t39 * t54 + t149 * t40 + t219 * t191 + (-t189 * t39 + t190 * t38) * qJD(4) + (-t88 + t166) * t50) * m(6) + (t125 / 0.2e1 + t105 / 0.2e1) * t18 + t282 * t166 + t270 * qJD(3) + (-mrSges(7,1) * t258 - mrSges(7,2) * t259) * t31 + (-qJD(4) * t262 + t191 * t261 - t199) * t190 + (qJD(4) * t263 + t191 * t260 + t200) * t189 + t139 * t304 + m(4) * (-t150 * qJ(3) - t157 * qJD(3)); -t212 * t22 + t139 * t23 + t259 * t53 + t258 * t52 + t261 * t190 + t260 * t189 + (-t270 - t282) * qJD(2) + (t243 * t253 + (-t189 * t262 - t190 * t263 + t141) * t193) * qJD(1) - m(4) * (-qJD(2) * t157 - t134 * t256) + (qJD(2) * t31 + t316) * m(7) + (-qJD(2) * t50 + t217 * t256 + t219) * m(6) + (-qJD(2) * t124 - t216 * t256 + t218) * m(5); t76 * t52 - t213 * t53 + t155 + t262 * t137 - t263 * t136 + (-t281 + (-mrSges(5,1) - mrSges(6,1)) * t190) * t238 - t9 + (-t213 * t5 + t6 * t76 - t28) * m(7) + (t136 * t39 - t137 * t38 + t40) * m(6) + (t136 * t46 + t137 * t45 + t116) * m(5); t192 * t23 + t194 * t22 + t154 + t282 * t137 + t215 * qJD(6) + (-mrSges(6,1) * t253 + (-t101 - t215) * t193) * qJD(1) + (t1 * t192 - t137 * t31 + t194 * t2 + t168 * (-t192 * t5 + t194 * t6)) * m(7) + (t137 * t50 - t256 * t39 + t25) * m(6); -Ifges(7,3) * t237 - t31 * (mrSges(7,1) * t213 + mrSges(7,2) * t76) + (Ifges(7,1) * t76 - t286) * t298 + t17 * t297 + (Ifges(7,5) * t76 - Ifges(7,6) * t213) * t291 - t5 * t52 + t6 * t53 + (t213 * t6 + t5 * t76) * mrSges(7,3) + (-Ifges(7,2) * t213 + t18 + t73) * t300 - t307;];
tauc  = t7(:);
