% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:03:01
% EndTime: 2019-03-09 06:03:29
% DurationCPUTime: 14.77s
% Computational Cost: add. (7491->629), mult. (17702->860), div. (0->0), fcn. (11310->8), ass. (0->246)
t334 = Ifges(6,1) + Ifges(7,1);
t333 = Ifges(7,4) + Ifges(6,5);
t332 = Ifges(7,5) - Ifges(6,4);
t202 = sin(qJ(5));
t203 = sin(qJ(4));
t205 = cos(qJ(5));
t206 = cos(qJ(4));
t166 = t202 * t203 - t205 * t206;
t325 = qJD(4) + qJD(5);
t115 = t325 * t166;
t207 = cos(qJ(3));
t212 = t166 * t207;
t139 = qJD(1) * t212;
t343 = t115 - t139;
t167 = t202 * t206 + t203 * t205;
t116 = t325 * t167;
t267 = qJD(1) * t207;
t138 = t167 * t267;
t271 = t116 - t138;
t331 = Ifges(7,2) + Ifges(6,3);
t330 = -Ifges(6,6) + Ifges(7,6);
t266 = qJD(3) * t203;
t204 = sin(qJ(3));
t268 = qJD(1) * t204;
t165 = t206 * t268 + t266;
t264 = qJD(3) * t206;
t214 = t203 * t268 - t264;
t110 = t202 * t165 + t205 * t214;
t105 = Ifges(6,4) * t110;
t191 = qJD(4) - t267;
t183 = qJD(5) + t191;
t210 = t205 * t165 - t202 * t214;
t285 = Ifges(7,5) * t110;
t329 = t333 * t183 + t334 * t210 - t105 + t285;
t342 = t329 / 0.2e1;
t256 = qJD(1) * qJD(3);
t240 = t204 * t256;
t255 = qJD(3) * qJD(4);
t261 = qJD(4) * t204;
t263 = qJD(3) * t207;
t327 = -t203 * t261 + t206 * t263;
t127 = qJD(1) * t327 + t206 * t255;
t244 = t203 * t263;
t260 = qJD(4) * t206;
t338 = t204 * t260 + t244;
t128 = -qJD(1) * t338 - t203 * t255;
t45 = -qJD(5) * t110 + t205 * t127 + t202 * t128;
t46 = qJD(5) * t210 + t202 * t127 - t205 * t128;
t341 = t333 * t240 + t332 * t46 + t334 * t45;
t193 = sin(pkin(10)) * pkin(1) + pkin(7);
t176 = t193 * qJD(1);
t142 = t204 * qJD(2) + t207 * t176;
t247 = t203 * t267;
t120 = pkin(4) * t247 + t142;
t262 = qJD(4) * t203;
t340 = pkin(4) * t262 + t271 * pkin(5) + t343 * qJ(6) - qJD(6) * t167 - t120;
t318 = t46 / 0.2e1;
t339 = Ifges(7,3) * t318;
t66 = pkin(5) * t210 + qJ(6) * t110;
t320 = t45 / 0.2e1;
t309 = t127 / 0.2e1;
t308 = t128 / 0.2e1;
t336 = t240 / 0.2e1;
t145 = t167 * t204;
t134 = qJD(3) * pkin(8) + t142;
t249 = -cos(pkin(10)) * pkin(1) - pkin(2);
t160 = -pkin(3) * t207 - t204 * pkin(8) + t249;
t137 = t160 * qJD(1);
t85 = t206 * t134 + t203 * t137;
t71 = -pkin(9) * t214 + t85;
t282 = t202 * t71;
t84 = -t134 * t203 + t206 * t137;
t70 = -pkin(9) * t165 + t84;
t64 = pkin(4) * t191 + t70;
t21 = t205 * t64 - t282;
t328 = qJD(6) - t21;
t141 = t207 * qJD(2) - t204 * t176;
t273 = t206 * t207;
t169 = t193 * t273;
t114 = t203 * t160 + t169;
t326 = t331 * t240 + t330 * t46 + t333 * t45;
t135 = t141 * qJD(3);
t233 = pkin(3) * t204 - pkin(8) * t207;
t172 = t233 * qJD(3);
t159 = qJD(1) * t172;
t38 = -t134 * t262 + t206 * t135 + t137 * t260 + t203 * t159;
t39 = -qJD(4) * t85 - t135 * t203 + t206 * t159;
t225 = -t203 * t39 + t206 * t38;
t133 = -qJD(3) * pkin(3) - t141;
t224 = t203 * t85 + t206 * t84;
t284 = Ifges(5,6) * t203;
t228 = Ifges(5,5) * t206 - t284;
t288 = Ifges(5,4) * t203;
t232 = Ifges(5,1) * t206 - t288;
t300 = t206 / 0.2e1;
t301 = t191 / 0.2e1;
t305 = t165 / 0.2e1;
t289 = Ifges(5,4) * t165;
t98 = -Ifges(5,2) * t214 + Ifges(5,6) * t191 + t289;
t161 = Ifges(5,4) * t214;
t99 = t165 * Ifges(5,1) + t191 * Ifges(5,5) - t161;
t324 = -t224 * mrSges(5,3) + t133 * (mrSges(5,1) * t203 + mrSges(5,2) * t206) + t232 * t305 + t228 * t301 - t203 * t98 / 0.2e1 + t99 * t300;
t280 = t205 * t71;
t22 = t202 * t64 + t280;
t26 = pkin(4) * t240 - pkin(9) * t127 + t39;
t28 = pkin(9) * t128 + t38;
t6 = -qJD(5) * t22 - t202 * t28 + t205 * t26;
t275 = t203 * t204;
t102 = -pkin(9) * t275 + t114;
t148 = t206 * t160;
t274 = t204 * t206;
t276 = t193 * t203;
t95 = -pkin(9) * t274 + t148 + (-pkin(4) - t276) * t207;
t292 = t205 * t102 + t202 * t95;
t217 = pkin(4) * t204 - pkin(9) * t273;
t265 = qJD(3) * t204;
t270 = t206 * t172 + t265 * t276;
t54 = t217 * qJD(3) + (-t169 + (pkin(9) * t204 - t160) * t203) * qJD(4) + t270;
t72 = t160 * t260 + t203 * t172 + (-t204 * t264 - t207 * t262) * t193;
t62 = -pkin(9) * t338 + t72;
t11 = -qJD(5) * t292 - t202 * t62 + t205 * t54;
t103 = pkin(4) * t214 + t133;
t20 = qJ(6) * t183 + t22;
t31 = t110 * pkin(5) - qJ(6) * t210 + t103;
t104 = Ifges(7,5) * t210;
t56 = t183 * Ifges(7,6) + t110 * Ifges(7,3) + t104;
t286 = Ifges(6,4) * t210;
t59 = -t110 * Ifges(6,2) + t183 * Ifges(6,6) + t286;
t323 = t20 * mrSges(7,2) + t22 * mrSges(6,3) - t56 / 0.2e1 + t59 / 0.2e1 - t103 * mrSges(6,1) - t31 * mrSges(7,1);
t322 = Ifges(7,5) * t320 + Ifges(7,6) * t336 + t339;
t321 = -t45 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t318 - Ifges(6,6) * t240 / 0.2e1;
t319 = -t46 / 0.2e1;
t316 = Ifges(5,1) * t309 + Ifges(5,4) * t308 + Ifges(5,5) * t336;
t315 = -pkin(9) - pkin(8);
t314 = -t110 / 0.2e1;
t313 = t110 / 0.2e1;
t311 = -t210 / 0.2e1;
t310 = t210 / 0.2e1;
t306 = -t165 / 0.2e1;
t303 = -t183 / 0.2e1;
t302 = t183 / 0.2e1;
t298 = m(6) * t103;
t297 = qJD(3) / 0.2e1;
t34 = -mrSges(7,2) * t46 + mrSges(7,3) * t240;
t37 = -mrSges(6,2) * t240 - mrSges(6,3) * t46;
t296 = t34 + t37;
t35 = mrSges(6,1) * t240 - mrSges(6,3) * t45;
t36 = -mrSges(7,1) * t240 + t45 * mrSges(7,2);
t295 = t36 - t35;
t170 = t233 * qJD(1);
t100 = -t203 * t141 + t206 * t170;
t83 = qJD(1) * t217 + t100;
t101 = t206 * t141 + t203 * t170;
t89 = -pkin(9) * t247 + t101;
t33 = t202 * t83 + t205 * t89;
t90 = -mrSges(7,2) * t110 + mrSges(7,3) * t183;
t291 = mrSges(6,3) * t110;
t91 = -mrSges(6,2) * t183 - t291;
t294 = t90 + t91;
t290 = mrSges(6,3) * t210;
t92 = mrSges(6,1) * t183 - t290;
t93 = -mrSges(7,1) * t183 + mrSges(7,2) * t210;
t293 = -t93 + t92;
t287 = Ifges(5,4) * t206;
t283 = Ifges(5,3) * t191;
t278 = Ifges(4,6) * qJD(3);
t136 = t142 * qJD(3);
t277 = t136 * t204;
t252 = mrSges(4,3) * t268;
t269 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t214 - t165 * mrSges(5,2) - t252;
t149 = pkin(4) * t275 + t204 * t193;
t259 = qJD(5) * t202;
t258 = qJD(5) * t205;
t251 = mrSges(4,3) * t267;
t250 = Ifges(5,5) * t127 + Ifges(5,6) * t128 + Ifges(5,3) * t240;
t118 = pkin(4) * t338 + t193 * t263;
t197 = -pkin(4) * t206 - pkin(3);
t248 = qJD(4) * t315;
t238 = t263 / 0.2e1;
t237 = -t261 / 0.2e1;
t236 = t206 * t248;
t178 = t249 * qJD(1);
t231 = Ifges(5,1) * t203 + t287;
t230 = -Ifges(5,2) * t203 + t287;
t229 = Ifges(5,2) * t206 + t288;
t32 = -t202 * t89 + t205 * t83;
t223 = t203 * t84 - t206 * t85;
t50 = -t202 * t102 + t205 * t95;
t107 = mrSges(5,1) * t240 - mrSges(5,3) * t127;
t108 = -mrSges(5,2) * t240 + mrSges(5,3) * t128;
t221 = -t203 * t107 + t206 * t108;
t130 = -t191 * mrSges(5,2) - mrSges(5,3) * t214;
t131 = mrSges(5,1) * t191 - mrSges(5,3) * t165;
t220 = -t203 * t130 - t206 * t131;
t219 = t135 * t207 + t277;
t181 = t315 * t203;
t182 = t315 * t206;
t218 = t205 * t181 + t182 * t202;
t124 = t181 * t202 - t182 * t205;
t5 = t202 * t26 + t205 * t28 + t64 * t258 - t259 * t71;
t10 = -t102 * t259 + t202 * t54 + t205 * t62 + t95 * t258;
t216 = t203 * t230;
t215 = t206 * t230;
t87 = -t128 * pkin(4) + t136;
t2 = qJ(6) * t240 + qJD(6) * t183 + t5;
t3 = -pkin(5) * t240 - t6;
t209 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t326;
t198 = Ifges(4,4) * t267;
t196 = -pkin(4) * t205 - pkin(5);
t194 = pkin(4) * t202 + qJ(6);
t187 = pkin(4) * t258 + qJD(6);
t179 = -qJD(3) * mrSges(4,2) + t251;
t171 = t203 * t248;
t153 = Ifges(4,1) * t268 + Ifges(4,5) * qJD(3) + t198;
t152 = t278 + (Ifges(4,4) * t204 + t207 * Ifges(4,2)) * qJD(1);
t146 = t166 * t204;
t113 = -t207 * t276 + t148;
t106 = pkin(5) * t166 - qJ(6) * t167 + t197;
t97 = Ifges(5,5) * t165 - Ifges(5,6) * t214 + t283;
t82 = -mrSges(5,1) * t128 + mrSges(5,2) * t127;
t80 = pkin(5) * t145 + qJ(6) * t146 + t149;
t79 = qJD(5) * t124 + t171 * t202 - t205 * t236;
t78 = qJD(5) * t218 + t205 * t171 + t202 * t236;
t76 = t127 * Ifges(5,4) + t128 * Ifges(5,2) + Ifges(5,6) * t240;
t75 = -t259 * t275 + (t274 * t325 + t244) * t205 + t327 * t202;
t74 = -qJD(3) * t212 - t145 * t325;
t73 = -qJD(4) * t114 + t270;
t68 = mrSges(6,1) * t110 + mrSges(6,2) * t210;
t67 = mrSges(7,1) * t110 - mrSges(7,3) * t210;
t58 = Ifges(7,4) * t210 + t183 * Ifges(7,2) + t110 * Ifges(7,6);
t57 = Ifges(6,5) * t210 - t110 * Ifges(6,6) + t183 * Ifges(6,3);
t53 = pkin(4) * t165 + t66;
t48 = pkin(5) * t207 - t50;
t47 = -qJ(6) * t207 + t292;
t30 = -pkin(5) * t268 - t32;
t29 = qJ(6) * t268 + t33;
t24 = t205 * t70 - t282;
t23 = t202 * t70 + t280;
t19 = -pkin(5) * t183 + t328;
t18 = pkin(5) * t75 - qJ(6) * t74 + qJD(6) * t146 + t118;
t17 = mrSges(6,1) * t46 + mrSges(6,2) * t45;
t16 = mrSges(7,1) * t46 - mrSges(7,3) * t45;
t9 = -pkin(5) * t265 - t11;
t8 = t46 * pkin(5) - t45 * qJ(6) - qJD(6) * t210 + t87;
t7 = qJ(6) * t265 - qJD(6) * t207 + t10;
t1 = [-(t250 + t326) * t207 / 0.2e1 + m(7) * (t18 * t31 + t19 * t9 + t2 * t47 + t20 * t7 + t3 * t48 + t8 * t80) - t214 * (-t229 * t261 + (Ifges(5,6) * t204 + t207 * t230) * qJD(3)) / 0.2e1 + (-t224 * t263 + (qJD(4) * t223 - t38 * t203 - t39 * t206) * t204) * mrSges(5,3) + (t203 * t237 + t206 * t238) * t99 + (Ifges(6,4) * t74 + Ifges(6,6) * t265) * t314 + (t204 * t82 + (-t204 * t179 - t207 * t269) * qJD(3) + m(5) * (t133 * t263 + t277) + m(4) * (-t141 * t263 - t142 * t265 + t219)) * t193 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t207 * t240 + (-Ifges(6,2) * t314 + Ifges(7,3) * t313 + t302 * t330 + t310 * t332 - t323) * t75 + 0.2e1 * t178 * (mrSges(4,1) * t204 + mrSges(4,2) * t207) * qJD(3) + m(6) * (t10 * t22 + t103 * t118 + t11 * t21 + t149 * t87 + t292 * t5 + t50 * t6) + t292 * t37 + ((t204 * t228 - t333 * t146 + t330 * t145 + (-Ifges(5,3) - t331) * t207) * qJD(1) + t97 + t58 + t57) * t265 / 0.2e1 + (t265 * t331 + t333 * t74) * t302 + (t265 * t333 + t334 * t74) * t310 + (-t146 * t334 - t207 * t333) * t320 + qJD(3) ^ 2 * (Ifges(4,5) * t207 - Ifges(4,6) * t204) / 0.2e1 + ((t133 * t264 + t38) * t207 + (-qJD(3) * t85 - t133 * t262 + t136 * t206) * t204) * mrSges(5,2) - t152 * t265 / 0.2e1 + t19 * (-mrSges(7,1) * t265 + mrSges(7,2) * t74) + t21 * (mrSges(6,1) * t265 - mrSges(6,3) * t74) + (t206 * t237 - t244 / 0.2e1) * t98 + (Ifges(7,5) * t74 + Ifges(7,6) * t265) * t313 + ((-t141 * t207 - t142 * t204) * qJD(3) + t219) * mrSges(4,3) - t76 * t275 / 0.2e1 + ((t133 * t266 - t39) * t207 + (qJD(3) * t84 + t133 * t260 + t136 * t203) * t204) * mrSges(5,1) + m(5) * (t39 * t113 + t38 * t114 + t85 * t72 + t84 * t73) + t149 * t17 + t72 * t130 + t73 * t131 + t118 * t68 + t113 * t107 + t114 * t108 + t11 * t92 + t9 * t93 + t7 * t90 + t10 * t91 + t80 * t16 + t18 * t67 + t47 * t34 + t48 * t36 + t50 * t35 + (mrSges(6,1) * t87 + mrSges(7,1) * t8 - mrSges(7,2) * t2 - mrSges(6,3) * t5 - Ifges(6,2) * t319 + t320 * t332 + t321 + t322 + t339) * t145 - t341 * t146 / 0.2e1 + (0.3e1 / 0.2e1 * t207 ^ 2 - 0.3e1 / 0.2e1 * t204 ^ 2) * Ifges(4,4) * t256 + t153 * t238 + (-Ifges(6,4) * t146 - Ifges(6,6) * t207) * t319 + (t146 * t8 - t2 * t207 + t20 * t265 - t31 * t74) * mrSges(7,3) + (-Ifges(7,5) * t146 - Ifges(7,6) * t207) * t318 + t6 * (-mrSges(6,1) * t207 + t146 * mrSges(6,3)) + t3 * (mrSges(7,1) * t207 - t146 * mrSges(7,2)) + (t103 * t74 - t146 * t87 + t207 * t5 - t22 * t265) * mrSges(6,2) + t74 * t342 + ((-Ifges(5,5) * t203 - Ifges(5,6) * t206) * t261 + (Ifges(5,3) * t204 + t207 * t228) * qJD(3)) * t301 + (-t231 * t261 + (Ifges(5,5) * t204 + t207 * t232) * qJD(3)) * t305 + (-Ifges(5,6) * t207 + t204 * t230) * t308 + (-Ifges(5,5) * t207 + t204 * t232) * t309 + t274 * t316; -t293 * t75 + t294 * t74 - t296 * t146 + t295 * t145 + (-t16 - t17 - t82) * t207 + (qJD(4) * t220 + t221) * t204 + ((t130 * t206 - t131 * t203 + t179 - t251) * t207 + (t67 + t68 - t252 - t269) * t204) * qJD(3) + m(4) * (t135 * t204 - t136 * t207 + (-t141 * t204 + t142 * t207) * qJD(3)) + m(7) * (t3 * t145 - t2 * t146 + t19 * t75 + t20 * t74 - t207 * t8 + t265 * t31) + m(6) * (t103 * t265 - t6 * t145 - t5 * t146 - t207 * t87 - t21 * t75 + t22 * t74) + m(5) * ((-qJD(3) * t223 - t136) * t207 + (qJD(3) * t133 - qJD(4) * t224 + t225) * t204); m(5) * (-pkin(3) * t136 + pkin(8) * t225) + (t215 * t297 + (t68 + t298) * t203 * pkin(4) + (-m(5) * t224 + t220) * pkin(8) + t324) * qJD(4) + t221 * pkin(8) + (-mrSges(5,1) * t206 + mrSges(5,2) * t203 - mrSges(4,1)) * t136 + (t59 - t56) * (t138 / 0.2e1 - t116 / 0.2e1) - t295 * t218 + (-t103 * t120 + t218 * t6 + t124 * t5 + t197 * t87 + (-t33 + t78) * t22 + (-t32 - t79) * t21) * m(6) + (-Ifges(6,4) * t115 - Ifges(7,5) * t139 - Ifges(6,2) * t116 + Ifges(7,3) * t138) * t314 + (-Ifges(6,4) * t139 - Ifges(7,5) * t115 - Ifges(6,2) * t138 + Ifges(7,3) * t116) * t313 + ((t141 * mrSges(4,3) - t153 / 0.2e1 - t178 * mrSges(4,2) - t198 / 0.2e1 + (Ifges(4,5) / 0.2e1 - t215 / 0.2e1) * qJD(3) - t324) * t207 + (-t58 / 0.2e1 - t57 / 0.2e1 + (t284 / 0.2e1 + Ifges(4,4) / 0.2e1) * t268 + t152 / 0.2e1 + t22 * mrSges(6,2) - t21 * mrSges(6,1) + t19 * mrSges(7,1) - t20 * mrSges(7,3) + t142 * mrSges(4,3) - t84 * mrSges(5,1) + t85 * mrSges(5,2) - t178 * mrSges(4,1) - t97 / 0.2e1 - t283 / 0.2e1 - t278 / 0.2e1 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t183 + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t210 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t110 - qJD(4) * t216 / 0.2e1 + (t216 / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t267 + (t266 / 0.2e1 + t306) * Ifges(5,5) + (t166 * t330 + t167 * t333) * t297) * t204) * qJD(1) + t225 * mrSges(5,3) - m(5) * (t100 * t84 + t101 * t85 + t133 * t142) - t293 * t79 + t294 * t78 + t296 * t124 + (mrSges(6,1) * t271 - mrSges(6,2) * t343) * t103 + (-t166 * t2 + t167 * t3 - t19 * t343 - t20 * t271) * mrSges(7,2) + (mrSges(7,1) * t271 + mrSges(7,3) * t343) * t31 + (-t166 * t5 - t167 * t6 + t21 * t343 - t22 * t271) * mrSges(6,3) + t269 * t142 + t197 * t17 - t141 * t179 + t8 * (mrSges(7,1) * t166 - mrSges(7,3) * t167) + t87 * (mrSges(6,1) * t166 + mrSges(6,2) * t167) - t101 * t130 - t100 * t131 - t135 * mrSges(4,2) - t120 * t68 + t106 * t16 - t32 * t92 - t30 * t93 - t29 * t90 - t33 * t91 - pkin(3) * t82 + (t106 * t8 - t218 * t3 + t124 * t2 + t340 * t31 + (-t29 + t78) * t20 + (-t30 + t79) * t19) * m(7) + t340 * t67 + t329 * (t139 / 0.2e1 - t115 / 0.2e1) + t341 * t167 / 0.2e1 + (-t115 * t333 + t116 * t330) * t302 + (t138 * t330 - t139 * t333) * t303 + (-t115 * t334 + t116 * t332) * t310 + (t138 * t332 - t139 * t334) * t311 + (t166 * t332 + t167 * t334) * t320 + t76 * t300 + t229 * t308 + t231 * t309 + t203 * t316 + (Ifges(7,5) * t167 + Ifges(7,3) * t166) * t318 + (Ifges(6,4) * t167 - Ifges(6,2) * t166) * t319 + t166 * t321 + t166 * t322; -t133 * (t165 * mrSges(5,1) - mrSges(5,2) * t214) - t191 * (-Ifges(5,5) * t214 - Ifges(5,6) * t165) / 0.2e1 + (t85 * t165 - t214 * t84) * mrSges(5,3) + (-Ifges(6,2) * t313 + Ifges(7,3) * t314 + t303 * t330 + t311 * t332 + t323) * t210 + t209 + (0.2e1 * t298 * t306 + m(7) * t19 * t259 + m(6) * (t202 * t5 + t205 * t6 - t21 * t259 + t22 * t258) - t165 * t68 + t202 * t37 + t205 * t35 + (-t202 * t293 + t205 * t91) * qJD(5)) * pkin(4) - m(6) * (-t21 * t23 + t22 * t24) + t293 * t23 - t294 * t24 + t250 + t194 * t34 + t196 * t36 + t187 * t90 - t84 * t130 + t85 * t131 - t53 * t67 + t39 * mrSges(5,1) - t38 * mrSges(5,2) + (-t19 * t23 + t194 * t2 + t196 * t3 - t31 * t53 + (t187 - t24) * t20) * m(7) + t98 * t305 + (-Ifges(5,1) * t214 - t289) * t306 + (-Ifges(5,2) * t165 - t161 + t99) * t214 / 0.2e1 + (t103 * mrSges(6,2) + t19 * mrSges(7,2) - t21 * mrSges(6,3) - t31 * mrSges(7,3) - Ifges(6,4) * t313 - Ifges(7,5) * t314 - t333 * t303 - t334 * t311 + t342) * t110; (t110 * t19 + t20 * t210) * mrSges(7,2) + t209 + (t290 + t293) * t22 + (-t291 - t294) * t21 - t103 * (mrSges(6,1) * t210 - mrSges(6,2) * t110) + t59 * t310 + (Ifges(7,3) * t210 - t285) * t314 - t31 * (mrSges(7,1) * t210 + mrSges(7,3) * t110) + qJD(6) * t90 - t66 * t67 + qJ(6) * t34 - pkin(5) * t36 + (-t110 * t333 + t210 * t330) * t303 + (-pkin(5) * t3 + qJ(6) * t2 - t19 * t22 + t20 * t328 - t31 * t66) * m(7) + (-Ifges(6,2) * t210 - t105 + t329) * t313 + (-t110 * t334 + t104 - t286 + t56) * t311; t210 * t67 - t183 * t90 + 0.2e1 * (t3 / 0.2e1 + t31 * t310 + t20 * t303) * m(7) + t36;];
tauc  = t1(:);
