% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:20
% EndTime: 2019-03-08 19:07:42
% DurationCPUTime: 9.96s
% Computational Cost: add. (10742->591), mult. (29788->888), div. (0->0), fcn. (26092->16), ass. (0->280)
t203 = cos(pkin(6));
t191 = qJD(1) * t203 + qJD(2);
t198 = sin(pkin(7));
t211 = cos(qJ(3));
t283 = t198 * t211;
t180 = t191 * t283;
t202 = cos(pkin(7));
t200 = cos(pkin(14));
t199 = sin(pkin(6));
t273 = qJD(1) * t199;
t258 = t200 * t273;
t196 = sin(pkin(14));
t207 = sin(qJ(3));
t289 = t196 * t207;
t122 = t202 * t211 * t258 - t273 * t289 + t180;
t284 = t198 * t207;
t282 = t200 * t202;
t341 = (t196 * t211 + t207 * t282) * t199;
t123 = qJD(1) * t341 + t191 * t284;
t197 = sin(pkin(8));
t206 = sin(qJ(4));
t287 = t197 * t206;
t192 = pkin(10) * t287;
t201 = cos(pkin(8));
t210 = cos(qJ(4));
t280 = t201 * t210;
t177 = pkin(3) * t280 - t192;
t281 = t201 * t206;
t342 = t177 * qJD(4) - t122 * t210 + t123 * t281;
t285 = t197 * t210;
t178 = pkin(3) * t281 + pkin(10) * t285;
t164 = pkin(11) * t201 + t178;
t245 = -pkin(4) * t210 - pkin(11) * t206;
t165 = (-pkin(3) + t245) * t197;
t244 = pkin(4) * t206 - pkin(11) * t210;
t224 = qJD(4) * t244;
t170 = t197 * t224;
t205 = sin(qJ(5));
t209 = cos(qJ(5));
t266 = qJD(5) * t209;
t267 = qJD(5) * t205;
t288 = t197 * t205;
t345 = -t123 * t288 - t164 * t267 + t165 * t266 + t205 * t170 + t342 * t209;
t343 = t178 * qJD(4) - t122 * t206 - t123 * t280;
t269 = qJD(4) * t197;
t255 = t206 * t269;
t354 = -pkin(12) * t255 - t345;
t175 = -t209 * t201 + t205 * t287;
t268 = qJD(4) * t210;
t253 = t209 * t268;
t138 = -qJD(5) * t175 + t197 * t253;
t286 = t197 * t209;
t176 = t201 * t205 + t206 * t286;
t254 = t205 * t268;
t139 = qJD(5) * t176 + t197 * t254;
t353 = pkin(5) * t139 - pkin(12) * t138 + t343;
t204 = sin(qJ(6));
t208 = cos(qJ(6));
t235 = Ifges(7,5) * t208 - Ifges(7,6) * t204;
t311 = Ifges(7,4) * t208;
t237 = -Ifges(7,2) * t204 + t311;
t312 = Ifges(7,4) * t204;
t239 = Ifges(7,1) * t208 - t312;
t270 = qJD(3) * t210;
t256 = t197 * t270;
t187 = qJD(5) - t256;
t190 = qJD(3) * t201 + qJD(4);
t118 = qJD(3) * pkin(3) + t122;
t153 = t191 * t202 - t198 * t258;
t232 = t118 * t201 + t153 * t197;
t272 = qJD(3) * t197;
t114 = pkin(10) * t272 + t123;
t292 = t114 * t210;
t57 = t206 * t232 + t292;
t50 = pkin(11) * t190 + t57;
t144 = t201 * t153;
t76 = t144 + (qJD(3) * t245 - t118) * t197;
t28 = -t205 * t50 + t209 * t76;
t24 = -pkin(5) * t187 - t28;
t240 = mrSges(7,1) * t204 + mrSges(7,2) * t208;
t29 = t205 * t76 + t209 * t50;
t25 = pkin(12) * t187 + t29;
t257 = t206 * t272;
t154 = t190 * t209 - t205 * t257;
t155 = t190 * t205 + t209 * t257;
t56 = -t206 * t114 + t232 * t210;
t49 = -pkin(4) * t190 - t56;
t39 = -pkin(5) * t154 - pkin(12) * t155 + t49;
t10 = t204 * t39 + t208 * t25;
t9 = -t204 * t25 + t208 * t39;
t242 = t10 * t204 + t208 * t9;
t321 = t208 / 0.2e1;
t322 = -t204 / 0.2e1;
t150 = qJD(6) - t154;
t325 = t150 / 0.2e1;
t121 = t155 * t208 + t187 * t204;
t329 = t121 / 0.2e1;
t120 = -t155 * t204 + t187 * t208;
t331 = t120 / 0.2e1;
t313 = Ifges(7,4) * t121;
t61 = Ifges(7,2) * t120 + Ifges(7,6) * t150 + t313;
t119 = Ifges(7,4) * t120;
t62 = Ifges(7,1) * t121 + Ifges(7,5) * t150 + t119;
t352 = -t242 * mrSges(7,3) + t235 * t325 + t237 * t331 + t239 * t329 + t24 * t240 + t321 * t62 + t322 * t61;
t340 = t209 * t164 + t205 * t165;
t344 = -qJD(5) * t340 - t123 * t286 + t170 * t209 - t205 * t342;
t310 = Ifges(6,5) * t187;
t149 = Ifges(6,4) * t154;
t316 = Ifges(6,1) * t155;
t94 = t149 + t310 + t316;
t219 = t28 * mrSges(6,3) - t94 / 0.2e1 - t310 / 0.2e1 - t49 * mrSges(6,2);
t351 = t219 - t352;
t135 = t203 * t284 + t341;
t350 = qJD(3) * t123;
t349 = -t190 * Ifges(5,6) / 0.2e1;
t163 = t192 + (-pkin(3) * t210 - pkin(4)) * t201;
t100 = pkin(5) * t175 - pkin(12) * t176 + t163;
t102 = -pkin(12) * t285 + t340;
t59 = t100 * t204 + t102 * t208;
t348 = -qJD(6) * t59 + t204 * t354 + t353 * t208;
t58 = t100 * t208 - t102 * t204;
t347 = qJD(6) * t58 + t353 * t204 - t208 * t354;
t346 = -pkin(5) * t255 - t344;
t291 = t350 * t197;
t274 = mrSges(5,1) * t190 + mrSges(6,1) * t154 - mrSges(6,2) * t155 - mrSges(5,3) * t257;
t275 = t210 * t211;
t279 = t206 * t207;
t339 = t201 * t275 - t279;
t128 = t190 * t266 + (-t206 * t267 + t253) * t272;
t129 = t190 * t267 + (t206 * t266 + t254) * t272;
t222 = (t211 * t282 - t289) * t199;
t115 = (qJD(1) * t222 + t180) * qJD(3);
t35 = qJD(4) * t57 + t115 * t206 + t280 * t350;
t17 = pkin(5) * t129 - pkin(12) * t128 + t35;
t252 = qJD(3) * t269;
t249 = t206 * t252;
t34 = qJD(4) * t56 + t115 * t210 - t281 * t350;
t88 = (qJD(3) * t224 + t350) * t197;
t7 = t205 * t88 + t209 * t34 + t76 * t266 - t267 * t50;
t5 = pkin(12) * t249 + t7;
t1 = qJD(6) * t9 + t17 * t204 + t208 * t5;
t2 = -qJD(6) * t10 + t17 * t208 - t204 * t5;
t338 = t1 * t208 - t2 * t204;
t8 = -qJD(5) * t29 - t205 * t34 + t209 * t88;
t337 = -t8 * mrSges(6,1) + t7 * mrSges(6,2) - Ifges(6,5) * t128 + Ifges(6,6) * t129;
t72 = qJD(6) * t120 + t128 * t208 + t204 * t249;
t73 = -qJD(6) * t121 - t128 * t204 + t208 * t249;
t33 = Ifges(7,1) * t72 + Ifges(7,4) * t73 + Ifges(7,5) * t129;
t336 = t33 / 0.2e1;
t335 = -t61 / 0.2e1;
t334 = t72 / 0.2e1;
t333 = t73 / 0.2e1;
t332 = -t120 / 0.2e1;
t330 = -t121 / 0.2e1;
t328 = t129 / 0.2e1;
t327 = -t149 / 0.2e1;
t326 = -t150 / 0.2e1;
t324 = -t175 / 0.2e1;
t323 = t176 / 0.2e1;
t68 = Ifges(7,5) * t72;
t67 = Ifges(7,6) * t73;
t320 = pkin(11) * t209;
t134 = t203 * t283 + t222;
t173 = -t198 * t199 * t200 + t202 * t203;
t290 = t135 * t206;
t65 = -t134 * t280 - t173 * t285 + t290;
t317 = t35 * t65;
t315 = Ifges(5,4) * t206;
t314 = Ifges(6,4) * t155;
t309 = Ifges(7,5) * t121;
t308 = Ifges(6,2) * t154;
t307 = Ifges(6,6) * t187;
t306 = Ifges(7,6) * t120;
t305 = Ifges(7,3) * t150;
t304 = t128 * Ifges(6,1);
t303 = t128 * Ifges(6,4);
t302 = t129 * Ifges(6,4);
t261 = t202 * t285;
t136 = -t198 * t339 - t261;
t301 = t136 * t35;
t104 = mrSges(6,1) * t249 - mrSges(6,3) * t128;
t36 = -mrSges(7,1) * t73 + mrSges(7,2) * t72;
t294 = t36 - t104;
t169 = t244 * t272;
t45 = t205 * t169 + t209 * t56;
t133 = mrSges(6,1) * t187 - mrSges(6,3) * t155;
t81 = -mrSges(7,1) * t120 + mrSges(7,2) * t121;
t293 = t81 - t133;
t278 = t206 * t211;
t277 = t207 * t210;
t276 = t209 * t210;
t271 = qJD(3) * t207;
t186 = -pkin(5) * t209 - pkin(12) * t205 - pkin(4);
t265 = qJD(6) * t186;
t264 = qJD(6) * t209;
t31 = Ifges(7,3) * t129 + t67 + t68;
t251 = t272 / 0.2e1;
t250 = t197 * t198 * t271;
t247 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t246 = -t35 * mrSges(5,1) - t34 * mrSges(5,2);
t243 = pkin(5) * t205 - pkin(12) * t209;
t241 = mrSges(7,1) * t208 - mrSges(7,2) * t204;
t238 = Ifges(7,1) * t204 + t311;
t236 = Ifges(7,2) * t208 + t312;
t234 = Ifges(7,5) * t204 + Ifges(7,6) * t208;
t231 = t134 * t201 + t173 * t197;
t66 = t135 * t210 + t206 * t231;
t97 = -t134 * t197 + t173 * t201;
t41 = t205 * t97 + t209 * t66;
t23 = t204 * t65 + t208 * t41;
t22 = -t204 * t41 + t208 * t65;
t40 = t205 * t66 - t97 * t209;
t225 = t201 * t278 + t277;
t137 = t198 * t225 + t202 * t287;
t174 = -t197 * t283 + t201 * t202;
t99 = t137 * t209 + t174 * t205;
t70 = t136 * t208 - t204 * t99;
t71 = t136 * t204 + t208 * t99;
t44 = t169 * t209 - t205 * t56;
t98 = t137 * t205 - t209 * t174;
t111 = -t164 * t205 + t165 * t209;
t140 = -t176 * t204 - t208 * t285;
t227 = -t176 * t208 + t204 * t285;
t189 = Ifges(5,4) * t256;
t86 = -t118 * t197 + t144;
t218 = -t56 * mrSges(5,3) + t86 * mrSges(5,2) + t190 * Ifges(5,5) + Ifges(5,1) * t257 / 0.2e1 + t189 / 0.2e1;
t217 = t28 * mrSges(6,1) + t86 * mrSges(5,1) + t349 - (Ifges(5,2) * t210 + t315) * t272 / 0.2e1 + t187 * Ifges(6,3) + t155 * Ifges(6,5) + t154 * Ifges(6,6) - t29 * mrSges(6,2) - t57 * mrSges(5,3);
t60 = t305 + t306 + t309;
t93 = t307 + t308 + t314;
t216 = -t305 / 0.2e1 - t306 / 0.2e1 - t309 / 0.2e1 + t307 / 0.2e1 + t93 / 0.2e1 - t60 / 0.2e1 - t49 * mrSges(6,1) - t9 * mrSges(7,1) + t10 * mrSges(7,2) + t29 * mrSges(6,3) + t314 / 0.2e1;
t215 = t308 / 0.2e1 + t216;
t185 = Ifges(5,5) * t210 * t252;
t184 = Ifges(6,3) * t249;
t182 = t243 * qJD(5);
t168 = (-mrSges(5,1) * t210 + mrSges(5,2) * t206) * t272;
t167 = -mrSges(5,2) * t190 + mrSges(5,3) * t256;
t161 = t186 * t204 + t208 * t320;
t160 = t186 * t208 - t204 * t320;
t159 = (mrSges(5,1) * t206 + mrSges(5,2) * t210) * t252;
t147 = (t204 * t206 + t208 * t276) * t272;
t146 = (-t204 * t276 + t206 * t208) * t272;
t132 = -mrSges(6,2) * t187 + mrSges(6,3) * t154;
t131 = t135 * qJD(3);
t130 = t134 * qJD(3);
t110 = pkin(5) * t155 - pkin(12) * t154;
t109 = -t204 * t265 + t182 * t208 + (t204 * t267 - t208 * t264) * pkin(11);
t108 = t208 * t265 + t182 * t204 + (-t204 * t264 - t208 * t267) * pkin(11);
t105 = -mrSges(6,2) * t249 - mrSges(6,3) * t129;
t101 = pkin(5) * t285 - t111;
t96 = t202 * t255 + (t225 * qJD(4) + (t201 * t277 + t278) * qJD(3)) * t198;
t95 = qJD(4) * t261 + (t339 * qJD(4) + (-t201 * t279 + t275) * qJD(3)) * t198;
t90 = mrSges(7,1) * t150 - mrSges(7,3) * t121;
t89 = -mrSges(7,2) * t150 + mrSges(7,3) * t120;
t85 = qJD(6) * t227 - t138 * t204 + t208 * t255;
t84 = qJD(6) * t140 + t138 * t208 + t204 * t255;
t83 = mrSges(6,1) * t129 + mrSges(6,2) * t128;
t80 = Ifges(6,5) * t249 - t302 + t304;
t79 = -t129 * Ifges(6,2) + Ifges(6,6) * t249 + t303;
t54 = -mrSges(7,2) * t129 + mrSges(7,3) * t73;
t53 = mrSges(7,1) * t129 - mrSges(7,3) * t72;
t52 = qJD(5) * t99 + t205 * t95 - t209 * t250;
t51 = -qJD(5) * t98 + t205 * t250 + t209 * t95;
t46 = t118 * t281 + t292 + (t153 * t206 + t243 * t270) * t197;
t43 = pkin(12) * t257 + t45;
t42 = -pkin(5) * t257 - t44;
t38 = -t131 * t281 + t130 * t210 + (t210 * t231 - t290) * qJD(4);
t37 = qJD(4) * t66 + t130 * t206 + t131 * t280;
t32 = Ifges(7,4) * t72 + Ifges(7,2) * t73 + Ifges(7,6) * t129;
t21 = t110 * t204 + t208 * t28;
t20 = t110 * t208 - t204 * t28;
t19 = -qJD(6) * t71 - t204 * t51 + t208 * t96;
t18 = qJD(6) * t70 + t204 * t96 + t208 * t51;
t14 = t204 * t46 + t208 * t43;
t13 = -t204 * t43 + t208 * t46;
t12 = -qJD(5) * t40 + t131 * t288 + t209 * t38;
t11 = qJD(5) * t41 - t131 * t286 + t205 * t38;
t6 = -pkin(5) * t249 - t8;
t4 = qJD(6) * t22 + t12 * t208 + t204 * t37;
t3 = -qJD(6) * t23 - t12 * t204 + t208 * t37;
t15 = [t131 * t197 * t168 + t41 * t105 + t12 * t132 + t97 * t159 + t38 * t167 + t22 * t53 + t23 * t54 + t3 * t90 + t4 * t89 + t65 * t83 + t294 * t40 - t274 * t37 + t293 * t11 + (-t131 * mrSges(4,1) - t130 * mrSges(4,2) + (-t206 * t66 + t210 * t65) * mrSges(5,3) * t269) * qJD(3) + m(6) * (-t11 * t28 + t12 * t29 + t37 * t49 - t40 * t8 + t41 * t7 + t317) + m(7) * (t1 * t23 + t10 * t4 + t11 * t24 + t2 * t22 + t3 * t9 + t40 * t6) + m(4) * (t115 * t135 - t122 * t131 + t123 * t130 - t134 * t350) + m(5) * (t34 * t66 + t317 - t37 * t56 + t38 * t57 + (t131 * t86 + t350 * t97) * t197); t99 * t105 + t51 * t132 + t136 * t83 + t174 * t159 + t95 * t167 + t18 * t89 + t19 * t90 + t70 * t53 + t71 * t54 + t294 * t98 - t274 * t96 + t293 * t52 + (t168 * t284 + (t136 * t210 - t137 * t206) * qJD(4) * mrSges(5,3)) * t272 + m(6) * (-t28 * t52 + t29 * t51 + t49 * t96 + t7 * t99 - t8 * t98 + t301) + m(7) * (t1 * t71 + t10 * t18 + t19 * t9 + t2 * t70 + t24 * t52 + t6 * t98) + m(5) * (t137 * t34 + t174 * t291 - t56 * t96 + t57 * t95 + t301) + ((-mrSges(4,1) * t207 - mrSges(4,2) * t211) * qJD(3) ^ 2 + 0.2e1 * m(5) * t207 * t86 * t251 + m(4) * (t115 * t207 - t122 * t271)) * t198; (t123 * mrSges(4,1) + t122 * mrSges(4,2) + ((Ifges(5,5) * t201 / 0.2e1 - t177 * mrSges(5,3) + 0.3e1 / 0.2e1 * Ifges(5,4) * t285) * t210 + (-Ifges(5,6) * t201 + Ifges(6,5) * t323 + Ifges(6,6) * t324 - t178 * mrSges(5,3) - 0.3e1 / 0.2e1 * Ifges(5,4) * t287 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - Ifges(6,3) / 0.2e1) * t285) * t206) * t269) * qJD(3) + t342 * t167 - t274 * t343 + (-t177 * t35 + t178 * t34 + t342 * t57 - t343 * t56 + (-pkin(3) * t291 - t123 * t86) * t197) * m(5) + (-pkin(3) * t159 - t123 * t168 + (mrSges(5,2) * t291 + mrSges(5,3) * t35) * t206 + (t34 * mrSges(5,3) - t184 / 0.2e1 - mrSges(5,1) * t291 + t337) * t210 + (t218 * t210 + (t349 + t217) * t206) * qJD(4)) * t197 + t340 * t105 + (t111 * t8 + t163 * t35 + t28 * t344 + t29 * t345 + t340 * t7 + t343 * t49) * m(6) + t6 * (-mrSges(7,1) * t140 - mrSges(7,2) * t227) + (-Ifges(7,5) * t227 + Ifges(7,6) * t140 + Ifges(7,3) * t175) * t328 + (-Ifges(7,4) * t227 + Ifges(7,2) * t140 + Ifges(7,6) * t175) * t333 + (-Ifges(7,1) * t227 + Ifges(7,4) * t140 + Ifges(7,5) * t175) * t334 + t2 * (mrSges(7,1) * t175 + mrSges(7,3) * t227) - t227 * t336 + t187 * (Ifges(6,5) * t138 - Ifges(6,6) * t139) / 0.2e1 - t350 * mrSges(4,1) + (t185 / 0.2e1 + t246) * t201 + t128 * (Ifges(6,1) * t176 - Ifges(6,4) * t175) / 0.2e1 + t155 * (Ifges(6,1) * t138 - Ifges(6,4) * t139) / 0.2e1 + (-t138 * t28 - t139 * t29 - t175 * t7 - t176 * t8) * mrSges(6,3) + t35 * (mrSges(6,1) * t175 + mrSges(6,2) * t176) + t175 * t31 / 0.2e1 + t1 * (-mrSges(7,2) * t175 + mrSges(7,3) * t140) + t163 * t83 + t10 * (-mrSges(7,2) * t139 + mrSges(7,3) * t85) + t139 * t60 / 0.2e1 + t49 * (mrSges(6,1) * t139 + mrSges(6,2) * t138) - t139 * t93 / 0.2e1 + t140 * t32 / 0.2e1 + t138 * t94 / 0.2e1 + t9 * (mrSges(7,1) * t139 - mrSges(7,3) * t84) - t115 * mrSges(4,2) + t111 * t104 + t101 * t36 + t24 * (-mrSges(7,1) * t85 + mrSges(7,2) * t84) + t85 * t61 / 0.2e1 + t84 * t62 / 0.2e1 + t58 * t53 + t59 * t54 + t344 * t133 + t345 * t132 + t346 * t81 + t347 * t89 + t348 * t90 + (t1 * t59 + t10 * t347 + t101 * t6 + t2 * t58 + t24 * t346 + t348 * t9) * m(7) - t129 * (Ifges(6,4) * t176 - Ifges(6,2) * t175) / 0.2e1 + t154 * (Ifges(6,4) * t138 - Ifges(6,2) * t139) / 0.2e1 + t80 * t323 + t79 * t324 + (Ifges(7,5) * t84 + Ifges(7,6) * t85 + Ifges(7,3) * t139) * t325 + (Ifges(7,1) * t84 + Ifges(7,4) * t85 + Ifges(7,5) * t139) * t329 + (Ifges(7,4) * t84 + Ifges(7,2) * t85 + Ifges(7,6) * t139) * t331; t246 + t274 * t57 + (t303 / 0.2e1 - t35 * mrSges(6,1) + t79 / 0.2e1 - t31 / 0.2e1 - t68 / 0.2e1 - t67 / 0.2e1 + t7 * mrSges(6,3) + (-Ifges(6,2) / 0.2e1 - Ifges(7,3) / 0.2e1) * t129 + (m(6) * t7 + t105) * pkin(11) + t247) * t209 + ((qJD(4) * (Ifges(6,5) * t205 + Ifges(6,6) * t209) / 0.2e1 + t251 * t315 + (t190 / 0.2e1 - qJD(4)) * Ifges(5,6) - t217) * t206 + (-t189 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t257 + (-t316 / 0.2e1 + t327 + t219) * t209 + t215 * t205 - t218) * t210) * t272 - m(7) * (t10 * t14 + t13 * t9 + t24 * t42) - m(6) * (t28 * t44 + t29 * t45 + t49 * t57) + (t108 - t14) * t89 + ((-t215 + (-m(6) * t29 - t132) * pkin(11)) * t205 + (t316 / 0.2e1 + t149 / 0.2e1 + (-m(6) * t28 + m(7) * t24 + t293) * pkin(11) - t351) * t209) * qJD(5) + (-t10 * t146 + t147 * t9) * mrSges(7,3) + (-t302 / 0.2e1 - t8 * mrSges(6,3) + t304 / 0.2e1 + t35 * mrSges(6,2) + t32 * t322 + t33 * t321 + t6 * t240 + t235 * t328 + t239 * t334 + t237 * t333 + t80 / 0.2e1 + (-t1 * t204 - t2 * t208) * mrSges(7,3) + (-m(6) * t8 + m(7) * t6 + t294) * pkin(11) + (t24 * t241 + t234 * t326 + t236 * t332 + t238 * t330 + t208 * t335 + t62 * t322 + (-t10 * t208 + t204 * t9) * mrSges(7,3)) * qJD(6)) * t205 + (t109 - t13) * t90 + m(7) * (t1 * t161 + t10 * t108 + t109 * t9 + t160 * t2) + t146 * t335 - t56 * t167 + t160 * t53 + t161 * t54 - t24 * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) - t147 * t62 / 0.2e1 - t45 * t132 - t44 * t133 - t42 * t81 + (-m(6) * t35 - t83) * pkin(4) + t185 + (Ifges(7,5) * t147 + Ifges(7,6) * t146) * t326 + (Ifges(7,1) * t147 + Ifges(7,4) * t146) * t330 + (Ifges(7,4) * t147 + Ifges(7,2) * t146) * t332; -t293 * t29 - t337 + t352 * qJD(6) + t216 * t155 + (t327 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t155 + t351) * t154 + t338 * mrSges(7,3) + (-t204 * t53 + t208 * t54 + m(7) * t338 + (-m(7) * t242 - t204 * t89 - t208 * t90) * qJD(6)) * pkin(12) - t6 * t241 + (-pkin(5) * t6 - t10 * t21 - t20 * t9 - t24 * t29) * m(7) + t204 * t336 - t28 * t132 - t21 * t89 - t20 * t90 - pkin(5) * t36 + t184 + t32 * t321 + t234 * t328 + t236 * t333 + t238 * t334; (Ifges(7,1) * t120 - t313) * t330 + t61 * t329 + (Ifges(7,5) * t120 - Ifges(7,6) * t121) * t326 - t9 * t89 - t24 * (mrSges(7,1) * t121 + mrSges(7,2) * t120) + t10 * t90 + (t10 * t121 + t120 * t9) * mrSges(7,3) - t247 + t31 + (-Ifges(7,2) * t121 + t119 + t62) * t332;];
tauc  = t15(:);
