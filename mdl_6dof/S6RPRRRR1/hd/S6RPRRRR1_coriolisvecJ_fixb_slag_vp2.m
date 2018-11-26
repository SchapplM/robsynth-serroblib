% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:32:24
% EndTime: 2018-11-23 16:32:29
% DurationCPUTime: 5.28s
% Computational Cost: add. (13451->510), mult. (32101->710), div. (0->0), fcn. (22815->10), ass. (0->237)
t213 = qJD(3) + qJD(4);
t211 = qJD(5) + t213;
t216 = sin(qJ(6));
t220 = cos(qJ(6));
t218 = sin(qJ(4));
t219 = sin(qJ(3));
t222 = cos(qJ(4));
t223 = cos(qJ(3));
t194 = -t218 * t219 + t222 * t223;
t189 = t194 * qJD(1);
t195 = t218 * t223 + t222 * t219;
t190 = t195 * qJD(1);
t217 = sin(qJ(5));
t221 = cos(qJ(5));
t242 = t189 * t217 + t221 * t190;
t117 = t211 * t220 - t216 * t242;
t329 = -t117 / 0.2e1;
t118 = t211 * t216 + t220 * t242;
t328 = -t118 / 0.2e1;
t258 = t221 * t189 - t190 * t217;
t128 = qJD(6) - t258;
t326 = -t128 / 0.2e1;
t153 = t213 * t194;
t141 = t153 * qJD(1);
t154 = t213 * t195;
t142 = t154 * qJD(1);
t71 = qJD(5) * t258 + t141 * t221 - t142 * t217;
t42 = qJD(6) * t117 + t220 * t71;
t43 = -qJD(6) * t118 - t216 * t71;
t14 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t316 = t141 * pkin(9);
t204 = sin(pkin(11)) * pkin(1) + pkin(7);
t197 = t204 * qJD(1);
t259 = pkin(8) * qJD(1) + t197;
t267 = t219 * qJD(2);
t161 = t223 * t259 + t267;
t157 = t222 * t161;
t212 = t223 * qJD(2);
t160 = -t259 * t219 + t212;
t296 = qJD(3) * pkin(3);
t158 = t160 + t296;
t107 = t158 * t218 + t157;
t317 = pkin(9) * t189;
t102 = t107 + t317;
t276 = t221 * t102;
t155 = t218 * t161;
t106 = t222 * t158 - t155;
t183 = t190 * pkin(9);
t101 = t106 - t183;
t94 = pkin(4) * t213 + t101;
t49 = t217 * t94 + t276;
t111 = t222 * t160 - t155;
t272 = qJD(4) * t222;
t273 = qJD(4) * t218;
t64 = qJD(3) * t111 + t158 * t272 - t161 * t273;
t50 = -pkin(9) * t142 + t64;
t174 = t197 * t223 + t267;
t274 = qJD(1) * t223;
t275 = qJD(1) * t219;
t341 = -t197 * t219 + t212;
t228 = (-t218 * (-pkin(8) * t275 + t341) + t222 * (-pkin(8) * t274 - t174)) * qJD(3);
t65 = -qJD(4) * t107 + t228;
t9 = t217 * t50 - t221 * (t65 - t316) + t49 * qJD(5);
t347 = m(7) * t9 + t14;
t254 = mrSges(7,1) * t216 + mrSges(7,2) * t220;
t285 = t102 * t217;
t48 = t221 * t94 - t285;
t44 = -pkin(5) * t211 - t48;
t237 = t44 * t254;
t249 = Ifges(7,5) * t220 - Ifges(7,6) * t216;
t303 = Ifges(7,4) * t220;
t251 = -Ifges(7,2) * t216 + t303;
t304 = Ifges(7,4) * t216;
t253 = Ifges(7,1) * t220 - t304;
t320 = t220 / 0.2e1;
t321 = -t216 / 0.2e1;
t327 = t118 / 0.2e1;
t305 = Ifges(7,4) * t118;
t58 = Ifges(7,2) * t117 + Ifges(7,6) * t128 + t305;
t116 = Ifges(7,4) * t117;
t59 = Ifges(7,1) * t118 + Ifges(7,5) * t128 + t116;
t346 = t117 * t251 / 0.2e1 + t253 * t327 + t128 * t249 / 0.2e1 + t237 + t58 * t321 + t59 * t320;
t306 = Ifges(6,4) * t242;
t127 = Ifges(6,4) * t258;
t290 = mrSges(6,1) * t211 + mrSges(7,1) * t117 - mrSges(7,2) * t118 - mrSges(6,3) * t242;
t103 = -t183 + t111;
t207 = pkin(3) * t222 + pkin(4);
t110 = -t160 * t218 - t157;
t239 = t110 - t317;
t270 = qJD(5) * t221;
t271 = qJD(5) * t217;
t277 = t218 * t221;
t344 = t103 * t217 - t221 * t239 - t207 * t271 - (t218 * t270 + (t217 * t222 + t277) * qJD(4)) * pkin(3);
t278 = t217 * t218;
t139 = t207 * t270 + (-t218 * t271 + (t221 * t222 - t278) * qJD(4)) * pkin(3);
t56 = t221 * t103 + t217 * t239;
t343 = -t56 + t139;
t311 = pkin(8) + t204;
t192 = t311 * t219;
t193 = t311 * t223;
t138 = -t218 * t192 + t222 * t193;
t45 = pkin(10) * t211 + t49;
t264 = -cos(pkin(11)) * pkin(1) - pkin(2);
t196 = -pkin(3) * t223 + t264;
t191 = qJD(1) * t196;
t148 = -t189 * pkin(4) + t191;
t74 = -pkin(5) * t258 - pkin(10) * t242 + t148;
t19 = -t216 * t45 + t220 * t74;
t20 = t216 * t74 + t220 * t45;
t340 = -t19 * t216 + t20 * t220;
t209 = pkin(3) * t275;
t123 = pkin(4) * t142 + qJD(3) * t209;
t72 = qJD(5) * t242 + t141 * t217 + t221 * t142;
t21 = pkin(5) * t72 - pkin(10) * t71 + t123;
t8 = t48 * qJD(5) + t221 * t50 + (-t158 * t273 - t161 * t272 + t228 - t316) * t217;
t2 = qJD(6) * t19 + t21 * t216 + t220 * t8;
t286 = qJD(6) * t20;
t3 = t21 * t220 - t216 * t8 - t286;
t339 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t42 + Ifges(7,6) * t43;
t338 = Ifges(6,6) * t211 / 0.2e1 + t306 / 0.2e1;
t95 = pkin(5) * t242 - pkin(10) * t258;
t333 = Ifges(7,5) * t328 + Ifges(7,6) * t329 + Ifges(7,3) * t326;
t246 = t19 * t220 + t20 * t216;
t302 = Ifges(6,5) * t211;
t309 = Ifges(6,1) * t242;
t88 = t127 + t302 + t309;
t337 = t148 * mrSges(6,2) + t88 / 0.2e1 + t302 / 0.2e1 + t127 / 0.2e1 - t246 * mrSges(7,3) + t346;
t335 = t42 / 0.2e1;
t334 = t43 / 0.2e1;
t332 = t72 / 0.2e1;
t300 = Ifges(6,2) * t258;
t331 = t300 / 0.2e1 + t338;
t137 = -t222 * t192 - t193 * t218;
t119 = -pkin(9) * t195 + t137;
t120 = pkin(9) * t194 + t138;
t243 = t221 * t119 - t120 * t217;
t330 = t243 * t9;
t324 = t189 / 0.2e1;
t323 = -t190 / 0.2e1;
t322 = t190 / 0.2e1;
t319 = m(5) * t191;
t318 = pkin(4) * t190;
t241 = t221 * t194 - t195 * t217;
t315 = t241 * t9;
t314 = t2 * t220;
t313 = t216 * t3;
t312 = t48 * mrSges(6,3);
t310 = mrSges(5,3) * t189;
t308 = Ifges(4,4) * t219;
t307 = Ifges(5,4) * t190;
t295 = t242 * t49;
t293 = t190 * mrSges(5,3);
t199 = t264 * qJD(1);
t292 = t199 * mrSges(4,2);
t288 = Ifges(4,5) * qJD(3);
t287 = Ifges(4,6) * qJD(3);
t284 = t258 * t216;
t283 = t258 * t220;
t281 = t194 * t141;
t280 = t195 * t142;
t185 = pkin(3) * t277 + t217 * t207;
t269 = qJD(6) * t216;
t268 = qJD(6) * t220;
t208 = Ifges(4,4) * t274;
t263 = t288 / 0.2e1;
t262 = -t287 / 0.2e1;
t261 = m(4) * t204 + mrSges(4,3);
t134 = pkin(4) * t154 + t219 * t296;
t260 = qJD(3) * t311;
t257 = t287 / 0.2e1 + (t223 * Ifges(4,2) + t308) * qJD(1) / 0.2e1 - t199 * mrSges(4,1);
t255 = mrSges(7,1) * t220 - mrSges(7,2) * t216;
t252 = Ifges(7,1) * t216 + t303;
t250 = Ifges(7,2) * t220 + t304;
t248 = Ifges(7,5) * t216 + Ifges(7,6) * t220;
t17 = mrSges(7,1) * t72 - mrSges(7,3) * t42;
t18 = -mrSges(7,2) * t72 + mrSges(7,3) * t43;
t247 = -t216 * t17 + t220 * t18;
t77 = t119 * t217 + t120 * t221;
t147 = t194 * t217 + t195 * t221;
t159 = -t194 * pkin(4) + t196;
t85 = -pkin(5) * t241 - t147 * pkin(10) + t159;
t30 = t216 * t85 + t220 * t77;
t29 = -t216 * t77 + t220 * t85;
t184 = -pkin(3) * t278 + t207 * t221;
t121 = -mrSges(6,2) * t211 + mrSges(6,3) * t258;
t83 = -mrSges(7,2) * t128 + mrSges(7,3) * t117;
t84 = mrSges(7,1) * t128 - mrSges(7,3) * t118;
t238 = -t216 * t84 + t220 * t83 + t121;
t180 = t219 * t260;
t181 = t223 * t260;
t96 = -t222 * t180 - t218 * t181 - t192 * t272 - t193 * t273;
t82 = t318 + t95;
t233 = -qJD(6) * t246 - t313;
t232 = t233 * mrSges(7,3);
t97 = -qJD(4) * t138 + t180 * t218 - t222 * t181;
t231 = t233 + t314;
t230 = -pkin(9) * t153 + t97;
t12 = t42 * Ifges(7,4) + t43 * Ifges(7,2) + t72 * Ifges(7,6);
t13 = t42 * Ifges(7,1) + t43 * Ifges(7,4) + t72 * Ifges(7,5);
t229 = -t8 * mrSges(6,2) + mrSges(7,3) * t314 + t216 * t13 / 0.2e1 + t12 * t320 + t252 * t335 + t250 * t334 + t248 * t332 - Ifges(6,6) * t72 + Ifges(6,5) * t71 + (-mrSges(6,1) - t255) * t9 + t346 * qJD(6);
t227 = -t148 * mrSges(6,1) - t19 * mrSges(7,1) + t20 * mrSges(7,2) + t331 + 0.2e1 * t333 + t338;
t226 = m(7) * (-t19 * t268 - t20 * t269 - t313 + t314) - t83 * t269 - t84 * t268 + t247;
t125 = t189 * Ifges(5,2) + t213 * Ifges(5,6) + t307;
t182 = Ifges(5,4) * t189;
t126 = t190 * Ifges(5,1) + t213 * Ifges(5,5) + t182;
t224 = t125 * t322 + (Ifges(5,1) * t189 - t307) * t323 + t106 * t310 + t229 - (-Ifges(5,2) * t190 + t126 + t182) * t189 / 0.2e1 - (-Ifges(6,2) * t242 + t127 + t88) * t258 / 0.2e1 - t64 * mrSges(5,2) + t65 * mrSges(5,1) + t258 * t312 + (Ifges(7,3) * t242 + t249 * t258) * t326 + (Ifges(7,5) * t242 + t253 * t258) * t328 + (Ifges(7,6) * t242 + t251 * t258) * t329 - t148 * (mrSges(6,1) * t242 + mrSges(6,2) * t258) - t211 * (Ifges(6,5) * t258 - Ifges(6,6) * t242) / 0.2e1 - t191 * (mrSges(5,1) * t190 + mrSges(5,2) * t189) + t242 * t331 + t242 * t333 - t20 * (-mrSges(7,2) * t242 - mrSges(7,3) * t284) - t19 * (mrSges(7,1) * t242 - mrSges(7,3) * t283) - t242 * (Ifges(6,1) * t258 - t306) / 0.2e1 - t258 * t237 + Ifges(5,5) * t141 - Ifges(5,6) * t142 - t59 * t283 / 0.2e1 + t58 * t284 / 0.2e1 - t213 * (Ifges(5,5) * t189 - Ifges(5,6) * t190) / 0.2e1;
t200 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t274;
t198 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t275;
t188 = Ifges(4,1) * t275 + t208 + t288;
t179 = pkin(10) + t185;
t178 = -pkin(5) - t184;
t166 = t174 * qJD(3);
t165 = t341 * qJD(3);
t164 = mrSges(5,1) * t213 - t293;
t163 = -mrSges(5,2) * t213 + t310;
t162 = t209 + t318;
t144 = -mrSges(5,1) * t189 + mrSges(5,2) * t190;
t93 = -mrSges(6,1) * t258 + mrSges(6,2) * t242;
t81 = qJD(5) * t147 + t153 * t217 + t221 * t154;
t80 = qJD(5) * t241 + t153 * t221 - t154 * t217;
t79 = t209 + t82;
t75 = -pkin(9) * t154 + t96;
t68 = Ifges(7,3) * t72;
t52 = t101 * t221 - t285;
t51 = t101 * t217 + t276;
t28 = pkin(5) * t81 - pkin(10) * t80 + t134;
t27 = t216 * t95 + t220 * t48;
t26 = -t216 * t48 + t220 * t95;
t25 = t216 * t79 + t220 * t56;
t24 = -t216 * t56 + t220 * t79;
t23 = t216 * t82 + t220 * t52;
t22 = -t216 * t52 + t220 * t82;
t16 = qJD(5) * t77 + t217 * t75 - t221 * t230;
t15 = qJD(5) * t243 + t217 * t230 + t221 * t75;
t5 = -qJD(6) * t30 - t15 * t216 + t220 * t28;
t4 = qJD(6) * t29 + t15 * t220 + t216 * t28;
t1 = [t15 * t121 + t191 * (mrSges(5,1) * t154 + mrSges(5,2) * t153) + t213 * (Ifges(5,5) * t153 - Ifges(5,6) * t154) / 0.2e1 + m(5) * (t106 * t97 + t107 * t96 + t137 * t65 + t138 * t64) + t4 * t83 + t5 * t84 - t243 * t14 + (-t243 * t71 - t48 * t80 - t49 * t81 - t72 * t77) * mrSges(6,3) - (-t8 * mrSges(6,3) + t68 / 0.2e1 - Ifges(6,4) * t71 + t123 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t72 + t339) * t241 + (t261 * t165 + (t188 / 0.2e1 - t204 * t198 + 0.3e1 / 0.2e1 * t208 + t263 - t261 * t341 + 0.2e1 * t292) * qJD(3)) * t223 + (Ifges(6,1) * t71 - Ifges(6,4) * t72 + t123 * mrSges(6,2) + t12 * t321 + t249 * t332 + t253 * t335 + t251 * t334 + t13 * t320 + (mrSges(6,3) + t254) * t9 + (-t2 * t216 - t220 * t3) * mrSges(7,3) + (t248 * t326 + t252 * t328 + t250 * t329 + t44 * t255 + t59 * t321 - t220 * t58 / 0.2e1 - t340 * mrSges(7,3)) * qJD(6)) * t147 + (-t106 * t153 - t107 * t154 - t137 * t141 - t138 * t142 + t194 * t64 - t195 * t65) * mrSges(5,3) + (-t194 * t142 - t154 * t324) * Ifges(5,2) + t196 * (mrSges(5,1) * t142 + mrSges(5,2) * t141) + (t153 * t324 - t154 * t322 - t280 + t281) * Ifges(5,4) + t29 * t17 + t30 * t18 + (t261 * t166 + (-t204 * t200 + t262 - t261 * t174 + (t264 * mrSges(4,1) - 0.3e1 / 0.2e1 * t308 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t223) * qJD(1) + (qJD(1) * (-mrSges(5,1) * t194 + mrSges(5,2) * t195) + 0.2e1 * t319 + t144) * pkin(3) - t257) * qJD(3)) * t219 + m(7) * (t16 * t44 + t19 * t5 + t2 * t30 + t20 * t4 + t29 * t3 - t330) + m(6) * (t123 * t159 + t134 * t148 + t15 * t49 - t16 * t48 + t77 * t8 - t330) + (t309 / 0.2e1 + t337) * t80 + t134 * t93 + t153 * t126 / 0.2e1 - t154 * t125 / 0.2e1 + t159 * (mrSges(6,1) * t72 + mrSges(6,2) * t71) + t96 * t163 + t97 * t164 - t290 * t16 + (-t300 / 0.2e1 - t227) * t81 + (t195 * t141 + t153 * t322) * Ifges(5,1); t153 * t163 - t154 * t164 - t290 * t81 - (t71 * mrSges(6,3) + t14) * t241 + (-t280 - t281) * mrSges(5,3) + t238 * t80 + (-t219 * t198 + t223 * t200 + (-t219 ^ 2 - t223 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t72 * mrSges(6,3) + (-t216 * t83 - t220 * t84) * qJD(6) + t247) * t147 + m(4) * (t165 * t219 - t166 * t223 + (t174 * t223 - t219 * t341) * qJD(3)) + m(5) * (-t106 * t154 + t107 * t153 + t194 * t65 + t195 * t64) + m(6) * (t147 * t8 - t48 * t81 + t49 * t80 - t315) + m(7) * (t147 * t231 + t340 * t80 + t44 * t81 - t315); ((t163 * t222 - t164 * t218) * qJD(4) + (-t141 * t222 - t142 * t218) * mrSges(5,3)) * pkin(3) + t107 * t293 + (t139 * t83 + t179 * t18 + (-t19 * mrSges(7,3) - t179 * t84) * qJD(6)) * t220 + t343 * t121 - t25 * t83 - t24 * t84 + t224 - t162 * t93 - t111 * t163 - t110 * t164 - t165 * mrSges(4,2) - t166 * mrSges(4,1) + (-t139 * t84 + (-qJD(6) * t83 - t17) * t179 + (-t3 - t286) * mrSges(7,3)) * t216 + t178 * t14 + (-t184 * t71 - t185 * t72 + t295) * mrSges(6,3) + t174 * t198 - t341 * t200 + ((t263 - t292 - t188 / 0.2e1 - t208 / 0.2e1 + t341 * mrSges(4,3)) * t223 + (t262 + t174 * mrSges(4,3) + (t308 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t223) * qJD(1) + (-t144 - t319) * pkin(3) + t257) * t219) * qJD(1) + t290 * t344 + (-t148 * t162 - t184 * t9 + t185 * t8 + t343 * t49 + t344 * t48) * m(6) + (-t106 * t110 - t107 * t111 + (t218 * t64 + t222 * t65 + (-t106 * t218 + t107 * t222) * qJD(4)) * pkin(3)) * m(5) + (t139 * t340 + t178 * t9 + t231 * t179 - t19 * t24 - t20 * t25 - t344 * t44) * m(7); -t52 * t121 + (t164 + t293) * t107 + t226 * (pkin(4) * t217 + pkin(10)) + (-t190 * t93 + (-t217 * t72 - t221 * t71) * mrSges(6,3) + (-t290 * t217 + t238 * t221 + m(7) * (t217 * t44 + t221 * t340)) * qJD(5) + (0.2e1 * t148 * t323 + t217 * t8 - t221 * t9 + (-t217 * t48 + t221 * t49) * qJD(5)) * m(6)) * pkin(4) + t232 + t290 * t51 - m(6) * (-t48 * t51 + t49 * t52) - t23 * t83 - t22 * t84 + t224 - m(7) * (t19 * t22 + t20 * t23 + t44 * t51) + mrSges(6,3) * t295 - t106 * t163 + t347 * (-pkin(4) * t221 - pkin(5)); -t48 * t121 + t232 + t290 * t49 - m(7) * (t19 * t26 + t20 * t27 + t44 * t49) + t229 - t27 * t83 - t26 * t84 + (t312 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t242 - t337) * t258 + t226 * pkin(10) + (t49 * mrSges(6,3) + t227) * t242 - t347 * pkin(5); t68 - t44 * (mrSges(7,1) * t118 + mrSges(7,2) * t117) + (Ifges(7,1) * t117 - t305) * t328 + t58 * t327 + (Ifges(7,5) * t117 - Ifges(7,6) * t118) * t326 - t19 * t83 + t20 * t84 + (t117 * t19 + t118 * t20) * mrSges(7,3) + (-Ifges(7,2) * t118 + t116 + t59) * t329 + t339;];
tauc  = t1(:);
