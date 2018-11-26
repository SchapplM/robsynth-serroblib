% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2018-11-23 16:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:47:30
% EndTime: 2018-11-23 16:47:39
% DurationCPUTime: 8.97s
% Computational Cost: add. (5374->570), mult. (12785->748), div. (0->0), fcn. (7908->6), ass. (0->250)
t207 = sin(pkin(9));
t208 = cos(pkin(9));
t265 = t208 * qJD(2);
t212 = cos(qJ(2));
t271 = qJD(1) * t212;
t156 = t207 * t271 - t265;
t210 = sin(qJ(5));
t266 = t207 * qJD(2);
t221 = t208 * t271 + t266;
t296 = cos(qJ(5));
t101 = t210 * t156 - t296 * t221;
t211 = sin(qJ(2));
t223 = -t207 * t296 - t210 * t208;
t219 = t211 * t223;
t216 = qJD(2) * t219;
t62 = -qJD(1) * t216 + qJD(5) * t101;
t314 = t62 / 0.2e1;
t215 = -t156 * t296 - t210 * t221;
t250 = t296 * t208;
t242 = t211 * t250;
t229 = qJD(1) * t242;
t263 = qJD(1) * qJD(2);
t246 = t211 * t263;
t241 = t207 * t246;
t63 = -qJD(2) * t229 + qJD(5) * t215 + t210 * t241;
t313 = -t63 / 0.2e1;
t329 = Ifges(6,1) + Ifges(7,1);
t328 = Ifges(6,4) - Ifges(7,5);
t327 = Ifges(7,4) + Ifges(6,5);
t310 = pkin(3) + pkin(7);
t343 = t310 * t212;
t275 = t210 * t207;
t252 = t211 * t275;
t122 = -qJD(1) * t252 + t229;
t247 = qJD(5) * t296;
t268 = qJD(5) * t210;
t144 = -t207 * t268 + t208 * t247;
t274 = t122 + t144;
t123 = qJD(1) * t219;
t143 = -t207 * t247 - t208 * t268;
t273 = -t123 - t143;
t342 = Ifges(6,6) / 0.2e1;
t341 = -Ifges(7,6) / 0.2e1;
t245 = t212 * t263;
t340 = t328 * t313 + t329 * t314 + t327 * t245 / 0.2e1;
t339 = -mrSges(3,1) + mrSges(4,2);
t272 = qJD(1) * t211;
t196 = pkin(7) * t272;
t251 = -pkin(4) * t208 - pkin(3);
t127 = t251 * t272 - t196;
t338 = qJD(3) - t127;
t209 = -pkin(2) - qJ(4);
t243 = -t211 * qJ(3) - pkin(1);
t155 = t209 * t212 + t243;
t124 = t155 * qJD(1);
t163 = -pkin(3) * t272 - t196;
t316 = qJD(3) - t163;
t129 = qJD(2) * t209 + t316;
t73 = -t124 * t207 + t208 * t129;
t48 = pkin(4) * t272 + pkin(8) * t156 + t73;
t74 = t208 * t124 + t207 * t129;
t51 = -pkin(8) * t221 + t74;
t12 = -t210 * t51 + t296 * t48;
t13 = t210 * t48 + t296 * t51;
t278 = t207 * t211;
t228 = pkin(4) * t212 - pkin(8) * t278;
t220 = t228 * qJD(2);
t134 = (qJD(1) * t343 - qJD(4)) * qJD(2);
t191 = pkin(2) * t246;
t230 = -qJ(3) * t212 + qJ(4) * t211;
t264 = t211 * qJD(3);
t214 = qJD(2) * t230 - qJD(4) * t212 - t264;
t98 = qJD(1) * t214 + t191;
t60 = t208 * t134 - t207 * t98;
t37 = qJD(1) * t220 + t60;
t240 = t208 * t246;
t61 = t207 * t134 + t208 * t98;
t43 = pkin(8) * t240 + t61;
t3 = t210 * t37 + t48 * t247 - t268 * t51 + t296 * t43;
t317 = t250 - t275;
t4 = -qJD(5) * t13 - t210 * t43 + t296 * t37;
t337 = t12 * t273 - t13 * t274 + t223 * t3 - t317 * t4;
t192 = qJD(5) + t272;
t1 = qJ(6) * t245 + t192 * qJD(6) + t3;
t318 = qJD(6) - t12;
t10 = -t192 * pkin(5) + t318;
t11 = t192 * qJ(6) + t13;
t2 = -pkin(5) * t245 - t4;
t336 = t1 * t223 - t10 * t273 - t11 * t274 + t2 * t317;
t335 = -Ifges(7,5) * t62 / 0.2e1 + t245 * t341 + Ifges(7,3) * t313;
t334 = Ifges(6,4) * t314 + Ifges(6,2) * t313 + t245 * t342;
t306 = -t101 / 0.2e1;
t307 = t101 / 0.2e1;
t299 = t192 / 0.2e1;
t304 = t215 / 0.2e1;
t332 = -qJD(1) / 0.2e1;
t331 = -qJD(2) / 0.2e1;
t330 = qJD(2) / 0.2e1;
t326 = Ifges(6,6) - Ifges(7,6);
t281 = Ifges(7,5) * t101;
t97 = Ifges(6,4) * t101;
t324 = t192 * t327 + t215 * t329 - t281 + t97;
t323 = pkin(5) * t274 + qJ(6) * t273 - qJD(6) * t317 + t338;
t181 = t310 * t211;
t160 = t208 * t181;
t84 = t211 * pkin(4) + t160 + (pkin(8) * t212 - t155) * t207;
t104 = t208 * t155 + t207 * t181;
t276 = t208 * t212;
t88 = -pkin(8) * t276 + t104;
t290 = t210 * t84 + t296 * t88;
t270 = qJD(2) * t211;
t200 = pkin(2) * t270;
t113 = t200 + t214;
t166 = qJD(2) * t343;
t76 = -t207 * t113 + t208 * t166;
t59 = t220 + t76;
t277 = t208 * t211;
t262 = pkin(8) * t277;
t77 = t208 * t113 + t207 * t166;
t66 = qJD(2) * t262 + t77;
t8 = -qJD(5) * t290 - t210 * t66 + t296 * t59;
t172 = -pkin(2) * t212 + t243;
t153 = t172 * qJD(1);
t171 = -qJD(2) * pkin(2) + qJD(3) + t196;
t195 = Ifges(3,4) * t271;
t256 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t257 = t342 + t341;
t258 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t315 = (m(4) * t171 + (mrSges(4,1) + mrSges(3,3)) * t272 + t339 * qJD(2)) * pkin(7) + t257 * t101 + t256 * t192 + t258 * t215 + t11 * mrSges(7,3) + t12 * mrSges(6,1) + t171 * mrSges(4,1) + t73 * mrSges(5,1) + Ifges(3,5) * t330 + t195 / 0.2e1 + Ifges(4,4) * t331 + (-Ifges(4,2) * t211 - Ifges(4,6) * t212) * t332 + Ifges(6,6) * t307 + Ifges(7,6) * t306 - Ifges(5,6) * t221 / 0.2e1 - Ifges(5,5) * t156 - t10 * mrSges(7,1) - t13 * mrSges(6,2) - t153 * mrSges(4,3) - t74 * mrSges(5,2) + (Ifges(3,1) + Ifges(5,3)) * t272 / 0.2e1 + t327 * t304 + (Ifges(6,3) + Ifges(7,2)) * t299;
t312 = t63 / 0.2e1;
t309 = pkin(1) * mrSges(3,1);
t308 = pkin(1) * mrSges(3,2);
t305 = -t215 / 0.2e1;
t233 = Ifges(5,4) * t207 + Ifges(5,2) * t208;
t303 = -(Ifges(5,6) * t212 + t211 * t233) * t263 / 0.2e1;
t300 = -t192 / 0.2e1;
t298 = -t207 / 0.2e1;
t297 = t208 / 0.2e1;
t295 = -pkin(8) + t209;
t44 = -t63 * mrSges(7,2) + mrSges(7,3) * t245;
t47 = -mrSges(6,2) * t245 - t63 * mrSges(6,3);
t294 = t44 + t47;
t45 = mrSges(6,1) * t245 - t62 * mrSges(6,3);
t46 = -mrSges(7,1) * t245 + t62 * mrSges(7,2);
t293 = t46 - t45;
t197 = pkin(2) * t272;
t135 = qJD(1) * t230 + t197;
t198 = pkin(7) * t271;
t164 = pkin(3) * t271 + t198;
t90 = -t207 * t135 + t208 * t164;
t69 = qJD(1) * t228 + t90;
t91 = t208 * t135 + t207 * t164;
t78 = qJD(1) * t262 + t91;
t24 = t210 * t69 + t296 * t78;
t289 = mrSges(6,3) * t101;
t80 = -mrSges(6,2) * t192 + t289;
t81 = mrSges(7,2) * t101 + mrSges(7,3) * t192;
t292 = t80 + t81;
t288 = mrSges(6,3) * t215;
t82 = mrSges(6,1) * t192 - t288;
t83 = -mrSges(7,1) * t192 + mrSges(7,2) * t215;
t291 = t83 - t82;
t287 = Ifges(5,1) * t156;
t286 = Ifges(5,1) * t207;
t285 = Ifges(5,4) * t208;
t284 = Ifges(6,4) * t215;
t282 = Ifges(5,5) * t207;
t280 = Ifges(5,6) * t208;
t279 = qJD(2) * mrSges(3,2);
t193 = t207 * pkin(4) + qJ(3);
t269 = qJD(2) * t212;
t267 = qJD(5) * t212;
t261 = -0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(4,6);
t260 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t259 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t107 = mrSges(5,1) * t221 - t156 * mrSges(5,2);
t179 = -mrSges(4,1) * t271 - qJD(2) * mrSges(4,3);
t42 = -mrSges(6,1) * t101 + mrSges(6,2) * t215;
t255 = -t179 + t107 + t42;
t254 = m(4) * pkin(7) + mrSges(4,1);
t142 = pkin(4) * t276 + t343;
t20 = t63 * mrSges(6,1) + t62 * mrSges(6,2);
t19 = t63 * mrSges(7,1) - t62 * mrSges(7,3);
t244 = -t272 / 0.2e1;
t165 = t310 * t270;
t206 = qJD(2) * qJ(3);
t141 = qJD(4) + t206 + t164;
t177 = -t198 - t206;
t238 = m(4) * t177 - mrSges(3,3) * t271 + t179 + t279;
t235 = mrSges(5,1) * t208 - mrSges(5,2) * t207;
t234 = t285 + t286;
t232 = t207 * t61 + t208 * t60;
t231 = t207 * t73 - t208 * t74;
t23 = -t210 * t78 + t296 * t69;
t28 = -t210 * t88 + t296 * t84;
t7 = t210 * t59 + t84 * t247 - t268 * t88 + t296 * t66;
t167 = t295 * t207;
t168 = t295 * t208;
t109 = t167 * t296 + t210 * t168;
t224 = -t210 * t167 + t168 * t296;
t128 = (-pkin(7) + t251) * t270;
t222 = -qJ(3) * t269 - t264;
t126 = -mrSges(5,1) * t240 + mrSges(5,2) * t241;
t218 = t153 * mrSges(4,2) + Ifges(3,6) * t330 + (Ifges(3,4) * t211 + t212 * Ifges(3,2)) * qJD(1) / 0.2e1 + Ifges(4,5) * t331 + (-Ifges(4,6) * t211 - t212 * Ifges(4,3)) * t332 - t177 * mrSges(4,1);
t217 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t205 = qJD(2) * qJD(3);
t110 = qJD(1) * t128 + t205;
t105 = pkin(4) * t221 + t141;
t189 = Ifges(7,2) * t245;
t188 = Ifges(6,3) * t245;
t170 = pkin(7) * t246 - t205;
t161 = (mrSges(4,2) * t212 - t211 * mrSges(4,3)) * qJD(1);
t138 = t200 + t222;
t137 = (-mrSges(5,2) * t212 + mrSges(5,3) * t277) * t263;
t136 = (mrSges(5,1) * t212 - mrSges(5,3) * t278) * t263;
t133 = -qJD(1) * t165 + t205;
t131 = t223 * t212;
t130 = t317 * t212;
t125 = qJD(1) * t222 + t191;
t121 = mrSges(5,1) * t272 + mrSges(5,3) * t156;
t120 = -mrSges(5,2) * t272 - mrSges(5,3) * t221;
t112 = (Ifges(5,5) * t212 + t211 * t234) * t263;
t103 = -t207 * t155 + t160;
t96 = Ifges(7,5) * t215;
t95 = -pkin(5) * t223 - qJ(6) * t317 + t193;
t94 = -Ifges(5,4) * t221 + Ifges(5,5) * t272 - t287;
t93 = -Ifges(5,4) * t156 - Ifges(5,2) * t221 + Ifges(5,6) * t272;
t87 = -t223 * t267 + (t242 - t252) * qJD(2);
t86 = -t267 * t317 - t216;
t72 = qJD(4) * t317 + qJD(5) * t109;
t71 = qJD(4) * t223 + qJD(5) * t224;
t65 = pkin(5) * t130 - qJ(6) * t131 + t142;
t58 = Ifges(7,4) * t62;
t57 = Ifges(6,5) * t62;
t56 = Ifges(6,6) * t63;
t55 = Ifges(7,6) * t63;
t41 = -mrSges(7,1) * t101 - mrSges(7,3) * t215;
t40 = pkin(5) * t215 - qJ(6) * t101;
t33 = Ifges(6,2) * t101 + Ifges(6,6) * t192 + t284;
t30 = Ifges(7,6) * t192 - Ifges(7,3) * t101 + t96;
t27 = -t211 * pkin(5) - t28;
t26 = qJ(6) * t211 + t290;
t25 = -pkin(5) * t101 - qJ(6) * t215 + t105;
t22 = -pkin(5) * t271 - t23;
t21 = qJ(6) * t271 + t24;
t18 = -pkin(5) * t87 - qJ(6) * t86 - qJD(6) * t131 + t128;
t9 = pkin(5) * t63 - qJ(6) * t62 - qJD(6) * t215 + t110;
t6 = -pkin(5) * t269 - t8;
t5 = qJ(6) * t269 + t211 * qJD(6) + t7;
t14 = [t110 * (mrSges(6,1) * t130 + mrSges(6,2) * t131) + t9 * (mrSges(7,1) * t130 - mrSges(7,3) * t131) + (Ifges(7,5) * t131 + Ifges(7,3) * t130) * t312 + (Ifges(6,4) * t131 - Ifges(6,2) * t130) * t313 + (-t130 * t328 + t131 * t329) * t314 + (t11 * mrSges(7,2) + t13 * mrSges(6,3) - t105 * mrSges(6,1) - Ifges(7,3) * t306 + Ifges(6,2) * t307 + t33 / 0.2e1 - t30 / 0.2e1 - t25 * mrSges(7,1) + t326 * t299 + t328 * t304) * t87 - t130 * t334 - t130 * t335 + (-t1 * t130 + t131 * t2) * mrSges(7,2) + (t125 * mrSges(4,2) + t133 * t235 + t112 * t298 + t208 * t303 - t254 * t170 + (t207 * t60 - t208 * t61) * mrSges(5,3) + (((-t280 - t282 / 0.2e1 - t261) * t212 - t172 * mrSges(4,3) - 0.2e1 * t308 + t258 * t131 - t257 * t130) * qJD(1) + (Ifges(5,6) * t298 + t260) * qJD(2) + t315) * qJD(2)) * t212 + (-t130 * t3 - t131 * t4) * mrSges(6,3) + m(4) * (t125 * t172 + t153 * t138) + m(5) * (t103 * t60 + t104 * t61 + t133 * t343 - t141 * t165 + t73 * t76 + t74 * t77) + t343 * t126 + (-t61 * mrSges(5,2) + t60 * mrSges(5,1) - t125 * mrSges(4,3) + t188 / 0.2e1 + t189 / 0.2e1 + t58 / 0.2e1 + t57 / 0.2e1 - t56 / 0.2e1 + t55 / 0.2e1 - t257 * t63 + t258 * t62 + (-t156 * t234 / 0.2e1 - t141 * t235 + (t233 * t298 + t259) * qJD(2) + t207 * t94 / 0.2e1 + t93 * t297 - t231 * mrSges(5,3) + t238 * pkin(7) + ((0.3e1 / 0.2e1 * t282 + 0.3e1 / 0.2e1 * t280 + t261) * t211 - t172 * mrSges(4,2) - 0.2e1 * t309 + (-Ifges(5,2) * t208 ^ 2 + 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t254 * pkin(7) + (-0.3e1 / 0.2e1 * t285 - t286 / 0.2e1) * t207 + t256) * t212) * qJD(1) - t218) * qJD(2) + t217) * t211 + m(7) * (t1 * t26 + t10 * t6 + t11 * t5 + t18 * t25 + t2 * t27 + t65 * t9) + m(6) * (t105 * t128 + t110 * t142 + t12 * t8 + t13 * t7 + t28 * t4 + t290 * t3) + t290 * t47 - t165 * t107 + t138 * t161 + t142 * t20 + t103 * t136 + t104 * t137 + t128 * t42 + t77 * t120 + t76 * t121 + t131 * t340 + (mrSges(6,2) * t105 + mrSges(7,2) * t10 - mrSges(6,3) * t12 - mrSges(7,3) * t25 + Ifges(6,4) * t307 + Ifges(7,5) * t306 + t324 / 0.2e1 + t327 * t299 + t329 * t304) * t86 + t18 * t41 + t26 * t44 + t28 * t45 + t27 * t46 + t65 * t19 + t7 * t80 + t5 * t81 + t8 * t82 + t6 * t83; (((t309 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t211) * qJD(1) + (t287 / 0.2e1 - t141 * mrSges(5,2) + t73 * mrSges(5,3) - t94 / 0.2e1 + Ifges(5,5) * t244) * t207 + (-qJ(3) * mrSges(4,1) + t259) * qJD(2) + (-t238 + t279) * pkin(7) + (-t74 * mrSges(5,3) + t141 * mrSges(5,1) - t93 / 0.2e1 + Ifges(5,1) * t266 / 0.2e1 + Ifges(5,6) * t244 + (t156 / 0.2e1 + t265 / 0.2e1) * Ifges(5,4)) * t208 + t218) * t211 + ((Ifges(5,5) * t297 - pkin(2) * mrSges(4,1) + t257 * t223 + t258 * t317 + (-m(4) * pkin(2) + t339) * pkin(7) + t260) * qJD(2) - t195 / 0.2e1 + (t233 * t297 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t272 + ((t280 / 0.2e1 - Ifges(4,6) / 0.2e1) * t212 + t308) * qJD(1) - t315) * t212) * qJD(1) + t110 * (-mrSges(6,1) * t223 + mrSges(6,2) * t317) + t9 * (-mrSges(7,1) * t223 - mrSges(7,3) * t317) + (Ifges(7,5) * t317 - Ifges(7,3) * t223) * t312 + (Ifges(6,4) * t317 + Ifges(6,2) * t223) * t313 + (t223 * t328 + t317 * t329) * t314 + t317 * t340 + (t224 * t4 + t109 * t3 + t110 * t193 + (-t24 + t71) * t13 + (-t23 - t72) * t12 + t338 * t105) * m(6) + t336 * mrSges(7,2) + t337 * mrSges(6,3) + (mrSges(7,1) * t274 + mrSges(7,3) * t273) * t25 + (mrSges(6,1) * t274 - mrSges(6,2) * t273) * t105 + t291 * t72 + t292 * t71 + t294 * t109 + m(4) * (-t170 * qJ(3) - t177 * qJD(3)) + t255 * qJD(3) - t293 * t224 + (t1 * t109 - t224 * t2 + t9 * t95 + t323 * t25 + (-t21 + t71) * t11 + (-t22 + t72) * t10) * m(7) + t223 * t334 + t223 * t335 + (Ifges(6,4) * t143 - Ifges(7,5) * t123 - Ifges(6,2) * t144 - Ifges(7,3) * t122) * t307 + (-Ifges(6,4) * t123 + Ifges(7,5) * t143 + Ifges(6,2) * t122 + Ifges(7,3) * t144) * t306 + t324 * (t143 / 0.2e1 + t123 / 0.2e1) + (t122 * t326 - t123 * t327) * t300 + (t122 * t328 - t123 * t329) * t305 + (t133 * mrSges(5,2) + t112 / 0.2e1 - t60 * mrSges(5,3) - qJD(4) * t121 + t209 * t136) * t208 + (t133 * mrSges(5,1) - t61 * mrSges(5,3) - qJD(4) * t120 + t209 * t137 + t303) * t207 + (-m(4) * t153 - t161) * (-qJ(3) * t271 + t197) + (t30 - t33) * (t144 / 0.2e1 + t122 / 0.2e1) + t193 * t20 - t170 * mrSges(4,3) - t163 * t107 + qJ(3) * t126 - t127 * t42 - t91 * t120 - t90 * t121 + (qJ(3) * t133 + t232 * t209 + (-t207 * t74 - t208 * t73) * qJD(4) - t73 * t90 - t74 * t91 + t316 * t141) * m(5) + t323 * t41 + (t143 * t327 - t144 * t326) * t299 + (t143 * t329 - t144 * t328) * t304 - t24 * t80 - t21 * t81 - t23 * t82 - t22 * t83 + t95 * t19; t208 * t136 + t207 * t137 - t294 * t223 - t293 * t317 + (t120 * t208 - t121 * t207 + t161) * t272 + (t254 * t271 - t255 - t41) * qJD(2) - m(4) * (-qJD(2) * t177 - t153 * t272) + t274 * t292 + t273 * t291 + (-qJD(2) * t25 - t336) * m(7) + (-qJD(2) * t105 - t337) * m(6) + (-qJD(2) * t141 - t231 * t272 + t232) * m(5); -t292 * t101 + t221 * t120 - t156 * t121 - t291 * t215 + t126 + t19 + t20 + (-t10 * t215 - t101 * t11 + t9) * m(7) + (-t101 * t13 + t12 * t215 + t110) * m(6) + (-t73 * t156 + t221 * t74 + t133) * m(5); t217 + (-t10 * t101 + t11 * t215) * mrSges(7,2) + (t288 - t291) * t13 + (t289 - t292) * t12 + t188 + t189 + t58 + t57 - t56 + t55 - t105 * (mrSges(6,1) * t215 + mrSges(6,2) * t101) - t25 * (mrSges(7,1) * t215 - mrSges(7,3) * t101) + t33 * t304 + (Ifges(7,3) * t215 + t281) * t307 - t40 * t41 + qJ(6) * t44 - pkin(5) * t46 + qJD(6) * t81 + (t101 * t327 - t215 * t326) * t300 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t13 + t11 * t318 - t25 * t40) * m(7) + (-Ifges(6,2) * t215 + t324 + t97) * t306 + (t101 * t329 - t284 + t30 + t96) * t305; t215 * t41 - t192 * t81 + 0.2e1 * (t2 / 0.2e1 + t25 * t304 + t11 * t300) * m(7) + t46;];
tauc  = t14(:);
