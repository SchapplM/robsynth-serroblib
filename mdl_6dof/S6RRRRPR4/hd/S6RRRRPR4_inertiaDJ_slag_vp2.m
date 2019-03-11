% Calculate time derivative of joint inertia matrix for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:16
% EndTime: 2019-03-09 22:06:25
% DurationCPUTime: 5.65s
% Computational Cost: add. (13690->560), mult. (30445->792), div. (0->0), fcn. (29555->10), ass. (0->241)
t273 = sin(qJ(2));
t350 = -pkin(8) - pkin(7);
t253 = t350 * t273;
t277 = cos(qJ(2));
t254 = t350 * t277;
t272 = sin(qJ(3));
t276 = cos(qJ(3));
t202 = t253 * t272 - t254 * t276;
t302 = qJD(2) * t350;
t246 = t273 * t302;
t294 = t277 * t302;
t135 = qJD(3) * t202 + t246 * t272 - t276 * t294;
t240 = t272 * t277 + t276 * t273;
t271 = sin(qJ(4));
t311 = qJD(4) * t271;
t301 = t240 * t311;
t239 = t272 * t273 - t276 * t277;
t359 = qJD(2) + qJD(3);
t190 = t359 * t239;
t275 = cos(qJ(4));
t323 = t190 * t275;
t284 = t301 + t323;
t310 = qJD(4) * t275;
t300 = t240 * t310;
t324 = t190 * t271;
t285 = t300 - t324;
t88 = mrSges(5,1) * t285 - mrSges(5,2) * t284;
t368 = m(5) * t135 + t88;
t360 = (t271 ^ 2 + t275 ^ 2) * t276;
t367 = Ifges(5,3) + Ifges(6,3);
t191 = t359 * t240;
t123 = pkin(2) * qJD(2) * t273 + pkin(3) * t191 + pkin(9) * t190;
t361 = t276 * t253 + t254 * t272;
t134 = qJD(3) * t361 + t276 * t246 + t272 * t294;
t260 = -pkin(2) * t277 - pkin(1);
t172 = t239 * pkin(3) - t240 * pkin(9) + t260;
t303 = t271 * t123 + t275 * t134 + t172 * t310;
t54 = -t202 * t311 + t303;
t195 = t275 * t202;
t128 = t271 * t172 + t195;
t296 = t275 * t123 - t134 * t271;
t55 = -t128 * qJD(4) + t296;
t366 = -t55 * t271 + t275 * t54;
t268 = sin(pkin(11));
t269 = cos(pkin(11));
t235 = t268 * t275 + t269 * t271;
t225 = t235 * qJD(4);
t287 = t268 * t271 - t269 * t275;
t226 = t287 * qJD(4);
t270 = sin(qJ(6));
t274 = cos(qJ(6));
t173 = -t235 * t270 - t274 * t287;
t116 = qJD(6) * t173 - t225 * t270 - t226 * t274;
t174 = t235 * t274 - t270 * t287;
t117 = -qJD(6) * t174 - t225 * t274 + t226 * t270;
t313 = Ifges(7,5) * t116 + Ifges(7,6) * t117;
t365 = Ifges(5,5) * t310 - Ifges(6,5) * t226 - Ifges(6,6) * t225 + t313;
t160 = t287 * t240;
t321 = t240 * t271;
t107 = -qJ(5) * t321 + t128;
t127 = t275 * t172 - t202 * t271;
t265 = t275 * qJ(5);
t94 = pkin(4) * t239 - t240 * t265 + t127;
t56 = -t107 * t268 + t269 * t94;
t42 = pkin(5) * t239 + pkin(10) * t160 + t56;
t159 = t235 * t240;
t57 = t269 * t107 + t268 * t94;
t47 = -pkin(10) * t159 + t57;
t16 = -t270 * t47 + t274 * t42;
t286 = qJ(5) * t190 - qJD(5) * t240;
t29 = pkin(4) * t191 + t286 * t275 + (-t195 + (qJ(5) * t240 - t172) * t271) * qJD(4) + t296;
t37 = -qJ(5) * t300 + (-qJD(4) * t202 + t286) * t271 + t303;
t13 = -t268 * t37 + t269 * t29;
t82 = t190 * t287 - t225 * t240;
t5 = pkin(5) * t191 - pkin(10) * t82 + t13;
t14 = t268 * t29 + t269 * t37;
t81 = t190 * t235 + t226 * t240;
t8 = pkin(10) * t81 + t14;
t3 = qJD(6) * t16 + t270 * t5 + t274 * t8;
t17 = t270 * t42 + t274 * t47;
t4 = -qJD(6) * t17 - t270 * t8 + t274 * t5;
t364 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t264 = t275 * qJD(5);
t257 = pkin(2) * t272 + pkin(9);
t314 = -qJ(5) - t257;
t295 = qJD(4) * t314;
t332 = pkin(2) * qJD(3);
t305 = t276 * t332;
t188 = t271 * t295 + t275 * t305 + t264;
t189 = (-qJD(5) - t305) * t271 + t275 * t295;
t125 = -t188 * t268 + t269 * t189;
t341 = pkin(10) * t226;
t101 = t125 + t341;
t126 = t269 * t188 + t268 * t189;
t219 = t225 * pkin(10);
t102 = -t219 + t126;
t231 = t314 * t271;
t232 = t257 * t275 + t265;
t164 = t269 * t231 - t232 * t268;
t340 = pkin(10) * t235;
t141 = t164 - t340;
t165 = t268 * t231 + t269 * t232;
t230 = t287 * pkin(10);
t142 = -t230 + t165;
t86 = t141 * t274 - t142 * t270;
t32 = qJD(6) * t86 + t101 * t270 + t102 * t274;
t87 = t141 * t270 + t142 * t274;
t33 = -qJD(6) * t87 + t101 * t274 - t102 * t270;
t363 = t33 * mrSges(7,1) - t32 * mrSges(7,2);
t337 = -qJ(5) - pkin(9);
t297 = qJD(4) * t337;
t222 = t271 * t297 + t264;
t223 = -qJD(5) * t271 + t275 * t297;
t153 = -t222 * t268 + t269 * t223;
t132 = t153 + t341;
t154 = t269 * t222 + t268 * t223;
t133 = -t219 + t154;
t248 = t337 * t271;
t250 = pkin(9) * t275 + t265;
t198 = t269 * t248 - t250 * t268;
t155 = t198 - t340;
t199 = t268 * t248 + t269 * t250;
t156 = -t230 + t199;
t98 = t155 * t274 - t156 * t270;
t50 = qJD(6) * t98 + t132 * t270 + t133 * t274;
t99 = t155 * t270 + t156 * t274;
t51 = -qJD(6) * t99 + t132 * t274 - t133 * t270;
t362 = t51 * mrSges(7,1) - t50 * mrSges(7,2);
t358 = 0.2e1 * m(5);
t357 = 2 * m(6);
t356 = 2 * m(7);
t64 = -t117 * mrSges(7,1) + t116 * mrSges(7,2);
t355 = 0.2e1 * t64;
t120 = -mrSges(7,1) * t173 + mrSges(7,2) * t174;
t354 = 0.2e1 * t120;
t353 = 0.2e1 * t135;
t166 = t225 * mrSges(6,1) - t226 * mrSges(6,2);
t352 = 0.2e1 * t166;
t351 = 0.2e1 * t260;
t344 = pkin(2) * t276;
t185 = mrSges(6,1) * t287 + mrSges(6,2) * t235;
t343 = pkin(4) * t185;
t342 = pkin(4) * t268;
t336 = mrSges(6,3) * t225;
t335 = Ifges(5,4) * t271;
t334 = Ifges(5,4) * t275;
t333 = Ifges(5,6) * t271;
t331 = t117 * mrSges(7,3);
t330 = t173 * mrSges(7,3);
t329 = t272 * mrSges(4,1);
t327 = t276 * mrSges(4,2);
t325 = t135 * t361;
t320 = t240 * t275;
t262 = pkin(4) * t311;
t263 = t272 * t332;
t245 = t263 + t262;
t319 = t245 * t185;
t249 = -mrSges(5,1) * t275 + mrSges(5,2) * t271;
t315 = t272 * t249;
t309 = 0.2e1 * mrSges(6,3);
t308 = 0.2e1 * mrSges(7,3);
t307 = 0.2e1 * t277;
t105 = -t159 * t274 + t160 * t270;
t25 = qJD(6) * t105 + t270 * t81 + t274 * t82;
t106 = -t159 * t270 - t160 * t274;
t26 = -qJD(6) * t106 - t270 * t82 + t274 * t81;
t306 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t191;
t304 = t269 * t226 * mrSges(6,3);
t259 = -pkin(4) * t275 - pkin(3);
t46 = -t81 * mrSges(6,1) + t82 * mrSges(6,2);
t11 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t299 = -t311 / 0.2e1;
t298 = -(2 * Ifges(4,4)) - t333;
t204 = pkin(5) * t225 + t262;
t293 = mrSges(5,3) * t360;
t152 = pkin(4) * t321 - t361;
t292 = mrSges(5,1) * t271 + mrSges(5,2) * t275;
t255 = pkin(4) * t269 + pkin(5);
t220 = t255 * t274 - t270 * t342;
t210 = t220 * qJD(6);
t221 = t255 * t270 + t274 * t342;
t211 = t221 * qJD(6);
t291 = -t211 * mrSges(7,1) - t210 * mrSges(7,2);
t290 = Ifges(5,1) * t275 - t335;
t289 = -Ifges(5,2) * t271 + t334;
t288 = Ifges(5,5) * t271 + Ifges(5,6) * t275;
t206 = pkin(5) * t287 + t259;
t283 = t166 + t64;
t282 = -Ifges(5,5) * t323 + Ifges(6,5) * t82 + Ifges(6,6) * t81 + t191 * t367 + t306;
t121 = Ifges(7,4) * t174 + Ifges(7,2) * t173;
t122 = Ifges(7,1) * t174 + Ifges(7,4) * t173;
t167 = -Ifges(6,4) * t226 - Ifges(6,2) * t225;
t168 = -Ifges(6,1) * t226 - Ifges(6,4) * t225;
t186 = Ifges(6,4) * t235 - Ifges(6,2) * t287;
t187 = Ifges(6,1) * t235 - Ifges(6,4) * t287;
t243 = t289 * qJD(4);
t244 = t290 * qJD(4);
t252 = Ifges(5,1) * t271 + t334;
t65 = Ifges(7,4) * t116 + Ifges(7,2) * t117;
t66 = Ifges(7,1) * t116 + Ifges(7,4) * t117;
t281 = t116 * t122 + t117 * t121 - t167 * t287 + t235 * t168 + t173 * t65 + t174 * t66 - t225 * t186 - t226 * t187 + t275 * t243 + t271 * t244 + t252 * t310;
t85 = pkin(4) * t285 + t135;
t280 = t210 * t330 - t336 * t342 + t221 * t331 + (-t116 * t220 + t174 * t211) * mrSges(7,3) + t365;
t103 = mrSges(5,1) * t191 + mrSges(5,3) * t284;
t104 = -mrSges(5,2) * t191 - mrSges(5,3) * t285;
t178 = -mrSges(5,2) * t239 - mrSges(5,3) * t321;
t179 = mrSges(5,1) * t239 - mrSges(5,3) * t320;
t279 = -t271 * t103 + t275 * t104 + m(5) * (-t127 * t310 - t128 * t311 + t366) - t178 * t311 - t179 * t310;
t10 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t191;
t108 = pkin(5) * t159 + t152;
t146 = Ifges(5,6) * t239 + t240 * t289;
t147 = Ifges(5,5) * t239 + t240 * t290;
t242 = t292 * qJD(4);
t251 = Ifges(5,2) * t275 + t335;
t40 = Ifges(6,4) * t82 + Ifges(6,2) * t81 + Ifges(6,6) * t191;
t41 = Ifges(6,1) * t82 + Ifges(6,4) * t81 + Ifges(6,5) * t191;
t48 = -pkin(5) * t81 + t85;
t60 = Ifges(7,4) * t106 + Ifges(7,2) * t105 + Ifges(7,6) * t239;
t61 = Ifges(7,1) * t106 + Ifges(7,4) * t105 + Ifges(7,5) * t239;
t76 = -Ifges(5,4) * t284 - Ifges(5,2) * t285 + Ifges(5,6) * t191;
t77 = -Ifges(5,1) * t284 - Ifges(5,4) * t285 + Ifges(5,5) * t191;
t9 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t191;
t95 = -Ifges(6,4) * t160 - Ifges(6,2) * t159 + Ifges(6,6) * t239;
t96 = -Ifges(6,1) * t160 - Ifges(6,4) * t159 + Ifges(6,5) * t239;
t278 = (t249 - mrSges(4,1)) * t135 + (-t116 * t16 - t174 * t4) * mrSges(7,3) + t146 * t299 - t361 * t242 + t275 * t76 / 0.2e1 + t271 * t77 / 0.2e1 + t235 * t41 / 0.2e1 - t225 * t95 / 0.2e1 - t226 * t96 / 0.2e1 + (t240 * t299 - t323 / 0.2e1) * t252 - Ifges(4,6) * t191 + (-t13 * t235 - t14 * t287 + t226 * t56) * mrSges(6,3) - t287 * t40 / 0.2e1 + (Ifges(6,5) * t235 + Ifges(7,5) * t174 - Ifges(6,6) * t287 + Ifges(7,6) * t173 + t288) * t191 / 0.2e1 + t85 * t185 + t81 * t186 / 0.2e1 + t82 * t187 / 0.2e1 - Ifges(4,5) * t190 - t159 * t167 / 0.2e1 - t160 * t168 / 0.2e1 + t173 * t9 / 0.2e1 + t174 * t10 / 0.2e1 + t152 * t166 - t134 * mrSges(4,2) + t116 * t61 / 0.2e1 + t117 * t60 / 0.2e1 + t48 * t120 + t26 * t121 / 0.2e1 + t25 * t122 / 0.2e1 + t108 * t64 + t105 * t65 / 0.2e1 + t106 * t66 / 0.2e1 + ((-t127 * t275 - t128 * t271) * qJD(4) + t366) * mrSges(5,3) - t285 * t251 / 0.2e1 + t3 * t330 + t17 * t331 - t57 * t336 + t147 * t310 / 0.2e1 + t244 * t320 / 0.2e1 - t243 * t321 / 0.2e1 + (-Ifges(5,6) * t311 + t365) * t239 / 0.2e1;
t258 = -pkin(3) - t344;
t247 = t259 - t344;
t203 = t206 - t344;
t192 = t204 + t263;
t169 = t292 * t240;
t138 = mrSges(6,1) * t239 + mrSges(6,3) * t160;
t137 = -mrSges(6,2) * t239 - mrSges(6,3) * t159;
t109 = mrSges(6,1) * t159 - mrSges(6,2) * t160;
t93 = mrSges(7,1) * t239 - mrSges(7,3) * t106;
t92 = -mrSges(7,2) * t239 + mrSges(7,3) * t105;
t68 = mrSges(6,1) * t191 - mrSges(6,3) * t82;
t67 = -mrSges(6,2) * t191 + mrSges(6,3) * t81;
t62 = -mrSges(7,1) * t105 + mrSges(7,2) * t106;
t21 = -mrSges(7,2) * t191 + mrSges(7,3) * t26;
t20 = mrSges(7,1) * t191 - mrSges(7,3) * t25;
t1 = [t191 * (Ifges(7,5) * t106 + Ifges(7,6) * t105) - 0.2e1 * t361 * t88 + 0.2e1 * (t190 * t361 - t191 * t202) * mrSges(4,3) + (mrSges(4,1) * t191 - mrSges(4,2) * t190) * t351 + (mrSges(4,3) * t353 - 0.2e1 * Ifges(4,1) * t190 - t271 * t76 + t275 * t77 + (Ifges(5,5) * t275 + t298) * t191 + (-t275 * t146 - t271 * t147 - t239 * t288) * qJD(4)) * t240 + (t127 * t55 + t128 * t54 - t325) * t358 + 0.2e1 * m(4) * (t134 * t202 - t325) + 0.2e1 * t54 * t178 + 0.2e1 * t55 * t179 - t159 * t40 - t160 * t41 + 0.2e1 * t152 * t46 + 0.2e1 * t14 * t137 + 0.2e1 * t13 * t138 + 0.2e1 * t127 * t103 + 0.2e1 * t128 * t104 + 0.2e1 * t108 * t11 + 0.2e1 * t85 * t109 + t105 * t9 + t106 * t10 + 0.2e1 * t3 * t92 + 0.2e1 * t4 * t93 + t81 * t95 + t82 * t96 + 0.2e1 * t57 * t67 + 0.2e1 * t56 * t68 + t26 * t60 + t25 * t61 + 0.2e1 * t48 * t62 + 0.2e1 * t16 * t20 + 0.2e1 * t17 * t21 + (-0.2e1 * t134 * mrSges(4,3) - t298 * t190 + ((2 * Ifges(4,2)) + Ifges(7,3) + t367) * t191 + t282) * t239 - t147 * t323 + t169 * t353 + (t108 * t48 + t16 * t4 + t17 * t3) * t356 + (t13 * t56 + t14 * t57 + t152 * t85) * t357 + t191 * (-Ifges(6,5) * t160 - Ifges(6,6) * t159) + t146 * t324 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t277) * t307 + (-0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t351 + 0.2e1 * pkin(2) * (mrSges(4,1) * t239 + mrSges(4,2) * t240) - 0.2e1 * Ifges(3,4) * t273 + (-Ifges(3,2) + Ifges(3,1)) * t307) * t273) * qJD(2); t279 * t257 + (m(4) * (t134 * t272 - t135 * t276) + (t190 * t276 - t191 * t272) * mrSges(4,3) + ((-t239 * mrSges(4,3) + t275 * t178 - t271 * t179 + m(4) * t202 + m(5) * (-t127 * t271 + t128 * t275)) * t276 + (t240 * mrSges(4,3) + t169 - (m(4) + m(5)) * t361) * t272) * qJD(3)) * pkin(2) + (Ifges(3,5) * t277 - Ifges(3,6) * t273 + (-mrSges(3,1) * t277 + mrSges(3,2) * t273) * pkin(7)) * qJD(2) + t245 * t109 + t247 * t46 + t278 + t203 * t11 + t192 * t62 + t164 * t68 + t165 * t67 + t126 * t137 + t125 * t138 + t87 * t21 + t32 * t92 + t33 * t93 + t86 * t20 + m(7) * (t108 * t192 + t16 * t33 + t17 * t32 + t203 * t48 + t3 * t87 + t4 * t86) + m(6) * (t125 * t56 + t126 * t57 + t13 * t164 + t14 * t165 + t152 * t245 + t247 * t85) + t368 * t258; (t192 * t203 + t32 * t87 + t33 * t86) * t356 + (t125 * t164 + t126 * t165 + t245 * t247) * t357 + ((t257 * t360 + t258 * t272) * t358 - 0.2e1 * t329 + 0.2e1 * t315 - 0.2e1 * t327 + 0.2e1 * t293) * t332 + (-t116 * t86 + t87 * t117 + t32 * t173 - t33 * t174) * t308 + (-t125 * t235 - t126 * t287 + t164 * t226 - t165 * t225) * t309 + t281 + 0.2e1 * t258 * t242 + 0.2e1 * t319 + t247 * t352 + t203 * t355 + t192 * t354 - t251 * t311; m(7) * (t108 * t204 + t16 * t51 + t17 * t50 + t206 * t48 + t3 * t99 + t4 * t98) + m(6) * (t13 * t198 + t14 * t199 + t152 * t262 + t153 * t56 + t154 * t57 + t259 * t85) + t279 * pkin(9) + t259 * t46 + t278 + t204 * t62 + t206 * t11 + t198 * t68 + t199 * t67 + t153 * t138 + t154 * t137 + t99 * t21 + t50 * t92 + t51 * t93 + t98 * t20 + t109 * t262 - t368 * pkin(3); (t259 + t247) * t166 + (t204 + t192) * t120 + (-t251 + t343) * t311 + (m(5) * (-pkin(3) * t272 + pkin(9) * t360) - t329 + t315 - t327 + t293) * t332 + ((-t125 - t153) * t235 - (t126 + t154) * t287 - (-t164 - t198) * t226 - (t165 + t199) * t225) * mrSges(6,3) + t281 + ((-t33 - t51) * t174 + (t32 + t50) * t173 + (t87 + t99) * t117 + (-t86 - t98) * t116) * mrSges(7,3) + t319 + m(7) * (t192 * t206 + t203 * t204 + t32 * t99 + t33 * t98 + t50 * t87 + t51 * t86) + m(6) * (t125 * t198 + t126 * t199 + t153 * t164 + t154 * t165 + t245 * t259 + t247 * t262) + (t258 - pkin(3)) * t242 + (t203 + t206) * t64; (t204 * t206 + t50 * t99 + t51 * t98) * t356 + (t153 * t198 + t154 * t199 + t259 * t262) * t357 + (-t116 * t98 + t99 * t117 + t50 * t173 - t51 * t174) * t308 + (-t153 * t235 - t154 * t287 + t198 * t226 - t199 * t225) * t309 + t281 + t259 * t352 - 0.2e1 * pkin(3) * t242 + (-t251 + 0.2e1 * t343) * t311 + t204 * t354 + t206 * t355; m(7) * (-t16 * t211 + t17 * t210 + t220 * t4 + t221 * t3) + (t268 * t67 + t269 * t68 + m(6) * (t13 * t269 + t14 * t268)) * pkin(4) - t285 * Ifges(5,6) + t220 * t20 + t221 * t21 + t210 * t92 - t211 * t93 - t54 * mrSges(5,2) + t55 * mrSges(5,1) - t14 * mrSges(6,2) + t13 * mrSges(6,1) + t282 - Ifges(5,5) * t301 + t364; t280 + m(7) * (t210 * t87 - t211 * t86 + t220 * t33 + t221 * t32) + (m(6) * (t125 * t269 + t126 * t268) + t304) * pkin(4) - t292 * t305 + (t249 * t257 - t333) * qJD(4) + t125 * mrSges(6,1) - t126 * mrSges(6,2) + t363; t280 + (m(6) * (t153 * t269 + t154 * t268) + t304) * pkin(4) + m(7) * (t210 * t99 - t211 * t98 + t220 * t51 + t221 * t50) + (pkin(9) * t249 - t333) * qJD(4) + t153 * mrSges(6,1) - t154 * mrSges(6,2) + t362; 0.2e1 * m(7) * (t210 * t221 - t211 * t220) + 0.2e1 * t291; m(6) * t85 + m(7) * t48 + t11 + t46; m(6) * t245 + m(7) * t192 + t283; m(6) * t262 + m(7) * t204 + t283; 0; 0; t306 + t364; t313 + t363; t313 + t362; t291; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
