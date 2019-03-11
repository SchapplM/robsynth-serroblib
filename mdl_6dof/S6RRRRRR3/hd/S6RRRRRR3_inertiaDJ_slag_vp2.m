% Calculate time derivative of joint inertia matrix for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:16
% EndTime: 2019-03-10 03:39:31
% DurationCPUTime: 6.58s
% Computational Cost: add. (18609->607), mult. (40609->863), div. (0->0), fcn. (39763->10), ass. (0->256)
t280 = sin(qJ(2));
t371 = -pkin(8) - pkin(7);
t261 = t371 * t280;
t285 = cos(qJ(2));
t263 = t371 * t285;
t279 = sin(qJ(3));
t284 = cos(qJ(3));
t214 = t261 * t279 - t263 * t284;
t312 = qJD(2) * t371;
t255 = t280 * t312;
t302 = t285 * t312;
t148 = t214 * qJD(3) + t255 * t279 - t284 * t302;
t244 = t279 * t285 + t284 * t280;
t278 = sin(qJ(4));
t327 = qJD(4) * t278;
t310 = t244 * t327;
t242 = t279 * t280 - t284 * t285;
t384 = qJD(2) + qJD(3);
t197 = t384 * t242;
t283 = cos(qJ(4));
t339 = t197 * t283;
t293 = t310 + t339;
t326 = qJD(4) * t283;
t340 = t197 * t278;
t294 = t244 * t326 - t340;
t90 = mrSges(5,1) * t294 - mrSges(5,2) * t293;
t392 = m(5) * t148 + t90;
t277 = sin(qJ(5));
t282 = cos(qJ(5));
t296 = t277 * t278 - t282 * t283;
t383 = qJD(4) + qJD(5);
t195 = t383 * t296;
t243 = t277 * t283 + t278 * t282;
t196 = t383 * t243;
t276 = sin(qJ(6));
t281 = cos(qJ(6));
t183 = -t243 * t276 - t281 * t296;
t98 = t183 * qJD(6) - t195 * t281 - t196 * t276;
t184 = t243 * t281 - t276 * t296;
t99 = -t184 * qJD(6) + t195 * t276 - t196 * t281;
t357 = Ifges(7,5) * t98 + Ifges(7,6) * t99;
t303 = -Ifges(6,5) * t195 - Ifges(6,6) * t196 + t357;
t391 = Ifges(5,5) * t326 + t303;
t385 = (t278 ^ 2 + t283 ^ 2) * t284;
t198 = t384 * t244;
t129 = pkin(2) * qJD(2) * t280 + pkin(3) * t198 + pkin(9) * t197;
t386 = t284 * t261 + t263 * t279;
t147 = qJD(3) * t386 + t284 * t255 + t279 * t302;
t269 = -pkin(2) * t285 - pkin(1);
t182 = t242 * pkin(3) - t244 * pkin(9) + t269;
t60 = t278 * t129 + t283 * t147 + t182 * t326 - t214 * t327;
t204 = t283 * t214;
t135 = t278 * t182 + t204;
t304 = t283 * t129 - t147 * t278;
t61 = -t135 * qJD(4) + t304;
t390 = -t61 * t278 + t283 * t60;
t265 = pkin(2) * t279 + pkin(9);
t358 = -pkin(10) - t265;
t306 = qJD(4) * t358;
t353 = pkin(2) * qJD(3);
t314 = t284 * t353;
t215 = t278 * t306 + t283 * t314;
t216 = -t278 * t314 + t283 * t306;
t238 = t358 * t278;
t273 = t283 * pkin(10);
t239 = t265 * t283 + t273;
t324 = qJD(5) * t282;
t325 = qJD(5) * t277;
t107 = t282 * t215 + t277 * t216 + t238 * t324 - t239 * t325;
t194 = t196 * pkin(11);
t82 = t107 - t194;
t177 = t277 * t238 + t282 * t239;
t108 = -qJD(5) * t177 - t215 * t277 + t282 * t216;
t363 = pkin(11) * t195;
t83 = t108 + t363;
t176 = t282 * t238 - t239 * t277;
t362 = pkin(11) * t243;
t153 = t176 - t362;
t236 = t296 * pkin(11);
t154 = -t236 + t177;
t91 = t153 * t281 - t154 * t276;
t28 = t91 * qJD(6) + t276 * t83 + t281 * t82;
t92 = t153 * t276 + t154 * t281;
t29 = -t92 * qJD(6) - t276 * t82 + t281 * t83;
t389 = t29 * mrSges(7,1) - t28 * mrSges(7,2);
t370 = -pkin(10) - pkin(9);
t311 = qJD(4) * t370;
t252 = t278 * t311;
t253 = t283 * t311;
t260 = t370 * t278;
t262 = pkin(9) * t283 + t273;
t145 = t282 * t252 + t277 * t253 + t260 * t324 - t262 * t325;
t111 = t145 - t194;
t213 = t277 * t260 + t282 * t262;
t146 = -qJD(5) * t213 - t252 * t277 + t282 * t253;
t112 = t146 + t363;
t211 = t282 * t260 - t262 * t277;
t162 = t211 - t362;
t163 = -t236 + t213;
t114 = t162 * t281 - t163 * t276;
t45 = t114 * qJD(6) + t111 * t281 + t112 * t276;
t115 = t162 * t276 + t163 * t281;
t46 = -t115 * qJD(6) - t111 * t276 + t112 * t281;
t388 = t46 * mrSges(7,1) - t45 * mrSges(7,2);
t170 = t296 * t244;
t134 = t283 * t182 - t214 * t278;
t336 = t244 * t283;
t105 = pkin(4) * t242 - pkin(10) * t336 + t134;
t337 = t244 * t278;
t118 = -pkin(10) * t337 + t135;
t63 = t277 * t105 + t282 * t118;
t387 = -Ifges(5,5) * t339 + Ifges(5,3) * t198;
t75 = -t196 * t244 + t296 * t197;
t76 = t170 * t383 + t243 * t197;
t382 = Ifges(6,5) * t75 + Ifges(6,6) * t76 + Ifges(6,3) * t198;
t381 = t108 * mrSges(6,1) - t107 * mrSges(6,2) + t389;
t380 = t146 * mrSges(6,1) - t145 * mrSges(6,2) + t388;
t379 = 0.2e1 * m(5);
t378 = 2 * m(6);
t377 = 2 * m(7);
t57 = -mrSges(7,1) * t99 + mrSges(7,2) * t98;
t376 = 0.2e1 * t57;
t130 = -mrSges(7,1) * t183 + mrSges(7,2) * t184;
t375 = 0.2e1 * t130;
t136 = mrSges(6,1) * t196 - mrSges(6,2) * t195;
t374 = 0.2e1 * t136;
t373 = 0.2e1 * t148;
t372 = 0.2e1 * t269;
t365 = pkin(2) * t284;
t201 = mrSges(6,1) * t296 + mrSges(6,2) * t243;
t364 = pkin(4) * t201;
t359 = t99 * mrSges(7,3);
t356 = Ifges(5,4) * t278;
t355 = Ifges(5,4) * t283;
t354 = Ifges(5,6) * t278;
t352 = pkin(4) * qJD(5);
t351 = pkin(5) * qJD(6);
t266 = pkin(4) * t282 + pkin(5);
t322 = qJD(6) * t281;
t323 = qJD(6) * t276;
t331 = t276 * t277;
t174 = t266 * t322 + (-t277 * t323 + (t281 * t282 - t331) * qJD(5)) * pkin(4);
t348 = t174 * mrSges(7,2);
t347 = t183 * mrSges(7,3);
t346 = t196 * mrSges(6,3);
t345 = t279 * mrSges(4,1);
t343 = t284 * mrSges(4,2);
t341 = t148 * t386;
t271 = pkin(4) * t327;
t272 = t279 * t353;
t254 = t272 + t271;
t335 = t254 * t201;
t330 = t277 * t281;
t257 = -mrSges(5,1) * t283 + mrSges(5,2) * t278;
t329 = t279 * t257;
t321 = 0.2e1 * mrSges(6,3);
t320 = 0.2e1 * mrSges(7,3);
t319 = 0.2e1 * t285;
t318 = mrSges(6,3) * t352;
t317 = mrSges(7,3) * t351;
t316 = t281 * t98 * mrSges(7,3);
t169 = t243 * t244;
t119 = -t169 * t281 + t170 * t276;
t25 = t119 * qJD(6) + t276 * t76 + t281 * t75;
t120 = -t169 * t276 - t170 * t281;
t26 = -t120 * qJD(6) - t276 * t75 + t281 * t76;
t315 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t198;
t313 = t282 * t195 * mrSges(6,3);
t268 = -pkin(4) * t283 - pkin(3);
t308 = -t327 / 0.2e1;
t307 = -(2 * Ifges(4,4)) - t354;
t168 = pkin(5) * t196 + t271;
t175 = -t266 * t323 + (-t277 * t322 + (-t276 * t282 - t330) * qJD(5)) * pkin(4);
t171 = t175 * mrSges(7,1);
t305 = t171 - t348;
t62 = t282 * t105 - t118 * t277;
t301 = mrSges(5,3) * t385;
t160 = pkin(4) * t337 - t386;
t300 = mrSges(5,1) * t278 + mrSges(5,2) * t283;
t299 = Ifges(5,1) * t283 - t356;
t298 = -Ifges(5,2) * t278 + t355;
t297 = Ifges(5,5) * t278 + Ifges(5,6) * t283;
t50 = pkin(5) * t242 + pkin(11) * t170 + t62;
t53 = -pkin(11) * t169 + t63;
t17 = -t276 * t53 + t281 * t50;
t18 = t276 * t50 + t281 * t53;
t218 = pkin(5) * t296 + t268;
t41 = pkin(10) * t339 + pkin(4) * t198 + (-t204 + (pkin(10) * t244 - t182) * t278) * qJD(4) + t304;
t49 = -pkin(10) * t294 + t60;
t15 = -qJD(5) * t63 - t277 * t49 + t282 * t41;
t5 = pkin(5) * t198 - pkin(11) * t75 + t15;
t14 = t105 * t324 - t118 * t325 + t277 * t41 + t282 * t49;
t6 = pkin(11) * t76 + t14;
t3 = qJD(6) * t17 + t276 * t5 + t281 * t6;
t4 = -qJD(6) * t18 - t276 * t6 + t281 * t5;
t295 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t315;
t292 = t281 * t183 * t317 + t303 + (pkin(5) * t359 + t184 * t317) * t276;
t291 = (-mrSges(6,1) * t277 - mrSges(6,2) * t282) * t352;
t131 = Ifges(7,4) * t184 + Ifges(7,2) * t183;
t132 = Ifges(7,1) * t184 + Ifges(7,4) * t183;
t137 = -Ifges(6,4) * t195 - Ifges(6,2) * t196;
t138 = -Ifges(6,1) * t195 - Ifges(6,4) * t196;
t202 = Ifges(6,4) * t243 - Ifges(6,2) * t296;
t203 = Ifges(6,1) * t243 - Ifges(6,4) * t296;
t250 = t298 * qJD(4);
t251 = t299 * qJD(4);
t259 = Ifges(5,1) * t278 + t355;
t58 = Ifges(7,4) * t98 + Ifges(7,2) * t99;
t59 = Ifges(7,1) * t98 + Ifges(7,4) * t99;
t290 = t99 * t131 + t98 * t132 - t137 * t296 + t243 * t138 + t183 * t58 + t184 * t59 - t195 * t203 - t196 * t202 + t283 * t250 + t278 * t251 + t259 * t326;
t289 = t15 * mrSges(6,1) - t14 * mrSges(6,2) + t295 + t382;
t87 = pkin(4) * t294 + t148;
t227 = -pkin(4) * t331 + t266 * t281;
t228 = pkin(4) * t330 + t266 * t276;
t288 = t174 * t347 - t282 * t296 * t318 + t228 * t359 + (-t175 * t184 - t227 * t98) * mrSges(7,3) + (-pkin(4) * t346 + t243 * t318) * t277 + t391;
t116 = mrSges(5,1) * t198 + mrSges(5,3) * t293;
t117 = -mrSges(5,2) * t198 - mrSges(5,3) * t294;
t185 = -mrSges(5,2) * t242 - mrSges(5,3) * t337;
t186 = mrSges(5,1) * t242 - mrSges(5,3) * t336;
t287 = m(5) * (-t134 * t326 - t135 * t327 + t390) + t283 * t117 - t278 * t116 - t185 * t327 - t186 * t326;
t10 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t198;
t109 = -Ifges(6,4) * t170 - Ifges(6,2) * t169 + Ifges(6,6) * t242;
t110 = -Ifges(6,1) * t170 - Ifges(6,4) * t169 + Ifges(6,5) * t242;
t121 = pkin(5) * t169 + t160;
t155 = Ifges(5,6) * t242 + t298 * t244;
t156 = Ifges(5,5) * t242 + t299 * t244;
t249 = t300 * qJD(4);
t258 = Ifges(5,2) * t283 + t356;
t32 = Ifges(6,4) * t75 + Ifges(6,2) * t76 + Ifges(6,6) * t198;
t33 = Ifges(6,1) * t75 + Ifges(6,4) * t76 + Ifges(6,5) * t198;
t47 = -pkin(5) * t76 + t87;
t66 = Ifges(7,4) * t120 + Ifges(7,2) * t119 + Ifges(7,6) * t242;
t67 = Ifges(7,1) * t120 + Ifges(7,4) * t119 + Ifges(7,5) * t242;
t80 = -Ifges(5,4) * t293 - Ifges(5,2) * t294 + Ifges(5,6) * t198;
t81 = -t293 * Ifges(5,1) - Ifges(5,4) * t294 + Ifges(5,5) * t198;
t9 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t198;
t286 = t18 * t359 - t63 * t346 + t3 * t347 + (-Ifges(5,6) * t327 + t391) * t242 / 0.2e1 + (t244 * t308 - t339 / 0.2e1) * t259 + (t257 - mrSges(4,1)) * t148 + (-t14 * t296 - t15 * t243 + t62 * t195) * mrSges(6,3) + ((-t134 * t283 - t135 * t278) * qJD(4) + t390) * mrSges(5,3) - t294 * t258 / 0.2e1 + (Ifges(6,5) * t243 + Ifges(7,5) * t184 - Ifges(6,6) * t296 + Ifges(7,6) * t183 + t297) * t198 / 0.2e1 - t296 * t32 / 0.2e1 + t283 * t80 / 0.2e1 + t278 * t81 / 0.2e1 + t243 * t33 / 0.2e1 + t76 * t202 / 0.2e1 + t75 * t203 / 0.2e1 - t196 * t109 / 0.2e1 - Ifges(4,5) * t197 - Ifges(4,6) * t198 + t87 * t201 - t195 * t110 / 0.2e1 + t183 * t9 / 0.2e1 + t184 * t10 / 0.2e1 - t169 * t137 / 0.2e1 - t170 * t138 / 0.2e1 + t160 * t136 - t147 * mrSges(4,2) + t47 * t130 + t26 * t131 / 0.2e1 + t25 * t132 / 0.2e1 + t119 * t58 / 0.2e1 + t120 * t59 / 0.2e1 + t121 * t57 + t98 * t67 / 0.2e1 + t99 * t66 / 0.2e1 + (-t17 * t98 - t4 * t184) * mrSges(7,3) - t386 * t249 + t155 * t308 + t156 * t326 / 0.2e1 + t251 * t336 / 0.2e1 - t250 * t337 / 0.2e1;
t267 = -pkin(3) - t365;
t256 = t268 - t365;
t237 = (-mrSges(7,1) * t276 - mrSges(7,2) * t281) * t351;
t217 = t218 - t365;
t179 = t300 * t244;
t159 = t168 + t272;
t152 = mrSges(6,1) * t242 + mrSges(6,3) * t170;
t151 = -mrSges(6,2) * t242 - mrSges(6,3) * t169;
t122 = mrSges(6,1) * t169 - mrSges(6,2) * t170;
t104 = mrSges(7,1) * t242 - mrSges(7,3) * t120;
t103 = -mrSges(7,2) * t242 + mrSges(7,3) * t119;
t72 = -mrSges(7,1) * t119 + mrSges(7,2) * t120;
t65 = -mrSges(6,2) * t198 + mrSges(6,3) * t76;
t64 = mrSges(6,1) * t198 - mrSges(6,3) * t75;
t35 = -mrSges(6,1) * t76 + mrSges(6,2) * t75;
t22 = -mrSges(7,2) * t198 + mrSges(7,3) * t26;
t21 = mrSges(7,1) * t198 - mrSges(7,3) * t25;
t11 = -mrSges(7,1) * t26 + mrSges(7,2) * t25;
t1 = [t179 * t373 + (t121 * t47 + t17 * t4 + t18 * t3) * t377 + (t14 * t63 + t15 * t62 + t160 * t87) * t378 - t156 * t339 + (mrSges(4,1) * t198 - mrSges(4,2) * t197) * t372 + (-0.2e1 * mrSges(4,3) * t147 - t307 * t197 + ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3) + Ifges(7,3)) * t198 + t315 + t382 + t387) * t242 + (mrSges(4,3) * t373 - 0.2e1 * Ifges(4,1) * t197 - t278 * t80 + t283 * t81 + (Ifges(5,5) * t283 + t307) * t198 + (-t283 * t155 - t278 * t156 - t242 * t297) * qJD(4)) * t244 + (t134 * t61 + t135 * t60 - t341) * t379 + 0.2e1 * m(4) * (t147 * t214 - t341) + t198 * (Ifges(7,5) * t120 + Ifges(7,6) * t119) + 0.2e1 * t60 * t185 + 0.2e1 * t61 * t186 - t169 * t32 - t170 * t33 + 0.2e1 * t160 * t35 + 0.2e1 * t14 * t151 + 0.2e1 * t15 * t152 + 0.2e1 * t135 * t117 + 0.2e1 * t134 * t116 + t119 * t9 + t120 * t10 + 0.2e1 * t121 * t11 + 0.2e1 * t87 * t122 + t76 * t109 + t75 * t110 + 0.2e1 * t3 * t103 + 0.2e1 * t4 * t104 + 0.2e1 * t47 * t72 + 0.2e1 * t62 * t64 + 0.2e1 * t63 * t65 + t26 * t66 + t25 * t67 + 0.2e1 * t17 * t21 + 0.2e1 * t18 * t22 - 0.2e1 * t386 * t90 + 0.2e1 * (t197 * t386 - t198 * t214) * mrSges(4,3) + t198 * (-Ifges(6,5) * t170 - Ifges(6,6) * t169) + t155 * t340 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t285) * t319 + (m(4) * pkin(2) * t372 - 0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (mrSges(4,1) * t242 + mrSges(4,2) * t244) - 0.2e1 * Ifges(3,4) * t280 + (-Ifges(3,2) + Ifges(3,1)) * t319) * t280) * qJD(2); t286 + m(6) * (t107 * t63 + t108 * t62 + t14 * t177 + t15 * t176 + t160 * t254 + t256 * t87) + m(7) * (t121 * t159 + t17 * t29 + t18 * t28 + t217 * t47 + t3 * t92 + t4 * t91) + t287 * t265 + (m(4) * (t147 * t279 - t148 * t284) + (t284 * t197 - t279 * t198) * mrSges(4,3) + ((-t242 * mrSges(4,3) + t283 * t185 - t278 * t186 + m(5) * (-t134 * t278 + t135 * t283) + m(4) * t214) * t284 + (t244 * mrSges(4,3) + t179 - (m(5) + m(4)) * t386) * t279) * qJD(3)) * pkin(2) + t254 * t122 + t256 * t35 + t217 * t11 + t176 * t64 + t177 * t65 + t159 * t72 + t107 * t151 + t108 * t152 + t28 * t103 + t29 * t104 + t91 * t21 + t92 * t22 + (Ifges(3,5) * t285 - Ifges(3,6) * t280 + (-mrSges(3,1) * t285 + mrSges(3,2) * t280) * pkin(7)) * qJD(2) + t392 * t267; ((t265 * t385 + t267 * t279) * t379 - 0.2e1 * t345 + 0.2e1 * t329 - 0.2e1 * t343 + 0.2e1 * t301) * t353 + (-t107 * t296 - t108 * t243 + t176 * t195 - t177 * t196) * t321 + (t183 * t28 - t184 * t29 - t91 * t98 + t92 * t99) * t320 + t290 + (t107 * t177 + t108 * t176 + t254 * t256) * t378 + (t159 * t217 + t28 * t92 + t29 * t91) * t377 - t258 * t327 + 0.2e1 * t267 * t249 + 0.2e1 * t335 + t256 * t374 + t217 * t376 + t159 * t375; t286 + t287 * pkin(9) + m(6) * (t14 * t213 + t145 * t63 + t146 * t62 + t15 * t211 + t160 * t271 + t268 * t87) + m(7) * (t114 * t4 + t115 * t3 + t121 * t168 + t17 * t46 + t18 * t45 + t218 * t47) + t268 * t35 + t218 * t11 + t211 * t64 + t213 * t65 + t168 * t72 + t145 * t151 + t146 * t152 + t114 * t21 + t115 * t22 + t45 * t103 + t46 * t104 + t122 * t271 - t392 * pkin(3); (m(5) * (-pkin(3) * t279 + pkin(9) * t385) - t345 + t329 - t343 + t301) * t353 + t290 + (-t258 + t364) * t327 + t335 + m(6) * (t107 * t213 + t108 * t211 + t145 * t177 + t146 * t176 + t254 * t268 + t256 * t271) + m(7) * (t114 * t29 + t115 * t28 + t159 * t218 + t168 * t217 + t45 * t92 + t46 * t91) + ((-t108 - t146) * t243 - (t107 + t145) * t296 - (t177 + t213) * t196 - (-t176 - t211) * t195) * mrSges(6,3) + (t217 + t218) * t57 + ((t115 + t92) * t99 + (-t114 - t91) * t98 + (-t29 - t46) * t184 + (t28 + t45) * t183) * mrSges(7,3) + (t168 + t159) * t130 + (t268 + t256) * t136 + (t267 - pkin(3)) * t249; (-t114 * t98 + t115 * t99 + t183 * t45 - t184 * t46) * t320 + (-t145 * t296 - t146 * t243 + t195 * t211 - t196 * t213) * t321 + t290 + (t145 * t213 + t146 * t211 + t268 * t271) * t378 + (t114 * t46 + t115 * t45 + t168 * t218) * t377 + t268 * t374 - 0.2e1 * pkin(3) * t249 + t218 * t376 + t168 * t375 + (-t258 + 0.2e1 * t364) * t327; t289 - t294 * Ifges(5,6) + m(7) * (t17 * t175 + t174 * t18 + t227 * t4 + t228 * t3) + (m(6) * (t14 * t277 + t15 * t282 + t324 * t63 - t325 * t62) + t282 * t64 + t277 * t65 - t152 * t325 + t151 * t324) * pkin(4) + t227 * t21 + t228 * t22 + t174 * t103 + t175 * t104 - t60 * mrSges(5,2) + t61 * mrSges(5,1) - Ifges(5,5) * t310 + t387; t288 + (m(6) * (t107 * t277 + t108 * t282 - t176 * t325 + t177 * t324) + t313) * pkin(4) + m(7) * (t174 * t92 + t175 * t91 + t227 * t29 + t228 * t28) + (t257 * t265 - t354) * qJD(4) - t300 * t314 + t381; t288 + m(7) * (t114 * t175 + t115 * t174 + t227 * t46 + t228 * t45) + (m(6) * (t145 * t277 + t146 * t282 - t211 * t325 + t213 * t324) + t313) * pkin(4) + (pkin(9) * t257 - t354) * qJD(4) + t380; -0.2e1 * t348 + 0.2e1 * t171 + (t174 * t228 + t175 * t227) * t377 + 0.2e1 * t291; (m(7) * (-t17 * t323 + t18 * t322 + t276 * t3 + t281 * t4) + t103 * t322 + t276 * t22 - t104 * t323 + t281 * t21) * pkin(5) + t289; (m(7) * (t276 * t28 + t281 * t29 + t322 * t92 - t323 * t91) - t316) * pkin(5) + t292 + t381; (m(7) * (-t114 * t323 + t115 * t322 + t276 * t45 + t281 * t46) - t316) * pkin(5) + t292 + t380; t291 + (m(7) * (t174 * t276 + t175 * t281 - t227 * t323 + t228 * t322) - mrSges(7,2) * t322 - mrSges(7,1) * t323) * pkin(5) + t305; 0.2e1 * t237; t295; t357 + t389; t357 + t388; t305; t237; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
