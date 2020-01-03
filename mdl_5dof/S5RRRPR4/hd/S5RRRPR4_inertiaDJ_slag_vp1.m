% Calculate time derivative of joint inertia matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:15
% DurationCPUTime: 9.19s
% Computational Cost: add. (12403->519), mult. (15663->714), div. (0->0), fcn. (14380->8), ass. (0->269)
t248 = sin(qJ(3));
t251 = cos(qJ(3));
t359 = Icges(5,5) * t251;
t298 = Icges(5,3) * t248 + t359;
t361 = Icges(4,4) * t251;
t301 = -Icges(4,2) * t248 + t361;
t360 = Icges(5,5) * t248;
t302 = Icges(5,1) * t251 + t360;
t362 = Icges(4,4) * t248;
t303 = Icges(4,1) * t251 - t362;
t411 = (t302 + t303) * t248 - (t298 - t301) * t251;
t410 = t360 - t362 + (-Icges(4,2) - Icges(5,3)) * t251;
t409 = -t359 + t361 + (Icges(4,1) + Icges(5,1)) * t248;
t398 = t410 * t248 + t409 * t251;
t246 = qJ(1) + qJ(2);
t242 = sin(t246);
t247 = sin(qJ(5));
t250 = cos(qJ(5));
t283 = t247 * t251 - t248 * t250;
t168 = t283 * t242;
t282 = t247 * t248 + t250 * t251;
t169 = t282 * t242;
t243 = cos(t246);
t126 = Icges(6,4) * t169 - Icges(6,2) * t168 + Icges(6,6) * t243;
t128 = Icges(6,1) * t169 - Icges(6,4) * t168 + Icges(6,5) * t243;
t141 = -Icges(6,5) * t283 - Icges(6,6) * t282;
t142 = -Icges(6,4) * t283 - Icges(6,2) * t282;
t143 = -Icges(6,1) * t283 - Icges(6,4) * t282;
t403 = -t126 * t282 - t128 * t283 + t141 * t243 - t142 * t168 + t143 * t169;
t170 = t283 * t243;
t171 = t282 * t243;
t127 = Icges(6,4) * t171 - Icges(6,2) * t170 - Icges(6,6) * t242;
t129 = Icges(6,1) * t171 - Icges(6,4) * t170 - Icges(6,5) * t242;
t402 = -t127 * t282 - t129 * t283 - t141 * t242 - t142 * t170 + t143 * t171;
t390 = qJD(3) - qJD(5);
t401 = t390 * t283;
t210 = Icges(4,5) * t248 + Icges(4,6) * t251;
t211 = Icges(5,4) * t248 - Icges(5,6) * t251;
t400 = t210 + t211;
t332 = qJD(3) * t251;
t322 = t242 * t332;
t245 = qJD(1) + qJD(2);
t349 = t245 * t248;
t397 = t243 * t349 + t322;
t299 = Icges(4,5) * t251 - Icges(4,6) * t248;
t300 = Icges(5,4) * t251 + Icges(5,6) * t248;
t396 = t398 * t245 + (-t299 - t300) * qJD(3);
t333 = qJD(3) * t248;
t323 = t242 * t333;
t354 = t242 * t245;
t138 = t390 * t282;
t86 = t138 * t242 - t170 * t245;
t87 = t171 * t245 + t242 * t401;
t41 = t87 * rSges(6,1) + t86 * rSges(6,2) - rSges(6,3) * t354;
t395 = -pkin(4) * t323 - pkin(8) * t354 + t41;
t277 = t301 * t243;
t155 = Icges(4,6) * t242 + t277;
t279 = t303 * t243;
t159 = Icges(4,5) * t242 + t279;
t286 = t155 * t248 - t159 * t251;
t394 = t242 * t286;
t274 = t298 * t243;
t149 = Icges(5,6) * t242 + t274;
t278 = t302 * t243;
t157 = Icges(5,4) * t242 + t278;
t290 = t149 * t248 + t157 * t251;
t393 = t242 * t290;
t154 = -Icges(4,6) * t243 + t242 * t301;
t158 = -Icges(4,5) * t243 + t242 * t303;
t288 = t154 * t248 - t158 * t251;
t392 = t243 * t288;
t148 = -Icges(5,6) * t243 + t242 * t298;
t156 = -Icges(5,4) * t243 + t242 * t302;
t292 = t148 * t248 + t156 * t251;
t391 = t243 * t292;
t267 = Icges(6,4) * t87 + Icges(6,2) * t86 - Icges(6,6) * t354;
t268 = Icges(6,1) * t87 + Icges(6,4) * t86 - Icges(6,5) * t354;
t352 = t243 * t245;
t84 = t138 * t243 + t283 * t354;
t85 = t243 * t401 - t282 * t354;
t38 = Icges(6,4) * t85 + Icges(6,2) * t84 - Icges(6,6) * t352;
t39 = Icges(6,1) * t85 + Icges(6,4) * t84 - Icges(6,5) * t352;
t94 = Icges(6,5) * t138 - Icges(6,6) * t401;
t95 = Icges(6,4) * t138 - Icges(6,2) * t401;
t96 = Icges(6,1) * t138 - Icges(6,4) * t401;
t387 = -(-t127 * t401 + t129 * t138 - t141 * t352 + t142 * t84 + t143 * t85 - t170 * t95 + t171 * t96 - t242 * t94 - t282 * t38 - t283 * t39) * t242 / 0.2e1 + (-t126 * t401 + t138 * t128 - t141 * t354 + t142 * t86 + t143 * t87 - t168 * t95 + t169 * t96 + t243 * t94 - t267 * t282 - t268 * t283) * t243 / 0.2e1;
t386 = 2 * m(3);
t385 = 2 * m(4);
t384 = 2 * m(5);
t383 = 2 * m(6);
t382 = m(5) / 0.2e1;
t381 = m(6) / 0.2e1;
t380 = -pkin(3) - pkin(4);
t377 = -rSges(5,1) - pkin(3);
t376 = -rSges(6,3) - pkin(8);
t369 = rSges(4,2) * t248;
t370 = rSges(4,1) * t251;
t201 = (-t369 + t370) * qJD(3);
t375 = m(4) * t201;
t223 = rSges(4,1) * t248 + rSges(4,2) * t251;
t374 = m(4) * t223;
t249 = sin(qJ(1));
t373 = pkin(1) * t249;
t372 = pkin(8) * t243;
t237 = t242 * pkin(7);
t371 = t85 * rSges(6,1) + t84 * rSges(6,2);
t368 = pkin(1) * qJD(1);
t233 = t242 * rSges(5,2);
t232 = t242 * rSges(4,3);
t363 = -rSges(5,3) - qJ(4);
t358 = qJ(4) * t248;
t125 = Icges(6,5) * t171 - Icges(6,6) * t170 - Icges(6,3) * t242;
t357 = t125 * t242;
t356 = t125 * t243;
t144 = -rSges(6,1) * t283 - rSges(6,2) * t282;
t355 = t144 * t245;
t353 = t242 * t251;
t351 = t243 * t248;
t350 = t243 * t251;
t307 = -rSges(6,1) * t169 + rSges(6,2) * t168;
t130 = rSges(6,3) * t243 - t307;
t348 = pkin(4) * t353 + t130 + t372;
t345 = t171 * rSges(6,1) - t170 * rSges(6,2);
t131 = -rSges(6,3) * t242 + t345;
t230 = pkin(4) * t350;
t347 = -pkin(8) * t242 + t131 + t230;
t306 = pkin(3) * t251 + t358;
t174 = t306 * t242;
t231 = pkin(3) * t350;
t175 = qJ(4) * t351 + t231;
t346 = t242 * t174 + t243 * t175;
t179 = qJD(3) * t306 - qJD(4) * t251;
t308 = rSges(5,1) * t251 + rSges(5,3) * t248;
t344 = -t308 * qJD(3) - t179;
t328 = t242 * t349;
t343 = rSges(4,2) * t328 + rSges(4,3) * t352;
t320 = t243 * t332;
t331 = qJD(4) * t248;
t342 = qJ(4) * t320 + t243 * t331;
t341 = rSges(5,2) * t352 + rSges(5,3) * t320;
t221 = pkin(3) * t248 - qJ(4) * t251;
t222 = rSges(5,1) * t248 - rSges(5,3) * t251;
t339 = -t221 - t222;
t338 = t243 * rSges(4,3) + t242 * t369;
t337 = t243 * pkin(2) + t237;
t336 = t242 ^ 2 + t243 ^ 2;
t335 = qJD(3) * t242;
t334 = qJD(3) * t243;
t330 = t249 * t368;
t252 = cos(qJ(1));
t329 = t252 * t368;
t208 = pkin(3) * t323;
t321 = t243 * t333;
t272 = -t245 * t353 - t321;
t326 = t242 * (qJ(4) * t397 + t245 * t231 + t242 * t331 - t208) + t243 * (pkin(3) * t272 - qJ(4) * t328 + t342) + t174 * t352;
t325 = -rSges(4,1) * t323 - rSges(4,2) * t397;
t220 = pkin(7) * t352;
t324 = t220 + t342;
t162 = rSges(5,1) * t350 + rSges(5,3) * t351 + t233;
t317 = -t352 / 0.2e1;
t315 = -pkin(2) - t370;
t124 = Icges(6,5) * t169 - Icges(6,6) * t168 + Icges(6,3) * t243;
t266 = Icges(6,5) * t87 + Icges(6,6) * t86 - Icges(6,3) * t354;
t270 = -t124 * t242 - t126 * t170 + t128 * t171;
t29 = -t127 * t170 + t129 * t171 - t357;
t37 = Icges(6,5) * t85 + Icges(6,6) * t84 - Icges(6,3) * t352;
t1 = t29 * t352 - (-t124 * t352 + t84 * t126 + t85 * t128 - t170 * t267 + t171 * t268 - t242 * t266) * t243 + (t84 * t127 + t85 * t129 - t170 * t38 + t171 * t39 - t242 * t37 + (t270 - t356) * t245) * t242;
t271 = t124 * t243 - t126 * t168 + t128 * t169;
t28 = -t127 * t168 + t129 * t169 + t356;
t19 = t242 * t28 - t243 * t271;
t314 = t245 * t19 + t1;
t2 = t28 * t352 - (-t124 * t354 + t86 * t126 + t87 * t128 - t168 * t267 + t169 * t268 + t243 * t266) * t243 + (t86 * t127 + t87 * t129 - t168 * t38 + t169 * t39 + t243 * t37 + (t271 - t357) * t245) * t242;
t20 = t242 * t29 - t243 * t270;
t313 = t245 * t20 - t2;
t312 = pkin(4) * t248 + t144;
t189 = t243 * rSges(3,1) - rSges(3,2) * t242;
t311 = t144 * t336;
t140 = t339 * t243;
t310 = t337 + t175;
t309 = -t221 - t312;
t173 = -rSges(3,1) * t352 + rSges(3,2) * t354;
t188 = -rSges(3,1) * t242 - rSges(3,2) * t243;
t238 = t243 * pkin(7);
t273 = t251 * t380 - pkin(2) - t358;
t258 = t242 * t273 + t243 * t376;
t66 = t238 + t258 + t307;
t64 = t66 - t373;
t244 = t252 * pkin(1);
t67 = t242 * t376 + t230 + t310 + t345;
t65 = t244 + t67;
t305 = t242 * t65 + t243 * t64;
t304 = t242 * t67 + t243 * t66;
t293 = -t148 * t251 + t156 * t248;
t291 = t149 * t251 - t157 * t248;
t289 = t154 * t251 + t158 * t248;
t287 = t155 * t251 + t159 * t248;
t97 = rSges(6,1) * t138 - rSges(6,2) * t401;
t281 = -pkin(4) * t332 - t179 - t97;
t163 = rSges(4,1) * t350 - rSges(4,2) * t351 + t232;
t109 = t309 * t243;
t172 = t188 * t245;
t280 = t387 - t403 * t354 / 0.2e1 + t402 * t317;
t276 = t300 * t243;
t275 = t299 * t243;
t123 = t310 + t162;
t135 = t163 + t337;
t134 = t242 * t315 + t238 + t338;
t269 = t248 * t363 + t251 * t377 - pkin(2);
t265 = t269 * t242;
t262 = Icges(5,2) * t245 - qJD(3) * t211;
t259 = Icges(4,3) * t245 - qJD(3) * t210;
t235 = t243 * rSges(5,2);
t122 = t235 + t238 + t265;
t89 = (t315 * t243 + (-rSges(4,3) - pkin(7)) * t242) * t245 - t325;
t254 = t411 * qJD(3) + t138 * t143 - t401 * t142 - t282 * t95 - t283 * t96 + t409 * t332 + t410 * t333;
t88 = rSges(4,1) * t272 - rSges(4,2) * t320 - pkin(2) * t354 + t220 + t343;
t253 = t291 * t317 - t387 + (-t396 * t242 + (-t286 + t290) * qJD(3) - t411 * t354) * t242 / 0.2e1 - (t396 * t243 + (-t288 + t292) * qJD(3) + ((-t274 + t277) * t251 + (t278 + t279) * t248) * t245) * t243 / 0.2e1 + (t242 * t400 + t243 * t398 + t287 + t402) * t352 / 0.2e1 + (t398 * t242 - t243 * t400 + t289 + t293 + t403) * t354 / 0.2e1;
t60 = t245 * t265 + t321 * t377 + t324 + t341;
t24 = t245 * t258 + t321 * t380 + t324 + t371;
t203 = rSges(5,1) * t323;
t61 = t203 + t208 + (t363 * t332 - t331) * t242 + ((-rSges(5,2) - pkin(7)) * t242 + t269 * t243) * t245;
t25 = t208 + (-qJ(4) * t332 - t331) * t242 + (t243 * t273 - t237) * t245 - t395;
t178 = t189 + t244;
t177 = t188 - t373;
t176 = t221 * t354;
t161 = rSges(4,1) * t353 - t338;
t160 = t242 * t308 - t235;
t153 = Icges(5,2) * t242 + t276;
t152 = -Icges(5,2) * t243 + t242 * t300;
t151 = Icges(4,3) * t242 + t275;
t150 = -Icges(4,3) * t243 + t242 * t299;
t147 = t173 - t329;
t146 = t172 - t330;
t139 = t339 * t242;
t133 = t135 + t244;
t132 = t134 - t373;
t115 = t242 * t262 + t245 * t276;
t114 = t243 * t262 - t300 * t354;
t113 = t242 * t259 + t245 * t275;
t112 = t243 * t259 - t299 * t354;
t108 = t309 * t242;
t107 = t244 + t123;
t106 = t122 - t373;
t83 = t140 * t245 + t242 * t344;
t82 = t222 * t354 + t243 * t344 + t176;
t77 = t89 - t329;
t76 = t88 - t330;
t75 = t151 * t242 - t286 * t243;
t74 = t150 * t242 - t392;
t73 = t153 * t242 + t290 * t243;
t72 = t152 * t242 + t391;
t71 = -t151 * t243 - t394;
t70 = -t150 * t243 - t288 * t242;
t69 = -t153 * t243 + t393;
t68 = -t152 * t243 + t292 * t242;
t63 = t160 * t242 + t162 * t243 + t346;
t62 = -t130 * t242 - t131 * t243;
t59 = t61 - t329;
t58 = t60 - t330;
t40 = -rSges(6,3) * t352 + t371;
t36 = t109 * t245 + t242 * t281;
t35 = t243 * t281 + t312 * t354 + t176;
t30 = t242 * t348 + t243 * t347 + t346;
t23 = t25 - t329;
t22 = t24 - t330;
t21 = t242 * (rSges(5,3) * t322 - t203) + t243 * (-rSges(5,1) * t321 + t341) + (t243 * t160 + (-t162 - t175 + t233) * t242) * t245 + t326;
t18 = (-t245 * t130 - t40) * t243 + (t131 * t245 - t41) * t242;
t5 = (-pkin(4) * t321 + t40 + (t348 - t372) * t245) * t243 + ((-t175 - t347) * t245 + t395) * t242 + t326;
t3 = [(t22 * t65 + t23 * t64) * t383 + (t106 * t59 + t107 * t58) * t384 + (t132 * t77 + t133 * t76) * t385 + (t146 * t178 + t147 * t177) * t386 + t254; t254 + m(6) * (t22 * t67 + t23 * t66 + t24 * t65 + t25 * t64) + m(5) * (t106 * t61 + t107 * t60 + t122 * t59 + t123 * t58) + m(4) * (t132 * t89 + t133 * t88 + t134 * t77 + t135 * t76) + m(3) * (t146 * t189 + t147 * t188 + t172 * t178 + t173 * t177); t254 + (t24 * t67 + t25 * t66) * t383 + (t122 * t61 + t123 * t60) * t384 + (t134 * t89 + t135 * t88) * t385 + (t172 * t189 + t173 * t188) * t386; (-t132 * t243 - t133 * t242) * t375 + t253 + ((-t133 * t245 - t77) * t243 + (t132 * t245 - t76) * t242) * t374 + m(5) * (t106 * t82 + t107 * t83 + t139 * t58 + t140 * t59) + m(6) * (t108 * t22 + t109 * t23 + t35 * t64 + t36 * t65); m(5) * (t122 * t82 + t123 * t83 + t139 * t60 + t140 * t61) + m(6) * (t108 * t24 + t109 * t25 + t35 * t66 + t36 * t67) + t253 + ((-t135 * t245 - t89) * t243 + (t134 * t245 - t88) * t242) * t374 + (-t134 * t243 - t135 * t242) * t375; (t108 * t36 + t109 * t35 + t30 * t5) * t383 + t242 * t1 - t243 * t2 + (t139 * t83 + t140 * t82 + t21 * t63) * t384 + t242 * ((t242 * t112 + (t74 + t394) * t245) * t242 + (t75 * t245 + (t154 * t332 + t158 * t333) * t243 + (-t287 * qJD(3) - t245 * t288 - t113) * t242) * t243) + ((t161 * t242 + t163 * t243) * (((-t163 + t232) * t245 + t325) * t242 + (t245 * t161 - t223 * t334 + t343) * t243) + t336 * t223 * t201) * t385 - t243 * ((t243 * t115 + (t69 - t391) * t245) * t243 + (t68 * t245 + (t149 * t332 - t157 * t333) * t242 + (t293 * qJD(3) + t245 * t290 - t114) * t243) * t242) - t243 * ((t243 * t113 + (t71 + t392) * t245) * t243 + (t70 * t245 + (-t155 * t332 - t159 * t333) * t242 + (t289 * qJD(3) - t245 * t286 - t112) * t243) * t242) + t242 * ((t242 * t114 + (t72 - t393) * t245) * t242 + (t73 * t245 + (-t148 * t332 + t156 * t333) * t243 + (t291 * qJD(3) + t245 * t292 - t115) * t242) * t243) + (t19 + (-t68 - t70) * t243 + (t69 + t71) * t242) * t354 + (t20 + (-t72 - t74) * t243 + (t73 + t75) * t242) * t352; 0.2e1 * (t305 * t381 + (t106 * t243 + t107 * t242) * t382) * t332 + 0.2e1 * ((t22 * t242 + t23 * t243 + t352 * t65 - t354 * t64) * t381 + (-t106 * t354 + t107 * t352 + t242 * t58 + t243 * t59) * t382) * t248; 0.2e1 * (t304 * t381 + (t122 * t243 + t123 * t242) * t382) * t332 + 0.2e1 * ((t24 * t242 + t243 * t25 + t352 * t67 - t354 * t66) * t381 + (-t122 * t354 + t123 * t352 + t242 * t60 + t243 * t61) * t382) * t248; 0.2e1 * ((t108 * t335 + t109 * t334 - t5) * t381 + (t139 * t335 + t140 * t334 - t21) * t382) * t251 + 0.2e1 * ((qJD(3) * t30 + t108 * t352 - t109 * t354 + t242 * t36 + t243 * t35) * t381 + (qJD(3) * t63 + t139 * t352 - t140 * t354 + t242 * t83 + t243 * t82) * t382) * t248; 0.4e1 * (t382 + t381) * (-0.1e1 + t336) * t248 * t332; m(6) * (t305 * t97 + ((t245 * t65 + t23) * t243 + (-t245 * t64 + t22) * t242) * t144) + t280; m(6) * (t304 * t97 + ((t245 * t67 + t25) * t243 + (-t245 * t66 + t24) * t242) * t144) + t280; m(6) * (t18 * t30 + t5 * t62) + (m(6) * (t108 * t355 + t109 * t97 + t144 * t35) - t313) * t243 + (m(6) * (t108 * t97 - t109 * t355 + t144 * t36) - t314) * t242; m(6) * (-t18 * t251 + t336 * t97 * t248 + (t248 * t62 + t251 * t311) * qJD(3)); (t62 * t18 + t311 * t97) * t383 + t313 * t243 + t314 * t242;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
