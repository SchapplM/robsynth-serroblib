% Calculate time derivative of joint inertia matrix for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:19
% EndTime: 2019-12-31 19:07:33
% DurationCPUTime: 7.37s
% Computational Cost: add. (23908->619), mult. (23514->889), div. (0->0), fcn. (22139->10), ass. (0->329)
t253 = pkin(9) + qJ(3);
t246 = qJ(4) + t253;
t239 = sin(t246);
t254 = qJD(3) + qJD(4);
t262 = cos(qJ(5));
t260 = sin(qJ(5));
t401 = Icges(6,4) * t262;
t303 = -Icges(6,2) * t260 + t401;
t240 = cos(t246);
t379 = t240 * t254;
t402 = Icges(6,4) * t260;
t110 = t303 * t379 + (Icges(6,6) * t254 + (-Icges(6,2) * t262 - t402) * qJD(5)) * t239;
t169 = -Icges(6,6) * t240 + t239 * t303;
t307 = Icges(6,1) * t262 - t402;
t170 = -Icges(6,5) * t240 + t239 * t307;
t435 = -t260 * t110 + (-t169 * t262 - t170 * t260) * qJD(5);
t261 = sin(qJ(1));
t263 = cos(qJ(1));
t244 = sin(t253);
t245 = cos(t253);
t405 = Icges(4,4) * t245;
t306 = -Icges(4,2) * t244 + t405;
t192 = Icges(4,6) * t261 + t263 * t306;
t406 = Icges(4,4) * t244;
t310 = Icges(4,1) * t245 - t406;
t194 = Icges(4,5) * t261 + t263 * t310;
t290 = t192 * t244 - t194 * t245;
t434 = t261 * t290;
t403 = Icges(5,4) * t240;
t304 = -Icges(5,2) * t239 + t403;
t181 = Icges(5,6) * t261 + t263 * t304;
t404 = Icges(5,4) * t239;
t308 = Icges(5,1) * t240 - t404;
t183 = Icges(5,5) * t261 + t263 * t308;
t292 = t181 * t239 - t183 * t240;
t433 = t261 * t292;
t191 = -Icges(4,6) * t263 + t261 * t306;
t193 = -Icges(4,5) * t263 + t261 * t310;
t291 = t191 * t244 - t193 * t245;
t432 = t263 * t291;
t180 = -Icges(5,6) * t263 + t261 * t304;
t182 = -Icges(5,5) * t263 + t261 * t308;
t293 = t180 * t239 - t182 * t240;
t431 = t263 * t293;
t300 = Icges(6,5) * t262 - Icges(6,6) * t260;
t109 = t300 * t379 + (Icges(6,3) * t254 + (-Icges(6,5) * t260 - Icges(6,6) * t262) * qJD(5)) * t239;
t392 = t169 * t260;
t430 = -t254 * t392 - t109;
t369 = t262 * t263;
t371 = t261 * t260;
t207 = -t240 * t371 - t369;
t370 = t261 * t262;
t373 = t260 * t263;
t208 = t240 * t370 - t373;
t317 = -t208 * rSges(6,1) - t207 * rSges(6,2);
t381 = t239 * t261;
t145 = rSges(6,3) * t381 - t317;
t209 = -t240 * t373 + t370;
t210 = t240 * t369 + t371;
t380 = t239 * t263;
t146 = t210 * rSges(6,1) + t209 * rSges(6,2) + rSges(6,3) * t380;
t429 = -t261 * t145 - t263 * t146;
t259 = -pkin(6) - qJ(2);
t258 = cos(pkin(9));
t241 = t258 * pkin(2) + pkin(1);
t410 = rSges(4,2) * t244;
t412 = rSges(4,1) * t245;
t320 = -t410 + t412;
t283 = -t241 - t320;
t162 = (rSges(4,3) - t259) * t263 + t283 * t261;
t251 = t261 * rSges(4,3);
t360 = t263 * t412 + t251;
t163 = -t261 * t259 + (t241 - t410) * t263 + t360;
t428 = t162 * t263 + t163 * t261;
t301 = Icges(5,5) * t240 - Icges(5,6) * t239;
t178 = -Icges(5,3) * t263 + t261 * t301;
t427 = qJD(1) * t178;
t302 = Icges(4,5) * t245 - Icges(4,6) * t244;
t189 = -Icges(4,3) * t263 + t261 * t302;
t323 = qJD(1) * t240 - qJD(5);
t375 = t254 * t263;
t342 = t239 * t375;
t426 = t261 * t323 + t342;
t216 = Icges(5,2) * t240 + t404;
t217 = Icges(5,1) * t239 + t403;
t289 = t216 * t239 - t217 * t240;
t424 = qJD(1) * t289 + t301 * t254;
t285 = rSges(3,1) * t258 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t407 = rSges(3,3) + qJ(2);
t177 = t261 * t407 + t263 * t285;
t252 = -pkin(7) + t259;
t359 = t252 - t259;
t227 = pkin(3) * t245 + t241;
t363 = t227 - t241;
t164 = t261 * t363 + t263 * t359;
t423 = 2 * m(4);
t422 = 2 * m(5);
t421 = 2 * m(6);
t255 = t261 ^ 2;
t256 = t263 ^ 2;
t420 = t261 / 0.2e1;
t419 = -t263 / 0.2e1;
t418 = -rSges(6,3) - pkin(8);
t224 = rSges(4,1) * t244 + rSges(4,2) * t245;
t417 = m(4) * t224;
t218 = rSges(5,1) * t239 + rSges(5,2) * t240;
t416 = m(5) * t218;
t415 = pkin(3) * t244;
t414 = pkin(4) * t239;
t413 = pkin(4) * t240;
t411 = rSges(5,1) * t240;
t139 = Icges(6,5) * t208 + Icges(6,6) * t207 + Icges(6,3) * t381;
t141 = Icges(6,4) * t208 + Icges(6,2) * t207 + Icges(6,6) * t381;
t143 = Icges(6,1) * t208 + Icges(6,4) * t207 + Icges(6,5) * t381;
t299 = -t141 * t260 + t143 * t262;
t324 = -qJD(5) * t240 + qJD(1);
t287 = t324 * t262;
t377 = t254 * t261;
t343 = t239 * t377;
t130 = t261 * t287 + (-t263 * t323 + t343) * t260;
t288 = t324 * t260;
t376 = t254 * t262;
t131 = t323 * t369 + (-t239 * t376 + t288) * t261;
t354 = qJD(1) * t263;
t271 = t239 * t354 + t240 * t377;
t71 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t271;
t73 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t271;
t75 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t271;
t19 = (t254 * t299 - t71) * t240 + (t139 * t254 - t260 * t73 + t262 * t75 + (-t141 * t262 - t143 * t260) * qJD(5)) * t239;
t409 = t19 * t263;
t140 = Icges(6,5) * t210 + Icges(6,6) * t209 + Icges(6,3) * t380;
t142 = Icges(6,4) * t210 + Icges(6,2) * t209 + Icges(6,6) * t380;
t144 = Icges(6,1) * t210 + Icges(6,4) * t209 + Icges(6,5) * t380;
t298 = -t142 * t260 + t144 * t262;
t128 = t260 * t426 + t263 * t287;
t129 = -t262 * t426 + t263 * t288;
t355 = qJD(1) * t261;
t335 = t239 * t355;
t341 = t240 * t375;
t270 = -t335 + t341;
t70 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t270;
t72 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t270;
t74 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t270;
t20 = (t254 * t298 - t70) * t240 + (t140 * t254 - t260 * t72 + t262 * t74 + (-t142 * t262 - t144 * t260) * qJD(5)) * t239;
t408 = t20 * t261;
t250 = t261 * rSges(5,3);
t330 = -t218 - t415;
t173 = t330 * t263;
t389 = t173 * t263;
t388 = t191 * t245;
t387 = t192 * t245;
t386 = t193 * t244;
t385 = t194 * t244;
t384 = t216 * t254;
t383 = t217 * t254;
t382 = t239 * t254;
t378 = t240 * t263;
t319 = -rSges(5,2) * t239 + t411;
t200 = t319 * t254;
t372 = t261 * t200;
t316 = rSges(6,1) * t262 - rSges(6,2) * t260;
t112 = t316 * t379 + (rSges(6,3) * t254 + (-rSges(6,1) * t260 - rSges(6,2) * t262) * qJD(5)) * t239;
t321 = pkin(8) * t239 + t413;
t368 = -t321 * t254 - t112;
t203 = pkin(4) * t378 + pkin(8) * t380;
t367 = -t146 - t203;
t220 = t263 * t227;
t165 = -t241 * t263 - t261 * t359 + t220;
t366 = t261 * t164 + t263 * t165;
t171 = -rSges(6,3) * t240 + t239 * t316;
t166 = t171 * t355;
t219 = -pkin(8) * t240 + t414;
t365 = t219 * t355 + t166;
t364 = -t171 - t219;
t184 = -rSges(5,3) * t263 + t261 * t319;
t185 = rSges(5,1) * t378 - rSges(5,2) * t380 + t250;
t119 = t261 * t184 + t263 * t185;
t362 = rSges(5,2) * t335 + rSges(5,3) * t354;
t353 = qJD(3) * t244;
t348 = pkin(3) * t353;
t361 = -t252 * t355 - t261 * t348;
t358 = t255 + t256;
t179 = Icges(5,3) * t261 + t263 * t301;
t357 = qJD(1) * t179;
t190 = Icges(4,3) * t261 + t263 * t302;
t356 = qJD(1) * t190;
t352 = qJD(3) * t245;
t351 = qJD(3) * t261;
t349 = t263 * t410;
t347 = pkin(3) * t352;
t62 = -t140 * t240 + t239 * t298;
t168 = -Icges(6,3) * t240 + t239 * t300;
t82 = t168 * t380 + t209 * t169 + t210 * t170;
t346 = t62 / 0.2e1 + t82 / 0.2e1;
t61 = -t139 * t240 + t239 * t299;
t81 = t168 * t381 + t169 * t207 + t170 * t208;
t345 = t81 / 0.2e1 + t61 / 0.2e1;
t111 = t307 * t379 + (Icges(6,5) * t254 + (-Icges(6,1) * t260 - t401) * qJD(5)) * t239;
t340 = t239 * t262 * t111 + t240 * t170 * t376 + t168 * t382;
t272 = -t240 * t355 - t342;
t281 = t218 * t254;
t339 = t261 * (-t261 * t281 + (t263 * t319 + t250) * qJD(1)) + t263 * (rSges(5,1) * t272 - rSges(5,2) * t341 + t362) + t184 * t354;
t338 = t129 * rSges(6,1) + t128 * rSges(6,2) + rSges(6,3) * t341;
t238 = t259 * t355;
t337 = t261 * (t354 * t363 + t238 + t361) + t263 * (-qJD(1) * t164 - t263 * t348) + t164 * t354;
t248 = qJD(2) * t263;
t336 = t248 - t361;
t334 = t244 * t355;
t332 = t355 / 0.2e1;
t331 = t354 / 0.2e1;
t138 = t364 * t263;
t122 = -qJD(1) * t180 - t263 * t384;
t329 = t183 * t254 + t122;
t123 = qJD(1) * t181 - t261 * t384;
t328 = t182 * t254 + t123;
t124 = -qJD(1) * t182 - t263 * t383;
t327 = -t181 * t254 + t124;
t125 = qJD(1) * t183 - t261 * t383;
t326 = t180 * t254 - t125;
t325 = -t261 * t252 + t220;
t202 = t321 * t261;
t65 = t261 * t202 + t263 * t203 - t429;
t322 = t364 - t415;
t318 = t131 * rSges(6,1) + t130 * rSges(6,2);
t54 = t139 * t381 + t141 * t207 + t143 * t208;
t55 = t140 * t381 + t142 * t207 + t144 * t208;
t39 = t55 * t261 - t263 * t54;
t315 = t261 * t54 + t263 * t55;
t56 = t139 * t380 + t209 * t141 + t210 * t143;
t57 = t140 * t380 + t209 * t142 + t210 * t144;
t40 = t57 * t261 - t263 * t56;
t314 = t261 * t56 + t263 * t57;
t313 = t62 * t261 - t263 * t61;
t312 = t261 * t61 + t263 * t62;
t215 = Icges(5,5) * t239 + Icges(5,6) * t240;
t277 = t254 * t215;
t120 = -t263 * t277 - t427;
t121 = -t261 * t277 + t357;
t13 = t128 * t141 + t129 * t143 + t139 * t270 + t209 * t73 + t210 * t75 + t380 * t71;
t14 = t128 * t142 + t129 * t144 + t140 * t270 + t209 * t72 + t210 * t74 + t380 * t70;
t8 = qJD(1) * t314 - t13 * t263 + t14 * t261;
t86 = -t178 * t263 - t261 * t293;
t87 = -t179 * t263 - t433;
t88 = t261 * t178 - t431;
t89 = t261 * t179 - t263 * t292;
t311 = t40 * t354 + t39 * t355 + (-t88 * t354 - t86 * t355) * t263 + (t8 + (t89 * qJD(1) + (t123 * t239 - t125 * t240 + t180 * t379 + t182 * t382 - t427) * t263) * t263 + t87 * t355 + t89 * t354 + ((t88 + t433) * qJD(1) + (-t121 + t327 * t240 - t329 * t239 + (t179 - t293) * qJD(1)) * t263 + t261 * t120) * t261) * t261;
t309 = Icges(4,1) * t244 + t405;
t305 = Icges(4,2) * t245 + t406;
t297 = t145 * t263 - t146 * t261;
t282 = -t227 - t319;
t153 = (rSges(5,3) - t252) * t263 + t282 * t261;
t154 = t185 + t325;
t294 = t153 * t263 + t154 * t261;
t286 = -t347 + t368;
t226 = pkin(8) * t341;
t76 = -rSges(6,3) * t335 + t338;
t77 = rSges(6,3) * t271 + t318;
t284 = (t145 + t202) * t354 + (pkin(4) * t272 - pkin(8) * t335 + t226 + t76) * t263 + (t77 + t271 * pkin(8) + (t240 * t354 - t343) * pkin(4)) * t261;
t114 = t322 * t263;
t280 = qJD(3) * t224;
t276 = qJD(3) * t309;
t275 = qJD(3) * t305;
t274 = qJD(3) * (-Icges(4,5) * t244 - Icges(4,6) * t245);
t273 = t239 * t418 - t227 - t413;
t12 = (t263 * t121 + (t87 + t431) * qJD(1)) * t263 + (t86 * qJD(1) + (-t122 * t239 + t124 * t240 - t181 * t379 - t183 * t382 + t357) * t261 + (-t120 + t326 * t240 + t328 * t239 + (-t178 - t292) * qJD(1)) * t263) * t261;
t15 = t130 * t141 + t131 * t143 + t139 * t271 + t207 * t73 + t208 * t75 + t381 * t71;
t16 = t130 * t142 + t131 * t144 + t140 * t271 + t207 * t72 + t208 * t74 + t381 * t70;
t9 = qJD(1) * t315 - t15 * t263 + t16 * t261;
t269 = (-t12 - t9) * t263 + t311;
t24 = t239 * t315 - t81 * t240;
t25 = t239 * t314 - t82 * t240;
t30 = t109 * t380 + t209 * t110 + t210 * t111 + t128 * t169 + t129 * t170 + t168 * t270;
t3 = (t254 * t314 - t30) * t240 + (-qJD(1) * t40 + t13 * t261 + t14 * t263 + t254 * t82) * t239;
t31 = t109 * t381 + t207 * t110 + t208 * t111 + t130 * t169 + t131 * t170 + t168 * t271;
t4 = (t254 * t315 - t31) * t240 + (-qJD(1) * t39 + t15 * t261 + t16 * t263 + t254 * t81) * t239;
t268 = t3 * t420 + t9 * t381 / 0.2e1 + t4 * t419 - t240 * (qJD(1) * t312 + t408 - t409) / 0.2e1 + t24 * t332 - t40 * t335 / 0.2e1 + t313 * t382 / 0.2e1 + t8 * t380 / 0.2e1 + (t261 * t39 + t263 * t40) * t379 / 0.2e1 + (t239 * t39 + t25) * t331;
t267 = rSges(4,2) * t334 + rSges(4,3) * t354 - t263 * t280;
t266 = -t252 * t263 + t261 * t273;
t176 = -t261 * t285 + t263 * t407;
t196 = t304 * t254;
t197 = t308 * t254;
t265 = qJD(1) * t215 + (t197 - t384) * t240 + (-t196 - t383) * t239;
t264 = -t409 / 0.2e1 + t408 / 0.2e1 + (t239 * t327 + t240 * t329 + t261 * t424 + t265 * t263 + t30) * t420 + (-t239 * t326 + t240 * t328 + t265 * t261 - t263 * t424 + t31) * t419 + (t180 * t240 + t182 * t239 - t215 * t263 - t261 * t289 + t61 + t81) * t332 + (t181 * t240 + t183 * t239 + t261 * t215 - t263 * t289 + t62 + t82) * t331;
t247 = qJD(2) * t261;
t231 = pkin(3) * t334;
t214 = t320 * qJD(3);
t199 = -t349 + t360;
t198 = -rSges(4,3) * t263 + t261 * t320;
t172 = t330 * t261;
t161 = -t177 * qJD(1) + t248;
t160 = qJD(1) * t176 + t247;
t148 = t261 * t274 + t356;
t147 = -qJD(1) * t189 + t263 * t274;
t137 = t364 * t261;
t113 = t322 * t261;
t108 = -t218 * t354 - t372 + (-t244 * t354 - t245 * t351) * pkin(3);
t107 = t218 * t355 + t231 + (-t200 - t347) * t263;
t106 = t238 + t248 + t224 * t351 + (t263 * t283 - t251) * qJD(1);
t105 = t247 + (-t259 * t263 + (-t241 - t412) * t261) * qJD(1) + t267;
t99 = t261 * t190 - t290 * t263;
t98 = t261 * t189 - t432;
t97 = -t190 * t263 - t434;
t96 = -t189 * t263 - t261 * t291;
t95 = t325 - t367;
t94 = t266 + t317;
t93 = -t240 * t146 - t171 * t380;
t92 = t145 * t240 + t171 * t381;
t91 = t218 * t377 + (t263 * t282 - t250) * qJD(1) + t336;
t90 = t247 + (-t227 - t411) * t355 + (-qJD(1) * t252 - t281 - t348) * t263 + t362;
t85 = -t168 * t240 + (t170 * t262 - t392) * t239;
t84 = t297 * t239;
t83 = t85 * t382;
t80 = t119 + t366;
t67 = qJD(1) * t138 + t261 * t368;
t66 = t263 * t368 + t365;
t64 = qJD(1) * t114 + t261 * t286;
t63 = t263 * t286 + t231 + t365;
t60 = -t185 * t355 + t339;
t47 = (t240 * t418 + t414) * t377 + t273 * t354 - t318 + t336;
t46 = t226 + t247 + (-pkin(4) * t382 - t348) * t263 + t266 * qJD(1) + t338;
t45 = t65 + t366;
t44 = (t171 * t377 + t77) * t240 + (t261 * t112 - t145 * t254 + t171 * t354) * t239;
t43 = (-t171 * t375 - t76) * t240 + (-t112 * t263 + t146 * t254 + t166) * t239;
t41 = (-t165 - t185) * t355 + t337 + t339;
t36 = t239 * t435 + t430 * t240 + t340;
t27 = t297 * t379 + (qJD(1) * t429 - t261 * t76 + t263 * t77) * t239;
t26 = t355 * t367 + t284;
t21 = (-t165 + t367) * t355 + t284 + t337;
t1 = [(t46 * t95 + t47 * t94) * t421 + (t153 * t91 + t154 * t90) * t422 + t217 * t379 - t216 * t382 + (t105 * t163 + t106 * t162) * t423 + 0.2e1 * m(3) * (t160 * t177 + t161 * t176) + t340 + (t310 - t305) * t353 + (t309 + t306) * t352 + (t196 + t430) * t240 + (t197 + t435) * t239; m(6) * (t261 * t47 - t263 * t46 + (t261 * t95 + t263 * t94) * qJD(1)) + m(5) * (qJD(1) * t294 + t261 * t91 - t263 * t90) + m(4) * (qJD(1) * t428 - t105 * t263 + t261 * t106) + m(3) * (-t160 * t263 + t261 * t161 + (t176 * t263 + t177 * t261) * qJD(1)); 0; t264 + (-qJD(3) * t290 + (-qJD(1) * t191 - t263 * t275) * t245 + (-qJD(1) * t193 - t263 * t276) * t244) * t420 + (-qJD(3) * t291 + (qJD(1) * t192 - t261 * t275) * t245 + (qJD(1) * t194 - t261 * t276) * t244) * t419 + m(4) * ((-t105 * t261 - t106 * t263) * t224 - t428 * t214) + m(6) * (t113 * t46 + t114 * t47 + t63 * t94 + t64 * t95) + m(5) * (t107 * t153 + t108 * t154 + t172 * t90 + t173 * t91) + (t255 / 0.2e1 + t256 / 0.2e1) * t302 * qJD(3) + ((-t163 * t417 + t387 / 0.2e1 + t385 / 0.2e1) * t263 + (t162 * t417 + t388 / 0.2e1 + t386 / 0.2e1) * t261) * qJD(1); m(5) * (t107 * t261 - t108 * t263 + (t172 * t261 + t389) * qJD(1)) + m(6) * (t63 * t261 - t263 * t64 + (t113 * t261 + t114 * t263) * qJD(1)); (t113 * t64 + t114 * t63 + t45 * t21) * t421 - t263 * t9 + (t107 * t173 + t108 * t172 + t41 * t80) * t422 - t263 * t12 + (t99 * t261 - t263 * t98) * t354 + t261 * ((t261 * t147 + (t98 + t434) * qJD(1)) * t261 + (t99 * qJD(1) + (t191 * t352 + t193 * t353) * t263 + (-t148 + (-t385 - t387) * qJD(3) + (t190 - t291) * qJD(1)) * t261) * t263) + (t97 * t261 - t263 * t96) * t355 - t263 * ((t263 * t148 + (t97 + t432) * qJD(1)) * t263 + (t96 * qJD(1) + (-t192 * t352 - t194 * t353 + t356) * t261 + (-t147 + (t386 + t388) * qJD(3) - t290 * qJD(1)) * t263) * t261) + ((t261 * t198 + t199 * t263) * ((qJD(1) * t198 + t267) * t263 + (-t261 * t280 + (-t199 - t349 + t251) * qJD(1)) * t261) + t358 * t224 * t214) * t423 + t311; m(6) * (t137 * t46 + t138 * t47 + t66 * t94 + t67 * t95) + t264 + (-t261 * t90 - t263 * t91 + (t153 * t261 - t154 * t263) * qJD(1)) * t416 - m(5) * t294 * t200; m(6) * (t66 * t261 - t263 * t67 + (t137 * t261 + t138 * t263) * qJD(1)); m(6) * (t113 * t67 + t114 * t66 + t137 * t64 + t138 * t63 + t65 * t21 + t26 * t45) + m(5) * (t119 * t41 - t172 * t372 - t200 * t389 + t60 * t80) + (-t107 * t263 - t108 * t261 + (-t172 * t263 + t173 * t261) * qJD(1)) * t416 + t269; (t200 * t218 * t358 + t119 * t60) * t422 + (t137 * t67 + t138 * t66 + t65 * t26) * t421 + t269; m(6) * (t43 * t95 + t44 * t94 + t46 * t93 + t47 * t92) + t83 + (-t36 + (t261 * t345 + t263 * t346) * t254) * t240 + ((t20 / 0.2e1 + t30 / 0.2e1) * t263 + (t31 / 0.2e1 + t19 / 0.2e1) * t261 + (-t261 * t346 + t263 * t345) * qJD(1)) * t239; m(6) * (t44 * t261 - t263 * t43 + (t261 * t93 + t263 * t92) * qJD(1)); t268 + m(6) * (t113 * t43 + t114 * t44 + t21 * t84 + t27 * t45 + t63 * t92 + t64 * t93); t268 + m(6) * (t137 * t43 + t138 * t44 + t26 * t84 + t27 * t65 + t66 * t92 + t67 * t93); (t27 * t84 + t43 * t93 + t44 * t92) * t421 + (t36 * t240 - t83 + (t261 * t24 - t240 * t312 + t263 * t25) * t254) * t240 + (t263 * t3 + t261 * t4 + t312 * t382 + (-t19 * t261 - t20 * t263 - t254 * t85) * t240 + (t263 * t24 + t240 * t313 - t261 * t25) * qJD(1)) * t239;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
