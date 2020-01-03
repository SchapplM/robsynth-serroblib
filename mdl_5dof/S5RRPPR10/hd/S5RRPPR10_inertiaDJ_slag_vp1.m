% Calculate time derivative of joint inertia matrix for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR10_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:03
% EndTime: 2019-12-31 19:43:30
% DurationCPUTime: 14.97s
% Computational Cost: add. (11356->792), mult. (31339->1086), div. (0->0), fcn. (32997->8), ass. (0->346)
t447 = Icges(5,2) + Icges(4,3);
t450 = t447 / 0.2e1;
t448 = Icges(5,4) + Icges(4,5);
t449 = -Icges(4,6) + Icges(5,6);
t446 = Icges(5,6) / 0.2e1 - Icges(4,6) / 0.2e1;
t268 = sin(qJ(1));
t419 = t268 / 0.2e1;
t271 = cos(qJ(1));
t417 = -t271 / 0.2e1;
t439 = -qJD(1) / 0.2e1;
t267 = sin(qJ(2));
t270 = cos(qJ(2));
t264 = sin(pkin(8));
t265 = cos(pkin(8));
t302 = Icges(5,5) * t265 + Icges(5,3) * t264;
t306 = Icges(4,4) * t265 - Icges(4,2) * t264;
t438 = ((t302 - t306) * t270 + t449 * t267) * qJD(2);
t309 = Icges(5,1) * t265 + Icges(5,5) * t264;
t310 = Icges(4,1) * t265 - Icges(4,4) * t264;
t437 = ((t309 + t310) * t270 + t448 * t267) * qJD(2);
t190 = -Icges(5,6) * t270 + t267 * t302;
t193 = -Icges(4,6) * t270 + t267 * t306;
t436 = t190 - t193;
t194 = -Icges(5,4) * t270 + t267 * t309;
t195 = -Icges(4,5) * t270 + t267 * t310;
t435 = t194 + t195;
t393 = t267 * t271;
t266 = sin(qJ(5));
t269 = cos(qJ(5));
t296 = t264 * t266 + t265 * t269;
t209 = t296 * t267;
t402 = Icges(3,4) * t270;
t308 = -Icges(3,2) * t267 + t402;
t202 = Icges(3,6) * t268 + t271 * t308;
t403 = Icges(3,4) * t267;
t312 = Icges(3,1) * t270 - t403;
t204 = Icges(3,5) * t268 + t271 * t312;
t298 = t202 * t267 - t204 * t270;
t288 = t298 * t268;
t201 = -Icges(3,6) * t271 + t268 * t308;
t203 = -Icges(3,5) * t271 + t268 * t312;
t299 = t201 * t267 - t203 * t270;
t289 = t299 * t271;
t392 = t268 * t270;
t221 = t264 * t392 + t265 * t271;
t222 = -t264 * t271 + t265 * t392;
t394 = t267 * t268;
t128 = Icges(5,5) * t222 + Icges(5,6) * t394 + Icges(5,3) * t221;
t134 = Icges(4,4) * t222 - Icges(4,2) * t221 + Icges(4,6) * t394;
t433 = t128 - t134;
t391 = t270 * t271;
t223 = t264 * t391 - t268 * t265;
t224 = t268 * t264 + t265 * t391;
t129 = Icges(5,5) * t224 + Icges(5,6) * t393 + Icges(5,3) * t223;
t135 = Icges(4,4) * t224 - Icges(4,2) * t223 + Icges(4,6) * t393;
t432 = t129 - t135;
t136 = Icges(5,1) * t222 + Icges(5,4) * t394 + Icges(5,5) * t221;
t138 = Icges(4,1) * t222 - Icges(4,4) * t221 + Icges(4,5) * t394;
t431 = t136 + t138;
t137 = Icges(5,1) * t224 + Icges(5,4) * t393 + Icges(5,5) * t223;
t139 = Icges(4,1) * t224 - Icges(4,4) * t223 + Icges(4,5) * t393;
t430 = t137 + t139;
t429 = t449 * t264 + t448 * t265;
t428 = -rSges(3,2) * t393 + t268 * rSges(3,3);
t304 = Icges(3,5) * t270 - Icges(3,6) * t267;
t199 = -Icges(3,3) * t271 + t268 * t304;
t427 = 2 * m(3);
t426 = 2 * m(4);
t425 = 2 * m(5);
t424 = 2 * m(6);
t262 = t268 ^ 2;
t263 = t271 ^ 2;
t423 = m(4) / 0.2e1;
t422 = m(5) / 0.2e1;
t421 = m(6) / 0.2e1;
t420 = -pkin(3) - pkin(4);
t418 = t270 / 0.2e1;
t415 = -rSges(5,1) - pkin(3);
t414 = -rSges(6,3) - pkin(7);
t413 = pkin(2) * t267;
t412 = pkin(2) * t270;
t259 = t268 * pkin(6);
t411 = qJD(1) / 0.2e1;
t410 = rSges(3,3) * t271;
t409 = rSges(5,3) * t221;
t372 = qJD(1) * t271;
t335 = t270 * t372;
t370 = qJD(2) * t268;
t340 = t267 * t370;
t373 = qJD(1) * t268;
t174 = -t265 * t373 + (t335 - t340) * t264;
t408 = t174 * rSges(5,3);
t407 = -rSges(5,2) - qJ(3);
t406 = -rSges(4,3) - qJ(3);
t157 = t221 * t269 - t222 * t266;
t158 = t221 * t266 + t222 * t269;
t320 = -rSges(6,1) * t158 - rSges(6,2) * t157;
t92 = -rSges(6,3) * t394 - t320;
t405 = pkin(4) * t222 - pkin(7) * t394 + t92;
t218 = t224 * pkin(4);
t364 = pkin(7) * t393;
t159 = t223 * t269 - t224 * t266;
t160 = t223 * t266 + t224 * t269;
t389 = t160 * rSges(6,1) + t159 * rSges(6,2);
t93 = -rSges(6,3) * t393 + t389;
t404 = t218 - t364 + t93;
t398 = t264 * t267;
t397 = t264 * t270;
t396 = t265 * t267;
t395 = t265 * t270;
t297 = t264 * t269 - t265 * t266;
t208 = t297 * t267;
t125 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t270;
t390 = pkin(4) * t396 + pkin(7) * t270 + t125;
t164 = t224 * pkin(3) + t223 * qJ(4);
t228 = pkin(2) * t391 + qJ(3) * t393;
t388 = -t164 - t228;
t387 = -t174 * qJ(4) - t221 * qJD(4);
t368 = qJD(2) * t271;
t338 = t267 * t368;
t173 = -qJD(1) * t222 - t265 * t338;
t343 = t267 * t373;
t386 = t173 * pkin(4) + pkin(7) * t343;
t319 = qJ(3) * t267 + t412;
t220 = qJD(2) * t319 - qJD(3) * t270;
t318 = pkin(3) * t265 + qJ(4) * t264;
t369 = qJD(2) * t270;
t385 = -qJD(4) * t398 - t318 * t369 - t220;
t322 = rSges(4,1) * t265 - rSges(4,2) * t264;
t384 = -(rSges(4,3) * t267 + t270 * t322) * qJD(2) - t220;
t225 = t318 * t267;
t239 = -qJ(3) * t270 + t413;
t229 = t239 * t373;
t383 = t225 * t373 + t229;
t198 = -rSges(4,3) * t270 + t267 * t322;
t382 = -t198 - t239;
t227 = t319 * t268;
t381 = t268 * t227 + t271 * t228;
t380 = -t225 - t239;
t336 = t270 * t368;
t367 = qJD(3) * t267;
t379 = qJ(3) * t336 + t271 * t367;
t378 = rSges(3,2) * t343 + rSges(3,3) * t372;
t377 = t271 * pkin(1) + t259;
t213 = t221 * qJ(4);
t260 = t271 * pkin(6);
t376 = t260 - t213;
t375 = t262 + t263;
t200 = Icges(3,3) * t268 + t271 * t304;
t374 = qJD(1) * t200;
t371 = qJD(2) * t267;
t366 = t422 + t421;
t365 = -qJ(3) - t414;
t172 = qJD(1) * t221 + t264 * t338;
t72 = -qJD(5) * t160 - t172 * t269 - t173 * t266;
t73 = qJD(5) * t159 - t172 * t266 + t173 * t269;
t363 = t73 * rSges(6,1) + t72 * rSges(6,2) + rSges(6,3) * t343;
t86 = Icges(6,5) * t158 + Icges(6,6) * t157 - Icges(6,3) * t394;
t88 = Icges(6,4) * t158 + Icges(6,2) * t157 - Icges(6,6) * t394;
t90 = Icges(6,1) * t158 + Icges(6,4) * t157 - Icges(6,5) * t394;
t31 = t208 * t88 + t209 * t90 + t270 * t86;
t122 = Icges(6,5) * t209 + Icges(6,6) * t208 + Icges(6,3) * t270;
t123 = Icges(6,4) * t209 + Icges(6,2) * t208 + Icges(6,6) * t270;
t124 = Icges(6,1) * t209 + Icges(6,4) * t208 + Icges(6,5) * t270;
t43 = -t122 * t394 + t123 * t157 + t124 * t158;
t361 = -t43 / 0.2e1 - t31 / 0.2e1;
t87 = Icges(6,5) * t160 + Icges(6,6) * t159 - Icges(6,3) * t393;
t89 = Icges(6,4) * t160 + Icges(6,2) * t159 - Icges(6,6) * t393;
t91 = Icges(6,1) * t160 + Icges(6,4) * t159 - Icges(6,5) * t393;
t32 = t208 * t89 + t209 * t91 + t270 * t87;
t44 = -t122 * t393 + t159 * t123 + t160 * t124;
t360 = t44 / 0.2e1 + t32 / 0.2e1;
t130 = Icges(4,5) * t222 - Icges(4,6) * t221 + Icges(4,3) * t394;
t359 = t130 * t394;
t358 = t130 * t393;
t131 = Icges(4,5) * t224 - Icges(4,6) * t223 + Icges(4,3) * t393;
t357 = t131 * t394;
t356 = t131 * t393;
t132 = Icges(5,4) * t222 + Icges(5,2) * t394 + Icges(5,6) * t221;
t355 = t132 * t394;
t354 = t132 * t393;
t133 = Icges(5,4) * t224 + Icges(5,2) * t393 + Icges(5,6) * t223;
t353 = t133 * t394;
t352 = t133 * t393;
t246 = pkin(2) * t340;
t337 = t268 * t369;
t281 = t267 * t372 + t337;
t282 = -t270 * t373 - t338;
t351 = t268 * (pkin(2) * t335 + qJ(3) * t281 + t268 * t367 - t246) + t271 * (pkin(2) * t282 - qJ(3) * t343 + t379) + t227 * t372;
t350 = t173 * rSges(4,1) + t172 * rSges(4,2) + rSges(4,3) * t336;
t349 = t173 * rSges(5,1) + rSges(5,2) * t336 - t172 * rSges(5,3);
t321 = rSges(5,1) * t265 + rSges(5,3) * t264;
t348 = -(rSges(5,2) * t267 + t270 * t321) * qJD(2) + t385;
t197 = -rSges(5,2) * t270 + t267 * t321;
t347 = -t197 + t380;
t142 = t224 * rSges(5,1) + rSges(5,2) * t393 + t223 * rSges(5,3);
t143 = t224 * rSges(4,1) - t223 * rSges(4,2) + rSges(4,3) * t393;
t257 = pkin(6) * t372;
t346 = t257 + t379;
t345 = t246 + t387;
t344 = -pkin(1) - t412;
t339 = t267 * t369;
t334 = -t193 / 0.2e1 + t190 / 0.2e1;
t333 = t195 / 0.2e1 + t194 / 0.2e1;
t332 = -t369 / 0.2e1;
t162 = t382 * t271;
t331 = pkin(2) * t338;
t144 = -qJD(5) * t209 + t297 * t369;
t145 = qJD(5) * t208 + t296 * t369;
t85 = rSges(6,1) * t145 + rSges(6,2) * t144 - rSges(6,3) * t371;
t330 = -(pkin(4) * t395 - pkin(7) * t267) * qJD(2) - t85 + t385;
t329 = t380 - t390;
t163 = pkin(3) * t222 + t213;
t328 = t268 * t163 + t271 * t164 + t381;
t327 = t377 + t228;
t119 = t347 * t271;
t175 = qJD(1) * t224 - t265 * t340;
t74 = -qJD(5) * t158 + t174 * t269 - t175 * t266;
t75 = qJD(5) * t157 + t174 * t266 + t175 * t269;
t326 = t75 * rSges(6,1) + t74 * rSges(6,2);
t325 = rSges(3,1) * t270 - rSges(3,2) * t267;
t240 = rSges(3,1) * t267 + rSges(3,2) * t270;
t324 = -t175 * rSges(4,1) + t174 * rSges(4,2);
t323 = -rSges(4,1) * t222 + rSges(4,2) * t221;
t27 = t157 * t88 + t158 * t90 - t394 * t86;
t28 = t157 * t89 + t158 * t91 - t394 * t87;
t17 = t28 * t268 - t27 * t271;
t317 = t268 * t27 + t271 * t28;
t29 = t159 * t88 + t160 * t90 - t393 * t86;
t30 = t159 * t89 + t160 * t91 - t393 * t87;
t18 = t30 * t268 - t271 * t29;
t316 = t268 * t29 + t271 * t30;
t315 = t32 * t268 - t31 * t271;
t314 = -t268 * t31 - t271 * t32;
t313 = t268 * t93 - t271 * t92;
t311 = Icges(3,1) * t267 + t402;
t307 = Icges(3,2) * t270 + t403;
t207 = rSges(3,1) * t391 + t428;
t81 = t329 * t271;
t295 = -pkin(1) - t325;
t291 = t173 * pkin(3) - t172 * qJ(4) + t223 * qJD(4);
t294 = t163 * t372 + t268 * (t175 * pkin(3) - t387) + t271 * t291 + t351;
t290 = qJD(2) * t240;
t287 = qJD(2) * t311;
t286 = qJD(2) * t307;
t285 = qJD(2) * (-Icges(3,5) * t267 - Icges(3,6) * t270);
t284 = t267 * t407 + t344;
t283 = t267 * t406 + t344;
t280 = -t336 + t343;
t279 = t164 + t327;
t278 = t267 * t365 + t344;
t82 = Icges(6,5) * t145 + Icges(6,6) * t144 - Icges(6,3) * t371;
t83 = Icges(6,4) * t145 + Icges(6,2) * t144 - Icges(6,6) * t371;
t84 = Icges(6,1) * t145 + Icges(6,4) * t144 - Icges(6,5) * t371;
t275 = -t122 * t371 + t144 * t123 + t145 * t124 + t208 * t83 + t209 * t84 + t270 * t82;
t274 = t284 * t268;
t273 = t283 * t268;
t272 = t291 + t346;
t233 = t325 * qJD(2);
t206 = t268 * t325 - t410;
t178 = t207 + t377;
t177 = t268 * t295 + t260 + t410;
t161 = t382 * t268;
t147 = t268 * t285 + t374;
t146 = -qJD(1) * t199 + t271 * t285;
t141 = rSges(4,3) * t394 - t323;
t140 = rSges(5,1) * t222 + rSges(5,2) * t394 + t409;
t121 = t240 * t370 + ((-rSges(3,3) - pkin(6)) * t268 + t295 * t271) * qJD(1);
t120 = rSges(3,1) * t282 - rSges(3,2) * t336 - pkin(1) * t373 + t257 + t378;
t118 = t347 * t268;
t117 = t327 + t143;
t116 = t260 + t273 + t323;
t115 = t268 * t200 - t271 * t298;
t114 = t268 * t199 - t289;
t113 = -t200 * t271 - t288;
t112 = -t199 * t271 - t268 * t299;
t111 = Icges(4,1) * t175 - Icges(4,4) * t174 + Icges(4,5) * t281;
t110 = Icges(4,1) * t173 + Icges(4,4) * t172 - Icges(4,5) * t280;
t109 = Icges(5,1) * t175 + Icges(5,4) * t281 + Icges(5,5) * t174;
t108 = Icges(5,1) * t173 - Icges(5,4) * t280 - Icges(5,5) * t172;
t107 = Icges(4,4) * t175 - Icges(4,2) * t174 + Icges(4,6) * t281;
t106 = Icges(4,4) * t173 + Icges(4,2) * t172 - Icges(4,6) * t280;
t101 = Icges(5,5) * t175 + Icges(5,6) * t281 + Icges(5,3) * t174;
t100 = Icges(5,5) * t173 - Icges(5,6) * t280 - Icges(5,3) * t172;
t97 = qJD(1) * t162 + t268 * t384;
t96 = t198 * t373 + t271 * t384 + t229;
t80 = t329 * t268;
t78 = t279 + t142;
t77 = t222 * t415 + t274 + t376 - t409;
t67 = t268 * t141 + t143 * t271 + t381;
t66 = t246 + (t369 * t406 - t367) * t268 + (t271 * t283 - t259) * qJD(1) + t324;
t65 = qJD(1) * t273 - t331 + t346 + t350;
t64 = qJD(1) * t119 + t268 * t348;
t63 = t197 * t373 + t271 * t348 + t383;
t62 = t125 * t393 + t270 * t93;
t61 = -t125 * t394 - t270 * t92;
t60 = -t223 * t135 + t224 * t139 + t356;
t59 = -t223 * t134 + t224 * t138 + t358;
t58 = t223 * t129 + t224 * t137 + t352;
t57 = t223 * t128 + t224 * t136 + t354;
t56 = -t135 * t221 + t139 * t222 + t357;
t55 = -t134 * t221 + t138 * t222 + t359;
t54 = t129 * t221 + t137 * t222 + t353;
t53 = t128 * t221 + t136 * t222 + t355;
t52 = t393 * t414 + t218 + t279 + t389;
t51 = t222 * t420 + t268 * t278 + t320 + t376;
t49 = t122 * t270 + t123 * t208 + t124 * t209;
t48 = t313 * t267;
t47 = t268 * t140 + t142 * t271 + t328;
t46 = -t408 + t415 * t175 + (t369 * t407 - t367) * t268 + (t271 * t284 - t259) * qJD(1) + t345;
t45 = qJD(1) * t274 + t272 - t331 + t349;
t42 = -rSges(6,3) * t281 + t326;
t41 = -rSges(6,3) * t336 + t363;
t40 = Icges(6,1) * t75 + Icges(6,4) * t74 - Icges(6,5) * t281;
t39 = Icges(6,1) * t73 + Icges(6,4) * t72 + Icges(6,5) * t280;
t38 = Icges(6,4) * t75 + Icges(6,2) * t74 - Icges(6,6) * t281;
t37 = Icges(6,4) * t73 + Icges(6,2) * t72 + Icges(6,6) * t280;
t36 = Icges(6,5) * t75 + Icges(6,6) * t74 - Icges(6,3) * t281;
t35 = Icges(6,5) * t73 + Icges(6,6) * t72 + Icges(6,3) * t280;
t34 = qJD(1) * t81 + t268 * t330;
t33 = t271 * t330 + t373 * t390 + t383;
t26 = t268 * t405 + t271 * t404 + t328;
t25 = t268 * (rSges(4,3) * t337 - t324) + t271 * t350 + (t271 * t141 + (-t143 - t228) * t268) * qJD(1) + t351;
t24 = t420 * t175 + (t365 * t369 - t367) * t268 + (t271 * t278 - t259) * qJD(1) - t326 + t345;
t23 = (t270 * t414 - t413) * t368 + (-pkin(1) - t319) * t373 + t272 + t363 + t386;
t22 = (-t125 * t370 - t42) * t270 + (qJD(2) * t92 - t125 * t372 - t268 * t85) * t267;
t21 = (t125 * t368 + t41) * t270 + (-qJD(2) * t93 - t125 * t373 + t271 * t85) * t267;
t20 = t275 * t270;
t19 = t268 * (t175 * rSges(5,1) + rSges(5,2) * t337 + t408) + t271 * t349 + (t271 * t140 + (-t142 + t388) * t268) * qJD(1) + t294;
t16 = -t122 * t281 + t74 * t123 + t75 * t124 + t157 * t83 + t158 * t84 - t394 * t82;
t15 = t122 * t280 + t72 * t123 + t73 * t124 + t159 * t83 + t160 * t84 - t393 * t82;
t14 = t313 * t369 + (t268 * t41 - t271 * t42 + (t268 * t92 + t271 * t93) * qJD(1)) * t267;
t13 = -t267 * t316 + t44 * t270;
t12 = -t267 * t317 + t43 * t270;
t11 = t144 * t89 + t145 * t91 + t208 * t37 + t209 * t39 + t270 * t35 - t371 * t87;
t10 = t144 * t88 + t145 * t90 + t208 * t38 + t209 * t40 + t270 * t36 - t371 * t86;
t9 = (-pkin(7) * t336 + t386 + t41) * t271 + (t175 * pkin(4) - pkin(7) * t337 + t42) * t268 + (t405 * t271 + (-t364 + t388 - t404) * t268) * qJD(1) + t294;
t8 = -t87 * t337 + t157 * t37 + t158 * t39 + t74 * t89 + t75 * t91 + (-t268 * t35 - t372 * t87) * t267;
t7 = -t86 * t337 + t157 * t38 + t158 * t40 + t74 * t88 + t75 * t90 + (-t268 * t36 - t372 * t86) * t267;
t6 = -t87 * t336 + t159 * t37 + t160 * t39 + t72 * t89 + t73 * t91 + (-t271 * t35 + t373 * t87) * t267;
t5 = -t86 * t336 + t159 * t38 + t160 * t40 + t72 * t88 + t73 * t90 + (-t271 * t36 + t373 * t86) * t267;
t4 = qJD(1) * t317 + t8 * t268 - t271 * t7;
t3 = qJD(1) * t316 + t6 * t268 - t271 * t5;
t2 = (-qJD(2) * t317 + t16) * t270 + (qJD(1) * t17 - qJD(2) * t43 - t268 * t7 - t271 * t8) * t267;
t1 = (-qJD(2) * t316 + t15) * t270 + (qJD(1) * t18 - qJD(2) * t44 - t268 * t5 - t271 * t6) * t267;
t50 = [t275 + (t120 * t178 + t121 * t177) * t427 + (t23 * t52 + t24 * t51) * t424 + (t45 * t78 + t46 * t77) * t425 + (t116 * t66 + t117 * t65) * t426 + t438 * t398 + t437 * t396 + (t267 * t429 - t270 * t447 - t307 + t312) * t371 + (t264 * t436 + t265 * t435 - t267 * t447 - t429 * t270 + t308 + t311) * t369; m(3) * ((-t120 * t268 - t121 * t271) * t240 + (-t177 * t271 - t178 * t268) * t233) + m(4) * (t116 * t96 + t117 * t97 + t161 * t65 + t162 * t66) + m(6) * (t23 * t80 + t24 * t81 + t33 * t51 + t34 * t52) + m(5) * (t118 * t45 + t119 * t46 + t63 * t77 + t64 * t78) + ((t202 * t439 + t286 * t419 + t448 * t175 / 0.2e1 + t281 * t450 + t446 * t174) * t271 + (t201 * t439 + t286 * t417 - t448 * t173 / 0.2e1 + t280 * t450 + t446 * t172) * t268) * t270 + ((t223 * t334 + t224 * t333 + t360) * t271 + (t221 * t334 + t222 * t333 - t361) * t268 + m(3) * (t177 * t268 - t178 * t271) * t240 + ((-t131 / 0.2e1 - t133 / 0.2e1 + t202 / 0.2e1) * t271 + (t201 / 0.2e1 - t130 / 0.2e1 - t132 / 0.2e1) * t268) * t270 + (t264 * t432 + t265 * t430 + t204) * t393 / 0.2e1) * qJD(1) + ((t262 / 0.2e1 + t263 / 0.2e1) * t304 - t288 / 0.2e1 + t289 / 0.2e1) * qJD(2) + (t11 + t15 + (t430 * t395 + t432 * t397) * qJD(2) + (-t271 * t287 + (t108 + t110) * t265 + (t100 - t106) * t264 + (t131 + t133) * qJD(2) + (t264 * t433 + t265 * t431) * qJD(1)) * t267 + t437 * t224 + t438 * t223 + t435 * t173 - t436 * t172) * t419 + (t10 + t16 + (qJD(1) * t204 - t268 * t287 + (t109 + t111) * t265 + (t101 - t107) * t264) * t267 + (t433 * t397 + t431 * t395 + (t130 + t132) * t267) * qJD(2) + t437 * t222 + t438 * t221 + t435 * t175 + t436 * t174) * t417; (t26 * t9 + t33 * t81 + t34 * t80) * t424 + t268 * t3 - t271 * t4 + (t118 * t64 + t119 * t63 + t19 * t47) * t425 + (t161 * t97 + t162 * t96 + t25 * t67) * t426 + t268 * ((-t223 * t106 + t224 * t110 + t172 * t135 + t173 * t139 + (t59 - t357) * qJD(1)) * t268 + (t223 * t107 - t224 * t111 - t172 * t134 - t173 * t138 + (t60 + t359) * qJD(1)) * t271) + t268 * ((t223 * t100 + t224 * t108 - t172 * t129 + t173 * t137 + (t57 - t353) * qJD(1)) * t268 + (-t223 * t101 - t224 * t109 + t172 * t128 - t173 * t136 + (t58 + t355) * qJD(1)) * t271) + t268 * ((t268 * t146 + (t114 + t288) * qJD(1)) * t268 + (t115 * qJD(1) + (t201 * t369 + t203 * t371) * t271 + (-t147 + (-t202 * t270 - t204 * t267) * qJD(2) + (t200 - t299) * qJD(1)) * t268) * t271) - t271 * ((t221 * t107 - t222 * t111 + t174 * t134 - t175 * t138 + (t56 - t358) * qJD(1)) * t271 + (-t221 * t106 + t222 * t110 - t174 * t135 + t175 * t139 + (t55 + t356) * qJD(1)) * t268) - t271 * ((-t221 * t101 - t222 * t109 - t174 * t128 - t175 * t136 + (t54 - t354) * qJD(1)) * t271 + (t221 * t100 + t222 * t108 + t174 * t129 + t175 * t137 + (t53 + t352) * qJD(1)) * t268) + ((t268 * t206 + t207 * t271) * ((qJD(1) * t206 - t271 * t290 + t378) * t271 + (-t268 * t290 + (-t207 + t428) * qJD(1)) * t268) + t375 * t240 * t233) * t427 - t271 * ((t271 * t147 + (t113 + t289) * qJD(1)) * t271 + (t112 * qJD(1) + (-t202 * t369 - t204 * t371 + t374) * t268 + (-t146 + (t201 * t270 + t203 * t267) * qJD(2) - t298 * qJD(1)) * t271) * t268) + (t17 + (-t112 - t53 - t55) * t271 + (t113 + t54 + t56) * t268) * t373 + (t18 + (-t114 - t57 - t59) * t271 + (t115 + t58 + t60) * t268) * t372; 0.2e1 * ((t268 * t52 + t271 * t51) * t421 + (t268 * t78 + t271 * t77) * t422 + (t116 * t271 + t117 * t268) * t423) * t369 + 0.2e1 * ((t23 * t268 + t24 * t271 + t372 * t52 - t373 * t51) * t421 + (t268 * t45 + t271 * t46 + t372 * t78 - t373 * t77) * t422 + (-t116 * t373 + t117 * t372 + t268 * t65 + t271 * t66) * t423) * t267; 0.2e1 * ((t368 * t81 + t370 * t80 - t9) * t421 + (t118 * t370 + t119 * t368 - t19) * t422 + (t161 * t370 + t162 * t368 - t25) * t423) * t270 + 0.2e1 * ((qJD(2) * t26 + t268 * t34 + t271 * t33 + t372 * t80 - t373 * t81) * t421 + (qJD(2) * t47 + t118 * t372 - t119 * t373 + t268 * t64 + t271 * t63) * t422 + (qJD(2) * t67 + t161 * t372 - t162 * t373 + t268 * t97 + t271 * t96) * t423) * t267; 0.4e1 * (t423 + t366) * (-0.1e1 + t375) * t339; m(6) * (-t172 * t51 + t174 * t52 + t221 * t23 + t223 * t24) + m(5) * (-t172 * t77 + t174 * t78 + t221 * t45 + t223 * t46); m(6) * (-t172 * t81 + t174 * t80 + t221 * t34 + t223 * t33 + (t26 * t369 + t267 * t9) * t264) + m(5) * (t118 * t174 - t119 * t172 + t221 * t64 + t223 * t63 + (t19 * t267 + t369 * t47) * t264); 0.2e1 * t366 * ((t268 * t221 + t271 * t223 - t397) * t369 + (t264 * t371 - t172 * t271 + t174 * t268 + (t271 * t221 - t223 * t268) * qJD(1)) * t267); 0.4e1 * t366 * (t264 ^ 2 * t339 - t172 * t223 + t174 * t221); m(6) * (t21 * t52 + t22 * t51 + t23 * t62 + t24 * t61) + t20 + (t268 * t361 - t271 * t360) * t369 + (-qJD(2) * t49 + (-t15 / 0.2e1 - t11 / 0.2e1) * t271 + (-t16 / 0.2e1 - t10 / 0.2e1) * t268 + (t268 * t360 + t271 * t361) * qJD(1)) * t267; m(6) * (t14 * t26 + t21 * t80 + t22 * t81 + t33 * t61 + t34 * t62 + t48 * t9) + (t13 * t411 - t2 / 0.2e1 + t18 * t332 + (t32 * qJD(1) - t10) * t418) * t271 + (t1 / 0.2e1 + t12 * t411 + t17 * t332 + (t31 * qJD(1) + t11) * t418) * t268 + (t3 * t417 - t268 * t4 / 0.2e1 - qJD(2) * t315 / 0.2e1 + (t17 * t417 + t18 * t419) * qJD(1)) * t267; m(6) * ((-t14 + (t268 * t62 + t271 * t61) * qJD(2)) * t270 + (qJD(2) * t48 + t21 * t268 + t22 * t271 + (-t268 * t61 + t271 * t62) * qJD(1)) * t267); m(6) * (-t172 * t61 + t174 * t62 + t21 * t221 + t22 * t223 + (t14 * t267 + t369 * t48) * t264); (t14 * t48 + t21 * t62 + t22 * t61) * t424 + (t20 + (-t268 * t12 - t271 * t13 + t270 * t314) * qJD(2)) * t270 + (-t271 * t1 - t268 * t2 + t270 * (-t10 * t268 - t11 * t271) + (-t267 * t314 - 0.2e1 * t49 * t270) * qJD(2) + (-t271 * t12 + t268 * t13 + t270 * t315) * qJD(1)) * t267;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t50(1), t50(2), t50(4), t50(7), t50(11); t50(2), t50(3), t50(5), t50(8), t50(12); t50(4), t50(5), t50(6), t50(9), t50(13); t50(7), t50(8), t50(9), t50(10), t50(14); t50(11), t50(12), t50(13), t50(14), t50(15);];
Mq = res;
