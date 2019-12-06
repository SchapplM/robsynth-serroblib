% Calculate time derivative of joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:03
% EndTime: 2019-12-05 19:00:18
% DurationCPUTime: 7.10s
% Computational Cost: add. (19037->566), mult. (12740->770), div. (0->0), fcn. (9780->10), ass. (0->338)
t264 = qJ(1) + qJ(2);
t255 = sin(t264);
t263 = qJ(3) + qJ(4);
t256 = cos(t263);
t425 = rSges(5,1) * t256;
t257 = cos(t264);
t254 = sin(t263);
t388 = t254 * t255;
t460 = rSges(5,2) * t388 + t257 * rSges(5,3);
t140 = -t255 * t425 + t460;
t261 = qJD(1) + qJD(2);
t267 = cos(qJ(3));
t364 = qJD(3) * t267;
t265 = sin(qJ(3));
t382 = t257 * t265;
t466 = t255 * t364 + t261 * t382;
t383 = t257 * t261;
t260 = qJD(3) + qJD(4);
t384 = t256 * t260;
t465 = t254 * t383 + t255 * t384;
t258 = qJ(5) + t263;
t248 = cos(t258);
t253 = qJD(5) + t260;
t389 = t248 * t253;
t247 = sin(t258);
t390 = t247 * t253;
t464 = rSges(6,1) * t390 + rSges(6,2) * t389;
t415 = Icges(4,4) * t267;
t315 = -Icges(4,2) * t265 + t415;
t290 = t315 * t255;
t155 = Icges(4,6) * t257 - t290;
t416 = Icges(4,4) * t265;
t318 = Icges(4,1) * t267 - t416;
t157 = Icges(4,5) * t257 - t255 * t318;
t304 = t155 * t265 - t157 * t267;
t463 = t257 * t304;
t413 = Icges(5,4) * t256;
t314 = -Icges(5,2) * t254 + t413;
t289 = t314 * t255;
t136 = Icges(5,6) * t257 - t289;
t414 = Icges(5,4) * t254;
t317 = Icges(5,1) * t256 - t414;
t138 = Icges(5,5) * t257 - t255 * t317;
t307 = t136 * t254 - t138 * t256;
t462 = t257 * t307;
t411 = Icges(6,4) * t248;
t313 = -Icges(6,2) * t247 + t411;
t288 = t313 * t255;
t127 = Icges(6,6) * t257 - t288;
t412 = Icges(6,4) * t247;
t316 = Icges(6,1) * t248 - t412;
t291 = t316 * t255;
t129 = Icges(6,5) * t257 - t291;
t309 = t127 * t247 - t129 * t248;
t461 = t257 * t309;
t421 = rSges(6,2) * t247;
t212 = t255 * t421;
t371 = t257 * rSges(6,3) + t212;
t200 = Icges(5,1) * t254 + t413;
t391 = t200 * t260;
t459 = -Icges(5,5) * t261 + t391;
t199 = Icges(5,2) * t256 + t414;
t392 = t199 * t260;
t458 = -Icges(5,6) * t261 + t392;
t180 = Icges(6,1) * t247 + t411;
t394 = t180 * t253;
t457 = -Icges(6,5) * t261 + t394;
t179 = Icges(6,2) * t248 + t412;
t395 = t179 * t253;
t456 = -Icges(6,6) * t261 + t395;
t137 = Icges(5,6) * t255 + t257 * t314;
t139 = Icges(5,5) * t255 + t257 * t317;
t306 = t137 * t254 - t139 * t256;
t283 = t306 * t255;
t81 = -t257 * t458 - t261 * t289;
t327 = t139 * t260 + t81;
t82 = t255 * t458 - t314 * t383;
t328 = t138 * t260 + t82;
t292 = t317 * t261;
t83 = -t255 * t292 - t257 * t459;
t329 = -t137 * t260 + t83;
t84 = t255 * t459 - t257 * t292;
t330 = -t136 * t260 + t84;
t311 = Icges(5,5) * t256 - Icges(5,6) * t254;
t134 = Icges(5,3) * t257 - t255 * t311;
t37 = t134 * t257 + t255 * t307;
t135 = Icges(5,3) * t255 + t257 * t311;
t38 = t135 * t257 + t283;
t387 = t254 * t260;
t39 = t134 * t255 - t462;
t40 = t135 * t255 - t257 * t306;
t286 = t311 * t261;
t198 = Icges(5,5) * t254 + Icges(5,6) * t256;
t453 = -Icges(5,3) * t261 + t198 * t260;
t79 = -t255 * t286 - t257 * t453;
t80 = t255 * t453 - t257 * t286;
t455 = (t255 * t40 + t257 * t39) * t383 + t255 * ((t255 * t79 + (-t39 + t283) * t261) * t255 + (t40 * t261 + (-t136 * t384 - t138 * t387 - t254 * t82 + t256 * t84) * t257 + (t80 + (-t138 * t261 + t329) * t256 + (t136 * t261 - t327) * t254) * t255) * t257) + t257 * ((t257 * t80 + (t38 + t462) * t261) * t257 + (-t37 * t261 + (t137 * t384 + t139 * t387 + t254 * t81 - t256 * t83) * t255 + (t79 + (-t139 * t261 - t330) * t256 + (t137 * t261 + t328) * t254) * t257) * t255);
t178 = Icges(6,5) * t247 + Icges(6,6) * t248;
t454 = -Icges(6,3) * t261 + t178 * t253;
t249 = t267 * pkin(3) + pkin(2);
t429 = pkin(2) - t249;
t452 = pkin(7) * t255 + t257 * t429;
t361 = rSges(5,1) * t387;
t422 = rSges(5,2) * t257;
t451 = t140 * t261 - t257 * t361 - t384 * t422;
t229 = Icges(4,5) * t265 + Icges(4,6) * t267;
t450 = -Icges(4,3) * t261 + qJD(3) * t229;
t230 = Icges(4,2) * t267 + t416;
t449 = -Icges(4,6) * t261 + qJD(3) * t230;
t231 = Icges(4,1) * t265 + t415;
t448 = -Icges(4,5) * t261 + qJD(3) * t231;
t205 = t315 * qJD(3);
t206 = t318 * qJD(3);
t447 = qJD(3) * (t230 * t267 + t231 * t265) + t205 * t265 - t206 * t267 - t229 * t261;
t148 = t313 * t253;
t149 = t316 * t253;
t446 = t247 * (t148 + t394) + t248 * (-t149 + t395) - t178 * t261;
t166 = t314 * t260;
t167 = t317 * t260;
t445 = t254 * (t166 + t391) + t256 * (-t167 + t392) - t198 * t261;
t444 = 2 * m(3);
t443 = 2 * m(4);
t442 = 2 * m(5);
t441 = 2 * m(6);
t269 = -pkin(8) - pkin(7);
t440 = t255 / 0.2e1;
t439 = t257 / 0.2e1;
t438 = -rSges(4,3) - pkin(7);
t423 = rSges(4,2) * t265;
t426 = rSges(4,1) * t267;
t209 = (-t423 + t426) * qJD(3);
t437 = m(4) * t209;
t236 = rSges(4,1) * t265 + rSges(4,2) * t267;
t436 = m(4) * t236;
t169 = (-rSges(5,2) * t254 + t425) * t260;
t435 = m(5) * t169;
t201 = rSges(5,1) * t254 + rSges(5,2) * t256;
t434 = m(5) * t201;
t268 = cos(qJ(1));
t433 = pkin(1) * t268;
t432 = pkin(3) * t265;
t266 = sin(qJ(1));
t430 = t266 * pkin(1);
t365 = qJD(3) * t265;
t181 = -pkin(3) * t365 - pkin(4) * t387;
t341 = t257 * t365;
t379 = t261 * t269;
t386 = t255 * t261;
t348 = pkin(3) * t341 + t249 * t386 + t257 * t379;
t216 = pkin(4) * t256 + t249;
t262 = -pkin(9) + t269;
t380 = t261 * t262;
t373 = t216 * t386 + t257 * t380;
t424 = rSges(6,1) * t248;
t359 = t255 * t424;
t349 = -t464 * t257 - t261 * t359;
t61 = t257 * (t261 * t371 + t349);
t428 = t257 * (t181 * t257 + t348 - t373) + t61;
t322 = (-t181 + t380) * t255;
t343 = t255 * t365;
t369 = pkin(3) * t343 + t255 * t379;
t370 = t216 - t249;
t419 = rSges(6,3) * t255;
t295 = -t257 * t424 - t419;
t213 = t257 * t421;
t350 = t261 * t213 + t464 * t255;
t72 = t261 * t295 + t350;
t427 = t370 * t383 - t322 + t369 - t72;
t420 = rSges(5,3) * t255;
t418 = pkin(1) * qJD(1);
t417 = t255 * rSges(4,3);
t245 = t257 * rSges(4,3);
t397 = t169 * t255;
t385 = t255 * t265;
t381 = t257 * t269;
t240 = t255 * t262;
t241 = t255 * t269;
t109 = t257 * t370 - t240 + t241;
t132 = -t213 - t295;
t117 = t257 * t132;
t378 = t257 * t109 + t117;
t131 = -t359 + t371;
t377 = -(-t262 + t269) * t257 + t370 * t255 - t131;
t376 = -t109 - t132;
t246 = t257 * pkin(7);
t123 = t255 * t429 - t246 - t381;
t375 = -t123 - t140;
t184 = rSges(6,1) * t247 + rSges(6,2) * t248;
t221 = pkin(4) * t388;
t374 = t184 * t386 + t261 * t221;
t372 = t466 * pkin(3);
t119 = t255 * t184 + t221;
t170 = rSges(3,1) * t386 + rSges(3,2) * t383;
t237 = rSges(4,2) * t385;
t367 = t237 + t245;
t366 = t255 ^ 2 + t257 ^ 2;
t128 = Icges(6,6) * t255 + t257 * t313;
t130 = Icges(6,5) * t255 + t257 * t316;
t308 = t128 * t247 - t130 * t248;
t284 = t308 * t255;
t310 = Icges(6,5) * t248 - Icges(6,6) * t247;
t285 = t310 * t255;
t125 = Icges(6,3) * t257 - t285;
t33 = t125 * t257 + t255 * t309;
t68 = -t257 * t456 - t261 * t288;
t331 = t130 * t253 + t68;
t69 = t255 * t456 - t313 * t383;
t332 = t129 * t253 + t69;
t70 = -t257 * t457 - t261 * t291;
t333 = -t128 * t253 + t70;
t71 = t255 * t457 - t316 * t383;
t334 = -t127 * t253 + t71;
t126 = Icges(6,3) * t255 + t257 * t310;
t34 = t126 * t257 + t284;
t35 = t125 * t255 - t461;
t36 = t126 * t255 - t257 * t308;
t66 = -t257 * t454 - t261 * t285;
t67 = t255 * t454 - t310 * t383;
t363 = t255 * ((t255 * t66 + (-t35 + t284) * t261) * t255 + (t36 * t261 + (-t127 * t389 - t129 * t390 - t247 * t69 + t248 * t71) * t257 + (t67 + (-t129 * t261 + t333) * t248 + (t127 * t261 - t331) * t247) * t255) * t257) + t257 * ((t257 * t67 + (t34 + t461) * t261) * t257 + (-t33 * t261 + (t128 * t389 + t130 * t390 + t247 * t68 - t248 * t70) * t255 + (t66 + (-t130 * t261 - t334) * t248 + (t128 * t261 + t332) * t247) * t257) * t255) + (t255 * t36 + t257 * t35) * t383;
t239 = pkin(3) * t385;
t362 = t255 * t426;
t357 = t268 * t418;
t356 = pkin(3) * t364;
t351 = -t123 + t377;
t347 = t465 * rSges(5,2) + t255 * t361;
t345 = -t257 * rSges(4,2) * t364 - rSges(4,1) * t341 - t261 * t362;
t344 = rSges(4,1) * t343 + t466 * rSges(4,2);
t340 = -t386 / 0.2e1;
t339 = t383 / 0.2e1;
t338 = -pkin(2) - t426;
t337 = -pkin(4) * t254 - t184;
t203 = -rSges(3,1) * t257 + t255 * rSges(3,2);
t336 = -t249 - t425;
t335 = -t216 - t424;
t152 = (-t421 + t424) * t253;
t58 = t465 * pkin(4) + t255 * t152 + t184 * t383;
t321 = -pkin(4) * t384 - t152;
t171 = -rSges(3,1) * t383 + rSges(3,2) * t386;
t202 = -rSges(3,1) * t255 - rSges(3,2) * t257;
t312 = Icges(4,5) * t267 - Icges(4,6) * t265;
t305 = t155 * t267 + t157 * t265;
t156 = Icges(4,6) * t255 + t257 * t315;
t158 = Icges(4,5) * t255 + t257 * t318;
t303 = t156 * t267 + t158 * t265;
t302 = t156 * t265 - t158 * t267;
t301 = t179 * t247 - t180 * t248;
t300 = t199 * t254 - t200 * t256;
t298 = t230 * t265 - t231 * t267;
t13 = t255 * t34 + t257 * t33;
t297 = -t13 * t386 + t363;
t296 = -t257 * t425 - t420;
t293 = t318 * t261;
t287 = t312 * t261;
t282 = t302 * t255;
t277 = t310 * t253 + t261 * t301;
t281 = (t247 * t333 + t248 * t331 + t277 * t255 - t257 * t446) * t440 + (t247 * t334 + t248 * t332 + t255 * t446 + t277 * t257) * t439 + (t127 * t248 + t129 * t247 + t178 * t257 + t255 * t301) * t340 + (t128 * t248 + t130 * t247 + t178 * t255 - t257 * t301) * t339;
t279 = t257 * t336 - t420;
t278 = t257 * t335 - t419;
t276 = t311 * t260 + t261 * t300;
t275 = t312 * qJD(3) + t261 * t298;
t114 = t255 * t338 + t246 + t367;
t274 = t255 * t438 + t257 * t338;
t17 = t255 * t38 + t257 * t37;
t273 = (-t13 - t17) * t386 + t363 + t455;
t238 = rSges(4,2) * t382;
t115 = t238 + t274;
t220 = t254 * t422;
t107 = t220 + t241 + t279;
t101 = t213 + t240 + t278;
t100 = t255 * t335 - t257 * t262 + t371;
t106 = t255 * t336 - t381 + t460;
t235 = pkin(2) * t386;
t73 = t235 + (t257 * t438 - t237) * t261 - t345;
t272 = t281 + (t254 * t329 + t276 * t255 + t256 * t327 - t257 * t445) * t440 + (t254 * t330 + t255 * t445 + t256 * t328 + t276 * t257) * t439 + (t136 * t256 + t138 * t254 + t198 * t257 + t255 * t300) * t340 + (t137 * t256 + t139 * t254 + t198 * t255 - t257 * t300) * t339;
t45 = t348 - t451;
t74 = t261 * t274 + t344;
t31 = -t261 * t212 + (-rSges(6,3) * t261 - t181) * t257 - t349 + t373;
t46 = t261 * t279 + t347 + t369;
t32 = t261 * t278 + t322 + t350;
t271 = t248 * t148 + t247 * t149 + t256 * t166 + t254 * t167 - t179 * t390 + t180 * t389 - t199 * t387 + t200 * t384 + t267 * t205 + t265 * t206 - t230 * t365 + t231 * t364;
t270 = t272 + (-qJD(3) * t302 + t275 * t255 - t257 * t447 + t265 * (-t255 * t293 - t257 * t448) + t267 * (-t257 * t449 - t261 * t290)) * t440 + (-qJD(3) * t304 + t255 * t447 + t275 * t257 + t265 * (t255 * t448 - t257 * t293) + t267 * (t255 * t449 - t315 * t383)) * t439 + (t229 * t257 + t255 * t298 + t305) * t340 + (t229 * t255 - t257 * t298 + t303) * t339;
t250 = t266 * t418;
t210 = t261 * t239;
t173 = t203 - t433;
t172 = t202 - t430;
t160 = t257 * t426 - t238 + t417;
t159 = -t362 + t367;
t154 = Icges(4,3) * t255 + t257 * t312;
t153 = Icges(4,3) * t257 - t255 * t312;
t146 = t171 - t357;
t145 = t250 + t170;
t143 = (-t201 - t432) * t257;
t142 = t201 * t255 + t239;
t141 = -t220 - t296;
t124 = -t241 - t452;
t120 = t337 * t257;
t118 = t257 * t141;
t116 = t257 * t124;
t113 = (t337 - t432) * t257;
t112 = t239 + t119;
t111 = t115 - t433;
t110 = t114 - t430;
t104 = t107 - t433;
t103 = t106 - t430;
t102 = t261 * t452 + t369;
t95 = t255 * t450 - t257 * t287;
t94 = -t255 * t287 - t257 * t450;
t93 = t257 * (-pkin(7) * t383 + t235 - t348);
t92 = t101 - t433;
t91 = t100 - t430;
t88 = t261 * t296 + t347;
t85 = -t140 * t255 + t118;
t78 = t201 * t383 + t372 + t397;
t77 = t201 * t386 + t210 + (-t169 - t356) * t257;
t76 = t257 * t451;
t75 = -t131 * t255 + t117;
t60 = t74 - t357;
t59 = t250 + t73;
t57 = t257 * t321 + t374;
t50 = t154 * t255 - t302 * t257;
t49 = t153 * t255 - t463;
t48 = t154 * t257 + t282;
t47 = t153 * t257 + t304 * t255;
t44 = t58 + t372;
t43 = t210 + (t321 - t356) * t257 + t374;
t42 = t46 - t357;
t41 = t250 + t45;
t30 = t32 - t357;
t29 = t250 + t31;
t28 = t255 * t375 + t116 + t118;
t25 = t255 * t377 + t378;
t22 = -t140 * t383 + t76 + (-t141 * t261 - t88) * t255;
t16 = -t131 * t383 + t61 + (-t132 * t261 - t72) * t255;
t14 = t255 * t351 + t116 + t378;
t7 = t76 + t93 + t375 * t383 + (-t102 - t88 + (-t124 - t141) * t261) * t255;
t6 = t377 * t383 + (t261 * t376 + t427) * t255 + t428;
t5 = t93 + t351 * t383 + (-t102 + (-t124 + t376) * t261 + t427) * t255 + t428;
t1 = [(t145 * t173 + t146 * t172) * t444 + (t110 * t60 + t111 * t59) * t443 + (t29 * t92 + t30 * t91) * t441 + (t103 * t42 + t104 * t41) * t442 + t271; m(6) * (t100 * t30 + t101 * t29 + t31 * t92 + t32 * t91) + m(5) * (t103 * t46 + t104 * t45 + t106 * t42 + t107 * t41) + m(3) * (t145 * t203 + t146 * t202 + t170 * t173 + t171 * t172) + m(4) * (t110 * t74 + t111 * t73 + t114 * t60 + t115 * t59) + t271; (t100 * t32 + t101 * t31) * t441 + (t106 * t46 + t107 * t45) * t442 + (t114 * t74 + t115 * t73) * t443 + (t170 * t203 + t171 * t202) * t444 + t271; ((t111 * t261 - t60) * t257 + (t110 * t261 + t59) * t255) * t436 + t270 + m(5) * (t103 * t77 + t104 * t78 + t142 * t41 + t143 * t42) + m(6) * (t112 * t29 + t113 * t30 + t43 * t91 + t44 * t92) + (-t110 * t257 + t111 * t255) * t437; m(5) * (t106 * t77 + t107 * t78 + t142 * t45 + t143 * t46) + m(6) * (t100 * t43 + t101 * t44 + t112 * t31 + t113 * t32) + t270 + (-t114 * t257 + t115 * t255) * t437 + ((t115 * t261 - t74) * t257 + (t114 * t261 + t73) * t255) * t436; (t112 * t44 + t113 * t43 + t14 * t5) * t441 + (t142 * t78 + t143 * t77 + t28 * t7) * t442 + (t255 * t50 + t257 * t49) * t383 + t255 * ((t255 * t94 + (-t49 + t282) * t261) * t255 + (t50 * t261 + (-t155 * t364 - t157 * t365) * t257 + (-t303 * qJD(3) + t304 * t261 + t95) * t255) * t257) + t257 * ((t257 * t95 + (t48 + t463) * t261) * t257 + (-t47 * t261 + (t156 * t364 + t158 * t365) * t255 + (t305 * qJD(3) + t302 * t261 + t94) * t257) * t255) + ((-t159 * t255 + t160 * t257) * (t257 * t345 - t255 * t344 + ((-t159 + t245) * t257 + (t417 - t160 + (t423 + t426) * t257) * t255) * t261) + t366 * t236 * t209) * t443 + t297 + (-t255 * t48 - t257 * t47 - t17) * t386 + t455; ((t104 * t261 - t42) * t257 + (t103 * t261 + t41) * t255) * t434 + m(6) * (t119 * t29 + t120 * t30 + t57 * t91 + t58 * t92) + (-t103 * t257 + t104 * t255) * t435 + t272; m(6) * (t100 * t57 + t101 * t58 + t119 * t31 + t120 * t32) + ((t107 * t261 - t46) * t257 + (t106 * t261 + t45) * t255) * t434 + (-t106 * t257 + t107 * t255) * t435 + t272; m(6) * (t112 * t58 + t113 * t57 + t119 * t44 + t120 * t43 + t14 * t6 + t25 * t5) + m(5) * (-t143 * t169 * t257 + t142 * t397 + t22 * t28 + t7 * t85) + ((t142 * t261 - t77) * t257 + (t143 * t261 + t78) * t255) * t434 + t273; (t169 * t201 * t366 + t22 * t85) * t442 + (t119 * t58 + t120 * t57 + t25 * t6) * t441 + t273; m(6) * ((t255 * t92 - t257 * t91) * t152 + ((t261 * t92 - t30) * t257 + (t261 * t91 + t29) * t255) * t184) + t281; m(6) * ((-t100 * t257 + t101 * t255) * t152 + ((t101 * t261 - t32) * t257 + (t100 * t261 + t31) * t255) * t184) + t281; m(6) * (t14 * t16 + t5 * t75 + (t112 * t255 - t113 * t257) * t152 + ((t112 * t261 - t43) * t257 + (t113 * t261 + t44) * t255) * t184) + t297; m(6) * (t16 * t25 + t6 * t75 + (t119 * t255 - t120 * t257) * t152 + ((t119 * t261 - t57) * t257 + (t120 * t261 + t58) * t255) * t184) + t297; (t152 * t184 * t366 + t16 * t75) * t441 + t297;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
