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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:14:47
% EndTime: 2020-01-03 12:15:02
% DurationCPUTime: 7.62s
% Computational Cost: add. (19037->590), mult. (12740->796), div. (0->0), fcn. (9780->10), ass. (0->345)
t266 = qJ(1) + qJ(2);
t255 = sin(t266);
t257 = cos(t266);
t265 = qJ(3) + qJ(4);
t258 = qJ(5) + t265;
t247 = sin(t258);
t248 = cos(t258);
t313 = Icges(6,5) * t248 - Icges(6,6) * t247;
t288 = t313 * t255;
t133 = -Icges(6,3) * t257 + t288;
t412 = Icges(6,4) * t248;
t316 = -Icges(6,2) * t247 + t412;
t291 = t316 * t255;
t135 = -Icges(6,6) * t257 + t291;
t413 = Icges(6,4) * t247;
t319 = Icges(6,1) * t248 - t413;
t294 = t319 * t255;
t137 = -Icges(6,5) * t257 + t294;
t312 = t135 * t247 - t137 * t248;
t287 = t312 * t257;
t35 = -t133 * t255 + t287;
t134 = -Icges(6,3) * t255 - t257 * t313;
t136 = -Icges(6,6) * t255 - t257 * t316;
t138 = -Icges(6,5) * t255 - t257 * t319;
t311 = t136 * t247 - t138 * t248;
t36 = -t134 * t255 + t257 * t311;
t13 = -t255 * t36 - t257 * t35;
t254 = sin(t265);
t256 = cos(t265);
t314 = Icges(5,5) * t256 - Icges(5,6) * t254;
t289 = t314 * t255;
t142 = -Icges(5,3) * t257 + t289;
t414 = Icges(5,4) * t256;
t317 = -Icges(5,2) * t254 + t414;
t292 = t317 * t255;
t144 = -Icges(5,6) * t257 + t292;
t415 = Icges(5,4) * t254;
t320 = Icges(5,1) * t256 - t415;
t295 = t320 * t255;
t146 = -Icges(5,5) * t257 + t295;
t310 = t144 * t254 - t146 * t256;
t286 = t310 * t257;
t39 = -t142 * t255 + t286;
t143 = -Icges(5,3) * t255 - t257 * t314;
t145 = -Icges(5,6) * t255 - t257 * t317;
t147 = -Icges(5,5) * t255 - t257 * t320;
t309 = t145 * t254 - t147 * t256;
t40 = -t143 * t255 + t257 * t309;
t17 = -t255 * t40 - t257 * t39;
t456 = -t13 - t17;
t201 = rSges(5,1) * t254 + rSges(5,2) * t256;
t267 = sin(qJ(3));
t428 = pkin(3) * t267;
t150 = (-t201 - t428) * t255;
t263 = qJD(1) + qJD(2);
t455 = t150 * t263;
t269 = cos(qJ(3));
t416 = Icges(4,4) * t269;
t318 = -Icges(4,2) * t267 + t416;
t164 = -Icges(4,6) * t255 - t257 * t318;
t417 = Icges(4,4) * t267;
t321 = Icges(4,1) * t269 - t417;
t166 = -Icges(4,5) * t255 - t257 * t321;
t305 = t164 * t267 - t166 * t269;
t454 = t255 * t305;
t453 = t255 * t309;
t452 = t255 * t311;
t451 = -t255 * rSges(3,1) - t257 * rSges(3,2);
t200 = Icges(5,1) * t254 + t414;
t262 = qJD(3) + qJD(4);
t389 = t200 * t262;
t450 = -Icges(5,5) * t263 + t389;
t199 = Icges(5,2) * t256 + t415;
t390 = t199 * t262;
t449 = -Icges(5,6) * t263 + t390;
t186 = Icges(6,1) * t247 + t412;
t253 = qJD(5) + t262;
t394 = t186 * t253;
t448 = -Icges(6,5) * t263 + t394;
t185 = Icges(6,2) * t248 + t413;
t395 = t185 * t253;
t447 = -Icges(6,6) * t263 + t395;
t184 = Icges(6,5) * t247 + Icges(6,6) * t248;
t446 = -Icges(6,3) * t263 + t184 * t253;
t198 = Icges(5,5) * t254 + Icges(5,6) * t256;
t445 = -Icges(5,3) * t263 + t198 * t262;
t219 = Icges(4,5) * t267 + Icges(4,6) * t269;
t444 = -Icges(4,3) * t263 + qJD(3) * t219;
t220 = Icges(4,2) * t269 + t417;
t443 = -Icges(4,6) * t263 + qJD(3) * t220;
t221 = Icges(4,1) * t267 + t416;
t442 = -Icges(4,5) * t263 + qJD(3) * t221;
t205 = t318 * qJD(3);
t206 = t321 * qJD(3);
t441 = t205 * t267 - t206 * t269 - t219 * t263 + (t220 * t269 + t221 * t267) * qJD(3);
t156 = t316 * t253;
t157 = t319 * t253;
t440 = t247 * (t156 + t394) + t248 * (-t157 + t395) - t184 * t263;
t174 = t317 * t262;
t175 = t320 * t262;
t439 = t254 * (t174 + t389) + t256 * (-t175 + t390) - t198 * t263;
t438 = 2 * m(3);
t437 = 2 * m(4);
t436 = 2 * m(5);
t435 = 2 * m(6);
t271 = -pkin(8) - pkin(7);
t434 = -t255 / 0.2e1;
t433 = -t257 / 0.2e1;
t422 = rSges(4,2) * t267;
t426 = rSges(4,1) * t269;
t209 = (-t422 + t426) * qJD(3);
t432 = m(4) * t209;
t232 = rSges(4,1) * t267 + rSges(4,2) * t269;
t431 = m(4) * t232;
t421 = rSges(5,2) * t254;
t425 = rSges(5,1) * t256;
t177 = (-t421 + t425) * t262;
t430 = m(5) * t177;
t429 = m(5) * t201;
t249 = t269 * pkin(3) + pkin(2);
t360 = qJD(3) * t267;
t355 = pkin(3) * t360;
t383 = t254 * t262;
t187 = -pkin(4) * t383 - t355;
t171 = t257 * t187;
t212 = pkin(4) * t256 + t249;
t345 = t257 * t360;
t326 = pkin(3) * t345;
t264 = -pkin(9) + t271;
t361 = t264 - t271;
t420 = rSges(6,2) * t247;
t356 = t255 * t420;
t378 = t257 * t263;
t371 = -rSges(6,3) * t378 - t263 * t356;
t382 = t255 * t263;
t385 = t253 * t257;
t419 = rSges(6,2) * t248;
t72 = t385 * t419 + (t247 * t385 + t248 * t382) * rSges(6,1) + t371;
t427 = t326 + t171 - (t361 * t257 + (t212 - t249) * t255) * t263 - t72;
t424 = rSges(6,1) * t247;
t423 = rSges(6,1) * t248;
t418 = pkin(1) * qJD(1);
t402 = t133 * t263;
t401 = t134 * t263;
t400 = t142 * t263;
t399 = t143 * t263;
t160 = (-t420 + t423) * t253;
t398 = t160 * t255;
t397 = t177 * t255;
t190 = t419 + t424;
t393 = t190 * t255;
t392 = t190 * t263;
t388 = t247 * t253;
t387 = t248 * t253;
t386 = t253 * t255;
t384 = t254 * t257;
t381 = t255 * t269;
t380 = t256 * t262;
t379 = t257 * t262;
t377 = t257 * t267;
t237 = t257 * t271;
t365 = t255 * t249 + t237;
t372 = t255 * t212 + t257 * t264;
t110 = -t365 + t372;
t139 = -rSges(6,3) * t257 + t255 * t423 - t356;
t125 = t255 * t139;
t376 = t255 * t110 + t125;
t183 = t257 * t212;
t214 = t257 * t249;
t111 = t255 * t361 - t183 + t214;
t211 = t257 * t423;
t140 = -t255 * rSges(6,3) + t257 * t420 - t211;
t375 = -t111 - t140;
t347 = t256 * t379;
t374 = pkin(4) * t347 + t257 * t160;
t373 = t255 * t187 + t212 * t378;
t370 = rSges(6,3) * t382 + t263 * t211;
t357 = t255 * t421;
t369 = rSges(5,3) * t378 + t263 * t357;
t216 = t257 * t425;
t368 = rSges(5,3) * t382 + t263 * t216;
t358 = t255 * t422;
t367 = rSges(4,3) * t378 + t263 * t358;
t234 = t257 * t426;
t366 = rSges(4,3) * t382 + t263 * t234;
t128 = pkin(4) * t384 + t257 * t190;
t364 = -pkin(2) * t378 - pkin(7) * t382;
t363 = t257 * pkin(2) + t255 * pkin(7);
t362 = t255 ^ 2 + t257 ^ 2;
t359 = qJD(3) * t269;
t268 = sin(qJ(1));
t354 = t268 * t418;
t245 = t255 * pkin(2);
t131 = pkin(7) * t257 - t245 + t365;
t327 = t255 * t271 - t214;
t132 = t327 + t363;
t193 = t249 * t378;
t298 = -t263 * t271 - t355;
t353 = t132 * t382 + t131 * t378 + t255 * (t255 * t298 + t193 + t364);
t349 = t247 * t378;
t352 = t140 * t382 + t139 * t378 + t255 * (-t386 * t424 + (-t248 * t386 - t349) * rSges(6,2) + t370);
t148 = -rSges(5,3) * t257 + t255 * t425 - t357;
t149 = rSges(5,2) * t384 - t255 * rSges(5,3) - t216;
t348 = t254 * t378;
t282 = -t255 * t380 - t348;
t351 = t149 * t382 + t148 * t378 + t255 * (-rSges(5,1) * t255 * t383 + rSges(5,2) * t282 + t368);
t350 = t190 * t382;
t68 = t257 * t447 + t263 * t291;
t336 = t138 * t253 + t68;
t70 = t257 * t448 + t263 * t294;
t338 = -t136 * t253 + t70;
t66 = t257 * t446 + t263 * t288;
t67 = -t255 * t446 + t313 * t378;
t69 = -t255 * t447 + t316 * t378;
t71 = -t255 * t448 + t319 * t378;
t1 = (t255 * t66 + (t35 + t452) * t263) * t255 + (-t36 * t263 + (-t135 * t387 - t137 * t388 - t247 * t69 + t248 * t71 + t402) * t257 + (t401 + t67 + (-t137 * t263 + t338) * t248 + (t135 * t263 - t336) * t247) * t255) * t257;
t33 = -t133 * t257 - t255 * t312;
t34 = -t134 * t257 - t452;
t8 = (-t255 * t34 - t257 * t33) * t382;
t346 = -t255 * t1 + t8;
t344 = t257 * t359;
t343 = t382 / 0.2e1;
t342 = -t378 / 0.2e1;
t337 = t137 * t253 + t69;
t339 = -t135 * t253 + t71;
t2 = (t257 * t67 + (-t34 + t287) * t263) * t257 + (t33 * t263 + (t136 * t387 + t138 * t388 + t247 * t68 - t248 * t70 - t401) * t255 + (-t402 + t66 + (-t138 * t263 - t339) * t248 + (t136 * t263 + t337) * t247) * t257) * t255;
t341 = -t263 * t13 - t2;
t340 = -pkin(4) * t254 - t190;
t203 = t257 * rSges(3,1) - rSges(3,2) * t255;
t84 = -t255 * t450 + t320 * t378;
t335 = -t144 * t262 + t84;
t83 = t257 * t450 + t263 * t295;
t334 = -t145 * t262 + t83;
t82 = -t255 * t449 + t317 * t378;
t333 = t146 * t262 + t82;
t81 = t257 * t449 + t263 * t292;
t332 = t147 * t262 + t81;
t179 = rSges(3,1) * t378 - rSges(3,2) * t382;
t325 = rSges(4,1) * t381 - t358;
t127 = t340 * t255;
t37 = -t142 * t257 - t255 * t310;
t38 = -t143 * t257 - t453;
t15 = (-t255 * t38 - t257 * t37) * t382;
t79 = t257 * t445 + t263 * t289;
t80 = -t255 * t445 + t314 * t378;
t3 = (t255 * t79 + (t39 + t453) * t263) * t255 + (-t40 * t263 + (-t144 * t380 - t146 * t383 - t254 * t82 + t256 * t84 + t400) * t257 + (t399 + t80 + (-t146 * t263 + t334) * t256 + (t144 * t263 - t332) * t254) * t255) * t257;
t322 = t15 + t8 + (-t1 - t3) * t255;
t315 = Icges(4,5) * t269 - Icges(4,6) * t267;
t293 = t318 * t255;
t163 = -Icges(4,6) * t257 + t293;
t296 = t321 * t255;
t165 = -Icges(4,5) * t257 + t296;
t308 = t163 * t269 + t165 * t267;
t307 = t163 * t267 - t165 * t269;
t306 = t164 * t269 + t166 * t267;
t304 = t185 * t247 - t186 * t248;
t303 = t199 * t254 - t200 * t256;
t301 = t220 * t267 - t221 * t269;
t168 = rSges(4,2) * t377 - t255 * rSges(4,3) - t234;
t300 = t340 - t428;
t299 = t111 * t382 + t110 * t378 + t255 * (-t193 + (-t263 * t361 + t355) * t255 + t373) + t352;
t178 = t451 * t263;
t297 = qJD(3) * t232;
t290 = t315 * t255;
t285 = t307 * t257;
t279 = -t313 * t253 - t263 * t304;
t284 = (t247 * t338 + t248 * t336 + t279 * t255 + t257 * t440) * t434 + (t247 * t339 + t248 * t337 - t255 * t440 + t279 * t257) * t433 + (t135 * t248 + t137 * t247 - t184 * t257 - t255 * t304) * t343 + (t136 * t248 + t138 * t247 - t184 * t255 + t257 * t304) * t342;
t283 = t300 * t263;
t281 = -t255 * t359 - t263 * t377;
t123 = -t168 + t363;
t278 = -t314 * t262 - t263 * t303;
t277 = -t315 * qJD(3) - t263 * t301;
t108 = t148 + t365;
t100 = t139 + t372;
t109 = -t149 - t327;
t101 = -t255 * t264 - t140 + t183;
t122 = t245 + (-rSges(4,3) - pkin(7)) * t257 + t325;
t276 = -t190 * t253 - t263 * t264;
t275 = -t201 * t262 + t298;
t274 = t284 + (t254 * t334 + t278 * t255 + t256 * t332 + t257 * t439) * t434 + (t254 * t335 - t255 * t439 + t256 * t333 + t278 * t257) * t433 + (t144 * t256 + t146 * t254 - t198 * t257 - t255 * t303) * t343 + (t145 * t256 + t147 * t254 - t198 * t255 + t257 * t303) * t342;
t74 = -rSges(4,1) * t255 * t360 + rSges(4,2) * t281 - t364 + t366;
t273 = t248 * t156 + t247 * t157 + t256 * t174 + t254 * t175 - t185 * t388 + t186 * t387 - t199 * t383 + t200 * t380 + t269 * t205 + t267 * t206 - t220 * t360 + t221 * t359;
t230 = pkin(7) * t378;
t73 = -rSges(4,2) * t344 - pkin(2) * t382 + t230 + (-t263 * t381 - t345) * rSges(4,1) + t367;
t32 = -rSges(6,2) * t349 + t255 * t276 + t370 + t373;
t272 = t274 + (-qJD(3) * t305 + t277 * t255 + t257 * t441 + t267 * (t257 * t442 + t263 * t296) + t269 * (t257 * t443 + t263 * t293)) * t434 + (-qJD(3) * t307 - t255 * t441 + t277 * t257 + t267 * (-t255 * t442 + t321 * t378) + t269 * (-t255 * t443 + t318 * t378)) * t433 + (-t219 * t257 - t255 * t301 + t308) * t343 + (-t219 * t255 + t257 * t301 + t306) * t342;
t31 = t171 + (-t212 - t423) * t382 + t276 * t257 - t371;
t46 = -rSges(5,2) * t348 + t255 * t275 + t193 + t368;
t45 = (-t249 - t425) * t382 + t275 * t257 + t369;
t270 = cos(qJ(1));
t261 = t270 * pkin(1);
t259 = t268 * pkin(1);
t250 = t270 * t418;
t235 = pkin(3) * t377;
t218 = pkin(3) * t344;
t181 = t203 + t261;
t180 = t259 - t451;
t167 = -rSges(4,3) * t257 + t325;
t162 = -Icges(4,3) * t255 - t257 * t315;
t161 = -Icges(4,3) * t257 + t290;
t154 = t179 + t250;
t153 = t178 - t354;
t151 = t201 * t257 + t235;
t126 = t255 * t148;
t124 = t255 * t131;
t115 = t235 + t128;
t114 = t300 * t255;
t113 = t123 + t261;
t112 = t259 + t122;
t106 = t109 + t261;
t105 = t108 + t259;
t102 = t326 + t230 + (t237 + (-pkin(2) + t249) * t255) * t263;
t95 = -t255 * t444 + t315 * t378;
t94 = t257 * t444 + t263 * t290;
t92 = t101 + t261;
t91 = t100 + t259;
t88 = rSges(5,2) * t347 + (t254 * t379 + t256 * t382) * rSges(5,1) - t369;
t85 = -t149 * t257 + t126;
t78 = pkin(3) * t281 - t201 * t378 - t397;
t77 = t177 * t257 + t218 + t455;
t75 = -t140 * t257 + t125;
t60 = t250 + t74;
t59 = t73 - t354;
t58 = pkin(4) * t282 - t190 * t378 - t398;
t57 = t127 * t263 + t374;
t50 = -t162 * t255 + t305 * t257;
t49 = -t161 * t255 + t285;
t48 = -t162 * t257 - t454;
t47 = -t161 * t257 - t307 * t255;
t44 = t257 * t283 + (-pkin(3) * t359 - pkin(4) * t380 - t160) * t255;
t43 = t255 * t283 + t218 + t374;
t42 = t250 + t46;
t41 = t45 - t354;
t30 = t250 + t32;
t29 = t31 - t354;
t28 = t124 + t126 + (-t132 - t149) * t257;
t25 = t257 * t375 + t376;
t22 = -t257 * t88 + t351;
t16 = -t257 * t72 + t352;
t14 = t124 + (-t132 + t375) * t257 + t376;
t7 = (-t102 - t88) * t257 + t351 + t353;
t6 = t257 * t427 + t299;
t5 = (-t102 + t427) * t257 + t299 + t353;
t4 = (t257 * t80 + (-t38 + t286) * t263) * t257 + (t37 * t263 + (t145 * t380 + t147 * t383 + t254 * t81 - t256 * t83 - t399) * t255 + (-t400 + t79 + (-t147 * t263 - t335) * t256 + (t145 * t263 + t333) * t254) * t257) * t255;
t9 = [(t29 * t92 + t30 * t91) * t435 + (t105 * t42 + t106 * t41) * t436 + (t112 * t60 + t113 * t59) * t437 + (t153 * t181 + t154 * t180) * t438 + t273; m(6) * (t100 * t30 + t101 * t29 + t31 * t92 + t32 * t91) + m(5) * (t105 * t46 + t106 * t45 + t108 * t42 + t109 * t41) + m(4) * (t112 * t74 + t113 * t73 + t122 * t60 + t123 * t59) + m(3) * (t153 * t203 - t154 * t451 + t178 * t181 + t179 * t180) + t273; (t100 * t32 + t101 * t31) * t435 + (t108 * t46 + t109 * t45) * t436 + (t122 * t74 + t123 * t73) * t437 + (t178 * t203 - t179 * t451) * t438 + t273; m(5) * (t105 * t77 + t106 * t78 + t150 * t41 + t151 * t42) + m(6) * (t114 * t29 + t115 * t30 + t43 * t91 + t44 * t92) + ((-t113 * t263 + t60) * t257 + (-t112 * t263 - t59) * t255) * t431 + (t112 * t257 - t113 * t255) * t432 + t272; ((-t123 * t263 + t74) * t257 + (-t122 * t263 - t73) * t255) * t431 + (t122 * t257 - t123 * t255) * t432 + t272 + m(5) * (t108 * t77 + t109 * t78 + t150 * t45 + t151 * t46) + m(6) * (t100 * t43 + t101 * t44 + t114 * t31 + t115 * t32); -t257 * t2 + (t114 * t44 + t115 * t43 + t14 * t5) * t435 + t15 - t257 * t4 - t255 * t3 + (t150 * t78 + t151 * t77 + t28 * t7) * t436 + ((t167 * t255 - t168 * t257) * ((t263 * t167 - t257 * t297 + t367) * t257 + (-t255 * t297 + (t168 + (-t422 - t426) * t257) * t263 + t366) * t255) + t362 * t232 * t209) * t437 - t255 * ((t255 * t94 + (t49 + t454) * t263) * t255 + (-t50 * t263 + (-t163 * t359 - t165 * t360) * t257 + (-t306 * qJD(3) + t307 * t263 + t95) * t255) * t257) + (-t255 * t48 - t47 * t257) * t382 - t257 * ((t257 * t95 + (-t48 + t285) * t263) * t257 + (t47 * t263 + (t164 * t359 + t166 * t360) * t255 + (t308 * qJD(3) + t263 * t305 + t94) * t257) * t255) + t346 + (t255 * t50 + t257 * t49 + t456) * t378; m(6) * (t127 * t29 + t128 * t30 + t57 * t91 + t58 * t92) + ((-t106 * t263 + t42) * t257 + (-t105 * t263 - t41) * t255) * t429 + (t105 * t257 - t106 * t255) * t430 + t274; m(6) * (t100 * t57 + t101 * t58 + t127 * t31 + t128 * t32) + ((-t109 * t263 + t46) * t257 + (-t108 * t263 - t45) * t255) * t429 + (t108 * t257 - t109 * t255) * t430 + t274; m(6) * (t114 * t58 + t115 * t57 + t127 * t44 + t128 * t43 + t14 * t6 + t25 * t5) + m(5) * (-t150 * t397 + t22 * t28 + t7 * t85 + (-t151 * t382 - t255 * t78) * t201) + (-t4 - t263 * t17 + m(5) * (t151 * t177 + (t77 - t455) * t201) + t341) * t257 + t322; (t177 * t201 * t362 + t22 * t85) * t436 + (t127 * t58 + t128 * t57 + t25 * t6) * t435 + (t263 * t456 - t2 - t4) * t257 + t322; m(6) * ((-t255 * t92 + t257 * t91) * t160 + ((-t263 * t92 + t30) * t257 + (-t263 * t91 - t29) * t255) * t190) + t284; m(6) * ((t100 * t257 - t101 * t255) * t160 + ((-t101 * t263 + t32) * t257 + (-t100 * t263 - t31) * t255) * t190) + t284; m(6) * (-t114 * t398 - t115 * t350 + t14 * t16 - t44 * t393 + t5 * t75) + (m(6) * (-t114 * t392 + t115 * t160 + t190 * t43) + t341) * t257 + t346; m(6) * (-t127 * t398 - t128 * t350 + t16 * t25 - t58 * t393 + t6 * t75) + (m(6) * (-t127 * t392 + t128 * t160 + t190 * t57) + t341) * t257 + t346; (t160 * t190 * t362 + t16 * t75) * t435 + t341 * t257 + t346;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
