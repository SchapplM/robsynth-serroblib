% Calculate time derivative of joint inertia matrix for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:41
% EndTime: 2019-12-31 21:50:56
% DurationCPUTime: 8.41s
% Computational Cost: add. (13500->523), mult. (10723->704), div. (0->0), fcn. (8170->8), ass. (0->309)
t250 = qJ(3) + qJ(4);
t245 = cos(t250);
t243 = sin(t250);
t385 = Icges(6,5) * t243;
t189 = -Icges(6,3) * t245 + t385;
t387 = Icges(5,4) * t243;
t192 = Icges(5,2) * t245 + t387;
t447 = t189 - t192;
t384 = Icges(6,5) * t245;
t193 = Icges(6,1) * t243 - t384;
t386 = Icges(5,4) * t245;
t194 = Icges(5,1) * t243 + t386;
t446 = t193 + t194;
t437 = rSges(6,1) + pkin(4);
t445 = rSges(6,3) + qJ(5);
t251 = qJ(1) + qJ(2);
t244 = sin(t251);
t246 = cos(t251);
t298 = Icges(6,3) * t243 + t384;
t276 = t298 * t246;
t126 = Icges(6,6) * t244 + t276;
t302 = -Icges(5,2) * t243 + t386;
t280 = t302 * t246;
t132 = Icges(5,6) * t244 + t280;
t444 = t126 - t132;
t304 = Icges(6,1) * t245 + t385;
t282 = t304 * t246;
t134 = Icges(6,4) * t244 + t282;
t305 = Icges(5,1) * t245 - t387;
t283 = t305 * t246;
t136 = Icges(5,5) * t244 + t283;
t443 = t134 + t136;
t248 = qJD(3) + qJD(4);
t442 = (-t298 + t302) * t248;
t441 = (t304 + t305) * t248;
t190 = Icges(5,5) * t243 + Icges(5,6) * t245;
t191 = Icges(6,4) * t243 - Icges(6,6) * t245;
t431 = t190 + t191;
t428 = t447 * t243 + t446 * t245;
t125 = -Icges(6,6) * t246 + t244 * t298;
t131 = -Icges(5,6) * t246 + t244 * t302;
t440 = t125 - t131;
t133 = -Icges(6,4) * t246 + t244 * t304;
t135 = -Icges(5,5) * t246 + t244 * t305;
t439 = t133 + t135;
t232 = t244 * rSges(6,2);
t368 = t245 * t246;
t438 = pkin(4) * t368 + t232;
t299 = Icges(5,5) * t245 - Icges(5,6) * t243;
t127 = -Icges(5,3) * t246 + t244 * t299;
t301 = Icges(6,4) * t245 + Icges(6,6) * t243;
t129 = -Icges(6,2) * t246 + t244 * t301;
t297 = t125 * t243 + t133 * t245;
t415 = t246 * t297;
t295 = t131 * t243 - t135 * t245;
t416 = t246 * t295;
t436 = -t415 + t416 + (-t127 - t129) * t244;
t277 = t299 * t246;
t128 = Icges(5,3) * t244 + t277;
t279 = t301 * t246;
t130 = Icges(6,2) * t244 + t279;
t294 = t132 * t243 - t136 * t245;
t296 = t126 * t243 + t134 * t245;
t435 = (-t294 + t296) * t246 + (t128 + t130) * t244;
t235 = t246 * rSges(6,2);
t310 = rSges(6,1) * t245 + rSges(6,3) * t243;
t434 = -t235 + (pkin(4) * t245 + qJ(5) * t243 + t310) * t244;
t372 = t243 * t246;
t433 = rSges(6,1) * t368 + t445 * t372 + t438;
t371 = t243 * t248;
t338 = t244 * t371;
t432 = t437 * t338;
t430 = t437 * t243 - t445 * t245;
t249 = qJD(1) + qJD(2);
t373 = t194 * t248;
t374 = t193 * t248;
t375 = t192 * t248;
t376 = t189 * t248;
t429 = t431 * t249 + (t376 - t375 + t441) * t245 + (-t374 - t373 - t442) * t243;
t268 = Icges(6,6) * t249 - t376;
t370 = t244 * t249;
t76 = t246 * t268 - t298 * t370;
t269 = Icges(5,6) * t249 - t375;
t82 = t246 * t269 - t302 * t370;
t427 = t248 * t443 - t76 + t82;
t272 = Icges(6,4) * t249 - t374;
t84 = t246 * t272 - t304 * t370;
t271 = Icges(5,5) * t249 - t373;
t86 = t246 * t271 - t305 * t370;
t426 = t248 * t444 + t84 + t86;
t254 = cos(qJ(3));
t346 = qJD(3) * t254;
t252 = sin(qJ(3));
t363 = t249 * t252;
t425 = -t244 * t346 - t246 * t363;
t367 = t245 * t248;
t335 = t244 * t367;
t365 = t246 * t249;
t424 = t243 * t365 + t335;
t366 = t246 * t248;
t333 = t245 * t366;
t345 = qJD(5) * t243;
t423 = rSges(6,2) * t365 + t246 * t345 + t445 * t333;
t393 = rSges(4,2) * t252;
t395 = rSges(4,1) * t254;
t422 = -t393 + t395;
t421 = t428 * t249 + (-t299 - t301) * t248;
t388 = Icges(4,4) * t254;
t303 = -Icges(4,2) * t252 + t388;
t281 = t303 * t246;
t154 = Icges(4,6) * t244 + t281;
t389 = Icges(4,4) * t252;
t306 = Icges(4,1) * t254 - t389;
t284 = t306 * t246;
t156 = Icges(4,5) * t244 + t284;
t290 = t154 * t252 - t156 * t254;
t420 = t244 * t290;
t419 = t244 * t294;
t418 = t244 * t296;
t153 = -Icges(4,6) * t246 + t244 * t303;
t155 = -Icges(4,5) * t246 + t244 * t306;
t292 = t153 * t252 - t155 * t254;
t417 = t246 * t292;
t230 = t244 * rSges(5,3);
t414 = rSges(5,1) * t368 + t230;
t220 = Icges(4,2) * t254 + t389;
t221 = Icges(4,1) * t252 + t388;
t287 = t220 * t252 - t221 * t254;
t300 = Icges(4,5) * t254 - Icges(4,6) * t252;
t411 = t300 * qJD(3) + t249 * t287;
t256 = -pkin(8) - pkin(7);
t364 = t246 * t256;
t240 = pkin(3) * t254 + pkin(2);
t396 = pkin(2) - t240;
t410 = t244 * t396 - t364;
t409 = 2 * m(3);
t408 = 2 * m(4);
t407 = 2 * m(5);
t406 = 2 * m(6);
t405 = t244 / 0.2e1;
t404 = -t246 / 0.2e1;
t205 = t422 * qJD(3);
t402 = m(4) * t205;
t227 = rSges(4,1) * t252 + rSges(4,2) * t254;
t401 = m(4) * t227;
t392 = rSges(5,2) * t243;
t394 = rSges(5,1) * t245;
t173 = (-t392 + t394) * t248;
t400 = m(5) * t173;
t197 = rSges(5,1) * t243 + rSges(5,2) * t245;
t399 = m(5) * t197;
t398 = pkin(3) * t252;
t237 = t244 * pkin(7);
t253 = sin(qJ(1));
t397 = t253 * pkin(1);
t391 = pkin(1) * qJD(1);
t231 = t244 * rSges(4,3);
t312 = -t430 - t398;
t106 = t312 * t246;
t383 = t106 * t249;
t111 = t430 * t246;
t382 = t111 * t249;
t381 = t127 * t249;
t380 = t128 * t249;
t379 = t129 * t249;
t378 = t130 * t249;
t377 = t173 * t244;
t369 = t244 * t254;
t362 = t249 * t256;
t238 = t246 * pkin(7);
t122 = t238 - t410;
t313 = t246 * t240 - t244 * t256;
t349 = -t246 * pkin(2) - t237;
t123 = t313 + t349;
t361 = t244 * t122 + t246 * t123;
t351 = t246 * rSges(5,3) + t244 * t392;
t138 = t244 * t394 - t351;
t140 = -rSges(5,2) * t372 + t414;
t88 = t244 * t138 + t246 * t140;
t359 = -qJ(5) * t371 - (pkin(4) * t248 - qJD(5)) * t245 - t310 * t248;
t358 = t430 * t370;
t337 = t243 * t370;
t355 = rSges(5,2) * t337 + rSges(5,3) * t365;
t334 = t244 * t363;
t353 = rSges(4,2) * t334 + rSges(4,3) * t365;
t347 = qJD(3) * t252;
t329 = t244 * t347;
t352 = pkin(3) * t329 + t244 * t362;
t350 = t246 * rSges(4,3) + t244 * t393;
t348 = t244 ^ 2 + t246 ^ 2;
t226 = pkin(7) * t365;
t327 = t246 * t347;
t344 = t122 * t365 + t244 * ((-t246 * t396 - t237) * t249 - t352) + t246 * (-pkin(3) * t327 + t410 * t249 - t226);
t275 = -t243 * t366 - t245 * t370;
t331 = -rSges(5,1) * t338 - rSges(5,2) * t424;
t343 = t138 * t365 + t244 * (t414 * t249 + t331) + t246 * (rSges(5,1) * t275 - rSges(5,2) * t333 + t355);
t255 = cos(qJ(1));
t342 = t255 * t391;
t341 = pkin(3) * t347;
t340 = pkin(3) * t346;
t339 = t253 * t391;
t330 = -rSges(4,1) * t329 + rSges(4,2) * t425;
t326 = t370 / 0.2e1;
t325 = t365 / 0.2e1;
t324 = -pkin(2) - t395;
t323 = -t197 - t398;
t322 = -t240 - t394;
t199 = t246 * rSges(3,1) - rSges(3,2) * t244;
t85 = t244 * t272 + t249 * t282;
t321 = t125 * t248 + t85;
t87 = t244 * t271 + t249 * t283;
t319 = -t131 * t248 + t87;
t77 = t244 * t268 + t249 * t276;
t317 = t133 * t248 - t77;
t83 = t244 * t269 + t249 * t280;
t315 = t135 * t248 + t83;
t33 = t434 * t244 + t433 * t246;
t175 = -rSges(3,1) * t365 + rSges(3,2) * t370;
t34 = -t129 * t246 + t244 * t297;
t35 = -t130 * t246 + t418;
t36 = -t127 * t246 - t244 * t295;
t37 = -t128 * t246 - t419;
t267 = Icges(5,3) * t249 - t190 * t248;
t78 = t246 * t267 - t299 * t370;
t79 = t244 * t267 + t249 * t277;
t270 = Icges(6,2) * t249 - t191 * t248;
t80 = t246 * t270 - t301 * t370;
t81 = t244 * t270 + t249 * t279;
t311 = ((-t34 - t36) * t370 + t436 * t365) * t246 + (((t80 + t78) * t244 + (-t418 + t419 - t436) * t249) * t244 + (t35 + t37) * t370 + t435 * t365 + ((t378 - t81 + t380 - t79) * t244 + (-t367 * t440 + t439 * t371 - t379 - t381) * t246 + t435 * t249 + ((t249 * t439 + t426) * t244 + (-t85 - t87) * t246) * t245 + ((t249 * t440 - t427) * t244 + (-t77 + t83) * t246) * t243) * t246) * t244;
t198 = -rSges(3,1) * t244 - rSges(3,2) * t246;
t308 = t322 * t244;
t219 = Icges(4,5) * t252 + Icges(4,6) * t254;
t293 = t153 * t254 + t155 * t252;
t291 = t154 * t254 + t156 * t252;
t286 = t434 * t365 + (t437 * t275 - t337 * t445 + t423) * t246 + (t424 * qJ(5) + rSges(6,3) * t335 + t244 * t345 + (t246 * t310 + t438) * t249 - t432) * t244;
t158 = t246 * t422 + t231;
t285 = -t340 + t359;
t174 = t198 * t249;
t278 = t300 * t246;
t115 = t158 - t349;
t3 = (t246 * t81 + (t35 - t415) * t249) * t246 + (t34 * t249 + (t126 * t367 - t134 * t371 + t243 * t76 + t245 * t84 + t378) * t244 + (-t379 - t80 + (t134 * t249 - t321) * t245 + (t126 * t249 + t317) * t243) * t246) * t244;
t4 = (t246 * t79 + (t37 + t416) * t249) * t246 + (t36 * t249 + (-t132 * t367 - t136 * t371 - t243 * t82 + t245 * t86 + t380) * t244 + (-t381 - t78 + (t136 * t249 - t319) * t245 + (-t132 * t249 + t315) * t243) * t246) * t244;
t274 = (-t4 - t3) * t246 + t311;
t114 = t244 * t324 + t238 + t350;
t75 = t313 + t433;
t273 = -t243 * t445 - t245 * t437 - t240;
t104 = t140 + t313;
t266 = Icges(4,5) * t249 - qJD(3) * t221;
t265 = Icges(4,6) * t249 - qJD(3) * t220;
t264 = Icges(4,3) * t249 - qJD(3) * t219;
t103 = t308 + t351 - t364;
t263 = t244 * t273 - t364;
t262 = (t426 * t243 - t421 * t244 + t427 * t245 + t429 * t246) * t405 + (t421 * t246 + (t315 + t317) * t245 + t429 * t244 + (t319 + t321) * t243) * t404 + (t439 * t243 + t428 * t244 - t245 * t440 - t431 * t246) * t326 + (t443 * t243 + t431 * t244 - t245 * t444 + t428 * t246) * t325;
t63 = (t324 * t246 + (-rSges(4,3) - pkin(7)) * t244) * t249 - t330;
t74 = t235 + t263;
t201 = t303 * qJD(3);
t202 = t306 * qJD(3);
t259 = -t201 * t252 + t202 * t254 + t219 * t249 + (-t220 * t254 - t221 * t252) * qJD(3);
t45 = (t246 * t322 - t230) * t249 - t331 + t352;
t258 = t254 * t201 + t252 * t202 - t220 * t347 + t221 * t346 + t441 * t243 + t442 * t245 + t446 * t367 + t447 * t371;
t62 = -rSges(4,2) * t246 * t346 - pkin(2) * t370 + t226 + (-t249 * t369 - t327) * rSges(4,1) + t353;
t257 = t262 + (-qJD(3) * t290 + t411 * t244 + t259 * t246 + t252 * (t246 * t266 - t306 * t370) + t254 * (t246 * t265 - t303 * t370)) * t405 + (-qJD(3) * t292 + (t244 * t266 + t249 * t284) * t252 + t259 * t244 - t411 * t246 + t254 * (t244 * t265 + t249 * t281)) * t404 + (-t219 * t246 - t244 * t287 + t293) * t326 + (t219 * t244 - t246 * t287 + t291) * t325;
t44 = t249 * t308 + (-t197 * t248 - t341 - t362) * t246 + t355;
t27 = (-t367 * t445 - t345) * t244 + (t246 * t273 - t232) * t249 + t352 + t432;
t26 = (-t371 * t437 - t341) * t246 + t263 * t249 + t423;
t247 = t255 * pkin(1);
t206 = pkin(3) * t334;
t177 = t199 + t247;
t176 = t198 - t397;
t157 = rSges(4,1) * t369 - t350;
t152 = Icges(4,3) * t244 + t278;
t151 = -Icges(4,3) * t246 + t244 * t300;
t147 = t175 - t342;
t146 = t174 - t339;
t142 = t323 * t246;
t141 = t323 * t244;
t110 = t430 * t244;
t108 = t115 + t247;
t107 = t114 - t397;
t105 = t312 * t244;
t102 = t104 + t247;
t101 = t103 - t397;
t96 = t244 * t264 + t249 * t278;
t95 = t246 * t264 - t300 * t370;
t73 = pkin(3) * t425 - t197 * t365 - t377;
t72 = t197 * t370 + t206 + (-t173 - t340) * t246;
t67 = t247 + t75;
t66 = t74 - t397;
t53 = t63 - t342;
t52 = t62 - t339;
t51 = t152 * t244 - t290 * t246;
t50 = t151 * t244 - t417;
t49 = -t152 * t246 - t420;
t48 = -t151 * t246 - t292 * t244;
t47 = t244 * t359 - t382;
t46 = t246 * t359 + t358;
t43 = t45 - t342;
t42 = t44 - t339;
t32 = t244 * t285 + t383;
t31 = t246 * t285 + t206 + t358;
t30 = t88 + t361;
t25 = t27 - t342;
t24 = t26 - t339;
t21 = t33 + t361;
t20 = -t140 * t370 + t343;
t7 = (-t123 - t140) * t370 + t343 + t344;
t6 = -t370 * t433 + t286;
t5 = (-t123 - t433) * t370 + t286 + t344;
t1 = [t258 + (t24 * t67 + t25 * t66) * t406 + (t101 * t43 + t102 * t42) * t407 + (t107 * t53 + t108 * t52) * t408 + (t146 * t177 + t147 * t176) * t409; t258 + m(6) * (t24 * t75 + t25 * t74 + t26 * t67 + t27 * t66) + m(5) * (t101 * t45 + t102 * t44 + t103 * t43 + t104 * t42) + m(4) * (t107 * t63 + t108 * t62 + t114 * t53 + t115 * t52) + m(3) * (t146 * t199 + t147 * t198 + t174 * t177 + t175 * t176); t258 + (t26 * t75 + t27 * t74) * t406 + (t103 * t45 + t104 * t44) * t407 + (t114 * t63 + t115 * t62) * t408 + (t174 * t199 + t175 * t198) * t409; t257 + m(6) * (t105 * t24 + t106 * t25 + t31 * t66 + t32 * t67) + m(5) * (t101 * t72 + t102 * t73 + t141 * t42 + t142 * t43) + ((-t108 * t249 - t53) * t246 + (t107 * t249 - t52) * t244) * t401 + (-t107 * t246 - t108 * t244) * t402; t257 + m(6) * (t105 * t26 + t106 * t27 + t31 * t74 + t32 * t75) + m(5) * (t103 * t72 + t104 * t73 + t141 * t44 + t142 * t45) + ((-t115 * t249 - t63) * t246 + (t114 * t249 - t62) * t244) * t401 + (-t114 * t246 - t115 * t244) * t402; (t105 * t32 + t106 * t31 + t21 * t5) * t406 - t246 * t3 - t246 * t4 + (t141 * t73 + t142 * t72 + t30 * t7) * t407 + (t244 * t51 - t246 * t50) * t365 + t244 * ((t244 * t95 + (t50 + t420) * t249) * t244 + (t51 * t249 + (t153 * t346 + t155 * t347) * t246 + (-t291 * qJD(3) - t249 * t292 - t96) * t244) * t246) + (t244 * t49 - t246 * t48) * t370 - t246 * ((t246 * t96 + (t49 + t417) * t249) * t246 + (t48 * t249 + (-t154 * t346 - t156 * t347) * t244 + (t293 * qJD(3) - t249 * t290 - t95) * t246) * t244) + ((t157 * t244 + t158 * t246) * (((-t158 + t231) * t249 + t330) * t244 + (-qJD(3) * t227 * t246 + t249 * t157 + t353) * t246) + t348 * t227 * t205) * t408 + t311; t262 + m(6) * (-t110 * t24 - t111 * t25 + t46 * t66 + t47 * t67) + ((-t102 * t249 - t43) * t246 + (t101 * t249 - t42) * t244) * t399 + (-t101 * t246 - t102 * t244) * t400; t262 + m(6) * (-t110 * t26 - t111 * t27 + t46 * t74 + t47 * t75) + ((-t104 * t249 - t45) * t246 + (t103 * t249 - t44) * t244) * t399 + (-t103 * t246 - t104 * t244) * t400; m(6) * (t105 * t47 + t106 * t46 - t110 * t32 - t111 * t31 + t21 * t6 + t33 * t5) + m(5) * (-t142 * t173 * t246 - t141 * t377 + t20 * t30 + t7 * t88) + ((-t141 * t249 - t72) * t246 + (t142 * t249 - t73) * t244) * t399 + t274; (t173 * t197 * t348 + t20 * t88) * t407 + (-t110 * t47 - t111 * t46 + t33 * t6) * t406 + t274; m(6) * ((t244 * t67 + t246 * t66) * t367 + ((t249 * t67 + t25) * t246 + (-t249 * t66 + t24) * t244) * t243); m(6) * ((t244 * t75 + t246 * t74) * t367 + ((t249 * t75 + t27) * t246 + (-t249 * t74 + t26) * t244) * t243); m(6) * ((-t5 + (t105 * t244 + t106 * t246) * t248) * t245 + (t21 * t248 + (t105 * t249 + t31) * t246 + (t32 - t383) * t244) * t243); m(6) * ((-t6 + (-t110 * t244 - t111 * t246) * t248) * t245 + (t248 * t33 + (-t110 * t249 + t46) * t246 + (t47 + t382) * t244) * t243); (-0.1e1 + t348) * t243 * t367 * t406;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
