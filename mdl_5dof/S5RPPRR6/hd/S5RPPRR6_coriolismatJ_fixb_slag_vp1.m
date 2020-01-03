% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:42
% EndTime: 2019-12-31 17:57:55
% DurationCPUTime: 10.40s
% Computational Cost: add. (42908->478), mult. (36917->720), div. (0->0), fcn. (40341->9), ass. (0->290)
t299 = pkin(9) + qJ(4);
t295 = sin(t299);
t300 = qJ(1) + pkin(8);
t296 = sin(t300);
t410 = t295 * t296;
t297 = cos(t299);
t302 = sin(qJ(5));
t304 = cos(qJ(5));
t334 = Icges(6,5) * t304 - Icges(6,6) * t302;
t233 = -Icges(6,3) * t297 + t295 * t334;
t420 = Icges(6,4) * t304;
t336 = -Icges(6,2) * t302 + t420;
t235 = -Icges(6,6) * t297 + t295 * t336;
t421 = Icges(6,4) * t302;
t337 = Icges(6,1) * t304 - t421;
t237 = -Icges(6,5) * t297 + t295 * t337;
t298 = cos(t300);
t400 = t298 * t304;
t407 = t296 * t302;
t255 = t297 * t407 + t400;
t401 = t298 * t302;
t406 = t296 * t304;
t256 = t297 * t406 - t401;
t129 = t233 * t410 - t235 * t255 + t237 * t256;
t177 = Icges(6,5) * t256 - Icges(6,6) * t255 + Icges(6,3) * t410;
t248 = Icges(6,4) * t256;
t180 = -Icges(6,2) * t255 + Icges(6,6) * t410 + t248;
t247 = Icges(6,4) * t255;
t184 = -Icges(6,1) * t256 - Icges(6,5) * t410 + t247;
t96 = t177 * t410 - t180 * t255 - t184 * t256;
t257 = -t297 * t401 + t406;
t258 = t297 * t400 + t407;
t409 = t295 * t298;
t179 = Icges(6,5) * t258 + Icges(6,6) * t257 + Icges(6,3) * t409;
t422 = Icges(6,4) * t258;
t182 = Icges(6,2) * t257 + Icges(6,6) * t409 + t422;
t249 = Icges(6,4) * t257;
t185 = Icges(6,1) * t258 + Icges(6,5) * t409 + t249;
t97 = t179 * t410 - t255 * t182 + t256 * t185;
t340 = t296 * t96 + t298 * t97;
t12 = t297 * t129 - t295 * t340;
t342 = t258 * rSges(6,1) + t257 * rSges(6,2);
t438 = -pkin(6) - qJ(3);
t442 = cos(qJ(1)) * pkin(1);
t352 = -t296 * t438 + t442;
t292 = cos(pkin(9)) * pkin(3) + pkin(2);
t441 = t297 * pkin(4);
t454 = rSges(6,3) + pkin(7);
t477 = t454 * t295 + t292 + t441;
t153 = t298 * t477 + t342 + t352;
t206 = -rSges(6,1) * t255 - rSges(6,2) * t256;
t207 = rSges(6,1) * t257 - rSges(6,2) * t258;
t440 = sin(qJ(1)) * pkin(1);
t324 = -t298 * t438 - t440;
t481 = -t256 * rSges(6,1) + t255 * rSges(6,2);
t489 = -t296 * t477 + t324 + t481;
t506 = m(6) * (t153 * t207 - t206 * t489);
t432 = m(6) * qJD(5);
t131 = t233 * t409 + t235 * t257 + t237 * t258;
t98 = t177 * t409 + t257 * t180 - t184 * t258;
t99 = t179 * t409 + t257 * t182 + t258 * t185;
t339 = t296 * t98 + t298 * t99;
t504 = -t297 * t131 + t295 * t339;
t48 = t296 * t97 - t298 * t96;
t189 = rSges(6,3) * t409 + t342;
t430 = rSges(6,1) * t304;
t341 = -rSges(6,2) * t302 + t430;
t239 = -rSges(6,3) * t297 + t295 * t341;
t156 = t189 * t297 + t239 * t409;
t187 = rSges(6,3) * t410 - t481;
t500 = t187 * t297 + t239 * t410;
t482 = t156 * t296 - t298 * t500;
t415 = t177 * t297;
t499 = t180 * t302 + t184 * t304;
t111 = t499 * t295 + t415;
t273 = pkin(4) * t295 - pkin(7) * t297;
t375 = t239 + t273;
t195 = t375 * t296;
t197 = t375 * t298;
t479 = t195 * t296 + t197 * t298;
t271 = rSges(5,1) * t295 + rSges(5,2) * t297;
t293 = t296 ^ 2;
t294 = t298 ^ 2;
t485 = (t293 + t294) * t271;
t498 = -m(6) / 0.2e1;
t389 = t479 * t498 - m(5) * t485 / 0.2e1;
t175 = (-t454 * t297 + (pkin(4) + t341) * t295) * t296;
t403 = t297 * t298;
t281 = pkin(7) * t403;
t372 = t295 * rSges(6,2) * t401 + rSges(6,3) * t403;
t176 = t281 + (-pkin(4) - t430) * t409 + t372;
t472 = m(6) / 0.2e1;
t250 = t271 * t296;
t251 = t271 * t298;
t190 = t250 * t296 + t251 * t298;
t496 = m(5) * t190;
t390 = (t175 * t296 - t176 * t298) * t472 + t496 / 0.2e1;
t45 = t390 - t389;
t502 = t45 * qJD(1);
t501 = t296 * t99 - t298 * t98;
t483 = -t153 * t296 - t298 * t489;
t408 = t296 * t297;
t228 = Icges(5,4) * t408 - Icges(5,2) * t410 - Icges(5,6) * t298;
t290 = Icges(5,4) * t297;
t418 = Icges(5,2) * t295;
t229 = Icges(5,6) * t296 + (t290 - t418) * t298;
t423 = Icges(5,4) * t295;
t269 = Icges(5,1) * t297 - t423;
t231 = Icges(5,5) * t296 + t269 * t298;
t211 = t231 * t408;
t265 = Icges(5,5) * t297 - Icges(5,6) * t295;
t227 = Icges(5,3) * t296 + t298 * t265;
t346 = t298 * t227 - t211;
t226 = Icges(5,5) * t408 - Icges(5,6) * t410 - Icges(5,3) * t298;
t277 = Icges(5,4) * t410;
t230 = Icges(5,1) * t408 - Icges(5,5) * t298 - t277;
t381 = -t296 * t226 - t230 * t403;
t490 = -t228 * t409 - t229 * t410 - t346 - t381;
t460 = t296 / 0.2e1;
t457 = -t298 / 0.2e1;
t455 = t298 / 0.2e1;
t486 = t239 * t296;
t260 = (-Icges(6,2) * t304 - t421) * t295;
t261 = (-Icges(6,1) * t302 - t420) * t295;
t476 = -(t237 / 0.2e1 + t260 / 0.2e1) * t302 + (t261 / 0.2e1 - t235 / 0.2e1) * t304;
t243 = -Icges(5,2) * t408 - t277;
t266 = Icges(5,2) * t297 + t423;
t244 = t266 * t298;
t424 = Icges(5,1) * t295;
t338 = -t290 - t424;
t245 = t338 * t296;
t246 = t338 * t298;
t475 = ((t228 - t245) * t298 + (-t229 + t246) * t296) * t297 + ((t230 + t243) * t298 + (-t231 + t244) * t296) * t295;
t474 = 0.2e1 * m(6);
t473 = 0.4e1 * qJD(1);
t471 = -t504 / 0.2e1;
t470 = t48 / 0.2e1;
t469 = t501 / 0.2e1;
t200 = -Icges(6,5) * t255 - Icges(6,6) * t256;
t384 = -Icges(6,2) * t256 - t184 - t247;
t386 = -Icges(6,1) * t255 - t180 - t248;
t84 = -t200 * t297 + (-t384 * t302 + t386 * t304) * t295;
t468 = t84 / 0.2e1;
t240 = t295 * rSges(6,3) + t297 * t341;
t117 = (t240 * t296 - t187) * t295;
t221 = -rSges(6,1) * t295 * t400 + t372;
t118 = (-t239 * t298 - t221) * t297 + (-t240 * t298 + t189) * t295;
t466 = m(6) * (t117 * t489 + t118 * t153 - t156 * t176 + t175 * t500);
t262 = (-rSges(6,1) * t302 - rSges(6,2) * t304) * t295;
t464 = m(6) * (-t195 * t207 + t197 * t206 + t483 * t262);
t462 = m(6) * (t153 * t176 + t175 * t489);
t461 = t295 / 0.2e1;
t459 = t296 / 0.4e1;
t458 = -t297 / 0.2e1;
t456 = -t298 / 0.4e1;
t425 = rSges(4,3) + qJ(3);
t453 = m(4) * ((t425 * t298 - t440) * t298 + (t425 * t296 + t442) * t296);
t431 = rSges(5,1) * t297;
t350 = t292 + t431;
t371 = rSges(5,2) * t410 + t298 * rSges(5,3);
t193 = -t350 * t296 + t324 + t371;
t349 = -rSges(5,2) * t409 + t296 * rSges(5,3);
t194 = t350 * t298 + t349 + t352;
t452 = m(5) * (t193 * t250 - t194 * t251);
t451 = m(5) * (t193 * t298 + t194 * t296);
t448 = m(6) * t483;
t447 = m(6) * t482;
t147 = t206 * t296 + t207 * t298;
t443 = m(6) * t147;
t439 = qJD(4) / 0.2e1;
t428 = t296 * t12;
t426 = t298 * t504;
t414 = t179 * t297;
t412 = t233 * t297;
t411 = t295 * t228;
t259 = (-Icges(6,5) * t302 - Icges(6,6) * t304) * t295;
t404 = t297 * t259;
t399 = t302 * t235;
t236 = Icges(6,6) * t295 + t297 * t336;
t398 = t302 * t236;
t397 = t304 * t237;
t238 = Icges(6,5) * t295 + t297 * t337;
t396 = t304 * t238;
t69 = t117 * t296 - t118 * t298;
t395 = t69 * qJD(3);
t394 = t69 * qJD(5);
t330 = t187 * t298 - t189 * t296;
t95 = t330 * t297 + (-t221 * t296 - t298 * t486) * t295;
t74 = (t95 / 0.4e1 - t147 / 0.4e1) * t474;
t393 = t74 * qJD(2);
t385 = Icges(6,1) * t257 - t182 - t422;
t383 = -Icges(6,2) * t258 + t185 + t249;
t380 = t296 * t227 + t231 * t403;
t377 = -t235 + t261;
t376 = t237 + t260;
t274 = t295 * pkin(7) + t441;
t374 = -t240 - t274;
t369 = qJD(1) * t295;
t368 = qJD(1) * t297;
t367 = qJD(5) * t295;
t366 = qJD(5) * t297;
t361 = t471 + t504 / 0.2e1;
t360 = t410 / 0.4e1;
t345 = t295 * t229 - t226;
t107 = -t376 * t255 + t377 * t256 + t259 * t410;
t108 = t376 * t257 + t377 * t258 + t259 * t409;
t201 = Icges(6,5) * t257 - Icges(6,6) * t258;
t85 = -t201 * t297 + (-t383 * t302 + t385 * t304) * t295;
t344 = t464 / 0.2e1 + (t108 + t85) * t459 + (t107 + t84) * t456;
t335 = -Icges(5,5) * t295 - Icges(5,6) * t297;
t331 = -t182 * t302 + t185 * t304;
t112 = t295 * t331 - t414;
t333 = -t111 * t296 + t112 * t298;
t328 = t397 - t399;
t136 = -t229 * t409 + t380;
t327 = -t48 / 0.2e1 + t470 + (t298 * t345 + t136 - t380) * t455 + (-(-t297 * t230 + t411) * t296 - t298 * t226) * t457 + (t345 * t296 + t346 + t490) * t460;
t326 = t469 - t501 / 0.2e1 + t136 * t460 - t380 * t296 / 0.2e1 + (-t211 + (t227 + t411) * t298 + t381 + t490) * t457;
t70 = t200 * t410 - t384 * t255 + t386 * t256;
t71 = t201 * t410 - t383 * t255 + t385 * t256;
t32 = t296 * t71 - t298 * t70;
t72 = t200 * t409 + t384 * t257 + t386 * t258;
t73 = t201 * t409 + t383 * t257 + t385 * t258;
t33 = t296 * t73 - t298 * t72;
t325 = t32 * t457 + t33 * t460;
t321 = -t233 * t296 + t499;
t320 = -t233 * t298 - t331;
t234 = Icges(6,3) * t295 + t297 * t334;
t319 = t234 - t328;
t316 = t12 * t459 + t504 * t456 - t428 / 0.4e1 + t426 / 0.4e1 + (t360 - t410 / 0.4e1) * t501;
t113 = -t319 * t297 + (t233 + t396 - t398) * t295;
t157 = t295 * t328 - t412;
t216 = t235 * t296;
t218 = t237 * t296;
t78 = -t321 * t297 + (t216 * t302 - t218 * t304 + t177) * t295;
t217 = t235 * t298;
t219 = t237 * t298;
t79 = -t320 * t297 + (t217 * t302 - t219 * t304 + t179) * t295;
t309 = t295 * t319 + t412;
t92 = -t236 * t255 + t238 * t256 + t309 * t296;
t93 = t236 * t257 + t238 * t258 + t298 * t309;
t312 = t113 * t458 + t157 * t461 + t466 / 0.2e1 + (t78 + t92) * t360 + (t79 + t93) * t409 / 0.4e1 + (-t111 + t129) * t408 / 0.4e1 + (t112 + t131) * t403 / 0.4e1;
t311 = t295 * t321 + t415;
t310 = t295 * t320 + t414;
t308 = t397 / 0.2e1 - t399 / 0.2e1 - t234 / 0.2e1 + t290 + t424 / 0.2e1 - t418 / 0.2e1;
t307 = t396 / 0.2e1 - t398 / 0.2e1 + t233 / 0.2e1 + t269 / 0.2e1 - t266 / 0.2e1;
t272 = -rSges(5,2) * t295 + t431;
t242 = t335 * t298;
t241 = t335 * t296;
t198 = t374 * t298;
t196 = t374 * t296;
t160 = -t207 * t297 - t262 * t409;
t159 = t206 * t297 + t262 * t410;
t143 = -t443 / 0.2e1;
t139 = (t206 * t298 - t207 * t296) * t295;
t127 = t330 * t295;
t121 = -t404 + (-t376 * t302 + t377 * t304) * t295;
t120 = (-pkin(4) * t409 + t221 + t281) * t298 + (-t273 * t296 - t486) * t296;
t115 = (t274 * t296 + t187) * t296 + (t274 * t298 + t189) * t298;
t101 = -t447 / 0.2e1;
t75 = (t147 + t95) * t472;
t68 = m(6) * t69 * t439;
t67 = -t217 * t257 - t219 * t258 + t298 * t310;
t66 = -t216 * t257 - t218 * t258 + t298 * t311;
t65 = t217 * t255 - t219 * t256 + t310 * t296;
t64 = t216 * t255 - t218 * t256 + t311 * t296;
t58 = t115 * t147 + t262 * t479;
t57 = -t448 + t451 + t453;
t54 = t506 - t404 / 0.2e1 + t476 * t295;
t46 = t389 + t390;
t42 = -t157 * t297 + t295 * t333;
t31 = t101 + t443 / 0.2e1;
t30 = t143 + t101;
t29 = t143 + t447 / 0.2e1;
t28 = t295 * t307 + t297 * t308 + t452 + t462;
t27 = t296 * t67 - t298 * t66;
t26 = t296 * t65 - t298 * t64;
t23 = t117 * t500 - t118 * t156 + t127 * t95;
t22 = -t108 * t297 + (t296 * t72 + t298 * t73) * t295;
t21 = -t107 * t297 + (t296 * t70 + t298 * t71) * t295;
t14 = (-t113 + t333) * t297 + (t78 * t296 + t79 * t298 + t157) * t295;
t9 = (t339 - t93) * t297 + (t296 * t66 + t298 * t67 + t131) * t295;
t8 = (t340 - t92) * t297 + (t296 * t64 + t298 * t65 + t129) * t295;
t7 = m(6) * t58 + t325;
t6 = t361 * t410;
t5 = t327 * t296 + t298 * t326;
t4 = m(6) * t23 + (t426 / 0.2e1 - t428 / 0.2e1 - t14 / 0.2e1) * t297 + (t9 * t455 + t8 * t460 + t42 / 0.2e1) * t295;
t3 = t312 + t344;
t2 = (-t157 / 0.2e1 + (-t93 / 0.4e1 - t79 / 0.4e1) * t298 + (-t92 / 0.4e1 - t78 / 0.4e1) * t296) * t295 + t316 - t466 / 0.2e1 + (t113 / 0.2e1 + (-t131 / 0.4e1 - t112 / 0.4e1) * t298 + (-t129 / 0.4e1 + t111 / 0.4e1) * t296) * t297 + t344;
t1 = (-t108 / 0.4e1 - t85 / 0.4e1) * t296 + (t84 / 0.4e1 + t107 / 0.4e1) * t298 + t316 - t464 / 0.2e1 + t312;
t10 = [t57 * qJD(3) + t28 * qJD(4) + t54 * qJD(5), 0, qJD(1) * t57 + qJD(4) * t46 + qJD(5) * t30, t28 * qJD(1) + t46 * qJD(3) + t3 * qJD(5) + (m(6) * (t153 * t196 - t175 * t197 - t176 * t195 + t198 * t489) + (m(5) * (-t193 * t272 - t250 * t271) - t78 / 0.2e1 - t92 / 0.2e1 + t265 * t455 + (-t230 / 0.2e1 - t243 / 0.2e1) * t297 + (t228 / 0.2e1 - t245 / 0.2e1) * t295 - t326) * t298 + (m(5) * (-t194 * t272 + t251 * t271) + t79 / 0.2e1 + t93 / 0.2e1 + t265 * t460 + (t231 / 0.2e1 - t244 / 0.2e1) * t297 + (-t229 / 0.2e1 + t246 / 0.2e1) * t295 - t327) * t296) * qJD(4), t54 * qJD(1) + t30 * qJD(3) + t3 * qJD(4) - t121 * t366 + (t160 * t153 - t156 * t207 + t159 * t489 - t206 * t500) * t432 + ((t85 / 0.2e1 + t108 / 0.2e1) * t298 + (t468 + t107 / 0.2e1 - t361) * t296) * t367; 0, 0, 0, 0.2e1 * (-t496 / 0.2e1 + t120 * t472) * qJD(4) + t75 * qJD(5), t75 * qJD(4) + t139 * t432; t45 * qJD(4) + t29 * qJD(5) + (t448 / 0.4e1 - t451 / 0.4e1 - t453 / 0.4e1) * t473, 0, 0, t502 + ((-t196 * t298 + t198 * t296) * t439 + t394 / 0.4e1) * t474, t29 * qJD(1) + t68 + (t159 * t296 - t160 * t298) * t432; -t45 * qJD(3) + t5 * qJD(4) + t2 * qJD(5) + (-t462 / 0.4e1 - t452 / 0.4e1) * t473 - t308 * t368 - t307 * t369, -qJD(5) * t74, t394 * t498 - t502, t5 * qJD(1) + (m(5) * (t272 * t485 - (t296 * (rSges(5,1) * t408 - t371) + t298 * (rSges(5,1) * t403 + t349)) * t190) + m(6) * (t115 * t120 - t195 * t196 - t197 * t198) + (t293 * t242 + (-t296 * t241 + t475) * t298 + t27) * t460 + (t294 * t241 + (-t298 * t242 + t475) * t296 + t26) * t457) * qJD(4) + t7 * qJD(5), t2 * qJD(1) - t393 + t7 * qJD(4) + (t21 * t457 + t22 * t460) * qJD(5) + (t14 / 0.2e1 + (t468 + t471) * t298 + (-t85 / 0.2e1 + t12 / 0.2e1) * t296) * t366 + (-t395 / 0.2e1 + (t139 * t115 + t127 * t147 - t159 * t197 - t160 * t195 + t262 * t482 - t23) * qJD(5)) * m(6) + (-t42 / 0.2e1 + (t33 / 0.2e1 - t9 / 0.2e1) * t298 + (t32 / 0.2e1 - t8 / 0.2e1) * t296) * t367; t259 * t368 / 0.2e1 + t31 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) - qJD(1) * t506 - t476 * t369, qJD(4) * t74, qJD(1) * t31 + t68, t1 * qJD(1) + t393 + (t9 * t460 + t403 * t469 + t27 * t409 / 0.2e1 + t8 * t457 + t408 * t470 + t26 * t410 / 0.2e1 + (t111 * t298 + t112 * t296) * t461 + (t79 * t296 - t78 * t298) * t458 - t325) * qJD(4) + t4 * qJD(5) + (t395 / 0.2e1 + (t115 * t95 - t117 * t197 - t118 * t195 + t120 * t127 - t156 * t196 + t198 * t500 - t58) * qJD(4)) * m(6), t6 * qJD(1) + t4 * qJD(4) + (m(6) * (t127 * t139 - t156 * t160 + t159 * t500) + t297 ^ 2 * t121 / 0.2e1 + (t22 * t455 + t21 * t460 + (t296 * t84 + t298 * t85) * t458) * t295) * qJD(5);];
Cq = t10;
