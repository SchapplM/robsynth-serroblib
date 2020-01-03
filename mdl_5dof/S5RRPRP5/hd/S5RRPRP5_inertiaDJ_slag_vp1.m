% Calculate time derivative of joint inertia matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:02
% EndTime: 2019-12-31 19:54:24
% DurationCPUTime: 12.63s
% Computational Cost: add. (9981->591), mult. (10053->796), div. (0->0), fcn. (7796->8), ass. (0->339)
t247 = qJ(2) + pkin(8);
t236 = qJ(4) + t247;
t228 = cos(t236);
t227 = sin(t236);
t415 = Icges(6,5) * t227;
t184 = -Icges(6,3) * t228 + t415;
t421 = Icges(5,4) * t227;
t187 = Icges(5,2) * t228 + t421;
t481 = t184 - t187;
t414 = Icges(6,5) * t228;
t188 = Icges(6,1) * t227 - t414;
t420 = Icges(5,4) * t228;
t189 = Icges(5,1) * t227 + t420;
t480 = t188 + t189;
t440 = rSges(6,1) + pkin(4);
t479 = rSges(6,3) + qJ(5);
t252 = sin(qJ(1));
t254 = cos(qJ(1));
t302 = Icges(6,3) * t227 + t414;
t129 = Icges(6,6) * t252 + t254 * t302;
t307 = -Icges(5,2) * t227 + t420;
t135 = Icges(5,6) * t252 + t254 * t307;
t478 = t129 - t135;
t303 = Icges(5,5) * t228 - Icges(5,6) * t227;
t131 = Icges(5,3) * t252 + t254 * t303;
t306 = Icges(6,4) * t228 + Icges(6,6) * t227;
t133 = Icges(6,2) * t252 + t254 * t306;
t477 = t131 + t133;
t312 = Icges(6,1) * t228 + t415;
t137 = Icges(6,4) * t252 + t254 * t312;
t313 = Icges(5,1) * t228 - t421;
t139 = Icges(5,5) * t252 + t254 * t313;
t476 = t137 + t139;
t246 = qJD(2) + qJD(4);
t475 = (-t302 + t307) * t246;
t474 = (t312 + t313) * t246;
t185 = Icges(5,5) * t227 + Icges(5,6) * t228;
t186 = Icges(6,4) * t227 - Icges(6,6) * t228;
t468 = t185 + t186;
t465 = t481 * t227 + t480 * t228;
t128 = -Icges(6,6) * t254 + t252 * t302;
t134 = -Icges(5,6) * t254 + t252 * t307;
t473 = -t128 + t134;
t136 = -Icges(6,4) * t254 + t252 * t312;
t138 = -Icges(5,5) * t254 + t252 * t313;
t472 = t136 + t138;
t130 = -Icges(5,3) * t254 + t252 * t303;
t132 = -Icges(6,2) * t254 + t252 * t306;
t298 = t128 * t227 + t136 * t228;
t453 = t254 * t298;
t296 = t134 * t227 - t138 * t228;
t454 = t254 * t296;
t471 = -t453 + t454 + (-t130 - t132) * t252;
t295 = t135 * t227 - t139 * t228;
t297 = t129 * t227 + t137 * t228;
t470 = (-t295 + t297) * t254 + t477 * t252;
t242 = t252 * rSges(6,2);
t388 = t228 * t254;
t390 = t227 * t254;
t469 = t388 * t440 + t479 * t390 + t242;
t191 = rSges(6,1) * t227 - rSges(6,3) * t228;
t467 = pkin(4) * t227 - qJ(5) * t228 + t191;
t392 = t189 * t246;
t393 = t188 * t246;
t394 = t187 * t246;
t395 = t184 * t246;
t466 = (t395 - t394 + t474) * t228 + (-t393 - t392 - t475) * t227 + t468 * qJD(1);
t71 = -qJD(1) * t128 - t254 * t395;
t77 = -qJD(1) * t134 - t254 * t394;
t464 = -t246 * t476 + t71 - t77;
t79 = -qJD(1) * t136 - t254 * t393;
t81 = -qJD(1) * t138 - t254 * t392;
t463 = t246 * t478 + t79 + t81;
t462 = (-t303 - t306) * t246 + t465 * qJD(1);
t192 = rSges(5,1) * t227 + rSges(5,2) * t228;
t280 = t192 * t246;
t251 = sin(qJ(2));
t253 = cos(qJ(2));
t424 = Icges(3,4) * t253;
t311 = -Icges(3,2) * t251 + t424;
t174 = Icges(3,6) * t252 + t254 * t311;
t425 = Icges(3,4) * t251;
t317 = Icges(3,1) * t253 - t425;
t176 = Icges(3,5) * t252 + t254 * t317;
t291 = t174 * t251 - t176 * t253;
t461 = t252 * t291;
t234 = sin(t247);
t235 = cos(t247);
t422 = Icges(4,4) * t235;
t309 = -Icges(4,2) * t234 + t422;
t154 = Icges(4,6) * t252 + t254 * t309;
t423 = Icges(4,4) * t234;
t315 = Icges(4,1) * t235 - t423;
t156 = Icges(4,5) * t252 + t254 * t315;
t293 = t154 * t234 - t156 * t235;
t460 = t252 * t293;
t459 = t252 * t295;
t458 = t252 * t297;
t243 = t253 * pkin(2);
t229 = t243 + pkin(1);
t434 = pkin(1) - t229;
t457 = t252 * t434;
t173 = -Icges(3,6) * t254 + t252 * t311;
t175 = -Icges(3,5) * t254 + t252 * t317;
t292 = t173 * t251 - t175 * t253;
t456 = t254 * t292;
t153 = -Icges(4,6) * t254 + t252 * t309;
t155 = -Icges(4,5) * t254 + t252 * t315;
t294 = t153 * t234 - t155 * t235;
t455 = t254 * t294;
t240 = t252 * rSges(4,3);
t428 = rSges(4,2) * t234;
t452 = -t254 * t428 + t240;
t437 = pkin(2) * t251;
t328 = -pkin(3) * t234 - t437;
t285 = -t192 + t328;
t114 = t285 * t252;
t115 = t285 * t254;
t451 = t114 * t252 + t115 * t254;
t450 = qJD(1) * t130;
t449 = qJD(1) * t132;
t304 = Icges(4,5) * t235 - Icges(4,6) * t234;
t151 = -Icges(4,3) * t254 + t252 * t304;
t305 = Icges(3,5) * t253 - Icges(3,6) * t251;
t171 = -Icges(3,3) * t254 + t252 * t305;
t250 = -qJ(3) - pkin(6);
t245 = -pkin(7) + t250;
t368 = t245 - t250;
t436 = pkin(3) * t235;
t204 = t229 + t436;
t375 = t204 - t229;
t112 = t252 * t375 + t254 * t368;
t446 = 2 * m(3);
t445 = 2 * m(4);
t444 = 2 * m(5);
t443 = 2 * m(6);
t248 = t252 ^ 2;
t249 = t254 ^ 2;
t442 = t252 / 0.2e1;
t441 = -t254 / 0.2e1;
t218 = rSges(3,1) * t251 + rSges(3,2) * t253;
t439 = m(3) * t218;
t438 = m(5) * t192;
t435 = t252 * pkin(6);
t244 = t254 * pkin(6);
t433 = -pkin(6) - t250;
t432 = rSges(3,1) * t253;
t431 = rSges(4,1) * t235;
t430 = rSges(5,1) * t228;
t429 = rSges(3,2) * t251;
t427 = rSges(3,3) * t254;
t241 = t252 * rSges(3,3);
t239 = t252 * rSges(5,3);
t403 = t153 * t235;
t402 = t154 * t235;
t401 = t155 * t234;
t400 = t156 * t234;
t399 = t173 * t253;
t398 = t174 * t253;
t397 = t175 * t251;
t396 = t176 * t251;
t391 = t227 * t246;
t389 = t228 * t246;
t387 = t245 * t254;
t386 = t246 * t252;
t385 = t246 * t254;
t193 = t254 * t204;
t220 = t254 * t229;
t113 = -t252 * t368 + t193 - t220;
t150 = -pkin(1) * t254 + t252 * t433 + t220;
t384 = -t113 - t150;
t323 = rSges(6,1) * t228 + rSges(6,3) * t227;
t340 = pkin(4) * t246 - qJD(5);
t383 = -qJ(5) * t391 - t228 * t340 - t323 * t246;
t324 = -rSges(5,2) * t227 + t430;
t142 = -rSges(5,3) * t254 + t252 * t324;
t144 = rSges(5,1) * t388 - rSges(5,2) * t390 + t239;
t70 = t252 * t142 + t254 * t144;
t149 = t250 * t254 + t244 - t457;
t382 = t252 * t149 + t254 * t150;
t362 = qJD(1) * t252;
t380 = t467 * t362;
t197 = t328 * qJD(2);
t179 = t254 * t197;
t237 = qJD(3) * t252;
t379 = t179 + t237;
t348 = t228 * t385;
t355 = qJD(5) * t227;
t377 = qJ(5) * t348 + t254 * t355;
t361 = qJD(1) * t254;
t376 = rSges(6,2) * t361 + rSges(6,3) * t348;
t346 = t227 * t362;
t374 = rSges(5,2) * t346 + rSges(5,3) * t361;
t345 = t234 * t362;
t373 = rSges(4,2) * t345 + rSges(4,3) * t361;
t344 = t251 * t362;
t224 = pkin(2) * t344;
t372 = pkin(3) * t345 + t224;
t222 = t245 * t362;
t238 = qJD(3) * t254;
t371 = t222 + t238;
t358 = qJD(2) * t251;
t350 = pkin(2) * t358;
t370 = t250 * t362 + t252 * t350;
t369 = t254 * t432 + t241;
t367 = t248 + t249;
t366 = qJD(1) * t131;
t365 = qJD(1) * t133;
t152 = Icges(4,3) * t252 + t254 * t304;
t364 = qJD(1) * t152;
t172 = Icges(3,3) * t252 + t254 * t305;
t363 = qJD(1) * t172;
t360 = qJD(2) * t234;
t359 = qJD(2) * t235;
t357 = qJD(2) * t252;
t356 = qJD(2) * t253;
t349 = t227 * t385;
t264 = -t228 * t362 - t349;
t354 = t142 * t361 + t252 * (-t252 * t280 + (t254 * t324 + t239) * qJD(1)) + t254 * (rSges(5,1) * t264 - rSges(5,2) * t348 + t374);
t330 = t254 * t350;
t347 = -t238 - t370;
t353 = t149 * t361 + t252 * ((-t254 * t434 - t435) * qJD(1) + t347) + t254 * (-t330 + t237 + (t254 * t433 + t457) * qJD(1));
t352 = t254 * t429;
t199 = rSges(4,1) * t234 + rSges(4,2) * t235;
t341 = -t199 - t437;
t80 = qJD(1) * t137 - t252 * t393;
t339 = t128 * t246 + t80;
t82 = qJD(1) * t139 - t252 * t392;
t337 = t134 * t246 - t82;
t72 = qJD(1) * t129 - t252 * t395;
t335 = t136 * t246 - t72;
t78 = qJD(1) * t135 - t252 * t394;
t333 = t138 * t246 + t78;
t109 = t467 * t254;
t331 = -t252 * t245 + t193;
t329 = t252 * t112 + t254 * t113 + t382;
t141 = -rSges(6,2) * t254 + t252 * t323;
t167 = (pkin(4) * t228 + qJ(5) * t227) * t252;
t27 = t469 * t254 + (t141 + t167) * t252;
t32 = -t132 * t254 + t252 * t298;
t33 = -t133 * t254 + t458;
t34 = -t130 * t254 - t252 * t296;
t35 = -t131 * t254 - t459;
t272 = t246 * t185;
t73 = -t254 * t272 - t450;
t74 = -t252 * t272 + t366;
t273 = t246 * t186;
t75 = -t254 * t273 - t449;
t76 = -t252 * t273 + t365;
t327 = ((-t32 - t34) * t362 + t471 * t361) * t254 + (((t75 + t73) * t252 + (-t458 + t459 - t471) * qJD(1)) * t252 + (t33 + t35) * t362 + t470 * t361 + ((-t74 - t76) * t252 + (t389 * t473 + t472 * t391 - t449 - t450) * t254 + (t463 * t252 + (-t80 - t82) * t254) * t228 + (t464 * t252 + (-t72 + t78) * t254) * t227 + ((t298 - t296 + t477) * t252 + t470) * qJD(1)) * t254) * t252;
t326 = -t429 + t432;
t325 = -t428 + t431;
t261 = -t227 * t479 - t228 * t440 - t204;
t258 = t261 * t252;
t62 = (rSges(6,2) - t245) * t254 + t258;
t63 = t331 + t469;
t320 = t252 * t63 + t254 * t62;
t277 = t328 - t467;
t83 = t277 * t252;
t84 = t277 * t254;
t319 = t252 * t83 + t254 * t84;
t283 = -t204 - t324;
t95 = (rSges(5,3) - t245) * t254 + t283 * t252;
t96 = t144 + t331;
t318 = t252 * t96 + t254 * t95;
t316 = Icges(3,1) * t251 + t424;
t314 = Icges(4,1) * t234 + t422;
t310 = Icges(3,2) * t253 + t425;
t308 = Icges(4,2) * t235 + t423;
t108 = t467 * t252;
t299 = -t108 * t252 - t109 * t254;
t288 = t252 * (t252 * t197 + t361 * t375 - t222 + t370) + t254 * (-qJD(1) * t112 + t179 + t330) + t112 * t361 + t353;
t287 = t141 * t361 + t167 * t361 + t252 * ((pkin(4) * t361 + qJ(5) * t386) * t228 + (qJ(5) * t361 - t252 * t340) * t227) + t254 * (pkin(4) * t264 - qJ(5) * t346 + t377) + t252 * (-t191 * t386 + (t254 * t323 + t242) * qJD(1)) + t254 * (rSges(6,1) * t264 - rSges(6,3) * t346 + t376);
t164 = t254 * t431 + t452;
t286 = -pkin(1) - t326;
t284 = -t229 - t325;
t281 = (-t243 - t436) * qJD(2);
t279 = qJD(2) * t218;
t278 = qJD(2) * t199;
t270 = qJD(2) * t316;
t269 = qJD(2) * t314;
t268 = qJD(2) * t310;
t267 = qJD(2) * t308;
t266 = qJD(2) * (-Icges(3,5) * t251 - Icges(3,6) * t253);
t265 = qJD(2) * (-Icges(4,5) * t234 - Icges(4,6) * t235);
t166 = t324 * t246;
t263 = -t166 + t281;
t3 = (t254 * t76 + (t33 - t453) * qJD(1)) * t254 + (t32 * qJD(1) + (t129 * t389 - t137 * t391 + t227 * t71 + t228 * t79 + t365) * t252 + (-t75 - t339 * t228 + t335 * t227 + (-t132 + t297) * qJD(1)) * t254) * t252;
t4 = (t254 * t74 + (t35 + t454) * qJD(1)) * t254 + (t34 * qJD(1) + (-t135 * t389 - t139 * t391 - t227 * t77 + t228 * t81 + t366) * t252 + (-t73 + t337 * t228 + t333 * t227 + (-t130 - t295) * qJD(1)) * t254) * t252;
t262 = (-t4 - t3) * t254 + t327;
t260 = t281 + t383;
t259 = rSges(3,2) * t344 + rSges(3,3) * t361 - t254 * t279;
t257 = (t227 * t463 - t228 * t464 - t462 * t252 + t254 * t466) * t442 + (t462 * t254 + t466 * t252 + (t333 + t335) * t228 + (-t337 + t339) * t227) * t441 + (t472 * t227 + t228 * t473 + t465 * t252 - t468 * t254) * t362 / 0.2e1 + (t227 * t476 - t228 * t478 + t468 * t252 + t465 * t254) * t361 / 0.2e1;
t205 = t326 * qJD(2);
t183 = t325 * qJD(2);
t178 = -t352 + t369;
t177 = t252 * t326 - t427;
t163 = -rSges(4,3) * t254 + t252 * t325;
t148 = t341 * t254;
t147 = t341 * t252;
t127 = t435 + (pkin(1) - t429) * t254 + t369;
t126 = t252 * t286 + t244 + t427;
t111 = -t252 * t250 + t164 + t220;
t110 = (rSges(4,3) - t250) * t254 + t284 * t252;
t103 = t252 * t266 + t363;
t102 = -qJD(1) * t171 + t254 * t266;
t90 = t252 * t265 + t364;
t89 = -qJD(1) * t151 + t254 * t265;
t88 = t218 * t357 + ((-rSges(3,3) - pkin(6)) * t252 + t286 * t254) * qJD(1);
t87 = (t244 + (-pkin(1) - t432) * t252) * qJD(1) + t259;
t86 = -t199 * t361 - t252 * t183 + (-t251 * t361 - t252 * t356) * pkin(2);
t85 = t199 * t362 + t224 + (-pkin(2) * t356 - t183) * t254;
t55 = t252 * t172 - t291 * t254;
t54 = t252 * t171 - t456;
t53 = -t172 * t254 - t461;
t52 = -t171 * t254 - t252 * t292;
t47 = t199 * t357 + (t254 * t284 - t240) * qJD(1) - t347;
t46 = t237 + (-t229 - t431) * t362 + (-qJD(1) * t250 + qJD(2) * t341) * t254 + t373;
t45 = qJD(1) * t115 + t252 * t263;
t44 = t192 * t362 + t254 * t263 + t372;
t43 = t252 * t152 - t293 * t254;
t42 = t252 * t151 - t455;
t41 = -t152 * t254 - t460;
t40 = -t151 * t254 - t252 * t294;
t31 = -qJD(1) * t109 + t252 * t383;
t30 = t254 * t383 + t380;
t29 = (-t197 + t280) * t252 + (t254 * t283 - t239) * qJD(1) + t371;
t28 = -t254 * t280 + (-t387 + (-t204 - t430) * t252) * qJD(1) + t374 + t379;
t26 = qJD(1) * t84 + t252 * t260;
t25 = t254 * t260 + t372 + t380;
t24 = t261 * t361 + (-rSges(6,2) * qJD(1) - t355 - t197 + (t227 * t440 - t228 * t479) * t246) * t252 + t371;
t23 = -t440 * t349 + (t258 - t387) * qJD(1) + t376 + t377 + t379;
t22 = -t144 * t362 + t354;
t21 = t329 + t70;
t8 = t27 + t329;
t7 = -t362 * t469 + t287;
t6 = (-t144 + t384) * t362 + t288 + t354;
t5 = (-t469 + t384) * t362 + t287 + t288;
t1 = [(t28 * t96 + t29 * t95) * t444 + (t23 * t63 + t24 * t62) * t443 + (t110 * t47 + t111 * t46) * t445 + (t126 * t88 + t127 * t87) * t446 + t481 * t391 + t480 * t389 + (t315 - t308) * t360 + (t314 + t309) * t359 + (t317 - t310) * t358 + (t316 + t311) * t356 + t475 * t228 + t474 * t227; t257 + m(3) * ((-t252 * t87 - t254 * t88) * t218 + (-t126 * t254 - t127 * t252) * t205) + ((t402 / 0.2e1 + t400 / 0.2e1 - t127 * t439 + t398 / 0.2e1 + t396 / 0.2e1) * t254 + (t126 * t439 + t403 / 0.2e1 + t401 / 0.2e1 + t399 / 0.2e1 + t397 / 0.2e1) * t252) * qJD(1) + m(6) * (t23 * t83 + t24 * t84 + t25 * t62 + t26 * t63) + m(5) * (t114 * t28 + t115 * t29 + t44 * t95 + t45 * t96) + m(4) * (t110 * t85 + t111 * t86 + t147 * t46 + t148 * t47) + ((-qJD(1) * t173 - t254 * t268) * t253 + (-qJD(1) * t175 - t254 * t270) * t251 + t234 * (-qJD(1) * t155 - t254 * t269) + t235 * (-qJD(1) * t153 - t254 * t267) + (-t291 - t293) * qJD(2)) * t442 + ((qJD(1) * t174 - t252 * t268) * t253 + (qJD(1) * t176 - t252 * t270) * t251 + t234 * (qJD(1) * t156 - t252 * t269) + t235 * (qJD(1) * t154 - t252 * t267) + (-t292 - t294) * qJD(2)) * t441 + (t305 + t304) * qJD(2) * (t249 / 0.2e1 + t248 / 0.2e1); t252 * ((t252 * t89 + (t42 + t460) * qJD(1)) * t252 + (t43 * qJD(1) + (t153 * t359 + t155 * t360) * t254 + (-t90 + (-t400 - t402) * qJD(2) + (t152 - t294) * qJD(1)) * t252) * t254) - t254 * ((t254 * t90 + (t41 + t455) * qJD(1)) * t254 + (t40 * qJD(1) + (-t154 * t359 - t156 * t360 + t364) * t252 + (-t89 + (t401 + t403) * qJD(2) - t293 * qJD(1)) * t254) * t252) + t327 + (t25 * t84 + t26 * t83 + t5 * t8) * t443 + (t114 * t45 + t115 * t44 + t21 * t6) * t444 + t252 * ((t252 * t102 + (t54 + t461) * qJD(1)) * t252 + (t55 * qJD(1) + (t173 * t356 + t175 * t358) * t254 + (-t103 + (-t396 - t398) * qJD(2) + (t172 - t292) * qJD(1)) * t252) * t254) - t254 * ((t254 * t103 + (t53 + t456) * qJD(1)) * t254 + (t52 * qJD(1) + (-t174 * t356 - t176 * t358 + t363) * t252 + (-t102 + (t397 + t399) * qJD(2) - t291 * qJD(1)) * t254) * t252) + (t148 * t85 + t147 * t86 + (t252 * t163 + t164 * t254 + t382) * ((qJD(1) * t163 - t254 * t278 + t373) * t254 + (-t252 * t278 + (-t150 - t164 + t452) * qJD(1)) * t252 + t353)) * t445 + ((t252 * t177 + t178 * t254) * ((qJD(1) * t177 + t259) * t254 + (-t252 * t279 + (-t178 - t352 + t241) * qJD(1)) * t252) + t367 * t218 * t205) * t446 - t254 * t4 - t254 * t3 + ((-t40 - t52) * t254 + (t41 + t53) * t252) * t362 + ((-t42 - t54) * t254 + (t43 + t55) * t252) * t361; m(5) * (qJD(1) * t318 + t252 * t29 - t254 * t28) + m(6) * (qJD(1) * t320 - t23 * t254 + t252 * t24) + m(4) * (t252 * t47 - t254 * t46 + (t110 * t254 + t111 * t252) * qJD(1)); m(6) * (qJD(1) * t319 + t252 * t25 - t254 * t26) + m(5) * (qJD(1) * t451 + t252 * t44 - t254 * t45) + m(4) * (t252 * t85 - t254 * t86 + (t147 * t252 + t148 * t254) * qJD(1)); 0; t257 + m(6) * (-t108 * t23 - t109 * t24 + t30 * t62 + t31 * t63) + (-t252 * t28 - t254 * t29 + (t252 * t95 - t254 * t96) * qJD(1)) * t438 - m(5) * t318 * t166; m(6) * (-t108 * t26 - t109 * t25 + t27 * t5 + t30 * t84 + t31 * t83 + t7 * t8) + m(5) * (-t166 * t451 + t22 * t21 + t70 * t6) + (-t252 * t45 - t254 * t44 + (-t114 * t254 + t115 * t252) * qJD(1)) * t438 + t262; m(6) * (qJD(1) * t299 + t30 * t252 - t254 * t31); (t166 * t192 * t367 + t22 * t70) * t444 + (-t108 * t31 - t109 * t30 + t27 * t7) * t443 + t262; m(6) * (t320 * t389 + (t23 * t252 + t24 * t254 + (-t252 * t62 + t254 * t63) * qJD(1)) * t227); m(6) * ((t246 * t319 - t5) * t228 + (t246 * t8 + t25 * t254 + t252 * t26 + (-t252 * t84 + t254 * t83) * qJD(1)) * t227); 0; m(6) * ((t246 * t299 - t7) * t228 + (t246 * t27 + t252 * t31 + t254 * t30 + (-t108 * t254 + t109 * t252) * qJD(1)) * t227); (-0.1e1 + t367) * t227 * t389 * t443;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
