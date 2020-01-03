% Calculate time derivative of joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:53
% EndTime: 2019-12-31 20:55:14
% DurationCPUTime: 12.46s
% Computational Cost: add. (11914->625), mult. (12191->820), div. (0->0), fcn. (9492->8), ass. (0->355)
t275 = qJ(2) + qJ(3);
t260 = pkin(8) + t275;
t254 = cos(t260);
t253 = sin(t260);
t440 = Icges(6,5) * t253;
t210 = -Icges(6,3) * t254 + t440;
t446 = Icges(5,4) * t253;
t213 = Icges(5,2) * t254 + t446;
t514 = -t213 + t210;
t439 = Icges(6,5) * t254;
t214 = Icges(6,1) * t253 - t439;
t445 = Icges(5,4) * t254;
t215 = Icges(5,1) * t253 + t445;
t513 = t215 + t214;
t466 = rSges(6,1) + pkin(4);
t512 = rSges(6,3) + qJ(5);
t277 = sin(qJ(1));
t279 = cos(qJ(1));
t332 = Icges(6,3) * t253 + t439;
t152 = Icges(6,6) * t277 + t279 * t332;
t337 = -Icges(5,2) * t253 + t445;
t158 = Icges(5,6) * t277 + t279 * t337;
t511 = t152 - t158;
t341 = Icges(6,1) * t254 + t440;
t160 = Icges(6,4) * t277 + t279 * t341;
t342 = Icges(5,1) * t254 - t446;
t162 = Icges(5,5) * t277 + t279 * t342;
t510 = t160 + t162;
t272 = qJD(2) + qJD(3);
t509 = (-t332 + t337) * t272;
t508 = (t341 + t342) * t272;
t333 = Icges(5,5) * t254 - Icges(5,6) * t253;
t154 = Icges(5,3) * t277 + t279 * t333;
t336 = Icges(6,4) * t254 + Icges(6,6) * t253;
t156 = Icges(6,2) * t277 + t279 * t336;
t263 = sin(t275);
t264 = cos(t275);
t334 = Icges(4,5) * t264 - Icges(4,6) * t263;
t183 = Icges(4,3) * t277 + t279 * t334;
t507 = t154 + t156 + t183;
t211 = Icges(5,5) * t253 + Icges(5,6) * t254;
t212 = Icges(6,4) * t253 - Icges(6,6) * t254;
t223 = Icges(4,5) * t263 + Icges(4,6) * t264;
t495 = t211 + t212 + t223;
t448 = Icges(4,4) * t263;
t224 = Icges(4,2) * t264 + t448;
t447 = Icges(4,4) * t264;
t225 = Icges(4,1) * t263 + t447;
t506 = t224 * t263 - t225 * t264 - t253 * t514 - t513 * t254;
t151 = -Icges(6,6) * t279 + t277 * t332;
t157 = -Icges(5,6) * t279 + t277 * t337;
t505 = -t151 + t157;
t159 = -Icges(6,4) * t279 + t277 * t341;
t161 = -Icges(5,5) * t279 + t277 * t342;
t504 = t159 + t161;
t266 = t277 * rSges(4,3);
t457 = rSges(4,1) * t264;
t349 = -rSges(4,2) * t263 + t457;
t191 = t279 * t349 + t266;
t268 = t277 * rSges(6,2);
t416 = t254 * t279;
t418 = t253 * t279;
t502 = t466 * t416 + t512 * t418 + t268;
t217 = rSges(6,1) * t253 - rSges(6,3) * t254;
t501 = pkin(4) * t253 - qJ(5) * t254 + t217;
t425 = t210 * t272;
t88 = -qJD(1) * t151 - t279 * t425;
t424 = t213 * t272;
t94 = -qJD(1) * t157 - t279 * t424;
t500 = -t510 * t272 + t88 - t94;
t423 = t214 * t272;
t96 = -qJD(1) * t159 - t279 * t423;
t422 = t215 * t272;
t98 = -qJD(1) * t161 - t279 * t422;
t499 = t511 * t272 + t96 + t98;
t280 = -pkin(7) - pkin(6);
t276 = sin(qJ(2));
t388 = qJD(2) * t276;
t382 = pkin(2) * t388;
t498 = qJD(1) * t280 + t382;
t153 = -Icges(5,3) * t279 + t277 * t333;
t155 = -Icges(6,2) * t279 + t277 * t336;
t182 = -Icges(4,3) * t279 + t277 * t334;
t327 = t151 * t253 + t159 * t254;
t482 = t279 * t327;
t325 = t157 * t253 - t161 * t254;
t483 = t279 * t325;
t338 = -Icges(4,2) * t263 + t447;
t184 = -Icges(4,6) * t279 + t277 * t338;
t343 = Icges(4,1) * t264 - t448;
t186 = -Icges(4,5) * t279 + t277 * t343;
t323 = t184 * t263 - t186 * t264;
t484 = t279 * t323;
t497 = -t482 + t483 + t484 + (-t153 - t155 - t182) * t277;
t185 = Icges(4,6) * t277 + t279 * t338;
t187 = Icges(4,5) * t277 + t279 * t343;
t322 = t185 * t263 - t187 * t264;
t324 = t158 * t253 - t162 * t254;
t326 = t152 * t253 + t160 * t254;
t496 = (-t322 - t324 + t326) * t279 + t507 * t277;
t198 = t338 * t272;
t199 = t343 * t272;
t420 = t225 * t272;
t421 = t224 * t272;
t494 = (t199 - t421) * t264 + (-t198 - t420) * t263 + (t425 - t424 + t508) * t254 + (-t423 - t422 - t509) * t253 + t495 * qJD(1);
t492 = (t333 + t334 + t336) * t272 + t506 * qJD(1);
t218 = rSges(5,1) * t253 + rSges(5,2) * t254;
t310 = t218 * t272;
t278 = cos(qJ(2));
t244 = rSges(3,1) * t276 + rSges(3,2) * t278;
t308 = qJD(2) * t244;
t491 = t277 * t308;
t449 = Icges(3,4) * t278;
t340 = -Icges(3,2) * t276 + t449;
t204 = Icges(3,6) * t277 + t279 * t340;
t450 = Icges(3,4) * t276;
t345 = Icges(3,1) * t278 - t450;
t206 = Icges(3,5) * t277 + t279 * t345;
t320 = t204 * t276 - t206 * t278;
t490 = t277 * t320;
t489 = t277 * t322;
t488 = t277 * t324;
t487 = t277 * t326;
t255 = t278 * pkin(2) + pkin(1);
t460 = pkin(1) - t255;
t486 = t277 * t460;
t203 = -Icges(3,6) * t279 + t277 * t340;
t205 = -Icges(3,5) * t279 + t277 * t345;
t321 = t203 * t276 - t205 * t278;
t485 = t279 * t321;
t347 = rSges(6,1) * t254 + rSges(6,3) * t253;
t367 = pkin(4) * t272 - qJD(5);
t419 = t253 * t272;
t481 = -qJ(5) * t419 - t254 * t367 - t347 * t272;
t478 = qJD(1) * t153;
t477 = qJD(1) * t155;
t476 = qJD(1) * t182;
t335 = Icges(3,5) * t278 - Icges(3,6) * t276;
t201 = -Icges(3,3) * t279 + t277 * t335;
t271 = -qJ(4) + t280;
t397 = t271 - t280;
t231 = pkin(3) * t264 + t255;
t403 = t231 - t255;
t132 = t277 * t403 + t279 * t397;
t472 = 2 * m(3);
t471 = 2 * m(4);
t470 = 2 * m(5);
t469 = 2 * m(6);
t273 = t277 ^ 2;
t274 = t279 ^ 2;
t468 = t277 / 0.2e1;
t467 = -t279 / 0.2e1;
t465 = m(3) * t244;
t453 = rSges(4,2) * t264;
t226 = rSges(4,1) * t263 + t453;
t464 = m(4) * t226;
t463 = pkin(2) * t276;
t462 = pkin(3) * t263;
t461 = t277 * pkin(6);
t270 = t279 * pkin(6);
t459 = -pkin(6) - t280;
t458 = rSges(3,1) * t278;
t456 = rSges(5,1) * t254;
t455 = rSges(3,2) * t276;
t452 = rSges(3,3) * t279;
t267 = t277 * rSges(3,3);
t265 = t277 * rSges(5,3);
t200 = t349 * t272;
t430 = t200 * t277;
t429 = t203 * t278;
t428 = t204 * t278;
t427 = t205 * t276;
t426 = t206 * t276;
t417 = t254 * t272;
t415 = t263 * t272;
t414 = t264 * t272;
t413 = t272 * t277;
t412 = t272 * t279;
t411 = t279 * t271;
t220 = t279 * t231;
t245 = t279 * t255;
t133 = -t277 * t397 + t220 - t245;
t410 = t277 * t132 + t279 * t133;
t169 = rSges(5,1) * t416 - rSges(5,2) * t418 + t265;
t409 = -t133 - t169;
t174 = t279 * t280 + t270 - t486;
t175 = -t279 * pkin(1) + t277 * t459 + t245;
t408 = t277 * t174 + t279 * t175;
t190 = -t279 * rSges(4,3) + t277 * t349;
t106 = t277 * t190 + t279 * t191;
t391 = qJD(1) * t277;
t373 = t263 * t391;
t240 = pkin(3) * t373;
t407 = t218 * t391 + t240;
t219 = -pkin(3) * t415 - t382;
t406 = qJD(4) * t277 + t279 * t219;
t379 = t254 * t412;
t386 = qJD(5) * t253;
t405 = qJ(5) * t379 + t279 * t386;
t390 = qJD(1) * t279;
t404 = rSges(6,2) * t390 + rSges(6,3) * t379;
t374 = t253 * t391;
t402 = rSges(5,2) * t374 + rSges(5,3) * t390;
t401 = rSges(4,2) * t373 + rSges(4,3) * t390;
t400 = qJD(4) * t279 + t271 * t391;
t399 = t498 * t277;
t398 = t279 * t458 + t267;
t396 = t273 + t274;
t395 = qJD(1) * t154;
t394 = qJD(1) * t156;
t393 = qJD(1) * t183;
t202 = Icges(3,3) * t277 + t279 * t335;
t392 = qJD(1) * t202;
t387 = qJD(2) * t278;
t385 = pkin(3) * t414;
t353 = t279 * t382;
t384 = t132 * t390 + t277 * (t219 * t277 + t390 * t403 + t399 - t400) + t279 * (-qJD(1) * t132 + t353 + t406);
t383 = t279 * t455;
t381 = pkin(2) * t387;
t380 = t253 * t412;
t309 = t226 * t272;
t378 = t277 * (t191 * qJD(1) - t277 * t309) + t279 * (-t412 * t453 + (-t263 * t412 - t264 * t391) * rSges(4,1) + t401) + t190 * t390;
t377 = t277 * ((-t279 * t460 - t461) * qJD(1) - t399) + t279 * (-t353 + (t279 * t459 + t486) * qJD(1)) + t174 * t390;
t376 = -t133 - t502;
t375 = t501 * t391 + t240;
t372 = t276 * t391;
t369 = -t226 - t463;
t368 = -t218 - t462;
t97 = qJD(1) * t160 - t277 * t423;
t366 = t151 * t272 + t97;
t99 = qJD(1) * t162 - t277 * t422;
t364 = t157 * t272 - t99;
t89 = qJD(1) * t152 - t277 * t425;
t362 = t159 * t272 - t89;
t95 = qJD(1) * t158 - t277 * t424;
t360 = t161 * t272 + t95;
t109 = -qJD(1) * t184 - t279 * t421;
t358 = t187 * t272 + t109;
t110 = qJD(1) * t185 - t277 * t421;
t357 = t186 * t272 + t110;
t111 = -qJD(1) * t186 - t279 * t420;
t356 = -t185 * t272 + t111;
t112 = qJD(1) * t187 - t277 * t420;
t355 = t184 * t272 - t112;
t354 = -t277 * t271 + t220;
t348 = -rSges(5,2) * t253 + t456;
t167 = -rSges(5,3) * t279 + t277 * t348;
t38 = t277 * t167 + t279 * t169 + t410;
t352 = -t462 - t501;
t351 = -t462 - t463;
t350 = -t455 + t458;
t289 = -t253 * t512 - t254 * t466 - t231;
t286 = t289 * t277;
t76 = (rSges(6,2) - t271) * t279 + t286;
t77 = t354 + t502;
t346 = t277 * t77 + t279 * t76;
t344 = Icges(3,1) * t276 + t449;
t339 = Icges(3,2) * t278 + t450;
t307 = t351 - t501;
t102 = t307 * t277;
t103 = t307 * t279;
t331 = t102 * t277 + t103 * t279;
t115 = t352 * t277;
t116 = t352 * t279;
t330 = t115 * t277 + t116 * t279;
t316 = -t385 + t481;
t292 = -t254 * t391 - t380;
t315 = t167 * t390 + t277 * (-t277 * t310 + (t279 * t348 + t265) * qJD(1)) + t279 * (rSges(5,1) * t292 - rSges(5,2) * t379 + t402) + t384;
t314 = -pkin(1) - t350;
t313 = -t218 + t351;
t312 = -t255 - t349;
t311 = -t231 - t348;
t166 = -rSges(6,2) * t279 + t277 * t347;
t192 = (pkin(4) * t254 + qJ(5) * t253) * t277;
t30 = t410 + t502 * t279 + (t166 + t192) * t277;
t301 = t272 * t212;
t300 = t272 * t223;
t299 = t272 * t211;
t297 = -t381 - t385;
t296 = qJD(2) * t344;
t295 = qJD(2) * t339;
t294 = qJD(2) * (-Icges(3,5) * t276 - Icges(3,6) * t278);
t135 = t313 * t279;
t107 = -t279 * t300 - t476;
t108 = -t277 * t300 + t393;
t44 = -t155 * t279 + t277 * t327;
t45 = -t156 * t279 + t487;
t46 = -t153 * t279 - t277 * t325;
t47 = -t154 * t279 - t488;
t52 = -t182 * t279 - t277 * t323;
t53 = -t183 * t279 - t489;
t90 = -t279 * t299 - t478;
t91 = -t277 * t299 + t395;
t92 = -t279 * t301 - t477;
t93 = -t277 * t301 + t394;
t293 = ((-t44 - t46 - t52) * t391 + t497 * t390) * t279 + (((t92 + t90 + t107) * t277 + (-t487 + t488 + t489 - t497) * qJD(1)) * t277 + (t45 + t47 + t53) * t391 + t496 * t390 + ((-t358 * t263 + t356 * t264 - t108 - t91 - t93) * t277 + (t110 * t263 - t112 * t264 + t184 * t414 + t186 * t415 + t505 * t417 + t504 * t419 - t476 - t477 - t478) * t279 + (t499 * t277 + (-t97 - t99) * t279) * t254 + (t500 * t277 + (-t89 + t95) * t279) * t253 + ((t327 - t325 - t323 + t507) * t277 + t496) * qJD(1)) * t279) * t277;
t189 = t348 * t272;
t291 = -t189 + t297;
t290 = t166 * t390 + t192 * t390 + t277 * ((pkin(4) * t390 + qJ(5) * t413) * t254 + (qJ(5) * t390 - t277 * t367) * t253) + t279 * (pkin(4) * t292 - qJ(5) * t374 + t405) + t277 * (-t217 * t413 + (t279 * t347 + t268) * qJD(1)) + t279 * (rSges(6,1) * t292 - rSges(6,3) * t374 + t404) + t384;
t288 = t297 + t481;
t287 = rSges(3,2) * t372 + rSges(3,3) * t390 - t279 * t308;
t3 = (t279 * t93 + (t45 - t482) * qJD(1)) * t279 + (t44 * qJD(1) + (t152 * t417 - t160 * t419 + t253 * t88 + t254 * t96 + t394) * t277 + (-t92 - t366 * t254 + t362 * t253 + (-t155 + t326) * qJD(1)) * t279) * t277;
t4 = (t279 * t91 + (t47 + t483) * qJD(1)) * t279 + (t46 * qJD(1) + (-t158 * t417 - t162 * t419 - t253 * t94 + t254 * t98 + t395) * t277 + (-t90 + t364 * t254 + t360 * t253 + (-t153 - t324) * qJD(1)) * t279) * t277;
t6 = (t279 * t108 + (t53 + t484) * qJD(1)) * t279 + (t52 * qJD(1) + (-t109 * t263 + t111 * t264 - t185 * t414 - t187 * t415 + t393) * t277 + (-t107 + t355 * t264 + t357 * t263 + (-t182 - t322) * qJD(1)) * t279) * t277;
t285 = (-t4 - t3 - t6) * t279 + t293;
t281 = (t499 * t253 - t500 * t254 + t263 * t356 + t264 * t358 + t492 * t277 + t494 * t279) * t468 + (-t263 * t355 + t264 * t357 + (t360 + t362) * t254 + (-t364 + t366) * t253 - t492 * t279 + t494 * t277) * t467 + (t184 * t264 + t186 * t263 + t504 * t253 + t505 * t254 - t277 * t506 - t495 * t279) * t391 / 0.2e1 + (t185 * t264 + t187 * t263 + t510 * t253 - t511 * t254 + t495 * t277 - t279 * t506) * t390 / 0.2e1;
t250 = pkin(2) * t372;
t230 = t350 * qJD(2);
t209 = -t383 + t398;
t208 = t277 * t350 - t452;
t173 = t369 * t279;
t172 = t369 * t277;
t148 = t461 + (pkin(1) - t455) * t279 + t398;
t147 = t277 * t314 + t270 + t452;
t144 = t368 * t279;
t143 = t368 * t277;
t134 = t313 * t277;
t131 = -t277 * t280 + t191 + t245;
t130 = (rSges(4,3) - t280) * t279 + t312 * t277;
t125 = t277 * t294 + t392;
t124 = -qJD(1) * t201 + t279 * t294;
t118 = t169 + t354;
t117 = (rSges(5,3) - t271) * t279 + t311 * t277;
t114 = t491 + ((-rSges(3,3) - pkin(6)) * t277 + t314 * t279) * qJD(1);
t113 = (t270 + (-pkin(1) - t458) * t277) * qJD(1) + t287;
t101 = -t226 * t390 - t430 + (-t276 * t390 - t277 * t387) * pkin(2);
t100 = t226 * t391 + t250 + (-t200 - t381) * t279;
t75 = -t218 * t390 - t189 * t277 + (-t263 * t390 - t264 * t413) * pkin(3);
t74 = (-t189 - t385) * t279 + t407;
t69 = t202 * t277 - t320 * t279;
t68 = t201 * t277 - t485;
t67 = -t202 * t279 - t490;
t66 = -t201 * t279 - t277 * t321;
t65 = t226 * t413 + (t279 * t312 - t266) * qJD(1) + t399;
t64 = (-t255 - t457) * t391 + (-t309 - t498) * t279 + t401;
t57 = qJD(1) * t135 + t277 * t291;
t56 = t279 * t291 + t250 + t407;
t43 = (-t219 + t310) * t277 + (t279 * t311 - t265) * qJD(1) + t400;
t42 = -t279 * t310 + (-t411 + (-t231 - t456) * t277) * qJD(1) + t402 + t406;
t41 = t106 + t408;
t40 = qJD(1) * t116 + t277 * t316;
t39 = t279 * t316 + t375;
t37 = qJD(1) * t103 + t277 * t288;
t36 = t279 * t288 + t250 + t375;
t35 = -t191 * t391 + t378;
t34 = t289 * t390 + (-rSges(6,2) * qJD(1) - t386 - t219 + (t253 * t466 - t254 * t512) * t272) * t277 + t400;
t33 = -t466 * t380 + (t286 - t411) * qJD(1) + t404 + t405 + t406;
t27 = t38 + t408;
t12 = t30 + t408;
t11 = (-t175 - t191) * t391 + t377 + t378;
t10 = t391 * t409 + t315;
t9 = (-t175 + t409) * t391 + t315 + t377;
t8 = t376 * t391 + t290;
t7 = (-t175 + t376) * t391 + t290 + t377;
t1 = [(t117 * t43 + t118 * t42) * t470 + (t33 * t77 + t34 * t76) * t469 + t225 * t414 + t263 * t199 - t224 * t415 + t264 * t198 + (t130 * t65 + t131 * t64) * t471 + (t113 * t148 + t114 * t147) * t472 + t514 * t419 + t513 * t417 + (t345 - t339) * t388 + (t344 + t340) * t387 + t509 * t254 + t508 * t253; m(6) * (t102 * t33 + t103 * t34 + t36 * t76 + t37 * t77) + m(5) * (t117 * t56 + t118 * t57 + t134 * t42 + t135 * t43) + m(4) * (t100 * t130 + t101 * t131 + t172 * t64 + t173 * t65) + t281 + m(3) * ((-t113 * t277 - t114 * t279) * t244 + (-t147 * t279 - t148 * t277) * t230) + ((t428 / 0.2e1 + t426 / 0.2e1 - t148 * t465) * t279 + (t147 * t465 + t429 / 0.2e1 + t427 / 0.2e1) * t277) * qJD(1) + (t274 / 0.2e1 + t273 / 0.2e1) * t335 * qJD(2) + (-qJD(2) * t320 + (-qJD(1) * t203 - t279 * t295) * t278 + (-qJD(1) * t205 - t279 * t296) * t276) * t468 + (-qJD(2) * t321 + (qJD(1) * t204 - t277 * t295) * t278 + (qJD(1) * t206 - t277 * t296) * t276) * t467; (t102 * t37 + t103 * t36 + t12 * t7) * t469 + (t134 * t57 + t135 * t56 + t27 * t9) * t470 + (t100 * t173 + t101 * t172 + t11 * t41) * t471 + t277 * ((t277 * t124 + (t68 + t490) * qJD(1)) * t277 + (t69 * qJD(1) + (t203 * t387 + t205 * t388) * t279 + (-t125 + (-t426 - t428) * qJD(2) + (t202 - t321) * qJD(1)) * t277) * t279) - t279 * t4 - t279 * t3 - t279 * t6 - t279 * ((t279 * t125 + (t67 + t485) * qJD(1)) * t279 + (t66 * qJD(1) + (-t204 * t387 - t206 * t388 + t392) * t277 + (-t124 + (t427 + t429) * qJD(2) - t320 * qJD(1)) * t279) * t277) + (t277 * t69 - t279 * t68) * t390 + (t277 * t67 - t279 * t66) * t391 + ((t208 * t277 + t209 * t279) * ((qJD(1) * t208 + t287) * t279 + (-t491 + (-t209 - t383 + t267) * qJD(1)) * t277) + t396 * t244 * t230) * t472 + t293; t281 + m(6) * (t115 * t33 + t116 * t34 + t39 * t76 + t40 * t77) + m(5) * (t117 * t74 + t118 * t75 + t143 * t42 + t144 * t43) + (-t277 * t64 - t279 * t65 + (t130 * t277 - t131 * t279) * qJD(1)) * t464 + m(4) * (-t130 * t279 - t131 * t277) * t200; (-t100 * t279 - t101 * t277 + (-t172 * t279 + t173 * t277) * qJD(1)) * t464 + m(4) * (-t173 * t200 * t279 + t106 * t11 - t172 * t430 + t35 * t41) + m(6) * (t102 * t40 + t103 * t39 + t115 * t37 + t116 * t36 + t12 * t8 + t30 * t7) + m(5) * (t10 * t27 + t134 * t75 + t135 * t74 + t143 * t57 + t144 * t56 + t38 * t9) + t285; (t10 * t38 + t143 * t75 + t144 * t74) * t470 + (t115 * t40 + t116 * t39 + t30 * t8) * t469 + (t200 * t226 * t396 + t106 * t35) * t471 + t285; m(5) * (t277 * t43 - t279 * t42 + (t117 * t279 + t118 * t277) * qJD(1)) + m(6) * (qJD(1) * t346 + t277 * t34 - t279 * t33); m(6) * (qJD(1) * t331 + t277 * t36 - t279 * t37) + m(5) * (t277 * t56 - t279 * t57 + (t134 * t277 + t135 * t279) * qJD(1)); m(5) * (t277 * t74 - t279 * t75 + (t143 * t277 + t144 * t279) * qJD(1)) + m(6) * (qJD(1) * t330 + t277 * t39 - t279 * t40); 0; m(6) * (t346 * t417 + (t277 * t33 + t279 * t34 + (-t277 * t76 + t279 * t77) * qJD(1)) * t253); m(6) * ((t272 * t331 - t7) * t254 + (t12 * t272 + t277 * t37 + t279 * t36 + (t102 * t279 - t103 * t277) * qJD(1)) * t253); m(6) * ((t272 * t330 - t8) * t254 + (t272 * t30 + t277 * t40 + t279 * t39 + (t115 * t279 - t116 * t277) * qJD(1)) * t253); 0; (-0.1e1 + t396) * t253 * t417 * t469;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
