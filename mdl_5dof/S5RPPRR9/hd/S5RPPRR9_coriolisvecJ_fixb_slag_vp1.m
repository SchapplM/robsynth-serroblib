% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR9_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:44
% DurationCPUTime: 20.07s
% Computational Cost: add. (13655->703), mult. (32711->964), div. (0->0), fcn. (37025->8), ass. (0->330)
t269 = sin(qJ(4));
t271 = cos(qJ(4));
t268 = sin(qJ(5));
t270 = cos(qJ(5));
t444 = sin(pkin(8));
t445 = cos(pkin(8));
t461 = sin(qJ(1));
t462 = cos(qJ(1));
t225 = -t444 * t461 - t445 * t462;
t226 = t462 * t444 - t461 * t445;
t422 = t268 * t271;
t161 = -t225 * t270 + t226 * t422;
t424 = t226 * t269;
t420 = t270 * t271;
t160 = t225 * t268 + t226 * t420;
t441 = Icges(6,4) * t160;
t91 = -Icges(6,2) * t161 + Icges(6,6) * t424 + t441;
t154 = Icges(6,4) * t161;
t94 = Icges(6,1) * t160 + Icges(6,5) * t424 - t154;
t511 = t268 * t91 - t270 * t94;
t88 = Icges(6,5) * t160 - Icges(6,6) * t161 + Icges(6,3) * t424;
t31 = -t511 * t269 - t271 * t88;
t393 = qJD(5) * t269;
t400 = qJD(4) * t225;
t179 = t226 * t393 + t400;
t399 = qJD(4) * t226;
t180 = -t225 * t393 + t399;
t426 = t225 * t269;
t163 = t225 * t422 + t226 * t270;
t164 = -t225 * t420 + t226 * t268;
t512 = t163 * t91 + t164 * t94;
t25 = t426 * t88 - t512;
t392 = qJD(5) * t271;
t251 = qJD(1) + t392;
t90 = Icges(6,5) * t164 + Icges(6,6) * t163 - Icges(6,3) * t426;
t440 = Icges(6,4) * t164;
t93 = Icges(6,2) * t163 - Icges(6,6) * t426 + t440;
t155 = Icges(6,4) * t163;
t96 = Icges(6,1) * t164 - Icges(6,5) * t426 + t155;
t26 = t163 * t93 + t164 * t96 - t90 * t426;
t439 = Icges(6,4) * t268;
t250 = t269 * t439;
t421 = t269 * t270;
t436 = Icges(6,5) * t271;
t197 = -Icges(6,1) * t421 + t250 + t436;
t330 = Icges(6,5) * t270 - Icges(6,6) * t268;
t285 = -Icges(6,3) * t271 + t269 * t330;
t438 = Icges(6,4) * t270;
t333 = -Icges(6,2) * t268 + t438;
t286 = -Icges(6,6) * t271 + t269 * t333;
t57 = -t163 * t286 + t164 * t197 + t285 * t426;
t11 = -t179 * t25 + t180 * t26 + t57 * t251;
t332 = Icges(5,5) * t271 - Icges(5,6) * t269;
t136 = Icges(5,3) * t225 + t226 * t332;
t442 = Icges(5,4) * t271;
t335 = -Icges(5,2) * t269 + t442;
t139 = Icges(5,6) * t225 + t226 * t335;
t443 = Icges(5,4) * t269;
t338 = Icges(5,1) * t271 - t443;
t142 = Icges(5,5) * t225 + t226 * t338;
t501 = t139 * t269 - t142 * t271;
t51 = t225 * t136 - t501 * t226;
t23 = t160 * t94 - t161 * t91 + t424 * t88;
t259 = qJD(2) * t461;
t346 = rSges(6,1) * t270 - rSges(6,2) * t268;
t288 = -rSges(6,3) * t271 + t269 * t346;
t459 = t269 * pkin(4);
t354 = -pkin(7) * t271 + t459;
t397 = qJD(4) * t354;
t347 = rSges(6,1) * t160 - t161 * rSges(6,2);
t98 = -rSges(6,3) * t424 - t347;
t508 = t179 * t288 + t225 * t397 - t251 * t98 + t259;
t491 = t160 * t197 + t161 * t286 - t285 * t424;
t507 = t179 * t23 + t491 * t251;
t138 = Icges(5,3) * t226 - t225 * t332;
t504 = t226 * t138;
t331 = Icges(5,5) * t269 + Icges(5,6) * t271;
t166 = t331 * t225;
t334 = Icges(5,2) * t271 + t443;
t337 = Icges(5,1) * t269 + t442;
t323 = -t269 * t334 + t271 * t337;
t109 = t226 * t323 + t166;
t503 = qJD(1) * t109;
t72 = -t139 * t271 - t142 * t269;
t469 = -t179 / 0.2e1;
t355 = pkin(4) * t271 + pkin(7) * t269;
t479 = t355 * qJD(4);
t500 = t225 * t479;
t499 = t226 * t479;
t165 = t331 * t226;
t214 = t225 * qJD(1);
t215 = t226 * qJD(1);
t498 = t215 * pkin(3) + t214 * pkin(6);
t266 = t462 * pkin(2);
t401 = t462 * pkin(1) + t461 * qJ(2);
t485 = t266 + t401;
t361 = -t225 * pkin(3) + t485;
t260 = qJD(2) * t462;
t497 = -t485 * qJD(1) + t260;
t395 = qJD(4) * t269;
t309 = t215 * t271 + t225 * t395;
t394 = qJD(4) * t271;
t429 = t215 * t269;
t310 = t225 * t394 - t429;
t87 = rSges(5,1) * t309 + rSges(5,2) * t310 + t214 * rSges(5,3);
t496 = -t225 * rSges(4,1) - t226 * rSges(4,2) + t485;
t492 = t160 * t96 - t161 * t93;
t490 = 0.2e1 * qJD(4);
t487 = t355 * t226;
t262 = t462 * qJ(2);
t388 = t461 * pkin(1);
t241 = t388 - t262;
t387 = t461 * pkin(2);
t294 = t226 * rSges(4,1) - t225 * rSges(4,2) - t387;
t486 = -t241 + t294;
t304 = -t388 - t387;
t484 = t304 + t387;
t302 = t462 * rSges(3,1) + t461 * rSges(3,3);
t483 = t401 + t302;
t233 = qJD(1) * t241;
t403 = t259 - t233;
t373 = qJD(1) * t462;
t404 = qJ(2) * t373 + t259;
t482 = t404 - t403;
t222 = t225 * pkin(6);
t481 = t226 * pkin(3) + t222;
t349 = rSges(5,1) * t271 - rSges(5,2) * t269;
t453 = t225 * rSges(5,3);
t145 = t226 * t349 + t453;
t210 = t215 * pkin(6);
t478 = -t210 + t497;
t477 = -qJD(1) * t481 + t404 + t498;
t415 = t334 * t226 - t142;
t417 = -t337 * t226 - t139;
t476 = t269 * t415 + t271 * t417;
t194 = -Icges(6,3) * t269 - t271 * t330;
t324 = -t197 * t270 - t268 * t286;
t339 = t268 * t93 - t270 * t96;
t475 = -t180 * (t285 * t225 + t339) + t179 * (t285 * t226 - t511) - t251 * (t194 + t324);
t218 = (Icges(6,1) * t268 + t438) * t269;
t278 = t180 * (Icges(6,1) * t163 - t440 - t93) - t179 * (Icges(6,1) * t161 + t441 + t91) - t251 * (-t286 - t218);
t272 = qJD(1) ^ 2;
t375 = t225 * t392;
t130 = t215 * t393 + (t214 - t375) * qJD(4);
t474 = t130 / 0.2e1;
t374 = t226 * t392;
t131 = -t214 * t393 + (t215 - t374) * qJD(4);
t473 = t131 / 0.2e1;
t472 = -t180 / 0.2e1;
t471 = t180 / 0.2e1;
t470 = t179 / 0.2e1;
t468 = t214 / 0.2e1;
t467 = t215 / 0.2e1;
t464 = -t251 / 0.2e1;
t463 = t251 / 0.2e1;
t430 = t214 * t269;
t312 = t226 * t394 + t430;
t311 = -t214 * t271 + t226 * t395;
t291 = -qJD(5) * t225 + t311;
t351 = t215 + t374;
t76 = -t268 * t291 + t270 * t351;
t77 = t268 * t351 + t270 * t291;
t39 = Icges(6,5) * t77 + Icges(6,6) * t76 - Icges(6,3) * t312;
t41 = Icges(6,4) * t77 + Icges(6,2) * t76 - Icges(6,6) * t312;
t43 = Icges(6,1) * t77 + Icges(6,4) * t76 - Icges(6,5) * t312;
t7 = (-qJD(4) * t511 + t39) * t271 + (qJD(4) * t88 + t268 * t41 - t270 * t43 + (-t268 * t94 - t270 * t91) * qJD(5)) * t269;
t458 = t7 * t179;
t290 = qJD(5) * t226 + t309;
t352 = t214 + t375;
t78 = -t268 * t290 + t270 * t352;
t79 = t268 * t352 + t270 * t290;
t40 = Icges(6,5) * t79 + Icges(6,6) * t78 - Icges(6,3) * t310;
t42 = Icges(6,4) * t79 + Icges(6,2) * t78 - Icges(6,6) * t310;
t44 = Icges(6,1) * t79 + Icges(6,4) * t78 - Icges(6,5) * t310;
t8 = (qJD(4) * t339 + t40) * t271 + (-qJD(4) * t90 + t268 * t42 - t270 * t44 + (t268 * t96 + t270 * t93) * qJD(5)) * t269;
t457 = t8 * t180;
t456 = rSges(6,3) * t269;
t454 = t215 * rSges(5,3);
t452 = t31 * t131;
t32 = t269 * t339 + t271 * t90;
t451 = t32 * t130;
t102 = pkin(4) * t309 - pkin(7) * t310;
t46 = t79 * rSges(6,1) + t78 * rSges(6,2) - rSges(6,3) * t310;
t450 = t102 + t46;
t447 = -t487 + t98;
t425 = t225 * t271;
t178 = -pkin(4) * t425 - pkin(7) * t426;
t99 = t164 * rSges(6,1) + t163 * rSges(6,2) - rSges(6,3) * t426;
t446 = t178 + t99;
t141 = Icges(5,6) * t226 - t225 * t335;
t431 = t141 * t269;
t423 = t226 * t271;
t144 = Icges(5,5) * t226 - t225 * t338;
t419 = -t225 * t138 - t144 * t423;
t418 = -t144 * t425 + t504;
t416 = -t337 * t225 + t141;
t414 = t334 * t225 + t144;
t200 = -t271 * t346 - t456;
t224 = (rSges(6,1) * t268 + rSges(6,2) * t270) * t269;
t153 = qJD(4) * t200 + qJD(5) * t224;
t413 = -t153 + t479;
t213 = qJD(1) * t401 - t260;
t412 = t214 * pkin(3) - t210 - t213;
t410 = -t288 - t354;
t372 = qJD(1) * t461;
t409 = (-pkin(1) * t372 + t259 + t404) * qJD(1);
t408 = t215 * rSges(4,1) - t214 * rSges(4,2);
t406 = -t334 + t338;
t405 = -t335 - t337;
t348 = rSges(5,1) * t269 + rSges(5,2) * t271;
t398 = qJD(4) * t348;
t391 = t332 * qJD(1);
t24 = -t90 * t424 - t492;
t390 = -t462 / 0.2e1;
t389 = t461 / 0.2e1;
t384 = t461 * rSges(3,1);
t371 = qJD(4) * t393;
t367 = t400 / 0.2e1;
t366 = -t399 / 0.2e1;
t364 = -t395 / 0.2e1;
t363 = -t394 / 0.2e1;
t360 = qJD(5) * t364;
t356 = t77 * rSges(6,1) + t76 * rSges(6,2);
t350 = t214 * rSges(4,1) + t215 * rSges(4,2);
t345 = -t225 * t24 - t226 * t23;
t344 = -t225 * t26 - t226 * t25;
t343 = -t225 * t32 - t226 * t31;
t321 = -t241 - t387 + t481;
t61 = t225 * t398 + t259 + (t145 + t321) * qJD(1);
t147 = -rSges(5,1) * t425 + rSges(5,2) * t426 + t226 * rSges(5,3);
t320 = pkin(6) * t226 + t361;
t62 = t226 * t398 - t260 + (t147 + t320) * qJD(1);
t342 = -t225 * t61 - t226 * t62;
t341 = -t225 * t98 + t226 * t99;
t336 = Icges(6,1) * t270 - t439;
t326 = -t144 * t271 + t431;
t325 = -t145 * t226 + t147 * t225;
t256 = qJD(1) * t260;
t322 = -t266 * t272 + t256;
t319 = pkin(3) + t349;
t53 = -t136 * t226 - t139 * t426 + t142 * t425;
t113 = t269 * t324 - t271 * t285;
t216 = (Icges(6,5) * t268 + Icges(6,6) * t270) * t269;
t148 = qJD(4) * t194 + qJD(5) * t216;
t196 = -Icges(6,6) * t269 - t271 * t333;
t149 = (Icges(6,2) * t270 + t439) * t393 + t196 * qJD(4);
t198 = -Icges(6,5) * t269 - t271 * t336;
t150 = qJD(4) * t198 + qJD(5) * t218;
t38 = (qJD(4) * t324 + t148) * t271 + (qJD(4) * t285 + t149 * t268 - t150 * t270 + (t197 * t268 - t270 * t286) * qJD(5)) * t269;
t318 = -t113 * t371 + t38 * t251;
t315 = -t269 * t39 + t394 * t88;
t314 = -t269 * t40 - t394 * t90;
t308 = -t272 * t387 + t409;
t52 = t141 * t424 + t419;
t307 = (-t225 * t51 + t226 * t52) * qJD(4);
t54 = t141 * t426 + t418;
t306 = (-t225 * t53 + t226 * t54) * qJD(4);
t303 = -t384 - t388;
t299 = qJD(1) * t498 + t308;
t298 = -t179 * t88 - t180 * t90 + t251 * t285;
t297 = -(Icges(6,5) * t161 + Icges(6,6) * t160) * t179 + (Icges(6,5) * t163 - Icges(6,6) * t164) * t180 + t216 * t251;
t296 = pkin(3) + t355 + t456;
t295 = t269 * t414 + t271 * t416;
t292 = t222 + t262 + t304;
t289 = t269 * t297;
t287 = t269 * t336 - t436;
t284 = (t269 * t405 + t271 * t406) * qJD(1);
t86 = rSges(5,1) * t311 + rSges(5,2) * t312 + t454;
t283 = -t145 * t214 - t147 * t215 + t225 * t87 + t226 * t86;
t279 = (-Icges(6,2) * t164 + t155 + t96) * t180 - (Icges(6,2) * t160 + t154 - t94) * t179 + (Icges(6,2) * t421 + t197 + t250) * t251;
t228 = t335 * qJD(4);
t229 = t338 * qJD(4);
t275 = -t228 * t269 + t229 * t271 + (-t269 * t337 - t271 * t334) * qJD(4);
t28 = t180 * t98 + t179 * t99 - qJD(3) + (t178 * t225 - t226 * t487) * qJD(4);
t36 = (t487 + t321) * qJD(1) + t508;
t37 = t226 * t397 + t180 * t288 + t251 * t99 - t260 + (t178 + t320) * qJD(1);
t274 = t28 * t341 - (t225 * t37 - t226 * t36) * t288;
t273 = t475 * t269;
t264 = t462 * rSges(3,3);
t258 = rSges(3,3) * t373;
t230 = t349 * qJD(4);
t227 = t332 * qJD(4);
t177 = t354 * t225;
t175 = t354 * t226;
t172 = t348 * t225;
t171 = t348 * t226;
t152 = -qJD(1) * t213 - t272 * t302 + t256;
t151 = qJD(1) * (-rSges(3,1) * t372 + t258) + t409;
t134 = qJD(1) * t486 + t259;
t133 = t288 * t225;
t132 = t288 * t226;
t129 = t287 * t225;
t128 = t287 * t226;
t127 = t286 * t225;
t126 = t286 * t226;
t115 = (-t213 + t350) * qJD(1) + t322;
t114 = qJD(1) * t408 + t308;
t112 = rSges(6,1) * t163 - rSges(6,2) * t164;
t111 = rSges(6,1) * t161 + rSges(6,2) * t160;
t110 = t225 * t323 - t165;
t101 = pkin(4) * t311 - pkin(7) * t312;
t100 = t110 * qJD(1);
t81 = Icges(5,5) * t309 + Icges(5,6) * t310 + Icges(5,3) * t214;
t80 = Icges(5,5) * t311 + Icges(5,6) * t312 + Icges(5,3) * t215;
t73 = t141 * t271 + t144 * t269;
t60 = qJD(4) * t325 - qJD(3);
t50 = (-t215 * t348 + t225 * t230) * qJD(4) + (-t86 + t412) * qJD(1) + t322;
t49 = qJD(1) * t87 + (t214 * t348 + t226 * t230) * qJD(4) + t299;
t48 = -t214 * t331 - t215 * t323 + t225 * t275 - t226 * t227;
t47 = t214 * t323 - t215 * t331 + t225 * t227 + t226 * t275;
t45 = -rSges(6,3) * t312 + t356;
t30 = qJD(4) * t326 - t269 * (Icges(5,1) * t309 + Icges(5,4) * t310 + Icges(5,5) * t214) - t271 * (Icges(5,4) * t309 + Icges(5,2) * t310 + Icges(5,6) * t214);
t29 = -qJD(4) * t501 - t269 * (Icges(5,1) * t311 + Icges(5,4) * t312 + Icges(5,5) * t215) - t271 * (Icges(5,4) * t311 + Icges(5,2) * t312 + Icges(5,6) * t215);
t27 = t283 * qJD(4);
t22 = t100 + t306;
t21 = t307 + t503;
t20 = -t148 * t426 + t149 * t163 + t150 * t164 + t197 * t79 + t285 * t310 - t286 * t78;
t19 = -t148 * t424 + t149 * t161 - t150 * t160 + t197 * t77 + t285 * t312 - t286 * t76;
t14 = -t131 * t288 - t179 * t153 - t251 * t45 + (-t215 * t354 + t393 * t98 + t500) * qJD(4) + (-t101 + t412) * qJD(1) + t322;
t13 = qJD(1) * t102 + t130 * t288 - t180 * t153 + t251 * t46 + (t214 * t354 - t393 * t99 + t499) * qJD(4) + t299;
t12 = t113 * t251 - t179 * t31 + t180 * t32;
t10 = t180 * t24 - t507;
t9 = t130 * t98 - t131 * t99 + t180 * t45 + t179 * t46 + (t101 * t226 + t102 * t225 - t178 * t215 - t214 * t487) * qJD(4);
t6 = t163 * t42 + t164 * t44 + t225 * t314 + t429 * t90 + t78 * t93 + t79 * t96;
t5 = t163 * t41 + t164 * t43 + t225 * t315 - t429 * t88 - t78 * t91 - t79 * t94;
t4 = -t160 * t44 + t161 * t42 + t226 * t314 - t430 * t90 + t76 * t93 + t77 * t96;
t3 = -t160 * t43 + t161 * t41 + t226 * t315 + t430 * t88 - t76 * t91 - t77 * t94;
t2 = t130 * t26 + t131 * t25 - t179 * t5 + t180 * t6 + t20 * t251 - t371 * t57;
t1 = t130 * t24 + t131 * t23 - t179 * t3 + t180 * t4 + t19 * t251 + t371 * t491;
t15 = [t451 / 0.2e1 + t452 / 0.2e1 + t100 * t367 - t458 / 0.2e1 + t457 / 0.2e1 + (t228 * t271 + t229 * t269) * qJD(1) + t11 * t470 + t20 * t471 - t491 * t473 + t57 * t474 + t318 + ((t25 + (-t225 * t88 + t226 * t90) * t269 + t492 + t512) * t180 + t10 + t507) * t472 + (t30 + t48) * t399 / 0.2e1 + (t14 * (t292 + t347) + t13 * (t361 + t446) + (t13 * pkin(6) + t14 * t296) * t226 + (-t356 + t296 * t214 + (-t459 + (rSges(6,3) + pkin(7)) * t271) * t399 + t478) * t36 + (t233 + t36 + t450 + (-t487 + t484) * qJD(1) + t477 - t508) * t37) * m(6) + (t50 * (t292 + t453) + t49 * (t147 + t361) + (t49 * pkin(6) + t319 * t50) * t226 + (t319 * t214 - t348 * t399 - t454 + t478) * t61 + (-t348 * t400 - t403 + t61 + (-t145 + t484) * qJD(1) + t477 + t87) * t62) * m(5) + (t115 * t486 + t114 * t496 + (t350 + t497) * t134 + (t134 + t408 + (-t294 + t304) * qJD(1) + t482) * (qJD(1) * t496 - t260)) * m(4) + (t152 * (t262 + t264 + t303) + t151 * t483 + (t258 + (t384 - t264 + t303) * qJD(1) + t482) * (qJD(1) * t483 - t260)) * m(3) - (t29 + t47 + t22) * t400 / 0.2e1 + (t21 - t503) * t366 + (t19 + t11) * t469 + (((t52 - t53 - t419) * t225 + t418 * t226) * t367 + t323 * qJD(1) + (t110 - t73) * t468 + (t109 - t72) * t467 + ((t53 + (t136 - t326) * t226) * t226 + (t54 + (t136 - t431) * t225 + t504 - t418) * t225) * t366) * qJD(4); 0.2e1 * (t13 * t390 + t14 * t389) * m(6) + 0.2e1 * (t389 * t50 + t390 * t49) * m(5) + 0.2e1 * (t114 * t390 + t115 * t389) * m(4) + 0.2e1 * (t151 * t390 + t152 * t389) * m(3); -m(5) * t27 - m(6) * t9; -qJD(1) * ((t406 * t269 - t405 * t271) * qJD(1) + ((t225 * t415 - t226 * t414) * t271 + (-t225 * t417 + t226 * t416) * t269) * qJD(4)) / 0.2e1 + ((t166 * t399 - t391) * t226 + (t284 + (-t476 * t225 + (-t165 + t295) * t226) * qJD(4)) * t225) * t366 + ((t165 * t400 + t391) * t225 + (t284 + (t295 * t226 + (-t166 - t476) * t225) * qJD(4)) * t226) * t367 + t12 * t393 / 0.2e1 + (-t225 * t31 + t226 * t32) * t360 + (t214 * t24 + t215 * t23 - t225 * t3 + t226 * t4) * t469 + ((t127 * t161 - t129 * t160) * t180 - (t126 * t161 - t128 * t160) * t179 + (-t160 * t198 + t161 * t196) * t251 + (-t24 * t425 + t269 * t491) * qJD(5) + ((-qJD(5) * t23 + t298) * t271 + t273) * t226) * t470 + (t214 * t26 + t215 * t25 - t225 * t5 + t226 * t6) * t471 + ((t127 * t163 + t129 * t164) * t180 - (t126 * t163 + t128 * t164) * t179 + (t163 * t196 + t164 * t198) * t251 + (-t25 * t423 - t269 * t57) * qJD(5) + ((-qJD(5) * t26 + t298) * t271 + t273) * t225) * t472 + (-t225 * t23 + t226 * t24) * t473 + (-t225 * t25 + t226 * t26) * t474 + (t214 * t32 + t215 * t31 - t225 * t7 + t226 * t8) * t463 + (((t127 * t268 - t129 * t270 - t90) * t180 - (t126 * t268 - t128 * t270 + t88) * t179 + (t196 * t268 - t198 * t270 + t285) * t251 - t113 * qJD(5)) * t269 + (qJD(5) * t343 - t475) * t271) * t464 + qJD(1) * (-t214 * t73 - t215 * t72 - t225 * t29 + t226 * t30) / 0.2e1 + (t226 * t10 + t225 * t11) * t392 / 0.2e1 + ((-t28 * t446 + t36 * t410) * t215 - (-t28 * t447 + t37 * t410) * t214 + (-t13 * t410 + t37 * t413 + t9 * t447 + t28 * (t101 + t45)) * t226 + (-t14 * t410 + t28 * t450 + t36 * t413 + t446 * t9) * t225 - t36 * (-qJD(1) * t175 - t132 * t251 - t179 * t200 + t500) - t37 * (qJD(1) * t177 + t133 * t251 - t180 * t200 + t499) - t28 * (t132 * t180 + t133 * t179 + t175 * t399 + t177 * t400) - ((t36 * t98 - t37 * t99) * t269 + t274 * t271) * qJD(5)) * m(6) + (t27 * t325 + t60 * t283 - t342 * t230 - (-t214 * t62 + t215 * t61 - t225 * t50 - t226 * t49) * t348 - (-t171 * t61 + t172 * t62) * qJD(1) - (t60 * (t171 * t226 + t172 * t225) - t342 * t349) * qJD(4)) * m(5) - (qJD(1) * t47 + t1 + (-(-t136 * t215 - t214 * t501 - t225 * t80) * t225 + (t138 * t215 + t214 * t326 - t225 * t81) * t226 + t214 * t52 + t215 * t51) * t490) * t225 / 0.2e1 + (qJD(1) * t48 + t2 + (-(-t136 * t214 + t215 * t501 + t226 * t80) * t225 + (t138 * t214 - t215 * t326 + t226 * t81) * t226 + t214 * t54 + t215 * t53) * t490) * t226 / 0.2e1 + (t306 + t22 + t11) * t468 + (t307 + t21 + t10) * t467; -t2 * t426 / 0.2e1 + (t269 * t344 + t271 * t57) * t474 + ((qJD(4) * t344 + t20) * t271 + (-qJD(4) * t57 - t214 * t25 + t215 * t26 - t225 * t6 - t226 * t5) * t269) * t471 - t1 * t424 / 0.2e1 + (t269 * t345 - t271 * t491) * t473 + ((qJD(4) * t345 + t19) * t271 + (qJD(4) * t491 - t214 * t23 + t215 * t24 - t225 * t4 - t226 * t3) * t269) * t469 + t12 * t364 + t271 * (t318 + t451 + t452 + t457 - t458) / 0.2e1 + (t113 * t271 + t269 * t343) * t360 + ((qJD(4) * t343 + t38) * t271 + (-qJD(4) * t113 - t214 * t31 + t215 * t32 - t225 * t8 - t226 * t7) * t269) * t463 + (t163 * t279 + t164 * t278 - t225 * t289) * t472 + (-t160 * t278 + t161 * t279 - t226 * t289) * t470 + (t297 * t271 + (t279 * t268 - t270 * t278) * t269) * t464 + (t429 / 0.2e1 + t225 * t363) * t11 + (-t430 / 0.2e1 + t226 * t363) * t10 + ((qJD(4) * t274 + t13 * t99 - t14 * t98 - t36 * t45 + t37 * t46) * t271 + (t36 * (qJD(4) * t98 - t153 * t226) + t37 * (-qJD(4) * t99 + t153 * t225) + t9 * t341 + t28 * (t214 * t99 + t215 * t98 - t225 * t45 + t226 * t46) - (t13 * t225 - t14 * t226 - t214 * t36 - t215 * t37) * t288) * t269 - t36 * (-t111 * t251 - t179 * t224) - t37 * (t112 * t251 - t180 * t224) - t28 * (t111 * t180 + t112 * t179)) * m(6);];
tauc = t15(:);
