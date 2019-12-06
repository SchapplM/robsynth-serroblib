% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:40
% EndTime: 2019-12-05 15:06:52
% DurationCPUTime: 8.40s
% Computational Cost: add. (20410->468), mult. (26023->682), div. (0->0), fcn. (27884->6), ass. (0->291)
t254 = pkin(8) + qJ(3);
t251 = cos(t254);
t256 = cos(pkin(7));
t259 = cos(qJ(4));
t353 = t256 * t259;
t255 = sin(pkin(7));
t258 = sin(qJ(4));
t356 = t255 * t258;
t235 = -t251 * t356 - t353;
t354 = t256 * t258;
t355 = t255 * t259;
t236 = t251 * t355 - t354;
t250 = sin(t254);
t368 = t250 * t255;
t128 = Icges(6,5) * t236 + Icges(6,6) * t235 + Icges(6,3) * t368;
t130 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t368;
t464 = t128 + t130;
t237 = -t251 * t354 + t355;
t238 = t251 * t353 + t356;
t367 = t250 * t256;
t129 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t367;
t131 = Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t367;
t463 = t129 + t131;
t380 = Icges(6,4) * t236;
t132 = Icges(6,2) * t235 + Icges(6,6) * t368 + t380;
t384 = Icges(5,4) * t236;
t134 = Icges(5,2) * t235 + Icges(5,6) * t368 + t384;
t478 = t132 + t134;
t379 = Icges(6,4) * t238;
t133 = Icges(6,2) * t237 + Icges(6,6) * t367 + t379;
t383 = Icges(5,4) * t238;
t135 = Icges(5,2) * t237 + Icges(5,6) * t367 + t383;
t477 = t133 + t135;
t221 = Icges(6,4) * t235;
t136 = Icges(6,1) * t236 + Icges(6,5) * t368 + t221;
t223 = Icges(5,4) * t235;
t138 = Icges(5,1) * t236 + Icges(5,5) * t368 + t223;
t476 = t136 + t138;
t222 = Icges(6,4) * t237;
t137 = Icges(6,1) * t238 + Icges(6,5) * t367 + t222;
t224 = Icges(5,4) * t237;
t139 = Icges(5,1) * t238 + Icges(5,5) * t367 + t224;
t475 = t137 + t139;
t377 = Icges(6,4) * t259;
t296 = -Icges(6,2) * t258 + t377;
t194 = -Icges(6,6) * t251 + t296 * t250;
t381 = Icges(5,4) * t259;
t297 = -Icges(5,2) * t258 + t381;
t196 = -Icges(5,6) * t251 + t297 * t250;
t461 = t194 + t196;
t378 = Icges(6,4) * t258;
t300 = Icges(6,1) * t259 - t378;
t198 = -Icges(6,5) * t251 + t300 * t250;
t382 = Icges(5,4) * t258;
t301 = Icges(5,1) * t259 - t382;
t200 = -Icges(5,5) * t251 + t301 * t250;
t460 = t198 + t200;
t151 = Icges(6,5) * t235 - Icges(6,6) * t236;
t153 = Icges(5,5) * t235 - Icges(5,6) * t236;
t474 = t151 + t153;
t152 = Icges(6,5) * t237 - Icges(6,6) * t238;
t154 = Icges(5,5) * t237 - Icges(5,6) * t238;
t473 = t152 + t154;
t293 = Icges(6,5) * t259 - Icges(6,6) * t258;
t190 = -Icges(6,3) * t251 + t250 * t293;
t294 = Icges(5,5) * t259 - Icges(5,6) * t258;
t192 = -Icges(5,3) * t251 + t250 * t294;
t462 = t190 + t192;
t342 = -Icges(5,2) * t238 + t139 + t224;
t344 = -Icges(6,2) * t238 + t137 + t222;
t472 = t342 + t344;
t343 = -Icges(5,2) * t236 + t138 + t223;
t345 = -Icges(6,2) * t236 + t136 + t221;
t471 = t343 + t345;
t346 = Icges(5,1) * t237 - t135 - t383;
t348 = Icges(6,1) * t237 - t133 - t379;
t470 = t346 + t348;
t347 = Icges(5,1) * t235 - t134 - t384;
t349 = Icges(6,1) * t235 - t132 - t380;
t469 = t347 + t349;
t447 = rSges(6,1) + pkin(4);
t468 = t477 * t235 + t475 * t236 + t463 * t368;
t467 = t478 * t237 + t476 * t238 + t464 * t367;
t306 = rSges(6,1) * t259 - rSges(6,2) * t258;
t397 = pkin(4) * t259;
t437 = (rSges(6,3) + qJ(5)) * t251 + (-t306 - t397) * t250;
t466 = t437 * t255;
t465 = t437 * t256;
t252 = t255 ^ 2;
t253 = t256 ^ 2;
t421 = t252 + t253;
t459 = -t461 * t258 + t460 * t259;
t458 = Icges(4,5) * t250 + Icges(4,6) * t251;
t457 = t471 * t235 + t469 * t236 + t474 * t368;
t456 = t472 * t235 + t470 * t236 + t473 * t368;
t455 = t471 * t237 + t469 * t238 + t474 * t367;
t454 = t472 * t237 + t470 * t238 + t473 * t367;
t307 = rSges(5,1) * t259 - rSges(5,2) * t258;
t203 = -t251 * rSges(5,3) + t307 * t250;
t183 = t203 * t255;
t185 = t203 * t256;
t451 = t460 + (-t378 - t382 + (-Icges(5,2) - Icges(6,2)) * t259) * t250;
t450 = t461 + (t377 + t381 + (Icges(5,1) + Icges(6,1)) * t258) * t250;
t227 = (-Icges(6,5) * t258 - Icges(6,6) * t259) * t250;
t228 = (-Icges(5,5) * t258 - Icges(5,6) * t259) * t250;
t449 = (-t227 - t228) * t251;
t432 = t462 * t251;
t174 = t194 * t255;
t176 = t196 * t255;
t178 = t198 * t255;
t180 = t200 * t255;
t289 = -t134 * t258 + t138 * t259;
t279 = -t192 * t255 - t289;
t291 = -t132 * t258 + t136 * t259;
t281 = -t190 * t255 - t291;
t446 = (-t281 - t279) * t251 + ((-t178 - t180) * t259 + (t174 + t176) * t258 + t464) * t250;
t175 = t194 * t256;
t177 = t196 * t256;
t179 = t198 * t256;
t181 = t200 * t256;
t288 = -t135 * t258 + t139 * t259;
t278 = -t192 * t256 - t288;
t290 = -t133 * t258 + t137 * t259;
t280 = -t190 * t256 - t290;
t445 = (-t280 - t278) * t251 + ((-t179 - t181) * t259 + (t175 + t177) * t258 + t463) * t250;
t444 = t478 * t235 + t476 * t236 + t464 * t368;
t443 = t477 * t237 + t475 * t238 + t463 * t367;
t442 = t461 * t235 + t460 * t236 + t462 * t368;
t441 = t461 * t237 + t460 * t238 + t462 * t367;
t440 = t459 * t250 - t432;
t439 = (-t296 - t297) * t251 + (-Icges(5,6) - Icges(6,6)) * t250;
t438 = (-t300 - t301) * t251 + (-Icges(5,5) - Icges(6,5)) * t250;
t436 = ((t294 + t293) * t251 + (Icges(5,3) + Icges(6,3)) * t250 - t459) * t251;
t366 = t251 * t128;
t92 = t250 * t291 - t366;
t365 = t251 * t129;
t93 = t250 * t290 - t365;
t364 = t251 * t130;
t94 = t250 * t289 - t364;
t363 = t251 * t131;
t95 = t250 * t288 - t363;
t433 = (t93 + t95) * t256 + (t92 + t94) * t255;
t431 = t467 * t255;
t430 = t468 * t256;
t338 = -rSges(6,2) * t238 + t447 * t237;
t339 = -rSges(6,2) * t236 + t447 * t235;
t101 = t339 * t255 + t338 * t256;
t164 = rSges(5,1) * t235 - rSges(5,2) * t236;
t166 = rSges(5,1) * t237 - rSges(5,2) * t238;
t118 = t164 * t255 + t166 * t256;
t242 = pkin(3) * t250 - pkin(6) * t251;
t314 = -t242 + t437;
t119 = t314 * t255;
t121 = t314 * t256;
t331 = -t203 - t242;
t145 = t331 * t255;
t147 = t331 * t256;
t308 = (rSges(6,2) * t259 + t447 * t258) * t250;
t168 = t308 * t255;
t169 = t308 * t256;
t234 = (-rSges(5,1) * t258 - rSges(5,2) * t259) * t250;
t243 = pkin(3) * t251 + pkin(6) * t250;
t326 = t421 * t243;
t187 = qJ(5) * t250 + t397 * t251;
t350 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t367 + pkin(4) * t356 + t187 * t256;
t351 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t368 - pkin(4) * t354 + t187 * t255;
t69 = t351 * t255 + t350 * t256 + t326;
t141 = rSges(5,1) * t236 + rSges(5,2) * t235 + rSges(5,3) * t368;
t143 = rSges(5,1) * t238 + rSges(5,2) * t237 + rSges(5,3) * t367;
t99 = t141 * t255 + t143 * t256 + t326;
t429 = -m(6) * (t101 * t69 + t119 * t168 + t121 * t169) - m(5) * (t118 * t99 + (-t145 * t255 - t147 * t256) * t234);
t428 = -t250 / 0.2e1;
t427 = t250 / 0.2e1;
t426 = -t251 / 0.2e1;
t425 = t251 / 0.2e1;
t424 = -t255 / 0.2e1;
t401 = t255 / 0.2e1;
t400 = -t256 / 0.2e1;
t399 = t256 / 0.2e1;
t423 = (-t451 * t235 + t450 * t236) * t251 + (t456 * t256 + (t449 + t457) * t255) * t250;
t422 = (-t451 * t237 + t450 * t238) * t251 + ((t449 + t454) * t256 + t455 * t255) * t250;
t419 = t458 * t421;
t417 = t251 ^ 2;
t416 = 2 * qJD(3);
t415 = 2 * qJD(4);
t414 = 4 * qJD(4);
t413 = m(5) / 0.2e1;
t412 = m(6) / 0.2e1;
t357 = t251 * t256;
t358 = t251 * t255;
t90 = -t250 * t466 + t351 * t251;
t91 = t250 * t465 - t350 * t251;
t395 = t90 * t357 + t91 * t358;
t273 = -t350 * t255 + t351 * t256;
t44 = t273 * t251 + (-t255 * t465 + t256 * t466) * t250;
t336 = t250 * rSges(6,3) + t306 * t251 + t187;
t46 = (t336 * t255 - t351) * t250;
t47 = (-t336 * t256 + t350) * t250;
t74 = t273 * t250;
t411 = m(6) * (-t251 * t44 + (t255 * t47 + t256 * t46 + t74) * t250 + t395);
t410 = m(6) * (t44 * t74 + t46 * t90 + t47 * t91);
t287 = t141 * t256 - t143 * t255;
t107 = t287 * t250;
t116 = t141 * t251 + t203 * t368;
t117 = -t143 * t251 - t203 * t367;
t79 = t287 * t251 + (-t183 * t256 + t185 * t255) * t250;
t205 = t250 * rSges(5,3) + t307 * t251;
t96 = (t205 * t255 - t141) * t250;
t97 = (-t205 * t256 + t143) * t250;
t409 = m(5) * (t107 * t79 + t116 * t96 + t117 * t97);
t239 = t421 * t250;
t407 = m(6) * (t239 * t74 + t395);
t405 = m(6) * (-t101 * t251 + (t168 * t255 + t169 * t256) * t250);
t404 = t239 / 0.2e1;
t394 = m(6) * qJD(3);
t393 = m(6) * qJD(5);
t369 = t250 * t251;
t352 = t119 * t358 + t121 * t357;
t330 = -t205 - t243;
t327 = t421 * t242;
t325 = t421 * t369;
t324 = qJD(4) * t250;
t189 = (t404 + t428) * m(6);
t323 = t189 * qJD(1);
t264 = t250 * t281 + t366;
t49 = -t235 * t174 - t236 * t178 + t255 * t264;
t263 = t250 * t280 + t365;
t50 = -t235 * t175 - t236 * t179 + t255 * t263;
t262 = t250 * t279 + t364;
t51 = -t235 * t176 - t236 * t180 + t255 * t262;
t261 = t250 * t278 + t363;
t52 = -t235 * t177 - t236 * t181 + t255 * t261;
t322 = ((t52 + t50) * t256 + (t51 + t49 - t436) * t255 + t442) * t427 + ((-t432 + t444) * t255 + t438 * t236 + t439 * t235 + t430) * t425;
t53 = -t237 * t174 - t238 * t178 + t256 * t264;
t54 = -t237 * t175 - t238 * t179 + t256 * t263;
t55 = -t237 * t176 - t238 * t180 + t256 * t262;
t56 = -t237 * t177 - t238 * t181 + t256 * t261;
t321 = ((t56 + t54 - t436) * t256 + (t55 + t53) * t255 + t441) * t427 + ((-t432 + t443) * t256 + t438 * t238 + t439 * t237 + t431) * t425;
t320 = (t445 * t256 + t446 * t255 + (-t258 * t439 + t259 * t438 - t462) * t251 + t440) * t428 + (t433 + t436) * t426;
t319 = t457 * t399 + t456 * t424;
t318 = t455 * t400 + t454 * t401;
t317 = (t255 * t444 + t430) * t427 + t442 * t426;
t316 = (t256 * t443 + t431) * t427 + t441 * t426;
t315 = t425 * t440 + t428 * t433;
t313 = -t243 - t336;
t240 = rSges(4,1) * t250 + rSges(4,2) * t251;
t260 = -m(5) * (t255 * t96 - t256 * t97) / 0.2e1 - m(6) * (t46 * t255 - t47 * t256) / 0.2e1;
t277 = (-t168 * t256 + t169 * t255) * t412;
t16 = t277 + t260;
t17 = 0.2e1 * (t44 / 0.4e1 - t101 / 0.4e1) * m(6) + 0.2e1 * (t79 / 0.4e1 - t118 / 0.4e1) * m(5);
t292 = -t17 * qJD(1) + t16 * qJD(2);
t284 = -t319 - t322;
t283 = -t318 + t321;
t282 = t308 * t250;
t216 = t458 * t256;
t215 = t458 * t255;
t188 = m(6) * t404 + t250 * t412;
t167 = t325 - t369;
t148 = t330 * t256;
t146 = t330 * t255;
t144 = t421 * t240;
t125 = -t166 * t251 - t234 * t367;
t124 = t164 * t251 + t234 * t368;
t122 = t313 * t256;
t120 = t313 * t255;
t112 = (t164 * t256 - t166 * t255) * t250;
t109 = -t338 * t251 + t282 * t256;
t108 = t339 * t251 - t282 * t255;
t106 = -t183 * t255 - t185 * t256 - t327;
t98 = (-t338 * t255 + t339 * t256) * t250;
t77 = t255 * t466 + t256 * t465 - t327;
t75 = t405 / 0.2e1;
t73 = -t251 * t154 + (-t342 * t258 + t346 * t259) * t250;
t72 = -t251 * t153 + (-t343 * t258 + t347 * t259) * t250;
t71 = -t251 * t152 + (-t344 * t258 + t348 * t259) * t250;
t70 = -t251 * t151 + (-t345 * t258 + t349 * t259) * t250;
t42 = t239 * t69 + t352;
t29 = t407 / 0.2e1;
t28 = t255 * t56 - t256 * t55;
t27 = t255 * t54 - t256 * t53;
t26 = t255 * t52 - t256 * t51;
t25 = t255 * t50 - t256 * t49;
t18 = (t118 + t79) * t413 + (t101 + t44) * t412;
t15 = t277 - t260;
t6 = t411 / 0.2e1;
t5 = t29 + t75 - t411 / 0.2e1;
t4 = t29 + t6 - t405 / 0.2e1;
t3 = t75 + t6 - t407 / 0.2e1;
t2 = t318 * t255 + t319 * t256 - t429;
t1 = t409 + t410 + (t317 * t255 + t316 * t256 + t320) * t251 + (t322 * t255 + t321 * t256 - t315) * t250;
t7 = [0, 0, t18 * qJD(4) + t188 * qJD(5) + (-m(4) * t144 / 0.2e1 + t106 * t413 + t77 * t412) * t416, t18 * qJD(3) + (t112 * t413 + t98 * t412) * t415, t188 * qJD(3); 0, 0, t15 * qJD(4) + ((-t146 * t256 + t148 * t255) * t413 + (-t120 * t256 + t122 * t255) * t412) * t416, t15 * qJD(3) + ((t124 * t255 - t125 * t256) * t413 + (t108 * t255 - t109 * t256) * t412) * t415, 0; -qJD(4) * t17 + qJD(5) * t189, t16 * qJD(4), t2 * qJD(4) + t42 * t393 + (m(6) * (t119 * t120 + t121 * t122 + t69 * t77) + m(5) * (t106 * t99 + t145 * t146 + t147 * t148) + m(4) * (-t144 + t240) * t421 * (rSges(4,1) * t251 - rSges(4,2) * t250) + (t28 + t27 - t252 * t216 + (t255 * t215 - t419) * t256) * t401 + (t25 + t26 - t253 * t215 + (t256 * t216 - t419) * t255) * t400) * qJD(3), t2 * qJD(3) + t5 * qJD(5) + (-t409 / 0.4e1 - t410 / 0.4e1) * t414 + ((t107 * t118 + t112 * t99 + t124 * t147 + t125 * t145 + (-t116 * t256 - t117 * t255) * t234) * t413 + (t101 * t74 + t108 * t121 + t109 * t119 + t168 * t91 + t169 * t90 + t69 * t98) * t412) * t415 + (t255 * t284 - t256 * t283 + t315) * t324 + t292 + (((t70 / 0.2e1 + t72 / 0.2e1 - t316) * t256 + (-t71 / 0.2e1 - t73 / 0.2e1 - t317) * t255 - t320) * t251 + t422 * t401 + t423 * t400) * qJD(4), t323 + t42 * t394 + t5 * qJD(4) + (-t239 * t251 - t167 + t325) * t393; t17 * qJD(3), -t16 * qJD(3), t1 * qJD(4) + t4 * qJD(5) + ((t106 * t107 + t116 * t148 + t117 * t146 + t145 * t97 + t147 * t96 + t79 * t99) * t413 + (t119 * t47 + t120 * t91 + t121 * t46 + t122 * t90 + t44 * t69 + t74 * t77) * t412) * t416 - t292 + (t284 * t256 + t283 * t255 + (t445 * t424 + (t468 * t255 - t444 * t256) * t401 + (t443 * t255 - t467 * t256 + t446) * t399) * t251 + ((-t92 / 0.2e1 - t94 / 0.2e1 + t28 / 0.2e1 + t27 / 0.2e1) * t256 + (t93 / 0.2e1 + t95 / 0.2e1 + t26 / 0.2e1 + t25 / 0.2e1) * t255) * t250 + t429) * qJD(3), t1 * qJD(3) + (m(6) * (t108 * t90 + t109 * t91 + t74 * t98) / 0.4e1 + m(5) * (t107 * t112 + t116 * t124 + t117 * t125) / 0.4e1) * t414 + (-t227 / 0.2e1 - t228 / 0.2e1) * qJD(4) * t251 * t417 + (((t71 + t73) * t256 + (t70 + t72) * t255 + (t451 * t258 + t450 * t259) * t251) * t426 + t423 * t401 + t422 * t399) * t324, t4 * qJD(3); -t189 * qJD(3), 0, -t323 + (-t251 * t77 + (t120 * t255 + t122 * t256 + t69) * t250 - t42 + t352) * t394 + t3 * qJD(4) + t167 * t393, t3 * qJD(3) + m(6) * (-t251 * t98 + (t108 * t256 + t109 * t255) * t250) * qJD(4), t167 * t394;];
Cq = t7;
