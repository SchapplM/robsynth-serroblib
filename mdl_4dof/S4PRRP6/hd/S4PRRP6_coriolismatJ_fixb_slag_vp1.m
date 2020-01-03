% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:25
% DurationCPUTime: 8.76s
% Computational Cost: add. (9856->491), mult. (26058->724), div. (0->0), fcn. (28608->6), ass. (0->300)
t444 = rSges(5,1) + pkin(3);
t437 = rSges(5,3) + qJ(4);
t250 = sin(pkin(6));
t251 = cos(pkin(6));
t254 = cos(qJ(3));
t252 = sin(qJ(3));
t255 = cos(qJ(2));
t347 = t255 * t252;
t219 = t250 * t347 + t251 * t254;
t348 = t254 * t255;
t358 = t251 * t252;
t220 = t250 * t348 - t358;
t253 = sin(qJ(2));
t359 = t250 * t253;
t129 = Icges(4,5) * t220 - Icges(4,6) * t219 + Icges(4,3) * t359;
t131 = Icges(5,4) * t220 + Icges(5,2) * t359 + Icges(5,6) * t219;
t474 = t129 + t131;
t221 = -t250 * t254 + t251 * t347;
t360 = t250 * t252;
t222 = t251 * t348 + t360;
t357 = t251 * t253;
t130 = Icges(4,5) * t222 - Icges(4,6) * t221 + Icges(4,3) * t357;
t132 = Icges(5,4) * t222 + Icges(5,2) * t357 + Icges(5,6) * t221;
t473 = t130 + t132;
t382 = Icges(4,4) * t254;
t291 = -Icges(4,2) * t252 + t382;
t201 = -Icges(4,6) * t255 + t291 * t253;
t349 = t253 * t254;
t246 = Icges(5,5) * t349;
t356 = t252 * t253;
t374 = Icges(5,6) * t255;
t459 = Icges(5,3) * t356 - t201 + t246 - t374;
t288 = Icges(4,5) * t254 - Icges(4,6) * t252;
t197 = -Icges(4,3) * t255 + t288 * t253;
t290 = Icges(5,4) * t254 + Icges(5,6) * t252;
t199 = -Icges(5,2) * t255 + t290 * t253;
t472 = t197 + t199;
t377 = Icges(5,5) * t252;
t294 = Icges(5,1) * t254 + t377;
t203 = -Icges(5,4) * t255 + t294 * t253;
t383 = Icges(4,4) * t252;
t295 = Icges(4,1) * t254 - t383;
t205 = -Icges(4,5) * t255 + t295 * t253;
t458 = t203 + t205;
t249 = t253 ^ 2;
t412 = t255 ^ 2;
t149 = -Icges(4,5) * t219 - Icges(4,6) * t220;
t151 = -Icges(5,4) * t219 + Icges(5,6) * t220;
t471 = t149 + t151;
t150 = -Icges(4,5) * t221 - Icges(4,6) * t222;
t152 = -Icges(5,4) * t221 + Icges(5,6) * t222;
t470 = t150 + t152;
t218 = Icges(4,4) * t221;
t138 = Icges(4,1) * t222 + Icges(4,5) * t357 - t218;
t339 = Icges(4,2) * t222 - t138 + t218;
t378 = Icges(5,5) * t221;
t136 = Icges(5,1) * t222 + Icges(5,4) * t357 + t378;
t341 = Icges(5,3) * t222 - t136 - t378;
t468 = t339 + t341;
t217 = Icges(4,4) * t219;
t137 = Icges(4,1) * t220 + Icges(4,5) * t359 - t217;
t340 = Icges(4,2) * t220 - t137 + t217;
t379 = Icges(5,5) * t219;
t135 = Icges(5,1) * t220 + Icges(5,4) * t359 + t379;
t342 = Icges(5,3) * t220 - t135 - t379;
t467 = t340 + t342;
t384 = Icges(4,4) * t222;
t134 = -Icges(4,2) * t221 + Icges(4,6) * t357 + t384;
t343 = -Icges(4,1) * t221 - t134 - t384;
t216 = Icges(5,5) * t222;
t128 = Icges(5,6) * t357 + Icges(5,3) * t221 + t216;
t345 = -Icges(5,1) * t221 + t128 + t216;
t466 = t343 + t345;
t385 = Icges(4,4) * t220;
t133 = -Icges(4,2) * t219 + Icges(4,6) * t359 + t385;
t344 = -Icges(4,1) * t219 - t133 - t385;
t215 = Icges(5,5) * t220;
t127 = Icges(5,6) * t359 + Icges(5,3) * t219 + t215;
t346 = -Icges(5,1) * t219 + t127 + t215;
t465 = t344 + t346;
t463 = t127 - t133;
t462 = t128 - t134;
t461 = t135 + t137;
t460 = t136 + t138;
t457 = t459 * t252 + t458 * t254;
t456 = t437 * t252 + t444 * t254;
t455 = t472 * t253;
t454 = t473 * t253;
t453 = t474 * t253;
t452 = t467 * t219 + t465 * t220 + t471 * t359;
t451 = t468 * t219 + t466 * t220 + t470 * t359;
t450 = t467 * t221 + t465 * t222 + t471 * t357;
t449 = t468 * t221 + t466 * t222 + t470 * t357;
t448 = t458 + (t377 - t383 + (-Icges(4,2) - Icges(5,3)) * t254) * t253;
t447 = -(-Icges(4,1) * t252 - t382) * t253 + Icges(5,1) * t356 - t246 - t459;
t230 = (-Icges(4,5) * t252 - Icges(4,6) * t254) * t253;
t231 = (-Icges(5,4) * t252 + Icges(5,6) * t254) * t253;
t446 = (-t230 - t231) * t255;
t426 = t472 * t255;
t445 = 0.2e1 * t253 * t255 * (Icges(3,1) - Icges(3,2)) + (0.2e1 * t412 - 0.2e1 * t249) * Icges(3,4);
t287 = Icges(5,5) * t254 + Icges(5,3) * t252;
t261 = -t287 * t253 + t374;
t173 = t261 * t250;
t179 = t201 * t250;
t181 = t203 * t250;
t183 = t205 * t250;
t283 = -t133 * t252 + t137 * t254;
t274 = -t197 * t250 - t283;
t285 = t127 * t252 + t135 * t254;
t276 = t199 * t250 + t285;
t443 = (t276 - t274) * t255 + ((-t181 - t183) * t254 + (t173 + t179) * t252 + t474) * t253;
t174 = t261 * t251;
t180 = t201 * t251;
t182 = t203 * t251;
t184 = t205 * t251;
t282 = -t134 * t252 + t138 * t254;
t273 = -t197 * t251 - t282;
t284 = t128 * t252 + t136 * t254;
t275 = t199 * t251 + t284;
t442 = (-t273 + t275) * t255 + ((-t182 - t184) * t254 + (t174 + t180) * t252 + t473) * t253;
t441 = t219 * t463 + t220 * t461 + t250 * t453;
t440 = t219 * t462 + t220 * t460 + t250 * t454;
t439 = t221 * t463 + t222 * t461 + t251 * t453;
t438 = t221 * t462 + t222 * t460 + t251 * t454;
t289 = -Icges(3,5) * t253 - Icges(3,6) * t255;
t223 = t289 * t250;
t224 = t289 * t251;
t436 = t219 * t459 + t220 * t458 + t250 * t455;
t435 = t221 * t459 + t222 * t458 + t251 * t455;
t434 = t457 * t253 - t426;
t433 = (t287 - t291) * t255 + (-Icges(4,6) + Icges(5,6)) * t253;
t432 = (t294 + t295) * t255 + (Icges(5,4) + Icges(4,5)) * t253;
t431 = rSges(5,2) * t255 - t456 * t253;
t430 = ((t288 + t290) * t255 + (Icges(5,2) + Icges(4,3)) * t253 - t457) * t255;
t370 = t131 * t255;
t88 = t285 * t253 - t370;
t369 = t132 * t255;
t89 = t284 * t253 - t369;
t372 = t129 * t255;
t90 = t283 * t253 - t372;
t371 = t130 * t255;
t91 = t282 * t253 - t371;
t427 = (t89 + t91) * t251 + (t88 + t90) * t250;
t165 = -rSges(4,1) * t219 - rSges(4,2) * t220;
t169 = -rSges(4,1) * t221 - rSges(4,2) * t222;
t107 = t165 * t250 + t169 * t251;
t242 = t253 * pkin(2) - pkin(5) * t255;
t309 = -t242 + t431;
t118 = t309 * t250;
t120 = t309 * t251;
t322 = (-t444 * t252 + t437 * t254) * t253;
t145 = t322 * t250;
t146 = t322 * t251;
t303 = rSges(4,1) * t254 - rSges(4,2) * t252;
t208 = -rSges(4,3) * t255 + t303 * t253;
t327 = -t208 - t242;
t159 = t327 * t250;
t161 = t327 * t251;
t237 = (-rSges(4,1) * t252 - rSges(4,2) * t254) * t253;
t243 = pkin(2) * t255 + t253 * pkin(5);
t410 = t251 ^ 2;
t411 = t250 ^ 2;
t416 = t410 + t411;
t323 = t416 * t243;
t337 = rSges(5,2) * t357 + t437 * t221 + t444 * t222;
t338 = rSges(5,2) * t359 + t437 * t219 + t444 * t220;
t72 = t338 * t250 + t337 * t251 + t323;
t335 = -t444 * t221 + t437 * t222;
t336 = -t444 * t219 + t437 * t220;
t86 = t336 * t250 + t335 * t251;
t140 = rSges(4,1) * t220 - rSges(4,2) * t219 + rSges(4,3) * t359;
t142 = rSges(4,1) * t222 - rSges(4,2) * t221 + rSges(4,3) * t357;
t94 = t140 * t250 + t142 * t251 + t323;
t425 = -m(5) * (-t118 * t145 - t120 * t146 + t72 * t86) - m(4) * (t107 * t94 + (-t159 * t250 - t161 * t251) * t237);
t424 = -t250 / 0.2e1;
t396 = t250 / 0.2e1;
t395 = -t251 / 0.2e1;
t394 = t251 / 0.2e1;
t423 = -t253 / 0.2e1;
t422 = t253 / 0.2e1;
t421 = -t255 / 0.2e1;
t419 = (t448 * t219 + t447 * t220) * t255 + (t451 * t251 + (t446 + t452) * t250) * t253;
t418 = (t448 * t221 + t447 * t222) * t255 + ((t446 + t449) * t251 + t450 * t250) * t253;
t143 = t219 * t250 + t221 * t251;
t115 = (-t143 + t347) * t356;
t306 = t347 / 0.2e1;
t122 = (t306 - t143 / 0.2e1) * m(5);
t390 = m(5) * qJD(4);
t417 = t122 * qJD(1) + t115 * t390;
t326 = t253 * rSges(5,2) + t255 * t456;
t414 = t326 * t253 - t255 * t431;
t409 = 2 * qJD(2);
t408 = 2 * qJD(3);
t407 = 4 * qJD(3);
t406 = m(4) / 0.2e1;
t405 = m(5) / 0.2e1;
t305 = t338 * t255;
t92 = -t359 * t431 + t305;
t304 = t337 * t255;
t93 = t357 * t431 - t304;
t298 = -t250 * t93 - t251 * t92;
t333 = t431 * t251;
t334 = t431 * t250;
t45 = (t334 * t253 + t305) * t251 + (-t333 * t253 - t304) * t250;
t66 = t414 * t250 - t338 * t253 + t334 * t255;
t67 = -t414 * t251 + t337 * t253 - t333 * t255;
t73 = (-t337 * t250 + t338 * t251) * t253;
t404 = m(5) * (t219 * t67 + t221 * t66 + (t255 * t73 + (t298 + t45) * t253) * t252);
t403 = m(5) * (t45 * t73 + t66 * t92 + t67 * t93);
t100 = (t140 * t251 - t142 * t250) * t253;
t368 = t140 * t255;
t111 = t208 * t359 + t368;
t367 = t142 * t255;
t112 = -t208 * t357 - t367;
t186 = t208 * t250;
t188 = t208 * t251;
t84 = (-t186 * t253 + t368) * t251 + (t188 * t253 - t367) * t250;
t210 = t253 * rSges(4,3) + t303 * t255;
t279 = t208 * t255 + t210 * t253;
t95 = -t253 * t140 - t186 * t255 + t279 * t250;
t96 = t253 * t142 + t188 * t255 - t279 * t251;
t402 = m(4) * (t100 * t84 + t111 * t95 + t112 * t96);
t125 = (t219 * t251 - t221 * t250) * t253;
t171 = t219 * t255 + t249 * t360;
t172 = -t221 * t255 - t249 * t358;
t400 = m(5) * (t118 * t172 + t120 * t171 + t125 * t72 + t143 * t73 + t298 * t356);
t399 = m(5) * (t118 * t220 + t120 * t222 - t145 * t219 - t146 * t221 + (t252 * t86 + t254 * t72) * t253);
t392 = m(5) * qJD(2);
t391 = m(5) * qJD(3);
t17 = 0.2e1 * (t45 / 0.4e1 - t86 / 0.4e1) * m(5) + 0.2e1 * (t84 / 0.4e1 - t107 / 0.4e1) * m(4);
t366 = t17 * qJD(1);
t325 = -t210 - t243;
t324 = t416 * t242;
t321 = qJD(3) * t253;
t307 = t349 / 0.2e1;
t116 = (t307 - t125 / 0.2e1) * m(5);
t320 = t116 * qJD(1);
t257 = -t276 * t253 + t370;
t46 = t219 * t173 - t220 * t181 + t250 * t257;
t256 = -t275 * t253 + t369;
t47 = t219 * t174 - t220 * t182 + t250 * t256;
t259 = t274 * t253 + t372;
t48 = t219 * t179 - t220 * t183 + t250 * t259;
t258 = t273 * t253 + t371;
t49 = t219 * t180 - t220 * t184 + t250 * t258;
t318 = t436 * t422 + (t219 * t433 + t220 * t432) * t421 + ((-t426 + t441) * t255 + (t46 + t48 - t430) * t253) * t396 + (t440 * t255 + (t47 + t49) * t253) * t394;
t50 = t221 * t173 - t222 * t181 + t251 * t257;
t51 = t221 * t174 - t222 * t182 + t251 * t256;
t52 = t221 * t179 - t222 * t183 + t251 * t259;
t53 = t221 * t180 - t222 * t184 + t251 * t258;
t316 = t435 * t422 + (t433 * t221 + t222 * t432) * t421 + (t439 * t255 + (t50 + t52) * t253) * t396 + ((-t426 + t438) * t255 + (t53 + t51 - t430) * t253) * t394;
t315 = ((-t252 * t433 - t254 * t432 - t472) * t255 + t442 * t251 + t443 * t250 + t434) * t423 + (t427 + t430) * t421;
t314 = t452 * t394 + t451 * t424;
t313 = t450 * t395 + t449 * t396;
t312 = (t441 * t250 + t440 * t251) * t422 + t436 * t421;
t311 = (t439 * t250 + t438 * t251) * t422 + t435 * t421;
t310 = t427 * t423 + t434 * t255 / 0.2e1;
t308 = -t243 - t326;
t240 = t253 * rSges(3,1) + rSges(3,2) * t255;
t286 = -t118 * t250 - t120 * t251;
t278 = -t314 - t318;
t277 = -t313 + t316;
t269 = t445 * t250 + t224;
t268 = -t445 * t251 + t223;
t162 = t325 * t251;
t160 = t325 * t250;
t144 = t416 * t240;
t123 = m(5) * t306 + t143 * t405;
t121 = t308 * t251;
t119 = t308 * t250;
t117 = m(5) * t307 + t125 * t405;
t114 = -t169 * t255 - t237 * t357;
t113 = t165 * t255 + t237 * t359;
t110 = t249 * t252 * t254 + t219 * t220 + t221 * t222;
t105 = (t165 * t251 - t169 * t250) * t253;
t99 = -t186 * t250 - t188 * t251 - t324;
t98 = -t253 * t146 - t335 * t255;
t97 = t336 * t255 + t322 * t359;
t87 = t334 * t250 + t333 * t251 - t324;
t83 = (-t335 * t250 + t336 * t251) * t253;
t71 = -t150 * t255 + (t339 * t252 + t343 * t254) * t253;
t70 = -t149 * t255 + (t340 * t252 + t344 * t254) * t253;
t69 = -t152 * t255 + (t341 * t252 + t345 * t254) * t253;
t68 = -t151 * t255 + (t342 * t252 + t346 * t254) * t253;
t42 = t143 * t72 + t286 * t356;
t31 = t125 * t73 + t171 * t92 + t172 * t93;
t29 = t250 * t53 - t251 * t52;
t28 = t250 * t51 - t251 * t50;
t27 = t250 * t49 - t251 * t48;
t26 = t250 * t47 - t251 * t46;
t24 = t399 / 0.2e1;
t18 = (t107 + t84) * t406 + (t45 + t86) * t405;
t15 = t400 / 0.2e1;
t6 = t404 / 0.2e1;
t5 = t24 + t6 - t400 / 0.2e1;
t4 = t15 + t24 - t404 / 0.2e1;
t3 = t15 + t6 - t399 / 0.2e1;
t2 = t313 * t250 + t314 * t251 - t425;
t1 = t402 + t403 + (t312 * t250 + t311 * t251 + t315) * t255 + (t318 * t250 + t316 * t251 - t310) * t253;
t7 = [0, t18 * qJD(3) + t123 * qJD(4) + (-m(3) * t144 / 0.2e1 + t99 * t406 + t87 * t405) * t409, t18 * qJD(2) + (t105 * t406 + t83 * t405) * t408 + t117 * qJD(4), qJD(2) * t123 + qJD(3) * t117; -qJD(3) * t17 - qJD(4) * t122, t2 * qJD(3) + t42 * t390 + (m(5) * (t118 * t119 + t120 * t121 + t72 * t87) + m(4) * (t159 * t160 + t161 * t162 + t94 * t99) + m(3) * (-t144 + t240) * t416 * (rSges(3,1) * t255 - t253 * rSges(3,2)) + (t29 + t28 + t411 * t224 + (t269 * t251 + (-t223 + t268) * t250) * t251) * t396 + (t26 + t27 + t410 * t223 + (t268 * t250 + (-t224 + t269) * t251) * t250) * t395) * qJD(2), -t366 + t2 * qJD(2) + t4 * qJD(4) + (-t403 / 0.4e1 - t402 / 0.4e1) * t407 + ((t100 * t107 + t105 * t94 + t113 * t161 + t114 * t159 + (-t111 * t251 - t112 * t250) * t237) * t406 + (t118 * t98 + t120 * t97 - t145 * t93 - t146 * t92 + t72 * t83 + t73 * t86) * t405) * t408 + (t250 * t278 - t251 * t277 + t310) * t321 + (((t70 / 0.2e1 + t68 / 0.2e1 - t311) * t251 + (-t71 / 0.2e1 - t69 / 0.2e1 - t312) * t250 - t315) * t255 + t418 * t396 + t419 * t395) * qJD(3), t4 * qJD(3) + t42 * t392 - t417; qJD(2) * t17 - qJD(4) * t116, t366 + t1 * qJD(3) + t3 * qJD(4) + ((t118 * t67 + t119 * t93 + t120 * t66 + t121 * t92 + t45 * t72 + t73 * t87) * t405 + (t100 * t99 + t111 * t162 + t112 * t160 + t159 * t96 + t161 * t95 + t84 * t94) * t406) * t409 + (t278 * t251 + t277 * t250 + (t442 * t424 + (t440 * t250 - t441 * t251) * t396 + (t438 * t250 - t439 * t251 + t443) * t394) * t255 + ((-t90 / 0.2e1 - t88 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1) * t251 + (t91 / 0.2e1 + t89 / 0.2e1 + t27 / 0.2e1 + t26 / 0.2e1) * t250) * t253 + t425) * qJD(2), t1 * qJD(2) + (m(4) * (t100 * t105 + t111 * t113 + t112 * t114) / 0.4e1 + m(5) * (t73 * t83 + t92 * t97 + t93 * t98) / 0.4e1) * t407 + t31 * t390 + (-t230 / 0.2e1 - t231 / 0.2e1) * qJD(3) * t255 * t412 + (t419 * t396 + t418 * t394 + ((t448 * t252 + t447 * t254) * t255 + (t71 + t69) * t251 + (t70 + t68) * t250) * t421) * t321, -t320 + t3 * qJD(2) + t31 * t391 + (t125 * t356 + t171 * t221 + t172 * t219 - t110) * t390; qJD(2) * t122 + qJD(3) * t116, (t219 * t119 + t221 * t121 - t42 + (t255 * t72 + (t286 + t87) * t253) * t252) * t392 + t5 * qJD(3) + t417, t320 + t5 * qJD(2) + (t219 * t98 + t220 * t93 + t221 * t97 + t222 * t92 + (t252 * t83 + t254 * t73) * t253 - t31) * t391 + t110 * t390, 0.4e1 * (t115 * qJD(2) / 0.4e1 + t110 * qJD(3) / 0.4e1) * m(5);];
Cq = t7;
