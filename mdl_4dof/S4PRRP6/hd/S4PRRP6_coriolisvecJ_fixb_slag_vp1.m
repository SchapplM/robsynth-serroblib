% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:39
% DurationCPUTime: 23.52s
% Computational Cost: add. (7276->561), mult. (20357->821), div. (0->0), fcn. (19918->6), ass. (0->284)
t496 = Icges(5,4) + Icges(4,5);
t495 = Icges(4,6) - Icges(5,6);
t236 = sin(pkin(6));
t240 = cos(qJ(3));
t241 = cos(qJ(2));
t376 = t240 * t241;
t237 = cos(pkin(6));
t238 = sin(qJ(3));
t383 = t237 * t238;
t197 = t236 * t376 - t383;
t192 = Icges(5,5) * t197;
t379 = t238 * t241;
t196 = t236 * t379 + t237 * t240;
t239 = sin(qJ(2));
t385 = t236 * t239;
t83 = Icges(5,6) * t385 + Icges(5,3) * t196 + t192;
t402 = Icges(4,4) * t197;
t89 = -Icges(4,2) * t196 + Icges(4,6) * t385 + t402;
t487 = t83 - t89;
t386 = t236 * t238;
t199 = t237 * t376 + t386;
t193 = Icges(5,5) * t199;
t198 = -t236 * t240 + t237 * t379;
t382 = t237 * t239;
t84 = Icges(5,6) * t382 + Icges(5,3) * t198 + t193;
t401 = Icges(4,4) * t199;
t90 = -Icges(4,2) * t198 + Icges(4,6) * t382 + t401;
t486 = t84 - t90;
t85 = Icges(4,5) * t197 - Icges(4,6) * t196 + Icges(4,3) * t385;
t87 = Icges(5,4) * t197 + Icges(5,2) * t385 + Icges(5,6) * t196;
t485 = t87 + t85;
t86 = Icges(4,5) * t199 - Icges(4,6) * t198 + Icges(4,3) * t382;
t88 = Icges(5,4) * t199 + Icges(5,2) * t382 + Icges(5,6) * t198;
t484 = t88 + t86;
t396 = Icges(5,5) * t196;
t91 = Icges(5,1) * t197 + Icges(5,4) * t385 + t396;
t194 = Icges(4,4) * t196;
t93 = Icges(4,1) * t197 + Icges(4,5) * t385 - t194;
t508 = t91 + t93;
t395 = Icges(5,5) * t198;
t92 = Icges(5,1) * t199 + Icges(5,4) * t382 + t395;
t195 = Icges(4,4) * t198;
t94 = Icges(4,1) * t199 + Icges(4,5) * t382 - t195;
t507 = t92 + t94;
t473 = Icges(4,2) + Icges(5,3);
t526 = Icges(5,2) + Icges(4,3);
t377 = t239 * t240;
t231 = Icges(5,5) * t377;
t380 = t238 * t239;
t390 = Icges(5,6) * t241;
t164 = Icges(5,3) * t380 + t231 - t390;
t399 = Icges(4,4) * t240;
t303 = -Icges(4,2) * t238 + t399;
t170 = -Icges(4,6) * t241 + t239 * t303;
t504 = t164 - t170;
t394 = Icges(5,5) * t238;
t306 = Icges(5,1) * t240 + t394;
t172 = -Icges(5,4) * t241 + t239 * t306;
t400 = Icges(4,4) * t238;
t307 = Icges(4,1) * t240 - t400;
t174 = -Icges(4,5) * t241 + t239 * t307;
t471 = t172 + t174;
t525 = (-t238 * t496 - t240 * t495) * t239;
t493 = t196 * t486 + t197 * t507 + t385 * t484;
t492 = t198 * t487 + t199 * t508 + t382 * t485;
t494 = t196 * t487 + t197 * t508 + t385 * t485;
t458 = t198 * t486 + t199 * t507 + t382 * t484;
t300 = Icges(4,5) * t240 - Icges(4,6) * t238;
t166 = -Icges(4,3) * t241 + t239 * t300;
t302 = Icges(5,4) * t240 + Icges(5,6) * t238;
t168 = -Icges(5,2) * t241 + t239 * t302;
t378 = t239 * t168;
t524 = t166 * t385 + t196 * t504 + t197 * t471 + t236 * t378;
t523 = t166 * t382 + t198 * t504 + t199 * t471 + t237 * t378;
t355 = qJD(3) * t241;
t341 = t240 * t355;
t358 = qJD(2) * t239;
t128 = -t236 * t341 + (qJD(3) * t237 + t236 * t358) * t238;
t343 = t240 * t358;
t129 = -qJD(3) * t196 - t236 * t343;
t357 = qJD(2) * t241;
t345 = t236 * t357;
t522 = t128 * t495 + t129 * t496 + t345 * t526;
t130 = -qJD(3) * t386 - t237 * t341 + t358 * t383;
t131 = -qJD(3) * t198 - t237 * t343;
t344 = t237 * t357;
t521 = t130 * t495 + t131 * t496 + t344 * t526;
t520 = Icges(4,1) + Icges(5,1);
t519 = Icges(4,4) - Icges(5,5);
t167 = Icges(4,3) * t239 + t241 * t300;
t169 = Icges(5,2) * t239 + t241 * t302;
t518 = t525 * qJD(3) + (t167 + t169) * qJD(2);
t393 = Icges(5,5) * t240;
t299 = Icges(5,3) * t238 + t393;
t165 = Icges(5,6) * t239 + t241 * t299;
t171 = Icges(4,6) * t239 + t241 * t303;
t451 = -t165 + t171;
t517 = -t166 - t168;
t173 = Icges(5,4) * t239 + t241 * t306;
t175 = Icges(4,5) * t239 + t241 * t307;
t450 = -t173 - t175;
t516 = (-t240 * t473 + t394 - t400) * t239;
t515 = t493 * t237;
t514 = t492 * t236;
t234 = t236 ^ 2;
t235 = t237 ^ 2;
t437 = t234 + t235;
t513 = m(3) * qJD(2) ^ 2 * (rSges(3,1) * t239 + rSges(3,2) * t241) * t437;
t512 = -t128 * t473 - t129 * t519 - t345 * t495;
t511 = -t130 * t473 - t131 * t519 - t344 * t495;
t510 = t128 * t519 + t129 * t520 + t345 * t496;
t509 = t130 * t519 + t131 * t520 + t344 * t496;
t506 = qJD(2) * t451 + qJD(3) * t516;
t213 = (-Icges(4,1) * t238 - t399) * t239;
t356 = qJD(3) * t239;
t505 = -(-Icges(5,1) * t238 + t393) * t356 - qJD(3) * t213 + t450 * qJD(2);
t503 = -t239 * t518 + t357 * t517;
t502 = t239 * t521 + t357 * t484;
t501 = t239 * t522 + t357 * t485;
t500 = t237 * t458 + t514;
t499 = t236 * t494 + t515;
t498 = t523 * t239;
t497 = t524 * t239;
t360 = qJD(2) * t236;
t219 = t237 * t356 + t360;
t425 = t219 / 0.2e1;
t359 = qJD(2) * t237;
t220 = t236 * t356 - t359;
t423 = t220 / 0.2e1;
t464 = -t128 * t487 + t129 * t508 + t196 * t512 + t197 * t510 + t236 * t501;
t463 = -t128 * t486 + t129 * t507 + t196 * t511 + t197 * t509 + t236 * t502;
t462 = -t130 * t487 + t131 * t508 + t198 * t512 + t199 * t510 + t237 * t501;
t461 = -t130 * t486 + t131 * t507 + t198 * t511 + t199 * t509 + t237 * t502;
t313 = t238 * t83 + t240 * t91;
t47 = t239 * t313 - t241 * t87;
t311 = -t238 * t89 + t240 * t93;
t49 = t239 * t311 - t241 * t85;
t491 = t47 + t49;
t312 = t238 * t84 + t240 * t92;
t48 = t239 * t312 - t241 * t88;
t310 = -t238 * t90 + t240 * t94;
t50 = t239 * t310 - t241 * t86;
t490 = t48 + t50;
t297 = t164 * t238 + t172 * t240;
t387 = t168 * t241;
t63 = t239 * t297 - t387;
t296 = -t170 * t238 + t174 * t240;
t388 = t166 * t241;
t64 = t239 * t296 - t388;
t476 = t63 + t64;
t325 = rSges(4,1) * t240 - rSges(4,2) * t238;
t177 = -rSges(4,3) * t241 + t239 * t325;
t481 = t177 * t219;
t480 = t177 * t220;
t479 = (t130 * t504 - t131 * t471 + t198 * t506 + t199 * t505 + t237 * t503) * t241 + (t241 * t500 + t498) * qJD(2);
t478 = (t128 * t504 - t129 * t471 + t196 * t506 + t197 * t505 + t236 * t503) * t241 + (t241 * t499 + t497) * qJD(2);
t477 = t462 * t423 + t461 * t425 + t479 * qJD(3) / 0.2e1;
t470 = t525 * t355 + (t196 * t496 + t197 * t495) * t220 + (t198 * t496 + t199 * t495) * t219;
t315 = t236 * t49 + t237 * t50;
t316 = t236 * t47 + t237 * t48;
t469 = t315 + t316;
t468 = (t219 * t493 + t220 * t494 - t355 * t524) * t236 + (t219 * t458 + t220 * t492 - t355 * t523) * t237;
t467 = qJD(3) * t478 + t219 * t463 + t220 * t464;
t441 = ((t311 + t313) * qJD(2) - t522) * t241 + (t510 * t240 + t512 * t238 + (-t238 * t508 + t240 * t487) * qJD(3) + t485 * qJD(2)) * t239;
t440 = ((t310 + t312) * qJD(2) - t521) * t241 + (t509 * t240 + t511 * t238 + (-t238 * t507 + t240 * t486) * qJD(3) + t484 * qJD(2)) * t239;
t465 = rSges(5,1) + pkin(3);
t460 = t219 * t490 + t220 * t491 - t355 * t476;
t457 = rSges(5,3) + qJ(4);
t261 = -t239 * t299 + t390;
t132 = t261 * t236;
t138 = t170 * t236;
t455 = t132 + t138;
t133 = t261 * t237;
t139 = t170 * t237;
t454 = t133 + t139;
t140 = t172 * t236;
t142 = t174 * t236;
t453 = -t140 - t142;
t141 = t172 * t237;
t143 = t174 * t237;
t452 = -t141 - t143;
t283 = t167 - t296;
t284 = -t169 + t297;
t429 = -(-t166 * t237 - t310) * t219 - (-t166 * t236 - t311) * t220;
t430 = (t168 * t237 + t312) * t219 + (t168 * t236 + t313) * t220;
t449 = (-t429 - t430 + (-t283 + t284) * t355) * t239;
t323 = pkin(3) * t240 + qJ(4) * t238;
t324 = rSges(5,1) * t240 + rSges(5,3) * t238;
t448 = rSges(5,2) * t241 + (-t323 - t324) * t239;
t447 = t219 * t484 + t220 * t485;
t227 = pkin(2) * t239 - pkin(5) * t241;
t282 = qJD(2) * t227;
t190 = t236 * t282;
t191 = t237 * t282;
t446 = -t190 * t236 - t191 * t237 + t282 * t437;
t445 = -t387 - t388;
t403 = Icges(3,4) * t241;
t404 = Icges(3,4) * t239;
t444 = -(-Icges(3,2) * t241 - t404) * t358 + (-Icges(3,1) * t239 - t403) * t357;
t436 = t476 * t358 + (t518 * t241 + (t505 * t240 + t506 * t238 + (t238 * t471 - t240 * t504) * qJD(3)) * t239 + ((-t296 - t297) * t241 + t517 * t239 + t469) * qJD(2)) * t241;
t435 = (Icges(5,1) * t380 - t213 - t231 - t504) * t355 + (-t196 * t520 + t192 - t402 + t487) * t220 + (-t198 * t520 + t193 - t401 + t486) * t219;
t434 = (t471 + t516) * t355 + (t197 * t473 + t194 - t396 - t508) * t220 + (t199 * t473 + t195 - t395 - t507) * t219;
t433 = t470 * t239;
t342 = t240 * t356;
t432 = t238 * t357 + t342;
t352 = qJD(4) * t238;
t230 = t239 * t352;
t228 = pkin(2) * t241 + pkin(5) * t239;
t206 = t228 * t236;
t207 = t228 * t237;
t338 = t206 * t360 + t207 * t359 + qJD(1);
t405 = rSges(5,2) * t382 + t198 * t457 + t199 * t465;
t406 = rSges(5,2) * t385 + t196 * t457 + t197 * t465;
t33 = t219 * t406 - t220 * t405 + t230 + t338;
t353 = qJD(4) * t198;
t418 = rSges(5,2) * t344 - t130 * t457 + t131 * t465 + t353;
t370 = -t190 * t360 - t191 * t359;
t354 = qJD(4) * t196;
t419 = rSges(5,2) * t345 - t128 * t457 + t129 * t465 + t354;
t5 = qJD(4) * t342 - t418 * t220 + t419 * t219 + (t352 + (-t236 * t405 + t237 * t406) * qJD(3)) * t357 + t370;
t431 = t33 * t418 + t405 * t5;
t426 = -t219 / 0.2e1;
t424 = -t220 / 0.2e1;
t178 = rSges(5,2) * t239 + t241 * t324;
t215 = (-rSges(5,1) * t238 + rSges(5,3) * t240) * t239;
t407 = qJD(2) * t178 + qJD(3) * t215 + t230 + t432 * qJ(4) + (-t238 * t356 + t240 * t357) * pkin(3);
t384 = t236 * t241;
t381 = t237 * t241;
t179 = rSges(4,3) * t239 + t241 * t325;
t216 = (-rSges(4,1) * t238 - rSges(4,2) * t240) * t239;
t107 = qJD(2) * t179 + qJD(3) * t216;
t222 = qJD(2) * t228;
t375 = -t107 - t222;
t374 = -t196 * t465 + t197 * t457;
t373 = t198 * t465 - t199 * t457;
t372 = t448 * t236;
t371 = t448 * t237;
t364 = -t177 - t227;
t362 = t206 * t236 + t207 * t237;
t361 = (-pkin(3) * t238 + qJ(4) * t240) * t239 + t215;
t351 = qJD(2) * qJD(3);
t350 = -t222 - t407;
t348 = -t227 + t448;
t347 = t227 * t360;
t346 = t227 * t359;
t334 = -t355 / 0.2e1;
t333 = t355 / 0.2e1;
t332 = t351 / 0.2e1;
t327 = t239 * t332;
t326 = t241 * t332;
t96 = rSges(4,1) * t197 - rSges(4,2) * t196 + rSges(4,3) * t385;
t98 = rSges(4,1) * t199 - rSges(4,2) * t198 + rSges(4,3) * t382;
t314 = -t236 * t98 + t237 * t96;
t309 = Icges(3,1) * t241 - t404;
t305 = -Icges(3,2) * t239 + t403;
t301 = -Icges(3,5) * t239 - Icges(3,6) * t241;
t294 = t236 * t326;
t293 = t237 * t326;
t46 = t219 * t96 - t220 * t98 + t338;
t281 = t46 * t314;
t273 = qJD(2) * t301;
t272 = t33 * t419 + t406 * t5;
t36 = t219 * t448 - t355 * t405 - t347 + t354;
t37 = -t220 * t448 + t355 * t406 - t346 + t353;
t271 = t36 * t405 - t37 * t406;
t257 = (-(-Icges(3,6) * t237 + t236 * t305) * t241 - (-Icges(3,5) * t237 + t236 * t309) * t239) * qJD(2) + t444 * t236;
t256 = (-(Icges(3,6) * t236 + t237 * t305) * t241 - (Icges(3,5) * t236 + t237 * t309) * t239) * qJD(2) + t444 * t237;
t245 = (t33 * t406 + t36 * t448) * t237 + (-t33 * t405 - t37 * t448) * t236;
t201 = t301 * t237;
t200 = t301 * t236;
t183 = t237 * t273;
t182 = t236 * t273;
t126 = -rSges(4,1) * t198 - rSges(4,2) * t199;
t122 = -rSges(4,1) * t196 - rSges(4,2) * t197;
t80 = rSges(4,1) * t131 + rSges(4,2) * t130 + rSges(4,3) * t344;
t78 = rSges(4,1) * t129 + rSges(4,2) * t128 + rSges(4,3) * t345;
t60 = t355 * t96 - t346 + t480;
t59 = -t355 * t98 - t347 - t481;
t35 = t78 * t355 + t107 * t220 + (-t222 * t237 + (t177 * t384 - t239 * t96) * qJD(3)) * qJD(2);
t34 = -t80 * t355 - t107 * t219 + (-t222 * t236 + (-t177 * t381 + t239 * t98) * qJD(3)) * qJD(2);
t26 = t241 * t314 * t351 + t219 * t78 - t220 * t80 + t370;
t15 = -t222 * t359 - qJD(4) * t130 + t407 * t220 + (t419 * t241 + (-t239 * t406 - t384 * t448) * qJD(2)) * qJD(3);
t14 = -t222 * t360 - qJD(4) * t128 - t407 * t219 + (-t418 * t241 + (t239 * t405 + t381 * t448) * qJD(2)) * qJD(3);
t1 = [m(4) * t26 + m(5) * t5 - t513; (t200 * qJD(2) * t235 - t201 * t237 * t360) * t359 / 0.2e1 - (t236 * (-t183 * t237 + t236 * t256) - t237 * (-t182 * t237 + t236 * t257)) * t359 - (qJD(2) * t201 * t234 - t200 * t236 * t359) * t360 / 0.2e1 + (t236 * (t183 * t236 + t237 * t256) - t237 * (t182 * t236 + t237 * t257)) * t360 + (((t198 * t451 + t199 * t450 + t514) * t241 + t498) * qJD(3) + (((t445 + t458) * qJD(3) + t447) * t241 + t449) * t237 + (t198 * t455 + t199 * t453) * t220 + (t198 * t454 + t199 * t452) * t219) * t426 + (t236 * t461 - t237 * t462) * t425 + (((t196 * t451 + t197 * t450 + t515) * t241 + t497) * qJD(3) + (((t445 + t494) * qJD(3) + t447) * t241 + t449) * t236 + (t196 * t455 + t197 * t453) * t220 + (t196 * t454 + t197 * t452) * t219) * t424 + (t236 * t463 - t237 * t464) * t423 + t236 * t477 - t467 * t237 / 0.2e1 - t460 * t356 / 0.2e1 + (((t139 * t238 - t143 * t240 + t86) * t219 + (t138 * t238 - t142 * t240 + t85) * t220 + t64 * qJD(3)) * t239 + ((t283 * t241 + (t171 * t238 - t175 * t240 - t166) * t239 + t315) * qJD(3) + t429) * t241 + ((t133 * t238 - t141 * t240 + t88) * t219 + (t132 * t238 - t140 * t240 + t87) * t220 + t63 * qJD(3)) * t239 + ((-t284 * t241 + (-t165 * t238 - t173 * t240 - t168) * t239 + t316) * qJD(3) + t430) * t241) * t333 + (t236 * t490 - t237 * t491) * t327 + (t236 * t493 - t237 * t494) * t294 + (t458 * t236 - t237 * t492) * t293 + (t5 * t362 + (t15 * t348 + t350 * t37 + t431) * t237 + (t14 * t348 + t350 * t36 + t272) * t236 - (t271 * t239 + (-t36 * t371 + t37 * t372 + t245) * t241) * qJD(3) - (t236 * t36 + t237 * t37) * (-t230 - t222) + (t219 * t36 - t220 * t37) * (t241 * t323 + t178) + (-t219 * t372 + t220 * t371 - t241 * t352 + t446) * t33) * m(5) + (t26 * t362 - t60 * (t179 * t220 - t228 * t359) - t59 * (-t179 * t219 - t228 * t360) - ((t59 * t98 - t60 * t96) * t239 + t281 * t241) * qJD(3) + t446 * t46 + (t26 * t98 + t35 * t364 + t375 * t60 + (-t480 + t80) * t46) * t237 + (t26 * t96 + t34 * t364 + t375 * t59 + (t481 + t78) * t46) * t236) * m(4) + (t236 * t440 - t237 * t441 + t468) * t334 + (-t437 + 0.1e1) * (rSges(3,1) * t241 - rSges(3,2) * t239) * t513; (t198 * t434 + t199 * t435 - t237 * t433) * t426 + ((t236 * t462 + t237 * t461) * t239 + t479) * t425 + (t196 * t434 + t197 * t435 - t236 * t433) * t424 + ((t236 * t464 + t237 * t463) * t239 + t478) * t423 - (qJD(3) * t436 + t219 * t440 + t220 * t441) * t241 / 0.2e1 + t467 * t385 / 0.2e1 + t382 * t477 + t460 * t358 / 0.2e1 + ((t236 * t441 + t237 * t440) * t239 + t436) * t334 + (t470 * t241 + (t238 * t434 + t240 * t435) * t239) * t333 + (t239 * t469 - t241 * t476) * t327 + (t239 * t499 - t241 * t524) * t294 + (t239 * t500 - t241 * t523) * t293 + (-(t197 * t36 + t199 * t37 + t33 * t377) * qJD(4) - (t33 * t373 + t361 * t37) * t220 - (t33 * t374 - t36 * t361) * t219 - (t36 * t373 + t37 * t374) * t355 + (qJD(2) * t245 - t14 * t405 + t15 * t406 - t36 * t418 + t37 * t419) * t241 + (t271 * qJD(2) + (t14 * t448 - t36 * t407 + t272) * t237 + (-t15 * t448 + t37 * t407 - t431) * t236) * t239) * m(5) + (-t60 * (t122 * t355 + t216 * t220) - t59 * (-t126 * t355 - t216 * t219) - t46 * (t122 * t219 - t126 * t220) + (-t34 * t98 + t35 * t96 - t59 * t80 + t60 * t78 + (t281 + (t236 * t60 - t237 * t59) * t177) * qJD(2)) * t241 + (t60 * (-qJD(2) * t96 + t107 * t236) + t59 * (qJD(2) * t98 - t107 * t237) + t26 * t314 + t46 * (-t236 * t80 + t237 * t78) + (t236 * t35 - t237 * t34) * t177) * t239) * m(4) + t468 * t357 / 0.2e1; (t5 * t380 + t14 * t196 + t15 * t198 + (-t196 * t355 - t220 * t380 - t130) * t37 + (t198 * t355 + t219 * t380 - t128) * t36 + (-t196 * t219 + t198 * t220 + t432) * t33) * m(5);];
tauc = t1(:);
