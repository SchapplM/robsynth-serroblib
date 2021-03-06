% Calculate vector of inverse dynamics joint torques for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:11
% EndTime: 2019-12-05 15:14:47
% DurationCPUTime: 29.11s
% Computational Cost: add. (23796->872), mult. (26351->1272), div. (0->0), fcn. (25456->8), ass. (0->416)
t332 = cos(pkin(8));
t564 = t332 ^ 2;
t331 = sin(pkin(8));
t565 = t331 ^ 2;
t577 = t564 + t565;
t328 = pkin(9) + qJ(3);
t324 = sin(t328);
t325 = cos(t328);
t405 = -Icges(4,5) * t324 - Icges(4,6) * t325;
t584 = t331 * t332;
t330 = qJ(4) + qJ(5);
t326 = sin(t330);
t537 = rSges(6,2) * t326;
t327 = cos(t330);
t539 = rSges(6,1) * t327;
t423 = -t537 + t539;
t197 = t324 * rSges(6,3) + t325 * t423;
t333 = sin(qJ(4));
t334 = cos(qJ(4));
t540 = rSges(5,1) * t334;
t424 = -rSges(5,2) * t333 + t540;
t211 = t324 * rSges(5,3) + t325 * t424;
t583 = 2 * qJD(3);
t582 = 2 * qJDD(3);
t469 = qJD(3) * qJD(4);
t579 = qJDD(4) * t324 + t325 * t469;
t319 = pkin(4) * t334 + pkin(3);
t542 = pkin(3) - t319;
t440 = t542 * t324;
t335 = -pkin(7) - pkin(6);
t541 = pkin(6) + t335;
t578 = -t325 * t541 + t440;
t497 = t332 * t333;
t499 = t331 * t334;
t278 = -t325 * t497 + t499;
t475 = qJD(3) * t324;
t474 = qJD(3) * t325;
t570 = t325 * pkin(3) + t324 * pkin(6);
t569 = g(1) * t332 + g(2) * t331;
t403 = Icges(6,5) * t327 - Icges(6,6) * t326;
t190 = -Icges(6,3) * t325 + t324 * t403;
t404 = Icges(5,5) * t334 - Icges(5,6) * t333;
t202 = -Icges(5,3) * t325 + t324 * t404;
t519 = Icges(6,4) * t327;
t407 = -Icges(6,2) * t326 + t519;
t192 = -Icges(6,6) * t325 + t324 * t407;
t523 = Icges(5,4) * t334;
t408 = -Icges(5,2) * t333 + t523;
t204 = -Icges(5,6) * t325 + t324 * t408;
t520 = Icges(6,4) * t326;
t411 = Icges(6,1) * t327 - t520;
t194 = -Icges(6,5) * t325 + t324 * t411;
t524 = Icges(5,4) * t333;
t412 = Icges(5,1) * t334 - t524;
t206 = -Icges(5,5) * t325 + t324 * t412;
t536 = pkin(4) * qJD(4);
t459 = t333 * t536;
t339 = qJD(3) * t578 - t325 * t459;
t458 = t334 * t536;
t102 = t331 * t339 - t332 * t458;
t103 = t331 * t458 + t332 * t339;
t503 = t327 * t332;
t505 = t326 * t331;
t252 = -t325 * t505 - t503;
t498 = t332 * t326;
t504 = t327 * t331;
t253 = t325 * t504 - t498;
t511 = t324 * t331;
t121 = rSges(6,1) * t253 + rSges(6,2) * t252 + rSges(6,3) * t511;
t254 = -t325 * t498 + t504;
t255 = t325 * t503 + t505;
t510 = t324 * t332;
t122 = rSges(6,1) * t255 + rSges(6,2) * t254 + rSges(6,3) * t510;
t359 = -t324 * t541 - t325 * t542;
t126 = -pkin(4) * t497 + t331 * t359;
t500 = t331 * t333;
t127 = pkin(4) * t500 + t332 * t359;
t320 = qJDD(3) * t331;
t198 = t332 * t579 + t320;
t468 = qJD(3) * qJD(5);
t363 = qJDD(5) * t324 + t325 * t468;
t139 = t332 * t363 + t198;
t466 = qJDD(3) * t332;
t199 = t331 * t579 - t466;
t140 = t331 * t363 + t199;
t322 = qJD(3) * t331;
t472 = qJD(4) * t324;
t283 = t332 * t472 + t322;
t470 = qJD(5) * t324;
t212 = t332 * t470 + t283;
t473 = qJD(3) * t332;
t284 = t331 * t472 - t473;
t213 = t331 * t470 + t284;
t291 = t324 * pkin(3) - t325 * pkin(6);
t376 = qJD(3) * t291;
t250 = t331 * t376;
t251 = t332 * t376;
t270 = t570 * t331;
t271 = t570 * t332;
t389 = -t250 * t322 - t251 * t473 + t270 * t320 + t271 * t466 + qJDD(1);
t329 = qJD(4) + qJD(5);
t501 = t329 * t332;
t372 = t322 * t324 + t501;
t502 = t329 * t331;
t455 = t325 * t502;
t153 = t326 * t372 - t327 * t455;
t154 = -t326 * t455 - t327 * t372;
t446 = t325 * t322;
t86 = rSges(6,1) * t154 + rSges(6,2) * t153 + rSges(6,3) * t446;
t373 = t324 * t473 - t502;
t454 = t325 * t501;
t155 = t326 * t373 - t327 * t454;
t156 = -t326 * t454 - t327 * t373;
t445 = t325 * t473;
t87 = rSges(6,1) * t156 + rSges(6,2) * t155 + rSges(6,3) * t445;
t13 = t102 * t283 - t103 * t284 + t121 * t139 - t122 * t140 + t126 * t198 - t127 * t199 + t212 * t86 - t213 * t87 + t389;
t442 = t270 * t322 + t271 * t473 + qJD(1);
t48 = t121 * t212 - t122 * t213 + t126 * t283 - t127 * t284 + t442;
t491 = t122 + t127;
t530 = t103 + t87;
t568 = -t13 * t491 - t48 * t530;
t496 = t332 * t334;
t279 = t325 * t496 + t500;
t525 = Icges(5,4) * t279;
t132 = Icges(5,2) * t278 + Icges(5,6) * t510 + t525;
t265 = Icges(5,4) * t278;
t134 = Icges(5,1) * t279 + Icges(5,5) * t510 + t265;
t398 = -t132 * t333 + t134 * t334;
t276 = -t325 * t500 - t496;
t277 = t325 * t499 - t497;
t526 = Icges(5,4) * t277;
t131 = Icges(5,2) * t276 + Icges(5,6) * t511 + t526;
t264 = Icges(5,4) * t276;
t133 = Icges(5,1) * t277 + Icges(5,5) * t511 + t264;
t399 = -t131 * t333 + t133 * t334;
t567 = -(-t202 * t332 - t398) * t283 - (-t202 * t331 - t399) * t284;
t191 = Icges(6,3) * t324 + t325 * t403;
t396 = -t192 * t326 + t194 * t327;
t521 = Icges(6,4) * t255;
t118 = Icges(6,2) * t254 + Icges(6,6) * t510 + t521;
t231 = Icges(6,4) * t254;
t120 = Icges(6,1) * t255 + Icges(6,5) * t510 + t231;
t401 = -t118 * t326 + t120 * t327;
t522 = Icges(6,4) * t253;
t117 = Icges(6,2) * t252 + Icges(6,6) * t511 + t522;
t230 = Icges(6,4) * t252;
t119 = Icges(6,1) * t253 + Icges(6,5) * t511 + t230;
t402 = -t117 * t326 + t119 * t327;
t429 = t325 * t329;
t566 = t212 * (-t190 * t332 - t401) + t213 * (-t190 * t331 - t402) - t429 * (t191 - t396);
t237 = (-Icges(6,2) * t327 - t520) * t324;
t344 = t212 * (-Icges(6,2) * t255 + t120 + t231) + t213 * (-Icges(6,2) * t253 + t119 + t230) - t429 * (t194 + t237);
t273 = (-Icges(5,2) * t334 - t524) * t324;
t471 = qJD(4) * t325;
t340 = t283 * (-Icges(5,2) * t279 + t134 + t265) + t284 * (-Icges(5,2) * t277 + t133 + t264) - t471 * (t206 + t273);
t563 = t139 / 0.2e1;
t562 = t140 / 0.2e1;
t311 = t324 * t469;
t189 = t324 * t468 + t311 + (-qJDD(4) - qJDD(5)) * t325;
t561 = t189 / 0.2e1;
t560 = t198 / 0.2e1;
t559 = t199 / 0.2e1;
t558 = -t212 / 0.2e1;
t557 = t212 / 0.2e1;
t556 = -t213 / 0.2e1;
t555 = t213 / 0.2e1;
t282 = -qJDD(4) * t325 + t311;
t554 = t282 / 0.2e1;
t553 = -t283 / 0.2e1;
t552 = t283 / 0.2e1;
t551 = -t284 / 0.2e1;
t550 = t284 / 0.2e1;
t549 = t429 / 0.2e1;
t548 = -t429 / 0.2e1;
t547 = -t325 / 0.2e1;
t129 = Icges(5,5) * t277 + Icges(5,6) * t276 + Icges(5,3) * t511;
t61 = t129 * t510 + t131 * t278 + t133 * t279;
t535 = t331 * t61;
t130 = Icges(5,5) * t279 + Icges(5,6) * t278 + Icges(5,3) * t510;
t60 = t130 * t511 + t132 * t276 + t134 * t277;
t534 = t332 * t60;
t196 = -rSges(6,3) * t325 + t324 * t423;
t476 = qJD(2) * t332;
t370 = -t291 * t322 - t476;
t52 = -t122 * t429 - t127 * t471 - t196 * t212 + t283 * t578 + t370;
t533 = t52 * t127;
t72 = t202 * t511 + t204 * t276 + t206 * t277;
t532 = t72 * t324;
t73 = t202 * t510 + t204 * t278 + t206 * t279;
t531 = t73 * t324;
t529 = t121 * t445 + t86 * t510;
t512 = t202 * t325;
t509 = t324 * t333;
t508 = t325 * t331;
t507 = t325 * t332;
t506 = t325 * t335;
t422 = -rSges(6,1) * t326 - rSges(6,2) * t327;
t247 = t422 * t324;
t113 = qJD(3) * t197 + t329 * t247;
t138 = qJD(3) * t359 - t324 * t459;
t495 = -t113 - t138;
t492 = t121 + t126;
t275 = (-rSges(5,1) * t333 - rSges(5,2) * t334) * t324;
t128 = qJD(3) * t211 + qJD(4) * t275;
t281 = t570 * qJD(3);
t490 = -t128 - t281;
t487 = t578 - t196;
t484 = -t331 * t250 - t332 * t251;
t210 = -rSges(5,3) * t325 + t324 * t424;
t483 = -t210 - t291;
t308 = pkin(6) * t508;
t309 = pkin(6) * t507;
t482 = (-pkin(3) * t511 + t308) * t322 + (-pkin(3) * t510 + t309) * t473;
t147 = t252 * rSges(6,1) - t253 * rSges(6,2);
t148 = t254 * rSges(6,1) - t255 * rSges(6,2);
t481 = t331 * t270 + t332 * t271;
t460 = t324 * t537;
t480 = rSges(6,3) * t508 + t331 * t460;
t479 = rSges(6,3) * t507 + t332 * t460;
t461 = rSges(5,2) * t509;
t478 = rSges(5,3) * t508 + t331 * t461;
t477 = rSges(5,3) * t507 + t332 * t461;
t467 = qJDD(2) * t332;
t464 = pkin(4) * t509;
t463 = t324 * t540;
t462 = t324 * t539;
t457 = t113 * t511 + t196 * t446 + t325 * t86;
t456 = -m(3) - m(4) - m(5) - m(6);
t452 = -t281 + t495;
t451 = -t291 + t487;
t450 = t570 * t322;
t449 = t570 * t473;
t448 = t333 * t475;
t447 = t334 * t475;
t444 = t511 / 0.2e1;
t443 = t510 / 0.2e1;
t438 = t475 / 0.2e1;
t437 = t473 / 0.2e1;
t436 = -t471 / 0.2e1;
t435 = t471 / 0.2e1;
t321 = qJDD(2) * t331;
t387 = -qJD(3) * t281 - qJDD(3) * t291;
t357 = t332 * t387 + t321;
t15 = t102 * t471 + t113 * t213 - t121 * t189 - t126 * t282 + t138 * t284 + t140 * t196 - t199 * t578 + t429 * t86 + t357;
t434 = t15 * (t325 * t121 + t196 * t511);
t432 = t212 * t147 - t148 * t213;
t431 = -t148 * t429 - t212 * t247;
t430 = t147 * t429 + t213 * t247;
t428 = t446 / 0.2e1;
t427 = t325 * t437;
t323 = qJD(2) * t331;
t425 = -t291 * t473 + t323;
t290 = rSges(4,1) * t325 - rSges(4,2) * t324;
t289 = rSges(4,1) * t324 + rSges(4,2) * t325;
t135 = rSges(5,1) * t277 + rSges(5,2) * t276 + rSges(5,3) * t511;
t183 = -qJD(4) * t277 + t331 * t448;
t184 = qJD(4) * t276 - t331 * t447;
t99 = rSges(5,1) * t184 + rSges(5,2) * t183 + rSges(5,3) * t446;
t49 = t128 * t284 - t135 * t282 + t199 * t210 + t471 * t99 + t357;
t185 = -qJD(4) * t279 + t332 * t448;
t186 = qJD(4) * t278 - t332 * t447;
t100 = rSges(5,1) * t186 + rSges(5,2) * t185 + rSges(5,3) * t445;
t136 = rSges(5,1) * t279 + rSges(5,2) * t278 + rSges(5,3) * t510;
t347 = t331 * t387 - t467;
t50 = -t100 * t471 - t128 * t283 + t136 * t282 - t198 * t210 + t347;
t421 = t331 * t49 - t332 * t50;
t115 = Icges(6,5) * t253 + Icges(6,6) * t252 + Icges(6,3) * t511;
t53 = t115 * t511 + t117 * t252 + t119 * t253;
t116 = Icges(6,5) * t255 + Icges(6,6) * t254 + Icges(6,3) * t510;
t54 = t116 * t511 + t118 * t252 + t120 * t253;
t420 = t331 * t53 + t332 * t54;
t55 = t115 * t510 + t117 * t254 + t119 * t255;
t56 = t116 * t510 + t118 * t254 + t120 * t255;
t419 = t331 * t55 + t332 * t56;
t57 = -t115 * t325 + t324 * t402;
t58 = -t116 * t325 + t324 * t401;
t418 = t331 * t57 + t332 * t58;
t59 = t129 * t511 + t131 * t276 + t133 * t277;
t417 = t331 * t59 + t534;
t62 = t130 * t510 + t132 * t278 + t134 * t279;
t416 = t332 * t62 + t535;
t64 = -t129 * t325 + t324 * t399;
t65 = -t130 * t325 + t324 * t398;
t415 = t64 * t331 + t65 * t332;
t406 = Icges(4,5) * t325 - Icges(4,6) * t324;
t400 = t129 * t284 + t130 * t283;
t397 = t135 * t332 - t136 * t331;
t395 = -t204 * t333 + t206 * t334;
t392 = -(-t289 * t473 + t323) * t332 - (-t289 * t322 - t476) * t331;
t391 = t577 * t290;
t390 = t577 * qJD(3) * t289;
t267 = t278 * pkin(4);
t280 = t290 * qJD(3);
t388 = -qJD(3) * t280 - qJDD(3) * t289;
t80 = Icges(6,5) * t154 + Icges(6,6) * t153 + Icges(6,3) * t446;
t386 = t115 * t474 + t324 * t80;
t81 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t445;
t385 = t116 * t474 + t324 * t81;
t93 = Icges(5,5) * t184 + Icges(5,6) * t183 + Icges(5,3) * t446;
t384 = t129 * t474 + t324 * t93;
t94 = Icges(5,5) * t186 + Icges(5,6) * t185 + Icges(5,3) * t445;
t383 = t130 * t474 + t324 * t94;
t203 = Icges(5,3) * t324 + t325 * t404;
t377 = t203 - t395;
t236 = (-Icges(6,5) * t326 - Icges(6,6) * t327) * t324;
t108 = qJD(3) * t191 + t236 * t329;
t375 = t108 * t324 + t190 * t474;
t272 = (-Icges(5,5) * t333 - Icges(5,6) * t334) * t324;
t123 = qJD(3) * t203 + qJD(4) * t272;
t374 = t123 * t324 + t202 * t474;
t371 = t440 - t506;
t63 = t135 * t283 - t136 * t284 + t442;
t368 = t63 * t397;
t274 = (-Icges(5,1) * t333 - t523) * t324;
t238 = (-Icges(6,1) * t326 - t519) * t324;
t360 = qJD(3) * t405;
t266 = t276 * pkin(4);
t358 = (Icges(6,5) * t252 - Icges(6,6) * t253) * t213 + (Icges(6,5) * t254 - Icges(6,6) * t255) * t212 - t236 * t429;
t356 = -(Icges(5,5) * t276 - Icges(5,6) * t277) * t284 - (Icges(5,5) * t278 - Icges(5,6) * t279) * t283 + t272 * t471;
t355 = t324 * t358;
t207 = Icges(5,5) * t324 + t325 * t412;
t195 = Icges(6,5) * t324 + t325 * t411;
t205 = Icges(5,6) * t324 + t325 * t408;
t193 = Icges(6,6) * t324 + t325 * t407;
t348 = t324 * t356;
t82 = Icges(6,4) * t154 + Icges(6,2) * t153 + Icges(6,6) * t446;
t84 = Icges(6,1) * t154 + Icges(6,4) * t153 + Icges(6,5) * t446;
t19 = (qJD(3) * t402 - t80) * t325 + (qJD(3) * t115 + (-t117 * t329 + t84) * t327 + (-t119 * t329 - t82) * t326) * t324;
t83 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t445;
t85 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t445;
t20 = (qJD(3) * t401 - t81) * t325 + (qJD(3) * t116 + (-t118 * t329 + t85) * t327 + (-t120 * t329 - t83) * t326) * t324;
t70 = t190 * t511 + t192 * t252 + t194 * t253;
t23 = t212 * t54 + t213 * t53 - t429 * t70;
t71 = t190 * t510 + t192 * t254 + t194 * t255;
t24 = t212 * t56 + t213 * t55 - t429 * t71;
t25 = t117 * t153 + t119 * t154 + t252 * t82 + t253 * t84 + t331 * t386;
t26 = t118 * t153 + t120 * t154 + t252 * t83 + t253 * t85 + t331 * t385;
t27 = t117 * t155 + t119 * t156 + t254 * t82 + t255 * t84 + t332 * t386;
t28 = t118 * t155 + t120 * t156 + t254 * t83 + t255 * t85 + t332 * t385;
t77 = -t190 * t325 + t324 * t396;
t30 = t212 * t58 + t213 * t57 - t429 * t77;
t345 = (Icges(6,1) * t254 - t118 - t521) * t212 + (Icges(6,1) * t252 - t117 - t522) * t213 - (-t192 + t238) * t429;
t109 = qJD(3) * t193 + t237 * t329;
t110 = qJD(3) * t195 + t238 * t329;
t42 = t109 * t252 + t110 * t253 + t153 * t192 + t154 * t194 + t331 * t375;
t4 = t139 * t54 + t140 * t53 + t189 * t70 + t212 * t26 + t213 * t25 - t42 * t429;
t43 = t109 * t254 + t110 * t255 + t155 * t192 + t156 * t194 + t332 * t375;
t44 = (qJD(3) * t396 - t108) * t325 + (qJD(3) * t190 + (-t192 * t329 + t110) * t327 + (-t194 * t329 - t109) * t326) * t324;
t5 = t139 * t56 + t140 * t55 + t189 * t71 + t212 * t28 + t213 * t27 - t429 * t43;
t346 = t4 * t444 + t5 * t443 + (t254 * t344 + t255 * t345 + t332 * t355) * t558 + (t252 * t344 + t253 * t345 + t331 * t355) * t556 + (-t358 * t325 + (-t326 * t344 + t345 * t327) * t324) * t549 + t23 * t428 + t24 * t427 + (t324 * t419 - t325 * t71) * t563 + (t324 * t420 - t325 * t70) * t562 + t30 * t438 + (t139 * t58 + t140 * t57 + t189 * t77 + t19 * t213 + t20 * t212 - t429 * t44) * t547 + (t324 * t418 - t325 * t77) * t561 + (-t325 * t43 + (t27 * t331 + t28 * t332) * t324 + (t324 * t71 + t325 * t419) * qJD(3)) * t557 + (-t325 * t42 + (t25 * t331 + t26 * t332) * t324 + (t324 * t70 + t325 * t420) * qJD(3)) * t555 + (-t325 * t44 + (t19 * t331 + t20 * t332) * t324 + (t324 * t77 + t325 * t418) * qJD(3)) * t548;
t341 = (Icges(5,1) * t278 - t132 - t525) * t283 + (Icges(5,1) * t276 - t131 - t526) * t284 - (-t204 + t274) * t471;
t338 = t577 * t405;
t337 = (-t377 * t471 - t567) * t324;
t336 = (t115 * t213 + t116 * t212 - t190 * t429) * t325 + t566 * t324;
t297 = t325 * t319;
t285 = t324 * t329;
t269 = t289 * t332;
t268 = t289 * t331;
t259 = t405 * t332;
t258 = t405 * t331;
t257 = t332 * t429;
t256 = t331 * t429;
t242 = t332 * t360;
t241 = t331 * t360;
t215 = Icges(4,3) * t331 + t332 * t406;
t214 = -Icges(4,3) * t332 + t331 * t406;
t188 = -t324 * t335 + t297 - t570;
t182 = -t332 * t463 + t477;
t181 = -t331 * t463 + t478;
t180 = t206 * t332;
t179 = t206 * t331;
t178 = t204 * t332;
t177 = t204 * t331;
t173 = -t332 * t462 + t479;
t172 = -t331 * t462 + t480;
t171 = t194 * t332;
t170 = t194 * t331;
t169 = t192 * t332;
t168 = t192 * t331;
t164 = rSges(5,1) * t278 - rSges(5,2) * t279;
t163 = rSges(5,1) * t276 - rSges(5,2) * t277;
t152 = t332 * t371 - t309;
t151 = t331 * t371 - t308;
t150 = t331 * t388 - t467;
t149 = t332 * t388 + t321;
t125 = qJD(3) * t207 + qJD(4) * t274;
t124 = qJD(3) * t205 + qJD(4) * t273;
t114 = qJD(3) * t391 + qJD(1);
t107 = t122 * t475;
t106 = t121 * t510;
t98 = Icges(5,1) * t186 + Icges(5,4) * t185 + Icges(5,5) * t445;
t97 = Icges(5,1) * t184 + Icges(5,4) * t183 + Icges(5,5) * t446;
t96 = Icges(5,4) * t186 + Icges(5,2) * t185 + Icges(5,6) * t445;
t95 = Icges(5,4) * t184 + Icges(5,2) * t183 + Icges(5,6) * t446;
t88 = t324 * t395 - t512;
t76 = -qJD(3) * t390 + qJDD(3) * t391 + qJDD(1);
t75 = -t136 * t471 - t210 * t283 + t370;
t74 = t135 * t471 + t210 * t284 + t425;
t51 = t121 * t429 + t126 * t471 + t196 * t213 - t284 * t578 + t425;
t47 = (qJD(3) * t395 - t123) * t325 + (qJD(3) * t202 - t124 * t333 + t125 * t334 + (-t204 * t334 - t206 * t333) * qJD(4)) * t324;
t46 = t124 * t278 + t125 * t279 + t185 * t204 + t186 * t206 + t332 * t374;
t45 = t124 * t276 + t125 * t277 + t183 * t204 + t184 * t206 + t331 * t374;
t41 = -t100 * t284 + t135 * t198 - t136 * t199 + t283 * t99 + t389;
t40 = t283 * t65 + t284 * t64 - t471 * t88;
t39 = t132 * t185 + t134 * t186 + t278 * t96 + t279 * t98 + t332 * t383;
t38 = t131 * t185 + t133 * t186 + t278 * t95 + t279 * t97 + t332 * t384;
t37 = t132 * t183 + t134 * t184 + t276 * t96 + t277 * t98 + t331 * t383;
t36 = t131 * t183 + t133 * t184 + t276 * t95 + t277 * t97 + t331 * t384;
t35 = t283 * t62 + t284 * t61 - t471 * t73;
t34 = t283 * t60 + t284 * t59 - t471 * t72;
t33 = (qJD(3) * t398 - t94) * t325 + (qJD(3) * t130 - t333 * t96 + t334 * t98 + (-t132 * t334 - t134 * t333) * qJD(4)) * t324;
t32 = (qJD(3) * t399 - t93) * t325 + (qJD(3) * t129 - t333 * t95 + t334 * t97 + (-t131 * t334 - t133 * t333) * qJD(4)) * t324;
t16 = -t103 * t471 - t113 * t212 + t122 * t189 + t127 * t282 - t138 * t283 - t139 * t196 + t198 * t578 - t429 * t87 + t347;
t10 = t198 * t62 + t199 * t61 + t282 * t73 + t283 * t39 + t284 * t38 - t46 * t471;
t9 = t198 * t60 + t199 * t59 + t282 * t72 + t283 * t37 + t284 * t36 - t45 * t471;
t1 = [(m(2) + m(3)) * qJDD(1) + m(4) * t76 + m(5) * t41 + m(6) * t13 + (-m(2) + t456) * g(3); t456 * (g(1) * t331 - g(2) * t332) + m(4) * (t149 * t331 - t150 * t332) + m(5) * t421 + m(6) * (t15 * t331 - t16 * t332) + m(3) * t577 * qJDD(2); ((-t169 * t254 - t171 * t255) * t212 + t56 * t257 + (-t168 * t254 - t170 * t255) * t213 + t55 * t256 - (t193 * t254 + t195 * t255) * t429 + t71 * t285 + t336 * t332) * t558 + (t57 * t256 + t58 * t257 + t77 * t285 - t566 * t325 + ((t169 * t326 - t171 * t327 + t116) * t212 + (t168 * t326 - t170 * t327 + t115) * t213 - (-t193 * t326 + t195 * t327 + t190) * t429) * t324) * t549 + (((t178 * t333 - t180 * t334 + t130) * t283 + (t177 * t333 - t179 * t334 + t129) * t284 + t88 * qJD(4)) * t324 + ((t377 * t325 + (t205 * t333 - t207 * t334 - t202) * t324 + t415) * qJD(4) + t567) * t325) * t435 + ((-t178 * t276 - t180 * t277) * t283 + (-t177 * t276 - t179 * t277) * t284 + (t532 + (-t205 * t276 - t207 * t277 + t534) * t325) * qJD(4) + (((t59 - t512) * qJD(4) + t400) * t325 + t337) * t331) * t551 + ((-t178 * t278 - t180 * t279) * t283 + (-t177 * t278 - t179 * t279) * t284 + (t531 + (-t205 * t278 - t207 * t279 + t535) * t325) * qJD(4) + (((t62 - t512) * qJD(4) + t400) * t325 + t337) * t332) * t553 + (t76 * t391 - t114 * t390 + (-t149 * t332 - t150 * t331) * t289 + t392 * t280 + g(1) * t269 + g(2) * t268 - g(3) * t290 - (t114 * (-t268 * t331 - t269 * t332) + t392 * t290) * qJD(3)) * m(4) + (-t51 * (-t121 * t285 + t172 * t429 + t188 * t284 + t196 * t256 + t197 * t213 - t449) - t52 * (t122 * t285 - t173 * t429 - t188 * t283 - t196 * t257 - t197 * t212 - t450) - t48 * (t121 * t257 - t122 * t256 + t151 * t283 - t152 * t284 + t172 * t212 - t173 * t213 + t482) - ((-t126 * t51 + t533) * t324 + (t51 * (-t331 * t578 + t151) + t52 * (t332 * t578 - t152) + t48 * (t126 * t332 - t127 * t331)) * t325) * qJD(4) - g(1) * (-t332 * t506 + t479) - g(2) * (-t331 * t506 + t480) - g(3) * (t197 + t297) - (-g(3) * t335 + t569 * (-t319 - t539)) * t324 + t13 * t481 + t48 * t484 + (t15 * t451 + t452 * t51 - t568) * t332 + (t16 * t451 + t52 * t452 + t13 * t492 + t48 * (t102 + t86)) * t331) * m(6) + ((-t32 + t35) * t332 + (t33 + t34) * t331) * t436 - t285 * t30 / 0.2e1 + ((-t169 * t252 - t171 * t253) * t212 + t54 * t257 + (-t168 * t252 - t170 * t253) * t213 + t53 * t256 - (t193 * t252 + t195 * t253) * t429 + t70 * t285 + t336 * t331) * t556 - ((t214 * t564 - t215 * t584) * t582 + (t241 * t564 - t242 * t584) * t583 + t9 + t4) * t332 / 0.2e1 + ((-t214 * t584 + t215 * t565) * t582 + (-t241 * t584 + t242 * t565) * t583 + t10 + t5) * t331 / 0.2e1 - t40 * t472 / 0.2e1 - t256 * t23 / 0.2e1 - t257 * t24 / 0.2e1 + (t331 * t65 - t332 * t64) * t554 + (-t25 * t332 + t26 * t331) * t555 + (-t27 * t332 + t28 * t331) * t557 + (t331 * t60 - t332 * t59) * t559 + (t331 * t62 - t332 * t61) * t560 + (t331 * t58 - t332 * t57) * t561 + (t331 * t54 - t332 * t53) * t562 + (t331 * t56 - t332 * t55) * t563 + (-t19 * t332 + t20 * t331) * t548 + (t331 * t37 - t332 * t36) * t550 + (t331 * t39 - t332 * t38) * t552 + (-t74 * (t211 * t284 - t449) - t75 * (-t211 * t283 - t450) - ((-t135 * t74 + t136 * t75) * t324 + (t74 * (t210 * t331 + t181) + t75 * (-t210 * t332 - t182) + t368) * t325) * qJD(4) + t41 * t481 + (t41 * t136 + t483 * t49 + t490 * t74) * t332 + (t41 * t135 + t483 * t50 + t490 * t75) * t331 - g(1) * (t309 + t477) - g(2) * (t308 + t478) - g(3) * (t211 + t570) - t569 * t324 * (-pkin(3) - t540) + (t100 * t332 - t181 * t283 + t182 * t284 + t99 * t331 - t482 + t484) * t63) * m(5) + (t258 * qJD(3) * t564 + (-t332 * t259 + t338) * t322) * t437 - (t259 * qJD(3) * t565 + (-t331 * t258 + t338) * t473) * t322 / 0.2e1; (t198 * t65 + t199 * t64 + t282 * t88 + t283 * t33 + t284 * t32 - t47 * t471) * t547 + (-t325 * t47 + (t32 * t331 + t33 * t332) * t324 + (t88 * t324 + t325 * t415) * qJD(3)) * t436 + t346 + t40 * t438 + (t324 * t415 - t325 * t88) * t554 + (t324 * t416 - t325 * t73) * t560 + (t324 * t417 - t325 * t72) * t559 + (-t325 * t46 + (t331 * t38 + t332 * t39) * t324 + (t325 * t416 + t531) * qJD(3)) * t552 + (-t325 * t45 + (t331 * t36 + t332 * t37) * t324 + (t325 * t417 + t532) * qJD(3)) * t550 + (t356 * t325 + (-t333 * t340 + t341 * t334) * t324) * t435 + (t278 * t340 + t279 * t341 - t332 * t348) * t553 + (t276 * t340 + t277 * t341 - t331 * t348) * t551 + t34 * t428 + t35 * t427 + t9 * t444 + t10 * t443 + (-g(1) * (t267 + t148) - g(2) * (t266 + t147) + t434 + t51 * t457 + t52 * t107 + t13 * t106 + t48 * t529 + (t15 * t126 + t51 * t102 - t16 * t491 - t52 * t530 + ((t48 * t126 + t487 * t52) * t332 + (-t48 * t491 - t51 * t578) * t331) * qJD(3)) * t325 - t51 * (t266 * t471 - t284 * t464 + t430) - t52 * (-t267 * t471 + t283 * t464 + t431) - t48 * (t266 * t283 - t267 * t284 + t432) + (-g(3) * (-pkin(4) * t333 + t422) + (-t492 * t51 + t533) * qJD(3) + (t48 * t102 + t13 * t126 + t16 * t487 + t495 * t52) * t332 + (t51 * t138 - t15 * t578 + t568) * t331) * t324) * m(6) + (-g(1) * t164 - g(2) * t163 - g(3) * t275 + (-t75 * t100 + t49 * t135 - t50 * t136 + t74 * t99 + (t368 + (t331 * t74 - t332 * t75) * t210) * qJD(3)) * t325 + (t74 * (-qJD(3) * t135 + t128 * t331) + t75 * (qJD(3) * t136 - t128 * t332) + t41 * t397 + t63 * (-t100 * t331 + t332 * t99) + t421 * t210) * t324 - t74 * (t163 * t471 + t275 * t284) - t75 * (-t164 * t471 - t275 * t283) - t63 * (t163 * t283 - t164 * t284)) * m(5); t346 + (t434 + t16 * (-t122 * t325 - t196 * t510) + t13 * (-t122 * t511 + t106) - g(1) * t148 - g(2) * t147 - g(3) * t247 + (-t113 * t510 + t107 + (-t196 * t473 - t87) * t325 - t431) * t52 + (-t121 * t475 - t430 + t457) * t51 + ((-t122 * t474 - t324 * t87) * t331 + t529 - t432) * t48) * m(6);];
tau = t1;
