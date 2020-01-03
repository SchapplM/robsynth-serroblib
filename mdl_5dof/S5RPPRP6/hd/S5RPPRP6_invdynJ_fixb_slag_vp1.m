% Calculate vector of inverse dynamics joint torques for
% S5RPPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:55:15
% DurationCPUTime: 14.39s
% Computational Cost: add. (6848->531), mult. (9604->636), div. (0->0), fcn. (7239->6), ass. (0->286)
t566 = Icges(6,4) + Icges(5,5);
t565 = Icges(5,6) - Icges(6,6);
t273 = sin(qJ(1));
t274 = cos(qJ(1));
t269 = pkin(7) + qJ(4);
t246 = cos(t269);
t245 = sin(t269);
t437 = Icges(5,4) * t245;
t318 = Icges(5,2) * t246 + t437;
t564 = t565 * t273 - t274 * t318;
t436 = Icges(5,4) * t246;
t320 = Icges(5,1) * t245 + t436;
t563 = -t566 * t273 + t274 * t320;
t418 = t246 * t274;
t420 = t245 * t274;
t559 = Icges(6,5) * t420 - Icges(6,3) * t418 + t564;
t213 = Icges(6,5) * t418;
t556 = -Icges(6,1) * t420 + t213 - t563;
t433 = Icges(6,5) * t245;
t178 = Icges(6,1) * t246 + t433;
t180 = Icges(5,1) * t246 - t437;
t552 = t178 + t180;
t316 = Icges(5,5) * t245 + Icges(5,6) * t246;
t109 = Icges(5,3) * t274 + t273 * t316;
t317 = Icges(6,4) * t245 - Icges(6,6) * t246;
t111 = Icges(6,2) * t274 + t273 * t317;
t562 = t109 + t111;
t530 = t556 * t245 + t559 * t246;
t561 = t530 * t274;
t421 = t245 * t273;
t212 = Icges(6,5) * t421;
t419 = t246 * t273;
t430 = Icges(6,6) * t274;
t107 = -Icges(6,3) * t419 + t212 + t430;
t113 = Icges(5,6) * t274 + t273 * t318;
t560 = t107 - t113;
t110 = -Icges(5,3) * t273 + t274 * t316;
t432 = Icges(6,2) * t273;
t112 = Icges(6,4) * t420 - Icges(6,6) * t418 - t432;
t558 = -t110 - t112;
t237 = Icges(6,5) * t246;
t319 = Icges(6,1) * t245 - t237;
t115 = Icges(6,4) * t274 + t273 * t319;
t214 = Icges(5,4) * t419;
t434 = Icges(5,5) * t274;
t117 = Icges(5,1) * t421 + t214 + t434;
t557 = t115 + t117;
t555 = t552 * t274;
t554 = -t565 * t245 + t566 * t246;
t176 = -Icges(5,2) * t245 + t436;
t492 = Icges(6,3) * t245 + t237;
t553 = -t176 + t492;
t314 = -Icges(6,3) * t246 + t433;
t551 = t314 - t318;
t550 = -t319 - t320;
t495 = -t107 * t418 - t273 * t562;
t304 = t246 * t176 + t245 * t180;
t541 = t245 * t178 - t246 * t492 + t304;
t310 = t113 * t246 + t117 * t245;
t491 = t115 * t245 + t310;
t439 = t274 * t111 + t115 * t421;
t37 = -t107 * t419 + t439;
t39 = t274 * t109 + t113 * t419 + t117 * t421;
t507 = t37 + t39;
t506 = t558 * t274 + t419 * t559 + t421 * t556;
t505 = -t115 * t420 - t310 * t274 - t495;
t504 = t273 * t558 - t561;
t137 = t176 * t274;
t375 = qJD(4) * t274;
t549 = qJD(4) * t137 - t492 * t375 + (t273 * t314 - t113 + t430) * qJD(1);
t130 = t492 * t273;
t376 = qJD(4) * t273;
t548 = qJD(4) * t130 - t176 * t376 + (t274 * t314 + t564) * qJD(1);
t547 = -t555 * qJD(4) + (t273 * t320 + t115 + t434) * qJD(1);
t140 = t180 * t273;
t546 = qJD(4) * t140 + t178 * t376 + (t274 * t319 + t563) * qJD(1);
t500 = -t245 * t559 + t246 * t556;
t501 = t245 * t560 + t246 * t557;
t521 = t554 * t273;
t545 = t551 * qJD(4);
t544 = t550 * qJD(4);
t543 = (Icges(6,1) * t419 + t140 + t212 + t560) * t274 + (-t555 - t559) * t273;
t542 = t245 * t553 + t246 * t552;
t425 = t107 * t246;
t540 = t425 - t491;
t539 = -t316 - t317;
t517 = t554 * t274;
t503 = t541 * t273 + t517;
t502 = t178 * t420 + t274 * t304 - t418 * t492 - t521;
t538 = (Icges(5,2) * t421 + t130 - t214 - t557) * t274 + (-Icges(6,3) * t420 + t137 - t213 - t556) * t273;
t537 = t562 * qJD(1);
t536 = t550 + t553;
t535 = t551 + t552;
t534 = t274 ^ 2;
t468 = rSges(6,1) + pkin(4);
t441 = rSges(6,3) + qJ(5);
t533 = t541 * qJD(1) + qJD(4) * t539;
t382 = qJD(1) * t110;
t532 = t382 + t521 * qJD(4) + (t274 * t317 - t432 - t540) * qJD(1);
t531 = qJD(1) * t530 - t517 * qJD(4) + t537;
t529 = t273 * t504 + t274 * t505;
t528 = t506 * t273 + t507 * t274;
t527 = -qJD(4) * t501 - t245 * t546 + t246 * t548 + t537;
t526 = -qJD(1) * t554 + t542 * qJD(4) + t544 * t245 + t545 * t246;
t525 = qJD(1) * t112 + qJD(4) * t500 + t245 * t547 - t246 * t549 + t382;
t524 = (t536 * t245 + t535 * t246) * qJD(1);
t523 = -t245 * t543 + t538 * t246;
t522 = t503 * qJD(1);
t184 = pkin(4) * t246 + qJ(5) * t245;
t185 = rSges(6,1) * t246 + rSges(6,3) * t245;
t394 = t184 + t185;
t236 = t246 * qJ(5);
t238 = t246 * rSges(6,3);
t330 = rSges(6,1) * t245 - t238;
t395 = -pkin(4) * t245 + t236 - t330;
t377 = qJD(1) * t274;
t520 = t245 * t377 + t246 * t376;
t519 = t539 * qJD(1);
t518 = t502 * qJD(1);
t515 = -qJ(3) * qJD(1) ^ 2 + qJDD(3);
t513 = t528 * qJD(4) + t522;
t512 = t529 * qJD(4) - t518;
t511 = t533 * t273 - t526 * t274;
t510 = t526 * t273 + t533 * t274;
t509 = qJD(4) * t540 + t245 * t548 + t246 * t546;
t508 = -t530 * qJD(4) + t245 * t549 + t246 * t547;
t499 = t245 * t468;
t270 = sin(pkin(7));
t417 = t270 * t273;
t271 = cos(pkin(7));
t451 = rSges(4,2) * t271;
t128 = rSges(4,1) * t417 + t274 * rSges(4,3) + t273 * t451;
t205 = t274 * pkin(1) + t273 * qJ(2);
t427 = qJ(3) * t274;
t345 = t205 + t427;
t494 = t128 + t345;
t391 = -t273 * rSges(6,2) - rSges(6,3) * t418;
t267 = t274 * rSges(6,2);
t493 = t468 * t421 + t267;
t206 = -rSges(3,2) * t274 + t273 * rSges(3,3);
t490 = t532 * t534 + (t525 * t273 + (-t527 + t531) * t274) * t273;
t489 = t527 * t534 + (t531 * t273 + (-t525 + t532) * t274) * t273;
t356 = t245 * t376;
t488 = -t441 * t356 - t468 * t520;
t487 = -qJ(5) * t418 + t391;
t233 = qJD(5) * t245;
t405 = qJD(4) * t395 + t233;
t480 = qJD(4) * (t233 + t405) - qJDD(5) * t246;
t477 = t245 * t441 + t246 * t468;
t475 = t273 ^ 2;
t369 = qJD(1) * qJD(4);
t189 = qJDD(4) * t273 + t274 * t369;
t473 = t189 / 0.2e1;
t190 = qJDD(4) * t274 - t273 * t369;
t472 = t190 / 0.2e1;
t471 = t273 / 0.2e1;
t470 = -t274 / 0.2e1;
t467 = rSges(3,2) - pkin(1);
t466 = -rSges(6,2) - pkin(1);
t465 = -rSges(5,3) - pkin(1);
t464 = pkin(3) * t270;
t462 = g(2) * t274;
t232 = pkin(3) * t417;
t272 = -pkin(6) - qJ(3);
t414 = qJ(3) + t272;
t155 = -t274 * t414 + t232;
t371 = qJD(1) * qJD(2);
t378 = qJD(1) * t273;
t242 = qJ(2) * t377;
t251 = qJD(2) * t273;
t387 = t242 + t251;
t364 = qJD(1) * (-pkin(1) * t378 + t387) + qJDD(1) * t205 + t273 * t371;
t370 = qJD(1) * qJD(3);
t277 = qJDD(1) * t427 + t273 * t515 + 0.2e1 * t274 * t370 + t364;
t357 = t270 * t377;
t389 = pkin(3) * t357 + t272 * t378;
t276 = qJD(1) * (qJ(3) * t378 + t389) + qJDD(1) * t155 + t277;
t373 = qJD(5) * t273;
t354 = t246 * t373;
t404 = -t419 * t441 + t493;
t455 = -(-qJ(5) * t377 - t373) * t246 - t391 * qJD(1) + t488;
t1 = -t394 * t190 + t404 * qJDD(1) + (-t354 - t455) * qJD(1) + (-qJDD(2) - t480) * t274 + t276;
t461 = t1 * t274;
t388 = qJDD(2) * t273 + t274 * t371;
t291 = -0.2e1 * t273 * t370 + t274 * t515 + t388;
t154 = t273 * t414 + t274 * t464;
t254 = t274 * qJ(2);
t202 = pkin(1) * t273 - t254;
t428 = qJ(3) * t273;
t346 = -t202 - t428;
t337 = t154 + t346;
t403 = t420 * t468 + t487;
t297 = t337 + t403;
t374 = qJD(5) * t246;
t252 = qJD(2) * t274;
t166 = qJD(1) * t205 - t252;
t229 = t272 * t377;
t402 = t229 - (t232 - t427) * qJD(1) - t166;
t148 = t185 * t274;
t456 = (pkin(4) * t378 - qJ(5) * t375) * t245 + (-qJ(5) * t378 + (-pkin(4) * qJD(4) + qJD(5)) * t274) * t246 - qJD(4) * t148 + (t273 * t330 + t267) * qJD(1);
t2 = t394 * t189 + t480 * t273 + t297 * qJDD(1) + (-t274 * t374 + t402 - t456) * qJD(1) + t291;
t460 = t2 * t273;
t457 = -pkin(1) - qJ(3);
t454 = rSges(5,1) * t246;
t450 = rSges(5,2) * t245;
t449 = rSges(5,2) * t246;
t448 = rSges(3,3) * t274;
t331 = rSges(5,1) * t245 + t449;
t168 = t331 * qJD(4);
t186 = -t450 + t454;
t261 = t273 * rSges(5,3);
t123 = t274 * t331 - t261;
t302 = t123 + t337;
t149 = t186 * t274;
t265 = t274 * rSges(5,3);
t82 = -qJD(4) * t149 + (t273 * t331 + t265) * qJD(1);
t16 = -t168 * t376 + t186 * t189 + (-t82 + t402) * qJD(1) + t302 * qJDD(1) + t291;
t447 = t16 * t273;
t121 = rSges(5,1) * t421 + rSges(5,2) * t419 + t265;
t361 = rSges(5,1) * t520 + t377 * t449;
t365 = qJD(4) * t450;
t84 = (-rSges(5,3) * qJD(1) - t365) * t273 + t361;
t17 = qJD(1) * t84 + qJDD(1) * t121 - t186 * t190 + (qJD(4) * t168 - qJDD(2)) * t274 + t276;
t446 = t17 * t274;
t153 = t186 * t376;
t250 = qJD(3) * t274;
t385 = t250 + t251;
t45 = qJD(1) * t302 + t153 + t385;
t444 = t274 * t45;
t33 = t233 + (-t273 * t404 - t274 * t403) * qJD(4);
t443 = t33 * t246;
t401 = t419 * t468 + t441 * t421;
t400 = -t184 * t274 - t148;
t203 = rSges(3,2) * t273 + t448;
t393 = -t202 + t203;
t157 = t205 + t206;
t390 = rSges(4,1) * t357 + t377 * t451;
t386 = rSges(3,2) * t378 + rSges(3,3) * t377;
t192 = qJD(1) * t202;
t384 = t251 - t192;
t372 = -m(4) - m(5) - m(6);
t368 = qJDD(2) * t274;
t366 = -rSges(4,3) + t457;
t360 = t242 + t385;
t359 = t250 + t384;
t351 = -t376 / 0.2e1;
t350 = t376 / 0.2e1;
t349 = -t375 / 0.2e1;
t348 = t375 / 0.2e1;
t347 = t246 * t441;
t344 = t112 - t425;
t343 = qJD(3) * t273 - t252;
t342 = g(1) * t273 - t462;
t332 = rSges(4,1) * t270 + t451;
t129 = -t273 * rSges(4,3) + t274 * t332;
t338 = t129 + t346;
t336 = t155 + t345;
t335 = t229 - t343;
t333 = t254 + (-pkin(1) + t272) * t273;
t207 = rSges(2,1) * t274 - rSges(2,2) * t273;
t204 = rSges(2,1) * t273 + rSges(2,2) * t274;
t300 = -t354 + t385;
t290 = t376 * t394 + t300;
t29 = qJD(1) * t297 + t290;
t30 = (-qJD(4) * t394 + t374) * t274 + (t336 + t404) * qJD(1) + t343;
t327 = t273 * t29 - t274 * t30;
t46 = -t186 * t375 + (t121 + t336) * qJD(1) + t343;
t322 = t273 * t45 - t274 * t46;
t321 = -t273 * t84 + t274 * t82;
t307 = -t121 * t273 - t123 * t274;
t301 = -t272 * t274 + t205 + t232;
t299 = t327 * t245;
t296 = t331 + t464;
t145 = t186 * t273;
t127 = qJD(1) * t154;
t99 = qJD(1) * t157 - t252;
t98 = qJD(1) * t393 + t251;
t80 = qJD(1) * t494 + t343;
t79 = qJD(1) * t338 + t385;
t61 = t307 * qJD(4);
t48 = qJD(1) * t386 + qJDD(1) * t206 + t364 - t368;
t47 = t393 * qJDD(1) + (-qJD(1) * t206 - t166) * qJD(1) + t388;
t32 = -t368 + qJDD(1) * t128 + qJD(1) * (-rSges(4,3) * t378 + t390) + t277;
t31 = t338 * qJDD(1) + (-qJD(1) * t128 - t166) * qJD(1) + t291;
t3 = qJDD(5) * t245 - t403 * t190 - t404 * t189 + (t273 * t455 + t274 * t456 + t374) * qJD(4);
t4 = [-m(2) * (-g(1) * t204 + g(2) * t207) - t502 * t189 / 0.2e1 + t500 * t473 + (((-t495 - t505 + t506) * t273 + (t439 + t39 + t561 + (t344 + t110 - t491) * t273 + t504) * t274) * qJD(4) + t522) * t351 + (-t541 * qJD(4) - t545 * t245 + t544 * t246) * qJD(1) + (t29 * t335 + t30 * (t242 + t300 + t389 - t488) + (qJD(4) * t477 - t374) * t274 * t29 + ((t29 * t466 - t30 * t347) * t274 + (t29 * (-qJ(2) - t464 + t395) + t30 * t466) * t273) * qJD(1) - (t127 - t192 - t29 + (t403 - t428) * qJD(1) + t290) * t30 + (t1 - g(2)) * (-t273 * t347 + t301 + t493) + (t2 - g(1)) * ((t464 + t499) * t274 + t333 + t487)) * m(6) + (t45 * (-t274 * t365 + t375 * t454 + t335) + t46 * (-rSges(5,2) * t356 + t360 + t361 + t389) + (t465 * t444 + (t45 * (-qJ(2) - t296) + t46 * t465) * t273) * qJD(1) - (t127 + t153 - t45 + (t123 - t428) * qJD(1) + t359) * t46 + (t17 - g(2)) * (t301 + t121) + (t16 - g(1)) * (t274 * t296 - t261 + t333)) * m(5) + (-t79 * t343 + t80 * (t360 + t390) + (t79 * t366 * t274 + (t79 * (-qJ(2) - t332) + t80 * t366) * t273) * qJD(1) - (-t79 + (t129 - t428) * qJD(1) + t359) * t80 + (t32 - g(2)) * t494 + (t31 - g(1)) * (t273 * t457 + t129 + t254)) * m(4) + (t98 * t252 + t99 * (t386 + t387) + (t98 * t467 * t274 + (t98 * (-rSges(3,3) - qJ(2)) - t99 * pkin(1)) * t273) * qJD(1) - (qJD(1) * t203 + t384 - t98) * t99 + (t48 - g(2)) * t157 + (t47 - g(1)) * (t273 * t467 + t254 + t448)) * m(3) + (t501 + t503) * t472 + (t509 + t510) * t348 + ((t475 * t110 + (t273 * t344 - t37 + t439) * t273 + ((t491 - t558) * t274 + t495 + t506) * t274) * qJD(4) + t512 + t518) * t349 + (t508 + t511 + t513) * t350 + (m(2) * (t204 ^ 2 + t207 ^ 2) + Icges(4,1) * t271 ^ 2 + (-0.2e1 * Icges(4,4) * t271 + Icges(4,2) * t270) * t270 + Icges(2,3) + Icges(3,1) + t542) * qJDD(1); (-m(3) + t372) * t342 + 0.2e1 * (-t461 / 0.2e1 + t460 / 0.2e1) * m(6) + 0.2e1 * (t447 / 0.2e1 - t446 / 0.2e1) * m(5) + 0.2e1 * (t31 * t471 + t32 * t470) * m(4) + 0.2e1 * (t47 * t471 + t470 * t48) * m(3); t372 * (g(1) * t274 + g(2) * t273) + m(4) * (t273 * t32 + t274 * t31) + m(5) * (t16 * t274 + t17 * t273) + m(6) * (t1 * t273 + t2 * t274); t529 * t473 + t528 * t472 + (qJD(1) * t511 + qJD(4) * t489 - t502 * qJDD(1) + t189 * t504 + t190 * t505) * t471 + (t510 * qJD(1) + t490 * qJD(4) + t503 * qJDD(1) + t506 * t189 + t507 * t190) * t274 / 0.2e1 - ((t538 * t245 + t246 * t543) * qJD(4) + (-t535 * t245 + t536 * t246) * qJD(1)) * qJD(1) / 0.2e1 + (t509 * t274 + t508 * t273 + (-t501 * t273 + t500 * t274) * qJD(1)) * qJD(1) / 0.2e1 + (t500 * t273 + t501 * t274) * qJDD(1) / 0.2e1 - t513 * t378 / 0.2e1 + t512 * t377 / 0.2e1 + ((-t376 * t517 + t519) * t273 + ((t273 * t521 + t523) * qJD(4) - t524) * t274) * t351 + ((-t505 * t273 + t504 * t274) * qJD(1) + t489) * t350 + ((t375 * t521 + t519) * t274 + ((-t274 * t517 - t523) * qJD(4) + t524) * t273) * t349 + ((-t273 * t507 + t506 * t274) * qJD(1) + t490) * t348 + (-g(1) * t401 - g(3) * (t236 + t238 - t499) + t477 * t462 + (-t1 * t394 - t30 * t405 - t3 * t403 + t33 * t456 + (t29 * t394 - t33 * t404) * qJD(1)) * t274 + (t2 * t394 + t29 * t405 - t3 * t404 + t33 * t455 + (t30 * t394 + t33 * t403) * qJD(1)) * t273 - (t299 + t443) * qJD(5) - (-t29 * t400 + t30 * t401) * qJD(1) - ((-t30 * t395 + t33 * t400) * t274 + (t29 * t395 - t33 * t401) * t273) * qJD(4)) * m(6) + ((qJD(4) * t321 - t121 * t189 - t123 * t190) * t307 + t61 * ((-t121 * t274 + t123 * t273) * qJD(1) + t321) - t322 * t168 + (t447 - t446 + (t273 * t46 + t444) * qJD(1)) * t186 - (t145 * t46 + t149 * t45) * qJD(1) - (t61 * (-t145 * t273 - t149 * t274) - t322 * t331) * qJD(4) - g(1) * t145 + g(2) * t149 + g(3) * t331) * m(5); (-((t475 + t534) * t443 + t299) * qJD(4) + (qJD(4) * t327 - g(3) + t3) * t245 + (qJD(4) * t33 + t342 - t460 + t461) * t246) * m(6);];
tau = t4;
