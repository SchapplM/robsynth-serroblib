% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:13
% EndTime: 2022-01-23 09:12:40
% DurationCPUTime: 24.04s
% Computational Cost: add. (12208->554), mult. (16242->722), div. (0->0), fcn. (15498->8), ass. (0->267)
t551 = Icges(5,4) + Icges(6,4);
t276 = qJ(1) + pkin(7);
t272 = sin(t276);
t278 = cos(pkin(8));
t282 = cos(qJ(4));
t402 = t278 * t282;
t273 = cos(t276);
t280 = sin(qJ(4));
t405 = t273 * t280;
t204 = t272 * t402 - t405;
t174 = Icges(6,4) * t204;
t403 = t278 * t280;
t203 = t272 * t403 + t273 * t282;
t277 = sin(pkin(8));
t411 = t272 * t277;
t113 = -Icges(6,2) * t203 + Icges(6,6) * t411 + t174;
t177 = Icges(5,4) * t204;
t116 = -Icges(5,2) * t203 + Icges(5,6) * t411 + t177;
t507 = t116 + t113;
t173 = Icges(6,4) * t203;
t120 = -Icges(6,1) * t204 - Icges(6,5) * t411 + t173;
t176 = Icges(5,4) * t203;
t123 = -Icges(5,1) * t204 - Icges(5,5) * t411 + t176;
t542 = t123 + t120;
t107 = Icges(6,5) * t204 - Icges(6,6) * t203 + Icges(6,3) * t411;
t110 = Icges(5,5) * t204 - Icges(5,6) * t203 + Icges(5,3) * t411;
t510 = t107 + t110;
t511 = Icges(5,1) + Icges(6,1);
t544 = Icges(5,5) + Icges(6,5);
t529 = Icges(5,2) + Icges(6,2);
t550 = Icges(5,6) + Icges(6,6);
t549 = Icges(5,3) + Icges(6,3);
t556 = t551 * t282;
t555 = t551 * t280;
t554 = t507 * t203 + t542 * t204;
t490 = t510 * t411 - t554;
t205 = t272 * t282 - t273 * t403;
t409 = t272 * t280;
t206 = t273 * t402 + t409;
t407 = t273 * t277;
t109 = Icges(6,5) * t206 + Icges(6,6) * t205 + Icges(6,3) * t407;
t112 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t407;
t509 = t109 + t112;
t416 = Icges(6,4) * t206;
t115 = Icges(6,2) * t205 + Icges(6,6) * t407 + t416;
t419 = Icges(5,4) * t206;
t118 = Icges(5,2) * t205 + Icges(5,6) * t407 + t419;
t508 = t115 + t118;
t175 = Icges(6,4) * t205;
t121 = Icges(6,1) * t206 + Icges(6,5) * t407 + t175;
t178 = Icges(5,4) * t205;
t124 = Icges(5,1) * t206 + Icges(5,5) * t407 + t178;
t506 = t121 + t124;
t541 = -t549 * t278 + (-t550 * t280 + t544 * t282) * t277;
t540 = t550 * t278 + (t529 * t280 - t556) * t277;
t539 = -t544 * t278 + (t511 * t282 - t555) * t277;
t553 = t510 * t407;
t552 = t507 * t205 - t542 * t206;
t489 = -t508 * t203 + t506 * t204 + t509 * t411;
t488 = t552 + t553;
t487 = t508 * t205 + t506 * t206 + t509 * t407;
t548 = (-t544 * t280 - t550 * t282) * t277;
t547 = (t529 * t282 + t555) * t277;
t546 = (t511 * t280 + t556) * t277;
t545 = t490 * t272;
t484 = t540 * t203 + t539 * t204 + t541 * t411;
t483 = -t540 * t205 + t539 * t206 + t541 * t407;
t137 = qJD(1) * t203 - qJD(4) * t206;
t138 = -qJD(1) * t204 + qJD(4) * t205;
t376 = qJD(1) * t277;
t356 = t272 * t376;
t518 = t550 * t137 + t544 * t138 - t549 * t356;
t139 = qJD(1) * t205 - qJD(4) * t204;
t293 = t203 * qJD(4);
t140 = qJD(1) * t206 - t293;
t355 = t273 * t376;
t517 = t550 * t139 + t544 * t140 + t549 * t355;
t516 = -t529 * t137 - t551 * t138 + t550 * t356;
t515 = t529 * t139 + t551 * t140 + t550 * t355;
t514 = t551 * t137 + t511 * t138 - t544 * t356;
t513 = t551 * t139 + t511 * t140 + t544 * t355;
t538 = t548 * qJD(4);
t537 = t547 * qJD(4);
t536 = t546 * qJD(4);
t535 = (t488 * t272 + t487 * t273) * t277;
t534 = (t489 * t273 + t545) * t277;
t270 = pkin(4) * t282 + pkin(3);
t367 = pkin(4) * t405;
t410 = t272 * t278;
t533 = -t204 * rSges(6,1) + t203 * rSges(6,2) - t270 * t410 + t367;
t373 = qJD(5) * t277;
t243 = t273 * t373;
t260 = -qJD(4) * t278 + qJD(1);
t375 = qJD(4) * t277;
t354 = t272 * t375;
t279 = -qJ(5) - pkin(6);
t440 = pkin(6) + t279;
t441 = pkin(3) - t270;
t390 = (t440 - rSges(6,3)) * t278 + (rSges(6,1) * t282 - rSges(6,2) * t280 - t441) * t277;
t445 = pkin(3) * t278;
t465 = t277 * t440;
t396 = (t445 + t465) * t272 - rSges(6,3) * t411 + t533;
t532 = t260 * t396 + t354 * t390 + t243;
t283 = cos(qJ(1));
t274 = t283 * pkin(1);
t531 = t483 * t260;
t233 = rSges(3,1) * t272 + rSges(3,2) * t273;
t281 = sin(qJ(1));
t446 = pkin(1) * t281;
t215 = -t233 - t446;
t528 = t484 * t260;
t388 = -t204 * rSges(5,1) + t203 * rSges(5,2);
t128 = rSges(5,3) * t411 - t388;
t202 = -rSges(5,3) * t278 + (rSges(5,1) * t282 - rSges(5,2) * t280) * t277;
t527 = -t128 * t260 + t202 * t354;
t369 = qJD(1) * qJD(4);
t171 = (qJDD(4) * t272 + t273 * t369) * t277;
t259 = -qJDD(4) * t278 + qJDD(1);
t234 = t273 * pkin(2) + t272 * qJ(3);
t263 = qJD(3) * t273;
t167 = qJD(1) * t234 - t263;
t370 = qJD(1) * qJD(3);
t284 = qJD(1) ^ 2;
t524 = t284 * t274;
t312 = qJDD(3) * t272 + t273 * t370 - t524;
t443 = pkin(6) * t277;
t332 = t443 + t445;
t207 = t332 * t272;
t265 = t273 * qJ(3);
t232 = pkin(2) * t272 - t265;
t347 = -t232 - t446;
t333 = -t207 + t347;
t377 = qJD(1) * t273;
t285 = (-t332 * t377 - t167) * qJD(1) + t333 * qJDD(1) + t312;
t252 = pkin(4) * t409;
t323 = rSges(6,1) * t140 + rSges(6,2) * t139;
t464 = t278 * t441;
t436 = rSges(6,3) * t355 + t323 + t272 * t373 - pkin(4) * t293 + (t252 + (-t464 - t465) * t273) * qJD(1);
t434 = rSges(6,2) * t282;
t306 = (-rSges(6,1) * t280 - t434) * t277;
t433 = pkin(4) * qJD(4);
t362 = t280 * t433;
t372 = qJD(5) * t278;
t383 = qJD(4) * t306 - t277 * t362 - t372;
t457 = -qJD(1) * qJD(5) + qJD(4) * t383;
t10 = -t436 * t260 + t396 * t259 + t390 * t171 + (qJDD(5) * t273 + t272 * t457) * t277 + t285;
t526 = -g(1) + t10;
t172 = (qJDD(4) * t273 - t272 * t369) * t277;
t406 = t273 * t278;
t208 = pkin(3) * t406 + pkin(6) * t407;
t336 = qJDD(1) * t274 - t284 * t446;
t378 = qJD(1) * t272;
t262 = qJD(3) * t272;
t380 = qJ(3) * t377 + t262;
t297 = qJD(1) * (-pkin(2) * t378 + t380) + qJDD(1) * t234 + t272 * t370 + t336;
t287 = qJDD(1) * t208 - t207 * t284 + t297;
t404 = t277 * t279;
t463 = t206 * rSges(6,1) + t205 * rSges(6,2) + rSges(6,3) * t407 + t270 * t406 - t273 * t404 + t252;
t395 = -t208 + t463;
t343 = t278 * t362;
t361 = t282 * t433;
t461 = t138 * rSges(6,1) + t137 * rSges(6,2) + qJD(1) * t367 + t272 * t361 - t273 * t343 + t279 * t356 + t243;
t437 = -rSges(6,3) * t356 + (t443 + t464) * t378 + t461;
t11 = qJDD(5) * t411 + t437 * t260 + t395 * t259 - t390 * t172 + (-t277 * t457 - qJDD(3)) * t273 + t287;
t525 = -g(2) + t11;
t523 = qJD(4) * t534 + t528;
t522 = qJD(4) * t535 + t531;
t521 = -t517 * t278 + (t513 * t282 - t515 * t280 + (t542 * t280 - t507 * t282) * qJD(4)) * t277;
t520 = -t518 * t278 + (t514 * t282 + t516 * t280 + (-t506 * t280 - t508 * t282) * qJD(4)) * t277;
t492 = (-t538 * t273 + t541 * t378) * t277 + t536 * t206 + t537 * t205 - t539 * t138 + t540 * t137;
t491 = (-t538 * t272 - t541 * t377) * t277 + t536 * t204 - t537 * t203 - t539 * t140 + t540 * t139;
t482 = t107 * t278;
t497 = t113 * t280 + t120 * t282;
t40 = -t277 * t497 - t482;
t480 = t110 * t278;
t494 = t116 * t280 + t123 * t282;
t42 = -t277 * t494 - t480;
t486 = t40 + t42;
t519 = t538 * t278 + (t536 * t282 - t537 * t280 + (t539 * t280 - t540 * t282) * qJD(4)) * t277;
t512 = t541 * t278 + (-t540 * t280 - t539 * t282) * t277;
t504 = (t544 * t205 - t206 * t550) * t273 + (-t544 * t203 - t204 * t550) * t272;
t503 = t488 - t553;
t502 = rSges(6,1) + pkin(4);
t476 = (t206 * t529 - t175 - t178 - t506) * t273 + (t204 * t529 + t173 + t176 + t542) * t272;
t41 = -t109 * t278 + (-t115 * t280 + t121 * t282) * t277;
t43 = -t112 * t278 + (-t118 * t280 + t124 * t282) * t277;
t485 = t41 + t43;
t163 = rSges(4,1) * t406 - rSges(4,2) * t407 + t272 * rSges(4,3);
t462 = t274 + t234;
t158 = t163 + t462;
t478 = t208 + t462;
t475 = (t511 * t205 - t416 - t419 - t508) * t273 + (-t511 * t203 - t174 - t177 - t507) * t272;
t474 = t504 * t277;
t473 = t521 * t272 + t520 * t273;
t472 = (t513 * t204 - t515 * t203 - t542 * t140 + t507 * t139 + (t517 * t272 + t510 * t377) * t277) * t272 + ((t518 * t272 + t509 * t377) * t277 + t514 * t204 + t516 * t203 + t506 * t140 + t508 * t139) * t273;
t471 = (t514 * t206 - t516 * t205 + t506 * t138 + t508 * t137 + (t518 * t273 - t509 * t378) * t277) * t273 + ((t517 * t273 - t510 * t378) * t277 + t513 * t206 + t515 * t205 - t542 * t138 + t507 * t137) * t272;
t470 = t539 - t547;
t469 = t540 - t546;
t468 = -t512 * t259 - t519 * t260;
t235 = t273 * rSges(3,1) - rSges(3,2) * t272;
t216 = t235 + t274;
t454 = t171 / 0.2e1;
t453 = t172 / 0.2e1;
t449 = t272 / 0.2e1;
t448 = -t273 / 0.2e1;
t444 = pkin(4) * t280;
t442 = g(1) * t272;
t228 = (-rSges(5,1) * t280 - rSges(5,2) * t282) * t277;
t218 = qJD(4) * t228;
t325 = -rSges(5,1) * t140 - rSges(5,2) * t139;
t77 = rSges(5,3) * t355 - t325;
t21 = -t128 * t259 + t171 * t202 + t218 * t354 - t260 * t77 + t285;
t432 = t21 * t272;
t132 = t206 * rSges(5,1) + t205 * rSges(5,2) + rSges(5,3) * t407;
t393 = t138 * rSges(5,1) + t137 * rSges(5,2);
t75 = -rSges(5,3) * t356 + t393;
t22 = t132 * t259 - t172 * t202 + t260 * t75 + (-t218 * t375 - qJDD(3)) * t273 + t287;
t431 = t22 * t273;
t296 = qJD(1) * t333 + t262;
t46 = t296 + t527;
t424 = t273 * t46;
t423 = t40 * t171;
t422 = t41 * t172;
t421 = t42 * t171;
t420 = t43 * t172;
t412 = t270 * t278;
t392 = -t204 * rSges(6,2) - t203 * t502;
t391 = -t206 * rSges(6,2) + t205 * t502;
t382 = rSges(4,2) * t356 + rSges(4,3) * t377;
t381 = rSges(4,2) * t411 + t273 * rSges(4,3);
t379 = -qJD(1) * t232 + t262;
t371 = -m(4) - m(5) - m(6);
t366 = rSges(4,1) * t410;
t351 = -rSges(4,1) * t278 - pkin(2);
t350 = -rSges(6,3) * t277 - pkin(2);
t349 = -t375 / 0.2e1;
t348 = t375 / 0.2e1;
t345 = t265 - t446;
t340 = t272 * t349;
t339 = t272 * t348;
t338 = t273 * t349;
t337 = t273 * t348;
t162 = t366 - t381;
t334 = -t162 + t347;
t327 = -pkin(2) + (-rSges(6,3) + t279) * t277;
t249 = rSges(2,1) * t283 - rSges(2,2) * t281;
t248 = rSges(2,1) * t281 + rSges(2,2) * t283;
t295 = qJD(1) * t478 - t263;
t47 = -t202 * t273 * t375 + t132 * t260 + t295;
t318 = t272 * t46 - t273 * t47;
t317 = -t272 * t75 + t273 * t77;
t316 = t128 * t273 - t132 * t272;
t301 = t379 + (-t207 - t446) * qJD(1);
t298 = -t445 - pkin(2) + (-rSges(5,3) - pkin(6)) * t277;
t156 = rSges(5,1) * t205 - rSges(5,2) * t206;
t154 = -rSges(5,1) * t203 - rSges(5,2) * t204;
t106 = qJD(1) * t158 - t263;
t105 = qJD(1) * t334 + t262;
t52 = t316 * t375 + qJD(2);
t51 = -qJDD(3) * t273 + qJDD(1) * t163 + qJD(1) * (-qJD(1) * t366 + t382) + t297;
t50 = t334 * qJDD(1) + (-qJD(1) * t163 - t167) * qJD(1) + t312;
t31 = t395 * t260 + (-qJD(4) * t273 * t390 + qJD(5) * t272) * t277 + t295;
t30 = t296 + t532;
t27 = -t372 + qJD(2) + (-t272 * t395 - t273 * t396) * t375;
t20 = t128 * t172 - t132 * t171 + t317 * t375 + qJDD(2);
t1 = -qJDD(5) * t278 + qJDD(2) - t396 * t172 - t395 * t171 + (-t272 * t437 + t273 * t436) * t375;
t2 = [t421 / 0.2e1 + t422 / 0.2e1 + t468 + t420 / 0.2e1 + t423 / 0.2e1 - m(2) * (-g(1) * t248 + g(2) * t249) + t484 * t454 + t483 * t453 + (((t487 + t490 + t554) * t273 + t503 * t272) * t375 + t531) * t340 + ((-t233 * t284 - g(2) + t336) * t216 + (-t524 + (-0.2e1 * t235 - t274 + t216) * t284 - g(1)) * t215) * m(3) + (m(2) * (t248 ^ 2 + t249 ^ 2) + Icges(4,2) * t278 ^ 2 + (Icges(4,1) * t277 + 0.2e1 * Icges(4,4) * t278) * t277 + m(3) * (t215 ^ 2 + t235 * t216) + Icges(2,3) + Icges(3,3)) * qJDD(1) + (t30 * (t273 * t361 + t263 - t323) + t31 * (t380 + t461) + (t10 * (t350 + t404) + t30 * (t343 - t373)) * t272 + ((-t281 * t31 - t283 * t30) * pkin(1) + t30 * (t327 - t412) * t273 + (t30 * (-qJ(3) - t444) + t31 * (t350 - t412)) * t272) * qJD(1) - t327 * t442 - (-t30 + t301 + t532) * t31 + t525 * (t462 + t463) + t526 * (t345 + t533)) * m(6) + (-(t301 - t46 + t527) * t47 + t46 * (t263 + t325) + t47 * (t380 + t393) + ((-t281 * t47 - t283 * t46) * pkin(1) + t298 * t424 + (-t46 * qJ(3) + t47 * (-rSges(5,3) * t277 - pkin(2) - t332)) * t272) * qJD(1) + (-g(2) + t22) * (t132 + t478) + (t21 - g(1)) * (t272 * t298 + t345 + t388)) * m(5) + (-(-t105 + (-t162 - t446) * qJD(1) + t379) * t106 + t105 * t263 + t106 * (t380 + t382) + ((-t105 * t283 - t106 * t281) * pkin(1) + t105 * (rSges(4,2) * t277 + t351) * t273 + (t105 * (-rSges(4,3) - qJ(3)) + t106 * t351) * t272) * qJD(1) + (t51 - g(2)) * t158 + (t50 - g(1)) * (t351 * t272 + t345 + t381)) * m(4) + (-t492 + t520) * t337 + (-t491 + t521 + t522) * t339 + ((t480 + t482 + (t494 + t497) * t277 + t486) * t260 + ((t503 - t489 - t552) * t273 - t545) * t375 + t523 - t528) * t338; (m(3) + m(4)) * qJDD(2) + m(5) * t20 + m(6) * t1 + (-m(3) + t371) * g(3); t371 * (-g(2) * t273 + t442) + 0.2e1 * (t10 * t449 + t11 * t448) * m(6) + 0.2e1 * (t432 / 0.2e1 - t431 / 0.2e1) * m(5) + 0.2e1 * (t448 * t51 + t449 * t50) * m(4); (-t278 * t484 + t534) * t454 + (-t278 * t483 + t535) * t453 + (t512 * t278 + (t272 * t486 + t273 * t485) * t277) * t259 / 0.2e1 - (((-t280 * t470 + t282 * t469) * t260 + ((t280 * t476 + t282 * t475) * t277 - t504 * t278) * qJD(4)) * t277 - t548 * t260 * t278) * t260 / 0.2e1 + (t519 * t278 + ((-t272 * t485 + t273 * t486) * qJD(1) + t473) * t277) * t260 / 0.2e1 - (t473 * t375 + t420 + t421 + t422 + t423 + t468) * t278 / 0.2e1 + (t490 * t171 + t489 * t172 + t484 * t259 - t491 * t260 + t472 * t375) * t411 / 0.2e1 - t522 * t356 / 0.2e1 + ((t203 * t476 + t204 * t475 + t272 * t474) * t375 + (-t203 * t470 + t204 * t469 + t411 * t548) * t260) * t340 + (t491 * t278 + ((-t272 * t489 + t273 * t490) * qJD(1) + t472) * t277) * t339 + ((-t205 * t476 + t475 * t206 + t474 * t273) * t375 + (t205 * t470 + t206 * t469 + t407 * t548) * t260) * t338 + (t492 * t278 + ((-t272 * t487 + t273 * t488) * qJD(1) + t471) * t277) * t337 + (-(-t30 * t392 + t31 * t391) * t260 - (t27 * (-t272 * t391 + t273 * t392) + (-t277 * t444 + t306) * (t272 * t30 - t273 * t31)) * t375 + (-t10 * t396 - t11 * t395 + t30 * t436 - t31 * t437) * t278 + ((-t11 * t390 - t31 * t383 - t1 * t396 + t27 * t436 + (-t27 * t395 + t30 * t390) * qJD(1)) * t273 + (t10 * t390 + t30 * t383 - t1 * t395 - t27 * t437 + (t27 * t396 + t31 * t390) * qJD(1)) * t272) * t277 - g(1) * t391 - g(2) * t392 - g(3) * (-t280 * t502 - t434) * t277) * m(6) + (-(-t154 * t46 + t156 * t47) * t260 - (t52 * (t154 * t273 - t156 * t272) + t318 * t228) * t375 + (t128 * t21 - t132 * t22 + t46 * t77 - t47 * t75) * t278 + (t20 * t316 + t52 * (-t128 * t378 - t132 * t377 + t317) + t318 * t218 + (t432 - t431 + (t272 * t47 + t424) * qJD(1)) * t202) * t277 - g(1) * t156 - g(2) * t154 - g(3) * t228) * m(5) + (t523 * qJD(1) + t488 * t171 + t487 * t172 + t483 * t259 - t492 * t260 + t471 * t375) * t407 / 0.2e1; ((-t1 + g(3)) * t278 + (t525 * t272 + t273 * t526) * t277) * m(6);];
tau = t2;
