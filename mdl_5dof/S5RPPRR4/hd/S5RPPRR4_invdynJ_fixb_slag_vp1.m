% Calculate vector of inverse dynamics joint torques for
% S5RPPRR4
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:05
% EndTime: 2019-12-05 17:44:51
% DurationCPUTime: 26.78s
% Computational Cost: add. (20580->833), mult. (23801->1105), div. (0->0), fcn. (22772->10), ass. (0->384)
t343 = cos(qJ(1));
t341 = -pkin(6) - qJ(3);
t338 = sin(pkin(8));
t513 = t338 * t343;
t303 = t341 * t513;
t337 = sin(pkin(9));
t342 = sin(qJ(1));
t506 = t342 * t337;
t447 = pkin(3) * t506;
t405 = -t303 + t447;
t339 = cos(pkin(9));
t320 = t339 * pkin(3) + pkin(2);
t335 = pkin(9) + qJ(4);
t324 = cos(t335);
t283 = pkin(4) * t324 + t320;
t323 = sin(t335);
t548 = pkin(4) * t323;
t549 = pkin(3) * t337;
t404 = t548 + t549;
t340 = cos(pkin(8));
t511 = t340 * t343;
t476 = -t283 * t511 - t342 * t404;
t333 = -pkin(7) + t341;
t515 = t333 * t338;
t517 = t320 * t340;
t110 = (t515 + t517) * t343 + t405 + t476;
t327 = qJ(5) + t335;
t318 = sin(t327);
t504 = t343 * t318;
t319 = cos(t327);
t509 = t342 * t319;
t234 = t340 * t504 - t509;
t235 = t318 * t342 + t319 * t511;
t392 = t235 * rSges(6,1) - t234 * rSges(6,2);
t134 = rSges(6,3) * t513 + t392;
t464 = t333 - t341;
t470 = t283 - t320;
t186 = t338 * t470 + t340 * t464;
t310 = qJDD(4) * t513;
t336 = qJD(4) + qJD(5);
t413 = qJD(1) * t336;
t194 = t310 + (qJDD(5) * t343 - t342 * t413) * t338;
t391 = -rSges(6,1) * t318 - rSges(6,2) * t319;
t236 = t391 * t338;
t201 = t336 * t236;
t208 = -rSges(6,3) * t340 + (rSges(6,1) * t319 - rSges(6,2) * t318) * t338;
t462 = qJD(1) * t342;
t435 = t338 * t462;
t268 = -qJD(4) * t435 + t310;
t275 = t336 * t513;
t448 = -qJDD(4) - qJDD(5);
t282 = t340 * t448 + qJDD(1);
t284 = -t336 * t340 + qJD(1);
t312 = -qJDD(4) * t340 + qJDD(1);
t315 = -qJD(4) * t340 + qJD(1);
t433 = t340 * t462;
t280 = t320 * t433;
t326 = qJDD(2) * t343;
t456 = qJD(3) * t343;
t431 = t338 * t456;
t460 = qJD(2) * t342;
t363 = t431 + t460;
t522 = qJ(3) * t338;
t365 = t522 + (pkin(2) - t320) * t340;
t349 = t343 * t365 - t447;
t182 = t303 + t349;
t550 = pkin(2) * t340;
t390 = t522 + t550;
t277 = t390 * t343;
t510 = t342 * qJ(2);
t296 = t343 * pkin(1) + t510;
t472 = -t277 - t296;
t438 = t182 + t472;
t467 = pkin(2) * t433 + qJ(3) * t435;
t502 = t343 * t337;
t514 = t338 * t342;
t468 = pkin(3) * t502 + t341 * t514;
t322 = pkin(1) * t462;
t416 = -t322 + t460;
t461 = qJD(1) * t343;
t271 = qJ(2) * t461 + t416;
t565 = -t431 + t467 - t271;
t346 = t326 + t438 * qJDD(1) + (-qJD(1) * t468 + t280 - t363 - t467 + t565) * qJD(1);
t417 = t338 ^ 2 * qJD(4) ^ 2 * t548;
t449 = qJDD(3) * t342;
t426 = t338 * t449;
t512 = t340 * t342;
t232 = t318 * t512 + t319 * t343;
t138 = qJD(1) * t232 - t235 * t336;
t233 = t340 * t509 - t504;
t139 = -qJD(1) * t233 - t234 * t336;
t393 = t139 * rSges(6,1) + t138 * rSges(6,2);
t80 = -rSges(6,3) * t435 + t393;
t503 = t343 * t323;
t507 = t342 * t324;
t253 = t340 * t503 - t507;
t356 = t253 * qJD(4);
t518 = t283 * t340;
t85 = t280 + (t338 * t464 - t518) * t462 + (t323 * t461 - t356) * pkin(4);
t12 = t312 * t110 - t282 * t134 + t268 * t186 + t194 * t208 + t275 * t201 - t284 * t80 - t315 * t85 - t343 * t417 + t346 - t426;
t588 = t12 * t343;
t453 = qJD(4) * t343;
t430 = t338 * t453;
t587 = t110 * t315 - t134 * t284 + t186 * t430 + t208 * t275;
t508 = t342 * t323;
t254 = t324 * t511 + t508;
t394 = t254 * rSges(5,1) - t253 * rSges(5,2);
t153 = rSges(5,3) * t513 + t394;
t213 = -rSges(5,3) * t340 + (rSges(5,1) * t324 - rSges(5,2) * t323) * t338;
t584 = -t153 * t315 + t213 * t430;
t124 = -Icges(6,5) * t233 + Icges(6,6) * t232 - Icges(6,3) * t514;
t125 = Icges(6,5) * t235 - Icges(6,6) * t234 + Icges(6,3) * t513;
t216 = Icges(6,4) * t235;
t128 = -Icges(6,2) * t234 + Icges(6,6) * t513 + t216;
t215 = Icges(6,4) * t234;
t132 = -Icges(6,1) * t235 - Icges(6,5) * t513 + t215;
t381 = t128 * t234 + t132 * t235;
t525 = Icges(6,4) * t233;
t127 = Icges(6,2) * t232 - Icges(6,6) * t514 - t525;
t214 = Icges(6,4) * t232;
t130 = -Icges(6,1) * t233 - Icges(6,5) * t514 + t214;
t501 = t232 * t127 - t233 * t130;
t583 = t381 + t501 + (-t124 * t342 - t125 * t343) * t338;
t274 = t336 * t514;
t47 = t124 * t513 - t234 * t127 + t235 * t130;
t205 = -Icges(6,3) * t340 + (Icges(6,5) * t319 - Icges(6,6) * t318) * t338;
t523 = Icges(6,4) * t319;
t206 = -Icges(6,6) * t340 + (-Icges(6,2) * t318 + t523) * t338;
t524 = Icges(6,4) * t318;
t207 = -Icges(6,5) * t340 + (Icges(6,1) * t319 - t524) * t338;
t68 = t205 * t513 - t206 * t234 + t207 * t235;
t582 = t274 * t47 - t284 * t68;
t144 = Icges(5,5) * t254 - Icges(5,6) * t253 + Icges(5,3) * t513;
t239 = Icges(5,4) * t254;
t147 = -Icges(5,2) * t253 + Icges(5,6) * t513 + t239;
t238 = Icges(5,4) * t253;
t151 = -Icges(5,1) * t254 - Icges(5,5) * t513 + t238;
t63 = -t144 * t340 - (t147 * t323 + t151 * t324) * t338;
t57 = -t125 * t340 - (t128 * t318 + t132 * t319) * t338;
t380 = -t147 * t253 - t151 * t254;
t251 = t324 * t343 + t340 * t508;
t252 = t340 * t507 - t503;
t528 = Icges(5,4) * t252;
t146 = Icges(5,2) * t251 - Icges(5,6) * t514 - t528;
t237 = Icges(5,4) * t251;
t149 = -Icges(5,1) * t252 - Icges(5,5) * t514 + t237;
t497 = -t251 * t146 + t252 * t149;
t579 = -t380 - t497;
t498 = t251 * t147 + t151 * t252;
t46 = -t125 * t514 + t232 * t128 + t132 * t233;
t505 = t342 * t339;
t564 = (t340 * t502 - t505) * rSges(4,2) - (t339 * t511 + t506) * rSges(4,1);
t572 = -rSges(4,3) * t513 + t564;
t475 = (-t340 * t505 + t502) * rSges(4,1) + (t339 * t343 + t340 * t506) * rSges(4,2);
t67 = -t205 * t514 + t206 * t232 - t207 * t233;
t567 = t275 * t46 + t67 * t284;
t162 = t232 * rSges(6,1) + t233 * rSges(6,2);
t434 = t338 * t461;
t563 = -t284 * t162 + t201 * t514 + t208 * t434 - t236 * t274;
t163 = -t234 * rSges(6,1) - t235 * rSges(6,2);
t562 = t163 * t284 + t201 * t513 - t275 * t236 + t340 * t80;
t143 = -Icges(5,5) * t252 + Icges(5,6) * t251 - Icges(5,3) * t514;
t52 = t143 * t513 - t253 * t146 + t254 * t149;
t353 = t342 * (Icges(5,2) * t252 + t149 + t237) - t343 * (-Icges(5,2) * t254 - t151 - t238);
t354 = t342 * (-Icges(5,1) * t251 + t146 - t528) - t343 * (Icges(5,1) * t253 + t147 + t239);
t230 = (-Icges(6,2) * t319 - t524) * t338;
t347 = t274 * (Icges(6,2) * t233 + t130 + t214) - t275 * (-Icges(6,2) * t235 - t132 - t215) - t284 * (t207 + t230);
t231 = (-Icges(6,1) * t318 - t523) * t338;
t348 = t274 * (-Icges(6,1) * t232 + t127 - t525) - t275 * (Icges(6,1) * t234 + t128 + t216) - t284 * (t206 - t231);
t561 = t340 ^ 2;
t193 = (t342 * t448 - t343 * t413) * t338;
t560 = t193 / 0.2e1;
t559 = t194 / 0.2e1;
t267 = (-qJD(1) * t453 - qJDD(4) * t342) * t338;
t558 = t267 / 0.2e1;
t557 = t268 / 0.2e1;
t556 = t274 / 0.2e1;
t555 = -t274 / 0.2e1;
t554 = -t275 / 0.2e1;
t553 = t275 / 0.2e1;
t551 = -t340 / 0.2e1;
t547 = g(2) * t343;
t229 = (-Icges(6,5) * t318 - Icges(6,6) * t319) * t338;
t198 = t336 * t229;
t199 = t336 * t230;
t200 = t336 * t231;
t58 = -t198 * t340 + ((-t206 * t336 + t200) * t319 + (-t207 * t336 - t199) * t318) * t338;
t84 = -t205 * t340 + (-t206 * t318 + t207 * t319) * t338;
t546 = t84 * t282 + t58 * t284;
t526 = Icges(5,4) * t324;
t211 = -Icges(5,6) * t340 + (-Icges(5,2) * t323 + t526) * t338;
t527 = Icges(5,4) * t323;
t212 = -Icges(5,5) * t340 + (Icges(5,1) * t324 - t527) * t338;
t248 = (-Icges(5,5) * t323 - Icges(5,6) * t324) * t338;
t225 = qJD(4) * t248;
t249 = (-Icges(5,2) * t324 - t527) * t338;
t226 = qJD(4) * t249;
t250 = (-Icges(5,1) * t323 - t526) * t338;
t227 = qJD(4) * t250;
t65 = -t225 * t340 + (-t226 * t323 + t227 * t324 + (-t211 * t324 - t212 * t323) * qJD(4)) * t338;
t210 = -Icges(5,3) * t340 + (Icges(5,5) * t324 - Icges(5,6) * t323) * t338;
t96 = -t210 * t340 + (-t211 * t323 + t212 * t324) * t338;
t545 = t96 * t312 + t65 * t315;
t140 = qJD(1) * t234 + t233 * t336;
t141 = -qJD(1) * t235 + t232 * t336;
t492 = t141 * rSges(6,1) + t140 * rSges(6,2);
t81 = -rSges(6,3) * t434 + t492;
t292 = t341 * t434;
t421 = t340 * t470;
t541 = pkin(4) * qJD(4);
t414 = t340 * t323 * t541;
t443 = t324 * t541;
t439 = t333 * t434 + t342 * t414 + t343 * t443;
t86 = -t292 + (-pkin(4) * t508 - t343 * t421) * qJD(1) + t439;
t544 = -t81 - t86;
t543 = rSges(3,1) * t340;
t542 = rSges(6,3) * t338;
t75 = Icges(6,5) * t141 + Icges(6,6) * t140 - Icges(6,3) * t434;
t77 = Icges(6,4) * t141 + Icges(6,2) * t140 - Icges(6,6) * t434;
t79 = Icges(6,1) * t141 + Icges(6,4) * t140 - Icges(6,5) * t434;
t29 = -t340 * t75 + ((-t127 * t336 + t79) * t319 + (-t130 * t336 - t77) * t318) * t338;
t538 = t29 * t274;
t74 = Icges(6,5) * t139 + Icges(6,6) * t138 - Icges(6,3) * t435;
t76 = Icges(6,4) * t139 + Icges(6,2) * t138 - Icges(6,6) * t435;
t78 = Icges(6,1) * t139 + Icges(6,4) * t138 - Icges(6,5) * t435;
t30 = -t340 * t74 + ((-t128 * t336 + t78) * t319 + (t132 * t336 - t76) * t318) * t338;
t537 = t30 * t275;
t73 = t210 * t513 - t211 * t253 + t212 * t254;
t536 = t315 * t73;
t56 = -t124 * t340 + (-t127 * t318 + t130 * t319) * t338;
t535 = t56 * t193;
t534 = t57 * t194;
t62 = -t143 * t340 + (-t146 * t323 + t149 * t324) * t338;
t533 = t62 * t267;
t532 = t63 * t268;
t531 = -rSges(3,3) - qJ(2);
t530 = -rSges(4,3) - qJ(3);
t519 = t143 * t342;
t471 = t333 * t514 + t343 * t404;
t109 = -t342 * t421 - t468 + t471;
t481 = -t233 * rSges(6,1) + t232 * rSges(6,2);
t133 = -rSges(6,3) * t514 + t481;
t500 = -t109 - t133;
t499 = t110 - t134;
t171 = qJD(1) * t253 + qJD(4) * t252;
t172 = -qJD(1) * t254 + qJD(4) * t251;
t487 = t172 * rSges(5,1) + t171 * rSges(5,2);
t483 = t211 - t250;
t482 = t212 + t249;
t479 = -t252 * rSges(5,1) + t251 * rSges(5,2);
t477 = t564 * qJD(1);
t276 = t390 * t342;
t330 = t343 * qJ(2);
t294 = -pkin(1) * t342 + t330;
t473 = -t276 + t294;
t240 = t251 * pkin(4);
t314 = rSges(3,2) * t513;
t374 = -rSges(3,1) * t511 - rSges(3,3) * t342;
t257 = -t314 - t374;
t469 = -t296 - t257;
t466 = rSges(3,2) * t514 + t343 * rSges(3,3);
t289 = qJD(1) * t296;
t328 = qJD(2) * t343;
t465 = t328 - t289;
t463 = qJD(1) * t338;
t459 = qJD(3) * t338;
t458 = qJD(3) * t340;
t457 = qJD(3) * t342;
t455 = qJD(4) * t338;
t454 = qJD(4) * t342;
t452 = -m(4) - m(5) - m(6);
t451 = qJD(1) * qJD(3);
t450 = qJDD(3) * t340;
t437 = t472 + t572;
t436 = -pkin(1) - t550;
t432 = t338 * t457;
t429 = -t514 / 0.2e1;
t428 = t513 / 0.2e1;
t427 = -pkin(1) - t543;
t425 = -t463 / 0.2e1;
t424 = -t455 / 0.2e1;
t423 = t455 / 0.2e1;
t422 = -qJ(2) - t549;
t420 = -t162 * t275 - t274 * t163;
t412 = qJDD(1) * t294 + qJDD(2) * t342 + (t328 + t465) * qJD(1);
t411 = t342 * t425;
t410 = t343 * t425;
t409 = t342 * t424;
t408 = t342 * t423;
t407 = t343 * t424;
t406 = t343 * t423;
t402 = t328 - t432;
t297 = rSges(2,1) * t343 - t342 * rSges(2,2);
t399 = rSges(2,1) * t342 + rSges(2,2) * t343;
t398 = -rSges(3,2) * t338 + t543;
t397 = t475 * qJD(1);
t169 = qJD(1) * t251 - qJD(4) * t254;
t170 = -qJD(1) * t252 - t356;
t395 = rSges(5,1) * t170 + rSges(5,2) * t169;
t88 = Icges(5,5) * t170 + Icges(5,6) * t169 - Icges(5,3) * t435;
t89 = Icges(5,5) * t172 + Icges(5,6) * t171 - Icges(5,3) * t434;
t90 = Icges(5,4) * t170 + Icges(5,2) * t169 - Icges(5,6) * t435;
t91 = Icges(5,4) * t172 + Icges(5,2) * t171 - Icges(5,6) * t434;
t92 = Icges(5,1) * t170 + Icges(5,4) * t169 - Icges(5,5) * t435;
t93 = Icges(5,1) * t172 + Icges(5,4) * t171 - Icges(5,5) * t434;
t389 = -(t146 * t169 + t149 * t170 - t253 * t91 + t254 * t93 + (-t143 * t462 + t343 * t89) * t338) * t342 + (t147 * t169 - t151 * t170 - t253 * t90 + t254 * t92 + (-t144 * t462 + t343 * t88) * t338) * t343;
t388 = -(t146 * t171 + t149 * t172 + t251 * t91 - t252 * t93 + (-t143 * t461 - t342 * t89) * t338) * t342 + (t147 * t171 - t151 * t172 + t251 * t90 - t252 * t92 + (-t144 * t461 - t342 * t88) * t338) * t343;
t32 = -t340 * t89 + (-t323 * t91 + t324 * t93 + (-t146 * t324 - t149 * t323) * qJD(4)) * t338;
t33 = -t340 * t88 + (-t323 * t90 + t324 * t92 + (-t147 * t324 + t151 * t323) * qJD(4)) * t338;
t387 = -t32 * t342 + t33 * t343;
t152 = -rSges(5,3) * t514 + t479;
t255 = (-rSges(5,1) * t323 - rSges(5,2) * t324) * t338;
t228 = qJD(4) * t255;
t181 = t342 * t365 + t468;
t359 = qJD(1) * (-t390 * t461 - t432) - qJDD(1) * t276 + qJDD(3) * t513 + t412;
t352 = qJD(1) * (qJD(1) * t349 + t292) + qJDD(1) * t181 + t359;
t95 = -rSges(5,3) * t434 + t487;
t34 = t152 * t312 - t213 * t267 + t315 * t95 + (qJD(4) * t228 - t451) * t514 + t352;
t94 = -rSges(5,3) * t435 + t395;
t35 = -t153 * t312 + t213 * t268 - t315 * t94 + (t228 * t453 - t449) * t338 + t346;
t386 = t34 * t342 + t35 * t343;
t50 = -t143 * t514 - t497;
t51 = -t144 * t514 + t498;
t385 = -t51 * t342 - t50 * t343;
t355 = t460 + (t181 + t473) * qJD(1);
t60 = t152 * t315 + (t213 * t454 + t456) * t338 + t355;
t351 = qJD(1) * t438 + t402;
t61 = t351 + t584;
t384 = t60 * t342 + t61 * t343;
t383 = -t342 * t94 - t343 * t95;
t382 = -t518 - t542;
t379 = -t152 * t343 - t153 * t342;
t378 = (Icges(5,5) * t251 + Icges(5,6) * t252) * t342 - (-Icges(5,5) * t253 - Icges(5,6) * t254) * t343;
t241 = t253 * pkin(4);
t377 = -qJ(2) - t404;
t376 = -rSges(5,3) * t338 - pkin(1) - t517;
t375 = -pkin(1) + t382;
t367 = (-t342 * t50 + t343 * t51) * t338;
t53 = t144 * t513 + t380;
t366 = (-t342 * t52 + t343 * t53) * t338;
t42 = -t458 - t133 * t275 - t134 * t274 + (-t109 * t343 + t110 * t342) * t455;
t364 = t12 * (t340 * t134 + t208 * t513) + t42 * t133 * t435;
t362 = -qJD(1) * t277 - t289 + t402;
t360 = t530 * t338 + t436;
t358 = -(Icges(6,5) * t232 + Icges(6,6) * t233) * t274 + (-Icges(6,5) * t234 - Icges(6,6) * t235) * t275 + t229 * t284;
t357 = qJD(1) * t182 + t362;
t13 = t127 * t138 + t130 * t139 - t234 * t77 + t235 * t79 + (-t124 * t462 + t343 * t75) * t338;
t14 = t128 * t138 - t132 * t139 - t234 * t76 + t235 * t78 + (-t125 * t462 + t343 * t74) * t338;
t15 = t127 * t140 + t130 * t141 + t232 * t77 - t233 * t79 + (-t124 * t461 - t342 * t75) * t338;
t16 = t128 * t140 - t132 * t141 + t232 * t76 - t233 * t78 + (-t125 * t461 - t342 * t74) * t338;
t45 = -t124 * t514 + t501;
t21 = -t274 * t45 + t567;
t48 = t125 * t513 - t381;
t22 = t275 * t48 - t582;
t37 = t138 * t206 + t139 * t207 - t199 * t234 + t200 * t235 + (t198 * t343 - t205 * t462) * t338;
t38 = t140 * t206 + t141 * t207 + t199 * t232 - t200 * t233 + (-t198 * t342 - t205 * t461) * t338;
t350 = (-t15 * t274 + t16 * t275 + t193 * t45 + t194 * t46 + t282 * t67 + t284 * t38) * t429 - (-t358 * t340 + (t318 * t347 + t319 * t348) * t338) * t284 / 0.2e1 + t22 * t411 + t21 * t410 + (-t340 * t67 + (-t342 * t45 + t343 * t46) * t338) * t560 + (-t13 * t274 + t14 * t275 + t193 * t47 + t194 * t48 + t282 * t68 + t284 * t37) * t428 + (-t340 * t68 + (-t342 * t47 + t343 * t48) * t338) * t559 + (-t340 * t38 + (-t15 * t342 + t16 * t343 + (-t342 * t46 - t343 * t45) * qJD(1)) * t338) * t555 + t282 * (-t340 * t84 + (-t342 * t56 + t343 * t57) * t338) / 0.2e1 + (-t340 * t37 + (-t13 * t342 + t14 * t343 + (-t342 * t48 - t343 * t47) * qJD(1)) * t338) * t553 + (t534 + t535 + t537 - t538 + t546) * t551 + t284 * (-t340 * t58 + (-t29 * t342 + t30 * t343 + (-t342 * t57 - t343 * t56) * qJD(1)) * t338) / 0.2e1 + (-t232 * t347 - t233 * t348 - t358 * t514) * t556 + (t234 * t347 + t235 * t348 + t358 * t513) * t554;
t306 = rSges(3,2) * t434;
t256 = -rSges(3,1) * t512 + t466;
t195 = t208 * t514;
t189 = qJD(1) * t469 + t328;
t188 = t460 + (t256 + t294) * qJD(1);
t183 = -rSges(4,3) * t514 + t475;
t180 = -rSges(5,1) * t253 - rSges(5,2) * t254;
t179 = rSges(5,1) * t251 + rSges(5,2) * t252;
t102 = qJD(1) * t437 + t402;
t101 = (t183 + t473) * qJD(1) + t363;
t99 = t326 + t469 * qJDD(1) + (-rSges(3,3) * t461 - t271 + (qJD(1) * t398 - qJD(2)) * t342) * qJD(1);
t98 = qJDD(1) * t256 + (qJD(1) * t374 + t306) * qJD(1) + t412;
t72 = -t210 * t514 + t211 * t251 - t212 * t252;
t71 = t379 * t455 - t458;
t69 = t72 * t315;
t55 = -t426 + t326 + t437 * qJDD(1) + (-t460 + (rSges(4,3) * t462 - t456) * t338 - t397 + t565) * qJD(1);
t54 = qJDD(1) * t183 + ((-rSges(4,3) * t461 - t457) * t338 + t477) * qJD(1) + t359;
t44 = t171 * t211 + t172 * t212 + t226 * t251 - t227 * t252 + (-t210 * t461 - t225 * t342) * t338;
t43 = t169 * t211 + t170 * t212 - t226 * t253 + t227 * t254 + (-t210 * t462 + t225 * t343) * t338;
t41 = t351 + t587;
t40 = t109 * t315 + t133 * t284 + t208 * t274 + (t186 * t454 + t456) * t338 + t355;
t36 = -t152 * t268 + t153 * t267 + t383 * t455 - t450;
t28 = qJD(4) * t366 + t536;
t27 = qJD(4) * t367 + t69;
t11 = t312 * t109 + t315 * t86 + t282 * t133 + t284 * t81 - t267 * t186 + t274 * t201 + t352 - t193 * t208 + (-t338 * t451 - t417) * t342;
t9 = -t450 - t109 * t268 - t110 * t267 - t133 * t194 + t134 * t193 - t274 * t80 - t275 * t81 + (-t342 * t85 - t343 * t86) * t455;
t1 = [(g(2) * t297 + g(3) * t399) * m(2) + (m(2) * (t297 ^ 2 + t399 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t561 + ((Icges(3,1) + Icges(4,1) * t339 ^ 2 + (-0.2e1 * Icges(4,4) * t339 + Icges(4,2) * t337) * t337) * t338 + 0.2e1 * (-Icges(4,5) * t339 + Icges(4,6) * t337 + Icges(3,4)) * t340) * t338) * qJDD(1) + (-(-pkin(1) + (-rSges(6,3) + t333) * t338) * t547 + (-pkin(1) + t515 - t542) * t588 + (-g(3) + t11) * (t342 * t375 + t330 + t471 + t481) + (-g(2) + t12) * (-t392 + t476 - t510) + (-t342 * t443 - t393 - t416 + (t414 - t459) * t343 + (t377 * t343 + (-t382 - t515) * t342) * qJD(1)) * t41 + (-t357 + t41 + t402 + t439 + t492 + (t377 * t342 + t375 * t343) * qJD(1) - t587) * t40) * m(6) + (t69 + (t498 * t343 + (t338 * t519 - t53 - t579) * t342) * t455) * t407 + (-t536 + (-(-t498 - t52) * t342 + t579 * t343 + (-t144 * t342 ^ 2 + (-t144 * t343 - t519) * t343) * t338 + t385) * t455 + t28) * t408 + (-(qJD(1) * t572 - t102 + t362) * t101 - g(2) * (t564 - t510) - t360 * t547 + t55 * t564 + t102 * (t322 - t397 + t467) + t101 * (t328 + t477) + (-t55 * qJ(2) + t102 * (rSges(4,3) * t463 - qJD(2)) + t101 * (-qJD(1) * qJ(2) - t459)) * t342 + (t55 * t436 + (-t102 * qJD(3) + t530 * t55) * t338 + (-t102 * qJ(2) + t101 * (-rSges(4,3) * t338 - pkin(1) - t390)) * qJD(1)) * t343 + (-g(3) + t54) * (t342 * t360 + t330 + t475)) * m(4) + (-(-qJD(1) * t257 - t189 + t465) * t188 - t189 * t416 + t188 * (t306 + t328) + ((t188 * t427 + t189 * t531) * t343 + (t188 * t531 + t189 * t398) * t342) * qJD(1) + (t99 - g(2)) * (t531 * t342 + t427 * t343 + t314) + (t98 - g(3)) * (t342 * t427 + t330 + t466)) * m(3) + ((-t45 + t583) * t275 + t22 + t582) * t556 + (-(t48 + t583) * t274 + t567) * t554 + t545 + t546 + (-t376 * t547 - (t357 + t584 - t61) * t60 + t61 * (t280 - t395 - t416) + t60 * (t292 + t402 + t487) + (t35 * t376 - t61 * t459) * t343 + ((t60 * t422 + t61 * (rSges(5,3) - t341) * t338) * t342 + (t376 * t60 + t422 * t61) * t343) * qJD(1) + (-g(3) + t34) * (t342 * t376 + t330 + t468 + t479) + (-g(2) + t35) * (-t405 - t510 - t394)) * m(5) - t538 / 0.2e1 + t537 / 0.2e1 + t532 / 0.2e1 + t533 / 0.2e1 + (t37 + t21) * t553 + (t33 + t43 + t27) * t406 + t38 * t555 + t73 * t557 + t72 * t558 + t68 * t559 + t67 * t560 + t535 / 0.2e1 + t534 / 0.2e1 + (t32 + t44) * t409; (-m(3) + t452) * (g(3) * t342 + t547) + m(3) * (t342 * t98 + t343 * t99) + m(4) * (t342 * t54 + t343 * t55) + m(5) * t386 + m(6) * (t11 * t342 + t588); t452 * (-g(1) * t340 + (-g(2) * t342 + g(3) * t343) * t338) + m(4) * (qJDD(3) * t561 + t513 * t54 - t514 * t55) + m(5) * (t34 * t513 - t340 * t36 - t35 * t514) + m(6) * (t11 * t513 - t12 * t514 - t340 * t9); t350 + t315 * (-t340 * t65 + ((-t63 * t342 - t343 * t62) * qJD(1) + t387) * t338) / 0.2e1 + (-t340 * t73 + t366) * t557 + t27 * t410 + t28 * t411 + (t267 * t50 + t268 * t51 + t312 * t72 + t315 * t44 + t388 * t455) * t429 + ((t248 * t513 - t253 * t482 - t254 * t483) * t315 + (t253 * t353 + t254 * t354 - t378 * t513) * t455) * t407 + ((-t248 * t514 + t251 * t482 + t252 * t483) * t315 + (-t251 * t353 - t252 * t354 + t378 * t514) * t455) * t408 + t312 * (-t340 * t96 + (-t342 * t62 + t343 * t63) * t338) / 0.2e1 - t315 * (-t340 * t248 * t315 + ((-t323 * t482 - t324 * t483) * t315 + ((t323 * t353 + t324 * t354) * t338 + t378 * t340) * qJD(4)) * t338) / 0.2e1 + (t387 * t455 + t532 + t533 + t545) * t551 + (-t340 * t43 + ((-t53 * t342 - t52 * t343) * qJD(1) + t389) * t338) * t406 + (t267 * t52 + t268 * t53 + t312 * t73 + t315 * t43 + t389 * t455) * t428 + (-t340 * t72 + t367) * t558 + (-t340 * t44 + (qJD(1) * t385 + t388) * t338) * t409 + (-g(2) * (t240 + t162) - g(3) * (-t241 + t163) + t11 * t195 + (t11 * t500 - t12 * t110) * t340 + t364 - t42 * ((-t240 * t343 + t241 * t342) * t455 + t420) + (-t241 * t315 + t340 * t85 + t562) * t41 + (-t240 * t315 + t340 * t544 + t563) * t40 + (-g(1) * (t391 - t548) + (t9 * t500 + t42 * t544 + t12 * t186 + (t40 * t186 + t42 * t499) * qJD(1)) * t343 + (t9 * t499 + t42 * (-t80 - t85) + t11 * t186 + (t42 * t109 + t41 * (-t186 - t208)) * qJD(1)) * t342) * t338) * m(6) + (-g(1) * t255 - g(2) * t179 - g(3) * t180 + (-t152 * t34 + t153 * t35 - t60 * t95 + t61 * t94) * t340 + (t36 * t379 + t71 * (t152 * t462 - t153 * t461 + t383) + t384 * t228 + ((-t61 * t342 + t60 * t343) * qJD(1) + t386) * t213) * t338 - (t179 * t60 - t180 * t61) * t315 - (t71 * (-t179 * t343 - t180 * t342) + t384 * t255) * t455) * m(5); t350 + (-g(1) * t236 - g(2) * t162 - g(3) * t163 + t11 * (-t133 * t340 + t195) + (t9 * (-t133 * t343 - t134 * t342) + t42 * (-t134 * t461 - t342 * t80 - t343 * t81)) * t338 + t364 - t42 * t420 + (-t208 * t435 + t562) * t41 + (-t340 * t81 + t563) * t40) * m(6);];
tau = t1;
