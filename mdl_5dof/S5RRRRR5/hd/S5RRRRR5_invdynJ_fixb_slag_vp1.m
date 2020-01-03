% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:03
% EndTime: 2020-01-03 12:13:19
% DurationCPUTime: 10.76s
% Computational Cost: add. (23835->652), mult. (14122->816), div. (0->0), fcn. (10979->10), ass. (0->376)
t316 = sin(qJ(1));
t504 = pkin(1) * qJD(1);
t420 = t316 * t504;
t314 = qJ(1) + qJ(2);
t300 = sin(t314);
t312 = qJD(1) + qJD(2);
t460 = t300 * t312;
t426 = pkin(2) * t460;
t354 = t420 + t426;
t298 = qJD(3) + t312;
t304 = qJ(3) + t314;
t289 = sin(t304);
t290 = cos(t304);
t546 = t289 * rSges(4,1) + t290 * rSges(4,2);
t560 = t546 * t298;
t143 = t354 + t560;
t315 = sin(qJ(4));
t317 = cos(qJ(4));
t272 = rSges(5,1) * t315 + rSges(5,2) * t317;
t197 = t272 * t289;
t198 = t272 * t290;
t469 = t289 * t315;
t422 = rSges(5,2) * t469;
t468 = t289 * t317;
t385 = rSges(5,1) * t468 - t422;
t167 = -rSges(5,3) * t290 + t385;
t280 = t289 * pkin(3);
t216 = -pkin(8) * t290 + t280;
t431 = qJD(4) * t290;
t347 = -t272 * t431 - t298 * (t167 + t216);
t77 = -t347 + t354;
t432 = qJD(4) * t289;
t412 = t272 * t432;
t302 = cos(t314);
t458 = t302 * t312;
t271 = pkin(2) * t458;
t318 = cos(qJ(1));
t296 = t318 * t504;
t434 = t271 + t296;
t464 = t290 * t315;
t259 = rSges(5,2) * t464;
t463 = t290 * t317;
t425 = rSges(5,1) * t463;
t168 = rSges(5,3) * t289 - t259 + t425;
t508 = pkin(3) * t290;
t217 = pkin(8) * t289 + t508;
t445 = t168 + t217;
t78 = t298 * t445 - t412 + t434;
t559 = -t197 * t77 - t198 * t78;
t325 = t559 * qJD(4);
t147 = t298 * t168;
t199 = t298 * t217;
t416 = t298 * t464;
t467 = t290 * t298;
t472 = t289 * t298;
t438 = pkin(3) * t467 + pkin(8) * t472;
t441 = rSges(5,3) * t472 + t298 * t425;
t540 = -rSges(5,2) * t416 - t147 - t199 + t412 + t438 + t441;
t566 = t540 * t77 + t325;
t544 = t300 * rSges(3,1) + t302 * rSges(3,2);
t204 = t544 * t312;
t183 = t420 + t204;
t313 = qJ(4) + qJ(5);
t299 = sin(t313);
t301 = cos(t313);
t226 = rSges(6,1) * t299 + rSges(6,2) * t301;
t285 = t301 * rSges(6,1);
t307 = t317 * pkin(4);
t291 = t307 + pkin(3);
t311 = qJD(4) + qJD(5);
t319 = -pkin(9) - pkin(8);
t429 = qJD(4) * t315;
t419 = pkin(4) * t429;
t362 = -t298 * t319 - t419;
t462 = t298 * t299;
t421 = rSges(6,2) * t462;
t461 = t299 * t311;
t423 = rSges(6,1) * t461;
t459 = t301 * t311;
t213 = t290 * t311;
t408 = t290 * t429;
t390 = pkin(4) * t408;
t366 = -t213 * t226 - t390;
t336 = -t366 + t426;
t276 = t290 * t319;
t437 = t289 * t291 + t276;
t145 = -t216 + t437;
t470 = t289 * t301;
t471 = t289 * t299;
t155 = rSges(6,1) * t470 - rSges(6,2) * t471 - rSges(6,3) * t290;
t456 = t145 + t155;
t415 = t216 + t456;
t58 = t298 * t415 + t336 + t420;
t212 = t289 * t311;
t410 = t289 * t429;
t352 = -pkin(4) * t410 - t212 * t226;
t249 = t290 * t291;
t146 = t508 - t249 + (pkin(8) + t319) * t289;
t465 = t290 * t301;
t424 = rSges(6,1) * t465;
t466 = t290 * t299;
t386 = -rSges(6,2) * t466 + t424;
t156 = rSges(6,3) * t289 + t386;
t455 = -t146 + t156;
t414 = t217 + t455;
t59 = t298 * t414 + t352 + t434;
t322 = (t58 * (-t226 * t311 - t419) + (t59 * (-t291 - t285) - t58 * t319) * t298) * t289 + (t59 * (-rSges(6,2) * t459 + t362 - t423) - t58 * t421) * t290;
t502 = t298 * t59;
t142 = t298 * t156;
t209 = t291 * t467;
t443 = rSges(6,3) * t472 + t298 * t424;
t541 = t298 * t146 - t142 - t199 + t209 - t352 + t443;
t565 = t415 * t502 + t541 * t58 + t322;
t497 = Icges(6,4) * t299;
t225 = Icges(6,1) * t301 - t497;
t359 = t225 * t289;
t153 = -Icges(6,5) * t290 + t359;
t247 = Icges(6,4) * t466;
t154 = Icges(6,1) * t465 + Icges(6,5) * t289 - t247;
t222 = Icges(6,2) * t301 + t497;
t282 = Icges(6,4) * t301;
t376 = -Icges(6,2) * t299 + t282;
t545 = Icges(6,1) * t299 + t282;
t558 = t545 + t376;
t327 = t212 * (-Icges(6,2) * t465 + t154 - t247) - t213 * (-t222 * t289 + t153) + t298 * t558;
t357 = t376 * t289;
t151 = -Icges(6,6) * t290 + t357;
t152 = Icges(6,4) * t465 - Icges(6,2) * t466 + Icges(6,6) * t289;
t525 = t212 * (t290 * t545 + t152) - t213 * (t289 * t545 + t151) + t298 * (t222 - t225);
t564 = t327 * t299 + t301 * t525;
t303 = Icges(5,4) * t317;
t377 = -Icges(5,2) * t315 + t303;
t543 = Icges(5,1) * t315 + t303;
t435 = t543 + t377;
t498 = Icges(5,4) * t315;
t265 = Icges(5,2) * t317 + t498;
t268 = Icges(5,1) * t317 - t498;
t436 = t265 - t268;
t561 = (t315 * t435 + t317 * t436) * t298;
t427 = -qJDD(4) - qJDD(5);
t134 = -t213 * t298 + t289 * t427;
t243 = pkin(8) * t467;
t173 = pkin(3) * t472 - t243;
t430 = qJD(4) * t298;
t187 = -qJDD(4) * t289 - t290 * t430;
t228 = -rSges(6,2) * t299 + t285;
t203 = t228 * t311;
t310 = qJDD(1) + qJDD(2);
t297 = qJDD(3) + t310;
t306 = t316 * pkin(1);
t321 = qJD(1) ^ 2;
t488 = pkin(1) * qJDD(1);
t387 = -t306 * t321 + t318 * t488;
t509 = pkin(2) * t310;
t510 = pkin(2) * t312 ^ 2;
t337 = -t300 * t510 + t302 * t509 + t387;
t457 = t317 * qJD(4) ^ 2;
t253 = rSges(6,2) * t465;
t444 = rSges(6,3) * t467 + t289 * t421;
t90 = t311 * t253 + (t290 * t461 + t298 * t470) * rSges(6,1) - t444;
t98 = t390 + t243 + (t276 + (-pkin(3) + t291) * t289) * t298;
t507 = -t90 - t98;
t18 = t134 * t226 - t203 * t212 + (t187 * t315 - t289 * t457) * pkin(4) + (-t173 + t507) * t298 + t414 * t297 + t337;
t557 = -g(2) + t18;
t417 = t298 * t468;
t428 = qJD(4) * t317;
t442 = -rSges(5,3) * t467 - t298 * t422;
t107 = rSges(5,2) * t290 * t428 + (t408 + t417) * rSges(5,1) + t442;
t506 = rSges(5,1) * t317;
t274 = -rSges(5,2) * t315 + t506;
t239 = t274 * qJD(4);
t45 = -t239 * t432 + t187 * t272 + (-t107 - t173) * t298 + t445 * t297 + t337;
t556 = -g(2) + t45;
t230 = t289 * t430;
t135 = qJD(5) * t472 + t290 * t427 + t230;
t188 = -qJDD(4) * t290 + t230;
t308 = t318 * pkin(1);
t433 = t321 * t308 + t316 * t488;
t388 = t300 * t509 + t302 * t510 + t433;
t363 = t297 * t216 + t298 * t438 + t388;
t91 = -t289 * t423 + (-t289 * t459 - t290 * t462) * rSges(6,2) + t443;
t99 = t289 * t362 + t209 - t438;
t17 = -t135 * t226 + t203 * t213 + (t91 + t99) * t298 + t456 * t297 + (-t188 * t315 + t290 * t457) * pkin(4) + t363;
t555 = -g(3) + t17;
t409 = t289 * t428;
t108 = -rSges(5,1) * t410 + (-t409 - t416) * rSges(5,2) + t441;
t44 = t108 * t298 + t167 * t297 - t188 * t272 + t239 * t431 + t363;
t554 = -g(3) + t44;
t172 = rSges(4,1) * t467 - rSges(4,2) * t472;
t553 = t172 * t298 + t297 * t546 - g(3) + t388;
t215 = rSges(4,1) * t290 - t289 * rSges(4,2);
t552 = t215 * t297 - t298 * t560 - g(2) + t337;
t205 = rSges(3,1) * t458 - rSges(3,2) * t460;
t551 = t205 * t312 + t310 * t544 - g(3) + t433;
t229 = rSges(3,1) * t302 - t300 * rSges(3,2);
t550 = -t204 * t312 + t229 * t310 - g(2) + t387;
t221 = Icges(6,5) * t301 - Icges(6,6) * t299;
t149 = -Icges(6,3) * t290 + t221 * t289;
t483 = t152 * t299;
t375 = -t154 * t301 + t483;
t361 = -t149 + t375;
t548 = t213 * t361;
t196 = t298 * t215;
t144 = t434 + t196;
t475 = t222 * t311;
t539 = -Icges(6,6) * t298 + t475;
t180 = t226 * t289;
t181 = rSges(6,1) * t466 + t253;
t52 = t155 * t212 + t156 * t213 + (t145 * t289 - t146 * t290) * qJD(4);
t538 = -t58 * (-t298 * t180 + t213 * t228) - t52 * (-t180 * t212 - t213 * t181) - t59 * (-t181 * t298 - t212 * t228);
t220 = Icges(6,5) * t299 + Icges(6,6) * t301;
t537 = -Icges(6,3) * t298 + t220 * t311;
t536 = -Icges(6,5) * t298 + t311 * t545;
t533 = -Icges(5,6) * t298 + qJD(4) * t265;
t103 = -t289 * t533 + t377 * t467;
t530 = -Icges(5,5) * t298 + qJD(4) * t543;
t105 = t268 * t467 - t289 * t530;
t358 = t377 * t289;
t163 = -Icges(5,6) * t290 + t358;
t360 = t268 * t289;
t165 = -Icges(5,5) * t290 + t360;
t109 = t163 * t317 + t165 * t315;
t264 = Icges(5,5) * t317 - Icges(5,6) * t315;
t356 = t264 * t289;
t161 = -Icges(5,3) * t290 + t356;
t535 = qJD(4) * t109 + t103 * t315 - t105 * t317 - t161 * t298;
t263 = Icges(5,5) * t315 + Icges(5,6) * t317;
t534 = -Icges(5,3) * t298 + qJD(4) * t263;
t232 = t377 * qJD(4);
t233 = t268 * qJD(4);
t369 = t265 * t317 + t315 * t543;
t532 = qJD(4) * t369 + t232 * t315 - t233 * t317 - t263 * t298;
t102 = t290 * t533 + t298 * t358;
t104 = t290 * t530 + t298 * t360;
t162 = Icges(5,5) * t463 - Icges(5,6) * t464 + Icges(5,3) * t289;
t164 = Icges(5,4) * t463 - Icges(5,2) * t464 + Icges(5,6) * t289;
t257 = Icges(5,4) * t464;
t166 = Icges(5,1) * t463 + Icges(5,5) * t289 - t257;
t373 = t164 * t317 + t315 * t166;
t531 = qJD(4) * t373 - t102 * t315 + t104 * t317 - t162 * t298;
t391 = -t225 * t311 + t475;
t392 = t558 * t311;
t528 = -t220 * t298 + t299 * t392 + t301 * t391;
t150 = Icges(6,5) * t465 - Icges(6,6) * t466 + Icges(6,3) * t289;
t396 = t154 * t311 - t290 * t539 - t298 * t357;
t398 = t152 * t311 + t290 * t536 + t298 * t359;
t527 = -t150 * t298 + t299 * t396 + t301 * t398;
t397 = t153 * t311 - t289 * t539 + t376 * t467;
t399 = t151 * t311 - t225 * t467 + t289 * t536;
t526 = -t149 * t298 + t299 * t397 + t301 * t399;
t524 = t134 / 0.2e1;
t523 = t135 / 0.2e1;
t522 = t187 / 0.2e1;
t521 = t188 / 0.2e1;
t520 = -t212 / 0.2e1;
t519 = t212 / 0.2e1;
t518 = t213 / 0.2e1;
t517 = -t213 / 0.2e1;
t516 = -t289 / 0.2e1;
t515 = -t290 / 0.2e1;
t514 = t297 / 0.2e1;
t513 = -t298 / 0.2e1;
t512 = t298 / 0.2e1;
t511 = rSges(5,3) + pkin(8);
t477 = t220 * t289;
t97 = -t222 * t466 + t465 * t545 + t477;
t501 = t97 * t298;
t474 = t263 * t289;
t117 = -t265 * t464 + t463 * t543 + t474;
t487 = t117 * t298;
t484 = t151 * t299;
t482 = t153 * t301;
t481 = t163 * t315;
t480 = t164 * t315;
t479 = t165 * t317;
t175 = t220 * t290;
t355 = t221 * t298;
t191 = t263 * t290;
t473 = t264 * t298;
t450 = t289 * t543 + t163;
t449 = t290 * t543 + t164;
t448 = -t265 * t289 + t165;
t447 = -Icges(5,2) * t463 + t166 - t257;
t418 = t155 * t467 + (-t142 + t91) * t289;
t287 = pkin(2) * t300;
t185 = t287 + t546;
t406 = t472 / 0.2e1;
t405 = -t467 / 0.2e1;
t404 = -t432 / 0.2e1;
t403 = t432 / 0.2e1;
t402 = -t431 / 0.2e1;
t401 = t431 / 0.2e1;
t400 = -pkin(4) * t315 - t226;
t395 = -t150 - t484;
t394 = -t150 + t482;
t184 = t229 * t312 + t296;
t288 = pkin(2) * t302;
t186 = t215 + t288;
t275 = rSges(2,1) * t318 - t316 * rSges(2,2);
t273 = rSges(2,1) * t316 + rSges(2,2) * t318;
t139 = t165 * t468;
t70 = -t161 * t290 - t163 * t469 + t139;
t140 = t166 * t468;
t71 = t162 * t290 + t164 * t469 - t140;
t382 = -t289 * t71 - t290 * t70;
t141 = t163 * t464;
t72 = -t161 * t289 - t165 * t463 + t141;
t372 = -t166 * t317 + t480;
t73 = t289 * t162 - t372 * t290;
t381 = -t289 * t73 - t290 * t72;
t380 = -t289 * t78 + t290 * t77;
t93 = -t152 * t301 - t154 * t299;
t374 = t479 - t481;
t371 = t167 * t289 + t168 * t290;
t370 = -t222 * t299 + t301 * t545;
t368 = -t265 * t315 + t317 * t543;
t365 = t172 + t271;
t349 = t175 * t212 - t213 * t477 - t355;
t346 = -t289 * t355 - t290 * t537 + t298 * t375;
t345 = -t290 * t355 + t537 * t289 + (t482 - t484) * t298;
t344 = -t290 * t534 + (-t356 + t372) * t298;
t343 = -t264 * t467 + t289 * t534 + t298 * t374;
t342 = -t221 * t311 + t298 * t370;
t341 = -t264 * qJD(4) + t298 * t368;
t340 = t315 * t448 + t317 * t450;
t339 = t315 * t447 + t317 * t449;
t122 = t155 + t437;
t132 = -t290 * t511 + t280 + t385;
t123 = t249 + (rSges(6,3) - t319) * t289 + t386;
t114 = t122 + t287;
t133 = -t259 + (pkin(3) + t506) * t290 + t511 * t289;
t124 = t287 + t132;
t115 = t288 + t123;
t13 = t345 * t289 + t290 * t526;
t14 = t346 * t289 - t290 * t527;
t15 = -t289 * t526 + t345 * t290;
t16 = t289 * t527 + t346 * t290;
t127 = t153 * t470;
t62 = -t149 * t290 - t151 * t471 + t127;
t128 = t154 * t470;
t63 = t150 * t290 + t152 * t471 - t128;
t96 = t289 * t370 - t175;
t94 = t96 * t298;
t30 = -t212 * t63 - t213 * t62 + t94;
t129 = t151 * t466;
t64 = -t149 * t289 - t153 * t465 + t129;
t65 = t150 * t289 - t290 * t375;
t31 = -t212 * t65 - t213 * t64 - t501;
t40 = -t299 * t399 + t301 * t397;
t41 = t299 * t398 - t301 * t396;
t46 = t342 * t289 + t290 * t528;
t47 = -t289 * t528 + t342 * t290;
t92 = t151 * t301 + t153 * t299;
t332 = (-t13 * t213 + t134 * t65 + t135 * t64 - t14 * t212 - t297 * t97 + t298 * t46) * t516 + (t349 * t289 + t564 * t290) * t519 + (-t564 * t289 + t349 * t290) * t518 + (t134 * t63 + t135 * t62 - t15 * t213 - t16 * t212 + t297 * t96 + t298 * t47) * t515 + (-t299 * t525 + t301 * t327) * t513 + t30 * t406 + t31 * t405 + ((-t298 * t65 - t13) * t290 + (t298 * t64 - t14) * t289) * t520 + (-t289 * t65 - t290 * t64) * t524 + (-t289 * t63 - t290 * t62) * t523 + ((-t298 * t63 - t15) * t290 + (t298 * t62 - t16) * t289) * t517 + (-t289 * t93 - t290 * t92) * t514 + ((-t298 * t93 - t40) * t290 + (t298 * t92 - t41) * t289) * t512;
t331 = -rSges(5,1) * t417 - t173 - t442;
t125 = t288 + t133;
t326 = t331 - t426;
t116 = t289 * t368 - t191;
t113 = t116 * t298;
t36 = qJD(4) * t382 + t113;
t37 = qJD(4) * t381 - t487;
t50 = qJD(4) * t374 + t103 * t317 + t105 * t315;
t51 = qJD(4) * t372 + t102 * t317 + t104 * t315;
t56 = t341 * t289 + t290 * t532;
t57 = -t289 * t532 + t341 * t290;
t324 = t93 * t524 - t373 * t522 + (t113 + ((t72 + t140 - t141 + (t161 - t480) * t289) * t289 + (-t139 - t73 + (t161 - t372) * t290 + (t479 + t481) * t289) * t290) * qJD(4)) * t403 + (t94 - (t289 * t395 + t127 + t65) * t213 + (t128 - t129 + t64 + (t149 - t483) * t289) * t212 + (t212 * t394 - t548) * t290) * t519 - t134 * t97 / 0.2e1 - t187 * t117 / 0.2e1 + (t92 + t96) * t523 + (t116 + t109) * t521 + (t31 + t501 - (-t129 + t63) * t213 + (-t127 + t62) * t212 + (-t212 * t361 - t213 * t394) * t290 + (-t212 * t395 + t548) * t289) * t518 + (t40 + t47) * t517 + (t50 + t57) * t402 + (t487 + ((t141 - t71 + (t162 - t479) * t290) * t290 + (-t139 + t70 + (t162 + t481) * t289) * t289) * qJD(4) + t37) * t401 + (t41 + t46 + t30) * t520 + (t51 + t56 + t36) * t404 + (qJD(4) * t368 + t232 * t317 + t233 * t315 - t299 * t391 + t301 * t392) * t298 + (t222 * t301 + t299 * t545 + Icges(4,3) + t369) * t297;
t323 = Icges(3,3) * t310 + t324;
t260 = pkin(4) * t464;
t138 = t289 * t155;
t95 = t371 * qJD(4);
t22 = t289 * t531 + t344 * t290;
t21 = -t289 * t535 + t343 * t290;
t20 = t344 * t289 - t290 * t531;
t19 = t343 * t289 + t290 * t535;
t10 = -t134 * t155 - t135 * t156 - t145 * t187 + t146 * t188 + t212 * t91 - t213 * t90 + (t289 * t99 - t290 * t98) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t323 + (t550 * (t229 + t308) + t551 * (t306 + t544) + (-t184 + t205 + t296) * t183) * m(3) + ((t273 ^ 2 + t275 ^ 2) * qJDD(1) - g(2) * t275 - g(3) * t273) * m(2) + (t59 * (-t354 + t444) + t322 + t557 * (t308 + t115) + t555 * (t114 + t306) + (t59 + t541) * t58) * m(6) + (t78 * (t326 - t420) + t325 + t556 * (t308 + t125) + t554 * (t306 + t124) + (t78 + t540) * t77) * m(5) + (t552 * (t186 + t308) + t553 * (t306 + t185) + (-t144 + t296 + t365) * t143) * m(4); t323 + ((t336 - t426 + t444) * t59 + t557 * t115 + t555 * t114 + t565) * m(6) + ((t326 - t347 + t426) * t78 + t556 * t125 + t554 * t124 + t566) * m(5) + (t552 * t186 + t553 * t185 + (t365 - t271 - t196) * t143) * m(4) + (t183 * t205 - t184 * t204 + (-t183 * t312 + t550) * t229 + (t184 * t312 + t551) * t544) * m(3); t324 + ((-t366 + t444) * t59 + t557 * t123 + t555 * t122 + t565) * m(6) + ((t331 - t347) * t78 + t556 * t133 + t554 * t132 + t566) * m(5) + (t143 * t172 - t144 * t560 + (-t143 * t298 + t552) * t215 + (t144 * t298 + t553) * t546) * m(4); ((-t431 * t474 - t473) * t290 + (-t561 + (-t339 * t289 + (t191 + t340) * t290) * qJD(4)) * t289) * t401 + ((-t315 * t436 + t317 * t435) * t298 + ((t289 * t447 - t290 * t448) * t317 + (-t289 * t449 + t290 * t450) * t315) * qJD(4)) * t513 + ((-t298 * t71 - t21) * t290 + (t298 * t70 - t22) * t289) * t402 + ((-t298 * t73 - t19) * t290 + (t298 * t72 - t20) * t289) * t404 + ((t191 * t432 - t473) * t289 + (t561 + (-t340 * t290 + (-t474 + t339) * t289) * qJD(4)) * t290) * t403 + t332 + ((t298 * t373 - t50) * t290 + (t109 * t298 - t51) * t289) * t512 + t36 * t406 + t37 * t405 + t381 * t522 + t382 * t521 + (-t109 * t290 + t289 * t373) * t514 + (t116 * t297 + t187 * t71 + t188 * t70 + t298 * t57 + (-t21 * t290 - t22 * t289) * qJD(4)) * t515 + (-t117 * t297 + t187 * t73 + t188 * t72 + t298 * t56 + (-t19 * t290 - t20 * t289) * qJD(4)) * t516 + (-g(1) * (t228 + t307) - g(3) * (t260 + t181) - g(2) * t400 * t289 + t10 * t138 + t52 * t418 + t17 * t260 + (t10 * t455 + t52 * t507 + t17 * t226 + t58 * t203 + (t52 * t145 + t400 * t59) * t298) * t290 + (t10 * t145 + t52 * t99 + t18 * t400 + t59 * (-pkin(4) * t428 - t203) + (t52 * t146 + t400 * t58) * t298) * t289 - (-t59 * t409 + ((-t289 * t58 - t290 * t59) * t298 + t52 * (-t289 ^ 2 - t290 ^ 2) * qJD(4)) * t315) * pkin(4) + t538) * m(6) + ((-t167 * t187 - t168 * t188 + (-t107 * t290 + t108 * t289) * qJD(4)) * t371 + t95 * ((t167 * t298 - t107) * t290 + (t108 - t147) * t289) + t380 * t239 + ((-t298 * t78 + t44) * t290 + (-t298 * t77 - t45) * t289) * t272 - t559 * t298 - (t95 * (-t197 * t289 - t198 * t290) + t380 * t274) * qJD(4) - g(1) * t274 + g(2) * t197 - g(3) * t198) * m(5); t332 + (t10 * (t156 * t290 + t138) + t52 * (-t290 * t90 + t418) + (-t289 * t59 + t290 * t58) * t203 + ((t17 - t502) * t290 + (-t298 * t58 - t18) * t289) * t226 - g(1) * t228 + g(2) * t180 - g(3) * t181 + t538) * m(6);];
tau = t1;
