% Calculate vector of inverse dynamics joint torques for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR7_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:39
% DurationCPUTime: 7.94s
% Computational Cost: add. (14670->638), mult. (12473->789), div. (0->0), fcn. (9599->8), ass. (0->363)
t320 = sin(qJ(4));
t432 = qJD(4) * t320;
t559 = rSges(5,2) * t432;
t317 = qJD(1) + qJD(2);
t322 = cos(qJ(4));
t486 = Icges(5,4) * t322;
t269 = -Icges(5,2) * t320 + t486;
t383 = Icges(5,1) * t320 + t486;
t439 = t269 + t383;
t487 = Icges(5,4) * t320;
t271 = Icges(5,1) * t322 - t487;
t381 = Icges(5,2) * t322 + t487;
t440 = -t381 + t271;
t558 = (t320 * t439 - t322 * t440) * t317;
t319 = qJ(1) + qJ(2);
t312 = cos(t319);
t318 = qJ(4) + qJ(5);
t309 = sin(t318);
t492 = rSges(6,2) * t309;
t255 = t312 * t492;
t311 = cos(t318);
t493 = rSges(6,1) * t311;
t191 = -t312 * t493 + t255;
t321 = sin(qJ(1));
t489 = pkin(1) * qJD(1);
t425 = t321 * t489;
t310 = sin(t319);
t232 = rSges(3,1) * t310 + rSges(3,2) * t312;
t468 = t232 * t317;
t180 = -t425 - t468;
t387 = rSges(6,1) * t309 + rSges(6,2) * t311;
t463 = t310 * t311;
t251 = Icges(6,4) * t463;
t465 = t309 * t310;
t160 = Icges(6,1) * t465 + Icges(6,5) * t312 + t251;
t484 = Icges(6,4) * t311;
t382 = Icges(6,1) * t309 + t484;
t161 = -Icges(6,5) * t310 + t312 * t382;
t316 = qJD(4) + qJD(5);
t221 = t310 * t316;
t222 = t312 * t316;
t226 = -Icges(6,2) * t309 + t484;
t485 = Icges(6,4) * t309;
t228 = Icges(6,1) * t311 - t485;
t380 = Icges(6,2) * t311 + t485;
t333 = t221 * (t226 * t312 + t161) - t222 * (-Icges(6,2) * t465 + t160 + t251) - t317 * (-t380 + t228);
t358 = t380 * t310;
t158 = Icges(6,6) * t312 + t358;
t357 = t380 * t312;
t159 = -Icges(6,6) * t310 + t357;
t334 = t221 * (-t228 * t312 + t159) - t222 * (-t228 * t310 + t158) - t317 * (t226 + t382);
t557 = t334 * t309 - t333 * t311;
t324 = -pkin(8) - pkin(7);
t457 = t317 * t324;
t250 = t312 * t457;
t431 = qJD(4) * t322;
t415 = t312 * t431;
t260 = pkin(4) * t415;
t461 = t310 * t320;
t283 = pkin(4) * t461;
t304 = t312 * pkin(7);
t129 = -t250 - t260 + (t283 - t304) * t317;
t433 = qJD(4) * t317;
t208 = qJDD(4) * t310 + t312 * t433;
t459 = t312 * t317;
t138 = qJD(5) * t459 + qJDD(5) * t310 + t208;
t305 = t312 * pkin(2);
t234 = qJ(3) * t310 + t305;
t289 = qJD(3) * t312;
t168 = t234 * t317 - t289;
t195 = t387 * t316;
t233 = -t492 + t493;
t315 = qJDD(1) + qJDD(2);
t314 = t317 ^ 2;
t323 = cos(qJ(1));
t326 = qJD(1) ^ 2;
t353 = (-qJDD(1) * t321 - t323 * t326) * pkin(1);
t436 = qJD(3) * t317;
t338 = qJDD(3) * t310 + t312 * t436 + t353;
t331 = -t304 * t314 + t338;
t298 = t310 * rSges(6,3);
t164 = t312 * t387 - t298;
t458 = t312 * t320;
t182 = pkin(4) * t458 + (pkin(7) + t324) * t310;
t291 = t312 * qJ(3);
t230 = pkin(2) * t310 - t291;
t495 = pkin(7) * t310;
t408 = -t230 - t495;
t367 = t164 + t182 + t408;
t456 = t320 * qJD(4) ^ 2;
t301 = t312 * rSges(6,3);
t366 = t191 * t316;
t97 = (t310 * t387 + t301) * t317 + t366;
t17 = t138 * t233 - t195 * t221 + (t208 * t322 - t310 * t456) * pkin(4) + (-t129 - t168 - t97) * t317 + t367 * t315 + t331;
t556 = t17 - g(1);
t417 = t310 * t431;
t258 = pkin(4) * t417;
t423 = t317 * t458;
t419 = pkin(4) * t423 + t310 * t457 + t258;
t462 = t310 * t317;
t428 = pkin(7) * t462;
t128 = t419 + t428;
t286 = qJDD(4) * t312;
t139 = qJDD(5) * t312 - t316 * t462 + t286;
t209 = -t310 * t433 + t286;
t313 = t323 * pkin(1);
t497 = pkin(1) * t321;
t394 = qJDD(1) * t313 - t326 * t497;
t265 = qJ(3) * t459;
t288 = qJD(3) * t310;
t441 = t265 + t288;
t351 = t317 * (-pkin(2) * t462 + t441) + t315 * t234 + t310 * t436 + t394;
t335 = t304 * t315 - t314 * t495 + t351;
t430 = qJDD(3) * t312;
t163 = rSges(6,1) * t465 + rSges(6,2) * t463 + t301;
t395 = -t312 * t324 + t283;
t183 = -t304 + t395;
t451 = t163 + t183;
t254 = rSges(6,1) * t463;
t421 = t316 * t254 + t387 * t459;
t426 = t316 * t492;
t98 = (-rSges(6,3) * t317 - t426) * t310 + t421;
t18 = -t430 - t139 * t233 + t195 * t222 + (t128 + t98) * t317 + t451 * t315 + (-t209 * t322 + t312 * t456) * pkin(4) + t335;
t555 = t18 - g(2);
t302 = t312 * rSges(5,3);
t364 = -rSges(5,1) * t415 + t312 * t559;
t548 = rSges(5,2) * t322;
t388 = rSges(5,1) * t320 + t548;
t117 = (t310 * t388 + t302) * t317 + t364;
t243 = t388 * qJD(4);
t278 = rSges(5,1) * t322 - rSges(5,2) * t320;
t179 = -rSges(5,3) * t310 + t312 * t388;
t392 = t179 + t408;
t435 = qJD(4) * t310;
t40 = -t243 * t435 + t208 * t278 + (-t117 - t168) * t317 + t392 * t315 + t331;
t554 = t40 - g(1);
t420 = t459 * t548 + (t417 + t423) * rSges(5,1);
t118 = (-rSges(5,3) * t317 - t559) * t310 + t420;
t460 = t310 * t322;
t178 = rSges(5,1) * t461 + rSges(5,2) * t460 + t302;
t41 = t118 * t317 + t178 * t315 - t209 * t278 + (qJD(4) * t243 - qJDD(3)) * t312 + t335;
t553 = t41 - g(2);
t275 = rSges(4,2) * t459;
t490 = rSges(4,3) * t312;
t231 = rSges(4,2) * t310 + t490;
t444 = -t230 + t231;
t64 = (-rSges(4,3) * t462 - t168 + t275) * t317 + t444 * t315 + t338;
t552 = t64 - g(1);
t235 = -rSges(4,2) * t312 + rSges(4,3) * t310;
t438 = rSges(4,2) * t462 + rSges(4,3) * t459;
t65 = t235 * t315 + t317 * t438 + t351 - t430;
t551 = t65 - g(2);
t197 = rSges(3,1) * t459 - rSges(3,2) * t462;
t550 = -t197 * t317 - t232 * t315 - g(1) + t353;
t236 = rSges(3,1) * t312 - rSges(3,2) * t310;
t549 = t236 * t315 - t317 * t468 - g(2) + t394;
t536 = t221 * t233 + t258;
t389 = t288 + t536;
t60 = t317 * t367 + t389 - t425;
t547 = t317 * t60;
t190 = -rSges(6,2) * t465 + t254;
t424 = t323 * t489;
t390 = -t289 + t424;
t407 = t234 + t304;
t529 = -t222 * t233 + t317 * (t407 + t451);
t61 = -t260 + t390 + t529;
t546 = t61 * (t190 * t317 + t222 * t387);
t378 = Icges(6,5) * t309 + Icges(6,6) * t311;
t157 = -Icges(6,3) * t310 + t312 * t378;
t67 = -t157 * t312 - t159 * t463 - t161 * t465;
t371 = t226 * t311 + t228 * t309;
t224 = Icges(6,5) * t311 - Icges(6,6) * t309;
t472 = t224 * t312;
t99 = t310 * t371 + t472;
t545 = t221 * t67 + t317 * t99;
t544 = t190 * t221 - t191 * t222 + t312 * t97;
t171 = t234 + t235;
t543 = t171 * t317;
t376 = t159 * t311 + t161 * t309;
t542 = t312 * t376;
t541 = t317 * t495;
t175 = -Icges(5,6) * t310 + t312 * t381;
t177 = -Icges(5,5) * t310 + t312 * t383;
t373 = t175 * t322 + t177 * t320;
t538 = t373 * t312;
t354 = t378 * t317;
t537 = t178 + t407;
t212 = t317 * t230;
t391 = t288 - t425;
t365 = -t212 + t391;
t535 = -t365 + t428 - t425;
t534 = -t231 * t317 + t438;
t469 = t228 * t316;
t533 = -Icges(6,5) * t317 + t469;
t470 = t226 * t316;
t532 = -Icges(6,6) * t317 + t470;
t148 = t317 * t164;
t167 = t317 * t182;
t531 = -t310 * t426 - t148 - t167 + t419 + t421 + t441;
t162 = t317 * t179;
t214 = t278 * t435;
t418 = t310 * t432;
t530 = -rSges(5,2) * t418 - t162 - t214 + t420 + t441;
t434 = qJD(4) * t312;
t528 = -t278 * t434 + t317 * t537;
t527 = -Icges(6,3) * t317 + t224 * t316;
t106 = t175 * t320 - t177 * t322;
t359 = t381 * t317;
t524 = -Icges(5,6) * t317 + qJD(4) * t269;
t111 = t310 * t359 - t312 * t524;
t361 = t383 * t317;
t523 = -Icges(5,5) * t317 + qJD(4) * t271;
t113 = t310 * t361 - t312 * t523;
t379 = Icges(5,5) * t320 + Icges(5,6) * t322;
t355 = t379 * t312;
t173 = -Icges(5,3) * t310 + t355;
t526 = qJD(4) * t106 + t111 * t322 + t113 * t320 + t173 * t317;
t267 = Icges(5,5) * t322 - Icges(5,6) * t320;
t525 = -Icges(5,3) * t317 + qJD(4) * t267;
t239 = t381 * qJD(4);
t240 = t383 * qJD(4);
t369 = t269 * t320 - t271 * t322;
t522 = qJD(4) * t369 + t239 * t322 + t240 * t320 + t267 * t317;
t112 = t310 * t524 + t312 * t359;
t114 = t310 * t523 + t312 * t361;
t356 = t379 * t310;
t172 = Icges(5,3) * t312 + t356;
t174 = Icges(5,6) * t312 + t310 * t381;
t280 = Icges(5,4) * t460;
t176 = Icges(5,1) * t461 + Icges(5,5) * t312 + t280;
t374 = t174 * t320 - t176 * t322;
t521 = qJD(4) * t374 - t112 * t322 - t114 * t320 + t172 * t317;
t447 = t269 * t312 + t177;
t449 = -t271 * t312 + t175;
t519 = t320 * t449 - t322 * t447;
t448 = -Icges(5,2) * t461 + t176 + t280;
t450 = -t271 * t310 + t174;
t518 = t320 * t450 - t322 * t448;
t396 = t316 * t382 + t470;
t397 = -t316 * t380 + t469;
t517 = t224 * t317 + t309 * t396 - t311 * t397;
t402 = t161 * t316 + t312 * t532 - t317 * t358;
t360 = t382 * t317;
t404 = t159 * t316 + t310 * t360 - t312 * t533;
t516 = t157 * t317 + t309 * t404 - t311 * t402;
t156 = Icges(6,3) * t312 + t310 * t378;
t403 = t160 * t316 + t310 * t532 + t317 * t357;
t405 = t158 * t316 - t310 * t533 - t312 * t360;
t515 = t156 * t317 + t309 * t405 - t311 * t403;
t514 = t310 ^ 2;
t513 = -pkin(2) - pkin(7);
t512 = t138 / 0.2e1;
t511 = t139 / 0.2e1;
t510 = t208 / 0.2e1;
t509 = t209 / 0.2e1;
t508 = -t221 / 0.2e1;
t507 = t221 / 0.2e1;
t506 = -t222 / 0.2e1;
t505 = t222 / 0.2e1;
t504 = t310 / 0.2e1;
t503 = -t312 / 0.2e1;
t502 = t312 / 0.2e1;
t501 = t315 / 0.2e1;
t500 = -t317 / 0.2e1;
t499 = t317 / 0.2e1;
t498 = -rSges(6,3) - pkin(2);
t496 = pkin(4) * t322;
t184 = t310 * t224;
t100 = t312 * t371 - t184;
t477 = t100 * t317;
t198 = t310 * t267;
t370 = t269 * t322 + t271 * t320;
t123 = t312 * t370 - t198;
t476 = t123 * t317;
t467 = t379 * t317;
t466 = t267 * t312;
t464 = t310 * t157;
t442 = t260 + t289;
t437 = t288 - t212;
t429 = -rSges(5,3) + t513;
t284 = pkin(4) * t460;
t66 = t156 * t312 + t158 * t463 + t160 * t465;
t75 = t172 * t312 + t174 * t460 + t176 * t461;
t76 = -t173 * t312 - t175 * t460 - t177 * t461;
t414 = -t462 / 0.2e1;
t413 = t459 / 0.2e1;
t412 = -t435 / 0.2e1;
t411 = t435 / 0.2e1;
t410 = -t434 / 0.2e1;
t409 = t434 / 0.2e1;
t406 = t317 * t61 + t17;
t401 = -t98 + t148;
t399 = -t191 * t317 - t221 * t387;
t279 = rSges(2,1) * t323 - rSges(2,2) * t321;
t277 = rSges(2,1) * t321 + rSges(2,2) * t323;
t386 = t310 * t76 + t312 * t75;
t151 = t310 * t172;
t375 = t174 * t322 + t176 * t320;
t77 = -t312 * t375 + t151;
t78 = -t310 * t173 + t538;
t385 = t310 * t78 + t312 * t77;
t81 = t317 * t392 + t214 + t391;
t82 = t390 + t528;
t384 = t310 * t81 - t312 * t82;
t377 = t158 * t311 + t160 * t309;
t90 = t159 * t309 - t161 * t311;
t372 = -t178 * t310 - t179 * t312;
t363 = t66 + t464;
t362 = t289 - t364;
t170 = t490 + t291 + (rSges(4,2) - pkin(2)) * t310;
t352 = pkin(4) * t320 + t387;
t348 = t184 * t222 - t221 * t472 - t354;
t346 = t310 * t354 - t312 * t527 - t317 * t376;
t345 = t310 * t527 + t312 * t354 + t317 * t377;
t344 = -t312 * t525 + (t356 - t373) * t317;
t343 = t310 * t525 + (t355 + t375) * t317;
t140 = t310 * t156;
t68 = -t312 * t377 + t140;
t342 = -t316 * t378 + t317 * t371;
t341 = -qJD(4) * t379 + t317 * t370;
t340 = t250 - t366 + t442;
t108 = t395 + t163 + t234;
t13 = t310 * t345 + t312 * t515;
t14 = t310 * t346 - t312 * t516;
t15 = -t310 * t515 + t312 * t345;
t16 = t310 * t516 + t312 * t346;
t29 = t222 * t66 + t545;
t69 = -t464 + t542;
t30 = t221 * t69 + t222 * t68 - t477;
t44 = -t309 * t403 - t311 * t405;
t45 = t309 * t402 + t311 * t404;
t46 = t310 * t342 + t312 * t517;
t47 = -t310 * t517 + t312 * t342;
t89 = -t158 * t309 + t160 * t311;
t336 = (-t100 * t315 + t13 * t222 + t138 * t69 + t139 * t68 + t14 * t221 + t317 * t46) * t504 + (t348 * t310 - t312 * t557) * t508 + (t310 * t557 + t348 * t312) * t506 + (t138 * t67 + t139 * t66 + t15 * t222 + t16 * t221 + t315 * t99 + t317 * t47) * t502 + (t309 * t333 + t311 * t334) * t500 + t29 * t414 + t30 * t413 + ((t317 * t69 + t13) * t312 + (-t317 * t68 + t14) * t310) * t507 + (t310 * t69 + t312 * t68) * t512 + (t310 * t67 + t312 * t66) * t511 + ((t317 * t67 + t15) * t312 + (-t317 * t66 + t16) * t310) * t505 + (t310 * t90 + t312 * t89) * t501 + ((t317 * t90 + t44) * t312 + (-t317 * t89 + t45) * t310) * t499;
t130 = t310 * t513 + t179 + t291;
t124 = t317 * t444 + t391;
t125 = t390 + t543;
t330 = (-t124 * t305 + (t124 * (-rSges(4,3) - qJ(3)) - t125 * pkin(2)) * t310) * t317;
t107 = t291 - t298 + (-pkin(2) + t324) * t310 + t352 * t312;
t329 = (t81 * t429 * t312 + (t81 * (-qJ(3) - t388) + t82 * t429) * t310) * t317;
t328 = (t60 * t498 * t312 + (t60 * (-qJ(3) - t352) + t61 * t498) * t310) * t317;
t122 = t310 * t370 + t466;
t119 = t122 * t317;
t36 = qJD(4) * t386 + t119;
t37 = qJD(4) * t385 - t476;
t50 = -qJD(4) * t375 - t112 * t320 + t114 * t322;
t51 = qJD(4) * t373 - t111 * t320 + t113 * t322;
t54 = t310 * t341 + t312 * t522;
t55 = -t310 * t522 + t312 * t341;
t327 = t106 * t510 + t90 * t512 + ((t363 + t69 - t542) * t222 + t545) * t508 - t208 * t123 / 0.2e1 - t138 * t100 / 0.2e1 + (t119 + ((-t77 + t151 + t76) * t310 + (t78 - t538 + (t173 - t375) * t310 + t75) * t312) * qJD(4)) * t412 + (t89 + t99) * t511 + (t122 - t374) * t509 + (t477 + (t376 * t310 - t140 + t67) * t222 + (t363 - t66) * t221 + ((t157 + t377) * t222 - t376 * t221) * t312 + t30) * t506 + (t44 + t47) * t505 + (t476 + (t514 * t173 + (-t151 + t76 + (t173 + t375) * t312) * t312) * qJD(4) + t37) * t410 + (t50 + t55) * t409 + (-qJD(4) * t370 + t239 * t320 - t240 * t322 - t309 * t397 - t311 * t396) * t317 + (t45 + t46 + t29) * t507 + (t51 + t54 + t36) * t411 + (-t226 * t309 + t228 * t311 + Icges(4,1) + Icges(3,3) - t369) * t315;
t206 = t278 * t312;
t205 = t278 * t310;
t181 = t236 * t317 + t424;
t143 = t312 * t164;
t102 = t372 * qJD(4);
t57 = -t163 * t221 - t164 * t222 + (-t182 * t312 - t183 * t310) * qJD(4);
t22 = t310 * t526 + t312 * t344;
t21 = -t310 * t521 + t312 * t343;
t20 = t310 * t344 - t312 * t526;
t19 = t310 * t343 + t312 * t521;
t12 = -t138 * t163 - t139 * t164 - t182 * t209 - t183 * t208 - t221 * t98 + t222 * t97 + (-t128 * t310 + t129 * t312) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t327 + (t549 * (t236 + t313) + t550 * (-t232 - t497) + (-t197 - t424 + t181) * t180) * m(3) + ((t277 ^ 2 + t279 ^ 2) * qJDD(1) + g(1) * t277 - g(2) * t279) * m(2) + (t60 * (t340 - t424) + t328 + (t60 + t531 + t535 - t536) * t61 + t555 * (t108 + t313) + t556 * (t107 - t497)) * m(6) + (t81 * (t362 - t424) + t329 + (t81 + t530 + t535) * t82 + t553 * (t313 + t537) + t554 * (t130 - t497)) * m(5) + (t124 * (t275 - t390) + t330 + t551 * (t171 + t313) + t552 * (t170 - t497) + (t265 + t391 + t124 - t365 + t534) * t125) * m(4); t327 + (t328 + (t212 - t389 + t531 + t541) * t61 + (t340 - t442 + t529) * t60 + t555 * t108 + t556 * t107) * m(6) + (t329 + (-t437 + t530 + t541) * t82 + (-t289 + t362 + t528) * t81 + t553 * t537 + t554 * t130) * m(5) + (t330 + t551 * t171 + t552 * t170 + (-t437 + t441 + t534) * t125 + (t275 + t543) * t124) * m(4) + (-t180 * t197 - t181 * t468 + (t180 * t317 + t549) * t236 + (t181 * t317 - t550) * t232) * m(3); (-m(4) - m(5) - m(6)) * (g(1) * t310 - g(2) * t312) + 0.2e1 * (t17 * t504 + t18 * t503) * m(6) + 0.2e1 * (t40 * t504 + t41 * t503) * m(5) + 0.2e1 * (t503 * t65 + t504 * t64) * m(4); t336 + ((t198 * t434 - t467) * t312 + (-t558 + (t519 * t310 + (-t466 - t518) * t312) * qJD(4)) * t310) * t410 + (-t123 * t315 + t208 * t78 + t209 * t77 + t317 * t54 + (t19 * t312 + t20 * t310) * qJD(4)) * t504 + ((t106 * t317 + t50) * t312 + (t317 * t374 + t51) * t310) * t499 + ((-t320 * t440 - t322 * t439) * t317 + ((t310 * t449 - t312 * t450) * t322 + (t310 * t447 - t312 * t448) * t320) * qJD(4)) * t500 + (t122 * t315 + t208 * t76 + t209 * t75 + t317 * t55 + (t21 * t312 + t22 * t310) * qJD(4)) * t502 + ((-t435 * t466 - t467) * t310 + (t558 + (t518 * t312 + (t198 - t519) * t310) * qJD(4)) * t312) * t412 + t385 * t510 + t386 * t509 + (t106 * t310 - t312 * t374) * t501 + ((t317 * t78 + t19) * t312 + (-t317 * t77 + t20) * t310) * t411 + ((t317 * t76 + t21) * t312 + (-t317 * t75 + t22) * t310) * t409 + t36 * t414 + t37 * t413 + (t17 * t284 - t12 * t143 + (t18 * (-t233 - t496) + t61 * t195 - t12 * t182 + t233 * t547) * t312 + (t60 * (-pkin(4) * t432 - t195) - t12 * t451 + t406 * t233) * t310 - t60 * (-pkin(4) * t418 + t399) - t546 - g(1) * (t190 + t284) - g(2) * (t255 + (-t493 - t496) * t312) + g(3) * t352 + ((-t317 * t451 + t129) * t312 + (-t128 + t401 + t167) * t310 - (-t312 ^ 2 - t514) * pkin(4) * t431 + t544) * t57) * m(6) + (-(t205 * t82 + t206 * t81) * t317 - (t102 * (-t205 * t310 - t206 * t312) - t384 * t388) * qJD(4) + (-t178 * t208 - t179 * t209 + (t117 * t312 - t118 * t310) * qJD(4)) * t372 + t102 * ((-t178 * t317 + t117) * t312 + (-t118 + t162) * t310) - t384 * t243 + ((t317 * t81 - t41) * t312 + (t317 * t82 + t40) * t310) * t278 - g(1) * t205 + g(2) * t206 + g(3) * t388) * m(5); t336 + (-t546 - t60 * t399 + t12 * (-t163 * t310 - t143) - (t310 * t60 - t312 * t61) * t195 + ((-t18 + t547) * t312 + t406 * t310) * t233 - g(1) * t190 - g(2) * t191 + g(3) * t387 + (-t163 * t459 + t310 * t401 + t544) * t57) * m(6);];
tau = t1;
