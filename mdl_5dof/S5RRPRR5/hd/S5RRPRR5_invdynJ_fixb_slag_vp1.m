% Calculate vector of inverse dynamics joint torques for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:52
% DurationCPUTime: 10.63s
% Computational Cost: add. (19881->644), mult. (13225->791), div. (0->0), fcn. (10283->10), ass. (0->384)
t321 = pkin(9) + qJ(4);
t313 = sin(t321);
t324 = qJ(1) + qJ(2);
t316 = sin(t324);
t457 = qJD(4) * t316;
t436 = t313 * t457;
t269 = pkin(4) * t436;
t540 = -pkin(7) - qJ(3);
t320 = -pkin(8) + t540;
t323 = qJD(1) + qJD(2);
t484 = t320 * t323;
t275 = t316 * t484;
t429 = t323 * t540;
t279 = t316 * t429;
t314 = cos(t321);
t305 = pkin(4) * t314;
t326 = cos(pkin(9));
t437 = pkin(3) * t326 + pkin(2);
t260 = t305 + t437;
t408 = -t260 + t437;
t317 = cos(t324);
t486 = t317 * t323;
t102 = t408 * t486 + t269 + t275 - t279;
t299 = qJDD(4) * t317;
t322 = qJD(4) + qJD(5);
t488 = t316 * t323;
t160 = qJDD(5) * t317 - t322 * t488 + t299;
t315 = qJ(5) + t321;
t310 = cos(t315);
t297 = t310 * rSges(6,1);
t309 = sin(t315);
t232 = -rSges(6,2) * t309 + t297;
t199 = t232 * t322;
t455 = qJD(4) * t323;
t215 = -t316 * t455 + t299;
t231 = rSges(6,1) * t309 + rSges(6,2) * t310;
t253 = t317 * t322;
t319 = qJDD(1) + qJDD(2);
t491 = t314 * qJD(4) ^ 2;
t543 = pkin(4) * t313;
t304 = qJD(3) * t317;
t327 = sin(qJ(1));
t328 = cos(qJ(1));
t330 = qJD(1) ^ 2;
t361 = (-qJDD(1) * t327 - t328 * t330) * pkin(1);
t397 = -t316 * pkin(2) + t317 * qJ(3);
t490 = t316 * qJ(3);
t255 = t317 * pkin(2) + t490;
t220 = t323 * t255;
t458 = t304 - t220;
t340 = t361 + qJDD(3) * t316 + t319 * t397 + (t304 + t458) * t323;
t485 = t317 * t326;
t371 = -pkin(3) * t485 + t490;
t487 = t316 * t326;
t379 = -pkin(3) * t487 + pkin(7) * t317;
t591 = t340 + t323 * (t323 * t371 + t279) + t319 * t379;
t341 = (-t320 + t540) * t317 + t408 * t316;
t497 = t310 * t316;
t451 = rSges(6,1) * t497;
t499 = t309 * t316;
t579 = rSges(6,2) * t499 + t317 * rSges(6,3);
t381 = t451 - t579;
t609 = t341 - t381;
t496 = t310 * t317;
t449 = rSges(6,1) * t496;
t374 = -rSges(6,3) * t316 - t449;
t192 = rSges(6,1) * t499 + rSges(6,2) * t497;
t498 = t309 * t317;
t274 = rSges(6,2) * t498;
t441 = t192 * t322 + t323 * t274;
t99 = t323 * t374 + t441;
t13 = -pkin(4) * t317 * t491 - t160 * t231 - t253 * t199 - t215 * t543 + (t102 + t99) * t323 + t609 * t319 + t591;
t611 = g(3) - t13;
t212 = t317 * t260;
t295 = t316 * t540;
t405 = t317 * t437;
t133 = t316 * t320 - t212 - t295 + t405;
t171 = -t274 - t374;
t610 = t133 - t171;
t399 = rSges(3,1) * t316 + rSges(3,2) * t317;
t210 = t399 * t323;
t532 = pkin(1) * qJD(1);
t453 = t327 * t532;
t194 = t210 + t453;
t608 = t323 * t192 - t231 * t488 - t232 * t253;
t252 = t316 * t322;
t503 = t231 * t317;
t595 = t231 * t253;
t440 = -t323 * t451 - t595;
t98 = t323 * t579 + t440;
t607 = t192 * t252 + t253 * t503 - t316 * t99 + t317 * t98;
t461 = t269 + t304;
t402 = t231 * t252 + t461;
t164 = t295 + t371;
t478 = t164 - t255;
t410 = t478 + t610;
t452 = t328 * t532;
t59 = t323 * t410 + t402 - t452;
t606 = (t199 * t316 + t231 * t486 - t252 * t232 - t323 * t503) * t59;
t267 = Icges(6,4) * t499;
t169 = -Icges(6,1) * t497 + Icges(6,5) * t317 + t267;
t268 = Icges(6,4) * t498;
t170 = Icges(6,1) * t496 + Icges(6,5) * t316 - t268;
t296 = Icges(6,4) * t310;
t390 = -Icges(6,2) * t309 + t296;
t577 = Icges(6,1) * t309 + t296;
t587 = t577 + t390;
t337 = t252 * (-Icges(6,2) * t496 + t170 - t268) + t253 * (Icges(6,2) * t497 + t169 + t267) + t323 * t587;
t364 = t390 * t316;
t167 = Icges(6,6) * t317 - t364;
t168 = Icges(6,4) * t496 - Icges(6,2) * t498 + Icges(6,6) * t316;
t522 = Icges(6,4) * t309;
t227 = Icges(6,2) * t310 + t522;
t230 = Icges(6,1) * t310 - t522;
t557 = t252 * (t317 * t577 + t168) + t253 * (-t316 * t577 + t167) + t323 * (t227 - t230);
t605 = t337 * t309 + t310 * t557;
t495 = t313 * t316;
t283 = Icges(5,4) * t495;
t493 = t314 * t316;
t178 = -Icges(5,1) * t493 + Icges(5,5) * t317 + t283;
t494 = t313 * t317;
t284 = Icges(5,4) * t494;
t492 = t314 * t317;
t179 = Icges(5,1) * t492 + Icges(5,5) * t316 - t284;
t348 = t316 * (-Icges(5,2) * t492 + t179 - t284) + t317 * (Icges(5,2) * t493 + t178 + t283);
t302 = Icges(5,4) * t314;
t391 = -Icges(5,2) * t313 + t302;
t365 = t391 * t316;
t176 = Icges(5,6) * t317 - t365;
t177 = Icges(5,4) * t492 - Icges(5,2) * t494 + Icges(5,6) * t316;
t561 = t176 * t317 + t177 * t316;
t604 = -t348 * t313 - t314 * t561;
t456 = qJD(4) * t317;
t435 = t313 * t456;
t411 = pkin(4) * t435;
t368 = -t260 * t488 - t317 * t484 - t411;
t466 = t317 * t429 + t437 * t488;
t101 = t368 + t466;
t214 = qJDD(4) * t316 + t317 * t455;
t159 = qJD(5) * t486 + qJDD(5) * t316 + t214;
t541 = t328 * pkin(1);
t542 = t327 * pkin(1);
t404 = -qJDD(1) * t541 + t330 * t542;
t378 = qJDD(3) * t317 + t404;
t291 = pkin(2) * t488;
t403 = qJ(3) * t486 - t291;
t303 = qJD(3) * t316;
t412 = -0.2e1 * t303 - t403;
t401 = t403 + t466 + t412;
t14 = t159 * t231 + t199 * t252 + (t214 * t313 + t316 * t491) * pkin(4) + t410 * t319 + (-t101 + t401 - t98) * t323 + t378;
t603 = -g(2) + t14;
t530 = t316 * rSges(5,3);
t375 = -rSges(5,1) * t492 - t530;
t438 = rSges(5,1) * t436 + rSges(5,2) * (t313 * t486 + t314 * t457);
t116 = t323 * t375 + t438;
t538 = rSges(5,1) * t314;
t251 = -rSges(5,2) * t313 + t538;
t224 = t251 * qJD(4);
t250 = rSges(5,1) * t313 + rSges(5,2) * t314;
t450 = rSges(5,1) * t493;
t578 = rSges(5,2) * t495 + t317 * rSges(5,3);
t382 = t450 - t578;
t34 = t323 * t116 - t215 * t250 - t224 * t456 - t319 * t382 + t591;
t601 = t34 - g(3);
t433 = t314 * t456;
t439 = -rSges(5,1) * t435 - rSges(5,2) * t433 - t323 * t450;
t115 = t323 * t578 + t439;
t286 = rSges(5,2) * t494;
t183 = -t286 - t375;
t442 = -t183 + t478;
t35 = t224 * t457 + t214 * t250 + t442 * t319 + (-t115 + t401) * t323 + t378;
t600 = t35 - g(2);
t325 = sin(pkin(9));
t536 = rSges(4,2) * t325;
t294 = t317 * t536;
t259 = t323 * t294;
t539 = rSges(4,1) * t326;
t398 = t536 - t539;
t533 = rSges(4,3) * t317;
t345 = t316 * t398 + t533;
t376 = -rSges(4,1) * t485 - rSges(4,3) * t316;
t62 = t319 * t345 + (t323 * t376 + t259) * t323 + t340;
t599 = t62 - g(3);
t258 = t323 * rSges(4,1) * t487;
t191 = -t294 - t376;
t462 = -t255 - t191;
t63 = t462 * t319 + (t258 + (-t316 * t536 - t533) * t323 + t412) * t323 + t378;
t598 = t63 - g(2);
t211 = -rSges(3,1) * t486 + rSges(3,2) * t488;
t597 = t211 * t323 - t319 * t399 - g(3) + t361;
t256 = rSges(3,1) * t317 - t316 * rSges(3,2);
t596 = t210 * t323 - t256 * t319 - g(2) + t404;
t576 = Icges(5,1) * t313 + t302;
t463 = t576 + t391;
t523 = Icges(5,4) * t313;
t246 = Icges(5,2) * t314 + t523;
t249 = Icges(5,1) * t314 - t523;
t464 = t246 - t249;
t594 = (t313 * t463 + t314 * t464) * t323;
t166 = Icges(6,5) * t496 - Icges(6,6) * t498 + Icges(6,3) * t316;
t593 = t317 * t166 - t170 * t497;
t592 = -t323 * t191 - t259;
t153 = t323 * t164;
t583 = t610 * t323;
t590 = t153 - t220 + t402 + t583 - t275 - t441 - t461;
t163 = t323 * t183;
t209 = t250 * t457;
t589 = -t279 - t304 - t438 - t163 + t153 + t209;
t387 = t177 * t313 - t179 * t314;
t586 = t317 * t387;
t581 = t609 * t323;
t580 = t258 + t291;
t226 = Icges(6,5) * t310 - Icges(6,6) * t309;
t362 = t226 * t316;
t165 = Icges(6,3) * t317 - t362;
t481 = -t316 * t165 - t169 * t496;
t575 = t481 - t593;
t505 = t227 * t322;
t573 = -Icges(6,6) * t323 + t505;
t162 = t323 * t382;
t218 = t323 * t397;
t469 = -t218 - t303;
t443 = -t323 * t379 + t469;
t356 = t250 * t456 + t162 + t443;
t70 = t356 + t453;
t406 = -t304 + t452;
t71 = t323 * t442 + t209 - t406;
t572 = t316 * t71 + t317 * t70;
t225 = Icges(6,5) * t309 + Icges(6,6) * t310;
t570 = -Icges(6,3) * t323 + t225 * t322;
t569 = -Icges(6,5) * t323 + t322 * t577;
t104 = t176 * t314 + t178 * t313;
t565 = -Icges(5,6) * t323 + qJD(4) * t246;
t112 = t316 * t565 - t391 * t486;
t366 = t249 * t323;
t563 = -Icges(5,5) * t323 + qJD(4) * t576;
t114 = t316 * t563 - t317 * t366;
t245 = Icges(5,5) * t314 - Icges(5,6) * t313;
t174 = Icges(5,3) * t317 - t245 * t316;
t568 = qJD(4) * t104 + t112 * t313 - t114 * t314 - t174 * t323;
t105 = t177 * t314 + t179 * t313;
t111 = -t317 * t565 - t323 * t365;
t113 = -t316 * t366 - t317 * t563;
t175 = Icges(5,5) * t492 - Icges(5,6) * t494 + Icges(5,3) * t316;
t567 = qJD(4) * t105 + t111 * t313 - t113 * t314 - t175 * t323;
t244 = Icges(5,5) * t313 + Icges(5,6) * t314;
t566 = -Icges(5,3) * t323 + qJD(4) * t244;
t222 = t391 * qJD(4);
t223 = t249 * qJD(4);
t385 = t246 * t314 + t313 * t576;
t564 = qJD(4) * t385 + t222 * t313 - t223 * t314 - t244 * t323;
t413 = -t230 * t322 + t505;
t414 = t587 * t322;
t560 = -t225 * t323 + t309 * t414 + t310 * t413;
t419 = t170 * t322 - t317 * t573 - t323 * t364;
t367 = t323 * t230;
t421 = t168 * t322 + t316 * t367 + t317 * t569;
t559 = -t166 * t323 + t309 * t419 + t310 * t421;
t420 = t169 * t322 + t316 * t573 - t390 * t486;
t422 = t167 * t322 - t316 * t569 + t317 * t367;
t558 = -t165 * t323 + t309 * t420 + t310 * t422;
t556 = t159 / 0.2e1;
t555 = t160 / 0.2e1;
t554 = t214 / 0.2e1;
t553 = t215 / 0.2e1;
t552 = -t252 / 0.2e1;
t551 = t252 / 0.2e1;
t550 = -t253 / 0.2e1;
t549 = t253 / 0.2e1;
t548 = t316 / 0.2e1;
t547 = t317 / 0.2e1;
t546 = t319 / 0.2e1;
t545 = -t323 / 0.2e1;
t544 = t323 / 0.2e1;
t531 = pkin(4) * qJD(4);
t185 = t225 * t316;
t386 = t227 * t309 - t310 * t577;
t89 = -t317 * t386 + t185;
t527 = t89 * t323;
t526 = rSges(4,3) + qJ(3);
t200 = t244 * t316;
t384 = t246 * t313 - t314 * t576;
t107 = -t317 * t384 + t200;
t513 = t107 * t323;
t510 = t176 * t313;
t509 = t178 * t314;
t507 = t225 * t317;
t504 = t231 * t316;
t502 = t244 * t317;
t501 = t245 * t323;
t206 = t250 * t316;
t500 = t250 * t317;
t489 = t316 * t166;
t482 = t317 * t165 + t167 * t499;
t480 = t317 * t174 + t176 * t495;
t479 = t316 * t174 + t178 * t492;
t287 = pkin(4) * t495;
t448 = t167 * t498;
t432 = -t488 / 0.2e1;
t431 = t486 / 0.2e1;
t430 = -pkin(2) - t539;
t428 = -t457 / 0.2e1;
t427 = t457 / 0.2e1;
t426 = -t456 / 0.2e1;
t425 = t456 / 0.2e1;
t416 = -t175 - t509;
t407 = -t303 + t453;
t293 = rSges(2,1) * t328 - t327 * rSges(2,2);
t400 = rSges(2,1) * t327 + rSges(2,2) * t328;
t73 = -t178 * t493 + t480;
t74 = t317 * t175 + t177 * t495 - t179 * t493;
t396 = t316 * t74 + t317 * t73;
t75 = -t176 * t494 + t479;
t76 = t175 * t316 - t586;
t395 = t316 * t76 + t317 * t75;
t91 = t168 * t310 + t170 * t309;
t389 = t168 * t309 - t170 * t310;
t388 = -t509 + t510;
t383 = -t437 - t538;
t380 = t220 + t406;
t369 = -t303 - t439 + t466;
t195 = -t256 * t323 - t452;
t358 = t185 * t253 + t226 * t323 - t252 * t507;
t355 = -t317 * t570 + (-t362 + t389) * t323;
t354 = -t226 * t486 + t570 * t316 + (t167 * t309 - t169 * t310) * t323;
t353 = -t501 * t316 - t317 * t566 + t323 * t387;
t352 = t316 * t566 - t501 * t317 + t323 * t388;
t351 = t226 * t322 + t323 * t386;
t350 = t245 * qJD(4) + t323 * t384;
t347 = t317 * t171 + t316 * t381;
t122 = -t449 - t212 + t274 + (-rSges(6,3) + t320) * t316;
t344 = t317 * t183 + t316 * t382;
t121 = -t317 * t320 + (-t260 - t297) * t316 + t579;
t146 = -t316 * t526 + t317 * t430 + t294;
t15 = t354 * t316 - t317 * t558;
t16 = t355 * t316 - t317 * t559;
t17 = t316 * t558 + t354 * t317;
t18 = t316 * t559 + t355 * t317;
t66 = -t169 * t497 + t482;
t137 = t168 * t499;
t67 = t137 + t593;
t88 = t316 * t386 + t507;
t86 = t88 * t323;
t30 = t252 * t67 + t253 * t66 + t86;
t68 = -t448 - t481;
t69 = -t317 * t389 + t489;
t31 = t252 * t69 + t253 * t68 + t527;
t44 = t351 * t316 - t317 * t560;
t45 = t316 * t560 + t351 * t317;
t46 = -t309 * t422 + t310 * t420;
t47 = -t309 * t421 + t310 * t419;
t90 = t167 * t310 + t169 * t309;
t343 = (t15 * t253 + t159 * t69 + t16 * t252 + t160 * t68 + t319 * t89 + t323 * t44) * t548 + (t358 * t316 - t317 * t605) * t552 + (t316 * t605 + t358 * t317) * t550 + (t159 * t67 + t160 * t66 + t17 * t253 + t18 * t252 + t319 * t88 + t323 * t45) * t547 + (-t309 * t557 + t310 * t337) * t545 + t30 * t432 + t31 * t431 + ((t323 * t69 + t15) * t317 + (-t323 * t68 + t16) * t316) * t551 + (t316 * t69 + t317 * t68) * t556 + (t316 * t67 + t317 * t66) * t555 + ((t323 * t67 + t17) * t317 + (-t323 * t66 + t18) * t316) * t549 + (t316 * t91 + t317 * t90) * t546 + ((t323 * t91 + t46) * t317 + (-t323 * t90 + t47) * t316) * t544;
t342 = -t303 - t368 - t440;
t145 = t526 * t317 + (-pkin(2) + t398) * t316;
t132 = t317 * t383 + t286 + t295 - t530;
t339 = t411 + t443 - t581 + t595;
t336 = t316 * t341;
t131 = t316 * t383 - t317 * t540 + t578;
t335 = t168 * t498 - t489 + (-t169 * t316 - t170 * t317) * t310 + t482;
t58 = t339 + t453;
t334 = (-t59 * t579 - t58 * (t374 - t212)) * t323;
t333 = (-t71 * t578 - t70 * (-t405 + t375)) * t323;
t180 = t323 * t345;
t119 = -t180 - t218 + t407;
t120 = t323 * t462 - t406;
t332 = ((t119 * t526 - t120 * t536) * t316 + (-t119 * t430 - t120 * t526) * t317) * t323;
t106 = t316 * t384 + t502;
t100 = t106 * t323;
t36 = qJD(4) * t396 + t100;
t37 = qJD(4) * t395 + t513;
t50 = -qJD(4) * t388 + t112 * t314 + t114 * t313;
t51 = -qJD(4) * t387 + t111 * t314 + t113 * t313;
t54 = t350 * t316 - t317 * t564;
t55 = t316 * t564 + t350 * t317;
t331 = (t100 + ((t480 + t76 + t586) * t317 + (-t75 + (t416 - t510) * t317 + t74 + t479) * t316) * qJD(4)) * t428 + (t86 + (t69 + t335) * t253 + (t137 - t68 - t448 - t575) * t252) * t552 + (t91 + t89) * t556 + (t90 + t88) * t555 + (t107 + t105) * t554 + (t106 + t104) * t553 + (t31 - t527 + (t67 + (t167 * t317 - t168 * t316) * t309 + t575) * t253 + (-t66 + t335) * t252) * t550 + (t46 + t45) * t549 + (t37 - t513 + ((t74 + (-t175 + t510) * t317 - t479) * t317 + (t316 * t416 + t480 - t73) * t316) * qJD(4)) * t426 + (t50 + t55) * t425 + (-qJD(4) * t384 + t222 * t314 + t223 * t313 - t309 * t413 + t310 * t414) * t323 + (t47 + t44 + t30) * t551 + (t51 + t54 + t36) * t427 + (t385 + Icges(4,2) * t326 ^ 2 + (Icges(4,1) * t325 + 0.2e1 * Icges(4,4) * t326) * t325 + Icges(3,3) + t227 * t310 + t577 * t309) * t319;
t103 = t344 * qJD(4);
t57 = -qJD(4) * t336 - t133 * t456 + t253 * t171 + t252 * t381;
t22 = t316 * t567 + t353 * t317;
t21 = t316 * t568 + t352 * t317;
t20 = t353 * t316 - t317 * t567;
t19 = t352 * t316 - t317 * t568;
t10 = t101 * t456 - t102 * t457 - t215 * t133 + t159 * t381 + t160 * t171 - t214 * t341 - t252 * t99 + t253 * t98;
t1 = [Icges(2,3) * qJDD(1) + t331 + (t596 * (-t256 - t541) + t597 * (-t399 - t542) + (t195 - t211 + t452) * t194) * m(3) + ((qJDD(1) * t400 + g(3)) * t400 + (qJDD(1) * t293 + g(2)) * t293) * m(2) + (t59 * (t342 + t453) + t334 + t603 * (t122 - t541) - t611 * (t121 - t542) + (-t59 + t590) * t58) * m(6) + (t71 * (t369 + t453) + t333 + (t452 - t380 - t71 + t589) * t70 + t600 * (t132 - t541) + t601 * (t131 - t542)) * m(5) + (t120 * (t407 + t580) + t332 + t598 * (t146 - t541) + t599 * (t145 - t542) + (t406 - t120 - t380 + t592) * t119) * m(4); t331 + (t334 + (-t339 + t342) * t59 + t590 * t58 + t603 * t122 - t611 * t121) * m(6) + (t333 + (t369 - t356) * t71 + (t458 + t589) * t70 + t600 * t132 + t601 * t131) * m(5) + (t332 + t598 * t146 + t599 * t145 + (t180 - t469 - t303 + t580) * t120 + (t458 - t304 + t592) * t119) * m(4) + (-t194 * t211 + t195 * t210 + (-t194 * t323 - t596) * t256 - (t195 * t323 + t597) * t399) * m(3); (-m(4) - m(5) - m(6)) * (g(2) * t317 + g(3) * t316) + m(4) * (t316 * t62 + t317 * t63) + m(5) * (t316 * t34 + t317 * t35) + m(6) * (t13 * t316 + t14 * t317); (t106 * t319 + t214 * t74 + t215 * t73 + t323 * t55 + (t21 * t317 + t22 * t316) * qJD(4)) * t547 + ((t105 * t323 + t50) * t317 + (-t104 * t323 + t51) * t316) * t544 + ((t200 * t456 + t501) * t317 + (t594 + (-t317 * t502 - t604) * qJD(4)) * t316) * t426 + (t107 * t319 + t214 * t76 + t215 * t75 + t323 * t54 + (t19 * t317 + t20 * t316) * qJD(4)) * t548 + ((t323 * t74 + t21) * t317 + (-t323 * t73 + t22) * t316) * t425 + ((-t313 * t464 + t314 * t463) * t323 + (-t313 * t561 + t314 * t348) * qJD(4)) * t545 + t343 + ((t323 * t76 + t19) * t317 + (-t323 * t75 + t20) * t316) * t427 + ((-t457 * t502 + t501) * t316 + (-t594 + (t316 * t200 + t604) * qJD(4)) * t317) * t428 + t36 * t432 + t37 * t431 + (t104 * t317 + t105 * t316) * t546 + t395 * t554 + t396 * t553 + (t10 * (-t336 + t347) + t14 * (t287 + t504) - g(1) * (t232 + t305) - g(2) * (t287 + t192) + (-t10 * t133 - t611 * (-t231 - t543)) * t317 + ((-t102 + t583) * t316 + (t101 - t581) * t317 - (-t316 ^ 2 - t317 ^ 2) * t313 * t531 + t607) * t57 + t606 + ((t314 * t531 + t199) * t317 - pkin(4) * t433 + t608) * t58) * m(6) + (-(-t206 * t70 + t500 * t71) * t323 - (t103 * (-t206 * t316 - t317 * t500) + t572 * t251) * qJD(4) + (t215 * t183 + t214 * t382 + (t317 * t115 - t316 * t116) * qJD(4)) * t344 + t103 * ((-t116 - t163) * t316 + (t115 + t162) * t317) + t35 * t206 - t34 * t500 + (t71 * t486 - t70 * t488) * t250 + t572 * t224 - g(1) * t251 - g(2) * t206 + g(3) * t500) * m(5); t343 + (-g(1) * t232 - g(2) * t192 + t10 * t347 + t14 * t504 + ((-t316 * t171 + t317 * t381) * t323 + t607) * t57 + t606 + (t199 * t317 + t608) * t58 + t611 * t503) * m(6);];
tau = t1;
