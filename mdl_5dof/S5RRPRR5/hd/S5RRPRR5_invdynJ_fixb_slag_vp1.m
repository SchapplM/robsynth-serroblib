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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:03:31
% EndTime: 2020-01-03 12:03:46
% DurationCPUTime: 10.20s
% Computational Cost: add. (19881->678), mult. (13225->843), div. (0->0), fcn. (10283->10), ass. (0->377)
t320 = qJD(1) + qJD(2);
t321 = qJ(1) + qJ(2);
t310 = sin(t321);
t311 = cos(t321);
t549 = t310 * rSges(3,1) + t311 * rSges(3,2);
t208 = t549 * t320;
t325 = sin(qJ(1));
t509 = pkin(1) * qJD(1);
t418 = t325 * t509;
t193 = t418 + t208;
t318 = pkin(9) + qJ(4);
t309 = qJ(5) + t318;
t298 = sin(t309);
t299 = cos(t309);
t500 = Icges(6,4) * t298;
t227 = Icges(6,1) * t299 - t500;
t361 = t227 * t310;
t163 = -Icges(6,5) * t311 + t361;
t476 = t298 * t311;
t255 = Icges(6,4) * t476;
t474 = t299 * t311;
t164 = Icges(6,1) * t474 + Icges(6,5) * t310 - t255;
t224 = Icges(6,2) * t299 + t500;
t319 = qJD(4) + qJD(5);
t243 = t310 * t319;
t467 = t311 * t319;
t290 = Icges(6,4) * t299;
t374 = -Icges(6,2) * t298 + t290;
t551 = Icges(6,1) * t298 + t290;
t566 = t551 + t374;
t333 = t243 * (-Icges(6,2) * t474 + t164 - t255) - t467 * (-t224 * t310 + t163) + t320 * t566;
t359 = t374 * t310;
t161 = -Icges(6,6) * t311 + t359;
t162 = Icges(6,4) * t474 - Icges(6,2) * t476 + Icges(6,6) * t310;
t530 = t243 * (t311 * t551 + t162) - t467 * (t310 * t551 + t161) + t320 * (t224 - t227);
t569 = t333 * t298 + t299 * t530;
t307 = sin(t318);
t308 = cos(t318);
t501 = Icges(5,4) * t307;
t240 = Icges(5,1) * t308 - t501;
t362 = t240 * t310;
t174 = -Icges(5,5) * t311 + t362;
t472 = t307 * t311;
t265 = Icges(5,4) * t472;
t470 = t308 * t311;
t175 = Icges(5,1) * t470 + Icges(5,5) * t310 - t265;
t237 = Icges(5,2) * t308 + t501;
t343 = t310 * (-Icges(5,2) * t470 + t175 - t265) - t311 * (-t237 * t310 + t174);
t292 = Icges(5,4) * t308;
t375 = -Icges(5,2) * t307 + t292;
t172 = -Icges(5,6) * t311 + t310 * t375;
t173 = Icges(5,4) * t470 - Icges(5,2) * t472 + Icges(5,6) * t310;
t534 = -t311 * t172 + t173 * t310;
t568 = -t343 * t307 - t308 * t534;
t550 = Icges(5,1) * t307 + t292;
t443 = t550 + t375;
t444 = t237 - t240;
t567 = (t307 * t443 + t308 * t444) * t320;
t294 = pkin(4) * t308;
t323 = cos(pkin(9));
t300 = t323 * pkin(3) + pkin(2);
t252 = t294 + t300;
t466 = t311 * t320;
t207 = t252 * t466;
t232 = t300 * t466;
t431 = qJD(4) * t307;
t417 = pkin(4) * t431;
t324 = -pkin(7) - qJ(3);
t317 = -pkin(8) + t324;
t434 = t317 - t324;
t102 = t207 - t232 + (-t320 * t434 - t417) * t310;
t427 = qJD(4) * t320;
t271 = t310 * t427;
t426 = -qJDD(4) - qJDD(5);
t468 = t310 * t320;
t154 = qJD(5) * t468 + t311 * t426 + t271;
t291 = t299 * rSges(6,1);
t511 = rSges(6,2) * t298;
t229 = t291 - t511;
t198 = t229 * t319;
t214 = -qJDD(4) * t311 + t271;
t510 = rSges(6,2) * t299;
t514 = rSges(6,1) * t298;
t228 = t510 + t514;
t316 = qJDD(1) + qJDD(2);
t301 = t310 * pkin(2);
t245 = -qJ(3) * t311 + t301;
t289 = t311 * t324;
t439 = t310 * t300 + t289;
t157 = -t245 + t439;
t433 = qJD(3) * t311;
t326 = cos(qJ(1));
t315 = t326 * pkin(1);
t328 = qJD(1) ^ 2;
t491 = pkin(1) * qJDD(1);
t435 = t328 * t315 + t325 * t491;
t438 = pkin(2) * t466 + qJ(3) * t468;
t356 = -qJDD(3) * t310 + t320 * (-t433 + t438) + t316 * t245 + t435;
t465 = t320 * t324;
t341 = t320 * (-t310 * t465 + t232 - t438) + t316 * t157 + t356;
t451 = t310 * t252 + t311 * t317;
t132 = -t439 + t451;
t475 = t299 * t310;
t477 = t298 * t310;
t165 = rSges(6,1) * t475 - rSges(6,2) * t477 - rSges(6,3) * t311;
t464 = t132 + t165;
t469 = t308 * qJD(4) ^ 2;
t423 = t319 * t514;
t422 = rSges(6,1) * t474;
t447 = rSges(6,3) * t468 + t320 * t422;
t99 = -t310 * t423 + (-t243 * t299 - t298 * t466) * rSges(6,2) + t447;
t13 = -t154 * t228 + t198 * t467 + t464 * t316 + (-t214 * t307 + t311 * t469) * pkin(4) + (t102 + t99 - t433) * t320 + t341;
t565 = t13 - g(3);
t153 = t310 * t426 - t320 * t467;
t213 = -qJDD(4) * t310 - t311 * t427;
t314 = t325 * pkin(1);
t387 = -t314 * t328 + t326 * t491;
t432 = qJD(3) * t320;
t349 = -qJDD(3) * t311 + t310 * t432 + t387;
t247 = pkin(2) * t311 + qJ(3) * t310;
t436 = -t311 * t300 + t310 * t324;
t158 = t247 + t436;
t442 = t247 - t158;
t211 = t311 * t252;
t133 = t310 * t317 - t211 - t436;
t386 = -rSges(6,2) * t476 + t422;
t166 = rSges(6,3) * t310 + t386;
t463 = -t133 + t166;
t388 = t442 + t463;
t274 = qJ(3) * t466;
t293 = qJD(3) * t310;
t437 = t274 + t293;
t181 = pkin(2) * t468 - t437;
t462 = -t274 - (t289 + (-pkin(2) + t300) * t310) * t320 - t181;
t428 = qJD(4) * t311;
t410 = t307 * t428;
t390 = pkin(4) * t410;
t101 = t390 + (t434 * t311 + (t252 - t300) * t310) * t320;
t260 = rSges(6,2) * t474;
t419 = t320 * t511;
t448 = rSges(6,3) * t466 + t310 * t419;
t98 = t319 * t260 + (t298 * t467 + t299 * t468) * rSges(6,1) - t448;
t504 = -t101 - t98;
t14 = t153 * t228 - t198 * t243 + (t213 * t307 - t310 * t469) * pkin(4) + (t462 + t504) * t320 + t388 * t316 + t349;
t564 = t14 - g(2);
t429 = qJD(4) * t310;
t409 = t308 * t429;
t411 = t307 * t429;
t424 = rSges(5,1) * t470;
t445 = rSges(5,3) * t468 + t320 * t424;
t116 = -rSges(5,1) * t411 + (-t307 * t466 - t409) * rSges(5,2) + t445;
t471 = t308 * t310;
t473 = t307 * t310;
t178 = rSges(5,1) * t471 - rSges(5,2) * t473 - rSges(5,3) * t311;
t512 = rSges(5,2) * t307;
t515 = rSges(5,1) * t308;
t242 = -t512 + t515;
t221 = t242 * qJD(4);
t557 = rSges(5,2) * t308;
t241 = rSges(5,1) * t307 + t557;
t34 = t116 * t320 + t178 * t316 - t214 * t241 + (qJD(4) * t221 - t432) * t311 + t341;
t563 = t34 - g(3);
t420 = t320 * t512;
t446 = rSges(5,3) * t466 + t310 * t420;
t115 = t428 * t557 + (t308 * t468 + t410) * rSges(5,1) - t446;
t179 = -rSges(5,2) * t472 + rSges(5,3) * t310 + t424;
t414 = t179 + t442;
t35 = -t221 * t429 + t213 * t241 + (-t115 + t462) * t320 + t414 * t316 + t349;
t562 = t35 - g(2);
t516 = rSges(4,1) * t323;
t285 = t310 * t516;
t322 = sin(pkin(9));
t513 = rSges(4,2) * t322;
t385 = -t310 * t513 + t285;
t188 = -rSges(4,3) * t311 + t385;
t421 = t320 * t513;
t425 = t311 * t516;
t342 = t320 * t425 + rSges(4,3) * t468 + (-qJD(3) - t421) * t311;
t62 = t316 * t188 + t320 * t342 + t356;
t561 = t62 - g(3);
t440 = rSges(4,3) * t466 + t310 * t421;
t286 = t311 * t513;
t189 = rSges(4,3) * t310 - t286 + t425;
t441 = t247 + t189;
t63 = (-t285 * t320 - t181 + t440) * t320 + t441 * t316 + t349;
t560 = t63 - g(2);
t248 = rSges(3,1) * t311 - t310 * rSges(3,2);
t559 = -t208 * t320 + t248 * t316 - g(2) + t387;
t209 = rSges(3,1) * t466 - rSges(3,2) * t468;
t558 = t209 * t320 + t316 * t549 - g(3) + t435;
t223 = Icges(6,5) * t299 - Icges(6,6) * t298;
t159 = -Icges(6,3) * t311 + t223 * t310;
t486 = t162 * t298;
t373 = -t164 * t299 + t486;
t365 = -t159 + t373;
t556 = t467 * t365;
t555 = t320 * (t188 + t245);
t217 = t320 * t247;
t552 = -t320 * t158 + t217;
t548 = -t320 * t189 - t217 + t342 + t433 + t438;
t156 = t320 * t179;
t364 = -t241 * t429 - t433;
t547 = t232 + t445 - t156 - t364 - t552;
t152 = t320 * t166;
t339 = -pkin(4) * t411 - t228 * t243 - t433;
t546 = t320 * t133 - t152 + t207 - t339 + t447 - t552;
t479 = t224 * t319;
t545 = -Icges(6,6) * t320 + t479;
t190 = t228 * t310;
t191 = rSges(6,1) * t476 + t260;
t57 = t165 * t243 + t166 * t467 + (t132 * t310 - t133 * t311) * qJD(4);
t363 = t228 * t467 - t293 + t390;
t461 = t157 + t245;
t389 = t461 + t464;
t58 = t320 * t389 + t363 + t418;
t306 = t326 * t509;
t59 = t320 * t388 + t306 + t339;
t544 = -(-t320 * t190 + t229 * t467) * t58 - t57 * (-t190 * t243 - t191 * t467) - t59 * (-t191 * t320 - t243 * t229);
t222 = Icges(6,5) * t298 + Icges(6,6) * t299;
t543 = -Icges(6,3) * t320 + t222 * t319;
t542 = -Icges(6,5) * t320 + t319 * t551;
t104 = t172 * t308 + t174 * t307;
t360 = t375 * t320;
t539 = -Icges(5,6) * t320 + qJD(4) * t237;
t112 = -t310 * t539 + t311 * t360;
t536 = -Icges(5,5) * t320 + qJD(4) * t550;
t114 = t240 * t466 - t310 * t536;
t236 = Icges(5,5) * t308 - Icges(5,6) * t307;
t170 = -Icges(5,3) * t311 + t236 * t310;
t541 = qJD(4) * t104 + t112 * t307 - t114 * t308 - t170 * t320;
t235 = Icges(5,5) * t307 + Icges(5,6) * t308;
t540 = -Icges(5,3) * t320 + qJD(4) * t235;
t219 = t375 * qJD(4);
t220 = t240 * qJD(4);
t367 = t237 * t308 + t307 * t550;
t538 = qJD(4) * t367 + t219 * t307 - t220 * t308 - t235 * t320;
t111 = t310 * t360 + t311 * t539;
t113 = t311 * t536 + t320 * t362;
t171 = Icges(5,5) * t470 - Icges(5,6) * t472 + Icges(5,3) * t310;
t371 = t173 * t308 + t175 * t307;
t537 = qJD(4) * t371 - t111 * t307 + t113 * t308 - t171 * t320;
t391 = -t227 * t319 + t479;
t392 = t566 * t319;
t533 = -t222 * t320 + t298 * t392 + t299 * t391;
t160 = Icges(6,5) * t474 - Icges(6,6) * t476 + Icges(6,3) * t310;
t396 = t164 * t319 - t311 * t545 - t320 * t359;
t398 = t162 * t319 + t311 * t542 + t320 * t361;
t532 = -t160 * t320 + t298 * t396 + t299 * t398;
t397 = t163 * t319 - t310 * t545 + t374 * t466;
t399 = t161 * t319 - t227 * t466 + t310 * t542;
t531 = -t159 * t320 + t298 * t397 + t299 * t399;
t529 = t153 / 0.2e1;
t528 = t154 / 0.2e1;
t527 = t213 / 0.2e1;
t526 = t214 / 0.2e1;
t525 = -t243 / 0.2e1;
t524 = t243 / 0.2e1;
t523 = t467 / 0.2e1;
t522 = -t467 / 0.2e1;
t521 = -t310 / 0.2e1;
t520 = -t311 / 0.2e1;
t519 = t316 / 0.2e1;
t518 = -t320 / 0.2e1;
t517 = t320 / 0.2e1;
t508 = t320 * t59;
t71 = t320 * t414 + t306 + t364;
t507 = t320 * t71;
t481 = t222 * t310;
t89 = -t224 * t476 + t474 * t551 + t481;
t506 = t89 * t320;
t505 = rSges(4,3) + qJ(3);
t478 = t235 * t310;
t107 = -t237 * t472 + t470 * t550 + t478;
t490 = t107 * t320;
t487 = t161 * t298;
t485 = t163 * t299;
t484 = t172 * t307;
t483 = t173 * t307;
t482 = t174 * t308;
t183 = t222 * t311;
t357 = t223 * t320;
t200 = t235 * t311;
t358 = t236 * t320;
t430 = qJD(4) * t308;
t416 = t165 * t466 + (-t152 + t99) * t310;
t415 = t178 + t461;
t407 = t468 / 0.2e1;
t406 = -t466 / 0.2e1;
t405 = pkin(2) + t516;
t404 = -t429 / 0.2e1;
t403 = t429 / 0.2e1;
t402 = -t428 / 0.2e1;
t401 = t428 / 0.2e1;
t400 = -pkin(4) * t307 - t228;
t395 = -t160 - t487;
t394 = -t160 + t485;
t194 = t248 * t320 + t306;
t384 = -t293 + t418;
t383 = -t241 * t428 + t293;
t284 = rSges(2,1) * t326 - t325 * rSges(2,2);
t283 = rSges(2,1) * t325 + rSges(2,2) * t326;
t70 = t320 * t415 - t383 + t418;
t380 = -t310 * t71 + t311 * t70;
t144 = t174 * t471;
t73 = -t170 * t311 - t172 * t473 + t144;
t145 = t175 * t471;
t74 = t171 * t311 + t173 * t473 - t145;
t379 = -t310 * t74 - t311 * t73;
t146 = t172 * t472;
t75 = -t170 * t310 - t174 * t470 + t146;
t370 = -t175 * t308 + t483;
t76 = t310 * t171 - t370 * t311;
t378 = -t310 * t76 - t311 * t75;
t91 = -t162 * t299 - t164 * t298;
t372 = t482 - t484;
t369 = t178 * t310 + t179 * t311;
t368 = -t224 * t298 + t299 * t551;
t366 = -t237 * t307 + t308 * t550;
t353 = t183 * t243 - t467 * t481 - t357;
t351 = -t310 * t357 - t311 * t543 + t320 * t373;
t350 = -t311 * t357 + t543 * t310 + (t485 - t487) * t320;
t348 = -t310 * t358 - t311 * t540 + t320 * t370;
t347 = t310 * t540 - t311 * t358 + t320 * t372;
t346 = -t223 * t319 + t320 * t368;
t345 = -t236 * qJD(4) + t320 * t366;
t131 = t179 - t436;
t130 = t178 + t439;
t121 = t165 + t451;
t122 = t211 + (rSges(6,3) - t317) * t310 + t386;
t147 = -t311 * t505 + t301 + t385;
t148 = t310 * t505 + t311 * t405 - t286;
t338 = -t405 * t468 + t437 + t440;
t15 = t350 * t310 + t311 * t531;
t16 = t351 * t310 - t311 * t532;
t17 = -t310 * t531 + t350 * t311;
t18 = t310 * t532 + t351 * t311;
t136 = t163 * t475;
t66 = -t159 * t311 - t161 * t477 + t136;
t137 = t164 * t475;
t67 = t160 * t311 + t162 * t477 - t137;
t88 = t310 * t368 - t183;
t86 = t88 * t320;
t30 = -t243 * t67 - t467 * t66 + t86;
t138 = t161 * t476;
t68 = -t159 * t310 - t163 * t474 + t138;
t69 = t160 * t310 - t311 * t373;
t31 = -t243 * t69 - t467 * t68 - t506;
t44 = t346 * t310 + t311 * t533;
t45 = -t310 * t533 + t346 * t311;
t46 = -t298 * t399 + t299 * t397;
t47 = t298 * t398 - t299 * t396;
t90 = t161 * t299 + t163 * t298;
t337 = (-t15 * t467 + t153 * t69 + t154 * t68 - t16 * t243 - t316 * t89 + t320 * t44) * t521 + (t353 * t310 + t569 * t311) * t524 + (-t569 * t310 + t353 * t311) * t523 + (t153 * t67 + t154 * t66 - t17 * t467 - t18 * t243 + t316 * t88 + t320 * t45) * t520 + (-t298 * t530 + t299 * t333) * t518 + t30 * t407 + t31 * t406 + ((-t320 * t69 - t15) * t311 + (t320 * t68 - t16) * t310) * t525 + (-t310 * t69 - t311 * t68) * t529 + (-t310 * t67 - t311 * t66) * t528 + ((-t320 * t67 - t17) * t311 + (t320 * t66 - t18) * t310) * t522 + (-t310 * t91 - t311 * t90) * t519 + ((-t320 * t91 - t46) * t311 + (t320 * t90 - t47) * t310) * t517;
t106 = t310 * t366 - t200;
t100 = t106 * t320;
t36 = qJD(4) * t379 + t100;
t37 = qJD(4) * t378 - t490;
t50 = qJD(4) * t372 + t112 * t308 + t114 * t307;
t51 = qJD(4) * t370 + t111 * t308 + t113 * t307;
t54 = t345 * t310 + t311 * t538;
t55 = -t310 * t538 + t345 * t311;
t331 = -t371 * t527 + t91 * t529 + (t100 + ((t75 + t145 - t146 + (t170 - t483) * t310) * t310 + (-t144 - t76 + (t170 - t370) * t311 + (t482 + t484) * t310) * t311) * qJD(4)) * t403 + (t86 - (t310 * t395 + t136 + t69) * t467 + (t137 - t138 + t68 + (t159 - t486) * t310) * t243 + (t243 * t394 - t556) * t311) * t524 - t153 * t89 / 0.2e1 - t213 * t107 / 0.2e1 + (t90 + t88) * t528 + (t106 + t104) * t526 + (t31 + t506 - (-t138 + t67) * t467 + (-t136 + t66) * t243 + (-t243 * t365 - t394 * t467) * t311 + (-t243 * t395 + t556) * t310) * t523 + (t46 + t45) * t522 + (t50 + t55) * t402 + (t490 + ((t146 - t74 + (t171 - t482) * t311) * t311 + (-t144 + t73 + (t171 + t484) * t310) * t310) * qJD(4) + t37) * t401 + (t47 + t44 + t30) * t525 + (t51 + t54 + t36) * t404 + (qJD(4) * t366 + t219 * t308 + t220 * t307 - t298 * t391 + t299 * t392) * t320 + (t367 + Icges(4,2) * t323 ^ 2 + (Icges(4,1) * t322 + 0.2e1 * Icges(4,4) * t323) * t322 + Icges(3,3) + t224 * t299 + t551 * t298) * t316;
t330 = (t71 * (-rSges(5,1) * t431 - rSges(5,2) * t430 - t465) + t70 * (-qJD(3) - t420)) * t311 + (-t70 * t241 * qJD(4) + (t71 * (-t300 - t515) - t70 * t324) * t320) * t310;
t329 = (t58 * (-t228 * t319 - t417) + (t59 * (-t252 - t291) - t58 * t317) * t320) * t310 + (t59 * (-t317 * t320 - t319 * t510 - t417 - t423) + t58 * (-qJD(3) - t419)) * t311;
t270 = pkin(4) * t472;
t206 = t241 * t311;
t205 = t241 * t310;
t150 = t310 * t165;
t120 = t320 * t441 + t306 - t433;
t119 = t384 + t555;
t103 = t369 * qJD(4);
t22 = t310 * t537 + t348 * t311;
t21 = -t310 * t541 + t347 * t311;
t20 = t348 * t310 - t311 * t537;
t19 = t347 * t310 + t311 * t541;
t10 = -t132 * t213 + t133 * t214 - t153 * t165 - t154 * t166 + t243 * t99 - t467 * t98 + (-t101 * t311 + t102 * t310) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t331 + (t559 * (t248 + t315) + t558 * (t314 + t549) + (-t194 + t209 + t306) * t193) * m(3) + (-g(2) * t284 - g(3) * t283 + (t283 ^ 2 + t284 ^ 2) * qJDD(1)) * m(2) + (t59 * (-t384 + t448) + t329 + t564 * (t315 + t122) + t565 * (t121 + t314) + (t59 + t546) * t58) * m(6) + (t71 * (-t384 + t446) + t330 + t562 * (t131 + t315) + t563 * (t130 + t314) + (t71 + t547) * t70) * m(5) + (t120 * (t338 - t418) + t560 * (t315 + t148) + t561 * (t314 + t147) + (t120 + t548) * t119) * m(4); t331 + (t389 * t508 + t329 + (t293 + t448 + t363) * t59 + t546 * t58 + t564 * t122 + t565 * t121) * m(6) + (t415 * t507 + t330 + (t293 + t446 - t383) * t71 + t547 * t70 + t562 * t131 + t563 * t130) * m(5) + (t560 * t148 + t561 * t147 + (-t293 + t338 + t555) * t120 + t548 * t119) * m(4) + (t193 * t209 - t194 * t208 + (-t193 * t320 + t559) * t248 + (t194 * t320 + t558) * t549) * m(3); (-m(4) - m(5) - m(6)) * (-g(2) * t311 - g(3) * t310) + m(4) * (-t310 * t62 - t311 * t63) + m(5) * (-t310 * t34 - t311 * t35) + m(6) * (-t13 * t310 - t14 * t311); (-t107 * t316 + t213 * t76 + t214 * t75 + t320 * t54 + (-t19 * t311 - t20 * t310) * qJD(4)) * t521 + (t106 * t316 + t213 * t74 + t214 * t73 + t320 * t55 + (-t21 * t311 - t22 * t310) * qJD(4)) * t520 + ((t320 * t371 - t50) * t311 + (t104 * t320 - t51) * t310) * t517 + ((-t428 * t478 - t358) * t311 + (-t567 + (t311 * t200 + t568) * qJD(4)) * t310) * t401 + t337 + ((-t307 * t444 + t308 * t443) * t320 + (-t307 * t534 + t308 * t343) * qJD(4)) * t518 + ((t200 * t429 - t358) * t310 + (t567 + (-t310 * t478 - t568) * qJD(4)) * t311) * t403 + ((-t320 * t74 - t21) * t311 + (t320 * t73 - t22) * t310) * t402 + ((-t320 * t76 - t19) * t311 + (t320 * t75 - t20) * t310) * t404 + t37 * t406 + t36 * t407 + (-t104 * t311 + t310 * t371) * t519 + t379 * t526 + t378 * t527 + (-(-t59 * t409 + ((-t310 * t58 - t311 * t59) * t320 + t57 * (-t310 ^ 2 - t311 ^ 2) * qJD(4)) * t307) * pkin(4) - g(1) * (t229 + t294) - g(3) * (t270 + t191) - g(2) * t400 * t310 + t10 * t150 + t57 * t416 + t13 * t270 + (t10 * t463 + t57 * t504 + t13 * t228 + t58 * t198 + (t132 * t57 + t400 * t59) * t320) * t311 + (t10 * t132 + t57 * t102 + t14 * t400 + t59 * (-pkin(4) * t430 - t198) + (t133 * t57 + t400 * t58) * t320) * t310 + t544) * m(6) + ((-t178 * t213 - t179 * t214 + (-t115 * t311 + t116 * t310) * qJD(4)) * t369 + t103 * ((t178 * t320 - t115) * t311 + (t116 - t156) * t310) + t380 * t221 + ((t34 - t507) * t311 + (-t320 * t70 - t35) * t310) * t241 - (-t205 * t70 - t206 * t71) * t320 - (t103 * (-t205 * t310 - t206 * t311) + t380 * t242) * qJD(4) - g(1) * t242 + g(2) * t205 - g(3) * t206) * m(5); t337 + (-g(1) * t229 + g(2) * t190 - g(3) * t191 + t10 * (t166 * t311 + t150) + t57 * (-t311 * t98 + t416) + (-t310 * t59 + t311 * t58) * t198 + ((t13 - t508) * t311 + (-t320 * t58 - t14) * t310) * t228 + t544) * m(6);];
tau = t1;
