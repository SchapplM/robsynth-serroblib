% Calculate vector of inverse dynamics joint torques for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:42
% DurationCPUTime: 12.33s
% Computational Cost: add. (13388->533), mult. (10177->626), div. (0->0), fcn. (7800->8), ass. (0->304)
t287 = qJ(1) + qJ(2);
t274 = pkin(8) + t287;
t268 = sin(t274);
t288 = sin(qJ(4));
t290 = cos(qJ(4));
t235 = Icges(5,5) * t290 - Icges(5,6) * t288;
t269 = cos(t274);
t324 = t235 * t269;
t133 = Icges(5,3) * t268 + t324;
t237 = Icges(6,4) * t290 + Icges(6,6) * t288;
t325 = t237 * t269;
t135 = Icges(6,2) * t268 + t325;
t561 = t133 + t135;
t278 = Icges(6,5) * t288;
t347 = Icges(6,1) * t290 + t278;
t138 = -Icges(6,4) * t269 + t268 * t347;
t429 = t268 * t288;
t223 = Icges(5,4) * t429;
t428 = t268 * t290;
t140 = Icges(5,1) * t428 - Icges(5,5) * t269 - t223;
t560 = -t138 - t140;
t327 = t347 * t269;
t139 = Icges(6,4) * t268 + t327;
t447 = Icges(5,4) * t288;
t243 = Icges(5,1) * t290 - t447;
t328 = t243 * t269;
t141 = Icges(5,5) * t268 + t328;
t553 = t139 + t141;
t238 = Icges(5,2) * t290 + t447;
t345 = Icges(6,3) * t290 - t278;
t559 = t238 + t345;
t446 = Icges(6,5) * t290;
t240 = Icges(6,1) * t288 - t446;
t279 = Icges(5,4) * t290;
t242 = Icges(5,1) * t288 + t279;
t558 = t240 + t242;
t233 = Icges(6,3) * t288 + t446;
t346 = -Icges(5,2) * t288 + t279;
t557 = t233 - t346;
t556 = t243 + t347;
t425 = t269 * t290;
t222 = Icges(6,5) * t425;
t426 = t269 * t288;
t131 = Icges(6,6) * t268 + Icges(6,3) * t426 + t222;
t555 = t131 * t426 + t268 * t561 + t425 * t553;
t134 = -Icges(6,2) * t269 + t237 * t268;
t122 = t268 * t134;
t130 = -Icges(6,6) * t269 + t233 * t268;
t132 = Icges(5,5) * t428 - Icges(5,6) * t429 - Icges(5,3) * t269;
t554 = -t130 * t426 - t268 * t132 + t425 * t560 - t122;
t234 = Icges(5,5) * t288 + Icges(5,6) * t290;
t236 = Icges(6,4) * t288 - Icges(6,6) * t290;
t552 = t234 + t236;
t286 = qJD(1) + qJD(2);
t551 = (-Icges(5,6) + Icges(6,6)) * t286 + t559 * qJD(4);
t550 = (Icges(6,4) + Icges(5,5)) * t286 - t558 * qJD(4);
t546 = -t288 * t559 + t290 * t558;
t136 = Icges(5,4) * t428 - Icges(5,2) * t429 - Icges(5,6) * t269;
t439 = t136 * t288;
t340 = -t140 * t290 + t439;
t440 = t134 * t269;
t343 = t130 * t288 + t138 * t290;
t489 = t268 * t343;
t50 = -t440 + t489;
t495 = -t132 * t269 - t268 * t340 + t50;
t493 = -t136 * t426 - t554;
t326 = t346 * t269;
t137 = Icges(5,6) * t268 + t326;
t492 = -t137 * t426 + t555;
t359 = -t131 * t429 + t135 * t269 - t139 * t428;
t113 = t141 * t428;
t364 = t133 * t269 - t113;
t53 = -t137 * t429 - t364;
t494 = t53 - t359;
t549 = t557 * qJD(4);
t548 = t556 * qJD(4);
t547 = t235 + t237;
t545 = -t288 * t558 - t290 * t559;
t524 = rSges(6,1) + pkin(4);
t430 = t268 * t286;
t544 = t269 * t551 - t430 * t557;
t427 = t269 * t286;
t543 = -t233 * t427 - t268 * t551 + t286 * t326;
t542 = t269 * t550 - t430 * t556;
t541 = (t327 + t328) * t286 + t550 * t268;
t531 = rSges(6,3) + qJ(5);
t540 = (t131 - t137) * t290 - t553 * t288;
t491 = (t130 - t136) * t290 + t560 * t288;
t486 = t531 * t428;
t432 = t236 * t269;
t435 = t234 * t269;
t516 = t268 * t546 - t432 - t435;
t433 = t236 * t268;
t436 = t234 * t268;
t515 = t269 * t546 + t433 + t436;
t539 = (-Icges(6,2) - Icges(5,3)) * t286 + t552 * qJD(4);
t438 = t137 * t288;
t538 = t131 * t288 + t290 * t553 - t438;
t537 = -t340 + t343;
t277 = t288 * qJ(5);
t357 = rSges(6,1) * t290 + rSges(6,3) * t288;
t536 = pkin(4) * t290 + t277 + t357;
t535 = qJD(4) * t545 + t286 * t552 + t288 * t549 + t290 * t548;
t534 = t492 * t268 - t269 * t493;
t533 = t268 * t494 - t495 * t269;
t532 = -qJD(4) * t547 + t546 * t286;
t530 = t515 * t286;
t262 = t269 * rSges(6,2);
t411 = t268 * t536 - t262;
t399 = t288 * t524 - t290 * t531;
t529 = t540 * qJD(4) + t286 * t561 + t544 * t288 + t542 * t290;
t528 = -t541 * t290 + t543 * t288 + (-t132 - t134) * t286 - t491 * qJD(4);
t527 = t516 * t286;
t526 = (-t325 - t324 + t537) * t286 + t539 * t268;
t525 = t269 * t539 + t286 * t538 + t430 * t547;
t523 = qJD(4) * t533 + t527;
t522 = qJD(4) * t534 + t530;
t521 = qJD(4) * t537 + t288 * t541 + t290 * t543;
t289 = sin(qJ(1));
t291 = cos(qJ(1));
t292 = qJD(1) ^ 2;
t323 = (-qJDD(1) * t289 - t291 * t292) * pkin(1);
t520 = t323 - g(1);
t519 = qJD(4) * t538 + t288 * t542 - t290 * t544;
t518 = -t268 * t532 + t269 * t535;
t517 = t268 * t535 + t269 * t532;
t229 = rSges(5,1) * t425;
t484 = rSges(5,3) * t268 + t229;
t145 = -rSges(5,2) * t426 + t484;
t266 = t269 * pkin(3);
t191 = pkin(7) * t268 + t266;
t276 = cos(t287);
t270 = pkin(2) * t276;
t483 = t270 + t191;
t110 = t145 + t483;
t265 = t269 * pkin(7);
t190 = pkin(3) * t268 - t265;
t215 = pkin(7) * t427;
t514 = t190 * t286 + t215;
t394 = qJD(4) * t288;
t380 = t268 * t394;
t513 = t524 * t380;
t512 = rSges(6,2) * t268 + pkin(4) * t425;
t511 = t440 + t555;
t393 = qJD(4) * t290;
t379 = t268 * t393;
t422 = t286 * t288;
t510 = t269 * t422 + t379;
t454 = pkin(1) * qJD(1);
t387 = t289 * t454;
t275 = sin(t287);
t194 = rSges(3,1) * t275 + rSges(3,2) * t276;
t437 = t194 * t286;
t150 = -t387 - t437;
t424 = t275 * t286;
t389 = pkin(2) * t424;
t509 = t389 - t387;
t410 = rSges(6,1) * t425 + rSges(6,3) * t426 + t269 * t277 + t512;
t74 = t483 + t410;
t392 = qJD(5) * t288;
t216 = t269 * t392;
t396 = qJD(4) * t269;
t508 = t396 * t399 - t216;
t507 = t268 * t526 + t269 * t528;
t506 = -t268 * t525 + t269 * t529;
t152 = t191 * t286;
t395 = qJD(4) * t286;
t176 = -qJDD(4) * t269 + t268 * t395;
t285 = qJDD(1) + qJDD(2);
t284 = t286 ^ 2;
t423 = t276 * t284;
t305 = -pkin(2) * t423 + t323;
t391 = qJD(5) * t290;
t407 = -qJD(4) * t536 + t391;
t312 = qJDD(5) * t288 + (t391 + t407) * qJD(4);
t459 = pkin(2) * t275;
t367 = -t190 - t459;
t333 = t367 - t411;
t376 = t268 * t392;
t482 = t376 - t513;
t456 = rSges(6,3) * t379 + t510 * qJ(5) + t482 + (t269 * t357 + t512) * t286;
t14 = t399 * t176 + (-t152 - t376 - t456) * t286 + t333 * t285 + t312 * t269 + t305;
t505 = t14 - g(1);
t175 = qJDD(4) * t268 + t269 * t395;
t283 = t291 * pkin(1);
t460 = pkin(1) * t289;
t362 = qJDD(1) * t283 - t292 * t460;
t316 = t270 * t285 - t284 * t459 + t362;
t304 = t286 * (-pkin(3) * t430 + t215) + t285 * t191 + t316;
t378 = t269 * t394;
t320 = -t286 * t428 - t378;
t386 = t268 * t422;
t377 = t269 * t393;
t481 = rSges(6,2) * t427 + t377 * t531 + t216;
t457 = t320 * t524 - t386 * t531 + t481;
t15 = t410 * t285 - t399 * t175 + (t216 + t457) * t286 + t312 * t268 + t304;
t504 = t15 - g(2);
t455 = rSges(5,1) * t290;
t251 = -rSges(5,2) * t288 + t455;
t203 = t251 * qJD(4);
t247 = rSges(5,1) * t288 + rSges(5,2) * t290;
t404 = rSges(5,2) * t429 + rSges(5,3) * t269;
t143 = rSges(5,1) * t428 - t404;
t361 = -t143 + t367;
t384 = rSges(5,1) * t380 + rSges(5,2) * t510;
t96 = t286 * t484 - t384;
t26 = -t203 * t396 + t176 * t247 + (-t152 - t96) * t286 + t361 * t285 + t305;
t503 = t26 - g(1);
t397 = qJD(4) * t268;
t329 = rSges(5,3) * t427 + (-t377 + t386) * rSges(5,2);
t94 = rSges(5,1) * t320 + t329;
t27 = t145 * t285 - t175 * t247 - t203 * t397 + t286 * t94 + t304;
t502 = t27 - g(2);
t188 = rSges(4,1) * t268 + rSges(4,2) * t269;
t212 = rSges(4,2) * t430;
t501 = -t285 * t188 - t286 * (rSges(4,1) * t427 - t212) + (-t275 * t285 - t423) * pkin(2) + t520;
t263 = t269 * rSges(4,1);
t189 = -rSges(4,2) * t268 + t263;
t500 = -t188 * t284 + t189 * t285 - g(2) + t316;
t499 = t268 * t528 - t269 * t526;
t267 = t276 * rSges(3,1);
t180 = -rSges(3,2) * t424 + t267 * t286;
t498 = -t180 * t286 - t194 * t285 + t520;
t195 = -rSges(3,2) * t275 + t267;
t497 = t195 * t285 - t286 * t437 - g(2) + t362;
t496 = t268 * t529 + t269 * t525;
t487 = t189 + t270;
t485 = t531 * t425;
t127 = t286 * t143;
t480 = -rSges(5,1) * t378 + t127 + t329 + t514;
t475 = t110 * t286 - t247 * t397;
t472 = t286 * t411 + t481 + t514;
t471 = t286 * t74 - t397 * t399;
t413 = -Icges(5,2) * t428 + t140 - t223;
t417 = t242 * t268 + t136;
t470 = -t288 * t413 - t290 * t417;
t415 = -t268 * t345 + t138;
t419 = -t240 * t268 + t130;
t469 = -t288 * t415 + t290 * t419;
t468 = t175 / 0.2e1;
t467 = t176 / 0.2e1;
t458 = g(2) * t268;
t315 = -t387 - t508;
t48 = t286 * t333 + t315;
t453 = t269 * t48;
t388 = t291 * t454;
t321 = t376 + t388;
t49 = t321 + t471;
t452 = t275 * t49;
t47 = -t391 + qJD(3) + (t268 * t411 + t269 * t410) * qJD(4);
t451 = t47 * t288;
t434 = t235 * t286;
t431 = t237 * t286;
t418 = -Icges(6,1) * t426 + t131 + t222;
t416 = -t242 * t269 - t137;
t414 = -t269 * t345 + t139;
t412 = -t238 * t269 + t141;
t409 = -t429 * t524 + t486;
t408 = -t426 * t524 + t485;
t403 = -t345 + t347;
t402 = t233 - t240;
t401 = -t238 + t243;
t400 = t242 + t346;
t63 = t388 + t475;
t390 = t63 * t459;
t382 = t269 * t524;
t381 = t247 * t396;
t373 = -pkin(3) - t455;
t372 = -t397 / 0.2e1;
t371 = t397 / 0.2e1;
t370 = -t396 / 0.2e1;
t369 = t396 / 0.2e1;
t153 = -t188 - t459;
t365 = t265 - t459;
t363 = -t132 + t438;
t252 = rSges(2,1) * t291 - rSges(2,2) * t289;
t248 = rSges(2,1) * t289 + rSges(2,2) * t291;
t355 = t268 * t49 + t453;
t322 = -t381 - t387;
t62 = t286 * t361 + t322;
t350 = -t268 * t63 - t269 * t62;
t338 = t143 * t268 + t145 * t269;
t330 = t355 * t290;
t319 = -t288 * t414 + t290 * t418;
t318 = -t288 * t412 + t290 * t416;
t317 = -t288 * t531 - t290 * t524 - pkin(3);
t314 = (t288 * t402 + t290 * t403) * t286;
t313 = (-t288 * t400 + t290 * t401) * t286;
t109 = t268 * t373 + t365 + t404;
t73 = t268 * t317 + t262 + t365;
t124 = t153 * t286 - t387;
t125 = t286 * t487 + t388;
t296 = (t124 * (-t263 - t270) + t125 * t153) * t286;
t295 = (t62 * (-t229 - t266 - t270) - t390 + (t62 * (-rSges(5,3) - pkin(7)) + t63 * t373) * t268) * t286;
t294 = (((t53 - t113 + (t133 + t439) * t269 + t554) * t269 + (t50 - t489 + t511) * t268) * qJD(4) + t530) * t369 + (qJD(4) * t546 + t288 * t548 - t290 * t549) * t286 + (t515 - t540) * t468 + (t516 - t491) * t467 + (t518 + t519) * t371 + (Icges(3,3) + Icges(4,3) - t545) * t285 + (((t269 * t363 + t492 - t511) * t269 + (t268 * t363 - t122 + t359 + t364 + t493) * t268) * qJD(4) + t523 - t527) * t372 + (t517 + t521 + t522) * t370;
t293 = (-t288 * t382 * t49 - t48 * t486) * qJD(4) + ((-t276 * t48 - t452) * pkin(2) + t317 * t453 + (t48 * (-rSges(6,2) - pkin(7)) + t49 * (-pkin(3) - t536)) * t268) * t286;
t178 = t286 * t188;
t173 = t247 * t269;
t169 = t247 * t268;
t151 = t195 * t286 + t388;
t70 = qJD(4) * t338 + qJD(3);
t25 = t143 * t175 - t145 * t176 + qJDD(3) + (t268 * t96 + t269 * t94) * qJD(4);
t5 = -qJDD(5) * t290 + qJDD(3) - t410 * t176 + t411 * t175 + (t268 * t456 + t269 * t457 + t392) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t294 + (t497 * (t195 + t283) + t498 * (-t194 - t460) + (-t180 - t388 + t151) * t150) * m(3) + ((t248 ^ 2 + t252 ^ 2) * qJDD(1) + g(1) * t248 - g(2) * t252) * m(2) + (t48 * (-t321 + t513) + t293 + t504 * (t283 + t74) + t505 * (t73 - t460) + (-t315 + t48 + t472 + t509) * t49) * m(6) + (t62 * (t384 - t388) + t295 + (-t322 + t62 + t480 + t509) * t63 + t502 * (t110 + t283) + t503 * (t109 - t460)) * m(5) + (t124 * (t212 - t388) + t296 + t500 * (t487 + t283) + t501 * (t153 - t460) + (t124 + t178 + t389) * t125) * m(4); t294 + (pkin(2) * t452 * t286 + t293 + t504 * t74 + t505 * t73 + (t472 + t508) * t49 + (t376 + t471 - t482) * t48) * m(6) + (t390 * t286 + t295 + (t381 + t480) * t63 + (t384 + t475) * t62 + t502 * t110 + t503 * t109) * m(5) + (t124 * t212 + t296 + t125 * t178 - (-t124 * t487 - t125 * t459) * t286 + t500 * t487 + t501 * t153) * m(4) + (-t150 * t180 - t151 * t437 + (t150 * t286 + t497) * t195 + (t151 * t286 - t498) * t194) * m(3); m(4) * qJDD(3) + m(5) * t25 + m(6) * t5 + (-m(4) - m(5) - m(6)) * g(3); t534 * t468 + t533 * t467 + (t518 * t286 + t515 * t285 + t493 * t176 + t492 * t175 + (t506 * t268 + t507 * t269) * qJD(4)) * t268 / 0.2e1 - (t517 * t286 + t516 * t285 + t495 * t176 + t494 * t175 + (t496 * t268 + t499 * t269) * qJD(4)) * t269 / 0.2e1 + (-t268 * t540 + t491 * t269) * t285 / 0.2e1 - (((t400 - t402) * t290 + (t401 + t403) * t288) * t286 + (((-t413 - t415) * t269 + (t412 + t414) * t268) * t290 + ((t417 - t419) * t269 + (t416 + t418) * t268) * t288) * qJD(4)) * t286 / 0.2e1 + ((-t286 * t540 - t521) * t269 + (-t286 * t491 + t519) * t268) * t286 / 0.2e1 + t523 * t430 / 0.2e1 + t522 * t427 / 0.2e1 + ((-t397 * t432 + t431) * t268 + (t314 + (-t469 * t269 + (t433 + t319) * t268) * qJD(4)) * t269 + (-t397 * t435 + t434) * t268 + (t313 + (-t470 * t269 + (t436 + t318) * t268) * qJD(4)) * t269) * t372 + ((t286 * t492 + t507) * t269 + (t286 * t493 + t506) * t268) * t371 + ((t286 * t494 + t499) * t269 + (t286 * t495 + t496) * t268) * t370 + ((-t396 * t433 - t431) * t269 + (t314 + (t319 * t268 + (t432 - t469) * t269) * qJD(4)) * t268 + (-t396 * t436 - t434) * t269 + (t313 + (t318 * t268 + (t435 - t470) * t269) * qJD(4)) * t268) * t369 + ((-t14 * t399 + t48 * t407 + t5 * t410 + t47 * t457 + (-t399 * t49 + t411 * t47) * t286) * t269 + (-t15 * t399 + t49 * t407 + t5 * t411 + t47 * t456 + (t399 * t48 - t410 * t47) * t286) * t268 - (t330 + t451) * qJD(5) - (t408 * t49 - t409 * t48) * t286 - ((t408 * t47 - t48 * t536) * t269 + (t409 * t47 - t49 * t536) * t268) * qJD(4) - g(1) * t485 - g(2) * t486 - g(3) * t536 - (-g(1) * t382 - t458 * t524) * t288) * m(6) + (t25 * t338 + t70 * ((t94 + t127) * t269 + (-t145 * t286 + t96) * t268) + t350 * t203 + ((-t286 * t63 - t26) * t269 + (t286 * t62 - t27) * t268) * t247 - (t169 * t62 - t173 * t63) * t286 - (t70 * (-t169 * t268 - t173 * t269) + t350 * t251) * qJD(4) + g(1) * t173 + g(2) * t169 - g(3) * t251) * m(5); (-(-t268 * t48 + t269 * t49) * t422 - (t330 + (t268 ^ 2 + t269 ^ 2) * t451) * qJD(4) + (qJD(4) * t355 + g(3) - t5) * t290 + (qJD(4) * t47 + (-t286 * t48 + t15) * t268 - t458 + (t286 * t49 + t505) * t269) * t288) * m(6);];
tau = t1;
