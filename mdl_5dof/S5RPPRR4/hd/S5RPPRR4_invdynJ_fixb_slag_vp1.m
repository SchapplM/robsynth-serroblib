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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:30:27
% EndTime: 2020-01-03 11:31:12
% DurationCPUTime: 21.41s
% Computational Cost: add. (20580->824), mult. (23801->1092), div. (0->0), fcn. (22772->10), ass. (0->369)
t367 = cos(qJ(1));
t362 = sin(pkin(8));
t516 = qJ(3) * t362;
t364 = cos(pkin(8));
t542 = pkin(2) * t364;
t408 = t516 + t542;
t286 = t408 * t367;
t366 = sin(qJ(1));
t517 = qJ(2) * t366;
t305 = pkin(1) * t367 + t517;
t301 = qJD(1) * t305;
t468 = qJD(2) * t367;
t469 = qJD(1) * t367;
t470 = qJD(1) * t366;
t413 = pkin(1) * t469 + qJ(2) * t470 - t468;
t597 = -qJD(1) * t286 - t301 + t413;
t359 = pkin(9) + qJ(4);
t350 = sin(t359);
t540 = pkin(4) * t350;
t361 = sin(pkin(9));
t541 = pkin(3) * t361;
t299 = t540 + t541;
t596 = t299 * t367;
t339 = -qJD(4) * t364 + qJD(1);
t351 = cos(t359);
t216 = -Icges(5,3) * t364 + (Icges(5,5) * t351 - Icges(5,6) * t350) * t362;
t521 = Icges(5,4) * t351;
t217 = -Icges(5,6) * t364 + (-Icges(5,2) * t350 + t521) * t362;
t522 = Icges(5,4) * t350;
t218 = -Icges(5,5) * t364 + (Icges(5,1) * t351 - t522) * t362;
t503 = t367 * t350;
t505 = t366 * t351;
t260 = t364 * t503 - t505;
t507 = t364 * t367;
t261 = t350 * t366 + t351 * t507;
t509 = t362 * t367;
t377 = t216 * t509 - t217 * t260 + t218 * t261;
t595 = t377 * t339;
t365 = -pkin(6) - qJ(3);
t314 = t365 * t509;
t363 = cos(pkin(9));
t344 = pkin(3) * t363 + pkin(2);
t512 = t361 * t366;
t447 = pkin(3) * t512 + t344 * t507 - t314;
t193 = t286 - t447;
t467 = qJD(3) * t362;
t327 = t366 * t467;
t594 = qJD(1) * t193 + t327 + t597;
t357 = -pkin(7) + t365;
t292 = pkin(4) * t351 + t344;
t484 = t292 * t507 + t299 * t366;
t118 = t357 * t509 + t447 - t484;
t352 = qJ(5) + t359;
t342 = sin(t352);
t504 = t367 * t342;
t343 = cos(t352);
t506 = t366 * t343;
t240 = t364 * t504 - t506;
t241 = t342 * t366 + t343 * t507;
t146 = rSges(6,1) * t241 - rSges(6,2) * t240 + rSges(6,3) * t509;
t472 = t357 - t365;
t480 = t292 - t344;
t197 = t362 * t480 + t364 * t472;
t214 = -rSges(6,3) * t364 + (rSges(6,1) * t343 - rSges(6,2) * t342) * t362;
t360 = qJD(4) + qJD(5);
t424 = t360 * t362;
t284 = t367 * t424;
t295 = -t360 * t364 + qJD(1);
t465 = qJD(4) * t362;
t593 = t295 * t146 + t327 + (-t197 * t465 - qJD(2)) * t367 - t339 * t118 - t214 * t284;
t156 = Icges(5,5) * t261 - Icges(5,6) * t260 + Icges(5,3) * t509;
t245 = Icges(5,4) * t261;
t160 = Icges(5,2) * t260 - Icges(5,6) * t509 - t245;
t244 = Icges(5,4) * t260;
t162 = Icges(5,1) * t261 + Icges(5,5) * t509 - t244;
t67 = t364 * t156 - (t160 * t350 + t162 * t351) * t362;
t137 = Icges(6,5) * t241 - Icges(6,6) * t240 + Icges(6,3) * t509;
t222 = Icges(6,4) * t241;
t141 = Icges(6,2) * t240 - Icges(6,6) * t509 - t222;
t221 = Icges(6,4) * t240;
t143 = Icges(6,1) * t241 + Icges(6,5) * t509 - t221;
t59 = t364 * t137 - (t141 * t342 + t143 * t343) * t362;
t399 = t160 * t260 + t162 * t261;
t508 = t364 * t366;
t570 = t350 * t508;
t258 = -t351 * t367 - t570;
t259 = t364 * t505 - t503;
t510 = t362 * t366;
t155 = Icges(5,5) * t259 + Icges(5,6) * t258 + Icges(5,3) * t510;
t523 = Icges(5,4) * t259;
t158 = Icges(5,2) * t258 + Icges(5,6) * t510 + t523;
t243 = Icges(5,4) * t258;
t161 = Icges(5,1) * t259 + Icges(5,5) * t510 + t243;
t52 = t155 * t510 + t158 * t258 + t161 * t259;
t592 = -t399 + t52;
t283 = t366 * t424;
t238 = -t342 * t508 - t343 * t367;
t239 = t364 * t506 - t504;
t136 = Icges(6,5) * t239 + Icges(6,6) * t238 + Icges(6,3) * t510;
t520 = Icges(6,4) * t239;
t139 = Icges(6,2) * t238 + Icges(6,6) * t510 + t520;
t220 = Icges(6,4) * t238;
t142 = Icges(6,1) * t239 + Icges(6,5) * t510 + t220;
t47 = t136 * t510 + t139 * t238 + t142 * t239;
t587 = -t141 * t238 + t239 * t143;
t48 = -t137 * t510 - t587;
t211 = -Icges(6,3) * t364 + (Icges(6,5) * t343 - Icges(6,6) * t342) * t362;
t518 = Icges(6,4) * t343;
t212 = -Icges(6,6) * t364 + (-Icges(6,2) * t342 + t518) * t362;
t519 = Icges(6,4) * t342;
t213 = -Icges(6,5) * t364 + (Icges(6,1) * t343 - t519) * t362;
t71 = t211 * t510 + t212 * t238 + t213 * t239;
t21 = t283 * t47 - t284 * t48 + t295 * t71;
t165 = rSges(5,1) * t261 - rSges(5,2) * t260 + rSges(5,3) * t509;
t219 = -rSges(5,3) * t364 + (rSges(5,1) * t351 - rSges(5,2) * t350) * t362;
t588 = t165 * t339 + t327 + (-t219 * t465 - qJD(2)) * t367;
t585 = -t160 * t258 + t162 * t259;
t50 = t137 * t509 + t141 * t240 + t143 * t241;
t511 = t361 * t367;
t577 = -(t363 * t508 - t511) * rSges(4,1) - (-t361 * t508 - t363 * t367) * rSges(4,2);
t412 = t577 * qJD(1);
t462 = qJD(1) * qJD(2);
t345 = t366 * t462;
t422 = -qJDD(2) * t367 + t345;
t442 = t362 * t470;
t579 = -(t361 * t507 - t363 * t366) * rSges(4,2) + (t363 * t507 + t512) * rSges(4,1);
t448 = rSges(4,3) * t509 + t579;
t482 = t286 + t305;
t445 = t448 + t482;
t353 = qJD(2) * t366;
t473 = qJ(2) * t469 + t353;
t279 = pkin(1) * t470 - t473;
t328 = t367 * t467;
t489 = -t408 * t470 - t279 + t328;
t434 = qJD(1) * t467;
t459 = qJDD(3) * t362;
t564 = t366 * t459 + t367 * t434;
t57 = t445 * qJDD(1) + (-rSges(4,3) * t442 + t412 + t489) * qJD(1) + t422 + t564;
t583 = t367 * t57;
t578 = -rSges(3,1) * t508 + rSges(3,2) * t510;
t378 = t211 * t509 - t212 * t240 + t213 * t241;
t576 = t284 * t50 + t378 * t295;
t569 = -t139 * t240 + t142 * t241;
t501 = -t158 * t260 + t161 * t261;
t285 = pkin(2) * t508 + qJ(3) * t510;
t293 = t344 * t508;
t415 = -t365 * t510 + t293;
t192 = -pkin(3) * t511 - t285 + t415;
t441 = t362 * t469;
t440 = t364 * t469;
t477 = pkin(2) * t440 + qJ(3) * t441;
t452 = qJD(1) * t541;
t481 = t344 * t440 + t366 * t452;
t568 = qJD(1) * (-t365 * t441 - t477 + t481) + qJDD(1) * t192;
t194 = rSges(4,3) * t510 - t577;
t356 = t366 * pkin(1);
t303 = -qJ(2) * t367 + t356;
t483 = t285 + t303;
t125 = t194 + t483;
t563 = rSges(3,3) * t367 + t578;
t174 = rSges(6,1) * t238 - rSges(6,2) * t239;
t175 = rSges(6,1) * t240 + rSges(6,2) * t241;
t150 = qJD(1) * t238 + t241 * t360;
t151 = qJD(1) * t239 + t240 * t360;
t410 = rSges(6,1) * t151 + rSges(6,2) * t150;
t86 = rSges(6,3) * t442 + t410;
t152 = -qJD(1) * t240 - t239 * t360;
t153 = qJD(1) * t241 + t238 * t360;
t87 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t441;
t562 = -t146 * t441 - t174 * t284 - t175 * t283 + t509 * t87 + t510 * t86;
t409 = -rSges(6,1) * t342 - rSges(6,2) * t343;
t242 = t409 * t362;
t561 = t175 * t295 + t214 * t442 + t242 * t284 + t364 * t86;
t379 = t366 * (-Icges(5,2) * t259 + t161 + t243) - t367 * (Icges(5,2) * t261 - t162 + t244);
t557 = t366 * (-Icges(5,1) * t258 + t158 + t523) - t367 * (-Icges(5,1) * t260 + t160 - t245);
t236 = (-Icges(6,2) * t343 - t519) * t362;
t371 = t283 * (-Icges(6,2) * t239 + t142 + t220) - t284 * (Icges(6,2) * t241 - t143 + t221) + t295 * (t213 + t236);
t237 = (-Icges(6,1) * t342 - t518) * t362;
t556 = t283 * (-Icges(6,1) * t238 + t139 + t520) - t284 * (-Icges(6,1) * t240 + t141 - t222) + t295 * (t212 - t237);
t555 = t364 ^ 2;
t433 = qJD(1) * t465;
t456 = qJDD(4) * t362;
t275 = t366 * t456 + t367 * t433;
t461 = qJD(1) * qJD(5);
t202 = (qJDD(5) * t366 + t367 * t461) * t362 + t275;
t554 = t202 / 0.2e1;
t310 = t366 * t433;
t455 = -qJDD(4) - qJDD(5);
t203 = t310 + (t366 * t461 + t367 * t455) * t362;
t553 = t203 / 0.2e1;
t552 = t275 / 0.2e1;
t276 = -t367 * t456 + t310;
t551 = t276 / 0.2e1;
t550 = -t283 / 0.2e1;
t549 = t283 / 0.2e1;
t548 = -t284 / 0.2e1;
t547 = t284 / 0.2e1;
t545 = -t364 / 0.2e1;
t544 = t366 / 0.2e1;
t543 = -t367 / 0.2e1;
t539 = g(2) * t367;
t291 = t364 * t455 + qJDD(1);
t235 = (-Icges(6,5) * t342 - Icges(6,6) * t343) * t362;
t204 = t360 * t235;
t205 = t360 * t236;
t206 = t360 * t237;
t60 = -t204 * t364 + ((-t212 * t360 + t206) * t343 + (-t213 * t360 - t205) * t342) * t362;
t90 = -t211 * t364 + (-t212 * t342 + t213 * t343) * t362;
t538 = t90 * t291 + t60 * t295;
t102 = -t216 * t364 + (-t217 * t350 + t218 * t351) * t362;
t329 = -qJDD(4) * t364 + qJDD(1);
t254 = (-Icges(5,5) * t350 - Icges(5,6) * t351) * t362;
t231 = qJD(4) * t254;
t255 = (-Icges(5,2) * t351 - t522) * t362;
t232 = qJD(4) * t255;
t256 = (-Icges(5,1) * t350 - t521) * t362;
t233 = qJD(4) * t256;
t69 = -t231 * t364 + (-t232 * t350 + t233 * t351 + (-t217 * t351 - t218 * t350) * qJD(4)) * t362;
t537 = t102 * t329 + t69 * t339;
t536 = rSges(3,2) * t362;
t535 = pkin(4) * qJD(4);
t81 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t441;
t83 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t441;
t85 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t441;
t29 = -t364 * t81 + ((-t139 * t360 + t85) * t343 + (-t142 * t360 - t83) * t342) * t362;
t532 = t29 * t283;
t80 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t442;
t82 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t442;
t84 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t442;
t30 = -t364 * t80 + ((-t141 * t360 + t84) * t343 + (t143 * t360 - t82) * t342) * t362;
t531 = t30 * t284;
t446 = -t193 + t482;
t414 = t446 * qJD(1);
t65 = t414 + t588;
t530 = t366 * t65;
t58 = -t136 * t364 + (-t139 * t342 + t142 * t343) * t362;
t529 = t58 * t202;
t528 = t59 * t203;
t66 = -t155 * t364 + (-t158 * t350 + t161 * t351) * t362;
t527 = t66 * t275;
t526 = t67 * t276;
t525 = rSges(3,3) + qJ(2);
t515 = t155 * t367;
t514 = t156 * t366;
t513 = t357 * t362;
t269 = t292 * t508;
t429 = t472 * t362;
t117 = t269 - t293 + (-t299 + t541) * t367 - t366 * t429;
t145 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t510;
t502 = -t117 - t145;
t492 = -t197 - t214;
t488 = t217 - t256;
t487 = t218 + t255;
t486 = t292 * t440 + t299 * t470;
t485 = qJD(1) * t413 + qJDD(1) * t303;
t479 = t365 * t442 + t367 * t452;
t338 = rSges(3,2) * t509;
t264 = rSges(3,1) * t507 + rSges(3,3) * t366 - t338;
t478 = t305 + t264;
t476 = rSges(3,1) * t440 + rSges(3,3) * t470;
t475 = -t328 - t353;
t471 = qJD(1) * t362;
t466 = qJD(3) * t364;
t464 = qJD(4) * t366;
t463 = -m(4) - m(5) - m(6);
t460 = qJDD(2) * t366;
t458 = qJDD(3) * t364;
t457 = qJDD(3) * t367;
t53 = -t156 * t510 - t585;
t182 = -qJD(1) * t260 - qJD(4) * t259;
t183 = qJD(1) * t261 + qJD(4) * t258;
t101 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t441;
t164 = rSges(5,1) * t259 + rSges(5,2) * t258 + rSges(5,3) * t510;
t449 = rSges(4,3) * t441 + qJD(1) * t579;
t444 = t327 + t477;
t443 = t328 + t473;
t439 = t362 * t464;
t438 = t510 / 0.2e1;
t437 = -t509 / 0.2e1;
t436 = rSges(3,1) * t364 + pkin(1);
t382 = t260 * qJD(4);
t91 = pkin(4) * t382 + (-t596 + (t364 * t480 - t513) * t366) * qJD(1) + t479;
t92 = -t535 * t570 + (-qJD(1) * t429 - t351 * t535) * t367 - t481 + t486;
t9 = -t458 - t117 * t276 + t118 * t275 - t145 * t203 - t146 * t202 + t283 * t86 + t284 * t87 + (t366 * t91 + t367 * t92) * t465;
t435 = t9 * (t145 * t509 - t146 * t510);
t432 = t471 / 0.2e1;
t431 = -t465 / 0.2e1;
t430 = t465 / 0.2e1;
t427 = t174 * t295 - t242 * t283;
t425 = t327 - t468;
t421 = t366 * t432;
t420 = t367 * t432;
t419 = t366 * t431;
t418 = t366 * t430;
t417 = t367 * t431;
t416 = t367 * t430;
t247 = t260 * pkin(4);
t306 = rSges(2,1) * t367 - rSges(2,2) * t366;
t304 = rSges(2,1) * t366 + rSges(2,2) * t367;
t180 = qJD(1) * t258 + qJD(4) * t261;
t181 = qJD(1) * t259 + t382;
t411 = rSges(5,1) * t181 + rSges(5,2) * t180;
t94 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t442;
t95 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t441;
t96 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t442;
t97 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t441;
t98 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t442;
t99 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t441;
t407 = (t158 * t180 + t161 * t181 + t260 * t97 - t261 * t99 + (t155 * t470 - t367 * t95) * t362) * t366 - (t160 * t180 - t162 * t181 + t260 * t96 - t261 * t98 + (-t156 * t470 - t367 * t94) * t362) * t367;
t406 = (t158 * t182 + t161 * t183 + t258 * t97 + t259 * t99 + (t155 * t469 + t366 * t95) * t362) * t366 - (t160 * t182 - t162 * t183 + t258 * t96 + t259 * t98 + (-t156 * t469 + t366 * t94) * t362) * t367;
t32 = -t364 * t95 + (-t350 * t97 + t351 * t99 + (-t158 * t351 - t161 * t350) * qJD(4)) * t362;
t33 = -t364 * t94 + (-t350 * t96 + t351 * t98 + (-t160 * t351 + t162 * t350) * qJD(4)) * t362;
t405 = t32 * t366 - t33 * t367;
t262 = (-rSges(5,1) * t350 - rSges(5,2) * t351) * t362;
t234 = qJD(4) * t262;
t401 = qJD(1) * t444 + qJDD(1) * t285 + t366 * t434 + t485;
t381 = t401 - t460;
t34 = -t367 * t462 + t101 * t339 + t164 * t329 - t219 * t275 + (-t234 * t464 - t457) * t362 + t381 + t568;
t100 = rSges(5,3) * t442 + t411;
t370 = t345 + (-(-t516 + (-pkin(2) + t344) * t364) * t470 + t479 + t489) * qJD(1) + t446 * qJDD(1) + t564;
t35 = -t100 * t339 + t165 * t329 + t219 * t276 + (-t234 * t465 - qJDD(2)) * t367 + t370;
t404 = -t34 * t366 - t35 * t367;
t403 = t366 * t53 + t367 * t52;
t384 = (t192 + t483) * qJD(1) + t475;
t64 = t164 * t339 - t219 * t439 + t384;
t402 = -t366 * t64 - t367 * t65;
t400 = t100 * t366 + t101 * t367;
t398 = t164 * t367 - t165 * t366;
t397 = (Icges(5,5) * t258 - Icges(5,6) * t259) * t366 - (Icges(5,5) * t260 + Icges(5,6) * t261) * t367;
t396 = qJD(4) ^ 2 * t362 ^ 2 * t540 - qJDD(2);
t49 = -t136 * t509 - t569;
t386 = (t366 * t52 - t367 * t53) * t362;
t54 = -t155 * t509 - t501;
t55 = t156 * t509 + t399;
t385 = (t366 * t54 - t367 * t55) * t362;
t246 = t258 * pkin(4);
t383 = (Icges(6,5) * t238 - Icges(6,6) * t239) * t283 - (Icges(6,5) * t240 + Icges(6,6) * t241) * t284 + t235 * t295;
t13 = t139 * t150 + t142 * t151 + t240 * t83 - t241 * t85 + (t136 * t470 - t367 * t81) * t362;
t14 = t141 * t150 - t143 * t151 + t240 * t82 - t241 * t84 + (-t137 * t470 - t367 * t80) * t362;
t15 = t139 * t152 + t142 * t153 + t238 * t83 + t239 * t85 + (t136 * t469 + t366 * t81) * t362;
t16 = t141 * t152 - t143 * t153 + t238 * t82 + t239 * t84 + (-t137 * t469 + t366 * t80) * t362;
t22 = t283 * t49 - t576;
t37 = t150 * t212 + t151 * t213 + t205 * t240 - t206 * t241 + (-t204 * t367 + t211 * t470) * t362;
t38 = t152 * t212 + t153 * t213 + t205 * t238 + t206 * t239 + (t204 * t366 + t211 * t469) * t362;
t374 = (t15 * t283 - t16 * t284 + t202 * t47 + t203 * t48 + t291 * t71 + t295 * t38) * t438 - (-t383 * t364 + (-t342 * t371 - t343 * t556) * t362) * t295 / 0.2e1 + t22 * t421 + t21 * t420 + (-t364 * t71 + (t366 * t47 - t367 * t48) * t362) * t554 + (t13 * t283 - t14 * t284 + t202 * t49 + t203 * t50 - t291 * t378 + t295 * t37) * t437 + (t364 * t378 + (t366 * t49 - t367 * t50) * t362) * t553 + (-t364 * t38 + (t15 * t366 - t16 * t367 + (t366 * t48 + t367 * t47) * qJD(1)) * t362) * t549 + t291 * (-t364 * t90 + (t366 * t58 - t367 * t59) * t362) / 0.2e1 + (-t364 * t37 + (t13 * t366 - t14 * t367 + (t366 * t50 + t367 * t49) * qJD(1)) * t362) * t548 + (t528 + t529 - t531 + t532 + t538) * t545 + t295 * (-t364 * t60 + (t29 * t366 - t30 * t367 + (t366 * t59 + t367 * t58) * qJD(1)) * t362) / 0.2e1 + (t238 * t371 - t239 * t556 + t383 * t510) * t550 + (t240 * t371 + t241 * t556 - t383 * t509) * t547;
t207 = t360 * t242;
t200 = qJD(1) * t478 - t468;
t199 = -t353 + (-t563 + t303) * qJD(1);
t191 = rSges(5,1) * t260 + rSges(5,2) * t261;
t190 = rSges(5,1) * t258 - rSges(5,2) * t259;
t133 = t364 * t146;
t110 = qJD(1) * t445 + t425;
t107 = t478 * qJDD(1) + (qJD(1) * t563 - t279) * qJD(1) + t422;
t106 = -qJDD(1) * t563 - t460 + ((-rSges(3,2) * t471 - qJD(2)) * t367 + t476) * qJD(1) + t485;
t78 = t216 * t510 + t217 * t258 + t218 * t259;
t77 = t398 * t465 - t466;
t73 = t78 * t339;
t56 = -t362 * t457 + qJDD(1) * t194 + (t449 - t468) * qJD(1) + t381;
t46 = t182 * t217 + t183 * t218 + t232 * t258 + t233 * t259 + (t216 * t469 + t231 * t366) * t362;
t45 = t180 * t217 + t181 * t218 + t232 * t260 - t233 * t261 + (t216 * t470 - t231 * t367) * t362;
t44 = -t466 + t145 * t284 - t146 * t283 + (t117 * t367 + t118 * t366) * t465;
t43 = t414 + t593;
t42 = t117 * t339 + t145 * t295 - t197 * t439 - t214 * t283 + t384;
t36 = -t164 * t276 - t165 * t275 + t400 * t465 - t458;
t28 = qJD(4) * t385 - t595;
t27 = qJD(4) * t386 + t73;
t12 = -t118 * t329 + t146 * t291 + t197 * t276 + t203 * t214 - t207 * t284 - t295 * t86 - t339 * t91 + t367 * t396 + t370;
t11 = t329 * t117 + t291 * t145 - t275 * t197 - t202 * t214 - t283 * t207 + t295 * t87 + t339 * t92 + (-t459 - t462) * t367 + t396 * t366 + t401 + t568;
t1 = [(m(2) * (t304 ^ 2 + t306 ^ 2) + Icges(2,3) + (Icges(4,3) + Icges(3,2)) * t555 + ((Icges(4,1) * t363 ^ 2 + (-0.2e1 * Icges(4,4) * t363 + Icges(4,2) * t361) * t361 + Icges(3,1)) * t362 + 0.2e1 * (-Icges(4,5) * t363 + Icges(4,6) * t361 + Icges(3,4)) * t364) * t362) * qJDD(1) - t378 * t553 - t377 * t551 + t78 * t552 + t71 * t554 + t38 * t549 + (t73 + ((t55 + t592) * t366 + (t54 + (-t514 + t515) * t362 - t53 + t501) * t367) * t465) * t416 + t538 + t532 / 0.2e1 + (t32 + t46) * t418 + (t37 + t21) * t548 + (t33 + t45 + t27) * t417 + (-(qJD(1) * t264 - t200 + t301 - t468) * t199 + t200 * t473 + t199 * (t413 + t476) + ((rSges(3,3) * t200 - t199 * t536) * t367 + t200 * (-t436 + t536) * t366) * qJD(1) + (-g(2) + t107) * (t366 * t525 + t367 * t436 - t338) + (-g(3) + t106) * (-t367 * t525 + t356 - t578)) * m(3) + ((-g(3) + t11) * (-t357 * t510 + t269 + t356 + (-qJ(2) - t299) * t367 + t145) + (-g(2) + t12) * (t517 + (pkin(1) - t513) * t367 + t484 + t146) + (-t410 + t443 - t260 * t535 + (t596 + (-t292 * t364 - pkin(1) + (-rSges(6,3) + t357) * t362) * t366) * qJD(1)) * t43 + (t258 * t535 - t469 * t513 + t43 + t486 - t593 + t594 + t87) * t42) * m(6) + t537 + t21 * t547 - t531 / 0.2e1 + (t28 + ((t501 + t585) * t366 - t592 * t367 + ((t514 + t515) * t366 + t156 * t367 ^ 2) * t362 + t403) * t465 + t595) * t419 + (t65 * (-t411 + t443 + t479) + (-rSges(5,3) * t362 - t344 * t364 - pkin(1)) * t530 * qJD(1) + (-qJD(1) * t314 + t101 + t481 - t588 + t594 + t65) * t64 + (-g(2) + t35) * (t305 + t447 + t165) + (-g(3) + t34) * (t356 + (-qJ(2) - t541) * t367 + t415 + t164)) * m(5) + ((t412 + t443 + (-t542 - pkin(1) + (-rSges(4,3) - qJ(3)) * t362) * t470) * t110 + (-t539 + t583) * (pkin(1) + t408) + (-g(3) + t56) * t125 + (-qJD(1) * t448 + t110 - t425 + t444 + t449 + t597) * (qJD(1) * t125 + t475) + (-g(2) + t57) * (t448 + t517)) * m(4) - m(2) * (g(2) * t306 + g(3) * t304) + t528 / 0.2e1 + t529 / 0.2e1 + t526 / 0.2e1 + t527 / 0.2e1 + (t22 + (t48 + (t136 * t367 + t137 * t366) * t362 + t569 + t587) * t283 + t576) * t550; (-m(3) + t463) * (-g(3) * t366 - t539) + m(3) * (-t106 * t366 - t107 * t367) + m(4) * (-t366 * t56 - t583) + m(5) * t404 + m(6) * (-t11 * t366 - t12 * t367); t463 * (-g(1) * t364 + (g(2) * t366 - g(3) * t367) * t362) + 0.2e1 * (t9 * t545 + (t11 * t543 + t12 * t544) * t362) * m(6) + 0.2e1 * (t36 * t545 + (t34 * t543 + t35 * t544) * t362) * m(5) + 0.2e1 * (qJDD(3) * t555 / 0.2e1 + (t543 * t56 + t544 * t57) * t362) * m(4); (-t364 * t45 + ((t55 * t366 + t367 * t54) * qJD(1) + t407) * t362) * t417 + (-t364 * t46 + (qJD(1) * t403 + t406) * t362) * t418 + t27 * t420 + t28 * t421 + t374 + (t275 * t52 + t276 * t53 + t329 * t78 + t339 * t46 + t406 * t465) * t438 + t339 * (-t364 * t69 + ((t366 * t67 + t367 * t66) * qJD(1) + t405) * t362) / 0.2e1 - t339 * (-t364 * t254 * t339 + ((-t350 * t487 - t351 * t488) * t339 + ((-t350 * t379 - t351 * t557) * t362 - t397 * t364) * qJD(4)) * t362) / 0.2e1 + t329 * (-t102 * t364 + (t366 * t66 - t367 * t67) * t362) / 0.2e1 + (t405 * t465 + t526 + t527 + t537) * t545 + ((-t254 * t509 + t260 * t487 + t261 * t488) * t339 + (t379 * t260 + t261 * t557 - t397 * t509) * t465) * t416 + ((t254 * t510 + t258 * t487 - t259 * t488) * t339 + (t258 * t379 - t259 * t557 + t397 * t510) * t465) * t419 + (-t364 * t78 + t386) * t552 + (t364 * t377 + t385) * t551 + (t275 * t54 + t276 * t55 - t329 * t377 + t339 * t45 + t407 * t465) * t437 + (-g(2) * (t246 + t174) - g(3) * (t247 + t175) - t42 * (t246 * t339 + t427) + t435 - t12 * t133 + (t12 * t118 + t11 * t502 + t42 * (-t87 - t92)) * t364 + (-(t246 * t367 + t247 * t366) * t465 + t562) * t44 + (t247 * t339 + t364 * t91 + t561) * t43 + (-g(1) * (t409 - t540) + (t9 * t117 + t44 * t92 + t12 * t492 - t43 * t207 + (t118 * t44 + t42 * t492) * qJD(1)) * t367 + (t9 * t118 + t44 * t91 + t11 * t492 - t42 * t207 + (t197 * t43 + t44 * t502) * qJD(1)) * t366) * t362) * m(6) + (-(t190 * t64 - t191 * t65) * t339 - (t77 * (t190 * t367 + t191 * t366) + t402 * t262) * t465 + (t100 * t65 - t101 * t64 - t164 * t34 - t165 * t35) * t364 + (t36 * t398 + t77 * (-t164 * t470 - t165 * t469 + t400) + t402 * t234 + ((-t367 * t64 + t530) * qJD(1) + t404) * t219) * t362 - g(1) * t262 - g(2) * t190 - g(3) * t191) * m(5); t374 + (t435 + t12 * (-t214 * t509 - t133) + t11 * (-t145 * t364 - t214 * t510) - g(1) * t242 - g(2) * t174 - g(3) * t175 + (-t145 * t442 + t562) * t44 + (-t207 * t509 + t561) * t43 + (-t364 * t87 + (-t207 * t366 - t214 * t469) * t362 - t427) * t42) * m(6);];
tau = t1;
