% Calculate vector of inverse dynamics joint torques for
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:42
% DurationCPUTime: 17.60s
% Computational Cost: add. (6746->526), mult. (17896->646), div. (0->0), fcn. (17749->6), ass. (0->241)
t525 = Icges(5,4) - Icges(6,5);
t502 = Icges(5,1) + Icges(6,1);
t501 = Icges(6,4) + Icges(5,5);
t500 = Icges(5,2) + Icges(6,3);
t499 = Icges(5,6) - Icges(6,6);
t297 = sin(qJ(1));
t296 = cos(pkin(7));
t423 = sin(pkin(7));
t440 = sin(qJ(4));
t441 = cos(qJ(4));
t245 = -t296 * t440 + t423 * t441;
t298 = cos(qJ(1));
t226 = t245 * t298;
t244 = t296 * t441 + t423 * t440;
t227 = t244 * t298;
t225 = t244 * t297;
t224 = t245 * t297;
t541 = t525 * t224;
t480 = t502 * t225 + t501 * t298 + t541;
t540 = t525 * t225;
t481 = t500 * t224 + t499 * t298 + t540;
t527 = t481 * t226 + t480 * t227;
t542 = Icges(6,2) + Icges(5,3);
t532 = t499 * t224 + t501 * t225 + t542 * t298;
t473 = -t532 * t297 + t527;
t539 = t525 * t226;
t538 = t525 * t227;
t498 = -t500 * t226 + t499 * t297 - t538;
t518 = t502 * t227 - t501 * t297 + t539;
t535 = t224 * t481 + t225 * t480;
t496 = t499 * t226 + t501 * t227 - t542 * t297;
t531 = t525 * t245;
t530 = t525 * t244;
t479 = t500 * t244 - t531;
t529 = -t499 * t244 + t501 * t245;
t478 = t502 * t245 - t530;
t528 = -t498 * t226 + t518 * t227;
t526 = -t498 * t224 + t518 * t225;
t475 = t532 * t298 + t535;
t474 = t496 * t298 + t526;
t472 = t496 * t297 - t528;
t290 = t298 * qJ(2);
t258 = t297 * pkin(1) - t290;
t287 = qJD(2) * t297;
t389 = qJD(1) * t298;
t392 = qJ(2) * t389 + t287;
t524 = qJD(1) * t258 - t287 + t392;
t491 = rSges(6,1) + pkin(4);
t484 = -t479 * t224 + t478 * t225 + t529 * t298;
t483 = -t479 * t226 + t478 * t227 - t529 * t297;
t482 = rSges(6,3) + qJ(5);
t360 = t423 * qJ(3);
t439 = pkin(2) * t296;
t315 = t360 + t439;
t234 = t315 * t297;
t522 = qJD(1) * t234 + t524;
t471 = -t244 * t481 + t245 * t480;
t386 = qJD(4) * t298;
t460 = t245 * qJD(1);
t134 = -t244 * t386 - t297 * t460;
t390 = qJD(1) * t297;
t477 = qJD(4) * t245;
t135 = -t244 * t390 + t298 * t477;
t508 = -t134 * t500 - t135 * t525 + t389 * t499;
t136 = qJD(4) * t225 - t298 * t460;
t137 = qJD(1) * t227 + t297 * t477;
t521 = t136 * t500 - t137 * t525 + t390 * t499;
t504 = -t134 * t525 - t135 * t502 + t389 * t501;
t503 = -t136 * t525 + t137 * t502 - t390 * t501;
t232 = t244 * qJD(4);
t517 = t232 * t525 + t477 * t500;
t516 = -t232 * t501 - t477 * t499;
t515 = -t232 * t502 - t477 * t525;
t513 = t472 * t297 + t298 * t473;
t512 = -t474 * t297 + t475 * t298;
t415 = t296 * t297;
t438 = t298 * pkin(6);
t248 = pkin(3) * t415 + t438;
t511 = qJD(1) * t248 + t522;
t510 = (t134 * t499 + t135 * t501 - t389 * t542) * t297 - (-t136 * t499 + t137 * t501 - t390 * t542) * t298;
t509 = t483 * qJD(1);
t120 = t225 * rSges(5,1) + t224 * rSges(5,2) + t298 * rSges(5,3);
t364 = t290 - t438;
t494 = -t120 + t364;
t401 = t482 * t244 + t245 * t491;
t493 = t484 * qJD(1);
t387 = qJD(4) * t297;
t383 = qJD(1) * qJD(4);
t250 = -qJDD(4) * t297 - t298 * t383;
t414 = t296 * t298;
t277 = pkin(3) * t414;
t249 = -pkin(6) * t297 + t277;
t235 = pkin(2) * t414 + t298 * t360;
t358 = qJD(3) * t423;
t271 = t298 * t358;
t260 = t298 * pkin(1) + t297 * qJ(2);
t384 = qJD(1) * qJD(2);
t318 = -qJDD(2) * t298 + qJD(1) * (-pkin(1) * t390 + t392) + qJDD(1) * t260 + t297 * t384;
t341 = qJD(1) * t358;
t356 = qJDD(3) * t423;
t304 = qJD(1) * (-t315 * t390 + t271) + qJDD(1) * t235 + t298 * t341 + t297 * t356 + t318;
t450 = qJD(1) ^ 2;
t301 = qJDD(1) * t249 - t450 * t248 + t304;
t385 = qJD(5) * t244;
t424 = -t232 * t491 + t477 * t482 + t385;
t357 = qJD(4) * t424;
t465 = -t482 * t226 + t491 * t227;
t411 = -rSges(6,2) * t297 + t465;
t195 = qJD(5) * t226;
t466 = -t482 * t134 + t491 * t135 - t195;
t435 = -rSges(6,2) * t389 + t466;
t3 = qJD(1) * t435 + qJD(5) * t136 + qJDD(1) * t411 - qJDD(5) * t224 - t250 * t401 + t297 * t357 + t301;
t492 = -g(2) + t3;
t490 = t512 * qJD(4) + t493;
t489 = t513 * qJD(4) + t509;
t488 = -t480 * t232 + t521 * t244 + t503 * t245 - t481 * t477;
t487 = t518 * t232 - t508 * t244 + t504 * t245 - t498 * t477;
t486 = -t479 * t134 + t478 * t135 - t517 * t226 + t515 * t227 - t516 * t297 - t389 * t529;
t485 = t479 * t136 + t478 * t137 - t517 * t224 + t515 * t225 + t516 * t298 - t390 * t529;
t470 = -t498 * t244 - t518 * t245;
t413 = t298 * rSges(6,2) - t482 * t224 + t491 * t225;
t476 = t472 - t475;
t194 = qJD(5) * t224;
t467 = -t482 * t136 - t491 * t137 + t194;
t464 = -t501 * t244 - t499 * t245;
t270 = t297 * t358;
t288 = qJD(2) * t298;
t396 = t270 - t288;
t463 = t194 - t396;
t462 = (t500 * t225 - t480 - t541) * t298 + (-t500 * t227 + t518 + t539) * t297;
t461 = (t502 * t224 - t481 - t540) * t298 + (-t502 * t226 - t498 + t538) * t297;
t459 = (t501 * t224 - t499 * t225) * t298 + (-t501 * t226 + t499 * t227) * t297;
t458 = (-t498 * t136 - t137 * t518 + t508 * t224 + t504 * t225 + t496 * t390) * t297 + (-t481 * t136 + t480 * t137 - t224 * t521 + t503 * t225 - t390 * t532 - t510) * t298;
t457 = (t481 * t134 + t480 * t135 - t226 * t521 + t503 * t227 - t389 * t532) * t298 + (t498 * t134 - t135 * t518 + t508 * t226 + t504 * t227 + t496 * t389 + t510) * t297;
t292 = t297 * rSges(4,2);
t362 = t298 * t423;
t220 = rSges(4,1) * t414 + rSges(4,3) * t362 + t292;
t353 = t235 + t260;
t133 = t353 + t220;
t456 = -t502 * t244 + t479 - t531;
t455 = t500 * t245 - t478 + t530;
t453 = t401 * t386 - t195;
t451 = t364 - t413;
t449 = t250 / 0.2e1;
t251 = qJDD(4) * t298 - t297 * t383;
t448 = t251 / 0.2e1;
t446 = t297 / 0.2e1;
t445 = -t298 / 0.2e1;
t443 = -rSges(6,2) - pkin(6);
t442 = -rSges(5,3) - pkin(6);
t434 = rSges(6,2) * t390 + t467;
t429 = rSges(3,1) * t296;
t428 = rSges(4,1) * t296;
t291 = t297 * rSges(3,3);
t410 = t135 * rSges(5,1) + t134 * rSges(5,2);
t408 = -t224 * t491 - t225 * t482;
t407 = t226 * t491 + t227 * t482;
t402 = -t244 * t491 + t245 * t482;
t399 = t227 * rSges(5,1) + t226 * rSges(5,2);
t221 = rSges(3,1) * t414 - rSges(3,2) * t362 + t291;
t191 = t221 + t260;
t398 = -t234 - t258;
t359 = qJD(1) * t423;
t346 = t297 * t359;
t397 = rSges(3,2) * t346 + rSges(3,3) * t389;
t395 = t271 + t287;
t363 = t297 * t423;
t394 = rSges(3,2) * t363 + t298 * rSges(3,3);
t393 = qJDD(2) * t297 + t298 * t384;
t388 = qJD(3) * t296;
t382 = qJDD(3) * t296;
t380 = rSges(3,1) * t415;
t379 = -t248 + t398;
t378 = t249 + t353;
t268 = t298 * t356;
t377 = t268 + t393;
t285 = pkin(6) * t390;
t375 = t285 - t396;
t372 = -pkin(1) - t429;
t371 = t423 * rSges(4,3);
t368 = -t387 / 0.2e1;
t367 = t387 / 0.2e1;
t366 = -t386 / 0.2e1;
t365 = t386 / 0.2e1;
t330 = t379 - t413;
t30 = qJD(1) * t330 + t395 + t453;
t361 = t30 * t401;
t355 = -t120 + t379;
t340 = -pkin(1) - t360;
t261 = rSges(2,1) * t298 - rSges(2,2) * t297;
t259 = rSges(2,1) * t297 + rSges(2,2) * t298;
t339 = t137 * rSges(5,1) - t136 * rSges(5,2);
t189 = rSges(5,1) * t245 - rSges(5,2) * t244;
t171 = t189 * t386;
t46 = qJD(1) * t355 + t171 + t395;
t124 = -rSges(5,3) * t297 + t399;
t47 = t189 * t387 + (t124 + t378) * qJD(1) + t396;
t332 = t297 * t47 + t298 * t46;
t75 = -rSges(5,3) * t389 + t410;
t77 = -rSges(5,3) * t390 + t339;
t331 = -t297 * t77 - t298 * t75;
t329 = t378 + t411;
t327 = t277 + t353;
t322 = -t120 * t297 - t124 * t298;
t317 = -rSges(3,2) * t423 + t429;
t316 = t371 + t428;
t230 = qJD(1) * t260 - t288;
t314 = -qJD(1) * t230 - qJDD(1) * t258 + t393;
t193 = t315 * t389 + t270;
t313 = -qJD(1) * t277 - t193 - t230 - t270 + t285;
t307 = -t371 + t340;
t306 = -pkin(3) * t296 - pkin(1) - t315;
t305 = t297 * ((-pkin(2) - pkin(3)) * t296 + t340);
t302 = g(1) * t305;
t300 = (-rSges(4,1) - pkin(2)) * t296 + t307;
t294 = t298 * rSges(4,2);
t284 = rSges(4,2) * t389;
t222 = g(1) * t362 + g(2) * t363 - g(3) * t296;
t219 = t380 - t394;
t218 = t297 * t316 - t294;
t186 = -rSges(5,1) * t244 - rSges(5,2) * t245;
t168 = -rSges(5,1) * t232 - rSges(5,2) * t477;
t160 = qJD(1) * t191 - t288;
t159 = t287 + (-t219 - t258) * qJD(1);
t157 = rSges(5,1) * t226 - rSges(5,2) * t227;
t152 = rSges(5,1) * t224 - rSges(5,2) * t225;
t97 = qJD(1) * t133 + t396;
t96 = (-t218 + t398) * qJD(1) + t395;
t79 = qJDD(1) * t221 + qJD(1) * (-qJD(1) * t380 + t397) + t318;
t78 = -qJDD(1) * t219 - t450 * (t298 * t317 + t291) + t314;
t61 = qJD(4) * t322 - t388;
t42 = qJDD(1) * t220 + qJD(1) * (-t316 * t390 + t284) + t304;
t41 = -qJD(1) * t193 + t268 - t297 * t341 - t450 * (t298 * t316 + t292) + t314 + (-t234 - t218) * qJDD(1);
t32 = -t388 + t385 + (-t297 * t413 - t298 * t411) * qJD(4);
t31 = t329 * qJD(1) + t387 * t401 - t463;
t26 = qJD(4) * t331 + t120 * t250 - t124 * t251 - t382;
t17 = qJD(1) * t75 + qJDD(1) * t124 + t168 * t387 - t189 * t250 + t301;
t16 = t168 * t386 + t251 * t189 + t355 * qJDD(1) + (-t77 + t313) * qJD(1) + t377;
t2 = -qJD(5) * t134 - qJDD(5) * t226 + t401 * t251 + t298 * t357 + t330 * qJDD(1) + (t313 + t434) * qJD(1) + t377;
t1 = qJD(5) * t477 - t382 + qJDD(5) * t244 - t411 * t251 + t413 * t250 + (t297 * t434 - t298 * t435) * qJD(4);
t4 = [-m(2) * (-g(1) * t259 + g(2) * t261) + ((t527 * t298 + (t476 + t535) * t297) * qJD(4) + t509) * t366 + (-g(1) * t451 - t302 + t361 * t387 + (t305 + t451) * t2 + t492 * (t297 * t443 + t327 + t465) + (-t453 + t466 + t511) * t31 + (-t463 + t375 + t467) * t30) * m(6) + (t46 * (-t339 + t375) - g(1) * t494 - t302 + (t305 + t494) * t16 + (t17 - g(2)) * (t297 * t442 + t327 + t399) + (t410 - t171 + t46 + t511) * t47) * m(5) + (-t96 * t396 + (t42 - g(2)) * t133 + (t41 - g(1)) * (t297 * t300 + t290 + t294) + (t284 + t96 + t522) * t97) * m(4) + (t159 * t288 + (t79 - g(2)) * t191 + (t78 - g(1)) * (t372 * t297 + t290 + t394) + (t397 + t159 + t524) * t160) * m(3) + (-t470 + t483) * t449 + (t471 + t484) * t448 + (t486 - t487) * t368 + (((t476 + t528) * t298 + t526 * t297) * qJD(4) + t490 - t493) * t367 + (t485 + t488 + t489) * t365 + (m(2) * (t259 ^ 2 + t261 ^ 2) + Icges(2,3) + (Icges(4,1) + Icges(3,1)) * t423 ^ 2 + t478 * t245 + t479 * t244 + ((Icges(4,3) + Icges(3,2)) * t296 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t423) * t296) * qJDD(1) + (t515 * t245 + t517 * t244 + t479 * t477 - t478 * t232 + ((t30 * t306 + t31 * t443) * t298 + (t30 * (rSges(6,2) - qJ(2)) + t31 * t306) * t297 + t30 * t329 + t31 * t413) * m(6) + ((t306 * t46 + t442 * t47) * t298 + (t46 * (rSges(5,3) - qJ(2)) + t47 * t306) * t297 + t120 * t47) * m(5) + (t96 * t300 * t298 + (t96 * (-rSges(4,2) - qJ(2)) + t97 * (t307 - t428 - t439)) * t297 + t218 * t97) * m(4) + (t159 * (-pkin(1) - t317) * t298 + (t159 * (-rSges(3,3) - qJ(2)) + t160 * t372) * t297 + t219 * t160) * m(3)) * qJD(1); (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t297 - g(2) * t298) + 0.2e1 * (t2 * t446 + t3 * t445) * m(6) + 0.2e1 * (t16 * t446 + t17 * t445) * m(5) + 0.2e1 * (t41 * t446 + t42 * t445) * m(4) + 0.2e1 * (t445 * t79 + t446 * t78) * m(3); (t2 * t362 - t1 * t296 + (t3 * t423 - t30 * t359) * t297 + t30 * t346 - t222) * m(6) + (t16 * t362 - t26 * t296 + (t17 * t423 - t359 * t46) * t297 + t346 * t46 - t222) * m(5) + (t41 * t362 + qJDD(3) * t296 ^ 2 + (-t359 * t96 + t42 * t423) * t297 + t346 * t96 - t222) * m(4); t513 * t449 + t512 * t448 - (qJD(1) * t486 + t457 * qJD(4) + qJDD(1) * t483 - t472 * t250 + t473 * t251) * t297 / 0.2e1 + (qJD(1) * t485 + t458 * qJD(4) + qJDD(1) * t484 + t474 * t250 + t475 * t251) * t298 / 0.2e1 - ((t244 * t462 + t245 * t461) * qJD(4) + (t455 * t244 + t456 * t245) * qJD(1)) * qJD(1) / 0.2e1 + (t488 * t298 + t487 * t297 + (-t471 * t297 + t470 * t298) * qJD(1)) * qJD(1) / 0.2e1 + (t470 * t297 + t471 * t298) * qJDD(1) / 0.2e1 - t490 * t390 / 0.2e1 - t489 * t389 / 0.2e1 + ((-t473 * t297 + t298 * t472) * qJD(1) + t457) * t368 + ((-t226 * t462 + t227 * t461 - t297 * t459) * qJD(4) + (-t226 * t455 + t227 * t456 - t297 * t464) * qJD(1)) * t367 + ((-t224 * t462 + t225 * t461 + t298 * t459) * qJD(4) + (-t224 * t455 + t225 * t456 + t298 * t464) * qJD(1)) * t366 + ((-t297 * t475 - t474 * t298) * qJD(1) + t458) * t365 + ((t2 * t401 + t30 * t424 - t1 * t411 - t32 * t435 + (t31 * t401 - t32 * t413) * qJD(1)) * t298 + (t3 * t401 + t31 * t424 - t1 * t413 + t32 * t434 + (t32 * t411 - t361) * qJD(1)) * t297 - g(1) * t407 + g(2) * t408 - g(3) * t402 - (t225 * t31 + t227 * t30 + t245 * t32) * qJD(5) - (t30 * t408 + t31 * t407) * qJD(1) - ((t30 * t402 - t32 * t407) * t298 + (t31 * t402 + t32 * t408) * t297) * qJD(4)) * m(6) + (-(-t152 * t46 + t157 * t47) * qJD(1) - (t61 * (-t152 * t297 - t157 * t298) + t332 * t186) * qJD(4) + t26 * t322 + t61 * ((-t120 * t298 + t124 * t297) * qJD(1) + t331) + t332 * t168 + (t16 * t298 + t17 * t297 + (-t297 * t46 + t298 * t47) * qJD(1)) * t189 - g(1) * t157 - g(2) * t152 - g(3) * t186) * m(5); (-t134 * t30 + t136 * t31 - (-qJD(1) * t31 - g(1) + t2) * t226 - (qJD(1) * t30 + t492) * t224 + (t477 - (t224 * t297 + t226 * t298) * qJD(4)) * t32 + (t1 - (t297 * t31 + t298 * t30) * qJD(4) - g(3)) * t244) * m(6);];
tau = t4;
