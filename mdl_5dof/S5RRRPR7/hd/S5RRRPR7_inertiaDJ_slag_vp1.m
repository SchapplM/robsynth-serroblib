% Calculate time derivative of joint inertia matrix for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:07
% EndTime: 2019-12-31 21:16:33
% DurationCPUTime: 13.67s
% Computational Cost: add. (29341->859), mult. (31613->1170), div. (0->0), fcn. (29615->10), ass. (0->427)
t339 = qJ(2) + qJ(3);
t331 = cos(t339);
t341 = cos(pkin(9));
t344 = sin(qJ(1));
t485 = t341 * t344;
t340 = sin(pkin(9));
t346 = cos(qJ(1));
t486 = t340 * t346;
t287 = t331 * t485 - t486;
t484 = t341 * t346;
t487 = t340 * t344;
t373 = t331 * t487 + t484;
t330 = sin(t339);
t495 = t330 * t344;
t195 = Icges(5,4) * t287 - Icges(5,2) * t373 + Icges(5,6) * t495;
t197 = Icges(5,1) * t287 - Icges(5,4) * t373 + Icges(5,5) * t495;
t397 = Icges(4,5) * t331 - Icges(4,6) * t330;
t252 = -Icges(4,3) * t346 + t397 * t344;
t288 = -t331 * t486 + t485;
t289 = t331 * t484 + t487;
t193 = Icges(5,5) * t287 - Icges(5,6) * t373 + Icges(5,3) * t495;
t494 = t330 * t346;
t455 = t193 * t494;
t515 = Icges(4,4) * t331;
t401 = -Icges(4,2) * t330 + t515;
t254 = -Icges(4,6) * t346 + t401 * t344;
t516 = Icges(4,4) * t330;
t406 = Icges(4,1) * t331 - t516;
t256 = -Icges(4,5) * t346 + t406 * t344;
t387 = t254 * t330 - t256 * t331;
t555 = t387 * t346;
t562 = -t195 * t288 - t197 * t289 - t252 * t344 - t455 + t555;
t196 = Icges(5,4) * t289 + Icges(5,2) * t288 + Icges(5,6) * t494;
t198 = Icges(5,1) * t289 + Icges(5,4) * t288 + Icges(5,5) * t494;
t253 = Icges(4,3) * t344 + t397 * t346;
t255 = Icges(4,6) * t344 + t401 * t346;
t257 = Icges(4,5) * t344 + t406 * t346;
t386 = t255 * t330 - t257 * t331;
t194 = Icges(5,5) * t289 + Icges(5,6) * t288 + Icges(5,3) * t494;
t453 = t194 * t494;
t561 = t196 * t288 + t198 * t289 + t253 * t344 - t386 * t346 + t453;
t347 = -pkin(7) - pkin(6);
t343 = sin(qJ(2));
t464 = qJD(2) * t343;
t458 = pkin(2) * t464;
t560 = qJD(1) * t347 + t458;
t345 = cos(qJ(2));
t309 = rSges(3,1) * t343 + rSges(3,2) * t345;
t370 = qJD(2) * t309;
t559 = t344 * t370;
t517 = Icges(3,4) * t345;
t403 = -Icges(3,2) * t343 + t517;
t276 = Icges(3,6) * t344 + t403 * t346;
t518 = Icges(3,4) * t343;
t408 = Icges(3,1) * t345 - t518;
t278 = Icges(3,5) * t344 + t408 * t346;
t384 = t276 * t343 - t278 * t345;
t558 = t384 * t344;
t275 = -Icges(3,6) * t346 + t403 * t344;
t277 = -Icges(3,5) * t346 + t408 * t344;
t385 = t275 * t343 - t277 * t345;
t557 = t385 * t346;
t556 = t386 * t344;
t342 = -pkin(8) - qJ(4);
t482 = qJ(4) + t342;
t554 = t482 * t331;
t324 = pkin(4) * t341 + pkin(3);
t528 = pkin(3) - t324;
t553 = t528 * t330;
t552 = t528 * t331;
t325 = pkin(2) * t345 + pkin(1);
t529 = pkin(1) - t325;
t551 = t529 * t344;
t335 = pkin(9) + qJ(5);
t328 = sin(t335);
t329 = cos(t335);
t336 = qJD(2) + qJD(3);
t413 = rSges(6,1) * t329 - rSges(6,2) * t328;
t493 = t331 * t336;
t150 = t413 * t493 + (rSges(6,3) * t336 + (-rSges(6,1) * t328 - rSges(6,2) * t329) * qJD(5)) * t330;
t356 = -t482 * t330 - t552;
t550 = -t356 * t336 - t150;
t491 = t331 * t344;
t265 = -t328 * t491 - t329 * t346;
t266 = -t328 * t346 + t329 * t491;
t414 = -t266 * rSges(6,1) - t265 * rSges(6,2);
t170 = rSges(6,3) * t495 - t414;
t490 = t331 * t346;
t267 = -t328 * t490 + t329 * t344;
t268 = t328 * t344 + t329 * t490;
t171 = t268 * rSges(6,1) + t267 * rSges(6,2) + rSges(6,3) * t494;
t549 = -t344 * t170 - t346 * t171;
t229 = -t553 + t554;
t234 = -rSges(6,3) * t331 + t413 * t330;
t548 = -t229 - t234;
t547 = qJD(1) * t252;
t398 = Icges(3,5) * t345 - Icges(3,6) * t343;
t273 = -Icges(3,3) * t346 + t398 * t344;
t423 = qJD(1) * t331 - qJD(5);
t488 = t336 * t346;
t450 = t330 * t488;
t546 = t423 * t344 + t450;
t224 = t373 * qJD(1) + t340 * t450;
t225 = -t287 * qJD(1) - t341 * t450;
t467 = qJD(1) * t344;
t437 = t330 * t467;
t449 = t331 * t488;
t360 = -t437 + t449;
t120 = Icges(5,5) * t225 + Icges(5,6) * t224 + t360 * Icges(5,3);
t489 = t336 * t344;
t451 = t330 * t489;
t226 = t288 * qJD(1) + t340 * t451;
t227 = t289 * qJD(1) - t341 * t451;
t466 = qJD(1) * t346;
t361 = t330 * t466 + t331 * t489;
t121 = Icges(5,5) * t227 + Icges(5,6) * t226 + t361 * Icges(5,3);
t545 = -(t121 * t330 + t193 * t493) * t346 + (t120 * t330 + t194 * t493) * t344;
t544 = t423 * t346 - t451;
t292 = Icges(4,2) * t331 + t516;
t293 = Icges(4,1) * t330 + t515;
t383 = t292 * t330 - t293 * t331;
t543 = t383 * qJD(1) + t397 * t336;
t542 = 2 * m(3);
t541 = 2 * m(4);
t540 = 2 * m(5);
t539 = 2 * m(6);
t337 = t344 ^ 2;
t338 = t346 ^ 2;
t538 = m(5) / 0.2e1;
t537 = m(6) / 0.2e1;
t536 = t344 / 0.2e1;
t535 = -t346 / 0.2e1;
t534 = m(3) * t309;
t295 = rSges(4,1) * t330 + rSges(4,2) * t331;
t533 = m(4) * t295;
t532 = pkin(2) * t343;
t531 = pkin(3) * t331;
t530 = t344 * pkin(6);
t334 = t346 * pkin(6);
t527 = -pkin(6) - t347;
t526 = rSges(3,1) * t345;
t525 = rSges(4,1) * t331;
t524 = rSges(3,2) * t343;
t523 = rSges(3,3) * t346;
t164 = Icges(6,5) * t266 + Icges(6,6) * t265 + Icges(6,3) * t495;
t166 = Icges(6,4) * t266 + Icges(6,2) * t265 + Icges(6,6) * t495;
t168 = Icges(6,1) * t266 + Icges(6,4) * t265 + Icges(6,5) * t495;
t394 = -t166 * t328 + t168 * t329;
t424 = -qJD(5) * t331 + qJD(1);
t381 = t344 * t424;
t153 = -t544 * t328 + t329 * t381;
t154 = t328 * t381 + t544 * t329;
t90 = Icges(6,5) * t154 + Icges(6,6) * t153 + t361 * Icges(6,3);
t92 = Icges(6,4) * t154 + Icges(6,2) * t153 + t361 * Icges(6,6);
t94 = Icges(6,1) * t154 + Icges(6,4) * t153 + t361 * Icges(6,5);
t23 = (t394 * t336 - t90) * t331 + (t164 * t336 - t328 * t92 + t329 * t94 + (-t166 * t329 - t168 * t328) * qJD(5)) * t330;
t522 = t23 * t346;
t165 = Icges(6,5) * t268 + Icges(6,6) * t267 + Icges(6,3) * t494;
t167 = Icges(6,4) * t268 + Icges(6,2) * t267 + Icges(6,6) * t494;
t169 = Icges(6,1) * t268 + Icges(6,4) * t267 + Icges(6,5) * t494;
t393 = -t167 * t328 + t169 * t329;
t380 = t346 * t424;
t151 = t546 * t328 + t329 * t380;
t152 = t328 * t380 - t546 * t329;
t89 = Icges(6,5) * t152 + Icges(6,6) * t151 + t360 * Icges(6,3);
t91 = Icges(6,4) * t152 + Icges(6,2) * t151 + t360 * Icges(6,6);
t93 = Icges(6,1) * t152 + Icges(6,4) * t151 + t360 * Icges(6,5);
t24 = (t393 * t336 - t89) * t331 + (t165 * t336 - t328 * t91 + t329 * t93 + (-t167 * t329 - t169 * t328) * qJD(5)) * t330;
t521 = t24 * t344;
t333 = t344 * rSges(3,3);
t332 = t344 * rSges(4,3);
t520 = -rSges(5,3) - qJ(4);
t519 = -rSges(6,3) + t342;
t514 = Icges(6,4) * t328;
t513 = Icges(6,4) * t329;
t506 = qJ(4) * t330;
t404 = Icges(6,1) * t329 - t514;
t233 = -Icges(6,5) * t331 + t404 * t330;
t505 = t233 * t329;
t419 = -rSges(4,2) * t330 + t525;
t272 = t419 * t336;
t504 = t272 * t344;
t503 = t275 * t345;
t502 = t276 * t345;
t501 = t277 * t343;
t500 = t278 * t343;
t499 = t292 * t336;
t498 = t293 * t336;
t497 = t324 * t330;
t496 = t330 * t336;
t492 = t331 * t342;
t483 = t346 * t347;
t200 = t289 * rSges(5,1) + t288 * rSges(5,2) + rSges(5,3) * t494;
t316 = pkin(3) * t490;
t284 = qJ(4) * t494 + t316;
t481 = -t200 - t284;
t416 = rSges(5,1) * t341 - rSges(5,2) * t340;
t223 = (rSges(5,3) * t330 + t416 * t331) * t336;
t241 = qJ(4) * t496 + (pkin(3) * t336 - qJD(4)) * t331;
t480 = -t223 - t241;
t247 = -rSges(5,3) * t331 + t416 * t330;
t294 = pkin(3) * t330 - qJ(4) * t331;
t285 = t294 * t467;
t479 = t247 * t467 + t285;
t250 = t334 + t483 - t551;
t312 = t346 * t325;
t251 = -t346 * pkin(1) + t527 * t344 + t312;
t478 = t344 * t250 + t346 * t251;
t261 = -t346 * rSges(4,3) + t419 * t344;
t262 = rSges(4,1) * t490 - rSges(4,2) * t494 + t332;
t178 = t344 * t261 + t346 * t262;
t477 = -t247 - t294;
t283 = (t506 + t531) * t344;
t476 = t344 * t283 + t346 * t284;
t299 = qJ(4) * t449;
t462 = qJD(4) * t330;
t310 = t346 * t462;
t475 = t299 + t310;
t460 = pkin(4) * t486;
t474 = qJD(1) * t460 + t342 * t437;
t473 = rSges(4,2) * t437 + rSges(4,3) * t466;
t472 = t560 * t344;
t471 = t346 * t526 + t333;
t470 = t337 + t338;
t469 = qJD(1) * t253;
t274 = Icges(3,3) * t344 + t398 * t346;
t468 = qJD(1) * t274;
t463 = qJD(2) * t345;
t461 = qJD(5) * t330;
t321 = pkin(4) * t487;
t459 = t346 * t524;
t457 = pkin(2) * t463;
t456 = t193 * t495;
t454 = t194 * t495;
t399 = -Icges(6,2) * t328 + t513;
t232 = -Icges(6,6) * t331 + t399 * t330;
t452 = t232 * t493;
t395 = Icges(6,5) * t329 - Icges(6,6) * t328;
t231 = -Icges(6,3) * t331 + t395 * t330;
t100 = t231 * t495 + t232 * t265 + t233 * t266;
t72 = -t164 * t331 + t394 * t330;
t448 = t100 / 0.2e1 + t72 / 0.2e1;
t101 = t231 * t494 + t232 * t267 + t233 * t268;
t73 = -t165 * t331 + t393 * t330;
t447 = t101 / 0.2e1 + t73 / 0.2e1;
t446 = t152 * rSges(6,1) + t151 * rSges(6,2) + rSges(6,3) * t449;
t445 = -t241 + t550;
t303 = pkin(3) * t451;
t362 = -t331 * t467 - t450;
t444 = t344 * (t361 * qJ(4) + qJD(1) * t316 + t344 * t462 - t303) + t346 * (t362 * pkin(3) - qJ(4) * t437 + t475) + t283 * t466;
t377 = t324 * t490 - t342 * t494 + t321;
t192 = t377 - t284;
t443 = -t171 - t192 - t284;
t372 = t295 * t336;
t442 = t344 * (-t344 * t372 + (t419 * t346 + t332) * qJD(1)) + t346 * (t362 * rSges(4,1) - rSges(4,2) * t449 + t473) + t261 * t466;
t441 = t344 * ((-t529 * t346 - t530) * qJD(1) - t472) + t346 * (-t346 * t458 + (t527 * t346 + t551) * qJD(1)) + t250 * t466;
t228 = t234 * t467;
t440 = t229 * t467 + t228 + t285;
t439 = t225 * rSges(5,1) + t224 * rSges(5,2) + rSges(5,3) * t449;
t438 = -t294 + t548;
t436 = t343 * t467;
t434 = t467 / 0.2e1;
t433 = t466 / 0.2e1;
t432 = -t294 - t532;
t431 = -t295 - t532;
t203 = t477 * t346;
t183 = -t254 * qJD(1) - t346 * t499;
t430 = t257 * t336 + t183;
t184 = t255 * qJD(1) - t344 * t499;
t429 = t256 * t336 + t184;
t185 = -t256 * qJD(1) - t346 * t498;
t428 = -t255 * t336 + t185;
t186 = t257 * qJD(1) - t344 * t498;
t427 = t254 * t336 - t186;
t426 = -t344 * t347 + t312;
t425 = -t324 * t331 - t325;
t417 = -t287 * rSges(5,1) + rSges(5,2) * t373;
t199 = rSges(5,3) * t495 - t417;
t99 = t344 * t199 + t346 * t200 + t476;
t422 = -t247 + t432;
t421 = -t241 - t457;
t133 = t438 * t346;
t420 = -t524 + t526;
t418 = rSges(5,1) * t227 + rSges(5,2) * t226;
t415 = rSges(6,1) * t154 + rSges(6,2) * t153;
t63 = t164 * t495 + t166 * t265 + t168 * t266;
t64 = t165 * t495 + t167 * t265 + t169 * t266;
t44 = t344 * t64 - t346 * t63;
t412 = t344 * t63 + t346 * t64;
t65 = t164 * t494 + t166 * t267 + t168 * t268;
t66 = t165 * t494 + t167 * t267 + t169 * t268;
t45 = t344 * t66 - t346 * t65;
t411 = t344 * t65 + t346 * t66;
t410 = t344 * t73 - t72 * t346;
t409 = t344 * t72 + t346 * t73;
t407 = Icges(3,1) * t343 + t517;
t405 = Icges(5,1) * t341 - Icges(5,4) * t340;
t402 = Icges(3,2) * t345 + t518;
t400 = Icges(5,4) * t341 - Icges(5,2) * t340;
t291 = Icges(4,5) * t330 + Icges(4,6) * t331;
t396 = Icges(5,5) * t341 - Icges(5,6) * t340;
t392 = t170 * t346 - t171 * t344;
t391 = -t195 * t340 + t197 * t341;
t390 = -t196 * t340 + t198 * t341;
t382 = t432 + t548;
t379 = -t223 + t421;
t378 = -pkin(1) - t420;
t180 = t422 * t346;
t376 = -t325 - t419;
t375 = t344 * (t361 * rSges(5,3) + t418) + t346 * (-rSges(5,3) * t437 + t439) + t199 * t466 + t444;
t191 = t356 * t344 - t460;
t56 = t344 * t191 + t346 * t192 + t476 - t549;
t371 = t421 + t550;
t367 = t336 * t291;
t366 = qJD(2) * t407;
t365 = qJD(2) * t402;
t364 = qJD(2) * (-Icges(3,5) * t343 - Icges(3,6) * t345);
t127 = t382 * t346;
t363 = t520 * t330 - t325 - t531;
t359 = t519 * t330 + t425;
t122 = Icges(5,4) * t225 + Icges(5,2) * t224 + t360 * Icges(5,6);
t123 = Icges(5,4) * t227 + Icges(5,2) * t226 + t361 * Icges(5,6);
t124 = Icges(5,1) * t225 + Icges(5,4) * t224 + t360 * Icges(5,5);
t125 = Icges(5,1) * t227 + Icges(5,4) * t226 + t361 * Icges(5,5);
t128 = -t252 * t346 - t387 * t344;
t129 = -t253 * t346 - t556;
t181 = -t346 * t367 - t547;
t182 = -t344 * t367 + t469;
t74 = -t195 * t373 + t197 * t287 + t456;
t75 = -t196 * t373 + t198 * t287 + t454;
t15 = t151 * t166 + t152 * t168 + t360 * t164 + t267 * t92 + t268 * t94 + t90 * t494;
t16 = t151 * t167 + t152 * t169 + t360 * t165 + t267 * t91 + t268 * t93 + t89 * t494;
t8 = t411 * qJD(1) - t15 * t346 + t16 * t344;
t358 = t44 * t467 + t45 * t466 + ((-t128 - t74) * t467 + t562 * t466) * t346 + (t8 + (t288 * t122 + t289 * t124 + t344 * t181 + t224 * t196 + t225 * t198 + (-t454 + t556 - t562) * qJD(1)) * t344 + (t75 + t129) * t467 + t561 * t466 + (-t288 * t123 - t289 * t125 - t224 * t195 - t225 * t197 + t545 + (-t430 * t330 + t428 * t331 - t182) * t344 + (t184 * t330 - t186 * t331 + t254 * t493 + t256 * t496 - t547) * t346 + (t456 + (t253 - t387) * t344 + t561) * qJD(1)) * t346) * t344;
t95 = -rSges(6,3) * t437 + t446;
t96 = t361 * rSges(6,3) + t415;
t355 = t444 + (t170 + t191) * t466 + (t95 - t299 + (-t492 + t553) * t488 + (t506 + t552) * t467 + t474) * t346 + (t96 + t303 + (-t497 - t554) * t489 + (t356 * t346 + t321) * qJD(1)) * t344;
t27 = -t100 * t331 + t412 * t330;
t28 = -t101 * t331 + t411 * t330;
t147 = t395 * t493 + (Icges(6,3) * t336 + (-Icges(6,5) * t328 - Icges(6,6) * t329) * qJD(5)) * t330;
t148 = t399 * t493 + (Icges(6,6) * t336 + (-Icges(6,2) * t329 - t514) * qJD(5)) * t330;
t149 = t404 * t493 + (Icges(6,5) * t336 + (-Icges(6,1) * t328 - t513) * qJD(5)) * t330;
t33 = t147 * t494 + t148 * t267 + t149 * t268 + t151 * t232 + t152 * t233 + t360 * t231;
t3 = (t411 * t336 - t33) * t331 + (-qJD(1) * t45 + t336 * t101 + t15 * t344 + t16 * t346) * t330;
t17 = t153 * t166 + t154 * t168 + t361 * t164 + t265 * t92 + t266 * t94 + t90 * t495;
t18 = t153 * t167 + t154 * t169 + t361 * t165 + t265 * t91 + t266 * t93 + t89 * t495;
t34 = t147 * t495 + t148 * t265 + t149 * t266 + t153 * t232 + t154 * t233 + t361 * t231;
t4 = (t412 * t336 - t34) * t331 + (-qJD(1) * t44 + t100 * t336 + t17 * t344 + t18 * t346) * t330;
t9 = t412 * qJD(1) - t17 * t346 + t18 * t344;
t354 = t3 * t536 + t4 * t535 + t9 * t495 / 0.2e1 - t331 * (t409 * qJD(1) + t521 - t522) / 0.2e1 + t27 * t434 - t45 * t437 / 0.2e1 + t410 * t496 / 0.2e1 + t8 * t494 / 0.2e1 + (t344 * t44 + t346 * t45) * t493 / 0.2e1 + (t330 * t44 + t28) * t433;
t353 = rSges(3,2) * t436 + rSges(3,3) * t466 - t346 * t370;
t352 = -t331 * t147 + t231 * t496 + t493 * t505 + (t149 * t330 - t232 * t461) * t329;
t351 = t363 * t344 - t483;
t12 = (t373 * t123 - t287 * t125 - t226 * t195 - t227 * t197 + (t75 - t455) * qJD(1)) * t346 + (-t373 * t122 + t287 * t124 + t226 * t196 + t227 * t198 + (t74 + t453) * qJD(1) + t545) * t344;
t20 = (t346 * t182 + (t129 + t555) * qJD(1)) * t346 + (t128 * qJD(1) + (-t183 * t330 + t185 * t331 - t255 * t493 - t257 * t496 + t469) * t344 + (-t181 + t427 * t331 + t429 * t330 + (-t252 - t386) * qJD(1)) * t346) * t344;
t350 = (-t9 - t12 - t20) * t346 + t358;
t270 = t401 * t336;
t271 = t406 * t336;
t349 = qJD(1) * t291 + (t271 - t499) * t331 + (-t270 - t498) * t330;
t220 = (Icges(5,3) * t330 + t396 * t331) * t336;
t221 = (Icges(5,6) * t330 + t400 * t331) * t336;
t222 = (Icges(5,5) * t330 + t405 * t331) * t336;
t244 = -Icges(5,3) * t331 + t396 * t330;
t245 = -Icges(5,6) * t331 + t400 * t330;
t246 = -Icges(5,5) * t331 + t405 * t330;
t348 = -t522 / 0.2e1 + t521 / 0.2e1 + (t220 * t494 + t221 * t288 + t222 * t289 + t224 * t245 + t225 * t246 + t360 * t244 + t543 * t344 + t349 * t346 + t33 + (t390 * t336 - t120 + t430) * t331 + (-t122 * t340 + t124 * t341 + t194 * t336 + t428) * t330) * t536 + (t220 * t495 - t221 * t373 + t222 * t287 + t226 * t245 + t227 * t246 + t361 * t244 + t349 * t344 - t543 * t346 + t34 + (t391 * t336 - t121 + t429) * t331 + (-t123 * t340 + t125 * t341 + t193 * t336 - t427) * t330) * t535 + (t244 * t495 - t245 * t373 + t246 * t287 - t291 * t346 - t383 * t344 + t100 + t72 + (-t193 + t254) * t331 + (t256 + t391) * t330) * t434 + (t244 * t494 + t245 * t288 + t246 * t289 + t291 * t344 - t383 * t346 + t101 + t73 + (-t194 + t255) * t331 + (t257 + t390) * t330) * t433;
t320 = pkin(2) * t436;
t302 = t420 * qJD(2);
t282 = -t459 + t471;
t281 = t344 * t420 - t523;
t249 = t431 * t346;
t248 = t431 * t344;
t238 = t530 + (pkin(1) - t524) * t346 + t471;
t237 = t378 * t344 + t334 + t523;
t216 = t262 + t426;
t215 = (rSges(4,3) - t347) * t346 + t376 * t344;
t209 = t344 * t364 + t468;
t208 = -qJD(1) * t273 + t346 * t364;
t202 = t477 * t344;
t190 = t559 + ((-rSges(3,3) - pkin(6)) * t344 + t378 * t346) * qJD(1);
t189 = (t334 + (-pkin(1) - t526) * t344) * qJD(1) + t353;
t179 = t422 * t344;
t161 = -t295 * t466 - t504 + (-t343 * t466 - t344 * t463) * pkin(2);
t160 = t295 * t467 + t320 + (-t272 - t457) * t346;
t142 = t274 * t344 - t384 * t346;
t141 = t273 * t344 - t557;
t140 = -t274 * t346 - t558;
t139 = -t273 * t346 - t385 * t344;
t137 = t295 * t489 + (t376 * t346 - t332) * qJD(1) + t472;
t136 = (-t325 - t525) * t467 + (-t372 - t560) * t346 + t473;
t135 = t426 - t481;
t134 = t351 + t417;
t132 = t438 * t344;
t126 = t382 * t344;
t119 = -t171 * t331 - t234 * t494;
t118 = t170 * t331 + t234 * t495;
t115 = t377 + t426 + t171;
t114 = (pkin(4) * t340 - t347) * t346 + t359 * t344 + t414;
t113 = t178 + t478;
t112 = -t231 * t331 + (-t232 * t328 + t505) * t330;
t109 = qJD(1) * t203 + t480 * t344;
t108 = t480 * t346 + t479;
t107 = t112 * t496;
t106 = t392 * t330;
t103 = qJD(1) * t180 + t379 * t344;
t102 = t379 * t346 + t320 + t479;
t86 = -t262 * t467 + t442;
t83 = t303 + (t520 * t493 - t462) * t344 + t363 * t466 - t418 + t472;
t82 = (-pkin(3) * t496 - t458) * t346 + t351 * qJD(1) + t439 + t475;
t69 = t99 + t478;
t62 = qJD(1) * t133 + t445 * t344;
t61 = t445 * t346 + t440;
t60 = qJD(1) * t127 + t371 * t344;
t59 = t371 * t346 + t320 + t440;
t58 = (-t462 + (t519 * t331 + t497) * t336) * t344 + (t359 * t346 - t321) * qJD(1) - t415 + t472;
t57 = t310 + (-t458 + (-t492 - t497) * t336) * t346 + (-t483 + (-rSges(6,3) * t330 + t425) * t344) * qJD(1) + t446 + t474;
t55 = (-t251 - t262) * t467 + t441 + t442;
t54 = t56 + t478;
t51 = (t234 * t489 + t96) * t331 + (t150 * t344 - t170 * t336 + t234 * t466) * t330;
t50 = (-t234 * t488 - t95) * t331 + (-t150 * t346 + t171 * t336 + t228) * t330;
t46 = (-t452 + (-qJD(5) * t233 - t148) * t330) * t328 + t352;
t37 = t481 * t467 + t375;
t30 = t392 * t493 + (t549 * qJD(1) - t344 * t95 + t346 * t96) * t330;
t29 = (-t251 + t481) * t467 + t375 + t441;
t14 = t443 * t467 + t355;
t13 = (-t251 + t443) * t467 + t355 + t441;
t1 = [t352 + (t134 * t83 + t135 * t82) * t540 + (t114 * t58 + t115 * t57) * t539 + (t136 * t216 + t137 * t215) * t541 + (t189 * t238 + t190 * t237) * t542 + (t244 - t292) * t496 + (-t402 + t408) * t464 + (t407 + t403) * t463 + (t270 - t220) * t331 + (-t233 * t461 - t452) * t328 + (-t340 * t245 + t341 * t246 + t293) * t493 + (-t328 * t148 - t340 * t221 + t341 * t222 + t271) * t330; t348 + (-qJD(2) * t385 + (t276 * qJD(1) - t344 * t365) * t345 + (t278 * qJD(1) - t344 * t366) * t343) * t535 + (-qJD(2) * t384 + (-t275 * qJD(1) - t346 * t365) * t345 + (-t277 * qJD(1) - t346 * t366) * t343) * t536 + ((-t238 * t534 + t502 / 0.2e1 + t500 / 0.2e1) * t346 + (t503 / 0.2e1 + t501 / 0.2e1 + t237 * t534) * t344) * qJD(1) + (t337 / 0.2e1 + t338 / 0.2e1) * t398 * qJD(2) + m(3) * ((-t189 * t344 - t190 * t346) * t309 + (-t237 * t346 - t238 * t344) * t302) + m(6) * (t114 * t59 + t115 * t60 + t126 * t57 + t127 * t58) + m(5) * (t102 * t134 + t103 * t135 + t179 * t82 + t180 * t83) + m(4) * (t136 * t248 + t137 * t249 + t160 * t215 + t161 * t216); t358 + t344 * ((t344 * t208 + (t141 + t558) * qJD(1)) * t344 + (t142 * qJD(1) + (t275 * t463 + t277 * t464) * t346 + (-t209 + (-t500 - t502) * qJD(2) + (t274 - t385) * qJD(1)) * t344) * t346) - t346 * ((t346 * t209 + (t140 + t557) * qJD(1)) * t346 + (t139 * qJD(1) + (-t276 * t463 - t278 * t464 + t468) * t344 + (-t208 + (t501 + t503) * qJD(2) - t384 * qJD(1)) * t346) * t344) - t346 * t9 - t346 * t12 - t346 * t20 + ((t281 * t344 + t282 * t346) * ((qJD(1) * t281 + t353) * t346 + (-t559 + (-t282 - t459 + t333) * qJD(1)) * t344) + t470 * t309 * t302) * t542 + (t102 * t180 + t103 * t179 + t29 * t69) * t540 + (t113 * t55 + t160 * t249 + t161 * t248) * t541 + (t126 * t60 + t127 * t59 + t13 * t54) * t539 + (-t139 * t346 + t140 * t344) * t467 + (-t141 * t346 + t142 * t344) * t466; t348 + (-t136 * t344 - t137 * t346 + (t215 * t344 - t216 * t346) * qJD(1)) * t533 + m(6) * (t114 * t61 + t115 * t62 + t132 * t57 + t133 * t58) + m(5) * (t108 * t134 + t109 * t135 + t202 * t82 + t203 * t83) + m(4) * (-t215 * t346 - t216 * t344) * t272; t350 + (-t160 * t346 - t161 * t344 + (-t248 * t346 + t249 * t344) * qJD(1)) * t533 + m(6) * (t126 * t62 + t127 * t61 + t13 * t56 + t132 * t60 + t133 * t59 + t14 * t54) + m(5) * (t102 * t203 + t103 * t202 + t108 * t180 + t109 * t179 + t29 * t99 + t37 * t69) + m(4) * (-t249 * t272 * t346 + t113 * t86 + t178 * t55 - t248 * t504); (t132 * t62 + t133 * t61 + t14 * t56) * t539 + (t108 * t203 + t109 * t202 + t37 * t99) * t540 + (t272 * t295 * t470 + t178 * t86) * t541 + t350; 0.2e1 * ((t114 * t346 + t115 * t344) * t537 + (t134 * t346 + t135 * t344) * t538) * t493 + 0.2e1 * ((-t114 * t467 + t115 * t466 + t344 * t57 + t346 * t58) * t537 + (-t134 * t467 + t135 * t466 + t344 * t82 + t346 * t83) * t538) * t330; 0.2e1 * ((t126 * t489 + t127 * t488 - t13) * t537 + (t179 * t489 + t180 * t488 - t29) * t538) * t331 + 0.2e1 * ((t126 * t466 - t127 * t467 + t336 * t54 + t344 * t60 + t346 * t59) * t537 + (t102 * t346 + t103 * t344 + t179 * t466 - t180 * t467 + t336 * t69) * t538) * t330; 0.2e1 * ((t132 * t489 + t133 * t488 - t14) * t537 + (t202 * t489 + t203 * t488 - t37) * t538) * t331 + 0.2e1 * ((t132 * t466 - t133 * t467 + t336 * t56 + t344 * t62 + t346 * t61) * t537 + (t108 * t346 + t109 * t344 + t202 * t466 - t203 * t467 + t336 * t99) * t538) * t330; 0.4e1 * (t538 + t537) * (-0.1e1 + t470) * t330 * t493; m(6) * (t114 * t51 + t115 * t50 + t118 * t58 + t119 * t57) + t107 + (-t46 + (t344 * t448 + t346 * t447) * t336) * t331 + ((t33 / 0.2e1 + t24 / 0.2e1) * t346 + (t34 / 0.2e1 + t23 / 0.2e1) * t344 + (-t344 * t447 + t346 * t448) * qJD(1)) * t330; m(6) * (t106 * t13 + t118 * t59 + t119 * t60 + t126 * t50 + t127 * t51 + t30 * t54) + t354; m(6) * (t106 * t14 + t118 * t61 + t119 * t62 + t132 * t50 + t133 * t51 + t30 * t56) + t354; m(6) * ((-t30 + (t118 * t346 + t119 * t344) * t336) * t331 + (t106 * t336 + t344 * t50 + t346 * t51 + (-t118 * t344 + t119 * t346) * qJD(1)) * t330); (t106 * t30 + t118 * t51 + t119 * t50) * t539 + (t46 * t331 - t107 + (t344 * t27 + t346 * t28 - t331 * t409) * t336) * t331 + (t346 * t3 + t344 * t4 + t409 * t496 + (-t112 * t336 - t23 * t344 - t24 * t346) * t331 + (t346 * t27 - t344 * t28 + t331 * t410) * qJD(1)) * t330;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;