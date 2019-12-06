% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:07
% EndTime: 2019-12-05 17:44:37
% DurationCPUTime: 12.68s
% Computational Cost: add. (40177->452), mult. (39793->669), div. (0->0), fcn. (43098->10), ass. (0->291)
t351 = cos(qJ(1));
t348 = sin(pkin(8));
t486 = -pkin(6) - qJ(3);
t395 = t348 * t486;
t336 = t351 * t395;
t349 = cos(pkin(8));
t477 = cos(pkin(9));
t386 = t477 * pkin(3) + pkin(2);
t374 = t349 * t386;
t419 = -pkin(7) + t486;
t389 = t348 * t419;
t350 = sin(qJ(1));
t476 = sin(pkin(9));
t392 = t350 * t476;
t347 = pkin(9) + qJ(4);
t503 = cos(t347);
t416 = pkin(4) * t503;
t359 = t349 * (t416 + t386);
t404 = t476 * pkin(3);
t343 = sin(t347);
t487 = pkin(4) * t343;
t370 = t404 + t487;
t424 = -t350 * t370 - t351 * t359;
t196 = pkin(3) * t392 - t336 + (t374 + t389) * t351 + t424;
t464 = t348 * t351;
t394 = qJ(5) + t347;
t385 = cos(t394);
t373 = t350 * t385;
t342 = sin(t394);
t458 = t351 * t342;
t295 = t349 * t458 - t373;
t372 = t351 * t385;
t460 = t350 * t342;
t296 = t349 * t372 + t460;
t541 = -t296 * rSges(6,1) + t295 * rSges(6,2);
t223 = rSges(6,3) * t464 - t541;
t211 = t349 * t223;
t354 = t348 * (-t349 * rSges(6,3) + (rSges(6,1) * t385 - rSges(6,2) * t342) * t348);
t265 = t351 * t354;
t362 = t348 * (-pkin(7) * t349 + t348 * t416);
t433 = t351 * t362 + t265;
t105 = -t196 * t349 + t211 + t433;
t568 = t196 - t223;
t456 = -t568 * t349 - t105 + t433;
t573 = m(6) * t456;
t344 = t351 * qJ(2);
t423 = t350 * t389 + t351 * t370;
t293 = t349 * t460 + t372;
t294 = t349 * t373 - t458;
t428 = t294 * rSges(6,1) - t293 * rSges(6,2);
t163 = -t344 + (t348 * rSges(6,3) + pkin(1) + t359) * t350 - t423 + t428;
t461 = t350 * qJ(2);
t164 = (-pkin(1) + (-rSges(6,3) + t419) * t348) * t351 - t461 + t424 + t541;
t244 = t293 * rSges(6,1) + t294 * rSges(6,2);
t245 = -t295 * rSges(6,1) - t296 * rSges(6,2);
t574 = m(6) * (-t163 * t244 - t164 * t245);
t474 = Icges(5,4) * t343;
t280 = -Icges(5,5) * t349 + (t503 * Icges(5,1) - t474) * t348;
t308 = (-t503 * Icges(5,2) - t474) * t348;
t407 = t503 * Icges(5,4);
t309 = (-Icges(5,1) * t343 - t407) * t348;
t405 = t503 * t309;
t279 = -Icges(5,6) * t349 + (-Icges(5,2) * t343 + t407) * t348;
t406 = t503 * t279;
t307 = (-Icges(5,5) * t343 - t503 * Icges(5,6)) * t348;
t462 = t349 * t307;
t572 = t462 / 0.2e1 - (-(t280 / 0.2e1 + t308 / 0.2e1) * t343 - t406 / 0.2e1 + t405 / 0.2e1) * t348;
t175 = t211 + t265;
t214 = Icges(6,5) * t296 - Icges(6,6) * t295 + Icges(6,3) * t464;
t465 = t348 * t350;
t473 = Icges(6,4) * t294;
t216 = Icges(6,2) * t293 - Icges(6,6) * t465 - t473;
t281 = Icges(6,4) * t293;
t219 = -Icges(6,1) * t294 - Icges(6,5) * t465 + t281;
t447 = -t293 * t216 + t294 * t219;
t571 = t214 * t464 - t447;
t534 = t350 ^ 2;
t570 = t348 * t534;
t525 = m(6) / 0.2e1;
t569 = t456 * t525;
t409 = t350 * t503;
t457 = t351 * t343;
t312 = t349 * t457 - t409;
t408 = t351 * t503;
t459 = t350 * t343;
t313 = t349 * t408 + t459;
t227 = Icges(5,5) * t313 - Icges(5,6) * t312 + Icges(5,3) * t464;
t310 = t349 * t459 + t408;
t388 = t349 * t409;
t311 = t388 - t457;
t475 = Icges(5,4) * t311;
t229 = Icges(5,2) * t310 - Icges(5,6) * t465 - t475;
t298 = Icges(5,4) * t310;
t232 = -Icges(5,1) * t311 - Icges(5,5) * t465 + t298;
t444 = -t310 * t229 + t311 * t232;
t567 = t227 * t464 - t444;
t540 = -t313 * rSges(5,1) + t312 * rSges(5,2);
t236 = rSges(5,3) * t464 - t540;
t356 = t348 * (-t349 * rSges(5,3) + (t503 * rSges(5,1) - rSges(5,2) * t343) * t348);
t183 = t236 * t349 + t351 * t356;
t549 = -t348 / 0.2e1;
t472 = Icges(6,4) * t342;
t277 = -Icges(6,5) * t349 + (Icges(6,1) * t385 - t472) * t348;
t431 = t277 + (-Icges(6,2) * t385 - t472) * t348;
t566 = t342 * t431;
t302 = t312 * pkin(4);
t205 = -t245 + t302;
t434 = -t310 * pkin(4) - t244;
t564 = t205 * t350 + t434 * t351;
t283 = Icges(6,4) * t296;
t217 = -Icges(6,2) * t295 + Icges(6,6) * t464 + t283;
t282 = Icges(6,4) * t295;
t221 = -Icges(6,1) * t296 - Icges(6,5) * t464 + t282;
t448 = t293 * t217 + t221 * t294;
t300 = Icges(5,4) * t313;
t230 = -Icges(5,2) * t312 + Icges(5,6) * t464 + t300;
t299 = Icges(5,4) * t312;
t234 = -Icges(5,1) * t313 - Icges(5,5) * t464 + t299;
t445 = t310 * t230 + t234 * t311;
t527 = m(5) / 0.2e1;
t547 = m(6) * t348;
t561 = t547 / 0.2e1;
t213 = -Icges(6,5) * t294 + Icges(6,6) * t293 - Icges(6,3) * t465;
t468 = t213 * t350;
t96 = -t213 * t465 - t447;
t97 = -t214 * t465 + t448;
t560 = (-t214 * t570 + (t448 - t97) * t350 + (-t96 + (-t214 * t351 - t468) * t348 + t571) * t351) * t348;
t358 = t348 * rSges(5,3) + pkin(1) + t374;
t181 = t336 + (-t404 - qJ(2)) * t350 - t358 * t351 + t540;
t264 = t350 * t354;
t222 = -rSges(6,3) * t465 - t428;
t390 = t351 * t476;
t422 = -pkin(3) * t390 - t350 * t395;
t446 = -pkin(4) * t388 + t222 + t422 + t423;
t103 = t446 * t349 - t350 * t362 - t264;
t426 = t311 * rSges(5,1) - t310 * rSges(5,2);
t235 = -rSges(5,3) * t465 - t426;
t182 = t235 * t349 - t350 * t356;
t478 = (-t182 * t351 - t183 * t350) * t527 + (-t103 * t351 - t105 * t350) * t525;
t484 = t350 * t573 / 0.2e1;
t552 = t478 - t484;
t548 = m(5) * t348;
t479 = (t182 * t350 - t183 * t351) * t548 / 0.2e1 + (t103 * t350 - t105 * t351) * t561;
t485 = t464 * t569;
t551 = t479 - t485;
t290 = (-Icges(6,5) * t342 - Icges(6,6) * t385) * t348;
t463 = t349 * t290;
t550 = -t463 / 0.2e1 + t566 * t549;
t504 = -t349 / 0.2e1;
t159 = (t244 * t351 + t245 * t350) * t348;
t297 = (-rSges(6,1) * t342 - rSges(6,2) * t385) * t348;
t269 = t297 * t465;
t185 = -t244 * t349 + t269;
t186 = t349 * t245 + t297 * t464;
t371 = t385 * Icges(6,4);
t276 = -Icges(6,6) * t349 + (-Icges(6,2) * t342 + t371) * t348;
t292 = (-Icges(6,1) * t342 - t371) * t348;
t432 = t276 - t292;
t100 = -t290 * t465 + t293 * t431 + t294 * t432;
t101 = t290 * t464 - t295 * t431 - t296 * t432;
t238 = Icges(6,5) * t293 + Icges(6,6) * t294;
t239 = -Icges(6,5) * t295 - Icges(6,6) * t296;
t398 = t464 / 0.2e1;
t401 = -t465 / 0.2e1;
t439 = -Icges(6,2) * t296 - t221 - t282;
t440 = Icges(6,2) * t294 + t219 + t281;
t441 = Icges(6,1) * t295 + t217 + t283;
t442 = -Icges(6,1) * t293 + t216 - t473;
t135 = -t463 + (-t432 * t385 - t566) * t348;
t471 = t135 * t349;
t80 = -t349 * t238 + (-t440 * t342 - t442 * t385) * t348;
t81 = -t349 * t239 + (-t439 * t342 - t441 * t385) * t348;
t418 = (-t471 + (-t350 * t80 + t351 * t81) * t348) * t504 + (-t100 * t349 - (-t238 * t465 + t440 * t293 + t442 * t294) * t465 + (-t239 * t465 + t439 * t293 + t441 * t294) * t464) * t401 + (-t101 * t349 - (t238 * t464 - t440 * t295 - t442 * t296) * t465 + (t239 * t464 - t439 * t295 - t441 * t296) * t464) * t398;
t86 = (t568 * t350 - t446 * t351) * t348;
t6 = t418 + m(6) * (-t103 * t185 + t105 * t186 - t86 * t159);
t544 = t6 * qJD(5);
t254 = rSges(5,1) * t310 + rSges(5,2) * t311;
t255 = -rSges(5,1) * t312 - rSges(5,2) * t313;
t169 = (t254 * t351 + t255 * t350) * t348;
t438 = -Icges(5,1) * t310 + t229 - t475;
t437 = Icges(5,1) * t312 + t230 + t300;
t453 = (t205 * t351 - t350 * t434) * t525 + (t350 * t254 - t255 * t351) * t527;
t454 = t169 * t527 - t561 * t564;
t536 = t349 * pkin(2) + pkin(1) + (rSges(4,3) + qJ(3)) * t348;
t328 = (t351 ^ 2 + t534) * t348;
t533 = 0.2e1 * t328;
t531 = 4 * qJD(1);
t530 = 2 * qJD(4);
t529 = 4 * qJD(4);
t526 = m(5) / 0.4e1;
t524 = m(6) / 0.4e1;
t180 = t350 * t358 - t344 + t422 + t426;
t520 = m(5) * (-t180 * t254 - t181 * t255);
t174 = t222 * t349 - t264;
t519 = t174 * t573;
t518 = t103 * t573;
t455 = -t185 * t163 + t186 * t164;
t515 = m(6) * (-t103 * t244 - t105 * t245 + t455);
t513 = m(6) * (t174 * t434 + t175 * t205 + t455);
t508 = m(6) * (t163 * t434 + t164 * t205);
t507 = (t163 * t350 - t164 * t351) * t547;
t506 = m(6) * (-t351 * t163 - t164 * t350);
t502 = m(3) * (-((-rSges(3,3) - qJ(2)) * t350 + rSges(3,2) * t464) * t350 - t351 * (-rSges(3,2) * t465 - t351 * rSges(3,3) - t344));
t391 = t351 * t477;
t393 = t350 * t477;
t207 = -(t349 * t392 + t391) * rSges(4,2) - (-t349 * t393 + t390) * rSges(4,1) - t344 + t536 * t350;
t353 = -(t349 * t391 + t392) * rSges(4,1) + (t349 * t390 - t393) * rSges(4,2) - t461 - t536 * t351;
t501 = m(4) * (t207 * t350 - t351 * t353) * t348;
t500 = m(4) * (-t351 * t207 - t350 * t353);
t499 = (t180 * t350 - t181 * t351) * t548;
t497 = m(5) * (-t351 * t180 - t181 * t350);
t493 = (t174 * t350 - t175 * t351) * t547;
t492 = m(6) * (-t174 * t351 - t175 * t350);
t489 = m(6) * t159;
t488 = m(6) * (t350 * t244 - t245 * t351);
t483 = m(6) * qJD(4);
t482 = m(6) * qJD(5);
t470 = (-(-Icges(6,3) * t349 + (Icges(6,5) * t385 - Icges(6,6) * t342) * t348) * t465 + t276 * t293 - t277 * t294) * t349;
t226 = -Icges(5,5) * t311 + Icges(5,6) * t310 - Icges(5,3) * t465;
t467 = t226 * t350;
t436 = Icges(5,2) * t311 + t232 + t298;
t435 = -Icges(5,2) * t313 - t234 - t299;
t430 = t279 - t309;
t429 = t280 + t308;
t167 = (m(4) / 0.2e1 + t527 + t525) * t533;
t420 = t167 * qJD(1);
t417 = t348 ^ 2 * t487;
t107 = -t226 * t465 - t444;
t108 = -t227 * t465 + t445;
t413 = (-t227 * t570 + (-t108 + t445) * t350 + (-t107 + (-t227 * t351 - t467) * t348 + t567) * t351) * t549 + (t504 + t349 / 0.2e1) * ((-Icges(5,3) * t349 + (t503 * Icges(5,5) - Icges(5,6) * t343) * t348) * t464 - t312 * t279 + t313 * t280);
t412 = (-t107 * t350 + t108 * t351) * t549 + (t445 * t351 + (t467 * t348 - t567) * t350) * t348 / 0.2e1;
t399 = -t464 / 0.2e1;
t21 = -t470 + (t448 * t351 + (t468 * t348 - t571) * t350) * t348;
t51 = -t470 + (-t350 * t96 + t351 * t97) * t348;
t387 = t21 * t398 + t51 * t399 + t560 * t401;
t375 = t348 * t385;
t364 = t375 / 0.2e1;
t365 = -t375 / 0.2e1;
t381 = t276 * t365 + t292 * t364 + t550;
t380 = -t519 / 0.2e1 + t387;
t363 = t21 * t399 + (t80 + t100) * t401 + t560 * t465 / 0.2e1 + (t51 + t81 + t101) * t398;
t360 = t276 * t364 + t292 * t365 - t550;
t357 = t363 - t471;
t314 = (-rSges(5,1) * t343 - t503 * rSges(5,2)) * t348;
t249 = -Icges(5,5) * t312 - Icges(5,6) * t313;
t248 = Icges(5,5) * t310 + Icges(5,6) * t311;
t194 = t349 * t255 + t314 * t464;
t193 = -t254 * t349 + t314 * t465;
t168 = (m(4) / 0.4e1 + t526 + t524) * t533 - (m(4) + m(5) + m(6)) * t328 / 0.2e1;
t161 = t488 / 0.2e1;
t157 = t489 / 0.2e1;
t154 = -t349 * t302 - t351 * t417 + t186;
t153 = t349 * t434 - t350 * t417 + t269;
t148 = (-t222 * t351 - t223 * t350) * t348;
t139 = -t462 + (-t343 * t429 + t405 - t406) * t348;
t134 = t564 * t348;
t133 = t185 * t350 + t186 * t351;
t130 = t133 * t482;
t116 = t492 / 0.2e1;
t115 = t307 * t464 - t312 * t429 - t313 * t430;
t114 = -t307 * t465 + t310 * t429 + t311 * t430;
t106 = t493 / 0.2e1;
t89 = -t349 * t249 + (-t435 * t343 - t437 * t503) * t348;
t88 = -t349 * t248 + (-t436 * t343 - t438 * t503) * t348;
t85 = t159 * t349 + (t185 * t351 - t186 * t350) * t348;
t84 = t85 * t482;
t62 = t381 + t574;
t54 = t499 + t501 + t507;
t45 = t497 + t500 + t502 + t506;
t44 = t116 - t488 / 0.2e1;
t43 = t161 + t116;
t42 = t161 - t492 / 0.2e1;
t40 = t513 / 0.2e1;
t38 = t106 - t489 / 0.2e1;
t37 = t157 + t106;
t36 = t157 - t493 / 0.2e1;
t32 = t515 / 0.2e1;
t31 = t520 + t508 + t381 - t572;
t14 = t478 + t484 - t453;
t13 = t453 + t552;
t12 = t453 - t552;
t11 = t479 + t485 - t454;
t10 = t454 + t551;
t9 = t454 - t551;
t8 = m(6) * (-t148 * t159 - t174 * t185 + t175 * t186) + t418;
t7 = t8 * qJD(5);
t4 = t40 - t515 / 0.2e1 + t380;
t3 = t32 - t513 / 0.2e1 + t380;
t2 = t32 + t40 + t519 / 0.2e1 + t357;
t1 = -t518 + (t350 * t413 + t351 * t412) * t348 + t387;
t5 = [t45 * qJD(2) + t54 * qJD(3) + t31 * qJD(4) + t62 * qJD(5), qJD(1) * t45 + qJD(3) * t168 + qJD(4) * t13 + qJD(5) * t43, qJD(1) * t54 + qJD(2) * t168 + qJD(4) * t10 + qJD(5) * t37, t31 * qJD(1) + t13 * qJD(2) + t10 * qJD(3) + t2 * qJD(5) + t518 * t529 / 0.4e1 + ((-t180 * t193 + t181 * t194 - t182 * t254 - t183 * t255) * t527 + (t103 * t434 + t105 * t205 - t153 * t163 + t154 * t164) * t525) * t530 + (t363 + (-t139 - t135) * t349 + ((t89 / 0.2e1 + t115 / 0.2e1 - t412) * t351 + (-t88 / 0.2e1 - t114 / 0.2e1 - t413) * t350) * t348) * qJD(4), t62 * qJD(1) + t43 * qJD(2) + t37 * qJD(3) + t2 * qJD(4) + ((-t174 * t244 - t175 * t245 + t455) * m(6) + t357) * qJD(5); -t167 * qJD(3) + t12 * qJD(4) + t42 * qJD(5) + (-t502 / 0.4e1 - t500 / 0.4e1 - t497 / 0.4e1 - t506 / 0.4e1) * t531, 0, -t420, t12 * qJD(1) + ((t193 * t350 + t194 * t351) * t527 + (t153 * t350 + t154 * t351) * t525) * t530 + t130, t42 * qJD(1) + t133 * t483 + t130; t167 * qJD(2) + t9 * qJD(4) + t36 * qJD(5) + (-t501 / 0.4e1 - t499 / 0.4e1 - t507 / 0.4e1) * t531, t420, 0, t9 * qJD(1) + ((t169 * t349 + (t193 * t351 - t194 * t350) * t348) * t527 + (-t134 * t349 + (t153 * t351 - t154 * t350) * t348) * t525) * t530 + t84, t36 * qJD(1) + t85 * t483 + t84; t14 * qJD(2) + t11 * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + (-t520 / 0.4e1 - t508 / 0.4e1) * t531 + (-0.2e1 * t163 * t569 + t360 + t572) * qJD(1), t14 * qJD(1), t11 * qJD(1), t1 * qJD(1) + ((-t139 * t349 + (-t350 * t88 + t351 * t89) * t348) * t504 + (-t114 * t349 - (-t248 * t465 + t436 * t310 + t438 * t311) * t465 + (-t249 * t465 + t435 * t310 + t437 * t311) * t464) * t401 + (-t115 * t349 - (t248 * t464 - t436 * t312 - t438 * t313) * t465 + (t249 * t464 - t435 * t312 - t437 * t313) * t464) * t398 + t418) * qJD(4) + t544 + ((-(-t235 * t351 - t236 * t350) * t348 * t169 - t182 * t193 + t183 * t194) * t526 + (-t103 * t153 + t105 * t154 + t134 * t86) * t524) * t529, t3 * qJD(1) + t6 * qJD(4) + t544; (t360 - t574) * qJD(1) + t44 * qJD(2) + t38 * qJD(3) + t4 * qJD(4) + t387 * qJD(5), t44 * qJD(1), t38 * qJD(1), t4 * qJD(1) + ((t134 * t148 - t153 * t174 + t154 * t175) * m(6) + t418) * qJD(4) + t7, qJD(1) * t387 + qJD(4) * t8 + t7;];
Cq = t5;
