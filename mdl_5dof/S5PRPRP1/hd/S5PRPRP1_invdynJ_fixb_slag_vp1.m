% Calculate vector of inverse dynamics joint torques for
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:18
% EndTime: 2019-12-05 15:28:38
% DurationCPUTime: 13.46s
% Computational Cost: add. (10066->460), mult. (8525->570), div. (0->0), fcn. (6585->6), ass. (0->243)
t503 = Icges(6,4) + Icges(5,5);
t242 = pkin(8) + qJ(4);
t240 = cos(t242);
t238 = sin(t242);
t388 = Icges(5,4) * t238;
t170 = Icges(5,1) * t240 - t388;
t223 = Icges(6,5) * t238;
t288 = Icges(6,1) * t240 + t223;
t488 = t170 + t288;
t243 = pkin(7) + qJ(2);
t239 = sin(t243);
t241 = cos(t243);
t502 = -t239 * t288 + t503 * t241;
t371 = t238 * t239;
t192 = Icges(5,4) * t371;
t369 = t239 * t240;
t493 = -Icges(5,1) * t369 + t192 + t502;
t492 = t503 * t239 + t488 * t241;
t286 = Icges(6,3) * t240 - t223;
t501 = Icges(5,2) * t240 + t286 + t388;
t385 = Icges(6,5) * t240;
t167 = Icges(6,1) * t238 - t385;
t226 = Icges(5,4) * t240;
t499 = Icges(5,1) * t238 + t167 + t226;
t162 = Icges(5,5) * t240 - Icges(5,6) * t238;
t107 = Icges(5,3) * t239 + t162 * t241;
t164 = Icges(6,4) * t240 + Icges(6,6) * t238;
t109 = Icges(6,2) * t239 + t164 * t241;
t500 = t107 + t109;
t379 = Icges(5,3) * t241;
t106 = Icges(5,5) * t369 - Icges(5,6) * t371 - t379;
t160 = Icges(6,3) * t238 + t385;
t104 = -Icges(6,6) * t241 + t160 * t239;
t382 = Icges(5,6) * t241;
t110 = Icges(5,4) * t369 - Icges(5,2) * t371 - t382;
t377 = t110 * t238;
t462 = -t104 * t238 + t493 * t240 + t377;
t498 = -t241 * t106 - t462 * t239;
t368 = t240 * t241;
t191 = Icges(6,5) * t368;
t370 = t238 * t241;
t381 = Icges(6,6) * t239;
t105 = Icges(6,3) * t370 + t191 + t381;
t497 = -t105 * t371 + t109 * t241 - t492 * t369;
t495 = t104 - t110;
t287 = -Icges(5,2) * t238 + t226;
t111 = Icges(5,6) * t239 + t241 * t287;
t494 = t105 - t111;
t491 = t160 - t287;
t490 = (Icges(5,6) - Icges(6,6)) * t240 + t503 * t238;
t487 = t501 * qJD(4);
t486 = t499 * qJD(4);
t485 = t105 * t370 + t500 * t239 + t492 * t368;
t108 = -Icges(6,2) * t241 + t164 * t239;
t98 = t239 * t108;
t484 = -t104 * t370 - t239 * t106 + t493 * t368 - t98;
t472 = -t238 * t501 + t499 * t240;
t483 = -t241 * t107 - t497;
t366 = t241 * t108;
t440 = -t366 + t498;
t438 = -t110 * t370 - t484;
t437 = -t111 * t370 + t485;
t482 = t487 * t241 + (t239 * t287 - t104 - t382) * qJD(2);
t481 = t487 * t239 + (t160 * t241 - t111 + t381) * qJD(2);
t480 = -t486 * t241 + (-t170 * t239 + t502) * qJD(2);
t479 = -qJD(2) * t492 + t239 * t486;
t478 = -t238 * t492 + t240 * t494;
t436 = t238 * t493 + t240 * t495;
t439 = -t111 * t371 + t483;
t477 = t491 * qJD(4);
t476 = t488 * qJD(4);
t475 = -t162 - t164;
t473 = t490 * qJD(4);
t471 = -t238 * t499 - t240 * t501;
t376 = t111 * t238;
t470 = t105 * t238 + t240 * t492 - t376;
t430 = t490 * t241;
t431 = t490 * t239;
t434 = t472 * t239 - t430;
t433 = t472 * t241 + t431;
t469 = t500 * qJD(2);
t468 = (-t241 * t501 + t492) * t239 - (-Icges(5,2) * t369 - t286 * t239 - t192 - t493) * t241;
t467 = t241 ^ 2;
t451 = rSges(6,3) + qJ(5);
t466 = qJD(2) * t490 + t471 * qJD(4) + t477 * t238 + t476 * t240;
t465 = t478 * qJD(4) + t482 * t238 + t480 * t240 + t469;
t420 = qJD(2) * t108;
t464 = -qJD(2) * t106 - t436 * qJD(4) - t481 * t238 + t479 * t240 - t420;
t463 = -t495 * t241 + (-Icges(6,1) * t370 + t167 * t241 + t191 + t494) * t239;
t461 = t437 * t239 - t438 * t241;
t460 = t439 * t239 - t440 * t241;
t459 = t499 - t491;
t458 = -t501 + t488;
t455 = t472 * qJD(2) + t475 * qJD(4);
t454 = t462 * qJD(2) - t473 * t239 + t469;
t453 = -t420 - t473 * t241 + (-t162 * t239 + t379 - t470) * qJD(2);
t409 = rSges(6,1) + pkin(4);
t452 = t433 * qJD(2);
t299 = t240 * rSges(6,1) + t238 * rSges(6,3);
t450 = t240 * pkin(4) + t238 * qJ(5) + t299;
t449 = t434 * qJD(2);
t448 = -m(2) - m(3);
t324 = qJD(2) * qJD(4);
t156 = -qJDD(4) * t241 + t239 * t324;
t327 = qJD(5) * t240;
t355 = -t450 * qJD(4) + t327;
t259 = qJDD(5) * t238 + (t327 + t355) * qJD(4);
t234 = t241 * rSges(6,2);
t357 = t450 * t239 - t234;
t218 = t241 * qJ(3);
t174 = pkin(2) * t239 - t218;
t246 = -pkin(6) - qJ(3);
t208 = t241 * t246;
t245 = cos(pkin(8));
t237 = pkin(3) * t245 + pkin(2);
t342 = -t239 * t237 - t208;
t102 = t174 + t342;
t363 = t102 - t174;
t305 = -t357 + t363;
t328 = qJD(5) * t238;
t317 = t239 * t328;
t325 = qJD(2) * qJD(3);
t338 = qJDD(3) * t239 + t241 * t325;
t172 = rSges(6,1) * t238 - rSges(6,3) * t240;
t348 = pkin(4) * t238 - qJ(5) * t240 + t172;
t217 = t239 * qJ(3);
t179 = t241 * pkin(2) + t217;
t215 = qJD(3) * t241;
t136 = t179 * qJD(2) - t215;
t332 = qJD(2) * t239;
t203 = t246 * t332;
t405 = pkin(2) - t237;
t391 = -t136 + t203 - (-t241 * t405 - t217) * qJD(2);
t230 = t239 * rSges(6,2);
t330 = qJD(4) * t239;
t331 = qJD(2) * t241;
t403 = (pkin(4) * t331 + qJ(5) * t330) * t240 + (qJ(5) * t331 + (-pkin(4) * qJD(4) + qJD(5)) * t239) * t238 - t172 * t330 + (t241 * t299 + t230) * qJD(2);
t1 = t348 * t156 + t259 * t241 + t305 * qJDD(2) + (-t317 + t391 - t403) * qJD(2) + t338;
t447 = -g(1) + t1;
t446 = t460 * qJD(4) + t449;
t445 = t461 * qJD(4) + t452;
t444 = t462 * qJD(4) + t479 * t238 + t481 * t240;
t443 = t470 * qJD(4) + t480 * t238 - t482 * t240;
t442 = -t455 * t239 + t466 * t241;
t441 = t466 * t239 + t455 * t241;
t158 = qJD(2) * t174;
t432 = qJD(2) * t102 - t158;
t329 = qJD(4) * t241;
t318 = t240 * t329;
t429 = rSges(6,2) * t331 + t451 * t318;
t428 = t451 * t369;
t427 = t451 * t368;
t426 = -t238 * t468 + t463 * t240;
t425 = (-t459 * t238 + t458 * t240) * qJD(2);
t424 = t366 + t485;
t423 = t454 * t467 + (t465 * t239 + (-t453 + t464) * t241) * t239;
t422 = t464 * t467 + (t453 * t239 + (-t454 + t465) * t241) * t239;
t421 = t475 * qJD(2);
t244 = sin(pkin(8));
t398 = rSges(4,2) * t244;
t400 = rSges(4,1) * t245;
t121 = t239 * rSges(4,3) + (-t398 + t400) * t241;
t155 = qJDD(4) * t239 + t241 * t324;
t413 = t155 / 0.2e1;
t412 = t156 / 0.2e1;
t411 = t239 / 0.2e1;
t410 = -t241 / 0.2e1;
t408 = g(2) * t239;
t185 = t241 * t328;
t263 = -t238 * t329 - t240 * t332;
t319 = t238 * t332;
t404 = t409 * t263 - t451 * t319 + t185 + t429;
t399 = rSges(5,1) * t240;
t173 = rSges(5,1) * t238 + rSges(5,2) * t240;
t143 = t173 * t241;
t228 = t239 * rSges(5,3);
t119 = rSges(5,1) * t368 - rSges(5,2) * t370 + t228;
t307 = t241 * t237 - t239 * t246;
t103 = t307 - t179;
t362 = t103 + t179;
t45 = -t173 * t330 - t215 + (t119 + t362) * qJD(2);
t397 = t143 * t45;
t214 = qJD(3) * t239;
t300 = -t173 * t329 + t214;
t117 = rSges(5,1) * t369 - rSges(5,2) * t371 - t241 * rSges(5,3);
t321 = -t117 + t363;
t44 = qJD(2) * t321 + t300;
t395 = t239 * t44;
t356 = t409 * t368 + t451 * t370 + t230;
t31 = -t327 + qJD(1) + (t239 * t357 + t241 * t356) * qJD(4);
t394 = t31 * t238;
t94 = t121 + t179;
t354 = -t409 * t371 + t428;
t353 = -t409 * t370 + t427;
t323 = t239 * t400;
t204 = t239 * t398;
t339 = t241 * rSges(4,3) + t204;
t120 = t323 - t339;
t347 = -t174 - t120;
t344 = rSges(5,2) * t319 + rSges(5,3) * t331;
t343 = t185 + t214;
t341 = rSges(4,3) * t331 + qJD(2) * t204;
t340 = t203 + t215;
t209 = qJ(3) * t331;
t337 = t209 + t214;
t326 = -m(4) - m(5) - m(6);
t320 = t409 * t241;
t316 = -pkin(2) - t400;
t313 = -t330 / 0.2e1;
t312 = t330 / 0.2e1;
t311 = -t329 / 0.2e1;
t310 = t329 / 0.2e1;
t308 = -t106 + t376;
t306 = t348 * qJD(4);
t180 = rSges(3,1) * t241 - rSges(3,2) * t239;
t175 = rSges(3,1) * t239 + rSges(3,2) * t241;
t178 = -rSges(5,2) * t238 + t399;
t264 = -t241 * t306 + t343;
t29 = qJD(2) * t305 + t264;
t30 = -t215 + (-t306 + t328) * t239 + (t356 + t362) * qJD(2);
t296 = t239 * t30 + t241 * t29;
t291 = -t239 * t45 - t241 * t44;
t74 = rSges(5,1) * t263 - rSges(5,2) * t318 + t344;
t139 = t173 * t239;
t76 = -qJD(4) * t139 + (t178 * t241 + t228) * qJD(2);
t290 = t239 * t76 + t241 * t74;
t279 = t117 * t239 + t119 * t241;
t273 = -qJDD(3) * t241 + qJD(2) * (-pkin(2) * t332 + t337) + qJDD(2) * t179 + t239 * t325;
t272 = t296 * t240;
t262 = qJDD(2) * t103 + t273 + qJD(2) * (-t209 + (t239 * t405 - t208) * qJD(2));
t256 = -t237 - t450;
t152 = t178 * qJD(4);
t80 = qJD(2) * t94 - t215;
t79 = qJD(2) * t347 + t214;
t46 = qJD(4) * t279 + qJD(1);
t35 = qJDD(2) * t121 + qJD(2) * (-qJD(2) * t323 + t341) + t273;
t34 = t347 * qJDD(2) + (-t121 * qJD(2) - t136) * qJD(2) + t338;
t18 = qJD(4) * t290 + t117 * t155 - t119 * t156 + qJDD(1);
t17 = qJD(2) * t74 + qJDD(2) * t119 - t152 * t330 - t155 * t173 + t262;
t16 = -t152 * t329 + t156 * t173 + t321 * qJDD(2) + (-t76 + t391) * qJD(2) + t338;
t3 = -qJDD(5) * t240 + qJDD(1) - t356 * t156 + t357 * t155 + (t239 * t403 + t241 * t404 + t328) * qJD(4);
t2 = -t348 * t155 + t356 * qJDD(2) + (t185 + t404) * qJD(2) + t259 * t239 + t262;
t4 = [m(5) * t18 + m(6) * t3 + (m(4) - t448) * qJDD(1) + (t326 + t448) * g(3); -m(3) * (-g(1) * t175 + g(2) * t180) + ((((t107 + t377) * t241 + t439 + t484 + t497) * t241 + (t424 + t440 - t498) * t239) * qJD(4) + t452) * t310 + (t472 * qJD(4) + t476 * t238 - t477 * t240) * qJD(2) + (t29 * (-t317 + t340) + t30 * (t343 + t429) + (-t30 * t238 * t320 + (t238 * t409 - t240 * t451) * t239 * t29) * qJD(4) + ((-t30 * t246 + t256 * t29) * t241 + (-t29 * rSges(6,2) + t256 * t30) * t239) * qJD(2) - (-qJD(2) * t357 + t264 - t29 + t432) * t30 + (t2 - g(2)) * (t307 + t356) + t447 * (t234 + (-t238 * t451 - t240 * t409) * t239 + t342)) * m(6) + (-(-qJD(2) * t117 + t300 + t432 - t44) * t45 + t44 * t340 + t45 * (t214 + t344) + (t173 * t395 - t397) * qJD(4) + ((-t44 * rSges(5,3) + t45 * (-t237 - t399)) * t239 + (t44 * (-t178 - t237) - t45 * t246) * t241) * qJD(2) + (t17 - g(2)) * (t119 + t307) + (t16 - g(1)) * (-t117 + t342)) * m(5) + (-(-qJD(2) * t120 - t158 + t214 - t79) * t80 + t79 * t215 + t80 * (t337 + t341) + (t79 * (t316 + t398) * t241 + (t79 * (-rSges(4,3) - qJ(3)) + t80 * t316) * t239) * qJD(2) + (t35 - g(2)) * t94 + (t34 - g(1)) * (t316 * t239 + t218 + t339)) * m(4) + (t433 - t478) * t413 + (t434 - t436) * t412 + (t442 + t443) * t312 + (((t241 * t308 - t424 + t437) * t241 + (t239 * t308 + t438 - t483 - t98) * t239) * qJD(4) + t446 - t449) * t313 + (m(3) * (t175 ^ 2 + t180 ^ 2) + Icges(4,2) * t245 ^ 2 + (Icges(4,1) * t244 + 0.2e1 * Icges(4,4) * t245) * t244 + Icges(3,3) - t471) * qJDD(2) + (t441 - t444 + t445) * t311; t326 * (g(1) * t239 - g(2) * t241) + 0.2e1 * (t1 * t411 + t2 * t410) * m(6) + 0.2e1 * (t16 * t411 + t17 * t410) * m(5) + 0.2e1 * (t34 * t411 + t35 * t410) * m(4); t461 * t413 + t460 * t412 + (t442 * qJD(2) + t422 * qJD(4) + t433 * qJDD(2) + t437 * t155 + t438 * t156) * t411 + (t441 * qJD(2) + t423 * qJD(4) + t434 * qJDD(2) + t439 * t155 + t440 * t156) * t410 - ((t463 * t238 + t240 * t468) * qJD(4) + (t458 * t238 + t459 * t240) * qJD(2)) * qJD(2) / 0.2e1 + (t444 * t241 + t443 * t239 + (-t436 * t239 - t241 * t478) * qJD(2)) * qJD(2) / 0.2e1 + (-t239 * t478 + t436 * t241) * qJDD(2) / 0.2e1 + t446 * t332 / 0.2e1 + t445 * t331 / 0.2e1 + ((-t430 * t330 - t421) * t239 + ((t431 * t239 + t426) * qJD(4) + t425) * t241) * t313 + ((t438 * t239 + t437 * t241) * qJD(2) + t422) * t312 + ((t440 * t239 + t439 * t241) * qJD(2) + t423) * t311 + ((-t431 * t329 + t421) * t241 + ((t430 * t241 + t426) * qJD(4) + t425) * t239) * t310 + (-g(1) * t427 - g(2) * t428 - g(3) * t450 - (-g(1) * t320 - t408 * t409) * t238 + (-t1 * t348 + t29 * t355 + t3 * t356 + t31 * t404 + (-t30 * t348 + t31 * t357) * qJD(2)) * t241 + (-t2 * t348 + t30 * t355 + t3 * t357 + t31 * t403 + (t29 * t348 - t31 * t356) * qJD(2)) * t239 - (t272 + t394) * qJD(5) - (-t29 * t354 + t30 * t353) * qJD(2) - ((-t29 * t450 + t31 * t353) * t241 + (-t30 * t450 + t31 * t354) * t239) * qJD(4)) * m(6) + (-(t139 * t44 - t397) * qJD(2) - (t46 * (-t139 * t239 - t143 * t241) + t291 * t178) * qJD(4) + t18 * t279 + t46 * ((t117 * t241 - t119 * t239) * qJD(2) + t290) + t291 * t152 + (-t16 * t241 - t17 * t239 + (-t241 * t45 + t395) * qJD(2)) * t173 + g(1) * t143 + g(2) * t139 - g(3) * t178) * m(5); (-(t272 + (t239 ^ 2 + t467) * t394) * qJD(4) + (qJD(4) * t296 + g(3) - t3) * t240 + (qJD(4) * t31 + t2 * t239 + t447 * t241 - t408) * t238) * m(6);];
tau = t4;
