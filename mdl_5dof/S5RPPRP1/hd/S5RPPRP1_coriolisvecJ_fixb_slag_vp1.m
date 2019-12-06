% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:36
% EndTime: 2019-12-05 17:36:20
% DurationCPUTime: 22.41s
% Computational Cost: add. (11133->464), mult. (15000->629), div. (0->0), fcn. (14318->8), ass. (0->232)
t498 = Icges(5,4) + Icges(6,4);
t499 = Icges(5,1) + Icges(6,1);
t465 = Icges(5,5) + Icges(6,5);
t464 = Icges(5,2) + Icges(6,2);
t463 = Icges(5,6) + Icges(6,6);
t500 = Icges(5,3) + Icges(6,3);
t235 = qJ(1) + pkin(7);
t231 = sin(t235);
t232 = cos(t235);
t241 = cos(qJ(4));
t237 = cos(pkin(8));
t239 = sin(qJ(4));
t374 = t237 * t239;
t180 = t231 * t374 + t232 * t241;
t510 = t498 * t180;
t182 = -t231 * t241 + t232 * t374;
t509 = t498 * t182;
t373 = t237 * t241;
t379 = t231 * t239;
t183 = t232 * t373 + t379;
t508 = t498 * t183;
t507 = t498 * t241;
t506 = t498 * t239;
t376 = t232 * t239;
t181 = t231 * t373 - t376;
t505 = t498 * t181;
t236 = sin(pkin(8));
t378 = t232 * t236;
t434 = t499 * t183 + t465 * t378 - t509;
t435 = -t464 * t182 + t463 * t378 + t508;
t491 = t435 * t180 - t434 * t181;
t493 = t182 * t435 - t183 * t434;
t380 = t231 * t236;
t462 = -t464 * t180 + t463 * t380 + t505;
t483 = -t499 * t181 - t465 * t380 + t510;
t467 = t463 * t180 - t465 * t181 - t500 * t380;
t486 = -t463 * t182 + t465 * t183 + t500 * t378;
t496 = (-t465 * t239 - t463 * t241) * t236;
t495 = (-t464 * t241 - t506) * t236;
t494 = (-t499 * t239 - t507) * t236;
t442 = -(-t239 * t435 + t241 * t434) * t236 + t486 * t237;
t460 = t500 * t237 + (t463 * t239 - t465 * t241) * t236;
t459 = -t463 * t237 + (-t464 * t239 + t507) * t236;
t458 = -t465 * t237 + (t499 * t241 - t506) * t236;
t492 = -t462 * t180 - t483 * t181;
t489 = -t486 * t380 + t491;
t477 = -t462 * t182 - t483 * t183 - t467 * t378;
t488 = t486 * t378 - t493;
t125 = qJD(1) * t180 - qJD(4) * t183;
t258 = t182 * qJD(4);
t126 = -qJD(1) * t181 - t258;
t346 = qJD(1) * t236;
t333 = t231 * t346;
t487 = t463 * t125 + t465 * t126 - t333 * t500;
t127 = qJD(1) * t182 + qJD(4) * t181;
t128 = -qJD(1) * t183 + qJD(4) * t180;
t332 = t232 * t346;
t472 = t463 * t127 + t465 * t128 - t332 * t500;
t471 = t464 * t125 + t498 * t126 - t463 * t333;
t470 = t464 * t127 + t498 * t128 - t463 * t332;
t469 = t498 * t125 + t499 * t126 - t465 * t333;
t468 = t498 * t127 + t499 * t128 - t465 * t332;
t457 = t496 * qJD(4);
t456 = t495 * qJD(4);
t455 = t494 * qJD(4);
t482 = t467 * t380 - t492;
t481 = t459 * t180 - t458 * t181 + t460 * t380;
t439 = t459 * t182 - t458 * t183 + t460 * t378;
t480 = t487 * t232;
t479 = t472 * t237 + (-t468 * t241 + t470 * t239 + (t483 * t239 - t462 * t241) * qJD(4)) * t236;
t478 = -t487 * t237 + (t469 * t241 - t471 * t239 + (-t434 * t239 - t435 * t241) * qJD(4)) * t236;
t475 = t467 * t237 + (-t462 * t239 - t483 * t241) * t236;
t474 = -t457 * t237 + (t455 * t241 - t456 * t239 + (-t458 * t239 - t459 * t241) * qJD(4)) * t236;
t454 = (t477 * t231 + t488 * t232) * t236;
t453 = (t482 * t231 + t489 * t232) * t236;
t452 = -t489 * t231 + t482 * t232;
t225 = -qJD(4) * t237 + qJD(1);
t451 = t481 * t225;
t238 = -qJ(5) - pkin(6);
t407 = pkin(6) + t238;
t229 = pkin(4) * t241 + pkin(3);
t408 = pkin(3) - t229;
t359 = (t407 - rSges(6,3)) * t237 + (rSges(6,1) * t241 - rSges(6,2) * t239 - t408) * t236;
t219 = pkin(4) * t379;
t377 = t232 * t237;
t450 = -rSges(6,1) * t183 + rSges(6,2) * t182 - t229 * t377 - t219;
t449 = t439 * t225;
t343 = qJD(5) * t236;
t329 = t231 * t343;
t345 = qJD(4) * t236;
t330 = t232 * t345;
t324 = t407 * t236;
t411 = pkin(3) * t237;
t364 = (t324 + t411) * t232 - rSges(6,3) * t378 + t450;
t448 = t225 * t364 + t359 * t330 - t329;
t447 = rSges(6,1) + pkin(4);
t446 = qJD(4) * t453 + t451;
t445 = qJD(4) * t454 - t449;
t348 = qJD(1) * t231;
t444 = (t232 * t457 + t348 * t460) * t236 + t455 * t183 - t456 * t182 + t458 * t126 + t459 * t125;
t347 = qJD(1) * t232;
t443 = (-t231 * t457 + t347 * t460) * t236 - t455 * t181 + t456 * t180 + t458 * t128 + t459 * t127;
t441 = t474 * t225;
t437 = qJ(3) * t232;
t432 = (t478 * t232 + t479 * t231 + (t231 * t442 + t232 * t475) * qJD(1)) * t236;
t431 = (t452 * qJD(1) + (-t236 * t347 * t486 + t127 * t435 + t128 * t434 + t180 * t471 - t181 * t469) * t232 + ((t231 * t472 + t347 * t467 - t480) * t236 + t468 * t181 - t470 * t180 - t483 * t128 + t462 * t127) * t231) * t236;
t430 = (((-t348 * t486 + t480) * t236 + t469 * t183 - t471 * t182 + t434 * t126 + t435 * t125 + t477 * qJD(1)) * t232 + ((-t232 * t472 + t348 * t467) * t236 - t468 * t183 + t470 * t182 - t483 * t126 + t462 * t125 - t488 * qJD(1)) * t231) * t236;
t429 = (t183 * t464 - t434 + t509) * t232 + (t181 * t464 + t483 + t510) * t231;
t428 = (-t182 * t499 - t435 - t508) * t232 + (-t180 * t499 - t462 - t505) * t231;
t427 = (t182 * t465 + t183 * t463) * t232 + (t180 * t465 + t181 * t463) * t231;
t426 = t458 + t495;
t425 = t459 - t494;
t298 = -rSges(5,1) * t183 + rSges(5,2) * t182;
t118 = rSges(5,3) * t378 - t298;
t179 = -rSges(5,3) * t237 + (rSges(5,1) * t241 - rSges(5,2) * t239) * t236;
t424 = -t118 * t225 + t179 * t330;
t390 = -rSges(6,3) + t238;
t341 = pkin(6) + t390;
t423 = t236 * t341;
t409 = pkin(6) * t236;
t305 = t409 + t411;
t184 = t305 * t232;
t325 = t408 * t237;
t215 = t232 * t343;
t422 = -rSges(6,1) * t126 - rSges(6,2) * t125 - t215;
t242 = cos(qJ(1));
t413 = pkin(1) * t242;
t304 = -rSges(3,1) * t232 - t413;
t421 = t231 * rSges(3,2) + t304;
t296 = -t181 * rSges(6,1) + t180 * rSges(6,2);
t257 = pkin(4) * t376 + t296;
t396 = pkin(4) * qJD(4);
t336 = t239 * t396;
t313 = t237 * t336;
t335 = t241 * t396;
t420 = t128 * rSges(6,1) + t127 * rSges(6,2) + t231 * t313 + t232 * t335 + t238 * t332;
t243 = qJD(1) ^ 2;
t240 = sin(qJ(1));
t414 = pkin(1) * t240;
t412 = pkin(1) * t243;
t410 = pkin(4) * t239;
t406 = -rSges(6,3) * t332 - t329 + (-t219 + (t325 + t409) * t232) * qJD(1) + t420;
t157 = -pkin(6) * t333 - t348 * t411;
t375 = t236 * t238;
t381 = t229 * t237;
t283 = t375 - t381;
t405 = rSges(6,3) * t333 - t283 * t348 - (t239 * t347 - t258) * pkin(4) + t157 + t422;
t399 = rSges(4,1) * t237;
t398 = rSges(6,3) * t236;
t397 = pkin(1) * qJD(1);
t391 = rSges(4,3) + qJ(3);
t383 = qJ(3) * t231;
t362 = t128 * rSges(5,1) + t127 * rSges(5,2);
t361 = rSges(6,2) * t181 + t180 * t447;
t360 = rSges(6,2) * t183 + t182 * t447;
t226 = qJD(3) * t231;
t349 = -pkin(2) * t348 + t226;
t156 = qJ(3) * t347 + t349;
t358 = -t156 - t157;
t357 = -t181 * rSges(5,1) + t180 * rSges(5,2);
t268 = t236 * (-rSges(6,1) * t239 - rSges(6,2) * t241);
t342 = qJD(5) * t237;
t352 = qJD(4) * t268 - t236 * t336 - t342;
t206 = pkin(2) * t232 + t383;
t227 = qJD(3) * t232;
t351 = qJD(1) * t206 - t227;
t339 = t242 * t412;
t338 = rSges(5,3) * t380;
t337 = t240 * t397;
t331 = t231 * t345;
t326 = -pkin(2) - t399;
t322 = -t345 / 0.2e1;
t321 = t345 / 0.2e1;
t320 = -t206 - t413;
t318 = -pkin(2) - t381;
t317 = qJ(3) + t410;
t316 = t359 * t232;
t315 = t352 * t232;
t310 = t231 * t322;
t309 = t231 * t321;
t308 = t232 * t322;
t307 = t232 * t321;
t306 = qJD(1) * t322;
t301 = rSges(3,1) * t231 + rSges(3,2) * t232;
t300 = -rSges(4,2) * t236 + t399;
t299 = rSges(5,1) * t126 + rSges(5,2) * t125;
t269 = t236 * (-rSges(5,1) * t239 - rSges(5,2) * t241);
t192 = qJD(4) * t269;
t230 = t240 * t412;
t70 = -rSges(5,3) * t333 + t299;
t29 = t192 * t330 - t225 * t70 + t230 + ((-t179 * t345 - qJD(3)) * t231 + t358) * qJD(1);
t280 = -t339 + (t227 - t351) * qJD(1);
t270 = -t184 * t243 + t280;
t72 = -rSges(5,3) * t332 + t362;
t30 = t225 * t72 + (t179 * t347 + t192 * t231) * t345 + t270;
t293 = t30 * t231 + t29 * t232;
t116 = -t338 + t357;
t276 = -qJD(1) * (-pkin(2) * t231 + t437) - t226 + t337;
t263 = t305 * t348 + t276;
t44 = -t225 * t116 - t179 * t331 + t263;
t260 = t227 + (-t184 + t320) * qJD(1);
t45 = t260 + t424;
t290 = -t231 * t44 + t232 * t45;
t289 = -t231 * t70 - t232 * t72;
t286 = -t116 * t232 - t118 * t231;
t278 = t231 * t306;
t277 = t232 * t306;
t275 = -rSges(4,1) * t377 - rSges(4,3) * t231;
t262 = qJD(1) * t184 + t242 * t397 + t351;
t259 = -rSges(5,3) * t236 - pkin(2) - t305;
t252 = qJD(1) * t316 + t231 * t352;
t245 = t236 * (-t257 * t232 + ((-t325 - t423) * t232 + t364) * t231);
t244 = t236 * ((qJD(1) * t364 - t406) * t232 + (((-t381 + t411 + t423) * t231 + t257) * qJD(1) + t405) * t231);
t223 = rSges(3,2) * t348;
t216 = rSges(4,2) * t378;
t211 = rSges(4,2) * t332;
t152 = -t216 - t275;
t144 = -rSges(5,1) * t182 - rSges(5,2) * t183;
t142 = rSges(5,1) * t180 + rSges(5,2) * t181;
t96 = t227 + (-t152 + t320) * qJD(1);
t95 = -qJD(1) * (rSges(4,3) * t232 - t300 * t231) + t276;
t88 = (qJD(1) * t275 + t211) * qJD(1) + t280;
t87 = t230 + (-rSges(4,3) * t347 - t156 + (qJD(1) * t300 - qJD(3)) * t231) * qJD(1);
t48 = t286 * t345 + qJD(2);
t28 = t260 + t448;
t27 = -t215 + t263 - t359 * t331 + (-(t325 + t324) * t231 + rSges(6,3) * t380 - t257) * t225;
t25 = qJD(4) * t245 + qJD(2) - t342;
t20 = ((t116 * t231 - t118 * t232) * qJD(1) + t289) * t345;
t19 = t406 * t225 + (qJD(4) * t252 - qJD(5) * t348) * t236 + t270;
t18 = t230 + t405 * t225 + t315 * t345 + (-t215 + (-t345 * t359 - qJD(3)) * t231 + t358) * qJD(1);
t1 = qJD(4) * t244;
t2 = [m(3) * ((t243 * t301 + t230) * t421 - (qJD(1) * t301 + t337) * (qJD(1) * t304 + t223) + (t339 + (rSges(3,1) * t347 + qJD(1) * t421 - t223) * qJD(1)) * (t301 + t414)) + ((t491 * t232 + (t482 - t488 - t493) * t231) * t345 + t451) * t308 + (-(t28 + t262 - t448) * t27 + t18 * (-t413 + t450) + t28 * (-t349 + t422) + t19 * (t296 - t414) - t27 * (t227 + t420) + (t18 * (-pkin(2) + t375 - t398) + t28 * t313 + t19 * t317) * t232 + (-t18 * qJ(3) - t28 * t335 + t19 * t318 + (t27 * qJD(5) + t19 * t390) * t236) * t231 + ((t240 * t28 + t242 * t27) * pkin(1) + (-t28 * t317 - t27 * (t318 - t398)) * t232 + (t28 * (-t283 + t398) + t27 * t317) * t231) * qJD(1)) * m(6) + (t29 * (t298 - t413) + t30 * (t357 - t414) + (t30 * qJ(3) + t29 * t259) * t232 + (-t29 * qJ(3) + t259 * t30) * t231 + (-t299 - t349 - t157 + (t338 + t414 - t437) * qJD(1)) * t45 + (-t262 - t45 - t227 - t362 + (-t232 * t259 + t383 + t413) * qJD(1) + t424) * t44) * m(5) + (-(t96 + (t152 + t413) * qJD(1) + t351) * t95 + t87 * (t216 - t413) - t96 * t349 - t88 * t414 - t95 * (t211 + t227) + (t326 * t87 + t391 * t88) * t232 + (-t87 * t391 + t88 * (-pkin(2) - t300)) * t231 + ((t240 * t96 + t242 * t95) * pkin(1) + (-t326 * t95 - t391 * t96) * t232 + (t300 * t96 + t391 * t95) * t231) * qJD(1)) * m(4) + (((t492 + t493) * t232 + (-t477 + t491) * t231 + (-t486 * t231 ^ 2 + (-t231 * t467 - t232 * t486) * t232) * t236 + t452) * t345 + t445 + t449) * t309 + (-t439 - t442) * t278 + (t481 - t475) * t277 + (t444 + t446 + t478) * t307 + t441 + (t443 - t479) * t310; m(5) * t20 + m(6) * t1; m(4) * (t231 * t88 + t232 * t87) + m(5) * t293 + m(6) * (t18 * t232 + t19 * t231); -(((-t239 * t426 - t241 * t425) * t225 + ((t239 * t429 + t241 * t428) * t236 + t427 * t237) * qJD(4)) * t236 - t496 * t225 * t237) * t225 / 0.2e1 + (-t237 * t474 + t432) * t225 / 0.2e1 - (t432 * qJD(4) + t441) * t237 / 0.2e1 - (t431 * qJD(4) + t443 * t225) * t380 / 0.2e1 + (t430 * qJD(4) + t444 * t225) * t378 / 0.2e1 + (-t237 * t443 + t431) * t310 + ((-t180 * t429 - t181 * t428 + t380 * t427) * t345 + (t180 * t426 + t181 * t425 - t380 * t496) * t225) * t309 + ((t182 * t429 + t183 * t428 - t378 * t427) * t345 + (-t182 * t426 - t183 * t425 + t378 * t496) * t225) * t308 + (-t237 * t444 + t430) * t307 + (t237 * t439 + t454) * t278 + (-t237 * t481 + t453) * t277 + (-(-t27 * t361 + t28 * t360) * t225 - (t25 * (t231 * t360 - t232 * t361) + (t236 * t410 - t268) * (t231 * t27 - t232 * t28)) * t345 + t1 * t245 + t25 * t244 + t18 * (t236 * t316 - t237 * t364) + t28 * (-t405 * t237 + (-t348 * t359 + t315) * t236) + t19 * (-t257 * t237 + (-t237 ^ 2 * t408 + (-t237 * t341 + t359) * t236) * t231) - t27 * (t236 * t252 - t237 * t406)) * m(6) + ((-t116 * t30 + t118 * t29 + t44 * t72 + t45 * t70) * t237 + (t20 * t286 + t48 * (t116 * t348 - t118 * t347 + t289) + t290 * t192 + ((-t231 * t45 - t232 * t44) * qJD(1) + t293) * t179) * t236 - (-t142 * t44 - t144 * t45) * t225 - (t48 * (-t142 * t232 - t144 * t231) + t290 * t269) * t345) * m(5) - (t231 * t445 + t232 * t446) * t346 / 0.2e1; (-t1 * t237 + (-t18 * t231 + t19 * t232) * t236) * m(6);];
tauc = t2(:);
