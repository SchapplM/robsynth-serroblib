% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:37
% EndTime: 2019-12-05 16:41:55
% DurationCPUTime: 10.64s
% Computational Cost: add. (11989->406), mult. (9136->503), div. (0->0), fcn. (7094->6), ass. (0->240)
t245 = cos(qJ(4));
t244 = sin(qJ(4));
t392 = Icges(5,4) * t244;
t208 = Icges(5,2) * t245 + t392;
t240 = Icges(6,5) * t244;
t301 = Icges(6,3) * t245 - t240;
t477 = t208 + t301;
t391 = Icges(6,5) * t245;
t210 = Icges(6,1) * t244 - t391;
t241 = Icges(5,4) * t245;
t212 = Icges(5,1) * t244 + t241;
t481 = t210 + t212;
t242 = pkin(8) + qJ(2);
t239 = qJ(3) + t242;
t235 = sin(t239);
t236 = cos(t239);
t303 = Icges(6,1) * t245 + t240;
t124 = -Icges(6,4) * t236 + t235 * t303;
t376 = t235 * t244;
t196 = Icges(5,4) * t376;
t375 = t235 * t245;
t126 = Icges(5,1) * t375 - Icges(5,5) * t236 - t196;
t480 = t124 + t126;
t281 = t303 * t236;
t125 = Icges(6,4) * t235 + t281;
t213 = Icges(5,1) * t245 - t392;
t282 = t213 * t236;
t127 = Icges(5,5) * t235 + t282;
t473 = t125 + t127;
t203 = Icges(6,3) * t244 + t391;
t302 = -Icges(5,2) * t244 + t241;
t479 = t203 - t302;
t478 = t213 + t303;
t468 = -t477 * t244 + t481 * t245;
t372 = t236 * t245;
t195 = Icges(6,5) * t372;
t373 = t236 * t244;
t117 = Icges(6,6) * t235 + Icges(6,3) * t373 + t195;
t205 = Icges(5,5) * t245 - Icges(5,6) * t244;
t278 = t205 * t236;
t119 = Icges(5,3) * t235 + t278;
t207 = Icges(6,4) * t245 + Icges(6,6) * t244;
t279 = t207 * t236;
t121 = Icges(6,2) * t235 + t279;
t475 = t117 * t373 + t473 * t372 + (t119 + t121) * t235;
t120 = -Icges(6,2) * t236 + t207 * t235;
t110 = t235 * t120;
t116 = -Icges(6,6) * t236 + t203 * t235;
t118 = Icges(5,5) * t375 - Icges(5,6) * t376 - Icges(5,3) * t236;
t474 = -t116 * t373 - t235 * t118 - t480 * t372 - t110;
t472 = t479 * qJD(4);
t471 = t478 * qJD(4);
t204 = Icges(5,5) * t244 + Icges(5,6) * t245;
t206 = Icges(6,4) * t244 - Icges(6,6) * t245;
t470 = t204 + t206;
t469 = -t205 - t207;
t243 = qJD(2) + qJD(3);
t467 = ((Icges(6,4) + Icges(5,5)) * t243 - t481 * qJD(4)) * t244 - ((-Icges(5,6) + Icges(6,6)) * t243 + t477 * qJD(4)) * t245;
t122 = Icges(5,4) * t375 - Icges(5,2) * t376 - Icges(5,6) * t236;
t386 = t122 * t244;
t296 = -t126 * t245 + t386;
t387 = t120 * t236;
t299 = t116 * t244 + t124 * t245;
t426 = t235 * t299;
t50 = -t387 + t426;
t466 = -t118 * t236 - t235 * t296 + t50;
t441 = -t122 * t373 - t474;
t280 = t302 * t236;
t123 = Icges(5,6) * t235 + t280;
t440 = -t123 * t373 + t475;
t379 = t206 * t236;
t382 = t204 * t236;
t465 = t468 * t235 - t379 - t382;
t380 = t206 * t235;
t383 = t204 * t235;
t464 = t468 * t236 + t380 + t383;
t310 = -t117 * t376 + t121 * t236 - t125 * t375;
t101 = t127 * t375;
t315 = t119 * t236 - t101;
t53 = -t123 * t376 - t315;
t463 = -t310 + t53;
t218 = rSges(6,1) * t245 + rSges(6,3) * t244;
t462 = pkin(4) * t245 + qJ(5) * t244 + t218;
t448 = rSges(6,1) + pkin(4);
t439 = rSges(6,3) + qJ(5);
t459 = t471 * t245 + t472 * t244 + t470 * t243 + (-t244 * t481 - t477 * t245) * qJD(4);
t458 = (Icges(6,2) + Icges(5,3)) * t243 - t470 * qJD(4);
t385 = t123 * t244;
t457 = -t117 * t244 - t473 * t245 + t385;
t456 = t296 - t299;
t455 = t469 * qJD(4) + t468 * t243;
t454 = t464 * t243;
t374 = t236 * t243;
t229 = t236 * rSges(6,2);
t360 = t462 * t235 - t229;
t453 = t360 * t243;
t452 = ((t278 + t279 + t456) * t243 + t458 * t235) * t236;
t451 = (t440 * t235 - t441 * t236) * qJD(4);
t450 = (t463 * t235 - t466 * t236) * qJD(4);
t449 = t465 * t243;
t447 = t449 + t450;
t446 = t451 + t454;
t445 = t203 * t374 * t245 + t456 * qJD(4) + (-t280 * t245 + (-t281 - t282) * t244) * t243 - t467 * t235;
t377 = t235 * t243;
t444 = -t457 * qJD(4) + (-t478 * t244 + t479 * t245) * t377 + t467 * t236;
t443 = -t455 * t235 + t459 * t236;
t442 = t459 * t235 + t455 * t236;
t438 = (-t117 + t123) * t245 + t473 * t244;
t437 = (-t116 + t122) * t245 + t480 * t244;
t347 = t448 * t244 - t439 * t245;
t316 = t347 * t236;
t232 = t236 * pkin(7);
t167 = pkin(3) * t235 - t232;
t190 = pkin(7) * t374;
t436 = t243 * t167 + t190;
t435 = t235 * rSges(6,2) + pkin(4) * t372;
t432 = t387 + t475;
t341 = qJD(4) * t245;
t328 = t235 * t341;
t371 = t243 * t244;
t431 = t236 * t371 + t328;
t237 = sin(t242);
t400 = pkin(2) * qJD(2);
t336 = t237 * t400;
t165 = rSges(4,1) * t235 + rSges(4,2) * t236;
t384 = t165 * t243;
t132 = -t336 - t384;
t429 = t458 * t236 + t457 * t243 + t469 * t377;
t428 = 0.2e1 * qJD(4);
t339 = qJD(5) * t245;
t359 = rSges(6,1) * t372 + t373 * t439 + t435;
t47 = -t339 + qJD(1) + (t360 * t235 + t359 * t236) * qJD(4);
t427 = qJD(4) * t47;
t168 = t236 * pkin(3) + t235 * pkin(7);
t425 = t168 + t359;
t342 = qJD(4) * t244;
t329 = t235 * t342;
t354 = t448 * t329;
t424 = rSges(5,1) * t372 + t235 * rSges(5,3);
t340 = qJD(5) * t244;
t191 = t236 * t340;
t326 = t236 * t341;
t423 = rSges(6,2) * t374 + t326 * t439 + t191;
t352 = rSges(5,2) * t376 + t236 * rSges(5,3);
t129 = rSges(5,1) * t375 - t352;
t115 = t243 * t129;
t334 = t235 * t371;
t283 = rSges(5,3) * t374 + (-t326 + t334) * rSges(5,2);
t327 = t236 * t342;
t422 = -rSges(5,1) * t327 + t115 + t283 + t436;
t216 = rSges(5,1) * t244 + rSges(5,2) * t245;
t344 = qJD(4) * t235;
t162 = t216 * t344;
t131 = -rSges(5,2) * t373 + t424;
t274 = t131 + t168;
t417 = t243 * t274 - t162;
t414 = t423 + t436 + t453;
t362 = -Icges(5,2) * t375 + t126 - t196;
t366 = t212 * t235 + t122;
t413 = -t244 * t362 - t245 * t366;
t364 = -t301 * t235 + t124;
t368 = -t210 * t235 + t116;
t412 = -t244 * t364 + t245 * t368;
t325 = t235 * t340;
t402 = rSges(6,3) * t328 + t325 + t431 * qJ(5) - t354 + (t218 * t236 + t435) * t243;
t275 = -t243 * t375 - t327;
t403 = t275 * t448 - t439 * t334 + t423;
t5 = (t340 + (t403 + t453) * t236 + (-t359 * t243 + t402) * t235) * qJD(4);
t411 = m(6) * t5;
t407 = t243 / 0.2e1;
t405 = pkin(2) * t237;
t404 = pkin(2) * qJD(2) ^ 2;
t401 = rSges(5,1) * t245;
t343 = qJD(4) * t236;
t330 = t216 * t343;
t276 = -t330 - t336;
t62 = (-t129 - t167) * t243 + t276;
t399 = t236 * t62;
t289 = -qJD(4) * t316 + t191;
t270 = t289 - t336;
t48 = (-t167 - t360) * t243 + t270;
t398 = t243 * t48;
t381 = t205 * t243;
t378 = t207 * t243;
t367 = -Icges(6,1) * t373 + t117 + t195;
t365 = -t212 * t236 - t123;
t363 = -t301 * t236 + t125;
t361 = -t208 * t236 + t127;
t358 = t347 * t235;
t356 = -qJD(4) * t462 + t339;
t351 = -t301 + t303;
t350 = t203 - t210;
t349 = -t208 + t213;
t348 = t212 + t302;
t345 = t229 + t232;
t338 = t237 * t404;
t238 = cos(t242);
t337 = t238 * t404;
t335 = t238 * t400;
t331 = rSges(5,1) * t329 + rSges(5,2) * t431;
t322 = -pkin(3) - t401;
t321 = -t344 / 0.2e1;
t318 = t343 / 0.2e1;
t166 = t236 * rSges(4,1) - rSges(4,2) * t235;
t314 = -t118 + t385;
t311 = t243 * (-pkin(3) * t377 + t190) - t338;
t137 = rSges(4,1) * t374 - rSges(4,2) * t377;
t309 = t339 + t356;
t307 = -rSges(5,2) * t244 + t401;
t63 = t335 + t417;
t306 = -t235 * t63 - t399;
t294 = t129 * t235 + t131 * t236;
t290 = -t344 * t347 + t325;
t273 = -t244 * t363 + t245 * t367;
t272 = -t244 * t361 + t245 * t365;
t271 = t235 * t322 + t232 + t352;
t269 = (t244 * t350 + t245 * t351) * t243;
t268 = (-t244 * t348 + t245 * t349) * t243;
t88 = rSges(5,1) * t275 + t283;
t90 = t243 * t424 - t331;
t261 = (t88 + t115) * t236 + (-t131 * t243 + t90) * t235;
t249 = (t322 * t399 + (t62 * (-rSges(5,3) - pkin(7)) + t63 * t322) * t235) * t243;
t248 = (((t53 - t101 + (t119 + t386) * t236 + t474) * t236 + (t50 - t426 + t432) * t235) * qJD(4) + t454) * t318 + (t468 * qJD(4) + t471 * t244 - t472 * t245) * t243 + (t443 + t444) * t344 / 0.2e1 + (((t236 * t314 - t432 + t440) * t236 + (t235 * t314 - t110 + t310 + t315 + t441) * t235) * qJD(4) + t447 - t449) * t321 - (t442 - t445 + t446) * t343 / 0.2e1 + ((t437 + t465) * t235 + (t438 + t464) * t236) * qJD(4) * t407;
t138 = t168 * t243;
t23 = -t337 + t309 * t343 + (-t138 + (qJD(4) * t347 - t340) * t235 - t402) * t243;
t49 = t243 * t425 + t290 + t335;
t247 = (-t49 * t448 * t342 + (-t244 * t439 - t245 * t448 - pkin(3)) * t398) * t236 + (-t23 * pkin(3) + (-t48 * qJD(5) - t23 * t439) * t244 + (-qJD(4) * t439 * t48 - t23 * t448) * t245 + (t48 * (-rSges(6,2) - pkin(7)) + t49 * (-pkin(3) - t462)) * t243) * t235;
t234 = pkin(2) * t238;
t180 = t307 * qJD(4);
t157 = t216 * t236;
t153 = t216 * t235;
t133 = t166 * t243 + t335;
t113 = -t137 * t243 - t337;
t112 = -t243 * t384 - t338;
t64 = qJD(4) * t294 + qJD(1);
t46 = -t337 - t180 * t343 + (-t138 - t90 + t162) * t243;
t45 = t243 * t88 + (-t180 * t235 - t216 * t374) * qJD(4) + t311;
t25 = t261 * qJD(4);
t22 = (t191 + t403) * t243 + (t235 * t309 - t243 * t316) * qJD(4) + t311;
t1 = [m(5) * t25 + t411; m(4) * (t113 * (-t165 - t405) + t112 * (t166 + t234) + (-t137 - t335 + t133) * t132) + t248 + (t23 * (t345 - t405) + t48 * (-t335 + t354) + t22 * (t234 + t425) + t247 + (t48 - t270 - t336 + t414) * t49) * m(6) + (t46 * (t271 - t405) + t62 * (t331 - t335) + t45 * (t234 + t274) + t249 + (-t276 + t62 - t336 + t422) * t63) * m(5); t248 + (t23 * t345 + t247 + (-t289 + t414) * t49 + (t290 + t354) * t48 + (t22 + t398) * t425) * m(6) + (t46 * t271 + t45 * t274 + t249 + (t330 + t422) * t63 + (t331 + t417) * t62) * m(5) + (t112 * t166 - t113 * t165 - t132 * t137 - t133 * t384 - (-t132 * t166 - t133 * t165) * t243) * m(4); -(((t348 - t350) * t245 + (t349 + t351) * t244) * t243 + (((-t362 - t364) * t236 + (t361 + t363) * t235) * t245 + ((t366 - t368) * t236 + (t365 + t367) * t235) * t244) * qJD(4)) * t243 / 0.2e1 + ((t243 * t438 + t445) * t236 + (t243 * t437 + t444) * t235) * t407 + ((-t344 * t379 + t378) * t235 + (t269 + (-t412 * t236 + (t380 + t273) * t235) * qJD(4)) * t236 + (-t344 * t382 + t381) * t235 + (t268 + (-t413 * t236 + (t383 + t272) * t235) * qJD(4)) * t236) * t321 + ((-t343 * t380 - t378) * t236 + (t269 + (t273 * t235 + (t379 - t412) * t236) * qJD(4)) * t235 + (-t343 * t383 - t381) * t236 + (t268 + (t272 * t235 + (t382 - t413) * t236) * qJD(4)) * t235) * t318 + (-(t244 * t47 + (t235 * t49 + t236 * t48) * t245) * qJD(5) - (-t316 * t49 + t358 * t48) * t243 - ((-t316 * t47 - t462 * t48) * t236 + (-t358 * t47 - t462 * t49) * t235) * qJD(4) + (-t23 * t347 + t48 * t356 + t5 * t359 + t47 * t403 + (-t347 * t49 + t360 * t47) * t243) * t236 + (-t22 * t347 + t49 * t356 + t5 * t360 + t47 * t402 + (t347 * t48 - t359 * t47) * t243) * t235) * m(6) + (t25 * t294 + t64 * t261 + t306 * t180 + ((-t243 * t63 - t46) * t236 + (t243 * t62 - t45) * t235) * t216 - (t153 * t62 - t157 * t63) * t243 - (t64 * (-t153 * t235 - t157 * t236) + t306 * t307) * qJD(4)) * m(5) + (t443 * t243 + (t440 * t374 + (t429 * t235 + t441 * t243 - t452) * t235) * t428) * t235 / 0.2e1 - (t442 * t243 + ((t243 * t463 + t452) * t236 + (-t429 * t236 + t243 * t466) * t235) * t428) * t236 / 0.2e1 + (t447 + t450) * t377 / 0.2e1 + (t446 + t451) * t374 / 0.2e1; -t245 * t411 + 0.2e1 * (m(6) * (t22 * t235 + t23 * t236 + t427) / 0.2e1 - m(6) * (t235 ^ 2 + t236 ^ 2) * t427 / 0.2e1) * t244;];
tauc = t1(:);
