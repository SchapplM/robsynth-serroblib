% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:04
% EndTime: 2022-01-20 10:48:16
% DurationCPUTime: 4.89s
% Computational Cost: add. (42777->342), mult. (23906->454), div. (0->0), fcn. (21722->10), ass. (0->225)
t322 = qJ(4) + qJ(5);
t319 = cos(t322);
t311 = Icges(6,4) * t319;
t317 = sin(t322);
t270 = -Icges(6,2) * t317 + t311;
t271 = Icges(6,1) * t317 + t311;
t497 = t270 + t271;
t326 = cos(qJ(4));
t321 = Icges(5,4) * t326;
t324 = sin(qJ(4));
t290 = -Icges(5,2) * t324 + t321;
t291 = Icges(5,1) * t324 + t321;
t496 = t290 + t291;
t323 = qJ(1) + qJ(2);
t316 = pkin(9) + t323;
t314 = cos(t316);
t308 = t314 * pkin(7);
t313 = sin(t316);
t433 = rSges(5,1) * t326;
t365 = pkin(3) + t433;
t403 = t313 * t324;
t368 = rSges(5,2) * t403 + t314 * rSges(5,3);
t318 = sin(t323);
t440 = pkin(2) * t318;
t181 = -t365 * t313 + t308 + t368 - t440;
t398 = t314 * t324;
t287 = rSges(5,2) * t398;
t320 = cos(t323);
t439 = pkin(2) * t320;
t182 = t439 - t287 + t365 * t314 + (rSges(5,3) + pkin(7)) * t313;
t293 = rSges(5,1) * t324 + rSges(5,2) * t326;
t248 = t293 * t313;
t249 = t293 * t314;
t102 = t181 * t248 - t182 * t249;
t404 = t313 * t319;
t405 = t313 * t317;
t210 = rSges(6,1) * t404 - rSges(6,2) * t405 - t314 * rSges(6,3);
t437 = pkin(4) * t326;
t315 = pkin(3) + t437;
t477 = -pkin(8) - pkin(7);
t369 = -t313 * t315 - t314 * t477;
t172 = -t210 + t369 - t440;
t297 = t313 * t477;
t400 = t314 * t317;
t362 = -rSges(6,2) * t400 + t313 * rSges(6,3);
t432 = rSges(6,1) * t319;
t173 = t439 - t297 + (t315 + t432) * t314 + t362;
t273 = rSges(6,1) * t317 + rSges(6,2) * t319;
t438 = pkin(4) * t324;
t342 = t273 + t438;
t488 = t342 * t314;
t489 = t342 * t313;
t91 = t172 * t489 - t173 * t488;
t495 = m(5) * t102 + m(6) * t91;
t479 = m(5) / 0.2e1;
t478 = m(6) / 0.2e1;
t452 = t313 / 0.2e1;
t451 = -t314 / 0.2e1;
t494 = t314 / 0.2e1;
t441 = cos(qJ(1)) * pkin(1);
t442 = sin(qJ(1)) * pkin(1);
t449 = m(3) * (t441 * (-rSges(3,1) * t318 - rSges(3,2) * t320) + (t320 * rSges(3,1) - t318 * rSges(3,2)) * t442);
t447 = m(4) * (t441 * (-rSges(4,1) * t313 - rSges(4,2) * t314 - t440) + (t314 * rSges(4,1) - t313 * rSges(4,2) + t439) * t442);
t309 = t313 ^ 2;
t310 = t314 ^ 2;
t367 = t309 + t310;
t453 = -t313 / 0.2e1;
t491 = t452 + t453;
t166 = t172 - t442;
t167 = t173 + t441;
t90 = t166 * t489 - t167 * t488;
t238 = t273 * t313;
t239 = t273 * t314;
t168 = -t313 * t238 - t314 * t239;
t275 = -rSges(6,2) * t317 + t432;
t351 = Icges(6,5) * t317 + Icges(6,6) * t319;
t232 = t351 * t313;
t233 = t314 * t351;
t422 = Icges(6,4) * t317;
t272 = Icges(6,1) * t319 - t422;
t209 = Icges(6,5) * t313 + t272 * t314;
t269 = Icges(6,2) * t319 + t422;
t374 = -t269 * t314 + t209;
t207 = Icges(6,6) * t313 + t270 * t314;
t376 = -t271 * t314 - t207;
t337 = -t374 * t317 + t376 * t319;
t279 = Icges(6,4) * t405;
t208 = Icges(6,1) * t404 - Icges(6,5) * t314 - t279;
t375 = -Icges(6,2) * t404 + t208 - t279;
t206 = Icges(6,4) * t404 - Icges(6,2) * t405 - Icges(6,6) * t314;
t377 = t271 * t313 + t206;
t338 = t375 * t317 + t377 * t319;
t436 = (-t309 * t233 + (t338 * t314 + (t232 + t337) * t313) * t314) * t452 + (-t310 * t232 + (t337 * t313 + (t233 + t338) * t314) * t313) * t451;
t399 = t314 * t319;
t149 = t313 * t210 + t314 * (rSges(6,1) * t399 + t362);
t95 = -t313 * (pkin(3) * t313 - t308 + t369) + (-t313 * pkin(7) - t297 + (-pkin(3) + t315) * t314) * t314 + t149;
t15 = t436 + m(6) * (t95 * t168 + (t313 * t489 + t314 * t488) * t275);
t490 = t15 * qJD(5);
t75 = -t173 * t166 + t167 * t172;
t178 = t181 - t442;
t179 = t182 + t441;
t88 = -t182 * t178 + t179 * t181;
t487 = qJD(1) + qJD(2);
t435 = (t90 - t91) * t478 + ((-t179 + t182) * t314 + (t178 - t181) * t313) * t293 * t479;
t99 = t178 * t248 - t179 * t249;
t486 = (t91 + t90) * t478 + (t102 + t99) * t479;
t423 = Icges(5,4) * t324;
t289 = Icges(5,2) * t326 + t423;
t292 = Icges(5,1) * t326 - t423;
t485 = t496 * t326 / 0.2e1 + (t292 / 0.2e1 - t289 / 0.2e1) * t324;
t357 = t497 * t319 / 0.2e1 + (-t269 / 0.2e1 + t272 / 0.2e1) * t317;
t187 = t209 * t404;
t268 = Icges(6,5) * t319 - Icges(6,6) * t317;
t411 = t268 * t314;
t205 = Icges(6,3) * t313 + t411;
t361 = t205 * t314 - t187;
t107 = -t207 * t405 - t361;
t204 = Icges(6,5) * t404 - Icges(6,6) * t405 - Icges(6,3) * t314;
t382 = -t313 * t204 - t208 * t399;
t108 = -t206 * t400 - t382;
t381 = t313 * t205 + t209 * t399;
t109 = -t207 * t400 + t381;
t359 = t207 * t317 - t204;
t415 = t206 * t317;
t363 = ((t107 - t187 + (t205 + t415) * t314 + t382) * t314 + t381 * t313) * t451 + (-t108 * t314 + t109 * t313) * t494 + ((t359 * t313 + t107 + t108 + t361) * t313 + (t109 - t381 + (-t208 * t319 + t415) * t313 + (t359 + t204) * t314) * t314) * t452;
t484 = 4 * qJD(1);
t482 = 2 * qJD(4);
t474 = m(5) * t88;
t472 = m(5) * t99;
t92 = t166 * t238 - t167 * t239;
t93 = t172 * t238 - t173 * t239;
t466 = m(6) * (t93 + t92);
t464 = m(6) * ((-t167 + t173) * t314 + (t166 - t172) * t313) * t273;
t346 = (-t166 * t314 - t167 * t313) * t275;
t383 = -t238 * t488 + t239 * t489;
t463 = m(6) * (t346 + t383);
t344 = (t313 * t488 - t314 * t489) * t273;
t462 = m(6) * (t346 + t344);
t345 = (-t172 * t314 - t173 * t313) * t275;
t461 = m(6) * (t345 + t383);
t460 = m(6) * (t345 + t344);
t459 = m(6) * t75;
t457 = m(6) * t90;
t455 = m(6) * t92;
t454 = m(6) * t93;
t402 = t313 * t326;
t223 = Icges(5,4) * t402 - Icges(5,2) * t403 - Icges(5,6) * t314;
t221 = Icges(5,5) * t402 - Icges(5,6) * t403 - Icges(5,3) * t314;
t285 = Icges(5,4) * t403;
t225 = Icges(5,1) * t402 - Icges(5,5) * t314 - t285;
t397 = t314 * t326;
t379 = -t313 * t221 - t225 * t397;
t120 = -t223 * t398 - t379;
t224 = Icges(5,6) * t313 + t290 * t314;
t288 = Icges(5,5) * t326 - Icges(5,6) * t324;
t408 = t288 * t314;
t222 = Icges(5,3) * t313 + t408;
t226 = Icges(5,5) * t313 + t292 * t314;
t378 = t313 * t222 + t226 * t397;
t121 = -t224 * t398 + t378;
t358 = t224 * t324 - t221;
t196 = t226 * t402;
t360 = t222 * t314 - t196;
t28 = (t358 * t314 + t121 - t378) * t314 + (t358 * t313 + t120 + t360) * t313;
t119 = -t224 * t403 - t360;
t413 = t223 * t324;
t29 = (t119 - t196 + (t222 + t413) * t314 + t379) * t314 + t378 * t313;
t84 = -(-(-t225 * t326 + t413) * t313 - t221 * t314) * t314 + t119 * t313;
t85 = -t120 * t314 + t121 * t313;
t2 = (t85 / 0.2e1 - t29 / 0.2e1) * t314 + (t28 / 0.2e1 + t84 / 0.2e1) * t313 + t363;
t450 = t2 * qJD(4);
t443 = m(6) * t168;
t427 = qJD(5) * t363;
t373 = t291 * t313 + t223;
t372 = -t291 * t314 - t224;
t371 = -Icges(5,2) * t402 + t225 - t285;
t370 = -t289 * t314 + t226;
t364 = -t275 - t437;
t355 = t466 / 0.2e1 + t357;
t352 = Icges(5,5) * t324 + Icges(5,6) * t326;
t343 = (-t248 * t314 + t249 * t313) * t293;
t334 = (-t269 + t272) * t319 - t497 * t317;
t341 = -t363 + (t268 * t313 + t334 * t314 + t376 * t317 + t374 * t319) * t452 + (t334 * t313 - t377 * t317 + t375 * t319 - t411) * t451;
t340 = -t357 + t491 * (t319 * t206 + t317 * t208);
t339 = t357 + t485;
t336 = t371 * t324 + t373 * t326;
t335 = -t370 * t324 + t372 * t326;
t333 = (-t289 + t292) * t326 - t496 * t324;
t331 = t339 + t486;
t330 = (-t238 * t314 + t239 * t313) * t273;
t329 = t340 - t485 + t491 * (t326 * t223 + t324 * t225);
t328 = (t29 * t494 + t341 + (t333 * t313 - t373 * t324 + t371 * t326 - t408 + t85) * t451 + (t28 + t84) * t453 + (t288 * t313 + t333 * t314 + t372 * t324 + t370 * t326) * t452) * qJD(4);
t295 = -rSges(5,2) * t324 + t433;
t243 = t314 * t352;
t242 = t352 * t313;
t220 = t364 * t314;
t218 = t364 * t313;
t174 = -t248 * t313 - t249 * t314;
t162 = qJD(5) * t443;
t148 = -t367 * t438 + t168;
t74 = t357 + t454;
t73 = t357 + t455;
t67 = t460 / 0.2e1;
t65 = t461 / 0.2e1;
t58 = t462 / 0.2e1;
t57 = t463 / 0.2e1;
t53 = t464 / 0.2e1;
t37 = t339 + t495;
t36 = t339 + t457 + t472;
t30 = t447 + t449 + t459 + t474;
t19 = -t464 / 0.2e1 + t355;
t18 = t53 + t355;
t17 = m(6) * (t367 * t273 * t275 + t149 * t168) + t436;
t16 = t17 * qJD(5);
t14 = t53 - t466 / 0.2e1 + t340;
t13 = t331 - t435;
t12 = t331 + t435;
t9 = t329 + t435 - t486;
t8 = t65 - t460 / 0.2e1 + t363;
t7 = t67 - t461 / 0.2e1 + t363;
t6 = t57 - t462 / 0.2e1 + t363;
t5 = t58 - t463 / 0.2e1 + t363;
t4 = t65 + t67 + t341;
t3 = t57 + t58 + t341;
t1 = [qJD(2) * t30 + qJD(4) * t36 + qJD(5) * t73, t30 * qJD(1) + t12 * qJD(4) + t18 * qJD(5) + 0.2e1 * (t447 / 0.2e1 + t449 / 0.2e1 + t75 * t478 + t88 * t479) * qJD(2), 0, t36 * qJD(1) + t12 * qJD(2) + t3 * qJD(5) + ((t166 * t220 + t167 * t218) * t478 + ((-t178 * t314 - t179 * t313) * t295 + t343) * t479) * t482 + t328, t73 * qJD(1) + t18 * qJD(2) + t3 * qJD(4) + ((t346 + t330) * m(6) + t341) * qJD(5); t13 * qJD(4) + t19 * qJD(5) + (-t459 / 0.4e1 - t474 / 0.4e1 - t447 / 0.4e1 - t449 / 0.4e1) * t484, qJD(4) * t37 + qJD(5) * t74, 0, t13 * qJD(1) + t37 * qJD(2) + t4 * qJD(5) + ((t172 * t220 + t173 * t218) * t478 + ((-t181 * t314 - t182 * t313) * t295 + t343) * t479) * t482 + t328, t19 * qJD(1) + t74 * qJD(2) + t4 * qJD(4) + ((t345 + t330) * m(6) + t341) * qJD(5); 0, 0, 0, (t148 * t478 + t174 * t479) * t482 + t162, qJD(4) * t443 + t162; t329 * qJD(1) + t9 * qJD(2) + t6 * qJD(5) + (-t457 / 0.4e1 - t472 / 0.4e1) * t484 + t450, t9 * qJD(1) + t8 * qJD(5) + t450 + (t329 - t495) * qJD(2), 0, (m(5) * (t367 * t295 * t293 + (t313 * (rSges(5,1) * t402 - t368) + t314 * (rSges(5,1) * t397 + t313 * rSges(5,3) - t287)) * t174) + (-t309 * t243 + (t336 * t314 + (t242 + t335) * t313) * t314) * t452 + (-t310 * t242 + (t335 * t313 + (t243 + t336) * t314) * t313) * t451 + m(6) * (t148 * t95 - t218 * t489 - t220 * t488) + t436) * qJD(4) + t490 + t487 * t2, t6 * qJD(1) + t8 * qJD(2) + t15 * qJD(4) + t490; (t340 - t455) * qJD(1) + t14 * qJD(2) + t5 * qJD(4) + t427, t14 * qJD(1) + (t340 - t454) * qJD(2) + t7 * qJD(4) + t427, 0, t5 * qJD(1) + t7 * qJD(2) + ((t148 * t149 + (-t218 * t313 - t220 * t314) * t273) * m(6) + t436) * qJD(4) + t16, qJD(4) * t17 + t363 * t487 + t16;];
Cq = t1;
