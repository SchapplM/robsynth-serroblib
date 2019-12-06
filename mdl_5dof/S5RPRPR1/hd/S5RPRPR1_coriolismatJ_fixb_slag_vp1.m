% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:04
% EndTime: 2019-12-05 17:47:20
% DurationCPUTime: 6.52s
% Computational Cost: add. (14561->418), mult. (15056->557), div. (0->0), fcn. (13640->8), ass. (0->248)
t308 = sin(qJ(3));
t310 = cos(qJ(3));
t275 = rSges(4,1) * t310 - rSges(4,2) * t308;
t309 = sin(qJ(1));
t249 = t275 * t309;
t311 = cos(qJ(1));
t250 = t275 * t311;
t142 = -t249 * t311 + t309 * t250;
t372 = t309 * t310;
t283 = pkin(3) * t372;
t304 = qJ(3) + pkin(8);
t288 = cos(t304);
t375 = t288 * t309;
t287 = sin(t304);
t376 = t287 * t309;
t346 = rSges(5,1) * t375 - rSges(5,2) * t376;
t179 = t283 + t346;
t443 = m(6) / 0.2e1;
t444 = m(5) / 0.2e1;
t445 = m(4) / 0.2e1;
t259 = rSges(5,1) * t288 - rSges(5,2) * t287;
t404 = pkin(3) * t310;
t344 = (t259 + t404) * t311;
t465 = t309 * t344;
t289 = qJ(5) + t304;
t285 = cos(t289);
t377 = t285 * t309;
t284 = sin(t289);
t380 = t284 * t309;
t215 = rSges(6,1) * t377 - rSges(6,2) * t380;
t403 = pkin(4) * t288;
t265 = t403 + t404;
t147 = t309 * t265 + t215;
t242 = rSges(6,1) * t285 - rSges(6,2) * t284;
t454 = (t242 + t265) * t311;
t472 = t309 * t454;
t85 = -t147 * t311 + t472;
t351 = (-t179 * t311 + t465) * t444 + t142 * t445 + t85 * t443;
t144 = t283 + (t242 + t403) * t309;
t196 = t259 * t309 + t283;
t394 = (t196 * t311 - t465) * t444 + (t144 * t311 - t472) * t443;
t21 = t394 - t351;
t479 = t21 * qJD(1);
t141 = t311 * t454;
t464 = t311 * t344;
t392 = (-t196 * t309 - t464) * t444 + (-t144 * t309 - t141) * t443;
t395 = (t309 * t179 + t464) * t444 + (t309 * t147 + t141) * t443;
t26 = t395 - t392;
t478 = t26 * qJD(1);
t386 = Icges(6,4) * t284;
t335 = Icges(6,2) * t285 + t386;
t174 = Icges(6,6) * t311 + t335 * t309;
t266 = Icges(6,4) * t377;
t176 = Icges(6,1) * t380 + Icges(6,5) * t311 + t266;
t327 = -t174 * t285 - t176 * t284;
t477 = t311 * t327;
t385 = Icges(6,4) * t285;
t238 = -Icges(6,2) * t284 + t385;
t338 = Icges(6,1) * t284 + t385;
t476 = t238 + t338;
t475 = Icges(4,5) * t310 + Icges(5,5) * t288 - Icges(4,6) * t308 - Icges(5,6) * t287;
t305 = t309 ^ 2;
t473 = -t309 / 0.2e1;
t417 = t309 / 0.2e1;
t415 = t311 / 0.2e1;
t306 = t311 ^ 2;
t279 = t305 + t306;
t341 = rSges(6,1) * t284 + rSges(6,2) * t285;
t139 = t279 * t341;
t329 = Icges(6,5) * t284 + Icges(6,6) * t285;
t471 = t311 * t329;
t470 = t329 * t309;
t383 = t341 * t309;
t342 = rSges(5,1) * t287 + rSges(5,2) * t288;
t405 = pkin(3) * t308;
t457 = t342 + t405;
t197 = t457 * t311;
t264 = pkin(4) * t287 + t405;
t319 = t341 + t264;
t177 = -Icges(6,5) * t309 + t311 * t338;
t367 = t238 * t311 + t177;
t368 = -Icges(6,2) * t380 + t176 + t266;
t175 = -Icges(6,6) * t309 + t311 * t335;
t240 = Icges(6,1) * t285 - t386;
t369 = -t240 * t311 + t175;
t370 = -t240 * t309 + t174;
t469 = -(t369 * t309 - t370 * t311) * t284 + (t367 * t309 - t368 * t311) * t285;
t387 = Icges(5,4) * t288;
t255 = -Icges(5,2) * t287 + t387;
t389 = Icges(4,4) * t310;
t270 = -Icges(4,2) * t308 + t389;
t339 = Icges(5,1) * t287 + t387;
t340 = Icges(4,1) * t308 + t389;
t468 = (t339 / 0.2e1 + t255 / 0.2e1) * t288 + (t340 / 0.2e1 + t270 / 0.2e1) * t310;
t330 = Icges(6,5) * t285 - Icges(6,6) * t284;
t209 = t309 * t330;
t210 = t330 * t311;
t402 = (-t305 * t210 + (t309 * t209 + t469) * t311) * t417 + (t306 * t209 + (-t311 * t210 - t469) * t309) * t415;
t216 = t242 * t311;
t194 = t311 * t216;
t460 = t309 * t215 + t194;
t313 = -t309 * rSges(6,3) + t311 * t341;
t166 = t311 * t313;
t178 = rSges(6,3) * t311 + t383;
t307 = -qJ(4) - pkin(6);
t353 = t311 * t405;
t218 = t311 * (t353 + (pkin(6) + t307) * t309);
t286 = t311 * t307;
t374 = t308 * t309;
t235 = pkin(3) * t374 - pkin(6) * t311 - t286;
t354 = -pkin(7) + t307;
t282 = t311 * t354;
t356 = t311 * t264 + t309 * t354;
t52 = t311 * (t309 * t307 + t353 - t356) - t166 - t218 + (-t235 + t282 - t286 - t178 + (-t264 + t405) * t309) * t309;
t6 = t402 + m(6) * (-t141 * t341 - t144 * t383 - t460 * t52);
t466 = t6 * qJD(5);
t390 = Icges(4,4) * t308;
t337 = Icges(4,2) * t310 + t390;
t222 = -Icges(4,6) * t309 + t311 * t337;
t224 = -Icges(4,5) * t309 + t311 * t340;
t463 = (t222 * t310 + t224 * t308) * t311;
t388 = Icges(5,4) * t287;
t336 = Icges(5,2) * t288 + t388;
t202 = -Icges(5,6) * t309 + t311 * t336;
t204 = -Icges(5,5) * t309 + t311 * t339;
t462 = (t202 * t288 + t204 * t287) * t311;
t461 = (t175 * t285 + t177 * t284) * t311;
t459 = t475 * t311;
t458 = t475 * t309;
t453 = t476 * t284 - (-t335 + t240) * t285;
t348 = -t476 * t285 / 0.2e1 + (t335 / 0.2e1 - t240 / 0.2e1) * t284;
t173 = -Icges(6,3) * t309 + t471;
t65 = t311 * (Icges(6,3) * t311 + t470) + t174 * t377 + t176 * t380;
t66 = -t311 * t173 - t175 * t377 - t177 * t380;
t68 = -t309 * t173 + t461;
t5 = ((t66 - t477) * t309 + (t68 - t461 + (t173 + t327) * t309 + t65) * t311) * t417 + (t66 * t309 + t311 * t65) * t473 + (t305 * t173 + t68 * t309 + (t66 + (t173 - t327) * t311 + t477) * t311) * t415;
t359 = t270 * t311 + t224;
t272 = Icges(4,1) * t310 - t390;
t361 = -t272 * t311 + t222;
t363 = t255 * t311 + t204;
t257 = Icges(5,1) * t288 - t388;
t365 = -t257 * t311 + t202;
t450 = t365 * t287 - t363 * t288 + t361 * t308 - t359 * t310;
t280 = Icges(4,4) * t372;
t223 = Icges(4,1) * t374 + Icges(4,5) * t311 + t280;
t360 = -Icges(4,2) * t374 + t223 + t280;
t221 = Icges(4,6) * t311 + t337 * t309;
t362 = -t272 * t309 + t221;
t277 = Icges(5,4) * t375;
t203 = Icges(5,1) * t376 + Icges(5,5) * t311 + t277;
t364 = -Icges(5,2) * t376 + t203 + t277;
t201 = Icges(5,6) * t311 + t336 * t309;
t366 = -t257 * t309 + t201;
t449 = t366 * t287 - t364 * t288 + t362 * t308 - t360 * t310;
t448 = 0.2e1 * t279;
t447 = 0.4e1 * qJD(1);
t446 = 2 * qJD(3);
t442 = -pkin(1) - pkin(6);
t290 = t311 * qJ(2);
t343 = rSges(4,1) * t308 + rSges(4,2) * t310;
t312 = -t309 * rSges(4,3) + t311 * t343;
t161 = t442 * t309 + t290 + t312;
t162 = (rSges(4,3) - t442) * t311 + (qJ(2) + t343) * t309;
t441 = m(4) * (t161 * t250 + t162 * t249);
t301 = t309 * rSges(5,3);
t133 = t290 - t301 + (-pkin(1) + t307) * t309 + t197;
t134 = -t286 + (rSges(5,3) + pkin(1)) * t311 + (qJ(2) + t457) * t309;
t440 = m(5) * (t133 * t344 + t134 * t179);
t439 = m(5) * (-t133 * t309 + t311 * t134);
t438 = m(5) * (t133 * t311 + t309 * t134);
t121 = -t282 + (rSges(6,3) + pkin(1)) * t311 + (qJ(2) + t319) * t309;
t110 = t311 * t121;
t120 = -t309 * pkin(1) + t290 + t313 + t356;
t345 = t110 * t341 - t120 * t383;
t434 = m(6) * (t144 * t216 - t215 * t454 + t345);
t433 = m(6) * (t242 * t85 + t345);
t432 = m(6) * (t120 * t454 + t121 * t147);
t431 = m(6) * (t120 * t216 + t121 * t215);
t430 = m(6) * (-t120 * t309 + t110);
t429 = m(6) * (t120 * t311 + t309 * t121);
t124 = -t215 * t311 + t309 * t216;
t422 = t124 / 0.2e1;
t421 = -t279 / 0.2e1;
t416 = -t311 / 0.2e1;
t414 = m(3) * ((rSges(3,3) * t311 + t290) * t311 + (rSges(3,3) + qJ(2)) * t305);
t413 = m(4) * (t161 * t311 + t162 * t309);
t401 = m(6) * qJD(5);
t317 = m(6) * t460;
t320 = m(6) * t242 * t421;
t69 = -t317 / 0.2e1 + t320;
t371 = t69 * qJD(1);
t128 = (t444 + t443) * t448;
t355 = t128 * qJD(1);
t352 = -m(6) * t124 / 0.2e1;
t331 = Icges(5,5) * t287 + Icges(5,6) * t288;
t199 = Icges(5,3) * t311 + t331 * t309;
t75 = t311 * t199 + t201 * t375 + t203 * t376;
t200 = -Icges(5,3) * t309 + t311 * t331;
t76 = -t311 * t200 - t202 * t375 - t204 * t376;
t333 = Icges(4,5) * t308 + Icges(4,6) * t310;
t219 = Icges(4,3) * t311 + t333 * t309;
t94 = t311 * t219 + t221 * t372 + t223 * t374;
t220 = -Icges(4,3) * t309 + t311 * t333;
t95 = -t311 * t220 - t222 * t372 - t224 * t374;
t349 = t279 * t343;
t143 = t319 * t309;
t145 = t319 * t311;
t328 = -t143 * t309 - t145 * t311;
t324 = -t201 * t288 - t203 * t287;
t322 = -t221 * t310 - t223 * t308;
t316 = -t5 + (t367 * t284 + t369 * t285 + t453 * t311 - t470) * t417 + (-t368 * t284 - t370 * t285 - t453 * t309 - t471) * t415;
t314 = -t348 + (t415 + t416) * (t175 * t284 - t177 * t285);
t206 = t309 * t219;
t195 = t457 * t309;
t169 = t309 * t199;
t130 = t139 * t401;
t129 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t448 + (m(5) + m(6)) * t421;
t119 = t401 * t422;
t113 = -t309 * t178 - t166;
t97 = -t309 * t220 + t463;
t96 = t311 * t322 + t206;
t78 = -t309 * t200 + t462;
t77 = t311 * t324 + t169;
t70 = t317 / 0.2e1 + t320;
t63 = -t306 * t265 - t194 + (t283 - t147) * t309 + (-t279 + t306) * t404;
t55 = t97 * t309 + t311 * t96;
t54 = t95 * t309 + t311 * t94;
t48 = t78 * t309 + t311 * t77;
t47 = t76 * t309 + t311 * t75;
t40 = t430 + t439;
t38 = t348 + t431;
t36 = t433 / 0.2e1;
t34 = t434 / 0.2e1;
t29 = t413 + t414 + t429 + t438;
t28 = t392 + t395;
t25 = t305 * t220 + (-t206 + t95 + (t220 - t322) * t311) * t311;
t24 = (-t96 + t206 + t95) * t309 + (t97 - t463 + (t220 + t322) * t309 + t94) * t311;
t22 = t351 + t394;
t19 = t305 * t200 + (-t169 + t76 + (t200 - t324) * t311) * t311;
t18 = (-t77 + t169 + t76) * t309 + (t78 - t462 + (t200 + t324) * t309 + t75) * t311;
t9 = (-t272 / 0.2e1 + t337 / 0.2e1) * t308 + (-t257 / 0.2e1 + t336 / 0.2e1) * t287 + t441 + t440 + t432 + t348 - t468;
t8 = m(6) * (-t113 * t460 - t242 * t139) + t402;
t7 = t8 * qJD(5);
t4 = t34 - t433 / 0.2e1 + t5;
t3 = t36 - t434 / 0.2e1 + t5;
t2 = t34 + t36 + t316;
t1 = (t48 / 0.2e1 + t25 / 0.2e1 + t55 / 0.2e1 + t19 / 0.2e1) * t311 + (t18 / 0.2e1 - t54 / 0.2e1 + t24 / 0.2e1 - t47 / 0.2e1) * t309 + t5;
t10 = [qJD(2) * t29 + qJD(3) * t9 + qJD(4) * t40 + qJD(5) * t38, qJD(1) * t29 + qJD(3) * t22 + qJD(4) * t129 + t119, t9 * qJD(1) + t22 * qJD(2) + t28 * qJD(4) + t2 * qJD(5) + ((t142 * t275 - (t161 * t309 - t162 * t311) * t343) * t445 + (-t133 * t195 + t134 * t197 + (-t179 + t196) * t344) * t444 + (-t120 * t143 + t121 * t145 + (t144 - t147) * t454) * t443) * t446 + ((-t333 - t331) * (t306 / 0.2e1 + t305 / 0.2e1) + t316 + (t18 + t24) * t473 + (-t287 * t364 - t288 * t366 - t308 * t360 - t310 * t362) * t415 + (t287 * t363 + t288 * t365 + t308 * t359 + t310 * t361 + t47 + t54) * t417 + (t48 + t25 + t55 + t19) * t416) * qJD(3), qJD(1) * t40 + qJD(2) * t129 + qJD(3) * t28 + qJD(5) * t70, t38 * qJD(1) + t2 * qJD(3) + t70 * qJD(4) + t316 * qJD(5) + (qJD(2) * t422 + (t124 * t242 + t345) * qJD(5)) * m(6); -t21 * qJD(3) - t128 * qJD(4) + t119 + (-t438 / 0.4e1 - t429 / 0.4e1 - t413 / 0.4e1 - t414 / 0.4e1) * t447, 0, -t479 - t130 + ((-t195 * t309 - t197 * t311) * t444 + t328 * t443 - t349 * t445) * t446, -t355, -t130 + 0.2e1 * (t124 * qJD(1) / 0.4e1 - t139 * qJD(3) / 0.2e1) * m(6); (t314 + (-t336 + t257) * t287 / 0.2e1 + (-t337 + t272) * t308 / 0.2e1 + t468) * qJD(1) + t21 * qJD(2) + t1 * qJD(3) - t26 * qJD(4) + t4 * qJD(5) + (-t441 / 0.4e1 - t432 / 0.4e1 - t440 / 0.4e1) * t447, t479, t1 * qJD(1) + (m(4) * ((-t311 * t312 + (-t311 * rSges(4,3) - t343 * t309) * t309) * (-t309 * t249 - t250 * t311) - t275 * t349) + m(6) * (-t143 * t144 - t145 * t454 + t52 * t63) + m(5) * (-t196 * t195 - t344 * t197 + (-t311 * (t311 * t342 - t301) - t218 + (-rSges(5,3) * t311 - t309 * t342 - t235) * t309) * (-t259 * t306 - t279 * t404 - t309 * t346)) + t402 + ((t449 * t311 + (-t450 + t458) * t309) * t311 - t459 * t305) * t417 + ((t450 * t309 + (-t449 - t459) * t311) * t309 + t458 * t306) * t415) * qJD(3) + t466, -t478, t4 * qJD(1) + t6 * qJD(3) + t466; t128 * qJD(2) + t26 * qJD(3) - t69 * qJD(5) + (-t439 / 0.4e1 - t430 / 0.4e1) * t447, t355, t478 + ((-t143 * t311 + t309 * t145) * t443 + (-t195 * t311 + t309 * t197) * t444) * t446, 0, -t371; (t314 - t431) * qJD(1) + qJD(2) * t352 + t3 * qJD(3) + t69 * qJD(4) + t5 * qJD(5), qJD(1) * t352, t3 * qJD(1) + ((t113 * t63 + t242 * t328) * m(6) + t402) * qJD(3) + t7, t371, qJD(1) * t5 + qJD(3) * t8 + t7;];
Cq = t10;
