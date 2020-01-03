% Calculate vector of inverse dynamics joint torques for
% S5RPPRP2
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:23
% DurationCPUTime: 15.18s
% Computational Cost: add. (10207->499), mult. (8829->597), div. (0->0), fcn. (6733->8), ass. (0->270)
t251 = pkin(8) + qJ(4);
t246 = sin(t251);
t248 = cos(t251);
t164 = Icges(5,5) * t248 - Icges(5,6) * t246;
t252 = qJ(1) + pkin(7);
t247 = sin(t252);
t249 = cos(t252);
t107 = Icges(5,3) * t247 + t164 * t249;
t166 = Icges(6,4) * t248 + Icges(6,6) * t246;
t109 = Icges(6,2) * t247 + t166 * t249;
t522 = t107 + t109;
t229 = Icges(6,5) * t246;
t303 = Icges(6,1) * t248 + t229;
t113 = Icges(6,4) * t247 + t249 * t303;
t407 = Icges(5,4) * t246;
t172 = Icges(5,1) * t248 - t407;
t115 = Icges(5,5) * t247 + t172 * t249;
t516 = t113 + t115;
t301 = Icges(6,3) * t248 - t229;
t520 = Icges(5,2) * t248 + t301 + t407;
t404 = Icges(6,5) * t248;
t169 = Icges(6,1) * t246 - t404;
t232 = Icges(5,4) * t248;
t519 = Icges(5,1) * t246 + t169 + t232;
t162 = Icges(6,3) * t246 + t404;
t302 = -Icges(5,2) * t246 + t232;
t521 = t162 - t302;
t513 = (Icges(5,6) - Icges(6,6)) * t248 + (Icges(6,4) + Icges(5,5)) * t246;
t518 = t172 + t303;
t385 = t248 * t249;
t195 = Icges(6,5) * t385;
t387 = t246 * t249;
t400 = Icges(6,6) * t247;
t105 = Icges(6,3) * t387 + t195 + t400;
t517 = t105 * t387 + t522 * t247 + t516 * t385;
t509 = -t246 * t520 + t248 * t519;
t386 = t247 * t248;
t388 = t246 * t247;
t398 = Icges(5,3) * t249;
t106 = Icges(5,5) * t386 - Icges(5,6) * t388 - t398;
t196 = Icges(5,4) * t388;
t405 = Icges(5,5) * t249;
t114 = Icges(5,1) * t386 - t196 - t405;
t401 = Icges(5,6) * t249;
t110 = Icges(5,4) * t386 - Icges(5,2) * t388 - t401;
t394 = t110 * t246;
t296 = -t114 * t248 + t394;
t108 = -Icges(6,2) * t249 + t166 * t247;
t395 = t108 * t249;
t104 = -Icges(6,6) * t249 + t162 * t247;
t112 = -Icges(6,4) * t249 + t247 * t303;
t299 = t104 * t246 + t112 * t248;
t454 = t247 * t299;
t36 = -t395 + t454;
t460 = -t106 * t249 - t247 * t296 + t36;
t111 = Icges(5,6) * t247 + t249 * t302;
t457 = -t111 * t387 + t517;
t321 = -t105 * t388 + t109 * t249 - t113 * t386;
t88 = t115 * t386;
t329 = t107 * t249 - t88;
t39 = -t111 * t388 - t329;
t459 = -t321 + t39;
t515 = t521 * qJD(4);
t514 = t518 * qJD(4);
t512 = -t164 - t166;
t511 = t520 * qJD(4);
t510 = t519 * qJD(4);
t508 = -t519 * t246 - t520 * t248;
t451 = t513 * t249;
t452 = t513 * t247;
t478 = t247 * t509 - t451;
t477 = t249 * t509 + t452;
t507 = t104 - t110;
t506 = t105 - t111;
t505 = -t112 - t114;
t255 = -pkin(6) - qJ(3);
t214 = t249 * t255;
t254 = cos(pkin(8));
t244 = pkin(3) * t254 + pkin(2);
t363 = -t247 * t244 - t214;
t256 = sin(qJ(1));
t428 = pkin(1) * t256;
t503 = t363 - t428;
t502 = t511 * t249 + (t247 * t302 - t104 - t401) * qJD(1);
t501 = t511 * t247 + (t162 * t249 - t111 + t400) * qJD(1);
t500 = -t510 * t249 + (-t172 * t247 - t112 + t405) * qJD(1);
t499 = -t516 * qJD(1) + t510 * t247;
t498 = t513 * qJD(1) + t508 * qJD(4) + t515 * t246 + t514 * t248;
t344 = rSges(5,1) * t386;
t497 = -t344 + t503;
t393 = t111 * t246;
t496 = t105 * t246 + t516 * t248 - t393;
t495 = t296 - t299;
t421 = -t247 * t106 - t114 * t385;
t42 = -t110 * t387 - t421;
t98 = t247 * t108;
t40 = t104 * t387 + t112 * t385 + t98;
t463 = t249 * t40;
t494 = t457 * t247 - t249 * t42 - t463;
t493 = t459 * t247 - t460 * t249;
t492 = t509 * qJD(1) + t512 * qJD(4);
t257 = cos(qJ(1));
t250 = t257 * pkin(1);
t491 = t477 * qJD(1);
t490 = -t246 * t516 + t506 * t248;
t456 = t505 * t246 + t507 * t248;
t177 = rSges(3,1) * t247 + rSges(3,2) * t249;
t147 = -t177 - t428;
t488 = t513 * qJD(4);
t487 = t478 * qJD(1);
t346 = qJD(1) * qJD(4);
t158 = -qJDD(4) * t249 + t247 * t346;
t349 = qJD(5) * t248;
t314 = t248 * rSges(6,1) + t246 * rSges(6,3);
t461 = t248 * pkin(4) + t246 * qJ(5) + t314;
t375 = -qJD(4) * t461 + t349;
t272 = qJDD(5) * t246 + (t349 + t375) * qJD(4);
t224 = t249 * qJ(3);
t176 = pkin(2) * t247 - t224;
t102 = t176 + t363;
t240 = t249 * rSges(6,2);
t377 = t247 * t461 - t240;
t323 = -t377 - t428;
t285 = t102 - t176 + t323;
t347 = qJD(1) * qJD(3);
t258 = qJD(1) ^ 2;
t485 = t258 * t250;
t288 = qJDD(3) * t247 + t249 * t347 - t485;
t350 = qJD(5) * t246;
t339 = t247 * t350;
t174 = rSges(6,1) * t246 - rSges(6,3) * t248;
t368 = pkin(4) * t246 - qJ(5) * t248 + t174;
t223 = t247 * qJ(3);
t181 = t249 * pkin(2) + t223;
t221 = qJD(3) * t249;
t136 = qJD(1) * t181 - t221;
t354 = qJD(1) * t247;
t207 = t255 * t354;
t424 = pkin(2) - t244;
t410 = -t136 + t207 - (-t249 * t424 - t223) * qJD(1);
t236 = t247 * rSges(6,2);
t352 = qJD(4) * t247;
t353 = qJD(1) * t249;
t422 = (pkin(4) * t353 + qJ(5) * t352) * t248 + (qJ(5) * t353 + (-pkin(4) * qJD(4) + qJD(5)) * t247) * t246 - t174 * t352 + (t249 * t314 + t236) * qJD(1);
t1 = t368 * t158 + t272 * t249 + t285 * qJDD(1) + (-t339 + t410 - t422) * qJD(1) + t288;
t486 = -g(1) + t1;
t484 = qJD(4) * t493 + t487;
t483 = qJD(4) * t494 + t491;
t482 = qJD(4) * t495 + t246 * t499 + t248 * t501;
t481 = t496 * qJD(4) + t500 * t246 - t248 * t502;
t480 = -t247 * t492 + t249 * t498;
t479 = t247 * t498 + t249 * t492;
t458 = t40 + t42;
t476 = t395 + t517;
t475 = t522 * qJD(1);
t474 = t247 * (-t249 * t520 + t516) - t249 * (-Icges(5,2) * t386 - t301 * t247 - t196 - t505);
t473 = t249 ^ 2;
t462 = rSges(6,3) + qJ(5);
t472 = t490 * qJD(4) + t246 * t502 + t500 * t248 + t475;
t441 = qJD(1) * t108;
t471 = -qJD(1) * t106 - t456 * qJD(4) - t246 * t501 + t248 * t499 - t441;
t470 = -t507 * t249 + (-Icges(6,1) * t387 + t169 * t249 + t195 + t506) * t247;
t469 = t519 - t521;
t468 = -t520 + t518;
t465 = qJD(1) * t495 - t488 * t247 + t475;
t464 = -t441 - t488 * t249 + (-t164 * t247 + t398 - t496) * qJD(1);
t429 = rSges(6,1) + pkin(4);
t160 = qJD(1) * t176;
t453 = qJD(1) * t102 - t160;
t253 = sin(pkin(8));
t417 = rSges(4,2) * t253;
t419 = rSges(4,1) * t254;
t121 = t247 * rSges(4,3) + (-t417 + t419) * t249;
t330 = t181 + t250;
t84 = t121 + t330;
t182 = t249 * rSges(3,1) - rSges(3,2) * t247;
t148 = t182 + t250;
t351 = qJD(4) * t249;
t340 = t248 * t351;
t450 = rSges(6,2) * t353 + t340 * t462;
t449 = t462 * t386;
t448 = t462 * t385;
t447 = -rSges(5,2) * t388 - t249 * rSges(5,3);
t446 = -t246 * t474 + t470 * t248;
t445 = (-t246 * t469 + t248 * t468) * qJD(1);
t444 = t465 * t473 + (t472 * t247 + (-t464 + t471) * t249) * t247;
t443 = t471 * t473 + (t464 * t247 + (-t465 + t472) * t249) * t247;
t442 = t512 * qJD(1);
t157 = qJDD(4) * t247 + t249 * t346;
t433 = t157 / 0.2e1;
t432 = t158 / 0.2e1;
t431 = t247 / 0.2e1;
t430 = -t249 / 0.2e1;
t427 = g(2) * t247;
t189 = t249 * t350;
t275 = -t246 * t351 - t248 * t354;
t341 = t246 * t354;
t423 = t275 * t429 - t341 * t462 + t189 + t450;
t175 = rSges(5,1) * t246 + rSges(5,2) * t248;
t143 = t175 * t249;
t234 = t247 * rSges(5,3);
t119 = rSges(5,1) * t385 - rSges(5,2) * t387 + t234;
t327 = t249 * t244 - t247 * t255;
t103 = t327 - t181;
t324 = t103 + t330;
t45 = -t175 * t352 - t221 + (t119 + t324) * qJD(1);
t416 = t143 * t45;
t117 = t344 + t447;
t331 = -t176 - t428;
t287 = t102 - t117 + t331;
t220 = qJD(3) * t247;
t316 = -t175 * t351 + t220;
t44 = qJD(1) * t287 + t316;
t414 = t247 * t44;
t376 = t385 * t429 + t387 * t462 + t236;
t31 = -t349 + qJD(2) + (t247 * t377 + t249 * t376) * qJD(4);
t413 = t31 * t246;
t374 = -t388 * t429 + t449;
t373 = -t387 * t429 + t448;
t365 = rSges(5,2) * t341 + rSges(5,3) * t353;
t364 = t189 + t220;
t208 = t247 * t417;
t362 = rSges(4,3) * t353 + qJD(1) * t208;
t361 = t207 + t221;
t360 = t249 * rSges(4,3) + t208;
t215 = qJ(3) * t353;
t359 = t215 + t220;
t348 = -m(4) - m(5) - m(6);
t345 = t247 * t419;
t342 = t429 * t249;
t338 = -pkin(2) - t419;
t335 = -t352 / 0.2e1;
t334 = t352 / 0.2e1;
t333 = -t351 / 0.2e1;
t332 = t351 / 0.2e1;
t328 = -t106 + t393;
t326 = t368 * qJD(4);
t325 = qJDD(1) * t250 - t258 * t428;
t120 = t345 - t360;
t322 = -t120 + t331;
t317 = t250 + t327;
t211 = rSges(2,1) * t257 - rSges(2,2) * t256;
t210 = rSges(2,1) * t256 + rSges(2,2) * t257;
t180 = rSges(5,1) * t248 - rSges(5,2) * t246;
t276 = -t249 * t326 + t364;
t29 = qJD(1) * t285 + t276;
t30 = -t221 + (-t326 + t350) * t247 + (t324 + t376) * qJD(1);
t311 = t247 * t30 + t249 * t29;
t306 = -t247 * t45 - t249 * t44;
t76 = rSges(5,1) * t275 - rSges(5,2) * t340 + t365;
t139 = t175 * t247;
t78 = -qJD(4) * t139 + (t180 * t249 + t234) * qJD(1);
t305 = t247 * t78 + t249 * t76;
t294 = t117 * t247 + t119 * t249;
t284 = t311 * t248;
t269 = -qJDD(3) * t249 + qJD(1) * (-pkin(2) * t354 + t359) + qJDD(1) * t181 + t247 * t347 + t325;
t268 = -t244 - t461;
t267 = qJDD(1) * t103 + t269 + qJD(1) * (-t215 + (t247 * t424 - t214) * qJD(1));
t156 = t180 * qJD(4);
t80 = qJD(1) * t84 - t221;
t79 = qJD(1) * t322 + t220;
t48 = qJD(4) * t294 + qJD(2);
t33 = qJDD(1) * t121 + qJD(1) * (-qJD(1) * t345 + t362) + t269;
t32 = t322 * qJDD(1) + (-qJD(1) * t121 - t136) * qJD(1) + t288;
t18 = qJD(4) * t305 + t117 * t157 - t119 * t158 + qJDD(2);
t17 = qJD(1) * t76 + qJDD(1) * t119 - t156 * t352 - t157 * t175 + t267;
t16 = -t156 * t351 + t158 * t175 + (-t78 + t410) * qJD(1) + t287 * qJDD(1) + t288;
t3 = -qJDD(5) * t248 + qJDD(2) - t376 * t158 + t377 * t157 + (t247 * t422 + t249 * t423 + t350) * qJD(4);
t2 = -t368 * t157 + t376 * qJDD(1) + (t189 + t423) * qJD(1) + t272 * t247 + t267;
t4 = [-m(2) * (-g(1) * t210 + g(2) * t211) + (((t39 - t88 + (t107 + t394) * t249 + t421) * t249 - t463 + (t36 - t454 + t476) * t247) * qJD(4) + t491) * t332 + (t509 * qJD(4) + t514 * t246 - t515 * t248) * qJD(1) + ((-t177 * t258 - g(2) + t325) * t148 + (-t485 + (-0.2e1 * t182 - t250 + t148) * t258 - g(1)) * t147) * m(3) + (t29 * (-t339 + t361) + t30 * (t364 + t450) + (-t30 * t246 * t342 + (t246 * t429 - t248 * t462) * t247 * t29) * qJD(4) + ((-t256 * t30 - t257 * t29) * pkin(1) + (-t30 * t255 + t268 * t29) * t249 + (-t29 * rSges(6,2) + t268 * t30) * t247) * qJD(1) - (qJD(1) * t323 + t276 - t29 + t453) * t30 + (-g(2) + t2) * (t317 + t376) + t486 * (t240 + (-t246 * t462 - t248 * t429) * t247 + t503)) * m(6) + ((t175 * t414 - t416) * qJD(4) + (t17 - g(2)) * (t119 + t317) + (t16 - g(1)) * (-t447 + t497) + (t361 + (-t234 - t250 + (-t180 - t244) * t249) * qJD(1)) * t44 + (t220 + t365 + t44 - t316 - t453 + (t117 + t428 + t497) * qJD(1)) * t45) * m(5) + (-(-t160 + t220 - t79 + (-t120 - t428) * qJD(1)) * t80 + t79 * t221 + t80 * (t359 + t362) + ((-t256 * t80 - t257 * t79) * pkin(1) + t79 * (t338 + t417) * t249 + (t79 * (-rSges(4,3) - qJ(3)) + t80 * t338) * t247) * qJD(1) + (t33 - g(2)) * t84 + (t32 - g(1)) * (t338 * t247 + t224 + t360 - t428)) * m(4) + (-t490 + t477) * t433 + (-t456 + t478) * t432 + (((t249 * t328 + t457 - t476) * t249 + (t247 * t328 + t321 + t329 + t458 - t98) * t247) * qJD(4) + t484 - t487) * t335 + (t480 + t481) * t334 + (m(3) * (t147 ^ 2 + t182 * t148) + m(2) * (t210 ^ 2 + t211 ^ 2) + Icges(4,2) * t254 ^ 2 + (Icges(4,1) * t253 + 0.2e1 * Icges(4,4) * t254) * t253 + Icges(2,3) + Icges(3,3) - t508) * qJDD(1) + (t479 - t482 + t483) * t333; (m(3) + m(4)) * qJDD(2) + m(5) * t18 + m(6) * t3 + (-m(3) + t348) * g(3); t348 * (g(1) * t247 - g(2) * t249) + 0.2e1 * (t1 * t431 + t2 * t430) * m(6) + 0.2e1 * (t16 * t431 + t17 * t430) * m(5) + 0.2e1 * (t32 * t431 + t33 * t430) * m(4); t494 * t433 + t493 * t432 + (qJD(1) * t480 + t443 * qJD(4) + qJDD(1) * t477 + t457 * t157 + t458 * t158) * t431 + (qJD(1) * t479 + t444 * qJD(4) + qJDD(1) * t478 + t459 * t157 + t460 * t158) * t430 - ((t470 * t246 + t248 * t474) * qJD(4) + (t468 * t246 + t469 * t248) * qJD(1)) * qJD(1) / 0.2e1 + (t482 * t249 + t481 * t247 + (-t456 * t247 - t249 * t490) * qJD(1)) * qJD(1) / 0.2e1 + (-t247 * t490 + t456 * t249) * qJDD(1) / 0.2e1 + t484 * t354 / 0.2e1 + t483 * t353 / 0.2e1 + ((-t352 * t451 - t442) * t247 + ((t247 * t452 + t446) * qJD(4) + t445) * t249) * t335 + ((t247 * t458 + t249 * t457) * qJD(1) + t443) * t334 + ((t247 * t460 + t459 * t249) * qJD(1) + t444) * t333 + ((-t351 * t452 + t442) * t249 + ((t249 * t451 + t446) * qJD(4) + t445) * t247) * t332 + (-g(1) * t448 - g(2) * t449 - g(3) * t461 - (-g(1) * t342 - t427 * t429) * t246 + (-t1 * t368 + t29 * t375 + t3 * t376 + t31 * t423 + (-t30 * t368 + t31 * t377) * qJD(1)) * t249 + (-t2 * t368 + t30 * t375 + t3 * t377 + t31 * t422 + (t29 * t368 - t31 * t376) * qJD(1)) * t247 - (t284 + t413) * qJD(5) - (-t29 * t374 + t30 * t373) * qJD(1) - ((-t29 * t461 + t31 * t373) * t249 + (-t30 * t461 + t31 * t374) * t247) * qJD(4)) * m(6) + (-(t139 * t44 - t416) * qJD(1) - (t48 * (-t139 * t247 - t143 * t249) + t306 * t180) * qJD(4) + t18 * t294 + t48 * ((t117 * t249 - t119 * t247) * qJD(1) + t305) + t306 * t156 + (-t16 * t249 - t17 * t247 + (-t249 * t45 + t414) * qJD(1)) * t175 + g(1) * t143 + g(2) * t139 - g(3) * t180) * m(5); (-(t284 + (t247 ^ 2 + t473) * t413) * qJD(4) + (qJD(4) * t311 + g(3) - t3) * t248 + (qJD(4) * t31 + t2 * t247 + t249 * t486 - t427) * t246) * m(6);];
tau = t4;
