% Calculate vector of inverse dynamics joint torques for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:45:00
% DurationCPUTime: 12.05s
% Computational Cost: add. (6186->450), mult. (8346->566), div. (0->0), fcn. (6481->6), ass. (0->239)
t499 = Icges(5,4) + Icges(4,5);
t239 = pkin(6) + qJ(3);
t220 = cos(t239);
t219 = sin(t239);
t385 = Icges(4,4) * t219;
t163 = Icges(4,1) * t220 - t385;
t206 = Icges(5,5) * t219;
t286 = Icges(5,1) * t220 + t206;
t484 = t163 + t286;
t243 = sin(qJ(1));
t244 = cos(qJ(1));
t498 = -t243 * t286 + t499 * t244;
t368 = t219 * t243;
t190 = Icges(4,4) * t368;
t366 = t220 * t243;
t489 = -Icges(4,1) * t366 + t190 + t498;
t488 = t499 * t243 + t484 * t244;
t284 = Icges(5,3) * t220 - t206;
t497 = Icges(4,2) * t220 + t284 + t385;
t382 = Icges(5,5) * t220;
t160 = Icges(5,1) * t219 - t382;
t207 = Icges(4,4) * t220;
t495 = Icges(4,1) * t219 + t160 + t207;
t155 = Icges(4,5) * t220 - Icges(4,6) * t219;
t106 = Icges(4,3) * t243 + t155 * t244;
t157 = Icges(5,4) * t220 + Icges(5,6) * t219;
t108 = Icges(5,2) * t243 + t157 * t244;
t496 = t106 + t108;
t376 = Icges(4,3) * t244;
t105 = Icges(4,5) * t366 - Icges(4,6) * t368 - t376;
t153 = Icges(5,3) * t219 + t382;
t103 = -Icges(5,6) * t244 + t153 * t243;
t379 = Icges(4,6) * t244;
t109 = Icges(4,4) * t366 - Icges(4,2) * t368 - t379;
t374 = t109 * t219;
t458 = -t103 * t219 + t489 * t220 + t374;
t494 = -t244 * t105 - t243 * t458;
t365 = t220 * t244;
t189 = Icges(5,5) * t365;
t367 = t219 * t244;
t378 = Icges(5,6) * t243;
t104 = Icges(5,3) * t367 + t189 + t378;
t493 = -t104 * t368 + t108 * t244 - t488 * t366;
t491 = t103 - t109;
t285 = -Icges(4,2) * t219 + t207;
t110 = Icges(4,6) * t243 + t244 * t285;
t490 = t104 - t110;
t487 = t153 - t285;
t486 = (Icges(4,6) - Icges(5,6)) * t220 + t499 * t219;
t483 = t497 * qJD(3);
t482 = t495 * qJD(3);
t481 = t104 * t367 + t243 * t496 + t488 * t365;
t107 = -Icges(5,2) * t244 + t157 * t243;
t97 = t243 * t107;
t480 = -t103 * t367 - t243 * t105 + t365 * t489 - t97;
t468 = -t219 * t497 + t220 * t495;
t479 = -t244 * t106 - t493;
t363 = t244 * t107;
t437 = -t363 + t494;
t435 = -t109 * t367 - t480;
t434 = -t110 * t367 + t481;
t478 = t483 * t244 + (t243 * t285 - t103 - t379) * qJD(1);
t477 = t483 * t243 + (t153 * t244 - t110 + t378) * qJD(1);
t476 = -t482 * t244 + (-t163 * t243 + t498) * qJD(1);
t475 = -t488 * qJD(1) + t482 * t243;
t474 = -t488 * t219 + t490 * t220;
t431 = t489 * t219 + t491 * t220;
t436 = -t110 * t368 + t479;
t473 = t487 * qJD(3);
t472 = t484 * qJD(3);
t471 = -t155 - t157;
t469 = t486 * qJD(3);
t467 = -t219 * t495 - t220 * t497;
t373 = t110 * t219;
t466 = t104 * t219 + t488 * t220 - t373;
t427 = t486 * t244;
t428 = t486 * t243;
t433 = t468 * t243 - t427;
t432 = t468 * t244 + t428;
t465 = t496 * qJD(1);
t464 = t243 * (-t244 * t497 + t488) - t244 * (-Icges(4,2) * t366 - t284 * t243 - t190 - t489);
t463 = t244 ^ 2;
t447 = rSges(5,3) + qJ(4);
t462 = t486 * qJD(1) + t467 * qJD(3) + t473 * t219 + t472 * t220;
t461 = qJD(3) * t474 + t219 * t478 + t220 * t476 + t465;
t417 = qJD(1) * t107;
t460 = -qJD(1) * t105 - qJD(3) * t431 - t219 * t477 + t220 * t475 - t417;
t459 = -t491 * t244 + (-Icges(5,1) * t367 + t160 * t244 + t189 + t490) * t243;
t457 = t243 * t434 - t244 * t435;
t456 = t436 * t243 - t244 * t437;
t455 = t495 - t487;
t454 = -t497 + t484;
t451 = t468 * qJD(1) + qJD(3) * t471;
t450 = qJD(1) * t458 - t243 * t469 + t465;
t449 = -t417 - t469 * t244 + (-t155 * t243 + t376 - t466) * qJD(1);
t406 = rSges(5,1) + pkin(3);
t448 = t432 * qJD(1);
t297 = t220 * rSges(5,1) + t219 * rSges(5,3);
t446 = t220 * pkin(3) + t219 * qJ(4) + t297;
t445 = t433 * qJD(1);
t322 = qJD(1) * qJD(3);
t173 = -qJDD(3) * t244 + t243 * t322;
t324 = qJD(4) * t220;
t354 = -t446 * qJD(3) + t324;
t257 = qJDD(4) * t219 + (t324 + t354) * qJD(3);
t237 = t244 * rSges(5,2);
t353 = t446 * t243 - t237;
t225 = t244 * qJ(2);
t180 = pkin(1) * t243 - t225;
t241 = cos(pkin(6));
t211 = pkin(2) * t241 + pkin(1);
t242 = -pkin(5) - qJ(2);
t214 = t244 * t242;
t339 = -t243 * t211 - t214;
t101 = t180 + t339;
t360 = t101 - t180;
t303 = -t353 + t360;
t325 = qJD(4) * t219;
t315 = t243 * t325;
t323 = qJD(1) * qJD(2);
t335 = qJDD(2) * t243 + t244 * t323;
t165 = rSges(5,1) * t219 - rSges(5,3) * t220;
t345 = pkin(3) * t219 - qJ(4) * t220 + t165;
t224 = t243 * qJ(2);
t182 = t244 * pkin(1) + t224;
t223 = qJD(2) * t244;
t149 = qJD(1) * t182 - t223;
t329 = qJD(1) * t243;
t202 = t242 * t329;
t402 = pkin(1) - t211;
t388 = -t149 + t202 - (-t244 * t402 - t224) * qJD(1);
t234 = t243 * rSges(5,2);
t327 = qJD(3) * t243;
t328 = qJD(1) * t244;
t400 = (pkin(3) * t328 + qJ(4) * t327) * t220 + (qJ(4) * t328 + (-pkin(3) * qJD(3) + qJD(4)) * t243) * t219 - t165 * t327 + (t244 * t297 + t234) * qJD(1);
t1 = t345 * t173 + t257 * t244 + t303 * qJDD(1) + (-t315 + t388 - t400) * qJD(1) + t335;
t444 = -g(1) + t1;
t443 = qJD(3) * t456 + t445;
t442 = qJD(3) * t457 + t448;
t441 = -t243 * t451 + t244 * t462;
t440 = t243 * t462 + t244 * t451;
t439 = t458 * qJD(3) + t219 * t475 + t220 * t477;
t438 = qJD(3) * t466 + t219 * t476 - t220 * t478;
t175 = qJD(1) * t180;
t429 = qJD(1) * t101 - t175;
t326 = qJD(3) * t244;
t316 = t220 * t326;
t426 = rSges(5,2) * t328 + t316 * t447;
t425 = t447 * t366;
t424 = t447 * t365;
t423 = -t219 * t464 + t459 * t220;
t422 = (-t219 * t455 + t220 * t454) * qJD(1);
t421 = t363 + t481;
t420 = t450 * t463 + (t461 * t243 + (-t449 + t460) * t244) * t243;
t419 = t460 * t463 + (t449 * t243 + (-t450 + t461) * t244) * t243;
t418 = t471 * qJD(1);
t240 = sin(pkin(6));
t395 = rSges(3,2) * t240;
t397 = rSges(3,1) * t241;
t121 = t243 * rSges(3,3) + (-t395 + t397) * t244;
t172 = qJDD(3) * t243 + t244 * t322;
t410 = t172 / 0.2e1;
t409 = t173 / 0.2e1;
t408 = t243 / 0.2e1;
t407 = -t244 / 0.2e1;
t405 = g(2) * t243;
t179 = t244 * t325;
t261 = -t219 * t326 - t220 * t329;
t317 = t219 * t329;
t401 = t406 * t261 - t447 * t317 + t179 + t426;
t396 = rSges(4,1) * t220;
t166 = rSges(4,1) * t219 + rSges(4,2) * t220;
t140 = t166 * t244;
t232 = t243 * rSges(4,3);
t119 = rSges(4,1) * t365 - rSges(4,2) * t367 + t232;
t305 = t244 * t211 - t242 * t243;
t102 = t305 - t182;
t359 = t102 + t182;
t44 = -t166 * t327 - t223 + (t119 + t359) * qJD(1);
t394 = t140 * t44;
t222 = qJD(2) * t243;
t298 = -t166 * t326 + t222;
t117 = rSges(4,1) * t366 - rSges(4,2) * t368 - t244 * rSges(4,3);
t319 = -t117 + t360;
t43 = qJD(1) * t319 + t298;
t392 = t243 * t43;
t352 = t406 * t365 + t447 * t367 + t234;
t30 = -t324 + (t243 * t353 + t244 * t352) * qJD(3);
t391 = t30 * t219;
t94 = t121 + t182;
t351 = -t406 * t368 + t425;
t350 = -t406 * t367 + t424;
t342 = rSges(4,2) * t317 + rSges(4,3) * t328;
t341 = t179 + t222;
t321 = t243 * t397;
t203 = t243 * t395;
t336 = t244 * rSges(3,3) + t203;
t120 = t321 - t336;
t340 = -t180 - t120;
t338 = rSges(3,3) * t328 + qJD(1) * t203;
t337 = t202 + t223;
t215 = qJ(2) * t328;
t334 = t215 + t222;
t318 = t406 * t244;
t314 = -pkin(1) - t397;
t311 = -t327 / 0.2e1;
t310 = t327 / 0.2e1;
t309 = -t326 / 0.2e1;
t308 = t326 / 0.2e1;
t306 = -t105 + t373;
t304 = t345 * qJD(3);
t183 = rSges(2,1) * t244 - rSges(2,2) * t243;
t181 = rSges(2,1) * t243 + rSges(2,2) * t244;
t169 = -rSges(4,2) * t219 + t396;
t262 = -t244 * t304 + t341;
t28 = qJD(1) * t303 + t262;
t29 = -t223 + (-t304 + t325) * t243 + (t352 + t359) * qJD(1);
t294 = t243 * t29 + t244 * t28;
t289 = -t243 * t44 - t244 * t43;
t73 = rSges(4,1) * t261 - rSges(4,2) * t316 + t342;
t136 = t166 * t243;
t75 = -qJD(3) * t136 + (t169 * t244 + t232) * qJD(1);
t288 = t243 * t75 + t244 * t73;
t277 = t117 * t243 + t119 * t244;
t271 = -qJDD(2) * t244 + qJD(1) * (-pkin(1) * t329 + t334) + qJDD(1) * t182 + t243 * t323;
t270 = t294 * t220;
t260 = t271 + qJD(1) * (-t215 + (t243 * t402 - t214) * qJD(1)) + qJDD(1) * t102;
t254 = -t211 - t446;
t151 = t169 * qJD(3);
t81 = qJD(1) * t94 - t223;
t80 = qJD(1) * t340 + t222;
t57 = t277 * qJD(3);
t42 = qJDD(1) * t121 + qJD(1) * (-qJD(1) * t321 + t338) + t271;
t41 = t340 * qJDD(1) + (-qJD(1) * t121 - t149) * qJD(1) + t335;
t17 = qJD(1) * t73 + qJDD(1) * t119 - t151 * t327 - t166 * t172 + t260;
t16 = -t151 * t326 + t166 * t173 + t319 * qJDD(1) + (-t75 + t388) * qJD(1) + t335;
t3 = -qJDD(4) * t220 - t352 * t173 + t353 * t172 + (t243 * t400 + t244 * t401 + t325) * qJD(3);
t2 = -t345 * t172 + t352 * qJDD(1) + (t179 + t401) * qJD(1) + t257 * t243 + t260;
t4 = [-m(2) * (-g(1) * t181 + g(2) * t183) + ((((t106 + t374) * t244 + t436 + t480 + t493) * t244 + (t421 + t437 - t494) * t243) * qJD(3) + t448) * t308 + (t468 * qJD(3) + t472 * t219 - t473 * t220) * qJD(1) + (t28 * (-t315 + t337) + t29 * (t341 + t426) + (-t29 * t219 * t318 + (t219 * t406 - t220 * t447) * t243 * t28) * qJD(3) + ((-t29 * t242 + t254 * t28) * t244 + (-t28 * rSges(5,2) + t254 * t29) * t243) * qJD(1) - (-qJD(1) * t353 + t262 - t28 + t429) * t29 + (t2 - g(2)) * (t305 + t352) + t444 * (t237 + (-t219 * t447 - t220 * t406) * t243 + t339)) * m(5) + (t43 * t337 + t44 * (t222 + t342) + (t166 * t392 - t394) * qJD(3) + ((-t43 * rSges(4,3) + t44 * (-t211 - t396)) * t243 + (t43 * (-t169 - t211) - t44 * t242) * t244) * qJD(1) - (-qJD(1) * t117 + t298 + t429 - t43) * t44 + (t17 - g(2)) * (t119 + t305) + (t16 - g(1)) * (-t117 + t339)) * m(4) + (-(-qJD(1) * t120 - t175 + t222 - t80) * t81 + t80 * t223 + t81 * (t334 + t338) + (t80 * (t314 + t395) * t244 + (t80 * (-rSges(3,3) - qJ(2)) + t81 * t314) * t243) * qJD(1) + (t42 - g(2)) * t94 + (t41 - g(1)) * (t243 * t314 + t225 + t336)) * m(3) + (-t474 + t432) * t410 + (-t431 + t433) * t409 + (t438 + t441) * t310 + (((t244 * t306 - t421 + t434) * t244 + (t243 * t306 + t435 - t479 - t97) * t243) * qJD(3) + t443 - t445) * t311 + (m(2) * (t181 ^ 2 + t183 ^ 2) + Icges(3,2) * t241 ^ 2 + (Icges(3,1) * t240 + 0.2e1 * Icges(3,4) * t241) * t240 + Icges(2,3) - t467) * qJDD(1) + (-t439 + t440 + t442) * t309; (-m(3) - m(4) - m(5)) * (g(1) * t243 - g(2) * t244) + 0.2e1 * (t1 * t408 + t2 * t407) * m(5) + 0.2e1 * (t16 * t408 + t17 * t407) * m(4) + 0.2e1 * (t407 * t42 + t408 * t41) * m(3); t457 * t410 + t456 * t409 + (qJD(1) * t441 + t419 * qJD(3) + t432 * qJDD(1) + t434 * t172 + t435 * t173) * t408 + (qJD(1) * t440 + t420 * qJD(3) + qJDD(1) * t433 + t172 * t436 + t173 * t437) * t407 - ((t459 * t219 + t220 * t464) * qJD(3) + (t454 * t219 + t455 * t220) * qJD(1)) * qJD(1) / 0.2e1 + (t439 * t244 + t438 * t243 + (-t243 * t431 - t244 * t474) * qJD(1)) * qJD(1) / 0.2e1 + (-t243 * t474 + t244 * t431) * qJDD(1) / 0.2e1 + t443 * t329 / 0.2e1 + t442 * t328 / 0.2e1 + ((-t327 * t427 - t418) * t243 + ((t243 * t428 + t423) * qJD(3) + t422) * t244) * t311 + ((t243 * t435 + t244 * t434) * qJD(1) + t419) * t310 + ((t243 * t437 + t244 * t436) * qJD(1) + t420) * t309 + ((-t326 * t428 + t418) * t244 + ((t244 * t427 + t423) * qJD(3) + t422) * t243) * t308 + (-g(1) * t424 - g(2) * t425 - g(3) * t446 - (-g(1) * t318 - t405 * t406) * t219 + (-t1 * t345 + t28 * t354 + t3 * t352 + t30 * t401 + (-t29 * t345 + t30 * t353) * qJD(1)) * t244 + (-t2 * t345 + t29 * t354 + t3 * t353 + t30 * t400 + (t28 * t345 - t30 * t352) * qJD(1)) * t243 - (t270 + t391) * qJD(4) - (-t28 * t351 + t29 * t350) * qJD(1) - ((-t28 * t446 + t30 * t350) * t244 + (-t29 * t446 + t30 * t351) * t243) * qJD(3)) * m(5) + (-(t136 * t43 - t394) * qJD(1) - (t57 * (-t136 * t243 - t140 * t244) + t289 * t169) * qJD(3) + (qJD(3) * t288 + t117 * t172 - t119 * t173) * t277 + t57 * ((t117 * t244 - t119 * t243) * qJD(1) + t288) + t289 * t151 + (-t16 * t244 - t17 * t243 + (-t244 * t44 + t392) * qJD(1)) * t166 + g(1) * t140 + g(2) * t136 - g(3) * t169) * m(4); (-(t270 + (t243 ^ 2 + t463) * t391) * qJD(3) + (qJD(3) * t294 + g(3) - t3) * t220 + (qJD(3) * t30 + t2 * t243 + t444 * t244 - t405) * t219) * m(5);];
tau = t4;
