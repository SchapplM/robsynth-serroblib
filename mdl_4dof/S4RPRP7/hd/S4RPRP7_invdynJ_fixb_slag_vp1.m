% Calculate vector of inverse dynamics joint torques for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:16
% DurationCPUTime: 11.58s
% Computational Cost: add. (3178->449), mult. (7968->551), div. (0->0), fcn. (6147->4), ass. (0->240)
t503 = Icges(5,4) + Icges(4,5);
t502 = Icges(4,6) - Icges(5,6);
t234 = sin(qJ(1));
t236 = cos(qJ(1));
t235 = cos(qJ(3));
t233 = sin(qJ(3));
t374 = Icges(4,4) * t233;
t271 = Icges(4,2) * t235 + t374;
t501 = t502 * t234 - t236 * t271;
t373 = Icges(4,4) * t235;
t273 = Icges(4,1) * t233 + t373;
t500 = -t503 * t234 + t236 * t273;
t355 = t235 * t236;
t357 = t233 * t236;
t495 = -Icges(5,5) * t357 + Icges(5,3) * t355 - t501;
t195 = Icges(5,5) * t355;
t499 = Icges(5,1) * t357 - t195 + t500;
t370 = Icges(5,5) * t233;
t165 = Icges(5,1) * t235 + t370;
t167 = Icges(4,1) * t235 - t374;
t488 = t165 + t167;
t269 = Icges(4,5) * t233 + Icges(4,6) * t235;
t97 = Icges(4,3) * t236 + t234 * t269;
t270 = Icges(5,4) * t233 - Icges(5,6) * t235;
t99 = Icges(5,2) * t236 + t234 * t270;
t498 = t97 + t99;
t369 = Icges(5,2) * t234;
t100 = Icges(5,4) * t357 - Icges(5,6) * t355 - t369;
t98 = -Icges(4,3) * t234 + t236 * t269;
t497 = -t100 - t98;
t101 = Icges(4,6) * t236 + t234 * t271;
t358 = t233 * t234;
t194 = Icges(5,5) * t358;
t356 = t234 * t235;
t367 = Icges(5,6) * t236;
t95 = -Icges(5,3) * t356 + t194 + t367;
t496 = -t101 + t95;
t466 = t499 * t233 + t495 * t235;
t494 = t466 * t236;
t222 = Icges(5,5) * t235;
t272 = Icges(5,1) * t233 - t222;
t103 = Icges(5,4) * t236 + t234 * t272;
t196 = Icges(4,4) * t356;
t371 = Icges(4,5) * t236;
t105 = Icges(4,1) * t358 + t196 + t371;
t493 = t103 + t105;
t491 = t488 * t236;
t490 = -t502 * t233 + t503 * t235;
t163 = -Icges(4,2) * t233 + t373;
t430 = Icges(5,3) * t233 + t222;
t489 = -t163 + t430;
t267 = -Icges(5,3) * t235 + t370;
t487 = t267 - t271;
t486 = -t272 - t273;
t436 = -t498 * t234 - t95 * t355;
t260 = t163 * t235 + t167 * t233;
t478 = t165 * t233 - t235 * t430 + t260;
t266 = t101 * t235 + t105 * t233;
t428 = -t103 * t233 - t266;
t393 = t103 * t358 + t236 * t99;
t31 = -t356 * t95 + t393;
t33 = t101 * t356 + t105 * t358 + t236 * t97;
t444 = t31 + t33;
t443 = t497 * t236 - t495 * t356 - t358 * t499;
t442 = -t103 * t357 - t266 * t236 - t436;
t441 = t234 * t497 + t494;
t124 = t163 * t236;
t323 = qJD(3) * t236;
t485 = qJD(3) * t124 - t430 * t323 + (t234 * t267 - t101 + t367) * qJD(1);
t117 = t430 * t234;
t324 = qJD(3) * t234;
t484 = qJD(3) * t117 - t163 * t324 + (t236 * t267 + t501) * qJD(1);
t483 = -t491 * qJD(3) + (t234 * t273 + t103 + t371) * qJD(1);
t127 = t167 * t234;
t482 = qJD(3) * t127 + t165 * t324 + (t236 * t272 + t500) * qJD(1);
t439 = t233 * t495 - t235 * t499;
t440 = t233 * t496 + t235 * t493;
t457 = t490 * t234;
t481 = t487 * qJD(3);
t480 = t486 * qJD(3);
t479 = t233 * t489 + t235 * t488;
t383 = t235 * t95;
t477 = t383 + t428;
t476 = -t269 - t270;
t453 = t490 * t236;
t438 = t478 * t234 + t453;
t437 = t165 * t357 + t236 * t260 - t355 * t430 - t457;
t475 = (Icges(4,2) * t358 + t117 - t196 - t493) * t236 + (-Icges(5,3) * t357 + t124 - t195 + t499) * t234;
t474 = (Icges(5,1) * t356 + t127 + t194 + t496) * t236 + (-t491 + t495) * t234;
t473 = t486 + t489;
t472 = t487 + t488;
t471 = t498 * qJD(1);
t470 = t236 ^ 2;
t404 = rSges(5,1) + pkin(3);
t378 = rSges(5,3) + qJ(4);
t469 = t478 * qJD(1) + t476 * qJD(3);
t363 = qJD(1) * t98;
t468 = t363 + t457 * qJD(3) + (t236 * t270 - t369 - t477) * qJD(1);
t467 = -t466 * qJD(1) - t453 * qJD(3) + t471;
t465 = t441 * t234 + t442 * t236;
t464 = t443 * t234 + t444 * t236;
t463 = -t440 * qJD(3) - t482 * t233 + t484 * t235 + t471;
t462 = -qJD(1) * t490 + t479 * qJD(3) + t480 * t233 + t481 * t235;
t461 = qJD(1) * t100 + t439 * qJD(3) + t483 * t233 - t485 * t235 + t363;
t460 = (t233 * t473 + t235 * t472) * qJD(1);
t333 = -t234 * rSges(5,2) - rSges(5,3) * t355;
t345 = -qJ(4) * t355 + t357 * t404 + t333;
t459 = t438 * qJD(1);
t229 = t236 * rSges(4,3);
t108 = rSges(4,1) * t358 + rSges(4,2) * t356 + t229;
t177 = t236 * pkin(1) + t234 * qJ(2);
t231 = t236 * pkin(5);
t429 = t231 + t177;
t458 = t108 + t429;
t174 = pkin(3) * t235 + qJ(4) * t233;
t175 = rSges(5,1) * t235 + rSges(5,3) * t233;
t335 = t174 + t175;
t216 = t235 * qJ(4);
t228 = t235 * rSges(5,3);
t286 = rSges(5,1) * t233 - t228;
t337 = -pkin(3) * t233 + t216 - t286;
t325 = qJD(1) * t236;
t456 = t233 * t325 + t235 * t324;
t455 = t476 * qJD(1);
t454 = t437 * qJD(1);
t451 = -t474 * t233 + t235 * t475;
t450 = qJD(3) * t464 + t459;
t449 = qJD(3) * t465 - t454;
t448 = t477 * qJD(3) + t484 * t233 + t482 * t235;
t447 = t466 * qJD(3) + t485 * t233 + t483 * t235;
t446 = t469 * t234 - t462 * t236;
t445 = t462 * t234 + t236 * t469;
t230 = t236 * rSges(5,2);
t431 = t404 * t358 + t230;
t178 = -rSges(3,2) * t236 + t234 * rSges(3,3);
t427 = t468 * t470 + (t461 * t234 + (-t463 + t467) * t236) * t234;
t426 = t463 * t470 + (t467 * t234 + (-t461 + t468) * t236) * t234;
t309 = t233 * t324;
t425 = t378 * t309 + t404 * t456;
t214 = qJD(2) * t236;
t115 = qJD(1) * t177 - t214;
t387 = rSges(4,2) * t235;
t287 = rSges(4,1) * t233 + t387;
t149 = t287 * qJD(3);
t320 = qJD(1) * qJD(3);
t152 = qJDD(3) * t234 + t236 * t320;
t388 = rSges(4,2) * t233;
t391 = rSges(4,1) * t235;
t176 = -t388 + t391;
t237 = qJD(1) ^ 2;
t321 = qJD(1) * qJD(2);
t332 = qJDD(2) * t234 + t236 * t321;
t258 = -t231 * t237 + t332;
t217 = t236 * qJ(2);
t171 = pkin(1) * t234 - t217;
t110 = -t234 * rSges(4,3) + t236 * t287;
t401 = pkin(5) * t234;
t300 = t110 - t401;
t295 = -t171 + t300;
t136 = t176 * t236;
t73 = -qJD(3) * t136 + (t234 * t287 + t229) * qJD(1);
t16 = -t149 * t324 + t152 * t176 + (-t115 - t73) * qJD(1) + t295 * qJDD(1) + t258;
t153 = qJDD(3) * t236 - t234 * t320;
t326 = qJD(1) * t234;
t207 = qJ(2) * t325;
t213 = qJD(2) * t234;
t331 = t207 + t213;
t315 = qJD(1) * (-pkin(1) * t326 + t331) + qJDD(1) * t177 + t234 * t321;
t365 = pkin(5) * qJDD(1);
t254 = t236 * t365 - t237 * t401 + t315;
t312 = t456 * rSges(4,1) + t325 * t387;
t316 = qJD(3) * t388;
t75 = (-rSges(4,3) * qJD(1) - t316) * t234 + t312;
t17 = qJD(1) * t75 + qJDD(1) * t108 - t153 * t176 + (qJD(3) * t149 - qJDD(2)) * t236 + t254;
t423 = t16 * t234 - t17 * t236;
t322 = qJD(4) * t235;
t135 = t175 * t236;
t395 = (pkin(3) * t326 - qJ(4) * t323) * t233 + (-qJ(4) * t326 + (-pkin(3) * qJD(3) + qJD(4)) * t236) * t235 - qJD(3) * t135 + (t234 * t286 + t230) * qJD(1);
t212 = qJD(4) * t233;
t344 = t337 * qJD(3) + t212;
t415 = -qJD(3) * (t212 + t344) + qJDD(4) * t235;
t2 = t335 * t152 + (-t171 + t345) * qJDD(1) + (-t236 * t322 - t115 - t395) * qJD(1) + (-t365 - t415) * t234 + t258;
t400 = g(2) * t236;
t293 = g(1) * t234 - t400;
t307 = t234 * t322;
t346 = -rSges(5,3) * t356 - t216 * t234 + t431;
t394 = -(-qJ(4) * t325 - qJD(4) * t234) * t235 - t333 * qJD(1) - t425;
t3 = -t335 * t153 + t346 * qJDD(1) + (-t307 - t394) * qJD(1) + (-qJDD(2) + t415) * t236 + t254;
t422 = t2 * t234 - t236 * t3 - t293;
t412 = t233 * t378 + t235 * t404;
t411 = t234 ^ 2;
t410 = -pkin(1) - pkin(5);
t408 = t152 / 0.2e1;
t407 = t153 / 0.2e1;
t406 = t234 / 0.2e1;
t403 = rSges(3,2) - pkin(1);
t386 = rSges(3,3) * t236;
t139 = t176 * t324;
t43 = qJD(1) * t295 + t139 + t213;
t381 = t236 * t43;
t30 = t212 + (-t234 * t346 - t236 * t345) * qJD(3);
t380 = t30 * t235;
t343 = t404 * t356 + t378 * t358;
t342 = -t174 * t236 - t135;
t172 = rSges(3,2) * t234 + t386;
t336 = -t171 + t172;
t113 = t177 + t178;
t330 = rSges(3,2) * t326 + rSges(3,3) * t325;
t155 = qJD(1) * t171;
t329 = t213 - t155;
t318 = -rSges(5,2) + t410;
t317 = -rSges(4,3) + t410;
t304 = -t324 / 0.2e1;
t303 = t324 / 0.2e1;
t302 = -t323 / 0.2e1;
t301 = t323 / 0.2e1;
t298 = t235 * t378;
t297 = t100 - t383;
t296 = t345 - t401;
t294 = t234 * t410 + t217;
t289 = t213 - t307;
t179 = rSges(2,1) * t236 - rSges(2,2) * t234;
t173 = rSges(2,1) * t234 + rSges(2,2) * t236;
t255 = t335 * t324 + t289;
t28 = (-t171 + t296) * qJD(1) + t255;
t29 = -t214 + (-qJD(3) * t335 + t322) * t236 + (t429 + t346) * qJD(1);
t283 = t234 * t28 - t236 * t29;
t44 = t458 * qJD(1) - t176 * t323 - t214;
t278 = t234 * t43 - t236 * t44;
t277 = -t234 * t75 + t236 * t73;
t263 = -t108 * t234 - t110 * t236;
t257 = t283 * t233;
t132 = t176 * t234;
t81 = qJD(1) * t113 - t214;
t80 = qJD(1) * t336 + t213;
t45 = t263 * qJD(3);
t40 = qJD(1) * t330 + qJDD(1) * t178 - qJDD(2) * t236 + t315;
t39 = t336 * qJDD(1) + (-qJD(1) * t178 - t115) * qJD(1) + t332;
t1 = qJDD(4) * t233 - t345 * t153 - t346 * t152 + (t234 * t394 + t236 * t395 + t322) * qJD(3);
t4 = [-m(2) * (-g(1) * t173 + g(2) * t179) - t437 * t152 / 0.2e1 + t439 * t408 + (((-t436 - t442 + t443) * t234 + (t33 + t393 - t494 + (t98 + t297 + t428) * t234 + t441) * t236) * qJD(3) + t459) * t304 + (-t478 * qJD(3) - t481 * t233 + t480 * t235) * qJD(1) + (t28 * t214 + t29 * (t207 + t289 + t425) + (qJD(3) * t412 - t322) * t236 * t28 + ((t28 * t318 - t29 * t298) * t236 + (t28 * (-qJ(2) + t337) + t29 * t318) * t234) * qJD(1) - (qJD(1) * t296 - t155 + t255 - t28) * t29 + (t3 - g(2)) * (-t234 * t298 + t429 + t431) + (t2 - g(1)) * (t294 + t345)) * m(5) + (-(qJD(1) * t300 + t139 + t329 - t43) * t44 + t43 * (-t236 * t316 + t323 * t391 + t214) + t44 * (-rSges(4,2) * t309 + t312 + t331) + (t317 * t381 + (t43 * (-qJ(2) - t287) + t44 * t317) * t234) * qJD(1) + (-g(2) + t17) * t458 + (-g(1) + t16) * (t110 + t294)) * m(4) + (t80 * t214 + t81 * (t330 + t331) + (t80 * t403 * t236 + (t80 * (-rSges(3,3) - qJ(2)) - t81 * pkin(1)) * t234) * qJD(1) - (qJD(1) * t172 + t329 - t80) * t81 + (t40 - g(2)) * t113 + (t39 - g(1)) * (t234 * t403 + t217 + t386)) * m(3) + (t438 + t440) * t407 + (t445 + t448) * t301 + ((t411 * t98 + (t234 * t297 - t31 + t393) * t234 + ((-t428 - t497) * t236 + t436 + t443) * t236) * qJD(3) + t449 + t454) * t302 + (m(2) * (t173 ^ 2 + t179 ^ 2) + Icges(2,3) + Icges(3,1) + t479) * qJDD(1) + (t446 + t447 + t450) * t303; (-t236 * t40 + 0.2e1 * t39 * t406 - t293) * m(3) + t422 * m(5) + (-t293 + t423) * m(4); t465 * t408 + t464 * t407 + (qJD(1) * t446 + qJD(3) * t426 - qJDD(1) * t437 + t152 * t441 + t153 * t442) * t406 + (qJD(1) * t445 + qJD(3) * t427 + qJDD(1) * t438 + t152 * t443 + t153 * t444) * t236 / 0.2e1 - ((t233 * t475 + t474 * t235) * qJD(3) + (-t233 * t472 + t235 * t473) * qJD(1)) * qJD(1) / 0.2e1 + (t448 * t236 + t447 * t234 + (-t440 * t234 + t439 * t236) * qJD(1)) * qJD(1) / 0.2e1 + (t439 * t234 + t440 * t236) * qJDD(1) / 0.2e1 - t450 * t326 / 0.2e1 + t449 * t325 / 0.2e1 + ((-t453 * t324 + t455) * t234 + ((t457 * t234 + t451) * qJD(3) - t460) * t236) * t304 + ((-t234 * t442 + t441 * t236) * qJD(1) + t426) * t303 + ((t457 * t323 + t455) * t236 + ((-t453 * t236 - t451) * qJD(3) + t460) * t234) * t302 + ((-t234 * t444 + t236 * t443) * qJD(1) + t427) * t301 + (-(t257 + t380) * qJD(4) - (-t28 * t342 + t29 * t343) * qJD(1) - ((-t29 * t337 + t30 * t342) * t236 + (t28 * t337 - t30 * t343) * t234) * qJD(3) - g(1) * t343 - g(3) * (-t233 * t404 + t216 + t228) + t412 * t400 + (-t3 * t335 - t29 * t344 - t1 * t345 + t30 * t395 + (t28 * t335 - t30 * t346) * qJD(1)) * t236 + (t2 * t335 + t28 * t344 - t1 * t346 + t30 * t394 + (t29 * t335 + t30 * t345) * qJD(1)) * t234) * m(5) + (-g(1) * t132 + g(2) * t136 + g(3) * t287 - (t132 * t44 + t136 * t43) * qJD(1) - (t45 * (-t132 * t234 - t136 * t236) - t278 * t287) * qJD(3) + (qJD(3) * t277 - t108 * t152 - t110 * t153) * t263 + t45 * ((-t108 * t236 + t110 * t234) * qJD(1) + t277) - t278 * t149 + ((t234 * t44 + t381) * qJD(1) + t423) * t176) * m(4); (-((t411 + t470) * t380 + t257) * qJD(3) + (qJD(3) * t283 - g(3) + t1) * t233 + (qJD(3) * t30 - t422) * t235) * m(5);];
tau = t4;
