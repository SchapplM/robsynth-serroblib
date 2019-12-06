% Calculate vector of inverse dynamics joint torques for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:51
% DurationCPUTime: 9.40s
% Computational Cost: add. (10466->632), mult. (11433->799), div. (0->0), fcn. (8761->8), ass. (0->343)
t289 = pkin(8) + qJ(4);
t268 = qJ(5) + t289;
t252 = sin(t268);
t253 = cos(t268);
t294 = sin(qJ(1));
t433 = t253 * t294;
t216 = Icges(6,4) * t433;
t434 = t252 * t294;
t295 = cos(qJ(1));
t441 = Icges(6,5) * t295;
t130 = Icges(6,1) * t434 + t216 + t441;
t290 = qJD(4) + qJD(5);
t226 = t290 * t294;
t227 = t290 * t295;
t444 = Icges(6,4) * t252;
t189 = Icges(6,1) * t253 - t444;
t332 = Icges(6,2) * t253 + t444;
t508 = -t332 + t189;
t443 = Icges(6,4) * t253;
t334 = Icges(6,1) * t252 + t443;
t131 = -Icges(6,5) * t294 + t295 * t334;
t187 = -Icges(6,2) * t252 + t443;
t510 = t187 * t295 + t131;
t485 = qJD(1) * t508 - t226 * t510 + t227 * (-Icges(6,2) * t434 + t130 + t216);
t509 = t187 + t334;
t129 = -Icges(6,6) * t294 + t295 * t332;
t511 = -t189 * t295 + t129;
t128 = Icges(6,6) * t295 + t294 * t332;
t512 = -t189 * t294 + t128;
t486 = qJD(1) * t509 - t226 * t511 + t227 * t512;
t514 = t252 * t486 - t253 * t485;
t261 = sin(t289);
t513 = rSges(5,2) * t261;
t342 = rSges(6,1) * t252 + rSges(6,2) * t253;
t262 = cos(t289);
t343 = rSges(5,1) * t261 + rSges(5,2) * t262;
t430 = t262 * t294;
t234 = Icges(5,4) * t430;
t432 = t261 * t294;
t442 = Icges(5,5) * t295;
t139 = Icges(5,1) * t432 + t234 + t442;
t445 = Icges(5,4) * t262;
t335 = Icges(5,1) * t261 + t445;
t140 = -Icges(5,5) * t294 + t295 * t335;
t198 = -Icges(5,2) * t261 + t445;
t167 = t198 * t295;
t305 = t294 * (t140 + t167) - t295 * (-Icges(5,2) * t432 + t139 + t234);
t446 = Icges(5,4) * t261;
t333 = Icges(5,2) * t262 + t446;
t137 = Icges(5,6) * t295 + t294 * t333;
t138 = -Icges(5,6) * t294 + t295 * t333;
t200 = Icges(5,1) * t262 - t446;
t168 = t200 * t294;
t169 = t200 * t295;
t306 = t294 * (t138 - t169) - t295 * (t137 - t168);
t507 = -t306 * t261 + t305 * t262;
t408 = t198 + t335;
t409 = -t333 + t200;
t506 = (t261 * t408 - t262 * t409) * qJD(1);
t505 = -qJ(3) * qJD(1) ^ 2 + qJDD(3);
t224 = rSges(6,1) * t433;
t379 = rSges(6,2) * t434;
t157 = t224 - t379;
t388 = qJD(4) * t295;
t367 = t262 * t388;
t219 = pkin(4) * t367;
t271 = qJD(2) * t295;
t349 = qJD(3) * t294 - t271;
t345 = -t219 + t349;
t291 = sin(pkin(8));
t428 = t291 * t294;
t248 = pkin(3) * t428;
t293 = -pkin(6) - qJ(3);
t424 = qJ(3) + t293;
t175 = -t295 * t424 + t248;
t231 = t295 * pkin(1) + t294 * qJ(2);
t439 = qJ(3) * t295;
t355 = t231 + t439;
t346 = t175 + t355;
t468 = pkin(4) * t261;
t469 = pkin(3) * t291;
t210 = t468 + t469;
t193 = t294 * t210;
t288 = -pkin(7) + t293;
t102 = t193 - t248 + (-t288 + t293) * t295;
t284 = t295 * rSges(6,3);
t132 = rSges(6,1) * t434 + rSges(6,2) * t433 + t284;
t422 = t102 + t132;
t456 = rSges(6,2) * t252;
t459 = rSges(6,1) * t253;
t191 = -t456 + t459;
t436 = t191 * t227;
t43 = -t436 + (t346 + t422) * qJD(1) + t345;
t503 = t43 * (qJD(1) * t157 + t227 * t342);
t330 = Icges(6,5) * t252 + Icges(6,6) * t253;
t127 = -Icges(6,3) * t294 + t295 * t330;
t48 = -t295 * t127 - t129 * t433 - t131 * t434;
t323 = t187 * t253 + t189 * t252;
t185 = Icges(6,5) * t253 - Icges(6,6) * t252;
t437 = t185 * t295;
t61 = t294 * t323 + t437;
t502 = t61 * qJD(1) + t226 * t48;
t225 = t295 * t456;
t158 = -t295 * t459 + t225;
t79 = -t436 + (t294 * t342 + t284) * qJD(1);
t501 = t157 * t226 - t227 * t158 + t295 * t79;
t328 = t129 * t253 + t131 * t252;
t500 = t295 * t328;
t325 = t138 * t262 + t140 * t261;
t497 = t325 * t295;
t392 = qJD(1) * t294;
t391 = qJD(1) * t295;
t371 = t291 * t391;
t404 = pkin(3) * t371 + t293 * t392;
t496 = qJD(1) * (qJ(3) * t392 + t404) + qJDD(1) * t175;
t292 = cos(pkin(8));
t457 = rSges(4,2) * t292;
t162 = rSges(4,1) * t428 + t295 * rSges(4,3) + t294 * t457;
t495 = t162 + t355;
t232 = -rSges(3,2) * t295 + t294 * rSges(3,3);
t494 = qJD(1) * t191;
t331 = Icges(5,5) * t261 + Icges(5,6) * t262;
t136 = -Icges(5,3) * t294 + t295 * t331;
t395 = qJD(1) * t136;
t72 = t138 * t261 - t140 * t262;
t91 = qJD(1) * t137 - qJD(4) * t167;
t93 = -qJD(4) * t169 + (t294 * t335 + t442) * qJD(1);
t493 = qJD(4) * t72 + t261 * t93 + t262 * t91 + t395;
t180 = t333 * qJD(4);
t181 = t335 * qJD(4);
t196 = Icges(5,5) * t262 - Icges(5,6) * t261;
t321 = t198 * t261 - t200 * t262;
t492 = qJD(1) * t196 + qJD(4) * t321 + t180 * t262 + t181 * t261;
t326 = t137 * t261 - t139 * t262;
t135 = Icges(5,3) * t295 + t294 * t331;
t396 = qJD(1) * t135;
t389 = qJD(4) * t294;
t92 = qJD(1) * t138 + t198 * t389;
t94 = qJD(1) * t140 + qJD(4) * t168;
t491 = qJD(4) * t326 - t261 * t94 - t262 * t92 + t396;
t352 = t509 * t290;
t353 = t508 * t290;
t489 = qJD(1) * t185 + t252 * t352 - t253 * t353;
t357 = -qJD(1) * t128 + t510 * t290;
t359 = (t294 * t334 + t441) * qJD(1) + t511 * t290;
t397 = qJD(1) * t127;
t488 = t252 * t359 - t253 * t357 + t397;
t358 = qJD(1) * t129 + t130 * t290 + t187 * t226;
t360 = -qJD(1) * t131 + t512 * t290;
t126 = Icges(6,3) * t295 + t294 * t330;
t398 = qJD(1) * t126;
t487 = t252 * t360 - t253 * t358 + t398;
t484 = t294 ^ 2;
t384 = qJD(1) * qJD(4);
t208 = qJDD(4) * t294 + t295 * t384;
t148 = qJD(5) * t391 + qJDD(5) * t294 + t208;
t483 = t148 / 0.2e1;
t264 = qJDD(4) * t295;
t149 = qJDD(5) * t295 - t290 * t392 + t264;
t482 = t149 / 0.2e1;
t481 = t208 / 0.2e1;
t209 = -t294 * t384 + t264;
t480 = t209 / 0.2e1;
t479 = -t226 / 0.2e1;
t478 = t226 / 0.2e1;
t477 = -t227 / 0.2e1;
t476 = t227 / 0.2e1;
t475 = t294 / 0.2e1;
t474 = -t295 / 0.2e1;
t473 = t295 / 0.2e1;
t472 = rSges(3,2) - pkin(1);
t471 = -rSges(5,3) - pkin(1);
t470 = -rSges(6,3) - pkin(1);
t467 = pkin(4) * t262;
t147 = t342 * t290;
t385 = qJD(1) * qJD(3);
t386 = qJD(1) * qJD(2);
t403 = qJDD(2) * t294 + t295 * t386;
t313 = -0.2e1 * t294 * t385 + t505 * t295 + t403;
t381 = t295 * t469;
t174 = t294 * t424 + t381;
t273 = t295 * qJ(2);
t228 = pkin(1) * t294 - t273;
t440 = qJ(3) * t294;
t356 = -t228 - t440;
t347 = t174 + t356;
t410 = t295 * t210 + t294 * t288;
t101 = t293 * t294 + t381 - t410;
t133 = -t294 * rSges(6,3) + t295 * t342;
t423 = -t101 + t133;
t317 = t347 + t423;
t182 = qJD(1) * t231 - t271;
t245 = t293 * t391;
t413 = t245 - (t248 - t439) * qJD(1) - t182;
t431 = t261 * qJD(4) ^ 2;
t429 = t288 * t295;
t88 = -t219 + t245 + (-t429 + (t210 - t469) * t294) * qJD(1);
t9 = -t147 * t226 + t148 * t191 + (t208 * t262 - t294 * t431) * pkin(4) + (-t79 - t88 + t413) * qJD(1) + t317 * qJDD(1) + t313;
t466 = t9 * t294;
t465 = -qJD(1) / 0.2e1;
t464 = qJD(1) / 0.2e1;
t463 = -pkin(1) - qJ(3);
t376 = t290 * t224 + t342 * t391;
t80 = (-rSges(6,3) * qJD(1) - t290 * t456) * t294 + t376;
t368 = t262 * t389;
t217 = pkin(4) * t368;
t377 = t210 * t391 + t288 * t392 + t217;
t87 = t377 - t404;
t462 = t80 + t87;
t454 = rSges(3,3) * t295;
t258 = qJ(2) * t391;
t270 = qJD(2) * t294;
t402 = t258 + t270;
t378 = qJD(1) * (-pkin(1) * t392 + t402) + qJDD(1) * t231 + t294 * t386;
t302 = qJDD(1) * t439 + t505 * t294 + 0.2e1 * t295 * t385 + t378;
t383 = qJDD(2) * t295;
t300 = t302 - t383;
t10 = t147 * t227 - t149 * t191 + t422 * qJDD(1) + t462 * qJD(1) + (-t209 * t262 + t295 * t431) * pkin(4) + t300 + t496;
t453 = t10 * t295;
t183 = t343 * qJD(4);
t202 = rSges(5,1) * t262 - t513;
t281 = t294 * rSges(5,3);
t146 = t295 * t343 - t281;
t320 = t146 + t347;
t171 = t202 * t295;
t285 = t295 * rSges(5,3);
t97 = -qJD(4) * t171 + (t294 * t343 + t285) * qJD(1);
t28 = -t183 * t389 + t202 * t208 + (-t97 + t413) * qJD(1) + t320 * qJDD(1) + t313;
t451 = t28 * t294;
t145 = rSges(5,1) * t432 + rSges(5,2) * t430 + t285;
t375 = rSges(5,1) * t368 + t343 * t391;
t390 = qJD(4) * t261;
t98 = (-rSges(5,2) * t390 - rSges(5,3) * qJD(1)) * t294 + t375;
t29 = qJD(1) * t98 + qJDD(1) * t145 - t202 * t209 + (qJD(4) * t183 - qJDD(2)) * t295 + t302 + t496;
t450 = t29 * t295;
t269 = qJD(3) * t295;
t400 = t269 + t270;
t318 = t226 * t191 + t217 + t400;
t42 = qJD(1) * t317 + t318;
t449 = t295 * t42;
t173 = t202 * t389;
t58 = qJD(1) * t320 + t173 + t400;
t448 = t295 * t58;
t447 = qJDD(1) / 0.2e1;
t435 = t196 * t295;
t427 = t294 * t127;
t151 = t294 * t185;
t164 = t294 * t196;
t62 = t295 * t323 - t151;
t426 = t62 * qJD(1);
t322 = t198 * t262 + t200 * t261;
t69 = t295 * t322 - t164;
t425 = t69 * qJD(1);
t229 = rSges(3,2) * t294 + t454;
t406 = -t228 + t229;
t177 = t231 + t232;
t405 = rSges(4,1) * t371 + t391 * t457;
t401 = rSges(3,2) * t392 + rSges(3,3) * t391;
t212 = qJD(1) * t228;
t399 = t270 - t212;
t393 = qJD(1) * t331;
t387 = -m(4) - m(5) - m(6);
t382 = -rSges(4,3) + t463;
t239 = pkin(4) * t430;
t47 = t295 * t126 + t128 * t433 + t130 * t434;
t54 = t295 * t135 + t137 * t430 + t139 * t432;
t55 = -t295 * t136 - t138 * t430 - t140 * t432;
t374 = t258 + t400;
t373 = t269 + t399;
t370 = t261 * t389;
t366 = -t392 / 0.2e1;
t365 = t391 / 0.2e1;
t364 = -t389 / 0.2e1;
t363 = t389 / 0.2e1;
t362 = -t388 / 0.2e1;
t361 = t388 / 0.2e1;
t350 = -qJD(1) * t158 - t226 * t342;
t344 = rSges(4,1) * t291 + t457;
t163 = -t294 * rSges(4,3) + t295 * t344;
t348 = t163 + t356;
t233 = rSges(2,1) * t295 - rSges(2,2) * t294;
t230 = rSges(2,1) * t294 + rSges(2,2) * t295;
t327 = t137 * t262 + t139 * t261;
t307 = qJD(1) * t327 + qJD(4) * t164 + t395;
t308 = -qJD(1) * t325 - qJD(4) * t435 + t396;
t341 = (t307 * t294 + t295 * t491) * t295 + (t308 * t294 - t295 * t493) * t294;
t340 = (-t294 * t491 + t307 * t295) * t295 + (t294 * t493 + t308 * t295) * t294;
t339 = t294 * t55 + t295 * t54;
t122 = t294 * t135;
t56 = -t327 * t295 + t122;
t57 = -t136 * t294 + t497;
t338 = t294 * t57 + t295 * t56;
t59 = -t202 * t388 + (t145 + t346) * qJD(1) + t349;
t337 = t294 * t58 - t295 * t59;
t336 = -t294 * t98 + t295 * t97;
t329 = t128 * t253 + t130 * t252;
t66 = t129 * t252 - t131 * t253;
t324 = -t145 * t294 - t146 * t295;
t319 = t47 + t427;
t316 = t343 + t469;
t312 = -qJD(1) * t330 + t151 * t227 - t226 * t437;
t116 = t294 * t126;
t49 = -t295 * t329 + t116;
t310 = -qJD(1) * t328 - t290 * t437 + t398;
t309 = qJD(1) * t329 + t151 * t290 + t397;
t304 = qJD(1) * t323 - t330 * t290;
t303 = t322 * qJD(1) - t331 * qJD(4);
t11 = t309 * t294 + t295 * t487;
t12 = t310 * t294 - t295 * t488;
t13 = -t294 * t487 + t309 * t295;
t14 = t294 * t488 + t310 * t295;
t22 = t227 * t47 + t502;
t50 = -t427 + t500;
t23 = t226 * t50 + t227 * t49 - t426;
t30 = t304 * t294 + t295 * t489;
t31 = -t294 * t489 + t304 * t295;
t32 = -t252 * t358 - t253 * t360;
t33 = t252 * t357 + t253 * t359;
t65 = -t128 * t252 + t130 * t253;
t301 = (qJD(1) * t30 - qJDD(1) * t62 + t11 * t227 + t12 * t226 + t148 * t50 + t149 * t49) * t475 + (-t252 * t485 - t253 * t486) * t465 + (qJD(1) * t31 + qJDD(1) * t61 + t13 * t227 + t14 * t226 + t148 * t48 + t149 * t47) * t473 + t22 * t366 + t23 * t365 + (t294 * t50 + t295 * t49) * t483 + (t294 * t48 + t295 * t47) * t482 + (t11 * t295 + t12 * t294 + (-t294 * t49 + t295 * t50) * qJD(1)) * t478 + (t13 * t295 + t14 * t294 + (-t294 * t47 + t295 * t48) * qJD(1)) * t476 + (t294 * t66 + t295 * t65) * t447 + (t294 * t33 + t295 * t32 + (-t294 * t65 + t295 * t66) * qJD(1)) * t464 + (t312 * t294 + t514 * t295) * t479 + (-t514 * t294 + t312 * t295) * t477;
t170 = t202 * t294;
t160 = qJD(1) * t174;
t121 = qJD(1) * t177 - t271;
t120 = qJD(1) * t406 + t270;
t119 = t295 * t133;
t96 = qJD(1) * t495 + t349;
t95 = qJD(1) * t348 + t400;
t81 = t324 * qJD(4);
t68 = t294 * t322 + t435;
t67 = t68 * qJD(1);
t64 = qJD(1) * t401 + qJDD(1) * t232 + t378 - t383;
t63 = t406 * qJDD(1) + (-qJD(1) * t232 - t182) * qJD(1) + t403;
t46 = qJDD(1) * t162 + qJD(1) * (-rSges(4,3) * t392 + t405) + t300;
t45 = t348 * qJDD(1) + (-qJD(1) * t162 - t182) * qJD(1) + t313;
t41 = -t132 * t226 - t133 * t227 + (t101 * t295 - t102 * t294) * qJD(4);
t37 = t325 * qJD(4) - t261 * t91 + t262 * t93;
t36 = -qJD(4) * t327 - t261 * t92 + t262 * t94;
t35 = -t294 * t492 + t303 * t295;
t34 = t303 * t294 + t295 * t492;
t25 = qJD(4) * t338 - t425;
t24 = qJD(4) * t339 + t67;
t8 = t101 * t209 - t102 * t208 - t132 * t148 - t133 * t149 - t226 * t80 + t227 * t79 + (-t294 * t87 + t295 * t88) * qJD(4);
t1 = [-m(2) * (-g(1) * t230 + g(2) * t233) - t208 * t69 / 0.2e1 - t148 * t62 / 0.2e1 + (t67 + ((-t56 + t122 + t55) * t294 + (t57 - t497 + (t136 - t327) * t294 + t54) * t295) * qJD(4)) * t364 + ((t319 + t50 - t500) * t227 + t502) * t479 + t72 * t481 + t66 * t483 + (t65 + t61) * t482 + (-t326 + t68) * t480 + (t426 + (t328 * t294 - t116 + t48) * t227 + (t319 - t47) * t226 + ((t127 + t329) * t227 - t328 * t226) * t295 + t23) * t477 + (t32 + t31) * t476 + (t425 + (t136 * t484 + (-t122 + t55 + (t136 + t327) * t295) * t295) * qJD(4) + t25) * t362 + (t36 + t35) * t361 + (-qJD(4) * t322 + t180 * t261 - t181 * t262 - t252 * t353 - t253 * t352) * qJD(1) + (t33 + t30 + t22) * t478 + (t37 + t34 + t24) * t363 + (-(t160 - t212 - t42 + (t423 - t440) * qJD(1) + t318) * t43 + t42 * (-t158 * t290 - t345) + t43 * (-t290 * t379 + t374 + t376 + t377) + ((t288 + t470) * t449 + (t42 * (-qJ(2) - t210 - t342) + t43 * t470) * t294) * qJD(1) + (-g(2) + t10) * (t132 + t193 + t231 - t429) + (-g(1) + t9) * (t133 - t228 + t410)) * m(6) + (-(t160 + t173 - t58 + (t146 - t440) * qJD(1) + t373) * t59 + t58 * (rSges(5,1) * t367 - t388 * t513 + t245 - t349) + t59 * (-rSges(5,2) * t370 + t374 + t375 + t404) + (t471 * t448 + (t58 * (-qJ(2) - t316) + t59 * t471) * t294) * qJD(1) + (-g(1) + t28) * (t273 - t281 + (-pkin(1) + t293) * t294 + t316 * t295) + (-g(2) + t29) * (-t293 * t295 + t145 + t231 + t248)) * m(5) + (-t95 * t349 + t96 * (t374 + t405) + (t95 * t382 * t295 + (t95 * (-qJ(2) - t344) + t96 * t382) * t294) * qJD(1) - (-t95 + (t163 - t440) * qJD(1) + t373) * t96 + (t46 - g(2)) * t495 + (t45 - g(1)) * (t294 * t463 + t163 + t273)) * m(4) + (-(qJD(1) * t229 - t120 + t399) * t121 + t120 * t271 + t121 * (t401 + t402) + (t120 * t472 * t295 + (t120 * (-rSges(3,3) - qJ(2)) - t121 * pkin(1)) * t294) * qJD(1) + (-g(2) + t64) * t177 + (-g(1) + t63) * (t294 * t472 + t273 + t454)) * m(3) + (-t187 * t252 + t189 * t253 - t321 + m(2) * (t230 ^ 2 + t233 ^ 2) + Icges(4,1) * t292 ^ 2 + (-0.2e1 * Icges(4,4) * t292 + Icges(4,2) * t291) * t291 + Icges(2,3) + Icges(3,1)) * qJDD(1); (-m(3) + t387) * (g(1) * t294 - g(2) * t295) + 0.2e1 * (-t453 / 0.2e1 + t466 / 0.2e1) * m(6) + 0.2e1 * (t451 / 0.2e1 - t450 / 0.2e1) * m(5) + 0.2e1 * (t45 * t475 + t46 * t474) * m(4) + 0.2e1 * (t474 * t64 + t475 * t63) * m(3); t387 * (g(1) * t295 + g(2) * t294) + m(4) * (t294 * t46 + t295 * t45) + m(5) * (t28 * t295 + t29 * t294) + m(6) * (t10 * t294 + t295 * t9); t301 + ((-t54 * t294 + t55 * t295) * qJD(1) + t340) * t361 + ((-t56 * t294 + t57 * t295) * qJD(1) + t341) * t363 + t338 * t481 + t339 * t480 + (t294 * t72 - t295 * t326) * t447 + (qJD(1) * t34 + qJD(4) * t341 - qJDD(1) * t69 + t208 * t57 + t209 * t56) * t475 + ((-t261 * t409 - t262 * t408) * qJD(1) + (t261 * t305 + t262 * t306) * qJD(4)) * t465 + ((t164 * t388 - t393) * t295 + (-t506 + (-t295 * t435 - t507) * qJD(4)) * t294) * t362 + (t294 * t37 + t295 * t36 + (t294 * t326 + t72 * t295) * qJD(1)) * t464 + ((-t389 * t435 - t393) * t294 + (t506 + (t294 * t164 + t507) * qJD(4)) * t295) * t364 + (qJD(1) * t35 + qJD(4) * t340 + qJDD(1) * t68 + t208 * t55 + t209 * t54) * t473 + t24 * t366 + t25 * t365 + (-g(1) * (t157 + t239) - g(2) * (t225 + (-t459 - t467) * t295) - g(3) * (-t342 - t468) + t9 * t239 - t8 * t119 + (t10 * (-t191 - t467) + t43 * t147 + t8 * t101 + t42 * t494) * t295 + (t9 * t191 + t42 * (-pkin(4) * t390 - t147) - t8 * t422 + t43 * t494) * t294 - t42 * (-pkin(4) * t370 + t350) - t503 + ((-qJD(1) * t422 + t88) * t295 + (t423 * qJD(1) - t462) * t294 - (-t295 ^ 2 - t484) * qJD(4) * t467 + t501) * t41) * m(6) + (-g(1) * t170 + g(2) * t171 + g(3) * t343 - (t170 * t59 + t171 * t58) * qJD(1) - (t81 * (-t170 * t294 - t171 * t295) - t337 * t343) * qJD(4) + (qJD(4) * t336 - t145 * t208 - t146 * t209) * t324 + t81 * ((-t145 * t295 + t146 * t294) * qJD(1) + t336) - t337 * t183 + (t451 - t450 + (t294 * t59 + t448) * qJD(1)) * t202) * m(5); t301 + (-g(1) * t157 - g(2) * t158 + g(3) * t342 + t8 * (-t132 * t294 - t119) - (t294 * t42 - t295 * t43) * t147 + (-t453 + t466 + (t294 * t43 + t449) * qJD(1)) * t191 - t42 * t350 - t503 + (-t294 * t80 + (-t132 * t295 + t133 * t294) * qJD(1) + t501) * t41) * m(6);];
tau = t1;
