% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:47:54
% EndTime: 2018-11-23 17:48:21
% DurationCPUTime: 27.57s
% Computational Cost: add. (17932->867), mult. (47425->1172), div. (0->0), fcn. (37170->10), ass. (0->357)
t337 = sin(qJ(2));
t332 = sin(pkin(6));
t396 = qJD(1) * t332;
t380 = t337 * t396;
t334 = cos(pkin(6));
t340 = cos(qJ(2));
t451 = pkin(1) * t340;
t387 = t334 * t451;
t288 = -pkin(8) * t380 + qJD(1) * t387;
t348 = (pkin(2) * t337 - pkin(9) * t340) * t332;
t289 = qJD(1) * t348;
t336 = sin(qJ(3));
t339 = cos(qJ(3));
t206 = t339 * t288 + t336 * t289;
t185 = qJ(4) * t380 + t206;
t452 = pkin(1) * t337;
t325 = t334 * t452;
t359 = pkin(3) * t336 - qJ(4) * t339;
t409 = t332 * t340;
t212 = (t325 + (pkin(8) + t359) * t409) * qJD(1);
t331 = sin(pkin(11));
t333 = cos(pkin(11));
t128 = -t331 * t185 + t333 * t212;
t294 = qJD(3) * t359 - qJD(4) * t336;
t392 = qJD(3) * t336;
t386 = pkin(9) * t392;
t232 = t333 * t294 + t331 * t386;
t532 = t232 - t128;
t406 = t339 * t340;
t253 = (t331 * t337 + t333 * t406) * t396;
t379 = t340 * t396;
t366 = t336 * t379;
t407 = t333 * t339;
t531 = -pkin(4) * t366 + t253 * pkin(10) + (pkin(4) * t336 - pkin(10) * t407) * qJD(3) + t532;
t129 = t333 * t185 + t331 * t212;
t252 = (-t331 * t406 + t333 * t337) * t396;
t279 = t331 * t294;
t408 = t333 * t336;
t411 = t331 * t339;
t530 = pkin(10) * t252 + t129 - t279 - (-pkin(9) * t408 - pkin(10) * t411) * qJD(3);
t508 = Ifges(6,1) + Ifges(7,1);
t507 = -Ifges(6,4) + Ifges(7,5);
t506 = Ifges(7,4) + Ifges(6,5);
t335 = sin(qJ(5));
t338 = cos(qJ(5));
t172 = -t338 * t252 + t253 * t335;
t306 = t331 * t338 + t333 * t335;
t388 = qJD(5) * t338;
t389 = qJD(5) * t336;
t391 = qJD(3) * t339;
t413 = t331 * t335;
t210 = t306 * t391 + t388 * t408 - t389 * t413;
t402 = t172 - t210;
t173 = t252 * t335 + t253 * t338;
t305 = -t338 * t333 + t413;
t209 = -t305 * t391 - t306 * t389;
t401 = t173 - t209;
t529 = t366 - t392;
t309 = -pkin(3) * t339 - qJ(4) * t336 - pkin(2);
t303 = t333 * t309;
t228 = -pkin(10) * t408 + t303 + (-pkin(9) * t331 - pkin(4)) * t339;
t267 = pkin(9) * t407 + t331 * t309;
t412 = t331 * t336;
t240 = -pkin(10) * t412 + t267;
t390 = qJD(5) * t335;
t523 = t228 * t388 - t240 * t390 + t335 * t531 - t530 * t338;
t493 = t335 * t228 + t338 * t240;
t522 = -qJD(5) * t493 + t530 * t335 + t338 * t531;
t205 = -t336 * t288 + t289 * t339;
t187 = -pkin(3) * t380 - t205;
t518 = pkin(4) * t252 - t187 + (pkin(4) * t331 + pkin(9)) * t391;
t368 = qJD(1) * t334 + qJD(2);
t492 = -t336 * t380 + t339 * t368;
t174 = t306 * t492;
t296 = t306 * qJD(5);
t400 = t174 - t296;
t175 = t305 * t492;
t295 = t305 * qJD(5);
t399 = -t175 + t295;
t269 = t336 * t368 + t339 * t380;
t395 = qJD(2) * t332;
t369 = qJD(1) * t395;
t364 = t340 * t369;
t225 = qJD(3) * t269 + t336 * t364;
t262 = qJD(5) - t492;
t224 = qJD(3) * t492 + t339 * t364;
t365 = t337 * t369;
t184 = t224 * t333 + t331 * t365;
t324 = pkin(8) * t409;
t250 = qJD(2) * pkin(9) + (t324 + (pkin(9) + t452) * t334) * qJD(1);
t283 = (-pkin(2) * t340 - pkin(9) * t337 - pkin(1)) * t332;
t260 = qJD(1) * t283;
t290 = qJD(2) * t348;
t276 = qJD(1) * t290;
t410 = t332 * t337;
t321 = pkin(8) * t410;
t299 = -t321 + t387;
t292 = t299 * qJD(2);
t277 = qJD(1) * t292;
t120 = -t250 * t392 + t260 * t391 + t336 * t276 + t339 * t277;
t315 = qJD(3) - t379;
t101 = qJ(4) * t365 + qJD(4) * t315 + t120;
t397 = t324 + t325;
t293 = t397 * qJD(2);
t278 = qJD(1) * t293;
t107 = t225 * pkin(3) - t224 * qJ(4) - t269 * qJD(4) + t278;
t44 = -t101 * t331 + t333 * t107;
t24 = pkin(4) * t225 - pkin(10) * t184 + t44;
t183 = -t224 * t331 + t333 * t365;
t45 = t333 * t101 + t331 * t107;
t32 = pkin(10) * t183 + t45;
t214 = t269 * t333 + t315 * t331;
t249 = -pkin(2) * t368 - t288;
t153 = -pkin(3) * t492 - t269 * qJ(4) + t249;
t170 = t339 * t250 + t336 * t260;
t158 = qJ(4) * t315 + t170;
t89 = t333 * t153 - t158 * t331;
t52 = -pkin(4) * t492 - pkin(10) * t214 + t89;
t350 = t269 * t331 - t315 * t333;
t90 = t331 * t153 + t333 * t158;
t72 = -pkin(10) * t350 + t90;
t3 = t335 * t24 + t338 * t32 + t52 * t388 - t390 * t72;
t1 = qJ(6) * t225 + qJD(6) * t262 + t3;
t467 = t225 / 0.2e1;
t513 = t338 * t214 - t335 * t350;
t69 = qJD(5) * t513 - t338 * t183 + t335 * t184;
t482 = t69 / 0.2e1;
t483 = -t69 / 0.2e1;
t140 = t214 * t335 + t338 * t350;
t68 = -qJD(5) * t140 + t335 * t183 + t338 * t184;
t484 = t68 / 0.2e1;
t510 = -t225 / 0.2e1;
t121 = -t250 * t391 - t260 * t392 + t276 * t339 - t336 * t277;
t108 = -pkin(3) * t365 - t121;
t76 = -pkin(4) * t183 + t108;
t9 = pkin(5) * t69 - qJ(6) * t68 - qJD(6) * t513 + t76;
t487 = mrSges(6,1) * t76 + mrSges(7,1) * t9 - mrSges(7,2) * t1 - mrSges(6,3) * t3 - t68 * Ifges(6,4) / 0.2e1 + Ifges(6,6) * t510 + Ifges(7,6) * t467 + 0.2e1 * Ifges(7,3) * t482 + (-t483 + t482) * Ifges(6,2) + (t507 + Ifges(7,5)) * t484;
t20 = t335 * t52 + t338 * t72;
t4 = -qJD(5) * t20 + t24 * t338 - t32 * t335;
t2 = -pkin(5) * t225 - t4;
t528 = -mrSges(6,2) * t76 - mrSges(7,2) * t2 + mrSges(6,3) * t4 + mrSges(7,3) * t9 - Ifges(6,4) * t483 - Ifges(7,5) * t482 - t506 * t467 - t508 * t484;
t502 = t506 * t225 + t507 * t69 + t508 * t68;
t527 = t502 / 0.2e1;
t505 = Ifges(6,6) - Ifges(7,6);
t136 = Ifges(7,5) * t513;
t62 = Ifges(7,6) * t262 + Ifges(7,3) * t140 + t136;
t434 = Ifges(6,4) * t513;
t65 = -Ifges(6,2) * t140 + Ifges(6,6) * t262 + t434;
t526 = t65 - t62;
t137 = Ifges(6,4) * t140;
t431 = Ifges(7,5) * t140;
t501 = t506 * t262 + t508 * t513 - t137 + t431;
t525 = -qJ(6) * t529 - qJD(6) * t339 + t523;
t524 = pkin(5) * t529 - t522;
t521 = Ifges(4,5) * t224;
t285 = t305 * t336;
t519 = -pkin(5) * t402 + qJ(6) * t401 + qJD(6) * t285 + t518;
t443 = pkin(10) + qJ(4);
t311 = t443 * t331;
t312 = t443 * t333;
t349 = -t338 * t311 - t312 * t335;
t169 = -t336 * t250 + t339 * t260;
t194 = pkin(3) * t269 - qJ(4) * t492;
t115 = -t169 * t331 + t333 * t194;
t85 = -pkin(10) * t333 * t492 + pkin(4) * t269 + t115;
t116 = t333 * t169 + t331 * t194;
t414 = t492 * t331;
t98 = -pkin(10) * t414 + t116;
t497 = -qJD(4) * t305 + qJD(5) * t349 - t335 * t85 - t338 * t98;
t245 = -t311 * t335 + t312 * t338;
t496 = -qJD(4) * t306 - qJD(5) * t245 + t335 * t98 - t338 * t85;
t517 = -t121 * mrSges(4,1) + t120 * mrSges(4,2);
t418 = t315 * Ifges(4,5);
t516 = -t418 / 0.2e1 - t249 * mrSges(4,2);
t19 = -t335 * t72 + t338 * t52;
t494 = qJD(6) - t19;
t17 = -pkin(5) * t262 + t494;
t18 = qJ(6) * t262 + t20;
t347 = Ifges(5,6) * t350;
t417 = t315 * Ifges(4,6);
t515 = -t249 * mrSges(4,1) - t89 * mrSges(5,1) - t19 * mrSges(6,1) + t17 * mrSges(7,1) + t90 * mrSges(5,2) + t20 * mrSges(6,2) - t18 * mrSges(7,3) + t347 / 0.2e1 + t417 / 0.2e1;
t514 = -t278 * mrSges(3,1) - t277 * mrSges(3,2);
t509 = -t350 / 0.2e1;
t376 = -Ifges(3,6) * qJD(2) / 0.2e1;
t504 = Ifges(6,3) + Ifges(7,2);
t12 = Ifges(6,5) * t68 - Ifges(6,6) * t69 + Ifges(6,3) * t225;
t13 = Ifges(7,4) * t68 + Ifges(7,2) * t225 + Ifges(7,6) * t69;
t503 = t12 + t13;
t145 = pkin(4) * t414 + t170;
t499 = -pkin(5) * t400 + qJ(6) * t399 - qJD(6) * t306 - t145;
t498 = -qJ(6) * t269 + t497;
t495 = pkin(5) * t269 - t496;
t432 = Ifges(5,5) * t214;
t123 = -Ifges(5,3) * t492 - t347 + t432;
t63 = Ifges(6,5) * t513 - t140 * Ifges(6,6) + t262 * Ifges(6,3);
t64 = Ifges(7,4) * t513 + t262 * Ifges(7,2) + t140 * Ifges(7,6);
t491 = t123 + t64 + t63;
t297 = -t334 * t339 + t336 * t410;
t393 = qJD(2) * t340;
t377 = t332 * t393;
t238 = -qJD(3) * t297 + t339 * t377;
t394 = qJD(2) * t337;
t378 = t332 * t394;
t204 = t238 * t333 + t331 * t378;
t298 = t334 * t336 + t339 * t410;
t239 = qJD(3) * t298 + t336 * t377;
t282 = pkin(9) * t334 + t397;
t134 = -t282 * t392 + t283 * t391 + t336 * t290 + t339 * t292;
t122 = (qJ(4) * t394 - qJD(4) * t340) * t332 + t134;
t131 = t239 * pkin(3) - t238 * qJ(4) - t298 * qJD(4) + t293;
t53 = -t122 * t331 + t333 * t131;
t36 = pkin(4) * t239 - pkin(10) * t204 + t53;
t203 = -t238 * t331 + t333 * t378;
t54 = t333 * t122 + t331 * t131;
t43 = pkin(10) * t203 + t54;
t281 = t321 + (-pkin(2) - t451) * t334;
t178 = t297 * pkin(3) - t298 * qJ(4) + t281;
t197 = t339 * t282 + t336 * t283;
t180 = -qJ(4) * t409 + t197;
t109 = t333 * t178 - t180 * t331;
t237 = t333 * t298 - t331 * t409;
t81 = pkin(4) * t297 - pkin(10) * t237 + t109;
t110 = t331 * t178 + t333 * t180;
t236 = -t331 * t298 - t333 * t409;
t91 = pkin(10) * t236 + t110;
t442 = t335 * t81 + t338 * t91;
t8 = -qJD(5) * t442 - t335 * t43 + t338 * t36;
t124 = Ifges(5,4) * t214 - t350 * Ifges(5,2) - Ifges(5,6) * t492;
t125 = Ifges(5,1) * t214 - t350 * Ifges(5,4) - Ifges(5,5) * t492;
t155 = -t315 * pkin(3) + qJD(4) - t169;
t435 = Ifges(5,4) * t333;
t361 = -Ifges(5,2) * t331 + t435;
t363 = mrSges(5,1) * t331 + mrSges(5,2) * t333;
t453 = t333 / 0.2e1;
t454 = -t331 / 0.2e1;
t489 = (-t331 * t90 - t333 * t89) * mrSges(5,3) + t155 * t363 + t361 * t509 + t125 * t453 + t124 * t454;
t382 = -Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1;
t383 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t384 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t385 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t420 = t269 * Ifges(4,4);
t165 = Ifges(4,2) * t492 + t417 + t420;
t472 = -t165 / 0.2e1;
t488 = t140 * t383 + t513 * t385 + t262 * t384 + t492 * t382 - t420 / 0.2e1 + t432 / 0.2e1 + t63 / 0.2e1 + t64 / 0.2e1 + t472 + t123 / 0.2e1 - t170 * mrSges(4,3) - t515;
t94 = t184 * Ifges(5,1) + t183 * Ifges(5,4) + t225 * Ifges(5,5);
t480 = t94 / 0.2e1;
t479 = pkin(1) * mrSges(3,1);
t478 = pkin(1) * mrSges(3,2);
t477 = -t140 / 0.2e1;
t476 = t140 / 0.2e1;
t475 = -t513 / 0.2e1;
t474 = t513 / 0.2e1;
t471 = t183 / 0.2e1;
t470 = t184 / 0.2e1;
t469 = -t214 / 0.2e1;
t468 = t214 / 0.2e1;
t465 = -t262 / 0.2e1;
t464 = t262 / 0.2e1;
t463 = t492 / 0.2e1;
t462 = -t492 / 0.2e1;
t461 = -t269 / 0.2e1;
t460 = t269 / 0.2e1;
t458 = -t297 / 0.2e1;
t456 = t298 / 0.2e1;
t330 = t336 * pkin(9);
t444 = qJD(2) / 0.2e1;
t441 = mrSges(4,3) * t492;
t440 = mrSges(4,3) * t269;
t439 = mrSges(6,3) * t140;
t438 = mrSges(6,3) * t513;
t437 = Ifges(3,4) * t337;
t261 = Ifges(4,4) * t492;
t436 = Ifges(5,4) * t331;
t433 = Ifges(3,5) * t334;
t430 = Ifges(3,6) * t334;
t427 = t183 * Ifges(5,6);
t426 = t184 * Ifges(5,5);
t425 = t224 * Ifges(4,1);
t424 = t224 * Ifges(4,4);
t423 = t225 * Ifges(4,4);
t421 = t269 * Ifges(4,1);
t111 = -mrSges(7,2) * t140 + mrSges(7,3) * t262;
t112 = -mrSges(6,2) * t262 - t439;
t405 = -t111 - t112;
t113 = mrSges(6,1) * t262 - t438;
t114 = -mrSges(7,1) * t262 + mrSges(7,2) * t513;
t404 = t113 - t114;
t146 = mrSges(5,1) * t350 + t214 * mrSges(5,2);
t227 = mrSges(4,1) * t315 - t440;
t403 = t146 - t227;
t398 = -mrSges(3,1) * t368 - mrSges(4,1) * t492 + mrSges(4,2) * t269 + mrSges(3,3) * t380;
t307 = pkin(4) * t412 + t330;
t381 = -Ifges(4,6) * t225 + Ifges(4,3) * t365 + t521;
t328 = -pkin(4) * t333 - pkin(3);
t22 = t69 * mrSges(6,1) + t68 * mrSges(6,2);
t21 = t69 * mrSges(7,1) - t68 * mrSges(7,3);
t47 = -t225 * mrSges(7,1) + t68 * mrSges(7,2);
t119 = -t183 * mrSges(5,1) + t184 * mrSges(5,2);
t196 = -t336 * t282 + t283 * t339;
t181 = pkin(3) * t409 - t196;
t362 = Ifges(5,1) * t333 - t436;
t360 = Ifges(5,5) * t333 - Ifges(5,6) * t331;
t357 = -t331 * t44 + t333 * t45;
t29 = -t335 * t91 + t338 * t81;
t156 = t228 * t338 - t240 * t335;
t351 = t338 * t236 - t237 * t335;
t160 = t236 * t335 + t237 * t338;
t135 = -t282 * t391 - t283 * t392 + t290 * t339 - t336 * t292;
t7 = t335 * t36 + t338 * t43 + t81 * t388 - t390 * t91;
t144 = -pkin(4) * t236 + t181;
t316 = Ifges(3,4) * t379;
t345 = -t288 * mrSges(3,3) + Ifges(3,1) * t380 / 0.2e1 + t316 / 0.2e1 + (t368 / 0.2e1 + t444) * Ifges(3,5);
t344 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t130 = -pkin(3) * t378 - t135;
t95 = -pkin(4) * t203 + t130;
t118 = pkin(4) * t350 + t155;
t166 = t261 + t418 + t421;
t343 = -t261 / 0.2e1 - t421 / 0.2e1 - t166 / 0.2e1 + t169 * mrSges(4,3) + t516;
t291 = t397 * qJD(1);
t342 = t169 * mrSges(4,1) + t315 * Ifges(4,3) + t269 * Ifges(4,5) + t492 * Ifges(4,6) + t376 - (t430 + (Ifges(3,2) * t340 + t437) * t332) * qJD(1) / 0.2e1 - t170 * mrSges(4,2) - t291 * mrSges(3,3);
t314 = Ifges(3,5) * t364;
t287 = -mrSges(3,2) * t368 + mrSges(3,3) * t379;
t284 = t306 * t336;
t266 = -pkin(9) * t411 + t303;
t233 = -t333 * t386 + t279;
t226 = -mrSges(4,2) * t315 + t441;
t223 = pkin(5) * t305 - qJ(6) * t306 + t328;
t191 = -mrSges(4,2) * t365 - mrSges(4,3) * t225;
t190 = mrSges(4,1) * t365 - mrSges(4,3) * t224;
t182 = pkin(5) * t284 + qJ(6) * t285 + t307;
t162 = -mrSges(5,1) * t492 - mrSges(5,3) * t214;
t161 = mrSges(5,2) * t492 - mrSges(5,3) * t350;
t152 = pkin(5) * t339 - t156;
t151 = -qJ(6) * t339 + t493;
t150 = mrSges(4,1) * t225 + mrSges(4,2) * t224;
t139 = Ifges(4,5) * t365 - t423 + t425;
t138 = -t225 * Ifges(4,2) + Ifges(4,6) * t365 + t424;
t133 = mrSges(5,1) * t225 - mrSges(5,3) * t184;
t132 = -mrSges(5,2) * t225 + mrSges(5,3) * t183;
t93 = t184 * Ifges(5,4) + t183 * Ifges(5,2) + t225 * Ifges(5,6);
t92 = t225 * Ifges(5,3) + t426 + t427;
t84 = qJD(5) * t160 - t338 * t203 + t204 * t335;
t83 = qJD(5) * t351 + t203 * t335 + t204 * t338;
t75 = mrSges(6,1) * t140 + mrSges(6,2) * t513;
t74 = mrSges(7,1) * t140 - mrSges(7,3) * t513;
t73 = pkin(5) * t513 + qJ(6) * t140;
t50 = -pkin(5) * t351 - qJ(6) * t160 + t144;
t49 = -mrSges(7,2) * t69 + mrSges(7,3) * t225;
t48 = -mrSges(6,2) * t225 - mrSges(6,3) * t69;
t46 = mrSges(6,1) * t225 - mrSges(6,3) * t68;
t37 = t140 * pkin(5) - qJ(6) * t513 + t118;
t26 = -pkin(5) * t297 - t29;
t25 = qJ(6) * t297 + t442;
t10 = pkin(5) * t84 - qJ(6) * t83 - qJD(6) * t160 + t95;
t6 = -pkin(5) * t239 - t8;
t5 = qJ(6) * t239 + qJD(6) * t297 + t7;
t11 = [(Ifges(5,5) * t204 + Ifges(5,6) * t203 + Ifges(5,3) * t239) * t462 + (Ifges(4,4) * t238 - Ifges(4,2) * t239) * t463 + m(6) * (t118 * t95 + t144 * t76 + t19 * t8 + t20 * t7 + t29 * t4 + t3 * t442) + t442 * t48 + m(3) * (t277 * t397 - t278 * t299 - t288 * t293 + t291 * t292) + ((-t299 * mrSges(3,3) + t433 + (-0.2e1 * t478 + 0.3e1 / 0.2e1 * Ifges(3,4) * t340) * t332) * t393 + (Ifges(4,5) * t456 + Ifges(4,6) * t458 - 0.3e1 / 0.2e1 * t430 - t397 * mrSges(3,3) + (-0.2e1 * t479 - 0.3e1 / 0.2e1 * t437 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1) * t340) * t332) * t394) * t396 + (mrSges(4,1) * t297 + mrSges(4,2) * t298 + mrSges(3,3) * t410) * t278 + (t1 * t297 - t160 * t9 + t18 * t239 - t37 * t83) * mrSges(7,3) + (-t120 * t297 - t121 * t298 - t169 * t238 - t170 * t239) * mrSges(4,3) + (t118 * t83 + t160 * t76 - t20 * t239 - t297 * t3) * mrSges(6,2) + (Ifges(4,4) * t298 - Ifges(4,2) * t297) * t510 + t249 * (mrSges(4,1) * t239 + mrSges(4,2) * t238) + t90 * (-mrSges(5,2) * t239 + mrSges(5,3) * t203) + t89 * (mrSges(5,1) * t239 - mrSges(5,3) * t204) + t19 * (mrSges(6,1) * t239 - mrSges(6,3) * t83) + t17 * (-mrSges(7,1) * t239 + mrSges(7,2) * t83) + t236 * t93 / 0.2e1 + t108 * (-mrSges(5,1) * t236 + mrSges(5,2) * t237) + t238 * t166 / 0.2e1 + t134 * t226 + t135 * t227 + t398 * t293 + (Ifges(7,5) * t83 + Ifges(7,6) * t239) * t476 + (Ifges(7,5) * t160 + Ifges(7,6) * t297) * t482 + (-t521 / 0.2e1 + t277 * mrSges(3,3) - t381 / 0.2e1 - Ifges(4,6) * t510 + t517) * t409 - t487 * t351 + (Ifges(5,5) * t237 + Ifges(5,6) * t236 + t506 * t160 + t505 * t351 + (Ifges(5,3) + t504) * t297) * t467 + m(4) * (t120 * t197 + t121 * t196 + t134 * t170 + t135 * t169 + t249 * t293 + t278 * t281) + m(5) * (t108 * t181 + t109 * t44 + t110 * t45 + t130 * t155 + t53 * t89 + t54 * t90) + m(7) * (t1 * t25 + t10 * t37 + t17 * t6 + t18 * t5 + t2 * t26 + t50 * t9) + t155 * (-mrSges(5,1) * t203 + mrSges(5,2) * t204) + t204 * t125 / 0.2e1 + t203 * t124 / 0.2e1 + t196 * t190 + t197 * t191 + t181 * t119 + t54 * t161 + t53 * t162 + t130 * t146 + t144 * t22 + t110 * t132 + t109 * t133 + t5 * t111 + t7 * t112 + t8 * t113 + t6 * t114 + t95 * t75 + t10 * t74 + t50 * t21 + t29 * t46 + t26 * t47 + t25 * t49 + t224 * (Ifges(4,1) * t298 - Ifges(4,4) * t297) / 0.2e1 + t315 * (Ifges(4,5) * t238 - Ifges(4,6) * t239) / 0.2e1 + t139 * t456 + t138 * t458 + (Ifges(4,1) * t238 - Ifges(4,4) * t239) * t460 + (Ifges(5,1) * t204 + Ifges(5,4) * t203 + Ifges(5,5) * t239) * t468 + (Ifges(5,1) * t237 + Ifges(5,4) * t236 + Ifges(5,5) * t297) * t470 + (Ifges(5,4) * t237 + Ifges(5,2) * t236 + Ifges(5,6) * t297) * t471 + t239 * t472 + t491 * t239 / 0.2e1 + t281 * t150 + t292 * t287 + t45 * (-mrSges(5,2) * t297 + mrSges(5,3) * t236) + t44 * (mrSges(5,1) * t297 - mrSges(5,3) * t237) + t4 * (mrSges(6,1) * t297 - mrSges(6,3) * t160) + t2 * (-mrSges(7,1) * t297 + mrSges(7,2) * t160) + t501 * t83 / 0.2e1 + (t92 + t503) * t297 / 0.2e1 + (t239 * t504 + t506 * t83) * t464 + (-t18 * mrSges(7,2) - t20 * mrSges(6,3) + t118 * mrSges(6,1) + t37 * mrSges(7,1) + t62 / 0.2e1 - t65 / 0.2e1 + Ifges(7,3) * t476 - Ifges(6,2) * t477 + t507 * t474 - t505 * t464) * t84 + (t239 * t506 + t508 * t83) * t474 + (t160 * t508 + t297 * t506) * t484 + (t314 / 0.2e1 + t514) * t334 + t237 * t480 + t160 * t527 + (Ifges(5,4) * t204 + Ifges(5,2) * t203 + Ifges(5,6) * t239) * t509 + (t345 * t340 + (t376 + t342) * t337) * t395 + (Ifges(6,4) * t83 + Ifges(6,6) * t239) * t477 + (Ifges(6,4) * t160 + Ifges(6,6) * t297) * t483; (Ifges(5,5) * t253 + Ifges(5,6) * t252) * t463 + t532 * t162 + (t233 - t129) * t161 + t522 * t113 + (t118 * t518 + t156 * t4 + t19 * t522 + t20 * t523 + t493 * t3 + t307 * t76) * m(6) + t523 * t112 + t524 * t114 + (t1 * t151 + t152 * t2 + t17 * t524 + t18 * t525 + t182 * t9 + t37 * t519) * m(7) + t525 * t111 + t501 * (t209 / 0.2e1 - t173 / 0.2e1) - t526 * (t210 / 0.2e1 - t172 / 0.2e1) + t314 + (((-m(4) * t170 - t226) * pkin(9) + t488) * t336 + (t362 * t468 + t360 * t462 - t343 + (-m(4) * t169 + m(5) * t155 + t403) * pkin(9) + t489) * t339) * qJD(3) + (-t502 / 0.2e1 + t528) * t285 + (Ifges(6,4) * t173 + Ifges(7,5) * t209 - Ifges(6,2) * t172 + Ifges(7,3) * t210) * t476 - t252 * t124 / 0.2e1 - t155 * (-mrSges(5,1) * t252 + mrSges(5,2) * t253) - t253 * t125 / 0.2e1 + t518 * t75 - t206 * t226 - t205 * t227 - t398 * t291 + (-t17 * t401 + t18 * t402) * mrSges(7,2) + (t19 * t401 + t20 * t402) * mrSges(6,3) + t514 - m(4) * (t169 * t205 + t170 * t206 + t249 * t291) - m(5) * (t128 * t89 + t129 * t90 + t155 * t187) + t519 * t74 + t493 * t48 + ((t376 + (Ifges(4,5) * t336 + Ifges(4,6) * t339) * t444 + (t430 / 0.2e1 + (t479 + t437 / 0.2e1) * t332) * qJD(1) - t342) * t337 + (-t316 / 0.2e1 + (-t433 / 0.2e1 + (t478 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t337) * t332) * qJD(1) + t343 * t339 - t488 * t336 - t345) * t340) * t396 + (-t172 * t505 + t173 * t506) * t465 + (t209 * t506 - t210 * t505) * t464 + (t209 * t508 + t210 * t507) * t474 + (t172 * t507 + t173 * t508) * t475 + (Ifges(6,4) * t209 + Ifges(7,5) * t173 - Ifges(6,2) * t210 + Ifges(7,3) * t172) * t477 + t350 * (Ifges(5,4) * t253 + Ifges(5,2) * t252) / 0.2e1 + (-mrSges(7,1) * t402 + mrSges(7,3) * t401) * t37 + (-mrSges(6,1) * t402 - mrSges(6,2) * t401) * t118 + (-t252 * t90 + t253 * t89) * mrSges(5,3) + t182 * t21 - t187 * t146 + t156 * t46 - pkin(2) * t150 + t151 * t49 + t152 * t47 + (Ifges(5,1) * t253 + Ifges(5,4) * t252) * t469 + t266 * t133 + t267 * t132 + (pkin(9) * t191 + t120 * mrSges(4,3) + t138 / 0.2e1 - t92 / 0.2e1 - t12 / 0.2e1 - t13 / 0.2e1 - t278 * mrSges(4,1) + t45 * mrSges(5,2) - t44 * mrSges(5,1) - t427 / 0.2e1 - t426 / 0.2e1 + t424 / 0.2e1 - t383 * t69 - t385 * t68 + (t382 - t384) * t225 + t344) * t339 - t288 * t287 + t307 * t22 + (-t467 * t505 + t487) * t284 + m(4) * (pkin(9) * t120 * t339 - pkin(2) * t278 - t121 * t330) + m(5) * (t108 * t330 + t232 * t89 + t233 * t90 + t266 * t44 + t267 * t45) + (t94 * t453 + t93 * t454 - t121 * mrSges(4,3) + t139 / 0.2e1 + t108 * t363 - t423 / 0.2e1 + t278 * mrSges(4,2) + t361 * t471 + t362 * t470 + t360 * t467 + t425 / 0.2e1 + (-t190 + t119) * pkin(9) + (-t331 * t45 - t333 * t44) * mrSges(5,3)) * t336; (-t226 + t441) * t169 + (-t17 * t399 + t18 * t400) * mrSges(7,2) + (t19 * t399 + t20 * t400) * mrSges(6,3) - (-Ifges(4,1) * t461 - t360 * t463 - t362 * t469 + t489 - t516) * t492 + t526 * (t174 / 0.2e1 - t296 / 0.2e1) + (-Ifges(6,4) * t175 - Ifges(7,5) * t295 - Ifges(6,2) * t174 + Ifges(7,3) * t296) * t476 + (-Ifges(6,4) * t295 - Ifges(7,5) * t175 - Ifges(6,2) * t296 + Ifges(7,3) * t174) * t477 + (-t174 * t505 - t175 * t506) * t465 + (t174 * t507 - t175 * t508) * t475 + t501 * (t175 / 0.2e1 - t295 / 0.2e1) + (t166 + t261) * t462 + (t132 * t333 - t133 * t331) * qJ(4) + (Ifges(5,5) * t331 + t333 * Ifges(5,6) - t305 * t505) * t467 + (t527 - t528) * t306 + t357 * mrSges(5,3) + (-pkin(3) * t108 + (-t331 * t89 + t333 * t90) * qJD(4) + t357 * qJ(4) - t115 * t89 - t116 * t90 - t155 * t170) * m(5) + (t48 + t49) * t245 + (t161 * t333 - t162 * t331) * qJD(4) - (t47 - t46) * t349 + (-t118 * t145 + t19 * t496 + t20 * t497 + t245 * t3 + t328 * t76 + t349 * t4) * m(6) + (t1 * t245 + t17 * t495 + t18 * t498 - t2 * t349 + t223 * t9 + t37 * t499) * m(7) + t487 * t305 + t381 + (-mrSges(7,1) * t400 + mrSges(7,3) * t399) * t37 + (-mrSges(6,1) * t400 - mrSges(6,2) * t399) * t118 - t116 * t161 - t115 * t162 - t145 * t75 - pkin(3) * t119 + t93 * t453 + t165 * t460 + (Ifges(5,1) * t331 + t435) * t470 + (Ifges(5,2) * t333 + t436) * t471 + (-t420 + t491) * t461 - t517 + t495 * t114 + t496 * t113 + t497 * t112 + t498 * t111 + t499 * t74 + (-t403 + t440) * t170 + t328 * t22 + (-t295 * t506 - t296 * t505) * t464 + (-t295 * t508 + t296 * t507) * t474 + t108 * (-mrSges(5,1) * t333 + mrSges(5,2) * t331) + (Ifges(5,5) * t469 - Ifges(4,2) * t462 + Ifges(6,6) * t476 + Ifges(7,6) * t477 + Ifges(5,3) * t463 + t504 * t465 + t506 * t475 + t515) * t269 + t331 * t480 + t223 * t21; t404 * t513 - t405 * t140 + t350 * t161 + t214 * t162 + t119 + t21 + t22 + (t140 * t18 - t17 * t513 + t9) * m(7) + (t140 * t20 + t19 * t513 + t76) * m(6) + (t214 * t89 + t350 * t90 + t108) * m(5); -t344 + (t404 + t438) * t20 + (t405 - t439) * t19 + (Ifges(7,3) * t513 - t431) * t477 - t37 * (mrSges(7,1) * t513 + mrSges(7,3) * t140) - t118 * (mrSges(6,1) * t513 - mrSges(6,2) * t140) + t65 * t474 + qJD(6) * t111 - t73 * t74 - pkin(5) * t47 + qJ(6) * t49 + (t140 * t17 + t18 * t513) * mrSges(7,2) + (-t140 * t506 - t505 * t513) * t465 + (-pkin(5) * t2 + qJ(6) * t1 - t17 * t20 + t18 * t494 - t37 * t73) * m(7) + (-Ifges(6,2) * t513 - t137 + t501) * t476 + (-t140 * t508 + t136 - t434 + t62) * t475 + t503; -t262 * t111 + t513 * t74 + 0.2e1 * (t2 / 0.2e1 + t37 * t474 + t18 * t465) * m(7) + t47;];
tauc  = t11(:);
