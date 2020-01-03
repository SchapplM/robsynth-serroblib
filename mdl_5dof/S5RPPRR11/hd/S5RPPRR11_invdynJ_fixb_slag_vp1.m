% Calculate vector of inverse dynamics joint torques for
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR11_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:32
% EndTime: 2019-12-31 18:05:53
% DurationCPUTime: 16.56s
% Computational Cost: add. (7722->808), mult. (19697->1077), div. (0->0), fcn. (18040->6), ass. (0->375)
t292 = cos(qJ(4));
t288 = sin(qJ(5));
t466 = rSges(6,2) * t288;
t291 = cos(qJ(5));
t468 = rSges(6,1) * t291;
t361 = -t466 + t468;
t511 = t292 * t361;
t290 = sin(qJ(1));
t293 = cos(qJ(1));
t289 = sin(qJ(4));
t432 = t293 * t291;
t437 = t290 * t288;
t179 = t289 * t437 - t432;
t436 = t290 * t291;
t180 = t288 * t293 + t289 * t436;
t172 = Icges(6,4) * t180;
t435 = t290 * t292;
t94 = -Icges(6,2) * t179 - Icges(6,6) * t435 + t172;
t171 = Icges(6,4) * t179;
t98 = -Icges(6,1) * t180 + Icges(6,5) * t435 + t171;
t359 = t179 * t94 + t180 * t98;
t438 = t289 * t293;
t181 = -t288 * t438 - t436;
t182 = t289 * t432 - t437;
t433 = t292 * t293;
t451 = Icges(6,4) * t182;
t96 = Icges(6,2) * t181 - Icges(6,6) * t433 + t451;
t173 = Icges(6,4) * t181;
t99 = Icges(6,1) * t182 - Icges(6,5) * t433 + t173;
t469 = t181 * t96 + t182 * t99;
t92 = -Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t435;
t93 = Icges(6,5) * t182 + Icges(6,6) * t181 - Icges(6,3) * t433;
t519 = t359 + t469 + (-t290 * t92 - t293 * t93) * t292;
t470 = t181 * t94 - t182 * t98;
t471 = -t179 * t96 + t180 * t99;
t518 = t470 + t292 * (-t290 * t93 + t293 * t92) + t471;
t100 = t180 * rSges(6,1) - t179 * rSges(6,2) - rSges(6,3) * t435;
t166 = rSges(6,3) * t289 + t511;
t400 = qJD(5) * t292;
t402 = qJD(4) * t293;
t203 = -t290 * t400 + t402;
t476 = t292 * pkin(4);
t234 = pkin(7) * t289 + t476;
t401 = qJD(5) * t289;
t251 = qJD(1) + t401;
t516 = -t100 * t251 + t166 * t203 + t234 * t402;
t353 = t288 * t94 + t291 * t98;
t37 = -t289 * t92 - t292 * t353;
t338 = Icges(6,5) * t291 - Icges(6,6) * t288;
t153 = Icges(6,3) * t289 + t292 * t338;
t449 = Icges(6,4) * t291;
t340 = -Icges(6,2) * t288 + t449;
t157 = Icges(6,6) * t289 + t292 * t340;
t450 = Icges(6,4) * t288;
t342 = Icges(6,1) * t291 - t450;
t161 = Icges(6,5) * t289 + t292 * t342;
t307 = t153 * t435 + t157 * t179 - t161 * t180;
t513 = t307 * t251;
t225 = t290 * rSges(4,2) + t293 * rSges(4,3);
t229 = t293 * pkin(1) + t290 * qJ(2);
t279 = t293 * qJ(3);
t508 = t279 + t229;
t512 = t225 + t508;
t383 = t289 * t402;
t407 = qJD(1) * t290;
t385 = t292 * t407;
t311 = t383 + t385;
t231 = -rSges(3,2) * t293 + t290 * rSges(3,3);
t286 = t292 * pkin(7);
t479 = pkin(4) * t289;
t507 = t286 - t479;
t339 = Icges(5,5) * t289 + Icges(5,6) * t292;
t154 = Icges(5,3) * t293 + t290 * t339;
t506 = qJD(1) * t154;
t453 = Icges(5,4) * t289;
t341 = Icges(5,2) * t292 + t453;
t158 = Icges(5,6) * t293 + t290 * t341;
t452 = Icges(5,4) * t292;
t343 = Icges(5,1) * t289 + t452;
t162 = Icges(5,5) * t293 + t290 * t343;
t362 = rSges(5,1) * t289 + rSges(5,2) * t292;
t165 = rSges(5,3) * t293 + t290 * t362;
t369 = qJD(1) * t289 + qJD(5);
t381 = t292 * t402;
t503 = t290 * t369 - t381;
t216 = Icges(5,5) * t292 - Icges(5,6) * t289;
t315 = qJD(4) * t216;
t159 = -Icges(5,6) * t290 + t293 * t341;
t244 = Icges(5,4) * t433;
t448 = Icges(5,5) * t290;
t163 = Icges(5,1) * t438 + t244 - t448;
t333 = t159 * t292 + t163 * t289;
t502 = qJD(1) * t333 + t293 * t315 - t506;
t334 = t158 * t292 + t162 * t289;
t155 = -Icges(5,3) * t290 + t293 * t339;
t409 = qJD(1) * t155;
t501 = qJD(1) * t334 + t290 * t315 + t409;
t218 = -Icges(5,2) * t289 + t452;
t220 = Icges(5,1) * t292 - t453;
t331 = t218 * t292 + t220 * t289;
t500 = t331 * qJD(1) - t339 * qJD(4);
t499 = t290 * (-Icges(5,2) * t438 + t163 + t244) - t293 * (t218 * t290 + t162);
t152 = Icges(6,3) * t292 - t289 * t338;
t404 = qJD(4) * t290;
t202 = -t293 * t400 - t404;
t335 = t157 * t288 - t161 * t291;
t352 = t288 * t96 - t291 * t99;
t498 = -t202 * (t153 * t293 + t352) - t203 * (t153 * t290 + t353) - t251 * (t152 + t335);
t186 = (-Icges(6,2) * t291 - t450) * t292;
t497 = t202 * (-Icges(6,2) * t182 + t173 + t99) + t203 * (-Icges(6,2) * t180 - t171 - t98) + t251 * (t161 + t186);
t306 = -qJDD(5) * t292 + (-qJD(1) + t401) * qJD(4);
t380 = qJD(1) * t400;
t117 = (-qJDD(4) + t380) * t290 + t306 * t293;
t496 = t117 / 0.2e1;
t271 = qJDD(4) * t293;
t118 = t290 * t306 - t293 * t380 + t271;
t495 = t118 / 0.2e1;
t494 = -t202 / 0.2e1;
t493 = t202 / 0.2e1;
t492 = -t203 / 0.2e1;
t491 = t203 / 0.2e1;
t396 = qJD(1) * qJD(4);
t211 = -qJDD(4) * t290 - t293 * t396;
t490 = t211 / 0.2e1;
t212 = -t290 * t396 + t271;
t489 = t212 / 0.2e1;
t488 = -t251 / 0.2e1;
t487 = t251 / 0.2e1;
t485 = t290 / 0.2e1;
t484 = -t293 / 0.2e1;
t482 = rSges(3,2) - pkin(1);
t481 = -rSges(5,3) - pkin(6);
t480 = -rSges(6,3) - pkin(7);
t478 = pkin(6) * t293;
t477 = g(1) * t290;
t406 = qJD(1) * t293;
t384 = t292 * t406;
t312 = t289 * t404 - t384;
t403 = qJD(4) * t292;
t382 = t290 * t403;
t84 = -t251 * t436 + (-t293 * t369 - t382) * t288;
t85 = t369 * t432 + (-t251 * t288 + t291 * t403) * t290;
t44 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t312;
t46 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t312;
t48 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t312;
t7 = (qJD(4) * t353 + t44) * t289 + (-qJD(4) * t92 - t288 * t46 + t291 * t48 + (t288 * t98 - t291 * t94) * qJD(5)) * t292;
t475 = t7 * t203;
t327 = t251 * t293;
t82 = t288 * t503 - t291 * t327;
t83 = -t288 * t327 - t291 * t503;
t43 = Icges(6,5) * t83 + Icges(6,6) * t82 + Icges(6,3) * t311;
t45 = Icges(6,4) * t83 + Icges(6,2) * t82 + Icges(6,6) * t311;
t47 = Icges(6,1) * t83 + Icges(6,4) * t82 + Icges(6,5) * t311;
t8 = (qJD(4) * t352 + t43) * t289 + (qJD(4) * t93 - t288 * t45 + t291 * t47 + (-t288 * t99 - t291 * t96) * qJD(5)) * t292;
t474 = t8 * t202;
t473 = -pkin(1) - qJ(3);
t183 = (-Icges(6,5) * t288 - Icges(6,6) * t291) * t292;
t103 = qJD(4) * t152 + qJD(5) * t183;
t156 = Icges(6,6) * t292 - t289 * t340;
t106 = qJD(4) * t156 + qJD(5) * t186;
t160 = Icges(6,5) * t292 - t289 * t342;
t189 = (-Icges(6,1) * t288 - t449) * t292;
t109 = qJD(4) * t160 + qJD(5) * t189;
t18 = (qJD(4) * t335 + t103) * t289 + (qJD(4) * t153 - t106 * t288 + t109 * t291 + (-t157 * t291 - t161 * t288) * qJD(5)) * t292;
t201 = qJD(4) * t400 + qJDD(5) * t289 + qJDD(1);
t59 = t153 * t289 - t292 * t335;
t472 = t18 * t251 + t59 * t201;
t465 = rSges(3,3) * t293;
t464 = rSges(5,3) * t290;
t193 = (-rSges(6,1) * t288 - rSges(6,2) * t291) * t292;
t284 = t292 * rSges(6,3);
t112 = qJD(5) * t193 + (-t289 * t361 + t284) * qJD(4);
t130 = t312 * pkin(7) + (t289 * t406 + t382) * pkin(4);
t277 = qJD(2) * t293;
t174 = qJD(1) * t229 - t277;
t255 = pkin(7) * t435;
t439 = t289 * t290;
t393 = pkin(4) * t439;
t195 = -t255 + t393;
t280 = t293 * qJ(2);
t224 = pkin(1) * t290 - t280;
t443 = qJ(3) * t290;
t360 = -t443 - t478;
t325 = -t224 + t360;
t313 = -t195 + t325;
t210 = t507 * qJD(4);
t294 = qJD(1) ^ 2;
t442 = qJ(3) * t294;
t328 = qJD(4) * t210 - t442;
t397 = qJD(1) * qJD(3);
t398 = qJD(1) * qJD(2);
t415 = qJDD(2) * t290 + t293 * t398;
t367 = qJDD(3) * t293 - 0.2e1 * t290 * t397 + t415;
t434 = t290 * t294;
t344 = pkin(6) * t434 + t367;
t365 = rSges(6,1) * t85 + rSges(6,2) * t84;
t50 = rSges(6,3) * t312 + t365;
t12 = -t100 * t201 + t112 * t203 + t118 * t166 + t212 * t234 - t251 * t50 + t328 * t293 + (-t130 - t174) * qJD(1) + t313 * qJDD(1) + t344;
t462 = t12 * t290;
t421 = t182 * rSges(6,1) + t181 * rSges(6,2);
t102 = -rSges(6,3) * t433 + t421;
t389 = pkin(4) * t381 + pkin(7) * t311;
t129 = -qJD(1) * t393 + t389;
t258 = pkin(4) * t438;
t197 = -pkin(7) * t433 + t258;
t263 = qJ(2) * t406;
t276 = qJD(2) * t290;
t414 = t263 + t276;
t390 = qJD(1) * (-pkin(1) * t407 + t414) + qJDD(1) * t229 + t290 * t398;
t322 = qJDD(1) * t279 + qJDD(3) * t290 + 0.2e1 * t293 * t397 + t390;
t303 = (-pkin(6) * t294 - qJDD(2)) * t293 + t322;
t444 = pkin(6) * qJDD(1);
t49 = t83 * rSges(6,1) + t82 * rSges(6,2) + rSges(6,3) * t311;
t13 = qJD(1) * t129 + qJDD(1) * t197 + t102 * t201 - t112 * t202 - t117 * t166 - t211 * t234 + t251 * t49 + (t328 - t444) * t290 + t303;
t461 = t13 * t293;
t228 = rSges(5,1) * t292 - rSges(5,2) * t289;
t199 = t228 * t402;
t314 = -t165 + t325;
t275 = qJD(3) * t293;
t411 = t275 + t276;
t66 = qJD(1) * t314 + t199 + t411;
t460 = t290 * t66;
t457 = t37 * t118;
t38 = t289 * t93 - t292 * t352;
t456 = t38 * t117;
t441 = t154 * t290;
t440 = t154 * t293;
t184 = t216 * t290;
t185 = t216 * t293;
t149 = t293 * t155;
t86 = t290 * t331 + t185;
t431 = t86 * qJD(1);
t430 = -t100 - t195;
t429 = t102 + t197;
t428 = t112 + t210;
t427 = t158 * t433 + t162 * t438;
t426 = t159 * t433 + t163 * t438;
t422 = t166 + t234;
t420 = -t341 + t220;
t419 = -t218 - t343;
t226 = rSges(3,2) * t290 + t465;
t418 = -t224 + t226;
t170 = t229 + t231;
t417 = t289 * t466 + t284;
t416 = rSges(5,1) * t438 + rSges(5,2) * t433;
t196 = pkin(4) * t435 + pkin(7) * t439;
t198 = pkin(4) * t433 + pkin(7) * t438;
t413 = rSges(3,2) * t407 + rSges(3,3) * t406;
t412 = pkin(6) * t407 + t277;
t410 = -qJD(1) * t224 + t276;
t408 = qJD(1) * t339;
t405 = qJD(4) * t289;
t399 = -m(4) - m(5) - m(6);
t395 = qJDD(2) * t293;
t394 = -rSges(4,3) + t473;
t58 = t159 * t435 + t163 * t439 + t149;
t388 = t263 + t411;
t387 = t275 + t410;
t379 = -t407 / 0.2e1;
t378 = -t406 / 0.2e1;
t377 = -t404 / 0.2e1;
t376 = t404 / 0.2e1;
t375 = -t402 / 0.2e1;
t374 = t402 / 0.2e1;
t373 = -t401 / 0.2e1;
t230 = rSges(4,2) * t293 - t290 * rSges(4,3);
t371 = t230 - t443;
t370 = -qJD(3) * t290 + t277;
t368 = qJD(4) * t228 + qJD(3);
t366 = t473 - t479;
t363 = -t224 + t371;
t232 = rSges(2,1) * t293 - rSges(2,2) * t290;
t227 = rSges(2,1) * t290 + rSges(2,2) * t293;
t316 = qJD(4) * t218;
t107 = -qJD(1) * t158 + t293 * t316;
t317 = qJD(4) * t220;
t110 = -qJD(1) * t162 + t293 * t317;
t81 = -t159 * t289 + t163 * t292;
t298 = qJD(4) * t81 + t107 * t292 + t110 * t289 - t409;
t108 = qJD(1) * t159 + t290 * t316;
t111 = t290 * t317 + (t293 * t343 - t448) * qJD(1);
t80 = -t158 * t289 + t162 * t292;
t299 = qJD(4) * t80 + t108 * t292 + t111 * t289 - t506;
t358 = (-t290 * t501 + t299 * t293) * t293 - (-t290 * t502 + t298 * t293) * t290;
t357 = (t299 * t290 + t293 * t501) * t293 - (t298 * t290 + t293 * t502) * t290;
t192 = t228 * t290;
t114 = qJD(4) * t192 + (t293 * t362 - t464) * qJD(1);
t207 = t362 * qJD(4);
t329 = -qJD(4) * t207 - t442;
t25 = t212 * t228 + t329 * t293 + (-t114 - t174) * qJD(1) + t314 * qJDD(1) + t344;
t326 = rSges(5,1) * t381 - rSges(5,2) * t383;
t113 = -t165 * qJD(1) + t326;
t167 = t416 - t464;
t26 = qJD(1) * t113 + qJDD(1) * t167 - t211 * t228 + (t329 - t444) * t290 + t303;
t356 = t25 * t293 + t26 * t290;
t27 = t435 * t92 - t359;
t28 = -t435 * t93 + t471;
t355 = t27 * t293 - t28 * t290;
t354 = t27 * t290 + t28 * t293;
t29 = t433 * t92 + t470;
t30 = -t433 * t93 + t469;
t351 = t29 * t293 - t290 * t30;
t350 = t29 * t290 + t293 * t30;
t349 = t290 * t38 - t293 * t37;
t348 = t290 * t37 + t293 * t38;
t57 = t290 * t334 + t440;
t347 = -t290 * t58 + t293 * t57;
t60 = t427 - t441;
t61 = -t155 * t290 + t426;
t346 = -t290 * t61 + t293 * t60;
t67 = t368 * t290 + (t167 + t508) * qJD(1) - t412;
t345 = t290 * t67 + t293 * t66;
t337 = t100 * t293 - t102 * t290;
t336 = -t113 * t293 - t114 * t290;
t332 = -t165 * t290 - t167 * t293;
t330 = -t218 * t289 + t220 * t292;
t141 = rSges(6,3) * t439 + t290 * t511;
t142 = rSges(6,3) * t438 + t293 * t511;
t318 = -t362 + t473;
t310 = t153 * t251 + t202 * t93 - t203 * t92;
t309 = (-Icges(6,5) * t179 - Icges(6,6) * t180) * t203 + (Icges(6,5) * t181 - Icges(6,6) * t182) * t202 + t183 * t251;
t308 = -t158 * t293 + t159 * t290;
t305 = t292 * t309;
t304 = (t289 * t419 + t292 * t420) * qJD(1);
t301 = (Icges(6,1) * t181 - t451 - t96) * t202 + (-Icges(6,1) * t179 - t172 - t94) * t203 + (-t157 + t189) * t251;
t34 = t100 * t202 - t102 * t203 + (-t195 * t290 - t197 * t293) * qJD(4);
t35 = qJD(1) * t313 + t411 + t516;
t36 = t102 * t251 - t166 * t202 + (qJD(4) * t234 + qJD(3)) * t290 + (t197 + t508) * qJD(1) - t412;
t300 = t34 * t337 + (t290 * t35 - t293 * t36) * t166;
t205 = t341 * qJD(4);
t206 = t343 * qJD(4);
t297 = -qJD(1) * t216 + qJD(4) * t330 - t205 * t292 - t206 * t289;
t296 = t308 * t289 - t292 * t499;
t295 = t498 * t292;
t268 = rSges(4,2) * t406;
t194 = t228 * t293;
t164 = -t289 * t468 + t417;
t140 = t161 * t293;
t139 = t161 * t290;
t138 = t157 * t293;
t137 = t157 * t290;
t134 = qJD(1) * t170 - t277;
t133 = qJD(1) * t418 + t276;
t128 = rSges(6,1) * t181 - rSges(6,2) * t182;
t127 = -rSges(6,1) * t179 - rSges(6,2) * t180;
t116 = qJD(1) * t512 - t370;
t115 = qJD(1) * t363 + t411;
t87 = t293 * t331 - t184;
t77 = t87 * qJD(1);
t76 = t332 * qJD(4);
t69 = qJD(1) * t413 + qJDD(1) * t231 + t390 - t395;
t68 = t418 * qJDD(1) + (-qJD(1) * t231 - t174) * qJD(1) + t415;
t55 = -t395 - qJ(3) * t434 + qJDD(1) * t225 + qJD(1) * (-rSges(4,3) * t407 + t268) + t322;
t54 = -t294 * t279 + t363 * qJDD(1) + (-qJD(1) * t225 - t174) * qJD(1) + t367;
t53 = -t153 * t433 + t157 * t181 + t161 * t182;
t51 = t53 * t251;
t42 = t297 * t290 + t293 * t500;
t41 = -t290 * t500 + t297 * t293;
t40 = -qJD(4) * t333 - t107 * t289 + t110 * t292;
t39 = -qJD(4) * t334 - t108 * t289 + t111 * t292;
t24 = qJD(4) * t346 + t77;
t23 = qJD(4) * t347 + t431;
t16 = -t103 * t435 - t106 * t179 + t109 * t180 + t153 * t312 + t157 * t84 + t161 * t85;
t15 = -t103 * t433 + t106 * t181 + t109 * t182 + t153 * t311 + t157 * t82 + t161 * t83;
t14 = t202 * t38 + t203 * t37 + t251 * t59;
t11 = t202 * t30 + t203 * t29 + t51;
t10 = t202 * t28 + t203 * t27 - t513;
t9 = t100 * t117 - t102 * t118 + t195 * t211 - t197 * t212 + t202 * t50 - t203 * t49 + (-t129 * t293 - t130 * t290) * qJD(4);
t6 = -t93 * t384 - t179 * t45 + t180 * t47 + t84 * t96 + t85 * t99 + (-t292 * t43 + t405 * t93) * t290;
t5 = t92 * t384 - t179 * t46 + t180 * t48 + t84 * t94 - t85 * t98 + (-t292 * t44 - t405 * t92) * t290;
t4 = t93 * t383 + t181 * t45 + t182 * t47 + t82 * t96 + t83 * t99 + (-t293 * t43 + t407 * t93) * t292;
t3 = -t92 * t383 + t181 * t46 + t182 * t48 + t82 * t94 - t83 * t98 + (-t293 * t44 - t407 * t92) * t292;
t2 = t117 * t28 + t118 * t27 + t16 * t251 - t201 * t307 + t202 * t6 + t203 * t5;
t1 = t117 * t30 + t118 * t29 + t15 * t251 + t201 * t53 + t202 * t4 + t203 * t3;
t17 = [(t133 * t277 + t134 * (t413 + t414) + (t133 * t482 * t293 + (t133 * (-rSges(3,3) - qJ(2)) - t134 * pkin(1)) * t290) * qJD(1) - (qJD(1) * t226 - t133 + t410) * t134 + (t69 - g(2)) * t170 + (t68 - g(1)) * (t290 * t482 + t280 + t465)) * m(3) - t307 * t495 + (t66 * t412 + t67 * (t326 + t388) - t368 * t460 + ((t318 * t66 + t481 * t67) * t293 + (t66 * (rSges(5,3) - qJ(2)) + t67 * t318) * t290) * qJD(1) - (t199 - t66 + (-t165 + t360) * qJD(1) + t387) * t67 + (t26 - g(2)) * (t290 * t481 + t416 + t508) + (t25 - g(1)) * (t318 * t290 + t481 * t293 + t280)) * m(5) + (t16 + t11) * t491 + (-qJD(4) * t331 + t205 * t289 - t206 * t292) * qJD(1) + (-(qJD(1) * t371 - t115 + t387) * t116 + t115 * t370 + t116 * (t268 + t388) + (t115 * t394 * t293 + (t115 * (-rSges(4,2) - qJ(2)) + t116 * t394) * t290) * qJD(1) + (t55 - g(2)) * t512 + (t54 - g(1)) * (t290 * t473 + t230 + t280)) * m(4) + (t40 + t41) * t377 + (t81 + t87) * t490 + t15 * t493 + t53 * t496 + (t39 + t42 + t24) * t374 + t456 / 0.2e1 + t457 / 0.2e1 + (m(2) * (t227 ^ 2 + t232 ^ 2) + t330 + Icges(2,3) + Icges(3,1) + Icges(4,1)) * qJDD(1) + (t77 + (t427 * t293 + (-t57 + (t155 + t334) * t290 - t426) * t290) * qJD(4)) * t375 + (t80 + t86) * t489 + t474 / 0.2e1 + t475 / 0.2e1 + (-t366 * t477 + t35 * (-t365 + t412) + t36 * (t49 + t388 + t389) + (t12 * t366 + (-qJD(3) + (t289 * t480 - t476) * qJD(4)) * t35) * t290 + ((-t35 * qJ(2) + t36 * t366) * t290 + (t35 * (t507 + t284 + t473) - t36 * pkin(6)) * t293) * qJD(1) - (-t35 + t387 + (-t195 + t360) * qJD(1) + t516) * t36 + (t13 - g(2)) * (-pkin(6) * t290 + t433 * t480 + t258 + t421 + t508) + (t12 - g(1)) * (-t100 + t255 + t280 - t478)) * m(6) + (t513 + (-t30 + t519) * t203 + (t29 - t518) * t202 + t10) * t494 + (t51 + (-t28 + t518) * t203 + (t27 + t519) * t202) * t492 + t472 + (-t431 + ((t426 - t61 - t440) * t293 + (-t149 + t58 - t60 - t441) * t290) * qJD(4) + t23) * t376 - m(2) * (-g(1) * t227 + g(2) * t232); (-m(3) + t399) * (-g(2) * t293 + t477) + 0.2e1 * (t462 / 0.2e1 - t461 / 0.2e1) * m(6) + 0.2e1 * (t25 * t485 + t26 * t484) * m(5) + 0.2e1 * (t484 * t55 + t485 * t54) * m(4) + 0.2e1 * (t484 * t69 + t485 * t68) * m(3); t399 * (g(1) * t293 + g(2) * t290) + m(4) * (t290 * t55 + t293 * t54) + m(5) * t356 + m(6) * (t12 * t293 + t13 * t290); ((-t138 * t179 + t140 * t180) * t202 + (-t137 * t179 + t139 * t180) * t203 + (-t156 * t179 + t160 * t180) * t251 + (t28 * t438 - t292 * t307) * qJD(5) + ((qJD(5) * t27 + t310) * t289 + t295) * t290) * t492 + (((-t138 * t288 + t140 * t291 + t93) * t202 + (-t137 * t288 + t139 * t291 - t92) * t203 + (-t156 * t288 + t160 * t291 + t153) * t251 + t59 * qJD(5)) * t292 + (qJD(5) * t348 - t498) * t289) * t488 + (-(-t192 * t66 + t194 * t67) * qJD(1) - (t76 * (-t192 * t290 - t194 * t293) - t345 * t362) * qJD(4) + (qJD(4) * t336 + t165 * t211 - t167 * t212) * t332 + t76 * ((-t165 * t293 + t167 * t290) * qJD(1) + t336) - t345 * t207 + ((t293 * t67 - t460) * qJD(1) + t356) * t228 - g(1) * t194 - g(2) * t192 + g(3) * t362) * m(5) + (-t35 * (-qJD(1) * t196 - t141 * t251 + t164 * t203 + t402 * t507) - t36 * (qJD(1) * t198 + t142 * t251 - t164 * t202 + t404 * t507) - t34 * (t141 * t202 - t142 * t203 - t196 * t404 - t198 * t402) - ((-t100 * t35 + t102 * t36) * t292 + t300 * t289) * qJD(5) - g(1) * (t142 + t198) - g(2) * (t141 + t196) - g(3) * (t286 + (-pkin(4) - t468) * t289 + t417) + (t12 * t422 + t35 * t428 - t9 * t429 + t34 * (-t129 - t49) + (t34 * t430 + t36 * t422) * qJD(1)) * t293 + (t13 * t422 + t36 * t428 + t9 * t430 + t34 * (-t130 - t50) + (t34 * t429 - t35 * t422) * qJD(1)) * t290) * m(6) + (qJD(1) * t42 + qJD(4) * t357 + qJDD(1) * t86 + t211 * t58 + t212 * t57 + t2) * t293 / 0.2e1 - (qJD(1) * t41 + qJD(4) * t358 + qJDD(1) * t87 + t211 * t61 + t212 * t60 + t1) * t290 / 0.2e1 + ((t185 * t404 + t408) * t290 + (t304 + (-t290 * t184 + t296) * qJD(4)) * t293) * t376 + ((t184 * t402 - t408) * t293 + (t304 + (-t293 * t185 + t296) * qJD(4)) * t290) * t375 - t201 * t349 / 0.2e1 + t347 * t489 + t346 * t490 + (-qJD(1) * t354 - t290 * t6 + t293 * t5) * t491 + (-qJD(1) * t350 - t290 * t4 + t293 * t3) * t493 + ((t138 * t181 + t140 * t182) * t202 + (t137 * t181 + t139 * t182) * t203 + (t156 * t181 + t160 * t182) * t251 + (t29 * t439 + t292 * t53) * qJD(5) + ((qJD(5) * t30 + t310) * t289 + t295) * t293) * t494 + t355 * t495 + t351 * t496 + (-qJD(1) * t348 - t290 * t8 + t293 * t7) * t487 + (t290 * t373 + t379) * t10 + qJD(1) * (-t290 * t40 + t293 * t39 + (-t80 * t290 - t293 * t81) * qJD(1)) / 0.2e1 + t24 * t378 + t23 * t379 + ((-t60 * t290 - t61 * t293) * qJD(1) + t358) * t377 - qJD(1) * ((-t420 * t289 + t419 * t292) * qJD(1) + (t289 * t499 + t308 * t292) * qJD(4)) / 0.2e1 + (t293 * t373 + t378) * t11 - t14 * t400 / 0.2e1 + ((-t57 * t290 - t58 * t293) * qJD(1) + t357) * t374 + qJDD(1) * (-t290 * t81 + t293 * t80) / 0.2e1; -t1 * t433 / 0.2e1 + (t289 * t53 - t292 * t350) * t496 + ((qJD(4) * t350 + t15) * t289 + (-qJD(1) * t351 + qJD(4) * t53 - t290 * t3 - t293 * t4) * t292) * t493 - t2 * t435 / 0.2e1 + (-t289 * t307 - t292 * t354) * t495 + ((qJD(4) * t354 + t16) * t289 + (-qJD(1) * t355 - qJD(4) * t307 - t290 * t5 - t293 * t6) * t292) * t491 + t14 * t403 / 0.2e1 + t289 * (t456 + t457 + t472 + t474 + t475) / 0.2e1 + t201 * (t289 * t59 - t292 * t348) / 0.2e1 + ((qJD(4) * t348 + t18) * t289 + (qJD(1) * t349 + qJD(4) * t59 - t290 * t7 - t293 * t8) * t292) * t487 + (t181 * t497 + t301 * t182 - t293 * t305) * t494 + (-t179 * t497 + t180 * t301 - t290 * t305) * t492 + (t309 * t289 + (-t288 * t497 + t291 * t301) * t292) * t488 + (t385 / 0.2e1 + t289 * t374) * t11 + (t289 * t376 + t292 * t378) * t10 + ((qJD(4) * t300 - t12 * t100 + t13 * t102 - t35 * t50 + t36 * t49) * t289 + (t35 * (-qJD(4) * t100 - t112 * t290) + t36 * (qJD(4) * t102 + t112 * t293) - t9 * t337 + t34 * (t100 * t407 + t102 * t406 + t290 * t49 - t293 * t50) + (-t462 + t461 + (-t290 * t36 - t293 * t35) * qJD(1)) * t166) * t292 - t35 * (-t127 * t251 + t193 * t203) - t36 * (t128 * t251 - t193 * t202) - t34 * (t127 * t202 - t128 * t203) - g(1) * t128 - g(2) * t127 - g(3) * t193) * m(6);];
tau = t17;
