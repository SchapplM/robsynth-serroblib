% Calculate vector of inverse dynamics joint torques for
% S5RPPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:58:06
% DurationCPUTime: 21.48s
% Computational Cost: add. (20724->785), mult. (19257->1042), div. (0->0), fcn. (17824->10), ass. (0->377)
t295 = qJ(1) + pkin(8);
t292 = cos(t295);
t298 = -pkin(6) - qJ(3);
t263 = t292 * t298;
t297 = cos(pkin(9));
t287 = pkin(3) * t297 + pkin(2);
t290 = sin(t295);
t423 = -t287 * t290 - t263;
t300 = sin(qJ(1));
t492 = pkin(1) * t300;
t530 = t423 - t492;
t302 = cos(qJ(1));
t293 = t302 * pkin(1);
t301 = cos(qJ(5));
t444 = t292 * t301;
t294 = pkin(9) + qJ(4);
t291 = cos(t294);
t299 = sin(qJ(5));
t448 = t291 * t299;
t192 = t290 * t448 + t444;
t445 = t292 * t299;
t447 = t291 * t301;
t193 = t290 * t447 - t445;
t367 = rSges(6,1) * t193 - rSges(6,2) * t192;
t289 = sin(t294);
t454 = t289 * t290;
t108 = -rSges(6,3) * t454 - t367;
t406 = qJD(1) * qJD(4);
t211 = -qJDD(4) * t292 + t290 * t406;
t412 = qJD(4) * t290;
t395 = t291 * t412;
t414 = qJD(1) * t292;
t325 = t289 * t414 + t395;
t405 = qJDD(5) * t289;
t113 = qJD(5) * t325 + t290 * t405 + t211;
t120 = t325 * pkin(7) + (-t289 * t412 + t291 * t414) * pkin(4);
t201 = (-rSges(6,1) * t299 - rSges(6,2) * t301) * t289;
t278 = t289 * rSges(6,3);
t481 = rSges(6,2) * t299;
t484 = rSges(6,1) * t301;
t366 = -t481 + t484;
t121 = qJD(5) * t201 + (t291 * t366 + t278) * qJD(4);
t167 = -rSges(6,3) * t291 + t289 * t366;
t410 = qJD(5) * t289;
t411 = qJD(4) * t292;
t200 = -t290 * t410 + t411;
t204 = qJD(4) * t410 - qJDD(5) * t291 + qJDD(1);
t285 = t291 * pkin(4);
t515 = pkin(7) * t289 + t285;
t209 = t515 * qJD(4);
t490 = t289 * pkin(4);
t228 = -pkin(7) * t291 + t490;
t409 = qJD(5) * t291;
t258 = qJD(1) - t409;
t187 = t515 * t290;
t273 = t292 * qJ(3);
t223 = pkin(2) * t290 - t273;
t149 = t223 + t423;
t384 = -t223 - t492;
t377 = t149 + t384;
t337 = -t187 + t377;
t407 = qJD(1) * qJD(3);
t303 = qJD(1) ^ 2;
t524 = t303 * t293;
t339 = qJDD(3) * t290 + t292 * t407 - t524;
t272 = t290 * qJ(3);
t226 = pkin(2) * t292 + t272;
t271 = qJD(3) * t292;
t179 = qJD(1) * t226 - t271;
t416 = qJD(1) * t290;
t252 = t298 * t416;
t487 = pkin(2) - t287;
t434 = t252 - (-t292 * t487 - t272) * qJD(1) - t179;
t413 = qJD(4) * t289;
t319 = t258 * t301 + t299 * t413;
t415 = qJD(1) * t291;
t379 = -qJD(5) + t415;
t88 = t290 * t319 - t379 * t445;
t318 = t258 * t299 - t301 * t413;
t89 = t290 * t318 + t379 * t444;
t373 = rSges(6,1) * t89 + rSges(6,2) * t88;
t50 = rSges(6,3) * t325 + t373;
t10 = -t209 * t411 + t108 * t204 + t113 * t167 - t121 * t200 + t211 * t228 - t258 * t50 + (-t120 + t434) * qJD(1) + t337 * qJDD(1) + t339;
t480 = t10 * t290;
t213 = qJD(1) * t223;
t270 = qJD(3) * t290;
t529 = -qJD(1) * t149 + t213 + t270;
t452 = t290 * t291;
t402 = rSges(5,1) * t452;
t528 = -t402 + t530;
t527 = -t258 * t108 + t167 * t200 + t228 * t411 - t270;
t224 = rSges(3,1) * t290 + rSges(3,2) * t292;
t202 = -t224 - t492;
t177 = Icges(6,4) * t193;
t101 = -Icges(6,2) * t192 + Icges(6,6) * t454 + t177;
t176 = Icges(6,4) * t192;
t105 = -Icges(6,1) * t193 - Icges(6,5) * t454 + t176;
t522 = t101 * t299 + t105 * t301;
t98 = Icges(6,5) * t193 - Icges(6,6) * t192 + Icges(6,3) * t454;
t41 = -t289 * t522 - t291 * t98;
t199 = t292 * t410 + t412;
t27 = -t101 * t192 - t105 * t193 + t454 * t98;
t450 = t290 * t301;
t194 = -t291 * t445 + t450;
t451 = t290 * t299;
t195 = t291 * t444 + t451;
t453 = t289 * t292;
t100 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t453;
t468 = Icges(6,4) * t195;
t103 = Icges(6,2) * t194 + Icges(6,6) * t453 + t468;
t178 = Icges(6,4) * t194;
t106 = Icges(6,1) * t195 + Icges(6,5) * t453 + t178;
t28 = t100 * t454 - t103 * t192 + t106 * t193;
t349 = Icges(6,5) * t301 - Icges(6,6) * t299;
t161 = -Icges(6,3) * t291 + t289 * t349;
t466 = Icges(6,4) * t301;
t350 = -Icges(6,2) * t299 + t466;
t163 = -Icges(6,6) * t291 + t289 * t350;
t467 = Icges(6,4) * t299;
t352 = Icges(6,1) * t301 - t467;
t165 = -Icges(6,5) * t291 + t289 * t352;
t52 = t161 * t454 - t163 * t192 + t165 * t193;
t12 = t199 * t28 - t200 * t27 + t258 * t52;
t29 = t101 * t194 - t105 * t195 + t453 * t98;
t30 = t100 * t453 + t103 * t194 + t106 * t195;
t53 = t161 * t453 + t163 * t194 + t165 * t195;
t13 = t199 * t30 - t200 * t29 + t258 * t53;
t296 = sin(pkin(9));
t482 = rSges(4,2) * t296;
t485 = rSges(4,1) * t297;
t160 = t290 * rSges(4,3) + (-t482 + t485) * t292;
t383 = t226 + t293;
t131 = t160 + t383;
t227 = rSges(3,1) * t292 - rSges(3,2) * t290;
t203 = t227 + t293;
t516 = -rSges(5,2) * t454 - rSges(5,3) * t292;
t277 = Icges(5,4) * t291;
t351 = -Icges(5,2) * t289 + t277;
t218 = Icges(5,1) * t289 + t277;
t215 = Icges(5,5) * t291 - Icges(5,6) * t289;
t214 = Icges(5,5) * t289 + Icges(5,6) * t291;
t327 = qJD(4) * t214;
t469 = Icges(5,4) * t289;
t219 = Icges(5,1) * t291 - t469;
t156 = Icges(5,5) * t290 + t219 * t292;
t154 = Icges(5,6) * t290 + t292 * t351;
t457 = t154 * t289;
t344 = -t156 * t291 + t457;
t460 = Icges(5,3) * t292;
t512 = -t292 * t327 + (-t215 * t290 + t344 + t460) * qJD(1);
t239 = Icges(5,4) * t454;
t465 = Icges(5,5) * t292;
t155 = Icges(5,1) * t452 - t239 - t465;
t462 = Icges(5,6) * t292;
t153 = Icges(5,4) * t452 - Icges(5,2) * t454 - t462;
t458 = t153 * t289;
t345 = -t155 * t291 + t458;
t152 = Icges(5,3) * t290 + t215 * t292;
t418 = qJD(1) * t152;
t511 = qJD(1) * t345 - t290 * t327 + t418;
t151 = Icges(5,5) * t452 - Icges(5,6) * t454 - t460;
t59 = -t151 * t292 - t290 * t345;
t216 = Icges(5,2) * t291 + t469;
t340 = t216 * t289 - t218 * t291;
t510 = t340 * qJD(1) + qJD(4) * t215;
t509 = t290 * (-t216 * t292 + t156) - t292 * (-Icges(5,2) * t452 + t155 - t239);
t162 = Icges(6,3) * t289 + t291 * t349;
t342 = -t163 * t299 + t165 * t301;
t347 = -t103 * t299 + t106 * t301;
t508 = t199 * (-t161 * t292 - t347) - t200 * (-t161 * t290 + t522) + t258 * (t162 - t342);
t197 = (-Icges(6,2) * t301 - t467) * t289;
t507 = t199 * (-Icges(6,2) * t195 + t106 + t178) - t200 * (-Icges(6,2) * t193 - t105 - t176) + t258 * (t165 + t197);
t210 = qJDD(4) * t290 + t292 * t406;
t394 = t291 * t411;
t397 = t289 * t416;
t324 = t394 - t397;
t112 = qJD(5) * t324 + t292 * t405 + t210;
t506 = t112 / 0.2e1;
t505 = t113 / 0.2e1;
t504 = -t199 / 0.2e1;
t503 = t199 / 0.2e1;
t502 = -t200 / 0.2e1;
t501 = t200 / 0.2e1;
t500 = t204 / 0.2e1;
t499 = t210 / 0.2e1;
t498 = t211 / 0.2e1;
t497 = -t258 / 0.2e1;
t496 = t258 / 0.2e1;
t495 = t290 / 0.2e1;
t494 = -t292 / 0.2e1;
t493 = -rSges(6,3) - pkin(7);
t491 = g(1) * t290;
t44 = Icges(6,5) * t89 + Icges(6,6) * t88 + Icges(6,3) * t325;
t46 = Icges(6,4) * t89 + Icges(6,2) * t88 + Icges(6,6) * t325;
t48 = Icges(6,1) * t89 + Icges(6,4) * t88 + Icges(6,5) * t325;
t7 = (-qJD(4) * t522 - t44) * t291 + (qJD(4) * t98 - t299 * t46 + t301 * t48 + (-t101 * t301 + t105 * t299) * qJD(5)) * t289;
t489 = t7 * t200;
t86 = t292 * t319 + t379 * t451;
t87 = t292 * t318 - t379 * t450;
t43 = Icges(6,5) * t87 + Icges(6,6) * t86 + Icges(6,3) * t324;
t45 = Icges(6,4) * t87 + Icges(6,2) * t86 + Icges(6,6) * t324;
t47 = Icges(6,1) * t87 + Icges(6,4) * t86 + Icges(6,5) * t324;
t8 = (qJD(4) * t347 - t43) * t291 + (qJD(4) * t100 - t299 * t45 + t301 * t47 + (-t103 * t301 - t106 * t299) * qJD(5)) * t289;
t488 = t8 * t199;
t196 = (-Icges(6,5) * t299 - Icges(6,6) * t301) * t289;
t116 = qJD(4) * t162 + qJD(5) * t196;
t164 = Icges(6,6) * t289 + t291 * t350;
t117 = qJD(4) * t164 + qJD(5) * t197;
t166 = Icges(6,5) * t289 + t291 * t352;
t198 = (-Icges(6,1) * t299 - t466) * t289;
t118 = qJD(4) * t166 + qJD(5) * t198;
t22 = (qJD(4) * t342 - t116) * t291 + (qJD(4) * t161 - t117 * t299 + t118 * t301 + (-t163 * t301 - t165 * t299) * qJD(5)) * t289;
t66 = -t161 * t291 + t289 * t342;
t486 = t66 * t204 + t22 * t258;
t109 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t453;
t233 = pkin(7) * t394;
t396 = t289 * t411;
t326 = -t290 * t415 - t396;
t119 = pkin(4) * t326 - pkin(7) * t397 + t233;
t449 = t291 * t292;
t189 = pkin(4) * t449 + pkin(7) * t453;
t380 = t287 * t292 - t290 * t298;
t150 = t380 - t226;
t264 = qJ(3) * t414;
t378 = qJDD(1) * t293 - t303 * t492;
t419 = t264 + t270;
t313 = -qJDD(3) * t292 + qJD(1) * (-pkin(2) * t416 + t419) + qJDD(1) * t226 + t290 * t407 + t378;
t312 = qJD(1) * (-t264 + (t290 * t487 - t263) * qJD(1)) + qJDD(1) * t150 + t313;
t404 = rSges(6,1) * t87 + rSges(6,2) * t86 + rSges(6,3) * t394;
t49 = -rSges(6,3) * t397 + t404;
t11 = qJD(1) * t119 + qJDD(1) * t189 + t109 * t204 - t112 * t167 - t121 * t199 - t209 * t412 - t210 * t228 + t258 * t49 + t312;
t479 = t11 * t292;
t222 = rSges(5,1) * t289 + rSges(5,2) * t291;
t181 = t222 * t292;
t279 = t290 * rSges(5,3);
t158 = rSges(5,1) * t449 - rSges(5,2) * t453 + t279;
t376 = t150 + t383;
t64 = -t222 * t412 - t271 + (t158 + t376) * qJD(1);
t478 = t181 * t64;
t35 = -t228 * t412 + t109 * t258 - t167 * t199 - t271 + (t189 + t376) * qJD(1);
t474 = t290 * t35;
t157 = t402 + t516;
t338 = -t157 + t377;
t370 = -t222 * t411 + t270;
t63 = qJD(1) * t338 + t370;
t473 = t290 * t63;
t472 = t41 * t113;
t42 = -t100 * t291 + t289 * t347;
t471 = t42 * t112;
t456 = t214 * t290;
t455 = t214 * t292;
t80 = -t290 * t340 - t455;
t442 = t80 * qJD(1);
t439 = -t108 + t187;
t438 = t109 + t189;
t437 = -t121 - t209;
t436 = -t151 * t290 - t155 * t449;
t435 = t152 * t290 + t156 * t449;
t430 = t167 + t228;
t428 = -t216 + t219;
t427 = t218 + t351;
t399 = t289 * t481;
t426 = rSges(6,3) * t452 + t290 * t399;
t425 = rSges(6,3) * t449 + t292 * t399;
t424 = rSges(5,2) * t397 + rSges(5,3) * t414;
t253 = t290 * t482;
t422 = rSges(4,3) * t414 + qJD(1) * t253;
t421 = t252 + t271;
t420 = rSges(4,3) * t292 + t253;
t417 = qJD(1) * t215;
t408 = -m(4) - m(5) - m(6);
t403 = t290 * t485;
t401 = t289 * t484;
t393 = -pkin(2) - t485;
t391 = t416 / 0.2e1;
t390 = t414 / 0.2e1;
t389 = -t412 / 0.2e1;
t388 = t412 / 0.2e1;
t387 = -t411 / 0.2e1;
t386 = t411 / 0.2e1;
t385 = -t409 / 0.2e1;
t133 = t156 * t452;
t382 = t292 * t152 - t133;
t381 = -t151 + t457;
t159 = t403 - t420;
t375 = -t159 + t384;
t371 = t293 + t380;
t256 = rSges(2,1) * t302 - rSges(2,2) * t300;
t255 = rSges(2,1) * t300 + rSges(2,2) * t302;
t225 = rSges(5,1) * t291 - rSges(5,2) * t289;
t75 = t154 * t291 + t156 * t289;
t328 = qJD(4) * t216;
t92 = -t292 * t328 + (-t290 * t351 + t462) * qJD(1);
t329 = qJD(4) * t218;
t94 = -t292 * t329 + (-t219 * t290 + t465) * qJD(1);
t307 = -qJD(4) * t75 - t289 * t92 + t291 * t94 + t418;
t74 = t153 * t291 + t155 * t289;
t93 = qJD(1) * t154 - t290 * t328;
t95 = qJD(1) * t156 - t290 * t329;
t308 = qJD(1) * t151 - qJD(4) * t74 - t289 * t93 + t291 * t95;
t365 = -(t290 * t511 + t292 * t308) * t292 + (t290 * t512 + t292 * t307) * t290;
t364 = -(t290 * t308 - t292 * t511) * t292 + (t290 * t307 - t292 * t512) * t290;
t363 = t27 * t292 - t28 * t290;
t362 = t27 * t290 + t28 * t292;
t361 = t29 * t292 - t290 * t30;
t360 = t29 * t290 + t292 * t30;
t359 = t290 * t42 - t292 * t41;
t358 = t290 * t41 + t292 * t42;
t60 = -t154 * t454 - t382;
t357 = t290 * t60 - t292 * t59;
t61 = -t153 * t453 - t436;
t62 = -t154 * t453 + t435;
t356 = t290 * t62 - t292 * t61;
t355 = -t290 * t64 - t292 * t63;
t96 = rSges(5,1) * t326 - rSges(5,2) * t394 + t424;
t180 = t222 * t290;
t97 = -qJD(4) * t180 + (t225 * t292 + t279) * qJD(1);
t354 = t290 * t97 + t292 * t96;
t346 = -t108 * t292 - t109 * t290;
t343 = t157 * t290 + t158 * t292;
t341 = t216 * t291 + t218 * t289;
t168 = rSges(6,1) * t447 - rSges(6,2) * t448 + t278;
t335 = t289 * t493 - t285;
t323 = t100 * t199 + t161 * t258 - t200 * t98;
t322 = (-Icges(6,5) * t192 - Icges(6,6) * t193) * t200 - (Icges(6,5) * t194 - Icges(6,6) * t195) * t199 - t196 * t258;
t321 = t153 * t292 - t154 * t290;
t320 = t289 * t322;
t314 = (-t289 * t427 + t291 * t428) * qJD(1);
t311 = (Icges(6,1) * t194 - t103 - t468) * t199 - (-Icges(6,1) * t192 - t101 - t177) * t200 + (-t163 + t198) * t258;
t34 = qJD(1) * t337 - t527;
t36 = -t108 * t199 + t109 * t200 + qJD(2) + (t187 * t290 + t189 * t292) * qJD(4);
t309 = t36 * t346 + (t290 * t34 - t292 * t35) * t167;
t206 = t351 * qJD(4);
t207 = t219 * qJD(4);
t306 = qJD(1) * t214 - qJD(4) * t341 - t206 * t289 + t207 * t291;
t305 = -t289 * t509 + t291 * t321;
t304 = t508 * t289;
t249 = pkin(7) * t449;
t247 = pkin(7) * t452;
t208 = t225 * qJD(4);
t188 = -pkin(4) * t453 + t249;
t186 = -pkin(4) * t454 + t247;
t143 = -t292 * t401 + t425;
t142 = -t290 * t401 + t426;
t141 = t165 * t292;
t140 = t165 * t290;
t139 = t163 * t292;
t138 = t163 * t290;
t129 = rSges(6,1) * t194 - rSges(6,2) * t195;
t128 = -rSges(6,1) * t192 - rSges(6,2) * t193;
t111 = qJD(1) * t131 - t271;
t110 = qJD(1) * t375 + t270;
t81 = -t292 * t340 + t456;
t76 = t81 * qJD(1);
t73 = qJD(4) * t343 + qJD(2);
t55 = qJDD(1) * t160 + qJD(1) * (-qJD(1) * t403 + t422) + t313;
t54 = t375 * qJDD(1) + (-qJD(1) * t160 - t179) * qJD(1) + t339;
t40 = t290 * t306 - t292 * t510;
t39 = t290 * t510 + t292 * t306;
t38 = -qJD(4) * t344 + t289 * t94 + t291 * t92;
t37 = -qJD(4) * t345 + t289 * t95 + t291 * t93;
t33 = qJD(4) * t354 + t157 * t210 - t158 * t211 + qJDD(2);
t26 = qJD(1) * t96 + qJDD(1) * t158 - t208 * t412 - t210 * t222 + t312;
t25 = -t208 * t411 + t211 * t222 + (-t97 + t434) * qJD(1) + t338 * qJDD(1) + t339;
t24 = qJD(4) * t356 + t76;
t23 = qJD(4) * t357 + t442;
t16 = t116 * t454 - t117 * t192 + t118 * t193 + t161 * t325 + t163 * t88 + t165 * t89;
t15 = t116 * t453 + t117 * t194 + t118 * t195 + t161 * t324 + t163 * t86 + t165 * t87;
t14 = t199 * t42 - t200 * t41 + t258 * t66;
t9 = -t108 * t112 - t109 * t113 + t187 * t210 - t189 * t211 + t199 * t50 + t200 * t49 + qJDD(2) + (t119 * t292 + t120 * t290) * qJD(4);
t6 = t100 * t325 + t103 * t88 + t106 * t89 - t192 * t45 + t193 * t47 + t43 * t454;
t5 = t98 * t395 + t101 * t88 - t105 * t89 - t192 * t46 + t193 * t48 + (t290 * t44 + t414 * t98) * t289;
t4 = t100 * t324 + t103 * t86 + t106 * t87 + t194 * t45 + t195 * t47 + t43 * t453;
t3 = t98 * t394 + t101 * t86 - t105 * t87 + t194 * t46 + t195 * t48 + (t292 * t44 - t416 * t98) * t289;
t2 = t112 * t28 + t113 * t27 + t16 * t258 + t199 * t6 - t200 * t5 + t204 * t52;
t1 = t112 * t30 + t113 * t29 + t15 * t258 + t199 * t4 - t200 * t3 + t204 * t53;
t17 = [(-qJD(4) * t340 + t206 * t291 + t207 * t289) * qJD(1) + t486 + (t76 + ((t60 - t133 + (t152 + t458) * t292 + t436) * t292 + t435 * t290) * qJD(4)) * t386 - m(2) * (-g(1) * t255 + g(2) * t256) + t471 / 0.2e1 + t472 / 0.2e1 - t489 / 0.2e1 + t488 / 0.2e1 + t16 * t502 + t15 * t503 + t52 * t505 + t53 * t506 + (t75 + t81) * t499 + (t74 + t80) * t498 + (t23 - t442 + ((t292 * t381 - t435 + t62) * t292 + (t290 * t381 + t382 + t61) * t290) * qJD(4)) * t389 + (t38 + t39) * t388 + (t501 + t502) * t13 + ((-t224 * t303 - g(2) + t378) * t203 + (-t524 + (-0.2e1 * t227 - t293 + t203) * t303 - g(1)) * t202) * m(3) + (t37 + t40 + t24) * t387 + (t335 * t480 + (-t287 + t335) * t474 * qJD(1) - t335 * t491 + (-g(2) + t11) * (t371 + t438) + (-g(1) + t10) * (-t367 + t530) + (-t373 + t421 + (t291 * t493 + t490) * t412 + (-t293 + (-t287 - t515 - t278) * t292) * qJD(1)) * t34 + (-pkin(4) * t396 + t233 + t34 + t404 + (t187 - t263) * qJD(1) + t527 + t529) * t35) * m(6) + ((t222 * t473 - t478) * qJD(4) + (t421 + (-t279 - t293 + (-t225 - t287) * t292) * qJD(1)) * t63 + (t63 - t370 + t424 + (t157 + t492 + t528) * qJD(1) + t529) * t64 + (-g(2) + t26) * (t158 + t371) + (-g(1) + t25) * (-t516 + t528)) * m(5) + (t110 * t271 + t111 * (t419 + t422) + ((-t110 * t302 - t111 * t300) * pkin(1) + t110 * (t393 + t482) * t292 + (t110 * (-rSges(4,3) - qJ(3)) + t111 * t393) * t290) * qJD(1) - (-t110 - t213 + t270 + (-t159 - t492) * qJD(1)) * t111 + (-g(2) + t55) * t131 + (-g(1) + t54) * (t290 * t393 + t273 + t420 - t492)) * m(4) + (t341 + m(2) * (t255 ^ 2 + t256 ^ 2) + m(3) * (t202 ^ 2 + t227 * t203) + Icges(4,2) * t297 ^ 2 + (Icges(4,1) * t296 + 0.2e1 * Icges(4,4) * t297) * t296 + Icges(2,3) + Icges(3,3)) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t33 + m(6) * t9 + (-m(3) + t408) * g(3); t408 * (-g(2) * t292 + t491) + 0.2e1 * (t480 / 0.2e1 - t479 / 0.2e1) * m(6) + 0.2e1 * (t25 * t495 + t26 * t494) * m(5) + 0.2e1 * (t494 * t55 + t495 * t54) * m(4); t24 * t390 + t23 * t391 + ((t290 * t61 + t62 * t292) * qJD(1) + t365) * t388 + ((t290 * t59 + t292 * t60) * qJD(1) + t364) * t387 + (t292 * t385 + t390) * t13 - t112 * t361 / 0.2e1 + (g(1) * t181 + g(2) * t180 - g(3) * t225 - (t180 * t63 - t478) * qJD(1) - (t73 * (-t180 * t290 - t181 * t292) + t355 * t225) * qJD(4) + t33 * t343 + t73 * ((t157 * t292 - t158 * t290) * qJD(1) + t354) + t355 * t208 + (-t25 * t292 - t26 * t290 + (-t292 * t64 + t473) * qJD(1)) * t222) * m(5) + ((-t412 * t455 + t417) * t290 + (t314 + (t290 * t456 + t305) * qJD(4)) * t292) * t389 + ((-t411 * t456 - t417) * t292 + (t314 + (t292 * t455 + t305) * qJD(4)) * t290) * t386 + ((t139 * t192 - t141 * t193) * t199 - (t138 * t192 - t140 * t193) * t200 + (-t164 * t192 + t166 * t193) * t258 + (t28 * t449 + t289 * t52) * qJD(5) + ((qJD(5) * t27 + t323) * t291 + t304) * t290) * t501 + ((-t139 * t194 - t141 * t195) * t199 - (-t138 * t194 - t140 * t195) * t200 + (t164 * t194 + t166 * t195) * t258 + (t289 * t53 + t29 * t452) * qJD(5) + ((qJD(5) * t30 + t323) * t291 + t304) * t292) * t504 + (t290 * t385 + t391) * t12 - qJD(1) * ((t289 * t428 + t291 * t427) * qJD(1) + (t321 * t289 + t291 * t509) * qJD(4)) / 0.2e1 + (qJD(1) * t40 + qJD(4) * t364 + qJDD(1) * t80 + t210 * t60 + t211 * t59 + t2) * t494 + (qJD(1) * t358 + t290 * t8 - t292 * t7) * t496 + t357 * t498 + t356 * t499 + t359 * t500 + (qJD(1) * t362 + t290 * t6 - t292 * t5) * t502 + (qJD(1) * t360 + t290 * t4 - t292 * t3) * t503 + (qJD(1) * t39 + qJD(4) * t365 + qJDD(1) * t81 + t210 * t62 + t211 * t61 + t1) * t495 + qJD(1) * (t290 * t38 - t292 * t37 + (t74 * t290 + t292 * t75) * qJD(1)) / 0.2e1 - t14 * t410 / 0.2e1 + (((t139 * t299 - t141 * t301 + t100) * t199 - (t138 * t299 - t140 * t301 + t98) * t200 + (-t164 * t299 + t166 * t301 + t161) * t258 + t66 * qJD(5)) * t289 + (qJD(5) * t358 - t508) * t291) * t497 - t113 * t363 / 0.2e1 + qJDD(1) * (t290 * t75 - t292 * t74) / 0.2e1 + (-g(1) * (t249 + t425) - g(2) * (t247 + t426) - g(3) * (t168 + t515) - (g(1) * t292 + g(2) * t290) * t289 * (-pkin(4) - t484) + (-t10 * t430 + t34 * t437 + t9 * t438 + t36 * (t119 + t49) + (-t35 * t430 + t36 * t439) * qJD(1)) * t292 + (-t11 * t430 + t35 * t437 + t9 * t439 + t36 * (t120 + t50) + (t34 * t430 - t36 * t438) * qJD(1)) * t290 - t34 * (-qJD(1) * t186 - t142 * t258 - t168 * t200 - t411 * t515) - t35 * (qJD(1) * t188 + t143 * t258 - t168 * t199 - t412 * t515) - t36 * (t142 * t199 + t143 * t200 + t186 * t412 + t188 * t411) - ((t108 * t34 + t109 * t35) * t289 + t309 * t291) * qJD(5)) * m(6); t1 * t453 / 0.2e1 + (t289 * t360 - t291 * t53) * t506 + ((qJD(4) * t360 - t15) * t291 + (qJD(1) * t361 + qJD(4) * t53 + t290 * t3 + t292 * t4) * t289) * t503 + t2 * t454 / 0.2e1 + (t289 * t362 - t291 * t52) * t505 + ((qJD(4) * t362 - t16) * t291 + (qJD(1) * t363 + qJD(4) * t52 + t290 * t5 + t292 * t6) * t289) * t502 + t14 * t413 / 0.2e1 - t291 * (t471 + t472 + t486 + t488 - t489) / 0.2e1 + (t289 * t358 - t291 * t66) * t500 + ((qJD(4) * t358 - t22) * t291 + (-qJD(1) * t359 + qJD(4) * t66 + t290 * t7 + t292 * t8) * t289) * t496 + (t194 * t507 + t311 * t195 - t292 * t320) * t504 + (-t192 * t507 + t193 * t311 - t290 * t320) * t501 + (t322 * t291 + (-t299 * t507 + t301 * t311) * t289) * t497 + (-t397 / 0.2e1 + t291 * t386) * t13 + (t289 * t390 + t291 * t388) * t12 + ((qJD(4) * t309 - t10 * t108 - t11 * t109 + t34 * t50 - t35 * t49) * t291 + (t34 * (qJD(4) * t108 + t121 * t290) + t35 * (qJD(4) * t109 - t121 * t292) + t9 * t346 + t36 * (t108 * t416 - t109 * t414 - t290 * t49 + t292 * t50) + (t480 - t479 + (t292 * t34 + t474) * qJD(1)) * t167) * t289 - t34 * (-t128 * t258 - t200 * t201) - t35 * (t129 * t258 - t199 * t201) - t36 * (t128 * t199 + t129 * t200) - g(1) * t129 - g(2) * t128 - g(3) * t201) * m(6);];
tau = t17;
