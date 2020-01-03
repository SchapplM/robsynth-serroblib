% Calculate vector of inverse dynamics joint torques for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:54:06
% DurationCPUTime: 12.23s
% Computational Cost: add. (18888->627), mult. (11804->804), div. (0->0), fcn. (9175->10), ass. (0->355)
t299 = qJ(1) + pkin(9);
t286 = sin(t299);
t287 = cos(t299);
t209 = rSges(3,1) * t287 - t286 * rSges(3,2);
t304 = cos(qJ(1));
t295 = t304 * pkin(1);
t538 = t209 + t295;
t298 = qJD(1) + qJD(3);
t288 = qJ(3) + t299;
t278 = sin(t288);
t279 = cos(t288);
t523 = t278 * rSges(4,1) + t279 * rSges(4,2);
t175 = t523 * t298;
t302 = sin(qJ(1));
t293 = t302 * pkin(1);
t522 = pkin(2) * t286 + t293;
t339 = t522 * qJD(1);
t140 = t339 + t175;
t300 = qJ(4) + qJ(5);
t289 = sin(t300);
t290 = cos(t300);
t473 = Icges(6,4) * t289;
t217 = Icges(6,1) * t290 - t473;
t146 = -Icges(6,5) * t279 + t217 * t278;
t444 = t279 * t289;
t225 = Icges(6,4) * t444;
t443 = t279 * t290;
t147 = Icges(6,1) * t443 + Icges(6,5) * t278 - t225;
t297 = qJD(4) + qJD(5);
t199 = t278 * t297;
t200 = t279 * t297;
t214 = Icges(6,2) * t290 + t473;
t276 = Icges(6,4) * t290;
t353 = -Icges(6,2) * t289 + t276;
t521 = Icges(6,1) * t289 + t276;
t533 = t521 + t353;
t311 = t199 * (-Icges(6,2) * t443 + t147 - t225) - t200 * (-t214 * t278 + t146) + t298 * t533;
t144 = -Icges(6,6) * t279 + t278 * t353;
t145 = Icges(6,4) * t443 - Icges(6,2) * t444 + Icges(6,6) * t278;
t500 = t199 * (t279 * t521 + t145) - t200 * (t278 * t521 + t144) + t298 * (t214 - t217);
t537 = t311 * t289 + t290 * t500;
t301 = sin(qJ(4));
t303 = cos(qJ(4));
t291 = Icges(5,4) * t303;
t354 = -Icges(5,2) * t301 + t291;
t520 = Icges(5,1) * t301 + t291;
t414 = t520 + t354;
t474 = Icges(5,4) * t301;
t253 = Icges(5,2) * t303 + t474;
t256 = Icges(5,1) * t303 - t474;
t415 = t253 - t256;
t536 = (t301 * t414 + t303 * t415) * t298;
t257 = rSges(5,1) * t301 + rSges(5,2) * t303;
t186 = t257 * t278;
t187 = t257 * t279;
t446 = t278 * t301;
t398 = rSges(5,2) * t446;
t445 = t278 * t303;
t364 = rSges(5,1) * t445 - t398;
t161 = -rSges(5,3) * t279 + t364;
t270 = t278 * pkin(3);
t203 = -pkin(7) * t279 + t270;
t406 = qJD(4) * t279;
t510 = t257 * t406 + t298 * (t161 + t203);
t78 = t339 + t510;
t407 = qJD(4) * t278;
t389 = t257 * t407;
t408 = qJD(1) * t287;
t481 = pkin(1) * qJD(1);
t413 = pkin(2) * t408 + t304 * t481;
t343 = -t389 + t413;
t441 = t279 * t301;
t249 = rSges(5,2) * t441;
t440 = t279 * t303;
t401 = rSges(5,1) * t440;
t162 = rSges(5,3) * t278 - t249 + t401;
t484 = pkin(3) * t279;
t204 = pkin(7) * t278 + t484;
t424 = t162 + t204;
t79 = t298 * t424 + t343;
t534 = -t186 * t78 - t187 * t79;
t402 = -qJDD(4) - qJDD(5);
t131 = -t200 * t298 + t278 * t402;
t442 = t279 * t298;
t243 = pkin(7) * t442;
t447 = t278 * t298;
t177 = pkin(3) * t447 - t243;
t405 = qJD(4) * t298;
t188 = -qJDD(4) * t278 - t279 * t405;
t277 = t290 * rSges(6,1);
t219 = -rSges(6,2) * t289 + t277;
t194 = t219 * t297;
t218 = rSges(6,1) * t289 + rSges(6,2) * t290;
t296 = qJDD(1) + qJDD(3);
t464 = pkin(1) * qJDD(1);
t282 = t304 * t464;
t307 = qJD(1) ^ 2;
t463 = pkin(2) * qJDD(1);
t317 = t287 * t463 - t307 * t522 + t282;
t294 = t303 * pkin(4);
t280 = t294 + pkin(3);
t227 = t279 * t280;
t305 = -pkin(8) - pkin(7);
t138 = t484 - t227 + (pkin(7) + t305) * t278;
t400 = rSges(6,1) * t443;
t365 = -rSges(6,2) * t444 + t400;
t150 = rSges(6,3) * t278 + t365;
t434 = -t138 + t150;
t391 = t204 + t434;
t436 = t303 * qJD(4) ^ 2;
t261 = t279 * t305;
t404 = qJD(4) * t301;
t385 = t279 * t404;
t367 = pkin(4) * t385;
t112 = t367 + t243 + (t261 + (-pkin(3) + t280) * t278) * t298;
t231 = rSges(6,2) * t443;
t438 = t289 * t298;
t397 = rSges(6,2) * t438;
t423 = rSges(6,3) * t442 + t278 * t397;
t439 = t289 * t297;
t94 = t297 * t231 + (t279 * t439 + t290 * t447) * rSges(6,1) - t423;
t477 = -t112 - t94;
t18 = t131 * t218 - t194 * t199 + (t188 * t301 - t278 * t436) * pkin(4) + (-t177 + t477) * t298 + t391 * t296 + t317;
t532 = -g(2) + t18;
t205 = t280 * t442;
t396 = pkin(4) * t404;
t341 = -t298 * t305 - t396;
t416 = pkin(3) * t442 + pkin(7) * t447;
t113 = t278 * t341 + t205 - t416;
t236 = t278 * t405;
t132 = qJD(5) * t447 + t279 * t402 + t236;
t189 = -qJDD(4) * t279 + t236;
t275 = pkin(2) * t287;
t410 = t307 * t295 + t302 * t464;
t366 = t307 * t275 + t286 * t463 + t410;
t342 = t296 * t203 + t298 * t416 + t366;
t417 = t278 * t280 + t261;
t137 = -t203 + t417;
t448 = t278 * t290;
t449 = t278 * t289;
t149 = rSges(6,1) * t448 - rSges(6,2) * t449 - rSges(6,3) * t279;
t435 = t137 + t149;
t399 = rSges(6,1) * t439;
t422 = rSges(6,3) * t447 + t298 * t400;
t437 = t290 * t297;
t95 = -t278 * t399 + (-t278 * t437 - t279 * t438) * rSges(6,2) + t422;
t17 = -t132 * t218 + t194 * t200 + (t113 + t95) * t298 + t435 * t296 + (-t189 * t301 + t279 * t436) * pkin(4) + t342;
t531 = -g(3) + t17;
t403 = qJD(4) * t303;
t386 = t278 * t403;
t387 = t278 * t404;
t393 = t298 * t441;
t420 = rSges(5,3) * t447 + t298 * t401;
t111 = -rSges(5,1) * t387 + (-t386 - t393) * rSges(5,2) + t420;
t483 = rSges(5,1) * t303;
t259 = -rSges(5,2) * t301 + t483;
t232 = t259 * qJD(4);
t47 = t111 * t298 + t161 * t296 - t189 * t257 + t232 * t406 + t342;
t530 = t47 - g(3);
t394 = t298 * t445;
t421 = -rSges(5,3) * t442 - t298 * t398;
t110 = rSges(5,2) * t279 * t403 + (t385 + t394) * rSges(5,1) + t421;
t48 = -t232 * t407 + t188 * t257 + (-t110 - t177) * t298 + t424 * t296 + t317;
t529 = t48 - g(2);
t176 = rSges(4,1) * t442 - rSges(4,2) * t447;
t528 = t176 * t298 + t296 * t523 - g(3) + t366;
t202 = rSges(4,1) * t279 - t278 * rSges(4,2);
t527 = -t175 * t298 + t202 * t296 - g(2) + t317;
t213 = Icges(6,5) * t290 - Icges(6,6) * t289;
t142 = -Icges(6,3) * t279 + t213 * t278;
t458 = t145 * t289;
t352 = -t147 * t290 + t458;
t340 = -t142 + t352;
t525 = t200 * t340;
t208 = t286 * rSges(3,1) + t287 * rSges(3,2);
t139 = t298 * t150;
t195 = t298 * t204;
t331 = -pkin(4) * t387 - t199 * t218;
t518 = t298 * t138 - t139 - t195 + t205 - t331 + t422;
t451 = t214 * t297;
t516 = -Icges(6,6) * t298 + t451;
t148 = t298 * t162;
t515 = -rSges(5,2) * t393 - t148 - t195 + t416 + t420;
t172 = t218 * t278;
t173 = rSges(6,1) * t444 + t231;
t53 = t149 * t199 + t150 * t200 + qJD(2) + (t137 * t278 - t138 * t279) * qJD(4);
t344 = -t200 * t218 - t367;
t392 = t203 + t435;
t60 = t298 * t392 + t339 - t344;
t61 = t298 * t391 + t331 + t413;
t514 = -(-t298 * t172 + t200 * t219) * t60 - t53 * (-t172 * t199 - t200 * t173) - t61 * (-t173 * t298 - t199 * t219);
t336 = t354 * t298;
t508 = -Icges(5,6) * t298 + qJD(4) * t253;
t107 = -t278 * t508 + t279 * t336;
t338 = t256 * t298;
t505 = -Icges(5,5) * t298 + qJD(4) * t520;
t109 = -t278 * t505 + t279 * t338;
t252 = Icges(5,5) * t303 - Icges(5,6) * t301;
t154 = -Icges(5,3) * t279 + t252 * t278;
t156 = -Icges(5,6) * t279 + t278 * t354;
t158 = -Icges(5,5) * t279 + t256 * t278;
t99 = t156 * t303 + t158 * t301;
t513 = qJD(4) * t99 + t107 * t301 - t109 * t303 - t154 * t298;
t212 = Icges(6,5) * t289 + Icges(6,6) * t290;
t512 = -Icges(6,3) * t298 + t212 * t297;
t511 = -Icges(6,5) * t298 + t297 * t521;
t251 = Icges(5,5) * t301 + Icges(5,6) * t303;
t509 = -Icges(5,3) * t298 + qJD(4) * t251;
t221 = t354 * qJD(4);
t222 = t256 * qJD(4);
t346 = t253 * t303 + t301 * t520;
t507 = qJD(4) * t346 + t221 * t301 - t222 * t303 - t251 * t298;
t106 = t278 * t336 + t279 * t508;
t108 = t278 * t338 + t279 * t505;
t155 = Icges(5,5) * t440 - Icges(5,6) * t441 + Icges(5,3) * t278;
t157 = Icges(5,4) * t440 - Icges(5,2) * t441 + Icges(5,6) * t278;
t247 = Icges(5,4) * t441;
t159 = Icges(5,1) * t440 + Icges(5,5) * t278 - t247;
t350 = t157 * t303 + t159 * t301;
t506 = qJD(4) * t350 - t106 * t301 + t108 * t303 - t155 * t298;
t368 = -t217 * t297 + t451;
t369 = t533 * t297;
t503 = -t212 * t298 + t289 * t369 + t290 * t368;
t143 = Icges(6,5) * t443 - Icges(6,6) * t444 + Icges(6,3) * t278;
t335 = t353 * t298;
t373 = t147 * t297 - t278 * t335 - t279 * t516;
t337 = t217 * t298;
t375 = t145 * t297 + t278 * t337 + t279 * t511;
t502 = -t143 * t298 + t289 * t373 + t290 * t375;
t374 = t146 * t297 - t278 * t516 + t279 * t335;
t376 = t144 * t297 + t278 * t511 - t279 * t337;
t501 = -t142 * t298 + t289 * t374 + t290 * t376;
t499 = m(3) + m(4);
t498 = t131 / 0.2e1;
t497 = t132 / 0.2e1;
t496 = t188 / 0.2e1;
t495 = t189 / 0.2e1;
t494 = -t199 / 0.2e1;
t493 = t199 / 0.2e1;
t492 = t200 / 0.2e1;
t491 = -t200 / 0.2e1;
t490 = -t278 / 0.2e1;
t489 = -t279 / 0.2e1;
t488 = t296 / 0.2e1;
t487 = -t298 / 0.2e1;
t486 = t298 / 0.2e1;
t485 = rSges(5,3) + pkin(7);
t479 = t298 * t61;
t453 = t212 * t278;
t98 = -t214 * t444 + t443 * t521 + t453;
t478 = t98 * t298;
t450 = t251 * t278;
t116 = -t253 * t441 + t440 * t520 + t450;
t462 = t116 * t298;
t459 = t144 * t289;
t457 = t146 * t290;
t456 = t156 * t301;
t455 = t157 * t301;
t454 = t158 * t303;
t167 = t212 * t279;
t333 = t213 * t298;
t181 = t251 * t279;
t334 = t252 * t298;
t429 = t278 * t520 + t156;
t428 = t279 * t520 + t157;
t427 = -t253 * t278 + t158;
t426 = -Icges(5,2) * t440 + t159 - t247;
t411 = t275 + t295;
t409 = qJD(1) * t286;
t395 = t149 * t442 + (-t139 + t95) * t278;
t383 = t447 / 0.2e1;
t382 = -t442 / 0.2e1;
t381 = -t407 / 0.2e1;
t380 = t407 / 0.2e1;
t379 = -t406 / 0.2e1;
t378 = t406 / 0.2e1;
t377 = -pkin(4) * t301 - t218;
t372 = -t143 - t459;
t371 = -t143 + t457;
t141 = t202 * t298 + t413;
t260 = rSges(2,1) * t304 - t302 * rSges(2,2);
t258 = rSges(2,1) * t302 + rSges(2,2) * t304;
t134 = t158 * t445;
t69 = -t154 * t279 - t156 * t446 + t134;
t135 = t159 * t445;
t70 = t155 * t279 + t157 * t446 - t135;
t359 = -t278 * t70 - t279 * t69;
t136 = t156 * t441;
t71 = -t154 * t278 - t158 * t440 + t136;
t349 = -t159 * t303 + t455;
t72 = t278 * t155 - t279 * t349;
t358 = -t278 * t72 - t279 * t71;
t357 = -t278 * t79 + t279 * t78;
t82 = -t145 * t290 - t147 * t289;
t351 = t454 - t456;
t348 = t161 * t278 + t162 * t279;
t347 = -t214 * t289 + t290 * t521;
t345 = -t253 * t301 + t303 * t520;
t332 = -pkin(2) * t409 - t302 * t481;
t328 = t167 * t199 - t200 * t453 - t333;
t326 = -t278 * t333 - t279 * t512 + t298 * t352;
t325 = -t279 * t333 + t512 * t278 + (t457 - t459) * t298;
t324 = -t278 * t334 - t279 * t509 + t298 * t349;
t323 = t278 * t509 - t279 * t334 + t298 * t351;
t322 = -t213 * t297 + t298 * t347;
t321 = -t252 * qJD(4) + t298 * t345;
t320 = t301 * t427 + t303 * t429;
t319 = t301 * t426 + t303 * t428;
t119 = t149 + t417;
t124 = -t279 * t485 + t270 + t364;
t120 = t227 + (rSges(6,3) - t305) * t278 + t365;
t125 = -t249 + (pkin(3) + t483) * t279 + t485 * t278;
t13 = t325 * t278 + t279 * t501;
t14 = t326 * t278 - t279 * t502;
t15 = -t278 * t501 + t325 * t279;
t16 = t278 * t502 + t326 * t279;
t121 = t146 * t448;
t63 = -t142 * t279 - t144 * t449 + t121;
t122 = t147 * t448;
t64 = t143 * t279 + t145 * t449 - t122;
t97 = t278 * t347 - t167;
t96 = t97 * t298;
t30 = -t199 * t64 - t200 * t63 + t96;
t123 = t144 * t444;
t65 = -t142 * t278 - t146 * t443 + t123;
t66 = t143 * t278 - t279 * t352;
t31 = -t199 * t66 - t200 * t65 - t478;
t40 = -t289 * t376 + t290 * t374;
t41 = t289 * t375 - t290 * t373;
t45 = t322 * t278 + t279 * t503;
t46 = -t278 * t503 + t322 * t279;
t81 = t144 * t290 + t146 * t289;
t315 = (-t13 * t200 + t131 * t66 + t132 * t65 - t14 * t199 - t296 * t98 + t298 * t45) * t490 + (t328 * t278 + t537 * t279) * t493 + (-t537 * t278 + t328 * t279) * t492 + (t131 * t64 + t132 * t63 - t15 * t200 - t16 * t199 + t296 * t97 + t298 * t46) * t489 + (-t289 * t500 + t290 * t311) * t487 + t30 * t383 + t31 * t382 + ((-t298 * t66 - t13) * t279 + (t298 * t65 - t14) * t278) * t494 + (-t278 * t66 - t279 * t65) * t498 + (-t278 * t64 - t279 * t63) * t497 + ((-t298 * t64 - t15) * t279 + (t298 * t63 - t16) * t278) * t491 + (-t278 * t82 - t279 * t81) * t488 + ((-t298 * t82 - t40) * t279 + (t298 * t81 - t41) * t278) * t486;
t314 = -rSges(5,1) * t394 - t177 - t421;
t310 = t534 * qJD(4);
t115 = t278 * t345 - t181;
t114 = t115 * t298;
t36 = qJD(4) * t359 + t114;
t37 = qJD(4) * t358 - t462;
t51 = qJD(4) * t351 + t107 * t303 + t109 * t301;
t52 = qJD(4) * t349 + t106 * t303 + t108 * t301;
t57 = t321 * t278 + t279 * t507;
t58 = -t278 * t507 + t321 * t279;
t309 = -t350 * t496 + t82 * t498 + (t114 + ((t71 + t135 - t136 + (t154 - t455) * t278) * t278 + (-t134 - t72 + (t154 - t349) * t279 + (t454 + t456) * t278) * t279) * qJD(4)) * t380 + (t96 - (t278 * t372 + t121 + t66) * t200 + (t122 - t123 + t65 + (t142 - t458) * t278) * t199 + (t199 * t371 - t525) * t279) * t493 - t131 * t98 / 0.2e1 - t188 * t116 / 0.2e1 + (t81 + t97) * t497 + (t115 + t99) * t495 + (t31 + t478 - (-t123 + t64) * t200 + (-t121 + t63) * t199 + (-t199 * t340 - t200 * t371) * t279 + (-t199 * t372 + t525) * t278) * t492 + (t40 + t46) * t491 + (t51 + t58) * t379 + (t37 + t462 + ((t136 - t70 + (t155 - t454) * t279) * t279 + (-t134 + t69 + (t155 + t456) * t278) * t278) * qJD(4)) * t378 + (t41 + t45 + t30) * t494 + (t52 + t57 + t36) * t381 + (qJD(4) * t345 + t221 * t303 + t222 * t301 - t289 * t368 + t290 * t369) * t298 + (t214 * t290 + t289 * t521 + Icges(4,3) + t346) * t296;
t308 = (t60 * (-t218 * t297 - t396) + (t61 * (-t280 - t277) - t60 * t305) * t298) * t278 + (t61 * (-rSges(6,2) * t437 + t341 - t399) - t60 * t397) * t279;
t250 = pkin(4) * t441;
t133 = t278 * t149;
t85 = qJD(4) * t348 + qJD(2);
t44 = -t161 * t188 - t162 * t189 + qJDD(2) + (-t110 * t279 + t111 * t278) * qJD(4);
t22 = t278 * t506 + t324 * t279;
t21 = -t278 * t513 + t323 * t279;
t20 = t324 * t278 - t279 * t506;
t19 = t323 * t278 + t279 * t513;
t10 = -t131 * t149 - t132 * t150 - t137 * t188 + t138 * t189 + t199 * t95 - t200 * t94 + qJDD(2) + (-t112 * t279 + t113 * t278) * qJD(4);
t1 = [t309 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t527 * (t202 + t411) + t528 * (t522 + t523) + (-t141 + t176 + t413) * t140) * m(4) + ((qJDD(1) * t209 - g(2) + t282) * t538 + (-t307 * t538 + qJDD(1) * t208 - g(3) + t410 + (0.2e1 * rSges(3,1) * t408 - 0.2e1 * rSges(3,2) * t409 - qJD(1) * t209) * qJD(1)) * (t293 + t208)) * m(3) + ((t258 ^ 2 + t260 ^ 2) * qJDD(1) - g(2) * t260 - g(3) * t258) * m(2) + (t61 * (t332 + t423) + t308 + t532 * (t120 + t411) + t531 * (t119 + t522) + (t61 + t518) * t60) * m(6) + (t79 * (t314 + t332) + t310 + (t413 - t343 + t79 + t515) * t78 + t529 * (t125 + t411) + t530 * (t124 + t522)) * m(5); t499 * qJDD(2) + m(5) * t44 + m(6) * t10 + (-m(5) - m(6) - t499) * g(1); t309 + (t392 * t479 + t308 + (-t344 + t423) * t61 + t518 * t60 + t532 * t120 + t531 * t119) * m(6) + (t310 + (t314 + t510) * t79 + (t389 + t515) * t78 + t529 * t125 + t530 * t124) * m(5) + (t140 * t176 - t141 * t175 + (-t140 * t298 + t527) * t202 + (t141 * t298 + t528) * t523) * m(4); ((-t406 * t450 - t334) * t279 + (-t536 + (-t319 * t278 + (t181 + t320) * t279) * qJD(4)) * t278) * t378 + t315 + ((-t301 * t415 + t303 * t414) * t298 + ((t278 * t426 - t279 * t427) * t303 + (-t278 * t428 + t279 * t429) * t301) * qJD(4)) * t487 + ((-t298 * t70 - t21) * t279 + (t298 * t69 - t22) * t278) * t379 + ((-t298 * t72 - t19) * t279 + (t298 * t71 - t20) * t278) * t381 + (-t116 * t296 + t188 * t72 + t189 * t71 + t298 * t57 + (-t19 * t279 - t20 * t278) * qJD(4)) * t490 + (t115 * t296 + t188 * t70 + t189 * t69 + t298 * t58 + (-t21 * t279 - t22 * t278) * qJD(4)) * t489 + ((t181 * t407 - t334) * t278 + (t536 + (-t320 * t279 + (-t450 + t319) * t278) * qJD(4)) * t279) * t380 + ((t298 * t350 - t51) * t279 + (t298 * t99 - t52) * t278) * t486 + t358 * t496 + t359 * t495 + (t278 * t350 - t279 * t99) * t488 + t36 * t383 + t37 * t382 + (t10 * t133 + t53 * t395 + t17 * t250 + (t10 * t434 + t53 * t477 + t17 * t218 + t60 * t194 + (t53 * t137 + t377 * t61) * t298) * t279 + (t10 * t137 + t53 * t113 + t18 * t377 + t61 * (-pkin(4) * t403 - t194) + (t53 * t138 + t377 * t60) * t298) * t278 - g(1) * (t219 + t294) - g(3) * (t250 + t173) - g(2) * t377 * t278 - (-t61 * t386 + ((-t278 * t60 - t279 * t61) * t298 + t53 * (-t278 ^ 2 - t279 ^ 2) * qJD(4)) * t301) * pkin(4) + t514) * m(6) + (-t534 * t298 - (t85 * (-t186 * t278 - t187 * t279) + t357 * t259) * qJD(4) + t44 * t348 + t85 * ((t161 * t298 - t110) * t279 + (t111 - t148) * t278) + t357 * t232 + ((-t298 * t79 + t47) * t279 + (-t298 * t78 - t48) * t278) * t257 - g(1) * t259 + g(2) * t186 - g(3) * t187) * m(5); t315 + (t10 * (t150 * t279 + t133) + t53 * (-t279 * t94 + t395) + (-t278 * t61 + t279 * t60) * t194 + ((t17 - t479) * t279 + (-t298 * t60 - t18) * t278) * t218 - g(1) * t219 + g(2) * t172 - g(3) * t173 + t514) * m(6);];
tau = t1;
