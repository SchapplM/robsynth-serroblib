% Calculate time derivative of joint inertia matrix for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:48:00
% EndTime: 2019-03-09 03:48:29
% DurationCPUTime: 18.21s
% Computational Cost: add. (38131->830), mult. (55670->1174), div. (0->0), fcn. (57852->10), ass. (0->389)
t291 = pkin(10) + qJ(3);
t285 = sin(t291);
t286 = cos(t291);
t456 = Icges(5,5) * t286;
t460 = Icges(4,4) * t286;
t518 = -t456 + t460 + (Icges(4,1) + Icges(5,1)) * t285;
t517 = t518 * qJD(3);
t477 = sin(qJ(5));
t478 = cos(qJ(5));
t242 = t285 * t478 - t286 * t477;
t298 = sin(qJ(1));
t300 = cos(qJ(1));
t324 = -t285 * t477 - t286 * t478;
t226 = t324 * t298;
t297 = sin(qJ(6));
t299 = cos(qJ(6));
t193 = t226 * t297 + t299 * t300;
t194 = -t226 * t299 + t297 * t300;
t225 = t242 * t298;
t112 = Icges(7,5) * t194 + Icges(7,6) * t193 - Icges(7,3) * t225;
t114 = Icges(7,4) * t194 + Icges(7,2) * t193 - Icges(7,6) * t225;
t116 = Icges(7,1) * t194 + Icges(7,4) * t193 - Icges(7,5) * t225;
t228 = t324 * t300;
t195 = t228 * t297 - t298 * t299;
t196 = -t228 * t299 - t297 * t298;
t227 = t242 * t300;
t48 = -t112 * t227 + t114 * t195 + t116 * t196;
t113 = Icges(7,5) * t196 + Icges(7,6) * t195 - Icges(7,3) * t227;
t115 = Icges(7,4) * t196 + Icges(7,2) * t195 - Icges(7,6) * t227;
t117 = Icges(7,1) * t196 + Icges(7,4) * t195 - Icges(7,5) * t227;
t49 = -t113 * t227 + t115 * t195 + t117 * t196;
t27 = t298 * t49 - t300 * t48;
t155 = -Icges(6,5) * t226 + Icges(6,6) * t225 + Icges(6,3) * t300;
t157 = -Icges(6,4) * t226 + Icges(6,2) * t225 + Icges(6,6) * t300;
t159 = -Icges(6,1) * t226 + Icges(6,4) * t225 + Icges(6,5) * t300;
t318 = -t155 * t298 + t157 * t227 - t159 * t228;
t158 = -Icges(6,4) * t228 + Icges(6,2) * t227 - Icges(6,6) * t298;
t160 = -Icges(6,1) * t228 + Icges(6,4) * t227 - Icges(6,5) * t298;
t156 = -Icges(6,5) * t228 + Icges(6,6) * t227 - Icges(6,3) * t298;
t450 = t156 * t298;
t69 = t158 * t227 - t160 * t228 - t450;
t498 = t298 * t69 - t300 * t318 + t27;
t46 = -t112 * t225 + t114 * t193 + t116 * t194;
t47 = -t113 * t225 + t115 * t193 + t117 * t194;
t26 = t298 * t47 - t300 * t46;
t319 = t155 * t300 + t157 * t225 - t159 * t226;
t449 = t156 * t300;
t68 = t158 * t225 - t160 * t226 + t449;
t499 = t298 * t68 - t300 * t319 + t26;
t516 = (t298 * t499 + t300 * t498) * qJD(1);
t515 = -t298 / 0.2e1;
t501 = t298 / 0.2e1;
t481 = -t300 / 0.2e1;
t514 = t300 / 0.2e1;
t513 = -qJD(1) / 0.2e1;
t512 = qJD(1) / 0.2e1;
t187 = rSges(6,1) * t242 + rSges(6,2) * t324;
t292 = t298 ^ 2;
t293 = t300 ^ 2;
t427 = t292 + t293;
t511 = t187 * t427;
t392 = qJD(5) * t477;
t393 = qJD(5) * t478;
t394 = qJD(3) * t477;
t395 = qJD(3) * t478;
t182 = (-t392 + t394) * t286 + (t393 - t395) * t285;
t509 = t182 / 0.2e1;
t185 = Icges(6,4) * t242 + Icges(6,2) * t324;
t508 = -t185 / 0.2e1;
t186 = Icges(6,1) * t242 + Icges(6,4) * t324;
t507 = t186 / 0.2e1;
t506 = -t225 / 0.2e1;
t505 = -t227 / 0.2e1;
t504 = t228 / 0.2e1;
t503 = -t324 / 0.2e1;
t502 = -t242 / 0.2e1;
t368 = -Icges(4,2) * t285 + t460;
t208 = Icges(4,6) * t298 + t300 * t368;
t461 = Icges(4,4) * t285;
t372 = Icges(4,1) * t286 - t461;
t212 = Icges(4,5) * t298 + t300 * t372;
t351 = t208 * t285 - t212 * t286;
t333 = t351 * t298;
t207 = -Icges(4,6) * t300 + t298 * t368;
t211 = -Icges(4,5) * t300 + t298 * t372;
t352 = t207 * t285 - t211 * t286;
t334 = t352 * t300;
t364 = Icges(5,3) * t285 + t456;
t202 = Icges(5,6) * t298 + t300 * t364;
t457 = Icges(5,5) * t285;
t370 = Icges(5,1) * t286 + t457;
t210 = Icges(5,4) * t298 + t300 * t370;
t353 = t202 * t285 + t210 * t286;
t335 = t353 * t298;
t201 = -Icges(5,6) * t300 + t298 * t364;
t209 = -Icges(5,4) * t300 + t298 * t370;
t354 = t201 * t285 + t209 * t286;
t336 = t354 * t300;
t377 = -rSges(7,1) * t194 - rSges(7,2) * t193;
t118 = -rSges(7,3) * t225 - t377;
t473 = t226 * pkin(5);
t441 = -pkin(9) * t225 + t118 - t473;
t497 = t441 * t300;
t290 = t298 * rSges(5,2);
t443 = t286 * t300;
t445 = t285 * t300;
t221 = rSges(5,1) * t443 + rSges(5,3) * t445 + t290;
t230 = pkin(3) * t443 + qJ(4) * t445;
t496 = -t221 - t230;
t495 = t242 * qJD(1);
t289 = t298 * rSges(4,3);
t494 = -rSges(4,2) * t445 + t289;
t365 = Icges(4,5) * t286 - Icges(4,6) * t285;
t203 = -Icges(4,3) * t300 + t298 * t365;
t366 = Icges(5,4) * t286 + Icges(5,6) * t285;
t205 = -Icges(5,2) * t300 + t298 * t366;
t304 = (-qJD(3) + qJD(5)) * t242;
t423 = qJD(1) * t298;
t140 = t300 * t304 + t324 * t423;
t422 = qJD(1) * t300;
t100 = -qJD(6) * t196 - t140 * t297 - t299 * t422;
t101 = qJD(6) * t195 + t140 * t299 - t297 * t422;
t310 = t324 * qJD(3);
t139 = -t298 * t495 - t300 * t310 - t392 * t445 - t393 * t443;
t53 = Icges(7,5) * t101 + Icges(7,6) * t100 - Icges(7,3) * t139;
t54 = Icges(7,4) * t101 + Icges(7,2) * t100 - Icges(7,6) * t139;
t55 = Icges(7,1) * t101 + Icges(7,4) * t100 - Icges(7,5) * t139;
t12 = t100 * t115 + t101 * t117 - t113 * t139 + t195 * t54 + t196 * t55 - t227 * t53;
t444 = t286 * t298;
t141 = -t285 * t298 * t394 - qJD(5) * t226 - t300 * t495 - t395 * t444;
t147 = -Icges(7,3) * t324 + (Icges(7,5) * t299 - Icges(7,6) * t297) * t242;
t148 = -Icges(7,6) * t324 + (Icges(7,4) * t299 - Icges(7,2) * t297) * t242;
t149 = -Icges(7,5) * t324 + (Icges(7,1) * t299 - Icges(7,4) * t297) * t242;
t415 = qJD(6) * t242;
t398 = t297 * t415;
t183 = qJD(5) * t324 - t310;
t446 = t183 * t299;
t339 = -t398 + t446;
t397 = t299 * t415;
t447 = t183 * t297;
t340 = -t397 - t447;
t79 = Icges(7,5) * t339 + Icges(7,6) * t340 + Icges(7,3) * t182;
t80 = Icges(7,4) * t339 + Icges(7,2) * t340 + Icges(7,6) * t182;
t81 = Icges(7,1) * t339 + Icges(7,4) * t340 + Icges(7,5) * t182;
t16 = t100 * t148 + t101 * t149 - t139 * t147 + t195 * t80 + t196 * t81 - t227 * t79;
t142 = -qJD(1) * t228 + t298 * t304;
t102 = -qJD(6) * t194 - t142 * t297 - t299 * t423;
t103 = qJD(6) * t193 + t142 * t299 - t297 * t423;
t320 = Icges(7,5) * t103 + Icges(7,6) * t102 + Icges(7,3) * t141;
t321 = Icges(7,4) * t103 + Icges(7,2) * t102 + Icges(7,6) * t141;
t322 = Icges(7,1) * t103 + Icges(7,4) * t102 + Icges(7,5) * t141;
t302 = t100 * t114 + t101 * t116 - t139 * t112 + t195 * t321 + t196 * t322 - t227 * t320;
t61 = -t147 * t227 + t148 * t195 + t149 * t196;
t1 = -t12 * t227 - t49 * t139 + t48 * t141 - t16 * t324 + t61 * t182 - t225 * t302;
t60 = -t147 * t225 + t148 * t193 + t149 * t194;
t19 = -t225 * t46 - t227 * t47 - t324 * t60;
t13 = t102 * t115 + t103 * t117 + t113 * t141 + t193 * t54 + t194 * t55 - t225 * t53;
t17 = t102 * t148 + t103 * t149 + t141 * t147 + t193 * t80 + t194 * t81 - t225 * t79;
t301 = t102 * t114 + t103 * t116 + t112 * t141 + t193 * t321 + t194 * t322 - t225 * t320;
t2 = -t13 * t227 - t139 * t47 + t141 * t46 - t17 * t324 + t182 * t60 - t225 * t301;
t20 = -t225 * t48 - t227 * t49 - t324 * t61;
t3 = t12 * t298 - t302 * t300 + (t298 * t48 + t300 * t49) * qJD(1);
t4 = t13 * t298 - t301 * t300 + (t298 * t46 + t300 * t47) * qJD(1);
t361 = -t114 * t297 + t116 * t299;
t50 = -t112 * t324 + t242 * t361;
t360 = -t115 * t297 + t117 * t299;
t51 = -t113 * t324 + t242 * t360;
t8 = t182 * t112 - t324 * t320 + t361 * t183 + (t299 * t322 - t297 * t321 + (-t114 * t299 - t116 * t297) * qJD(6)) * t242;
t9 = t113 * t182 - t324 * t53 + t360 * t183 + (-t297 * t54 + t299 * t55 + (-t115 * t299 - t117 * t297) * qJD(6)) * t242;
t493 = (t19 * t515 + t20 * t481) * qJD(1) + t139 * t27 / 0.2e1 - t141 * t26 / 0.2e1 - (t298 * t51 - t300 * t50) * t509 - t4 * t506 - t3 * t505 - (t9 * t298 - t8 * t300 + (t298 * t50 + t300 * t51) * qJD(1)) * t503 - t1 * t501 - t2 * t481;
t295 = cos(pkin(10));
t346 = rSges(3,1) * t295 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t464 = rSges(3,3) + qJ(2);
t200 = t298 * t464 + t300 * t346;
t119 = rSges(7,1) * t196 + rSges(7,2) * t195 - rSges(7,3) * t227;
t440 = -pkin(5) * t228 - pkin(9) * t227 + t119;
t58 = -t298 * t441 - t300 * t440;
t56 = rSges(7,1) * t101 + rSges(7,2) * t100 - rSges(7,3) * t139;
t470 = pkin(5) * t140 - t139 * pkin(9) + t56;
t474 = t142 * pkin(5);
t378 = -t103 * rSges(7,1) - t102 * rSges(7,2);
t57 = t141 * rSges(7,3) - t378;
t492 = -(t141 * pkin(9) + t474 + t57) * t298 - t470 * t300;
t491 = 2 * m(4);
t490 = 2 * m(5);
t489 = 2 * m(6);
t488 = 2 * m(7);
t487 = m(5) / 0.2e1;
t486 = m(6) / 0.2e1;
t485 = m(7) / 0.2e1;
t484 = -pkin(3) - pkin(4);
t311 = Icges(6,5) * t142 - Icges(6,6) * t141 - Icges(6,3) * t423;
t312 = -Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t423;
t313 = Icges(6,1) * t142 - Icges(6,4) * t141 - Icges(6,5) * t423;
t85 = Icges(6,5) * t140 + Icges(6,6) * t139 - Icges(6,3) * t422;
t86 = Icges(6,4) * t140 + Icges(6,2) * t139 - Icges(6,6) * t422;
t87 = Icges(6,1) * t140 + Icges(6,4) * t139 - Icges(6,5) * t422;
t6 = (t139 * t158 + t140 * t160 + t227 * t86 - t228 * t87 - t298 * t85) * t298 - (t139 * t157 + t140 * t159 - t155 * t422 - t227 * t312 - t228 * t313 - t298 * t311) * t300 + (t69 * t300 + (t318 - t449) * t298) * qJD(1);
t483 = t6 + t3;
t7 = (-t141 * t158 + t142 * t160 + t225 * t86 - t226 * t87 + t300 * t85) * t298 - (-t141 * t157 + t142 * t159 - t155 * t423 - t225 * t312 - t226 * t313 + t300 * t311) * t300 + (t68 * t300 + (t319 - t450) * t298) * qJD(1);
t482 = -t7 - t4;
t480 = -rSges(5,1) - pkin(3);
t479 = -rSges(7,3) - pkin(9);
t255 = rSges(4,1) * t285 + rSges(4,2) * t286;
t476 = m(4) * t255;
t475 = m(6) * t187;
t296 = -pkin(7) - qJ(2);
t472 = -pkin(8) - t296;
t325 = t242 * t299 * t81 + t147 * t182 - t148 * t447 + t149 * t446 - t324 * t79;
t465 = t297 * t80;
t64 = -t147 * t324 + (-t148 * t297 + t149 * t299) * t242;
t471 = -((-t465 + (-t148 * t299 - t149 * t297) * qJD(6)) * t242 + t325) * t324 + t64 * t182;
t468 = rSges(4,1) * t286;
t467 = rSges(5,1) * t285;
t466 = rSges(6,3) * t298;
t463 = -rSges(5,3) - qJ(4);
t82 = rSges(7,1) * t339 + rSges(7,2) * t340 + rSges(7,3) * t182;
t462 = pkin(5) * t183 + pkin(9) * t182 + t82;
t451 = qJ(4) * t285;
t379 = t226 * rSges(6,1) - t225 * rSges(6,2);
t161 = rSges(6,3) * t300 - t379;
t448 = t161 * t300;
t442 = t296 * t300;
t439 = rSges(6,1) * t140 + rSges(6,2) * t139;
t150 = -rSges(7,3) * t324 + (rSges(7,1) * t299 - rSges(7,2) * t297) * t242;
t438 = pkin(5) * t242 - pkin(9) * t324 + t150;
t437 = -rSges(6,1) * t228 + rSges(6,2) * t227;
t376 = pkin(3) * t286 + t451;
t218 = qJD(3) * t376 - qJD(4) * t286;
t381 = rSges(5,1) * t286 + rSges(5,3) * t285;
t436 = -qJD(3) * t381 - t218;
t229 = t376 * t298;
t435 = t229 * t298 + t230 * t300;
t278 = pkin(4) * t443;
t252 = -pkin(8) * t298 + t278;
t434 = -t230 - t252;
t253 = pkin(3) * t285 - qJ(4) * t286;
t231 = t253 * t423;
t403 = t285 * t423;
t433 = pkin(4) * t403 + t231;
t254 = -rSges(5,3) * t286 + t467;
t432 = -t253 - t254;
t417 = qJD(3) * t300;
t399 = t286 * t417;
t416 = qJD(4) * t285;
t431 = qJ(4) * t399 + t300 * t416;
t430 = rSges(5,2) * t422 + rSges(5,3) * t399;
t418 = qJD(3) * t298;
t401 = t285 * t418;
t429 = pkin(4) * t401 + pkin(8) * t423;
t288 = qJD(2) * t300;
t428 = t296 * t423 + t288;
t426 = qJD(1) * t187;
t204 = Icges(4,3) * t298 + t300 * t365;
t425 = qJD(1) * t204;
t206 = Icges(5,2) * t298 + t300 * t366;
t424 = qJD(1) * t206;
t420 = qJD(3) * t285;
t419 = qJD(3) * t286;
t414 = -rSges(6,3) + t472;
t413 = t16 / 0.2e1 + t9 / 0.2e1;
t412 = t8 / 0.2e1 + t17 / 0.2e1;
t410 = t50 / 0.2e1 + t60 / 0.2e1;
t409 = -t61 / 0.2e1 - t51 / 0.2e1;
t268 = pkin(3) * t401;
t400 = t285 * t417;
t323 = -t286 * t423 - t400;
t402 = t286 * t422;
t408 = t298 * (pkin(3) * t402 + t298 * t416 - t268 + (t285 * t422 + t286 * t418) * qJ(4)) + t300 * (pkin(3) * t323 - qJ(4) * t403 + t431) + t229 * t422;
t287 = qJD(2) * t298;
t407 = t287 + t431;
t406 = t268 + t428;
t396 = (t365 + t366) * qJD(3) / 0.2e1;
t391 = -pkin(4) * t285 - t253;
t190 = t432 * t300;
t281 = pkin(2) * t295 + pkin(1);
t274 = t300 * t281;
t390 = -t296 * t298 + t274;
t389 = qJD(1) * t438;
t251 = pkin(4) * t444 + pkin(8) * t300;
t388 = t251 * t298 + t252 * t300 + t435;
t387 = t274 + t278 + t230;
t384 = -t187 + t391;
t383 = -pkin(4) * t419 - t218;
t382 = -rSges(4,2) * t285 + t468;
t380 = -t142 * rSges(6,1) + t141 * rSges(6,2);
t70 = t118 * t324 - t150 * t225;
t71 = -t119 * t324 + t150 * t227;
t375 = t298 * t71 + t300 * t70;
t326 = t286 * t484 - t281 - t451;
t317 = t326 * t298;
t307 = t300 * t472 + t317;
t77 = -t225 * t479 + t307 + t377 + t473;
t78 = t298 * t472 + t387 + t440;
t374 = t298 * t78 + t300 * t77;
t373 = -t298 * (-rSges(6,3) * t423 - t380) - t300 * (-rSges(6,3) * t422 + t439);
t367 = Icges(4,2) * t286 + t461;
t363 = -Icges(5,3) * t286 + t457;
t305 = t300 * t414 + t317;
t109 = t305 + t379;
t110 = t298 * t414 + t387 + t437;
t362 = t109 * t300 + t110 * t298;
t315 = t285 * t463 + t286 * t480 - t281;
t308 = t315 * t298;
t151 = (rSges(5,2) - t296) * t300 + t308;
t152 = t390 - t496;
t359 = t151 * t300 + t152 * t298;
t162 = t437 - t466;
t104 = -t161 * t298 - t162 * t300;
t350 = t391 - t438;
t222 = rSges(4,1) * t443 + t494;
t133 = rSges(6,1) * t183 - rSges(6,2) * t182;
t349 = -t133 + t383;
t130 = Icges(6,5) * t183 - Icges(6,6) * t182;
t131 = Icges(6,4) * t183 - Icges(6,2) * t182;
t132 = Icges(6,1) * t183 - Icges(6,4) * t182;
t184 = Icges(6,5) * t242 + Icges(6,6) * t324;
t348 = t130 * t501 + t131 * t505 + t132 * t504 + t139 * t508 - t140 * t186 / 0.2e1 + t184 * t422 / 0.2e1 + t158 * t509 - t160 * t183 / 0.2e1 + t86 * t503 + t87 * t502 - t413;
t347 = t130 * t514 + t225 * t131 / 0.2e1 - t226 * t132 / 0.2e1 + t141 * t508 + t142 * t507 - t184 * t423 / 0.2e1 - t182 * t157 / 0.2e1 + t183 * t159 / 0.2e1 + t312 * t503 + t242 * t313 / 0.2e1 + t412;
t345 = t157 * t503 + t159 * t502 + t184 * t481 + t185 * t506 + t226 * t507 - t410;
t344 = t158 * t503 + t160 * t502 + t184 * t501 + t185 * t505 + t186 * t504 + t409;
t146 = t384 * t300;
t343 = -t281 - t382;
t341 = t298 * (pkin(4) * t402 - t429) + t300 * (pkin(4) * t323 - pkin(8) * t422) + t251 * t422 + t408;
t338 = t383 - t462;
t337 = qJD(3) * t255;
t330 = qJD(3) * t367;
t329 = qJD(3) * (-Icges(5,4) * t285 + Icges(5,6) * t286);
t328 = qJD(3) * (-Icges(4,5) * t285 - Icges(4,6) * t286);
t327 = qJD(3) * t363;
t95 = t350 * t300;
t316 = t326 * t300;
t314 = t400 * t484 + t407;
t309 = rSges(4,2) * t403 + rSges(4,3) * t422 - t300 * t337;
t199 = -t298 * t346 + t300 * t464;
t306 = (-qJ(4) * t419 - t416) * t298 + t406 + t429;
t244 = t382 * qJD(3);
t220 = -rSges(4,3) * t300 + t298 * t382;
t219 = -rSges(5,2) * t300 + t298 * t381;
t189 = t432 * t298;
t181 = t222 + t390;
t180 = (rSges(4,3) - t296) * t300 + t343 * t298;
t179 = -qJD(1) * t200 + t288;
t178 = qJD(1) * t199 + t287;
t168 = t298 * t329 + t424;
t167 = -qJD(1) * t205 + t300 * t329;
t166 = t298 * t328 + t425;
t165 = -qJD(1) * t203 + t300 * t328;
t145 = t384 * t298;
t144 = t255 * t418 + (t300 * t343 - t289) * qJD(1) + t428;
t143 = t287 + (-t442 + (-t281 - t468) * t298) * qJD(1) + t309;
t129 = qJD(1) * t190 + t298 * t436;
t128 = t254 * t423 + t300 * t436 + t231;
t127 = t204 * t298 - t300 * t351;
t126 = t203 * t298 - t334;
t125 = t206 * t298 + t300 * t353;
t124 = t205 * t298 + t336;
t123 = -t204 * t300 - t333;
t122 = -t203 * t300 - t298 * t352;
t121 = -t206 * t300 + t335;
t120 = -t205 * t300 + t298 * t354;
t111 = t219 * t298 + t221 * t300 + t435;
t106 = t438 * t300;
t105 = t438 * t298;
t97 = (-t416 + (t286 * t463 + t467) * qJD(3)) * t298 + (t300 * t315 - t290) * qJD(1) + t406;
t96 = t480 * t400 + (t308 - t442) * qJD(1) + t407 + t430;
t94 = t350 * t298;
t75 = -t104 + t388;
t73 = qJD(1) * t146 + t298 * t349;
t72 = t187 * t423 + t300 * t349 + t433;
t65 = -t118 * t227 + t119 * t225;
t63 = (t316 + t466) * qJD(1) + t306 + t380;
t62 = qJD(1) * t305 + t314 + t439;
t52 = t300 * t430 + (-t254 * t292 - t293 * t467) * qJD(3) + (t300 * t219 + (t290 + t496) * t298) * qJD(1) + t408;
t43 = t298 * t462 + t300 * t389;
t42 = t300 * t462 - t423 * t438;
t41 = t388 - t58;
t38 = (t162 * t298 - t448) * qJD(1) + t373;
t37 = qJD(1) * t95 + t298 * t338;
t36 = t298 * t389 + t300 * t338 + t433;
t31 = qJD(1) * t316 + t141 * t479 + t306 + t378 - t474;
t30 = qJD(1) * t307 + t314 + t470;
t28 = (t448 + (-t162 + t434) * t298) * qJD(1) + t341 - t373;
t25 = -t118 * t182 + t141 * t150 - t225 * t82 + t324 * t57;
t24 = t119 * t182 + t139 * t150 + t227 * t82 - t324 * t56;
t21 = -t118 * t139 - t119 * t141 + t225 * t56 - t227 * t57;
t15 = (t298 * t440 - t497) * qJD(1) + t492;
t14 = (t497 + (t434 - t440) * t298) * qJD(1) + t341 - t492;
t5 = [t325 + t324 * t131 - t182 * t185 + t183 * t186 - t148 * t397 - t149 * t398 + 0.2e1 * m(3) * (t178 * t200 + t179 * t199) + (t143 * t181 + t144 * t180) * t491 + (t151 * t97 + t152 * t96) * t490 + (t109 * t63 + t110 * t62) * t489 + (t30 * t78 + t31 * t77) * t488 + (-t465 + t132) * t242 + (t363 + t370 - t367 + t372) * t420 + (-t364 + t368 + t518) * t419; m(7) * (qJD(1) * t374 + t298 * t31 - t30 * t300) + m(6) * (qJD(1) * t362 + t298 * t63 - t300 * t62) + m(5) * (qJD(1) * t359 + t298 * t97 - t300 * t96) + m(4) * (-t143 * t300 + t298 * t144 + (t180 * t300 + t181 * t298) * qJD(1)) + m(3) * (-t178 * t300 + t298 * t179 + (t199 * t300 + t200 * t298) * qJD(1)); 0; m(5) * (t128 * t151 + t129 * t152 + t189 * t96 + t190 * t97) + m(6) * (t109 * t72 + t110 * t73 + t145 * t62 + t146 * t63) + m(7) * (t30 * t94 + t31 * t95 + t36 * t77 + t37 * t78) + (m(4) * (-t144 * t255 - t180 * t244) + t396 * t300 + (t202 * t512 + t208 * t513 + t327 * t515 + t330 * t501) * t286 + ((t210 + t212) * t513 + t517 * t501) * t285 - t347) * t300 + (m(4) * (-t143 * t255 - t181 * t244) + t396 * t298 + (t201 * t512 + t207 * t513 + t327 * t514 + t330 * t481) * t286 + ((t209 + t211) * t513 + t517 * t481) * t285 - t348) * t298 + (-t336 / 0.2e1 + t334 / 0.2e1 - t333 / 0.2e1 + t335 / 0.2e1) * qJD(3) + ((-t181 * t476 + (-t202 / 0.2e1 + t208 / 0.2e1) * t286 + (t210 / 0.2e1 + t212 / 0.2e1) * t285 - t344) * t300 + (t180 * t476 + (-t201 / 0.2e1 + t207 / 0.2e1) * t286 + (t209 / 0.2e1 + t211 / 0.2e1) * t285 - t345) * t298) * qJD(1); m(5) * (t128 * t298 - t129 * t300 + (t189 * t298 + t190 * t300) * qJD(1)) + m(6) * (t72 * t298 - t300 * t73 + (t145 * t298 + t146 * t300) * qJD(1)) + m(7) * (t36 * t298 - t300 * t37 + (t298 * t94 + t300 * t95) * qJD(1)); (t14 * t41 + t36 * t95 + t37 * t94) * t488 - t300 * t4 + t298 * t3 + (t145 * t73 + t146 * t72 + t28 * t75) * t489 + t298 * t6 - t300 * t7 + (t111 * t52 + t128 * t190 + t129 * t189) * t490 + ((t220 * t298 + t222 * t300) * ((qJD(1) * t220 + t309) * t300 + (-t298 * t337 + (-t222 + t494) * qJD(1)) * t298) + t427 * t255 * t244) * t491 - t300 * ((t300 * t168 + (t121 - t336) * qJD(1)) * t300 + (t120 * qJD(1) + (t202 * t419 - t210 * t420 + t424) * t298 + (-t167 + (-t201 * t286 + t209 * t285) * qJD(3) + t353 * qJD(1)) * t300) * t298) + t298 * ((t298 * t165 + (t126 + t333) * qJD(1)) * t298 + (t127 * qJD(1) + (t207 * t419 + t211 * t420) * t300 + (-t166 + (-t208 * t286 - t212 * t285) * qJD(3) + (t204 - t352) * qJD(1)) * t298) * t300) - t300 * ((t300 * t166 + (t123 + t334) * qJD(1)) * t300 + (t122 * qJD(1) + (-t208 * t419 - t212 * t420 + t425) * t298 + (-t165 + (t207 * t286 + t211 * t285) * qJD(3) - t351 * qJD(1)) * t300) * t298) + t298 * ((t298 * t167 + (t124 - t335) * qJD(1)) * t298 + (t125 * qJD(1) + (-t201 * t419 + t209 * t420) * t300 + (-t168 + (t202 * t286 - t210 * t285) * qJD(3) + (t206 + t354) * qJD(1)) * t298) * t300) + ((-t120 - t122) * t300 + (t121 + t123) * t298 + t499) * t423 + ((-t124 - t126) * t300 + (t125 + t127) * t298 + t498) * t422; 0.2e1 * (t359 * t487 + t362 * t486 + t374 * t485) * t419 + 0.2e1 * ((t298 * t30 + t300 * t31 + t422 * t78 - t423 * t77) * t485 + (-t109 * t423 + t110 * t422 + t298 * t62 + t300 * t63) * t486 + (-t151 * t423 + t152 * t422 + t298 * t96 + t300 * t97) * t487) * t285; 0; 0.2e1 * ((t417 * t95 + t418 * t94 - t14) * t485 + (t145 * t418 + t146 * t417 - t28) * t486 + (t189 * t418 + t190 * t417 - t52) * t487) * t286 + 0.2e1 * ((qJD(3) * t41 + t298 * t37 + t300 * t36 + t422 * t94 - t423 * t95) * t485 + (qJD(3) * t75 + t145 * t422 - t146 * t423 + t298 * t73 + t300 * t72) * t486 + (qJD(3) * t111 + t128 * t300 + t129 * t298 + t189 * t422 - t190 * t423) * t487) * t285; 0.4e1 * (t487 + t486 + t485) * (-0.1e1 + t427) * t285 * t419; m(7) * (t105 * t30 + t106 * t31 + t42 * t77 + t43 * t78) + (m(6) * (t109 * t133 + t187 * t63) + (t110 * t475 + t344) * qJD(1) + t347) * t300 + (m(6) * (t110 * t133 + t187 * t62) + (-t109 * t475 + t345) * qJD(1) + t348) * t298; m(7) * (t42 * t298 - t300 * t43 + (t105 * t298 + t106 * t300) * qJD(1)); m(7) * (t105 * t37 + t106 * t36 + t14 * t58 + t15 * t41 + t42 * t95 + t43 * t94) + m(6) * (t104 * t28 + t38 * t75) + (m(6) * (t133 * t146 + t145 * t426 + t187 * t72) - t482) * t300 + (m(6) * (t133 * t145 - t146 * t426 + t187 * t73) - t483) * t298 - t516; 0.2e1 * ((qJD(3) * t511 - t38) * t486 + (t105 * t418 + t106 * t417 - t15) * t485) * t286 + 0.2e1 * ((qJD(3) * t104 + t133 * t427) * t486 + (qJD(3) * t58 + t105 * t422 - t106 * t423 + t298 * t43 + t300 * t42) * t485) * t285; t482 * t300 + t483 * t298 + (t104 * t38 + t133 * t511) * t489 + (t105 * t43 + t106 * t42 + t15 * t58) * t488 + t516; m(7) * (t24 * t78 + t25 * t77 + t30 * t71 + t31 * t70) - t413 * t227 - t412 * t225 + t410 * t141 + t409 * t139 + t471; m(7) * (qJD(1) * t375 - t24 * t300 + t25 * t298); m(7) * (t14 * t65 + t21 * t41 + t24 * t94 + t25 * t95 + t36 * t70 + t37 * t71) - t493; m(7) * ((qJD(3) * t375 - t21) * t286 + (qJD(3) * t65 + t24 * t298 + t25 * t300 + (-t298 * t70 + t300 * t71) * qJD(1)) * t285); m(7) * (t105 * t24 + t106 * t25 + t15 * t65 + t21 * t58 + t42 * t70 + t43 * t71) + t493; -t139 * t20 - t227 * t1 + t141 * t19 - t225 * t2 + t182 * (-t225 * t50 - t227 * t51 - t324 * t64) - t324 * (-t139 * t51 + t141 * t50 - t225 * t8 - t227 * t9 + t471) + (t21 * t65 + t24 * t71 + t25 * t70) * t488;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
