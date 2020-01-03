% Calculate time derivative of joint inertia matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR13_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR13_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:14:54
% DurationCPUTime: 11.26s
% Computational Cost: add. (22487->785), mult. (37136->1112), div. (0->0), fcn. (35880->8), ass. (0->393)
t302 = sin(qJ(1));
t305 = cos(qJ(1));
t301 = sin(qJ(3));
t416 = qJD(3) * t301;
t379 = -t416 / 0.2e1;
t304 = cos(qJ(3));
t418 = qJD(1) * t304;
t380 = -t418 / 0.2e1;
t507 = t302 * t380 + t305 * t379;
t378 = t416 / 0.2e1;
t506 = t302 * t378 + t305 * t380;
t413 = qJD(3) * t305;
t385 = t301 * t413;
t391 = t302 * t418;
t505 = t385 + t391;
t415 = qJD(3) * t302;
t386 = t301 * t415;
t417 = qJD(1) * t305;
t390 = t304 * t417;
t315 = t386 - t390;
t414 = qJD(3) * t304;
t389 = t302 * t414;
t504 = t301 * t417 + t389;
t300 = sin(qJ(4));
t303 = cos(qJ(4));
t466 = Icges(5,4) * t303;
t340 = -Icges(5,2) * t300 + t466;
t225 = Icges(5,6) * t301 + t304 * t340;
t467 = Icges(5,4) * t300;
t344 = Icges(5,1) * t303 - t467;
t228 = Icges(5,5) * t301 + t304 * t344;
t503 = -t225 * t300 + t228 * t303;
t299 = qJ(4) + qJ(5);
t289 = sin(t299);
t290 = cos(t299);
t464 = Icges(6,4) * t290;
t339 = -Icges(6,2) * t289 + t464;
t217 = Icges(6,6) * t301 + t304 * t339;
t465 = Icges(6,4) * t289;
t343 = Icges(6,1) * t290 - t465;
t218 = Icges(6,5) * t301 + t304 * t343;
t502 = -t217 * t289 + t218 * t290;
t468 = Icges(4,4) * t304;
t345 = Icges(4,1) * t301 + t468;
t229 = Icges(4,5) * t305 + t302 * t345;
t448 = t229 * t301;
t469 = Icges(4,4) * t301;
t341 = Icges(4,2) * t304 + t469;
t226 = Icges(4,6) * t305 + t302 * t341;
t452 = t226 * t304;
t327 = t448 + t452;
t501 = t305 * t327;
t359 = rSges(4,1) * t301 + rSges(4,2) * t304;
t319 = t305 * t359;
t306 = -pkin(8) - pkin(7);
t478 = -pkin(7) - t306;
t286 = pkin(4) * t303 + pkin(3);
t479 = pkin(3) - t286;
t215 = t478 * t301 - t479 * t304;
t353 = rSges(6,1) * t290 - rSges(6,2) * t289;
t219 = rSges(6,3) * t301 + t304 * t353;
t426 = t215 + t219;
t375 = t426 * t302;
t500 = t479 * t301;
t398 = t504 * rSges(4,1) + rSges(4,2) * t390;
t488 = -pkin(1) - pkin(6);
t410 = -rSges(4,3) + t488;
t423 = qJ(2) * t417 + qJD(2) * t302;
t142 = (-rSges(4,2) * t416 + t410 * qJD(1)) * t302 + t398 + t423;
t477 = rSges(4,2) * t301;
t265 = rSges(4,1) * t304 - t477;
t288 = qJD(2) * t305;
t143 = t288 + t265 * t413 + (t410 * t305 + (-qJ(2) - t359) * t302) * qJD(1);
t499 = -t142 * t305 + t143 * t302;
t420 = qJD(1) * t301;
t366 = qJD(4) + t420;
t388 = t304 * t413;
t498 = t366 * t302 - t388;
t296 = qJD(4) + qJD(5);
t372 = t296 + t420;
t497 = t372 * t302 - t388;
t444 = t286 * t301;
t473 = rSges(6,3) - t306;
t496 = -t473 * t304 + t444;
t337 = Icges(4,5) * t301 + Icges(4,6) * t304;
t495 = -Icges(4,3) * t302 + t305 * t337;
t494 = -Icges(4,6) * t302 + t305 * t341;
t493 = -Icges(4,5) * t302 + t305 * t345;
t492 = t372 * t305 + t389;
t491 = 2 * m(4);
t490 = 2 * m(5);
t489 = 2 * m(6);
t297 = t302 ^ 2;
t298 = t305 ^ 2;
t487 = -t301 / 0.2e1;
t486 = t301 / 0.2e1;
t485 = -t302 / 0.2e1;
t484 = t304 / 0.2e1;
t483 = t305 / 0.2e1;
t482 = rSges(3,2) - pkin(1);
t481 = m(4) * t265;
t480 = pkin(3) * t301;
t476 = rSges(5,3) * t304;
t475 = pkin(4) * qJD(4);
t474 = t302 * rSges(4,3);
t293 = t305 * rSges(4,3);
t440 = t300 * t305;
t283 = pkin(4) * t440;
t408 = t303 * t475;
t369 = t302 * t408;
t409 = t300 * t475;
t370 = t301 * t409;
t396 = pkin(3) * t388 + t505 * pkin(7);
t432 = t304 * t306;
t434 = t304 * t286;
t373 = t296 * t301 + qJD(1);
t324 = t373 * t305;
t147 = -t497 * t289 + t290 * t324;
t148 = t289 * t324 + t497 * t290;
t355 = t148 * rSges(6,1) + t147 * rSges(6,2);
t99 = -rSges(6,3) * t505 + t355;
t472 = t369 + (t370 + (t301 * t306 - t434) * qJD(3)) * t305 + (t283 + (t432 - t500) * t302) * qJD(1) + t396 + t99;
t335 = Icges(6,5) * t290 - Icges(6,6) * t289;
t216 = Icges(6,3) * t301 + t304 * t335;
t125 = t216 * t301 + t502 * t304;
t443 = t296 * t304;
t159 = (-Icges(6,2) * t290 - t465) * t443 + (Icges(6,6) * t304 - t301 * t339) * qJD(3);
t158 = (-Icges(6,5) * t289 - Icges(6,6) * t290) * t443 + (Icges(6,3) * t304 - t301 * t335) * qJD(3);
t160 = (-Icges(6,1) * t289 - t464) * t443 + (Icges(6,5) * t304 - t301 * t343) * qJD(3);
t311 = t304 * t290 * t160 + t301 * t158 + t216 * t414 - t502 * t416;
t405 = t296 * t290 * t217;
t471 = t125 * t414 + ((-t405 + (-t218 * t296 - t159) * t289) * t304 + t311) * t301;
t336 = Icges(5,5) * t303 - Icges(5,6) * t300;
t222 = Icges(5,3) * t301 + t304 * t336;
t133 = t222 * t301 + t503 * t304;
t411 = qJD(4) * t304;
t189 = (-Icges(5,5) * t300 - Icges(5,6) * t303) * t411 + (Icges(5,3) * t304 - t301 * t336) * qJD(3);
t195 = (-Icges(5,1) * t300 - t466) * t411 + (Icges(5,5) * t304 - t301 * t344) * qJD(3);
t310 = t304 * t303 * t195 + t301 * t189 + t222 * t414 - t503 * t416;
t192 = (-Icges(5,2) * t303 - t467) * t411 + (Icges(5,6) * t304 - t301 * t340) * qJD(3);
t442 = t300 * t192;
t454 = t225 * t303;
t470 = t133 * t414 + ((-t442 + (-t228 * t300 - t454) * qJD(4)) * t304 + t310) * t301;
t437 = t302 * t303;
t438 = t301 * t305;
t245 = t300 * t438 + t437;
t435 = t303 * t305;
t441 = t300 * t302;
t246 = -t301 * t435 + t441;
t357 = -rSges(5,1) * t246 - rSges(5,2) * t245;
t433 = t304 * t305;
t186 = rSges(5,3) * t433 - t357;
t458 = t186 * t305;
t453 = t226 * t301;
t451 = t494 * t301;
t450 = t494 * t304;
t447 = t229 * t304;
t446 = t493 * t301;
t445 = t493 * t304;
t439 = t301 * t302;
t436 = t302 * t304;
t325 = t302 * t373;
t149 = -t492 * t289 - t290 * t325;
t150 = -t289 * t325 + t492 * t290;
t100 = t150 * rSges(6,1) + t149 * rSges(6,2) + t315 * rSges(6,3);
t397 = t504 * pkin(3) + pkin(7) * t386;
t313 = pkin(7) * t390 - t397;
t367 = qJD(4) * t301 + qJD(1);
t323 = t367 * t300;
t365 = t504 * t286 + t305 * t408 + t306 * t390;
t412 = qJD(3) * t306;
t118 = (-pkin(4) * t323 - t301 * t412) * t302 + t313 + t365;
t431 = -t100 - t118;
t162 = (-rSges(6,1) * t289 - rSges(6,2) * t290) * t443 + (rSges(6,3) * t304 - t301 * t353) * qJD(3);
t381 = t478 * t304;
t387 = t300 * t411;
t199 = -pkin(4) * t387 + (t381 + t500) * qJD(3);
t430 = t162 + t199;
t234 = -t289 * t439 + t290 * t305;
t235 = t289 * t305 + t290 * t439;
t169 = t235 * rSges(6,1) + t234 * rSges(6,2) - rSges(6,3) * t436;
t282 = pkin(3) * t439;
t248 = -pkin(7) * t436 + t282;
t399 = t286 * t439 + t302 * t432 + t283;
t187 = -t248 + t399;
t429 = -t169 - t187;
t236 = t289 * t438 + t290 * t302;
t237 = t289 * t302 - t290 * t438;
t354 = -t237 * rSges(6,1) - t236 * rSges(6,2);
t170 = rSges(6,3) * t433 - t354;
t284 = pkin(3) * t438;
t188 = pkin(4) * t441 + t284 + (t381 - t444) * t305;
t428 = t170 + t188;
t243 = -t300 * t439 + t435;
t244 = t301 * t437 + t440;
t425 = t244 * rSges(5,1) + t243 * rSges(5,2);
t185 = -rSges(5,3) * t436 + t425;
t427 = -t185 - t248;
t128 = t301 * t169 + t219 * t436;
t259 = (pkin(7) * t304 - t480) * qJD(3);
t266 = t304 * pkin(3) + t301 * pkin(7);
t424 = t302 * t259 + t266 * t417;
t422 = t305 * pkin(1) + t302 * qJ(2);
t223 = Icges(4,3) * t305 + t302 * t337;
t421 = qJD(1) * t223;
t419 = qJD(1) * t302;
t176 = -t367 * t437 + (-t366 * t305 - t389) * t300;
t177 = t366 * t435 + (t303 * t414 - t323) * t302;
t109 = Icges(5,5) * t177 + Icges(5,6) * t176 + Icges(5,3) * t315;
t111 = Icges(5,4) * t177 + Icges(5,2) * t176 + Icges(5,6) * t315;
t113 = Icges(5,1) * t177 + Icges(5,4) * t176 + Icges(5,5) * t315;
t179 = Icges(5,5) * t244 + Icges(5,6) * t243 - Icges(5,3) * t436;
t181 = Icges(5,4) * t244 + Icges(5,2) * t243 - Icges(5,6) * t436;
t183 = Icges(5,1) * t244 + Icges(5,4) * t243 - Icges(5,5) * t436;
t332 = t181 * t300 - t183 * t303;
t30 = (qJD(3) * t332 + t109) * t301 + (qJD(3) * t179 - t111 * t300 + t113 * t303 + (-t181 * t303 - t183 * t300) * qJD(4)) * t304;
t57 = t176 * t225 + t177 * t228 - t189 * t436 + t192 * t243 + t195 * t244 + t222 * t315;
t407 = -t30 / 0.2e1 - t57 / 0.2e1;
t322 = t367 * t305;
t174 = -t498 * t300 + t303 * t322;
t175 = t300 * t322 + t498 * t303;
t108 = Icges(5,5) * t175 + Icges(5,6) * t174 - Icges(5,3) * t505;
t110 = Icges(5,4) * t175 + Icges(5,2) * t174 - Icges(5,6) * t505;
t112 = Icges(5,1) * t175 + Icges(5,4) * t174 - Icges(5,5) * t505;
t180 = Icges(5,5) * t246 + Icges(5,6) * t245 + Icges(5,3) * t433;
t182 = Icges(5,4) * t246 + Icges(5,2) * t245 + Icges(5,6) * t433;
t184 = Icges(5,1) * t246 + Icges(5,4) * t245 + Icges(5,5) * t433;
t331 = t182 * t300 - t184 * t303;
t31 = (qJD(3) * t331 + t108) * t301 + (qJD(3) * t180 - t110 * t300 + t112 * t303 + (-t182 * t303 - t184 * t300) * qJD(4)) * t304;
t56 = t174 * t225 + t175 * t228 + t189 * t433 + t192 * t245 + t195 * t246 - t222 * t505;
t406 = t31 / 0.2e1 + t56 / 0.2e1;
t119 = -t222 * t436 + t225 * t243 + t228 * t244;
t91 = t179 * t301 - t304 * t332;
t404 = t91 / 0.2e1 + t119 / 0.2e1;
t120 = t222 * t433 + t225 * t245 + t228 * t246;
t92 = t180 * t301 - t304 * t331;
t403 = -t92 / 0.2e1 - t120 / 0.2e1;
t402 = t505 * t169 + t170 * t386;
t401 = -t248 + t429;
t400 = t177 * rSges(5,1) + t176 * rSges(5,2) + rSges(5,3) * t386;
t231 = rSges(4,1) * t439 + rSges(4,2) * t436 + t293;
t395 = t305 * pkin(6) + t422;
t394 = (-rSges(5,3) - pkin(7)) * t304;
t356 = rSges(5,1) * t303 - rSges(5,2) * t300;
t232 = rSges(5,3) * t301 + t356 * t304;
t393 = t232 * t419;
t384 = -t436 / 0.2e1;
t383 = t433 / 0.2e1;
t163 = Icges(6,5) * t235 + Icges(6,6) * t234 - Icges(6,3) * t436;
t165 = Icges(6,4) * t235 + Icges(6,2) * t234 - Icges(6,6) * t436;
t167 = Icges(6,1) * t235 + Icges(6,4) * t234 - Icges(6,5) * t436;
t334 = t165 * t289 - t167 * t290;
t94 = Icges(6,5) * t150 + Icges(6,6) * t149 + Icges(6,3) * t315;
t96 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t315;
t98 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t315;
t24 = (qJD(3) * t334 + t94) * t301 + (qJD(3) * t163 + (-t165 * t296 + t98) * t290 + (-t167 * t296 - t96) * t289) * t304;
t164 = Icges(6,5) * t237 + Icges(6,6) * t236 + Icges(6,3) * t433;
t166 = Icges(6,4) * t237 + Icges(6,2) * t236 + Icges(6,6) * t433;
t168 = Icges(6,1) * t237 + Icges(6,4) * t236 + Icges(6,5) * t433;
t333 = t166 * t289 - t168 * t290;
t93 = Icges(6,5) * t148 + Icges(6,6) * t147 - Icges(6,3) * t505;
t95 = Icges(6,4) * t148 + Icges(6,2) * t147 - Icges(6,6) * t505;
t97 = Icges(6,1) * t148 + Icges(6,4) * t147 - Icges(6,5) * t505;
t25 = (qJD(3) * t333 + t93) * t301 + (qJD(3) * t164 + (-t166 * t296 + t97) * t290 + (-t168 * t296 - t95) * t289) * t304;
t80 = t163 * t301 - t304 * t334;
t81 = t164 * t301 - t304 * t333;
t349 = t80 * t302 - t81 * t305;
t350 = t302 * t81 + t305 * t80;
t106 = -t216 * t436 + t217 * t234 + t218 * t235;
t70 = -t163 * t436 + t165 * t234 + t167 * t235;
t71 = -t164 * t436 + t166 * t234 + t168 * t235;
t352 = t302 * t70 - t305 * t71;
t37 = t106 * t301 - t304 * t352;
t107 = t216 * t433 + t217 * t236 + t218 * t237;
t18 = t147 * t165 + t148 * t167 - t163 * t505 + t236 * t96 + t237 * t98 + t433 * t94;
t19 = t147 * t166 + t148 * t168 - t164 * t505 + t236 * t95 + t237 * t97 + t433 * t93;
t72 = t163 * t433 + t165 * t236 + t167 * t237;
t73 = t164 * t433 + t166 * t236 + t168 * t237;
t351 = t302 * t72 - t305 * t73;
t46 = t147 * t217 + t148 * t218 + t158 * t433 + t159 * t236 + t160 * t237 - t216 * t505;
t54 = t302 * t73 + t305 * t72;
t4 = (qJD(3) * t351 + t46) * t301 + (-qJD(1) * t54 + qJD(3) * t107 - t18 * t302 + t19 * t305) * t304;
t382 = t4 * t433 + t37 * t386 + (t125 * t301 - t304 * t349) * t414 + t301 * (t349 * t416 + (-qJD(1) * t350 - t24 * t302 + t25 * t305) * t304 + t471);
t84 = -t179 * t436 + t181 * t243 + t183 * t244;
t85 = -t180 * t436 + t182 * t243 + t184 * t244;
t348 = t302 * t84 - t305 * t85;
t41 = t301 * t119 - t304 * t348;
t377 = t91 * t301 + t41;
t258 = t359 * qJD(3);
t376 = t258 * (t297 + t298);
t374 = qJD(1) * t426;
t371 = -pkin(4) * t300 + t488;
t368 = t301 * t100 + t162 * t436 + t169 * t414 + t219 * t390;
t38 = t107 * t301 - t304 * t351;
t86 = t179 * t433 + t181 * t245 + t183 * t246;
t87 = t180 * t433 + t182 * t245 + t184 * t246;
t347 = t302 * t86 - t305 * t87;
t42 = t301 * t120 - t304 * t347;
t364 = -t92 * t301 - t38 - t42;
t358 = rSges(5,1) * t175 + rSges(5,2) * t174;
t53 = t302 * t71 + t305 * t70;
t59 = t302 * t85 + t305 * t84;
t60 = t302 * t87 + t305 * t86;
t346 = Icges(4,1) * t304 - t469;
t342 = -Icges(4,2) * t301 + t468;
t338 = Icges(4,5) * t304 - Icges(4,6) * t301;
t330 = t185 * t305 + t186 * t302;
t326 = -t446 - t450;
t321 = rSges(3,3) * t305 + t482 * t302;
t198 = (-rSges(5,1) * t300 - rSges(5,2) * t303) * t411 + (-t356 * t301 + t476) * qJD(3);
t320 = t198 * t302 + t232 * t417;
t318 = t326 * t302;
t317 = qJD(3) * t346;
t316 = qJD(3) * t342;
t312 = t488 * t302 + t305 * t394;
t11 = -qJD(1) * t351 + t18 * t305 + t19 * t302;
t20 = t149 * t165 + t150 * t167 + t163 * t315 + t234 * t96 + t235 * t98 - t436 * t94;
t21 = t149 * t166 + t150 * t168 + t164 * t315 + t234 * t95 + t235 * t97 - t436 * t93;
t12 = -qJD(1) * t352 + t20 * t305 + t21 * t302;
t47 = t149 * t217 + t150 * t218 - t158 * t436 + t159 * t234 + t160 * t235 + t216 * t315;
t5 = (qJD(3) * t352 + t47) * t301 + (-qJD(1) * t53 + qJD(3) * t106 - t20 * t302 + t21 * t305) * t304;
t309 = t11 * t383 + (-qJD(1) * t349 + t24 * t305 + t25 * t302) * t486 + t302 * t4 / 0.2e1 + t5 * t483 - t37 * t419 / 0.2e1 + t38 * t417 / 0.2e1 + t350 * t414 / 0.2e1 + t12 * t384 + t507 * t54 + t506 * t53;
t308 = t471 + (t24 + t47) * t384 + (t25 + t46) * t383 + (t107 + t81) * t507 + (t106 + t80) * t506;
t307 = (-t302 * t5 + (-t302 * t38 - t305 * t37) * qJD(1)) * t304 - t38 * t385 + t382;
t292 = t305 * qJ(2);
t257 = t302 * t266;
t250 = t266 * t419;
t249 = pkin(7) * t433 - t284;
t240 = t305 * t249;
t239 = -rSges(3,2) * t305 + rSges(3,3) * t302 + t422;
t238 = t292 + t321;
t233 = t474 - t319;
t211 = t288 + (t482 * t305 + (-rSges(3,3) - qJ(2)) * t302) * qJD(1);
t210 = qJD(1) * t321 + t423;
t209 = t219 * t433;
t206 = t395 + t231;
t205 = t410 * t302 + t292 + t319;
t202 = (-t232 - t266) * t305;
t201 = t232 * t302 + t257;
t200 = t305 * (qJD(1) * t282 - t396);
t191 = t495 * qJD(1) + t338 * t415;
t190 = -t338 * t413 + t421;
t152 = t162 * t433;
t141 = (-t266 - t426) * t305;
t140 = t257 + t375;
t139 = t302 * t394 + t282 + t395 + t425;
t138 = t284 + t292 + t312 + t357;
t137 = -t186 * t301 + t232 * t433;
t136 = t185 * t301 + t232 * t436;
t135 = -t302 * t495 - t305 * t326;
t134 = t223 * t302 - t501;
t132 = -t305 * t495 + t318;
t131 = t223 * t305 + t302 * t327;
t129 = -t170 * t301 + t209;
t127 = t169 + t395 + t399;
t126 = t371 * t302 + t496 * t305 + t292 + t354;
t123 = t330 * t304;
t122 = t320 + t424;
t121 = t393 + t250 + (-t198 - t259) * t305;
t116 = (-t169 * t305 - t170 * t302) * t304;
t115 = -rSges(5,3) * t390 + t400;
t114 = -rSges(5,3) * t505 + t358;
t101 = t427 * t302 + t240 + t458;
t89 = t215 * t433 - t428 * t301 + t209;
t88 = t187 * t301 + t215 * t436 + t128;
t83 = rSges(5,3) * t385 + t288 + (t488 * t305 + (-qJ(2) + t476 - t480) * t302) * qJD(1) - t358 + t396;
t82 = qJD(1) * t312 + t397 + t400 + t423;
t75 = t430 * t302 + t305 * t374 + t424;
t74 = t250 + t302 * t374 + (-t259 - t430) * t305;
t69 = (-t428 * t302 + t429 * t305) * t304;
t68 = t401 * t302 + t428 * t305 + t240;
t67 = -t369 + t288 + (-t370 + (t473 * t301 + t434) * qJD(3)) * t305 + (t371 * t305 + (-qJ(2) - t496) * t302) * qJD(1) - t355;
t66 = ((-t409 - t412) * t301 + t371 * qJD(1)) * t302 + t100 + t365 + t423;
t65 = (-t232 * t415 + t115) * t301 + (qJD(3) * t185 + t320) * t304;
t64 = (-t232 * t413 - t114) * t301 + (-qJD(3) * t186 + t198 * t305 - t393) * t304;
t62 = -t219 * t386 + t368;
t61 = -t219 * t391 - t301 * t99 + t152 + (-t170 * t304 - t219 * t438) * qJD(3);
t48 = t330 * t416 + (-t114 * t302 - t115 * t305 + (t185 * t302 - t458) * qJD(1)) * t304;
t45 = t114 * t305 + t200 + (-t115 + t313) * t302 + (t427 * t305 + (-t186 - t249) * t302) * qJD(1);
t39 = (-t302 * t99 + (-qJD(1) * t170 - t100) * t305) * t304 + t402;
t33 = t118 * t301 + (t199 * t302 + t215 * t417) * t304 + (t187 * t304 - t301 * t375) * qJD(3) + t368;
t32 = t152 + (-t426 * t413 - t472) * t301 + (-qJD(1) * t375 - t428 * qJD(3) + t199 * t305) * t304;
t29 = -t108 * t436 + t110 * t243 + t112 * t244 + t176 * t182 + t177 * t184 + t180 * t315;
t28 = -t109 * t436 + t111 * t243 + t113 * t244 + t176 * t181 + t177 * t183 + t179 * t315;
t27 = t108 * t433 + t110 * t245 + t112 * t246 + t174 * t182 + t175 * t184 - t180 * t505;
t26 = t109 * t433 + t111 * t245 + t113 * t246 + t174 * t181 + t175 * t183 - t179 * t505;
t17 = t200 + t472 * t305 + (t313 + t431) * t302 + (t401 * t305 + (-t249 - t428) * t302) * qJD(1);
t16 = (t187 * t305 + t188 * t302) * t416 + ((qJD(1) * t187 - t472) * t302 + (-qJD(1) * t428 + t431) * t305) * t304 + t402;
t15 = -qJD(1) * t348 + t28 * t305 + t29 * t302;
t14 = -qJD(1) * t347 + t26 * t305 + t27 * t302;
t8 = (qJD(3) * t348 + t57) * t301 + (-qJD(1) * t59 + qJD(3) * t119 - t28 * t302 + t29 * t305) * t304;
t7 = (qJD(3) * t347 + t56) * t301 + (-qJD(1) * t60 + qJD(3) * t120 - t26 * t302 + t27 * t305) * t304;
t1 = [t341 * t416 - t345 * t414 + 0.2e1 * m(3) * (t210 * t239 + t211 * t238) + (t142 * t206 + t143 * t205) * t491 + (t138 * t83 + t139 * t82) * t490 + (t126 * t67 + t127 * t66) * t489 - t301 * t317 + t311 + t310 - t228 * t387 - t411 * t454 - t289 * t218 * t443 + (-t289 * t159 - t316 - t405 - t442) * t304; m(6) * (t302 * t67 - t305 * t66 + (t126 * t305 + t127 * t302) * qJD(1)) + m(5) * (t302 * t83 - t305 * t82 + (t138 * t305 + t139 * t302) * qJD(1)) + m(4) * ((t205 * t305 + t206 * t302) * qJD(1) + t499) + m(3) * (-t210 * t305 + t211 * t302 + (t238 * t305 + t239 * t302) * qJD(1)); 0; ((t494 * qJD(1) + t302 * t316) * t487 + (t493 * qJD(1) + t302 * t317) * t484 + t47 / 0.2e1 + t24 / 0.2e1 + (-t452 / 0.2e1 - t448 / 0.2e1) * qJD(3) - t407) * t305 + ((qJD(1) * t226 - t342 * t413) * t487 + (qJD(1) * t229 - t346 * t413) * t484 + t25 / 0.2e1 + t46 / 0.2e1 + (t450 / 0.2e1 + t446 / 0.2e1) * qJD(3) + t406) * t302 + m(5) * (t121 * t139 + t122 * t138 + t201 * t83 + t202 * t82) + m(6) * (t126 * t75 + t127 * t74 + t140 * t67 + t141 * t66) + m(4) * (t499 * t265 - (t205 * t302 - t206 * t305) * t258) - (t297 / 0.2e1 + t298 / 0.2e1) * t337 * qJD(3) + ((-t80 / 0.2e1 - t106 / 0.2e1 + t453 / 0.2e1 - t447 / 0.2e1 + t206 * t481 - t404) * t302 + (t81 / 0.2e1 + t107 / 0.2e1 + t205 * t481 + t451 / 0.2e1 - t445 / 0.2e1 - t403) * t305) * qJD(1); m(5) * (-t121 * t305 + t122 * t302 + (t201 * t305 + t202 * t302) * qJD(1)) + m(6) * (t302 * t75 - t305 * t74 + (t140 * t305 + t141 * t302) * qJD(1)) - m(4) * t376; (t140 * t75 + t141 * t74 + t17 * t68) * t489 + t305 * t12 + t302 * t11 + t302 * t14 + (t101 * t45 + t121 * t202 + t122 * t201) * t490 + t305 * t15 + t305 * ((t305 * t191 + (t132 + t501) * qJD(1)) * t305 + (-t131 * qJD(1) + (-t414 * t493 + t416 * t494) * t302 + (t190 + (t447 - t453) * qJD(3) + (-t223 + t326) * qJD(1)) * t305) * t302) + t302 * ((t302 * t190 + (-t134 + t318) * qJD(1)) * t302 + (t135 * qJD(1) + (t226 * t416 - t229 * t414 + t421) * t305 + (t191 + (t445 - t451) * qJD(3) + t327 * qJD(1)) * t302) * t305) + ((-t231 * t302 + t233 * t305) * (-t302 * t398 + (-t265 * t298 + t297 * t477) * qJD(3) + ((-t231 + t293) * t305 + (-t233 + t319 + t474) * t302) * qJD(1)) - t265 * t376) * t491 + (-t131 * t305 - t132 * t302 - t53 - t59) * t419 + (t134 * t305 + t135 * t302 + t54 + t60) * t417; t308 + (t406 * t305 + t407 * t302 + (t302 * t403 - t305 * t404) * qJD(1)) * t304 + m(5) * (t136 * t82 + t137 * t83 + t138 * t64 + t139 * t65) + m(6) * (t126 * t32 + t127 * t33 + t66 * t88 + t67 * t89) + (t302 * t404 + t305 * t403) * t416 + t470; m(5) * (t302 * t64 - t305 * t65 + (t136 * t302 + t137 * t305) * qJD(1)) + m(6) * (t302 * t32 - t305 * t33 + (t302 * t88 + t305 * t89) * qJD(1)); m(5) * (t101 * t48 + t121 * t136 + t122 * t137 - t123 * t45 + t201 * t64 + t202 * t65) + m(6) * (t140 * t32 + t141 * t33 + t16 * t68 + t17 * t69 + t74 * t88 + t75 * t89) + (t15 * t485 + t14 * t483 + qJD(3) * (t302 * t92 + t305 * t91) / 0.2e1 + (-t305 * t59 / 0.2e1 + t60 * t485) * qJD(1)) * t304 + t309 + ((qJD(1) * t92 + t30) * t486 + t8 / 0.2e1 + t60 * t379 + qJD(1) * t42 / 0.2e1) * t305 + ((-qJD(1) * t91 + t31) * t486 + t7 / 0.2e1 - qJD(1) * t41 / 0.2e1 + t59 * t378) * t302; (t16 * t69 + t32 * t89 + t33 * t88) * t489 + (-t123 * t48 + t136 * t65 + t137 * t64) * t490 + ((t302 * t377 + t305 * t364) * qJD(3) + t470) * t301 + ((t301 * t31 + t7) * t305 + (-t30 * t301 - t5 - t8) * t302 + (t133 * t301 + (-t302 * t91 + t305 * t92) * t304) * qJD(3) + ((-t37 - t377) * t305 + t364 * t302) * qJD(1)) * t304 + t382; t308 + m(6) * (t126 * t61 + t127 * t62 + t128 * t66 + t129 * t67); m(6) * (t302 * t61 - t305 * t62 + (t128 * t302 + t129 * t305) * qJD(1)); m(6) * (t116 * t17 + t128 * t74 + t129 * t75 + t140 * t61 + t141 * t62 + t39 * t68) + t309; m(6) * (t116 * t16 + t128 * t33 + t129 * t32 + t39 * t69 + t61 * t89 + t62 * t88) + t307; (t116 * t39 + t128 * t62 + t129 * t61) * t489 + t307;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
