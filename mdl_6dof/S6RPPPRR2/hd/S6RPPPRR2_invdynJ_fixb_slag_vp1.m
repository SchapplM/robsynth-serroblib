% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:39
% EndTime: 2019-03-09 01:32:01
% DurationCPUTime: 18.84s
% Computational Cost: add. (22144->856), mult. (20733->1094), div. (0->0), fcn. (18752->10), ass. (0->410)
t323 = cos(qJ(1));
t314 = t323 * pkin(1);
t316 = qJ(1) + pkin(9);
t311 = sin(t316);
t313 = cos(t316);
t239 = rSges(3,1) * t311 + rSges(3,2) * t313;
t321 = sin(qJ(1));
t537 = pkin(1) * t321;
t216 = -t239 - t537;
t315 = pkin(10) + qJ(5);
t312 = cos(t315);
t493 = t312 * t313;
t437 = rSges(7,3) * t493;
t310 = sin(t315);
t320 = sin(qJ(6));
t492 = t313 * t320;
t322 = cos(qJ(6));
t494 = t311 * t322;
t205 = t310 * t492 + t494;
t491 = t313 * t322;
t495 = t311 * t320;
t206 = t310 * t491 - t495;
t471 = t206 * rSges(7,1) - t205 * rSges(7,2);
t116 = -t471 + t437;
t520 = rSges(7,2) * t320;
t525 = rSges(7,1) * t322;
t397 = -t520 + t525;
t518 = rSges(7,3) * t310;
t177 = t312 * t397 + t518;
t448 = qJD(6) * t313;
t421 = t312 * t448;
t453 = qJD(5) * t311;
t213 = t421 + t453;
t531 = t312 * pkin(5);
t245 = pkin(8) * t310 + t531;
t450 = qJD(6) * t310;
t276 = qJD(1) + t450;
t584 = t116 * t276 - t177 * t213 - t245 * t453;
t324 = qJD(1) ^ 2;
t583 = t324 * t314;
t536 = pkin(2) * t311;
t203 = -t310 * t495 + t491;
t204 = t310 * t494 + t492;
t497 = t311 * t312;
t105 = Icges(7,5) * t204 + Icges(7,6) * t203 - Icges(7,3) * t497;
t107 = -Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t493;
t189 = Icges(7,4) * t206;
t110 = Icges(7,2) * t205 + Icges(7,6) * t493 - t189;
t188 = Icges(7,4) * t205;
t112 = Icges(7,1) * t206 - Icges(7,5) * t493 - t188;
t371 = -t205 * t110 - t206 * t112;
t509 = Icges(7,4) * t204;
t108 = Icges(7,2) * t203 - Icges(7,6) * t497 + t509;
t187 = Icges(7,4) * t203;
t111 = Icges(7,1) * t204 - Icges(7,5) * t497 + t187;
t526 = t203 * t108 + t204 * t111;
t582 = t371 + t526 + (-t105 * t311 - t107 * t313) * t312;
t319 = -pkin(7) - qJ(4);
t488 = qJ(4) + t319;
t317 = sin(pkin(10));
t535 = pkin(4) * t317;
t179 = t311 * t488 + t313 * t535;
t292 = qJD(4) * t313;
t296 = t313 * qJ(3);
t237 = -t296 + t536;
t293 = qJD(3) * t311;
t461 = -qJD(1) * t237 + t293;
t432 = t292 + t461;
t455 = qJD(1) * t313;
t283 = qJ(3) * t455;
t462 = t292 + t293;
t433 = t283 + t462;
t428 = t317 * t455;
t454 = qJD(1) * t319;
t467 = pkin(4) * t428 + t311 * t454;
t581 = -qJD(1) * t179 - t432 + t433 + t467;
t395 = -qJ(4) * t311 - t537;
t580 = -t395 - t537;
t449 = qJD(6) * t312;
t451 = qJD(5) * t313;
t214 = -t311 * t449 + t451;
t32 = t105 * t493 + t205 * t108 - t206 * t111;
t374 = Icges(7,5) * t322 - Icges(7,6) * t320;
t171 = Icges(7,3) * t310 + t312 * t374;
t507 = Icges(7,4) * t322;
t376 = -Icges(7,2) * t320 + t507;
t173 = Icges(7,6) * t310 + t312 * t376;
t508 = Icges(7,4) * t320;
t378 = Icges(7,1) * t322 - t508;
t175 = Icges(7,5) * t310 + t312 * t378;
t56 = t171 * t493 + t173 * t205 - t175 * t206;
t579 = -t214 * t32 - t276 * t56;
t242 = -rSges(4,2) * t313 + t311 * rSges(4,3);
t308 = t313 * pkin(2);
t564 = t311 * qJ(3) + t308;
t431 = t314 + t564;
t156 = t242 + t431;
t264 = pkin(8) * t493;
t404 = -t536 - t537;
t577 = t311 * t319 - t264 + t296 + t404 + t471;
t370 = t110 * t320 + t112 * t322;
t43 = t107 * t310 - t312 * t370;
t31 = -t107 * t497 + t203 * t110 - t112 * t204;
t452 = qJD(5) * t312;
t574 = t310 * t455 + t311 * t452;
t257 = Icges(6,4) * t497;
t499 = t310 * t311;
t506 = Icges(6,5) * t313;
t161 = Icges(6,1) * t499 + t257 + t506;
t510 = Icges(6,4) * t312;
t379 = Icges(6,1) * t310 + t510;
t162 = -Icges(6,5) * t311 + t313 * t379;
t231 = -Icges(6,2) * t310 + t510;
t184 = t231 * t313;
t334 = t311 * (t162 + t184) - t313 * (-Icges(6,2) * t499 + t161 + t257);
t511 = Icges(6,4) * t310;
t377 = Icges(6,2) * t312 + t511;
t159 = Icges(6,6) * t313 + t311 * t377;
t160 = -Icges(6,6) * t311 + t313 * t377;
t233 = Icges(6,1) * t312 - t511;
t185 = t233 * t311;
t186 = t233 * t313;
t335 = t311 * (t160 - t186) - t313 * (t159 - t185);
t573 = -t335 * t310 + t334 * t312;
t469 = t231 + t379;
t470 = -t377 + t233;
t572 = (t310 * t469 - t312 * t470) * qJD(1);
t55 = -t171 * t497 + t173 * t203 + t175 * t204;
t570 = t213 * t31 + t55 * t276;
t538 = rSges(7,3) + pkin(8);
t569 = t310 * t538;
t170 = Icges(7,3) * t312 - t310 * t374;
t364 = t173 * t320 - t175 * t322;
t372 = t108 * t320 - t111 * t322;
t325 = t213 * (-t171 * t313 + t370) + t214 * (t171 * t311 + t372) + t276 * (t170 + t364);
t568 = t325 * t312;
t366 = t160 * t312 + t162 * t310;
t565 = t366 * t313;
t243 = t313 * rSges(3,1) - rSges(3,2) * t311;
t217 = t243 + t314;
t307 = t312 * pkin(8);
t534 = pkin(5) * t310;
t244 = t307 - t534;
t446 = qJD(1) * qJD(3);
t457 = qJD(1) * t311;
t464 = t283 + t293;
t562 = qJD(1) * (-pkin(2) * t457 + t464) + qJDD(1) * t564 + t311 * t446;
t496 = t311 * t317;
t318 = cos(pkin(10));
t522 = rSges(5,2) * t318;
t168 = rSges(5,1) * t496 + t313 * rSges(5,3) + t311 * t522;
t101 = -qJD(5) * t186 + (t311 * t379 + t506) * qJD(1);
t375 = Icges(6,5) * t310 + Icges(6,6) * t312;
t158 = -Icges(6,3) * t311 + t313 * t375;
t459 = qJD(1) * t158;
t82 = t160 * t310 - t162 * t312;
t99 = qJD(1) * t159 - qJD(5) * t184;
t561 = qJD(5) * t82 + t101 * t310 + t312 * t99 + t459;
t559 = t276 * t320 - t322 * t452;
t558 = t276 * t322 + t320 * t452;
t220 = t377 * qJD(5);
t221 = t379 * qJD(5);
t229 = Icges(6,5) * t312 - Icges(6,6) * t310;
t362 = t231 * t310 - t233 * t312;
t557 = qJD(1) * t229 + qJD(5) * t362 + t220 * t312 + t221 * t310;
t100 = qJD(1) * t160 + t231 * t453;
t102 = qJD(1) * t162 + qJD(5) * t185;
t367 = t159 * t310 - t161 * t312;
t157 = Icges(6,3) * t313 + t311 * t375;
t460 = qJD(1) * t157;
t556 = qJD(5) * t367 - t100 * t312 - t102 * t310 + t460;
t211 = (-Icges(7,2) * t322 - t508) * t312;
t329 = t213 * (Icges(7,2) * t206 - t112 + t188) + t214 * (-Icges(7,2) * t204 + t111 + t187) + t276 * (t175 + t211);
t212 = (-Icges(7,1) * t320 - t507) * t312;
t554 = t213 * (-Icges(7,1) * t205 + t110 - t189) + t214 * (-Icges(7,1) * t203 + t108 + t509) + t276 * (t173 - t212);
t444 = qJD(1) * qJD(5);
t224 = qJDD(5) * t311 + t313 * t444;
t426 = t310 * t451;
t456 = qJD(1) * t312;
t345 = -t311 * t456 - t426;
t442 = qJDD(6) * t312;
t117 = qJD(6) * t345 + t313 * t442 + t224;
t553 = t117 / 0.2e1;
t288 = qJDD(5) * t313;
t118 = -qJD(1) * t421 + t288 + (-t442 + (-qJD(1) + t450) * qJD(5)) * t311;
t552 = t118 / 0.2e1;
t551 = -t213 / 0.2e1;
t550 = t213 / 0.2e1;
t549 = -t214 / 0.2e1;
t548 = t214 / 0.2e1;
t218 = qJD(5) * t449 + qJDD(6) * t310 + qJDD(1);
t547 = t218 / 0.2e1;
t546 = t224 / 0.2e1;
t225 = -t311 * t444 + t288;
t545 = t225 / 0.2e1;
t544 = -t276 / 0.2e1;
t543 = t276 / 0.2e1;
t542 = t310 / 0.2e1;
t541 = t311 / 0.2e1;
t540 = -t313 / 0.2e1;
t533 = g(1) * t313;
t532 = g(2) * t313;
t427 = t310 * t453;
t429 = t312 * t455;
t346 = t427 - t429;
t409 = qJD(1) * t310 + qJD(6);
t361 = t320 * t409;
t95 = -t311 * t558 - t313 * t361;
t360 = t409 * t322;
t96 = -t311 * t559 + t313 * t360;
t47 = Icges(7,5) * t96 + Icges(7,6) * t95 + Icges(7,3) * t346;
t49 = Icges(7,4) * t96 + Icges(7,2) * t95 + Icges(7,6) * t346;
t51 = Icges(7,1) * t96 + Icges(7,4) * t95 + Icges(7,5) * t346;
t7 = (qJD(5) * t372 + t47) * t310 + (qJD(5) * t105 - t320 * t49 + t322 * t51 + (-t108 * t322 - t111 * t320) * qJD(6)) * t312;
t530 = t7 * t214;
t93 = -t311 * t361 + t313 * t558;
t94 = t311 * t360 + t313 * t559;
t46 = Icges(7,5) * t94 + Icges(7,6) * t93 + Icges(7,3) * t345;
t48 = Icges(7,4) * t94 + Icges(7,2) * t93 + Icges(7,6) * t345;
t50 = Icges(7,1) * t94 + Icges(7,4) * t93 + Icges(7,5) * t345;
t8 = (qJD(5) * t370 + t46) * t310 + (qJD(5) * t107 - t320 * t48 + t322 * t50 + (-t110 * t322 + t112 * t320) * qJD(6)) * t312;
t529 = t8 * t213;
t528 = -pkin(2) - qJ(4);
t210 = (-Icges(7,5) * t320 - Icges(7,6) * t322) * t312;
t121 = qJD(5) * t170 + qJD(6) * t210;
t172 = Icges(7,6) * t312 - t310 * t376;
t122 = qJD(5) * t172 + qJD(6) * t211;
t174 = Icges(7,5) * t312 - t310 * t378;
t123 = qJD(5) * t174 + qJD(6) * t212;
t22 = (qJD(5) * t364 + t121) * t310 + (qJD(5) * t171 - t122 * t320 + t123 * t322 + (-t173 * t322 - t175 * t320) * qJD(6)) * t312;
t70 = t171 * t310 - t312 * t364;
t527 = t70 * t218 + t22 * t276;
t521 = rSges(6,2) * t310;
t519 = rSges(4,3) * t313;
t240 = rSges(6,1) * t312 - t521;
t192 = t240 * t313;
t304 = t313 * rSges(6,3);
t399 = rSges(6,1) * t310 + rSges(6,2) * t312;
t103 = -qJD(5) * t192 + (t311 * t399 + t304) * qJD(1);
t222 = t399 * qJD(5);
t502 = qJ(4) * t313;
t394 = -t502 - t314;
t445 = qJD(1) * qJD(4);
t465 = qJDD(3) * t311 + t313 * t446;
t331 = qJDD(4) * t313 - 0.2e1 * t311 * t445 + t324 * t394 + t465;
t300 = t311 * rSges(6,3);
t164 = t313 * t399 - t300;
t359 = -t237 + t395;
t350 = t179 + t359;
t339 = t164 + t350;
t294 = qJD(3) * t313;
t190 = qJD(1) * t564 - t294;
t269 = t313 * t454;
t274 = pkin(4) * t496;
t476 = t269 - (t274 - t502) * qJD(1) - t190;
t25 = -t222 * t453 + t224 * t240 + (-t103 + t476) * qJD(1) + t339 * qJDD(1) + t331;
t516 = t25 * t311;
t435 = rSges(6,1) * t574 + rSges(6,2) * t429;
t104 = (-rSges(6,3) * qJD(1) - qJD(5) * t521) * t311 + t435;
t163 = rSges(6,1) * t499 + rSges(6,2) * t497 + t304;
t180 = -t313 * t488 + t274;
t309 = qJDD(1) * t314;
t328 = qJDD(1) * t502 + qJDD(4) * t311 + 0.2e1 * t313 * t445 + t324 * t395 + t309 + t562;
t326 = qJD(1) * (qJ(4) * t457 + t467) + qJDD(1) * t180 + t328;
t26 = qJD(1) * t104 + qJDD(1) * t163 - t225 * t240 + (qJD(5) * t222 - qJDD(3)) * t313 + t326;
t515 = t26 * t313;
t303 = t312 * rSges(7,3);
t42 = t105 * t310 - t312 * t372;
t513 = t42 * t118;
t512 = t43 * t117;
t500 = t229 * t313;
t498 = t310 * t313;
t181 = t311 * t229;
t363 = t231 * t312 + t233 * t310;
t88 = t313 * t363 - t181;
t489 = t88 * qJD(1);
t472 = t204 * rSges(7,1) + t203 * rSges(7,2);
t114 = -rSges(7,3) * t497 + t472;
t262 = pkin(5) * t499;
t197 = -pkin(8) * t497 + t262;
t483 = -t114 - t197;
t199 = pkin(5) * t498 - t264;
t482 = t116 - t199;
t215 = (-rSges(7,1) * t320 - rSges(7,2) * t322) * t312;
t126 = qJD(6) * t215 + (-t310 * t397 + t303) * qJD(5);
t223 = t244 * qJD(5);
t481 = t126 + t223;
t473 = t177 + t245;
t198 = pkin(5) * t497 + pkin(8) * t499;
t468 = rSges(5,1) * t428 + t455 * t522;
t466 = t310 * t520 + t303;
t463 = rSges(4,2) * t457 + rSges(4,3) * t455;
t458 = qJD(1) * t375;
t447 = -m(5) - m(6) - m(7);
t443 = qJDD(3) * t313;
t441 = -rSges(5,3) + t528;
t440 = t96 * rSges(7,1) + t95 * rSges(7,2) + rSges(7,3) * t427;
t439 = t312 * t525;
t438 = t312 * t520;
t436 = -m(4) + t447;
t63 = t313 * t157 + t159 * t497 + t161 * t499;
t64 = -t313 * t158 - t160 * t497 - t162 * t499;
t434 = pkin(5) * t574 + pkin(8) * t427;
t424 = t312 * t451;
t420 = -pkin(5) - t525;
t419 = -t457 / 0.2e1;
t418 = -t456 / 0.2e1;
t417 = t455 / 0.2e1;
t416 = -t453 / 0.2e1;
t415 = t453 / 0.2e1;
t414 = -t451 / 0.2e1;
t413 = t451 / 0.2e1;
t412 = rSges(4,2) * t311 + t519 - t537;
t411 = t296 - t537;
t410 = qJD(4) * t311 - t294;
t407 = -t324 * t537 + t309;
t405 = rSges(7,1) * t94 + rSges(7,2) * t93;
t403 = -t314 - t308;
t402 = t269 - t410;
t273 = rSges(2,1) * t323 - rSges(2,2) * t321;
t272 = rSges(2,1) * t321 + rSges(2,2) * t323;
t400 = rSges(5,1) * t317 + t522;
t125 = -pkin(8) * t429 + t434;
t53 = -rSges(7,3) * t429 + t440;
t10 = t218 * t114 - t214 * t126 + (-qJD(5) * t223 - qJDD(3)) * t313 + t326 + t276 * t53 - t225 * t245 + qJDD(1) * t197 + qJD(1) * t125 - t118 * t177;
t124 = t345 * pkin(8) + (t310 * t457 - t424) * pkin(5);
t338 = t199 + t350;
t52 = rSges(7,3) * t345 + t405;
t11 = t223 * t453 - t116 * t218 + t117 * t177 + t126 * t213 + t224 * t245 - t276 * t52 + (-t124 + t476) * qJD(1) + t338 * qJDD(1) + t331;
t392 = t10 * t311 + t11 * t313;
t368 = t159 * t312 + t161 * t310;
t336 = qJD(1) * t368 + qJD(5) * t181 + t459;
t337 = -qJD(1) * t366 - qJD(5) * t500 + t460;
t391 = (t336 * t311 + t313 * t556) * t313 + (t337 * t311 - t313 * t561) * t311;
t390 = (-t311 * t556 + t336 * t313) * t313 + (t311 * t561 + t337 * t313) * t311;
t30 = -t105 * t497 + t526;
t389 = t30 * t313 + t31 * t311;
t388 = t30 * t311 - t31 * t313;
t33 = t107 * t493 - t371;
t387 = t311 * t33 + t313 * t32;
t386 = t311 * t32 - t313 * t33;
t385 = t311 * t43 + t313 * t42;
t384 = t311 * t42 - t313 * t43;
t207 = t240 * t453;
t60 = qJD(1) * t339 + t207 + t462;
t358 = t564 - t394;
t349 = t180 + t358;
t61 = -t240 * t451 + (t163 + t349) * qJD(1) + t410;
t383 = t311 * t60 - t313 * t61;
t382 = t311 * t64 + t313 * t63;
t151 = t311 * t157;
t65 = -t313 * t368 + t151;
t66 = -t311 * t158 + t565;
t381 = t311 * t66 + t313 * t65;
t373 = t103 * t313 - t104 * t311;
t369 = t114 * t313 + t116 * t311;
t365 = -t163 * t311 - t164 * t313;
t149 = rSges(7,3) * t499 + (-t438 + t439) * t311;
t352 = -t313 * t319 + t274 + t431;
t169 = -t311 * rSges(5,3) + t313 * t400;
t351 = t169 + t359;
t348 = -t303 + t534 + t535;
t347 = t399 + t535;
t343 = t105 * t214 + t107 * t213 + t171 * t276;
t342 = (Icges(7,5) * t203 - Icges(7,6) * t204) * t214 + (Icges(7,5) * t205 + Icges(7,6) * t206) * t213 + t210 * t276;
t333 = t363 * qJD(1) - t375 * qJD(5);
t34 = qJD(1) * t338 + t462 - t584;
t35 = -t245 * t451 + t114 * t276 - t177 * t214 + (t197 + t349) * qJD(1) + t410;
t37 = -t114 * t213 + t116 * t214 + qJD(2) + (-t197 * t311 - t199 * t313) * qJD(5);
t327 = t37 * t369 + (-t311 * t35 - t313 * t34) * t177;
t235 = t313 * t438;
t200 = t245 * t313;
t191 = t240 * t311;
t176 = -t310 * t525 + t466;
t150 = t235 + (-t439 - t518) * t313;
t148 = t175 * t313;
t147 = t175 * t311;
t146 = t173 * t313;
t145 = t173 * t311;
t136 = rSges(7,1) * t205 + rSges(7,2) * t206;
t135 = rSges(7,1) * t203 - rSges(7,2) * t204;
t87 = t311 * t363 + t500;
t85 = t87 * qJD(1);
t84 = (t168 + t358) * qJD(1) + t410;
t83 = qJD(1) * t351 + t462;
t80 = qJD(5) * t365 + qJD(2);
t69 = qJD(1) * t463 + qJDD(1) * t242 + t407 - t443 + t562;
t68 = -t583 + (-t237 + t412) * qJDD(1) + (-qJD(1) * t242 - t190) * qJD(1) + t465;
t45 = -t443 + qJDD(1) * t168 + qJD(1) * (-rSges(5,3) * t457 + t468) + t328;
t44 = t351 * qJDD(1) + (-qJD(1) * t168 - t190) * qJD(1) + t331;
t41 = -t311 * t557 + t333 * t313;
t40 = t333 * t311 + t313 * t557;
t39 = qJD(5) * t366 + t101 * t312 - t310 * t99;
t38 = -qJD(5) * t368 - t100 * t310 + t102 * t312;
t36 = qJD(5) * t373 - t163 * t224 - t164 * t225 + qJDD(2);
t24 = qJD(5) * t381 - t489;
t23 = qJD(5) * t382 + t85;
t16 = -t121 * t497 + t122 * t203 + t123 * t204 + t171 * t346 + t173 * t95 + t175 * t96;
t15 = t121 * t493 + t122 * t205 - t123 * t206 + t171 * t345 + t173 * t93 + t175 * t94;
t14 = t213 * t43 + t214 * t42 + t276 * t70;
t13 = t213 * t33 - t579;
t12 = t214 * t30 + t570;
t9 = -t114 * t117 + t116 * t118 - t197 * t224 - t199 * t225 - t213 * t53 + t214 * t52 + qJDD(2) + (t124 * t313 - t125 * t311) * qJD(5);
t6 = t107 * t346 + t110 * t95 - t112 * t96 + t203 * t48 + t204 * t50 - t46 * t497;
t5 = t105 * t346 + t108 * t95 + t111 * t96 + t203 * t49 + t204 * t51 - t47 * t497;
t4 = t107 * t345 + t110 * t93 - t112 * t94 + t205 * t48 - t206 * t50 + t46 * t493;
t3 = t105 * t345 + t108 * t93 + t111 * t94 + t205 * t49 - t206 * t51 + t47 * t493;
t2 = t117 * t31 + t118 * t30 + t16 * t276 + t213 * t6 + t214 * t5 + t218 * t55;
t1 = t117 * t33 + t118 * t32 + t15 * t276 + t213 * t4 + t214 * t3 + t218 * t56;
t17 = [(-t83 * t410 + t84 * (t433 + t468) + ((-t321 * t84 - t323 * t83) * pkin(1) + t83 * t441 * t313 + (t83 * (-qJ(3) - t400) + t84 * t441) * t311) * qJD(1) - (-t83 + (t169 + t395) * qJD(1) + t432) * t84 + (t45 - g(2)) * (t168 + t431 + t502) + (t44 - g(1)) * (t311 * t528 + t169 + t411)) * m(5) + (t15 + t12) * t550 + t527 - t224 * t88 / 0.2e1 + t529 / 0.2e1 + (t489 + (t311 ^ 2 * t158 + (-t151 + t64 + (t158 + t368) * t313) * t313) * qJD(5) + t24) * t414 + ((t68 - g(1)) * (t519 + (rSges(4,2) - pkin(2)) * t311 + t411) + (t69 - g(2)) * t156 + (-t461 + t463 + t464 + (t404 - t412) * qJD(1)) * (qJD(1) * t156 - t294)) * m(4) + (-t367 + t87) * t545 + (t38 + t41) * t413 + t82 * t546 + t16 * t548 + t55 * t552 + t56 * t553 + (-qJD(5) * t363 + t220 * t310 - t221 * t312) * qJD(1) - m(2) * (-g(1) * t272 + g(2) * t273) + t530 / 0.2e1 + t512 / 0.2e1 + t513 / 0.2e1 + (-g(1) * t577 - t348 * t533 + (t10 - g(2)) * (-t497 * t538 + t262 + t352 + t472) + (t313 * t348 + t577) * t11 + (t402 - t405 + (t531 + t569) * t451 + (t403 + (-qJ(3) - t348 + t307) * t311) * qJD(1)) * t34 + (t34 + t434 + t440 + (-t199 - t264 - t437 + t580 - t536) * qJD(1) + t581 + t584) * t35) * m(7) + (m(3) * (t216 ^ 2 + t243 * t217) + Icges(5,1) * t318 ^ 2 + (-0.2e1 * Icges(5,4) * t318 + Icges(5,2) * t317) * t317 + m(2) * (t272 ^ 2 + t273 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1) - t362) * qJDD(1) + ((t26 - g(2)) * (t352 + t163) + (t25 - g(1)) * (-t300 + (-pkin(2) + t319) * t311 + t347 * t313 + t411) + (rSges(6,1) * t424 - rSges(6,2) * t426 + t402 + (t403 - t304 + (-qJ(3) - t347) * t311) * qJD(1)) * t60 + (-rSges(6,2) * t427 - t207 + t435 + t60 + ((-rSges(6,3) - pkin(2)) * t311 - t164 + t580) * qJD(1) + t581) * t61) * m(6) + ((t33 + t582) * t214 + t570) * t551 + (t13 + (-t30 + t582) * t213 + t579) * t549 + (t85 + ((-t65 + t151 + t64) * t311 + (t66 - t565 + (t158 - t368) * t311 + t63) * t313) * qJD(5)) * t416 + ((-t239 * t324 - g(2) + t407) * t217 + (-t583 + (-0.2e1 * t243 - t314 + t217) * t324 - g(1)) * t216) * m(3) + (t39 + t40 + t23) * t415; m(6) * t36 + m(7) * t9 + (m(3) + m(4) + m(5)) * qJDD(2) + (-m(3) + t436) * g(3); t436 * (g(1) * t311 - t532) + 0.2e1 * (t10 * t540 + t11 * t541) * m(7) + 0.2e1 * (t516 / 0.2e1 - t515 / 0.2e1) * m(6) + 0.2e1 * (t44 * t541 + t45 * t540) * m(5) + 0.2e1 * (t540 * t69 + t541 * t68) * m(4); t447 * (g(2) * t311 + t533) + m(5) * (t311 * t45 + t313 * t44) + m(6) * (t25 * t313 + t26 * t311) + m(7) * t392; (t448 * t542 + t417) * t13 + (-t311 * t450 / 0.2e1 + t419) * t12 - qJD(1) * ((-t310 * t470 - t312 * t469) * qJD(1) + (t310 * t334 + t312 * t335) * qJD(5)) / 0.2e1 - t14 * t449 / 0.2e1 + t24 * t417 + t23 * t419 + ((t145 * t205 - t147 * t206) * t214 + (-t146 * t205 + t148 * t206) * t213 + (t172 * t205 - t174 * t206) * t276 + (t312 * t56 + t32 * t499) * qJD(6) + ((-qJD(6) * t33 - t343) * t310 + t568) * t313) * t551 + (((-t145 * t320 + t147 * t322 + t105) * t214 + (t146 * t320 - t148 * t322 + t107) * t213 + (-t172 * t320 + t174 * t322 + t171) * t276 + t70 * qJD(6)) * t312 + (qJD(6) * t384 + t325) * t310) * t544 + ((t145 * t203 + t147 * t204) * t214 + (-t146 * t203 - t148 * t204) * t213 + (t172 * t203 + t174 * t204) * t276 + (-t31 * t498 + t312 * t55) * qJD(6) + ((qJD(6) * t30 + t343) * t310 - t568) * t311) * t549 + ((-t65 * t311 + t66 * t313) * qJD(1) + t391) * t415 + qJDD(1) * (t311 * t82 - t313 * t367) / 0.2e1 + qJD(1) * (t311 * t39 + t313 * t38 + (t311 * t367 + t82 * t313) * qJD(1)) / 0.2e1 + (-g(1) * t191 + g(2) * t192 + g(3) * t399 - (t191 * t61 + t192 * t60) * qJD(1) - (t80 * (-t191 * t311 - t192 * t313) - t383 * t399) * qJD(5) + t36 * t365 + t80 * ((-t163 * t313 + t164 * t311) * qJD(1) + t373) - t383 * t222 + (t516 - t515 + (t311 * t61 + t313 * t60) * qJD(1)) * t240) * m(6) + (-qJD(1) * t384 + t311 * t8 + t313 * t7) * t543 + t382 * t545 + t381 * t546 + t385 * t547 + (-qJD(1) * t388 + t311 * t6 + t313 * t5) * t548 + (-qJD(1) * t386 + t3 * t313 + t311 * t4) * t550 + t389 * t552 + t387 * t553 + (qJD(1) * t41 + qJD(5) * t390 + qJDD(1) * t87 + t224 * t64 + t225 * t63 + t2) * t313 / 0.2e1 + ((-t63 * t311 + t64 * t313) * qJD(1) + t390) * t413 + (qJD(1) * t40 + qJD(5) * t391 - qJDD(1) * t88 + t224 * t66 + t225 * t65 + t1) * t541 + (-g(1) * (t149 + t198) - g(2) * t235 - g(3) * (t310 * t420 + t307 + t466) - (t312 * t420 - t569) * t532 - t34 * (qJD(1) * t200 - t150 * t276 + t176 * t213 + t244 * t453) - t35 * (qJD(1) * t198 + t149 * t276 - t176 * t214 - t244 * t451) - t37 * (-t149 * t213 + t150 * t214 - t198 * t453 - t200 * t451) - ((t114 * t35 - t116 * t34) * t312 + t327 * t310) * qJD(6) + (-t10 * t473 - t35 * t481 + t9 * t482 + t37 * (t124 + t52) + (t34 * t473 + t37 * t483) * qJD(1)) * t313 + (t11 * t473 + t34 * t481 + t9 * t483 + t37 * (-t125 - t53) + (t35 * t473 - t37 * t482) * qJD(1)) * t311) * m(7) + ((t181 * t451 - t458) * t313 + (-t572 + (-t313 * t500 - t573) * qJD(5)) * t311) * t414 + ((-t453 * t500 - t458) * t311 + (t572 + (t311 * t181 + t573) * qJD(5)) * t313) * t416; -t2 * t497 / 0.2e1 + (t310 * t55 - t312 * t388) * t552 + ((qJD(5) * t388 + t16) * t310 + (-qJD(1) * t389 + qJD(5) * t55 - t311 * t5 + t313 * t6) * t312) * t548 + t1 * t493 / 0.2e1 + (t310 * t56 - t312 * t386) * t553 + ((qJD(5) * t386 + t15) * t310 + (-qJD(1) * t387 + qJD(5) * t56 - t3 * t311 + t313 * t4) * t312) * t550 + t14 * t452 / 0.2e1 + (t512 + t513 + t527 + t529 + t530) * t542 + (t310 * t70 - t312 * t384) * t547 + ((qJD(5) * t384 + t22) * t310 + (-qJD(1) * t385 + qJD(5) * t70 - t311 * t7 + t313 * t8) * t312) * t543 + (t203 * t329 - t204 * t554 - t342 * t497) * t549 + (t329 * t205 + t206 * t554 + t342 * t493) * t551 + (t342 * t310 + (-t320 * t329 - t322 * t554) * t312) * t544 + (t310 * t414 + t311 * t418) * t13 + (t310 * t415 + t313 * t418) * t12 + ((qJD(5) * t327 + t10 * t114 - t11 * t116 - t34 * t52 + t35 * t53) * t310 + (t34 * (-qJD(5) * t116 + t126 * t313) + t35 * (qJD(5) * t114 + t126 * t311) - t9 * t369 + t37 * (t114 * t457 - t116 * t455 - t311 * t52 - t313 * t53) + ((-t311 * t34 + t313 * t35) * qJD(1) + t392) * t177) * t312 - t34 * (-t136 * t276 + t213 * t215) - t35 * (t135 * t276 - t214 * t215) - t37 * (-t135 * t213 + t136 * t214) - g(1) * t135 - g(2) * t136 - g(3) * t215) * m(7);];
tau  = t17;
