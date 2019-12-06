% Calculate vector of inverse dynamics joint torques for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:31:03
% DurationCPUTime: 21.72s
% Computational Cost: add. (12071->525), mult. (15946->698), div. (0->0), fcn. (15354->6), ass. (0->250)
t525 = Icges(5,4) + Icges(6,4);
t267 = pkin(7) + qJ(2);
t264 = sin(t267);
t265 = cos(t267);
t272 = cos(qJ(4));
t269 = cos(pkin(8));
t271 = sin(qJ(4));
t386 = t269 * t271;
t201 = t264 * t386 + t265 * t272;
t385 = t269 * t272;
t388 = t265 * t271;
t202 = t264 * t385 - t388;
t171 = Icges(6,4) * t201;
t268 = sin(pkin(8));
t394 = t264 * t268;
t118 = -Icges(6,1) * t202 - Icges(6,5) * t394 + t171;
t174 = Icges(5,4) * t201;
t121 = -Icges(5,1) * t202 - Icges(5,5) * t394 + t174;
t461 = t118 + t121;
t172 = Icges(6,4) * t202;
t111 = -Icges(6,2) * t201 + Icges(6,6) * t394 + t172;
t175 = Icges(5,4) * t202;
t114 = -Icges(5,2) * t201 + Icges(5,6) * t394 + t175;
t462 = t111 + t114;
t533 = t201 * t462 + t202 * t461;
t497 = Icges(5,1) + Icges(6,1);
t518 = Icges(5,5) + Icges(6,5);
t509 = Icges(5,2) + Icges(6,2);
t524 = Icges(5,6) + Icges(6,6);
t105 = Icges(6,5) * t202 - Icges(6,6) * t201 + Icges(6,3) * t394;
t108 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t394;
t496 = t105 + t108;
t530 = t525 * t272;
t529 = t525 * t271;
t523 = Icges(5,3) + Icges(6,3);
t203 = t264 * t272 - t265 * t386;
t392 = t264 * t271;
t204 = t265 * t385 + t392;
t528 = t203 * t462 - t204 * t461;
t472 = t394 * t496 - t533;
t390 = t265 * t268;
t107 = Icges(6,5) * t204 + Icges(6,6) * t203 + Icges(6,3) * t390;
t110 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t390;
t495 = t107 + t110;
t399 = Icges(6,4) * t204;
t113 = Icges(6,2) * t203 + Icges(6,6) * t390 + t399;
t402 = Icges(5,4) * t204;
t116 = Icges(5,2) * t203 + Icges(5,6) * t390 + t402;
t494 = t113 + t116;
t173 = Icges(6,4) * t203;
t119 = Icges(6,1) * t204 + Icges(6,5) * t390 + t173;
t176 = Icges(5,4) * t203;
t122 = Icges(5,1) * t204 + Icges(5,5) * t390 + t176;
t493 = t119 + t122;
t522 = (-t271 * t518 - t272 * t524) * t268;
t521 = (-t272 * t509 - t529) * t268;
t520 = (t271 * t497 + t530) * t268;
t492 = -t523 * t269 + (-t524 * t271 + t518 * t272) * t268;
t516 = -t524 * t269 + (-t509 * t271 + t530) * t268;
t515 = -t518 * t269 + (t497 * t272 - t529) * t268;
t519 = t496 * t390;
t471 = -t494 * t201 + t493 * t202 + t495 * t394;
t470 = t519 + t528;
t469 = t494 * t203 + t493 * t204 + t495 * t390;
t137 = qJD(2) * t201 - qJD(4) * t204;
t138 = -qJD(2) * t202 + qJD(4) * t203;
t356 = qJD(2) * t268;
t335 = t264 * t356;
t504 = t524 * t137 + t518 * t138 - t523 * t335;
t139 = qJD(2) * t203 - qJD(4) * t202;
t281 = t201 * qJD(4);
t140 = qJD(2) * t204 - t281;
t334 = t265 * t356;
t503 = t524 * t139 + t518 * t140 + t523 * t334;
t502 = -t509 * t137 - t525 * t138 + t524 * t335;
t501 = t509 * t139 + t525 * t140 + t524 * t334;
t500 = t525 * t137 + t497 * t138 - t518 * t335;
t499 = t525 * t139 + t497 * t140 + t518 * t334;
t514 = t522 * qJD(4);
t488 = t521 * qJD(4);
t487 = t520 * qJD(4);
t263 = pkin(4) * t272 + pkin(3);
t347 = pkin(4) * t388;
t393 = t264 * t269;
t513 = -t202 * rSges(6,1) + t201 * rSges(6,2) - t263 * t393 + t347;
t512 = t472 * t264;
t353 = qJD(5) * t268;
t239 = t265 * t353;
t254 = -qJD(4) * t269 + qJD(2);
t355 = qJD(4) * t268;
t333 = t264 * t355;
t270 = -qJ(5) - pkin(6);
t422 = pkin(6) + t270;
t423 = pkin(3) - t263;
t374 = (t422 - rSges(6,3)) * t269 + (rSges(6,1) * t272 - rSges(6,2) * t271 - t423) * t268;
t427 = pkin(3) * t269;
t445 = t268 * t422;
t380 = (t427 + t445) * t264 - rSges(6,3) * t394 + t513;
t511 = t254 * t380 + t333 * t374 + t239;
t466 = -t516 * t201 + t515 * t202 + t492 * t394;
t465 = t516 * t203 + t515 * t204 + t492 * t390;
t349 = qJD(2) * qJD(4);
t169 = (qJDD(4) * t264 + t265 * t349) * t268;
t253 = -qJDD(4) * t269 + qJDD(2);
t232 = t265 * pkin(2) + t264 * qJ(3);
t257 = qJD(3) * t265;
t167 = qJD(2) * t232 - t257;
t425 = pkin(6) * t268;
t317 = t425 + t427;
t357 = qJD(2) * t265;
t350 = qJD(2) * qJD(3);
t361 = qJDD(3) * t264 + t265 * t350;
t205 = t317 * t264;
t259 = t265 * qJ(3);
t230 = pkin(2) * t264 - t259;
t367 = -t205 - t230;
t274 = (-t317 * t357 - t167) * qJD(2) + t367 * qJDD(2) + t361;
t246 = pkin(4) * t392;
t309 = rSges(6,1) * t140 + rSges(6,2) * t139;
t444 = t269 * t423;
t418 = rSges(6,3) * t334 + t309 + t264 * t353 - pkin(4) * t281 + (t246 + (-t444 - t445) * t265) * qJD(2);
t417 = rSges(6,2) * t272;
t289 = t268 * (-rSges(6,1) * t271 - t417);
t416 = pkin(4) * qJD(4);
t342 = t271 * t416;
t352 = qJD(5) * t269;
t365 = qJD(4) * t289 - t268 * t342 - t352;
t438 = -qJD(2) * qJD(5) + qJD(4) * t365;
t10 = -t418 * t254 + t380 * t253 + t374 * t169 + (qJDD(5) * t265 + t264 * t438) * t268 + t274;
t508 = -g(1) + t10;
t507 = -t503 * t269 + (t499 * t272 - t501 * t271 + (t271 * t461 - t272 * t462) * qJD(4)) * t268;
t506 = -t504 * t269 + (t500 * t272 + t502 * t271 + (-t271 * t493 - t272 * t494) * qJD(4)) * t268;
t505 = t514 * t269 + (t487 * t272 + t488 * t271 + (t271 * t515 + t272 * t516) * qJD(4)) * t268;
t498 = t492 * t269 + (t271 * t516 - t272 * t515) * t268;
t486 = (t264 * t470 + t265 * t469) * t268;
t485 = (t265 * t471 + t512) * t268;
t484 = (t203 * t518 - t204 * t524) * t265 + (-t201 * t518 - t202 * t524) * t264;
t442 = -t202 * rSges(5,1) + t201 * rSges(5,2);
t126 = rSges(5,3) * t394 - t442;
t200 = -rSges(5,3) * t269 + (rSges(5,1) * t272 - rSges(5,2) * t271) * t268;
t483 = -t126 * t254 + t200 * t333;
t482 = rSges(6,1) + pkin(4);
t481 = t465 * t254;
t479 = t466 * t254;
t459 = (t204 * t509 - t173 - t176 - t493) * t265 + (t202 * t509 + t171 + t174 + t461) * t264;
t170 = (qJDD(4) * t265 - t264 * t349) * t268;
t389 = t265 * t269;
t206 = pkin(3) * t389 + pkin(6) * t390;
t358 = qJD(2) * t264;
t256 = qJD(3) * t264;
t360 = qJ(3) * t357 + t256;
t337 = qJD(2) * (-pkin(2) * t358 + t360) + qJDD(2) * t232 + t264 * t350;
t302 = -qJD(2) ^ 2 * t205 + qJDD(2) * t206 + t337;
t387 = t268 * t270;
t443 = t204 * rSges(6,1) + t203 * rSges(6,2) + rSges(6,3) * t390 + t263 * t389 - t265 * t387 + t246;
t379 = -t206 + t443;
t325 = t269 * t342;
t341 = t272 * t416;
t440 = t138 * rSges(6,1) + t137 * rSges(6,2) + qJD(2) * t347 + t264 * t341 - t265 * t325 + t270 * t335 + t239;
t419 = -rSges(6,3) * t335 + (t425 + t444) * t358 + t440;
t11 = qJDD(5) * t394 + t419 * t254 + t379 * t253 - t374 * t170 + (-t268 * t438 - qJDD(3)) * t265 + t302;
t477 = -g(2) + t11;
t476 = t485 * qJD(4) + t479;
t475 = t486 * qJD(4) + t481;
t474 = (-t265 * t514 + t492 * t358) * t268 + t487 * t204 - t488 * t203 - t515 * t138 - t516 * t137;
t473 = (-t264 * t514 - t492 * t357) * t268 + t487 * t202 + t488 * t201 - t515 * t140 - t516 * t139;
t40 = -t105 * t269 + (-t111 * t271 - t118 * t272) * t268;
t42 = -t108 * t269 + (-t114 * t271 - t121 * t272) * t268;
t468 = t40 + t42;
t41 = -t107 * t269 + (-t113 * t271 + t119 * t272) * t268;
t43 = -t110 * t269 + (-t116 * t271 + t122 * t272) * t268;
t467 = t41 + t43;
t458 = (t497 * t203 - t399 - t402 - t494) * t265 + (-t497 * t201 - t172 - t175 - t462) * t264;
t457 = t484 * t268;
t456 = t507 * t264 + t506 * t265;
t455 = (t499 * t202 - t501 * t201 - t461 * t140 + t462 * t139 + (t503 * t264 + t496 * t357) * t268) * t264 + ((t504 * t264 + t495 * t357) * t268 + t500 * t202 + t502 * t201 + t493 * t140 + t494 * t139) * t265;
t454 = (t500 * t204 - t502 * t203 + t493 * t138 + t494 * t137 + (t504 * t265 - t495 * t358) * t268) * t265 + ((t503 * t265 - t496 * t358) * t268 + t499 * t204 + t501 * t203 - t461 * t138 + t462 * t137) * t264;
t453 = t515 + t521;
t452 = -t516 - t520;
t451 = -t498 * t253 - t505 * t254;
t449 = t470 - t519;
t448 = -m(2) - m(3);
t163 = rSges(4,1) * t389 - rSges(4,2) * t390 + t264 * rSges(4,3);
t435 = t169 / 0.2e1;
t434 = t170 / 0.2e1;
t430 = t264 / 0.2e1;
t429 = -t265 / 0.2e1;
t426 = pkin(4) * t271;
t424 = g(1) * t264;
t226 = (-rSges(5,1) * t271 - rSges(5,2) * t272) * t268;
t214 = qJD(4) * t226;
t311 = -rSges(5,1) * t140 - rSges(5,2) * t139;
t76 = rSges(5,3) * t334 - t311;
t21 = -t126 * t253 + t169 * t200 + t214 * t333 - t254 * t76 + t274;
t415 = t21 * t264;
t130 = t204 * rSges(5,1) + t203 * rSges(5,2) + rSges(5,3) * t390;
t377 = t138 * rSges(5,1) + t137 * rSges(5,2);
t74 = -rSges(5,3) * t335 + t377;
t22 = t130 * t253 - t170 * t200 + t254 * t74 + (-t214 * t355 - qJDD(3)) * t265 + t302;
t414 = t22 * t265;
t297 = qJD(2) * t367 + t256;
t46 = t297 + t483;
t407 = t265 * t46;
t406 = t40 * t169;
t405 = t41 * t170;
t404 = t42 * t169;
t403 = t43 * t170;
t395 = t263 * t269;
t376 = -t202 * rSges(6,2) - t201 * t482;
t375 = -t204 * rSges(6,2) + t203 * t482;
t158 = t163 + t232;
t366 = t206 + t232;
t346 = rSges(4,1) * t393;
t362 = rSges(4,2) * t394 + t265 * rSges(4,3);
t162 = t346 - t362;
t364 = -t230 - t162;
t363 = rSges(4,2) * t335 + rSges(4,3) * t357;
t359 = -qJD(2) * t230 + t256;
t351 = -m(4) - m(5) - m(6);
t336 = -qJD(2) * t205 + t359;
t330 = -rSges(4,1) * t269 - pkin(2);
t329 = -rSges(6,3) * t268 - pkin(2);
t328 = -t355 / 0.2e1;
t327 = t355 / 0.2e1;
t322 = t264 * t328;
t321 = t264 * t327;
t320 = t265 * t328;
t319 = t265 * t327;
t312 = -pkin(2) + (-rSges(6,3) + t270) * t268;
t233 = rSges(3,1) * t265 - rSges(3,2) * t264;
t231 = rSges(3,1) * t264 + rSges(3,2) * t265;
t296 = qJD(2) * t366 - t257;
t47 = -t200 * t265 * t355 + t130 * t254 + t296;
t304 = t264 * t46 - t265 * t47;
t303 = -t264 * t74 + t265 * t76;
t301 = t126 * t265 - t130 * t264;
t282 = -t427 - pkin(2) + (-rSges(5,3) - pkin(6)) * t268;
t156 = rSges(5,1) * t203 - rSges(5,2) * t204;
t154 = -rSges(5,1) * t201 - rSges(5,2) * t202;
t132 = qJD(2) * t158 - t257;
t131 = qJD(2) * t364 + t256;
t56 = -qJDD(3) * t265 + qJDD(2) * t163 + qJD(2) * (-qJD(2) * t346 + t363) + t337;
t55 = t364 * qJDD(2) + (-qJD(2) * t163 - t167) * qJD(2) + t361;
t50 = t301 * t355 + qJD(1);
t31 = t379 * t254 + (-qJD(4) * t265 * t374 + qJD(5) * t264) * t268 + t296;
t30 = t297 + t511;
t27 = -t352 + qJD(1) + (-t264 * t379 - t265 * t380) * t355;
t20 = t126 * t170 - t130 * t169 + t303 * t355 + qJDD(1);
t1 = -qJDD(5) * t269 + qJDD(1) - t380 * t170 - t379 * t169 + (-t264 * t419 + t265 * t418) * t355;
t2 = [m(5) * t20 + m(6) * t1 + (m(4) - t448) * qJDD(1) + (t351 + t448) * g(3); t406 / 0.2e1 + t404 / 0.2e1 + t405 / 0.2e1 + t403 / 0.2e1 - m(3) * (-g(1) * t231 + g(2) * t233) + t466 * t435 + t465 * t434 + (((t469 + t472 + t533) * t265 + t449 * t264) * t355 + t481) * t322 + (m(3) * (t231 ^ 2 + t233 ^ 2) + Icges(4,2) * t269 ^ 2 + (Icges(4,1) * t268 + 0.2e1 * Icges(4,4) * t269) * t268 + Icges(3,3)) * qJDD(2) + (t30 * (t265 * t341 + t257 - t309) + t31 * (t360 + t440) + (t10 * (t329 + t387) + t30 * (t325 - t353)) * t264 + (t30 * (t312 - t395) * t265 + (t30 * (-qJ(3) - t426) + t31 * (t329 - t395)) * t264) * qJD(2) - t312 * t424 - (-t30 + t336 + t511) * t31 + t477 * (t232 + t443) + t508 * (t259 + t513)) * m(6) + (t46 * (t257 + t311) + t47 * (t360 + t377) + (t282 * t407 + (-t46 * qJ(3) + t47 * (-rSges(5,3) * t268 - pkin(2) - t317)) * t264) * qJD(2) - (t336 - t46 + t483) * t47 + (-g(2) + t22) * (t130 + t366) + (-g(1) + t21) * (t264 * t282 + t259 + t442)) * m(5) + (t131 * t257 + t132 * (t360 + t363) + (t131 * (rSges(4,2) * t268 + t330) * t265 + (t131 * (-rSges(4,3) - qJ(3)) + t132 * t330) * t264) * qJD(2) - (-qJD(2) * t162 - t131 + t359) * t132 + (t56 - g(2)) * t158 + (t55 - g(1)) * (t330 * t264 + t259 + t362)) * m(4) + (-t474 + t506) * t319 + (-t473 + t475 + t507) * t321 + t451 + ((t496 * t269 + (t271 * t462 + t272 * t461) * t268 + t468) * t254 + ((t449 - t471 - t528) * t265 - t512) * t355 + t476 - t479) * t320; t351 * (-g(2) * t265 + t424) + 0.2e1 * (t10 * t430 + t11 * t429) * m(6) + 0.2e1 * (t415 / 0.2e1 - t414 / 0.2e1) * m(5) + 0.2e1 * (t429 * t56 + t430 * t55) * m(4); (-t269 * t466 + t485) * t435 + (-t269 * t465 + t486) * t434 + (t498 * t269 + (t264 * t468 + t265 * t467) * t268) * t253 / 0.2e1 - (((-t271 * t453 + t272 * t452) * t254 + ((t271 * t459 + t272 * t458) * t268 - t484 * t269) * qJD(4)) * t268 - t522 * t254 * t269) * t254 / 0.2e1 + (t505 * t269 + ((-t264 * t467 + t265 * t468) * qJD(2) + t456) * t268) * t254 / 0.2e1 - (t456 * t355 + t403 + t404 + t405 + t406 + t451) * t269 / 0.2e1 + (t472 * t169 + t471 * t170 + t466 * t253 - t473 * t254 + t455 * t355) * t394 / 0.2e1 - t475 * t335 / 0.2e1 + ((t201 * t459 + t202 * t458 + t264 * t457) * t355 + (-t201 * t453 + t202 * t452 + t394 * t522) * t254) * t322 + (t473 * t269 + ((-t264 * t471 + t265 * t472) * qJD(2) + t455) * t268) * t321 + ((-t203 * t459 + t458 * t204 + t457 * t265) * t355 + (t203 * t453 + t204 * t452 + t390 * t522) * t254) * t320 + (t474 * t269 + ((-t264 * t469 + t265 * t470) * qJD(2) + t454) * t268) * t319 + ((-t10 * t380 - t11 * t379 + t30 * t418 - t31 * t419) * t269 + ((-t11 * t374 - t31 * t365 - t1 * t380 + t27 * t418 + (-t27 * t379 + t30 * t374) * qJD(2)) * t265 + (t10 * t374 + t30 * t365 - t1 * t379 - t27 * t419 + (t27 * t380 + t31 * t374) * qJD(2)) * t264) * t268 - g(1) * t375 - g(2) * t376 - g(3) * (-t271 * t482 - t417) * t268 - (-t30 * t376 + t31 * t375) * t254 - (t27 * (-t264 * t375 + t265 * t376) + (-t268 * t426 + t289) * (t264 * t30 - t265 * t31)) * t355) * m(6) + ((t126 * t21 - t130 * t22 + t46 * t76 - t47 * t74) * t269 + (t20 * t301 + t50 * (-t126 * t358 - t130 * t357 + t303) + t304 * t214 + (t415 - t414 + (t264 * t47 + t407) * qJD(2)) * t200) * t268 - (-t154 * t46 + t156 * t47) * t254 - (t50 * (t154 * t265 - t156 * t264) + t304 * t226) * t355 - g(1) * t156 - g(2) * t154 - g(3) * t226) * m(5) + (qJD(2) * t476 + t169 * t470 + t170 * t469 + t253 * t465 - t254 * t474 + t355 * t454) * t390 / 0.2e1; ((-t1 + g(3)) * t269 + (t477 * t264 + t508 * t265) * t268) * m(6);];
tau = t2;
