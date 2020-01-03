% Calculate vector of inverse dynamics joint torques for
% S5RPRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:41
% EndTime: 2019-12-31 18:40:56
% DurationCPUTime: 12.03s
% Computational Cost: add. (13137->507), mult. (10055->610), div. (0->0), fcn. (7736->8), ass. (0->286)
t279 = qJ(1) + pkin(8);
t269 = qJ(3) + t279;
t263 = sin(t269);
t280 = sin(qJ(4));
t282 = cos(qJ(4));
t229 = Icges(5,5) * t282 - Icges(5,6) * t280;
t264 = cos(t269);
t313 = t229 * t264;
t131 = Icges(5,3) * t263 + t313;
t231 = Icges(6,4) * t282 + Icges(6,6) * t280;
t314 = t231 * t264;
t133 = Icges(6,2) * t263 + t314;
t538 = t131 + t133;
t271 = Icges(6,5) * t280;
t338 = Icges(6,1) * t282 + t271;
t136 = -Icges(6,4) * t264 + t263 * t338;
t413 = t263 * t280;
t218 = Icges(5,4) * t413;
t412 = t263 * t282;
t138 = Icges(5,1) * t412 - Icges(5,5) * t264 - t218;
t537 = -t136 - t138;
t316 = t338 * t264;
t137 = Icges(6,4) * t263 + t316;
t431 = Icges(5,4) * t280;
t237 = Icges(5,1) * t282 - t431;
t317 = t237 * t264;
t139 = Icges(5,5) * t263 + t317;
t530 = t137 + t139;
t232 = Icges(5,2) * t282 + t431;
t336 = Icges(6,3) * t282 - t271;
t536 = t232 + t336;
t430 = Icges(6,5) * t282;
t234 = Icges(6,1) * t280 - t430;
t272 = Icges(5,4) * t282;
t236 = Icges(5,1) * t280 + t272;
t535 = t234 + t236;
t227 = Icges(6,3) * t280 + t430;
t337 = -Icges(5,2) * t280 + t272;
t534 = t227 - t337;
t533 = t237 + t338;
t409 = t264 * t282;
t217 = Icges(6,5) * t409;
t410 = t264 * t280;
t129 = Icges(6,6) * t263 + Icges(6,3) * t410 + t217;
t532 = t129 * t410 + t538 * t263 + t530 * t409;
t132 = -Icges(6,2) * t264 + t231 * t263;
t120 = t263 * t132;
t128 = -Icges(6,6) * t264 + t227 * t263;
t130 = Icges(5,5) * t412 - Icges(5,6) * t413 - Icges(5,3) * t264;
t531 = -t128 * t410 - t263 * t130 + t537 * t409 - t120;
t228 = Icges(5,5) * t280 + Icges(5,6) * t282;
t230 = Icges(6,4) * t280 - Icges(6,6) * t282;
t529 = t228 + t230;
t278 = qJD(1) + qJD(3);
t528 = (-Icges(5,6) + Icges(6,6)) * t278 + t536 * qJD(4);
t527 = (Icges(6,4) + Icges(5,5)) * t278 - t535 * qJD(4);
t523 = -t536 * t280 + t535 * t282;
t134 = Icges(5,4) * t412 - Icges(5,2) * t413 - Icges(5,6) * t264;
t422 = t134 * t280;
t331 = -t138 * t282 + t422;
t423 = t132 * t264;
t334 = t128 * t280 + t136 * t282;
t472 = t263 * t334;
t50 = -t423 + t472;
t479 = -t130 * t264 - t263 * t331 + t50;
t477 = -t134 * t410 - t531;
t315 = t337 * t264;
t135 = Icges(5,6) * t263 + t315;
t476 = -t135 * t410 + t532;
t350 = -t129 * t413 + t133 * t264 - t137 * t412;
t111 = t139 * t412;
t354 = t131 * t264 - t111;
t53 = -t135 * t413 - t354;
t478 = t53 - t350;
t526 = t534 * qJD(4);
t525 = t533 * qJD(4);
t524 = t229 + t231;
t522 = -t535 * t280 - t536 * t282;
t414 = t263 * t278;
t521 = t528 * t264 - t534 * t414;
t411 = t264 * t278;
t520 = -t227 * t411 - t528 * t263 + t278 * t315;
t519 = t527 * t264 - t533 * t414;
t518 = (t316 + t317) * t278 + t527 * t263;
t517 = (t129 - t135) * t282 - t530 * t280;
t475 = (t128 - t134) * t282 + t537 * t280;
t503 = rSges(6,3) + qJ(5);
t468 = t503 * t412;
t416 = t230 * t264;
t419 = t228 * t264;
t489 = t523 * t263 - t416 - t419;
t417 = t230 * t263;
t420 = t228 * t263;
t488 = t523 * t264 + t417 + t420;
t516 = (-Icges(6,2) - Icges(5,3)) * t278 + t529 * qJD(4);
t421 = t135 * t280;
t515 = t129 * t280 + t530 * t282 - t421;
t514 = -t331 + t334;
t270 = t280 * qJ(5);
t347 = t282 * rSges(6,1) + t280 * rSges(6,3);
t513 = t282 * pkin(4) + t270 + t347;
t496 = rSges(6,1) + pkin(4);
t186 = t264 * pkin(3) + t263 * pkin(7);
t150 = t186 * t278;
t378 = qJD(4) * t278;
t174 = -qJDD(4) * t264 + t263 * t378;
t277 = qJDD(1) + qJDD(3);
t267 = sin(t279);
t268 = cos(t279);
t284 = qJD(1) ^ 2;
t281 = sin(qJ(1));
t283 = cos(qJ(1));
t312 = (-qJDD(1) * t281 - t283 * t284) * pkin(1);
t294 = (-qJDD(1) * t267 - t268 * t284) * pkin(2) + t312;
t374 = qJD(5) * t282;
t391 = -qJD(4) * t513 + t374;
t304 = qJDD(5) * t280 + (t374 + t391) * qJD(4);
t375 = qJD(5) * t280;
t362 = t263 * t375;
t259 = t264 * pkin(7);
t185 = pkin(3) * t263 - t259;
t256 = t264 * rSges(6,2);
t396 = t513 * t263 - t256;
t371 = -t185 - t396;
t383 = t496 * t280 - t503 * t282;
t376 = qJD(4) * t282;
t366 = t263 * t376;
t377 = qJD(4) * t280;
t367 = t263 * t377;
t469 = t367 * t496 - t362;
t407 = t278 * t280;
t485 = t264 * t407 + t366;
t487 = t263 * rSges(6,2) + pkin(4) * t409;
t441 = rSges(6,3) * t366 + t485 * qJ(5) - t469 + (t264 * t347 + t487) * t278;
t14 = t383 * t174 + t371 * t277 + (-t150 - t362 - t441) * t278 + t304 * t264 + t294;
t482 = t14 - g(1);
t173 = qJDD(4) * t263 + t264 * t378;
t211 = t264 * t375;
t210 = pkin(7) * t411;
t262 = pkin(2) * t268;
t276 = t283 * pkin(1);
t265 = qJDD(1) * t276;
t444 = pkin(1) * t281;
t352 = -pkin(2) * t267 - t444;
t307 = qJDD(1) * t262 + t284 * t352 + t265;
t296 = t278 * (-pkin(3) * t414 + t210) + t277 * t186 + t307;
t394 = rSges(6,1) * t409 + rSges(6,3) * t410 + t264 * t270 + t487;
t364 = t264 * t377;
t311 = -t278 * t412 - t364;
t373 = t263 * t407;
t363 = t264 * t376;
t464 = rSges(6,2) * t411 + t363 * t503 + t211;
t442 = t496 * t311 - t373 * t503 + t464;
t15 = t394 * t277 - t383 * t173 + (t211 + t442) * t278 + t304 * t263 + t296;
t512 = t15 - g(2);
t440 = rSges(5,1) * t282;
t244 = -rSges(5,2) * t280 + t440;
t198 = t244 * qJD(4);
t240 = rSges(5,1) * t280 + rSges(5,2) * t282;
t379 = qJD(4) * t264;
t388 = rSges(5,2) * t413 + t264 * rSges(5,3);
t141 = rSges(5,1) * t412 - t388;
t395 = -t141 - t185;
t369 = rSges(5,1) * t367 + t485 * rSges(5,2);
t466 = rSges(5,1) * t409 + t263 * rSges(5,3);
t96 = t278 * t466 - t369;
t26 = -t198 * t379 + t174 * t240 + (-t150 - t96) * t278 + t395 * t277 + t294;
t511 = t26 - g(1);
t143 = -rSges(5,2) * t410 + t466;
t380 = qJD(4) * t263;
t318 = rSges(5,3) * t411 + (-t363 + t373) * rSges(5,2);
t94 = rSges(5,1) * t311 + t318;
t27 = t143 * t277 - t173 * t240 - t198 * t380 + t278 * t94 + t296;
t510 = t27 - g(2);
t207 = rSges(4,2) * t414;
t149 = rSges(4,1) * t411 - t207;
t183 = rSges(4,1) * t263 + rSges(4,2) * t264;
t509 = -t149 * t278 - t183 * t277 - g(1) + t294;
t257 = t264 * rSges(4,1);
t184 = -rSges(4,2) * t263 + t257;
t408 = t278 * t183;
t508 = t184 * t277 - t278 * t408 - g(2) + t307;
t507 = t522 * qJD(4) + t529 * t278 + t526 * t280 + t525 * t282;
t506 = t476 * t263 - t477 * t264;
t505 = t478 * t263 - t479 * t264;
t504 = -t524 * qJD(4) + t523 * t278;
t502 = t488 * t278;
t501 = t517 * qJD(4) + t538 * t278 + t521 * t280 + t519 * t282;
t500 = -t518 * t282 + t520 * t280 + (-t130 - t132) * t278 - t475 * qJD(4);
t499 = t489 * t278;
t498 = (-t314 - t313 + t514) * t278 + t516 * t263;
t497 = t516 * t264 + t515 * t278 + t524 * t414;
t495 = t505 * qJD(4) + t499;
t494 = t506 * qJD(4) + t502;
t493 = t514 * qJD(4) + t518 * t280 + t520 * t282;
t492 = t515 * qJD(4) + t519 * t280 - t521 * t282;
t491 = -t504 * t263 + t507 * t264;
t490 = t507 * t263 + t504 * t264;
t486 = t423 + t532;
t484 = t263 * t498 + t264 * t500;
t483 = -t263 * t497 + t264 * t501;
t481 = t263 * t500 - t264 * t498;
t480 = t263 * t501 + t264 * t497;
t322 = t352 * qJD(1);
t365 = t240 * t379;
t297 = t322 - t365;
t62 = t278 * t395 + t297;
t473 = t278 * t62;
t125 = t278 * t141;
t176 = t278 * t185;
t471 = -t125 - t176;
t470 = t186 + t394;
t188 = t268 * rSges(3,1) - rSges(3,2) * t267;
t178 = t188 + t276;
t467 = t503 * t409;
t465 = t262 + t276;
t462 = -t396 * t278 - t176;
t398 = -Icges(5,2) * t412 + t138 - t218;
t402 = t236 * t263 + t134;
t455 = -t280 * t398 - t282 * t402;
t400 = -t336 * t263 + t136;
t404 = -t234 * t263 + t128;
t454 = -t280 * t400 + t282 * t404;
t453 = m(3) + m(4);
t452 = t173 / 0.2e1;
t451 = t174 / 0.2e1;
t443 = g(2) * t263;
t320 = -t379 * t383 + t211;
t295 = t322 + t320;
t48 = t278 * t371 + t295;
t438 = t264 * t48;
t437 = t264 * t62;
t436 = t278 * t48;
t47 = -t374 + qJD(2) + (t396 * t263 + t394 * t264) * qJD(4);
t435 = t47 * t280;
t321 = t465 * qJD(1);
t123 = t184 * t278 + t321;
t425 = t123 * t183;
t418 = t229 * t278;
t415 = t231 * t278;
t403 = -Icges(6,1) * t410 + t129 + t217;
t401 = -t236 * t264 - t135;
t399 = -t336 * t264 + t137;
t397 = -t232 * t264 + t139;
t108 = t143 + t186;
t393 = -t496 * t413 + t468;
t392 = -t496 * t410 + t467;
t387 = -t336 + t338;
t386 = t227 - t234;
t385 = -t232 + t237;
t384 = t236 + t337;
t368 = t264 * t496;
t359 = -pkin(3) - t440;
t358 = -t380 / 0.2e1;
t357 = t380 / 0.2e1;
t356 = -t379 / 0.2e1;
t355 = t379 / 0.2e1;
t353 = -t130 + t421;
t245 = rSges(2,1) * t283 - rSges(2,2) * t281;
t241 = rSges(2,1) * t281 + rSges(2,2) * t283;
t187 = rSges(3,1) * t267 + rSges(3,2) * t268;
t323 = -t383 * t380 + t362;
t49 = t278 * t470 + t321 + t323;
t345 = t263 * t49 + t438;
t181 = t240 * t380;
t63 = t108 * t278 - t181 + t321;
t340 = -t263 * t63 - t437;
t329 = t141 * t263 + t143 * t264;
t319 = t345 * t282;
t310 = -t280 * t399 + t282 * t403;
t309 = -t280 * t397 + t282 * t401;
t107 = t263 * t359 + t259 + t388;
t308 = -t280 * t503 - t282 * t496 - pkin(3);
t122 = t322 - t408;
t306 = (t280 * t386 + t282 * t387) * t278;
t305 = (-t280 * t384 + t282 * t385) * t278;
t77 = t263 * t308 + t256 + t259;
t287 = (((t53 - t111 + (t131 + t422) * t264 + t531) * t264 + (t50 - t472 + t486) * t263) * qJD(4) + t502) * t355 + (t523 * qJD(4) + t525 * t280 - t526 * t282) * t278 + (Icges(4,3) - t522) * t277 + (t488 - t517) * t452 + (t489 - t475) * t451 + (t491 + t492) * t357 + (((t264 * t353 + t476 - t486) * t264 + (t263 * t353 - t120 + t350 + t354 + t477) * t263) * qJD(4) + t495 - t499) * t358 + (t490 + t493 + t494) * t356;
t286 = t62 * t369 + t63 * (-rSges(5,1) * t364 + t210 + t318) + (t359 * t437 + (t62 * (-rSges(5,3) - pkin(7)) + t63 * t359) * t263) * t278;
t285 = t48 * t469 + t49 * (t210 + t464) + (t308 * t438 + (t48 * (-rSges(6,2) - pkin(7)) + t49 * (-pkin(3) - t513)) * t263) * t278 + (-t280 * t368 * t49 - t48 * t468) * qJD(4);
t171 = t240 * t264;
t167 = t240 * t263;
t68 = qJD(4) * t329 + qJD(2);
t25 = t141 * t173 - t143 * t174 + qJDD(2) + (t263 * t96 + t264 * t94) * qJD(4);
t5 = -qJDD(5) * t282 + qJDD(2) - t394 * t174 + t396 * t173 + (t263 * t441 + t264 * t442 + t375) * qJD(4);
t1 = [t287 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t122 * t207 + (-t122 * t257 - t425) * t278 + (-t122 * t465 + t123 * t352) * qJD(1) + t508 * (t184 + t465) + t509 * (-t183 + t352)) * m(4) + ((qJDD(1) * t188 - g(2) + t265) * t178 + (-qJDD(1) * t187 + t312 + (-0.2e1 * t188 + 0.2e1 * t178 - t276) * t284 - g(1)) * (-t187 - t444)) * m(3) + ((t241 ^ 2 + t245 ^ 2) * qJDD(1) + g(1) * t241 - g(2) * t245) * m(2) + (-(-t48 + t295 + t462) * t49 + (t352 * t49 - t465 * t48) * qJD(1) + t285 + t512 * (t470 + t465) + t482 * (t352 + t77)) * m(6) + (-(-t62 + t297 + t471) * t63 + (t352 * t63 - t465 * t62) * qJD(1) + t286 + t510 * (t108 + t465) + t511 * (t107 + t352)) * m(5); t453 * qJDD(2) + m(5) * t25 + m(6) * t5 + (-m(5) - m(6) - t453) * g(3); t287 + (t285 + t48 * t323 - t49 * (t320 + t462) + t482 * t77 + (t436 + t512) * t470) * m(6) + (t286 - t62 * t181 - t63 * (-t365 + t471) + (t473 + t510) * t108 + t511 * t107) * m(5) + (-t122 * t149 - t123 * t408 + t425 * t278 + (t122 * t278 + t508) * t184 - t509 * t183) * m(4); t506 * t452 + t505 * t451 + (t491 * t278 + t488 * t277 + t477 * t174 + t476 * t173 + (t483 * t263 + t484 * t264) * qJD(4)) * t263 / 0.2e1 - (t490 * t278 + t489 * t277 + t479 * t174 + t478 * t173 + (t480 * t263 + t481 * t264) * qJD(4)) * t264 / 0.2e1 + (-t263 * t517 + t264 * t475) * t277 / 0.2e1 - (((t384 - t386) * t282 + (t385 + t387) * t280) * t278 + (((-t398 - t400) * t264 + (t397 + t399) * t263) * t282 + ((t402 - t404) * t264 + (t401 + t403) * t263) * t280) * qJD(4)) * t278 / 0.2e1 + ((-t278 * t517 - t493) * t264 + (-t278 * t475 + t492) * t263) * t278 / 0.2e1 + t495 * t414 / 0.2e1 + t494 * t411 / 0.2e1 + ((-t380 * t416 + t415) * t263 + (t306 + (-t454 * t264 + (t417 + t310) * t263) * qJD(4)) * t264 + (-t380 * t419 + t418) * t263 + (t305 + (-t455 * t264 + (t420 + t309) * t263) * qJD(4)) * t264) * t358 + ((t278 * t476 + t484) * t264 + (t278 * t477 + t483) * t263) * t357 + ((t278 * t478 + t481) * t264 + (t278 * t479 + t480) * t263) * t356 + ((-t379 * t417 - t415) * t264 + (t306 + (t310 * t263 + (t416 - t454) * t264) * qJD(4)) * t263 + (-t379 * t420 - t418) * t264 + (t305 + (t309 * t263 + (t419 - t455) * t264) * qJD(4)) * t263) * t355 + (-g(1) * t467 - g(2) * t468 - g(3) * t513 - (-g(1) * t368 - t443 * t496) * t280 + (-t14 * t383 + t48 * t391 + t5 * t394 + t47 * t442 + (-t383 * t49 + t396 * t47) * t278) * t264 + (-t15 * t383 + t49 * t391 + t5 * t396 + t47 * t441 + (t383 * t48 - t394 * t47) * t278) * t263 - (t319 + t435) * qJD(5) - (t392 * t49 - t393 * t48) * t278 - ((t392 * t47 - t48 * t513) * t264 + (t393 * t47 - t49 * t513) * t263) * qJD(4)) * m(6) + (t25 * t329 + t68 * ((t94 + t125) * t264 + (-t143 * t278 + t96) * t263) + t340 * t198 + ((-t278 * t63 - t26) * t264 + (-t27 + t473) * t263) * t240 - (t167 * t62 - t171 * t63) * t278 - (t68 * (-t167 * t263 - t171 * t264) + t340 * t244) * qJD(4) + g(1) * t171 + g(2) * t167 - g(3) * t244) * m(5); (-(-t263 * t48 + t264 * t49) * t407 - (t319 + (t263 ^ 2 + t264 ^ 2) * t435) * qJD(4) + (qJD(4) * t345 + g(3) - t5) * t282 + (qJD(4) * t47 + (t15 - t436) * t263 - t443 + (t278 * t49 + t482) * t264) * t280) * m(6);];
tau = t1;
