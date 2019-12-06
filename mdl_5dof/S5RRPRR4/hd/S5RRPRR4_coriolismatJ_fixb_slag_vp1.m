% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:50
% EndTime: 2019-12-05 18:31:58
% DurationCPUTime: 4.64s
% Computational Cost: add. (42777->363), mult. (23906->459), div. (0->0), fcn. (21722->10), ass. (0->239)
t329 = qJ(1) + qJ(2);
t322 = pkin(9) + t329;
t319 = sin(t322);
t320 = cos(t322);
t328 = qJ(4) + qJ(5);
t325 = cos(t328);
t317 = Icges(6,4) * t325;
t323 = sin(t328);
t496 = Icges(6,2) * t323 - t317;
t215 = Icges(6,6) * t320 + t319 * t496;
t399 = t320 * t323;
t509 = t215 * t399;
t398 = t320 * t325;
t214 = Icges(6,5) * t398 - Icges(6,6) * t399 + Icges(6,3) * t319;
t405 = t319 * t323;
t277 = Icges(6,5) * t325 - Icges(6,6) * t323;
t412 = t277 * t319;
t379 = t320 * (Icges(6,3) * t320 - t412) + t215 * t405;
t216 = Icges(6,4) * t398 - Icges(6,2) * t399 + Icges(6,6) * t319;
t289 = Icges(6,4) * t399;
t218 = Icges(6,1) * t398 + Icges(6,5) * t319 - t289;
t499 = (t216 * t323 - t218 * t325) * t320;
t508 = t214 * t319 + t379 - t499;
t332 = cos(qJ(4));
t442 = pkin(4) * t332;
t321 = pkin(3) + t442;
t365 = -rSges(6,2) * t405 - rSges(6,3) * t320;
t334 = -pkin(8) - pkin(7);
t395 = t320 * t334;
t435 = rSges(6,1) * t325;
t324 = sin(t329);
t445 = pkin(2) * t324;
t178 = t445 + t395 + (t321 + t435) * t319 + t365;
t441 = sin(qJ(1)) * pkin(1);
t172 = t178 + t441;
t507 = -t172 + t178;
t330 = sin(qJ(4));
t397 = t320 * t330;
t300 = rSges(5,2) * t397;
t436 = rSges(5,1) * t332;
t360 = pkin(3) + t436;
t326 = cos(t329);
t444 = pkin(2) * t326;
t190 = -t444 + t300 - t360 * t320 + (-rSges(5,3) - pkin(7)) * t319;
t307 = rSges(5,1) * t330 + rSges(5,2) * t332;
t258 = t307 * t320;
t314 = t320 * pkin(7);
t403 = t319 * t330;
t364 = -rSges(5,2) * t403 - rSges(5,3) * t320;
t189 = t319 * t360 - t314 + t364 + t445;
t257 = t307 * t319;
t419 = t189 * t257;
t102 = t190 * t258 - t419;
t290 = t320 * t321;
t354 = -rSges(6,1) * t398 + rSges(6,2) * t399;
t179 = -t444 - t290 + (-rSges(6,3) + t334) * t319 + t354;
t282 = rSges(6,1) * t323 + rSges(6,2) * t325;
t443 = pkin(4) * t330;
t352 = (t282 + t443) * t320;
t136 = t352 * t179;
t404 = t319 * t325;
t247 = rSges(6,1) * t405 + rSges(6,2) * t404;
t301 = pkin(4) * t403;
t211 = t301 + t247;
t91 = -t178 * t211 + t136;
t506 = -m(5) * t102 - m(6) * t91;
t479 = m(5) / 0.2e1;
t478 = m(6) / 0.2e1;
t505 = -t319 / 0.2e1;
t455 = t319 / 0.2e1;
t453 = t320 / 0.2e1;
t440 = cos(qJ(1)) * pkin(1);
t451 = m(3) * (-t440 * (rSges(3,1) * t324 + rSges(3,2) * t326) - (-rSges(3,1) * t326 + rSges(3,2) * t324) * t441);
t449 = m(4) * (-t440 * (rSges(4,1) * t319 + rSges(4,2) * t320 + t445) - (-rSges(4,1) * t320 + rSges(4,2) * t319 - t444) * t441);
t315 = t319 ^ 2;
t316 = t320 ^ 2;
t361 = t315 + t316;
t454 = -t320 / 0.2e1;
t504 = t453 + t454;
t248 = t282 * t320;
t93 = -t178 * t247 + t179 * t248;
t503 = m(6) * t93;
t174 = -t247 * t319 - t248 * t320;
t410 = t282 * t319;
t226 = t301 + t410;
t284 = -rSges(6,2) * t323 + t435;
t408 = t284 * t320;
t409 = t284 * t319;
t347 = Icges(6,5) * t323 + Icges(6,6) * t325;
t241 = t319 * t347;
t242 = t347 * t320;
t288 = Icges(6,4) * t405;
t217 = -Icges(6,1) * t404 + Icges(6,5) * t320 + t288;
t373 = Icges(6,2) * t404 + t217 + t288;
t498 = Icges(6,1) * t323 + t317;
t375 = -t319 * t498 + t215;
t490 = -t323 * t373 - t325 * t375;
t372 = -Icges(6,2) * t398 + t218 - t289;
t374 = t320 * t498 + t216;
t491 = t323 * t372 + t325 * t374;
t439 = (t316 * t241 + (t491 * t319 + (-t242 - t490) * t320) * t319) * t453 + (-t315 * t242 + (t490 * t320 + (t241 - t491) * t319) * t320) * t455;
t202 = t320 * (rSges(6,3) * t319 - t354);
t219 = -rSges(6,1) * t404 - t365;
t96 = -t320 * (t320 * pkin(3) - t290 + (pkin(7) + t334) * t319) + t202 + (t395 - t219 + t314 + (-pkin(3) + t321) * t319) * t319;
t15 = t439 + m(6) * (t174 * t96 + t226 * t409 + t352 * t408);
t502 = t15 * qJD(5);
t309 = -rSges(5,2) * t330 + t436;
t501 = t309 * t479;
t396 = t320 * t332;
t233 = Icges(5,4) * t396 - Icges(5,2) * t397 + Icges(5,6) * t319;
t298 = Icges(5,4) * t397;
t235 = Icges(5,1) * t396 + Icges(5,5) * t319 - t298;
t500 = (t233 * t330 - t235 * t332) * t320;
t173 = t179 - t440;
t134 = t352 * t173;
t75 = -t172 * t179 + t173 * t178;
t185 = t189 + t441;
t186 = t190 - t440;
t88 = -t185 * t190 + t186 * t189;
t327 = Icges(5,4) * t332;
t497 = Icges(5,1) * t330 + t327;
t495 = Icges(5,2) * t330 - t327;
t494 = qJD(1) + qJD(2);
t420 = t185 * t257;
t438 = (t226 * t507 + t134 - t136) * t478 + (-t420 + t419 + (t186 - t190) * t258) * t479;
t90 = -t172 * t211 + t134;
t99 = t186 * t258 - t420;
t493 = (t91 + t90) * t478 + (t102 + t99) * t479;
t426 = Icges(6,4) * t323;
t278 = Icges(6,2) * t325 + t426;
t281 = Icges(6,1) * t325 - t426;
t492 = (-t496 + t498) * t323 + (t278 - t281) * t325;
t427 = Icges(5,4) * t330;
t303 = Icges(5,2) * t332 + t427;
t306 = Icges(5,1) * t332 - t427;
t489 = (-t495 + t497) * t330 + (t303 - t306) * t332;
t368 = -Icges(5,2) * t396 + t235 - t298;
t370 = t320 * t497 + t233;
t488 = t330 * t368 + t332 * t370;
t297 = Icges(5,4) * t403;
t402 = t319 * t332;
t234 = -Icges(5,1) * t402 + Icges(5,5) * t320 + t297;
t369 = Icges(5,2) * t402 + t234 + t297;
t232 = Icges(5,6) * t320 + t319 * t495;
t371 = -t319 * t497 + t232;
t487 = -t330 * t369 - t332 * t371;
t486 = (t497 / 0.2e1 - t495 / 0.2e1) * t332 + (t306 / 0.2e1 - t303 / 0.2e1) * t330;
t356 = (-t496 / 0.2e1 + t498 / 0.2e1) * t325 + (-t278 / 0.2e1 + t281 / 0.2e1) * t323;
t106 = -t217 * t404 + t379;
t107 = t320 * t214 + t216 * t405 - t218 * t404;
t358 = -t217 * t325 - t214;
t416 = t215 * t323;
t11 = ((t499 + t508) * t320 + ((t358 - t416) * t320 + t107 + t509) * t319) * t455 + (t106 * t320 + t107 * t319) * t505 + ((t107 + (-t214 + t416) * t320 - t509) * t320 + (t319 * t358 - t106 + t508) * t319) * t453;
t485 = 4 * qJD(1);
t482 = 2 * qJD(4);
t475 = m(5) * t88;
t473 = m(5) * t99;
t92 = -t172 * t247 + t173 * t248;
t468 = m(6) * (t93 + t92);
t466 = m(6) * ((t173 - t179) * t248 + t507 * t410);
t150 = t173 * t409;
t380 = t226 * t248 - t247 * t352;
t465 = m(6) * (t172 * t408 + t150 + t380);
t183 = t352 * t410;
t417 = t211 * t282;
t423 = t172 * t284;
t464 = m(6) * (t150 + t183 + (-t417 + t423) * t320);
t152 = t179 * t409;
t463 = m(6) * (t178 * t408 + t152 + t380);
t422 = t178 * t284;
t462 = m(6) * (t152 + t183 + (-t417 + t422) * t320);
t461 = m(6) * t75;
t459 = m(6) * t90;
t457 = m(6) * t92;
t231 = Icges(5,5) * t396 - Icges(5,6) * t397 + Icges(5,3) * t319;
t119 = t320 * t231 + t233 * t403 - t235 * t402;
t302 = Icges(5,5) * t332 - Icges(5,6) * t330;
t406 = t302 * t319;
t230 = Icges(5,3) * t320 - t406;
t376 = t230 * t319 + t234 * t396;
t120 = -t232 * t397 + t376;
t121 = t231 * t319 - t500;
t357 = -t234 * t332 - t231;
t377 = t230 * t320 + t232 * t403;
t415 = t232 * t330;
t28 = (t121 + t377 + t500) * t320 + (-t120 + (t357 - t415) * t320 + t119 + t376) * t319;
t118 = -t234 * t402 + t377;
t29 = (t119 + (-t231 + t415) * t320 - t376) * t320 + (t319 * t357 - t118 + t377) * t319;
t84 = t118 * t320 + t119 * t319;
t85 = t120 * t320 + t121 * t319;
t2 = (t29 / 0.2e1 + t85 / 0.2e1) * t320 + (-t84 / 0.2e1 + t28 / 0.2e1) * t319 + t11;
t452 = t2 * qJD(4);
t446 = m(6) * t174;
t430 = t11 * qJD(5);
t414 = t247 * t282;
t381 = (-t211 + t226) * t352;
t359 = t284 + t442;
t351 = t468 / 0.2e1 + t356;
t348 = Icges(5,5) * t330 + Icges(5,6) * t332;
t251 = t319 * t348;
t341 = -t11 + (-t320 * t492 - t323 * t374 + t325 * t372 + t412) * t455 + (t277 * t320 + t319 * t492 - t323 * t375 + t325 * t373) * t453;
t340 = -t356 + t504 * (t325 * t216 + t323 * t218);
t339 = t356 + t486;
t337 = t339 + t493;
t336 = t340 - t486 + t504 * (t233 * t332 + t330 * t235);
t335 = (t28 * t505 + t341 + (t85 + t29) * t454 + (t302 * t320 + t319 * t489 - t330 * t371 + t332 * t369) * t453 + (-t320 * t489 - t330 * t370 + t332 * t368 + t406 + t84) * t455) * qJD(4);
t252 = t348 * t320;
t229 = t359 * t320;
t227 = t359 * t319;
t195 = t248 * t410;
t180 = -t257 * t319 - t258 * t320;
t166 = qJD(5) * t446;
t154 = -t319 * t219 + t202;
t148 = -t361 * t443 + t174;
t74 = t356 + t503;
t73 = t356 + t457;
t67 = t462 / 0.2e1;
t65 = t463 / 0.2e1;
t58 = t464 / 0.2e1;
t57 = t465 / 0.2e1;
t53 = t466 / 0.2e1;
t37 = t339 - t506;
t36 = t339 + t459 + t473;
t30 = t449 + t451 + t461 + t475;
t19 = -t466 / 0.2e1 + t351;
t18 = t53 + t351;
t17 = m(6) * (t282 * t284 * t361 + t154 * t174) + t439;
t16 = t17 * qJD(5);
t14 = t53 - t468 / 0.2e1 + t340;
t13 = t337 + t438;
t12 = t337 - t438;
t9 = t336 + t438 - t493;
t8 = t65 - t462 / 0.2e1 + t11;
t7 = t67 - t463 / 0.2e1 + t11;
t6 = t57 - t464 / 0.2e1 + t11;
t5 = t58 - t465 / 0.2e1 + t11;
t4 = t65 + t67 + t341;
t3 = t57 + t58 + t341;
t1 = [qJD(2) * t30 + qJD(4) * t36 + qJD(5) * t73, t30 * qJD(1) + t13 * qJD(4) + t18 * qJD(5) + 0.2e1 * (t449 / 0.2e1 + t451 / 0.2e1 + t75 * t478 + t88 * t479) * qJD(2), 0, t36 * qJD(1) + t13 * qJD(2) + t3 * qJD(5) + ((t172 * t229 + t173 * t227 + t381) * t478 + (t185 * t320 + t186 * t319) * t501) * t482 + t335, t73 * qJD(1) + t18 * qJD(2) + t3 * qJD(4) + ((t150 + t195 + (-t414 + t423) * t320) * m(6) + t341) * qJD(5); t12 * qJD(4) + t19 * qJD(5) + (-t451 / 0.4e1 - t449 / 0.4e1 - t475 / 0.4e1 - t461 / 0.4e1) * t485, qJD(4) * t37 + qJD(5) * t74, 0, t12 * qJD(1) + t37 * qJD(2) + t4 * qJD(5) + ((t178 * t229 + t179 * t227 + t381) * t478 + (t189 * t320 + t190 * t319) * t501) * t482 + t335, t19 * qJD(1) + t74 * qJD(2) + t4 * qJD(4) + ((t152 + t195 + (-t414 + t422) * t320) * m(6) + t341) * qJD(5); 0, 0, 0, (t148 * t478 + t180 * t479) * t482 + t166, qJD(4) * t446 + t166; t9 * qJD(2) + t6 * qJD(5) + (-t473 / 0.4e1 - t459 / 0.4e1) * t485 + t452 + t336 * qJD(1), t9 * qJD(1) + t8 * qJD(5) + t452 + (t336 + t506) * qJD(2), 0, (m(5) * (t361 * t309 * t307 + (t320 * (rSges(5,1) * t396 + rSges(5,3) * t319 - t300) - t319 * (-rSges(5,1) * t402 - t364)) * t180) + (t316 * t251 + (t488 * t319 + (-t252 - t487) * t320) * t319) * t453 + (-t315 * t252 + (t487 * t320 + (t251 - t488) * t319) * t320) * t455 + m(6) * (t148 * t96 + t226 * t227 + t229 * t352) + t439) * qJD(4) + t502 + t494 * t2, t6 * qJD(1) + t8 * qJD(2) + t15 * qJD(4) + t502; (t340 - t457) * qJD(1) + t14 * qJD(2) + t5 * qJD(4) + t430, t14 * qJD(1) + (t340 - t503) * qJD(2) + t7 * qJD(4) + t430, 0, t5 * qJD(1) + t7 * qJD(2) + ((t148 * t154 + (t227 * t319 + t229 * t320) * t282) * m(6) + t439) * qJD(4) + t16, qJD(4) * t17 + t11 * t494 + t16;];
Cq = t1;
