% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:07
% EndTime: 2019-12-05 18:16:19
% DurationCPUTime: 4.57s
% Computational Cost: add. (42067->354), mult. (23472->450), div. (0->0), fcn. (21350->10), ass. (0->234)
t319 = qJ(1) + pkin(9);
t315 = qJ(3) + t319;
t310 = sin(t315);
t311 = cos(t315);
t320 = qJ(4) + qJ(5);
t317 = cos(t320);
t309 = Icges(6,4) * t317;
t316 = sin(t320);
t481 = Icges(6,2) * t316 - t309;
t211 = Icges(6,6) * t311 + t481 * t310;
t392 = t311 * t316;
t494 = t211 * t392;
t391 = t311 * t317;
t210 = Icges(6,5) * t391 - Icges(6,6) * t392 + Icges(6,3) * t310;
t398 = t310 * t316;
t271 = Icges(6,5) * t317 - Icges(6,6) * t316;
t404 = t271 * t310;
t372 = t311 * (Icges(6,3) * t311 - t404) + t211 * t398;
t212 = Icges(6,4) * t391 - Icges(6,2) * t392 + Icges(6,6) * t310;
t281 = Icges(6,4) * t392;
t214 = Icges(6,1) * t391 + Icges(6,5) * t310 - t281;
t484 = (t212 * t316 - t214 * t317) * t311;
t493 = t210 * t310 + t372 - t484;
t323 = cos(qJ(4));
t432 = pkin(4) * t323;
t312 = pkin(3) + t432;
t358 = -rSges(6,2) * t398 - t311 * rSges(6,3);
t325 = -pkin(8) - pkin(7);
t388 = t311 * t325;
t427 = rSges(6,1) * t317;
t184 = t388 + (t312 + t427) * t310 + t358;
t346 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t319);
t171 = t184 + t346;
t492 = -t171 + t184;
t321 = sin(qJ(4));
t390 = t311 * t321;
t292 = rSges(5,2) * t390;
t428 = rSges(5,1) * t323;
t353 = pkin(3) + t428;
t192 = t292 - t353 * t311 + (-rSges(5,3) - pkin(7)) * t310;
t299 = rSges(5,1) * t321 + rSges(5,2) * t323;
t250 = t299 * t311;
t306 = t311 * pkin(7);
t396 = t310 * t321;
t357 = -rSges(5,2) * t396 - t311 * rSges(5,3);
t191 = t353 * t310 - t306 + t357;
t249 = t299 * t310;
t411 = t191 * t249;
t103 = t192 * t250 - t411;
t282 = t311 * t312;
t347 = -rSges(6,1) * t391 + rSges(6,2) * t392;
t185 = -t282 + (-rSges(6,3) + t325) * t310 + t347;
t276 = rSges(6,1) * t316 + rSges(6,2) * t317;
t433 = pkin(4) * t321;
t343 = (t276 + t433) * t311;
t144 = t343 * t185;
t397 = t310 * t317;
t241 = rSges(6,1) * t398 + rSges(6,2) * t397;
t293 = pkin(4) * t396;
t207 = t293 + t241;
t91 = -t184 * t207 + t144;
t491 = -m(5) * t103 - m(6) * t91;
t465 = m(5) / 0.2e1;
t464 = m(6) / 0.2e1;
t490 = -t310 / 0.2e1;
t441 = t310 / 0.2e1;
t439 = t311 / 0.2e1;
t345 = -cos(qJ(1)) * pkin(1) - pkin(2) * cos(t319);
t437 = m(4) * (t345 * (rSges(4,1) * t310 + rSges(4,2) * t311) - (-rSges(4,1) * t311 + t310 * rSges(4,2)) * t346);
t307 = t310 ^ 2;
t308 = t311 ^ 2;
t354 = t307 + t308;
t440 = -t311 / 0.2e1;
t489 = t439 + t440;
t242 = t276 * t311;
t96 = -t184 * t241 + t185 * t242;
t488 = m(6) * t96;
t173 = -t241 * t310 - t311 * t242;
t402 = t276 * t310;
t222 = t293 + t402;
t277 = -rSges(6,2) * t316 + t427;
t400 = t277 * t311;
t401 = t277 * t310;
t338 = Icges(6,5) * t316 + Icges(6,6) * t317;
t235 = t310 * t338;
t236 = t338 * t311;
t280 = Icges(6,4) * t398;
t213 = -Icges(6,1) * t397 + Icges(6,5) * t311 + t280;
t366 = Icges(6,2) * t397 + t213 + t280;
t483 = Icges(6,1) * t316 + t309;
t368 = -t310 * t483 + t211;
t475 = -t366 * t316 - t368 * t317;
t365 = -Icges(6,2) * t391 + t214 - t281;
t367 = t311 * t483 + t212;
t476 = t365 * t316 + t367 * t317;
t431 = (t308 * t235 + (t476 * t310 + (-t236 - t475) * t311) * t310) * t439 + (-t307 * t236 + (t475 * t311 + (t235 - t476) * t310) * t311) * t441;
t198 = t311 * (t310 * rSges(6,3) - t347);
t217 = -rSges(6,1) * t397 - t358;
t95 = -t311 * (t311 * pkin(3) - t282 + (pkin(7) + t325) * t310) + t198 + (t388 - t217 + t306 + (-pkin(3) + t312) * t310) * t310;
t15 = t431 + m(6) * (t95 * t173 + t222 * t401 + t343 * t400);
t487 = t15 * qJD(5);
t301 = -rSges(5,2) * t321 + t428;
t486 = t301 * t465;
t389 = t311 * t323;
t229 = Icges(5,4) * t389 - Icges(5,2) * t390 + Icges(5,6) * t310;
t290 = Icges(5,4) * t390;
t231 = Icges(5,1) * t389 + Icges(5,5) * t310 - t290;
t485 = (t229 * t321 - t231 * t323) * t311;
t172 = t185 + t345;
t129 = t343 * t172;
t84 = -t185 * t171 + t172 * t184;
t182 = t191 + t346;
t183 = t192 + t345;
t89 = -t192 * t182 + t183 * t191;
t318 = Icges(5,4) * t323;
t482 = Icges(5,1) * t321 + t318;
t480 = Icges(5,2) * t321 - t318;
t479 = qJD(1) + qJD(3);
t414 = t182 * t249;
t430 = (t492 * t222 + t129 - t144) * t464 + (-t414 + t411 + (t183 - t192) * t250) * t465;
t88 = -t171 * t207 + t129;
t97 = t183 * t250 - t414;
t478 = (t91 + t88) * t464 + (t103 + t97) * t465;
t418 = Icges(6,4) * t316;
t272 = Icges(6,2) * t317 + t418;
t275 = Icges(6,1) * t317 - t418;
t477 = (-t481 + t483) * t316 + (t272 - t275) * t317;
t419 = Icges(5,4) * t321;
t295 = Icges(5,2) * t323 + t419;
t298 = Icges(5,1) * t323 - t419;
t474 = (-t480 + t482) * t321 + (t295 - t298) * t323;
t361 = -Icges(5,2) * t389 + t231 - t290;
t363 = t311 * t482 + t229;
t473 = t361 * t321 + t363 * t323;
t289 = Icges(5,4) * t396;
t395 = t310 * t323;
t230 = -Icges(5,1) * t395 + Icges(5,5) * t311 + t289;
t362 = Icges(5,2) * t395 + t230 + t289;
t228 = Icges(5,6) * t311 + t480 * t310;
t364 = -t310 * t482 + t228;
t472 = -t362 * t321 - t364 * t323;
t471 = (t482 / 0.2e1 - t480 / 0.2e1) * t323 + (t298 / 0.2e1 - t295 / 0.2e1) * t321;
t349 = (-t481 / 0.2e1 + t483 / 0.2e1) * t317 + (-t272 / 0.2e1 + t275 / 0.2e1) * t316;
t104 = -t213 * t397 + t372;
t105 = t311 * t210 + t212 * t398 - t214 * t397;
t351 = -t213 * t317 - t210;
t408 = t211 * t316;
t11 = ((t484 + t493) * t311 + ((t351 - t408) * t311 + t105 + t494) * t310) * t441 + (t104 * t311 + t105 * t310) * t490 + ((t105 + (-t210 + t408) * t311 - t494) * t311 + (t351 * t310 - t104 + t493) * t310) * t439;
t470 = 4 * qJD(1);
t467 = 2 * qJD(4);
t461 = m(5) * t89;
t459 = m(5) * t97;
t90 = -t171 * t241 + t172 * t242;
t454 = m(6) * (t96 + t90);
t452 = m(6) * ((t172 - t185) * t242 + t492 * t402);
t147 = t172 * t401;
t373 = t222 * t242 - t241 * t343;
t451 = m(6) * (t171 * t400 + t147 + t373);
t180 = t343 * t402;
t409 = t207 * t276;
t415 = t171 * t277;
t450 = m(6) * (t147 + t180 + (-t409 + t415) * t311);
t160 = t185 * t401;
t449 = m(6) * (t184 * t400 + t160 + t373);
t413 = t184 * t277;
t448 = m(6) * (t160 + t180 + (-t409 + t413) * t311);
t447 = m(6) * t84;
t445 = m(6) * t88;
t444 = m(6) * t90;
t227 = Icges(5,5) * t389 - Icges(5,6) * t390 + Icges(5,3) * t310;
t117 = t311 * t227 + t229 * t396 - t231 * t395;
t294 = Icges(5,5) * t323 - Icges(5,6) * t321;
t399 = t294 * t310;
t226 = Icges(5,3) * t311 - t399;
t369 = t310 * t226 + t230 * t389;
t118 = -t228 * t390 + t369;
t119 = t227 * t310 - t485;
t350 = -t230 * t323 - t227;
t370 = t311 * t226 + t228 * t396;
t407 = t228 * t321;
t28 = (t119 + t370 + t485) * t311 + (-t118 + (t350 - t407) * t311 + t117 + t369) * t310;
t116 = -t230 * t395 + t370;
t29 = (t117 + (-t227 + t407) * t311 - t369) * t311 + (t350 * t310 - t116 + t370) * t310;
t81 = t116 * t311 + t117 * t310;
t82 = t118 * t311 + t119 * t310;
t2 = (t29 / 0.2e1 + t82 / 0.2e1) * t311 + (-t81 / 0.2e1 + t28 / 0.2e1) * t310 + t11;
t438 = t2 * qJD(4);
t434 = m(6) * t173;
t422 = t11 * qJD(5);
t406 = t241 * t276;
t374 = (-t207 + t222) * t343;
t352 = t277 + t432;
t342 = t454 / 0.2e1 + t349;
t339 = Icges(5,5) * t321 + Icges(5,6) * t323;
t243 = t310 * t339;
t332 = -t11 + (-t477 * t311 - t367 * t316 + t365 * t317 + t404) * t441 + (t271 * t311 + t477 * t310 - t368 * t316 + t366 * t317) * t439;
t331 = -t349 + t489 * (t317 * t212 + t316 * t214);
t330 = t349 + t471;
t328 = t330 + t478;
t327 = t331 - t471 + t489 * (t323 * t229 + t321 * t231);
t326 = (t332 + t28 * t490 + (t29 + t82) * t440 + (t294 * t311 + t474 * t310 - t364 * t321 + t362 * t323) * t439 + (-t474 * t311 - t363 * t321 + t361 * t323 + t399 + t81) * t441) * qJD(4);
t244 = t339 * t311;
t225 = t352 * t311;
t223 = t352 * t310;
t187 = t242 * t402;
t177 = -t249 * t310 - t250 * t311;
t164 = qJD(5) * t434;
t151 = -t310 * t217 + t198;
t145 = -t354 * t433 + t173;
t83 = t349 + t488;
t70 = t448 / 0.2e1;
t69 = t349 + t444;
t68 = t449 / 0.2e1;
t58 = t450 / 0.2e1;
t56 = t451 / 0.2e1;
t53 = t452 / 0.2e1;
t39 = t330 - t491;
t34 = t330 + t445 + t459;
t33 = t437 + t447 + t461;
t19 = -t452 / 0.2e1 + t342;
t18 = t53 + t342;
t17 = m(6) * (t354 * t276 * t277 + t151 * t173) + t431;
t16 = t17 * qJD(5);
t14 = t53 - t454 / 0.2e1 + t331;
t13 = t328 + t430;
t12 = t328 - t430;
t9 = t327 + t430 - t478;
t8 = t68 - t448 / 0.2e1 + t11;
t7 = t70 - t449 / 0.2e1 + t11;
t6 = t56 - t450 / 0.2e1 + t11;
t5 = t58 - t451 / 0.2e1 + t11;
t4 = t68 + t70 + t332;
t3 = t56 + t58 + t332;
t1 = [qJD(3) * t33 + qJD(4) * t34 + qJD(5) * t69, 0, t33 * qJD(1) + t13 * qJD(4) + t18 * qJD(5) + 0.2e1 * (t437 / 0.2e1 + t84 * t464 + t89 * t465) * qJD(3), t34 * qJD(1) + t13 * qJD(3) + t3 * qJD(5) + ((t182 * t311 + t183 * t310) * t486 + (t171 * t225 + t172 * t223 + t374) * t464) * t467 + t326, t69 * qJD(1) + t18 * qJD(3) + t3 * qJD(4) + ((t147 + t187 + (-t406 + t415) * t311) * m(6) + t332) * qJD(5); 0, 0, 0, (t145 * t464 + t177 * t465) * t467 + t164, qJD(4) * t434 + t164; t12 * qJD(4) + t19 * qJD(5) + (-t437 / 0.4e1 - t461 / 0.4e1 - t447 / 0.4e1) * t470, 0, qJD(4) * t39 + qJD(5) * t83, t12 * qJD(1) + t39 * qJD(3) + t4 * qJD(5) + ((t184 * t225 + t185 * t223 + t374) * t464 + (t191 * t311 + t192 * t310) * t486) * t467 + t326, t19 * qJD(1) + t83 * qJD(3) + t4 * qJD(4) + ((t160 + t187 + (-t406 + t413) * t311) * m(6) + t332) * qJD(5); t9 * qJD(3) + t6 * qJD(5) + (-t445 / 0.4e1 - t459 / 0.4e1) * t470 + t438 + t327 * qJD(1), 0, t9 * qJD(1) + t8 * qJD(5) + t438 + (t327 + t491) * qJD(3), (m(5) * (t354 * t301 * t299 + (t311 * (rSges(5,1) * t389 + rSges(5,3) * t310 - t292) - t310 * (-rSges(5,1) * t395 - t357)) * t177) + (t308 * t243 + (t473 * t310 + (-t244 - t472) * t311) * t310) * t439 + (-t307 * t244 + (t472 * t311 + (t243 - t473) * t310) * t311) * t441 + m(6) * (t145 * t95 + t222 * t223 + t225 * t343) + t431) * qJD(4) + t487 + t479 * t2, t6 * qJD(1) + t8 * qJD(3) + t15 * qJD(4) + t487; (t331 - t444) * qJD(1) + t14 * qJD(3) + t5 * qJD(4) + t422, 0, t14 * qJD(1) + (t331 - t488) * qJD(3) + t7 * qJD(4) + t422, t5 * qJD(1) + t7 * qJD(3) + ((t145 * t151 + (t223 * t310 + t225 * t311) * t276) * m(6) + t431) * qJD(4) + t16, qJD(4) * t17 + t479 * t11 + t16;];
Cq = t1;
