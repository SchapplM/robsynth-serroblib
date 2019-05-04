% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR10V2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:39
% EndTime: 2019-04-11 14:42:44
% DurationCPUTime: 21.95s
% Computational Cost: add. (17343->808), mult. (40829->1157), div. (0->0), fcn. (31393->10), ass. (0->368)
t295 = sin(qJ(4));
t392 = qJD(4) * t295;
t296 = sin(qJ(3));
t297 = sin(qJ(2));
t301 = cos(qJ(3));
t302 = cos(qJ(2));
t264 = t296 * t297 - t301 * t302;
t258 = t264 * qJD(1);
t415 = t258 * t295;
t543 = t415 + t392;
t292 = qJD(2) + qJD(3);
t229 = t292 * t264;
t217 = t229 * qJD(1);
t265 = t296 * t302 + t297 * t301;
t259 = t265 * qJD(1);
t300 = cos(qJ(4));
t235 = t259 * t300 + t292 * t295;
t162 = qJD(4) * t235 - t217 * t295;
t491 = -t162 / 0.2e1;
t234 = -t259 * t295 + t292 * t300;
t161 = qJD(4) * t234 - t217 * t300;
t250 = qJD(4) + t258;
t294 = sin(qJ(5));
t299 = cos(qJ(5));
t188 = t235 * t299 + t250 * t294;
t230 = t292 * t265;
t218 = t230 * qJD(1);
t83 = qJD(5) * t188 + t161 * t294 - t218 * t299;
t500 = t83 / 0.2e1;
t225 = pkin(3) * t259 + pkin(5) * t258;
t432 = pkin(2) * qJD(2);
t382 = t301 * t432;
t196 = t225 * t295 + t300 * t382;
t386 = qJD(5) * t300;
t391 = qJD(4) * t299;
t315 = -t294 * t386 - t295 * t391;
t383 = t296 * t432;
t387 = qJD(5) * t299;
t529 = -pkin(3) * t387 + pkin(5) * t315 - t196 * t299 - t294 * t383;
t560 = t543 * pkin(6);
t402 = t299 * t300;
t201 = -t258 * t402 + t259 * t294;
t388 = qJD(5) * t295;
t366 = t294 * t388;
t559 = (t201 + t366) * pkin(6);
t470 = pkin(2) * t296;
t286 = pkin(5) + t470;
t469 = pkin(2) * t301;
t287 = -pkin(3) - t469;
t431 = pkin(2) * qJD(3);
t380 = t301 * t431;
t357 = t300 * t380;
t381 = t296 * t431;
t181 = t286 * t315 + t287 * t387 + t294 * t381 + t299 * t357;
t397 = qJD(1) * t297;
t384 = pkin(2) * t397;
t205 = t225 + t384;
t405 = t295 * t299;
t558 = -t205 * t405 + t181;
t187 = -t235 * t294 + t250 * t299;
t82 = qJD(5) * t187 + t161 * t299 + t218 * t294;
t557 = Ifges(6,6) * t491 - t82 * Ifges(6,4) / 0.2e1;
t195 = -t225 * t300 + t295 * t382;
t466 = pkin(6) * t299;
t374 = pkin(5) - t466;
t390 = qJD(4) * t300;
t556 = t374 * t390 - t195 + t559;
t358 = t295 * t380;
t359 = t286 - t466;
t417 = t205 * t300;
t555 = t359 * t390 + t358 + t417 + t559;
t554 = -t560 - t529;
t553 = -t558 - t560;
t293 = sin(qJ(6));
t298 = cos(qJ(6));
t406 = t295 * t298;
t164 = -t201 * t293 - t258 * t406;
t355 = -qJD(6) + t391;
t356 = -qJD(6) * t299 + qJD(4);
t199 = t356 * t406 + (-t300 * t355 + t366) * t293;
t552 = t164 - t199;
t411 = t293 * t295;
t165 = t201 * t298 - t258 * t411;
t389 = qJD(5) * t294;
t403 = t298 * t300;
t198 = t355 * t403 + (t293 * t356 - t298 * t389) * t295;
t551 = t165 - t198;
t284 = pkin(5) * t402;
t369 = t294 * t392;
t550 = -pkin(3) * t389 - pkin(5) * t369 + qJD(5) * t284 - t196 * t294 + t299 * t383;
t270 = pkin(5) * t292 + t383;
t288 = -pkin(2) * t302 - pkin(1);
t275 = t288 * qJD(1);
t309 = pkin(3) * t258 - pkin(5) * t259 + t275;
t172 = t270 * t300 + t295 * t309;
t464 = t292 * pkin(3);
t326 = t382 + t464;
t318 = t299 * t326;
t379 = qJD(2) * t431;
t354 = t296 * t379;
t385 = qJD(1) * qJD(2);
t364 = t297 * t385;
t311 = pkin(2) * t364 + pkin(3) * t218 + pkin(5) * t217;
t353 = t301 * t379;
t171 = t270 * t295 - t300 * t309;
t394 = qJD(4) * t171;
t87 = t295 * t311 + t300 * t353 - t394;
t43 = -qJD(5) * t318 - t172 * t389 + t294 * t354 + t299 * t87;
t29 = pkin(6) * t162 + t43;
t393 = qJD(4) * t172;
t88 = t295 * t353 - t300 * t311 + t393;
t30 = -pkin(6) * t82 + t88;
t113 = -pkin(6) * t188 + t171;
t147 = t172 * t299 - t294 * t326;
t233 = qJD(5) - t234;
t116 = pkin(6) * t233 + t147;
t53 = t113 * t298 - t116 * t293;
t3 = qJD(6) * t53 + t29 * t298 + t293 * t30;
t54 = t113 * t293 + t116 * t298;
t4 = -qJD(6) * t54 - t29 * t293 + t298 * t30;
t352 = -mrSges(7,1) * t4 + mrSges(7,2) * t3;
t549 = -t43 * mrSges(6,3) - t352;
t141 = t188 * t298 + t233 * t293;
t28 = -qJD(6) * t141 + t162 * t298 - t293 * t82;
t25 = Ifges(7,6) * t28;
t140 = -t188 * t293 + t233 * t298;
t27 = qJD(6) * t140 + t162 * t293 + t298 * t82;
t26 = Ifges(7,5) * t27;
t8 = Ifges(7,3) * t83 + t25 + t26;
t548 = Ifges(6,2) * t500 + t8 / 0.2e1 + t557;
t540 = t43 * mrSges(6,2);
t424 = -mrSges(6,1) * t233 - mrSges(7,1) * t140 + mrSges(7,2) * t141 + mrSges(6,3) * t188;
t428 = t235 * mrSges(5,3);
t399 = -mrSges(5,1) * t250 - mrSges(6,1) * t187 + mrSges(6,2) * t188 + t428;
t409 = t294 * t300;
t200 = -t258 * t409 - t259 * t299;
t368 = t294 * t390;
t314 = t295 * t387 + t368;
t546 = t200 - t314;
t545 = -t88 * mrSges(6,1) - t549;
t426 = t88 * t295;
t324 = t171 * t390 + t426;
t544 = -t172 * t392 + t324;
t508 = t28 / 0.2e1;
t509 = t27 / 0.2e1;
t542 = Ifges(7,5) * t509 + Ifges(7,6) * t508 + Ifges(7,3) * t500 + t548;
t539 = t88 * mrSges(5,1);
t458 = -mrSges(6,1) * t162 - mrSges(7,1) * t28 + mrSges(7,2) * t27 + mrSges(6,3) * t82;
t245 = t286 * t402 + t287 * t294;
t465 = pkin(6) * t300;
t236 = t245 - t465;
t251 = t359 * t295;
t190 = -t236 * t293 + t251 * t298;
t538 = qJD(6) * t190 + t293 * t555 - t298 * t553;
t191 = t236 * t298 + t251 * t293;
t537 = -qJD(6) * t191 + t293 * t553 + t298 * t555;
t536 = t147 * mrSges(6,2);
t425 = -mrSges(5,1) * t218 + mrSges(6,1) * t83 + mrSges(6,2) * t82 + mrSges(5,3) * t161;
t267 = -pkin(3) * t294 + t284;
t254 = t267 - t465;
t268 = t374 * t295;
t221 = -t254 * t293 + t268 * t298;
t535 = qJD(6) * t221 + t293 * t556 - t298 * t554;
t222 = t254 * t298 + t268 * t293;
t534 = -qJD(6) * t222 + t293 * t554 + t298 * t556;
t146 = t172 * t294 + t318;
t346 = mrSges(7,1) * t293 + mrSges(7,2) * t298;
t325 = mrSges(6,3) + t346;
t533 = t325 * t146;
t532 = t399 * t300;
t414 = t258 * t300;
t531 = mrSges(5,1) * t259 - mrSges(6,1) * t200 - mrSges(6,2) * t201 + mrSges(5,3) * t414;
t456 = mrSges(4,3) * t259;
t530 = -mrSges(4,1) * t292 - mrSges(5,1) * t234 + mrSges(5,2) * t235 + t456;
t185 = Ifges(6,4) * t187;
t440 = Ifges(6,5) * t233;
t454 = Ifges(6,1) * t188;
t110 = t185 + t440 + t454;
t497 = -t110 / 0.2e1;
t528 = t497 - t440 / 0.2e1 - t171 * mrSges(6,2) - t185 / 0.2e1;
t266 = pkin(3) * t299 + pkin(5) * t409;
t44 = qJD(5) * t147 + t294 * t87 - t299 * t354;
t527 = t146 * t550 + t266 * t44;
t526 = t286 * t390 + t358;
t400 = t300 * t229;
t321 = -t265 * t392 - t400;
t524 = -m(4) * t275 - mrSges(4,1) * t258 - mrSges(4,2) * t259;
t108 = Ifges(6,5) * t188 + Ifges(6,6) * t187 + Ifges(6,3) * t233;
t430 = t162 * Ifges(6,5);
t461 = t83 * Ifges(6,4);
t463 = t82 * Ifges(6,1);
t21 = t430 - t461 + t463;
t522 = t88 * mrSges(6,2) + t21 / 0.2e1;
t332 = t293 * t54 + t298 * t53;
t335 = Ifges(7,5) * t298 - Ifges(7,6) * t293;
t443 = Ifges(7,4) * t298;
t339 = -Ifges(7,2) * t293 + t443;
t444 = Ifges(7,4) * t293;
t343 = Ifges(7,1) * t298 - t444;
t473 = t298 / 0.2e1;
t474 = -t293 / 0.2e1;
t186 = qJD(6) - t187;
t487 = t186 / 0.2e1;
t493 = t141 / 0.2e1;
t495 = t140 / 0.2e1;
t445 = Ifges(7,4) * t141;
t56 = Ifges(7,2) * t140 + Ifges(7,6) * t186 + t445;
t134 = Ifges(7,4) * t140;
t57 = Ifges(7,1) * t141 + Ifges(7,5) * t186 + t134;
t521 = -t332 * mrSges(7,3) + t335 * t487 + t339 * t495 + t343 * t493 + t473 * t57 + t474 * t56;
t436 = Ifges(6,6) * t233;
t438 = Ifges(6,2) * t187;
t448 = Ifges(6,4) * t188;
t109 = t436 + t438 + t448;
t152 = Ifges(5,5) * t235 + Ifges(5,6) * t234 + Ifges(5,3) * t250;
t451 = Ifges(5,4) * t235;
t153 = Ifges(5,2) * t234 + Ifges(5,6) * t250 + t451;
t232 = Ifges(5,4) * t234;
t154 = Ifges(5,1) * t235 + Ifges(5,5) * t250 + t232;
t19 = Ifges(6,5) * t82 - Ifges(6,6) * t83 + Ifges(6,3) * t162;
t452 = Ifges(4,4) * t259;
t206 = -Ifges(4,2) * t258 + Ifges(4,6) * t292 + t452;
t249 = Ifges(4,4) * t258;
t207 = Ifges(4,1) * t259 + t292 * Ifges(4,5) - t249;
t260 = -t293 * t405 - t403;
t404 = t298 * t299;
t261 = -t293 * t300 + t295 * t404;
t316 = t299 * t390 - t366;
t336 = Ifges(6,5) * t299 - Ifges(6,6) * t294;
t337 = Ifges(5,5) * t300 - Ifges(5,6) * t295;
t446 = Ifges(6,4) * t299;
t340 = -Ifges(6,2) * t294 + t446;
t449 = Ifges(5,4) * t300;
t341 = -Ifges(5,2) * t295 + t449;
t447 = Ifges(6,4) * t294;
t344 = Ifges(6,1) * t299 - t447;
t450 = Ifges(5,4) * t295;
t345 = Ifges(5,1) * t300 - t450;
t348 = mrSges(5,1) * t295 + mrSges(5,2) * t300;
t361 = t390 / 0.2e1;
t362 = -t392 / 0.2e1;
t410 = t294 * t295;
t427 = t87 * t300;
t457 = mrSges(4,3) * t258;
t472 = t300 / 0.2e1;
t475 = t259 / 0.2e1;
t476 = -t259 / 0.2e1;
t477 = t258 / 0.2e1;
t478 = -t250 / 0.2e1;
t480 = -t235 / 0.2e1;
t481 = -t234 / 0.2e1;
t482 = t233 / 0.2e1;
t484 = t218 / 0.2e1;
t485 = t188 / 0.2e1;
t486 = t187 / 0.2e1;
t488 = -t186 / 0.2e1;
t490 = t162 / 0.2e1;
t492 = t161 / 0.2e1;
t494 = -t141 / 0.2e1;
t496 = -t140 / 0.2e1;
t498 = t109 / 0.2e1;
t499 = -t109 / 0.2e1;
t501 = -t83 / 0.2e1;
t502 = t82 / 0.2e1;
t503 = t57 / 0.2e1;
t504 = -t57 / 0.2e1;
t505 = t56 / 0.2e1;
t506 = -t56 / 0.2e1;
t433 = Ifges(7,3) * t186;
t435 = Ifges(7,6) * t140;
t439 = Ifges(7,5) * t141;
t55 = t433 + t435 + t439;
t507 = -t55 / 0.2e1;
t10 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t83;
t512 = t10 / 0.2e1;
t9 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t83;
t513 = t9 / 0.2e1;
t68 = Ifges(5,4) * t161 - Ifges(5,2) * t162 + Ifges(5,6) * t218;
t69 = Ifges(5,1) * t161 - Ifges(5,4) * t162 + Ifges(5,5) * t218;
t519 = -t172 * (-mrSges(5,2) * t259 + mrSges(5,3) * t415) - t233 * (Ifges(6,5) * t201 - Ifges(6,6) * t200 - Ifges(6,3) * t415) / 0.2e1 - t187 * (Ifges(6,4) * t201 - Ifges(6,2) * t200 - Ifges(6,6) * t415) / 0.2e1 - t188 * (Ifges(6,1) * t201 - Ifges(6,4) * t200 - Ifges(6,5) * t415) / 0.2e1 - t275 * (mrSges(4,1) * t259 - mrSges(4,2) * t258) - t292 * (-Ifges(4,5) * t258 - Ifges(4,6) * t259) / 0.2e1 + (Ifges(5,3) * t259 - t258 * t337) * t478 + (Ifges(5,5) * t259 - t258 * t345) * t480 + (Ifges(5,6) * t259 - t258 * t341) * t481 + (Ifges(7,4) * t261 + Ifges(7,2) * t260) * t508 + (-Ifges(4,1) * t258 + t152 - t452) * t476 + (-t415 / 0.2e1 + t362) * t153 + (Ifges(7,1) * t261 + Ifges(7,4) * t260) * t509 + (-Ifges(4,2) * t259 + t207 - t249) * t477 + (-mrSges(5,1) * t300 + mrSges(5,2) * t295) * t354 + t154 * t414 / 0.2e1 + (t110 * t299 + t294 * t55 + t154) * t361 + t300 * t540 + (Ifges(7,5) * t261 + Ifges(7,6) * t260) * t500 + t543 * t108 / 0.2e1 + t295 * t69 / 0.2e1 - Ifges(4,6) * t218 - Ifges(4,5) * t217 - t300 * t19 / 0.2e1 + (t261 * t44 + t54 * t546) * mrSges(7,2) + (-mrSges(6,2) * t543 + mrSges(6,3) * t546) * t147 + (-t260 * t44 - t53 * t546) * mrSges(7,1) - t44 * (-mrSges(6,1) * t300 - mrSges(6,3) * t405) - (t109 * t299 + t110 * t294) * t388 / 0.2e1 + (qJD(5) * t55 + t21) * t405 / 0.2e1 + (t234 * t341 + t235 * t345 + t250 * t337) * qJD(4) / 0.2e1 - t326 * t348 * qJD(4) + (t542 + t549) * t410 + ((-t201 + t316) * mrSges(6,3) - t551 * mrSges(7,2) + t552 * mrSges(7,1) - t543 * mrSges(6,1)) * t146 + (t260 * t3 - t261 * t4 + t53 * t551 - t54 * t552) * mrSges(7,3) + (mrSges(6,1) * t294 + mrSges(6,2) * t299) * t426 + mrSges(5,3) * t427 - t382 * t457 + t68 * t472 + t206 * t475 + ((-Ifges(6,5) * t294 - Ifges(6,6) * t299) * t388 + (Ifges(6,3) * t295 + t300 * t336) * qJD(4)) * t482 + (Ifges(5,5) * t295 + Ifges(5,6) * t300) * t484 + ((-Ifges(6,1) * t294 - t446) * t388 + (Ifges(6,5) * t295 + t300 * t344) * qJD(4)) * t485 + ((-Ifges(6,2) * t299 - t447) * t388 + (Ifges(6,6) * t295 + t300 * t340) * qJD(4)) * t486 + (Ifges(7,5) * t198 + Ifges(7,6) * t199 + Ifges(7,3) * t314) * t487 + (Ifges(7,5) * t165 + Ifges(7,6) * t164 + Ifges(7,3) * t200) * t488 + (-Ifges(6,3) * t300 + t295 * t336) * t490 + (Ifges(5,2) * t300 + t450) * t491 + (Ifges(5,1) * t295 + t449) * t492 + (Ifges(7,1) * t198 + Ifges(7,4) * t199 + Ifges(7,5) * t314) * t493 + (Ifges(7,1) * t165 + Ifges(7,4) * t164 + Ifges(7,5) * t200) * t494 + (Ifges(7,4) * t198 + Ifges(7,2) * t199 + Ifges(7,6) * t314) * t495 + (Ifges(7,4) * t165 + Ifges(7,2) * t164 + Ifges(7,6) * t200) * t496 + t201 * t497 + t200 * t498 + t368 * t499 + (-Ifges(6,6) * t300 + t295 * t340) * t501 + (-Ifges(6,5) * t300 + t295 * t344) * t502 + t198 * t503 + t165 * t504 + t199 * t505 + t164 * t506 + t200 * t507 + t261 * t512 + t260 * t513 + t171 * (mrSges(6,1) * t314 + mrSges(6,2) * t316);
t518 = -mrSges(6,1) * t171 - mrSges(7,1) * t53 + mrSges(7,2) * t54 + mrSges(6,3) * t147;
t517 = -0.2e1 * pkin(1);
t515 = m(6) / 0.2e1;
t483 = -t229 / 0.2e1;
t479 = t235 / 0.2e1;
t468 = pkin(6) * t187;
t460 = t87 * mrSges(5,2);
t455 = mrSges(5,3) * t234;
t453 = Ifges(3,4) * t297;
t423 = Ifges(3,5) * qJD(2);
t422 = Ifges(3,6) * qJD(2);
t420 = t146 * t294;
t419 = t171 * t195;
t418 = t171 * t300;
t416 = t234 * t299;
t413 = t265 * t295;
t192 = -mrSges(5,2) * t250 + t455;
t408 = t295 * t192;
t407 = t295 * t229;
t118 = -mrSges(5,2) * t218 - mrSges(5,3) * t162;
t401 = t300 * t118;
t396 = qJD(1) * t302;
t395 = qJD(2) * t297;
t378 = t146 * t410;
t377 = t205 * t410;
t375 = Ifges(5,5) * t161 - Ifges(5,6) * t162 + Ifges(5,3) * t218;
t373 = qJD(4) * t420;
t351 = -t293 * t4 + t298 * t3;
t350 = -t293 * t3 - t298 * t4;
t111 = (-t265 * t386 + t230) * t294 + (qJD(5) * t264 + t321) * t299;
t349 = qJD(6) * t413 + t111;
t347 = mrSges(7,1) * t298 - mrSges(7,2) * t293;
t342 = Ifges(7,1) * t293 + t443;
t338 = Ifges(7,2) * t298 + t444;
t334 = Ifges(7,5) * t293 + Ifges(7,6) * t298;
t226 = pkin(3) * t264 - pkin(5) * t265 + t288;
t333 = pkin(6) * t265 + t226 * t299;
t331 = t293 * t53 - t298 * t54;
t182 = qJD(5) * t245 - t286 * t369 + t294 * t357 - t299 * t381;
t244 = t286 * t409 - t287 * t299;
t330 = t146 * t182 + t244 * t44;
t220 = t264 * t294 + t265 * t402;
t155 = -pkin(6) * t220 - t226 * t300;
t179 = t333 * t295;
t102 = t155 * t298 - t179 * t293;
t103 = t155 * t293 + t179 * t298;
t328 = -t172 * t295 + t418;
t327 = t171 * t295 + t172 * t300;
t323 = t146 * t387 + t294 * t44;
t322 = t265 * t390 - t407;
t320 = m(5) * t328;
t142 = -mrSges(6,2) * t233 + mrSges(6,3) * t187;
t313 = t142 * t299 + t294 * t424 + t192;
t312 = -qJD(6) * t220 + t322;
t310 = qJD(4) * t328 + t426;
t308 = -t454 / 0.2e1 + t528;
t307 = -t433 / 0.2e1 - t435 / 0.2e1 - t439 / 0.2e1 + t436 / 0.2e1 + t507 + t498 + t448 / 0.2e1 + t518;
t305 = -t438 / 0.2e1 - t307;
t98 = -mrSges(7,2) * t186 + mrSges(7,3) * t140;
t99 = mrSges(7,1) * t186 - mrSges(7,3) * t141;
t304 = (-m(7) * t332 - t293 * t98 - t298 * t99) * pkin(6) + t521;
t289 = Ifges(3,4) * t396;
t257 = Ifges(3,1) * t397 + t289 + t423;
t256 = t422 + (t302 * Ifges(3,2) + t453) * qJD(1);
t238 = -mrSges(4,2) * t292 - t457;
t208 = t348 * t258;
t176 = t220 * t298 + t265 * t411;
t175 = -t220 * t293 + t265 * t406;
t169 = t234 * t404 + t235 * t293;
t168 = t235 * t298 - t293 * t416;
t160 = pkin(2) * t395 + pkin(3) * t230 + pkin(5) * t229;
t131 = -pkin(6) * t416 + t172;
t127 = pkin(6) * t235 - t171 * t299;
t105 = -t146 * t298 - t293 * t468;
t104 = t146 * t293 - t298 * t468;
t101 = mrSges(5,1) * t162 + mrSges(5,2) * t161;
t63 = t127 * t298 + t131 * t293;
t62 = -t127 * t293 + t131 * t298;
t58 = t333 * t390 + (-pkin(6) * t229 + t160 * t299 - t226 * t389) * t295;
t51 = -mrSges(6,2) * t162 - mrSges(6,3) * t83;
t47 = -pkin(6) * t111 - t160 * t300 + t226 * t392;
t46 = -t293 * t349 + t298 * t312;
t45 = t293 * t312 + t298 * t349;
t15 = -mrSges(7,2) * t83 + mrSges(7,3) * t28;
t14 = mrSges(7,1) * t83 - mrSges(7,3) * t27;
t13 = -qJD(6) * t103 - t293 * t58 + t298 * t47;
t12 = qJD(6) * t102 + t293 * t47 + t298 * t58;
t1 = [(t217 * t264 - t218 * t265 + t230 * t476 - t258 * t483) * Ifges(4,4) + (-t217 * t265 - t229 * t475) * Ifges(4,1) + t288 * (mrSges(4,1) * t218 - mrSges(4,2) * t217) + (t229 * t301 - t230 * t296 + (-t264 * t301 + t265 * t296) * qJD(3)) * mrSges(4,3) * t432 + (-t532 + t313 * t295 - t320 + m(6) * (t147 * t405 + t378 - t418) + m(7) * t378) * t160 + (-t256 / 0.2e1 - t422 / 0.2e1 - t524 * pkin(2) + (mrSges(3,1) * t517 - 0.3e1 / 0.2e1 * t453 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t302 + (m(4) * t288 + mrSges(4,1) * t264 + mrSges(4,2) * t265) * pkin(2)) * qJD(1)) * t395 + (t55 / 0.2e1 - Ifges(6,6) * t482 - Ifges(6,4) * t485 - Ifges(6,2) * t486 + Ifges(7,3) * t487 + Ifges(7,5) * t493 + Ifges(7,6) * t495 + t499 - t518) * (qJD(5) * t220 - t230 * t299 + t294 * t321) + (t44 * mrSges(6,3) + Ifges(6,1) * t502 + Ifges(6,4) * t501 + Ifges(6,5) * t490 + t522) * t220 + (-mrSges(5,1) * t230 + mrSges(6,2) * t111) * t171 + (t218 * t264 + t230 * t477) * Ifges(4,2) + t292 * (-Ifges(4,5) * t229 - Ifges(4,6) * t230) / 0.2e1 + t275 * (mrSges(4,1) * t230 - mrSges(4,2) * t229) + (-t328 * t229 + (-qJD(4) * t327 - t295 * t87 + t300 * t88) * t265) * mrSges(5,3) + (-t400 / 0.2e1 + t265 * t362) * t154 + (t108 - t153) * (-t407 / 0.2e1 + t265 * t361) + m(7) * (t102 * t4 + t103 * t3 + t12 * t54 + t13 * t53) - t172 * mrSges(5,2) * t230 + (Ifges(7,5) * t45 + Ifges(7,6) * t46) * t487 + (Ifges(7,5) * t176 + Ifges(7,6) * t175) * t500 + (Ifges(7,4) * t45 + Ifges(7,2) * t46) * t495 + (Ifges(7,4) * t176 + Ifges(7,2) * t175) * t508 + ((t313 * qJD(4) + m(5) * (-t88 + t393) + m(6) * (t147 * t391 + t373 - t88) + m(7) * t373 - t425) * t300 + (t299 * t51 + t118 + t458 * t294 + t399 * qJD(4) + (-t142 * t294 + t299 * t424) * qJD(5) + m(5) * (t87 + t394) + m(6) * (-t147 * t389 + t299 * t43 + t323 + t394) + m(7) * t323) * t295) * t226 - t264 * t460 + t264 * t375 / 0.2e1 + (Ifges(6,5) * t111 + Ifges(6,3) * t322) * t482 + (-Ifges(6,4) * t502 - Ifges(6,2) * t501 - Ifges(6,6) * t490 + t542 - t545) * (-t264 * t299 + t265 * t409) + (Ifges(6,4) * t111 + Ifges(6,6) * t322) * t486 + t230 * t152 / 0.2e1 - t230 * t206 / 0.2e1 + t44 * (-mrSges(7,1) * t175 + mrSges(7,2) * t176) + (Ifges(7,1) * t45 + Ifges(7,4) * t46) * t493 + (Ifges(7,1) * t176 + Ifges(7,4) * t175) * t509 + t111 * t110 / 0.2e1 + t12 * t98 + t13 * t99 + t102 * t14 + t103 * t15 - t326 * (mrSges(5,1) * t322 + mrSges(5,2) * t321) + (-mrSges(6,1) * t322 - mrSges(7,1) * t46 + mrSges(7,2) * t45 + mrSges(6,3) * t111) * t146 + t250 * (Ifges(5,5) * t321 - Ifges(5,6) * t322 + Ifges(5,3) * t230) / 0.2e1 + t234 * (Ifges(5,4) * t321 - Ifges(5,2) * t322 + Ifges(5,6) * t230) / 0.2e1 + (t175 * t3 - t176 * t4 - t45 * t53 + t46 * t54) * mrSges(7,3) - t322 * t536 + (Ifges(5,1) * t321 - Ifges(5,4) * t322 + Ifges(5,5) * t230) * t479 + t207 * t483 + (Ifges(5,3) * t264 + t265 * t337) * t484 + (Ifges(5,6) * t264 + t265 * t341) * t491 + (Ifges(5,5) * t264 + t265 * t345) * t492 + t45 * t503 + t46 * t505 + t176 * t512 + t175 * t513 - t264 * t539 + (-t68 / 0.2e1 + t19 / 0.2e1 - t540 - mrSges(6,1) * t44 + Ifges(6,3) * t490 + Ifges(6,6) * t501 + Ifges(6,5) * t502) * t413 + t265 * t69 * t472 + t348 * t265 * t354 + (t257 / 0.2e1 + t423 / 0.2e1 + (mrSges(3,2) * t517 + 0.3e1 / 0.2e1 * Ifges(3,4) * t302) * qJD(1)) * qJD(2) * t302 + (Ifges(6,1) * t111 + Ifges(6,5) * t322) * t485; (t217 * t469 - t218 * t470) * mrSges(4,3) + (-t408 + t320 - m(6) * (-t418 + (t147 * t299 + t420) * t295)) * t205 - (-Ifges(3,2) * t397 + t257 + t289) * t396 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t297 + mrSges(3,2) * t302) - t297 * (Ifges(3,1) * t302 - t453) / 0.2e1) * qJD(1) ^ 2 + t530 * t381 + (m(6) * t526 + t531) * t171 + t524 * t384 + t385 * Ifges(3,5) * t302 / 0.2e1 + ((m(6) * t88 + t425) * t295 + t401 + m(5) * (t310 + t427) - t392 * t192) * t286 - Ifges(3,6) * t364 / 0.2e1 + t256 * t397 / 0.2e1 + t192 * t357 - mrSges(4,2) * t353 + t424 * (-t377 + t182) + m(6) * (t147 * t181 + t245 * t43 + t330) + t544 * mrSges(5,3) + t287 * t101 + t245 * t51 + t190 * t14 + t191 * t15 - t326 * t208 + t399 * (t417 + t526) + t558 * t142 + m(5) * (t327 * t301 + (-t464 + (t287 - t469) * qJD(2)) * t296) * t431 + t519 + t238 * t380 + t537 * t99 + t538 * t98 + t383 * t456 + (-t146 * t377 + t190 * t4 + t191 * t3 + t53 * t537 + t538 * t54 + t330) * m(7) + t458 * t244 - mrSges(4,1) * t354; -t399 * t195 + ((-mrSges(4,2) * qJD(3) - t208 - t238) * t301 + (-mrSges(4,1) * qJD(3) + t456 - t530) * t296) * t432 + t531 * t171 + t529 * t142 + t458 * t266 + t267 * t51 + t221 * t14 + t222 * t15 - t196 * t192 + t310 * mrSges(5,3) + t535 * t98 + t534 * t99 + (-t208 * t292 - t101) * pkin(3) + t424 * t550 + (t221 * t4 + t222 * t3 + t53 * t534 + t535 * t54 + t527) * m(7) + (t147 * t529 + t267 * t43 - t419 + t527) * m(6) + (-pkin(3) * t354 - t172 * t196 + t326 * t383 - t419) * m(5) + (m(5) * (t427 + t544) + 0.2e1 * t324 * t515 + t401 + t425 * t295 + (-t408 + t532) * qJD(4)) * pkin(5) + t519; -m(7) * (t53 * t62 + t54 * t63) + (t192 - t455) * t171 + (mrSges(6,1) * t235 + mrSges(7,1) * t168 - mrSges(7,2) * t169) * t146 + (-m(6) * t171 - t399 + t428) * t172 + (-t26 / 0.2e1 - t25 / 0.2e1 + (-Ifges(7,3) / 0.2e1 - Ifges(6,2) / 0.2e1) * t83 + (m(6) * t147 + t142) * t171 + (-t146 * mrSges(6,3) + t308) * t234 + t545 - t548 - t557) * t299 + (-t168 * t54 + t169 * t53) * mrSges(7,3) + (t343 * t509 + t339 * t508 + t335 * t500 + t430 / 0.2e1 - t461 / 0.2e1 + t463 / 0.2e1 + t9 * t474 + t10 * t473 + t325 * t44 + t350 * mrSges(7,3) + (0.2e1 * (m(7) / 0.2e1 + t515) * t146 + t424) * t171 + (m(7) * t350 - t298 * t14 - t293 * t15) * pkin(6) - t305 * t234 + (t334 * t488 + t338 * t496 + t342 * t494 + t57 * t474 + t298 * t506 + t146 * t347 + t331 * mrSges(7,3) + (m(7) * t331 + t293 * t99 - t298 * t98) * pkin(6)) * qJD(6) + t522) * t294 - t460 - t63 * t98 - t62 * t99 - t539 + t235 * t536 + t326 * (mrSges(5,1) * t235 + mrSges(5,2) * t234) + (t305 * t294 + (t304 - t308 + t533) * t299) * qJD(5) + t375 + (Ifges(5,5) * t234 - Ifges(5,6) * t235) * t478 + t153 * t479 + (Ifges(7,5) * t169 + Ifges(7,6) * t168) * t488 + (Ifges(7,1) * t169 + Ifges(7,4) * t168) * t494 + (Ifges(7,4) * t169 + Ifges(7,2) * t168) * t496 + t169 * t504 + t168 * t506 + (-Ifges(5,2) * t235 + t154 + t232) * t481 + (Ifges(5,1) * t234 + 0.2e1 * t108 - t451) * t480; -m(7) * (t104 * t53 + t105 * t54) + t19 + t9 * t473 + t307 * t188 + t293 * t512 + t146 * t142 - t105 * t98 - t104 * t99 - t540 + t351 * mrSges(7,3) + (m(7) * t351 - t14 * t293 + t15 * t298) * pkin(6) + (t146 * t346 + t304) * qJD(6) + (-mrSges(6,1) - t347) * t44 + t342 * t509 + t334 * t500 + t338 * t508 + ((Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t188 - t533 - t521 + t528) * t187 + (-m(7) * t146 - t424) * t147; -t146 * (mrSges(7,1) * t141 + mrSges(7,2) * t140) + (Ifges(7,1) * t140 - t445) * t494 + t56 * t493 + (Ifges(7,5) * t140 - Ifges(7,6) * t141) * t488 - t53 * t98 + t54 * t99 + (t140 * t53 + t141 * t54) * mrSges(7,3) - t352 + t8 + (-Ifges(7,2) * t141 + t134 + t57) * t496;];
tauc  = t1(:);
