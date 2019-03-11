% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:00:06
% EndTime: 2019-03-08 20:00:18
% DurationCPUTime: 6.03s
% Computational Cost: add. (7388->492), mult. (18760->698), div. (0->0), fcn. (18501->10), ass. (0->267)
t293 = sin(qJ(5));
t289 = t293 ^ 2;
t296 = cos(qJ(5));
t290 = t296 ^ 2;
t508 = t289 + t290;
t297 = cos(qJ(4));
t450 = -t297 / 0.2e1;
t496 = Ifges(7,4) + Ifges(6,5);
t421 = t296 * mrSges(6,2);
t424 = t293 * mrSges(6,1);
t243 = t421 + t424;
t456 = t243 / 0.2e1;
t420 = t296 * mrSges(7,3);
t423 = t293 * mrSges(7,1);
t242 = -t420 + t423;
t457 = t242 / 0.2e1;
t507 = t456 + t457;
t292 = sin(pkin(11));
t409 = sin(pkin(6));
t410 = cos(pkin(11));
t340 = t410 * t409;
t295 = sin(qJ(2));
t348 = t295 * t409;
t447 = cos(qJ(2));
t187 = t292 * t348 - t340 * t447;
t341 = t447 * t409;
t188 = t292 * t341 + t295 * t340;
t384 = t296 * t297;
t96 = -t187 * t384 + t188 * t293;
t413 = t96 * t296;
t389 = t293 * t297;
t95 = -t187 * t389 - t188 * t296;
t414 = t95 * t293;
t334 = t413 + t414;
t321 = t334 * pkin(9);
t335 = pkin(5) * t296 + qJ(6) * t293;
t237 = -pkin(4) - t335;
t294 = sin(qJ(4));
t396 = t237 * t294;
t400 = t187 * t294;
t475 = -m(7) / 0.2e1;
t478 = -m(6) / 0.2e1;
t497 = mrSges(7,2) + mrSges(6,3);
t506 = -t497 * (t413 / 0.2e1 + t414 / 0.2e1) + (pkin(4) * t400 + t321) * t478 + (-t187 * t396 + t321) * t475;
t477 = m(6) / 0.2e1;
t411 = cos(pkin(6));
t140 = t188 * t297 + t294 * t411;
t83 = t140 * t296 + t187 * t293;
t415 = t83 * t296;
t505 = t140 - t415;
t338 = t296 * mrSges(7,1) + t293 * mrSges(7,3);
t339 = t296 * mrSges(6,1) - t293 * mrSges(6,2);
t504 = -t339 - t338;
t284 = Ifges(7,5) * t293;
t492 = Ifges(7,1) * t296 + t284;
t193 = -t297 * Ifges(7,4) + t294 * t492;
t432 = Ifges(6,4) * t293;
t252 = Ifges(6,1) * t296 - t432;
t195 = -t297 * Ifges(6,5) + t294 * t252;
t370 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t343 = t370 * t297;
t503 = -t343 + t193 / 0.2e1 + t195 / 0.2e1;
t502 = t508 * pkin(9);
t484 = m(6) / 0.4e1 + m(7) / 0.4e1;
t500 = m(6) + m(7);
t499 = mrSges(6,1) + mrSges(7,1);
t498 = mrSges(6,2) - mrSges(7,3);
t495 = Ifges(7,2) + Ifges(6,3);
t276 = pkin(2) * t292 + pkin(8);
t393 = t276 * t296;
t347 = -qJ(6) + t393;
t277 = -pkin(2) * t410 - pkin(3);
t223 = -t297 * pkin(4) - t294 * pkin(9) + t277;
t391 = t293 * t223;
t124 = t297 * t347 + t391;
t137 = t276 * t384 + t391;
t383 = -t124 + t137;
t349 = t276 * t293 + pkin(5);
t385 = t296 * t223;
t127 = t297 * t349 - t385;
t136 = -t276 * t389 + t385;
t382 = t127 + t136;
t390 = t293 * t294;
t229 = t297 * mrSges(6,2) - mrSges(6,3) * t390;
t282 = t297 * mrSges(7,3);
t236 = -mrSges(7,2) * t390 - t282;
t494 = t229 + t236;
t493 = m(7) * t237 - t338;
t491 = Ifges(7,6) * t293 + t496 * t296;
t287 = Ifges(6,4) * t296;
t490 = -Ifges(6,2) * t293 + t287;
t251 = Ifges(6,1) * t293 + t287;
t438 = t294 * pkin(4);
t442 = pkin(9) * t297;
t258 = t438 - t442;
t395 = t258 * t296;
t149 = t276 * t390 + t395;
t228 = t293 * t258;
t387 = t294 * t296;
t150 = -t276 * t387 + t228;
t486 = -t149 * t293 + t150 * t296;
t138 = -t294 * t347 + t228;
t141 = -t294 * t349 - t395;
t329 = t138 * t296 + t141 * t293;
t485 = Ifges(7,3) * t296 - t284;
t369 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t281 = m(7) * qJ(6) + mrSges(7,3);
t483 = t297 * t484;
t470 = -mrSges(7,3) / 0.2e1;
t471 = mrSges(6,2) / 0.2e1;
t371 = t471 + t470;
t372 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t482 = -t372 * t293 - t371 * t296 + t507;
t217 = t297 * t242;
t218 = t297 * t243;
t355 = t236 / 0.2e1 + t229 / 0.2e1;
t231 = -t297 * mrSges(6,1) - mrSges(6,3) * t387;
t232 = t297 * mrSges(7,1) + mrSges(7,2) * t387;
t358 = t231 / 0.2e1 - t232 / 0.2e1;
t306 = t293 * t358 - t296 * t355 + t217 / 0.2e1 + t218 / 0.2e1;
t139 = t188 * t294 - t297 * t411;
t386 = t296 * qJ(6);
t439 = t293 * pkin(5);
t241 = -t386 + t439;
t325 = t241 + t276;
t153 = t325 * t294;
t154 = t325 * t297;
t320 = t124 * t296 + t127 * t293 - t154;
t330 = -t136 * t293 + t137 * t296;
t230 = -t294 * mrSges(6,2) - mrSges(6,3) * t389;
t235 = -mrSges(7,2) * t389 + t294 * mrSges(7,3);
t356 = t235 / 0.2e1 + t230 / 0.2e1;
t233 = t294 * mrSges(6,1) - mrSges(6,3) * t384;
t374 = mrSges(7,2) * t384;
t422 = t294 * mrSges(7,1);
t234 = t374 - t422;
t461 = t234 / 0.2e1;
t357 = -t233 / 0.2e1 + t461;
t215 = t242 * t294;
t216 = t243 * t294;
t360 = t215 / 0.2e1 + t216 / 0.2e1;
t394 = t276 * t294;
t474 = m(7) / 0.2e1;
t82 = t140 * t293 - t187 * t296;
t481 = t360 * t140 + t357 * t82 + t356 * t83 + (t140 * t394 - t149 * t82 + t150 * t83) * t477 + (t138 * t83 + t140 * t153 + t141 * t82) * t474 + ((t276 * t297 - t330) * t477 + t320 * t475 + t306) * t139;
t480 = 0.2e1 * m(7);
t479 = 2 * qJD(4);
t472 = -mrSges(6,1) / 0.2e1;
t467 = -t82 / 0.2e1;
t466 = t83 / 0.2e1;
t464 = -qJ(6) / 0.2e1;
t463 = t137 / 0.2e1;
t462 = t141 / 0.2e1;
t460 = t237 / 0.2e1;
t459 = -t338 / 0.2e1;
t458 = t241 / 0.2e1;
t419 = t297 * mrSges(5,2);
t244 = t294 * mrSges(5,1) + t419;
t455 = t244 / 0.2e1;
t454 = t293 / 0.2e1;
t453 = t294 / 0.2e1;
t452 = -t296 / 0.2e1;
t451 = t296 / 0.2e1;
t445 = m(7) * t241;
t444 = m(7) * t297;
t319 = m(6) * t486;
t291 = t297 ^ 2;
t392 = t291 * t276;
t16 = t392 * t478 + (t356 * t296 + t357 * t293 + (t153 + t329) * t474 + t319 / 0.2e1 + t394 * t477 + t360) * t294 + (t320 * t474 + t330 * t477 - t306) * t297;
t214 = t335 * t294;
t425 = t290 * mrSges(6,3);
t426 = t290 * mrSges(7,2);
t427 = t289 * mrSges(6,3);
t428 = t289 * mrSges(7,2);
t21 = t214 * t444 / 0.2e1 + ((t383 * t293 + t382 * t296) * t475 + t231 * t451 + t232 * t452 + (t426 / 0.2e1 + t425 / 0.2e1 + t427 / 0.2e1 + t428 / 0.2e1) * t294 + t494 * t454 + t504 * t450) * t294;
t437 = t16 * qJD(4) - t21 * qJD(5);
t436 = m(7) * qJD(5);
t435 = m(7) * qJD(6);
t431 = Ifges(7,5) * t296;
t418 = t297 * Ifges(6,6);
t417 = t297 * t83;
t416 = t82 * t293;
t412 = -t339 - mrSges(5,1);
t398 = t187 * t297;
t24 = 0.2e1 * t484 * (t334 + t398) * t294;
t407 = qJD(1) * t24;
t405 = t139 * t293;
t399 = t187 * t294 ^ 2;
t397 = t21 * qJD(2);
t388 = t294 * t139;
t381 = t502 * t139;
t380 = t502 * t297;
t377 = qJD(4) * t297;
t376 = qJD(5) * t294;
t375 = -0.1e1 + t508;
t373 = m(7) * t462;
t368 = t139 * t387;
t364 = -t387 / 0.2e1;
t272 = Ifges(7,5) * t387;
t189 = -Ifges(7,6) * t297 + Ifges(7,3) * t390 + t272;
t191 = t294 * t490 - t418;
t362 = t189 / 0.2e1 - t191 / 0.2e1;
t247 = Ifges(6,2) * t296 + t432;
t353 = -t485 / 0.2e1 - t247 / 0.2e1;
t249 = Ifges(7,1) * t293 - t431;
t352 = t249 / 0.2e1 + t251 / 0.2e1;
t342 = qJD(4) * (-t338 + t412);
t246 = Ifges(7,3) * t293 + t431;
t34 = (t95 / 0.4e1 + t368 / 0.4e1 + t417 / 0.4e1) * t480;
t46 = t215 * t387 - m(7) * (-t124 * t297 - t153 * t387) + t297 * t236;
t333 = qJD(1) * t34 + qJD(2) * t46;
t12 = m(5) * (-t140 * t297 + t188 - t388) * t187 + t500 * (-t187 * t388 + t82 * t95 + t83 * t96);
t332 = t12 * qJD(1) + t24 * qJD(3);
t13 = t500 * (t505 * t139 - t405 * t82);
t14 = -0.2e1 * t484 * t375 * t388 + 0.2e1 * (t416 - t505) * t483;
t331 = t13 * qJD(1) + t14 * qJD(3);
t145 = t493 * t293;
t305 = (-t293 * t153 + (-t396 - t442) * t296) * t474 - t293 * t215 / 0.2e1;
t39 = t374 + (-mrSges(7,1) / 0.2e1 - t338 * t451) * t294 + t373 - t305;
t327 = qJD(2) * t39 + qJD(4) * t145;
t55 = -t282 + (t391 / 0.4e1 - t137 / 0.4e1 + (t393 / 0.4e1 + t464) * t297) * t480;
t326 = qJD(2) * t55 + qJD(5) * t281;
t318 = t139 * t339;
t317 = t139 * t338;
t316 = (-pkin(5) * t95 + qJ(6) * t96) * t474;
t2 = (-t419 / 0.2e1 + t455 + (-mrSges(5,1) / 0.2e1 - t339 / 0.2e1 + t459) * t294) * t187 + t481 + t506;
t190 = Ifges(7,6) * t294 + t246 * t297;
t192 = Ifges(6,6) * t294 + t297 * t490;
t194 = Ifges(7,4) * t294 + t297 * t492;
t196 = Ifges(6,5) * t294 + t252 * t297;
t7 = t124 * t235 + t127 * t234 + t136 * t233 + t137 * t230 + t138 * t236 + t141 * t232 + t149 * t231 + t150 * t229 + t153 * t217 + t154 * t215 + t277 * t244 + m(6) * (t136 * t149 + t137 * t150) + m(7) * (t124 * t138 + t127 * t141 + t153 * t154) + (t276 * t216 + Ifges(5,4) * t297 + t503 * t296 + (-t297 * t369 + t362) * t293) * t297 + (-Ifges(5,4) * t294 + t276 * t218 + (t194 / 0.2e1 + t196 / 0.2e1 + t370 * t294) * t296 + (t190 / 0.2e1 - t192 / 0.2e1 + t369 * t294) * t293 + (m(6) * t276 ^ 2 + Ifges(5,1) - Ifges(5,2) - t495) * t297) * t294;
t314 = t2 * qJD(1) + t7 * qJD(2) + t16 * qJD(3);
t308 = m(7) * (t139 * t214 + t382 * t83 + t383 * t82);
t3 = -t371 * t96 - t372 * t95 + t358 * t83 + t355 * t82 + t316 - t308 / 0.2e1 + (-t317 / 0.2e1 - t318 / 0.2e1 + t497 * (t416 / 0.2e1 + t415 / 0.2e1)) * t294;
t219 = t485 * t294;
t220 = t294 * t247;
t221 = -Ifges(7,1) * t390 + t272;
t222 = t294 * t251;
t271 = Ifges(7,6) * t387;
t9 = t214 * t215 + m(7) * (t124 * t136 + t127 * t137 + t153 * t214) + t136 * t236 + t136 * t229 - t137 * t231 + t137 * t232 + t271 * t450 + ((t153 * mrSges(7,1) - t124 * mrSges(7,2) + t418 / 0.2e1 - t137 * mrSges(6,3) - t222 / 0.2e1 + t221 / 0.2e1 + mrSges(6,1) * t394 + t362) * t296 + (t153 * mrSges(7,3) + t136 * mrSges(6,3) - t127 * mrSges(7,2) + t219 / 0.2e1 + t220 / 0.2e1 - mrSges(6,2) * t394 - t503) * t293) * t294;
t313 = -t3 * qJD(1) + t9 * qJD(2) - t21 * qJD(3);
t81 = 0.4e1 * t375 * t294 * t483;
t311 = t14 * qJD(1) + t16 * qJD(2) + t81 * qJD(3);
t10 = ((t439 / 0.4e1 - t386 / 0.4e1 - t241 / 0.4e1) * t480 - t482) * t139;
t30 = -pkin(4) * t243 - t241 * t338 + (t242 + t445) * t237 + (-t246 / 0.2e1 + t490 / 0.2e1 + t352) * t296 + (t492 / 0.2e1 + t252 / 0.2e1 + t353) * t293;
t50 = ((-t439 / 0.2e1 + t386 / 0.2e1 + t458) * m(7) + t482) * t297;
t298 = t276 * t456 + (t252 / 0.4e1 + t492 / 0.4e1 - t247 / 0.4e1 - t485 / 0.4e1 + pkin(4) * t472 + mrSges(7,1) * t460) * t296 + (-t251 / 0.4e1 - t249 / 0.4e1 - t490 / 0.4e1 + t246 / 0.4e1 + pkin(4) * t471 + mrSges(7,3) * t460) * t293 + t497 * pkin(9) * (-t290 / 0.2e1 - t289 / 0.2e1);
t300 = (-pkin(5) * t141 + qJ(6) * t138) * t475 + pkin(5) * t461 + t235 * t464 + t138 * t470 + mrSges(7,1) * t462 + t149 * t472 + t150 * t471;
t302 = (t153 * t241 + t214 * t237) * t474 + t153 * t457 + t214 * t459 + t215 * t458 - t491 * t297 / 0.4e1;
t303 = (t463 - t124 / 0.2e1) * mrSges(7,2) + (t383 * t474 - t355) * pkin(9) + t189 / 0.4e1 - t191 / 0.4e1 + t221 / 0.4e1 - t222 / 0.4e1;
t304 = -t220 / 0.4e1 - t219 / 0.4e1 + t195 / 0.4e1 + t193 / 0.4e1 + (t136 / 0.2e1 + t127 / 0.2e1) * mrSges(7,2) + (t382 * t474 - t358) * pkin(9);
t6 = (-t343 + t304) * t296 + ((0.3e1 / 0.4e1 * Ifges(6,6) - Ifges(7,6) / 0.2e1) * t297 + t303) * t293 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1 + t298) * t294 + t300 + t302;
t309 = t10 * qJD(1) - t6 * qJD(2) + t50 * qJD(3) - t30 * qJD(4);
t307 = -t424 / 0.2e1 - t423 / 0.2e1 - t421 / 0.2e1 + t420 / 0.2e1;
t238 = (m(7) * pkin(9) + mrSges(7,2)) * t296;
t142 = t276 * t399;
t69 = m(7) * t405;
t53 = (t391 + (-0.2e1 * qJ(6) + t393) * t297) * t474 + m(7) * t463 + t236;
t51 = -t241 * t444 / 0.2e1 + (-t445 / 0.2e1 + t307) * t297 + (t243 + t242) * t450;
t43 = -t338 * t364 + t373 - t422 / 0.2e1 + t305;
t35 = (-t368 - t417 + t95) * t474;
t11 = (t445 - t307 + t507) * t139;
t8 = qJD(2) * t24 + qJD(4) * t14;
t5 = t298 * t294 + (t418 / 0.4e1 + t303) * t293 + t304 * t296 - t300 + t302 + t495 * t453 + t369 * t389 + t496 * t384 / 0.2e1;
t4 = -t83 * t231 / 0.2e1 + t232 * t466 + t308 / 0.2e1 + t316 + (-mrSges(6,2) / 0.2e1 + mrSges(7,3) / 0.2e1) * t96 + t494 * t467 - t499 * t95 / 0.2e1 + (t317 + t318) * t453 + t497 * (t83 * t364 + t390 * t467);
t1 = mrSges(5,2) * t398 / 0.2e1 + t187 * t455 + t481 + (mrSges(5,1) / 0.2e1 - t504 / 0.2e1) * t400 - t506;
t15 = [qJD(2) * t12 + qJD(4) * t13 (-mrSges(3,2) * t341 - mrSges(3,1) * t348 - m(6) * t142 + m(5) * (-t187 * t392 + t188 * t277 - t142) + t188 * (-t297 * mrSges(5,1) + t294 * mrSges(5,2)) + m(4) * (-t187 * t292 - t188 * t410) * pkin(2) + t187 * mrSges(4,2) - t188 * mrSges(4,1) + (m(6) * t137 + m(7) * t124 + t494) * t96 + (-m(6) * t136 + m(7) * t127 - t231 + t232) * t95 + (-m(7) * t153 - t215 - t216) * t400 + (-t187 * t291 - t399) * mrSges(5,3)) * qJD(2) + t1 * qJD(4) + t4 * qJD(5) + t35 * qJD(6) + t332, t8, t1 * qJD(2) + t11 * qJD(5) - t69 * qJD(6) + t140 * t342 + ((-pkin(4) * t140 - t381) * t477 + (t140 * t237 - t381) * t474) * t479 + (-t497 * t508 + mrSges(5,2)) * qJD(4) * t139 + t331, t4 * qJD(2) + t11 * qJD(4) + (t498 * t82 - t499 * t83) * qJD(5) + ((-pkin(5) * t83 - qJ(6) * t82) * qJD(5) / 0.2e1 + qJD(6) * t466) * t480, t35 * qJD(2) - t69 * qJD(4) + t436 * t83; qJD(4) * t2 - qJD(5) * t3 - qJD(6) * t34 - t332, qJD(4) * t7 + qJD(5) * t9 - qJD(6) * t46, -t407 + t437, t5 * qJD(5) + t43 * qJD(6) + (Ifges(5,5) + t352 * t296 + t353 * t293 + (-m(6) * pkin(4) + t412) * t276) * t377 + t314 + (mrSges(5,2) * t394 + t190 * t452 + t192 * t451 - Ifges(5,6) * t294 + t237 * t217 - pkin(4) * t218 + ((t230 + t235) * t296 + (-t233 + t234) * t293 + t319 + m(7) * t329) * pkin(9) + (t194 + t196) * t454 + ((Ifges(6,6) - Ifges(7,6)) * t296 + t496 * t293) * t453 + t493 * t154 + t486 * mrSges(6,3) + t329 * mrSges(7,2)) * qJD(4), t5 * qJD(4) + (t271 + (-m(7) * pkin(5) - t499) * t137 + (-mrSges(6,2) + t281) * t136) * qJD(5) + t53 * qJD(6) + ((-qJ(6) * mrSges(7,2) - Ifges(6,6)) * t296 + (pkin(5) * mrSges(7,2) - t496) * t293) * t376 + t313, qJD(4) * t43 + qJD(5) * t53 - t333; t8, t407 + t437, t81 * qJD(4), t51 * qJD(5) + t294 * t342 + ((t380 - t438) * t477 + (t380 + t396) * t474) * t479 + ((-mrSges(5,2) + t425 + t426 + t427 + t428) * qJD(4) + t293 * t435) * t297 + t311, -t397 + t51 * qJD(4) - t214 * t436 + (t498 * qJD(5) * t293 + (-t499 * qJD(5) + t435) * t296) * t294 (t293 * t377 + t296 * t376) * m(7); -qJD(2) * t2 - qJD(5) * t10 - t331, qJD(5) * t6 - qJD(6) * t39 - t314, -qJD(5) * t50 - t311, qJD(5) * t30 - qJD(6) * t145, t238 * qJD(6) - t309 + (-Ifges(6,6) * t293 + (-m(7) * t335 + t504) * pkin(9) - t335 * mrSges(7,2) + t491) * qJD(5), qJD(5) * t238 - t327; qJD(2) * t3 + qJD(4) * t10, -qJD(4) * t6 + qJD(6) * t55 - t313, qJD(4) * t50 + t397, t309, t281 * qJD(6), t326; qJD(2) * t34, qJD(4) * t39 - qJD(5) * t55 + t333, 0, t327, -t326, 0;];
Cq  = t15;
