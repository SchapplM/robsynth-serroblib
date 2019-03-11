% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:14
% EndTime: 2019-03-08 23:00:51
% DurationCPUTime: 21.32s
% Computational Cost: add. (13222->730), mult. (30805->1003), div. (0->0), fcn. (23955->18), ass. (0->346)
t332 = sin(qJ(3));
t338 = -pkin(9) - pkin(8);
t295 = t338 * t332;
t336 = cos(qJ(3));
t296 = t338 * t336;
t331 = sin(qJ(4));
t335 = cos(qJ(4));
t219 = t331 * t295 - t335 * t296;
t275 = t331 * t336 + t332 * t335;
t398 = qJD(3) * t338;
t282 = t332 * t398;
t283 = t336 * t398;
t337 = cos(qJ(2));
t327 = sin(pkin(6));
t425 = qJD(1) * t327;
t396 = t337 * t425;
t512 = -qJD(4) * t219 + t275 * t396 - t282 * t331 + t335 * t283;
t274 = -t331 * t332 + t335 * t336;
t416 = qJD(4) * t335;
t417 = qJD(4) * t331;
t511 = -t274 * t396 + t335 * t282 + t331 * t283 + t295 * t416 + t296 * t417;
t352 = t275 * qJD(4);
t205 = -qJD(3) * t275 - t352;
t545 = -qJ(5) * t205 - qJD(5) * t274 - t511;
t351 = t274 * qJD(4);
t204 = qJD(3) * t274 + t351;
t544 = -qJ(5) * t204 - qJD(5) * t275 + t512;
t543 = -mrSges(7,3) + mrSges(6,2);
t330 = sin(qJ(6));
t334 = cos(qJ(6));
t375 = -mrSges(7,1) * t334 + mrSges(7,2) * t330;
t541 = t375 - mrSges(6,1);
t325 = sin(pkin(12));
t328 = cos(pkin(12));
t542 = -t544 * t325 + t328 * t545;
t128 = t204 * t325 - t328 * t205;
t129 = t204 * t328 + t205 * t325;
t419 = qJD(3) * t332;
t315 = pkin(3) * t419;
t181 = -pkin(4) * t205 + t315;
t333 = sin(qJ(2));
t397 = t333 * t425;
t540 = pkin(5) * t128 - pkin(10) * t129 + t181 - t397;
t323 = qJD(3) + qJD(4);
t286 = qJD(2) * pkin(8) + t397;
t380 = pkin(9) * qJD(2) + t286;
t329 = cos(pkin(6));
t424 = qJD(1) * t329;
t395 = t332 * t424;
t211 = t336 * t380 + t395;
t201 = t335 * t211;
t304 = t336 * t424;
t210 = -t332 * t380 + t304;
t202 = qJD(3) * pkin(3) + t210;
t131 = t202 * t331 + t201;
t266 = t274 * qJD(2);
t454 = qJ(5) * t266;
t114 = t131 + t454;
t101 = t328 * t114;
t199 = t331 * t211;
t130 = t335 * t202 - t199;
t267 = t275 * qJD(2);
t254 = t267 * qJ(5);
t113 = t130 - t254;
t99 = pkin(4) * t323 + t113;
t50 = t325 * t99 + t101;
t45 = pkin(10) * t323 + t50;
t481 = pkin(3) * t336;
t312 = pkin(2) + t481;
t256 = -qJD(2) * t312 - t396;
t187 = -pkin(4) * t266 + qJD(5) + t256;
t366 = t266 * t325 + t328 * t267;
t379 = t328 * t266 - t267 * t325;
t75 = -pkin(5) * t379 - pkin(10) * t366 + t187;
t19 = -t330 * t45 + t334 * t75;
t539 = t19 * mrSges(7,1);
t20 = t330 * t75 + t334 * t45;
t538 = t20 * mrSges(7,2);
t517 = t325 * t545 + t544 * t328;
t152 = t323 * t334 - t330 * t366;
t153 = t323 * t330 + t334 * t366;
t456 = mrSges(6,1) * t323 + mrSges(7,1) * t152 - mrSges(7,2) * t153 - mrSges(6,3) * t366;
t455 = cos(pkin(11));
t382 = t455 * t337;
t326 = sin(pkin(11));
t441 = t326 * t333;
t260 = -t329 * t441 + t382;
t324 = qJ(3) + qJ(4);
t319 = sin(t324);
t320 = cos(t324);
t442 = t326 * t327;
t537 = -t260 * t319 + t320 * t442;
t439 = t327 * t333;
t536 = -t319 * t439 + t320 * t329;
t413 = qJD(2) * qJD(3);
t284 = qJDD(2) * t336 - t332 * t413;
t285 = qJDD(2) * t332 + t336 * t413;
t155 = -qJD(2) * t352 + t284 * t335 - t285 * t331;
t423 = qJD(2) * t327;
t390 = qJD(1) * t423;
t297 = t333 * t390;
t412 = qJDD(1) * t327;
t246 = t337 * t412 - t297;
t234 = -qJDD(2) * pkin(2) - t246;
t195 = -pkin(3) * t284 + t234;
t112 = -pkin(4) * t155 + qJDD(5) + t195;
t154 = qJD(2) * t351 + t284 * t331 + t285 * t335;
t88 = -t154 * t325 + t155 * t328;
t89 = t154 * t328 + t155 * t325;
t22 = -pkin(5) * t88 - pkin(10) * t89 + t112;
t321 = qJDD(3) + qJDD(4);
t227 = t286 * t336 + t395;
t298 = t337 * t390;
t247 = t333 * t412 + t298;
t235 = qJDD(2) * pkin(8) + t247;
t411 = qJDD(1) * t329;
t141 = -t227 * qJD(3) - t235 * t332 + t336 * t411;
t126 = qJDD(3) * pkin(3) - pkin(9) * t285 + t141;
t140 = qJD(3) * t304 + t336 * t235 - t286 * t419 + t332 * t411;
t127 = pkin(9) * t284 + t140;
t35 = -qJD(4) * t131 + t335 * t126 - t127 * t331;
t23 = pkin(4) * t321 - qJ(5) * t154 - qJD(5) * t267 + t35;
t34 = t331 * t126 + t335 * t127 + t202 * t416 - t211 * t417;
t25 = qJ(5) * t155 + qJD(5) * t266 + t34;
t9 = t325 * t23 + t328 * t25;
t6 = pkin(10) * t321 + t9;
t2 = qJD(6) * t19 + t22 * t330 + t334 * t6;
t3 = -qJD(6) * t20 + t22 * t334 - t330 * t6;
t535 = t2 * t334 - t3 * t330;
t374 = mrSges(7,1) * t330 + mrSges(7,2) * t334;
t501 = -m(4) * pkin(8) + m(5) * t338 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - t374;
t318 = pkin(12) + t324;
t306 = sin(t318);
t307 = cos(t318);
t377 = -mrSges(4,1) * t336 + mrSges(4,2) * t332;
t534 = -m(4) * pkin(2) - m(5) * t312 - t320 * mrSges(5,1) + t319 * mrSges(5,2) - mrSges(3,1) + t377 + (-m(7) * pkin(5) + t541) * t307 + (-m(7) * pkin(10) + t543) * t306;
t63 = qJD(6) * t152 + t321 * t330 + t334 * t89;
t496 = t63 / 0.2e1;
t64 = -qJD(6) * t153 + t321 * t334 - t330 * t89;
t495 = t64 / 0.2e1;
t86 = qJDD(6) - t88;
t493 = t86 / 0.2e1;
t533 = m(6) + m(7);
t532 = t284 / 0.2e1;
t531 = t285 / 0.2e1;
t486 = t321 / 0.2e1;
t530 = -t366 / 0.2e1;
t529 = t366 / 0.2e1;
t528 = -t379 / 0.2e1;
t16 = -mrSges(7,1) * t64 + mrSges(7,2) * t63;
t77 = mrSges(6,1) * t321 - mrSges(6,3) * t89;
t527 = t16 - t77;
t526 = Ifges(6,4) * t366;
t525 = Ifges(6,4) * t379;
t524 = Ifges(5,5) * t275;
t523 = Ifges(5,6) * t274;
t522 = t336 * Ifges(4,2);
t453 = t114 * t325;
t49 = t328 * t99 - t453;
t44 = -pkin(5) * t323 - t49;
t521 = t374 * t44;
t197 = -t328 * t274 + t275 * t325;
t198 = t274 * t325 + t275 * t328;
t238 = -pkin(4) * t274 - t312;
t100 = pkin(5) * t197 - pkin(10) * t198 + t238;
t166 = qJ(5) * t274 + t219;
t218 = t335 * t295 + t296 * t331;
t357 = -qJ(5) * t275 + t218;
t104 = t328 * t166 + t325 * t357;
t48 = t100 * t330 + t104 * t334;
t519 = -qJD(6) * t48 + t330 * t542 + t334 * t540;
t47 = t100 * t334 - t104 * t330;
t518 = qJD(6) * t47 + t330 * t540 - t334 * t542;
t498 = m(5) * pkin(3);
t514 = -t498 - mrSges(4,1);
t156 = -mrSges(6,2) * t323 + mrSges(6,3) * t379;
t163 = qJD(6) - t379;
t96 = -mrSges(7,2) * t163 + mrSges(7,3) * t152;
t97 = mrSges(7,1) * t163 - mrSges(7,3) * t153;
t367 = -t330 * t97 + t334 * t96;
t513 = -t156 - t367;
t414 = qJD(6) * t334;
t356 = t129 * t330 + t198 * t414;
t510 = t140 * t336 - t141 * t332;
t509 = 0.2e1 * t486;
t120 = -mrSges(6,1) * t379 + mrSges(6,2) * t366;
t193 = -mrSges(5,1) * t266 + mrSges(5,2) * t267;
t506 = t377 * qJD(2) + t120 + t193;
t230 = -t306 * t439 + t307 * t329;
t231 = t306 * t329 + t307 * t439;
t505 = -t536 * mrSges(5,1) - (-t319 * t329 - t320 * t439) * mrSges(5,2) + t543 * t231 + t541 * t230;
t184 = -t260 * t306 + t307 * t442;
t185 = t260 * t307 + t306 * t442;
t504 = -t537 * mrSges(5,1) - (-t260 * t320 - t319 * t442) * mrSges(5,2) + t543 * t185 + t541 * t184;
t383 = t455 * t333;
t440 = t326 * t337;
t258 = t329 * t383 + t440;
t384 = t327 * t455;
t182 = -t258 * t306 - t307 * t384;
t183 = t258 * t307 - t306 * t384;
t347 = -t258 * t319 - t320 * t384;
t503 = -t347 * mrSges(5,1) - (-t258 * t320 + t319 * t384) * mrSges(5,2) + t543 * t183 + t541 * t182;
t502 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t30 = mrSges(7,1) * t86 - mrSges(7,3) * t63;
t31 = -mrSges(7,2) * t86 + mrSges(7,3) * t64;
t415 = qJD(6) * t330;
t499 = m(7) * ((-t19 * t334 - t20 * t330) * qJD(6) + t535) - t97 * t414 - t96 * t415 + t31 * t334 - t30 * t330;
t339 = qJD(2) ^ 2;
t497 = Ifges(7,1) * t496 + Ifges(7,4) * t495 + Ifges(7,5) * t493;
t492 = -t152 / 0.2e1;
t491 = -t153 / 0.2e1;
t490 = t153 / 0.2e1;
t489 = -t163 / 0.2e1;
t487 = t267 / 0.2e1;
t485 = -t323 / 0.2e1;
t483 = t334 / 0.2e1;
t482 = pkin(3) * t335;
t480 = pkin(4) * t267;
t479 = pkin(4) * t325;
t478 = pkin(4) * t328;
t474 = t49 * mrSges(6,3);
t473 = t50 * mrSges(6,3);
t469 = Ifges(4,4) * t332;
t468 = Ifges(4,4) * t336;
t467 = Ifges(7,4) * t330;
t466 = Ifges(7,4) * t334;
t465 = pkin(3) * qJD(4);
t464 = t130 * mrSges(5,3);
t463 = t131 * mrSges(5,3);
t462 = t153 * Ifges(7,4);
t461 = t267 * Ifges(5,4);
t449 = t379 * t330;
t448 = t379 * t334;
t447 = t198 * t330;
t446 = t198 * t334;
t443 = t325 * t331;
t438 = t327 * t336;
t437 = t327 * t337;
t436 = t328 * t331;
t133 = t335 * t210 - t199;
t288 = -pkin(3) * t332 - pkin(4) * t319;
t289 = pkin(4) * t320 + t481;
t429 = t260 * t288 + t289 * t442;
t426 = t288 * t439 + t329 * t289;
t311 = pkin(4) + t482;
t252 = pkin(3) * t436 + t325 * t311;
t422 = qJD(2) * t332;
t421 = qJD(2) * t333;
t420 = qJD(2) * t336;
t418 = qJD(3) * t336;
t410 = Ifges(7,5) * t63 + Ifges(7,6) * t64 + Ifges(7,3) * t86;
t314 = pkin(3) * t422;
t407 = mrSges(4,3) * t422;
t406 = mrSges(4,3) * t420;
t401 = t330 * t437;
t400 = t334 * t437;
t149 = Ifges(7,4) * t152;
t69 = t153 * Ifges(7,1) + t163 * Ifges(7,5) + t149;
t399 = t69 * t483;
t394 = t327 * t421;
t393 = t337 * t423;
t40 = -t88 * mrSges(6,1) + t89 * mrSges(6,2);
t388 = -t415 / 0.2e1;
t387 = t182 * pkin(5) + t183 * pkin(10);
t386 = t184 * pkin(5) + pkin(10) * t185;
t385 = t230 * pkin(5) + pkin(10) * t231;
t381 = t413 / 0.2e1;
t132 = -t210 * t331 - t201;
t378 = t537 * pkin(4);
t373 = Ifges(7,1) * t334 - t467;
t372 = t469 + t522;
t371 = -Ifges(7,2) * t330 + t466;
t370 = Ifges(4,5) * t336 - Ifges(4,6) * t332;
t369 = Ifges(7,5) * t334 - Ifges(7,6) * t330;
t368 = -t19 * t330 + t20 * t334;
t8 = t23 * t328 - t25 * t325;
t262 = t329 * t336 - t332 * t439;
t263 = t329 * t332 + t333 * t438;
t167 = t262 * t335 - t263 * t331;
t168 = t262 * t331 + t263 * t335;
t365 = t536 * pkin(4);
t251 = -pkin(3) * t443 + t311 * t328;
t361 = t258 * t288 - t289 * t384;
t110 = t167 * t325 + t168 * t328;
t93 = -t110 * t330 - t400;
t360 = -t110 * t334 + t401;
t358 = t132 - t454;
t355 = -t129 * t334 + t198 * t415;
t287 = -qJD(2) * pkin(2) - t396;
t354 = t287 * (mrSges(4,1) * t332 + mrSges(4,2) * t336);
t353 = t332 * (Ifges(4,1) * t336 - t469);
t95 = pkin(5) * t366 - pkin(10) * t379 + t480;
t257 = -t329 * t382 + t441;
t259 = t329 * t440 + t383;
t345 = -g(1) * t259 - g(2) * t257 + g(3) * t437;
t344 = t347 * pkin(4);
t106 = Ifges(6,2) * t379 + t323 * Ifges(6,6) + t526;
t107 = Ifges(6,1) * t366 + t323 * Ifges(6,5) + t525;
t14 = t63 * Ifges(7,4) + t64 * Ifges(7,2) + t86 * Ifges(7,6);
t164 = t266 * Ifges(5,2) + t323 * Ifges(5,6) + t461;
t255 = Ifges(5,4) * t266;
t165 = t267 * Ifges(5,1) + t323 * Ifges(5,5) + t255;
t5 = -pkin(5) * t321 - t8;
t67 = t153 * Ifges(7,5) + t152 * Ifges(7,6) + t163 * Ifges(7,3);
t68 = t152 * Ifges(7,2) + t163 * Ifges(7,6) + t462;
t340 = (Ifges(5,3) + Ifges(6,3)) * t321 + (t449 / 0.2e1 + t388) * t68 + (t521 + t399) * qJD(6) - t256 * (mrSges(5,1) * t267 + mrSges(5,2) * t266) + (t152 * t371 + t153 * t373 + t163 * t369) * qJD(6) / 0.2e1 - (-Ifges(5,2) * t267 + t165 + t255) * t266 / 0.2e1 + t5 * t375 - t69 * t448 / 0.2e1 + (-t187 * mrSges(6,2) + Ifges(6,1) * t530 + Ifges(6,5) * t485 + t369 * t489 + t371 * t492 + t373 * t491 + t474 - t521) * t379 + t106 * t529 + (-t187 * mrSges(6,1) + Ifges(7,5) * t491 - Ifges(6,2) * t528 - Ifges(6,6) * t485 + Ifges(7,6) * t492 + Ifges(7,3) * t489 + t473 + t538 - t539) * t366 + Ifges(5,5) * t154 + Ifges(5,6) * t155 + Ifges(6,6) * t88 + Ifges(6,5) * t89 - t34 * mrSges(5,2) + t35 * mrSges(5,1) + (-t526 + t67) * t530 + t8 * mrSges(6,1) - t9 * mrSges(6,2) + t267 * t463 + t266 * t464 + t14 * t483 + (Ifges(5,5) * t266 - Ifges(5,6) * t267) * t485 + t164 * t487 + (Ifges(7,5) * t330 + Ifges(7,6) * t334) * t493 + (Ifges(7,2) * t334 + t467) * t495 + (Ifges(7,1) * t330 + t466) * t496 + t330 * t497 + ((-t415 + t449) * t20 + (-t414 + t448) * t19 + t535) * mrSges(7,3) + (t525 + t107) * t528 - t267 * (Ifges(5,1) * t266 - t461) / 0.2e1;
t322 = -qJ(5) + t338;
t313 = Ifges(4,4) * t420;
t309 = -pkin(5) - t478;
t291 = -qJD(3) * mrSges(4,2) + t406;
t290 = qJD(3) * mrSges(4,1) - t407;
t281 = pkin(2) + t289;
t265 = Ifges(4,1) * t422 + Ifges(4,5) * qJD(3) + t313;
t264 = Ifges(4,6) * qJD(3) + qJD(2) * t372;
t250 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t285;
t249 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t284;
t241 = -pkin(5) - t251;
t229 = mrSges(5,1) * t323 - mrSges(5,3) * t267;
t228 = -mrSges(5,2) * t323 + mrSges(5,3) * t266;
t226 = -t286 * t332 + t304;
t225 = t314 + t480;
t214 = -mrSges(4,1) * t284 + mrSges(4,2) * t285;
t207 = qJD(3) * t262 + t336 * t393;
t206 = -qJD(3) * t263 - t332 * t393;
t145 = -mrSges(5,2) * t321 + mrSges(5,3) * t155;
t144 = mrSges(5,1) * t321 - mrSges(5,3) * t154;
t118 = -t254 + t133;
t109 = -t328 * t167 + t168 * t325;
t103 = t166 * t325 - t328 * t357;
t92 = -mrSges(5,1) * t155 + mrSges(5,2) * t154;
t90 = t314 + t95;
t82 = -qJD(4) * t168 + t206 * t335 - t207 * t331;
t81 = qJD(4) * t167 + t206 * t331 + t207 * t335;
t76 = -mrSges(6,2) * t321 + mrSges(6,3) * t88;
t59 = t328 * t118 + t325 * t358;
t56 = t113 * t328 - t453;
t55 = t113 * t325 + t101;
t37 = t325 * t82 + t328 * t81;
t36 = t325 * t81 - t328 * t82;
t29 = t330 * t90 + t334 * t59;
t28 = -t330 * t59 + t334 * t90;
t27 = t330 * t95 + t334 * t56;
t26 = -t330 * t56 + t334 * t95;
t18 = qJD(6) * t360 - t330 * t37 + t334 * t394;
t17 = qJD(6) * t93 + t330 * t394 + t334 * t37;
t1 = [m(2) * qJDD(1) + t110 * t76 + t167 * t144 + t168 * t145 + t37 * t156 + t17 * t96 + t18 * t97 + t206 * t290 + t207 * t291 + t81 * t228 + t82 * t229 + t263 * t249 + t262 * t250 + t93 * t30 - t360 * t31 - t456 * t36 + t527 * t109 + (-m(2) - m(3) - m(4) - m(5) - t533) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t339 - t214 - t40 - t92) * t337 + (-mrSges(3,1) * t339 - mrSges(3,2) * qJDD(2) + qJD(2) * t506) * t333) * t327 + m(3) * (qJDD(1) * t329 ^ 2 + (t246 * t337 + t247 * t333) * t327) + m(6) * (-t109 * t8 + t110 * t9 - t36 * t49 + t37 * t50 + (-t112 * t337 + t187 * t421) * t327) + m(5) * (t130 * t82 + t131 * t81 + t167 * t35 + t168 * t34 + (-t195 * t337 + t256 * t421) * t327) + m(4) * (t140 * t263 + t141 * t262 + t206 * t226 + t207 * t227 + (-t234 * t337 + t287 * t421) * t327) + m(7) * (t109 * t5 + t17 * t20 + t18 * t19 - t2 * t360 + t3 * t93 + t36 * t44); (Ifges(5,5) * t204 + Ifges(6,5) * t129 + Ifges(5,6) * t205 - Ifges(6,6) * t128) * t323 / 0.2e1 + t379 * (Ifges(6,4) * t129 - Ifges(6,2) * t128) / 0.2e1 + (t112 * mrSges(6,2) - t8 * mrSges(6,3) + Ifges(6,1) * t89 + Ifges(6,4) * t88 + Ifges(6,5) * t509 + t369 * t493 + t371 * t495 + t373 * t496 + t5 * t374 + t69 * t388) * t198 + (t410 / 0.2e1 - Ifges(6,4) * t89 - Ifges(6,2) * t88 + t112 * mrSges(6,1) + Ifges(7,3) * t493 + Ifges(7,6) * t495 + Ifges(7,5) * t496 - t9 * mrSges(6,3) - t509 * Ifges(6,6) + t502) * t197 + (-m(5) * t256 - m(6) * t187 - t506) * t397 + (t246 + t297) * mrSges(3,1) + (-t533 * (-t259 * t281 - t260 * t322) + t501 * t260 - t534 * t259) * g(1) + (-t533 * (-t257 * t281 - t258 * t322) + t501 * t258 - t534 * t257) * g(2) + qJDD(3) * (Ifges(4,5) * t332 + Ifges(4,6) * t336) + t234 * t377 + t265 * t418 / 0.2e1 - t264 * t419 / 0.2e1 + (-t103 * t8 + t104 * t9 + t112 * t238 + t181 * t187 + t49 * t517 - t50 * t542) * m(6) - t542 * t156 - t14 * t447 / 0.2e1 + (-t247 + t298) * mrSges(3,2) + t44 * (mrSges(7,1) * t356 - mrSges(7,2) * t355) + t163 * (-Ifges(7,5) * t355 - Ifges(7,6) * t356 + Ifges(7,3) * t128) / 0.2e1 + t152 * (-Ifges(7,4) * t355 - Ifges(7,2) * t356 + Ifges(7,6) * t128) / 0.2e1 - t128 * t473 - t129 * t474 + (-(t287 * t333 + (-t226 * t332 + t227 * t336) * t337) * t425 - pkin(2) * t234) * m(4) + t193 * t315 + (Ifges(6,1) * t129 - Ifges(6,4) * t128) * t529 + t468 * t531 + t372 * t532 - t128 * t538 - t312 * t92 + (t354 + t370 * qJD(3) / 0.2e1) * qJD(3) + t195 * (-mrSges(5,1) * t274 + mrSges(5,2) * t275) + t154 * (Ifges(5,1) * t275 + Ifges(5,4) * t274) + t155 * (Ifges(5,4) * t275 + Ifges(5,2) * t274) + Ifges(3,3) * qJDD(2) + t266 * (Ifges(5,4) * t204 + Ifges(5,2) * t205) / 0.2e1 + t256 * (-mrSges(5,1) * t205 + mrSges(5,2) * t204) + t238 * t40 + (t19 * t355 - t2 * t447 - t20 * t356 - t3 * t446) * mrSges(7,3) + t218 * t144 + t219 * t145 - pkin(2) * t214 + t204 * t165 / 0.2e1 + t205 * t164 / 0.2e1 + t187 * (mrSges(6,1) * t128 + mrSges(6,2) * t129) + t181 * t120 + t129 * t107 / 0.2e1 + t128 * t67 / 0.2e1 - t128 * t106 / 0.2e1 + t104 * t76 + t47 * t30 + t48 * t31 + (t274 * t34 - t275 * t35) * mrSges(5,3) + t517 * t456 + t205 * t463 + t128 * t539 + t353 * t381 + (-t290 * t418 - t291 * t419 + m(4) * ((-t226 * t336 - t227 * t332) * qJD(3) + t510) + t336 * t249 - t332 * t250) * pkin(8) + (-t226 * t418 - t227 * t419 + t510) * mrSges(4,3) - t356 * t68 / 0.2e1 + t511 * t228 + t512 * t229 + (t130 * t512 + t131 * t511 - t195 * t312 + t218 * t35 + t219 * t34 + t256 * t315) * m(5) + t518 * t96 + t519 * t97 + (t103 * t5 + t19 * t519 + t2 * t48 + t20 * t518 + t3 * t47 - t44 * t517) * m(7) + (t523 / 0.2e1 + t524 / 0.2e1) * t321 + (t523 + t524) * t486 + t527 * t103 + t129 * t399 + (Ifges(5,1) * t204 + Ifges(5,4) * t205) * t487 + (-Ifges(7,1) * t355 - Ifges(7,4) * t356 + Ifges(7,5) * t128) * t490 + t446 * t497 + (Ifges(4,4) * t531 + Ifges(4,2) * t532 - t291 * t396 + t468 * t381) * t336 + (Ifges(4,1) * t285 + Ifges(4,4) * t532 + t290 * t396 - t381 * t522) * t332 + (-t533 * t281 * t437 + (t534 * t337 + (t322 * t533 + t501) * t333) * t327) * g(3) - t204 * t464; -(-Ifges(4,2) * t422 + t265 + t313) * t420 / 0.2e1 + t499 * (pkin(10) + t252) + (m(6) * t50 + m(7) * t368 - t513) * (t328 * t335 - t443) * t465 + (-t19 * t28 - t20 * t29 + t241 * t5) * m(7) + (-t187 * t225 + t251 * t8 + t252 * t9 - t50 * t59) * m(6) - t370 * t413 / 0.2e1 - m(5) * (t130 * t132 + t131 * t133 + t256 * t314) + (t407 + t290) * t227 + (t331 * t145 + t228 * t416 - t229 * t417) * pkin(3) + t264 * t422 / 0.2e1 - qJD(2) * t354 + t340 - t339 * t353 / 0.2e1 + (t406 - t291) * t226 + Ifges(4,6) * t284 + Ifges(4,5) * t285 + Ifges(4,3) * qJDD(3) + t251 * t77 + t252 * t76 + t241 * t16 - t132 * t229 - t225 * t120 - t133 * t228 - t59 * t156 + t141 * mrSges(4,1) - t140 * mrSges(4,2) - t28 * t97 - t29 * t96 + (-m(6) * t49 + m(7) * t44 - t456) * (-t118 * t325 + t328 * t358 + (t325 * t335 + t436) * t465) + (-m(6) * t361 - m(7) * (t361 + t387) - (-t258 * t336 + t332 * t384) * mrSges(4,2) + t514 * (-t258 * t332 - t336 * t384) + t503) * g(2) + (-(-t260 * t336 - t332 * t442) * mrSges(4,2) - m(6) * t429 - m(7) * (t386 + t429) + t514 * (-t260 * t332 + t326 * t438) + t504) * g(1) + (-m(7) * (t385 + t426) - m(6) * t426 + mrSges(4,2) * t263 + t514 * t262 + t505) * g(3) + t144 * t482 + (t331 * t34 + t335 * t35 + (-t130 * t331 + t131 * t335) * qJD(4)) * t498 - t193 * t314; -t120 * t480 - t130 * t228 + t131 * t229 - t56 * t156 + t309 * t16 - t26 * t97 - t27 * t96 + t77 * t478 + t76 * t479 + t340 + t456 * t55 + (-t19 * t26 - t20 * t27 + t309 * t5 - t44 * t55) * m(7) + (-t187 * t480 + t49 * t55 - t50 * t56 + (t325 * t9 + t328 * t8) * pkin(4)) * m(6) + (-m(6) * t365 - m(7) * (t365 + t385) + t505) * g(3) + (-m(6) * t344 - m(7) * (t344 + t387) + t503) * g(2) + (-m(6) * t378 - m(7) * (t378 + t386) + t504) * g(1) + t499 * (pkin(10) + t479); t334 * t30 + t330 * t31 + t456 * t366 + t367 * qJD(6) + t513 * t379 + t40 + (t163 * t368 + t2 * t330 + t3 * t334 - t366 * t44 + t345) * m(7) + (t366 * t49 - t379 * t50 + t112 + t345) * m(6); -t44 * (mrSges(7,1) * t153 + mrSges(7,2) * t152) + (Ifges(7,1) * t152 - t462) * t491 + t68 * t490 + (Ifges(7,5) * t152 - Ifges(7,6) * t153) * t489 - t19 * t96 + t20 * t97 - g(1) * ((-t185 * t330 + t259 * t334) * mrSges(7,1) + (-t185 * t334 - t259 * t330) * mrSges(7,2)) - g(2) * ((-t183 * t330 + t257 * t334) * mrSges(7,1) + (-t183 * t334 - t257 * t330) * mrSges(7,2)) - g(3) * ((-t231 * t330 - t400) * mrSges(7,1) + (-t231 * t334 + t401) * mrSges(7,2)) + (t152 * t19 + t153 * t20) * mrSges(7,3) + t410 + (-Ifges(7,2) * t153 + t149 + t69) * t492 + t502;];
tau  = t1;
