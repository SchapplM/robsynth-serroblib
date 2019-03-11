% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:55
% EndTime: 2019-03-09 16:31:33
% DurationCPUTime: 24.80s
% Computational Cost: add. (16963->749), mult. (39455->942), div. (0->0), fcn. (29187->14), ass. (0->362)
t640 = Ifges(6,4) + Ifges(7,4);
t364 = cos(qJ(5));
t641 = mrSges(7,1) + mrSges(6,1);
t643 = t364 * t641;
t602 = Ifges(6,1) + Ifges(7,1);
t600 = Ifges(6,5) + Ifges(7,5);
t629 = Ifges(6,2) + Ifges(7,2);
t598 = Ifges(7,6) + Ifges(6,6);
t361 = sin(qJ(3));
t362 = sin(qJ(2));
t365 = cos(qJ(3));
t366 = cos(qJ(2));
t295 = -t361 * t362 + t365 * t366;
t274 = t295 * qJD(1);
t296 = t361 * t366 + t362 * t365;
t275 = t296 * qJD(1);
t357 = sin(pkin(10));
t358 = cos(pkin(10));
t415 = t358 * t274 - t275 * t357;
t605 = -t415 / 0.2e1;
t642 = Ifges(5,2) * t605;
t453 = qJD(1) * qJD(2);
t306 = qJDD(1) * t366 - t362 * t453;
t307 = qJDD(1) * t362 + t366 * t453;
t379 = t295 * qJD(3);
t207 = qJD(1) * t379 + t306 * t361 + t307 * t365;
t380 = t296 * qJD(3);
t208 = -qJD(1) * t380 + t306 * t365 - t307 * t361;
t139 = t207 * t358 + t208 * t357;
t355 = qJD(2) + qJD(3);
t360 = sin(qJ(5));
t390 = t274 * t357 + t358 * t275;
t203 = t355 * t364 - t360 * t390;
t353 = qJDD(2) + qJDD(3);
t89 = qJD(5) * t203 + t139 * t364 + t353 * t360;
t555 = t89 / 0.2e1;
t204 = t355 * t360 + t364 * t390;
t90 = -qJD(5) * t204 - t139 * t360 + t353 * t364;
t554 = t90 / 0.2e1;
t138 = -t207 * t357 + t208 * t358;
t133 = qJDD(5) - t138;
t553 = t133 / 0.2e1;
t634 = -mrSges(6,3) - mrSges(7,3);
t597 = Ifges(7,3) + Ifges(6,3);
t639 = t640 * t203;
t456 = qJD(5) * t360;
t492 = t415 * t360;
t638 = t456 - t492;
t637 = t640 * t204;
t539 = t355 / 0.2e1;
t592 = t390 * Ifges(5,4);
t636 = -Ifges(5,6) * t539 + t642 - t592 / 0.2e1;
t635 = t600 * t553 + t640 * t554 + t602 * t555;
t603 = mrSges(6,2) + mrSges(7,2);
t215 = qJD(5) - t415;
t624 = t203 * t629 + t215 * t598 + t637;
t623 = t204 * t602 + t600 * t215 + t639;
t633 = t640 * t364;
t632 = t640 * t360;
t356 = qJ(2) + qJ(3);
t347 = pkin(10) + t356;
t333 = sin(t347);
t630 = t333 * t603;
t628 = t133 * t598 + t629 * t90 + t640 * t89;
t626 = t203 * t598 + t204 * t600 + t215 * t597;
t478 = t357 * t361;
t503 = pkin(2) * qJD(3);
t265 = (t358 * t365 - t478) * t503;
t368 = -pkin(8) - pkin(7);
t324 = t368 * t366;
t303 = qJD(1) * t324;
t279 = t365 * t303;
t323 = t368 * t362;
t302 = qJD(1) * t323;
t234 = -t302 * t361 + t279;
t496 = qJ(4) * t274;
t205 = t234 - t496;
t276 = t361 * t303;
t235 = t365 * t302 + t276;
t268 = t275 * qJ(4);
t206 = -t268 + t235;
t137 = t205 * t357 + t206 * t358;
t533 = pkin(3) * t275;
t148 = pkin(4) * t390 - pkin(9) * t415 + t533;
t461 = qJD(1) * t362;
t343 = pkin(2) * t461;
t140 = t148 + t343;
t55 = t364 * t137 + t360 * t140;
t625 = t364 * t265 - t55;
t476 = t358 * t361;
t580 = -t358 * t205 + t206 * t357 - (t357 * t365 + t476) * t503;
t622 = qJ(6) * t492 + t364 * qJD(6);
t621 = -t360 * t598 + t364 * t600;
t620 = -t360 * t629 + t633;
t619 = t364 * t602 - t632;
t334 = cos(t347);
t349 = sin(t356);
t350 = cos(t356);
t618 = mrSges(4,1) * t349 + mrSges(5,1) * t333 + mrSges(4,2) * t350 + mrSges(5,2) * t334;
t615 = t333 * t643;
t614 = t638 * pkin(5);
t455 = qJD(5) * t364;
t491 = t415 * t364;
t613 = -t455 + t491;
t283 = qJD(2) * pkin(2) + t302;
t227 = t365 * t283 + t276;
t196 = t227 - t268;
t183 = pkin(3) * t355 + t196;
t228 = t283 * t361 - t279;
t197 = t228 + t496;
t477 = t358 * t197;
t115 = t357 * t183 + t477;
t111 = pkin(9) * t355 + t115;
t352 = t366 * pkin(2);
t340 = t352 + pkin(1);
t322 = t340 * qJD(1);
t240 = -pkin(3) * t274 + qJD(4) - t322;
t124 = -pkin(4) * t415 - pkin(9) * t390 + t240;
t290 = t307 * pkin(7);
t239 = qJDD(2) * pkin(2) - pkin(8) * t307 - t290;
t289 = t306 * pkin(7);
t241 = pkin(8) * t306 + t289;
t147 = -qJD(3) * t228 + t365 * t239 - t241 * t361;
t78 = pkin(3) * t353 - qJ(4) * t207 - qJD(4) * t275 + t147;
t457 = qJD(3) * t365;
t458 = qJD(3) * t361;
t146 = t361 * t239 + t365 * t241 + t283 * t457 + t303 * t458;
t88 = qJ(4) * t208 + qJD(4) * t274 + t146;
t27 = t357 * t78 + t358 * t88;
t24 = pkin(9) * t353 + t27;
t495 = qJDD(1) * pkin(1);
t269 = -pkin(2) * t306 - t495;
t171 = -pkin(3) * t208 + qJDD(4) + t269;
t38 = -pkin(4) * t138 - pkin(9) * t139 + t171;
t5 = -t111 * t456 + t124 * t455 + t364 * t24 + t360 * t38;
t47 = t111 * t364 + t124 * t360;
t6 = -qJD(5) * t47 - t24 * t360 + t364 * t38;
t611 = -t360 * t6 + t364 * t5;
t363 = sin(qJ(1));
t367 = cos(qJ(1));
t567 = g(1) * t367 + g(2) * t363;
t188 = t357 * t197;
t114 = t183 * t358 - t188;
t110 = -pkin(4) * t355 - t114;
t517 = mrSges(7,2) * t364;
t401 = mrSges(7,1) * t360 + t517;
t402 = mrSges(6,1) * t360 + mrSges(6,2) * t364;
t77 = -pkin(5) * t203 + qJD(6) + t110;
t610 = t110 * t402 + t401 * t77;
t609 = t240 * mrSges(5,2) - t114 * mrSges(5,3);
t608 = t350 * mrSges(4,1) + t334 * mrSges(5,1) - mrSges(4,2) * t349 + (-mrSges(5,2) - t634) * t333;
t46 = -t111 * t360 + t364 * t124;
t34 = -qJ(6) * t204 + t46;
t31 = pkin(5) * t215 + t34;
t35 = qJ(6) * t203 + t47;
t607 = -t240 * mrSges(5,1) - t46 * mrSges(6,1) - t31 * mrSges(7,1) + t47 * mrSges(6,2) + t35 * mrSges(7,2) + t115 * mrSges(5,3) - t636;
t556 = m(7) * pkin(5);
t606 = t306 / 0.2e1;
t541 = t353 / 0.2e1;
t536 = t366 / 0.2e1;
t545 = -t390 / 0.2e1;
t596 = Ifges(4,5) * t296;
t595 = Ifges(4,6) * t295;
t594 = t366 * Ifges(3,2);
t593 = t390 * Ifges(5,1);
t590 = t415 * Ifges(5,4);
t534 = pkin(2) * t365;
t339 = pkin(3) + t534;
t267 = pkin(2) * t476 + t357 * t339;
t262 = pkin(9) + t267;
t470 = -qJ(6) - t262;
t413 = qJD(5) * t470;
t588 = t360 * t413 + t622 + t625;
t351 = t364 * qJ(6);
t408 = pkin(5) * t390 - t351 * t415;
t54 = -t137 * t360 + t364 * t140;
t587 = (-qJD(6) - t265) * t360 + t364 * t413 - t408 - t54;
t531 = pkin(3) * t357;
t335 = pkin(9) + t531;
t469 = -qJ(6) - t335;
t412 = qJD(5) * t469;
t122 = t196 * t358 - t188;
t52 = -t122 * t360 + t364 * t148;
t586 = -qJD(6) * t360 + t364 * t412 - t408 - t52;
t53 = t364 * t122 + t360 * t148;
t585 = t360 * t412 - t53 + t622;
t584 = t556 + mrSges(7,1);
t121 = t196 * t357 + t477;
t583 = -t121 + t614;
t582 = t614 - t580;
t579 = -t137 + t265;
t578 = -mrSges(5,1) * t355 - mrSges(6,1) * t203 + mrSges(6,2) * t204 + mrSges(5,3) * t390;
t244 = t365 * t323 + t324 * t361;
t218 = -qJ(4) * t296 + t244;
t245 = t361 * t323 - t365 * t324;
t219 = qJ(4) * t295 + t245;
t159 = t218 * t357 + t219 * t358;
t156 = t364 * t159;
t230 = -t358 * t295 + t296 * t357;
t231 = t295 * t357 + t296 * t358;
t253 = -pkin(3) * t295 - t340;
t157 = pkin(4) * t230 - pkin(9) * t231 + t253;
t64 = t360 * t157 + t156;
t527 = pkin(5) * t364;
t338 = pkin(4) + t527;
t359 = -qJ(6) - pkin(9);
t389 = -t333 * t359 + t334 * t338;
t409 = t334 * pkin(4) + t333 * pkin(9);
t388 = -t333 * t338 - t334 * t359;
t529 = pkin(4) * t333;
t532 = pkin(3) * t349;
t576 = -m(7) * (t388 - t532) - m(6) * (-t529 - t532) + t615;
t575 = -m(7) * t388 + t615;
t572 = t133 * t597 + t598 * t90 + t600 * t89;
t460 = qJD(1) * t366;
t525 = pkin(7) * t366;
t526 = pkin(7) * t362;
t571 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t461) * t525 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t460) * t526;
t474 = t360 * t367;
t480 = t334 * t367;
t570 = -t474 * t630 + t480 * t634;
t475 = t360 * t363;
t482 = t334 * t363;
t569 = -t475 * t630 + t482 * t634;
t568 = t289 * t366 + t290 * t362;
t566 = m(5) + m(6) + m(7);
t565 = 0.2e1 * t541;
t564 = mrSges(6,1) + t584;
t562 = -m(3) * pkin(7) + m(4) * t368 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t561 = -t608 + (t360 * t603 - t643) * t334;
t514 = mrSges(6,3) * t203;
t151 = -mrSges(6,2) * t215 + t514;
t513 = mrSges(6,3) * t204;
t153 = mrSges(6,1) * t215 - t513;
t40 = mrSges(6,1) * t133 - mrSges(6,3) * t89;
t560 = m(6) * ((-t360 * t47 - t364 * t46) * qJD(5) + t611) - t153 * t455 - t151 * t456 - t360 * t40;
t320 = -mrSges(3,1) * t366 + mrSges(3,2) * t362;
t559 = m(3) * pkin(1) + m(4) * t340 + mrSges(2,1) - t320 + t608;
t3 = qJ(6) * t90 + qJD(6) * t203 + t5;
t558 = t6 * mrSges(6,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t551 = -t203 / 0.2e1;
t550 = t203 / 0.2e1;
t549 = -t204 / 0.2e1;
t548 = t204 / 0.2e1;
t547 = -t215 / 0.2e1;
t546 = t215 / 0.2e1;
t542 = t275 / 0.2e1;
t540 = -t355 / 0.2e1;
t535 = pkin(2) * t362;
t337 = pkin(3) * t350;
t530 = pkin(3) * t358;
t516 = mrSges(4,3) * t274;
t515 = mrSges(4,3) * t275;
t512 = mrSges(7,3) * t203;
t511 = mrSges(7,3) * t204;
t510 = Ifges(3,4) * t362;
t509 = Ifges(3,4) * t366;
t508 = Ifges(4,4) * t275;
t236 = qJD(2) * t295 + t379;
t237 = -qJD(2) * t296 - t380;
t173 = t236 * t358 + t237 * t357;
t494 = t173 * t360;
t493 = t173 * t364;
t490 = t231 * t360;
t489 = t231 * t364;
t488 = t262 * t364;
t479 = t335 * t364;
t473 = t363 * t364;
t471 = t364 * t367;
t462 = t337 + t352;
t459 = qJD(2) * t362;
t438 = qJD(2) * t368;
t304 = t362 * t438;
t305 = t366 * t438;
t179 = t365 * t304 + t361 * t305 + t323 * t457 + t324 * t458;
t134 = qJ(4) * t237 + qJD(4) * t295 + t179;
t180 = -qJD(3) * t245 - t304 * t361 + t365 * t305;
t135 = -qJ(4) * t236 - qJD(4) * t296 + t180;
t51 = t134 * t358 + t135 * t357;
t172 = t236 * t357 - t358 * t237;
t344 = pkin(2) * t459;
t224 = -pkin(3) * t237 + t344;
t75 = pkin(4) * t172 - pkin(9) * t173 + t224;
t449 = t157 * t455 + t360 * t75 + t364 * t51;
t439 = t337 + t409;
t435 = t231 * t455;
t28 = -t90 * mrSges(7,1) + t89 * mrSges(7,2);
t422 = -t456 / 0.2e1;
t26 = -t357 * t88 + t358 * t78;
t420 = -t360 * t51 + t364 * t75;
t419 = t453 / 0.2e1;
t418 = -t138 * mrSges(5,1) + t139 * mrSges(5,2);
t50 = t134 * t357 - t358 * t135;
t63 = t364 * t157 - t159 * t360;
t158 = -t358 * t218 + t219 * t357;
t266 = -pkin(2) * t478 + t339 * t358;
t261 = -pkin(4) - t266;
t407 = t337 + t389;
t406 = mrSges(3,1) * t362 + mrSges(3,2) * t366;
t398 = t510 + t594;
t395 = Ifges(3,5) * t366 - Ifges(3,6) * t362;
t391 = -t360 * t46 + t364 * t47;
t387 = -qJ(6) * t173 - qJD(6) * t231;
t23 = -pkin(4) * t353 - t26;
t384 = pkin(1) * t406;
t259 = -t334 * t474 + t473;
t257 = t334 * t475 + t471;
t383 = t435 + t494;
t382 = t231 * t456 - t493;
t381 = t362 * (Ifges(3,1) * t366 - t510);
t1 = pkin(5) * t133 - qJ(6) * t89 - qJD(6) * t204 + t6;
t12 = -pkin(5) * t90 + qJDD(6) + t23;
t162 = t355 * Ifges(5,5) + t590 + t593;
t216 = Ifges(4,2) * t274 + Ifges(4,6) * t355 + t508;
t270 = Ifges(4,4) * t274;
t217 = Ifges(4,1) * t275 + Ifges(4,5) * t355 + t270;
t369 = (-t1 * t360 + t3 * t364 + t31 * t613 - t35 * t638) * mrSges(7,3) + (t46 * t613 - t47 * t638 + t611) * mrSges(6,3) - t275 * (Ifges(4,1) * t274 - t508) / 0.2e1 + (Ifges(4,5) * t274 - Ifges(4,6) * t275) * t540 + t216 * t542 + (-t592 + t626) * t545 + t628 * t364 / 0.2e1 + (t360 * t600 + t364 * t598) * t553 + t228 * t515 + t227 * t516 + t322 * (mrSges(4,1) * t275 + mrSges(4,2) * t274) + (Ifges(4,3) + Ifges(5,3)) * t353 - (-Ifges(4,2) * t275 + t217 + t270) * t274 / 0.2e1 + (Ifges(5,1) * t545 + Ifges(5,5) * t540 + t547 * t621 + t549 * t619 + t551 * t620 - t609 - t610) * t415 + t610 * qJD(5) + (t203 * t620 + t204 * t619 + t215 * t621) * qJD(5) / 0.2e1 + t623 * (t455 / 0.2e1 - t491 / 0.2e1) + t624 * (t422 + t492 / 0.2e1) + (t162 + t590) * t605 + (t364 * t629 + t632) * t554 + (t360 * t602 + t633) * t555 + t360 * t635 - t146 * mrSges(4,2) + t147 * mrSges(4,1) + Ifges(5,5) * t139 + Ifges(5,6) * t138 + t26 * mrSges(5,1) - t27 * mrSges(5,2) + Ifges(4,5) * t207 + Ifges(4,6) * t208 + t12 * (-mrSges(7,1) * t364 + mrSges(7,2) * t360) + t23 * (-mrSges(6,1) * t364 + mrSges(6,2) * t360) + (-Ifges(5,6) * t540 + t547 * t597 + t549 * t600 + t551 * t598 + t607 - t642) * t390;
t354 = -qJ(4) + t368;
t342 = Ifges(3,4) * t460;
t336 = -pkin(4) - t530;
t317 = pkin(9) * t480;
t316 = pkin(9) * t482;
t315 = -t338 - t530;
t308 = -t532 - t535;
t301 = pkin(1) + t462;
t287 = t367 * t308;
t286 = t363 * t308;
t285 = t351 + t479;
t284 = t469 * t360;
t273 = Ifges(3,1) * t461 + Ifges(3,5) * qJD(2) + t342;
t272 = Ifges(3,6) * qJD(2) + qJD(1) * t398;
t260 = t334 * t471 + t475;
t258 = -t334 * t473 + t474;
t250 = t261 - t527;
t249 = mrSges(4,1) * t355 - t515;
t248 = -mrSges(4,2) * t355 + t516;
t247 = t343 + t533;
t243 = t351 + t488;
t242 = t470 * t360;
t226 = -mrSges(4,1) * t274 + mrSges(4,2) * t275;
t209 = -mrSges(5,2) * t355 + mrSges(5,3) * t415;
t185 = -mrSges(4,2) * t353 + mrSges(4,3) * t208;
t184 = mrSges(4,1) * t353 - mrSges(4,3) * t207;
t170 = -mrSges(5,1) * t415 + mrSges(5,2) * t390;
t152 = mrSges(7,1) * t215 - t511;
t150 = -mrSges(7,2) * t215 + t512;
t141 = -mrSges(7,1) * t203 + mrSges(7,2) * t204;
t119 = mrSges(5,1) * t353 - mrSges(5,3) * t139;
t118 = -mrSges(5,2) * t353 + mrSges(5,3) * t138;
t107 = pkin(5) * t490 + t158;
t48 = -qJ(6) * t490 + t64;
t45 = pkin(5) * t230 - t231 * t351 + t63;
t42 = -mrSges(6,2) * t133 + mrSges(6,3) * t90;
t41 = -mrSges(7,2) * t133 + mrSges(7,3) * t90;
t39 = mrSges(7,1) * t133 - mrSges(7,3) * t89;
t30 = pkin(5) * t383 + t50;
t29 = -mrSges(6,1) * t90 + mrSges(6,2) * t89;
t10 = -qJD(5) * t64 + t420;
t9 = -t159 * t456 + t449;
t8 = -qJ(6) * t435 + (-qJD(5) * t159 + t387) * t360 + t449;
t7 = pkin(5) * t172 + t387 * t364 + (-t156 + (qJ(6) * t231 - t157) * t360) * qJD(5) + t420;
t2 = [-t320 * t495 + (-t382 * t640 - t383 * t629) * t550 + (-t382 * t602 - t383 * t640) * t548 + (-t475 * t556 - t566 * (t367 * t301 - t354 * t363) - t641 * t260 - t603 * t259 + t562 * t363 + (-m(6) * t409 - m(7) * t389 - t559) * t367) * g(2) + (-t641 * t258 - t603 * t257 + (t354 * t566 - t360 * t556 + t562) * t367 + (-m(7) * (-t301 - t389) - m(6) * (-t301 - t409) + m(5) * t301 + t559) * t363) * g(1) + m(5) * (t115 * t51 + t159 * t27 + t171 * t253 + t224 * t240) + (Ifges(4,5) * t236 + Ifges(4,6) * t237) * t539 + (Ifges(4,1) * t236 + Ifges(4,4) * t237) * t542 - t628 * t490 / 0.2e1 + (t366 * t509 + t381) * t419 + (t572 / 0.2e1 + t171 * mrSges(5,1) + t1 * mrSges(7,1) - t27 * mrSges(5,3) - Ifges(5,4) * t139 - Ifges(5,2) * t138 - Ifges(5,6) * t565 + t553 * t597 + t554 * t598 + t555 * t600 + t558) * t230 + (t382 * t46 - t383 * t47 - t489 * t6 - t490 * t5) * mrSges(6,3) + t307 * t509 / 0.2e1 + t398 * t606 + (-m(5) * t26 + m(6) * t23 - t119 + t29) * t158 + (t146 * t295 - t147 * t296 - t227 * t236 + t228 * t237) * mrSges(4,3) + m(4) * (t146 * t245 + t147 * t244 + t179 * t228 + t180 * t227 - t269 * t340 - t322 * t344) + (-mrSges(4,2) * t340 + Ifges(4,1) * t296 + Ifges(4,4) * t295) * t207 + (t306 * t525 + t307 * t526 + t568) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t568) + (t273 * t536 + t395 * qJD(2) / 0.2e1 - t571) * qJD(2) + t623 * (t231 * t422 + t493 / 0.2e1) + t624 * (-t435 / 0.2e1 - t494 / 0.2e1) + (t171 * mrSges(5,2) - t26 * mrSges(5,3) + Ifges(5,1) * t139 + Ifges(5,4) * t138 + t565 * Ifges(5,5) + t12 * t401 + t23 * t402 + t553 * t621 + t554 * t620 + t555 * t619) * t231 + m(7) * (t1 * t45 + t107 * t12 + t3 * t48 + t30 * t77 + t31 * t7 + t35 * t8) - t272 * t459 / 0.2e1 + (mrSges(4,1) * t340 + Ifges(4,4) * t296 + Ifges(4,2) * t295) * t208 - t384 * t453 + (t626 / 0.2e1 + t598 * t550 + t600 * t548 + t597 * t546 - t607 + t636) * t172 + (-m(5) * t114 + m(6) * t110 + t578) * t50 + (Ifges(5,5) * t539 + t162 / 0.2e1 + t590 / 0.2e1 + t593 / 0.2e1 + t609) * t173 + (-t1 * t489 - t3 * t490 + t31 * t382 - t35 * t383) * mrSges(7,3) + m(6) * (t10 * t46 + t47 * t9 + t5 * t64 + t6 * t63) + t489 * t635 + t7 * t152 + t10 * t153 + t159 * t118 + t8 * t150 + t9 * t151 + t30 * t141 + Ifges(2,3) * qJDD(1) + t107 * t28 + t63 * t40 + t64 * t42 + t48 * t41 + t45 * t39 + t51 * t209 + t253 * t418 + t224 * t170 + (t595 / 0.2e1 + t596 / 0.2e1) * t353 + (t595 + t596) * t541 + t236 * t217 / 0.2e1 + t237 * t216 / 0.2e1 + t226 * t344 + t244 * t184 + t245 * t185 + t179 * t248 + t180 * t249 + (-t382 * t600 - t383 * t598) * t546 + t274 * (Ifges(4,4) * t236 + Ifges(4,2) * t237) / 0.2e1 + (-mrSges(3,1) * t526 - mrSges(3,2) * t525 + 0.2e1 * Ifges(3,6) * t536) * qJDD(2) + (Ifges(3,1) * t307 + Ifges(3,4) * t606 + Ifges(3,5) * qJDD(2) - t419 * t594) * t362 + t269 * (-mrSges(4,1) * t295 + mrSges(4,2) * t296) - pkin(1) * (-mrSges(3,1) * t306 + mrSges(3,2) * t307) - t322 * (-mrSges(4,1) * t237 + mrSges(4,2) * t236) + (Ifges(3,4) * t307 + Ifges(3,2) * t306) * t536 + t77 * (mrSges(7,1) * t383 - mrSges(7,2) * t382) + t110 * (mrSges(6,1) * t383 - mrSges(6,2) * t382); t184 * t534 + (-t265 * t360 - t54) * t153 + t625 * t151 + t272 * t461 / 0.2e1 + (t571 + (t384 - t381 / 0.2e1) * qJD(1)) * qJD(1) - (-Ifges(3,2) * t461 + t273 + t342) * t460 / 0.2e1 + t560 * t262 + t369 - m(4) * (t227 * t234 + t228 * t235 - t322 * t343) + (t575 * t363 + t569) * g(2) + (t367 * t575 + t570) * g(1) - t580 * t578 + t567 * (m(4) * t535 - m(5) * t308 + t406 + t618) - t226 * t343 + (m(4) * (t146 * t361 + t147 * t365 + (-t227 * t361 + t228 * t365) * qJD(3)) + t248 * t457 + t361 * t185 - t249 * t458) * pkin(2) + (-m(5) * t462 - m(6) * (t352 + t439) - m(7) * (t352 + t407) - m(4) * t352 + t320 + t561) * g(3) - t395 * t453 / 0.2e1 + t579 * t209 + (t114 * t580 + t115 * t579 - t240 * t247 + t26 * t266 + t267 * t27) * m(5) + (-t46 * t54 - t47 * t55 + t23 * t261 + t391 * t265 - g(1) * (-t367 * t529 + t287 + t317) - g(2) * (-t363 * t529 + t286 + t316) - t580 * t110) * m(6) + Ifges(3,3) * qJDD(2) + t582 * t141 + t587 * t152 + t588 * t150 + (-g(1) * t287 - g(2) * t286 + t1 * t242 + t12 * t250 + t243 * t3 + t31 * t587 + t35 * t588 + t582 * t77) * m(7) + t242 * t39 + t243 * t41 - t247 * t170 - t235 * t248 - t234 * t249 + t250 * t28 + t261 * t29 + t266 * t119 + t267 * t118 - t289 * mrSges(3,2) - t290 * mrSges(3,1) + Ifges(3,6) * t306 + Ifges(3,5) * t307 + t42 * t488; t42 * t479 + t119 * t530 + t118 * t531 + ((t26 * t358 + t27 * t357) * pkin(3) + t114 * t121 - t115 * t122 - t240 * t533) * m(5) + t560 * t335 + t369 - t170 * t533 + (t363 * t576 + t569) * g(2) + (t367 * t576 + t570) * g(1) + (m(5) * t532 + t618) * t567 + (-m(5) * t337 - m(6) * t439 - m(7) * t407 + t561) * g(3) - t578 * t121 + (-g(1) * t317 - g(2) * t316 - t110 * t121 + t23 * t336 - t46 * t52 - t47 * t53) * m(6) - t52 * t153 - t53 * t151 - t122 * t209 + t583 * t141 + t585 * t150 + t586 * t152 + (t1 * t284 + t12 * t315 + t285 * t3 + t31 * t586 + t35 * t585 + t583 * t77) * m(7) - t227 * t248 + t228 * t249 + t284 * t39 + t285 * t41 + t315 * t28 + t336 * t29; -t415 * t209 + (-t141 - t578) * t390 + (t39 + t40 + t215 * (t150 + t151)) * t364 + (t41 + t42 - t215 * (t152 + t153)) * t360 + t418 + (-g(1) * t363 + g(2) * t367) * t566 + (t1 * t364 - t390 * t77 + t3 * t360 + t215 * (-t31 * t360 + t35 * t364)) * m(7) + (-t110 * t390 + t215 * t391 + t360 * t5 + t364 * t6) * m(6) + (t114 * t390 - t115 * t415 + t171) * m(5); (t203 * t600 - t204 * t598) * t547 + (t203 * t602 - t637) * t549 + (t511 + t152 - m(7) * (-t31 + t34)) * t35 + t624 * t548 + t584 * t1 + t31 * t512 + (t360 * t584 + t402 + t517) * g(3) * t333 + (-t204 * t629 + t623 + t639) * t551 + (-t259 * t564 + t260 * t603) * g(1) + (t513 + t153) * t47 + (t514 - t151) * t46 + t558 - t34 * t150 - t110 * (mrSges(6,1) * t204 + mrSges(6,2) * t203) - t77 * (mrSges(7,1) * t204 + mrSges(7,2) * t203) + (t257 * t564 - t258 * t603) * g(2) + t572 + ((-m(7) * t77 - t141) * t204 + t39) * pkin(5); -t203 * t150 + t204 * t152 + (g(3) * t334 - t35 * t203 + t31 * t204 - t333 * t567 + t12) * m(7) + t28;];
tau  = t2;
