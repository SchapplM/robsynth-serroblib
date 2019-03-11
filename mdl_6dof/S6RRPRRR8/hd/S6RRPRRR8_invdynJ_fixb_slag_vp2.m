% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR8
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:22
% EndTime: 2019-03-09 14:03:33
% DurationCPUTime: 44.95s
% Computational Cost: add. (29209->992), mult. (66064->1320), div. (0->0), fcn. (50580->18), ass. (0->437)
t391 = sin(pkin(11));
t392 = cos(pkin(11));
t396 = sin(qJ(4));
t401 = cos(qJ(4));
t336 = t391 * t401 + t392 * t396;
t402 = cos(qJ(2));
t421 = t336 * t402;
t285 = qJD(1) * t421;
t312 = t336 * qJD(4);
t639 = t285 - t312;
t427 = t391 * t396 - t392 * t401;
t420 = t427 * t402;
t286 = qJD(1) * t420;
t311 = t427 * qJD(4);
t638 = t286 - t311;
t397 = sin(qJ(2));
t431 = pkin(2) * t397 - qJ(3) * t402;
t339 = t431 * qJD(1);
t482 = qJD(1) * t397;
t462 = t391 * t482;
t265 = pkin(7) * t462 + t392 * t339;
t492 = t392 * t402;
t426 = pkin(3) * t397 - pkin(8) * t492;
t225 = qJD(1) * t426 + t265;
t313 = t391 * t339;
t493 = t392 * t397;
t494 = t391 * t402;
t422 = -pkin(7) * t493 - pkin(8) * t494;
t248 = qJD(1) * t422 + t313;
t521 = pkin(8) + qJ(3);
t346 = t521 * t391;
t347 = t521 * t392;
t259 = -t396 * t346 + t401 * t347;
t614 = -t336 * qJD(3) - qJD(4) * t259 - t401 * t225 + t248 * t396;
t474 = qJD(4) * t401;
t477 = qJD(3) * t392;
t478 = qJD(3) * t391;
t613 = -t346 * t474 + (-t248 + t477) * t401 + (-qJD(4) * t347 - t225 - t478) * t396;
t648 = -pkin(4) * t482 - pkin(9) * t638 + t614;
t647 = -pkin(9) * t639 - t613;
t395 = sin(qJ(5));
t400 = cos(qJ(5));
t247 = t336 * t400 - t395 * t427;
t177 = -qJD(5) * t247 + t311 * t395 - t312 * t400;
t198 = -t285 * t400 + t286 * t395;
t646 = t177 - t198;
t258 = -t401 * t346 - t347 * t396;
t217 = -pkin(9) * t336 + t258;
t218 = -pkin(9) * t427 + t259;
t472 = qJD(5) * t400;
t473 = qJD(5) * t395;
t623 = t217 * t472 - t218 * t473 + t395 * t648 - t647 * t400;
t150 = t395 * t217 + t400 * t218;
t622 = -qJD(5) * t150 + t647 * t395 + t400 * t648;
t645 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t329 = qJD(2) * t392 - t462;
t460 = t392 * t482;
t330 = qJD(2) * t391 + t460;
t245 = t329 * t396 + t330 * t401;
t443 = t401 * t329 - t330 * t396;
t175 = t245 * t400 + t395 * t443;
t394 = sin(qJ(6));
t399 = cos(qJ(6));
t630 = -t245 * t395 + t400 * t443;
t644 = -t175 * t394 + t399 * t630;
t110 = t175 * t399 + t394 * t630;
t643 = pkin(10) * t646 + t623;
t246 = -t336 * t395 - t400 * t427;
t176 = qJD(5) * t246 - t311 * t400 - t312 * t395;
t199 = -t285 * t395 - t286 * t400;
t642 = -pkin(5) * t482 + t622 + (-t176 + t199) * pkin(10);
t376 = pkin(7) * t482;
t344 = -qJD(2) * pkin(2) + qJD(3) + t376;
t264 = -pkin(3) * t329 + t344;
t190 = -pkin(4) * t443 + t264;
t120 = -pkin(5) * t630 + t190;
t469 = qJD(1) * qJD(2);
t342 = qJDD(1) * t397 + t402 * t469;
t287 = qJDD(2) * t392 - t342 * t391;
t288 = qJDD(2) * t391 + t342 * t392;
t160 = qJD(4) * t443 + t287 * t396 + t288 * t401;
t161 = -qJD(4) * t245 + t287 * t401 - t288 * t396;
t74 = qJD(5) * t630 + t160 * t400 + t161 * t395;
t75 = -qJD(5) * t175 - t160 * t395 + t161 * t400;
t29 = qJD(6) * t644 + t394 * t75 + t399 * t74;
t30 = -qJD(6) * t110 - t394 * t74 + t399 * t75;
t381 = t402 * qJDD(1);
t454 = t397 * t469;
t341 = -t381 + t454;
t334 = qJDD(4) + t341;
t319 = qJDD(5) + t334;
t310 = qJDD(6) + t319;
t467 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t310;
t481 = qJD(1) * t402;
t365 = qJD(4) - t481;
t357 = qJD(5) + t365;
t448 = -qJ(3) * t397 - pkin(1);
t345 = -pkin(2) * t402 + t448;
t317 = t345 * qJD(1);
t377 = pkin(7) * t481;
t350 = qJD(2) * qJ(3) + t377;
t250 = t392 * t317 - t350 * t391;
t197 = -pkin(3) * t481 - pkin(8) * t330 + t250;
t251 = t391 * t317 + t392 * t350;
t200 = pkin(8) * t329 + t251;
t130 = t401 * t197 - t200 * t396;
t118 = -pkin(9) * t245 + t130;
t114 = pkin(4) * t365 + t118;
t131 = t197 * t396 + t200 * t401;
t119 = pkin(9) * t443 + t131;
t115 = t395 * t119;
t61 = t400 * t114 - t115;
t634 = pkin(10) * t175;
t47 = t61 - t634;
t45 = pkin(5) * t357 + t47;
t117 = t400 * t119;
t62 = t114 * t395 + t117;
t627 = pkin(10) * t630;
t48 = t62 + t627;
t509 = t394 * t48;
t16 = t399 * t45 - t509;
t476 = qJD(3) * t397;
t504 = qJDD(1) * pkin(1);
t233 = pkin(2) * t341 - qJ(3) * t342 - qJD(1) * t476 - t504;
t374 = pkin(7) * t381;
t292 = qJDD(2) * qJ(3) + t374 + (qJD(3) - t376) * qJD(2);
t186 = t392 * t233 - t292 * t391;
t140 = pkin(3) * t341 - pkin(8) * t288 + t186;
t187 = t391 * t233 + t392 * t292;
t154 = pkin(8) * t287 + t187;
t58 = -qJD(4) * t131 + t401 * t140 - t154 * t396;
t42 = pkin(4) * t334 - pkin(9) * t160 + t58;
t475 = qJD(4) * t396;
t57 = t396 * t140 + t401 * t154 + t197 * t474 - t200 * t475;
t44 = pkin(9) * t161 + t57;
t13 = -qJD(5) * t62 - t395 * t44 + t400 * t42;
t6 = pkin(5) * t319 - pkin(10) * t74 + t13;
t12 = t114 * t472 - t119 * t473 + t395 * t42 + t400 * t44;
t9 = pkin(10) * t75 + t12;
t2 = qJD(6) * t16 + t394 * t6 + t399 * t9;
t505 = t399 * t48;
t17 = t394 * t45 + t505;
t3 = -qJD(6) * t17 - t394 * t9 + t399 * t6;
t595 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t425 = t467 - t595;
t466 = Ifges(6,5) * t74 + Ifges(6,6) * t75 + Ifges(6,3) * t319;
t532 = mrSges(7,3) * t17;
t533 = mrSges(7,3) * t16;
t348 = qJD(6) + t357;
t542 = -t348 / 0.2e1;
t563 = -t110 / 0.2e1;
t565 = -t644 / 0.2e1;
t102 = Ifges(7,4) * t644;
t55 = Ifges(7,1) * t110 + Ifges(7,5) * t348 + t102;
t575 = -t55 / 0.2e1;
t510 = Ifges(7,4) * t110;
t54 = Ifges(7,2) * t644 + Ifges(7,6) * t348 + t510;
t577 = -t54 / 0.2e1;
t593 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t637 = t425 + t466 - t593 + (-mrSges(7,2) * t120 + Ifges(7,1) * t563 + Ifges(7,4) * t565 + Ifges(7,5) * t542 + t533 + t575) * t644 - (mrSges(7,1) * t120 + Ifges(7,4) * t563 + Ifges(7,2) * t565 + Ifges(7,6) * t542 - t532 + t577) * t110;
t390 = pkin(11) + qJ(4);
t380 = cos(t390);
t367 = pkin(4) * t380;
t370 = pkin(3) * t392 + pkin(2);
t338 = t367 + t370;
t382 = qJ(5) + t390;
t369 = cos(t382);
t359 = pkin(5) * t369;
t293 = t338 + t359;
t371 = qJ(6) + t382;
t360 = sin(t371);
t361 = cos(t371);
t368 = sin(t382);
t379 = sin(t390);
t636 = -m(5) * t370 - m(6) * t338 - m(7) * t293 - mrSges(5,1) * t380 - mrSges(6,1) * t369 - mrSges(7,1) * t361 + mrSges(5,2) * t379 + mrSges(6,2) * t368 + mrSges(7,2) * t360;
t389 = -pkin(9) - t521;
t383 = -pkin(10) + t389;
t635 = -m(5) * t521 + m(6) * t389 + m(7) * t383 - t645;
t512 = Ifges(5,4) * t245;
t163 = Ifges(5,2) * t443 + t365 * Ifges(5,6) + t512;
t559 = t163 / 0.2e1;
t239 = Ifges(5,4) * t443;
t164 = t245 * Ifges(5,1) + t365 * Ifges(5,5) + t239;
t558 = t164 / 0.2e1;
t528 = pkin(5) * t175;
t626 = qJD(2) / 0.2e1;
t398 = sin(qJ(1));
t403 = cos(qJ(1));
t631 = g(1) * t403 + g(2) * t398;
t461 = t391 * t481;
t320 = pkin(3) * t461 + t377;
t603 = -pkin(4) * t639 - t320;
t581 = t29 / 0.2e1;
t580 = t30 / 0.2e1;
t573 = t74 / 0.2e1;
t572 = t75 / 0.2e1;
t561 = t160 / 0.2e1;
t560 = t161 / 0.2e1;
t548 = t287 / 0.2e1;
t547 = t288 / 0.2e1;
t546 = t310 / 0.2e1;
t545 = t319 / 0.2e1;
t544 = t334 / 0.2e1;
t629 = -t341 / 0.2e1;
t543 = t341 / 0.2e1;
t628 = t342 / 0.2e1;
t149 = t400 * t217 - t218 * t395;
t124 = -pkin(10) * t247 + t149;
t125 = pkin(10) * t246 + t150;
t70 = t124 * t394 + t125 * t399;
t625 = -qJD(6) * t70 - t394 * t643 + t642 * t399;
t69 = t124 * t399 - t125 * t394;
t624 = qJD(6) * t69 + t642 * t394 + t399 * t643;
t529 = pkin(4) * t400;
t372 = pkin(5) + t529;
t470 = qJD(6) * t399;
t471 = qJD(6) * t394;
t490 = t395 * t399;
t63 = -t118 * t395 - t117;
t51 = t63 - t627;
t64 = t400 * t118 - t115;
t52 = t64 - t634;
t621 = t394 * t52 - t399 * t51 - t372 * t471 + (-t395 * t470 + (-t394 * t400 - t490) * qJD(5)) * pkin(4);
t491 = t394 * t395;
t620 = -t394 * t51 - t399 * t52 + t372 * t470 + (-t395 * t471 + (t399 * t400 - t491) * qJD(5)) * pkin(4);
t583 = m(6) * pkin(4);
t619 = t583 + mrSges(5,1);
t440 = -mrSges(4,1) * t392 + mrSges(4,2) * t391;
t419 = m(4) * pkin(2) - t440;
t616 = t402 * t419;
t326 = t392 * t345;
t249 = -pkin(8) * t493 + t326 + (-pkin(7) * t391 - pkin(3)) * t402;
t282 = pkin(7) * t492 + t391 * t345;
t495 = t391 * t397;
t257 = -pkin(8) * t495 + t282;
t181 = t401 * t249 - t257 * t396;
t302 = t427 * t397;
t147 = -pkin(4) * t402 + pkin(9) * t302 + t181;
t182 = t396 * t249 + t401 * t257;
t301 = t336 * t397;
t152 = -pkin(9) * t301 + t182;
t93 = t395 * t147 + t400 * t152;
t615 = -pkin(5) * t646 + t603;
t610 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t329 - mrSges(4,2) * t330 - mrSges(3,3) * t482;
t513 = Ifges(4,4) * t392;
t434 = -Ifges(4,2) * t391 + t513;
t514 = Ifges(4,4) * t391;
t436 = Ifges(4,1) * t392 - t514;
t609 = t329 * (Ifges(4,6) * t397 + t402 * t434) + t330 * (Ifges(4,5) * t397 + t402 * t436);
t437 = -mrSges(7,1) * t360 - mrSges(7,2) * t361;
t608 = mrSges(6,1) * t368 + mrSges(6,2) * t369 - t437;
t488 = t402 * t403;
t279 = -t368 * t488 + t398 * t369;
t280 = t398 * t368 + t369 * t488;
t269 = -t360 * t488 + t398 * t361;
t270 = t398 * t360 + t361 * t488;
t486 = t269 * mrSges(7,1) - t270 * mrSges(7,2);
t607 = -t279 * mrSges(6,1) + t280 * mrSges(6,2) - t486;
t489 = t398 * t402;
t277 = t368 * t489 + t369 * t403;
t278 = t368 * t403 - t369 * t489;
t267 = t360 * t489 + t361 * t403;
t268 = t360 * t403 - t361 * t489;
t487 = -t267 * mrSges(7,1) + t268 * mrSges(7,2);
t606 = t277 * mrSges(6,1) - t278 * mrSges(6,2) - t487;
t327 = -pkin(7) * t454 + t374;
t328 = t342 * pkin(7);
t605 = t327 * t402 + t328 * t397;
t604 = -t186 * t391 + t187 * t392;
t602 = m(4) + m(5) + m(6) + m(7);
t439 = mrSges(4,1) * t391 + t392 * mrSges(4,2);
t601 = -t344 * t402 * t439 - t251 * (-mrSges(4,2) * t397 - mrSges(4,3) * t494) - t250 * (mrSges(4,1) * t397 - mrSges(4,3) * t492);
t600 = -mrSges(6,1) * t190 + mrSges(6,3) * t62;
t599 = mrSges(6,2) * t190 - mrSges(6,3) * t61;
t598 = -m(3) - t602;
t375 = Ifges(3,4) * t481;
t597 = t392 * (t330 * Ifges(4,1) + t329 * Ifges(4,4) - Ifges(4,5) * t481) + Ifges(3,1) * t482 + Ifges(3,5) * qJD(2) + t375;
t596 = t330 * Ifges(4,5) + t245 * Ifges(5,5) + t175 * Ifges(6,5) + t110 * Ifges(7,5) + t329 * Ifges(4,6) + Ifges(5,6) * t443 + Ifges(6,6) * t630 + Ifges(7,6) * t644 - Ifges(4,3) * t481 + t365 * Ifges(5,3) + t357 * Ifges(6,3) + t348 * Ifges(7,3);
t594 = -t58 * mrSges(5,1) + t57 * mrSges(5,2);
t442 = mrSges(3,1) * t402 - mrSges(3,2) * t397;
t592 = t397 * t645 + mrSges(2,1) + t442;
t366 = pkin(4) * t379;
t527 = pkin(5) * t368;
t322 = -t366 - t527;
t384 = t391 * pkin(3);
t590 = -m(5) * t384 - m(6) * (t366 + t384) - m(7) * (-t322 + t384) + mrSges(2,2) - mrSges(3,3) - t439;
t516 = Ifges(3,4) * t397;
t435 = Ifges(3,2) * t402 + t516;
t587 = t131 * mrSges(5,2) + t17 * mrSges(7,2) + t62 * mrSges(6,2) + Ifges(3,6) * t626 + qJD(1) * t435 / 0.2e1 - t130 * mrSges(5,1) - t16 * mrSges(7,1) - t61 * mrSges(6,1);
t585 = Ifges(7,4) * t581 + Ifges(7,2) * t580 + Ifges(7,6) * t546;
t584 = Ifges(7,1) * t581 + Ifges(7,4) * t580 + Ifges(7,5) * t546;
t582 = m(7) * pkin(5);
t579 = Ifges(6,4) * t573 + Ifges(6,2) * t572 + Ifges(6,6) * t545;
t578 = Ifges(6,1) * t573 + Ifges(6,4) * t572 + Ifges(6,5) * t545;
t576 = t54 / 0.2e1;
t574 = t55 / 0.2e1;
t571 = Ifges(5,4) * t561 + Ifges(5,2) * t560 + Ifges(5,6) * t544;
t570 = Ifges(5,1) * t561 + Ifges(5,4) * t560 + Ifges(5,5) * t544;
t511 = Ifges(6,4) * t175;
t98 = Ifges(6,2) * t630 + Ifges(6,6) * t357 + t511;
t569 = -t98 / 0.2e1;
t568 = t98 / 0.2e1;
t169 = Ifges(6,4) * t630;
t99 = Ifges(6,1) * t175 + Ifges(6,5) * t357 + t169;
t567 = -t99 / 0.2e1;
t566 = t99 / 0.2e1;
t564 = t644 / 0.2e1;
t562 = t110 / 0.2e1;
t557 = -t630 / 0.2e1;
t556 = t630 / 0.2e1;
t555 = -t175 / 0.2e1;
t554 = t175 / 0.2e1;
t553 = Ifges(4,1) * t547 + Ifges(4,4) * t548 + Ifges(4,5) * t543;
t552 = -t443 / 0.2e1;
t551 = t443 / 0.2e1;
t550 = -t245 / 0.2e1;
t549 = t245 / 0.2e1;
t541 = t348 / 0.2e1;
t540 = -t357 / 0.2e1;
t539 = t357 / 0.2e1;
t538 = -t365 / 0.2e1;
t537 = t365 / 0.2e1;
t531 = pkin(4) * t245;
t524 = g(3) * t397;
t385 = t397 * pkin(7);
t520 = mrSges(5,3) * t443;
t519 = mrSges(5,3) * t245;
t518 = mrSges(6,3) * t630;
t517 = mrSges(6,3) * t175;
t515 = Ifges(3,4) * t402;
t303 = -qJDD(2) * pkin(2) + qJDD(3) + t328;
t499 = t303 * t397;
t305 = qJD(2) * t431 - t476;
t480 = qJD(2) * t397;
t465 = pkin(7) * t480;
t255 = t392 * t305 + t391 * t465;
t479 = qJD(2) * t402;
t378 = pkin(7) * t479;
t459 = t391 * t479;
t321 = pkin(3) * t459 + t378;
t340 = pkin(3) * t495 + t385;
t463 = Ifges(5,5) * t160 + Ifges(5,6) * t161 + Ifges(5,3) * t334;
t458 = m(4) * qJ(3) + mrSges(4,3);
t37 = -t75 * mrSges(6,1) + t74 * mrSges(6,2);
t10 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t447 = -t469 / 0.2e1;
t204 = -t287 * mrSges(4,1) + t288 * mrSges(4,2);
t94 = -t161 * mrSges(5,1) + t160 * mrSges(5,2);
t92 = t400 * t147 - t152 * t395;
t220 = -qJD(2) * t421 + t311 * t397;
t193 = -pkin(4) * t220 + t321;
t253 = pkin(4) * t301 + t340;
t441 = mrSges(3,1) * t397 + mrSges(3,2) * t402;
t433 = Ifges(3,5) * t402 - Ifges(3,6) * t397;
t432 = Ifges(4,5) * t392 - Ifges(4,6) * t391;
t213 = -t301 * t395 - t302 * t400;
t76 = -pkin(5) * t402 - pkin(10) * t213 + t92;
t212 = -t301 * t400 + t302 * t395;
t77 = pkin(10) * t212 + t93;
t38 = -t394 * t77 + t399 * t76;
t39 = t394 * t76 + t399 * t77;
t144 = t212 * t399 - t213 * t394;
t145 = t212 * t394 + t213 * t399;
t178 = t246 * t399 - t247 * t394;
t179 = t246 * t394 + t247 * t399;
t430 = t293 * t402 - t383 * t397;
t429 = t338 * t402 - t389 * t397;
t428 = t370 * t402 + t397 * t521;
t289 = pkin(4) * t427 - t370;
t424 = pkin(1) * t441;
t296 = -t379 * t488 + t398 * t380;
t294 = t379 * t489 + t380 * t403;
t423 = t397 * (Ifges(3,1) * t402 - t516);
t214 = qJD(2) * t426 + t255;
t290 = t391 * t305;
t226 = qJD(2) * t422 + t290;
t101 = -qJD(4) * t182 + t401 * t214 - t226 * t396;
t219 = -qJD(2) * t420 - t312 * t397;
t86 = pkin(4) * t480 - pkin(9) * t219 + t101;
t100 = t396 * t214 + t401 * t226 + t249 * t474 - t257 * t475;
t88 = pkin(9) * t220 + t100;
t31 = t147 * t472 - t152 * t473 + t395 * t86 + t400 * t88;
t222 = -pkin(3) * t287 + t303;
t123 = -pkin(4) * t161 + t222;
t32 = -qJD(5) * t93 - t395 * t88 + t400 * t86;
t410 = t402 * (Ifges(4,3) * t397 + t402 * t432);
t408 = t397 * t458 + t616;
t351 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t481;
t323 = t359 + t367;
t308 = pkin(4) * t490 + t372 * t394;
t307 = -pkin(4) * t491 + t372 * t399;
t297 = t398 * t379 + t380 * t488;
t295 = t379 * t403 - t380 * t489;
t284 = -mrSges(4,1) * t481 - mrSges(4,3) * t330;
t283 = mrSges(4,2) * t481 + mrSges(4,3) * t329;
t281 = -pkin(7) * t494 + t326;
t266 = -pkin(7) * t460 + t313;
t256 = -t392 * t465 + t290;
t235 = t330 * Ifges(4,4) + t329 * Ifges(4,2) - Ifges(4,6) * t481;
t231 = mrSges(4,1) * t341 - mrSges(4,3) * t288;
t230 = -mrSges(4,2) * t341 + mrSges(4,3) * t287;
t211 = mrSges(5,1) * t365 - t519;
t210 = -mrSges(5,2) * t365 + t520;
t196 = -pkin(5) * t246 + t289;
t191 = t288 * Ifges(4,4) + t287 * Ifges(4,2) + t341 * Ifges(4,6);
t180 = -mrSges(5,1) * t443 + mrSges(5,2) * t245;
t168 = -pkin(5) * t212 + t253;
t156 = mrSges(6,1) * t357 - t517;
t155 = -mrSges(6,2) * t357 + t518;
t138 = -mrSges(5,2) * t334 + mrSges(5,3) * t161;
t137 = mrSges(5,1) * t334 - mrSges(5,3) * t160;
t133 = t198 * t394 + t199 * t399;
t132 = t198 * t399 - t199 * t394;
t126 = t531 + t528;
t122 = -qJD(5) * t213 - t219 * t395 + t220 * t400;
t121 = qJD(5) * t212 + t219 * t400 + t220 * t395;
t113 = -mrSges(6,1) * t630 + mrSges(6,2) * t175;
t96 = mrSges(7,1) * t348 - mrSges(7,3) * t110;
t95 = -mrSges(7,2) * t348 + mrSges(7,3) * t644;
t89 = -pkin(5) * t122 + t193;
t79 = -qJD(6) * t179 - t176 * t394 + t177 * t399;
t78 = qJD(6) * t178 + t176 * t399 + t177 * t394;
t66 = -mrSges(6,2) * t319 + mrSges(6,3) * t75;
t65 = mrSges(6,1) * t319 - mrSges(6,3) * t74;
t56 = -mrSges(7,1) * t644 + mrSges(7,2) * t110;
t50 = -qJD(6) * t145 - t121 * t394 + t122 * t399;
t49 = qJD(6) * t144 + t121 * t399 + t122 * t394;
t46 = -pkin(5) * t75 + t123;
t25 = -mrSges(7,2) * t310 + mrSges(7,3) * t30;
t24 = mrSges(7,1) * t310 - mrSges(7,3) * t29;
t21 = pkin(10) * t122 + t31;
t20 = pkin(5) * t480 - pkin(10) * t121 + t32;
t19 = t399 * t47 - t509;
t18 = -t394 * t47 - t505;
t5 = -qJD(6) * t39 + t20 * t399 - t21 * t394;
t4 = qJD(6) * t38 + t20 * t394 + t21 * t399;
t1 = [t597 * t479 / 0.2e1 + (-t186 * t493 - t187 * t495) * mrSges(4,3) + (Ifges(7,5) * t49 + Ifges(7,6) * t50) * t541 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t546 + (Ifges(6,5) * t121 + Ifges(6,6) * t122) * t539 + (Ifges(6,5) * t213 + Ifges(6,6) * t212) * t545 + (Ifges(5,5) * t219 + Ifges(5,6) * t220) * t537 + (Ifges(7,4) * t49 + Ifges(7,2) * t50) * t564 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t580 + (t596 / 0.2e1 + Ifges(5,5) * t549 + Ifges(6,5) * t554 + Ifges(7,5) * t562 + Ifges(5,6) * t551 + Ifges(6,6) * t556 + Ifges(7,6) * t564 + Ifges(5,3) * t537 + Ifges(6,3) * t539 + Ifges(7,3) * t541 - t587) * t480 + (-t297 * mrSges(5,1) - t280 * mrSges(6,1) - t270 * mrSges(7,1) - t296 * mrSges(5,2) - t279 * mrSges(6,2) - t269 * mrSges(7,2) + t598 * (t403 * pkin(1) + t398 * pkin(7)) + t590 * t398 + (-m(5) * t428 - m(6) * t429 - m(7) * t430 - t408 - t592) * t403) * g(2) + (-pkin(7) * t341 * t402 + t342 * t385 + t605) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t605) + (-mrSges(3,1) * t385 + Ifges(3,5) * t397 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * t402) * qJDD(2) + (Ifges(5,1) * t219 + Ifges(5,4) * t220) * t549 - t610 * t378 + (Ifges(5,4) * t219 + Ifges(5,2) * t220) * t551 + (Ifges(6,4) * t121 + Ifges(6,2) * t122) * t556 + (Ifges(6,4) * t213 + Ifges(6,2) * t212) * t572 + t120 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t442 * t504 + (-Ifges(5,5) * t302 - Ifges(5,6) * t301) * t544 + (-Ifges(5,1) * t302 - Ifges(5,4) * t301) * t561 + (-Ifges(5,4) * t302 - Ifges(5,2) * t301) * t560 + (-t130 * t219 + t131 * t220 - t301 * t57 + t302 * t58) * mrSges(5,3) + t222 * (mrSges(5,1) * t301 - mrSges(5,2) * t302) + (Ifges(6,1) * t121 + Ifges(6,4) * t122) * t554 + (Ifges(6,1) * t213 + Ifges(6,4) * t212) * t573 + t168 * t10 - (Ifges(4,5) * t288 + Ifges(4,6) * t287 + Ifges(4,3) * t341 + t463 + t466 + t467) * t402 / 0.2e1 + t4 * t95 + t5 * t96 + t89 * t56 + t92 * t65 + t93 * t66 + Ifges(2,3) * qJDD(1) + t31 * t155 + t32 * t156 + t46 * (-mrSges(7,1) * t144 + mrSges(7,2) * t145) + t39 * t25 + t38 * t24 + (t144 * t2 - t145 * t3 - t16 * t49 + t17 * t50) * mrSges(7,3) + t515 * t628 + t435 * t629 + t609 * t626 + (t12 * t212 - t121 * t61 + t122 * t62 - t13 * t213) * mrSges(6,3) + m(4) * (t186 * t281 + t187 * t282 + t250 * t255 + t251 * t256 + (t344 * t479 + t499) * pkin(7)) - t191 * t495 / 0.2e1 + (-t295 * mrSges(5,1) - t278 * mrSges(6,1) - t268 * mrSges(7,1) - t294 * mrSges(5,2) - t277 * mrSges(6,2) - t267 * mrSges(7,2) + (-m(4) * t448 + t397 * mrSges(4,3) + t616 + m(3) * pkin(1) - m(5) * (-pkin(1) - t428) - m(6) * (-pkin(1) - t429) - m(7) * (-pkin(1) - t430) + t592) * t398 + (pkin(7) * t598 + t590) * t403) * g(1) - pkin(1) * (mrSges(3,1) * t341 + mrSges(3,2) * t342) + t340 * t94 + (Ifges(7,1) * t49 + Ifges(7,4) * t50) * t562 + (Ifges(7,1) * t145 + Ifges(7,4) * t144) * t581 + t204 * t385 + (t402 * (-Ifges(3,2) * t397 + t515) + t423) * t469 / 0.2e1 + t181 * t137 + t182 * t138 + m(5) * (t100 * t131 + t101 * t130 + t181 * t58 + t182 * t57 + t222 * t340 + t264 * t321) + m(7) * (t120 * t89 + t16 * t5 + t168 * t46 + t17 * t4 + t2 * t39 + t3 * t38) + m(6) * (t12 * t93 + t123 * t253 + t13 * t92 + t190 * t193 + t31 * t62 + t32 * t61) - t235 * t459 / 0.2e1 + t190 * (-mrSges(6,1) * t122 + mrSges(6,2) * t121) + t193 * t113 + (t433 * t626 - t601) * qJD(2) + t100 * t210 + t101 * t211 + t123 * (-mrSges(6,1) * t212 + mrSges(6,2) * t213) + (t593 - Ifges(4,3) * t543 - Ifges(5,3) * t544 - Ifges(6,3) * t545 - Ifges(7,3) * t546 - Ifges(4,5) * t547 - Ifges(4,6) * t548 - Ifges(5,6) * t560 - Ifges(5,5) * t561 - Ifges(6,6) * t572 - Ifges(6,5) * t573 - Ifges(7,6) * t580 - Ifges(7,5) * t581 - t186 * mrSges(4,1) + t187 * mrSges(4,2) + Ifges(3,2) * t629 + Ifges(3,4) * t628 + t594 + t595) * t402 + (Ifges(3,1) * t342 + Ifges(3,4) * t629 + t432 * t543 + t434 * t548 + t436 * t547) * t397 + t410 * t447 + t493 * t553 + t219 * t558 + t220 * t559 + t121 * t566 + t122 * t568 - t302 * t570 - t301 * t571 + t49 * t574 + t50 * t576 + t213 * t578 + t212 * t579 + t145 * t584 + t144 * t585 + t253 * t37 + t264 * (-mrSges(5,1) * t220 + mrSges(5,2) * t219) - t351 * t465 + t281 * t231 + t282 * t230 + t256 * t283 + t255 * t284 + t439 * t499 - t424 * t469 + t321 * t180; -t596 * t482 / 0.2e1 - (-Ifges(3,2) * t482 + t375 + t597) * t481 / 0.2e1 + (-mrSges(5,1) * t639 + mrSges(5,2) * t638) * t264 + (-t130 * t638 + t131 * t639 - t336 * t58 - t427 * t57) * mrSges(5,3) + (-Ifges(5,4) * t286 - Ifges(5,2) * t285) * t552 + (Ifges(6,5) * t199 + Ifges(6,6) * t198) * t540 + (Ifges(7,5) * t133 + Ifges(7,6) * t132) * t542 + t631 * t441 + (-t442 - t408) * g(3) + (-t265 - t478) * t284 + (Ifges(7,4) * t562 + Ifges(7,2) * t564 + Ifges(7,6) * t541 + t532 + t576) * t79 + (Ifges(6,4) * t199 + Ifges(6,2) * t198) * t557 + (Ifges(6,1) * t554 + Ifges(6,4) * t556 + Ifges(6,5) * t539 + t566 + t599) * t176 + (Ifges(6,4) * t554 + Ifges(6,2) * t556 + Ifges(6,6) * t539 + t568 + t600) * t177 + t613 * t210 + (t130 * t614 + t131 * t613 - t222 * t370 + t258 * t58 + t259 * t57 - t264 * t320) * m(5) + t614 * t211 + t615 * t56 + t603 * t113 + (-pkin(2) * t303 + (-t250 * t391 + t251 * t392) * qJD(3) + t604 * qJ(3) - t250 * t265 - t251 * t266 - t344 * t377) * m(4) + t604 * mrSges(4,3) + t351 * t376 + (Ifges(7,1) * t562 + Ifges(7,4) * t564 + Ifges(7,5) * t541 - t533 + t574) * t78 + (-Ifges(5,5) * t311 - Ifges(5,6) * t312) * t537 + (-Ifges(5,1) * t311 - Ifges(5,4) * t312) * t549 + (-Ifges(5,4) * t311 - Ifges(5,2) * t312) * t551 + (-t132 * t17 + t133 * t16 + t178 * t2 - t179 * t3) * mrSges(7,3) + (t230 * t392 - t231 * t391) * qJ(3) + t610 * t377 + (Ifges(7,4) * t133 + Ifges(7,2) * t132) * t565 + (Ifges(6,1) * t199 + Ifges(6,4) * t198) * t555 + (Ifges(5,5) * t336 - Ifges(5,6) * t427) * t544 + (Ifges(5,4) * t336 - Ifges(5,2) * t427) * t560 + (Ifges(5,1) * t336 - Ifges(5,4) * t427) * t561 + t222 * (mrSges(5,1) * t427 + mrSges(5,2) * t336) - t427 * t571 + (-Ifges(5,5) * t286 - Ifges(5,6) * t285) * t538 + (-Ifges(5,1) * t286 - Ifges(5,4) * t285) * t550 + (Ifges(7,1) * t133 + Ifges(7,4) * t132) * t563 + (Ifges(5,5) * t550 + Ifges(6,5) * t555 + Ifges(7,5) * t563 + Ifges(5,6) * t552 + Ifges(6,6) * t557 + Ifges(7,6) * t565 + Ifges(5,3) * t538 + Ifges(6,3) * t540 + Ifges(7,3) * t542 + t587) * t482 + t303 * t440 + (t601 - t609 / 0.2e1 + (t424 + t410 / 0.2e1 - t423 / 0.2e1) * qJD(1)) * qJD(1) + ((-t133 + t78) * mrSges(7,2) + (t132 - t79) * mrSges(7,1)) * t120 + (t12 * t246 - t13 * t247 - t198 * t62 + t199 * t61) * mrSges(6,3) + t69 * t24 + t70 * t25 + t149 * t65 + t150 * t66 + Ifges(3,3) * qJDD(2) + t392 * t191 / 0.2e1 - t370 * t94 - Ifges(3,6) * t341 + Ifges(3,5) * t342 + (t631 * (t419 - t636) + t635 * g(3)) * t397 + (t631 * (-t458 + t635) + t636 * g(3)) * t402 + t46 * (-mrSges(7,1) * t178 + mrSges(7,2) * t179) + t638 * t558 + t639 * t559 + t196 * t10 + t622 * t156 + t623 * t155 + t235 * t461 / 0.2e1 + (t12 * t150 + t123 * t289 + t13 * t149 + t190 * t603 + t61 * t622 + t62 * t623) * m(6) + t624 * t95 + t625 * t96 + (t120 * t615 + t16 * t625 + t17 * t624 + t196 * t46 + t2 * t70 + t3 * t69) * m(7) - t190 * (-mrSges(6,1) * t198 + mrSges(6,2) * t199) - pkin(2) * t204 + t123 * (-mrSges(6,1) * t246 + mrSges(6,2) * t247) + t433 * t447 + (Ifges(4,5) * t391 + Ifges(4,6) * t392) * t543 + (Ifges(6,5) * t247 + Ifges(6,6) * t246) * t545 + (Ifges(7,5) * t179 + Ifges(7,6) * t178) * t546 + (Ifges(4,1) * t391 + t513) * t547 + (Ifges(4,2) * t392 + t514) * t548 + t391 * t553 + t199 * t567 + t198 * t569 + t336 * t570 + (Ifges(6,4) * t247 + Ifges(6,2) * t246) * t572 + (Ifges(6,1) * t247 + Ifges(6,4) * t246) * t573 + t133 * t575 + t132 * t577 + t247 * t578 + t246 * t579 + (Ifges(7,4) * t179 + Ifges(7,2) * t178) * t580 + (Ifges(7,1) * t179 + Ifges(7,4) * t178) * t581 + t179 * t584 + t178 * t585 + t258 * t137 + t259 * t138 + t289 * t37 + (t477 - t266) * t283 - t320 * t180 - t327 * mrSges(3,2) - t328 * mrSges(3,1); t110 * t96 - t644 * t95 - t630 * t155 + t175 * t156 - t443 * t210 + t245 * t211 - t329 * t283 + t330 * t284 + t10 + t204 + t37 + t94 + (t110 * t16 - t17 * t644 + t46) * m(7) + (t175 * t61 - t62 * t630 + t123) * m(6) + (t130 * t245 - t131 * t443 + t222) * m(5) + (t250 * t330 - t251 * t329 + t303) * m(4) + (t402 * g(3) - t397 * t631) * t602; t637 + (Ifges(6,1) * t555 + Ifges(6,4) * t557 + Ifges(6,5) * t540 + t567 - t599) * t630 - t594 + (-m(7) * (t322 * t489 - t323 * t403) - mrSges(5,2) * t295 + t619 * t294 + t606) * g(2) + (-m(7) * (t322 * t488 + t398 * t323) + mrSges(5,2) * t297 - t619 * t296 + t607) * g(1) + t620 * t95 + (-t120 * t126 + t16 * t621 + t17 * t620 + t2 * t308 + t3 * t307) * m(7) + t621 * t96 + t463 + (m(6) * t366 - m(7) * t322 + mrSges(5,1) * t379 + mrSges(5,2) * t380 + t608) * t524 - t126 * t56 + (t520 - t210) * t130 - (Ifges(6,4) * t555 + Ifges(6,2) * t557 + Ifges(6,6) * t540 + t569 - t600) * t175 + (Ifges(5,5) * t443 - Ifges(5,6) * t245) * t538 + (Ifges(5,1) * t443 - t512) * t550 - t264 * (mrSges(5,1) * t245 + mrSges(5,2) * t443) + (t155 * t472 - t156 * t473 + t395 * t66) * pkin(4) - t64 * t155 - t63 * t156 + (t519 + t211) * t131 - t113 * t531 - m(6) * (t190 * t531 + t61 * t63 + t62 * t64) + (-Ifges(5,2) * t245 + t164 + t239) * t552 + t65 * t529 + t163 * t549 + (t12 * t395 + t13 * t400 + (-t395 * t61 + t400 * t62) * qJD(5)) * t583 + t307 * t24 + t308 * t25; -t19 * t95 - t18 * t96 - t56 * t528 - m(7) * (t120 * t528 + t16 * t18 + t17 * t19) - t190 * (mrSges(6,1) * t175 + mrSges(6,2) * t630) + (Ifges(6,5) * t630 - Ifges(6,6) * t175) * t540 + t98 * t554 + (Ifges(6,1) * t630 - t511) * t555 + (t2 * t394 + t3 * t399 + (-t16 * t394 + t17 * t399) * qJD(6)) * t582 + (t517 + t156) * t62 + (t518 - t155) * t61 + (-Ifges(6,2) * t175 + t169 + t99) * t557 + (m(7) * t527 + t608) * t524 + (t277 * t582 + t606) * g(2) + (-t279 * t582 + t607) * g(1) + (t24 * t399 + t25 * t394 + t470 * t95 - t471 * t96) * pkin(5) + t637; -t120 * (mrSges(7,1) * t110 + mrSges(7,2) * t644) + (Ifges(7,1) * t644 - t510) * t563 + t54 * t562 + (Ifges(7,5) * t644 - Ifges(7,6) * t110) * t542 - t16 * t95 + t17 * t96 - g(1) * t486 - g(2) * t487 - t437 * t524 + (t110 * t17 + t16 * t644) * mrSges(7,3) + t425 + (-Ifges(7,2) * t110 + t102 + t55) * t565;];
tau  = t1;
