% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:36
% EndTime: 2019-03-09 12:56:49
% DurationCPUTime: 46.31s
% Computational Cost: add. (13100->932), mult. (31162->1203), div. (0->0), fcn. (23376->10), ass. (0->417)
t629 = Ifges(6,4) + Ifges(7,4);
t325 = cos(qJ(2));
t317 = cos(pkin(6));
t469 = qJD(1) * t317;
t451 = pkin(1) * t469;
t295 = t325 * t451;
t653 = qJD(3) - t295;
t316 = sin(pkin(6));
t321 = sin(qJ(2));
t490 = t316 * t321;
t431 = qJD(2) * t490;
t488 = t316 * t325;
t242 = qJD(1) * t431 - qJDD(1) * t488;
t454 = qJDD(1) * t317;
t299 = qJDD(2) + t454;
t320 = sin(qJ(4));
t324 = cos(qJ(4));
t301 = qJD(2) + t469;
t470 = qJD(1) * t316;
t432 = t325 * t470;
t336 = -t301 * t324 + t320 * t432;
t126 = qJD(4) * t336 + t242 * t324 - t299 * t320;
t123 = qJDD(5) - t126;
t560 = t123 / 0.2e1;
t216 = -t301 * t320 - t324 * t432;
t125 = qJD(4) * t216 + t242 * t320 + t299 * t324;
t433 = t321 * t470;
t273 = qJD(4) + t433;
t319 = sin(qJ(5));
t323 = cos(qJ(5));
t154 = t273 * t319 - t323 * t336;
t466 = qJD(2) * t325;
t243 = (qJD(1) * t466 + qJDD(1) * t321) * t316;
t228 = qJDD(4) + t243;
t62 = -qJD(5) * t154 - t125 * t319 + t228 * t323;
t566 = t62 / 0.2e1;
t153 = t273 * t323 + t319 * t336;
t61 = qJD(5) * t153 + t125 * t323 + t228 * t319;
t567 = t61 / 0.2e1;
t626 = Ifges(6,5) + Ifges(7,5);
t630 = Ifges(6,1) + Ifges(7,1);
t640 = t560 * t626 + t566 * t629 + t567 * t630;
t625 = Ifges(6,2) + Ifges(7,2);
t624 = Ifges(6,6) + Ifges(7,6);
t623 = Ifges(6,3) + Ifges(7,3);
t382 = pkin(4) * t324 + pkin(10) * t320;
t562 = pkin(3) + pkin(8);
t652 = -(-t382 - t562) * t433 + qJD(4) * t382 + t653;
t563 = pkin(2) + pkin(9);
t146 = -t301 * t563 + t433 * t562 + t653;
t407 = -qJ(3) * t321 - pkin(1);
t169 = (-t325 * t563 + t407) * t470;
t461 = qJD(4) * t324;
t463 = qJD(4) * t320;
t293 = pkin(8) * t432;
t302 = pkin(8) * t490;
t537 = pkin(1) * t317;
t450 = qJD(2) * t537;
t396 = qJD(1) * t450;
t444 = pkin(1) * t454;
t156 = -qJD(2) * t293 - qJDD(1) * t302 - t321 * t396 + t325 * t444;
t333 = qJDD(3) - t156;
t92 = pkin(3) * t243 - t299 * t563 + t333;
t465 = qJD(3) * t321;
t328 = -qJ(3) * t243 + (-pkin(1) * qJDD(1) - qJD(1) * t465) * t316;
t98 = t242 * t563 + t328;
t22 = t146 * t461 - t169 * t463 + t320 * t92 + t324 * t98;
t17 = pkin(10) * t228 + t22;
t155 = -pkin(8) * t242 + t321 * t444 + t325 * t396;
t124 = -t299 * qJ(3) - t301 * qJD(3) - t155;
t95 = -pkin(3) * t242 - t124;
t29 = -pkin(4) * t126 - pkin(10) * t125 + t95;
t81 = t146 * t320 + t169 * t324;
t73 = pkin(10) * t273 + t81;
t239 = t321 * t451 + t293;
t202 = pkin(3) * t432 + t239;
t282 = t301 * qJ(3);
t163 = t282 + t202;
t84 = -pkin(4) * t216 + pkin(10) * t336 + t163;
t33 = t319 * t84 + t323 * t73;
t4 = -qJD(5) * t33 - t17 * t319 + t323 * t29;
t1 = pkin(5) * t123 - qJ(6) * t61 - qJD(6) * t154 + t4;
t543 = t228 / 0.2e1;
t558 = t126 / 0.2e1;
t458 = qJD(5) * t323;
t459 = qJD(5) * t319;
t3 = t323 * t17 + t319 * t29 + t84 * t458 - t459 * t73;
t2 = qJ(6) * t62 + qJD(6) * t153 + t3;
t575 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t651 = t575 + t1 * mrSges(7,1) + t560 * t623 + t566 * t624 + t567 * t626 - t125 * Ifges(5,4) / 0.2e1 + (-t543 - t228 / 0.2e1) * Ifges(5,6) + (-t558 - t126 / 0.2e1) * Ifges(5,2);
t650 = -m(5) - m(7);
t621 = t123 * t624 + t61 * t629 + t62 * t625;
t649 = t621 / 0.2e1;
t622 = t123 * t623 + t61 * t626 + t62 * t624;
t648 = t622 / 0.2e1;
t647 = -mrSges(3,1) + mrSges(4,2);
t646 = t629 * t153;
t484 = t320 * t321;
t218 = (-t319 * t484 + t323 * t325) * t316;
t209 = qJD(1) * t218;
t457 = qJD(5) * t324;
t425 = t323 * t457;
t429 = t319 * t463;
t602 = t209 + t425 - t429;
t480 = t321 * t323;
t219 = (t319 * t325 + t320 * t480) * t316;
t210 = qJD(1) * t219;
t462 = qJD(4) * t323;
t645 = t320 * t462 + t210;
t644 = t629 * t154;
t312 = pkin(5) * t323 + pkin(4);
t318 = -qJ(6) - pkin(10);
t376 = mrSges(5,1) * t320 + mrSges(5,2) * t324;
t643 = -m(7) * (t312 * t320 + t318 * t324) + t324 * mrSges(7,3) - t376;
t559 = t125 / 0.2e1;
t632 = Ifges(5,1) * t559 + Ifges(5,5) * t543;
t568 = Ifges(5,4) * t558 + t632;
t213 = qJD(5) - t216;
t619 = t153 * t624 + t154 * t626 + t213 * t623;
t641 = t619 / 0.2e1;
t618 = t153 * t625 + t213 * t624 + t644;
t617 = t154 * t630 + t626 * t213 + t646;
t290 = pkin(2) * t433;
t353 = pkin(9) * t321 - qJ(3) * t325;
t200 = t353 * t470 + t290;
t115 = t324 * t200 + t320 * t202;
t103 = pkin(10) * t432 + t115;
t381 = pkin(4) * t320 - pkin(10) * t324;
t266 = qJ(3) + t381;
t482 = t320 * t563;
t215 = t319 * t266 - t323 * t482;
t639 = -qJD(5) * t215 + t103 * t319 + t323 * t652;
t460 = qJD(4) * t563;
t426 = t324 * t460;
t638 = t266 * t458 + (-t103 - t426) * t323 + t652 * t319;
t114 = -t320 * t200 + t202 * t324;
t102 = -pkin(4) * t432 - t114;
t427 = t320 * t460;
t637 = -t102 - t427;
t32 = -t319 * t73 + t323 * t84;
t379 = t3 * t323 - t319 * t4;
t636 = -t32 * t458 - t33 * t459 + t379;
t634 = t629 * t323;
t633 = t629 * t319;
t570 = m(7) * pkin(5);
t526 = -mrSges(6,1) - mrSges(7,1);
t631 = mrSges(6,2) + mrSges(7,2);
t524 = mrSges(4,3) - mrSges(3,2);
t628 = Ifges(3,5) - Ifges(4,4);
t627 = Ifges(4,5) - Ifges(3,6);
t394 = t324 * t433;
t405 = t319 * t563 + pkin(5);
t456 = qJD(6) * t323;
t616 = pkin(5) * t394 + (qJD(4) * t405 - t456) * t324 + t639 + (t324 * t459 + t645) * qJ(6);
t20 = -mrSges(6,1) * t62 + mrSges(6,2) * t61;
t88 = mrSges(5,1) * t228 - mrSges(5,3) * t125;
t615 = t88 - t20;
t614 = -qJ(6) * t209 - qJ(6) * t425 + (-qJD(6) * t324 + (qJ(6) * qJD(4) + qJD(5) * t563) * t320) * t319 + t638;
t613 = t273 * Ifges(5,5);
t612 = t273 * Ifges(5,6);
t438 = t319 * t482;
t611 = qJD(5) * t438 + t638;
t610 = t319 * t426 + t639;
t402 = qJD(5) * t318;
t145 = -pkin(4) * t336 - pkin(10) * t216;
t80 = t146 * t324 - t320 * t169;
t47 = t323 * t145 - t319 * t80;
t495 = t216 * t323;
t609 = pkin(5) * t336 + qJ(6) * t495 - qJD(6) * t319 + t323 * t402 - t47;
t48 = t319 * t145 + t323 * t80;
t496 = t216 * t319;
t608 = qJ(6) * t496 + t319 * t402 + t456 - t48;
t607 = -t81 + (t459 - t496) * pkin(5);
t606 = t570 + mrSges(7,1);
t605 = pkin(5) * t602 + t637;
t521 = mrSges(5,3) * t336;
t165 = mrSges(5,1) * t273 + t521;
t86 = -mrSges(6,1) * t153 + mrSges(6,2) * t154;
t604 = t86 - t165;
t144 = -mrSges(5,1) * t216 - mrSges(5,2) * t336;
t399 = mrSges(4,1) * t432;
t234 = -mrSges(4,3) * t301 - t399;
t603 = t144 - t234;
t601 = t319 * t457 + t645;
t398 = mrSges(3,3) * t433;
t400 = mrSges(4,1) * t433;
t600 = t301 * t647 + t398 + t400;
t371 = -mrSges(7,1) * t323 + mrSges(7,2) * t319;
t374 = -mrSges(6,1) * t323 + mrSges(6,2) * t319;
t599 = m(6) * pkin(4) + m(7) * t312 - t371 - t374;
t370 = mrSges(7,1) * t319 + mrSges(7,2) * t323;
t373 = mrSges(6,1) * t319 + mrSges(6,2) * t323;
t72 = -pkin(4) * t273 - t80;
t49 = -pkin(5) * t153 + qJD(6) + t72;
t598 = t49 * t370 + t72 * t373;
t597 = t319 * t626 + t323 * t624;
t596 = -t319 * t624 + t323 * t626;
t595 = t323 * t625 + t633;
t594 = -t319 * t625 + t634;
t593 = t319 * t630 + t634;
t592 = t323 * t630 - t633;
t591 = t394 + t461;
t590 = -m(6) * pkin(10) + m(7) * t318 - mrSges(6,3) - mrSges(7,3);
t23 = -t146 * t463 - t169 * t461 - t320 * t98 + t324 * t92;
t588 = t22 * t320 + t23 * t324;
t455 = m(6) - t650;
t516 = Ifges(3,4) * t321;
t587 = -t321 * (Ifges(3,1) * t325 - t516) / 0.2e1 + pkin(1) * (mrSges(3,1) * t321 + mrSges(3,2) * t325);
t586 = mrSges(5,3) - t647;
t585 = -mrSges(6,1) - t606;
t584 = mrSges(5,1) + t599;
t339 = mrSges(5,2) + t590;
t533 = pkin(5) * t319;
t453 = m(7) * t533;
t581 = pkin(9) * t455 + t453 + t586;
t446 = m(4) + t455;
t580 = qJ(3) * t446 + t524;
t322 = sin(qJ(1));
t478 = t322 * t325;
t326 = cos(qJ(1));
t479 = t321 * t326;
t255 = t317 * t479 + t478;
t475 = t325 * t326;
t481 = t321 * t322;
t254 = -t317 * t475 + t481;
t487 = t316 * t326;
t346 = -t254 * t320 + t324 * t487;
t579 = t255 * t323 + t319 * t346;
t578 = -t255 * t319 + t323 * t346;
t577 = t23 * mrSges(5,1) - t22 * mrSges(5,2) + Ifges(5,5) * t125 + Ifges(5,6) * t126 + Ifges(5,3) * t228;
t24 = -qJ(6) * t154 + t32;
t21 = pkin(5) * t213 + t24;
t25 = qJ(6) * t153 + t33;
t576 = -t1 * t319 + t2 * t323 - t21 * t458 - t25 * t459;
t498 = t324 * mrSges(6,3);
t574 = -m(6) * t266 + t498 - t524 + (-m(4) + t650) * qJ(3) + t643;
t573 = m(5) / 0.2e1;
t572 = m(6) / 0.2e1;
t571 = m(7) / 0.2e1;
t515 = Ifges(5,4) * t336;
t119 = t216 * Ifges(5,2) - t515 + t612;
t561 = -t119 / 0.2e1;
t557 = -t153 / 0.2e1;
t556 = t153 / 0.2e1;
t555 = -t154 / 0.2e1;
t554 = t154 / 0.2e1;
t548 = -t213 / 0.2e1;
t547 = t213 / 0.2e1;
t546 = -t216 / 0.2e1;
t545 = t336 / 0.2e1;
t544 = -t336 / 0.2e1;
t538 = pkin(1) * t316;
t536 = pkin(1) * t325;
t535 = pkin(5) * t154;
t439 = t320 * t488;
t253 = t317 * t324 - t439;
t184 = -t253 * t319 + t316 * t480;
t534 = pkin(5) * t184;
t523 = Ifges(3,4) + Ifges(4,6);
t522 = mrSges(5,3) * t216;
t520 = mrSges(6,3) * t153;
t519 = mrSges(6,3) * t154;
t518 = mrSges(7,3) * t153;
t517 = mrSges(7,3) * t154;
t514 = Ifges(5,4) * t320;
t513 = Ifges(5,4) * t324;
t508 = Ifges(4,6) * t321;
t507 = Ifges(4,6) * t325;
t18 = -pkin(4) * t228 - t23;
t504 = t18 * t324;
t503 = t216 * Ifges(5,6);
t502 = t336 * Ifges(5,5);
t499 = t273 * Ifges(5,3);
t434 = -pkin(2) - t536;
t172 = pkin(3) * t490 + t302 + (-pkin(9) + t434) * t317;
t472 = pkin(2) * t488 + qJ(3) * t490;
t222 = -t472 - t538;
t305 = pkin(9) * t488;
t194 = t222 - t305;
t106 = t320 * t172 + t324 * t194;
t100 = pkin(10) * t490 + t106;
t311 = t321 * t537;
t260 = pkin(8) * t488 + t311;
t221 = -t317 * qJ(3) - t260;
t193 = pkin(3) * t488 - t221;
t252 = t317 * t320 + t324 * t488;
t113 = pkin(4) * t252 - pkin(10) * t253 + t193;
t46 = t323 * t100 + t319 * t113;
t494 = t254 * t319;
t256 = t317 * t478 + t479;
t491 = t256 * t319;
t489 = t316 * t322;
t486 = t319 * t320;
t485 = t319 * t324;
t483 = t320 * t323;
t477 = t323 * t324;
t107 = -mrSges(7,2) * t213 + t518;
t108 = -mrSges(6,2) * t213 + t520;
t474 = t107 + t108;
t109 = mrSges(7,1) * t213 - t517;
t110 = mrSges(6,1) * t213 - t519;
t473 = -t109 - t110;
t471 = t326 * pkin(1) + pkin(8) * t489;
t467 = qJD(1) ^ 2 * t316 ^ 2;
t464 = qJD(4) * t319;
t448 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t447 = -Ifges(3,6) / 0.2e1 + Ifges(4,5) / 0.2e1;
t257 = -t317 * t481 + t475;
t437 = t257 * pkin(2) + t471;
t430 = t316 * t466;
t19 = -t62 * mrSges(7,1) + t61 * mrSges(7,2);
t417 = t470 / 0.2e1;
t415 = -t463 / 0.2e1;
t408 = -t457 / 0.2e1;
t406 = -pkin(1) * t322 + pkin(8) * t487;
t181 = t243 * mrSges(4,1) + t299 * mrSges(4,2);
t45 = -t100 * t319 + t323 * t113;
t105 = t172 * t324 - t320 * t194;
t401 = t562 * t490;
t397 = mrSges(3,3) * t432;
t395 = pkin(3) * t489 + t437;
t390 = t321 * t417;
t384 = -t255 * pkin(2) + t406;
t378 = mrSges(5,1) * t252 + mrSges(5,2) * t253;
t185 = t253 * t323 + t319 * t490;
t375 = mrSges(6,1) * t184 - mrSges(6,2) * t185;
t372 = -t184 * mrSges(7,1) + t185 * mrSges(7,2);
t369 = mrSges(4,2) * t325 - mrSges(4,3) * t321;
t368 = Ifges(5,1) * t320 + t513;
t363 = Ifges(5,2) * t324 + t514;
t358 = Ifges(5,5) * t320 + Ifges(5,6) * t324;
t187 = t256 * t320 + t324 * t489;
t132 = -t187 * t319 + t257 * t323;
t349 = pkin(3) * t487 + t384;
t238 = pkin(8) * t433 - t295;
t296 = t325 * t450;
t240 = -pkin(8) * t431 + t296;
t292 = pkin(2) * t431;
t167 = t292 + (qJD(2) * t353 - t465) * t316;
t203 = (t488 * t562 + t311) * qJD(2);
t53 = -t320 * t167 - t172 * t463 - t194 * t461 + t203 * t324;
t188 = t254 * t324 + t320 * t487;
t342 = t321 * (-Ifges(4,2) * t325 + t508);
t341 = t325 * (Ifges(4,3) * t321 - t507);
t52 = t324 * t167 + t172 * t461 - t194 * t463 + t320 * t203;
t43 = pkin(10) * t430 + t52;
t313 = t317 * qJD(3);
t171 = -qJD(2) * t401 + t296 + t313;
t182 = -qJD(4) * t252 + t320 * t431;
t183 = -qJD(4) * t439 + t317 * t461 - t324 * t431;
t77 = pkin(4) * t183 - pkin(10) * t182 + t171;
t8 = -t100 * t459 + t113 * t458 + t319 * t77 + t323 * t43;
t99 = -pkin(4) * t490 - t105;
t142 = -pkin(2) * t299 + t333;
t331 = t156 * mrSges(3,1) - t155 * mrSges(3,2) + t142 * mrSges(4,2) - t124 * mrSges(4,3);
t44 = -pkin(4) * t430 - t53;
t9 = -qJD(5) * t46 - t319 * t43 + t323 * t77;
t289 = Ifges(3,4) * t432;
t281 = Ifges(4,1) * t299;
t280 = Ifges(3,3) * t299;
t276 = t318 * t323;
t275 = t318 * t319;
t265 = (t563 + t533) * t324;
t263 = t323 * t266;
t259 = t317 * t536 - t302;
t258 = (-mrSges(3,1) * t325 + mrSges(3,2) * t321) * t316;
t246 = t256 * pkin(2);
t244 = t254 * pkin(2);
t241 = t260 * qJD(2);
t237 = -qJ(3) * t432 + t290;
t236 = t369 * t470;
t233 = -mrSges(3,2) * t301 + t397;
t227 = Ifges(4,4) * t243;
t226 = Ifges(3,5) * t243;
t225 = Ifges(4,5) * t242;
t224 = Ifges(3,6) * t242;
t223 = t317 * t434 + t302;
t214 = t263 + t438;
t212 = -t240 - t313;
t211 = Ifges(5,4) * t216;
t208 = t301 * t319 - t323 * t394;
t207 = t301 * t323 + t319 * t394;
t206 = (-pkin(2) * t325 + t407) * t470;
t204 = t292 + (-qJ(3) * t466 - t465) * t316;
t201 = -qJD(1) * t401 + t295;
t199 = -t282 - t239;
t198 = t301 * Ifges(4,4) + (-Ifges(4,2) * t321 - t507) * t470;
t197 = t301 * Ifges(4,5) + (-t325 * Ifges(4,3) - t508) * t470;
t196 = Ifges(3,1) * t433 + Ifges(3,5) * t301 + t289;
t195 = t301 * Ifges(3,6) + (t325 * Ifges(3,2) + t516) * t470;
t192 = -pkin(2) * t301 + qJD(3) + t238;
t186 = -t256 * t324 + t320 * t489;
t180 = mrSges(4,1) * t242 - mrSges(4,3) * t299;
t177 = -qJ(6) * t485 + t215;
t168 = -qJ(6) * t477 + t320 * t405 + t263;
t164 = -mrSges(5,2) * t273 + t522;
t133 = t187 * t323 + t257 * t319;
t127 = pkin(2) * t242 + t328;
t120 = -Ifges(5,1) * t336 + t211 + t613;
t118 = t499 - t502 + t503;
t97 = qJD(5) * t184 + t182 * t323 + t319 * t430;
t96 = -qJD(5) * t185 - t182 * t319 + t323 * t430;
t89 = -mrSges(5,2) * t228 + mrSges(5,3) * t126;
t85 = -mrSges(7,1) * t153 + mrSges(7,2) * t154;
t71 = t99 - t534;
t63 = -mrSges(5,1) * t126 + mrSges(5,2) * t125;
t38 = qJ(6) * t184 + t46;
t37 = -mrSges(6,2) * t123 + mrSges(6,3) * t62;
t36 = -mrSges(7,2) * t123 + mrSges(7,3) * t62;
t35 = mrSges(6,1) * t123 - mrSges(6,3) * t61;
t34 = mrSges(7,1) * t123 - mrSges(7,3) * t61;
t30 = pkin(5) * t252 - qJ(6) * t185 + t45;
t26 = -pkin(5) * t96 + t44;
t7 = -pkin(5) * t62 + qJDD(6) + t18;
t6 = qJ(6) * t96 + qJD(6) * t184 + t8;
t5 = pkin(5) * t183 - qJ(6) * t97 - qJD(6) * t185 + t9;
t10 = [(-t22 * mrSges(5,3) - Ifges(5,4) * t559 + t648 + t651) * t252 + (-mrSges(5,3) * t23 + 0.2e1 * t568) * t253 + (-mrSges(6,3) * t4 - mrSges(7,3) * t1 + 0.2e1 * t640) * t185 + m(4) * (t124 * t221 + t127 * t222 + t142 * t223 + t192 * t241 + t199 * t212 + t204 * t206) + m(7) * (t1 * t30 + t2 * t38 + t21 * t5 + t25 * t6 + t26 * t49 + t7 * t71) + m(6) * (t18 * t99 + t3 * t46 + t32 * t9 + t33 * t8 + t4 * t45 + t44 * t72) + m(5) * (t105 * t23 + t106 * t22 + t163 * t171 + t193 * t95 + t52 * t81 + t53 * t80) + t182 * t120 / 0.2e1 + t32 * (mrSges(6,1) * t183 - mrSges(6,3) * t97) + t21 * (mrSges(7,1) * t183 - mrSges(7,3) * t97) + t33 * (-mrSges(6,2) * t183 + mrSges(6,3) * t96) + t25 * (-mrSges(7,2) * t183 + mrSges(7,3) * t96) + t163 * (mrSges(5,1) * t183 + mrSges(5,2) * t182) + t7 * t372 - t18 * t375 + t95 * t378 + t600 * t241 + t193 * t63 + t216 * (Ifges(5,4) * t182 - Ifges(5,2) * t183) / 0.2e1 + (t226 / 0.2e1 - t224 / 0.2e1 + t280 / 0.2e1 + t281 / 0.2e1 - t227 / 0.2e1 + t225 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t299 + t448 * t243 + t447 * t242 + t331) * t317 + t183 * t561 + (-t182 * t80 - t183 * t81) * mrSges(5,3) + t183 * t641 + m(3) * (t155 * t260 + t156 * t259 + t238 * t241 + t239 * t240) + (-m(3) * t406 - m(4) * t384 + mrSges(2,1) * t322 + mrSges(2,2) * t326 - m(6) * (pkin(4) * t346 + t349) - m(5) * t349 - t346 * mrSges(5,1) - m(7) * (t312 * t346 + t349) + t580 * t254 + t526 * t578 + t631 * t579 + t339 * t188 + t581 * t255) * g(1) + (-m(3) * t471 - m(4) * t437 - mrSges(2,1) * t326 + mrSges(2,2) * t322 - m(6) * (pkin(4) * t187 + t395) - m(5) * t395 - t187 * mrSges(5,1) - m(7) * (t187 * t312 + t395) - t580 * t256 + t526 * t133 - t631 * t132 + t339 * t186 - t581 * t257) * g(2) + t52 * t164 + t53 * t165 + t171 * t144 + (Ifges(5,1) * t182 - Ifges(5,4) * t183) * t544 + t6 * t107 + t8 * t108 + t5 * t109 + t9 * t110 + t105 * t88 + t106 * t89 + t72 * (-mrSges(6,1) * t96 + mrSges(6,2) * t97) + t49 * (-mrSges(7,1) * t96 + mrSges(7,2) * t97) + t99 * t20 + Ifges(2,3) * qJDD(1) + t26 * t85 + t44 * t86 + t71 * t19 + t45 * t35 + t46 * t37 + t38 * t36 + t30 * t34 + ((-mrSges(3,1) * t242 - mrSges(3,2) * t243 + (m(3) * t538 - t258) * qJDD(1)) * pkin(1) + (-t124 * mrSges(4,1) + t127 * mrSges(4,2) + t155 * mrSges(3,3) - t627 * t299 + t523 * t243 + (-Ifges(3,2) - Ifges(4,3)) * t242) * t325 + (-t127 * mrSges(4,3) - t156 * mrSges(3,3) + t142 * mrSges(4,1) + t628 * t299 + (Ifges(3,1) + Ifges(4,2)) * t243 - t523 * t242 + t577) * t321 + ((-t239 * mrSges(3,3) + t199 * mrSges(4,1) - t206 * mrSges(4,2) - t195 / 0.2e1 + t197 / 0.2e1 + t447 * t301) * t321 + (t238 * mrSges(3,3) + t192 * mrSges(4,1) - t81 * mrSges(5,2) + t80 * mrSges(5,1) + t503 / 0.2e1 - t502 / 0.2e1 + t499 / 0.2e1 - t206 * mrSges(4,3) + t196 / 0.2e1 - t198 / 0.2e1 + t118 / 0.2e1 + t448 * t301) * t325 + (-t341 / 0.2e1 + t325 * (Ifges(3,4) * t325 - Ifges(3,2) * t321) / 0.2e1 - t342 / 0.2e1 - t587) * t470) * qJD(2) + (mrSges(3,3) + mrSges(4,1)) * (-g(1) * t326 - g(2) * t322)) * t316 + t273 * (Ifges(5,5) * t182 - Ifges(5,6) * t183) / 0.2e1 + t221 * t180 + t223 * t181 + t212 * t234 + t204 * t236 + t240 * t233 + t617 * t97 / 0.2e1 + t618 * t96 / 0.2e1 + t222 * (-mrSges(4,2) * t242 - mrSges(4,3) * t243) + (t183 * t623 + t624 * t96 + t626 * t97) * t547 + (t183 * t624 + t625 * t96 + t629 * t97) * t556 + (t183 * t626 + t629 * t96 + t630 * t97) * t554 + (t3 * mrSges(6,3) + t2 * mrSges(7,3) + t624 * t560 + t625 * t566 + t629 * t567 + t649) * t184 + t259 * (mrSges(3,1) * t299 - mrSges(3,3) * t243) + t260 * (-mrSges(3,2) * t299 - mrSges(3,3) * t242); (t324 * t619 + t195) * t390 + t513 * t558 + (-t81 * (mrSges(5,3) * t321 * t324 - mrSges(5,2) * t325) - t206 * (-mrSges(4,2) * t321 - mrSges(4,3) * t325) - t80 * (mrSges(5,1) * t325 - mrSges(5,3) * t484)) * t470 + (t342 + t341) * t467 / 0.2e1 + t651 * t320 + (-pkin(2) * t142 - qJ(3) * t124 - qJD(3) * t199 - t206 * t237) * m(4) - t615 * t324 * t563 + (-m(4) * t199 + t233 - t234 - t397) * t238 + (-t3 * t485 - t4 * t477) * mrSges(6,3) + t273 * (mrSges(5,1) * t324 - mrSges(5,2) * t320) * t163 - (t216 * t363 + t273 * t358 - t336 * t368) * qJD(4) / 0.2e1 - ((-Ifges(3,2) * t433 + t118 + t196 + t289) * t325 + t273 * (Ifges(5,3) * t325 + t321 * t358) - t336 * (Ifges(5,5) * t325 + t321 * t368) + t216 * (Ifges(5,6) * t325 + t321 * t363) + (t324 * t119 + t320 * t120 + t197) * t321 + (t321 * t627 + t325 * t628) * t301) * t470 / 0.2e1 + (t4 * t214 + t3 * t215 - (t463 * t72 - t504) * t563 - t102 * t72 + t611 * t33 + t610 * t32) * m(6) + (t95 * qJ(3) - ((-t320 * t80 + t324 * t81) * qJD(4) + t588) * t563 - t114 * t80 - t115 * t81 + (qJD(3) - t201) * t163) * m(5) + (-t115 - t426) * t164 + (t427 - t114) * t165 + t177 * t36 - pkin(2) * t181 + t95 * t376 + t610 * t110 + t611 * t108 + t605 * t85 + (-m(4) * t192 + t398 - t600) * t239 + (mrSges(7,1) * t591 + mrSges(7,3) * t601) * t21 + (mrSges(6,1) * t591 + mrSges(6,3) * t601) * t32 + (mrSges(7,1) * t602 - mrSges(7,2) * t601) * t49 + (mrSges(6,1) * t602 - mrSges(6,2) * t601) * t72 + (-mrSges(6,2) * t591 - mrSges(6,3) * t602) * t33 + (-mrSges(7,2) * t591 - mrSges(7,3) * t602) * t25 + t603 * qJD(3) + t637 * t86 - t201 * t144 + t226 - t227 - t224 + t225 + t614 * t107 + t616 * t109 + (t1 * t168 + t177 * t2 + t21 * t616 + t25 * t614 + t265 * t7 + t49 * t605) * m(7) - t621 * t485 / 0.2e1 + (-t180 + t63) * qJ(3) + (-m(4) * t472 - (m(6) * t381 - t498) * t490 + t258 + t526 * t219 - t631 * t218 - t455 * (t305 + t472) + (t369 + (-t453 - mrSges(5,3)) * t325 + t643 * t321) * t316) * g(3) + (t494 * t570 + m(4) * t244 - t455 * (-pkin(9) * t254 - t244) + t526 * (t255 * t483 - t494) - t631 * (-t254 * t323 - t255 * t486) + t586 * t254 + t574 * t255) * g(2) + (t491 * t570 + m(4) * t246 - t455 * (-pkin(9) * t256 - t246) + t526 * (t257 * t483 - t491) - t631 * (-t256 * t323 - t257 * t486) + t586 * t256 + t574 * t257) * g(1) - t89 * t482 + t280 + t281 + t120 * t415 + t331 + t477 * t640 + t373 * t504 - t514 * t559 + (-t597 * t457 + (-t320 * t596 + t324 * t623) * qJD(4)) * t547 + (-t595 * t457 + (-t320 * t594 + t324 * t624) * qJD(4)) * t556 - t192 * t399 - t199 * t400 + t168 * t34 + (-t1 * t477 - t2 * t485) * mrSges(7,3) + (-t593 * t457 + (-t320 * t592 + t324 * t626) * qJD(4)) * t554 + (t209 * t624 + t210 * t626 - t394 * t623) * t548 + (t209 * t625 + t210 * t629 - t394 * t624) * t557 + (t209 * t629 + t210 * t630 - t394 * t626) * t555 + t320 * t648 + t214 * t35 + t215 * t37 + t325 * t198 * t417 - t237 * t236 + t617 * (t319 * t408 + t323 * t415 - t210 / 0.2e1) + t618 * (t323 * t408 + t429 / 0.2e1 - t209 / 0.2e1) + (-t81 * mrSges(5,3) + t561 + t641) * t461 + t265 * t19 + t587 * t467 + (t463 * t80 - t588) * mrSges(5,3) + (t7 * t370 + t560 * t596 + t566 * t594 + t567 * t592 + t568 + t632) * t324; -t474 * t208 + t473 * t207 + t236 * t433 + (t164 * t433 - t19 + (t319 * t473 + t323 * t474 + t164) * qJD(4) + t615) * t324 + (t89 + (t36 + t37) * t323 + (-t34 - t35) * t319 + (-t319 * t474 + t323 * t473) * qJD(5) + t273 * (t85 + t604)) * t320 - m(7) * (t207 * t21 + t208 * t25) - m(6) * (t207 * t32 + t208 * t33) + 0.2e1 * ((-t21 * t464 + t25 * t462 - t7) * t571 + (-t32 * t464 + t33 * t462 - t18) * t572 + m(5) * t81 * t390 + (qJD(4) * t81 + t23) * t573) * t324 + 0.2e1 * ((qJD(4) * t49 + t576) * t571 + (qJD(4) * t72 + t636) * t572 + (-qJD(4) * t80 + t22) * t573 + (t49 * t571 + t72 * t572 - m(5) * t80 / 0.2e1) * t433) * t320 + t181 + (-g(1) * t256 - g(2) * t254 + g(3) * t488) * t446 + (-m(5) * t163 - t603) * t301 + (t199 * t301 + t206 * t433 + t142) * m(4); (-t459 / 0.2e1 + t496 / 0.2e1) * t618 + (-pkin(4) * t18 - t32 * t47 - t33 * t48) * m(6) + (t153 * t594 + t154 * t592 + t213 * t596) * qJD(5) / 0.2e1 + (t522 - t164) * t80 + t593 * t567 + t595 * t566 + t597 * t560 + t598 * qJD(5) - (-t21 * mrSges(7,1) - t32 * mrSges(6,1) + t25 * mrSges(7,2) + t33 * mrSges(6,2) - Ifges(5,2) * t546 - t163 * mrSges(5,1) + t612 / 0.2e1 + t624 * t557 + t626 * t555 + t623 * t548) * t336 + t7 * t371 + t18 * t374 + (Ifges(5,1) * t545 - t163 * mrSges(5,2) - t613 / 0.2e1 + t594 * t557 + t592 * t555 + t596 * t548 - t598) * t216 + t607 * t85 + t608 * t107 + (t1 * t275 - t2 * t276 + t21 * t609 + t25 * t608 - t312 * t7 + t49 * t607) * m(7) + t609 * t109 + (-m(6) * t72 - t521 - t604) * t81 + (t252 * t599 + t253 * t590 + t378) * g(3) + (t32 * t495 + t33 * t496 + t636) * mrSges(6,3) + (t323 * t37 + m(6) * ((-t319 * t33 - t32 * t323) * qJD(5) + t379) - t319 * t35 - t110 * t458 - t108 * t459) * pkin(10) + (t515 + t619) * t545 + (t211 + t120) * t546 + t319 * t640 + t577 + (t186 * t584 + t187 * t339) * g(1) + (-t188 * t584 - t339 * t346) * g(2) + (t21 * t495 + t25 * t496 + t576) * mrSges(7,3) + (t458 / 0.2e1 - t495 / 0.2e1) * t617 + t119 * t544 - t48 * t108 - t47 * t110 - pkin(4) * t20 + t323 * t649 + t275 * t34 - t276 * t36 - t312 * t19; (-t578 * t631 + t579 * t585) * g(2) + (-t154 * t625 + t617 + t646) * t557 + (t132 * t585 + t133 * t631) * g(1) + (t153 * t626 - t154 * t624) * t548 + (t153 * t630 - t644) * t555 + t575 - t85 * t535 + t21 * t518 - t72 * (mrSges(6,1) * t154 + mrSges(6,2) * t153) - t24 * t107 + t618 * t554 + pkin(5) * t34 + (-m(7) * t535 - mrSges(7,1) * t154 - mrSges(7,2) * t153) * t49 + (t519 + t110) * t33 + (t520 - t108) * t32 + (-m(7) * (-t21 + t24) + t517 + t109) * t25 + (-m(7) * t534 + t372 - t375) * g(3) + t606 * t1 + t622; -t153 * t107 + t154 * t109 + (-g(1) * t186 + g(2) * t188 - g(3) * t252 - t25 * t153 + t21 * t154 + t7) * m(7) + t19;];
tau  = t10;
