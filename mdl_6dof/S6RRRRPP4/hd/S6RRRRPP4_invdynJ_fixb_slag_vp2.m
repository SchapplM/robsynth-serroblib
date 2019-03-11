% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:40
% EndTime: 2019-03-09 20:59:46
% DurationCPUTime: 42.03s
% Computational Cost: add. (18636->901), mult. (41606->1154), div. (0->0), fcn. (29606->14), ass. (0->412)
t379 = sin(qJ(2));
t383 = cos(qJ(2));
t431 = pkin(2) * t379 - pkin(8) * t383;
t316 = t431 * qJD(1);
t382 = cos(qJ(3));
t378 = sin(qJ(3));
t479 = qJD(1) * t379;
t454 = t378 * t479;
t222 = pkin(7) * t454 + t382 * t316;
t491 = t382 * t383;
t409 = pkin(3) * t379 - pkin(9) * t491;
t385 = -pkin(9) - pkin(8);
t456 = qJD(3) * t385;
t683 = -qJD(1) * t409 + t382 * t456 - t222;
t285 = t378 * t316;
t495 = t379 * t382;
t497 = t378 * t383;
t682 = t285 + (-pkin(7) * t495 - pkin(9) * t497) * qJD(1) - t378 * t456;
t377 = sin(qJ(4));
t381 = cos(qJ(4));
t311 = t377 * t382 + t378 * t381;
t614 = qJD(3) + qJD(4);
t216 = t614 * t311;
t402 = t311 * t383;
t249 = qJD(1) * t402;
t681 = t216 - t249;
t671 = mrSges(6,1) + mrSges(7,1);
t670 = mrSges(6,2) - mrSges(7,3);
t337 = t385 * t378;
t338 = t385 * t382;
t220 = t377 * t337 - t381 * t338;
t625 = -qJD(4) * t220 + t682 * t377 + t381 * t683;
t470 = qJD(4) * t381;
t471 = qJD(4) * t377;
t624 = t337 * t470 + t338 * t471 + t377 * t683 - t682 * t381;
t478 = qJD(1) * t383;
t349 = qJD(3) - t478;
t341 = qJD(4) + t349;
t547 = t341 / 0.2e1;
t476 = qJD(2) * t382;
t308 = -t454 + t476;
t452 = t382 * t479;
t309 = qJD(2) * t378 + t452;
t207 = t308 * t377 + t309 * t381;
t376 = sin(pkin(10));
t434 = t381 * t308 - t309 * t377;
t509 = cos(pkin(10));
t653 = t207 * t509 + t376 * t434;
t563 = t653 / 0.2e1;
t642 = Ifges(7,4) + Ifges(6,5);
t644 = Ifges(6,1) + Ifges(7,1);
t680 = t642 * t547 + t563 * t644;
t137 = t207 * t376 - t434 * t509;
t567 = t137 / 0.2e1;
t679 = mrSges(4,3) + mrSges(5,3);
t643 = -Ifges(6,4) + Ifges(7,5);
t678 = t643 + Ifges(7,5);
t469 = qJD(1) * qJD(2);
t320 = qJDD(1) * t383 - t379 * t469;
t303 = t320 * pkin(7);
t677 = t303 * t383;
t411 = t377 * t378 - t381 * t382;
t215 = t614 * t411;
t401 = t411 * t383;
t250 = qJD(1) * t401;
t676 = -pkin(4) * t479 - qJD(5) * t311 + t625 + (t215 - t250) * qJ(5);
t675 = -qJ(5) * t681 - qJD(5) * t411 + t624;
t674 = pkin(5) * t653 + qJ(6) * t137;
t548 = -t341 / 0.2e1;
t305 = qJDD(3) - t320;
t293 = qJDD(4) + t305;
t321 = qJDD(1) * t379 + t383 * t469;
t194 = qJD(3) * t308 + qJDD(2) * t378 + t321 * t382;
t195 = -qJD(3) * t309 + qJDD(2) * t382 - t321 * t378;
t101 = qJD(4) * t434 + t194 * t381 + t195 * t377;
t432 = pkin(2) * t383 + pkin(8) * t379;
t327 = -pkin(1) - t432;
t294 = t327 * qJD(1);
t364 = pkin(7) * t478;
t336 = qJD(2) * pkin(8) + t364;
t210 = t382 * t294 - t336 * t378;
t171 = -pkin(9) * t309 + t210;
t159 = pkin(3) * t349 + t171;
t211 = t294 * t378 + t336 * t382;
t172 = pkin(9) * t308 + t211;
t168 = t381 * t172;
t106 = t159 * t377 + t168;
t508 = qJDD(1) * pkin(1);
t214 = -pkin(2) * t320 - pkin(8) * t321 - t508;
t271 = qJDD(2) * pkin(8) + t303;
t122 = -qJD(3) * t211 + t382 * t214 - t271 * t378;
t87 = pkin(3) * t305 - pkin(9) * t194 + t122;
t472 = qJD(3) * t382;
t474 = qJD(3) * t378;
t121 = t378 * t214 + t382 * t271 + t294 * t472 - t336 * t474;
t96 = pkin(9) * t195 + t121;
t25 = -qJD(4) * t106 - t377 * t96 + t381 * t87;
t13 = pkin(4) * t293 - qJ(5) * t101 - qJD(5) * t207 + t25;
t102 = -qJD(4) * t207 - t194 * t377 + t195 * t381;
t24 = t159 * t470 - t172 * t471 + t377 * t87 + t381 * t96;
t15 = qJ(5) * t102 + qJD(5) * t434 + t24;
t6 = t376 * t13 + t509 * t15;
t2 = qJ(6) * t293 + qJD(6) * t341 + t6;
t5 = t13 * t509 - t376 * t15;
t3 = -t293 * pkin(5) + qJDD(6) - t5;
t599 = -t25 * mrSges(5,1) - t5 * mrSges(6,1) + t3 * mrSges(7,1) + t24 * mrSges(5,2) + t6 * mrSges(6,2) - t2 * mrSges(7,3);
t51 = t101 * t376 - t102 * t509;
t52 = t101 * t509 + t376 * t102;
t613 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t641 = -Ifges(6,6) + Ifges(7,6);
t612 = Ifges(5,5) * t101 + Ifges(5,6) * t102 + t293 * t613 + t51 * t641 + t52 * t642;
t564 = -t653 / 0.2e1;
t568 = -t137 / 0.2e1;
t363 = pkin(7) * t479;
t335 = -qJD(2) * pkin(2) + t363;
t236 = -pkin(3) * t308 + t335;
t154 = -pkin(4) * t434 + qJD(5) + t236;
t166 = t377 * t172;
t105 = t381 * t159 - t166;
t660 = qJ(5) * t207;
t85 = t105 - t660;
t78 = pkin(4) * t341 + t85;
t629 = qJ(5) * t434;
t86 = t106 + t629;
t81 = t509 * t86;
t29 = t376 * t78 + t81;
t27 = qJ(6) * t341 + t29;
t61 = pkin(5) * t137 - qJ(6) * t653 + t154;
t666 = mrSges(6,1) * t154 + mrSges(7,1) * t61 - Ifges(6,4) * t563 - Ifges(7,5) * t564 - Ifges(6,6) * t547 - Ifges(7,6) * t548 - (Ifges(6,2) + Ifges(7,3)) * t568 - mrSges(7,2) * t27 - mrSges(6,3) * t29;
t654 = -Ifges(6,2) * t567 + Ifges(7,3) * t568 + t564 * t643 - t666;
t514 = t376 * t86;
t28 = t509 * t78 - t514;
t26 = -t341 * pkin(5) + qJD(6) - t28;
t665 = mrSges(6,2) * t154 + mrSges(7,2) * t26 - mrSges(6,3) * t28 - mrSges(7,3) * t61 + t567 * t643 + t680;
t656 = Ifges(6,4) * t567 + Ifges(7,5) * t568 + t642 * t548 + t564 * t644 - t665;
t673 = -t656 * t137 + (t548 * t641 + t654) * t653 - t599 + t612;
t582 = t51 / 0.2e1;
t581 = t52 / 0.2e1;
t649 = m(6) + m(7);
t552 = t293 / 0.2e1;
t669 = -mrSges(6,3) - mrSges(7,2);
t453 = t378 * t478;
t622 = -t364 + (-t453 + t474) * pkin(3);
t668 = -m(4) * pkin(8) + m(5) * t385 - t679;
t126 = mrSges(6,1) * t341 - mrSges(6,3) * t653;
t127 = -mrSges(7,1) * t341 + mrSges(7,2) * t653;
t626 = t126 - t127;
t667 = -m(7) * t26 + t626;
t375 = qJ(3) + qJ(4);
t367 = pkin(10) + t375;
t353 = sin(t367);
t354 = cos(t367);
t537 = pkin(3) * t382;
t359 = pkin(2) + t537;
t368 = sin(t375);
t369 = cos(t375);
t427 = -mrSges(4,1) * t382 + mrSges(4,2) * t378;
t663 = m(4) * pkin(2) + m(5) * t359 + t369 * mrSges(5,1) - t368 * mrSges(5,2) - t670 * t353 + t354 * t671 - t427;
t536 = pkin(4) * t207;
t637 = -t376 * t675 + t509 * t676;
t636 = t376 * t676 + t509 * t675;
t661 = m(7) * qJ(6) + mrSges(7,3);
t111 = -t171 * t377 - t168;
t407 = t111 - t629;
t463 = pkin(3) * t470;
t464 = pkin(3) * t471;
t112 = t381 * t171 - t166;
t90 = t112 - t660;
t631 = (t463 - t90) * t509 + (-t407 - t464) * t376;
t623 = pkin(4) * t681 + t622;
t475 = qJD(2) * t383;
t396 = t378 * t475 + t379 * t472;
t380 = sin(qJ(1));
t384 = cos(qJ(1));
t490 = t383 * t384;
t257 = -t368 * t490 + t369 * t380;
t659 = g(1) * t384 + g(2) * t380;
t655 = -Ifges(6,2) * t568 + Ifges(7,3) * t567 + t563 * t643 + t666;
t583 = -t51 / 0.2e1;
t304 = t321 * pkin(7);
t272 = -qJDD(2) * pkin(2) + t304;
t160 = -pkin(3) * t195 + t272;
t70 = -pkin(4) * t102 + qJDD(5) + t160;
t7 = pkin(5) * t51 - qJ(6) * t52 - qJD(6) * t653 + t70;
t651 = mrSges(6,2) * t70 + mrSges(7,2) * t3 - mrSges(6,3) * t5 - mrSges(7,3) * t7 + Ifges(6,4) * t583 + 0.2e1 * t642 * t552 + 0.2e1 * t581 * t644 + t582 * t678;
t650 = Ifges(6,4) * t568 + Ifges(7,5) * t567 + t665 + t680;
t588 = m(5) * pkin(3);
t575 = t101 / 0.2e1;
t574 = t102 / 0.2e1;
t559 = t194 / 0.2e1;
t558 = t195 / 0.2e1;
t551 = t305 / 0.2e1;
t648 = t320 / 0.2e1;
t647 = t321 / 0.2e1;
t366 = pkin(7) * t475;
t645 = -mrSges(3,3) + mrSges(2,2);
t639 = -qJ(6) * t479 + t636;
t638 = pkin(5) * t479 - t637;
t145 = -t215 * t376 + t216 * t509;
t146 = -t215 * t509 - t376 * t216;
t164 = t249 * t509 - t250 * t376;
t165 = -t376 * t249 - t250 * t509;
t203 = t311 * t509 - t376 * t411;
t635 = -qJD(6) * t203 + t623 + (-t146 + t165) * qJ(6) + (t145 - t164) * pkin(5);
t633 = qJD(6) + t631;
t630 = t588 + mrSges(4,1);
t265 = t411 * t379;
t307 = t382 * t327;
t209 = -pkin(9) * t495 + t307 + (-pkin(7) * t378 - pkin(3)) * t383;
t351 = pkin(7) * t491;
t244 = t378 * t327 + t351;
t499 = t378 * t379;
t218 = -pkin(9) * t499 + t244;
t148 = t377 * t209 + t381 * t218;
t621 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t308 + mrSges(4,2) * t309 + mrSges(3,3) * t479;
t620 = mrSges(5,1) * t368 + mrSges(6,1) * t353 + mrSges(5,2) * t369 + mrSges(6,2) * t354;
t298 = Ifges(4,4) * t308;
t183 = t309 * Ifges(4,1) + t349 * Ifges(4,5) + t298;
t362 = Ifges(3,4) * t478;
t619 = Ifges(3,1) * t479 + Ifges(3,5) * qJD(2) + t382 * t183 + t362;
t618 = t304 * t379 + t677;
t617 = t121 * t382 - t122 * t378;
t616 = t669 * t379;
t615 = -m(4) - m(5) - m(3);
t489 = t384 * t353;
t240 = -t380 * t354 + t383 * t489;
t241 = t353 * t380 + t354 * t490;
t258 = t368 * t380 + t369 * t490;
t611 = -t257 * mrSges(5,1) + t258 * mrSges(5,2) + t240 * t671 + t241 * t670;
t493 = t380 * t383;
t238 = t353 * t493 + t354 * t384;
t239 = t354 * t493 - t489;
t255 = t368 * t493 + t369 * t384;
t256 = t368 * t384 - t369 * t493;
t610 = t255 * mrSges(5,1) - t256 * mrSges(5,2) + t238 * t671 + t239 * t670;
t609 = -mrSges(5,1) * t236 + mrSges(5,3) * t106;
t608 = mrSges(5,2) * t236 - mrSges(5,3) * t105;
t607 = t309 * Ifges(4,5) + t207 * Ifges(5,5) + t308 * Ifges(4,6) + Ifges(5,6) * t434 + t349 * Ifges(4,3) + t137 * t641 + t341 * t613 + t642 * t653;
t606 = m(7) * pkin(5) + t671;
t429 = t383 * mrSges(3,1) - mrSges(3,2) * t379;
t605 = t379 * t679 + mrSges(2,1) + t429;
t604 = -mrSges(6,2) + t661;
t601 = -t122 * mrSges(4,1) + t121 * mrSges(4,2);
t523 = Ifges(3,4) * t379;
t419 = t383 * Ifges(3,2) + t523;
t594 = t106 * mrSges(5,2) + t26 * mrSges(7,1) + t29 * mrSges(6,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t419 / 0.2e1 - t105 * mrSges(5,1) - t27 * mrSges(7,3) - t28 * mrSges(6,1);
t591 = mrSges(6,1) * t70 + mrSges(7,1) * t7 - mrSges(7,2) * t2 - mrSges(6,3) * t6 + 0.2e1 * Ifges(7,3) * t582 - t52 * Ifges(6,4) / 0.2e1 - t293 * Ifges(6,6) / 0.2e1 + Ifges(7,6) * t552 + t678 * t581 + (-t583 + t582) * Ifges(6,2);
t585 = Ifges(5,4) * t575 + Ifges(5,2) * t574 + Ifges(5,6) * t552;
t584 = Ifges(5,1) * t575 + Ifges(5,4) * t574 + Ifges(5,5) * t552;
t573 = Ifges(4,1) * t559 + Ifges(4,4) * t558 + Ifges(4,5) * t551;
t519 = Ifges(5,4) * t207;
t130 = Ifges(5,2) * t434 + t341 * Ifges(5,6) + t519;
t572 = -t130 / 0.2e1;
t571 = t130 / 0.2e1;
t198 = Ifges(5,4) * t434;
t131 = t207 * Ifges(5,1) + t341 * Ifges(5,5) + t198;
t570 = -t131 / 0.2e1;
t569 = t131 / 0.2e1;
t556 = -t434 / 0.2e1;
t555 = t434 / 0.2e1;
t554 = -t207 / 0.2e1;
t553 = t207 / 0.2e1;
t549 = t309 / 0.2e1;
t541 = pkin(3) * t309;
t540 = pkin(3) * t377;
t539 = pkin(3) * t378;
t538 = pkin(3) * t381;
t535 = pkin(4) * t368;
t534 = pkin(4) * t376;
t533 = pkin(5) * t353;
t530 = g(3) * t379;
t370 = t379 * pkin(7);
t152 = -qJD(2) * t401 - t216 * t379;
t477 = qJD(2) * t379;
t319 = t431 * qJD(2);
t462 = pkin(7) * t477;
t481 = t382 * t319 + t378 * t462;
t144 = t409 * qJD(2) + (-t351 + (pkin(9) * t379 - t327) * t378) * qJD(3) + t481;
t169 = t378 * t319 + t327 * t472 + (-t379 * t476 - t383 * t474) * pkin(7);
t151 = -pkin(9) * t396 + t169;
t63 = -qJD(4) * t148 + t381 * t144 - t151 * t377;
t31 = pkin(4) * t477 - qJ(5) * t152 + qJD(5) * t265 + t63;
t153 = -qJD(2) * t402 + t265 * t614;
t264 = t311 * t379;
t62 = t377 * t144 + t381 * t151 + t209 * t470 - t218 * t471;
t35 = qJ(5) * t153 - qJD(5) * t264 + t62;
t11 = t376 * t31 + t509 * t35;
t525 = mrSges(5,3) * t434;
t524 = mrSges(5,3) * t207;
t522 = Ifges(3,4) * t383;
t521 = Ifges(4,4) * t378;
t520 = Ifges(4,4) * t382;
t518 = t210 * mrSges(4,3);
t517 = t211 * mrSges(4,3);
t516 = t309 * Ifges(4,4);
t515 = t353 * mrSges(7,1);
t374 = -qJ(5) + t385;
t500 = t374 * t379;
t498 = t378 * t380;
t496 = t378 * t384;
t494 = t379 * t384;
t324 = pkin(4) * t369 + t537;
t315 = pkin(2) + t324;
t291 = t383 * t315;
t323 = t535 + t539;
t301 = t384 * t323;
t147 = t381 * t209 - t218 * t377;
t114 = -pkin(4) * t383 + qJ(5) * t265 + t147;
t123 = -qJ(5) * t264 + t148;
t69 = t376 * t114 + t509 * t123;
t482 = -t383 * t301 + t380 * t324;
t358 = pkin(4) + t538;
t439 = t509 * t377;
t268 = pkin(3) * t439 + t376 * t358;
t322 = pkin(3) * t499 + t370;
t480 = t384 * pkin(1) + t380 * pkin(7);
t473 = qJD(3) * t379;
t458 = Ifges(4,5) * t194 + Ifges(4,6) * t195 + Ifges(4,3) * t305;
t237 = pkin(3) * t396 + t366;
t455 = t509 * pkin(4);
t182 = t308 * Ifges(4,2) + t349 * Ifges(4,6) + t516;
t449 = -t378 * t182 / 0.2e1;
t21 = t51 * mrSges(6,1) + t52 * mrSges(6,2);
t20 = t51 * mrSges(7,1) - t52 * mrSges(7,3);
t441 = t661 * t354 * t379;
t41 = -t293 * mrSges(7,1) + t52 * mrSges(7,2);
t438 = t469 / 0.2e1;
t437 = -t238 * pkin(5) + qJ(6) * t239;
t436 = -t240 * pkin(5) + qJ(6) * t241;
t435 = -t323 * t493 - t324 * t384;
t219 = t381 * t337 + t338 * t377;
t212 = pkin(4) * t264 + t322;
t163 = t541 + t536;
t428 = mrSges(3,1) * t379 + mrSges(3,2) * t383;
t426 = mrSges(4,1) * t378 + mrSges(4,2) * t382;
t421 = Ifges(4,1) * t382 - t521;
t420 = Ifges(4,1) * t378 + t520;
t418 = -Ifges(4,2) * t378 + t520;
t417 = Ifges(4,2) * t382 + t521;
t416 = Ifges(3,5) * t383 - Ifges(3,6) * t379;
t415 = Ifges(4,5) * t382 - Ifges(4,6) * t378;
t414 = Ifges(4,5) * t378 + Ifges(4,6) * t382;
t413 = pkin(5) * t354 + qJ(6) * t353;
t412 = t383 * t359 - t379 * t385;
t410 = t257 * pkin(4);
t251 = pkin(4) * t411 - t359;
t128 = -pkin(4) * t153 + t237;
t408 = pkin(1) * t428;
t282 = -t378 * t490 + t380 * t382;
t280 = t378 * t493 + t382 * t384;
t406 = -qJ(5) * t311 + t219;
t10 = t31 * t509 - t376 * t35;
t405 = t335 * t426;
t404 = t379 * (Ifges(3,1) * t383 - t523);
t68 = t114 * t509 - t376 * t123;
t267 = t358 * t509 - t376 * t540;
t398 = t255 * pkin(4);
t397 = -t378 * t473 + t382 * t475;
t394 = Ifges(4,5) * t379 + t383 * t421;
t393 = Ifges(4,6) * t379 + t383 * t418;
t392 = Ifges(4,3) * t379 + t383 * t415;
t372 = t384 * pkin(7);
t355 = -t455 - pkin(5);
t352 = qJ(6) + t534;
t329 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t478;
t295 = t426 * t379;
t283 = t382 * t490 + t498;
t281 = -t380 * t491 + t496;
t260 = -pkin(5) - t267;
t259 = qJ(6) + t268;
t243 = -pkin(7) * t497 + t307;
t233 = mrSges(4,1) * t349 - mrSges(4,3) * t309;
t232 = -mrSges(4,2) * t349 + mrSges(4,3) * t308;
t223 = -pkin(7) * t452 + t285;
t202 = t311 * t376 + t411 * t509;
t180 = -qJ(5) * t411 + t220;
t176 = -t376 * t264 - t265 * t509;
t175 = t264 * t509 - t265 * t376;
t174 = mrSges(5,1) * t341 - t524;
t173 = -mrSges(5,2) * t341 + t525;
t170 = -qJD(3) * t244 + t481;
t162 = -mrSges(4,2) * t305 + mrSges(4,3) * t195;
t161 = mrSges(4,1) * t305 - mrSges(4,3) * t194;
t143 = -mrSges(5,1) * t434 + mrSges(5,2) * t207;
t134 = -mrSges(4,1) * t195 + mrSges(4,2) * t194;
t125 = -mrSges(6,2) * t341 - mrSges(6,3) * t137;
t124 = -mrSges(7,2) * t137 + mrSges(7,3) * t341;
t120 = t180 * t509 + t376 * t406;
t119 = t180 * t376 - t406 * t509;
t118 = pkin(5) * t202 - qJ(6) * t203 + t251;
t115 = t194 * Ifges(4,4) + t195 * Ifges(4,2) + t305 * Ifges(4,6);
t97 = pkin(5) * t175 - qJ(6) * t176 + t212;
t94 = t152 * t509 + t376 * t153;
t93 = t152 * t376 - t153 * t509;
t92 = -mrSges(5,2) * t293 + mrSges(5,3) * t102;
t91 = mrSges(5,1) * t293 - mrSges(5,3) * t101;
t80 = mrSges(6,1) * t137 + mrSges(6,2) * t653;
t79 = mrSges(7,1) * t137 - mrSges(7,3) * t653;
t67 = t536 + t674;
t66 = t383 * pkin(5) - t68;
t65 = -qJ(6) * t383 + t69;
t64 = t163 + t674;
t53 = -mrSges(5,1) * t102 + mrSges(5,2) * t101;
t40 = mrSges(6,1) * t293 - mrSges(6,3) * t52;
t39 = -mrSges(6,2) * t293 - mrSges(6,3) * t51;
t38 = -mrSges(7,2) * t51 + mrSges(7,3) * t293;
t34 = t509 * t85 - t514;
t33 = t376 * t85 + t81;
t22 = pkin(5) * t93 - qJ(6) * t94 - qJD(6) * t176 + t128;
t9 = -pkin(5) * t477 - t10;
t8 = qJ(6) * t477 - qJD(6) * t383 + t11;
t1 = [t650 * t94 + (-t105 * t152 + t106 * t153 - t24 * t264 + t25 * t265) * mrSges(5,3) + (-Ifges(5,1) * t265 - Ifges(5,4) * t264) * t575 + (-Ifges(5,4) * t265 - Ifges(5,2) * t264) * t574 + t160 * (mrSges(5,1) * t264 - mrSges(5,2) * t265) - (t382 * t182 + t378 * t183) * t473 / 0.2e1 + t429 * t508 + (-Ifges(5,5) * t265 - Ifges(5,6) * t264 + t175 * t641 - t383 * t613) * t552 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t618) + t651 * t176 + (Ifges(5,5) * t152 + Ifges(5,6) * t153 + t641 * t93) * t547 + (t607 / 0.2e1 + t210 * mrSges(4,1) - t211 * mrSges(4,2) + Ifges(5,5) * t553 + Ifges(5,6) * t555 + Ifges(6,6) * t568 + Ifges(7,6) * t567 + t563 * t642 - t594 + t613 * t547) * t477 + t522 * t647 + t419 * t648 + t591 * t175 + t655 * t93 + (qJD(2) * t394 - t420 * t473) * t549 + (t619 / 0.2e1 + t449) * t475 + (-mrSges(3,1) * t370 + Ifges(3,5) * t379 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * t383) * qJDD(2) + t152 * t569 + t335 * (mrSges(4,1) * t396 + mrSges(4,2) * t397) - (t458 + t612) * t383 / 0.2e1 + t404 * t438 + (Ifges(5,1) * t152 + Ifges(5,4) * t153) * t553 + t134 * t370 + t153 * t571 + t495 * t573 - t265 * t584 - t264 * t585 + m(7) * (t2 * t65 + t22 * t61 + t26 * t9 + t27 * t8 + t3 * t66 + t7 * t97) + m(6) * (t10 * t28 + t11 * t29 + t128 * t154 + t212 * t70 + t5 * t68 + t6 * t69) + m(5) * (t105 * t63 + t106 * t62 + t147 * t25 + t148 * t24 + t160 * t322 + t236 * t237) + t621 * t366 + (-t121 * t499 - t122 * t495 - t210 * t397 - t211 * t396) * mrSges(4,3) + qJD(2) ^ 2 * t416 / 0.2e1 + (Ifges(5,4) * t152 + Ifges(5,2) * t153) * t555 + Ifges(2,3) * qJDD(1) + t308 * (qJD(2) * t393 - t417 * t473) / 0.2e1 + t349 * (qJD(2) * t392 - t414 * t473) / 0.2e1 - t408 * t469 - t329 * t462 - pkin(1) * (-mrSges(3,1) * t320 + mrSges(3,2) * t321) + t322 * t53 + t272 * t295 + t243 * t161 + t244 * t162 + t170 * t233 + t236 * (-mrSges(5,1) * t153 + mrSges(5,2) * t152) + t237 * t143 + t169 * t232 + t212 * t21 + t62 * t173 + t63 * t174 + (-t498 * t588 - t283 * mrSges(4,1) - t258 * mrSges(5,1) - t282 * mrSges(4,2) - t257 * mrSges(5,2) + t669 * t494 + t615 * t480 - t649 * (t315 * t490 + t380 * t323 - t374 * t494 + t480) + t645 * t380 - t606 * t241 - t604 * t240 + (-m(4) * t432 - m(5) * t412 - t605) * t384) * g(2) + t65 * t38 + t66 * t41 + t68 * t40 + t69 * t39 + m(4) * (t121 * t244 + t122 * t243 + t169 * t211 + t170 * t210 + t335 * t366) + t22 * t79 + (t321 * t370 + t618 + t677) * mrSges(3,3) + (t599 - t642 * t581 + Ifges(3,4) * t647 + Ifges(3,2) * t648 - Ifges(6,6) * t583 - Ifges(4,3) * t551 - Ifges(4,6) * t558 - Ifges(4,5) * t559 + (-Ifges(3,2) * t379 + t522) * t438 - Ifges(5,6) * t574 - Ifges(5,5) * t575 - Ifges(7,6) * t582 + t601) * t383 + (m(4) * t272 * pkin(7) + Ifges(3,1) * t321 + Ifges(3,4) * t648 + t415 * t551 + t418 * t558 + t421 * t559) * t379 + (-t496 * t588 - t281 * mrSges(4,1) - t256 * mrSges(5,1) - t280 * mrSges(4,2) - t255 * mrSges(5,2) - t649 * (t380 * t500 + t301 + t372) + t645 * t384 + t615 * t372 + t606 * t239 + t604 * t238 + (-m(4) * t327 - m(5) * (-pkin(1) - t412) + m(3) * pkin(1) - t649 * (-pkin(1) - t291) + t605 - t616) * t380) * g(1) + t97 * t20 - t115 * t499 / 0.2e1 + t8 * t124 + t11 * t125 + t10 * t126 + t9 * t127 + t128 * t80 + t147 * t91 + t148 * t92; t650 * t146 + (t308 * t418 + t309 * t421 + t349 * t415) * qJD(3) / 0.2e1 - (t308 * t393 + t309 * t394 + t349 * t392) * qJD(1) / 0.2e1 + (-t404 / 0.2e1 + t408) * qJD(1) ^ 2 + (Ifges(5,5) * t311 - Ifges(5,6) * t411) * t552 + (Ifges(5,4) * t311 - Ifges(5,2) * t411) * t574 + (Ifges(5,1) * t311 - Ifges(5,4) * t411) * t575 + t160 * (mrSges(5,1) * t411 + mrSges(5,2) * t311) - t411 * t585 + (m(4) * ((-t210 * t382 - t211 * t378) * qJD(3) + t617) + t382 * t162 - t233 * t472 - t232 * t474 - t378 * t161) * pkin(8) + t617 * mrSges(4,3) - (-Ifges(3,2) * t479 + t362 + t619) * t478 / 0.2e1 + (-t518 + t183 / 0.2e1) * t472 + (-t429 - t649 * (t291 - t500) + t668 * t379 + (-m(7) * t413 - t663) * t383 + t616) * g(3) + (-t211 * (-mrSges(4,2) * t379 - mrSges(4,3) * t497) - t210 * (mrSges(4,1) * t379 - mrSges(4,3) * t491)) * qJD(1) + (t405 + t449) * qJD(3) + t651 * t203 - (Ifges(5,1) * t553 + Ifges(5,4) * t555 + Ifges(5,5) * t547 + t569 + t608) * t215 - (Ifges(5,4) * t553 + Ifges(5,2) * t555 + Ifges(5,6) * t547 + t571 + t609) * t216 + (-t105 * t250 + t106 * t249 - t24 * t411 - t25 * t311) * mrSges(5,3) + (-Ifges(5,4) * t250 - Ifges(5,2) * t249) * t556 + (-Ifges(5,5) * t250 - Ifges(5,6) * t249 + t164 * t641 + t479 * t613) * t548 + (-Ifges(5,1) * t250 - Ifges(5,4) * t249) * t554 - t236 * (mrSges(5,1) * t249 - mrSges(5,2) * t250) + (t41 - t40) * t119 + (t39 + t38) * t120 + (-pkin(2) * t272 - t210 * t222 - t211 * t223 - t335 * t364) * m(4) + t654 * t164 + (t547 * t641 + t655) * t145 + t414 * t551 - t474 * t517 + t417 * t558 + t420 * t559 - t250 * t570 - t249 * t572 + t378 * t573 + t311 * t584 + t272 * t427 - t607 * t479 / 0.2e1 + t656 * t165 - t621 * t364 + t622 * t143 + t623 * t80 + t624 * t173 + (t105 * t625 + t106 * t624 - t160 * t359 + t219 * t25 + t220 * t24 + t236 * t622) * m(5) + t625 * t174 - t405 * t478 - t416 * t469 / 0.2e1 + t182 * t453 / 0.2e1 + t382 * t115 / 0.2e1 + Ifges(3,3) * qJDD(2) - t359 * t53 + Ifges(3,6) * t320 + Ifges(3,5) * t321 - t303 * mrSges(3,2) - t304 * mrSges(3,1) + t251 * t21 - t222 * t233 - t223 * t232 + t219 * t91 + t220 * t92 + t635 * t79 + t636 * t125 + t637 * t126 + (-t119 * t5 + t120 * t6 + t154 * t623 + t251 * t70 + t28 * t637 + t29 * t636) * m(6) + t638 * t127 + t639 * t124 + (t118 * t7 + t119 * t3 + t120 * t2 + t26 * t638 + t27 * t639 + t61 * t635) * m(7) + (t552 * t641 + t591) * t202 + t659 * (t428 + (t374 * t649 + t668 + t669) * t383 + (m(6) * t315 - m(7) * (-t315 - t413) + t663) * t379) + (Ifges(5,5) * t554 + Ifges(5,6) * t556 + Ifges(6,6) * t567 + Ifges(7,6) * t568 + t564 * t642 + t594) * t479 + t329 * t363 + t118 * t20 - pkin(2) * t134; (-(m(7) * (-t323 - t533) - t515) * t379 - t441 + t295) * g(3) - (-Ifges(4,2) * t309 + t183 + t298) * t308 / 0.2e1 - (Ifges(5,4) * t554 + Ifges(5,2) * t556 + Ifges(5,6) * t548 + t572 - t609) * t207 + (-m(6) * t28 - t667) * (-t376 * t90 + t407 * t509 + (t376 * t381 + t439) * qJD(4) * pkin(3)) + (-t154 * t163 + t267 * t5 + t268 * t6 + t29 * t631 + t323 * t530) * m(6) + (Ifges(5,1) * t554 + Ifges(5,4) * t556 + Ifges(5,5) * t548 + t570 - t608) * t434 + t182 * t549 + t309 * t517 + t308 * t518 - t309 * (Ifges(4,1) * t308 - t516) / 0.2e1 - t143 * t541 - m(5) * (t105 * t111 + t106 * t112 + t236 * t541) + (-m(7) * (t435 + t437) - m(6) * t435 - mrSges(4,2) * t281 + t630 * t280 + t610) * g(2) + (-m(7) * (t436 + t482) - m(6) * t482 + mrSges(4,2) * t283 - t630 * t282 + t611) * g(1) + t631 * t125 + t633 * t124 + (-t464 - t111) * t174 + (t2 * t259 + t260 * t3 + t27 * t633 - t61 * t64) * m(7) + t91 * t538 + t92 * t540 + t673 + (m(5) * t539 + t620) * t530 - t601 + (-t112 + t463) * t173 + (t24 * t377 + t25 * t381 + (-t105 * t377 + t106 * t381) * qJD(4)) * t588 - t349 * (Ifges(4,5) * t308 - Ifges(4,6) * t309) / 0.2e1 - t335 * (mrSges(4,1) * t309 + mrSges(4,2) * t308) + t267 * t40 + t268 * t39 + t259 * t38 + t260 * t41 + t211 * t233 - t210 * t232 + t458 - t64 * t79 - t163 * t80; (Ifges(5,5) * t434 - Ifges(5,6) * t207) * t548 + (Ifges(5,1) * t434 - t519) * t554 - t236 * (mrSges(5,1) * t207 + mrSges(5,2) * t434) + t667 * t33 + (-m(7) * (-t398 + t437) + m(6) * t398 + t610) * g(2) + (-m(7) * (t410 + t436) - m(6) * t410 + t611) * g(1) + t40 * t455 + t130 * t553 + (t525 - t173) * t105 - g(3) * ((m(7) * (-t533 - t535) - t515) * t379 + t441) - t80 * t536 + (-Ifges(5,2) * t207 + t131 + t198) * t556 + t39 * t534 + m(7) * (qJD(6) * t27 + t2 * t352 + t3 * t355) + t673 - m(7) * (t27 * t34 + t61 * t67) + t620 * t530 + (t535 * t530 - t154 * t536 + t28 * t33 - t29 * t34 + (t376 * t6 + t5 * t509) * pkin(4)) * m(6) + (t524 + t174) * t106 + (qJD(6) - t34) * t124 + t352 * t38 + t355 * t41 - t67 * t79 - t34 * t125; -(-t124 - t125) * t137 + t626 * t653 + t20 + t21 + (t383 * g(3) - t379 * t659) * t649 + (t137 * t27 - t26 * t653 + t7) * m(7) + (t137 * t29 + t28 * t653 + t70) * m(6); -t341 * t124 + t653 * t79 + (-g(1) * t240 - g(2) * t238 - t27 * t341 - t353 * t530 + t61 * t653 + t3) * m(7) + t41;];
tau  = t1;
