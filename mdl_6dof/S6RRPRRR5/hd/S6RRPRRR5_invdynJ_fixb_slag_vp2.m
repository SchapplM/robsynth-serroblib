% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:05
% EndTime: 2019-03-09 13:39:37
% DurationCPUTime: 52.02s
% Computational Cost: add. (31684->1075), mult. (86686->1464), div. (0->0), fcn. (71696->16), ass. (0->476)
t394 = sin(qJ(2));
t399 = cos(qJ(2));
t537 = sin(pkin(12));
t538 = cos(pkin(12));
t343 = t394 * t537 - t399 * t538;
t389 = sin(pkin(6));
t651 = t389 * t343;
t307 = qJD(1) * t651;
t638 = t307 + qJD(4);
t390 = cos(pkin(6));
t575 = pkin(1) * t390;
t379 = t399 * t575;
t372 = qJD(1) * t379;
t567 = pkin(8) + qJ(3);
t469 = t567 * t394;
t450 = t389 * t469;
t294 = -qJD(1) * t450 + t372;
t525 = t390 * t394;
t378 = pkin(1) * t525;
t527 = t389 * t399;
t295 = (t527 * t567 + t378) * qJD(1);
t460 = t538 * t295;
t212 = t294 * t537 + t460;
t393 = sin(qJ(4));
t398 = cos(qJ(4));
t694 = -t212 + t638 * (pkin(4) * t393 - pkin(10) * t398);
t410 = t394 * t538 + t399 * t537;
t506 = qJD(1) * t389;
t308 = t410 * t506;
t397 = cos(qJ(5));
t392 = sin(qJ(5));
t523 = t392 * t398;
t236 = t307 * t523 + t308 * t397;
t501 = qJD(4) * t398;
t693 = t392 * t501 + t236;
t498 = qJD(5) * t397;
t647 = t393 * t498 + t693;
t282 = t537 * t295;
t213 = t294 * t538 - t282;
t479 = t394 * t506;
t455 = pkin(2) * t479;
t229 = pkin(3) * t308 + pkin(9) * t307 + t455;
t133 = t398 * t213 + t393 * t229;
t115 = pkin(10) * t308 + t133;
t481 = t538 * pkin(2);
t383 = -t481 - pkin(3);
t341 = -t398 * pkin(4) - t393 * pkin(10) + t383;
t480 = t537 * pkin(2);
t382 = t480 + pkin(9);
t500 = qJD(5) * t392;
t502 = qJD(4) * t393;
t657 = -t397 * t115 + t341 * t498 + (-t397 * t502 - t398 * t500) * t382 + t694 * t392;
t475 = t382 * t502;
t692 = t694 * t397 + (t115 + t475) * t392;
t375 = qJD(1) * t390 + qJD(2);
t274 = pkin(2) * t375 + t294;
t191 = t537 * t274 + t460;
t180 = pkin(9) * t375 + t191;
t385 = pkin(2) * t399 + pkin(1);
t331 = -t385 * t506 + qJD(3);
t205 = pkin(3) * t307 - pkin(9) * t308 + t331;
t113 = t180 * t398 + t205 * t393;
t303 = t343 * t506 + qJD(4);
t100 = pkin(10) * t303 + t113;
t190 = t274 * t538 - t282;
t179 = -t375 * pkin(3) - t190;
t272 = -t393 * t308 + t375 * t398;
t273 = t308 * t398 + t375 * t393;
t107 = -t272 * pkin(4) - t273 * pkin(10) + t179;
t55 = -t100 * t392 + t397 * t107;
t56 = t100 * t397 + t107 * t392;
t691 = t55 * mrSges(6,1) - t56 * mrSges(6,2) - t113 * mrSges(5,3);
t496 = qJD(1) * qJD(2);
t327 = (qJDD(1) * t399 - t394 * t496) * t389;
t328 = (qJDD(1) * t394 + t399 * t496) * t389;
t251 = t327 * t538 - t328 * t537;
t250 = qJDD(4) - t251;
t507 = pkin(8) * t527 + t378;
t322 = t507 * qJD(2);
t494 = qJDD(1) * t390;
t487 = pkin(1) * t494;
t370 = t399 * t487;
t374 = qJDD(2) + t494;
t529 = t389 * t394;
t476 = qJD(3) * t529;
t495 = qJDD(1) * t389;
t486 = pkin(8) * t495;
t186 = -t394 * t486 + pkin(2) * t374 - qJ(3) * t328 + t370 + (-t322 - t476) * qJD(1);
t493 = qJD(2) * t575;
t451 = qJD(1) * t493;
t484 = t394 * t487 + (t451 + t486) * t399;
t503 = qJD(3) * t399;
t504 = qJD(2) * t394;
t198 = qJ(3) * t327 + (-pkin(8) * t504 + t503) * t506 + t484;
t111 = t537 * t186 + t538 * t198;
t105 = pkin(9) * t374 + t111;
t252 = t327 * t537 + t328 * t538;
t290 = -pkin(1) * t495 - pkin(2) * t327 + qJDD(3);
t140 = -pkin(3) * t251 - pkin(9) * t252 + t290;
t44 = t398 * t105 + t393 * t140 - t180 * t502 + t205 * t501;
t39 = pkin(10) * t250 + t44;
t110 = t186 * t538 - t537 * t198;
t104 = -t374 * pkin(3) - t110;
t165 = t272 * qJD(4) + t252 * t398 + t374 * t393;
t166 = -t393 * t252 - t308 * t501 + t374 * t398 - t375 * t502;
t53 = -t166 * pkin(4) - t165 * pkin(10) + t104;
t12 = -qJD(5) * t56 - t39 * t392 + t397 * t53;
t666 = t12 * mrSges(6,1);
t11 = -t100 * t500 + t107 * t498 + t397 * t39 + t392 * t53;
t667 = t11 * mrSges(6,2);
t690 = t666 - t667;
t270 = qJD(5) - t272;
t264 = qJD(6) + t270;
t589 = t264 / 0.2e1;
t195 = -t273 * t392 + t303 * t397;
t196 = t273 * t397 + t303 * t392;
t391 = sin(qJ(6));
t396 = cos(qJ(6));
t127 = t195 * t391 + t196 * t396;
t602 = t127 / 0.2e1;
t458 = t396 * t195 - t196 * t391;
t604 = t458 / 0.2e1;
t550 = t273 * Ifges(5,4);
t685 = t303 * Ifges(5,6);
t686 = t272 * Ifges(5,2);
t155 = t550 + t685 + t686;
t663 = t196 * Ifges(6,5) + t127 * Ifges(7,5) + t195 * Ifges(6,6) + Ifges(7,6) * t458 + t270 * Ifges(6,3) + t264 * Ifges(7,3);
t688 = t155 / 0.2e1 - t663 / 0.2e1;
t689 = Ifges(7,5) * t602 + Ifges(7,6) * t604 + Ifges(7,3) * t589 - t688 + t691;
t442 = mrSges(5,1) * t398 - mrSges(5,2) * t393;
t401 = -pkin(11) - pkin(10);
t643 = -m(6) * pkin(10) + m(7) * t401 - mrSges(6,3) - mrSges(7,3);
t384 = pkin(5) * t397 + pkin(4);
t388 = qJ(5) + qJ(6);
t386 = sin(t388);
t387 = cos(t388);
t440 = -mrSges(6,1) * t397 + mrSges(6,2) * t392;
t644 = m(6) * pkin(4) + m(7) * t384 + mrSges(7,1) * t387 - mrSges(7,2) * t386 - t440;
t629 = t393 * t643 - t398 * t644 - mrSges(4,1) - t442;
t112 = -t393 * t180 + t205 * t398;
t553 = t112 * mrSges(5,3);
t687 = t179 * mrSges(5,2);
t517 = t397 * t398;
t237 = -t307 * t517 + t308 * t392;
t348 = t382 * t517;
t533 = t307 * t393;
t684 = pkin(5) * t533 + pkin(11) * t237 + (pkin(5) * t393 - pkin(11) * t517) * qJD(4) + (-t348 + (pkin(11) * t393 - t341) * t392) * qJD(5) + t692;
t683 = -pkin(11) * t647 + t657;
t482 = qJD(5) * t401;
t535 = t272 * t392;
t185 = pkin(4) * t273 - pkin(10) * t272;
t82 = t397 * t112 + t392 * t185;
t682 = pkin(11) * t535 + t392 * t482 - t82;
t534 = t272 * t397;
t81 = -t112 * t392 + t397 * t185;
t681 = -pkin(5) * t273 + pkin(11) * t534 + t397 * t482 - t81;
t678 = t500 - t535;
t588 = -t270 / 0.2e1;
t590 = -t264 / 0.2e1;
t594 = -t196 / 0.2e1;
t596 = -t195 / 0.2e1;
t603 = -t127 / 0.2e1;
t605 = -t458 / 0.2e1;
t677 = Ifges(6,5) * t594 + Ifges(7,5) * t603 + Ifges(6,6) * t596 + Ifges(7,6) * t605 + Ifges(6,3) * t588 + Ifges(7,3) * t590 - t691;
t48 = -pkin(11) * t196 + t55;
t41 = pkin(5) * t270 + t48;
t49 = pkin(11) * t195 + t56;
t546 = t391 * t49;
t16 = t396 * t41 - t546;
t164 = qJDD(5) - t166;
t85 = qJD(5) * t195 + t165 * t397 + t250 * t392;
t6 = pkin(5) * t164 - pkin(11) * t85 + t12;
t86 = -qJD(5) * t196 - t165 * t392 + t250 * t397;
t7 = pkin(11) * t86 + t11;
t2 = qJD(6) * t16 + t391 * t6 + t396 * t7;
t542 = t396 * t49;
t17 = t391 * t41 + t542;
t3 = -qJD(6) * t17 - t391 * t7 + t396 * t6;
t676 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t675 = mrSges(5,1) * t179 + mrSges(7,1) * t16 - mrSges(7,2) * t17;
t591 = t250 / 0.2e1;
t597 = t166 / 0.2e1;
t598 = t165 / 0.2e1;
t612 = Ifges(5,1) * t598 + Ifges(5,4) * t597 + Ifges(5,5) * t591;
t159 = qJDD(6) + t164;
t600 = t159 / 0.2e1;
t34 = -qJD(6) * t127 - t391 * t85 + t396 * t86;
t619 = t34 / 0.2e1;
t33 = qJD(6) * t458 + t391 * t86 + t396 * t85;
t620 = t33 / 0.2e1;
t621 = Ifges(7,1) * t620 + Ifges(7,4) * t619 + Ifges(7,5) * t600;
t623 = Ifges(7,4) * t620 + Ifges(7,2) * t619 + Ifges(7,6) * t600;
t35 = Ifges(6,5) * t85 + Ifges(6,6) * t86 + Ifges(6,3) * t164;
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t159;
t672 = t35 + t8;
t144 = t236 * t391 + t237 * t396;
t346 = t391 * t397 + t392 * t396;
t637 = qJD(5) + qJD(6);
t286 = t637 * t346;
t426 = t391 * t392 - t396 * t397;
t214 = -t286 * t393 - t426 * t501;
t649 = t214 - t144;
t143 = t236 * t396 - t237 * t391;
t324 = t426 * t393;
t215 = t324 * t637 - t346 * t501;
t648 = t215 - t143;
t438 = t386 * mrSges(7,1) + t387 * mrSges(7,2);
t541 = t397 * mrSges(6,2);
t439 = t392 * mrSges(6,1) + t541;
t568 = mrSges(5,3) - mrSges(4,2);
t572 = pkin(5) * t392;
t671 = -m(7) * (pkin(9) + t572) - t438 - m(6) * pkin(9) - t439 - t568;
t132 = -t393 * t213 + t229 * t398;
t114 = -pkin(4) * t308 - t132;
t474 = t382 * t501;
t670 = t474 - t114;
t400 = cos(qJ(1));
t516 = t399 * t400;
t395 = sin(qJ(1));
t521 = t394 * t395;
t669 = t390 * t516 - t521;
t642 = t502 + t533;
t611 = t85 / 0.2e1;
t610 = t86 / 0.2e1;
t668 = -m(6) - m(7);
t599 = t164 / 0.2e1;
t326 = t397 * t341;
t522 = t393 * t397;
t253 = -pkin(11) * t522 + t326 + (-t382 * t392 - pkin(5)) * t398;
t281 = t392 * t341 + t348;
t524 = t392 * t393;
t271 = -pkin(11) * t524 + t281;
t168 = t253 * t391 + t271 * t396;
t665 = -qJD(6) * t168 - t391 * t683 + t396 * t684;
t167 = t253 * t396 - t271 * t391;
t664 = qJD(6) * t167 + t391 * t684 + t396 * t683;
t662 = t112 * mrSges(5,1);
t661 = t113 * mrSges(5,2);
t622 = m(7) * pkin(5);
t659 = -mrSges(6,1) - t622;
t117 = mrSges(5,1) * t250 - mrSges(5,3) * t165;
t42 = -mrSges(6,1) * t86 + mrSges(6,2) * t85;
t658 = -t117 + t42;
t656 = -qJD(5) * t281 + t692;
t362 = t401 * t392;
t363 = t401 * t397;
t299 = t362 * t396 + t363 * t391;
t655 = qJD(6) * t299 + t391 * t681 + t396 * t682;
t300 = t362 * t391 - t363 * t396;
t654 = -qJD(6) * t300 - t391 * t682 + t396 * t681;
t653 = pkin(5) * t647 + t670;
t652 = pkin(5) * t678 - t113;
t291 = pkin(2) * t390 + t379 - t450;
t304 = qJ(3) * t527 + t507;
t225 = t537 * t291 + t538 * t304;
t207 = pkin(9) * t390 + t225;
t315 = t410 * t389;
t376 = pkin(2) * t527;
t510 = -pkin(3) * t651 + t376;
t447 = t315 * pkin(9) + t510;
t576 = pkin(1) * t389;
t235 = -t447 - t576;
t136 = t398 * t207 + t393 * t235;
t122 = pkin(10) * t651 + t136;
t224 = t291 * t538 - t537 * t304;
t206 = -t390 * pkin(3) - t224;
t287 = t315 * t393 - t390 * t398;
t288 = t315 * t398 + t390 * t393;
t134 = t287 * pkin(4) - t288 * pkin(10) + t206;
t70 = t397 * t122 + t392 * t134;
t128 = -mrSges(6,1) * t195 + mrSges(6,2) * t196;
t204 = mrSges(5,1) * t303 - mrSges(5,3) * t273;
t650 = t204 - t128;
t499 = qJD(5) * t393;
t646 = t392 * t499 - t397 * t501 + t237;
t565 = mrSges(4,3) * t308;
t645 = mrSges(4,1) * t375 + mrSges(5,1) * t272 - mrSges(5,2) * t273 - t565;
t45 = -t393 * t105 + t140 * t398 - t180 * t501 - t205 * t502;
t641 = -t393 * t45 + t398 * t44;
t62 = mrSges(6,1) * t164 - mrSges(6,3) * t85;
t63 = -mrSges(6,2) * t164 + mrSges(6,3) * t86;
t640 = -t392 * t62 + t397 * t63;
t639 = t11 * t397 - t12 * t392;
t497 = m(5) - t668;
t562 = Ifges(3,4) * t394;
t636 = pkin(1) * (mrSges(3,1) * t394 + mrSges(3,2) * t399) - t394 * (Ifges(3,1) * t399 - t562) / 0.2e1;
t634 = mrSges(5,1) + t644;
t415 = mrSges(5,2) + t643;
t99 = -pkin(4) * t303 - t112;
t80 = -pkin(5) * t195 + t99;
t633 = -mrSges(7,1) * t80 + mrSges(7,3) * t17;
t632 = mrSges(7,2) * t80 - t16 * mrSges(7,3);
t631 = -t497 * pkin(9) - t568;
t508 = t410 * t390;
t263 = -t400 * t343 - t395 * t508;
t630 = t395 * t343 - t400 * t508;
t628 = t45 * mrSges(5,1) - t44 * mrSges(5,2) + Ifges(5,5) * t165 + Ifges(5,6) * t166 + Ifges(5,3) * t250;
t477 = t389 * t504;
t454 = pkin(8) * t477;
t267 = -qJD(1) * t454 + t484;
t268 = -pkin(8) * t328 - t394 * t451 + t370;
t626 = t268 * mrSges(3,1) + t110 * mrSges(4,1) - t267 * mrSges(3,2) - t111 * mrSges(4,2) + Ifges(3,5) * t328 + Ifges(4,5) * t252 + Ifges(3,6) * t327 + Ifges(4,6) * t251;
t625 = m(5) * pkin(9) - t671;
t36 = Ifges(6,4) * t85 + Ifges(6,2) * t86 + Ifges(6,6) * t164;
t618 = t36 / 0.2e1;
t617 = Ifges(6,1) * t611 + Ifges(6,4) * t610 + Ifges(6,5) * t599;
t556 = Ifges(7,4) * t127;
t59 = Ifges(7,2) * t458 + Ifges(7,6) * t264 + t556;
t616 = -t59 / 0.2e1;
t615 = t59 / 0.2e1;
t123 = Ifges(7,4) * t458;
t60 = Ifges(7,1) * t127 + Ifges(7,5) * t264 + t123;
t614 = -t60 / 0.2e1;
t613 = t60 / 0.2e1;
t559 = Ifges(6,4) * t196;
t97 = t195 * Ifges(6,2) + t270 * Ifges(6,6) + t559;
t609 = -t97 / 0.2e1;
t608 = t97 / 0.2e1;
t192 = Ifges(6,4) * t195;
t98 = t196 * Ifges(6,1) + t270 * Ifges(6,5) + t192;
t607 = -t98 / 0.2e1;
t606 = t98 / 0.2e1;
t595 = t195 / 0.2e1;
t593 = t196 / 0.2e1;
t587 = t270 / 0.2e1;
t586 = -t272 / 0.2e1;
t585 = -t273 / 0.2e1;
t584 = t273 / 0.2e1;
t582 = -t303 / 0.2e1;
t580 = t308 / 0.2e1;
t574 = pkin(5) * t196;
t231 = -t288 * t392 + t397 * t651;
t573 = pkin(5) * t231;
t566 = mrSges(4,3) * t307;
t564 = mrSges(6,3) * t195;
t563 = mrSges(6,3) * t196;
t561 = Ifges(5,4) * t393;
t560 = Ifges(5,4) * t398;
t558 = Ifges(6,4) * t392;
t557 = Ifges(6,4) * t397;
t549 = t308 * Ifges(4,4);
t548 = t375 * Ifges(3,5);
t547 = t375 * Ifges(3,6);
t40 = -pkin(4) * t250 - t45;
t544 = t393 * t40;
t532 = t307 * t398;
t528 = t389 * t395;
t526 = t389 * t400;
t520 = t394 * t400;
t518 = t395 * t399;
t239 = -t393 * t526 - t398 * t630;
t407 = t390 * t343;
t259 = -t395 * t410 - t400 * t407;
t514 = (-t239 * t386 - t259 * t387) * mrSges(7,1) + (-t239 * t387 + t259 * t386) * mrSges(7,2);
t243 = t263 * t398 + t393 * t528;
t262 = t395 * t407 - t400 * t410;
t152 = -t243 * t386 - t262 * t387;
t153 = t243 * t387 - t262 * t386;
t513 = t152 * mrSges(7,1) - t153 * mrSges(7,2);
t511 = (-t288 * t386 + t387 * t651) * mrSges(7,1) + (-t288 * t387 - t386 * t651) * mrSges(7,2);
t490 = m(4) + t497;
t489 = -m(3) * pkin(1) - mrSges(2,1);
t478 = t399 * t506;
t470 = t567 * t389;
t466 = t501 / 0.2e1;
t465 = -t499 / 0.2e1;
t459 = -t251 * mrSges(4,1) + t252 * mrSges(4,2);
t69 = -t122 * t392 + t397 * t134;
t135 = -t393 * t207 + t235 * t398;
t238 = t393 * t630 - t398 * t526;
t329 = pkin(2) * t525 - t470;
t457 = -t395 * t329 + t400 * t385;
t453 = mrSges(3,3) * t479;
t452 = mrSges(3,3) * t478;
t448 = t669 * pkin(2);
t444 = t263 * pkin(3) + t457;
t443 = mrSges(5,1) * t287 + mrSges(5,2) * t288;
t232 = t288 * t397 + t392 * t651;
t441 = mrSges(6,1) * t231 - mrSges(6,2) * t232;
t437 = Ifges(5,1) * t398 - t561;
t436 = Ifges(6,1) * t397 - t558;
t435 = Ifges(6,1) * t392 + t557;
t434 = -Ifges(5,2) * t393 + t560;
t433 = -Ifges(6,2) * t392 + t557;
t432 = Ifges(6,2) * t397 + t558;
t431 = Ifges(5,5) * t398 - Ifges(5,6) * t393;
t430 = Ifges(6,5) * t397 - Ifges(6,6) * t392;
t429 = Ifges(6,5) * t392 + Ifges(6,6) * t397;
t50 = pkin(5) * t287 - pkin(11) * t232 + t69;
t57 = pkin(11) * t231 + t70;
t22 = -t391 * t57 + t396 * t50;
t23 = t391 * t50 + t396 * t57;
t373 = t399 * t493;
t275 = t373 + (-qJD(2) * t469 + t503) * t389;
t276 = -t476 + (-t399 * t470 - t378) * qJD(2);
t183 = t275 * t537 - t538 * t276;
t138 = t231 * t396 - t232 * t391;
t139 = t231 * t391 + t232 * t396;
t160 = -t243 * t392 - t262 * t397;
t184 = t275 * t538 + t276 * t537;
t309 = qJD(2) * t315;
t310 = qJD(2) * t651;
t230 = pkin(2) * t477 + pkin(3) * t309 + pkin(9) * t310;
t75 = -t393 * t184 - t207 * t501 + t230 * t398 - t235 * t502;
t424 = t676 + t8;
t121 = -pkin(4) * t651 - t135;
t334 = -t390 * t518 - t520;
t422 = t99 * t439;
t74 = t398 * t184 - t207 * t502 + t393 * t230 + t235 * t501;
t66 = pkin(10) * t309 + t74;
t222 = qJD(4) * t288 - t310 * t393;
t223 = -qJD(4) * t287 - t310 * t398;
t92 = t222 * pkin(4) - t223 * pkin(10) + t183;
t20 = -t122 * t500 + t134 * t498 + t392 * t92 + t397 * t66;
t413 = t334 * pkin(2);
t67 = -pkin(4) * t309 - t75;
t21 = -qJD(5) * t70 - t392 * t66 + t397 * t92;
t404 = (-t392 * t56 - t397 * t55) * qJD(5) + t639;
t371 = Ifges(3,4) * t478;
t361 = Ifges(3,3) * t374;
t360 = Ifges(4,3) * t374;
t347 = -t376 - t576;
t339 = -pkin(8) * t529 + t379;
t336 = (-mrSges(3,1) * t399 + mrSges(3,2) * t394) * t389;
t335 = -t390 * t521 + t516;
t333 = -t390 * t520 - t518;
t330 = (t382 + t572) * t393;
t323 = t346 * t393;
t321 = t373 - t454;
t320 = t507 * qJD(1);
t319 = -pkin(8) * t479 + t372;
t318 = -mrSges(3,2) * t375 + t452;
t317 = mrSges(3,1) * t375 - t453;
t302 = Ifges(4,4) * t307;
t293 = Ifges(3,1) * t479 + t371 + t548;
t292 = t547 + (t399 * Ifges(3,2) + t562) * t506;
t280 = -t382 * t523 + t326;
t277 = -mrSges(4,2) * t375 - t566;
t266 = Ifges(5,4) * t272;
t244 = mrSges(4,1) * t307 + mrSges(4,2) * t308;
t242 = t263 * t393 - t398 * t528;
t234 = mrSges(4,1) * t374 - mrSges(4,3) * t252;
t233 = -mrSges(4,2) * t374 + mrSges(4,3) * t251;
t228 = t308 * Ifges(4,1) + t375 * Ifges(4,5) - t302;
t227 = -t307 * Ifges(4,2) + t375 * Ifges(4,6) + t549;
t203 = -mrSges(5,2) * t303 + mrSges(5,3) * t272;
t170 = t426 * t272;
t169 = t346 * t272;
t161 = t243 * t397 - t262 * t392;
t156 = t273 * Ifges(5,1) + t303 * Ifges(5,5) + t266;
t154 = t273 * Ifges(5,5) + t272 * Ifges(5,6) + t303 * Ifges(5,3);
t147 = mrSges(6,1) * t270 - t563;
t146 = -mrSges(6,2) * t270 + t564;
t120 = qJD(5) * t231 + t223 * t397 + t309 * t392;
t119 = -qJD(5) * t232 - t223 * t392 + t309 * t397;
t118 = -mrSges(5,2) * t250 + mrSges(5,3) * t166;
t95 = mrSges(7,1) * t264 - mrSges(7,3) * t127;
t94 = -mrSges(7,2) * t264 + mrSges(7,3) * t458;
t91 = -mrSges(5,1) * t166 + mrSges(5,2) * t165;
t88 = t121 - t573;
t78 = t165 * Ifges(5,4) + t166 * Ifges(5,2) + t250 * Ifges(5,6);
t68 = -mrSges(7,1) * t458 + mrSges(7,2) * t127;
t47 = -qJD(6) * t139 + t119 * t396 - t120 * t391;
t46 = qJD(6) * t138 + t119 * t391 + t120 * t396;
t43 = -pkin(5) * t119 + t67;
t26 = -mrSges(7,2) * t159 + mrSges(7,3) * t34;
t25 = mrSges(7,1) * t159 - mrSges(7,3) * t33;
t24 = -pkin(5) * t86 + t40;
t19 = t396 * t48 - t546;
t18 = -t391 * t48 - t542;
t15 = pkin(11) * t119 + t20;
t14 = pkin(5) * t222 - pkin(11) * t120 + t21;
t13 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t5 = -qJD(6) * t23 + t14 * t396 - t15 * t391;
t4 = qJD(6) * t22 + t14 * t391 + t15 * t396;
t1 = [t223 * t687 + t272 * (Ifges(5,4) * t223 + Ifges(5,6) * t309) / 0.2e1 + (t138 * t2 - t139 * t3 - t16 * t46 + t17 * t47) * mrSges(7,3) + (Ifges(6,4) * t232 + Ifges(6,2) * t231) * t610 + (Ifges(6,4) * t120 + Ifges(6,2) * t119) * t595 + (Ifges(7,4) * t46 + Ifges(7,2) * t47) * t604 + ((mrSges(3,1) * t327 - mrSges(3,2) * t328 + (m(3) * t576 - t336) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t267 + Ifges(3,4) * t328 + Ifges(3,2) * t327 + Ifges(3,6) * t374) * t399 + (-mrSges(3,3) * t268 + Ifges(3,1) * t328 + Ifges(3,4) * t327 + Ifges(3,5) * t374) * t394 + ((t548 / 0.2e1 + t293 / 0.2e1 - t319 * mrSges(3,3)) * t399 + (-t547 / 0.2e1 - t292 / 0.2e1 - t320 * mrSges(3,3) + (m(4) * t331 + t244) * pkin(2)) * t394 + (t399 * (Ifges(3,4) * t399 - Ifges(3,2) * t394) / 0.2e1 - t636) * t506) * qJD(2) + (g(1) * t400 + g(2) * t395) * (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3))) * t389 + (-m(4) * t457 - t263 * mrSges(4,1) - t335 * mrSges(3,1) - t334 * mrSges(3,2) + mrSges(2,2) * t395 - m(7) * (t243 * t384 + t444) - t153 * mrSges(7,1) - t152 * mrSges(7,2) - m(6) * (pkin(4) * t243 + t444) - t161 * mrSges(6,1) - t160 * mrSges(6,2) - m(5) * t444 - t243 * mrSges(5,1) + t489 * t400 + (m(7) * t572 - t631) * t262 + t415 * t242) * g(2) + t303 * (Ifges(5,5) * t223 + Ifges(5,3) * t309) / 0.2e1 + (Ifges(7,4) * t139 + Ifges(7,2) * t138) * t619 + (t190 * t310 - t191 * t309) * mrSges(4,3) + t375 * (-Ifges(4,5) * t310 - Ifges(4,6) * t309) / 0.2e1 + t331 * (mrSges(4,1) * t309 - mrSges(4,2) * t310) - t307 * (-Ifges(4,4) * t310 - Ifges(4,2) * t309) / 0.2e1 + (-Ifges(4,1) * t310 - Ifges(4,4) * t309) * t580 + (Ifges(5,1) * t223 + Ifges(5,5) * t309) * t584 - t309 * t661 + (-t333 * mrSges(3,1) + t669 * mrSges(3,2) + (t329 * t490 + mrSges(2,2)) * t400 + (t385 * t490 - t489) * t395 + t634 * t239 + (t392 * t659 - t438 - t541 + t631) * t259 + t415 * t238 - (pkin(3) * t497 + mrSges(4,1)) * t630) * g(1) + t309 * t662 + (t11 * t231 + t119 * t56 - t12 * t232 - t120 * t55) * mrSges(6,3) + m(5) * (t104 * t206 + t112 * t75 + t113 * t74 + t135 * t45 + t136 * t44 + t179 * t183) + m(6) * (t11 * t70 + t12 * t69 + t121 * t40 + t20 * t56 + t21 * t55 + t67 * t99) + m(7) * (t16 * t5 + t17 * t4 + t2 * t23 + t22 * t3 + t24 * t88 + t43 * t80) + (t290 * mrSges(4,1) - t111 * mrSges(4,3) - Ifges(4,4) * t252 - Ifges(4,2) * t251 - Ifges(4,6) * t374 + t628) * t651 + (Ifges(7,5) * t139 + Ifges(7,6) * t138) * t600 + (Ifges(7,5) * t46 + Ifges(7,6) * t47) * t589 + (mrSges(4,2) * t290 - mrSges(4,3) * t110 + Ifges(4,1) * t252 + Ifges(4,4) * t251 + Ifges(4,5) * t374) * t315 + t507 * (-mrSges(3,2) * t374 + mrSges(3,3) * t327) + m(3) * (t267 * t507 + t268 * t339 - t319 * t322 + t320 * t321) + (Ifges(6,5) * t232 + Ifges(6,6) * t231) * t599 + (Ifges(6,5) * t120 + Ifges(6,6) * t119) * t587 + (Ifges(7,1) * t139 + Ifges(7,4) * t138) * t620 + (Ifges(7,1) * t46 + Ifges(7,4) * t47) * t602 + (Ifges(6,1) * t232 + Ifges(6,4) * t231) * t611 + (Ifges(6,1) * t120 + Ifges(6,4) * t119) * t593 - t645 * t183 + t347 * t459 + (-t685 / 0.2e1 - t686 / 0.2e1 + Ifges(6,6) * t595 - Ifges(5,4) * t584 + Ifges(6,3) * t587 + Ifges(6,5) * t593 + t675 + t689) * t222 + (t672 / 0.2e1 - t44 * mrSges(5,3) - t78 / 0.2e1 + Ifges(7,6) * t619 + Ifges(7,5) * t620 + Ifges(6,6) * t610 + Ifges(6,5) * t611 - Ifges(5,2) * t597 - Ifges(5,4) * t598 + Ifges(6,3) * t599 + Ifges(7,3) * t600 - Ifges(5,6) * t591 + t676 + t690) * t287 + (-mrSges(5,3) * t45 + 0.2e1 * t612) * t288 - t223 * t553 - t40 * t441 + t104 * t443 + t22 * t25 + t23 * t26 + Ifges(2,3) * qJDD(1) + t339 * (mrSges(3,1) * t374 - mrSges(3,3) * t328) + t321 * t318 - t322 * t317 - t310 * t228 / 0.2e1 + t309 * t154 / 0.2e1 - t309 * t227 / 0.2e1 + t184 * t277 + t225 * t233 + t224 * t234 + t223 * t156 / 0.2e1 + t74 * t203 + t75 * t204 + t206 * t91 + t139 * t621 + t138 * t623 + t46 * t613 + t47 * t615 + t232 * t617 + t231 * t618 + t120 * t606 + t119 * t608 + t43 * t68 + t69 * t62 + t70 * t63 + t80 * (-mrSges(7,1) * t47 + mrSges(7,2) * t46) + t88 * t13 + t4 * t94 + t5 * t95 + (t360 / 0.2e1 + t361 / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t374 + t626) * t390 + m(4) * (t110 * t224 + t111 * t225 - t183 * t190 + t184 * t191 + t290 * t347) + t99 * (-mrSges(6,1) * t119 + mrSges(6,2) * t120) + t121 * t42 + t67 * t128 + t135 * t117 + t136 * t118 + t24 * (-mrSges(7,1) * t138 + mrSges(7,2) * t139) + t20 * t146 + t21 * t147; (Ifges(6,1) * t237 + Ifges(6,4) * t236) * t594 + t664 * t94 + (t16 * t665 + t167 * t3 + t168 * t2 + t17 * t664 + t24 * t330 + t653 * t80) * m(7) + t665 * t95 + (t11 * t281 - t114 * t99 + t12 * t280 + t656 * t55 + t657 * t56) * m(6) + (t104 * t383 - t112 * t132 - t113 * t133 - t179 * t212) * m(5) + (t658 * t393 + (t501 * t99 + t544) * m(6) + ((-t112 * t398 - t113 * t393) * qJD(4) + t641) * m(5) + t398 * t118) * t382 + t645 * t212 + (mrSges(6,1) * t647 - mrSges(6,2) * t646) * t99 + (-t11 * t524 - t12 * t522 + t55 * t646 - t56 * t647) * mrSges(6,3) + (-mrSges(7,2) * t642 + mrSges(7,3) * t648) * t17 + (-mrSges(7,1) * t648 + mrSges(7,2) * t649) * t80 + (mrSges(7,1) * t642 - mrSges(7,3) * t649) * t16 + t638 * t179 * (mrSges(5,1) * t393 + mrSges(5,2) * t398) + t636 * qJD(1) ^ 2 * t389 ^ 2 + (t532 / 0.2e1 + t466) * t156 + (-t112 * t532 + t641) * mrSges(5,3) - (-Ifges(4,1) * t307 + t154 - t549) * t308 / 0.2e1 - t375 * (-Ifges(4,5) * t307 - Ifges(4,6) * t308) / 0.2e1 - t331 * (mrSges(4,1) * t308 - mrSges(4,2) * t307) + (Ifges(5,3) * t308 - t307 * t431) * t582 + (Ifges(5,5) * t308 - t307 * t437) * t585 + (Ifges(5,6) * t308 - t307 * t434) * t586 + ((t110 * t538 + t111 * t537) * pkin(2) + t190 * t212 - t191 * t213 - t331 * t455) * m(4) + (-t475 - t133) * t203 + t626 + t653 * t68 + t656 * t147 + t657 * t146 + t308 * t661 + (-m(4) * t448 - mrSges(3,1) * t669 - mrSges(3,2) * t333 - t497 * (t259 * pkin(3) + t448) + t629 * t259 + t625 * t630) * g(2) - t308 * t662 + t398 * t667 + (-t474 - t132) * t204 + t439 * t544 - t672 * t398 / 0.2e1 + (Ifges(7,1) * t214 + Ifges(7,4) * t215) * t602 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t603 - t398 * t666 + (Ifges(6,4) * t237 + Ifges(6,2) * t236) * t596 + (Ifges(6,5) * t237 + Ifges(6,6) * t236) * t588 + (-m(4) * t413 - mrSges(3,1) * t334 + mrSges(3,2) * t335 - t497 * (t262 * pkin(3) + t413) + t629 * t262 - t625 * t263) * g(1) + (t465 * t97 + t466 * t98) * t397 - (t375 * (Ifges(3,5) * t399 - Ifges(3,6) * t394) + (-Ifges(3,2) * t479 + t293 + t371) * t399) * t506 / 0.2e1 + (Ifges(7,4) * t214 + Ifges(7,2) * t215) * t604 + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t605 + t392 * t98 * t465 + (t272 * t434 + t273 * t437 + t303 * t431) * qJD(4) / 0.2e1 + t191 * t565 + t360 + t361 + t292 * t479 / 0.2e1 + t670 * t128 + (-Ifges(4,2) * t308 + t228 - t302) * t307 / 0.2e1 + t2 * (mrSges(7,2) * t398 - mrSges(7,3) * t323) + t24 * (mrSges(7,1) * t323 - mrSges(7,2) * t324) + (-Ifges(7,1) * t324 - Ifges(7,4) * t323 - Ifges(7,5) * t398) * t620 + (-Ifges(7,4) * t324 - Ifges(7,2) * t323 - Ifges(7,6) * t398) * t619 + (-Ifges(7,5) * t324 - Ifges(7,6) * t323 - Ifges(7,3) * t398) * t600 + t3 * (-mrSges(7,1) * t398 + mrSges(7,3) * t324) - t244 * t455 + t233 * t480 + t234 * t481 - (t677 + t688) * t533 + (-m(4) * t376 - m(5) * t447 + t671 * t315 + t668 * t510 - t629 * t651 + t336) * g(3) + t689 * t502 - t104 * t442 + t398 * t78 / 0.2e1 + t383 * t91 + t330 * t13 + (t453 + t317) * t320 + t280 * t62 + t281 * t63 - t213 * t277 + (t452 - t318) * t319 + t167 * t25 + t168 * t26 - t324 * t621 - t323 * t623 + t144 * t614 + t215 * t615 + t143 * t616 + t522 * t617 + t237 * t607 + (-Ifges(6,6) * t398 + t393 * t433) * t610 + (-Ifges(6,5) * t398 + t393 * t436) * t611 + t393 * t612 + t214 * t613 + (-t432 * t499 + (Ifges(6,6) * t393 + t398 * t433) * qJD(4)) * t595 + (Ifges(5,2) * t398 + t561) * t597 + (Ifges(5,1) * t393 + t560) * t598 + (-Ifges(6,3) * t398 + t393 * t430) * t599 + (-t429 * t499 + (Ifges(6,3) * t393 + t398 * t430) * qJD(4)) * t587 + (Ifges(5,5) * t393 + Ifges(5,6) * t398) * t591 + (-t435 * t499 + (Ifges(6,5) * t393 + t398 * t436) * qJD(4)) * t593 - t190 * t566 - t501 * t553 + t693 * t609 + t227 * t580 + (Ifges(7,5) * t214 + Ifges(7,6) * t215) * t589 + (Ifges(7,5) * t144 + Ifges(7,6) * t143) * t590 - t36 * t524 / 0.2e1; t645 * t308 + t459 + t649 * t94 + t648 * t95 - t323 * t25 - t324 * t26 + t307 * t277 + (t118 + (-t392 * t146 - t397 * t147) * qJD(5) + t638 * (t68 - t650) + t640) * t393 - t236 * t147 - t237 * t146 + (t307 * t203 - t13 + (t146 * t397 - t147 * t392 + t203) * qJD(4) - t658) * t398 + (-t390 * g(3) + (-g(1) * t395 + g(2) * t400) * t389) * t490 + (t16 * t648 + t17 * t649 - t2 * t324 - t24 * t398 - t3 * t323 + t642 * t80) * m(7) + (-t236 * t55 - t237 * t56 + t99 * t533 + (-t40 + (-t392 * t55 + t397 * t56) * qJD(4)) * t398 + (qJD(4) * t99 + t404) * t393) * m(6) + (-t179 * t308 + t393 * t44 + t398 * t45 + t638 * (-t112 * t393 + t113 * t398)) * m(5) + (t190 * t308 + t191 * t307 + t290) * m(4); (-t550 + t663) * t585 + t650 * t113 + (t287 * t644 + t288 * t643 + t443) * g(3) + (t242 * t634 + t243 * t415) * g(1) + (-t238 * t634 + t239 * t415) * g(2) + (-pkin(4) * t40 - t113 * t99 - t55 * t81 - t56 * t82) * m(6) + t628 + t652 * t68 + t654 * t95 + (t16 * t654 + t17 * t655 + t2 * t300 - t24 * t384 + t299 * t3 + t652 * t80) * m(7) + t655 * t94 + (-t16 * t170 + t169 * t17) * mrSges(7,3) + (Ifges(5,1) * t585 + Ifges(5,5) * t582 + t430 * t588 + t433 * t596 + t436 * t594 - t422 + t553 - t687) * t272 + (t266 + t156) * t586 - pkin(4) * t42 + (m(6) * t404 - t146 * t500 - t147 * t498 + t640) * pkin(10) + qJD(5) * t422 + (t195 * t433 + t196 * t436 + t270 * t430) * qJD(5) / 0.2e1 + (-Ifges(7,4) * t170 - Ifges(7,2) * t169) * t605 + (-Ifges(7,5) * t170 - Ifges(7,6) * t169) * t590 + (-Ifges(7,1) * t170 - Ifges(7,4) * t169) * t603 - t80 * (mrSges(7,1) * t169 - mrSges(7,2) * t170) - (Ifges(7,4) * t602 + Ifges(7,2) * t604 + Ifges(7,6) * t589 + t615 + t633) * t286 + (-(Ifges(7,1) * t602 + Ifges(7,4) * t604 + Ifges(7,5) * t589 + t613 + t632) * t637 - t2 * mrSges(7,3) + t24 * mrSges(7,1) - 0.2e1 * t623) * t426 + (mrSges(7,2) * t24 - mrSges(7,3) * t3 + 0.2e1 * t621) * t346 + t40 * t440 - t384 * t13 + t299 * t25 + t300 * t26 - t112 * t203 - t170 * t614 - t169 * t616 + t392 * t617 + t397 * t618 + t498 * t606 + t534 * t607 + t535 * t608 + t500 * t609 + t432 * t610 + t435 * t611 + t429 * t599 + t155 * t584 + (-Ifges(5,2) * t586 - Ifges(5,6) * t582 - t675 + t677) * t273 + (-t678 * t56 + (-t498 + t534) * t55 + t639) * mrSges(6,3) - t82 * t146 - t81 * t147; (-(-t239 * t397 + t259 * t392) * mrSges(6,2) - t514 + t659 * (-t239 * t392 - t259 * t397)) * g(2) + (mrSges(6,2) * t161 + t160 * t659 - t513) * g(1) + (-Ifges(6,2) * t196 + t192 + t98) * t596 + t424 + ((-t391 * t95 + t396 * t94) * qJD(6) + t25 * t396 + t26 * t391) * pkin(5) + t690 + t35 - t68 * t574 - m(7) * (t16 * t18 + t17 * t19 + t574 * t80) + (t564 - t146) * t55 + (Ifges(7,1) * t603 + Ifges(7,4) * t605 + Ifges(7,5) * t590 + t614 - t632) * t458 - (Ifges(7,4) * t603 + Ifges(7,2) * t605 + Ifges(7,6) * t590 + t616 - t633) * t127 + (-m(7) * t573 - t441 - t511) * g(3) - t99 * (mrSges(6,1) * t196 + mrSges(6,2) * t195) + (t2 * t391 + t3 * t396 + (-t16 * t391 + t17 * t396) * qJD(6)) * t622 + (Ifges(6,1) * t195 - t559) * t594 + (Ifges(6,5) * t195 - Ifges(6,6) * t196) * t588 + t97 * t593 + (t563 + t147) * t56 - t19 * t94 - t18 * t95; -t80 * (mrSges(7,1) * t127 + mrSges(7,2) * t458) + (Ifges(7,5) * t458 - Ifges(7,6) * t127) * t590 + (Ifges(7,1) * t458 - t556) * t603 + t59 * t602 - t16 * t94 + t17 * t95 - g(1) * t513 - g(2) * t514 - g(3) * t511 + (t127 * t17 + t16 * t458) * mrSges(7,3) + t424 + (-Ifges(7,2) * t127 + t123 + t60) * t605;];
tau  = t1;
