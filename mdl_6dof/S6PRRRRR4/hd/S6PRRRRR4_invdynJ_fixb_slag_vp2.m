% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:55:03
% EndTime: 2019-03-09 00:56:11
% DurationCPUTime: 40.84s
% Computational Cost: add. (23222->985), mult. (56953->1399), div. (0->0), fcn. (48084->18), ass. (0->455)
t389 = sin(pkin(6));
t396 = sin(qJ(2));
t543 = t389 * t396;
t506 = qJD(1) * t543;
t388 = sin(pkin(7));
t530 = qJD(2) * t388;
t344 = pkin(9) * t530 + t506;
t395 = sin(qJ(3));
t400 = cos(qJ(3));
t596 = cos(qJ(2));
t510 = t389 * t596;
t477 = qJD(1) * t510;
t356 = qJD(2) * pkin(2) + t477;
t390 = cos(pkin(7));
t391 = cos(pkin(6));
t531 = qJD(1) * t391;
t507 = t388 * t531;
t701 = t356 * t390 + t507;
t226 = -t395 * t344 + t400 * t701;
t438 = (pkin(3) * t395 - pkin(10) * t400) * t388;
t323 = qJD(2) * t438;
t394 = sin(qJ(4));
t399 = cos(qJ(4));
t167 = -t226 * t394 + t399 * t323;
t401 = -pkin(11) - pkin(10);
t511 = qJD(4) * t401;
t536 = t399 * t400;
t705 = -(pkin(4) * t395 - pkin(11) * t536) * t530 - t167 + t399 * t511;
t168 = t399 * t226 + t394 * t323;
t504 = t400 * t530;
t478 = t394 * t504;
t704 = -pkin(11) * t478 - t394 * t511 + t168;
t527 = qJD(4) * t394;
t703 = -t527 + t478;
t694 = mrSges(6,2) - mrSges(7,3);
t542 = t390 * t395;
t544 = t388 * t400;
t343 = pkin(2) * t542 + pkin(9) * t544;
t311 = pkin(10) * t390 + t343;
t461 = -pkin(3) * t400 - pkin(10) * t395;
t312 = (-pkin(2) + t461) * t388;
t229 = t399 * t311 + t394 * t312;
t324 = qJD(3) * t438;
t546 = t388 * t395;
t377 = pkin(9) * t546;
t541 = t390 * t400;
t342 = pkin(2) * t541 - t377;
t325 = t342 * qJD(3);
t147 = -qJD(4) * t229 + t399 * t324 - t325 * t394;
t339 = t390 * t399 - t394 * t546;
t528 = qJD(3) * t400;
t500 = t388 * t528;
t268 = qJD(4) * t339 + t399 * t500;
t529 = qJD(3) * t395;
t501 = t388 * t529;
t101 = pkin(4) * t501 - pkin(11) * t268 + t147;
t526 = qJD(4) * t399;
t146 = -t311 * t527 + t312 * t526 + t394 * t324 + t399 * t325;
t545 = t388 * t399;
t340 = t390 * t394 + t395 * t545;
t269 = -qJD(4) * t340 - t394 * t500;
t112 = pkin(11) * t269 + t146;
t228 = -t311 * t394 + t399 * t312;
t169 = -pkin(4) * t544 - pkin(11) * t340 + t228;
t180 = pkin(11) * t339 + t229;
t508 = t596 * t400;
t539 = t395 * t396;
t424 = -t390 * t539 + t508;
t300 = t424 * t389;
t284 = qJD(1) * t300;
t480 = t388 * t506;
t239 = -t284 * t394 + t399 * t480;
t240 = t284 * t399 + t394 * t480;
t393 = sin(qJ(5));
t398 = cos(qJ(5));
t524 = qJD(5) * t398;
t525 = qJD(5) * t393;
t672 = t169 * t524 - t180 * t525 + (t112 - t240) * t398 + (t101 - t239) * t393;
t366 = t401 * t394;
t367 = t401 * t399;
t439 = t398 * t366 + t367 * t393;
t671 = qJD(5) * t439 + t393 * t705 - t704 * t398;
t227 = t400 * t344 + t356 * t542 + t395 * t507;
t665 = -pkin(4) * t703 - t227;
t503 = qJD(2) * t543;
t476 = qJD(1) * t503;
t331 = qJDD(1) * t510 - t476;
t308 = qJDD(2) * pkin(2) + t331;
t521 = qJDD(1) * t391;
t496 = t388 * t521;
t369 = qJD(2) * t477;
t332 = qJDD(1) * t543 + t369;
t520 = qJDD(2) * t388;
t687 = pkin(9) * t520 + qJD(3) * t701 + t332;
t123 = t400 * (t308 * t390 + t496) - t344 * t528 - t687 * t395;
t375 = qJDD(2) * t390 + qJDD(3);
t109 = -pkin(3) * t375 - t123;
t376 = qJD(2) * t390 + qJD(3);
t505 = t395 * t530;
t298 = t376 * t394 + t399 * t505;
t522 = qJD(2) * qJD(3);
t330 = (qJDD(2) * t395 + t400 * t522) * t388;
t204 = -qJD(4) * t298 - t330 * t394 + t375 * t399;
t74 = -pkin(4) * t204 + t109;
t297 = t376 * t399 - t394 * t505;
t203 = qJD(4) * t297 + t330 * t399 + t375 * t394;
t487 = t398 * t297 - t298 * t393;
t84 = qJD(5) * t487 + t203 * t398 + t204 * t393;
t441 = t297 * t393 + t398 * t298;
t85 = -qJD(5) * t441 - t203 * t393 + t204 * t398;
t20 = -pkin(5) * t85 - pkin(12) * t84 + t74;
t392 = sin(qJ(6));
t397 = cos(qJ(6));
t362 = qJD(4) - t504;
t357 = qJD(5) + t362;
t206 = pkin(10) * t376 + t227;
t372 = t390 * t531;
t243 = t372 + (qJD(2) * t461 - t356) * t388;
t130 = t206 * t399 + t243 * t394;
t106 = pkin(11) * t297 + t130;
t537 = t398 * t106;
t129 = -t206 * t394 + t399 * t243;
t105 = -pkin(11) * t298 + t129;
t94 = pkin(4) * t362 + t105;
t56 = t393 * t94 + t537;
t48 = pkin(12) * t357 + t56;
t205 = -pkin(3) * t376 - t226;
t157 = -pkin(4) * t297 + t205;
t75 = -pkin(5) * t487 - pkin(12) * t441 + t157;
t21 = -t392 * t48 + t397 * t75;
t329 = -qJD(2) * t501 + t400 * t520;
t315 = qJDD(4) - t329;
t307 = qJDD(5) + t315;
t122 = t308 * t542 - t344 * t529 + t395 * t496 + t400 * t687;
t108 = pkin(10) * t375 + t122;
t258 = -t308 * t388 + t390 * t521;
t166 = -pkin(3) * t329 - pkin(10) * t330 + t258;
t45 = -qJD(4) * t130 - t108 * t394 + t399 * t166;
t34 = pkin(4) * t315 - pkin(11) * t203 + t45;
t44 = t399 * t108 + t394 * t166 - t206 * t527 + t243 * t526;
t40 = pkin(11) * t204 + t44;
t8 = -t106 * t525 + t393 * t34 + t398 * t40 + t94 * t524;
t5 = pkin(12) * t307 + t8;
t2 = qJD(6) * t21 + t20 * t392 + t397 * t5;
t22 = t392 * t75 + t397 * t48;
t3 = -qJD(6) * t22 + t20 * t397 - t392 * t5;
t450 = t21 * t397 + t22 * t392;
t702 = -qJD(6) * t450 + t2 * t397 - t3 * t392;
t458 = -mrSges(7,1) * t397 + mrSges(7,2) * t392;
t688 = -mrSges(6,1) + t458;
t700 = m(7) * pkin(5) - t688;
t649 = -m(7) * pkin(12) + t694;
t699 = -pkin(12) * t501 - t672;
t698 = -pkin(12) * t505 + t671;
t440 = t398 * t339 - t340 * t393;
t125 = qJD(5) * t440 + t268 * t398 + t269 * t393;
t248 = t339 * t393 + t340 * t398;
t126 = qJD(5) * t248 + t268 * t393 - t398 * t269;
t326 = t343 * qJD(3);
t224 = -pkin(4) * t269 + t326;
t509 = t596 * t395;
t538 = t396 * t400;
t423 = t390 * t538 + t509;
t299 = t423 * t389;
t283 = qJD(1) * t299;
t697 = pkin(5) * t126 - pkin(12) * t125 + t224 - t283;
t347 = t393 * t394 - t398 * t399;
t655 = qJD(4) + qJD(5);
t264 = t655 * t347;
t348 = t393 * t399 + t394 * t398;
t265 = t655 * t348;
t276 = t348 * t504;
t277 = t347 * t504;
t696 = t665 + (t264 - t277) * pkin(12) + (t265 - t276) * pkin(5);
t681 = t22 * mrSges(7,2);
t682 = t21 * mrSges(7,1);
t695 = -t157 * mrSges(6,1) + t56 * mrSges(6,3) + t681 - t682;
t279 = t366 * t393 - t367 * t398;
t670 = -qJD(5) * t279 + t704 * t393 + t398 * t705;
t172 = t357 * t397 - t392 * t441;
t57 = qJD(6) * t172 + t307 * t392 + t397 * t84;
t173 = t357 * t392 + t397 * t441;
t58 = -qJD(6) * t173 + t307 * t397 - t392 * t84;
t19 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t9 = -qJD(5) * t56 + t34 * t398 - t393 * t40;
t6 = -pkin(5) * t307 - t9;
t668 = m(7) * t6 + t19;
t83 = qJDD(6) - t85;
t16 = t57 * Ifges(7,1) + t58 * Ifges(7,4) + t83 * Ifges(7,5);
t216 = qJD(6) - t487;
t451 = Ifges(7,5) * t397 - Ifges(7,6) * t392;
t429 = t216 * t451;
t579 = Ifges(7,4) * t392;
t455 = Ifges(7,1) * t397 - t579;
t431 = t173 * t455;
t578 = Ifges(7,4) * t397;
t453 = -Ifges(7,2) * t392 + t578;
t432 = t172 * t453;
t457 = mrSges(7,1) * t392 + mrSges(7,2) * t397;
t564 = t106 * t393;
t55 = t398 * t94 - t564;
t47 = -pkin(5) * t357 - t55;
t436 = t47 * t457;
t523 = qJD(6) * t392;
t493 = -t523 / 0.2e1;
t170 = Ifges(7,4) * t172;
t70 = Ifges(7,1) * t173 + Ifges(7,5) * t216 + t170;
t569 = t397 * t70;
t514 = t569 / 0.2e1;
t518 = Ifges(6,5) * t84 + Ifges(6,6) * t85 + Ifges(6,3) * t307;
t597 = t392 / 0.2e1;
t625 = t83 / 0.2e1;
t629 = t58 / 0.2e1;
t630 = t57 / 0.2e1;
t15 = t57 * Ifges(7,4) + t58 * Ifges(7,2) + t83 * Ifges(7,6);
t634 = t15 / 0.2e1;
t580 = Ifges(7,4) * t173;
t69 = Ifges(7,2) * t172 + Ifges(7,6) * t216 + t580;
t693 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + t16 * t597 + t397 * t634 + t6 * t458 + (Ifges(7,1) * t392 + t578) * t630 + (Ifges(7,2) * t397 + t579) * t629 + t518 + t69 * t493 + (Ifges(7,5) * t392 + Ifges(7,6) * t397) * t625 + (t436 + t514) * qJD(6) + (t432 + t431 + t429) * qJD(6) / 0.2e1 + t702 * mrSges(7,3);
t422 = t390 * t509 + t538;
t261 = t389 * t422 + t391 * t546;
t336 = -t388 * t510 + t391 * t390;
t211 = -t261 * t394 + t336 * t399;
t566 = cos(pkin(13));
t463 = t566 * t596;
t565 = sin(pkin(13));
t488 = t565 * t396;
t338 = -t391 * t488 + t463;
t462 = t565 * t596;
t489 = t566 * t396;
t412 = t391 * t462 + t489;
t490 = t389 * t565;
t465 = t388 * t490;
t210 = t338 * t400 + (-t390 * t412 + t465) * t395;
t263 = t388 * t412 + t390 * t490;
t692 = -t210 * t394 + t263 * t399;
t337 = t391 * t489 + t462;
t411 = -t391 * t463 + t488;
t409 = t411 * t395;
t491 = t389 * t566;
t466 = t388 * t491;
t208 = t337 * t400 - t390 * t409 - t395 * t466;
t262 = t388 * t411 - t390 * t491;
t691 = -t208 * t394 + t262 * t399;
t179 = mrSges(6,1) * t357 - mrSges(6,3) * t441;
t93 = -mrSges(7,1) * t172 + mrSges(7,2) * t173;
t567 = t179 - t93;
t690 = -m(6) * t55 + m(7) * t47 - t567;
t574 = Ifges(6,6) * t357;
t575 = Ifges(6,2) * t487;
t581 = Ifges(6,4) * t441;
t117 = t574 + t575 + t581;
t598 = t357 / 0.2e1;
t604 = t441 / 0.2e1;
t605 = t487 / 0.2e1;
t606 = t216 / 0.2e1;
t613 = t173 / 0.2e1;
t615 = t172 / 0.2e1;
t572 = Ifges(7,3) * t216;
t573 = Ifges(7,6) * t172;
t576 = Ifges(7,5) * t173;
t68 = t572 + t573 + t576;
t689 = -Ifges(6,4) * t604 + Ifges(7,5) * t613 - Ifges(6,2) * t605 - Ifges(6,6) * t598 + Ifges(7,6) * t615 + Ifges(7,3) * t606 - t117 / 0.2e1 + t68 / 0.2e1 - t695;
t686 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - t457;
t600 = t307 / 0.2e1;
t623 = t85 / 0.2e1;
t624 = t84 / 0.2e1;
t631 = Ifges(6,1) * t624 + Ifges(6,4) * t623 + Ifges(6,5) * t600;
t684 = -m(7) - m(6);
t610 = t203 / 0.2e1;
t609 = t204 / 0.2e1;
t599 = t315 / 0.2e1;
t683 = m(6) * t157;
t384 = pkin(4) * t399 + pkin(3);
t255 = pkin(5) * t347 - pkin(12) * t348 - t384;
t176 = t255 * t392 + t279 * t397;
t679 = -qJD(6) * t176 - t392 * t698 + t397 * t696;
t175 = t255 * t397 - t279 * t392;
t678 = qJD(6) * t175 + t392 * t696 + t397 * t698;
t677 = t122 * mrSges(4,2);
t676 = t123 * mrSges(4,1);
t310 = t377 + (-pkin(2) * t400 - pkin(3)) * t390;
t250 = -pkin(4) * t339 + t310;
t115 = -pkin(5) * t440 - pkin(12) * t248 + t250;
t667 = t393 * t169 + t398 * t180;
t89 = -pkin(12) * t544 + t667;
t51 = t115 * t397 - t392 * t89;
t674 = qJD(6) * t51 + t392 * t697 - t397 * t699;
t52 = t115 * t392 + t397 * t89;
t673 = -qJD(6) * t52 + t392 * t699 + t397 * t697;
t669 = pkin(5) * t505 - t670;
t666 = t298 * Ifges(5,5) + Ifges(6,5) * t441 + t297 * Ifges(5,6) + Ifges(6,6) * t487 + t362 * Ifges(5,3) + t357 * Ifges(6,3);
t486 = mrSges(4,3) * t505;
t664 = -mrSges(4,1) * t376 - mrSges(5,1) * t297 + mrSges(5,2) * t298 + t486;
t663 = -t239 + t147;
t662 = -t240 + t146;
t244 = t277 * t392 + t397 * t505;
t551 = t348 * t397;
t434 = qJD(6) * t551 - t392 * t264;
t661 = t244 + t434;
t245 = -t277 * t397 + t392 * t505;
t433 = t397 * t264 + t348 * t523;
t660 = t245 + t433;
t212 = t261 * t399 + t336 * t394;
t120 = t211 * t393 + t212 * t398;
t657 = t390 * t508 - t539;
t260 = -t389 * t657 - t391 * t544;
t256 = t260 * t397;
t86 = -t120 * t392 + t256;
t659 = -t284 + t325;
t215 = Ifges(6,4) * t487;
t577 = Ifges(6,5) * t357;
t586 = Ifges(6,1) * t441;
t118 = t215 + t577 + t586;
t618 = -t118 / 0.2e1;
t651 = -t157 * mrSges(6,2) + t55 * mrSges(6,3);
t658 = t450 * mrSges(7,3) + t618 - t577 / 0.2e1 - t432 / 0.2e1 - t431 / 0.2e1 - t429 / 0.2e1 + t69 * t597 - t569 / 0.2e1 - t436 + t651 - t215 / 0.2e1;
t656 = -t394 * t45 + t399 * t44;
t287 = -t356 * t388 + t372;
t654 = (t376 * (Ifges(4,5) * t400 - Ifges(4,6) * t395) / 0.2e1 + t287 * (mrSges(4,1) * t395 + mrSges(4,2) * t400)) * t388;
t66 = mrSges(6,1) * t307 - mrSges(6,3) * t84;
t653 = m(6) * t9 + t66;
t28 = -qJD(5) * t667 + t101 * t398 - t112 * t393;
t387 = qJ(4) + qJ(5);
t385 = sin(t387);
t386 = cos(t387);
t460 = -mrSges(5,1) * t399 + mrSges(5,2) * t394;
t648 = m(5) * pkin(3) - t649 * t385 + t386 * t700 + mrSges(4,1) - t460;
t646 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t645 = -m(5) * t205 - t664;
t128 = pkin(5) * t441 - pkin(12) * t487;
t110 = -mrSges(7,2) * t216 + mrSges(7,3) * t172;
t111 = mrSges(7,1) * t216 - mrSges(7,3) * t173;
t23 = mrSges(7,1) * t83 - mrSges(7,3) * t57;
t24 = -mrSges(7,2) * t83 + mrSges(7,3) * t58;
t644 = (-t110 * t392 - t111 * t397) * qJD(6) + m(7) * t702 - t392 * t23 + t397 * t24;
t127 = -mrSges(6,1) * t487 + mrSges(6,2) * t441;
t640 = m(4) * t226 - t127 + t645 - t683;
t14 = Ifges(7,5) * t57 + Ifges(7,6) * t58 + Ifges(7,3) * t83;
t638 = mrSges(6,1) * t74 - mrSges(6,3) * t8 + Ifges(7,5) * t630 + Ifges(7,6) * t629 + Ifges(7,3) * t625 + t14 / 0.2e1 + t646 + (-t600 - t307 / 0.2e1) * Ifges(6,6) + (-t623 - t85 / 0.2e1) * Ifges(6,2) + (-t624 - t84 / 0.2e1) * Ifges(6,4);
t637 = t388 ^ 2;
t402 = qJD(2) ^ 2;
t633 = t16 / 0.2e1;
t628 = -t68 / 0.2e1;
t622 = Ifges(5,4) * t610 + Ifges(5,2) * t609 + Ifges(5,6) * t599;
t621 = Ifges(5,1) * t610 + Ifges(5,4) * t609 + Ifges(5,5) * t599;
t619 = t117 / 0.2e1;
t617 = t118 / 0.2e1;
t616 = -t172 / 0.2e1;
t614 = -t173 / 0.2e1;
t571 = t298 * Ifges(5,4);
t191 = t297 * Ifges(5,2) + t362 * Ifges(5,6) + t571;
t612 = t191 / 0.2e1;
t288 = Ifges(5,4) * t297;
t192 = t298 * Ifges(5,1) + t362 * Ifges(5,5) + t288;
t611 = t192 / 0.2e1;
t607 = -t216 / 0.2e1;
t602 = -t298 / 0.2e1;
t601 = t298 / 0.2e1;
t595 = pkin(2) * t388;
t585 = Ifges(4,4) * t395;
t584 = Ifges(4,4) * t400;
t583 = Ifges(5,4) * t394;
t582 = Ifges(5,4) * t399;
t560 = t260 * t392;
t554 = t337 * t388;
t553 = t338 * t388;
t552 = t348 * t392;
t549 = t385 * t388;
t548 = t386 * t388;
t547 = t388 * t394;
t515 = t388 * t543;
t532 = pkin(2) * t510 + pkin(9) * t515;
t513 = Ifges(5,5) * t203 + Ifges(5,6) * t204 + Ifges(5,3) * t315;
t512 = Ifges(4,5) * t330 + Ifges(4,6) * t329 + Ifges(4,3) * t375;
t498 = t546 / 0.2e1;
t143 = -t210 * t385 + t263 * t386;
t144 = t210 * t386 + t263 * t385;
t492 = t143 * pkin(5) + pkin(12) * t144;
t485 = mrSges(4,3) * t504;
t484 = t394 * t515;
t479 = t388 * t503;
t471 = t691 * pkin(4);
t470 = t692 * pkin(4);
t469 = t211 * pkin(4);
t327 = t411 * pkin(2);
t468 = pkin(9) * t554 - t327;
t328 = t412 * pkin(2);
t467 = pkin(9) * t553 - t328;
t464 = (pkin(4) * t394 + pkin(9)) * t388;
t456 = Ifges(5,1) * t399 - t583;
t454 = -Ifges(5,2) * t394 + t582;
t452 = Ifges(5,5) * t399 - Ifges(5,6) * t394;
t87 = t120 * t397 + t560;
t91 = t169 * t398 - t180 * t393;
t442 = t398 * t211 - t212 * t393;
t213 = -t248 * t392 - t397 * t544;
t437 = -t248 * t397 + t392 * t544;
t428 = t388 * (-mrSges(4,1) * t400 + mrSges(4,2) * t395);
t427 = (Ifges(4,2) * t400 + t585) * t388;
t419 = t395 * t637 * (Ifges(4,1) * t400 - t585);
t141 = -t208 * t385 + t262 * t386;
t142 = t208 * t386 + t262 * t385;
t418 = t688 * t141 + t142 * t694;
t417 = t688 * t143 + t144 * t694;
t195 = -t261 * t385 + t336 * t386;
t196 = t261 * t386 + t336 * t385;
t416 = t688 * t195 + t196 * t694;
t410 = t412 * t400;
t408 = t411 * t400;
t404 = -t572 / 0.2e1 - t573 / 0.2e1 - t576 / 0.2e1 + t574 / 0.2e1 + t628 + t619 + t581 / 0.2e1 + t695;
t370 = Ifges(4,4) * t504;
t322 = qJD(2) * t428;
t321 = -mrSges(4,2) * t376 + t485;
t272 = Ifges(4,1) * t505 + t376 * Ifges(4,5) + t370;
t271 = t376 * Ifges(4,6) + qJD(2) * t427;
t267 = mrSges(4,1) * t375 - mrSges(4,3) * t330;
t266 = -mrSges(4,2) * t375 + mrSges(4,3) * t329;
t252 = mrSges(5,1) * t362 - mrSges(5,3) * t298;
t251 = -mrSges(5,2) * t362 + mrSges(5,3) * t297;
t246 = -mrSges(4,1) * t329 + mrSges(4,2) * t330;
t238 = -t338 * t542 - t410;
t237 = t338 * t541 - t395 * t412;
t236 = -t337 * t542 - t408;
t235 = t337 * t541 - t409;
t209 = t338 * t395 + t390 * t410 - t400 * t465;
t207 = t337 * t395 + t390 * t408 + t400 * t466;
t194 = t391 * t500 + (t424 * qJD(2) + qJD(3) * t657) * t389;
t193 = t391 * t501 + (qJD(2) * t423 + qJD(3) * t422) * t389;
t187 = t195 * pkin(5);
t178 = -mrSges(6,2) * t357 + mrSges(6,3) * t487;
t156 = -mrSges(5,2) * t315 + mrSges(5,3) * t204;
t155 = mrSges(5,1) * t315 - mrSges(5,3) * t203;
t149 = -t398 * t239 + t240 * t393;
t139 = t141 * pkin(5);
t121 = -mrSges(5,1) * t204 + mrSges(5,2) * t203;
t104 = pkin(4) * t298 + t128;
t100 = qJD(4) * t211 + t194 * t399 + t394 * t479;
t99 = -qJD(4) * t212 - t194 * t394 + t399 * t479;
t88 = pkin(5) * t544 - t91;
t80 = qJD(6) * t437 - t125 * t392 + t397 * t501;
t79 = qJD(6) * t213 + t125 * t397 + t392 * t501;
t67 = -mrSges(6,2) * t307 + mrSges(6,3) * t85;
t61 = t105 * t398 - t564;
t60 = t105 * t393 + t537;
t41 = -mrSges(6,1) * t85 + mrSges(6,2) * t84;
t36 = t128 * t392 + t397 * t55;
t35 = t128 * t397 - t392 * t55;
t32 = qJD(5) * t442 + t100 * t398 + t393 * t99;
t30 = t104 * t392 + t397 * t61;
t29 = t104 * t397 - t392 * t61;
t26 = -pkin(5) * t501 - t28;
t18 = -qJD(6) * t87 + t193 * t397 - t32 * t392;
t17 = qJD(6) * t86 + t193 * t392 + t32 * t397;
t1 = [m(6) * (t120 * t8 + t32 * t56) + m(7) * (t17 * t22 + t18 * t21 + t2 * t87 + t3 * t86) + m(4) * (t122 * t261 + t194 * t227 + t258 * t336 + t287 * t479) + m(5) * (t100 * t130 + t129 * t99 + t211 * t45 + t212 * t44) + (-m(4) * t123 + m(5) * t109 + m(6) * t74 + t121 - t267 + t41) * t260 + (-qJDD(2) * t543 - t402 * t510) * mrSges(3,2) + (qJDD(2) * t510 - t402 * t543) * mrSges(3,1) - t640 * t193 + m(3) * (t391 ^ 2 * qJDD(1) + (t331 * t596 + t332 * t396) * t389) + t690 * (qJD(5) * t120 + t100 * t393 - t398 * t99) + t336 * t246 + t194 * t321 + t261 * t266 + t100 * t251 + t99 * t252 + t212 * t156 + t211 * t155 + (-m(2) - m(3) - m(4) - m(5) + t684) * g(3) + t32 * t178 - (-t653 + t668) * t442 + t86 * t23 + t87 * t24 + t17 * t110 + t18 * t111 + t120 * t67 + t322 * t479 + m(2) * qJDD(1); t667 * t67 + (t157 * t224 + t250 * t74 + t8 * t667 + t9 * t91 + t672 * t56 + (t149 + t28) * t55) * m(6) + (-Ifges(7,1) * t437 + Ifges(7,4) * t213) * t630 + (-Ifges(7,4) * t437 + Ifges(7,2) * t213) * t629 + (t2 * t213 - t21 * t79 + t22 * t80 + t3 * t437) * mrSges(7,3) + (-Ifges(7,5) * t437 + Ifges(7,6) * t213) * t625 + t6 * (-mrSges(7,1) * t213 - mrSges(7,2) * t437) - t437 * t633 - t638 * t440 + (mrSges(3,1) * t411 + t337 * mrSges(3,2) - m(4) * t468 - t236 * mrSges(4,1) - mrSges(4,3) * t554 - m(5) * (pkin(3) * t236 + t468) - (t236 * t399 + t337 * t547) * mrSges(5,1) - (-t236 * t394 + t337 * t545) * mrSges(5,2) + t684 * (-t235 * t401 + t236 * t384 + t337 * t464 - t327) + t649 * (t236 * t385 - t337 * t548) - t700 * (t236 * t386 + t337 * t549) + t686 * t235) * g(2) + (mrSges(3,1) * t412 + t338 * mrSges(3,2) - m(4) * t467 - t238 * mrSges(4,1) - mrSges(4,3) * t553 - m(5) * (pkin(3) * t238 + t467) - (t238 * t399 + t338 * t547) * mrSges(5,1) - (-t238 * t394 + t338 * t545) * mrSges(5,2) + t684 * (-t237 * t401 + t238 * t384 + t338 * t464 - t328) + t649 * (t238 * t385 - t338 * t548) - t700 * (t238 * t386 + t338 * t549) + t686 * t237) * g(1) + (-(mrSges(3,1) * t596 - mrSges(3,2) * t396) * t389 - m(4) * t532 - t300 * mrSges(4,1) - mrSges(4,3) * t515 - m(5) * (pkin(3) * t300 + t532) - (t300 * t399 + t484) * mrSges(5,1) - (-t300 * t394 + t399 * t515) * mrSges(5,2) + t684 * (pkin(4) * t484 - t299 * t401 + t300 * t384 + t532) + t649 * (t300 * t385 - t386 * t515) - t700 * (t300 * t386 + t385 * t515) + t686 * t299) * g(3) + (-t332 + t369) * mrSges(3,2) + (t498 * t666 + t654) * qJD(3) + t640 * t283 + (t637 * qJD(2) * (-Ifges(4,2) * t395 + t584) + t388 * t272) * t528 / 0.2e1 - (t518 + t513) * t544 / 0.2e1 + (t125 * t157 + t248 * t74 - t501 * t56 + t544 * t8) * mrSges(6,2) + t375 * (Ifges(4,3) * t390 + (Ifges(4,5) * t395 + Ifges(4,6) * t400) * t388) / 0.2e1 + t258 * t428 - t246 * t595 + (Ifges(6,5) * t125 + Ifges(6,3) * t501) * t598 + t330 * (Ifges(4,5) * t390 + (Ifges(4,1) * t395 + t584) * t388) / 0.2e1 + t689 * t126 + (Ifges(6,1) * t248 - Ifges(6,5) * t544) * t624 + (Ifges(6,1) * t125 + Ifges(6,5) * t501) * t604 + (t331 + t476) * mrSges(3,1) + (t122 * t544 - t123 * t546 - t226 * t500 - t227 * t501) * mrSges(4,3) + (Ifges(7,5) * t79 + Ifges(7,6) * t80) * t606 + (Ifges(4,4) * t330 + Ifges(4,2) * t329 + Ifges(4,6) * t375) * t544 / 0.2e1 + t9 * (-mrSges(6,1) * t544 - mrSges(6,3) * t248) + t45 * (-mrSges(5,1) * t544 - mrSges(5,3) * t340) + t44 * (mrSges(5,2) * t544 + mrSges(5,3) * t339) + Ifges(3,3) * qJDD(2) + t419 * t522 / 0.2e1 + t390 * t512 / 0.2e1 - t271 * t501 / 0.2e1 + t130 * (-mrSges(5,2) * t501 + mrSges(5,3) * t269) + t362 * (Ifges(5,5) * t268 + Ifges(5,6) * t269 + Ifges(5,3) * t501) / 0.2e1 + t297 * (Ifges(5,4) * t268 + Ifges(5,2) * t269 + Ifges(5,6) * t501) / 0.2e1 + t55 * (mrSges(6,1) * t501 - mrSges(6,3) * t125) + t129 * (mrSges(5,1) * t501 - mrSges(5,3) * t268) + (Ifges(6,4) * t248 - Ifges(6,6) * t544) * t623 + (Ifges(6,4) * t125 + Ifges(6,6) * t501) * t605 + t109 * (-mrSges(5,1) * t339 + mrSges(5,2) * t340) + t342 * t267 + t343 * t266 + t310 * t121 + t205 * (-mrSges(5,1) * t269 + mrSges(5,2) * t268) + t250 * t41 + t228 * t155 + t229 * t156 + t224 * t127 + t248 * t631 + t213 * t634 + t340 * t621 + t339 * t622 + (Ifges(5,4) * t340 + Ifges(5,2) * t339 - Ifges(5,6) * t544) * t609 + (Ifges(5,1) * t340 + Ifges(5,4) * t339 - Ifges(5,5) * t544) * t610 + t268 * t611 + t269 * t612 + t125 * t617 + (Ifges(5,5) * t340 + Ifges(5,6) * t339 - Ifges(5,3) * t544) * t599 + (Ifges(5,1) * t268 + Ifges(5,4) * t269 + Ifges(5,5) * t501) * t601 + t390 * t676 + t28 * t179 + (Ifges(7,1) * t79 + Ifges(7,4) * t80) * t613 - t390 * t677 + (Ifges(4,1) * t330 + Ifges(4,4) * t329 + Ifges(4,5) * t375) * t498 + t51 * t23 + t52 * t24 - t322 * t480 + t79 * t70 / 0.2e1 + t80 * t69 / 0.2e1 + t47 * (-mrSges(7,1) * t80 + mrSges(7,2) * t79) + (Ifges(6,5) * t248 - Ifges(6,3) * t544) * t600 + t88 * t19 + t91 * t66 + (t122 * t343 + t123 * t342 - t226 * t326 + t227 * t659 - t258 * t595 - t287 * t480) * m(4) + t659 * t321 + t26 * t93 + t662 * t251 + t663 * t252 + (t109 * t310 + t129 * t663 + t130 * t662 + t205 * t326 + t228 * t45 + t229 * t44) * m(5) + t664 * t326 + t329 * (Ifges(4,6) * t390 + t427) / 0.2e1 + (Ifges(7,4) * t79 + Ifges(7,2) * t80) * t615 + t567 * t149 + t672 * t178 + t673 * t111 + t674 * t110 + (t2 * t52 + t3 * t51 + t6 * t88 + (-t149 + t26) * t47 + t674 * t22 + t673 * t21) * m(7); (-Ifges(7,5) * t433 - Ifges(7,6) * t434) * t606 - (t19 - t66) * t439 + (t175 * t3 + t176 * t2 + t21 * t679 + t22 * t678 - t439 * t6 + t47 * t669) * m(7) + (t157 * t665 + t279 * t8 - t384 * t74 + t439 * t9 + t55 * t670 + t56 * t671) * m(6) + t638 * t347 - t487 * (-Ifges(6,4) * t277 - Ifges(6,2) * t276 + Ifges(6,6) * t505) / 0.2e1 - t441 * (-Ifges(6,1) * t277 - Ifges(6,4) * t276 + Ifges(6,5) * t505) / 0.2e1 - t357 * (-Ifges(6,5) * t277 - Ifges(6,6) * t276 + Ifges(6,3) * t505) / 0.2e1 - t55 * (mrSges(6,1) * t505 + mrSges(6,3) * t277) - t157 * (mrSges(6,1) * t276 - mrSges(6,2) * t277) + (t297 * t454 + t298 * t456 + t362 * t452) * qJD(4) / 0.2e1 + (-Ifges(7,4) * t433 - Ifges(7,2) * t434) * t615 + (-t321 + t485) * t226 + (t684 * (-t260 * t384 - t261 * t401) + t686 * t261 + t648 * t260) * g(3) + (t684 * (-t209 * t384 - t210 * t401) + t686 * t210 + t648 * t209) * g(1) + (t684 * (-t207 * t384 - t208 * t401) + t686 * t208 + t648 * t207) * g(2) + t109 * t460 + t512 - ((-Ifges(4,2) * t505 + t399 * t192 + t272 + t370) * t400 + t362 * (Ifges(5,3) * t395 + t400 * t452) + t298 * (Ifges(5,5) * t395 + t400 * t456) + t297 * (Ifges(5,6) * t395 + t400 * t454) + t666 * t395) * t530 / 0.2e1 + t676 - t677 + (-Ifges(7,1) * t433 - Ifges(7,4) * t434) * t613 - t402 * t419 / 0.2e1 + t689 * t265 + (mrSges(6,2) * t74 - mrSges(6,3) * t9 + t451 * t625 + t453 * t629 + t455 * t630 + t457 * t6 + t493 * t70 + 0.2e1 * t631) * t348 + (t703 * t130 + (t530 * t536 - t526) * t129 + t656) * mrSges(5,3) - t15 * t552 / 0.2e1 - t191 * t527 / 0.2e1 - t56 * (-mrSges(6,2) * t505 - mrSges(6,3) * t276) - t384 * t41 + (-pkin(3) * t109 - t129 * t167 - t130 * t168) * m(5) + t279 * t67 - t168 * t251 - t167 * t252 - t245 * t70 / 0.2e1 + t551 * t633 + t276 * t619 + t394 * t621 + t399 * t622 + t276 * t628 + (Ifges(7,5) * t245 + Ifges(7,6) * t244 + Ifges(7,3) * t276) * t607 + (Ifges(5,2) * t399 + t583) * t609 + (Ifges(5,1) * t394 + t582) * t610 + t526 * t611 + t478 * t612 + (Ifges(7,1) * t245 + Ifges(7,4) * t244 + Ifges(7,5) * t276) * t614 + (Ifges(7,4) * t245 + Ifges(7,2) * t244 + Ifges(7,6) * t276) * t616 - t277 * t618 + (Ifges(5,5) * t394 + Ifges(5,6) * t399) * t599 + t276 * t681 + (-mrSges(5,1) * t129 + mrSges(5,2) * t130) * t505 - t276 * t682 + t175 * t23 + t176 * t24 + t678 * t110 + t679 * t111 + t362 * t205 * (mrSges(5,1) * t394 + mrSges(5,2) * t399) + (t486 + t645) * t227 - (Ifges(6,1) * t604 + Ifges(6,4) * t605 + Ifges(6,5) * t598 + t514 + t617 - t651) * t264 + (t271 * t498 - t654) * qJD(2) + (-t252 * t526 - t251 * t527 + m(5) * ((-t129 * t399 - t130 * t394) * qJD(4) + t656) - t394 * t155 + t399 * t156) * pkin(10) - t661 * t69 / 0.2e1 + (mrSges(7,1) * t661 - mrSges(7,2) * t660) * t47 + (-t2 * t552 + t21 * t660 - t22 * t661 - t3 * t551) * mrSges(7,3) + t665 * t127 - pkin(3) * t121 + t669 * t93 + t670 * t179 + t671 * t178; t693 - (-Ifges(5,2) * t298 + t192 + t288) * t297 / 0.2e1 + (t575 / 0.2e1 + t404) * t441 + (-m(7) * (pkin(12) * t196 + t187 + t469) - m(6) * t469 - mrSges(5,1) * t211 + mrSges(5,2) * t212 + t416) * g(3) + t513 + (t129 * t297 + t130 * t298) * mrSges(5,3) + (0.2e1 * t602 * t683 - t127 * t298 + t653 * t398 + (m(6) * t8 + t67) * t393 + (t690 * t393 + (t110 * t397 - t111 * t392 + t178 + m(6) * t56 + m(7) * (-t21 * t392 + t22 * t397)) * t398) * qJD(5)) * pkin(4) + (-t691 * mrSges(5,1) - (-t208 * t399 - t262 * t394) * mrSges(5,2) - m(7) * (pkin(12) * t142 + t139 + t471) - m(6) * t471 + t418) * g(2) + (-t692 * mrSges(5,1) - (-t210 * t399 - t263 * t394) * mrSges(5,2) - m(7) * (t470 + t492) - m(6) * t470 + t417) * g(1) + t567 * t60 - t362 * (Ifges(5,5) * t297 - Ifges(5,6) * t298) / 0.2e1 - t205 * (mrSges(5,1) * t298 + mrSges(5,2) * t297) - t129 * t251 + t130 * t252 - m(6) * (-t55 * t60 + t56 * t61) + t191 * t601 + (Ifges(5,1) * t297 - t571) * t602 - t61 * t178 + t644 * (pkin(4) * t393 + pkin(12)) - t44 * mrSges(5,2) + t45 * mrSges(5,1) - m(7) * (t21 * t29 + t22 * t30 + t47 * t60) + (-t586 / 0.2e1 + t658) * t487 - t30 * t110 - t29 * t111 + t668 * (-pkin(4) * t398 - pkin(5)); (-m(7) * t139 + t418) * g(2) + (-m(7) * t492 + t417) * g(1) + (-m(7) * t187 + t416) * g(3) - m(7) * (t21 * t35 + t22 * t36 + t47 * t56) + ((-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t441 + t658) * t487 + t567 * t56 + t404 * t441 - t55 * t178 - t36 * t110 - t35 * t111 + ((-g(2) * t142 - g(3) * t196) * m(7) + t644) * pkin(12) - t668 * pkin(5) + t693; -t47 * (mrSges(7,1) * t173 + mrSges(7,2) * t172) + (Ifges(7,1) * t172 - t580) * t614 + t69 * t613 + (Ifges(7,5) * t172 - Ifges(7,6) * t173) * t607 - t21 * t110 + t22 * t111 - g(1) * ((-t144 * t392 + t209 * t397) * mrSges(7,1) + (-t144 * t397 - t209 * t392) * mrSges(7,2)) - g(2) * ((-t142 * t392 + t207 * t397) * mrSges(7,1) + (-t142 * t397 - t207 * t392) * mrSges(7,2)) - g(3) * ((-t196 * t392 + t256) * mrSges(7,1) + (-t196 * t397 - t560) * mrSges(7,2)) + (t172 * t21 + t173 * t22) * mrSges(7,3) + t14 + (-Ifges(7,2) * t173 + t170 + t70) * t616 + t646;];
tau  = t1;
