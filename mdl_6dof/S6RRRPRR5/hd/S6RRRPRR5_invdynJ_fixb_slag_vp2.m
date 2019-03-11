% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:39
% EndTime: 2019-03-09 18:21:22
% DurationCPUTime: 25.55s
% Computational Cost: add. (16691->863), mult. (35158->1084), div. (0->0), fcn. (24899->14), ass. (0->421)
t687 = -mrSges(4,1) + mrSges(5,2);
t392 = sin(qJ(5));
t393 = sin(qJ(3));
t497 = qJD(3) * t393;
t484 = pkin(2) * t497;
t576 = cos(qJ(3));
t577 = cos(qJ(2));
t449 = t576 * t577;
t394 = sin(qJ(2));
t499 = qJD(1) * t394;
t281 = -qJD(1) * t449 + t393 * t499;
t307 = t393 * t577 + t394 * t576;
t282 = t307 * qJD(1);
t221 = pkin(3) * t282 + qJ(4) * t281;
t372 = pkin(2) * t499;
t201 = t221 + t372;
t276 = t282 * pkin(9);
t155 = t201 + t276;
t400 = -pkin(8) - pkin(7);
t331 = t400 * t394;
t309 = qJD(1) * t331;
t489 = t577 * pkin(7);
t332 = pkin(8) * t577 + t489;
t310 = t332 * qJD(1);
t473 = t576 * t310;
t229 = t309 * t393 + t473;
t573 = pkin(4) * t281;
t188 = t229 - t573;
t397 = cos(qJ(5));
t88 = t397 * t155 + t392 * t188;
t686 = t392 * t484 - t88;
t391 = sin(qJ(6));
t494 = qJD(6) * t391;
t496 = qJD(5) * t392;
t396 = cos(qJ(6));
t508 = t396 * t397;
t623 = qJD(5) + qJD(6);
t232 = -t391 * t496 - t392 * t494 + t508 * t623;
t430 = t391 * t392 - t508;
t644 = t430 * t282;
t675 = t232 - t644;
t431 = t391 * t397 + t396 * t392;
t233 = t623 * t431;
t417 = t431 * t282;
t685 = t233 + t417;
t684 = t281 * pkin(5) + pkin(10) * t496;
t659 = -m(7) - m(5);
t176 = t397 * t188;
t568 = pkin(10) * t282;
t488 = t576 * pkin(2);
t367 = -t488 - pkin(3);
t353 = -pkin(9) + t367;
t631 = -t353 * t496 + t397 * t484;
t683 = t631 - t176 - (-t155 - t568) * t392 + t684;
t528 = t282 * t397;
t262 = pkin(10) * t528;
t559 = -pkin(10) + t353;
t293 = t559 * t397;
t682 = qJD(5) * t293 - t262 + t686;
t291 = qJD(2) * pkin(2) + t309;
t227 = t393 * t291 + t473;
t179 = t227 - t573;
t165 = t397 * t179;
t171 = t221 + t276;
t601 = pkin(3) + pkin(9);
t470 = t601 * t496;
t681 = t470 - t165 - (-t171 - t568) * t392 + t684;
t558 = -pkin(10) - t601;
t323 = t558 * t397;
t92 = t397 * t171 + t392 * t179;
t680 = qJD(5) * t323 - t262 - t92;
t390 = qJ(2) + qJ(3);
t384 = cos(t390);
t389 = qJ(5) + qJ(6);
t383 = cos(t389);
t519 = t383 * t384;
t381 = sin(t389);
t524 = t381 * t384;
t555 = mrSges(6,2) * t397;
t679 = -mrSges(7,1) * t524 - mrSges(7,2) * t519 - t384 * t555;
t388 = qJD(2) + qJD(3);
t380 = t388 * qJ(4);
t157 = t179 + t380;
t244 = t281 * t397 - t388 * t392;
t120 = -pkin(5) * t244 + t157;
t245 = t281 * t392 + t388 * t397;
t147 = t244 * t391 + t245 * t396;
t277 = qJD(5) + t282;
t267 = t576 * t291;
t511 = t393 * t310;
t226 = -t267 + t511;
t564 = t282 * pkin(4);
t424 = t226 + t564;
t666 = t424 + qJD(4);
t139 = -t388 * t601 + t666;
t386 = t577 * pkin(2);
t368 = t386 + pkin(1);
t330 = t368 * qJD(1);
t407 = -t282 * qJ(4) - t330;
t143 = t281 * t601 + t407;
t77 = t397 * t139 - t143 * t392;
t61 = -pkin(10) * t245 + t77;
t50 = pkin(5) * t277 + t61;
t78 = t139 * t392 + t143 * t397;
t62 = pkin(10) * t244 + t78;
t537 = t391 * t62;
t19 = t396 * t50 - t537;
t535 = t396 * t62;
t20 = t391 * t50 + t535;
t456 = t396 * t244 - t245 * t391;
t465 = qJD(2) * t577;
t312 = qJD(1) * t465 + t394 * qJDD(1);
t492 = qJD(1) * qJD(2);
t671 = t577 * qJDD(1) - t394 * t492;
t186 = qJD(3) * t281 - t576 * t312 - t393 * t671;
t184 = qJDD(5) - t186;
t174 = qJDD(6) + t184;
t406 = t307 * qJD(3);
t187 = qJD(1) * t406 + t393 * t312 - t576 * t671;
t387 = qJDD(2) + qJDD(3);
t114 = qJD(5) * t244 + t187 * t392 + t387 * t397;
t115 = -qJD(5) * t245 + t187 * t397 - t387 * t392;
t39 = qJD(6) * t456 + t114 * t396 + t115 * t391;
t40 = -qJD(6) * t147 - t114 * t391 + t115 * t396;
t490 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t174;
t545 = Ifges(7,4) * t147;
t265 = qJD(6) + t277;
t588 = -t265 / 0.2e1;
t596 = -t147 / 0.2e1;
t533 = qJDD(1) * pkin(1);
t272 = -pkin(2) * t671 - t533;
t403 = t186 * qJ(4) - t282 * qJD(4) + t272;
t48 = t187 * t601 + t403;
t296 = t312 * pkin(7);
t238 = qJDD(2) * pkin(2) - t312 * pkin(8) - t296;
t295 = t671 * pkin(7);
t243 = pkin(8) * t671 + t295;
t464 = qJD(3) * t576;
t97 = t238 * t576 - t393 * t243 - t291 * t497 - t310 * t464;
t410 = qJDD(4) - t97;
t58 = -t186 * pkin(4) - t387 * t601 + t410;
t15 = -qJD(5) * t78 - t392 * t48 + t397 * t58;
t7 = pkin(5) * t184 - pkin(10) * t114 + t15;
t495 = qJD(5) * t397;
t14 = t139 * t495 - t143 * t496 + t392 * t58 + t397 * t48;
t8 = pkin(10) * t115 + t14;
t3 = qJD(6) * t19 + t391 * t7 + t396 * t8;
t4 = -qJD(6) * t20 - t391 * t8 + t396 * t7;
t616 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t678 = t490 + t616 + (Ifges(7,5) * t456 - Ifges(7,6) * t147) * t588 + (t147 * t20 + t19 * t456) * mrSges(7,3) - t120 * (mrSges(7,1) * t147 + mrSges(7,2) * t456) + (Ifges(7,1) * t456 - t545) * t596;
t677 = -Ifges(4,5) + Ifges(5,4);
t676 = -Ifges(4,6) + Ifges(5,5);
t474 = t576 * t309;
t230 = t474 - t511;
t451 = pkin(2) * t464;
t346 = t451 + qJD(4);
t673 = -t346 + t230;
t570 = pkin(5) * t397;
t366 = pkin(4) + t570;
t672 = pkin(5) * t495 + t282 * t366 + t511;
t305 = t393 * t394 - t449;
t670 = t484 - t229;
t382 = sin(t390);
t668 = t687 * t384 + (mrSges(4,2) - mrSges(5,3)) * t382;
t548 = Ifges(6,4) * t245;
t128 = Ifges(6,2) * t244 + Ifges(6,6) * t277 + t548;
t438 = mrSges(6,1) * t397 - mrSges(6,2) * t392;
t667 = t157 * t438 - t128 * t397 / 0.2e1;
t274 = Ifges(4,4) * t281;
t653 = Ifges(6,3) * t277;
t654 = Ifges(6,6) * t244;
t665 = Ifges(4,1) * t282 + Ifges(4,5) * t388 + Ifges(6,5) * t245 + Ifges(7,5) * t147 + Ifges(7,6) * t456 + Ifges(7,3) * t265 - t274 + t653 + t654;
t140 = Ifges(7,4) * t456;
t664 = -Ifges(7,2) * t147 + t140;
t193 = t281 * pkin(3) + t407;
t663 = t330 * mrSges(4,1) + t193 * mrSges(5,2);
t662 = t19 * t685 - t20 * t675 - t3 * t431 + t4 * t430;
t661 = t77 * mrSges(6,1) + t19 * mrSges(7,1) - t330 * mrSges(4,2) - t78 * mrSges(6,2) - t20 * mrSges(7,2) - t193 * mrSges(5,3);
t607 = t39 / 0.2e1;
t606 = t40 / 0.2e1;
t600 = t114 / 0.2e1;
t599 = t115 / 0.2e1;
t594 = t174 / 0.2e1;
t593 = t184 / 0.2e1;
t580 = t387 / 0.2e1;
t658 = t671 / 0.2e1;
t610 = m(7) * pkin(5);
t651 = -mrSges(6,1) - t610;
t292 = t559 * t392;
t218 = -t292 * t391 + t293 * t396;
t650 = qJD(6) * t218 + t391 * t683 + t396 * t682;
t219 = t292 * t396 + t293 * t391;
t649 = -qJD(6) * t219 - t391 * t682 + t396 * t683;
t322 = t558 * t392;
t240 = -t322 * t391 + t323 * t396;
t648 = qJD(6) * t240 + t391 * t681 + t396 * t680;
t241 = t322 * t396 + t323 * t391;
t647 = -qJD(6) * t241 - t391 * t680 + t396 * t681;
t161 = mrSges(5,1) * t187 - mrSges(5,3) * t387;
t49 = -mrSges(6,1) * t115 + mrSges(6,2) * t114;
t646 = t49 - t161;
t214 = t430 * t305;
t385 = t392 * pkin(5);
t360 = qJ(4) + t385;
t156 = -mrSges(6,1) * t244 + mrSges(6,2) * t245;
t557 = mrSges(5,1) * t281;
t252 = -mrSges(5,3) * t388 + t557;
t643 = t156 - t252;
t444 = -qJ(4) * t307 - t368;
t185 = t305 * t601 + t444;
t248 = -t576 * t331 + t332 * t393;
t211 = pkin(4) * t307 + t248;
t198 = t392 * t211;
t103 = t397 * t185 + t198;
t554 = mrSges(4,3) * t281;
t250 = -mrSges(4,2) * t388 - t554;
t642 = t250 - t252;
t553 = mrSges(4,3) * t282;
t556 = mrSges(5,1) * t282;
t641 = t388 * t687 + t553 + t556;
t640 = t346 - t474 + t672;
t395 = sin(qJ(1));
t513 = t392 * t395;
t481 = t384 * t513;
t399 = -pkin(10) - pkin(9);
t520 = t382 * t399;
t639 = pkin(5) * t481 + t395 * t520;
t398 = cos(qJ(1));
t512 = t392 * t398;
t480 = t384 * t512;
t638 = pkin(5) * t480 + t398 * t520;
t515 = t384 * t398;
t521 = t382 * t398;
t637 = pkin(3) * t515 + qJ(4) * t521;
t523 = t382 * t392;
t636 = pkin(5) * t523 - t384 * t399;
t635 = t564 - t673;
t634 = qJD(4) - t267 + t672;
t630 = t577 * t295 + t296 * t394;
t466 = qJD(1) * t577;
t569 = pkin(7) * t394;
t629 = (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t466) * t569 + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t499) * t489;
t72 = mrSges(6,1) * t184 - mrSges(6,3) * t114;
t73 = -mrSges(6,2) * t184 + mrSges(6,3) * t115;
t628 = t392 * t73 + t397 * t72;
t626 = t14 * t392 + t15 * t397;
t625 = g(1) * t398 + g(2) * t395;
t624 = -t226 - qJD(4);
t622 = 0.2e1 * t580;
t423 = mrSges(3,1) * t394 + mrSges(3,2) * t577;
t550 = Ifges(3,4) * t394;
t621 = pkin(1) * t423 - t394 * (Ifges(3,1) * t577 - t550) / 0.2e1;
t335 = qJ(4) * t515;
t620 = -m(6) * t335 - mrSges(6,1) * t480 - mrSges(5,2) * t521 - mrSges(5,3) * t515 + t398 * t679;
t516 = t384 * t395;
t333 = qJ(4) * t516;
t522 = t382 * t395;
t619 = -m(6) * t333 - mrSges(6,1) * t481 - mrSges(5,2) * t522 - mrSges(5,3) * t516 + t395 * t679;
t329 = -mrSges(3,1) * t577 + t394 * mrSges(3,2);
t618 = -m(3) * pkin(1) - mrSges(2,1) + t329 + t668;
t356 = t384 * mrSges(7,3);
t617 = -mrSges(6,1) * t523 - t384 * mrSges(6,3) - t356 + t668 + (-mrSges(7,1) * t381 - mrSges(7,2) * t383 - t555) * t382;
t450 = -m(6) * t601 - mrSges(6,3);
t425 = t450 * t382;
t574 = pkin(2) * t394;
t615 = m(6) * t574 + t382 * mrSges(7,3) - t425 + t659 * (-pkin(3) * t382 - t574);
t614 = t15 * mrSges(6,1) - t14 * mrSges(6,2);
t552 = mrSges(6,3) * t244;
t190 = -mrSges(6,2) * t277 + t552;
t433 = t392 * t77 - t397 * t78;
t404 = -qJD(5) * t433 + t626;
t613 = m(6) * t404 + t190 * t495 + t628;
t612 = -m(3) * pkin(7) - m(6) * (pkin(4) - t400) - m(7) * t366 - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t609 = Ifges(7,4) * t607 + Ifges(7,2) * t606 + Ifges(7,6) * t594;
t608 = Ifges(7,1) * t607 + Ifges(7,4) * t606 + Ifges(7,5) * t594;
t605 = Ifges(6,1) * t600 + Ifges(6,4) * t599 + Ifges(6,5) * t593;
t69 = Ifges(7,2) * t456 + Ifges(7,6) * t265 + t545;
t604 = -t69 / 0.2e1;
t70 = Ifges(7,1) * t147 + Ifges(7,5) * t265 + t140;
t603 = -t70 / 0.2e1;
t602 = t70 / 0.2e1;
t598 = -t456 / 0.2e1;
t597 = t456 / 0.2e1;
t595 = t147 / 0.2e1;
t591 = -t244 / 0.2e1;
t590 = -t245 / 0.2e1;
t589 = t245 / 0.2e1;
t587 = t265 / 0.2e1;
t586 = -t277 / 0.2e1;
t585 = -t281 / 0.2e1;
t584 = t281 / 0.2e1;
t583 = -t282 / 0.2e1;
t582 = t282 / 0.2e1;
t579 = -t388 / 0.2e1;
t578 = t388 / 0.2e1;
t575 = pkin(2) * t393;
t572 = pkin(5) * t245;
t571 = pkin(5) * t396;
t565 = g(3) * t384;
t362 = t384 * pkin(3);
t551 = mrSges(6,3) * t245;
t549 = Ifges(4,4) * t282;
t547 = Ifges(6,4) * t392;
t546 = Ifges(6,4) * t397;
t544 = Ifges(5,6) * t282;
t235 = qJD(2) * t307 + t406;
t532 = t235 * t392;
t531 = t235 * t397;
t529 = t282 * t392;
t526 = t305 * t392;
t525 = t305 * t397;
t354 = t382 * qJ(4);
t518 = t383 * t395;
t517 = t383 * t398;
t509 = t395 * t397;
t507 = t397 * t398;
t258 = -t381 * t395 + t382 * t517;
t259 = t381 * t521 + t518;
t506 = t258 * mrSges(7,1) - t259 * mrSges(7,2);
t260 = t381 * t398 + t382 * t518;
t261 = -t381 * t522 + t517;
t505 = t260 * mrSges(7,1) + t261 * mrSges(7,2);
t500 = t362 + t354;
t498 = qJD(2) * t394;
t375 = pkin(2) * t498;
t482 = Ifges(3,4) * t577;
t479 = Ifges(6,5) * t114 + Ifges(6,6) * t115 + Ifges(6,3) * t184;
t476 = t386 + t500;
t467 = t525 / 0.2e1;
t462 = -t496 / 0.2e1;
t359 = qJ(4) + t575;
t461 = -pkin(10) * t305 - t185;
t162 = -t186 * mrSges(5,1) + t387 * mrSges(5,2);
t457 = -t368 - t354;
t336 = t398 * t368;
t455 = -t395 * t400 + t336;
t448 = -pkin(3) * t522 + t333;
t447 = -pkin(3) * t521 + t335;
t439 = mrSges(4,1) * t382 + mrSges(4,2) * t384;
t436 = Ifges(6,1) * t392 + t546;
t435 = Ifges(6,2) * t397 + t547;
t434 = Ifges(6,5) * t392 + Ifges(6,6) * t397;
t199 = t397 * t211;
t81 = pkin(5) * t307 + t392 * t461 + t199;
t89 = pkin(10) * t525 + t103;
t37 = -t391 * t89 + t396 * t81;
t38 = t391 * t81 + t396 * t89;
t191 = mrSges(6,1) * t277 - t551;
t432 = t190 * t397 - t191 * t392;
t427 = t500 + t636;
t422 = t577 * Ifges(3,2) + t550;
t421 = Ifges(3,5) * t577 - Ifges(3,6) * t394;
t268 = t382 * t507 - t513;
t270 = t382 * t509 + t512;
t249 = t393 * t331 + t332 * t576;
t420 = t305 * t495 + t532;
t419 = t305 * t496 - t531;
t234 = t305 * t388;
t418 = qJ(4) * t234 - qJD(4) * t307 + t375;
t311 = qJD(2) * t331;
t409 = qJD(2) * t332;
t149 = qJD(3) * t249 + t393 * t311 + t576 * t409;
t111 = -t234 * pkin(4) + t149;
t86 = t235 * t601 + t418;
t23 = t392 * t111 - t185 * t496 + t211 * t495 + t397 * t86;
t96 = t393 * t238 + t576 * t243 + t291 * t464 - t310 * t497;
t148 = -t576 * t311 - t331 * t464 + t332 * t497 + t393 * t409;
t90 = -t387 * qJ(4) - t388 * qJD(4) - t96;
t110 = -pkin(4) * t235 - t148;
t59 = -pkin(4) * t187 - t90;
t242 = Ifges(6,4) * t244;
t129 = Ifges(6,1) * t245 + Ifges(6,5) * t277 + t242;
t200 = -pkin(3) * t388 - t624;
t206 = Ifges(5,5) * t388 + Ifges(5,3) * t281 - t544;
t273 = Ifges(5,6) * t281;
t207 = Ifges(5,4) * t388 - Ifges(5,2) * t282 + t273;
t208 = -Ifges(4,2) * t281 + Ifges(4,6) * t388 + t549;
t213 = -t380 - t227;
t34 = -pkin(5) * t115 + t59;
t43 = t114 * Ifges(6,4) + t115 * Ifges(6,2) + t184 * Ifges(6,6);
t94 = -t387 * pkin(3) + t410;
t402 = (Ifges(5,3) * t585 + t434 * t586 + t435 * t591 + t436 * t590 + t579 * t676 + t663 + t667) * t282 + t676 * t187 + t677 * t186 + (-Ifges(4,1) * t583 - Ifges(6,5) * t590 - Ifges(7,5) * t596 + Ifges(5,2) * t582 - Ifges(6,6) * t591 - Ifges(7,6) * t598 - Ifges(6,3) * t586 - Ifges(7,3) * t588 + t579 * t677 + t661) * t281 + (-t549 + t206) * t583 + t417 * t603 + (t273 + t207) * t585 + t675 * t604 + t667 * qJD(5) + (mrSges(7,1) * t675 - mrSges(7,2) * t685) * t120 + (-t626 + (-t495 - t528) * t78 + (t496 + t529) * t77) * mrSges(6,3) + (-t529 / 0.2e1 + t462) * t129 - (t244 * t435 + t245 * t436 + t277 * t434) * qJD(5) / 0.2e1 + (Ifges(7,5) * t417 - Ifges(7,6) * t644) * t588 + (t544 + t208) * t582 + (Ifges(5,1) + Ifges(4,3)) * t387 + (-Ifges(4,2) * t282 - t274 + t665) * t584 - t213 * t556 + t227 * t553 + t226 * t554 + t200 * t557 + t59 * (mrSges(6,1) * t392 + t555) - t392 * t43 / 0.2e1 - t96 * mrSges(4,2) + t97 * mrSges(4,1) - t90 * mrSges(5,3) + t94 * mrSges(5,2) + t397 * t605 + (Ifges(7,1) * t417 - Ifges(7,4) * t644) * t596 + (Ifges(7,4) * t417 - Ifges(7,2) * t644) * t598 - t431 * t609 - t430 * t608 + t34 * (mrSges(7,1) * t431 - mrSges(7,2) * t430) + (-Ifges(7,4) * t430 - Ifges(7,2) * t431) * t606 + (-Ifges(7,1) * t430 - Ifges(7,4) * t431) * t607 + (-Ifges(7,5) * t430 - Ifges(7,6) * t431) * t594 + (-Ifges(7,5) * t233 - Ifges(7,6) * t232) * t587 + (-Ifges(7,1) * t233 - Ifges(7,4) * t232) * t595 + (-Ifges(7,4) * t233 - Ifges(7,2) * t232) * t597 + (Ifges(6,5) * t397 - Ifges(6,6) * t392) * t593 + (-Ifges(6,2) * t392 + t546) * t599 + (Ifges(6,1) * t397 - t547) * t600 - t233 * t602 + t662 * mrSges(7,3);
t370 = Ifges(3,4) * t466;
t361 = t384 * pkin(9);
t324 = t359 + t385;
t321 = mrSges(7,2) * t524;
t280 = Ifges(3,1) * t499 + Ifges(3,5) * qJD(2) + t370;
t279 = Ifges(3,6) * qJD(2) + qJD(1) * t422;
t271 = -t382 * t513 + t507;
t269 = t382 * t512 + t509;
t225 = pkin(3) * t305 + t444;
t224 = -mrSges(5,2) * t281 - mrSges(5,3) * t282;
t223 = mrSges(4,1) * t281 + mrSges(4,2) * t282;
t215 = t431 * t305;
t212 = -t305 * pkin(4) + t249;
t160 = -mrSges(4,2) * t387 - mrSges(4,3) * t187;
t159 = mrSges(4,1) * t387 + mrSges(4,3) * t186;
t158 = -t305 * t366 + t249;
t123 = mrSges(7,1) * t265 - mrSges(7,3) * t147;
t122 = -mrSges(7,2) * t265 + mrSges(7,3) * t456;
t118 = pkin(3) * t235 + t418;
t107 = t397 * t111;
t102 = -t185 * t392 + t199;
t91 = -t171 * t392 + t165;
t87 = -t155 * t392 + t176;
t84 = -mrSges(7,1) * t456 + mrSges(7,2) * t147;
t75 = t187 * pkin(3) + t403;
t71 = pkin(5) * t419 + t110;
t64 = -t233 * t305 - t235 * t430;
t63 = -t214 * t623 + t431 * t235;
t26 = -mrSges(7,2) * t174 + mrSges(7,3) * t40;
t25 = mrSges(7,1) * t174 - mrSges(7,3) * t39;
t24 = -qJD(5) * t103 - t392 * t86 + t107;
t22 = t396 * t61 - t537;
t21 = -t391 * t61 - t535;
t18 = -pkin(10) * t419 + t23;
t16 = -pkin(5) * t234 + t107 + (-pkin(10) * t235 - t86) * t392 + (t397 * t461 - t198) * qJD(5);
t13 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t6 = -qJD(6) * t38 + t16 * t396 - t18 * t391;
t5 = qJD(6) * t37 + t16 * t391 + t18 * t396;
t1 = [t312 * t482 / 0.2e1 + (t206 / 0.2e1 - t208 / 0.2e1 + mrSges(5,1) * t213 - mrSges(4,3) * t227 - Ifges(4,4) * t582 + Ifges(5,6) * t583 + Ifges(5,3) * t584 - Ifges(4,2) * t585 + t676 * t578 - t663) * t235 + (t272 * mrSges(4,1) + t90 * mrSges(5,1) - t75 * mrSges(5,2) - t96 * mrSges(4,3) + t128 * t462 + t434 * t593 + t435 * t599 + t436 * t600 - t59 * t438 + (Ifges(5,3) + Ifges(4,2)) * t187 + (Ifges(5,6) + Ifges(4,4)) * t186 + t676 * t622) * t305 + (Ifges(7,5) * t215 - Ifges(7,6) * t214) * t594 + t34 * (mrSges(7,1) * t214 + mrSges(7,2) * t215) + m(4) * (-t272 * t368 - t330 * t375) + (Ifges(7,1) * t215 - Ifges(7,4) * t214) * t607 + (-t19 * t63 + t20 * t64 - t214 * t3 - t215 * t4) * mrSges(7,3) + (Ifges(7,4) * t215 - Ifges(7,2) * t214) * t606 + (qJD(1) * (-Ifges(3,2) * t394 + t482) + t280) * t465 / 0.2e1 + (t14 * t525 - t15 * t526 - t419 * t78 - t420 * t77) * mrSges(6,3) + t422 * t658 + (m(4) * t226 + m(5) * t200 + t641) * t149 + (-t271 * mrSges(6,1) - t261 * mrSges(7,1) + t270 * mrSges(6,2) + t260 * mrSges(7,2) + ((m(4) - t659) * t400 + t612) * t398 + (-m(7) * (-t360 * t382 - t368) - m(5) * (t457 - t362) - m(6) * t457 + m(4) * t368 + (-m(7) * (-pkin(3) + t399) + mrSges(7,3) - t450) * t384 - t618) * t395) * g(1) + (-m(6) * (pkin(9) * t515 + t336 + t637) - t269 * mrSges(6,1) - t268 * mrSges(6,2) - mrSges(6,3) * t515 - t259 * mrSges(7,1) - t258 * mrSges(7,2) - m(4) * t455 + t659 * (t455 + t637) + t612 * t395 + (-m(7) * t636 - t356 + t618) * t398) * g(2) + (Ifges(3,1) * t312 + Ifges(3,4) * t658) * t394 + (-t653 / 0.2e1 - t661 - Ifges(4,1) * t582 + Ifges(5,2) * t583 + Ifges(5,6) * t584 - Ifges(4,4) * t585 - Ifges(7,3) * t587 - Ifges(6,5) * t589 - Ifges(7,5) * t595 - Ifges(7,6) * t597 - mrSges(4,3) * t226 + t677 * t578 - mrSges(5,1) * t200 - t654 / 0.2e1 + t207 / 0.2e1 - t665 / 0.2e1) * t234 + (-Ifges(4,4) * t187 + Ifges(4,5) * t387 + t479 + t490) * t307 / 0.2e1 + (t94 * mrSges(5,1) + t272 * mrSges(4,2) - t97 * mrSges(4,3) - t75 * mrSges(5,3) + Ifges(4,5) * t580 + Ifges(6,5) * t600 + Ifges(7,5) * t607 + Ifges(6,6) * t599 + Ifges(7,6) * t606 + Ifges(6,3) * t593 + Ifges(7,3) * t594 + (-Ifges(5,6) - Ifges(4,4) / 0.2e1) * t187 - Ifges(5,2) * t186 - t622 * Ifges(5,4) + t614 + t616) * t307 + m(5) * (t118 * t193 + t225 * t75) + (-m(4) * t227 + m(5) * t213 - t642) * t148 + (-m(4) * t97 + m(5) * t94 - t159 + t162) * t248 - t621 * t492 - t186 * Ifges(4,1) * t307 + t244 * (Ifges(6,4) * t420 - Ifges(6,2) * t419) / 0.2e1 + (t421 * qJD(2) / 0.2e1 - t629) * qJD(2) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t630) + (Ifges(6,1) * t420 - Ifges(6,4) * t419) * t589 + m(6) * (t102 * t15 + t103 * t14 + t110 * t157 + t212 * t59 + t23 * t78 + t24 * t77) + m(7) * (t120 * t71 + t158 * t34 + t19 * t6 + t20 * t5 + t3 * t38 + t37 * t4) + (t312 * t569 + t489 * t671 + t630) * mrSges(3,3) + t577 * (Ifges(3,4) * t312 + Ifges(3,2) * t671) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t671 + t312 * mrSges(3,2)) + (-mrSges(3,1) * t569 - mrSges(3,2) * t489 + Ifges(3,5) * t394 + Ifges(3,6) * t577) * qJDD(2) + t277 * (Ifges(6,5) * t420 - Ifges(6,6) * t419) / 0.2e1 + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t595 - t279 * t498 / 0.2e1 + t223 * t375 + (m(4) * t96 - m(5) * t90 + t160 - t161) * t249 - t329 * t533 + t128 * t531 / 0.2e1 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t597 + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t587 - t368 * (mrSges(4,1) * t187 - mrSges(4,2) * t186) + Ifges(2,3) * qJDD(1) + t118 * t224 + t225 * (-mrSges(5,2) * t187 + mrSges(5,3) * t186) + t212 * t49 + t23 * t190 + t24 * t191 + t158 * t13 + t110 * t156 + t120 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t5 * t122 + t6 * t123 + t102 * t72 + t103 * t73 + t71 * t84 + t64 * t69 / 0.2e1 + t37 * t25 + t38 * t26 + t526 * t605 + t215 * t608 - t214 * t609 + t63 * t602 + t157 * (mrSges(6,1) * t419 + mrSges(6,2) * t420) + t43 * t467 + (qJD(5) * t467 + t532 / 0.2e1) * t129; t686 * t190 + (-t226 * t229 - t227 * t230 + t330 * t372 + (t576 * t97 + t393 * t96 + (t226 * t393 + t227 * t576) * qJD(3)) * pkin(2)) * m(4) - (-Ifges(3,2) * t499 + t280 + t370) * t466 / 0.2e1 + (t359 * t59 + (t392 * t78 + t397 * t77) * t484 - t77 * t87 - t78 * t88 + t635 * t157) * m(6) + t635 * t156 + (-m(7) * (t335 + t638) - m(5) * t335 + t615 * t398 + t620) * g(1) + (-m(7) * (t333 + t639) - m(5) * t333 + t615 * t395 + t619) * g(2) + t640 * t84 + (-t193 * t201 + t200 * t670 + t213 * t673 - t359 * t90 + t367 * t94) * m(5) + t670 * t641 + (-m(4) * t386 - m(5) * t476 + t329 - m(6) * (t361 + t476) - m(7) * (t386 + t427) + t617) * g(3) + t613 * t353 - t642 * t230 + (-t87 + t631) * t191 + t649 * t123 + (t120 * t640 + t19 * t649 + t20 * t650 + t218 * t4 + t219 * t3 + t324 * t34) * m(7) + t650 * t122 + t646 * t359 + Ifges(3,6) * t671 + (m(4) * t574 + t423 + t439) * t625 + t402 + (qJD(1) * t621 + t629) * qJD(1) + t279 * t499 / 0.2e1 - t421 * t492 / 0.2e1 - t223 * t372 + t367 * t162 - t346 * t252 + t324 * t13 + Ifges(3,5) * t312 - t295 * mrSges(3,2) - t296 * mrSges(3,1) - t201 * t224 + t218 * t25 + t219 * t26 + Ifges(3,3) * qJDD(2) + t159 * t488 + t160 * t575 + t250 * t451; -pkin(3) * t162 + t360 * t13 + t424 * t156 - t92 * t190 - t221 * t224 + t240 * t25 + t241 * t26 + t402 + t634 * t84 + t625 * t439 - t641 * t227 + t642 * t226 + (t470 - t91) * t191 + t647 * t123 + t648 * t122 + t643 * qJD(4) + t646 * qJ(4) + (t120 * t634 + t19 * t647 + t20 * t648 + t240 * t4 + t241 * t3 + t34 * t360) * m(7) + (t59 * qJ(4) + t666 * t157 - t77 * t91 - t78 * t92) * m(6) + (-pkin(3) * t94 - qJ(4) * t90 - t193 * t221 - t200 * t227 + t213 * t624) * m(5) - t613 * t601 + (-m(5) * t448 - m(7) * (t448 + t639) + mrSges(7,3) * t522 - t395 * t425 + t619) * g(2) + (-m(5) * t447 - m(7) * (t447 + t638) + mrSges(7,3) * t521 - t398 * t425 + t620) * g(1) + (-m(6) * (t361 + t500) - m(5) * t500 - m(7) * t427 + t617) * g(3); -t430 * t25 + t431 * t26 - t685 * t123 + t675 * t122 + t432 * qJD(5) + (-t84 - t643) * t388 + (t224 + t432) * t282 + t162 + (-t120 * t388 - t662) * m(7) + (-t157 * t388 - t282 * t433 + t404) * m(6) + (t193 * t282 + t213 * t388 + t94) * m(5) + (-t382 * t625 + t565) * (m(6) - t659) + t628; t456 * t603 + t664 * t598 + t614 + (t552 - t190) * t77 - t147 * t604 + (-Ifges(6,2) * t245 + t129 + t242) * t591 + (t551 + t191) * t78 + (-mrSges(6,2) * t271 + t270 * t651 - t505) * g(2) + (mrSges(6,2) * t269 + t268 * t651 - t506) * g(1) + (-t123 * t494 + t26 * t391) * pkin(5) + (qJD(6) * t571 - t22) * t122 + t438 * t565 + t25 * t571 + t678 - g(3) * (t321 + (-m(7) * t570 - mrSges(7,1) * t383) * t384) + t479 - t157 * (mrSges(6,1) * t245 + mrSges(6,2) * t244) - t21 * t123 + (t3 * t391 + t396 * t4 + (-t19 * t391 + t20 * t396) * qJD(6)) * t610 + (Ifges(6,5) * t244 - Ifges(6,6) * t245) * t586 + t128 * t589 + (Ifges(6,1) * t244 - t548) * t590 - t84 * t572 - m(7) * (t120 * t572 + t19 * t21 + t20 * t22); t69 * t595 - t19 * t122 + t20 * t123 - g(1) * t506 - g(2) * t505 - g(3) * (-mrSges(7,1) * t519 + t321) + (t664 + t70) * t598 + t678;];
tau  = t1;
