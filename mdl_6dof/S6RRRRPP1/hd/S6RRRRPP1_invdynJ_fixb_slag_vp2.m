% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP1
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
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:44:04
% EndTime: 2019-03-09 20:44:49
% DurationCPUTime: 25.77s
% Computational Cost: add. (17999->803), mult. (39145->1004), div. (0->0), fcn. (28237->14), ass. (0->390)
t413 = sin(qJ(2));
t417 = cos(qJ(2));
t503 = qJD(1) * qJD(2);
t345 = qJDD(1) * t417 - t413 * t503;
t346 = qJDD(1) * t413 + t417 * t503;
t412 = sin(qJ(3));
t416 = cos(qJ(3));
t338 = t412 * t417 + t413 * t416;
t435 = t338 * qJD(3);
t223 = -qJD(1) * t435 + t345 * t416 - t346 * t412;
t221 = qJDD(4) - t223;
t598 = t221 / 0.2e1;
t728 = 0.2e1 * t598;
t337 = t412 * t413 - t416 * t417;
t434 = t337 * qJD(3);
t222 = -qJD(1) * t434 + t345 * t412 + t346 * t416;
t322 = t338 * qJD(1);
t406 = qJD(2) + qJD(3);
t411 = sin(qJ(4));
t415 = cos(qJ(4));
t287 = -t322 * t411 + t406 * t415;
t405 = qJDD(2) + qJDD(3);
t159 = qJD(4) * t287 + t222 * t415 + t405 * t411;
t288 = t322 * t415 + t406 * t411;
t160 = -qJD(4) * t288 - t222 * t411 + t405 * t415;
t409 = sin(pkin(10));
t549 = cos(pkin(10));
t85 = t159 * t549 + t160 * t409;
t609 = t85 / 0.2e1;
t727 = 0.2e1 * t609;
t321 = t337 * qJD(1);
t263 = pkin(3) * t322 + pkin(9) * t321;
t511 = qJD(1) * t413;
t495 = pkin(2) * t511;
t239 = t263 + t495;
t419 = -pkin(8) - pkin(7);
t363 = t419 * t417;
t341 = qJD(1) * t363;
t323 = t412 * t341;
t362 = t419 * t413;
t340 = qJD(1) * t362;
t272 = t340 * t416 + t323;
t174 = t239 * t411 + t272 * t415;
t507 = qJD(3) * t416;
t493 = pkin(2) * t507;
t726 = t415 * t493 - t174;
t546 = t321 * t411;
t725 = qJ(5) * t546 - qJD(5) * t415;
t505 = qJD(4) * t411;
t724 = t505 + t546;
t327 = qJD(2) * pkin(2) + t340;
t267 = t327 * t416 + t323;
t244 = -pkin(3) * t406 - t267;
t184 = -pkin(4) * t287 + qJD(5) + t244;
t314 = qJD(4) + t321;
t404 = t417 * pkin(2);
t392 = t404 + pkin(1);
t361 = t392 * qJD(1);
t232 = pkin(3) * t321 - pkin(9) * t322 - t361;
t520 = t416 * t341;
t268 = t327 * t412 - t520;
t245 = pkin(9) * t406 + t268;
t155 = t232 * t415 - t245 * t411;
t120 = -qJ(5) * t288 + t155;
t109 = pkin(4) * t314 + t120;
t156 = t232 * t411 + t245 * t415;
t121 = qJ(5) * t287 + t156;
t115 = t549 * t121;
t40 = t109 * t409 + t115;
t35 = qJ(6) * t314 + t40;
t590 = t314 / 0.2e1;
t591 = -t314 / 0.2e1;
t437 = t287 * t409 + t288 * t549;
t599 = t437 / 0.2e1;
t600 = -t437 / 0.2e1;
t197 = -t287 * t549 + t288 * t409;
t603 = -t197 / 0.2e1;
t70 = pkin(5) * t197 - qJ(6) * t437 + t184;
t723 = -mrSges(6,1) * t184 - mrSges(7,1) * t70 + Ifges(6,4) * t599 + Ifges(7,5) * t600 + Ifges(6,6) * t590 + Ifges(7,6) * t591 + (Ifges(6,2) + Ifges(7,3)) * t603 + mrSges(7,2) * t35 + mrSges(6,3) * t40;
t532 = t409 * t121;
t39 = t109 * t549 - t532;
t34 = -pkin(5) * t314 + qJD(6) - t39;
t722 = -mrSges(6,2) * t184 - mrSges(7,2) * t34 + mrSges(6,3) * t39 + mrSges(7,3) * t70;
t712 = -mrSges(6,1) - mrSges(7,1);
t623 = m(7) * pkin(5) - t712;
t84 = t159 * t409 - t160 * t549;
t611 = -t84 / 0.2e1;
t720 = mrSges(6,3) + mrSges(7,2);
t665 = Ifges(6,1) + Ifges(7,1);
t664 = Ifges(6,4) - Ifges(7,5);
t663 = Ifges(7,4) + Ifges(6,5);
t662 = Ifges(6,6) - Ifges(7,6);
t173 = t239 * t415 - t272 * t411;
t403 = t415 * qJ(5);
t460 = pkin(4) * t322 + t321 * t403;
t577 = pkin(2) * t412;
t389 = pkin(9) + t577;
t517 = -qJ(5) - t389;
t464 = qJD(4) * t517;
t719 = -t173 - t460 + (-qJD(5) - t493) * t411 + t415 * t464;
t176 = t263 * t415 - t267 * t411;
t410 = -qJ(5) - pkin(9);
t474 = qJD(4) * t410;
t718 = -qJD(5) * t411 + t415 * t474 - t176 - t460;
t717 = -t411 * t464 + t725 - t726;
t177 = t263 * t411 + t267 * t415;
t716 = -t411 * t474 + t177 + t725;
t602 = t197 / 0.2e1;
t715 = -Ifges(6,2) * t603 + Ifges(7,3) * t602 - t590 * t662 - t599 * t664 - t723;
t714 = Ifges(6,4) * t602 + Ifges(7,5) * t603 + t600 * t665 + t722;
t713 = -Ifges(6,2) * t602 + Ifges(7,3) * t603 - t600 * t664 + t723;
t630 = Ifges(6,3) + Ifges(5,3) + Ifges(7,2);
t711 = mrSges(6,2) - mrSges(7,3);
t671 = -m(7) - m(6);
t695 = pkin(4) * t671;
t709 = mrSges(5,1) - t695;
t691 = t724 * pkin(4);
t408 = qJ(2) + qJ(3);
t401 = sin(t408);
t407 = qJ(4) + pkin(10);
t398 = sin(t407);
t540 = t398 * t401;
t563 = mrSges(5,2) * t411;
t708 = -mrSges(6,2) * t540 - t401 * t563;
t657 = t406 * Ifges(4,5);
t707 = -t361 * mrSges(4,2) - t267 * mrSges(4,3) + t657 / 0.2e1;
t656 = t406 * Ifges(4,6);
t706 = t361 * mrSges(4,1) - t155 * mrSges(5,1) - t39 * mrSges(6,1) + t34 * mrSges(7,1) + t156 * mrSges(5,2) + t40 * mrSges(6,2) - t35 * mrSges(7,3) + t656 / 0.2e1;
t705 = -mrSges(5,3) - t720;
t548 = qJDD(1) * pkin(1);
t310 = -pkin(2) * t345 - t548;
t122 = -pkin(3) * t223 - pkin(9) * t222 + t310;
t333 = t346 * pkin(7);
t279 = qJDD(2) * pkin(2) - pkin(8) * t346 - t333;
t332 = t345 * pkin(7);
t286 = pkin(8) * t345 + t332;
t508 = qJD(3) * t412;
t137 = t279 * t412 + t286 * t416 + t327 * t507 + t341 * t508;
t131 = pkin(9) * t405 + t137;
t30 = -qJD(4) * t156 + t122 * t415 - t131 * t411;
t14 = pkin(4) * t221 - qJ(5) * t159 - qJD(5) * t288 + t30;
t504 = qJD(4) * t415;
t29 = t122 * t411 + t131 * t415 + t232 * t504 - t245 * t505;
t16 = qJ(5) * t160 + qJD(5) * t287 + t29;
t6 = t14 * t409 + t16 * t549;
t1 = qJ(6) * t221 + qJD(6) * t314 + t6;
t610 = t84 / 0.2e1;
t138 = t279 * t416 - t286 * t412 - t327 * t508 + t341 * t507;
t132 = -pkin(3) * t405 - t138;
t69 = -pkin(4) * t160 + qJDD(5) + t132;
t9 = pkin(5) * t84 - qJ(6) * t85 - qJD(6) * t437 + t69;
t701 = t662 * t728 + t664 * t727 - mrSges(6,1) * t69 - mrSges(7,1) * t9 + mrSges(7,2) * t1 + mrSges(6,3) * t6 - 0.2e1 * Ifges(7,3) * t610 + (-t610 + t611) * Ifges(6,2);
t700 = Ifges(6,4) * t603 + Ifges(7,5) * t602 + t590 * t663 + t599 * t665 - t722;
t5 = t14 * t549 - t16 * t409;
t3 = -pkin(5) * t221 + qJDD(6) - t5;
t699 = mrSges(6,2) * t69 + mrSges(7,2) * t3 - mrSges(6,3) * t5 - mrSges(7,3) * t9 + Ifges(7,5) * t610 + t663 * t728 + t665 * t727 + (Ifges(6,4) + t664) * t611;
t103 = Ifges(7,1) * t437 + Ifges(7,4) * t314 + Ifges(7,5) * t197;
t104 = Ifges(6,1) * t437 - Ifges(6,4) * t197 + Ifges(6,5) * t314;
t698 = t104 / 0.2e1 + t103 / 0.2e1;
t402 = cos(t408);
t638 = pkin(3) * t402 + pkin(9) * t401;
t696 = m(5) * t638;
t655 = t409 * t717 + t549 * t719;
t652 = t409 * t719 - t549 * t717;
t650 = t409 * t716 + t549 * t718;
t647 = t409 * t718 - t549 * t716;
t473 = t549 * t411;
t335 = t409 * t415 + t473;
t228 = t335 * t321;
t472 = t549 * t415;
t531 = t409 * t411;
t436 = t472 - t531;
t229 = t436 * t321;
t317 = t335 * qJD(4);
t318 = t436 * qJD(4);
t692 = -qJD(6) * t335 + t691 + (-t318 - t229) * qJ(6) + (t228 + t317) * pkin(5);
t552 = t322 * mrSges(4,3);
t641 = mrSges(4,1) * t406 + mrSges(5,1) * t287 - mrSges(5,2) * t288 - t552;
t690 = t720 * t401;
t689 = -t402 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t401;
t271 = t340 * t412 - t520;
t688 = pkin(2) * t508 - t271;
t275 = -qJD(2) * t337 - t434;
t481 = t338 * t504;
t440 = t275 * t411 + t481;
t414 = sin(qJ(1));
t524 = t414 * t415;
t418 = cos(qJ(1));
t526 = t411 * t418;
t308 = -t402 * t526 + t524;
t687 = t29 * t415 - t30 * t411;
t564 = mrSges(5,1) * t415;
t686 = t563 - t564;
t658 = t287 * Ifges(5,6);
t682 = t288 * Ifges(5,5) - t197 * t662 + t314 * t630 + t437 * t663 + t658;
t390 = pkin(4) * t415 + pkin(3);
t399 = cos(t407);
t681 = (t564 - m(7) * (-qJ(6) * t398 - t390) + t398 * mrSges(7,3) + t623 * t399) * t401;
t605 = t159 / 0.2e1;
t604 = t160 / 0.2e1;
t666 = t345 / 0.2e1;
t582 = t417 / 0.2e1;
t45 = -mrSges(7,2) * t84 + mrSges(7,3) * t221;
t46 = -mrSges(6,2) * t221 - mrSges(6,3) * t84;
t661 = t45 + t46;
t47 = mrSges(6,1) * t221 - mrSges(6,3) * t85;
t48 = -mrSges(7,1) * t221 + mrSges(7,2) * t85;
t660 = t48 - t47;
t659 = Ifges(3,2) * t417;
t566 = t322 * pkin(5);
t654 = t566 - t655;
t305 = t322 * qJ(6);
t653 = -t305 + t652;
t651 = -t268 + t692;
t649 = t566 - t650;
t648 = -t305 + t647;
t455 = mrSges(5,1) * t411 + mrSges(5,2) * t415;
t646 = t244 * t455;
t169 = mrSges(6,1) * t314 - mrSges(6,3) * t437;
t170 = -mrSges(7,1) * t314 + mrSges(7,2) * t437;
t644 = t169 - t170;
t643 = t688 + t692;
t642 = t688 + t691;
t266 = pkin(3) * t337 - pkin(9) * t338 - t392;
t290 = t362 * t412 - t363 * t416;
t282 = t415 * t290;
t188 = t266 * t411 + t282;
t640 = t362 * t416 + t363 * t412;
t637 = -t268 + t691;
t42 = t120 * t549 - t532;
t635 = -t42 + qJD(6);
t510 = qJD(1) * t417;
t570 = pkin(7) * t417;
t571 = pkin(7) * t413;
t634 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t511) * t570 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t510) * t571;
t633 = t332 * t417 + t333 * t413;
t632 = g(1) * t418 + g(2) * t414;
t629 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t534 = t402 * t414;
t628 = t414 * t708 + t534 * t705;
t533 = t402 * t418;
t627 = t418 * t708 + t533 * t705;
t535 = t402 * t410;
t448 = -t390 * t401 - t535;
t576 = pkin(2) * t413;
t626 = -m(7) * (-t535 - t576) - m(6) * (t448 - t576) - m(5) * (-pkin(3) * t401 - t576) + t681;
t625 = -m(6) * t448 + m(7) * t535 + t681;
t624 = Ifges(5,5) * t159 + Ifges(5,6) * t160 + t221 * t630 - t662 * t84 + t663 * t85;
t359 = -mrSges(3,1) * t417 + mrSges(3,2) * t413;
t622 = m(3) * pkin(1) + mrSges(2,1) - t359 - t689;
t538 = t399 * t402;
t539 = t398 * t402;
t621 = t402 * t686 + t538 * t712 + t539 * t711 + t689 - t690;
t620 = m(7) * qJ(6) - t711;
t106 = mrSges(5,1) * t221 - mrSges(5,3) * t159;
t107 = -mrSges(5,2) * t221 + mrSges(5,3) * t160;
t562 = mrSges(5,3) * t287;
t225 = -mrSges(5,2) * t314 + t562;
t561 = mrSges(5,3) * t288;
t226 = mrSges(5,1) * t314 - t561;
t618 = m(5) * ((-t155 * t415 - t156 * t411) * qJD(4) + t687) - t226 * t504 - t225 * t505 + t415 * t107 - t411 * t106;
t616 = t30 * mrSges(5,1) + t5 * mrSges(6,1) - t3 * mrSges(7,1) - t29 * mrSges(5,2) - t6 * mrSges(6,2) + t1 * mrSges(7,3);
t612 = Ifges(5,1) * t605 + Ifges(5,4) * t604 + Ifges(5,5) * t598;
t594 = -t287 / 0.2e1;
t593 = -t288 / 0.2e1;
t592 = t288 / 0.2e1;
t587 = -t321 / 0.2e1;
t585 = t322 / 0.2e1;
t575 = pkin(2) * t416;
t574 = pkin(4) * t288;
t573 = pkin(4) * t409;
t276 = qJD(2) * t338 + t435;
t447 = -qJ(5) * t275 - qJD(5) * t338;
t509 = qJD(2) * t413;
t494 = pkin(2) * t509;
t183 = pkin(3) * t276 - pkin(9) * t275 + t494;
t484 = qJD(2) * t419;
t343 = t413 * t484;
t344 = t417 * t484;
t200 = qJD(3) * t640 + t343 * t416 + t344 * t412;
t468 = t183 * t415 - t200 * t411;
t32 = pkin(4) * t276 + t447 * t415 + (-t282 + (qJ(5) * t338 - t266) * t411) * qJD(4) + t468;
t486 = t183 * t411 + t200 * t415 + t266 * t504;
t38 = -qJ(5) * t481 + (-qJD(4) * t290 + t447) * t411 + t486;
t12 = t32 * t409 + t38 * t549;
t560 = Ifges(3,4) * t413;
t559 = Ifges(3,4) * t417;
t558 = Ifges(5,4) * t288;
t557 = Ifges(5,4) * t411;
t556 = Ifges(5,4) * t415;
t551 = t322 * Ifges(4,4);
t543 = t338 * t411;
t542 = t338 * t415;
t537 = t401 * t414;
t536 = t401 * t418;
t354 = t402 * t390;
t179 = Ifges(5,2) * t287 + Ifges(5,6) * t314 + t558;
t528 = t411 * t179;
t527 = t411 * t414;
t525 = t414 * t399;
t283 = Ifges(5,4) * t287;
t180 = Ifges(5,1) * t288 + Ifges(5,5) * t314 + t283;
t522 = t415 * t180;
t521 = t415 * t418;
t519 = t418 * t398;
t518 = t418 * t419;
t187 = t266 * t415 - t290 * t411;
t147 = pkin(4) * t337 - t338 * t403 + t187;
t164 = -qJ(5) * t543 + t188;
t72 = t147 * t409 + t164 * t549;
t506 = qJD(4) * t338;
t483 = t549 * pkin(4);
t482 = t338 * t505;
t478 = t522 / 0.2e1;
t27 = mrSges(6,1) * t84 + mrSges(6,2) * t85;
t26 = mrSges(7,1) * t84 - mrSges(7,3) * t85;
t475 = -t505 / 0.2e1;
t471 = t503 / 0.2e1;
t469 = t517 * t411;
t467 = -t401 * t410 + t354;
t466 = t392 * t418 - t414 * t419;
t233 = pkin(4) * t543 - t640;
t458 = mrSges(3,1) * t413 + mrSges(3,2) * t417;
t456 = mrSges(4,1) * t401 + mrSges(4,2) * t402;
t454 = Ifges(5,1) * t415 - t557;
t453 = t560 + t659;
t452 = -Ifges(5,2) * t411 + t556;
t451 = Ifges(3,5) * t417 - Ifges(3,6) * t413;
t450 = Ifges(5,5) * t415 - Ifges(5,6) * t411;
t443 = pkin(5) * t538 + qJ(6) * t539 + t467;
t441 = pkin(1) * t458;
t306 = t402 * t527 + t521;
t11 = t32 * t549 - t38 * t409;
t439 = -t275 * t415 + t482;
t438 = t413 * (Ifges(3,1) * t417 - t560);
t71 = t147 * t549 - t164 * t409;
t260 = -pkin(5) * t436 - qJ(6) * t335 - t390;
t201 = qJD(3) * t290 + t343 * t412 - t344 * t416;
t129 = pkin(4) * t440 + t201;
t240 = -t321 * Ifges(4,2) + t551 + t656;
t311 = Ifges(4,4) * t321;
t241 = t322 * Ifges(4,1) - t311 + t657;
t53 = Ifges(5,4) * t159 + Ifges(5,2) * t160 + Ifges(5,6) * t221;
t420 = (-t724 * t156 + (-t321 * t415 - t504) * t155 + t687) * mrSges(5,3) + (-t452 * t594 - t454 * t593 + t646 + t707) * t321 + (Ifges(5,5) * t593 + Ifges(5,6) * t594 + Ifges(6,6) * t602 + Ifges(7,6) * t603 + t663 * t600 + t706) * t322 + (t478 + t646) * qJD(4) + (-Ifges(4,2) * t322 + t241 - t311 + t522) * t321 / 0.2e1 + (t228 * t662 - t229 * t663 - t321 * t450 + t322 * t630) * t591 + Ifges(4,6) * t223 + t699 * t335 + t700 * t318 - t713 * t228 - t714 * t229 + t715 * t317 + (t287 * t452 + t288 * t454 + t314 * t450) * qJD(4) / 0.2e1 + t179 * t475 - (-Ifges(4,1) * t321 - t551 + t682) * t322 / 0.2e1 + t132 * t686 + t701 * t436 + Ifges(4,5) * t222 - t137 * mrSges(4,2) + t138 * mrSges(4,1) + t268 * t552 + t240 * t585 + t528 * t587 + (Ifges(5,5) * t411 + Ifges(5,6) * t415) * t598 + (Ifges(5,2) * t415 + t557) * t604 + (Ifges(5,1) * t411 + t556) * t605 + t411 * t612 + Ifges(4,3) * t405 + t415 * t53 / 0.2e1 + (t104 + t103) * (t318 / 0.2e1 + t229 / 0.2e1);
t394 = Ifges(3,4) * t510;
t391 = -pkin(3) - t575;
t385 = -t483 - pkin(5);
t379 = qJ(6) + t573;
t373 = pkin(9) * t533;
t372 = pkin(9) * t534;
t360 = pkin(9) * t415 + t403;
t356 = -t390 - t575;
t328 = t389 * t415 + t403;
t320 = Ifges(3,1) * t511 + Ifges(3,5) * qJD(2) + t394;
t319 = Ifges(3,6) * qJD(2) + qJD(1) * t453;
t309 = t402 * t521 + t527;
t307 = -t402 * t524 + t526;
t296 = t398 * t414 + t399 * t533;
t295 = t402 * t519 - t525;
t294 = t402 * t525 - t519;
t293 = t398 * t534 + t399 * t418;
t291 = -mrSges(4,2) * t406 - mrSges(4,3) * t321;
t285 = t360 * t549 + t410 * t531;
t284 = t360 * t409 - t410 * t473;
t262 = mrSges(4,1) * t321 + mrSges(4,2) * t322;
t257 = t328 * t549 + t409 * t469;
t256 = t328 * t409 - t469 * t549;
t248 = t436 * t338;
t247 = t335 * t338;
t243 = t260 - t575;
t207 = -mrSges(4,2) * t405 + mrSges(4,3) * t223;
t206 = mrSges(4,1) * t405 - mrSges(4,3) * t222;
t168 = -mrSges(6,2) * t314 - mrSges(6,3) * t197;
t167 = -mrSges(7,2) * t197 + mrSges(7,3) * t314;
t127 = t275 * t436 - t335 * t506;
t126 = -t275 * t335 + t409 * t482 - t472 * t506;
t114 = pkin(5) * t247 - qJ(6) * t248 + t233;
t113 = mrSges(6,1) * t197 + mrSges(6,2) * t437;
t112 = mrSges(7,1) * t197 - mrSges(7,3) * t437;
t97 = pkin(5) * t437 + qJ(6) * t197 + t574;
t86 = -mrSges(5,1) * t160 + mrSges(5,2) * t159;
t68 = -pkin(5) * t337 - t71;
t67 = -qJD(4) * t188 + t468;
t66 = -t290 * t505 + t486;
t65 = qJ(6) * t337 + t72;
t41 = t120 * t409 + t115;
t25 = -pkin(5) * t126 - qJ(6) * t127 - qJD(6) * t248 + t129;
t10 = -pkin(5) * t276 - t11;
t8 = qJ(6) * t276 + qJD(6) * t337 + t12;
t2 = [(-t309 * mrSges(5,1) - t308 * mrSges(5,2) - t720 * t536 + (-m(4) - m(5)) * t466 + t671 * (pkin(4) * t527 + t390 * t533 - t410 * t536 + t466) - t623 * t296 - t620 * t295 + t629 * t414 + (-t622 - t696) * t418) * g(2) + (-t528 / 0.2e1 + t241 / 0.2e1 + t478 + Ifges(4,1) * t585 + Ifges(4,4) * t587 + t707) * t275 + (t682 / 0.2e1 + t658 / 0.2e1 - Ifges(4,4) * t585 - Ifges(4,2) * t587 + Ifges(5,5) * t592 + Ifges(7,6) * t602 + Ifges(6,6) * t603 - t240 / 0.2e1 + t663 * t599 + t630 * t590 - t268 * mrSges(4,3) - t706) * t276 + (-Ifges(5,5) * t439 - Ifges(5,6) * t440) * t590 + (t310 * mrSges(4,1) - t137 * mrSges(4,3) - Ifges(4,4) * t222 + Ifges(5,5) * t605 - Ifges(4,2) * t223 - Ifges(4,6) * t405 + Ifges(5,6) * t604 + Ifges(6,6) * t611 + Ifges(7,6) * t610 + t630 * t598 + t663 * t609 + t616 + t624 / 0.2e1) * t337 + t287 * (-Ifges(5,4) * t439 - Ifges(5,2) * t440) / 0.2e1 + (Ifges(3,4) * t346 + Ifges(3,2) * t345) * t582 + m(4) * (t137 * t290 + t200 * t268 - t310 * t392 - t361 * t494) + (t155 * t439 - t156 * t440 - t29 * t543 - t30 * t542) * mrSges(5,3) + m(5) * (t155 * t67 + t156 * t66 + t187 * t30 + t188 * t29) + (mrSges(4,2) * t310 - mrSges(4,3) * t138 + Ifges(4,1) * t222 + Ifges(4,4) * t223 + Ifges(4,5) * t405 + t132 * t455 + t180 * t475 + t450 * t598 + t452 * t604 + t454 * t605) * t338 + t233 * t27 + t66 * t225 + t67 * t226 + t262 * t494 + t346 * t559 / 0.2e1 + t699 * t248 + (t698 + t700) * t127 - t715 * t126 + t200 * t291 + t290 * t207 + (-Ifges(5,1) * t439 - Ifges(5,4) * t440) * t592 - t319 * t509 / 0.2e1 - t441 * t503 + (m(5) * t518 - t307 * mrSges(5,1) - t306 * mrSges(5,2) + t671 * (pkin(4) * t526 + t410 * t537 - t518) + t623 * t294 + t620 * t293 + (m(4) * t419 + t629) * t418 + (-m(5) * (-t392 - t638) + m(4) * t392 + t671 * (-t392 - t354) + t622 + t690) * t414) * g(1) - t179 * t481 / 0.2e1 + m(7) * (t1 * t65 + t10 * t34 + t114 * t9 + t25 * t70 + t3 * t68 + t35 * t8) + m(6) * (t11 * t39 + t12 * t40 + t129 * t184 + t233 * t69 + t5 * t71 + t6 * t72) - t53 * t543 / 0.2e1 - t701 * t247 + t187 * t106 + t188 * t107 + t10 * t170 + t8 * t167 + t12 * t168 + t11 * t169 + (t417 * t559 + t438) * t471 + Ifges(2,3) * qJDD(1) + t129 * t113 + t25 * t112 + t114 * t26 + t71 * t47 + t72 * t46 + t68 * t48 + t65 * t45 + (-mrSges(3,1) * t571 - mrSges(3,2) * t570 + 0.2e1 * Ifges(3,6) * t582) * qJDD(2) + (Ifges(3,1) * t346 + Ifges(3,4) * t666 + Ifges(3,5) * qJDD(2) - t471 * t659) * t413 + t244 * (mrSges(5,1) * t440 - mrSges(5,2) * t439) + (-m(4) * t267 + m(5) * t244 - t641) * t201 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t633) + (t320 * t582 + t451 * qJD(2) / 0.2e1 - t634) * qJD(2) + (t345 * t570 + t346 * t571 + t633) * mrSges(3,3) - pkin(1) * (-mrSges(3,1) * t345 + mrSges(3,2) * t346) + t453 * t666 + t542 * t612 - t392 * (-mrSges(4,1) * t223 + mrSges(4,2) * t222) - t359 * t548 - (-m(4) * t138 + m(5) * t132 - t206 + t86) * t640; t726 * t225 - (-Ifges(3,2) * t511 + t320 + t394) * t510 / 0.2e1 + (m(4) * t576 + t456 + t458) * t632 + ((t137 * t412 + t138 * t416 + (-t267 * t412 + t268 * t416) * qJD(3)) * pkin(2) + t267 * t271 - t268 * t272 + t361 * t495) * m(4) + t243 * t26 + (t634 + (t441 - t438 / 0.2e1) * qJD(1)) * qJD(1) + (-m(4) * t404 - m(6) * (t404 + t467) - m(5) * (t404 + t638) - m(7) * (t404 + t443) + t359 + t621) * g(3) + t319 * t511 / 0.2e1 - t451 * t503 / 0.2e1 - t641 * t688 + (-t411 * t493 - t173) * t226 + t420 + (t493 - t272) * t291 - t262 * t495 + Ifges(3,3) * qJDD(2) + t660 * t256 + t661 * t257 + t652 * t168 + t653 * t167 + t654 * t170 + (t1 * t257 + t243 * t9 + t256 * t3 + t34 * t654 + t35 * t653 + t643 * t70) * m(7) + t655 * t169 + (t642 * t184 - t256 * t5 + t257 * t6 + t356 * t69 + t655 * t39 + t652 * t40) * m(6) + t642 * t113 + t643 * t112 + (-m(5) * t373 + t418 * t626 + t627) * g(1) + (-m(5) * t372 + t414 * t626 + t628) * g(2) + t618 * t389 + (-t155 * t173 - t156 * t174 - t244 * t271 + t132 * t391 + (t244 * t412 + (-t155 * t411 + t156 * t415) * t416) * qJD(3) * pkin(2)) * m(5) - t332 * mrSges(3,2) - t333 * mrSges(3,1) + Ifges(3,6) * t345 + Ifges(3,5) * t346 + t356 * t27 + t206 * t575 + t207 * t577 + t391 * t86; -pkin(3) * t86 - t176 * t226 - t177 * t225 + t260 * t26 - t267 * t291 - t390 * t27 + t420 + t632 * t456 + t661 * t285 + t660 * t284 + t641 * t268 + t649 * t170 + t650 * t169 + t647 * t168 + t648 * t167 + t637 * t113 + t651 * t112 + (t1 * t285 + t260 * t9 + t284 * t3 + t34 * t649 + t35 * t648 + t651 * t70) * m(7) + (t184 * t637 - t284 * t5 + t285 * t6 + t39 * t650 - t390 * t69 + t40 * t647) * m(6) + (-pkin(3) * t132 - t155 * t176 - t156 * t177 - t244 * t268) * m(5) + (-m(5) * (-pkin(3) * t537 + t372) + t625 * t414 + t628) * g(2) + (-m(5) * (-pkin(3) * t536 + t373) + t625 * t418 + t627) * g(1) + (-m(6) * t467 - m(7) * t443 + t621 - t696) * g(3) + t618 * pkin(9); (-t307 * mrSges(5,2) + t293 * t623 - t294 * t620 + t306 * t709) * g(2) + (t309 * mrSges(5,2) + t623 * t295 - t620 * t296 - t308 * t709) * g(1) + t47 * t483 + t624 + (-t591 * t662 + t713) * t437 + (-t591 * t663 + t698 - t714) * t197 - t244 * (mrSges(5,1) * t288 + mrSges(5,2) * t287) + t616 + (-t225 + t562) * t155 - t113 * t574 + (t398 * t623 - t399 * t620 - t411 * t695 + t455) * g(3) * t401 + (-t184 * t574 + t39 * t41 - t40 * t42 + (t409 * t6 + t5 * t549) * pkin(4)) * m(6) - t42 * t168 - t97 * t112 + t644 * t41 + (t1 * t379 + t3 * t385 - t34 * t41 + t35 * t635 - t70 * t97) * m(7) + t635 * t167 + (-Ifges(5,2) * t288 + t180 + t283) * t594 + t46 * t573 + (Ifges(5,5) * t287 - Ifges(5,6) * t288) * t591 + t179 * t592 + (Ifges(5,1) * t287 - t558) * t593 + t379 * t45 + t385 * t48 + (t226 + t561) * t156; -(-t167 - t168) * t197 + t644 * t437 + t26 + t27 + (-g(3) * t402 + t401 * t632) * t671 + (t197 * t35 - t34 * t437 + t9) * m(7) + (t197 * t40 + t39 * t437 + t69) * m(6); t437 * t112 - t314 * t167 + (-g(1) * t295 - g(2) * t293 - g(3) * t540 - t314 * t35 + t437 * t70 + t3) * m(7) + t48;];
tau  = t2;
