% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:40
% EndTime: 2019-03-09 16:13:09
% DurationCPUTime: 60.64s
% Computational Cost: add. (17399->1100), mult. (42618->1469), div. (0->0), fcn. (34030->12), ass. (0->464)
t395 = sin(qJ(3));
t399 = cos(qJ(2));
t392 = sin(pkin(6));
t528 = qJD(1) * t392;
t501 = t399 * t528;
t468 = t395 * t501;
t523 = qJD(3) * t395;
t647 = t468 - t523;
t684 = Ifges(5,1) + Ifges(6,1);
t683 = -Ifges(5,4) + Ifges(6,5);
t682 = Ifges(6,4) + Ifges(5,5);
t681 = Ifges(5,6) - Ifges(6,6);
t680 = Ifges(5,3) + Ifges(6,2);
t553 = cos(pkin(6));
t481 = t553 * qJD(1);
t467 = pkin(1) * t481;
t396 = sin(qJ(2));
t502 = t396 * t528;
t319 = -pkin(8) * t502 + t399 * t467;
t425 = (pkin(2) * t396 - pkin(9) * t399) * t392;
t320 = qJD(1) * t425;
t398 = cos(qJ(3));
t206 = t398 * t319 + t395 * t320;
t184 = qJ(4) * t502 + t206;
t506 = pkin(1) * t553;
t386 = t396 * t506;
t441 = pkin(3) * t395 - qJ(4) * t398;
t536 = t392 * t399;
t215 = (t386 + (pkin(8) + t441) * t536) * qJD(1);
t391 = sin(pkin(11));
t393 = cos(pkin(11));
t115 = t393 * t184 + t391 * t215;
t713 = -qJ(5) * t647 - qJD(5) * t398 - t115;
t333 = qJD(3) * t441 - qJD(4) * t395;
t698 = t391 * t184 + (-t215 + t333) * t393;
t433 = t481 + qJD(2);
t285 = t395 * t433 + t398 * t502;
t517 = qJD(1) * qJD(2);
t327 = (qJDD(1) * t396 + t399 * t517) * t392;
t475 = t553 * qJDD(1);
t376 = t475 + qJDD(2);
t170 = qJD(3) * t285 + t395 * t327 - t398 * t376;
t601 = t170 / 0.2e1;
t284 = t395 * t502 - t398 * t433;
t524 = qJD(3) * t284;
t169 = t398 * t327 + t395 * t376 - t524;
t326 = (-qJDD(1) * t399 + t396 * t517) * t392;
t309 = qJDD(3) + t326;
t131 = t169 * t393 + t309 * t391;
t605 = t131 / 0.2e1;
t712 = t682 * t601 + t684 * t605;
t130 = t169 * t391 - t393 * t309;
t606 = t130 / 0.2e1;
t602 = -t170 / 0.2e1;
t532 = t398 * t399;
t281 = (t391 * t396 + t393 * t532) * t392;
t264 = qJD(1) * t281;
t507 = -pkin(9) * t391 - pkin(4);
t534 = t393 * t398;
t614 = pkin(4) + pkin(5);
t711 = t264 * pkin(10) + t468 * t614 + (-pkin(10) * t534 + (-pkin(5) + t507) * t395) * qJD(3) - t698;
t510 = t392 * t532;
t473 = t391 * t510;
t263 = qJD(1) * t473 - t393 * t502;
t297 = t391 * t333;
t535 = t393 * t395;
t538 = t391 * t398;
t710 = -pkin(10) * t263 + t297 + (-pkin(9) * t535 + pkin(10) * t538) * qJD(3) + t713;
t549 = qJ(5) * t391;
t575 = pkin(4) * t393;
t709 = -t549 - t575;
t363 = qJD(3) - t501;
t437 = t393 * t285 + t363 * t391;
t258 = -pkin(2) * t433 - t319;
t139 = t284 * pkin(3) - t285 * qJ(4) + t258;
t530 = pkin(8) * t536 + t386;
t304 = t553 * pkin(9) + t530;
t259 = qJD(2) * pkin(9) + qJD(1) * t304;
t268 = (-pkin(2) * t399 - pkin(9) * t396 - pkin(1)) * t528;
t162 = t259 * t398 + t268 * t395;
t143 = qJ(4) * t363 + t162;
t73 = t139 * t393 - t391 * t143;
t458 = qJD(5) - t73;
t36 = -pkin(10) * t437 - t284 * t614 + t458;
t394 = sin(qJ(6));
t397 = cos(qJ(6));
t217 = t285 * t391 - t393 * t363;
t74 = t391 * t139 + t393 * t143;
t63 = t284 * qJ(5) + t74;
t44 = pkin(10) * t217 + t63;
t13 = t36 * t397 - t394 * t44;
t14 = t36 * t394 + t397 * t44;
t559 = t162 * mrSges(4,3);
t62 = -pkin(4) * t284 + t458;
t674 = t363 * Ifges(4,6);
t708 = t559 - t258 * mrSges(4,1) - t73 * mrSges(5,1) + t62 * mrSges(6,1) + t13 * mrSges(7,1) + t74 * mrSges(5,2) - t14 * mrSges(7,2) - t63 * mrSges(6,3) + t674 / 0.2e1;
t662 = t683 * t217 + t682 * t284 + t437 * t684;
t707 = t662 / 0.2e1;
t706 = t606 * t683 + t712;
t699 = -pkin(8) * t392 * t517 + pkin(1) * t475;
t516 = qJDD(1) * t392;
t700 = pkin(8) * t516 + qJD(2) * t467;
t221 = t396 * t699 + t399 * t700;
t192 = pkin(9) * t376 + t221;
t513 = pkin(1) * t516;
t411 = pkin(2) * t326 - pkin(9) * t327 - t513;
t522 = qJD(3) * t398;
t70 = t398 * t192 - t259 * t523 + t268 * t522 + t395 * t411;
t53 = qJ(4) * t309 + qJD(4) * t363 + t70;
t222 = -t396 * t700 + t399 * t699;
t193 = -t376 * pkin(2) - t222;
t55 = t170 * pkin(3) - t169 * qJ(4) - t285 * qJD(4) + t193;
t18 = t391 * t55 + t393 * t53;
t12 = t170 * qJ(5) + t284 * qJD(5) + t18;
t705 = t12 * mrSges(6,3);
t704 = t18 * mrSges(5,2);
t703 = m(5) * qJD(4);
t702 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t121 = t217 * t397 - t394 * t437;
t274 = qJD(6) - t284;
t568 = Ifges(4,4) * t285;
t696 = t217 * t394 + t397 * t437;
t673 = Ifges(7,5) * t696 - t284 * Ifges(4,2) + t121 * Ifges(7,6) + t274 * Ifges(7,3) + t568 + t674;
t663 = -t217 * t681 + t284 * t680 + t437 * t682;
t205 = -t395 * t319 + t398 * t320;
t422 = pkin(3) * t502 + t205;
t701 = qJ(5) * t264 - qJD(5) * t535 + t422;
t349 = t391 * t397 - t393 * t394;
t436 = t391 * t394 + t393 * t397;
t697 = mrSges(7,1) * t436 + mrSges(7,2) * t349;
t579 = cos(qJ(1));
t461 = t553 * t579;
t578 = sin(qJ(1));
t340 = t396 * t461 + t399 * t578;
t505 = t392 * t579;
t246 = t340 * t398 - t395 * t505;
t339 = t396 * t578 - t399 * t461;
t178 = t246 * t391 - t339 * t393;
t179 = t246 * t393 + t339 * t391;
t607 = -t130 / 0.2e1;
t695 = Ifges(5,2) * t607 - Ifges(6,3) * t606 + t681 * t601;
t17 = -t391 * t53 + t393 * t55;
t459 = qJDD(5) - t17;
t5 = -pkin(10) * t131 - t170 * t614 + t459;
t9 = pkin(10) * t130 + t12;
t1 = qJD(6) * t13 + t394 * t5 + t397 * t9;
t2 = -qJD(6) * t14 - t394 * t9 + t397 * t5;
t694 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t693 = Ifges(6,5) * t605 + Ifges(6,6) * t601 - t131 * Ifges(5,4) / 0.2e1 + Ifges(5,6) * t602 - mrSges(6,2) * t12 + (Ifges(6,3) + Ifges(5,2)) * t606;
t15 = -pkin(4) * t170 + t459;
t692 = t15 * mrSges(6,2) + t706 + t712;
t690 = m(7) * pkin(5);
t33 = qJD(6) * t121 + t130 * t394 + t131 * t397;
t623 = t33 / 0.2e1;
t34 = -qJD(6) * t696 + t130 * t397 - t131 * t394;
t622 = t34 / 0.2e1;
t689 = -m(7) - m(6);
t160 = qJDD(6) - t170;
t604 = t160 / 0.2e1;
t603 = t169 / 0.2e1;
t585 = t309 / 0.2e1;
t6 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t160;
t688 = t169 * Ifges(4,4) - t170 * Ifges(4,2) + t309 * Ifges(4,6) + t6;
t687 = mrSges(5,2) - mrSges(6,3);
t685 = mrSges(4,3) - mrSges(3,2);
t679 = -t130 * t681 + t131 * t682 + t170 * t680;
t576 = pkin(3) * t398;
t356 = -qJ(4) * t395 - pkin(2) - t576;
t377 = pkin(9) * t538;
t389 = t398 * pkin(4);
t223 = pkin(5) * t398 + t377 + t389 + (-pkin(10) * t395 - t356) * t393;
t283 = pkin(9) * t534 + t391 * t356;
t260 = -qJ(5) * t398 + t283;
t539 = t391 * t395;
t227 = pkin(10) * t539 + t260;
t136 = t223 * t397 - t227 * t394;
t677 = qJD(6) * t136 + t394 * t711 + t397 * t710;
t137 = t223 * t394 + t227 * t397;
t676 = -qJD(6) * t137 - t394 * t710 + t397 * t711;
t675 = t363 * Ifges(4,5);
t572 = -pkin(10) + qJ(4);
t357 = t572 * t391;
t358 = t572 * t393;
t253 = t357 * t394 + t358 * t397;
t161 = -t395 * t259 + t398 * t268;
t150 = t391 * t161;
t194 = pkin(3) * t285 + qJ(4) * t284;
t58 = t150 + (pkin(10) * t284 - t194) * t393 - t614 * t285;
t546 = t284 * t391;
t97 = t393 * t161 + t391 * t194;
t84 = t285 * qJ(5) + t97;
t69 = -pkin(10) * t546 + t84;
t672 = qJD(4) * t349 - qJD(6) * t253 + t394 * t69 - t397 * t58;
t252 = t357 * t397 - t358 * t394;
t671 = qJD(4) * t436 + qJD(6) * t252 - t394 * t58 - t397 * t69;
t456 = mrSges(4,1) * t398 - mrSges(4,2) * t395;
t670 = t456 + mrSges(3,1);
t519 = qJD(5) * t391;
t548 = qJ(5) * t393;
t637 = t391 * t614 - t548;
t669 = -t284 * t637 + t162 + t519;
t420 = -pkin(9) - t637;
t668 = t263 * t614 + t420 * t522 - t701;
t515 = pkin(9) * t523;
t239 = -t393 * t515 + t297;
t667 = t239 + t713;
t664 = pkin(4) * t468 + t507 * t523 - t698;
t440 = -pkin(4) * t391 + t548;
t427 = pkin(9) - t440;
t661 = -pkin(4) * t263 + t427 * t522 + t701;
t165 = t263 * t397 - t264 * t394;
t334 = t436 * qJD(6);
t208 = -t334 * t395 + t349 * t522;
t660 = t165 - t208;
t166 = t263 * t394 + t264 * t397;
t314 = t349 * t395;
t207 = qJD(6) * t314 + t436 * t522;
t659 = t166 - t207;
t167 = t349 * t284;
t335 = t349 * qJD(6);
t658 = t167 - t335;
t168 = t436 * t284;
t657 = t168 - t334;
t245 = t340 * t395 + t398 * t505;
t656 = t709 * t245;
t460 = t553 * t578;
t342 = -t396 * t460 + t399 * t579;
t504 = t392 * t578;
t249 = t342 * t395 - t398 * t504;
t655 = t709 * t249;
t654 = t391 * t515 + t698;
t653 = t239 - t115;
t537 = t392 * t396;
t336 = t395 * t537 - t398 * t553;
t325 = t336 * qJ(5);
t652 = -t391 * t325 - t336 * t575;
t344 = -pkin(8) * t537 + t399 * t506;
t650 = -t391 * t681 + t393 * t682;
t562 = Ifges(6,5) * t391;
t565 = Ifges(5,4) * t391;
t649 = -t393 * t684 - t562 + t565;
t428 = pkin(3) * t363 - qJD(4) + t161;
t71 = -qJD(3) * t162 - t192 * t395 + t398 * t411;
t60 = -t309 * pkin(3) + qJDD(4) - t71;
t644 = t395 * t60 - t428 * t522;
t643 = -t395 * t71 + t398 * t70;
t642 = -t17 * t391 + t18 * t393;
t570 = Ifges(3,4) * t396;
t627 = t392 ^ 2;
t641 = (t396 * (Ifges(3,1) * t399 - t570) / 0.2e1 - pkin(1) * (mrSges(3,1) * t396 + mrSges(3,2) * t399)) * t627;
t640 = -m(7) * t572 + t702;
t453 = -mrSges(6,1) * t393 - mrSges(6,3) * t391;
t455 = -mrSges(5,1) * t393 + mrSges(5,2) * t391;
t639 = t393 * t690 - t453 - t455 + t697;
t638 = mrSges(4,2) + t640;
t636 = -mrSges(5,1) - mrSges(6,1) - t690;
t635 = mrSges(4,1) + t639;
t472 = mrSges(3,3) * t502;
t634 = -m(4) * t258 + mrSges(3,1) * t433 - mrSges(4,1) * t284 - mrSges(4,2) * t285 - t472;
t452 = mrSges(6,1) * t391 - mrSges(6,3) * t393;
t109 = Ifges(5,4) * t437 - t217 * Ifges(5,2) + t284 * Ifges(5,6);
t612 = -t109 / 0.2e1;
t412 = qJ(5) * t437 + t428;
t72 = pkin(4) * t217 - t412;
t633 = -t161 * mrSges(4,3) + t391 * t612 + t72 * t452;
t632 = m(7) * pkin(10) + t702;
t630 = -mrSges(7,1) * t394 - mrSges(7,2) * t397 + t687;
t629 = -mrSges(5,3) * t18 + t693;
t628 = -mrSges(7,1) * t397 + mrSges(7,2) * t394 + t636;
t625 = Ifges(7,4) * t623 + Ifges(7,2) * t622 + Ifges(7,6) * t604;
t624 = Ifges(7,1) * t623 + Ifges(7,4) * t622 + Ifges(7,5) * t604;
t563 = Ifges(7,4) * t696;
t48 = t121 * Ifges(7,2) + t274 * Ifges(7,6) + t563;
t619 = -t48 / 0.2e1;
t618 = t48 / 0.2e1;
t120 = Ifges(7,4) * t121;
t49 = Ifges(7,1) * t696 + t274 * Ifges(7,5) + t120;
t617 = -t49 / 0.2e1;
t616 = t49 / 0.2e1;
t615 = Ifges(4,1) * t603 + Ifges(4,4) * t602 + Ifges(4,5) * t585;
t106 = Ifges(6,5) * t437 + t284 * Ifges(6,6) + t217 * Ifges(6,3);
t613 = t106 / 0.2e1;
t611 = -t121 / 0.2e1;
t610 = t121 / 0.2e1;
t609 = -t696 / 0.2e1;
t608 = t696 / 0.2e1;
t599 = -t217 / 0.2e1;
t598 = t217 / 0.2e1;
t597 = -t437 / 0.2e1;
t596 = t437 / 0.2e1;
t591 = -t274 / 0.2e1;
t590 = t274 / 0.2e1;
t589 = -t284 / 0.2e1;
t588 = t284 / 0.2e1;
t587 = -t285 / 0.2e1;
t586 = t285 / 0.2e1;
t577 = pkin(1) * t392;
t569 = Ifges(3,4) * t399;
t567 = Ifges(4,4) * t395;
t566 = Ifges(4,4) * t398;
t564 = Ifges(5,4) * t393;
t561 = Ifges(6,5) * t393;
t552 = qJ(4) * t249;
t547 = t245 * qJ(4);
t545 = t284 * t393;
t541 = t339 * t395;
t341 = t396 * t579 + t399 * t460;
t540 = t341 * t395;
t533 = t395 * t399;
t531 = pkin(2) * t536 + pkin(9) * t537;
t306 = -t531 - t577;
t321 = qJD(2) * t425;
t323 = t344 * qJD(2);
t118 = -t304 * t523 + t306 * t522 + t395 * t321 + t398 * t323;
t527 = qJD(2) * t396;
t105 = (qJ(4) * t527 - qJD(4) * t399) * t392 + t118;
t337 = t395 * t553 + t398 * t537;
t526 = qJD(2) * t399;
t499 = t392 * t526;
t243 = qJD(3) * t337 + t395 * t499;
t244 = -qJD(3) * t336 + t398 * t499;
t324 = t530 * qJD(2);
t117 = t243 * pkin(3) - t244 * qJ(4) - t337 * qJD(4) + t324;
t46 = t393 * t105 + t391 * t117;
t303 = -pkin(2) * t553 - t344;
t328 = t336 * pkin(3);
t478 = t337 * qJ(4) - t328;
t173 = t303 - t478;
t197 = t398 * t304 + t395 * t306;
t174 = -qJ(4) * t536 + t197;
t91 = t391 * t173 + t393 * t174;
t196 = -t395 * t304 + t398 * t306;
t529 = t579 * pkin(1) + pkin(8) * t504;
t525 = qJD(3) * t217;
t521 = qJD(4) * t391;
t520 = qJD(4) * t393;
t514 = pkin(9) * t522;
t511 = t392 * t533;
t509 = Ifges(4,5) * t169 - Ifges(4,6) * t170 + Ifges(4,3) * t309;
t87 = t325 + t91;
t175 = pkin(3) * t536 - t196;
t508 = Ifges(3,5) * t327 - Ifges(3,6) * t326 + Ifges(3,3) * t376;
t500 = t392 * t527;
t494 = t537 / 0.2e1;
t492 = -t528 / 0.2e1;
t485 = t522 / 0.2e1;
t484 = pkin(3) + t549;
t483 = -t339 * pkin(2) + pkin(9) * t340;
t482 = -t341 * pkin(2) + pkin(9) * t342;
t65 = t130 * mrSges(5,1) + t131 * mrSges(5,2);
t80 = -t170 * mrSges(6,1) + t131 * mrSges(6,2);
t64 = t130 * mrSges(6,1) - t131 * mrSges(6,3);
t45 = -t391 * t105 + t117 * t393;
t234 = t245 * pkin(3);
t480 = qJ(4) * t246 - t234;
t236 = t249 * pkin(3);
t250 = t342 * t398 + t395 * t504;
t479 = qJ(4) * t250 - t236;
t96 = t194 * t393 - t150;
t90 = t173 * t393 - t391 * t174;
t282 = t356 * t393 - t377;
t471 = mrSges(3,3) * t501;
t32 = t243 * qJ(5) + t336 * qJD(5) + t46;
t119 = -t304 * t522 - t306 * t523 + t398 * t321 - t395 * t323;
t470 = pkin(3) * t510 + qJ(4) * t511 + t531;
t466 = t399 * t492;
t462 = -pkin(1) * t578 + pkin(8) * t505;
t10 = -t34 * mrSges(7,1) + t33 * mrSges(7,2);
t457 = mrSges(4,1) * t336 + mrSges(4,2) * t337;
t454 = mrSges(5,1) * t391 + mrSges(5,2) * t393;
t241 = t337 * t391 + t393 * t536;
t242 = t337 * t393 - t391 * t536;
t144 = t241 * t397 - t242 * t394;
t145 = t241 * t394 + t242 * t397;
t451 = mrSges(7,1) * t144 - mrSges(7,2) * t145;
t450 = Ifges(4,1) * t398 - t567;
t447 = -Ifges(4,2) * t395 + t566;
t446 = -Ifges(5,2) * t391 + t564;
t444 = Ifges(4,5) * t398 - Ifges(4,6) * t395;
t442 = Ifges(6,3) * t391 + t561;
t59 = -pkin(10) * t242 - t336 * t614 - t90;
t66 = pkin(10) * t241 + t87;
t19 = -t394 * t66 + t397 * t59;
t20 = t394 * t59 + t397 * t66;
t94 = -mrSges(7,2) * t274 + mrSges(7,3) * t121;
t95 = mrSges(7,1) * t274 - mrSges(7,3) * t696;
t438 = -t394 * t95 + t397 * t94;
t432 = -qJ(4) * t541 - t339 * t576 + t483;
t431 = -qJ(4) * t540 - t341 * t576 + t482;
t430 = t342 * pkin(2) + pkin(9) * t341 + t529;
t426 = qJ(5) * t242 - t175;
t424 = t250 * pkin(3) + t430;
t423 = t258 * (mrSges(4,1) * t395 + mrSges(4,2) * t398);
t182 = t250 * t391 - t341 * t393;
t419 = -g(1) * t182 - g(2) * t178 - g(3) * t241;
t418 = -g(1) * t249 - g(2) * t245 - g(3) * t336;
t414 = -t340 * pkin(2) - t339 * pkin(9) + t462;
t410 = pkin(3) * t500 + t119;
t409 = -pkin(3) * t246 + t414;
t183 = t250 * t393 + t341 * t391;
t405 = t183 * pkin(4) + qJ(5) * t182 + t424;
t404 = Ifges(3,6) * t553 + (t399 * Ifges(3,2) + t570) * t392;
t403 = t392 * t433 * (Ifges(3,5) * t399 - Ifges(3,6) * t396);
t402 = -pkin(4) * t179 - qJ(5) * t178 + t409;
t202 = t244 * t393 + t391 * t500;
t401 = qJ(5) * t202 + qJD(5) * t242 + t410;
t16 = t130 * pkin(4) - t131 * qJ(5) - qJD(5) * t437 + t60;
t374 = Ifges(3,4) * t501;
t351 = -t484 - t575;
t343 = (-mrSges(3,1) * t399 + mrSges(3,2) * t396) * t392;
t338 = t393 * t614 + t484;
t322 = t530 * qJD(1);
t317 = -mrSges(3,2) * t433 + t471;
t315 = t436 * t395;
t305 = t427 * t395;
t280 = -t393 * t537 + t473;
t273 = Ifges(4,4) * t284;
t261 = -t282 + t389;
t257 = t420 * t395;
t255 = Ifges(3,1) * t502 + Ifges(3,5) * t433 + t374;
t254 = Ifges(3,6) * qJD(2) + qJD(1) * t404;
t226 = mrSges(4,1) * t363 - mrSges(4,3) * t285;
t225 = -mrSges(4,2) * t363 - mrSges(4,3) * t284;
t212 = -t341 * t534 + t342 * t391;
t211 = -t341 * t538 - t342 * t393;
t210 = -t339 * t534 + t340 * t391;
t209 = -t339 * t538 - t340 * t393;
t201 = t244 * t391 - t393 * t500;
t155 = t285 * Ifges(4,1) - t273 + t675;
t153 = t285 * Ifges(4,5) - t284 * Ifges(4,6) + t363 * Ifges(4,3);
t149 = -mrSges(6,1) * t284 + mrSges(6,2) * t437;
t148 = mrSges(5,1) * t284 - mrSges(5,3) * t437;
t147 = -mrSges(5,2) * t284 - mrSges(5,3) * t217;
t146 = -mrSges(6,2) * t217 + mrSges(6,3) * t284;
t135 = -mrSges(4,2) * t309 - mrSges(4,3) * t170;
t134 = mrSges(4,1) * t309 - mrSges(4,3) * t169;
t133 = mrSges(5,1) * t217 + mrSges(5,2) * t437;
t132 = mrSges(6,1) * t217 - mrSges(6,3) * t437;
t104 = t284 * t440 + t162;
t103 = t182 * t394 + t183 * t397;
t102 = t182 * t397 - t183 * t394;
t93 = pkin(4) * t241 - t426;
t92 = mrSges(4,1) * t170 + mrSges(4,2) * t169;
t88 = -pkin(4) * t336 - t90;
t86 = -pkin(4) * t285 - t96;
t79 = mrSges(5,1) * t170 - mrSges(5,3) * t131;
t78 = -mrSges(5,2) * t170 - mrSges(5,3) * t130;
t77 = -mrSges(6,2) * t130 + mrSges(6,3) * t170;
t76 = -t241 * t614 + t426;
t68 = qJD(6) * t144 + t201 * t394 + t202 * t397;
t67 = -qJD(6) * t145 + t201 * t397 - t202 * t394;
t61 = -mrSges(7,1) * t121 + mrSges(7,2) * t696;
t52 = -t217 * t614 + t412;
t43 = pkin(4) * t201 - t401;
t35 = -pkin(4) * t243 - t45;
t29 = -t201 * t614 + t401;
t26 = pkin(10) * t201 + t32;
t25 = -pkin(10) * t202 - t243 * t614 - t45;
t24 = -mrSges(7,2) * t160 + mrSges(7,3) * t34;
t23 = mrSges(7,1) * t160 - mrSges(7,3) * t33;
t11 = pkin(5) * t130 + t16;
t4 = -qJD(6) * t20 + t25 * t397 - t26 * t394;
t3 = qJD(6) * t19 + t25 * t394 + t26 * t397;
t7 = [(-m(3) * t319 - t634) * t324 + m(7) * (t1 * t20 - t11 * t76 + t13 * t4 + t14 * t3 + t19 * t2 + t29 * t52) + m(6) * (t12 * t87 + t15 * t88 + t16 * t93 + t32 * t63 + t35 * t62 + t43 * t72) + (Ifges(4,1) * t244 + Ifges(4,5) * t500) * t586 + (Ifges(4,1) * t337 - Ifges(4,5) * t536) * t603 + (Ifges(7,5) * t68 + Ifges(7,6) * t67) * t590 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t604 + (-m(7) * t405 - t103 * mrSges(7,1) - t102 * mrSges(7,2) - mrSges(2,1) * t579 + mrSges(2,2) * t578 - m(4) * t430 - t250 * mrSges(4,1) - m(5) * (t424 + t552) - m(6) * (t405 + t552) - m(3) * t529 - t342 * mrSges(3,1) - mrSges(3,3) * t504 - t685 * t341 + t636 * t183 + t687 * t182 + t638 * t249) * g(2) + (t403 / 0.2e1 + t153 * t494) * qJD(2) + t641 * t517 + (-t428 * mrSges(5,2) + t62 * mrSges(6,2) - t73 * mrSges(5,3) - t72 * mrSges(6,3) + Ifges(5,4) * t599 + Ifges(6,5) * t598 + t682 * t588 + t684 * t596 + t707) * t202 + (Ifges(4,4) * t244 + Ifges(4,6) * t500) * t589 + (Ifges(4,4) * t337 - Ifges(4,6) * t536) * t602 + (Ifges(7,4) * t68 + Ifges(7,2) * t67) * t610 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t622 + (Ifges(7,1) * t68 + Ifges(7,4) * t67) * t608 + (Ifges(7,1) * t145 + Ifges(7,4) * t144) * t623 + t193 * t457 + t11 * t451 + t363 * (Ifges(4,5) * t244 + Ifges(4,3) * t500) / 0.2e1 + (Ifges(4,5) * t337 - Ifges(4,3) * t536) * t585 + (t663 / 0.2e1 - Ifges(4,4) * t586 - Ifges(4,2) * t589 - Ifges(7,3) * t590 - Ifges(7,5) * t608 - Ifges(7,6) * t610 + t680 * t588 + t682 * t596 + Ifges(6,6) * t598 + Ifges(5,6) * t599 - t673 / 0.2e1 - t708) * t243 + (t705 + t680 * t601 + t682 * t605 + Ifges(6,6) * t606 + Ifges(5,6) * t607 + t679 / 0.2e1 - t688 / 0.2e1 - t70 * mrSges(4,3) - t704 + t17 * mrSges(5,1) - t15 * mrSges(6,1) - Ifges(4,6) * t585 - Ifges(4,2) * t602 - Ifges(4,4) * t603 - Ifges(7,3) * t604 - Ifges(7,6) * t622 - Ifges(7,5) * t623 - t694) * t336 + (t60 * mrSges(5,2) - t17 * mrSges(5,3) - t16 * mrSges(6,3) + Ifges(5,4) * t607 + Ifges(6,5) * t606 + t692) * t242 + t76 * t10 + t376 * (Ifges(3,3) * t553 + (Ifges(3,5) * t396 + Ifges(3,6) * t399) * t392) / 0.2e1 + t553 * t508 / 0.2e1 - t509 * t536 / 0.2e1 - t343 * t513 - t410 * t133 - t326 * t404 / 0.2e1 + m(4) * (t118 * t162 + t119 * t161 + t193 * t303 + t196 * t71 + t197 * t70) + (t1 * t144 - t13 * t68 + t14 * t67 - t145 * t2) * mrSges(7,3) + t327 * (Ifges(3,5) * t553 + (t396 * Ifges(3,1) + t569) * t392) / 0.2e1 + t90 * t79 + t91 * t78 + t87 * t77 + t88 * t80 + (-mrSges(5,1) * t428 + mrSges(6,1) * t72 - mrSges(6,2) * t63 - mrSges(5,3) * t74 - Ifges(5,2) * t599 + Ifges(6,3) * t598 - t588 * t681 + t596 * t683 + t612 + t613) * t201 + m(5) * (t17 * t90 + t175 * t60 + t18 * t91 + t410 * t428 + t45 * t73 + t46 * t74) + t52 * (-mrSges(7,1) * t67 + mrSges(7,2) * t68) + (-m(7) * t402 - m(4) * t414 + t246 * mrSges(4,1) + mrSges(2,1) * t578 + mrSges(2,2) * t579 - m(6) * (t402 - t547) - m(5) * (t409 - t547) - m(3) * t462 + t340 * mrSges(3,1) - mrSges(3,3) * t505 + t685 * t339 - t628 * t179 - t630 * t178 - t638 * t245) * g(1) + t71 * (-mrSges(4,1) * t536 - t337 * mrSges(4,3)) + t29 * t61 - t254 * t500 / 0.2e1 + t161 * (mrSges(4,1) * t500 - mrSges(4,3) * t244) + t20 * t24 + Ifges(2,3) * qJDD(1) + t19 * t23 + (t221 * t536 - t222 * t537 - t319 * t499 - t322 * t500 - t326 * t530 - t327 * t344) * mrSges(3,3) + (Ifges(3,1) * t327 - Ifges(3,4) * t326 + Ifges(3,5) * t376) * t494 + (t222 * t553 - t326 * t577 + t344 * t376) * mrSges(3,1) + (Ifges(3,4) * t327 - Ifges(3,2) * t326 + Ifges(3,6) * t376) * t536 / 0.2e1 + (t627 * qJD(1) * (-Ifges(3,2) * t396 + t569) + t392 * t255) * t526 / 0.2e1 + (-t221 * t553 - t327 * t577 - t376 * t530) * mrSges(3,2) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t627 + t221 * t530 + t222 * t344 + t322 * t323) + t323 * t317 + t303 * t92 + (-t162 * t500 + t244 * t258 + t536 * t70) * mrSges(4,2) + t244 * t155 / 0.2e1 + (mrSges(5,1) * t60 + mrSges(6,1) * t16 + t605 * t683 + t629 - t695) * t241 + t93 * t64 + t3 * t94 + t4 * t95 + t43 * t132 + t32 * t146 + t46 * t147 + t45 * t148 + t35 * t149 + t175 * t65 + t196 * t134 + t197 * t135 + t118 * t225 + t119 * t226 + t337 * t615 + t68 * t616 + t67 * t618 + t145 * t624 + t144 * t625; (t472 + t634) * t322 + (-t263 * t681 + t264 * t682 + t468 * t680) * t589 + (-t395 * t649 - t398 * t682) * t605 + (t263 * t683 + t264 * t684 + t468 * t682) * t597 + t662 * (t393 * t485 - t264 / 0.2e1) + t643 * mrSges(4,3) + t676 * t95 + (t1 * t137 - t11 * t257 + t13 * t676 + t136 * t2 + t14 * t677 + t52 * t668) * m(7) + t677 * t94 + t398 * t704 + t535 * t706 + (-t514 - t205) * t226 + (t62 * (-mrSges(6,1) * t395 + mrSges(6,2) * t534) + t74 * (-mrSges(5,2) * t395 - mrSges(5,3) * t538) + t63 * (-mrSges(6,2) * t538 + mrSges(6,3) * t395) + t73 * (mrSges(5,1) * t395 - mrSges(5,3) * t534) + t423) * qJD(3) + (-t515 - t206) * t225 + (-Ifges(3,2) * t502 + t398 * t155 + t255 + t374) * t466 + t644 * t454 + (mrSges(7,1) * t647 + mrSges(7,3) * t659) * t13 + (mrSges(7,1) * t660 - mrSges(7,2) * t659) * t52 + (-mrSges(7,2) * t647 - mrSges(7,3) * t660) * t14 + t661 * t132 + t664 * t149 + t663 * (t395 * t466 + t523 / 0.2e1) + (-t403 / 0.2e1 - t641 * qJD(1)) * qJD(1) + t673 * (-t523 / 0.2e1 + t468 / 0.2e1) + (-t161 * (mrSges(4,1) * t396 - mrSges(4,3) * t532) - t162 * (-mrSges(4,2) * t396 - mrSges(4,3) * t533)) * t528 + (-t317 + t471) * t319 + (t391 * t485 - t263 / 0.2e1) * t106 - t193 * t456 + (t12 * t260 + t15 * t261 + t16 * t305 + t62 * t664 + t63 * t667 + t661 * t72) * m(6) + t667 * t146 + t668 * t61 + t16 * t452 * t395 - t679 * t398 / 0.2e1 + (t395 * t680 + t398 * t650) * t524 / 0.2e1 + (t395 * t650 - t398 * t680) * t601 + t633 * t522 - t398 * t705 + (-pkin(2) * t193 - t161 * t205 - t162 * t206) * m(4) + t653 * t147 + t654 * t148 + t17 * (-mrSges(5,1) * t398 - mrSges(5,3) * t535) + t15 * (mrSges(6,1) * t398 + mrSges(6,2) * t535) + (Ifges(6,6) * t395 + t398 * t442) * t525 / 0.2e1 - (Ifges(5,6) * t395 + t398 * t446) * t525 / 0.2e1 + t688 * t398 / 0.2e1 + (t363 * (Ifges(4,3) * t396 + t399 * t444) + t285 * (Ifges(4,5) * t396 + t399 * t450) + t396 * t153) * t492 + ((t395 * t682 - t398 * t649) * t437 + t285 * t450 + t363 * t444) * qJD(3) / 0.2e1 + t428 * (mrSges(5,1) * t263 + mrSges(5,2) * t264) + t155 * t485 + (-m(4) * t482 - m(5) * t431 + t689 * (t212 * pkin(4) + qJ(5) * t211 + t431) - t685 * t342 + t670 * t341 + t628 * t212 + t630 * t211 - t632 * t540) * g(1) + (-m(4) * t483 - m(5) * t432 + t689 * (t210 * pkin(4) + qJ(5) * t209 + t432) - t685 * t340 + t670 * t339 + t628 * t210 + t630 * t209 - t632 * t541) * g(2) + (-m(4) * t531 - (t396 * mrSges(4,3) + t399 * t456) * t392 - m(5) * t470 + t343 + t689 * (t281 * pkin(4) + t280 * qJ(5) + t470) + t628 * t281 + t630 * t280 + t632 * t511) * g(3) + ((-t134 + t65) * t395 + t644 * m(5) + t398 * t135 + ((-t161 * t398 - t162 * t395) * qJD(3) + t643) * m(4)) * pkin(9) + (t17 * t282 + t18 * t283 - t422 * t428 + t653 * t74 + t654 * t73) * m(5) + (t514 + t422) * t133 + t508 - t523 * t559 - t447 * t524 / 0.2e1 - t423 * t501 + (t284 * (Ifges(4,6) * t396 + t399 * t447) + t396 * t254) * t528 / 0.2e1 + t1 * (-mrSges(7,2) * t398 + mrSges(7,3) * t314) + t2 * (mrSges(7,1) * t398 - mrSges(7,3) * t315) - t11 * (-mrSges(7,1) * t314 + mrSges(7,2) * t315) + t305 * t64 + t283 * t78 + t282 * t79 - t63 * (-t263 * mrSges(6,2) + mrSges(6,3) * t468) - t73 * (mrSges(5,1) * t468 - t264 * mrSges(5,3)) + t263 * t109 / 0.2e1 - t74 * (-mrSges(5,2) * t468 - t263 * mrSges(5,3)) - t72 * (mrSges(6,1) * t263 - mrSges(6,3) * t264) - t62 * (-mrSges(6,1) * t468 + t264 * mrSges(6,2)) + t260 * t77 + t261 * t80 + t257 * t10 - pkin(2) * t92 + t136 * t23 + t137 * t24 - t221 * mrSges(3,2) + t222 * mrSges(3,1) + t629 * t539 + (Ifges(4,5) * t395 + Ifges(4,6) * t398) * t585 + (Ifges(7,5) * t207 + Ifges(7,6) * t208 - Ifges(7,3) * t523) * t590 + (Ifges(7,5) * t166 + Ifges(7,6) * t165 - Ifges(7,3) * t468) * t591 + (Ifges(5,4) * t264 - Ifges(5,2) * t263 + Ifges(5,6) * t468) * t598 + (Ifges(6,5) * t264 + Ifges(6,6) * t468 + Ifges(6,3) * t263) * t599 + (Ifges(4,2) * t398 + t567) * t602 + (Ifges(4,1) * t395 + t566) * t603 + (Ifges(7,5) * t315 + Ifges(7,6) * t314 + Ifges(7,3) * t398) * t604 + (-Ifges(6,6) * t398 + t395 * t442) * t606 + (-Ifges(5,6) * t398 + t395 * t446) * t607 + (Ifges(7,1) * t207 + Ifges(7,4) * t208 - Ifges(7,5) * t523) * t608 + (Ifges(7,1) * t166 + Ifges(7,4) * t165 - Ifges(7,5) * t468) * t609 + (Ifges(7,4) * t207 + Ifges(7,2) * t208 - Ifges(7,6) * t523) * t610 + (Ifges(7,4) * t166 + Ifges(7,2) * t165 - Ifges(7,6) * t468) * t611 + t395 * t615 + t207 * t616 + t166 * t617 + t208 * t618 + t165 * t619 + (Ifges(7,4) * t315 + Ifges(7,2) * t314 + Ifges(7,6) * t398) * t622 + (Ifges(7,1) * t315 + Ifges(7,4) * t314 + Ifges(7,5) * t398) * t623 + t315 * t624 + t314 * t625; t562 * t606 + t565 * t607 + (-t545 * t73 - t546 * t74 + t642) * mrSges(5,3) + (-t521 - t96) * t148 + t545 * t707 + (-t104 * t72 + t16 * t351 - t62 * t86 - t63 * t84) * m(6) - t11 * t697 + (-t568 + t663) * t587 + (-t86 + t521) * t149 + (-t73 * t703 + (qJ(4) * t15 + qJD(4) * t62 - qJD(5) * t72) * m(6) + (-t79 + t80) * qJ(4) + t692) * t391 + (t74 * t703 + (qJ(4) * t12 + qJD(4) * t63) * m(6) + (t78 + t77) * qJ(4) - t693 + t695) * t393 + t16 * t453 + t60 * t455 + t669 * t61 + t671 * t94 + (t1 * t253 - t11 * t338 + t13 * t672 + t14 * t671 + t2 * t252 + t52 * t669) * m(7) + t672 * t95 + t673 * t586 + (-t97 + t520) * t147 + (-Ifges(7,5) * t609 - Ifges(4,2) * t588 + Ifges(5,6) * t598 + Ifges(6,6) * t599 - Ifges(7,6) * t611 - Ifges(7,3) * t591 + t680 * t589 + t682 * t597 + t708) * t285 + (-t133 + t226) * t162 + (t545 * t62 - t546 * t63) * mrSges(6,2) + (t457 - m(7) * (-t328 + t652) - m(6) * (t478 + t652) - m(5) * t478 + t640 * t337 + t639 * t336) * g(3) + (-m(7) * (-t236 + t655) - m(6) * (t479 + t655) - m(5) * t479 + t638 * t250 + t635 * t249) * g(1) + (-m(7) * (-t234 + t656) - m(6) * (t480 + t656) - m(5) * t480 + t638 * t246 + t635 * t245) * g(2) + (-mrSges(7,1) * t658 + mrSges(7,2) * t657) * t52 - t70 * mrSges(4,2) + t71 * mrSges(4,1) + (-t428 * t454 + t675 / 0.2e1 + t258 * mrSges(4,2) - Ifges(4,1) * t587 - t446 * t598 - t442 * t599 + t649 * t597 - t650 * t589 + t633) * t284 + (-pkin(3) * t60 + t642 * qJ(4) + t162 * t428 - t73 * t96 - t74 * t97) * m(5) + (-t1 * t436 - t13 * t657 + t14 * t658 - t2 * t349) * mrSges(7,3) + (Ifges(7,5) * t349 - Ifges(7,6) * t436) * t604 + (Ifges(7,4) * t349 - Ifges(7,2) * t436) * t622 + (Ifges(7,1) * t349 - Ifges(7,4) * t436) * t623 - t436 * t625 + (-t273 + t155) * t588 + t509 + (-t561 + t564) * t605 - pkin(3) * t65 + (-t519 - t104) * t132 + (-Ifges(7,5) * t334 - Ifges(7,6) * t335) * t590 + (-Ifges(7,1) * t334 - Ifges(7,4) * t335) * t608 + (-Ifges(7,4) * t334 - Ifges(7,2) * t335) * t610 + (-Ifges(7,1) * t168 - Ifges(7,4) * t167) * t609 + (-Ifges(7,4) * t168 - Ifges(7,2) * t167) * t611 + (-Ifges(7,5) * t168 - Ifges(7,6) * t167) * t591 + t351 * t64 + t338 * t10 + t252 * t23 + t253 * t24 + (-t84 + t520) * t146 - t161 * t225 + t546 * t613 - t334 * t616 - t168 * t617 - t335 * t618 - t167 * t619 + t349 * t624; -t696 * t95 + t121 * t94 - (-t146 - t147) * t217 + (t148 - t149) * t437 - t10 + t64 + t65 + (t121 * t14 - t13 * t696 + t11 + t418) * m(7) + (t217 * t63 - t437 * t62 + t16 + t418) * m(6) + (t217 * t74 + t437 * t73 + t418 + t60) * m(5); t397 * t23 + t394 * t24 + (t132 - t61) * t437 + t438 * qJD(6) + (-t146 - t438) * t284 + t80 + (t1 * t394 + t2 * t397 - t437 * t52 + t419 + t274 * (-t13 * t394 + t14 * t397)) * m(7) + (-t284 * t63 + t437 * t72 + t15 + t419) * m(6); -t52 * (mrSges(7,1) * t696 + mrSges(7,2) * t121) + (Ifges(7,1) * t121 - t563) * t609 + t48 * t608 + (Ifges(7,5) * t121 - Ifges(7,6) * t696) * t591 - t13 * t94 + t14 * t95 - g(1) * (mrSges(7,1) * t102 - mrSges(7,2) * t103) - g(2) * ((t178 * t397 - t179 * t394) * mrSges(7,1) + (-t178 * t394 - t179 * t397) * mrSges(7,2)) - g(3) * t451 + (t121 * t13 + t14 * t696) * mrSges(7,3) + t6 + (-Ifges(7,2) * t696 + t120 + t49) * t611 + t694;];
tau  = t7;
