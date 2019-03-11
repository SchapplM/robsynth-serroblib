% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:12:03
% EndTime: 2019-03-09 22:12:58
% DurationCPUTime: 31.66s
% Computational Cost: add. (17376->902), mult. (37261->1133), div. (0->0), fcn. (26586->12), ass. (0->420)
t415 = sin(qJ(4));
t539 = qJD(4) * t415;
t416 = sin(qJ(3));
t421 = cos(qJ(3));
t422 = cos(qJ(2));
t544 = qJD(1) * t422;
t417 = sin(qJ(2));
t545 = qJD(1) * t417;
t328 = -t416 * t545 + t421 * t544;
t572 = t328 * t415;
t754 = t539 - t572;
t731 = Ifges(5,1) + Ifges(6,1);
t702 = Ifges(6,4) + Ifges(5,5);
t536 = qJD(1) * qJD(2);
t354 = qJDD(1) * t422 - t417 * t536;
t355 = qJDD(1) * t417 + t422 * t536;
t346 = t416 * t422 + t417 * t421;
t446 = t346 * qJD(3);
t206 = -qJD(1) * t446 + t354 * t421 - t355 * t416;
t204 = qJDD(4) - t206;
t634 = t204 / 0.2e1;
t344 = t416 * t417 - t421 * t422;
t445 = t344 * qJD(3);
t205 = -qJD(1) * t445 + t354 * t416 + t355 * t421;
t412 = qJDD(2) + qJDD(3);
t420 = cos(qJ(4));
t329 = t346 * qJD(1);
t534 = qJD(2) + qJD(3);
t283 = t329 * t415 - t420 * t534;
t540 = qJD(4) * t283;
t134 = t420 * t205 + t415 * t412 - t540;
t642 = t134 / 0.2e1;
t753 = t634 * t702 + t642 * t731;
t284 = t329 * t420 + t415 * t534;
t135 = qJD(4) * t284 + t205 * t415 - t412 * t420;
t640 = t135 / 0.2e1;
t752 = mrSges(6,2) + mrSges(5,3);
t751 = -Ifges(5,4) + Ifges(6,5);
t730 = Ifges(5,6) - Ifges(6,6);
t700 = Ifges(5,3) + Ifges(6,2);
t252 = pkin(3) * t329 - pkin(9) * t328;
t527 = pkin(2) * t545;
t229 = t252 + t527;
t424 = -pkin(8) - pkin(7);
t368 = t424 * t422;
t349 = qJD(1) * t368;
t330 = t416 * t349;
t366 = t424 * t417;
t348 = qJD(1) * t366;
t262 = t348 * t421 + t330;
t152 = t229 * t415 + t262 * t420;
t316 = t329 * qJ(5);
t111 = t316 + t152;
t541 = qJD(3) * t421;
t524 = pkin(2) * t541;
t488 = t420 * t524;
t750 = -t111 + t488;
t749 = t754 * pkin(10);
t411 = t422 * pkin(2);
t399 = t411 + pkin(1);
t364 = t399 * qJD(1);
t219 = -pkin(3) * t328 - pkin(9) * t329 - t364;
t331 = t421 * t349;
t334 = qJD(2) * pkin(2) + t348;
t257 = t334 * t416 - t331;
t235 = pkin(9) * t534 + t257;
t127 = t219 * t420 - t235 * t415;
t747 = qJD(5) - t127;
t706 = -m(6) - m(7);
t746 = t640 * t751 + t753;
t745 = -mrSges(5,1) - mrSges(6,1);
t704 = mrSges(5,2) - mrSges(6,3);
t613 = pkin(2) * t416;
t397 = pkin(9) + t613;
t743 = -t397 * t539 + t749 + t750;
t253 = t415 * t262;
t600 = -pkin(10) + t397;
t340 = t600 * t420;
t489 = t415 * t524;
t425 = -pkin(4) - pkin(5);
t521 = t425 * t329;
t604 = pkin(10) * t328;
t742 = qJD(4) * t340 + t489 - t253 - (-t229 - t604) * t420 - t521;
t256 = t334 * t421 + t330;
t157 = t252 * t415 + t256 * t420;
t120 = t316 + t157;
t741 = -pkin(9) * t539 - t120 + t749;
t246 = t415 * t256;
t643 = pkin(9) - pkin(10);
t367 = t643 * t420;
t740 = qJD(4) * t367 - t246 - (-t252 - t604) * t420 - t521;
t739 = -pkin(10) * t284 + t747;
t281 = Ifges(5,4) * t283;
t324 = qJD(4) - t328;
t590 = Ifges(6,5) * t283;
t725 = t284 * t731 + t324 * t702 - t281 + t590;
t599 = mrSges(6,2) * t283;
t209 = mrSges(6,3) * t324 - t599;
t597 = mrSges(5,3) * t283;
t210 = -mrSges(5,2) * t324 - t597;
t738 = -t210 - t209;
t211 = mrSges(5,1) * t324 - mrSges(5,3) * t284;
t598 = mrSges(6,2) * t284;
t212 = -mrSges(6,1) * t324 + t598;
t737 = t211 - t212;
t261 = t348 * t416 - t331;
t542 = qJD(3) * t416;
t736 = pkin(2) * t542 - t261;
t413 = qJ(2) + qJ(3);
t406 = sin(t413);
t414 = sin(qJ(6));
t419 = cos(qJ(6));
t460 = t414 * t420 - t415 * t419;
t297 = t460 * t406;
t563 = t406 * t415;
t459 = t414 * t415 + t419 * t420;
t683 = t459 * t406;
t735 = mrSges(7,1) * t683 - mrSges(5,2) * t563 - mrSges(7,2) * t297;
t128 = t219 * t415 + t235 * t420;
t69 = t324 * t425 + t739;
t307 = t324 * qJ(5);
t89 = pkin(10) * t283 + t128;
t79 = t307 + t89;
t29 = -t414 * t79 + t419 * t69;
t30 = t414 * t69 + t419 * t79;
t692 = t534 * Ifges(4,6);
t96 = -pkin(4) * t324 + t747;
t99 = t307 + t128;
t734 = t364 * mrSges(4,1) - t127 * mrSges(5,1) + t96 * mrSges(6,1) + t29 * mrSges(7,1) + t128 * mrSges(5,2) - t30 * mrSges(7,2) - t99 * mrSges(6,3) + t692 / 0.2e1;
t693 = t534 * Ifges(4,5);
t733 = t364 * mrSges(4,2) - t693 / 0.2e1;
t407 = cos(t413);
t423 = cos(qJ(1));
t552 = t420 * t423;
t418 = sin(qJ(1));
t554 = t415 * t418;
t317 = t407 * t554 + t552;
t551 = t423 * t415;
t553 = t418 * t420;
t318 = t407 * t553 - t551;
t461 = t317 * t414 + t318 * t419;
t670 = t317 * t419 - t318 * t414;
t732 = mrSges(7,1) * t670 - t461 * mrSges(7,2);
t200 = qJDD(6) - t204;
t635 = t200 / 0.2e1;
t462 = t283 * t414 + t284 * t419;
t41 = -qJD(6) * t462 - t134 * t414 + t135 * t419;
t648 = t41 / 0.2e1;
t178 = t283 * t419 - t284 * t414;
t40 = qJD(6) * t178 + t134 * t419 + t135 * t414;
t649 = t40 / 0.2e1;
t650 = Ifges(7,1) * t649 + Ifges(7,4) * t648 + Ifges(7,5) * t635;
t651 = Ifges(7,4) * t649 + Ifges(7,2) * t648 + Ifges(7,6) * t635;
t314 = qJD(6) - t324;
t594 = Ifges(4,4) * t329;
t727 = Ifges(7,5) * t462 + Ifges(4,2) * t328 + Ifges(7,6) * t178 + Ifges(7,3) * t314 + t594 + t692;
t726 = -t283 * t730 + t284 * t702 + t324 * t700;
t580 = t329 * mrSges(4,3);
t681 = mrSges(4,1) * t534 - mrSges(5,1) * t283 - mrSges(5,2) * t284 - t580;
t538 = qJD(4) * t420;
t323 = pkin(4) * t539 - qJ(5) * t538 - qJD(5) * t415;
t724 = -pkin(5) * t754 - t323;
t723 = t752 * t406;
t721 = -t415 * t730 + t420 * t702;
t589 = Ifges(6,5) * t415;
t593 = Ifges(5,4) * t415;
t720 = t420 * t731 + t589 - t593;
t571 = t328 * t420;
t549 = pkin(4) * t572 - qJ(5) * t571;
t718 = -t549 + t736;
t717 = t538 - t571;
t342 = t355 * pkin(7);
t272 = qJDD(2) * pkin(2) - pkin(8) * t355 - t342;
t341 = t354 * pkin(7);
t282 = pkin(8) * t354 + t341;
t104 = t272 * t416 + t282 * t421 + t334 * t541 + t349 * t542;
t101 = pkin(9) * t412 + t104;
t575 = qJDD(1) * pkin(1);
t321 = -pkin(2) * t354 - t575;
t90 = -pkin(3) * t206 - pkin(9) * t205 + t321;
t22 = t101 * t420 + t219 * t538 - t235 * t539 + t415 * t90;
t23 = -t101 * t415 - t219 * t539 - t235 * t538 + t420 * t90;
t716 = t22 * t420 - t23 * t415;
t16 = qJ(5) * t204 + qJD(5) * t324 + t22;
t447 = qJDD(5) - t23;
t18 = -pkin(4) * t204 + t447;
t715 = t16 * t420 + t18 * t415;
t234 = -pkin(3) * t534 - t256;
t437 = qJ(5) * t284 - t234;
t117 = pkin(4) * t283 - t437;
t578 = t420 * mrSges(6,3);
t477 = mrSges(6,1) * t415 - t578;
t479 = mrSges(5,1) * t415 + mrSges(5,2) * t420;
t714 = t117 * t477 + t234 * t479;
t713 = Ifges(6,5) * t642 + Ifges(6,6) * t634 - Ifges(5,4) * t134 / 0.2e1 - Ifges(5,6) * t204 / 0.2e1 + (Ifges(6,3) + Ifges(5,2)) * t640;
t408 = t415 * qJ(5);
t453 = t420 * t425 - t408;
t478 = -mrSges(6,1) * t420 - mrSges(6,3) * t415;
t562 = t406 * t420;
t466 = pkin(4) * t420 + t408;
t699 = -pkin(3) - t466;
t712 = mrSges(5,1) * t562 + (-m(6) * t699 - t478 - m(7) * (-pkin(3) + t453)) * t406;
t75 = -mrSges(6,2) * t135 + mrSges(6,3) * t204;
t76 = mrSges(5,1) * t204 - mrSges(5,3) * t134;
t77 = -mrSges(6,1) * t204 + mrSges(6,2) * t134;
t78 = -mrSges(5,2) * t204 - mrSges(5,3) * t135;
t711 = m(5) * ((-t127 * t420 - t128 * t415) * qJD(4) + t716) + m(6) * ((-t415 * t99 + t420 * t96) * qJD(4) + t715) + t738 * t539 - t737 * t538 + t420 * (t75 + t78) + t415 * (-t76 + t77);
t86 = t283 * t425 + t437;
t710 = -mrSges(7,1) * t86 + mrSges(7,3) * t30;
t709 = mrSges(7,2) * t86 - mrSges(7,3) * t29;
t708 = mrSges(7,1) * t459 - mrSges(7,2) * t460;
t705 = t354 / 0.2e1;
t616 = t422 / 0.2e1;
t697 = Ifges(3,2) * t422;
t339 = t600 * t415;
t248 = t339 * t419 - t340 * t414;
t691 = qJD(6) * t248 + t414 * t742 + t419 * t743;
t249 = t339 * t414 + t340 * t419;
t690 = -qJD(6) * t249 - t414 * t743 + t419 * t742;
t365 = t643 * t415;
t285 = t365 * t419 - t367 * t414;
t689 = qJD(6) * t285 + t414 * t740 + t419 * t741;
t287 = t365 * t414 + t367 * t419;
t688 = -qJD(6) * t287 - t414 * t741 + t419 * t740;
t358 = -qJ(5) * t414 + t419 * t425;
t687 = qJD(6) * t358 - t414 * t89 + t419 * t739;
t359 = qJ(5) * t419 + t414 * t425;
t686 = -qJD(6) * t359 - t414 * t739 - t419 * t89;
t239 = t460 * t346;
t682 = -t718 + t724;
t164 = t257 + t549;
t680 = t164 + t724;
t679 = t323 + t718;
t678 = t323 - t164;
t677 = t366 * t421 + t368 * t416;
t676 = mrSges(4,1) * t407 - mrSges(4,2) * t406;
t674 = t134 * t702 - t135 * t730 + t204 * t700;
t607 = pkin(7) * t422;
t608 = pkin(7) * t417;
t672 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t545) * t607 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t544) * t608;
t671 = t341 * t422 + t342 * t417;
t669 = g(1) * t423 + g(2) * t418;
t667 = qJD(4) - qJD(6);
t666 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t363 = -mrSges(3,1) * t422 + mrSges(3,2) * t417;
t665 = -m(3) * pkin(1) - mrSges(2,1) + t363 - t676;
t664 = -m(7) * pkin(10) - mrSges(7,3);
t610 = pkin(3) * t406;
t612 = pkin(2) * t417;
t663 = -m(7) * (-pkin(10) * t407 - t612) + t407 * mrSges(7,3) + m(6) * t612 - m(5) * (-t610 - t612) + t712;
t662 = -t407 * t664 + t712;
t557 = t407 * t423;
t383 = pkin(9) * t557;
t661 = t383 * t706 + t423 * t735 - t557 * t752;
t288 = t366 * t416 - t368 * t421;
t516 = qJD(2) * t424;
t352 = t417 * t516;
t353 = t422 * t516;
t185 = qJD(3) * t288 + t352 * t416 - t353 * t421;
t658 = m(7) * pkin(5) - t745;
t558 = t407 * t420;
t657 = t406 * mrSges(7,3) + t558 * t745 - t676 - t723 + (t415 * t704 - t708) * t407;
t10 = -pkin(10) * t134 + t204 * t425 + t447;
t12 = pkin(10) * t135 + t16;
t2 = qJD(6) * t29 + t10 * t414 + t12 * t419;
t3 = -qJD(6) * t30 + t10 * t419 - t12 * t414;
t656 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t655 = (t735 + (-t752 + (-m(5) + t706) * pkin(9)) * t407) * t418;
t653 = t23 * mrSges(5,1) - t18 * mrSges(6,1) - t22 * mrSges(5,2) + t16 * mrSges(6,3);
t591 = Ifges(7,4) * t462;
t73 = Ifges(7,2) * t178 + Ifges(7,6) * t314 + t591;
t647 = -t73 / 0.2e1;
t646 = t73 / 0.2e1;
t174 = Ifges(7,4) * t178;
t74 = Ifges(7,1) * t462 + Ifges(7,5) * t314 + t174;
t645 = -t74 / 0.2e1;
t644 = t74 / 0.2e1;
t641 = -t135 / 0.2e1;
t639 = -t178 / 0.2e1;
t638 = t178 / 0.2e1;
t637 = -t462 / 0.2e1;
t636 = t462 / 0.2e1;
t631 = -t283 / 0.2e1;
t630 = t283 / 0.2e1;
t629 = -t284 / 0.2e1;
t628 = t284 / 0.2e1;
t627 = -t314 / 0.2e1;
t626 = t314 / 0.2e1;
t625 = -t324 / 0.2e1;
t624 = t324 / 0.2e1;
t622 = t328 / 0.2e1;
t620 = t329 / 0.2e1;
t611 = pkin(2) * t421;
t609 = pkin(4) * t329;
t391 = t406 * pkin(9);
t392 = t407 * pkin(3);
t596 = Ifges(3,4) * t417;
t595 = Ifges(3,4) * t422;
t592 = Ifges(5,4) * t420;
t588 = Ifges(6,5) * t420;
t587 = t128 * mrSges(5,3);
t582 = t284 * Ifges(5,4);
t581 = t328 * mrSges(4,3);
t577 = qJ(5) * t283;
t576 = qJ(5) * t420;
t268 = -qJD(2) * t344 - t445;
t574 = t268 * t420;
t568 = t346 * t415;
t567 = t346 * t420;
t561 = t406 * t423;
t280 = Ifges(6,5) * t284;
t158 = Ifges(6,6) * t324 + Ifges(6,3) * t283 + t280;
t556 = t415 * t158;
t161 = -Ifges(5,2) * t283 + Ifges(5,6) * t324 + t582;
t555 = t415 * t161;
t550 = t423 * t424;
t255 = pkin(3) * t344 - pkin(9) * t346 - t399;
t169 = t255 * t415 + t288 * t420;
t546 = t392 + t391;
t543 = qJD(2) * t417;
t537 = qJD(5) * t420;
t529 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t200;
t526 = pkin(2) * t543;
t520 = t425 * t415;
t145 = qJ(5) * t344 + t169;
t398 = -pkin(3) - t611;
t515 = t346 * t538;
t507 = t556 / 0.2e1;
t502 = -t539 / 0.2e1;
t501 = t538 / 0.2e1;
t498 = -t399 - t392;
t497 = t536 / 0.2e1;
t151 = t229 * t420 - t253;
t156 = t252 * t420 - t246;
t273 = t415 * t288;
t168 = t255 * t420 - t273;
t490 = t399 * t423 - t418 * t424;
t105 = t272 * t421 - t282 * t416 - t334 * t542 + t349 * t541;
t487 = pkin(4) * t558 + t407 * t408 + t546;
t482 = mrSges(3,1) * t417 + mrSges(3,2) * t422;
t480 = mrSges(4,1) * t406 + mrSges(4,2) * t407;
t319 = t407 * t551 - t553;
t320 = t407 * t552 + t554;
t222 = t319 * t419 - t320 * t414;
t223 = t319 * t414 + t320 * t419;
t476 = mrSges(7,1) * t222 - mrSges(7,2) * t223;
t475 = -mrSges(7,1) * t297 - mrSges(7,2) * t683;
t472 = t596 + t697;
t471 = -Ifges(5,2) * t415 + t592;
t469 = Ifges(3,5) * t422 - Ifges(3,6) * t417;
t467 = Ifges(6,3) * t415 + t588;
t465 = pkin(4) * t415 - t576;
t110 = pkin(10) * t568 + t145;
t95 = t273 + (-pkin(10) * t346 - t255) * t420 + t425 * t344;
t55 = t110 * t419 + t414 * t95;
t54 = -t110 * t414 + t419 * t95;
t336 = t398 - t466;
t143 = -mrSges(7,2) * t314 + mrSges(7,3) * t178;
t144 = mrSges(7,1) * t314 - mrSges(7,3) * t462;
t463 = t143 * t419 - t144 * t414;
t269 = qJD(2) * t346 + t446;
t167 = pkin(3) * t269 - pkin(9) * t268 + t526;
t184 = qJD(3) * t677 + t352 * t421 + t353 * t416;
t57 = t167 * t420 - t184 * t415 - t255 * t539 - t288 * t538;
t456 = pkin(3) * t557 + pkin(9) * t561 + t490;
t455 = pkin(1) * t482;
t454 = t520 + t576;
t451 = t268 * t415 + t515;
t450 = t346 * t539 - t574;
t449 = pkin(3) * t412 + t105;
t448 = t417 * (Ifges(3,1) * t422 - t596);
t56 = t167 * t415 + t184 * t420 + t255 * t538 - t288 * t539;
t443 = t529 - t656;
t442 = pkin(5) * t558 - pkin(10) * t406 + t487;
t438 = -g(1) * t319 - g(2) * t317 - g(3) * t563;
t35 = qJ(5) * t269 + qJD(5) * t344 + t56;
t433 = qJ(5) * t134 + qJD(5) * t284 + t449;
t266 = t667 * t459;
t14 = t135 * t425 + t433;
t217 = t460 * t328;
t218 = t459 * t328;
t322 = Ifges(4,4) * t328;
t231 = Ifges(4,1) * t329 + t322 + t693;
t25 = pkin(4) * t135 - t433;
t426 = (t217 * t30 + t218 * t29) * mrSges(7,3) + (-mrSges(5,2) * t449 + t746 + t753) * t415 + (Ifges(7,4) * t218 - Ifges(7,2) * t217) * t639 + (-mrSges(7,3) * t2 - 0.2e1 * t651) * t459 + (-(Ifges(7,4) * t636 + Ifges(7,2) * t638 + Ifges(7,6) * t626 + t646 + t710) * t667 + t3 * mrSges(7,3) - 0.2e1 * t650) * t460 - t86 * (mrSges(7,1) * t217 + mrSges(7,2) * t218) + t725 * (-t571 / 0.2e1 + t501) + t589 * t640 + t593 * t641 + (-t588 + t592) * t642 + (t284 * t720 + t324 * t721) * qJD(4) / 0.2e1 + (t717 * t96 - t754 * t99 + t715) * mrSges(6,2) + t25 * t478 + (-t471 / 0.2e1 + t467 / 0.2e1) * t540 - (-Ifges(4,2) * t329 + t231 + t322 + t556) * t328 / 0.2e1 + t14 * t708 + (t467 * t631 + t471 * t630 + t721 * t625 + t720 * t629 - t714 + t733) * t328 + (-Ifges(7,5) * t637 + Ifges(5,6) * t630 + Ifges(6,6) * t631 - Ifges(7,6) * t639 - Ifges(7,3) * t627 + t700 * t625 + t702 * t629 + t734) * t329 + (Ifges(7,5) * t218 - Ifges(7,6) * t217) * t627 + t218 * t645 - t217 * t647 - t539 * t587 + (Ifges(7,1) * t218 - Ifges(7,4) * t217) * t637 + Ifges(4,5) * t205 + Ifges(4,6) * t206 - t104 * mrSges(4,2) + t105 * mrSges(4,1) + t555 * t622 + t257 * t580 + t256 * t581 + Ifges(4,3) * t412 + (Ifges(7,1) * t636 + Ifges(7,4) * t638 + Ifges(7,5) * t626 + t644 + t709) * t266 + (t507 + t714) * qJD(4) + (-t127 * t717 + t128 * t572 + t716) * mrSges(5,3) - (Ifges(4,1) * t328 - t594 + t726) * t329 / 0.2e1 + t161 * t502 + t727 * t620 + (mrSges(5,1) * t449 + Ifges(5,2) * t641 - Ifges(6,3) * t640 + t634 * t730 - t713) * t420;
t409 = t420 * pkin(5);
t401 = Ifges(3,4) * t544;
t337 = t409 - t699;
t327 = Ifges(3,1) * t545 + Ifges(3,5) * qJD(2) + t401;
t326 = Ifges(3,6) * qJD(2) + qJD(1) * t472;
t306 = -t336 + t409;
t293 = -mrSges(4,2) * t534 + t581;
t251 = -mrSges(4,1) * t328 + mrSges(4,2) * t329;
t240 = t459 * t346;
t192 = -mrSges(4,2) * t412 + mrSges(4,3) * t206;
t191 = mrSges(4,1) * t412 - mrSges(4,3) * t205;
t189 = mrSges(6,1) * t283 - mrSges(6,3) * t284;
t188 = pkin(4) * t284 + t577;
t181 = t346 * t465 - t677;
t153 = t346 * t454 + t677;
t150 = -pkin(4) * t344 - t168;
t141 = t284 * t425 - t577;
t122 = -t156 - t609;
t112 = -t151 - t609;
t85 = -mrSges(7,1) * t178 + mrSges(7,2) * t462;
t71 = t266 * t346 - t268 * t460;
t70 = t239 * t667 + t268 * t459;
t63 = mrSges(5,1) * t135 + mrSges(5,2) * t134;
t62 = mrSges(6,1) * t135 - mrSges(6,3) * t134;
t61 = t465 * t268 + (qJD(4) * t466 - t537) * t346 + t185;
t49 = t454 * t268 + (qJD(4) * t453 + t537) * t346 - t185;
t44 = -pkin(4) * t269 - t57;
t28 = -mrSges(7,2) * t200 + mrSges(7,3) * t41;
t27 = mrSges(7,1) * t200 - mrSges(7,3) * t40;
t26 = pkin(10) * t451 + t35;
t24 = pkin(10) * t450 + t269 * t425 - t57;
t11 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t5 = -qJD(6) * t55 + t24 * t419 - t26 * t414;
t4 = qJD(6) * t54 + t24 * t414 + t26 * t419;
t1 = [m(4) * (t104 * t288 + t184 * t257 - t321 * t399 - t364 * t526) + (Ifges(7,4) * t240 - Ifges(7,2) * t239) * t648 + (Ifges(7,1) * t240 - Ifges(7,4) * t239) * t649 + t14 * (mrSges(7,1) * t239 + mrSges(7,2) * t240) + (-t2 * t239 - t240 * t3 - t29 * t70 + t30 * t71) * mrSges(7,3) + (Ifges(7,5) * t240 - Ifges(7,6) * t239) * t635 + (-Ifges(6,5) * t450 + Ifges(6,3) * t451) * t630 + (-Ifges(5,4) * t450 - Ifges(5,2) * t451) * t631 + (Ifges(3,4) * t355 + Ifges(3,2) * t354) * t616 + (t422 * t595 + t448) * t497 + (-m(4) * t490 - m(5) * t456 - t223 * mrSges(7,1) - t222 * mrSges(7,2) + t706 * (pkin(4) * t320 + qJ(5) * t319 + t456) - t658 * t320 + t704 * t319 + (-t664 - t752) * t561 + t665 * t423 + t666 * t418) * g(2) - (-m(4) * t105 - m(5) * t449 - t191 + t63) * t677 + (m(5) * t550 + t461 * mrSges(7,1) + t670 * mrSges(7,2) + t706 * (-pkin(4) * t318 - qJ(5) * t317 - t550) + t658 * t318 - t704 * t317 + (m(4) * t424 + t666) * t423 + (m(4) * t399 - m(7) * t498 - (-m(7) * t643 + mrSges(7,3)) * t406 + (-m(5) - m(6)) * (t498 - t391) - t665 + t723) * t418) * g(1) + (t18 * t567 - t450 * t96 - t451 * t99) * mrSges(6,2) + (Ifges(7,5) * t70 + Ifges(7,6) * t71) * t626 + (t127 * t450 - t128 * t451 - t23 * t567) * mrSges(5,3) + m(6) * (t117 * t61 + t145 * t16 + t150 * t18 + t181 * t25 + t35 * t99 + t44 * t96) + m(7) * (t14 * t153 + t2 * t55 + t29 * t5 + t3 * t54 + t30 * t4 + t49 * t86) + (-t529 / 0.2e1 - Ifges(4,6) * t412 - Ifges(4,2) * t206 - Ifges(4,4) * t205 + Ifges(6,6) * t640 + Ifges(5,6) * t641 - Ifges(7,6) * t648 - Ifges(7,5) * t649 + t321 * mrSges(4,1) - Ifges(7,3) * t635 - t104 * mrSges(4,3) + t702 * t642 + t700 * t634 + t653 + t656 + t674 / 0.2e1) * t344 + t117 * (mrSges(6,1) * t451 + mrSges(6,3) * t450) + t234 * (mrSges(5,1) * t451 - mrSges(5,2) * t450) - t363 * t575 + (-t555 / 0.2e1 - t256 * mrSges(4,3) + t231 / 0.2e1 + Ifges(4,1) * t620 + Ifges(4,4) * t622 + t507 - t733) * t268 + (-Ifges(7,5) * t636 - Ifges(7,6) * t638 + Ifges(6,6) * t630 + Ifges(5,6) * t631 - Ifges(4,4) * t620 - Ifges(4,2) * t622 - Ifges(7,3) * t626 + t700 * t624 + t702 * t628 - t257 * mrSges(4,3) + t726 / 0.2e1 - t727 / 0.2e1 - t734) * t269 + m(5) * (t127 * t57 + t128 * t56 + t168 * t23 + t169 * t22) + (Ifges(7,4) * t70 + Ifges(7,2) * t71) * t638 + (-m(4) * t256 + m(5) * t234 - t681) * t185 + t70 * t644 + t71 * t646 + t240 * t650 - t239 * t651 - t326 * t543 / 0.2e1 + (-t450 * t702 - t451 * t730) * t624 - t399 * (-mrSges(4,1) * t206 + mrSges(4,2) * t205) + t472 * t705 + t567 * t746 + (t354 * t607 + t355 * t608 + t671) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t671) + (t327 * t616 + t469 * qJD(2) / 0.2e1 - t672) * qJD(2) + (Ifges(7,1) * t70 + Ifges(7,4) * t71) * t636 + t251 * t526 - t161 * t515 / 0.2e1 - pkin(1) * (-mrSges(3,1) * t354 + mrSges(3,2) * t355) - t455 * t536 + t184 * t293 + Ifges(2,3) * qJDD(1) + t288 * t192 + t355 * t595 / 0.2e1 + t35 * t209 + t56 * t210 + t57 * t211 + t44 * t212 + t61 * t189 + t181 * t62 + t168 * t76 + t169 * t78 + t153 * t11 + t5 * t144 + t145 * t75 + t150 * t77 + t4 * t143 + t49 * t85 + t86 * (-mrSges(7,1) * t71 + mrSges(7,2) * t70) + t54 * t27 + t55 * t28 + (-t16 * mrSges(6,2) - t22 * mrSges(5,3) + t713) * t568 + (t321 * mrSges(4,2) - t105 * mrSges(4,3) + Ifges(4,1) * t205 + Ifges(4,4) * t206 + Ifges(4,5) * t412 + t158 * t501 + t25 * t477 - t449 * t479 + t467 * t640 + t471 * t641 + t634 * t721 + t642 * t720) * t346 + (-t450 * t731 + t451 * t751) * t628 + t725 * (t346 * t502 + t574 / 0.2e1) + (-mrSges(3,1) * t608 - mrSges(3,2) * t607 + 0.2e1 * Ifges(3,6) * t616) * qJDD(2) + (Ifges(3,1) * t355 + Ifges(3,4) * t705 + Ifges(3,5) * qJDD(2) - t497 * t697) * t417; (t256 * t261 - t257 * t262 + t364 * t527 + (t104 * t416 + t105 * t421 + (-t256 * t416 + t257 * t421) * qJD(3)) * pkin(2)) * m(4) - (-Ifges(3,2) * t545 + t327 + t401) * t544 / 0.2e1 + (t524 - t262) * t293 + (m(4) * t612 + t480 + t482) * t669 + (-t112 + t489) * t212 - t681 * t736 + (-t152 + t488) * t210 + t750 * t209 + (-t449 * t398 + (t234 * t416 + (-t127 * t415 + t128 * t420) * t421) * qJD(3) * pkin(2) - t127 * t151 - t128 * t152 - t234 * t261) * m(5) + t690 * t144 + (t14 * t306 + t2 * t249 + t248 * t3 + t29 * t690 + t30 * t691 + t682 * t86) * m(7) + t691 * t143 + (-t489 - t151) * t211 + t682 * t85 + t326 * t545 / 0.2e1 + t398 * t63 - t251 * t527 + Ifges(3,6) * t354 + Ifges(3,5) * t355 + (-m(5) * t383 + t423 * t663 + t661) * g(1) + (t418 * t663 + t655) * g(2) - t469 * t536 / 0.2e1 + (t25 * t336 + (t415 * t96 + t420 * t99) * t524 - t111 * t99 - t112 * t96 + t679 * t117) * m(6) + t679 * t189 + (t363 - m(7) * (t411 + t442) - m(6) * (t411 + t487) - m(5) * (t411 + t546) - m(4) * t411 + t657) * g(3) - t341 * mrSges(3,2) - t342 * mrSges(3,1) + t336 * t62 + t306 * t11 + t248 * t27 + t249 * t28 + Ifges(3,3) * qJDD(2) + t191 * t611 + t192 * t613 + t426 + t711 * t397 + (t672 + (-t448 / 0.2e1 + t455) * qJD(1)) * qJD(1); (pkin(3) * t449 - t127 * t156 - t128 * t157 - t234 * t257) * m(5) + t699 * t62 + t689 * t143 + (t14 * t337 + t2 * t287 + t285 * t3 + t29 * t688 + t30 * t689 + t680 * t86) * m(7) + t680 * t85 + t681 * t257 + t669 * t480 + t688 * t144 + ((m(5) * t610 + t662) * t418 + t655) * g(2) + (-m(5) * (-pkin(3) * t561 + t383) + t662 * t423 + t661) * g(1) + (t117 * t678 - t120 * t99 - t122 * t96 + t25 * t699) * m(6) + t678 * t189 + (-m(5) * t546 - m(6) * t487 - m(7) * t442 + t657) * g(3) + t337 * t11 - t256 * t293 + t285 * t27 + t287 * t28 - t120 * t209 - t157 * t210 - t156 * t211 - t122 * t212 - pkin(3) * t63 + t426 + t711 * pkin(9); (-m(6) * t99 - t597 + t738) * t127 + t674 + t737 * t128 + t358 * t27 + t359 * t28 + (t706 * (-pkin(4) * t317 + qJ(5) * t318) + t704 * t318 + t658 * t317 + t732) * g(2) + (-pkin(4) * t18 + qJ(5) * t16 + qJD(5) * t99 - t117 * t188 - t128 * t96) * m(6) + t653 + (-t283 * t702 - t284 * t730) * t625 + (-t283 * t731 + t158 + t280 - t582) * t629 + (t475 + t706 * qJ(5) * t562 + (-t578 - (-m(6) * pkin(4) - mrSges(6,1)) * t415 - m(7) * t520 + t479) * t406) * g(3) + (t476 + t706 * (-pkin(4) * t319 + qJ(5) * t320) + t704 * t320 + t658 * t319) * g(1) + t686 * t144 + t687 * t143 + (-t141 * t86 + t2 * t359 + t29 * t686 + t3 * t358 + t30 * t687) * m(7) - t443 - t117 * (mrSges(6,1) * t284 + mrSges(6,3) * t283) - t234 * (mrSges(5,1) * t284 - mrSges(5,2) * t283) + qJD(5) * t209 - t188 * t189 - t141 * t85 + qJ(5) * t75 - pkin(4) * t77 + t161 * t628 + (Ifges(6,3) * t284 - t590) * t631 + t99 * t598 + t96 * t599 + t284 * t587 - (Ifges(7,1) * t637 + Ifges(7,4) * t639 + Ifges(7,5) * t627 + t645 - t709) * t178 + (Ifges(7,4) * t637 + Ifges(7,2) * t639 + Ifges(7,6) * t627 + t647 - t710) * t462 + (-Ifges(5,2) * t284 - t281 + t725) * t630; t419 * t27 + t414 * t28 + (t189 - t85) * t284 + t463 * qJD(6) + (-t209 - t463) * t324 + t77 + (t2 * t414 - t284 * t86 + t3 * t419 + t438 + t314 * (-t29 * t414 + t30 * t419)) * m(7) + (t117 * t284 - t324 * t99 + t18 + t438) * m(6); -t86 * (mrSges(7,1) * t462 + mrSges(7,2) * t178) + t73 * t636 + (Ifges(7,5) * t178 - Ifges(7,6) * t462) * t627 + (Ifges(7,1) * t178 - t591) * t637 - t29 * t143 + t30 * t144 - g(1) * t476 - g(2) * t732 - g(3) * t475 + (t178 * t29 + t30 * t462) * mrSges(7,3) + t443 + (-Ifges(7,2) * t462 + t174 + t74) * t639;];
tau  = t1;
