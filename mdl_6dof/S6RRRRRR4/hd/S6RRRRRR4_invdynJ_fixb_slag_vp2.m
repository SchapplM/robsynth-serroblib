% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:45
% EndTime: 2019-03-10 03:48:22
% DurationCPUTime: 60.32s
% Computational Cost: add. (39621->1111), mult. (87802->1484), div. (0->0), fcn. (64826->18), ass. (0->491)
t445 = sin(qJ(2));
t451 = cos(qJ(2));
t495 = pkin(2) * t445 - pkin(8) * t451;
t382 = t495 * qJD(1);
t450 = cos(qJ(3));
t444 = sin(qJ(3));
t539 = qJD(1) * t445;
t515 = t444 * t539;
t294 = pkin(7) * t515 + t450 * t382;
t549 = t450 * t451;
t475 = pkin(3) * t445 - pkin(9) * t549;
t453 = -pkin(9) - pkin(8);
t516 = qJD(3) * t453;
t744 = -qJD(1) * t475 + t450 * t516 - t294;
t353 = t444 * t382;
t552 = t445 * t450;
t554 = t444 * t451;
t743 = t353 + (-pkin(7) * t552 - pkin(9) * t554) * qJD(1) - t444 * t516;
t443 = sin(qJ(4));
t449 = cos(qJ(4));
t377 = t443 * t450 + t444 * t449;
t681 = qJD(3) + qJD(4);
t285 = t681 * t377;
t470 = t377 * t451;
t318 = qJD(1) * t470;
t742 = t285 - t318;
t399 = t453 * t444;
t400 = t453 * t450;
t293 = t443 * t399 - t449 * t400;
t692 = -qJD(4) * t293 + t743 * t443 + t449 * t744;
t530 = qJD(4) * t449;
t531 = qJD(4) * t443;
t691 = t399 * t530 + t400 * t531 + t443 * t744 - t743 * t449;
t536 = qJD(2) * t450;
t374 = -t515 + t536;
t513 = t450 * t539;
t375 = qJD(2) * t444 + t513;
t274 = t374 * t443 + t375 * t449;
t442 = sin(qJ(5));
t448 = cos(qJ(5));
t497 = t449 * t374 - t375 * t443;
t202 = t274 * t448 + t442 * t497;
t441 = sin(qJ(6));
t447 = cos(qJ(6));
t710 = -t274 * t442 + t448 * t497;
t120 = t202 * t447 + t441 * t710;
t426 = pkin(7) * t539;
t397 = -qJD(2) * pkin(2) + t426;
t302 = -pkin(3) * t374 + t397;
t218 = -pkin(4) * t497 + t302;
t137 = -pkin(5) * t710 + t218;
t499 = -t202 * t441 + t447 * t710;
t525 = qJD(1) * qJD(2);
t387 = qJDD(1) * t445 + t451 * t525;
t264 = qJD(3) * t374 + qJDD(2) * t444 + t387 * t450;
t265 = -qJD(3) * t375 + qJDD(2) * t450 - t387 * t444;
t149 = qJD(4) * t497 + t264 * t449 + t265 * t443;
t150 = -qJD(4) * t274 - t264 * t443 + t265 * t449;
t69 = qJD(5) * t710 + t149 * t448 + t150 * t442;
t70 = -qJD(5) * t202 - t149 * t442 + t150 * t448;
t22 = qJD(6) * t499 + t441 * t70 + t447 * t69;
t23 = -qJD(6) * t120 - t441 * t69 + t447 * t70;
t386 = qJDD(1) * t451 - t445 * t525;
t371 = qJDD(3) - t386;
t360 = qJDD(4) + t371;
t344 = qJDD(5) + t360;
t333 = qJDD(6) + t344;
t522 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t333;
t538 = qJD(1) * t451;
t412 = qJD(3) - t538;
t402 = qJD(4) + t412;
t393 = qJD(5) + t402;
t496 = pkin(2) * t451 + pkin(8) * t445;
t392 = -pkin(1) - t496;
t361 = t392 * qJD(1);
t427 = pkin(7) * t538;
t398 = qJD(2) * pkin(8) + t427;
t279 = t450 * t361 - t398 * t444;
t233 = -pkin(9) * t375 + t279;
t220 = pkin(3) * t412 + t233;
t280 = t361 * t444 + t398 * t450;
t234 = pkin(9) * t374 + t280;
t226 = t443 * t234;
t154 = t449 * t220 - t226;
t714 = pkin(10) * t274;
t128 = t154 - t714;
t121 = pkin(4) * t402 + t128;
t228 = t449 * t234;
t155 = t220 * t443 + t228;
t708 = pkin(10) * t497;
t129 = t155 + t708;
t123 = t442 * t129;
t71 = t448 * t121 - t123;
t713 = pkin(11) * t202;
t55 = t71 - t713;
t49 = pkin(5) * t393 + t55;
t707 = pkin(11) * t710;
t125 = t448 * t129;
t72 = t121 * t442 + t125;
t56 = t72 + t707;
t574 = t441 * t56;
t24 = t447 * t49 - t574;
t568 = qJDD(1) * pkin(1);
t283 = -pkin(2) * t386 - pkin(8) * t387 - t568;
t369 = t386 * pkin(7);
t337 = qJDD(2) * pkin(8) + t369;
t180 = -qJD(3) * t280 + t450 * t283 - t337 * t444;
t131 = pkin(3) * t371 - pkin(9) * t264 + t180;
t532 = qJD(3) * t450;
t534 = qJD(3) * t444;
t179 = t444 * t283 + t450 * t337 + t361 * t532 - t398 * t534;
t141 = pkin(9) * t265 + t179;
t52 = -qJD(4) * t155 + t449 * t131 - t141 * t443;
t41 = pkin(4) * t360 - pkin(10) * t149 + t52;
t51 = t443 * t131 + t449 * t141 + t220 * t530 - t234 * t531;
t43 = pkin(10) * t150 + t51;
t13 = -qJD(5) * t72 + t448 * t41 - t43 * t442;
t6 = pkin(5) * t344 - pkin(11) * t69 + t13;
t528 = qJD(5) * t448;
t529 = qJD(5) * t442;
t12 = t121 * t528 - t129 * t529 + t442 * t41 + t448 * t43;
t7 = pkin(11) * t70 + t12;
t2 = qJD(6) * t24 + t441 * t6 + t447 * t7;
t569 = t447 * t56;
t25 = t441 * t49 + t569;
t3 = -qJD(6) * t25 - t441 * t7 + t447 * t6;
t719 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t474 = t522 - t719;
t521 = Ifges(6,5) * t69 + Ifges(6,6) * t70 + Ifges(6,3) * t344;
t604 = mrSges(7,3) * t25;
t605 = mrSges(7,3) * t24;
t389 = qJD(6) + t393;
t614 = -t389 / 0.2e1;
t639 = -t120 / 0.2e1;
t641 = -t499 / 0.2e1;
t111 = Ifges(7,4) * t499;
t65 = Ifges(7,1) * t120 + Ifges(7,5) * t389 + t111;
t651 = -t65 / 0.2e1;
t578 = Ifges(7,4) * t120;
t64 = Ifges(7,2) * t499 + Ifges(7,6) * t389 + t578;
t653 = -t64 / 0.2e1;
t717 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t741 = t474 + t521 - t717 - (mrSges(7,1) * t137 + Ifges(7,4) * t639 + Ifges(7,2) * t641 + Ifges(7,6) * t614 - t604 + t653) * t120 + (-mrSges(7,2) * t137 + Ifges(7,1) * t639 + Ifges(7,4) * t641 + Ifges(7,5) * t614 + t605 + t651) * t499;
t476 = t443 * t444 - t449 * t450;
t284 = t681 * t476;
t469 = t476 * t451;
t319 = qJD(1) * t469;
t740 = -pkin(4) * t539 + t692 + (t284 - t319) * pkin(10);
t739 = t742 * pkin(10) - t691;
t518 = Ifges(5,5) * t149 + Ifges(5,6) * t150 + Ifges(5,3) * t360;
t612 = -t393 / 0.2e1;
t628 = -t202 / 0.2e1;
t630 = -t710 / 0.2e1;
t196 = Ifges(6,4) * t710;
t109 = t202 * Ifges(6,1) + t393 * Ifges(6,5) + t196;
t643 = -t109 / 0.2e1;
t579 = Ifges(6,4) * t202;
t108 = Ifges(6,2) * t710 + t393 * Ifges(6,6) + t579;
t645 = -t108 / 0.2e1;
t674 = mrSges(6,2) * t218 - mrSges(6,3) * t71;
t676 = -mrSges(6,1) * t218 + mrSges(6,3) * t72;
t718 = -t52 * mrSges(5,1) + t51 * mrSges(5,2);
t738 = -(Ifges(6,4) * t628 + Ifges(6,2) * t630 + Ifges(6,6) * t612 + t645 - t676) * t202 + (Ifges(6,1) * t628 + Ifges(6,4) * t630 + Ifges(6,5) * t612 + t643 - t674) * t710 + t518 - t718 + t741;
t276 = t377 * t448 - t442 * t476;
t161 = -qJD(5) * t276 + t284 * t442 - t285 * t448;
t229 = -t318 * t448 + t319 * t442;
t737 = t161 - t229;
t597 = pkin(5) * t202;
t292 = t449 * t399 + t400 * t443;
t247 = -pkin(10) * t377 + t292;
t248 = -pkin(10) * t476 + t293;
t704 = t247 * t528 - t248 * t529 + t442 * t740 - t739 * t448;
t185 = t442 * t247 + t448 * t248;
t703 = -qJD(5) * t185 + t739 * t442 + t448 * t740;
t163 = -t233 * t443 - t228;
t135 = t163 - t708;
t164 = t449 * t233 - t226;
t136 = t164 - t714;
t601 = pkin(3) * t449;
t420 = pkin(4) + t601;
t559 = t442 * t443;
t700 = -t442 * t135 - t448 * t136 + t420 * t528 + (-t443 * t529 + (t448 * t449 - t559) * qJD(4)) * pkin(3);
t557 = t443 * t448;
t699 = -t448 * t135 + t136 * t442 - t420 * t529 + (-t443 * t528 + (-t442 * t449 - t557) * qJD(4)) * pkin(3);
t734 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t732 = m(4) * pkin(7);
t731 = t387 / 0.2e1;
t730 = pkin(11) * t737 + t704;
t275 = -t377 * t442 - t448 * t476;
t160 = qJD(5) * t275 - t284 * t448 - t285 * t442;
t230 = -t318 * t442 - t319 * t448;
t729 = -pkin(5) * t539 + t703 + (-t160 + t230) * pkin(11);
t728 = t700 + t713;
t727 = t707 + t699;
t514 = t444 * t538;
t688 = -t427 + (-t514 + t534) * pkin(3);
t493 = mrSges(3,1) * t445 + mrSges(3,2) * t451;
t584 = Ifges(3,4) * t445;
t723 = t445 * (Ifges(3,1) * t451 - t584) / 0.2e1 - pkin(1) * t493;
t440 = qJ(3) + qJ(4);
t431 = cos(t440);
t416 = pkin(4) * t431;
t436 = t450 * pkin(3);
t391 = t416 + t436;
t432 = qJ(5) + t440;
t418 = cos(t432);
t408 = pkin(5) * t418;
t330 = t408 + t391;
t324 = pkin(2) + t330;
t381 = pkin(2) + t391;
t425 = qJ(6) + t432;
t409 = sin(t425);
t410 = cos(t425);
t417 = sin(t432);
t421 = t436 + pkin(2);
t430 = sin(t440);
t492 = -mrSges(4,1) * t450 + mrSges(4,2) * t444;
t721 = m(4) * pkin(2) + m(5) * t421 + m(6) * t381 + m(7) * t324 + mrSges(5,1) * t431 + mrSges(6,1) * t418 + mrSges(7,1) * t410 - mrSges(5,2) * t430 - mrSges(6,2) * t417 - mrSges(7,2) * t409 - t492;
t439 = -pkin(10) + t453;
t433 = -pkin(11) + t439;
t720 = -m(4) * pkin(8) + m(5) * t453 + m(6) * t439 + m(7) * t433 - t734;
t716 = -t180 * mrSges(4,1) + t179 * mrSges(4,2);
t485 = t451 * Ifges(3,2) + t584;
t715 = t155 * mrSges(5,2) + t25 * mrSges(7,2) + t72 * mrSges(6,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t485 / 0.2e1 - t154 * mrSges(5,1) - t24 * mrSges(7,1) - t71 * mrSges(6,1);
t600 = pkin(4) * t274;
t690 = t742 * pkin(4) + t688;
t535 = qJD(2) * t451;
t463 = t444 * t535 + t445 * t532;
t660 = m(5) * pkin(3);
t657 = t22 / 0.2e1;
t656 = t23 / 0.2e1;
t649 = t69 / 0.2e1;
t648 = t70 / 0.2e1;
t637 = t149 / 0.2e1;
t636 = t150 / 0.2e1;
t626 = t264 / 0.2e1;
t625 = t265 / 0.2e1;
t620 = t333 / 0.2e1;
t619 = t344 / 0.2e1;
t618 = t360 / 0.2e1;
t617 = t371 / 0.2e1;
t184 = t448 * t247 - t248 * t442;
t145 = -pkin(11) * t276 + t184;
t146 = pkin(11) * t275 + t185;
t90 = t145 * t441 + t146 * t447;
t706 = -qJD(6) * t90 - t441 * t730 + t447 * t729;
t89 = t145 * t447 - t146 * t441;
t705 = qJD(6) * t89 + t441 * t729 + t447 * t730;
t598 = pkin(4) * t448;
t419 = pkin(5) + t598;
t526 = qJD(6) * t447;
t527 = qJD(6) * t441;
t560 = t441 * t442;
t73 = -t128 * t442 - t125;
t57 = t73 - t707;
t74 = t448 * t128 - t123;
t58 = t74 - t713;
t702 = -t441 * t57 - t447 * t58 + t419 * t526 + (-t442 * t527 + (t447 * t448 - t560) * qJD(5)) * pkin(4);
t558 = t442 * t447;
t701 = t441 * t58 - t447 * t57 - t419 * t527 + (-t442 * t526 + (-t441 * t448 - t558) * qJD(5)) * pkin(4);
t340 = -pkin(3) * t559 + t448 * t420;
t336 = pkin(5) + t340;
t342 = pkin(3) * t557 + t420 * t442;
t244 = t336 * t441 + t342 * t447;
t698 = -qJD(6) * t244 - t441 * t728 + t447 * t727;
t243 = t336 * t447 - t342 * t441;
t697 = qJD(6) * t243 + t441 * t727 + t447 * t728;
t696 = t660 + mrSges(4,1);
t335 = t476 * t445;
t693 = -pkin(5) * t737 + t690;
t373 = t450 * t392;
t278 = -pkin(9) * t552 + t373 + (-pkin(7) * t444 - pkin(3)) * t451;
t414 = pkin(7) * t549;
t313 = t444 * t392 + t414;
t556 = t444 * t445;
t287 = -pkin(9) * t556 + t313;
t209 = t449 * t278 - t287 * t443;
t175 = -pkin(4) * t451 + pkin(10) * t335 + t209;
t210 = t443 * t278 + t449 * t287;
t334 = t377 * t445;
t186 = -pkin(10) * t334 + t210;
t101 = t442 * t175 + t448 * t186;
t689 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t374 - mrSges(4,2) * t375 - mrSges(3,3) * t539;
t488 = -mrSges(7,1) * t409 - mrSges(7,2) * t410;
t687 = mrSges(6,1) * t417 + mrSges(6,2) * t418 - t488;
t446 = sin(qJ(1));
t452 = cos(qJ(1));
t548 = t451 * t452;
t310 = -t417 * t548 + t418 * t446;
t311 = t417 * t446 + t418 * t548;
t298 = -t409 * t548 + t410 * t446;
t299 = t409 * t446 + t410 * t548;
t546 = t298 * mrSges(7,1) - t299 * mrSges(7,2);
t686 = -t310 * mrSges(6,1) + t311 * mrSges(6,2) - t546;
t551 = t446 * t451;
t308 = t417 * t551 + t418 * t452;
t309 = t417 * t452 - t418 * t551;
t296 = t409 * t551 + t410 * t452;
t297 = t409 * t452 - t410 * t551;
t547 = -t296 * mrSges(7,1) + t297 * mrSges(7,2);
t685 = t308 * mrSges(6,1) - t309 * mrSges(6,2) - t547;
t365 = Ifges(4,4) * t374;
t251 = t375 * Ifges(4,1) + t412 * Ifges(4,5) + t365;
t424 = Ifges(3,4) * t538;
t684 = Ifges(3,1) * t539 + Ifges(3,5) * qJD(2) + t450 * t251 + t424;
t370 = t387 * pkin(7);
t683 = t369 * t451 + t370 * t445;
t682 = t179 * t450 - t180 * t444;
t680 = mrSges(5,1) * t430 + mrSges(5,2) * t431 + t687;
t327 = -t430 * t548 + t431 * t446;
t328 = t430 * t446 + t431 * t548;
t679 = -t327 * mrSges(5,1) + t328 * mrSges(5,2) + t686;
t325 = t430 * t551 + t431 * t452;
t326 = t430 * t452 - t431 * t551;
t678 = t325 * mrSges(5,1) - t326 * mrSges(5,2) + t685;
t677 = -mrSges(5,1) * t302 + mrSges(5,3) * t155;
t675 = mrSges(5,2) * t302 - mrSges(5,3) * t154;
t673 = -m(7) - m(6) - m(5) - m(4) - m(3);
t672 = t375 * Ifges(4,5) + t274 * Ifges(5,5) + t202 * Ifges(6,5) + t120 * Ifges(7,5) + t374 * Ifges(4,6) + Ifges(5,6) * t497 + Ifges(6,6) * t710 + Ifges(7,6) * t499 + t412 * Ifges(4,3) + t402 * Ifges(5,3) + t393 * Ifges(6,3) + t389 * Ifges(7,3);
t596 = pkin(5) * t417;
t599 = pkin(4) * t430;
t366 = -t596 - t599;
t602 = pkin(3) * t444;
t329 = -t366 + t602;
t390 = t599 + t602;
t671 = -m(6) * t390 - m(7) * t329;
t494 = mrSges(3,1) * t451 - mrSges(3,2) * t445;
t669 = t445 * t734 + mrSges(2,1) + t494;
t668 = mrSges(2,2) - mrSges(3,3) + t671;
t662 = Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t620;
t661 = Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t620;
t659 = m(6) * pkin(4);
t658 = m(7) * pkin(5);
t655 = Ifges(6,4) * t649 + Ifges(6,2) * t648 + Ifges(6,6) * t619;
t654 = Ifges(6,1) * t649 + Ifges(6,4) * t648 + Ifges(6,5) * t619;
t652 = t64 / 0.2e1;
t650 = t65 / 0.2e1;
t647 = Ifges(5,4) * t637 + Ifges(5,2) * t636 + Ifges(5,6) * t618;
t646 = Ifges(5,1) * t637 + Ifges(5,4) * t636 + Ifges(5,5) * t618;
t644 = t108 / 0.2e1;
t642 = t109 / 0.2e1;
t640 = t499 / 0.2e1;
t638 = t120 / 0.2e1;
t635 = Ifges(4,1) * t626 + Ifges(4,4) * t625 + Ifges(4,5) * t617;
t580 = Ifges(5,4) * t274;
t189 = Ifges(5,2) * t497 + t402 * Ifges(5,6) + t580;
t634 = -t189 / 0.2e1;
t633 = t189 / 0.2e1;
t268 = Ifges(5,4) * t497;
t190 = t274 * Ifges(5,1) + t402 * Ifges(5,5) + t268;
t632 = -t190 / 0.2e1;
t631 = t190 / 0.2e1;
t629 = t710 / 0.2e1;
t627 = t202 / 0.2e1;
t624 = -t497 / 0.2e1;
t623 = t497 / 0.2e1;
t622 = -t274 / 0.2e1;
t621 = t274 / 0.2e1;
t615 = t375 / 0.2e1;
t613 = t389 / 0.2e1;
t611 = t393 / 0.2e1;
t610 = -t402 / 0.2e1;
t609 = t402 / 0.2e1;
t603 = pkin(3) * t375;
t593 = g(3) * t445;
t434 = t445 * pkin(7);
t588 = mrSges(5,3) * t497;
t587 = mrSges(5,3) * t274;
t586 = mrSges(6,3) * t710;
t585 = mrSges(6,3) * t202;
t583 = Ifges(3,4) * t451;
t582 = Ifges(4,4) * t444;
t581 = Ifges(4,4) * t450;
t577 = t279 * mrSges(4,3);
t576 = t280 * mrSges(4,3);
t575 = t375 * Ifges(4,4);
t555 = t444 * t446;
t553 = t444 * t452;
t385 = t495 * qJD(2);
t537 = qJD(2) * t445;
t520 = pkin(7) * t537;
t541 = t450 * t385 + t444 * t520;
t388 = pkin(3) * t556 + t434;
t533 = qJD(3) * t445;
t429 = pkin(7) * t535;
t517 = Ifges(4,5) * t264 + Ifges(4,6) * t265 + Ifges(4,3) * t371;
t307 = pkin(3) * t463 + t429;
t250 = t374 * Ifges(4,2) + t412 * Ifges(4,6) + t575;
t510 = -t444 * t250 / 0.2e1;
t100 = t448 * t175 - t186 * t442;
t320 = pkin(4) * t476 - t421;
t281 = pkin(4) * t334 + t388;
t225 = t603 + t600;
t338 = -qJDD(2) * pkin(2) + t370;
t491 = mrSges(4,1) * t444 + mrSges(4,2) * t450;
t487 = Ifges(4,1) * t450 - t582;
t486 = Ifges(4,1) * t444 + t581;
t484 = -Ifges(4,2) * t444 + t581;
t483 = Ifges(4,2) * t450 + t582;
t482 = Ifges(3,5) * t451 - Ifges(3,6) * t445;
t481 = Ifges(4,5) * t450 - Ifges(4,6) * t444;
t480 = Ifges(4,5) * t444 + Ifges(4,6) * t450;
t239 = -t334 * t442 - t335 * t448;
t93 = -pkin(5) * t451 - pkin(11) * t239 + t100;
t238 = -t334 * t448 + t335 * t442;
t95 = pkin(11) * t238 + t101;
t45 = -t441 * t95 + t447 * t93;
t46 = t441 * t93 + t447 * t95;
t171 = t238 * t447 - t239 * t441;
t172 = t238 * t441 + t239 * t447;
t203 = t275 * t447 - t276 * t441;
t204 = t275 * t441 + t276 * t447;
t479 = t324 * t451 - t433 * t445;
t478 = t381 * t451 - t439 * t445;
t477 = t451 * t421 - t445 * t453;
t214 = -qJD(2) * t470 + t335 * t681;
t187 = -pkin(4) * t214 + t307;
t350 = -t444 * t548 + t446 * t450;
t348 = t444 * t551 + t450 * t452;
t472 = t397 * t491;
t213 = -qJD(2) * t469 - t285 * t445;
t208 = t475 * qJD(2) + (-t414 + (pkin(9) * t445 - t392) * t444) * qJD(3) + t541;
t231 = t444 * t385 + t392 * t532 + (-t445 * t536 - t451 * t534) * pkin(7);
t212 = -pkin(9) * t463 + t231;
t99 = -qJD(4) * t210 + t449 * t208 - t212 * t443;
t81 = pkin(4) * t537 - pkin(10) * t213 + t99;
t98 = t443 * t208 + t449 * t212 + t278 * t530 - t287 * t531;
t85 = pkin(10) * t214 + t98;
t28 = t175 * t528 - t186 * t529 + t442 * t81 + t448 * t85;
t222 = -pkin(3) * t265 + t338;
t464 = -t444 * t533 + t450 * t535;
t106 = -pkin(4) * t150 + t222;
t462 = Ifges(4,5) * t445 + t451 * t487;
t461 = Ifges(4,6) * t445 + t451 * t484;
t460 = Ifges(4,3) * t445 + t451 * t481;
t29 = -qJD(5) * t101 - t442 * t85 + t448 * t81;
t395 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t538;
t367 = t408 + t416;
t362 = t491 * t445;
t351 = t450 * t548 + t555;
t349 = -t446 * t549 + t553;
t341 = pkin(4) * t558 + t419 * t441;
t339 = -pkin(4) * t560 + t419 * t447;
t312 = -pkin(7) * t554 + t373;
t301 = mrSges(4,1) * t412 - mrSges(4,3) * t375;
t300 = -mrSges(4,2) * t412 + mrSges(4,3) * t374;
t295 = -pkin(7) * t513 + t353;
t236 = mrSges(5,1) * t402 - t587;
t235 = -mrSges(5,2) * t402 + t588;
t232 = -qJD(3) * t313 + t541;
t224 = -mrSges(4,2) * t371 + mrSges(4,3) * t265;
t223 = mrSges(4,1) * t371 - mrSges(4,3) * t264;
t221 = -pkin(5) * t275 + t320;
t207 = -mrSges(5,1) * t497 + mrSges(5,2) * t274;
t195 = -mrSges(4,1) * t265 + mrSges(4,2) * t264;
t194 = -pkin(5) * t238 + t281;
t183 = mrSges(6,1) * t393 - t585;
t182 = -mrSges(6,2) * t393 + t586;
t176 = t264 * Ifges(4,4) + t265 * Ifges(4,2) + t371 * Ifges(4,6);
t158 = t229 * t441 + t230 * t447;
t157 = t229 * t447 - t230 * t441;
t151 = t600 + t597;
t142 = t225 + t597;
t139 = -mrSges(5,2) * t360 + mrSges(5,3) * t150;
t138 = mrSges(5,1) * t360 - mrSges(5,3) * t149;
t122 = -mrSges(6,1) * t710 + mrSges(6,2) * t202;
t105 = mrSges(7,1) * t389 - mrSges(7,3) * t120;
t104 = -mrSges(7,2) * t389 + mrSges(7,3) * t499;
t103 = -qJD(5) * t239 - t213 * t442 + t214 * t448;
t102 = qJD(5) * t238 + t213 * t448 + t214 * t442;
t94 = -mrSges(5,1) * t150 + mrSges(5,2) * t149;
t88 = -pkin(5) * t103 + t187;
t76 = -qJD(6) * t204 - t160 * t441 + t161 * t447;
t75 = qJD(6) * t203 + t160 * t447 + t161 * t441;
t66 = -mrSges(7,1) * t499 + mrSges(7,2) * t120;
t62 = -mrSges(6,2) * t344 + mrSges(6,3) * t70;
t61 = mrSges(6,1) * t344 - mrSges(6,3) * t69;
t48 = -qJD(6) * t172 - t102 * t441 + t103 * t447;
t47 = qJD(6) * t171 + t102 * t447 + t103 * t441;
t44 = -pkin(5) * t70 + t106;
t36 = -mrSges(6,1) * t70 + mrSges(6,2) * t69;
t27 = t447 * t55 - t574;
t26 = -t441 * t55 - t569;
t19 = -mrSges(7,2) * t333 + mrSges(7,3) * t23;
t18 = mrSges(7,1) * t333 - mrSges(7,3) * t22;
t17 = pkin(11) * t103 + t28;
t16 = pkin(5) * t537 - pkin(11) * t102 + t29;
t10 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t5 = -qJD(6) * t46 + t16 * t447 - t17 * t441;
t4 = qJD(6) * t45 + t16 * t441 + t17 * t447;
t1 = [(t584 + t485) * t386 / 0.2e1 + (-t553 * t660 - t349 * mrSges(4,1) - t326 * mrSges(5,1) - t309 * mrSges(6,1) - t297 * mrSges(7,1) - t348 * mrSges(4,2) - t325 * mrSges(5,2) - t308 * mrSges(6,2) - t296 * mrSges(7,2) + (-m(7) * (-pkin(1) - t479) - m(6) * (-pkin(1) - t478) - m(5) * (-pkin(1) - t477) - m(4) * t392 + m(3) * pkin(1) + t669) * t446 + (pkin(7) * t673 + t668) * t452) * g(1) + (Ifges(6,1) * t102 + Ifges(6,4) * t103) * t627 + (Ifges(6,1) * t239 + Ifges(6,4) * t238) * t649 + (t171 * t2 - t172 * t3 - t24 * t47 + t25 * t48) * mrSges(7,3) + t494 * t568 + (-t555 * t660 - t351 * mrSges(4,1) - t328 * mrSges(5,1) - t311 * mrSges(6,1) - t299 * mrSges(7,1) - t350 * mrSges(4,2) - t327 * mrSges(5,2) - t310 * mrSges(6,2) - t298 * mrSges(7,2) + t673 * (t452 * pkin(1) + t446 * pkin(7)) + t668 * t446 + (-m(4) * t496 - m(5) * t477 - m(6) * t478 - m(7) * t479 - t669) * t452) * g(2) + (-t280 * mrSges(4,2) + Ifges(5,3) * t609 + Ifges(6,3) * t611 + Ifges(7,3) * t613 + t279 * mrSges(4,1) + Ifges(5,5) * t621 + Ifges(5,6) * t623 + Ifges(6,5) * t627 + Ifges(6,6) * t629 + Ifges(7,5) * t638 + Ifges(7,6) * t640 - t715 + t672 / 0.2e1) * t537 + t683 * mrSges(3,3) + (-t154 * t213 + t155 * t214 - t334 * t51 + t335 * t52) * mrSges(5,3) + t723 * t525 - t689 * t429 + (Ifges(6,5) * t102 + Ifges(6,6) * t103) * t611 + (Ifges(6,5) * t239 + Ifges(6,6) * t238) * t619 - (t250 * t450 + t251 * t444) * t533 / 0.2e1 + m(4) * (t179 * t313 + t180 * t312 + t231 * t280 + t232 * t279) + (-t102 * t71 + t103 * t72 + t12 * t238 - t13 * t239) * mrSges(6,3) + (Ifges(6,4) * t102 + Ifges(6,2) * t103) * t629 + (Ifges(6,4) * t239 + Ifges(6,2) * t238) * t648 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t640 + (Ifges(7,4) * t172 + Ifges(7,2) * t171) * t656 + (-Ifges(5,4) * t335 - Ifges(5,2) * t334) * t636 + (Ifges(5,4) * t213 + Ifges(5,2) * t214) * t623 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t613 + (Ifges(7,5) * t172 + Ifges(7,6) * t171) * t620 + (-Ifges(5,5) * t335 - Ifges(5,6) * t334) * t618 + (Ifges(5,5) * t213 + Ifges(5,6) * t214) * t609 + (-t179 * t556 - t180 * t552 - t279 * t464 - t280 * t463) * mrSges(4,3) + t583 * t731 + t222 * (mrSges(5,1) * t334 - mrSges(5,2) * t335) + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t638 + (Ifges(7,1) * t172 + Ifges(7,4) * t171) * t657 + t88 * t66 - t176 * t556 / 0.2e1 + t374 * (qJD(2) * t461 - t483 * t533) / 0.2e1 + t412 * (qJD(2) * t460 - t480 * t533) / 0.2e1 + m(5) * (t154 * t99 + t155 * t98 + t209 * t52 + t210 * t51 + t222 * t388 + t302 * t307) + m(6) * (t100 * t13 + t101 * t12 + t106 * t281 + t187 * t218 + t28 * t72 + t29 * t71) + m(7) * (t137 * t88 + t194 * t44 + t2 * t46 + t24 * t5 + t25 * t4 + t3 * t45) + t397 * (mrSges(4,1) * t463 + mrSges(4,2) * t464) + (-Ifges(5,1) * t335 - Ifges(5,4) * t334) * t637 + (Ifges(5,1) * t213 + Ifges(5,4) * t214) * t621 - t395 * t520 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t683) + Ifges(2,3) * qJDD(1) + t46 * t19 + t45 * t18 + qJD(2) ^ 2 * t482 / 0.2e1 - pkin(1) * (-mrSges(3,1) * t386 + mrSges(3,2) * t387) + t388 * t94 + t338 * t362 + t312 * t223 + t313 * t224 + t307 * t207 + t231 * t300 + t232 * t301 + t302 * (-mrSges(5,1) * t214 + mrSges(5,2) * t213) + t281 * t36 + t106 * (-mrSges(6,1) * t238 + mrSges(6,2) * t239) + t99 * t236 + (-qJDD(2) * mrSges(3,1) + mrSges(3,3) * t387 + t195) * t434 + t98 * t235 + t100 * t61 + t101 * t62 + t4 * t104 + t5 * t105 + t137 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + (t716 + t719 + t717 + t718 + (pkin(7) * mrSges(3,3) + Ifges(3,2) / 0.2e1) * t386 + (-Ifges(3,2) * t445 + t583) * t525 / 0.2e1 - t517 / 0.2e1 - t518 / 0.2e1 - t522 / 0.2e1 - t521 / 0.2e1 - Ifges(5,3) * t618 - Ifges(5,6) * t636 - Ifges(5,5) * t637 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * qJDD(2) - Ifges(4,3) * t617 - Ifges(6,3) * t619 - Ifges(7,3) * t620 - Ifges(4,6) * t625 - Ifges(4,5) * t626 - Ifges(6,6) * t648 - Ifges(6,5) * t649 - Ifges(7,6) * t656 - Ifges(7,5) * t657 + Ifges(3,4) * t731) * t451 + (t397 * t732 + t510 + t684 / 0.2e1) * t535 + (Ifges(3,1) * t387 + Ifges(3,5) * qJDD(2) + t338 * t732 + t481 * t617 + t484 * t625 + t487 * t626) * t445 + (qJD(2) * t462 - t486 * t533) * t615 + t213 * t631 + t214 * t633 + t552 * t635 + t102 * t642 + t103 * t644 - t335 * t646 - t334 * t647 + t47 * t650 + t48 * t652 + t239 * t654 + t238 * t655 + t172 * t661 + t171 * t662 + t44 * (-mrSges(7,1) * t171 + mrSges(7,2) * t172) + t28 * t182 + t29 * t183 + t187 * t122 + t194 * t10 + t209 * t138 + t210 * t139 + t218 * (-mrSges(6,1) * t103 + mrSges(6,2) * t102); (Ifges(7,4) * t158 + Ifges(7,2) * t157) * t641 + (Ifges(6,4) * t230 + Ifges(6,2) * t229) * t630 + t222 * (mrSges(5,1) * t476 + mrSges(5,2) * t377) + (g(1) * t452 + g(2) * t446) * (t445 * t721 + t451 * t720 + t493) + (t445 * t720 - t451 * t721 - t494) * g(3) + (-t157 * t25 + t158 * t24 + t2 * t203 - t204 * t3) * mrSges(7,3) + (Ifges(7,1) * t158 + Ifges(7,4) * t157) * t639 + (Ifges(6,1) * t230 + Ifges(6,4) * t229) * t628 + (Ifges(5,5) * t377 - Ifges(5,6) * t476) * t618 + (Ifges(5,4) * t377 - Ifges(5,2) * t476) * t636 + (Ifges(5,1) * t377 - Ifges(5,4) * t476) * t637 - t476 * t647 + (t12 * t275 - t13 * t276 - t229 * t72 + t230 * t71) * mrSges(6,3) + (Ifges(5,5) * t622 + Ifges(6,5) * t628 + Ifges(7,5) * t639 + Ifges(5,6) * t624 + Ifges(6,6) * t630 + Ifges(7,6) * t641 + Ifges(5,3) * t610 + Ifges(6,3) * t612 + Ifges(7,3) * t614 + t715) * t539 - t723 * qJD(1) ^ 2 + t693 * t66 + t688 * t207 + t689 * t427 + t690 * t122 + t691 * t235 + t692 * t236 + (t154 * t692 + t155 * t691 - t222 * t421 + t292 * t52 + t293 * t51 + t302 * t688) * m(5) + ((-t158 + t75) * mrSges(7,2) + (t157 - t76) * mrSges(7,1)) * t137 + (-pkin(2) * t338 - t279 * t294 - t280 * t295 - t397 * t427) * m(4) + (Ifges(6,5) * t230 + Ifges(6,6) * t229) * t612 + (Ifges(7,5) * t158 + Ifges(7,6) * t157) * t614 + (t374 * t484 + t375 * t487 + t412 * t481) * qJD(3) / 0.2e1 - (t374 * t461 + t375 * t462 + t412 * t460) * qJD(1) / 0.2e1 + (-t280 * (-mrSges(4,2) * t445 - mrSges(4,3) * t554) - t279 * (mrSges(4,1) * t445 - mrSges(4,3) * t549)) * qJD(1) + (t510 + t472) * qJD(3) + t703 * t183 + t704 * t182 + (t106 * t320 + t12 * t185 + t13 * t184 + t218 * t690 + t703 * t71 + t704 * t72) * m(6) + t705 * t104 + (t137 * t693 + t2 * t90 + t221 * t44 + t24 * t706 + t25 * t705 + t3 * t89) * m(7) + t706 * t105 + (-Ifges(5,4) * t319 - Ifges(5,2) * t318) * t624 + (-t154 * t319 + t155 * t318 - t377 * t52 - t476 * t51) * mrSges(5,3) + (-Ifges(5,5) * t319 - Ifges(5,6) * t318) * t610 + (-Ifges(5,1) * t319 - Ifges(5,4) * t318) * t622 - t302 * (mrSges(5,1) * t318 - mrSges(5,2) * t319) + t89 * t18 + t90 * t19 + (t251 / 0.2e1 - t577) * t532 - t672 * t539 / 0.2e1 + (Ifges(7,1) * t638 + Ifges(7,4) * t640 + Ifges(7,5) * t613 - t605 + t650) * t75 - t482 * t525 / 0.2e1 + (Ifges(6,1) * t627 + Ifges(6,4) * t629 + Ifges(6,5) * t611 + t642 + t674) * t160 - (Ifges(5,1) * t621 + Ifges(5,4) * t623 + Ifges(5,5) * t609 + t631 + t675) * t284 + t395 * t426 - t534 * t576 + (Ifges(7,4) * t638 + Ifges(7,2) * t640 + Ifges(7,6) * t613 + t604 + t652) * t76 + t250 * t514 / 0.2e1 + (-t444 * t223 + t450 * t224 - t301 * t532 - t300 * t534 + m(4) * ((-t279 * t450 - t280 * t444) * qJD(3) + t682)) * pkin(8) + t682 * mrSges(4,3) - (-Ifges(3,2) * t539 + t424 + t684) * t538 / 0.2e1 + t338 * t492 + (Ifges(6,4) * t627 + Ifges(6,2) * t629 + Ifges(6,6) * t611 + t644 + t676) * t161 - (Ifges(5,4) * t621 + Ifges(5,2) * t623 + Ifges(5,6) * t609 + t633 + t677) * t285 + t450 * t176 / 0.2e1 + Ifges(3,3) * qJDD(2) - t421 * t94 + Ifges(3,6) * t386 + Ifges(3,5) * t387 - t369 * mrSges(3,2) - t370 * mrSges(3,1) + t320 * t36 - t295 * t300 - t294 * t301 + t292 * t138 + t293 * t139 + t106 * (-mrSges(6,1) * t275 + mrSges(6,2) * t276) - t218 * (-mrSges(6,1) * t229 + mrSges(6,2) * t230) - t472 * t538 + t480 * t617 + (Ifges(6,5) * t276 + Ifges(6,6) * t275) * t619 + (Ifges(7,5) * t204 + Ifges(7,6) * t203) * t620 + t483 * t625 + t486 * t626 - t319 * t632 - t318 * t634 + t444 * t635 + t230 * t643 + t229 * t645 + t377 * t646 + (Ifges(6,4) * t276 + Ifges(6,2) * t275) * t648 + (Ifges(6,1) * t276 + Ifges(6,4) * t275) * t649 + t158 * t651 + t157 * t653 + t276 * t654 + t275 * t655 + (Ifges(7,4) * t204 + Ifges(7,2) * t203) * t656 + (Ifges(7,1) * t204 + Ifges(7,4) * t203) * t657 + t204 * t661 + t203 * t662 + t184 * t61 + t185 * t62 - pkin(2) * t195 + t44 * (-mrSges(7,1) * t203 + mrSges(7,2) * t204) + t221 * t10; -(Ifges(5,4) * t622 + Ifges(5,2) * t624 + Ifges(5,6) * t610 + t634 - t677) * t274 - t716 - (-Ifges(4,2) * t375 + t251 + t365) * t374 / 0.2e1 + (t139 * t443 + t235 * t530 - t236 * t531) * pkin(3) + t375 * t576 + t374 * t577 + (-m(6) * (-t390 * t551 - t391 * t452) - m(7) * (-t329 * t551 - t330 * t452) - mrSges(4,2) * t349 + t696 * t348 + t678) * g(2) + (-m(6) * (-t390 * t548 + t391 * t446) - m(7) * (-t329 * t548 + t330 * t446) + mrSges(4,2) * t351 - t696 * t350 + t679) * g(1) + t697 * t104 + (-t137 * t142 + t2 * t244 + t24 * t698 + t243 * t3 + t25 * t697) * m(7) + t698 * t105 + t699 * t183 + (t12 * t342 + t13 * t340 - t218 * t225 + t699 * t71 + t700 * t72) * m(6) + t700 * t182 + t738 - t375 * (Ifges(4,1) * t374 - t575) / 0.2e1 + (Ifges(5,1) * t622 + Ifges(5,4) * t624 + Ifges(5,5) * t610 + t632 - t675) * t497 + t517 - t207 * t603 - m(5) * (t154 * t163 + t155 * t164 + t302 * t603) + (m(5) * t602 - t671 + t680) * t593 + t138 * t601 - t412 * (Ifges(4,5) * t374 - Ifges(4,6) * t375) / 0.2e1 - t397 * (mrSges(4,1) * t375 + mrSges(4,2) * t374) + g(3) * t362 + t340 * t61 + t342 * t62 - t279 * t300 + t280 * t301 + t243 * t18 + t244 * t19 - t163 * t236 - t164 * t235 - t225 * t122 - t142 * t66 + t250 * t615 + (t443 * t51 + t449 * t52 + (-t154 * t443 + t155 * t449) * qJD(4)) * t660; (t182 * t528 - t183 * t529 + t442 * t62) * pkin(4) - t302 * (mrSges(5,1) * t274 + mrSges(5,2) * t497) + (Ifges(5,5) * t497 - Ifges(5,6) * t274) * t610 + (Ifges(5,1) * t497 - t580) * t622 + (t588 - t235) * t154 + (t587 + t236) * t155 + t701 * t105 + (-t137 * t151 + t2 * t341 + t24 * t701 + t25 * t702 + t3 * t339) * m(7) + t702 * t104 + t738 - t122 * t600 - m(6) * (t218 * t600 + t71 * t73 + t72 * t74) + (-m(7) * (t366 * t551 - t367 * t452) + t325 * t659 + t678) * g(2) + (-t327 * t659 - m(7) * (t366 * t548 + t367 * t446) + t679) * g(1) + (m(6) * t599 - m(7) * t366 + t680) * t593 + t61 * t598 + (-Ifges(5,2) * t274 + t190 + t268) * t624 + t339 * t18 + t341 * t19 - t151 * t66 + t189 * t621 + (t12 * t442 + t13 * t448 + (-t442 * t71 + t448 * t72) * qJD(5)) * t659 - t74 * t182 - t73 * t183; -t66 * t597 - m(7) * (t137 * t597 + t24 * t26 + t25 * t27) + (Ifges(6,5) * t710 - Ifges(6,6) * t202) * t612 - t27 * t104 - t26 * t105 + t108 * t627 + (Ifges(6,1) * t710 - t579) * t628 + (t2 * t441 + t3 * t447 + (-t24 * t441 + t25 * t447) * qJD(6)) * t658 - t218 * (mrSges(6,1) * t202 + mrSges(6,2) * t710) + (t585 + t183) * t72 + (t586 - t182) * t71 + (-Ifges(6,2) * t202 + t109 + t196) * t630 + (m(7) * t596 + t687) * t593 + (t308 * t658 + t685) * g(2) + (-t310 * t658 + t686) * g(1) + (t104 * t526 - t105 * t527 + t18 * t447 + t19 * t441) * pkin(5) + t741; -t137 * (mrSges(7,1) * t120 + mrSges(7,2) * t499) + (Ifges(7,1) * t499 - t578) * t639 + t64 * t638 + (Ifges(7,5) * t499 - Ifges(7,6) * t120) * t614 - t24 * t104 + t25 * t105 - g(1) * t546 - g(2) * t547 - t488 * t593 + (t120 * t25 + t24 * t499) * mrSges(7,3) + t474 + (-Ifges(7,2) * t120 + t111 + t65) * t641;];
tau  = t1;
