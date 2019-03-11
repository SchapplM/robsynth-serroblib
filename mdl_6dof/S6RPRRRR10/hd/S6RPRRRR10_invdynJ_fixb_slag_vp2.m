% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:47
% EndTime: 2019-03-09 07:29:09
% DurationCPUTime: 49.30s
% Computational Cost: add. (49807->997), mult. (157351->1383), div. (0->0), fcn. (136915->18), ass. (0->460)
t415 = sin(pkin(13));
t417 = sin(pkin(6));
t418 = cos(pkin(13));
t428 = cos(qJ(3));
t419 = cos(pkin(7));
t424 = sin(qJ(3));
t566 = t419 * t424;
t450 = (-t415 * t566 + t418 * t428) * t417;
t351 = qJD(1) * t450;
t416 = sin(pkin(7));
t553 = qJD(3) * t428;
t524 = t416 * t553;
t748 = t351 - t524;
t569 = t417 * t418;
t401 = qJ(2) * t569;
t420 = cos(pkin(6));
t625 = pkin(1) * t420;
t544 = qJD(1) * t625;
t362 = qJD(1) * t401 + t415 * t544;
t568 = t417 * t419;
t572 = t416 * t420;
t467 = t418 * t568 + t572;
t448 = t467 * qJD(1);
t301 = pkin(9) * t448 + t362;
t623 = pkin(9) * t415;
t353 = (-pkin(2) * t418 - t416 * t623 - pkin(1)) * t417;
t344 = qJD(1) * t353 + qJD(2);
t399 = t418 * t544;
t574 = t415 * t417;
t624 = pkin(2) * t420;
t446 = t624 + (-pkin(9) * t419 - qJ(2)) * t574;
t311 = qJD(1) * t446 + t399;
t585 = t311 * t419;
t191 = -t424 * t301 + t428 * (t344 * t416 + t585);
t555 = qJD(1) * t417;
t527 = t415 * t555;
t388 = t424 * t527;
t315 = t428 * t448 - t388;
t468 = t415 * t428 + t418 * t566;
t571 = t416 * t424;
t334 = t417 * t468 + t420 * t571;
t318 = t334 * qJD(1);
t245 = pkin(3) * t318 - pkin(10) * t315;
t423 = sin(qJ(4));
t427 = cos(qJ(4));
t138 = -t191 * t423 + t427 * t245;
t430 = -pkin(11) - pkin(10);
t530 = qJD(4) * t430;
t583 = t315 * t427;
t747 = -pkin(4) * t318 + pkin(11) * t583 + t427 * t530 - t138;
t139 = t427 * t191 + t423 * t245;
t584 = t315 * t423;
t746 = -pkin(11) * t584 - t423 * t530 + t139;
t374 = t419 * t423 + t427 * t571;
t507 = t416 * t527;
t716 = -qJD(4) * t374 + t423 * t748 - t427 * t507;
t373 = t419 * t427 - t423 * t571;
t742 = -qJD(4) * t373 + t423 * t507 + t427 * t748;
t552 = qJD(4) * t423;
t745 = t552 - t584;
t626 = sin(qJ(1));
t528 = t626 * t418;
t429 = cos(qJ(1));
t563 = t429 * t415;
t367 = t420 * t563 + t528;
t406 = t626 * t415;
t562 = t429 * t418;
t501 = t420 * t562 - t406;
t567 = t417 * t429;
t739 = -t416 * t567 + t419 * t501;
t279 = t367 * t428 + t424 * t739;
t338 = t501 * t416 + t419 * t567;
t580 = t338 * t423;
t744 = -t279 * t427 + t580;
t725 = mrSges(6,2) - mrSges(7,3);
t421 = sin(qJ(6));
t425 = cos(qJ(6));
t496 = -mrSges(7,1) * t425 + mrSges(7,2) * t421;
t710 = mrSges(6,1) - t496;
t422 = sin(qJ(5));
t426 = cos(qJ(5));
t394 = t430 * t423;
t395 = t430 * t427;
t474 = t426 * t394 + t395 * t422;
t697 = qJD(5) * t474 + t422 * t747 - t746 * t426;
t475 = t426 * t373 - t374 * t422;
t718 = qJD(5) * t475 + t422 * t716 - t426 * t742;
t565 = t419 * t428;
t451 = (t415 * t565 + t418 * t424) * t417;
t350 = qJD(1) * t451;
t554 = qJD(3) * t424;
t525 = t416 * t554;
t741 = -t350 + t525;
t469 = -t415 * t424 + t418 * t565;
t547 = qJD(1) * qJD(3);
t521 = t428 * t547;
t243 = (qJDD(1) * t424 + t521) * t572 + (qJDD(1) * t468 + t469 * t547) * t417;
t366 = -t416 * t569 + t419 * t420;
t358 = qJD(1) * t366 + qJD(3);
t253 = -t318 * t423 + t358 * t427;
t357 = qJDD(1) * t366 + qJDD(3);
t155 = qJD(4) * t253 + t243 * t427 + t357 * t423;
t254 = t318 * t427 + t358 * t423;
t156 = -qJD(4) * t254 - t243 * t423 + t357 * t427;
t515 = t426 * t253 - t254 * t422;
t86 = qJD(5) * t515 + t155 * t426 + t156 * t422;
t481 = t253 * t422 + t426 * t254;
t87 = -qJD(5) * t481 - t155 * t422 + t156 * t426;
t377 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t417;
t540 = qJDD(1) * t625;
t345 = -t377 * t415 + t418 * t540;
t291 = (-t568 * t623 + t624) * qJDD(1) + t345;
t340 = qJDD(1) * t353 + qJDD(2);
t346 = t418 * t377 + t415 * t540;
t447 = t467 * qJDD(1);
t717 = pkin(9) * t447 + qJD(3) * t585 + t346;
t124 = t428 * (t291 * t419 + t340 * t416) - t301 * t553 - t344 * t525 - t717 * t424;
t118 = -pkin(3) * t357 - t124;
t89 = -pkin(4) * t156 + t118;
t22 = -pkin(5) * t87 - pkin(12) * t86 + t89;
t449 = t467 * t428;
t310 = -qJD(1) * t449 + qJD(4) + t388;
t304 = qJD(5) + t310;
t246 = -t311 * t416 + t419 * t344;
t166 = -pkin(3) * t315 - pkin(10) * t318 + t246;
t192 = t428 * t301 + t311 * t566 + t344 * t571;
t168 = pkin(10) * t358 + t192;
t115 = t166 * t423 + t168 * t427;
t98 = pkin(11) * t253 + t115;
t592 = t426 * t98;
t114 = t427 * t166 - t168 * t423;
t97 = -pkin(11) * t254 + t114;
t94 = pkin(4) * t310 + t97;
t49 = t422 * t94 + t592;
t47 = pkin(12) * t304 + t49;
t167 = -pkin(3) * t358 - t191;
t135 = -pkin(4) * t253 + t167;
t88 = -pkin(5) * t515 - pkin(12) * t481 + t135;
t23 = -t421 * t47 + t425 * t88;
t546 = qJDD(1) * t417;
t520 = t415 * t546;
t511 = -t521 * t574 + (-t467 * t547 - t520) * t424;
t241 = -qJDD(1) * t449 + qJDD(4) - t511;
t123 = t291 * t566 - t301 * t554 + t340 * t571 + t344 * t524 + t428 * t717;
t117 = pkin(10) * t357 + t123;
t233 = -t291 * t416 + t419 * t340;
t244 = t428 * t447 + t511;
t133 = -pkin(3) * t244 - pkin(10) * t243 + t233;
t44 = -qJD(4) * t115 - t117 * t423 + t427 * t133;
t32 = pkin(4) * t241 - pkin(11) * t155 + t44;
t551 = qJD(4) * t427;
t43 = t427 * t117 + t423 * t133 + t166 * t551 - t168 * t552;
t38 = pkin(11) * t156 + t43;
t549 = qJD(5) * t426;
t550 = qJD(5) * t422;
t10 = t422 * t32 + t426 * t38 + t94 * t549 - t550 * t98;
t236 = qJDD(5) + t241;
t7 = pkin(12) * t236 + t10;
t2 = qJD(6) * t23 + t22 * t421 + t425 * t7;
t24 = t421 * t88 + t425 * t47;
t3 = -qJD(6) * t24 + t22 * t425 - t421 * t7;
t488 = t23 * t425 + t24 * t421;
t740 = -qJD(6) * t488 + t2 * t425 - t3 * t421;
t686 = pkin(4) * t745 - t192;
t169 = Ifges(6,4) * t515;
t604 = Ifges(6,5) * t304;
t611 = Ifges(6,1) * t481;
t113 = t169 + t604 + t611;
t637 = t304 / 0.2e1;
t643 = t481 / 0.2e1;
t644 = t515 / 0.2e1;
t595 = t422 * t98;
t48 = t426 * t94 - t595;
t683 = -t135 * mrSges(6,2) + t48 * mrSges(6,3);
t738 = -t683 + Ifges(6,1) * t643 + Ifges(6,4) * t644 + Ifges(6,5) * t637 + t113 / 0.2e1;
t414 = qJ(4) + qJ(5);
t411 = sin(t414);
t412 = cos(t414);
t219 = t279 * t412 - t338 * t411;
t218 = -t279 * t411 - t338 * t412;
t737 = t279 * t423 + t338 * t427;
t309 = Ifges(4,4) * t315;
t736 = Ifges(4,2) * t315;
t735 = t253 * Ifges(5,6);
t734 = t310 * Ifges(5,3);
t733 = m(7) * pkin(5) + t710;
t513 = -m(7) * pkin(12) + t725;
t732 = -pkin(12) * t318 + t697;
t379 = t422 * t427 + t423 * t426;
t231 = t379 * t315;
t378 = t422 * t423 - t426 * t427;
t232 = t378 * t315;
t684 = qJD(4) + qJD(5);
t336 = t684 * t378;
t337 = t684 * t379;
t731 = t686 + (-t232 + t336) * pkin(12) + (-t231 + t337) * pkin(5);
t601 = Ifges(6,6) * t304;
t602 = Ifges(6,2) * t515;
t608 = Ifges(6,4) * t481;
t112 = t601 + t602 + t608;
t170 = qJD(6) - t515;
t645 = t170 / 0.2e1;
t144 = t304 * t421 + t425 * t481;
t652 = t144 / 0.2e1;
t143 = t304 * t425 - t421 * t481;
t654 = t143 / 0.2e1;
t702 = t24 * mrSges(7,2);
t703 = t23 * mrSges(7,1);
t677 = -t135 * mrSges(6,1) + t49 * mrSges(6,3) + t702 - t703;
t599 = Ifges(7,3) * t170;
t600 = Ifges(7,6) * t143;
t603 = Ifges(7,5) * t144;
t76 = t599 + t600 + t603;
t730 = -t677 - Ifges(6,4) * t643 + Ifges(7,5) * t652 - Ifges(6,2) * t644 - Ifges(6,6) * t637 + Ifges(7,6) * t654 + Ifges(7,3) * t645 - t112 / 0.2e1 + t76 / 0.2e1;
t54 = qJD(6) * t143 + t236 * t421 + t425 * t86;
t55 = -qJD(6) * t144 + t236 * t425 - t421 * t86;
t85 = qJDD(6) - t87;
t14 = Ifges(7,5) * t54 + Ifges(7,6) * t55 + Ifges(7,3) * t85;
t642 = t236 / 0.2e1;
t662 = t87 / 0.2e1;
t663 = t86 / 0.2e1;
t664 = t85 / 0.2e1;
t668 = t55 / 0.2e1;
t669 = t54 / 0.2e1;
t679 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t729 = t679 + mrSges(6,1) * t89 - mrSges(6,3) * t10 + Ifges(7,5) * t669 + Ifges(7,6) * t668 + Ifges(7,3) * t664 + t14 / 0.2e1 + (-t642 - t236 / 0.2e1) * Ifges(6,6) + (-t662 - t87 / 0.2e1) * Ifges(6,2) + (-t663 - t86 / 0.2e1) * Ifges(6,4);
t728 = t124 * mrSges(4,1) - t123 * mrSges(4,2) + Ifges(4,5) * t243 + Ifges(4,6) * t244 + Ifges(4,3) * t357;
t276 = -t334 * t423 + t366 * t427;
t503 = pkin(4) * t276;
t11 = -qJD(5) * t49 + t32 * t426 - t38 * t422;
t726 = t11 * mrSges(6,1);
t724 = t114 * mrSges(5,1);
t723 = t115 * mrSges(5,2);
t349 = t394 * t422 - t395 * t426;
t696 = -qJD(5) * t349 + t746 * t422 + t426 * t747;
t147 = mrSges(6,1) * t304 - mrSges(6,3) * t481;
t96 = -mrSges(7,1) * t143 + mrSges(7,2) * t144;
t590 = t96 - t147;
t722 = t254 * Ifges(5,5) + Ifges(6,5) * t481 + Ifges(6,6) * t515 + t304 * Ifges(6,3) + t734 + t735;
t303 = t373 * t422 + t374 * t426;
t570 = t416 * t428;
t470 = -t303 * t425 + t421 * t570;
t721 = qJD(6) * t470 - t421 * t718 + t425 * t741;
t284 = -t303 * t421 - t425 * t570;
t720 = qJD(6) * t284 + t421 * t741 + t425 * t718;
t719 = -qJD(5) * t303 + t422 * t742 + t426 * t716;
t370 = t415 * t625 + t401;
t324 = pkin(9) * t467 + t370;
t405 = t418 * t625;
t335 = t405 + t446;
t477 = t335 * t419 + t353 * t416;
t199 = -t424 * t324 + t477 * t428;
t16 = t54 * Ifges(7,1) + t55 * Ifges(7,4) + t85 * Ifges(7,5);
t489 = Ifges(7,5) * t425 - Ifges(7,6) * t421;
t460 = t170 * t489;
t606 = Ifges(7,4) * t421;
t493 = Ifges(7,1) * t425 - t606;
t462 = t144 * t493;
t605 = Ifges(7,4) * t425;
t491 = -Ifges(7,2) * t421 + t605;
t463 = t143 * t491;
t46 = -pkin(5) * t304 - t48;
t495 = mrSges(7,1) * t421 + mrSges(7,2) * t425;
t466 = t46 * t495;
t548 = qJD(6) * t421;
t519 = -t548 / 0.2e1;
t142 = Ifges(7,4) * t143;
t78 = t144 * Ifges(7,1) + t170 * Ifges(7,5) + t142;
t593 = t425 * t78;
t534 = t593 / 0.2e1;
t545 = Ifges(6,5) * t86 + Ifges(6,6) * t87 + Ifges(6,3) * t236;
t627 = t421 / 0.2e1;
t15 = t54 * Ifges(7,4) + t55 * Ifges(7,2) + t85 * Ifges(7,6);
t673 = t15 / 0.2e1;
t607 = Ifges(7,4) * t144;
t77 = Ifges(7,2) * t143 + Ifges(7,6) * t170 + t607;
t8 = -pkin(5) * t236 - t11;
t714 = t726 - t10 * mrSges(6,2) + t16 * t627 + t425 * t673 + t8 * t496 + (Ifges(7,1) * t421 + t605) * t669 + (Ifges(7,2) * t425 + t606) * t668 + t545 + t77 * t519 + (Ifges(7,5) * t421 + Ifges(7,6) * t425) * t664 + (t466 + t534) * qJD(6) + (t463 + t462 + t460) * qJD(6) / 0.2e1 + t740 * mrSges(7,3);
t368 = -t406 * t420 + t562;
t455 = t420 * t528 + t563;
t529 = t417 * t626;
t687 = -t416 * t529 + t455 * t419;
t283 = t368 * t428 - t424 * t687;
t713 = t455 * t416 + t419 * t529;
t226 = -t283 * t423 + t427 * t713;
t712 = t417 * t469 + t420 * t570;
t280 = -t367 * t424 + t428 * t739;
t709 = t44 * mrSges(5,1) - t43 * mrSges(5,2);
t707 = -t11 * mrSges(6,3) + 0.2e1 * Ifges(6,1) * t663 + 0.2e1 * Ifges(6,4) * t662 + 0.2e1 * Ifges(6,5) * t642;
t704 = -m(6) - m(7);
t649 = t155 / 0.2e1;
t648 = t156 / 0.2e1;
t641 = t241 / 0.2e1;
t21 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t68 = mrSges(6,1) * t236 - mrSges(6,3) * t86;
t701 = t21 - t68;
t410 = pkin(4) * t427 + pkin(3);
t325 = pkin(5) * t378 - pkin(12) * t379 - t410;
t250 = t325 * t425 - t349 * t421;
t699 = qJD(6) * t250 + t421 * t731 + t425 * t732;
t251 = t325 * t421 + t349 * t425;
t698 = -qJD(6) * t251 - t421 * t732 + t425 * t731;
t695 = pkin(5) * t318 - t696;
t694 = m(7) * t8 + t21;
t256 = -t335 * t416 + t419 * t353;
t182 = -pkin(3) * t712 - pkin(10) * t334 + t256;
t308 = t428 * t324;
t200 = t335 * t566 + t353 * t571 + t308;
t190 = pkin(10) * t366 + t200;
t125 = t427 * t182 - t190 * t423;
t277 = t334 * t427 + t366 * t423;
t101 = -pkin(4) * t712 - pkin(11) * t277 + t125;
t126 = t423 * t182 + t427 * t190;
t110 = pkin(11) * t276 + t126;
t691 = t422 * t101 + t426 * t110;
t183 = t232 * t421 + t318 * t425;
t575 = t379 * t425;
t465 = qJD(6) * t575 - t421 * t336;
t690 = t183 + t465;
t184 = -t232 * t425 + t318 * t421;
t464 = t425 * t336 + t379 * t548;
t689 = t184 + t464;
t612 = mrSges(4,3) * t318;
t558 = -mrSges(4,1) * t358 - mrSges(5,1) * t253 + mrSges(5,2) * t254 + t612;
t194 = t276 * t422 + t277 * t426;
t322 = t712 * t425;
t160 = -t194 * t421 - t322;
t657 = -t113 / 0.2e1;
t688 = t488 * mrSges(7,3) + t657 - t604 / 0.2e1 - t463 / 0.2e1 - t462 / 0.2e1 - t460 / 0.2e1 + t77 * t627 - t593 / 0.2e1 - t466 + t683 - t169 / 0.2e1;
t500 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t681 = -t495 + t500;
t685 = -t423 * t44 + t427 * t43;
t682 = t415 * ((-mrSges(4,1) * t315 + mrSges(4,2) * t318) * t416 - (mrSges(3,1) * t420 - mrSges(3,3) * t574) * qJD(1)) + (-mrSges(3,2) * t420 + mrSges(3,3) * t569) * qJD(1) * t418;
t316 = t712 * qJD(3);
t209 = qJD(4) * t276 + t316 * t427;
t317 = t334 * qJD(3);
t176 = qJD(2) * t450 + qJD(3) * t199;
t506 = qJD(2) * t416 * t574;
t230 = pkin(3) * t317 - pkin(10) * t316 + t506;
t82 = -qJD(4) * t126 - t176 * t423 + t427 * t230;
t63 = pkin(4) * t317 - pkin(11) * t209 + t82;
t208 = -qJD(4) * t277 - t316 * t423;
t81 = t427 * t176 + t182 * t551 - t190 * t552 + t423 * t230;
t65 = pkin(11) * t208 + t81;
t20 = -qJD(5) * t691 - t422 * t65 + t426 * t63;
t498 = -mrSges(5,1) * t427 + mrSges(5,2) * t423;
t680 = m(5) * pkin(3) - t513 * t411 + t412 * t733 + mrSges(4,1) - t498;
t122 = pkin(5) * t481 - pkin(12) * t515;
t107 = -mrSges(7,2) * t170 + mrSges(7,3) * t143;
t108 = mrSges(7,1) * t170 - mrSges(7,3) * t144;
t25 = mrSges(7,1) * t85 - mrSges(7,3) * t54;
t26 = -mrSges(7,2) * t85 + mrSges(7,3) * t55;
t678 = (-t107 * t421 - t108 * t425) * qJD(6) + m(7) * t740 - t421 * t25 + t425 * t26;
t676 = 0.2e1 * t420;
t672 = t16 / 0.2e1;
t667 = -t76 / 0.2e1;
t661 = Ifges(5,4) * t649 + Ifges(5,2) * t648 + Ifges(5,6) * t641;
t660 = Ifges(5,1) * t649 + Ifges(5,4) * t648 + Ifges(5,5) * t641;
t658 = t112 / 0.2e1;
t655 = -t143 / 0.2e1;
t653 = -t144 / 0.2e1;
t598 = t254 * Ifges(5,4);
t153 = t253 * Ifges(5,2) + t310 * Ifges(5,6) + t598;
t651 = t153 / 0.2e1;
t252 = Ifges(5,4) * t253;
t154 = t254 * Ifges(5,1) + t310 * Ifges(5,5) + t252;
t650 = t154 / 0.2e1;
t646 = -t170 / 0.2e1;
t640 = -t253 / 0.2e1;
t639 = -t254 / 0.2e1;
t638 = t254 / 0.2e1;
t636 = -t310 / 0.2e1;
t631 = t318 / 0.2e1;
t629 = -t358 / 0.2e1;
t616 = mrSges(3,1) * t418;
t613 = mrSges(4,3) * t315;
t610 = Ifges(5,4) * t423;
t609 = Ifges(5,4) * t427;
t597 = t318 * Ifges(4,4);
t596 = t358 * Ifges(4,5);
t582 = t712 * t421;
t576 = t379 * t421;
t556 = t429 * pkin(1) + qJ(2) * t529;
t533 = Ifges(5,5) * t155 + Ifges(5,6) * t156 + Ifges(5,3) * t241;
t121 = -mrSges(6,1) * t515 + mrSges(6,2) * t481;
t532 = -t121 - t558;
t222 = t283 * t411 - t412 * t713;
t223 = t283 * t412 + t411 * t713;
t518 = -t222 * pkin(5) + pkin(12) * t223;
t505 = t737 * pkin(4);
t504 = t226 * pkin(4);
t502 = -pkin(1) * t626 + qJ(2) * t567;
t499 = -mrSges(5,1) * t276 + mrSges(5,2) * t277;
t494 = Ifges(5,1) * t427 - t610;
t492 = -Ifges(5,2) * t423 + t609;
t490 = Ifges(5,5) * t427 - Ifges(5,6) * t423;
t59 = -pkin(12) * t712 + t691;
t189 = -pkin(3) * t366 - t199;
t145 = t189 - t503;
t480 = t426 * t276 - t277 * t422;
t92 = -pkin(5) * t480 - pkin(12) * t194 + t145;
t28 = t421 * t92 + t425 * t59;
t27 = -t421 * t59 + t425 * t92;
t60 = t101 * t426 - t110 * t422;
t161 = t194 * t425 - t582;
t476 = -(-qJ(2) * t527 + t399) * t415 + t362 * t418;
t19 = t101 * t549 - t110 * t550 + t422 * t63 + t426 * t65;
t454 = -t218 * t710 + t219 * t725;
t453 = t222 * t710 + t223 * t725;
t267 = -t334 * t411 + t366 * t412;
t268 = t334 * t412 + t366 * t411;
t452 = -t267 * t710 + t268 * t725;
t441 = -t367 * pkin(2) + pkin(9) * t338 + t502;
t438 = t368 * pkin(2) + pkin(9) * t713 + t556;
t437 = t713 * t423;
t282 = t368 * t424 + t428 * t687;
t433 = pkin(4) * t437 - t282 * t430 + t283 * t410 + t438;
t177 = qJD(2) * t451 + (t424 * t477 + t308) * qJD(3);
t134 = -pkin(4) * t208 + t177;
t432 = t601 / 0.2e1 - t599 / 0.2e1 - t600 / 0.2e1 - t603 / 0.2e1 + t667 + t658 + t608 / 0.2e1 + t677;
t400 = -pkin(1) * t546 + qJDD(2);
t391 = mrSges(3,2) * t520;
t369 = -qJ(2) * t574 + t405;
t265 = -mrSges(4,2) * t358 + t613;
t264 = t267 * pkin(5);
t227 = t283 * t427 + t437;
t216 = t218 * pkin(5);
t207 = t318 * Ifges(4,1) + t309 + t596;
t206 = t358 * Ifges(4,6) + t597 + t736;
t202 = -mrSges(4,2) * t357 + mrSges(4,3) * t244;
t201 = mrSges(4,1) * t357 - mrSges(4,3) * t243;
t198 = mrSges(5,1) * t310 - mrSges(5,3) * t254;
t197 = -mrSges(5,2) * t310 + mrSges(5,3) * t253;
t162 = -mrSges(4,1) * t244 + mrSges(4,2) * t243;
t159 = t223 * t425 + t282 * t421;
t158 = -t223 * t421 + t282 * t425;
t146 = -mrSges(6,2) * t304 + mrSges(6,3) * t515;
t128 = -mrSges(5,2) * t241 + mrSges(5,3) * t156;
t127 = mrSges(5,1) * t241 - mrSges(5,3) * t155;
t106 = qJD(5) * t194 - t426 * t208 + t209 * t422;
t105 = qJD(5) * t480 + t208 * t422 + t209 * t426;
t103 = pkin(4) * t254 + t122;
t102 = -mrSges(5,1) * t156 + mrSges(5,2) * t155;
t73 = -qJD(6) * t161 - t105 * t421 + t317 * t425;
t72 = qJD(6) * t160 + t105 * t425 + t317 * t421;
t69 = -mrSges(6,2) * t236 + mrSges(6,3) * t87;
t58 = pkin(5) * t712 - t60;
t57 = t426 * t97 - t595;
t56 = t422 * t97 + t592;
t42 = pkin(5) * t106 - pkin(12) * t105 + t134;
t39 = -mrSges(6,1) * t87 + mrSges(6,2) * t86;
t34 = t122 * t421 + t425 * t48;
t33 = t122 * t425 - t421 * t48;
t31 = t103 * t421 + t425 * t57;
t30 = t103 * t425 - t421 * t57;
t18 = -pkin(5) * t317 - t20;
t17 = pkin(12) * t317 + t19;
t5 = -qJD(6) * t28 - t17 * t421 + t42 * t425;
t4 = qJD(6) * t27 + t17 * t425 + t42 * t421;
t1 = [(Ifges(7,5) * t161 + Ifges(7,6) * t160) * t664 + (Ifges(7,5) * t72 + Ifges(7,6) * t73) * t645 + (mrSges(4,2) * t233 - mrSges(4,3) * t124 + Ifges(4,1) * t243 + Ifges(4,4) * t244 + Ifges(4,5) * t357) * t334 + (Ifges(2,3) + (mrSges(3,1) * t369 - mrSges(3,2) * t370 + Ifges(3,3) * t420) * t420 + ((-mrSges(3,3) * t369 + Ifges(3,1) * t574 + Ifges(3,5) * t676) * t415 + (t370 * mrSges(3,3) + Ifges(3,6) * t676 + (mrSges(3,1) * pkin(1) + 0.2e1 * Ifges(3,4) * t415 + Ifges(3,2) * t418) * t417) * t418) * t417) * qJDD(1) + (t596 / 0.2e1 + t207 / 0.2e1 - t191 * mrSges(4,3) + t246 * mrSges(4,2) + t309 / 0.2e1 + Ifges(4,1) * t631) * t316 + t161 * t672 + t160 * t673 + m(6) * (t10 * t691 + t11 * t60 + t134 * t135 + t145 * t89 + t19 * t49 + t20 * t48) + t691 * t69 + (Ifges(5,1) * t277 + Ifges(5,4) * t276) * t649 + (Ifges(5,5) * t277 + Ifges(5,6) * t276) * t641 + t707 * t194 + (Ifges(7,1) * t161 + Ifges(7,4) * t160) * t669 + (Ifges(7,1) * t72 + Ifges(7,4) * t73) * t652 + (-m(4) * t441 + t279 * mrSges(4,1) - t338 * mrSges(4,3) - m(3) * t502 + t367 * mrSges(3,1) + t501 * mrSges(3,2) - mrSges(3,3) * t567 - m(5) * (-pkin(3) * t279 + t441) - t744 * mrSges(5,1) - t737 * mrSges(5,2) + t626 * mrSges(2,1) + t429 * mrSges(2,2) + t513 * t218 + t733 * t219 + t681 * t280 + t704 * (pkin(4) * t580 - t279 * t410 - t280 * t430 + t441)) * g(1) + (Ifges(5,1) * t209 + Ifges(5,4) * t208) * t638 + (Ifges(5,4) * t277 + Ifges(5,2) * t276) * t648 + (t724 + t722 / 0.2e1 + t246 * mrSges(4,1) + Ifges(6,6) * t644 + Ifges(6,3) * t637 + Ifges(6,5) * t643 - t192 * mrSges(4,3) - t206 / 0.2e1 + t629 * Ifges(4,6) + t734 / 0.2e1 + t735 / 0.2e1 + t48 * mrSges(6,1) - t736 / 0.2e1 - t49 * mrSges(6,2) - t723 + Ifges(5,5) * t638 - Ifges(4,4) * t631) * t317 + (mrSges(3,1) * t345 - mrSges(3,2) * t346) * t420 + (t160 * t2 - t161 * t3 - t23 * t72 + t24 * t73) * mrSges(7,3) + t310 * (Ifges(5,5) * t209 + Ifges(5,6) * t208) / 0.2e1 + t558 * t177 + m(5) * (t114 * t82 + t115 * t81 + t118 * t189 + t125 * t44 + t126 * t43 + t167 * t177) + m(7) * (t18 * t46 + t2 * t28 + t23 * t5 + t24 * t4 + t27 * t3 + t58 * t8) + (-pkin(1) * t391 + t400 * (mrSges(3,2) * t415 - t616) + (-t345 * t415 + t346 * t418) * mrSges(3,3) + t682 * qJD(2)) * t417 + t728 * t366 + (Ifges(7,4) * t161 + Ifges(7,2) * t160) * t668 + (Ifges(7,4) * t72 + Ifges(7,2) * t73) * t654 + t253 * (Ifges(5,4) * t209 + Ifges(5,2) * t208) / 0.2e1 - (t545 + t533) * t712 / 0.2e1 - (mrSges(4,1) * t233 - mrSges(4,3) * t123 - Ifges(4,4) * t243 + Ifges(5,5) * t649 + Ifges(6,5) * t663 - Ifges(4,2) * t244 - Ifges(4,6) * t357 + Ifges(5,6) * t648 + Ifges(6,6) * t662 + Ifges(5,3) * t641 + Ifges(6,3) * t642 + t709 + t726) * t712 + (t10 * t712 + t194 * t89) * mrSges(6,2) + (-m(5) * (t283 * pkin(3) + t438) - t227 * mrSges(5,1) - t226 * mrSges(5,2) - m(6) * t433 - t223 * mrSges(6,1) - m(7) * (t223 * pkin(5) + t433) - t159 * mrSges(7,1) - t158 * mrSges(7,2) - m(3) * t556 - t368 * mrSges(3,1) + t455 * mrSges(3,2) - mrSges(3,3) * t529 - m(4) * t438 - t283 * mrSges(4,1) - t713 * mrSges(4,3) - t429 * mrSges(2,1) + t626 * mrSges(2,2) + t513 * t222 + t500 * t282) * g(2) - t729 * t480 + t730 * t106 + (-t114 * t209 + t115 * t208 + t276 * t43 - t277 * t44) * mrSges(5,3) + t176 * t265 + m(4) * (t123 * t200 + t124 * t199 + t176 * t192 - t177 * t191 + t233 * t256 + t246 * t506) + t256 * t162 + t167 * (-mrSges(5,1) * t208 + mrSges(5,2) * t209) + t81 * t197 + t82 * t198 + t199 * t201 + t200 * t202 + t118 * t499 + t189 * t102 + t8 * (-mrSges(7,1) * t160 + mrSges(7,2) * t161) + t277 * t660 + t276 * t661 + t208 * t651 + t209 * t650 + t19 * t146 + t20 * t147 + m(3) * (t345 * t369 + t346 * t370 + (-pkin(1) * t400 + qJD(2) * t476) * t417) + t145 * t39 + t134 * t121 + t125 * t127 + t126 * t128 + t4 * t107 + t5 * t108 + t18 * t96 + t73 * t77 / 0.2e1 + t72 * t78 / 0.2e1 + t46 * (-mrSges(7,1) * t73 + mrSges(7,2) * t72) + t60 * t68 + t58 * t21 + t27 * t25 + t28 * t26 + t738 * t105; (-qJD(1) * t682 - qJDD(1) * t616) * t417 + t532 * t350 + (t424 * t202 + (-t102 + t201 - t39) * t428 + (t265 * t428 - t424 * t532) * qJD(3)) * t416 + t716 * t198 - t742 * t197 + t303 * t69 + t284 * t25 - t470 * t26 + t718 * t146 - t351 * t265 + t373 * t127 + t374 * t128 + t391 + t419 * t162 + t720 * t107 - t701 * t475 + t721 * t108 - t590 * t719 + ((-g(1) * t626 + g(2) * t429) * t417 - t420 * g(3)) * (m(3) + m(4) + m(5) - t704) + (-t2 * t470 + t23 * t721 + t24 * t720 + t284 * t3 - t46 * t719 - t475 * t8) * m(7) + (t10 * t303 + t11 * t475 + (t135 * t554 - t428 * t89) * t416 - t135 * t350 + t718 * t49 + t719 * t48) * m(6) + (-t167 * t350 + t373 * t44 + t374 * t43 + (-t118 * t428 + t167 * t554) * t416 - t742 * t115 + t716 * t114) * m(5) + (t191 * t350 - t192 * t351 - t246 * t507 + t233 * t419 + (t123 * t424 + t124 * t428 + (-t191 * t424 + t192 * t428) * qJD(3)) * t416) * m(4) + (-t476 * t555 + t400) * m(3); (t253 * t492 + t254 * t494 + t310 * t490) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t318 + t207 + t309) * t315 / 0.2e1 + t318 * t723 - (Ifges(4,1) * t315 - t597 + t722) * t318 / 0.2e1 + t575 * t672 + (t10 * t349 + t11 * t474 + t135 * t686 - t410 * t89 + t48 * t696 + t49 * t697) * m(6) + (t2 * t251 + t23 * t698 + t24 * t699 + t250 * t3 + t46 * t695 - t474 * t8) * m(7) - t701 * t474 - t515 * (-Ifges(6,4) * t232 - Ifges(6,2) * t231 + Ifges(6,6) * t318) / 0.2e1 - t481 * (-Ifges(6,1) * t232 - Ifges(6,4) * t231 + Ifges(6,5) * t318) / 0.2e1 - t304 * (-Ifges(6,5) * t232 - Ifges(6,6) * t231 + Ifges(6,3) * t318) / 0.2e1 - t48 * (mrSges(6,1) * t318 + mrSges(6,3) * t232) + (t89 * mrSges(6,2) + t489 * t664 + t491 * t668 + t493 * t669 + t495 * t8 + t519 * t78 + t707) * t379 - t135 * (mrSges(6,1) * t231 - mrSges(6,2) * t232) + (-pkin(3) * t118 - t114 * t138 - t115 * t139) * m(5) + (-t265 + t613) * t191 + (-t315 + qJD(4)) * t167 * (mrSges(5,1) * t423 + mrSges(5,2) * t427) + (-Ifges(7,1) * t464 - Ifges(7,4) * t465) * t652 - t15 * t576 / 0.2e1 - t49 * (-mrSges(6,2) * t318 - mrSges(6,3) * t231) - t246 * (mrSges(4,1) * t318 + mrSges(4,2) * t315) - t153 * t552 / 0.2e1 - t154 * t583 / 0.2e1 + t728 + (t704 * (-t279 * t430 + t280 * t410) + t681 * t279 - t680 * t280) * g(2) + (t704 * (-t334 * t430 + t410 * t712) + t681 * t334 - t680 * t712) * g(3) + t729 * t378 + t730 * t337 - t318 * t724 + t250 * t25 + t251 * t26 - t139 * t197 - t138 * t198 + t118 * t498 - t184 * t78 / 0.2e1 + t423 * t660 + t427 * t661 + t231 * t667 + t551 * t650 + t584 * t651 + (Ifges(7,1) * t184 + Ifges(7,4) * t183 + Ifges(7,5) * t231) * t653 + (Ifges(7,4) * t184 + Ifges(7,2) * t183 + Ifges(7,6) * t231) * t655 - t232 * t657 + t231 * t658 + (Ifges(5,6) * t318 + t315 * t492) * t640 + (Ifges(5,5) * t423 + Ifges(5,6) * t427) * t641 + (Ifges(7,5) * t184 + Ifges(7,6) * t183 + Ifges(7,3) * t231) * t646 + (Ifges(5,2) * t427 + t610) * t648 + (Ifges(5,1) * t423 + t609) * t649 + (Ifges(5,3) * t318 + t315 * t490) * t636 + (Ifges(5,5) * t318 + t315 * t494) * t639 + (Ifges(4,5) * t315 - Ifges(4,6) * t318) * t629 + t206 * t631 + (-Ifges(7,5) * t464 - Ifges(7,6) * t465) * t645 - pkin(3) * t102 + (-Ifges(7,4) * t464 - Ifges(7,2) * t465) * t654 + (-t198 * t551 - t197 * t552 + m(5) * ((-t114 * t427 - t115 * t423) * qJD(4) + t685) - t423 * t127 + t427 * t128) * pkin(10) + t686 * t121 + (-m(5) * t167 - t558 + t612) * t192 + (-t745 * t115 + (-t551 + t583) * t114 + t685) * mrSges(5,3) - t690 * t77 / 0.2e1 + (mrSges(7,1) * t690 - mrSges(7,2) * t689) * t46 + (-t2 * t576 + t23 * t689 - t24 * t690 - t3 * t575) * mrSges(7,3) + t349 * t69 + t695 * t96 + t696 * t147 + t697 * t146 + t698 * t108 + t699 * t107 - t410 * t39 + t231 * t702 - (t534 + t738) * t336 - t231 * t703 + (t704 * (-t282 * t410 - t283 * t430) + t681 * t283 + t680 * t282) * g(1); (m(6) * t505 - m(7) * (pkin(12) * t219 + t216 - t505) + t737 * mrSges(5,1) - t744 * mrSges(5,2) + t454) * g(2) + (t602 / 0.2e1 + t432) * t481 - m(6) * (-t48 * t56 + t49 * t57) + t678 * (pkin(4) * t422 + pkin(12)) - m(7) * (t23 * t30 + t24 * t31 + t46 * t56) + t709 + (-Ifges(5,2) * t254 + t154 + t252) * t640 + t714 + (t114 * t253 + t115 * t254) * mrSges(5,3) + (0.2e1 * t135 * t639 * m(6) - t121 * t254 + (m(6) * t11 + t68) * t426 + (m(6) * t10 + t69) * t422 + ((-m(6) * t48 + m(7) * t46 + t590) * t422 + (m(6) * t49 + m(7) * (-t23 * t421 + t24 * t425) + t107 * t425 - t108 * t421 + t146) * t426) * qJD(5)) * pkin(4) + t533 + (-mrSges(5,1) * t226 + mrSges(5,2) * t227 - m(6) * t504 - m(7) * (t504 + t518) + t453) * g(1) - t590 * t56 - t167 * (mrSges(5,1) * t254 + mrSges(5,2) * t253) + (-m(6) * t503 - m(7) * (pkin(12) * t268 + t264 + t503) + t452 + t499) * g(3) - t114 * t197 + t115 * t198 + (Ifges(5,5) * t253 - Ifges(5,6) * t254) * t636 + t153 * t638 + (Ifges(5,1) * t253 - t598) * t639 - t57 * t146 - t31 * t107 - t30 * t108 + (-t611 / 0.2e1 + t688) * t515 + t694 * (-pkin(4) * t426 - pkin(5)); t432 * t481 + (-m(7) * t216 + t454) * g(2) + (-m(7) * t518 + t453) * g(1) + (-m(7) * t264 + t452) * g(3) - m(7) * (t23 * t33 + t24 * t34 + t46 * t49) - t590 * t49 + ((-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t481 + t688) * t515 - t48 * t146 - t34 * t107 - t33 * t108 + ((-g(2) * t219 - g(3) * t268) * m(7) + t678) * pkin(12) - t694 * pkin(5) + t714; -t46 * (mrSges(7,1) * t144 + mrSges(7,2) * t143) + (Ifges(7,1) * t143 - t607) * t653 + t77 * t652 + (Ifges(7,5) * t143 - Ifges(7,6) * t144) * t646 - t23 * t107 + t24 * t108 - g(1) * (mrSges(7,1) * t158 - mrSges(7,2) * t159) - g(2) * ((-t219 * t421 - t280 * t425) * mrSges(7,1) + (-t219 * t425 + t280 * t421) * mrSges(7,2)) - g(3) * ((-t268 * t421 - t322) * mrSges(7,1) + (-t268 * t425 + t582) * mrSges(7,2)) + (t143 * t23 + t144 * t24) * mrSges(7,3) + t14 + (-Ifges(7,2) * t144 + t142 + t78) * t655 + t679;];
tau  = t1;
