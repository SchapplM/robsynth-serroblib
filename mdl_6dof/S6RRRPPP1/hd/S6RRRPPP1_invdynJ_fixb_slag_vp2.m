% Calculate vector of inverse dynamics joint torques for
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:30
% EndTime: 2019-03-09 15:16:13
% DurationCPUTime: 77.82s
% Computational Cost: add. (14460->947), mult. (34493->1214), div. (0->0), fcn. (25643->10), ass. (0->458)
t381 = sin(qJ(3));
t382 = sin(qJ(2));
t521 = qJD(1) * t382;
t488 = t381 * t521;
t384 = cos(qJ(3));
t519 = qJD(2) * t384;
t317 = -t488 + t519;
t385 = cos(qJ(2));
t509 = qJD(1) * qJD(2);
t331 = qJDD(1) * t382 + t385 * t509;
t199 = qJD(3) * t317 + qJDD(2) * t381 + t331 * t384;
t487 = t384 * t521;
t318 = qJD(2) * t381 + t487;
t200 = -qJD(3) * t318 + qJDD(2) * t384 - t331 * t381;
t330 = qJDD(1) * t385 - t382 * t509;
t310 = qJDD(3) - t330;
t377 = sin(pkin(10));
t379 = cos(pkin(10));
t380 = cos(pkin(6));
t546 = t379 * t380;
t378 = sin(pkin(6));
t551 = t378 * t379;
t99 = t199 * t377 - t200 * t546 - t310 * t551;
t616 = -t99 / 0.2e1;
t436 = t200 * t380 + t310 * t378;
t100 = t199 * t379 + t377 * t436;
t613 = t100 / 0.2e1;
t142 = -t200 * t378 + t310 * t380;
t609 = t142 / 0.2e1;
t460 = pkin(2) * t382 - pkin(9) * t385;
t321 = t460 * qJD(1);
t236 = pkin(8) * t488 + t384 * t321;
t542 = t380 * t385;
t501 = qJ(4) * t542;
t414 = pkin(3) * t382 - t384 * t501;
t181 = qJD(1) * t414 + t236;
t570 = qJ(4) * t380;
t468 = pkin(9) + t570;
t442 = qJD(3) * t468;
t513 = qJD(4) * t381;
t481 = t380 * t513;
t243 = -t384 * t442 - t481;
t734 = t181 - t243;
t510 = t385 * qJD(1);
t371 = pkin(8) * t510;
t549 = t378 * t384;
t583 = pkin(3) * t381;
t427 = -qJ(4) * t549 + t583;
t234 = t427 * t510 + t371;
t241 = qJD(3) * t427 - t378 * t513;
t733 = -t234 + t241;
t294 = t381 * t321;
t512 = qJD(4) * t384;
t480 = t380 * t512;
t537 = t382 * t384;
t539 = t381 * t385;
t550 = t378 * t382;
t700 = -t380 * t539 + t550;
t732 = t381 * t442 - t480 + t294 + (-pkin(8) * t537 + qJ(4) * t700) * qJD(1);
t466 = -qJD(3) + t510;
t557 = t317 * t378;
t229 = t380 * t466 + t557;
t569 = qJDD(1) * pkin(1);
t217 = -pkin(2) * t330 - pkin(9) * t331 - t569;
t307 = t330 * pkin(8);
t278 = qJDD(2) * pkin(9) + t307;
t372 = t382 * pkin(9);
t374 = t385 * pkin(2);
t490 = -pkin(1) - t374;
t432 = t490 - t372;
t299 = t432 * qJD(1);
t348 = qJD(2) * pkin(9) + t371;
t515 = qJD(3) * t384;
t517 = qJD(3) * t381;
t111 = t381 * t217 + t384 * t278 + t299 * t515 - t348 * t517;
t431 = t466 * t378;
t556 = t317 * t380;
t230 = -t431 + t556;
t47 = qJ(4) * t436 + qJD(4) * t230 + t111;
t552 = t377 * t380;
t553 = t377 * t378;
t212 = t381 * t299 + t384 * t348;
t112 = -qJD(3) * t212 + t384 * t217 - t278 * t381;
t433 = -qJ(4) * t199 - qJD(4) * t318;
t60 = pkin(3) * t310 + t380 * t433 + t112;
t308 = t331 * pkin(8);
t279 = -qJDD(2) * pkin(2) + t308;
t89 = -pkin(3) * t200 + t378 * t433 + t279;
t8 = t379 * t47 + t60 * t552 + t89 * t553;
t4 = -qJ(5) * t142 + qJD(5) * t229 - t8;
t3 = -pkin(5) * t99 + qJDD(6) - t4;
t614 = -t100 / 0.2e1;
t642 = -Ifges(6,6) + Ifges(7,6);
t643 = Ifges(7,2) + Ifges(6,3);
t645 = Ifges(7,4) + Ifges(6,5);
t730 = -t142 / 0.2e1;
t731 = -t4 * mrSges(6,1) + t3 * mrSges(7,1) + t8 * mrSges(5,3) + Ifges(5,4) * t613 + Ifges(5,6) * t609 + t642 * t614 + t645 * t730 + (Ifges(5,2) + t643) * t616;
t644 = Ifges(5,5) + Ifges(7,5);
t729 = -Ifges(6,4) + t644;
t728 = -Ifges(5,6) + t645;
t647 = Ifges(5,1) + Ifges(7,3);
t727 = t647 * t613;
t533 = t384 * t385;
t216 = t377 * t700 + t379 * t533;
t203 = t216 * qJD(1);
t267 = t379 * t515 - t517 * t552;
t672 = t203 - t267;
t16 = -t378 * t60 + t380 * t89 + qJDD(4);
t555 = t318 * t377;
t159 = -t317 * t546 + t379 * t431 + t555;
t547 = t379 * t318;
t566 = t230 * t377;
t160 = t547 + t566;
t390 = -qJ(5) * t100 - qJD(5) * t160 + t16;
t579 = pkin(4) + qJ(6);
t2 = qJD(6) * t159 + t579 * t99 + t390;
t5 = pkin(4) * t99 + t390;
t615 = t99 / 0.2e1;
t7 = -t377 * t47 + t379 * (t378 * t89 + t380 * t60);
t397 = qJDD(5) - t7;
t1 = pkin(5) * t100 + qJD(6) * t229 - t142 * t579 + t397;
t6 = -pkin(4) * t142 + t397;
t646 = Ifges(5,4) - Ifges(7,6);
t696 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t7 * mrSges(5,3) + Ifges(6,4) * t730 + Ifges(6,2) * t613 + t644 * t609 + t727 + (Ifges(6,6) + t646) * t616;
t726 = t16 * mrSges(5,2) - t2 * mrSges(7,2) - t5 * mrSges(6,3) + Ifges(5,4) * t616 - Ifges(6,2) * t614 + t609 * t729 + t642 * t615 + t696 + t727;
t725 = t16 * mrSges(5,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) - Ifges(5,2) * t616 + Ifges(6,6) * t614 + t728 * t609 - t646 * t613 + t643 * t615 - t731;
t675 = -t732 * t379 - t552 * t734 + t733 * t553;
t724 = -t379 * (t181 * t380 + t234 * t378) + t732 * t377;
t674 = t378 * t734 + t733 * t380;
t545 = t379 * t381;
t423 = t377 * t384 + t380 * t545;
t493 = t379 * t550;
t215 = t423 * t385 - t493;
t202 = t215 * qJD(1);
t266 = t423 * qJD(3);
t673 = t202 - t266;
t426 = qJ(4) * t378 * t381 + pkin(3) * t384;
t311 = -pkin(2) - t426;
t723 = t382 * t311 + t501;
t629 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t594 = t229 / 0.2e1;
t605 = t160 / 0.2e1;
t606 = -t160 / 0.2e1;
t607 = t159 / 0.2e1;
t608 = -t159 / 0.2e1;
t720 = Ifges(6,4) * t605 + Ifges(5,6) * t607 + t629 * t594 + t644 * t606 + t645 * t608;
t595 = -t229 / 0.2e1;
t719 = Ifges(6,4) * t606 + Ifges(5,6) * t608 + t629 * t595 + t644 * t605 + t645 * t607;
t72 = -t229 * Ifges(6,4) - t160 * Ifges(6,2) + t159 * Ifges(6,6);
t718 = Ifges(5,4) * t607 - Ifges(6,2) * t605 + t594 * t729 + t647 * t606 + t642 * t608 + t72 / 0.2e1;
t717 = Ifges(5,4) * t608 - Ifges(6,2) * t606 + t729 * t595 + t647 * t605 + t642 * t607 - t72 / 0.2e1;
t73 = t160 * Ifges(5,4) - t159 * Ifges(5,2) - t229 * Ifges(5,6);
t716 = -Ifges(5,2) * t607 + Ifges(6,6) * t605 + t594 * t728 - t646 * t606 + t643 * t608 + t73 / 0.2e1;
t715 = -Ifges(5,2) * t608 + Ifges(6,6) * t606 + t728 * t595 - t646 * t605 + t643 * t607 - t73 / 0.2e1;
t712 = -t7 * mrSges(5,1) + t8 * mrSges(5,2) - t6 * mrSges(6,2) - t3 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3) - t613 * t644 - t615 * t645 + (-t609 + t730) * t629 + (-Ifges(5,6) + t728) * t616 + (-Ifges(6,4) + t729) * t614;
t711 = -t317 / 0.2e1;
t710 = -t318 / 0.2e1;
t650 = t466 / 0.2e1;
t658 = t159 * t728 + t160 * t729 - t229 * t629;
t709 = -t658 / 0.2e1;
t682 = -t159 * t646 + t160 * t647 - t229 * t644;
t708 = -t682 / 0.2e1;
t683 = t159 * t643 + t160 * t642 - t229 * t645;
t707 = -t683 / 0.2e1;
t500 = t378 * t539;
t544 = t380 * t382;
t287 = t500 + t544;
t269 = t287 * qJD(1);
t678 = qJ(5) * t269 - t378 * (qJ(5) * t517 - qJD(5) * t384) - t675;
t703 = -t243 * t546 - t724;
t514 = qJD(4) * t378;
t211 = t384 * t299 - t381 * t348;
t161 = -t318 * t570 + t211;
t162 = -qJ(4) * t556 - t212;
t208 = pkin(3) * t318 - qJ(4) * t557;
t67 = t379 * t161 + t162 * t552 + t208 * t553;
t702 = t379 * t514 - t67;
t543 = t380 * t384;
t281 = t377 * t543 + t545;
t701 = qJ(5) * t672 - qJD(5) * t281 + t674;
t540 = t381 * t382;
t422 = t378 * t385 + t380 * t540;
t649 = -mrSges(6,2) + mrSges(5,1);
t699 = m(7) * qJ(6) + mrSges(7,3) + t649;
t286 = t378 * t540 - t542;
t453 = mrSges(4,1) * t384 - mrSges(4,2) * t381;
t454 = mrSges(3,1) * t382 + mrSges(3,2) * t385;
t648 = -mrSges(5,3) - mrSges(6,1);
t655 = m(7) * pkin(5) + mrSges(7,1) - t648;
t698 = t286 * t655 + t454 - t385 * mrSges(4,3) - (-m(4) * pkin(2) - t453) * t382;
t135 = qJ(4) * t230 + t212;
t141 = -pkin(3) * t466 + t161;
t370 = pkin(8) * t521;
t347 = -qJD(2) * pkin(2) + t370;
t554 = t318 * t378;
t182 = -pkin(3) * t317 - qJ(4) * t554 + t347;
t45 = -t377 * t135 + t379 * (t141 * t380 + t182 * t378);
t392 = qJD(5) - t45;
t14 = pkin(5) * t160 + t229 * t579 + t392;
t46 = t379 * t135 + t141 * t552 + t182 * t553;
t37 = qJ(5) * t229 - t46;
t27 = -pkin(5) * t159 + qJD(6) - t37;
t35 = pkin(4) * t229 + t392;
t695 = -t45 * mrSges(5,1) + t46 * mrSges(5,2) - t35 * mrSges(6,2) - t27 * mrSges(7,2) + t37 * mrSges(6,3) + t14 * mrSges(7,3) + t709;
t90 = -t141 * t378 + t380 * t182 + qJD(4);
t410 = -qJ(5) * t160 + t90;
t26 = t159 * t579 + t410;
t31 = pkin(4) * t159 + t410;
t694 = -t35 * mrSges(6,1) - t14 * mrSges(7,1) - t90 * mrSges(5,2) + t26 * mrSges(7,2) + t45 * mrSges(5,3) + t31 * mrSges(6,3) + t708;
t693 = -t90 * mrSges(5,1) - t37 * mrSges(6,1) + t27 * mrSges(7,1) + t31 * mrSges(6,2) + t46 * mrSges(5,3) - t26 * mrSges(7,3) + t707;
t691 = -m(7) - m(6);
t690 = t330 / 0.2e1;
t689 = t331 / 0.2e1;
t518 = qJD(2) * t385;
t688 = pkin(8) * t518;
t687 = -mrSges(3,3) + mrSges(2,2);
t496 = t379 * t543;
t280 = t377 * t381 - t496;
t684 = qJD(6) * t280 - t579 * t673 + t701;
t565 = t241 * t379;
t681 = t269 * t579 + (qJD(6) * t384 - t517 * t579 - t565) * t378 + t703 - t672 * pkin(5);
t680 = pkin(5) * t673 - t678;
t679 = -pkin(4) * t673 + t701;
t677 = pkin(4) * t269 + (-pkin(4) * t517 - t565) * t378 + t703;
t676 = (t241 * t378 + t243 * t380) * t379 + t724;
t639 = qJ(5) * t554 - qJD(5) * t380 - t702;
t486 = t378 * t517;
t671 = t269 - t486;
t304 = Ifges(4,4) * t317;
t193 = Ifges(4,1) * t318 - Ifges(4,5) * t466 + t304;
t369 = Ifges(3,4) * t510;
t670 = Ifges(3,1) * t521 + Ifges(3,5) * qJD(2) + t384 * t193 + t369;
t669 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t317 + mrSges(4,2) * t318 + mrSges(3,3) * t521;
t383 = sin(qJ(1));
t386 = cos(qJ(1));
t531 = t386 * t381;
t291 = t383 * t533 - t531;
t375 = t386 * pkin(8);
t535 = t383 * t385;
t290 = t381 * t535 + t384 * t386;
t564 = t290 * t378;
t668 = t383 * (-t382 * t468 + t490) - t291 * pkin(3) - qJ(4) * t564 + t375;
t532 = t385 * t386;
t363 = pkin(9) * t532;
t667 = t386 * t723 + t363;
t361 = pkin(9) * t535;
t666 = t383 * t723 + t361;
t455 = mrSges(3,1) * t385 - mrSges(3,2) * t382;
t664 = -t382 * mrSges(4,3) - t455;
t663 = t307 * t385 + t308 * t382;
t662 = t111 * t384 - t112 * t381;
t630 = -mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t578 = Ifges(3,4) * t382;
t449 = t385 * Ifges(3,2) + t578;
t657 = Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t449 / 0.2e1 + Ifges(4,5) * t710 + Ifges(4,6) * t711 + Ifges(4,3) * t650;
t654 = -t112 * mrSges(4,1) + t111 * mrSges(4,2);
t602 = t199 / 0.2e1;
t601 = t200 / 0.2e1;
t587 = t310 / 0.2e1;
t651 = -t466 / 0.2e1;
t196 = t317 * t377 + t318 * t546;
t108 = -t162 * t378 + t380 * t208;
t197 = t317 * t379 - t318 * t552;
t421 = -qJ(5) * t197 + t108;
t511 = qJD(5) * t377;
t641 = -t196 * t579 + (-qJD(6) * t379 - t511) * t378 - t421;
t483 = t377 * t514;
t157 = t377 * t161;
t529 = -pkin(4) * t554 + t157;
t640 = -qJD(6) * t380 + t483 + t162 * t546 - pkin(5) * t197 - (-qJ(6) * t318 - t208 * t379) * t378 - t529;
t638 = pkin(5) * t196 - t639;
t637 = -t378 * (m(7) * (pkin(5) + qJ(4)) + mrSges(7,1)) + mrSges(4,2);
t319 = t468 * t381;
t636 = t379 * (t311 * t378 - t319 * t380);
t635 = t379 * (t162 * t380 + t208 * t378);
t523 = t374 + t372;
t343 = -pkin(1) - t523;
t360 = pkin(8) * t533;
t255 = t381 * t343 + t360;
t631 = t378 ^ 2 + t380 ^ 2;
t214 = -t377 * t422 + t379 * t537;
t322 = t460 * qJD(2);
t527 = t381 * t322 + t343 * t515;
t103 = (-t514 + (-pkin(8) * qJD(3) - qJD(2) * t570) * t381) * t385 + (-pkin(8) * t519 - t481 + (qJD(2) * t378 - t380 * t515) * qJ(4)) * t382 + t527;
t356 = qJ(4) * t544;
t520 = qJD(2) * t382;
t505 = pkin(8) * t520;
t526 = t384 * t322 + t381 * t505;
t115 = -t382 * t480 + t414 * qJD(2) + (-t360 + (-t343 + t356) * t381) * qJD(3) + t526;
t143 = (qJD(3) * t426 - t378 * t512) * t382 + (pkin(8) + t427) * t518;
t28 = -t377 * t103 + t379 * (t115 * t380 + t143 * t378);
t190 = -qJ(4) * t422 + t255;
t316 = t384 * t343;
t201 = -t384 * t356 + t316 + (-pkin(8) * t381 - pkin(3)) * t385;
t498 = t378 * t537;
t339 = qJ(4) * t498;
t244 = -t339 + (pkin(8) + t583) * t382;
t87 = -t377 * t190 + t379 * (t201 * t380 + t244 * t378);
t612 = Ifges(4,1) * t602 + Ifges(4,4) * t601 + Ifges(4,5) * t587;
t585 = t318 / 0.2e1;
t582 = pkin(8) * t382;
t581 = -qJD(1) / 0.2e1;
t580 = qJD(3) / 0.2e1;
t577 = Ifges(3,4) * t385;
t576 = Ifges(4,4) * t318;
t575 = Ifges(4,4) * t381;
t574 = Ifges(4,4) * t384;
t573 = t317 * mrSges(4,3);
t572 = t318 * mrSges(4,3);
t563 = t291 * t377;
t562 = t291 * t378;
t292 = t383 * t384 - t385 * t531;
t561 = t292 * t378;
t293 = t383 * t381 + t384 * t532;
t560 = t293 * t378;
t192 = Ifges(4,2) * t317 - Ifges(4,6) * t466 + t576;
t541 = t381 * t192;
t538 = t382 * t383;
t536 = t382 * t386;
t119 = mrSges(6,1) * t159 + mrSges(6,3) * t229;
t120 = -mrSges(7,1) * t159 - mrSges(7,2) * t229;
t530 = -t119 + t120;
t320 = t468 * t384;
t284 = t377 * t320;
t528 = pkin(4) * t549 + t284;
t283 = pkin(3) * t552 + qJ(4) * t551;
t522 = t386 * pkin(1) + t383 * pkin(8);
t516 = qJD(3) * t382;
t506 = pkin(3) * t540;
t503 = qJ(4) * t562;
t502 = qJ(4) * t560;
t499 = t378 * t538;
t29 = t379 * t103 + t115 * t552 + t143 * t553;
t88 = t379 * t190 + t201 * t552 + t244 * t553;
t491 = Ifges(4,5) * t199 + Ifges(4,6) * t200 + Ifges(4,3) * t310;
t145 = t311 * t553 - t319 * t552 + t379 * t320;
t489 = -pkin(3) * t379 - pkin(4);
t485 = t381 * t516;
t484 = t382 * t515;
t475 = -t541 / 0.2e1;
t57 = -t99 * mrSges(7,1) + t142 * mrSges(7,2);
t32 = -t100 * mrSges(7,2) + t99 * mrSges(7,3);
t55 = t100 * mrSges(7,1) - t142 * mrSges(7,3);
t469 = -qJ(5) * t377 - pkin(3);
t467 = t509 / 0.2e1;
t58 = t100 * mrSges(6,1) + t142 * mrSges(6,2);
t63 = -t115 * t378 + t380 * t143;
t128 = -t201 * t378 + t380 * t244;
t204 = t380 * t311 + t319 * t378;
t463 = pkin(2) * t532 + pkin(9) * t536 + t522;
t462 = t380 * t484;
t177 = -t290 * t377 + t291 * t546;
t178 = -t290 * t379 - t291 * t552;
t274 = t290 * pkin(3);
t458 = t178 * pkin(4) + qJ(5) * t177 - t274;
t179 = t292 * t377 + t293 * t546;
t180 = t292 * t379 - t293 * t552;
t276 = t292 * pkin(3);
t457 = t180 * pkin(4) + qJ(5) * t179 + t276;
t249 = t377 * t540 - t382 * t496;
t250 = t281 * t382;
t456 = -t250 * pkin(4) - qJ(5) * t249 + t339;
t253 = -qJ(5) * t380 - t283;
t452 = mrSges(4,1) * t381 + mrSges(4,2) * t384;
t451 = Ifges(4,1) * t384 - t575;
t450 = Ifges(4,1) * t381 + t574;
t448 = -Ifges(4,2) * t381 + t574;
t447 = Ifges(4,2) * t384 + t575;
t446 = Ifges(3,5) * t385 - Ifges(3,6) * t382;
t445 = Ifges(4,5) * t384 - Ifges(4,6) * t381;
t444 = Ifges(4,5) * t381 + Ifges(4,6) * t384;
t441 = pkin(3) * t533 + qJ(4) * t500 + t356 + t523;
t78 = -qJ(5) * t286 - t88;
t424 = pkin(1) * t454;
t226 = t380 * t538 + t564;
t419 = -qJ(5) * t214 + t128;
t418 = -qJ(5) * t281 + t204;
t417 = t347 * t452;
t416 = t382 * (Ifges(3,1) * t385 - t578);
t132 = qJ(5) * t549 - t145;
t153 = -t292 * t546 + t293 * t377 - t386 * t493;
t213 = t377 * t537 + t379 * t422;
t413 = -g(1) * t153 - g(2) * (t563 + (t290 * t380 - t499) * t379) - g(3) * t213;
t227 = t380 * t536 - t561;
t412 = -g(1) * t227 - g(2) * t226 - g(3) * t286;
t152 = t290 * t552 - t291 * t379 - t377 * t499;
t409 = t384 * t518 - t485;
t408 = t381 * t518 + t484;
t222 = qJD(2) * t287 + t378 * t484;
t12 = -qJ(5) * t222 - qJD(5) * t286 - t29;
t138 = qJD(2) * t216 - t377 * t462 - t379 * t485;
t401 = -qJ(5) * t138 - qJD(5) * t214 + t63;
t395 = Ifges(4,5) * t382 + t385 * t451;
t394 = Ifges(4,6) * t382 + t385 * t448;
t393 = Ifges(4,3) * t382 + t385 * t445;
t391 = t293 * pkin(3) - qJ(4) * t561 + t386 * t356 + t463;
t354 = qJ(4) * t553;
t345 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t510;
t300 = t452 * t382;
t282 = pkin(3) * t546 - t354;
t257 = (-pkin(4) * t379 + t469) * t378;
t256 = t380 * t489 + t354;
t254 = -pkin(8) * t539 + t316;
t239 = -mrSges(4,1) * t466 - t572;
t238 = mrSges(4,2) * t466 + t573;
t237 = -pkin(8) * t487 + t294;
t224 = (-t379 * t579 + t469) * t378;
t223 = pkin(5) * t551 - t253;
t207 = pkin(5) * t553 + t354 + (-qJ(6) + t489) * t380;
t189 = t214 * t386;
t188 = t213 * t386;
t187 = t214 * t383;
t186 = t213 * t383;
t168 = -qJD(3) * t255 + t526;
t167 = (-t382 * t519 - t385 * t517) * pkin(8) + t527;
t156 = -mrSges(4,2) * t310 + mrSges(4,3) * t200;
t155 = mrSges(4,1) * t310 - mrSges(4,3) * t199;
t154 = t293 * t379 + (t292 * t380 + t378 * t536) * t377;
t151 = -t290 * t546 + t383 * t493 - t563;
t144 = -t284 + t636;
t137 = qJD(2) * t215 - t377 * t485 + t379 * t462;
t133 = t528 - t636;
t127 = pkin(4) * t280 + t418;
t126 = t230 * t379 - t555 * t631;
t125 = t547 * t631 + t566;
t124 = -mrSges(4,1) * t200 + mrSges(4,2) * t199;
t122 = -pkin(5) * t280 - t132;
t121 = mrSges(6,1) * t160 - mrSges(6,2) * t229;
t118 = mrSges(7,1) * t160 + mrSges(7,3) * t229;
t117 = -mrSges(5,1) * t229 - mrSges(5,3) * t160;
t116 = mrSges(5,2) * t229 - mrSges(5,3) * t159;
t110 = t280 * t579 + t418;
t109 = t319 * t546 + pkin(5) * t281 + (qJ(6) * t384 - t311 * t379) * t378 + t528;
t106 = t199 * Ifges(4,4) + t200 * Ifges(4,2) + t310 * Ifges(4,6);
t97 = t100 * mrSges(6,3);
t95 = t100 * mrSges(5,2);
t93 = -mrSges(6,2) * t159 - mrSges(6,3) * t160;
t92 = mrSges(5,1) * t159 + mrSges(5,2) * t160;
t91 = -mrSges(7,2) * t160 + mrSges(7,3) * t159;
t80 = pkin(4) * t213 + t419;
t79 = -pkin(4) * t286 - t87;
t66 = -t157 + t635;
t59 = t213 * t579 + t419;
t56 = mrSges(6,1) * t99 - mrSges(6,3) * t142;
t54 = mrSges(5,1) * t142 - mrSges(5,3) * t100;
t53 = -mrSges(5,2) * t142 - mrSges(5,3) * t99;
t52 = -pkin(5) * t213 - t78;
t51 = t529 - t635;
t49 = pkin(4) * t196 + t421;
t42 = pkin(5) * t214 - t286 * t579 - t87;
t34 = -t99 * mrSges(6,2) - t97;
t33 = t99 * mrSges(5,1) + t95;
t15 = -pkin(4) * t222 - t28;
t13 = pkin(4) * t137 + t401;
t11 = -pkin(5) * t137 - t12;
t10 = pkin(5) * t138 - qJD(6) * t286 - t222 * t579 - t28;
t9 = qJD(6) * t213 + t137 * t579 + t401;
t17 = [(-qJDD(2) * mrSges(3,1) + t124) * t582 + t455 * t569 + (qJD(2) * t393 - t444 * t516) * t651 + t347 * (mrSges(4,1) * t408 + mrSges(4,2) * t409) - (t192 * t384 + t193 * t381) * t516 / 0.2e1 - t106 * t540 / 0.2e1 + (t211 * mrSges(4,1) - t212 * mrSges(4,2) - t657) * t520 + t669 * t688 + t577 * t689 + t449 * t690 + (-t668 * m(5) + t291 * mrSges(4,1) - t290 * mrSges(4,2) + t691 * (t152 * pkin(4) + qJ(5) * t151 + t668) + t687 * t386 + (-m(3) - m(4)) * t375 + (m(3) * pkin(1) - m(4) * t432 + mrSges(2,1) - t664) * t383 + t655 * t226 - t699 * t152 + t630 * t151) * g(1) + (-m(3) * t522 - m(4) * t463 - m(5) * t391 - t293 * mrSges(4,1) - t292 * mrSges(4,2) - mrSges(4,3) * t536 + t691 * (t154 * pkin(4) + qJ(5) * t153 + t391) + (-t455 - mrSges(2,1)) * t386 + t687 * t383 - t655 * t227 - t699 * t154 + t630 * t153) * g(2) + t416 * t467 + (t670 / 0.2e1 + t475) * t518 + qJD(2) ^ 2 * t446 / 0.2e1 + t317 * (qJD(2) * t394 - t447 * t516) / 0.2e1 - t424 * t509 + (-t695 + t719) * t222 + m(7) * (t1 * t42 + t10 * t14 + t11 * t27 + t2 * t59 + t26 * t9 + t3 * t52) + m(6) * (t12 * t37 + t13 * t31 + t15 * t35 + t4 * t78 + t5 * t80 + t6 * t79) + m(5) * (t128 * t16 + t28 * t45 + t29 * t46 + t63 * t90 + t7 * t87 + t8 * t88) - t345 * t505 + t725 * t213 + t726 * t214 - t712 * t286 + (t331 * t582 + t663) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(8) * t663) + (-t693 + t715) * t137 + (-t694 + t717) * t138 + (-t111 * t540 - t112 * t537 - t211 * t409 - t212 * t408) * mrSges(4,3) + t537 * t612 + Ifges(2,3) * qJDD(1) - pkin(1) * (-mrSges(3,1) * t330 + mrSges(3,2) * t331) + t279 * t300 + t254 * t155 + t255 * t156 + t167 * t238 + t168 * t239 + t128 * t33 + t29 * t116 + t28 * t117 + t10 * t118 + t12 * t119 + t11 * t120 + t15 * t121 + t87 * t54 + t88 * t53 + t9 * t91 + t63 * t92 + t13 * t93 + t78 * t56 + t79 * t58 + t80 * t34 + t52 * t57 + t59 * t32 + t42 * t55 + m(4) * (t111 * t255 + t112 * t254 + t167 * t212 + t168 * t211 + t347 * t688) + ((-Ifges(3,2) * t382 + t577) * t467 - t491 / 0.2e1 + pkin(8) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t330) + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t689 + Ifges(3,2) * t690 - Ifges(4,3) * t587 - Ifges(4,6) * t601 - Ifges(4,5) * t602 + t654) * t385 + (m(4) * t279 * pkin(8) + Ifges(3,1) * t331 + Ifges(3,4) * t690 + Ifges(3,5) * qJDD(2) + t445 * t587 + t448 * t601 + t451 * t602) * t382 + (qJD(2) * t395 - t450 * t516) * t585; t712 * t549 + (t393 * t650 - t211 * (mrSges(4,1) * t382 - mrSges(4,3) * t533) - t212 * (-mrSges(4,2) * t382 - mrSges(4,3) * t539) + (t424 - t416 / 0.2e1) * qJD(1)) * qJD(1) - t669 * t371 - (-Ifges(3,2) * t521 + t369 + t670) * t510 / 0.2e1 + t657 * t521 + (t445 * t651 + t417 + t475) * qJD(3) + t675 * t116 + (t144 * t7 + t145 * t8 + t16 * t204 + t45 * t676 + t46 * t675 + t674 * t90) * m(5) + t676 * t117 + t677 * t121 + (-t667 * m(5) - m(4) * t363 + t691 * (-t189 * pkin(4) - t188 * qJ(5) + t667) + t699 * t189 - t630 * t188 + t698 * t386) * g(1) + (-t666 * m(5) - m(4) * t361 + t691 * (-t187 * pkin(4) - qJ(5) * t186 + t666) + t699 * t187 - t630 * t186 + t698 * t383) * g(2) + (-m(4) * t523 - m(5) * t441 - t385 * t453 + t691 * (t216 * pkin(4) + qJ(5) * t215 + t441) - t655 * t287 - t699 * t216 + t630 * t215 + t664) * g(3) - t279 * t453 + t193 * t515 / 0.2e1 - t446 * t509 / 0.2e1 + (t658 / 0.2e1 + t719) * t486 + (t709 + t720) * t269 + (t395 * t581 + t451 * t580) * t318 + t725 * t280 + t726 * t281 + (m(4) * ((-t211 * t384 - t212 * t381) * qJD(3) + t662) - t239 * t515 - t238 * t517 + t384 * t156 - t381 * t155) * pkin(9) + (-t211 * t515 - t212 * t517 + t662) * mrSges(4,3) + (t394 * t581 + t448 * t580) * t317 + (-pkin(2) * t279 - t211 * t236 - t212 * t237 - t347 * t371) * m(4) + (t541 / 0.2e1 - t417) * t510 + (t683 / 0.2e1 + t715) * t266 + (t707 + t716) * t202 + (t682 / 0.2e1 + t717) * t267 + (t708 + t718) * t203 + t381 * t612 + (-mrSges(5,1) * t671 + mrSges(5,3) * t672) * t45 + (-mrSges(7,1) * t672 + mrSges(7,3) * t671) * t14 + (-mrSges(6,1) * t672 - mrSges(6,2) * t671) * t35 + (mrSges(5,2) * t671 + mrSges(5,3) * t673) * t46 + (mrSges(7,1) * t673 - mrSges(7,2) * t671) * t27 + (mrSges(6,2) * t673 + mrSges(6,3) * t672) * t31 + (-mrSges(6,1) * t673 + mrSges(6,3) * t671) * t37 + (mrSges(7,2) * t672 - mrSges(7,3) * t673) * t26 + (-mrSges(5,1) * t673 - mrSges(5,2) * t672) * t90 + t674 * t92 + t345 * t370 + Ifges(3,6) * t330 + Ifges(3,5) * t331 - t308 * mrSges(3,1) - t307 * mrSges(3,2) - t237 * t238 - t236 * t239 + Ifges(3,3) * qJDD(2) + t204 * t33 + t144 * t54 + t145 * t53 + t132 * t56 + t133 * t58 - pkin(2) * t124 + t127 * t34 + t122 * t57 + t109 * t55 + t110 * t32 + t384 * t106 / 0.2e1 + t678 * t119 + t679 * t93 + (t127 * t5 + t132 * t4 + t133 * t6 + t31 * t679 + t35 * t677 + t37 * t678) * m(6) + t680 * t120 + t681 * t118 + t684 * t91 + (t1 * t109 + t110 * t2 + t122 * t3 + t14 * t681 + t26 * t684 + t27 * t680) * m(7) + t444 * t587 + t447 * t601 + t450 * t602; t731 * t551 + (Ifges(4,5) * t317 - Ifges(4,6) * t318) * t650 + (-t483 - t66) * t117 + (-t238 + t573) * t211 + (m(6) * (qJD(4) * t35 - qJD(5) * t31) + t696) * t553 + (-t108 * t90 + t282 * t7 + t283 * t8 - t45 * t66 - t46 * t67) * m(5) + (-m(7) * t458 - m(5) * (-t274 + t503) - m(6) * (t458 + t503) + mrSges(4,1) * t290 + t648 * t562 + t637 * t291 - t699 * t178 + t630 * t177) * g(2) + (-m(7) * t457 - m(5) * (t276 + t502) - m(6) * (t457 + t502) - mrSges(4,1) * t292 + t648 * t560 + t637 * t293 - t699 * t180 + t630 * t179) * g(1) + (-m(5) * (t339 - t506) - m(6) * (t456 - t506) - m(7) * t456 - (m(7) * (pkin(5) * t549 - t583) + mrSges(7,1) * t549) * t382 + t300 + t648 * t498 + t699 * t250 - t630 * t249) * g(3) + (t239 + t572) * t212 - t654 + (t695 + t720) * t554 + (Ifges(4,1) * t317 - t576) * t710 + (-Ifges(4,2) * t318 + t193 + t304) * t711 - t712 * t380 + (t5 * (mrSges(6,2) * t379 - mrSges(6,3) * t377) + t2 * (-mrSges(7,2) * t377 - mrSges(7,3) * t379) + t16 * (-mrSges(5,1) * t379 + mrSges(5,2) * t377) + (Ifges(5,4) * t377 + Ifges(5,2) * t379) * t616 + m(5) * (-pkin(3) * t16 + (-t377 * t45 + t379 * t46) * qJD(4)) - t93 * t511 + (-Ifges(6,2) * t377 - Ifges(6,6) * t379) * t614 - pkin(3) * t33 + (t377 * t642 - t379 * t643) * t615 + (t377 * t647 + t379 * t646) * t613 + (t377 * t729 - t379 * t728) * t609) * t378 + (-t51 + t483) * t121 + (t693 + t716) * t196 + (t694 + t718) * t197 - t347 * (mrSges(4,1) * t318 + mrSges(4,2) * t317) + t282 * t54 + t283 * t53 + t256 * t58 + t257 * t34 + t253 * t56 + t223 * t57 + t224 * t32 + t207 * t55 - t108 * t92 - t49 * t93 + t192 * t585 + t638 * t120 + t639 * t119 + (t253 * t4 + t256 * t6 + t257 * t5 - t31 * t49 - t35 * t51 + t37 * t639) * m(6) + t640 * t118 + t641 * t91 + (t1 * t207 + t14 * t640 + t2 * t224 + t223 * t3 + t26 * t641 + t27 * t638) * m(7) + t702 * t116 + t491; t95 - t97 + t649 * t99 + (-t116 - t530) * t126 + (t117 - t118 - t121) * t125 + t32 + (-t125 * t14 - t126 * t27 + t2 + t412) * m(7) + (-t125 * t35 + t126 * t37 + t412 + t5) * m(6) + (t125 * t45 - t126 * t46 + t16 + t412) * m(5); t530 * t229 + (t91 + t93) * t160 + t58 + t55 + (t160 * t26 + t229 * t27 + t1 + t413) * m(7) + (t160 * t31 - t229 * t37 + t413 + t6) * m(6); -t229 * t118 - t159 * t91 + (-g(1) * t154 + g(2) * t152 - g(3) * t214 - t14 * t229 - t159 * t26 + t3) * m(7) + t57;];
tau  = t17;
