% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:09
% EndTime: 2019-03-09 18:26:39
% DurationCPUTime: 55.86s
% Computational Cost: add. (30067->1044), mult. (67451->1400), div. (0->0), fcn. (50452->18), ass. (0->454)
t418 = sin(qJ(2));
t423 = cos(qJ(2));
t464 = pkin(2) * t418 - pkin(8) * t423;
t357 = t464 * qJD(1);
t422 = cos(qJ(3));
t417 = sin(qJ(3));
t506 = qJD(1) * t418;
t484 = t417 * t506;
t277 = pkin(7) * t484 + t422 * t357;
t515 = t422 * t423;
t444 = pkin(3) * t418 - qJ(4) * t515;
t414 = -qJ(4) - pkin(8);
t469 = qJD(3) * t414;
t686 = -qJD(1) * t444 - qJD(4) * t417 + t422 * t469 - t277;
t332 = t417 * t357;
t498 = qJD(4) * t422;
t519 = t418 * t422;
t521 = t417 * t423;
t685 = t332 + (-pkin(7) * t519 - qJ(4) * t521) * qJD(1) - t417 * t469 - t498;
t412 = sin(pkin(11));
t413 = cos(pkin(11));
t348 = t412 * t422 + t413 * t417;
t440 = t423 * t348;
t291 = qJD(1) * t440;
t324 = t348 * qJD(3);
t677 = -t291 + t324;
t446 = t412 * t417 - t413 * t422;
t439 = t446 * t423;
t292 = qJD(1) * t439;
t325 = t446 * qJD(3);
t684 = t292 - t325;
t641 = t685 * t412 + t413 * t686;
t640 = t412 * t686 - t685 * t413;
t683 = -pkin(4) * t506 - pkin(9) * t684 + t641;
t682 = pkin(9) * t677 - t640;
t503 = qJD(2) * t422;
t354 = -t484 + t503;
t482 = t422 * t506;
t355 = qJD(2) * t417 + t482;
t256 = t354 * t412 + t355 * t413;
t416 = sin(qJ(5));
t421 = cos(qJ(5));
t466 = t413 * t354 - t355 * t412;
t189 = t256 * t421 + t416 * t466;
t415 = sin(qJ(6));
t420 = cos(qJ(6));
t656 = -t256 * t416 + t421 * t466;
t111 = t189 * t420 + t415 * t656;
t397 = pkin(7) * t506;
t372 = -qJD(2) * pkin(2) + t397;
t272 = -pkin(3) * t354 + qJD(4) + t372;
t198 = -pkin(4) * t466 + t272;
t124 = -pkin(5) * t656 + t198;
t675 = -t189 * t415 + t420 * t656;
t494 = qJD(1) * qJD(2);
t360 = qJDD(1) * t418 + t423 * t494;
t247 = qJD(3) * t354 + qJDD(2) * t417 + t360 * t422;
t248 = -qJD(3) * t355 + qJDD(2) * t422 - t360 * t417;
t172 = -t247 * t412 + t248 * t413;
t174 = t247 * t413 + t248 * t412;
t74 = qJD(5) * t656 + t172 * t416 + t174 * t421;
t75 = -qJD(5) * t189 + t172 * t421 - t174 * t416;
t29 = qJD(6) * t675 + t415 * t75 + t420 * t74;
t30 = -qJD(6) * t111 - t415 * t74 + t420 * t75;
t359 = qJDD(1) * t423 - t418 * t494;
t346 = qJDD(3) - t359;
t338 = qJDD(5) + t346;
t322 = qJDD(6) + t338;
t490 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t322;
t505 = qJD(1) * t423;
t384 = qJD(3) - t505;
t374 = qJD(5) + t384;
t465 = pkin(2) * t423 + pkin(8) * t418;
t364 = -pkin(1) - t465;
t339 = t364 * qJD(1);
t398 = pkin(7) * t505;
t373 = qJD(2) * pkin(8) + t398;
t261 = t422 * t339 - t373 * t417;
t218 = -qJ(4) * t355 + t261;
t205 = pkin(3) * t384 + t218;
t262 = t339 * t417 + t373 * t422;
t219 = qJ(4) * t354 + t262;
t211 = t412 * t219;
t133 = t413 * t205 - t211;
t658 = pkin(9) * t256;
t115 = pkin(4) * t384 + t133 - t658;
t523 = t413 * t219;
t134 = t412 * t205 + t523;
t654 = pkin(9) * t466;
t116 = t134 + t654;
t57 = t421 * t115 - t116 * t416;
t678 = pkin(10) * t189;
t48 = t57 - t678;
t43 = pkin(5) * t374 + t48;
t58 = t115 * t416 + t116 * t421;
t674 = pkin(10) * t656;
t49 = t58 + t674;
t535 = t415 * t49;
t16 = t420 * t43 - t535;
t529 = qJDD(1) * pkin(1);
t263 = -pkin(2) * t359 - pkin(8) * t360 - t529;
t344 = t359 * pkin(7);
t319 = qJDD(2) * pkin(8) + t344;
t161 = -qJD(3) * t262 + t422 * t263 - t319 * t417;
t113 = pkin(3) * t346 - qJ(4) * t247 - qJD(4) * t355 + t161;
t499 = qJD(3) * t422;
t501 = qJD(3) * t417;
t160 = t417 * t263 + t422 * t319 + t339 * t499 - t373 * t501;
t118 = qJ(4) * t248 + qJD(4) * t354 + t160;
t59 = t413 * t113 - t118 * t412;
t42 = pkin(4) * t346 - pkin(9) * t174 + t59;
t60 = t412 * t113 + t413 * t118;
t45 = pkin(9) * t172 + t60;
t13 = -qJD(5) * t58 - t416 * t45 + t421 * t42;
t6 = pkin(5) * t338 - pkin(10) * t74 + t13;
t496 = qJD(5) * t421;
t497 = qJD(5) * t416;
t12 = t115 * t496 - t116 * t497 + t416 * t42 + t421 * t45;
t7 = pkin(10) * t75 + t12;
t2 = qJD(6) * t16 + t415 * t6 + t420 * t7;
t530 = t420 * t49;
t17 = t415 * t43 + t530;
t3 = -qJD(6) * t17 - t415 * t7 + t420 * t6;
t663 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t445 = t490 - t663;
t489 = Ifges(6,5) * t74 + Ifges(6,6) * t75 + Ifges(6,3) * t338;
t561 = mrSges(7,3) * t17;
t562 = mrSges(7,3) * t16;
t365 = qJD(6) + t374;
t571 = -t365 / 0.2e1;
t595 = -t111 / 0.2e1;
t597 = -t675 / 0.2e1;
t102 = Ifges(7,4) * t675;
t55 = Ifges(7,1) * t111 + Ifges(7,5) * t365 + t102;
t607 = -t55 / 0.2e1;
t540 = Ifges(7,4) * t111;
t54 = Ifges(7,2) * t675 + Ifges(7,6) * t365 + t540;
t609 = -t54 / 0.2e1;
t662 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t681 = -(mrSges(7,1) * t124 + Ifges(7,4) * t595 + Ifges(7,2) * t597 + Ifges(7,6) * t571 - t561 + t609) * t111 + (-mrSges(7,2) * t124 + Ifges(7,1) * t595 + Ifges(7,4) * t597 + Ifges(7,5) * t571 + t562 + t607) * t675 + t445 + t489 - t662;
t250 = t348 * t421 - t416 * t446;
t181 = -qJD(5) * t250 - t324 * t421 + t325 * t416;
t209 = -t291 * t421 + t292 * t416;
t680 = t181 - t209;
t556 = pkin(5) * t189;
t368 = t414 * t417;
t370 = t414 * t422;
t266 = t413 * t368 + t370 * t412;
t235 = -pkin(9) * t348 + t266;
t267 = t412 * t368 - t413 * t370;
t236 = -pkin(9) * t446 + t267;
t155 = t416 * t235 + t421 * t236;
t651 = -qJD(5) * t155 + t682 * t416 + t421 * t683;
t650 = t235 * t496 - t236 * t497 + t416 * t683 - t682 * t421;
t141 = -t218 * t412 - t523;
t126 = t141 - t654;
t142 = t413 * t218 - t211;
t127 = t142 - t658;
t558 = pkin(3) * t413;
t390 = pkin(4) + t558;
t559 = pkin(3) * t412;
t316 = t390 * t416 + t421 * t559;
t646 = -t316 * qJD(5) - t421 * t126 + t127 * t416;
t315 = t421 * t390 - t416 * t559;
t645 = t315 * qJD(5) - t416 * t126 - t421 * t127;
t676 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t673 = Ifges(4,3) + Ifges(5,3);
t672 = pkin(10) * t680 + t650;
t249 = -t348 * t416 - t421 * t446;
t180 = qJD(5) * t249 - t324 * t416 - t325 * t421;
t210 = -t291 * t416 - t292 * t421;
t671 = -pkin(5) * t506 + t651 + (-t180 + t210) * pkin(10);
t670 = t645 + t678;
t669 = t674 + t646;
t483 = t417 * t505;
t638 = -t398 + (-t483 + t501) * pkin(3);
t462 = mrSges(3,1) * t418 + mrSges(3,2) * t423;
t545 = Ifges(3,4) * t418;
t666 = pkin(1) * t462 - t418 * (Ifges(3,1) * t423 - t545) / 0.2e1;
t411 = qJ(3) + pkin(11);
t402 = cos(t411);
t407 = t422 * pkin(3);
t363 = pkin(4) * t402 + t407;
t403 = qJ(5) + t411;
t389 = cos(t403);
t309 = pkin(5) * t389 + t363;
t301 = pkin(2) + t309;
t356 = pkin(2) + t363;
t392 = qJ(6) + t403;
t381 = sin(t392);
t382 = cos(t392);
t388 = sin(t403);
t393 = t407 + pkin(2);
t401 = sin(t411);
t461 = -mrSges(4,1) * t422 + mrSges(4,2) * t417;
t665 = m(4) * pkin(2) + m(5) * t393 + m(6) * t356 + m(7) * t301 + mrSges(5,1) * t402 + mrSges(6,1) * t389 + mrSges(7,1) * t382 - mrSges(5,2) * t401 - mrSges(6,2) * t388 - mrSges(7,2) * t381 - t461;
t410 = -pkin(9) + t414;
t404 = -pkin(10) + t410;
t664 = -m(4) * pkin(8) + m(5) * t414 + m(6) * t410 + m(7) * t404 - t676;
t661 = -t261 * mrSges(4,1) + t262 * mrSges(4,2);
t660 = -t161 * mrSges(4,1) - t59 * mrSges(5,1) + t160 * mrSges(4,2) + t60 * mrSges(5,2);
t455 = t423 * Ifges(3,2) + t545;
t659 = t134 * mrSges(5,2) + t17 * mrSges(7,2) + t58 * mrSges(6,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t455 / 0.2e1 - t133 * mrSges(5,1) - t16 * mrSges(7,1) - t57 * mrSges(6,1);
t639 = pkin(4) * t677 + t638;
t502 = qJD(2) * t423;
t433 = t417 * t502 + t418 * t499;
t419 = sin(qJ(1));
t424 = cos(qJ(1));
t657 = g(1) * t424 + g(2) * t419;
t615 = m(5) * pkin(3);
t613 = t29 / 0.2e1;
t612 = t30 / 0.2e1;
t605 = t74 / 0.2e1;
t604 = t75 / 0.2e1;
t588 = t172 / 0.2e1;
t587 = t174 / 0.2e1;
t582 = t247 / 0.2e1;
t581 = t248 / 0.2e1;
t576 = t322 / 0.2e1;
t575 = t338 / 0.2e1;
t574 = t346 / 0.2e1;
t154 = t421 * t235 - t236 * t416;
t128 = -pkin(10) * t250 + t154;
t129 = pkin(10) * t249 + t155;
t70 = t128 * t415 + t129 * t420;
t653 = -qJD(6) * t70 - t415 * t672 + t420 * t671;
t69 = t128 * t420 - t129 * t415;
t652 = qJD(6) * t69 + t415 * t671 + t420 * t672;
t312 = pkin(5) + t315;
t228 = t312 * t415 + t316 * t420;
t649 = -qJD(6) * t228 - t415 * t670 + t420 * t669;
t227 = t312 * t420 - t316 * t415;
t648 = qJD(6) * t227 + t415 * t669 + t420 * t670;
t647 = t615 + mrSges(4,1);
t350 = t422 * t364;
t258 = -qJ(4) * t519 + t350 + (-pkin(7) * t417 - pkin(3)) * t423;
t386 = pkin(7) * t515;
t294 = t417 * t364 + t386;
t522 = t417 * t418;
t265 = -qJ(4) * t522 + t294;
t192 = t413 * t258 - t265 * t412;
t311 = t446 * t418;
t152 = -pkin(4) * t423 + pkin(9) * t311 + t192;
t193 = t412 * t258 + t413 * t265;
t310 = t348 * t418;
t159 = -pkin(9) * t310 + t193;
t91 = t416 * t152 + t421 * t159;
t642 = -pkin(5) * t680 + t639;
t637 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t354 + mrSges(4,2) * t355 + mrSges(3,3) * t506;
t458 = -mrSges(7,1) * t381 - mrSges(7,2) * t382;
t636 = mrSges(6,1) * t388 + mrSges(6,2) * t389 - t458;
t635 = Ifges(4,5) * t247 + Ifges(5,5) * t174 + Ifges(4,6) * t248 + Ifges(5,6) * t172 + t346 * t673;
t514 = t423 * t424;
t289 = -t388 * t514 + t419 * t389;
t290 = t419 * t388 + t389 * t514;
t275 = -t381 * t514 + t419 * t382;
t276 = t419 * t381 + t382 * t514;
t512 = t275 * mrSges(7,1) - t276 * mrSges(7,2);
t634 = -t289 * mrSges(6,1) + t290 * mrSges(6,2) - t512;
t517 = t419 * t423;
t287 = t388 * t517 + t389 * t424;
t288 = t388 * t424 - t389 * t517;
t273 = t381 * t517 + t382 * t424;
t274 = t381 * t424 - t382 * t517;
t513 = -t273 * mrSges(7,1) + t274 * mrSges(7,2);
t633 = t287 * mrSges(6,1) - t288 * mrSges(6,2) - t513;
t343 = Ifges(4,4) * t354;
t240 = t355 * Ifges(4,1) + t384 * Ifges(4,5) + t343;
t396 = Ifges(3,4) * t505;
t632 = Ifges(3,1) * t506 + Ifges(3,5) * qJD(2) + t422 * t240 + t396;
t345 = t360 * pkin(7);
t631 = t344 * t423 + t345 * t418;
t630 = t160 * t422 - t161 * t417;
t629 = m(5) + m(6) + m(7);
t628 = -mrSges(6,1) * t198 + mrSges(6,3) * t58;
t627 = mrSges(6,2) * t198 - mrSges(6,3) * t57;
t626 = -m(4) - m(3) - t629;
t625 = t355 * Ifges(4,5) + t256 * Ifges(5,5) + t189 * Ifges(6,5) + t111 * Ifges(7,5) + t354 * Ifges(4,6) + Ifges(5,6) * t466 + t656 * Ifges(6,6) + Ifges(7,6) * t675 + t374 * Ifges(6,3) + t365 * Ifges(7,3) + t384 * t673;
t557 = pkin(3) * t417;
t362 = pkin(4) * t401 + t557;
t555 = pkin(5) * t388;
t308 = t362 + t555;
t624 = -m(6) * t362 - m(7) * t308;
t463 = mrSges(3,1) * t423 - mrSges(3,2) * t418;
t622 = t418 * t676 + mrSges(2,1) + t463;
t621 = mrSges(2,2) - mrSges(3,3) + t624;
t617 = Ifges(7,4) * t613 + Ifges(7,2) * t612 + Ifges(7,6) * t576;
t616 = Ifges(7,1) * t613 + Ifges(7,4) * t612 + Ifges(7,5) * t576;
t614 = m(7) * pkin(5);
t611 = Ifges(6,4) * t605 + Ifges(6,2) * t604 + Ifges(6,6) * t575;
t610 = Ifges(6,1) * t605 + Ifges(6,4) * t604 + Ifges(6,5) * t575;
t608 = t54 / 0.2e1;
t606 = t55 / 0.2e1;
t603 = Ifges(5,4) * t587 + Ifges(5,2) * t588 + Ifges(5,6) * t574;
t602 = Ifges(5,1) * t587 + Ifges(5,4) * t588 + Ifges(5,5) * t574;
t541 = Ifges(6,4) * t189;
t97 = Ifges(6,2) * t656 + t374 * Ifges(6,6) + t541;
t601 = -t97 / 0.2e1;
t600 = t97 / 0.2e1;
t179 = Ifges(6,4) * t656;
t98 = t189 * Ifges(6,1) + t374 * Ifges(6,5) + t179;
t599 = -t98 / 0.2e1;
t598 = t98 / 0.2e1;
t596 = t675 / 0.2e1;
t594 = t111 / 0.2e1;
t593 = Ifges(4,1) * t582 + Ifges(4,4) * t581 + Ifges(4,5) * t574;
t170 = t256 * Ifges(5,4) + Ifges(5,2) * t466 + t384 * Ifges(5,6);
t592 = -t170 / 0.2e1;
t591 = t170 / 0.2e1;
t171 = t256 * Ifges(5,1) + Ifges(5,4) * t466 + t384 * Ifges(5,5);
t590 = -t171 / 0.2e1;
t589 = t171 / 0.2e1;
t586 = -t656 / 0.2e1;
t585 = t656 / 0.2e1;
t584 = -t189 / 0.2e1;
t583 = t189 / 0.2e1;
t580 = -t466 / 0.2e1;
t579 = t466 / 0.2e1;
t578 = -t256 / 0.2e1;
t577 = t256 / 0.2e1;
t572 = t355 / 0.2e1;
t570 = t365 / 0.2e1;
t569 = -t374 / 0.2e1;
t568 = t374 / 0.2e1;
t567 = -t384 / 0.2e1;
t566 = t384 / 0.2e1;
t560 = pkin(3) * t355;
t552 = g(3) * t418;
t405 = t418 * pkin(7);
t549 = mrSges(5,3) * t133;
t548 = mrSges(5,3) * t134;
t547 = mrSges(6,3) * t656;
t546 = mrSges(6,3) * t189;
t544 = Ifges(3,4) * t423;
t543 = Ifges(4,4) * t417;
t542 = Ifges(4,4) * t422;
t538 = t261 * mrSges(4,3);
t537 = t262 * mrSges(4,3);
t536 = t355 * Ifges(4,4);
t520 = t417 * t424;
t518 = t419 * t417;
t358 = t464 * qJD(2);
t504 = qJD(2) * t418;
t488 = pkin(7) * t504;
t508 = t422 * t358 + t417 * t488;
t178 = -t418 * t498 + t444 * qJD(2) + (-t386 + (qJ(4) * t418 - t364) * t417) * qJD(3) + t508;
t509 = t417 * t358 + t364 * t499;
t190 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t519 + (-qJD(4) * t418 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t423) * t417 + t509;
t105 = t412 * t178 + t413 * t190;
t361 = pkin(3) * t522 + t405;
t500 = qJD(3) * t418;
t400 = pkin(7) * t502;
t286 = pkin(3) * t433 + t400;
t239 = t354 * Ifges(4,2) + t384 * Ifges(4,6) + t536;
t479 = -t417 * t239 / 0.2e1;
t37 = -t75 * mrSges(6,1) + t74 * mrSges(6,2);
t10 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t99 = -t172 * mrSges(5,1) + t174 * mrSges(5,2);
t90 = t421 * t152 - t159 * t416;
t104 = t413 * t178 - t190 * t412;
t295 = pkin(4) * t446 - t393;
t259 = pkin(4) * t310 + t361;
t208 = pkin(4) * t256 + t560;
t320 = -qJDD(2) * pkin(2) + t345;
t460 = mrSges(4,1) * t417 + mrSges(4,2) * t422;
t457 = Ifges(4,1) * t422 - t543;
t456 = Ifges(4,1) * t417 + t542;
t454 = -Ifges(4,2) * t417 + t542;
t453 = Ifges(4,2) * t422 + t543;
t452 = Ifges(3,5) * t423 - Ifges(3,6) * t418;
t451 = Ifges(4,5) * t422 - Ifges(4,6) * t417;
t450 = Ifges(4,5) * t417 + Ifges(4,6) * t422;
t226 = -t310 * t416 - t311 * t421;
t76 = -pkin(5) * t423 - pkin(10) * t226 + t90;
t225 = -t310 * t421 + t311 * t416;
t79 = pkin(10) * t225 + t91;
t38 = -t415 * t79 + t420 * t76;
t39 = t415 * t76 + t420 * t79;
t146 = t225 * t420 - t226 * t415;
t147 = t225 * t415 + t226 * t420;
t182 = t249 * t420 - t250 * t415;
t183 = t249 * t415 + t250 * t420;
t449 = t301 * t423 - t404 * t418;
t448 = t356 * t423 - t410 * t418;
t447 = t393 * t423 - t414 * t418;
t229 = -qJD(2) * t440 + t446 * t500;
t194 = -pkin(4) * t229 + t286;
t330 = -t417 * t514 + t419 * t422;
t328 = t417 * t517 + t422 * t424;
t442 = t372 * t460;
t230 = -qJD(2) * t439 - t324 * t418;
t86 = pkin(4) * t504 - pkin(9) * t230 + t104;
t89 = pkin(9) * t229 + t105;
t31 = t152 * t496 - t159 * t497 + t416 * t86 + t421 * t89;
t434 = -t417 * t500 + t422 * t502;
t204 = -pkin(3) * t248 + qJDD(4) + t320;
t432 = Ifges(4,5) * t418 + t423 * t457;
t431 = Ifges(4,6) * t418 + t423 * t454;
t430 = Ifges(4,3) * t418 + t423 * t451;
t125 = -pkin(4) * t172 + t204;
t32 = -qJD(5) * t91 - t416 * t89 + t421 * t86;
t367 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t505;
t340 = t460 * t418;
t331 = t422 * t514 + t518;
t329 = -t419 * t515 + t520;
t306 = t419 * t401 + t402 * t514;
t305 = -t401 * t514 + t419 * t402;
t304 = t401 * t424 - t402 * t517;
t303 = t401 * t517 + t402 * t424;
t293 = -pkin(7) * t521 + t350;
t285 = mrSges(4,1) * t384 - mrSges(4,3) * t355;
t284 = -mrSges(4,2) * t384 + mrSges(4,3) * t354;
t278 = -pkin(7) * t482 + t332;
t224 = mrSges(5,1) * t384 - mrSges(5,3) * t256;
t223 = -mrSges(5,2) * t384 + mrSges(5,3) * t466;
t216 = -qJD(3) * t294 + t508;
t215 = (-t418 * t503 - t423 * t501) * pkin(7) + t509;
t207 = -mrSges(4,2) * t346 + mrSges(4,3) * t248;
t206 = mrSges(4,1) * t346 - mrSges(4,3) * t247;
t203 = -pkin(5) * t249 + t295;
t191 = -mrSges(5,1) * t466 + mrSges(5,2) * t256;
t184 = -mrSges(4,1) * t248 + mrSges(4,2) * t247;
t176 = -pkin(5) * t225 + t259;
t163 = mrSges(6,1) * t374 - t546;
t162 = -mrSges(6,2) * t374 + t547;
t156 = t247 * Ifges(4,4) + t248 * Ifges(4,2) + t346 * Ifges(4,6);
t145 = mrSges(5,1) * t346 - mrSges(5,3) * t174;
t144 = -mrSges(5,2) * t346 + mrSges(5,3) * t172;
t137 = t209 * t415 + t210 * t420;
t136 = t209 * t420 - t210 * t415;
t130 = t208 + t556;
t123 = -qJD(5) * t226 + t229 * t421 - t230 * t416;
t122 = qJD(5) * t225 + t229 * t416 + t230 * t421;
t114 = -mrSges(6,1) * t656 + mrSges(6,2) * t189;
t95 = mrSges(7,1) * t365 - mrSges(7,3) * t111;
t94 = -mrSges(7,2) * t365 + mrSges(7,3) * t675;
t87 = -pkin(5) * t123 + t194;
t78 = -qJD(6) * t183 - t180 * t415 + t181 * t420;
t77 = qJD(6) * t182 + t180 * t420 + t181 * t415;
t68 = -mrSges(6,2) * t338 + mrSges(6,3) * t75;
t67 = mrSges(6,1) * t338 - mrSges(6,3) * t74;
t56 = -mrSges(7,1) * t675 + mrSges(7,2) * t111;
t50 = -pkin(5) * t75 + t125;
t47 = -qJD(6) * t147 - t122 * t415 + t123 * t420;
t46 = qJD(6) * t146 + t122 * t420 + t123 * t415;
t25 = -mrSges(7,2) * t322 + mrSges(7,3) * t30;
t24 = mrSges(7,1) * t322 - mrSges(7,3) * t29;
t21 = pkin(10) * t123 + t31;
t20 = pkin(5) * t504 - pkin(10) * t122 + t32;
t19 = t420 * t48 - t535;
t18 = -t415 * t48 - t530;
t5 = -qJD(6) * t39 + t20 * t420 - t21 * t415;
t4 = qJD(6) * t38 + t20 * t415 + t21 * t420;
t1 = [(t12 * t225 - t122 * t57 + t123 * t58 - t13 * t226) * mrSges(6,3) + (Ifges(5,5) * t230 + Ifges(5,6) * t229 + qJD(2) * t430 - t450 * t500) * t566 + t631 * mrSges(3,3) + (t663 + t662 + t660 - t635 / 0.2e1 - t490 / 0.2e1 - t489 / 0.2e1 + Ifges(3,4) * t360 / 0.2e1 + (Ifges(3,2) / 0.2e1 + pkin(7) * mrSges(3,3)) * t359 + (-Ifges(3,2) * t418 + t544) * t494 / 0.2e1 - t673 * t574 - Ifges(5,5) * t587 - Ifges(5,6) * t588 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * qJDD(2) - Ifges(6,6) * t604 - Ifges(6,5) * t605 - Ifges(7,6) * t612 - Ifges(7,5) * t613 - Ifges(6,3) * t575 - Ifges(7,3) * t576 - Ifges(4,6) * t581 - Ifges(4,5) * t582) * t423 - (t422 * t239 + t417 * t240) * t500 / 0.2e1 + (t146 * t2 - t147 * t3 - t16 * t46 + t17 * t47) * mrSges(7,3) - t666 * t494 + (-Ifges(5,1) * t311 - Ifges(5,4) * t310) * t587 + (Ifges(5,1) * t230 + Ifges(5,4) * t229) * t577 + (Ifges(6,5) * t122 + Ifges(6,6) * t123) * t568 + (Ifges(6,5) * t226 + Ifges(6,6) * t225) * t575 + (t545 + t455) * t359 / 0.2e1 + (-t133 * t230 + t134 * t229 - t310 * t60 + t311 * t59) * mrSges(5,3) + (Ifges(7,5) * t46 + Ifges(7,6) * t47) * t570 + (Ifges(7,5) * t147 + Ifges(7,6) * t146) * t576 + (Ifges(7,1) * t147 + Ifges(7,4) * t146) * t613 + (Ifges(7,1) * t46 + Ifges(7,4) * t47) * t594 + (Ifges(7,4) * t147 + Ifges(7,2) * t146) * t612 + (Ifges(7,4) * t46 + Ifges(7,2) * t47) * t596 - t367 * t488 + m(7) * (t124 * t87 + t16 * t5 + t17 * t4 + t176 * t50 + t2 * t39 + t3 * t38) + m(6) * (t12 * t91 + t125 * t259 + t13 * t90 + t194 * t198 + t31 * t58 + t32 * t57) + m(5) * (t104 * t133 + t105 * t134 + t192 * t59 + t193 * t60 + t204 * t361 + t272 * t286) + t124 * (-mrSges(7,1) * t47 + mrSges(7,2) * t46) + t204 * (mrSges(5,1) * t310 - mrSges(5,2) * t311) + (Ifges(6,1) * t226 + Ifges(6,4) * t225) * t605 + (Ifges(6,1) * t122 + Ifges(6,4) * t123) * t583 + m(4) * (t160 * t294 + t161 * t293 + t215 * t262 + t216 * t261 + (t320 * t418 + t372 * t502) * pkin(7)) + t354 * (qJD(2) * t431 - t453 * t500) / 0.2e1 + (t405 * mrSges(3,3) + t544 / 0.2e1 + t418 * Ifges(3,1) - pkin(1) * mrSges(3,2)) * t360 + t122 * t598 + t123 * t600 - t311 * t602 - t310 * t603 + t46 * t606 + t47 * t608 + t226 * t610 + t225 * t611 + (qJD(2) * t432 - t456 * t500) * t572 + (t625 / 0.2e1 - t659 + Ifges(5,3) * t566 + Ifges(6,3) * t568 + Ifges(7,3) * t570 + Ifges(5,5) * t577 + Ifges(5,6) * t579 + Ifges(6,5) * t583 + Ifges(6,6) * t585 + Ifges(7,5) * t594 + Ifges(7,6) * t596 - t661) * t504 + (-Ifges(5,4) * t311 - Ifges(5,2) * t310) * t588 + (Ifges(5,4) * t230 + Ifges(5,2) * t229) * t579 + (Ifges(6,4) * t226 + Ifges(6,2) * t225) * t604 + (Ifges(6,4) * t122 + Ifges(6,2) * t123) * t585 + (-t160 * t522 - t161 * t519 - t261 * t434 - t262 * t433) * mrSges(4,3) + t479 * t502 + t230 * t589 + t229 * t591 + t519 * t593 + t463 * t529 + t4 * t94 + t5 * t95 + t90 * t67 + t91 * t68 + t87 * t56 + Ifges(2,3) * qJDD(1) + t38 * t24 + t39 * t25 + qJD(2) ^ 2 * t452 / 0.2e1 + (-t520 * t615 - t329 * mrSges(4,1) - t304 * mrSges(5,1) - t288 * mrSges(6,1) - t274 * mrSges(7,1) - t328 * mrSges(4,2) - t303 * mrSges(5,2) - t287 * mrSges(6,2) - t273 * mrSges(7,2) + (-m(7) * (-pkin(1) - t449) - m(6) * (-pkin(1) - t448) - m(5) * (-pkin(1) - t447) - m(4) * t364 + m(3) * pkin(1) + t622) * t419 + (pkin(7) * t626 + t621) * t424) * g(1) + t372 * (mrSges(4,1) * t433 + mrSges(4,2) * t434) + t361 * t99 + t320 * t340 + t184 * t405 + t147 * t616 + t146 * t617 + t418 * t454 * t581 + t418 * t457 * t582 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t631) + t31 * t162 + t32 * t163 + t632 * t502 / 0.2e1 + t637 * t400 + t176 * t10 + t192 * t145 + t193 * t144 + t194 * t114 + t198 * (-mrSges(6,1) * t123 + mrSges(6,2) * t122) + t105 * t223 + t104 * t224 + t125 * (-mrSges(6,1) * t225 + mrSges(6,2) * t226) + pkin(1) * mrSges(3,1) * t359 + t259 * t37 + t272 * (-mrSges(5,1) * t229 + mrSges(5,2) * t230) + t50 * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) + t215 * t284 + t216 * t285 + t286 * t191 + (-Ifges(5,5) * t311 - Ifges(5,6) * t310 + t418 * t451) * t574 + (-mrSges(3,1) * t405 + Ifges(3,5) * t418) * qJDD(2) + t293 * t206 + t294 * t207 + (-t518 * t615 - t331 * mrSges(4,1) - t306 * mrSges(5,1) - t290 * mrSges(6,1) - t276 * mrSges(7,1) - t330 * mrSges(4,2) - t305 * mrSges(5,2) - t289 * mrSges(6,2) - t275 * mrSges(7,2) + t626 * (t424 * pkin(1) + t419 * pkin(7)) + t621 * t419 + (-m(4) * t465 - m(5) * t447 - m(6) * t448 - m(7) * t449 - t622) * t424) * g(2) - t156 * t522 / 0.2e1; (-t538 + t240 / 0.2e1) * t499 - t625 * t506 / 0.2e1 + (g(3) * t664 + qJD(1) * t661 + t657 * t665) * t418 + (-g(3) * t665 + t657 * t664) * t423 + (t354 * t454 + t355 * t457 + t384 * t451) * qJD(3) / 0.2e1 - (t354 * t431 + t355 * t432 + t384 * t430) * qJD(1) / 0.2e1 + t666 * qJD(1) ^ 2 + ((t261 * t515 + t262 * t521) * qJD(1) + t630) * mrSges(4,3) + (Ifges(6,1) * t583 + Ifges(6,4) * t585 + Ifges(6,5) * t568 + t598 + t627) * t180 + (Ifges(6,4) * t583 + Ifges(6,2) * t585 + Ifges(6,6) * t568 + t600 + t628) * t181 + (t442 + t479) * qJD(3) - t463 * g(3) + t657 * t462 + (Ifges(7,1) * t594 + Ifges(7,4) * t596 + Ifges(7,5) * t570 - t562 + t606) * t77 + (Ifges(6,5) * t210 + Ifges(6,6) * t209) * t569 + (Ifges(7,5) * t137 + Ifges(7,6) * t136) * t571 + (Ifges(5,5) * t578 + Ifges(6,5) * t584 + Ifges(7,5) * t595 + Ifges(5,6) * t580 + Ifges(6,6) * t586 + Ifges(7,6) * t597 + Ifges(5,3) * t567 + Ifges(6,3) * t569 + Ifges(7,3) * t571 + t659) * t506 + t239 * t483 / 0.2e1 + (Ifges(7,1) * t137 + Ifges(7,4) * t136) * t595 + (-pkin(2) * t320 - t261 * t277 - t262 * t278 - t372 * t398) * m(4) + (Ifges(7,4) * t137 + Ifges(7,2) * t136) * t597 - t501 * t537 - (Ifges(5,1) * t577 + Ifges(5,4) * t579 + Ifges(5,5) * t566 - t549 + t589) * t325 + (Ifges(6,1) * t210 + Ifges(6,4) * t209) * t584 - t452 * t494 / 0.2e1 + t210 * t599 + t209 * t601 + t348 * t602 + (Ifges(6,4) * t250 + Ifges(6,2) * t249) * t604 + (Ifges(6,1) * t250 + Ifges(6,4) * t249) * t605 + t137 * t607 + t136 * t609 + t250 * t610 + t249 * t611 + (Ifges(7,4) * t183 + Ifges(7,2) * t182) * t612 + (Ifges(7,1) * t183 + Ifges(7,4) * t182) * t613 + (Ifges(6,5) * t250 + Ifges(6,6) * t249) * t575 + (Ifges(7,5) * t183 + Ifges(7,6) * t182) * t576 + ((-t137 + t77) * mrSges(7,2) + (t136 - t78) * mrSges(7,1)) * t124 + (Ifges(7,4) * t594 + Ifges(7,2) * t596 + Ifges(7,6) * t570 + t561 + t608) * t78 + (-Ifges(5,5) * t292 - Ifges(5,6) * t291) * t567 + (-Ifges(5,4) * t292 - Ifges(5,2) * t291) * t580 + (-Ifges(5,1) * t292 - Ifges(5,4) * t291) * t578 + (-t133 * t292 + t134 * t291 - t348 * t59 - t446 * t60) * mrSges(5,3) - t442 * t505 + (Ifges(6,4) * t210 + Ifges(6,2) * t209) * t586 + t453 * t581 + t456 * t582 - t292 * t590 - t291 * t592 + t417 * t593 + (Ifges(5,5) * t348 - Ifges(5,6) * t446 + t450) * t574 - t446 * t603 + (Ifges(5,1) * t348 - Ifges(5,4) * t446) * t587 + (Ifges(5,4) * t348 - Ifges(5,2) * t446) * t588 + t204 * (mrSges(5,1) * t446 + mrSges(5,2) * t348) + t69 * t24 + t70 * t25 + t320 * t461 + Ifges(3,3) * qJDD(2) + t422 * t156 / 0.2e1 - t393 * t99 + Ifges(3,5) * t360 + Ifges(3,6) * t359 - t344 * mrSges(3,2) - t345 * mrSges(3,1) + t183 * t616 + t182 * t617 + (-t284 * t501 + m(4) * ((-t261 * t422 - t262 * t417) * qJD(3) + t630) - t285 * t499 + t422 * t207 - t417 * t206) * pkin(8) - (-Ifges(3,2) * t506 + t396 + t632) * t505 / 0.2e1 - t637 * t398 + t638 * t191 + t639 * t114 + t640 * t223 + t641 * t224 + (t133 * t641 + t134 * t640 - t204 * t393 + t266 * t59 + t267 * t60 + t272 * t638) * m(5) + t642 * t56 + (t12 * t249 - t13 * t250 - t209 * t58 + t210 * t57) * mrSges(6,3) + t50 * (-mrSges(7,1) * t182 + mrSges(7,2) * t183) - pkin(2) * t184 + t367 * t397 + t203 * t10 + (t677 * mrSges(5,1) + mrSges(5,2) * t684) * t272 - t198 * (-mrSges(6,1) * t209 + mrSges(6,2) * t210) + t125 * (-mrSges(6,1) * t249 + mrSges(6,2) * t250) + t650 * t162 + t651 * t163 + (t12 * t155 + t125 * t295 + t13 * t154 + t198 * t639 + t57 * t651 + t58 * t650) * m(6) + t652 * t94 + t653 * t95 + (t124 * t642 + t16 * t653 + t17 * t652 + t2 * t70 + t203 * t50 + t3 * t69) * m(7) + t266 * t145 + t267 * t144 - t278 * t284 - t277 * t285 + t295 * t37 + (-t136 * t17 + t137 * t16 + t182 * t2 - t183 * t3) * mrSges(7,3) - (Ifges(5,4) * t577 + Ifges(5,2) * t579 + Ifges(5,6) * t566 + t548 + t591) * t324 + t154 * t67 + t155 * t68; -(Ifges(6,4) * t584 + Ifges(6,2) * t586 + Ifges(6,6) * t569 + t601 - t628) * t189 - (-Ifges(4,2) * t355 + t240 + t343) * t354 / 0.2e1 - t660 + (Ifges(6,1) * t584 + Ifges(6,4) * t586 + Ifges(6,5) * t569 + t599 - t627) * t656 + t635 - t191 * t560 - m(5) * (t133 * t141 + t134 * t142 + t272 * t560) - t355 * (Ifges(4,1) * t354 - t536) / 0.2e1 + t355 * t537 + (Ifges(4,5) * t354 - Ifges(4,6) * t355) * t567 + t239 * t572 + (-mrSges(5,2) * t272 + Ifges(5,1) * t578 + Ifges(5,4) * t580 + Ifges(5,5) * t567 + t549 + t590) * t466 - (mrSges(5,1) * t272 + Ifges(5,4) * t578 + Ifges(5,2) * t580 + Ifges(5,6) * t567 - t548 + t592) * t256 + t681 - t372 * (mrSges(4,1) * t355 + mrSges(4,2) * t354) + g(3) * t340 + t354 * t538 + t145 * t558 + t144 * t559 + (t412 * t60 + t413 * t59) * t615 + (m(5) * t557 + mrSges(5,1) * t401 + mrSges(5,2) * t402 - t624 + t636) * t552 - t208 * t114 + t645 * t162 + t646 * t163 + (t12 * t316 + t13 * t315 - t198 * t208 + t57 * t646 + t58 * t645) * m(6) + (t303 * mrSges(5,1) - t304 * mrSges(5,2) - m(6) * (-t362 * t517 - t363 * t424) - m(7) * (-t308 * t517 - t309 * t424) - mrSges(4,2) * t329 + t647 * t328 + t633) * g(2) + (-t305 * mrSges(5,1) + t306 * mrSges(5,2) - m(6) * (-t362 * t514 + t419 * t363) - m(7) * (-t308 * t514 + t419 * t309) + mrSges(4,2) * t331 - t647 * t330 + t634) * g(1) + t648 * t94 + t649 * t95 + (-t124 * t130 + t16 * t649 + t17 * t648 + t2 * t228 + t227 * t3) * m(7) - t142 * t223 - t141 * t224 + t227 * t24 + t228 * t25 - t130 * t56 - t261 * t284 + t262 * t285 + t315 * t67 + t316 * t68; t111 * t95 - t675 * t94 - t656 * t162 + t189 * t163 - t466 * t223 + t256 * t224 + t10 + t37 + t99 + (t111 * t16 - t17 * t675 + t50) * m(7) + (t189 * t57 - t58 * t656 + t125) * m(6) + (t133 * t256 - t134 * t466 + t204) * m(5) + (t423 * g(3) - t418 * t657) * t629; -t56 * t556 - m(7) * (t124 * t556 + t16 * t18 + t17 * t19) + (t2 * t415 + t3 * t420 + (-t16 * t415 + t17 * t420) * qJD(6)) * t614 + (Ifges(6,5) * t656 - Ifges(6,6) * t189) * t569 + t97 * t583 + (Ifges(6,1) * t656 - t541) * t584 - t19 * t94 - t18 * t95 - t198 * (mrSges(6,1) * t189 + mrSges(6,2) * t656) + (-Ifges(6,2) * t189 + t179 + t98) * t586 + (t546 + t163) * t58 + (t547 - t162) * t57 + (m(7) * t555 + t636) * t552 + (t287 * t614 + t633) * g(2) + (-t289 * t614 + t634) * g(1) + ((-t415 * t95 + t420 * t94) * qJD(6) + t420 * t24 + t415 * t25) * pkin(5) + t681; -t124 * (mrSges(7,1) * t111 + mrSges(7,2) * t675) + (Ifges(7,1) * t675 - t540) * t595 + t54 * t594 + (Ifges(7,5) * t675 - Ifges(7,6) * t111) * t571 - t16 * t94 + t17 * t95 - g(1) * t512 - g(2) * t513 - t458 * t552 + (t111 * t17 + t16 * t675) * mrSges(7,3) + t445 + (-Ifges(7,2) * t111 + t102 + t55) * t597;];
tau  = t1;
