% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR11
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:20
% EndTime: 2019-03-09 07:39:50
% DurationCPUTime: 53.25s
% Computational Cost: add. (48772->1093), mult. (154245->1487), div. (0->0), fcn. (134155->18), ass. (0->489)
t405 = sin(pkin(6));
t408 = cos(pkin(6));
t403 = sin(pkin(13));
t407 = cos(pkin(7));
t412 = sin(qJ(3));
t406 = cos(pkin(13));
t595 = cos(qJ(3));
t507 = t595 * t406;
t441 = -t403 * t412 + t407 * t507;
t404 = sin(pkin(7));
t514 = t404 * t595;
t425 = t405 * t441 + t408 * t514;
t301 = t425 * qJD(1);
t296 = -t301 + qJD(4);
t547 = t407 * t412;
t434 = t405 * (-t403 * t547 + t507);
t333 = qJD(1) * t434;
t498 = qJD(3) * t595;
t479 = t404 * t498;
t722 = t479 - t333;
t535 = qJD(1) * t405;
t504 = t406 * t535;
t591 = pkin(1) * t408;
t524 = qJD(1) * t591;
t346 = qJ(2) * t504 + t403 * t524;
t548 = t405 * t407;
t551 = t404 * t408;
t440 = (t406 * t548 + t551) * pkin(9);
t290 = qJD(1) * t440 + t346;
t587 = pkin(9) * t403;
t335 = (-pkin(2) * t406 - t404 * t587 - pkin(1)) * t405;
t322 = qJD(1) * t335 + qJD(2);
t387 = t406 * t524;
t553 = t403 * t405;
t590 = pkin(2) * t408;
t435 = t590 + (-pkin(9) * t407 - qJ(2)) * t553;
t297 = qJD(1) * t435 + t387;
t557 = t297 * t407;
t201 = t595 * t290 + (t322 * t404 + t557) * t412;
t411 = sin(qJ(4));
t415 = cos(qJ(4));
t725 = -t201 + t296 * (pkin(4) * t411 - pkin(11) * t415);
t508 = t595 * t403;
t443 = t406 * t547 + t508;
t550 = t404 * t412;
t313 = t405 * t443 + t408 * t550;
t304 = t313 * qJD(1);
t414 = cos(qJ(5));
t410 = sin(qJ(5));
t544 = t410 * t415;
t229 = -t301 * t544 + t304 * t414;
t528 = qJD(5) * t414;
t531 = qJD(4) * t415;
t670 = t410 * t531 + t411 * t528 + t229;
t358 = -t415 * t407 + t411 * t550;
t505 = t403 * t535;
t481 = t404 * t505;
t666 = -qJD(4) * t358 - t411 * t481 + t415 * t722;
t433 = t405 * (t406 * t412 + t407 * t508);
t332 = qJD(1) * t433;
t533 = qJD(3) * t412;
t502 = t404 * t533;
t724 = t502 - t332;
t511 = t407 * t595;
t487 = t297 * t511;
t200 = -t412 * t290 + t322 * t514 + t487;
t240 = pkin(3) * t304 - pkin(10) * t301;
t142 = t415 * t200 + t411 * t240;
t121 = pkin(11) * t304 + t142;
t378 = -pkin(4) * t415 - pkin(11) * t411 - pkin(3);
t530 = qJD(5) * t410;
t532 = qJD(4) * t411;
t684 = -t414 * t121 + t378 * t528 + (-t414 * t532 - t415 * t530) * pkin(10) + t725 * t410;
t522 = pkin(10) * t532;
t723 = t725 * t414 + (t121 + t522) * t410;
t594 = sin(qJ(1));
t396 = t594 * t403;
t596 = cos(qJ(1));
t509 = t596 * t406;
t456 = -t408 * t509 + t396;
t513 = t405 * t596;
t721 = t404 * t513 + t407 * t456;
t506 = t594 * t406;
t510 = t596 * t403;
t352 = t408 * t510 + t506;
t273 = -t352 * t595 + t412 * t721;
t317 = t456 * t404 - t407 * t513;
t720 = t273 * t415 - t317 * t411;
t719 = t273 * t411 + t317 * t415;
t541 = t414 * t415;
t230 = t301 * t541 + t304 * t410;
t395 = pkin(10) * t541;
t556 = t301 * t411;
t718 = -pkin(5) * t556 + pkin(12) * t230 + (pkin(5) * t411 - pkin(12) * t541) * qJD(4) + (-t395 + (pkin(12) * t411 - t378) * t410) * qJD(5) + t723;
t717 = -pkin(12) * t670 + t684;
t416 = -pkin(12) - pkin(11);
t515 = qJD(5) * t416;
t549 = t405 * t406;
t452 = t404 * t549 - t407 * t408;
t341 = -qJD(1) * t452 + qJD(3);
t254 = -t411 * t304 + t341 * t415;
t559 = t254 * t410;
t255 = t304 * t415 + t341 * t411;
t190 = pkin(4) * t255 - pkin(11) * t254;
t245 = -t297 * t404 + t407 * t322;
t175 = -pkin(3) * t301 - pkin(10) * t304 + t245;
t178 = t341 * pkin(10) + t201;
t95 = t175 * t415 - t411 * t178;
t78 = t410 * t190 + t414 * t95;
t716 = pkin(12) * t559 + t410 * t515 - t78;
t558 = t254 * t414;
t77 = t414 * t190 - t410 * t95;
t715 = -pkin(5) * t255 + pkin(12) * t558 + t414 * t515 - t77;
t359 = t407 * t411 + t415 * t550;
t442 = -t414 * t359 + t410 * t514;
t668 = qJD(5) * t442 - t410 * t666 + t414 * t724;
t323 = -t410 * t359 - t414 * t514;
t714 = -qJD(5) * t323 - t410 * t724 - t414 * t666;
t398 = pkin(5) * t414 + pkin(4);
t402 = qJ(5) + qJ(6);
t399 = sin(t402);
t400 = cos(t402);
t471 = -mrSges(6,1) * t414 + mrSges(6,2) * t410;
t662 = m(6) * pkin(4) + m(7) * t398 + mrSges(7,1) * t400 - mrSges(7,2) * t399 - t471;
t660 = -m(6) * pkin(11) + m(7) * t416 - mrSges(6,3) - mrSges(7,3);
t713 = t530 - t559;
t527 = qJD(1) * qJD(3);
t238 = (qJD(1) * t498 + qJDD(1) * t412) * t551 + (qJDD(1) * t443 + t441 * t527) * t405;
t340 = -qJDD(1) * t452 + qJDD(3);
t158 = -t411 * t238 - t304 * t531 + t340 * t415 - t341 * t532;
t153 = qJDD(5) - t158;
t150 = qJDD(6) + t153;
t621 = t150 / 0.2e1;
t204 = -t255 * t410 + t296 * t414;
t205 = t255 * t414 + t296 * t410;
t409 = sin(qJ(6));
t413 = cos(qJ(6));
t133 = t204 * t409 + t205 * t413;
t157 = t254 * qJD(4) + t238 * t415 + t340 * t411;
t239 = t405 * (qJDD(1) * t441 - t443 * t527) - (-qJDD(1) * t595 + t412 * t527) * t551;
t236 = qJDD(4) - t239;
t85 = qJD(5) * t204 + t157 * t414 + t236 * t410;
t86 = -qJD(5) * t205 - t157 * t410 + t236 * t414;
t34 = -qJD(6) * t133 - t409 * t85 + t413 * t86;
t636 = t34 / 0.2e1;
t491 = t413 * t204 - t205 * t409;
t33 = qJD(6) * t491 + t409 * t86 + t413 * t85;
t637 = t33 / 0.2e1;
t638 = Ifges(7,1) * t637 + Ifges(7,4) * t636 + Ifges(7,5) * t621;
t639 = Ifges(7,4) * t637 + Ifges(7,2) * t636 + Ifges(7,6) * t621;
t712 = mrSges(4,2) - mrSges(5,3);
t177 = -t341 * pkin(3) - t200;
t711 = t177 * mrSges(5,2);
t710 = t296 * Ifges(5,5);
t709 = t296 * Ifges(5,6);
t251 = qJD(5) - t254;
t247 = qJD(6) + t251;
t689 = t205 * Ifges(6,5) + t133 * Ifges(7,5) + t204 * Ifges(6,6) + Ifges(7,6) * t491 + t251 * Ifges(6,3) + t247 * Ifges(7,3);
t580 = mrSges(4,3) * t304;
t708 = -mrSges(4,1) * t341 - mrSges(5,1) * t254 + mrSges(5,2) * t255 + t580;
t665 = qJD(4) * t359 + t411 * t722 + t415 * t481;
t436 = t408 * t506 + t510;
t512 = t405 * t594;
t422 = t436 * t404 + t407 * t512;
t470 = t410 * mrSges(6,1) + t414 * mrSges(6,2);
t516 = pkin(5) * t410 + pkin(10);
t700 = m(7) * t516;
t707 = -t399 * mrSges(7,1) - t400 * mrSges(7,2) - t470 - t700;
t141 = -t411 * t200 + t240 * t415;
t120 = -pkin(4) * t304 - t141;
t521 = pkin(10) * t531;
t706 = t521 - t120;
t473 = -mrSges(5,1) * t415 + mrSges(5,2) * t411;
t652 = -t411 * t660 + t415 * t662 - t473;
t701 = m(6) + m(7);
t705 = mrSges(4,1) + t652 + pkin(3) * (m(5) + t701);
t270 = t352 * t412 + t595 * t721;
t111 = -t254 * pkin(4) - t255 * pkin(11) + t177;
t363 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t405;
t525 = qJDD(1) * t408;
t519 = pkin(1) * t525;
t327 = t406 * t363 + t403 * t519;
t279 = qJDD(1) * t440 + t327;
t326 = -t363 * t403 + t406 * t519;
t282 = (-t548 * t587 + t590) * qJDD(1) + t326;
t318 = qJDD(1) * t335 + qJDD(2);
t112 = qJD(3) * t487 + t595 * t279 + t282 * t547 - t290 * t533 + t318 * t550 + t322 * t479;
t107 = pkin(10) * t340 + t112;
t228 = -t282 * t404 + t407 * t318;
t128 = -pkin(3) * t239 - pkin(10) * t238 + t228;
t41 = t415 * t107 + t411 * t128 + t175 * t531 - t178 * t532;
t36 = pkin(11) * t236 + t41;
t113 = -t412 * t279 + t282 * t511 - t290 * t498 + t318 * t514 - t322 * t502 - t533 * t557;
t108 = -t340 * pkin(3) - t113;
t54 = -t158 * pkin(4) - t157 * pkin(11) + t108;
t96 = t175 * t411 + t178 * t415;
t89 = pkin(11) * t296 + t96;
t8 = t111 * t528 + t414 * t36 + t410 * t54 - t530 * t89;
t56 = t111 * t410 + t414 * t89;
t9 = -qJD(5) * t56 - t36 * t410 + t414 * t54;
t704 = t9 * mrSges(6,1) - t8 * mrSges(6,2);
t55 = t414 * t111 - t410 * t89;
t46 = -pkin(12) * t205 + t55;
t43 = pkin(5) * t251 + t46;
t47 = pkin(12) * t204 + t56;
t564 = t409 * t47;
t16 = t413 * t43 - t564;
t6 = pkin(5) * t153 - pkin(12) * t85 + t9;
t7 = pkin(12) * t86 + t8;
t2 = qJD(6) * t16 + t409 * t6 + t413 * t7;
t561 = t413 * t47;
t17 = t409 * t43 + t561;
t3 = -qJD(6) * t17 - t409 * t7 + t413 * t6;
t703 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t611 = t236 / 0.2e1;
t617 = t158 / 0.2e1;
t618 = t157 / 0.2e1;
t629 = Ifges(5,1) * t618 + Ifges(5,4) * t617 + Ifges(5,5) * t611;
t702 = mrSges(5,1) * t177 + t55 * mrSges(6,1) + mrSges(7,1) * t16 - t56 * mrSges(6,2) - mrSges(7,2) * t17;
t628 = t85 / 0.2e1;
t627 = t86 / 0.2e1;
t620 = t153 / 0.2e1;
t699 = t95 * mrSges(5,1);
t698 = t96 * mrSges(5,2);
t10 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t150;
t38 = Ifges(6,5) * t85 + Ifges(6,6) * t86 + Ifges(6,3) * t153;
t697 = t38 + t10;
t295 = Ifges(4,4) * t301;
t696 = Ifges(4,5) * t238;
t695 = Ifges(4,6) * t239;
t694 = Ifges(4,3) * t340;
t693 = t296 * Ifges(5,3);
t692 = t301 * Ifges(4,2);
t691 = t341 * Ifges(4,5);
t690 = t341 * Ifges(4,6);
t366 = t414 * t378;
t543 = t411 * t414;
t311 = -pkin(12) * t543 + t366 + (-pkin(10) * t410 - pkin(5)) * t415;
t343 = t410 * t378 + t395;
t545 = t410 * t411;
t325 = -pkin(12) * t545 + t343;
t248 = t311 * t413 - t325 * t409;
t688 = qJD(6) * t248 + t409 * t718 + t717 * t413;
t249 = t311 * t409 + t325 * t413;
t687 = -qJD(6) * t249 - t717 * t409 + t413 * t718;
t382 = t416 * t410;
t383 = t416 * t414;
t330 = t382 * t413 + t383 * t409;
t686 = qJD(6) * t330 + t409 * t715 + t413 * t716;
t331 = t382 * t409 - t383 * t413;
t685 = -qJD(6) * t331 - t409 * t716 + t413 * t715;
t683 = -qJD(5) * t343 + t723;
t682 = pkin(5) * t670 + t706;
t116 = mrSges(5,1) * t236 - mrSges(5,3) * t157;
t45 = -mrSges(6,1) * t86 + mrSges(6,2) * t85;
t681 = t45 - t116;
t680 = pkin(5) * t713 - t96;
t640 = m(7) * pkin(5);
t679 = -t640 - mrSges(6,1);
t458 = t409 * t410 - t413 * t414;
t348 = t458 * t411;
t159 = t229 * t413 - t230 * t409;
t368 = t409 * t414 + t410 * t413;
t655 = qJD(5) + qJD(6);
t261 = t348 * t655 - t368 * t531;
t675 = t159 - t261;
t160 = t229 * t409 + t230 * t413;
t316 = t655 * t368;
t260 = -t316 * t411 - t458 * t531;
t674 = t160 - t260;
t253 = t323 * t409 - t413 * t442;
t673 = -qJD(6) * t253 + t409 * t714 + t413 * t668;
t252 = t323 * t413 + t409 * t442;
t672 = qJD(6) * t252 + t409 * t668 - t413 * t714;
t136 = -mrSges(6,1) * t204 + mrSges(6,2) * t205;
t578 = mrSges(5,3) * t255;
t207 = mrSges(5,1) * t296 - t578;
t671 = t207 - t136;
t529 = qJD(5) * t411;
t669 = t410 * t529 - t414 * t531 + t230;
t355 = qJ(2) * t549 + t403 * t591;
t305 = t440 + t355;
t394 = t406 * t591;
t314 = t394 + t435;
t208 = -t412 * t305 + t314 * t511 + t335 * t514;
t664 = -t404 * t512 + t436 * t407;
t663 = -m(6) * pkin(10) + t707;
t658 = t532 - t556;
t42 = -t411 * t107 + t128 * t415 - t175 * t532 - t178 * t531;
t657 = t41 * t415 - t411 * t42;
t656 = -t410 * t9 + t414 * t8;
t654 = mrSges(5,1) + t662;
t641 = mrSges(5,2) + t660;
t651 = -g(1) * t512 + g(2) * t513 - g(3) * t408;
t88 = -pkin(4) * t296 - t95;
t76 = -pkin(5) * t204 + t88;
t650 = -mrSges(7,1) * t76 + mrSges(7,3) * t17;
t649 = mrSges(7,2) * t76 - mrSges(7,3) * t16;
t646 = -m(5) * t177 - t708;
t645 = t42 * mrSges(5,1) - t41 * mrSges(5,2);
t644 = -t113 * mrSges(4,1) + t112 * mrSges(4,2);
t643 = -m(5) * pkin(10) + mrSges(4,2) + t663;
t39 = Ifges(6,4) * t85 + Ifges(6,2) * t86 + Ifges(6,6) * t153;
t635 = t39 / 0.2e1;
t634 = Ifges(6,1) * t628 + Ifges(6,4) * t627 + Ifges(6,5) * t620;
t570 = Ifges(7,4) * t133;
t68 = Ifges(7,2) * t491 + t247 * Ifges(7,6) + t570;
t633 = -t68 / 0.2e1;
t632 = t68 / 0.2e1;
t129 = Ifges(7,4) * t491;
t69 = t133 * Ifges(7,1) + t247 * Ifges(7,5) + t129;
t631 = -t69 / 0.2e1;
t630 = t69 / 0.2e1;
t203 = Ifges(6,4) * t204;
t106 = t205 * Ifges(6,1) + t251 * Ifges(6,5) + t203;
t626 = t106 / 0.2e1;
t625 = -t491 / 0.2e1;
t624 = t491 / 0.2e1;
t623 = -t133 / 0.2e1;
t622 = t133 / 0.2e1;
t567 = t255 * Ifges(5,4);
t155 = t254 * Ifges(5,2) + t567 + t709;
t619 = -t155 / 0.2e1;
t616 = -t204 / 0.2e1;
t615 = t204 / 0.2e1;
t614 = -t205 / 0.2e1;
t613 = t205 / 0.2e1;
t610 = -t247 / 0.2e1;
t609 = t247 / 0.2e1;
t608 = -t251 / 0.2e1;
t607 = t251 / 0.2e1;
t606 = -t254 / 0.2e1;
t605 = t254 / 0.2e1;
t604 = -t255 / 0.2e1;
t603 = t255 / 0.2e1;
t601 = -t296 / 0.2e1;
t598 = t304 / 0.2e1;
t589 = pkin(5) * t205;
t269 = t313 * t415 - t411 * t452;
t218 = -t269 * t410 - t414 * t425;
t588 = pkin(5) * t218;
t579 = mrSges(5,3) * t254;
t577 = mrSges(6,3) * t204;
t576 = mrSges(6,3) * t205;
t575 = Ifges(5,4) * t411;
t574 = Ifges(5,4) * t415;
t573 = Ifges(6,4) * t205;
t572 = Ifges(6,4) * t410;
t571 = Ifges(6,4) * t414;
t568 = t200 * mrSges(4,3);
t566 = t304 * Ifges(4,4);
t37 = -pkin(4) * t236 - t42;
t565 = t37 * t411;
t257 = -t314 * t404 + t407 * t335;
t309 = t425 * pkin(3);
t492 = t313 * pkin(10) + t309;
t192 = t257 - t492;
t294 = t595 * t305;
t209 = t314 * t547 + t335 * t550 + t294;
t199 = -pkin(10) * t452 + t209;
t115 = t411 * t192 + t415 * t199;
t101 = -pkin(11) * t425 + t115;
t198 = pkin(3) * t452 - t208;
t268 = t313 * t411 + t415 * t452;
t125 = t268 * pkin(4) - t269 * pkin(11) + t198;
t62 = t414 * t101 + t410 * t125;
t555 = t301 * t415;
t105 = t204 * Ifges(6,2) + t251 * Ifges(6,6) + t573;
t546 = t410 * t105;
t542 = t414 * t106;
t540 = (t270 * t400 + t399 * t720) * mrSges(7,1) + (-t270 * t399 + t400 * t720) * mrSges(7,2);
t353 = -t396 * t408 + t509;
t275 = t353 * t595 - t412 * t664;
t225 = t275 * t415 + t411 * t422;
t274 = t353 * t412 + t595 * t664;
t165 = -t225 * t399 + t274 * t400;
t166 = t225 * t400 + t274 * t399;
t539 = t165 * mrSges(7,1) - t166 * mrSges(7,2);
t538 = (-t269 * t399 - t400 * t425) * mrSges(7,1) + (-t269 * t400 + t399 * t425) * mrSges(7,2);
t536 = t596 * pkin(1) + qJ(2) * t512;
t534 = qJD(2) * t405;
t526 = qJDD(1) * t405;
t518 = Ifges(5,5) * t157 + Ifges(5,6) * t158 + Ifges(5,3) * t236;
t517 = t694 + t695 + t696;
t503 = t403 * t534;
t500 = -t546 / 0.2e1;
t497 = t403 * t526;
t496 = t406 * t526;
t494 = t531 / 0.2e1;
t493 = -t529 / 0.2e1;
t61 = -t101 * t410 + t414 * t125;
t114 = t192 * t415 - t411 * t199;
t480 = t404 * t503;
t478 = -pkin(1) * t594 + qJ(2) * t513;
t475 = -mrSges(4,1) * t425 + mrSges(4,2) * t313;
t474 = mrSges(5,1) * t268 + mrSges(5,2) * t269;
t219 = t269 * t414 - t410 * t425;
t472 = mrSges(6,1) * t218 - mrSges(6,2) * t219;
t469 = Ifges(5,1) * t415 - t575;
t468 = Ifges(6,1) * t414 - t572;
t467 = Ifges(6,1) * t410 + t571;
t466 = -Ifges(5,2) * t411 + t574;
t465 = -Ifges(6,2) * t410 + t571;
t464 = Ifges(6,2) * t414 + t572;
t463 = Ifges(5,5) * t415 - Ifges(5,6) * t411;
t462 = Ifges(6,5) * t414 - Ifges(6,6) * t410;
t461 = Ifges(6,5) * t410 + Ifges(6,6) * t414;
t48 = pkin(5) * t268 - pkin(12) * t219 + t61;
t53 = pkin(12) * t218 + t62;
t22 = -t409 * t53 + t413 * t48;
t23 = t409 * t48 + t413 * t53;
t144 = t218 * t413 - t219 * t409;
t145 = t218 * t409 + t219 * t413;
t170 = -t225 * t410 + t274 * t414;
t459 = -(-qJ(2) * t505 + t387) * t403 + t346 * t406;
t185 = qJD(2) * t434 + qJD(3) * t208;
t302 = t425 * qJD(3);
t303 = t313 * qJD(3);
t227 = pkin(3) * t303 - pkin(10) * t302 + t480;
t72 = -t411 * t185 - t192 * t532 - t199 * t531 + t227 * t415;
t457 = t10 + t703;
t455 = -mrSges(3,1) * t496 + mrSges(3,2) * t497;
t454 = mrSges(3,1) * t408 - mrSges(3,3) * t553;
t453 = -mrSges(3,2) * t408 + mrSges(3,3) * t549;
t100 = pkin(4) * t425 - t114;
t451 = t88 * t470;
t71 = t415 * t185 + t192 * t531 - t199 * t532 + t411 * t227;
t65 = pkin(11) * t303 + t71;
t186 = qJD(2) * t433 + (t294 + (t314 * t407 + t335 * t404) * t412) * qJD(3);
t216 = qJD(4) * t269 + t302 * t411;
t217 = -qJD(4) * t268 + t302 * t415;
t92 = t216 * pkin(4) - t217 * pkin(11) + t186;
t20 = -t101 * t530 + t125 * t528 + t410 * t92 + t414 * t65;
t66 = -pkin(4) * t303 - t72;
t21 = -qJD(5) * t62 - t410 * t65 + t414 * t92;
t423 = -t352 * pkin(2) - pkin(9) * t317 + t478;
t421 = t353 * pkin(2) + pkin(9) * t422 + t536;
t420 = t273 * pkin(3) + t423;
t419 = t275 * pkin(3) + t421;
t417 = t274 * pkin(10) + t419;
t388 = -pkin(1) * t526 + qJDD(2);
t375 = t516 * t411;
t357 = t453 * qJD(1);
t356 = t454 * qJD(1);
t354 = -qJ(2) * t553 + t394;
t347 = t368 * t411;
t342 = -pkin(10) * t544 + t366;
t258 = -mrSges(4,2) * t341 + mrSges(4,3) * t301;
t250 = Ifges(5,4) * t254;
t237 = -mrSges(4,1) * t301 + mrSges(4,2) * t304;
t224 = t275 * t411 - t415 * t422;
t215 = t304 * Ifges(4,1) + t295 + t691;
t214 = t566 + t690 + t692;
t211 = -mrSges(4,2) * t340 + mrSges(4,3) * t239;
t210 = mrSges(4,1) * t340 - mrSges(4,3) * t238;
t206 = -mrSges(5,2) * t296 + t579;
t180 = t458 * t254;
t179 = t368 * t254;
t171 = t225 * t414 + t274 * t410;
t169 = -mrSges(4,1) * t239 + mrSges(4,2) * t238;
t156 = t255 * Ifges(5,1) + t250 + t710;
t154 = t255 * Ifges(5,5) + t254 * Ifges(5,6) + t693;
t147 = mrSges(6,1) * t251 - t576;
t146 = -mrSges(6,2) * t251 + t577;
t127 = qJD(5) * t218 + t217 * t414 + t303 * t410;
t126 = -qJD(5) * t219 - t217 * t410 + t303 * t414;
t117 = -mrSges(5,2) * t236 + mrSges(5,3) * t158;
t99 = mrSges(7,1) * t247 - mrSges(7,3) * t133;
t98 = -mrSges(7,2) * t247 + mrSges(7,3) * t491;
t87 = -mrSges(5,1) * t158 + mrSges(5,2) * t157;
t81 = t100 - t588;
t79 = t157 * Ifges(5,4) + t158 * Ifges(5,2) + t236 * Ifges(5,6);
t73 = -mrSges(7,1) * t491 + mrSges(7,2) * t133;
t60 = -mrSges(6,2) * t153 + mrSges(6,3) * t86;
t59 = mrSges(6,1) * t153 - mrSges(6,3) * t85;
t50 = -qJD(6) * t145 + t126 * t413 - t127 * t409;
t49 = qJD(6) * t144 + t126 * t409 + t127 * t413;
t44 = -pkin(5) * t126 + t66;
t26 = -mrSges(7,2) * t150 + mrSges(7,3) * t34;
t25 = mrSges(7,1) * t150 - mrSges(7,3) * t33;
t24 = -pkin(5) * t86 + t37;
t19 = t413 * t46 - t564;
t18 = -t409 * t46 - t561;
t15 = pkin(12) * t126 + t20;
t14 = pkin(5) * t216 - pkin(12) * t127 + t21;
t13 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t5 = -qJD(6) * t23 + t14 * t413 - t15 * t409;
t4 = qJD(6) * t22 + t14 * t409 + t15 * t413;
t1 = [(Ifges(7,1) * t145 + Ifges(7,4) * t144) * t637 + (Ifges(7,1) * t49 + Ifges(7,4) * t50) * t622 + (Ifges(6,1) * t127 + Ifges(6,4) * t126) * t613 + (Ifges(6,1) * t219 + Ifges(6,4) * t218) * t628 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t636 + (Ifges(7,4) * t49 + Ifges(7,2) * t50) * t624 + (t126 * t56 - t127 * t55 + t218 * t8 - t219 * t9) * mrSges(6,3) + (Ifges(7,5) * t49 + Ifges(7,6) * t50) * t609 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t621 + (-mrSges(4,3) * t113 + Ifges(4,1) * t238 + Ifges(4,4) * t239 + Ifges(4,5) * t340) * t313 - t37 * t472 + t108 * t474 + t228 * t475 + (Ifges(6,4) * t127 + Ifges(6,2) * t126) * t615 + (Ifges(6,4) * t219 + Ifges(6,2) * t218) * t627 + (-mrSges(5,3) * t42 + 0.2e1 * t629) * t269 + m(4) * (t112 * t209 + t113 * t208 + t185 * t201 + t228 * t257 + t245 * t480) + (t144 * t2 - t145 * t3 - t16 * t49 + t17 * t50) * mrSges(7,3) + t406 * t357 * t534 + t237 * t480 + m(6) * (t100 * t37 + t20 * t56 + t21 * t55 + t61 * t9 + t62 * t8 + t66 * t88) + m(7) * (t16 * t5 + t17 * t4 + t2 * t23 + t22 * t3 + t24 * t81 + t44 * t76) - t356 * t503 + m(3) * (t326 * t354 + t327 * t355) + ((Ifges(3,5) * t403 + Ifges(3,6) * t406) * t525 + m(3) * (-pkin(1) * t388 + qJD(2) * t459) + (Ifges(3,1) * t403 + Ifges(3,4) * t406) * t497 + (Ifges(3,4) * t403 + Ifges(3,2) * t406) * t496 - pkin(1) * t455 + t388 * (-mrSges(3,1) * t406 + mrSges(3,2) * t403)) * t405 + (Ifges(3,5) * t497 + Ifges(3,6) * t496 + Ifges(3,3) * t525) * t408 + (Ifges(6,5) * t127 + Ifges(6,6) * t126) * t607 + (Ifges(6,5) * t219 + Ifges(6,6) * t218) * t620 + (t354 * t454 + t355 * t453 + Ifges(2,3)) * qJDD(1) + m(5) * (t108 * t198 + t114 * t42 + t115 * t41 + t71 * t96 + t72 * t95) + t327 * t453 + t326 * t454 - (t518 / 0.2e1 - Ifges(4,6) * t340 - Ifges(4,4) * t238 - Ifges(4,2) * t239 - t112 * mrSges(4,3) + Ifges(5,3) * t611 + Ifges(5,6) * t617 + Ifges(5,5) * t618 + t645) * t425 + t22 * t25 + t23 * t26 + (-Ifges(4,4) * t598 + Ifges(5,5) * t603 + Ifges(5,6) * t605 - t690 / 0.2e1 + t245 * mrSges(4,1) - t214 / 0.2e1 + t693 / 0.2e1 - t698 + t699 + t154 / 0.2e1 - t692 / 0.2e1 - mrSges(4,3) * t201) * t303 + (Ifges(4,1) * t598 - t568 + t691 / 0.2e1 + t245 * mrSges(4,2) + t215 / 0.2e1 + t295 / 0.2e1) * t302 + (-t517 / 0.2e1 - t694 / 0.2e1 - t696 / 0.2e1 - t695 / 0.2e1 + t644) * t452 + t185 * t258 + t257 * t169 + (-m(4) * t200 - t646) * t186 + t208 * t210 + t209 * t211 + t71 * t206 + t72 * t207 + (-m(3) * t478 - m(4) * t423 - m(7) * t420 + mrSges(2,1) * t594 + t352 * mrSges(3,1) - t273 * mrSges(4,1) + mrSges(2,2) * t596 - mrSges(3,2) * t456 - mrSges(3,3) * t513 + mrSges(4,3) * t317 + (-m(6) - m(5)) * (-pkin(10) * t270 + t420) - t654 * t720 - (t707 + t712) * t270 + t641 * t719) * g(1) + t198 * t87 + t50 * t632 + t219 * t634 + t218 * t635 + t145 * t638 + t144 * t639 + t20 * t146 + t21 * t147 + t24 * (-mrSges(7,1) * t144 + mrSges(7,2) * t145) + t66 * t136 + t126 * t105 / 0.2e1 + t88 * (-mrSges(6,1) * t126 + mrSges(6,2) * t127) + t114 * t116 + t115 * t117 + (-t96 * mrSges(5,3) + Ifges(6,3) * t607 - Ifges(5,2) * t605 - t709 / 0.2e1 - Ifges(5,4) * t603 + t689 / 0.2e1 + Ifges(7,3) * t609 + Ifges(6,5) * t613 + Ifges(6,6) * t615 + t619 + Ifges(7,5) * t622 + Ifges(7,6) * t624 + t702) * t216 + (-t95 * mrSges(5,3) + Ifges(5,4) * t605 + t710 / 0.2e1 + Ifges(5,1) * t603 + t711 + t156 / 0.2e1) * t217 + t61 * t59 + t62 * t60 + (-m(7) * (t225 * t398 + t419) - t166 * mrSges(7,1) - t165 * mrSges(7,2) - mrSges(2,1) * t596 + mrSges(2,2) * t594 - m(4) * t421 - t275 * mrSges(4,1) - mrSges(4,3) * t422 - m(3) * t536 - t353 * mrSges(3,1) + mrSges(3,2) * t436 - mrSges(3,3) * t512 - m(6) * (t225 * pkin(4) + t417) - t171 * mrSges(6,1) - t170 * mrSges(6,2) - m(5) * t417 - t225 * mrSges(5,1) + (-t700 + t712) * t274 + t641 * t224) * g(2) + t44 * t73 + t76 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t127 * t626 + t49 * t630 + t81 * t13 + t4 * t98 + t5 * t99 + t100 * t45 + (-t41 * mrSges(5,3) - Ifges(5,2) * t617 - Ifges(5,6) * t611 - Ifges(5,4) * t618 + t697 / 0.2e1 - t79 / 0.2e1 + Ifges(7,6) * t636 + Ifges(7,5) * t637 + Ifges(6,3) * t620 + Ifges(7,3) * t621 + Ifges(6,6) * t627 + Ifges(6,5) * t628 + t703 + t704) * t268; (t323 * t9 + t358 * t37 - t442 * t8 + t55 * t668 - t56 * t714 + t665 * t88 + t651) * m(6) - t714 * t146 - t442 * t60 + t665 * (t73 - t671) - t237 * t481 + (-t87 + t210) * t514 + t211 * t550 + t356 * t505 - t357 * t504 + t455 + t407 * t169 + t359 * t117 + (t13 + t681) * t358 + t323 * t59 + t708 * t724 + t672 * t98 + t673 * t99 + (t16 * t673 + t17 * t672 + t2 * t253 + t24 * t358 + t252 * t3 + t665 * t76 + t651) * m(7) + t666 * t206 + (-t42 * t358 + t41 * t359 + (-t108 * t595 + t177 * t533) * t404 - t177 * t332 + t666 * t96 - t665 * t95 + t651) * m(5) + t668 * t147 + (t200 * t332 - t201 * t333 - t245 * t481 + t228 * t407 + (t595 * t113 + t112 * t412 + (-t200 * t412 + t201 * t595) * qJD(3)) * t404 + t651) * m(4) + (-t459 * t535 + t388 + t651) * m(3) + t252 * t25 + t253 * t26 + t722 * t258; t3 * (-mrSges(7,1) * t415 + mrSges(7,3) * t348) + t2 * (mrSges(7,2) * t415 - mrSges(7,3) * t347) + t24 * (mrSges(7,1) * t347 - mrSges(7,2) * t348) + (-Ifges(7,4) * t348 - Ifges(7,2) * t347 - Ifges(7,6) * t415) * t636 + (-Ifges(7,1) * t348 - Ifges(7,4) * t347 - Ifges(7,5) * t415) * t637 + (-Ifges(7,5) * t348 - Ifges(7,6) * t347 - Ifges(7,3) * t415) * t621 + (t254 * t466 + t255 * t469 + t296 * t463) * qJD(4) / 0.2e1 - (Ifges(4,1) * t301 + t154 - t566) * t304 / 0.2e1 - (-Ifges(4,2) * t304 + t215 + t295) * t301 / 0.2e1 + t296 * t177 * (mrSges(5,1) * t411 + mrSges(5,2) * t415) + t304 * t698 + (t274 * t705 + t275 * t643) * g(1) + (t410 * t493 - t230 / 0.2e1) * t106 + t155 * t556 / 0.2e1 + t706 * t136 + t494 * t542 - t39 * t545 / 0.2e1 + t8 * (mrSges(6,2) * t415 - mrSges(6,3) * t545) + t517 + t214 * t598 + (Ifges(5,3) * t304 + t301 * t463) * t601 + (Ifges(5,5) * t304 + t301 * t469) * t604 + (Ifges(5,6) * t304 + t301 * t466) * t606 + (-t461 * t529 + (Ifges(6,3) * t411 + t415 * t462) * qJD(4)) * t607 + (Ifges(6,5) * t230 + Ifges(6,6) * t229 + Ifges(6,3) * t556) * t608 + t108 * t473 + t470 * t565 + t301 * t568 + (-t120 * t88 + t342 * t9 + t343 * t8 + t683 * t55 + t684 * t56) * m(6) + (-pkin(3) * t108 - t141 * t95 - t142 * t96) * m(5) + ((t531 * t88 + t565) * m(6) + ((-t411 * t96 - t415 * t95) * qJD(4) + t657) * m(5) + t415 * t117 + t681 * t411) * pkin(10) + (t414 * t493 - t229 / 0.2e1) * t105 + (-t555 / 0.2e1 + t494) * t156 - t644 + t500 * t531 + (-t522 - t142) * t206 + t9 * (-mrSges(6,1) * t415 - mrSges(6,3) * t543) + (t270 * t705 - t273 * t643) * g(2) + (-g(1) * t275 + g(2) * t273 - g(3) * t313 - t658 * t96 + (-t531 + t555) * t95 + t657) * mrSges(5,3) + (-m(5) * t492 - t309 * t701 + t313 * t663 - t425 * t652 + t475) * g(3) + t415 * t79 / 0.2e1 - t697 * t415 / 0.2e1 - t304 * t699 + t682 * t73 + t683 * t147 + t684 * t146 + t687 * t99 + t688 * t98 + (t16 * t687 + t17 * t688 + t2 * t249 + t24 * t375 + t248 * t3 + t682 * t76) * m(7) + t375 * t13 - t341 * (Ifges(4,5) * t301 - Ifges(4,6) * t304) / 0.2e1 + t342 * t59 + t343 * t60 + (mrSges(7,1) * t658 + mrSges(7,3) * t674) * t16 + (mrSges(7,1) * t675 - mrSges(7,2) * t674) * t76 + (-mrSges(7,2) * t658 - mrSges(7,3) * t675) * t17 + (mrSges(6,1) * t658 + mrSges(6,3) * t669) * t55 + (mrSges(6,1) * t670 - mrSges(6,2) * t669) * t88 + (-mrSges(6,2) * t658 - mrSges(6,3) * t670) * t56 - t245 * (mrSges(4,1) * t304 + mrSges(4,2) * t301) - t200 * t258 + (t580 + t646) * t201 + t248 * t25 + t249 * t26 + t160 * t631 + t261 * t632 + t159 * t633 + t543 * t634 - t348 * t638 - t347 * t639 + (-t521 - t141) * t207 + (Ifges(7,5) * t260 + Ifges(7,6) * t261 + Ifges(7,3) * t532) * t609 + (Ifges(7,5) * t160 + Ifges(7,6) * t159 + Ifges(7,3) * t556) * t610 + (Ifges(5,5) * t411 + Ifges(5,6) * t415) * t611 + (-t467 * t529 + (Ifges(6,5) * t411 + t415 * t468) * qJD(4)) * t613 + (Ifges(6,1) * t230 + Ifges(6,4) * t229 + Ifges(6,5) * t556) * t614 + (-t464 * t529 + (Ifges(6,6) * t411 + t415 * t465) * qJD(4)) * t615 + (Ifges(6,4) * t230 + Ifges(6,2) * t229 + Ifges(6,6) * t556) * t616 + t689 * (t532 / 0.2e1 - t556 / 0.2e1) + (Ifges(5,2) * t415 + t575) * t617 + (Ifges(5,1) * t411 + t574) * t618 + t532 * t619 + (-Ifges(6,3) * t415 + t411 * t462) * t620 + (Ifges(7,1) * t260 + Ifges(7,4) * t261 + Ifges(7,5) * t532) * t622 + (Ifges(7,1) * t160 + Ifges(7,4) * t159 + Ifges(7,5) * t556) * t623 + (Ifges(7,4) * t260 + Ifges(7,2) * t261 + Ifges(7,6) * t532) * t624 + (Ifges(7,4) * t160 + Ifges(7,2) * t159 + Ifges(7,6) * t556) * t625 + (-Ifges(6,6) * t415 + t411 * t465) * t627 + (-Ifges(6,5) * t415 + t411 * t468) * t628 + t411 * t629 + t260 * t630 - pkin(3) * t87; (-Ifges(7,4) * t180 - Ifges(7,2) * t179) * t625 - t76 * (mrSges(7,1) * t179 - mrSges(7,2) * t180) + (-Ifges(7,1) * t180 - Ifges(7,4) * t179) * t623 + (-Ifges(7,5) * t180 - Ifges(7,6) * t179) * t610 + (t204 * t465 + t205 * t468 + t251 * t462) * qJD(5) / 0.2e1 + (-t2 * mrSges(7,3) + t24 * mrSges(7,1) - 0.2e1 * t639 - (Ifges(7,1) * t622 + Ifges(7,4) * t624 + Ifges(7,5) * t609 + t630 + t649) * t655) * t458 + (mrSges(7,2) * t24 - mrSges(7,3) * t3 + 0.2e1 * t638) * t368 + (-t641 * t720 - t654 * t719) * g(2) + (-t16 * t180 + t17 * t179) * mrSges(7,3) + t518 + t155 * t603 + t546 * t605 + t37 * t471 + (-Ifges(5,2) * t255 + t156 + t250 + t542) * t606 + (t500 + t451) * qJD(5) + (Ifges(6,5) * t614 + Ifges(7,5) * t623 - Ifges(5,6) * t601 + Ifges(6,6) * t616 + Ifges(7,6) * t625 + Ifges(6,3) * t608 + Ifges(7,3) * t610 - t702) * t255 + t645 - pkin(4) * t45 + (t579 - t206) * t95 + (-pkin(4) * t37 - t55 * t77 - t56 * t78) * m(6) - t398 * t13 + t685 * t99 + t686 * t98 + (t16 * t685 + t17 * t686 + t2 * t331 - t24 * t398 + t3 * t330 + t680 * t76) * m(7) + (Ifges(5,1) * t254 - t567 + t689) * t604 + t680 * t73 + t330 * t25 + t331 * t26 + (-m(6) * t88 + t578 + t671) * t96 + (t268 * t662 + t269 * t660 + t474) * g(3) + (-t146 * t530 - t147 * t528 + t414 * t60 + m(6) * ((-t410 * t56 - t414 * t55) * qJD(5) + t656) - t410 * t59) * pkin(11) + (t224 * t654 + t225 * t641) * g(1) - (Ifges(7,4) * t622 + Ifges(7,2) * t624 + Ifges(7,6) * t609 + t632 + t650) * t316 - t180 * t631 - t179 * t633 + t410 * t634 + t414 * t635 - t78 * t146 - t77 * t147 + (-t713 * t56 + (-t528 + t558) * t55 + t656) * mrSges(6,3) + (Ifges(5,5) * t601 + t462 * t608 + t465 * t616 + t468 * t614 - t451 - t711) * t254 + t461 * t620 + t528 * t626 + t464 * t627 + t467 * t628; (-m(7) * t588 - t472 - t538) * g(3) - t73 * t589 - m(7) * (t16 * t18 + t17 * t19 + t589 * t76) + t457 + t704 + t38 + (-t540 - (-t270 * t410 + t414 * t720) * mrSges(6,2) + t679 * (t270 * t414 + t410 * t720)) * g(2) + (mrSges(6,2) * t171 + t170 * t679 - t539) * g(1) + (-Ifges(6,2) * t205 + t106 + t203) * t616 + (Ifges(7,1) * t623 + Ifges(7,4) * t625 + Ifges(7,5) * t610 + t631 - t649) * t491 - (Ifges(7,4) * t623 + Ifges(7,2) * t625 + Ifges(7,6) * t610 + t633 - t650) * t133 + (t576 + t147) * t56 + (t577 - t146) * t55 - t88 * (mrSges(6,1) * t205 + mrSges(6,2) * t204) + (t2 * t409 + t3 * t413 + (-t16 * t409 + t17 * t413) * qJD(6)) * t640 + (Ifges(6,5) * t204 - Ifges(6,6) * t205) * t608 + t105 * t613 + (Ifges(6,1) * t204 - t573) * t614 + ((-t409 * t99 + t413 * t98) * qJD(6) + t25 * t413 + t26 * t409) * pkin(5) - t19 * t98 - t18 * t99; -t76 * (mrSges(7,1) * t133 + mrSges(7,2) * t491) + t17 * t99 + (Ifges(7,1) * t491 - t570) * t623 + t68 * t622 + (Ifges(7,5) * t491 - Ifges(7,6) * t133) * t610 - t16 * t98 - g(1) * t539 - g(2) * t540 - g(3) * t538 + (t133 * t17 + t16 * t491) * mrSges(7,3) + t457 + (-Ifges(7,2) * t133 + t129 + t69) * t625;];
tau  = t1;
