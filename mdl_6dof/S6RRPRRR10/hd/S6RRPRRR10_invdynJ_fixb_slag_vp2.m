% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR10
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:22
% EndTime: 2019-03-09 14:19:55
% DurationCPUTime: 56.94s
% Computational Cost: add. (35525->1070), mult. (85481->1456), div. (0->0), fcn. (71220->18), ass. (0->463)
t418 = sin(pkin(12));
t420 = cos(pkin(12));
t425 = sin(qJ(4));
t430 = cos(qJ(4));
t367 = t418 * t430 + t420 * t425;
t419 = sin(pkin(6));
t431 = cos(qJ(2));
t530 = t419 * t431;
t437 = t367 * t530;
t303 = qJD(1) * t437;
t357 = t367 * qJD(4);
t711 = -t303 + t357;
t429 = cos(qJ(5));
t410 = pkin(5) * t429 + pkin(4);
t417 = qJ(5) + qJ(6);
t413 = sin(t417);
t414 = cos(t417);
t424 = sin(qJ(5));
t469 = -mrSges(6,1) * t429 + mrSges(6,2) * t424;
t647 = m(6) * pkin(4) + m(7) * t410 + mrSges(7,1) * t414 - mrSges(7,2) * t413 + mrSges(5,1) - t469;
t426 = sin(qJ(2));
t461 = pkin(2) * t426 - qJ(3) * t431;
t519 = qJD(1) * t419;
t347 = t461 * t519;
t492 = t426 * t519;
t421 = cos(pkin(6));
t518 = qJD(1) * t421;
t507 = pkin(1) * t518;
t348 = -pkin(8) * t492 + t431 * t507;
t258 = t420 * t347 - t348 * t418;
t528 = t420 * t431;
t439 = (pkin(3) * t426 - pkin(9) * t528) * t419;
t212 = qJD(1) * t439 + t258;
t259 = t418 * t347 + t420 * t348;
t491 = t431 * t519;
t477 = t418 * t491;
t236 = -pkin(9) * t477 + t259;
t366 = t418 * t425 - t430 * t420;
t571 = pkin(9) + qJ(3);
t378 = t571 * t418;
t379 = t571 * t420;
t656 = -t430 * t378 - t379 * t425;
t664 = -t366 * qJD(3) + qJD(4) * t656 - t425 * t212 - t430 * t236;
t432 = -pkin(11) - pkin(10);
t440 = -m(6) * pkin(10) + m(7) * t432 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t710 = pkin(10) * t492 - t664;
t396 = pkin(8) * t491;
t349 = t426 * t507 + t396;
t296 = pkin(3) * t477 + t349;
t436 = t366 * t530;
t304 = qJD(1) * t436;
t356 = t366 * qJD(4);
t709 = -t296 + (-t304 + t356) * pkin(10) + t711 * pkin(4);
t260 = t304 * t424 + t429 * t492;
t510 = qJD(5) * t429;
t447 = -t356 * t424 + t367 * t510;
t661 = t260 + t447;
t385 = qJD(4) - t491;
t585 = t385 / 0.2e1;
t401 = qJD(2) + t518;
t326 = t401 * t420 - t418 * t492;
t327 = t401 * t418 + t420 * t492;
t455 = t326 * t425 + t430 * t327;
t596 = t455 / 0.2e1;
t657 = t430 * t326 - t425 * t327;
t598 = t657 / 0.2e1;
t244 = qJD(5) - t657;
t600 = t244 / 0.2e1;
t240 = qJD(6) + t244;
t602 = t240 / 0.2e1;
t206 = t385 * t424 + t429 * t455;
t604 = t206 / 0.2e1;
t205 = t385 * t429 - t424 * t455;
t606 = t205 / 0.2e1;
t423 = sin(qJ(6));
t428 = cos(qJ(6));
t129 = t205 * t423 + t206 * t428;
t615 = t129 / 0.2e1;
t483 = t428 * t205 - t206 * t423;
t617 = t483 / 0.2e1;
t307 = qJ(3) * t401 + t349;
t335 = (-pkin(2) * t431 - qJ(3) * t426 - pkin(1)) * t419;
t317 = qJD(1) * t335;
t220 = -t307 * t418 + t420 * t317;
t173 = -pkin(3) * t491 - pkin(9) * t327 + t220;
t221 = t420 * t307 + t418 * t317;
t186 = pkin(9) * t326 + t221;
t103 = t173 * t425 + t186 * t430;
t100 = pkin(10) * t385 + t103;
t298 = -pkin(2) * t401 + qJD(3) - t348;
t245 = -pkin(3) * t326 + t298;
t114 = -pkin(4) * t657 - pkin(10) * t455 + t245;
t54 = -t100 * t424 + t429 * t114;
t48 = -pkin(11) * t206 + t54;
t45 = pkin(5) * t244 + t48;
t55 = t100 * t429 + t114 * t424;
t49 = pkin(11) * t205 + t55;
t546 = t423 * t49;
t15 = t428 * t45 - t546;
t545 = t428 * t49;
t16 = t423 * t45 + t545;
t689 = t245 * mrSges(5,1) + t54 * mrSges(6,1) + t15 * mrSges(7,1) - t55 * mrSges(6,2) - t16 * mrSges(7,2);
t563 = Ifges(5,4) * t455;
t158 = Ifges(5,2) * t657 + t385 * Ifges(5,6) + t563;
t678 = t206 * Ifges(6,5) + t129 * Ifges(7,5) + t205 * Ifges(6,6) + Ifges(7,6) * t483 + t244 * Ifges(6,3) + t240 * Ifges(7,3);
t707 = mrSges(5,3) * t103 - t678 / 0.2e1 + t158 / 0.2e1;
t708 = -Ifges(5,4) * t596 + Ifges(6,5) * t604 + Ifges(7,5) * t615 - Ifges(5,2) * t598 - Ifges(5,6) * t585 + Ifges(6,6) * t606 + Ifges(7,6) * t617 + Ifges(6,3) * t600 + Ifges(7,3) * t602 + t689 - t707;
t416 = pkin(12) + qJ(4);
t411 = sin(t416);
t412 = cos(t416);
t472 = -mrSges(4,1) * t420 + mrSges(4,2) * t418;
t442 = m(4) * pkin(2) - t472;
t706 = t440 * t411 - t647 * t412 - t442;
t705 = t424 * t710 + t429 * t709;
t409 = pkin(3) * t420 + pkin(2);
t278 = pkin(4) * t366 - pkin(10) * t367 - t409;
t302 = -t378 * t425 + t379 * t430;
t511 = qJD(5) * t424;
t670 = t278 * t510 - t302 * t511 + t424 * t709 - t429 * t710;
t243 = Ifges(5,4) * t657;
t159 = Ifges(5,1) * t455 + t385 * Ifges(5,5) + t243;
t609 = t159 / 0.2e1;
t102 = t173 * t430 - t425 * t186;
t643 = t245 * mrSges(5,2) - t102 * mrSges(5,3);
t704 = Ifges(5,1) * t596 + Ifges(5,4) * t598 + Ifges(5,5) * t585 + t609 + t643;
t532 = t419 * t426;
t402 = pkin(8) * t532;
t578 = pkin(1) * t421;
t506 = qJD(2) * t578;
t479 = qJD(1) * t506;
t509 = qJDD(1) * t421;
t502 = pkin(1) * t509;
t274 = -qJD(2) * t396 - qJDD(1) * t402 - t426 * t479 + t431 * t502;
t400 = qJDD(2) + t509;
t251 = -pkin(2) * t400 + qJDD(3) - t274;
t515 = qJD(2) * t431;
t353 = (qJD(1) * t515 + qJDD(1) * t426) * t419;
t282 = -t353 * t418 + t400 * t420;
t185 = -pkin(3) * t282 + t251;
t516 = qJD(2) * t419;
t490 = t426 * t516;
t684 = -qJD(1) * t490 + qJDD(1) * t530;
t273 = pkin(8) * t684 + t426 * t502 + t431 * t479;
t231 = qJ(3) * t400 + qJD(3) * t401 + t273;
t514 = qJD(3) * t426;
t239 = -pkin(2) * t684 - qJ(3) * t353 + (-pkin(1) * qJDD(1) - qJD(1) * t514) * t419;
t155 = -t231 * t418 + t420 * t239;
t283 = t353 * t420 + t400 * t418;
t108 = -pkin(3) * t684 - pkin(9) * t283 + t155;
t156 = t420 * t231 + t418 * t239;
t118 = pkin(9) * t282 + t156;
t512 = qJD(4) * t430;
t513 = qJD(4) * t425;
t44 = t108 * t430 - t425 * t118 - t173 * t513 - t186 * t512;
t340 = qJDD(4) - t684;
t591 = t340 / 0.2e1;
t144 = t282 * t430 - t425 * t283 - t326 * t513 - t327 * t512;
t611 = t144 / 0.2e1;
t143 = qJD(4) * t657 + t282 * t425 + t283 * t430;
t612 = t143 / 0.2e1;
t701 = mrSges(5,2) * t185 - mrSges(5,3) * t44 + 0.2e1 * Ifges(5,1) * t612 + 0.2e1 * Ifges(5,4) * t611 + 0.2e1 * Ifges(5,5) * t591;
t142 = qJDD(5) - t144;
t87 = qJD(5) * t205 + t143 * t429 + t340 * t424;
t88 = -qJD(5) * t206 - t143 * t424 + t340 * t429;
t35 = Ifges(6,5) * t87 + Ifges(6,6) * t88 + Ifges(6,3) * t142;
t43 = t425 * t108 + t430 * t118 + t173 * t512 - t186 * t513;
t613 = t142 / 0.2e1;
t138 = qJDD(6) + t142;
t614 = t138 / 0.2e1;
t622 = t88 / 0.2e1;
t623 = t87 / 0.2e1;
t34 = -qJD(6) * t129 - t423 * t87 + t428 * t88;
t631 = t34 / 0.2e1;
t33 = qJD(6) * t483 + t423 * t88 + t428 * t87;
t632 = t33 / 0.2e1;
t41 = pkin(10) * t340 + t43;
t64 = -pkin(4) * t144 - pkin(10) * t143 + t185;
t11 = -t100 * t511 + t114 * t510 + t429 * t41 + t424 * t64;
t12 = -qJD(5) * t55 - t41 * t424 + t429 * t64;
t641 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t6 = pkin(5) * t142 - pkin(11) * t87 + t12;
t7 = pkin(11) * t88 + t11;
t2 = qJD(6) * t15 + t423 * t6 + t428 * t7;
t3 = -qJD(6) * t16 - t423 * t7 + t428 * t6;
t642 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t138;
t700 = t642 + t641 + mrSges(5,1) * t185 + Ifges(6,5) * t623 + Ifges(7,5) * t632 + Ifges(6,6) * t622 + Ifges(7,6) * t631 + Ifges(6,3) * t613 + Ifges(7,3) * t614 + t35 / 0.2e1 + t8 / 0.2e1 - mrSges(5,3) * t43 + (-t340 / 0.2e1 - t591) * Ifges(5,6) + (-t144 / 0.2e1 - t611) * Ifges(5,2) + (-t143 / 0.2e1 - t612) * Ifges(5,4);
t261 = -t304 * t429 + t424 * t492;
t290 = t429 * t302;
t537 = t356 * t429;
t698 = pkin(11) * t261 + pkin(11) * t537 + (-t290 + (pkin(11) * t367 - t278) * t424) * qJD(5) + t705 + t711 * pkin(5);
t697 = pkin(11) * t661 - t670;
t496 = qJD(5) * t432;
t540 = t657 * t424;
t169 = pkin(4) * t455 - pkin(10) * t657;
t73 = t429 * t102 + t424 * t169;
t696 = pkin(11) * t540 + t424 * t496 - t73;
t539 = t657 * t429;
t72 = -t102 * t424 + t429 * t169;
t695 = -pkin(5) * t455 + pkin(11) * t539 + t429 * t496 - t72;
t665 = -qJD(3) * t367 - qJD(4) * t302 - t212 * t430 + t425 * t236;
t692 = t511 - t540;
t634 = m(7) * pkin(5);
t685 = t273 * mrSges(3,2);
t369 = t423 * t429 + t424 * t428;
t160 = t369 * t657;
t453 = t423 * t424 - t428 * t429;
t162 = t453 * t657;
t668 = pkin(4) * t492 - t665;
t683 = -t44 * mrSges(5,1) + t43 * mrSges(5,2);
t467 = t413 * mrSges(7,1) + t414 * mrSges(7,2);
t544 = t429 * mrSges(6,2);
t468 = mrSges(6,1) * t424 + t544;
t508 = t424 * t634;
t676 = m(4) * qJ(3) + mrSges(4,3) + mrSges(5,3);
t681 = -t467 - t508 - t468 - t676;
t196 = t429 * t278 - t302 * t424;
t535 = t367 * t429;
t170 = pkin(5) * t366 - pkin(11) * t535 + t196;
t197 = t424 * t278 + t290;
t536 = t367 * t424;
t180 = -pkin(11) * t536 + t197;
t95 = t170 * t423 + t180 * t428;
t680 = -qJD(6) * t95 + t423 * t697 + t428 * t698;
t94 = t170 * t428 - t180 * t423;
t679 = qJD(6) * t94 + t423 * t698 - t428 * t697;
t387 = t432 * t424;
t388 = t432 * t429;
t313 = t387 * t423 - t388 * t428;
t675 = -qJD(6) * t313 - t423 * t696 + t428 * t695;
t312 = t387 * t428 + t388 * t423;
t674 = qJD(6) * t312 + t423 * t695 + t428 * t696;
t673 = pkin(5) * t692 - t103;
t672 = -t634 - mrSges(6,1);
t671 = -qJD(5) * t197 + t705;
t268 = t453 * t367;
t669 = pkin(5) * t661 + t668;
t364 = pkin(8) * t530 + t426 * t578;
t334 = qJ(3) * t421 + t364;
t253 = -t334 * t418 + t420 * t335;
t355 = t418 * t421 + t420 * t532;
t192 = -pkin(3) * t530 - pkin(9) * t355 + t253;
t254 = t420 * t334 + t418 * t335;
t354 = -t418 * t532 + t420 * t421;
t209 = pkin(9) * t354 + t254;
t125 = t425 * t192 + t430 * t209;
t120 = -pkin(10) * t530 + t125;
t263 = t354 * t425 + t355 * t430;
t577 = pkin(1) * t431;
t337 = t402 + (-pkin(2) - t577) * t421;
t275 = -pkin(3) * t354 + t337;
t454 = t430 * t354 - t355 * t425;
t154 = -pkin(4) * t454 - pkin(10) * t263 + t275;
t76 = t429 * t120 + t424 * t154;
t649 = qJD(5) + qJD(6);
t293 = t649 * t369;
t145 = -t293 * t367 + t356 * t453;
t175 = t260 * t423 + t261 * t428;
t667 = t145 - t175;
t146 = t268 * t649 + t369 * t356;
t174 = t260 * t428 - t261 * t423;
t666 = t146 - t174;
t130 = -mrSges(6,1) * t205 + mrSges(6,2) * t206;
t570 = mrSges(5,3) * t455;
t211 = mrSges(5,1) * t385 - t570;
t663 = t211 - t130;
t481 = mrSges(3,3) * t492;
t662 = -mrSges(3,1) * t401 - mrSges(4,1) * t326 + mrSges(4,2) * t327 + t481;
t446 = t367 * t511 + t537;
t660 = t261 + t446;
t292 = t649 * t453;
t659 = -t292 + t162;
t658 = -t293 + t160;
t652 = -t155 * t418 + t156 * t420;
t651 = t11 * t429 - t12 * t424;
t650 = m(7) + m(6) + m(5);
t567 = Ifges(3,4) * t426;
t648 = pkin(1) * (mrSges(3,1) * t426 + mrSges(3,2) * t431) - t426 * (Ifges(3,1) * t431 - t567) / 0.2e1;
t99 = -pkin(4) * t385 - t102;
t80 = -pkin(5) * t205 + t99;
t645 = -mrSges(7,1) * t80 + mrSges(7,3) * t16;
t644 = mrSges(7,2) * t80 - mrSges(7,3) * t15;
t640 = mrSges(3,2) + t681;
t639 = mrSges(3,1) - t706;
t586 = -t385 / 0.2e1;
t599 = -t657 / 0.2e1;
t601 = -t244 / 0.2e1;
t603 = -t240 / 0.2e1;
t605 = -t206 / 0.2e1;
t607 = -t205 / 0.2e1;
t616 = -t129 / 0.2e1;
t618 = -t483 / 0.2e1;
t637 = Ifges(6,5) * t605 + Ifges(7,5) * t616 - Ifges(5,2) * t599 - Ifges(5,6) * t586 + Ifges(6,6) * t607 + Ifges(7,6) * t618 + Ifges(6,3) * t601 + Ifges(7,3) * t603 - t689;
t636 = t419 ^ 2;
t635 = Ifges(7,4) * t632 + Ifges(7,2) * t631 + Ifges(7,6) * t614;
t633 = Ifges(7,1) * t632 + Ifges(7,4) * t631 + Ifges(7,5) * t614;
t630 = Ifges(6,1) * t623 + Ifges(6,4) * t622 + Ifges(6,5) * t613;
t559 = Ifges(7,4) * t129;
t60 = Ifges(7,2) * t483 + Ifges(7,6) * t240 + t559;
t629 = -t60 / 0.2e1;
t628 = t60 / 0.2e1;
t123 = Ifges(7,4) * t483;
t61 = Ifges(7,1) * t129 + Ifges(7,5) * t240 + t123;
t627 = -t61 / 0.2e1;
t626 = t61 / 0.2e1;
t562 = Ifges(6,4) * t206;
t97 = t205 * Ifges(6,2) + t244 * Ifges(6,6) + t562;
t620 = t97 / 0.2e1;
t200 = Ifges(6,4) * t205;
t98 = t206 * Ifges(6,1) + t244 * Ifges(6,5) + t200;
t619 = -t98 / 0.2e1;
t597 = -t455 / 0.2e1;
t594 = t282 / 0.2e1;
t593 = t283 / 0.2e1;
t590 = t354 / 0.2e1;
t589 = t355 / 0.2e1;
t584 = t421 / 0.2e1;
t583 = t429 / 0.2e1;
t582 = t431 / 0.2e1;
t581 = cos(qJ(1));
t576 = pkin(5) * t206;
t569 = mrSges(6,3) * t205;
t568 = mrSges(6,3) * t206;
t566 = Ifges(3,4) * t431;
t565 = Ifges(4,4) * t418;
t564 = Ifges(4,4) * t420;
t561 = Ifges(6,4) * t424;
t560 = Ifges(6,4) * t429;
t558 = Ifges(4,3) * t426;
t553 = t657 * Ifges(5,6);
t552 = t455 * Ifges(5,5);
t551 = t326 * Ifges(4,6);
t550 = t327 * Ifges(4,5);
t549 = t385 * Ifges(5,3);
t548 = t401 * Ifges(3,5);
t547 = t401 * Ifges(3,6);
t533 = t418 * t431;
t427 = sin(qJ(1));
t531 = t419 * t427;
t529 = t420 * (Ifges(4,1) * t327 + Ifges(4,4) * t326 - Ifges(4,5) * t491);
t527 = t426 * t427;
t526 = t427 * t431;
t494 = t581 * t426;
t359 = t421 * t494 + t526;
t495 = t419 * t581;
t285 = t359 * t412 - t411 * t495;
t493 = t581 * t431;
t358 = -t421 * t493 + t527;
t525 = (-t285 * t413 + t358 * t414) * mrSges(7,1) + (-t285 * t414 - t358 * t413) * mrSges(7,2);
t361 = -t421 * t527 + t493;
t289 = t361 * t412 + t411 * t531;
t360 = t421 * t526 + t494;
t218 = -t289 * t413 + t360 * t414;
t219 = t289 * t414 + t360 * t413;
t524 = t218 * mrSges(7,1) - t219 * mrSges(7,2);
t333 = t411 * t421 + t412 * t532;
t523 = (-t333 * t413 - t414 * t530) * mrSges(7,1) + (-t333 * t414 + t413 * t530) * mrSges(7,2);
t311 = (qJD(2) * t461 - t514) * t419;
t350 = -pkin(8) * t490 + t431 * t506;
t321 = qJD(3) * t421 + t350;
t233 = t418 * t311 + t420 * t321;
t520 = t581 * pkin(1) + pkin(8) * t531;
t501 = t424 * t530;
t500 = t429 * t530;
t499 = t98 * t583;
t498 = Ifges(5,5) * t143 + Ifges(5,6) * t144 + Ifges(5,3) * t340;
t497 = Ifges(3,5) * t353 + Ifges(3,6) * t684 + Ifges(3,3) * t400;
t485 = -t511 / 0.2e1;
t484 = -t427 * pkin(1) + pkin(8) * t495;
t199 = -t282 * mrSges(4,1) + t283 * mrSges(4,2);
t81 = -t144 * mrSges(5,1) + t143 * mrSges(5,2);
t75 = -t120 * t424 + t429 * t154;
t124 = t192 * t430 - t425 * t209;
t232 = t420 * t311 - t321 * t418;
t284 = -t359 * t411 - t412 * t495;
t480 = mrSges(3,3) * t491;
t478 = t418 * t495;
t476 = t418 * t419 * t515;
t119 = pkin(4) * t530 - t124;
t471 = mrSges(4,1) * t418 + mrSges(4,2) * t420;
t466 = Ifges(4,1) * t420 - t565;
t465 = Ifges(6,1) * t429 - t561;
t464 = -Ifges(4,2) * t418 + t564;
t463 = -Ifges(6,2) * t424 + t560;
t462 = Ifges(6,5) * t429 - Ifges(6,6) * t424;
t450 = -t263 * t429 + t501;
t53 = -pkin(5) * t454 + pkin(11) * t450 + t75;
t241 = -t263 * t424 - t500;
t65 = pkin(11) * t241 + t76;
t27 = -t423 * t65 + t428 * t53;
t28 = t423 * t53 + t428 * t65;
t459 = -mrSges(3,2) + t676;
t458 = t418 * pkin(3) * t531 + t360 * t571 + t361 * t409 + t520;
t150 = -mrSges(6,2) * t244 + t569;
t151 = mrSges(6,1) * t244 - t568;
t457 = t150 * t429 - t151 * t424;
t164 = t241 * t428 + t423 * t450;
t165 = t241 * t423 - t428 * t450;
t234 = -t289 * t424 + t360 * t429;
t190 = qJD(2) * t439 + t232;
t207 = -pkin(9) * t476 + t233;
t71 = t190 * t430 - t192 * t513 - t425 * t207 - t209 * t512;
t452 = t642 + t8;
t448 = t99 * t468;
t444 = t298 * t471;
t201 = -qJD(2) * t436 + qJD(4) * t454;
t202 = qJD(2) * t437 + qJD(4) * t263;
t351 = t364 * qJD(2);
t297 = pkin(3) * t476 + t351;
t109 = pkin(4) * t202 - pkin(10) * t201 + t297;
t70 = t425 * t190 + t192 * t512 + t430 * t207 - t209 * t513;
t67 = pkin(10) * t490 + t70;
t20 = t424 * t109 - t120 * t511 + t154 * t510 + t429 * t67;
t42 = -pkin(4) * t340 - t44;
t68 = -pkin(4) * t490 - t71;
t21 = -qJD(5) * t76 + t429 * t109 - t424 * t67;
t394 = Ifges(3,4) * t491;
t363 = t421 * t577 - t402;
t362 = (-mrSges(3,1) * t431 + mrSges(3,2) * t426) * t419;
t346 = -mrSges(3,2) * t401 + t480;
t306 = Ifges(3,1) * t492 + t394 + t548;
t305 = t547 + (t431 * Ifges(3,2) + t567) * t519;
t288 = t361 * t411 - t412 * t531;
t281 = -mrSges(4,1) * t491 - mrSges(4,3) * t327;
t280 = mrSges(4,2) * t491 + mrSges(4,3) * t326;
t267 = t369 * t367;
t257 = pkin(5) * t536 - t656;
t238 = -mrSges(4,1) * t684 - mrSges(4,3) * t283;
t237 = mrSges(4,2) * t684 + mrSges(4,3) * t282;
t235 = t289 * t429 + t360 * t424;
t229 = Ifges(4,4) * t327 + Ifges(4,2) * t326 - Ifges(4,6) * t491;
t228 = -Ifges(4,3) * t491 + t550 + t551;
t210 = -mrSges(5,2) * t385 + mrSges(5,3) * t657;
t179 = t283 * Ifges(4,1) + t282 * Ifges(4,4) - Ifges(4,5) * t684;
t178 = t283 * Ifges(4,4) + t282 * Ifges(4,2) - Ifges(4,6) * t684;
t168 = -mrSges(5,1) * t657 + mrSges(5,2) * t455;
t157 = t549 + t552 + t553;
t132 = qJD(5) * t450 - t201 * t424 + t429 * t490;
t131 = qJD(5) * t241 + t201 * t429 + t424 * t490;
t122 = -mrSges(5,2) * t340 + mrSges(5,3) * t144;
t121 = mrSges(5,1) * t340 - mrSges(5,3) * t143;
t92 = mrSges(7,1) * t240 - mrSges(7,3) * t129;
t91 = -mrSges(7,2) * t240 + mrSges(7,3) * t483;
t90 = -pkin(5) * t241 + t119;
t74 = -mrSges(7,1) * t483 + mrSges(7,2) * t129;
t57 = -mrSges(6,2) * t142 + mrSges(6,3) * t88;
t56 = mrSges(6,1) * t142 - mrSges(6,3) * t87;
t52 = -qJD(6) * t165 - t131 * t423 + t132 * t428;
t51 = qJD(6) * t164 + t131 * t428 + t132 * t423;
t47 = -pkin(5) * t132 + t68;
t46 = -mrSges(6,1) * t88 + mrSges(6,2) * t87;
t36 = t87 * Ifges(6,4) + t88 * Ifges(6,2) + t142 * Ifges(6,6);
t26 = -mrSges(7,2) * t138 + mrSges(7,3) * t34;
t25 = mrSges(7,1) * t138 - mrSges(7,3) * t33;
t22 = -pkin(5) * t88 + t42;
t19 = t428 * t48 - t546;
t18 = -t423 * t48 - t545;
t17 = pkin(11) * t132 + t20;
t14 = pkin(5) * t202 - pkin(11) * t131 + t21;
t13 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t5 = -qJD(6) * t28 + t14 * t428 - t17 * t423;
t4 = qJD(6) * t27 + t14 * t423 + t17 * t428;
t1 = [t662 * t351 + (Ifges(7,5) * t165 + Ifges(7,6) * t164) * t614 + (Ifges(7,5) * t51 + Ifges(7,6) * t52) * t602 + t156 * mrSges(4,3) * t354 + t42 * (-mrSges(6,1) * t241 - mrSges(6,2) * t450) - t450 * t630 + (Ifges(6,4) * t131 + Ifges(6,2) * t132) * t606 + ((t327 * t466 / 0.2e1 + t326 * t464 / 0.2e1 + t548 / 0.2e1 + t444 + t306 / 0.2e1 - t348 * mrSges(3,3) - t418 * t229 / 0.2e1 + t529 / 0.2e1 + (-t220 * t420 - t221 * t418) * mrSges(4,3)) * t431 + (t220 * mrSges(4,1) + t550 / 0.2e1 + t551 / 0.2e1 - t221 * mrSges(4,2) - t547 / 0.2e1 + t553 / 0.2e1 + t552 / 0.2e1 + t549 / 0.2e1 + t102 * mrSges(5,1) - t103 * mrSges(5,2) + t157 / 0.2e1 + t228 / 0.2e1 - t305 / 0.2e1 - t349 * mrSges(3,3)) * t426 + (-t431 * (Ifges(4,5) * t528 - Ifges(4,6) * t533 + t558) / 0.2e1 + (-Ifges(3,2) * t426 + t566) * t582 - t648) * t519) * t516 + (-Ifges(5,6) * t611 - Ifges(5,5) * t612 - Ifges(5,3) * t591 - Ifges(4,6) * t282 - t283 * Ifges(4,5) - t498 / 0.2e1 - t155 * mrSges(4,1) + t273 * mrSges(3,3) + t156 * mrSges(4,2) + t683) * t530 + (Ifges(7,4) * t165 + Ifges(7,2) * t164) * t631 + t708 * t202 + (Ifges(6,5) * t131 + Ifges(6,6) * t132) * t600 + (-Ifges(6,5) * t450 + Ifges(6,6) * t241) * t613 - (-Ifges(3,6) * t421 / 0.2e1 + Ifges(4,5) * t589 + Ifges(4,6) * t590 - t364 * mrSges(3,3) + (-pkin(1) * mrSges(3,1) - t567 + (-Ifges(3,2) - Ifges(4,3)) * t431) * t419) * t684 + (-Ifges(6,4) * t450 + Ifges(6,2) * t241) * t622 + (Ifges(6,1) * t131 + Ifges(6,4) * t132) * t604 - t700 * t454 + (Ifges(7,4) * t51 + Ifges(7,2) * t52) * t617 + (t11 * t241 + t12 * t450 - t131 * t54 + t132 * t55) * mrSges(6,3) + t701 * t263 + t704 * t201 + (Ifges(7,1) * t165 + Ifges(7,4) * t164) * t632 + (Ifges(7,1) * t51 + Ifges(7,4) * t52) * t615 + (-Ifges(6,1) * t450 + Ifges(6,4) * t241) * t623 + (-m(4) * (-pkin(2) * t359 + t484) - (-t359 * t420 + t478) * mrSges(4,1) - (t359 * t418 + t420 * t495) * mrSges(4,2) - m(3) * t484 + t359 * mrSges(3,1) - mrSges(3,3) * t495 + t427 * mrSges(2,1) + t581 * mrSges(2,2) + t647 * t285 + (-t424 * t672 + t459 + t467 + t544) * t358 + t440 * t284 + t650 * (-pkin(3) * t478 + t358 * t571 + t359 * t409 - t484)) * g(1) + (-t581 * mrSges(2,1) - m(7) * (t289 * t410 + t458) - t219 * mrSges(7,1) - t218 * mrSges(7,2) - m(5) * t458 - t289 * mrSges(5,1) - m(6) * (pkin(4) * t289 + t458) - t235 * mrSges(6,1) - t234 * mrSges(6,2) + (-mrSges(3,1) - t442) * t361 + (mrSges(2,2) + (-mrSges(3,3) - t471) * t419) * t427 + (-t459 - t508) * t360 + t440 * t288 + (-m(4) - m(3)) * t520) * g(2) + (-t15 * t51 + t16 * t52 + t164 * t2 - t165 * t3) * mrSges(7,3) + m(4) * (t155 * t253 + t156 * t254 + t220 * t232 + t221 * t233 + t251 * t337 + t298 * t351) + m(5) * (t102 * t71 + t103 * t70 + t124 * t44 + t125 * t43 + t185 * t275 + t245 * t297) + m(6) * (t11 * t76 + t119 * t42 + t12 * t75 + t20 * t55 + t21 * t54 + t68 * t99) + m(7) * (t15 * t5 + t16 * t4 + t2 * t28 + t22 * t90 + t27 * t3 + t47 * t80) + (-pkin(1) * t362 * t419 + Ifges(2,3)) * qJDD(1) - t155 * mrSges(4,3) * t355 + t28 * t26 + t27 * t25 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t636 + t273 * t364 + t274 * t363 - t348 * t351 + t349 * t350) + (Ifges(3,3) * t584 - t364 * mrSges(3,2) + t363 * mrSges(3,1) + (Ifges(3,5) * t426 + Ifges(3,6) * t431) * t419) * t400 + (Ifges(3,5) * t584 - t363 * mrSges(3,3) + (-pkin(1) * mrSges(3,2) + t426 * Ifges(3,1) + t566) * t419) * t353 + t251 * (-mrSges(4,1) * t354 + mrSges(4,2) * t355) + t350 * t346 + t337 * t199 + t297 * t168 + t233 * t280 + t232 * t281 + t275 * t81 + t164 * t635 + t165 * t633 + t51 * t626 + t52 * t628 + t132 * t620 + t497 * t584 + t179 * t589 + t178 * t590 + (Ifges(4,1) * t355 + Ifges(4,4) * t354) * t593 + (Ifges(4,4) * t355 + Ifges(4,2) * t354) * t594 - t421 * t685 + t47 * t74 + t75 * t56 + t76 * t57 + t80 * (-mrSges(7,1) * t52 + mrSges(7,2) * t51) + t90 * t13 + t4 * t91 + t5 * t92 + t119 * t46 + t124 * t121 + t125 * t122 + t68 * t130 + t131 * t98 / 0.2e1 + t99 * (-mrSges(6,1) * t132 + mrSges(6,2) * t131) + t20 * t150 + t21 * t151 + t22 * (-mrSges(7,1) * t164 + mrSges(7,2) * t165) + t274 * (mrSges(3,1) * t421 - mrSges(3,3) * t532) + t70 * t210 + t71 * t211 + t241 * t36 / 0.2e1 + t253 * t238 + t254 * t237; (mrSges(6,1) * t661 - mrSges(6,2) * t660) * t99 - t661 * t97 / 0.2e1 + (-t11 * t536 - t12 * t535 + t54 * t660 - t55 * t661) * mrSges(6,3) + (-m(4) * t298 + t481 - t662) * t349 + t664 * t210 + t665 * t211 + (-mrSges(7,1) * t666 + mrSges(7,2) * t667) * t80 + (-t15 * t667 + t16 * t666 - t2 * t267 + t268 * t3) * mrSges(7,3) + t668 * t130 + (-Ifges(6,4) * t446 - Ifges(6,2) * t447) * t606 + (Ifges(6,4) * t261 + Ifges(6,2) * t260) * t607 - ((-Ifges(3,2) * t492 + t306 + t394 + t529) * t431 + (t228 + t157) * t426 + t327 * (Ifges(4,5) * t426 + t431 * t466) + t326 * (Ifges(4,6) * t426 + t431 * t464) + t401 * (Ifges(3,5) * t431 - Ifges(3,6) * t426)) * t519 / 0.2e1 + (-Ifges(5,4) * t597 + t637 + t707) * t303 + (-Ifges(6,1) * t446 - Ifges(6,4) * t447) * t604 + (Ifges(6,1) * t261 + Ifges(6,4) * t260) * t605 + (t237 * t420 - t238 * t418) * qJ(3) + (t103 * t492 + t245 * t304) * mrSges(5,2) + (-Ifges(5,5) * t304 + Ifges(5,3) * t492) * t586 + (-Ifges(5,4) * t304 + Ifges(5,6) * t492) * t599 - t102 * (mrSges(5,1) * t492 + mrSges(5,3) * t304) + (-Ifges(5,1) * t304 + Ifges(5,5) * t492) * t597 + (-t346 + t480) * t348 + (Ifges(7,4) * t145 + Ifges(7,2) * t146) * t617 + (Ifges(7,4) * t175 + Ifges(7,2) * t174) * t618 + (-t220 * t258 - t221 * t259 - pkin(2) * t251 + (-t220 * t418 + t221 * t420) * qJD(3) + t652 * qJ(3)) * m(4) + t652 * mrSges(4,3) + (Ifges(7,5) * t145 + Ifges(7,6) * t146) * t602 + (Ifges(7,5) * t175 + Ifges(7,6) * t174) * t603 + (-Ifges(6,5) * t446 - Ifges(6,6) * t447) * t600 + (Ifges(6,5) * t261 + Ifges(6,6) * t260) * t601 + (t362 - t650 * t409 * t530 + (t706 * t431 + (-t571 * t650 + t681) * t426) * t419) * g(3) + t708 * t357 - t684 * (Ifges(4,5) * t418 + Ifges(4,6) * t420) / 0.2e1 + ((t558 + (Ifges(4,5) * t420 - Ifges(4,6) * t418) * t431) * t582 + t648) * qJD(1) ^ 2 * t636 + t700 * t366 + (t42 * t468 + t462 * t613 + t463 * t622 + t465 * t623 + t485 * t98 + t701) * t367 - (t499 + t704) * t356 + (Ifges(7,1) * t145 + Ifges(7,4) * t146) * t615 + (Ifges(7,1) * t175 + Ifges(7,4) * t174) * t616 + (-Ifges(7,4) * t268 - Ifges(7,2) * t267) * t631 + (-Ifges(7,5) * t268 - Ifges(7,6) * t267) * t614 + (-Ifges(7,1) * t268 - Ifges(7,4) * t267) * t632 + t22 * (mrSges(7,1) * t267 - mrSges(7,2) * t268) + t304 * t609 + (-t650 * (-t360 * t409 + t361 * t571) + t640 * t361 + t639 * t360) * g(1) + (-t650 * (-t358 * t409 + t359 * t571) + t640 * t359 + t639 * t358) * g(2) + t229 * t477 / 0.2e1 + (-t221 * (-mrSges(4,2) * t426 - mrSges(4,3) * t533) - t220 * (mrSges(4,1) * t426 - mrSges(4,3) * t528)) * t519 + (t420 * t280 - t418 * t281) * qJD(3) + t305 * t492 / 0.2e1 + t251 * t472 - (t46 - t121) * t656 + (t102 * t665 + t103 * t664 - t185 * t409 - t245 * t296 + t302 * t43 + t44 * t656) * m(5) + (t11 * t197 + t12 * t196 - t42 * t656 + t54 * t671 + t55 * t670 + t668 * t99) * m(6) + t497 + t420 * t178 / 0.2e1 + t418 * t179 / 0.2e1 - t409 * t81 - t296 * t168 + t302 * t122 - t259 * t280 - t258 * t281 + t274 * mrSges(3,1) - t267 * t635 + t535 * t630 - t268 * t633 + t145 * t626 + t175 * t627 + t146 * t628 + t174 * t629 + t261 * t619 + (Ifges(4,1) * t418 + t564) * t593 + (Ifges(4,2) * t420 + t565) * t594 - t685 + t669 * t74 + t670 * t150 + t671 * t151 + t94 * t25 + t95 * t26 + t679 * t91 + t680 * t92 + (t15 * t680 + t16 * t679 + t2 * t95 + t22 * t257 + t3 * t94 + t669 * t80) * m(7) - t444 * t491 - t36 * t536 / 0.2e1 + t196 * t56 + t197 * t57 - pkin(2) * t199 + t257 * t13; t81 + (-t74 + t663) * t455 + (-t210 - t457) * t657 + t457 * qJD(5) + t658 * t92 + t429 * t56 + t424 * t57 - t453 * t25 + t369 * t26 + t327 * t281 - t326 * t280 + t199 + t659 * t91 + (t15 * t658 + t16 * t659 + t2 * t369 - t3 * t453 - t455 * t80) * m(7) + (t11 * t424 + t12 * t429 - t455 * t99 + t244 * (-t424 * t54 + t429 * t55)) * m(6) + (t102 * t455 - t103 * t657 + t185) * m(5) + (t220 * t327 - t221 * t326 + t251) * m(4) + (-g(1) * t360 - g(2) * t358 + g(3) * t530) * (m(4) + t650); -(Ifges(7,1) * t615 + Ifges(7,4) * t617 + Ifges(7,5) * t602 + t626 + t644) * t292 - (Ifges(7,4) * t615 + Ifges(7,2) * t617 + Ifges(7,6) * t602 + t628 + t645) * t293 + (-m(6) * t99 + t570 + t663) * t103 + (Ifges(5,5) * t586 + t462 * t601 + t463 * t607 + t465 * t605 - t448 - t643) * t657 + (-t692 * t55 + (-t510 + t539) * t54 + t651) * mrSges(6,3) + (-t151 * t510 - t150 * t511 - t424 * t56 + t429 * t57 + m(6) * ((-t424 * t55 - t429 * t54) * qJD(5) + t651)) * pkin(10) + (t243 + t159) * t599 + t637 * t455 - t683 - t80 * (mrSges(7,1) * t160 - mrSges(7,2) * t162) + (-Ifges(7,5) * t162 - Ifges(7,6) * t160) * t603 + (-Ifges(7,4) * t162 - Ifges(7,2) * t160) * t618 + (-Ifges(7,1) * t162 - Ifges(7,4) * t160) * t616 + (t448 + t499) * qJD(5) + (-pkin(4) * t42 - t54 * t72 - t55 * t73) * m(6) + (Ifges(5,1) * t657 - t563 + t678) * t597 + (-t15 * t162 + t16 * t160 - t2 * t453 - t3 * t369) * mrSges(7,3) + t22 * (mrSges(7,1) * t453 + mrSges(7,2) * t369) + (Ifges(7,4) * t369 - Ifges(7,2) * t453) * t631 + (Ifges(7,1) * t369 - Ifges(7,4) * t453) * t632 + (Ifges(7,5) * t369 - Ifges(7,6) * t453) * t614 - t453 * t635 + t97 * t485 + (t205 * t463 + t206 * t465 + t244 * t462) * qJD(5) / 0.2e1 + (t288 * t647 + t289 * t440) * g(1) + (-t284 * t647 + t285 * t440) * g(2) + (t440 * t333 - t647 * (-t411 * t532 + t412 * t421)) * g(3) - pkin(4) * t46 + t498 - t410 * t13 + t312 * t25 + t313 * t26 + t369 * t633 + (Ifges(6,2) * t429 + t561) * t622 + (Ifges(6,1) * t424 + t560) * t623 - t162 * t627 - t160 * t629 + t424 * t630 + (Ifges(6,5) * t424 + Ifges(6,6) * t429) * t613 + t539 * t619 + t540 * t620 + t158 * t596 + t36 * t583 + t42 * t469 + t673 * t74 + t674 * t91 + t675 * t92 + (t15 * t675 + t16 * t674 + t2 * t313 - t22 * t410 + t3 * t312 + t673 * t80) * m(7) - t73 * t150 - t72 * t151 - t102 * t210; (mrSges(6,2) * t235 + t234 * t672 - t524) * g(1) + t641 + t35 + (-t523 - (-t333 * t429 + t501) * mrSges(6,2) + t672 * (-t333 * t424 - t500)) * g(3) + (-(-t285 * t429 - t358 * t424) * mrSges(6,2) - t525 + t672 * (-t285 * t424 + t358 * t429)) * g(2) + (-Ifges(6,2) * t206 + t200 + t98) * t607 + t452 - (Ifges(7,4) * t616 + Ifges(7,2) * t618 + Ifges(7,6) * t603 + t629 - t645) * t129 + (Ifges(7,1) * t616 + Ifges(7,4) * t618 + Ifges(7,5) * t603 + t627 - t644) * t483 - t74 * t576 - m(7) * (t15 * t18 + t16 * t19 + t576 * t80) + (t2 * t423 + t3 * t428 + (-t15 * t423 + t16 * t428) * qJD(6)) * t634 + (Ifges(6,5) * t205 - Ifges(6,6) * t206) * t601 + t97 * t604 + (Ifges(6,1) * t205 - t562) * t605 - t19 * t91 - t18 * t92 + (t569 - t150) * t54 - t99 * (mrSges(6,1) * t206 + mrSges(6,2) * t205) + (t568 + t151) * t55 + ((-t423 * t92 + t428 * t91) * qJD(6) + t25 * t428 + t26 * t423) * pkin(5); -t80 * (mrSges(7,1) * t129 + mrSges(7,2) * t483) + (Ifges(7,1) * t483 - t559) * t616 + t60 * t615 + (Ifges(7,5) * t483 - Ifges(7,6) * t129) * t603 - t15 * t91 + t16 * t92 - g(1) * t524 - g(2) * t525 - g(3) * t523 + (t129 * t16 + t15 * t483) * mrSges(7,3) + t452 + (-Ifges(7,2) * t129 + t123 + t61) * t618;];
tau  = t1;
