% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:57
% EndTime: 2019-03-09 17:05:27
% DurationCPUTime: 58.92s
% Computational Cost: add. (22918->1010), mult. (55594->1318), div. (0->0), fcn. (44683->14), ass. (0->436)
t400 = cos(qJ(2));
t396 = sin(qJ(2));
t390 = sin(pkin(6));
t514 = qJD(1) * t390;
t484 = t396 * t514;
t392 = cos(pkin(6));
t513 = qJD(1) * t392;
t501 = pkin(1) * t513;
t313 = -pkin(8) * t484 + t400 * t501;
t417 = (pkin(2) * t396 - pkin(9) * t400) * t390;
t314 = qJD(1) * t417;
t395 = sin(qJ(3));
t399 = cos(qJ(3));
t216 = -t395 * t313 + t399 * t314;
t393 = -qJ(4) - pkin(9);
t466 = qJD(3) * t393;
t520 = t399 * t400;
t730 = -(pkin(3) * t396 - qJ(4) * t520) * t514 - t216 - qJD(4) * t395 + t399 * t466;
t217 = t399 * t313 + t395 * t314;
t483 = t400 * t514;
t454 = t395 * t483;
t729 = -qJ(4) * t454 - qJD(4) * t399 - t395 * t466 + t217;
t389 = sin(pkin(11));
t391 = cos(pkin(11));
t337 = t389 * t399 + t391 * t395;
t253 = t337 * t483;
t324 = t337 * qJD(3);
t718 = -t324 + t253;
t353 = qJD(3) - t483;
t591 = t353 / 0.2e1;
t373 = qJD(2) + t513;
t282 = t373 * t399 - t395 * t484;
t283 = t373 * t395 + t399 * t484;
t465 = t391 * t282 - t283 * t389;
t603 = t465 / 0.2e1;
t728 = -Ifges(5,2) * t603 - Ifges(5,6) * t591;
t672 = Ifges(6,1) + Ifges(7,1);
t670 = Ifges(7,4) + Ifges(6,5);
t652 = t389 * t730 - t729 * t391;
t316 = pkin(8) * t483 + t396 * t501;
t510 = qJD(3) * t395;
t647 = -t316 + (-t454 + t510) * pkin(3);
t260 = -t373 * pkin(2) - t313;
t195 = -t282 * pkin(3) + qJD(4) + t260;
t194 = qJD(5) - t465;
t394 = sin(qJ(5));
t398 = cos(qJ(5));
t261 = pkin(9) * t373 + t316;
t297 = (-pkin(2) * t400 - pkin(9) * t396 - pkin(1)) * t390;
t265 = qJD(1) * t297;
t181 = -t261 * t395 + t399 * t265;
t152 = -qJ(4) * t283 + t181;
t142 = pkin(3) * t353 + t152;
t182 = t261 * t399 + t265 * t395;
t153 = qJ(4) * t282 + t182;
t526 = t391 * t153;
t82 = t389 * t142 + t526;
t78 = pkin(10) * t353 + t82;
t421 = t282 * t389 + t391 * t283;
t92 = -pkin(4) * t465 - pkin(10) * t421 + t195;
t28 = -t394 * t78 + t398 * t92;
t656 = qJD(6) - t28;
t24 = -pkin(5) * t194 + t656;
t29 = t394 * t92 + t398 * t78;
t25 = qJ(6) * t194 + t29;
t661 = Ifges(5,4) * t421;
t727 = t195 * mrSges(5,1) + t28 * mrSges(6,1) - t24 * mrSges(7,1) - t29 * mrSges(6,2) + t25 * mrSges(7,3) - t82 * mrSges(5,3) - t661 / 0.2e1 + t728;
t505 = qJD(1) * qJD(2);
t320 = (qJDD(1) * t396 + t400 * t505) * t390;
t503 = qJDD(1) * t392;
t372 = qJDD(2) + t503;
t183 = qJD(3) * t282 + t320 * t399 + t372 * t395;
t184 = -qJD(3) * t283 - t320 * t395 + t372 * t399;
t133 = t183 * t391 + t184 * t389;
t319 = (-qJDD(1) * t400 + t396 * t505) * t390;
t301 = qJDD(3) + t319;
t162 = -t398 * t353 + t394 * t421;
t508 = qJD(5) * t162;
t68 = t133 * t398 + t301 * t394 - t508;
t635 = t68 / 0.2e1;
t163 = t353 * t394 + t398 * t421;
t69 = qJD(5) * t163 + t133 * t394 - t398 * t301;
t634 = -t69 / 0.2e1;
t132 = -t183 * t389 + t184 * t391;
t131 = qJDD(5) - t132;
t624 = t131 / 0.2e1;
t671 = -Ifges(6,4) + Ifges(7,5);
t669 = Ifges(7,2) + Ifges(6,3);
t711 = Ifges(6,6) - Ifges(7,6);
t726 = Ifges(4,3) + Ifges(5,3);
t725 = pkin(10) * t484 - t652;
t336 = t389 * t395 - t391 * t399;
t254 = t336 * t483;
t325 = t336 * qJD(3);
t724 = t647 + (-t254 + t325) * pkin(10) - t718 * pkin(4);
t592 = -t353 / 0.2e1;
t604 = -t465 / 0.2e1;
t607 = -t194 / 0.2e1;
t615 = -t163 / 0.2e1;
t616 = t162 / 0.2e1;
t617 = -t162 / 0.2e1;
t723 = -Ifges(5,2) * t604 - Ifges(5,6) * t592 + Ifges(6,6) * t616 + Ifges(7,6) * t617 + t607 * t669 + t615 * t670 - t727;
t601 = t421 / 0.2e1;
t606 = t194 / 0.2e1;
t614 = t163 / 0.2e1;
t662 = -t162 * t711 + t163 * t670 + t194 * t669;
t722 = -Ifges(5,4) * t601 + Ifges(6,6) * t617 + Ifges(7,6) * t616 + t606 * t669 + t614 * t670 + t662 / 0.2e1 + t727 + t728;
t148 = t389 * t153;
t81 = t142 * t391 - t148;
t721 = mrSges(5,2) * t195 - t81 * mrSges(5,3);
t675 = m(7) + m(6);
t710 = t131 * t670 + t671 * t69 + t672 * t68;
t719 = t710 / 0.2e1;
t159 = Ifges(6,4) * t162;
t559 = Ifges(7,5) * t162;
t707 = t163 * t672 + t194 * t670 - t159 + t559;
t654 = t729 * t389 + t391 * t730;
t660 = Ifges(5,4) * t465;
t137 = Ifges(5,1) * t421 + t353 * Ifges(5,5) + t660;
t717 = Ifges(5,1) * t601 + Ifges(5,4) * t603 + Ifges(5,5) * t591 + t137 / 0.2e1 + t721;
t716 = -t68 * Ifges(7,5) / 0.2e1 - t131 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t635 + Ifges(6,6) * t624 + (Ifges(7,3) + Ifges(6,2)) * t634;
t633 = t69 / 0.2e1;
t715 = Ifges(6,2) * t634 - Ifges(7,3) * t633 + t711 * t624;
t504 = qJDD(1) * t390;
t697 = pkin(8) * t504 + qJD(2) * t501;
t698 = -pkin(8) * t390 * t505 + pkin(1) * t503;
t228 = -t396 * t697 + t400 * t698;
t205 = -t372 * pkin(2) - t228;
t139 = -t184 * pkin(3) + qJDD(4) + t205;
t227 = t396 * t698 + t400 * t697;
t204 = pkin(9) * t372 + t227;
t209 = -pkin(1) * t504 + pkin(2) * t319 - pkin(9) * t320;
t91 = -qJD(3) * t182 - t204 * t395 + t399 * t209;
t52 = pkin(3) * t301 - qJ(4) * t183 - qJD(4) * t283 + t91;
t509 = qJD(3) * t399;
t90 = t399 * t204 + t395 * t209 - t261 * t510 + t265 * t509;
t59 = qJ(4) * t184 + qJD(4) * t282 + t90;
t21 = t389 * t52 + t391 * t59;
t595 = t301 / 0.2e1;
t622 = t133 / 0.2e1;
t623 = t132 / 0.2e1;
t19 = pkin(10) * t301 + t21;
t31 = -t132 * pkin(4) - t133 * pkin(10) + t139;
t506 = qJD(5) * t398;
t507 = qJD(5) * t394;
t3 = t398 * t19 + t394 * t31 + t92 * t506 - t507 * t78;
t1 = qJ(6) * t131 + qJD(6) * t194 + t3;
t4 = -qJD(5) * t29 - t19 * t394 + t31 * t398;
t2 = -pkin(5) * t131 + qJDD(6) - t4;
t642 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t667 = t131 * t669 + t670 * t68 - t69 * t711;
t713 = t667 / 0.2e1 + t624 * t669 + t635 * t670 + t642 + t139 * mrSges(5,1) - t21 * mrSges(5,3) + Ifges(6,6) * t634 + Ifges(7,6) * t633 + (-t595 - t301 / 0.2e1) * Ifges(5,6) + (-t623 - t132 / 0.2e1) * Ifges(5,2) + (-t622 - t133 / 0.2e1) * Ifges(5,4);
t530 = t390 * t396;
t326 = t392 * t399 - t395 * t530;
t420 = t326 * pkin(3);
t32 = mrSges(6,1) * t131 - mrSges(6,3) * t68;
t33 = -t131 * mrSges(7,1) + t68 * mrSges(7,2);
t709 = t33 - t32;
t34 = -mrSges(6,2) * t131 - mrSges(6,3) * t69;
t35 = -mrSges(7,2) * t69 + mrSges(7,3) * t131;
t708 = t35 + t34;
t384 = pkin(3) * t399 + pkin(2);
t230 = pkin(4) * t336 - pkin(10) * t337 - t384;
t357 = t393 * t395;
t358 = t393 * t399;
t259 = t357 * t389 - t358 * t391;
t649 = t394 * t230 + t398 * t259;
t664 = -qJD(5) * t649 + t394 * t725 + t398 * t724;
t663 = t230 * t506 - t259 * t507 + t394 * t724 - t398 * t725;
t705 = t227 * mrSges(3,2);
t430 = pkin(5) * t394 - qJ(6) * t398;
t85 = t152 * t389 + t526;
t704 = -qJD(6) * t394 + t194 * t430 - t85;
t107 = mrSges(6,1) * t162 + mrSges(6,2) * t163;
t170 = mrSges(5,1) * t353 - mrSges(5,3) * t421;
t703 = t170 - t107;
t653 = pkin(4) * t484 - t654;
t77 = -pkin(4) * t353 - t81;
t36 = pkin(5) * t162 - qJ(6) * t163 + t77;
t442 = mrSges(7,1) * t394 - mrSges(7,3) * t398;
t444 = mrSges(6,1) * t394 + mrSges(6,2) * t398;
t701 = t36 * t442 + t77 * t444;
t700 = -t394 * t711 + t398 * t670;
t558 = Ifges(7,5) * t394;
t561 = Ifges(6,4) * t394;
t699 = t398 * t672 + t558 - t561;
t696 = Ifges(4,5) * t183 + Ifges(5,5) * t133 + Ifges(4,6) * t184 + Ifges(5,6) * t132 + t301 * t726;
t588 = cos(qJ(1));
t485 = t588 * t400;
t397 = sin(qJ(1));
t524 = t396 * t397;
t331 = -t392 * t524 + t485;
t528 = t390 * t399;
t250 = -t331 * t395 + t397 * t528;
t544 = t465 * t398;
t695 = t506 - t544;
t545 = t465 * t394;
t694 = -t507 + t545;
t693 = t3 * t398 - t394 * t4;
t692 = t1 * t398 + t2 * t394;
t158 = Ifges(7,5) * t163;
t71 = Ifges(7,6) * t194 + Ifges(7,3) * t162 + t158;
t562 = Ifges(6,4) * t163;
t74 = -Ifges(6,2) * t162 + Ifges(6,6) * t194 + t562;
t691 = t74 / 0.2e1 - t71 / 0.2e1;
t431 = pkin(5) * t398 + qJ(6) * t394;
t443 = -t398 * mrSges(7,1) - t394 * mrSges(7,3);
t445 = mrSges(6,1) * t398 - mrSges(6,2) * t394;
t690 = m(7) * t431 + mrSges(5,1) - t443 + t445;
t486 = t588 * t396;
t523 = t397 * t400;
t329 = t392 * t486 + t523;
t388 = qJ(3) + pkin(11);
t385 = sin(t388);
t386 = cos(t388);
t487 = t390 * t588;
t240 = t329 * t386 - t385 * t487;
t328 = -t392 * t485 + t524;
t185 = t240 * t394 - t328 * t398;
t689 = t240 * t398 + t328 * t394;
t20 = -t389 * t59 + t391 * t52;
t688 = -t91 * mrSges(4,1) - t20 * mrSges(5,1) + t90 * mrSges(4,2) + t21 * mrSges(5,2);
t687 = t77 * mrSges(6,1) + t36 * mrSges(7,1);
t686 = t77 * mrSges(6,2) - t36 * mrSges(7,3);
t684 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 - t716;
t529 = t390 * t397;
t244 = t331 * t386 + t385 * t529;
t294 = t385 * t392 + t386 * t530;
t683 = -g(1) * t244 - g(2) * t240 - g(3) * t294;
t682 = t624 * t670 + t635 * t672;
t678 = mrSges(5,2) * t139 - mrSges(5,3) * t20 + 0.2e1 * Ifges(5,1) * t622 + 0.2e1 * Ifges(5,4) * t623 + 0.2e1 * Ifges(5,5) * t595;
t610 = t183 / 0.2e1;
t609 = t184 / 0.2e1;
t674 = -t421 / 0.2e1;
t407 = t329 * t395 + t399 * t487;
t403 = t407 * pkin(3);
t673 = mrSges(6,3) + mrSges(7,2);
t666 = -qJ(6) * t718 + qJD(6) * t336 + t663;
t665 = pkin(5) * t718 - t664;
t659 = m(4) * pkin(9) + mrSges(4,3) + mrSges(5,3);
t221 = -t254 * t394 - t398 * t484;
t222 = -t254 * t398 + t394 * t484;
t658 = -pkin(5) * t221 + qJ(6) * t222 - t430 * t325 + (qJD(5) * t431 - qJD(6) * t398) * t337 + t653;
t327 = t392 * t395 + t396 * t528;
t224 = -t391 * t326 + t327 * t389;
t225 = t326 * t389 + t327 * t391;
t374 = pkin(8) * t530;
t587 = pkin(1) * t400;
t295 = t374 + (-pkin(2) - t587) * t392;
t229 = t295 - t420;
t134 = t224 * pkin(4) - t225 * pkin(10) + t229;
t527 = t390 * t400;
t334 = t392 * t396 * pkin(1) + pkin(8) * t527;
t296 = pkin(9) * t392 + t334;
t207 = -t395 * t296 + t399 * t297;
t157 = -pkin(3) * t527 - t327 * qJ(4) + t207;
t208 = t399 * t296 + t395 * t297;
t171 = qJ(4) * t326 + t208;
t101 = t389 * t157 + t391 * t171;
t96 = -pkin(10) * t527 + t101;
t657 = t394 * t134 + t398 * t96;
t462 = mrSges(3,3) * t484;
t655 = -mrSges(3,1) * t373 - mrSges(4,1) * t282 + mrSges(4,2) * t283 + t462;
t414 = -t325 * t394 + t337 * t506;
t651 = t221 - t414;
t522 = t398 * t325;
t413 = t337 * t507 + t522;
t650 = t222 + t413;
t333 = t392 * t587 - t374;
t447 = -mrSges(4,1) * t399 + mrSges(4,2) * t395;
t648 = -m(4) * pkin(2) - t386 * mrSges(5,1) + mrSges(5,2) * t385 + t447;
t646 = -t395 * t91 + t399 * t90;
t645 = mrSges(3,2) - t659;
t644 = mrSges(3,1) - t648;
t315 = qJD(2) * t417;
t317 = t333 * qJD(2);
t143 = -t296 * t510 + t297 * t509 + t395 * t315 + t399 * t317;
t511 = qJD(2) * t390;
t481 = t400 * t511;
t247 = -qJD(3) * t327 - t395 * t481;
t104 = qJ(4) * t247 + qJD(4) * t326 + t143;
t144 = -qJD(3) * t208 + t399 * t315 - t317 * t395;
t248 = qJD(3) * t326 + t399 * t481;
t482 = t396 * t511;
t97 = pkin(3) * t482 - qJ(4) * t248 - qJD(4) * t327 + t144;
t45 = t391 * t104 + t389 * t97;
t43 = pkin(10) * t482 + t45;
t167 = -t391 * t247 + t248 * t389;
t168 = t247 * t389 + t248 * t391;
t318 = t334 * qJD(2);
t203 = -t247 * pkin(3) + t318;
t80 = t167 * pkin(4) - t168 * pkin(10) + t203;
t9 = -qJD(5) * t657 - t394 * t43 + t398 * t80;
t463 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t456 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t640 = t390 ^ 2;
t631 = -t74 / 0.2e1;
t629 = pkin(1) * mrSges(3,1);
t628 = pkin(1) * mrSges(3,2);
t627 = Ifges(4,4) * t610 + Ifges(4,2) * t609 + Ifges(4,6) * t595;
t626 = Ifges(4,1) * t610 + Ifges(4,4) * t609 + Ifges(4,5) * t595;
t619 = -t137 / 0.2e1;
t550 = t283 * Ifges(4,4);
t175 = t282 * Ifges(4,2) + t353 * Ifges(4,6) + t550;
t612 = t175 / 0.2e1;
t268 = Ifges(4,4) * t282;
t176 = t283 * Ifges(4,1) + t353 * Ifges(4,5) + t268;
t611 = t176 / 0.2e1;
t596 = t283 / 0.2e1;
t590 = t392 / 0.2e1;
t589 = t394 / 0.2e1;
t586 = pkin(3) * t283;
t585 = pkin(3) * t389;
t584 = pkin(3) * t391;
t583 = pkin(4) * t386;
t568 = mrSges(6,3) * t162;
t567 = mrSges(6,3) * t163;
t566 = Ifges(3,4) * t396;
t565 = Ifges(3,4) * t400;
t564 = Ifges(4,4) * t395;
t563 = Ifges(4,4) * t399;
t560 = Ifges(6,4) * t398;
t557 = Ifges(7,5) * t398;
t556 = Ifges(3,6) * t373;
t555 = t181 * mrSges(4,3);
t554 = t182 * mrSges(4,3);
t553 = t465 * Ifges(5,6);
t552 = t421 * Ifges(5,5);
t551 = t282 * Ifges(4,6);
t549 = t283 * Ifges(4,5);
t548 = t373 * Ifges(3,5);
t117 = pkin(4) * t421 - pkin(10) * t465 + t586;
t86 = t152 * t391 - t148;
t41 = t394 * t117 + t398 * t86;
t542 = t328 * t385;
t330 = t392 * t523 + t486;
t538 = t330 * t385;
t535 = t337 * t398;
t532 = t386 * t394;
t531 = t386 * t398;
t521 = t398 * t400;
t120 = -mrSges(7,2) * t162 + mrSges(7,3) * t194;
t121 = -mrSges(6,2) * t194 - t568;
t519 = t120 + t121;
t122 = mrSges(6,1) * t194 - t567;
t123 = -mrSges(7,1) * t194 + mrSges(7,2) * t163;
t518 = t122 - t123;
t517 = -t328 * t384 - t329 * t393;
t516 = -t330 * t384 - t331 * t393;
t515 = t588 * pkin(1) + pkin(8) * t529;
t376 = pkin(4) * t527;
t496 = t385 * t527;
t494 = t395 * t529;
t492 = t390 * t521;
t362 = t394 * t527;
t491 = t71 * t589;
t488 = Ifges(3,5) * t320 - Ifges(3,6) * t319 + Ifges(3,3) * t372;
t383 = -pkin(4) - t584;
t382 = pkin(10) + t585;
t479 = t382 * t507;
t478 = t382 * t506;
t470 = -t507 / 0.2e1;
t469 = t506 / 0.2e1;
t468 = -pkin(1) * t397 + pkin(8) * t487;
t44 = -t389 * t104 + t391 * t97;
t58 = -t132 * mrSges(5,1) + t133 * mrSges(5,2);
t100 = t157 * t391 - t389 * t171;
t239 = -t329 * t385 - t386 * t487;
t363 = t395 * t487;
t464 = -t329 * t399 + t363;
t258 = -t391 * t357 - t358 * t389;
t461 = mrSges(3,3) * t483;
t450 = t250 * pkin(3);
t95 = -t100 + t376;
t448 = mrSges(4,1) * t326 - mrSges(4,2) * t327;
t441 = Ifges(4,1) * t399 - t564;
t438 = Ifges(3,2) * t400 + t566;
t437 = -Ifges(4,2) * t395 + t563;
t436 = -Ifges(6,2) * t394 + t560;
t434 = Ifges(4,5) * t399 - Ifges(4,6) * t395;
t432 = Ifges(7,3) * t394 + t557;
t425 = pkin(3) * t494 - t330 * t393 + t331 * t384 + t515;
t40 = t117 * t398 - t394 * t86;
t47 = t134 * t398 - t394 * t96;
t160 = t230 * t398 - t259 * t394;
t18 = -pkin(4) * t301 - t20;
t419 = -pkin(10) * t675 + mrSges(5,2) - t673;
t192 = t394 * t225 + t492;
t8 = t134 * t506 + t394 * t80 + t398 * t43 - t507 * t96;
t412 = pkin(3) * t363 + t328 * t393 - t329 * t384 + t468;
t42 = -pkin(4) * t482 - t44;
t404 = -g(1) * t330 - g(2) * t328 + g(3) * t527;
t368 = Ifges(3,4) * t483;
t342 = t384 * t527;
t335 = t383 - t431;
t332 = (-mrSges(3,1) * t400 + mrSges(3,2) * t396) * t390;
t312 = -t373 * mrSges(3,2) + t461;
t293 = -t385 * t530 + t386 * t392;
t256 = Ifges(3,1) * t484 + t368 + t548;
t255 = t438 * t514 + t556;
t251 = t331 * t399 + t494;
t243 = t331 * t385 - t386 * t529;
t237 = t294 * t394 + t492;
t232 = mrSges(4,1) * t353 - mrSges(4,3) * t283;
t231 = -mrSges(4,2) * t353 + mrSges(4,3) * t282;
t193 = t225 * t398 - t362;
t190 = t244 * t398 + t330 * t394;
t189 = t244 * t394 - t330 * t398;
t174 = t353 * Ifges(4,3) + t549 + t551;
t173 = t337 * t430 + t258;
t169 = -mrSges(5,2) * t353 + mrSges(5,3) * t465;
t155 = -mrSges(4,2) * t301 + mrSges(4,3) * t184;
t154 = mrSges(4,1) * t301 - mrSges(4,3) * t183;
t151 = -pkin(5) * t336 - t160;
t147 = qJ(6) * t336 + t649;
t140 = -mrSges(5,1) * t465 + mrSges(5,2) * t421;
t138 = -mrSges(4,1) * t184 + mrSges(4,2) * t183;
t135 = t353 * Ifges(5,3) + t552 + t553;
t113 = -qJD(5) * t362 + t168 * t394 + t225 * t506 - t398 * t482;
t112 = -qJD(5) * t192 + t398 * t168 + t394 * t482;
t106 = mrSges(7,1) * t162 - mrSges(7,3) * t163;
t105 = pkin(5) * t163 + qJ(6) * t162;
t103 = mrSges(5,1) * t301 - mrSges(5,3) * t133;
t102 = -mrSges(5,2) * t301 + mrSges(5,3) * t132;
t53 = pkin(5) * t192 - qJ(6) * t193 + t95;
t38 = -pkin(5) * t224 - t47;
t37 = qJ(6) * t224 + t657;
t27 = -pkin(5) * t421 - t40;
t26 = qJ(6) * t421 + t41;
t23 = mrSges(6,1) * t69 + mrSges(6,2) * t68;
t22 = mrSges(7,1) * t69 - mrSges(7,3) * t68;
t10 = pkin(5) * t113 - qJ(6) * t112 - qJD(6) * t193 + t42;
t7 = -pkin(5) * t167 - t9;
t6 = qJ(6) * t167 + qJD(6) * t224 + t8;
t5 = pkin(5) * t69 - qJ(6) * t68 - qJD(6) * t163 + t18;
t11 = [(t707 / 0.2e1 + mrSges(7,2) * t24 - mrSges(6,3) * t28 + Ifges(6,4) * t617 + Ifges(7,5) * t616 + t606 * t670 + t614 * t672 + t686) * t112 + t228 * (mrSges(3,1) * t392 - mrSges(3,3) * t530) + m(4) * (t143 * t182 + t144 * t181 + t205 * t295 + t207 * t91 + t208 * t90 + t260 * t318) + m(5) * (t100 * t20 + t101 * t21 + t139 * t229 + t195 * t203 + t44 * t81 + t45 * t82) + m(7) * (t1 * t37 + t10 * t36 + t2 * t38 + t24 * t7 + t25 * t6 + t5 * t53) - t205 * t448 + (Ifges(4,5) * t327 + Ifges(4,6) * t326) * t595 + t488 * t590 + t282 * (Ifges(4,4) * t248 + Ifges(4,2) * t247) / 0.2e1 + (Ifges(3,5) * t590 - t333 * mrSges(3,3) + (t396 * Ifges(3,1) + t565 - t628) * t390) * t320 + ((t548 / 0.2e1 - t313 * mrSges(3,3) + t256 / 0.2e1 + (-t628 + t565 / 0.2e1) * t514) * t400 + (-t556 / 0.2e1 - t255 / 0.2e1 + t174 / 0.2e1 + t135 / 0.2e1 + t551 / 0.2e1 + t549 / 0.2e1 - t182 * mrSges(4,2) + t181 * mrSges(4,1) + t552 / 0.2e1 + t553 / 0.2e1 + t81 * mrSges(5,1) - t82 * mrSges(5,2) - t316 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t353 + (-t629 - t566 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t400) * t514) * t396) * t511 - (Ifges(3,6) * t590 + t334 * mrSges(3,3) + (t438 + t629) * t390) * t319 + (-t181 * t248 + t182 * t247 + t326 * t90 - t327 * t91) * mrSges(4,3) + m(6) * (t18 * t95 + t28 * t9 + t29 * t8 + t3 * t657 + t4 * t47 + t42 * t77) + t657 * t34 + t678 * t225 + (Ifges(4,4) * t327 + Ifges(4,2) * t326) * t609 + (-m(4) * (-pkin(2) * t329 + t468) - t464 * mrSges(4,1) - t407 * mrSges(4,2) - m(3) * t468 + t329 * mrSges(3,1) - mrSges(3,3) * t487 + t397 * mrSges(2,1) + t588 * mrSges(2,2) - m(5) * t412 + t240 * mrSges(5,1) + t463 * t689 - t456 * t185 - t645 * t328 + t419 * t239 + t675 * (pkin(4) * t240 - t412)) * g(1) + (-t588 * mrSges(2,1) - m(3) * t515 - t331 * mrSges(3,1) - m(5) * t425 - t244 * mrSges(5,1) - m(4) * (pkin(2) * t331 + t515) - t251 * mrSges(4,1) - t250 * mrSges(4,2) + (-mrSges(3,3) * t390 + mrSges(2,2)) * t397 - t463 * t190 + t456 * t189 + t645 * t330 + t419 * t243 - t675 * (t244 * pkin(4) + t425)) * g(2) + (Ifges(3,3) * t590 + t333 * mrSges(3,1) - t334 * mrSges(3,2) + (Ifges(3,5) * t396 + Ifges(3,6) * t400) * t390) * t372 + t713 * t224 + (Ifges(4,1) * t327 + Ifges(4,4) * t326) * t610 + t37 * t35 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t640 + t227 * t334 + t228 * t333 - t313 * t318 + t316 * t317) + (-pkin(1) * t332 * t390 + Ifges(2,3)) * qJDD(1) + t317 * t312 + t295 * t138 + t260 * (-mrSges(4,1) * t247 + mrSges(4,2) * t248) + t143 * t231 + t144 * t232 + t229 * t58 + t203 * t140 + (mrSges(6,2) * t18 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t5 + Ifges(6,4) * t634 + Ifges(7,5) * t633 + t682 + t719) * t193 + t207 * t154 + t208 * t155 + t45 * t169 + t44 * t170 + (Ifges(4,5) * t248 + Ifges(4,6) * t247) * t591 + (Ifges(4,1) * t248 + Ifges(4,4) * t247) * t596 + t38 * t33 + t47 * t32 + t53 * t22 + t722 * t167 + t655 * t318 + t248 * t611 + t247 * t612 + t327 * t626 + t326 * t627 + t95 * t23 + t101 * t102 + t100 * t103 + (mrSges(6,1) * t18 + mrSges(7,1) * t5 + t635 * t671 + t684 - t715) * t192 - t392 * t705 + t10 * t106 + t42 * t107 + t6 * t120 + t8 * t121 + t9 * t122 + t7 * t123 + t717 * t168 + (-t25 * mrSges(7,2) - t29 * mrSges(6,3) + Ifges(7,3) * t616 - Ifges(6,2) * t617 + t631 + t71 / 0.2e1 + t671 * t614 - t711 * t606 + t687) * t113 + (mrSges(3,3) * t227 - Ifges(4,5) * t610 - Ifges(5,5) * t622 - Ifges(4,6) * t609 - Ifges(5,6) * t623 - t726 * t595 + t688 - t696 / 0.2e1) * t527; (t18 * t444 + t394 * t684 + t432 * t633 + t436 * t634 + t442 * t5 + t469 * t71 + t624 * t700 + t635 * t699 + t678) * t337 - t81 * (mrSges(5,1) * t484 + mrSges(5,3) * t254) + (t195 * t254 + t484 * t82) * mrSges(5,2) + (-Ifges(5,5) * t254 + Ifges(5,3) * t484) * t592 + (-Ifges(5,4) * t254 + Ifges(5,6) * t484) * t604 + (pkin(1) * (mrSges(3,1) * t396 + mrSges(3,2) * t400) - t396 * (Ifges(3,1) * t400 - t566) / 0.2e1) * qJD(1) ^ 2 * t640 + (-Ifges(7,5) * t413 + Ifges(7,3) * t414) * t616 + t255 * t484 / 0.2e1 + (-t181 * (mrSges(4,1) * t396 - mrSges(4,3) * t520) - t182 * (-mrSges(4,3) * t395 * t400 - mrSges(4,2) * t396)) * t514 + (-t555 + t611) * t509 + (-t175 / 0.2e1 - t554) * t510 + t205 * t447 + (-pkin(2) * t205 - t181 * t216 - t182 * t217) * m(4) + (-t103 + t23) * t258 + (t461 - t312) * t313 + (t282 * t437 + t283 * t441 + t353 * t434) * qJD(3) / 0.2e1 + t649 * t34 + (t160 * t4 + t18 * t258 + t28 * t664 + t29 * t663 + t3 * t649 + t653 * t77) * m(6) + (-Ifges(5,1) * t254 + Ifges(5,5) * t484) * t674 - t705 - ((-Ifges(3,2) * t484 + t399 * t176 + t256 + t368) * t400 + (t174 + t135) * t396 + t353 * (Ifges(4,3) * t396 + t400 * t434) + t283 * (Ifges(4,5) * t396 + t400 * t441) + t282 * (Ifges(4,6) * t396 + t400 * t437) + t373 * (Ifges(3,5) * t400 - Ifges(3,6) * t396)) * t514 / 0.2e1 + t713 * t336 + (-Ifges(6,4) * t413 - Ifges(6,2) * t414) * t617 + t36 * (mrSges(7,1) * t414 + mrSges(7,3) * t413) + t77 * (mrSges(6,1) * t414 - mrSges(6,2) * t413) + (-t413 * t670 - t414 * t711) * t606 + (-Ifges(6,2) * t616 + Ifges(7,3) * t617 - t607 * t711 + t615 * t671 - t687 + t691) * t221 + (Ifges(6,4) * t616 + Ifges(7,5) * t617 + t607 * t670 + t615 * t672 - t686) * t222 + t535 * t719 + t488 - t384 * t58 + (-m(5) * t342 + t332 - t673 * t496 - t675 * (pkin(10) * t496 + t386 * t376 - t393 * t530 + t342) + t456 * (t362 * t386 - t398 * t530) + (t648 * t400 + (m(5) * t393 - t659) * t396 - t463 * (t386 * t521 + t394 * t396)) * t390) * g(3) + t259 * t102 - t217 * t231 - t216 * t232 + t228 * mrSges(3,1) + (-m(5) * t516 + t673 * t538 - t675 * (-pkin(10) * t538 - t330 * t583 + t516) - t463 * (-t330 * t531 + t331 * t394) + t456 * (-t330 * t532 - t331 * t398) + t645 * t331 + t644 * t330) * g(1) + (-m(5) * t517 + t673 * t542 - t675 * (-pkin(10) * t542 - t328 * t583 + t517) - t463 * (-t328 * t531 + t329 * t394) + t456 * (-t328 * t532 - t329 * t398) + t645 * t329 + t644 * t328) * g(2) + t707 * (t337 * t470 - t522 / 0.2e1 - t222 / 0.2e1) + t173 * t22 + (-t413 * t672 + t414 * t671) * t614 + t666 * t120 + (t1 * t147 + t151 * t2 + t173 * t5 + t24 * t665 + t25 * t666 + t36 * t658) * m(7) + t353 * t260 * (mrSges(4,1) * t395 + mrSges(4,2) * t399) + (Ifges(4,5) * t395 + Ifges(4,6) * t399) * t595 + (m(4) * ((-t181 * t399 - t182 * t395) * qJD(3) + t646) - t232 * t509 - t231 * t510 - t395 * t154 + t399 * t155) * pkin(9) + t646 * mrSges(4,3) + t647 * t140 + t414 * t631 + t722 * t324 + (-Ifges(5,4) * t674 - t662 / 0.2e1 + t723) * t253 + (t28 * t650 + t29 * t651 - t4 * t535) * mrSges(6,3) + (t2 * t535 - t24 * t650 + t25 * t651) * mrSges(7,2) + t652 * t169 + t653 * t107 + t654 * t170 + (-t139 * t384 + t195 * t647 - t20 * t258 + t21 * t259 + t652 * t82 + t654 * t81) * m(5) + (-m(4) * t260 + t462 - t655) * t316 + (Ifges(4,2) * t399 + t564) * t609 + (Ifges(4,1) * t395 + t563) * t610 + t454 * t612 - t254 * t619 + t395 * t626 + t399 * t627 + t658 * t106 - (t491 + t717) * t325 - pkin(2) * t138 + t147 * t35 + t151 * t33 + t160 * t32 + t663 * t121 + t664 * t122 + t665 * t123; (Ifges(5,1) * t465 - t661 + t662) * t674 - t283 * (Ifges(4,1) * t282 - t550) / 0.2e1 + t74 * t470 + ((t20 * t391 + t21 * t389) * pkin(3) - t195 * t586 + t81 * t85 - t82 * t86) * m(5) - t18 * t445 + t660 * t604 + t103 * t584 + t102 * t585 + (-t544 / 0.2e1 + t469) * t707 + (-t436 / 0.2e1 + t432 / 0.2e1) * t508 - (-Ifges(4,2) * t283 + t176 + t268) * t282 / 0.2e1 + t5 * t443 + t682 * t394 - t140 * t586 + (t163 * t699 + t194 * t700) * qJD(5) / 0.2e1 + (t18 * t383 - t28 * t40 - t29 * t41 - t77 * t85) * m(6) + (-t24 * t27 - t25 * t26 + t335 * t5 + t704 * t36) * m(7) + (t708 * t398 + t709 * t394 + ((-t28 * t398 - t29 * t394) * qJD(5) + t693) * m(6) + ((t24 * t398 - t25 * t394) * qJD(5) + t692) * m(7)) * t382 - t688 + t696 + t283 * t554 + t282 * t555 + (t478 - t27) * t123 + (-t557 + t560) * t635 + (-t478 - t40) * t122 + (-t479 - t41) * t121 + (-t479 - t26) * t120 + t383 * t23 + t335 * t22 - t260 * (mrSges(4,1) * t283 + mrSges(4,2) * t282) - t181 * t231 + t182 * t232 + (m(5) * t403 + mrSges(4,1) * t407 - mrSges(4,2) * t464 + t240 * mrSges(5,2) - t239 * t690 + t675 * (-t239 * pkin(4) - t240 * pkin(10) + t403)) * g(2) - t86 * t169 + (Ifges(4,5) * t282 - Ifges(4,6) * t283) * t592 + t175 * t596 + t691 * t545 + (t24 * t695 + t25 * t694 + t683 + t692) * mrSges(7,2) + (-t28 * t695 + t29 * t694 + t683 + t693) * mrSges(6,3) + (t491 + t701) * qJD(5) + t703 * t85 + (Ifges(5,5) * t592 + t432 * t617 + t436 * t616 + t607 * t700 + t615 * t699 + t619 - t701 - t721) * t465 + t723 * t421 + t704 * t106 + t558 * t633 + t561 * t634 + (t715 + t716) * t398 + t710 * t589 + (-m(5) * t450 - mrSges(4,1) * t250 + mrSges(4,2) * t251 + t244 * mrSges(5,2) - t675 * (-t243 * pkin(4) + pkin(10) * t244 + t450) + t690 * t243) * g(1) + (-m(5) * t420 + t294 * mrSges(5,2) - t448 - t675 * (t293 * pkin(4) + pkin(10) * t294 + t420) - t690 * t293) * g(3); -t465 * t169 + (-t106 + t703) * t421 + (t194 * t519 - t709) * t398 + (-t194 * t518 + t708) * t394 + t58 + (t1 * t394 - t421 * t36 - t2 * t398 + t404 + t194 * (t24 * t394 + t25 * t398)) * m(7) + (-t421 * t77 + t3 * t394 + t4 * t398 + t404 + t194 * (-t28 * t394 + t29 * t398)) * m(6) + (t421 * t81 - t465 * t82 + t139 + t404) * m(5); (t189 * t463 + t190 * t456) * g(1) + (t456 * (t294 * t398 - t362) + t463 * t237) * g(3) + (-t162 * t670 - t163 * t711) * t607 + (t463 * t185 + t456 * t689) * g(2) + (-pkin(5) * t2 + qJ(6) * t1 - t105 * t36 - t24 * t29 + t25 * t656) * m(7) + (-t162 * t672 + t158 - t562 + t71) * t615 + (t162 * t24 + t163 * t25) * mrSges(7,2) + qJ(6) * t35 + (-Ifges(6,2) * t163 - t159 + t707) * t616 + (t518 + t567) * t29 + (-t519 - t568) * t28 - pkin(5) * t33 - t77 * (mrSges(6,1) * t163 - mrSges(6,2) * t162) - t36 * (mrSges(7,1) * t163 + mrSges(7,3) * t162) + t642 + t74 * t614 + (Ifges(7,3) * t163 - t559) * t617 - t105 * t106 + qJD(6) * t120 + t667; t163 * t106 - t194 * t120 + (-g(1) * t189 - g(2) * t185 - g(3) * t237 + t163 * t36 - t194 * t25 + t2) * m(7) + t33;];
tau  = t11;
