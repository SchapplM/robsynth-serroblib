% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:53
% EndTime: 2019-03-09 21:43:31
% DurationCPUTime: 61.98s
% Computational Cost: add. (15565->1051), mult. (37792->1319), div. (0->0), fcn. (29564->10), ass. (0->452)
t373 = sin(qJ(3));
t376 = cos(qJ(3));
t370 = cos(pkin(6));
t471 = qJD(1) * t370 + qJD(2);
t374 = sin(qJ(2));
t369 = sin(pkin(6));
t538 = qJD(1) * t369;
t506 = t374 * t538;
t256 = t373 * t471 + t376 * t506;
t377 = cos(qJ(2));
t529 = qJD(1) * qJD(2);
t295 = (qJDD(1) * t374 + t377 * t529) * t369;
t527 = qJDD(1) * t370;
t351 = qJDD(2) + t527;
t153 = -qJD(3) * t256 - t373 * t295 + t376 * t351;
t619 = t153 / 0.2e1;
t372 = sin(qJ(4));
t375 = cos(qJ(4));
t255 = -t373 * t506 + t376 * t471;
t380 = qJD(3) * t255 + t376 * t295 + t373 * t351;
t294 = (-qJDD(1) * t377 + t374 * t529) * t369;
t387 = qJDD(3) + t294;
t505 = t377 * t538;
t422 = -qJD(3) + t505;
t200 = t375 * t256 - t372 * t422;
t535 = qJD(4) * t200;
t76 = t372 * t380 - t375 * t387 + t535;
t628 = t76 / 0.2e1;
t629 = -t76 / 0.2e1;
t406 = t375 * t422;
t534 = qJD(4) * t372;
t75 = qJD(4) * t406 + t256 * t534 - t372 * t387 - t375 * t380;
t630 = t75 / 0.2e1;
t631 = -t75 / 0.2e1;
t700 = Ifges(7,5) + Ifges(5,5);
t701 = Ifges(6,5) + Ifges(7,4);
t707 = t387 / 0.2e1;
t148 = qJDD(4) - t153;
t248 = qJD(4) - t255;
t551 = t369 * t377;
t357 = pkin(8) * t551;
t597 = pkin(1) * t374;
t236 = qJD(2) * pkin(9) + (t357 + (pkin(9) + t597) * t370) * qJD(1);
t246 = (-pkin(2) * t377 - pkin(9) * t374 - pkin(1)) * t538;
t150 = t376 * t236 + t373 * t246;
t123 = -pkin(10) * t422 + t150;
t528 = qJDD(1) * t369;
t720 = pkin(1) * t370 * t529 + pkin(8) * t528;
t721 = -pkin(8) * t369 * t529 + pkin(1) * t527;
t202 = -t374 * t720 + t377 * t721;
t178 = -t351 * pkin(2) - t202;
t50 = -t153 * pkin(3) - pkin(10) * t380 + t178;
t532 = qJD(4) * t375;
t552 = t369 * t374;
t352 = pkin(8) * t552;
t596 = pkin(1) * t377;
t275 = t352 + (-pkin(2) - t596) * t370;
t235 = -qJD(2) * pkin(2) + qJD(1) * t275;
t381 = -t255 * pkin(3) - t256 * pkin(10) + t235;
t201 = t374 * t721 + t377 * t720;
t177 = pkin(9) * t351 + t201;
t521 = pkin(1) * t528;
t185 = pkin(2) * t294 - pkin(9) * t295 - t521;
t530 = t376 * qJD(3);
t531 = t373 * qJD(3);
t52 = t376 * t177 + t373 * t185 - t236 * t531 + t246 * t530;
t725 = pkin(10) * t387 + qJD(4) * t381 + t52;
t7 = -t123 * t532 - t372 * t725 + t375 * t50;
t403 = qJDD(5) - t7;
t583 = pkin(4) + qJ(6);
t1 = -pkin(5) * t75 - qJD(6) * t248 - t148 * t583 + t403;
t6 = -t123 * t534 + t372 * t50 + t375 * t725;
t4 = -qJ(5) * t148 - qJD(5) * t248 - t6;
t2 = -pkin(5) * t76 + qJDD(6) - t4;
t5 = -pkin(4) * t148 + t403;
t711 = t7 * mrSges(5,1) - t6 * mrSges(5,2) + t5 * mrSges(6,2) + t2 * mrSges(7,2) - t4 * mrSges(6,3) - t1 * mrSges(7,3);
t708 = t380 / 0.2e1;
t745 = Ifges(4,4) * t708;
t747 = Ifges(6,4) * t630 - 0.2e1 * Ifges(4,2) * t619 - 0.2e1 * Ifges(4,6) * t707 + Ifges(5,6) * t629 + t701 * t628 + t700 * t631 + t711 - t745;
t359 = t370 * t597;
t453 = pkin(3) * t373 - pkin(10) * t376;
t746 = (t359 + (pkin(8) + t453) * t551) * qJD(1) - t453 * qJD(3);
t620 = t148 / 0.2e1;
t703 = Ifges(5,1) + Ifges(7,3);
t699 = Ifges(7,2) + Ifges(6,3);
t744 = Ifges(6,4) - t700;
t743 = -Ifges(5,6) + t701;
t545 = t376 * t377;
t258 = (t372 * t374 + t375 * t545) * t369;
t245 = qJD(1) * t258;
t533 = qJD(4) * t373;
t742 = t372 * t533 + t245;
t558 = t255 * t372;
t659 = t534 - t558;
t740 = -t148 / 0.2e1;
t702 = -Ifges(5,4) + Ifges(7,6);
t738 = Ifges(6,6) - Ifges(7,6);
t666 = t370 * t596 - t352;
t288 = t666 * qJD(1);
t413 = (pkin(2) * t374 - pkin(9) * t377) * t369;
t289 = qJD(1) * t413;
t187 = -t373 * t288 + t289 * t376;
t164 = -pkin(3) * t506 - t187;
t365 = pkin(9) * t530;
t737 = -t164 + t365;
t557 = t255 * t375;
t660 = t532 - t557;
t188 = t376 * t288 + t373 * t289;
t165 = pkin(10) * t506 + t188;
t595 = pkin(3) * t376;
t327 = -pkin(10) * t373 - pkin(2) - t595;
t735 = -t375 * t165 + t327 * t532 - t372 * t746;
t559 = qJ(5) * t372;
t594 = pkin(4) * t375;
t734 = t559 + t594;
t733 = mrSges(7,1) * t2 + Ifges(5,6) * t620 + t701 * t740 + (Ifges(5,4) + t738) * t631 + (Ifges(5,2) + t699) * t629;
t55 = t123 * t372 - t375 * t381;
t414 = pkin(5) * t200 + t55;
t717 = qJD(5) + t414;
t29 = -t248 * t583 + t717;
t56 = t375 * t123 + t372 * t381;
t47 = -t248 * qJ(5) - t56;
t199 = t256 * t372 + t406;
t716 = t199 * pkin(5) - qJD(6);
t30 = -t47 - t716;
t676 = -qJD(5) - t55;
t46 = -pkin(4) * t248 - t676;
t690 = t422 * Ifges(4,6);
t732 = -t235 * mrSges(4,1) + t55 * mrSges(5,1) + t56 * mrSges(5,2) - t46 * mrSges(6,2) - t30 * mrSges(7,2) + t47 * mrSges(6,3) + t29 * mrSges(7,3) - t690 / 0.2e1;
t651 = Ifges(5,3) + Ifges(7,1) + Ifges(6,1);
t53 = -t373 * t177 + t376 * t185 - t236 * t530 - t246 * t531;
t45 = -pkin(3) * t387 - t53;
t378 = t75 * qJ(5) - t200 * qJD(5) + t45;
t3 = t199 * qJD(6) + t583 * t76 + t378;
t729 = t3 * mrSges(7,3);
t649 = t199 * t743 - t200 * t744 + t248 * t651;
t728 = t649 / 0.2e1;
t195 = Ifges(7,6) * t200;
t567 = t200 * Ifges(6,6);
t695 = t199 * t699 + t248 * t701 + t195 - t567;
t197 = Ifges(5,4) * t199;
t571 = Ifges(7,6) * t199;
t694 = t200 * t703 + t248 * t700 - t197 + t571;
t727 = Ifges(4,2) * t255;
t446 = mrSges(6,2) * t375 - mrSges(6,3) * t372;
t448 = mrSges(5,1) * t375 - mrSges(5,2) * t372;
t493 = m(7) * qJ(6) + mrSges(7,3);
t564 = t372 * mrSges(7,2);
t646 = t375 * t493 - t446 + t448 + t564;
t726 = mrSges(4,1) + t646;
t463 = t373 * t505;
t724 = -qJ(5) * t463 - qJD(5) * t376 + t735;
t546 = t375 * t376;
t362 = pkin(9) * t546;
t723 = qJD(4) * t362 - t372 * t165 + t327 * t534 + t375 * t746;
t718 = t372 * t530 + t373 * t532;
t722 = pkin(4) * t718 + qJ(5) * t742 + t737;
t719 = pkin(4) * t659 - qJD(5) * t372 - t150;
t599 = sin(qJ(1));
t507 = t599 * t377;
t600 = cos(qJ(1));
t510 = t600 * t374;
t305 = t370 * t510 + t507;
t512 = t369 * t600;
t224 = t305 * t376 - t373 * t512;
t508 = t599 * t374;
t509 = t600 * t377;
t304 = -t370 * t509 + t508;
t166 = t224 * t372 - t304 * t375;
t167 = t224 * t375 + t304 * t372;
t714 = mrSges(7,1) * t1 + Ifges(6,4) * t740 + Ifges(6,6) * t629 + t700 * t620 + t702 * t628 + (Ifges(6,2) + t703) * t631;
t713 = mrSges(6,1) * t4 - mrSges(5,3) * t6 - t733;
t712 = mrSges(6,1) * t5 - mrSges(5,3) * t7 + t714;
t625 = m(6) + m(7);
t706 = pkin(9) * (t375 * t531 + t376 * t534);
t705 = -mrSges(6,1) - mrSges(5,3);
t704 = -mrSges(4,3) + mrSges(3,2);
t693 = Ifges(4,3) * t422;
t691 = t422 * Ifges(4,5);
t624 = pkin(5) + pkin(10);
t461 = m(7) * t624 + mrSges(7,1);
t689 = mrSges(4,2) - t461;
t513 = -pkin(9) * t372 - pkin(4);
t688 = t463 * t583 + qJD(6) * t376 + (-qJ(6) + t513) * t531 + t723 + (qJD(3) * t546 - t742) * pkin(5);
t517 = t372 * t551;
t468 = t376 * t517;
t244 = qJD(1) * t468 - t375 * t506;
t549 = t372 * t376;
t361 = pkin(9) * t549;
t548 = t373 * t375;
t687 = t244 * pkin(5) + (-pkin(5) * t548 - t361) * qJD(4) + (-pkin(5) * t549 + (-pkin(9) * t375 + qJ(5)) * t373) * qJD(3) + t724;
t420 = -qJ(5) * t375 + qJ(6) * t372;
t686 = -t244 * t583 + t420 * t530 + (qJD(6) * t372 + (qJ(6) * qJD(4) - qJD(5)) * t375) * t373 + t722;
t685 = -qJ(5) * t531 + t706 - t724;
t684 = pkin(4) * t463 + t513 * t531 + t723;
t683 = -pkin(4) * t244 + (-qJ(5) * t530 - qJD(5) * t373) * t375 + t722;
t682 = -qJD(6) * t375 + t248 * t420 + t719;
t681 = -qJ(5) * t660 + t719;
t149 = -t373 * t236 + t376 * t246;
t182 = pkin(3) * t256 - pkin(10) * t255;
t83 = t375 * t149 + t372 * t182;
t63 = -qJ(5) * t256 - t83;
t680 = pkin(5) * t558 - t624 * t534 + t63;
t130 = t372 * t149;
t335 = t624 * t375;
t679 = qJD(4) * t335 - t130 - (pkin(5) * t255 - t182) * t375 + t583 * t256;
t449 = mrSges(4,1) * t376 - mrSges(4,2) * t373;
t678 = t449 + mrSges(3,1);
t677 = t422 * (Ifges(4,5) * t376 - Ifges(4,6) * t373);
t675 = -t56 + t716;
t674 = -t706 + t735;
t524 = pkin(9) * t531;
t673 = t372 * t524 - t723;
t473 = -t305 * t373 - t376 * t512;
t672 = t734 * t473;
t307 = -t370 * t508 + t509;
t511 = t369 * t599;
t227 = t307 * t373 - t376 * t511;
t671 = t734 * t227;
t302 = -t370 * t376 + t373 * t552;
t537 = qJD(2) * t377;
t503 = t369 * t537;
t220 = -qJD(3) * t302 + t376 * t503;
t670 = -qJD(4) * t551 + t220;
t669 = t244 - t718;
t668 = -t375 * t530 + t742;
t667 = t734 * t302;
t421 = mrSges(5,1) - mrSges(6,2) + t493;
t665 = pkin(4) * t625 + t421;
t569 = Ifges(7,6) * t375;
t572 = Ifges(6,6) * t375;
t664 = t372 * t699 + t569 - t572;
t570 = Ifges(7,6) * t372;
t573 = Ifges(6,6) * t372;
t663 = t375 * t699 - t570 + t573;
t575 = Ifges(5,4) * t372;
t662 = t375 * t703 + t570 - t575;
t661 = t463 - t531;
t658 = -t373 * t53 + t376 * t52;
t657 = -t372 * t7 + t375 * t6;
t656 = t372 * t5 - t375 * t4;
t580 = Ifges(3,4) * t374;
t636 = t369 ^ 2;
t655 = (pkin(1) * (mrSges(3,1) * t374 + mrSges(3,2) * t377) - t374 * (Ifges(3,1) * t377 - t580) / 0.2e1) * t636;
t654 = -mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t650 = t148 * t651 + t743 * t76 + t744 * t75;
t122 = pkin(3) * t422 - t149;
t383 = -t200 * qJ(5) + t122;
t31 = t199 * t583 + t383;
t444 = -mrSges(7,2) * t375 + mrSges(7,3) * t372;
t445 = -mrSges(6,2) * t372 - mrSges(6,3) * t375;
t447 = mrSges(5,1) * t372 + mrSges(5,2) * t375;
t54 = t199 * pkin(4) + t383;
t648 = t122 * t447 + t31 * t444 + t54 * t445;
t647 = t372 * t743 - t375 * t744;
t643 = t689 + t705;
t467 = mrSges(3,3) * t506;
t642 = -m(4) * t235 + mrSges(3,1) * t471 + mrSges(4,1) * t255 - mrSges(4,2) * t256 - t467;
t641 = m(7) * pkin(5) + mrSges(7,1) - t705;
t639 = Ifges(4,1) * t708 + Ifges(4,5) * t707;
t228 = t307 * t376 + t373 * t511;
t303 = t370 * t373 + t376 * t552;
t638 = -g(1) * t228 - g(2) * t224 - g(3) * t303;
t632 = Ifges(4,4) * t619 + t639;
t196 = Ifges(6,6) * t199;
t91 = Ifges(6,4) * t248 - Ifges(6,2) * t200 + t196;
t627 = t91 / 0.2e1;
t568 = t200 * Ifges(5,4);
t92 = -Ifges(5,2) * t199 + Ifges(5,6) * t248 + t568;
t626 = t92 / 0.2e1;
t578 = Ifges(4,4) * t256;
t141 = t578 - t690 + t727;
t621 = -t141 / 0.2e1;
t618 = -t199 / 0.2e1;
t617 = t199 / 0.2e1;
t616 = -t200 / 0.2e1;
t615 = t200 / 0.2e1;
t609 = -t248 / 0.2e1;
t608 = t248 / 0.2e1;
t607 = -t255 / 0.2e1;
t606 = -t256 / 0.2e1;
t605 = t256 / 0.2e1;
t598 = pkin(1) * t369;
t593 = pkin(10) * t227;
t589 = t473 * pkin(10);
t366 = t373 * pkin(9);
t582 = mrSges(5,3) * t199;
t581 = mrSges(5,3) * t200;
t579 = Ifges(3,4) * t377;
t577 = Ifges(4,4) * t373;
t576 = Ifges(4,4) * t376;
t574 = Ifges(5,4) * t375;
t566 = t255 * mrSges(4,3);
t565 = t256 * mrSges(4,3);
t563 = t373 * t45;
t560 = qJ(5) * t199;
t555 = t304 * t373;
t306 = t370 * t507 + t510;
t553 = t306 * t373;
t550 = t372 * t373;
t547 = t373 * t377;
t125 = mrSges(6,1) * t199 - mrSges(6,3) * t248;
t126 = -mrSges(7,1) * t199 + mrSges(7,2) * t248;
t544 = t126 - t125;
t296 = t302 * pkin(3);
t474 = t303 * pkin(10) - t296;
t156 = t275 - t474;
t541 = t357 + t359;
t276 = pkin(9) * t370 + t541;
t542 = pkin(2) * t551 + pkin(9) * t552;
t277 = -t542 - t598;
t184 = t376 * t276 + t373 * t277;
t158 = -pkin(10) * t551 + t184;
t80 = t372 * t156 + t375 * t158;
t540 = pkin(4) * t550 + t366;
t267 = t372 * t327 + t362;
t539 = t600 * pkin(1) + pkin(8) * t511;
t355 = pkin(3) * t551;
t523 = pkin(10) * t534;
t522 = pkin(10) * t532;
t518 = t369 * t547;
t516 = Ifges(4,5) * t380 + Ifges(4,6) * t153 + Ifges(4,3) * t387;
t514 = Ifges(3,5) * t295 - Ifges(3,6) * t294 + Ifges(3,3) * t351;
t504 = qJD(2) * t552;
t496 = t552 / 0.2e1;
t491 = t538 / 0.2e1;
t486 = -t534 / 0.2e1;
t485 = t534 / 0.2e1;
t484 = -t532 / 0.2e1;
t483 = t532 / 0.2e1;
t480 = t530 / 0.2e1;
t37 = -t75 * mrSges(6,1) + t148 * mrSges(6,2);
t36 = -t76 * mrSges(7,1) + t148 * mrSges(7,2);
t34 = -t75 * mrSges(7,1) - t148 * mrSges(7,3);
t479 = -pkin(3) - t559;
t478 = -t304 * pkin(2) + pkin(9) * t305;
t477 = -t306 * pkin(2) + pkin(9) * t307;
t213 = t473 * pkin(3);
t476 = pkin(10) * t224 + t213;
t215 = t227 * pkin(3);
t475 = pkin(10) * t228 - t215;
t82 = t182 * t375 - t130;
t79 = t156 * t375 - t372 * t158;
t183 = -t373 * t276 + t277 * t376;
t266 = t327 * t375 - t361;
t466 = mrSges(3,3) * t505;
t465 = pkin(10) * t518 + t376 * t355 + t542;
t455 = -pkin(1) * t599 + pkin(8) * t512;
t65 = -qJ(5) * t302 - t80;
t239 = qJ(5) * t376 - t267;
t157 = -t183 + t355;
t450 = mrSges(4,1) * t302 + mrSges(4,2) * t303;
t443 = Ifges(4,1) * t376 - t577;
t441 = Ifges(5,1) * t372 + t574;
t440 = -Ifges(4,2) * t373 + t576;
t439 = -Ifges(5,2) * t372 + t574;
t438 = Ifges(5,2) * t375 + t575;
t436 = Ifges(6,4) * t372 + Ifges(6,5) * t375;
t435 = Ifges(7,4) * t375 - Ifges(7,5) * t372;
t431 = Ifges(5,5) * t372 + Ifges(5,6) * t375;
t430 = -Ifges(6,2) * t375 + t573;
t429 = Ifges(6,2) * t372 + t572;
t424 = -Ifges(7,3) * t372 + t569;
t417 = -pkin(10) * t555 - t304 * t595 + t478;
t416 = -pkin(10) * t553 - t306 * t595 + t477;
t415 = t307 * pkin(2) + pkin(9) * t306 + t539;
t219 = qJD(3) * t303 + t373 * t503;
t293 = t541 * qJD(2);
t113 = t219 * pkin(3) - t220 * pkin(10) + t293;
t290 = qJD(2) * t413;
t292 = t666 * qJD(2);
t103 = -t276 * t531 + t277 * t530 + t373 * t290 + t376 * t292;
t99 = pkin(10) * t504 + t103;
t24 = t113 * t375 - t156 * t534 - t158 * t532 - t372 * t99;
t104 = -t276 * t530 - t277 * t531 + t290 * t376 - t373 * t292;
t412 = -qJ(5) * t625 + t654;
t409 = t228 * pkin(3) + t415;
t407 = t235 * (mrSges(4,1) * t373 + mrSges(4,2) * t376);
t23 = t372 * t113 + t156 * t532 - t158 * t534 + t375 * t99;
t170 = t228 * t372 - t306 * t375;
t221 = t303 * t372 + t375 * t551;
t401 = -g(1) * t170 - g(2) * t166 - g(3) * t221;
t393 = -t305 * pkin(2) - t304 * pkin(9) + t455;
t212 = t221 * pkin(4);
t222 = t303 * t375 - t517;
t81 = -qJ(5) * t222 + t157 + t212;
t392 = -pkin(3) * t224 + t393;
t389 = Ifges(3,6) * t370 + (Ifges(3,2) * t377 + t580) * t369;
t171 = t228 * t375 + t306 * t372;
t388 = t171 * pkin(4) + qJ(5) * t170 + t409;
t100 = -pkin(3) * t504 - t104;
t385 = t369 * t471 * (Ifges(3,5) * t377 - Ifges(3,6) * t374);
t12 = -qJ(5) * t219 - qJD(5) * t302 - t23;
t384 = -pkin(4) * t167 - qJ(5) * t166 + t392;
t118 = -t303 * t534 + t372 * t504 + t375 * t670;
t382 = -qJ(5) * t118 - qJD(5) * t222 + t100;
t367 = t376 * pkin(4);
t346 = Ifges(3,4) * t505;
t334 = t624 * t372;
t322 = t479 - t594;
t312 = -t375 * t583 + t479;
t308 = (-mrSges(3,1) * t377 + mrSges(3,2) * t374) * t369;
t291 = t541 * qJD(1);
t287 = -mrSges(3,2) * t471 + t466;
t280 = -qJ(5) * t548 + t540;
t257 = -t375 * t552 + t468;
t247 = Ifges(4,4) * t255;
t240 = -t266 + t367;
t237 = t373 * t420 + t540;
t233 = Ifges(3,1) * t506 + Ifges(3,5) * t471 + t346;
t232 = Ifges(3,6) * qJD(2) + qJD(1) * t389;
t210 = -pkin(5) * t550 - t239;
t205 = -mrSges(4,1) * t422 - t565;
t204 = mrSges(4,2) * t422 + t566;
t203 = qJ(6) * t376 + t361 + t367 + (pkin(5) * t373 - t327) * t375;
t194 = -t306 * t546 + t307 * t372;
t193 = -t306 * t549 - t307 * t375;
t192 = -t304 * t546 + t305 * t372;
t191 = -t304 * t549 - t305 * t375;
t142 = Ifges(4,1) * t256 + t247 - t691;
t140 = Ifges(4,5) * t256 + Ifges(4,6) * t255 - t693;
t129 = mrSges(5,1) * t248 - t581;
t128 = -mrSges(5,2) * t248 - t582;
t127 = mrSges(6,1) * t200 + mrSges(6,2) * t248;
t124 = mrSges(7,1) * t200 - mrSges(7,3) * t248;
t117 = t303 * t532 + t372 * t670 - t375 * t504;
t112 = -mrSges(4,2) * t387 + t153 * mrSges(4,3);
t111 = mrSges(4,1) * t387 - mrSges(4,3) * t380;
t109 = -mrSges(6,2) * t199 - mrSges(6,3) * t200;
t108 = mrSges(5,1) * t199 + mrSges(5,2) * t200;
t107 = pkin(4) * t200 + t560;
t106 = -mrSges(7,2) * t200 + mrSges(7,3) * t199;
t78 = -t153 * mrSges(4,1) + mrSges(4,2) * t380;
t77 = t200 * t583 + t560;
t66 = -pkin(4) * t302 - t79;
t64 = -pkin(4) * t256 - t82;
t58 = qJ(6) * t221 + t81;
t49 = -pkin(5) * t221 - t65;
t43 = pkin(5) * t222 - t302 * t583 - t79;
t39 = -mrSges(5,2) * t148 - mrSges(5,3) * t76;
t38 = mrSges(5,1) * t148 + mrSges(5,3) * t75;
t35 = mrSges(6,1) * t76 - mrSges(6,3) * t148;
t28 = mrSges(7,2) * t75 + mrSges(7,3) * t76;
t27 = mrSges(5,1) * t76 - mrSges(5,2) * t75;
t26 = -mrSges(6,2) * t76 + mrSges(6,3) * t75;
t25 = pkin(4) * t117 + t382;
t13 = -pkin(4) * t219 - t24;
t11 = qJD(6) * t221 + t117 * t583 + t382;
t10 = -pkin(5) * t117 - t12;
t9 = pkin(5) * t118 - qJD(6) * t302 - t219 * t583 - t24;
t8 = t76 * pkin(4) + t378;
t14 = [(t385 / 0.2e1 + t140 * t496) * qJD(2) + (Ifges(4,5) * t303 - Ifges(4,3) * t551) * t707 - t422 * (Ifges(4,5) * t220 + Ifges(4,3) * t504) / 0.2e1 + t255 * (Ifges(4,4) * t220 + Ifges(4,6) * t504) / 0.2e1 + (Ifges(4,4) * t303 - Ifges(4,6) * t551) * t619 + (-t150 * t504 + t220 * t235 + t52 * t551) * mrSges(4,2) + (-t201 * t370 - t295 * t598 - t351 * t541) * mrSges(3,2) + (Ifges(3,1) * t295 - Ifges(3,4) * t294 + Ifges(3,5) * t351) * t496 + (Ifges(3,4) * t295 - Ifges(3,2) * t294 + Ifges(3,6) * t351) * t551 / 0.2e1 + (-m(3) * t288 - t642) * t293 + (t369 * t233 + t636 * qJD(1) * (-Ifges(3,2) * t374 + t579)) * t537 / 0.2e1 + (-t745 + t650 / 0.2e1 + t651 * t620 - t52 * mrSges(4,3) + t747) * t302 + (-m(7) * t388 - m(4) * t415 - t228 * mrSges(4,1) - mrSges(2,1) * t600 + mrSges(2,2) * t599 - m(5) * (t409 + t593) - m(6) * (t388 + t593) - m(3) * t539 - t307 * mrSges(3,1) - mrSges(3,3) * t511 + t704 * t306 - t421 * t171 + t654 * t170 + t643 * t227) * g(2) - t232 * t504 / 0.2e1 + (Ifges(4,1) * t303 - Ifges(4,5) * t551) * t708 + (Ifges(4,1) * t220 + Ifges(4,5) * t504) * t605 + t149 * (mrSges(4,1) * t504 - mrSges(4,3) * t220) - t308 * t521 + (t700 * t615 + t651 * t608 - Ifges(4,4) * t605 + Ifges(6,4) * t616 + Ifges(5,6) * t618 + t701 * t617 + t728 + t621 - t150 * mrSges(4,3) - t727 / 0.2e1 - t732) * t219 + m(6) * (t12 * t47 + t13 * t46 + t25 * t54 + t4 * t65 + t5 * t66 + t8 * t81) + m(7) * (t1 * t43 + t10 * t30 + t11 * t31 + t2 * t49 + t29 * t9 + t3 * t58) + m(5) * (t100 * t122 + t157 * t45 + t23 * t56 - t24 * t55 + t6 * t80 + t7 * t79) + t178 * t450 + (-m(3) * t455 + t305 * mrSges(3,1) - mrSges(3,3) * t512 - m(4) * t393 + t224 * mrSges(4,1) - m(7) * t384 + mrSges(2,1) * t599 + mrSges(2,2) * t600 - m(5) * (t392 + t589) - m(6) * (t384 + t589) - t704 * t304 + t421 * t167 - t654 * t166 + t643 * t473) * g(1) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t636 + t201 * t541 + t202 * t666 + t291 * t292) + (t201 * t551 - t202 * t552 - t288 * t503 - t291 * t504 - t294 * t541 - t295 * t666) * mrSges(3,3) + (t202 * t370 - t294 * t598 + t351 * t666) * mrSges(3,1) + (t45 * mrSges(5,1) - t8 * mrSges(6,2) - Ifges(5,2) * t629 + Ifges(6,6) * t630 + t620 * t743 + t699 * t628 + t702 * t631 + t713 + t729) * t221 + (t47 * mrSges(6,1) - t30 * mrSges(7,1) - t56 * mrSges(5,3) + t122 * mrSges(5,1) - t54 * mrSges(6,2) - t92 / 0.2e1 + t31 * mrSges(7,3) + t699 * t617 + t702 * t615 + t695 / 0.2e1 + t743 * t608 + Ifges(6,6) * t616 - Ifges(5,2) * t618) * t117 + (t55 * mrSges(5,3) + t29 * mrSges(7,1) + t46 * mrSges(6,1) + t122 * mrSges(5,2) - t54 * mrSges(6,3) - t91 / 0.2e1 - t31 * mrSges(7,2) - t738 * t617 + t703 * t615 + t694 / 0.2e1 - t744 * t608 - Ifges(6,2) * t616 + Ifges(5,4) * t618) * t118 + (t45 * mrSges(5,2) - t3 * mrSges(7,2) - t8 * mrSges(6,3) + Ifges(5,4) * t629 - Ifges(6,2) * t630 - t620 * t744 - t628 * t738 + t703 * t631 + t712) * t222 + t370 * t514 / 0.2e1 + m(4) * (t103 * t150 + t104 * t149 + t178 * t275 + t183 * t53 + t184 * t52) + t351 * (Ifges(3,3) * t370 + (Ifges(3,5) * t374 + Ifges(3,6) * t377) * t369) / 0.2e1 + t292 * t287 + t275 * t78 + Ifges(2,3) * qJDD(1) + t220 * t142 / 0.2e1 + t103 * t204 + t104 * t205 + t183 * t111 + t184 * t112 + t157 * t27 + t23 * t128 + t24 * t129 + t10 * t126 + t13 * t127 + t9 * t124 + t12 * t125 + t100 * t108 + t25 * t109 + t11 * t106 + t80 * t39 + t81 * t26 + t79 * t38 + t65 * t35 + t66 * t37 + t58 * t28 + t49 * t36 + t43 * t34 + t295 * (Ifges(3,5) * t370 + (t374 * Ifges(3,1) + t579) * t369) / 0.2e1 - t655 * t529 - t516 * t551 / 0.2e1 + t53 * (-mrSges(4,1) * t551 - t303 * mrSges(4,3)) - t294 * t389 / 0.2e1 + t303 * t632; (-t524 - t188) * t204 + (-t385 / 0.2e1 + t655 * qJD(1)) * qJD(1) + (t27 - t111) * t366 + (t467 + t642) * t291 + (-t287 + t466) * t288 + (-t365 - t187) * t205 + t737 * t108 - (t372 * t92 + t375 * t91) * t530 / 0.2e1 + (t255 * t440 + t256 * t443) * qJD(3) / 0.2e1 + (pkin(9) * t112 - t747) * t376 - ((-Ifges(3,2) * t506 + t376 * t142 + t373 * t649 + t233 + t346) * t377 + t256 * (Ifges(4,5) * t374 + t377 * t443) + t255 * (Ifges(4,6) * t374 + t377 * t440) + t374 * t140) * t538 / 0.2e1 + (-m(4) * t477 - m(5) * t416 - t625 * (t194 * pkin(4) + qJ(5) * t193 + t416) + t704 * t307 + t678 * t306 + t641 * t553 - t421 * t194 + t654 * t193) * g(1) + t577 * t619 + (-m(4) * t478 - m(5) * t417 - t625 * (t192 * pkin(4) + qJ(5) * t191 + t417) + t704 * t305 + t678 * t304 + t641 * t555 - t421 * t192 + t654 * t191) * g(2) + (-m(5) * t465 + t308 - m(4) * t542 - (t374 * mrSges(4,3) + t377 * t449) * t369 - t625 * (t258 * pkin(4) + t257 * qJ(5) + t465) - t641 * t518 - t421 * t258 + t654 * t257) * g(3) + t576 * t708 + t514 + t377 * t491 * t677 - t407 * t505 + t142 * t480 + t447 * t563 - t178 * t449 + t531 * t728 + (t244 * t743 - t245 * t744 + t463 * t651) * t609 + (-t438 * t533 + (Ifges(5,6) * t373 + t376 * t439) * qJD(3) + t701 * t463 - t738 * t245 + t699 * t244) * t618 + t694 * (t373 * t486 + t375 * t480 - t245 / 0.2e1) + t695 * (t372 * t480 + t373 * t483 - t244 / 0.2e1) + t141 * t463 / 0.2e1 + t280 * t26 + t266 * t38 + t267 * t39 + t239 * t35 + t240 * t37 + t237 * t28 + t210 * t36 + t712 * t548 + t713 * t550 - t201 * mrSges(3,2) + t202 * mrSges(3,1) + t203 * t34 - pkin(2) * t78 + (Ifges(5,4) * t245 - Ifges(5,2) * t244 + Ifges(5,6) * t463 + t663 * t533 + (t373 * t701 + t376 * t664) * qJD(3)) * t617 + (t429 * t533 + (Ifges(6,4) * t373 + t376 * t430) * qJD(3) + t700 * t463 + t703 * t245 + t702 * t244) * t616 + (Ifges(6,4) * t463 - Ifges(6,2) * t245 + Ifges(6,6) * t244 + (-t441 + t424) * t533 + (t373 * t700 + t376 * t662) * qJD(3)) * t615 - t650 * t376 / 0.2e1 + ((t436 + t435 - t431) * t533 + (t373 * t651 + t376 * t647) * qJD(3)) * t608 + (t373 * t647 - t376 * t651) * t620 + ((t538 * t547 - t531) * t150 + (t538 * t545 - t530) * t149 + t658) * mrSges(4,3) + (-t149 * t187 - t150 * t188 - pkin(2) * t178 + ((-t149 * t376 - t150 * t373) * qJD(3) + t658) * pkin(9)) * m(4) + (t3 * t444 + t430 * t630 + t439 * t629 + t8 * t445 + t92 * t484 + t91 * t485 + t628 * t664 + t631 * t662 + t632 + t639) * t373 + (mrSges(5,1) * t661 - mrSges(5,3) * t668) * t55 + (-mrSges(7,1) * t668 + mrSges(7,3) * t661) * t29 + (-mrSges(6,1) * t668 - mrSges(6,2) * t661) * t46 + (mrSges(5,2) * t661 + mrSges(5,3) * t669) * t56 + (mrSges(7,1) * t669 - mrSges(7,2) * t661) * t30 + (mrSges(6,2) * t669 + mrSges(6,3) * t668) * t54 + (-mrSges(6,1) * t669 + mrSges(6,3) * t661) * t47 + (mrSges(7,2) * t668 - mrSges(7,3) * t669) * t31 + (-mrSges(5,1) * t669 - mrSges(5,2) * t668) * t122 + t673 * t129 + t674 * t128 + (-t122 * t164 + t266 * t7 + t267 * t6 + (t122 * t530 + t563) * pkin(9) + t674 * t56 - t673 * t55) * m(5) + t531 * t621 + t244 * t626 + t245 * t627 + (t407 - t677 / 0.2e1) * qJD(3) + t683 * t109 + t684 * t127 + t685 * t125 + (t239 * t4 + t240 * t5 + t280 * t8 + t46 * t684 + t47 * t685 + t54 * t683) * m(6) + t686 * t106 + t687 * t126 + t688 * t124 + (t1 * t203 + t2 * t210 + t237 * t3 + t29 * t688 + t30 * t687 + t31 * t686) * m(7) + ((-mrSges(4,1) * t149 + mrSges(4,2) * t150) * t538 + (t232 + t693) * t491) * t374; (-t439 / 0.2e1 + t664 / 0.2e1) * qJD(4) * t199 + (t608 * t647 + t648) * qJD(4) + (-t695 / 0.2e1 + t626) * t558 + (-t694 / 0.2e1 + t627) * t557 + (t441 + t429) * t631 + (-t430 / 0.2e1 + t662 / 0.2e1) * t535 + (t247 + t142) * t607 + t8 * t446 + t516 + t92 * t486 + t91 * t484 + (t566 - t204) * t149 + (Ifges(6,4) * t615 - Ifges(4,2) * t607 + Ifges(5,6) * t617 + t651 * t609 + t700 * t616 + t701 * t618 + t732) * t256 + (-t729 + t733) * t375 + (-m(5) * t122 - t108 + t205 + t565) * t150 + (-t435 / 0.2e1 - t436 / 0.2e1) * t148 + (-t522 - t82) * t129 + t334 * t34 + t335 * t36 + t322 * t26 + (-t64 + t522) * t127 - t45 * t448 + (-m(6) * (t474 - t667) - m(5) * t474 + t450 - m(7) * (-t296 - t667) - t461 * t303 + t646 * t302) * g(3) + (-m(6) * (t475 - t671) - m(5) * t475 - m(7) * (-t215 - t671) + t689 * t228 + t726 * t227) * g(1) - t3 * t564 + (-t63 + t523) * t125 + (-m(6) * (t476 + t672) - m(5) * t476 - m(7) * (t213 + t672) + t689 * t224 - t726 * t473) * g(2) + (-t523 - t83) * t128 + t312 * t28 + t714 * t372 + (t29 * t660 - t30 * t659) * mrSges(7,1) - t52 * mrSges(4,2) + t53 * mrSges(4,1) - pkin(3) * t27 + t694 * t483 + t695 * t485 + (-pkin(3) * t45 + t55 * t82 - t56 * t83) * m(5) + (t322 * t8 - t46 * t64 - t47 * t63 + t681 * t54) * m(6) + ((t39 - t35) * t375 + (t37 - t38) * t372 + ((-t372 * t56 + t375 * t55) * qJD(4) + t657) * m(5) + ((t372 * t47 + t375 * t46) * qJD(4) + t656) * m(6)) * pkin(10) + (-t578 + t649) * t606 + (t46 * t660 + t47 * t659 + t638 + t656) * mrSges(6,1) + (t55 * t660 - t56 * t659 + t638 + t657) * mrSges(5,3) + (t438 + t663) * t629 + t141 * t605 + t431 * t620 + t679 * t124 + t680 * t126 + t681 * t109 + t682 * t106 + (t1 * t334 + t2 * t335 + t29 * t679 + t3 * t312 + t30 * t680 + t31 * t682) * m(7) + (t691 / 0.2e1 - t235 * mrSges(4,2) + Ifges(4,1) * t606 + t430 * t615 + t439 * t617 + t664 * t618 + t662 * t616 + t647 * t609 - t648) * t255 + t424 * t630; t711 + (t199 * t46 - t200 * t47) * mrSges(6,1) + t414 * t126 - t583 * t34 + (t199 * t744 + t200 * t743) * t609 + t650 + (Ifges(6,2) * t199 + t567 + t92) * t615 + (t36 - t35) * qJ(5) + (t199 * t29 + t200 * t30) * mrSges(7,1) + (qJ(5) * t2 - t1 * t583 + t675 * t29 + t30 * t717 - t31 * t77) * m(7) - t122 * (mrSges(5,1) * t200 - mrSges(5,2) * t199) - t54 * (-mrSges(6,2) * t200 + mrSges(6,3) * t199) - t31 * (mrSges(7,2) * t199 + mrSges(7,3) * t200) - t107 * t109 - t77 * t106 - pkin(4) * t37 + (t212 * t625 + t221 * t421 + t222 * t412) * g(3) + (-t127 + t129 + t581) * t56 + (-t125 + t128 + t582) * t55 + (-t199 * t703 + t195 - t568 + t695) * t616 + (t200 * t699 + t196 - t571 + t91) * t618 + (-Ifges(5,2) * t200 - t197 + t694) * t617 + t544 * qJD(5) + (t166 * t665 + t167 * t412) * g(2) + (t170 * t665 + t171 * t412) * g(1) + t675 * t124 + (-pkin(4) * t5 - qJ(5) * t4 - t107 * t54 - t46 * t56 + t47 * t676) * m(6); -t544 * t248 + (t106 + t109) * t200 + t34 + t37 + (t200 * t31 - t248 * t30 + t1 + t401) * m(7) + (t200 * t54 + t248 * t47 + t401 + t5) * m(6); -t199 * t106 + t248 * t124 + (-g(1) * t171 - g(2) * t167 - g(3) * t222 - t31 * t199 + t29 * t248 + t2) * m(7) + t36;];
tau  = t14;
