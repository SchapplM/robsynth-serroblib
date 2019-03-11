% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:35
% EndTime: 2019-03-10 01:27:52
% DurationCPUTime: 45.05s
% Computational Cost: add. (23290->952), mult. (51134->1214), div. (0->0), fcn. (36553->14), ass. (0->439)
t383 = sin(qJ(2));
t387 = cos(qJ(2));
t441 = pkin(2) * t383 - pkin(8) * t387;
t317 = t441 * qJD(1);
t386 = cos(qJ(3));
t382 = sin(qJ(3));
t494 = qJD(1) * t383;
t465 = t382 * t494;
t225 = pkin(7) * t465 + t386 * t317;
t506 = t386 * t387;
t420 = pkin(3) * t383 - pkin(9) * t506;
t389 = -pkin(9) - pkin(8);
t467 = qJD(3) * t389;
t710 = -qJD(1) * t420 + t386 * t467 - t225;
t286 = t382 * t317;
t510 = t383 * t386;
t512 = t382 * t387;
t709 = t286 + (-pkin(7) * t510 - pkin(9) * t512) * qJD(1) - t382 * t467;
t381 = sin(qJ(4));
t385 = cos(qJ(4));
t312 = t381 * t386 + t382 * t385;
t635 = qJD(3) + qJD(4);
t220 = t635 * t312;
t408 = t312 * t387;
t253 = qJD(1) * t408;
t703 = t220 - t253;
t698 = mrSges(6,1) + mrSges(7,1);
t697 = mrSges(6,2) - mrSges(7,3);
t337 = t389 * t382;
t338 = t389 * t386;
t485 = qJD(4) * t385;
t486 = qJD(4) * t381;
t650 = t337 * t485 + t338 * t486 + t381 * t710 - t709 * t385;
t224 = t381 * t337 - t385 * t338;
t649 = -qJD(4) * t224 + t709 * t381 + t385 * t710;
t708 = mrSges(4,3) + mrSges(5,3);
t668 = -Ifges(6,4) + Ifges(7,5);
t707 = t668 + Ifges(7,5);
t483 = qJD(1) * qJD(2);
t321 = qJDD(1) * t387 - t383 * t483;
t304 = t321 * pkin(7);
t706 = t304 * t387;
t422 = t381 * t382 - t385 * t386;
t219 = t635 * t422;
t407 = t422 * t387;
t254 = qJD(1) * t407;
t705 = -pkin(4) * t494 + t649 + (t219 - t254) * pkin(10);
t704 = pkin(10) * t703 - t650;
t491 = qJD(2) * t386;
t309 = -t465 + t491;
t463 = t386 * t494;
t310 = qJD(2) * t382 + t463;
t208 = t309 * t385 - t310 * t381;
t367 = pkin(7) * t494;
t334 = -qJD(2) * pkin(2) + t367;
t230 = -t309 * pkin(3) + t334;
t158 = -pkin(4) * t208 + t230;
t493 = qJD(1) * t387;
t352 = qJD(3) - t493;
t343 = qJD(4) + t352;
t330 = qJD(5) + t343;
t380 = sin(qJ(5));
t442 = pkin(2) * t387 + pkin(8) * t383;
t329 = -pkin(1) - t442;
t295 = t329 * qJD(1);
t368 = pkin(7) * t493;
t335 = qJD(2) * pkin(8) + t368;
t214 = t386 * t295 - t335 * t382;
t172 = -pkin(9) * t310 + t214;
t160 = pkin(3) * t352 + t172;
t215 = t295 * t382 + t335 * t386;
t173 = pkin(9) * t309 + t215;
t167 = t385 * t173;
t103 = t381 * t160 + t167;
t554 = pkin(10) * t208;
t86 = t103 + t554;
t528 = t380 * t86;
t568 = cos(qJ(5));
t165 = t381 * t173;
t102 = t385 * t160 - t165;
t209 = t309 * t381 + t310 * t385;
t690 = pkin(10) * t209;
t85 = t102 - t690;
t78 = pkin(4) * t343 + t85;
t39 = t568 * t78 - t528;
t687 = -t39 + qJD(6);
t30 = -t330 * pkin(5) + t687;
t140 = -t568 * t208 + t380 * t209;
t680 = t208 * t380 + t209 * t568;
t61 = t140 * pkin(5) - qJ(6) * t680 + t158;
t137 = Ifges(6,4) * t140;
t533 = Ifges(7,5) * t140;
t667 = Ifges(7,4) + Ifges(6,5);
t669 = Ifges(6,1) + Ifges(7,1);
t658 = t330 * t667 + t669 * t680 - t137 + t533;
t702 = mrSges(6,2) * t158 - mrSges(7,3) * t61 + mrSges(7,2) * t30 - mrSges(6,3) * t39 + t658 / 0.2e1;
t79 = pkin(5) * t680 + qJ(6) * t140;
t306 = qJDD(3) - t321;
t294 = qJDD(4) + t306;
t277 = qJDD(5) + t294;
t322 = qJDD(1) * t383 + t387 * t483;
t523 = qJDD(1) * pkin(1);
t218 = -pkin(2) * t321 - pkin(8) * t322 - t523;
t272 = qJDD(2) * pkin(8) + t304;
t121 = -qJD(3) * t215 + t386 * t218 - t272 * t382;
t199 = qJD(3) * t309 + qJDD(2) * t382 + t322 * t386;
t87 = pkin(3) * t306 - pkin(9) * t199 + t121;
t487 = qJD(3) * t386;
t489 = qJD(3) * t382;
t120 = t382 * t218 + t386 * t272 + t295 * t487 - t335 * t489;
t200 = -qJD(3) * t310 + qJDD(2) * t386 - t322 * t382;
t94 = pkin(9) * t200 + t120;
t25 = -qJD(4) * t103 - t381 * t94 + t385 * t87;
t99 = qJD(4) * t208 + t199 * t385 + t200 * t381;
t19 = pkin(4) * t294 - pkin(10) * t99 + t25;
t100 = -qJD(4) * t209 - t199 * t381 + t200 * t385;
t24 = t160 * t485 - t173 * t486 + t381 * t87 + t385 * t94;
t21 = pkin(10) * t100 + t24;
t457 = qJD(5) * t568;
t484 = qJD(5) * t380;
t5 = t380 * t19 + t568 * t21 + t78 * t457 - t484 * t86;
t2 = qJ(6) * t277 + qJD(6) * t330 + t5;
t471 = t568 * t86;
t40 = t380 * t78 + t471;
t6 = -qJD(5) * t40 + t19 * t568 - t380 * t21;
t3 = -t277 * pkin(5) + qJDD(6) - t6;
t623 = -t6 * mrSges(6,1) + t3 * mrSges(7,1) + t5 * mrSges(6,2) - t2 * mrSges(7,3);
t37 = -qJD(5) * t140 + t380 * t100 + t568 * t99;
t38 = qJD(5) * t680 - t100 * t568 + t380 * t99;
t665 = Ifges(6,3) + Ifges(7,2);
t666 = -Ifges(6,6) + Ifges(7,6);
t642 = t277 * t665 + t37 * t667 + t38 * t666;
t396 = -t623 + t642;
t477 = Ifges(5,5) * t99 + Ifges(5,6) * t100 + Ifges(5,3) * t294;
t203 = Ifges(5,4) * t208;
t132 = t209 * Ifges(5,1) + t343 * Ifges(5,5) + t203;
t593 = -t132 / 0.2e1;
t573 = -t330 / 0.2e1;
t587 = -t680 / 0.2e1;
t590 = t140 / 0.2e1;
t591 = -t140 / 0.2e1;
t31 = t330 * qJ(6) + t40;
t136 = Ifges(7,5) * t680;
t71 = t330 * Ifges(7,6) + t140 * Ifges(7,3) + t136;
t535 = Ifges(6,4) * t680;
t74 = -t140 * Ifges(6,2) + t330 * Ifges(6,6) + t535;
t694 = -mrSges(7,2) * t31 - mrSges(6,3) * t40 + mrSges(6,1) * t158 + mrSges(7,1) * t61 + t71 / 0.2e1 - t74 / 0.2e1;
t615 = -Ifges(6,2) * t590 + Ifges(7,3) * t591 + t573 * t666 + t587 * t668 - t694;
t628 = -t25 * mrSges(5,1) + t24 * mrSges(5,2);
t683 = Ifges(6,4) * t590 + Ifges(7,5) * t591 + t573 * t667 + t587 * t669 - t702;
t701 = -t683 * t140 + t208 * t593 + t615 * t680 + t396 + t477 - t628;
t609 = t37 / 0.2e1;
t607 = t38 / 0.2e1;
t700 = m(6) + m(7);
t578 = t277 / 0.2e1;
t696 = -mrSges(6,3) - mrSges(7,2);
t464 = t382 * t493;
t646 = -t368 + (-t464 + t489) * pkin(3);
t695 = -m(4) * pkin(8) + m(5) * t389 - t708;
t379 = qJ(3) + qJ(4);
t373 = qJ(5) + t379;
t359 = sin(t373);
t360 = cos(t373);
t559 = pkin(3) * t386;
t363 = pkin(2) + t559;
t371 = sin(t379);
t372 = cos(t379);
t438 = -mrSges(4,1) * t386 + mrSges(4,2) * t382;
t692 = m(4) * pkin(2) + m(5) * t363 + t372 * mrSges(5,1) - t371 * mrSges(5,2) - t697 * t359 + t360 * t698 - t438;
t558 = pkin(4) * t209;
t223 = t385 * t337 + t338 * t381;
t181 = -pkin(10) * t312 + t223;
t182 = -pkin(10) * t422 + t224;
t414 = t181 * t568 - t380 * t182;
t661 = qJD(5) * t414 + t380 * t705 - t568 * t704;
t127 = t380 * t181 + t182 * t568;
t660 = -qJD(5) * t127 + t380 * t704 + t568 * t705;
t689 = m(7) * qJ(6) + mrSges(7,3);
t560 = pkin(3) * t385;
t362 = pkin(4) + t560;
t109 = -t172 * t381 - t167;
t419 = t109 - t554;
t475 = pkin(3) * t485;
t562 = pkin(3) * t381;
t480 = t380 * t562;
t110 = t385 * t172 - t165;
t90 = t110 - t690;
t655 = -t380 * t419 + t362 * t457 + (-qJD(4) - qJD(5)) * t480 + (-t90 + t475) * t568;
t542 = mrSges(6,3) * t680;
t124 = mrSges(6,1) * t330 - t542;
t125 = -mrSges(7,1) * t330 + mrSges(7,2) * t680;
t688 = t124 - t125;
t648 = pkin(4) * t703 + t646;
t490 = qJD(2) * t387;
t402 = t382 * t490 + t383 * t487;
t384 = sin(qJ(1));
t388 = cos(qJ(1));
t505 = t387 * t388;
t260 = -t371 * t505 + t372 * t384;
t608 = -t38 / 0.2e1;
t305 = t322 * pkin(7);
t273 = -qJDD(2) * pkin(2) + t305;
t161 = -pkin(3) * t200 + t273;
t70 = -pkin(4) * t100 + t161;
t7 = pkin(5) * t38 - qJ(6) * t37 - qJD(6) * t680 + t70;
t684 = mrSges(6,2) * t70 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t7 + Ifges(6,4) * t608 + 0.2e1 * t578 * t667 + t607 * t707 + 0.2e1 * t609 * t669;
t572 = t330 / 0.2e1;
t586 = t680 / 0.2e1;
t682 = -Ifges(6,2) * t591 + Ifges(7,3) * t590 + t572 * t666 + t586 * t668 + t694;
t678 = Ifges(6,4) * t591 + Ifges(7,5) * t590 + t572 * t667 + t586 * t669 + t702;
t612 = m(5) * pkin(3);
t599 = t99 / 0.2e1;
t598 = t100 / 0.2e1;
t583 = t199 / 0.2e1;
t582 = t200 / 0.2e1;
t676 = -t208 / 0.2e1;
t675 = t208 / 0.2e1;
t577 = t294 / 0.2e1;
t576 = t306 / 0.2e1;
t674 = t321 / 0.2e1;
t673 = t322 / 0.2e1;
t370 = pkin(7) * t490;
t670 = mrSges(2,2) - mrSges(3,3);
t412 = -t380 * t312 - t422 * t568;
t106 = qJD(5) * t412 - t219 * t568 - t380 * t220;
t211 = t312 * t568 - t380 * t422;
t107 = qJD(5) * t211 - t380 * t219 + t220 * t568;
t168 = t253 * t568 - t254 * t380;
t169 = -t380 * t253 - t254 * t568;
t663 = -qJD(6) * t211 + t648 + (-t106 + t169) * qJ(6) + (t107 - t168) * pkin(5);
t662 = -qJ(6) * t494 + t661;
t659 = pkin(5) * t494 - t660;
t657 = t208 * Ifges(5,6);
t656 = qJD(6) + t655;
t652 = t612 + mrSges(4,1);
t267 = t422 * t383;
t308 = t386 * t329;
t213 = -pkin(9) * t510 + t308 + (-pkin(7) * t382 - pkin(3)) * t387;
t355 = pkin(7) * t506;
t248 = t382 * t329 + t355;
t514 = t382 * t383;
t222 = -pkin(9) * t514 + t248;
t148 = t385 * t213 - t222 * t381;
t115 = -pkin(4) * t387 + pkin(10) * t267 + t148;
t149 = t381 * t213 + t385 * t222;
t266 = t312 * t383;
t128 = -pkin(10) * t266 + t149;
t651 = t380 * t115 + t568 * t128;
t299 = Ifges(4,4) * t309;
t185 = t310 * Ifges(4,1) + t352 * Ifges(4,5) + t299;
t366 = Ifges(3,4) * t493;
t647 = Ifges(3,1) * t494 + Ifges(3,5) * qJD(2) + t386 * t185 + t366;
t645 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t309 - mrSges(4,2) * t310 - mrSges(3,3) * t494;
t547 = mrSges(6,2) * t360;
t643 = mrSges(5,1) * t371 + mrSges(6,1) * t359 + mrSges(5,2) * t372 + t547;
t504 = t388 * t359;
t244 = -t384 * t360 + t387 * t504;
t245 = t359 * t384 + t360 * t505;
t641 = t244 * t698 + t245 * t697;
t508 = t384 * t387;
t242 = t359 * t508 + t360 * t388;
t243 = t360 * t508 - t504;
t640 = t242 * t698 + t243 * t697;
t639 = t305 * t383 + t706;
t638 = t120 * t386 - t121 * t382;
t637 = t696 * t383;
t636 = -m(5) - m(4) - m(3);
t261 = t371 * t384 + t372 * t505;
t634 = -t260 * mrSges(5,1) + t261 * mrSges(5,2) + t641;
t258 = t371 * t508 + t372 * t388;
t259 = t371 * t388 - t372 * t508;
t633 = t258 * mrSges(5,1) - t259 * mrSges(5,2) + t640;
t632 = t310 * Ifges(4,5) + Ifges(5,5) * t209 + t309 * Ifges(4,6) + t352 * Ifges(4,3) + Ifges(5,3) * t343 + t140 * t666 + t330 * t665 + t667 * t680 + t657;
t153 = -qJD(2) * t407 - t220 * t383;
t492 = qJD(2) * t383;
t320 = t441 * qJD(2);
t473 = pkin(7) * t492;
t496 = t386 * t320 + t382 * t473;
t147 = t420 * qJD(2) + (-t355 + (pkin(9) * t383 - t329) * t382) * qJD(3) + t496;
t170 = t382 * t320 + t329 * t487 + (-t383 * t491 - t387 * t489) * pkin(7);
t152 = -pkin(9) * t402 + t170;
t60 = -qJD(4) * t149 + t385 * t147 - t152 * t381;
t47 = pkin(4) * t492 - pkin(10) * t153 + t60;
t154 = -qJD(2) * t408 + t267 * t635;
t59 = t381 * t147 + t385 * t152 + t213 * t485 - t222 * t486;
t51 = pkin(10) * t154 + t59;
t11 = -qJD(5) * t651 - t380 * t51 + t47 * t568;
t631 = m(7) * pkin(5) + t698;
t440 = mrSges(3,1) * t387 - mrSges(3,2) * t383;
t630 = t383 * t708 + mrSges(2,1) + t440;
t629 = -mrSges(6,2) + t689;
t625 = -t121 * mrSges(4,1) + t120 * mrSges(4,2);
t540 = Ifges(3,4) * t383;
t430 = t387 * Ifges(3,2) + t540;
t620 = t103 * mrSges(5,2) + t30 * mrSges(7,1) + t40 * mrSges(6,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t430 / 0.2e1 - t102 * mrSges(5,1) - t657 / 0.2e1 - t31 * mrSges(7,3) - t39 * mrSges(6,1);
t614 = mrSges(6,1) * t70 + mrSges(7,1) * t7 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t607 - t37 * Ifges(6,4) / 0.2e1 - t277 * Ifges(6,6) / 0.2e1 + t707 * t609 + (t666 + Ifges(7,6)) * t578 + (-t608 + t607) * Ifges(6,2);
t606 = Ifges(5,4) * t599 + Ifges(5,2) * t598 + Ifges(5,6) * t577;
t605 = Ifges(5,1) * t599 + Ifges(5,4) * t598 + Ifges(5,5) * t577;
t596 = Ifges(4,1) * t583 + Ifges(4,4) * t582 + Ifges(4,5) * t576;
t536 = Ifges(5,4) * t209;
t131 = t208 * Ifges(5,2) + Ifges(5,6) * t343 + t536;
t595 = -t131 / 0.2e1;
t594 = t131 / 0.2e1;
t592 = t132 / 0.2e1;
t581 = -t209 / 0.2e1;
t580 = t209 / 0.2e1;
t574 = t310 / 0.2e1;
t571 = -t343 / 0.2e1;
t570 = t343 / 0.2e1;
t563 = pkin(3) * t310;
t561 = pkin(3) * t382;
t557 = pkin(4) * t371;
t556 = pkin(4) * t380;
t555 = pkin(5) * t359;
t551 = g(3) * t383;
t374 = t383 * pkin(7);
t548 = mrSges(5,2) * t208;
t546 = mrSges(5,3) * t102;
t545 = mrSges(5,3) * t103;
t544 = mrSges(5,3) * t209;
t543 = mrSges(6,3) * t140;
t541 = Ifges(5,1) * t208;
t539 = Ifges(3,4) * t387;
t538 = Ifges(4,4) * t382;
t537 = Ifges(4,4) * t386;
t534 = Ifges(5,5) * t208;
t532 = t214 * mrSges(4,3);
t531 = t215 * mrSges(4,3);
t530 = t310 * Ifges(4,4);
t529 = t359 * mrSges(7,1);
t378 = -pkin(10) + t389;
t515 = t378 * t383;
t513 = t382 * t384;
t511 = t382 * t388;
t509 = t383 * t388;
t325 = pkin(4) * t372 + t559;
t316 = pkin(2) + t325;
t292 = t387 * t316;
t324 = t557 + t561;
t302 = t388 * t324;
t497 = -t387 * t302 + t384 * t325;
t466 = t568 * t381;
t275 = pkin(3) * t466 + t380 * t362;
t323 = pkin(3) * t514 + t374;
t495 = t388 * pkin(1) + t384 * pkin(7);
t488 = qJD(3) * t383;
t476 = t568 * pkin(4);
t474 = pkin(4) * t484;
t469 = Ifges(4,5) * t199 + Ifges(4,6) * t200 + Ifges(4,3) * t306;
t241 = pkin(3) * t402 + t370;
t184 = t309 * Ifges(4,2) + t352 * Ifges(4,6) + t530;
t459 = -t382 * t184 / 0.2e1;
t450 = t689 * t360 * t383;
t28 = -t277 * mrSges(7,1) + t37 * mrSges(7,2);
t448 = t483 / 0.2e1;
t447 = -t242 * pkin(5) + qJ(6) * t243;
t446 = -t244 * pkin(5) + qJ(6) * t245;
t445 = -t324 * t508 - t325 * t388;
t443 = pkin(4) * t457;
t216 = pkin(4) * t266 + t323;
t164 = t563 + t558;
t439 = mrSges(3,1) * t383 + mrSges(3,2) * t387;
t437 = mrSges(4,1) * t382 + mrSges(4,2) * t386;
t432 = Ifges(4,1) * t386 - t538;
t431 = Ifges(4,1) * t382 + t537;
t429 = -Ifges(4,2) * t382 + t537;
t428 = Ifges(4,2) * t386 + t538;
t427 = Ifges(3,5) * t387 - Ifges(3,6) * t383;
t426 = Ifges(4,5) * t386 - Ifges(4,6) * t382;
t425 = Ifges(4,5) * t382 + Ifges(4,6) * t386;
t424 = pkin(5) * t360 + qJ(6) * t359;
t423 = t387 * t363 - t383 * t389;
t421 = t260 * pkin(4);
t255 = pkin(4) * t422 - t363;
t129 = -pkin(4) * t154 + t241;
t416 = pkin(1) * t439;
t283 = -t382 * t505 + t384 * t386;
t281 = t382 * t508 + t386 * t388;
t66 = t115 * t568 - t380 * t128;
t413 = -t266 * t568 + t380 * t267;
t177 = -t380 * t266 - t267 * t568;
t411 = t334 * t437;
t410 = t383 * (Ifges(3,1) * t387 - t540);
t10 = t115 * t457 - t128 * t484 + t380 * t47 + t568 * t51;
t409 = t208 * mrSges(5,3);
t274 = t362 * t568 - t480;
t404 = t258 * pkin(4);
t403 = -t382 * t488 + t386 * t490;
t400 = Ifges(4,5) * t383 + t387 * t432;
t399 = Ifges(4,6) * t383 + t387 * t429;
t398 = Ifges(4,3) * t383 + t387 * t426;
t376 = t388 * pkin(7);
t361 = -t476 - pkin(5);
t356 = qJ(6) + t556;
t350 = t443 + qJD(6);
t332 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t493;
t296 = t437 * t383;
t284 = t386 * t505 + t513;
t282 = -t384 * t506 + t511;
t269 = -pkin(5) - t274;
t265 = qJ(6) + t275;
t247 = -pkin(7) * t512 + t308;
t229 = mrSges(4,1) * t352 - mrSges(4,3) * t310;
t228 = -mrSges(4,2) * t352 + mrSges(4,3) * t309;
t226 = -pkin(7) * t463 + t286;
t175 = mrSges(5,1) * t343 - t544;
t174 = -t343 * mrSges(5,2) + t409;
t171 = -qJD(3) * t248 + t496;
t163 = -mrSges(4,2) * t306 + mrSges(4,3) * t200;
t162 = mrSges(4,1) * t306 - mrSges(4,3) * t199;
t146 = -mrSges(5,1) * t208 + t209 * mrSges(5,2);
t135 = -mrSges(4,1) * t200 + mrSges(4,2) * t199;
t123 = -mrSges(6,2) * t330 - t543;
t122 = -mrSges(7,2) * t140 + mrSges(7,3) * t330;
t119 = -pkin(5) * t412 - qJ(6) * t211 + t255;
t116 = t199 * Ifges(4,4) + t200 * Ifges(4,2) + t306 * Ifges(4,6);
t95 = -pkin(5) * t413 - qJ(6) * t177 + t216;
t92 = -mrSges(5,2) * t294 + mrSges(5,3) * t100;
t91 = mrSges(5,1) * t294 - mrSges(5,3) * t99;
t81 = mrSges(6,1) * t140 + mrSges(6,2) * t680;
t80 = mrSges(7,1) * t140 - mrSges(7,3) * t680;
t69 = qJD(5) * t177 + t380 * t153 - t154 * t568;
t68 = qJD(5) * t413 + t153 * t568 + t380 * t154;
t65 = t558 + t79;
t64 = t387 * pkin(5) - t66;
t63 = -qJ(6) * t387 + t651;
t62 = t164 + t79;
t54 = -mrSges(5,1) * t100 + mrSges(5,2) * t99;
t42 = t568 * t85 - t528;
t41 = t380 * t85 + t471;
t29 = -mrSges(6,2) * t277 - mrSges(6,3) * t38;
t27 = mrSges(6,1) * t277 - mrSges(6,3) * t37;
t26 = -mrSges(7,2) * t38 + mrSges(7,3) * t277;
t22 = pkin(5) * t69 - qJ(6) * t68 - qJD(6) * t177 + t129;
t17 = mrSges(6,1) * t38 + mrSges(6,2) * t37;
t16 = mrSges(7,1) * t38 - mrSges(7,3) * t37;
t9 = -pkin(5) * t492 - t11;
t8 = qJ(6) * t492 - qJD(6) * t387 + t10;
t1 = [m(4) * (t120 * t248 + t121 * t247 + t170 * t215 + t171 * t214 + t334 * t370) + m(6) * (t10 * t40 + t11 * t39 + t129 * t158 + t216 * t70 + t5 * t651 + t6 * t66) + t651 * t29 - t614 * t413 + t440 * t523 + t678 * t68 + (t214 * mrSges(4,1) - t215 * mrSges(4,2) + Ifges(5,5) * t580 + Ifges(6,6) * t591 + Ifges(7,6) * t590 + Ifges(5,3) * t570 + t665 * t572 + t667 * t586 - t620 + t632 / 0.2e1) * t492 + (-mrSges(3,1) * t374 + Ifges(3,5) * t383 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * t387) * qJDD(2) - (t184 * t386 + t185 * t382) * t488 / 0.2e1 + t539 * t673 + t430 * t674 + (Ifges(5,4) * t153 + Ifges(5,2) * t154) * t675 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t639) - (t477 + t469 + t642) * t387 / 0.2e1 - t645 * t370 + (t647 / 0.2e1 + t459) * t490 + t410 * t448 + (-t120 * t514 - t121 * t510 - t214 * t403 - t215 * t402) * mrSges(4,3) + t334 * (mrSges(4,1) * t402 + mrSges(4,2) * t403) + qJD(2) ^ 2 * t427 / 0.2e1 + t682 * t69 + (Ifges(5,1) * t153 + Ifges(5,4) * t154) * t580 + m(5) * (t102 * t60 + t103 * t59 + t148 * t25 + t149 * t24 + t161 * t323 + t230 * t241) + m(7) * (t2 * t63 + t22 * t61 + t3 * t64 + t30 * t9 + t31 * t8 + t7 * t95) + (-Ifges(5,4) * t267 - Ifges(5,2) * t266) * t598 + (-t102 * t153 + t103 * t154 - t24 * t266 + t25 * t267) * mrSges(5,3) + (-Ifges(5,1) * t267 - Ifges(5,4) * t266) * t599 + t161 * (mrSges(5,1) * t266 - mrSges(5,2) * t267) + (-Ifges(5,5) * t267 - Ifges(5,6) * t266) * t577 + t684 * t177 + (-t511 * t612 - t282 * mrSges(4,1) - t259 * mrSges(5,1) - t281 * mrSges(4,2) - t258 * mrSges(5,2) - t700 * (t384 * t515 + t302 + t376) + t670 * t388 + t636 * t376 + t631 * t243 + t629 * t242 + (m(3) * pkin(1) - m(5) * (-pkin(1) - t423) - m(4) * t329 - t700 * (-pkin(1) - t292) + t630 - t637) * t384) * g(1) + (-t513 * t612 - t284 * mrSges(4,1) - t261 * mrSges(5,1) - t283 * mrSges(4,2) - t260 * mrSges(5,2) + t696 * t509 + t636 * t495 - t700 * (t316 * t505 + t384 * t324 - t378 * t509 + t495) + t670 * t384 - t631 * t245 - t629 * t244 + (-m(4) * t442 - m(5) * t423 - t630) * t388) * g(2) + (t322 * t374 + t639 + t706) * mrSges(3,3) + Ifges(2,3) * qJDD(1) + t510 * t596 - t267 * t605 - t266 * t606 + t153 * t592 + t154 * t594 - pkin(1) * (-mrSges(3,1) * t321 + mrSges(3,2) * t322) + t323 * t54 + t273 * t296 + t247 * t162 + t248 * t163 + t241 * t146 + t170 * t228 + t171 * t229 + t230 * (-mrSges(5,1) * t154 + mrSges(5,2) * t153) + t216 * t17 - t116 * t514 / 0.2e1 + t59 * t174 + t60 * t175 + (-Ifges(5,6) * t598 - Ifges(5,5) * t599 - Ifges(7,6) * t607 - Ifges(6,6) * t608 + (-Ifges(3,2) * t383 + t539) * t448 - Ifges(4,3) * t576 - Ifges(5,3) * t577 - Ifges(4,6) * t582 - Ifges(4,5) * t583 + t623 - t665 * t578 - t667 * t609 + Ifges(3,2) * t674 + Ifges(3,4) * t673 + t625 + t628) * t387 + (m(4) * t273 * pkin(7) + Ifges(3,1) * t322 + Ifges(3,4) * t674 + t426 * t576 + t429 * t582 + t432 * t583) * t383 + t63 * t26 + t64 * t28 + t66 * t27 + t22 * t80 + t95 * t16 - t332 * t473 + t135 * t374 + (qJD(2) * t400 - t431 * t488) * t574 + (Ifges(5,5) * t153 + Ifges(5,6) * t154) * t570 + t8 * t122 + t10 * t123 + t11 * t124 + t9 * t125 + t129 * t81 - t416 * t483 + t309 * (qJD(2) * t399 - t428 * t488) / 0.2e1 + t352 * (qJD(2) * t398 - t425 * t488) / 0.2e1 + t148 * t91 + t149 * t92; (Ifges(5,5) * t581 + Ifges(6,6) * t590 + Ifges(7,6) * t591 + Ifges(5,3) * t571 + t573 * t665 + t587 * t667 + t620) * t494 + (-mrSges(5,2) * t230 - Ifges(5,1) * t580 - Ifges(5,4) * t675 - Ifges(5,5) * t570 + t546 - t592) * t219 + t683 * t169 - (t28 - t27) * t414 + (t127 * t5 + t158 * t648 + t255 * t70 + t39 * t660 + t40 * t661 + t414 * t6) * m(6) + (t119 * t7 + t127 * t2 - t3 * t414 + t30 * t659 + t31 * t662 + t61 * t663) * m(7) - t614 * t412 + t161 * (mrSges(5,1) * t422 + mrSges(5,2) * t312) + (Ifges(5,5) * t312 - Ifges(5,6) * t422) * t577 + (Ifges(5,4) * t312 - Ifges(5,2) * t422) * t598 + (Ifges(5,1) * t312 - Ifges(5,4) * t422) * t599 - t422 * t606 + (-t532 + t185 / 0.2e1) * t487 + (t26 + t29) * t127 - t632 * t494 / 0.2e1 + t678 * t106 + (t309 * t429 + t310 * t432 + t352 * t426) * qJD(3) / 0.2e1 - (t309 * t399 + t310 * t400 + t352 * t398) * qJD(1) / 0.2e1 + (-t410 / 0.2e1 + t416) * qJD(1) ^ 2 + (t386 * t163 - t382 * t162 - t228 * t489 - t229 * t487 + m(4) * ((-t214 * t386 - t215 * t382) * qJD(3) + t638)) * pkin(8) + t638 * mrSges(4,3) + t645 * t368 + t646 * t146 - (-Ifges(3,2) * t494 + t366 + t647) * t493 / 0.2e1 + t648 * t81 + t649 * t175 + t650 * t174 + (t102 * t649 + t103 * t650 - t161 * t363 + t223 * t25 + t224 * t24 + t230 * t646) * m(5) + (-Ifges(5,4) * t580 - Ifges(5,2) * t675 - Ifges(5,6) * t570 - t545 - t594) * t220 + t332 * t367 + t184 * t464 / 0.2e1 + t615 * t168 + t682 * t107 + (-pkin(2) * t273 - t214 * t225 - t215 * t226 - t334 * t368) * m(4) + (t459 + t411) * qJD(3) + t273 * t438 + (-t440 - t700 * (t292 - t515) + t695 * t383 + (-m(7) * t424 - t692) * t387 + t637) * g(3) + t684 * t211 + (-t215 * (-mrSges(4,2) * t383 - mrSges(4,3) * t512) - t214 * (mrSges(4,1) * t383 - mrSges(4,3) * t506)) * qJD(1) + (t703 * mrSges(5,1) + mrSges(5,2) * t254) * t230 + t382 * t596 + t312 * t605 - t254 * t593 - t253 * t595 + t386 * t116 / 0.2e1 + Ifges(3,3) * qJDD(2) - t363 * t54 + Ifges(3,6) * t321 + Ifges(3,5) * t322 - t304 * mrSges(3,2) - t305 * mrSges(3,1) - t489 * t531 + t255 * t17 - t226 * t228 - t225 * t229 + t223 * t91 + t224 * t92 + (-t102 * t254 + t103 * t253 - t24 * t422 - t25 * t312) * mrSges(5,3) + (-Ifges(5,4) * t254 - Ifges(5,2) * t253) * t676 + (-Ifges(5,5) * t254 - Ifges(5,6) * t253) * t571 + (-Ifges(5,1) * t254 - Ifges(5,4) * t253) * t581 + t659 * t125 + t660 * t124 + t661 * t123 + t662 * t122 + t663 * t80 + (g(1) * t388 + g(2) * t384) * (t439 + (t378 * t700 + t695 + t696) * t387 + (m(6) * t316 - m(7) * (-t316 - t424) + t692) * t383) + t425 * t576 + t428 * t582 + t431 * t583 + t119 * t16 - pkin(2) * t135 - t427 * t483 / 0.2e1 - t411 * t493; t701 + (-mrSges(5,1) * t230 - Ifges(5,4) * t581 - Ifges(5,2) * t676 - Ifges(5,6) * t571 + t545 - t595) * t209 - t625 + t655 * t123 + t656 * t122 - (-Ifges(4,2) * t310 + t185 + t299) * t309 / 0.2e1 + (-mrSges(4,2) * t282 - m(7) * (t445 + t447) - m(6) * t445 + t652 * t281 + t633) * g(2) + (mrSges(4,2) * t284 - m(7) * (t446 + t497) - m(6) * t497 - t652 * t283 + t634) * g(1) + (m(5) * t561 + t643) * t551 + t469 + (-t110 + t475) * t174 + (-pkin(3) * t486 - t109) * t175 + (-t158 * t164 + t274 * t6 + t275 * t5 + t324 * t551 + t40 * t655) * m(6) + (-m(6) * t39 + m(7) * t30 - t688) * (-t380 * t90 + t419 * t568 + t362 * t484 + (t381 * t457 + (t380 * t385 + t466) * qJD(4)) * pkin(3)) + (-(m(7) * (-t324 - t555) - t529) * t383 - t450 + t296) * g(3) + t541 * t581 + t203 * t676 - t230 * t548 + t534 * t571 + (t24 * t381 + t25 * t385 + (-t102 * t381 + t103 * t385) * qJD(4)) * t612 - t146 * t563 - m(5) * (t102 * t109 + t103 * t110 + t230 * t563) - t352 * (Ifges(4,5) * t309 - Ifges(4,6) * t310) / 0.2e1 - t334 * (mrSges(4,1) * t310 + mrSges(4,2) * t309) + t269 * t28 + t274 * t27 + t275 * t29 - t310 * (Ifges(4,1) * t309 - t530) / 0.2e1 + t265 * t26 - t214 * t228 + t215 * t229 - t164 * t81 + (t2 * t265 + t269 * t3 + t31 * t656 - t61 * t62) * m(7) - t62 * t80 + t310 * t531 + t309 * t532 + t208 * t546 + t91 * t560 + t92 * t562 + t184 * t574; t701 + (-m(7) * (-t404 + t447) + m(6) * t404 + t633) * g(2) + (-m(7) * (t421 + t446) - m(6) * t421 + t634) * g(1) + (-Ifges(5,2) * t209 + t203) * t676 + t643 * t551 + (t409 - t174) * t102 + (t175 + t544) * t103 + t688 * (-t474 + t41) - m(7) * (t30 * t41 + t31 * t42 + t61 * t65) + t27 * t476 + (t350 - t42) * t122 + ((t568 * t6 + t380 * t5 + (-t380 * t39 + t40 * t568) * qJD(5)) * pkin(4) + t557 * t551 - t158 * t558 + t39 * t41 - t40 * t42) * m(6) + (t443 - t42) * t123 - g(3) * ((m(7) * (-t555 - t557) - t529) * t383 + t450) - t81 * t558 + t361 * t28 + t356 * t26 - t230 * (t209 * mrSges(5,1) + t548) - t65 * t80 + m(7) * (t2 * t356 + t3 * t361 + t30 * t474 + t31 * t350) + t29 * t556 + (-Ifges(5,6) * t209 + t534) * t571 + t131 * t580 + (-t536 + t541) * t581; t396 + qJ(6) * t26 - pkin(5) * t28 + (Ifges(7,3) * t680 - t533) * t591 + ((t359 * t631 + t547) * t383 - t450) * g(3) + (t542 + t688) * t40 + (-t122 - t123 - t543) * t39 + t641 * g(1) + t640 * g(2) - t158 * (mrSges(6,1) * t680 - mrSges(6,2) * t140) - t79 * t80 + t74 * t586 + qJD(6) * t122 + (t140 * t30 + t31 * t680) * mrSges(7,2) - t61 * (mrSges(7,1) * t680 + mrSges(7,3) * t140) + (-t140 * t667 + t666 * t680) * t573 + (-Ifges(6,2) * t680 - t137 + t658) * t590 + (-t140 * t669 + t136 - t535 + t71) * t587 + (-pkin(5) * t3 - t446 * g(1) - t447 * g(2) + qJ(6) * t2 - t30 * t40 + t31 * t687 - t61 * t79) * m(7); -t330 * t122 + t680 * t80 + (-g(1) * t244 - g(2) * t242 - t31 * t330 - t359 * t551 + t61 * t680 + t3) * m(7) + t28;];
tau  = t1;
