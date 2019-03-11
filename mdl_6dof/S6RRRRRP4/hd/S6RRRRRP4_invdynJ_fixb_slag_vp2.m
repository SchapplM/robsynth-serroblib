% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP4
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:07
% EndTime: 2019-03-10 01:12:57
% DurationCPUTime: 28.52s
% Computational Cost: add. (22274->887), mult. (47721->1116), div. (0->0), fcn. (34663->14), ass. (0->424)
t708 = Ifges(6,1) + Ifges(7,1);
t706 = Ifges(7,4) + Ifges(6,5);
t707 = -Ifges(6,4) + Ifges(7,5);
t764 = t707 + Ifges(7,5);
t431 = sin(qJ(4));
t535 = qJD(4) * t431;
t432 = sin(qJ(3));
t437 = cos(qJ(3));
t438 = cos(qJ(2));
t539 = qJD(1) * t438;
t433 = sin(qJ(2));
t540 = qJD(1) * t433;
t332 = -t432 * t540 + t437 * t539;
t578 = t332 * t431;
t763 = t535 - t578;
t351 = t432 * t438 + t433 * t437;
t333 = t351 * qJD(1);
t427 = qJD(2) + qJD(3);
t436 = cos(qJ(4));
t288 = t333 * t436 + t427 * t431;
t430 = sin(qJ(5));
t435 = cos(qJ(5));
t472 = t333 * t431 - t427 * t436;
t450 = t435 * t288 - t430 * t472;
t199 = Ifges(7,5) * t450;
t203 = t430 * t288 + t435 * t472;
t327 = qJD(4) - t332;
t315 = qJD(5) + t327;
t101 = Ifges(7,6) * t315 + Ifges(7,3) * t203 + t199;
t589 = Ifges(6,4) * t450;
t104 = -Ifges(6,2) * t203 + Ifges(6,6) * t315 + t589;
t441 = -pkin(8) - pkin(7);
t381 = t441 * t438;
t355 = qJD(1) * t381;
t334 = t432 * t355;
t379 = t441 * t433;
t354 = qJD(1) * t379;
t340 = qJD(2) * pkin(2) + t354;
t269 = t437 * t340 + t334;
t247 = -t427 * pkin(3) - t269;
t191 = pkin(4) * t472 + t247;
t425 = t438 * pkin(2);
t414 = t425 + pkin(1);
t377 = t414 * qJD(1);
t238 = -pkin(3) * t332 - pkin(9) * t333 - t377;
t552 = t437 * t355;
t270 = t432 * t340 - t552;
t248 = pkin(9) * t427 + t270;
t159 = t436 * t238 - t248 * t431;
t127 = -pkin(10) * t288 + t159;
t111 = pkin(4) * t327 + t127;
t160 = t431 * t238 + t436 * t248;
t128 = -pkin(10) * t472 + t160;
t554 = t435 * t128;
t48 = t111 * t430 + t554;
t42 = qJ(6) * t315 + t48;
t79 = t203 * pkin(5) - qJ(6) * t450 + t191;
t762 = -mrSges(7,2) * t42 - mrSges(6,3) * t48 + t191 * mrSges(6,1) + t79 * mrSges(7,1) + t101 / 0.2e1 - t104 / 0.2e1;
t531 = qJD(1) * qJD(2);
t360 = qJDD(1) * t438 - t433 * t531;
t361 = qJDD(1) * t433 + t438 * t531;
t349 = t432 * t433 - t437 * t438;
t457 = t349 * qJD(3);
t228 = -qJD(1) * t457 + t360 * t432 + t361 * t437;
t426 = qJDD(2) + qJDD(3);
t458 = t472 * qJD(4);
t163 = t228 * t436 + t426 * t431 - t458;
t164 = -qJD(4) * t288 - t228 * t431 + t426 * t436;
t57 = -qJD(5) * t203 + t435 * t163 + t430 * t164;
t658 = t57 / 0.2e1;
t58 = qJD(5) * t450 + t430 * t163 - t435 * t164;
t656 = t58 / 0.2e1;
t229 = -qJD(3) * t333 + t360 * t437 - t432 * t361;
t227 = qJDD(4) - t229;
t222 = qJDD(5) + t227;
t641 = t222 / 0.2e1;
t754 = -mrSges(6,3) - mrSges(7,2);
t705 = Ifges(7,2) + Ifges(6,3);
t704 = -Ifges(6,6) + Ifges(7,6);
t264 = pkin(3) * t333 - pkin(9) * t332;
t183 = t436 * t264 - t269 * t431;
t577 = t332 * t436;
t487 = t333 * pkin(4) - pkin(10) * t577;
t440 = -pkin(10) - pkin(9);
t509 = qJD(4) * t440;
t761 = t436 * t509 - t183 - t487;
t184 = t431 * t264 + t436 * t269;
t528 = pkin(10) * t578;
t760 = -t431 * t509 + t184 - t528;
t581 = t128 * t430;
t47 = t111 * t435 - t581;
t735 = qJD(6) - t47;
t41 = -pkin(5) * t315 + t735;
t758 = t191 * mrSges(6,2) + mrSges(7,2) * t41 - mrSges(6,3) * t47 - t79 * mrSges(7,3);
t633 = -t315 / 0.2e1;
t643 = -t450 / 0.2e1;
t645 = t203 / 0.2e1;
t646 = -t203 / 0.2e1;
t757 = -Ifges(6,2) * t645 + Ifges(7,3) * t646 + t633 * t704 + t643 * t707 - t762;
t756 = mrSges(6,1) + mrSges(7,1);
t755 = mrSges(6,2) - mrSges(7,3);
t200 = Ifges(6,4) * t203;
t588 = Ifges(7,5) * t203;
t694 = t706 * t315 + t450 * t708 - t200 + t588;
t731 = t763 * pkin(4);
t429 = qJ(2) + qJ(3);
t421 = sin(t429);
t428 = qJ(4) + qJ(5);
t420 = sin(t428);
t570 = t420 * t421;
t601 = mrSges(5,2) * t431;
t753 = -mrSges(6,2) * t570 - t421 * t601;
t751 = -Ifges(6,2) * t646 + Ifges(7,3) * t645 + t762;
t750 = -mrSges(5,3) + t754;
t345 = t361 * pkin(7);
t281 = qJDD(2) * pkin(2) - pkin(8) * t361 - t345;
t344 = t360 * pkin(7);
t287 = pkin(8) * t360 + t344;
t536 = qJD(3) * t437;
t537 = qJD(3) * t432;
t143 = t281 * t437 - t432 * t287 - t340 * t537 + t355 * t536;
t137 = -pkin(3) * t426 - t143;
t78 = -pkin(4) * t164 + t137;
t12 = pkin(5) * t58 - qJ(6) * t57 - qJD(6) * t450 + t78;
t583 = qJDD(1) * pkin(1);
t324 = -pkin(2) * t360 - t583;
t129 = -pkin(3) * t229 - pkin(9) * t228 + t324;
t142 = t432 * t281 + t437 * t287 + t340 * t536 + t355 * t537;
t136 = pkin(9) * t426 + t142;
t32 = -qJD(4) * t160 + t436 * t129 - t136 * t431;
t26 = pkin(4) * t227 - pkin(10) * t163 + t32;
t534 = qJD(4) * t436;
t31 = t431 * t129 + t436 * t136 + t238 * t534 - t248 * t535;
t29 = pkin(10) * t164 + t31;
t532 = qJD(5) * t435;
t533 = qJD(5) * t430;
t7 = t111 * t532 - t128 * t533 + t430 * t26 + t435 * t29;
t2 = qJ(6) * t222 + qJD(6) * t315 + t7;
t657 = -t58 / 0.2e1;
t748 = mrSges(6,1) * t78 + mrSges(7,1) * t12 - Ifges(6,4) * t57 / 0.2e1 - Ifges(6,6) * t222 / 0.2e1 + 0.2e1 * Ifges(7,3) * t656 - mrSges(6,3) * t7 - mrSges(7,2) * t2 + (-t657 + t656) * Ifges(6,2) + t764 * t658 + (t704 + Ifges(7,6)) * t641;
t8 = -qJD(5) * t48 + t26 * t435 - t29 * t430;
t4 = -pkin(5) * t222 + qJDD(6) - t8;
t746 = mrSges(6,2) * t78 + mrSges(7,2) * t4 - mrSges(6,3) * t8 - mrSges(7,3) * t12 + Ifges(6,4) * t657 + 0.2e1 * t641 * t706 + t656 * t764 + 0.2e1 * t658 * t708;
t745 = -Ifges(6,4) * t645 - Ifges(7,5) * t646 - t706 * t633 - t708 * t643 + t758;
t744 = t694 / 0.2e1;
t423 = cos(t429);
t683 = t423 * pkin(3) + t421 * pkin(9);
t743 = m(5) * t683;
t741 = Ifges(4,5) * t427;
t740 = Ifges(4,6) * t427;
t738 = t159 * mrSges(5,1);
t737 = t160 * mrSges(5,2);
t736 = -m(7) * qJ(6) - mrSges(7,3);
t378 = t440 * t431;
t424 = t436 * pkin(10);
t608 = pkin(9) * t436;
t380 = t424 + t608;
t470 = t435 * t378 - t380 * t430;
t699 = qJD(5) * t470 + t430 * t761 - t435 * t760;
t291 = t378 * t430 + t380 * t435;
t698 = -qJD(5) * t291 + t430 * t760 + t435 * t761;
t561 = t430 * t436;
t350 = t431 * t435 + t561;
t236 = t350 * t332;
t348 = t430 * t431 - t435 * t436;
t237 = t348 * t332;
t675 = qJD(4) + qJD(5);
t275 = t675 * t348;
t276 = t675 * t350;
t734 = -qJD(6) * t350 + t731 + (t275 - t237) * qJ(6) + (-t236 + t276) * pkin(5);
t596 = mrSges(6,3) * t450;
t176 = mrSges(6,1) * t315 - t596;
t177 = -mrSges(7,1) * t315 + mrSges(7,2) * t450;
t689 = t177 - t176;
t585 = t333 * mrSges(4,3);
t688 = mrSges(4,1) * t427 - mrSges(5,1) * t472 - t288 * mrSges(5,2) - t585;
t733 = t754 * t421;
t732 = -t423 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t421;
t272 = t354 * t432 - t552;
t730 = -pkin(2) * t537 + t272;
t277 = -qJD(2) * t349 - t457;
t580 = t277 * t431;
t462 = t351 * t534 + t580;
t434 = sin(qJ(1));
t555 = t434 * t436;
t439 = cos(qJ(1));
t559 = t431 * t439;
t322 = -t423 * t559 + t555;
t729 = t31 * t436 - t32 * t431;
t602 = mrSges(5,1) * t436;
t728 = t601 - t602;
t725 = Ifges(5,5) * t288 - t472 * Ifges(5,6) + Ifges(5,3) * t327 + t203 * t704 + t315 * t705 + t450 * t706;
t412 = pkin(4) * t436 + pkin(3);
t422 = cos(t428);
t568 = t421 * t422;
t724 = mrSges(6,1) * t568 + (t602 - m(7) * (-pkin(5) * t422 - qJ(6) * t420 - t412) + t422 * mrSges(7,1) + t420 * mrSges(7,3)) * t421;
t119 = pkin(5) * t450 + qJ(6) * t203;
t717 = -m(6) - m(7);
t649 = t163 / 0.2e1;
t648 = t164 / 0.2e1;
t640 = t227 / 0.2e1;
t712 = t360 / 0.2e1;
t624 = t438 / 0.2e1;
t711 = t472 / 0.2e1;
t35 = mrSges(6,1) * t222 - mrSges(6,3) * t57;
t36 = -t222 * mrSges(7,1) + t57 * mrSges(7,2);
t703 = t36 - t35;
t34 = -mrSges(7,2) * t58 + mrSges(7,3) * t222;
t37 = -mrSges(6,2) * t222 - mrSges(6,3) * t58;
t702 = t37 + t34;
t701 = t438 * Ifges(3,2);
t319 = t333 * qJ(6);
t700 = -t319 + t699;
t611 = pkin(5) * t333;
t697 = t611 - t698;
t481 = mrSges(5,1) * t431 + t436 * mrSges(5,2);
t695 = t247 * t481;
t253 = t350 * t351;
t693 = -t730 + t734;
t268 = pkin(3) * t349 - pkin(9) * t351 - t414;
t292 = t379 * t432 - t381 * t437;
t192 = t436 * t268 - t292 * t431;
t573 = t351 * t436;
t148 = pkin(4) * t349 - pkin(10) * t573 + t192;
t284 = t436 * t292;
t193 = t431 * t268 + t284;
t574 = t351 * t431;
t171 = -pkin(10) * t574 + t193;
t692 = t430 * t148 + t435 * t171;
t691 = -t270 + t734;
t174 = -mrSges(7,2) * t203 + mrSges(7,3) * t315;
t597 = mrSges(6,3) * t203;
t175 = -mrSges(6,2) * t315 - t597;
t690 = t174 + t175;
t687 = -t270 + t731;
t686 = -t730 + t731;
t685 = t437 * t379 + t381 * t432;
t681 = t222 * t705 + t57 * t706 + t58 * t704;
t609 = pkin(7) * t438;
t610 = pkin(7) * t433;
t680 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t540) * t609 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t539) * t610;
t551 = t439 * t420;
t556 = t434 * t422;
t308 = t423 * t551 - t556;
t563 = t423 * t439;
t309 = t420 * t434 + t422 * t563;
t679 = t308 * t756 + t755 * t309;
t564 = t423 * t434;
t306 = t420 * t564 + t422 * t439;
t307 = t423 * t556 - t551;
t678 = t306 * t756 + t755 * t307;
t677 = t344 * t438 + t345 * t433;
t676 = g(1) * t439 + g(2) * t434;
t674 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t673 = t439 * t753 + t563 * t750;
t562 = t423 * t440;
t469 = -t412 * t421 - t562;
t616 = pkin(3) * t421;
t618 = pkin(2) * t433;
t672 = -m(7) * (-t562 - t618) - m(6) * (t469 - t618) - m(5) * (-t616 - t618) + t724;
t671 = -m(6) * t469 + m(7) * t562 + t724;
t278 = t427 * t351;
t538 = qJD(2) * t433;
t520 = pkin(2) * t538;
t190 = pkin(3) * t278 - pkin(9) * t277 + t520;
t510 = qJD(2) * t441;
t358 = t433 * t510;
t359 = t438 * t510;
t208 = qJD(3) * t685 + t358 * t437 + t359 * t432;
t494 = t436 * t190 - t208 * t431;
t579 = t277 * t436;
t40 = -pkin(10) * t579 + pkin(4) * t278 + (-t284 + (pkin(10) * t351 - t268) * t431) * qJD(4) + t494;
t74 = t431 * t190 + t436 * t208 + t268 * t534 - t292 * t535;
t50 = -pkin(10) * t462 + t74;
t14 = -qJD(5) * t692 + t40 * t435 - t430 * t50;
t670 = m(7) * pkin(5) + t756;
t376 = -mrSges(3,1) * t438 + mrSges(3,2) * t433;
t669 = m(3) * pkin(1) + mrSges(2,1) - t376 - t732;
t565 = t422 * t423;
t569 = t420 * t423;
t668 = t423 * t728 - t565 * t756 + t755 * t569 + t732 + t733;
t667 = -mrSges(6,2) - t736;
t666 = t32 * mrSges(5,1) - t31 * mrSges(5,2);
t664 = t753 * t434 + (-m(5) * pkin(9) + t750) * t564;
t107 = mrSges(5,1) * t227 - mrSges(5,3) * t163;
t459 = t472 * mrSges(5,3);
t231 = -t327 * mrSges(5,2) - t459;
t598 = mrSges(5,3) * t288;
t232 = mrSges(5,1) * t327 - t598;
t663 = m(5) * ((-t159 * t436 - t160 * t431) * qJD(4) + t729) - t232 * t534 - t231 * t535 - t107 * t431;
t662 = t8 * mrSges(6,1) - t4 * mrSges(7,1) - t7 * mrSges(6,2) + t2 * mrSges(7,3);
t655 = Ifges(5,1) * t649 + Ifges(5,4) * t648 + Ifges(5,5) * t640;
t642 = t450 / 0.2e1;
t635 = -t288 / 0.2e1;
t634 = t288 / 0.2e1;
t632 = t315 / 0.2e1;
t631 = -t327 / 0.2e1;
t630 = -t332 / 0.2e1;
t628 = t333 / 0.2e1;
t625 = t436 / 0.2e1;
t619 = pkin(2) * t432;
t617 = pkin(2) * t437;
t615 = pkin(4) * t288;
t614 = pkin(4) * t430;
t613 = pkin(4) * t431;
t612 = pkin(4) * t435;
t605 = g(3) * t421;
t410 = pkin(9) + t619;
t603 = -pkin(10) - t410;
t600 = mrSges(6,2) * t422;
t599 = mrSges(4,3) * t332;
t595 = Ifges(3,4) * t433;
t594 = Ifges(3,4) * t438;
t593 = Ifges(4,4) * t333;
t592 = Ifges(5,4) * t288;
t591 = Ifges(5,4) * t431;
t590 = Ifges(5,4) * t436;
t571 = t410 * t436;
t567 = t421 * t439;
t566 = t421 * t440;
t366 = t423 * t412;
t560 = t431 * t434;
t553 = t436 * t439;
t550 = t439 * t441;
t521 = pkin(2) * t540;
t243 = t264 + t521;
t273 = t354 * t437 + t334;
t180 = t436 * t243 - t273 * t431;
t132 = t180 + t487;
t181 = t431 * t243 + t436 * t273;
t147 = -t528 + t181;
t71 = t430 * t132 + t435 * t147;
t519 = pkin(2) * t536;
t518 = pkin(4) * t533;
t517 = pkin(4) * t532;
t512 = Ifges(5,5) * t163 + Ifges(5,6) * t164 + Ifges(5,3) * t227;
t506 = t351 * t535;
t286 = Ifges(5,4) * t472;
t187 = Ifges(5,1) * t288 + Ifges(5,5) * t327 - t286;
t504 = t187 * t625;
t501 = -t535 / 0.2e1;
t500 = t736 * t568;
t499 = qJD(4) * t603;
t498 = t531 / 0.2e1;
t496 = -t306 * pkin(5) + qJ(6) * t307;
t495 = -t308 * pkin(5) + qJ(6) * t309;
t493 = t366 - t566;
t492 = t439 * t414 - t434 * t441;
t490 = t431 * t519;
t489 = t436 * t519;
t239 = pkin(4) * t574 - t685;
t484 = mrSges(3,1) * t433 + mrSges(3,2) * t438;
t482 = mrSges(4,1) * t421 + mrSges(4,2) * t423;
t480 = Ifges(5,1) * t436 - t591;
t479 = t595 + t701;
t478 = -Ifges(5,2) * t431 + t590;
t477 = Ifges(3,5) * t438 - Ifges(3,6) * t433;
t476 = Ifges(5,5) * t436 - Ifges(5,6) * t431;
t70 = t132 * t435 - t147 * t430;
t84 = t148 * t435 - t171 * t430;
t341 = t603 * t431;
t342 = t424 + t571;
t471 = t435 * t341 - t342 * t430;
t261 = t341 * t430 + t342 * t435;
t468 = t322 * pkin(4);
t465 = pkin(5) * t565 + qJ(6) * t569 + t493;
t463 = pkin(1) * t484;
t320 = t423 * t560 + t553;
t461 = t506 - t579;
t460 = t433 * (Ifges(3,1) * t438 - t595);
t13 = t148 * t532 - t171 * t533 + t430 * t40 + t435 * t50;
t455 = t320 * pkin(4);
t265 = pkin(5) * t348 - qJ(6) * t350 - t412;
t452 = t436 * t499 - t490;
t209 = qJD(3) * t292 + t358 * t432 - t437 * t359;
t448 = t662 + t681;
t134 = pkin(4) * t462 + t209;
t186 = -t472 * Ifges(5,2) + Ifges(5,6) * t327 + t592;
t244 = Ifges(4,2) * t332 + t593 + t740;
t325 = Ifges(4,4) * t332;
t245 = Ifges(4,1) * t333 + t325 + t741;
t67 = Ifges(5,4) * t163 + Ifges(5,2) * t164 + Ifges(5,6) * t227;
t442 = -(Ifges(4,1) * t332 - t593 + t725) * t333 / 0.2e1 + t757 * t236 - (Ifges(6,4) * t646 + Ifges(7,5) * t645 + t632 * t706 + t642 * t708 + t758) * t275 + (-t763 * t160 + (-t534 + t577) * t159 + t729) * mrSges(5,3) + (-t47 * mrSges(6,1) + t41 * mrSges(7,1) + Ifges(6,6) * t645 + Ifges(7,6) * t646 + t48 * mrSges(6,2) - t42 * mrSges(7,3) + t737 - t738 + t377 * mrSges(4,1) - Ifges(4,2) * t630 + Ifges(5,3) * t631 + Ifges(5,5) * t635 + t740 / 0.2e1 + Ifges(5,6) * t711 + t706 * t643 + t705 * t633) * t333 + t694 * (t237 / 0.2e1 - t275 / 0.2e1) + (t325 + t245) * t630 + (t377 * mrSges(4,2) + t476 * t631 + t480 * t635 - t695 - t741 / 0.2e1 + t478 * t711) * t332 + t745 * t237 + (Ifges(5,5) * t431 + Ifges(5,6) * t436) * t640 + (Ifges(5,2) * t436 + t591) * t648 + t67 * t625 + t244 * t628 + t270 * t585 + t269 * t599 + (Ifges(5,1) * t431 + t590) * t649 + t431 * t655 - t187 * t577 / 0.2e1 + t137 * t728 + t746 * t350 + (t501 + t578 / 0.2e1) * t186 + (t504 + t695) * qJD(4) + (t632 * t704 + t642 * t707 + t751) * t276 + t748 * t348 + Ifges(4,6) * t229 + Ifges(4,5) * t228 + (t288 * t480 + t327 * t476) * qJD(4) / 0.2e1 - t142 * mrSges(4,2) + t143 * mrSges(4,1) - t478 * t458 / 0.2e1 + Ifges(4,3) * t426;
t416 = Ifges(3,4) * t539;
t413 = -pkin(3) - t617;
t411 = -pkin(5) - t612;
t405 = qJ(6) + t614;
t396 = qJD(6) + t517;
t393 = pkin(9) * t563;
t373 = -t412 - t617;
t331 = Ifges(3,1) * t540 + Ifges(3,5) * qJD(2) + t416;
t330 = Ifges(3,6) * qJD(2) + qJD(1) * t479;
t323 = t423 * t553 + t560;
t321 = -t423 * t555 + t559;
t294 = -mrSges(4,2) * t427 + t599;
t293 = t431 * t499 + t489;
t263 = -mrSges(4,1) * t332 + mrSges(4,2) * t333;
t254 = t348 * t351;
t250 = t265 - t617;
t214 = -mrSges(4,2) * t426 + mrSges(4,3) * t229;
t213 = mrSges(4,1) * t426 - mrSges(4,3) * t228;
t153 = qJD(5) * t261 + t293 * t430 - t435 * t452;
t152 = qJD(5) * t471 + t435 * t293 + t430 * t452;
t122 = pkin(5) * t253 + qJ(6) * t254 + t239;
t121 = mrSges(6,1) * t203 + mrSges(6,2) * t450;
t120 = mrSges(7,1) * t203 - mrSges(7,3) * t450;
t108 = -mrSges(5,2) * t227 + mrSges(5,3) * t164;
t99 = t119 + t615;
t98 = t277 * t561 - t430 * t506 - t533 * t574 + (t573 * t675 + t580) * t435;
t97 = -t253 * t675 - t348 * t277;
t86 = -mrSges(5,1) * t164 + mrSges(5,2) * t163;
t77 = -pkin(5) * t349 - t84;
t76 = qJ(6) * t349 + t692;
t75 = -qJD(4) * t193 + t494;
t62 = -t70 - t611;
t61 = t319 + t71;
t60 = t127 * t435 - t581;
t59 = t127 * t430 + t554;
t27 = pkin(5) * t98 - qJ(6) * t97 + qJD(6) * t254 + t134;
t24 = mrSges(6,1) * t58 + mrSges(6,2) * t57;
t23 = mrSges(7,1) * t58 - mrSges(7,3) * t57;
t11 = -pkin(5) * t278 - t14;
t9 = qJ(6) * t278 + qJD(6) * t349 + t13;
t1 = [t725 * t278 / 0.2e1 + (t191 * t97 - t278 * t48) * mrSges(6,2) + (-t269 * t277 - t270 * t278) * mrSges(4,3) + (Ifges(7,5) * t97 + Ifges(7,6) * t278) * t645 + (t438 * t594 + t460) * t498 - (-m(4) * t143 + m(5) * t137 - t213 + t86) * t685 + (Ifges(6,4) * t97 + Ifges(6,6) * t278) * t646 - t278 * t737 + t361 * t594 / 0.2e1 + t278 * t738 + t97 * t744 + (t159 * t461 - t160 * t462 - t31 * t574 - t32 * t573) * mrSges(5,3) + (Ifges(4,1) * t277 - Ifges(4,4) * t278) * t628 + (-Ifges(5,1) * t461 - Ifges(5,4) * t462 + Ifges(5,5) * t278) * t634 + t239 * t24 + t277 * t504 + (t324 * mrSges(4,2) - t143 * mrSges(4,3) + Ifges(4,1) * t228 + Ifges(4,4) * t229 + Ifges(4,5) * t426 + t137 * t481 + t187 * t501 + t476 * t640 + t478 * t648 + t480 * t649) * t351 + t573 * t655 - t67 * t574 / 0.2e1 - t376 * t583 - t330 * t538 / 0.2e1 + (Ifges(3,4) * t361 + Ifges(3,2) * t360) * t624 + m(7) * (t11 * t41 + t12 * t122 + t2 * t76 + t27 * t79 + t4 * t77 + t42 * t9) + (t278 * t42 - t79 * t97) * mrSges(7,3) - t746 * t254 - t463 * t531 + t751 * t98 + t748 * t253 + Ifges(2,3) * qJDD(1) + t74 * t231 + t75 * t232 + m(6) * (t13 * t48 + t134 * t191 + t14 * t47 + t239 * t78 + t692 * t7 + t8 * t84) + t692 * t37 + m(4) * (t142 * t292 + t208 * t270 - t324 * t414 - t377 * t520) + t192 * t107 + t193 * t108 + t9 * t174 + t13 * t175 + t14 * t176 + t11 * t177 + t134 * t121 + t122 * t23 + t27 * t120 + m(5) * (t159 * t75 + t160 * t74 + t192 * t32 + t193 * t31) + t84 * t35 + t76 * t34 + t77 * t36 - t472 * (-Ifges(5,4) * t461 - Ifges(5,2) * t462 + Ifges(5,6) * t278) / 0.2e1 + (-mrSges(3,1) * t610 - mrSges(3,2) * t609 + 0.2e1 * Ifges(3,6) * t624) * qJDD(2) + (Ifges(3,1) * t361 + Ifges(3,4) * t712 + Ifges(3,5) * qJDD(2) - t498 * t701) * t433 + t247 * (mrSges(5,1) * t462 - mrSges(5,2) * t461) + t327 * (-Ifges(5,5) * t461 - Ifges(5,6) * t462 + Ifges(5,3) * t278) / 0.2e1 + (t278 * t705 + t704 * t98 + t706 * t97) * t632 + (t324 * mrSges(4,1) - t142 * mrSges(4,3) - Ifges(4,4) * t228 + Ifges(5,5) * t649 - Ifges(4,2) * t229 - Ifges(4,6) * t426 + Ifges(5,6) * t648 + Ifges(6,6) * t657 + Ifges(7,6) * t656 + Ifges(5,3) * t640 + t641 * t705 + t658 * t706 + t662 + t666) * t349 + (t278 * t706 + t707 * t98 + t708 * t97) * t642 + (-m(4) * t269 + m(5) * t247 - t688) * t209 - t462 * t186 / 0.2e1 + t427 * (Ifges(4,5) * t277 - Ifges(4,6) * t278) / 0.2e1 + (t360 * t609 + t361 * t610 + t677) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t677) + (t331 * t624 + t477 * qJD(2) / 0.2e1 - t680) * qJD(2) + (t512 + t681) * t349 / 0.2e1 + (-t323 * mrSges(5,1) - t322 * mrSges(5,2) + t754 * t567 + (-m(4) - m(5)) * t492 + t717 * (pkin(4) * t560 + t412 * t563 - t439 * t566 + t492) - t670 * t309 - t667 * t308 + t674 * t434 + (-t669 - t743) * t439) * g(2) + t277 * t245 / 0.2e1 + t47 * (mrSges(6,1) * t278 - mrSges(6,3) * t97) + t41 * (-mrSges(7,1) * t278 + mrSges(7,2) * t97) - t278 * t244 / 0.2e1 + t292 * t214 + t208 * t294 + t263 * t520 + t332 * (Ifges(4,4) * t277 - Ifges(4,2) * t278) / 0.2e1 + t479 * t712 + (m(5) * t550 - t321 * mrSges(5,1) - t320 * mrSges(5,2) + t717 * (pkin(4) * t559 + t434 * t566 - t550) + t670 * t307 + t667 * t306 + (m(4) * t441 + t674) * t439 + (-m(5) * (-t414 - t683) + m(4) * t414 + t717 * (-t414 - t366) + t669 - t733) * t434) * g(1) - pkin(1) * (-mrSges(3,1) * t360 + mrSges(3,2) * t361) - t377 * (mrSges(4,1) * t278 + mrSges(4,2) * t277) - t414 * (-mrSges(4,1) * t229 + mrSges(4,2) * t228); (t680 + (-t460 / 0.2e1 + t463) * qJD(1)) * qJD(1) + (-t273 + t519) * t294 - (-Ifges(3,2) * t540 + t331 + t416) * t539 / 0.2e1 + (-t159 * t180 - t160 * t181 - t247 * t272 + t137 * t413 + (t247 * t432 + (-t159 * t431 + t160 * t436) * t437) * qJD(3) * pkin(2)) * m(5) + t688 * t730 + (m(4) * t618 + t482 + t484) * t676 + (-t490 - t180) * t232 + t213 * t617 + t214 * t619 + t108 * t571 + t442 + t330 * t540 / 0.2e1 - t263 * t521 + (t489 - t181) * t231 + t250 * t23 - t477 * t531 / 0.2e1 - t703 * t471 + (t12 * t250 + t2 * t261 - t471 * t4 + t693 * t79 + (t152 - t61) * t42 + (t153 - t62) * t41) * m(7) + (t471 * t8 + t261 * t7 + t373 * t78 + (t152 - t71) * t48 + (-t153 - t70) * t47 + t686 * t191) * m(6) + (-m(4) * t425 - m(7) * (t425 + t465) - m(6) * (t425 + t493) - m(5) * (t425 + t683) + t376 + t668) * g(3) + (t269 * t272 - t270 * t273 + t377 * t521 + (t142 * t432 + t143 * t437 + (-t269 * t432 + t270 * t437) * qJD(3)) * pkin(2)) * m(4) - t61 * t174 - t71 * t175 - t70 * t176 - t62 * t177 + Ifges(3,3) * qJDD(2) + t702 * t261 + t693 * t120 + t686 * t121 + t689 * t153 + t690 * t152 + t663 * t410 - t344 * mrSges(3,2) - t345 * mrSges(3,1) + Ifges(3,6) * t360 + Ifges(3,5) * t361 + t373 * t24 + t413 * t86 + (t434 * t672 + t664) * g(2) + (-m(5) * t393 + t439 * t672 + t673) * g(1); t108 * t608 + t442 + t265 * t23 + (-m(6) * t493 - m(7) * t465 + t668 - t743) * g(3) + (-pkin(3) * t137 - t159 * t183 - t160 * t184 - t247 * t270) * m(5) - t184 * t231 - t183 * t232 - t703 * t470 + (t191 * t687 + t291 * t7 - t412 * t78 + t47 * t698 + t470 * t8 + t48 * t699) * m(6) + (t12 * t265 + t2 * t291 - t4 * t470 + t41 * t697 + t42 * t700 + t691 * t79) * m(7) - pkin(3) * t86 + t702 * t291 + t697 * t177 + t698 * t176 + t699 * t175 + t700 * t174 + t691 * t120 + t687 * t121 + t688 * t270 + t676 * t482 - t269 * t294 + t663 * pkin(9) - t412 * t24 + ((m(5) * t616 + t671) * t434 + t664) * g(2) + (-m(5) * (-pkin(3) * t567 + t393) + t671 * t439 + t673) * g(1); (-t191 * t615 + t47 * t59 - t48 * t60 + t613 * t605 + (t430 * t7 + t435 * t8 + (-t430 * t47 + t435 * t48) * qJD(5)) * pkin(4)) * m(6) + (-Ifges(5,2) * t288 + t187 - t286) * t711 + t757 * t450 + (t744 + t745) * t203 + t689 * (-t59 + t518) + m(7) * (t2 * t405 + t396 * t42 + t4 * t411 + t41 * t518) + t666 + (-Ifges(5,5) * t472 - Ifges(5,6) * t288) * t631 + t186 * t634 + (-Ifges(5,1) * t472 - t592) * t635 + t35 * t612 + t37 * t614 + (mrSges(6,1) * t420 + t481 + t600) * t605 - t121 * t615 - g(3) * ((m(7) * (-pkin(5) * t420 - t613) - t420 * mrSges(7,1)) * t421 - t500) + (-t60 + t517) * t175 + t448 - m(7) * (t41 * t59 + t42 * t60 + t79 * t99) + (t598 + t232) * t160 + t512 - t99 * t120 + (-t60 + t396) * t174 - t247 * (t288 * mrSges(5,1) - mrSges(5,2) * t472) + (-m(7) * (-t455 + t496) + m(6) * t455 + mrSges(5,1) * t320 - mrSges(5,2) * t321 + t678) * g(2) + (-m(7) * (t468 + t495) - m(6) * t468 - mrSges(5,1) * t322 + mrSges(5,2) * t323 + t679) * g(1) + (-t231 - t459) * t159 + t405 * t34 + t411 * t36; t104 * t642 + (Ifges(7,3) * t450 - t588) * t646 + (t596 - t689) * t48 + (-t597 - t690) * t47 + (t203 * t41 + t42 * t450) * mrSges(7,2) + t448 + t678 * g(2) + ((t420 * t670 + t600) * t421 + t500) * g(3) + t679 * g(1) - t191 * (mrSges(6,1) * t450 - mrSges(6,2) * t203) - t79 * (mrSges(7,1) * t450 + mrSges(7,3) * t203) + qJD(6) * t174 - t119 * t120 - pkin(5) * t36 + qJ(6) * t34 + (-t203 * t706 + t450 * t704) * t633 + (-Ifges(6,2) * t450 - t200 + t694) * t645 + (-t203 * t708 + t101 + t199 - t589) * t643 + (-pkin(5) * t4 - t495 * g(1) - t496 * g(2) + qJ(6) * t2 - t119 * t79 - t41 * t48 + t42 * t735) * m(7); t450 * t120 - t315 * t174 + (-g(1) * t308 - g(2) * t306 - g(3) * t570 - t42 * t315 + t450 * t79 + t4) * m(7) + t36;];
tau  = t1;
