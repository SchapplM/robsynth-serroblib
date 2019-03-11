% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:58
% EndTime: 2019-03-09 12:09:09
% DurationCPUTime: 46.43s
% Computational Cost: add. (19884->990), mult. (54965->1320), div. (0->0), fcn. (44970->12), ass. (0->440)
t346 = sin(pkin(6));
t345 = sin(pkin(11));
t350 = sin(qJ(2));
t354 = cos(qJ(2));
t525 = cos(pkin(11));
t376 = -t350 * t345 + t354 * t525;
t365 = t346 * t376;
t278 = qJD(1) * t365;
t622 = t278 - qJD(4);
t663 = Ifges(6,1) + Ifges(7,1);
t661 = Ifges(7,4) + Ifges(6,5);
t483 = qJD(1) * qJD(2);
t297 = (qJDD(1) * t354 - t350 * t483) * t346;
t298 = (qJDD(1) * t350 + t354 * t483) * t346;
t219 = t345 * t297 + t298 * t525;
t437 = t350 * t525;
t494 = qJD(1) * t346;
t460 = t354 * t494;
t279 = -t345 * t460 - t437 * t494;
t347 = cos(pkin(6));
t335 = qJD(1) * t347 + qJD(2);
t349 = sin(qJ(4));
t353 = cos(qJ(4));
t241 = t279 * t349 + t335 * t353;
t481 = qJDD(1) * t347;
t334 = qJDD(2) + t481;
t128 = qJD(4) * t241 + t219 * t353 + t334 * t349;
t583 = t128 / 0.2e1;
t687 = Ifges(5,4) * t583;
t559 = pkin(1) * t347;
t339 = t354 * t559;
t332 = qJD(1) * t339;
t547 = pkin(8) + qJ(3);
t449 = t547 * t350;
t423 = t346 * t449;
t263 = -qJD(1) * t423 + t332;
t508 = t347 * t350;
t338 = pkin(1) * t508;
t510 = t346 * t354;
t264 = (t510 * t547 + t338) * qJD(1);
t435 = t525 * t264;
t172 = t263 * t345 + t435;
t686 = -t172 - t622 * (pkin(4) * t349 - pkin(10) * t353);
t386 = t279 * t353 - t335 * t349;
t129 = qJD(4) * t386 - t219 * t349 + t334 * t353;
t127 = qJDD(5) - t129;
t584 = t127 / 0.2e1;
t218 = t297 * t525 - t345 * t298;
t213 = qJDD(4) - t218;
t348 = sin(qJ(5));
t352 = cos(qJ(5));
t488 = qJD(5) * t348;
t361 = -t376 * t494 + qJD(4);
t637 = qJD(5) * t361 + t128;
t62 = t348 * t213 + t352 * t637 + t386 * t488;
t593 = t62 / 0.2e1;
t685 = t584 * t661 + t593 * t663;
t486 = qJD(5) * t352;
t63 = -t352 * t213 + t348 * t637 - t386 * t486;
t592 = -t63 / 0.2e1;
t591 = t63 / 0.2e1;
t662 = -Ifges(6,4) + Ifges(7,5);
t684 = Ifges(6,6) - Ifges(7,6);
t683 = Ifges(3,3) + Ifges(4,3);
t659 = Ifges(6,3) + Ifges(7,2);
t240 = qJD(5) - t241;
t682 = t662 * t591 + t685;
t157 = -t348 * t386 - t352 * t361;
t153 = Ifges(6,4) * t157;
t158 = t348 * t361 - t352 * t386;
t535 = Ifges(7,5) * t157;
t654 = t158 * t663 + t661 * t240 - t153 + t535;
t582 = t129 / 0.2e1;
t681 = Ifges(5,4) * t582;
t462 = t525 * pkin(2);
t343 = -t462 - pkin(3);
t551 = t353 * pkin(4);
t311 = -t349 * pkin(10) + t343 - t551;
t254 = t345 * t264;
t173 = t263 * t525 - t254;
t461 = t350 * t494;
t430 = pkin(2) * t461;
t185 = -pkin(3) * t279 - pkin(9) * t278 + t430;
t102 = t353 * t173 + t349 * t185;
t88 = -pkin(10) * t279 + t102;
t680 = t311 * t486 + t348 * t686 - t352 * t88;
t558 = pkin(2) * t345;
t342 = pkin(9) + t558;
t500 = t352 * t353;
t635 = t348 * t311 + t342 * t500;
t679 = -qJD(5) * t635 + t348 * t88 + t352 * t686;
t243 = pkin(2) * t335 + t263;
t151 = t345 * t243 + t435;
t141 = pkin(9) * t335 + t151;
t344 = pkin(2) * t354 + pkin(1);
t301 = -t344 * t494 + qJD(3);
t167 = -pkin(3) * t278 + pkin(9) * t279 + t301;
t86 = t353 * t141 + t349 * t167;
t74 = pkin(10) * t361 + t86;
t150 = t243 * t525 - t254;
t140 = -t335 * pkin(3) - t150;
t80 = -t241 * pkin(4) + pkin(10) * t386 + t140;
t28 = -t348 * t74 + t352 * t80;
t29 = t348 * t80 + t352 * t74;
t482 = qJDD(1) * t346;
t474 = pkin(1) * t482;
t259 = -pkin(2) * t297 + qJDD(3) - t474;
t107 = -pkin(3) * t218 - pkin(9) * t219 + t259;
t484 = t353 * qJD(4);
t485 = t349 * qJD(4);
t495 = pkin(8) * t510 + t338;
t295 = t495 * qJD(2);
t473 = pkin(1) * t481;
t330 = t354 * t473;
t512 = t346 * t350;
t458 = qJD(3) * t512;
t472 = pkin(8) * t482;
t146 = -t350 * t472 + pkin(2) * t334 - qJ(3) * t298 + t330 + (-t295 - t458) * qJD(1);
t478 = qJD(2) * t559;
t425 = qJD(1) * t478;
t463 = t350 * t473 + (t425 + t472) * t354;
t492 = qJD(3) * t354;
t493 = qJD(2) * t350;
t160 = qJ(3) * t297 + (-pkin(8) * t493 + t492) * t494 + t463;
t84 = t345 * t146 + t525 * t160;
t79 = pkin(9) * t334 + t84;
t22 = t349 * t107 - t141 * t485 + t167 * t484 + t353 * t79;
t18 = pkin(10) * t213 + t22;
t83 = t146 * t525 - t345 * t160;
t78 = -t334 * pkin(3) - t83;
t27 = -t129 * pkin(4) - t128 * pkin(10) + t78;
t3 = t352 * t18 + t348 * t27 + t80 * t486 - t488 * t74;
t4 = -qJD(5) * t29 - t18 * t348 + t27 * t352;
t415 = t3 * t352 - t348 * t4;
t678 = -t28 * t486 - t29 * t488 + t415;
t638 = qJD(6) - t28;
t24 = -pkin(5) * t240 + t638;
t25 = qJ(6) * t240 + t29;
t1 = qJ(6) * t127 + qJD(6) * t240 + t3;
t2 = -pkin(5) * t127 + qJDD(6) - t4;
t416 = t1 * t352 + t2 * t348;
t677 = t24 * t486 - t25 * t488 + t416;
t355 = cos(qJ(1));
t499 = t354 * t355;
t351 = sin(qJ(1));
t502 = t351 * t350;
t676 = t347 * t499 - t502;
t636 = mrSges(4,1) * t335 + mrSges(5,1) * t241 + mrSges(5,2) * t386 + mrSges(4,3) * t279;
t675 = m(5) * t140 - t636;
t674 = mrSges(7,2) * t2 - mrSges(6,3) * t4 + t682;
t673 = Ifges(6,4) * t593 + Ifges(6,6) * t584 - Ifges(7,5) * t62 / 0.2e1 - Ifges(7,6) * t127 / 0.2e1 + (Ifges(6,2) + Ifges(7,3)) * t592;
t377 = -t354 * t345 - t437;
t288 = t377 * t347;
t231 = -t288 * t355 + t351 * t376;
t509 = t346 * t355;
t199 = t231 * t353 - t349 * t509;
t287 = t376 * t347;
t230 = t355 * t287 + t351 * t377;
t119 = t199 * t348 + t230 * t352;
t672 = t199 * t352 - t230 * t348;
t574 = t213 / 0.2e1;
t606 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t671 = -t659 * t584 - t661 * t593 - Ifges(6,6) * t592 - Ifges(7,6) * t591 + t687 - t606 - (-t213 / 0.2e1 - t574) * Ifges(5,6) - (-t129 / 0.2e1 - t582) * Ifges(5,2);
t670 = -m(6) - m(7);
t85 = -t349 * t141 + t353 * t167;
t73 = -pkin(4) * t361 - t85;
t669 = m(6) * t73;
t564 = t334 / 0.2e1;
t667 = t85 * mrSges(5,1);
t666 = t86 * mrSges(5,2);
t665 = mrSges(4,2) - mrSges(5,3);
t664 = mrSges(6,3) + mrSges(7,2);
t658 = t127 * t659 + t62 * t661 - t63 * t684;
t21 = mrSges(6,1) * t63 + mrSges(6,2) * t62;
t90 = mrSges(5,1) * t213 - mrSges(5,3) * t128;
t656 = t21 - t90;
t655 = -t157 * t684 + t158 * t661 + t240 * t659;
t653 = Ifges(3,5) * t298;
t652 = Ifges(4,5) * t219;
t651 = Ifges(5,5) * t361;
t650 = Ifges(3,6) * t297;
t649 = Ifges(4,6) * t218;
t648 = Ifges(5,6) * t361;
t412 = mrSges(5,1) * t353 - mrSges(5,2) * t349;
t647 = -mrSges(4,1) - t412;
t527 = t386 * mrSges(5,3);
t166 = mrSges(5,1) * t361 + t527;
t98 = mrSges(6,1) * t157 + mrSges(6,2) * t158;
t646 = t166 - t98;
t506 = t348 * t353;
t196 = t278 * t506 + t352 * t279;
t197 = t278 * t500 - t279 * t348;
t391 = pkin(5) * t348 - qJ(6) * t352;
t384 = t342 + t391;
t392 = pkin(5) * t352 + qJ(6) * t348;
t101 = -t349 * t173 + t185 * t353;
t87 = pkin(4) * t279 - t101;
t645 = -pkin(5) * t196 + qJ(6) * t197 + (qJD(5) * t392 - qJD(6) * t352) * t349 + t384 * t484 - t87;
t644 = -qJD(6) * t348 + t240 * t391 - t86;
t518 = t278 * t349;
t643 = -qJ(6) * t518 + (-t342 * t488 - qJD(6)) * t353 + (-t342 * t352 + qJ(6)) * t485 + t680;
t441 = t342 * t348 + pkin(5);
t642 = pkin(5) * t518 - t441 * t485 - t679;
t457 = t342 * t485;
t641 = t348 * t457 + t679;
t640 = (-t352 * t485 - t353 * t488) * t342 + t680;
t260 = pkin(2) * t347 + t339 - t423;
t274 = qJ(3) * t510 + t495;
t180 = t260 * t525 - t345 * t274;
t168 = -t347 * pkin(3) - t180;
t286 = t377 * t346;
t257 = -t286 * t349 - t347 * t353;
t258 = -t286 * t353 + t347 * t349;
t438 = -t257 * pkin(4) + t258 * pkin(10);
t103 = t168 - t438;
t181 = t345 * t260 + t525 * t274;
t169 = pkin(9) * t347 + t181;
t336 = pkin(2) * t510;
t418 = pkin(3) * t365 - pkin(9) * t286 + t336;
t560 = pkin(1) * t346;
t195 = -t418 - t560;
t105 = t353 * t169 + t349 * t195;
t95 = -pkin(10) * t365 + t105;
t639 = t348 * t103 + t352 * t95;
t409 = -t352 * mrSges(7,1) - t348 * mrSges(7,3);
t411 = mrSges(6,1) * t352 - mrSges(6,2) * t348;
t634 = m(7) * t392 - t409 + t411;
t34 = t157 * pkin(5) - t158 * qJ(6) + t73;
t408 = mrSges(7,1) * t348 - mrSges(7,3) * t352;
t410 = mrSges(6,1) * t348 + mrSges(6,2) * t352;
t633 = t34 * t408 + t73 * t410;
t632 = t348 * t661 + t352 * t684;
t631 = -t348 * t684 + t352 * t661;
t534 = Ifges(7,5) * t348;
t537 = Ifges(6,4) * t348;
t630 = t352 * t663 + t534 - t537;
t533 = Ifges(7,5) * t352;
t536 = Ifges(6,4) * t352;
t629 = t348 * t663 - t533 + t536;
t628 = t334 * t683 + t649 + t650 + t652 + t653;
t23 = t107 * t353 - t141 * t484 - t167 * t485 - t349 * t79;
t625 = t22 * t353 - t23 * t349;
t152 = Ifges(7,5) * t158;
t67 = Ifges(7,6) * t240 + Ifges(7,3) * t157 + t152;
t538 = Ifges(6,4) * t158;
t70 = -Ifges(6,2) * t157 + Ifges(6,6) * t240 + t538;
t624 = t70 / 0.2e1 - t67 / 0.2e1;
t588 = -t70 / 0.2e1;
t623 = t588 + t67 / 0.2e1;
t543 = Ifges(3,4) * t350;
t601 = t346 ^ 2;
t621 = (t350 * (Ifges(3,1) * t354 - t543) / 0.2e1 - pkin(1) * (mrSges(3,1) * t350 + mrSges(3,2) * t354)) * t601;
t620 = 0.2e1 * t564;
t619 = mrSges(5,2) - t664;
t618 = mrSges(5,1) + t634;
t617 = m(3) * pkin(1) + mrSges(2,1);
t30 = -mrSges(7,2) * t63 + mrSges(7,3) * t127;
t33 = -mrSges(6,2) * t127 - mrSges(6,3) * t63;
t615 = (t30 + t33) * t352;
t31 = mrSges(6,1) * t127 - mrSges(6,3) * t62;
t32 = -t127 * mrSges(7,1) + t62 * mrSges(7,2);
t614 = (-t31 + t32) * t348;
t280 = qJD(2) * t286;
t333 = t354 * t478;
t246 = t333 + (-qJD(2) * t449 + t492) * t346;
t450 = t547 * t346;
t247 = -t458 + (-t354 * t450 - t338) * qJD(2);
t144 = t246 * t525 + t345 * t247;
t281 = qJD(2) * t365;
t459 = t346 * t493;
t429 = pkin(2) * t459;
t186 = -pkin(3) * t280 - pkin(9) * t281 + t429;
t48 = t353 * t144 - t169 * t485 + t349 * t186 + t195 * t484;
t40 = -pkin(10) * t280 + t48;
t143 = t246 * t345 - t525 * t247;
t178 = qJD(4) * t258 + t281 * t349;
t179 = -qJD(4) * t257 + t281 * t353;
t66 = pkin(4) * t178 - pkin(10) * t179 + t143;
t9 = -qJD(5) * t639 - t348 * t40 + t352 * t66;
t612 = -m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3);
t431 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t424 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t610 = t23 * mrSges(5,1) - t22 * mrSges(5,2);
t609 = 0.2e1 * Ifges(5,1) * t583 + 0.2e1 * Ifges(5,5) * t574 + t681;
t608 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 - t673;
t232 = -t351 * t288 - t355 * t376;
t511 = t346 * t351;
t203 = -t232 * t353 + t349 * t511;
t607 = -g(1) * t203 - g(2) * t199 - g(3) * t258;
t605 = -t73 * mrSges(6,1) - t34 * mrSges(7,1) + t25 * mrSges(7,2) + t29 * mrSges(6,3);
t428 = pkin(8) * t459;
t237 = -qJD(1) * t428 + t463;
t238 = -pkin(8) * t298 - t350 * t425 + t330;
t604 = t238 * mrSges(3,1) + t83 * mrSges(4,1) - t237 * mrSges(3,2) - t84 * mrSges(4,2);
t573 = -t240 / 0.2e1;
t579 = -t158 / 0.2e1;
t580 = t157 / 0.2e1;
t581 = -t157 / 0.2e1;
t602 = -t28 * mrSges(6,1) + t24 * mrSges(7,1) + t29 * mrSges(6,2) - t25 * mrSges(7,3) + Ifges(6,6) * t580 + Ifges(7,6) * t581 + t573 * t659 + t579 * t661;
t599 = m(5) / 0.2e1;
t598 = m(6) / 0.2e1;
t597 = m(7) / 0.2e1;
t541 = Ifges(5,4) * t386;
t115 = Ifges(5,2) * t241 - t541 + t648;
t585 = -t115 / 0.2e1;
t578 = t158 / 0.2e1;
t572 = t240 / 0.2e1;
t571 = -t241 / 0.2e1;
t570 = t386 / 0.2e1;
t569 = -t386 / 0.2e1;
t567 = -t278 / 0.2e1;
t566 = -t279 / 0.2e1;
t550 = qJD(2) / 0.2e1;
t546 = mrSges(4,3) * t278;
t545 = mrSges(6,3) * t157;
t544 = mrSges(6,3) * t158;
t542 = Ifges(3,4) * t354;
t540 = Ifges(5,4) * t349;
t539 = Ifges(5,4) * t353;
t532 = t151 * mrSges(4,3);
t19 = -pkin(4) * t213 - t23;
t531 = t19 * t349;
t528 = t241 * mrSges(5,3);
t526 = t279 * Ifges(4,4);
t145 = -pkin(4) * t386 - pkin(10) * t241;
t55 = t348 * t145 + t352 * t85;
t523 = t230 * t349;
t233 = -t287 * t351 + t355 * t377;
t521 = t233 * t349;
t520 = t241 * t348;
t519 = t241 * t352;
t517 = t278 * t353;
t516 = t365 * t349;
t515 = t311 * t352;
t503 = t350 * t355;
t501 = t351 * t354;
t110 = -mrSges(7,2) * t157 + mrSges(7,3) * t240;
t111 = -mrSges(6,2) * t240 - t545;
t498 = t110 + t111;
t112 = mrSges(6,1) * t240 - t544;
t113 = -mrSges(7,1) * t240 + mrSges(7,2) * t158;
t497 = -t113 + t112;
t491 = qJD(4) * t348;
t490 = qJD(4) * t352;
t487 = qJD(5) * t349;
t477 = pkin(10) * t488;
t476 = pkin(10) * t486;
t466 = Ifges(5,5) * t128 + Ifges(5,6) * t129 + Ifges(5,3) * t213;
t456 = t342 * t484;
t455 = t348 * t484;
t452 = t512 / 0.2e1;
t448 = t354 * t550;
t445 = -t487 / 0.2e1;
t444 = t486 / 0.2e1;
t442 = t484 / 0.2e1;
t434 = -t218 * mrSges(4,1) + t219 * mrSges(4,2);
t104 = -t349 * t169 + t195 * t353;
t433 = -t231 * t349 - t353 * t509;
t299 = pkin(2) * t508 - t450;
t432 = -t299 * t351 + t355 * t344;
t427 = mrSges(3,3) * t461;
t426 = mrSges(3,3) * t460;
t419 = t676 * pkin(2);
t414 = -mrSges(4,1) * t365 - mrSges(4,2) * t286;
t413 = mrSges(5,1) * t257 + mrSges(5,2) * t258;
t407 = Ifges(5,1) * t353 - t540;
t402 = -Ifges(5,2) * t349 + t539;
t401 = -Ifges(6,2) * t348 + t536;
t400 = Ifges(6,2) * t352 + t537;
t397 = Ifges(5,5) * t353 - Ifges(5,6) * t349;
t394 = Ifges(7,3) * t348 + t533;
t393 = -Ifges(7,3) * t352 + t534;
t42 = t103 * t352 - t348 * t95;
t54 = t145 * t352 - t348 * t85;
t188 = t258 * t352 - t348 * t365;
t187 = t258 * t348 + t352 * t365;
t385 = -t299 * t355 - t351 * t344;
t49 = -t349 * t144 - t169 * t484 + t186 * t353 - t195 * t485;
t94 = pkin(4) * t365 - t104;
t305 = -t347 * t501 - t503;
t8 = t103 * t486 + t348 * t66 + t352 * t40 - t488 * t95;
t378 = (t354 * Ifges(3,2) + t543) * t346;
t373 = -pkin(3) * t232 - pkin(9) * t233 + t432;
t371 = t305 * pkin(2);
t370 = t335 * t346 * (Ifges(3,5) * t354 - Ifges(3,6) * t350);
t368 = -t348 * t487 + t352 * t484;
t367 = t349 * t486 + t455;
t366 = t230 * pkin(3) + pkin(9) * t231 + t419;
t41 = pkin(4) * t280 - t49;
t364 = -t231 * pkin(3) + t230 * pkin(9) + t385;
t360 = t233 * pkin(3) - t232 * pkin(9) + t371;
t331 = Ifges(3,4) * t460;
t318 = -pkin(4) - t392;
t313 = -t336 - t560;
t309 = -pkin(8) * t512 + t339;
t307 = (-mrSges(3,1) * t354 + mrSges(3,2) * t350) * t346;
t306 = -t347 * t502 + t499;
t304 = -t347 * t503 - t501;
t294 = t333 - t428;
t293 = t495 * qJD(1);
t292 = -pkin(8) * t461 + t332;
t291 = -mrSges(3,2) * t335 + t426;
t290 = mrSges(3,1) * t335 - t427;
t271 = Ifges(4,4) * t278;
t270 = t384 * t349;
t262 = Ifges(3,1) * t461 + Ifges(3,5) * t335 + t331;
t261 = t335 * Ifges(3,6) + qJD(1) * t378;
t252 = -t342 * t506 + t515;
t248 = -mrSges(4,2) * t335 + t546;
t245 = t353 * t441 - t515;
t244 = -qJ(6) * t353 + t635;
t236 = Ifges(5,4) * t241;
t207 = -mrSges(4,1) * t278 - mrSges(4,2) * t279;
t202 = -t232 * t349 - t353 * t511;
t190 = mrSges(4,1) * t334 - mrSges(4,3) * t219;
t189 = -mrSges(4,2) * t334 + mrSges(4,3) * t218;
t184 = -Ifges(4,1) * t279 + Ifges(4,5) * t335 + t271;
t183 = t278 * Ifges(4,2) + t335 * Ifges(4,6) - t526;
t165 = -mrSges(5,2) * t361 + t528;
t124 = t203 * t352 - t233 * t348;
t123 = t203 * t348 + t233 * t352;
t116 = -Ifges(5,1) * t386 + t236 + t651;
t114 = -Ifges(5,5) * t386 + Ifges(5,6) * t241 + t361 * Ifges(5,3);
t97 = mrSges(7,1) * t157 - mrSges(7,3) * t158;
t96 = pkin(5) * t158 + qJ(6) * t157;
t93 = -qJD(5) * t187 + t179 * t352 - t280 * t348;
t92 = qJD(5) * t188 + t179 * t348 + t280 * t352;
t91 = -mrSges(5,2) * t213 + mrSges(5,3) * t129;
t65 = -mrSges(5,1) * t129 + mrSges(5,2) * t128;
t51 = pkin(5) * t187 - qJ(6) * t188 + t94;
t45 = pkin(5) * t386 - t54;
t44 = -qJ(6) * t386 + t55;
t36 = -pkin(5) * t257 - t42;
t35 = qJ(6) * t257 + t639;
t20 = mrSges(7,1) * t63 - mrSges(7,3) * t62;
t10 = pkin(5) * t92 - qJ(6) * t93 - qJD(6) * t188 + t41;
t7 = -pkin(5) * t178 - t9;
t6 = qJ(6) * t178 + qJD(6) * t257 + t8;
t5 = pkin(5) * t63 - qJ(6) * t62 - qJD(6) * t158 + t19;
t11 = [(-m(4) * t385 - m(5) * t364 - t304 * mrSges(3,1) + t231 * mrSges(4,1) + t199 * mrSges(5,1) + mrSges(2,2) * t355 + t676 * mrSges(3,2) + t617 * t351 + t612 * t509 + t670 * (-pkin(4) * t199 + pkin(10) * t433 + t364) + t665 * t230 + t431 * t672 - t424 * t119 + t619 * t433) * g(1) + (-m(4) * t432 - m(5) * t373 - t306 * mrSges(3,1) + t232 * mrSges(4,1) - t203 * mrSges(5,1) + t351 * mrSges(2,2) - t305 * mrSges(3,2) - t617 * t355 + t612 * t511 + t670 * (t203 * pkin(4) + pkin(10) * t202 + t373) - t665 * t233 - t431 * t124 + t424 * t123 + t619 * t202) * g(2) + (Ifges(4,1) * t281 + Ifges(4,4) * t280) * t566 + (Ifges(5,1) * t179 - Ifges(5,4) * t178 - Ifges(5,5) * t280) * t569 + t361 * (Ifges(5,5) * t179 - Ifges(5,6) * t178 - Ifges(5,3) * t280) / 0.2e1 + t335 * (Ifges(4,5) * t281 + Ifges(4,6) * t280) / 0.2e1 + t301 * (-mrSges(4,1) * t280 + mrSges(4,2) * t281) + t241 * (Ifges(5,4) * t179 - Ifges(5,2) * t178 - Ifges(5,6) * t280) / 0.2e1 + t278 * (Ifges(4,4) * t281 + Ifges(4,2) * t280) / 0.2e1 + m(6) * (t19 * t94 + t28 * t9 + t29 * t8 + t3 * t639 + t4 * t42 + t41 * t73) + t639 * t33 - (Ifges(5,6) * t582 + Ifges(5,5) * t583 + Ifges(5,3) * t574 + t466 / 0.2e1 - Ifges(4,2) * t218 - Ifges(4,4) * t219 - t84 * mrSges(4,3) - t620 * Ifges(4,6) + t610) * t365 + ((Ifges(3,5) * t350 + Ifges(3,6) * t354) * t564 + t262 * t448 + t298 * (Ifges(3,1) * t350 + t542) / 0.2e1) * t346 + (Ifges(7,5) * t93 + Ifges(7,6) * t178) * t580 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t601 + t237 * t495 + t238 * t309 - t292 * t295 + t293 * t294) + (-m(4) * t150 + t675) * t143 - t280 * t667 + t297 * t378 / 0.2e1 + m(5) * (t104 * t23 + t105 * t22 + t168 * t78 + t48 * t86 + t49 * t85) + ((-qJD(2) * t292 + t237) * t510 - t238 * t512 - t293 * t459 + t297 * t495 - t298 * t309) * mrSges(3,3) + m(7) * (t1 * t35 + t10 * t34 + t2 * t36 + t24 * t7 + t25 * t6 + t5 * t51) + t280 * t532 + (Ifges(3,4) * t298 + Ifges(3,2) * t297 + Ifges(3,6) * t334) * t510 / 0.2e1 + t178 * t585 + (Ifges(3,1) * t298 + Ifges(3,4) * t297 + Ifges(3,5) * t334) * t452 + (t658 / 0.2e1 - t687 - t22 * mrSges(5,3) - t671) * t257 - (-mrSges(3,1) * t297 + mrSges(3,2) * t298) * t560 + (t178 * t661 + t662 * t92 + t663 * t93) * t578 - t307 * t474 - t261 * t459 / 0.2e1 + t207 * t429 + t654 * t93 / 0.2e1 + t655 * t178 / 0.2e1 + (mrSges(6,2) * t19 - mrSges(7,3) * t5 + Ifges(6,4) * t592 + Ifges(7,5) * t591 + t674 + t685) * t188 + m(4) * (t144 * t151 + t180 * t83 + t181 * t84 + t259 * t313 + t301 * t429) + t313 * t434 + (Ifges(6,4) * t93 + Ifges(6,6) * t178) * t581 + t280 * t666 - t495 * mrSges(3,2) * t334 + (-Ifges(6,2) * t581 + Ifges(7,3) * t580 - t605 + t623) * t92 + t621 * t483 - (-t83 * mrSges(4,3) + Ifges(4,1) * t219 + Ifges(4,4) * t218 + Ifges(4,5) * t620) * t286 + t259 * t414 + t78 * t413 + Ifges(2,3) * qJDD(1) + t51 * t20 + t42 * t31 + t35 * t30 + t36 * t32 + (t178 * t25 - t34 * t93) * mrSges(7,3) - t150 * t281 * mrSges(4,3) + (-t178 * t29 + t73 * t93) * mrSges(6,2) + t601 * qJD(1) * (-Ifges(3,2) * t350 + t542) * t448 + t294 * t291 - t295 * t290 - t280 * t114 / 0.2e1 + t280 * t183 / 0.2e1 + t281 * t184 / 0.2e1 + t144 * t248 + t309 * mrSges(3,1) * t334 + t94 * t21 + t10 * t97 + t41 * t98 + (t609 + t681) * t258 + t104 * t90 + t105 * t91 + t6 * t110 + t8 * t111 + t9 * t112 + t7 * t113 + (t650 / 0.2e1 + t649 / 0.2e1 + t652 / 0.2e1 + t653 / 0.2e1 + t683 * t564 + t604 + t628 / 0.2e1) * t347 + (mrSges(6,1) * t19 + mrSges(7,1) * t5 - Ifges(6,2) * t592 + Ifges(7,3) * t591 - t584 * t684 + t593 * t662 + t608) * t187 + (t178 * t659 + t661 * t93 - t684 * t92) * t572 + (-t178 * t86 - t179 * t85 - t23 * t258) * mrSges(5,3) + t48 * t165 + t49 * t166 + t168 * t65 + t370 * t550 + t28 * (mrSges(6,1) * t178 - mrSges(6,3) * t93) + t24 * (-mrSges(7,1) * t178 + mrSges(7,2) * t93) + t179 * t116 / 0.2e1 + t140 * (mrSges(5,1) * t178 + mrSges(5,2) * t179) + t181 * t189 + t180 * t190; (-m(4) * t419 - m(5) * t366 - mrSges(3,1) * t676 - mrSges(3,2) * t304 - t664 * t523 + t670 * (pkin(10) * t523 + t230 * t551 + t366) + t647 * t230 + t665 * t231 - t431 * (t230 * t500 + t231 * t348) + t424 * (t230 * t506 - t231 * t352)) * g(2) + (t241 * t402 + t361 * t397 - t386 * t407) * qJD(4) / 0.2e1 + (t252 * t4 + t635 * t3 + (t484 * t73 + t531) * t342 - t73 * t87 + t640 * t29 + t641 * t28) * m(6) + t635 * t33 + (t261 * t452 - t370 / 0.2e1 - t621 * qJD(1)) * qJD(1) + t674 * t349 * t352 + (-m(4) * t336 - m(5) * t418 + t286 * mrSges(5,3) - t365 * t412 + t307 + t414 - t664 * t516 + t670 * (pkin(10) * t516 + t365 * t551 + t418) - t431 * (-t286 * t348 + t365 * t500) + t424 * (t286 * t352 + t365 * t506)) * g(3) + (-t457 - t102) * t165 + (t426 - t291) * t292 + (-t456 - t101) * t166 + (Ifges(4,2) * t279 + t184 + t271) * t567 + t279 * t667 - (-Ifges(3,2) * t461 + t262 + t331) * t460 / 0.2e1 + (Ifges(4,1) * t278 + t114 + t526) * t279 / 0.2e1 + t539 * t583 - t207 * t430 + t189 * t558 + t183 * t566 + (t456 - t87) * t98 + t352 * t70 * t445 + t348 * t67 * t442 + t540 * t582 - t279 * t532 + t485 * t585 + t455 * t588 + t410 * t531 - t622 * t140 * (mrSges(5,1) * t349 + mrSges(5,2) * t353) + (t427 + t290) * t293 + t656 * t342 * t349 + (-t73 * mrSges(6,2) - t24 * mrSges(7,2) + t28 * mrSges(6,3) + t34 * mrSges(7,3) + Ifges(6,4) * t580 + Ifges(7,5) * t581 + t573 * t661 + t579 * t663) * t197 + (-t632 * t487 + (t349 * t659 + t353 * t631) * qJD(4)) * t572 + (-t629 * t487 + (t349 * t661 + t353 * t630) * qJD(4)) * t578 + (t342 * t91 + (t394 * t580 + t401 * t581) * qJD(4) + t671) * t353 + (t150 * t172 - t151 * t173 - t301 * t430 + (t345 * t84 + t525 * t83) * pkin(2)) * m(4) + (-Ifges(5,5) * t279 + t278 * t407) * t570 + (-Ifges(5,6) * t279 + t278 * t402) * t571 + t29 * (-mrSges(6,2) * t485 - mrSges(6,3) * t367) + t24 * (-mrSges(7,1) * t485 + mrSges(7,2) * t368) + t25 * (-mrSges(7,2) * t367 + mrSges(7,3) * t485) + t28 * (mrSges(6,1) * t485 - mrSges(6,3) * t368) + (t485 / 0.2e1 - t518 / 0.2e1) * t655 + (-m(4) * t371 - m(5) * t360 - mrSges(3,1) * t305 + mrSges(3,2) * t306 - t664 * t521 + t670 * (pkin(10) * t521 + t233 * t551 + t360) + t647 * t233 - t665 * t232 - t431 * (-t232 * t348 + t233 * t500) + t424 * (t232 * t352 + t233 * t506)) * g(1) - t658 * t353 / 0.2e1 + t604 + (-t393 * t580 - t400 * t581) * t487 + t636 * t172 + t640 * t111 + t641 * t112 + t642 * t113 + t643 * t110 + t645 * t97 + (t1 * t244 + t2 * t245 + t24 * t642 + t25 * t643 + t270 * t5 + t34 * t645) * m(7) - t279 * t666 + (t394 * t591 + t401 * t592 + t5 * t408 + t67 * t444 + t630 * t593 + t631 * t584 + (Ifges(6,6) * t581 + Ifges(7,6) * t580) * qJD(4) + t609) * t349 + ((-t485 + t518) * t86 + (-t484 + t517) * t85 + t625) * mrSges(5,3) + (t343 * t78 + ((-t349 * t86 - t353 * t85) * qJD(4) + t625) * t342 - t101 * t85 - t102 * t86 - t140 * t172) * m(5) - t78 * t412 - t361 * (-Ifges(5,3) * t279 + t278 * t397) / 0.2e1 + (t115 / 0.2e1 + t602) * t518 + t608 * t348 * t349 + (-t517 / 0.2e1 + t442) * t116 + t73 * (mrSges(6,1) * t367 + mrSges(6,2) * t368) + t34 * (mrSges(7,1) * t367 - mrSges(7,3) * t368) - t335 * (Ifges(4,5) * t278 + Ifges(4,6) * t279) / 0.2e1 + t343 * t65 - t301 * (-mrSges(4,1) * t279 + mrSges(4,2) * t278) + t270 * t20 + t252 * t31 + t244 * t30 + t245 * t32 - t173 * t248 + t628 + t654 * (t348 * t445 + t352 * t442 - t197 / 0.2e1) + (-Ifges(6,2) * t580 + Ifges(7,3) * t581 - t573 * t684 + t579 * t662 + t605 + t624) * t196 + t150 * t546 + t190 * t462; -t278 * t248 - t498 * t197 + t497 * t196 + (-t278 * t165 - t20 + (-t348 * t497 + t352 * t498 + t165) * qJD(4) - t656) * t353 + (t91 + t615 + t614 + (-t348 * t498 - t352 * t497) * qJD(5) + t622 * (-t97 + t646)) * t349 - m(7) * (t196 * t24 + t197 * t25) - m(6) * (-t196 * t28 + t197 * t29) + 0.2e1 * ((t24 * t491 + t25 * t490 - t5) * t597 + (-t28 * t491 + t29 * t490 - t19) * t598 + m(5) * t86 * t567 + (qJD(4) * t86 + t23) * t599) * t353 + 0.2e1 * ((qJD(4) * t34 + t677) * t597 + (qJD(4) * t73 + t678) * t598 + (-qJD(4) * t85 + t22) * t599 + (t85 * t599 - m(7) * t34 / 0.2e1 - t669 / 0.2e1) * t278) * t349 + t434 + (-t347 * g(3) + (-g(1) * t351 + g(2) * t355) * t346) * (m(4) + m(5) - t670) + t675 * t279 + (-t150 * t279 - t151 * t278 + t259) * m(4); -(-Ifges(5,2) * t571 + t648 / 0.2e1 - t140 * mrSges(5,1) + t602) * t386 + (t158 * t630 + t240 * t631) * qJD(5) / 0.2e1 + (t444 - t519 / 0.2e1) * t654 + (-t24 * t519 + t25 * t520 + t607 + t677) * mrSges(7,2) + (t28 * t519 + t29 * t520 + t607 + t678) * mrSges(6,3) + (-t28 * t54 - t29 * t55 - pkin(4) * t19 + ((-t28 * t352 - t29 * t348) * qJD(5) + t415) * pkin(10)) * m(6) + (t476 - t45) * t113 + (-t476 - t54) * t112 + (-t477 - t55) * t111 + (-t477 - t44) * t110 + pkin(10) * t614 + pkin(10) * t615 + t673 * t352 + (-t401 / 0.2e1 + t394 / 0.2e1) * qJD(5) * t157 + t348 * t682 + (t401 * t580 + t394 * t581 + Ifges(5,1) * t570 - t651 / 0.2e1 - t140 * mrSges(5,2) + t630 * t579 + t631 * t573 - t633) * t241 + t115 * t569 + (-t527 + t646 - t669) * t86 + (t634 * t257 + t438 * t670 + t413) * g(3) + (mrSges(5,2) * t199 + t670 * (pkin(4) * t433 + pkin(10) * t199) - t618 * t433) * g(2) + (mrSges(5,2) * t203 + t670 * (-t202 * pkin(4) + pkin(10) * t203) + t618 * t202) * g(1) + (t541 + t655) * t570 + t466 + (t528 - t165) * t85 + (t236 + t116) * t571 + t610 + (-t24 * t45 - t25 * t44 + t318 * t5 + ((t24 * t352 - t25 * t348) * qJD(5) + t416) * pkin(10) + t644 * t34) * m(7) + t644 * t97 + t629 * t593 + t632 * t584 + t633 * qJD(5) + t623 * t488 + t624 * t520 - t19 * t411 + t5 * t409 - pkin(4) * t21 + t318 * t20 + t393 * t591 + t400 * t592; (-t157 * t663 + t152 - t538 + t67) * t579 + (-pkin(5) * t2 + qJ(6) * t1 - t24 * t29 + t25 * t638 - t34 * t96) * m(7) + (-t157 * t661 - t158 * t684) * t573 + (t187 * t431 + t188 * t424) * g(3) + (t123 * t431 + t124 * t424) * g(1) + (t497 + t544) * t29 + (-t498 - t545) * t28 + (Ifges(7,3) * t158 - t535) * t581 + (t431 * t119 + t424 * t672) * g(2) + t70 * t578 + t606 + (t157 * t24 + t158 * t25) * mrSges(7,2) + (-Ifges(6,2) * t158 - t153 + t654) * t580 + qJ(6) * t30 - pkin(5) * t32 - t96 * t97 + qJD(6) * t110 - t34 * (mrSges(7,1) * t158 + mrSges(7,3) * t157) - t73 * (mrSges(6,1) * t158 - mrSges(6,2) * t157) + t658; -t240 * t110 + t158 * t97 + (-g(1) * t123 - g(2) * t119 - g(3) * t187 + t34 * t158 - t25 * t240 + t2) * m(7) + t32;];
tau  = t11;
