% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR9
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:22
% EndTime: 2019-03-09 14:09:49
% DurationCPUTime: 52.98s
% Computational Cost: add. (38229->1029), mult. (92820->1397), div. (0->0), fcn. (77757->18), ass. (0->479)
t421 = cos(pkin(6));
t542 = qJD(1) * t421;
t401 = qJD(2) + t542;
t418 = sin(pkin(12));
t420 = cos(pkin(12));
t426 = sin(qJ(2));
t419 = sin(pkin(6));
t543 = qJD(1) * t419;
t508 = t426 * t543;
t319 = t401 * t420 - t418 * t508;
t320 = t401 * t418 + t420 * t508;
t425 = sin(qJ(4));
t430 = cos(qJ(4));
t240 = t319 * t425 + t320 * t430;
t424 = sin(qJ(5));
t429 = cos(qJ(5));
t495 = t430 * t319 - t320 * t425;
t717 = t429 * t240 + t424 * t495;
t641 = -t717 / 0.2e1;
t431 = cos(qJ(2));
t507 = t431 * t543;
t383 = qJD(4) - t507;
t376 = qJD(5) + t383;
t629 = -t376 / 0.2e1;
t718 = -t240 * t424 + t429 * t495;
t643 = -t718 / 0.2e1;
t740 = Ifges(6,4) * t643 + Ifges(6,5) * t629;
t745 = Ifges(6,1) * t641 + t740;
t367 = t418 * t430 + t420 * t425;
t556 = t419 * t431;
t440 = t367 * t556;
t298 = qJD(1) * t440;
t354 = t367 * qJD(4);
t733 = t298 - t354;
t628 = t376 / 0.2e1;
t640 = t717 / 0.2e1;
t743 = Ifges(6,4) * t640 + Ifges(6,6) * t628;
t160 = -t718 + qJD(6);
t645 = -t160 / 0.2e1;
t423 = sin(qJ(6));
t428 = cos(qJ(6));
t143 = t376 * t423 + t428 * t717;
t651 = -t143 / 0.2e1;
t142 = t376 * t428 - t423 * t717;
t653 = -t142 / 0.2e1;
t742 = Ifges(7,5) * t651 + Ifges(7,6) * t653 + Ifges(7,3) * t645;
t737 = mrSges(6,2) - mrSges(7,3);
t642 = t718 / 0.2e1;
t741 = Ifges(6,2) * t642;
t480 = -mrSges(7,1) * t428 + mrSges(7,2) * t423;
t736 = mrSges(6,1) - t480;
t674 = m(7) * pkin(5);
t739 = -t674 - t736;
t673 = m(7) * pkin(11);
t493 = -t673 + t737;
t473 = pkin(2) * t426 - qJ(3) * t431;
t341 = t473 * t543;
t527 = pkin(1) * t542;
t344 = -pkin(8) * t508 + t431 * t527;
t248 = t420 * t341 - t418 * t344;
t554 = t420 * t431;
t443 = (pkin(3) * t426 - pkin(9) * t554) * t419;
t209 = qJD(1) * t443 + t248;
t249 = t418 * t341 + t420 * t344;
t488 = t418 * t507;
t228 = -pkin(9) * t488 + t249;
t149 = t430 * t209 - t228 * t425;
t366 = -t418 * t425 + t420 * t430;
t439 = t366 * t556;
t299 = qJD(1) * t439;
t130 = pkin(4) * t508 - pkin(10) * t299 + t149;
t150 = t425 * t209 + t430 * t228;
t133 = -pkin(10) * t298 + t150;
t607 = pkin(9) + qJ(3);
t378 = t607 * t418;
t379 = t607 * t420;
t534 = qJD(4) * t430;
t537 = qJD(3) * t420;
t538 = qJD(3) * t418;
t245 = -t378 * t534 + t430 * t537 + (-qJD(4) * t379 - t538) * t425;
t210 = -pkin(10) * t354 + t245;
t297 = -t425 * t378 + t430 * t379;
t246 = -t367 * qJD(3) - qJD(4) * t297;
t353 = t366 * qJD(4);
t432 = -pkin(10) * t353 + t246;
t296 = -t430 * t378 - t379 * t425;
t251 = -pkin(10) * t367 + t296;
t252 = pkin(10) * t366 + t297;
t463 = t429 * t251 - t252 * t424;
t706 = qJD(5) * t463 + (-t133 + t210) * t429 + (-t130 + t432) * t424;
t393 = pkin(8) * t507;
t345 = t426 * t527 + t393;
t287 = pkin(3) * t488 + t345;
t686 = -pkin(4) * t733 - t287;
t293 = -t401 * pkin(2) + qJD(3) - t344;
t236 = -t319 * pkin(3) + t293;
t168 = -pkin(4) * t495 + t236;
t684 = -t168 * mrSges(6,2) + t745;
t303 = qJ(3) * t401 + t345;
t328 = (-pkin(2) * t431 - qJ(3) * t426 - pkin(1)) * t419;
t311 = qJD(1) * t328;
t214 = -t418 * t303 + t420 * t311;
t177 = -pkin(3) * t507 - t320 * pkin(9) + t214;
t215 = t420 * t303 + t418 * t311;
t183 = pkin(9) * t319 + t215;
t124 = t177 * t425 + t183 * t430;
t109 = pkin(10) * t495 + t124;
t550 = t429 * t109;
t123 = t430 * t177 - t183 * t425;
t108 = -pkin(10) * t240 + t123;
t99 = pkin(4) * t383 + t108;
t57 = t424 * t99 + t550;
t55 = pkin(11) * t376 + t57;
t87 = -pkin(5) * t718 - pkin(11) * t717 + t168;
t25 = -t423 * t55 + t428 * t87;
t26 = t423 * t87 + t428 * t55;
t738 = t742 + t743;
t732 = -t168 * mrSges(6,1) - t25 * mrSges(7,1) + t26 * mrSges(7,2) + t57 * mrSges(6,3) + t738 + t741;
t735 = -pkin(11) * t508 + t706;
t461 = t429 * t366 - t367 * t424;
t193 = qJD(5) * t461 + t353 * t429 - t354 * t424;
t276 = t366 * t424 + t367 * t429;
t194 = qJD(5) * t276 + t353 * t424 + t429 * t354;
t211 = t429 * t298 + t299 * t424;
t212 = -t298 * t424 + t299 * t429;
t734 = t686 + (-t193 + t212) * pkin(11) + (t194 - t211) * pkin(5);
t539 = qJD(2) * t431;
t349 = (qJD(1) * t539 + qJDD(1) * t426) * t419;
t529 = qJDD(1) * t421;
t400 = qJDD(2) + t529;
t279 = -t349 * t418 + t400 * t420;
t280 = t349 * t420 + t400 * t418;
t147 = qJD(4) * t495 + t279 * t425 + t280 * t430;
t540 = qJD(2) * t419;
t506 = t426 * t540;
t724 = -qJD(1) * t506 + qJDD(1) * t556;
t337 = qJDD(4) - t724;
t620 = pkin(1) * t421;
t526 = qJD(2) * t620;
t490 = qJD(1) * t526;
t521 = pkin(1) * t529;
t257 = pkin(8) * t724 + t426 * t521 + t431 * t490;
t225 = qJ(3) * t400 + qJD(3) * t401 + t257;
t536 = qJD(3) * t426;
t234 = -pkin(2) * t724 - qJ(3) * t349 + (-pkin(1) * qJDD(1) - qJD(1) * t536) * t419;
t154 = -t225 * t418 + t420 * t234;
t127 = -pkin(3) * t724 - pkin(9) * t280 + t154;
t155 = t420 * t225 + t418 * t234;
t134 = pkin(9) * t279 + t155;
t52 = -qJD(4) * t124 + t430 * t127 - t134 * t425;
t35 = pkin(4) * t337 - pkin(10) * t147 + t52;
t148 = -qJD(4) * t240 + t279 * t430 - t280 * t425;
t535 = qJD(4) * t425;
t51 = t425 * t127 + t430 * t134 + t177 * t534 - t183 * t535;
t39 = pkin(10) * t148 + t51;
t532 = qJD(5) * t429;
t533 = qJD(5) * t424;
t10 = -t109 * t533 + t424 * t35 + t429 * t39 + t99 * t532;
t558 = t419 * t426;
t402 = pkin(8) * t558;
t258 = -qJD(2) * t393 - qJDD(1) * t402 - t426 * t490 + t431 * t521;
t241 = -t400 * pkin(2) + qJDD(3) - t258;
t182 = -t279 * pkin(3) + t241;
t114 = -t148 * pkin(4) + t182;
t326 = qJDD(5) + t337;
t80 = qJD(5) * t718 + t147 * t429 + t148 * t424;
t47 = qJD(6) * t142 + t326 * t423 + t428 * t80;
t48 = -qJD(6) * t143 + t326 * t428 - t423 * t80;
t81 = -qJD(5) * t717 - t147 * t424 + t148 * t429;
t75 = qJDD(6) - t81;
t14 = Ifges(7,5) * t47 + Ifges(7,6) * t48 + Ifges(7,3) * t75;
t633 = t326 / 0.2e1;
t660 = t81 / 0.2e1;
t661 = t80 / 0.2e1;
t665 = t75 / 0.2e1;
t666 = t48 / 0.2e1;
t667 = t47 / 0.2e1;
t22 = -t81 * pkin(5) - t80 * pkin(11) + t114;
t7 = pkin(11) * t326 + t10;
t2 = qJD(6) * t25 + t22 * t423 + t428 * t7;
t3 = -qJD(6) * t26 + t22 * t428 - t423 * t7;
t683 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t731 = t683 + mrSges(6,1) * t114 - mrSges(6,3) * t10 + Ifges(7,5) * t667 + Ifges(7,6) * t666 + Ifges(7,3) * t665 + t14 / 0.2e1 + (-t633 - t326 / 0.2e1) * Ifges(6,6) + (-t660 - t81 / 0.2e1) * Ifges(6,2) + (-t661 - t80 / 0.2e1) * Ifges(6,4);
t644 = t160 / 0.2e1;
t650 = t143 / 0.2e1;
t652 = t142 / 0.2e1;
t729 = Ifges(7,5) * t650 + Ifges(7,6) * t652 + Ifges(7,3) * t644 - t732 - t741 - t743;
t568 = t109 * t424;
t56 = t429 * t99 - t568;
t610 = t56 * mrSges(6,3);
t728 = Ifges(6,1) * t640 + Ifges(6,4) * t642 + Ifges(6,5) * t628 - t610 - t684;
t595 = Ifges(5,4) * t240;
t157 = Ifges(5,2) * t495 + t383 * Ifges(5,6) + t595;
t647 = t157 / 0.2e1;
t11 = -qJD(5) * t57 + t35 * t429 - t39 * t424;
t727 = t11 * mrSges(6,1) - t10 * mrSges(6,2);
t726 = t257 * mrSges(3,2);
t725 = -m(4) * qJ(3) - m(5) * t607 - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t153 = mrSges(6,1) * t376 - mrSges(6,3) * t717;
t95 = -mrSges(7,1) * t142 + mrSges(7,2) * t143;
t571 = -t95 + t153;
t409 = pkin(3) * t420 + pkin(2);
t482 = -mrSges(4,1) * t420 + mrSges(4,2) * t418;
t723 = -m(4) * pkin(2) - m(5) * t409 + t482;
t622 = cos(qJ(1));
t509 = t622 * t431;
t427 = sin(qJ(1));
t552 = t426 * t427;
t358 = -t421 * t552 + t509;
t417 = pkin(12) + qJ(4);
t412 = sin(t417);
t413 = cos(t417);
t557 = t419 * t427;
t281 = -t358 * t412 + t413 * t557;
t722 = -t412 * t558 + t413 * t421;
t721 = -Ifges(6,4) * t641 - Ifges(6,2) * t643 - Ifges(6,6) * t629 + t732 + t742;
t720 = t52 * mrSges(5,1) - t51 * mrSges(5,2);
t414 = qJ(5) + t417;
t407 = sin(t414);
t408 = cos(t414);
t719 = -t413 * mrSges(5,1) + t412 * mrSges(5,2) + t407 * t493 + t408 * t739 + t723;
t479 = mrSges(7,1) * t423 + mrSges(7,2) * t428;
t54 = -pkin(5) * t376 - t56;
t454 = t54 * t479;
t141 = Ifges(7,4) * t142;
t78 = Ifges(7,1) * t143 + Ifges(7,5) * t160 + t141;
t572 = t428 * t78;
t624 = t423 / 0.2e1;
t593 = Ifges(7,4) * t143;
t77 = Ifges(7,2) * t142 + Ifges(7,6) * t160 + t593;
t716 = -t572 / 0.2e1 + t77 * t624 + t610 - t454;
t714 = t114 * mrSges(6,2) - mrSges(6,3) * t11 + 0.2e1 * Ifges(6,1) * t661 + 0.2e1 * Ifges(6,4) * t660 + 0.2e1 * Ifges(6,5) * t633;
t711 = -m(4) - m(5);
t710 = m(7) + m(6);
t649 = t147 / 0.2e1;
t648 = t148 / 0.2e1;
t632 = t337 / 0.2e1;
t175 = t251 * t424 + t252 * t429;
t324 = -pkin(4) * t366 - t409;
t176 = -pkin(5) * t461 - pkin(11) * t276 + t324;
t120 = t175 * t428 + t176 * t423;
t709 = -qJD(6) * t120 - t423 * t735 + t428 * t734;
t119 = -t175 * t423 + t176 * t428;
t708 = qJD(6) * t119 + t423 * t734 + t428 * t735;
t361 = pkin(8) * t556 + t426 * t620;
t327 = qJ(3) * t421 + t361;
t243 = -t418 * t327 + t420 * t328;
t351 = t418 * t421 + t420 * t558;
t195 = -pkin(3) * t556 - t351 * pkin(9) + t243;
t244 = t420 * t327 + t418 * t328;
t350 = -t418 * t558 + t420 * t421;
t206 = pkin(9) * t350 + t244;
t137 = t430 * t195 - t425 * t206;
t254 = t350 * t425 + t351 * t430;
t117 = -pkin(4) * t556 - t254 * pkin(10) + t137;
t138 = t425 * t195 + t430 * t206;
t253 = t350 * t430 - t351 * t425;
t121 = pkin(10) * t253 + t138;
t703 = t424 * t117 + t429 * t121;
t702 = -t149 + t246;
t701 = -t150 + t245;
t152 = -mrSges(6,2) * t376 + mrSges(6,3) * t718;
t100 = -mrSges(7,2) * t160 + mrSges(7,3) * t142;
t101 = mrSges(7,1) * t160 - mrSges(7,3) * t143;
t467 = t100 * t428 - t101 * t423;
t700 = -t152 - t467;
t191 = -t212 * t423 + t428 * t508;
t530 = qJD(6) * t428;
t453 = t423 * t193 + t276 * t530;
t699 = t191 + t453;
t192 = t212 * t428 + t423 * t508;
t531 = qJD(6) * t423;
t452 = -t428 * t193 + t276 * t531;
t698 = t192 + t452;
t492 = mrSges(3,3) * t508;
t697 = -mrSges(3,1) * t401 - mrSges(4,1) * t319 + mrSges(4,2) * t320 + t492;
t695 = t299 - t353;
t438 = mrSges(3,2) + t725;
t694 = -t438 + t479;
t316 = -t407 * t558 + t408 * t421;
t317 = t407 * t421 + t408 * t558;
t692 = -t736 * t316 + t317 * t737;
t272 = t358 * t407 - t408 * t557;
t273 = t358 * t408 + t407 * t557;
t691 = t736 * t272 + t273 * t737;
t510 = t622 * t426;
t551 = t427 * t431;
t356 = t421 * t510 + t551;
t511 = t419 * t622;
t268 = -t356 * t407 - t408 * t511;
t269 = t356 * t408 - t407 * t511;
t690 = -t736 * t268 + t269 * t737;
t688 = -t154 * t418 + t155 * t420;
t23 = mrSges(7,1) * t75 - mrSges(7,3) * t47;
t24 = -mrSges(7,2) * t75 + mrSges(7,3) * t48;
t687 = -t423 * t23 + t428 * t24;
t599 = Ifges(3,4) * t426;
t685 = pkin(1) * (mrSges(3,1) * t426 + mrSges(3,2) * t431) - t426 * (Ifges(3,1) * t431 - t599) / 0.2e1;
t199 = qJD(2) * t439 + qJD(4) * t253;
t305 = (qJD(2) * t473 - t536) * t419;
t346 = -pkin(8) * t506 + t431 * t526;
t315 = qJD(3) * t421 + t346;
t226 = t420 * t305 - t418 * t315;
t188 = qJD(2) * t443 + t226;
t227 = t418 * t305 + t420 * t315;
t505 = t419 * t539;
t487 = t418 * t505;
t202 = -pkin(9) * t487 + t227;
t91 = -qJD(4) * t138 + t430 * t188 - t202 * t425;
t70 = pkin(4) * t506 - pkin(10) * t199 + t91;
t200 = -qJD(2) * t440 - qJD(4) * t254;
t90 = t425 * t188 + t195 * t534 + t430 * t202 - t206 * t535;
t82 = pkin(10) * t200 + t90;
t20 = -qJD(5) * t703 - t424 * t82 + t429 * t70;
t113 = pkin(5) * t717 - pkin(11) * t718;
t681 = mrSges(3,1) - t719;
t678 = t684 + t745;
t587 = Ifges(3,6) * t401;
t677 = t124 * mrSges(5,2) + t57 * mrSges(6,2) + t587 / 0.2e1 + (t431 * Ifges(3,2) + t599) * t543 / 0.2e1 - t123 * mrSges(5,1) - t56 * mrSges(6,1);
t675 = t419 ^ 2;
t15 = t47 * Ifges(7,4) + t48 * Ifges(7,2) + t75 * Ifges(7,6);
t671 = t15 / 0.2e1;
t16 = t47 * Ifges(7,1) + t48 * Ifges(7,4) + t75 * Ifges(7,5);
t670 = t16 / 0.2e1;
t659 = Ifges(5,4) * t649 + Ifges(5,2) * t648 + Ifges(5,6) * t632;
t658 = Ifges(5,1) * t649 + Ifges(5,4) * t648 + Ifges(5,5) * t632;
t235 = Ifges(5,4) * t495;
t158 = t240 * Ifges(5,1) + t383 * Ifges(5,5) + t235;
t646 = t158 / 0.2e1;
t639 = -t495 / 0.2e1;
t638 = t495 / 0.2e1;
t637 = -t240 / 0.2e1;
t636 = t240 / 0.2e1;
t635 = t279 / 0.2e1;
t634 = t280 / 0.2e1;
t631 = t350 / 0.2e1;
t630 = t351 / 0.2e1;
t627 = -t383 / 0.2e1;
t626 = t383 / 0.2e1;
t625 = t421 / 0.2e1;
t623 = t431 / 0.2e1;
t619 = pkin(1) * t431;
t618 = pkin(4) * t240;
t616 = pkin(4) * t424;
t615 = pkin(4) * t429;
t605 = mrSges(4,2) * t420;
t603 = mrSges(5,3) * t495;
t602 = mrSges(5,3) * t240;
t601 = mrSges(7,3) * t423;
t600 = mrSges(7,3) * t428;
t598 = Ifges(3,4) * t431;
t597 = Ifges(4,4) * t418;
t596 = Ifges(4,4) * t420;
t592 = Ifges(7,4) * t423;
t591 = Ifges(7,4) * t428;
t590 = Ifges(4,5) * t320;
t586 = Ifges(4,6) * t319;
t583 = Ifges(4,3) * t426;
t581 = t718 * Ifges(6,6);
t580 = t717 * Ifges(6,5);
t579 = t495 * Ifges(5,6);
t578 = t240 * Ifges(5,5);
t577 = t376 * Ifges(6,3);
t576 = t383 * Ifges(5,3);
t575 = t401 * Ifges(3,5);
t564 = t276 * t423;
t563 = t276 * t428;
t559 = t418 * t431;
t555 = t420 * (Ifges(4,1) * t320 + Ifges(4,4) * t319 - Ifges(4,5) * t507);
t347 = pkin(8) * t505 + t426 * t526;
t544 = t622 * pkin(1) + pkin(8) * t557;
t528 = Ifges(6,5) * t80 + Ifges(6,6) * t81 + Ifges(6,3) * t326;
t518 = t423 * t556;
t517 = t428 * t556;
t514 = t572 / 0.2e1;
t513 = Ifges(5,5) * t147 + Ifges(5,6) * t148 + Ifges(5,3) * t337;
t512 = Ifges(3,5) * t349 + Ifges(3,6) * t724 + Ifges(3,3) * t400;
t288 = pkin(3) * t487 + t347;
t31 = -t81 * mrSges(6,1) + t80 * mrSges(6,2);
t500 = -t543 / 0.2e1;
t499 = -t531 / 0.2e1;
t498 = -t427 * pkin(1) + pkin(8) * t511;
t497 = -t272 * pkin(5) + pkin(11) * t273;
t198 = -t279 * mrSges(4,1) + t280 * mrSges(4,2);
t97 = -t148 * mrSges(5,1) + t147 * mrSges(5,2);
t494 = -t356 * t413 + t412 * t511;
t491 = mrSges(3,3) * t507;
t489 = t418 * t511;
t484 = t431 * t500;
t483 = t281 * pkin(4);
t478 = Ifges(4,1) * t420 - t597;
t477 = Ifges(7,1) * t428 - t592;
t476 = -Ifges(4,2) * t418 + t596;
t475 = -Ifges(7,2) * t423 + t591;
t474 = Ifges(7,5) * t428 - Ifges(7,6) * t423;
t472 = t25 * t428 + t26 * t423;
t471 = -t25 * t423 + t26 * t428;
t63 = -pkin(11) * t556 + t703;
t181 = t253 * t424 + t254 * t429;
t333 = t402 + (-pkin(2) - t619) * t421;
t259 = -t350 * pkin(3) + t333;
t186 = -t253 * pkin(4) + t259;
t462 = t429 * t253 - t254 * t424;
t96 = -pkin(5) * t462 - t181 * pkin(11) + t186;
t37 = t423 * t96 + t428 * t63;
t36 = -t423 * t63 + t428 * t96;
t357 = t421 * t551 + t510;
t368 = pkin(4) * t413 + t409;
t373 = pkin(3) * t418 + pkin(4) * t412;
t416 = -pkin(10) - t607;
t468 = -t357 * t416 + t358 * t368 + t373 * t557 + t544;
t64 = t429 * t117 - t424 * t121;
t88 = t130 * t429 - t133 * t424;
t460 = t722 * pkin(4);
t171 = -pkin(4) * t200 + t288;
t169 = -t423 * t181 - t517;
t455 = -t428 * t181 + t518;
t450 = t142 * t475;
t449 = t143 * t477;
t448 = t160 * t474;
t447 = t293 * (mrSges(4,1) * t418 + t605);
t19 = t117 * t532 - t121 * t533 + t424 * t70 + t429 * t82;
t441 = t356 * t412 + t413 * t511;
t437 = t441 * pkin(4);
t436 = -qJD(6) * t472 - t3 * t423;
t435 = t2 * t428 + t436;
t8 = -pkin(5) * t326 - t11;
t433 = t16 * t624 + t2 * t600 + t428 * t671 + t8 * t480 + (Ifges(7,1) * t423 + t591) * t667 + (Ifges(7,2) * t428 + t592) * t666 + t528 + (Ifges(7,5) * t423 + Ifges(7,6) * t428) * t665 + t77 * t499 + (t454 + t514) * qJD(6) + (t450 + t449 + t448) * qJD(6) / 0.2e1 + t727;
t411 = -pkin(5) - t615;
t391 = Ifges(3,4) * t507;
t360 = t421 * t619 - t402;
t359 = (-mrSges(3,1) * t431 + mrSges(3,2) * t426) * t419;
t355 = -t421 * t509 + t552;
t340 = -t401 * mrSges(3,2) + t491;
t310 = t316 * pkin(5);
t302 = Ifges(3,1) * t508 + t391 + t575;
t282 = t358 * t413 + t412 * t557;
t278 = -mrSges(4,1) * t507 - t320 * mrSges(4,3);
t277 = mrSges(4,2) * t507 + t319 * mrSges(4,3);
t266 = t268 * pkin(5);
t230 = -mrSges(4,1) * t724 - mrSges(4,3) * t280;
t229 = mrSges(4,2) * t724 + mrSges(4,3) * t279;
t223 = Ifges(4,4) * t320 + Ifges(4,2) * t319 - Ifges(4,6) * t507;
t222 = -Ifges(4,3) * t507 + t586 + t590;
t218 = t273 * t428 + t357 * t423;
t217 = -t273 * t423 + t357 * t428;
t208 = mrSges(5,1) * t383 - t602;
t207 = -mrSges(5,2) * t383 + t603;
t179 = t280 * Ifges(4,1) + t279 * Ifges(4,4) - Ifges(4,5) * t724;
t178 = t280 * Ifges(4,4) + t279 * Ifges(4,2) - Ifges(4,6) * t724;
t167 = -mrSges(5,1) * t495 + mrSges(5,2) * t240;
t156 = t576 + t578 + t579;
t136 = -mrSges(5,2) * t337 + mrSges(5,3) * t148;
t135 = mrSges(5,1) * t337 - mrSges(5,3) * t147;
t112 = -mrSges(6,1) * t718 + mrSges(6,2) * t717;
t111 = qJD(5) * t175 + t210 * t424 - t429 * t432;
t107 = qJD(5) * t181 + t199 * t424 - t429 * t200;
t106 = qJD(5) * t462 + t199 * t429 + t200 * t424;
t102 = t577 + t580 + t581;
t94 = t113 + t618;
t85 = -pkin(5) * t508 - t88;
t84 = qJD(6) * t455 - t423 * t106 + t428 * t506;
t83 = qJD(6) * t169 + t428 * t106 + t423 * t506;
t62 = pkin(5) * t556 - t64;
t61 = -mrSges(6,2) * t326 + mrSges(6,3) * t81;
t60 = mrSges(6,1) * t326 - mrSges(6,3) * t80;
t59 = t108 * t429 - t568;
t58 = t108 * t424 + t550;
t42 = pkin(5) * t107 - pkin(11) * t106 + t171;
t33 = t113 * t423 + t428 * t56;
t32 = t113 * t428 - t423 * t56;
t30 = t423 * t94 + t428 * t59;
t29 = -t423 * t59 + t428 * t94;
t21 = -mrSges(7,1) * t48 + mrSges(7,2) * t47;
t18 = -pkin(5) * t506 - t20;
t17 = pkin(11) * t506 + t19;
t5 = -qJD(6) * t37 - t17 * t423 + t42 * t428;
t4 = qJD(6) * t36 + t17 * t428 + t42 * t423;
t1 = [-(-Ifges(3,6) * t421 / 0.2e1 + Ifges(4,5) * t630 + Ifges(4,6) * t631 - t361 * mrSges(3,3) + (-pkin(1) * mrSges(3,1) - t599 + (-Ifges(3,2) - Ifges(4,3)) * t431) * t419) * t724 + (Ifges(5,5) * t254 + Ifges(5,6) * t253) * t632 + ((t319 * t476 / 0.2e1 + t320 * t478 / 0.2e1 + t575 / 0.2e1 + t447 + t302 / 0.2e1 + t555 / 0.2e1 - t418 * t223 / 0.2e1 - t344 * mrSges(3,3) + (-t214 * t420 - t215 * t418) * mrSges(4,3)) * t431 + (-t677 + t586 / 0.2e1 - t587 / 0.2e1 + t222 / 0.2e1 + t580 / 0.2e1 + t581 / 0.2e1 + t590 / 0.2e1 + t156 / 0.2e1 + t578 / 0.2e1 + t579 / 0.2e1 + t576 / 0.2e1 + t577 / 0.2e1 + t214 * mrSges(4,1) - t215 * mrSges(4,2) - t345 * mrSges(3,3) + t102 / 0.2e1) * t426 + ((-Ifges(3,2) * t426 + t598) * t623 - t431 * (Ifges(4,5) * t554 - Ifges(4,6) * t559 + t583) / 0.2e1 - t685) * t543) * t540 + t155 * t350 * mrSges(4,3) - t154 * t351 * mrSges(4,3) - t731 * t462 + m(6) * (t10 * t703 + t11 * t64 + t114 * t186 + t168 * t171 + t19 * t57 + t20 * t56) + t703 * t61 + t8 * (-mrSges(7,1) * t169 - mrSges(7,2) * t455) + (-Ifges(7,1) * t455 + Ifges(7,4) * t169) * t667 + (-Ifges(7,5) * t455 + Ifges(7,6) * t169) * t665 + (-Ifges(7,4) * t455 + Ifges(7,2) * t169) * t666 + (t169 * t2 - t25 * t83 + t26 * t84 + t3 * t455) * mrSges(7,3) - t455 * t670 - (t528 + t513) * t556 / 0.2e1 + t728 * t106 + t729 * t107 + m(4) * (t154 * t243 + t155 * t244 + t214 * t226 + t215 * t227 + t241 * t333 + t293 * t347) + m(5) * (t123 * t91 + t124 * t90 + t137 * t52 + t138 * t51 + t182 * t259 + t236 * t288) + m(7) * (t18 * t54 + t2 * t37 + t25 * t5 + t26 * t4 + t3 * t36 + t62 * t8) + t714 * t181 + (-pkin(1) * t359 * t419 + Ifges(2,3)) * qJDD(1) + (Ifges(7,1) * t83 + Ifges(7,4) * t84) * t650 + t19 * t152 + t20 * t153 + t171 * t112 + (Ifges(3,3) * t625 + t360 * mrSges(3,1) - t361 * mrSges(3,2) + (Ifges(3,5) * t426 + Ifges(3,6) * t431) * t419) * t400 + (Ifges(3,5) * t625 - t360 * mrSges(3,3) + (-pkin(1) * mrSges(3,2) + t426 * Ifges(3,1) + t598) * t419) * t349 + (Ifges(7,5) * t83 + Ifges(7,6) * t84) * t644 + t258 * (mrSges(3,1) * t421 - mrSges(3,3) * t558) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t675 + t257 * t361 + t258 * t360 - t344 * t347 + t345 * t346) + (Ifges(5,4) * t254 + Ifges(5,2) * t253) * t648 + t137 * t135 + t138 * t136 + (Ifges(5,1) * t254 + Ifges(5,4) * t253) * t649 + (-t123 * t199 + t124 * t200 + t253 * t51 - t254 * t52) * mrSges(5,3) + (-m(4) * (-pkin(2) * t356 + t498) - (-t356 * t420 + t489) * mrSges(4,1) - (t356 * t418 + t420 * t511) * mrSges(4,2) - m(5) * (pkin(3) * t489 - t356 * t409 + t498) - t494 * mrSges(5,1) - t441 * mrSges(5,2) - m(3) * t498 + t356 * mrSges(3,1) - mrSges(3,3) * t511 + t427 * mrSges(2,1) + t622 * mrSges(2,2) + t493 * t268 - t739 * t269 + t694 * t355 + t710 * (-t355 * t416 + t356 * t368 - t373 * t511 - t498)) * g(1) + t4 * t100 + t5 * t101 + t18 * t95 + t54 * (-mrSges(7,1) * t84 + mrSges(7,2) * t83) + t84 * t77 / 0.2e1 + t83 * t78 / 0.2e1 + t62 * t21 + t64 * t60 + t36 * t23 + t37 * t24 + t186 * t31 + t697 * t347 + (Ifges(7,4) * t83 + Ifges(7,2) * t84) * t652 + t90 * t207 + t91 * t208 + (-t622 * mrSges(2,1) - t282 * mrSges(5,1) - t281 * mrSges(5,2) - m(7) * (pkin(5) * t273 + t468) - t218 * mrSges(7,1) - t217 * mrSges(7,2) - m(6) * t468 - t273 * mrSges(6,1) + (-mrSges(3,1) + t723) * t358 + (mrSges(2,2) + (-t605 - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t418) * t419) * t427 + t493 * t272 + t438 * t357 + (-m(3) + t711) * t544) * g(2) + t236 * (-mrSges(5,1) * t200 + mrSges(5,2) * t199) + t243 * t230 + t244 * t229 + t512 * t625 + (Ifges(5,5) * t199 + Ifges(5,6) * t200) * t626 + t179 * t630 + t178 * t631 + (Ifges(4,1) * t351 + Ifges(4,4) * t350) * t634 + (Ifges(4,4) * t351 + Ifges(4,2) * t350) * t635 + (Ifges(5,1) * t199 + Ifges(5,4) * t200) * t636 + (Ifges(5,4) * t199 + Ifges(5,2) * t200) * t638 + t199 * t646 + t200 * t647 + t254 * t658 + t253 * t659 + t169 * t671 + t182 * (-mrSges(5,1) * t253 + mrSges(5,2) * t254) + t259 * t97 + t227 * t277 + t226 * t278 - t421 * t726 + t288 * t167 + (-mrSges(4,1) * t154 + mrSges(4,2) * t155 + mrSges(3,3) * t257 - Ifges(4,5) * t280 - Ifges(5,5) * t649 - Ifges(6,5) * t661 - Ifges(4,6) * t279 - Ifges(5,6) * t648 - Ifges(6,6) * t660 - Ifges(5,3) * t632 - Ifges(6,3) * t633 - t720 - t727) * t556 + t333 * t198 + t346 * t340 + t241 * (-mrSges(4,1) * t350 + mrSges(4,2) * t351); (-Ifges(7,4) * t452 - Ifges(7,2) * t453) * t652 + (Ifges(7,4) * t192 + Ifges(7,2) * t191) * t653 + ((t583 + (Ifges(4,5) * t420 - Ifges(4,6) * t418) * t431) * t623 + t685) * qJD(1) ^ 2 * t675 - t724 * (Ifges(4,5) * t418 + Ifges(4,6) * t420) / 0.2e1 + (-Ifges(7,1) * t452 - Ifges(7,4) * t453) * t650 + (Ifges(7,1) * t192 + Ifges(7,4) * t191) * t651 - t731 * t461 + t733 * t647 - (t21 - t60) * t463 + (t10 * t175 + t11 * t463 + t114 * t324 + t706 * t57 + (-t111 - t88) * t56 + t686 * t168) * m(6) + (t119 * t3 + t120 * t2 - t463 * t8 + (t111 - t85) * t54 + t708 * t26 + t709 * t25) * m(7) + (Ifges(5,5) * t299 - Ifges(5,6) * t298) * t627 + (Ifges(5,4) * t299 - Ifges(5,2) * t298) * t639 + (Ifges(5,1) * t299 - Ifges(5,4) * t298) * t637 + (t514 + t728) * t193 + t729 * t194 + (Ifges(5,5) * t637 + Ifges(6,5) * t641 - Ifges(3,2) * t484 + Ifges(5,6) * t639 + Ifges(6,6) * t643 + Ifges(5,3) * t627 + Ifges(6,3) * t629 + t677) * t508 + (Ifges(7,5) * t192 + Ifges(7,6) * t191) * t645 + (-t538 - t248) * t278 + (-Ifges(7,5) * t452 - Ifges(7,6) * t453) * t644 + (-t710 * (-t355 * t368 - t356 * t416) - t694 * t356 + t681 * t355) * g(2) + (-t710 * (-t357 * t368 - t358 * t416) - t694 * t358 + t681 * t357) * g(1) + (t474 * t665 + t475 * t666 + t477 * t667 + t479 * t8 + t499 * t78 + t714) * t276 + (t391 + t555 + t302) * t484 + t223 * t488 / 0.2e1 - t88 * t153 + (t537 - t249) * t277 + (Ifges(5,5) * t353 - Ifges(5,6) * t354) * t626 + (Ifges(5,1) * t353 - Ifges(5,4) * t354) * t636 + (Ifges(5,4) * t353 - Ifges(5,2) * t354) * t638 + ((t222 + t156 + t102) * t426 + t320 * (Ifges(4,5) * t426 + t431 * t478) + t319 * (Ifges(4,6) * t426 + t431 * t476) + t401 * (Ifges(3,5) * t431 - Ifges(3,6) * t426)) * t500 + (t491 - t340) * t344 + t241 * t482 - t15 * t564 / 0.2e1 - t447 * t507 + (t229 * t420 - t230 * t418) * qJ(3) + (-mrSges(5,1) * t733 - mrSges(5,2) * t695) * t236 + (t123 * t695 + t124 * t733 + t366 * t51 - t367 * t52) * mrSges(5,3) + t119 * t23 + t120 * t24 - t85 * t95 + t420 * t178 / 0.2e1 + t418 * t179 / 0.2e1 - t726 + t175 * t61 + t686 * t112 + t688 * mrSges(4,3) + (-t214 * t248 - t215 * t249 - pkin(2) * t241 + (-t214 * t418 + t215 * t420) * qJD(3) + t688 * qJ(3)) * m(4) + (-m(4) * t293 + t492 - t697) * t345 - t699 * t77 / 0.2e1 + (mrSges(7,1) * t699 - mrSges(7,2) * t698) * t54 + (-t2 * t564 + t25 * t698 - t26 * t699 - t3 * t563) * mrSges(7,3) + t701 * t207 - t192 * t78 / 0.2e1 + t702 * t208 + (t123 * t702 + t124 * t701 - t182 * t409 - t236 * t287 + t296 * t52 + t297 * t51) * m(5) - pkin(2) * t198 + t721 * t211 - t571 * t111 + t706 * t152 + (Ifges(5,5) * t367 + Ifges(5,6) * t366) * t632 + (Ifges(4,1) * t418 + t596) * t634 + (Ifges(4,2) * t420 + t597) * t635 + t353 * t646 + (Ifges(5,4) * t367 + Ifges(5,2) * t366) * t648 + (Ifges(5,1) * t367 + Ifges(5,4) * t366) * t649 + t367 * t658 + t366 * t659 + t563 * t670 + t258 * mrSges(3,1) + (-t710 * t368 * t556 + t359 + (t719 * t431 + (t416 * t710 - t479 + t725) * t426) * t419) * g(3) + t708 * t100 + t709 * t101 + (t678 + t610) * t212 - t287 * t167 + t296 * t135 + t297 * t136 - t299 * t158 / 0.2e1 + t324 * t31 + (-t215 * (-mrSges(4,2) * t426 - mrSges(4,3) * t559) - t214 * (mrSges(4,1) * t426 - mrSges(4,3) * t554)) * t543 + t512 + t182 * (-mrSges(5,1) * t366 + mrSges(5,2) * t367) - t409 * t97; t700 * t718 + t467 * qJD(6) + t31 + t97 + t571 * t717 + t198 + t428 * t23 + t423 * t24 + t240 * t208 - t495 * t207 - t319 * t277 + t320 * t278 + (t160 * t471 + t2 * t423 + t3 * t428 - t717 * t54) * m(7) + (t56 * t717 - t57 * t718 + t114) * m(6) + (t123 * t240 - t124 * t495 + t182) * m(5) + (t214 * t320 - t215 * t319 + t241) * m(4) + (-g(1) * t357 - g(2) * t355 + g(3) * t556) * (t710 - t711); -t700 * pkin(4) * t532 + (t411 * t8 + (t424 * t54 + t429 * t471) * qJD(5) * pkin(4) - t25 * t29 - t26 * t30 - t54 * t58) * m(7) + (-t168 * t618 + t56 * t58 - t57 * t59 + (t10 * t424 + t11 * t429 + (-t424 * t56 + t429 * t57) * qJD(5)) * pkin(4)) * m(6) + t721 * t717 + (t25 * t600 + t26 * t601 + t474 * t645 + t475 * t653 + t477 * t651 + t678 + t716) * t718 - t236 * (mrSges(5,1) * t240 + mrSges(5,2) * t495) + (Ifges(5,5) * t495 - Ifges(5,6) * t240) * t627 + (Ifges(5,1) * t495 - t595) * t637 - t3 * t601 - t112 * t618 - t59 * t152 + (t208 + t602) * t124 + t720 + (-Ifges(5,2) * t240 + t158 + t235) * t639 + (-t207 + t603) * t123 - t30 * t100 - t29 * t101 + (-t25 * t530 - t26 * t531) * mrSges(7,3) + (m(7) * t435 - t100 * t531 - t101 * t530 + t687) * (pkin(11) + t616) + (m(6) * t437 - m(7) * (t269 * pkin(11) + t266 - t437) + mrSges(5,1) * t441 - mrSges(5,2) * t494 + t690) * g(2) + (-m(6) * t483 - m(7) * (t483 + t497) - mrSges(5,1) * t281 + mrSges(5,2) * t282 + t691) * g(1) + (-t722 * mrSges(5,1) - (-t412 * t421 - t413 * t558) * mrSges(5,2) - m(7) * (pkin(11) * t317 + t310 + t460) - m(6) * t460 + t692) * g(3) + t60 * t615 + t61 * t616 + t157 * t636 + t571 * (-pkin(4) * t533 + t58) + t433 + t513 + t411 * t21; -t56 * t152 - m(7) * (t25 * t32 + t26 * t33 + t54 * t57) + (-m(7) * t266 + t690) * g(2) + (-m(7) * t497 + t691) * g(1) + (-m(7) * t310 + t692) * g(3) + (t732 + t738) * t717 - t33 * t100 - t32 * t101 - pkin(5) * t21 + (-t450 / 0.2e1 - t449 / 0.2e1 - t448 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t717 + t472 * mrSges(7,3) + t684 + t716 + t740) * t718 + t436 * mrSges(7,3) + ((-t100 * t423 - t101 * t428) * qJD(6) + (-g(2) * t269 - g(3) * t317) * m(7) + t687) * pkin(11) + t571 * t57 - t8 * t674 + t435 * t673 + t433; -t54 * (mrSges(7,1) * t143 + mrSges(7,2) * t142) + (Ifges(7,1) * t142 - t593) * t651 + t77 * t650 + (Ifges(7,5) * t142 - Ifges(7,6) * t143) * t645 - t25 * t100 + t26 * t101 - g(1) * (mrSges(7,1) * t217 - mrSges(7,2) * t218) - g(2) * ((-t269 * t423 + t355 * t428) * mrSges(7,1) + (-t269 * t428 - t355 * t423) * mrSges(7,2)) - g(3) * ((-t317 * t423 - t517) * mrSges(7,1) + (-t317 * t428 + t518) * mrSges(7,2)) + (t142 * t25 + t143 * t26) * mrSges(7,3) + t14 + (-Ifges(7,2) * t143 + t141 + t78) * t653 + t683;];
tau  = t1;
