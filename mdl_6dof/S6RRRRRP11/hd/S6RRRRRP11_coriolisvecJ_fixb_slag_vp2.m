% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:41
% EndTime: 2019-03-10 02:34:17
% DurationCPUTime: 48.71s
% Computational Cost: add. (37561->988), mult. (111923->1353), div. (0->0), fcn. (92634->12), ass. (0->444)
t409 = cos(qJ(2));
t401 = cos(pkin(6));
t562 = pkin(1) * t401;
t394 = t409 * t562;
t384 = qJD(1) * t394;
t405 = sin(qJ(2));
t399 = sin(pkin(6));
t400 = cos(pkin(7));
t474 = pkin(10) * t400 + pkin(9);
t452 = t399 * t474;
t427 = t405 * t452;
t303 = -qJD(1) * t427 + t384;
t393 = t405 * t562;
t416 = -t409 * t452 - t393;
t304 = t416 * qJD(1);
t398 = sin(pkin(7));
t560 = pkin(10) * t398;
t418 = t399 * (pkin(2) * t405 - t409 * t560);
t337 = qJD(1) * t418;
t404 = sin(qJ(3));
t511 = t398 * t404;
t388 = pkin(10) * t511;
t408 = cos(qJ(3));
t506 = t400 * t408;
t363 = pkin(2) * t506 - t388;
t507 = t400 * t404;
t624 = t363 * qJD(3) - t408 * t303 - t304 * t507 - t337 * t511;
t249 = -t304 * t398 + t400 * t337;
t500 = t405 * t408;
t501 = t404 * t409;
t424 = t400 * t500 + t501;
t489 = qJD(1) * t399;
t322 = t424 * t489;
t498 = t408 * t409;
t502 = t404 * t405;
t422 = -t400 * t502 + t498;
t323 = t422 * t489;
t694 = -pkin(3) * t322 + pkin(11) * t323 - t249 + (pkin(3) * t404 - pkin(11) * t408) * t398 * qJD(3);
t473 = t405 * t489;
t455 = t398 * t473;
t693 = pkin(11) * t455 - t624;
t387 = qJD(1) * t401 + qJD(2);
t425 = t400 * t498 - t502;
t510 = t398 * t408;
t274 = t387 * t510 + t425 * t489;
t272 = qJD(4) - t274;
t403 = sin(qJ(4));
t407 = cos(qJ(4));
t267 = t323 * t403 - t407 * t455;
t360 = t400 * t403 + t407 * t511;
t486 = qJD(3) * t408;
t469 = t398 * t486;
t308 = qJD(4) * t360 + t403 * t469;
t494 = t267 - t308;
t268 = t323 * t407 + t403 * t455;
t359 = -t407 * t400 + t403 * t511;
t307 = -qJD(4) * t359 + t407 * t469;
t493 = t268 - t307;
t365 = pkin(2) * t507 + pkin(10) * t510;
t692 = -t365 * qJD(3) + t404 * t303;
t487 = qJD(3) * t404;
t470 = t398 * t487;
t691 = t322 - t470;
t645 = Ifges(6,4) + Ifges(7,4);
t423 = t400 * t501 + t500;
t417 = t423 * t399;
t275 = qJD(1) * t417 + t387 * t511;
t472 = t409 * t489;
t327 = t400 * t387 - t398 * t472 + qJD(3);
t238 = -t275 * t403 + t327 * t407;
t414 = (qJD(2) * t422 + qJD(3) * t425) * t399;
t240 = qJD(1) * t414 + t387 * t469;
t488 = qJD(2) * t399;
t460 = qJD(1) * t488;
t453 = t405 * t460;
t428 = t398 * t453;
t154 = qJD(4) * t238 + t240 * t407 + t403 * t428;
t601 = t154 / 0.2e1;
t690 = Ifges(5,4) * t601;
t345 = pkin(11) * t400 + t365;
t346 = (-pkin(3) * t408 - pkin(11) * t404 - pkin(2)) * t398;
t484 = qJD(4) * t407;
t485 = qJD(4) * t403;
t627 = -t345 * t485 + t346 * t484 + t403 * t694 - t693 * t407;
t508 = t399 * t409;
t269 = t387 * t560 + (t474 * t508 + t393) * qJD(1);
t273 = pkin(2) * t387 + t303;
t333 = (-pkin(2) * t409 - t405 * t560 - pkin(1)) * t399;
t318 = qJD(1) * t333;
t430 = t273 * t400 + t318 * t398;
t184 = t269 * t408 + t404 * t430;
t689 = -t184 + t272 * (pkin(4) * t403 - pkin(12) * t407);
t626 = t304 * t506 - (-pkin(3) * t473 - t337 * t408) * t398 - t692;
t239 = t275 * t407 + t327 * t403;
t155 = qJD(4) * t239 + t240 * t403 - t407 * t428;
t599 = t155 / 0.2e1;
t402 = sin(qJ(5));
t406 = cos(qJ(5));
t192 = t239 * t406 + t272 * t402;
t413 = (qJD(2) * t424 + qJD(3) * t423) * t399;
t241 = qJD(1) * t413 + t387 * t470;
t74 = -qJD(5) * t192 - t154 * t402 + t241 * t406;
t609 = t74 / 0.2e1;
t191 = -t239 * t402 + t272 * t406;
t73 = qJD(5) * t191 + t154 * t406 + t241 * t402;
t610 = t73 / 0.2e1;
t646 = Ifges(6,1) + Ifges(7,1);
t672 = -Ifges(7,5) - Ifges(6,5);
t673 = -t599 * t672 + t645 * t609 + t646 * t610;
t600 = -t155 / 0.2e1;
t579 = t241 / 0.2e1;
t643 = Ifges(6,2) + Ifges(7,2);
t642 = Ifges(6,6) + Ifges(7,6);
t671 = -Ifges(7,3) - Ifges(6,3);
t688 = pkin(12) * t691 - t627;
t687 = -pkin(4) * t494 + t493 * pkin(12) + t626;
t306 = t416 * qJD(2);
t289 = qJD(1) * t306;
t338 = qJD(2) * t418;
t328 = qJD(1) * t338;
t379 = qJD(2) * t384;
t420 = qJD(2) * t427;
t288 = -qJD(1) * t420 + t379;
t468 = t400 * t487;
t457 = -t269 * t486 - t273 * t468 - t404 * t288 - t318 * t470;
t110 = -t289 * t506 + (-pkin(3) * t453 - t328 * t408) * t398 - t457;
t467 = t400 * t486;
t116 = -t269 * t487 + t273 * t467 + t408 * t288 + t289 * t507 + t318 * t469 + t328 * t511;
t109 = pkin(11) * t428 + t116;
t244 = -t289 * t398 + t400 * t328;
t133 = pkin(3) * t241 - pkin(11) * t240 + t244;
t231 = -t273 * t398 + t400 * t318;
t162 = -pkin(3) * t274 - pkin(11) * t275 + t231;
t165 = pkin(11) * t327 + t184;
t21 = t407 * t109 + t403 * t133 + t162 * t484 - t165 * t485;
t18 = pkin(12) * t241 + t21;
t82 = t162 * t403 + t165 * t407;
t76 = pkin(12) * t272 + t82;
t183 = -t404 * t269 + t408 * t430;
t164 = -pkin(3) * t327 - t183;
t95 = -pkin(4) * t238 - pkin(12) * t239 + t164;
t34 = t402 * t95 + t406 * t76;
t41 = pkin(4) * t155 - pkin(12) * t154 + t110;
t6 = -qJD(5) * t34 - t18 * t402 + t406 * t41;
t1 = pkin(5) * t155 - qJ(6) * t73 - qJD(6) * t192 + t6;
t482 = qJD(5) * t406;
t483 = qJD(5) * t402;
t5 = t406 * t18 + t402 * t41 + t95 * t482 - t483 * t76;
t2 = qJ(6) * t74 + qJD(6) * t191 + t5;
t656 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t686 = -t110 * mrSges(5,1) + 0.2e1 * Ifges(5,2) * t600 + 0.2e1 * Ifges(5,6) * t579 + t671 * t599 - t642 * t609 + t672 * t610 - t656 + t690;
t234 = qJD(5) - t238;
t679 = t645 * t192;
t634 = t191 * t643 + t234 * t642 + t679;
t685 = -t634 / 0.2e1;
t224 = pkin(3) * t275 - pkin(11) * t274;
t130 = t407 * t183 + t403 * t224;
t108 = pkin(12) * t275 + t130;
t559 = pkin(11) * t402;
t684 = t108 * t402 + t406 * t689 + t485 * t559;
t378 = -pkin(4) * t407 - pkin(12) * t403 - pkin(3);
t683 = -t406 * t108 + t378 * t482 + t402 * t689;
t628 = -t345 * t484 - t346 * t485 + t693 * t403 + t407 * t694;
t682 = t645 * t191;
t309 = -t360 * t402 - t406 * t510;
t227 = qJD(5) * t309 + t307 * t406 + t402 * t470;
t230 = t268 * t406 + t322 * t402;
t496 = t227 - t230;
t426 = -t360 * t406 + t402 * t510;
t228 = qJD(5) * t426 - t307 * t402 + t406 * t470;
t229 = -t268 * t402 + t322 * t406;
t495 = t228 - t229;
t681 = t645 * t406;
t680 = t645 * t402;
t233 = Ifges(5,4) * t238;
t525 = t272 * Ifges(5,5);
t33 = -t402 * t76 + t406 * t95;
t23 = -qJ(6) * t192 + t33;
t20 = pkin(5) * t234 + t23;
t24 = qJ(6) * t191 + t34;
t431 = t33 * t406 + t34 * t402;
t445 = mrSges(7,1) * t402 + mrSges(7,2) * t406;
t447 = mrSges(6,1) * t402 + mrSges(6,2) * t406;
t81 = t162 * t407 - t403 * t165;
t75 = -pkin(4) * t272 - t81;
t54 = -pkin(5) * t191 + qJD(6) + t75;
t564 = t406 / 0.2e1;
t567 = -t402 / 0.2e1;
t584 = t234 / 0.2e1;
t592 = t192 / 0.2e1;
t594 = t191 / 0.2e1;
t633 = t192 * t646 - t234 * t672 + t682;
t661 = t406 * t646 - t680;
t662 = -t402 * t643 + t681;
t663 = -t402 * t642 - t406 * t672;
t615 = t431 * mrSges(6,3) + (t20 * t406 + t24 * t402) * mrSges(7,3) - t445 * t54 - t447 * t75 - t662 * t594 - t661 * t592 - t663 * t584 - t634 * t567 - t633 * t564;
t678 = t615 - t164 * mrSges(5,2) - t233 / 0.2e1 - t525 / 0.2e1 + t81 * mrSges(5,3);
t676 = t240 / 0.2e1;
t675 = -t241 / 0.2e1;
t674 = t327 / 0.2e1;
t640 = -t155 * t671 + t642 * t74 - t672 * t73;
t639 = t155 * t642 + t643 * t74 + t645 * t73;
t83 = t192 * Ifges(7,5) + t191 * Ifges(7,6) + t234 * Ifges(7,3);
t84 = t192 * Ifges(6,5) + t191 * Ifges(6,6) + t234 * Ifges(6,3);
t635 = t84 + t83;
t344 = t388 + (-pkin(2) * t408 - pkin(3)) * t400;
t256 = pkin(4) * t359 - pkin(12) * t360 + t344;
t262 = t407 * t345 + t403 * t346;
t258 = -pkin(12) * t510 + t262;
t632 = t256 * t482 - t258 * t483 + t402 * t687 - t406 * t688;
t190 = t402 * t256 + t406 * t258;
t631 = -qJD(5) * t190 + t402 * t688 + t406 * t687;
t670 = Ifges(5,4) * t600;
t669 = t231 * mrSges(4,2);
t499 = t406 * t407;
t216 = t274 * t499 + t275 * t402;
t395 = pkin(11) * t499;
t481 = qJD(6) * t406;
t512 = t274 * t403;
t668 = -pkin(5) * t512 + qJ(6) * t216 - t403 * t481 + (pkin(5) * t403 - qJ(6) * t499) * qJD(4) + (-t395 + (qJ(6) * t403 - t378) * t402) * qJD(5) + t684;
t504 = t402 * t407;
t215 = -t274 * t504 + t275 * t406;
t503 = t403 * t406;
t667 = -qJ(6) * t215 + (-pkin(11) * qJD(4) - qJ(6) * qJD(5)) * t503 + (-qJD(6) * t403 + (-pkin(11) * qJD(5) - qJ(6) * qJD(4)) * t407) * t402 + t683;
t340 = t402 * t378 + t395;
t666 = -qJD(5) * t340 + t684;
t665 = (-t406 * t485 - t407 * t483) * pkin(11) + t683;
t129 = -t403 * t183 + t224 * t407;
t107 = -pkin(4) * t275 - t129;
t664 = pkin(11) * t484 - t107 + (t402 * t484 + t403 * t482 + t215) * pkin(5);
t629 = pkin(4) * t691 - t628;
t660 = -t231 * mrSges(4,1) - t81 * mrSges(5,1) + t82 * mrSges(5,2);
t570 = -t327 / 0.2e1;
t574 = -t275 / 0.2e1;
t659 = Ifges(4,1) * t574 + Ifges(4,5) * t570;
t658 = Ifges(5,1) * t601 + Ifges(5,5) * t579;
t611 = t658 + t670;
t657 = t110 * mrSges(5,2) + t611 + t658;
t655 = -t33 * mrSges(6,1) - t20 * mrSges(7,1) + t34 * mrSges(6,2) + t24 * mrSges(7,2);
t575 = -t274 / 0.2e1;
t577 = -t272 / 0.2e1;
t581 = -t239 / 0.2e1;
t583 = -t238 / 0.2e1;
t654 = Ifges(5,5) * t581 - Ifges(4,2) * t575 - Ifges(4,6) * t570 + Ifges(5,6) * t583 + Ifges(5,3) * t577;
t653 = t164 * mrSges(5,1) - t82 * mrSges(5,3) - t655;
t22 = -t403 * t109 + t133 * t407 - t162 * t485 - t165 * t484;
t19 = -pkin(4) * t241 - t22;
t9 = -pkin(5) * t74 + t19;
t651 = -t19 * mrSges(6,1) - t9 * mrSges(7,1) + t639 / 0.2e1 + t642 * t599 + t643 * t609 + t645 * t610 + t5 * mrSges(6,3) + t2 * mrSges(7,3);
t650 = t19 * mrSges(6,2) + t9 * mrSges(7,2) - t6 * mrSges(6,3) - t1 * mrSges(7,3) + 0.2e1 * t673;
t649 = -t690 - t21 * mrSges(5,3) + t640 / 0.2e1 - t686;
t648 = Ifges(4,3) * t674;
t647 = -Ifges(3,6) * t387 / 0.2e1;
t637 = -pkin(5) * t494 - qJ(6) * t496 + qJD(6) * t426 + t631;
t636 = qJ(6) * t495 + qJD(6) * t309 + t632;
t630 = -pkin(5) * t495 + t629;
t490 = pkin(9) * t508 + t393;
t290 = (t398 * t401 + t400 * t508) * pkin(10) + t490;
t302 = pkin(2) * t401 + t394 - t427;
t205 = -t404 * t290 + t408 * (t302 * t400 + t333 * t398);
t358 = -t398 * t508 + t401 * t400;
t185 = -pkin(3) * t358 - t205;
t297 = t401 * t511 + t417;
t254 = t297 * t403 - t358 * t407;
t255 = t297 * t407 + t358 * t403;
t118 = pkin(4) * t254 - pkin(12) * t255 + t185;
t245 = -t302 * t398 + t400 * t333;
t296 = -t399 * t425 - t401 * t510;
t179 = pkin(3) * t296 - pkin(11) * t297 + t245;
t206 = t408 * t290 + t302 * t507 + t333 * t511;
t186 = pkin(11) * t358 + t206;
t101 = t403 * t179 + t407 * t186;
t98 = pkin(12) * t296 + t101;
t46 = t402 * t118 + t406 * t98;
t625 = -(t304 * t400 + t337 * t398) * t408 + t692;
t623 = -t402 * t672 + t406 * t642;
t622 = t406 * t643 + t680;
t621 = t402 * t646 + t681;
t620 = t21 * t407 - t22 * t403;
t619 = -t402 * t6 + t406 * t5;
t618 = Ifges(4,1) * t676 + Ifges(4,4) * t675;
t617 = t183 * mrSges(4,3) - t669;
t614 = t184 * mrSges(4,3) + t660;
t524 = t272 * Ifges(5,6);
t548 = Ifges(5,4) * t239;
t532 = t238 * Ifges(5,2);
t144 = t524 + t532 + t548;
t603 = t144 / 0.2e1;
t613 = (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t191 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t192 - (-Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1) * t234 - t603 + t83 / 0.2e1 + t84 / 0.2e1 - t548 / 0.2e1 - t524 / 0.2e1 + t653;
t606 = pkin(1) * mrSges(3,1);
t605 = pkin(1) * mrSges(3,2);
t523 = t272 * Ifges(5,3);
t529 = t239 * Ifges(5,5);
t531 = t238 * Ifges(5,6);
t143 = t523 + t529 + t531;
t604 = t143 / 0.2e1;
t530 = t239 * Ifges(5,1);
t145 = t233 + t525 + t530;
t602 = t145 / 0.2e1;
t236 = Ifges(4,6) * t241;
t237 = Ifges(4,5) * t240;
t158 = Ifges(4,3) * t428 - t236 + t237;
t598 = t158 / 0.2e1;
t597 = Ifges(4,5) * t428 / 0.2e1 + t618;
t595 = -t191 / 0.2e1;
t593 = -t192 / 0.2e1;
t516 = t327 * Ifges(4,6);
t519 = t275 * Ifges(4,4);
t522 = t274 * Ifges(4,2);
t202 = t516 + t519 + t522;
t591 = -t202 / 0.2e1;
t271 = Ifges(4,4) * t274;
t517 = t327 * Ifges(4,5);
t520 = t275 * Ifges(4,1);
t203 = t271 + t517 + t520;
t590 = t203 / 0.2e1;
t585 = -t234 / 0.2e1;
t582 = t238 / 0.2e1;
t580 = t239 / 0.2e1;
t576 = t272 / 0.2e1;
t573 = t275 / 0.2e1;
t568 = t401 / 0.2e1;
t561 = pkin(5) * t402;
t553 = qJD(2) / 0.2e1;
t552 = -qJ(6) - pkin(12);
t551 = mrSges(4,3) * t274;
t550 = mrSges(4,3) * t275;
t549 = Ifges(3,4) * t405;
t547 = Ifges(5,4) * t403;
t546 = Ifges(5,4) * t407;
t153 = Ifges(5,5) * t154;
t152 = Ifges(5,6) * t155;
t538 = t116 * mrSges(4,2);
t117 = (t289 * t400 + t328 * t398) * t408 + t457;
t537 = t117 * mrSges(4,1);
t535 = t19 * t403;
t527 = t240 * Ifges(4,4);
t170 = pkin(4) * t239 - pkin(12) * t238;
t58 = t402 * t170 + t406 * t81;
t513 = t238 * t402;
t509 = t399 * t405;
t505 = t402 * t403;
t115 = -mrSges(6,1) * t191 + mrSges(6,2) * t192;
t195 = mrSges(5,1) * t272 - mrSges(5,3) * t239;
t497 = -t115 + t195;
t63 = Ifges(5,3) * t241 - t152 + t153;
t471 = t405 * t488;
t28 = -t74 * mrSges(7,1) + t73 * mrSges(7,2);
t459 = qJD(5) * t552;
t45 = t406 * t118 - t402 * t98;
t57 = t406 * t170 - t402 * t81;
t100 = t179 * t407 - t403 * t186;
t189 = t406 * t256 - t258 * t402;
t250 = -t306 * t398 + t400 * t338;
t261 = -t403 * t345 + t346 * t407;
t351 = t490 * qJD(1);
t458 = t647 - (Ifges(3,2) * t409 + t549) * t489 / 0.2e1 - t351 * mrSges(3,3);
t385 = qJD(2) * t394;
t305 = t385 - t420;
t456 = -t290 * t486 - t302 * t468 - t404 * t305 - t333 * t470;
t454 = t398 * t471;
t451 = -t22 * mrSges(5,1) + t21 * mrSges(5,2);
t257 = pkin(4) * t510 - t261;
t449 = t537 - t538;
t448 = mrSges(6,1) * t406 - mrSges(6,2) * t402;
t446 = mrSges(7,1) * t406 - mrSges(7,2) * t402;
t212 = t255 * t406 + t296 * t402;
t211 = -t255 * t402 + t296 * t406;
t134 = -t290 * t487 + t302 * t467 + t408 * t305 + t306 * t507 + t333 * t469 + t338 * t511;
t125 = pkin(11) * t454 + t134;
t247 = t401 * t470 + t413;
t248 = t401 * t469 + t414;
t147 = pkin(3) * t247 - pkin(11) * t248 + t250;
t32 = -t403 * t125 + t147 * t407 - t179 * t485 - t186 * t484;
t97 = -pkin(4) * t296 - t100;
t31 = t407 * t125 + t403 * t147 + t179 * t484 - t186 * t485;
t26 = pkin(12) * t247 + t31;
t126 = -t306 * t506 + (-pkin(3) * t471 - t338 * t408) * t398 - t456;
t172 = qJD(4) * t255 + t248 * t403 - t407 * t454;
t173 = -qJD(4) * t254 + t248 * t407 + t403 * t454;
t53 = pkin(4) * t172 - pkin(12) * t173 + t126;
t7 = t118 * t482 + t406 * t26 + t402 * t53 - t483 * t98;
t349 = -pkin(9) * t473 + t384;
t383 = Ifges(3,4) * t472;
t419 = Ifges(3,1) * t473 / 0.2e1 + t383 / 0.2e1 + t387 * Ifges(3,5) - t349 * mrSges(3,3);
t27 = -pkin(4) * t247 - t32;
t355 = t490 * qJD(2);
t8 = -qJD(5) * t46 - t26 * t402 + t406 * t53;
t415 = t183 * mrSges(4,1) - t184 * mrSges(4,2) + t275 * Ifges(4,5) + t274 * Ifges(4,6) + t648;
t397 = -pkin(5) * t406 - pkin(4);
t381 = t552 * t406;
t380 = t552 * t402;
t377 = Ifges(3,5) * t409 * t460;
t376 = (pkin(11) + t561) * t403;
t372 = t406 * t378;
t364 = -pkin(9) * t509 + t394;
t357 = -qJD(6) * t402 + t406 * t459;
t356 = t402 * t459 + t481;
t354 = -pkin(9) * t471 + t385;
t348 = -t387 * mrSges(3,2) + mrSges(3,3) * t472;
t347 = mrSges(3,1) * t387 - mrSges(3,3) * t473;
t342 = qJD(1) * t355;
t341 = -pkin(9) * t453 + t379;
t339 = -pkin(11) * t504 + t372;
t312 = -qJ(6) * t505 + t340;
t295 = -qJ(6) * t503 + t372 + (-pkin(5) - t559) * t407;
t243 = mrSges(4,1) * t327 - t550;
t242 = -mrSges(4,2) * t327 + t551;
t223 = -mrSges(4,1) * t274 + mrSges(4,2) * t275;
t219 = -mrSges(4,2) * t428 - mrSges(4,3) * t241;
t218 = mrSges(4,1) * t428 - mrSges(4,3) * t240;
t217 = -pkin(5) * t309 + t257;
t194 = -mrSges(5,2) * t272 + mrSges(5,3) * t238;
t171 = mrSges(4,1) * t241 + mrSges(4,2) * t240;
t169 = -mrSges(5,1) * t238 + mrSges(5,2) * t239;
t166 = qJ(6) * t309 + t190;
t159 = -t241 * Ifges(4,2) + Ifges(4,6) * t428 + t527;
t149 = pkin(5) * t359 + qJ(6) * t426 + t189;
t142 = mrSges(6,1) * t234 - mrSges(6,3) * t192;
t141 = mrSges(7,1) * t234 - mrSges(7,3) * t192;
t140 = -mrSges(6,2) * t234 + mrSges(6,3) * t191;
t139 = -mrSges(7,2) * t234 + mrSges(7,3) * t191;
t135 = (t306 * t400 + t338 * t398) * t408 + t456;
t120 = -mrSges(5,2) * t241 - mrSges(5,3) * t155;
t119 = mrSges(5,1) * t241 - mrSges(5,3) * t154;
t114 = -mrSges(7,1) * t191 + mrSges(7,2) * t192;
t92 = qJD(5) * t211 + t173 * t406 + t247 * t402;
t91 = -qJD(5) * t212 - t173 * t402 + t247 * t406;
t79 = mrSges(5,1) * t155 + mrSges(5,2) * t154;
t67 = -pkin(5) * t211 + t97;
t66 = pkin(5) * t513 + t82;
t50 = -mrSges(6,2) * t155 + mrSges(6,3) * t74;
t49 = -mrSges(7,2) * t155 + mrSges(7,3) * t74;
t48 = mrSges(6,1) * t155 - mrSges(6,3) * t73;
t47 = mrSges(7,1) * t155 - mrSges(7,3) * t73;
t43 = -qJ(6) * t513 + t58;
t36 = -qJ(6) * t238 * t406 + pkin(5) * t239 + t57;
t35 = qJ(6) * t211 + t46;
t30 = pkin(5) * t254 - qJ(6) * t212 + t45;
t29 = -mrSges(6,1) * t74 + mrSges(6,2) * t73;
t10 = -pkin(5) * t91 + t27;
t4 = qJ(6) * t91 + qJD(6) * t211 + t7;
t3 = pkin(5) * t172 - qJ(6) * t92 - qJD(6) * t212 + t8;
t11 = [(Ifges(5,5) * t173 + Ifges(5,3) * t247) * t576 + (Ifges(5,5) * t255 + Ifges(5,3) * t296) * t579 + (t110 * t255 + t164 * t173 - t21 * t296 - t247 * t82) * mrSges(5,2) + m(3) * (t341 * t490 - t342 * t364 - t349 * t355 + t351 * t354) + (-Ifges(5,6) * t576 - Ifges(5,4) * t580 - Ifges(5,2) * t582 - t144 / 0.2e1 + t635 / 0.2e1 - t671 * t584 + t642 * t594 - t672 * t592 + t653) * t172 + (t75 * mrSges(6,2) + t54 * mrSges(7,2) - t33 * mrSges(6,3) - t20 * mrSges(7,3) + t633 / 0.2e1 - t672 * t584 + t645 * t594 + t646 * t592) * t92 + (Ifges(5,1) * t173 + Ifges(5,5) * t247) * t580 + (Ifges(5,1) * t255 + Ifges(5,5) * t296) * t601 + ((-t364 * mrSges(3,3) + Ifges(3,5) * t568 + (-0.2e1 * t605 + 0.3e1 / 0.2e1 * Ifges(3,4) * t409) * t399) * t409 + (t398 * (Ifges(4,5) * t297 - Ifges(4,6) * t296 + Ifges(4,3) * t358) / 0.2e1 - t490 * mrSges(3,3) - Ifges(3,6) * t401 + (-0.2e1 * t606 - 0.3e1 / 0.2e1 * t549 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t409) * t399) * t405) * t460 + (-t75 * mrSges(6,1) - t54 * mrSges(7,1) + t34 * mrSges(6,3) + t24 * mrSges(7,3) + t634 / 0.2e1 + t642 * t584 + t643 * t594 + t645 * t592) * t91 + t341 * (-t401 * mrSges(3,2) + mrSges(3,3) * t508) - t342 * (mrSges(3,1) * t401 - mrSges(3,3) * t509) + m(7) * (t1 * t30 + t10 * t54 + t2 * t35 + t20 * t3 + t24 * t4 + t67 * t9) + m(6) * (t19 * t97 + t27 * t75 + t33 * t8 + t34 * t7 + t45 * t6 + t46 * t5) + m(5) * (t100 * t22 + t101 * t21 + t110 * t185 + t126 * t164 + t31 * t82 + t32 * t81) + m(4) * (t116 * t206 + t117 * t205 + t134 * t184 + t135 * t183 + t231 * t250 + t244 * t245) + t97 * t29 - t358 * t538 + t10 * t114 + t27 * t115 + t274 * (Ifges(4,4) * t248 - Ifges(4,2) * t247) / 0.2e1 + t100 * t119 + (-t116 * t296 - t117 * t297 - t183 * t248 - t184 * t247) * mrSges(4,3) + t67 * t28 + t35 * t49 + t46 * t50 + t30 * t47 + t45 * t48 + (Ifges(4,1) * t248 - Ifges(4,4) * t247) * t573 + t377 * t568 + (Ifges(5,4) * t173 + Ifges(5,6) * t247) * t582 + (Ifges(5,4) * t255 + Ifges(5,6) * t296) * t600 + t358 * t537 + t248 * t590 + t247 * t591 + t101 * t120 + t4 * t139 + t7 * t140 + t3 * t141 + t8 * t142 + t126 * t169 + t185 * t79 + t297 * t597 + t358 * t598 + t173 * t602 + t247 * t604 + t31 * t194 + t32 * t195 + t205 * t218 + t206 * t219 + (t419 * t409 + (t647 + (t648 + t415) * t398 + t458) * t405) * t488 + t134 * t242 + t135 * t243 + t245 * t171 + t81 * (mrSges(5,1) * t247 - mrSges(5,3) * t173) + t231 * (mrSges(4,1) * t247 + mrSges(4,2) * t248) + t250 * t223 + t296 * t63 / 0.2e1 + t22 * (mrSges(5,1) * t296 - mrSges(5,3) * t255) - t296 * t159 / 0.2e1 + t244 * (mrSges(4,1) * t296 + mrSges(4,2) * t297) + t354 * t348 - t355 * t347 + t255 * t611 + (Ifges(4,5) * t248 - Ifges(4,6) * t247) * t674 + (Ifges(4,4) * t297 - Ifges(4,2) * t296 + Ifges(4,6) * t358) * t675 + (Ifges(4,1) * t297 - Ifges(4,4) * t296 + Ifges(4,5) * t358) * t676 + t649 * t254 + t650 * t212 + t651 * t211; (-pkin(2) * t171 + (t244 * mrSges(4,2) - t117 * mrSges(4,3) + t597 + t618) * t404 + (-t63 / 0.2e1 + t159 / 0.2e1 - t244 * mrSges(4,1) + t527 / 0.2e1 - t153 / 0.2e1 + t152 / 0.2e1 + t116 * mrSges(4,3) + (-Ifges(4,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t241 + t451) * t408 + ((t271 / 0.2e1 + t520 / 0.2e1 + t517 / 0.2e1 + t590 - t617) * t408 + (-t522 / 0.2e1 - t519 / 0.2e1 - t516 / 0.2e1 + t531 / 0.2e1 + t529 / 0.2e1 + t523 / 0.2e1 + t604 + t591 - t614) * t404) * qJD(3)) * t398 + (t229 * t642 - t230 * t672 - t267 * t671) * t585 + (t229 * t645 + t230 * t646 - t267 * t672) * t593 + (t227 * t646 + t228 * t645 - t308 * t672) * t592 + (-t227 * t672 + t228 * t642 - t308 * t671) * t584 + (Ifges(5,5) * t268 - Ifges(5,6) * t267) * t577 + (t493 * t81 + t494 * t82) * mrSges(5,3) + (Ifges(5,1) * t268 - Ifges(5,4) * t267) * t581 + (Ifges(5,4) * t268 - Ifges(5,2) * t267) * t583 + (-t268 / 0.2e1 + t307 / 0.2e1) * t145 + (t598 + t237 / 0.2e1 - t236 / 0.2e1 + t449) * t400 + ((-t383 / 0.2e1 + t489 * t605 - t419) * t409 + ((t606 + t549 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t409) * t489 + (-qJD(2) + t387 / 0.2e1) * Ifges(3,6) + ((Ifges(4,5) * t404 + Ifges(4,6) * t408) * t398 * t553 + (t400 * t553 + t570) * Ifges(4,3) - t415) * t398 - t458) * t405) * t489 + t377 + (Ifges(5,5) * t307 - Ifges(5,6) * t308) * t576 + (Ifges(5,1) * t307 - Ifges(5,4) * t308) * t580 + (Ifges(5,4) * t307 - Ifges(5,2) * t308) * t582 + (-mrSges(5,3) * t22 + t657 + t670) * t360 + (Ifges(4,4) * t575 - t203 / 0.2e1 + t617 + t659) * t323 + t624 * t242 + t625 * t243 + (-pkin(2) * t244 * t398 + t116 * t365 + t117 * t363 + t183 * t625 + t184 * t624 - t231 * t249) * m(4) + t626 * t169 + t627 * t194 + (t110 * t344 + t164 * t626 + t21 * t262 + t22 * t261 + t627 * t82 + t628 * t81) * m(5) + t628 * t195 + t629 * t115 + t149 * t47 + t166 * t49 + t630 * t114 + t189 * t48 + t190 * t50 + (-mrSges(5,1) * t494 - mrSges(5,2) * t493) * t164 + (mrSges(7,2) * t494 + mrSges(7,3) * t495) * t24 + (mrSges(6,2) * t494 + mrSges(6,3) * t495) * t34 + (-mrSges(7,1) * t495 + mrSges(7,2) * t496) * t54 + (-mrSges(6,1) * t495 + mrSges(6,2) * t496) * t75 + (-mrSges(7,1) * t494 - mrSges(7,3) * t496) * t20 + (-mrSges(6,1) * t494 - mrSges(6,3) * t496) * t33 + t631 * t142 + t632 * t140 + (t189 * t6 + t19 * t257 + t190 * t5 + t33 * t631 + t34 * t632 + t629 * t75) * m(6) + t633 * (t227 / 0.2e1 - t230 / 0.2e1) + t634 * (t228 / 0.2e1 - t229 / 0.2e1) + (t144 - t635) * (t267 / 0.2e1 - t308 / 0.2e1) + t636 * t139 + t637 * t141 + (t1 * t149 + t166 * t2 + t20 * t637 + t217 * t9 + t24 * t636 + t54 * t630) * m(7) + t217 * t28 + (t227 * t645 + t228 * t643 + t308 * t642) * t594 + (t229 * t643 + t230 * t645 + t267 * t642) * t595 - t249 * t223 + t257 * t29 + t261 * t119 + t262 * t120 - t341 * mrSges(3,2) - t342 * mrSges(3,1) + t344 * t79 - t349 * t348 + t351 * t347 + t363 * t218 + t365 * t219 + t649 * t359 - t650 * t426 + t651 * t309 + (-Ifges(4,4) * t574 - t143 / 0.2e1 + t202 / 0.2e1 + t614 + t654) * t322; (-pkin(3) * t110 + pkin(11) * t620 - t129 * t81 - t130 * t82 - t164 * t184) * m(5) + (-t519 + t143) * t574 + (-t671 * t585 - t672 * t593 + t642 * t595 + t603 + t655 - t635 / 0.2e1) * t512 + (-t1 * t503 - t2 * t505 + t20 * t216 - t215 * t24) * mrSges(7,3) + (-t34 * t215 + t33 * t216 - t5 * t505 - t6 * t503) * mrSges(6,3) + t215 * t685 + t449 + ((t602 + t530 / 0.2e1 + (-m(5) * t81 + m(6) * t75 - t497) * pkin(11) - t678) * qJD(4) + t120 * pkin(11) + t686) * t407 - t107 * t115 + t546 * t601 + t547 * t600 + t620 * mrSges(5,3) + t158 - pkin(3) * t79 + t202 * t573 + (t215 * t642 - t216 * t672) * t585 + (t215 * t643 + t216 * t645) * t595 + (t215 * t645 + t216 * t646) * t593 - t639 * t505 / 0.2e1 - (t145 * t274 + t640) * t407 / 0.2e1 - t633 * t216 / 0.2e1 + t447 * t535 + (-t164 * (mrSges(5,1) * t403 + mrSges(5,2) * t407) + (t403 * t82 + t407 * t81) * mrSges(5,3) + (Ifges(5,5) * t407 - Ifges(5,6) * t403) * t577 + (Ifges(5,1) * t407 - t547) * t581 + (-Ifges(5,2) * t403 + t546) * t583 - t669 + t659) * t274 + t664 * t114 + t665 * t140 + t666 * t142 + (pkin(11) * t535 - t107 * t75 + t33 * t666 + t339 * t6 + t34 * t665 + t340 * t5) * m(6) + t667 * t139 + t668 * t141 + (t1 * t295 + t2 * t312 + t20 * t668 + t24 * t667 + t376 * t9 + t54 * t664) * m(7) + (t271 + t203) * t575 + (t654 + t660) * t275 + (-t169 + t243 + t550) * t184 + (-t242 + t551) * t183 + ((-t532 / 0.2e1 + (-m(5) * t82 - t194) * pkin(11) + t613) * qJD(4) + (-t119 + t29) * pkin(11) + t9 * t445 + (t75 * t448 + t54 * t446 + (t20 * t402 - t24 * t406) * mrSges(7,3) + (t33 * t402 - t34 * t406) * mrSges(6,3) + t622 * t595 + t621 * t593 + t623 * t585 + t633 * t567 + t406 * t685) * qJD(5) + t661 * t610 + t662 * t609 + t663 * t599 + t657) * t403 - t130 * t194 - t129 * t195 - t75 * (-mrSges(6,1) * t215 + mrSges(6,2) * t216) - t54 * (-mrSges(7,1) * t215 + mrSges(7,2) * t216) + t295 * t47 + t312 * t49 + t339 * t48 + t340 * t50 + t376 * t28 + t503 * t673; -t613 * t239 + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t239 - t145 / 0.2e1 + t678) * t238 + m(7) * (t1 * t380 - t2 * t381 + t20 * t357 + t24 * t356 + t397 * t9) - t66 * t114 + t63 - m(7) * (t20 * t36 + t24 * t43 + t54 * t66) + ((m(7) * t54 + t114) * t561 - t615) * qJD(5) + t497 * t82 - pkin(4) * t29 - t451 + (-pkin(4) * t19 - t33 * t57 - t34 * t58 - t75 * t82) * m(6) + (m(6) * t619 + (-m(6) * t431 - t402 * t140 - t406 * t142) * qJD(5) - t402 * t48 + t406 * t50) * pkin(12) + t619 * mrSges(6,3) + t621 * t610 + t622 * t609 + t623 * t599 - t58 * t140 - t57 * t142 - t81 * t194 + t639 * t564 + (-t43 + t356) * t139 + t380 * t47 - t381 * t49 + t397 * t28 + t402 * t673 + (-t1 * t402 + t2 * t406) * mrSges(7,3) - t9 * t446 - t19 * t448 + (-t36 + t357) * t141; t640 + (-(-t20 + t23) * t24 + (-t192 * t54 + t1) * pkin(5)) * m(7) + (-t114 * t192 + t47) * pkin(5) + (t191 * t20 + t192 * t24) * mrSges(7,3) + (t191 * t33 + t192 * t34) * mrSges(6,3) - t23 * t139 - t33 * t140 + t24 * t141 + t34 * t142 - t75 * (mrSges(6,1) * t192 + mrSges(6,2) * t191) - t54 * (mrSges(7,1) * t192 + mrSges(7,2) * t191) + (t191 * t646 - t679) * t593 + t634 * t592 + (-t191 * t672 - t192 * t642) * t585 + (-t192 * t643 + t633 + t682) * t595 + t656; -t191 * t139 + t192 * t141 + 0.2e1 * (t9 / 0.2e1 + t24 * t595 + t20 * t592) * m(7) + t28;];
tauc  = t11(:);
