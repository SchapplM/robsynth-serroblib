% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRP12
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:41
% EndTime: 2019-03-10 03:01:43
% DurationCPUTime: 62.71s
% Computational Cost: add. (37524->1062), mult. (111875->1441), div. (0->0), fcn. (92347->12), ass. (0->454)
t393 = cos(pkin(6));
t380 = qJD(1) * t393 + qJD(2);
t391 = sin(pkin(6));
t392 = cos(pkin(7));
t397 = sin(qJ(2));
t400 = cos(qJ(3));
t514 = t397 * t400;
t396 = sin(qJ(3));
t401 = cos(qJ(2));
t515 = t396 * t401;
t418 = t392 * t515 + t514;
t407 = t418 * t391;
t390 = sin(pkin(7));
t526 = t390 * t396;
t270 = qJD(1) * t407 + t380 * t526;
t582 = -t270 / 0.2e1;
t697 = Ifges(5,3) * t582;
t569 = pkin(1) * t393;
t387 = t401 * t569;
t378 = qJD(1) * t387;
t485 = pkin(10) * t392 + pkin(9);
t460 = t391 * t485;
t424 = t397 * t460;
t298 = -qJD(1) * t424 + t378;
t386 = t397 * t569;
t406 = -t401 * t460 - t386;
t299 = t406 * qJD(1);
t567 = pkin(10) * t390;
t410 = t391 * (pkin(2) * t397 - t401 * t567);
t333 = qJD(1) * t410;
t381 = pkin(10) * t526;
t521 = t392 * t400;
t358 = pkin(2) * t521 - t381;
t522 = t392 * t396;
t628 = t358 * qJD(3) - t400 * t298 - t299 * t522 - t333 * t526;
t523 = t391 * t401;
t264 = t380 * t567 + (t485 * t523 + t386) * qJD(1);
t268 = pkin(2) * t380 + t298;
t329 = (-pkin(2) * t401 - t397 * t567 - pkin(1)) * t391;
t313 = qJD(1) * t329;
t427 = t268 * t392 + t313 * t390;
t176 = -t396 * t264 + t400 * t427;
t502 = qJD(1) * t391;
t483 = t401 * t502;
t324 = t392 * t380 - t390 * t483 + qJD(3);
t157 = -pkin(3) * t324 - t176;
t395 = sin(qJ(4));
t399 = cos(qJ(4));
t511 = t400 * t401;
t488 = t392 * t511;
t467 = t391 * t488;
t484 = t397 * t502;
t525 = t390 * t400;
t489 = t380 * t525;
t405 = -qJD(1) * t467 + t396 * t484 - t489;
t222 = -t268 * t390 + t392 * t313;
t516 = t396 * t397;
t420 = t488 - t516;
t635 = t420 * t502;
t269 = t489 + t635;
t155 = -pkin(3) * t269 - pkin(11) * t270 + t222;
t177 = t264 * t400 + t396 * t427;
t158 = pkin(11) * t324 + t177;
t78 = t399 * t155 - t395 * t158;
t79 = t395 * t155 + t399 * t158;
t431 = t395 * t79 + t399 * t78;
t543 = Ifges(5,6) * t395;
t548 = Ifges(5,5) * t399;
t696 = t431 * mrSges(5,3) - t157 * (mrSges(5,1) * t395 + mrSges(5,2) * t399) - t405 * (-t543 + t548) / 0.2e1;
t241 = -t299 * t390 + t392 * t333;
t419 = t392 * t514 + t515;
t319 = t419 * t502;
t417 = -t392 * t516 + t511;
t320 = t417 * t502;
t499 = qJD(3) * t390;
t695 = -pkin(3) * t319 + pkin(11) * t320 - t241 + (pkin(3) * t396 - pkin(11) * t400) * t499;
t464 = t390 * t484;
t694 = pkin(11) * t464 - t628;
t693 = qJD(4) - t269;
t229 = -t270 * t395 + t324 * t399;
t403 = (qJD(2) * t417 + qJD(3) * t420) * t391;
t497 = qJD(3) * t400;
t478 = t390 * t497;
t231 = qJD(1) * t403 + t380 * t478;
t500 = qJD(2) * t391;
t470 = qJD(1) * t500;
t462 = t397 * t470;
t425 = t390 * t462;
t147 = qJD(4) * t229 + t231 * t399 + t395 * t425;
t606 = t147 / 0.2e1;
t692 = Ifges(5,4) * t606;
t262 = t320 * t395 - t399 * t464;
t356 = t392 * t395 + t399 * t526;
t303 = qJD(4) * t356 + t395 * t478;
t505 = t262 - t303;
t263 = t320 * t399 + t395 * t464;
t355 = -t399 * t392 + t395 * t526;
t302 = -qJD(4) * t355 + t399 * t478;
t504 = t263 - t302;
t360 = pkin(2) * t522 + pkin(10) * t525;
t691 = -t360 * qJD(3) + t396 * t298;
t498 = qJD(3) * t396;
t481 = t390 * t498;
t685 = t319 - t481;
t341 = pkin(11) * t392 + t360;
t342 = (-pkin(3) * t400 - pkin(11) * t396 - pkin(2)) * t390;
t495 = qJD(4) * t399;
t496 = qJD(4) * t395;
t631 = -t341 * t496 + t342 * t495 + t395 * t695 - t694 * t399;
t630 = t299 * t521 - (-pkin(3) * t484 - t333 * t400) * t390 - t691;
t216 = pkin(3) * t270 - pkin(11) * t269;
t126 = t399 * t176 + t395 * t216;
t690 = pkin(11) * t496 + pkin(12) * t270 + t126;
t689 = -pkin(11) * qJD(5) * t399 - t177 + t693 * (pkin(4) * t395 - pkin(12) * t399);
t650 = Ifges(6,1) + Ifges(7,1);
t649 = -Ifges(6,4) + Ifges(7,5);
t648 = Ifges(7,4) + Ifges(6,5);
t688 = Ifges(6,6) - Ifges(7,6);
t678 = -Ifges(6,3) - Ifges(7,2);
t687 = pkin(12) * t685 - t631;
t686 = -pkin(4) * t505 + pkin(12) * t504 + t630;
t528 = t269 * t395;
t684 = t496 - t528;
t225 = qJD(5) - t229;
t394 = sin(qJ(5));
t398 = cos(qJ(5));
t457 = qJD(4) - t489;
t404 = -t457 + t635;
t74 = -pkin(12) * t404 + t79;
t230 = t270 * t399 + t324 * t395;
t91 = -pkin(4) * t229 - pkin(12) * t230 + t157;
t32 = t394 * t91 + t398 * t74;
t24 = qJ(6) * t225 + t32;
t184 = t394 * t230 + t398 * t404;
t185 = t398 * t230 - t394 * t404;
t73 = pkin(4) * t404 - t78;
t33 = t184 * pkin(5) - t185 * qJ(6) + t73;
t180 = Ifges(7,5) * t185;
t80 = Ifges(7,6) * t225 + Ifges(7,3) * t184 + t180;
t552 = Ifges(6,4) * t185;
t83 = -Ifges(6,2) * t184 + Ifges(6,6) * t225 + t552;
t668 = -t80 / 0.2e1 + t83 / 0.2e1;
t683 = t73 * mrSges(6,1) + t33 * mrSges(7,1) - t24 * mrSges(7,2) - t32 * mrSges(6,3) - t668;
t301 = t406 * qJD(2);
t285 = qJD(1) * t301;
t334 = qJD(2) * t410;
t325 = qJD(1) * t334;
t375 = qJD(2) * t378;
t411 = qJD(2) * t424;
t284 = -qJD(1) * t411 + t375;
t480 = t392 * t498;
t466 = -t264 * t497 - t268 * t480 - t396 * t284 - t313 * t481;
t107 = -t285 * t521 + (-pkin(3) * t462 - t325 * t400) * t390 - t466;
t402 = (qJD(2) * t419 + qJD(3) * t418) * t391;
t232 = qJD(1) * t402 + t380 * t481;
t586 = t232 / 0.2e1;
t587 = -t232 / 0.2e1;
t148 = qJD(4) * t230 + t231 * t395 - t399 * t425;
t604 = t148 / 0.2e1;
t605 = -t148 / 0.2e1;
t72 = qJD(5) * t185 + t394 * t147 - t398 * t232;
t610 = t72 / 0.2e1;
t611 = -t72 / 0.2e1;
t71 = -qJD(5) * t184 + t398 * t147 + t394 * t232;
t612 = t71 / 0.2e1;
t479 = t392 * t497;
t113 = -t264 * t498 + t268 * t479 + t400 * t284 + t285 * t522 + t313 * t478 + t325 * t526;
t106 = pkin(11) * t425 + t113;
t235 = -t285 * t390 + t392 * t325;
t129 = pkin(3) * t232 - pkin(11) * t231 + t235;
t20 = t399 * t106 + t395 * t129 + t155 * t495 - t158 * t496;
t18 = pkin(12) * t232 + t20;
t37 = pkin(4) * t148 - pkin(12) * t147 + t107;
t493 = qJD(5) * t398;
t494 = qJD(5) * t394;
t3 = t398 * t18 + t394 * t37 + t91 * t493 - t494 * t74;
t1 = qJ(6) * t148 + qJD(6) * t225 + t3;
t4 = -qJD(5) * t32 - t18 * t394 + t37 * t398;
t2 = -pkin(5) * t148 - t4;
t621 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t645 = -t148 * t678 + t648 * t71 - t688 * t72;
t682 = t678 * t604 - t621 - Ifges(6,6) * t611 - Ifges(7,6) * t610 - t107 * mrSges(5,1) + t692 - (t587 - t586) * Ifges(5,6) - (t604 - t605) * Ifges(5,2) - t612 * t648 - t645 / 0.2e1;
t374 = -pkin(4) * t399 - pkin(12) * t395 - pkin(3);
t675 = -t374 * t494 + t394 * t690 + t398 * t689;
t674 = t374 * t493 + t394 * t689 - t398 * t690;
t632 = -t341 * t495 - t342 * t496 + t694 * t395 + t399 * t695;
t304 = t356 * t394 + t398 * t525;
t218 = -qJD(5) * t304 + t302 * t398 + t394 * t481;
t221 = t263 * t398 + t319 * t394;
t507 = t218 - t221;
t490 = t394 * t525;
t219 = -qJD(5) * t490 + t302 * t394 + t356 * t493 - t398 * t481;
t220 = t263 * t394 - t398 * t319;
t506 = t219 - t220;
t31 = -t394 * t74 + t398 * t91;
t634 = qJD(6) - t31;
t23 = -pkin(5) * t225 + t634;
t661 = -t31 * mrSges(6,1) + t23 * mrSges(7,1) + t32 * mrSges(6,2) - t24 * mrSges(7,3);
t181 = Ifges(6,4) * t184;
t547 = Ifges(7,5) * t184;
t640 = t185 * t650 + t648 * t225 - t181 + t547;
t680 = t640 / 0.2e1;
t644 = t648 * t148 + t649 * t72 + t650 * t71;
t679 = t644 / 0.2e1;
t641 = -t184 * t688 + t185 * t648 - t225 * t678;
t340 = t381 + (-pkin(2) * t400 - pkin(3)) * t392;
t248 = pkin(4) * t355 - pkin(12) * t356 + t340;
t256 = t399 * t341 + t395 * t342;
t250 = -pkin(12) * t525 + t256;
t639 = t248 * t493 - t250 * t494 + t394 * t686 - t398 * t687;
t627 = t394 * t248 + t398 * t250;
t638 = -qJD(5) * t627 + t394 * t687 + t398 * t686;
t677 = qJ(6) * t684 - qJD(6) * t399 + t674;
t676 = -pkin(5) * t684 - t675;
t125 = -t395 * t176 + t216 * t399;
t104 = -pkin(4) * t270 - t125;
t519 = t394 * t399;
t209 = t269 * t519 - t398 * t270;
t513 = t398 * t399;
t210 = t269 * t513 + t270 * t394;
t436 = pkin(5) * t394 - qJ(6) * t398;
t423 = pkin(11) + t436;
t437 = pkin(5) * t398 + qJ(6) * t394;
t673 = -pkin(5) * t209 + qJ(6) * t210 - t104 + (qJD(5) * t437 - qJD(6) * t398) * t395 + t423 * t495;
t633 = pkin(4) * t685 - t632;
t672 = -t394 * t688 + t398 * t648;
t671 = t394 * t648 + t398 * t688;
t546 = Ifges(7,5) * t394;
t551 = Ifges(6,4) * t394;
t670 = t398 * t650 + t546 - t551;
t545 = Ifges(7,5) * t398;
t550 = Ifges(6,4) * t398;
t669 = t394 * t650 - t545 + t550;
t21 = -t395 * t106 + t129 * t399 - t155 * t496 - t158 * t495;
t19 = -pkin(4) * t232 - t21;
t9 = pkin(5) * t72 - qJ(6) * t71 - qJD(6) * t185 + t19;
t667 = mrSges(6,2) * t19 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t9 + Ifges(6,4) * t611 + Ifges(7,5) * t610 + t604 * t648 + t612 * t650 + t679;
t666 = -mrSges(5,3) * t20 - t682 - t692;
t11 = t71 * Ifges(7,5) + t148 * Ifges(7,6) + t72 * Ifges(7,3);
t14 = t71 * Ifges(6,4) - t72 * Ifges(6,2) + t148 * Ifges(6,6);
t665 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 - t14 / 0.2e1 + t11 / 0.2e1;
t664 = -t157 * mrSges(5,1) + t661;
t432 = t31 * t398 + t32 * t394;
t434 = t23 * t398 - t24 * t394;
t453 = mrSges(7,1) * t394 - mrSges(7,3) * t398;
t455 = mrSges(6,1) * t394 + mrSges(6,2) * t398;
t663 = t434 * mrSges(7,2) - t432 * mrSges(6,3) + t33 * t453 + t73 * t455;
t662 = t107 * mrSges(5,2) + Ifges(5,1) * t606 + Ifges(5,5) * t586;
t659 = -t73 * mrSges(6,2) - mrSges(7,2) * t23 + mrSges(6,3) * t31 + t33 * mrSges(7,3);
t593 = t225 / 0.2e1;
t597 = t185 / 0.2e1;
t599 = t184 / 0.2e1;
t600 = -t184 / 0.2e1;
t658 = Ifges(6,6) * t600 + Ifges(7,6) * t599 - t593 * t678 + t648 * t597;
t594 = -t225 / 0.2e1;
t598 = -t185 / 0.2e1;
t657 = Ifges(6,6) * t599 + Ifges(7,6) * t600 - t594 * t678 + t648 * t598;
t588 = t231 / 0.2e1;
t655 = t404 / 0.2e1;
t643 = -qJ(6) * t505 + qJD(6) * t355 + t639;
t642 = pkin(5) * t505 - t638;
t305 = t356 * t398 - t490;
t637 = pkin(5) * t506 - qJ(6) * t507 - qJD(6) * t305 + t633;
t503 = pkin(9) * t523 + t386;
t286 = (t390 * t393 + t392 * t523) * pkin(10) + t503;
t297 = pkin(2) * t393 + t387 - t424;
t199 = -t396 * t286 + t400 * (t297 * t392 + t329 * t390);
t354 = -t390 * t523 + t393 * t392;
t178 = -pkin(3) * t354 - t199;
t292 = t393 * t526 + t407;
t246 = t292 * t395 - t354 * t399;
t247 = t292 * t399 + t354 * t395;
t115 = pkin(4) * t246 - pkin(12) * t247 + t178;
t291 = t391 * t516 - t393 * t525 - t467;
t237 = -t297 * t390 + t392 * t329;
t172 = pkin(3) * t291 - pkin(11) * t292 + t237;
t200 = t400 * t286 + t297 * t522 + t329 * t526;
t179 = pkin(11) * t354 + t200;
t98 = t395 * t172 + t399 * t179;
t94 = pkin(12) * t291 + t98;
t636 = t394 * t115 + t398 * t94;
t629 = -(t299 * t392 + t333 * t390) * t400 + t691;
t574 = t394 / 0.2e1;
t575 = -t394 / 0.2e1;
t626 = t80 * t574 + t83 * t575;
t625 = t20 * t399 - t21 * t395;
t624 = t3 * t398 - t394 * t4;
t623 = t1 * t398 + t2 * t394;
t239 = t393 * t481 + t402;
t379 = qJD(2) * t387;
t300 = t379 - t411;
t130 = -t286 * t498 + t297 * t479 + t400 * t300 + t301 * t522 + t329 * t478 + t334 * t526;
t482 = t397 * t500;
t463 = t390 * t482;
t122 = pkin(11) * t463 + t130;
t240 = t393 * t478 + t403;
t242 = -t301 * t390 + t392 * t334;
t142 = pkin(3) * t239 - pkin(11) * t240 + t242;
t29 = t399 * t122 + t395 * t142 + t172 * t495 - t179 * t496;
t25 = pkin(12) * t239 + t29;
t465 = -t286 * t497 - t297 * t480 - t396 * t300 - t329 * t481;
t123 = -t301 * t521 + (-pkin(3) * t482 - t334 * t400) * t390 - t465;
t165 = qJD(4) * t247 + t240 * t395 - t399 * t463;
t166 = -qJD(4) * t246 + t240 * t399 + t395 * t463;
t51 = pkin(4) * t165 - pkin(12) * t166 + t123;
t8 = -qJD(5) * t636 - t25 * t394 + t398 * t51;
t622 = -t21 * mrSges(5,1) + t20 * mrSges(5,2);
t619 = mrSges(6,1) * t19 + mrSges(7,1) * t9 - Ifges(6,2) * t611 + Ifges(7,3) * t610 - t604 * t688 + t612 * t649 + t665;
t617 = -0.2e1 * pkin(1);
t64 = t147 * Ifges(5,1) - t148 * Ifges(5,4) + t232 * Ifges(5,5);
t613 = t64 / 0.2e1;
t151 = Ifges(4,5) * t231 - Ifges(4,6) * t232 + Ifges(4,3) * t425;
t603 = t151 / 0.2e1;
t602 = Ifges(4,1) * t588 + Ifges(4,4) * t587 + Ifges(4,5) * t425 / 0.2e1;
t592 = -t229 / 0.2e1;
t591 = t229 / 0.2e1;
t590 = -t230 / 0.2e1;
t589 = t230 / 0.2e1;
t584 = -t269 / 0.2e1;
t583 = t269 / 0.2e1;
t581 = t270 / 0.2e1;
t579 = -t324 / 0.2e1;
t577 = t390 / 0.2e1;
t576 = t393 / 0.2e1;
t573 = t395 / 0.2e1;
t572 = -t398 / 0.2e1;
t571 = t398 / 0.2e1;
t568 = pkin(2) * t390;
t562 = mrSges(4,3) * t269;
t561 = mrSges(4,3) * t270;
t560 = mrSges(6,3) * t184;
t559 = mrSges(6,3) * t185;
t558 = Ifges(3,4) * t397;
t557 = Ifges(4,4) * t396;
t556 = Ifges(4,4) * t400;
t555 = Ifges(5,4) * t230;
t554 = Ifges(5,4) * t395;
t553 = Ifges(5,4) * t399;
t549 = Ifges(3,5) * t401;
t544 = Ifges(3,6) * t380;
t541 = t113 * mrSges(4,2);
t114 = (t285 * t392 + t325 * t390) * t400 + t466;
t540 = t114 * mrSges(4,1);
t539 = t19 * t395;
t536 = t229 * mrSges(5,3);
t535 = t230 * mrSges(5,3);
t534 = t269 * Ifges(4,6);
t533 = t270 * Ifges(4,4);
t532 = t270 * Ifges(4,5);
t531 = t324 * Ifges(4,3);
t530 = t380 * Ifges(3,5);
t161 = pkin(4) * t230 - pkin(12) * t229;
t56 = t394 * t161 + t398 * t78;
t527 = t374 * t398;
t524 = t391 * t397;
t139 = Ifges(5,2) * t229 - Ifges(5,6) * t404 + t555;
t518 = t395 * t139;
t517 = t395 * t398;
t224 = Ifges(5,4) * t229;
t140 = Ifges(5,1) * t230 - Ifges(5,5) * t404 + t224;
t512 = t399 * t140;
t112 = mrSges(6,1) * t184 + mrSges(6,2) * t185;
t189 = -mrSges(5,1) * t404 - t535;
t510 = t112 - t189;
t134 = -mrSges(7,2) * t184 + mrSges(7,3) * t225;
t135 = -mrSges(6,2) * t225 - t560;
t509 = -t134 - t135;
t136 = mrSges(6,1) * t225 - t559;
t137 = -mrSges(7,1) * t225 + mrSges(7,2) * t185;
t508 = -t136 + t137;
t336 = pkin(11) * t513 + t394 * t374;
t62 = Ifges(5,5) * t147 - Ifges(5,6) * t148 + Ifges(5,3) * t232;
t468 = t499 / 0.2e1;
t43 = -t148 * mrSges(7,1) + t71 * mrSges(7,2);
t97 = t172 * t399 - t395 * t179;
t255 = -t395 * t341 + t342 * t399;
t249 = pkin(4) * t525 - t255;
t456 = mrSges(6,1) * t398 - mrSges(6,2) * t394;
t454 = mrSges(7,1) * t398 + mrSges(7,3) * t394;
t452 = Ifges(5,1) * t399 - t554;
t447 = -Ifges(5,2) * t395 + t553;
t446 = -Ifges(6,2) * t394 + t550;
t445 = Ifges(6,2) * t398 + t551;
t439 = Ifges(7,3) * t394 + t545;
t438 = -Ifges(7,3) * t398 + t546;
t40 = t115 * t398 - t394 * t94;
t55 = t161 * t398 - t394 * t78;
t206 = t247 * t398 + t291 * t394;
t205 = t247 * t394 - t291 * t398;
t182 = t248 * t398 - t250 * t394;
t30 = -t395 * t122 + t142 * t399 - t172 * t496 - t179 * t495;
t93 = -pkin(4) * t291 - t97;
t7 = t115 * t493 + t398 * t25 + t394 * t51 - t494 * t94;
t409 = Ifges(5,5) * t166 - Ifges(5,6) * t165 + Ifges(5,3) * t239;
t26 = -pkin(4) * t239 - t30;
t352 = t503 * qJD(2);
t377 = Ifges(3,4) * t483;
t373 = t470 * t549;
t370 = -pkin(4) - t437;
t359 = -pkin(9) * t524 + t387;
t353 = qJD(5) * t436 - qJD(6) * t394;
t351 = -pkin(9) * t482 + t379;
t348 = t503 * qJD(1);
t346 = -pkin(9) * t484 + t378;
t345 = -t380 * mrSges(3,2) + mrSges(3,3) * t483;
t344 = mrSges(3,1) * t380 - mrSges(3,3) * t484;
t343 = t423 * t395;
t338 = qJD(1) * t352;
t337 = -pkin(9) * t462 + t375;
t335 = -pkin(11) * t519 + t527;
t316 = -t527 + (pkin(11) * t394 + pkin(5)) * t399;
t315 = -qJ(6) * t399 + t336;
t310 = Ifges(3,1) * t484 + t377 + t530;
t309 = t544 + (Ifges(3,2) * t401 + t558) * t502;
t267 = Ifges(4,4) * t269;
t234 = mrSges(4,1) * t324 - t561;
t233 = -mrSges(4,2) * t324 + t562;
t215 = -mrSges(4,1) * t269 + mrSges(4,2) * t270;
t212 = -mrSges(4,2) * t425 - mrSges(4,3) * t232;
t211 = mrSges(4,1) * t425 - mrSges(4,3) * t231;
t197 = t270 * Ifges(4,1) + t324 * Ifges(4,5) + t267;
t196 = t269 * Ifges(4,2) + t324 * Ifges(4,6) + t533;
t195 = t531 + t532 + t534;
t188 = mrSges(5,2) * t404 + t536;
t186 = pkin(5) * t304 - qJ(6) * t305 + t249;
t164 = -pkin(5) * t355 - t182;
t163 = qJ(6) * t355 + t627;
t162 = mrSges(4,1) * t232 + mrSges(4,2) * t231;
t160 = -mrSges(5,1) * t229 + mrSges(5,2) * t230;
t152 = Ifges(4,4) * t231 - Ifges(4,2) * t232 + Ifges(4,6) * t425;
t138 = Ifges(5,5) * t230 + Ifges(5,6) * t229 - t404 * Ifges(5,3);
t131 = (t301 * t392 + t334 * t390) * t400 + t465;
t117 = -mrSges(5,2) * t232 - mrSges(5,3) * t148;
t116 = mrSges(5,1) * t232 - mrSges(5,3) * t147;
t111 = mrSges(7,1) * t184 - mrSges(7,3) * t185;
t110 = pkin(5) * t185 + qJ(6) * t184;
t89 = -qJD(5) * t205 + t166 * t398 + t239 * t394;
t88 = qJD(5) * t206 + t166 * t394 - t239 * t398;
t76 = mrSges(5,1) * t148 + mrSges(5,2) * t147;
t57 = t229 * t436 + t79;
t52 = pkin(5) * t205 - qJ(6) * t206 + t93;
t45 = -mrSges(7,2) * t72 + mrSges(7,3) * t148;
t44 = -mrSges(6,2) * t148 - mrSges(6,3) * t72;
t42 = mrSges(6,1) * t148 - mrSges(6,3) * t71;
t39 = -pkin(5) * t230 - t55;
t38 = qJ(6) * t230 + t56;
t35 = -pkin(5) * t246 - t40;
t34 = qJ(6) * t246 + t636;
t28 = mrSges(6,1) * t72 + mrSges(6,2) * t71;
t27 = mrSges(7,1) * t72 - mrSges(7,3) * t71;
t10 = pkin(5) * t88 - qJ(6) * t89 - qJD(6) * t206 + t26;
t6 = -pkin(5) * t165 - t8;
t5 = qJ(6) * t165 + qJD(6) * t246 + t7;
t12 = [(Ifges(5,5) * t247 + Ifges(5,3) * t291) * t586 + (-Ifges(6,2) * t600 + Ifges(7,3) * t599 - t593 * t688 + t649 * t597 + t683) * t88 + (Ifges(5,4) * t247 + Ifges(5,6) * t291) * t605 + (Ifges(5,4) * t166 + Ifges(5,6) * t239) * t591 + (Ifges(5,1) * t247 + Ifges(5,5) * t291) * t606 + (Ifges(5,1) * t166 + Ifges(5,5) * t239) * t589 + t636 * t44 + m(6) * (t19 * t93 + t26 * t73 + t3 * t636 + t31 * t8 + t32 * t7 + t4 * t40) + m(3) * (t337 * t503 - t338 * t359 - t346 * t352 + t348 * t351) + t97 * t116 + t98 * t117 + t10 * t111 + t26 * t112 + (-t113 * t291 - t114 * t292 - t176 * t240 - t177 * t239) * mrSges(4,3) + t457 * t409 / 0.2e1 + (t107 * t247 + t157 * t166 - t20 * t291 - t239 * t79) * mrSges(5,2) + (-t420 * t409 / 0.2e1 + ((-t359 * mrSges(3,3) + Ifges(3,5) * t576 + (mrSges(3,2) * t617 + 0.3e1 / 0.2e1 * Ifges(3,4) * t401) * t391) * t401 + (-t503 * mrSges(3,3) + (Ifges(4,5) * t292 - Ifges(4,6) * t291 + Ifges(4,3) * t354) * t577 - Ifges(3,6) * t393 + (mrSges(3,1) * t617 - 0.3e1 / 0.2e1 * t558 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t401) * t391) * t397) * qJD(2)) * t502 + t324 * (Ifges(4,5) * t240 - Ifges(4,6) * t239) / 0.2e1 + t93 * t28 + t52 * t27 + t34 * t45 + t35 * t43 + t40 * t42 + ((t310 / 0.2e1 - t346 * mrSges(3,3) + t530 / 0.2e1) * t401 + (-t309 / 0.2e1 - t348 * mrSges(3,3) - t544 / 0.2e1 + (t195 / 0.2e1 - t177 * mrSges(4,2) + t176 * mrSges(4,1) + t534 / 0.2e1 + t532 / 0.2e1 + t531 / 0.2e1) * t390) * t397) * t500 - t354 * t541 + t337 * (-t393 * mrSges(3,2) + mrSges(3,3) * t523) - t338 * (mrSges(3,1) * t393 - mrSges(3,3) * t524) + m(7) * (t1 * t34 + t10 * t33 + t2 * t35 + t23 * t6 + t24 * t5 + t52 * t9) + m(5) * (t107 * t178 + t123 * t157 + t20 * t98 + t21 * t97 + t29 * t79 + t30 * t78) + m(4) * (t113 * t200 + t114 * t199 + t130 * t177 + t131 * t176 + t222 * t242 + t235 * t237) + (Ifges(6,4) * t600 + Ifges(7,5) * t599 + t593 * t648 + t597 * t650 - t659 + t680) * t89 + t247 * t613 + t292 * t602 + t354 * t603 + (Ifges(4,1) * t240 - Ifges(4,4) * t239) * t581 + (Ifges(4,4) * t240 - Ifges(4,2) * t239) * t583 + (Ifges(4,4) * t292 - Ifges(4,2) * t291 + Ifges(4,6) * t354) * t587 + (Ifges(4,1) * t292 - Ifges(4,4) * t291 + Ifges(4,5) * t354) * t588 + t373 * t576 + t354 * t540 + t619 * t205 + t5 * t134 + t7 * t135 + t8 * t136 + t6 * t137 + t123 * t160 + t166 * t140 / 0.2e1 + t178 * t76 + t29 * t188 + t30 * t189 + t199 * t211 + t200 * t212 + t130 * t233 + t131 * t234 + t237 * t162 + t239 * t138 / 0.2e1 + t78 * (mrSges(5,1) * t239 - mrSges(5,3) * t166) - t239 * t196 / 0.2e1 + t240 * t197 / 0.2e1 + t222 * (mrSges(4,1) * t239 + mrSges(4,2) * t240) + t242 * t215 + t291 * t62 / 0.2e1 + t21 * (mrSges(5,1) * t291 - mrSges(5,3) * t247) - t291 * t152 / 0.2e1 + t235 * (mrSges(4,1) * t291 + mrSges(4,2) * t292) + t351 * t345 - t352 * t344 + (t641 / 0.2e1 - Ifges(5,2) * t591 - Ifges(5,4) * t589 - t139 / 0.2e1 - t79 * mrSges(5,3) + t658 - t664) * t165 + t666 * t246 + t667 * t206; (t324 * (Ifges(4,5) * t400 - Ifges(4,6) * t396) + t269 * (-Ifges(4,2) * t396 + t556) + t270 * (Ifges(4,1) * t400 - t557) + t400 * t197) * t468 + (-mrSges(5,1) * t685 + mrSges(5,3) * t504) * t78 + (mrSges(5,2) * t685 + mrSges(5,3) * t505) * t79 + (-t220 * t688 + t221 * t648 - t262 * t678) * t594 + (t218 * t648 - t219 * t688 - t303 * t678) * t593 + (Ifges(6,4) * t221 + Ifges(7,5) * t218 - Ifges(6,2) * t220 + Ifges(6,6) * t262 + Ifges(7,6) * t303 + Ifges(7,3) * t219) * t599 - ((-Ifges(3,2) * t484 + t310 + t377) * t401 + t380 * (-Ifges(3,6) * t397 + t549)) * t502 / 0.2e1 + (t182 * t4 + t19 * t249 + t3 * t627 + t31 * t638 + t32 * t639 + t633 * t73) * m(6) + t627 * t44 + (-t397 * (Ifges(3,1) * t401 - t558) / 0.2e1 + pkin(1) * (mrSges(3,1) * t397 + mrSges(3,2) * t401)) * qJD(1) ^ 2 * t391 ^ 2 + (Ifges(5,5) * t263 - Ifges(5,6) * t262 + Ifges(5,3) * t319) * t655 + (Ifges(6,4) * t218 + Ifges(7,5) * t221 - Ifges(6,2) * t219 + Ifges(6,6) * t303 + Ifges(7,6) * t262 + Ifges(7,3) * t220) * t600 + t628 * t233 + t629 * t234 + (t113 * t360 + t114 * t358 + t176 * t629 + t177 * t628 - t222 * t241 - t235 * t568) * m(4) + t630 * t160 + t631 * t188 + (t107 * t340 + t157 * t630 + t20 * t256 + t21 * t255 + t631 * t79 + t632 * t78) * m(5) + t632 * t189 + t633 * t112 + (t220 * t649 + t221 * t650 + t262 * t648) * t598 + (t218 * t650 + t219 * t649 + t303 * t648) * t597 + t373 + (t138 - t196) * (t396 * t468 - t319 / 0.2e1) + (t80 - t83) * (t219 / 0.2e1 - t220 / 0.2e1) + (t152 / 0.2e1 - t62 / 0.2e1 - Ifges(5,5) * t606 - Ifges(5,6) * t605 - Ifges(5,3) * t586 + t622) * t525 + (t176 * t320 + t177 * t319 + (t113 * t400 - t114 * t396 + (-t176 * t400 - t177 * t396) * qJD(3)) * t390) * mrSges(4,3) + t637 * t111 + t638 * t136 + t639 * t135 + t640 * (t218 / 0.2e1 - t221 / 0.2e1) + (t139 - t641) * (t262 / 0.2e1 - t303 / 0.2e1) + t642 * t137 + t643 * t134 + (t1 * t163 + t164 * t2 + t186 * t9 + t23 * t642 + t24 * t643 + t33 * t637) * m(7) - t195 * t464 / 0.2e1 + ((Ifges(4,3) * t392 + (Ifges(4,5) * t396 + Ifges(4,6) * t400) * t390) * t577 - Ifges(3,6)) * t462 - t162 * t568 + (-t263 / 0.2e1 + t302 / 0.2e1) * t140 + (mrSges(7,1) * t505 + mrSges(7,2) * t507) * t23 + (mrSges(6,1) * t506 + mrSges(6,2) * t507) * t73 + (-mrSges(6,1) * t505 - mrSges(6,3) * t507) * t31 + (mrSges(7,1) * t506 - mrSges(7,3) * t507) * t33 + (-mrSges(5,1) * t505 - mrSges(5,2) * t504) * t157 + (mrSges(6,2) * t505 - mrSges(6,3) * t506) * t32 + (-mrSges(7,2) * t506 - mrSges(7,3) * t505) * t24 + (-t113 * t392 - t222 * t320 + (t177 * t484 + t222 * t497 + t235 * t396) * t390) * mrSges(4,2) + (t114 * t392 - t222 * t319 + (-t176 * t484 + t222 * t498 - t235 * t400) * t390) * mrSges(4,1) + t309 * t484 / 0.2e1 - t404 * (Ifges(5,5) * t302 - Ifges(5,6) * t303 + Ifges(5,3) * t481) / 0.2e1 + t526 * t602 + t392 * t603 + (Ifges(5,4) * t302 - Ifges(5,2) * t303 + Ifges(5,6) * t481) * t591 + (Ifges(5,4) * t263 - Ifges(5,2) * t262 + Ifges(5,6) * t319) * t592 + (Ifges(4,1) * t320 - Ifges(4,4) * t319 + Ifges(4,5) * t464) * t582 + (Ifges(4,4) * t320 - Ifges(4,2) * t319 + Ifges(4,6) * t464) * t584 + (Ifges(4,6) * t392 + (Ifges(4,2) * t400 + t557) * t390) * t587 + (Ifges(4,5) * t392 + (Ifges(4,1) * t396 + t556) * t390) * t588 + (Ifges(5,1) * t302 - Ifges(5,4) * t303 + Ifges(5,5) * t481) * t589 + (Ifges(5,1) * t263 - Ifges(5,4) * t262 + Ifges(5,5) * t319) * t590 + (Ifges(4,5) * t320 - Ifges(4,6) * t319 + Ifges(4,3) * t464) * t579 + t619 * t304 + t163 * t45 + t164 * t43 + t182 * t42 + (t346 * t401 + t348 * t397) * mrSges(3,3) * t502 + t186 * t27 - t241 * t215 + t249 * t28 + t255 * t116 + t256 * t117 - t320 * t197 / 0.2e1 - t337 * mrSges(3,2) - t338 * mrSges(3,1) + t340 * t76 - t346 * t345 + t348 * t344 + t358 * t211 + t360 * t212 + (-mrSges(5,3) * t21 + Ifges(5,4) * t605 + t613 + t662) * t356 + t666 * t355 + t667 * t305; (Ifges(6,4) * t599 + Ifges(7,5) * t600 + t648 * t594 + t598 * t650 + t659 - t640 / 0.2e1) * t210 + (-Ifges(6,2) * t599 + Ifges(7,3) * t600 - t594 * t688 + t598 * t649 - t683) * t209 + (t665 * t394 + (t80 * t571 + t83 * t572 + t73 * t456 + t33 * t454 + t438 * t600 + t445 * t599 + (t31 * t394 - t32 * t398) * mrSges(6,3) + (-t23 * t394 - t24 * t398) * mrSges(7,2) + t669 * t598 + t671 * t594 + t640 * t575) * qJD(5) + t9 * t453 + (-pkin(11) * t188 + t658 - t661) * qJD(4) + (-t116 + t28) * pkin(11) + t439 * t610 + t446 * t611 + t670 * t612 + t672 * t604 + t662) * t395 + (((m(6) * t73 + t510) * pkin(11) + t626) * qJD(4) + t117 * pkin(11) + t682) * t399 + t625 * mrSges(5,3) + t554 * t605 + (-Ifges(4,2) * t270 + t197 + t267 + t512) * t584 + (-mrSges(4,1) * t222 - mrSges(5,1) * t78 + mrSges(5,2) * t79 + Ifges(5,5) * t590 - Ifges(4,6) * t579 + Ifges(5,6) * t592) * t270 + t553 * t606 + (-t4 * t517 + (-t31 * t513 - t32 * t519) * qJD(4)) * mrSges(6,3) + (-pkin(3) * t107 - t125 * t78 - t126 * t79 - t157 * t177 + (-qJD(4) * t431 + t625) * pkin(11)) * m(5) + (-t222 * mrSges(4,2) + Ifges(4,1) * t582 + Ifges(4,5) * t579 + t447 * t592 + t452 * t590 + t696) * t269 + (t138 - t533) * t582 - t104 * t112 + (t2 * t517 + (t23 * t513 - t24 * t519) * qJD(4)) * mrSges(7,2) + t151 - pkin(3) * t76 + t517 * t679 + (-t160 + t234 + t561) * t177 + (-t233 + t562) * t176 + t196 * t581 + t518 * t583 + t64 * t573 + t455 * t539 + t675 * t136 + (pkin(11) * t539 - t104 * t73 + t3 * t336 + t31 * t675 + t32 * t674 + t335 * t4) * m(6) + t676 * t137 + t677 * t134 + (t1 * t315 + t2 * t316 + t23 * t676 + t24 * t677 + t33 * t673 + t343 * t9) * m(7) + t540 - t541 - t126 * t188 - t125 * t189 + (t512 / 0.2e1 - t518 / 0.2e1 + t33 * (mrSges(7,1) * t519 - mrSges(7,3) * t513) + t73 * (mrSges(6,1) * t519 + mrSges(6,2) * t513) + (Ifges(6,4) * t513 - Ifges(6,2) * t519) * t600 + (Ifges(7,5) * t513 + Ifges(7,3) * t519) * t599 + t447 * t591 + t452 * t589 + t697 + (t513 * t650 + t519 * t649) * t597 + (t513 * t648 - t519 * t688) * t593 + t641 * t573 + t513 * t680 + t693 * (t548 / 0.2e1 - t543 / 0.2e1) - t696) * qJD(4) + t315 * t45 + t316 * t43 + t335 * t42 + t336 * t44 + t343 * t27 + t405 * t697 + t673 * t111 + (-t641 / 0.2e1 + t657 + t661) * t528 + t674 * t135; -t622 - m(7) * (t23 * t39 + t24 * t38 + t33 * t57) - m(6) * (t31 * t55 + t32 * t56 + t73 * t79) - t19 * t456 - t9 * t454 + t623 * mrSges(7,2) + t624 * mrSges(6,3) + m(7) * (pkin(12) * t623 + t33 * t353 + t370 * t9) + m(6) * (-pkin(4) * t19 + pkin(12) * t624) + ((t44 + t45) * t398 + (-t42 + t43) * t394) * pkin(12) - pkin(4) * t28 + (-t157 * mrSges(5,2) + Ifges(5,1) * t590 + Ifges(5,5) * t655 + t394 * t668 + t439 * t600 + t446 * t599 + t572 * t640 + t594 * t672 + t598 * t670 - t663) * t229 + (-t510 + t535) * t79 + (-t188 + t536) * t78 + t438 * t610 + t445 * t611 + t139 * t589 + t14 * t571 + t11 * t572 + (-t555 + t641) * t590 + t644 * t574 - t38 * t134 - t56 * t135 - t55 * t136 - t39 * t137 + t62 + (t224 + t140) * t592 + (-t57 + t353) * t111 + t370 * t27 + (-Ifges(5,2) * t592 - Ifges(5,6) * t655 + t657 + t664) * t230 + t669 * t612 + t671 * t604 + (t439 * t599 + t446 * t600 + (-m(6) * t432 + m(7) * t434 + t394 * t509 + t398 * t508) * pkin(12) + t670 * t597 + t672 * t593 + t640 * t571 + t626 + t663) * qJD(5); t621 - t110 * t111 + (-Ifges(6,2) * t185 - t181 + t640) * t599 + (-pkin(5) * t2 + qJ(6) * t1 - t110 * t33 - t23 * t32 + t24 * t634) * m(7) + (-t184 * t648 - t185 * t688) * t594 - pkin(5) * t43 + qJ(6) * t45 + (-t508 + t559) * t32 + (t509 - t560) * t31 + (t184 * t23 + t185 * t24) * mrSges(7,2) + (Ifges(7,3) * t185 - t547) * t600 + t83 * t597 + qJD(6) * t134 - t33 * (mrSges(7,1) * t185 + mrSges(7,3) * t184) + (-t184 * t650 + t180 - t552 + t80) * t598 - t73 * (mrSges(6,1) * t185 - mrSges(6,2) * t184) + t645; t185 * t111 - t225 * t134 + 0.2e1 * (t2 / 0.2e1 + t33 * t597 + t24 * t594) * m(7) + t43;];
tauc  = t12(:);
