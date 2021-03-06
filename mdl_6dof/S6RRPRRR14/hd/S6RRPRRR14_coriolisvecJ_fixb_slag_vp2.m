% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:57:27
% EndTime: 2019-03-09 14:58:59
% DurationCPUTime: 46.76s
% Computational Cost: add. (72801->1039), mult. (229239->1501), div. (0->0), fcn. (199291->16), ass. (0->466)
t425 = sin(pkin(14));
t428 = sin(pkin(6));
t429 = cos(pkin(14));
t440 = cos(qJ(2));
t431 = cos(pkin(7));
t436 = sin(qJ(2));
t558 = t431 * t436;
t456 = (-t425 * t440 - t429 * t558) * t428;
t373 = qJD(1) * t456;
t455 = (-t425 * t558 + t429 * t440) * t428;
t376 = qJD(1) * t455;
t435 = sin(qJ(4));
t439 = cos(qJ(4));
t426 = sin(pkin(8));
t427 = sin(pkin(7));
t539 = qJD(1) * t428;
t521 = t436 * t539;
t504 = t427 * t521;
t481 = t426 * t504;
t430 = cos(pkin(8));
t559 = t430 * t439;
t270 = -t373 * t559 + t376 * t435 - t439 * t481;
t560 = t430 * t435;
t569 = t426 * t435;
t348 = t431 * t569 + (t425 * t439 + t429 * t560) * t427;
t336 = t348 * qJD(4);
t546 = t270 - t336;
t453 = t373 * t430 + t481;
t271 = t376 * t439 + t435 * t453;
t567 = t426 * t439;
t690 = t427 * (-t425 * t435 + t429 * t559) + t431 * t567;
t335 = t690 * qJD(4);
t545 = t271 - t335;
t432 = cos(pkin(6));
t619 = pkin(1) * t432;
t421 = t440 * t619;
t413 = qJD(1) * t421;
t509 = qJ(3) * t431 + pkin(10);
t496 = t509 * t436;
t467 = t428 * t496;
t354 = -qJD(1) * t467 + t413;
t420 = t436 * t619;
t563 = t428 * t440;
t675 = t509 * t563 + t420;
t355 = qJD(1) * t675;
t577 = qJ(3) * t427;
t466 = pkin(2) * t436 - t440 * t577;
t383 = t466 * t539;
t562 = t429 * t431;
t566 = t427 * t429;
t257 = -t354 * t425 - t355 * t562 + t383 * t566;
t617 = pkin(11) * t430;
t219 = pkin(3) * t504 - t376 * t617 + t257;
t287 = t355 * t427 + t431 * t383;
t618 = pkin(11) * t426;
t249 = -pkin(3) * t373 - t376 * t618 + t287;
t148 = -t219 * t426 + t430 * t249;
t537 = qJD(3) * t427;
t573 = t425 * t426;
t702 = t537 * t573 - t148;
t571 = t425 * t431;
t572 = t425 * t427;
t258 = t429 * t354 - t355 * t571 + t383 * t572;
t217 = pkin(11) * t453 + t258;
t401 = pkin(2) * t571 + qJ(3) * t566;
t338 = (t426 * t431 + t430 * t566) * pkin(11) + t401;
t417 = pkin(2) * t562;
t351 = pkin(3) * t431 + t417 + (-qJ(3) - t617) * t572;
t378 = (-pkin(3) * t429 - pkin(11) * t573 - pkin(2)) * t427;
t470 = t351 * t430 + t378 * t426;
t253 = -t435 * t338 + t470 * t439;
t551 = -t439 * t217 - t219 * t560 - t249 * t569 + (-t425 * t560 + t429 * t439) * t537 + t253 * qJD(4);
t313 = -t373 * t426 + t430 * t504;
t701 = pkin(12) * t313 - t551;
t700 = -pkin(4) * t546 + t545 * pkin(12) + t702;
t418 = qJD(1) * t432 + qJD(2);
t557 = t431 * t440;
t561 = t429 * t436;
t321 = t418 * t572 + (t425 * t557 + t561) * t539;
t457 = (-t425 * t436 + t429 * t557) * t428;
t320 = qJD(1) * t457 + t418 * t566;
t520 = t440 * t539;
t381 = t431 * t418 - t427 * t520;
t472 = t320 * t430 + t381 * t426;
t248 = t321 * t439 + t435 * t472;
t374 = qJD(2) * t456;
t364 = qJD(1) * t374;
t375 = qJD(2) * t455;
t365 = qJD(1) * t375;
t538 = qJD(2) * t428;
t510 = qJD(1) * t538;
t501 = t436 * t510;
t469 = t427 * t501;
t462 = t426 * t469;
t171 = qJD(4) * t248 - t364 * t559 + t365 * t435 - t439 * t462;
t640 = -t171 / 0.2e1;
t328 = t439 * t338;
t550 = -(t425 * t559 + t429 * t435) * t537 - (t435 * t470 + t328) * qJD(4) + t435 * t217 - t439 * (t219 * t430 + t249 * t426);
t434 = sin(qJ(5));
t438 = cos(qJ(5));
t234 = t271 * t434 - t438 * t313;
t396 = -t426 * t566 + t430 * t431;
t291 = t348 * t438 + t396 * t434;
t260 = qJD(5) * t291 + t335 * t434;
t699 = t234 - t260;
t235 = t271 * t438 + t313 * t434;
t290 = t348 * t434 - t438 * t396;
t259 = -qJD(5) * t290 + t335 * t438;
t547 = t235 - t259;
t574 = t418 * t427;
t382 = t431 * t520 + t574;
t468 = (-t427 ^ 2 - t431 ^ 2) * t521;
t292 = -t382 * t425 + t429 * t468;
t293 = t382 * t429 + t425 * t468;
t245 = t292 * t560 + t293 * t439;
t534 = qJD(4) * t439;
t518 = t426 * t534;
t698 = -t245 + t518;
t300 = -t364 * t426 + t430 * t469;
t410 = qJD(2) * t413;
t447 = (-qJD(2) * t496 + qJD(3) * t557) * t428;
t288 = qJD(1) * t447 + t418 * t537 + t410;
t536 = qJD(3) * t436;
t318 = -t428 * t431 * t536 - qJD(2) * t675;
t307 = t318 * qJD(1);
t344 = (qJD(2) * t466 - t427 * t536) * t428;
t334 = qJD(1) * t344;
t210 = t429 * t288 + t307 * t571 + t334 * t572;
t449 = t364 * t430 + t462;
t166 = pkin(11) * t449 + t210;
t209 = -t288 * t425 + t307 * t562 + t334 * t566;
t172 = pkin(3) * t469 - t365 * t617 + t209;
t317 = qJ(3) * t574 + t355;
t322 = pkin(2) * t418 + t354;
t379 = (-pkin(2) * t440 - t436 * t577 - pkin(1)) * t428;
t370 = qJD(1) * t379;
t238 = t429 * t317 + t322 * t571 + t370 * t572;
t179 = pkin(11) * t472 + t238;
t237 = -t317 * t425 + t322 * t562 + t370 * t566;
t182 = pkin(3) * t381 - t321 * t617 + t237;
t272 = -t322 * t427 + t431 * t370 + qJD(3);
t206 = -pkin(3) * t320 - t321 * t618 + t272;
t265 = -t307 * t427 + t431 * t334;
t212 = -pkin(3) * t364 - t365 * t618 + t265;
t535 = qJD(4) * t435;
t517 = t430 * t535;
t519 = t426 * t535;
t41 = t439 * (t172 * t430 + t212 * t426) - t435 * t166 - t179 * t534 - t182 * t517 - t206 * t519;
t39 = -pkin(4) * t300 - t41;
t639 = t171 / 0.2e1;
t247 = -t321 * t435 + t439 * t472;
t170 = qJD(4) * t247 + t365 * t439 + t435 * t449;
t278 = -t320 * t426 + t381 * t430 + qJD(4);
t181 = t248 * t438 + t278 * t434;
t119 = qJD(5) * t181 + t170 * t434 - t438 * t300;
t647 = -t119 / 0.2e1;
t180 = -t248 * t434 + t278 * t438;
t118 = qJD(5) * t180 + t170 * t438 + t300 * t434;
t648 = t118 / 0.2e1;
t478 = t182 * t430 + t206 * t426;
t89 = t179 * t439 + t435 * t478;
t84 = pkin(12) * t278 + t89;
t130 = -t182 * t426 + t430 * t206;
t86 = -pkin(4) * t247 - pkin(12) * t248 + t130;
t36 = t434 * t86 + t438 * t84;
t516 = t430 * t534;
t40 = t439 * t166 + t172 * t560 - t179 * t535 + t182 * t516 + t206 * t518 + t212 * t569;
t38 = pkin(12) * t300 + t40;
t129 = -t172 * t426 + t430 * t212;
t76 = pkin(4) * t171 - pkin(12) * t170 + t129;
t8 = -qJD(5) * t36 - t38 * t434 + t438 * t76;
t665 = t39 * mrSges(6,2) - t8 * mrSges(6,3) + 0.2e1 * Ifges(6,1) * t648 + 0.2e1 * Ifges(6,4) * t647 + 0.2e1 * Ifges(6,5) * t639;
t246 = qJD(5) - t247;
t646 = t119 / 0.2e1;
t281 = -t351 * t426 + t430 * t378;
t233 = -pkin(4) * t690 - pkin(12) * t348 + t281;
t254 = t351 * t560 + t378 * t569 + t328;
t240 = pkin(12) * t396 + t254;
t532 = qJD(5) * t438;
t533 = qJD(5) * t434;
t582 = t233 * t532 - t240 * t533 + t700 * t434 - t438 * t701;
t549 = pkin(4) * t313 - t550;
t34 = pkin(13) * t246 + t36;
t433 = sin(qJ(6));
t437 = cos(qJ(6));
t88 = -t435 * t179 + t439 * t478;
t83 = -pkin(4) * t278 - t88;
t52 = -pkin(5) * t180 - pkin(13) * t181 + t83;
t16 = -t34 * t433 + t437 * t52;
t17 = t34 * t437 + t433 * t52;
t35 = -t434 * t84 + t438 * t86;
t33 = -pkin(5) * t246 - t35;
t489 = Ifges(7,5) * t437 - Ifges(7,6) * t433;
t606 = Ifges(7,4) * t437;
t491 = -Ifges(7,2) * t433 + t606;
t607 = Ifges(7,4) * t433;
t493 = Ifges(7,1) * t437 - t607;
t494 = mrSges(7,1) * t433 + mrSges(7,2) * t437;
t620 = t437 / 0.2e1;
t621 = -t433 / 0.2e1;
t178 = qJD(6) - t180;
t636 = t178 / 0.2e1;
t135 = t181 * t437 + t246 * t433;
t642 = t135 / 0.2e1;
t134 = -t181 * t433 + t246 * t437;
t644 = t134 / 0.2e1;
t608 = Ifges(7,4) * t135;
t70 = Ifges(7,2) * t134 + Ifges(7,6) * t178 + t608;
t133 = Ifges(7,4) * t134;
t71 = Ifges(7,1) * t135 + Ifges(7,5) * t178 + t133;
t697 = -t33 * t494 - t489 * t636 - t491 * t644 - t493 * t642 - t620 * t71 - t621 * t70 + (t16 * t437 + t17 * t433) * mrSges(7,3);
t696 = Ifges(6,6) * t640 - t118 * Ifges(6,4) / 0.2e1;
t625 = t364 / 0.2e1;
t624 = t365 / 0.2e1;
t695 = t381 / 0.2e1;
t694 = pkin(13) * t546 - t582;
t693 = -pkin(5) * t699 + pkin(13) * t547 + t549;
t244 = -t292 * t559 + t293 * t435;
t692 = -t244 + t519;
t402 = -t438 * t430 + t434 * t569;
t542 = t292 * t426 * t434 - qJD(5) * t402 + t438 * t698;
t605 = Ifges(6,5) * t246;
t177 = Ifges(6,4) * t180;
t612 = Ifges(6,1) * t181;
t101 = t177 + t605 + t612;
t669 = t35 * mrSges(6,3) - t101 / 0.2e1 - t83 * mrSges(6,2);
t450 = -t605 / 0.2e1 + t669;
t691 = t450 + t697;
t137 = mrSges(6,1) * t246 - mrSges(6,3) * t181;
t82 = -mrSges(7,1) * t134 + mrSges(7,2) * t135;
t522 = m(7) * t33 + t82;
t681 = -m(6) * t35 - t137 + t522;
t689 = -pkin(12) * qJD(6) * t438 - t89 + t246 * (pkin(5) * t434 - pkin(13) * t438);
t409 = -pkin(5) * t438 - pkin(13) * t434 - pkin(4);
t160 = pkin(4) * t248 - pkin(12) * t247;
t73 = t434 * t160 + t438 * t88;
t688 = pkin(12) * t533 + pkin(13) * t248 - qJD(6) * t409 + t73;
t56 = -qJD(6) * t135 - t118 * t433 + t171 * t437;
t53 = Ifges(7,6) * t56;
t55 = qJD(6) * t134 + t118 * t437 + t171 * t433;
t54 = Ifges(7,5) * t55;
t13 = Ifges(7,3) * t119 + t53 + t54;
t18 = pkin(5) * t119 - pkin(13) * t118 + t39;
t7 = t438 * t38 + t434 * t76 + t86 * t532 - t533 * t84;
t5 = pkin(13) * t171 + t7;
t1 = qJD(6) * t16 + t18 * t433 + t437 * t5;
t2 = -qJD(6) * t17 + t18 * t437 - t433 * t5;
t500 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t687 = t13 / 0.2e1 + Ifges(6,2) * t646 - t500 + t39 * mrSges(6,1) - t7 * mrSges(6,3) + t696;
t641 = t170 / 0.2e1;
t626 = t300 / 0.2e1;
t686 = t469 / 0.2e1;
t685 = Ifges(4,3) * t695;
t684 = -Ifges(3,6) * t418 / 0.2e1;
t565 = t427 * t432;
t464 = t428 * t557 + t565;
t540 = pkin(10) * t563 + t420;
t342 = qJ(3) * t464 + t540;
t353 = pkin(2) * t432 + t421 - t467;
t255 = -t342 * t425 + t353 * t562 + t379 * t566;
t346 = t425 * t464 + t428 * t561;
t397 = -t427 * t563 + t432 * t431;
t207 = pkin(3) * t397 - t346 * t617 + t255;
t283 = -t353 * t427 + t431 * t379;
t345 = t429 * t565 + t457;
t232 = -pkin(3) * t345 - t346 * t618 + t283;
t143 = -t207 * t426 + t430 * t232;
t575 = t346 * t435;
t261 = -t345 * t559 - t397 * t567 + t575;
t471 = t345 * t430 + t397 * t426;
t262 = t346 * t439 + t435 * t471;
t104 = pkin(4) * t261 - pkin(12) * t262 + t143;
t289 = -t345 * t426 + t397 * t430;
t256 = t429 * t342 + t353 * t571 + t379 * t572;
t203 = pkin(11) * t471 + t256;
t524 = t439 * t203 + t207 * t560 + t232 * t569;
t96 = pkin(12) * t289 + t524;
t683 = t434 * t104 + t438 * t96;
t682 = t434 * t233 + t438 * t240;
t21 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t6 = -pkin(5) * t171 - t8;
t527 = m(7) * t6 + t21;
t80 = mrSges(6,1) * t171 - mrSges(6,3) * t118;
t679 = -m(6) * t8 + t527 - t80;
t678 = Ifges(4,5) * t624 + Ifges(4,6) * t625;
t677 = -t435 * t203 + t439 * (t207 * t430 + t232 * t426);
t515 = t436 * t538;
t503 = t427 * t515;
t480 = t426 * t503;
t452 = t374 * t430 + t480;
t414 = qJD(2) * t421;
t299 = t432 * t537 + t414 + t447;
t523 = t429 * t299 + t318 * t571 + t344 * t572;
t188 = pkin(11) * t452 + t523;
t497 = -t299 * t425 + t318 * t562 + t344 * t566;
t191 = pkin(3) * t503 - t375 * t617 + t497;
t508 = -t318 * t427 + t431 * t344;
t226 = -pkin(3) * t374 - t375 * t618 + t508;
t674 = t439 * (t191 * t430 + t226 * t426) - t435 * t188 - t203 * t534 - t207 * t517 - t232 * t519;
t672 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t581 = -qJD(5) * t682 + t434 * t701 + t700 * t438;
t668 = t83 * mrSges(6,1) + t16 * mrSges(7,1) - t17 * mrSges(7,2) - t36 * mrSges(6,3);
t601 = Ifges(6,6) * t246;
t603 = Ifges(6,2) * t180;
t609 = Ifges(6,4) * t181;
t100 = t601 + t603 + t609;
t599 = Ifges(7,3) * t178;
t600 = Ifges(7,6) * t134;
t604 = Ifges(7,5) * t135;
t69 = t599 + t600 + t604;
t667 = t100 / 0.2e1 - t69 / 0.2e1 - t668;
t589 = t278 * Ifges(5,6);
t610 = Ifges(5,4) * t248;
t139 = t247 * Ifges(5,2) + t589 + t610;
t591 = t246 * Ifges(6,3);
t592 = t181 * Ifges(6,5);
t593 = t180 * Ifges(6,6);
t99 = t591 + t592 + t593;
t666 = t36 * mrSges(6,2) + t89 * mrSges(5,3) + t139 / 0.2e1 - t99 / 0.2e1 - t130 * mrSges(5,1) - t35 * mrSges(6,1);
t655 = t56 / 0.2e1;
t656 = t55 / 0.2e1;
t664 = -Ifges(6,4) * t648 + Ifges(7,5) * t656 - Ifges(6,2) * t647 - Ifges(6,6) * t639 + Ifges(7,6) * t655 + Ifges(7,3) * t646 + t687;
t44 = Ifges(6,5) * t118 - Ifges(6,6) * t119 + Ifges(6,3) * t171;
t663 = mrSges(5,1) * t129 - mrSges(5,3) * t40 + Ifges(6,5) * t648 + Ifges(6,6) * t647 + Ifges(6,3) * t639 + t44 / 0.2e1 + t672 + (-t626 - t300 / 0.2e1) * Ifges(5,6) + (-t640 + t639) * Ifges(5,2) + (-t641 - t170 / 0.2e1) * Ifges(5,4);
t14 = Ifges(7,4) * t55 + Ifges(7,2) * t56 + Ifges(7,6) * t119;
t661 = t14 / 0.2e1;
t15 = Ifges(7,1) * t55 + Ifges(7,4) * t56 + Ifges(7,5) * t119;
t660 = t15 / 0.2e1;
t654 = -t70 / 0.2e1;
t653 = pkin(1) * mrSges(3,1);
t652 = pkin(1) * mrSges(3,2);
t110 = Ifges(5,5) * t170 - Ifges(5,6) * t171 + Ifges(5,3) * t300;
t651 = t110 / 0.2e1;
t649 = Ifges(5,1) * t641 + Ifges(5,4) * t640 + Ifges(5,5) * t626;
t645 = -t134 / 0.2e1;
t643 = -t135 / 0.2e1;
t638 = -t177 / 0.2e1;
t637 = -t178 / 0.2e1;
t635 = t180 / 0.2e1;
t634 = t181 / 0.2e1;
t633 = t246 / 0.2e1;
t632 = t247 / 0.2e1;
t631 = t248 / 0.2e1;
t630 = Ifges(4,3) * t686 + t678;
t629 = Ifges(4,4) * t624 + Ifges(4,2) * t625 + Ifges(4,6) * t686;
t628 = Ifges(4,1) * t624 + Ifges(4,4) * t625 + Ifges(4,5) * t686;
t627 = t278 / 0.2e1;
t623 = -t381 / 0.2e1;
t622 = t432 / 0.2e1;
t616 = qJD(2) / 0.2e1;
t145 = -pkin(13) * t690 + t682;
t239 = -pkin(4) * t396 - t253;
t156 = pkin(5) * t290 - pkin(13) * t291 + t239;
t91 = t145 * t437 + t156 * t433;
t615 = -qJD(6) * t91 + t433 * t694 + t437 * t693;
t90 = -t145 * t433 + t156 * t437;
t614 = qJD(6) * t90 + t433 * t693 - t437 * t694;
t613 = pkin(5) * t546 - t581;
t611 = Ifges(3,4) * t436;
t590 = t278 * Ifges(5,5);
t580 = t433 * t689 - t437 * t688;
t579 = t433 * t688 + t437 * t689;
t568 = t426 * t438;
t564 = t428 * t436;
t556 = t433 * t438;
t555 = t437 * t438;
t125 = -mrSges(6,1) * t180 + mrSges(6,2) * t181;
t184 = mrSges(5,1) * t278 - mrSges(5,3) * t248;
t554 = t125 - t184;
t263 = -t291 * t433 - t437 * t690;
t157 = qJD(6) * t263 + t259 * t437 + t336 * t433;
t162 = t235 * t437 + t270 * t433;
t553 = t157 - t162;
t264 = t291 * t437 - t433 * t690;
t158 = -qJD(6) * t264 - t259 * t433 + t336 * t437;
t161 = -t235 * t433 + t270 * t437;
t552 = t158 - t161;
t403 = t430 * t434 + t435 * t568;
t362 = -t403 * t433 - t437 * t567;
t544 = qJD(6) * t362 + t433 * t692 + t437 * t542;
t465 = -t403 * t437 + t433 * t567;
t543 = qJD(6) * t465 - t433 * t542 + t437 * t692;
t81 = -mrSges(6,2) * t171 - mrSges(6,3) * t119;
t529 = m(6) * t7 + t81;
t136 = -mrSges(6,2) * t246 + mrSges(6,3) * t180;
t512 = m(6) * t36 + t136;
t285 = -t364 * mrSges(4,1) + t365 * mrSges(4,2);
t132 = -t191 * t426 + t430 * t226;
t392 = t540 * qJD(1);
t507 = t684 - (Ifges(3,2) * t440 + t611) * t539 / 0.2e1 - t392 * mrSges(3,3);
t498 = t1 * t437 - t2 * t433;
t495 = mrSges(7,1) * t437 - mrSges(7,2) * t433;
t492 = Ifges(7,1) * t433 + t606;
t490 = Ifges(7,2) * t437 + t607;
t488 = Ifges(7,5) * t433 + Ifges(7,6) * t437;
t43 = pkin(13) * t261 + t683;
t204 = t262 * t434 - t289 * t438;
t205 = t262 * t438 + t289 * t434;
t95 = -pkin(4) * t289 - t677;
t66 = pkin(5) * t204 - pkin(13) * t205 + t95;
t20 = t43 * t437 + t433 * t66;
t485 = -t43 * t433 + t437 * t66;
t72 = t160 * t438 - t434 * t88;
t147 = t205 * t437 + t261 * t433;
t146 = -t205 * t433 + t261 * t437;
t149 = t233 * t438 - t240 * t434;
t314 = -t374 * t426 + t430 * t503;
t451 = t439 * t188 + t191 * t560 - t203 * t535 + t207 * t516 + t226 * t569 + t232 * t518;
t48 = pkin(12) * t314 + t451;
t192 = qJD(4) * t262 - t374 * t559 + t375 * t435 - t439 * t480;
t193 = t375 * t439 + t452 * t435 + (t439 * t471 - t575) * qJD(4);
t79 = pkin(4) * t192 - pkin(12) * t193 + t132;
t461 = t104 * t532 + t434 * t79 + t438 * t48 - t533 * t96;
t391 = -pkin(10) * t521 + t413;
t412 = Ifges(3,4) * t520;
t458 = -t391 * mrSges(3,3) + t418 * Ifges(3,5) + Ifges(3,1) * t521 / 0.2e1 + t412 / 0.2e1;
t454 = qJD(2) * t540;
t446 = t237 * mrSges(4,1) - t238 * mrSges(4,2) + t321 * Ifges(4,5) + t320 * Ifges(4,6) + t685;
t445 = -pkin(4) * t314 - t674;
t444 = -t599 / 0.2e1 - t600 / 0.2e1 - t604 / 0.2e1 + t601 / 0.2e1 + t609 / 0.2e1 + t667;
t443 = t603 / 0.2e1 + t444;
t408 = Ifges(3,5) * t440 * t510;
t404 = -pkin(10) * t564 + t421;
t400 = -qJ(3) * t572 + t417;
t390 = -t418 * mrSges(3,2) + mrSges(3,3) * t520;
t389 = mrSges(3,1) * t418 - mrSges(3,3) * t521;
t387 = qJD(1) * t454;
t386 = -pkin(10) * t501 + t410;
t385 = pkin(12) * t555 + t409 * t433;
t384 = -pkin(12) * t556 + t409 * t437;
t312 = mrSges(4,1) * t469 - mrSges(4,3) * t365;
t311 = -mrSges(4,2) * t469 + mrSges(4,3) * t364;
t280 = mrSges(4,1) * t381 - mrSges(4,3) * t321;
t279 = -mrSges(4,2) * t381 + mrSges(4,3) * t320;
t268 = -mrSges(4,1) * t320 + mrSges(4,2) * t321;
t252 = t321 * Ifges(4,1) + t320 * Ifges(4,4) + t381 * Ifges(4,5);
t251 = t321 * Ifges(4,4) + t320 * Ifges(4,2) + t381 * Ifges(4,6);
t243 = Ifges(5,4) * t247;
t183 = -mrSges(5,2) * t278 + mrSges(5,3) * t247;
t159 = -mrSges(5,1) * t247 + mrSges(5,2) * t248;
t154 = t247 * t555 + t248 * t433;
t153 = -t247 * t556 + t248 * t437;
t152 = -mrSges(5,2) * t300 - mrSges(5,3) * t171;
t151 = mrSges(5,1) * t300 - mrSges(5,3) * t170;
t144 = pkin(5) * t690 - t149;
t140 = t248 * Ifges(5,1) + t243 + t590;
t138 = t248 * Ifges(5,5) + t247 * Ifges(5,6) + t278 * Ifges(5,3);
t128 = -qJD(5) * t204 + t193 * t438 + t314 * t434;
t127 = qJD(5) * t205 + t193 * t434 - t314 * t438;
t126 = pkin(5) * t181 - pkin(13) * t180;
t124 = mrSges(5,1) * t171 + mrSges(5,2) * t170;
t98 = mrSges(7,1) * t178 - mrSges(7,3) * t135;
t97 = -mrSges(7,2) * t178 + mrSges(7,3) * t134;
t68 = qJD(6) * t146 + t128 * t437 + t192 * t433;
t67 = -qJD(6) * t147 - t128 * t433 + t192 * t437;
t63 = mrSges(6,1) * t119 + mrSges(6,2) * t118;
t61 = -pkin(5) * t248 - t72;
t30 = -mrSges(7,2) * t119 + mrSges(7,3) * t56;
t29 = mrSges(7,1) * t119 - mrSges(7,3) * t55;
t28 = t126 * t433 + t35 * t437;
t27 = t126 * t437 - t35 * t433;
t22 = pkin(5) * t127 - pkin(13) * t128 + t445;
t9 = pkin(13) * t192 + t461;
t4 = -qJD(6) * t20 + t22 * t437 - t433 * t9;
t3 = qJD(6) * t485 + t22 * t433 + t437 * t9;
t10 = [-t679 * (t104 * t438 - t434 * t96) + ((Ifges(3,5) * t622 - t404 * mrSges(3,3) + (-0.2e1 * t652 + 0.3e1 / 0.2e1 * Ifges(3,4) * t440) * t428) * t440 + (-Ifges(3,6) * t432 - t540 * mrSges(3,3) + t427 * (Ifges(4,5) * t346 + Ifges(4,6) * t345 + Ifges(4,3) * t397) / 0.2e1 + (-0.2e1 * t653 - 0.3e1 / 0.2e1 * t611 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t440) * t428) * t436) * t510 + (m(3) * t540 - t432 * mrSges(3,2) + mrSges(3,3) * t563) * t386 + (m(5) * t88 + t184) * t674 + t681 * (qJD(5) * t683 + t434 * t48 - t438 * t79) + (m(5) * t41 + t151) * t677 + (-Ifges(6,4) * t634 + Ifges(7,5) * t642 - Ifges(6,2) * t635 - Ifges(6,6) * t633 + Ifges(7,6) * t644 + Ifges(7,3) * t636 - t667) * t127 + (-pkin(5) * t527 + t663) * t261 + t664 * t204 + t665 * t205 + (-Ifges(5,4) * t631 + Ifges(6,5) * t634 - Ifges(5,2) * t632 - Ifges(5,6) * t627 + Ifges(6,6) * t635 + Ifges(6,3) * t633 - t522 * pkin(5) - t666) * t192 + (Ifges(4,5) * t375 + Ifges(4,6) * t374) * t695 + (m(4) * t283 - mrSges(4,1) * t345 + mrSges(4,2) * t346) * t265 + t321 * (Ifges(4,1) * t375 + Ifges(4,4) * t374) / 0.2e1 + (m(6) * t83 + t125) * t445 + (m(4) * t255 + mrSges(4,1) * t397) * t209 + (Ifges(5,4) * t262 + Ifges(5,6) * t289) * t640 + (Ifges(5,4) * t193 + Ifges(5,6) * t314) * t632 + (Ifges(5,5) * t262 + Ifges(5,3) * t289) * t626 + (Ifges(5,5) * t193 + Ifges(5,3) * t314) * t627 + (Ifges(5,1) * t262 + Ifges(5,5) * t289) * t641 + (Ifges(5,1) * t193 + Ifges(5,5) * t314) * t631 + (Ifges(7,1) * t68 + Ifges(7,4) * t67) * t642 + (Ifges(6,1) * t634 + Ifges(6,4) * t635 + Ifges(6,5) * t633 - t669) * t128 + (m(6) * t39 + t63) * t95 + (Ifges(7,4) * t68 + Ifges(7,2) * t67) * t644 + (t129 * t262 + t130 * t193 - t289 * t40 - t314 * t89) * mrSges(5,2) + (m(7) * t1 + t30) * t20 + (m(4) * t256 - mrSges(4,2) * t397) * t210 + (m(7) * t3 + mrSges(7,3) * t67) * t17 + (Ifges(7,5) * t68 + Ifges(7,6) * t67) * t636 + (mrSges(7,2) * t6 - mrSges(7,3) * t2 + Ifges(7,1) * t656 + Ifges(7,4) * t655 + Ifges(7,5) * t646 + t660) * t147 + t374 * t251 / 0.2e1 + t272 * (-mrSges(4,1) * t374 + mrSges(4,2) * t375) + t375 * t252 / 0.2e1 + t88 * (mrSges(5,1) * t314 - mrSges(5,3) * t193) + t314 * t138 / 0.2e1 + t256 * t311 + t255 * t312 - (m(3) * t404 + mrSges(3,1) * t432 - mrSges(3,3) * t564) * t387 + (-mrSges(7,1) * t6 + mrSges(7,3) * t1 + Ifges(7,4) * t656 + Ifges(7,2) * t655 + Ifges(7,6) * t646 + t661) * t146 + t41 * (mrSges(5,1) * t289 - mrSges(5,3) * t262) + t283 * t285 + t408 * t622 + (Ifges(4,1) * t346 + Ifges(4,4) * t345 + Ifges(4,5) * t397) * t624 + (Ifges(4,4) * t346 + Ifges(4,2) * t345 + Ifges(4,6) * t397) * t625 + t346 * t628 + t345 * t629 + t397 * t630 + (m(4) * t238 + t279) * t523 + (m(5) * t40 + t152) * t524 + t193 * t140 / 0.2e1 + (m(5) * t129 + t124) * t143 + (m(3) * t392 + t390) * (-pkin(10) * t515 + t414) + t512 * t461 + (-t209 * t346 + t210 * t345 - t237 * t375 + t238 * t374) * mrSges(4,3) + (m(4) * t272 + t268) * t508 + t3 * t97 + t4 * t98 + (m(5) * t130 + t159) * t132 + t262 * t649 + t289 * t651 + t320 * (Ifges(4,4) * t375 + Ifges(4,2) * t374) / 0.2e1 + t67 * t70 / 0.2e1 + t68 * t71 / 0.2e1 + t33 * (-mrSges(7,1) * t67 + mrSges(7,2) * t68) + t529 * t683 + (t458 * t440 + (t684 + (t685 + t446) * t427 + t507) * t436) * t538 + (m(5) * t89 + t183) * t451 - (m(3) * t391 + t389) * t454 + (m(7) * t4 - mrSges(7,3) * t68) * t16 + (m(7) * t2 + t29) * t485 + (m(4) * t237 + t280) * t497; (t209 * mrSges(4,1) - t210 * mrSges(4,2) + t630 + t678) * t431 - t668 * t699 + (Ifges(5,1) * t348 + Ifges(5,5) * t396) * t641 + (-t348 * t41 + t545 * t88 + t546 * t89) * mrSges(5,3) + (t99 - t139) * (t336 / 0.2e1 - t270 / 0.2e1) + (t69 - t100) * (t260 / 0.2e1 - t234 / 0.2e1) + (Ifges(5,4) * t348 + Ifges(5,6) * t396) * t640 + t664 * t290 + t665 * t291 + (Ifges(5,5) * t348 + Ifges(5,3) * t396) * t626 + (Ifges(7,4) * t264 + Ifges(7,2) * t263) * t655 - t320 * (Ifges(4,4) * t376 + Ifges(4,2) * t373) / 0.2e1 - t663 * t690 + (t1 * t263 - t16 * t553 + t17 * t552 - t2 * t264) * mrSges(7,3) + (t36 * t546 - t547 * t83) * mrSges(6,2) + (Ifges(7,1) * t264 + Ifges(7,4) * t263) * t656 + (t209 * t400 + t210 * t401 - t237 * t257 - t238 * t258 - t272 * t287 + (-pkin(2) * t265 + (-t237 * t425 + t238 * t429) * qJD(3)) * t427) * m(4) + t408 + (t335 / 0.2e1 - t271 / 0.2e1) * t140 + ((-t412 / 0.2e1 + t539 * t652 - t458) * t440 + ((t653 + t611 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t440) * t539 + (t418 / 0.2e1 - qJD(2)) * Ifges(3,6) + ((Ifges(4,5) * t425 + Ifges(4,6) * t429) * t427 * t616 + (t431 * t616 + t623) * Ifges(4,3) - t446) * t427 - t507) * t436) * t539 + ((Ifges(4,4) * t425 + Ifges(4,2) * t429) * t625 + (Ifges(4,1) * t425 + Ifges(4,4) * t429) * t624 + t265 * (-mrSges(4,1) * t429 + mrSges(4,2) * t425) + t429 * t629 + t425 * t628 - pkin(2) * t285 + (-t209 * t425 + t210 * t429) * mrSges(4,3) + (t429 * t279 + (t159 * t426 - t280) * t425) * qJD(3)) * t427 + t613 * t82 + t614 * t97 + t615 * t98 + (t1 * t91 + t144 * t6 + t16 * t615 + t17 * t614 + t2 * t90 + t33 * t613) * m(7) + t401 * t311 + t400 * t312 - t391 * t390 + t392 * t389 - t386 * mrSges(3,2) - t387 * mrSges(3,1) - t272 * (-mrSges(4,1) * t373 + mrSges(4,2) * t376) - t376 * t252 / 0.2e1 + t581 * t137 + t582 * t136 - t373 * t251 / 0.2e1 - t248 * (Ifges(5,1) * t271 - Ifges(5,4) * t270 + Ifges(5,5) * t313) / 0.2e1 - t313 * t138 / 0.2e1 - t278 * (Ifges(5,5) * t271 - Ifges(5,6) * t270 + Ifges(5,3) * t313) / 0.2e1 - t247 * (Ifges(5,4) * t271 - Ifges(5,2) * t270 + Ifges(5,6) * t313) / 0.2e1 - t287 * t268 + (-t162 / 0.2e1 + t157 / 0.2e1) * t71 - t258 * t279 - t257 * t280 + t281 * t124 + t550 * t184 - t246 * (Ifges(6,5) * t235 - Ifges(6,6) * t234 + Ifges(6,3) * t270) / 0.2e1 - t180 * (Ifges(6,4) * t235 - Ifges(6,2) * t234 + Ifges(6,6) * t270) / 0.2e1 - t181 * (Ifges(6,1) * t235 - Ifges(6,4) * t234 + Ifges(6,5) * t270) / 0.2e1 + t551 * t183 + (-mrSges(7,1) * t552 + mrSges(7,2) * t553) * t33 + (-mrSges(6,1) * t546 + mrSges(6,3) * t547) * t35 + t549 * t125 + t6 * (-mrSges(7,1) * t263 + mrSges(7,2) * t264) + (t129 * t348 - t130 * t545 + t313 * t89 - t396 * t40) * mrSges(5,2) + (-t130 * t546 - t313 * t88 + t396 * t41) * mrSges(5,1) + t253 * t151 + t254 * t152 + t239 * t63 + (Ifges(6,5) * t259 - Ifges(6,6) * t260 + Ifges(6,3) * t336) * t633 + (Ifges(6,1) * t259 - Ifges(6,4) * t260 + Ifges(6,5) * t336) * t634 + (Ifges(6,4) * t259 - Ifges(6,2) * t260 + Ifges(6,6) * t336) * t635 + (Ifges(7,5) * t157 + Ifges(7,6) * t158 + Ifges(7,3) * t260) * t636 + (Ifges(7,5) * t162 + Ifges(7,6) * t161 + Ifges(7,3) * t234) * t637 + (Ifges(7,1) * t157 + Ifges(7,4) * t158 + Ifges(7,5) * t260) * t642 + (Ifges(7,1) * t162 + Ifges(7,4) * t161 + Ifges(7,5) * t234) * t643 + (Ifges(4,5) * t376 + Ifges(4,6) * t373) * t623 + (Ifges(5,5) * t335 - Ifges(5,6) * t336) * t627 + (Ifges(5,1) * t335 - Ifges(5,4) * t336) * t631 + (Ifges(5,4) * t335 - Ifges(5,2) * t336) * t632 + (-t161 / 0.2e1 + t158 / 0.2e1) * t70 - t148 * t159 + t149 * t80 + t144 * t21 + (t129 * t281 + t130 * t702 + t253 * t41 + t254 * t40 + t550 * t88 + t551 * t89) * m(5) + t91 * t30 + t90 * t29 + (Ifges(7,4) * t157 + Ifges(7,2) * t158 + Ifges(7,6) * t260) * t644 + (Ifges(7,4) * t162 + Ifges(7,2) * t161 + Ifges(7,6) * t234) * t645 + t348 * t649 + t396 * t651 - t321 * (Ifges(4,1) * t376 + Ifges(4,4) * t373) / 0.2e1 + (t149 * t8 + t239 * t39 + t35 * t581 + t36 * t582 + t549 * t83 + t682 * t7) * m(6) + t682 * t81 + t264 * t660 + t263 * t661 + (t259 / 0.2e1 - t235 / 0.2e1) * t101 + (Ifges(7,5) * t264 + Ifges(7,6) * t263) * t646 + (t237 * t376 - t238 * t373) * mrSges(4,3); t430 * t124 - t245 * t183 - t293 * t279 - t292 * t280 + t362 * t29 - t465 * t30 + t403 * t81 + t543 * t98 + t544 * t97 + (t21 - t80) * t402 - t554 * t244 + t542 * t136 + (t435 * t152 + t292 * t159 + (t151 - t63) * t439 + (t183 * t439 + t435 * t554) * qJD(4)) * t426 + (-t237 * t292 - t238 * t293 + t265) * m(4) + (t129 * t430 + t88 * t244 - t89 * t245 + (t130 * t292 + t40 * t435 + t41 * t439 + (-t435 * t88 + t439 * t89) * qJD(4)) * t426) * m(5) + (-t1 * t465 + t16 * t543 + t17 * t544 + t2 * t362 + t6 * t402) * m(7) + (-t244 * t83 - t8 * t402 + t7 * t403 + (-t39 * t439 + t535 * t83) * t426 + t542 * t36) * m(6) + t285 + t681 * (qJD(5) * t403 - t292 * t568 + t434 * t698); (t493 * t656 + t491 * t655 + t6 * t494 + t15 * t620 + t14 * t621 + t489 * t646 + (-t1 * t433 - t2 * t437) * mrSges(7,3) + t679 * pkin(12) + t443 * t247 + (t71 * t621 + t437 * t654 + t33 * t495 + t488 * t637 + t490 * t645 + t492 * t643 + (t16 * t433 - t17 * t437) * mrSges(7,3)) * qJD(6) + t665) * t434 + (-t54 / 0.2e1 - t53 / 0.2e1 + (-Ifges(7,3) / 0.2e1 - Ifges(6,2) / 0.2e1) * t119 + t529 * pkin(12) + (t638 - t612 / 0.2e1 + t450) * t247 - t687 - t696) * t438 + ((-pkin(12) * t512 - t443) * t434 + (t681 * pkin(12) + t612 / 0.2e1 + t177 / 0.2e1 - t691) * t438) * qJD(5) + (-t591 / 0.2e1 - t593 / 0.2e1 - t592 / 0.2e1 + t589 / 0.2e1 + t610 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t247 + t666) * t248 + (-t153 * t17 + t154 * t16) * mrSges(7,3) + t110 + (-pkin(4) * t39 - t35 * t72 - t36 * t73 - t83 * t89) * m(6) + (-t243 / 0.2e1 - t590 / 0.2e1 - t140 / 0.2e1 - t130 * mrSges(5,2) + t88 * mrSges(5,3)) * t247 + t385 * t30 + t579 * t98 + t580 * t97 + (t1 * t385 + t16 * t579 + t17 * t580 + t2 * t384 - t33 * t61) * m(7) + t384 * t29 - t554 * t89 + (Ifges(7,5) * t154 + Ifges(7,6) * t153) * t637 + (Ifges(7,1) * t154 + Ifges(7,4) * t153) * t643 - t88 * t183 - t33 * (-mrSges(7,1) * t153 + mrSges(7,2) * t154) - t154 * t71 / 0.2e1 - t72 * t137 - t73 * t136 + (Ifges(7,4) * t154 + Ifges(7,2) * t153) * t645 + t153 * t654 - t61 * t82 - pkin(4) * t63 - t40 * mrSges(5,2) + t41 * mrSges(5,1); t444 * t181 + t44 + t14 * t620 + t433 * t660 + (-t82 + t137) * t36 - t35 * t136 - t28 * t97 - t27 * t98 + (-t29 * t433 + t30 * t437) * pkin(13) + (((-m(7) * t16 - t98) * t437 + (-m(7) * t17 - t97) * t433) * pkin(13) - t697) * qJD(6) - pkin(5) * t21 + (t638 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t181 + t691) * t180 + t488 * t646 + t490 * t655 + t492 * t656 - t6 * t495 + t498 * mrSges(7,3) + (-pkin(5) * t6 + pkin(13) * t498 - t16 * t27 - t17 * t28 - t33 * t36) * m(7) + t672; t17 * t98 - t33 * (mrSges(7,1) * t135 + mrSges(7,2) * t134) + (Ifges(7,5) * t134 - Ifges(7,6) * t135) * t637 + (Ifges(7,1) * t134 - t608) * t643 + t70 * t642 - t16 * t97 + (t134 * t16 + t135 * t17) * mrSges(7,3) - t500 + t13 + (-Ifges(7,2) * t135 + t133 + t71) * t645;];
tauc  = t10(:);
