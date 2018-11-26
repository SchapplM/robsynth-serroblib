% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:43:44
% EndTime: 2018-11-23 18:44:30
% DurationCPUTime: 46.89s
% Computational Cost: add. (62251->1121), mult. (180725->1593), div. (0->0), fcn. (152192->14), ass. (0->509)
t451 = cos(qJ(2));
t441 = cos(pkin(6));
t623 = pkin(1) * t441;
t433 = t451 * t623;
t425 = qJD(1) * t433;
t446 = sin(qJ(2));
t439 = sin(pkin(6));
t440 = cos(pkin(7));
t524 = pkin(10) * t440 + pkin(9);
t501 = t439 * t524;
t477 = t446 * t501;
t357 = -qJD(1) * t477 + t425;
t432 = t446 * t623;
t462 = -t451 * t501 - t432;
t358 = t462 * qJD(1);
t438 = sin(pkin(7));
t622 = pkin(10) * t438;
t465 = (pkin(2) * t446 - t451 * t622) * t439;
t387 = qJD(1) * t465;
t445 = sin(qJ(3));
t560 = t438 * t445;
t428 = pkin(10) * t560;
t450 = cos(qJ(3));
t555 = t440 * t450;
t408 = pkin(2) * t555 - t428;
t556 = t440 * t445;
t698 = t408 * qJD(3) - t450 * t357 - t358 * t556 - t387 * t560;
t295 = -t358 * t438 + t440 * t387;
t551 = t446 * t450;
t552 = t445 * t451;
t475 = t440 * t551 + t552;
t537 = qJD(1) * t439;
t375 = t475 * t537;
t548 = t450 * t451;
t553 = t445 * t446;
t473 = -t440 * t553 + t548;
t376 = t473 * t537;
t739 = -pkin(3) * t375 + pkin(11) * t376 - t295 + (pkin(3) * t445 - pkin(11) * t450) * t438 * qJD(3);
t522 = t446 * t537;
t505 = t438 * t522;
t738 = pkin(11) * t505 - t698;
t444 = sin(qJ(4));
t449 = cos(qJ(4));
t313 = -t376 * t444 + t449 * t505;
t406 = t440 * t444 + t449 * t560;
t534 = qJD(3) * t450;
t518 = t438 * t534;
t362 = -qJD(4) * t406 - t444 * t518;
t540 = t313 - t362;
t314 = t376 * t449 + t444 * t505;
t405 = t440 * t449 - t444 * t560;
t361 = qJD(4) * t405 + t449 * t518;
t539 = t314 - t361;
t535 = qJD(3) * t445;
t519 = t438 * t535;
t732 = t375 - t519;
t559 = t438 * t450;
t410 = pkin(2) * t556 + pkin(10) * t559;
t392 = pkin(11) * t440 + t410;
t393 = (-pkin(3) * t450 - pkin(11) * t445 - pkin(2)) * t438;
t532 = qJD(4) * t449;
t533 = qJD(4) * t444;
t702 = -t392 * t533 + t393 * t532 + t444 * t739 - t738 * t449;
t304 = t449 * t392 + t444 * t393;
t701 = -qJD(4) * t304 + t738 * t444 + t449 * t739;
t737 = -pkin(4) * t732 + pkin(12) * t539 + t701;
t736 = pkin(12) * t540 - t702;
t735 = -t410 * qJD(3) + t445 * t357;
t443 = sin(qJ(5));
t448 = cos(qJ(5));
t480 = t448 * t405 - t406 * t443;
t246 = qJD(5) * t480 + t361 * t448 + t362 * t443;
t258 = t313 * t443 + t314 * t448;
t734 = t246 - t258;
t321 = t405 * t443 + t406 * t448;
t247 = qJD(5) * t321 + t361 * t443 - t448 * t362;
t257 = -t448 * t313 + t314 * t443;
t733 = t247 - t257;
t521 = t451 * t537;
t427 = qJD(1) * t441 + qJD(2);
t561 = t427 * t438;
t323 = -t445 * t522 + t450 * (t440 * t521 + t561);
t414 = t443 * t449 + t444 * t448;
t259 = t414 * t323;
t690 = qJD(4) + qJD(5);
t354 = t690 * t414;
t542 = t259 - t354;
t413 = t443 * t444 - t448 * t449;
t260 = t413 * t323;
t353 = t690 * t413;
t541 = -t260 + t353;
t696 = t358 * t555 - (-pkin(3) * t522 - t387 * t450) * t438 - t735;
t303 = -t392 * t444 + t449 * t393;
t272 = -pkin(4) * t559 - pkin(12) * t406 + t303;
t289 = pkin(12) * t405 + t304;
t530 = qJD(5) * t448;
t531 = qJD(5) * t443;
t710 = t272 * t530 - t289 * t531 + t443 * t737 - t736 * t448;
t664 = -pkin(12) - pkin(11);
t523 = qJD(4) * t664;
t416 = t444 * t523;
t422 = t664 * t444;
t423 = t664 * t449;
t479 = t448 * t422 + t423 * t443;
t506 = t449 * t523;
t291 = qJD(5) * t479 + t448 * t416 + t443 * t506;
t557 = t439 * t451;
t315 = pkin(10) * t561 + (t524 * t557 + t432) * qJD(1);
t322 = pkin(2) * t427 + t357;
t383 = (-pkin(2) * t451 - t446 * t622 - pkin(1)) * t439;
t369 = qJD(1) * t383;
t228 = -t445 * t315 + t450 * (t322 * t440 + t369 * t438);
t474 = t440 * t552 + t551;
t463 = t474 * t439;
t324 = qJD(1) * t463 + t427 * t560;
t266 = pkin(3) * t324 - pkin(11) * t323;
t159 = -t228 * t444 + t449 * t266;
t136 = -pkin(12) * t323 * t449 + pkin(4) * t324 + t159;
t160 = t449 * t228 + t444 * t266;
t562 = t323 * t444;
t145 = -pkin(12) * t562 + t160;
t82 = t443 * t136 + t448 * t145;
t731 = t291 - t82;
t699 = pkin(4) * t540 + t696;
t360 = t462 * qJD(2);
t340 = qJD(1) * t360;
t388 = qJD(2) * t465;
t379 = qJD(1) * t388;
t536 = qJD(2) * t439;
t513 = qJD(1) * t536;
t503 = t446 * t513;
t420 = qJD(2) * t425;
t467 = qJD(2) * t477;
t339 = -qJD(1) * t467 + t420;
t517 = t440 * t535;
t508 = -t315 * t534 - t322 * t517 - t445 * t339 - t369 * t519;
t148 = -t340 * t555 + (-pkin(3) * t503 - t379 * t450) * t438 - t508;
t378 = t440 * t427 - t438 * t521 + qJD(3);
t280 = t324 * t449 + t378 * t444;
t694 = t440 * t548 - t553;
t458 = (t473 * qJD(2) + qJD(3) * t694) * t439;
t281 = qJD(1) * t458 + t427 * t518;
t478 = t438 * t503;
t190 = -qJD(4) * t280 - t281 * t444 + t449 * t478;
t102 = -pkin(4) * t190 + t148;
t279 = -t324 * t444 + t378 * t449;
t189 = qJD(4) * t279 + t281 * t449 + t444 * t478;
t457 = (qJD(2) * t475 + qJD(3) * t474) * t439;
t282 = qJD(1) * t457 + t427 * t519;
t269 = -t322 * t438 + t440 * t369;
t198 = -pkin(3) * t323 - pkin(11) * t324 + t269;
t229 = t450 * t315 + t322 * t556 + t369 * t560;
t203 = pkin(11) * t378 + t229;
t130 = t198 * t444 + t203 * t449;
t516 = t440 * t534;
t149 = -t315 * t535 + t322 * t516 + t450 * t339 + t340 * t556 + t369 * t518 + t379 * t560;
t147 = pkin(11) * t478 + t149;
t288 = -t340 * t438 + t440 * t379;
t165 = pkin(3) * t282 - pkin(11) * t281 + t288;
t58 = -qJD(4) * t130 - t147 * t444 + t449 * t165;
t38 = pkin(4) * t282 - pkin(12) * t189 + t58;
t57 = t449 * t147 + t444 * t165 + t198 * t532 - t203 * t533;
t47 = pkin(12) * t190 + t57;
t113 = pkin(12) * t279 + t130;
t549 = t448 * t113;
t129 = t449 * t198 - t203 * t444;
t112 = -pkin(12) * t280 + t129;
t319 = qJD(4) - t323;
t98 = pkin(4) * t319 + t112;
t54 = t443 * t98 + t549;
t11 = -qJD(5) * t54 + t38 * t448 - t443 * t47;
t638 = t282 / 0.2e1;
t484 = t279 * t443 + t448 * t280;
t96 = qJD(5) * t484 + t189 * t443 - t448 * t190;
t666 = -t96 / 0.2e1;
t510 = t448 * t279 - t280 * t443;
t95 = qJD(5) * t510 + t189 * t448 + t190 * t443;
t667 = t95 / 0.2e1;
t671 = Ifges(6,1) * t667 + Ifges(6,4) * t666 + Ifges(6,5) * t638;
t730 = t102 * mrSges(6,2) - t11 * mrSges(6,3) + 0.2e1 * t671;
t729 = pkin(13) * t732 - t710;
t728 = -pkin(13) * t324 + t731;
t316 = qJD(5) + t319;
t442 = sin(qJ(6));
t447 = cos(qJ(6));
t171 = t316 * t447 - t442 * t484;
t65 = qJD(6) * t171 + t282 * t442 + t447 * t95;
t172 = t316 * t442 + t447 * t484;
t66 = -qJD(6) * t172 + t282 * t447 - t442 * t95;
t21 = -mrSges(7,1) * t66 + mrSges(7,2) * t65;
t8 = -pkin(5) * t282 - t11;
t727 = -m(7) * t8 - t21;
t184 = pkin(4) * t562 + t229;
t726 = pkin(4) * t533 - pkin(5) * t542 + pkin(13) * t541 - t184;
t725 = t733 * pkin(5) - pkin(13) * t734 + t699;
t724 = t281 / 0.2e1;
t723 = -t282 / 0.2e1;
t722 = t378 / 0.2e1;
t700 = t443 * t272 + t448 * t289;
t709 = -qJD(5) * t700 + t736 * t443 + t448 * t737;
t372 = t422 * t443 - t423 * t448;
t292 = qJD(5) * t372 + t416 * t443 - t448 * t506;
t81 = t136 * t448 - t145 * t443;
t721 = -t292 - t81;
t10 = -t113 * t531 + t443 * t38 + t448 * t47 + t98 * t530;
t63 = Ifges(7,6) * t66;
t64 = Ifges(7,5) * t65;
t18 = Ifges(7,3) * t96 + t63 + t64;
t52 = pkin(13) * t316 + t54;
t202 = -pkin(3) * t378 - t228;
t155 = -pkin(4) * t279 + t202;
t90 = -pkin(5) * t510 - pkin(13) * t484 + t155;
t22 = -t442 * t52 + t447 * t90;
t24 = pkin(5) * t96 - pkin(13) * t95 + t102;
t7 = pkin(13) * t282 + t10;
t2 = qJD(6) * t22 + t24 * t442 + t447 * t7;
t23 = t442 * t90 + t447 * t52;
t3 = -qJD(6) * t23 + t24 * t447 - t442 * t7;
t502 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t665 = t96 / 0.2e1;
t688 = Ifges(6,6) * t723 - t95 * Ifges(6,4) / 0.2e1;
t720 = t502 + t102 * mrSges(6,1) - t10 * mrSges(6,3) + Ifges(6,2) * t665 + t18 / 0.2e1 + t688;
t569 = t378 * Ifges(4,6);
t572 = t324 * Ifges(4,4);
t575 = t323 * Ifges(4,2);
t249 = t569 + t572 + t575;
t564 = t113 * t443;
t53 = t448 * t98 - t564;
t719 = -t269 * mrSges(4,1) - t129 * mrSges(5,1) - t53 * mrSges(6,1) + t130 * mrSges(5,2) + t54 * mrSges(6,2) + t229 * mrSges(4,3) + t249 / 0.2e1 + t572 / 0.2e1;
t318 = Ifges(4,4) * t323;
t718 = t269 * mrSges(4,2) - t228 * mrSges(4,3) + t318 / 0.2e1;
t655 = t189 / 0.2e1;
t654 = t190 / 0.2e1;
t717 = Ifges(4,3) * t722;
t715 = -Ifges(3,6) * t427 / 0.2e1;
t201 = -pkin(13) * t559 + t700;
t391 = t428 + (-pkin(2) * t450 - pkin(3)) * t440;
t328 = -pkin(4) * t405 + t391;
t235 = -pkin(5) * t480 - pkin(13) * t321 + t328;
t137 = -t201 * t442 + t235 * t447;
t714 = qJD(6) * t137 + t725 * t442 - t447 * t729;
t138 = t201 * t447 + t235 * t442;
t713 = -qJD(6) * t138 + t442 * t729 + t725 * t447;
t712 = -t11 * mrSges(6,1) + t10 * mrSges(6,2);
t573 = t324 * Ifges(4,1);
t711 = pkin(5) * t732 - t709;
t437 = -pkin(4) * t449 - pkin(3);
t338 = pkin(5) * t413 - pkin(13) * t414 + t437;
t284 = t338 * t447 - t372 * t442;
t708 = qJD(6) * t284 + t442 * t726 + t447 * t728;
t285 = t338 * t442 + t372 * t447;
t707 = -qJD(6) * t285 - t442 * t728 + t447 * t726;
t706 = pkin(5) * t324 - t721;
t597 = Ifges(6,6) * t316;
t599 = Ifges(6,2) * t510;
t605 = Ifges(6,4) * t484;
t123 = t597 + t599 + t605;
t205 = qJD(6) - t510;
t595 = Ifges(7,3) * t205;
t596 = Ifges(7,6) * t171;
t600 = Ifges(7,5) * t172;
t87 = t595 + t596 + t600;
t705 = t87 - t123;
t356 = pkin(2) * t441 + t433 - t477;
t290 = -t356 * t438 + t440 * t383;
t348 = -t439 * t694 - t441 * t559;
t349 = t441 * t560 + t463;
t226 = pkin(3) * t348 - pkin(11) * t349 + t290;
t538 = pkin(9) * t557 + t432;
t342 = (t438 * t441 + t440 * t557) * pkin(10) + t538;
t252 = t450 * t342 + t356 * t556 + t383 * t560;
t404 = -t438 * t557 + t441 * t440;
t234 = pkin(11) * t404 + t252;
t143 = t449 * t226 - t234 * t444;
t299 = t349 * t449 + t404 * t444;
t121 = pkin(4) * t348 - pkin(12) * t299 + t143;
t144 = t444 * t226 + t449 * t234;
t298 = -t349 * t444 + t404 * t449;
t128 = pkin(12) * t298 + t144;
t704 = t443 * t121 + t448 * t128;
t577 = t316 * Ifges(6,3);
t589 = t484 * Ifges(6,5);
t590 = t510 * Ifges(6,6);
t122 = t577 + t589 + t590;
t576 = t319 * Ifges(5,3);
t583 = t280 * Ifges(5,5);
t585 = t279 * Ifges(5,6);
t175 = t576 + t583 + t585;
t703 = t175 + t122;
t697 = -(t358 * t440 + t387 * t438) * t450 + t735;
t204 = Ifges(6,4) * t510;
t601 = Ifges(6,5) * t316;
t609 = Ifges(6,1) * t484;
t124 = t204 + t601 + t609;
t493 = Ifges(7,5) * t447 - Ifges(7,6) * t442;
t469 = t205 * t493;
t603 = Ifges(7,4) * t442;
t497 = Ifges(7,1) * t447 - t603;
t470 = t172 * t497;
t602 = Ifges(7,4) * t447;
t495 = -Ifges(7,2) * t442 + t602;
t471 = t171 * t495;
t498 = mrSges(7,1) * t442 + mrSges(7,2) * t447;
t51 = -pkin(5) * t316 - t53;
t472 = t51 * t498;
t626 = -t447 / 0.2e1;
t627 = t442 / 0.2e1;
t685 = t155 * mrSges(6,2) - t53 * mrSges(6,3);
t604 = Ifges(7,4) * t172;
t88 = Ifges(7,2) * t171 + Ifges(7,6) * t205 + t604;
t170 = Ifges(7,4) * t171;
t89 = Ifges(7,1) * t172 + Ifges(7,5) * t205 + t170;
t695 = -t124 / 0.2e1 - t601 / 0.2e1 - t471 / 0.2e1 - t470 / 0.2e1 - t469 / 0.2e1 + t88 * t627 + t89 * t626 - t472 - t685 - t204 / 0.2e1;
t693 = -t444 * t58 + t449 * t57;
t692 = -t22 * t442 + t23 * t447;
t691 = -t58 * mrSges(5,1) + t57 * mrSges(5,2);
t570 = t378 * Ifges(4,5);
t250 = t318 + t570 + t573;
t643 = t250 / 0.2e1;
t689 = t643 + t570 / 0.2e1 + t718;
t687 = Ifges(4,1) * t724 + Ifges(4,4) * t723;
t584 = t280 * Ifges(5,4);
t176 = t279 * Ifges(5,2) + t319 * Ifges(5,6) + t584;
t273 = Ifges(5,4) * t279;
t177 = t280 * Ifges(5,1) + t319 * Ifges(5,5) + t273;
t487 = t129 * t449 + t130 * t444;
t606 = Ifges(5,4) * t449;
t607 = Ifges(5,4) * t444;
t624 = t449 / 0.2e1;
t633 = t319 / 0.2e1;
t639 = t280 / 0.2e1;
t641 = t279 / 0.2e1;
t686 = -t487 * mrSges(5,3) + t202 * (mrSges(5,1) * t444 + mrSges(5,2) * t449) + (-Ifges(5,2) * t444 + t606) * t641 + (Ifges(5,1) * t449 - t607) * t639 + (Ifges(5,5) * t449 - Ifges(5,6) * t444) * t633 - t444 * t176 / 0.2e1 + t177 * t624;
t251 = -t445 * t342 + t450 * (t356 * t440 + t383 * t438);
t293 = t441 * t518 + t458;
t520 = t446 * t536;
t504 = t438 * t520;
t218 = qJD(4) * t298 + t293 * t449 + t444 * t504;
t294 = t441 * t519 + t457;
t426 = qJD(2) * t433;
t359 = t426 - t467;
t168 = -t342 * t535 + t356 * t516 + t450 * t359 + t360 * t556 + t383 * t518 + t388 * t560;
t157 = pkin(11) * t504 + t168;
t296 = -t360 * t438 + t440 * t388;
t180 = pkin(3) * t294 - pkin(11) * t293 + t296;
t72 = -qJD(4) * t144 - t157 * t444 + t449 * t180;
t50 = pkin(4) * t294 - pkin(12) * t218 + t72;
t219 = -qJD(4) * t299 - t293 * t444 + t449 * t504;
t71 = t449 * t157 + t444 * t180 + t226 * t532 - t234 * t533;
t56 = pkin(12) * t219 + t71;
t15 = -qJD(5) * t704 - t443 * t56 + t448 * t50;
t135 = pkin(5) * t484 - pkin(13) * t510;
t681 = t155 * mrSges(6,1) + t22 * mrSges(7,1) - t23 * mrSges(7,2) - t54 * mrSges(6,3);
t679 = t123 / 0.2e1 - t87 / 0.2e1 - t681;
t669 = t66 / 0.2e1;
t670 = t65 / 0.2e1;
t678 = -Ifges(6,4) * t667 + Ifges(7,5) * t670 - Ifges(6,2) * t666 - Ifges(6,6) * t638 + Ifges(7,6) * t669 + Ifges(7,3) * t665 + t720;
t676 = Ifges(6,2) / 0.2e1;
t19 = t65 * Ifges(7,4) + t66 * Ifges(7,2) + t96 * Ifges(7,6);
t674 = t19 / 0.2e1;
t20 = t65 * Ifges(7,1) + t66 * Ifges(7,4) + t96 * Ifges(7,5);
t673 = t20 / 0.2e1;
t668 = t89 / 0.2e1;
t663 = pkin(1) * mrSges(3,1);
t662 = pkin(1) * mrSges(3,2);
t108 = t189 * Ifges(5,4) + t190 * Ifges(5,2) + t282 * Ifges(5,6);
t661 = t108 / 0.2e1;
t660 = Ifges(5,1) * t655 + Ifges(5,4) * t654 + Ifges(5,5) * t638;
t659 = -t171 / 0.2e1;
t658 = t171 / 0.2e1;
t657 = -t172 / 0.2e1;
t656 = t172 / 0.2e1;
t276 = Ifges(4,6) * t282;
t277 = Ifges(4,5) * t281;
t195 = Ifges(4,3) * t478 - t276 + t277;
t653 = t195 / 0.2e1;
t652 = Ifges(4,5) * t478 / 0.2e1 + t687;
t650 = -t205 / 0.2e1;
t649 = t205 / 0.2e1;
t648 = -t510 / 0.2e1;
t647 = t510 / 0.2e1;
t646 = -t484 / 0.2e1;
t645 = t484 / 0.2e1;
t642 = -t279 / 0.2e1;
t640 = -t280 / 0.2e1;
t636 = -t316 / 0.2e1;
t635 = t316 / 0.2e1;
t634 = -t319 / 0.2e1;
t630 = -t378 / 0.2e1;
t629 = t441 / 0.2e1;
t628 = -t442 / 0.2e1;
t625 = t447 / 0.2e1;
t94 = Ifges(6,5) * t95;
t93 = Ifges(6,6) * t96;
t620 = t2 * t447;
t619 = t3 * t442;
t611 = qJD(2) / 0.2e1;
t608 = Ifges(3,4) * t446;
t188 = Ifges(5,5) * t189;
t187 = Ifges(5,6) * t190;
t592 = t149 * mrSges(4,2);
t150 = (t340 * t440 + t379 * t438) * t450 + t508;
t591 = t150 * mrSges(4,1);
t587 = t22 * t447;
t581 = t281 * Ifges(4,4);
t558 = t439 * t446;
t554 = t442 * t353;
t550 = t447 * t353;
t110 = -mrSges(7,1) * t171 + mrSges(7,2) * t172;
t174 = mrSges(6,1) * t316 - mrSges(6,3) * t484;
t547 = t110 - t174;
t300 = -t321 * t442 - t447 * t559;
t182 = qJD(6) * t300 + t246 * t447 + t442 * t519;
t232 = t258 * t447 + t375 * t442;
t546 = t182 - t232;
t476 = -t321 * t447 + t442 * t559;
t183 = qJD(6) * t476 - t246 * t442 + t447 * t519;
t231 = -t258 * t442 + t375 * t447;
t545 = t183 - t231;
t529 = qJD(6) * t442;
t528 = qJD(6) * t447;
t34 = Ifges(6,3) * t282 - t93 + t94;
t107 = Ifges(5,3) * t282 + t187 + t188;
t216 = t260 * t442 + t324 * t447;
t512 = t216 - t554;
t217 = -t260 * t447 + t324 * t442;
t511 = -t217 - t550;
t399 = t538 * qJD(1);
t509 = t715 - (Ifges(3,2) * t451 + t608) * t537 / 0.2e1 - t399 * mrSges(3,3);
t507 = -t342 * t534 - t356 * t517 - t445 * t359 - t383 * t519;
t500 = t591 - t592;
t499 = mrSges(7,1) * t447 - mrSges(7,2) * t442;
t496 = Ifges(7,1) * t442 + t602;
t494 = Ifges(7,2) * t447 + t603;
t492 = Ifges(7,5) * t442 + Ifges(7,6) * t447;
t491 = t23 * t442 + t587;
t233 = -pkin(3) * t404 - t251;
t181 = -pkin(4) * t298 + t233;
t237 = t298 * t443 + t299 * t448;
t483 = t448 * t298 - t299 * t443;
t101 = -pkin(5) * t483 - pkin(13) * t237 + t181;
t68 = pkin(13) * t348 + t704;
t30 = t101 * t447 - t442 * t68;
t31 = t101 * t442 + t447 * t68;
t69 = t121 * t448 - t128 * t443;
t193 = t237 * t447 + t348 * t442;
t192 = -t237 * t442 + t348 * t447;
t210 = t272 * t448 - t289 * t443;
t14 = t121 * t530 - t128 * t531 + t443 * t50 + t448 * t56;
t397 = -pkin(9) * t522 + t425;
t424 = Ifges(3,4) * t521;
t466 = Ifges(3,1) * t522 / 0.2e1 + t424 / 0.2e1 + t427 * Ifges(3,5) - t397 * mrSges(3,3);
t403 = t538 * qJD(2);
t460 = t228 * mrSges(4,1) - t229 * mrSges(4,2) + t324 * Ifges(4,5) + t323 * Ifges(4,6) + t717;
t158 = -t360 * t555 + (-pkin(3) * t520 - t388 * t450) * t438 - t507;
t459 = mrSges(7,3) * t620 + qJD(6) * t472 + t19 * t625 + t20 * t627 + t494 * t669 + t496 * t670 - t8 * t499 + t34 - t88 * t529 / 0.2e1 + t528 * t668 + t492 * t665 - t712 + (t471 + t470 + t469) * qJD(6) / 0.2e1;
t115 = -pkin(4) * t219 + t158;
t119 = -mrSges(7,2) * t205 + mrSges(7,3) * t171;
t120 = mrSges(7,1) * t205 - mrSges(7,3) * t172;
t25 = mrSges(7,1) * t96 - mrSges(7,3) * t65;
t26 = -mrSges(7,2) * t96 + mrSges(7,3) * t66;
t456 = -t120 * t528 - t119 * t529 - t442 * t25 + m(7) * (-t22 * t528 - t23 * t529 - t619 + t620) + t447 * t26;
t455 = -t595 / 0.2e1 - t600 / 0.2e1 - t596 / 0.2e1 + t597 / 0.2e1 + t605 / 0.2e1 + t679;
t453 = -t577 / 0.2e1 - t583 / 0.2e1 - t585 / 0.2e1 - t576 / 0.2e1 - t590 / 0.2e1 - t589 / 0.2e1 + t569 / 0.2e1 - t122 / 0.2e1 - t175 / 0.2e1 + t719;
t419 = Ifges(3,5) * t451 * t513;
t409 = -pkin(9) * t558 + t433;
t402 = -pkin(9) * t520 + t426;
t396 = -t427 * mrSges(3,2) + mrSges(3,3) * t521;
t395 = mrSges(3,1) * t427 - mrSges(3,3) * t522;
t390 = qJD(1) * t403;
t389 = -pkin(9) * t503 + t420;
t287 = mrSges(4,1) * t378 - mrSges(4,3) * t324;
t286 = -mrSges(4,2) * t378 + mrSges(4,3) * t323;
t265 = -mrSges(4,1) * t323 + mrSges(4,2) * t324;
t262 = -mrSges(4,2) * t478 - mrSges(4,3) * t282;
t261 = mrSges(4,1) * t478 - mrSges(4,3) * t281;
t241 = mrSges(5,1) * t319 - mrSges(5,3) * t280;
t240 = -mrSges(5,2) * t319 + mrSges(5,3) * t279;
t213 = mrSges(4,1) * t282 + mrSges(4,2) * t281;
t212 = -mrSges(5,1) * t279 + mrSges(5,2) * t280;
t200 = pkin(5) * t559 - t210;
t196 = -t282 * Ifges(4,2) + Ifges(4,6) * t478 + t581;
t173 = -mrSges(6,2) * t316 + mrSges(6,3) * t510;
t169 = (t360 * t440 + t388 * t438) * t450 + t507;
t152 = -mrSges(5,2) * t282 + mrSges(5,3) * t190;
t151 = mrSges(5,1) * t282 - mrSges(5,3) * t189;
t134 = -mrSges(6,1) * t510 + mrSges(6,2) * t484;
t126 = -mrSges(5,1) * t190 + mrSges(5,2) * t189;
t114 = pkin(4) * t280 + t135;
t106 = qJD(5) * t237 + t218 * t443 - t448 * t219;
t105 = qJD(5) * t483 + t218 * t448 + t219 * t443;
t86 = -mrSges(6,2) * t282 - mrSges(6,3) * t96;
t85 = mrSges(6,1) * t282 - mrSges(6,3) * t95;
t74 = -qJD(6) * t193 - t105 * t442 + t294 * t447;
t73 = qJD(6) * t192 + t105 * t447 + t294 * t442;
t67 = -pkin(5) * t348 - t69;
t60 = t112 * t448 - t564;
t59 = t112 * t443 + t549;
t41 = mrSges(6,1) * t96 + mrSges(6,2) * t95;
t33 = t135 * t442 + t447 * t53;
t32 = t135 * t447 - t442 * t53;
t29 = t114 * t442 + t447 * t60;
t28 = t114 * t447 - t442 * t60;
t27 = pkin(5) * t106 - pkin(13) * t105 + t115;
t13 = -pkin(5) * t294 - t15;
t12 = pkin(13) * t294 + t14;
t5 = -qJD(6) * t31 - t12 * t442 + t27 * t447;
t4 = qJD(6) * t30 + t12 * t447 + t27 * t442;
t1 = [(Ifges(7,4) * t73 + Ifges(7,2) * t74) * t658 + (Ifges(7,4) * t193 + Ifges(7,2) * t192) * t669 + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t649 + (Ifges(7,5) * t193 + Ifges(7,6) * t192) * t665 + (-Ifges(6,4) * t645 + Ifges(7,5) * t656 - Ifges(6,2) * t647 - Ifges(6,6) * t635 + Ifges(7,6) * t658 + Ifges(7,3) * t649 - t679) * t106 + ((Ifges(3,5) * t629 - t409 * mrSges(3,3) + (-0.2e1 * t662 + 0.3e1 / 0.2e1 * Ifges(3,4) * t451) * t439) * t451 + (-Ifges(3,6) * t441 - t538 * mrSges(3,3) + t438 * (Ifges(4,5) * t349 - Ifges(4,6) * t348 + Ifges(4,3) * t404) / 0.2e1 + (-0.2e1 * t663 - 0.3e1 / 0.2e1 * t608 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t451) * t439) * t446) * t513 + (Ifges(5,5) * t299 + Ifges(6,5) * t237 + Ifges(5,6) * t298 + (Ifges(5,3) + Ifges(6,3)) * t348) * t638 + (t192 * t2 - t193 * t3 - t22 * t73 + t23 * t74) * mrSges(7,3) + (-t10 * t348 + t102 * t237 + t105 * t155 - t294 * t54) * mrSges(6,2) + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t656 + (Ifges(7,1) * t193 + Ifges(7,4) * t192) * t670 + (Ifges(6,4) * t105 + Ifges(6,6) * t294) * t647 + (Ifges(6,4) * t237 + Ifges(6,6) * t348) * t666 + (Ifges(6,5) * t105 + Ifges(6,3) * t294) * t635 + m(6) * (t10 * t704 + t102 * t181 + t11 * t69 + t115 * t155 + t14 * t54 + t15 * t53) + t704 * t86 - t678 * t483 + t349 * t652 + t404 * t653 + (Ifges(5,4) * t299 + Ifges(5,2) * t298 + Ifges(5,6) * t348) * t654 + (Ifges(6,1) * t105 + Ifges(6,5) * t294) * t645 + (Ifges(6,1) * t237 + Ifges(6,5) * t348) * t667 + (-t149 * t348 - t150 * t349 - t228 * t293 - t229 * t294) * mrSges(4,3) + t389 * (-t441 * mrSges(3,2) + mrSges(3,3) * t557) - t390 * (mrSges(3,1) * t441 - mrSges(3,3) * t558) + (Ifges(4,5) * t293 - Ifges(4,6) * t294) * t722 + (Ifges(4,4) * t349 - Ifges(4,2) * t348 + Ifges(4,6) * t404) * t723 + (Ifges(4,1) * t349 - Ifges(4,4) * t348 + Ifges(4,5) * t404) * t724 + m(4) * (t149 * t252 + t150 * t251 + t168 * t229 + t169 * t228 + t269 * t296 + t288 * t290) + m(5) * (t129 * t72 + t130 * t71 + t143 * t58 + t144 * t57 + t148 * t233 + t158 * t202) + m(7) * (t13 * t51 + t2 * t31 + t22 * t5 + t23 * t4 + t3 * t30 + t67 * t8) + (Ifges(5,1) * t218 + Ifges(5,4) * t219 + Ifges(5,5) * t294) * t639 + (Ifges(5,4) * t218 + Ifges(5,2) * t219 + Ifges(5,6) * t294) * t641 + t293 * t643 + m(3) * (t389 * t538 - t390 * t409 - t397 * t403 + t399 * t402) + t323 * (Ifges(4,4) * t293 - Ifges(4,2) * t294) / 0.2e1 - t404 * t592 + t402 * t396 - t403 * t395 + t419 * t629 + (Ifges(5,5) * t218 + Ifges(5,6) * t219 + Ifges(5,3) * t294) * t633 + t404 * t591 + t57 * (-mrSges(5,2) * t348 + mrSges(5,3) * t298) + t58 * (mrSges(5,1) * t348 - mrSges(5,3) * t299) + t11 * (mrSges(6,1) * t348 - mrSges(6,3) * t237) - t348 * t196 / 0.2e1 + t288 * (mrSges(4,1) * t348 + mrSges(4,2) * t349) + t148 * (-mrSges(5,1) * t298 + mrSges(5,2) * t299) + t269 * (mrSges(4,1) * t294 + mrSges(4,2) * t293) - t294 * t249 / 0.2e1 + t296 * t265 + t129 * (mrSges(5,1) * t294 - mrSges(5,3) * t218) + t130 * (-mrSges(5,2) * t294 + mrSges(5,3) * t219) + t53 * (mrSges(6,1) * t294 - mrSges(6,3) * t105) + t290 * t213 + t168 * t286 + t169 * t287 + t252 * t262 + t251 * t261 + t72 * t241 + t71 * t240 + t233 * t126 + t218 * t177 / 0.2e1 + t202 * (-mrSges(5,1) * t219 + mrSges(5,2) * t218) + t219 * t176 / 0.2e1 + t158 * t212 + t8 * (-mrSges(7,1) * t192 + mrSges(7,2) * t193) + t15 * t174 + t181 * t41 + t14 * t173 + t143 * t151 + t144 * t152 + t115 * t134 + t105 * t124 / 0.2e1 + t5 * t120 + t4 * t119 + t13 * t110 + (Ifges(5,1) * t299 + Ifges(5,4) * t298 + Ifges(5,5) * t348) * t655 + t299 * t660 + t298 * t661 + t703 * t294 / 0.2e1 + t73 * t668 + t237 * t671 + t193 * t673 + t192 * t674 + (t107 + t34) * t348 / 0.2e1 + (t466 * t451 + (t715 + (t717 + t460) * t438 + t509) * t446) * t536 + t30 * t25 + t31 * t26 + t324 * (Ifges(4,1) * t293 - Ifges(4,4) * t294) / 0.2e1 + t67 * t21 + t51 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t69 * t85 + t74 * t88 / 0.2e1; (-t703 / 0.2e1 + Ifges(6,6) * t648 + t575 / 0.2e1 + Ifges(5,3) * t634 + Ifges(6,3) * t636 + Ifges(5,5) * t640 + Ifges(5,6) * t642 + Ifges(6,5) * t646 - Ifges(4,6) * t630 + t719) * t375 + (Ifges(6,4) * t258 - Ifges(6,2) * t257) * t648 + (Ifges(5,4) * t314 + Ifges(5,2) * t313) * t642 + (Ifges(6,5) * t258 - Ifges(6,6) * t257) * t636 + (Ifges(5,5) * t314 + Ifges(5,6) * t313) * t634 + t700 * t86 + (t10 * t700 + t102 * t328 + t11 * t210 + t155 * t699 + t53 * t709 + t54 * t710) * m(6) + (Ifges(6,4) * t246 - Ifges(6,2) * t247) * t647 + (Ifges(7,5) * t182 + Ifges(7,6) * t183 + Ifges(7,3) * t247) * t649 + (Ifges(7,5) * t232 + Ifges(7,6) * t231 + Ifges(7,3) * t257) * t650 + (Ifges(5,4) * t406 + Ifges(5,2) * t405) * t654 + (t361 / 0.2e1 - t314 / 0.2e1) * t177 + (-mrSges(7,1) * t545 + mrSges(7,2) * t546) * t51 + (-t231 / 0.2e1 + t183 / 0.2e1) * t88 + (-t258 / 0.2e1 + t246 / 0.2e1) * t124 + (-t232 / 0.2e1 + t182 / 0.2e1) * t89 + (t653 - t276 / 0.2e1 + t277 / 0.2e1 + t500) * t440 + (Ifges(6,1) * t258 - Ifges(6,4) * t257) * t646 + t730 * t321 + t681 * t733 + t685 * t734 + (Ifges(5,1) * t314 + Ifges(5,4) * t313) * t640 + (Ifges(6,5) * t246 - Ifges(6,6) * t247) * t635 + (Ifges(5,5) * t406 + Ifges(5,6) * t405) * t638 + (Ifges(5,1) * t361 + Ifges(5,4) * t362) * t639 + (Ifges(5,4) * t361 + Ifges(5,2) * t362) * t641 + (Ifges(6,1) * t246 - Ifges(6,4) * t247) * t645 + (-Ifges(7,4) * t476 + Ifges(7,2) * t300) * t669 + (-Ifges(7,1) * t476 + Ifges(7,4) * t300) * t670 + (t2 * t300 - t22 * t546 + t23 * t545 + t3 * t476) * mrSges(7,3) + (-Ifges(7,5) * t476 + Ifges(7,6) * t300) * t665 + t8 * (-mrSges(7,1) * t300 - mrSges(7,2) * t476) - t476 * t673 - t678 * t480 + (mrSges(5,1) * t540 - mrSges(5,2) * t539) * t202 + (t129 * t539 - t130 * t540 + t405 * t57 - t406 * t58) * mrSges(5,3) + t148 * (-mrSges(5,1) * t405 + mrSges(5,2) * t406) + t408 * t261 + t410 * t262 - t389 * mrSges(3,2) - t390 * mrSges(3,1) + t391 * t126 - t397 * t396 + t399 * t395 + (Ifges(5,5) * t361 + Ifges(5,6) * t362) * t633 + (-t250 / 0.2e1 + Ifges(4,5) * t630 - t573 / 0.2e1 - t718) * t376 + (t362 / 0.2e1 - t313 / 0.2e1) * t176 + t328 * t41 + t303 * t151 + t304 * t152 - t295 * t265 + t210 * t85 + t200 * t21 + t138 * t26 + t137 * t25 + (Ifges(5,1) * t406 + Ifges(5,4) * t405) * t655 + (Ifges(7,1) * t182 + Ifges(7,4) * t183 + Ifges(7,5) * t247) * t656 + (Ifges(7,1) * t232 + Ifges(7,4) * t231 + Ifges(7,5) * t257) * t657 + (Ifges(7,4) * t182 + Ifges(7,2) * t183 + Ifges(7,6) * t247) * t658 + (Ifges(7,4) * t232 + Ifges(7,2) * t231 + Ifges(7,6) * t257) * t659 + t406 * t660 + t405 * t661 + t419 + t696 * t212 + t697 * t287 + t698 * t286 + (-pkin(2) * t288 * t438 + t149 * t410 + t150 * t408 + t228 * t697 + t229 * t698 - t269 * t295) * m(4) + t699 * t134 + t701 * t241 + t702 * t240 + (t129 * t701 + t130 * t702 + t148 * t391 + t202 * t696 + t303 * t58 + t304 * t57) * m(5) + t300 * t674 - t705 * (t257 / 0.2e1 - t247 / 0.2e1) + t709 * t174 + t710 * t173 + t711 * t110 + (-pkin(2) * t213 + (t288 * mrSges(4,2) - t150 * mrSges(4,3) + t652 + t687) * t445 + (t149 * mrSges(4,3) + t93 / 0.2e1 - t94 / 0.2e1 - t188 / 0.2e1 - t187 / 0.2e1 + t581 / 0.2e1 - t288 * mrSges(4,1) + t196 / 0.2e1 - t34 / 0.2e1 - t107 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1) * t282 + t691 + t712) * t450 + ((t573 / 0.2e1 + t689) * t450 + (-t575 / 0.2e1 - t453) * t445) * qJD(3)) * t438 + t713 * t120 + t714 * t119 + (t137 * t3 + t138 * t2 + t200 * t8 + t22 * t713 + t23 * t714 + t51 * t711) * m(7) + ((-t424 / 0.2e1 + t537 * t662 - t466) * t451 + ((t663 + t608 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t451) * t537 + (-qJD(2) + t427 / 0.2e1) * Ifges(3,6) + ((Ifges(4,5) * t445 + Ifges(4,6) * t450) * t438 * t611 + (t440 * t611 + t630) * Ifges(4,3) - t460) * t438 - t509) * t446) * t537; (t64 / 0.2e1 + t63 / 0.2e1 + (t676 + Ifges(7,3) / 0.2e1) * t96 + t688 + t720) * t413 + (-Ifges(6,4) * t260 - Ifges(6,2) * t259) * t648 + (-Ifges(6,5) * t260 - Ifges(6,6) * t259) * t636 + (-Ifges(6,1) * t260 - Ifges(6,4) * t259) * t646 + (-pkin(3) * t148 - t129 * t159 - t130 * t160 - t202 * t229) * m(5) + (Ifges(7,5) * t217 + Ifges(7,6) * t216 + Ifges(7,3) * t259) * t650 + (Ifges(5,2) * t449 + t607) * t654 + (t53 * t541 + t54 * t542) * mrSges(6,3) + (mrSges(7,2) * t542 - mrSges(7,3) * t512) * t23 + (-mrSges(6,1) * t542 - mrSges(6,2) * t541) * t155 + (-mrSges(7,1) * t542 - mrSges(7,3) * t511) * t22 - m(6) * (t155 * t184 + t53 * t81 + t54 * t82) + t195 + (t497 * t670 + t495 * t669 + t493 * t665 + t8 * t498 + t19 * t628 + t20 * t625 + (-t2 * t442 - t3 * t447) * mrSges(7,3) + (-mrSges(7,3) * t692 + t492 * t650 + t494 * t659 + t496 * t657 + t51 * t499 + t88 * t626 + t89 * t628) * qJD(6) + t730) * t414 + t500 + t731 * t173 + (Ifges(5,5) * t444 + Ifges(5,6) * t449) * t638 + (mrSges(7,1) * t512 + mrSges(7,2) * t511) * t51 + m(6) * (t10 * t372 + t102 * t437 + t11 * t479 + t291 * t54 - t292 * t53) - (t21 - t85) * t479 + (t2 * t285 + t22 * t707 + t23 * t708 + t284 * t3 - t479 * t8 + t51 * t706) * m(7) + (-t353 / 0.2e1 + t260 / 0.2e1) * t124 + (-Ifges(6,5) * t353 - Ifges(6,6) * t354) * t635 + (-Ifges(6,1) * t353 - Ifges(6,4) * t354) * t645 + (-Ifges(6,4) * t353 - Ifges(6,2) * t354) * t647 + (-t217 / 0.2e1 - t550 / 0.2e1) * t89 + (-t216 / 0.2e1 + t554 / 0.2e1) * t88 + (-Ifges(7,5) * t550 + Ifges(7,6) * t554 + Ifges(7,3) * t354) * t649 + (-Ifges(7,1) * t550 + Ifges(7,4) * t554 + Ifges(7,5) * t354) * t656 + (-Ifges(7,4) * t550 + Ifges(7,2) * t554 + Ifges(7,6) * t354) * t658 + t437 * t41 + t372 * t86 + t108 * t624 + t284 * t25 + t285 * t26 - t228 * t286 - t160 * t240 - t159 * t241 - t184 * t134 - pkin(3) * t126 + (t287 - t212) * t229 + (Ifges(5,1) * t444 + t606) * t655 + (Ifges(7,1) * t217 + Ifges(7,4) * t216 + Ifges(7,5) * t259) * t657 + (Ifges(7,4) * t217 + Ifges(7,2) * t216 + Ifges(7,6) * t259) * t659 + t444 * t660 + t148 * (-mrSges(5,1) * t449 + mrSges(5,2) * t444) + t721 * t174 + ((m(6) * t155 + t134) * t444 * pkin(4) + t686) * qJD(4) + ((Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t324 - t686 - t689) * t323 + (m(5) * t693 + (-m(5) * t487 - t444 * t240 - t449 * t241) * qJD(4) - t151 * t444 + t152 * t449) * pkin(11) + t693 * mrSges(5,3) + t705 * (t354 / 0.2e1 - t259 / 0.2e1) + t706 * t110 + t707 * t120 + t708 * t119 + t453 * t324; -t547 * t59 - m(7) * (t22 * t28 + t23 * t29 + t51 * t59) - t691 + (Ifges(5,5) * t279 - Ifges(5,6) * t280) * t634 + t176 * t639 + (Ifges(5,1) * t279 - t584) * t640 + (t129 * t279 + t130 * t280) * mrSges(5,3) + (-t205 * t587 + (-t205 * t23 - t3) * t442) * mrSges(7,3) + t107 - t202 * (mrSges(5,1) * t280 + mrSges(5,2) * t279) - t129 * t240 + t130 * t241 - t60 * t173 - t29 * t119 - t28 * t120 + t459 + (-t609 / 0.2e1 + t695) * t510 - m(6) * (-t53 * t59 + t54 * t60) + (t599 / 0.2e1 + t455) * t484 + (-t280 * t134 + t443 * t86 + t448 * t85 + ((m(7) * t51 + t547) * t443 + (m(7) * t692 + t119 * t447 - t120 * t442 + t173) * t448) * qJD(5) + (0.2e1 * t155 * t640 + t10 * t443 + t11 * t448 + (-t443 * t53 + t448 * t54) * qJD(5)) * m(6)) * pkin(4) + t456 * (pkin(4) * t443 + pkin(13)) + (-Ifges(5,2) * t280 + t177 + t273) * t642 - t727 * (-pkin(4) * t448 - pkin(5)); -t547 * t54 - m(7) * (t22 * t32 + t23 * t33 + t51 * t54) - t53 * t173 - t33 * t119 - t32 * t120 + t459 + (-qJD(6) * t491 - t619) * mrSges(7,3) + ((-Ifges(6,1) / 0.2e1 + t676) * t484 + t491 * mrSges(7,3) + t695) * t510 + t455 * t484 + t456 * pkin(13) + t727 * pkin(5); t88 * t656 - t22 * t119 + (Ifges(7,1) * t171 - t604) * t657 + (Ifges(7,5) * t171 - Ifges(7,6) * t172) * t650 + t23 * t120 - t51 * (mrSges(7,1) * t172 + mrSges(7,2) * t171) + (t171 * t22 + t172 * t23) * mrSges(7,3) + t502 + t18 + (-Ifges(7,2) * t172 + t170 + t89) * t659;];
tauc  = t1(:);
