% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR15_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:35
% EndTime: 2019-03-09 20:26:05
% DurationCPUTime: 47.48s
% Computational Cost: add. (27840->978), mult. (83113->1294), div. (0->0), fcn. (67934->12), ass. (0->459)
t366 = sin(qJ(2));
t509 = cos(pkin(6));
t482 = pkin(1) * t509;
t357 = t366 * t482;
t362 = sin(pkin(6));
t508 = cos(pkin(7));
t539 = cos(qJ(2));
t448 = t508 * t539;
t373 = -t357 + (-pkin(9) * t539 - pkin(10) * t448) * t362;
t268 = t373 * qJD(1);
t361 = sin(pkin(7));
t537 = pkin(10) * t361;
t392 = t362 * (pkin(2) * t366 - t537 * t539);
t300 = qJD(1) * t392;
t208 = -t268 * t361 + t508 * t300;
t480 = t362 * t539;
t538 = cos(qJ(3));
t421 = t538 * t480;
t403 = qJD(1) * t421;
t365 = sin(qJ(3));
t467 = t365 * t508;
t502 = t362 * t366;
t415 = t467 * t502;
t293 = -qJD(1) * t415 + t403;
t493 = qJD(3) * t365;
t474 = t361 * t493;
t679 = pkin(3) * t474 + qJ(4) * t293 - t208;
t411 = t362 * t448;
t395 = t538 * t411;
t417 = qJD(1) * t509 + qJD(2);
t399 = t361 * t417;
t495 = qJD(1) * t362;
t476 = t366 * t495;
t460 = t365 * t476;
t228 = -qJD(1) * t395 - t538 * t399 + t460;
t458 = qJD(1) * t480;
t294 = t361 * t458 - t417 * t508 - qJD(3);
t364 = sin(qJ(5));
t368 = cos(qJ(5));
t192 = t228 * t368 + t294 * t364;
t561 = -t192 / 0.2e1;
t672 = Ifges(6,2) * t561;
t356 = pkin(2) * t467;
t358 = t539 * t482;
t350 = qJD(1) * t358;
t397 = (-pkin(10) * t508 - pkin(9)) * t502;
t267 = qJD(1) * t397 + t350;
t447 = t508 * t538;
t396 = -t365 * t267 + t268 * t447;
t573 = pkin(3) + pkin(11);
t464 = t573 * t502;
t478 = t538 * t300;
t481 = t361 * t538;
t635 = pkin(4) + pkin(10);
t678 = -t293 * pkin(4) - (-qJD(1) * t464 - t478) * t361 + t396 + (t481 * t635 + t356) * qJD(3);
t382 = t365 * t539 + t366 * t447;
t292 = t382 * t495;
t492 = qJD(4) * t365;
t677 = t292 * t573 - (-t492 + (pkin(11) * t365 - qJ(4) * t538) * qJD(3)) * t361 - t679;
t383 = t365 * t448 + t366 * t538;
t378 = t383 * t362;
t394 = t365 * t399;
t229 = qJD(1) * t378 + t394;
t226 = qJD(5) + t229;
t552 = -t226 / 0.2e1;
t193 = t228 * t364 - t294 * t368;
t559 = -t193 / 0.2e1;
t676 = -Ifges(6,4) * t559 - Ifges(6,6) * t552;
t186 = qJD(6) - t192;
t562 = t186 / 0.2e1;
t363 = sin(qJ(6));
t367 = cos(qJ(6));
t147 = t193 * t367 + t226 * t363;
t565 = t147 / 0.2e1;
t146 = -t193 * t363 + t226 * t367;
t567 = t146 / 0.2e1;
t675 = Ifges(7,5) * t565 + Ifges(7,6) * t567 + Ifges(7,3) * t562;
t503 = t361 * t365;
t162 = t538 * t267 + t268 * t467 + t300 * t503;
t461 = t361 * t476;
t674 = qJ(4) * t461 - t508 * qJD(4) + t162;
t220 = -t368 * t292 + t364 * t461;
t319 = -t364 * t481 + t368 * t508;
t272 = qJD(5) * t319 - t368 * t474;
t617 = t220 - t272;
t221 = t292 * t364 + t368 * t461;
t389 = -t364 * t508 - t368 * t481;
t271 = qJD(5) * t389 + t364 * t474;
t673 = t221 - t271;
t470 = qJD(3) * t538;
t456 = t361 * t470;
t612 = t293 - t456;
t496 = pkin(9) * t480 + t357;
t312 = t496 * qJD(1);
t222 = t312 + (qJD(1) * t411 + t399) * pkin(10);
t379 = pkin(2) * t509 + t397;
t227 = qJD(2) * pkin(2) + qJD(1) * t379 + t350;
t299 = (-pkin(2) * t539 - t366 * t537 - pkin(1)) * t362;
t284 = qJD(1) * t299;
t139 = t222 * t365 - t227 * t447 - t284 * t481;
t645 = -qJD(4) - t139;
t445 = pkin(5) * t368 + pkin(12) * t364;
t671 = qJD(5) * t445 - (-pkin(4) - t445) * t229 - t645;
t490 = qJD(5) * t368;
t140 = t538 * t222 + t227 * t467 + t284 * t503;
t103 = -pkin(4) * t228 + t140;
t507 = qJ(4) * t228;
t144 = t229 * t573 + t507;
t68 = t364 * t103 + t368 * t144;
t670 = -pkin(12) * t228 + t573 * t490 + t68;
t320 = pkin(2) * t447 - pkin(10) * t503;
t306 = -pkin(3) * t508 - t320;
t257 = pkin(4) * t503 - pkin(11) * t508 + t306;
t401 = -pkin(3) * t538 - qJ(4) * t365 - pkin(2);
t279 = (-pkin(11) * t538 + t401) * t361;
t491 = qJD(5) * t364;
t622 = t257 * t490 - t279 * t491 + t364 * t678 - t368 * t677;
t410 = qJD(3) * t447;
t349 = pkin(2) * t410;
t616 = pkin(4) * t292 - t474 * t635 + t349 - t674;
t398 = pkin(4) * t229 + t139;
t644 = qJD(4) + t398;
t86 = t294 * t573 + t644;
t181 = -t227 * t361 + t508 * t284;
t408 = -qJ(4) * t229 + t181;
t91 = t228 * t573 + t408;
t34 = t364 * t86 + t368 * t91;
t667 = t675 - t676;
t648 = t672 + t667;
t289 = t294 * qJ(4);
t94 = t103 - t289;
t598 = -t94 * mrSges(6,1) + t34 * mrSges(6,3) - t648;
t194 = (t460 * t508 - t403) * qJD(2) + t228 * qJD(3);
t345 = qJD(2) * t350;
t393 = qJD(2) * t397;
t247 = qJD(1) * t393 + t345;
t270 = t373 * qJD(2);
t248 = qJD(1) * t270;
t446 = qJD(3) * t467;
t381 = -t222 * t470 - t227 * t446 - t365 * t247 + t248 * t447 - t284 * t474;
t426 = qJD(2) * t464;
t301 = qJD(2) * t392;
t295 = qJD(1) * t301;
t479 = t538 * t295;
t50 = -t194 * pkin(4) + (-qJD(1) * t426 - t479) * t361 - t381;
t372 = (qJD(2) * t382 + qJD(3) * t383) * t362;
t195 = qJD(1) * t372 + qJD(3) * t394;
t202 = -t248 * t361 + t508 * t295;
t387 = qJ(4) * t194 - qJD(4) * t229 + t202;
t56 = t195 * t573 + t387;
t7 = t364 * t50 + t368 * t56 + t86 * t490 - t491 * t91;
t8 = -qJD(5) * t34 - t364 * t56 + t368 * t50;
t669 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t543 = t294 / 0.2e1;
t547 = t229 / 0.2e1;
t548 = -t229 / 0.2e1;
t33 = -t364 * t91 + t368 * t86;
t631 = t33 * mrSges(6,1);
t658 = Ifges(5,4) - Ifges(4,5);
t668 = -Ifges(4,1) * t548 - Ifges(6,5) * t559 + Ifges(5,2) * t547 - Ifges(6,6) * t561 - Ifges(6,3) * t552 + t543 * t658 + t631;
t557 = -t194 / 0.2e1;
t555 = -t195 / 0.2e1;
t630 = Ifges(5,1) + Ifges(4,3);
t628 = -Ifges(4,6) + Ifges(5,5);
t666 = -pkin(12) * t612 + t622;
t665 = -pkin(5) * t617 + pkin(12) * t673 + t616;
t475 = qJD(2) * t502;
t455 = qJD(1) * t475;
t420 = t361 * t455;
t108 = qJD(5) * t192 + t195 * t364 + t368 * t420;
t109 = qJD(5) * t193 - t368 * t195 + t364 * t420;
t44 = qJD(6) * t146 + t108 * t367 - t194 * t363;
t45 = -qJD(6) * t147 - t108 * t363 - t194 * t367;
t13 = Ifges(7,5) * t44 + Ifges(7,6) * t45 + Ifges(7,3) * t109;
t79 = -t222 * t493 + t227 * t410 + t538 * t247 + t248 * t467 + t284 * t456 + t295 * t503;
t72 = -qJ(4) * t420 + t294 * qJD(4) - t79;
t46 = -pkin(4) * t195 - t72;
t556 = t194 / 0.2e1;
t570 = t109 / 0.2e1;
t571 = -t109 / 0.2e1;
t584 = t45 / 0.2e1;
t585 = t44 / 0.2e1;
t32 = pkin(12) * t226 + t34;
t51 = -pkin(5) * t192 - pkin(12) * t193 + t94;
t16 = -t32 * t363 + t367 * t51;
t19 = pkin(5) * t109 - pkin(12) * t108 + t46;
t5 = -pkin(12) * t194 + t7;
t1 = qJD(6) * t16 + t19 * t363 + t367 * t5;
t17 = t32 * t367 + t363 * t51;
t2 = -qJD(6) * t17 + t19 * t367 - t363 * t5;
t602 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t664 = t602 + t46 * mrSges(6,1) + Ifges(7,5) * t585 + Ifges(7,6) * t584 + Ifges(7,3) * t570 - t108 * Ifges(6,4) / 0.2e1 + t13 / 0.2e1 + (-t557 + t556) * Ifges(6,6) + (-t571 + t570) * Ifges(6,2);
t647 = t16 * mrSges(7,1) - t17 * mrSges(7,2);
t35 = Ifges(6,5) * t108 - Ifges(6,6) * t109 - Ifges(6,3) * t194;
t572 = t108 / 0.2e1;
t75 = (-pkin(3) * t455 - t479) * t361 - t381;
t80 = t361 * t479 + t381;
t663 = Ifges(6,5) * t572 + Ifges(6,6) * t571 + t75 * mrSges(5,1) + t35 / 0.2e1 - t80 * mrSges(4,3) + (Ifges(5,2) + Ifges(4,1)) * t557 + (Ifges(5,6) + Ifges(4,4)) * t555 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t420 + t669;
t118 = pkin(3) * t228 + t408;
t119 = pkin(3) * t294 - t645;
t662 = mrSges(5,1) * t119 + t181 * mrSges(4,2) - t34 * mrSges(6,2) - t118 * mrSges(5,3);
t549 = t228 / 0.2e1;
t550 = -t228 / 0.2e1;
t661 = -Ifges(4,2) * t549 + Ifges(5,3) * t550 + t543 * t628;
t659 = t46 * mrSges(6,2);
t225 = Ifges(4,4) * t228;
t653 = t229 * Ifges(4,1) - t294 * Ifges(4,5) + t193 * Ifges(6,5) + t192 * Ifges(6,6) + t226 * Ifges(6,3) - t225;
t340 = pkin(5) * t364 - pkin(12) * t368 + qJ(4);
t500 = t364 * t573;
t298 = t363 * t340 - t367 * t500;
t652 = -qJD(6) * t298 + t363 * t670 + t367 * t671;
t297 = t367 * t340 + t363 * t500;
t651 = qJD(6) * t297 + t363 * t671 - t367 * t670;
t615 = t364 * t257 + t368 * t279;
t623 = -qJD(5) * t615 + t364 * t677 + t368 * t678;
t322 = pkin(10) * t481 + t356;
t650 = t322 * qJD(3) + t396;
t534 = t33 * mrSges(6,3);
t185 = Ifges(6,4) * t192;
t515 = t226 * Ifges(6,5);
t516 = t193 * Ifges(6,1);
t97 = t185 + t515 + t516;
t575 = -t97 / 0.2e1;
t649 = t534 + t575;
t431 = t33 * t364 - t34 * t368;
t646 = qJD(5) * t431;
t459 = t361 * t475;
t418 = t459 / 0.2e1;
t402 = qJD(1) * t418;
t554 = t195 / 0.2e1;
t597 = t80 * mrSges(4,1) - t79 * mrSges(4,2) + t75 * mrSges(5,2) - t72 * mrSges(5,3);
t621 = t194 * t658 + t195 * t628 + t420 * t630;
t643 = t621 / 0.2e1 + Ifges(5,4) * t556 + Ifges(4,5) * t557 + Ifges(5,5) * t554 + Ifges(4,6) * t555 + t402 * t630 + t597;
t642 = Ifges(6,1) * t572 + Ifges(6,5) * t557;
t121 = t289 - t140;
t641 = -t181 * mrSges(4,1) - t121 * mrSges(5,1) + t118 * mrSges(5,2);
t586 = Ifges(6,4) * t571 + t642;
t563 = -t186 / 0.2e1;
t566 = -t147 / 0.2e1;
t568 = -t146 / 0.2e1;
t640 = -Ifges(7,5) * t566 - Ifges(7,6) * t568 - Ifges(7,3) * t563;
t639 = -mrSges(6,3) * t8 + 0.2e1 * t586;
t180 = pkin(12) * t503 + t615;
t304 = -t508 * qJ(4) - t322;
t278 = pkin(4) * t481 - t304;
t203 = -pkin(5) * t389 - pkin(12) * t319 + t278;
t126 = t180 * t367 + t203 * t363;
t627 = -qJD(6) * t126 - t363 * t666 + t367 * t665;
t125 = -t180 * t363 + t203 * t367;
t626 = qJD(6) * t125 + t363 * t665 + t367 * t666;
t434 = t16 * t367 + t17 * t363;
t405 = t434 * mrSges(7,3);
t624 = pkin(5) * t612 - t623;
t620 = t228 * t628 - t229 * t658 - t294 * t630;
t273 = -t319 * t363 + t367 * t503;
t175 = qJD(6) * t273 + t367 * t271 + t363 * t456;
t178 = t221 * t367 + t293 * t363;
t619 = t175 - t178;
t274 = t319 * t367 + t363 * t503;
t176 = -qJD(6) * t274 - t363 * t271 + t367 * t456;
t177 = -t221 * t363 + t293 * t367;
t618 = t176 - t177;
t614 = -pkin(3) * t292 + (-qJ(4) * t470 - t492) * t361 + t679;
t613 = t292 - t474;
t313 = -pkin(10) * t474 + t349;
t611 = -t313 + t674;
t610 = t313 - t162;
t609 = -t361 * t478 - t650;
t608 = -(-pkin(3) * t476 - t478) * t361 + t650;
t465 = t509 * t361;
t607 = -t365 * t502 + t465 * t538 + t395;
t606 = t431 * t229 - t364 * t7 - t368 * t8;
t444 = t1 * t367 - t2 * t363;
t512 = t229 * Ifges(5,6);
t153 = -t294 * Ifges(5,5) + t228 * Ifges(5,3) - t512;
t513 = t229 * Ifges(4,4);
t156 = -t228 * Ifges(4,2) - t294 * Ifges(4,6) + t513;
t604 = -t153 / 0.2e1 + t156 / 0.2e1;
t31 = -pkin(5) * t226 - t33;
t441 = mrSges(7,1) * t363 + mrSges(7,2) * t367;
t540 = t367 / 0.2e1;
t541 = -t363 / 0.2e1;
t524 = Ifges(7,4) * t147;
t60 = t146 * Ifges(7,2) + t186 * Ifges(7,6) + t524;
t145 = Ifges(7,4) * t146;
t61 = t147 * Ifges(7,1) + t186 * Ifges(7,5) + t145;
t603 = -t31 * t441 - t540 * t61 - t541 * t60;
t449 = t365 * t465;
t259 = t449 + t378;
t317 = -t361 * t480 + t508 * t509;
t251 = (t411 + t465) * pkin(10) + t496;
t265 = t358 + t379;
t377 = t365 * t251 - t265 * t447 - t299 * t481;
t100 = t259 * pkin(4) - t317 * t573 + t377;
t205 = -t265 * t361 + t508 * t299;
t407 = -qJ(4) * t259 + t205;
t110 = -t573 * t607 + t407;
t511 = t364 * t100 + t368 * t110;
t207 = t607 * qJD(3) + (-t415 + t421) * qJD(2);
t351 = qJD(2) * t358;
t269 = t351 + t393;
t380 = -t251 * t470 - t265 * t446 - t365 * t269 + t270 * t447 - t299 * t474;
t477 = t538 * t301;
t70 = t207 * pkin(4) + (-t477 - t426) * t361 - t380;
t206 = qJD(3) * t449 + t372;
t209 = -t270 * t361 + t508 * t301;
t386 = -qJ(4) * t207 - qJD(4) * t259 + t209;
t73 = t206 * t573 + t386;
t12 = -qJD(5) * t511 - t364 * t73 + t368 * t70;
t302 = -pkin(9) * t455 + t345;
t316 = t496 * qJD(2);
t303 = qJD(1) * t316;
t600 = -t303 * mrSges(3,1) - t302 * mrSges(3,2);
t599 = t405 - t515 / 0.2e1 - t185 / 0.2e1 - t94 * mrSges(6,2) + t649;
t596 = -t139 * mrSges(4,1) - t140 * mrSges(4,2) + t119 * mrSges(5,2) - t121 * mrSges(5,3);
t551 = t226 / 0.2e1;
t558 = t193 / 0.2e1;
t560 = t192 / 0.2e1;
t595 = -Ifges(6,4) * t558 - Ifges(6,2) * t560 - Ifges(6,6) * t551 - t598 + t675;
t593 = -mrSges(6,3) * t7 - Ifges(6,4) * t572 + t664;
t592 = t362 ^ 2;
t14 = t44 * Ifges(7,4) + t45 * Ifges(7,2) + t109 * Ifges(7,6);
t589 = t14 / 0.2e1;
t588 = Ifges(7,1) * t585 + Ifges(7,4) * t584 + Ifges(7,5) * t570;
t581 = -t60 / 0.2e1;
t580 = t60 / 0.2e1;
t579 = -t61 / 0.2e1;
t578 = t61 / 0.2e1;
t574 = t97 / 0.2e1;
t544 = -t294 / 0.2e1;
t18 = -mrSges(7,1) * t45 + mrSges(7,2) * t44;
t82 = -mrSges(6,1) * t194 - mrSges(6,3) * t108;
t532 = -t18 + t82;
t531 = mrSges(4,3) * t228;
t530 = mrSges(4,3) * t229;
t529 = Ifges(3,4) * t366;
t528 = Ifges(4,4) * t365;
t526 = Ifges(6,4) * t364;
t525 = Ifges(6,4) * t368;
t523 = Ifges(7,4) * t363;
t522 = Ifges(7,4) * t367;
t521 = Ifges(5,6) * t365;
t149 = mrSges(6,1) * t226 - mrSges(6,3) * t193;
t78 = -mrSges(7,1) * t146 + mrSges(7,2) * t147;
t510 = -t78 + t149;
t506 = t202 * t361;
t505 = t229 * t364;
t501 = t363 * t368;
t499 = t367 * t368;
t498 = t368 * t573;
t127 = -mrSges(6,1) * t192 + mrSges(6,2) * t193;
t199 = mrSges(5,1) * t228 + mrSges(5,3) * t294;
t497 = t127 - t199;
t494 = qJD(3) * t361;
t488 = t539 * Ifges(3,4);
t487 = t538 * Ifges(4,4);
t486 = t538 * Ifges(5,6);
t160 = t538 * t251 + t265 * t467 + t299 * t503;
t471 = qJD(2) * t539;
t169 = -t194 * mrSges(5,1) + mrSges(5,2) * t420;
t463 = mrSges(3,3) * t476;
t142 = -t317 * qJ(4) - t160;
t457 = t362 * t471;
t454 = -t481 / 0.2e1;
t453 = t481 / 0.2e1;
t450 = qJD(1) * t471;
t442 = mrSges(7,1) * t367 - mrSges(7,2) * t363;
t440 = Ifges(7,1) * t367 - t523;
t439 = Ifges(7,1) * t363 + t522;
t438 = -Ifges(7,2) * t363 + t522;
t437 = Ifges(7,2) * t367 + t523;
t436 = Ifges(7,5) * t367 - Ifges(7,6) * t363;
t435 = Ifges(7,5) * t363 + Ifges(7,6) * t367;
t433 = t16 * t363 - t17 * t367;
t29 = mrSges(7,1) * t109 - mrSges(7,3) * t44;
t30 = -mrSges(7,2) * t109 + mrSges(7,3) * t45;
t432 = -t363 * t29 + t367 * t30;
t41 = pkin(12) * t259 + t511;
t117 = pkin(4) * t607 - t142;
t211 = t317 * t368 - t364 * t607;
t422 = -t317 * t364 - t368 * t607;
t74 = -pkin(5) * t422 - pkin(12) * t211 + t117;
t21 = t363 * t74 + t367 * t41;
t20 = -t363 * t41 + t367 * t74;
t92 = -mrSges(7,2) * t186 + mrSges(7,3) * t146;
t93 = mrSges(7,1) * t186 - mrSges(7,3) * t147;
t430 = -t363 * t92 - t367 * t93;
t427 = mrSges(3,3) * t458;
t47 = t100 * t368 - t110 * t364;
t67 = t103 * t368 - t144 * t364;
t164 = t211 * t367 + t259 * t363;
t163 = -t211 * t363 + t259 * t367;
t196 = t257 * t368 - t279 * t364;
t416 = t362 * t450;
t11 = t100 * t490 - t110 * t491 + t364 * t70 + t368 * t73;
t400 = Ifges(3,5) * t416 - Ifges(3,6) * t455;
t391 = pkin(1) * t592 * (mrSges(3,1) * t366 + mrSges(3,2) * t539);
t390 = t366 * t592 * (Ifges(3,1) * t539 - t529);
t89 = -t251 * t493 + t265 * t410 + t538 * t269 + t270 * t467 + t299 * t456 + t301 * t503;
t81 = -qJ(4) * t459 - t317 * qJD(4) - t89;
t64 = -pkin(4) * t206 - t81;
t376 = t362 * t417 * (Ifges(3,5) * t539 - Ifges(3,6) * t366);
t375 = (Ifges(3,6) * t509 + (Ifges(3,2) * t539 + t529) * t362) * qJD(1);
t371 = t598 - t647 - t667;
t370 = t436 * t562 + t438 * t567 + t440 * t565 - t603;
t347 = Ifges(3,4) * t458;
t321 = -pkin(9) * t502 + t358;
t315 = -pkin(9) * t475 + t351;
t311 = -pkin(9) * t476 + t350;
t310 = -mrSges(3,2) * t417 + t427;
t309 = mrSges(3,1) * t417 - t463;
t305 = t401 * t361;
t281 = Ifges(3,1) * t476 + Ifges(3,5) * t417 + t347;
t280 = Ifges(3,6) * qJD(2) + t375;
t224 = Ifges(5,6) * t228;
t201 = -mrSges(4,1) * t294 - t530;
t200 = mrSges(5,1) * t229 - mrSges(5,2) * t294;
t198 = mrSges(4,2) * t294 - t531;
t184 = -t229 * t499 - t294 * t363;
t183 = t229 * t501 - t294 * t367;
t179 = -pkin(5) * t503 - t196;
t174 = -mrSges(5,2) * t228 - mrSges(5,3) * t229;
t173 = mrSges(4,1) * t228 + mrSges(4,2) * t229;
t172 = pkin(3) * t229 + t507;
t171 = -mrSges(4,2) * t420 - mrSges(4,3) * t195;
t170 = mrSges(4,1) * t420 + mrSges(4,3) * t194;
t168 = mrSges(5,1) * t195 - mrSges(5,3) * t420;
t167 = -t228 * t363 + t367 * t505;
t166 = -t228 * t367 - t363 * t505;
t155 = -t294 * Ifges(5,4) - t229 * Ifges(5,2) + t224;
t148 = -mrSges(6,2) * t226 + mrSges(6,3) * t192;
t143 = -t317 * pkin(3) + t377;
t138 = -pkin(3) * t607 + t407;
t137 = qJD(5) * t422 + t206 * t364 + t368 * t459;
t136 = qJD(5) * t211 - t206 * t368 + t364 * t459;
t130 = mrSges(4,1) * t195 - mrSges(4,2) * t194;
t129 = -mrSges(5,2) * t195 + mrSges(5,3) * t194;
t128 = pkin(5) * t193 - pkin(12) * t192;
t115 = -Ifges(4,4) * t194 - Ifges(4,2) * t195 + Ifges(4,6) * t420;
t111 = Ifges(5,5) * t420 + Ifges(5,6) * t194 + Ifges(5,3) * t195;
t90 = t361 * t477 + t380;
t88 = pkin(3) * t206 + t386;
t87 = (-pkin(3) * t475 - t477) * t361 - t380;
t83 = mrSges(6,2) * t194 - mrSges(6,3) * t109;
t77 = pkin(3) * t195 + t387;
t63 = qJD(6) * t163 + t137 * t367 + t207 * t363;
t62 = -qJD(6) * t164 - t137 * t363 + t207 * t367;
t54 = mrSges(6,1) * t109 + mrSges(6,2) * t108;
t52 = pkin(5) * t228 - t67;
t40 = -pkin(5) * t259 - t47;
t26 = pkin(5) * t136 - pkin(12) * t137 + t64;
t25 = t128 * t363 + t33 * t367;
t24 = t128 * t367 - t33 * t363;
t10 = -pkin(5) * t207 - t12;
t9 = pkin(12) * t207 + t11;
t6 = pkin(5) * t194 - t8;
t4 = -qJD(6) * t21 + t26 * t367 - t363 * t9;
t3 = qJD(6) * t20 + t26 * t363 + t367 * t9;
t15 = [m(4) * (-t139 * t90 + t140 * t89 + t160 * t79 + t181 * t209 + t202 * t205 - t377 * t80) - t377 * t170 - t593 * t422 + (t302 * t480 + t303 * t502 - t311 * t457 - t312 * t475 - t321 * t416 - t455 * t496) * mrSges(3,3) + m(3) * (t302 * t496 - t303 * t321 - t311 * t316 + t312 * t315) + ((t390 - 0.2e1 * t391) * qJD(1) + t376 / 0.2e1) * qJD(2) + t639 * t211 + (t653 / 0.2e1 + Ifges(6,3) * t551 + t139 * mrSges(4,3) + Ifges(6,6) * t560 + Ifges(6,5) * t558 + Ifges(4,1) * t547 - Ifges(5,2) * t548 - Ifges(5,6) * t549 + Ifges(4,4) * t550 + t631 - t155 / 0.2e1 - t658 * t544 + t662) * t207 + (t202 * mrSges(4,2) - t77 * mrSges(5,3) + Ifges(4,4) * t555 - Ifges(5,2) * t556 - Ifges(5,6) * t554 + (Ifges(4,1) + Ifges(6,3)) * t557 - t658 * t402 + t663) * t259 + (t281 + qJD(1) * (Ifges(3,5) * t509 + (Ifges(3,1) * t366 + t488) * t362)) * t457 / 0.2e1 - (t375 + t280) * t475 / 0.2e1 + m(6) * (t11 * t34 + t117 * t46 + t12 * t33 + t47 * t8 + t511 * t7 + t64 * t94) + t511 * t83 + (t137 * t94 + t211 * t46) * mrSges(6,2) + (Ifges(7,5) * t63 + Ifges(7,6) * t62) * t562 + (Ifges(7,5) * t164 + Ifges(7,6) * t163) * t570 + t63 * t578 + t62 * t580 + (t1 * t163 - t16 * t63 - t164 * t2 + t17 * t62) * mrSges(7,3) + m(7) * (t1 * t21 + t10 * t31 + t16 * t4 + t17 * t3 + t2 * t20 + t40 * t6) + m(5) * (t118 * t88 + t119 * t87 + t121 * t81 + t138 * t77 + t142 * t72 + t143 * t75) + t643 * t317 - (Ifges(5,3) * t554 - Ifges(4,2) * t555 + Ifges(5,6) * t556 - Ifges(4,4) * t557 + t72 * mrSges(5,1) - t79 * mrSges(4,3) - t77 * mrSges(5,2) + t111 / 0.2e1 - t115 / 0.2e1 + t202 * mrSges(4,1) + t628 * t402) * t607 + (t595 + t647) * t136 + (Ifges(7,4) * t63 + Ifges(7,2) * t62) * t567 + (Ifges(7,4) * t164 + Ifges(7,2) * t163) * t584 + (-t140 * mrSges(4,3) - Ifges(4,4) * t547 - Ifges(4,2) * t550 + Ifges(5,6) * t548 + Ifges(5,3) * t549 + t544 * t628 - t604 - t641) * t206 + t164 * t588 + t163 * t589 + (Ifges(7,1) * t63 + Ifges(7,4) * t62) * t565 + (Ifges(7,1) * t164 + Ifges(7,4) * t163) * t585 + t592 * (-Ifges(3,2) * t366 + t488) * t450 + (Ifges(6,1) * t558 + Ifges(6,4) * t560 + Ifges(6,5) * t551 - t534 + t574) * t137 + (Ifges(5,4) * t548 + Ifges(4,5) * t547 + Ifges(5,5) * t549 + Ifges(4,6) * t550 + t544 * t630 + t596) * t459 + t620 * t418 + (t400 / 0.2e1 + t600) * t509 + t315 * t310 - t316 * t309 + t20 * t29 + t21 * t30 + t40 * t18 + t31 * (-mrSges(7,1) * t62 + mrSges(7,2) * t63) + t10 * t78 + t47 * t82 + t3 * t92 + t4 * t93 + t117 * t54 + t64 * t127 + t138 * t129 + t11 * t148 + t12 * t149 + t6 * (-mrSges(7,1) * t163 + mrSges(7,2) * t164) + t142 * t168 + t143 * t169 + t160 * t171 + t88 * t174 + t89 * t198 + t81 * t199 + t87 * t200 + t90 * t201 + t205 * t130 + t209 * t173; (t196 * t8 + t278 * t46 + t33 * t623 + t34 * t622 + t615 * t7 + t616 * t94) * m(6) + t615 * t83 - t593 * t389 + (Ifges(7,1) * t175 + Ifges(7,4) * t176) * t565 + (Ifges(7,1) * t178 + Ifges(7,4) * t177) * t566 + (Ifges(7,1) * t274 + Ifges(7,4) * t273) * t585 + t653 * (qJD(3) * t453 - t293 / 0.2e1) + ((-mrSges(5,2) * t365 - mrSges(5,3) * t538) * t494 + mrSges(5,2) * t292 + mrSges(5,3) * t293) * t118 + (Ifges(6,1) * t271 + Ifges(6,5) * t456) * t558 + (-t139 * t612 + t140 * t613 + t481 * t79) * mrSges(4,3) + (-t119 * t612 - t121 * t613 - t481 * t72) * mrSges(5,1) + (Ifges(6,5) * t271 + Ifges(6,3) * t456) * t551 - (t229 * (-Ifges(5,2) * t538 + t521) + t228 * (-Ifges(4,2) * t365 + t487) + t365 * t156 + (t365 * t628 - t538 * t658) * t294) * t494 / 0.2e1 + ((-Ifges(5,3) * t538 - t521) * t554 + (Ifges(4,2) * t538 + t528) * t555 + (-Ifges(5,2) * t365 - t486) * t556 + (Ifges(4,1) * t365 + t487) * t557 + t77 * (mrSges(5,2) * t538 - mrSges(5,3) * t365) - pkin(2) * t130 + (-t365 * t658 - t538 * t628) * t402) * t361 + (-t376 / 0.2e1 + (-t390 / 0.2e1 + t391) * qJD(1)) * qJD(1) + (Ifges(6,4) * t271 + Ifges(6,6) * t456) * t560 - (-Ifges(3,2) * t476 + t281 + t347) * t458 / 0.2e1 + (t229 * (Ifges(4,1) * t538 - t528) + t228 * (Ifges(5,3) * t365 - t486) + t365 * t153) * t494 / 0.2e1 + (Ifges(4,4) * t549 - Ifges(5,6) * t550 - t668) * t293 + (Ifges(7,4) * t175 + Ifges(7,2) * t176) * t567 + (Ifges(7,4) * t178 + Ifges(7,2) * t177) * t568 + (Ifges(7,4) * t274 + Ifges(7,2) * t273) * t584 + t33 * (mrSges(6,1) * t456 - t271 * mrSges(6,3)) + t271 * t574 + t175 * t578 + t178 * t579 + t176 * t580 + (-mrSges(4,1) * t538 + mrSges(4,2) * t365) * t506 + (t309 + t463) * t312 + t400 + t643 * t508 + (t1 * t273 - t2 * t274) * mrSges(7,3) + t280 * t476 / 0.2e1 + (Ifges(6,1) * t559 + Ifges(6,4) * t561 + Ifges(6,5) * t552 + t649) * t221 + (-t310 + t427) * t311 + (t293 / 0.2e1 + qJD(3) * t454) * t155 + t177 * t581 + t274 * t588 + t273 * t589 + (-Ifges(4,4) * t548 + Ifges(5,6) * t547 + t604 + t661) * t292 + (Ifges(5,4) * t547 + Ifges(4,5) * t548 + Ifges(5,5) * t550 + Ifges(4,6) * t549 + t543 * t630 - t596 - t620 / 0.2e1) * t461 + ((mrSges(4,1) * t365 + mrSges(4,2) * t538) * t494 - mrSges(4,1) * t292 - mrSges(4,2) * t293) * t181 + (Ifges(7,5) * t175 + Ifges(7,6) * t176) * t562 + (Ifges(7,5) * t178 + Ifges(7,6) * t177) * t563 + (Ifges(7,5) * t274 + Ifges(7,6) * t273) * t570 + t600 + t626 * t92 + t627 * t93 + (t1 * t126 + t125 * t2 + t16 * t627 + t17 * t626 + t179 * t6 + t31 * t624) * m(7) + t622 * t148 + t623 * t149 + t624 * t78 + t614 * t174 + (t118 * t614 + t119 * t608 + t121 * t611 + t304 * t72 + t305 * t77 + t306 * t75) * m(5) + t616 * t127 + (mrSges(7,2) * t617 + mrSges(7,3) * t618) * t17 + (-mrSges(7,1) * t618 + mrSges(7,2) * t619) * t31 + (-mrSges(7,1) * t617 - mrSges(7,3) * t619) * t16 + t608 * t200 + t609 * t201 + t610 * t198 + (-pkin(2) * t506 - t139 * t609 + t140 * t610 - t181 * t208 + t320 * t80 + t322 * t79) * m(4) + t611 * t199 + t322 * t171 + t320 * t170 + t304 * t168 + t305 * t129 + t306 * t169 + t278 * t54 + t6 * (-mrSges(7,1) * t273 + mrSges(7,2) * t274) + t595 * t272 + (t639 + t659) * t319 + (Ifges(6,3) * t557 + t663) * t503 + t115 * t453 + t111 * t454 + (t612 * t34 - t673 * t94) * mrSges(6,2) + (t598 - t640 + t676 - t672) * t220 + t125 * t29 + t126 * t30 + t179 * t18 + t196 * t82 - t208 * t173; t398 * t127 + t652 * t93 + (t1 * t298 + t16 * t652 + t17 * t651 + t2 * t297 - t31 * t52 + t498 * t6) * m(7) + (t436 * t570 - t532 * t573 + t438 * t584 + t440 * t585 + t586 + (mrSges(7,3) * t433 + t31 * t442 + t367 * t581 + t435 * t563 + t437 * t568 + t439 * t566 + t541 * t61) * qJD(6) + t659 + (-t148 * t573 - t371 + t672) * qJD(5) + t6 * t441 + t642) * t368 + t597 + (t662 + t668) * t228 + (t198 - t199 + t531) * t139 - t14 * t501 / 0.2e1 + t497 * qJD(4) + (-t1 * t501 + t16 * t167 - t166 * t17 - t2 * t499) * mrSges(7,3) + (t512 + t156) * t547 + t505 * t575 + t167 * t579 + t621 + (-pkin(3) * t75 - qJ(4) * t72 - t118 * t172 - t119 * t140 + t121 * t645) * m(5) + (t46 * qJ(4) - t33 * t67 - t34 * t68 - t498 * t8 - t500 * t7 + t573 * t646 + t644 * t94) * m(6) + (t224 + t155) * t550 + (-t513 + t153) * t548 + (Ifges(7,5) * t167 + Ifges(7,6) * t166) * t563 + t651 * t92 + t166 * t581 + t499 * t588 + (Ifges(7,4) * t167 + Ifges(7,2) * t166) * t568 + (-t225 + t653) * t549 + t525 * t571 + (-t200 + t201 + t530) * t140 + (t54 - t168) * qJ(4) + (Ifges(7,1) * t167 + Ifges(7,4) * t166) * t566 - t526 * t572 + t606 * mrSges(6,3) + t297 * t29 + t298 * t30 + ((t640 + t648 + t647) * t368 + (Ifges(6,5) * t364 + Ifges(6,6) * t368) * t552 + (Ifges(6,1) * t364 + t525) * t559 + (Ifges(6,2) * t368 + t526) * t561 - t94 * (-mrSges(6,1) * t368 + mrSges(6,2) * t364) + t641 + t661) * t229 + (-t83 * t573 + (t438 * t568 + t440 * t566 + t436 * t563 - t516 / 0.2e1 - (m(7) * t31 - t510) * t573 + t599 + t603) * qJD(5) + t664) * t364 - t52 * t78 - t68 * t148 - t67 * t149 - t31 * (-mrSges(7,1) * t166 + mrSges(7,2) * t167) - pkin(3) * t169 - t172 * t174; t229 * t174 - t183 * t93 - t184 * t92 + t497 * t294 + (t229 * t148 + (-t363 * t93 + t367 * t92 + t148) * qJD(5) + t532) * t368 + (qJD(6) * t430 - t226 * t510 + t432 + t83) * t364 + t169 + ((-qJD(5) * t433 - t6) * t368 + (qJD(5) * t31 - qJD(6) * t434 + t444) * t364 - t16 * t183 - t17 * t184 + t31 * t505) * m(7) + (t294 * t94 - t606 - t646) * m(6) + (t118 * t229 - t121 * t294 + t75) * m(5); (t370 - t405) * qJD(6) + t435 * t570 + t444 * mrSges(7,3) + t510 * t34 + t437 * t584 + t439 * t585 - t6 * t442 + t35 + ((Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t193 - t370 + t599) * t192 + t371 * t193 + t14 * t540 + t363 * t588 - pkin(5) * t18 - t25 * t92 - t24 * t93 - t33 * t148 + (-pkin(5) * t6 - t16 * t24 - t17 * t25 - t31 * t34) * m(7) + (t432 + (-m(7) * t434 + t430) * qJD(6) + m(7) * t444) * pkin(12) + t669; -t31 * (mrSges(7,1) * t147 + mrSges(7,2) * t146) + (Ifges(7,1) * t146 - t524) * t566 + t60 * t565 + (Ifges(7,5) * t146 - Ifges(7,6) * t147) * t563 - t16 * t92 + t17 * t93 + (t146 * t16 + t147 * t17) * mrSges(7,3) + t13 + (-Ifges(7,2) * t147 + t145 + t61) * t568 + t602;];
tauc  = t15(:);
