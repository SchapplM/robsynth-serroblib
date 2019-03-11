% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRR9
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:46
% EndTime: 2019-03-10 05:18:37
% DurationCPUTime: 59.20s
% Computational Cost: add. (60735->1163), mult. (176805->1629), div. (0->0), fcn. (148527->14), ass. (0->519)
t470 = cos(qJ(2));
t460 = cos(pkin(6));
t622 = pkin(1) * t460;
t453 = t470 * t622;
t443 = qJD(1) * t453;
t465 = sin(qJ(2));
t458 = sin(pkin(6));
t459 = cos(pkin(7));
t528 = pkin(10) * t459 + pkin(9);
t509 = t458 * t528;
t491 = t465 * t509;
t356 = -qJD(1) * t491 + t443;
t452 = t465 * t622;
t477 = -t470 * t509 - t452;
t357 = t477 * qJD(1);
t457 = sin(pkin(7));
t620 = pkin(10) * t457;
t480 = (pkin(2) * t465 - t470 * t620) * t458;
t392 = qJD(1) * t480;
t464 = sin(qJ(3));
t566 = t457 * t464;
t447 = pkin(10) * t566;
t469 = cos(qJ(3));
t561 = t459 * t469;
t418 = pkin(2) * t561 - t447;
t562 = t459 * t464;
t686 = t418 * qJD(3) - t469 * t356 - t357 * t562 - t392 * t566;
t299 = -t357 * t457 + t459 * t392;
t555 = t465 * t469;
t556 = t464 * t470;
t487 = t459 * t555 + t556;
t539 = qJD(1) * t458;
t377 = t487 * t539;
t553 = t469 * t470;
t557 = t464 * t465;
t485 = -t459 * t557 + t553;
t378 = t485 * t539;
t753 = -pkin(3) * t377 + pkin(11) * t378 - t299 + (pkin(3) * t464 - pkin(11) * t469) * t457 * qJD(3);
t526 = t465 * t539;
t512 = t457 * t526;
t752 = pkin(11) * t512 - t686;
t446 = qJD(1) * t460 + qJD(2);
t488 = t459 * t553 - t557;
t565 = t457 * t469;
t325 = t446 * t565 + t488 * t539;
t323 = qJD(4) - t325;
t463 = sin(qJ(4));
t468 = cos(qJ(4));
t317 = t378 * t463 - t468 * t512;
t415 = t459 * t463 + t468 * t566;
t536 = qJD(3) * t469;
t522 = t457 * t536;
t361 = qJD(4) * t415 + t463 * t522;
t742 = t317 - t361;
t318 = t378 * t468 + t463 * t512;
t414 = -t468 * t459 + t463 * t566;
t360 = -qJD(4) * t414 + t468 * t522;
t751 = t318 - t360;
t420 = pkin(2) * t562 + pkin(10) * t565;
t750 = -t420 * qJD(3) + t464 * t356;
t537 = qJD(3) * t464;
t523 = t457 * t537;
t749 = t377 - t523;
t563 = t458 * t470;
t319 = t446 * t620 + (t528 * t563 + t452) * qJD(1);
t324 = pkin(2) * t446 + t356;
t388 = (-pkin(2) * t470 - t465 * t620 - pkin(1)) * t458;
t371 = qJD(1) * t388;
t495 = t324 * t459 + t371 * t457;
t225 = t319 * t469 + t464 * t495;
t748 = -t225 + t323 * (pkin(4) * t463 - pkin(12) * t468);
t400 = pkin(11) * t459 + t420;
t401 = (-pkin(3) * t469 - pkin(11) * t464 - pkin(2)) * t457;
t534 = qJD(4) * t468;
t535 = qJD(4) * t463;
t688 = -t400 * t535 + t401 * t534 + t463 * t753 - t752 * t468;
t684 = t357 * t561 - (-pkin(3) * t526 - t392 * t469) * t457 - t750;
t747 = pkin(12) * t749 - t688;
t746 = -pkin(4) * t742 + t751 * pkin(12) + t684;
t486 = t459 * t556 + t555;
t478 = t486 * t458;
t326 = qJD(1) * t478 + t446 * t566;
t467 = cos(qJ(5));
t462 = sin(qJ(5));
t559 = t462 * t468;
t256 = -t325 * t559 + t326 * t467;
t532 = qJD(5) * t467;
t745 = t462 * t534 + t463 * t532 + t256;
t525 = t470 * t539;
t382 = t459 * t446 - t457 * t525 + qJD(3);
t279 = -t463 * t326 + t382 * t468;
t475 = (qJD(2) * t485 + qJD(3) * t488) * t458;
t281 = qJD(1) * t475 + t446 * t522;
t538 = qJD(2) * t458;
t517 = qJD(1) * t538;
t510 = t465 * t517;
t492 = t457 * t510;
t191 = qJD(4) * t279 + t281 * t468 + t463 * t492;
t655 = t191 / 0.2e1;
t744 = Ifges(5,4) * t655;
t224 = -t464 * t319 + t469 * t495;
t265 = pkin(3) * t326 - pkin(11) * t325;
t163 = t468 * t224 + t463 * t265;
t135 = pkin(12) * t326 + t163;
t437 = -pkin(4) * t468 - pkin(12) * t463 - pkin(3);
t533 = qJD(5) * t462;
t724 = -t467 * t135 + t437 * t532 + (-t467 * t535 - t468 * t533) * pkin(11) + t748 * t462;
t619 = pkin(11) * t462;
t743 = t135 * t462 + t467 * t748 + t535 * t619;
t362 = -t415 * t462 - t467 * t565;
t267 = qJD(5) * t362 + t360 * t467 + t462 * t523;
t270 = t318 * t467 + t377 * t462;
t545 = t267 - t270;
t489 = -t415 * t467 + t462 * t565;
t268 = qJD(5) * t489 - t360 * t462 + t467 * t523;
t269 = -t318 * t462 + t377 * t467;
t544 = t268 - t269;
t280 = t326 * t468 + t382 * t463;
t192 = qJD(4) * t280 + t281 * t463 - t468 * t492;
t654 = -t192 / 0.2e1;
t474 = (qJD(2) * t487 + qJD(3) * t486) * t458;
t282 = qJD(1) * t474 + t446 * t523;
t635 = t282 / 0.2e1;
t554 = t467 * t468;
t257 = t325 * t554 + t326 * t462;
t454 = pkin(11) * t554;
t568 = t325 * t463;
t741 = -pkin(5) * t568 + pkin(13) * t257 + (pkin(5) * t463 - pkin(13) * t554) * qJD(4) + (-t454 + (pkin(13) * t463 - t437) * t462) * qJD(5) + t743;
t740 = -t745 * pkin(13) + t724;
t399 = t447 + (-pkin(2) * t469 - pkin(3)) * t459;
t306 = pkin(4) * t414 - pkin(12) * t415 + t399;
t312 = t468 * t400 + t463 * t401;
t308 = -pkin(12) * t565 + t312;
t231 = t462 * t306 + t467 * t308;
t692 = -qJD(5) * t231 + t747 * t462 + t746 * t467;
t691 = t306 * t532 - t308 * t533 + t746 * t462 - t747 * t467;
t739 = -pkin(5) * t742 - pkin(13) * t545 + t692;
t738 = -pkin(13) * t544 - t691;
t667 = -pkin(13) - pkin(12);
t527 = qJD(5) * t667;
t569 = t279 * t462;
t271 = -t324 * t457 + t459 * t371;
t201 = -pkin(3) * t325 - pkin(11) * t326 + t271;
t206 = pkin(11) * t382 + t225;
t111 = t201 * t468 - t463 * t206;
t211 = pkin(4) * t280 - pkin(12) * t279;
t87 = t467 * t111 + t462 * t211;
t737 = pkin(13) * t569 + t462 * t527 - t87;
t86 = -t111 * t462 + t467 * t211;
t736 = -pkin(5) * t280 - t86 + (pkin(13) * t279 + t527) * t467;
t687 = -t400 * t534 - t401 * t535 + t752 * t463 + t468 * t753;
t101 = -pkin(4) * t323 - t111;
t232 = -t280 * t462 + t323 * t467;
t274 = qJD(5) - t279;
t233 = t280 * t467 + t323 * t462;
t611 = Ifges(6,4) * t233;
t114 = Ifges(6,2) * t232 + Ifges(6,6) * t274 + t611;
t229 = Ifges(6,4) * t232;
t115 = Ifges(6,1) * t233 + Ifges(6,5) * t274 + t229;
t112 = t201 * t463 + t206 * t468;
t102 = pkin(12) * t323 + t112;
t205 = -pkin(3) * t382 - t224;
t122 = -pkin(4) * t279 - pkin(12) * t280 + t205;
t60 = -t102 * t462 + t467 * t122;
t61 = t102 * t467 + t122 * t462;
t497 = t462 * t61 + t467 * t60;
t499 = Ifges(6,5) * t467 - Ifges(6,6) * t462;
t609 = Ifges(6,4) * t467;
t501 = -Ifges(6,2) * t462 + t609;
t610 = Ifges(6,4) * t462;
t503 = Ifges(6,1) * t467 - t610;
t504 = mrSges(6,1) * t462 + mrSges(6,2) * t467;
t624 = t467 / 0.2e1;
t625 = -t462 / 0.2e1;
t640 = t274 / 0.2e1;
t646 = t233 / 0.2e1;
t648 = t232 / 0.2e1;
t735 = t497 * mrSges(6,3) - t101 * t504 - t114 * t625 - t115 * t624 - t499 * t640 - t501 * t648 - t503 * t646;
t359 = t477 * qJD(2);
t340 = qJD(1) * t359;
t393 = qJD(2) * t480;
t383 = qJD(1) * t393;
t438 = qJD(2) * t443;
t482 = qJD(2) * t491;
t339 = -qJD(1) * t482 + t438;
t521 = t459 * t537;
t514 = -t319 * t536 - t324 * t521 - t464 * t339 - t371 * t523;
t137 = -t340 * t561 + (-pkin(3) * t510 - t383 * t469) * t457 - t514;
t653 = t192 / 0.2e1;
t100 = -qJD(5) * t233 - t191 * t462 + t282 * t467;
t663 = t100 / 0.2e1;
t99 = qJD(5) * t232 + t191 * t467 + t282 * t462;
t668 = t99 / 0.2e1;
t461 = sin(qJ(6));
t466 = cos(qJ(6));
t142 = t232 * t461 + t233 * t466;
t39 = -qJD(6) * t142 + t100 * t466 - t461 * t99;
t673 = t39 / 0.2e1;
t516 = t466 * t232 - t233 * t461;
t38 = qJD(6) * t516 + t100 * t461 + t466 * t99;
t674 = t38 / 0.2e1;
t520 = t459 * t536;
t147 = -t319 * t537 + t324 * t520 + t469 * t339 + t340 * t562 + t371 * t522 + t383 * t566;
t136 = pkin(11) * t492 + t147;
t287 = -t340 * t457 + t459 * t383;
t168 = pkin(3) * t282 - pkin(11) * t281 + t287;
t49 = t468 * t136 + t463 * t168 + t201 * t534 - t206 * t535;
t44 = pkin(12) * t282 + t49;
t68 = pkin(4) * t192 - pkin(12) * t191 + t137;
t12 = -t102 * t533 + t122 * t532 + t467 * t44 + t462 * t68;
t13 = -qJD(5) * t61 - t44 * t462 + t467 * t68;
t718 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t51 = -pkin(13) * t233 + t60;
t46 = pkin(5) * t274 + t51;
t52 = pkin(13) * t232 + t61;
t574 = t461 * t52;
t18 = t46 * t466 - t574;
t6 = pkin(5) * t192 - pkin(13) * t99 + t13;
t7 = pkin(13) * t100 + t12;
t2 = qJD(6) * t18 + t461 * t6 + t466 * t7;
t572 = t466 * t52;
t19 = t46 * t461 + t572;
t3 = -qJD(6) * t19 - t461 * t7 + t466 * t6;
t720 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t734 = -t137 * mrSges(5,1) - Ifges(6,5) * t668 - Ifges(7,5) * t674 + 0.2e1 * Ifges(5,2) * t654 + 0.2e1 * Ifges(5,6) * t635 - Ifges(6,6) * t663 - Ifges(7,6) * t673 - t718 - t720 + (-Ifges(7,3) - Ifges(6,3)) * t653 + t744;
t138 = Ifges(7,4) * t516;
t608 = Ifges(7,4) * t142;
t272 = qJD(6) + t274;
t643 = -t272 / 0.2e1;
t659 = -t142 / 0.2e1;
t661 = -t516 / 0.2e1;
t71 = Ifges(7,1) * t142 + Ifges(7,5) * t272 + t138;
t83 = -pkin(5) * t232 + t101;
t9 = Ifges(7,5) * t38 + Ifges(7,6) * t39 + Ifges(7,3) * t192;
t733 = t720 + t9 + (Ifges(7,5) * t516 - Ifges(7,6) * t142) * t643 + (t142 * t19 + t18 * t516) * mrSges(7,3) + (-Ifges(7,2) * t142 + t138 + t71) * t661 - t83 * (mrSges(7,1) * t142 + mrSges(7,2) * t516) + (Ifges(7,1) * t516 - t608) * t659;
t732 = t281 / 0.2e1;
t731 = -t282 / 0.2e1;
t730 = t382 / 0.2e1;
t40 = Ifges(6,5) * t99 + Ifges(6,6) * t100 + Ifges(6,3) * t192;
t704 = t40 + t9;
t728 = Ifges(5,4) * t654;
t727 = t271 * mrSges(4,2);
t595 = t274 * Ifges(6,3);
t597 = t233 * Ifges(6,5);
t598 = t232 * Ifges(6,6);
t113 = t595 + t597 + t598;
t596 = t272 * Ifges(7,3);
t601 = t142 * Ifges(7,5);
t602 = t516 * Ifges(7,6);
t69 = t596 + t601 + t602;
t696 = t113 + t69;
t427 = t467 * t437;
t558 = t463 * t467;
t346 = -pkin(13) * t558 + t427 + (-pkin(5) - t619) * t468;
t395 = t462 * t437 + t454;
t560 = t462 * t463;
t365 = -pkin(13) * t560 + t395;
t283 = t346 * t466 - t365 * t461;
t726 = qJD(6) * t283 + t461 * t741 + t466 * t740;
t284 = t346 * t461 + t365 * t466;
t725 = -qJD(6) * t284 - t461 * t740 + t466 * t741;
t723 = -qJD(5) * t395 + t743;
t689 = pkin(4) * t749 - t687;
t162 = -t463 * t224 + t265 * t468;
t134 = -pkin(4) * t326 - t162;
t722 = t745 * pkin(5) + pkin(11) * t534 - t134;
t273 = Ifges(5,4) * t279;
t587 = t323 * Ifges(5,5);
t592 = t280 * Ifges(5,1);
t180 = t273 + t587 + t592;
t717 = -t205 * mrSges(5,2) + t111 * mrSges(5,3);
t721 = -t180 / 0.2e1 - t587 / 0.2e1 - t273 / 0.2e1 + t717 + t735;
t719 = -t271 * mrSges(4,1) - t111 * mrSges(5,1) + t112 * mrSges(5,2);
t628 = -t382 / 0.2e1;
t630 = -t326 / 0.2e1;
t716 = Ifges(4,1) * t630 + Ifges(4,5) * t628;
t715 = Ifges(5,1) * t655 + Ifges(5,5) * t635;
t669 = t715 + t728;
t714 = t137 * mrSges(5,2) + t669 + t715;
t713 = t60 * mrSges(6,1) + t18 * mrSges(7,1) - t61 * mrSges(6,2) - t19 * mrSges(7,2);
t631 = -t325 / 0.2e1;
t633 = -t323 / 0.2e1;
t637 = -t280 / 0.2e1;
t639 = -t279 / 0.2e1;
t712 = Ifges(5,5) * t637 - Ifges(4,2) * t631 - Ifges(4,6) * t628 + Ifges(5,6) * t639 + Ifges(5,3) * t633;
t711 = t205 * mrSges(5,1) - t112 * mrSges(5,3) + t713;
t708 = t704 / 0.2e1 - t49 * mrSges(5,3) - t744 - t734;
t70 = Ifges(7,2) * t516 + Ifges(7,6) * t272 + t608;
t707 = t70 / 0.2e1;
t706 = Ifges(4,3) * t730;
t705 = -Ifges(3,6) * t446 / 0.2e1;
t230 = t467 * t306 - t308 * t462;
t186 = pkin(5) * t414 + pkin(13) * t489 + t230;
t207 = pkin(13) * t362 + t231;
t106 = t186 * t466 - t207 * t461;
t702 = qJD(6) * t106 + t461 * t739 - t738 * t466;
t107 = t186 * t461 + t207 * t466;
t701 = -qJD(6) * t107 + t738 * t461 + t466 * t739;
t440 = t667 * t462;
t441 = t667 * t467;
t373 = t440 * t466 + t441 * t461;
t695 = qJD(6) * t373 + t461 * t736 + t466 * t737;
t374 = t440 * t461 - t441 * t466;
t694 = -qJD(6) * t374 - t461 * t737 + t466 * t736;
t82 = -mrSges(7,1) * t516 + mrSges(7,2) * t142;
t693 = m(7) * t83 + t82;
t493 = t461 * t462 - t466 * t467;
t412 = t493 * t463;
t355 = pkin(2) * t460 + t453 - t491;
t293 = -t355 * t457 + t459 * t388;
t347 = -t458 * t488 - t460 * t565;
t348 = t460 * t566 + t478;
t220 = pkin(3) * t347 - pkin(11) * t348 + t293;
t540 = pkin(9) * t563 + t452;
t341 = (t457 * t460 + t459 * t563) * pkin(10) + t540;
t247 = t469 * t341 + t355 * t562 + t388 * t566;
t413 = -t457 * t563 + t460 * t459;
t227 = pkin(11) * t413 + t247;
t128 = t463 * t220 + t468 * t227;
t125 = pkin(12) * t347 + t128;
t246 = -t464 * t341 + t469 * (t355 * t459 + t388 * t457);
t226 = -pkin(3) * t413 - t246;
t304 = t348 * t463 - t413 * t468;
t305 = t348 * t468 + t413 * t463;
t149 = pkin(4) * t304 - pkin(12) * t305 + t226;
t76 = t467 * t125 + t462 * t149;
t690 = -pkin(5) * t544 + t689;
t685 = -(t357 * t459 + t392 * t457) * t469 + t750;
t50 = -t463 * t136 + t168 * t468 - t201 * t535 - t206 * t534;
t683 = -t463 * t50 + t468 * t49;
t682 = t12 * t467 - t13 * t462;
t681 = qJD(5) + qJD(6);
t679 = Ifges(4,1) * t732 + Ifges(4,4) * t731;
t678 = -t224 * mrSges(4,3) + t727;
t677 = t225 * mrSges(4,3) + t719;
t676 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t653;
t675 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t653;
t41 = t99 * Ifges(6,4) + t100 * Ifges(6,2) + t192 * Ifges(6,6);
t672 = t41 / 0.2e1;
t671 = Ifges(6,1) * t668 + Ifges(6,4) * t663 + Ifges(6,5) * t653;
t665 = pkin(1) * mrSges(3,1);
t664 = pkin(1) * mrSges(3,2);
t662 = -t114 / 0.2e1;
t660 = t516 / 0.2e1;
t658 = t142 / 0.2e1;
t585 = t323 * Ifges(5,3);
t591 = t280 * Ifges(5,5);
t593 = t279 * Ifges(5,6);
t178 = t585 + t591 + t593;
t657 = t178 / 0.2e1;
t586 = t323 * Ifges(5,6);
t594 = t279 * Ifges(5,2);
t614 = Ifges(5,4) * t280;
t179 = t586 + t594 + t614;
t656 = t179 / 0.2e1;
t276 = Ifges(4,6) * t282;
t277 = Ifges(4,5) * t281;
t195 = Ifges(4,3) * t492 - t276 + t277;
t652 = t195 / 0.2e1;
t651 = Ifges(4,5) * t492 / 0.2e1 + t679;
t649 = -t232 / 0.2e1;
t647 = -t233 / 0.2e1;
t578 = t382 * Ifges(4,6);
t581 = t326 * Ifges(4,4);
t584 = t325 * Ifges(4,2);
t243 = t578 + t581 + t584;
t645 = -t243 / 0.2e1;
t322 = Ifges(4,4) * t325;
t579 = t382 * Ifges(4,5);
t582 = t326 * Ifges(4,1);
t244 = t322 + t579 + t582;
t644 = t244 / 0.2e1;
t642 = t272 / 0.2e1;
t641 = -t274 / 0.2e1;
t638 = t279 / 0.2e1;
t636 = t280 / 0.2e1;
t632 = t323 / 0.2e1;
t629 = t326 / 0.2e1;
t626 = t460 / 0.2e1;
t621 = pkin(5) * t462;
t618 = qJD(2) / 0.2e1;
t617 = mrSges(4,3) * t325;
t616 = mrSges(4,3) * t326;
t615 = Ifges(3,4) * t465;
t613 = Ifges(5,4) * t463;
t612 = Ifges(5,4) * t468;
t190 = Ifges(5,5) * t191;
t189 = Ifges(5,6) * t192;
t600 = t147 * mrSges(4,2);
t148 = (t340 * t459 + t383 * t457) * t469 + t514;
t599 = t148 * mrSges(4,1);
t589 = t281 * Ifges(4,4);
t45 = -pkin(4) * t282 - t50;
t575 = t45 * t463;
t564 = t458 * t465;
t288 = t362 * t466 + t461 * t489;
t154 = qJD(6) * t288 + t267 * t466 + t268 * t461;
t199 = t269 * t461 + t270 * t466;
t552 = t154 - t199;
t289 = t362 * t461 - t466 * t489;
t155 = -qJD(6) * t289 - t267 * t461 + t268 * t466;
t198 = t269 * t466 - t270 * t461;
t551 = t155 - t198;
t183 = t256 * t466 - t257 * t461;
t430 = t461 * t467 + t462 * t466;
t292 = t412 * t681 - t430 * t534;
t550 = t183 - t292;
t184 = t256 * t461 + t257 * t466;
t353 = t681 * t430;
t291 = -t353 * t463 - t493 * t534;
t549 = t184 - t291;
t203 = t430 * t279;
t548 = -t203 + t353;
t204 = t493 * t279;
t352 = t681 * t493;
t547 = -t204 + t352;
t146 = -mrSges(6,1) * t232 + mrSges(6,2) * t233;
t236 = mrSges(5,1) * t323 - mrSges(5,3) * t280;
t546 = t236 - t146;
t91 = Ifges(5,3) * t282 - t189 + t190;
t524 = t465 * t538;
t75 = -t125 * t462 + t467 * t149;
t127 = t220 * t468 - t463 * t227;
t300 = -t359 * t457 + t459 * t393;
t311 = -t463 * t400 + t401 * t468;
t406 = t540 * qJD(1);
t515 = t705 - (Ifges(3,2) * t470 + t615) * t539 / 0.2e1 - t406 * mrSges(3,3);
t444 = qJD(2) * t453;
t358 = t444 - t482;
t513 = -t341 * t536 - t355 * t521 - t464 * t358 - t388 * t523;
t511 = t457 * t524;
t508 = -t50 * mrSges(5,1) + t49 * mrSges(5,2);
t307 = pkin(4) * t565 - t311;
t506 = t599 - t600;
t505 = mrSges(6,1) * t467 - mrSges(6,2) * t462;
t502 = Ifges(6,1) * t462 + t609;
t500 = Ifges(6,2) * t467 + t610;
t498 = Ifges(6,5) * t462 + Ifges(6,6) * t467;
t253 = t305 * t467 + t347 * t462;
t57 = pkin(5) * t304 - pkin(13) * t253 + t75;
t252 = -t305 * t462 + t347 * t467;
t62 = pkin(13) * t252 + t76;
t22 = -t461 * t62 + t466 * t57;
t23 = t461 * t57 + t466 * t62;
t496 = t111 * t468 + t112 * t463;
t174 = t252 * t466 - t253 * t461;
t175 = t252 * t461 + t253 * t466;
t169 = -t341 * t537 + t355 * t520 + t469 * t358 + t359 * t562 + t388 * t522 + t393 * t566;
t158 = pkin(11) * t511 + t169;
t297 = t460 * t523 + t474;
t298 = t460 * t522 + t475;
t182 = pkin(3) * t297 - pkin(11) * t298 + t300;
t59 = -t463 * t158 + t182 * t468 - t220 * t535 - t227 * t534;
t124 = -pkin(4) * t347 - t127;
t58 = t468 * t158 + t463 * t182 + t220 * t534 - t227 * t535;
t54 = pkin(12) * t297 + t58;
t159 = -t359 * t561 + (-pkin(3) * t524 - t393 * t469) * t457 - t513;
t213 = qJD(4) * t305 + t298 * t463 - t468 * t511;
t214 = -qJD(4) * t304 + t298 * t468 + t463 * t511;
t81 = pkin(4) * t213 - pkin(12) * t214 + t159;
t16 = -t125 * t533 + t149 * t532 + t462 * t81 + t467 * t54;
t404 = -pkin(9) * t526 + t443;
t442 = Ifges(3,4) * t525;
t481 = t446 * Ifges(3,5) + Ifges(3,1) * t526 / 0.2e1 + t442 / 0.2e1 - t404 * mrSges(3,3);
t55 = -pkin(4) * t297 - t59;
t410 = t540 * qJD(2);
t17 = -qJD(5) * t76 - t462 * t54 + t467 * t81;
t476 = t224 * mrSges(4,1) - t225 * mrSges(4,2) + t326 * Ifges(4,5) + t325 * Ifges(4,6) + t706;
t472 = -t595 / 0.2e1 - t597 / 0.2e1 - t598 / 0.2e1 + t586 / 0.2e1 + t656 - t596 / 0.2e1 - t601 / 0.2e1 - t602 / 0.2e1 - t69 / 0.2e1 - t113 / 0.2e1 + t614 / 0.2e1 - t711;
t456 = -pkin(5) * t467 - pkin(4);
t436 = Ifges(3,5) * t470 * t517;
t435 = (pkin(11) + t621) * t463;
t419 = -pkin(9) * t564 + t453;
t411 = t430 * t463;
t409 = -pkin(9) * t524 + t444;
t403 = -t446 * mrSges(3,2) + mrSges(3,3) * t525;
t402 = mrSges(3,1) * t446 - mrSges(3,3) * t526;
t397 = qJD(1) * t410;
t396 = -pkin(9) * t510 + t438;
t394 = -pkin(11) * t559 + t427;
t286 = mrSges(4,1) * t382 - t616;
t285 = -mrSges(4,2) * t382 + t617;
t264 = -mrSges(4,1) * t325 + mrSges(4,2) * t326;
t260 = -mrSges(4,2) * t492 - mrSges(4,3) * t282;
t259 = mrSges(4,1) * t492 - mrSges(4,3) * t281;
t258 = -pkin(5) * t362 + t307;
t235 = -mrSges(5,2) * t323 + mrSges(5,3) * t279;
t212 = mrSges(4,1) * t282 + mrSges(4,2) * t281;
t210 = -mrSges(5,1) * t279 + mrSges(5,2) * t280;
t196 = -t282 * Ifges(4,2) + Ifges(4,6) * t492 + t589;
t177 = mrSges(6,1) * t274 - mrSges(6,3) * t233;
t176 = -mrSges(6,2) * t274 + mrSges(6,3) * t232;
t170 = (t359 * t459 + t393 * t457) * t469 + t513;
t151 = -mrSges(5,2) * t282 - mrSges(5,3) * t192;
t150 = mrSges(5,1) * t282 - mrSges(5,3) * t191;
t119 = qJD(5) * t252 + t214 * t467 + t297 * t462;
t118 = -qJD(5) * t253 - t214 * t462 + t297 * t467;
t109 = mrSges(7,1) * t272 - mrSges(7,3) * t142;
t108 = -mrSges(7,2) * t272 + mrSges(7,3) * t516;
t105 = mrSges(5,1) * t192 + mrSges(5,2) * t191;
t96 = -pkin(5) * t252 + t124;
t95 = pkin(5) * t569 + t112;
t78 = -mrSges(6,2) * t192 + mrSges(6,3) * t100;
t77 = mrSges(6,1) * t192 - mrSges(6,3) * t99;
t56 = -mrSges(6,1) * t100 + mrSges(6,2) * t99;
t48 = -qJD(6) * t175 + t118 * t466 - t119 * t461;
t47 = qJD(6) * t174 + t118 * t461 + t119 * t466;
t35 = -pkin(5) * t118 + t55;
t32 = -mrSges(7,2) * t192 + mrSges(7,3) * t39;
t31 = mrSges(7,1) * t192 - mrSges(7,3) * t38;
t24 = -pkin(5) * t100 + t45;
t21 = t466 * t51 - t574;
t20 = -t461 * t51 - t572;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t14 = pkin(13) * t118 + t16;
t8 = pkin(5) * t213 - pkin(13) * t119 + t17;
t5 = -qJD(6) * t23 - t14 * t461 + t466 * t8;
t4 = qJD(6) * t22 + t14 * t466 + t461 * t8;
t1 = [(Ifges(7,4) * t175 + Ifges(7,2) * t174) * t673 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t660 + (Ifges(6,5) * t119 + Ifges(6,6) * t118) * t640 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t642 + (-t179 / 0.2e1 + Ifges(7,6) * t660 + Ifges(6,6) * t648 + Ifges(7,5) * t658 - Ifges(5,2) * t638 + Ifges(6,3) * t640 + Ifges(7,3) * t642 + Ifges(6,5) * t646 - Ifges(5,6) * t632 - Ifges(5,4) * t636 + t696 / 0.2e1 + t711) * t213 + m(3) * (t396 * t540 - t397 * t419 - t404 * t410 + t406 * t409) + (t118 * t61 - t119 * t60 + t12 * t252 - t13 * t253) * mrSges(6,3) + (Ifges(5,4) * t305 + Ifges(5,6) * t347) * t654 + (Ifges(5,4) * t214 + Ifges(5,6) * t297) * t638 + (Ifges(6,1) * t253 + Ifges(6,4) * t252) * t668 + (Ifges(6,1) * t119 + Ifges(6,4) * t118) * t646 + (-t112 * t297 + t137 * t305 + t205 * t214 - t347 * t49) * mrSges(5,2) + (Ifges(6,4) * t253 + Ifges(6,2) * t252) * t663 + (Ifges(6,4) * t119 + Ifges(6,2) * t118) * t648 + (-t147 * t347 - t148 * t348 - t224 * t298 - t225 * t297) * mrSges(4,3) + (Ifges(7,1) * t175 + Ifges(7,4) * t174) * t674 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t658 - t413 * t600 + ((Ifges(3,5) * t626 - t419 * mrSges(3,3) + (-0.2e1 * t664 + 0.3e1 / 0.2e1 * Ifges(3,4) * t470) * t458) * t470 + (-Ifges(3,6) * t460 - t540 * mrSges(3,3) + t457 * (Ifges(4,5) * t348 - Ifges(4,6) * t347 + Ifges(4,3) * t413) / 0.2e1 + (-0.2e1 * t665 - 0.3e1 / 0.2e1 * t615 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t470) * t458) * t465) * t517 + t396 * (-t460 * mrSges(3,2) + mrSges(3,3) * t563) - t397 * (mrSges(3,1) * t460 - mrSges(3,3) * t564) + m(4) * (t147 * t247 + t148 * t246 + t169 * t225 + t170 * t224 + t271 * t300 + t287 * t293) + m(5) * (t111 * t59 + t112 * t58 + t127 * t50 + t128 * t49 + t137 * t226 + t159 * t205) + m(6) * (t101 * t55 + t12 * t76 + t124 * t45 + t13 * t75 + t16 * t61 + t17 * t60) + m(7) * (t18 * t5 + t19 * t4 + t2 * t23 + t22 * t3 + t24 * t96 + t35 * t83) + (t174 * t2 - t175 * t3 - t18 * t47 + t19 * t48) * mrSges(7,3) + t325 * (Ifges(4,4) * t298 - Ifges(4,2) * t297) / 0.2e1 + t708 * t304 + (Ifges(5,1) * t305 + Ifges(5,5) * t347) * t655 + (Ifges(5,1) * t214 + Ifges(5,5) * t297) * t636 + (Ifges(6,5) * t253 + Ifges(7,5) * t175 + Ifges(6,6) * t252 + Ifges(7,6) * t174) * t653 + (Ifges(4,5) * t298 - Ifges(4,6) * t297) * t730 + (Ifges(4,4) * t348 - Ifges(4,2) * t347 + Ifges(4,6) * t413) * t731 + (Ifges(4,1) * t348 - Ifges(4,4) * t347 + Ifges(4,5) * t413) * t732 + t409 * t403 - t410 * t402 - t347 * t196 / 0.2e1 + t287 * (mrSges(4,1) * t347 + mrSges(4,2) * t348) + t50 * (mrSges(5,1) * t347 - mrSges(5,3) * t305) + t347 * t91 / 0.2e1 + t111 * (mrSges(5,1) * t297 - mrSges(5,3) * t214) + t271 * (mrSges(4,1) * t297 + mrSges(4,2) * t298) + t300 * t264 + t293 * t212 + t169 * t285 + t170 * t286 + t246 * t259 + t247 * t260 + t45 * (-mrSges(6,1) * t252 + mrSges(6,2) * t253) + t59 * t236 + t58 * t235 + t226 * t105 + (t481 * t470 + (t705 + (t706 + t476) * t457 + t515) * t465) * t538 + t214 * t180 / 0.2e1 + t159 * t210 + t305 * t669 + t253 * t671 + t252 * t672 + t175 * t675 + t174 * t676 + t348 * t651 + t413 * t652 + t297 * t657 + t298 * t644 + t297 * t645 + (Ifges(4,1) * t298 - Ifges(4,4) * t297) * t629 + t436 * t626 + t413 * t599 + t17 * t177 + t24 * (-mrSges(7,1) * t174 + mrSges(7,2) * t175) + t16 * t176 + t55 * t146 + t127 * t150 + t128 * t151 + (Ifges(5,5) * t214 + Ifges(5,3) * t297) * t632 + (Ifges(5,5) * t305 + Ifges(5,3) * t347) * t635 + t22 * t31 + t23 * t32 + t48 * t707 + t47 * t71 / 0.2e1 + t75 * t77 + t76 * t78 + t35 * t82 + t83 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + t96 * t15 + t4 * t108 + t5 * t109 + t118 * t114 / 0.2e1 + t119 * t115 / 0.2e1 + t101 * (-mrSges(6,1) * t118 + mrSges(6,2) * t119) + t124 * t56; -t711 * t742 + t717 * t751 + (Ifges(5,4) * t318 - Ifges(5,2) * t317) * t639 + (t243 / 0.2e1 - t178 / 0.2e1 - Ifges(4,4) * t630 + t677 + t712) * t377 + (Ifges(7,1) * t289 + Ifges(7,4) * t288) * t674 + (Ifges(5,5) * t318 - Ifges(5,6) * t317) * t633 + (-mrSges(6,1) * t544 + mrSges(6,2) * t545) * t101 + (t277 / 0.2e1 - t276 / 0.2e1 + t652 + t506) * t459 + (-t18 * t552 + t19 * t551 + t2 * t288 - t289 * t3) * mrSges(7,3) + (-t198 / 0.2e1 + t155 / 0.2e1) * t70 + (t360 / 0.2e1 - t318 / 0.2e1) * t180 + (-mrSges(7,1) * t551 + mrSges(7,2) * t552) * t83 + (-t269 / 0.2e1 + t268 / 0.2e1) * t114 + (Ifges(5,1) * t318 - Ifges(5,4) * t317) * t637 + ((-t442 / 0.2e1 + t539 * t664 - t481) * t470 + ((t665 + t615 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t470) * t539 + (t446 / 0.2e1 - qJD(2)) * Ifges(3,6) + ((Ifges(4,5) * t464 + Ifges(4,6) * t469) * t457 * t618 + (t459 * t618 + t628) * Ifges(4,3) - t476) * t457 - t515) * t465) * t539 + (Ifges(7,4) * t289 + Ifges(7,2) * t288) * t673 + (-Ifges(6,1) * t489 + Ifges(6,4) * t362) * t668 + (t12 * t362 + t13 * t489 + t544 * t61 - t545 * t60) * mrSges(6,3) + t45 * (-mrSges(6,1) * t362 - mrSges(6,2) * t489) - t489 * t671 + (-Ifges(6,4) * t489 + Ifges(6,2) * t362) * t663 + (-Ifges(6,5) * t489 + Ifges(7,5) * t289 + Ifges(6,6) * t362 + Ifges(7,6) * t288) * t653 + t436 + (-t270 / 0.2e1 + t267 / 0.2e1) * t115 + t708 * t414 + (-t199 / 0.2e1 + t154 / 0.2e1) * t71 + t418 * t259 + t420 * t260 - t404 * t403 + t406 * t402 - t396 * mrSges(3,2) - t397 * mrSges(3,1) + t399 * t105 + t307 * t56 + t311 * t150 + t312 * t151 - t299 * t264 + t24 * (-mrSges(7,1) * t288 + mrSges(7,2) * t289) + t258 * t15 + t230 * t77 + t231 * t78 + t701 * t109 + t702 * t108 + (t106 * t3 + t107 * t2 + t18 * t701 + t19 * t702 + t24 * t258 + t690 * t83) * m(7) + t362 * t672 + t289 * t675 + t288 * t676 + (Ifges(7,1) * t199 + Ifges(7,4) * t198 + Ifges(7,5) * t317) * t659 + (Ifges(7,4) * t154 + Ifges(7,2) * t155 + Ifges(7,6) * t361) * t660 + (Ifges(7,4) * t199 + Ifges(7,2) * t198 + Ifges(7,6) * t317) * t661 + (Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t317) * t649 + (Ifges(7,1) * t154 + Ifges(7,4) * t155 + Ifges(7,5) * t361) * t658 + (Ifges(5,4) * t360 - Ifges(5,2) * t361) * t638 + (Ifges(6,5) * t267 + Ifges(6,6) * t268 + Ifges(6,3) * t361) * t640 + (Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t317) * t641 + (Ifges(7,5) * t154 + Ifges(7,6) * t155 + Ifges(7,3) * t361) * t642 + (Ifges(7,5) * t199 + Ifges(7,6) * t198 + Ifges(7,3) * t317) * t643 + (Ifges(6,1) * t267 + Ifges(6,4) * t268 + Ifges(6,5) * t361) * t646 + (Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t317) * t647 + (Ifges(6,4) * t267 + Ifges(6,2) * t268 + Ifges(6,6) * t361) * t648 + (Ifges(5,5) * t360 - Ifges(5,6) * t361) * t632 + (Ifges(5,1) * t360 - Ifges(5,4) * t361) * t636 + (-t179 + t696) * (t361 / 0.2e1 - t317 / 0.2e1) + t691 * t176 + t692 * t177 + (t101 * t689 + t12 * t231 + t13 * t230 + t307 * t45 + t60 * t692 + t61 * t691) * m(6) + (-t244 / 0.2e1 + Ifges(4,4) * t631 - t678 + t716) * t378 + (-mrSges(5,3) * t50 + t714 + t728) * t415 + (-pkin(2) * t212 + (t287 * mrSges(4,2) - t148 * mrSges(4,3) + t651 + t679) * t464 + (t147 * mrSges(4,3) + t189 / 0.2e1 - t190 / 0.2e1 + t589 / 0.2e1 - t287 * mrSges(4,1) + t196 / 0.2e1 - t91 / 0.2e1 + (-Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1) * t282 + t508) * t469 + ((t644 + t579 / 0.2e1 + t582 / 0.2e1 + t322 / 0.2e1 + t678) * t469 + (t585 / 0.2e1 + t593 / 0.2e1 + t591 / 0.2e1 + t645 + t657 - t578 / 0.2e1 - t581 / 0.2e1 - t584 / 0.2e1 - t677) * t464) * qJD(3)) * t457 + t106 * t31 + t107 * t32 + t684 * t210 + t685 * t286 + t686 * t285 + (-pkin(2) * t287 * t457 + t147 * t420 + t148 * t418 + t224 * t685 + t225 * t686 - t271 * t299) * m(4) + t687 * t236 + t688 * t235 + (t111 * t687 + t112 * t688 + t137 * t399 + t205 * t684 + t311 * t50 + t312 * t49) * m(5) + t689 * t146 + t690 * t82; ((t498 * t641 + t502 * t647 + t500 * t649 + t101 * t505 + t467 * t662 + t115 * t625 + (t462 * t60 - t467 * t61) * mrSges(6,3)) * qJD(5) + t501 * t663 + t503 * t668 + t499 * t653 + (-t150 + t56) * pkin(11) + (-t594 / 0.2e1 - t472 - t235 * pkin(11)) * qJD(4) + t714) * t463 + (-t581 + t178) * t630 + (t24 * t411 + t550 * t83) * mrSges(7,1) + (-Ifges(7,1) * t412 - Ifges(7,4) * t411) * t674 + (-Ifges(7,4) * t412 - Ifges(7,2) * t411) * t673 + (-t210 + t286 + t616) * t225 + (-t285 + t617) * t224 + (Ifges(7,1) * t184 + Ifges(7,4) * t183) * t659 + (Ifges(7,4) * t184 + Ifges(7,2) * t183) * t661 + (Ifges(6,5) * t257 + Ifges(6,6) * t256) * t641 + (t292 / 0.2e1 - t183 / 0.2e1) * t70 + (-t12 * t560 - t13 * t558 - t256 * t61 + t257 * t60) * mrSges(6,3) + (t151 * pkin(11) + (t592 / 0.2e1 + (m(6) * t101 - t546) * pkin(11) - t721) * qJD(4) + t734) * t468 + (-pkin(3) * t137 - t111 * t162 - t112 * t163 - t205 * t225 + (-qJD(4) * t496 + t683) * pkin(11)) * m(5) + (Ifges(6,4) * t257 + Ifges(6,2) * t256) * t649 + (t322 + t244) * t631 + (Ifges(6,5) * t647 + Ifges(7,5) * t659 + Ifges(6,6) * t649 + Ifges(7,6) * t661 + Ifges(6,3) * t641 + Ifges(7,3) * t643 + t656 - t713 - t696 / 0.2e1) * t568 + t506 + t612 * t655 + (-t24 * t412 - t549 * t83) * mrSges(7,2) + (-Ifges(7,5) * t412 - Ifges(7,6) * t411) * t653 + (t291 / 0.2e1 - t184 / 0.2e1) * t71 - t41 * t560 / 0.2e1 + (t18 * t549 - t19 * t550 - t2 * t411 + t3 * t412) * mrSges(7,3) + (Ifges(6,1) * t257 + Ifges(6,4) * t256) * t647 + t435 * t15 + t394 * t77 + t395 * t78 + t283 * t31 + t284 * t32 - t101 * (-mrSges(6,1) * t256 + mrSges(6,2) * t257) - t257 * t115 / 0.2e1 - t162 * t236 - t163 * t235 + t558 * t671 - t412 * t675 - t411 * t676 + (Ifges(7,4) * t291 + Ifges(7,2) * t292) * t660 + t256 * t662 + (Ifges(7,1) * t291 + Ifges(7,4) * t292) * t658 + (Ifges(7,5) * t291 + Ifges(7,6) * t292) * t642 + t243 * t629 + t504 * t575 - t134 * t146 + t195 + (Ifges(7,5) * t184 + Ifges(7,6) * t183) * t643 + t683 * mrSges(5,3) + (t712 + t719) * t326 + t722 * t82 + t723 * t177 + t724 * t176 + (pkin(11) * t575 - t101 * t134 + t12 * t395 + t13 * t394 + t60 * t723 + t61 * t724) * m(6) + t725 * t109 + t726 * t108 + (t18 * t725 + t19 * t726 + t2 * t284 + t24 * t435 + t283 * t3 + t722 * t83) * m(7) + (-t205 * (mrSges(5,1) * t463 + mrSges(5,2) * t468) - t727 + (Ifges(5,1) * t468 - t613) * t637 + (-Ifges(5,2) * t463 + t612) * t639 + (Ifges(5,5) * t468 - Ifges(5,6) * t463) * t633 + t496 * mrSges(5,3) + t716) * t325 + t613 * t654 - (t180 * t325 + t704) * t468 / 0.2e1 - pkin(3) * t105; (-pkin(4) * t45 - t101 * t112 - t60 * t86 - t61 * t87) * m(6) + (-Ifges(7,1) * t352 - Ifges(7,4) * t353) * t658 + (-Ifges(7,5) * t352 - Ifges(7,6) * t353) * t642 + (-t353 / 0.2e1 + t203 / 0.2e1) * t70 + t546 * t112 + (mrSges(7,1) * t548 - mrSges(7,2) * t547) * t83 + t91 - t45 * t505 + (t621 * t693 - t735) * qJD(5) + t472 * t280 + (t18 * t547 - t19 * t548 - t2 * t493 - t3 * t430) * mrSges(7,3) + (Ifges(7,5) * t430 - Ifges(7,6) * t493 + t498) * t653 + t24 * (mrSges(7,1) * t493 + mrSges(7,2) * t430) + (Ifges(7,4) * t430 - Ifges(7,2) * t493) * t673 + (Ifges(7,1) * t430 - Ifges(7,4) * t493) * t674 - t493 * t676 + (-t352 / 0.2e1 + t204 / 0.2e1) * t71 + (-Ifges(7,4) * t352 - Ifges(7,2) * t353) * t660 + (-Ifges(7,1) * t204 - Ifges(7,4) * t203) * t659 + (-Ifges(7,4) * t204 - Ifges(7,2) * t203) * t661 + (-Ifges(7,5) * t204 - Ifges(7,6) * t203) * t643 - t508 + t456 * t15 + t373 * t31 + t374 * t32 - t111 * t235 + t462 * t671 + t430 * t675 + t500 * t663 + t502 * t668 + t41 * t624 - t86 * t177 - t87 * t176 + t695 * t108 + (t18 * t694 + t19 * t695 + t2 * t374 + t24 * t456 + t3 * t373 - t83 * t95) * m(7) + t694 * t109 + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t280 + t721) * t279 - pkin(4) * t56 - t95 * t82 + (-t462 * t77 + t467 * t78 + m(6) * t682 + (-m(6) * t497 - t462 * t176 - t467 * t177) * qJD(5)) * pkin(12) + t682 * mrSges(6,3); -m(7) * (t18 * t20 + t19 * t21) + (-Ifges(6,2) * t233 + t115 + t229) * t649 + t142 * t707 + t40 - t101 * (mrSges(6,1) * t233 + mrSges(6,2) * t232) + (Ifges(6,5) * t232 - Ifges(6,6) * t233) * t641 + t114 * t646 + (Ifges(6,1) * t232 - t611) * t647 + t61 * t177 - t60 * t176 + (t466 * t31 + t461 * t32 + m(7) * (t2 * t461 + t3 * t466) - t693 * t233 + (t108 * t466 - t109 * t461 + m(7) * (-t18 * t461 + t19 * t466)) * qJD(6)) * pkin(5) + (t232 * t60 + t233 * t61) * mrSges(6,3) - t21 * t108 - t20 * t109 + t718 + t733; -t18 * t108 + t19 * t109 + t70 * t658 + t733;];
tauc  = t1(:);
