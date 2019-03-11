% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:17
% EndTime: 2019-03-10 03:34:07
% DurationCPUTime: 28.79s
% Computational Cost: add. (35536->912), mult. (80707->1186), div. (0->0), fcn. (61383->18), ass. (0->465)
t440 = sin(qJ(5));
t539 = qJD(5) * t440;
t442 = sin(qJ(3));
t443 = sin(qJ(2));
t448 = cos(qJ(3));
t449 = cos(qJ(2));
t366 = -t442 * t443 + t448 * t449;
t344 = t366 * qJD(1);
t368 = t442 * t449 + t443 * t448;
t345 = t368 * qJD(1);
t441 = sin(qJ(4));
t447 = cos(qJ(4));
t692 = t447 * t344 - t441 * t345;
t734 = t692 * t440;
t757 = t539 - t734;
t756 = -mrSges(6,3) - mrSges(7,3);
t446 = cos(qJ(5));
t540 = qJD(4) * t447;
t528 = pkin(3) * t540;
t452 = -pkin(8) - pkin(7);
t400 = t452 * t449;
t376 = qJD(1) * t400;
t349 = t448 * t376;
t398 = t452 * t443;
t375 = qJD(1) * t398;
t353 = qJD(2) * pkin(2) + t375;
t286 = t353 * t442 - t349;
t620 = pkin(9) * t344;
t251 = t286 + t620;
t240 = t441 * t251;
t346 = t442 * t376;
t285 = t448 * t353 + t346;
t338 = t345 * pkin(9);
t250 = t285 - t338;
t167 = t250 * t447 - t240;
t482 = t344 * t441 + t447 * t345;
t218 = pkin(4) * t482 - pkin(10) * t692;
t630 = pkin(3) * t345;
t189 = t218 + t630;
t99 = t446 * t167 + t440 * t189;
t755 = t446 * t528 - t99;
t98 = -t167 * t440 + t446 * t189;
t754 = -t440 * t528 - t98;
t438 = qJ(2) + qJ(3);
t431 = qJ(4) + t438;
t412 = sin(t431);
t437 = qJ(5) + qJ(6);
t429 = cos(t437);
t611 = mrSges(7,1) * t429;
t612 = mrSges(6,1) * t446;
t753 = (t611 + t612) * t412;
t427 = sin(t437);
t608 = mrSges(7,2) * t427;
t609 = mrSges(6,2) * t440;
t752 = (-t608 - t609) * t412;
t292 = -t375 * t442 + t349;
t258 = t292 - t620;
t293 = t448 * t375 + t346;
t259 = -t338 + t293;
t180 = t258 * t441 + t259 * t447;
t546 = qJD(1) * t443;
t421 = pkin(2) * t546;
t181 = t189 + t421;
t101 = t446 * t180 + t440 * t181;
t631 = pkin(2) * t448;
t417 = pkin(3) + t631;
t541 = qJD(4) * t441;
t559 = t441 * t442;
t279 = t417 * t540 + (-t442 * t541 + (t447 * t448 - t559) * qJD(3)) * pkin(2);
t751 = t279 * t446 - t101;
t100 = -t180 * t440 + t446 * t181;
t750 = -t279 * t440 - t100;
t535 = pkin(11) * t734;
t558 = t442 * t447;
t340 = pkin(2) * t558 + t441 * t417;
t335 = pkin(10) + t340;
t614 = -pkin(11) - t335;
t509 = qJD(5) * t614;
t749 = t440 * t509 + t535 + t751;
t733 = t692 * t446;
t497 = pkin(5) * t482 - pkin(11) * t733;
t748 = t446 * t509 - t497 + t750;
t628 = pkin(3) * t441;
t414 = pkin(10) + t628;
t613 = -pkin(11) - t414;
t508 = qJD(5) * t613;
t747 = t440 * t508 + t535 + t755;
t746 = t446 * t508 - t497 + t754;
t435 = qJD(2) + qJD(3);
t238 = pkin(3) * t435 + t250;
t161 = t238 * t447 - t240;
t109 = t446 * t161 + t440 * t218;
t451 = -pkin(11) - pkin(10);
t521 = qJD(5) * t451;
t745 = t440 * t521 - t109 + t535;
t108 = -t161 * t440 + t446 * t218;
t744 = t446 * t521 - t108 - t497;
t413 = cos(t431);
t623 = pkin(5) * t446;
t415 = pkin(4) + t623;
t479 = -t412 * t415 - t413 * t451;
t743 = -m(7) * t479 + t753;
t426 = qJD(4) + t435;
t636 = -t426 / 0.2e1;
t642 = -t692 / 0.2e1;
t273 = qJD(5) - t692;
t644 = -t273 / 0.2e1;
t265 = qJD(6) + t273;
t646 = -t265 / 0.2e1;
t253 = t426 * t440 + t446 * t482;
t648 = -t253 / 0.2e1;
t252 = t426 * t446 - t440 * t482;
t649 = -t252 / 0.2e1;
t439 = sin(qJ(6));
t445 = cos(qJ(6));
t172 = t252 * t439 + t253 * t445;
t652 = -t172 / 0.2e1;
t504 = t445 * t252 - t253 * t439;
t654 = -t504 / 0.2e1;
t241 = t447 * t251;
t162 = t238 * t441 + t241;
t157 = pkin(10) * t426 + t162;
t433 = t449 * pkin(2);
t418 = t433 + pkin(1);
t396 = t418 * qJD(1);
t301 = -pkin(3) * t344 - t396;
t177 = -pkin(4) * t692 - pkin(10) * t482 + t301;
t93 = t157 * t446 + t177 * t440;
t75 = pkin(11) * t252 + t93;
t588 = t439 * t75;
t92 = -t157 * t440 + t446 * t177;
t74 = -pkin(11) * t253 + t92;
t68 = pkin(5) * t273 + t74;
t26 = t445 * t68 - t588;
t586 = t445 * t75;
t27 = t439 * t68 + t586;
t670 = t301 * mrSges(5,1) + t92 * mrSges(6,1) + t26 * mrSges(7,1) - t93 * mrSges(6,2) - t27 * mrSges(7,2);
t742 = Ifges(6,5) * t648 + Ifges(7,5) * t652 - Ifges(5,2) * t642 - Ifges(5,6) * t636 + Ifges(6,6) * t649 + Ifges(7,6) * t654 + Ifges(6,3) * t644 + Ifges(7,3) * t646 - t670;
t156 = -pkin(4) * t426 - t161;
t605 = mrSges(5,3) * t482;
t697 = -mrSges(5,1) * t426 - mrSges(6,1) * t252 + mrSges(6,2) * t253 + t605;
t675 = -m(6) * t156 - t697;
t741 = -m(5) * t161 - t675;
t269 = Ifges(5,4) * t692;
t708 = t426 * Ifges(5,5);
t208 = Ifges(5,1) * t482 + t269 + t708;
t679 = t301 * mrSges(5,2) - t161 * mrSges(5,3);
t740 = -t208 / 0.2e1 - t679;
t731 = t757 * pkin(5);
t537 = qJD(1) * qJD(2);
t382 = qJDD(1) * t449 - t443 * t537;
t383 = qJDD(1) * t443 + t449 * t537;
t466 = t366 * qJD(3);
t256 = qJD(1) * t466 + t382 * t442 + t383 * t448;
t467 = t368 * qJD(3);
t257 = -qJD(1) * t467 + t382 * t448 - t383 * t442;
t149 = -t441 * t256 + t257 * t447 - t344 * t541 - t345 * t540;
t147 = qJDD(5) - t149;
t140 = qJDD(6) + t147;
t656 = t140 / 0.2e1;
t148 = qJD(4) * t692 + t256 * t447 + t257 * t441;
t434 = qJDD(2) + qJDD(3);
t425 = qJDD(4) + t434;
t104 = qJD(5) * t252 + t148 * t446 + t425 * t440;
t105 = -qJD(5) * t253 - t148 * t440 + t425 * t446;
t43 = -qJD(6) * t172 - t104 * t439 + t105 * t445;
t663 = t43 / 0.2e1;
t42 = qJD(6) * t504 + t104 * t445 + t105 * t439;
t664 = t42 / 0.2e1;
t666 = Ifges(7,1) * t664 + Ifges(7,4) * t663 + Ifges(7,5) * t656;
t667 = Ifges(7,4) * t664 + Ifges(7,2) * t663 + Ifges(7,6) * t656;
t121 = -pkin(5) * t252 + t156;
t598 = Ifges(6,4) * t253;
t142 = Ifges(6,2) * t252 + Ifges(6,6) * t273 + t598;
t599 = Ifges(5,4) * t482;
t711 = Ifges(5,6) * t426;
t712 = Ifges(5,2) * t692;
t207 = t599 + t711 + t712;
t362 = t383 * pkin(7);
t299 = qJDD(2) * pkin(2) - pkin(8) * t383 - t362;
t361 = t382 * pkin(7);
t300 = pkin(8) * t382 + t361;
t188 = -qJD(3) * t286 + t448 * t299 - t300 * t442;
t126 = pkin(3) * t434 - pkin(9) * t256 + t188;
t542 = qJD(3) * t448;
t543 = qJD(3) * t442;
t187 = t442 * t299 + t448 * t300 + t353 * t542 + t376 * t543;
t137 = pkin(9) * t257 + t187;
t54 = t126 * t447 - t441 * t137 - t238 * t541 - t251 * t540;
t51 = -pkin(4) * t425 - t54;
t25 = -pkin(5) * t105 + t51;
t367 = t439 * t446 + t440 * t445;
t682 = qJD(5) + qJD(6);
t295 = t682 * t367;
t53 = t441 * t126 + t447 * t137 + t238 * t540 - t251 * t541;
t50 = pkin(10) * t425 + t53;
t585 = qJDD(1) * pkin(1);
t336 = -pkin(2) * t382 - t585;
t223 = -pkin(3) * t257 + t336;
t63 = -pkin(4) * t149 - pkin(10) * t148 + t223;
t16 = -qJD(5) * t93 - t440 * t50 + t446 * t63;
t7 = pkin(5) * t147 - pkin(11) * t104 + t16;
t538 = qJD(5) * t446;
t15 = -t157 * t539 + t177 * t538 + t440 * t63 + t446 * t50;
t8 = pkin(11) * t105 + t15;
t3 = qJD(6) * t26 + t439 * t7 + t445 * t8;
t32 = Ifges(6,4) * t104 + Ifges(6,2) * t105 + Ifges(6,6) * t147;
t4 = -qJD(6) * t27 - t439 * t8 + t445 * t7;
t490 = mrSges(6,1) * t440 + mrSges(6,2) * t446;
t471 = t156 * t490;
t478 = t439 * t440 - t445 * t446;
t484 = Ifges(6,5) * t446 - Ifges(6,6) * t440;
t596 = Ifges(6,4) * t446;
t486 = -Ifges(6,2) * t440 + t596;
t597 = Ifges(6,4) * t440;
t488 = Ifges(6,1) * t446 - t597;
t511 = -t539 / 0.2e1;
t246 = Ifges(6,4) * t252;
t143 = Ifges(6,1) * t253 + Ifges(6,5) * t273 + t246;
t556 = t446 * t143;
t515 = t556 / 0.2e1;
t640 = t482 / 0.2e1;
t645 = t265 / 0.2e1;
t651 = t172 / 0.2e1;
t653 = t504 / 0.2e1;
t655 = t147 / 0.2e1;
t657 = t105 / 0.2e1;
t658 = t104 / 0.2e1;
t165 = Ifges(7,4) * t504;
t83 = Ifges(7,1) * t172 + Ifges(7,5) * t265 + t165;
t659 = t83 / 0.2e1;
t660 = -t83 / 0.2e1;
t595 = Ifges(7,4) * t172;
t82 = Ifges(7,2) * t504 + Ifges(7,6) * t265 + t595;
t661 = t82 / 0.2e1;
t662 = -t82 / 0.2e1;
t665 = Ifges(6,1) * t658 + Ifges(6,4) * t657 + Ifges(6,5) * t655;
t709 = Ifges(6,3) * t273;
t710 = Ifges(6,6) * t252;
t703 = Ifges(6,5) * t253 + Ifges(7,5) * t172 + Ifges(7,6) * t504 + Ifges(7,3) * t265 + t709 + t710;
t716 = -t482 / 0.2e1;
t722 = mrSges(7,2) * t121 - mrSges(7,3) * t26;
t723 = -mrSges(7,1) * t121 + mrSges(7,3) * t27;
t725 = t609 - t612;
t726 = t15 * t446 - t16 * t440;
t735 = t478 * t692;
t736 = t367 * t692;
t739 = (t252 * t486 + t253 * t488 + t273 * t484) * qJD(5) / 0.2e1 + (-t599 + t703) * t716 + (t471 + t515) * qJD(5) + t51 * t725 - (Ifges(7,4) * t651 + Ifges(7,2) * t653 + Ifges(7,6) * t645 + t661 + t723) * t295 + (t726 - t757 * t93 + (-t538 + t733) * t92) * mrSges(6,3) - t736 * t662 + (-Ifges(7,5) * t735 - Ifges(7,6) * t736) * t646 + (-Ifges(7,4) * t735 - Ifges(7,2) * t736) * t654 + (-Ifges(7,1) * t735 - Ifges(7,4) * t736) * t652 - t121 * (mrSges(7,1) * t736 - mrSges(7,2) * t735) + (-t26 * t735 + t27 * t736) * mrSges(7,3) - t735 * t660 + t54 * mrSges(5,1) - t53 * mrSges(5,2) + (Ifges(6,2) * t446 + t597) * t657 + (Ifges(6,1) * t440 + t596) * t658 + t440 * t665 + (Ifges(6,5) * t440 + Ifges(6,6) * t446) * t655 + t207 * t640 + t446 * t32 / 0.2e1 + Ifges(5,3) * t425 + Ifges(5,5) * t148 + Ifges(5,6) * t149 + t142 * t511 + (t25 * mrSges(7,1) - 0.2e1 * t667 - (Ifges(7,1) * t651 + Ifges(7,4) * t653 + Ifges(7,5) * t645 + t659 + t722) * t682 - t3 * mrSges(7,3)) * t478 + (mrSges(7,2) * t25 - mrSges(7,3) * t4 + 0.2e1 * t666) * t367;
t668 = m(7) * pkin(5);
t732 = t447 * t258 - t259 * t441 + t417 * t541 + (t442 * t540 + (t441 * t448 + t558) * qJD(3)) * pkin(2);
t428 = sin(t438);
t430 = cos(t438);
t491 = mrSges(5,1) * t412 + mrSges(5,2) * t413;
t730 = mrSges(4,1) * t428 + mrSges(4,2) * t430 + t491;
t166 = t250 * t441 + t241;
t729 = pkin(3) * t541 - t166;
t450 = cos(qJ(1));
t568 = t413 * t450;
t728 = t752 * t450 + t568 * t756;
t444 = sin(qJ(1));
t569 = t413 * t444;
t727 = t752 * t444 + t569 * t756;
t684 = g(1) * t450 + g(2) * t444;
t724 = -t413 * mrSges(5,1) + (mrSges(5,2) + t756) * t412;
t688 = t413 * pkin(4) + t412 * pkin(10);
t691 = -t412 * t451 + t413 * t415;
t720 = -m(6) * t688 - m(7) * t691;
t562 = t440 * t142;
t719 = Ifges(5,1) * t716 + Ifges(5,5) * t636 + t484 * t644 + t486 * t649 + t488 * t648 - t471 - t556 / 0.2e1 + t562 / 0.2e1;
t717 = t382 / 0.2e1;
t635 = t449 / 0.2e1;
t715 = t440 * t668;
t707 = t449 * Ifges(3,2);
t706 = mrSges(6,1) + t668;
t302 = t614 * t440;
t432 = t446 * pkin(11);
t573 = t335 * t446;
t303 = t432 + t573;
t232 = t302 * t445 - t303 * t439;
t705 = qJD(6) * t232 + t439 * t748 + t445 * t749;
t233 = t302 * t439 + t303 * t445;
t704 = -qJD(6) * t233 - t439 * t749 + t445 * t748;
t358 = t613 * t440;
t567 = t414 * t446;
t359 = t432 + t567;
t282 = t358 * t439 + t359 * t445;
t702 = -qJD(6) * t282 - t439 * t747 + t445 * t746;
t281 = t358 * t445 - t359 * t439;
t701 = qJD(6) * t281 + t439 * t746 + t445 * t747;
t397 = t451 * t440;
t619 = pkin(10) * t446;
t399 = t432 + t619;
t306 = t397 * t439 + t399 * t445;
t700 = -qJD(6) * t306 - t439 * t745 + t445 * t744;
t304 = t397 * t445 - t399 * t439;
t699 = qJD(6) * t304 + t439 * t744 + t445 * t745;
t291 = t366 * t441 + t368 * t447;
t220 = t478 * t291;
t698 = -t162 + t731;
t305 = t448 * t398 + t400 * t442;
t270 = -pkin(9) * t368 + t305;
t307 = t442 * t398 - t448 * t400;
t271 = pkin(9) * t366 + t307;
t210 = t270 * t441 + t271 * t447;
t200 = t446 * t210;
t324 = -pkin(3) * t366 - t418;
t481 = t447 * t366 - t368 * t441;
t205 = -pkin(4) * t481 - pkin(10) * t291 + t324;
t114 = t440 * t205 + t200;
t696 = t447 * t270 - t271 * t441;
t695 = t731 + t732;
t694 = t279 - t180;
t690 = t729 + t731;
t506 = t430 * mrSges(4,1) - mrSges(4,2) * t428;
t626 = pkin(4) * t412;
t629 = pkin(3) * t428;
t687 = -m(7) * (t479 - t629) - m(6) * (-t626 - t629) + t753;
t545 = qJD(1) * t449;
t621 = pkin(7) * t449;
t622 = pkin(7) * t443;
t686 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t546) * t621 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t545) * t622;
t685 = t361 * t449 + t362 * t443;
t683 = -m(7) - m(6) - m(5);
t680 = t724 + (t608 - t611 + t725) * t413;
t678 = -m(3) * pkin(7) + m(4) * t452 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t677 = -t506 + t680;
t676 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t674 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t604 = mrSges(6,3) * t252;
t192 = -mrSges(6,2) * t273 + t604;
t603 = mrSges(6,3) * t253;
t193 = mrSges(6,1) * t273 - t603;
t61 = mrSges(6,1) * t147 - mrSges(6,3) * t104;
t673 = m(6) * ((-t440 * t93 - t446 * t92) * qJD(5) + t726) - t193 * t538 - t192 * t539 - t440 * t61;
t395 = -mrSges(3,1) * t449 + mrSges(3,2) * t443;
t672 = m(3) * pkin(1) + m(4) * t418 + mrSges(2,1) - t395 + t506 - t724;
t671 = (t450 * t743 + t728) * g(1) + (t444 * t743 + t727) * g(2);
t647 = t253 / 0.2e1;
t637 = t345 / 0.2e1;
t632 = pkin(2) * t443;
t411 = pkin(3) * t430;
t627 = pkin(3) * t447;
t625 = pkin(5) * t253;
t616 = g(3) * t412;
t607 = mrSges(4,3) * t345;
t606 = mrSges(5,3) * t162;
t602 = Ifges(3,4) * t443;
t601 = Ifges(3,4) * t449;
t600 = Ifges(4,4) * t345;
t594 = pkin(5) * qJD(6);
t590 = t285 * mrSges(4,3);
t296 = qJD(2) * t366 + t466;
t297 = -qJD(2) * t368 - t467;
t182 = qJD(4) * t481 + t296 * t447 + t297 * t441;
t584 = t182 * t446;
t575 = t291 * t440;
t574 = t291 * t446;
t566 = t427 * t444;
t565 = t427 * t450;
t564 = t429 * t444;
t563 = t429 * t450;
t561 = t440 * t444;
t560 = t440 * t450;
t557 = t444 * t446;
t555 = t446 * t450;
t317 = t413 * t566 + t563;
t318 = -t413 * t564 + t565;
t554 = -t317 * mrSges(7,1) + t318 * mrSges(7,2);
t319 = -t413 * t565 + t564;
t320 = t413 * t563 + t566;
t553 = t319 * mrSges(7,1) - t320 * mrSges(7,2);
t547 = t411 + t433;
t544 = qJD(2) * t443;
t532 = Ifges(7,5) * t42 + Ifges(7,6) * t43 + Ifges(7,3) * t140;
t423 = pkin(2) * t544;
t524 = Ifges(6,5) * t104 + Ifges(6,6) * t105 + Ifges(6,3) * t147;
t523 = t411 + t688;
t522 = qJD(2) * t452;
t518 = t291 * t538;
t278 = -pkin(3) * t297 + t423;
t380 = t443 * t522;
t381 = t449 * t522;
t230 = t448 * t380 + t442 * t381 + t398 * t542 + t400 * t543;
t198 = pkin(9) * t297 + t230;
t231 = -qJD(3) * t307 - t380 * t442 + t448 * t381;
t199 = -pkin(9) * t296 + t231;
t78 = qJD(4) * t696 + t198 * t447 + t199 * t441;
t183 = qJD(4) * t291 + t296 * t441 - t447 * t297;
t89 = pkin(4) * t183 - pkin(10) * t182 + t278;
t510 = -t440 * t78 + t446 * t89;
t507 = t537 / 0.2e1;
t113 = t446 * t205 - t210 * t440;
t339 = -pkin(2) * t559 + t417 * t447;
t391 = pkin(10) * t569;
t499 = -t444 * t626 + t391;
t392 = pkin(10) * t568;
t498 = -t450 * t626 + t392;
t334 = -pkin(4) - t339;
t495 = t411 + t691;
t494 = mrSges(3,1) * t443 + mrSges(3,2) * t449;
t489 = -mrSges(7,1) * t427 - mrSges(7,2) * t429;
t487 = t602 + t707;
t485 = Ifges(3,5) * t449 - Ifges(3,6) * t443;
t86 = -pkin(5) * t481 - pkin(11) * t574 + t113;
t97 = -pkin(11) * t575 + t114;
t46 = -t439 * t97 + t445 * t86;
t47 = t439 * t86 + t445 * t97;
t483 = -t440 * t92 + t446 * t93;
t475 = t532 + t676;
t474 = pkin(1) * t494;
t331 = -t413 * t560 + t557;
t329 = t413 * t561 + t555;
t473 = t182 * t440 + t518;
t472 = t291 * t539 - t584;
t470 = t443 * (Ifges(3,1) * t449 - t602);
t19 = t205 * t538 - t210 * t539 + t440 * t89 + t446 * t78;
t79 = qJD(4) * t210 + t198 * t441 - t447 * t199;
t267 = Ifges(4,2) * t344 + Ifges(4,6) * t435 + t600;
t337 = Ifges(4,4) * t344;
t268 = t345 * Ifges(4,1) + t435 * Ifges(4,5) + t337;
t453 = -(-Ifges(4,2) * t345 + t268 + t337) * t344 / 0.2e1 + (Ifges(5,4) * t642 + t719 + t740) * t692 - t345 * (Ifges(4,1) * t344 - t600) / 0.2e1 + t396 * (mrSges(4,1) * t345 + mrSges(4,2) * t344) + (t606 + t742) * t482 + t739 + t267 * t637 - t435 * (Ifges(4,5) * t344 - Ifges(4,6) * t345) / 0.2e1 + Ifges(4,3) * t434 + t344 * t590 - t187 * mrSges(4,2) + t188 * mrSges(4,1) + t286 * t607 + Ifges(4,5) * t256 + Ifges(4,6) * t257;
t436 = -pkin(9) + t452;
t420 = Ifges(3,4) * t545;
t416 = -pkin(4) - t627;
t390 = -t415 - t627;
t384 = -t629 - t632;
t374 = pkin(1) + t547;
t357 = t450 * t384;
t356 = t444 * t384;
t343 = Ifges(3,1) * t546 + Ifges(3,5) * qJD(2) + t420;
t342 = Ifges(3,6) * qJD(2) + qJD(1) * t487;
t332 = t413 * t555 + t561;
t330 = -t413 * t557 + t560;
t321 = t334 - t623;
t316 = mrSges(4,1) * t435 - t607;
t315 = -mrSges(4,2) * t435 + mrSges(4,3) * t344;
t314 = t421 + t630;
t284 = -mrSges(4,1) * t344 + mrSges(4,2) * t345;
t260 = -mrSges(5,2) * t426 + mrSges(5,3) * t692;
t237 = -mrSges(4,2) * t434 + mrSges(4,3) * t257;
t236 = mrSges(4,1) * t434 - mrSges(4,3) * t256;
t219 = t367 * t291;
t217 = -mrSges(5,1) * t692 + mrSges(5,2) * t482;
t155 = pkin(5) * t575 - t696;
t133 = mrSges(7,1) * t265 - mrSges(7,3) * t172;
t132 = -mrSges(7,2) * t265 + mrSges(7,3) * t504;
t129 = -mrSges(5,2) * t425 + mrSges(5,3) * t149;
t128 = mrSges(5,1) * t425 - mrSges(5,3) * t148;
t96 = -mrSges(7,1) * t504 + mrSges(7,2) * t172;
t62 = -mrSges(6,2) * t147 + mrSges(6,3) * t105;
t58 = -t367 * t182 + t220 * t682;
t57 = -t182 * t478 - t291 * t295;
t56 = pkin(5) * t473 + t79;
t55 = -mrSges(6,1) * t105 + mrSges(6,2) * t104;
t29 = t445 * t74 - t588;
t28 = -t439 * t74 - t586;
t22 = -mrSges(7,2) * t140 + mrSges(7,3) * t43;
t21 = mrSges(7,1) * t140 - mrSges(7,3) * t42;
t20 = -qJD(5) * t114 + t510;
t18 = -pkin(11) * t473 + t19;
t17 = -pkin(11) * t584 + pkin(5) * t183 + (-t200 + (pkin(11) * t291 - t205) * t440) * qJD(5) + t510;
t13 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t6 = -qJD(6) * t47 + t17 * t445 - t18 * t439;
t5 = qJD(6) * t46 + t17 * t439 + t18 * t445;
t1 = [(-t561 * t668 - t332 * mrSges(6,1) - t320 * mrSges(7,1) - t331 * mrSges(6,2) - t319 * mrSges(7,2) + t683 * (t450 * t374 - t436 * t444) + t678 * t444 + (-t672 + t720) * t450) * g(2) + (Ifges(7,1) * t57 + Ifges(7,4) * t58) * t651 + (Ifges(7,5) * t57 + Ifges(7,6) * t58) * t645 + (t223 * mrSges(5,2) - t54 * mrSges(5,3) + Ifges(5,1) * t148 + Ifges(5,4) * t149 + Ifges(5,5) * t425 + t143 * t511 + t484 * t655 + t486 * t657 + t488 * t658 + t490 * t51) * t291 + t383 * t601 / 0.2e1 + t741 * t79 + (t449 * t601 + t470) * t507 + m(4) * (t187 * t307 + t188 * t305 + t230 * t286 + t231 * t285 - t336 * t418 - t396 * t423) - (mrSges(5,1) * t223 - mrSges(5,3) * t53 - Ifges(5,4) * t148 + Ifges(6,5) * t658 + Ifges(7,5) * t664 - Ifges(5,2) * t149 - Ifges(5,6) * t425 + Ifges(6,6) * t657 + Ifges(7,6) * t663 + Ifges(6,3) * t655 + Ifges(7,3) * t656 + t674 + t676) * t481 - (t532 + t524) * t481 / 0.2e1 + (-t330 * mrSges(6,1) - t318 * mrSges(7,1) - t329 * mrSges(6,2) - t317 * mrSges(7,2) + (-t436 * t683 + t678 - t715) * t450 + (-m(6) * (-t374 - t688) - m(7) * (-t374 - t691) + m(5) * t374 + t672) * t444) * g(1) + (-t15 * t575 - t16 * t574 + t472 * t92 - t473 * t93) * mrSges(6,3) + t252 * (-Ifges(6,4) * t472 - Ifges(6,2) * t473) / 0.2e1 + (Ifges(7,4) * t57 + Ifges(7,2) * t58) * t653 + (-t606 + t710 / 0.2e1 + t709 / 0.2e1 + Ifges(7,6) * t653 + Ifges(7,3) * t645 + Ifges(6,5) * t647 + Ifges(7,5) * t651 - Ifges(5,4) * t640 - t711 / 0.2e1 - t712 / 0.2e1 - t207 / 0.2e1 + t670 + t703 / 0.2e1) * t183 + m(5) * (t162 * t78 + t210 * t53 + t223 * t324 + t278 * t301) + t273 * (-Ifges(6,5) * t472 - Ifges(6,6) * t473) / 0.2e1 + (Ifges(3,4) * t383 + Ifges(3,2) * t382) * t635 + (-t562 / 0.2e1 + Ifges(5,1) * t640 + t708 / 0.2e1 + t269 / 0.2e1 + t515 - t740) * t182 + (mrSges(4,2) * t336 - mrSges(4,3) * t188 + Ifges(4,1) * t256 + Ifges(4,4) * t257 + Ifges(4,5) * t434) * t368 + (-mrSges(4,1) * t336 + mrSges(4,3) * t187 + Ifges(4,4) * t256 + Ifges(4,2) * t257 + Ifges(4,6) * t434) * t366 - t395 * t585 + t487 * t717 - t32 * t575 / 0.2e1 + (-Ifges(6,1) * t472 - Ifges(6,4) * t473) * t647 + t121 * (-mrSges(7,1) * t58 + mrSges(7,2) * t57) + t156 * (mrSges(6,1) * t473 - mrSges(6,2) * t472) - t296 * t590 + (-Ifges(7,1) * t220 - Ifges(7,4) * t219) * t664 + (-Ifges(7,5) * t220 - Ifges(7,6) * t219) * t656 + (-t219 * t3 + t220 * t4 - t26 * t57 + t27 * t58) * mrSges(7,3) + (-Ifges(7,4) * t220 - Ifges(7,2) * t219) * t663 + t25 * (mrSges(7,1) * t219 - mrSges(7,2) * t220) + (t382 * t621 + t383 * t622 + t685) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t685) + (t343 * t635 + t485 * qJD(2) / 0.2e1 - t686) * qJD(2) + t113 * t61 + t114 * t62 + t56 * t96 - t342 * t544 / 0.2e1 - (-m(5) * t54 + m(6) * t51 - t128 + t55) * t696 - t474 * t537 + Ifges(2,3) * qJDD(1) + t46 * t21 + t47 * t22 + m(7) * (t121 * t56 + t155 * t25 + t26 * t6 + t27 * t5 + t3 * t47 + t4 * t46) - t142 * t518 / 0.2e1 + t57 * t659 + t58 * t661 + t574 * t665 - t220 * t666 - t219 * t667 + (Ifges(4,1) * t296 + Ifges(4,4) * t297) * t637 + t435 * (Ifges(4,5) * t296 + Ifges(4,6) * t297) / 0.2e1 - t418 * (-mrSges(4,1) * t257 + mrSges(4,2) * t256) - t396 * (-mrSges(4,1) * t297 + mrSges(4,2) * t296) - pkin(1) * (-mrSges(3,1) * t382 + mrSges(3,2) * t383) + t344 * (Ifges(4,4) * t296 + Ifges(4,2) * t297) / 0.2e1 + t324 * (-mrSges(5,1) * t149 + mrSges(5,2) * t148) + t230 * t315 + t231 * t316 + t305 * t236 + t307 * t237 + t296 * t268 / 0.2e1 + t297 * t267 / 0.2e1 + t278 * t217 + m(6) * (t113 * t16 + t114 * t15 + t19 * t93 + t20 * t92) + t5 * t132 + t6 * t133 + t155 * t13 + t284 * t423 + (-mrSges(3,1) * t622 - mrSges(3,2) * t621 + 0.2e1 * Ifges(3,6) * t635) * qJDD(2) + (Ifges(3,1) * t383 + Ifges(3,4) * t717 + Ifges(3,5) * qJDD(2) - t507 * t707) * t443 + t19 * t192 + t20 * t193 + t210 * t129 + t286 * t297 * mrSges(4,3) + t78 * t260; (t162 * t694 - t301 * t314 + t339 * t54 + t340 * t53) * m(5) + (-t100 * t92 - t101 * t93 + t279 * t483 + t334 * t51 - g(1) * (t357 + t498) - g(2) * (t356 + t499)) * m(6) + t741 * t732 + t673 * t335 - m(4) * (t285 * t292 + t286 * t293 - t396 * t421) + t671 + t684 * (m(4) * t632 - m(5) * t384 + t494 + t730) + t694 * t260 + t695 * t96 + (-m(4) * t433 - m(5) * t547 - m(6) * (t433 + t523) + t395 - m(7) * (t433 + t495) + t677) * g(3) + t62 * t573 + t453 - (-Ifges(3,2) * t546 + t343 + t420) * t545 / 0.2e1 + t751 * t192 + (t686 + (-t470 / 0.2e1 + t474) * qJD(1)) * qJD(1) + t750 * t193 + (-t316 * t543 + m(4) * (t187 * t442 + t188 * t448 + (-t285 * t442 + t286 * t448) * qJD(3)) + t315 * t542 + t442 * t237) * pkin(2) + t342 * t546 / 0.2e1 - t485 * t537 / 0.2e1 - t284 * t421 + Ifges(3,3) * qJDD(2) + t236 * t631 + Ifges(3,6) * t382 + Ifges(3,5) * t383 - t361 * mrSges(3,2) - t362 * mrSges(3,1) + t339 * t128 + t340 * t129 + t334 * t55 - t314 * t217 - t293 * t315 - t292 * t316 + t321 * t13 + t232 * t21 + t233 * t22 + t704 * t133 + t705 * t132 + (-g(1) * t357 - g(2) * t356 + t121 * t695 + t232 * t4 + t233 * t3 + t25 * t321 + t26 * t704 + t27 * t705) * m(7); (t416 * t51 + (t156 * t441 + t447 * t483) * qJD(4) * pkin(3) - t156 * t166 - t92 * t98 - t93 * t99 - g(1) * t392 - g(2) * t391) * m(6) + t673 * t414 + ((t441 * t53 + t447 * t54 + (-t161 * t441 + t162 * t447) * qJD(4)) * pkin(3) + t161 * t166 - t162 * t167 - t301 * t630) * m(5) + (t444 * t687 + t727) * g(2) + (t450 * t687 + t728) * g(1) + t729 * t697 + (m(5) * t629 + t730) * t684 + (t528 - t167) * t260 + t701 * t132 + (-m(5) * t411 - m(6) * t523 - m(7) * t495 + t677) * g(3) + t62 * t567 + t453 + t755 * t192 + t690 * t96 - t217 * t630 + t128 * t627 + t129 * t628 + t416 * t55 + t390 * t13 - t285 * t315 + t286 * t316 + t281 * t21 + t282 * t22 + t754 * t193 + (t121 * t690 + t25 * t390 + t26 * t702 + t27 * t701 + t281 * t4 + t282 * t3) * m(7) + t702 * t133; (t680 + t720) * g(3) + t673 * pkin(10) + (-pkin(4) * t51 - g(1) * t498 - g(2) * t499 - t108 * t92 - t109 * t93) * m(6) + t671 + t698 * t96 + t699 * t132 + (t121 * t698 - t25 * t415 + t26 * t700 + t27 * t699 + t3 * t306 + t304 * t4) * m(7) + t700 * t133 + (t605 + t675) * t162 + (t269 + t208) * t642 + t684 * t491 + t739 + t742 * t482 + (-t679 + t719) * t692 - pkin(4) * t55 - t415 * t13 + t304 * t21 + t306 * t22 + t62 * t619 - t109 * t192 - t108 * t193 - t161 * t260; t674 + (mrSges(6,2) * t332 - t331 * t706 - t553) * g(1) + (-mrSges(6,2) * t330 + t329 * t706 - t554) * g(2) + (t604 - t192) * t92 + t524 + (-Ifges(6,2) * t253 + t143 + t246) * t649 + (t603 + t193) * t93 + (-t439 * t594 - t28) * t133 + (t21 * t445 + t22 * t439) * pkin(5) - (Ifges(7,4) * t652 + Ifges(7,2) * t654 + Ifges(7,6) * t646 + t662 - t723) * t172 + (Ifges(7,1) * t652 + Ifges(7,4) * t654 + Ifges(7,5) * t646 + t660 - t722) * t504 + (t445 * t594 - t29) * t132 - t96 * t625 - m(7) * (t121 * t625 + t26 * t28 + t27 * t29) + t475 + (-t489 + t490 + t715) * t616 + (t3 * t439 + t4 * t445 + (-t26 * t439 + t27 * t445) * qJD(6)) * t668 + (Ifges(6,5) * t252 - Ifges(6,6) * t253) * t644 + t142 * t647 + (Ifges(6,1) * t252 - t598) * t648 - t156 * (mrSges(6,1) * t253 + mrSges(6,2) * t252); -t121 * (mrSges(7,1) * t172 + mrSges(7,2) * t504) + (Ifges(7,1) * t504 - t595) * t652 + t82 * t651 + (Ifges(7,5) * t504 - Ifges(7,6) * t172) * t646 - t26 * t132 + t27 * t133 - g(1) * t553 - g(2) * t554 - t489 * t616 + (t172 * t27 + t26 * t504) * mrSges(7,3) + t475 + (-Ifges(7,2) * t172 + t165 + t83) * t654;];
tau  = t1;
