% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:26
% EndTime: 2019-03-09 11:05:38
% DurationCPUTime: 45.35s
% Computational Cost: add. (17652->937), mult. (43020->1222), div. (0->0), fcn. (35056->14), ass. (0->440)
t349 = sin(qJ(4));
t345 = cos(pkin(11));
t526 = cos(qJ(4));
t437 = t526 * t345;
t344 = sin(pkin(6));
t353 = cos(qJ(2));
t477 = t344 * t353;
t393 = t437 * t477;
t343 = sin(pkin(11));
t465 = qJD(1) * t344;
t435 = t353 * t465;
t418 = t343 * t435;
t227 = qJD(1) * t393 - t349 * t418;
t430 = qJD(4) * t526;
t458 = qJD(4) * t349;
t281 = t343 * t458 - t345 * t430;
t651 = t227 + t281;
t293 = t343 * t526 + t349 * t345;
t363 = t293 * t477;
t226 = qJD(1) * t363;
t282 = t293 * qJD(4);
t641 = t226 - t282;
t346 = cos(pkin(6));
t464 = qJD(1) * t346;
t330 = qJD(2) + t464;
t350 = sin(qJ(2));
t436 = t350 * t465;
t369 = -t330 * t345 + t343 * t436;
t234 = t526 * t369;
t247 = t330 * t343 + t345 * t436;
t177 = t247 * t349 + t234;
t310 = -qJD(4) + t435;
t358 = t247 * t526 - t349 * t369;
t501 = t358 * Ifges(6,6);
t101 = -t310 * Ifges(6,5) + t177 * Ifges(6,3) - t501;
t502 = t358 * Ifges(5,4);
t104 = -t177 * Ifges(5,2) - t310 * Ifges(5,6) + t502;
t320 = pkin(8) * t435;
t453 = pkin(1) * t464;
t275 = t350 * t453 + t320;
t230 = qJ(3) * t330 + t275;
t261 = (-pkin(2) * t353 - qJ(3) * t350 - pkin(1)) * t344;
t236 = qJD(1) * t261;
t152 = t345 * t230 + t343 * t236;
t121 = -pkin(9) * t369 + t152;
t151 = -t230 * t343 + t345 * t236;
t361 = -pkin(3) * t435 - pkin(9) * t247 + t151;
t359 = t349 * t361;
t60 = t121 * t526 + t359;
t58 = t310 * qJ(5) - t60;
t650 = -t58 * mrSges(6,1) + t60 * mrSges(5,3) - t101 / 0.2e1 + t104 / 0.2e1;
t461 = qJD(2) * t353;
t279 = (qJD(1) * t461 + qJDD(1) * t350) * t344;
t454 = qJDD(1) * t346;
t329 = qJDD(2) + t454;
t209 = t279 * t345 + t329 * t343;
t392 = t279 * t343 - t329 * t345;
t92 = qJD(4) * t234 - t526 * t209 + t247 * t458 + t349 * t392;
t560 = -t92 / 0.2e1;
t93 = qJD(4) * t358 + t349 * t209 + t392 * t526;
t558 = -t93 / 0.2e1;
t462 = qJD(2) * t350;
t434 = t344 * t462;
t628 = -qJD(1) * t434 + qJDD(1) * t477;
t266 = qJDD(4) - t628;
t537 = t266 / 0.2e1;
t219 = pkin(3) * t418 + t275;
t649 = qJ(5) * t651 - qJD(5) * t293 - t219;
t608 = mrSges(5,2) - mrSges(6,3);
t615 = m(7) + m(6);
t391 = -qJ(5) * t615 + t608;
t348 = sin(qJ(6));
t352 = cos(qJ(6));
t406 = mrSges(7,1) * t348 + mrSges(7,2) * t352;
t581 = t391 - t406;
t399 = pkin(2) * t350 - qJ(3) * t353;
t273 = t399 * t465;
t274 = -pkin(8) * t436 + t353 * t453;
t187 = t345 * t273 - t274 * t343;
t475 = t345 * t353;
t371 = (pkin(3) * t350 - pkin(9) * t475) * t344;
t147 = qJD(1) * t371 + t187;
t188 = t343 * t273 + t345 * t274;
t165 = -pkin(9) * t418 + t188;
t518 = pkin(9) + qJ(3);
t303 = t518 * t343;
t304 = t518 * t345;
t225 = -t349 * t303 + t304 * t526;
t600 = -qJD(3) * t293 - qJD(4) * t225 - t147 * t526 + t349 * t165;
t171 = Ifges(6,6) * t177;
t103 = -t310 * Ifges(6,4) - Ifges(6,2) * t358 + t171;
t115 = t526 * t361;
t59 = t121 * t349 - t115;
t57 = pkin(4) * t310 + qJD(5) + t59;
t648 = t103 / 0.2e1 - mrSges(6,1) * t57 - mrSges(5,3) * t59;
t557 = t93 / 0.2e1;
t647 = t557 - t558;
t559 = t92 / 0.2e1;
t646 = t559 - t560;
t533 = -t310 / 0.2e1;
t541 = t358 / 0.2e1;
t542 = -t358 / 0.2e1;
t543 = t177 / 0.2e1;
t544 = -t177 / 0.2e1;
t221 = -t330 * pkin(2) + qJD(3) - t274;
t174 = pkin(3) * t369 + t221;
t356 = -qJ(5) * t358 + t174;
t67 = t177 * pkin(4) + t356;
t583 = -t174 * mrSges(5,1) + t67 * mrSges(6,2);
t606 = Ifges(6,5) - Ifges(5,6);
t645 = -Ifges(5,4) * t541 - Ifges(5,2) * t544 + Ifges(6,6) * t542 + Ifges(6,3) * t543 + t533 * t606 - t583 - t650;
t644 = m(5) + t615;
t519 = -mrSges(6,2) + mrSges(5,1);
t138 = t177 * t352 + t310 * t348;
t139 = t177 * t348 - t310 * t352;
t172 = Ifges(5,4) * t177;
t173 = qJD(6) + t358;
t629 = Ifges(5,1) * t358 - t310 * Ifges(5,5) + t139 * Ifges(7,5) + t138 * Ifges(7,6) + t173 * Ifges(7,3) - t172;
t479 = t344 * t350;
t555 = pkin(4) + pkin(10);
t425 = t555 * t479;
t643 = -pkin(5) * t651 + qJD(1) * t425 - t600;
t460 = qJD(3) * t343;
t599 = qJD(3) * t437 - t526 * t165 - t303 * t430 + (-qJD(4) * t304 - t147 - t460) * t349;
t642 = t555 * t641 - t649;
t48 = qJD(6) * t138 + t266 * t352 + t348 * t93;
t91 = qJDD(6) - t92;
t25 = mrSges(7,1) * t91 - mrSges(7,3) * t48;
t49 = -qJD(6) * t139 - t266 * t348 + t352 * t93;
t26 = -mrSges(7,2) * t91 + mrSges(7,3) * t49;
t96 = -mrSges(7,2) * t173 + mrSges(7,3) * t138;
t97 = mrSges(7,1) * t173 - mrSges(7,3) * t139;
t396 = -t348 * t97 + t352 * t96;
t640 = t396 * qJD(6) + t352 * t25 + t348 * t26;
t524 = pkin(1) * t346;
t452 = qJD(2) * t524;
t422 = qJD(1) * t452;
t447 = pkin(1) * t454;
t195 = pkin(8) * t628 + t350 * t447 + t353 * t422;
t160 = qJ(3) * t329 + qJD(3) * t330 + t195;
t459 = qJD(3) * t350;
t168 = -pkin(2) * t628 - qJ(3) * t279 + (-pkin(1) * qJDD(1) - qJD(1) * t459) * t344;
t99 = -t160 * t343 + t345 * t168;
t62 = -pkin(3) * t628 - pkin(9) * t209 + t99;
t100 = t345 * t160 + t343 * t168;
t70 = -pkin(9) * t392 + t100;
t13 = qJD(4) * t115 - t121 * t458 + t349 * t62 + t526 * t70;
t10 = -qJ(5) * t266 + qJD(5) * t310 - t13;
t331 = pkin(8) * t479;
t196 = -qJD(2) * t320 - qJDD(1) * t331 - t350 * t422 + t353 * t447;
t181 = -t329 * pkin(2) + qJDD(3) - t196;
t120 = pkin(3) * t392 + t181;
t355 = t92 * qJ(5) - qJD(5) * t358 + t120;
t20 = t93 * pkin(4) + t355;
t614 = -t266 / 0.2e1;
t637 = mrSges(5,1) * t120 + mrSges(6,1) * t10 - mrSges(6,2) * t20 - mrSges(5,3) * t13 + Ifges(5,4) * t646 + Ifges(5,2) * t647 + Ifges(5,6) * t614 + 0.2e1 * Ifges(6,6) * t559 + 0.2e1 * Ifges(6,3) * t557 + (Ifges(6,5) + t606) * t537;
t14 = -qJD(4) * t359 - t121 * t430 - t349 * t70 + t526 * t62;
t367 = qJDD(5) - t14;
t11 = -t266 * pkin(4) + t367;
t561 = t91 / 0.2e1;
t566 = t49 / 0.2e1;
t567 = t48 / 0.2e1;
t12 = t555 * t93 + t355;
t386 = pkin(5) * t358 + t59;
t624 = t386 + qJD(5);
t40 = t310 * t555 + t624;
t51 = t177 * t555 + t356;
t15 = -t348 * t51 + t352 * t40;
t5 = -t92 * pkin(5) - t266 * t555 + t367;
t1 = qJD(6) * t15 + t12 * t352 + t348 * t5;
t16 = t348 * t40 + t352 * t51;
t2 = -qJD(6) * t16 - t12 * t348 + t352 * t5;
t584 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t607 = Ifges(6,4) - Ifges(5,5);
t7 = Ifges(7,5) * t48 + Ifges(7,6) * t49 + Ifges(7,3) * t91;
t635 = t584 + mrSges(6,1) * t11 + mrSges(5,2) * t120 - mrSges(5,3) * t14 - mrSges(6,3) * t20 + 0.2e1 * Ifges(5,1) * t560 + 0.2e1 * Ifges(5,4) * t558 + Ifges(6,4) * t614 + Ifges(7,5) * t567 + Ifges(7,6) * t566 + Ifges(7,3) * t561 - t647 * Ifges(6,6) - t646 * Ifges(6,2) + t7 / 0.2e1 + (-t607 + Ifges(5,5)) * t537;
t545 = t173 / 0.2e1;
t547 = t139 / 0.2e1;
t549 = t138 / 0.2e1;
t577 = t15 * mrSges(7,1) + t174 * mrSges(5,2) - t16 * mrSges(7,2) - t67 * mrSges(6,3);
t634 = -Ifges(5,1) * t541 - Ifges(5,4) * t544 - Ifges(7,5) * t547 + Ifges(6,2) * t542 + Ifges(6,6) * t543 - Ifges(7,6) * t549 - Ifges(7,3) * t545 + t533 * t607 - t577 + t648;
t630 = t196 * mrSges(3,1);
t602 = qJ(5) * t436 - t599;
t616 = m(7) * pkin(10);
t627 = pkin(4) * t615 + t519 + t616;
t623 = Ifges(3,5) * t279 + Ifges(3,3) * t329;
t532 = t310 / 0.2e1;
t546 = -t173 / 0.2e1;
t548 = -t139 / 0.2e1;
t550 = -t138 / 0.2e1;
t622 = -Ifges(5,1) * t542 - Ifges(7,5) * t548 + Ifges(6,2) * t541 - Ifges(7,6) * t550 - Ifges(7,3) * t546 + t532 * t607 + t577;
t613 = -t392 / 0.2e1;
t292 = t343 * t349 - t437;
t338 = pkin(3) * t345 + pkin(2);
t388 = -qJ(5) * t293 - t338;
t180 = t292 * t555 + t388;
t224 = t526 * t303 + t304 * t349;
t189 = pkin(5) * t293 + t224;
t111 = -t180 * t348 + t189 * t352;
t605 = qJD(6) * t111 + t348 * t643 - t352 * t642;
t112 = t180 * t352 + t189 * t348;
t604 = -qJD(6) * t112 + t348 * t642 + t352 * t643;
t517 = mrSges(6,1) * t177;
t143 = mrSges(6,3) * t310 + t517;
t79 = -mrSges(7,1) * t138 + mrSges(7,2) * t139;
t494 = -t143 + t79;
t603 = pkin(5) * t641 - t602;
t601 = pkin(4) * t436 - t600;
t598 = -pkin(4) * t641 + t649;
t516 = mrSges(6,1) * t358;
t144 = -mrSges(6,2) * t310 + t516;
t514 = mrSges(5,3) * t358;
t146 = -mrSges(5,1) * t310 - t514;
t597 = t146 - t144;
t424 = mrSges(3,3) * t436;
t596 = -mrSges(3,1) * t330 + mrSges(4,1) * t369 + t247 * mrSges(4,2) + t424;
t191 = t226 * t352 - t348 * t436;
t457 = qJD(6) * t348;
t471 = t352 * t282;
t380 = t292 * t457 - t471;
t595 = t191 + t380;
t192 = t226 * t348 + t352 * t436;
t456 = qJD(6) * t352;
t488 = t282 * t348;
t381 = t292 * t456 + t488;
t594 = t192 - t381;
t280 = t343 * t346 + t345 * t479;
t445 = t343 * t479;
t382 = -t345 * t346 + t445;
t364 = t526 * t382;
t193 = t280 * t349 + t364;
t474 = t348 * t353;
t313 = t344 * t474;
t169 = t193 * t352 + t313;
t593 = t100 * t345 - t343 * t99;
t512 = Ifges(3,4) * t350;
t590 = pkin(1) * (mrSges(3,1) * t350 + mrSges(3,2) * t353) - t350 * (Ifges(3,1) * t353 - t512) / 0.2e1;
t588 = -mrSges(7,3) - t616;
t342 = pkin(11) + qJ(4);
t339 = sin(t342);
t340 = cos(t342);
t410 = -mrSges(4,1) * t345 + mrSges(4,2) * t343;
t373 = m(4) * pkin(2) - t410;
t587 = t339 * t608 - t340 * t519 - t373;
t407 = mrSges(7,1) * t352 - mrSges(7,2) * t348;
t521 = t177 * pkin(5);
t41 = -t58 - t521;
t503 = t139 * Ifges(7,4);
t54 = t138 * Ifges(7,2) + t173 * Ifges(7,6) + t503;
t565 = -t54 / 0.2e1;
t586 = t352 * t565 + t41 * t407;
t233 = (qJD(2) * t399 - t459) * t344;
t276 = -pkin(8) * t434 + t353 * t452;
t240 = qJD(3) * t346 + t276;
t161 = t345 * t233 - t240 * t343;
t125 = qJD(2) * t371 + t161;
t289 = pkin(8) * t477 + t350 * t524;
t260 = qJ(3) * t346 + t289;
t183 = -t260 * t343 + t345 * t261;
t127 = -pkin(3) * t477 - pkin(9) * t280 + t183;
t162 = t343 * t233 + t345 * t240;
t433 = t344 * t461;
t417 = t343 * t433;
t140 = -pkin(9) * t417 + t162;
t184 = t345 * t260 + t343 * t261;
t142 = -pkin(9) * t382 + t184;
t374 = -t349 * t125 - t127 * t430 - t526 * t140 + t142 * t458;
t27 = -t344 * (qJ(5) * t462 - qJD(5) * t353) + t374;
t582 = -m(4) * qJ(3) - m(7) * pkin(5) - mrSges(6,1) - mrSges(4,3) - mrSges(5,3);
t527 = cos(qJ(1));
t439 = t527 * t350;
t351 = sin(qJ(1));
t472 = t351 * t353;
t284 = t346 * t439 + t472;
t440 = t344 * t527;
t210 = t284 * t339 + t340 * t440;
t438 = t527 * t353;
t473 = t350 * t351;
t286 = -t346 * t473 + t438;
t478 = t344 * t351;
t214 = t286 * t339 - t340 * t478;
t258 = t339 * t479 - t346 * t340;
t580 = g(1) * t214 + g(2) * t210 + g(3) * t258;
t579 = t339 * t406 + mrSges(3,1) - t587;
t578 = t59 * mrSges(5,1) + t60 * mrSges(5,2) - t57 * mrSges(6,2) + t58 * mrSges(6,3);
t576 = -t407 + t582;
t575 = mrSges(3,2) + t576;
t574 = Ifges(6,3) * t544 + t532 * t606 + t583;
t572 = t344 ^ 2;
t571 = Ifges(7,1) * t567 + Ifges(7,4) * t566 + Ifges(7,5) * t561;
t564 = t54 / 0.2e1;
t132 = Ifges(7,4) * t138;
t55 = t139 * Ifges(7,1) + t173 * Ifges(7,5) + t132;
t563 = -t55 / 0.2e1;
t562 = t55 / 0.2e1;
t539 = t209 / 0.2e1;
t536 = t280 / 0.2e1;
t531 = -t343 / 0.2e1;
t530 = t345 / 0.2e1;
t528 = t353 / 0.2e1;
t525 = pkin(1) * t344;
t523 = pkin(1) * t353;
t522 = t1 * t348;
t515 = mrSges(5,3) * t177;
t513 = mrSges(7,3) * t352;
t511 = Ifges(4,4) * t343;
t510 = Ifges(4,4) * t345;
t509 = Ifges(7,4) * t348;
t508 = Ifges(7,4) * t352;
t507 = Ifges(4,5) * t247;
t506 = Ifges(4,6) * t345;
t505 = Ifges(4,6) * t350;
t504 = Ifges(4,3) * t350;
t500 = t343 * Ifges(4,6);
t493 = qJ(5) * t177;
t492 = qJ(5) * t339;
t490 = t358 * t348;
t283 = -t346 * t438 + t473;
t487 = t283 * t340;
t285 = t346 * t472 + t439;
t486 = t285 * t340;
t485 = t292 * t348;
t484 = t292 * t352;
t480 = t343 * t353;
t476 = t345 * (Ifges(4,1) * t247 - Ifges(4,4) * t369 - Ifges(4,5) * t435);
t470 = t352 * t353;
t78 = t349 * t127 + t526 * t142;
t469 = -t283 * t338 + t284 * t518;
t468 = -t285 * t338 + t286 * t518;
t277 = pkin(8) * t433 + t350 * t452;
t466 = t527 * pkin(1) + pkin(8) * t478;
t333 = pkin(4) * t477;
t451 = -Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t450 = Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t449 = -Ifges(5,3) / 0.2e1 - Ifges(6,1) / 0.2e1;
t446 = qJ(5) * t477;
t444 = t344 * t470;
t442 = Ifges(3,6) * t628 + t623;
t220 = pkin(3) * t417 + t277;
t429 = -t465 / 0.2e1;
t427 = -t457 / 0.2e1;
t76 = -t92 * mrSges(6,1) + t266 * mrSges(6,2);
t426 = -pkin(1) * t351 + pkin(8) * t440;
t211 = t284 * t340 - t339 * t440;
t423 = mrSges(3,3) * t435;
t421 = -pkin(4) * t487 - t283 * t492 + t469;
t420 = -pkin(4) * t486 - t285 * t492 + t468;
t419 = t343 * t440;
t414 = t353 * t429;
t413 = t506 / 0.2e1 - Ifges(3,6) / 0.2e1;
t412 = -t519 + t588;
t77 = t127 * t526 - t349 * t142;
t409 = mrSges(4,1) * t343 + mrSges(4,2) * t345;
t404 = Ifges(4,1) * t345 - t511;
t403 = Ifges(7,1) * t348 + t508;
t402 = -Ifges(4,2) * t343 + t510;
t401 = Ifges(7,2) * t352 + t509;
t400 = Ifges(7,5) * t348 + Ifges(7,6) * t352;
t397 = t15 * t348 - t16 * t352;
t194 = t280 * t526 - t349 * t382;
t72 = t333 - t77;
t52 = t194 * pkin(5) + pkin(10) * t477 + t72;
t197 = pkin(3) * t445 + t331 + (-t338 - t523) * t346;
t357 = -t194 * qJ(5) + t197;
t65 = t193 * t555 + t357;
t23 = -t348 * t65 + t352 * t52;
t24 = t348 * t52 + t352 * t65;
t395 = -t348 * t96 - t352 * t97;
t394 = t343 * pkin(3) * t478 + t285 * t518 + t286 * t338 + t466;
t71 = t446 - t78;
t383 = -t193 * t348 + t444;
t379 = pkin(3) * t419 - t283 * t518 - t284 * t338 + t426;
t378 = t221 * t409;
t375 = -t125 * t526 + t127 * t458 + t349 * t140 + t142 * t430;
t366 = mrSges(3,2) + t582;
t131 = mrSges(4,1) * t392 + t209 * mrSges(4,2);
t133 = qJD(4) * t364 - qJD(2) * t393 + (qJD(4) * t280 + t417) * t349;
t365 = qJ(5) * t133 - qJD(5) * t194 + t220;
t362 = -t14 * mrSges(5,1) + t13 * mrSges(5,2) - t11 * mrSges(6,2) + t10 * mrSges(6,3);
t360 = -qJD(6) * t397 + t2 * t352 + t522;
t318 = Ifges(3,4) * t435;
t297 = t338 * t477;
t288 = t346 * t523 - t331;
t287 = (-mrSges(3,1) * t353 + mrSges(3,2) * t350) * t344;
t272 = -mrSges(3,2) * t330 + t423;
t263 = t331 + (-pkin(2) - t523) * t346;
t257 = Ifges(6,1) * t266;
t256 = Ifges(5,3) * t266;
t229 = Ifges(3,1) * t436 + t330 * Ifges(3,5) + t318;
t228 = t330 * Ifges(3,6) + (Ifges(3,2) * t353 + t512) * t465;
t215 = t286 * t340 + t339 * t478;
t208 = -mrSges(4,1) * t435 - mrSges(4,3) * t247;
t207 = mrSges(4,2) * t435 - mrSges(4,3) * t369;
t201 = pkin(4) * t292 + t388;
t190 = -t292 * pkin(5) + t225;
t167 = -mrSges(4,1) * t628 - mrSges(4,3) * t209;
t166 = mrSges(4,2) * t628 - mrSges(4,3) * t392;
t164 = t214 * t348 + t285 * t352;
t163 = t214 * t352 - t285 * t348;
t158 = Ifges(4,4) * t247 - Ifges(4,2) * t369 - Ifges(4,6) * t435;
t157 = -Ifges(4,6) * t369 - Ifges(4,3) * t435 + t507;
t145 = mrSges(5,2) * t310 - t515;
t134 = qJD(2) * t363 + qJD(4) * t194;
t117 = Ifges(4,1) * t209 - Ifges(4,4) * t392 - Ifges(4,5) * t628;
t116 = Ifges(4,4) * t209 - t392 * Ifges(4,2) - Ifges(4,6) * t628;
t110 = -mrSges(6,2) * t177 - mrSges(6,3) * t358;
t109 = mrSges(5,1) * t177 + mrSges(5,2) * t358;
t108 = pkin(4) * t358 + t493;
t105 = -t310 * Ifges(6,1) - Ifges(6,4) * t358 + t177 * Ifges(6,5);
t102 = Ifges(5,5) * t358 - t177 * Ifges(5,6) - t310 * Ifges(5,3);
t98 = t193 * pkin(4) + t357;
t90 = Ifges(6,4) * t92;
t89 = Ifges(5,5) * t92;
t88 = Ifges(6,5) * t93;
t87 = Ifges(5,6) * t93;
t86 = t92 * mrSges(5,2);
t85 = t92 * mrSges(6,3);
t81 = qJD(6) * t169 + t134 * t348 + t352 * t434;
t80 = qJD(6) * t383 + t134 * t352 - t348 * t434;
t75 = mrSges(6,1) * t93 - mrSges(6,3) * t266;
t74 = -mrSges(5,2) * t266 - mrSges(5,3) * t93;
t73 = mrSges(5,1) * t266 + mrSges(5,3) * t92;
t66 = t358 * t555 + t493;
t56 = -pkin(5) * t193 - t71;
t50 = pkin(4) * t134 + t365;
t45 = t60 - t521;
t39 = t93 * mrSges(5,1) - t86;
t38 = -t93 * mrSges(6,2) + t85;
t33 = t134 * t555 + t365;
t28 = -pkin(4) * t434 + t375;
t22 = -pkin(5) * t134 - t27;
t21 = -t133 * pkin(5) - qJD(2) * t425 + t375;
t19 = t348 * t45 + t352 * t66;
t18 = -t348 * t66 + t352 * t45;
t17 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t8 = Ifges(7,4) * t48 + Ifges(7,2) * t49 + Ifges(7,6) * t91;
t6 = -pkin(5) * t93 - t10;
t4 = -qJD(6) * t24 + t21 * t352 - t33 * t348;
t3 = qJD(6) * t23 + t21 * t348 + t33 * t352;
t9 = [(Ifges(4,4) * t280 - Ifges(4,2) * t382) * t613 + t637 * t193 + (-t629 / 0.2e1 + t634) * t133 + (-Ifges(7,5) * t383 + Ifges(7,6) * t169) * t561 + (-Ifges(7,1) * t383 + Ifges(7,4) * t169) * t567 + (t1 * t169 - t15 * t81 + t16 * t80 + t2 * t383) * mrSges(7,3) + t6 * (-mrSges(7,1) * t169 - mrSges(7,2) * t383) + (-Ifges(7,4) * t383 + Ifges(7,2) * t169) * t566 - t383 * t571 + (t351 * mrSges(2,1) + t527 * mrSges(2,2) - m(4) * (-pkin(2) * t284 + t426) - (-t284 * t345 + t419) * mrSges(4,1) - (t284 * t343 + t345 * t440) * mrSges(4,2) - m(3) * t426 + t284 * mrSges(3,1) - mrSges(3,3) * t440 - m(5) * t379 - t412 * t211 - t581 * t210 + (-t366 + t407) * t283 + t615 * (pkin(4) * t211 - t379)) * g(1) + t645 * t134 + t635 * t194 - (Ifges(4,5) * t536 - t289 * mrSges(3,3) + t413 * t346 + (-pkin(1) * mrSges(3,1) + (-Ifges(4,3) - Ifges(3,2)) * t353 + (-Ifges(3,4) - t500 / 0.2e1) * t350) * t344) * t628 + t596 * t277 + (-t287 * t525 + Ifges(2,3)) * qJDD(1) + (-t195 * t346 - t279 * t525 - t289 * t329) * mrSges(3,2) + t346 * t630 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t572 + t195 * t289 + t196 * t288 - t274 * t277 + t275 * t276) + (-t100 * t382 - t99 * t280) * mrSges(4,3) + t181 * (mrSges(4,1) * t382 + t280 * mrSges(4,2)) - t382 * t116 / 0.2e1 + m(4) * (t100 * t184 + t151 * t161 + t152 * t162 + t181 * t263 + t183 * t99 + t221 * t277) + m(7) * (t1 * t24 + t15 * t4 + t16 * t3 + t2 * t23 + t22 * t41 + t56 * t6) + m(6) * (t10 * t71 + t11 * t72 + t20 * t98 + t27 * t58 + t28 * t57 + t50 * t67) + ((t158 * t531 - t274 * mrSges(3,3) + t476 / 0.2e1 + t229 / 0.2e1 + t247 * t404 / 0.2e1 + t378 + (t402 * t530 + Ifges(3,5) / 0.2e1) * t330 + (-t151 * t345 - t152 * t343) * mrSges(4,3)) * t353 + (-t275 * mrSges(3,3) - t228 / 0.2e1 + t157 / 0.2e1 + t102 / 0.2e1 + t105 / 0.2e1 + t507 / 0.2e1 + t151 * mrSges(4,1) - t152 * mrSges(4,2) + t413 * t330 + t449 * t310 + t451 * t358 + t450 * t177 - t578) * t350 + (t350 * (Ifges(4,4) * t475 - Ifges(4,2) * t480 + t505) * t531 - t353 * (Ifges(4,5) * t475 - Ifges(4,6) * t480 + t504) / 0.2e1 + (Ifges(3,4) * t353 - Ifges(3,2) * t350) * t528 - t590) * t465) * t344 * qJD(2) - t374 * t145 + m(5) * (t120 * t197 + t13 * t78 + t14 * t77 + t174 * t220 - t374 * t60 + t375 * t59) - t375 * t146 + t98 * t38 + t3 * t96 + t4 * t97 + (Ifges(3,6) * t329 - Ifges(4,5) * t209 + t195 * mrSges(3,3) + t89 / 0.2e1 + t87 / 0.2e1 - t256 / 0.2e1 - t257 / 0.2e1 - t90 / 0.2e1 - t88 / 0.2e1 - t99 * mrSges(4,1) + Ifges(4,6) * t392 + t100 * mrSges(4,2) + t279 * Ifges(3,4) - t450 * t93 + t451 * t92 + t449 * t266 + t362) * t477 + t41 * (-mrSges(7,1) * t80 + mrSges(7,2) * t81) + (-mrSges(3,3) * t196 + Ifges(3,1) * t279 + Ifges(3,5) * t329) * t479 + (Ifges(7,5) * t81 + Ifges(7,6) * t80) * t545 + (Ifges(7,1) * t81 + Ifges(7,4) * t80) * t547 + t71 * t75 + t72 * t76 + t77 * t73 + t78 * t74 + t22 * t79 + t56 * t17 + t23 * t25 + t24 * t26 + (t442 + t623) * t346 / 0.2e1 + t288 * (mrSges(3,1) * t329 - mrSges(3,3) * t279) + t276 * t272 + t263 * t131 + t117 * t536 + (Ifges(7,4) * t81 + Ifges(7,2) * t80) * t549 + t50 * t110 + (-t527 * mrSges(2,1) - m(5) * t394 - t164 * mrSges(7,1) - t163 * mrSges(7,2) + (-mrSges(3,1) - t373) * t286 + (mrSges(2,2) + (-mrSges(3,3) - t409) * t344) * t351 + t391 * t214 + t412 * t215 + t366 * t285 + (-m(4) - m(3)) * t466 - t615 * (t215 * pkin(4) + t394)) * g(2) + t27 * t143 + t28 * t144 + t169 * t8 / 0.2e1 + (Ifges(4,1) * t280 - Ifges(4,4) * t382) * t539 + t81 * t562 + t80 * t564 + t183 * t167 + t184 * t166 + t197 * t39 + t162 * t207 + t161 * t208 + t220 * t109; (Ifges(7,1) * t192 + Ifges(7,4) * t191) * t548 + (Ifges(7,1) * t381 - Ifges(7,4) * t380) * t547 + (Ifges(4,2) * t345 + t511) * t613 + (t423 - t272) * t274 + t604 * t97 + (t1 * t112 + t111 * t2 + t15 * t604 + t16 * t605 + t190 * t6 + t41 * t603) * m(7) + t605 * t96 + t629 * (-t281 / 0.2e1 - t227 / 0.2e1) + (t400 * t561 + t401 * t566 + t403 * t567 - t407 * t6 + t427 * t54 + t637) * t292 + (Ifges(7,4) * t381 - Ifges(7,2) * t380) * t549 + (Ifges(7,4) * t192 + Ifges(7,2) * t191) * t550 + (-t152 * (-mrSges(4,2) * t350 - mrSges(4,3) * t480) - t151 * (mrSges(4,1) * t350 - mrSges(4,3) * t475)) * t465 + t634 * t281 + t601 * t144 + (-t10 * t225 + t11 * t224 + t20 * t201 + t57 * t601 + t58 * t602 + t598 * t67) * m(6) + t602 * t143 + t603 * t79 + (Ifges(5,4) * t543 - Ifges(6,6) * t544 - t622 + t648) * t227 + t645 * t282 + (t318 + t229 + t476) * t414 + (-t460 - t187) * t208 + t635 * t293 + (Ifges(7,5) * t381 - Ifges(7,6) * t380) * t545 + (Ifges(7,5) * t192 + Ifges(7,6) * t191) * t546 + (-m(5) * t297 + t287 - t615 * (t340 * t333 + t339 * t446 + t297) + ((-mrSges(7,1) * t474 - mrSges(7,2) * t470) * t339 + (-t518 * t644 + t576) * t350 + (t340 * t588 + t587) * t353) * t344) * g(3) - t628 * (Ifges(4,5) * t343 + t506) / 0.2e1 + (t166 * t345 - t167 * t343) * qJ(3) + t599 * t145 + t600 * t146 + (-t120 * t338 + t13 * t225 - t14 * t224 - t174 * t219 - t59 * t600 + t599 * t60) * m(5) + (-m(4) * t221 + t424 - t596) * t275 + t598 * t110 + (-pkin(2) * t181 + (-t151 * t343 + t152 * t345) * qJD(3) + t593 * qJ(3) - t151 * t187 - t152 * t188) * m(4) + t593 * mrSges(4,3) + (mrSges(7,1) * t595 - mrSges(7,2) * t594) * t41 + (t1 * t484 + t15 * t594 - t16 * t595 - t2 * t485) * mrSges(7,3) + (t369 * (t353 * t402 + t505) + t350 * t228) * t465 / 0.2e1 + (qJD(6) * t55 + t8) * t484 / 0.2e1 - t378 * t435 + t442 + ((t504 + (Ifges(4,5) * t345 - t500) * t353) * t528 + t590) * qJD(1) ^ 2 * t572 + t181 * t410 + (t74 - t75) * t225 + (t76 - t73) * t224 + (qJD(3) * t345 - t188) * t207 + t158 * t418 / 0.2e1 + ((t157 + t105 + t102) * t350 + t247 * (Ifges(4,5) * t350 + t353 * t404) + t330 * (Ifges(3,5) * t353 - Ifges(3,6) * t350)) * t429 + (-Ifges(5,4) * t542 - Ifges(5,2) * t543 + Ifges(6,6) * t541 + t574 + t650) * t226 + t343 * t117 / 0.2e1 - t338 * t39 + t116 * t530 + t111 * t25 + t112 * t26 - pkin(2) * t131 + t630 + (Ifges(6,4) * t541 + Ifges(5,5) * t542 + Ifges(6,5) * t544 - Ifges(3,2) * t414 + Ifges(5,6) * t543 + (Ifges(6,1) + Ifges(5,3)) * t532 + t578) * t436 + (-m(6) * t421 - m(7) * (-pkin(10) * t487 + t421) + mrSges(7,3) * t487 - m(5) * t469 + t575 * t284 + t579 * t283) * g(2) + (-m(6) * t420 - m(7) * (-pkin(10) * t486 + t420) + mrSges(7,3) * t486 - m(5) * t468 + t575 * t286 + t579 * t285) * g(1) + (Ifges(4,1) * t343 + t510) * t539 + t488 * t562 + t192 * t563 + t471 * t564 + t191 * t565 + t485 * t571 + t190 * t17 - t195 * mrSges(3,2) + t201 * t38 - t219 * t109; -(-t145 - t494) * t177 + (t395 + t597) * t358 + t395 * qJD(6) + t131 + t519 * t93 - t86 + t85 + t352 * t26 + t369 * t207 - t348 * t25 + t247 * t208 + (t1 * t352 + t177 * t41 - t2 * t348 - t173 * (t15 * t352 + t16 * t348)) * m(7) + (-t177 * t58 - t358 * t57 + t20) * m(6) + (t177 * t60 - t358 * t59 + t120) * m(5) + (t151 * t247 + t152 * t369 + t181) * m(4) + (-g(1) * t285 - g(2) * t283 + g(3) * t477) * (m(4) + t644); (-t502 + t101) * t542 + (t501 + t104) * t541 + (t171 + t103) * t544 + t494 * qJD(5) + t386 * t79 + (-t16 * t513 + t400 * t546 + t401 * t550 + t403 * t548 + t574 + t586) * t358 + t586 * qJD(6) + (-t16 * t456 - t522 + (t457 + t490) * t15 + t580) * mrSges(7,3) + (t210 * t627 + t211 * t581) * g(2) + (t214 * t627 + t215 * t581) * g(1) + (t581 * (t339 * t346 + t340 * t479) + t627 * t258) * g(3) + (-t75 + t17) * qJ(5) + t622 * t177 + (-m(6) * t57 + t514 + t597) * t60 + t90 + (-pkin(4) * t11 - qJ(5) * t10 - qJD(5) * t58 - t108 * t67) * m(6) - (t138 * t401 + t139 * t403 + t173 * t400) * qJD(6) / 0.2e1 + t257 + t256 + t57 * t517 + (t6 * qJ(5) - t15 * t18 - t16 * t19 + t41 * t624) * m(7) + (-m(6) * t58 - t143 + t145 + t515) * t59 + t55 * t427 + t6 * t406 - t87 - t58 * t516 - t2 * t513 - t19 * t96 - t18 * t97 - t89 + t88 - t362 + (-m(7) * t360 - t640) * t555 - pkin(4) * t76 - t348 * t8 / 0.2e1 - t108 * t110 + (-Ifges(5,2) * t358 - t172 + t629) * t543 + (Ifges(7,5) * t352 - Ifges(7,6) * t348) * t561 + t490 * t563 + (-Ifges(7,2) * t348 + t508) * t566 + (Ifges(7,1) * t352 - t509) * t567 + t352 * t571; t494 * t310 + (t110 + t396) * t358 + t76 + (t310 * t41 - t358 * t397 + t360 - t580) * m(7) + (-t310 * t58 + t358 * t67 + t11 - t580) * m(6) + t640; -t41 * (mrSges(7,1) * t139 + mrSges(7,2) * t138) + (Ifges(7,1) * t138 - t503) * t548 + t54 * t547 + (Ifges(7,5) * t138 - Ifges(7,6) * t139) * t546 - t15 * t96 + t16 * t97 - g(1) * (mrSges(7,1) * t163 - mrSges(7,2) * t164) - g(2) * ((t210 * t352 - t283 * t348) * mrSges(7,1) + (-t210 * t348 - t283 * t352) * mrSges(7,2)) - g(3) * ((t258 * t352 + t313) * mrSges(7,1) + (-t258 * t348 + t444) * mrSges(7,2)) + (t138 * t15 + t139 * t16) * mrSges(7,3) + t7 + (-Ifges(7,2) * t139 + t132 + t55) * t550 + t584;];
tau  = t9;
