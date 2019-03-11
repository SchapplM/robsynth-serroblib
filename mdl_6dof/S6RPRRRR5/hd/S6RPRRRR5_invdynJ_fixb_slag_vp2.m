% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:13
% EndTime: 2019-03-09 07:08:54
% DurationCPUTime: 26.19s
% Computational Cost: add. (25297->803), mult. (61767->1043), div. (0->0), fcn. (49356->18), ass. (0->382)
t339 = sin(qJ(5));
t421 = qJD(5) * t339;
t335 = sin(pkin(11));
t341 = sin(qJ(3));
t336 = cos(pkin(11));
t346 = cos(qJ(3));
t439 = t336 * t346;
t277 = -t335 * t341 + t439;
t264 = t277 * qJD(1);
t278 = t335 * t346 + t336 * t341;
t265 = t278 * qJD(1);
t340 = sin(qJ(4));
t345 = cos(qJ(4));
t548 = t345 * t264 - t340 * t265;
t578 = t548 * t339;
t600 = t421 - t578;
t595 = -mrSges(6,3) - mrSges(7,3);
t422 = qJD(4) * t345;
t410 = pkin(3) * t422;
t483 = pkin(7) + qJ(2);
t300 = t483 * t335;
t281 = qJD(1) * t300;
t301 = t483 * t336;
t282 = qJD(1) * t301;
t224 = -t281 * t341 + t282 * t346;
t195 = pkin(8) * t264 + t224;
t189 = t340 * t195;
t448 = t282 * t341;
t223 = -t346 * t281 - t448;
t194 = -pkin(8) * t265 + t223;
t141 = t194 * t345 - t189;
t368 = t264 * t340 + t345 * t265;
t171 = pkin(4) * t368 - pkin(9) * t548;
t498 = pkin(3) * t265;
t150 = t171 + t498;
t344 = cos(qJ(5));
t83 = -t141 * t339 + t344 * t150;
t599 = -t339 * t410 - t83;
t84 = t344 * t141 + t339 * t150;
t598 = t344 * t410 - t84;
t418 = qJD(1) * qJD(2);
t308 = qJ(2) * qJDD(1) + t418;
t532 = m(7) * pkin(5);
t333 = qJD(3) + qJD(4);
t196 = t333 * t344 - t339 * t368;
t208 = qJD(5) - t548;
t197 = t333 * t339 + t344 * t368;
t462 = t197 * Ifges(6,4);
t109 = t196 * Ifges(6,2) + t208 * Ifges(6,6) + t462;
t597 = -t109 / 0.2e1;
t596 = pkin(10) * t578;
t495 = pkin(3) * t340;
t317 = pkin(9) + t495;
t482 = -pkin(10) - t317;
t393 = qJD(5) * t482;
t594 = t339 * t393 + t596 + t598;
t577 = t548 * t344;
t588 = pkin(5) * t368 - pkin(10) * t577;
t593 = t344 * t393 - t588 + t599;
t348 = -pkin(10) - pkin(9);
t406 = qJD(5) * t348;
t191 = qJD(3) * pkin(3) + t194;
t136 = t191 * t345 - t189;
t86 = t344 * t136 + t339 * t171;
t592 = t339 * t406 + t596 - t86;
t85 = -t136 * t339 + t344 * t171;
t591 = t344 * t406 - t588 - t85;
t425 = t335 ^ 2 + t336 ^ 2;
t332 = pkin(11) + qJ(3);
t324 = qJ(4) + t332;
t314 = sin(t324);
t334 = qJ(5) + qJ(6);
t325 = sin(t334);
t477 = mrSges(7,2) * t325;
t478 = mrSges(6,2) * t339;
t590 = (-t477 - t478) * t314;
t589 = t600 * pkin(5);
t129 = -pkin(4) * t333 - t136;
t100 = -pkin(5) * t196 + t129;
t266 = t277 * qJD(3);
t219 = qJD(1) * t266 + qJDD(1) * t278;
t267 = t278 * qJD(3);
t220 = -qJD(1) * t267 + qJDD(1) * t277;
t127 = qJD(4) * t548 + t219 * t345 + t220 * t340;
t423 = qJD(4) * t340;
t128 = -t340 * t219 + t220 * t345 - t264 * t423 - t265 * t422;
t343 = cos(qJ(6));
t338 = sin(qJ(6));
t190 = t345 * t195;
t137 = t191 * t340 + t190;
t130 = pkin(9) * t333 + t137;
t316 = pkin(2) * t336 + pkin(1);
t292 = -qJD(1) * t316 + qJD(2);
t227 = -pkin(3) * t264 + t292;
t138 = -pkin(4) * t548 - pkin(9) * t368 + t227;
t76 = t130 * t344 + t138 * t339;
t67 = pkin(10) * t196 + t76;
t460 = t338 * t67;
t75 = -t130 * t339 + t344 * t138;
t66 = -pkin(10) * t197 + t75;
t55 = pkin(5) * t208 + t66;
t20 = t343 * t55 - t460;
t458 = t343 * t67;
t21 = t338 * t55 + t458;
t364 = t338 * t339 - t343 * t344;
t537 = qJD(5) + qJD(6);
t225 = t537 * t364;
t285 = t338 * t344 + t339 * t343;
t226 = t537 * t285;
t420 = qJD(5) * t344;
t328 = qJDD(3) + qJDD(4);
t388 = pkin(7) * qJDD(1) + t308;
t258 = t388 * t335;
t259 = t388 * t336;
t175 = -qJD(3) * t224 - t346 * t258 - t259 * t341;
t126 = qJDD(3) * pkin(3) - pkin(8) * t219 + t175;
t424 = qJD(3) * t346;
t174 = -qJD(3) * t448 - t341 * t258 + t346 * t259 - t281 * t424;
t133 = pkin(8) * t220 + t174;
t52 = t340 * t126 + t345 * t133 + t191 * t422 - t195 * t423;
t49 = pkin(9) * t328 + t52;
t291 = -qJDD(1) * t316 + qJDD(2);
t198 = -pkin(3) * t220 + t291;
t65 = -pkin(4) * t128 - pkin(9) * t127 + t198;
t15 = -t130 * t421 + t138 * t420 + t339 * t65 + t344 * t49;
t94 = -qJD(5) * t197 - t127 * t339 + t328 * t344;
t12 = pkin(10) * t94 + t15;
t123 = qJDD(5) - t128;
t16 = -qJD(5) * t76 - t339 * t49 + t344 * t65;
t93 = qJD(5) * t196 + t127 * t344 + t328 * t339;
t7 = pkin(5) * t123 - pkin(10) * t93 + t16;
t3 = qJD(6) * t20 + t12 * t343 + t338 * t7;
t53 = t126 * t345 - t340 * t133 - t191 * t423 - t195 * t422;
t50 = -pkin(4) * t328 - t53;
t31 = -pkin(5) * t94 + t50;
t375 = mrSges(6,1) * t339 + mrSges(6,2) * t344;
t358 = t129 * t375;
t371 = Ifges(6,5) * t344 - Ifges(6,6) * t339;
t469 = Ifges(6,4) * t344;
t372 = -Ifges(6,2) * t339 + t469;
t470 = Ifges(6,4) * t339;
t373 = Ifges(6,1) * t344 - t470;
t395 = -t421 / 0.2e1;
t193 = Ifges(6,4) * t196;
t110 = t197 * Ifges(6,1) + t208 * Ifges(6,5) + t193;
t434 = t344 * t110;
t398 = t434 / 0.2e1;
t4 = -qJD(6) * t21 - t12 * t338 + t343 * t7;
t40 = t93 * Ifges(6,4) + t94 * Ifges(6,2) + t123 * Ifges(6,6);
t206 = qJD(6) + t208;
t509 = t206 / 0.2e1;
t510 = -t206 / 0.2e1;
t146 = t196 * t338 + t197 * t343;
t515 = t146 / 0.2e1;
t516 = -t146 / 0.2e1;
t387 = t343 * t196 - t197 * t338;
t517 = t387 / 0.2e1;
t518 = -t387 / 0.2e1;
t519 = t123 / 0.2e1;
t119 = qJDD(6) + t123;
t520 = t119 / 0.2e1;
t521 = t94 / 0.2e1;
t522 = t93 / 0.2e1;
t139 = Ifges(7,4) * t387;
t72 = Ifges(7,1) * t146 + Ifges(7,5) * t206 + t139;
t523 = t72 / 0.2e1;
t524 = -t72 / 0.2e1;
t468 = Ifges(7,4) * t146;
t71 = Ifges(7,2) * t387 + Ifges(7,6) * t206 + t468;
t525 = t71 / 0.2e1;
t526 = -t71 / 0.2e1;
t45 = -qJD(6) * t146 - t338 * t93 + t343 * t94;
t527 = t45 / 0.2e1;
t44 = qJD(6) * t387 + t338 * t94 + t343 * t93;
t528 = t44 / 0.2e1;
t529 = Ifges(6,1) * t522 + Ifges(6,4) * t521 + Ifges(6,5) * t519;
t530 = Ifges(7,1) * t528 + Ifges(7,4) * t527 + Ifges(7,5) * t520;
t531 = Ifges(7,4) * t528 + Ifges(7,2) * t527 + Ifges(7,6) * t520;
t568 = mrSges(7,2) * t100 - t20 * mrSges(7,3);
t569 = -mrSges(7,1) * t100 + t21 * mrSges(7,3);
t481 = mrSges(6,1) * t344;
t573 = t478 - t481;
t574 = t15 * t344 - t16 * t339;
t579 = t364 * t548;
t580 = t285 * t548;
t587 = t31 * (mrSges(7,1) * t364 + mrSges(7,2) * t285) - t580 * t526 + (-Ifges(7,1) * t579 - Ifges(7,4) * t580) * t516 + (-Ifges(7,4) * t579 - Ifges(7,2) * t580) * t518 + t50 * t573 + (-Ifges(7,5) * t579 - Ifges(7,6) * t580) * t510 - t100 * (mrSges(7,1) * t580 - mrSges(7,2) * t579) - t579 * t524 + t339 * t529 + t285 * t530 + (t196 * t372 + t197 * t373 + t208 * t371) * qJD(5) / 0.2e1 - t364 * t531 + (Ifges(7,1) * t285 - Ifges(7,4) * t364) * t528 + (Ifges(7,5) * t285 - Ifges(7,6) * t364) * t520 + (Ifges(7,4) * t285 - Ifges(7,2) * t364) * t527 + (Ifges(6,5) * t339 + Ifges(6,6) * t344) * t519 + (Ifges(6,2) * t344 + t470) * t521 + (Ifges(6,1) * t339 + t469) * t522 + (t358 + t398) * qJD(5) - (Ifges(7,1) * t515 + Ifges(7,4) * t517 + Ifges(7,5) * t509 + t523 + t568) * t225 - (Ifges(7,4) * t515 + Ifges(7,2) * t517 + Ifges(7,6) * t509 + t525 + t569) * t226 + t344 * t40 / 0.2e1 + Ifges(5,3) * t328 + Ifges(5,5) * t127 + Ifges(5,6) * t128 + (t574 - t600 * t76 + (-t420 + t577) * t75) * mrSges(6,3) - t52 * mrSges(5,2) + t53 * mrSges(5,1) + t109 * t395 + (-t20 * t579 + t21 * t580 - t285 * t4 - t3 * t364) * mrSges(7,3);
t586 = t339 * t532;
t207 = Ifges(5,4) * t548;
t585 = t196 * Ifges(6,6);
t584 = t208 * Ifges(6,3);
t583 = t333 * Ifges(5,5);
t582 = t333 * Ifges(5,6);
t581 = t548 * Ifges(5,2);
t557 = t197 * Ifges(6,5) + t146 * Ifges(7,5) + Ifges(7,6) * t387 + t206 * Ifges(7,3) + t584 + t585;
t140 = t194 * t340 + t190;
t576 = -pkin(3) * t423 + t140;
t552 = mrSges(5,1) * t333 + mrSges(6,1) * t196 - mrSges(6,2) * t197 - mrSges(5,3) * t368;
t326 = cos(t334);
t480 = mrSges(7,1) * t326;
t575 = (t480 + t481) * t314;
t572 = -m(7) - m(6) - m(5);
t315 = cos(t324);
t571 = -t315 * mrSges(5,1) + (mrSges(5,2) + t595) * t314;
t397 = m(3) * qJ(2) + mrSges(3,3);
t567 = -m(4) * t483 + mrSges(2,2) - mrSges(4,3) - mrSges(5,3) - t397;
t566 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t544 = t315 * pkin(4) + t314 * pkin(9);
t318 = pkin(5) * t344 + pkin(4);
t547 = -t314 * t348 + t315 * t318;
t565 = -m(6) * t544 - m(7) * t547;
t564 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t563 = t227 * mrSges(5,2) - t136 * mrSges(5,3) + t339 * t597;
t322 = sin(t332);
t323 = cos(t332);
t379 = mrSges(4,1) * t323 - mrSges(4,2) * t322;
t380 = -mrSges(3,1) * t336 + mrSges(3,2) * t335;
t562 = m(3) * pkin(1) + m(4) * t316 + mrSges(2,1) + t379 - t380 - t571;
t561 = mrSges(5,1) * t227 + mrSges(6,1) * t75 + mrSges(7,1) * t20 - mrSges(6,2) * t76 - mrSges(7,2) * t21 - t137 * mrSges(5,3);
t505 = -t368 / 0.2e1;
t558 = mrSges(6,1) + t532;
t272 = t482 * t339;
t327 = t344 * pkin(10);
t444 = t317 * t344;
t273 = t327 + t444;
t217 = t272 * t343 - t273 * t338;
t556 = qJD(6) * t217 + t338 * t593 + t343 * t594;
t218 = t272 * t338 + t273 * t343;
t555 = -qJD(6) * t218 - t338 * t594 + t343 * t593;
t303 = t348 * t339;
t490 = pkin(9) * t344;
t304 = t327 + t490;
t230 = t303 * t343 - t304 * t338;
t554 = qJD(6) * t230 + t338 * t591 + t343 * t592;
t231 = t303 * t338 + t304 * t343;
t553 = -qJD(6) * t231 - t338 * t592 + t343 * t591;
t222 = t277 * t340 + t278 * t345;
t173 = t364 * t222;
t228 = -t346 * t300 - t301 * t341;
t204 = -pkin(8) * t278 + t228;
t229 = -t341 * t300 + t346 * t301;
t205 = pkin(8) * t277 + t229;
t160 = t204 * t340 + t205 * t345;
t152 = t344 * t160;
t247 = -pkin(3) * t277 - t316;
t367 = t345 * t277 - t278 * t340;
t161 = -pkin(4) * t367 - pkin(9) * t222 + t247;
t90 = t339 * t161 + t152;
t551 = t345 * t204 - t205 * t340;
t550 = -t225 + t579;
t549 = -t226 + t580;
t546 = -t576 + t589;
t365 = -t314 * t318 - t315 * t348;
t493 = pkin(4) * t314;
t496 = pkin(3) * t322;
t543 = -m(7) * (t365 - t496) - m(6) * (-t493 - t496) + t575;
t542 = -m(7) * t365 + t575;
t541 = -t137 + t589;
t347 = cos(qJ(1));
t445 = t315 * t347;
t540 = t347 * t590 + t445 * t595;
t342 = sin(qJ(1));
t446 = t315 * t342;
t539 = t342 * t590 + t446 * t595;
t538 = g(1) * t347 + g(2) * t342;
t536 = t571 + (t477 - t480 + t573) * t315;
t474 = mrSges(6,3) * t196;
t153 = -mrSges(6,2) * t208 + t474;
t473 = mrSges(6,3) * t197;
t154 = mrSges(6,1) * t208 - t473;
t58 = mrSges(6,1) * t123 - mrSges(6,3) * t93;
t535 = m(6) * ((-t339 * t76 - t344 * t75) * qJD(5) + t574) - t154 * t420 - t153 * t421 - t339 * t58;
t500 = -t333 / 0.2e1;
t508 = -t208 / 0.2e1;
t512 = -t197 / 0.2e1;
t513 = -t196 / 0.2e1;
t534 = Ifges(5,1) * t505 + Ifges(5,5) * t500 + t371 * t508 + t372 * t513 + t373 * t512 - t358 - t434 / 0.2e1 - t563;
t506 = -t548 / 0.2e1;
t533 = Ifges(6,5) * t512 + Ifges(7,5) * t516 - Ifges(5,2) * t506 - Ifges(5,6) * t500 + Ifges(6,6) * t513 + Ifges(7,6) * t518 + Ifges(6,3) * t508 + Ifges(7,3) * t510 - t561;
t511 = t197 / 0.2e1;
t504 = t368 / 0.2e1;
t501 = t265 / 0.2e1;
t497 = pkin(3) * t267;
t312 = pkin(3) * t323;
t494 = pkin(3) * t345;
t492 = pkin(5) * t197;
t487 = g(3) * t314;
t475 = mrSges(4,3) * t264;
t472 = Ifges(4,4) * t265;
t471 = Ifges(5,4) * t368;
t467 = pkin(5) * qJD(6);
t461 = t224 * mrSges(4,3);
t176 = qJD(4) * t367 + t266 * t345 - t267 * t340;
t457 = t176 * t344;
t450 = t222 * t339;
t449 = t222 * t344;
t443 = t325 * t342;
t442 = t325 * t347;
t441 = t326 * t342;
t440 = t326 * t347;
t437 = t339 * t342;
t436 = t339 * t347;
t435 = t342 * t344;
t433 = t344 * t347;
t238 = t315 * t443 + t440;
t239 = -t315 * t441 + t442;
t432 = -t238 * mrSges(7,1) + t239 * mrSges(7,2);
t240 = -t315 * t442 + t441;
t241 = t315 * t440 + t443;
t431 = t240 * mrSges(7,1) - t241 * mrSges(7,2);
t417 = qJDD(1) * t335;
t416 = qJDD(1) * t336;
t415 = Ifges(7,5) * t44 + Ifges(7,6) * t45 + Ifges(7,3) * t119;
t414 = Ifges(6,5) * t93 + Ifges(6,6) * t94 + Ifges(6,3) * t123;
t401 = t222 * t420;
t201 = -t300 * t424 + qJD(2) * t439 + (-qJD(2) * t335 - qJD(3) * t301) * t341;
t181 = -pkin(8) * t267 + t201;
t202 = -t278 * qJD(2) - qJD(3) * t229;
t182 = -pkin(8) * t266 + t202;
t80 = qJD(4) * t551 + t181 * t345 + t182 * t340;
t177 = qJD(4) * t222 + t266 * t340 + t345 * t267;
t99 = pkin(4) * t177 - pkin(9) * t176 + t497;
t394 = -t339 * t80 + t344 * t99;
t392 = -t220 * mrSges(4,1) + t219 * mrSges(4,2);
t391 = -t128 * mrSges(5,1) + t127 * mrSges(5,2);
t89 = -t160 * t339 + t344 * t161;
t381 = -mrSges(3,1) * t416 + mrSges(3,2) * t417;
t376 = mrSges(5,1) * t314 + mrSges(5,2) * t315;
t374 = -mrSges(7,1) * t325 - mrSges(7,2) * t326;
t69 = -pkin(5) * t367 - pkin(10) * t449 + t89;
t79 = -pkin(10) * t450 + t90;
t36 = -t338 * t79 + t343 * t69;
t37 = t338 * t69 + t343 * t79;
t370 = -t339 * t75 + t344 * t76;
t369 = t153 * t344 - t154 * t339;
t361 = t415 + t566;
t256 = -t315 * t436 + t435;
t254 = t315 * t437 + t433;
t360 = t176 * t339 + t401;
t359 = t222 * t421 - t457;
t24 = -t160 * t421 + t161 * t420 + t339 * t99 + t344 * t80;
t81 = qJD(4) * t160 + t181 * t340 - t345 * t182;
t329 = -pkin(8) - t483;
t321 = -qJDD(1) * pkin(1) + qJDD(2);
t319 = -pkin(4) - t494;
t302 = -t318 - t494;
t299 = pkin(9) * t445;
t298 = pkin(9) * t446;
t283 = t312 + t316;
t260 = Ifges(4,4) * t264;
t257 = t315 * t433 + t437;
t255 = -t315 * t435 + t436;
t243 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t265;
t242 = -qJD(3) * mrSges(4,2) + t475;
t210 = t265 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t260;
t209 = t264 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t472;
t199 = -mrSges(5,2) * t333 + mrSges(5,3) * t548;
t172 = t285 * t222;
t170 = -mrSges(5,1) * t548 + mrSges(5,2) * t368;
t167 = Ifges(5,1) * t368 + t207 + t583;
t166 = t471 + t581 + t582;
t115 = -mrSges(5,2) * t328 + mrSges(5,3) * t128;
t114 = mrSges(5,1) * t328 - mrSges(5,3) * t127;
t112 = pkin(5) * t450 - t551;
t103 = mrSges(7,1) * t206 - mrSges(7,3) * t146;
t102 = -mrSges(7,2) * t206 + mrSges(7,3) * t387;
t82 = -mrSges(7,1) * t387 + mrSges(7,2) * t146;
t59 = -mrSges(6,2) * t123 + mrSges(6,3) * t94;
t57 = t173 * t537 - t285 * t176;
t56 = -t176 * t364 - t222 * t226;
t54 = pkin(5) * t360 + t81;
t47 = -mrSges(6,1) * t94 + mrSges(6,2) * t93;
t27 = t343 * t66 - t460;
t26 = -t338 * t66 - t458;
t25 = -qJD(5) * t90 + t394;
t23 = -mrSges(7,2) * t119 + mrSges(7,3) * t45;
t22 = mrSges(7,1) * t119 - mrSges(7,3) * t44;
t18 = -pkin(10) * t360 + t24;
t17 = -pkin(10) * t457 + pkin(5) * t177 + (-t152 + (pkin(10) * t222 - t161) * t339) * qJD(5) + t394;
t13 = -mrSges(7,1) * t45 + mrSges(7,2) * t44;
t6 = -qJD(6) * t37 + t17 * t343 - t18 * t338;
t5 = qJD(6) * t36 + t17 * t338 + t18 * t343;
t1 = [(-t437 * t532 - t257 * mrSges(6,1) - t241 * mrSges(7,1) - t256 * mrSges(6,2) - t240 * mrSges(7,2) + t572 * (t347 * t283 - t329 * t342) + t567 * t342 + (-t562 + t565) * t347) * g(2) + (-m(5) * t136 + m(6) * t129 - t552) * t81 - (-m(5) * t53 + m(6) * t50 - t114 + t47) * t551 + m(5) * (t137 * t80 + t160 * t52 + t198 * t247 + t227 * t497) + (-t172 * t3 + t173 * t4 - t20 * t56 + t21 * t57) * mrSges(7,3) + t449 * t529 - t173 * t530 - t172 * t531 + (Ifges(7,1) * t56 + Ifges(7,4) * t57) * t515 + (-Ifges(7,1) * t173 - Ifges(7,4) * t172) * t528 - (t415 + t414) * t367 / 0.2e1 + 0.2e1 * t425 * t308 * mrSges(3,3) + m(3) * (-pkin(1) * t321 + (t308 + t418) * qJ(2) * t425) + t196 * (-Ifges(6,4) * t359 - Ifges(6,2) * t360) / 0.2e1 + (t198 * mrSges(5,2) - t53 * mrSges(5,3) + Ifges(5,1) * t127 + Ifges(5,4) * t128 + Ifges(5,5) * t328 + t110 * t395 + t371 * t519 + t372 * t521 + t373 * t522 + t375 * t50) * t222 + (Ifges(3,4) * t335 + Ifges(3,2) * t336) * t416 + (Ifges(3,1) * t335 + Ifges(3,4) * t336) * t417 + m(6) * (t15 * t90 + t16 * t89 + t24 * t76 + t25 * t75) + (Ifges(7,4) * t56 + Ifges(7,2) * t57) * t517 + (t557 / 0.2e1 + Ifges(7,3) * t509 + Ifges(6,5) * t511 + Ifges(7,5) * t515 + Ifges(7,6) * t517 - Ifges(5,4) * t504 - t582 / 0.2e1 - t581 / 0.2e1 + t585 / 0.2e1 + t584 / 0.2e1 - t166 / 0.2e1 + t561) * t177 + (-Ifges(7,4) * t173 - Ifges(7,2) * t172) * t527 - t40 * t450 / 0.2e1 + (-t15 * t450 - t16 * t449 + t359 * t75 - t360 * t76) * mrSges(6,3) + (Ifges(7,5) * t56 + Ifges(7,6) * t57) * t509 + (-mrSges(5,1) * t198 + mrSges(5,3) * t52 + Ifges(5,4) * t127 - Ifges(6,5) * t522 - Ifges(7,5) * t528 + Ifges(5,2) * t128 + Ifges(5,6) * t328 - Ifges(6,6) * t521 - Ifges(7,6) * t527 - Ifges(6,3) * t519 - Ifges(7,3) * t520 - t564 - t566) * t367 + t56 * t523 + t57 * t525 + (Ifges(5,1) * t504 + t583 / 0.2e1 + t207 / 0.2e1 + t167 / 0.2e1 + t398 + t563) * t176 + t208 * (-Ifges(6,5) * t359 - Ifges(6,6) * t360) / 0.2e1 + (-Ifges(6,1) * t359 - Ifges(6,4) * t360) * t511 - t267 * t461 + t401 * t597 + Ifges(2,3) * qJDD(1) + (mrSges(4,2) * t291 - mrSges(4,3) * t175 + Ifges(4,1) * t219 + Ifges(4,4) * t220 + Ifges(4,5) * qJDD(3)) * t278 + t247 * t391 - t316 * t392 + (-mrSges(4,1) * t291 + mrSges(4,3) * t174 + Ifges(4,4) * t219 + Ifges(4,2) * t220 + Ifges(4,6) * qJDD(3)) * t277 - t223 * t266 * mrSges(4,3) + t266 * t210 / 0.2e1 - t267 * t209 / 0.2e1 - pkin(1) * t381 + t321 * t380 + t201 * t242 + t202 * t243 + t228 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t219) + t229 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t220) + t80 * t199 + t129 * (mrSges(6,1) * t360 - mrSges(6,2) * t359) + t160 * t115 + t24 * t153 + t25 * t154 + (-Ifges(7,5) * t173 - Ifges(7,6) * t172) * t520 + t31 * (mrSges(7,1) * t172 - mrSges(7,2) * t173) + (Ifges(4,1) * t266 - Ifges(4,4) * t267) * t501 + t292 * (mrSges(4,1) * t267 + mrSges(4,2) * t266) + t264 * (Ifges(4,4) * t266 - Ifges(4,2) * t267) / 0.2e1 + qJD(3) * (Ifges(4,5) * t266 - Ifges(4,6) * t267) / 0.2e1 + t36 * t22 + t37 * t23 + (-t255 * mrSges(6,1) - t239 * mrSges(7,1) - t254 * mrSges(6,2) - t238 * mrSges(7,2) + (-t329 * t572 + t567 - t586) * t347 + (-m(6) * (-t283 - t544) + m(5) * t283 - m(7) * (-t283 - t547) + t562) * t342) * g(1) + t54 * t82 + t89 * t58 + t90 * t59 + m(4) * (t174 * t229 + t175 * t228 + t201 * t224 + t202 * t223 - t291 * t316) + m(7) * (t100 * t54 + t112 * t31 + t20 * t6 + t21 * t5 + t3 * t37 + t36 * t4) + t100 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t5 * t102 + t6 * t103 + t112 * t13 + t170 * t497; t391 + t392 + (-t82 + t552) * t368 + t549 * t103 + t550 * t102 + m(3) * t321 + t381 + t344 * t58 + t339 * t59 - t364 * t22 + t285 * t23 + t265 * t243 - t264 * t242 + (-t199 - t369) * t548 + t369 * qJD(5) + (-g(1) * t342 + g(2) * t347) * (m(3) + m(4) - t572) - t397 * t425 * qJD(1) ^ 2 + (-t100 * t368 + t20 * t549 + t21 * t550 + t285 * t3 - t364 * t4) * m(7) + (-t129 * t368 + t15 * t339 + t16 * t344 + t208 * t370) * m(6) + (t136 * t368 - t137 * t548 + t198) * m(5) + (t223 * t265 - t224 * t264 + t291) * m(4); (-m(6) * (t312 + t544) - m(7) * (t312 + t547) - t379 - m(5) * t312 + t536) * g(3) + t587 + t598 * t153 + t599 * t154 + t576 * t552 - (-Ifges(4,2) * t265 + t210 + t260) * t264 / 0.2e1 + ((t340 * t52 + t345 * t53 + (-t136 * t340 + t137 * t345) * qJD(4)) * pkin(3) + t136 * t140 - t137 * t141 - t227 * t498) * m(5) + (Ifges(5,4) * t506 - t167 / 0.2e1 + t534) * t548 + (t475 - t242) * t223 + (t319 * t50 + (t129 * t340 + t345 * t370) * qJD(4) * pkin(3) - t129 * t140 - t75 * t83 - t76 * t84 - g(1) * t299 - g(2) * t298) * m(6) + (m(5) * t496 + mrSges(4,1) * t322 + mrSges(4,2) * t323 + t376) * t538 + (t410 - t141) * t199 + t114 * t494 + t115 * t495 + t209 * t501 + (t533 - Ifges(5,4) * t505 + t166 / 0.2e1) * t368 + t535 * t317 + t59 * t444 + t265 * t461 + Ifges(4,3) * qJDD(3) + t319 * t47 - t292 * (mrSges(4,1) * t265 + mrSges(4,2) * t264) + t302 * t13 - qJD(3) * (Ifges(4,5) * t264 - Ifges(4,6) * t265) / 0.2e1 + t224 * t243 + Ifges(4,5) * t219 + Ifges(4,6) * t220 + t217 * t22 + t218 * t23 - t174 * mrSges(4,2) + t175 * mrSges(4,1) + t557 * t505 + (t342 * t543 + t539) * g(2) + (t347 * t543 + t540) * g(1) + t546 * t82 + t555 * t103 + t556 * t102 + (t100 * t546 + t20 * t555 + t21 * t556 + t217 * t4 + t218 * t3 + t302 * t31) * m(7) - t265 * (Ifges(4,1) * t264 - t472) / 0.2e1 - t170 * t498; t587 + t534 * t548 + (-t129 * t137 - t75 * t85 - t76 * t86 - pkin(4) * t50 - g(1) * (-t347 * t493 + t299) - g(2) * (-t342 * t493 + t298)) * m(6) + t59 * t490 + t166 * t504 + t533 * t368 + (t536 + t565) * g(3) + t535 * pkin(9) + (t167 + t207) * t506 - t318 * t13 + t230 * t22 + t231 * t23 - t136 * t199 - t86 * t153 - t85 * t154 + t538 * t376 + t541 * t82 + (t342 * t542 + t539) * g(2) + (t347 * t542 + t540) * g(1) + t552 * t137 - pkin(4) * t47 + t553 * t103 + t554 * t102 + (t100 * t541 + t20 * t553 + t21 * t554 + t230 * t4 + t231 * t3 - t31 * t318) * m(7) + (-t471 + t557) * t505; t564 + (t3 * t338 + t343 * t4 + (-t20 * t338 + t21 * t343) * qJD(6)) * t532 + (t343 * t467 - t27) * t102 + t414 + (Ifges(6,5) * t196 - Ifges(6,6) * t197) * t508 + t109 * t511 + (Ifges(6,1) * t196 - t462) * t512 + (-t338 * t467 - t26) * t103 + t361 - (Ifges(7,4) * t516 + Ifges(7,2) * t518 + Ifges(7,6) * t510 + t526 - t569) * t146 + (Ifges(7,1) * t516 + Ifges(7,4) * t518 + Ifges(7,5) * t510 + t524 - t568) * t387 - t129 * (mrSges(6,1) * t197 + mrSges(6,2) * t196) + (t473 + t154) * t76 + (t474 - t153) * t75 + (-t374 + t375 + t586) * t487 + (-mrSges(6,2) * t255 + t254 * t558 - t432) * g(2) + (mrSges(6,2) * t257 - t256 * t558 - t431) * g(1) + (t22 * t343 + t23 * t338) * pkin(5) + (-Ifges(6,2) * t197 + t110 + t193) * t513 - t82 * t492 - m(7) * (t100 * t492 + t20 * t26 + t21 * t27); -t100 * (mrSges(7,1) * t146 + mrSges(7,2) * t387) + (Ifges(7,1) * t387 - t468) * t516 + t71 * t515 + (Ifges(7,5) * t387 - Ifges(7,6) * t146) * t510 - t20 * t102 + t21 * t103 - g(1) * t431 - g(2) * t432 - t374 * t487 + (t146 * t21 + t20 * t387) * mrSges(7,3) + t361 + (-Ifges(7,2) * t146 + t139 + t72) * t518;];
tau  = t1;
