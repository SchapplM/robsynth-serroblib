% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:50
% EndTime: 2019-03-09 22:49:15
% DurationCPUTime: 44.65s
% Computational Cost: add. (31609->995), mult. (81174->1402), div. (0->0), fcn. (65156->12), ass. (0->425)
t404 = cos(qJ(2));
t401 = sin(qJ(2));
t395 = sin(pkin(6));
t461 = qJD(1) * t395;
t445 = t401 * t461;
t397 = cos(pkin(6));
t460 = qJD(1) * t397;
t452 = pkin(1) * t460;
t345 = -pkin(8) * t445 + t404 * t452;
t413 = (pkin(2) * t401 - pkin(9) * t404) * t395;
t346 = qJD(1) * t413;
t400 = sin(qJ(3));
t403 = cos(qJ(3));
t277 = -t400 * t345 + t403 * t346;
t567 = -pkin(10) - pkin(9);
t446 = qJD(3) * t567;
t468 = t403 * t404;
t624 = -(pkin(3) * t401 - pkin(10) * t468) * t461 - t277 + t403 * t446;
t278 = t403 * t345 + t400 * t346;
t444 = t404 * t461;
t431 = t400 * t444;
t623 = -pkin(10) * t431 - t400 * t446 + t278;
t399 = sin(qJ(4));
t524 = cos(qJ(4));
t411 = -t399 * t400 + t403 * t524;
t580 = qJD(3) + qJD(4);
t300 = t580 * t411;
t315 = t411 * t444;
t463 = t300 - t315;
t364 = t399 * t403 + t400 * t524;
t301 = t580 * t364;
t314 = t364 * t444;
t462 = t301 - t314;
t377 = t567 * t400;
t378 = t567 * t403;
t582 = t524 * t377 + t399 * t378;
t599 = qJD(4) * t582 + t399 * t624 - t623 * t524;
t348 = pkin(8) * t444 + t401 * t452;
t305 = pkin(3) * t431 + t348;
t457 = qJD(3) * t400;
t622 = pkin(3) * t457 - t305;
t621 = qJ(5) * t445 - t599;
t620 = t462 * pkin(4) - t463 * qJ(5) - qJD(5) * t364 + t622;
t394 = sin(pkin(12));
t396 = cos(pkin(12));
t280 = -t315 * t394 + t396 * t445;
t480 = t300 * t394;
t435 = t280 + t480;
t281 = t315 * t396 + t394 * t445;
t479 = t300 * t396;
t619 = t281 - t479;
t433 = qJD(2) + t460;
t408 = t400 * t445 - t403 * t433;
t459 = qJD(2) * t395;
t438 = qJD(1) * t459;
t429 = t404 * t438;
t294 = -qJD(3) * t408 + t403 * t429;
t329 = t400 * t433 + t403 * t445;
t295 = -qJD(3) * t329 - t400 * t429;
t406 = t329 * t524 - t399 * t408;
t167 = qJD(4) * t406 + t399 * t294 - t295 * t524;
t265 = t399 * t329 + t524 * t408;
t166 = -qJD(4) * t265 + t294 * t524 + t399 * t295;
t430 = t401 * t438;
t145 = -t166 * t394 + t396 * t430;
t146 = t166 * t396 + t394 * t430;
t376 = qJD(3) - t444;
t368 = qJD(4) + t376;
t235 = t368 * t394 + t396 * t406;
t398 = sin(qJ(6));
t402 = cos(qJ(6));
t436 = t396 * t368 - t394 * t406;
t151 = t235 * t402 + t398 * t436;
t53 = -qJD(6) * t151 + t145 * t402 - t146 * t398;
t50 = Ifges(7,6) * t53;
t595 = -t235 * t398 + t402 * t436;
t52 = qJD(6) * t595 + t145 * t398 + t146 * t402;
t51 = Ifges(7,5) * t52;
t10 = Ifges(7,3) * t167 + t50 + t51;
t313 = pkin(9) * t433 + t348;
t341 = (-pkin(2) * t404 - pkin(9) * t401 - pkin(1)) * t395;
t324 = qJD(1) * t341;
t256 = t403 * t313 + t400 * t324;
t347 = qJD(2) * t413;
t336 = qJD(1) * t347;
t475 = t395 * t401;
t384 = pkin(8) * t475;
t523 = pkin(1) * t404;
t356 = t397 * t523 - t384;
t349 = t356 * qJD(2);
t337 = qJD(1) * t349;
t198 = -qJD(3) * t256 + t403 * t336 - t337 * t400;
t135 = pkin(3) * t430 - pkin(10) * t294 + t198;
t456 = qJD(3) * t403;
t197 = -t313 * t457 + t324 * t456 + t400 * t336 + t403 * t337;
t152 = pkin(10) * t295 + t197;
t255 = -t313 * t400 + t403 * t324;
t225 = -pkin(10) * t329 + t255;
t207 = pkin(3) * t376 + t225;
t226 = -pkin(10) * t408 + t256;
t439 = qJD(4) * t524;
t455 = qJD(4) * t399;
t44 = t399 * t135 + t524 * t152 + t207 * t439 - t226 * t455;
t40 = qJ(5) * t430 + qJD(5) * t368 + t44;
t474 = t395 * t404;
t357 = t397 * t401 * pkin(1) + pkin(8) * t474;
t350 = t357 * qJD(2);
t338 = qJD(1) * t350;
t254 = -t295 * pkin(3) + t338;
t69 = t167 * pkin(4) - t166 * qJ(5) - qJD(5) * t406 + t254;
t19 = -t394 * t40 + t396 * t69;
t20 = t394 * t69 + t396 * t40;
t224 = t524 * t226;
t124 = t399 * t207 + t224;
t116 = t368 * qJ(5) + t124;
t312 = -pkin(2) * t433 - t345;
t268 = pkin(3) * t408 + t312;
t136 = t265 * pkin(4) - qJ(5) * t406 + t268;
t76 = -t116 * t394 + t396 * t136;
t48 = pkin(5) * t265 - pkin(11) * t235 + t76;
t77 = t396 * t116 + t394 * t136;
t59 = pkin(11) * t436 + t77;
t15 = -t398 * t59 + t402 * t48;
t6 = pkin(5) * t167 - pkin(11) * t146 + t19;
t9 = pkin(11) * t145 + t20;
t2 = qJD(6) * t15 + t398 * t6 + t402 * t9;
t16 = t398 * t48 + t402 * t59;
t3 = -qJD(6) * t16 - t398 * t9 + t402 * t6;
t427 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t493 = t166 * Ifges(5,4);
t560 = t167 / 0.2e1;
t565 = t146 / 0.2e1;
t566 = t145 / 0.2e1;
t618 = t427 + t19 * mrSges(6,1) + t254 * mrSges(5,1) + t10 / 0.2e1 + Ifges(6,3) * t560 - t20 * mrSges(6,2) + 0.2e1 * Ifges(6,5) * t565 + 0.2e1 * Ifges(6,6) * t566 - t493 / 0.2e1;
t601 = t394 * t621 + t620 * t396;
t600 = t620 * t394 - t396 * t621;
t362 = t394 * t402 + t396 * t398;
t352 = t362 * qJD(6);
t603 = t362 * t265;
t617 = t352 + t603;
t414 = t394 * t398 - t396 * t402;
t351 = t414 * qJD(6);
t602 = t414 * t265;
t616 = t602 + t351;
t615 = pkin(5) * t462 + pkin(11) * t619 + t601;
t614 = pkin(11) * t435 - t600;
t317 = t399 * t377 - t378 * t524;
t597 = -qJD(4) * t317 + t623 * t399 + t524 * t624;
t604 = t265 * t394;
t613 = pkin(5) * t604;
t612 = pkin(11) * t604;
t596 = pkin(4) * t445 - t597;
t223 = t399 * t226;
t123 = t207 * t524 - t223;
t115 = -t368 * pkin(4) + qJD(5) - t123;
t102 = -pkin(5) * t436 + t115;
t45 = t135 * t524 - t399 * t152 - t207 * t455 - t226 * t439;
t42 = -pkin(4) * t430 - t45;
t31 = -t145 * pkin(5) + t42;
t448 = Ifges(5,5) * t166 - Ifges(5,6) * t167 + Ifges(5,3) * t430;
t490 = t20 * t396;
t511 = Ifges(6,4) * t396;
t512 = Ifges(6,4) * t394;
t529 = t396 / 0.2e1;
t263 = qJD(6) + t265;
t549 = t263 / 0.2e1;
t550 = -t263 / 0.2e1;
t561 = t151 / 0.2e1;
t562 = -t151 / 0.2e1;
t563 = t595 / 0.2e1;
t564 = -t595 / 0.2e1;
t147 = Ifges(7,4) * t595;
t80 = Ifges(7,1) * t151 + Ifges(7,5) * t263 + t147;
t568 = t80 / 0.2e1;
t57 = t146 * Ifges(6,4) + t145 * Ifges(6,2) + t167 * Ifges(6,6);
t510 = Ifges(7,4) * t151;
t79 = Ifges(7,2) * t595 + Ifges(7,6) * t263 + t510;
t570 = t79 / 0.2e1;
t58 = t146 * Ifges(6,1) + t145 * Ifges(6,4) + t167 * Ifges(6,5);
t572 = t58 / 0.2e1;
t573 = t53 / 0.2e1;
t574 = t52 / 0.2e1;
t575 = Ifges(7,1) * t574 + Ifges(7,4) * t573 + Ifges(7,5) * t560;
t576 = Ifges(7,4) * t574 + Ifges(7,2) * t573 + Ifges(7,6) * t560;
t609 = (mrSges(7,1) * t617 - mrSges(7,2) * t616) * t102 + (t15 * t616 - t16 * t617 - t2 * t414 - t3 * t362) * mrSges(7,3) + (-Ifges(7,1) * t351 - Ifges(7,4) * t352) * t561 + (-Ifges(7,5) * t351 - Ifges(7,6) * t352) * t549 + (-Ifges(7,4) * t351 - Ifges(7,2) * t352) * t563 - t414 * t576 + t31 * (mrSges(7,1) * t414 + mrSges(7,2) * t362) + (Ifges(7,1) * t362 - Ifges(7,4) * t414) * t574 + (Ifges(7,4) * t362 - Ifges(7,2) * t414) * t573 + (Ifges(6,5) * t394 + Ifges(7,5) * t362 + Ifges(6,6) * t396 - Ifges(7,6) * t414) * t560 + t394 * t572 + t362 * t575 + t448 - t602 * t80 / 0.2e1 - t603 * t79 / 0.2e1 + (Ifges(7,1) * t602 + Ifges(7,4) * t603) * t562 + (Ifges(7,4) * t602 + Ifges(7,2) * t603) * t564 + (Ifges(7,5) * t602 + Ifges(7,6) * t603) * t550 + t42 * (-mrSges(6,1) * t396 + mrSges(6,2) * t394) + (Ifges(6,1) * t394 + t511) * t565 + (Ifges(6,2) * t396 + t512) * t566 - t351 * t568 - t352 * t570 + t57 * t529 + mrSges(6,3) * t490 - t44 * mrSges(5,2) + t45 * mrSges(5,1);
t113 = t235 * Ifges(6,4) + Ifges(6,2) * t436 + t265 * Ifges(6,6);
t608 = t113 / 0.2e1;
t546 = t265 / 0.2e1;
t392 = -pkin(3) * t403 - pkin(2);
t298 = -pkin(4) * t411 - qJ(5) * t364 + t392;
t239 = t396 * t298 - t317 * t394;
t393 = t396 * pkin(11);
t203 = -pkin(5) * t411 - t364 * t393 + t239;
t240 = t394 * t298 + t396 * t317;
t478 = t364 * t394;
t221 = -pkin(11) * t478 + t240;
t117 = t203 * t402 - t221 * t398;
t607 = qJD(6) * t117 + t615 * t398 - t614 * t402;
t118 = t203 * t398 + t221 * t402;
t606 = -qJD(6) * t118 + t614 * t398 + t615 * t402;
t605 = t167 * Ifges(5,2);
t502 = Ifges(6,6) * t436;
t508 = Ifges(6,5) * t235;
t112 = Ifges(6,3) * t265 + t502 + t508;
t499 = Ifges(7,3) * t263;
t500 = Ifges(7,6) * t595;
t506 = Ifges(7,5) * t151;
t78 = t499 + t500 + t506;
t591 = t112 + t78;
t598 = t435 * pkin(5) + t596;
t201 = pkin(4) * t406 + qJ(5) * t265;
t503 = Ifges(5,6) * t368;
t513 = Ifges(5,4) * t406;
t185 = -Ifges(5,2) * t265 + t503 + t513;
t578 = -t268 * mrSges(5,1) - t76 * mrSges(6,1) - t15 * mrSges(7,1) + t77 * mrSges(6,2) + t16 * mrSges(7,2) + t124 * mrSges(5,3) + t185 / 0.2e1;
t544 = -t406 / 0.2e1;
t592 = pkin(5) * t406;
t521 = pkin(3) * t399;
t388 = qJ(5) + t521;
t358 = (-pkin(11) - t388) * t394;
t477 = t388 * t396;
t359 = t393 + t477;
t291 = t358 * t402 - t359 * t398;
t432 = pkin(3) * t439;
t383 = t432 + qJD(5);
t481 = t265 * t396;
t128 = t225 * t524 - t223;
t522 = pkin(3) * t329;
t171 = t201 + t522;
t86 = -t128 * t394 + t396 * t171;
t63 = pkin(11) * t481 + t592 + t86;
t87 = t396 * t128 + t394 * t171;
t75 = t87 + t612;
t590 = qJD(6) * t291 - t383 * t414 - t398 * t63 - t402 * t75;
t292 = t358 * t398 + t359 * t402;
t589 = -qJD(6) * t292 - t362 * t383 + t398 * t75 - t402 * t63;
t371 = (-pkin(11) - qJ(5)) * t394;
t484 = qJ(5) * t396;
t372 = t393 + t484;
t307 = t371 * t402 - t372 * t398;
t92 = -t123 * t394 + t396 * t201;
t66 = t265 * t393 + t592 + t92;
t93 = t396 * t123 + t394 * t201;
t81 = t93 + t612;
t588 = -qJD(5) * t414 + qJD(6) * t307 - t398 * t66 - t402 * t81;
t308 = t371 * t398 + t372 * t402;
t587 = -qJD(5) * t362 - qJD(6) * t308 + t398 * t81 - t402 * t66;
t156 = -mrSges(6,1) * t436 + mrSges(6,2) * t235;
t245 = mrSges(5,1) * t368 - mrSges(5,3) * t406;
t584 = t245 - t156;
t489 = t265 * Ifges(5,6);
t581 = t406 * Ifges(5,5) + t368 * Ifges(5,3);
t583 = Ifges(4,5) * t329 - Ifges(4,6) * t408 + Ifges(4,3) * t376 - t489 + t581;
t340 = pkin(9) * t397 + t357;
t274 = t403 * t340 + t400 * t341;
t426 = mrSges(6,1) * t394 + mrSges(6,2) * t396;
t579 = -t268 * mrSges(5,2) + t123 * mrSges(5,3) - t115 * t426 + t394 * t608;
t577 = -0.2e1 * pkin(1);
t262 = Ifges(5,4) * t265;
t509 = Ifges(5,5) * t368;
t186 = Ifges(5,1) * t406 - t262 + t509;
t558 = t186 / 0.2e1;
t556 = -t436 / 0.2e1;
t555 = t436 / 0.2e1;
t554 = -t235 / 0.2e1;
t553 = t235 / 0.2e1;
t516 = Ifges(4,4) * t329;
t252 = -Ifges(4,2) * t408 + Ifges(4,6) * t376 + t516;
t552 = t252 / 0.2e1;
t325 = Ifges(4,4) * t408;
t253 = t329 * Ifges(4,1) + t376 * Ifges(4,5) - t325;
t551 = t253 / 0.2e1;
t547 = -t265 / 0.2e1;
t543 = t406 / 0.2e1;
t470 = t400 * t401;
t472 = t397 * t403;
t353 = -t395 * t470 + t472;
t354 = t397 * t400 + t403 * t475;
t412 = t353 * t524 - t399 * t354;
t542 = t412 / 0.2e1;
t284 = t399 * t353 + t354 * t524;
t540 = t284 / 0.2e1;
t539 = t294 / 0.2e1;
t538 = t295 / 0.2e1;
t537 = -t329 / 0.2e1;
t536 = t329 / 0.2e1;
t535 = t353 / 0.2e1;
t534 = t354 / 0.2e1;
t533 = -t368 / 0.2e1;
t532 = t368 / 0.2e1;
t531 = -t376 / 0.2e1;
t530 = t376 / 0.2e1;
t528 = t400 / 0.2e1;
t526 = t401 / 0.2e1;
t525 = t403 / 0.2e1;
t519 = t396 * pkin(5);
t458 = qJD(2) * t401;
t209 = -qJD(3) * t274 + t403 * t347 - t349 * t400;
t442 = t404 * t459;
t304 = qJD(3) * t353 + t403 * t442;
t443 = t395 * t458;
t170 = pkin(3) * t443 - pkin(10) * t304 + t209;
t208 = -t340 * t457 + t341 * t456 + t400 * t347 + t403 * t349;
t303 = -qJD(3) * t354 - t400 * t442;
t180 = pkin(10) * t303 + t208;
t273 = -t400 * t340 + t403 * t341;
t231 = -pkin(3) * t474 - t354 * pkin(10) + t273;
t246 = pkin(10) * t353 + t274;
t64 = t399 * t170 + t524 * t180 + t231 * t439 - t246 * t455;
t61 = (qJ(5) * t458 - qJD(5) * t404) * t395 + t64;
t199 = qJD(4) * t412 + t399 * t303 + t304 * t524;
t200 = qJD(4) * t284 - t303 * t524 + t399 * t304;
t271 = -t303 * pkin(3) + t350;
t85 = t200 * pkin(4) - t199 * qJ(5) - t284 * qJD(5) + t271;
t27 = t394 * t85 + t396 * t61;
t518 = mrSges(4,3) * t329;
t517 = Ifges(3,4) * t401;
t515 = Ifges(4,4) * t400;
t514 = Ifges(4,4) * t403;
t507 = Ifges(6,5) * t396;
t505 = Ifges(3,2) * t404;
t504 = Ifges(3,6) * t397;
t501 = Ifges(6,6) * t394;
t494 = t166 * Ifges(5,1);
t492 = t167 * Ifges(5,4);
t491 = t19 * t394;
t487 = t337 * mrSges(3,2);
t485 = t400 * Ifges(4,6);
t114 = t235 * Ifges(6,1) + Ifges(6,4) * t436 + Ifges(6,5) * t265;
t473 = t396 * t114;
t469 = t400 * t404;
t379 = Ifges(3,4) * t444;
t467 = t404 * (Ifges(3,1) * t445 + Ifges(3,5) * t433 + t379);
t158 = t399 * t231 + t524 * t246;
t140 = -qJ(5) * t474 + t158;
t339 = t384 + (-pkin(2) - t523) * t397;
t293 = -t353 * pkin(3) + t339;
t183 = -pkin(4) * t412 - t284 * qJ(5) + t293;
t99 = t396 * t140 + t394 * t183;
t187 = -t300 * t414 - t352 * t364;
t211 = t280 * t398 + t281 * t402;
t466 = t187 - t211;
t188 = -t300 * t362 + t351 * t364;
t210 = t280 * t402 - t281 * t398;
t465 = t188 - t210;
t464 = -mrSges(3,1) * t433 + mrSges(4,1) * t408 + t329 * mrSges(4,2) + mrSges(3,3) * t445;
t453 = t524 * pkin(3);
t449 = -Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t447 = Ifges(4,5) * t294 + Ifges(4,6) * t295 + Ifges(4,3) * t430;
t17 = -t53 * mrSges(7,1) + t52 * mrSges(7,2);
t26 = -t394 * t61 + t396 * t85;
t89 = -t145 * mrSges(6,1) + t146 * mrSges(6,2);
t98 = -t140 * t394 + t396 * t183;
t127 = t225 * t399 + t224;
t391 = -t453 - pkin(4);
t157 = t231 * t524 - t399 * t246;
t425 = Ifges(6,1) * t396 - t512;
t424 = Ifges(4,4) * t304 + Ifges(4,2) * t303;
t423 = -Ifges(4,2) * t400 + t514;
t422 = -Ifges(6,2) * t394 + t511;
t421 = Ifges(3,5) * t404 - Ifges(3,6) * t401;
t420 = -t501 + t507;
t419 = t490 - t491;
t418 = -t394 * t76 + t396 * t77;
t261 = t396 * t284 - t394 * t474;
t70 = -pkin(5) * t412 - pkin(11) * t261 + t98;
t260 = -t394 * t284 - t396 * t474;
t82 = pkin(11) * t260 + t99;
t28 = -t398 * t82 + t402 * t70;
t29 = t398 * t70 + t402 * t82;
t417 = t197 * t403 - t198 * t400;
t416 = -t255 * t403 - t256 * t400;
t193 = t260 * t402 - t261 * t398;
t194 = t260 * t398 + t261 * t402;
t415 = t345 * t404 + t348 * t401;
t142 = pkin(4) * t474 - t157;
t410 = t403 * t423;
t65 = t170 * t524 - t399 * t180 - t231 * t455 - t246 * t439;
t409 = Ifges(4,4) * t468 - Ifges(4,2) * t469 + Ifges(4,6) * t401;
t62 = -pkin(4) * t443 - t65;
t407 = t408 * mrSges(4,3);
t389 = -pkin(4) - t519;
t375 = Ifges(3,5) * t429;
t370 = t391 - t519;
t344 = -mrSges(3,2) * t433 + mrSges(3,3) * t444;
t309 = Ifges(3,6) * qJD(2) + (t504 + (t505 + t517) * t395) * qJD(1);
t297 = mrSges(4,1) * t376 - t518;
t296 = -t376 * mrSges(4,2) - t407;
t286 = t414 * t364;
t285 = t362 * t364;
t279 = pkin(5) * t478 - t582;
t270 = -mrSges(4,2) * t430 + mrSges(4,3) * t295;
t269 = mrSges(4,1) * t430 - mrSges(4,3) * t294;
t244 = -mrSges(5,2) * t368 - mrSges(5,3) * t265;
t227 = -mrSges(4,1) * t295 + mrSges(4,2) * t294;
t218 = t294 * Ifges(4,1) + t295 * Ifges(4,4) + Ifges(4,5) * t430;
t217 = t294 * Ifges(4,4) + t295 * Ifges(4,2) + Ifges(4,6) * t430;
t202 = mrSges(5,1) * t265 + mrSges(5,2) * t406;
t182 = t199 * t396 + t394 * t443;
t181 = -t199 * t394 + t396 * t443;
t177 = mrSges(6,1) * t265 - mrSges(6,3) * t235;
t176 = -mrSges(6,2) * t265 + mrSges(6,3) * t436;
t154 = -mrSges(5,2) * t430 - mrSges(5,3) * t167;
t153 = mrSges(5,1) * t430 - mrSges(5,3) * t166;
t110 = mrSges(7,1) * t263 - mrSges(7,3) * t151;
t109 = -mrSges(7,2) * t263 + mrSges(7,3) * t595;
t108 = -t260 * pkin(5) + t142;
t106 = t127 - t613;
t103 = t124 - t613;
t101 = mrSges(5,1) * t167 + mrSges(5,2) * t166;
t97 = Ifges(5,5) * t430 - t492 + t494;
t96 = Ifges(5,6) * t430 + t493 - t605;
t95 = mrSges(6,1) * t167 - mrSges(6,3) * t146;
t94 = -mrSges(6,2) * t167 + mrSges(6,3) * t145;
t91 = -mrSges(7,1) * t595 + mrSges(7,2) * t151;
t74 = -qJD(6) * t194 + t181 * t402 - t182 * t398;
t73 = qJD(6) * t193 + t181 * t398 + t182 * t402;
t38 = -t181 * pkin(5) + t62;
t35 = -mrSges(7,2) * t167 + mrSges(7,3) * t53;
t34 = mrSges(7,1) * t167 - mrSges(7,3) * t52;
t21 = pkin(11) * t181 + t27;
t13 = pkin(5) * t200 - pkin(11) * t182 + t26;
t5 = -qJD(6) * t29 + t13 * t402 - t21 * t398;
t4 = qJD(6) * t28 + t13 * t398 + t21 * t402;
t1 = [((-Ifges(6,3) - Ifges(7,3)) * t560 - Ifges(7,6) * t573 - Ifges(7,5) * t574 - t605 / 0.2e1 + t44 * mrSges(5,3) - t618) * t412 + t166 * (Ifges(5,1) * t284 - Ifges(5,5) * t474) / 0.2e1 + (Ifges(5,1) * t199 + Ifges(5,5) * t443) * t543 + (Ifges(5,5) * t199 + Ifges(5,3) * t443) * t532 + (Ifges(7,5) * t561 + Ifges(7,6) * t563 + Ifges(6,5) * t553 + Ifges(6,6) * t555 + Ifges(6,3) * t546 - Ifges(5,2) * t547 + Ifges(7,3) * t549 - Ifges(5,4) * t543 - Ifges(5,6) * t532 - t578 + t591 / 0.2e1) * t200 + (Ifges(7,1) * t194 + Ifges(7,4) * t193) * t574 + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t561 - t167 * (Ifges(5,4) * t284 - Ifges(5,6) * t474) / 0.2e1 + (Ifges(5,4) * t199 + Ifges(5,6) * t443) * t547 + (Ifges(6,5) * t261 + Ifges(7,5) * t194 + Ifges(6,6) * t260 + Ifges(7,6) * t193) * t560 + (Ifges(6,4) * t261 + Ifges(6,2) * t260) * t566 + (Ifges(6,4) * t182 + Ifges(6,2) * t181) * t555 - (t447 + t448) * t474 / 0.2e1 + (Ifges(7,4) * t194 + Ifges(7,2) * t193) * t573 + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t563 + m(7) * (t102 * t38 + t108 * t31 + t15 * t5 + t16 * t4 + t2 * t29 + t28 * t3) + m(6) * (t115 * t62 + t142 * t42 + t19 * t98 + t20 * t99 + t26 * t76 + t27 * t77) + m(5) * (t123 * t65 + t124 * t64 + t157 * t45 + t158 * t44 + t254 * t293 + t268 * t271) + m(4) * (t197 * t274 + t198 * t273 + t208 * t256 + t209 * t255 + t312 * t350 + t338 * t339) + m(3) * (t337 * t357 - t338 * t356 - t345 * t350 + t348 * t349) + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t549 + (-qJD(2) * t415 + t337 * t404 + t338 * t401) * t395 * mrSges(3,3) + (Ifges(6,5) * t182 + Ifges(6,6) * t181) * t546 + (-t15 * t73 + t16 * t74 + t193 * t2 - t194 * t3) * mrSges(7,3) + t74 * t570 + t261 * t572 + t194 * t575 + t193 * t576 + (Ifges(6,1) * t261 + Ifges(6,4) * t260) * t565 + (Ifges(6,1) * t182 + Ifges(6,4) * t181) * t553 + (t401 * t583 + t467) * t459 / 0.2e1 + (t181 * t77 - t182 * t76 - t19 * t261 + t20 * t260) * mrSges(6,3) + (-t124 * t443 + t199 * t268 + t254 * t284 + t44 * t474) * mrSges(5,2) + (t424 * t535 + ((-t356 * mrSges(3,3) + Ifges(3,5) * t397 + (mrSges(3,2) * t577 + 0.3e1 / 0.2e1 * Ifges(3,4) * t404) * t395) * t404 + (-0.3e1 / 0.2e1 * t504 + Ifges(4,5) * t534 + Ifges(5,5) * t540 + Ifges(5,6) * t542 - t357 * mrSges(3,3) + (t472 / 0.2e1 + t535) * Ifges(4,6) + (mrSges(3,1) * t577 + (-0.3e1 / 0.2e1 * Ifges(3,4) - t485 / 0.2e1) * t401) * t395 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 - Ifges(5,3) / 0.2e1) * t474) * t401) * t459) * qJD(1) + qJD(2) * (Ifges(4,6) * t443 + t424) * t525 + qJD(2) ^ 2 * t421 * t395 / 0.2e1 + (-t487 - t338 * mrSges(3,1) + t375 / 0.2e1) * t397 + t45 * (-mrSges(5,1) * t474 - t284 * mrSges(5,3)) + t198 * (-mrSges(4,1) * t474 - t354 * mrSges(4,3)) + t197 * (mrSges(4,2) * t474 + t353 * mrSges(4,3)) + t464 * t350 - t309 * t443 / 0.2e1 + t256 * (-mrSges(4,2) * t443 + mrSges(4,3) * t303) + t123 * (mrSges(5,1) * t443 - mrSges(5,3) * t199) + t255 * (mrSges(4,1) * t443 - mrSges(4,3) * t304) + t73 * t568 + t199 * t558 + t304 * t551 + t303 * t552 + t217 * t535 + (Ifges(4,1) * t304 + Ifges(4,4) * t303 + Ifges(4,5) * t443) * t536 + (Ifges(4,4) * t354 + Ifges(4,2) * t353 - Ifges(4,6) * t474) * t538 + (Ifges(4,1) * t354 + Ifges(4,4) * t353 - Ifges(4,5) * t474) * t539 + t97 * t540 + t96 * t542 + (Ifges(4,5) * t304 + Ifges(4,6) * t303 + Ifges(4,3) * t443) * t530 + t218 * t534 + t338 * (-mrSges(4,1) * t353 + mrSges(4,2) * t354) + t349 * t344 + t339 * t227 + t312 * (-mrSges(4,1) * t303 + mrSges(4,2) * t304) + t293 * t101 + t208 * t296 + t209 * t297 + t271 * t202 + t273 * t269 + t274 * t270 + t42 * (-mrSges(6,1) * t260 + mrSges(6,2) * t261) + t260 * t57 / 0.2e1 + t64 * t244 + t65 * t245 + t31 * (-mrSges(7,1) * t193 + mrSges(7,2) * t194) + t115 * (-mrSges(6,1) * t181 + mrSges(6,2) * t182) + t182 * t114 / 0.2e1 + t27 * t176 + t26 * t177 + t62 * t156 + t157 * t153 + t158 * t154 + t142 * t89 + t28 * t34 + t29 * t35 + t181 * t608 + t38 * t91 + t98 * t95 + t99 * t94 + t102 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t108 * t17 + t4 * t109 + t5 * t110; (mrSges(6,1) * t462 + mrSges(6,3) * t619) * t76 + (mrSges(6,1) * t435 - mrSges(6,2) * t619) * t115 - (t50 / 0.2e1 - t96 / 0.2e1 + t51 / 0.2e1 + (Ifges(7,3) / 0.2e1 - t449) * t167 + t618) * t411 + (-t102 * t465 + t15 * t462 + t285 * t31) * mrSges(7,1) + (-Ifges(7,4) * t286 - Ifges(7,2) * t285) * t573 + (-Ifges(7,1) * t286 - Ifges(7,4) * t285) * t574 + (-t15 * t466 + t16 * t465 - t2 * t285 + t286 * t3) * mrSges(7,3) + (-Ifges(7,5) * t286 - Ifges(7,6) * t285) * t560 + (t102 * t466 - t16 * t462 - t286 * t31) * mrSges(7,2) + (-t269 * t400 + t270 * t403) * pkin(9) + (-t123 * t463 - t124 * t462 - t364 * t45 + t411 * t44) * mrSges(5,3) + ((Ifges(4,1) * t468 - Ifges(4,4) * t469 + Ifges(4,5) * t401) * t537 - t255 * (mrSges(4,1) * t401 - mrSges(4,3) * t468) - t312 * (mrSges(4,1) * t469 + mrSges(4,2) * t468) - t256 * (-mrSges(4,2) * t401 - mrSges(4,3) * t469) + (Ifges(4,5) * t468 - Ifges(4,6) * t469 + Ifges(4,3) * t401) * t531 - t404 * t379 / 0.2e1 - t467 / 0.2e1 + t309 * t526 + (-t397 * t421 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t401 + mrSges(3,2) * t404) + t505 * t526 + t409 * t470 / 0.2e1) * t395 - t409 * t472 / 0.2e1) * qJD(1) - t123 * mrSges(5,1) * t401 + t124 * mrSges(5,2) * t401 + t489 * t526 - t253 * t468 / 0.2e1 + t469 * t552 + t415 * mrSges(3,3) + ((-Ifges(3,5) / 0.2e1 - t410 / 0.2e1) * t404 + (-Ifges(3,6) / 0.2e1 + Ifges(4,5) * t528 + Ifges(5,5) * t364 / 0.2e1 + Ifges(5,6) * t411 / 0.2e1) * t401) * qJD(2) - ((Ifges(3,1) * t404 - t517) * t461 + t423 * t457 + t581 + t583) * t401 / 0.2e1) * t461 + (t425 * t565 + t42 * t426 + t420 * t560 + t422 * t566 + t254 * mrSges(5,2) + t494 / 0.2e1 - t492 / 0.2e1 + t97 / 0.2e1 - t394 * t57 / 0.2e1 + t58 * t529 + (-t19 * t396 - t20 * t394) * mrSges(6,3)) * t364 - t286 * t575 - t285 * t576 + t375 + (-t315 / 0.2e1 + t300 / 0.2e1) * t186 + (t433 * t410 / 0.2e1 - t485 * t530 - t515 * t536 + t253 * t525 + t416 * mrSges(4,3) + (mrSges(4,2) * t312 + Ifges(4,1) * t536 + Ifges(4,5) * t530 - pkin(9) * t297) * t403 + (pkin(3) * t202 + t312 * mrSges(4,1) - t252 / 0.2e1 - pkin(9) * t296) * t400) * qJD(3) - t487 - (t89 - t153) * t582 + (t115 * t596 + t19 * t239 + t20 * t240 - t42 * t582 + t600 * t77 + t601 * t76) * m(6) + (t597 * t123 + t599 * t124 + t254 * t392 + t268 * t622 + t317 * t44 + t582 * t45) * m(5) + (-t281 / 0.2e1 + t479 / 0.2e1) * t114 + (-t280 / 0.2e1 - t480 / 0.2e1) * t113 + (mrSges(5,1) * t462 + mrSges(5,2) * t463) * t268 - t464 * t348 + (-mrSges(6,2) * t462 - mrSges(6,3) * t435) * t77 + (-mrSges(4,1) * t403 + mrSges(4,2) * t400 - mrSges(3,1)) * t338 + t392 * t101 + (Ifges(7,1) * t187 + Ifges(7,4) * t188 + Ifges(7,5) * t301) * t561 + (Ifges(7,1) * t211 + Ifges(7,4) * t210 + Ifges(7,5) * t314) * t562 + (Ifges(7,4) * t187 + Ifges(7,2) * t188 + Ifges(7,6) * t301) * t563 + (Ifges(7,4) * t211 + Ifges(7,2) * t210 + Ifges(7,6) * t314) * t564 + (Ifges(6,5) * t301 + t300 * t425) * t553 + (Ifges(6,1) * t281 + Ifges(6,4) * t280 + Ifges(6,5) * t314) * t554 + (Ifges(6,6) * t301 + t300 * t422) * t555 + (Ifges(6,4) * t281 + Ifges(6,2) * t280 + Ifges(6,6) * t314) * t556 + (Ifges(5,1) * t300 - Ifges(5,4) * t301) * t543 + (Ifges(5,1) * t315 - Ifges(5,4) * t314) * t544 + (Ifges(7,5) * t187 + Ifges(7,6) * t188 + Ifges(7,3) * t301) * t549 + (Ifges(7,5) * t211 + Ifges(7,6) * t210 + Ifges(7,3) * t314) * t550 + (Ifges(4,2) * t403 + t515) * t538 + (Ifges(4,1) * t400 + t514) * t539 + t217 * t525 + t218 * t528 + (Ifges(5,5) * t300 - Ifges(5,6) * t301) * t532 + (Ifges(5,5) * t315 - Ifges(5,6) * t314) * t533 - t345 * t344 + (-t211 / 0.2e1 + t187 / 0.2e1) * t80 + t317 * t154 - t305 * t202 - t278 * t296 - t277 * t297 + t279 * t17 + t239 * t95 + t240 * t94 - pkin(2) * t227 + (-t210 / 0.2e1 + t188 / 0.2e1) * t79 + t417 * mrSges(4,3) + t596 * t156 + t597 * t245 + t598 * t91 + t599 * t244 + t600 * t176 + t601 * t177 + (t185 - t591) * (t314 / 0.2e1 - t301 / 0.2e1) + t606 * t110 + t607 * t109 + (t102 * t598 + t117 * t3 + t118 * t2 + t15 * t606 + t16 * t607 + t279 * t31) * m(7) + (Ifges(5,4) * t315 - Ifges(5,2) * t314 + Ifges(6,3) * t301 + t300 * t420) * t546 + (Ifges(5,4) * t300 + Ifges(6,5) * t281 - Ifges(5,2) * t301 + Ifges(6,6) * t280 + Ifges(6,3) * t314) * t547 + ((qJD(3) * t416 + t417) * pkin(9) - pkin(2) * t338 - t255 * t277 - t256 * t278 - t312 * t348) * m(4) + t117 * t34 + t118 * t35; (m(6) * t115 + m(7) * t102 - t584 + t91) * pkin(3) * t455 + t609 + (-Ifges(5,4) * t544 + Ifges(6,5) * t554 + Ifges(7,5) * t562 - Ifges(5,2) * t546 - Ifges(5,6) * t533 + Ifges(6,6) * t556 + Ifges(7,6) * t564 + Ifges(6,3) * t547 + Ifges(7,3) * t550 + t578) * t406 + (-t296 - t407) * t255 + t408 * (-Ifges(4,2) * t329 - t325) / 0.2e1 - (Ifges(5,1) * t544 + Ifges(5,4) * t546 + Ifges(5,5) * t533 + t420 * t547 + t422 * t556 + t425 * t554 + t579) * t265 + (t518 + t297) * t256 - t312 * (t329 * mrSges(4,1) - mrSges(4,2) * t408) + (t123 * t127 - t124 * t128 - t268 * t522 + (t524 * t45 + t399 * t44 + (-t123 * t399 + t124 * t524) * qJD(4)) * pkin(3)) * m(5) + t447 + t589 * t110 + t590 * t109 + (-t102 * t106 + t15 * t589 + t16 * t590 + t2 * t292 + t291 * t3 + t31 * t370) * m(7) + (t432 - t128) * t244 + (-t481 * t76 - t604 * t77 - t491) * mrSges(6,3) + (-t383 * t394 - t86) * t177 + t584 * t127 + (t383 * t396 - t87) * t176 - t202 * t522 + t391 * t89 + t370 * t17 + t408 * t551 + t252 * t536 + (-Ifges(4,1) * t408 - t516) * t537 + (-Ifges(4,5) * t408 - Ifges(4,6) * t329) * t531 + t154 * t521 + t591 * t544 + t94 * t477 + t291 * t34 + t292 * t35 + (-t115 * t127 + t383 * t418 + t388 * t419 + t391 * t42 - t76 * t86 - t77 * t87) * m(6) - t197 * mrSges(4,2) + t198 * mrSges(4,1) - t394 * t388 * t95 + (t473 + t186) * t546 + t153 * t453 - t106 * t91; t609 + (-t508 / 0.2e1 - t502 / 0.2e1 + t503 / 0.2e1 - t499 / 0.2e1 - t506 / 0.2e1 - t500 / 0.2e1 - t78 / 0.2e1 - t112 / 0.2e1 + t513 / 0.2e1 + t578) * t406 + (qJD(5) * t396 - t93) * t176 + (-t19 * mrSges(6,3) - qJ(5) * t95 - qJD(5) * t177) * t394 + (t425 * t553 + t422 * t555 + t509 / 0.2e1 - t262 / 0.2e1 + t558 + t473 / 0.2e1 + (t507 / 0.2e1 - t501 / 0.2e1) * t265 + (-t394 * t77 - t396 * t76) * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + t449) * t406 - t579) * t265 + t587 * t110 + t588 * t109 + (-t102 * t103 + t15 * t587 + t16 * t588 + t2 * t308 + t3 * t307 + t31 * t389) * m(7) + t584 * t124 + (-pkin(4) * t42 + qJ(5) * t419 + qJD(5) * t418 - t115 * t124 - t76 * t92 - t77 * t93) * m(6) + t389 * t17 + t94 * t484 + t307 * t34 + t308 * t35 - t123 * t244 - t92 * t177 - pkin(4) * t89 - t103 * t91; -t595 * t109 + t151 * t110 - t436 * t176 + t235 * t177 + t17 + t89 + (t15 * t151 - t16 * t595 + t31) * m(7) + (t235 * t76 - t436 * t77 + t42) * m(6); -t102 * (mrSges(7,1) * t151 + mrSges(7,2) * t595) + (Ifges(7,1) * t595 - t510) * t562 + t79 * t561 + (Ifges(7,5) * t595 - Ifges(7,6) * t151) * t550 - t15 * t109 + t16 * t110 + (t15 * t595 + t151 * t16) * mrSges(7,3) + t427 + t10 + (-Ifges(7,2) * t151 + t147 + t80) * t564;];
tauc  = t1(:);
