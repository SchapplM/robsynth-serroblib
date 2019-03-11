% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:18
% EndTime: 2019-03-09 18:07:03
% DurationCPUTime: 25.27s
% Computational Cost: add. (27282->879), mult. (63092->1146), div. (0->0), fcn. (48017->18), ass. (0->433)
t683 = -mrSges(6,3) - mrSges(7,3);
t414 = sin(qJ(5));
t498 = qJD(5) * t414;
t415 = sin(qJ(3));
t416 = sin(qJ(2));
t420 = cos(qJ(3));
t421 = cos(qJ(2));
t346 = -t415 * t416 + t420 * t421;
t324 = t346 * qJD(1);
t348 = t415 * t421 + t416 * t420;
t325 = t348 * qJD(1);
t411 = sin(pkin(11));
t412 = cos(pkin(11));
t642 = t412 * t324 - t411 * t325;
t536 = t642 * t414;
t687 = t498 - t536;
t408 = qJD(2) + qJD(3);
t588 = t408 / 0.2e1;
t447 = t324 * t411 + t412 * t325;
t654 = t447 * Ifges(5,4);
t684 = t642 * Ifges(5,2) / 0.2e1;
t686 = -Ifges(5,6) * t588 - t654 / 0.2e1 - t684;
t419 = cos(qJ(5));
t239 = t408 * t419 - t414 * t447;
t252 = qJD(5) - t642;
t240 = t408 * t414 + t419 * t447;
t548 = t240 * Ifges(6,4);
t129 = t239 * Ifges(6,2) + t252 * Ifges(6,6) + t548;
t685 = -t129 / 0.2e1;
t521 = t411 * t415;
t553 = pkin(2) * qJD(3);
t313 = (t412 * t420 - t521) * t553;
t583 = pkin(2) * t420;
t391 = pkin(3) + t583;
t519 = t412 * t415;
t315 = pkin(2) * t519 + t411 * t391;
t310 = pkin(9) + t315;
t568 = -pkin(10) - t310;
t472 = qJD(5) * t568;
t493 = pkin(10) * t536;
t424 = -pkin(8) - pkin(7);
t376 = t424 * t421;
t356 = qJD(1) * t376;
t329 = t420 * t356;
t375 = t424 * t416;
t355 = qJD(1) * t375;
t273 = -t355 * t415 + t329;
t540 = qJ(4) * t324;
t241 = t273 - t540;
t326 = t415 * t356;
t274 = t420 * t355 + t326;
t316 = t325 * qJ(4);
t242 = -t316 + t274;
t173 = t241 * t411 + t242 * t412;
t582 = pkin(3) * t325;
t183 = pkin(4) * t447 - pkin(9) * t642 + t582;
t503 = qJD(1) * t416;
t395 = pkin(2) * t503;
t176 = t183 + t395;
t92 = t419 * t173 + t414 * t176;
t682 = t313 * t419 + t414 * t472 + t493 - t92;
t535 = t642 * t419;
t463 = pkin(5) * t447 - pkin(10) * t535;
t91 = -t173 * t414 + t419 * t176;
t681 = -t313 * t414 + t419 * t472 - t463 - t91;
t580 = pkin(3) * t411;
t387 = pkin(9) + t580;
t567 = -pkin(10) - t387;
t471 = qJD(5) * t567;
t335 = qJD(2) * pkin(2) + t355;
t267 = t335 * t415 - t329;
t232 = t267 + t540;
t222 = t411 * t232;
t266 = t420 * t335 + t326;
t231 = t266 - t316;
t151 = t231 * t412 - t222;
t90 = t419 * t151 + t414 * t183;
t680 = t414 * t471 + t493 - t90;
t89 = -t151 * t414 + t419 * t183;
t679 = t419 * t471 - t463 - t89;
t410 = qJ(2) + qJ(3);
t399 = pkin(11) + t410;
t385 = sin(t399);
t543 = t414 * mrSges(6,2);
t409 = qJ(5) + qJ(6);
t400 = sin(t409);
t563 = mrSges(7,2) * t400;
t678 = (-t543 - t563) * t385;
t219 = pkin(3) * t408 + t231;
t143 = t219 * t412 - t222;
t405 = t421 * pkin(2);
t392 = t405 + pkin(1);
t374 = t392 * qJD(1);
t281 = -pkin(3) * t324 + qJD(4) - t374;
t658 = t642 * Ifges(5,4);
t677 = t281 * mrSges(5,2) - t143 * mrSges(5,3) + t414 * t685 + t658 / 0.2e1;
t620 = m(7) * pkin(5);
t413 = sin(qJ(6));
t418 = cos(qJ(6));
t165 = t239 * t413 + t240 * t418;
t249 = qJD(6) + t252;
t465 = t418 * t239 - t240 * t413;
t659 = t252 * Ifges(6,3);
t660 = t239 * Ifges(6,6);
t675 = t240 * Ifges(6,5) + t165 * Ifges(7,5) + Ifges(7,6) * t465 + t249 * Ifges(7,3) + t659 + t660;
t444 = t413 * t414 - t418 * t419;
t191 = t444 * t642;
t630 = qJD(5) + qJD(6);
t275 = t630 * t444;
t645 = -t275 + t191;
t347 = t413 * t419 + t414 * t418;
t190 = t347 * t642;
t276 = t630 * t347;
t674 = t276 - t190;
t641 = t412 * t241 - t242 * t411 + (t411 * t420 + t519) * t553;
t386 = cos(t399);
t401 = sin(t410);
t403 = cos(t410);
t673 = mrSges(4,1) * t401 + mrSges(5,1) * t385 + mrSges(4,2) * t403 + mrSges(5,2) * t386;
t495 = qJD(1) * qJD(2);
t359 = qJDD(1) * t421 - t416 * t495;
t402 = cos(t409);
t565 = mrSges(7,1) * t402;
t566 = mrSges(6,1) * t419;
t672 = (t565 + t566) * t385;
t671 = t687 * pkin(5);
t670 = t543 - t566;
t520 = t412 * t232;
t144 = t411 * t219 + t520;
t141 = pkin(9) * t408 + t144;
t154 = -pkin(4) * t642 - pkin(9) * t447 + t281;
t406 = qJDD(2) + qJDD(3);
t360 = qJDD(1) * t416 + t421 * t495;
t342 = t360 * pkin(7);
t280 = qJDD(2) * pkin(2) - pkin(8) * t360 - t342;
t341 = t359 * pkin(7);
t282 = pkin(8) * t359 + t341;
t182 = -qJD(3) * t267 + t420 * t280 - t282 * t415;
t435 = t346 * qJD(3);
t243 = qJD(1) * t435 + t359 * t415 + t360 * t420;
t107 = pkin(3) * t406 - qJ(4) * t243 - qJD(4) * t325 + t182;
t499 = qJD(3) * t420;
t500 = qJD(3) * t415;
t181 = t415 * t280 + t420 * t282 + t335 * t499 + t356 * t500;
t436 = t348 * qJD(3);
t244 = -qJD(1) * t436 + t359 * t420 - t360 * t415;
t118 = qJ(4) * t244 + qJD(4) * t324 + t181;
t54 = t411 * t107 + t412 * t118;
t49 = pkin(9) * t406 + t54;
t497 = qJD(5) * t419;
t174 = -t411 * t243 + t244 * t412;
t175 = t243 * t412 + t244 * t411;
t539 = qJDD(1) * pkin(1);
t319 = -pkin(2) * t359 - t539;
t207 = -pkin(3) * t244 + qJDD(4) + t319;
t71 = -pkin(4) * t174 - pkin(9) * t175 + t207;
t15 = -t141 * t498 + t154 * t497 + t414 * t71 + t419 * t49;
t83 = t141 * t419 + t154 * t414;
t16 = -qJD(5) * t83 - t414 * t49 + t419 * t71;
t669 = t15 * t419 - t16 * t414;
t417 = sin(qJ(1));
t422 = cos(qJ(1));
t668 = g(1) * t422 + g(2) * t417;
t667 = -t403 * mrSges(4,1) - t386 * mrSges(5,1) + mrSges(4,2) * t401 + (mrSges(5,2) + t683) * t385;
t68 = pkin(10) * t239 + t83;
t544 = t413 * t68;
t82 = -t141 * t414 + t419 * t154;
t67 = -pkin(10) * t240 + t82;
t59 = pkin(5) * t252 + t67;
t20 = t418 * t59 - t544;
t541 = t418 * t68;
t21 = t413 * t59 + t541;
t665 = -t281 * mrSges(5,1) - t82 * mrSges(6,1) - t20 * mrSges(7,1) + t83 * mrSges(6,2) + t21 * mrSges(7,2) + mrSges(5,3) * t144 - t686;
t119 = qJD(5) * t239 + t175 * t419 + t406 * t414;
t120 = -qJD(5) * t240 - t175 * t414 + t406 * t419;
t42 = qJD(6) * t465 + t119 * t418 + t120 * t413;
t617 = t42 / 0.2e1;
t43 = -qJD(6) * t165 - t119 * t413 + t120 * t418;
t616 = t43 / 0.2e1;
t610 = t119 / 0.2e1;
t609 = t120 / 0.2e1;
t169 = qJDD(5) - t174;
t159 = qJDD(6) + t169;
t608 = t159 / 0.2e1;
t603 = t169 / 0.2e1;
t664 = t359 / 0.2e1;
t590 = t406 / 0.2e1;
t587 = t421 / 0.2e1;
t663 = t414 * t620;
t662 = Ifges(4,5) * t348;
t661 = Ifges(4,6) * t346;
t656 = t421 * Ifges(3,2);
t655 = t447 * Ifges(5,1);
t653 = mrSges(6,1) + t620;
t283 = t568 * t414;
t404 = t419 * pkin(10);
t532 = t310 * t419;
t284 = t404 + t532;
t212 = t283 * t418 - t284 * t413;
t652 = qJD(6) * t212 + t413 * t681 + t418 * t682;
t213 = t283 * t413 + t284 * t418;
t651 = -qJD(6) * t213 - t413 * t682 + t418 * t681;
t336 = t567 * t414;
t526 = t387 * t419;
t337 = t404 + t526;
t262 = t336 * t418 - t337 * t413;
t650 = qJD(6) * t262 + t413 * t679 + t418 * t680;
t263 = t336 * t413 + t337 * t418;
t649 = -qJD(6) * t263 - t413 * t680 + t418 * t679;
t140 = -pkin(4) * t408 - t143;
t456 = mrSges(6,1) * t414 + mrSges(6,2) * t419;
t648 = t140 * t456;
t270 = t346 * t411 + t348 * t412;
t206 = t444 * t270;
t285 = t420 * t375 + t376 * t415;
t255 = -qJ(4) * t348 + t285;
t286 = t415 * t375 - t420 * t376;
t256 = qJ(4) * t346 + t286;
t194 = t255 * t411 + t256 * t412;
t189 = t419 * t194;
t269 = -t412 * t346 + t348 * t411;
t300 = -pkin(3) * t346 - t392;
t192 = pkin(4) * t269 - pkin(9) * t270 + t300;
t97 = t414 * t192 + t189;
t245 = -mrSges(5,2) * t408 + mrSges(5,3) * t642;
t560 = mrSges(6,3) * t239;
t185 = -mrSges(6,2) * t252 + t560;
t559 = mrSges(6,3) * t240;
t186 = mrSges(6,1) * t252 - t559;
t448 = t419 * t185 - t414 * t186;
t647 = -t245 - t448;
t646 = mrSges(5,1) * t408 + mrSges(6,1) * t239 - mrSges(6,2) * t240 - mrSges(5,3) * t447;
t643 = t671 + t641;
t575 = pkin(5) * t419;
t390 = pkin(4) + t575;
t423 = -pkin(10) - pkin(9);
t446 = -t385 * t423 + t386 * t390;
t462 = t386 * pkin(4) + t385 * pkin(9);
t445 = -t385 * t390 - t386 * t423;
t578 = pkin(4) * t385;
t581 = pkin(3) * t401;
t639 = -m(7) * (t445 - t581) - m(6) * (-t578 - t581) + t672;
t638 = -m(7) * t445 + t672;
t150 = t231 * t411 + t520;
t637 = -t150 + t671;
t502 = qJD(1) * t421;
t573 = pkin(7) * t421;
t574 = pkin(7) * t416;
t636 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t503) * t573 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t502) * t574;
t527 = t386 * t422;
t635 = t422 * t678 + t527 * t683;
t528 = t386 * t417;
t634 = t417 * t678 + t528 * t683;
t633 = t341 * t421 + t342 * t416;
t631 = m(5) + m(6) + m(7);
t629 = 0.2e1 * t590;
t627 = -m(3) * pkin(7) + m(4) * t424 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t626 = t667 + (t563 - t565 + t670) * t386;
t7 = pkin(5) * t169 - pkin(10) * t119 + t16;
t8 = pkin(10) * t120 + t15;
t3 = qJD(6) * t20 + t413 * t7 + t418 * t8;
t4 = -qJD(6) * t21 - t413 * t8 + t418 * t7;
t625 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t624 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t72 = mrSges(6,1) * t169 - mrSges(6,3) * t119;
t623 = m(6) * ((-t414 * t83 - t419 * t82) * qJD(5) + t669) - t186 * t497 - t185 * t498 - t414 * t72;
t372 = -mrSges(3,1) * t421 + mrSges(3,2) * t416;
t622 = m(3) * pkin(1) + m(4) * t392 + mrSges(2,1) - t372 - t667;
t619 = Ifges(7,4) * t617 + Ifges(7,2) * t616 + Ifges(7,6) * t608;
t618 = Ifges(7,1) * t617 + Ifges(7,4) * t616 + Ifges(7,5) * t608;
t615 = Ifges(6,1) * t610 + Ifges(6,4) * t609 + Ifges(6,5) * t603;
t554 = Ifges(7,4) * t165;
t76 = Ifges(7,2) * t465 + Ifges(7,6) * t249 + t554;
t614 = -t76 / 0.2e1;
t613 = t76 / 0.2e1;
t155 = Ifges(7,4) * t465;
t77 = t165 * Ifges(7,1) + t249 * Ifges(7,5) + t155;
t612 = -t77 / 0.2e1;
t611 = t77 / 0.2e1;
t607 = -t465 / 0.2e1;
t606 = t465 / 0.2e1;
t605 = -t165 / 0.2e1;
t604 = t165 / 0.2e1;
t601 = -t239 / 0.2e1;
t600 = -t240 / 0.2e1;
t599 = t240 / 0.2e1;
t598 = -t249 / 0.2e1;
t597 = t249 / 0.2e1;
t596 = -t252 / 0.2e1;
t591 = t325 / 0.2e1;
t589 = -t408 / 0.2e1;
t586 = mrSges(7,3) * t20;
t585 = mrSges(7,3) * t21;
t584 = pkin(2) * t416;
t389 = pkin(3) * t403;
t579 = pkin(3) * t412;
t577 = pkin(5) * t240;
t570 = g(3) * t385;
t562 = mrSges(4,3) * t324;
t558 = Ifges(3,4) * t416;
t557 = Ifges(3,4) * t421;
t556 = Ifges(6,4) * t414;
t555 = Ifges(6,4) * t419;
t552 = pkin(5) * qJD(6);
t547 = t325 * mrSges(4,3);
t546 = t325 * Ifges(4,4);
t277 = qJD(2) * t346 + t435;
t278 = -qJD(2) * t348 - t436;
t209 = t277 * t412 + t278 * t411;
t538 = t209 * t419;
t534 = t270 * t414;
t533 = t270 * t419;
t525 = t400 * t417;
t524 = t400 * t422;
t523 = t402 * t417;
t522 = t402 * t422;
t516 = t414 * t417;
t515 = t414 * t422;
t514 = t417 * t419;
t233 = Ifges(6,4) * t239;
t130 = t240 * Ifges(6,1) + t252 * Ifges(6,5) + t233;
t513 = t419 * t130;
t511 = t419 * t422;
t295 = t386 * t525 + t522;
t296 = -t386 * t523 + t524;
t510 = -t295 * mrSges(7,1) + t296 * mrSges(7,2);
t297 = -t386 * t524 + t523;
t298 = t386 * t522 + t525;
t509 = t297 * mrSges(7,1) - t298 * mrSges(7,2);
t504 = t389 + t405;
t501 = qJD(2) * t416;
t490 = Ifges(7,5) * t42 + Ifges(7,6) * t43 + Ifges(7,3) * t159;
t396 = pkin(2) * t501;
t483 = Ifges(6,5) * t119 + Ifges(6,6) * t120 + Ifges(6,3) * t169;
t482 = t389 + t462;
t481 = qJD(2) * t424;
t478 = t270 * t497;
t477 = t513 / 0.2e1;
t473 = -t498 / 0.2e1;
t261 = -pkin(3) * t278 + t396;
t470 = t495 / 0.2e1;
t469 = -t174 * mrSges(5,1) + t175 * mrSges(5,2);
t208 = t277 * t411 - t412 * t278;
t105 = pkin(4) * t208 - pkin(9) * t209 + t261;
t357 = t416 * t481;
t358 = t421 * t481;
t215 = t420 * t357 + t415 * t358 + t375 * t499 + t376 * t500;
t170 = qJ(4) * t278 + qJD(4) * t346 + t215;
t216 = -qJD(3) * t286 - t357 * t415 + t420 * t358;
t171 = -qJ(4) * t277 - qJD(4) * t348 + t216;
t87 = t170 * t412 + t171 * t411;
t466 = t419 * t105 - t414 * t87;
t53 = t107 * t412 - t411 * t118;
t86 = t170 * t411 - t412 * t171;
t96 = t419 * t192 - t194 * t414;
t193 = -t412 * t255 + t256 * t411;
t314 = -pkin(2) * t521 + t391 * t412;
t309 = -pkin(4) - t314;
t461 = t389 + t446;
t460 = mrSges(3,1) * t416 + mrSges(3,2) * t421;
t455 = -mrSges(7,1) * t400 - mrSges(7,2) * t402;
t454 = Ifges(6,1) * t419 - t556;
t453 = t558 + t656;
t452 = -Ifges(6,2) * t414 + t555;
t451 = Ifges(3,5) * t421 - Ifges(3,6) * t416;
t450 = Ifges(6,5) * t419 - Ifges(6,6) * t414;
t80 = pkin(5) * t269 - pkin(10) * t533 + t96;
t84 = -pkin(10) * t534 + t97;
t34 = -t413 * t84 + t418 * t80;
t35 = t413 * t80 + t418 * t84;
t449 = -t414 * t82 + t419 * t83;
t441 = t490 + t625;
t48 = -pkin(4) * t406 - t53;
t440 = pkin(1) * t460;
t307 = -t386 * t515 + t514;
t305 = t386 * t516 + t511;
t439 = t209 * t414 + t478;
t438 = t270 * t498 - t538;
t437 = t416 * (Ifges(3,1) * t421 - t558);
t22 = t414 * t105 + t192 * t497 - t194 * t498 + t419 * t87;
t110 = -pkin(5) * t239 + t140;
t197 = t408 * Ifges(5,5) + t655 + t658;
t253 = t324 * Ifges(4,2) + t408 * Ifges(4,6) + t546;
t320 = Ifges(4,4) * t324;
t254 = t325 * Ifges(4,1) + t408 * Ifges(4,5) + t320;
t33 = -pkin(5) * t120 + t48;
t44 = t119 * Ifges(6,4) + t120 * Ifges(6,2) + t169 * Ifges(6,6);
t425 = (Ifges(5,3) + Ifges(4,3)) * t406 + (-t687 * t83 + (-t497 + t535) * t82 + t669) * mrSges(6,3) + (-Ifges(7,1) * t191 - Ifges(7,4) * t190) * t605 + t266 * t562 - (t513 + t197) * t642 / 0.2e1 + t33 * (mrSges(7,1) * t444 + mrSges(7,2) * t347) + (Ifges(7,5) * t347 - Ifges(7,6) * t444) * t608 + (Ifges(7,4) * t347 - Ifges(7,2) * t444) * t616 + (Ifges(7,1) * t347 - Ifges(7,4) * t444) * t617 + (t190 * t21 - t191 * t20 - t3 * t444 - t347 * t4) * mrSges(7,3) - t444 * t619 - (Ifges(7,4) * t604 + Ifges(7,2) * t606 + Ifges(7,6) * t597 + t585 + t613) * t276 + t48 * t670 - (Ifges(7,1) * t604 + Ifges(7,4) * t606 + Ifges(7,5) * t597 - t586 + t611) * t275 + (Ifges(4,5) * t324 - Ifges(4,6) * t325) * t589 + t253 * t591 + (t239 * t452 + t240 * t454 + t252 * t450) * qJD(5) / 0.2e1 - (-Ifges(4,2) * t325 + t254 + t320) * t324 / 0.2e1 + (-Ifges(7,4) * t191 - Ifges(7,2) * t190) * t607 + (mrSges(7,1) * t674 + mrSges(7,2) * t645) * t110 + t374 * (mrSges(4,1) * t325 + mrSges(4,2) * t324) + (-Ifges(7,5) * t191 - Ifges(7,6) * t190) * t598 + (Ifges(6,5) * t600 + Ifges(7,5) * t605 - Ifges(5,6) * t589 + Ifges(6,6) * t601 + Ifges(7,6) * t607 + Ifges(6,3) * t596 + Ifges(7,3) * t598 + t665 + t684) * t447 + t129 * t473 - (Ifges(5,1) * t642 - t654 + t675) * t447 / 0.2e1 + (t648 + t477) * qJD(5) + (Ifges(5,5) * t589 + t450 * t596 + t452 * t601 + t454 * t600 - t648 - t677) * t642 - t325 * (Ifges(4,1) * t324 - t546) / 0.2e1 + t419 * t44 / 0.2e1 + Ifges(4,5) * t243 + Ifges(4,6) * t244 - t181 * mrSges(4,2) + t182 * mrSges(4,1) + Ifges(5,6) * t174 + Ifges(5,5) * t175 + (Ifges(6,5) * t414 + Ifges(6,6) * t419) * t603 + (Ifges(6,2) * t419 + t556) * t609 + (Ifges(6,1) * t414 + t555) * t610 - t191 * t612 - t190 * t614 + t414 * t615 + t347 * t618 + t53 * mrSges(5,1) - t54 * mrSges(5,2) + t267 * t547;
t407 = -qJ(4) + t424;
t394 = Ifges(3,4) * t502;
t388 = -pkin(4) - t579;
t369 = pkin(9) * t527;
t368 = pkin(9) * t528;
t367 = -t390 - t579;
t361 = -t581 - t584;
t354 = pkin(1) + t504;
t339 = t422 * t361;
t338 = t417 * t361;
t323 = Ifges(3,1) * t503 + Ifges(3,5) * qJD(2) + t394;
t322 = Ifges(3,6) * qJD(2) + qJD(1) * t453;
t308 = t386 * t511 + t516;
t306 = -t386 * t514 + t515;
t299 = t309 - t575;
t294 = mrSges(4,1) * t408 - t547;
t293 = -mrSges(4,2) * t408 + t562;
t292 = t395 + t582;
t265 = -mrSges(4,1) * t324 + mrSges(4,2) * t325;
t221 = -mrSges(4,2) * t406 + mrSges(4,3) * t244;
t220 = mrSges(4,1) * t406 - mrSges(4,3) * t243;
t205 = t347 * t270;
t202 = -mrSges(5,1) * t642 + mrSges(5,2) * t447;
t148 = mrSges(5,1) * t406 - mrSges(5,3) * t175;
t147 = -mrSges(5,2) * t406 + mrSges(5,3) * t174;
t137 = pkin(5) * t534 + t193;
t124 = mrSges(7,1) * t249 - mrSges(7,3) * t165;
t123 = -mrSges(7,2) * t249 + mrSges(7,3) * t465;
t88 = -mrSges(7,1) * t465 + mrSges(7,2) * t165;
t73 = -mrSges(6,2) * t169 + mrSges(6,3) * t120;
t61 = t206 * t630 - t347 * t209;
t60 = -t209 * t444 - t270 * t276;
t58 = pkin(5) * t439 + t86;
t55 = -mrSges(6,1) * t120 + mrSges(6,2) * t119;
t31 = -mrSges(7,2) * t159 + mrSges(7,3) * t43;
t30 = mrSges(7,1) * t159 - mrSges(7,3) * t42;
t25 = t418 * t67 - t544;
t24 = -t413 * t67 - t541;
t23 = -qJD(5) * t97 + t466;
t18 = -pkin(10) * t439 + t22;
t17 = -pkin(10) * t538 + pkin(5) * t208 + (-t189 + (pkin(10) * t270 - t192) * t414) * qJD(5) + t466;
t13 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t6 = -qJD(6) * t35 + t17 * t418 - t18 * t413;
t5 = qJD(6) * t34 + t17 * t413 + t18 * t418;
t1 = [(-t306 * mrSges(6,1) - t296 * mrSges(7,1) - t305 * mrSges(6,2) - t295 * mrSges(7,2) + (t407 * t631 + t627 - t663) * t422 + (-m(7) * (-t354 - t446) - m(6) * (-t354 - t462) + m(5) * t354 + t622) * t417) * g(1) + (-mrSges(3,1) * t574 - mrSges(3,2) * t573 + 0.2e1 * Ifges(3,6) * t587) * qJDD(2) + (Ifges(3,1) * t360 + Ifges(3,4) * t664 + Ifges(3,5) * qJDD(2) - t470 * t656) * t416 + (-t516 * t620 - t308 * mrSges(6,1) - t298 * mrSges(7,1) - t307 * mrSges(6,2) - t297 * mrSges(7,2) - t631 * (t422 * t354 - t407 * t417) + t627 * t417 + (-m(6) * t462 - m(7) * t446 - t622) * t422) * g(2) + t478 * t685 + (-t20 * t60 - t205 * t3 + t206 * t4 + t21 * t61) * mrSges(7,3) + (-Ifges(7,1) * t206 - Ifges(7,4) * t205) * t617 + (-Ifges(7,4) * t206 - Ifges(7,2) * t205) * t616 + (-Ifges(7,5) * t206 - Ifges(7,6) * t205) * t608 + t33 * (mrSges(7,1) * t205 - mrSges(7,2) * t206) + m(4) * (t181 * t286 + t182 * t285 + t215 * t267 + t216 * t266 - t319 * t392 - t374 * t396) + (t323 * t587 + t451 * qJD(2) / 0.2e1 - t636) * qJD(2) + t453 * t664 + (t490 + t483) * t269 / 0.2e1 + m(7) * (t110 * t58 + t137 * t33 + t20 * t6 + t21 * t5 + t3 * t35 + t34 * t4) + (-mrSges(4,2) * t392 + Ifges(4,1) * t348 + Ifges(4,4) * t346) * t243 + (-t15 * t534 - t16 * t533 + t438 * t82 - t439 * t83) * mrSges(6,3) + t265 * t396 + (t359 * t573 + t360 * t574 + t633) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t633) + (t662 / 0.2e1 + t661 / 0.2e1) * t406 + t239 * (-Ifges(6,4) * t438 - Ifges(6,2) * t439) / 0.2e1 + (t181 * t346 - t182 * t348 - t266 * t277 + t267 * t278) * mrSges(4,3) + (t207 * mrSges(5,2) - t53 * mrSges(5,3) + Ifges(5,1) * t175 + Ifges(5,4) * t174 + Ifges(5,5) * t629 + t130 * t473 + t450 * t603 + t452 * t609 + t454 * t610 + t48 * t456) * t270 + (t207 * mrSges(5,1) - t54 * mrSges(5,3) - Ifges(5,4) * t175 + Ifges(6,5) * t610 + Ifges(7,5) * t617 - Ifges(5,2) * t174 - Ifges(5,6) * t629 + Ifges(6,6) * t609 + Ifges(7,6) * t616 + Ifges(6,3) * t603 + Ifges(7,3) * t608 + t624 + t625) * t269 + (Ifges(4,5) * t277 + Ifges(4,6) * t278) * t588 + (Ifges(4,1) * t277 + Ifges(4,4) * t278) * t591 + (Ifges(7,3) * t597 + Ifges(6,5) * t599 + Ifges(7,5) * t604 + Ifges(7,6) * t606 + t659 / 0.2e1 + t660 / 0.2e1 - t665 + t675 / 0.2e1 + t686) * t208 + (-m(5) * t143 + m(6) * t140 - t646) * t86 + (mrSges(4,1) * t392 + Ifges(4,4) * t348 + Ifges(4,2) * t346) * t244 + m(5) * (t144 * t87 + t194 * t54 + t207 * t300 + t261 * t281) + (-Ifges(6,1) * t438 - Ifges(6,4) * t439) * t599 + m(6) * (t15 * t97 + t16 * t96 + t22 * t83 + t23 * t82) + t252 * (-Ifges(6,5) * t438 - Ifges(6,6) * t439) / 0.2e1 + (Ifges(7,1) * t60 + Ifges(7,4) * t61) * t604 + t360 * t557 / 0.2e1 + (Ifges(3,4) * t360 + Ifges(3,2) * t359) * t587 + (Ifges(7,4) * t60 + Ifges(7,2) * t61) * t606 + (-m(5) * t53 + m(6) * t48 - t148 + t55) * t193 + (t661 + t662) * t590 + (Ifges(5,5) * t588 + t477 + t655 / 0.2e1 + t197 / 0.2e1 + t677) * t209 + (Ifges(7,5) * t60 + Ifges(7,6) * t61) * t597 + Ifges(2,3) * qJDD(1) - t372 * t539 - t44 * t534 / 0.2e1 - t322 * t501 / 0.2e1 - t440 * t495 + (t421 * t557 + t437) * t470 - t374 * (-mrSges(4,1) * t278 + mrSges(4,2) * t277) - pkin(1) * (-mrSges(3,1) * t359 + mrSges(3,2) * t360) + t319 * (-mrSges(4,1) * t346 + mrSges(4,2) * t348) + t300 * t469 + t324 * (Ifges(4,4) * t277 + Ifges(4,2) * t278) / 0.2e1 + t285 * t220 + t286 * t221 + t215 * t293 + t216 * t294 + t278 * t253 / 0.2e1 + t277 * t254 / 0.2e1 + t261 * t202 + t87 * t245 + t194 * t147 + t22 * t185 + t23 * t186 + t137 * t13 + t5 * t123 + t6 * t124 + t110 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) + t96 * t72 + t97 * t73 + t58 * t88 + t140 * (mrSges(6,1) * t439 - mrSges(6,2) * t438) + t60 * t611 + t61 * t613 + t533 * t615 - t206 * t618 - t205 * t619 + t34 * t30 + t35 * t31; (-t294 * t500 + m(4) * (t181 * t415 + t182 * t420 + (-t266 * t415 + t267 * t420) * qJD(3)) + t293 * t499 + t415 * t221) * pkin(2) - m(4) * (t266 * t273 + t267 * t274 - t374 * t395) + (t638 * t417 + t634) * g(2) + (t638 * t422 + t635) * g(1) + (t636 + (t440 - t437 / 0.2e1) * qJD(1)) * qJD(1) - (-Ifges(3,2) * t503 + t323 + t394) * t502 / 0.2e1 + (-m(4) * t405 - m(5) * t504 - m(6) * (t405 + t482) + t372 - m(7) * (t405 + t461) + t626) * g(3) + t668 * (m(4) * t584 - m(5) * t361 + t460 + t673) - t641 * t646 + t220 * t583 + t425 + (t309 * t48 - t82 * t91 - t83 * t92 - g(1) * (-t422 * t578 + t339 + t369) - g(2) * (-t417 * t578 + t338 + t368) + t641 * t140) * m(6) + (-t143 * t641 - t144 * t173 - t281 * t292 + t314 * t53 + t315 * t54) * m(5) + t643 * t88 + (m(5) * t144 + m(6) * t449 - t647) * t313 + t651 * t124 + t652 * t123 + (-g(1) * t339 - g(2) * t338 + t110 * t643 + t20 * t651 + t21 * t652 + t212 * t4 + t213 * t3 + t299 * t33) * m(7) + t623 * t310 + t322 * t503 / 0.2e1 - t451 * t495 / 0.2e1 - t265 * t395 + Ifges(3,6) * t359 + Ifges(3,5) * t360 - t341 * mrSges(3,2) - t342 * mrSges(3,1) + t309 * t55 + Ifges(3,3) * qJDD(2) + t314 * t148 + t315 * t147 + t299 * t13 - t292 * t202 - t274 * t293 - t273 * t294 - t173 * t245 + t212 * t30 + t213 * t31 - t92 * t185 - t91 * t186 + t73 * t532; t637 * t88 + (t417 * t639 + t634) * g(2) + (t422 * t639 + t635) * g(1) + (-m(5) * t389 - m(6) * t482 - m(7) * t461 + t626) * g(3) + (m(5) * t581 + t673) * t668 + t148 * t579 + t147 * t580 + t425 + t646 * t150 + ((t411 * t54 + t412 * t53) * pkin(3) + t143 * t150 - t144 * t151 - t281 * t582) * m(5) + (-g(1) * t369 - g(2) * t368 - t140 * t150 + t388 * t48 - t82 * t89 - t83 * t90) * m(6) + t623 * t387 + t649 * t124 + (t110 * t637 + t20 * t649 + t21 * t650 + t262 * t4 + t263 * t3 + t33 * t367) * m(7) + t650 * t123 - t202 * t582 + t73 * t526 + t388 * t55 + t367 * t13 - t266 * t293 + t267 * t294 + t262 * t30 + t263 * t31 - t151 * t245 - t90 * t185 - t89 * t186; -t444 * t30 + t347 * t31 + t414 * t73 + t419 * t72 - t674 * t124 + t645 * t123 + t448 * qJD(5) + t647 * t642 + (-t88 + t646) * t447 + t469 + (-g(1) * t417 + g(2) * t422) * t631 + (-t110 * t447 - t20 * t674 + t21 * t645 + t3 * t347 - t4 * t444) * m(7) + (-t140 * t447 + t15 * t414 + t16 * t419 + t252 * t449) * m(6) + (t143 * t447 - t144 * t642 + t207) * m(5); (t418 * t552 - t25) * t123 + (Ifges(6,5) * t239 - Ifges(6,6) * t240) * t596 - (mrSges(7,1) * t110 + Ifges(7,4) * t605 + Ifges(7,2) * t607 + Ifges(7,6) * t598 - t585 + t614) * t165 + (-mrSges(7,2) * t110 + Ifges(7,1) * t605 + Ifges(7,4) * t607 + Ifges(7,5) * t598 + t586 + t612) * t465 + (t30 * t418 + t31 * t413) * pkin(5) + (-t413 * t552 - t24) * t124 + (t560 - t185) * t82 + (-mrSges(6,2) * t306 + t305 * t653 - t510) * g(2) + t624 + (-t455 + t456 + t663) * t570 + (-Ifges(6,2) * t240 + t130 + t233) * t601 - t88 * t577 - m(7) * (t110 * t577 + t20 * t24 + t21 * t25) + (mrSges(6,2) * t308 - t307 * t653 - t509) * g(1) + t483 + (t559 + t186) * t83 - t140 * (mrSges(6,1) * t240 + mrSges(6,2) * t239) + t129 * t599 + (Ifges(6,1) * t239 - t548) * t600 + (t3 * t413 + t4 * t418 + (-t20 * t413 + t21 * t418) * qJD(6)) * t620 + t441; -t110 * (mrSges(7,1) * t165 + mrSges(7,2) * t465) + (Ifges(7,1) * t465 - t554) * t605 + t76 * t604 + (Ifges(7,5) * t465 - Ifges(7,6) * t165) * t598 - t20 * t123 + t21 * t124 - g(1) * t509 - g(2) * t510 - t455 * t570 + (t165 * t21 + t20 * t465) * mrSges(7,3) + t441 + (-Ifges(7,2) * t165 + t155 + t77) * t607;];
tau  = t1;
