% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:56
% EndTime: 2019-03-09 16:49:02
% DurationCPUTime: 42.68s
% Computational Cost: add. (17467->894), mult. (38952->1149), div. (0->0), fcn. (28094->14), ass. (0->403)
t366 = sin(qJ(2));
t369 = cos(qJ(2));
t415 = pkin(2) * t366 - pkin(8) * t369;
t306 = t415 * qJD(1);
t368 = cos(qJ(3));
t365 = sin(qJ(3));
t461 = qJD(1) * t366;
t438 = t365 * t461;
t213 = pkin(7) * t438 + t368 * t306;
t473 = t368 * t369;
t397 = pkin(3) * t366 - qJ(4) * t473;
t363 = -qJ(4) - pkin(8);
t422 = qJD(3) * t363;
t655 = -qJD(1) * t397 - qJD(4) * t365 + t368 * t422 - t213;
t278 = t365 * t306;
t453 = qJD(4) * t368;
t478 = t366 * t368;
t480 = t365 * t369;
t654 = t278 + (-pkin(7) * t478 - qJ(4) * t480) * qJD(1) - t365 * t422 - t453;
t361 = sin(pkin(10));
t362 = cos(pkin(10));
t297 = t361 * t368 + t362 * t365;
t385 = t297 * t369;
t235 = qJD(1) * t385;
t269 = t297 * qJD(3);
t653 = t235 - t269;
t643 = mrSges(6,1) + mrSges(7,1);
t642 = mrSges(6,2) - mrSges(7,3);
t598 = t654 * t361 + t362 * t655;
t597 = t361 * t655 - t654 * t362;
t458 = qJD(2) * t368;
t303 = -t438 + t458;
t348 = pkin(7) * t461;
t326 = -qJD(2) * pkin(2) + t348;
t211 = -t303 * pkin(3) + qJD(4) + t326;
t416 = pkin(2) * t369 + pkin(8) * t366;
t316 = -pkin(1) - t416;
t285 = t316 * qJD(1);
t460 = qJD(1) * t369;
t349 = pkin(7) * t460;
t327 = qJD(2) * pkin(8) + t349;
t204 = t368 * t285 - t327 * t365;
t436 = t368 * t461;
t304 = qJD(2) * t365 + t436;
t163 = -qJ(4) * t304 + t204;
t336 = qJD(3) - t460;
t151 = pkin(3) * t336 + t163;
t205 = t285 * t365 + t327 * t368;
t164 = qJ(4) * t303 + t205;
t157 = t361 * t164;
t91 = t362 * t151 - t157;
t652 = -t211 * mrSges(5,2) + mrSges(5,3) * t91;
t651 = mrSges(4,3) + mrSges(5,3);
t650 = t653 * pkin(9) + t597;
t398 = t361 * t365 - t362 * t368;
t384 = t398 * t369;
t236 = qJD(1) * t384;
t270 = t398 * qJD(3);
t649 = -pkin(4) * t461 + t598 + (-t236 + t270) * pkin(9);
t615 = -Ifges(6,4) + Ifges(7,5);
t648 = t615 + Ifges(7,5);
t450 = qJD(1) * qJD(2);
t308 = qJDD(1) * t369 - t366 * t450;
t293 = t308 * pkin(7);
t647 = t293 * t369;
t198 = t303 * t362 - t304 * t361;
t199 = t303 * t361 + t304 * t362;
t364 = sin(qJ(5));
t525 = cos(qJ(5));
t137 = -t525 * t198 + t364 * t199;
t626 = t198 * t364 + t199 * t525;
t76 = pkin(5) * t626 + qJ(6) * t137;
t309 = qJDD(1) * t366 + t369 * t450;
t191 = qJD(3) * t303 + qJDD(2) * t365 + t309 * t368;
t192 = -qJD(3) * t304 + qJDD(2) * t368 - t309 * t365;
t125 = -t191 * t361 + t192 * t362;
t127 = t191 * t362 + t192 * t361;
t42 = -qJD(5) * t137 + t364 * t125 + t127 * t525;
t564 = t42 / 0.2e1;
t43 = qJD(5) * t626 - t125 * t525 + t364 * t127;
t562 = t43 / 0.2e1;
t645 = m(6) + m(7);
t295 = qJDD(3) - t308;
t284 = qJDD(5) + t295;
t534 = t284 / 0.2e1;
t641 = -mrSges(6,3) - mrSges(7,2);
t616 = Ifges(6,1) + Ifges(7,1);
t614 = Ifges(7,4) + Ifges(6,5);
t613 = -Ifges(6,6) + Ifges(7,6);
t640 = -Ifges(5,3) - Ifges(4,3);
t612 = Ifges(6,3) + Ifges(7,2);
t437 = t365 * t460;
t456 = qJD(3) * t365;
t593 = -t349 + (-t437 + t456) * pkin(3);
t639 = -m(4) * pkin(8) + m(5) * t363 - t651;
t146 = -pkin(4) * t198 + t211;
t328 = qJD(5) + t336;
t635 = pkin(9) * t199;
t79 = pkin(4) * t336 - t635 + t91;
t513 = pkin(9) * t198;
t482 = t362 * t164;
t92 = t361 * t151 + t482;
t80 = t92 + t513;
t25 = -t364 * t80 + t525 * t79;
t633 = qJD(6) - t25;
t23 = -t328 * pkin(5) + t633;
t51 = t137 * pkin(5) - qJ(6) * t626 + t146;
t131 = Ifges(6,4) * t137;
t498 = Ifges(7,5) * t137;
t606 = t328 * t614 + t616 * t626 - t131 + t498;
t638 = mrSges(6,2) * t146 - mrSges(7,3) * t51 - mrSges(6,3) * t25 + mrSges(7,2) * t23 + t606 / 0.2e1;
t26 = t364 * t79 + t525 * t80;
t24 = t328 * qJ(6) + t26;
t130 = Ifges(7,5) * t626;
t63 = Ifges(7,6) * t328 + Ifges(7,3) * t137 + t130;
t499 = Ifges(6,4) * t626;
t66 = -Ifges(6,2) * t137 + Ifges(6,6) * t328 + t499;
t637 = -mrSges(7,2) * t24 - mrSges(6,3) * t26 + mrSges(6,1) * t146 + mrSges(7,1) * t51 + t63 / 0.2e1 - t66 / 0.2e1;
t360 = qJ(3) + pkin(10);
t354 = qJ(5) + t360;
t340 = sin(t354);
t341 = cos(t354);
t514 = pkin(3) * t368;
t344 = pkin(2) + t514;
t352 = sin(t360);
t353 = cos(t360);
t412 = -mrSges(4,1) * t368 + mrSges(4,2) * t365;
t636 = m(4) * pkin(2) + m(5) * t344 + t353 * mrSges(5,1) - t352 * mrSges(5,2) - t642 * t340 + t341 * t643 - t412;
t535 = t199 / 0.2e1;
t320 = t363 * t365;
t322 = t363 * t368;
t209 = t362 * t320 + t322 * t361;
t177 = -pkin(9) * t297 + t209;
t210 = t361 * t320 - t362 * t322;
t178 = -pkin(9) * t398 + t210;
t106 = t364 * t177 + t178 * t525;
t608 = -qJD(5) * t106 - t364 * t650 + t525 * t649;
t390 = t177 * t525 - t364 * t178;
t607 = qJD(5) * t390 + t364 * t649 + t525 * t650;
t634 = -m(7) * qJ(6) - mrSges(7,3);
t516 = pkin(3) * t362;
t342 = pkin(4) + t516;
t96 = -t163 * t361 - t482;
t396 = t96 - t513;
t431 = qJD(5) * t525;
t517 = pkin(3) * t361;
t447 = t364 * t517;
t97 = t362 * t163 - t157;
t88 = t97 - t635;
t601 = -qJD(5) * t447 + t342 * t431 - t364 * t396 - t525 * t88;
t505 = mrSges(6,3) * t626;
t114 = mrSges(6,1) * t328 - t505;
t115 = -mrSges(7,1) * t328 + mrSges(7,2) * t626;
t470 = t114 - t115;
t596 = -pkin(4) * t653 + t593;
t454 = qJD(3) * t368;
t457 = qJD(2) * t369;
t380 = t365 * t457 + t366 * t454;
t367 = sin(qJ(1));
t370 = cos(qJ(1));
t632 = g(1) * t370 + g(2) * t367;
t490 = qJDD(1) * pkin(1);
t206 = -pkin(2) * t308 - pkin(8) * t309 - t490;
t264 = qJDD(2) * pkin(8) + t293;
t111 = -qJD(3) * t205 + t368 * t206 - t264 * t365;
t75 = pkin(3) * t295 - qJ(4) * t191 - qJD(4) * t304 + t111;
t110 = t365 * t206 + t368 * t264 + t285 * t454 - t327 * t456;
t83 = qJ(4) * t192 + qJD(4) * t303 + t110;
t27 = -t361 * t83 + t362 * t75;
t19 = pkin(4) * t295 - pkin(9) * t127 + t27;
t28 = t361 * t75 + t362 * t83;
t22 = pkin(9) * t125 + t28;
t6 = -qJD(5) * t26 + t19 * t525 - t364 * t22;
t3 = -t284 * pkin(5) + qJDD(6) - t6;
t563 = -t43 / 0.2e1;
t294 = t309 * pkin(7);
t265 = -qJDD(2) * pkin(2) + t294;
t150 = -pkin(3) * t192 + qJDD(4) + t265;
t87 = -pkin(4) * t125 + t150;
t7 = pkin(5) * t43 - qJ(6) * t42 - qJD(6) * t626 + t87;
t630 = mrSges(6,2) * t87 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t7 + Ifges(6,4) * t563 + 0.2e1 * t534 * t614 + t562 * t648 + 0.2e1 * t564 * t616;
t529 = t328 / 0.2e1;
t542 = t626 / 0.2e1;
t545 = t137 / 0.2e1;
t546 = -t137 / 0.2e1;
t629 = -Ifges(6,2) * t546 + Ifges(7,3) * t545 + t529 * t613 + t542 * t615 + t637;
t530 = -t328 / 0.2e1;
t543 = -t626 / 0.2e1;
t628 = -Ifges(6,2) * t545 + Ifges(7,3) * t546 + t530 * t613 + t543 * t615 - t637;
t624 = Ifges(6,4) * t545 + Ifges(7,5) * t546 + t530 * t614 + t543 * t616 - t638;
t623 = Ifges(6,4) * t546 + Ifges(7,5) * t545 + t529 * t614 + t542 * t616 + t638;
t567 = m(5) * pkin(3);
t549 = t125 / 0.2e1;
t548 = t127 / 0.2e1;
t539 = t191 / 0.2e1;
t538 = t192 / 0.2e1;
t621 = t198 / 0.2e1;
t533 = t295 / 0.2e1;
t620 = t308 / 0.2e1;
t619 = t309 / 0.2e1;
t351 = pkin(7) * t457;
t617 = -mrSges(3,3) + mrSges(2,2);
t610 = -qJ(6) * t461 + t607;
t609 = pkin(5) * t461 - t608;
t388 = -t364 * t297 - t398 * t525;
t132 = qJD(5) * t388 - t364 * t269 - t270 * t525;
t194 = t297 * t525 - t364 * t398;
t133 = qJD(5) * t194 + t269 * t525 - t364 * t270;
t155 = t235 * t525 - t236 * t364;
t156 = -t364 * t235 - t236 * t525;
t605 = -qJD(6) * t194 + t596 + (-t132 + t156) * qJ(6) + (t133 - t155) * pkin(5);
t604 = t198 * Ifges(5,6);
t602 = qJD(6) + t601;
t600 = t567 + mrSges(4,1);
t299 = t368 * t316;
t201 = -qJ(4) * t478 + t299 + (-pkin(7) * t365 - pkin(3)) * t369;
t338 = pkin(7) * t473;
t238 = t365 * t316 + t338;
t481 = t365 * t366;
t208 = -qJ(4) * t481 + t238;
t142 = t362 * t201 - t208 * t361;
t254 = t398 * t366;
t102 = -pkin(4) * t369 + pkin(9) * t254 + t142;
t143 = t361 * t201 + t362 * t208;
t253 = t297 * t366;
t109 = -pkin(9) * t253 + t143;
t599 = t364 * t102 + t525 * t109;
t292 = Ifges(4,4) * t303;
t182 = t304 * Ifges(4,1) + t336 * Ifges(4,5) + t292;
t347 = Ifges(3,4) * t460;
t595 = Ifges(3,1) * t461 + Ifges(3,5) * qJD(2) + t368 * t182 + t347;
t594 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t303 - mrSges(4,2) * t304 - mrSges(3,3) * t461;
t259 = t364 * t342 + t525 * t517;
t592 = Ifges(4,5) * t191 + Ifges(5,5) * t127 + Ifges(4,6) * t192 + Ifges(5,6) * t125 - t295 * t640;
t591 = t284 * t612 + t42 * t614 + t43 * t613;
t472 = t369 * t370;
t231 = t340 * t472 - t367 * t341;
t232 = t367 * t340 + t341 * t472;
t590 = t231 * t643 + t232 * t642;
t475 = t367 * t369;
t229 = t340 * t475 + t341 * t370;
t230 = -t340 * t370 + t341 * t475;
t589 = t229 * t643 + t230 * t642;
t588 = t294 * t366 + t647;
t587 = t110 * t368 - t111 * t365;
t586 = t641 * t366;
t585 = -m(4) - m(5) - m(3);
t584 = t304 * Ifges(4,5) + Ifges(5,5) * t199 + t303 * Ifges(4,6) + t137 * t613 + t328 * t612 - t336 * t640 + t614 * t626 + t604;
t172 = -qJD(2) * t384 - t269 * t366;
t459 = qJD(2) * t366;
t307 = t415 * qJD(2);
t444 = pkin(7) * t459;
t463 = t368 * t307 + t365 * t444;
t129 = -t366 * t453 + t397 * qJD(2) + (-t338 + (qJ(4) * t366 - t316) * t365) * qJD(3) + t463;
t464 = t365 * t307 + t316 * t454;
t140 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t478 + (-qJD(4) * t366 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t369) * t365 + t464;
t73 = t362 * t129 - t140 * t361;
t53 = pkin(4) * t459 - pkin(9) * t172 + t73;
t455 = qJD(3) * t366;
t171 = -qJD(2) * t385 + t398 * t455;
t74 = t361 * t129 + t362 * t140;
t56 = pkin(9) * t171 + t74;
t11 = -qJD(5) * t599 - t364 * t56 + t525 * t53;
t583 = m(7) * pkin(5) + t643;
t414 = t369 * mrSges(3,1) - mrSges(3,2) * t366;
t582 = t366 * t651 + mrSges(2,1) + t414;
t581 = -mrSges(6,2) - t634;
t452 = qJD(5) * t364;
t5 = t364 * t19 + t525 * t22 + t79 * t431 - t452 * t80;
t2 = qJ(6) * t284 + qJD(6) * t328 + t5;
t577 = -t6 * mrSges(6,1) + t3 * mrSges(7,1) + t5 * mrSges(6,2) - t2 * mrSges(7,3);
t576 = -t111 * mrSges(4,1) - t27 * mrSges(5,1) + t110 * mrSges(4,2) + t28 * mrSges(5,2);
t504 = Ifges(3,4) * t366;
t406 = t369 * Ifges(3,2) + t504;
t573 = t23 * mrSges(7,1) + t26 * mrSges(6,2) + t92 * mrSges(5,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t406 / 0.2e1 - t604 / 0.2e1 - t24 * mrSges(7,3) - t25 * mrSges(6,1) - t91 * mrSges(5,1);
t569 = mrSges(6,1) * t87 + mrSges(7,1) * t7 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t562 - t42 * Ifges(6,4) / 0.2e1 - t284 * Ifges(6,6) / 0.2e1 + t648 * t564 + (t613 + Ifges(7,6)) * t534 + (-t563 + t562) * Ifges(6,2);
t561 = Ifges(5,4) * t548 + Ifges(5,2) * t549 + Ifges(5,6) * t533;
t560 = Ifges(5,1) * t548 + Ifges(5,4) * t549 + Ifges(5,5) * t533;
t554 = Ifges(4,1) * t539 + Ifges(4,4) * t538 + Ifges(4,5) * t533;
t123 = Ifges(5,4) * t199 + Ifges(5,2) * t198 + Ifges(5,6) * t336;
t553 = -t123 / 0.2e1;
t552 = t123 / 0.2e1;
t500 = Ifges(5,4) * t198;
t124 = Ifges(5,1) * t199 + Ifges(5,5) * t336 + t500;
t551 = -t124 / 0.2e1;
t550 = t124 / 0.2e1;
t536 = -t199 / 0.2e1;
t531 = t304 / 0.2e1;
t528 = -t336 / 0.2e1;
t527 = t336 / 0.2e1;
t521 = mrSges(5,3) * t92;
t518 = pkin(3) * t304;
t515 = pkin(3) * t365;
t510 = g(3) * t366;
t355 = t366 * pkin(7);
t507 = mrSges(6,2) * t341;
t506 = mrSges(6,3) * t137;
t503 = Ifges(3,4) * t369;
t502 = Ifges(4,4) * t365;
t501 = Ifges(4,4) * t368;
t497 = t204 * mrSges(4,3);
t496 = t205 * mrSges(4,3);
t495 = t304 * Ifges(4,4);
t359 = -pkin(9) + t363;
t483 = t359 * t366;
t479 = t365 * t370;
t477 = t366 * t370;
t476 = t367 * t365;
t312 = pkin(4) * t353 + t514;
t305 = pkin(2) + t312;
t271 = t369 * t305;
t311 = pkin(4) * t352 + t515;
t290 = t370 * t311;
t112 = -mrSges(7,2) * t137 + mrSges(7,3) * t328;
t113 = -mrSges(6,2) * t328 - t506;
t471 = -t112 - t113;
t465 = -t369 * t290 + t367 * t312;
t310 = pkin(3) * t481 + t355;
t462 = t370 * pkin(1) + t367 * pkin(7);
t228 = pkin(3) * t380 + t351;
t181 = t303 * Ifges(4,2) + t336 * Ifges(4,6) + t495;
t433 = -t365 * t181 / 0.2e1;
t17 = t43 * mrSges(6,1) + t42 * mrSges(6,2);
t16 = t43 * mrSges(7,1) - t42 * mrSges(7,3);
t424 = t634 * t341 * t366;
t32 = -t284 * mrSges(7,1) + t42 * mrSges(7,2);
t421 = t450 / 0.2e1;
t69 = -t125 * mrSges(5,1) + t127 * mrSges(5,2);
t420 = -t229 * pkin(5) + t230 * qJ(6);
t419 = -t231 * pkin(5) + qJ(6) * t232;
t418 = -t311 * t475 - t312 * t370;
t202 = pkin(4) * t253 + t310;
t154 = pkin(4) * t199 + t518;
t413 = mrSges(3,1) * t366 + mrSges(3,2) * t369;
t411 = mrSges(4,1) * t365 + mrSges(4,2) * t368;
t408 = Ifges(4,1) * t368 - t502;
t407 = Ifges(4,1) * t365 + t501;
t405 = -Ifges(4,2) * t365 + t501;
t404 = Ifges(4,2) * t368 + t502;
t403 = Ifges(3,5) * t369 - Ifges(3,6) * t366;
t402 = Ifges(4,5) * t368 - Ifges(4,6) * t365;
t401 = Ifges(4,5) * t365 + Ifges(4,6) * t368;
t400 = pkin(5) * t341 + qJ(6) * t340;
t399 = t344 * t369 - t363 * t366;
t239 = pkin(4) * t398 - t344;
t144 = -pkin(4) * t171 + t228;
t392 = pkin(1) * t413;
t276 = -t365 * t472 + t367 * t368;
t274 = t365 * t475 + t368 * t370;
t59 = t102 * t525 - t364 * t109;
t389 = -t253 * t525 + t364 * t254;
t170 = -t364 * t253 - t254 * t525;
t10 = t102 * t431 - t109 * t452 + t364 * t53 + t525 * t56;
t387 = t326 * t411;
t386 = t366 * (Ifges(3,1) * t369 - t504);
t258 = t342 * t525 - t447;
t381 = -t365 * t455 + t368 * t457;
t378 = Ifges(4,5) * t366 + t369 * t408;
t377 = Ifges(4,6) * t366 + t369 * t405;
t376 = Ifges(4,3) * t366 + t369 * t402;
t375 = -t577 + t591;
t357 = t370 * pkin(7);
t318 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t460;
t286 = t411 * t366;
t277 = t368 * t472 + t476;
t275 = -t367 * t473 + t479;
t255 = -pkin(5) - t258;
t252 = qJ(6) + t259;
t249 = t367 * t352 + t353 * t472;
t248 = -t352 * t472 + t367 * t353;
t247 = t352 * t370 - t353 * t475;
t246 = t352 * t475 + t353 * t370;
t237 = -pkin(7) * t480 + t299;
t225 = mrSges(4,1) * t336 - mrSges(4,3) * t304;
t224 = -mrSges(4,2) * t336 + mrSges(4,3) * t303;
t214 = -pkin(7) * t436 + t278;
t168 = mrSges(5,1) * t336 - mrSges(5,3) * t199;
t167 = -t336 * mrSges(5,2) + mrSges(5,3) * t198;
t162 = -qJD(3) * t238 + t463;
t161 = (-t366 * t458 - t369 * t456) * pkin(7) + t464;
t153 = -mrSges(4,2) * t295 + mrSges(4,3) * t192;
t152 = mrSges(4,1) * t295 - mrSges(4,3) * t191;
t141 = -mrSges(5,1) * t198 + t199 * mrSges(5,2);
t134 = -mrSges(4,1) * t192 + mrSges(4,2) * t191;
t107 = t191 * Ifges(4,4) + t192 * Ifges(4,2) + t295 * Ifges(4,6);
t103 = -pkin(5) * t388 - qJ(6) * t194 + t239;
t100 = mrSges(5,1) * t295 - mrSges(5,3) * t127;
t99 = -mrSges(5,2) * t295 + mrSges(5,3) * t125;
t89 = -pkin(5) * t389 - qJ(6) * t170 + t202;
t86 = qJD(5) * t170 - t171 * t525 + t364 * t172;
t85 = qJD(5) * t389 + t364 * t171 + t172 * t525;
t78 = mrSges(6,1) * t137 + mrSges(6,2) * t626;
t77 = mrSges(7,1) * t137 - mrSges(7,3) * t626;
t58 = t369 * pkin(5) - t59;
t57 = -qJ(6) * t369 + t599;
t55 = t154 + t76;
t34 = -mrSges(7,2) * t43 + mrSges(7,3) * t284;
t33 = -mrSges(6,2) * t284 - mrSges(6,3) * t43;
t31 = mrSges(6,1) * t284 - mrSges(6,3) * t42;
t20 = pkin(5) * t86 - qJ(6) * t85 - qJD(6) * t170 + t144;
t9 = -pkin(5) * t459 - t11;
t8 = qJ(6) * t459 - qJD(6) * t369 + t10;
t1 = [m(6) * (t10 * t26 + t11 * t25 + t144 * t146 + t202 * t87 + t5 * t599 + t59 * t6) + t599 * t33 - t569 * t389 + t414 * t490 + t623 * t85 + t629 * t86 + t630 * t170 + (t433 + t595 / 0.2e1) * t457 + t150 * (mrSges(5,1) * t253 - mrSges(5,2) * t254) + (-Ifges(5,4) * t254 - Ifges(5,2) * t253) * t549 + (t171 * t92 - t172 * t91 - t253 * t28 + t254 * t27) * mrSges(5,3) + (-Ifges(5,5) * t254 - Ifges(5,6) * t253) * t533 + (-Ifges(5,1) * t254 - Ifges(5,4) * t253) * t548 - (t181 * t368 + t182 * t365) * t455 / 0.2e1 + (m(4) * t265 * pkin(7) + Ifges(3,1) * t309 + Ifges(3,4) * t620 + t402 * t533 + t405 * t538 + t408 * t539) * t366 + (t584 / 0.2e1 + t204 * mrSges(4,1) - t205 * mrSges(4,2) + Ifges(5,5) * t535 + Ifges(6,6) * t546 + Ifges(7,6) * t545 + Ifges(5,3) * t527 + t529 * t612 + t542 * t614 - t573) * t459 - t392 * t450 + t503 * t619 + t406 * t620 + (t309 * t355 + t588 + t647) * mrSges(3,3) + t303 * (qJD(2) * t377 - t404 * t455) / 0.2e1 + (Ifges(5,4) * t172 + Ifges(5,2) * t171) * t621 + (-t110 * t481 - t111 * t478 - t204 * t381 - t205 * t380) * mrSges(4,3) - t594 * t351 + m(7) * (t2 * t57 + t20 * t51 + t23 * t9 + t24 * t8 + t3 * t58 + t7 * t89) + m(5) * (t142 * t27 + t143 * t28 + t150 * t310 + t211 * t228 + t73 * t91 + t74 * t92) + (-t476 * t567 - t277 * mrSges(4,1) - t249 * mrSges(5,1) - t276 * mrSges(4,2) - t248 * mrSges(5,2) + t641 * t477 + t585 * t462 - t645 * (t305 * t472 + t367 * t311 - t359 * t477 + t462) + t617 * t367 - t583 * t232 - t581 * t231 + (-m(4) * t416 - m(5) * t399 - t582) * t370) * g(2) + (-t479 * t567 - t275 * mrSges(4,1) - t247 * mrSges(5,1) - t274 * mrSges(4,2) - t246 * mrSges(5,2) - t645 * (t367 * t483 + t290 + t357) + t617 * t370 + t585 * t357 + t583 * t230 + t581 * t229 + (-m(4) * t316 - m(5) * (-pkin(1) - t399) + m(3) * pkin(1) - t645 * (-pkin(1) - t271) + t582 - t586) * t367) * g(1) + (-mrSges(3,1) * t355 + Ifges(3,5) * t366 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * t369) * qJDD(2) + (t576 - t614 * t564 + (-Ifges(3,2) * t366 + t503) * t421 - Ifges(4,6) * t538 - Ifges(4,5) * t539 - Ifges(5,5) * t548 - Ifges(5,6) * t549 - Ifges(7,6) * t562 - Ifges(6,6) * t563 + Ifges(3,4) * t619 + Ifges(3,2) * t620 - t612 * t534 + t640 * t533 + t577) * t369 + t74 * t167 + t73 * t168 + t386 * t421 + m(4) * (t110 * t238 + t111 * t237 + t161 * t205 + t162 * t204 + t326 * t351) - t107 * t481 / 0.2e1 + (Ifges(5,1) * t172 + Ifges(5,4) * t171) * t535 + (Ifges(5,5) * t172 + Ifges(5,6) * t171 + qJD(2) * t376 - t401 * t455) * t527 + t142 * t100 + t143 * t99 + t144 * t78 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t588) - (t591 + t592) * t369 / 0.2e1 + t8 * t112 + t10 * t113 + t11 * t114 + t9 * t115 + Ifges(2,3) * qJDD(1) + t89 * t16 + qJD(2) ^ 2 * t403 / 0.2e1 + t20 * t77 + t59 * t31 + t57 * t34 + t58 * t32 - t318 * t444 + t326 * (mrSges(4,1) * t380 + mrSges(4,2) * t381) + t202 * t17 + t211 * (-mrSges(5,1) * t171 + mrSges(5,2) * t172) + t161 * t224 + t162 * t225 + t228 * t141 + t237 * t152 + t238 * t153 + t134 * t355 + (qJD(2) * t378 - t407 * t455) * t531 + t172 * t550 + t171 * t552 + t478 * t554 - t254 * t560 - t253 * t561 + t265 * t286 - pkin(1) * (-mrSges(3,1) * t308 + mrSges(3,2) * t309) + t310 * t69; (-Ifges(5,1) * t236 - Ifges(5,4) * t235) * t536 - t198 * (-Ifges(5,4) * t236 - Ifges(5,2) * t235) / 0.2e1 + (t235 * t92 - t236 * t91 - t27 * t297 - t28 * t398) * mrSges(5,3) + (-Ifges(5,5) * t236 - Ifges(5,6) * t235) * t528 + t624 * t156 + (t106 * t5 + t146 * t596 + t239 * t87 + t25 * t608 + t26 * t607 + t390 * t6) * m(6) + (t103 * t7 + t106 * t2 + t23 * t609 + t24 * t610 - t3 * t390 + t51 * t605) * m(7) - (t32 - t31) * t390 - t569 * t388 + (Ifges(5,5) * t297 - Ifges(5,6) * t398 + t401) * t533 + (Ifges(5,1) * t297 - Ifges(5,4) * t398) * t548 + (Ifges(5,4) * t297 - Ifges(5,2) * t398) * t549 + t150 * (mrSges(5,1) * t398 + mrSges(5,2) * t297) - t398 * t561 + t623 * t132 + t628 * t155 + t629 * t133 + t630 * t194 - t584 * t461 / 0.2e1 + t605 * t77 + t607 * t113 + t608 * t114 + t609 * t115 + t610 * t112 + (t392 - t386 / 0.2e1) * qJD(1) ^ 2 + (t303 * t405 + t304 * t408 + t336 * t402) * qJD(3) / 0.2e1 - (t303 * t377 + t304 * t378 + t336 * t376) * qJD(1) / 0.2e1 + (Ifges(5,5) * t536 + Ifges(6,6) * t545 + Ifges(7,6) * t546 + Ifges(5,3) * t528 + t530 * t612 + t543 * t614 + t573) * t461 + (t182 / 0.2e1 - t497) * t454 - t403 * t450 / 0.2e1 + (t433 + t387) * qJD(3) + t181 * t437 / 0.2e1 + (-Ifges(5,4) * t535 - Ifges(5,2) * t621 - Ifges(5,6) * t527 - t521 - t552) * t269 + t593 * t141 + t594 * t349 - (-Ifges(3,2) * t461 + t347 + t595) * t460 / 0.2e1 + t596 * t78 + t597 * t167 + t598 * t168 + (-t150 * t344 + t209 * t27 + t210 * t28 + t211 * t593 + t597 * t92 + t598 * t91) * m(5) + (-t414 - t645 * (t271 - t483) + t639 * t366 + (-m(7) * t400 - t636) * t369 + t586) * g(3) + t632 * (t413 + (t359 * t645 + t639 + t641) * t369 + (m(6) * t305 - m(7) * (-t305 - t400) + t636) * t366) + t318 * t348 - t387 * t460 + (t33 + t34) * t106 - t456 * t496 + (-t204 * (mrSges(4,1) * t366 - mrSges(4,3) * t473) - t205 * (-mrSges(4,2) * t366 - mrSges(4,3) * t480)) * qJD(1) + (-Ifges(5,1) * t535 - Ifges(5,4) * t621 - Ifges(5,5) * t527 - t550 + t652) * t270 - pkin(2) * t134 + (m(4) * ((-t204 * t368 - t205 * t365) * qJD(3) + t587) - t225 * t454 - t224 * t456 + t368 * t153 - t365 * t152) * pkin(8) + t587 * mrSges(4,3) + t103 * t16 + t265 * t412 + Ifges(3,3) * qJDD(2) + t368 * t107 / 0.2e1 - t344 * t69 + t209 * t100 + t210 * t99 - t214 * t224 - t213 * t225 + t239 * t17 + t404 * t538 + t407 * t539 - t236 * t551 - t235 * t553 + t365 * t554 + t297 * t560 - t293 * mrSges(3,2) - t294 * mrSges(3,1) + Ifges(3,6) * t308 + Ifges(3,5) * t309 + (-mrSges(5,1) * t653 + mrSges(5,2) * t236) * t211 + (-pkin(2) * t265 - t204 * t213 - t205 * t214 - t326 * t349) * m(4); (t246 * mrSges(5,1) - t247 * mrSges(5,2) - m(7) * (t418 + t420) - m(6) * t418 - mrSges(4,2) * t275 + t600 * t274 + t589) * g(2) + (-t248 * mrSges(5,1) + t249 * mrSges(5,2) - m(7) * (t419 + t465) - m(6) * t465 + mrSges(4,2) * t277 - t600 * t276 + t590) * g(1) + t601 * t113 + t602 * t112 + t592 - (t211 * mrSges(5,1) + Ifges(5,4) * t536 + Ifges(5,6) * t528 - t521 + t553) * t199 - t576 + t628 * t626 + (-(m(7) * (-pkin(5) * t340 - t311) - t340 * mrSges(7,1)) * t366 + t424 + t286) * g(3) - (-Ifges(4,2) * t304 + t182 + t292) * t303 / 0.2e1 - t624 * t137 + t375 - t304 * (Ifges(4,1) * t303 - t495) / 0.2e1 + (-t146 * t154 + t258 * t6 + t259 * t5 + t26 * t601 + t311 * t510) * m(6) + (t2 * t252 + t24 * t602 + t255 * t3 - t51 * t55) * m(7) - t97 * t167 - t96 * t168 + (m(6) * t25 - m(7) * t23 + t470) * (-t259 * qJD(5) + t364 * t88 - t396 * t525) + (m(5) * t515 + mrSges(5,1) * t352 + mrSges(6,1) * t340 + mrSges(5,2) * t353 + t507) * t510 + (Ifges(5,2) * t535 - t500 / 0.2e1 + Ifges(5,5) * t528 + Ifges(5,1) * t536 + t551 + t652) * t198 - t154 * t78 - t55 * t77 - t141 * t518 - m(5) * (t211 * t518 + t91 * t96 + t92 * t97) - t204 * t224 + t205 * t225 + t252 * t34 + t255 * t32 + t258 * t31 + t259 * t33 + t304 * t496 + t303 * t497 + t100 * t516 + t99 * t517 + (Ifges(4,5) * t303 - Ifges(4,6) * t304) * t528 + t181 * t531 + (t27 * t362 + t28 * t361) * t567 - t326 * (mrSges(4,1) * t304 + mrSges(4,2) * t303); t470 * t626 - t471 * t137 - t198 * t167 + t199 * t168 + t16 + t17 + t69 + (t137 * t24 - t23 * t626 + t7) * m(7) + (t137 * t26 + t25 * t626 + t87) * m(6) + (-t198 * t92 + t199 * t91 + t150) * m(5) + (t369 * g(3) - t366 * t632) * (m(5) + t645); (t137 * t23 + t24 * t626) * mrSges(7,2) + t375 + t590 * g(1) + t589 * g(2) + ((t340 * t583 + t507) * t366 + t424) * g(3) + (t470 + t505) * t26 + (t471 - t506) * t25 - t146 * (mrSges(6,1) * t626 - mrSges(6,2) * t137) - t51 * (mrSges(7,1) * t626 + mrSges(7,3) * t137) + qJD(6) * t112 - t76 * t77 - pkin(5) * t32 + qJ(6) * t34 + t66 * t542 + (Ifges(7,3) * t626 - t498) * t546 + (-t137 * t614 + t613 * t626) * t530 + (-Ifges(6,2) * t626 - t131 + t606) * t545 + (-t137 * t616 + t130 - t499 + t63) * t543 + (-pkin(5) * t3 - t419 * g(1) - t420 * g(2) + qJ(6) * t2 - t23 * t26 + t24 * t633 - t51 * t76) * m(7); -t328 * t112 + t626 * t77 + (-g(1) * t231 - g(2) * t229 - t24 * t328 - t340 * t510 + t51 * t626 + t3) * m(7) + t32;];
tau  = t1;
