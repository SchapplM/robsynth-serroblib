% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:24
% EndTime: 2019-03-09 04:12:25
% DurationCPUTime: 35.50s
% Computational Cost: add. (36725->939), mult. (117512->1294), div. (0->0), fcn. (102878->18), ass. (0->431)
t368 = sin(pkin(7));
t374 = sin(qJ(3));
t369 = sin(pkin(6));
t511 = cos(pkin(12));
t512 = cos(pkin(7));
t436 = t512 * t511;
t535 = cos(qJ(3));
t403 = t535 * t436;
t398 = t369 * t403;
t513 = cos(pkin(6));
t443 = t513 * t535;
t412 = qJD(1) * t443;
t367 = sin(pkin(12));
t488 = qJD(1) * t369;
t469 = t367 * t488;
t275 = -qJD(1) * t398 - t368 * t412 + t374 * t469;
t366 = sin(pkin(13));
t370 = cos(pkin(13));
t373 = sin(qJ(5));
t376 = cos(qJ(5));
t335 = t366 * t376 + t370 * t373;
t194 = t335 * t275;
t327 = t335 * qJD(5);
t644 = t194 + t327;
t334 = t366 * t373 - t376 * t370;
t195 = t334 * t275;
t326 = t334 * qJD(5);
t643 = t195 + t326;
t459 = t367 * t513;
t357 = pkin(1) * t459;
t457 = t369 * t511;
t441 = qJD(1) * t457;
t315 = qJ(2) * t441 + qJD(1) * t357;
t455 = t513 * t368;
t386 = (t369 * t436 + t455) * pkin(9);
t261 = qJD(1) * t386 + t315;
t437 = t513 * t511;
t358 = pkin(1) * t437;
t352 = qJD(1) * t358;
t470 = t513 * pkin(2);
t497 = t367 * t369;
t389 = t470 + (-pkin(9) * t512 - qJ(2)) * t497;
t269 = qJD(1) * t389 + t352;
t307 = (-pkin(9) * t367 * t368 - pkin(2) * t511 - pkin(1)) * t369;
t297 = qJD(1) * t307 + qJD(2);
t442 = t512 * t535;
t474 = t368 * t535;
t170 = -t374 * t261 + t269 * t442 + t297 * t474;
t392 = t367 * t535 + t374 * t436;
t290 = t369 * t392 + t374 * t455;
t278 = t290 * qJD(1);
t203 = pkin(3) * t278 + qJ(4) * t275;
t117 = t370 * t170 + t366 * t203;
t505 = t275 * t366;
t111 = pkin(10) * t505 + t117;
t530 = pkin(10) + qJ(4);
t342 = t530 * t366;
t343 = t530 * t370;
t415 = -t376 * t342 - t343 * t373;
t116 = -t170 * t366 + t370 * t203;
t504 = t275 * t370;
t98 = pkin(4) * t278 + pkin(10) * t504 + t116;
t604 = -qJD(4) * t334 + qJD(5) * t415 - t376 * t111 - t373 * t98;
t456 = t374 * t512;
t304 = (-t367 * t456 + t511 * t535) * t488;
t449 = t368 * t469;
t252 = -t304 * t366 + t370 * t449;
t253 = t304 * t370 + t366 * t449;
t495 = t368 * t374;
t320 = -t366 * t495 + t370 * t512;
t321 = t366 * t512 + t370 * t495;
t416 = t376 * t320 - t321 * t373;
t471 = t376 * t535;
t472 = t373 * t535;
t487 = qJD(3) * t368;
t630 = -t252 * t373 - t253 * t376 + t416 * qJD(5) + (-t366 * t472 + t370 * t471) * t487;
t388 = t369 * (t367 * t442 + t374 * t511);
t303 = qJD(1) * t388;
t486 = qJD(3) * t374;
t467 = t368 * t486;
t642 = -t303 + t467;
t268 = qJD(5) + t275;
t319 = -t368 * t457 + t512 * t513;
t312 = qJD(1) * t319 + qJD(3);
t216 = t278 * t370 + t312 * t366;
t207 = -t269 * t368 + t512 * t297;
t143 = pkin(3) * t275 - qJ(4) * t278 + t207;
t254 = t269 * t456;
t171 = t535 * t261 + t297 * t495 + t254;
t145 = qJ(4) * t312 + t171;
t92 = t370 * t143 - t145 * t366;
t64 = pkin(4) * t275 - pkin(10) * t216 + t92;
t453 = -t278 * t366 + t370 * t312;
t93 = t366 * t143 + t370 * t145;
t66 = pkin(10) * t453 + t93;
t33 = -t373 * t66 + t376 * t64;
t28 = -pkin(5) * t268 - t33;
t622 = t216 * t376 + t373 * t453;
t527 = mrSges(6,3) * t622;
t125 = mrSges(6,1) * t268 - t527;
t372 = sin(qJ(6));
t375 = cos(qJ(6));
t122 = t268 * t375 - t372 * t622;
t123 = t268 * t372 + t375 * t622;
t67 = -mrSges(7,1) * t122 + mrSges(7,2) * t123;
t608 = t125 - t67;
t641 = -m(7) * t28 + t608;
t377 = cos(qJ(1));
t534 = sin(qJ(1));
t322 = t377 * t459 + t511 * t534;
t494 = t369 * t377;
t399 = t534 * t367 - t377 * t437;
t598 = t399 * t512;
t236 = -t322 * t535 + t374 * t598 + t494 * t495;
t458 = t369 * t512;
t294 = t399 * t368 - t377 * t458;
t365 = pkin(13) + qJ(5);
t362 = sin(t365);
t363 = cos(t365);
t640 = t236 * t363 - t294 * t362;
t639 = t236 * t362 + t294 * t363;
t431 = -mrSges(7,1) * t375 + mrSges(7,2) * t372;
t609 = m(7) * pkin(5) + mrSges(6,1) - t431;
t587 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t637 = -pkin(11) * t278 + t604;
t132 = -pkin(4) * t505 + t171;
t636 = pkin(5) * t644 + t643 * pkin(11) - t132;
t544 = t268 / 0.2e1;
t550 = t622 / 0.2e1;
t623 = -t216 * t373 + t376 * t453;
t552 = t623 / 0.2e1;
t144 = -t312 * pkin(3) + qJD(4) - t170;
t113 = -pkin(4) * t453 + t144;
t614 = t113 * mrSges(6,2);
t148 = Ifges(6,4) * t623;
t89 = Ifges(6,1) * t622 + t268 * Ifges(6,5) + t148;
t635 = Ifges(6,1) * t550 + Ifges(6,4) * t552 + Ifges(6,5) * t544 + t89 / 0.2e1 + t614;
t149 = qJD(6) - t623;
t277 = t290 * qJD(3);
t481 = qJDD(1) * t369;
t461 = t367 * t481;
t414 = t368 * t443;
t628 = -t398 - t414;
t206 = qJD(1) * t277 + qJDD(1) * t628 + t374 * t461;
t202 = qJDD(5) + t206;
t496 = t367 * t374;
t393 = t403 - t496;
t452 = qJDD(1) * t513;
t205 = (qJD(3) * t412 + t374 * t452) * t368 + (qJD(1) * qJD(3) * t393 + qJDD(1) * t392) * t369;
t311 = qJDD(1) * t319 + qJDD(3);
t178 = -t205 * t366 + t311 * t370;
t179 = t205 * t370 + t311 * t366;
t75 = qJD(5) * t623 + t178 * t373 + t179 * t376;
t45 = qJD(6) * t122 + t202 * t372 + t375 * t75;
t46 = -qJD(6) * t123 + t202 * t375 - t372 * t75;
t18 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t58 = mrSges(6,1) * t202 - mrSges(6,3) * t75;
t634 = t18 - t58;
t112 = -t178 * mrSges(5,1) + t179 * mrSges(5,2);
t76 = -qJD(5) * t622 + t178 * t376 - t179 * t373;
t36 = -t76 * mrSges(6,1) + t75 * mrSges(6,2);
t633 = -t112 - t36;
t299 = -t342 * t373 + t343 * t376;
t603 = -qJD(4) * t335 - qJD(5) * t299 + t111 * t373 - t376 * t98;
t258 = t320 * t373 + t321 * t376;
t400 = -t375 * t258 + t372 * t474;
t632 = qJD(6) * t400 - t372 * t630 + t375 * t642;
t239 = -t372 * t258 - t375 * t474;
t631 = qJD(6) * t239 + t372 * t642 + t375 * t630;
t595 = mrSges(4,1) * t312 + mrSges(5,1) * t453 - mrSges(5,2) * t216 - mrSges(4,3) * t278;
t394 = t377 * t367 + t437 * t534;
t381 = t394 * t368 + t534 * t458;
t180 = -mrSges(5,2) * t275 + mrSges(5,3) * t453;
t181 = mrSges(5,1) * t275 - mrSges(5,3) * t216;
t627 = t180 * t370 - t181 * t366;
t54 = t123 * Ifges(7,5) + t122 * Ifges(7,6) + t149 * Ifges(7,3);
t524 = Ifges(6,4) * t622;
t88 = Ifges(6,2) * t623 + t268 * Ifges(6,6) + t524;
t626 = -t54 / 0.2e1 + t88 / 0.2e1;
t473 = t369 * t534;
t624 = -g(1) * t473 + g(2) * t494 - g(3) * t513;
t233 = t322 * t374 + t474 * t494 + t535 * t598;
t571 = t45 / 0.2e1;
t570 = t46 / 0.2e1;
t73 = qJDD(6) - t76;
t565 = t73 / 0.2e1;
t564 = t75 / 0.2e1;
t563 = t76 / 0.2e1;
t621 = -m(7) - m(6);
t549 = t178 / 0.2e1;
t548 = t179 / 0.2e1;
t547 = t202 / 0.2e1;
t546 = t206 / 0.2e1;
t620 = t33 * mrSges(6,1);
t619 = t92 * mrSges(5,1);
t618 = t93 * mrSges(5,2);
t617 = Ifges(4,5) * t205;
t616 = Ifges(4,6) * t206;
t615 = Ifges(4,3) * t311;
t613 = t216 * Ifges(5,5);
t612 = t312 * Ifges(4,5);
t611 = t312 * Ifges(4,6);
t610 = t453 * Ifges(5,6);
t607 = Ifges(6,5) * t622 + Ifges(6,6) * t623 + t275 * Ifges(5,3) + t268 * Ifges(6,3) + t610 + t613;
t361 = pkin(4) * t370 + pkin(3);
t273 = pkin(5) * t334 - pkin(11) * t335 - t361;
t208 = t273 * t375 - t299 * t372;
t606 = qJD(6) * t208 + t372 * t636 + t375 * t637;
t209 = t273 * t372 + t299 * t375;
t605 = -qJD(6) * t209 - t372 * t637 + t375 * t636;
t602 = pkin(5) * t278 - t603;
t601 = t144 * (mrSges(5,1) * t366 + mrSges(5,2) * t370);
t160 = -t195 * t372 + t278 * t375;
t482 = qJD(6) * t375;
t409 = -t372 * t326 + t335 * t482;
t597 = t160 + t409;
t161 = t195 * t375 + t278 * t372;
t483 = qJD(6) * t372;
t408 = t375 * t326 + t335 * t483;
t596 = t161 + t408;
t231 = -t290 * t366 + t319 * t370;
t232 = t290 * t370 + t319 * t366;
t173 = t231 * t373 + t232 * t376;
t289 = t369 * t496 + t628;
t283 = t289 * t375;
t133 = -t173 * t372 + t283;
t594 = -t368 * t473 + t394 * t512;
t430 = mrSges(7,1) * t372 + mrSges(7,2) * t375;
t463 = m(5) * qJ(4) + mrSges(5,3);
t593 = -t430 - t463;
t592 = -t366 * t92 + t370 * t93;
t331 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t369;
t300 = qJDD(1) * t358 - t367 * t331;
t446 = t367 * t458;
t248 = (-pkin(9) * t446 + t470) * qJDD(1) + t300;
t295 = qJDD(1) * t307 + qJDD(2);
t196 = -t248 * t368 + t512 * t295;
t105 = pkin(3) * t206 - qJ(4) * t205 - qJD(4) * t278 + t196;
t301 = qJDD(1) * t357 + t511 * t331;
t244 = qJDD(1) * t386 + t301;
t411 = qJD(3) * t442;
t462 = qJD(3) * t535;
t447 = t368 * t462;
t106 = t535 * t244 + t248 * t456 - t261 * t486 + t269 * t411 + t295 * t495 + t297 * t447;
t91 = qJ(4) * t311 + qJD(4) * t312 + t106;
t47 = t370 * t105 - t366 * t91;
t48 = t366 * t105 + t370 * t91;
t591 = -t366 * t47 + t370 * t48;
t34 = t373 * t64 + t376 * t66;
t29 = pkin(11) * t268 + t34;
t57 = -pkin(5) * t623 - pkin(11) * t622 + t113;
t16 = -t29 * t372 + t375 * t57;
t107 = -qJD(3) * t254 - t374 * t244 + t248 * t442 - t261 * t462 + t295 * t474 - t297 * t467;
t96 = -t311 * pkin(3) + qJDD(4) - t107;
t65 = -t178 * pkin(4) + t96;
t21 = -t76 * pkin(5) - t75 * pkin(11) + t65;
t35 = pkin(4) * t206 - pkin(10) * t179 + t47;
t40 = pkin(10) * t178 + t48;
t484 = qJD(5) * t376;
t485 = qJD(5) * t373;
t7 = t373 * t35 + t376 * t40 + t64 * t484 - t485 * t66;
t5 = pkin(11) * t202 + t7;
t1 = qJD(6) * t16 + t21 * t372 + t375 * t5;
t17 = t29 * t375 + t372 * t57;
t2 = -qJD(6) * t17 + t21 * t375 - t372 * t5;
t590 = t1 * t375 - t2 * t372;
t589 = mrSges(4,2) + t593;
t434 = -mrSges(5,1) * t370 + mrSges(5,2) * t366;
t588 = m(5) * pkin(3) - t362 * t587 + t363 * t609 - t434;
t8 = -qJD(5) * t34 + t35 * t376 - t373 * t40;
t291 = t358 + t389;
t217 = -t291 * t368 + t512 * t307;
t159 = pkin(3) * t289 - qJ(4) * t290 + t217;
t325 = qJ(2) * t457 + t357;
t286 = t386 + t325;
t266 = t535 * t286;
t454 = t512 * t291;
t177 = t307 * t495 + t374 * t454 + t266;
t166 = qJ(4) * t319 + t177;
t108 = t370 * t159 - t166 * t366;
t79 = pkin(4) * t289 - pkin(10) * t232 + t108;
t109 = t366 * t159 + t370 * t166;
t86 = pkin(10) * t231 + t109;
t529 = t373 * t79 + t376 * t86;
t276 = (t369 * t393 + t414) * qJD(3);
t502 = t276 * t370;
t440 = qJD(2) * t457;
t155 = t291 * t411 + t307 * t447 + t535 * t440 + (-qJD(2) * t446 - qJD(3) * t286) * t374;
t139 = t319 * qJD(4) + t155;
t468 = qJD(2) * t497;
t448 = t368 * t468;
t174 = pkin(3) * t277 - qJ(4) * t276 - qJD(4) * t290 + t448;
t99 = -t139 * t366 + t370 * t174;
t74 = pkin(4) * t277 - pkin(10) * t502 + t99;
t100 = t370 * t139 + t366 * t174;
t503 = t276 * t366;
t85 = -pkin(10) * t503 + t100;
t15 = -qJD(5) * t529 - t373 * t85 + t376 * t74;
t586 = mrSges(4,1) + t588;
t585 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t584 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t583 = t107 * mrSges(4,1) - t106 * mrSges(4,2);
t582 = mrSges(4,2) - mrSges(6,3) - t463;
t581 = -mrSges(6,1) * t113 - mrSges(7,1) * t16 + mrSges(7,2) * t17;
t579 = t65 * mrSges(6,2) + 0.2e1 * Ifges(6,1) * t564 + 0.2e1 * Ifges(6,4) * t563 + 0.2e1 * Ifges(6,5) * t547;
t545 = -t268 / 0.2e1;
t553 = -t623 / 0.2e1;
t555 = -t149 / 0.2e1;
t557 = -t123 / 0.2e1;
t559 = -t122 / 0.2e1;
t578 = Ifges(7,5) * t557 - Ifges(6,2) * t553 - Ifges(6,6) * t545 + Ifges(7,6) * t559 + Ifges(7,3) * t555 + t581;
t9 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t73;
t577 = mrSges(6,1) * t65 + Ifges(7,5) * t571 + Ifges(7,6) * t570 + Ifges(7,3) * t565 + t9 / 0.2e1 + t584 + (-t547 - t202 / 0.2e1) * Ifges(6,6) + (-t563 - t76 / 0.2e1) * Ifges(6,2) + (-t564 - t75 / 0.2e1) * Ifges(6,4);
t554 = t149 / 0.2e1;
t556 = t123 / 0.2e1;
t558 = t122 / 0.2e1;
t576 = -Ifges(6,4) * t550 + Ifges(7,5) * t556 - Ifges(6,2) * t552 - Ifges(6,6) * t544 + Ifges(7,6) * t558 + Ifges(7,3) * t554 - t581 - t626;
t574 = Ifges(7,1) * t571 + Ifges(7,4) * t570 + Ifges(7,5) * t565;
t523 = Ifges(7,4) * t123;
t55 = t122 * Ifges(7,2) + t149 * Ifges(7,6) + t523;
t567 = t55 / 0.2e1;
t120 = Ifges(7,4) * t122;
t56 = t123 * Ifges(7,1) + t149 * Ifges(7,5) + t120;
t566 = -t56 / 0.2e1;
t560 = Ifges(5,1) * t548 + Ifges(5,4) * t549 + Ifges(5,5) * t546;
t551 = -t622 / 0.2e1;
t543 = -t275 / 0.2e1;
t542 = t275 / 0.2e1;
t539 = t278 / 0.2e1;
t537 = t370 / 0.2e1;
t536 = t375 / 0.2e1;
t528 = mrSges(6,3) * t623;
t526 = Ifges(5,4) * t366;
t525 = Ifges(5,4) * t370;
t522 = Ifges(7,4) * t372;
t521 = Ifges(7,4) * t375;
t520 = t170 * mrSges(4,3);
t519 = t171 * mrSges(4,3);
t518 = t278 * Ifges(4,4);
t510 = t623 * t372;
t509 = t623 * t375;
t501 = t289 * t372;
t500 = t294 * t366;
t499 = t335 * t372;
t498 = t335 * t375;
t489 = t377 * pkin(1) + qJ(2) * t473;
t479 = Ifges(6,5) * t75 + Ifges(6,6) * t76 + Ifges(6,3) * t202;
t477 = t56 * t536;
t476 = t615 - t616 + t617;
t465 = -t366 * (t216 * Ifges(5,4) + Ifges(5,2) * t453 + Ifges(5,6) * t275) / 0.2e1;
t464 = (t216 * Ifges(5,1) + Ifges(5,4) * t453 + Ifges(5,5) * t275) * t537;
t460 = -t483 / 0.2e1;
t444 = -pkin(1) * t534 + qJ(2) * t494;
t439 = qJDD(1) * t457;
t435 = mrSges(4,1) * t289 + mrSges(4,2) * t290;
t429 = Ifges(5,1) * t370 - t526;
t428 = Ifges(7,1) * t375 - t522;
t427 = -Ifges(5,2) * t366 + t525;
t426 = -Ifges(7,2) * t372 + t521;
t425 = Ifges(5,5) * t370 - Ifges(5,6) * t366;
t424 = Ifges(7,5) * t375 - Ifges(7,6) * t372;
t38 = pkin(11) * t289 + t529;
t176 = -t374 * t286 + t291 * t442 + t307 * t474;
t169 = -t319 * pkin(3) - t176;
t121 = -t231 * pkin(4) + t169;
t417 = t376 * t231 - t232 * t373;
t60 = -pkin(5) * t417 - t173 * pkin(11) + t121;
t20 = t372 * t60 + t375 * t38;
t19 = -t372 * t38 + t375 * t60;
t82 = -mrSges(7,2) * t149 + mrSges(7,3) * t122;
t83 = mrSges(7,1) * t149 - mrSges(7,3) * t123;
t421 = -t372 * t83 + t375 * t82;
t41 = -t373 * t86 + t376 * t79;
t134 = t173 * t375 + t501;
t14 = t373 * t74 + t376 * t85 + t79 * t484 - t485 * t86;
t410 = t28 * t430;
t407 = -(-qJ(2) * t469 + t352) * t367 + t315 * t511;
t323 = t377 * t511 - t459 * t534;
t237 = t323 * t374 + t535 * t594;
t404 = -g(1) * t237 - g(2) * t233 - g(3) * t289;
t402 = -mrSges(3,1) * t439 + mrSges(3,2) * t461;
t401 = mrSges(3,1) * t513 - mrSges(3,3) * t497;
t396 = -mrSges(3,2) * t513 + mrSges(3,3) * t457;
t383 = -t322 * pkin(2) - pkin(9) * t294 + t444;
t382 = t323 * pkin(2) + pkin(9) * t381 + t489;
t379 = t381 * t366;
t238 = t323 * t535 - t374 * t594;
t378 = pkin(4) * t379 + t237 * t530 + t238 * t361 + t382;
t156 = qJD(2) * t388 + (t266 + (t307 * t368 + t454) * t374) * qJD(3);
t126 = pkin(4) * t503 + t156;
t353 = -pkin(1) * t481 + qJDD(2);
t329 = t396 * qJD(1);
t328 = t401 * qJD(1);
t324 = -qJ(2) * t497 + t358;
t267 = Ifges(4,4) * t275;
t225 = t290 * t363 + t319 * t362;
t222 = -mrSges(4,2) * t312 - mrSges(4,3) * t275;
t204 = mrSges(4,1) * t275 + mrSges(4,2) * t278;
t191 = t238 * t363 + t362 * t381;
t190 = t238 * t362 - t363 * t381;
t185 = t278 * Ifges(4,1) - t267 + t612;
t184 = -t275 * Ifges(4,2) + t518 + t611;
t183 = -mrSges(4,2) * t311 - mrSges(4,3) * t206;
t182 = mrSges(4,1) * t311 - mrSges(4,3) * t205;
t138 = mrSges(4,1) * t206 + mrSges(4,2) * t205;
t136 = t191 * t375 + t237 * t372;
t135 = -t191 * t372 + t237 * t375;
t124 = -mrSges(6,2) * t268 + t528;
t119 = mrSges(5,1) * t206 - mrSges(5,3) * t179;
t118 = -mrSges(5,2) * t206 + mrSges(5,3) * t178;
t115 = qJD(5) * t173 + t276 * t335;
t114 = qJD(5) * t417 - t276 * t334;
t104 = pkin(5) * t622 - pkin(11) * t623;
t103 = -mrSges(6,1) * t623 + mrSges(6,2) * t622;
t94 = t179 * Ifges(5,4) + t178 * Ifges(5,2) + t206 * Ifges(5,6);
t62 = -qJD(6) * t134 - t114 * t372 + t277 * t375;
t61 = qJD(6) * t133 + t114 * t375 + t277 * t372;
t59 = -mrSges(6,2) * t202 + mrSges(6,3) * t76;
t49 = t115 * pkin(5) - t114 * pkin(11) + t126;
t37 = -pkin(5) * t289 - t41;
t25 = -mrSges(7,2) * t73 + mrSges(7,3) * t46;
t24 = mrSges(7,1) * t73 - mrSges(7,3) * t45;
t23 = t104 * t372 + t33 * t375;
t22 = t104 * t375 - t33 * t372;
t13 = -pkin(5) * t277 - t15;
t12 = pkin(11) * t277 + t14;
t10 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t73 * Ifges(7,6);
t6 = -pkin(5) * t202 - t8;
t4 = -qJD(6) * t20 - t12 * t372 + t375 * t49;
t3 = qJD(6) * t19 + t12 * t375 + t372 * t49;
t11 = [(Ifges(7,4) * t134 + Ifges(7,2) * t133) * t570 + (Ifges(3,5) * t461 + Ifges(3,6) * t439 + Ifges(3,3) * t452) * t513 + (mrSges(5,1) * t47 - mrSges(5,2) * t48 - mrSges(4,3) * t106 - Ifges(4,4) * t205 + Ifges(5,5) * t548 + Ifges(6,5) * t564 + Ifges(4,2) * t206 - Ifges(4,6) * t311 + Ifges(5,6) * t549 + Ifges(6,6) * t563 + Ifges(5,3) * t546 + Ifges(6,3) * t547 + t585) * t289 + (Ifges(5,5) * t232 + Ifges(5,6) * t231) * t546 + (-t520 + t453 * t427 / 0.2e1 + t216 * t429 / 0.2e1 + t464 + t465 + t612 / 0.2e1 + t207 * mrSges(4,2) + t185 / 0.2e1 + Ifges(4,1) * t539 + t425 * t542 + Ifges(4,4) * t543 + t601) * t276 - t328 * t468 + t196 * t435 + (t324 * t401 + t325 * t396 + Ifges(2,3)) * qJDD(1) + m(3) * (t300 * t324 + t301 * t325) + ((Ifges(3,1) * t367 + Ifges(3,4) * t511) * t461 + (Ifges(3,5) * t367 + Ifges(3,6) * t511) * t452 + (Ifges(3,4) * t367 + Ifges(3,2) * t511) * t439 + m(3) * (-pkin(1) * t353 + qJD(2) * t407) + t353 * (-mrSges(3,1) * t511 + mrSges(3,2) * t367) - pkin(1) * t402) * t369 + (t1 * t133 - t134 * t2 - t16 * t61 + t17 * t62) * mrSges(7,3) + (Ifges(5,4) * t232 + Ifges(5,2) * t231) * t549 + (-m(3) * t489 - t323 * mrSges(3,1) + mrSges(3,2) * t394 - mrSges(3,3) * t473 - t377 * mrSges(2,1) + mrSges(2,2) * t534 - m(5) * (t238 * pkin(3) + t382) - (t238 * t370 + t379) * mrSges(5,1) - (-t238 * t366 + t370 * t381) * mrSges(5,2) - m(4) * t382 - t238 * mrSges(4,1) - mrSges(4,3) * t381 - m(7) * (t191 * pkin(5) + t378) - t136 * mrSges(7,1) - t135 * mrSges(7,2) - m(6) * t378 - t191 * mrSges(6,1) + t587 * t190 + t582 * t237) * g(2) + (-m(3) * t444 + t322 * mrSges(3,1) - mrSges(3,2) * t399 - mrSges(3,3) * t494 - m(5) * (t236 * pkin(3) + t383) - (t236 * t370 - t500) * mrSges(5,1) - (-t236 * t366 - t294 * t370) * mrSges(5,2) + mrSges(2,1) * t534 + t377 * mrSges(2,2) - m(4) * t383 - t236 * mrSges(4,1) + t294 * mrSges(4,3) + t621 * (-pkin(4) * t500 - t233 * t530 + t236 * t361 + t383) + t587 * t639 - t609 * t640 - (-t430 + t582) * t233) * g(1) + t300 * t401 + (Ifges(7,4) * t61 + Ifges(7,2) * t62) * t558 + t204 * t448 + (t231 * t48 - t232 * t47 - t502 * t92 - t503 * t93) * mrSges(5,3) + m(5) * (t100 * t93 + t108 * t47 + t109 * t48 + t169 * t96 + t92 * t99) + (-m(4) * t170 + m(5) * t144 - t595) * t156 + (Ifges(5,1) * t232 + Ifges(5,4) * t231) * t548 + (t476 / 0.2e1 + t615 / 0.2e1 - t616 / 0.2e1 + t617 / 0.2e1 + t583) * t319 + (-mrSges(6,3) * t34 + t576) * t115 + (-t8 * mrSges(6,3) + t579) * t173 + (-t33 * mrSges(6,3) + t635) * t114 + m(4) * (t106 * t177 + t107 * t176 + t155 * t171 + t196 * t217 + t207 * t448) - (-mrSges(6,3) * t7 + t577) * t417 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t571 + (Ifges(7,1) * t61 + Ifges(7,4) * t62) * t556 + t301 * t396 + (Ifges(5,5) * t179 + Ifges(5,6) * t178 + Ifges(5,3) * t206 + t479) * t289 / 0.2e1 + m(6) * (t113 * t126 + t121 * t65 + t14 * t34 + t15 * t33 + t41 * t8 + t529 * t7) + t529 * t59 + m(7) * (t1 * t20 + t13 * t28 + t16 * t4 + t17 * t3 + t19 * t2 + t37 * t6) + (Ifges(7,5) * t134 + Ifges(7,6) * t133) * t565 + (Ifges(7,5) * t61 + Ifges(7,6) * t62) * t554 + t37 * t18 + t20 * t25 + t19 * t24 + (-t519 + t610 / 0.2e1 + Ifges(6,5) * t550 + Ifges(6,6) * t552 + t619 - t618 + Ifges(6,3) * t544 + t613 / 0.2e1 - t611 / 0.2e1 + t620 + t207 * mrSges(4,1) - t184 / 0.2e1 - Ifges(4,4) * t539 + Ifges(5,3) * t542 - Ifges(4,2) * t543 - t34 * mrSges(6,2) + t607 / 0.2e1) * t277 + t134 * t574 + t232 * t560 + t62 * t567 + t96 * (-mrSges(5,1) * t231 + mrSges(5,2) * t232) + t231 * t94 / 0.2e1 + t155 * t222 + t217 * t138 + t100 * t180 + t99 * t181 + t176 * t182 + t177 * t183 + t169 * t112 + t133 * t10 / 0.2e1 + t6 * (-mrSges(7,1) * t133 + mrSges(7,2) * t134) + t126 * t103 + t14 * t124 + t15 * t125 + t329 * t440 + t41 * t58 + t61 * t56 / 0.2e1 + t28 * (-mrSges(7,1) * t62 + mrSges(7,2) * t61) + t13 * t67 + t3 * t82 + t4 * t83 + t109 * t118 + t108 * t119 + t121 * t36 + (-mrSges(4,3) * t107 + Ifges(4,1) * t205 - Ifges(4,4) * t206 + Ifges(4,5) * t311) * t290; t512 * t138 - t204 * t449 + (t222 + t627) * t447 + t630 * t124 + t631 * t82 + t632 * t83 + (t170 * t303 - t171 * t304 - t207 * t449 + t196 * t512 + (t535 * t107 + t106 * t374 + (-t170 * t374 + t171 * t535) * qJD(3)) * t368 + t624) * m(4) + (t47 * t320 + t48 * t321 + (t144 * t486 + (qJD(3) * t592 - t96) * t535) * t368 - t144 * t303 - t252 * t92 - t253 * t93 + t624) * m(5) + (-t407 * t488 + t353 + t624) * m(3) - t400 * t25 + (t8 * t416 + t7 * t258 + (t113 * t486 - t535 * t65) * t368 - t113 * t303 + t630 * t34 + t624) * m(6) + (-t1 * t400 + t16 * t632 + t17 * t631 + t2 * t239 - t416 * t6 + t624) * m(7) + t402 - t329 * t441 + t321 * t118 + t320 * t119 - t304 * t222 - t252 * t181 - t253 * t180 + t258 * t59 + t239 * t24 + (-m(6) * t33 - t641) * (t376 * t252 - t253 * t373 + t258 * qJD(5) + (t366 * t471 + t370 * t472) * t487) + t642 * (t103 - t595) + t328 * t469 + (t182 + t633) * t474 + t183 * t495 - t634 * t416; (-Ifges(7,1) * t408 - Ifges(7,4) * t409) * t556 + (Ifges(7,1) * t161 + Ifges(7,4) * t160) * t557 - t275 * t520 + (Ifges(6,4) * t195 + Ifges(6,6) * t278) * t553 + t96 * t434 + t602 * t67 + t603 * t125 + t604 * t124 + t605 * t83 + t606 * t82 + (-Ifges(7,5) * t408 - Ifges(7,6) * t409) * t554 + (Ifges(7,5) * t161 + Ifges(7,6) * t160) * t555 + (-g(1) * t238 + g(2) * t236 - g(3) * t290 + t643 * t33 - t334 * t7 - t335 * t8 - t34 * t644) * mrSges(6,3) + (-t113 * t195 + t278 * t34) * mrSges(6,2) - t278 * t619 - t278 * t620 + (-t504 * t92 - t505 * t93 + t591) * mrSges(5,3) + (-pkin(3) * t96 + qJ(4) * t591 + qJD(4) * t592 - t116 * t92 - t117 * t93 - t144 * t171) * m(5) - t10 * t499 / 0.2e1 + (t621 * (-t233 * t361 - t236 * t530) - t589 * t236 + t586 * t233) * g(2) + t627 * qJD(4) + t595 * t171 - t597 * t55 / 0.2e1 + (mrSges(7,1) * t597 - mrSges(7,2) * t596) * t28 + (-t1 * t499 + t16 * t596 - t17 * t597 - t2 * t498) * mrSges(7,3) - (-Ifges(6,4) * t551 + t578 + t626) * t194 + t576 * t327 + t577 * t334 + (t424 * t565 + t426 * t570 + t428 * t571 + t430 * t6 + t460 * t56 + t579) * t335 + t476 + (Ifges(6,1) * t195 + Ifges(6,5) * t278) * t551 - (t477 + t635) * t326 + (-t113 * t132 + t299 * t7 + t33 * t603 + t34 * t604 - t361 * t65 + t415 * t8) * m(6) + (t1 * t209 + t16 * t605 + t17 * t606 + t2 * t208 + t28 * t602 - t415 * t6) * m(7) + (t621 * (-t237 * t361 + t238 * t530) + t589 * t238 + t586 * t237) * g(1) + (t435 + t621 * (-t289 * t361 + t290 * t530) + t593 * t290 + t588 * t289) * g(3) + t583 - t453 * (Ifges(5,6) * t278 - t275 * t427) / 0.2e1 + (t118 * t370 - t119 * t366) * qJ(4) + t275 * t601 + t278 * t618 - t216 * (Ifges(5,5) * t278 - t275 * t429) / 0.2e1 + (Ifges(7,4) * t161 + Ifges(7,2) * t160) * t559 + (-Ifges(7,4) * t408 - Ifges(7,2) * t409) * t558 + (-Ifges(4,2) * t278 + t185 - t267) * t542 + t275 * t464 + t275 * t465 + (Ifges(6,5) * t195 + Ifges(6,3) * t278) * t545 - t361 * t36 - t312 * (-Ifges(4,5) * t275 - Ifges(4,6) * t278) / 0.2e1 + t299 * t59 - t207 * (mrSges(4,1) * t278 - mrSges(4,2) * t275) + t498 * t574 + t366 * t560 + t161 * t566 + (Ifges(5,1) * t366 + t525) * t548 + (Ifges(5,2) * t370 + t526) * t549 + t184 * t539 + (Ifges(5,3) * t278 - t275 * t425) * t543 + (Ifges(5,5) * t366 + Ifges(5,6) * t370) * t546 + t94 * t537 + t278 * t519 - t170 * t222 + t208 * t24 + t209 * t25 - t195 * t89 / 0.2e1 - t117 * t180 - t116 * t181 - t132 * t103 - t634 * t415 - (-Ifges(4,1) * t275 - t518 + t607) * t278 / 0.2e1 - pkin(3) * t112; -t453 * t180 + t216 * t181 + t375 * t24 + t372 * t25 + t608 * t622 + t421 * qJD(6) + (-t124 - t421) * t623 + (t1 * t372 - t622 * t28 + t2 * t375 + t404 + t149 * (-t16 * t372 + t17 * t375)) * m(7) + (t33 * t622 - t34 * t623 + t404 + t65) * m(6) + (t216 * t92 - t453 * t93 + t404 + t96) * m(5) - t633; t585 + (t477 + t410) * qJD(6) + (-pkin(5) * t6 - t16 * t22 - t17 * t23) * m(7) + ((-t483 + t510) * t17 + (-t482 + t509) * t16 + t590) * mrSges(7,3) + (m(7) * ((-t16 * t375 - t17 * t372) * qJD(6) + t590) + t375 * t25 - t372 * t24 - t83 * t482 - t82 * t483) * pkin(11) + t479 + (t148 + t89) * t553 + (-t524 + t54) * t551 + (Ifges(6,1) * t551 + Ifges(6,5) * t545 + t424 * t555 + t426 * t559 + t428 * t557 - t410 - t614) * t623 + t578 * t622 + (-t587 * t640 - t609 * t639) * g(2) + (t122 * t426 + t123 * t428 + t149 * t424) * qJD(6) / 0.2e1 + t6 * t431 + t55 * t460 - pkin(5) * t18 + (t528 - t124) * t33 + (Ifges(7,2) * t375 + t522) * t570 + (Ifges(7,1) * t372 + t521) * t571 + t372 * t574 + (Ifges(7,5) * t372 + Ifges(7,6) * t375) * t565 + t509 * t566 + t510 * t567 + t88 * t550 + t10 * t536 + (t527 + t641) * t34 - t23 * t82 - t22 * t83 + (t190 * t609 + t191 * t587) * g(1) + (t587 * t225 - t609 * (-t290 * t362 + t319 * t363)) * g(3); -t28 * (mrSges(7,1) * t123 + mrSges(7,2) * t122) - t16 * t82 + t17 * t83 + (Ifges(7,1) * t122 - t523) * t557 + t55 * t556 + (Ifges(7,5) * t122 - Ifges(7,6) * t123) * t555 - g(1) * (mrSges(7,1) * t135 - mrSges(7,2) * t136) - g(2) * ((t233 * t375 + t372 * t640) * mrSges(7,1) + (-t233 * t372 + t375 * t640) * mrSges(7,2)) - g(3) * ((-t225 * t372 + t283) * mrSges(7,1) + (-t225 * t375 - t501) * mrSges(7,2)) + (t122 * t16 + t123 * t17) * mrSges(7,3) + t9 + (-Ifges(7,2) * t123 + t120 + t56) * t559 + t584;];
tau  = t11;
