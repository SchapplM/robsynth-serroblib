% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:58
% EndTime: 2019-03-09 21:43:15
% DurationCPUTime: 45.29s
% Computational Cost: add. (12139->857), mult. (32061->1106), div. (0->0), fcn. (23888->8), ass. (0->363)
t294 = cos(qJ(3));
t288 = cos(pkin(6));
t372 = qJD(1) * t288 + qJD(2);
t258 = t294 * t372;
t291 = sin(qJ(3));
t292 = sin(qJ(2));
t287 = sin(pkin(6));
t408 = qJD(1) * t287;
t394 = t292 * t408;
t215 = -t291 * t394 + t258;
t579 = t215 / 0.2e1;
t295 = cos(qJ(2));
t393 = t295 * t408;
t329 = -qJD(3) + t393;
t544 = -t329 / 0.2e1;
t576 = Ifges(4,6) * t544;
t578 = Ifges(4,2) * t579 + t576;
t216 = t291 * t372 + t294 * t394;
t577 = t216 / 0.2e1;
t575 = Ifges(4,4) * t577;
t456 = pkin(1) * t292;
t278 = t288 * t456;
t356 = pkin(3) * t291 - pkin(10) * t294;
t420 = t287 * t295;
t574 = (t278 + (pkin(8) + t356) * t420) * qJD(1) - t356 * qJD(3);
t277 = pkin(8) * t420;
t196 = qJD(2) * pkin(9) + (t277 + (pkin(9) + t456) * t288) * qJD(1);
t231 = (-pkin(2) * t295 - pkin(9) * t292 - pkin(1)) * t287;
t209 = qJD(1) * t231;
t139 = t294 * t196 + t291 * t209;
t121 = -pkin(10) * t329 + t139;
t290 = sin(qJ(4));
t293 = cos(qJ(4));
t421 = t287 * t292;
t275 = pkin(8) * t421;
t455 = pkin(1) * t295;
t229 = t275 + (-pkin(2) - t455) * t288;
t195 = -qJD(2) * pkin(2) + qJD(1) * t229;
t299 = -t215 * pkin(3) - t216 * pkin(10) + t195;
t41 = t121 * t290 - t293 * t299;
t42 = t293 * t121 + t290 * t299;
t491 = t41 * mrSges(5,1) + t42 * mrSges(5,2) + t575 + t578;
t540 = Ifges(5,1) + Ifges(7,3);
t538 = -Ifges(6,5) - Ifges(7,4);
t537 = -Ifges(7,5) - Ifges(5,5);
t536 = Ifges(7,2) + Ifges(6,3);
t573 = Ifges(6,4) + t537;
t572 = -Ifges(5,6) - t538;
t415 = t294 * t295;
t205 = (t290 * t292 + t293 * t415) * t408;
t404 = qJD(4) * t291;
t571 = t290 * t404 + t205;
t168 = t216 * t290 + t293 * t329;
t169 = t293 * t216 - t290 * t329;
t211 = qJD(4) - t215;
t75 = t169 * Ifges(5,5) - t168 * Ifges(5,6) + t211 * Ifges(5,3);
t79 = t211 * Ifges(7,1) + t168 * Ifges(7,4) + t169 * Ifges(7,5);
t80 = t211 * Ifges(6,1) - t169 * Ifges(6,4) + t168 * Ifges(6,5);
t497 = t80 + t79 + t75;
t570 = t497 / 0.2e1 - t491;
t402 = t291 * qJD(3);
t406 = qJD(2) * t295;
t309 = -t292 * t402 + t294 * t406;
t541 = qJD(2) / 0.2e1;
t380 = t287 * t541;
t358 = qJD(1) * t380;
t323 = t292 * t358;
t407 = qJD(2) * t292;
t400 = qJD(1) * qJD(2);
t384 = t295 * t400;
t365 = t287 * t384;
t178 = qJD(3) * t216 + t291 * t365;
t470 = -t178 / 0.2e1;
t383 = t292 * t400;
t366 = t287 * t383;
t403 = qJD(4) * t293;
t298 = qJD(3) * t215 + t294 * t365;
t510 = -qJD(4) * t329 + t298;
t99 = t216 * t403 + t290 * t510 - t293 * t366;
t481 = t99 / 0.2e1;
t482 = -t99 / 0.2e1;
t405 = qJD(4) * t290;
t98 = t216 * t405 - t290 * t366 - t293 * t510;
t483 = t98 / 0.2e1;
t484 = -t98 / 0.2e1;
t543 = t408 / 0.2e1;
t321 = (pkin(2) * t292 - pkin(9) * t295) * t287;
t237 = qJD(2) * t321;
t225 = qJD(1) * t237;
t511 = t288 * t455 - t275;
t239 = t511 * qJD(2);
t226 = qJD(1) * t239;
t401 = t294 * qJD(3);
t71 = -t196 * t402 + t209 * t401 + t291 * t225 + t294 * t226;
t559 = pkin(10) * t366 + qJD(4) * t299 + t71;
t90 = pkin(1) * t288 * t383 + t178 * pkin(3) + pkin(8) * t365 - pkin(10) * t298;
t10 = -t121 * t403 - t290 * t559 + t293 * t90;
t449 = pkin(4) + qJ(6);
t1 = -pkin(5) * t98 - qJD(6) * t211 - t178 * t449 - t10;
t9 = -t121 * t405 + t290 * t90 + t293 * t559;
t5 = -qJ(5) * t178 - qJD(5) * t211 - t9;
t2 = -pkin(5) * t99 - t5;
t7 = -pkin(4) * t178 - t10;
t547 = t10 * mrSges(5,1) - t9 * mrSges(5,2) + t7 * mrSges(6,2) + t2 * mrSges(7,2) - t5 * mrSges(6,3) - t1 * mrSges(7,3);
t316 = qJD(3) * t258;
t568 = t316 / 0.2e1;
t569 = t547 - 0.2e1 * Ifges(4,2) * t470 - Ifges(4,6) * t323 - t537 * t484 - t538 * t481 - Ifges(4,4) * t568 - (Ifges(4,4) * t309 + Ifges(4,6) * t407) * t543 + Ifges(6,4) * t483 + Ifges(5,6) * t482;
t539 = -Ifges(5,4) + Ifges(7,6);
t567 = Ifges(6,6) - Ifges(7,6);
t235 = t511 * qJD(1);
t236 = qJD(1) * t321;
t162 = t294 * t235 + t291 * t236;
t149 = pkin(10) * t394 + t162;
t266 = -pkin(3) * t294 - pkin(10) * t291 - pkin(2);
t566 = -t293 * t149 + t266 * t403 - t290 * t574;
t161 = -t291 * t235 + t236 * t294;
t148 = -pkin(3) * t394 - t161;
t284 = pkin(9) * t401;
t565 = -t148 + t284;
t501 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t322 = pkin(5) * t169 + t41;
t553 = qJD(5) + t322;
t17 = -t211 * t449 + t553;
t37 = -t211 * qJ(5) - t42;
t552 = t168 * pkin(5) - qJD(6);
t27 = -t37 - t552;
t516 = -qJD(5) - t41;
t36 = -pkin(4) * t211 - t516;
t433 = t139 * mrSges(4,3);
t561 = -t195 * mrSges(4,1) - t36 * mrSges(6,2) - t27 * mrSges(7,2) + t37 * mrSges(6,3) + t17 * mrSges(7,3) + t433 + t575;
t165 = Ifges(5,4) * t168;
t432 = t168 * Ifges(7,6);
t532 = t169 * t540 - t537 * t211 - t165 + t432;
t163 = Ifges(7,6) * t169;
t430 = t169 * Ifges(6,6);
t531 = t168 * t536 - t211 * t538 + t163 - t430;
t368 = t291 * t393;
t558 = -qJ(5) * t368 - qJD(5) * t294 + t566;
t416 = t293 * t294;
t281 = pkin(9) * t416;
t557 = qJD(4) * t281 - t290 * t149 + t266 * t405 + t293 * t574;
t554 = t290 * t401 + t291 * t403;
t556 = pkin(4) * t554 + qJ(5) * t571 + t565;
t422 = t215 * t290;
t555 = -qJD(5) * t290 - t139 + (t405 - t422) * pkin(4);
t23 = -t98 * Ifges(5,4) - t99 * Ifges(5,2) + t178 * Ifges(5,6);
t534 = -t178 * t538 + t536 * t99 + t567 * t98;
t550 = mrSges(6,1) * t5 - mrSges(7,1) * t2 - mrSges(5,3) * t9 - t23 / 0.2e1 + t534 / 0.2e1;
t486 = Ifges(6,4) * t470 + Ifges(6,2) * t484 + Ifges(6,6) * t482;
t533 = -t178 * t537 + t539 * t99 - t540 * t98;
t549 = mrSges(6,1) * t7 + mrSges(7,1) * t1 - mrSges(5,3) * t10 + t486 + t533 / 0.2e1;
t545 = t298 / 0.2e1;
t542 = pkin(9) * (t293 * t402 + t294 * t405);
t530 = t226 * mrSges(3,2);
t395 = -pkin(9) * t290 - pkin(4);
t529 = t368 * t449 + qJD(6) * t294 + (-qJ(6) + t395) * t402 + t557 + (qJD(3) * t416 - t571) * pkin(5);
t418 = t290 * t294;
t204 = -t293 * t394 + t393 * t418;
t280 = pkin(9) * t418;
t417 = t291 * t293;
t528 = t204 * pkin(5) + (-pkin(5) * t417 - t280) * qJD(4) + (-pkin(5) * t418 + (-pkin(9) * t293 + qJ(5)) * t291) * qJD(3) + t558;
t423 = qJ(5) * t293;
t325 = qJ(6) * t290 - t423;
t527 = -t204 * t449 + t325 * t401 + (qJD(6) * t290 + (qJ(6) * qJD(4) - qJD(5)) * t293) * t291 + t556;
t526 = -qJ(5) * t402 + t542 - t558;
t525 = pkin(4) * t368 + t395 * t402 + t557;
t524 = -pkin(4) * t204 + (-qJ(5) * t401 - qJD(5) * t291) * t293 + t556;
t523 = -t542 + t566;
t522 = -qJD(6) * t293 + t211 * t325 + t555;
t521 = -qJ(5) * t403 + t215 * t423 + t555;
t138 = -t291 * t196 + t294 * t209;
t157 = pkin(3) * t216 - pkin(10) * t215;
t67 = t293 * t138 + t290 * t157;
t47 = -qJ(5) * t216 - t67;
t480 = pkin(5) + pkin(10);
t520 = pkin(5) * t422 - t480 * t405 + t47;
t128 = t290 * t138;
t268 = t480 * t293;
t519 = qJD(4) * t268 - t128 - (pkin(5) * t215 - t157) * t293 + t449 * t216;
t398 = pkin(9) * t402;
t518 = -t290 * t398 + t557;
t517 = t292 * (Ifges(4,5) * t216 + Ifges(4,6) * t215 - Ifges(4,3) * t329);
t515 = -t42 + t552;
t242 = -t288 * t294 + t291 * t421;
t387 = t287 * t406;
t186 = -qJD(3) * t242 + t294 * t387;
t514 = -qJD(4) * t420 + t186;
t513 = t204 - t554;
t512 = -t293 * t401 + t571;
t436 = Ifges(7,6) * t293;
t438 = Ifges(6,6) * t293;
t509 = t290 * t536 + t436 - t438;
t437 = Ifges(7,6) * t290;
t439 = Ifges(6,6) * t290;
t508 = t293 * t536 - t437 + t439;
t507 = -t538 * t293 + (Ifges(6,4) - Ifges(7,5)) * t290;
t441 = Ifges(5,4) * t290;
t506 = t293 * t540 + t437 - t441;
t505 = t368 - t402;
t72 = -t196 * t401 - t209 * t402 + t225 * t294 - t291 * t226;
t504 = -t291 * t72 + t294 * t71;
t503 = -t10 * t290 + t293 * t9;
t502 = t290 * t7 - t293 * t5;
t498 = t178 * t501 + t572 * t99 + t573 * t98;
t496 = t290 * t572 - t293 * t573;
t371 = mrSges(3,3) * t394;
t493 = -m(4) * t195 + mrSges(3,1) * t372 + mrSges(4,1) * t215 - mrSges(4,2) * t216 - t371;
t120 = pkin(3) * t329 - t138;
t302 = -t169 * qJ(5) + t120;
t28 = t168 * t449 + t302;
t326 = t290 * t42 - t293 * t41;
t327 = t290 * t37 + t293 * t36;
t337 = -Ifges(6,2) * t293 + t439;
t440 = Ifges(5,4) * t293;
t346 = -Ifges(5,2) * t290 + t440;
t351 = -mrSges(7,2) * t293 + mrSges(7,3) * t290;
t352 = -mrSges(6,2) * t290 - mrSges(6,3) * t293;
t353 = mrSges(5,1) * t290 + mrSges(5,2) * t293;
t40 = t168 * pkin(4) + t302;
t164 = Ifges(6,6) * t168;
t77 = t211 * Ifges(6,4) - t169 * Ifges(6,2) + t164;
t426 = t293 * t77;
t431 = t169 * Ifges(5,4);
t78 = -Ifges(5,2) * t168 + Ifges(5,6) * t211 + t431;
t429 = t290 * t78;
t458 = t293 / 0.2e1;
t460 = t290 / 0.2e1;
t462 = t211 / 0.2e1;
t471 = t169 / 0.2e1;
t472 = -t169 / 0.2e1;
t473 = t168 / 0.2e1;
t474 = -t168 / 0.2e1;
t490 = t327 * mrSges(6,1) + (t17 * t293 - t27 * t290) * mrSges(7,1) - t326 * mrSges(5,3) + t120 * t353 + t28 * t351 + t40 * t352 + t346 * t474 + t337 * t472 - t429 / 0.2e1 - t426 / 0.2e1 + t509 * t473 + t506 * t471 + t531 * t460 + t532 * t458 + t496 * t462;
t488 = t287 ^ 2;
t479 = Ifges(4,1) * t568 + Ifges(4,4) * t470 + (Ifges(4,1) * t309 + Ifges(4,5) * t407) * t543;
t469 = t178 / 0.2e1;
t463 = -t211 / 0.2e1;
t285 = t291 * pkin(9);
t448 = mrSges(5,3) * t168;
t447 = mrSges(5,3) * t169;
t446 = Ifges(3,4) * t292;
t445 = Ifges(3,4) * t295;
t443 = Ifges(4,4) * t291;
t442 = Ifges(4,4) * t294;
t434 = t138 * mrSges(4,3);
t62 = -pkin(3) * t366 - t72;
t428 = t291 * t62;
t424 = qJ(5) * t168;
t419 = t290 * t291;
t123 = mrSges(6,1) * t168 - mrSges(6,3) * t211;
t126 = -mrSges(5,2) * t211 - t448;
t414 = t123 - t126;
t124 = -mrSges(7,1) * t168 + mrSges(7,2) * t211;
t413 = t124 - t123;
t125 = mrSges(6,1) * t169 + mrSges(6,2) * t211;
t127 = mrSges(5,1) * t211 - t447;
t412 = t125 - t127;
t243 = t288 * t291 + t294 * t421;
t144 = t242 * pkin(3) - t243 * pkin(10) + t229;
t410 = t277 + t278;
t230 = pkin(9) * t288 + t410;
t159 = t294 * t230 + t291 * t231;
t146 = -pkin(10) * t420 + t159;
t64 = t290 * t144 + t293 * t146;
t409 = pkin(4) * t419 + t285;
t224 = t290 * t266 + t281;
t397 = Ifges(4,5) * t298 - Ifges(4,6) * t178 + Ifges(4,3) * t366;
t392 = t287 * t407;
t379 = -t404 / 0.2e1;
t378 = t404 / 0.2e1;
t375 = t401 / 0.2e1;
t56 = -t98 * mrSges(6,1) + t178 * mrSges(6,2);
t55 = -t99 * mrSges(7,1) + t178 * mrSges(7,2);
t53 = -t98 * mrSges(7,1) - t178 * mrSges(7,3);
t374 = -qJ(5) * t290 - pkin(3);
t66 = t157 * t293 - t128;
t63 = t144 * t293 - t290 * t146;
t158 = -t291 * t230 + t231 * t294;
t223 = t266 * t293 - t280;
t370 = mrSges(3,3) * t393;
t49 = -qJ(5) * t242 - t64;
t199 = qJ(5) * t294 - t224;
t145 = pkin(3) * t420 - t158;
t350 = Ifges(4,1) * t294 - t443;
t348 = Ifges(5,1) * t290 + t440;
t347 = -Ifges(4,2) * t291 + t442;
t345 = Ifges(5,2) * t293 + t441;
t340 = Ifges(4,5) * t294 - Ifges(4,6) * t291;
t338 = Ifges(5,5) * t290 + Ifges(5,6) * t293;
t336 = Ifges(6,2) * t290 + t438;
t331 = -Ifges(7,3) * t290 + t436;
t185 = qJD(3) * t243 + t291 * t387;
t240 = t410 * qJD(2);
t110 = t185 * pkin(3) - t186 * pkin(10) + t240;
t100 = -t230 * t402 + t231 * t401 + t291 * t237 + t294 * t239;
t85 = pkin(10) * t392 + t100;
t15 = t110 * t293 - t144 * t405 - t146 * t403 - t290 * t85;
t101 = -t230 * t401 - t231 * t402 + t237 * t294 - t291 * t239;
t320 = Ifges(4,5) * t329;
t318 = t195 * (mrSges(4,1) * t291 + mrSges(4,2) * t294);
t14 = t290 * t110 + t144 * t403 - t146 * t405 + t293 * t85;
t313 = Ifges(3,5) * t365 - Ifges(3,6) * t366;
t312 = pkin(1) * t488 * (mrSges(3,1) * t292 + mrSges(3,2) * t295);
t188 = t243 * t293 - t290 * t420;
t311 = -qJ(5) * t188 + t145;
t310 = t292 * t488 * (Ifges(3,1) * t295 - t446);
t86 = -pkin(3) * t392 - t101;
t304 = t287 * t372 * (Ifges(3,5) * t295 - Ifges(3,6) * t292);
t303 = (Ifges(3,6) * t288 + (Ifges(3,2) * t295 + t446) * t287) * qJD(1);
t11 = -qJ(5) * t185 - qJD(5) * t242 - t14;
t116 = -t243 * t405 + t290 * t392 + t293 * t514;
t301 = -qJ(5) * t116 - qJD(5) * t188 + t86;
t300 = qJ(5) * t98 - qJD(5) * t169 + t62;
t286 = t294 * pkin(4);
t270 = Ifges(3,4) * t393;
t267 = t480 * t290;
t259 = -pkin(4) * t293 + t374;
t247 = -t293 * t449 + t374;
t238 = t410 * qJD(1);
t234 = -mrSges(3,2) * t372 + t370;
t232 = -qJ(5) * t417 + t409;
t227 = qJD(1) * t240;
t210 = Ifges(4,4) * t215;
t200 = -t223 + t286;
t197 = t291 * t325 + t409;
t193 = Ifges(3,1) * t394 + Ifges(3,5) * t372 + t270;
t192 = Ifges(3,6) * qJD(2) + t303;
t187 = t243 * t290 + t293 * t420;
t181 = -pkin(5) * t419 - t199;
t180 = -mrSges(4,1) * t329 - t216 * mrSges(4,3);
t179 = mrSges(4,2) * t329 + t215 * mrSges(4,3);
t177 = qJ(6) * t294 + t280 + t286 + (pkin(5) * t291 - t266) * t293;
t152 = -mrSges(4,2) * t366 - mrSges(4,3) * t178;
t151 = -mrSges(4,3) * t316 + (mrSges(4,1) * t407 - mrSges(4,3) * t309) * t408;
t134 = Ifges(4,1) * t216 + t210 - t320;
t122 = mrSges(7,1) * t169 - mrSges(7,3) * t211;
t115 = t243 * t403 + t290 * t514 - t293 * t392;
t113 = t178 * mrSges(4,1) + mrSges(4,2) * t298;
t108 = -mrSges(6,2) * t168 - mrSges(6,3) * t169;
t107 = mrSges(5,1) * t168 + mrSges(5,2) * t169;
t106 = pkin(4) * t169 + t424;
t105 = -mrSges(7,2) * t169 + mrSges(7,3) * t168;
t65 = pkin(4) * t187 + t311;
t58 = t169 * t449 + t424;
t54 = mrSges(6,1) * t99 - mrSges(6,3) * t178;
t52 = -mrSges(5,2) * t178 - mrSges(5,3) * t99;
t51 = mrSges(5,1) * t178 + mrSges(5,3) * t98;
t50 = -pkin(4) * t242 - t63;
t48 = -pkin(4) * t216 - t66;
t44 = t187 * t449 + t311;
t38 = -pkin(5) * t187 - t49;
t35 = pkin(5) * t188 - t242 * t449 - t63;
t33 = mrSges(7,2) * t98 + mrSges(7,3) * t99;
t32 = mrSges(5,1) * t99 - mrSges(5,2) * t98;
t31 = -mrSges(6,2) * t99 + mrSges(6,3) * t98;
t16 = pkin(4) * t115 + t301;
t13 = -pkin(4) * t185 - t15;
t12 = pkin(4) * t99 + t300;
t8 = qJD(6) * t187 + t115 * t449 + t301;
t6 = -pkin(5) * t115 - t11;
t4 = pkin(5) * t116 - qJD(6) * t242 - t185 * t449 - t15;
t3 = qJD(6) * t168 + t449 * t99 + t300;
t18 = [(t313 / 0.2e1 - t530) * t288 + (-t139 * t392 + t186 * t195 + t227 * t243 + t420 * t71) * mrSges(4,2) + (t498 / 0.2e1 - t71 * mrSges(4,3) - Ifges(4,4) * t545 + t227 * mrSges(4,1) + t501 * t469 + t569) * t242 + t380 * t517 + (-m(3) * t511 + m(4) * t229 - mrSges(3,1) * t288) * t227 + (t226 * t420 + t227 * t421 - t235 * t387 - t238 * t392 - t365 * t511 - t366 * t410) * mrSges(3,3) + (-m(3) * t235 - t493) * t240 + (Ifges(4,5) * t186 + Ifges(4,3) * t392) * t544 + (Ifges(4,5) * t243 - Ifges(4,3) * t420) * t323 + t304 * t541 + (t310 - 0.2e1 * t312) * t400 + (Ifges(4,1) * t243 - Ifges(4,5) * t420) * t545 + m(5) * (t10 * t63 + t120 * t86 + t14 * t42 + t145 * t62 - t15 * t41 + t64 * t9) + m(7) * (t1 * t35 + t17 * t4 + t2 * t38 + t27 * t6 + t28 * t8 + t3 * t44) + m(6) * (t11 * t37 + t12 * t65 + t13 * t36 + t16 * t40 + t49 * t5 + t50 * t7) + (t62 * mrSges(5,1) - t12 * mrSges(6,2) + t3 * mrSges(7,3) - Ifges(5,2) * t482 + Ifges(6,6) * t483 + t469 * t572 + t536 * t481 + t539 * t484 + t550) * t187 + (t572 * t462 + Ifges(6,6) * t472 - Ifges(5,2) * t474 + t37 * mrSges(6,1) - t27 * mrSges(7,1) - t42 * mrSges(5,3) - t40 * mrSges(6,2) + t120 * mrSges(5,1) - t78 / 0.2e1 + t28 * mrSges(7,3) + t531 / 0.2e1 + t536 * t473 + t539 * t471) * t115 + (t62 * mrSges(5,2) - t3 * mrSges(7,2) - t12 * mrSges(6,3) + Ifges(5,4) * t482 - Ifges(6,2) * t483 - t469 * t573 - t481 * t567 + t540 * t484 + t549) * t188 + (-t573 * t462 - Ifges(6,2) * t472 + Ifges(5,4) * t474 + t36 * mrSges(6,1) + t41 * mrSges(5,3) + t17 * mrSges(7,1) - t40 * mrSges(6,3) - t77 / 0.2e1 + t120 * mrSges(5,2) - t28 * mrSges(7,2) + t532 / 0.2e1 - t567 * t473 + t540 * t471) * t116 + m(4) * (t100 * t139 + t101 * t138 + t158 * t72 + t159 * t71) + (Ifges(4,4) * t243 - Ifges(4,6) * t420) * t470 + (t193 * t380 + (Ifges(3,5) * t288 + (t292 * Ifges(3,1) + t445) * t287) * t358) * t295 + t239 * t234 + t229 * t113 + t186 * t134 / 0.2e1 + t100 * t179 + t101 * t180 + t159 * t152 + t158 * t151 + t145 * t32 + t11 * t123 + t72 * (-mrSges(4,1) * t420 - t243 * mrSges(4,3)) + t6 * t124 + t13 * t125 + t14 * t126 + t15 * t127 - t397 * t420 / 0.2e1 + t4 * t122 + t8 * t105 + t86 * t107 + t16 * t108 - (t303 + t192) * t392 / 0.2e1 + m(3) * (t226 * t410 + t238 * t239) + t64 * t52 + t65 * t31 + t63 * t51 + t35 * t53 + t49 * t54 + t38 * t55 + t50 * t56 + t44 * t33 + t138 * (mrSges(4,1) * t392 - mrSges(4,3) * t186) + t243 * t479 + (Ifges(4,4) * t186 + Ifges(4,6) * t392) * t579 + (Ifges(6,4) * t472 + Ifges(5,6) * t474 + t501 * t462 - t537 * t471 - t538 * t473 - t561 + t570 - t578) * t185 + t488 * (-Ifges(3,2) * t292 + t445) * t384 + (Ifges(4,1) * t186 + Ifges(4,5) * t392) * t577; t524 * t108 + t525 * t125 + (t12 * t232 + t199 * t5 + t200 * t7 + t36 * t525 + t37 * t526 + t40 * t524) * m(6) + t526 * t123 + t527 * t105 + t528 * t124 + (t1 * t177 + t17 * t529 + t181 * t2 + t197 * t3 + t27 * t528 + t28 * t527) * m(7) + t529 * t122 + (pkin(9) * t152 - t569) * t294 + (-t433 + t570) * t402 + (mrSges(7,1) * t513 - mrSges(7,2) * t505) * t27 + (mrSges(6,2) * t513 + mrSges(6,3) * t512) * t40 + (-t41 * t512 + t42 * t513) * mrSges(5,3) + (-mrSges(6,1) * t513 + mrSges(6,3) * t505) * t37 + (mrSges(7,2) * t512 - mrSges(7,3) * t513) * t28 + (-t513 * mrSges(5,1) - t512 * mrSges(5,2)) * t120 + t532 * (t290 * t379 + t293 * t375 - t205 / 0.2e1) + t565 * t107 + (-t138 * t161 - t139 * t162 - pkin(2) * t227 + ((-t138 * t294 - t139 * t291) * qJD(3) + t504) * pkin(9)) * m(4) + t504 * mrSges(4,3) + (-t304 / 0.2e1 + (t312 - t310 / 0.2e1) * qJD(1)) * qJD(1) + t549 * t417 + t550 * t419 + t491 * t368 - t518 * t127 + (-t120 * t148 + t10 * t223 + t224 * t9 + (t120 * t401 + t428) * pkin(9) + t523 * t42 + t518 * t41) * m(5) + t523 * t126 + (-mrSges(4,1) * t294 + mrSges(4,2) * t291 - mrSges(3,1)) * t227 + (t371 + t493) * t238 + t531 * (t290 * t375 + t293 * t378 - t204 / 0.2e1) + (t139 * mrSges(4,2) * t292 - t138 * (mrSges(4,1) * t292 - mrSges(4,3) * t415)) * t408 + (t291 * t496 - t294 * t501) * t469 + (t329 * (Ifges(4,3) * t292 + t295 * t340) + t292 * t192) * t543 + t442 * t545 + (t204 / 0.2e1 + t293 * t379) * t78 + t353 * t428 + (t204 * t572 - t205 * t573 + t368 * t501) * t463 + (-t345 * t404 + (Ifges(5,6) * t291 + t294 * t346) * qJD(3) - t538 * t368 - t567 * t205 + t536 * t204) * t474 - t530 + ((-t338 + t507) * t404 + (t291 * t501 + t294 * t496) * qJD(3)) * t462 + (-mrSges(7,1) * t512 + mrSges(7,3) * t505) * t17 + (-mrSges(6,1) * t512 - mrSges(6,2) * t505) * t36 + (-t151 + t32) * t285 + t313 + (t205 / 0.2e1 + t290 * t378) * t77 + t134 * t375 - ((-Ifges(3,2) * t394 + t294 * t134 + t291 * t497 + t193 + t270) * t295 + t216 * (Ifges(4,5) * t292 + t295 * t350) + t215 * (Ifges(4,6) * t292 + t295 * t347) + t517) * t408 / 0.2e1 + t443 * t470 + (-t162 - t398) * t179 + (-t161 - t284) * t180 - t498 * t294 / 0.2e1 + t232 * t31 + t223 * t51 + t224 * t52 + t199 * t54 + t200 * t56 + t197 * t33 + t181 * t55 + t177 * t53 - t401 * t434 + (-t234 + t370) * t235 - pkin(2) * t113 - (t426 + t429) * t401 / 0.2e1 + (t215 * t347 + t216 * t350) * qJD(3) / 0.2e1 - t318 * t393 + (Ifges(6,4) * t368 - Ifges(6,2) * t205 + Ifges(6,6) * t204 + (-t348 + t331) * t404 + (-t291 * t537 + t294 * t506) * qJD(3)) * t471 + (Ifges(5,4) * t205 - Ifges(5,2) * t204 + Ifges(5,6) * t368 + t508 * t404 + (-t291 * t538 + t294 * t509) * qJD(3)) * t473 + (t336 * t404 + (Ifges(6,4) * t291 + t294 * t337) * qJD(3) - t537 * t368 + t540 * t205 + t539 * t204) * t472 + (t340 * t544 + t318) * qJD(3) + (Ifges(4,1) * t545 + Ifges(4,5) * t323 + t12 * t352 + t3 * t351 + t337 * t483 + t346 * t482 + t393 * t433 + t481 * t509 + t484 * t506 + t479) * t291; t502 * mrSges(6,1) + ((t52 - t54) * t293 + (-t51 + t56) * t290 + m(6) * t502 + m(5) * t503 + (-m(5) * t326 + m(6) * t327 + t290 * t414 + t293 * t412) * qJD(4)) * pkin(10) + t503 * mrSges(5,3) + t519 * t122 + t520 * t124 + (t12 * t259 - t36 * t48 - t37 * t47 + t40 * t521) * m(6) + t521 * t108 + (t1 * t267 + t17 * t519 + t2 * t268 + t247 * t3 + t27 * t520 + t28 * t522) * m(7) + t522 * t105 + (t576 - t79 / 0.2e1 - t80 / 0.2e1 - t75 / 0.2e1 + (-Ifges(7,1) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(5,3) / 0.2e1) * t211 + (-Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t169 + (Ifges(5,6) / 0.2e1 - Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t168 + t491 + t561) * t216 + (t180 - t107) * t139 + t397 + (-pkin(3) * t62 - t120 * t139 + t41 * t66 - t42 * t67) * m(5) + t338 * t469 + t507 * t470 + (t345 + t508) * t482 + (t348 + t336) * t484 + (t1 * t290 + t2 * t293) * mrSges(7,1) + t62 * (-mrSges(5,1) * t293 + mrSges(5,2) * t290) + t12 * (mrSges(6,2) * t293 - mrSges(6,3) * t290) + t3 * (-mrSges(7,2) * t290 - mrSges(7,3) * t293) + t259 * t31 + t267 * t53 + t268 * t55 + t247 * t33 - t138 * t179 - t48 * t125 - t67 * t126 - t66 * t127 - t47 * t123 - t71 * mrSges(4,2) + t72 * mrSges(4,1) - pkin(3) * t32 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t216 - t134 / 0.2e1 - t210 / 0.2e1 + t320 / 0.2e1 + t434 - t195 * mrSges(4,2) - t490) * t215 + t490 * qJD(4) + t331 * t483 + t290 * t486 + t533 * t460 - t534 * t293 / 0.2e1 + t23 * t458; t515 * t122 + (t169 * t536 + t164 - t432 + t77) * t474 + t498 + t547 + (-pkin(4) * t7 - qJ(5) * t5 - t106 * t40 - t36 * t42 + t37 * t516) * m(6) + (-t168 * t540 + t163 - t431 + t531) * t472 + (-t54 + t55) * qJ(5) + (t168 * t573 + t169 * t572) * t463 + (t168 * t17 + t169 * t27) * mrSges(7,1) + (-Ifges(5,2) * t169 - t165 + t532) * t473 - t449 * t53 + (-t414 + t448) * t41 + (-t412 + t447) * t42 - t120 * (mrSges(5,1) * t169 - mrSges(5,2) * t168) - t28 * (mrSges(7,2) * t168 + mrSges(7,3) * t169) - t40 * (-mrSges(6,2) * t169 + mrSges(6,3) * t168) + t322 * t124 - t58 * t105 - t106 * t108 + (qJ(5) * t2 - t1 * t449 + t515 * t17 + t27 * t553 - t28 * t58) * m(7) + t413 * qJD(5) + (t168 * t36 - t169 * t37) * mrSges(6,1) - pkin(4) * t56 + (Ifges(6,2) * t168 + t430 + t78) * t471; -t413 * t211 + (t105 + t108) * t169 + t53 + t56 + (t169 * t28 - t211 * t27 + t1) * m(7) + (t169 * t40 + t211 * t37 + t7) * m(6); -t168 * t105 + t211 * t122 + 0.2e1 * (t2 / 0.2e1 + t28 * t474 + t17 * t462) * m(7) + t55;];
tauc  = t18(:);
