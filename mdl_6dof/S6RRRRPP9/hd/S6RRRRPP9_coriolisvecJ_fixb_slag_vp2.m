% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:10:42
% EndTime: 2018-11-23 18:11:16
% DurationCPUTime: 35.31s
% Computational Cost: add. (12139->862), mult. (32061->1121), div. (0->0), fcn. (23888->8), ass. (0->370)
t287 = cos(pkin(6));
t291 = sin(qJ(2));
t462 = pkin(1) * t291;
t277 = t287 * t462;
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t460 = pkin(10) * t293;
t360 = pkin(3) * t290 - t460;
t286 = sin(pkin(6));
t294 = cos(qJ(2));
t425 = t286 * t294;
t579 = (t277 + (pkin(8) + t360) * t425) * qJD(1) - t360 * qJD(3);
t545 = Ifges(5,1) + Ifges(7,3);
t543 = -Ifges(6,5) - Ifges(7,4);
t542 = -Ifges(7,5) - Ifges(5,5);
t541 = Ifges(7,2) + Ifges(6,3);
t578 = Ifges(6,4) + t542;
t577 = -Ifges(5,6) - t543;
t289 = sin(qJ(4));
t292 = cos(qJ(4));
t413 = qJD(1) * t286;
t420 = t293 * t294;
t205 = (t289 * t291 + t292 * t420) * t413;
t409 = qJD(4) * t290;
t576 = t289 * t409 + t205;
t376 = qJD(1) * t287 + qJD(2);
t257 = t293 * t376;
t398 = t291 * t413;
t215 = -t290 * t398 + t257;
t211 = qJD(4) - t215;
t454 = pkin(4) + qJ(6);
t216 = t290 * t376 + t293 * t398;
t397 = t294 * t413;
t333 = -qJD(3) + t397;
t169 = t292 * t216 - t289 * t333;
t276 = pkin(8) * t425;
t196 = qJD(2) * pkin(9) + (t276 + (pkin(9) + t462) * t287) * qJD(1);
t230 = (-pkin(2) * t294 - pkin(9) * t291 - pkin(1)) * t286;
t209 = qJD(1) * t230;
t139 = t293 * t196 + t290 * t209;
t121 = -pkin(10) * t333 + t139;
t426 = t286 * t291;
t274 = pkin(8) * t426;
t461 = pkin(1) * t294;
t228 = t274 + (-pkin(2) - t461) * t287;
t195 = -qJD(2) * pkin(2) + qJD(1) * t228;
t298 = -t215 * pkin(3) - t216 * pkin(10) + t195;
t41 = t121 * t289 - t292 * t298;
t326 = pkin(5) * t169 + t41;
t557 = qJD(5) + t326;
t17 = -t211 * t454 + t557;
t42 = t292 * t121 + t289 * t298;
t37 = -t211 * qJ(5) - t42;
t168 = t216 * t289 + t292 * t333;
t556 = t168 * pkin(5) - qJD(6);
t27 = -t37 - t556;
t521 = -qJD(5) - t41;
t36 = -pkin(4) * t211 - t521;
t449 = Ifges(4,4) * t216;
t438 = t139 * mrSges(4,3);
t499 = t42 * mrSges(5,2) + t438;
t575 = t499 - t195 * mrSges(4,1) + t41 * mrSges(5,1) - t36 * mrSges(6,2) - t27 * mrSges(7,2) + t37 * mrSges(6,3) + t17 * mrSges(7,3) + t449 / 0.2e1;
t406 = t294 * qJD(2);
t412 = qJD(3) * t290;
t313 = -t291 * t412 + t293 * t406;
t407 = t291 * qJD(2);
t302 = (Ifges(4,4) * t313 + Ifges(4,6) * t407) * t286;
t546 = qJD(2) / 0.2e1;
t384 = t286 * t546;
t362 = qJD(1) * t384;
t327 = t291 * t362;
t405 = qJD(1) * qJD(2);
t387 = t294 * t405;
t369 = t286 * t387;
t178 = qJD(3) * t216 + t290 * t369;
t476 = -t178 / 0.2e1;
t388 = t291 * t405;
t370 = t286 * t388;
t408 = qJD(4) * t292;
t297 = qJD(3) * t215 + t293 * t369;
t516 = -qJD(4) * t333 + t297;
t99 = t216 * t408 + t289 * t516 - t292 * t370;
t487 = t99 / 0.2e1;
t488 = -t99 / 0.2e1;
t410 = qJD(4) * t289;
t98 = t216 * t410 - t289 * t370 - t292 * t516;
t489 = t98 / 0.2e1;
t490 = -t98 / 0.2e1;
t325 = (pkin(2) * t291 - pkin(9) * t294) * t286;
t236 = qJD(2) * t325;
t225 = qJD(1) * t236;
t517 = t287 * t461 - t274;
t238 = t517 * qJD(2);
t226 = qJD(1) * t238;
t411 = qJD(3) * t293;
t71 = -t196 * t412 + t209 * t411 + t290 * t225 + t293 * t226;
t563 = pkin(10) * t370 + qJD(4) * t298 + t71;
t90 = pkin(1) * t287 * t388 + t178 * pkin(3) + pkin(8) * t369 - pkin(10) * t297;
t10 = -t121 * t408 - t289 * t563 + t292 * t90;
t1 = -pkin(5) * t98 - qJD(6) * t211 - t178 * t454 - t10;
t9 = -t121 * t410 + t289 * t90 + t292 * t563;
t5 = -qJ(5) * t178 - qJD(5) * t211 - t9;
t2 = -pkin(5) * t99 - t5;
t7 = -pkin(4) * t178 - t10;
t552 = t10 * mrSges(5,1) - t9 * mrSges(5,2) + t7 * mrSges(6,2) + t2 * mrSges(7,2) - t5 * mrSges(6,3) - t1 * mrSges(7,3);
t572 = qJD(1) / 0.2e1;
t320 = qJD(3) * t257;
t573 = t320 / 0.2e1;
t574 = -Ifges(4,4) * t573 + Ifges(6,4) * t489 - 0.2e1 * Ifges(4,2) * t476 - Ifges(4,6) * t327 + Ifges(5,6) * t488 - t302 * t572 - t543 * t487 - t542 * t490 + t552;
t544 = -Ifges(5,4) + Ifges(7,6);
t540 = Ifges(6,6) - Ifges(7,6);
t234 = t517 * qJD(1);
t235 = qJD(1) * t325;
t162 = t293 * t234 + t290 * t235;
t149 = pkin(10) * t398 + t162;
t265 = -pkin(3) * t293 - pkin(10) * t290 - pkin(2);
t571 = -t292 * t149 + t265 * t408 - t289 * t579;
t161 = -t290 * t234 + t235 * t293;
t148 = -pkin(3) * t398 - t161;
t283 = pkin(9) * t411;
t570 = -t148 + t283;
t507 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t163 = Ifges(7,6) * t169;
t435 = t169 * Ifges(6,6);
t537 = t168 * t541 - t211 * t543 + t163 - t435;
t165 = Ifges(5,4) * t168;
t437 = t168 * Ifges(7,6);
t536 = t169 * t545 - t542 * t211 - t165 + t437;
t564 = Ifges(4,2) * t215;
t372 = t290 * t397;
t562 = -qJ(5) * t372 - qJD(5) * t293 + t571;
t421 = t292 * t293;
t280 = pkin(9) * t421;
t561 = qJD(4) * t280 - t289 * t149 + t265 * t410 + t292 * t579;
t558 = t289 * t411 + t290 * t408;
t560 = pkin(4) * t558 + qJ(5) * t576 + t570;
t427 = t215 * t289;
t559 = -qJD(5) * t289 - t139 + (t410 - t427) * pkin(4);
t323 = Ifges(4,6) * t333;
t133 = -t323 + t449 + t564;
t75 = t169 * Ifges(5,5) - t168 * Ifges(5,6) + t211 * Ifges(5,3);
t79 = t211 * Ifges(7,1) + t168 * Ifges(7,4) + t169 * Ifges(7,5);
t80 = t211 * Ifges(6,1) - t169 * Ifges(6,4) + t168 * Ifges(6,5);
t503 = t80 + t79 + t75;
t555 = -t133 / 0.2e1 + t503 / 0.2e1;
t23 = -Ifges(5,4) * t98 - Ifges(5,2) * t99 + Ifges(5,6) * t178;
t539 = -t178 * t543 + t540 * t98 + t541 * t99;
t554 = t539 / 0.2e1 + mrSges(6,1) * t5 - mrSges(7,1) * t2 - mrSges(5,3) * t9 - t23 / 0.2e1;
t493 = Ifges(6,4) * t476 + Ifges(6,2) * t490 + Ifges(6,6) * t488;
t538 = -t178 * t542 + t544 * t99 - t545 * t98;
t553 = t538 / 0.2e1 + mrSges(6,1) * t7 + mrSges(7,1) * t1 - mrSges(5,3) * t10 + t493;
t549 = t297 / 0.2e1;
t548 = -t333 / 0.2e1;
t547 = pkin(9) * (t292 * t412 + t293 * t410);
t535 = t226 * mrSges(3,2);
t399 = -pkin(9) * t289 - pkin(4);
t534 = t372 * t454 + qJD(6) * t293 + (-qJ(6) + t399) * t412 + t561 + (qJD(3) * t421 - t576) * pkin(5);
t423 = t289 * t293;
t204 = -t292 * t398 + t397 * t423;
t279 = pkin(9) * t423;
t422 = t290 * t292;
t533 = t204 * pkin(5) + (-pkin(5) * t422 - t279) * qJD(4) + (-pkin(5) * t423 + (-pkin(9) * t292 + qJ(5)) * t290) * qJD(3) + t562;
t428 = qJ(5) * t292;
t329 = qJ(6) * t289 - t428;
t532 = -t204 * t454 + t329 * t411 + (qJD(6) * t289 + (qJ(6) * qJD(4) - qJD(5)) * t292) * t290 + t560;
t531 = -qJ(5) * t412 + t547 - t562;
t530 = pkin(4) * t372 + t399 * t412 + t561;
t529 = -pkin(4) * t204 + (-qJ(5) * t411 - qJD(5) * t290) * t292 + t560;
t528 = -t547 + t571;
t527 = -qJD(6) * t292 + t211 * t329 + t559;
t526 = -qJ(5) * t408 + t215 * t428 + t559;
t138 = -t290 * t196 + t293 * t209;
t157 = pkin(3) * t216 - pkin(10) * t215;
t67 = t292 * t138 + t289 * t157;
t47 = -qJ(5) * t216 - t67;
t486 = pkin(5) + pkin(10);
t525 = pkin(5) * t427 - t486 * t410 + t47;
t128 = t289 * t138;
t267 = t486 * t292;
t524 = qJD(4) * t267 - t128 - (pkin(5) * t215 - t157) * t292 + t454 * t216;
t403 = pkin(9) * t412;
t523 = -t289 * t403 + t561;
t522 = t291 * (Ifges(4,5) * t216 + Ifges(4,6) * t215 - Ifges(4,3) * t333);
t520 = -t42 + t556;
t519 = t204 - t558;
t518 = -t292 * t411 + t576;
t442 = Ifges(7,6) * t289;
t446 = Ifges(5,4) * t289;
t515 = t292 * t545 + t442 - t446;
t441 = Ifges(7,6) * t292;
t443 = Ifges(6,6) * t292;
t514 = t289 * t541 + t441 - t443;
t444 = Ifges(6,6) * t289;
t513 = t292 * t541 - t442 + t444;
t512 = -t543 * t292 + (Ifges(6,4) - Ifges(7,5)) * t289;
t511 = t372 - t412;
t72 = -t196 * t411 - t209 * t412 + t225 * t293 - t290 * t226;
t510 = -t290 * t72 + t293 * t71;
t509 = -t10 * t289 + t292 * t9;
t508 = t289 * t7 - t292 * t5;
t504 = t178 * t507 + t577 * t99 + t578 * t98;
t502 = t289 * t577 - t292 * t578;
t375 = mrSges(3,3) * t398;
t498 = -m(4) * t195 + mrSges(3,1) * t376 + mrSges(4,1) * t215 - mrSges(4,2) * t216 - t375;
t120 = pkin(3) * t333 - t138;
t305 = -t169 * qJ(5) + t120;
t28 = t168 * t454 + t305;
t330 = t289 * t42 - t292 * t41;
t331 = t289 * t37 + t292 * t36;
t341 = -Ifges(6,2) * t292 + t444;
t445 = Ifges(5,4) * t292;
t350 = -Ifges(5,2) * t289 + t445;
t355 = -mrSges(7,2) * t292 + mrSges(7,3) * t289;
t356 = -mrSges(6,2) * t289 - mrSges(6,3) * t292;
t357 = mrSges(5,1) * t289 + mrSges(5,2) * t292;
t40 = t168 * pkin(4) + t305;
t164 = Ifges(6,6) * t168;
t77 = t211 * Ifges(6,4) - t169 * Ifges(6,2) + t164;
t431 = t292 * t77;
t436 = t169 * Ifges(5,4);
t78 = -Ifges(5,2) * t168 + Ifges(5,6) * t211 + t436;
t434 = t289 * t78;
t464 = t292 / 0.2e1;
t466 = t289 / 0.2e1;
t468 = t211 / 0.2e1;
t477 = t169 / 0.2e1;
t478 = -t169 / 0.2e1;
t479 = t168 / 0.2e1;
t480 = -t168 / 0.2e1;
t497 = t331 * mrSges(6,1) + (t17 * t292 - t27 * t289) * mrSges(7,1) - t330 * mrSges(5,3) - t431 / 0.2e1 - t434 / 0.2e1 + t120 * t357 + t40 * t356 + t28 * t355 + t341 * t478 + t350 * t480 + t514 * t479 + t515 * t477 + t537 * t466 + t536 * t464 + t502 * t468;
t495 = t286 ^ 2;
t491 = t77 / 0.2e1;
t303 = (Ifges(4,1) * t313 + Ifges(4,5) * t407) * t286;
t485 = Ifges(4,1) * t573 + Ifges(4,4) * t476 + t303 * t572;
t481 = t133 / 0.2e1;
t475 = t178 / 0.2e1;
t469 = -t211 / 0.2e1;
t284 = t290 * pkin(9);
t453 = mrSges(5,3) * t168;
t452 = mrSges(5,3) * t169;
t451 = Ifges(3,4) * t291;
t450 = Ifges(3,4) * t294;
t448 = Ifges(4,4) * t290;
t447 = Ifges(4,4) * t293;
t439 = t138 * mrSges(4,3);
t62 = -pkin(3) * t370 - t72;
t433 = t290 * t62;
t429 = qJ(5) * t168;
t424 = t289 * t290;
t123 = mrSges(6,1) * t168 - mrSges(6,3) * t211;
t124 = -mrSges(7,1) * t168 + mrSges(7,2) * t211;
t419 = t123 - t124;
t126 = -mrSges(5,2) * t211 - t453;
t418 = t123 - t126;
t125 = mrSges(6,1) * t169 + mrSges(6,2) * t211;
t127 = mrSges(5,1) * t211 - t452;
t417 = t125 - t127;
t241 = -t287 * t293 + t290 * t426;
t242 = t287 * t290 + t293 * t426;
t144 = t241 * pkin(3) - t242 * pkin(10) + t228;
t415 = t276 + t277;
t229 = pkin(9) * t287 + t415;
t159 = t293 * t229 + t290 * t230;
t146 = -pkin(10) * t425 + t159;
t64 = t289 * t144 + t292 * t146;
t414 = pkin(4) * t424 + t284;
t224 = t289 * t265 + t280;
t402 = t289 * t425;
t401 = Ifges(4,5) * t297 - Ifges(4,6) * t178 + Ifges(4,3) * t370;
t393 = t286 * t407;
t392 = t286 * t406;
t391 = t287 * t411;
t381 = t411 / 0.2e1;
t380 = -t409 / 0.2e1;
t379 = t409 / 0.2e1;
t56 = -t98 * mrSges(6,1) + t178 * mrSges(6,2);
t55 = -t99 * mrSges(7,1) + t178 * mrSges(7,2);
t53 = -t98 * mrSges(7,1) - t178 * mrSges(7,3);
t378 = -qJ(5) * t289 - pkin(3);
t66 = t157 * t292 - t128;
t63 = t144 * t292 - t289 * t146;
t158 = -t290 * t229 + t230 * t293;
t223 = t265 * t292 - t279;
t374 = mrSges(3,3) * t397;
t49 = -qJ(5) * t241 - t64;
t199 = qJ(5) * t293 - t224;
t145 = pkin(3) * t425 - t158;
t354 = Ifges(4,1) * t293 - t448;
t352 = Ifges(5,1) * t289 + t445;
t351 = -Ifges(4,2) * t290 + t447;
t349 = Ifges(5,2) * t292 + t446;
t344 = Ifges(4,5) * t293 - Ifges(4,6) * t290;
t342 = Ifges(5,5) * t289 + Ifges(5,6) * t292;
t340 = Ifges(6,2) * t289 + t443;
t335 = -Ifges(7,3) * t289 + t441;
t186 = qJD(3) * t242 + t290 * t392;
t309 = t241 * qJD(3);
t110 = t186 * pkin(3) + pkin(10) * t309 + (t277 + (pkin(8) - t460) * t425) * qJD(2);
t100 = -t229 * t412 + t230 * t411 + t290 * t236 + t293 * t238;
t85 = pkin(10) * t393 + t100;
t15 = t110 * t292 - t144 * t410 - t146 * t408 - t289 * t85;
t101 = -t229 * t411 - t230 * t412 + t236 * t293 - t290 * t238;
t187 = t289 * t242 + t292 * t425;
t324 = Ifges(4,5) * t333;
t322 = t195 * (mrSges(4,1) * t290 + mrSges(4,2) * t293);
t14 = t289 * t110 + t144 * t408 - t146 * t410 + t292 * t85;
t317 = Ifges(3,5) * t369 - Ifges(3,6) * t370;
t316 = pkin(1) * t495 * (mrSges(3,1) * t291 + mrSges(3,2) * t294);
t188 = t242 * t292 - t402;
t315 = -qJ(5) * t188 + t145;
t314 = t291 * t495 * (Ifges(3,1) * t294 - t451);
t239 = t415 * qJD(2);
t86 = -pkin(3) * t393 - t101;
t307 = t286 * t376 * (Ifges(3,5) * t294 - Ifges(3,6) * t291);
t306 = (Ifges(3,6) * t287 + (t294 * Ifges(3,2) + t451) * t286) * qJD(1);
t11 = -qJ(5) * t186 - qJD(5) * t241 - t14;
t304 = t293 * t392 - t309;
t301 = (mrSges(4,1) * t407 - mrSges(4,3) * t313) * t286;
t115 = qJD(4) * t187 - t289 * t393 - t292 * t304;
t300 = qJ(5) * t115 - qJD(5) * t188 + t86;
t299 = qJ(5) * t98 - qJD(5) * t169 + t62;
t285 = t293 * pkin(4);
t269 = Ifges(3,4) * t397;
t266 = t486 * t289;
t258 = -pkin(4) * t292 + t378;
t246 = -t292 * t454 + t378;
t237 = t415 * qJD(1);
t233 = -mrSges(3,2) * t376 + t374;
t231 = -qJ(5) * t422 + t414;
t227 = qJD(1) * t239;
t210 = Ifges(4,4) * t215;
t200 = -t223 + t285;
t197 = t290 * t329 + t414;
t193 = Ifges(3,1) * t398 + Ifges(3,5) * t376 + t269;
t192 = Ifges(3,6) * qJD(2) + t306;
t182 = -pkin(5) * t424 - t199;
t180 = -mrSges(4,1) * t333 - t216 * mrSges(4,3);
t179 = mrSges(4,2) * t333 + t215 * mrSges(4,3);
t177 = qJ(6) * t293 + t279 + t285 + (pkin(5) * t290 - t265) * t292;
t152 = -mrSges(4,2) * t370 - mrSges(4,3) * t178;
t151 = -mrSges(4,3) * t320 + qJD(1) * t301;
t134 = Ifges(4,1) * t216 + t210 - t324;
t122 = mrSges(7,1) * t169 - mrSges(7,3) * t211;
t116 = -qJD(4) * t402 + t242 * t408 + t289 * t304 - t292 * t393;
t113 = t178 * mrSges(4,1) + mrSges(4,2) * t297;
t108 = -mrSges(6,2) * t168 - mrSges(6,3) * t169;
t107 = mrSges(5,1) * t168 + mrSges(5,2) * t169;
t106 = pkin(4) * t169 + t429;
t105 = -mrSges(7,2) * t169 + mrSges(7,3) * t168;
t65 = pkin(4) * t187 + t315;
t58 = t169 * t454 + t429;
t54 = mrSges(6,1) * t99 - mrSges(6,3) * t178;
t52 = -mrSges(5,2) * t178 - mrSges(5,3) * t99;
t51 = mrSges(5,1) * t178 + mrSges(5,3) * t98;
t50 = -pkin(4) * t241 - t63;
t48 = -pkin(4) * t216 - t66;
t44 = t187 * t454 + t315;
t38 = -pkin(5) * t187 - t49;
t35 = pkin(5) * t188 - t241 * t454 - t63;
t33 = mrSges(7,2) * t98 + mrSges(7,3) * t99;
t32 = mrSges(5,1) * t99 - mrSges(5,2) * t98;
t31 = -mrSges(6,2) * t99 + mrSges(6,3) * t98;
t16 = pkin(4) * t116 + t300;
t13 = -pkin(4) * t186 - t15;
t12 = pkin(4) * t99 + t299;
t8 = qJD(6) * t187 + t116 * t454 + t300;
t6 = -pkin(5) * t116 - t11;
t4 = -pkin(5) * t115 - qJD(6) * t241 - t186 * t454 - t15;
t3 = qJD(6) * t168 + t454 * t99 + t299;
t18 = [t384 * t522 + (-m(3) * t234 - t498) * t239 + (mrSges(4,1) * t227 + t507 * t475 + t504 / 0.2e1 - Ifges(4,4) * t549 - t71 * mrSges(4,3) + t574) * t241 + (-t564 / 0.2e1 + Ifges(6,4) * t478 + t507 * t468 - t543 * t479 - Ifges(4,6) * t548 - t542 * t477 + Ifges(5,6) * t480 + t555 - t575) * t186 + t307 * t546 + (-m(3) * t517 + m(4) * t228 - mrSges(3,1) * t287) * t227 + (t226 * t425 + t227 * t426 - t234 * t392 - t237 * t393 - t369 * t517 - t370 * t415) * mrSges(3,3) + (t317 / 0.2e1 - t535) * t287 + (-t139 * t393 + t195 * t304 + t227 * t242 + t425 * t71) * mrSges(4,2) + (t62 * mrSges(5,1) - t12 * mrSges(6,2) + t3 * mrSges(7,3) - Ifges(5,2) * t488 + Ifges(6,6) * t489 + t475 * t577 + t541 * t487 + t544 * t490 + t554) * t187 + (t541 * t479 + t544 * t477 + t537 / 0.2e1 - Ifges(5,2) * t480 + Ifges(6,6) * t478 + t577 * t468 + mrSges(6,1) * t37 - mrSges(7,1) * t27 - mrSges(5,3) * t42 + mrSges(5,1) * t120 + mrSges(7,3) * t28 - mrSges(6,2) * t40 - t78 / 0.2e1) * t116 + (t62 * mrSges(5,2) - t3 * mrSges(7,2) - t12 * mrSges(6,3) + Ifges(5,4) * t488 - Ifges(6,2) * t489 - t475 * t578 - t540 * t487 + t545 * t490 + t553) * t188 + (t540 * t479 - t545 * t477 - t536 / 0.2e1 - Ifges(5,4) * t480 + t491 + Ifges(6,2) * t478 + t578 * t468 - mrSges(5,3) * t41 - mrSges(7,1) * t17 - mrSges(6,1) * t36 - mrSges(5,2) * t120 + mrSges(7,2) * t28 + mrSges(6,3) * t40) * t115 + (Ifges(4,1) * t242 - Ifges(4,5) * t425) * t549 + t216 * (Ifges(4,1) * t391 + t303) / 0.2e1 + (Ifges(4,5) * t391 + (Ifges(4,5) * t313 + Ifges(4,3) * t407) * t286) * t548 + (Ifges(4,5) * t242 - Ifges(4,3) * t425) * t327 + (t193 * t384 + (Ifges(3,5) * t287 + (t291 * Ifges(3,1) + t450) * t286) * t362) * t294 - (t306 + t192) * t393 / 0.2e1 + m(3) * (t226 * t415 + t237 * t238) + t304 * t134 / 0.2e1 + m(4) * (t100 * t139 + t101 * t138 + t158 * t72 + t159 * t71) + t242 * t485 + m(5) * (t10 * t63 + t120 * t86 + t14 * t42 + t145 * t62 - t15 * t41 + t64 * t9) + m(6) * (t11 * t37 + t12 * t65 + t13 * t36 + t16 * t40 + t49 * t5 + t50 * t7) + m(7) * (t1 * t35 + t17 * t4 + t2 * t38 + t27 * t6 + t28 * t8 + t3 * t44) + (t314 - 0.2e1 * t316) * t405 + t238 * t233 + t228 * t113 + t101 * t180 + t100 * t179 + t72 * (-mrSges(4,1) * t425 - t242 * mrSges(4,3)) - t401 * t425 / 0.2e1 + t159 * t152 + t158 * t151 + t145 * t32 + t6 * t124 + t13 * t125 + t14 * t126 + t15 * t127 + t4 * t122 + t11 * t123 + t86 * t107 + t16 * t108 + t8 * t105 + t64 * t52 + t65 * t31 + t63 * t51 + t35 * t53 + t49 * t54 + t38 * t55 + t50 * t56 + t44 * t33 + t138 * (-mrSges(4,3) * t391 + t301) + t495 * (-Ifges(3,2) * t291 + t450) * t387 + (Ifges(4,4) * t242 - Ifges(4,6) * t425) * t476 + t215 * (Ifges(4,4) * t391 + t302) / 0.2e1; (Ifges(6,4) * t372 - Ifges(6,2) * t205 + Ifges(6,6) * t204 + (-t352 + t335) * t409 + (-t290 * t542 + t293 * t515) * qJD(3)) * t477 + (-t349 * t409 + (Ifges(5,6) * t290 + t293 * t350) * qJD(3) - t543 * t372 - t540 * t205 + t541 * t204) * t480 + (Ifges(5,4) * t205 - Ifges(5,2) * t204 + Ifges(5,6) * t372 + t513 * t409 + (-t290 * t543 + t293 * t514) * qJD(3)) * t479 + (t375 + t498) * t237 + (pkin(9) * t152 - t574) * t293 + t447 * t549 + (-t499 + t555) * t412 + (-mrSges(4,1) * t293 + mrSges(4,2) * t290 - mrSges(3,1)) * t227 + (-t120 * t148 + t10 * t223 + t224 * t9 + (t120 * t411 + t433) * pkin(9) + t528 * t42 + t523 * t41) * m(5) + t528 * t126 + t529 * t108 + t530 * t125 + (t12 * t231 + t199 * t5 + t200 * t7 + t36 * t530 + t37 * t531 + t40 * t529) * m(6) + t531 * t123 + t532 * t105 + t533 * t124 + (t1 * t177 + t17 * t534 + t182 * t2 + t197 * t3 + t27 * t533 + t28 * t532) * m(7) + t534 * t122 - t523 * t127 + (t204 * t577 - t205 * t578 + t372 * t507) * t469 - (t215 * (Ifges(4,6) * t291 + t294 * t351) + t522 + t216 * (Ifges(4,5) * t291 + t294 * t354) + (-Ifges(3,2) * t398 + t293 * t134 + t290 * t503 + t193 + t269) * t294) * t413 / 0.2e1 + t536 * (t289 * t380 + t292 * t381 - t205 / 0.2e1) + t537 * (t289 * t381 + t292 * t379 - t204 / 0.2e1) + t289 * t77 * t379 + t553 * t422 + t554 * t424 + (-mrSges(7,1) * t518 + mrSges(7,3) * t511) * t17 + (mrSges(5,1) * t511 - mrSges(5,3) * t518) * t41 + (-mrSges(6,1) * t518 - mrSges(6,2) * t511) * t36 + (mrSges(6,2) * t519 + mrSges(6,3) * t518) * t40 + (mrSges(5,2) * t372 + mrSges(5,3) * t519) * t42 + (mrSges(7,1) * t519 - mrSges(7,2) * t511) * t27 + (mrSges(7,2) * t518 - mrSges(7,3) * t519) * t28 + (-t519 * mrSges(5,1) - t518 * mrSges(5,2)) * t120 + (-mrSges(6,1) * t519 + mrSges(6,3) * t511) * t37 + (-t161 - t283) * t180 + (t340 * t409 + (Ifges(6,4) * t290 + t293 * t341) * qJD(3) - t542 * t372 + t545 * t205 + t544 * t204) * t478 - (t431 + t434) * t411 / 0.2e1 + (t333 * (Ifges(4,3) * t291 + t294 * t344) + t291 * t192) * t413 / 0.2e1 + (t215 * t351 + t216 * t354) * qJD(3) / 0.2e1 + ((-t342 + t512) * t409 + (t290 * t507 + t293 * t502) * qJD(3)) * t468 + t372 * t481 + t205 * t491 + (-t307 / 0.2e1 + (t316 - t314 / 0.2e1) * qJD(1)) * qJD(1) + (-t233 + t374) * t234 + t448 * t476 + (t344 * t548 + t322) * qJD(3) + (Ifges(4,1) * t549 + Ifges(4,5) * t327 + t12 * t356 + t3 * t355 + t341 * t489 + t350 * t488 + t397 * t438 + t487 * t514 + t490 * t515 + t485) * t290 - t504 * t293 / 0.2e1 + (t290 * t502 - t293 * t507) * t475 + t317 + t134 * t381 + (-t138 * t161 - t139 * t162 - pkin(2) * t227 + ((-t138 * t293 - t139 * t290) * qJD(3) + t510) * pkin(9)) * m(4) + t510 * mrSges(4,3) + (-t162 - t403) * t179 + t223 * t51 + t224 * t52 + t231 * t31 + t199 * t54 + t200 * t56 - t411 * t439 + t197 * t33 + t182 * t55 + t177 * t53 + t357 * t433 - pkin(2) * t113 - t322 * t397 + (t292 * t380 + t204 / 0.2e1) * t78 + t570 * t107 + (t139 * mrSges(4,2) * t291 - t138 * (mrSges(4,1) * t291 - mrSges(4,3) * t420)) * t413 - t535 + (t32 - t151) * t284; (-t323 / 0.2e1 - t79 / 0.2e1 - t80 / 0.2e1 - t75 / 0.2e1 + t481 + (-Ifges(5,3) / 0.2e1 - Ifges(7,1) / 0.2e1 - Ifges(6,1) / 0.2e1) * t211 + (Ifges(6,4) / 0.2e1 - Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * t169 + (Ifges(5,6) / 0.2e1 - Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t168 + t575) * t216 + (-t210 / 0.2e1 - t134 / 0.2e1 - t195 * mrSges(4,2) + t439 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t216 + t324 / 0.2e1 - t497) * t215 + t497 * qJD(4) + t524 * t122 + t525 * t124 + t526 * t108 + (t12 * t258 - t36 * t48 - t37 * t47 + t40 * t526) * m(6) + (t1 * t266 + t17 * t524 + t2 * t267 + t246 * t3 + t27 * t525 + t28 * t527) * m(7) + t527 * t105 + t401 + t538 * t466 - t539 * t292 / 0.2e1 + t512 * t476 + (t349 + t513) * t488 + t23 * t464 + t335 * t489 + t289 * t493 + (t180 - t107) * t139 + (-pkin(3) * t62 - t120 * t139 + t41 * t66 - t42 * t67) * m(5) + t342 * t475 + (t1 * t289 + t2 * t292) * mrSges(7,1) + (t352 + t340) * t490 + t508 * mrSges(6,1) + (m(6) * t508 + m(5) * t509 + (-m(5) * t330 + m(6) * t331 + t289 * t418 + t292 * t417) * qJD(4) + (t52 - t54) * t292 + (-t51 + t56) * t289) * pkin(10) + t509 * mrSges(5,3) + t266 * t53 + t267 * t55 + t258 * t31 + t246 * t33 - t138 * t179 - t48 * t125 - t67 * t126 - t66 * t127 - t47 * t123 - t71 * mrSges(4,2) + t72 * mrSges(4,1) - pkin(3) * t32 + t62 * (-mrSges(5,1) * t292 + mrSges(5,2) * t289) + t12 * (mrSges(6,2) * t292 - mrSges(6,3) * t289) + t3 * (-mrSges(7,2) * t289 - mrSges(7,3) * t292); t552 + (qJ(5) * t2 - t1 * t454 + t520 * t17 + t27 * t557 - t28 * t58) * m(7) + (-Ifges(5,2) * t169 - t165 + t536) * t479 + t520 * t122 - t454 * t53 + t504 + (-t54 + t55) * qJ(5) + (t168 * t17 + t169 * t27) * mrSges(7,1) + (t168 * t36 - t169 * t37) * mrSges(6,1) + (-t417 + t452) * t42 + (-t418 + t453) * t41 - t120 * (mrSges(5,1) * t169 - mrSges(5,2) * t168) - t28 * (mrSges(7,2) * t168 + mrSges(7,3) * t169) - t40 * (-mrSges(6,2) * t169 + mrSges(6,3) * t168) - t419 * qJD(5) + t326 * t124 - t106 * t108 - t58 * t105 - pkin(4) * t56 + (-pkin(4) * t7 - qJ(5) * t5 - t106 * t40 - t36 * t42 + t37 * t521) * m(6) + (t169 * t541 + t164 - t437 + t77) * t480 + (t168 * t578 + t169 * t577) * t469 + (-t168 * t545 + t163 - t436 + t537) * t478 + (Ifges(6,2) * t168 + t435 + t78) * t477; t419 * t211 + (t105 + t108) * t169 + t53 + t56 + (t169 * t28 - t211 * t27 + t1) * m(7) + (t169 * t40 + t211 * t37 + t7) * m(6); -t168 * t105 + t211 * t122 + 0.2e1 * (t2 / 0.2e1 + t28 * t480 + t17 * t468) * m(7) + t55;];
tauc  = t18(:);
