% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:08
% EndTime: 2019-03-09 23:52:12
% DurationCPUTime: 32.87s
% Computational Cost: add. (18738->944), mult. (47931->1261), div. (0->0), fcn. (36745->10), ass. (0->403)
t344 = sin(qJ(2));
t339 = sin(pkin(6));
t430 = qJD(1) * t339;
t413 = t344 * t430;
t340 = cos(pkin(6));
t348 = cos(qJ(2));
t483 = pkin(1) * t348;
t421 = t340 * t483;
t287 = -pkin(8) * t413 + qJD(1) * t421;
t371 = (pkin(2) * t344 - pkin(9) * t348) * t339;
t288 = qJD(1) * t371;
t343 = sin(qJ(3));
t347 = cos(qJ(3));
t203 = t347 * t287 + t343 * t288;
t190 = pkin(10) * t413 + t203;
t333 = t340 * t344 * pkin(1);
t403 = pkin(3) * t343 - pkin(10) * t347;
t444 = t339 * t348;
t206 = (t333 + (pkin(8) + t403) * t444) * qJD(1);
t342 = sin(qJ(4));
t346 = cos(qJ(4));
t113 = t346 * t190 + t342 * t206;
t412 = t348 * t430;
t406 = t343 * t412;
t427 = qJD(3) * t343;
t582 = -t113 + (-t406 + t427) * qJ(5);
t311 = t403 * qJD(3);
t323 = -pkin(3) * t347 - pkin(10) * t343 - pkin(2);
t442 = t346 * t347;
t335 = pkin(9) * t442;
t425 = qJD(4) * t342;
t579 = -qJD(4) * t335 + t342 * t190 - t323 * t425 + (-t206 + t311) * t346;
t553 = Ifges(5,1) + Ifges(6,1);
t552 = Ifges(6,4) + Ifges(5,5);
t441 = t346 * t348;
t249 = (t342 * t344 + t347 * t441) * t430;
t482 = pkin(9) * t342;
t414 = -pkin(4) - t482;
t480 = pkin(11) * t343;
t420 = t342 * t480;
t518 = pkin(4) + pkin(5);
t581 = pkin(11) * t249 + t406 * t518 + qJD(4) * t420 + (-pkin(11) * t442 + (-pkin(5) + t414) * t343) * qJD(3) - t579;
t443 = t342 * t347;
t248 = -t346 * t413 + t412 * t443;
t424 = qJD(4) * t346;
t432 = t342 * t311 + t323 * t424;
t580 = -pkin(11) * t248 + (-pkin(9) * qJD(3) + pkin(11) * qJD(4)) * t346 * t343 + (-qJD(5) + (-pkin(9) * qJD(4) + pkin(11) * qJD(3)) * t342) * t347 + t432 + t582;
t551 = Ifges(6,5) - Ifges(5,4);
t578 = Ifges(5,6) - Ifges(6,6);
t329 = qJD(1) * t340 + qJD(2);
t235 = -pkin(2) * t329 - t287;
t263 = t347 * t329 - t343 * t413;
t264 = t329 * t343 + t347 * t413;
t149 = -pkin(3) * t263 - pkin(10) * t264 + t235;
t431 = pkin(8) * t444 + t333;
t290 = t431 * qJD(1);
t236 = pkin(9) * t329 + t290;
t280 = (-pkin(2) * t348 - pkin(9) * t344 - pkin(1)) * t339;
t254 = qJD(1) * t280;
t176 = t347 * t236 + t343 * t254;
t382 = -qJD(3) + t412;
t153 = -pkin(10) * t382 + t176;
t69 = t346 * t149 - t342 * t153;
t534 = qJD(5) - t69;
t334 = pkin(9) * t443;
t338 = t347 * pkin(4);
t215 = pkin(5) * t347 + t334 + t338 + (-t323 - t480) * t346;
t274 = t342 * t323 + t335;
t240 = -qJ(5) * t347 + t274;
t220 = t240 + t420;
t341 = sin(qJ(6));
t345 = cos(qJ(6));
t142 = t215 * t345 - t220 * t341;
t577 = qJD(6) * t142 + t341 * t581 + t345 * t580;
t143 = t215 * t341 + t220 * t345;
t576 = -qJD(6) * t143 - t341 * t580 + t345 * t581;
t481 = pkin(11) * t263;
t517 = pkin(10) - pkin(11);
t175 = -t343 * t236 + t347 * t254;
t198 = pkin(3) * t264 - pkin(10) * t263;
t95 = t346 * t175 + t342 * t198;
t74 = t264 * qJ(5) + t95;
t575 = t342 * t481 + t517 * t425 + t74;
t161 = t342 * t175;
t325 = t517 * t346;
t574 = qJD(4) * t325 - t161 - (-t198 - t481) * t346 + t518 * t264;
t202 = -t343 * t287 + t347 * t288;
t366 = pkin(3) * t413 + t202;
t361 = qJ(5) * t249 + t366;
t446 = qJ(5) * t346;
t370 = -t342 * t518 + t446;
t364 = -pkin(9) + t370;
t422 = qJD(5) * t346;
t426 = qJD(3) * t347;
t447 = qJ(5) * t342;
t525 = -t346 * t518 - t447;
t573 = t248 * t518 + (qJD(4) * t525 + t422) * t343 + t364 * t426 - t361;
t196 = (-t346 * t427 - t347 * t425) * pkin(9) + t432;
t572 = -qJD(5) * t347 + t196 + t582;
t571 = pkin(4) * t406 + t414 * t427 - t579;
t208 = t346 * t264 - t342 * t382;
t570 = -pkin(11) * t208 + t534;
t207 = t342 * t264 + t346 * t382;
t205 = Ifges(5,4) * t207;
t257 = qJD(4) - t263;
t468 = Ifges(6,5) * t207;
t533 = t208 * t553 + t552 * t257 - t205 + t468;
t383 = pkin(4) * t342 - t446;
t373 = pkin(9) + t383;
t384 = pkin(4) * t346 + t447;
t569 = -pkin(4) * t248 + (qJD(4) * t384 - t422) * t343 + t373 * t426 + t361;
t568 = t427 * t482 + t579;
t567 = -t113 + t196;
t129 = t207 * t345 - t208 * t341;
t127 = Ifges(7,4) * t129;
t377 = t207 * t341 + t208 * t345;
t566 = Ifges(7,2) * t377 - t127;
t465 = Ifges(4,2) * t263;
t102 = t208 * Ifges(5,5) - t207 * Ifges(5,6) + t257 * Ifges(5,3);
t103 = t208 * Ifges(6,4) + t257 * Ifges(6,2) + t207 * Ifges(6,6);
t41 = -t257 * t518 + t570;
t243 = t257 * qJ(5);
t70 = t342 * t149 + t346 * t153;
t58 = pkin(11) * t207 + t70;
t52 = t243 + t58;
t11 = -t341 * t52 + t345 * t41;
t12 = t341 * t41 + t345 * t52;
t368 = Ifges(4,6) * t382;
t473 = Ifges(4,4) * t264;
t170 = -t368 + t465 + t473;
t417 = Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t418 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t419 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t244 = qJD(6) - t257;
t451 = t244 * Ifges(7,3);
t460 = t377 * Ifges(7,5);
t461 = t129 * Ifges(7,6);
t53 = t451 + t460 + t461;
t63 = -pkin(4) * t257 + t534;
t64 = t243 + t70;
t523 = t418 * t207 + t419 * t208 + t417 * t257 - t11 * mrSges(7,1) - t176 * mrSges(4,3) - t63 * mrSges(6,1) - t70 * mrSges(5,2) + t102 / 0.2e1 + t103 / 0.2e1 - t170 / 0.2e1 - t53 / 0.2e1 - t473 / 0.2e1 + t12 * mrSges(7,2) - t461 / 0.2e1 - t460 / 0.2e1 + t235 * mrSges(4,1) - t451 / 0.2e1 + t64 * mrSges(6,3) + t69 * mrSges(5,1);
t565 = t465 / 0.2e1 - t523;
t497 = -t244 / 0.2e1;
t428 = qJD(2) * t348;
t410 = t343 * t428;
t217 = t329 * t427 + (t344 * t426 + t410) * t430;
t409 = t347 * t428;
t216 = t329 * t426 + (-t344 * t427 + t409) * t430;
t429 = qJD(2) * t339;
t408 = qJD(1) * t429;
t405 = t344 * t408;
t121 = -qJD(4) * t207 + t346 * t216 + t342 * t405;
t122 = qJD(4) * t208 + t342 * t216 - t346 * t405;
t40 = -qJD(6) * t377 - t121 * t341 + t122 * t345;
t37 = Ifges(7,6) * t40;
t39 = qJD(6) * t129 + t121 * t345 + t122 * t341;
t38 = Ifges(7,5) * t39;
t5 = -Ifges(7,3) * t217 + t37 + t38;
t152 = pkin(3) * t382 - t175;
t359 = t208 * qJ(5) - t152;
t56 = -t207 * t518 + t359;
t564 = (t11 * t129 + t12 * t377) * mrSges(7,3) + (Ifges(7,5) * t129 - Ifges(7,6) * t377) * t497 - t56 * (mrSges(7,1) * t377 + mrSges(7,2) * t129) + t5;
t204 = Ifges(6,5) * t208;
t101 = t257 * Ifges(6,6) + t207 * Ifges(6,3) + t204;
t456 = t208 * Ifges(5,4);
t104 = -t207 * Ifges(5,2) + t257 * Ifges(5,6) + t456;
t68 = t207 * pkin(4) - t359;
t563 = t152 * mrSges(5,1) + t68 * mrSges(6,1) + t101 / 0.2e1 - t104 / 0.2e1 - t64 * mrSges(6,2) - t70 * mrSges(5,3);
t547 = t121 * t553 + t551 * t122 + t552 * t217;
t562 = qJD(5) * t342 + t176;
t561 = -t342 * t578 + t346 * t552;
t467 = Ifges(6,5) * t342;
t472 = Ifges(5,4) * t342;
t560 = t346 * t553 + t467 - t472;
t470 = Ifges(7,4) * t377;
t558 = Ifges(7,1) * t129 - t470;
t520 = t39 / 0.2e1;
t519 = t40 / 0.2e1;
t55 = Ifges(7,1) * t377 + Ifges(7,5) * t244 + t127;
t555 = t55 / 0.2e1;
t502 = -t217 / 0.2e1;
t509 = -t377 / 0.2e1;
t407 = -t430 / 0.2e1;
t554 = -t329 * Ifges(3,6) / 0.2e1;
t550 = Ifges(6,2) + Ifges(5,3);
t45 = Ifges(5,5) * t121 - Ifges(5,6) * t122 + Ifges(5,3) * t217;
t46 = Ifges(6,4) * t121 + Ifges(6,2) * t217 + Ifges(6,6) * t122;
t548 = t45 + t46;
t545 = Ifges(4,5) * t216;
t445 = t339 * t344;
t330 = pkin(8) * t445;
t300 = -t330 + t421;
t291 = t300 * qJD(2);
t276 = qJD(1) * t291;
t541 = t276 * mrSges(3,2);
t315 = t345 * qJ(5) - t341 * t518;
t540 = -qJD(6) * t315 - t341 * t570 - t345 * t58;
t314 = -t341 * qJ(5) - t345 * t518;
t539 = qJD(6) * t314 - t341 * t58 + t345 * t570;
t324 = t517 * t342;
t239 = t324 * t341 + t325 * t345;
t538 = -qJD(6) * t239 + t341 * t575 + t345 * t574;
t238 = t324 * t345 - t325 * t341;
t537 = qJD(6) * t238 + t341 * t574 - t345 * t575;
t536 = t257 * t370 + t562;
t535 = t257 * t383 - t562;
t376 = t341 * t346 - t342 * t345;
t293 = t376 * t343;
t256 = Ifges(4,4) * t263;
t369 = Ifges(4,5) * t382;
t475 = Ifges(4,1) * t264;
t171 = t256 - t369 + t475;
t532 = t175 * mrSges(4,3) - t171 / 0.2e1 - t235 * mrSges(4,2) - t256 / 0.2e1;
t531 = t342 * t552 + t346 * t578;
t466 = Ifges(6,5) * t346;
t471 = Ifges(5,4) * t346;
t530 = t342 * t553 - t466 + t471;
t292 = t431 * qJD(2);
t277 = qJD(1) * t292;
t115 = pkin(3) * t217 - pkin(10) * t216 + t277;
t289 = qJD(2) * t371;
t275 = qJD(1) * t289;
t99 = -t236 * t427 + t254 * t426 + t343 * t275 + t347 * t276;
t87 = pkin(10) * t405 + t99;
t20 = t342 * t115 + t149 * t424 - t153 * t425 + t346 * t87;
t21 = t115 * t346 - t149 * t425 - t153 * t424 - t342 * t87;
t529 = t20 * t346 - t21 * t342;
t15 = t217 * qJ(5) + t257 * qJD(5) + t20;
t17 = -pkin(4) * t217 - t21;
t528 = t15 * t346 + t17 * t342;
t100 = -t236 * t426 - t254 * t427 + t347 * t275 - t343 * t276;
t527 = -t100 * mrSges(4,1) + t99 * mrSges(4,2);
t526 = qJD(4) - qJD(6);
t378 = t342 * t70 + t346 * t69;
t379 = t342 * t64 - t346 * t63;
t386 = Ifges(6,3) * t342 + t466;
t393 = -Ifges(5,2) * t342 + t471;
t398 = mrSges(6,1) * t342 - mrSges(6,3) * t346;
t400 = mrSges(5,1) * t342 + mrSges(5,2) * t346;
t484 = t346 / 0.2e1;
t486 = t342 / 0.2e1;
t487 = -t342 / 0.2e1;
t492 = t257 / 0.2e1;
t503 = t208 / 0.2e1;
t505 = t207 / 0.2e1;
t506 = -t207 / 0.2e1;
t524 = t379 * mrSges(6,2) + t378 * mrSges(5,3) - t101 * t486 - t104 * t487 - t152 * t400 - t386 * t505 - t393 * t506 - t398 * t68 - t533 * t484 - t492 * t561 - t503 * t560;
t522 = Ifges(7,4) * t520 + Ifges(7,2) * t519 + Ifges(7,6) * t502;
t521 = Ifges(7,1) * t520 + Ifges(7,4) * t519 + Ifges(7,5) * t502;
t516 = pkin(1) * mrSges(3,1);
t515 = pkin(1) * mrSges(3,2);
t514 = t121 / 0.2e1;
t513 = -t122 / 0.2e1;
t512 = t122 / 0.2e1;
t511 = -t129 / 0.2e1;
t510 = t129 / 0.2e1;
t508 = t377 / 0.2e1;
t504 = -t208 / 0.2e1;
t501 = t217 / 0.2e1;
t496 = t244 / 0.2e1;
t493 = -t257 / 0.2e1;
t297 = -t340 * t347 + t343 * t445;
t491 = -t297 / 0.2e1;
t298 = t340 * t343 + t347 * t445;
t489 = t298 / 0.2e1;
t488 = t340 / 0.2e1;
t485 = -t346 / 0.2e1;
t478 = qJD(3) / 0.2e1;
t477 = mrSges(5,3) * t207;
t476 = mrSges(5,3) * t208;
t474 = Ifges(3,4) * t344;
t463 = Ifges(4,3) * t344;
t454 = t216 * Ifges(4,1);
t453 = t216 * Ifges(4,4);
t452 = t217 * Ifges(4,4);
t448 = qJ(5) * t207;
t154 = -mrSges(6,2) * t207 + mrSges(6,3) * t257;
t155 = -mrSges(5,2) * t257 - t477;
t440 = -t154 - t155;
t156 = mrSges(5,1) * t257 - t476;
t157 = -mrSges(6,1) * t257 + mrSges(6,2) * t208;
t439 = -t156 + t157;
t375 = t341 * t342 + t345 * t346;
t163 = t293 * t526 + t375 * t426;
t182 = t248 * t341 + t249 * t345;
t438 = t163 - t182;
t223 = t526 * t375;
t164 = t223 * t343 - t376 * t426;
t181 = t248 * t345 - t249 * t341;
t437 = t164 - t181;
t278 = t330 + (-pkin(2) - t483) * t340;
t185 = pkin(3) * t297 - pkin(10) * t298 + t278;
t279 = pkin(9) * t340 + t431;
t200 = t347 * t279 + t343 * t280;
t187 = -pkin(10) * t444 + t200;
t90 = t342 * t185 + t346 * t187;
t179 = t376 * t263;
t224 = t526 * t376;
t436 = -t179 + t224;
t180 = t375 * t263;
t435 = t180 - t223;
t434 = -mrSges(3,1) * t329 - mrSges(4,1) * t263 + mrSges(4,2) * t264 + mrSges(3,3) * t413;
t135 = mrSges(5,1) * t207 + mrSges(5,2) * t208;
t219 = -mrSges(4,1) * t382 - t264 * mrSges(4,3);
t433 = -t219 + t135;
t199 = -t343 * t279 + t347 * t280;
t416 = t342 * t444;
t77 = t297 * qJ(5) + t90;
t415 = -Ifges(4,6) * t217 + Ifges(4,3) * t405 + t545;
t186 = pkin(3) * t444 - t199;
t411 = t344 * t429;
t80 = -t217 * mrSges(6,1) + t121 * mrSges(6,2);
t94 = t198 * t346 - t161;
t89 = t185 * t346 - t342 * t187;
t273 = t323 * t346 - t334;
t124 = -t279 * t426 - t280 * t427 + t347 * t289 - t343 * t291;
t10 = pkin(11) * t122 + t15;
t9 = -pkin(11) * t121 - t217 * t518 - t21;
t1 = qJD(6) * t11 + t10 * t345 + t341 * t9;
t2 = -qJD(6) * t12 - t10 * t341 + t345 * t9;
t404 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t401 = mrSges(5,1) * t346 - mrSges(5,2) * t342;
t399 = mrSges(6,1) * t346 + mrSges(6,3) * t342;
t392 = Ifges(5,2) * t346 + t472;
t226 = qJD(3) * t298 + t339 * t410;
t227 = -qJD(3) * t297 + t339 * t409;
t389 = Ifges(4,5) * t227 - Ifges(4,6) * t226;
t385 = -Ifges(6,3) * t346 + t467;
t229 = t298 * t346 - t416;
t62 = -pkin(11) * t229 - t297 * t518 - t89;
t228 = t298 * t342 + t339 * t441;
t66 = pkin(11) * t228 + t77;
t24 = -t341 * t66 + t345 * t62;
t26 = t341 * t62 + t345 * t66;
t92 = -mrSges(7,2) * t244 + mrSges(7,3) * t129;
t93 = mrSges(7,1) * t244 - mrSges(7,3) * t377;
t380 = -t341 * t93 + t345 * t92;
t158 = t228 * t345 - t229 * t341;
t159 = t228 * t341 + t229 * t345;
t374 = t412 - qJD(3) / 0.2e1;
t372 = qJ(5) * t229 - t186;
t123 = -t279 * t427 + t280 * t426 + t343 * t289 + t347 * t291;
t110 = pkin(10) * t411 + t123;
t136 = pkin(3) * t226 - pkin(10) * t227 + t292;
t33 = -t342 * t110 + t136 * t346 - t185 * t425 - t187 * t424;
t32 = t346 * t110 + t342 * t136 + t185 * t424 - t187 * t425;
t326 = Ifges(3,4) * t412;
t363 = -t287 * mrSges(3,3) + Ifges(3,1) * t413 / 0.2e1 + t326 / 0.2e1 + t329 * Ifges(3,5);
t362 = pkin(3) * t411 + t124;
t23 = t226 * qJ(5) + t297 * qJD(5) + t32;
t360 = pkin(3) * t405 + t100;
t358 = -t475 / 0.2e1 + t532;
t357 = -t21 * mrSges(5,1) + t17 * mrSges(6,1) + t20 * mrSges(5,2) - t15 * mrSges(6,3) + t404;
t146 = -qJD(4) * t228 + t227 * t346 + t342 * t411;
t356 = qJ(5) * t146 + qJD(5) * t229 + t362;
t355 = qJ(5) * t121 + qJD(5) * t208 + t360;
t354 = t175 * mrSges(4,1) - Ifges(4,3) * t382 / 0.2e1 + Ifges(4,6) * t263 + Ifges(4,5) * t264 + t554 + (Ifges(3,2) * t348 + t474) * t407 - t176 * mrSges(4,2) - t290 * mrSges(3,3);
t322 = Ifges(3,5) * t348 * t408;
t317 = -pkin(3) - t384;
t304 = pkin(3) - t525;
t294 = t375 * t343;
t285 = -mrSges(3,2) * t329 + mrSges(3,3) * t412;
t281 = t373 * t343;
t241 = -t273 + t338;
t237 = t364 * t343;
t218 = mrSges(4,2) * t382 + t263 * mrSges(4,3);
t193 = -mrSges(4,2) * t405 - mrSges(4,3) * t217;
t192 = mrSges(4,1) * t405 - mrSges(4,3) * t216;
t145 = -qJD(4) * t416 + t227 * t342 + t298 * t424 - t346 * t411;
t139 = mrSges(4,1) * t217 + mrSges(4,2) * t216;
t134 = mrSges(6,1) * t207 - mrSges(6,3) * t208;
t133 = pkin(4) * t208 + t448;
t126 = Ifges(4,5) * t405 - t452 + t454;
t125 = -t217 * Ifges(4,2) + Ifges(4,6) * t405 + t453;
t91 = pkin(4) * t228 - t372;
t84 = -t208 * t518 - t448;
t82 = -mrSges(6,2) * t122 + mrSges(6,3) * t217;
t81 = -mrSges(5,2) * t217 - mrSges(5,3) * t122;
t79 = mrSges(5,1) * t217 - mrSges(5,3) * t121;
t78 = -pkin(4) * t297 - t89;
t76 = -pkin(4) * t264 - t94;
t72 = -t228 * t518 + t372;
t65 = -mrSges(7,1) * t129 + mrSges(7,2) * t377;
t60 = mrSges(5,1) * t122 + mrSges(5,2) * t121;
t59 = mrSges(6,1) * t122 - mrSges(6,3) * t121;
t54 = Ifges(7,2) * t129 + Ifges(7,6) * t244 + t470;
t51 = qJD(6) * t158 + t145 * t341 + t146 * t345;
t50 = -qJD(6) * t159 + t145 * t345 - t146 * t341;
t47 = Ifges(5,4) * t121 - Ifges(5,2) * t122 + Ifges(5,6) * t217;
t44 = Ifges(6,5) * t121 + Ifges(6,6) * t217 + Ifges(6,3) * t122;
t36 = pkin(4) * t145 - t356;
t31 = mrSges(7,2) * t217 + mrSges(7,3) * t40;
t30 = -mrSges(7,1) * t217 - mrSges(7,3) * t39;
t29 = -pkin(4) * t226 - t33;
t28 = pkin(4) * t122 - t355;
t22 = -t145 * t518 + t356;
t16 = pkin(11) * t145 + t23;
t14 = -pkin(11) * t146 - t226 * t518 - t33;
t13 = -t122 * t518 + t355;
t8 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t4 = -qJD(6) * t26 + t14 * t345 - t16 * t341;
t3 = qJD(6) * t24 + t14 * t341 + t16 * t345;
t6 = [(-Ifges(5,2) * t506 + Ifges(6,3) * t505 - t492 * t578 + t551 * t503 + t563) * t145 + (Ifges(6,3) * t512 - Ifges(5,2) * t513 + t44 / 0.2e1 - t47 / 0.2e1 + t28 * mrSges(6,1) - t360 * mrSges(5,1) - t15 * mrSges(6,2) - t20 * mrSges(5,3) + t551 * t514 - t578 * t501) * t228 + t216 * (Ifges(4,1) * t298 - Ifges(4,4) * t297) / 0.2e1 + (t103 + t102) * t226 / 0.2e1 - (t170 + t53) * t226 / 0.2e1 + (-t348 * t389 / 0.2e1 + ((Ifges(3,5) * t488 - t300 * mrSges(3,3) + (-0.2e1 * t515 + 0.3e1 / 0.2e1 * Ifges(3,4) * t348) * t339) * t348 + (Ifges(4,5) * t489 + Ifges(4,6) * t491 - Ifges(3,6) * t340 - t431 * mrSges(3,3) + (-0.2e1 * t516 - 0.3e1 / 0.2e1 * t474 + (-Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t348) * t339) * t344) * qJD(2)) * t430 + m(3) * (t276 * t431 - t277 * t300 - t287 * t292 + t290 * t291) + (Ifges(4,4) * t298 + Ifges(7,5) * t159 + Ifges(7,6) * t158 + (-Ifges(4,2) - Ifges(7,3)) * t297) * t502 + (Ifges(6,5) * t146 + Ifges(6,6) * t226) * t505 + (Ifges(6,5) * t229 + Ifges(6,6) * t297) * t512 + (Ifges(5,4) * t146 + Ifges(5,6) * t226) * t506 + (Ifges(5,4) * t229 + Ifges(5,6) * t297) * t513 + t36 * t134 + (-t146 * t68 + t15 * t297 + t226 * t64 - t229 * t28) * mrSges(6,3) + t125 * t491 + t5 * t491 + (Ifges(7,5) * t51 + Ifges(7,6) * t50 - Ifges(7,3) * t226) * t496 + (t146 * t552 + t226 * t550) * t492 + (t229 * t552 + t297 * t550) * t501 + (t146 * t553 + t226 * t552) * t503 + (t229 * t553 + t297 * t552) * t514 + t389 * t478 + t322 * t488 + t126 * t489 + t264 * (Ifges(4,1) * t227 - Ifges(4,4) * t226) / 0.2e1 + (-Ifges(4,6) * t502 - t415 / 0.2e1 - t545 / 0.2e1 + t276 * mrSges(3,3) + t527) * t444 + (t146 * t152 - t20 * t297 - t226 * t70 - t229 * t360) * mrSges(5,2) - t362 * t135 + m(5) * (-t152 * t362 - t186 * t360 + t20 * t90 + t21 * t89 + t32 * t70 + t33 * t69) + (t363 * t348 + (Ifges(4,3) * t478 + t354 + t554) * t344) * t429 + (Ifges(7,1) * t51 + Ifges(7,4) * t50 - Ifges(7,5) * t226) * t508 + (Ifges(7,4) * t51 + Ifges(7,2) * t50 - Ifges(7,6) * t226) * t510 + t533 * t146 / 0.2e1 + t89 * t79 + t90 * t81 + t91 * t59 + t3 * t92 + t4 * t93 - t340 * t541 + t547 * t229 / 0.2e1 + t548 * t297 / 0.2e1 + t263 * (Ifges(4,4) * t227 - Ifges(4,2) * t226) / 0.2e1 + t77 * t82 + t78 * t80 + t72 * t8 + t22 * t65 + t50 * t54 / 0.2e1 + t56 * (-mrSges(7,1) * t50 + mrSges(7,2) * t51) + t26 * t31 + t24 * t30 + (-t100 * t298 - t175 * t227 - t176 * t226 - t297 * t99) * mrSges(4,3) + m(4) * (t100 * t199 + t123 * t176 + t124 * t175 + t200 * t99 + t235 * t292 + t277 * t278) + m(7) * (t1 * t26 + t11 * t4 + t12 * t3 + t13 * t72 + t2 * t24 + t22 * t56) + m(6) * (t15 * t77 + t17 * t78 + t23 * t64 + t28 * t91 + t29 * t63 + t36 * t68) + (-mrSges(3,1) * t340 + mrSges(4,1) * t297 + mrSges(4,2) * t298 + mrSges(3,3) * t445) * t277 + t434 * t292 + t23 * t154 + t32 * t155 + t33 * t156 + t29 * t157 + t13 * (-mrSges(7,1) * t158 + mrSges(7,2) * t159) + t186 * t60 + t199 * t192 + t200 * t193 + t123 * t218 + t124 * t219 + (Ifges(7,4) * t159 + Ifges(7,2) * t158 - Ifges(7,6) * t297) * t519 + (Ifges(7,1) * t159 + Ifges(7,4) * t158 - Ifges(7,5) * t297) * t520 + t159 * t521 + t158 * t522 + t69 * (mrSges(5,1) * t226 - mrSges(5,3) * t146) + t63 * (-mrSges(6,1) * t226 + mrSges(6,2) * t146) + t11 * (-mrSges(7,1) * t226 - mrSges(7,3) * t51) + t12 * (mrSges(7,2) * t226 + mrSges(7,3) * t50) + t227 * t171 / 0.2e1 + t235 * (mrSges(4,1) * t226 + mrSges(4,2) * t227) + t278 * t139 + t291 * t285 + t21 * (mrSges(5,1) * t297 - mrSges(5,3) * t229) + t17 * (-mrSges(6,1) * t297 + mrSges(6,2) * t229) + t1 * (mrSges(7,2) * t297 + mrSges(7,3) * t158) + t2 * (-mrSges(7,1) * t297 - mrSges(7,3) * t159) + t51 * t555; (-pkin(2) * t277 - t175 * t202 - t176 * t203 - t235 * t290) * m(4) + (Ifges(7,5) * t294 - Ifges(7,6) * t293) * t502 + (-t1 * t293 - t11 * t438 + t12 * t437 - t2 * t294) * mrSges(7,3) + (Ifges(7,4) * t294 - Ifges(7,2) * t293) * t519 + (Ifges(7,1) * t294 - Ifges(7,4) * t293) * t520 + t13 * (mrSges(7,1) * t293 + mrSges(7,2) * t294) + t322 - t541 + (Ifges(7,5) * t163 + Ifges(7,6) * t164) * t496 + (Ifges(7,5) * t182 + Ifges(7,6) * t181) * t497 + t366 * t135 + ((-Ifges(7,3) / 0.2e1 - Ifges(4,2) / 0.2e1 - t417) * t217 - t418 * t122 - t46 / 0.2e1 - t45 / 0.2e1 - t419 * t121 + t357 - t277 * mrSges(4,1) + t99 * mrSges(4,3) + t37 / 0.2e1 + t38 / 0.2e1 + t125 / 0.2e1 + t5 / 0.2e1 + t453 / 0.2e1) * t347 + (t163 / 0.2e1 - t182 / 0.2e1) * t55 + (Ifges(7,1) * t163 + Ifges(7,4) * t164) * t508 + (Ifges(7,1) * t182 + Ifges(7,4) * t181) * t509 + (Ifges(7,4) * t163 + Ifges(7,2) * t164) * t510 + (Ifges(7,4) * t182 + Ifges(7,2) * t181) * t511 + (t164 / 0.2e1 - t181 / 0.2e1) * t54 + (-t533 / 0.2e1 - mrSges(5,2) * t152 - t63 * mrSges(6,2) + t69 * mrSges(5,3) + mrSges(6,3) * t68 + Ifges(5,4) * t505 + Ifges(6,5) * t506 + t552 * t493 + t553 * t504) * t249 - pkin(2) * t139 + t142 * t30 + t143 * t31 + (-mrSges(7,1) * t437 + mrSges(7,2) * t438) * t56 + (t386 * t512 + t393 * t513 - t100 * mrSges(4,3) + t44 * t486 - t360 * t400 + t28 * t398 + t454 / 0.2e1 - t452 / 0.2e1 + t277 * mrSges(4,2) + t47 * t487 + t126 / 0.2e1 + (-t20 * t342 - t21 * t346) * mrSges(5,3) + (-t15 * t342 + t17 * t346) * mrSges(6,2) + (t385 * t506 + t392 * t505 + t68 * t399 + t152 * t401 + t104 * t485 + (t342 * t69 - t346 * t70) * mrSges(5,3) + (-t342 * t63 - t346 * t64) * mrSges(6,2) + t530 * t504 + t531 * t493 + t533 * t487) * qJD(4) + t560 * t514 + t561 * t501 + (qJD(4) * t101 + t547) * t484) * t343 + (m(4) * (-t100 * t343 - t175 * t426 - t176 * t427 + t347 * t99) + m(5) * (t152 * t426 - t343 * t360) + (-t218 * t343 + t347 * t433) * qJD(3) + t193 * t347 + (t60 - t192) * t343) * pkin(9) - t434 * t290 + ((qJD(2) * (Ifges(4,5) * t343 + Ifges(4,6) * t347) / 0.2e1 + (t516 + t474 / 0.2e1) * t430 + (t329 / 0.2e1 - qJD(2)) * Ifges(3,6) - t354) * t344 + (-t326 / 0.2e1 + (t515 + t463 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t344) * t430 + (Ifges(4,5) * t412 / 0.2e1 + t358) * t347 + (Ifges(4,6) * t348 * t407 + t565) * t343 - t363) * t348) * t430 + (t407 * t463 + (t374 * Ifges(4,6) - t565) * t343 + (-Ifges(4,5) * t374 - t358 - t524) * t347) * qJD(3) + t567 * t155 + t568 * t156 + (t152 * t366 + t20 * t274 + t21 * t273 + t567 * t70 + t568 * t69) * m(5) + t569 * t134 + t571 * t157 + t572 * t154 + (t15 * t240 + t17 * t241 + t28 * t281 + t569 * t68 + t571 * t63 + t572 * t64) * m(6) + t573 * t65 - t203 * t218 - t202 * t219 + t294 * t521 - t293 * t522 + t576 * t93 + t577 * t92 + (t1 * t143 + t11 * t576 + t12 * t577 + t13 * t237 + t142 * t2 + t56 * t573) * m(7) + (-Ifges(5,2) * t505 + Ifges(6,3) * t506 - t493 * t578 + t504 * t551 - t563) * t248 + t237 * t8 + t240 * t82 + t241 * t80 + t273 * t79 + t274 * t81 - t277 * mrSges(3,1) + t281 * t59 - t287 * t285; (Ifges(7,5) * t180 - Ifges(7,6) * t179) * t497 + (Ifges(7,1) * t180 - Ifges(7,4) * t179) * t509 + (Ifges(7,4) * t180 - Ifges(7,2) * t179) * t511 + (pkin(3) * t360 - t152 * t176 - t69 * t94 - t70 * t95) * m(5) + t360 * t401 - t527 + t47 * t484 + t44 * t485 - t524 * qJD(4) + (-t523 - t368 / 0.2e1) * t264 + t415 + t535 * t134 + (t28 * t317 + t535 * t68 - t63 * t76 - t64 * t74) * m(6) + t536 * t65 + t537 * t92 + t538 * t93 + (t1 * t239 + t11 * t538 + t12 * t537 + t13 * t304 + t2 * t238 + t536 * t56) * m(7) - t375 * t522 + (-Ifges(7,4) * t376 - Ifges(7,2) * t375) * t519 + (-Ifges(7,1) * t376 - Ifges(7,4) * t375) * t520 + t13 * (mrSges(7,1) * t375 - mrSges(7,2) * t376) + (-Ifges(7,5) * t376 - Ifges(7,6) * t375) * t502 + (-t1 * t375 + t11 * t435 - t12 * t436 + t2 * t376) * mrSges(7,3) - t376 * t521 + (t179 / 0.2e1 - t224 / 0.2e1) * t54 + (Ifges(7,5) * t223 - Ifges(7,6) * t224) * t496 + (Ifges(7,1) * t223 - Ifges(7,4) * t224) * t508 + (Ifges(7,4) * t223 - Ifges(7,2) * t224) * t510 + (-t180 / 0.2e1 + t223 / 0.2e1) * t55 + t531 * t501 + ((Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t264 + t369 / 0.2e1 + t524 + t532) * t263 + t547 * t486 + t528 * mrSges(6,2) + ((t81 + t82) * t346 + (-t79 + t80) * t342 + m(5) * t529 + m(6) * t528 + (-m(5) * t378 - m(6) * t379 + t342 * t440 + t346 * t439) * qJD(4)) * pkin(10) + t529 * mrSges(5,3) + t530 * t514 - pkin(3) * t60 - t433 * t176 + (mrSges(7,1) * t436 - mrSges(7,2) * t435) * t56 - t74 * t154 - t95 * t155 - t94 * t156 - t76 * t157 - t175 * t218 + t385 * t512 + t392 * t513 + t238 * t30 + t239 * t31 - t28 * t399 + t304 * t8 + t317 * t59; (-t207 * t552 - t208 * t578) * t493 + (t54 - t558) * t509 - t133 * t134 + (-t207 * t553 + t101 + t204 - t456) * t504 + t539 * t92 + (t1 * t315 + t11 * t540 + t12 * t539 + t2 * t314 - t56 * t84) * m(7) + t540 * t93 + t129 * t555 + t104 * t503 + (Ifges(6,3) * t208 - t468) * t506 + (-Ifges(5,2) * t208 - t205 + t533) * t505 + (-pkin(4) * t17 + qJ(5) * t15 - t133 * t68 + t534 * t64 - t63 * t70) * m(6) - t357 - t84 * t65 - pkin(4) * t80 + qJ(5) * t82 + (t207 * t63 + t208 * t64) * mrSges(6,2) + (-t439 + t476) * t70 + (t440 - t477) * t69 + qJD(5) * t154 + t566 * t511 - t564 + t548 - t68 * (mrSges(6,1) * t208 + mrSges(6,3) * t207) - t152 * (mrSges(5,1) * t208 - mrSges(5,2) * t207) + t314 * t30 + t315 * t31; t345 * t30 + t341 * t31 + (t134 - t65) * t208 + t380 * qJD(6) + (-t154 - t380) * t257 + t80 + (t1 * t341 + t2 * t345 - t208 * t56 + t244 * (-t11 * t341 + t12 * t345)) * m(7) + (t208 * t68 - t257 * t64 + t17) * m(6); t558 * t509 + t54 * t508 - t11 * t92 + t12 * t93 + t404 + (t55 - t566) * t511 + t564;];
tauc  = t6(:);
