% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:18
% EndTime: 2019-03-09 09:54:37
% DurationCPUTime: 9.84s
% Computational Cost: add. (13038->680), mult. (28612->862), div. (0->0), fcn. (28446->6), ass. (0->317)
t526 = m(7) / 0.2e1;
t574 = -0.2e1 * t526;
t566 = Ifges(5,5) + Ifges(7,5);
t543 = Ifges(6,4) - t566;
t567 = Ifges(7,4) + Ifges(6,5);
t573 = Ifges(5,6) - t567;
t528 = m(6) / 0.2e1;
t344 = sin(pkin(9));
t345 = cos(pkin(9));
t489 = sin(qJ(4));
t490 = cos(qJ(4));
t301 = t489 * t344 - t490 * t345;
t572 = -t301 / 0.2e1;
t348 = cos(qJ(2));
t428 = t348 * t301;
t556 = t428 / 0.2e1;
t571 = -Ifges(5,1) - Ifges(7,3);
t564 = Ifges(7,2) + Ifges(6,3);
t346 = pkin(4) + qJ(6);
t303 = t344 * t490 + t345 * t489;
t328 = -pkin(3) * t345 - pkin(2);
t368 = -qJ(5) * t303 + t328;
t123 = t301 * t346 + t368;
t204 = -mrSges(7,2) * t303 + mrSges(7,3) * t301;
t401 = m(7) * t123 + t204;
t483 = pkin(8) + qJ(3);
t313 = t483 * t345;
t402 = t483 * t344;
t219 = t313 * t489 + t490 * t402;
t220 = t313 * t490 - t402 * t489;
t347 = sin(qJ(2));
t277 = t303 * t347;
t279 = t301 * t347;
t545 = t220 * t279 / 0.2e1 - t219 * t277 / 0.2e1;
t341 = t348 * pkin(4);
t312 = -pkin(2) * t348 - t347 * qJ(3) - pkin(1);
t430 = t345 * t348;
t240 = pkin(7) * t430 + t344 * t312;
t433 = t344 * t347;
t214 = -pkin(8) * t433 + t240;
t358 = (-pkin(7) * t344 - pkin(3)) * t348 + (-pkin(8) * t347 + t312) * t345;
t92 = t214 * t489 - t490 * t358;
t81 = t341 + t92;
t570 = -(t92 / 0.2e1 - t81 / 0.2e1) * t301 - t545;
t568 = Ifges(5,4) + Ifges(6,6);
t565 = Ifges(6,2) + Ifges(5,1);
t184 = pkin(4) * t301 + t368;
t513 = t184 / 0.2e1;
t563 = mrSges(6,2) * t513;
t339 = t347 * pkin(7);
t310 = pkin(3) * t433 + t339;
t442 = qJ(5) * t279;
t390 = t310 + t442;
t121 = pkin(4) * t277 + t390;
t441 = qJ(5) * t301;
t152 = t303 * t346 + t441;
t472 = Ifges(7,6) * t277;
t154 = -t348 * Ifges(7,5) - Ifges(7,3) * t279 + t472;
t474 = Ifges(6,6) * t279;
t156 = -t348 * Ifges(6,5) + Ifges(6,3) * t277 + t474;
t263 = Ifges(7,6) * t279;
t158 = -t348 * Ifges(7,4) + Ifges(7,2) * t277 - t263;
t264 = Ifges(6,6) * t277;
t160 = -t348 * Ifges(6,4) + Ifges(6,2) * t279 + t264;
t478 = Ifges(5,4) * t279;
t161 = -Ifges(5,2) * t277 - t348 * Ifges(5,6) - t478;
t271 = Ifges(5,4) * t277;
t163 = -Ifges(5,1) * t279 - t348 * Ifges(5,5) - t271;
t166 = t303 * pkin(5) + t219;
t168 = -t301 * pkin(5) + t220;
t443 = qJ(5) * t277;
t170 = -pkin(4) * t279 + t443;
t205 = pkin(4) * t303 + t441;
t445 = t348 * mrSges(7,3);
t462 = t279 * mrSges(7,1);
t225 = t445 - t462;
t261 = t277 * mrSges(5,2);
t397 = t301 * mrSges(7,2) + t303 * mrSges(7,3);
t293 = t301 * mrSges(6,3);
t399 = -t303 * mrSges(6,2) + t293;
t294 = t301 * mrSges(5,2);
t400 = t303 * mrSges(5,1) - t294;
t495 = t328 / 0.2e1;
t496 = t303 / 0.4e1;
t498 = -t303 / 0.4e1;
t500 = t301 / 0.4e1;
t502 = -t301 / 0.4e1;
t337 = t348 * mrSges(7,2);
t468 = t277 * mrSges(7,1);
t227 = -t337 - t468;
t510 = -t227 / 0.2e1;
t206 = -mrSges(6,2) * t301 - mrSges(6,3) * t303;
t511 = -t206 / 0.2e1;
t512 = -t204 / 0.2e1;
t172 = -mrSges(6,2) * t277 + mrSges(6,3) * t279;
t514 = -t172 / 0.2e1;
t169 = mrSges(7,2) * t279 + mrSges(7,3) * t277;
t515 = -t169 / 0.2e1;
t527 = -m(7) / 0.2e1;
t529 = -m(6) / 0.2e1;
t429 = t348 * qJ(5);
t93 = t490 * t214 + t489 * t358;
t78 = -t93 + t429;
t549 = t78 + t93;
t63 = t279 * pkin(5) - t92;
t360 = -t341 + t63;
t48 = t348 * qJ(6) - t360;
t550 = t48 + t63;
t486 = t277 * pkin(5);
t56 = -t78 - t486;
t64 = t93 - t486;
t82 = t277 * t346 + t390;
t99 = -t279 * t346 + t443;
t561 = t261 * t495 + (t121 * t205 + t170 * t184 + (t81 - t92) * t220 + t549 * t219) * t529 + (t123 * t99 + t152 * t82 + t168 * t550) * t527 - t121 * t399 / 0.2e1 + t152 * t515 - t168 * t225 / 0.2e1 + t170 * t511 + t205 * t514 + t160 * t502 + t161 * t496 - t310 * t400 / 0.2e1 - t82 * t397 / 0.2e1 + t99 * t512 + (t154 + t163) * t500 + (t156 + t158) * t498 + (t301 * t543 - t303 * t573) * t348 / 0.4e1 + ((-t56 + t64) * t527 - t510) * t166;
t560 = -0.2e1 * t303;
t524 = -mrSges(5,1) / 0.2e1;
t559 = m(6) + m(7);
t278 = t348 * t303;
t557 = t278 / 0.2e1;
t555 = -t428 / 0.2e1;
t554 = -t347 / 0.2e1;
t553 = t347 / 0.2e1;
t485 = -mrSges(5,1) + mrSges(6,2);
t552 = mrSges(6,1) + mrSges(7,1);
t551 = mrSges(5,2) - mrSges(6,3);
t477 = Ifges(5,4) * t303;
t211 = -t301 * Ifges(5,2) + t477;
t473 = Ifges(6,6) * t303;
t546 = Ifges(6,2) * t301 + t211 + t473;
t314 = t347 * pkin(2) - qJ(3) * t348;
t258 = pkin(7) * t433 + t345 * t314;
t431 = t345 * t347;
t259 = -pkin(7) * t431 + t344 * t314;
t544 = -t258 * t344 + t259 * t345;
t297 = Ifges(6,6) * t301;
t210 = -t303 * Ifges(6,2) + t297;
t471 = Ifges(7,6) * t301;
t540 = t564 * t303 + t210 + t297 - t471;
t207 = t303 * Ifges(7,3) + t471;
t298 = Ifges(5,4) * t301;
t212 = t303 * Ifges(5,1) - t298;
t539 = -Ifges(5,2) * t303 + t207 + t212 - t298;
t538 = t277 * t543 + t573 * t279;
t228 = -t279 * mrSges(6,1) - mrSges(6,2) * t348;
t461 = t279 * mrSges(5,3);
t231 = -mrSges(5,1) * t348 + t461;
t537 = t228 / 0.2e1 - t231 / 0.2e1;
t336 = t348 * mrSges(6,3);
t469 = t277 * mrSges(6,1);
t226 = t336 + t469;
t467 = t277 * mrSges(5,3);
t229 = mrSges(5,2) * t348 - t467;
t536 = -t229 / 0.2e1 + t226 / 0.2e1;
t208 = t301 * Ifges(6,3) - t473;
t296 = Ifges(7,6) * t303;
t209 = t301 * Ifges(7,2) + t296;
t535 = t571 * t301 + t208 + t209 + t296 - t477;
t391 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1 - Ifges(5,6) / 0.2e1;
t392 = -Ifges(5,5) / 0.2e1 - Ifges(7,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t534 = -t391 * t278 - t392 * t428;
t343 = t345 ^ 2;
t533 = -0.2e1 * t279;
t532 = 2 * qJD(4);
t531 = m(4) / 0.2e1;
t530 = m(5) / 0.2e1;
t525 = m(4) * pkin(7);
t523 = -mrSges(7,1) / 0.2e1;
t522 = -mrSges(5,2) / 0.2e1;
t521 = mrSges(7,2) / 0.2e1;
t199 = t347 * pkin(3) - pkin(8) * t430 + t258;
t432 = t344 * t348;
t218 = -pkin(8) * t432 + t259;
t388 = -t199 * t490 + t489 * t218;
t50 = -pkin(5) * t428 - t346 * t347 + t388;
t520 = -t50 / 0.2e1;
t519 = t63 / 0.2e1;
t518 = -t64 / 0.2e1;
t517 = -t93 / 0.2e1;
t516 = -t123 / 0.2e1;
t508 = t277 / 0.2e1;
t501 = t301 / 0.2e1;
t499 = -t303 / 0.2e1;
t497 = t303 / 0.2e1;
t494 = -t344 / 0.2e1;
t493 = t345 / 0.2e1;
t487 = m(6) * t205;
t340 = t348 * pkin(7);
t482 = m(7) * qJD(5);
t481 = m(7) * qJD(6);
t480 = Ifges(4,4) * t344;
t479 = Ifges(4,4) * t345;
t476 = Ifges(4,5) * t345;
t475 = Ifges(4,6) * t344;
t239 = -pkin(7) * t432 + t345 * t312;
t306 = t348 * mrSges(4,2) - mrSges(4,3) * t433;
t308 = -t348 * mrSges(4,1) - mrSges(4,3) * t431;
t425 = t228 - t231;
t426 = t226 - t229;
t8 = -(t225 + t425) * t279 + (-t227 + t426) * t277 + m(7) * (-t277 * t56 - t279 * t48) + m(6) * (t277 * t78 - t279 * t81) + m(5) * (-t277 * t93 - t279 * t92) + (-t306 * t344 - t308 * t345 + m(4) * (-t239 * t345 - t240 * t344)) * t347;
t470 = qJD(1) * t8;
t466 = t278 * mrSges(5,1);
t465 = t278 * mrSges(7,1);
t464 = t278 * mrSges(6,2);
t463 = t278 * mrSges(7,3);
t460 = t428 * mrSges(6,1);
t459 = t428 * mrSges(7,1);
t458 = t428 * mrSges(5,2);
t457 = t428 * mrSges(7,2);
t456 = t428 * mrSges(6,3);
t311 = pkin(3) * t432 + t340;
t389 = qJ(5) * t428 + t311;
t122 = pkin(4) * t278 + t389;
t171 = t457 + t463;
t173 = -t458 + t466;
t174 = t456 - t464;
t446 = t347 * mrSges(7,3);
t221 = -t446 - t459;
t222 = mrSges(6,1) * t278 - mrSges(6,3) * t347;
t447 = t347 * mrSges(7,2);
t223 = t447 - t465;
t448 = t347 * mrSges(6,2);
t224 = t448 - t460;
t230 = -mrSges(5,2) * t347 - mrSges(5,3) * t278;
t232 = mrSges(5,1) * t347 + mrSges(5,3) * t428;
t450 = t344 * Ifges(4,2);
t272 = Ifges(4,6) * t347 + (-t450 + t479) * t348;
t273 = Ifges(4,5) * t347 + (t345 * Ifges(4,1) - t480) * t348;
t449 = t345 * mrSges(4,2);
t451 = t344 * mrSges(4,1);
t292 = (t449 + t451) * t348;
t307 = -t347 * mrSges(4,2) - mrSges(4,3) * t432;
t309 = t347 * mrSges(4,1) - mrSges(4,3) * t430;
t362 = t311 * mrSges(5,2) + Ifges(6,4) * t554 + Ifges(7,6) * t557 - t568 * t278 / 0.2e1 + t566 * t553 + (Ifges(7,3) + t565) * t555;
t363 = t311 * mrSges(5,1) + Ifges(5,6) * t554 + Ifges(7,6) * t555 + t568 * t556 + t567 * t553 + (Ifges(5,2) + t564) * t557;
t385 = -t161 / 0.2e1 + t156 / 0.2e1 + t158 / 0.2e1;
t386 = -t160 / 0.2e1 + t154 / 0.2e1 + t163 / 0.2e1;
t95 = t489 * t199 + t490 * t218;
t84 = -qJ(5) * t347 - t95;
t58 = -pkin(5) * t278 - t84;
t83 = t278 * t346 + t389;
t85 = -t347 * pkin(4) + t388;
t3 = (-pkin(1) * mrSges(3,2) + (t343 * Ifges(4,1) / 0.2e1 - Ifges(4,3) - Ifges(6,1) - Ifges(7,1) - Ifges(5,3) + Ifges(3,1) - Ifges(3,2) + (t449 + t525) * pkin(7) + (pkin(7) * mrSges(4,1) - t479 + t450 / 0.2e1) * t344) * t347 + (Ifges(3,4) + t475 - t476) * t348 + t534) * t348 + (-pkin(1) * mrSges(3,1) + t272 * t494 + t273 * t493 + pkin(7) * t292 + (-Ifges(3,4) + t476 / 0.2e1 - t475 / 0.2e1) * t347 + t392 * t279 + t391 * t277) * t347 + m(5) * (t310 * t311 + t388 * t92 + t93 * t95) - t388 * t231 - t386 * t428 + m(6) * (t121 * t122 + t78 * t84 + t81 * t85) + m(7) * (t48 * t50 + t56 * t58 + t82 * t83) + t259 * t306 + t240 * t307 + t258 * t308 + t239 * t309 + t310 * t173 + t84 * t226 + t58 * t227 + t85 * t228 + t95 * t229 + t93 * t230 - t92 * t232 + t48 * t221 + t78 * t222 + t56 * t223 + t81 * t224 + t50 * t225 + m(4) * (t239 * t258 + t240 * t259) + t82 * t171 + t122 * t172 + t121 * t174 + t83 * t169 - t362 * t279 + t363 * t277 + t385 * t278;
t455 = t3 * qJD(1);
t452 = t303 * mrSges(7,1);
t260 = t277 * mrSges(6,3);
t262 = t279 * mrSges(7,3);
t413 = Ifges(5,4) / 0.2e1 + Ifges(6,6) / 0.2e1;
t4 = -t310 * t261 + t121 * t260 - t82 * t262 + t63 * t227 + t64 * t225 + t170 * t172 + t99 * t169 + t425 * t93 + t426 * t92 + m(6) * (t121 * t170 + t78 * t92 + t81 * t93) + m(7) * (t48 * t64 + t56 * t63 + t82 * t99) + (t271 / 0.2e1 + t264 / 0.2e1 + t82 * mrSges(7,2) - t48 * mrSges(7,1) - t92 * mrSges(5,3) - t81 * mrSges(6,1) - t472 / 0.2e1 - t386) * t277 - (t310 * mrSges(5,1) - t263 / 0.2e1 - t121 * mrSges(6,2) - t93 * mrSges(5,3) - t56 * mrSges(7,1) + t78 * mrSges(6,1) + t413 * t279 + (-Ifges(5,1) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1 - Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t277 + t385) * t279 - t538 * t348 / 0.2e1;
t444 = t4 * qJD(1);
t13 = (t226 - t227) * t348 - (-t169 - t172) * t279 + m(7) * (t82 * t279 - t348 * t56) + m(6) * (t121 * t279 + t348 * t78);
t440 = t13 * qJD(1);
t439 = t168 * t279;
t22 = t348 * t225 + t277 * t169 + m(7) * (t82 * t277 + t348 * t48);
t437 = t22 * qJD(1);
t420 = -mrSges(6,1) + t523;
t419 = t527 + t529;
t418 = -m(7) / 0.4e1 - m(6) / 0.4e1;
t416 = m(7) * t519;
t415 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t414 = mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t412 = -t468 / 0.2e1;
t411 = mrSges(7,1) * t572;
t405 = t348 * t497;
t398 = t277 * mrSges(7,2) - t262;
t393 = mrSges(5,3) / 0.2e1 + t415;
t351 = (-t222 / 0.2e1 + t223 / 0.2e1) * qJ(5) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t347 + (-pkin(4) * t85 - qJ(5) * t84) * t528 + (qJ(5) * t58 - t346 * t50) * t526 - pkin(4) * t224 / 0.2e1 - t346 * t221 / 0.2e1 + mrSges(7,3) * t520 + t58 * t521 - t84 * mrSges(6,3) / 0.2e1 + t85 * mrSges(6,2) / 0.2e1 + t388 * t524 + t95 * t522 - t534;
t365 = -Ifges(5,2) / 0.4e1 - Ifges(7,2) / 0.4e1 - Ifges(6,3) / 0.4e1 + Ifges(7,3) / 0.4e1;
t2 = -t184 * t260 / 0.2e1 + t351 - t262 * t516 + ((-t78 / 0.2e1 + t517) * t303 + t570) * mrSges(6,1) + (-t439 / 0.2e1 + t166 * t508 + (t518 + t56 / 0.2e1) * t303 + (t48 / 0.2e1 + t519) * t301) * mrSges(7,1) - (t328 * t524 - t296 / 0.4e1 + t211 / 0.4e1 - t209 / 0.4e1 - t208 / 0.4e1 + t563 + t413 * t303 + t365 * t301) * t279 + (-t461 / 0.2e1 - t537) * t220 + (t467 / 0.2e1 - t536) * t219 + (-t298 / 0.4e1 + t471 / 0.2e1 - t297 / 0.4e1 + t207 / 0.4e1 + mrSges(7,2) * t516 - t210 / 0.4e1 + t212 / 0.4e1 + t365 * t303) * t277 - t263 * t498 + (t271 + t264) * t502 + t565 * (t277 * t496 - t279 * t500) + t561;
t5 = t328 * t400 + t123 * t397 + t205 * t206 + t546 * t499 + (t399 + t487) * t184 + t539 * t572 + t540 * t501 + t535 * t497 + t401 * t152;
t381 = -t2 * qJD(1) + t5 * qJD(2);
t16 = m(7) * (t166 * t303 - t168 * t301) + (m(4) * qJ(3) + mrSges(4,3)) * (t344 ^ 2 + t343) + (t301 ^ 2 + t303 ^ 2) * (mrSges(5,3) + t552) + (m(6) + m(5)) * (t219 * t303 - t220 * t301);
t375 = -t219 * t279 - t220 * t277;
t350 = (t277 * t393 + t510 + t536) * t301 + (t225 / 0.2e1 - t393 * t279 + t537) * t303 + (-t239 * t344 + t240 * t345) * t531 + (-t301 * t93 + t303 * t92 + t375) * t530 + (t301 * t78 + t303 * t81 + t375) * t528 + (-t166 * t279 - t168 * t277 - t301 * t56 + t303 * t48) * t526 + t308 * t494 + t306 * t493;
t359 = -m(5) * t311 / 0.2e1 + t122 * t529 + t83 * t527;
t7 = (-t449 / 0.2e1 - t451 / 0.2e1 - t525 / 0.2e1) * t348 - (t522 + t521 + mrSges(6,3) / 0.2e1) * t428 + (t524 - t414) * t278 + t350 + t359;
t380 = qJD(1) * t7 + qJD(2) * t16;
t30 = (-m(6) * t184 - t206 - t401) * t303;
t352 = (t515 + t514) * t303 - (t512 + t511) * t279 + (-t303 * t121 + t184 * t279 - t220 * t348) * t528 + (t123 * t279 - t168 * t348 - t303 * t82) * t526;
t372 = m(7) * t520 + t529 * t85;
t9 = t347 * t414 + 0.2e1 * t552 * t556 + t352 + t372;
t379 = -qJD(1) * t9 - qJD(2) * t30;
t354 = (t123 * t277 + t166 * t348 + t301 * t82) * t526 + t204 * t508 + t169 * t501;
t366 = t58 * t527 - t447 / 0.2e1;
t15 = (t405 + t557) * mrSges(7,1) + t354 + t366;
t38 = t401 * t301;
t378 = qJD(1) * t15 + qJD(2) * t38;
t21 = -0.2e1 * t170 * t528 - t279 * t485 + t574 * t99 - t260 + t261 - t398;
t25 = t152 * t574 + t303 * t485 - t293 + t294 - t397 - t487;
t377 = qJD(1) * t21 + qJD(2) * t25;
t60 = t419 * t533;
t88 = t419 * t560;
t376 = qJD(1) * t60 - qJD(2) * t88;
t335 = -0.2e1 * t429;
t355 = -t336 - t337 + (t335 + t93) * t528 + (t335 + t64) * t526;
t371 = m(6) * t517 + m(7) * t518;
t19 = t355 + t371;
t321 = qJ(5) * t559 + mrSges(7,2) + mrSges(6,3);
t373 = qJD(1) * t19 + qJD(4) * t321;
t149 = m(7) * t277;
t186 = m(7) * t301;
t370 = qJD(1) * t149 + qJD(2) * t186;
t353 = t445 + ((-qJ(6) - t346) * t348 + t360) * t527;
t29 = t416 + t353;
t318 = m(7) * t346 + mrSges(7,3);
t361 = -qJD(1) * t29 + qJD(4) * t318;
t315 = (qJD(1) * t348 - qJD(4)) * m(7);
t89 = t418 * t560 + t499 * t559;
t59 = t418 * t533 - t559 * t279 / 0.2e1;
t41 = -t452 + (-t526 + t527) * t166;
t36 = m(6) * t220 + m(7) * t168 + t301 * t420 + t411;
t26 = -t353 + t416 + t462;
t17 = t277 * t420 + t355 - t371 + t412;
t14 = mrSges(7,1) * t405 - t465 / 0.2e1 + t354 - t366;
t10 = -t446 / 0.2e1 - t459 / 0.2e1 + t448 / 0.2e1 - t460 / 0.2e1 + t415 * t428 + t352 - t372;
t6 = t350 - t458 / 0.2e1 + t457 / 0.2e1 + t456 / 0.2e1 - t464 / 0.2e1 + t463 / 0.2e1 + t466 / 0.2e1 + mrSges(4,2) * t430 / 0.2e1 + mrSges(4,1) * t432 / 0.2e1 + t340 * t531 - t359;
t1 = t166 * t412 + t351 + t123 * t398 / 0.2e1 + t260 * t513 - t439 * t523 - t474 * t498 - t271 * t502 + (t264 - t472) * t500 + (-t263 + t478) * t496 + (t64 / 0.2e1 - t56 / 0.2e1) * t452 + t550 * t411 + t537 * t220 + t536 * t219 + t545 * mrSges(5,3) - (mrSges(5,1) * t495 - t563 - Ifges(5,2) * t502 - t546 / 0.4e1 + t564 * t500 + t535 / 0.4e1) * t279 + (Ifges(6,2) * t498 + t571 * t496 - t539 / 0.4e1 + t540 / 0.4e1) * t277 + (t497 * t549 - t570) * mrSges(6,1) - t561;
t11 = [qJD(2) * t3 + qJD(3) * t8 + qJD(4) * t4 + qJD(5) * t13 + qJD(6) * t22, t6 * qJD(3) + t1 * qJD(4) + t10 * qJD(5) + t14 * qJD(6) + t455 + (t544 * mrSges(4,3) + 0.2e1 * (-pkin(2) * t340 + qJ(3) * t544) * t531 + (t85 * mrSges(6,1) + t50 * mrSges(7,1) + mrSges(5,3) * t388 - t347 * t392 + t362) * t303 + 0.2e1 * (t219 * t388 + t220 * t95 + t311 * t328) * t530 - (-t210 / 0.2e1 + t212 / 0.2e1 + t207 / 0.2e1) * t428 + mrSges(3,2) * t339 + t272 * t493 + (Ifges(4,5) * t344 + Ifges(4,6) * t345) * t553 + (t84 * mrSges(6,1) - t58 * mrSges(7,1) - t95 * mrSges(5,3) + t347 * t391 + t363) * t301 + (-t211 / 0.2e1 + t208 / 0.2e1 + t209 / 0.2e1) * t278 + (Ifges(3,5) + (Ifges(4,1) * t344 + t479) * t493 + (Ifges(4,2) * t345 + t480) * t494 + (-mrSges(4,1) * t345 + mrSges(4,2) * t344 - mrSges(3,1)) * pkin(7)) * t348 - Ifges(3,6) * t347 + t344 * t273 / 0.2e1 + t328 * t173 - pkin(2) * t292 + t220 * t230 - t219 * t232 + t166 * t221 - t220 * t222 + t168 * t223 + t219 * t224 + t83 * t204 + t122 * t206 + t184 * t174 + t123 * t171 - t344 * qJ(3) * t309 + t345 * qJ(3) * t307 + 0.2e1 * (t123 * t83 + t166 * t50 + t168 * t58) * t526 + 0.2e1 * (t122 * t184 + t219 * t85 - t220 * t84) * t528) * qJD(2), qJD(2) * t6 + qJD(5) * t59 + t470, t444 + t1 * qJD(2) + (t63 * mrSges(7,2) - t64 * mrSges(7,3) + pkin(4) * t469 + t346 * t468 + t442 * t552 + t485 * t93 + t551 * t92 + t538) * qJD(4) + t17 * qJD(5) + t26 * qJD(6) + ((qJ(5) * t63 - t346 * t64) * t526 + (-pkin(4) * t93 - qJ(5) * t92) * t528) * t532, qJD(2) * t10 + qJD(3) * t59 + qJD(4) * t17 + t440, qJD(2) * t14 + qJD(4) * t26 + t437; qJD(3) * t7 - qJD(4) * t2 + qJD(5) * t9 + qJD(6) * t15 - t455, qJD(3) * t16 + qJD(4) * t5 + qJD(5) * t30 + qJD(6) * t38, qJD(5) * t89 + t380, t36 * qJD(5) + t41 * qJD(6) + ((-pkin(4) * t220 - qJ(5) * t219) * t528 + (-qJ(5) * t166 - t168 * t346) * t526) * t532 + t381 + (-t166 * mrSges(7,2) - t168 * mrSges(7,3) + (-qJ(5) * t552 - t573) * t303 + (mrSges(6,1) * pkin(4) + mrSges(7,1) * t346 + t543) * t301 + t485 * t220 + t551 * t219) * qJD(4), qJD(3) * t89 + qJD(4) * t36 - t379, qJD(4) * t41 + t378; -qJD(2) * t7 - qJD(4) * t21 + qJD(5) * t60 + qJD(6) * t149 - t470, -qJD(4) * t25 - qJD(5) * t88 + qJD(6) * t186 - t380, 0, -t377, t376, t370; qJD(2) * t2 + qJD(3) * t21 + qJD(5) * t19 - qJD(6) * t29 - t444, qJD(3) * t25 - t381, t377, qJD(5) * t321 + qJD(6) * t318, t373, t361; -t9 * qJD(2) - t60 * qJD(3) - t19 * qJD(4) + t348 * t481 - t440, qJD(3) * t88 + t379, -t376, -t373 - t481, 0, t315; -t15 * qJD(2) - t149 * qJD(3) + t29 * qJD(4) - t348 * t482 - t437, -qJD(3) * t186 - t378, -t370, -t361 + t482, -t315, 0;];
Cq  = t11;
