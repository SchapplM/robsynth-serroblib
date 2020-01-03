% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:55
% EndTime: 2019-12-31 20:24:14
% DurationCPUTime: 8.60s
% Computational Cost: add. (19808->633), mult. (51082->910), div. (0->0), fcn. (55079->10), ass. (0->340)
t368 = sin(qJ(2));
t476 = sin(pkin(5));
t522 = cos(qJ(2));
t415 = t522 * t476;
t478 = cos(pkin(5));
t438 = pkin(1) * t478;
t309 = pkin(7) * t415 + t368 * t438;
t352 = t522 * t438;
t421 = t368 * t476;
t577 = -(-pkin(7) * t421 + t352) * mrSges(3,2) + Ifges(3,5) * t415 - t309 * mrSges(3,1);
t366 = sin(qJ(5));
t532 = t366 / 0.2e1;
t369 = cos(qJ(5));
t529 = -t369 / 0.2e1;
t365 = sin(pkin(10));
t477 = cos(pkin(10));
t413 = t477 * t476;
t294 = t365 * t415 + t368 * t413;
t367 = sin(qJ(4));
t370 = cos(qJ(4));
t252 = t294 * t370 + t367 * t478;
t493 = t252 * mrSges(5,3);
t293 = t365 * t421 - t413 * t522;
t196 = -t252 * t366 + t293 * t369;
t197 = t252 * t369 + t293 * t366;
t104 = -mrSges(6,1) * t196 + mrSges(6,2) * t197;
t207 = mrSges(5,1) * t293 - t493;
t576 = -t104 + t207;
t358 = Ifges(6,5) * t369;
t502 = Ifges(6,6) * t366;
t575 = t358 - t502;
t360 = Ifges(6,4) * t369;
t574 = -Ifges(6,2) * t366 + t360;
t336 = Ifges(6,1) * t366 + t360;
t251 = t294 * t367 - t370 * t478;
t130 = -mrSges(6,2) * t251 + t196 * mrSges(6,3);
t131 = mrSges(6,1) * t251 - mrSges(6,3) * t197;
t534 = -t366 / 0.2e1;
t393 = t130 * t534 + t131 * t529;
t518 = t367 * pkin(4);
t343 = -pkin(9) * t370 + t518;
t521 = pkin(2) * t365;
t355 = pkin(8) + t521;
t459 = t366 * t367;
t262 = t343 * t369 + t355 * t459;
t457 = t367 * t369;
t263 = t343 * t366 - t355 * t457;
t573 = -t262 * t366 + t263 * t369;
t275 = qJ(3) * t415 + t309;
t259 = t365 * t275;
t274 = t352 + (-pkin(7) - qJ(3)) * t421;
t204 = t274 * t477 - t259;
t416 = pkin(2) * t421;
t217 = pkin(3) * t294 + pkin(8) * t293 + t416;
t110 = -t204 * t367 + t217 * t370;
t111 = t204 * t370 + t217 * t367;
t572 = -t110 * t367 + t111 * t370;
t361 = Ifges(5,4) * t370;
t571 = Ifges(5,2) * t367 - t361;
t570 = t196 * t532 + t197 * t529;
t569 = -mrSges(6,3) * t570 - t393;
t328 = -mrSges(6,1) * t369 + mrSges(6,2) * t366;
t568 = -m(6) * pkin(4) - mrSges(5,1) + t328;
t567 = t367 ^ 2;
t566 = m(5) / 0.2e1;
t565 = -m(6) / 0.2e1;
t564 = m(6) / 0.2e1;
t563 = pkin(4) / 0.2e1;
t562 = -pkin(9) / 0.2e1;
t561 = m(4) * pkin(2);
t560 = -mrSges(6,1) / 0.2e1;
t559 = -mrSges(6,3) / 0.2e1;
t558 = mrSges(6,3) / 0.2e1;
t256 = pkin(2) * t478 + t274;
t420 = t477 * t275;
t192 = t256 * t365 + t420;
t185 = pkin(8) * t478 + t192;
t323 = -pkin(1) * t476 - pkin(2) * t415;
t200 = pkin(3) * t293 - pkin(8) * t294 + t323;
t91 = -t185 * t367 + t200 * t370;
t72 = -pkin(4) * t293 - t91;
t557 = -t72 / 0.2e1;
t556 = t72 / 0.2e1;
t92 = t185 * t370 + t200 * t367;
t555 = t92 / 0.2e1;
t554 = mrSges(6,3) * pkin(9);
t243 = Ifges(5,4) * t251;
t509 = Ifges(5,5) * t293;
t138 = Ifges(5,1) * t252 - t243 + t509;
t553 = -t138 / 0.2e1;
t552 = t196 / 0.2e1;
t551 = t197 / 0.2e1;
t550 = t197 / 0.4e1;
t458 = t366 * t370;
t218 = t293 * t458 + t294 * t369;
t549 = t218 / 0.2e1;
t452 = t369 * t370;
t219 = -t293 * t452 + t294 * t366;
t548 = t219 / 0.2e1;
t437 = t477 * pkin(2);
t356 = -t437 - pkin(3);
t316 = -pkin(4) * t370 - pkin(9) * t367 + t356;
t247 = t316 * t369 - t355 * t458;
t547 = t247 / 0.2e1;
t546 = t251 / 0.4e1;
t545 = t252 / 0.2e1;
t544 = t263 / 0.2e1;
t543 = t293 / 0.2e1;
t542 = t294 / 0.2e1;
t324 = mrSges(6,2) * t370 - mrSges(6,3) * t459;
t541 = t324 / 0.2e1;
t326 = -mrSges(6,1) * t370 - mrSges(6,3) * t457;
t540 = -t326 / 0.2e1;
t539 = t326 / 0.2e1;
t538 = -t328 / 0.2e1;
t332 = Ifges(6,5) * t366 + Ifges(6,6) * t369;
t537 = t332 / 0.2e1;
t536 = -t336 / 0.4e1;
t535 = t355 / 0.2e1;
t533 = -t366 / 0.4e1;
t531 = -t367 / 0.2e1;
t530 = t367 / 0.2e1;
t528 = -t369 / 0.4e1;
t527 = t369 / 0.2e1;
t526 = t369 / 0.4e1;
t525 = -t370 / 0.2e1;
t524 = -t370 / 0.4e1;
t523 = t370 / 0.2e1;
t520 = pkin(9) * t366;
t73 = pkin(9) * t293 + t92;
t191 = t256 * t477 - t259;
t184 = -pkin(3) * t478 - t191;
t88 = pkin(4) * t251 - pkin(9) * t252 + t184;
t36 = -t366 * t73 + t369 * t88;
t519 = t36 * mrSges(6,3);
t37 = t366 * t88 + t369 * t73;
t517 = t37 * mrSges(6,3);
t248 = t316 * t366 + t355 * t452;
t485 = t369 * mrSges(6,2);
t487 = t366 * mrSges(6,1);
t330 = t485 + t487;
t311 = t330 * t367;
t312 = t370 * t330;
t325 = -mrSges(6,2) * t367 - mrSges(6,3) * t458;
t327 = mrSges(6,1) * t367 - mrSges(6,3) * t452;
t430 = t457 / 0.2e1;
t432 = -t459 / 0.2e1;
t456 = t567 * t355;
t364 = t370 ^ 2;
t462 = t364 * t355;
t39 = t311 * t530 + t327 * t432 + (-t262 * t459 + t263 * t457 + t456 - t462) * t564 + t325 * t430 + (-t312 / 0.2e1 + t326 * t534 + (-t247 * t366 + t248 * t369) * t564 + t324 * t527) * t370;
t310 = -mrSges(6,1) * t457 + mrSges(6,2) * t459;
t362 = t366 ^ 2;
t363 = t369 ^ 2;
t423 = t362 / 0.2e1 + t363 / 0.2e1;
t68 = t310 * t525 + (mrSges(6,3) * t367 * t423 + t324 * t532 + t326 * t527) * t367;
t516 = qJD(4) * t39 - qJD(5) * t68;
t515 = m(6) * qJD(3);
t514 = mrSges(6,3) * t251;
t512 = Ifges(5,4) * t252;
t511 = Ifges(5,4) * t367;
t510 = Ifges(6,4) * t366;
t508 = Ifges(6,5) * t219;
t507 = Ifges(6,5) * t370;
t504 = Ifges(5,6) * t293;
t503 = Ifges(6,6) * t218;
t501 = Ifges(6,3) * t252;
t500 = Ifges(6,3) * t367;
t499 = t197 * Ifges(6,4);
t498 = t218 * mrSges(6,1);
t497 = t219 * mrSges(6,2);
t496 = t247 * mrSges(6,3);
t495 = t248 * mrSges(6,3);
t494 = t251 * mrSges(5,3);
t492 = t293 * mrSges(5,2);
t491 = t293 * mrSges(4,3);
t490 = t294 * mrSges(4,3);
t117 = t497 - t498;
t466 = t293 * t367;
t150 = mrSges(6,2) * t466 + mrSges(6,3) * t218;
t151 = -mrSges(6,1) * t466 - mrSges(6,3) * t219;
t159 = Ifges(5,6) * t294 + t293 * t571;
t338 = Ifges(5,1) * t370 - t511;
t160 = Ifges(5,5) * t294 - t293 * t338;
t203 = t274 * t365 + t420;
t206 = -t492 - t494;
t331 = t367 * mrSges(5,1) + t370 * mrSges(5,2);
t220 = t331 * t293;
t222 = -mrSges(5,2) * t294 + mrSges(5,3) * t466;
t465 = t293 * t370;
t223 = mrSges(5,1) * t294 + mrSges(5,3) * t465;
t257 = -mrSges(4,2) * t478 - t491;
t277 = t293 * mrSges(4,2);
t446 = t476 ^ 2;
t414 = t446 * t522;
t422 = Ifges(3,4) * t446;
t137 = -Ifges(5,2) * t251 + t504 + t512;
t69 = Ifges(6,5) * t197 + Ifges(6,6) * t196 + Ifges(6,3) * t251;
t439 = t137 / 0.2e1 - t69 / 0.2e1;
t440 = -t478 / 0.2e1;
t448 = -Ifges(4,5) * t293 - Ifges(4,6) * t294;
t450 = -mrSges(4,1) * t478 + mrSges(5,1) * t251 + mrSges(5,2) * t252 + t490;
t118 = -t293 * t343 + t203;
t78 = pkin(9) * t294 + t111;
t46 = t118 * t369 - t366 * t78;
t47 = t118 * t366 + t369 * t78;
t70 = Ifges(6,2) * t196 + Ifges(6,6) * t251 + t499;
t190 = Ifges(6,4) * t196;
t71 = Ifges(6,1) * t197 + Ifges(6,5) * t251 + t190;
t77 = -pkin(4) * t294 - t110;
t85 = -Ifges(6,3) * t466 + t503 + t508;
t86 = Ifges(6,4) * t219 + Ifges(6,2) * t218 - Ifges(6,6) * t466;
t87 = Ifges(6,1) * t219 + Ifges(6,4) * t218 - Ifges(6,5) * t466;
t3 = -t323 * t277 + (t448 / 0.2e1 + t577) * t478 + (t85 / 0.2e1 - t159 / 0.2e1) * t251 + (t191 * mrSges(4,3) + Ifges(4,4) * t293 + Ifges(4,5) * t440 + (t553 - t509 / 0.2e1) * t370 + (t504 / 0.2e1 + t439) * t367 + (Ifges(5,3) + Ifges(4,2) - Ifges(4,1)) * t294) * t293 + (Ifges(5,5) * t545 - Ifges(5,6) * t251 / 0.2e1 - t192 * mrSges(4,3) + t323 * mrSges(4,1) + Ifges(4,6) * t440 - Ifges(4,4) * t294) * t294 + t522 ^ 2 * t422 + t450 * t203 + t204 * t257 - t184 * t220 + t92 * t222 + t91 * t223 + t111 * t206 + t110 * t207 + m(4) * (-t191 * t203 + t192 * t204) + t37 * t150 + t36 * t151 + t47 * t130 + t46 * t131 + t72 * t117 + t77 * t104 - pkin(1) * t414 * mrSges(3,2) + (-t446 * pkin(1) * mrSges(3,1) + (-t478 * Ifges(3,6) + (m(4) * t323 + mrSges(4,1) * t293 + mrSges(4,2) * t294) * pkin(2)) * t476 + (Ifges(3,1) - Ifges(3,2)) * t414 - t422 * t368) * t368 + m(5) * (t110 * t91 + t111 * t92 + t184 * t203) + m(6) * (t36 * t46 + t37 * t47 + t72 * t77) + t160 * t545 + t71 * t548 + t70 * t549 + t87 * t551 + t86 * t552;
t489 = t3 * qJD(1);
t486 = t366 * t70;
t484 = t369 * t71;
t483 = t370 * Ifges(6,6);
t119 = -t251 * t575 + t501;
t120 = Ifges(6,6) * t252 - t251 * t574;
t412 = Ifges(6,1) * t369 - t510;
t121 = Ifges(6,5) * t252 - t251 * t412;
t161 = t330 * t251;
t168 = -mrSges(6,2) * t252 + t366 * t514;
t169 = mrSges(6,1) * t252 + t369 * t514;
t175 = mrSges(5,1) * t252 - mrSges(5,2) * t251;
t177 = -Ifges(5,2) * t252 - t243;
t178 = -Ifges(5,1) * t251 - t512;
t449 = -Ifges(5,5) * t251 - Ifges(5,6) * t252;
t179 = pkin(4) * t252 + pkin(9) * t251;
t53 = t179 * t369 - t366 * t91;
t54 = t179 * t366 + t369 * t91;
t4 = t449 * t543 + t91 * t206 - t72 * t161 + t121 * t551 + t120 * t552 + t54 * t130 + t37 * t168 + t53 * t131 + t36 * t169 + m(6) * (t36 * t53 + t37 * t54) + t184 * t175 + (t178 / 0.2e1 - t439) * t252 + (t553 - t177 / 0.2e1 + t119 / 0.2e1 - t484 / 0.2e1 + t486 / 0.2e1 + t91 * mrSges(5,3)) * t251 + (m(6) * t72 - t493 - t576) * t92;
t482 = t4 * qJD(1);
t481 = t46 * t366;
t480 = t47 * t369;
t103 = mrSges(6,1) * t197 + mrSges(6,2) * t196;
t105 = Ifges(6,5) * t196 - Ifges(6,6) * t197;
t106 = -Ifges(6,2) * t197 + t190;
t107 = Ifges(6,1) * t196 - t499;
t7 = t72 * t103 + t251 * t105 / 0.2e1 + t36 * t130 - t37 * t131 + (t107 / 0.2e1 - t70 / 0.2e1 - t517) * t197 + (t71 / 0.2e1 + t106 / 0.2e1 - t519) * t196;
t479 = t7 * qJD(1);
t253 = t293 * t456;
t329 = -mrSges(5,1) * t370 + mrSges(5,2) * t367;
t436 = -t466 / 0.2e1;
t372 = (-t293 * t462 + t294 * t356 - t253) * t566 + (t218 * t247 + t219 * t248 - t253) * t564 + t218 * t539 + t219 * t541 + t329 * t542 + (-t293 * t365 - t294 * t477) * t561 / 0.2e1 + t311 * t436 - (t364 + t567) * mrSges(5,3) * t293 / 0.2e1;
t409 = t480 - t481;
t428 = t150 * t527;
t373 = (t110 * t370 + t111 * t367) * t566 + (t367 * t409 - t370 * t77) * t564 + t222 * t530 + t117 * t525 + t223 * t523 + t151 * t432 + t367 * t428 + m(4) * t416 / 0.2e1;
t13 = -t294 * mrSges(4,1) + t277 + t372 - t373;
t475 = qJD(1) * t13;
t451 = t370 * t206;
t12 = t219 * t130 + t218 * t131 + t450 * t294 + (t576 * t367 - t257 - t451) * t293 + m(6) * (t218 * t36 + t219 * t37 - t466 * t72) + m(5) * (t184 * t294 + (t367 * t91 - t370 * t92) * t293) + m(4) * (-t191 * t294 - t192 * t293);
t472 = t12 * qJD(1);
t470 = t218 * t366;
t469 = t219 * t369;
t463 = t355 * t370;
t461 = t366 * t169;
t297 = t367 * t574 - t483;
t460 = t366 * t297;
t455 = t369 * t168;
t299 = t367 * t412 - t507;
t453 = t369 * t299;
t447 = t362 + t363;
t445 = m(6) * t555;
t443 = -t358 / 0.2e1;
t442 = t492 / 0.2e1;
t441 = mrSges(5,3) * t531;
t435 = t466 / 0.2e1;
t434 = -t465 / 0.2e1;
t431 = -t458 / 0.2e1;
t427 = t452 / 0.2e1;
t426 = -t161 / 0.2e1 - t206 / 0.2e1;
t315 = t336 * t367;
t425 = -t297 / 0.4e1 - t315 / 0.4e1;
t333 = Ifges(6,2) * t369 + t510;
t314 = t333 * t367;
t424 = t299 / 0.4e1 - t314 / 0.4e1;
t417 = t447 * t370;
t408 = -t366 * t53 + t369 * t54;
t379 = t103 * t525 - t367 * t569;
t401 = t498 / 0.2e1 - t497 / 0.2e1;
t18 = t379 - t401;
t407 = -qJD(1) * t18 + qJD(2) * t68;
t406 = t469 - t470;
t404 = -t36 * t366 + t369 * t37 - t92;
t403 = t408 + t72;
t402 = t310 * t563 + t358 * t524;
t400 = mrSges(6,2) * t544 + t262 * t560;
t399 = -t486 / 0.4e1 + t484 / 0.4e1;
t398 = -t485 / 0.2e1 - t487 / 0.2e1;
t397 = t324 * t562 + t425;
t396 = pkin(9) * t540 + t424;
t395 = t104 / 0.2e1 - t207 / 0.2e1 - t493 / 0.2e1;
t394 = -t494 / 0.2e1 + t426;
t391 = t333 * t532 + t336 * t529;
t313 = t367 * t332;
t390 = t406 * pkin(9);
t295 = -Ifges(6,3) * t370 + t367 * t575;
t298 = Ifges(6,6) * t367 + t370 * t574;
t300 = Ifges(6,5) * t367 + t370 * t412;
t334 = Ifges(5,2) * t370 + t511;
t359 = Ifges(5,5) * t370;
t371 = (t338 / 0.4e1 - t334 / 0.4e1 + t295 / 0.4e1) * t252 + (t247 * t53 + t248 * t54 + t262 * t36 + t263 * t37) * t564 + t184 * t331 / 0.2e1 + t196 * t298 / 0.4e1 + t300 * t550 + t169 * t547 + t248 * t168 / 0.2e1 + t262 * t131 / 0.2e1 + t130 * t544 + t293 * t359 / 0.4e1 + t356 * t175 / 0.2e1 + t36 * t327 / 0.2e1 + t37 * t325 / 0.2e1 + t53 * t539 + t54 * t541 + t312 * t556 + t311 * t555;
t376 = (-pkin(4) * t77 + pkin(9) * t409) * t565 - Ifges(5,3) * t294 / 0.2e1 + t117 * t563 - t110 * mrSges(5,1) / 0.2e1 + t111 * mrSges(5,2) / 0.2e1 - t218 * t333 / 0.4e1 + t219 * t536 + t77 * t538;
t378 = t177 / 0.4e1 + t138 / 0.4e1 - t119 / 0.4e1 + (m(6) * t556 + t395) * t355 + t399;
t296 = t370 * t575 + t500;
t337 = Ifges(5,1) * t367 + t361;
t382 = -t453 / 0.4e1 + t460 / 0.4e1 + t571 / 0.4e1 - t337 / 0.4e1 + t296 / 0.4e1;
t383 = -t137 / 0.4e1 + t178 / 0.4e1 + t69 / 0.4e1 + t120 * t533 + t121 * t526;
t1 = t371 + t382 * t251 + ((-0.3e1 / 0.4e1 * Ifges(5,6) + t332 / 0.4e1) * t293 + (t445 + t394) * t355 + t383) * t367 + (t509 / 0.2e1 + t378) * t370 + (t150 * t562 + t47 * t559 - t86 / 0.4e1) * t369 + (t46 * t558 + pkin(9) * t151 / 0.2e1 - t87 / 0.4e1) * t366 + t376;
t24 = t262 * t326 + t247 * t327 + m(6) * (t247 * t262 + t248 * t263) + t263 * t324 + t248 * t325 + t356 * t331 + (t355 * t311 - t460 / 0.2e1 + t453 / 0.2e1 - t296 / 0.2e1 + t337 / 0.2e1 - t571 / 0.2e1) * t370 + (t298 * t534 + t300 * t527 + t295 / 0.2e1 + t338 / 0.2e1 - t334 / 0.2e1 + (m(6) * t463 + t312) * t355) * t367;
t389 = t1 * qJD(1) + t24 * qJD(2) + t39 * qJD(3);
t27 = -t313 * t525 - t248 * t326 + t247 * t324 + (-t355 * t310 + (-t297 / 0.2e1 - t315 / 0.2e1 - t495) * t369 + (t314 / 0.2e1 - t299 / 0.2e1 + t496) * t366) * t367;
t375 = (-t496 / 0.2e1 + t424) * t196 + (-t495 / 0.2e1 + t425) * t197 + t130 * t547 - t248 * t131 / 0.2e1 - t313 * t546 + t36 * t541 + t37 * t540 + t105 * t524 + t310 * t557;
t380 = t103 * t535 + (t107 / 0.4e1 - t70 / 0.4e1 - t517 / 0.2e1) * t369 + (-t106 / 0.4e1 - t71 / 0.4e1 + t519 / 0.2e1) * t366;
t381 = -t508 / 0.2e1 - t503 / 0.2e1 + t46 * t560 + t47 * mrSges(6,2) / 0.2e1;
t5 = (Ifges(6,3) * t543 + t380) * t367 + t375 + t381;
t388 = qJD(1) * t5 + qJD(2) * t27 - qJD(3) * t68;
t387 = (t469 / 0.2e1 - t470 / 0.2e1) * mrSges(6,3);
t10 = t390 * t564 + t387 + (t403 * t565 - t455 / 0.2e1 + t461 / 0.2e1 + (t538 + m(6) * t563 + mrSges(5,1) / 0.2e1) * t293 - t395) * t367 + (t130 * t529 + t131 * t532 + t404 * t565 + t394 + t442) * t370;
t271 = (-0.1e1 + t447) * t370 * t367;
t386 = -qJD(1) * t10 + qJD(2) * t39 + t271 * t515;
t385 = t501 / 0.2e1 + t53 * mrSges(6,1) / 0.2e1 - t54 * mrSges(6,2) / 0.2e1;
t170 = pkin(4) * t330 + t412 * t534 + t529 * t574 + t391;
t249 = (t330 / 0.2e1 + t398) * t370;
t377 = t330 * t535 + t333 * t528 + t412 * t526 - t423 * t554 + (t336 + t574) * t533;
t25 = (-t507 / 0.2e1 + t396) * t369 + (0.3e1 / 0.4e1 * t483 + t397) * t366 + (-Ifges(6,3) / 0.2e1 + t377) * t367 + t400 + t402;
t374 = t103 * t563 + t333 * t550 - t197 * t412 / 0.4e1 + t107 * t533 + t106 * t528 + t330 * t557 - t399 + (t536 - t574 / 0.4e1) * t196;
t8 = t374 + t569 * pkin(9) + (t443 + 0.3e1 / 0.4e1 * t502 - t358 / 0.4e1) * t251 + t385;
t384 = t8 * qJD(1) - t25 * qJD(2) + t249 * qJD(3) + t170 * qJD(4);
t250 = t330 * t525 + t370 * t398;
t26 = Ifges(6,6) * t431 + Ifges(6,5) * t427 + t500 / 0.2e1 + (t483 / 0.4e1 + t397) * t366 + t396 * t369 + t377 * t367 - t400 + t402;
t19 = t379 + t401;
t14 = t372 + t373;
t11 = t104 * t530 - t161 * t525 + t130 * t427 + t168 * t430 + t131 * t431 + t169 * t432 + t451 / 0.2e1 + t252 * t441 + t207 * t531 + t494 * t523 + t328 * t436 + t370 * t442 + mrSges(5,1) * t435 + t387 + (pkin(4) * t466 + t367 * t403 + t370 * t404 + t390) * t564;
t9 = -t374 + (t443 + t502 / 0.2e1) * t251 + t575 * t546 + t385 + t570 * t554 + t393 * pkin(9);
t6 = Ifges(6,3) * t436 + t367 * t380 + t375 - t381;
t2 = -t376 + t371 + (t355 * t441 + t382) * t251 + Ifges(5,5) * t434 + pkin(9) * t428 + t480 * t558 + Ifges(5,6) * t435 - t332 * t466 / 0.4e1 + t481 * t559 - t151 * t520 / 0.2e1 + t86 * t526 + t366 * t87 / 0.4e1 + (-t504 / 0.4e1 + (t445 + t426) * t355 + t383) * t367 + t378 * t370;
t15 = [qJD(2) * t3 + qJD(3) * t12 + qJD(4) * t4 + qJD(5) * t7, t489 + (t577 + (t365 * t561 - mrSges(4,2)) * t204 + (m(5) * t356 - t477 * t561 - mrSges(4,1) + t329) * t203 + t448 + t86 * t432 + t337 * t434 + t334 * t435 + t295 * t436 + t87 * t430 + t85 * t525 + t159 * t523 + t160 * t530 + m(6) * (t247 * t46 + t248 * t47) + t437 * t491 + t222 * t463 - t490 * t521 - t356 * t220 + t47 * t324 + t46 * t326 + t77 * t311 + t248 * t150 + t247 * t151 - Ifges(3,6) * t421 + (m(5) * t572 + (m(6) * t77 + t117 - t223) * t367) * t355 + t572 * mrSges(5,3) + (Ifges(5,5) * t367 + Ifges(5,6) * t370) * t542 + t299 * t548 + t297 * t549) * qJD(2) + t14 * qJD(3) + t2 * qJD(4) + t6 * qJD(5), t472 + t14 * qJD(2) + t11 * qJD(4) + t19 * qJD(5) + (t406 + t465) * t367 * t515, t482 + t2 * qJD(2) + t11 * qJD(3) + (-t91 * mrSges(5,2) + t408 * mrSges(6,3) + pkin(4) * t161 + t120 * t527 + t121 * t532 + t391 * t251 + t252 * t537 + t449 + t568 * t92 + (m(6) * t408 + t455 - t461) * pkin(9)) * qJD(4) + t9 * qJD(5), t479 + t6 * qJD(2) + t19 * qJD(3) + t9 * qJD(4) + (-mrSges(6,1) * t37 - mrSges(6,2) * t36 + t105) * qJD(5); qJD(3) * t13 + qJD(4) * t1 + qJD(5) * t5 - t489, qJD(4) * t24 + qJD(5) * t27, t475 + t516, t26 * qJD(5) + t389 + (t333 * t431 + t336 * t427 + t300 * t532 + t298 * t527 - pkin(4) * t312 - t327 * t520 + t359 + (mrSges(5,2) * t355 - Ifges(5,6) + t537) * t367 + (m(6) * t573 + t369 * t325) * pkin(9) + t568 * t463 + t573 * mrSges(6,3)) * qJD(4), t26 * qJD(4) + (-mrSges(6,1) * t248 - mrSges(6,2) * t247 - t313) * qJD(5) + t388; -qJD(2) * t13 - qJD(4) * t10 + qJD(5) * t18 - t472, -t475 + t516, m(6) * t271 * qJD(4), (t367 * t328 + m(6) * (pkin(9) * t417 - t518) + mrSges(6,3) * t417 - t331) * qJD(4) + t250 * qJD(5) + t386, qJD(4) * t250 + qJD(5) * t310 - t407; -qJD(2) * t1 + qJD(3) * t10 - qJD(5) * t8 - t482, qJD(5) * t25 - t389, -qJD(5) * t249 - t386, -t170 * qJD(5), (pkin(9) * t328 + t575) * qJD(5) - t384; -qJD(2) * t5 - qJD(3) * t18 + qJD(4) * t8 - t479, -qJD(4) * t25 - t388, qJD(4) * t249 + t407, t384, 0;];
Cq = t15;
