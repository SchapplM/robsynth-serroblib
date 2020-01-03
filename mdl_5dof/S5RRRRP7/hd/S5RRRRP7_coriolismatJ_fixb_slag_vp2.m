% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:11
% EndTime: 2019-12-31 21:56:28
% DurationCPUTime: 7.83s
% Computational Cost: add. (11674->506), mult. (23838->674), div. (0->0), fcn. (23853->6), ass. (0->283)
t335 = sin(qJ(4));
t333 = t335 ^ 2;
t336 = cos(qJ(4));
t334 = t336 ^ 2;
t594 = t334 + t333;
t503 = Ifges(6,5) * t336;
t308 = Ifges(6,1) * t335 - t503;
t332 = Ifges(5,4) * t336;
t310 = Ifges(5,1) * t335 + t332;
t593 = t310 + t308;
t528 = -t335 / 0.2e1;
t527 = t335 / 0.2e1;
t525 = -t336 / 0.2e1;
t521 = cos(qJ(3));
t592 = t594 * t521;
t303 = Ifges(6,3) * t335 + t503;
t329 = Ifges(6,5) * t335;
t309 = Ifges(6,1) * t336 + t329;
t386 = Ifges(5,2) * t335 - t332;
t500 = Ifges(6,3) * t336;
t302 = t329 - t500;
t504 = Ifges(5,4) * t335;
t306 = Ifges(5,2) * t336 + t504;
t560 = t302 * t528 + t306 * t527 + t593 * t525;
t591 = t309 * t527 + (t303 + t386) * t525 - t560;
t511 = Ifges(6,4) + Ifges(5,5);
t582 = Ifges(5,6) - Ifges(6,6);
t590 = t511 * t335 + t582 * t336;
t519 = sin(qJ(3));
t520 = sin(qJ(2));
t522 = cos(qJ(2));
t288 = t519 * t520 - t521 * t522;
t273 = t288 * mrSges(6,3);
t289 = -t519 * t522 - t520 * t521;
t466 = t289 * t335;
t213 = mrSges(6,2) * t466 + t273;
t429 = mrSges(5,3) * t466;
t369 = -mrSges(5,2) * t288 + t429;
t566 = t213 + t369;
t311 = Ifges(5,1) * t336 - t504;
t564 = -t311 - t309;
t465 = t289 * t336;
t370 = t288 * mrSges(6,1) + mrSges(6,2) * t465;
t580 = mrSges(5,3) * t465;
t371 = t288 * mrSges(5,1) + t580;
t563 = t371 + t370;
t297 = -t336 * mrSges(6,1) - t335 * mrSges(6,3);
t467 = t288 * t336;
t468 = t288 * t335;
t431 = t520 * pkin(6);
t316 = -pkin(7) * t520 - t431;
t434 = t522 * pkin(6);
t317 = pkin(7) * t522 + t434;
t551 = t519 * t316 + t521 * t317;
t554 = -pkin(4) * t468 + qJ(5) * t467 + t551;
t568 = t554 * t297;
t238 = -t521 * t316 + t519 * t317;
t579 = t238 * mrSges(4,2);
t588 = t568 / 0.2e1 + t579 / 0.2e1;
t298 = -t336 * mrSges(5,1) + t335 * mrSges(5,2);
t569 = t551 * t298;
t578 = t551 * mrSges(4,1);
t587 = t569 / 0.2e1 - t578 / 0.2e1;
t586 = t289 * t594;
t524 = t336 / 0.2e1;
t534 = -t289 / 0.2e1;
t352 = t590 * t534 + (-Ifges(5,6) * t524 - Ifges(6,6) * t525 - t511 * t527 + Ifges(4,6)) * t289 + (-t303 * t525 + t386 * t524 + t564 * t527 - Ifges(4,5) + t560) * t288;
t585 = t352 + t568 + t569 + t579 - t578;
t426 = mrSges(6,2) * t467;
t493 = t289 * mrSges(6,1);
t212 = -t426 + t493;
t584 = t212 / 0.2e1;
t583 = pkin(3) * t551;
t556 = mrSges(6,2) + mrSges(5,3);
t581 = Ifges(5,3) + Ifges(6,2);
t390 = -t522 * pkin(2) - pkin(1);
t205 = t288 * pkin(3) + t289 * pkin(8) + t390;
t97 = t205 * t336 - t335 * t551;
t74 = -t288 * pkin(4) - t97;
t508 = t74 + t97;
t299 = t335 * pkin(4) - qJ(5) * t336;
t127 = -t289 * t299 + t238;
t577 = t127 * t554;
t576 = t238 * t519;
t473 = t238 * t551;
t458 = t335 * qJ(5);
t385 = -t336 * pkin(4) - t458;
t295 = -pkin(3) + t385;
t433 = t521 * pkin(2);
t271 = -t433 + t295;
t575 = t271 * t554;
t574 = t295 * t554;
t325 = -t433 - pkin(3);
t573 = t325 * t551;
t572 = t335 * t238;
t571 = t336 * t238;
t570 = t336 * t508;
t265 = Ifges(6,5) * t465;
t203 = Ifges(6,1) * t466 - t265;
t204 = t310 * t289;
t567 = t203 + t204;
t330 = Ifges(5,5) * t336;
t501 = Ifges(5,6) * t335;
t304 = t330 - t501;
t328 = Ifges(6,6) * t335;
t331 = Ifges(6,4) * t336;
t305 = t331 + t328;
t565 = t305 + t304;
t562 = t592 * pkin(2);
t561 = -Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t507 = mrSges(6,3) * t336;
t300 = t335 * mrSges(6,1) - t507;
t198 = t300 * t288;
t301 = t335 * mrSges(5,1) + mrSges(5,2) * t336;
t199 = t301 * t288;
t200 = t300 * t289;
t210 = mrSges(5,2) * t289 + mrSges(5,3) * t468;
t211 = -t289 * mrSges(5,1) + mrSges(5,3) * t467;
t214 = mrSges(6,2) * t468 - mrSges(6,3) * t289;
t367 = t301 * t289;
t372 = -Ifges(4,4) + t330 + t331;
t469 = t288 * qJ(5);
t448 = t336 * t551;
t455 = t335 * t205;
t98 = t448 + t455;
t73 = t98 + t469;
t558 = t127 * t198 + t238 * t199 - t98 * t210 - t97 * t211 - t74 * t212 - t73 * t214 - t390 * (-mrSges(4,1) * t289 - mrSges(4,2) * t288) + t554 * t200 + (-t328 - t372 + t501) * t289 ^ 2 + t551 * t367;
t557 = mrSges(5,1) + mrSges(6,1);
t555 = t210 + t214;
t442 = t271 + t295;
t553 = t297 + t298;
t513 = t288 * pkin(8);
t215 = -t289 * pkin(3) + t513;
t432 = t520 * pkin(2);
t207 = t432 + t215;
t112 = t335 * t207 - t571;
t481 = t112 * t336;
t111 = t207 * t336 + t572;
t482 = t111 * t335;
t377 = t481 - t482;
t512 = t289 * pkin(4);
t77 = -t111 + t512;
t489 = t77 * t335;
t272 = t289 * qJ(5);
t76 = -t272 + t112;
t490 = t76 * t336;
t384 = t489 + t490;
t550 = -qJ(5) * t214 / 0.2e1 + pkin(4) * t584;
t326 = m(6) * qJ(5) + mrSges(6,3);
t430 = t519 * pkin(2);
t324 = t430 + pkin(8);
t461 = t324 * t336;
t549 = t555 * t461;
t548 = t563 * t525 + t566 * t528;
t197 = t385 * t289;
t547 = t299 * t200 / 0.2e1 - t238 * t301 / 0.2e1 - t197 * t297 / 0.2e1 - t127 * t300 / 0.2e1;
t449 = t335 * t311;
t279 = t449 / 0.2e1;
t546 = t279 + t591;
t464 = t295 * t198;
t517 = pkin(3) * t199;
t541 = m(6) / 0.2e1;
t543 = m(5) / 0.2e1;
t545 = (pkin(8) * t377 - t583) * t543 + (pkin(8) * t384 + t574) * t541 + t517 / 0.2e1 - t464 / 0.2e1 + (t481 / 0.2e1 - t482 / 0.2e1) * mrSges(5,3) + (t490 / 0.2e1 + t489 / 0.2e1) * mrSges(6,2) + t588;
t544 = 2 * qJD(3);
t542 = -m(6) / 0.2e1;
t540 = mrSges(5,1) / 0.2e1;
t539 = mrSges(6,1) / 0.2e1;
t538 = -mrSges(5,2) / 0.2e1;
t120 = t215 * t336 + t572;
t81 = -t120 + t512;
t537 = t81 / 0.2e1;
t536 = t271 / 0.2e1;
t533 = t295 / 0.2e1;
t531 = -t324 / 0.2e1;
t530 = t324 / 0.2e1;
t529 = t325 / 0.2e1;
t518 = m(6) * t299;
t516 = pkin(3) * t301;
t514 = pkin(8) * t336;
t509 = -t73 + t98;
t342 = (-t265 / 0.2e1 - t582 * t288 + ((-Ifges(5,2) - Ifges(6,3)) * t335 + (0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(6,5)) * t336) * t289) * t335 + (-Ifges(4,1) + Ifges(4,2) + (-Ifges(5,1) - Ifges(6,1)) * t334 + t581) * t289;
t494 = t288 * mrSges(4,3);
t2 = -t76 * t213 + (t238 * mrSges(4,3) + t288 * t372 + t342) * t288 + pkin(1) * (mrSges(3,1) * t520 + mrSges(3,2) * t522) - t238 * t494 - m(5) * (t111 * t97 + t112 * t98 + t473) - m(6) * (t73 * t76 + t74 * t77 + t577) - t112 * t369 - t111 * t371 + t77 * t370 + (t520 ^ 2 - t522 ^ 2) * Ifges(3,4) + (Ifges(3,2) - Ifges(3,1)) * t522 * t520 + t558 + (-m(4) * t390 - mrSges(4,1) * t288 + mrSges(4,2) * t289) * t432;
t499 = t2 * qJD(1);
t121 = t335 * t215 - t571;
t80 = -t272 + t121;
t3 = -t80 * t213 + ((t336 * t511 - Ifges(4,4)) * t288 + t342) * t288 - m(5) * (t120 * t97 + t121 * t98 + t473) - m(6) * (t73 * t80 + t74 * t81 + t577) - t121 * t369 - t120 * t371 + t81 * t370 + t558;
t491 = t3 * qJD(1);
t327 = t336 * mrSges(6,2);
t488 = t80 * t336;
t487 = t81 * t335;
t201 = t302 * t289;
t202 = t306 * t289;
t366 = t289 * t297;
t368 = t298 * t289;
t9 = -t238 * t368 + t197 * t200 + (t265 * t525 + t201 * t527 + (t386 * t525 + t279) * t289 + (-t335 * t74 - t336 * t73) * mrSges(6,2) - t590 * t288 + (t202 + (-t309 + t500) * t289) * t528 + t567 * t524) * t289 + (-m(6) * t197 - t366) * t127 + (-m(6) * t74 + t563 - t580) * t98 + (-m(6) * t73 + t429 - t566) * t97;
t486 = t9 * qJD(1);
t485 = t98 * t335;
t26 = t288 * t213 + m(6) * (t127 * t465 + t73 * t288) - t200 * t465;
t483 = qJD(1) * t26;
t480 = t120 * t335;
t479 = t121 * t336;
t470 = t271 * t198;
t463 = t299 * t127;
t460 = t325 * t199;
t459 = t325 * t301;
t457 = t335 * t127;
t456 = t335 * t200;
t454 = t335 * t211;
t453 = t335 * t212;
t447 = t336 * t297;
t439 = mrSges(6,2) * t488;
t438 = t77 * t541;
t437 = m(6) * t537;
t436 = m(5) * t519;
t435 = m(6) * t519;
t428 = mrSges(5,3) * t480;
t427 = mrSges(5,3) * t479;
t424 = t324 * t454;
t423 = t324 * t453;
t422 = mrSges(5,3) * t530;
t421 = mrSges(6,2) * t528;
t420 = t327 / 0.2e1;
t417 = t485 / 0.2e1;
t414 = t335 * t521;
t413 = t336 * t521;
t409 = -t467 / 0.2e1;
t398 = t300 + t518;
t395 = -t433 / 0.2e1;
t394 = t433 / 0.2e1;
t383 = t487 + t488;
t380 = t336 * t395;
t251 = t299 * t297;
t356 = t251 + t546;
t28 = t271 * t398 + t356 + t459;
t362 = t511 * t409 + t561 * t468 + t581 * t534 - t550;
t344 = (-pkin(4) * t77 + qJ(5) * t76) * t541 + t111 * t540 + t112 * t538 + t76 * mrSges(6,3) / 0.2e1 - t77 * mrSges(6,1) / 0.2e1 + t362;
t375 = t271 * t197 + t463;
t5 = ((t386 / 0.4e1 + t303 / 0.4e1 - t310 / 0.4e1 - t308 / 0.4e1 + t325 * t538 + mrSges(6,3) * t536 + (t422 + Ifges(6,3) / 0.4e1 + Ifges(5,2) / 0.4e1) * t335) * t335 + (t311 / 0.4e1 + t309 / 0.4e1 - t306 / 0.4e1 + t302 / 0.4e1 + mrSges(5,1) * t529 + mrSges(6,1) * t536 + (mrSges(6,2) * t530 + t422 + Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1) * t336 + (-Ifges(5,4) / 0.2e1 + Ifges(6,5) / 0.4e1) * t335) * t336 + t556 * t324 * (-t333 / 0.2e1 - t334 / 0.2e1)) * t289 + (-t304 / 0.4e1 - t305 / 0.4e1 + (mrSges(5,2) * t531 - Ifges(6,6) / 0.4e1 + Ifges(5,6) / 0.4e1) * t335 + (-Ifges(5,5) / 0.4e1 - Ifges(6,4) / 0.4e1 + (t539 + t540) * t324) * t336) * t288 + (-t202 / 0.4e1 + t201 / 0.4e1 + m(6) * t508 * t531 + (-t97 / 0.2e1 - t74 / 0.2e1) * mrSges(6,2)) * t336 + (t265 / 0.4e1 - t204 / 0.4e1 - t203 / 0.4e1 + (-t98 / 0.2e1 + t73 / 0.2e1) * mrSges(6,2) + (t509 * t542 + t213 / 0.2e1) * t324) * t335 + t375 * t542 + t344 + t547;
t379 = -t5 * qJD(1) + t28 * qJD(2);
t341 = -t521 * mrSges(4,2) + t556 * t592 + (-mrSges(4,1) + t553) * t519;
t35 = (t271 * t435 + t325 * t436 + t341) * pkin(2) + (m(6) + m(5)) * t562 * t324;
t376 = t479 - t480;
t337 = -m(5) * (t573 + t376 * t324 + (t413 * t98 - t414 * t97 + t576) * pkin(2)) / 0.2e1 + (t575 + t383 * t324 + (t127 * t519 + t413 * t73 + t414 * t74) * pkin(2)) * t542 + t470 / 0.2e1 + t460 / 0.2e1 + t428 / 0.2e1 - t427 / 0.2e1 + t424 / 0.2e1 - t423 / 0.2e1 - t439 / 0.2e1 + t81 * t421 - t549 / 0.2e1 - (-t200 - t367) * t430 / 0.2e1 + t563 * t335 * t394 + t566 * t380 - t587 - t588;
t4 = t337 + t555 * t514 / 0.2e1 + (t453 / 0.2e1 - t454 / 0.2e1) * pkin(8) + t545 + t587;
t378 = -t4 * qJD(1) + t35 * qJD(2);
t350 = (t539 - t447 / 0.2e1) * t289 - t456 / 0.2e1 - t426;
t353 = m(6) * (-t457 + (t271 * t289 + t288 * t324) * t336);
t18 = t438 - t353 / 0.2e1 + t350;
t209 = (m(6) * t271 + t297) * t335;
t374 = qJD(1) * t18 + qJD(2) * t209;
t32 = -t273 + 0.2e1 * (t98 / 0.4e1 - t469 / 0.2e1 - t455 / 0.4e1 - t448 / 0.4e1) * m(6);
t373 = qJD(1) * t32 - qJD(4) * t326;
t339 = (-pkin(4) * t414 + qJ(5) * t413) * pkin(2) * t541 + mrSges(5,2) * t380 + t394 * t507 + t557 * t335 * t395;
t345 = t251 - t516 / 0.2e1 + t459 / 0.2e1 + (t518 / 0.2e1 + t300 / 0.2e1) * t442;
t21 = t339 - t449 / 0.2e1 - t345 - t591;
t30 = t295 * t398 + t356 - t516;
t343 = -t547 + t306 * t465 / 0.4e1 - t303 * t466 / 0.4e1 - t335 * (Ifges(5,6) * t288 + t289 * t386) / 0.4e1 - t336 * t201 / 0.4e1 + mrSges(6,2) * t417 + t73 * t421 + t565 * t288 / 0.4e1 + t508 * t420 + (-t485 / 0.2e1 + t417) * mrSges(5,3) + (Ifges(6,6) * t288 - Ifges(6,3) * t466 - t265 + t567) * t335 / 0.4e1 + (t511 * t288 + t564 * t289 + t202) * t336 / 0.4e1 + (-t386 + t593) * t466 / 0.4e1 - (t302 - t564) * t465 / 0.4e1;
t338 = ((-t73 * t335 + t485 + t570) * t541 + (mrSges(6,2) / 0.2e1 + mrSges(5,3) / 0.2e1) * t586 + t548) * pkin(8) + t343 + (t295 * t197 + t463) * t541;
t346 = (-pkin(4) * t81 + qJ(5) * t80) * t542 - t120 * mrSges(5,1) / 0.2e1 + t121 * mrSges(5,2) / 0.2e1 - t80 * mrSges(6,3) / 0.2e1 + mrSges(6,1) * t537;
t348 = -pkin(3) * t298 / 0.2e1 + t297 * t533;
t8 = t338 + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t348) * t289 + ((Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t336 - t561 * t335) * t288 + t346 + t550;
t363 = t8 * qJD(1) - t21 * qJD(2) + t30 * qJD(3);
t129 = (t297 + (t394 + t536 + t533) * m(6)) * t335;
t222 = (m(6) * t295 + t297) * t335;
t354 = m(6) * (-t457 + (t289 * t295 + t513) * t336);
t23 = t437 - t354 / 0.2e1 + t350;
t361 = qJD(1) * t23 + qJD(2) * t129 + qJD(3) * t222;
t355 = t456 / 0.2e1 + t289 * t447 / 0.2e1 + t288 * t420 + t493 / 0.2e1 + mrSges(6,2) * t409;
t349 = (-mrSges(6,2) * t458 - pkin(4) * t327 + t565) * qJD(4);
t347 = qJD(4) * (m(6) * t385 + t553);
t296 = m(6) * t514 + t327;
t270 = m(6) * t461 + t327;
t130 = m(6) * t442 * t528 + (m(6) * t394 - t297) * t335;
t29 = 0.2e1 * t541 * t73 + t213;
t24 = t354 / 0.2e1 + t437 + t355;
t22 = t339 + t345 + t546;
t19 = t353 / 0.2e1 + t438 + t355;
t7 = t289 * t348 + t338 - t346 + t362;
t6 = t366 * t536 + t368 * t529 + t375 * t541 + t343 + t344 + ((t335 * t509 + t570) * t541 + t548) * t324 + t556 * t530 * t586;
t1 = t352 - t337 + (-mrSges(4,1) / 0.2e1 + t298 / 0.2e1) * t551 + ((t214 / 0.2e1 + t210 / 0.2e1) * t336 + (t584 - t211 / 0.2e1) * t335) * pkin(8) + t545;
t10 = [-qJD(2) * t2 - qJD(3) * t3 - qJD(4) * t9 + qJD(5) * t26, -t499 + (m(6) * (t324 * t384 + t575) + m(5) * (t324 * t377 + t573) + m(4) * (-t521 * t551 - t576) * pkin(2) + t585 + t549 - mrSges(3,1) * t434 - t460 + t430 * t289 * mrSges(4,3) + t433 * t494 + mrSges(3,2) * t431 - t470 - t424 + t423 - Ifges(3,6) * t520 + Ifges(3,5) * t522 + t384 * mrSges(6,2) + t377 * mrSges(5,3)) * qJD(2) + t1 * qJD(3) + t6 * qJD(4) + t19 * qJD(5), -t491 + t1 * qJD(2) + t7 * qJD(4) + t24 * qJD(5) + ((pkin(8) * t383 + t574) * t541 + (pkin(8) * t376 - t583) * t543) * t544 + (mrSges(6,2) * t487 + t427 - t428 + t439 - t464 + t517 + (t555 * t336 + (-t211 + t212) * t335) * pkin(8) + t585) * qJD(3), t6 * qJD(2) + t7 * qJD(3) + t29 * qJD(5) - t486 + ((-m(6) * pkin(4) - t557) * t98 + (-mrSges(5,2) + t326) * t97 + ((mrSges(6,2) * qJ(5) + t582) * t336 + (-pkin(4) * mrSges(6,2) + t511) * t335) * t289) * qJD(4), qJD(2) * t19 + qJD(3) * t24 + qJD(4) * t29 + t483; -qJD(3) * t4 - qJD(4) * t5 - qJD(5) * t18 + t499, qJD(3) * t35 + qJD(4) * t28 - qJD(5) * t209, t22 * qJD(4) + t130 * qJD(5) + (-pkin(3) * t436 + t295 * t435 + t341) * qJD(3) * pkin(2) + t378 + (t541 + t543) * t544 * t562 * pkin(8), t22 * qJD(3) + t270 * qJD(5) + t324 * t347 + t349 + t379, qJD(3) * t130 + qJD(4) * t270 - t374; qJD(2) * t4 + qJD(4) * t8 - qJD(5) * t23 + t491, -qJD(4) * t21 - qJD(5) * t129 - t378, qJD(4) * t30 - qJD(5) * t222, pkin(8) * t347 + t296 * qJD(5) + t349 + t363, qJD(4) * t296 - t361; qJD(2) * t5 - qJD(3) * t8 - qJD(5) * t32 + t486, qJD(3) * t21 - t379, -t363, t326 * qJD(5), -t373; qJD(2) * t18 + qJD(3) * t23 + qJD(4) * t32 - t483, qJD(3) * t129 + t374, t361, t373, 0;];
Cq = t10;
