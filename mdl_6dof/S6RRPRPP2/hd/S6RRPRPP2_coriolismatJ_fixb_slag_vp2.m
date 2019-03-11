% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPP2
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:51
% EndTime: 2019-03-09 09:50:03
% DurationCPUTime: 8.16s
% Computational Cost: add. (12225->610), mult. (24341->773), div. (0->0), fcn. (25201->6), ass. (0->291)
t539 = -m(7) / 0.2e1;
t512 = sin(qJ(2));
t439 = t512 * pkin(2);
t567 = m(4) * t439;
t556 = Ifges(6,4) + Ifges(5,5);
t473 = sin(pkin(9));
t474 = cos(pkin(9));
t513 = cos(qJ(2));
t322 = t473 * t512 - t474 * t513;
t500 = mrSges(4,3) * t322;
t366 = sin(qJ(4));
t528 = pkin(4) + pkin(5);
t560 = t528 * t366;
t566 = t322 * t560;
t364 = t366 ^ 2;
t367 = cos(qJ(4));
t365 = t367 ^ 2;
t445 = t364 + t365;
t324 = -t473 * t513 - t474 * t512;
t462 = t324 * t366;
t298 = mrSges(7,3) * t462;
t488 = t322 * mrSges(7,2);
t224 = -t298 + t488;
t461 = t324 * t367;
t229 = -mrSges(7,1) * t322 + mrSges(7,3) * t461;
t514 = t367 / 0.2e1;
t516 = t366 / 0.2e1;
t451 = t224 * t514 + t229 * t516;
t413 = -pkin(2) * t513 - pkin(1);
t216 = t322 * pkin(3) + t324 * pkin(8) + t413;
t438 = t512 * pkin(7);
t329 = -qJ(3) * t512 - t438;
t440 = t513 * pkin(7);
t334 = qJ(3) * t513 + t440;
t550 = t473 * t329 + t474 * t334;
t100 = t216 * t366 + t367 * t550;
t310 = t322 * qJ(5);
t80 = t100 + t310;
t475 = t100 - t80;
t447 = -t367 * t216 + t366 * t550;
t81 = -pkin(4) * t322 + t447;
t501 = t81 - t447;
t457 = t366 * qJ(6);
t295 = t324 * t457;
t62 = -t295 + t80;
t77 = t100 - t295;
t502 = t62 - t77;
t393 = qJ(6) * t461 + t447;
t51 = -t322 * t528 + t393;
t503 = t51 - t393;
t541 = -m(6) / 0.2e1;
t565 = t451 - (t366 * t501 - t367 * t475) * t541 - (t366 * t503 + t367 * t502) * t539;
t455 = t367 * qJ(5);
t564 = t560 - t455;
t427 = t473 * pkin(2);
t353 = t427 + pkin(8);
t460 = t353 * t366;
t315 = -t457 + t460;
t316 = (-qJ(6) + t353) * t367;
t404 = -t366 * t51 - t367 * t62;
t563 = t451 + ((t315 * t367 - t316 * t366) * t324 + t404) * t539;
t428 = t474 * pkin(2);
t354 = -t428 - pkin(3);
t458 = t366 * qJ(5);
t385 = t354 - t458;
t507 = t367 * pkin(4);
t314 = t385 - t507;
t327 = -t367 * mrSges(6,1) - t366 * mrSges(6,3);
t562 = m(6) * t314 + t327;
t561 = mrSges(6,2) + mrSges(5,3);
t218 = -t324 * pkin(3) + t322 * pkin(8) + t439;
t252 = -t474 * t329 + t334 * t473;
t234 = t366 * t252;
t113 = t218 * t367 + t234;
t114 = t366 * t218 - t252 * t367;
t559 = -t113 * t366 + t114 * t367;
t83 = -t324 * qJ(5) + t114;
t84 = pkin(4) * t324 - t113;
t558 = t366 * t84 + t367 * t83;
t529 = m(6) + m(7);
t328 = mrSges(7,1) * t367 + mrSges(7,2) * t366;
t206 = t324 * t328;
t557 = t206 / 0.2e1;
t504 = Ifges(5,6) + Ifges(7,6);
t555 = mrSges(7,3) * t445;
t484 = t366 * mrSges(5,1);
t333 = t367 * mrSges(5,2) + t484;
t554 = mrSges(4,3) + t333;
t483 = t366 * mrSges(6,1);
t331 = -t367 * mrSges(6,3) + t483;
t208 = t331 * t324;
t356 = t366 * mrSges(7,1);
t481 = t367 * mrSges(7,2);
t408 = -t481 + t356;
t209 = t408 * t324;
t553 = t208 + t209;
t464 = t322 * t366;
t222 = -mrSges(7,2) * t324 - mrSges(7,3) * t464;
t233 = mrSges(6,2) * t464 - mrSges(6,3) * t324;
t552 = t222 + t233;
t225 = -mrSges(5,2) * t322 + mrSges(5,3) * t462;
t312 = t322 * mrSges(6,3);
t437 = mrSges(6,2) * t462;
t232 = t312 + t437;
t449 = t225 + t232;
t358 = Ifges(6,5) * t366;
t335 = -t367 * Ifges(6,3) + t358;
t344 = Ifges(6,1) * t367 + t358;
t551 = Ifges(6,6) * t366 + t556 * t367;
t360 = Ifges(7,4) * t366;
t337 = -t367 * Ifges(7,2) + t360;
t342 = Ifges(7,1) * t367 + t360;
t362 = Ifges(5,4) * t367;
t345 = Ifges(5,1) * t366 + t362;
t330 = -t366 * pkin(4) + t455;
t405 = Ifges(5,2) * t366 - t362;
t548 = Ifges(7,5) - t556;
t547 = Ifges(7,3) + Ifges(6,2) + Ifges(5,3);
t506 = mrSges(6,2) - mrSges(7,3);
t546 = qJ(5) * t506 + t504;
t297 = t322 * t455;
t535 = mrSges(7,1) / 0.2e1;
t536 = mrSges(6,1) / 0.2e1;
t435 = t536 + t535;
t441 = pkin(4) * t464;
t532 = -mrSges(6,3) / 0.2e1;
t533 = -mrSges(7,2) / 0.2e1;
t534 = mrSges(5,2) / 0.2e1;
t538 = m(7) / 0.2e1;
t540 = m(6) / 0.2e1;
t545 = (-t297 + t441) * t540 + (-t297 + t566) * t538 + ((t534 + t532 + t533) * t367 + (mrSges(5,1) / 0.2e1 + t435) * t366) * t322;
t544 = 0.2e1 * m(7);
t543 = 2 * qJD(4);
t542 = m(5) / 0.2e1;
t537 = m(4) * pkin(2);
t531 = -mrSges(7,3) / 0.2e1;
t530 = Ifges(6,6) / 0.2e1;
t526 = t100 / 0.2e1;
t296 = t324 * t458;
t431 = t528 * t367;
t137 = t324 * t431 + t296;
t525 = t137 / 0.2e1;
t524 = t208 / 0.2e1;
t523 = t209 / 0.2e1;
t463 = t322 * t367;
t485 = t324 * mrSges(7,1);
t226 = mrSges(7,3) * t463 + t485;
t522 = -t226 / 0.2e1;
t521 = -t315 / 0.2e1;
t520 = -t322 / 0.2e1;
t518 = -t324 / 0.2e1;
t517 = t327 / 0.2e1;
t515 = -t367 / 0.2e1;
t511 = m(6) * t330;
t321 = -pkin(5) * t366 + t330;
t510 = m(7) * t321;
t509 = m(7) * t324;
t508 = m(7) * t366;
t505 = mrSges(7,2) + mrSges(6,3);
t499 = mrSges(4,3) * t324;
t495 = Ifges(5,4) * t366;
t494 = Ifges(7,4) * t367;
t493 = Ifges(6,5) * t367;
t126 = t324 * t330 + t252;
t230 = mrSges(5,1) * t322 + mrSges(5,3) * t461;
t231 = -mrSges(6,1) * t322 - mrSges(6,2) * t461;
t448 = -t230 + t231;
t467 = t252 * t324;
t87 = t564 * t324 - t252;
t5 = (t324 * t554 + t553) * t324 + (t500 + (-t224 - t449) * t367 + (-t229 - t448) * t366) * t322 + m(6) * (-t126 * t324 + (-t366 * t81 - t367 * t80) * t322) + m(7) * (t322 * t404 + t324 * t87) + m(5) * (-t467 + (-t100 * t367 - t366 * t447) * t322) + m(4) * (-t322 * t550 - t467);
t491 = qJD(1) * t5;
t223 = t324 * mrSges(5,2) + mrSges(5,3) * t464;
t227 = -mrSges(5,1) * t324 + mrSges(5,3) * t463;
t486 = t324 * mrSges(6,1);
t228 = -mrSges(6,2) * t463 + t486;
t205 = t324 * t327;
t409 = -t367 * mrSges(5,1) + t366 * mrSges(5,2);
t207 = t324 * t409;
t266 = t431 - t385;
t401 = -t315 * t366 - t316 * t367;
t465 = t322 * t353;
t411 = t445 * t465;
t466 = t314 * t324;
t369 = (-t324 * t354 - t411) * t542 + (-t411 - t466) * t540 + (t266 * t324 + t322 * t401) * t538 - t205 / 0.2e1 + t557 - t207 / 0.2e1 + (-t322 * t473 + t324 * t474) * t537 / 0.2e1 + t322 * t555 / 0.2e1 + t561 * t445 * t520;
t52 = -t234 + (qJ(6) * t322 - t218) * t367 + t528 * t324;
t63 = -t322 * t457 + t83;
t375 = (t113 * t367 + t114 * t366) * t542 + (t366 * t83 - t367 * t84) * t540 + (t366 * t63 - t367 * t52) * t538 + t567 / 0.2e1;
t414 = -t324 * mrSges(4,1) - t322 * mrSges(4,2);
t8 = t227 * t514 - t369 + t375 + t414 + (t226 + t228) * t515 + (t223 + t552) * t516;
t490 = qJD(1) * t8;
t390 = -t297 - t550;
t125 = -t390 - t441;
t487 = t322 * Ifges(7,5);
t144 = -t324 * t342 - t487;
t145 = t322 * Ifges(6,4) - t324 * t344;
t346 = Ifges(5,1) * t367 - t495;
t146 = t322 * Ifges(5,5) - t324 * t346;
t434 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t383 = t146 / 0.2e1 + t145 / 0.2e1 + t144 / 0.2e1 + (-Ifges(7,5) / 0.2e1 + t434) * t322;
t386 = t408 * t322;
t387 = t331 * t322;
t388 = t333 * t322;
t300 = Ifges(6,5) * t461;
t141 = Ifges(6,6) * t322 - Ifges(6,3) * t462 - t300;
t301 = Ifges(7,4) * t461;
t142 = -Ifges(7,2) * t462 - t322 * Ifges(7,6) - t301;
t143 = t322 * Ifges(5,6) + t324 * t405;
t410 = t141 / 0.2e1 + t142 / 0.2e1 - t143 / 0.2e1;
t433 = Ifges(5,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t86 = t390 + t566;
t3 = (-mrSges(4,1) * t439 + t410 * t366 + t383 * t367 + (-Ifges(4,4) + (t530 - t433) * t366) * t322) * t322 - m(5) * (t100 * t114 - t113 * t447 + t252 * t550) + (-t414 - t567) * t413 + t252 * t388 + t126 * t387 - t87 * t386 + (-Ifges(3,1) + Ifges(3,2)) * t513 * t512 + pkin(1) * (mrSges(3,1) * t512 + mrSges(3,2) * t513) + (t512 ^ 2 - t513 ^ 2) * Ifges(3,4) - m(6) * (t125 * t126 + t80 * t83 + t81 * t84) - m(7) * (t51 * t52 + t62 * t63 + t86 * t87) + t447 * t227 + (t548 * t461 + (-Ifges(6,6) + t504) * t462 + (Ifges(4,2) - Ifges(4,1) + (-Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(7,1) / 0.2e1) * t365 + ((-Ifges(5,2) / 0.2e1 - Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t366 + (Ifges(5,4) - Ifges(7,4) - Ifges(6,5)) * t367) * t366 + t547) * t322 + mrSges(4,2) * t439 + Ifges(4,4) * t324 + t554 * t550) * t324 - t550 * t499 + t125 * t208 - t86 * t209 - t62 * t222 - t100 * t223 - t63 * t224 - t114 * t225 - t51 * t226 - t81 * t228 - t52 * t229 - t113 * t230 - t84 * t231 - t83 * t232 - t80 * t233;
t489 = t3 * qJD(1);
t204 = -pkin(4) * t461 - t296;
t210 = t324 * t335;
t211 = t324 * t337;
t339 = t367 * Ifges(5,2) + t495;
t212 = t324 * t339;
t213 = Ifges(7,1) * t462 - t301;
t214 = Ifges(6,1) * t462 - t300;
t215 = t324 * t345;
t299 = Ifges(6,6) * t461;
t4 = t299 * t520 + t126 * t205 + t87 * t206 - t204 * t208 + t137 * t209 - t393 * t224 + t77 * t229 + t252 * t207 - t449 * t447 + t448 * t100 + m(6) * (t100 * t81 + t126 * t204 - t447 * t80) + m(7) * (t137 * t87 - t393 * t62 + t51 * t77) + ((-t213 / 0.2e1 - t214 / 0.2e1 - t215 / 0.2e1 - t62 * mrSges(7,3) + t80 * mrSges(6,2) + t100 * mrSges(5,3) + t433 * t322 - t410) * t367 + (-t210 / 0.2e1 - t211 / 0.2e1 + t212 / 0.2e1 - t51 * mrSges(7,3) + t81 * mrSges(6,2) + t447 * mrSges(5,3) + t383) * t366) * t324;
t477 = t4 * qJD(1);
t418 = -t232 / 0.2e1 - t225 / 0.2e1;
t419 = -t231 / 0.2e1 + t230 / 0.2e1;
t6 = t419 * t366 + t418 * t367 + t545 - t565;
t476 = t6 * qJD(1);
t14 = (-m(6) * t126 + m(7) * t87 + t553) * t461 + (-m(6) * t80 - m(7) * t62 - t224 - t232) * t322;
t472 = qJD(1) * t14;
t23 = (-t366 * t224 + m(7) * (-t62 * t366 + t51 * t367) + t367 * t229) * t324;
t471 = qJD(1) * t23;
t347 = t529 * t366;
t95 = t322 * t347;
t452 = t95 * qJD(1);
t444 = qJD(4) * t367;
t129 = 0.2e1 * (0.1e1 / 0.4e1 + t364 / 0.4e1 + t365 / 0.4e1) * t509;
t443 = t129 * qJD(1);
t442 = m(7) * t461;
t430 = t528 * t461;
t426 = -t464 / 0.2e1;
t425 = t464 / 0.2e1;
t424 = -t463 / 0.2e1;
t423 = t463 / 0.2e1;
t420 = -t222 / 0.2e1 - t233 / 0.2e1;
t417 = t517 - t328 / 0.2e1;
t338 = Ifges(7,2) * t366 + t494;
t336 = Ifges(6,3) * t366 + t493;
t370 = (-t126 * t330 + t204 * t314) * t540 + (t137 * t266 - t315 * t502 + t316 * t503 + t321 * t87) * t538 + t126 * t331 / 0.2e1 + t328 * t525 + t204 * t517 + t252 * t333 / 0.2e1 + t266 * t557 + t314 * t205 / 0.2e1 + t224 * t521 + t316 * t229 / 0.2e1 + t321 * t523 + t330 * t524 + t354 * t207 / 0.2e1 - t87 * t408 / 0.2e1 + t551 * t322 / 0.4e1;
t371 = (-pkin(4) * t84 + qJ(5) * t83) * t541 + (qJ(5) * t63 - t52 * t528) * t539 + pkin(4) * t228 / 0.2e1 - t113 * mrSges(5,1) / 0.2e1 + t114 * t534 - t528 * t522 + t52 * t535 + t63 * t533 + t83 * t532 + t84 * t536;
t372 = (-t77 / 0.2e1 + t62 / 0.2e1) * mrSges(7,3) + (-t80 / 0.2e1 + t526) * mrSges(6,2) + (t475 * t540 + t418) * t353 + t141 / 0.4e1 + t142 / 0.4e1 - t143 / 0.4e1 + t213 / 0.4e1 + t214 / 0.4e1 + t215 / 0.4e1;
t373 = (t393 / 0.2e1 - t51 / 0.2e1) * mrSges(7,3) + (-t447 / 0.2e1 + t81 / 0.2e1) * mrSges(6,2) + (t501 * t540 - t419) * t353 + t144 / 0.4e1 + t145 / 0.4e1 + t146 / 0.4e1 - t210 / 0.4e1 - t211 / 0.4e1 + t212 / 0.4e1;
t380 = -t346 / 0.4e1 - t344 / 0.4e1 - t342 / 0.4e1 + t339 / 0.4e1 - t337 / 0.4e1 - t335 / 0.4e1 + t316 * t531;
t341 = Ifges(7,1) * t366 - t494;
t343 = Ifges(6,1) * t366 - t493;
t381 = t345 / 0.4e1 + t343 / 0.4e1 + t341 / 0.4e1 - t405 / 0.4e1 - t338 / 0.4e1 - t336 / 0.4e1 + mrSges(7,3) * t521;
t382 = t561 * (t364 / 0.2e1 + t365 / 0.2e1) * t353;
t2 = t370 + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t382) * t324 + t420 * qJ(5) + ((-0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(7,6) + t530) * t322 + t381 * t324 + t372) * t366 + ((-0.3e1 / 0.4e1 * Ifges(7,5) + t434) * t322 + t380 * t324 + t373) * t367 + t371;
t20 = t321 * t328 - t330 * t327 + t354 * t333 + (t331 - t511) * t314 + (-t408 + t510) * t266 + (t345 / 0.2e1 - t405 / 0.2e1 + t343 / 0.2e1 - t336 / 0.2e1 + t341 / 0.2e1 - t338 / 0.2e1) * t367 + (t346 / 0.2e1 - t339 / 0.2e1 + t344 / 0.2e1 + t335 / 0.2e1 + t342 / 0.2e1 + t337 / 0.2e1) * t366;
t403 = t2 * qJD(1) + t20 * qJD(2);
t374 = (t524 + t523) * t366 + (-t126 * t366 + (t465 + t466) * t367) * t540 + (-t266 * t461 + t316 * t322 + t366 * t87) * t538;
t396 = t52 * t539 + t84 * t541;
t11 = t506 * t463 + (t367 * t417 - t435) * t324 + t374 + t396;
t115 = (m(7) * t266 + t328 - t562) * t366;
t400 = qJD(1) * t11 + qJD(2) * t115;
t147 = m(7) * t401 + t555;
t384 = (-t481 / 0.2e1 + t356 / 0.2e1) * t322 + t86 * t538;
t15 = t384 + t563;
t399 = -qJD(1) * t15 + qJD(2) * t147;
t164 = (t321 / 0.4e1 + t455 / 0.4e1 - t560 / 0.4e1) * t544 - t408;
t41 = -t206 + (-t430 / 0.4e1 - t296 / 0.4e1 - t137 / 0.4e1) * t544;
t398 = qJD(1) * t41 - qJD(2) * t164;
t389 = t100 + 0.2e1 * t310;
t376 = -t312 + t389 * t541 + (-t295 + t389) * t539 - t488;
t394 = m(6) * t526 + t538 * t77;
t22 = t376 + t394;
t348 = qJ(5) * t529 + t505;
t397 = qJD(1) * t22 - qJD(4) * t348;
t395 = -mrSges(6,2) * pkin(4) + mrSges(7,3) * t528 - Ifges(7,5);
t255 = (-qJD(1) * t461 + qJD(2) * t366) * m(7);
t219 = -t564 * t539 + t510 / 0.2e1;
t181 = m(7) * t316 + (m(6) * t353 + t506) * t367;
t130 = (-t445 + 0.1e1) * t509 / 0.2e1;
t94 = (t540 + t538) * t464 + t529 * t426;
t65 = (-t296 - t430) * t538 + m(7) * t525;
t19 = -t298 - t376 + t394 + t437;
t16 = t384 - t563;
t10 = mrSges(7,3) * t423 + t485 / 0.2e1 + mrSges(6,2) * t424 + t486 / 0.2e1 + (t417 * t324 + (mrSges(6,2) / 0.2e1 + t531) * t322) * t367 + t374 - t396;
t9 = (t223 / 0.2e1 - t420) * t366 + (t522 - t228 / 0.2e1 + t227 / 0.2e1) * t367 + t369 + t375;
t7 = t231 * t516 - t366 * t230 / 0.2e1 + t449 * t514 + t545 + t565;
t1 = (t366 * t381 + t367 * t380 + t382) * t324 + t370 + Ifges(7,5) * t423 + Ifges(6,6) * t426 + (-t487 / 0.4e1 + t373) * t367 + ((-Ifges(7,6) / 0.4e1 - Ifges(5,6) / 0.4e1) * t322 + t372) * t366 - t371 + t552 * qJ(5) / 0.2e1 + t504 * t425 + t556 * t424 + t547 * t518;
t12 = [-qJD(2) * t3 + qJD(3) * t5 + qJD(4) * t4 - qJD(5) * t14 + qJD(6) * t23, t9 * qJD(3) + t1 * qJD(4) + t10 * qJD(5) + t16 * qJD(6) - t489 + ((t548 * t324 + (-t346 - t342 - t344) * t322) * t516 + ((-Ifges(6,6) + Ifges(7,6)) * t324 + (-t336 - t338) * t322) * t515 + t339 * t425 - mrSges(3,1) * t440 + ((Ifges(5,6) - Ifges(6,6)) * t367 + t556 * t366) * t518 + t558 * mrSges(6,2) + t559 * mrSges(5,3) + (-t227 + t228) * t460 + m(7) * (t266 * t86 + t315 * t52 + t316 * t63) + (t337 + t335) * t426 + (t345 + t343 + t341) * t424 + (-Ifges(5,6) * t324 + t322 * t405) * t514 + t427 * t499 + t428 * t500 + t562 * t125 + (-t366 * t52 - t367 * t63) * mrSges(7,3) + (m(5) * t559 + m(6) * t558 + (t233 + t223) * t367) * t353 + mrSges(3,2) * t438 + t266 * t386 - t314 * t387 - t354 * t388 - Ifges(3,6) * t512 + Ifges(3,5) * t513 + t324 * (Ifges(7,5) * t366 - Ifges(7,6) * t367) / 0.2e1 + (m(5) * t354 - t474 * t537 - mrSges(4,1) + t409) * t550 + (-t473 * t537 + mrSges(4,2)) * t252 + t315 * t226 + t316 * t222 - Ifges(4,5) * t322 + Ifges(4,6) * t324 + t86 * t328) * qJD(2), qJD(2) * t9 + qJD(4) * t7 + qJD(5) * t94 + qJD(6) * t130 + t491, t477 + t1 * qJD(2) + t7 * qJD(3) + t19 * qJD(5) + t65 * qJD(6) + ((-pkin(4) * t100 - qJ(5) * t447) * t540 + (-qJ(5) * t393 - t528 * t77) * t538) * t543 + (-t77 * mrSges(7,1) - t393 * mrSges(7,2) - t299 + (t546 * t367 + (t395 + t556) * t366) * t324 - (-mrSges(5,2) + mrSges(6,3)) * t447 + (-mrSges(5,1) - mrSges(6,1)) * t100) * qJD(4), qJD(2) * t10 + qJD(3) * t94 + qJD(4) * t19 - t472, qJD(2) * t16 + qJD(3) * t130 + qJD(4) * t65 + t471; -qJD(3) * t8 + qJD(4) * t2 + qJD(5) * t11 - qJD(6) * t15 + t489, qJD(4) * t20 + qJD(5) * t115 + qJD(6) * t147, -t490, t181 * qJD(5) + t219 * qJD(6) + t395 * t444 + t403 + (-t316 * mrSges(7,1) - t315 * mrSges(7,2) + m(7) * (-qJ(5) * t315 - t316 * t528) - t546 * t366 + (m(6) * (-t458 - t507) + t327 + t409) * t353 + t551) * qJD(4), qJD(4) * t181 + t400, qJD(4) * t219 + t399; qJD(2) * t8 - qJD(4) * t6 + qJD(5) * t95 - qJD(6) * t129 - t491, t490, 0, -t476 + (-t356 - t483 - t484) * qJD(4) + t347 * qJD(5) + (-mrSges(5,2) + t505) * t444 + (t511 / 0.2e1 - t564 * t538) * t543, qJD(4) * t347 + t452, -t443; -qJD(2) * t2 + qJD(3) * t6 - qJD(5) * t22 + qJD(6) * t41 - t477, -t164 * qJD(6) - t403, t476, t348 * qJD(5), -t397, t398; -qJD(2) * t11 - qJD(3) * t95 + qJD(4) * t22 + qJD(6) * t442 + t472, -qJD(6) * t508 - t400, -t452, t397, 0, -t255; qJD(2) * t15 + qJD(3) * t129 - qJD(4) * t41 - qJD(5) * t442 - t471, qJD(4) * t164 + qJD(5) * t508 - t399, t443, -t398, t255, 0;];
Cq  = t12;
