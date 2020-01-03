% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:27
% EndTime: 2019-12-31 22:21:46
% DurationCPUTime: 9.95s
% Computational Cost: add. (23104->480), mult. (45719->644), div. (0->0), fcn. (51397->8), ass. (0->294)
t317 = cos(qJ(5));
t487 = t317 * mrSges(6,2);
t315 = sin(qJ(5));
t491 = t315 * mrSges(6,1);
t294 = t487 + t491;
t527 = sin(qJ(3));
t528 = sin(qJ(2));
t530 = cos(qJ(3));
t531 = cos(qJ(2));
t290 = -t527 * t531 - t530 * t528;
t316 = sin(qJ(4));
t346 = t527 * t528 - t530 * t531;
t529 = cos(qJ(4));
t333 = -t529 * t290 - t316 * t346;
t163 = t294 * t333;
t242 = t290 * t316 - t529 * t346;
t308 = -t531 * pkin(2) - pkin(1);
t264 = t346 * pkin(3) + t308;
t127 = -pkin(4) * t242 - pkin(9) * t333 + t264;
t438 = t528 * pkin(6);
t364 = -t528 * pkin(7) - t438;
t440 = t531 * pkin(6);
t365 = t531 * pkin(7) + t440;
t255 = t527 * t364 + t530 * t365;
t220 = -t346 * pkin(8) + t255;
t557 = t530 * t364 - t527 * t365;
t589 = t290 * pkin(8) + t557;
t614 = t529 * t220 + t316 * t589;
t48 = t127 * t317 - t315 * t614;
t592 = t294 * t242;
t593 = t242 * t317;
t601 = mrSges(6,1) * t333 - mrSges(6,3) * t593;
t615 = -t316 * t220 + t529 * t589;
t642 = t614 * t163 + t48 * t601 - t615 * t592;
t238 = Ifges(5,5) * t242;
t310 = Ifges(6,4) * t317;
t298 = Ifges(6,1) * t315 + t310;
t444 = t317 * t298;
t518 = Ifges(6,4) * t315;
t297 = Ifges(6,2) * t317 + t518;
t451 = t315 * t297;
t368 = t451 / 0.2e1 - t444 / 0.2e1;
t295 = Ifges(6,5) * t315 + t317 * Ifges(6,6);
t571 = t333 * t295;
t574 = Ifges(5,6) * t333;
t513 = Ifges(6,6) * t333;
t566 = -Ifges(6,2) * t315 + t310;
t91 = t242 * t566 + t513;
t595 = t317 * t91;
t386 = Ifges(6,1) * t317 - t518;
t517 = Ifges(6,5) * t333;
t93 = t386 * t242 + t517;
t596 = t315 * t93;
t600 = t571 / 0.2e1 - t574 + t595 / 0.2e1 + t596 / 0.2e1;
t352 = -Ifges(4,5) * t346 + Ifges(4,6) * t290 - t242 * t368 + t238 + t600;
t594 = t557 * mrSges(4,2);
t624 = t255 * mrSges(4,1);
t387 = mrSges(6,1) * t317 - t315 * mrSges(6,2);
t630 = t614 * t387;
t632 = t615 * mrSges(5,2);
t633 = t614 * mrSges(5,1);
t627 = -t630 - t633 - t632;
t640 = t627 + t352 - t594 - t624;
t639 = m(5) * t264 - mrSges(5,1) * t242 + mrSges(5,2) * t333;
t550 = pkin(4) / 0.2e1;
t586 = t571 / 0.4e1 - t574 / 0.2e1;
t351 = -t238 / 0.2e1 + t592 * t550 - t586;
t587 = -t595 / 0.4e1 - t596 / 0.4e1;
t467 = t242 * t315;
t575 = mrSges(6,2) * t333;
t167 = -mrSges(6,3) * t467 - t575;
t447 = t317 * t167;
t618 = t315 * t601;
t603 = t447 / 0.2e1 - t618 / 0.2e1;
t637 = t603 * pkin(9) - t351 - t587;
t549 = m(5) * pkin(3);
t621 = t615 * t316;
t636 = (-t529 * t614 + t621) * t549;
t49 = t315 * t127 + t317 * t614;
t532 = t317 / 0.2e1;
t533 = -t315 / 0.2e1;
t309 = Ifges(6,5) * t317;
t511 = Ifges(6,6) * t315;
t567 = t309 - t511;
t581 = -t242 / 0.2e1;
t580 = t333 / 0.2e1;
t583 = t264 * mrSges(5,1) + t567 * t580 - Ifges(5,4) * t333 + (Ifges(5,2) + Ifges(6,3)) * t581;
t366 = t386 * t333;
t516 = Ifges(6,5) * t242;
t94 = t366 - t516;
t485 = t317 * t94;
t512 = Ifges(6,6) * t242;
t92 = t333 * t566 - t512;
t489 = t315 * t92;
t585 = -t264 * mrSges(5,2) - Ifges(5,4) * t242 - Ifges(5,1) * t580 + t489 / 0.2e1 - t485 / 0.2e1;
t602 = Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t635 = t308 * (-t290 * mrSges(4,1) - t346 * mrSges(4,2)) + t49 * t167 + t346 ^ 2 * Ifges(4,4) + (t567 * t581 - t585) * t242 + ((Ifges(6,1) * t593 + t517) * t532 + (-Ifges(6,2) * t467 + t513) * t533 - t593 * t518 + t583 + t602 * t242) * t333 + t642;
t555 = t632 / 0.2e1 + t630 / 0.2e1 + t633 / 0.2e1;
t634 = pkin(4) * t614;
t501 = t242 * mrSges(5,3);
t307 = t530 * pkin(2) + pkin(3);
t415 = t527 * t316;
t275 = -pkin(2) * t415 + t529 * t307;
t270 = -pkin(4) - t275;
t612 = t270 * t614;
t439 = t529 * pkin(3);
t306 = -t439 - pkin(4);
t611 = t306 * t614;
t619 = t315 * t615;
t629 = t614 * t615;
t393 = t527 * t529;
t276 = pkin(2) * t393 + t316 * t307;
t622 = t615 * t276;
t620 = t615 * t317;
t331 = t368 + t386 * t533 - t317 * t566 / 0.2e1;
t626 = -t594 / 0.2e1 - t624 / 0.2e1;
t313 = t315 ^ 2;
t314 = t317 ^ 2;
t443 = t313 + t314;
t573 = t443 * mrSges(6,3);
t599 = -t615 / 0.2e1;
t598 = -t242 / 0.4e1;
t570 = (t309 / 0.2e1 - t511 / 0.2e1) * t242;
t281 = (t530 * t529 - t415) * pkin(2);
t402 = t443 * t281;
t565 = -t444 / 0.4e1 + t451 / 0.4e1;
t591 = t565 * t242;
t177 = pkin(4) * t333 - pkin(9) * t242;
t534 = -t306 / 0.2e1;
t510 = Ifges(6,3) * t333;
t579 = -t510 / 0.2e1;
t464 = t333 * t315;
t437 = mrSges(6,3) * t464;
t168 = mrSges(6,2) * t242 - t437;
t463 = t333 * t317;
t171 = -mrSges(6,1) * t242 - mrSges(6,3) * t463;
t543 = -t171 / 0.2e1;
t369 = t168 * t533 + t317 * t543;
t562 = t447 - t618;
t60 = t177 * t317 - t619;
t61 = t315 * t177 + t620;
t380 = -t60 * t315 + t61 * t317;
t526 = pkin(3) * t290;
t135 = t177 - t526;
t312 = t528 * pkin(2);
t128 = t135 + t312;
t50 = t128 * t317 - t619;
t51 = t315 * t128 + t620;
t382 = -t50 * t315 + t51 * t317;
t559 = (-t527 * mrSges(4,1) - t530 * mrSges(4,2)) * pkin(2);
t558 = mrSges(6,3) * t402;
t522 = mrSges(6,3) * t333;
t556 = (t314 / 0.2e1 + t313 / 0.2e1) * t522 - t369;
t554 = t451 * t598 + t242 * t444 / 0.4e1 + t637;
t164 = t297 * t333;
t165 = t298 * t333;
t553 = -t297 * t463 / 0.4e1 - t315 * t165 / 0.4e1 + t567 * t598 + t294 * t599 + t485 / 0.4e1 - t489 / 0.4e1 - (t566 + t298) * t464 / 0.4e1 + (t366 / 0.4e1 - t164 / 0.4e1) * t317;
t552 = -m(6) / 0.2e1;
t551 = m(6) / 0.2e1;
t548 = m(6) * pkin(3);
t547 = -mrSges(6,2) / 0.2e1;
t546 = -Ifges(5,5) / 0.2e1;
t545 = -t60 / 0.2e1;
t544 = t61 / 0.2e1;
t542 = -t270 / 0.2e1;
t541 = t270 / 0.2e1;
t540 = t275 / 0.2e1;
t280 = (t530 * t316 + t393) * pkin(2);
t539 = -t280 / 0.2e1;
t538 = t281 / 0.2e1;
t537 = t297 / 0.4e1;
t536 = -t298 / 0.4e1;
t525 = pkin(3) * t316;
t305 = pkin(9) + t525;
t535 = -t305 / 0.2e1;
t523 = pkin(4) * t294;
t327 = -Ifges(4,4) * t290 + (Ifges(4,1) - Ifges(4,2)) * t346;
t342 = pkin(2) * t346;
t1 = (-mrSges(4,2) * t312 + t327) * t290 + m(6) * (t48 * t50 + t49 * t51 - t629) + m(4) * t308 * t312 + t51 * t168 + t50 * t171 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t531) * t531 + (-pkin(1) * mrSges(3,1) + mrSges(4,1) * t342 + (Ifges(3,1) - Ifges(3,2)) * t531 - Ifges(3,4) * t528) * t528 + t635 + t639 * (t312 - t526);
t509 = t1 * qJD(1);
t502 = t333 * mrSges(5,3);
t496 = t275 * mrSges(5,2);
t495 = t276 * mrSges(5,1);
t494 = t280 * mrSges(5,1);
t493 = t281 * mrSges(5,2);
t54 = t135 * t317 - t619;
t55 = t315 * t135 + t620;
t3 = m(6) * (t48 * t54 + t49 * t55 - t629) + (-t639 * pkin(3) + t327) * t290 + t55 * t168 + t54 * t171 + t635;
t492 = t3 * qJD(1);
t490 = t315 * mrSges(6,3);
t482 = t54 * t315;
t481 = t55 * t317;
t166 = -t242 * t490 - t575;
t9 = m(6) * (t48 * t60 + t49 * t61 - t629) + t61 * t168 + t49 * t166 + t60 * t171 + (t93 * t532 + t91 * t533 + t583) * t333 - (-t602 * t333 + t570 + t585) * t242 + t642;
t478 = t9 * qJD(1);
t160 = t387 * t333;
t425 = -t512 / 0.2e1;
t10 = t615 * t160 + t49 * t171 + ((t94 / 0.2e1 - t164 / 0.2e1 - t516 / 0.2e1) * t315 + (t165 / 0.2e1 + t92 / 0.2e1 + t425 + t49 * mrSges(6,3)) * t317) * t333 + (-t168 - t437) * t48;
t477 = t10 * qJD(1);
t462 = t270 * t592;
t461 = t270 * t294;
t460 = t276 * t387;
t459 = t280 * t615;
t458 = t280 * t387;
t457 = t306 * t160;
t456 = t306 * t592;
t455 = t306 * t294;
t448 = t317 * t166;
t446 = t317 * t168;
t442 = mrSges(5,1) * t525;
t441 = mrSges(6,3) * t481;
t436 = t275 * t501;
t435 = t276 * t502;
t434 = t525 / 0.2e1;
t423 = -t490 / 0.2e1;
t422 = t490 / 0.2e1;
t421 = mrSges(6,3) * t532;
t419 = t510 / 0.2e1 + t570;
t418 = t315 * t529;
t417 = t317 * t529;
t403 = t443 * t275;
t401 = t443 * t305;
t400 = t502 * t525;
t398 = mrSges(5,2) * t439;
t397 = t439 / 0.2e1;
t392 = -t418 / 0.2e1;
t391 = t439 * t501;
t390 = t55 * t421 + t54 * t423 - t555;
t389 = t51 * t421 + t50 * t423 - t555;
t383 = -t315 * t48 + t317 * t49;
t381 = t481 - t482;
t378 = -t523 / 0.2e1 - t331;
t271 = pkin(9) + t276;
t334 = mrSges(6,3) * t403 - t460 - t495 - t496;
t31 = m(6) * (t270 * t276 + t271 * t403) + t334;
t339 = (-Ifges(5,6) / 0.2e1 + t295 / 0.4e1) * t333 + t592 * t541 + t276 * t163 / 0.2e1;
t344 = m(6) * (t382 * pkin(9) - t634);
t357 = t93 / 0.4e1 - t271 * t601 / 0.2e1 + t275 * t543;
t358 = t91 / 0.4e1 + t271 * t166 / 0.2e1 + t168 * t540;
t375 = t612 - t622;
t5 = (t599 + t615 / 0.2e1) * mrSges(5,2) + (-pkin(9) * t167 / 0.2e1 - t513 / 0.4e1 + (t536 - t310 / 0.4e1) * t242 + t358) * t317 + (pkin(9) * t601 / 0.2e1 - t517 / 0.4e1 + (t537 + t518 / 0.4e1 + (Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1) * t317) * t242 + t357) * t315 - t344 / 0.2e1 + (t380 * t271 + t383 * t275 + t375) * t551 - (t546 + t565) * t242 + ((t544 - t51 / 0.2e1) * t317 + (t545 + t50 / 0.2e1) * t315) * mrSges(6,3) + t339 + t351;
t377 = t5 * qJD(1) + t31 * qJD(2);
t32 = (-mrSges(5,1) - t387) * t280 + t559 + (-mrSges(5,2) + t573) * t281 + m(6) * (t270 * t280 + t271 * t402) + m(5) * (-t275 * t280 + t276 * t281);
t320 = -m(5) * (t622 - t459 + (-t275 + t281) * t614) / 0.2e1 + (t383 * t281 - t459 + t612) * t552 - t462 / 0.2e1 + t163 * t539 + t436 / 0.2e1 + t435 / 0.2e1 + t315 * t171 * t538 - t281 * t446 / 0.2e1 + (t381 * t552 - t603) * t271 + (-t242 * t538 + t333 * t539) * mrSges(5,3) - t626;
t322 = t551 * t611 + t456 / 0.2e1 + t636 / 0.2e1 - t400 / 0.2e1 - t391 / 0.2e1 + t389 + (t382 * t551 + t603) * t305 + t626;
t332 = t54 * t422 - t441 / 0.2e1 + t555;
t6 = t322 + t320 + t332;
t376 = -t6 * qJD(1) + t32 * qJD(2);
t338 = t553 + (t422 + t423) * t49;
t324 = t160 * t541 - t271 * t556 + t338;
t373 = -t50 * mrSges(6,1) / 0.2e1 + t51 * mrSges(6,2) / 0.2e1;
t12 = t324 + t373 + t579 - t570;
t145 = -t331 + t461;
t374 = t12 * qJD(1) + t145 * qJD(2);
t370 = -t487 / 0.2e1 - t491 / 0.2e1;
t362 = t443 * t529;
t361 = t370 * t275;
t360 = t370 * t281;
t330 = -t387 * t525 + t439 * t573 - t398 - t442;
t146 = (t362 * t305 + t306 * t316) * t548 + t330;
t321 = (t306 * t276 + t275 * t401 + (t270 * t316 + t362 * t271) * pkin(3)) * t551 - t496 / 0.2e1 - t495 / 0.2e1 - t460 / 0.2e1 - t442 / 0.2e1 - t387 * t434 - t398 / 0.2e1 + (t397 + t540) * t573;
t326 = (-pkin(4) * t280 + pkin(9) * t402) * t552 + t494 / 0.2e1 + t458 / 0.2e1 + t493 / 0.2e1 - t558 / 0.2e1;
t22 = t321 + t326;
t336 = -t242 * t546 - t555;
t318 = (t611 + t380 * t305 + (t49 * t417 - t48 * t418 - t621) * pkin(3)) * t551 - t592 * t534 + t618 * t535 + t163 * t434 + t305 * t448 / 0.2e1 + t60 * t423 + t61 * t421 + pkin(3) * t171 * t392 + t397 * t446 + t336 - t591 + t586 - t587;
t345 = m(6) * (t381 * pkin(9) - t634);
t8 = t318 - t345 / 0.2e1 + t332 + t591 - t637;
t356 = t8 * qJD(1) + t22 * qJD(2) + t146 * qJD(3);
t353 = t54 * mrSges(6,1) / 0.2e1 + t55 * t547 + t419;
t13 = -t457 / 0.2e1 + t556 * t305 + t353 - t553;
t172 = -t331 + t455;
t72 = (t534 + t542) * t294 + t360 + t331;
t355 = -t13 * qJD(1) - t72 * qJD(2) + t172 * qJD(3);
t354 = mrSges(6,1) * t545 + mrSges(6,2) * t544 + t579;
t325 = -t556 * pkin(9) + t338 - pkin(4) * t160 / 0.2e1;
t16 = t325 + t354 - t570;
t178 = t331 + t523;
t74 = (t550 + t542) * t294 + t361 + t331;
t337 = (mrSges(6,1) * t392 + t417 * t547) * pkin(3);
t99 = (t550 + t534) * t294 + t337 + t331;
t343 = t16 * qJD(1) - t74 * qJD(2) - t99 * qJD(3) - t178 * qJD(4);
t274 = t455 / 0.2e1;
t245 = t461 / 0.2e1;
t100 = t274 + t337 + t378;
t75 = t245 + t361 + t378;
t73 = t274 + t245 + t360 - t331;
t21 = t321 - t326;
t15 = t325 + Ifges(6,5) * t593 / 0.2e1 + t315 * t425 - t354;
t14 = t338 + t457 / 0.2e1 + t353 + t443 * t522 * t535 + t369 * t305;
t11 = t324 - t373 + t419;
t7 = t318 + t345 / 0.2e1 + t390 + t554;
t4 = t554 + t339 + t344 / 0.2e1 + t336 + ((-t271 * t60 - t275 * t48) * t551 + mrSges(6,3) * t545 - t242 * t537 + t357) * t315 + ((t271 * t61 + t275 * t49) * t551 - t242 * t536 + mrSges(6,3) * t544 + t358) * t317 + t389 + t375 * t551;
t2 = t322 - t320 + t352 + t390;
t17 = [qJD(2) * t1 + qJD(3) * t3 + qJD(4) * t9 - qJD(5) * t10, t509 + (mrSges(3,2) * t438 - mrSges(3,1) * t440 + Ifges(3,5) * t531 - Ifges(3,6) * t528 + m(5) * (-t275 * t614 + t622) - t435 - t436 + m(6) * t612 + m(4) * (-t255 * t530 + t527 * t557) * pkin(2) + t462 + (m(6) * t382 + t562) * t271 + t382 * mrSges(6,3) + (t527 * pkin(2) * t290 + t530 * t342) * mrSges(4,3) + t640) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t11 * qJD(5), t492 + t2 * qJD(2) + (t441 - mrSges(6,3) * t482 + m(6) * t611 + t636 - t391 - t400 + t456 + (m(6) * t381 + t562) * t305 + t640) * qJD(3) + t7 * qJD(4) + t14 * qJD(5), t4 * qJD(2) + t7 * qJD(3) + t15 * qJD(5) + t478 + (-(-Ifges(5,5) + t368) * t242 + (-m(6) * t614 - t592) * pkin(4) + (m(6) * t380 + t448 - t618) * pkin(9) + t380 * mrSges(6,3) + t600 + t627) * qJD(4), -t477 + t11 * qJD(2) + t14 * qJD(3) + t15 * qJD(4) + (-t49 * mrSges(6,1) - t48 * mrSges(6,2) - t571) * qJD(5); -qJD(3) * t6 + qJD(4) * t5 + qJD(5) * t12 - t509, qJD(3) * t32 + qJD(4) * t31 + qJD(5) * t145, (m(6) * (t280 * t306 + t281 * t401) - t458 - t494 + (-t529 * t280 + t281 * t316) * t549 - t493 + t558 + t559) * qJD(3) + t21 * qJD(4) + t73 * qJD(5) + t376, t21 * qJD(3) + (m(6) * (-pkin(4) * t276 + pkin(9) * t403) + t334) * qJD(4) + t75 * qJD(5) + t377, t73 * qJD(3) + t75 * qJD(4) + (-t271 * t387 + t567) * qJD(5) + t374; qJD(2) * t6 + qJD(4) * t8 - qJD(5) * t13 - t492, qJD(4) * t22 - qJD(5) * t72 - t376, qJD(4) * t146 + qJD(5) * t172, ((-pkin(4) * t316 + pkin(9) * t362) * t548 + t330) * qJD(4) + t100 * qJD(5) + t356, t100 * qJD(4) + (-t305 * t387 + t567) * qJD(5) + t355; -qJD(2) * t5 - qJD(3) * t8 + qJD(5) * t16 - t478, -qJD(3) * t22 - qJD(5) * t74 - t377, -qJD(5) * t99 - t356, -t178 * qJD(5), (-pkin(9) * t387 + t567) * qJD(5) + t343; -qJD(2) * t12 + qJD(3) * t13 - qJD(4) * t16 + t477, qJD(3) * t72 + qJD(4) * t74 - t374, qJD(4) * t99 - t355, -t343, 0;];
Cq = t17;
