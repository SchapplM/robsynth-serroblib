% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP6
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:40
% EndTime: 2019-03-09 00:29:49
% DurationCPUTime: 41.00s
% Computational Cost: add. (14449->946), mult. (36430->1280), div. (0->0), fcn. (30826->14), ass. (0->431)
t350 = cos(pkin(7));
t340 = qJD(2) * t350 + qJD(3);
t353 = sin(qJ(4));
t357 = cos(qJ(4));
t354 = sin(qJ(3));
t348 = sin(pkin(7));
t497 = qJD(2) * t348;
t466 = t354 * t497;
t264 = t340 * t357 - t353 * t466;
t358 = cos(qJ(3));
t485 = qJD(2) * qJD(3);
t297 = (qJDD(2) * t354 + t358 * t485) * t348;
t339 = qJDD(2) * t350 + qJDD(3);
t170 = qJD(4) * t264 + t297 * t357 + t339 * t353;
t588 = t170 / 0.2e1;
t672 = Ifges(5,4) * t588;
t296 = (-qJDD(2) * t358 + t354 * t485) * t348;
t282 = qJDD(4) + t296;
t574 = t282 / 0.2e1;
t265 = t340 * t353 + t357 * t466;
t171 = -qJD(4) * t265 - t297 * t353 + t339 * t357;
t587 = t171 / 0.2e1;
t165 = qJDD(5) - t171;
t589 = t165 / 0.2e1;
t352 = sin(qJ(5));
t356 = cos(qJ(5));
t496 = qJD(2) * t358;
t461 = t348 * t496;
t396 = -qJD(4) + t461;
t214 = t356 * t265 - t352 * t396;
t67 = qJD(5) * t214 + t352 * t170 - t356 * t282;
t597 = t67 / 0.2e1;
t598 = -t67 / 0.2e1;
t213 = t352 * t265 + t356 * t396;
t491 = qJD(5) * t213;
t66 = t356 * t170 + t352 * t282 - t491;
t599 = t66 / 0.2e1;
t643 = Ifges(6,3) + Ifges(7,2);
t645 = Ifges(6,5) + Ifges(7,4);
t255 = qJD(5) - t264;
t355 = sin(qJ(2));
t359 = cos(qJ(2));
t349 = sin(pkin(6));
t499 = qJD(1) * t349;
t456 = qJD(2) * t499;
t484 = qJDD(1) * t349;
t298 = -t355 * t456 + t359 * t484;
t273 = qJDD(2) * pkin(2) + t298;
t467 = t355 * t499;
t310 = pkin(9) * t497 + t467;
t351 = cos(pkin(6));
t483 = qJDD(1) * t351;
t455 = t348 * t483;
t494 = qJD(3) * t358;
t299 = t355 * t484 + t359 * t456;
t322 = qJD(2) * pkin(2) + t359 * t499;
t498 = qJD(1) * t351;
t378 = t322 * t350 + t348 * t498;
t652 = pkin(9) * qJDD(2) * t348 + qJD(3) * t378 + t299;
t78 = t358 * (t273 * t350 + t455) - t310 * t494 - t652 * t354;
t65 = -pkin(3) * t339 - t78;
t29 = -pkin(4) * t171 - pkin(11) * t170 + t65;
t488 = qJD(5) * t356;
t490 = qJD(5) * t352;
t526 = t310 * t358;
t193 = t354 * t378 + t526;
t173 = pkin(10) * t340 + t193;
t336 = t350 * t498;
t421 = -pkin(3) * t358 - pkin(10) * t354;
t208 = t336 + (qJD(2) * t421 - t322) * t348;
t93 = t357 * t173 + t353 * t208;
t82 = -pkin(11) * t396 + t93;
t223 = -t273 * t348 + t350 * t483;
t124 = pkin(3) * t296 - pkin(10) * t297 + t223;
t492 = qJD(4) * t357;
t493 = qJD(4) * t353;
t495 = qJD(3) * t354;
t516 = t350 * t354;
t77 = t273 * t516 - t310 * t495 + t354 * t455 + t358 * t652;
t64 = pkin(10) * t339 + t77;
t17 = t353 * t124 - t173 * t493 + t208 * t492 + t357 * t64;
t9 = pkin(11) * t282 + t17;
t192 = -t354 * t310 + t358 * t378;
t172 = -pkin(3) * t340 - t192;
t95 = -pkin(4) * t264 - pkin(11) * t265 + t172;
t3 = t352 * t29 + t356 * t9 + t95 * t488 - t490 * t82;
t1 = qJ(6) * t165 + qJD(6) * t255 + t3;
t31 = t352 * t95 + t356 * t82;
t4 = -qJD(5) * t31 + t29 * t356 - t352 * t9;
t2 = -pkin(5) * t165 + qJDD(6) - t4;
t664 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t671 = -0.2e1 * Ifges(5,2) * t587 - 0.2e1 * Ifges(5,6) * t574 + Ifges(6,6) * t598 + Ifges(7,6) * t597 + t643 * t589 + t645 * t599 + t664 - t672;
t646 = Ifges(6,1) + Ifges(7,1);
t661 = Ifges(6,4) * t598 + Ifges(7,5) * t597 + t589 * t645 + t599 * t646;
t667 = Ifges(6,6) - Ifges(7,6);
t385 = t348 * (pkin(3) * t354 - pkin(10) * t358);
t290 = qJD(2) * t385;
t126 = t357 * t192 + t353 * t290;
t670 = pkin(10) * t493 + pkin(11) * t466 + t126;
t420 = pkin(4) * t353 - pkin(11) * t357;
t669 = -pkin(10) * qJD(5) * t357 - t322 * t516 - t526 - (t354 * t498 + t420 * t496) * t348 + t420 * qJD(4);
t428 = t353 * t461;
t666 = t428 - t493;
t665 = Ifges(6,4) * t599 + Ifges(6,6) * t589 - t66 * Ifges(7,5) / 0.2e1 - t165 * Ifges(7,6) / 0.2e1 + (Ifges(6,2) + Ifges(7,3)) * t598;
t650 = Ifges(5,1) * t588 + Ifges(5,5) * t574;
t600 = Ifges(5,4) * t587 + t650;
t638 = -t213 * t667 + t214 * t645 + t255 * t643;
t662 = t638 / 0.2e1;
t660 = t77 * mrSges(4,2);
t659 = -Ifges(6,4) + Ifges(7,5);
t642 = t165 * t643 + t645 * t66 - t667 * t67;
t210 = Ifges(6,4) * t213;
t546 = Ifges(7,5) * t213;
t637 = t214 * t646 + t255 * t645 - t210 + t546;
t569 = pkin(4) * t357;
t329 = -pkin(11) * t353 - pkin(3) - t569;
t630 = t329 * t488 + t352 * t669 - t356 * t670;
t629 = -t329 * t490 + t352 * t670 + t356 * t669;
t519 = t348 * t358;
t308 = pkin(2) * t516 + pkin(9) * t519;
t278 = pkin(10) * t350 + t308;
t279 = (-pkin(2) + t421) * t348;
t291 = qJD(3) * t385;
t522 = t348 * t354;
t341 = pkin(9) * t522;
t515 = t350 * t358;
t307 = pkin(2) * t515 - t341;
t292 = t307 * qJD(3);
t106 = -t278 * t493 + t279 * t492 + t353 * t291 + t357 * t292;
t506 = t358 * t359;
t511 = t354 * t355;
t381 = -t350 * t511 + t506;
t267 = t381 * t349;
t245 = qJD(1) * t267;
t430 = t348 * t467;
t207 = t245 * t357 + t353 * t430;
t658 = -t207 + t106;
t657 = -t245 + t292;
t435 = mrSges(4,3) * t466;
t621 = -mrSges(4,1) * t340 - mrSges(5,1) * t264 + mrSges(5,2) * t265 + t435;
t509 = t355 * t358;
t510 = t354 * t359;
t383 = t350 * t509 + t510;
t266 = t383 * t349;
t244 = qJD(1) * t266;
t293 = t308 * qJD(3);
t656 = t293 - t244;
t18 = t124 * t357 - t173 * t492 - t208 * t493 - t353 * t64;
t655 = -t18 * mrSges(5,1) + t17 * mrSges(5,2);
t550 = Ifges(6,4) * t214;
t89 = -Ifges(6,2) * t213 + Ifges(6,6) * t255 + t550;
t594 = -t89 / 0.2e1;
t209 = Ifges(7,5) * t214;
t86 = Ifges(7,6) * t255 + Ifges(7,3) * t213 + t209;
t654 = t594 + t86 / 0.2e1;
t92 = -t353 * t173 + t357 * t208;
t81 = pkin(4) * t396 - t92;
t32 = t213 * pkin(5) - t214 * qJ(6) + t81;
t651 = -mrSges(6,1) * t81 - mrSges(7,1) * t32;
t592 = -m(6) - m(7);
t649 = -t396 / 0.2e1;
t648 = -mrSges(5,3) + mrSges(4,2);
t647 = mrSges(6,3) + mrSges(7,2);
t33 = -mrSges(7,2) * t67 + mrSges(7,3) * t165;
t36 = -mrSges(6,2) * t165 - mrSges(6,3) * t67;
t640 = t33 + t36;
t34 = mrSges(6,1) * t165 - mrSges(6,3) * t66;
t35 = -t165 * mrSges(7,1) + t66 * mrSges(7,2);
t639 = t35 - t34;
t636 = Ifges(5,5) * t396;
t635 = Ifges(5,6) * t396;
t115 = mrSges(5,1) * t282 - mrSges(5,3) * t170;
t25 = mrSges(6,1) * t67 + mrSges(6,2) * t66;
t634 = -t115 + t25;
t633 = -qJ(6) * t666 - qJD(6) * t357 + t630;
t632 = pkin(5) * t666 - t629;
t125 = -t353 * t192 + t290 * t357;
t113 = -pkin(4) * t466 - t125;
t513 = t352 * t357;
t249 = -t356 * t466 + t461 * t513;
t507 = t357 * t358;
t250 = (t352 * t354 + t356 * t507) * t497;
t397 = pkin(5) * t352 - qJ(6) * t356;
t388 = pkin(10) + t397;
t398 = pkin(5) * t356 + qJ(6) * t352;
t631 = -pkin(5) * t249 + qJ(6) * t250 - t113 + (qJD(5) * t398 - qJD(6) * t356) * t353 + t388 * t492;
t628 = -qJD(6) * t352 + t255 * t397 - t93;
t418 = -mrSges(5,1) * t357 + mrSges(5,2) * t353;
t627 = -t418 + mrSges(4,1);
t30 = -t352 * t82 + t356 * t95;
t626 = qJD(6) - t30;
t277 = t341 + (-pkin(2) * t358 - pkin(3)) * t350;
t304 = -t357 * t350 + t353 * t522;
t520 = t348 * t357;
t305 = t350 * t353 + t354 * t520;
t180 = pkin(4) * t304 - pkin(11) * t305 + t277;
t195 = t357 * t278 + t353 * t279;
t182 = -pkin(11) * t519 + t195;
t625 = t352 * t180 + t356 * t182;
t112 = mrSges(6,1) * t213 + mrSges(6,2) * t214;
t537 = t265 * mrSges(5,3);
t218 = -mrSges(5,1) * t396 - t537;
t624 = t218 - t112;
t462 = t352 * t492;
t369 = t353 * t488 + t462;
t623 = t249 - t369;
t489 = qJD(5) * t353;
t370 = -t352 * t489 + t356 * t492;
t622 = t250 - t370;
t414 = mrSges(7,1) * t352 - mrSges(7,3) * t356;
t416 = mrSges(6,1) * t352 + mrSges(6,2) * t356;
t620 = t32 * t414 + t81 * t416;
t619 = t352 * t645 + t356 * t667;
t618 = -t352 * t667 + t356 * t645;
t544 = Ifges(7,5) * t356;
t548 = Ifges(6,4) * t356;
t617 = t352 * t646 - t544 + t548;
t545 = Ifges(7,5) * t352;
t549 = Ifges(6,4) * t352;
t616 = t356 * t646 + t545 - t549;
t615 = t350 * t506 - t511;
t527 = t264 * t356;
t614 = t488 - t527;
t528 = t264 * t352;
t613 = -t490 + t528;
t612 = t17 * t357 - t18 * t353;
t611 = t3 * t356 - t352 * t4;
t610 = t1 * t356 + t2 * t352;
t609 = t89 / 0.2e1 - t86 / 0.2e1;
t415 = -t356 * mrSges(7,1) - t352 * mrSges(7,3);
t417 = mrSges(6,1) * t356 - mrSges(6,2) * t352;
t608 = -m(7) * t398 - mrSges(5,1) + t415 - t417;
t463 = t348 * t494;
t231 = -qJD(4) * t304 + t357 * t463;
t232 = qJD(4) * t305 + t353 * t463;
t117 = pkin(4) * t232 - pkin(11) * t231 + t293;
t464 = t348 * t495;
t96 = pkin(11) * t464 + t106;
t20 = -qJD(5) * t625 + t117 * t356 - t352 * t96;
t436 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t431 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t533 = sin(pkin(12));
t437 = t533 * t359;
t534 = cos(pkin(12));
t440 = t534 * t355;
t302 = t351 * t440 + t437;
t438 = t533 * t355;
t439 = t534 * t359;
t365 = -t351 * t439 + t438;
t363 = t365 * t354;
t442 = t349 * t534;
t423 = t348 * t442;
t175 = t302 * t358 - t350 * t363 - t354 * t423;
t227 = t348 * t365 - t350 * t442;
t103 = t175 * t357 + t227 * t353;
t303 = -t351 * t438 + t439;
t366 = t351 * t437 + t440;
t441 = t349 * t533;
t422 = t348 * t441;
t177 = t303 * t358 + (-t350 * t366 + t422) * t354;
t228 = t348 * t366 + t350 * t441;
t105 = t177 * t357 + t228 * t353;
t382 = t350 * t510 + t509;
t226 = t349 * t382 + t351 * t522;
t517 = t350 * t351;
t518 = t349 * t359;
t384 = -t348 * t518 + t517;
t179 = t226 * t357 + t353 * t384;
t606 = -g(1) * t105 - g(2) * t103 - g(3) * t179;
t26 = -pkin(5) * t255 + t626;
t27 = qJ(6) * t255 + t31;
t604 = -t30 * mrSges(6,1) + t26 * mrSges(7,1) + t31 * mrSges(6,2) - t27 * mrSges(7,3);
t579 = -t255 / 0.2e1;
t584 = -t214 / 0.2e1;
t585 = t213 / 0.2e1;
t586 = -t213 / 0.2e1;
t603 = Ifges(6,6) * t585 + Ifges(7,6) * t586 + t579 * t643 + t584 * t645 + t604;
t360 = qJD(2) ^ 2;
t553 = Ifges(5,4) * t265;
t152 = Ifges(5,2) * t264 + t553 - t635;
t590 = -t152 / 0.2e1;
t583 = t214 / 0.2e1;
t578 = t255 / 0.2e1;
t577 = -t264 / 0.2e1;
t576 = -t265 / 0.2e1;
t575 = t265 / 0.2e1;
t572 = t350 / 0.2e1;
t559 = qJD(4) / 0.2e1;
t557 = mrSges(6,3) * t213;
t556 = mrSges(6,3) * t214;
t555 = Ifges(4,4) * t354;
t554 = Ifges(4,4) * t358;
t552 = Ifges(5,4) * t353;
t551 = Ifges(5,4) * t357;
t547 = Ifges(5,5) * t265;
t543 = Ifges(5,6) * t264;
t542 = Ifges(5,3) * t354;
t10 = -pkin(4) * t282 - t18;
t541 = t10 * t353;
t538 = t264 * mrSges(5,3);
t536 = t340 * Ifges(4,5);
t535 = t340 * Ifges(4,6);
t191 = pkin(4) * t265 - pkin(11) * t264;
t42 = t352 * t191 + t356 * t92;
t362 = t365 * t358;
t174 = t302 * t354 + t350 * t362 + t358 * t423;
t532 = t174 * t353;
t364 = t366 * t358;
t176 = t303 * t354 + t350 * t364 - t358 * t422;
t531 = t176 * t353;
t211 = mrSges(4,1) * t296 + mrSges(4,2) * t297;
t530 = t211 * t348;
t225 = -t349 * t615 - t351 * t519;
t529 = t225 * t353;
t524 = t329 * t356;
t523 = t348 * t353;
t521 = t348 * t355;
t514 = t352 * t353;
t512 = t353 * t356;
t508 = t356 * t357;
t138 = -mrSges(7,2) * t213 + mrSges(7,3) * t255;
t139 = -mrSges(6,2) * t255 - t557;
t505 = t138 + t139;
t140 = mrSges(6,1) * t255 - t556;
t141 = -mrSges(7,1) * t255 + mrSges(7,2) * t214;
t504 = t141 - t140;
t203 = -t302 * t516 - t362;
t294 = t365 * pkin(2);
t503 = t203 * pkin(3) - t294;
t205 = -t303 * t516 - t364;
t295 = t366 * pkin(2);
t502 = t205 * pkin(3) - t295;
t477 = t349 * t521;
t500 = pkin(2) * t518 + pkin(9) * t477;
t275 = pkin(10) * t508 + t352 * t329;
t486 = -m(5) + t592;
t481 = pkin(10) * t492;
t480 = pkin(11) * t490;
t479 = pkin(11) * t488;
t478 = -m(4) + t486;
t476 = t352 * t519;
t472 = Ifges(5,5) * t170 + Ifges(5,6) * t171 + Ifges(5,3) * t282;
t111 = mrSges(7,1) * t213 - mrSges(7,3) * t214;
t471 = -t111 + t624;
t470 = t267 * pkin(3) + t500;
t469 = Ifges(4,5) * t297 - Ifges(4,6) * t296 + Ifges(4,3) * t339;
t451 = t492 / 0.2e1;
t448 = -t489 / 0.2e1;
t447 = t488 / 0.2e1;
t446 = -t174 * pkin(3) + pkin(10) * t175;
t445 = -t176 * pkin(3) + pkin(10) * t177;
t444 = -t225 * pkin(3) + pkin(10) * t226;
t194 = -t353 * t278 + t279 * t357;
t434 = mrSges(4,3) * t461;
t429 = qJD(2) * t477;
t181 = pkin(4) * t519 - t194;
t413 = Ifges(5,1) * t357 - t552;
t408 = -Ifges(5,2) * t353 + t551;
t407 = -Ifges(6,2) * t352 + t548;
t406 = Ifges(6,2) * t356 + t549;
t403 = Ifges(5,5) * t357 - Ifges(5,6) * t353;
t400 = Ifges(7,3) * t352 + t544;
t399 = -Ifges(7,3) * t356 + t545;
t41 = t191 * t356 - t352 * t92;
t101 = t179 * t356 + t225 * t352;
t100 = t179 * t352 - t225 * t356;
t79 = t180 * t356 - t182 * t352;
t107 = -t278 * t492 - t279 * t493 + t291 * t357 - t353 * t292;
t387 = pkin(10) * t486 + t648;
t386 = pkin(11) * t592 + mrSges(5,2) - t647;
t233 = t305 * t352 + t356 * t519;
t376 = t172 * (mrSges(5,1) * t353 + mrSges(5,2) * t357);
t375 = t354 * (Ifges(4,1) * t358 - t555);
t374 = (-mrSges(4,1) * t358 + mrSges(4,2) * t354) * t348;
t373 = (Ifges(4,2) * t358 + t555) * t348;
t19 = t352 * t117 + t180 * t488 - t182 * t490 + t356 * t96;
t367 = mrSges(3,2) + (pkin(9) * t478 - mrSges(4,3)) * t348;
t97 = -pkin(4) * t464 - t107;
t178 = t226 * t353 - t357 * t384;
t334 = Ifges(4,4) * t461;
t324 = -pkin(4) - t398;
t289 = qJD(2) * t374;
t288 = -mrSges(4,2) * t340 + t434;
t286 = t388 * t353;
t274 = -pkin(10) * t513 + t524;
t253 = Ifges(5,4) * t264;
t252 = -t322 * t348 + t336;
t242 = -t524 + (pkin(10) * t352 + pkin(5)) * t357;
t241 = -qJ(6) * t357 + t275;
t237 = Ifges(4,1) * t466 + t334 + t536;
t236 = qJD(2) * t373 + t535;
t234 = t305 * t356 - t476;
t230 = mrSges(4,1) * t339 - mrSges(4,3) * t297;
t229 = -mrSges(4,2) * t339 - mrSges(4,3) * t296;
t217 = mrSges(5,2) * t396 + t538;
t216 = t267 * t357 + t353 * t477;
t206 = t245 * t353 - t357 * t430;
t204 = t303 * t515 - t354 * t366;
t202 = t302 * t515 - t363;
t156 = t351 * t463 + (t381 * qJD(2) + qJD(3) * t615) * t349;
t155 = t351 * t464 + (qJD(2) * t383 + qJD(3) * t382) * t349;
t153 = Ifges(5,1) * t265 + t253 - t636;
t151 = -Ifges(5,3) * t396 + t543 + t547;
t134 = t205 * t357 + t303 * t523;
t132 = t203 * t357 + t302 * t523;
t122 = t207 * t356 + t244 * t352;
t121 = t207 * t352 - t356 * t244;
t120 = -qJD(5) * t476 + t231 * t352 + t305 * t488 - t356 * t464;
t119 = -qJD(5) * t233 + t231 * t356 + t352 * t464;
t116 = -mrSges(5,2) * t282 + mrSges(5,3) * t171;
t110 = pkin(5) * t214 + qJ(6) * t213;
t104 = -t177 * t353 + t228 * t357;
t102 = -t175 * t353 + t227 * t357;
t83 = pkin(5) * t233 - qJ(6) * t234 + t181;
t76 = -mrSges(5,1) * t171 + mrSges(5,2) * t170;
t55 = -qJD(4) * t178 + t156 * t357 + t353 * t429;
t54 = qJD(4) * t179 + t156 * t353 - t357 * t429;
t53 = -pkin(5) * t304 - t79;
t52 = qJ(6) * t304 + t625;
t48 = t105 * t352 - t176 * t356;
t46 = t103 * t352 - t174 * t356;
t38 = -pkin(5) * t265 - t41;
t37 = qJ(6) * t265 + t42;
t24 = mrSges(7,1) * t67 - mrSges(7,3) * t66;
t23 = pkin(5) * t120 - qJ(6) * t119 - qJD(6) * t234 + t97;
t22 = -qJD(5) * t100 + t155 * t352 + t356 * t55;
t21 = qJD(5) * t101 - t155 * t356 + t352 * t55;
t7 = -pkin(5) * t232 - t20;
t6 = qJ(6) * t232 + qJD(6) * t304 + t19;
t5 = pkin(5) * t67 - qJ(6) * t66 - qJD(6) * t214 + t10;
t8 = [t211 * t517 + m(2) * qJDD(1) + t179 * t116 + t156 * t288 + t55 * t217 + t226 * t229 + (-t230 + t76) * t225 + t505 * t22 + t504 * t21 + t621 * t155 + t640 * t101 + t639 * t100 - t471 * t54 + (t24 + t634) * t178 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t360 - t530) * t359 + (-mrSges(3,1) * t360 - mrSges(3,2) * qJDD(2) + t289 * t497) * t355) * t349 + (-m(2) - m(3) + t478) * g(3) + m(5) * (t155 * t172 + t17 * t179 - t178 * t18 + t225 * t65 - t54 * t92 + t55 * t93) + m(7) * (t1 * t101 + t100 * t2 + t178 * t5 + t21 * t26 + t22 * t27 + t32 * t54) + m(6) * (t10 * t178 - t100 * t4 + t101 * t3 - t21 * t30 + t22 * t31 + t54 * t81) + m(3) * (qJDD(1) * t351 ^ 2 + (t298 * t359 + t299 * t355) * t349) + m(4) * (-t192 * t155 + t193 * t156 + t223 * t384 - t78 * t225 + t77 * t226 + t252 * t429); (-t231 * t92 - t232 * t93) * mrSges(5,3) + (t10 * t181 + t3 * t625 + t4 * t79 + (-t206 + t97) * t81 + (-t122 + t19) * t31 + (t121 + t20) * t30) * m(6) + t625 * t36 + (t119 * t646 + t232 * t645) * t583 - t350 * t660 + t232 * t662 + (Ifges(5,5) * t231 - Ifges(5,6) * t232) * t649 + (t119 * t645 + t232 * t643) * t578 + (t1 * t52 + t2 * t53 + t5 * t83 + (-t206 + t23) * t32 + (-t122 + t6) * t27 + (-t121 + t7) * t26) * m(7) + (Ifges(6,4) * t119 + Ifges(6,6) * t232) * t586 + (-t289 * t521 + (mrSges(3,1) * t355 + mrSges(3,2) * t359) * qJD(2)) * t499 + t78 * (mrSges(4,1) * t350 - mrSges(4,3) * t522) + (mrSges(6,2) * t10 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t5 + 0.2e1 * t661) * t234 - t504 * t121 - t505 * t122 + (t65 * mrSges(5,2) - t18 * mrSges(5,3) + 0.2e1 * t600) * t305 + t637 * t119 / 0.2e1 + (Ifges(7,5) * t119 + Ifges(7,6) * t232) * t585 + t223 * t374 + (-t472 / 0.2e1 + t77 * mrSges(4,3) - Ifges(5,3) * t574 - Ifges(5,6) * t587 - Ifges(5,5) * t588 + t655) * t519 + t621 * t656 + (-pkin(2) * t223 * t348 + g(1) * t295 + g(2) * t294 - g(3) * t500 - t192 * t656 + t193 * t657 - t252 * t430 + t307 * t78 + t308 * t77) * m(4) + t657 * t288 + t658 * t217 + (-g(1) * t502 - g(2) * t503 - g(3) * t470 + t17 * t195 + t18 * t194 + t277 * t65 + t658 * t93 + (t107 + t206) * t92 + t656 * t172) * m(5) + (t365 * mrSges(3,1) - t203 * mrSges(4,1) - t132 * mrSges(5,1) - t436 * (t132 * t356 + t202 * t352) + t431 * (t132 * t352 - t202 * t356) + t367 * t302 + t387 * t202 + t386 * (t203 * t353 - t302 * t520) + t592 * (t132 * pkin(4) + t503)) * g(2) + (t366 * mrSges(3,1) - t205 * mrSges(4,1) - t134 * mrSges(5,1) - t436 * (t134 * t356 + t204 * t352) + t431 * (t134 * t352 - t204 * t356) + t367 * t303 + t387 * t204 + t386 * (t205 * t353 - t303 * t520) + t592 * (t134 * pkin(4) + t502)) * g(1) + (-t267 * mrSges(4,1) - t216 * mrSges(5,1) + t387 * t266 - t436 * (t216 * t356 + t266 * t352) + t431 * (t216 * t352 - t266 * t356) + (-mrSges(3,1) * t359 + (-mrSges(4,3) * t348 + mrSges(3,2)) * t355) * t349 + t386 * (t267 * t353 - t357 * t477) + t592 * (t216 * pkin(4) + t470)) * g(3) + (-t119 * t32 + t232 * t27) * mrSges(7,3) + ((-t192 * mrSges(4,3) + t252 * mrSges(4,2) + t536 / 0.2e1 + t237 / 0.2e1) * t358 + (-t193 * mrSges(4,3) + t252 * mrSges(4,1) + Ifges(5,3) * t559 - t535 / 0.2e1 + t547 / 0.2e1 + t543 / 0.2e1 + t92 * mrSges(5,1) - t93 * mrSges(5,2) - t236 / 0.2e1 + t151 / 0.2e1) * t354 + (-t358 * t542 / 0.2e1 + t375 / 0.2e1 + t358 * (-Ifges(4,2) * t354 + t554) / 0.2e1) * t497) * t348 * qJD(3) - pkin(2) * t530 + t469 * t572 + (Ifges(5,1) * t231 - Ifges(5,4) * t232) * t575 + t471 * t206 + (Ifges(4,3) * t572 + (Ifges(4,5) * t354 + Ifges(4,6) * t358) * t348) * t339 - (Ifges(4,6) * t572 + t373) * t296 + (Ifges(4,5) * t572 + (Ifges(4,1) * t354 + t554) * t348) * t297 + (t642 / 0.2e1 - t17 * mrSges(5,3) + t65 * mrSges(5,1) - t672 + t671) * t304 + (-mrSges(7,2) * t27 - mrSges(6,3) * t31 - Ifges(6,2) * t586 + Ifges(7,3) * t585 - t578 * t667 + t583 * t659 - t651 + t654) * t120 + (mrSges(6,1) * t10 + mrSges(7,1) * t5 - mrSges(7,2) * t1 - mrSges(6,3) * t3 - Ifges(6,2) * t598 + Ifges(7,3) * t597 - t589 * t667 + t599 * t659 - t665) * t233 + Ifges(3,3) * qJDD(2) + t307 * t230 + t308 * t229 + t298 * mrSges(3,1) - t299 * mrSges(3,2) + t277 * t76 + t231 * t153 / 0.2e1 + t30 * (mrSges(6,1) * t232 - mrSges(6,3) * t119) + t26 * (-mrSges(7,1) * t232 + mrSges(7,2) * t119) + t172 * (mrSges(5,1) * t232 + mrSges(5,2) * t231) + t107 * t218 + t194 * t115 + t195 * t116 + t181 * t25 + t6 * t138 + t19 * t139 + t20 * t140 + t7 * t141 + t23 * t111 + t97 * t112 + t264 * (Ifges(5,4) * t231 - Ifges(5,2) * t232) / 0.2e1 + t83 * t24 + t79 * t34 + t52 * t33 + t53 * t35 + t232 * t590 + (t119 * t81 - t232 * t31) * mrSges(6,2); -((-Ifges(4,2) * t466 + t153 * t357 + t353 * t638 + t237 + t334) * t358 + t265 * (Ifges(5,5) * t354 + t358 * t413) + t264 * (Ifges(5,6) * t354 + t358 * t408) + t340 * (Ifges(4,5) * t358 - Ifges(4,6) * t354) + t354 * t151) * t497 / 0.2e1 + (t400 * t597 + t407 * t598 + t5 * t414 + t86 * t447 + t589 * t618 + t599 * t616 + t600 + t650) * t353 + t551 * t588 + (-t481 - t125) * t218 + t512 * t661 + t552 * t587 + (t396 * (t358 * t403 + t542) + t354 * t236) * t497 / 0.2e1 + (t264 * t408 + t265 * t413) * t559 - t348 ^ 2 * t360 * t375 / 0.2e1 + (-t93 * mrSges(5,3) - pkin(10) * t217 + t590 - t604 + t662) * t493 + t416 * t541 + (-t93 * (-mrSges(5,3) * t353 * t358 - mrSges(5,2) * t354) - t252 * (mrSges(4,1) * t354 + mrSges(4,2) * t358) - t92 * (mrSges(5,1) * t354 - mrSges(5,3) * t507)) * t497 + (t434 - t288) * t192 + t637 * (t352 * t448 + t356 * t451 - t250 / 0.2e1) + t81 * (mrSges(6,1) * t369 + mrSges(6,2) * t370) + t32 * (mrSges(7,1) * t369 - mrSges(7,3) * t370) + t356 * t89 * t448 + (t352 * t86 + t153) * t451 + t65 * t418 + (t152 / 0.2e1 + t603) * t428 + (pkin(10) * t116 - t671) * t357 + (t481 - t113) * t112 + (-Ifges(6,2) * t585 + Ifges(7,3) * t586 - t579 * t667 + t584 * t659 + t609 + t651) * t249 - t665 * t514 + t286 * t24 + t274 * t34 + t275 * t36 + t241 * t33 + t242 * t35 - t126 * t217 - pkin(3) * t76 + t78 * mrSges(4,1) + t634 * pkin(10) * t353 + (-t125 * t92 - t126 * t93 - pkin(3) * t65 + ((-t353 * t93 - t357 * t92) * qJD(4) + t612) * pkin(10)) * m(5) + (-t492 * t92 + t612) * mrSges(5,3) + t469 + (-t399 * t585 - t406 * t586 - t578 * t619 - t583 * t617) * t489 + (-m(5) * t172 + t435 - t621) * t193 + (-t3 * t514 + t30 * t622 + t31 * t623 - t4 * t512) * mrSges(6,3) + (-t1 * t514 + t2 * t512 - t26 * t622 + t27 * t623) * mrSges(7,2) + t462 * t594 + t629 * t140 + t630 * t139 + (-t113 * t81 + t274 * t4 + t275 * t3 + (t492 * t81 + t541) * pkin(10) + t630 * t31 + t629 * t30) * m(6) + t631 * t111 + t632 * t141 + t633 * t138 + (t1 * t241 + t2 * t242 + t26 * t632 + t27 * t633 + t286 * t5 + t32 * t631) * m(7) - t642 * t357 / 0.2e1 + (-t81 * mrSges(6,2) + t32 * mrSges(7,3) + Ifges(6,4) * t585 + Ifges(7,5) * t586 + t579 * t645 + t584 * t646) * t250 - t660 + (t376 + t403 * t649 + (Ifges(7,6) * t353 + t357 * t400) * t585 + (Ifges(6,6) * t353 + t357 * t407) * t586 + (t353 * t645 + t357 * t616) * t583 + (t353 * t643 + t357 * t618) * t578) * qJD(4) + (-m(5) * t444 + t647 * t529 + t592 * (-pkin(11) * t529 - t225 * t569 + t444) + t648 * t226 + t627 * t225 - t436 * (-t225 * t508 + t226 * t352) + t431 * (-t225 * t513 - t226 * t356)) * g(3) + (-m(5) * t445 + t592 * (-pkin(11) * t531 - t176 * t569 + t445) - t436 * (-t176 * t508 + t177 * t352) + t431 * (-t176 * t513 - t177 * t356) + t647 * t531 + t648 * t177 + t627 * t176) * g(1) + (-m(5) * t446 + t592 * (-pkin(11) * t532 - t174 * t569 + t446) - t436 * (-t174 * t508 + t175 * t352) + t431 * (-t174 * t513 - t175 * t356) + t647 * t532 + t648 * t175 + t627 * t174) * g(2) - t376 * t461; (-t480 - t37) * t138 + (-t407 / 0.2e1 + t400 / 0.2e1) * t491 + (t214 * t616 + t255 * t618) * qJD(5) / 0.2e1 + t352 * t661 + (t447 - t527 / 0.2e1) * t637 - t655 + (t479 - t38) * t141 + (-t479 - t41) * t140 + (-t480 - t42) * t139 + (t538 - t217) * t92 + t654 * t490 + t152 * t575 + t5 * t415 - t10 * t417 + (t253 + t153) * t577 + t324 * t24 - pkin(4) * t25 + t639 * pkin(11) * t352 + t640 * pkin(11) * t356 + t609 * t528 + (-pkin(4) * t10 + ((-t30 * t356 - t31 * t352) * qJD(5) + t611) * pkin(11) - t30 * t41 - t31 * t42) * m(6) + t472 + (t26 * t614 + t27 * t613 + t606 + t610) * mrSges(7,2) + (-t30 * t614 + t31 * t613 + t606 + t611) * mrSges(6,3) + t617 * t599 + t619 * t589 + t620 * qJD(5) + (-m(6) * t81 + t537 + t624) * t93 + t665 * t356 + t399 * t597 + t406 * t598 + t628 * t111 + (t324 * t5 + ((t26 * t356 - t27 * t352) * qJD(5) + t610) * pkin(11) - t26 * t38 - t27 * t37 + t628 * t32) * m(7) + (-t635 / 0.2e1 - Ifges(5,2) * t577 - t172 * mrSges(5,1) + t603) * t265 + (t636 / 0.2e1 + Ifges(5,1) * t576 - t172 * mrSges(5,2) + t407 * t585 + t400 * t586 + t616 * t584 + t618 * t579 - t620) * t264 + (-t553 + t638) * t576 + (mrSges(5,2) * t179 + t592 * (-t178 * pkin(4) + pkin(11) * t179) - t608 * t178) * g(3) + (mrSges(5,2) * t105 + t592 * (t104 * pkin(4) + pkin(11) * t105) + t608 * t104) * g(1) + (mrSges(5,2) * t103 + t592 * (t102 * pkin(4) + pkin(11) * t103) + t608 * t102) * g(2); (-t213 * t646 + t209 - t550 + t86) * t584 + (t213 * t26 + t214 * t27) * mrSges(7,2) + (-t213 * t645 - t214 * t667) * t579 + (-pkin(5) * t2 + qJ(6) * t1 - t110 * t32 - t26 * t31 + t27 * t626) * m(7) + (t431 * (t103 * t356 + t174 * t352) + t436 * t46) * g(2) + (t431 * (t105 * t356 + t176 * t352) + t436 * t48) * g(1) + (t100 * t436 + t101 * t431) * g(3) + t89 * t583 + (-t504 + t556) * t31 + (-t505 - t557) * t30 + t664 + (-Ifges(6,2) * t214 - t210 + t637) * t585 - t32 * (mrSges(7,1) * t214 + mrSges(7,3) * t213) - t81 * (mrSges(6,1) * t214 - mrSges(6,2) * t213) + qJD(6) * t138 - t110 * t111 + qJ(6) * t33 - pkin(5) * t35 + (Ifges(7,3) * t214 - t546) * t586 + t642; t214 * t111 - t255 * t138 + (-g(1) * t48 - g(2) * t46 - g(3) * t100 + t32 * t214 - t255 * t27 + t2) * m(7) + t35;];
tau  = t8;
