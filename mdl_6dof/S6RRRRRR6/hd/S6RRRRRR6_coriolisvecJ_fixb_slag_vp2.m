% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:48
% EndTime: 2019-03-10 04:13:16
% DurationCPUTime: 46.96s
% Computational Cost: add. (40472->1022), mult. (101268->1434), div. (0->0), fcn. (81604->12), ass. (0->447)
t424 = cos(qJ(2));
t419 = sin(qJ(2));
t413 = sin(pkin(6));
t488 = qJD(1) * t413;
t469 = t419 * t488;
t414 = cos(pkin(6));
t487 = qJD(1) * t414;
t477 = pkin(1) * t487;
t364 = -pkin(8) * t469 + t424 * t477;
t437 = (pkin(2) * t419 - pkin(9) * t424) * t413;
t365 = qJD(1) * t437;
t418 = sin(qJ(3));
t423 = cos(qJ(3));
t289 = -t364 * t418 + t365 * t423;
t590 = -pkin(10) - pkin(9);
t471 = qJD(3) * t590;
t664 = -(-pkin(10) * t423 * t424 + pkin(3) * t419) * t488 - t289 + t423 * t471;
t290 = t364 * t423 + t365 * t418;
t468 = t424 * t488;
t455 = t418 * t468;
t663 = -pkin(10) * t455 - t418 * t471 + t290;
t417 = sin(qJ(4));
t422 = cos(qJ(4));
t378 = t417 * t418 - t422 * t423;
t602 = qJD(3) + qJD(4);
t315 = t602 * t378;
t328 = t378 * t468;
t662 = t315 - t328;
t380 = t417 * t423 + t418 * t422;
t316 = t602 * t380;
t327 = t380 * t468;
t489 = t316 - t327;
t395 = t590 * t418;
t397 = t590 * t423;
t605 = t395 * t422 + t397 * t417;
t622 = qJD(4) * t605 + t417 * t664 - t422 * t663;
t367 = pkin(8) * t468 + t419 * t477;
t319 = pkin(3) * t455 + t367;
t484 = qJD(3) * t418;
t661 = pkin(3) * t484 - t319;
t660 = pkin(11) * t469 - t622;
t659 = pkin(4) * t489 + pkin(11) * t662 + t661;
t416 = sin(qJ(5));
t421 = cos(qJ(5));
t292 = t328 * t416 + t421 * t469;
t500 = t315 * t416;
t457 = t292 - t500;
t479 = qJD(5) * t421;
t658 = t380 * t479 + t457;
t408 = pkin(3) * t417 + pkin(11);
t541 = -pkin(12) - t408;
t460 = qJD(5) * t541;
t481 = qJD(4) * t422;
t475 = pkin(3) * t481;
t402 = qJD(2) + t487;
t347 = t402 * t423 - t418 * t469;
t348 = t402 * t418 + t423 * t469;
t606 = t347 * t422 - t348 * t417;
t628 = t606 * t416;
t647 = pkin(12) * t628;
t326 = pkin(9) * t402 + t367;
t360 = (-pkin(2) * t424 - pkin(9) * t419 - pkin(1)) * t413;
t342 = qJD(1) * t360;
t269 = t326 * t423 + t342 * t418;
t235 = pkin(10) * t347 + t269;
t229 = t417 * t235;
t268 = -t326 * t418 + t342 * t423;
t234 = -pkin(10) * t348 + t268;
t143 = t234 * t422 - t229;
t440 = t347 * t417 + t348 * t422;
t212 = pkin(4) * t440 - pkin(11) * t606;
t182 = pkin(3) * t348 + t212;
t91 = t143 * t421 + t182 * t416;
t657 = t416 * t460 + t421 * t475 + t647 - t91;
t412 = t421 * pkin(12);
t619 = pkin(5) * t440 - t412 * t606;
t90 = -t143 * t416 + t182 * t421;
t656 = -t416 * t475 + t421 * t460 - t619 - t90;
t411 = -pkin(3) * t423 - pkin(2);
t311 = pkin(4) * t378 - pkin(11) * t380 + t411;
t332 = t395 * t417 - t397 * t422;
t480 = qJD(5) * t416;
t626 = t311 * t479 - t332 * t480 + t416 * t659 - t421 * t660;
t655 = t416 * t660 + t421 * t659;
t485 = qJD(2) * t424;
t467 = t418 * t485;
t483 = qJD(3) * t423;
t306 = -t402 * t484 + (-t419 * t483 - t467) * t488;
t495 = t413 * t424;
t374 = pkin(1) * t414 * t419 + pkin(8) * t495;
t369 = t374 * qJD(2);
t357 = qJD(1) * t369;
t267 = -pkin(3) * t306 + t357;
t466 = t423 * t485;
t305 = t402 * t483 + (-t419 * t484 + t466) * t488;
t177 = qJD(4) * t606 + t305 * t422 + t306 * t417;
t526 = t177 * Ifges(5,4);
t366 = qJD(2) * t437;
t355 = qJD(1) * t366;
t496 = t413 * t419;
t403 = pkin(8) * t496;
t546 = pkin(1) * t424;
t373 = t414 * t546 - t403;
t368 = t373 * qJD(2);
t356 = qJD(1) * t368;
t206 = -qJD(3) * t269 + t355 * t423 - t356 * t418;
t486 = qJD(2) * t413;
t462 = qJD(1) * t486;
t454 = t419 * t462;
t148 = pkin(3) * t454 - pkin(10) * t305 + t206;
t205 = -t326 * t484 + t342 * t483 + t355 * t418 + t356 * t423;
t158 = pkin(10) * t306 + t205;
t393 = qJD(3) - t468;
t218 = pkin(3) * t393 + t234;
t482 = qJD(4) * t417;
t58 = t148 * t417 + t158 * t422 + t218 * t481 - t235 * t482;
t230 = t422 * t235;
t137 = t218 * t417 + t230;
t387 = qJD(4) + t393;
t132 = pkin(11) * t387 + t137;
t325 = -pkin(2) * t402 - t364;
t282 = -pkin(3) * t347 + t325;
t149 = -pkin(4) * t606 - pkin(11) * t440 + t282;
t55 = pkin(11) * t454 + t58;
t178 = qJD(4) * t440 + t305 * t417 - t306 * t422;
t89 = pkin(4) * t178 - pkin(11) * t177 + t267;
t16 = -t132 * t480 + t149 * t479 + t416 * t89 + t421 * t55;
t80 = t132 * t421 + t149 * t416;
t17 = -qJD(5) * t80 - t416 * t55 + t421 * t89;
t420 = cos(qJ(6));
t415 = sin(qJ(6));
t246 = t387 * t421 - t416 * t440;
t65 = pkin(12) * t246 + t80;
t507 = t415 * t65;
t277 = qJD(5) - t606;
t247 = t387 * t416 + t421 * t440;
t79 = -t132 * t416 + t149 * t421;
t64 = -pkin(12) * t247 + t79;
t60 = pkin(5) * t277 + t64;
t23 = t420 * t60 - t507;
t116 = qJD(5) * t246 + t177 * t421 + t416 * t454;
t7 = pkin(5) * t178 - pkin(12) * t116 + t17;
t117 = -qJD(5) * t247 - t177 * t416 + t421 * t454;
t8 = pkin(12) * t117 + t16;
t3 = qJD(6) * t23 + t415 * t7 + t420 * t8;
t505 = t420 * t65;
t24 = t415 * t60 + t505;
t4 = -qJD(6) * t24 - t415 * t8 + t420 * t7;
t634 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t615 = t17 * mrSges(6,1) - t16 * mrSges(6,2) + t634;
t654 = -t615 - t267 * mrSges(5,1) + t58 * mrSges(5,3) + t526 / 0.2e1;
t438 = t415 * t416 - t420 * t421;
t601 = qJD(5) + qJD(6);
t313 = t601 * t438;
t629 = t438 * t606;
t653 = -t629 + t313;
t379 = t415 * t421 + t416 * t420;
t314 = t601 * t379;
t630 = t379 * t606;
t652 = t630 - t314;
t293 = -t328 * t421 + t416 * t469;
t320 = t421 * t332;
t499 = t315 * t421;
t651 = pkin(12) * t293 + pkin(12) * t499 + (-t320 + (pkin(12) * t380 - t311) * t416) * qJD(5) + t655 + t489 * pkin(5);
t650 = pkin(12) * t658 - t626;
t623 = -qJD(4) * t332 + t417 * t663 + t422 * t664;
t648 = pkin(5) * t628;
t646 = t178 * Ifges(5,2);
t375 = t541 * t416;
t376 = t408 * t421 + t412;
t310 = t375 * t415 + t376 * t420;
t645 = -qJD(6) * t310 - t415 * t657 + t420 * t656;
t309 = t375 * t420 - t376 * t415;
t644 = qJD(6) * t309 + t415 * t656 + t420 * t657;
t589 = -pkin(12) - pkin(11);
t470 = qJD(5) * t589;
t136 = t218 * t422 - t229;
t95 = t136 * t421 + t212 * t416;
t643 = t416 * t470 + t647 - t95;
t94 = -t136 * t416 + t212 * t421;
t642 = t421 * t470 - t619 - t94;
t59 = t148 * t422 - t158 * t417 - t218 * t482 - t235 * t481;
t56 = -pkin(4) * t454 - t59;
t63 = -mrSges(6,1) * t117 + mrSges(6,2) * t116;
t641 = -m(6) * t56 - t63;
t142 = t234 * t417 + t230;
t474 = pkin(5) * t480;
t640 = pkin(3) * t482 - t142 + t474 - t648;
t624 = pkin(4) * t469 - t623;
t131 = -pkin(4) * t387 - t136;
t106 = -pkin(5) * t246 + t131;
t162 = t246 * t415 + t247 * t420;
t458 = t246 * t420 - t247 * t415;
t533 = Ifges(7,4) * t162;
t273 = qJD(6) + t277;
t571 = -t273 / 0.2e1;
t579 = -t162 / 0.2e1;
t637 = (Ifges(7,5) * t458 - Ifges(7,6) * t162) * t571 + (t162 * t24 + t23 * t458) * mrSges(7,3) - t106 * (mrSges(7,1) * t162 + mrSges(7,2) * t458) + (Ifges(7,1) * t458 - t533) * t579;
t511 = t387 * Ifges(5,6);
t520 = t606 * Ifges(5,2);
t537 = Ifges(5,4) * t440;
t196 = t511 + t520 + t537;
t636 = t137 * mrSges(5,3) + t24 * mrSges(7,2) + t80 * mrSges(6,2) + t196 / 0.2e1 - t23 * mrSges(7,1) - t282 * mrSges(5,1) - t79 * mrSges(6,1);
t536 = Ifges(6,4) * t247;
t129 = Ifges(6,2) * t246 + Ifges(6,6) * t277 + t536;
t36 = -pkin(5) * t117 + t56;
t446 = Ifges(6,5) * t421 - Ifges(6,6) * t416;
t431 = t277 * t446;
t535 = Ifges(6,4) * t416;
t450 = Ifges(6,1) * t421 - t535;
t432 = t247 * t450;
t534 = Ifges(6,4) * t421;
t448 = -Ifges(6,2) * t416 + t534;
t433 = t246 * t448;
t451 = mrSges(6,1) * t416 + mrSges(6,2) * t421;
t434 = t131 * t451;
t445 = Ifges(6,5) * t416 + Ifges(6,6) * t421;
t447 = Ifges(6,2) * t421 + t535;
t449 = Ifges(6,1) * t416 + t534;
t452 = mrSges(6,1) * t421 - mrSges(6,2) * t416;
t473 = Ifges(5,5) * t177 - Ifges(5,6) * t178 + Ifges(5,3) * t454;
t49 = Ifges(6,4) * t116 + Ifges(6,2) * t117 + Ifges(6,6) * t178;
t50 = Ifges(6,1) * t116 + Ifges(6,4) * t117 + Ifges(6,5) * t178;
t548 = t421 / 0.2e1;
t550 = t416 / 0.2e1;
t570 = t273 / 0.2e1;
t577 = t178 / 0.2e1;
t578 = t162 / 0.2e1;
t580 = t458 / 0.2e1;
t581 = -t458 / 0.2e1;
t244 = Ifges(6,4) * t246;
t130 = Ifges(6,1) * t247 + Ifges(6,5) * t277 + t244;
t582 = t130 / 0.2e1;
t583 = t117 / 0.2e1;
t584 = t116 / 0.2e1;
t157 = Ifges(7,4) * t458;
t83 = Ifges(7,1) * t162 + Ifges(7,5) * t273 + t157;
t591 = t83 / 0.2e1;
t592 = -t83 / 0.2e1;
t82 = Ifges(7,2) * t458 + Ifges(7,6) * t273 + t533;
t593 = t82 / 0.2e1;
t594 = -t82 / 0.2e1;
t44 = -qJD(6) * t162 - t116 * t415 + t117 * t420;
t595 = t44 / 0.2e1;
t43 = qJD(6) * t458 + t116 * t420 + t117 * t415;
t596 = t43 / 0.2e1;
t597 = Ifges(7,1) * t596 + Ifges(7,4) * t595 + Ifges(7,5) * t577;
t598 = Ifges(7,4) * t596 + Ifges(7,2) * t595 + Ifges(7,6) * t577;
t621 = t16 * t421 - t17 * t416;
t635 = (-Ifges(7,5) * t313 - Ifges(7,6) * t314) * t570 + (-Ifges(7,1) * t313 - Ifges(7,4) * t314) * t578 + (-Ifges(7,4) * t313 - Ifges(7,2) * t314) * t580 + (t23 * t653 + t24 * t652 - t3 * t438 - t4 * t379) * mrSges(7,3) + (-mrSges(7,1) * t652 - mrSges(7,2) * t653) * t106 + (t433 + t432 + t431) * qJD(5) / 0.2e1 - t56 * t452 - t630 * t594 + (-Ifges(7,5) * t629 - Ifges(7,6) * t630) * t571 + (-Ifges(7,1) * t629 - Ifges(7,4) * t630) * t579 + (-Ifges(7,4) * t629 - Ifges(7,2) * t630) * t581 - t629 * t592 + t473 + (Ifges(7,1) * t379 - Ifges(7,4) * t438) * t596 + t36 * (mrSges(7,1) * t438 + mrSges(7,2) * t379) + (Ifges(7,4) * t379 - Ifges(7,2) * t438) * t595 + (Ifges(7,5) * t379 - Ifges(7,6) * t438 + t445) * t577 - t438 * t598 + t621 * mrSges(6,3) - t58 * mrSges(5,2) + t59 * mrSges(5,1) - t129 * t480 / 0.2e1 + qJD(5) * t434 + t49 * t548 + t50 * t550 + t479 * t582 + t447 * t583 + t449 * t584 - t313 * t591 - t314 * t593 + t379 * t597;
t251 = t311 * t421 - t332 * t416;
t213 = pkin(5) * t378 - t380 * t412 + t251;
t252 = t311 * t416 + t320;
t498 = t380 * t416;
t231 = -pkin(12) * t498 + t252;
t134 = t213 * t415 + t231 * t420;
t633 = -qJD(6) * t134 + t415 * t650 + t420 * t651;
t133 = t213 * t420 - t231 * t415;
t632 = qJD(6) * t133 + t415 * t651 - t420 * t650;
t521 = t277 * Ifges(6,3);
t523 = t247 * Ifges(6,5);
t524 = t246 * Ifges(6,6);
t128 = t521 + t523 + t524;
t522 = t273 * Ifges(7,3);
t529 = t162 * Ifges(7,5);
t530 = t458 * Ifges(7,6);
t81 = t522 + t529 + t530;
t631 = t128 + t81;
t627 = -qJD(5) * t252 + t655;
t625 = pkin(5) * t658 + t624;
t620 = -Ifges(7,2) * t162 + t157;
t41 = Ifges(7,6) * t44;
t42 = Ifges(7,5) * t43;
t11 = Ifges(7,3) * t178 + t41 + t42;
t114 = Ifges(6,6) * t117;
t115 = Ifges(6,5) * t116;
t48 = Ifges(6,3) * t178 + t114 + t115;
t614 = t48 + t11;
t443 = t416 * t80 + t421 * t79;
t609 = t443 * mrSges(6,3);
t394 = t589 * t416;
t396 = pkin(11) * t421 + t412;
t331 = t394 * t415 + t396 * t420;
t608 = -qJD(6) * t331 - t415 * t643 + t420 * t642;
t329 = t394 * t420 - t396 * t415;
t607 = qJD(6) * t329 + t415 * t642 + t420 * t643;
t299 = t438 * t380;
t359 = pkin(9) * t414 + t374;
t287 = -t359 * t418 + t360 * t423;
t371 = t414 * t418 + t423 * t496;
t243 = -pkin(3) * t495 - pkin(10) * t371 + t287;
t288 = t359 * t423 + t360 * t418;
t370 = t414 * t423 - t418 * t496;
t257 = pkin(10) * t370 + t288;
t167 = t243 * t417 + t257 * t422;
t154 = -pkin(11) * t495 + t167;
t296 = t370 * t417 + t371 * t422;
t358 = t403 + (-pkin(2) - t546) * t414;
t304 = -pkin(3) * t370 + t358;
t439 = t370 * t422 - t371 * t417;
t194 = -pkin(4) * t439 - pkin(11) * t296 + t304;
t100 = t154 * t421 + t194 * t416;
t276 = Ifges(5,4) * t606;
t512 = t387 * Ifges(5,5);
t518 = t440 * Ifges(5,1);
t197 = t276 + t512 + t518;
t549 = -t421 / 0.2e1;
t604 = t136 * mrSges(5,3) + t609 - t197 / 0.2e1 - t434 - t433 / 0.2e1 - t432 / 0.2e1 - t431 / 0.2e1 - t282 * mrSges(5,2) - t512 / 0.2e1 + t129 * t550 + t130 * t549 - t276 / 0.2e1;
t603 = -t416 * t79 + t421 * t80;
t515 = t348 * Ifges(4,4);
t265 = Ifges(4,2) * t347 + Ifges(4,6) * t393 + t515;
t343 = Ifges(4,4) * t347;
t266 = Ifges(4,1) * t348 + Ifges(4,5) * t393 + t343;
t441 = t268 * t423 + t269 * t418;
t538 = Ifges(4,4) * t423;
t539 = Ifges(4,4) * t418;
t547 = t423 / 0.2e1;
t552 = t393 / 0.2e1;
t556 = t348 / 0.2e1;
t558 = t347 / 0.2e1;
t600 = -t441 * mrSges(4,3) + t325 * (mrSges(4,1) * t418 + mrSges(4,2) * t423) + (-Ifges(4,2) * t418 + t538) * t558 + (Ifges(4,1) * t423 - t539) * t556 + (Ifges(4,5) * t423 - Ifges(4,6) * t418) * t552 - t418 * t265 / 0.2e1 + t266 * t547;
t599 = Ifges(5,2) / 0.2e1;
t587 = pkin(1) * mrSges(3,1);
t586 = pkin(1) * mrSges(3,2);
t575 = -t246 / 0.2e1;
t574 = t246 / 0.2e1;
t573 = -t247 / 0.2e1;
t572 = t247 / 0.2e1;
t568 = -t277 / 0.2e1;
t567 = t277 / 0.2e1;
t566 = t606 / 0.2e1;
t565 = t440 / 0.2e1;
t564 = t439 / 0.2e1;
t562 = t296 / 0.2e1;
t561 = t305 / 0.2e1;
t560 = t306 / 0.2e1;
t557 = -t348 / 0.2e1;
t555 = t370 / 0.2e1;
t554 = t371 / 0.2e1;
t553 = t387 / 0.2e1;
t551 = -t416 / 0.2e1;
t545 = pkin(3) * t422;
t540 = Ifges(3,4) * t419;
t532 = Ifges(3,5) * t424;
t527 = t177 * Ifges(5,1);
t525 = t178 * Ifges(5,4);
t519 = t606 * Ifges(5,6);
t517 = t440 * Ifges(5,5);
t516 = t347 * Ifges(4,6);
t514 = t348 * Ifges(4,5);
t513 = t356 * mrSges(3,2);
t510 = t387 * Ifges(5,3);
t509 = t393 * Ifges(4,3);
t508 = t402 * Ifges(3,5);
t151 = -t314 * t380 + t315 * t438;
t222 = t292 * t415 + t293 * t420;
t494 = t151 - t222;
t152 = t299 * t601 + t379 * t315;
t221 = t292 * t420 - t293 * t415;
t493 = t152 - t221;
t168 = -mrSges(6,1) * t246 + mrSges(6,2) * t247;
t256 = mrSges(5,1) * t387 - mrSges(5,3) * t440;
t492 = t168 - t256;
t491 = -mrSges(3,1) * t402 - mrSges(4,1) * t347 + mrSges(4,2) * t348 + mrSges(3,3) * t469;
t472 = Ifges(4,5) * t305 + Ifges(4,6) * t306 + Ifges(4,3) * t454;
t410 = -pkin(5) * t421 - pkin(4);
t465 = t419 * t486;
t99 = -t154 * t416 + t194 * t421;
t166 = t243 * t422 - t257 * t417;
t456 = -t293 - t499;
t153 = pkin(4) * t495 - t166;
t436 = -t296 * t421 + t416 * t495;
t75 = -pkin(5) * t439 + pkin(12) * t436 + t99;
t274 = -t296 * t416 - t421 * t495;
t85 = pkin(12) * t274 + t100;
t37 = -t415 * t85 + t420 * t75;
t38 = t415 * t75 + t420 * t85;
t77 = mrSges(6,1) * t178 - mrSges(6,3) * t116;
t78 = -mrSges(6,2) * t178 + mrSges(6,3) * t117;
t444 = -t416 * t77 + t421 * t78;
t442 = t205 * t423 - t206 * t418;
t202 = t274 * t420 + t415 * t436;
t203 = t274 * t415 - t420 * t436;
t220 = -qJD(3) * t288 + t366 * t423 - t368 * t418;
t318 = qJD(3) * t370 + t413 * t466;
t181 = pkin(3) * t465 - pkin(10) * t318 + t220;
t219 = -t359 * t484 + t360 * t483 + t366 * t418 + t368 * t423;
t317 = -qJD(3) * t371 - t413 * t467;
t190 = pkin(10) * t317 + t219;
t71 = t181 * t422 - t190 * t417 - t243 * t482 - t257 * t481;
t209 = qJD(4) * t439 + t317 * t417 + t318 * t422;
t210 = qJD(4) * t296 - t317 * t422 + t318 * t417;
t285 = -pkin(3) * t317 + t369;
t105 = pkin(4) * t210 - pkin(11) * t209 + t285;
t70 = t181 * t417 + t190 * t422 + t243 * t481 - t257 * t482;
t67 = pkin(11) * t465 + t70;
t20 = t105 * t416 - t154 * t480 + t194 * t479 + t421 * t67;
t68 = -pkin(4) * t465 - t71;
t21 = -qJD(5) * t100 + t105 * t421 - t416 * t67;
t429 = m(6) * (-qJD(5) * t443 + t621);
t426 = t636 + t537 / 0.2e1 - t524 / 0.2e1 - t529 / 0.2e1 - t521 / 0.2e1 - t81 / 0.2e1 - t522 / 0.2e1 - t128 / 0.2e1 - t530 / 0.2e1 - t523 / 0.2e1 + t511 / 0.2e1;
t398 = Ifges(3,4) * t468;
t392 = t410 - t545;
t391 = t462 * t532;
t363 = -mrSges(3,2) * t402 + mrSges(3,3) * t468;
t323 = Ifges(3,1) * t469 + t398 + t508;
t322 = Ifges(3,6) * t402 + (Ifges(3,2) * t424 + t540) * t488;
t308 = mrSges(4,1) * t393 - mrSges(4,3) * t348;
t307 = -mrSges(4,2) * t393 + mrSges(4,3) * t347;
t298 = t379 * t380;
t291 = pkin(5) * t498 - t605;
t284 = -mrSges(4,2) * t454 + mrSges(4,3) * t306;
t283 = mrSges(4,1) * t454 - mrSges(4,3) * t305;
t264 = t509 + t514 + t516;
t255 = -mrSges(5,2) * t387 + mrSges(5,3) * t606;
t236 = -mrSges(4,1) * t306 + mrSges(4,2) * t305;
t226 = Ifges(4,1) * t305 + Ifges(4,4) * t306 + Ifges(4,5) * t454;
t225 = Ifges(4,4) * t305 + Ifges(4,2) * t306 + Ifges(4,6) * t454;
t211 = -mrSges(5,1) * t606 + mrSges(5,2) * t440;
t195 = t510 + t517 + t519;
t186 = mrSges(6,1) * t277 - mrSges(6,3) * t247;
t185 = -mrSges(6,2) * t277 + mrSges(6,3) * t246;
t164 = -mrSges(5,2) * t454 - mrSges(5,3) * t178;
t163 = mrSges(5,1) * t454 - mrSges(5,3) * t177;
t140 = qJD(5) * t436 - t209 * t416 + t421 * t465;
t139 = qJD(5) * t274 + t209 * t421 + t416 * t465;
t125 = mrSges(7,1) * t273 - mrSges(7,3) * t162;
t124 = -mrSges(7,2) * t273 + mrSges(7,3) * t458;
t123 = -pkin(5) * t274 + t153;
t108 = t137 + t648;
t104 = mrSges(5,1) * t178 + mrSges(5,2) * t177;
t98 = Ifges(5,5) * t454 - t525 + t527;
t97 = Ifges(5,6) * t454 + t526 - t646;
t93 = -mrSges(7,1) * t458 + mrSges(7,2) * t162;
t62 = -qJD(6) * t203 - t139 * t415 + t140 * t420;
t61 = qJD(6) * t202 + t139 * t420 + t140 * t415;
t45 = -pkin(5) * t140 + t68;
t35 = -mrSges(7,2) * t178 + mrSges(7,3) * t44;
t34 = mrSges(7,1) * t178 - mrSges(7,3) * t43;
t26 = t420 * t64 - t507;
t25 = -t415 * t64 - t505;
t19 = pkin(12) * t140 + t20;
t18 = -mrSges(7,1) * t44 + mrSges(7,2) * t43;
t14 = pkin(5) * t210 - pkin(12) * t139 + t21;
t6 = -qJD(6) * t38 + t14 * t420 - t19 * t415;
t5 = qJD(6) * t37 + t14 * t415 + t19 * t420;
t1 = [-t178 * (Ifges(5,4) * t296 - Ifges(5,6) * t495) / 0.2e1 + (Ifges(5,4) * t209 + Ifges(5,6) * t465) * t566 + (Ifges(7,5) * t61 + Ifges(7,6) * t62) * t570 + (Ifges(7,4) * t203 + Ifges(7,2) * t202) * t595 + (Ifges(7,4) * t61 + Ifges(7,2) * t62) * t580 + t177 * (Ifges(5,1) * t296 - Ifges(5,5) * t495) / 0.2e1 + (Ifges(5,1) * t209 + Ifges(5,5) * t465) * t565 + (-Ifges(6,4) * t436 + Ifges(6,2) * t274) * t583 - (t472 + t473) * t495 / 0.2e1 + (t202 * t3 - t203 * t4 - t23 * t61 + t24 * t62) * mrSges(7,3) + (Ifges(6,4) * t139 + Ifges(6,2) * t140) * t574 + (Ifges(5,5) * t209 + Ifges(5,3) * t465) * t553 + ((-Ifges(6,3) - Ifges(7,3)) * t577 - Ifges(6,6) * t583 - Ifges(6,5) * t584 - t646 / 0.2e1 - Ifges(7,6) * t595 - Ifges(7,5) * t596 - t614 / 0.2e1 + t654) * t439 + (-Ifges(6,1) * t436 + Ifges(6,4) * t274) * t584 + (Ifges(6,5) * t139 + Ifges(6,6) * t140) * t567 + ((t264 + t195) * t419 + t402 * (-Ifges(3,6) * t419 + t532) + t424 * t323) * t486 / 0.2e1 + t56 * (-mrSges(6,1) * t274 - mrSges(6,2) * t436) - t436 * t50 / 0.2e1 + (-t137 * t465 + t209 * t282 + t267 * t296 + t495 * t58) * mrSges(5,2) + (-t636 + t631 / 0.2e1 - Ifges(5,6) * t553 - Ifges(5,4) * t565 - Ifges(5,2) * t566 + Ifges(6,3) * t567 + Ifges(6,5) * t572 + Ifges(6,6) * t574 + Ifges(7,3) * t570 + Ifges(7,5) * t578 + Ifges(7,6) * t580) * t210 + (-Ifges(6,5) * t436 + Ifges(7,5) * t203 + Ifges(6,6) * t274 + Ifges(7,6) * t202) * t577 + (Ifges(6,1) * t139 + Ifges(6,4) * t140) * t572 + (Ifges(7,1) * t203 + Ifges(7,4) * t202) * t596 + (Ifges(7,1) * t61 + Ifges(7,4) * t62) * t578 + (-t139 * t79 + t140 * t80 + t16 * t274 + t17 * t436) * mrSges(6,3) + t37 * t34 + t38 * t35 + (t356 * t424 + t357 * t419 + (-t364 * t424 - t367 * t419) * qJD(2)) * t413 * mrSges(3,3) + t45 * t93 + t99 * t77 + t100 * t78 + t106 * (-mrSges(7,1) * t62 + mrSges(7,2) * t61) - t322 * t465 / 0.2e1 + t269 * (-mrSges(4,2) * t465 + mrSges(4,3) * t317) + t136 * (mrSges(5,1) * t465 - mrSges(5,3) * t209) + t268 * (mrSges(4,1) * t465 - mrSges(4,3) * t318) + t123 * t18 + t5 * t124 + t6 * t125 + t140 * t129 / 0.2e1 + t131 * (-mrSges(6,1) * t140 + mrSges(6,2) * t139) + t153 * t63 + t166 * t163 + t167 * t164 + t68 * t168 + m(7) * (t106 * t45 + t123 * t36 + t23 * t6 + t24 * t5 + t3 * t38 + t37 * t4) + m(6) * (t100 * t16 + t131 * t68 + t153 * t56 + t17 * t99 + t20 * t80 + t21 * t79) + m(5) * (t136 * t71 + t137 * t70 + t166 * t59 + t167 * t58 + t267 * t304 + t282 * t285) + m(4) * (t205 * t288 + t206 * t287 + t219 * t269 + t220 * t268 + t325 * t369 + t357 * t358) + m(3) * (t356 * t374 - t357 * t373 - t364 * t369 + t367 * t368) + t20 * t185 + t21 * t186 + t491 * t369 + t36 * (-mrSges(7,1) * t202 + mrSges(7,2) * t203) + t209 * t197 / 0.2e1 + t59 * (-mrSges(5,1) * t495 - mrSges(5,3) * t296) + t206 * (-mrSges(4,1) * t495 - mrSges(4,3) * t371) + t205 * (mrSges(4,2) * t495 + mrSges(4,3) * t370) + t70 * t255 + t71 * t256 + (t391 / 0.2e1 - t513 - t357 * mrSges(3,1)) * t414 + t274 * t49 / 0.2e1 + t285 * t211 + t287 * t283 + t288 * t284 + t304 * t104 + t219 * t307 + t220 * t308 + t317 * t265 / 0.2e1 + t318 * t266 / 0.2e1 + t325 * (-mrSges(4,1) * t317 + mrSges(4,2) * t318) + ((Ifges(3,5) * t414 / 0.2e1 - t373 * mrSges(3,3) + (-0.2e1 * t586 + 0.3e1 / 0.2e1 * Ifges(3,4) * t424) * t413) * t424 + (Ifges(5,5) * t562 + Ifges(5,6) * t564 + Ifges(4,5) * t554 + Ifges(4,6) * t555 - Ifges(3,6) * t414 - t374 * mrSges(3,3) + (-0.2e1 * t587 - 0.3e1 / 0.2e1 * t540) * t413 + (-Ifges(5,3) / 0.2e1 - Ifges(4,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t495) * t419) * t462 + t358 * t236 + t368 * t363 + t357 * (-mrSges(4,1) * t370 + mrSges(4,2) * t371) + (Ifges(4,5) * t318 + Ifges(4,6) * t317 + Ifges(4,3) * t465) * t552 + t226 * t554 + t225 * t555 + (Ifges(4,1) * t318 + Ifges(4,4) * t317 + Ifges(4,5) * t465) * t556 + (Ifges(4,4) * t318 + Ifges(4,2) * t317 + Ifges(4,6) * t465) * t558 + (Ifges(4,4) * t371 + Ifges(4,2) * t370 - Ifges(4,6) * t495) * t560 + (Ifges(4,1) * t371 + Ifges(4,4) * t370 - Ifges(4,5) * t495) * t561 + t98 * t562 + t97 * t564 + t139 * t582 + t61 * t591 + t62 * t593 + t203 * t597 + t202 * t598; (-t268 * t289 - t269 * t290 - t325 * t367 + (-qJD(3) * t441 + t442) * pkin(9) - pkin(2) * t357) * m(4) + (t418 * pkin(3) * t211 + (-t307 * t418 - t308 * t423) * pkin(9) + t600) * qJD(3) + ((-t398 / 0.2e1 - t323 / 0.2e1 - t508 / 0.2e1 + t364 * mrSges(3,3) + t488 * t586 - t600) * t424 + (-t268 * mrSges(4,1) - t516 / 0.2e1 - t514 / 0.2e1 - t509 / 0.2e1 + t269 * mrSges(4,2) + t322 / 0.2e1 - t264 / 0.2e1 - t195 / 0.2e1 + t137 * mrSges(5,2) - t519 / 0.2e1 - t517 / 0.2e1 - t510 / 0.2e1 + t367 * mrSges(3,3) - t136 * mrSges(5,1) + (t587 + t540 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t424) * t488 + (t402 / 0.2e1 - qJD(2)) * Ifges(3,6) + (Ifges(4,5) * t418 + Ifges(5,5) * t380 + Ifges(4,6) * t423 - Ifges(5,6) * t378) * qJD(2) / 0.2e1) * t419) * t488 + (-t23 * t494 + t24 * t493 - t298 * t3 + t299 * t4) * mrSges(7,3) + (-Ifges(7,5) * t299 - Ifges(7,6) * t298) * t577 + (-Ifges(7,4) * t299 - Ifges(7,2) * t298) * t595 + (-Ifges(7,1) * t299 - Ifges(7,4) * t298) * t596 + (t106 * t494 - t24 * t489 - t299 * t36) * mrSges(7,2) - t440 * (-Ifges(5,1) * t328 - Ifges(5,4) * t327) / 0.2e1 - t387 * (-Ifges(5,5) * t328 - Ifges(5,6) * t327) / 0.2e1 + (mrSges(5,1) * t489 - mrSges(5,2) * t662) * t282 + (t136 * t662 - t137 * t489) * mrSges(5,3) + (t50 * t548 + t49 * t551 + t56 * t451 + t450 * t584 + t448 * t583 + t446 * t577 + t98 / 0.2e1 + t527 / 0.2e1 - t525 / 0.2e1 + t267 * mrSges(5,2) - t59 * mrSges(5,3) + (-t16 * t416 - t17 * t421) * mrSges(6,3) + (-mrSges(6,3) * t603 + t129 * t549 + t130 * t551 + t131 * t452 + t445 * t568 + t447 * t575 + t449 * t573) * qJD(5)) * t380 + (-t315 / 0.2e1 + t328 / 0.2e1) * t197 + (-Ifges(5,5) * t315 - Ifges(5,6) * t316) * t553 + (-Ifges(5,1) * t315 - Ifges(5,4) * t316) * t565 + (-Ifges(5,4) * t315 - Ifges(5,2) * t316) * t566 + (-t499 / 0.2e1 - t293 / 0.2e1) * t130 + (-Ifges(6,5) * t499 + Ifges(6,6) * t500 + Ifges(6,3) * t316) * t567 + (-Ifges(6,1) * t499 + Ifges(6,4) * t500 + Ifges(6,5) * t316) * t572 + (-Ifges(6,4) * t499 + Ifges(6,2) * t500 + Ifges(6,6) * t316) * t574 + (t500 / 0.2e1 - t292 / 0.2e1) * t129 + (mrSges(6,1) * t457 + mrSges(6,2) * t456) * t131 + (t48 / 0.2e1 + t11 / 0.2e1 + t41 / 0.2e1 + t42 / 0.2e1 + t115 / 0.2e1 + t114 / 0.2e1 - t97 / 0.2e1 + (t599 + Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t178 - t654) * t378 + (t151 / 0.2e1 - t222 / 0.2e1) * t83 + (-t283 * t418 + t284 * t423) * pkin(9) + (t131 * t624 + t16 * t252 + t17 * t251 - t56 * t605 + t626 * t80 + t627 * t79) * m(6) - (t63 - t163) * t605 + t391 + (-t106 * t493 + t23 * t489 + t298 * t36) * mrSges(7,1) + t622 * t255 + t623 * t256 + t624 * t168 + t625 * t93 + t626 * t185 + t627 * t186 - t513 + (-mrSges(4,1) * t423 + mrSges(4,2) * t418 - mrSges(3,1)) * t357 + (-t196 + t631) * (t316 / 0.2e1 - t327 / 0.2e1) + (t152 / 0.2e1 - t221 / 0.2e1) * t82 - t606 * (-Ifges(5,4) * t328 - Ifges(5,2) * t327) / 0.2e1 + t442 * mrSges(4,3) + t632 * t124 + t633 * t125 + (t106 * t625 + t133 * t4 + t134 * t3 + t23 * t633 + t24 * t632 + t291 * t36) * m(7) + (t623 * t136 + t622 * t137 + t267 * t411 + t282 * t661 + t332 * t58 + t605 * t59) * m(5) + t133 * t34 + t134 * t35 + (mrSges(6,1) * t489 - mrSges(6,3) * t456) * t79 + (-mrSges(6,2) * t489 - mrSges(6,3) * t457) * t80 - t491 * t367 - pkin(2) * t236 + t251 * t77 + t252 * t78 + t291 * t18 - t290 * t307 - t289 * t308 - t319 * t211 + t332 * t164 - t364 * t363 + t225 * t547 + (Ifges(4,2) * t423 + t539) * t560 + (Ifges(4,1) * t418 + t538) * t561 + (Ifges(6,5) * t293 + Ifges(6,6) * t292 + Ifges(6,3) * t327) * t568 + (Ifges(7,5) * t151 + Ifges(7,6) * t152 + Ifges(7,3) * t316) * t570 + (Ifges(7,5) * t222 + Ifges(7,6) * t221 + Ifges(7,3) * t327) * t571 + (Ifges(6,1) * t293 + Ifges(6,4) * t292 + Ifges(6,5) * t327) * t573 + (Ifges(6,4) * t293 + Ifges(6,2) * t292 + Ifges(6,6) * t327) * t575 + (Ifges(7,1) * t151 + Ifges(7,4) * t152 + Ifges(7,5) * t316) * t578 + (Ifges(7,1) * t222 + Ifges(7,4) * t221 + Ifges(7,5) * t327) * t579 + (Ifges(7,4) * t151 + Ifges(7,2) * t152 + Ifges(7,6) * t316) * t580 + (Ifges(7,4) * t222 + Ifges(7,2) * t221 + Ifges(7,6) * t327) * t581 - t299 * t597 - t298 * t598 + t411 * t104 + t418 * t226 / 0.2e1; -t609 * qJD(5) + t635 + (t426 + t520 / 0.2e1) * t440 + t640 * t93 - (-Ifges(4,2) * t348 + t266 + t343) * t347 / 0.2e1 - m(6) * (t131 * t142 + t79 * t90 + t80 * t91) + t472 + (-t518 / 0.2e1 + t604) * t606 + (t268 * t347 + t269 * t348) * mrSges(4,3) + ((-t185 * t416 - t186 * t421) * qJD(5) + t444 + t429) * t408 - t641 * (-pkin(4) - t545) + t644 * t124 + t645 * t125 + (t106 * t640 + t23 * t645 + t24 * t644 + t3 * t310 + t309 * t4 + t36 * t392) * m(7) - m(5) * (-t136 * t142 + t137 * t143) - t91 * t185 - t90 * t186 - t492 * t142 - t205 * mrSges(4,2) + t206 * mrSges(4,1) - t143 * t255 - t268 * t307 + t269 * t308 + t309 * t34 + t310 * t35 + (t422 * t163 + t417 * t164 - t348 * t211 + ((m(6) * t131 + t492) * t417 + (m(6) * t603 + t185 * t421 - t186 * t416 + t255) * t422) * qJD(4) + (t417 * t58 + t422 * t59 + 0.2e1 * t282 * t557 + (-t136 * t417 + t137 * t422) * qJD(4)) * m(5)) * pkin(3) - t325 * (mrSges(4,1) * t348 + mrSges(4,2) * t347) + t265 * t556 + (Ifges(4,1) * t347 - t515) * t557 + t392 * t18 - t393 * (Ifges(4,5) * t347 - Ifges(4,6) * t348) / 0.2e1; t635 + t607 * t124 + (t3 * t331 + t329 * t4 + t36 * t410 + t607 * t24 + t608 * t23 + (-t108 + t474) * t106) * m(7) + t608 * t125 + t426 * t440 + ((-mrSges(6,3) * t79 - pkin(11) * t186) * t421 + (-t80 * mrSges(6,3) + pkin(5) * t93 - pkin(11) * t185) * t416) * qJD(5) + t641 * pkin(4) - m(6) * (t131 * t137 + t79 * t94 + t80 * t95) + t444 * pkin(11) + ((t599 - Ifges(5,1) / 0.2e1) * t440 + t604) * t606 - t108 * t93 - t95 * t185 - t94 * t186 - t492 * t137 - t136 * t255 + t329 * t34 + t331 * t35 + pkin(11) * t429 + t410 * t18; -t162 * t594 + t458 * t592 + t614 + t615 + (-t247 * t93 + t420 * t34 + t415 * t35 + (t124 * t420 - t125 * t415) * qJD(6) + (-t106 * t247 + t3 * t415 + t4 * t420 + (-t23 * t415 + t24 * t420) * qJD(6)) * m(7)) * pkin(5) + (-Ifges(6,2) * t247 + t130 + t244) * t575 + (t246 * t79 + t247 * t80) * mrSges(6,3) - m(7) * (t23 * t25 + t24 * t26) - t26 * t124 - t25 * t125 - t79 * t185 + t80 * t186 - t131 * (mrSges(6,1) * t247 + mrSges(6,2) * t246) + (Ifges(6,5) * t246 - Ifges(6,6) * t247) * t568 + t129 * t572 + (Ifges(6,1) * t246 - t536) * t573 + t620 * t581 + t637; -t23 * t124 + t24 * t125 + t82 * t578 + t11 + (t620 + t83) * t581 + t634 + t637;];
tauc  = t1(:);
