% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:29:01
% EndTime: 2019-03-09 10:30:16
% DurationCPUTime: 47.53s
% Computational Cost: add. (24496->1014), mult. (67944->1377), div. (0->0), fcn. (56247->16), ass. (0->446)
t375 = sin(qJ(2));
t379 = cos(qJ(2));
t510 = sin(pkin(11));
t511 = cos(pkin(11));
t327 = t375 * t510 - t379 * t511;
t369 = sin(pkin(6));
t609 = t369 * t327;
t287 = qJD(1) * t609;
t591 = t287 + qJD(4);
t371 = cos(pkin(6));
t539 = pkin(1) * t371;
t358 = t379 * t539;
t352 = qJD(1) * t358;
t535 = pkin(8) + qJ(3);
t447 = t535 * t375;
t426 = t369 * t447;
t278 = -qJD(1) * t426 + t352;
t493 = t371 * t375;
t357 = pkin(1) * t493;
t497 = t369 * t379;
t279 = (t497 * t535 + t357) * qJD(1);
t438 = t511 * t279;
t194 = t278 * t510 + t438;
t374 = sin(qJ(4));
t378 = cos(qJ(4));
t652 = -qJD(5) * t374 - t194 + t591 * (pkin(4) * t374 - qJ(5) * t378);
t266 = t510 * t279;
t195 = t278 * t511 - t266;
t389 = t375 * t511 + t379 * t510;
t482 = qJD(1) * t369;
t288 = t389 * t482;
t457 = t375 * t482;
t433 = pkin(2) * t457;
t209 = pkin(3) * t288 + pkin(9) * t287 + t433;
t127 = t378 * t195 + t374 * t209;
t111 = qJ(5) * t288 + t127;
t368 = sin(pkin(12));
t370 = cos(pkin(12));
t458 = t510 * pkin(2);
t361 = t458 + pkin(9);
t477 = qJD(4) * t374;
t453 = t361 * t477;
t612 = t652 * t370 + (t111 + t453) * t368;
t651 = -t370 * t111 + t368 * t652;
t362 = pkin(5) * t370 + pkin(4);
t418 = -mrSges(6,1) * t370 + mrSges(6,2) * t368;
t650 = -m(6) * pkin(4) - m(7) * t362 + t418;
t393 = t371 * pkin(2) - t426;
t256 = qJD(2) * pkin(2) + qJD(1) * t393 + t352;
t175 = t256 * t511 - t266;
t431 = qJD(1) * t371 + qJD(2);
t166 = -pkin(3) * t431 - t175;
t254 = t374 * t288 - t378 * t431;
t255 = t378 * t288 + t374 * t431;
t101 = t254 * pkin(4) - t255 * qJ(5) + t166;
t176 = t510 * t256 + t438;
t167 = pkin(9) * t431 + t176;
t364 = pkin(2) * t379 + pkin(1);
t314 = -t364 * t482 + qJD(3);
t191 = pkin(3) * t287 - pkin(9) * t288 + t314;
t108 = t167 * t378 + t191 * t374;
t283 = t327 * t482 + qJD(4);
t94 = qJ(5) * t283 + t108;
t49 = t370 * t101 - t368 * t94;
t50 = t368 * t101 + t370 * t94;
t522 = t108 * mrSges(5,3);
t649 = t49 * mrSges(6,1) - t50 * mrSges(6,2) - t522;
t419 = mrSges(5,1) * t378 - mrSges(5,2) * t374;
t534 = pkin(10) + qJ(5);
t599 = -m(6) * qJ(5) - m(7) * t534 - mrSges(6,3) - mrSges(7,3);
t367 = pkin(12) + qJ(6);
t365 = sin(t367);
t366 = cos(t367);
t592 = -mrSges(7,1) * t366 + mrSges(7,2) * t365;
t600 = -t592 - t650;
t586 = t374 * t599 - t378 * t600 - mrSges(4,1) - t419;
t520 = t255 * Ifges(5,4);
t642 = t283 * Ifges(5,6);
t145 = -t254 * Ifges(5,2) + t520 + t642;
t182 = t255 * t370 + t283 * t368;
t373 = sin(qJ(6));
t377 = cos(qJ(6));
t436 = -t255 * t368 + t370 * t283;
t122 = t182 * t377 + t373 * t436;
t253 = qJD(6) + t254;
t625 = -t182 * t373 + t377 * t436;
t644 = t182 * Ifges(6,5);
t647 = Ifges(6,6) * t436;
t619 = t122 * Ifges(7,5) + Ifges(7,6) * t625 + t254 * Ifges(6,3) + t253 * Ifges(7,3) + t644 + t647;
t648 = t145 / 0.2e1 - t619 / 0.2e1;
t107 = -t374 * t167 + t191 * t378;
t646 = t107 * mrSges(5,3);
t645 = t166 * mrSges(5,2);
t643 = t283 * Ifges(5,5);
t494 = t370 * t378;
t214 = -t287 * t494 + t288 * t368;
t504 = t287 * t374;
t641 = pkin(5) * t504 + pkin(10) * t214 + (pkin(5) * t374 - pkin(10) * t494) * qJD(4) + t612;
t500 = t368 * t378;
t213 = t287 * t500 + t288 * t370;
t495 = t370 * t374;
t640 = -pkin(10) * t213 + (-pkin(10) * t500 - t361 * t495) * qJD(4) + t651;
t551 = -t254 / 0.2e1;
t553 = -t253 / 0.2e1;
t556 = -t182 / 0.2e1;
t557 = -t436 / 0.2e1;
t564 = -t122 / 0.2e1;
t566 = -t625 / 0.2e1;
t637 = Ifges(6,5) * t556 + Ifges(7,5) * t564 + Ifges(6,6) * t557 + Ifges(7,6) * t566 + Ifges(6,3) * t551 + Ifges(7,3) * t553 - t649;
t36 = pkin(5) * t254 - pkin(10) * t182 + t49;
t39 = pkin(10) * t436 + t50;
t13 = t36 * t377 - t373 * t39;
t474 = qJD(1) * qJD(2);
t309 = (qJDD(1) * t379 - t375 * t474) * t369;
t310 = (qJDD(1) * t375 + t379 * t474) * t369;
t232 = t309 * t511 - t310 * t510;
t230 = qJDD(4) - t232;
t233 = t309 * t510 + t310 * t511;
t473 = qJDD(1) * t369;
t271 = -pkin(1) * t473 - pkin(2) * t309 + qJDD(3);
t136 = -pkin(3) * t232 - pkin(9) * t233 + t271;
t476 = qJD(4) * t378;
t483 = pkin(8) * t497 + t357;
t308 = t483 * qJD(2);
t472 = qJDD(1) * t371;
t466 = pkin(1) * t472;
t350 = t379 * t466;
t354 = qJDD(2) + t472;
t499 = t369 * t375;
t454 = qJD(3) * t499;
t465 = pkin(8) * t473;
t172 = -t375 * t465 + pkin(2) * t354 - qJ(3) * t310 + t350 + (-t308 - t454) * qJD(1);
t471 = qJD(2) * t539;
t430 = qJD(1) * t471;
t461 = t375 * t466 + (t430 + t465) * t379;
t479 = qJD(3) * t379;
t480 = qJD(2) * t375;
t184 = qJ(3) * t309 + (-pkin(8) * t480 + t479) * t482 + t461;
t106 = t510 * t172 + t511 * t184;
t98 = pkin(9) * t354 + t106;
t37 = t374 * t136 - t167 * t477 + t191 * t476 + t378 * t98;
t29 = qJ(5) * t230 + qJD(5) * t283 + t37;
t478 = qJD(4) * t254;
t151 = t378 * t233 + t374 * t354 - t478;
t152 = qJD(4) * t255 + t374 * t233 - t378 * t354;
t105 = t172 * t511 - t510 * t184;
t97 = -t354 * pkin(3) - t105;
t45 = t152 * pkin(4) - t151 * qJ(5) - t255 * qJD(5) + t97;
t10 = -t29 * t368 + t370 * t45;
t110 = t151 * t370 + t230 * t368;
t5 = pkin(5) * t152 - pkin(10) * t110 + t10;
t109 = -t151 * t368 + t230 * t370;
t11 = t370 * t29 + t368 * t45;
t6 = pkin(10) * t109 + t11;
t1 = qJD(6) * t13 + t373 * t5 + t377 * t6;
t14 = t36 * t373 + t377 * t39;
t2 = -qJD(6) * t14 - t373 * t6 + t377 * t5;
t636 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t634 = t166 * mrSges(5,1) + t13 * mrSges(7,1) - t14 * mrSges(7,2);
t554 = t230 / 0.2e1;
t559 = -t152 / 0.2e1;
t560 = t151 / 0.2e1;
t571 = Ifges(5,1) * t560 + Ifges(5,4) * t559 + Ifges(5,5) * t554;
t552 = t253 / 0.2e1;
t563 = t122 / 0.2e1;
t565 = t625 / 0.2e1;
t633 = Ifges(7,5) * t563 + Ifges(7,6) * t565 + Ifges(7,3) * t552 - t648;
t150 = qJDD(6) + t152;
t34 = qJD(6) * t625 + t109 * t373 + t110 * t377;
t35 = -qJD(6) * t122 + t109 * t377 - t110 * t373;
t7 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t150;
t631 = t110 * Ifges(6,5) + t109 * Ifges(6,6) + t152 * Ifges(6,3) + t7;
t630 = Ifges(4,4) * t287;
t629 = Ifges(4,2) * t287;
t628 = t283 * Ifges(5,3);
t416 = t365 * mrSges(7,1) + t366 * mrSges(7,2);
t514 = t370 * mrSges(6,2);
t417 = t368 * mrSges(6,1) + t514;
t538 = pkin(5) * t368;
t622 = mrSges(5,3) - mrSges(4,2);
t627 = -m(7) * (pkin(9) + t538) - t416 - m(6) * pkin(9) - t417 - t622;
t380 = cos(qJ(1));
t488 = t379 * t380;
t376 = sin(qJ(1));
t490 = t376 * t375;
t626 = t371 * t488 - t490;
t579 = t34 / 0.2e1;
t578 = t35 / 0.2e1;
t623 = -m(7) - m(6);
t568 = t109 / 0.2e1;
t567 = t110 / 0.2e1;
t561 = t150 / 0.2e1;
t558 = t152 / 0.2e1;
t459 = t511 * pkin(2);
t363 = -t459 - pkin(3);
t324 = -t378 * pkin(4) - t374 * qJ(5) + t363;
t302 = t370 * t324;
t234 = -pkin(10) * t495 + t302 + (-t361 * t368 - pkin(5)) * t378;
t264 = t368 * t324 + t361 * t494;
t501 = t368 * t374;
t248 = -pkin(10) * t501 + t264;
t154 = t234 * t373 + t248 * t377;
t621 = -qJD(6) * t154 - t373 * t640 + t377 * t641;
t153 = t234 * t377 - t248 * t373;
t620 = qJD(6) * t153 + t373 * t641 + t377 * t640;
t618 = t107 * mrSges(5,1);
t617 = t108 * mrSges(5,2);
t114 = mrSges(5,1) * t230 - mrSges(5,3) * t151;
t61 = -t109 * mrSges(6,1) + t110 * mrSges(6,2);
t615 = -t114 + t61;
t336 = t534 * t368;
t337 = t534 * t370;
t272 = -t336 * t377 - t337 * t373;
t403 = t368 * t373 - t370 * t377;
t505 = t254 * t370;
t168 = pkin(4) * t255 + qJ(5) * t254;
t76 = -t107 * t368 + t370 * t168;
t57 = pkin(5) * t255 + pkin(10) * t505 + t76;
t506 = t254 * t368;
t77 = t370 * t107 + t368 * t168;
t66 = pkin(10) * t506 + t77;
t614 = -qJD(5) * t403 + qJD(6) * t272 - t373 * t57 - t377 * t66;
t273 = -t336 * t373 + t337 * t377;
t328 = t368 * t377 + t370 * t373;
t613 = -qJD(5) * t328 - qJD(6) * t273 + t373 * t66 - t377 * t57;
t611 = -t370 * t453 + t651;
t126 = -t374 * t195 + t209 * t378;
t112 = -pkin(4) * t288 - t126;
t443 = t361 + t538;
t610 = pkin(5) * t213 + t443 * t476 - t112;
t157 = t328 * t254;
t316 = t328 * qJD(6);
t606 = t157 + t316;
t158 = t403 * t254;
t315 = t403 * qJD(6);
t605 = t158 + t315;
t123 = -mrSges(6,1) * t436 + mrSges(6,2) * t182;
t190 = mrSges(5,1) * t283 - mrSges(5,3) * t255;
t604 = t190 - t123;
t138 = t213 * t373 + t214 * t377;
t235 = -t316 * t374 - t403 * t476;
t603 = t235 - t138;
t137 = t213 * t377 - t214 * t373;
t236 = t315 * t374 - t328 * t476;
t602 = t236 - t137;
t518 = t288 * mrSges(4,3);
t601 = mrSges(4,1) * t431 - mrSges(5,1) * t254 - mrSges(5,2) * t255 - t518;
t91 = t182 * Ifges(6,4) + Ifges(6,2) * t436 + Ifges(6,6) * t254;
t570 = -t91 / 0.2e1;
t598 = t368 * t570 - t646;
t38 = t136 * t378 - t167 * t476 - t191 * t477 - t374 * t98;
t30 = -pkin(4) * t230 + qJDD(5) - t38;
t93 = -pkin(4) * t283 + qJD(5) - t107;
t597 = t30 * t374 + t93 * t476;
t596 = t477 + t504;
t140 = -mrSges(6,2) * t254 + mrSges(6,3) * t436;
t141 = mrSges(6,1) * t254 - mrSges(6,3) * t182;
t595 = t140 * t370 - t141 * t368;
t71 = -mrSges(6,2) * t152 + mrSges(6,3) * t109;
t72 = mrSges(6,1) * t152 - mrSges(6,3) * t110;
t594 = -t368 * t72 + t370 * t71;
t593 = t37 * t378 - t374 * t38;
t406 = -t10 * t368 + t11 * t370;
t475 = m(5) - t623;
t531 = Ifges(3,4) * t375;
t590 = -t375 * (Ifges(3,1) * t379 - t531) / 0.2e1 + pkin(1) * (mrSges(3,1) * t375 + mrSges(3,2) * t379);
t588 = mrSges(5,1) + t600;
t391 = mrSges(5,2) + t599;
t484 = t389 * t371;
t246 = -t380 * t327 - t376 * t484;
t587 = t376 * t327 - t380 * t484;
t585 = t38 * mrSges(5,1) - t37 * mrSges(5,2) + Ifges(5,5) * t151 - Ifges(5,6) * t152 + Ifges(5,3) * t230;
t432 = pkin(8) * t457;
t250 = -qJD(2) * t432 + t461;
t251 = -pkin(8) * t310 - t375 * t430 + t350;
t584 = t251 * mrSges(3,1) + t105 * mrSges(4,1) - t250 * mrSges(3,2) - t106 * mrSges(4,2) + Ifges(3,5) * t310 + Ifges(4,5) * t233 + Ifges(3,6) * t309 + Ifges(4,6) * t232;
t583 = m(5) * pkin(9) - t627;
t582 = pkin(9) * t475 + t368 * (m(7) * pkin(5) + mrSges(6,1)) + t514 + t622;
t581 = Ifges(7,4) * t579 + Ifges(7,2) * t578 + Ifges(7,6) * t561;
t580 = Ifges(7,1) * t579 + Ifges(7,4) * t578 + Ifges(7,5) * t561;
t41 = Ifges(6,4) * t110 + Ifges(6,2) * t109 + Ifges(6,6) * t152;
t577 = t41 / 0.2e1;
t576 = Ifges(6,1) * t567 + Ifges(6,4) * t568 + Ifges(6,5) * t558;
t525 = Ifges(7,4) * t122;
t53 = Ifges(7,2) * t625 + t253 * Ifges(7,6) + t525;
t575 = -t53 / 0.2e1;
t574 = t53 / 0.2e1;
t116 = Ifges(7,4) * t625;
t54 = t122 * Ifges(7,1) + t253 * Ifges(7,5) + t116;
t573 = -t54 / 0.2e1;
t572 = t54 / 0.2e1;
t92 = t182 * Ifges(6,1) + Ifges(6,4) * t436 + Ifges(6,5) * t254;
t569 = t92 / 0.2e1;
t550 = t254 / 0.2e1;
t549 = -t255 / 0.2e1;
t548 = t255 / 0.2e1;
t546 = -t283 / 0.2e1;
t544 = t288 / 0.2e1;
t540 = pkin(1) * t369;
t537 = qJD(2) / 0.2e1;
t298 = t389 * t369;
t289 = qJD(2) * t298;
t353 = t379 * t471;
t259 = t353 + (-qJD(2) * t447 + t479) * t369;
t448 = t535 * t369;
t260 = -t454 + (-t379 * t448 - t357) * qJD(2);
t171 = t259 * t511 + t260 * t510;
t274 = t358 + t393;
t284 = qJ(3) * t497 + t483;
t203 = t510 * t274 + t511 * t284;
t193 = pkin(9) * t371 + t203;
t290 = qJD(2) * t609;
t455 = t369 * t480;
t210 = pkin(2) * t455 + pkin(3) * t289 + pkin(9) * t290;
t355 = pkin(2) * t497;
t485 = -pkin(3) * t609 + t355;
t423 = pkin(9) * t298 + t485;
t215 = -t423 - t540;
t69 = t378 * t171 - t193 * t477 + t374 * t210 + t215 * t476;
t58 = qJ(5) * t289 + qJD(5) * t609 + t69;
t170 = t259 * t510 - t511 * t260;
t270 = t298 * t378 + t371 * t374;
t200 = qJD(4) * t270 - t290 * t374;
t269 = t298 * t374 - t371 * t378;
t201 = -qJD(4) * t269 - t290 * t378;
t80 = t200 * pkin(4) - t201 * qJ(5) - t270 * qJD(5) + t170;
t27 = t368 * t80 + t370 * t58;
t530 = Ifges(4,4) * t288;
t529 = Ifges(5,4) * t374;
t528 = Ifges(5,4) * t378;
t527 = Ifges(6,4) * t368;
t526 = Ifges(6,4) * t370;
t519 = t287 * mrSges(4,3);
t509 = Ifges(3,6) * qJD(2);
t503 = t287 * t378;
t498 = t369 * t376;
t496 = t369 * t380;
t492 = t375 * t380;
t489 = t376 * t379;
t130 = t378 * t193 + t374 * t215;
t117 = qJ(5) * t609 + t130;
t202 = t274 * t511 - t510 * t284;
t192 = -t371 * pkin(3) - t202;
t128 = t269 * pkin(4) - t270 * qJ(5) + t192;
t65 = t370 * t117 + t368 * t128;
t469 = m(4) + t475;
t468 = -m(3) * pkin(1) - mrSges(2,1);
t456 = t379 * t482;
t452 = t361 * t476;
t12 = -t35 * mrSges(7,1) + t34 * mrSges(7,2);
t444 = t476 / 0.2e1;
t26 = -t368 * t58 + t370 * t80;
t437 = -t232 * mrSges(4,1) + t233 * mrSges(4,2);
t64 = -t117 * t368 + t370 * t128;
t129 = -t374 * t193 + t215 * t378;
t217 = -t374 * t496 - t378 * t587;
t311 = pkin(2) * t493 - t448;
t435 = -t311 * t376 + t380 * t364;
t429 = mrSges(3,3) * t457;
t428 = mrSges(3,3) * t456;
t424 = t626 * pkin(2);
t420 = mrSges(5,1) * t269 + mrSges(5,2) * t270;
t415 = Ifges(5,1) * t378 - t529;
t414 = Ifges(6,1) * t370 - t527;
t413 = -Ifges(5,2) * t374 + t528;
t412 = -Ifges(6,2) * t368 + t526;
t411 = Ifges(3,5) * t379 - Ifges(3,6) * t375;
t410 = -Ifges(4,5) * t290 - Ifges(4,6) * t289;
t409 = Ifges(5,5) * t378 - Ifges(5,6) * t374;
t408 = Ifges(6,5) * t370 - Ifges(6,6) * t368;
t405 = -t368 * t49 + t370 * t50;
t208 = t270 * t370 + t368 * t609;
t47 = pkin(5) * t269 - pkin(10) * t208 + t64;
t207 = -t270 * t368 + t370 * t609;
t51 = pkin(10) * t207 + t65;
t15 = -t373 * t51 + t377 * t47;
t16 = t373 * t47 + t377 * t51;
t131 = t207 * t377 - t208 * t373;
t132 = t207 * t373 + t208 * t377;
t70 = -t374 * t171 - t193 * t476 + t210 * t378 - t215 * t477;
t118 = -pkin(4) * t609 - t129;
t216 = -t374 * t587 + t378 * t496;
t319 = -t371 * t489 - t492;
t400 = t166 * (mrSges(5,1) * t374 + mrSges(5,2) * t378);
t220 = t246 * t374 - t378 * t498;
t394 = -g(1) * t220 - g(2) * t216 - g(3) * t269;
t390 = t319 * pkin(2);
t63 = -pkin(4) * t289 - t70;
t386 = t371 * t327;
t383 = -mrSges(5,1) + t650;
t351 = Ifges(3,4) * t456;
t344 = Ifges(3,3) * t354;
t343 = Ifges(4,3) * t354;
t330 = -t355 - t540;
t322 = -pkin(8) * t499 + t358;
t321 = (-mrSges(3,1) * t379 + mrSges(3,2) * t375) * t369;
t320 = -t371 * t490 + t488;
t318 = -t371 * t492 - t489;
t312 = t443 * t374;
t307 = -pkin(8) * t455 + t353;
t306 = t483 * qJD(1);
t305 = t352 - t432;
t304 = -mrSges(3,2) * t431 + t428;
t303 = mrSges(3,1) * t431 - t429;
t300 = t403 * t374;
t299 = t328 * t374;
t276 = Ifges(3,1) * t457 + Ifges(3,5) * t431 + t351;
t275 = t509 + (Ifges(3,6) * t371 + (t379 * Ifges(3,2) + t531) * t369) * qJD(1);
t263 = -t361 * t500 + t302;
t261 = -mrSges(4,2) * t431 - t519;
t249 = Ifges(5,4) * t254;
t245 = t376 * t386 - t380 * t389;
t242 = -t376 * t389 - t380 * t386;
t223 = mrSges(4,1) * t287 + mrSges(4,2) * t288;
t221 = t246 * t378 + t374 * t498;
t212 = mrSges(4,1) * t354 - mrSges(4,3) * t233;
t211 = -mrSges(4,2) * t354 + mrSges(4,3) * t232;
t206 = Ifges(4,1) * t288 + t431 * Ifges(4,5) - t630;
t205 = t431 * Ifges(4,6) + t530 - t629;
t189 = -mrSges(5,2) * t283 - mrSges(5,3) * t254;
t156 = t201 * t370 + t289 * t368;
t155 = -t201 * t368 + t289 * t370;
t146 = t255 * Ifges(5,1) - t249 + t643;
t144 = t255 * Ifges(5,5) - t254 * Ifges(5,6) + t628;
t143 = t221 * t366 - t245 * t365;
t142 = -t221 * t365 - t245 * t366;
t115 = -mrSges(5,2) * t230 - mrSges(5,3) * t152;
t88 = mrSges(7,1) * t253 - mrSges(7,3) * t122;
t87 = -mrSges(7,2) * t253 + mrSges(7,3) * t625;
t86 = -pkin(5) * t506 + t108;
t83 = mrSges(5,1) * t152 + mrSges(5,2) * t151;
t82 = -pkin(5) * t207 + t118;
t74 = t151 * Ifges(5,4) - t152 * Ifges(5,2) + t230 * Ifges(5,6);
t73 = -pkin(5) * t436 + t93;
t62 = -mrSges(7,1) * t625 + mrSges(7,2) * t122;
t60 = -qJD(6) * t132 + t155 * t377 - t156 * t373;
t59 = qJD(6) * t131 + t155 * t373 + t156 * t377;
t46 = -pkin(5) * t155 + t63;
t25 = -mrSges(7,2) * t150 + mrSges(7,3) * t35;
t24 = mrSges(7,1) * t150 - mrSges(7,3) * t34;
t21 = -pkin(5) * t109 + t30;
t20 = pkin(10) * t155 + t27;
t17 = pkin(5) * t200 - pkin(10) * t156 + t26;
t4 = -qJD(6) * t16 + t17 * t377 - t20 * t373;
t3 = qJD(6) * t15 + t17 * t373 + t20 * t377;
t8 = [(-Ifges(5,4) * t548 + t647 / 0.2e1 + Ifges(6,3) * t550 - t642 / 0.2e1 - Ifges(5,2) * t551 + t644 / 0.2e1 + t633 + t634 + t649) * t200 + (-t318 * mrSges(3,1) + t626 * mrSges(3,2) + (t311 * t469 + mrSges(2,2)) * t380 + (t364 * t469 - t468) * t376 - (t383 + t592) * t217 + (-t416 - t582) * t242 - t391 * t216 - (pkin(3) * t475 + mrSges(4,1)) * t587) * g(1) + (t343 / 0.2e1 + t344 / 0.2e1 + qJD(1) * t410 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t354 + t584) * t371 + (t1 * t131 - t13 * t59 - t132 * t2 + t14 * t60) * mrSges(7,3) + (Ifges(7,5) * t59 + Ifges(7,6) * t60) * t552 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t561 + (-mrSges(5,3) * t38 + 0.2e1 * t571) * t270 + (Ifges(6,5) * t156 + Ifges(6,6) * t155) * t550 + (Ifges(6,5) * t208 + Ifges(6,6) * t207) * t558 + m(3) * (t250 * t483 + t251 * t322 - t305 * t308 + t306 * t307) + (Ifges(6,1) * t208 + Ifges(6,4) * t207) * t567 + t182 * (Ifges(6,1) * t156 + Ifges(6,4) * t155) / 0.2e1 + t483 * (-mrSges(3,2) * t354 + mrSges(3,3) * t309) + (-t10 * t208 + t11 * t207 + t155 * t50 - t156 * t49) * mrSges(6,3) + t410 * t537 + t59 * t572 + t60 * t574 + t208 * t576 + t207 * t577 + t132 * t580 + t131 * t581 + t330 * t437 + m(4) * (t105 * t202 + t106 * t203 - t170 * t175 + t171 * t176 + t271 * t330) + m(5) * (t107 * t70 + t108 * t69 + t129 * t38 + t130 * t37 + t166 * t170 + t192 * t97) + m(6) * (t10 * t64 + t11 * t65 + t118 * t30 + t26 * t49 + t27 * t50 + t63 * t93) + m(7) * (t1 * t16 + t13 * t4 + t14 * t3 + t15 * t2 + t21 * t82 + t46 * t73) + (Ifges(7,4) * t132 + Ifges(7,2) * t131) * t578 + (Ifges(7,4) * t59 + Ifges(7,2) * t60) * t565 + (t271 * mrSges(4,1) - t106 * mrSges(4,3) - Ifges(4,4) * t233 - Ifges(4,2) * t232 - Ifges(4,6) * t354 + t585) * t609 + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t579 + (Ifges(7,1) * t59 + Ifges(7,4) * t60) * t563 + (mrSges(4,2) * t271 - mrSges(4,3) * t105 + Ifges(4,1) * t233 + Ifges(4,4) * t232 + Ifges(4,5) * t354) * t298 + t436 * (Ifges(6,4) * t156 + Ifges(6,2) * t155) / 0.2e1 + (Ifges(6,4) * t208 + Ifges(6,2) * t207) * t568 - t601 * t170 + t156 * t569 + Ifges(2,3) * qJDD(1) + t322 * (mrSges(3,1) * t354 - mrSges(3,3) * t310) + t307 * t304 - t308 * t303 + t171 * t261 + t97 * t420 + t203 * t211 + t202 * t212 + t30 * (-mrSges(6,1) * t207 + mrSges(6,2) * t208) + t69 * t189 + t70 * t190 + t192 * t83 + t155 * t91 / 0.2e1 + t93 * (-mrSges(6,1) * t155 + mrSges(6,2) * t156) + t27 * t140 + t26 * t141 + t129 * t114 + t130 * t115 + t21 * (-mrSges(7,1) * t131 + mrSges(7,2) * t132) + t63 * t123 + t118 * t61 + t4 * t88 + (-t37 * mrSges(5,3) - Ifges(5,6) * t554 + Ifges(7,6) * t578 + Ifges(7,5) * t579 + Ifges(6,3) * t558 - Ifges(5,2) * t559 - Ifges(5,4) * t560 + Ifges(7,3) * t561 + Ifges(6,5) * t567 + Ifges(6,6) * t568 - t74 / 0.2e1 - t11 * mrSges(6,2) + t10 * mrSges(6,1) + t631 / 0.2e1 + t636) * t269 + (-t176 * mrSges(4,3) - Ifges(4,4) * t544 + Ifges(5,5) * t548 + Ifges(5,6) * t551 + t314 * mrSges(4,1) + t629 / 0.2e1 + t628 / 0.2e1 + t144 / 0.2e1 - t205 / 0.2e1 + t618 - t617) * t289 - (-t175 * mrSges(4,3) + Ifges(4,1) * t544 + t314 * mrSges(4,2) - t630 / 0.2e1 + t206 / 0.2e1) * t290 + ((mrSges(3,1) * t309 - mrSges(3,2) * t310 + (m(3) * t540 - t321) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t250 + Ifges(3,4) * t310 + Ifges(3,2) * t309 + Ifges(3,6) * t354) * t379 + (-mrSges(3,3) * t251 + Ifges(3,1) * t310 + Ifges(3,4) * t309 + Ifges(3,5) * t354) * t375 + ((-t305 * mrSges(3,3) + t276 / 0.2e1 + Ifges(3,5) * t537) * t379 + (-t306 * mrSges(3,3) - t275 / 0.2e1 - t509 / 0.2e1 + (m(4) * t314 + t223) * pkin(2)) * t375 + (t371 * t411 / 0.2e1 + (t379 * (Ifges(3,4) * t379 - Ifges(3,2) * t375) / 0.2e1 - t590) * t369) * qJD(1)) * qJD(2) + (g(1) * t380 + g(2) * t376) * (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3))) * t369 + (-m(4) * t435 - t246 * mrSges(4,1) + t376 * mrSges(2,2) - t143 * mrSges(7,1) - t142 * mrSges(7,2) - t320 * mrSges(3,1) - t319 * mrSges(3,2) + t468 * t380 + t383 * t221 + t582 * t245 + t391 * t220 + t475 * (-t246 * pkin(3) - t435)) * g(2) + t15 * t24 + t16 * t25 + (Ifges(5,1) * t548 - t646 + t643 / 0.2e1 + Ifges(5,4) * t551 + t645 + t146 / 0.2e1) * t201 + t46 * t62 + t65 * t71 + t64 * t72 + t73 * (-mrSges(7,1) * t60 + mrSges(7,2) * t59) + t82 * t12 + t3 * t87; -(t637 + t648) * t504 + (-m(4) * t355 - m(5) * t423 + t627 * t298 + t623 * t485 - t586 * t609 + t321) * g(3) + (-Ifges(4,2) * t288 + t206 - t630) * t287 / 0.2e1 + (-m(4) * t424 - mrSges(3,1) * t626 - mrSges(3,2) * t318 - t475 * (t242 * pkin(3) + t424) + t586 * t242 + t583 * t587) * g(2) + (Ifges(7,5) * t235 + Ifges(7,6) * t236) * t552 + (t400 + t49 * (mrSges(6,1) * t374 - mrSges(6,3) * t494) + t50 * (-mrSges(6,2) * t374 - mrSges(6,3) * t500)) * qJD(4) + t21 * (mrSges(7,1) * t299 - mrSges(7,2) * t300) + (-Ifges(7,4) * t300 - Ifges(7,2) * t299 - Ifges(7,6) * t378) * t578 + (-Ifges(7,1) * t300 - Ifges(7,4) * t299 - Ifges(7,5) * t378) * t579 + (-Ifges(7,5) * t300 - Ifges(7,6) * t299 - Ifges(7,3) * t378) * t561 + t2 * (-mrSges(7,1) * t378 + mrSges(7,3) * t300) + (t182 * (Ifges(6,5) * t374 + t378 * t414) + t436 * (Ifges(6,6) * t374 + t378 * t412) + t283 * t409 + t255 * t415) * qJD(4) / 0.2e1 + t1 * (mrSges(7,2) * t378 - mrSges(7,3) * t299) + (t175 * t194 - t176 * t195 - t314 * t433 + (t105 * t511 + t106 * t510) * pkin(2)) * m(4) + (-t452 - t126) * t190 + t584 - t223 * t433 - t175 * t519 + (-t112 + t452) * t123 + (-t522 + t633) * t477 + (t428 - t304) * t305 + t343 + t344 + t205 * t544 + (Ifges(5,5) * t374 + Ifges(5,6) * t378) * t554 + t374 * t571 + t235 * t572 + t138 * t573 + t236 * t574 + t137 * t575 + t495 * t576 - t300 * t580 - t299 * t581 + (-t453 - t127) * t189 + (mrSges(6,1) * t93 - mrSges(6,3) * t50 + Ifges(6,4) * t556 + Ifges(6,2) * t557 + Ifges(6,6) * t551 + t570) * t213 + t211 * t458 + t212 * t459 + t176 * t518 + (Ifges(7,4) * t138 + Ifges(7,2) * t137) * t566 - t41 * t501 / 0.2e1 + t11 * (mrSges(6,2) * t378 - mrSges(6,3) * t501) + (t429 + t303) * t306 + (-t92 / 0.2e1 + Ifges(6,5) * t551 + t49 * mrSges(6,3) + Ifges(6,1) * t556 + Ifges(6,4) * t557 - t93 * mrSges(6,2)) * t214 + (Ifges(7,4) * t235 + Ifges(7,2) * t236) * t565 - ((-Ifges(3,2) * t457 + t276 + t351) * t379 + t431 * t411) * t482 / 0.2e1 + t10 * (-mrSges(6,1) * t378 - mrSges(6,3) * t495) + (Ifges(6,3) * t374 + t378 * t408) * t478 / 0.2e1 - t413 * t478 / 0.2e1 + (-t107 * t126 - t108 * t127 - t166 * t194 + t363 * t97) * m(5) + (t10 * t263 + t11 * t264 - t112 * t93 + t49 * t612 + t50 * t611) * m(6) + (t615 * t374 + ((-t107 * t378 - t108 * t374) * qJD(4) + t593) * m(5) + t378 * t115 + t597 * m(6)) * t361 + (Ifges(7,1) * t235 + Ifges(7,4) * t236) * t563 + (Ifges(7,1) * t138 + Ifges(7,4) * t137) * t564 + (-Ifges(6,3) * t378 + t374 * t408) * t558 + (Ifges(5,2) * t378 + t529) * t559 + (Ifges(5,1) * t374 + t528) * t560 + (-Ifges(6,5) * t378 + t374 * t414) * t567 + (-Ifges(6,6) * t378 + t374 * t412) * t568 + t275 * t457 / 0.2e1 + t378 * t74 / 0.2e1 + t363 * t83 + t312 * t12 - t195 * t261 + t263 * t72 + t264 * t71 - t97 * t419 + t153 * t24 + t154 * t25 + t287 * t400 + t590 * qJD(1) ^ 2 * t369 ^ 2 - t431 * (-Ifges(4,5) * t287 - Ifges(4,6) * t288) / 0.2e1 + (Ifges(5,3) * t288 - t287 * t409) * t546 + (Ifges(5,5) * t288 - t287 * t415) * t549 + (Ifges(5,6) * t288 - t287 * t413) * t550 - t314 * (mrSges(4,1) * t288 - mrSges(4,2) * t287) + (-t107 * t503 + t593) * mrSges(5,3) + (t503 / 0.2e1 + t444) * t146 - (-Ifges(4,1) * t287 + t144 - t530) * t288 / 0.2e1 + (-m(4) * t390 - mrSges(3,1) * t319 + mrSges(3,2) * t320 - t475 * (t245 * pkin(3) + t390) + t586 * t245 - t583 * t246) * g(1) + (Ifges(7,5) * t138 + Ifges(7,6) * t137) * t553 + t370 * t92 * t444 + t597 * t417 + t598 * t476 + t601 * t194 + (-mrSges(7,2) * t596 + mrSges(7,3) * t602) * t14 + (-mrSges(7,1) * t602 + mrSges(7,2) * t603) * t73 + (mrSges(7,1) * t596 - mrSges(7,3) * t603) * t13 - t631 * t378 / 0.2e1 - t288 * t618 + t610 * t62 + t611 * t140 + t612 * t141 + t288 * t617 + t620 * t87 + t621 * t88 + (t1 * t154 + t13 * t621 + t14 * t620 + t153 * t2 + t21 * t312 + t610 * t73) * m(7); t602 * t88 + (t287 * t189 - t12 + (t189 + t595) * qJD(4) - t615) * t378 + t437 + t601 * t288 + (t115 + t591 * (t62 - t604) + t594) * t374 + t603 * t87 - t299 * t24 - t300 * t25 + t287 * t261 - t213 * t141 - t214 * t140 + (-t1 * t300 + t13 * t602 + t14 * t603 - t2 * t299 - t21 * t378 + t596 * t73) * m(7) + (-t213 * t49 - t214 * t50 + t93 * t504 - t30 * t378 + t406 * t374 + (t374 * t93 + t378 * t405) * qJD(4)) * m(6) + (-t166 * t288 + t37 * t374 + t378 * t38 + t591 * (-t107 * t374 + t108 * t378)) * m(5) + (t175 * t288 + t176 * t287 + t271) * m(4) + (-t371 * g(3) + (-g(1) * t376 + g(2) * t380) * t369) * t469; (-Ifges(7,5) * t315 - Ifges(7,6) * t316) * t552 + (-Ifges(7,1) * t315 - Ifges(7,4) * t316) * t563 + (-Ifges(7,4) * t315 - Ifges(7,2) * t316) * t565 + (Ifges(7,1) * t158 + Ifges(7,4) * t157) * t564 + t145 * t548 + t585 + (Ifges(7,4) * t158 + Ifges(7,2) * t157) * t566 - t315 * t572 + t158 * t573 - t316 * t574 + t157 * t575 + t368 * t576 + t370 * t577 + t328 * t580 + (t146 - t249) * t550 + (Ifges(7,5) * t158 + Ifges(7,6) * t157) * t553 + (-pkin(4) * t30 + qJ(5) * t406 + qJD(5) * t405 - t108 * t93 - t49 * t76 - t50 * t77) * m(6) + (Ifges(7,4) * t328 - Ifges(7,2) * t403) * t578 + (Ifges(7,1) * t328 - Ifges(7,4) * t403) * t579 + (Ifges(7,5) * t328 - Ifges(7,6) * t403) * t561 + t21 * (mrSges(7,1) * t403 + mrSges(7,2) * t328) - t403 * t581 + (-t1 * t403 + t13 * t605 - t14 * t606 - t2 * t328) * mrSges(7,3) + (Ifges(6,5) * t368 + Ifges(6,6) * t370) * t558 + (Ifges(6,1) * t368 + t526) * t567 + (Ifges(6,2) * t370 + t527) * t568 + t505 * t569 - t362 * t12 + t272 * t24 + t273 * t25 + t30 * t418 - t107 * t189 - t77 * t140 - t76 * t141 + (t220 * t588 + t221 * t391) * g(1) + (t216 * t588 + t217 * t391) * g(2) + (-t49 * t505 - t50 * t506 + t406) * mrSges(6,3) + t594 * qJ(5) + t595 * qJD(5) + (-Ifges(5,2) * t550 - Ifges(5,6) * t546 - t634 + t637) * t255 + (t269 * t600 + t270 * t599 + t420) * g(3) + t604 * t108 + (mrSges(7,1) * t606 - mrSges(7,2) * t605) * t73 + t613 * t88 + t614 * t87 + (t1 * t273 + t13 * t613 + t14 * t614 + t2 * t272 - t21 * t362 - t73 * t86) * m(7) + (-Ifges(5,1) * t549 - Ifges(5,5) * t546 - t408 * t551 - t412 * t557 - t414 * t556 + t417 * t93 + t598 + t645) * t254 - pkin(4) * t61 + (-t520 + t619) * t549 - t86 * t62; t122 * t88 - t625 * t87 - t436 * t140 + t182 * t141 + t12 + t61 + (t122 * t13 - t14 * t625 + t21 + t394) * m(7) + (t182 * t49 - t436 * t50 + t30 + t394) * m(6); -t73 * (mrSges(7,1) * t122 + mrSges(7,2) * t625) + t53 * t563 + (Ifges(7,5) * t625 - Ifges(7,6) * t122) * t553 + (Ifges(7,1) * t625 - t525) * t564 - t13 * t87 + t14 * t88 - g(1) * (mrSges(7,1) * t142 - mrSges(7,2) * t143) - g(2) * ((-t217 * t365 - t242 * t366) * mrSges(7,1) + (-t217 * t366 + t242 * t365) * mrSges(7,2)) - g(3) * ((-t270 * t365 + t366 * t609) * mrSges(7,1) + (-t270 * t366 - t365 * t609) * mrSges(7,2)) + (t122 * t14 + t13 * t625) * mrSges(7,3) + t7 + (-Ifges(7,2) * t122 + t116 + t54) * t566 + t636;];
tau  = t8;
