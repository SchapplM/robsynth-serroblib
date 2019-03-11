% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:05
% EndTime: 2019-03-09 09:39:56
% DurationCPUTime: 32.27s
% Computational Cost: add. (15894->903), mult. (37825->1179), div. (0->0), fcn. (29539->14), ass. (0->423)
t402 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t574 = -m(7) - m(6);
t575 = m(5) + m(4);
t595 = -t575 + t574;
t596 = -qJ(3) * t595 - mrSges(3,2) + mrSges(4,3);
t324 = sin(pkin(11));
t326 = cos(pkin(11));
t331 = sin(qJ(5));
t335 = cos(qJ(5));
t266 = t324 * t335 + t326 * t331;
t325 = sin(pkin(6));
t332 = sin(qJ(2));
t459 = t325 * t332;
t345 = t266 * t459;
t196 = qJD(1) * t345;
t434 = qJD(5) * t335;
t435 = qJD(5) * t331;
t255 = -t324 * t434 - t326 * t435;
t557 = t255 - t196;
t443 = qJD(1) * t325;
t415 = t332 * t443;
t394 = t326 * t415;
t423 = t324 * t459;
t395 = t331 * t423;
t195 = -qJD(1) * t395 + t335 * t394;
t256 = -t324 * t435 + t326 * t434;
t556 = t256 + t195;
t330 = sin(qJ(6));
t334 = cos(qJ(6));
t383 = -mrSges(7,1) * t334 + mrSges(7,2) * t330;
t567 = m(7) * pkin(5) + mrSges(6,1) - t383;
t292 = pkin(2) * t415;
t336 = cos(qJ(2));
t368 = -qJ(3) * t336 + qJ(4) * t332;
t198 = t368 * t443 + t292;
t414 = t336 * t443;
t295 = pkin(8) * t414;
t327 = cos(pkin(6));
t442 = qJD(1) * t327;
t430 = pkin(1) * t442;
t244 = t332 * t430 + t295;
t205 = pkin(3) * t414 + t244;
t122 = -t198 * t324 + t326 * t205;
t460 = t324 * t332;
t350 = (pkin(4) * t336 - pkin(9) * t460) * t325;
t106 = qJD(1) * t350 + t122;
t123 = t326 * t198 + t324 * t205;
t111 = pkin(9) * t394 + t123;
t495 = pkin(2) + qJ(4);
t494 = -pkin(9) - t495;
t271 = t494 * t324;
t272 = t494 * t326;
t186 = t271 * t331 - t335 * t272;
t566 = -qJD(4) * t266 - qJD(5) * t186 - t331 * t106 - t335 * t111;
t323 = pkin(11) + qJ(5);
t318 = sin(t323);
t319 = cos(t323);
t385 = t324 * mrSges(5,1) + t326 * mrSges(5,2);
t594 = -t318 * t567 - t319 * t402 - t385;
t413 = qJD(2) * t459;
t457 = t325 * t336;
t247 = qJD(1) * t413 - qJDD(1) * t457;
t431 = qJDD(1) * t327;
t304 = qJDD(2) + t431;
t173 = t247 * t326 - t304 * t324;
t174 = t247 * t324 + t304 * t326;
t306 = qJD(2) + t442;
t210 = -t306 * t326 + t324 * t414;
t347 = -t306 * t324 - t326 * t414;
t582 = t210 * t331 + t335 * t347;
t71 = qJD(5) * t582 + t173 * t331 + t174 * t335;
t527 = t71 / 0.2e1;
t581 = -t335 * t210 + t331 * t347;
t72 = -qJD(5) * t581 + t173 * t335 - t174 * t331;
t526 = t72 / 0.2e1;
t439 = qJD(2) * t336;
t248 = (qJD(1) * t439 + qJDD(1) * t332) * t325;
t232 = qJDD(5) + t248;
t509 = t232 / 0.2e1;
t593 = -mrSges(3,1) + mrSges(4,2);
t592 = -m(5) * qJ(4) - mrSges(5,3) - mrSges(6,3);
t591 = -pkin(10) * t414 + t566;
t297 = t336 * t430;
t316 = pkin(4) * t326 + pkin(3);
t359 = (-pkin(8) - t316) * t459;
t161 = qJD(1) * t359 + t297;
t590 = t556 * pkin(5) - pkin(10) * t557 + qJD(3) - t161;
t70 = qJDD(6) - t72;
t528 = t70 / 0.2e1;
t278 = qJD(5) + t415;
t117 = t278 * t330 + t334 * t581;
t40 = -qJD(6) * t117 + t232 * t334 - t330 * t71;
t534 = t40 / 0.2e1;
t116 = t278 * t334 - t330 * t581;
t39 = qJD(6) * t116 + t232 * t330 + t334 * t71;
t535 = t39 / 0.2e1;
t522 = pkin(3) + pkin(8);
t149 = -t306 * t495 + t415 * t522 + qJD(3) - t297;
t404 = -qJ(3) * t332 - pkin(1);
t175 = (-t336 * t495 + t404) * t443;
t99 = t326 * t149 - t175 * t324;
t79 = pkin(4) * t415 + pkin(9) * t210 + t99;
t100 = t324 * t149 + t326 * t175;
t81 = pkin(9) * t347 + t100;
t34 = t331 * t79 + t335 * t81;
t32 = pkin(10) * t278 + t34;
t284 = t306 * qJ(3);
t163 = t284 + qJD(4) + t205;
t124 = -pkin(4) * t347 + t163;
t57 = -pkin(5) * t582 - pkin(10) * t581 + t124;
t13 = -t32 * t330 + t334 * t57;
t504 = pkin(1) * t327;
t429 = qJD(2) * t504;
t396 = qJD(1) * t429;
t424 = pkin(1) * t431;
t159 = -pkin(8) * t247 + t332 * t424 + t336 * t396;
t128 = -t304 * qJ(3) - t306 * qJD(3) - t159;
t107 = -pkin(3) * t247 + qJDD(4) - t128;
t80 = -pkin(4) * t173 + t107;
t19 = -pkin(5) * t72 - pkin(10) * t71 + t80;
t438 = qJD(3) * t332;
t365 = -qJD(4) * t336 - t438;
t468 = qJ(3) * t248;
t469 = pkin(1) * qJDD(1);
t103 = -t468 + t495 * t247 + (qJD(1) * t365 - t469) * t325;
t307 = pkin(8) * t459;
t160 = -qJD(2) * t295 - qJDD(1) * t307 - t332 * t396 + t336 * t424;
t341 = qJDD(3) - t160;
t98 = pkin(3) * t248 - qJD(4) * t306 - t304 * t495 + t341;
t51 = -t103 * t324 + t326 * t98;
t36 = pkin(4) * t248 - pkin(9) * t174 + t51;
t52 = t326 * t103 + t324 * t98;
t42 = pkin(9) * t173 + t52;
t7 = t331 * t36 + t335 * t42 + t79 * t434 - t435 * t81;
t5 = pkin(10) * t232 + t7;
t1 = qJD(6) * t13 + t19 * t330 + t334 * t5;
t14 = t32 * t334 + t330 * t57;
t2 = -qJD(6) * t14 + t19 * t334 - t330 * t5;
t546 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t70;
t588 = t546 + mrSges(6,1) * t80 + Ifges(7,5) * t535 + Ifges(7,6) * t534 + Ifges(7,3) * t528 + t9 / 0.2e1 + (-t509 - t232 / 0.2e1) * Ifges(6,6) + (-t526 - t72 / 0.2e1) * Ifges(6,2) + (-t527 - t71 / 0.2e1) * Ifges(6,4);
t506 = t278 / 0.2e1;
t512 = t581 / 0.2e1;
t514 = t582 / 0.2e1;
t137 = qJD(6) - t582;
t516 = t137 / 0.2e1;
t518 = t117 / 0.2e1;
t520 = t116 / 0.2e1;
t43 = t117 * Ifges(7,5) + t116 * Ifges(7,6) + t137 * Ifges(7,3);
t487 = Ifges(6,4) * t581;
t77 = Ifges(6,2) * t582 + t278 * Ifges(6,6) + t487;
t552 = -t77 / 0.2e1 + t43 / 0.2e1;
t587 = -Ifges(6,4) * t512 + Ifges(7,5) * t518 - Ifges(6,2) * t514 - Ifges(6,6) * t506 + Ifges(7,6) * t520 + Ifges(7,3) * t516 + t552;
t586 = -t332 / 0.2e1;
t585 = mrSges(6,2) * t124;
t187 = t271 * t335 + t272 * t331;
t436 = qJD(4) * t326;
t437 = qJD(4) * t324;
t565 = -qJD(5) * t187 + (-t106 - t436) * t335 + (t111 + t437) * t331;
t147 = -mrSges(5,1) * t347 - mrSges(5,2) * t210;
t82 = -mrSges(6,1) * t582 + mrSges(6,2) * t581;
t584 = -t147 - t82;
t364 = t592 + t593;
t382 = mrSges(7,1) * t330 + mrSges(7,2) * t334;
t543 = -t364 + t382;
t580 = -mrSges(6,1) * t124 - mrSges(7,1) * t13 + mrSges(7,2) * t14;
t136 = Ifges(6,4) * t582;
t78 = Ifges(6,1) * t581 + t278 * Ifges(6,5) + t136;
t579 = Ifges(6,1) * t512 + Ifges(6,4) * t514 + Ifges(6,5) * t506 + t78 / 0.2e1;
t578 = t80 * mrSges(6,2) + 0.2e1 * Ifges(6,1) * t527 + 0.2e1 * Ifges(6,4) * t526 + 0.2e1 * Ifges(6,5) * t509;
t511 = t173 / 0.2e1;
t510 = t174 / 0.2e1;
t508 = t248 / 0.2e1;
t572 = Ifges(3,5) - Ifges(4,4);
t571 = Ifges(4,5) - Ifges(3,6);
t12 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t60 = mrSges(6,1) * t232 - mrSges(6,3) * t71;
t570 = t12 - t60;
t265 = t324 * t331 - t335 * t326;
t321 = t324 * pkin(4);
t314 = qJ(3) + t321;
t165 = pkin(5) * t266 + pkin(10) * t265 + t314;
t109 = t165 * t334 - t187 * t330;
t569 = qJD(6) * t109 + t330 * t590 + t334 * t591;
t110 = t165 * t330 + t187 * t334;
t568 = -qJD(6) * t110 - t330 * t591 + t334 * t590;
t564 = pkin(5) * t414 - t565;
t491 = mrSges(6,3) * t581;
t121 = mrSges(6,1) * t278 - t491;
t62 = -mrSges(7,1) * t116 + mrSges(7,2) * t117;
t563 = -t62 + t121;
t155 = -t196 * t330 + t334 * t414;
t432 = qJD(6) * t334;
t358 = -t330 * t255 + t265 * t432;
t560 = t155 - t358;
t156 = t196 * t334 + t330 * t414;
t433 = qJD(6) * t330;
t449 = t334 * t255;
t357 = t265 * t433 + t449;
t559 = t156 - t357;
t398 = mrSges(3,3) * t415;
t400 = mrSges(4,1) * t415;
t558 = t306 * t593 + t398 + t400;
t253 = -t324 * t327 - t326 * t457;
t254 = -t324 * t457 + t326 * t327;
t158 = t253 * t331 + t254 * t335;
t286 = t334 * t459;
t134 = -t158 * t330 + t286;
t131 = -mrSges(5,2) * t248 + mrSges(5,3) * t173;
t132 = mrSges(5,1) * t248 - mrSges(5,3) * t174;
t555 = t324 * t131 + t326 * t132;
t20 = mrSges(7,1) * t70 - mrSges(7,3) * t39;
t21 = -mrSges(7,2) * t70 + mrSges(7,3) * t40;
t554 = -t330 * t20 + t334 * t21;
t374 = t324 * t52 + t326 * t51;
t553 = t1 * t334 - t2 * t330;
t490 = Ifges(3,4) * t332;
t551 = pkin(1) * (mrSges(3,1) * t332 + mrSges(3,2) * t336) + (Ifges(3,1) * t336 - t490) * t586;
t203 = -t284 - t244;
t399 = mrSges(4,1) * t414;
t239 = -mrSges(4,3) * t306 - t399;
t549 = -m(4) * t203 - t239;
t548 = t385 + t596;
t8 = -qJD(5) * t34 - t331 * t42 + t335 * t36;
t547 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t71 + Ifges(6,6) * t72 + Ifges(6,3) * t232;
t503 = pkin(1) * t336;
t416 = -pkin(2) - t503;
t182 = pkin(3) * t459 + t307 + (-qJ(4) + t416) * t327;
t445 = pkin(2) * t457 + qJ(3) * t459;
t194 = (-qJ(4) * t336 - pkin(1)) * t325 - t445;
t119 = t324 * t182 + t326 * t194;
t101 = pkin(9) * t253 + t119;
t118 = t326 * t182 - t194 * t324;
t91 = pkin(4) * t459 - pkin(9) * t254 + t118;
t493 = t335 * t101 + t331 * t91;
t294 = pkin(2) * t413;
t152 = t294 + (qJD(2) * t368 + t365) * t325;
t313 = t332 * t504;
t184 = -qJD(4) * t327 + (t457 * t522 + t313) * qJD(2);
t104 = -t152 * t324 + t326 * t184;
t85 = qJD(2) * t350 + t104;
t105 = t326 * t152 + t324 * t184;
t393 = t326 * t413;
t90 = pkin(9) * t393 + t105;
t18 = -qJD(5) * t493 - t331 * t90 + t335 * t85;
t33 = -t331 * t81 + t335 * t79;
t545 = t33 * mrSges(6,1) - t34 * mrSges(6,2);
t507 = -t278 / 0.2e1;
t513 = -t581 / 0.2e1;
t544 = Ifges(6,1) * t513 + Ifges(6,5) * t507;
t542 = t594 - t596;
t541 = t265 * t8 - t266 * t7 - t33 * t557 - t34 * t556;
t515 = -t582 / 0.2e1;
t517 = -t137 / 0.2e1;
t519 = -t117 / 0.2e1;
t521 = -t116 / 0.2e1;
t540 = -Ifges(7,5) * t519 + Ifges(6,2) * t515 + Ifges(6,6) * t507 - Ifges(7,6) * t521 - Ifges(7,3) * t517;
t10 = t39 * Ifges(7,4) + t40 * Ifges(7,2) + t70 * Ifges(7,6);
t538 = t10 / 0.2e1;
t11 = t39 * Ifges(7,1) + t40 * Ifges(7,4) + t70 * Ifges(7,5);
t537 = t11 / 0.2e1;
t486 = Ifges(7,4) * t117;
t44 = Ifges(7,2) * t116 + Ifges(7,6) * t137 + t486;
t532 = -t44 / 0.2e1;
t531 = t44 / 0.2e1;
t113 = Ifges(7,4) * t116;
t45 = Ifges(7,1) * t117 + Ifges(7,5) * t137 + t113;
t530 = -t45 / 0.2e1;
t529 = t45 / 0.2e1;
t523 = Ifges(5,1) * t510 + Ifges(5,4) * t511 + Ifges(5,5) * t508;
t505 = pkin(1) * t325;
t6 = -pkin(5) * t232 - t8;
t500 = t265 * t6;
t496 = Ifges(3,4) + Ifges(4,6);
t492 = mrSges(6,3) * t582;
t489 = Ifges(5,4) * t324;
t488 = Ifges(5,4) * t326;
t485 = Ifges(7,4) * t330;
t484 = Ifges(7,4) * t334;
t483 = Ifges(4,6) * t332;
t482 = Ifges(4,6) * t336;
t481 = Ifges(5,3) * t336;
t480 = t582 * Ifges(6,6);
t479 = t581 * Ifges(6,5);
t478 = t347 * Ifges(5,6);
t477 = t210 * Ifges(5,5);
t474 = t278 * Ifges(6,3);
t467 = t582 * t330;
t466 = t582 * t334;
t464 = t265 * t330;
t463 = t265 * t334;
t462 = t324 * (-Ifges(5,1) * t210 + Ifges(5,4) * t347 + Ifges(5,5) * t415);
t333 = sin(qJ(1));
t458 = t325 * t333;
t337 = cos(qJ(1));
t456 = t325 * t337;
t455 = t326 * (-Ifges(5,4) * t210 + Ifges(5,2) * t347 + Ifges(5,6) * t415);
t453 = t326 * t332;
t452 = t332 * t333;
t451 = t332 * t337;
t450 = t333 * t336;
t448 = t336 * t337;
t298 = t336 * t429;
t317 = t327 * qJD(3);
t446 = t298 + t317;
t263 = pkin(8) * t457 + t313;
t444 = t337 * pkin(1) + pkin(8) * t458;
t440 = qJD(1) ^ 2 * t325 ^ 2;
t428 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t427 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t422 = t330 * t459;
t421 = t330 * t532;
t260 = -t327 * t452 + t448;
t418 = t260 * pkin(2) + t444;
t224 = -t327 * qJ(3) - t263;
t412 = t325 * t439;
t30 = -t72 * mrSges(6,1) + t71 * mrSges(6,2);
t405 = t432 / 0.2e1;
t403 = -t333 * pkin(1) + pkin(8) * t456;
t192 = t248 * mrSges(4,1) + t304 * mrSges(4,2);
t112 = -t173 * mrSges(5,1) + t174 * mrSges(5,2);
t401 = t522 * t459;
t397 = mrSges(3,3) * t414;
t197 = pkin(3) * t457 - t224;
t258 = t327 * t451 + t450;
t387 = t258 * pkin(2) - t403;
t386 = -mrSges(5,1) * t326 + mrSges(5,2) * t324;
t381 = mrSges(4,2) * t336 - mrSges(4,3) * t332;
t380 = Ifges(5,1) * t324 + t488;
t379 = Ifges(7,1) * t334 - t485;
t378 = Ifges(5,2) * t326 + t489;
t377 = -Ifges(7,2) * t330 + t484;
t376 = Ifges(7,5) * t334 - Ifges(7,6) * t330;
t375 = -t13 * t330 + t14 * t334;
t47 = pkin(10) * t459 + t493;
t150 = -pkin(4) * t253 + t197;
t366 = t335 * t253 - t254 * t331;
t66 = -pkin(5) * t366 - pkin(10) * t158 + t150;
t25 = t330 * t66 + t334 * t47;
t24 = -t330 * t47 + t334 * t66;
t73 = -mrSges(7,2) * t137 + mrSges(7,3) * t116;
t74 = mrSges(7,1) * t137 - mrSges(7,3) * t117;
t373 = -t330 * t74 + t334 * t73;
t370 = t100 * t326 - t324 * t99;
t49 = -t101 * t331 + t335 * t91;
t243 = pkin(8) * t415 - t297;
t245 = -pkin(8) * t413 + t298;
t120 = -mrSges(6,2) * t278 + t492;
t363 = -t120 - t373;
t259 = t327 * t450 + t451;
t328 = -pkin(9) - qJ(4);
t362 = t259 * t321 - t260 * t328 + t316 * t458 + t418;
t135 = t158 * t334 + t422;
t257 = -t327 * t448 + t452;
t180 = -t257 * t318 + t319 * t456;
t178 = t257 * t319 + t318 * t456;
t31 = -pkin(5) * t278 - t33;
t360 = t31 * t382;
t17 = -t101 * t435 + t331 * t85 + t335 * t90 + t91 * t434;
t356 = t163 * t386;
t354 = t332 * (-Ifges(4,2) * t336 + t483);
t353 = t336 * (Ifges(4,3) * t332 - t482);
t344 = -g(1) * t260 - g(2) * t258 - g(3) * t459;
t146 = -pkin(2) * t304 + t341;
t340 = t160 * mrSges(3,1) - t159 * mrSges(3,2) + t146 * mrSges(4,2) - t128 * mrSges(4,3);
t151 = qJD(2) * t359 + t446;
t339 = (-t13 * t334 - t14 * t330) * qJD(6) + t553;
t291 = Ifges(3,4) * t414;
t283 = Ifges(4,1) * t304;
t282 = Ifges(3,3) * t304;
t262 = t327 * t503 - t307;
t261 = (-mrSges(3,1) * t336 + mrSges(3,2) * t332) * t325;
t251 = t259 * pkin(2);
t249 = t257 * pkin(2);
t246 = t263 * qJD(2);
t242 = -qJ(3) * t414 + t292;
t241 = t381 * t443;
t238 = -mrSges(3,2) * t306 + t397;
t231 = Ifges(4,4) * t248;
t230 = Ifges(3,5) * t248;
t229 = Ifges(4,5) * t247;
t228 = Ifges(3,6) * t247;
t227 = t327 * t416 + t307;
t225 = -t445 - t505;
t219 = -t318 * t457 + t319 * t327;
t209 = -t245 - t317;
t208 = (-pkin(2) * t336 + t404) * t443;
t206 = t294 + (-qJ(3) * t439 - t438) * t325;
t204 = -qJD(1) * t401 + t297;
t202 = t306 * Ifges(4,4) + (-t332 * Ifges(4,2) - t482) * t443;
t201 = t306 * Ifges(4,5) + (-t336 * Ifges(4,3) - t483) * t443;
t200 = Ifges(3,1) * t415 + t306 * Ifges(3,5) + t291;
t199 = t306 * Ifges(3,6) + (t336 * Ifges(3,2) + t490) * t443;
t193 = -pkin(2) * t306 + qJD(3) + t243;
t191 = mrSges(4,1) * t247 - mrSges(4,3) * t304;
t183 = -qJD(2) * t401 + t446;
t177 = t259 * t318 + t319 * t458;
t176 = -t259 * t319 + t318 * t458;
t172 = mrSges(5,1) * t415 + mrSges(5,3) * t210;
t171 = -mrSges(5,2) * t415 + mrSges(5,3) * t347;
t154 = -t195 * t334 + t306 * t330;
t153 = t195 * t330 + t306 * t334;
t133 = pkin(2) * t247 - t468 + (-qJD(1) * t438 - t469) * t325;
t130 = t177 * t334 + t260 * t330;
t129 = -t177 * t330 + t260 * t334;
t125 = Ifges(5,3) * t415 - t477 + t478;
t115 = qJD(2) * t395 + qJD(5) * t158 - t335 * t393;
t114 = qJD(2) * t345 + qJD(5) * t366;
t94 = t174 * Ifges(5,4) + t173 * Ifges(5,2) + t248 * Ifges(5,6);
t83 = pkin(5) * t581 - pkin(10) * t582;
t76 = t474 + t479 + t480;
t64 = -qJD(6) * t135 - t114 * t330 + t334 * t412;
t63 = qJD(6) * t134 + t114 * t334 + t330 * t412;
t61 = -mrSges(6,2) * t232 + mrSges(6,3) * t72;
t48 = pkin(5) * t115 - pkin(10) * t114 + t151;
t46 = -pkin(5) * t459 - t49;
t23 = t33 * t334 + t330 * t83;
t22 = -t33 * t330 + t334 * t83;
t16 = -pkin(5) * t412 - t18;
t15 = pkin(10) * t412 + t17;
t4 = -qJD(6) * t25 - t15 * t330 + t334 * t48;
t3 = qJD(6) * t24 + t15 * t334 + t330 * t48;
t26 = [(t253 * t52 - t254 * t51) * mrSges(5,3) + m(7) * (t1 * t25 + t13 * t4 + t14 * t3 + t16 * t31 + t2 * t24 + t46 * t6) + m(5) * (t100 * t105 + t104 * t99 + t107 * t197 + t118 * t51 + t119 * t52 + t163 * t183) + m(4) * (t128 * t224 + t133 * t225 + t146 * t227 + t193 * t246 + t203 * t209 + t206 * t208) + (-t8 * mrSges(6,3) + t578) * t158 + (Ifges(7,5) * t135 + Ifges(7,6) * t134) * t528 + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t516 + (-mrSges(6,3) * t34 - t580 + t587) * t115 + (mrSges(6,3) * t7 - t588) * t366 + (Ifges(7,4) * t135 + Ifges(7,2) * t134) * t534 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t520 + (t230 / 0.2e1 - t228 / 0.2e1 + t282 / 0.2e1 + t283 / 0.2e1 - t231 / 0.2e1 + t229 / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t304 + t428 * t248 + t427 * t247 + t340) * t327 + ((-mrSges(3,1) * t247 - mrSges(3,2) * t248 + (m(3) * t505 - t261) * qJDD(1)) * pkin(1) + (-t128 * mrSges(4,1) + t133 * mrSges(4,2) + t159 * mrSges(3,3) - t571 * t304 + t496 * t248 + (-Ifges(3,2) - Ifges(4,3)) * t247) * t336 + (Ifges(5,5) * t174 + Ifges(5,6) * t173 - t160 * mrSges(3,3) + t146 * mrSges(4,1) - t52 * mrSges(5,2) + t51 * mrSges(5,1) - t133 * mrSges(4,3) + t572 * t304 - t496 * t247 + (Ifges(3,1) + Ifges(4,2) + Ifges(5,3)) * t248 + t547) * t332 + ((t455 / 0.2e1 + t462 / 0.2e1 - t244 * mrSges(3,3) + t203 * mrSges(4,1) - t199 / 0.2e1 + t201 / 0.2e1 + t347 * t378 / 0.2e1 - t210 * t380 / 0.2e1 + t356 - t208 * mrSges(4,2) + t427 * t306 + t370 * mrSges(5,3)) * t332 + (t474 / 0.2e1 + t479 / 0.2e1 + t480 / 0.2e1 + t243 * mrSges(3,3) + t193 * mrSges(4,1) + t200 / 0.2e1 - t202 / 0.2e1 + t125 / 0.2e1 + t76 / 0.2e1 - t100 * mrSges(5,2) + t478 / 0.2e1 - t477 / 0.2e1 + t99 * mrSges(5,1) - t208 * mrSges(4,3) + t428 * t306 + t545) * t336 + (-t354 / 0.2e1 - t353 / 0.2e1 + t336 * (Ifges(3,4) * t336 - Ifges(3,2) * t332) / 0.2e1 + t332 * (Ifges(5,5) * t460 + Ifges(5,6) * t453 + t481) / 0.2e1 - t551) * t443) * qJD(2) + (g(1) * t337 + g(2) * t333) * (-m(5) * pkin(3) - mrSges(4,1) - mrSges(3,3) + t386)) * t325 + m(3) * (t159 * t263 + t160 * t262 + t243 * t246 + t244 * t245) + (Ifges(7,1) * t135 + Ifges(7,4) * t134) * t535 + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t518 + (-m(3) * t403 + t333 * mrSges(2,1) + mrSges(2,2) * t337 + t402 * t178 - t567 * t180 + t543 * t258 + t548 * t257 + t575 * t387 - t574 * (t257 * t321 - t258 * t328 - t316 * t456 + t387)) * g(1) + (t1 * t134 - t13 * t63 - t135 * t2 + t14 * t64) * mrSges(7,3) + m(6) * (t124 * t151 + t150 * t80 + t17 * t34 + t18 * t33 + t49 * t8 + t493 * t7) + t493 * t61 + t558 * t246 + (-mrSges(6,3) * t33 + t579 + t585) * t114 + Ifges(2,3) * qJDD(1) + t262 * (mrSges(3,1) * t304 - mrSges(3,3) * t248) + t263 * (-mrSges(3,2) * t304 - mrSges(3,3) * t247) + t253 * t94 / 0.2e1 + t107 * (-mrSges(5,1) * t253 + mrSges(5,2) * t254) + t225 * (-mrSges(4,2) * t247 - mrSges(4,3) * t248) + t245 * t238 + t209 * t239 + t206 * t241 + t224 * t191 + t227 * t192 + t197 * t112 + t183 * t147 + t105 * t171 + t104 * t172 + t150 * t30 + t151 * t82 + t6 * (-mrSges(7,1) * t134 + mrSges(7,2) * t135) + t119 * t131 + t118 * t132 + t17 * t120 + t18 * t121 + t4 * t74 + t3 * t73 + t31 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t49 * t60 + t16 * t62 + t46 * t12 + t24 * t20 + t25 * t21 + t64 * t531 + t135 * t537 + t134 * t538 + t63 * t529 + t254 * t523 + (-m(3) * t444 - mrSges(2,1) * t337 + t333 * mrSges(2,2) - m(7) * (pkin(5) * t177 + t362) - t130 * mrSges(7,1) - t129 * mrSges(7,2) - m(6) * t362 - t177 * mrSges(6,1) + t402 * t176 + t364 * t260 - t548 * t259 - t575 * t418) * g(2) + (Ifges(5,5) * t254 + Ifges(5,6) * t253) * t508 + (Ifges(5,1) * t254 + Ifges(5,4) * t253) * t510 + (Ifges(5,4) * t254 + Ifges(5,2) * t253) * t511; t564 * t62 + t565 * t121 + (-t124 * t161 - t186 * t8 + t187 * t7 + t314 * t80 + t33 * t565 + t34 * t566) * m(6) + t566 * t120 + (-t376 * t528 - t377 * t534 - t379 * t535 + t44 * t405 - t578) * t265 - t382 * t500 + (Ifges(6,4) * t513 + t540 + t552) * t195 - t374 * mrSges(5,3) + (-pkin(2) * t146 - qJ(3) * t128 - t208 * t242) * m(4) + (-t99 * (mrSges(5,1) * t336 - mrSges(5,3) * t460) - t100 * (-mrSges(5,2) * t336 + mrSges(5,3) * t453) - t208 * (-mrSges(4,2) * t332 - mrSges(4,3) * t336)) * t443 + t541 * mrSges(6,3) + (-t437 - t123) * t171 + (-t436 - t122) * t172 + t587 * t256 + t588 * t266 + (m(5) * t163 + m(6) * t124 + t549 - t584) * qJD(3) + (t238 - t397 + t549) * t243 + (Ifges(7,4) * t357 + Ifges(7,2) * t358) * t520 + (Ifges(7,4) * t156 + Ifges(7,2) * t155) * t521 - t203 * t400 + (t421 + t579) * t255 - t193 * t399 + (-t191 + t112) * qJ(3) + (Ifges(6,5) * t513 + Ifges(6,6) * t515 + Ifges(6,3) * t507 - t545) * t414 - t356 * t415 + (-t78 / 0.2e1 + Ifges(6,4) * t515 + t544) * t196 - (-t210 * (Ifges(5,5) * t336 + t332 * t380) + t347 * (Ifges(5,6) * t336 + t332 * t378) + (t455 + t462 + t201) * t332 + (-Ifges(3,2) * t415 + t125 + t200 + t291 + t76) * t336 + (t332 * t571 + t336 * t572) * t306) * t443 / 0.2e1 - t555 * t495 + (qJ(3) * t107 - t374 * t495 + (-t100 * t324 - t326 * t99) * qJD(4) - t100 * t123 - t122 * t99 - t163 * t204) * m(5) - t11 * t463 / 0.2e1 + (t353 + t354) * t440 / 0.2e1 + (t332 * t199 + t336 * t202) * t443 / 0.2e1 + (qJD(6) * t45 + t10) * t464 / 0.2e1 + t282 + t283 + t340 + (t261 - t575 * t445 + t574 * (pkin(4) * t423 + t445) + (t381 + (-t328 * t574 - t382 + t592) * t336 + t594 * t332) * t325) * g(3) + (t1 * t464 + t2 * t463) * mrSges(7,3) + t230 - t231 - t228 + t229 + (mrSges(6,1) * t556 + mrSges(6,2) * t557) * t124 + (-m(4) * t193 + t398 - t558) * t244 + (mrSges(7,1) * t556 + mrSges(7,3) * t559) * t13 + (mrSges(7,1) * t560 - mrSges(7,2) * t559) * t31 + (-mrSges(7,2) * t556 - mrSges(7,3) * t560) * t14 + (Ifges(7,5) * t357 + Ifges(7,6) * t358) * t516 + (Ifges(7,5) * t156 + Ifges(7,6) * t155) * t517 + t107 * t385 + (Ifges(7,1) * t156 + Ifges(7,4) * t155) * t519 + (Ifges(7,1) * t357 + Ifges(7,4) * t358) * t518 - t324 * t94 / 0.2e1 + t314 * t30 - t242 * t241 - t204 * t147 + t187 * t61 - pkin(2) * t192 - t161 * t82 + t109 * t20 + t110 * t21 + t156 * t530 + t155 * t532 + t449 * t529 + t326 * t523 + t568 * t74 + t569 * t73 + (t1 * t110 + t109 * t2 + t13 * t568 + t14 * t569 + t186 * t6 + t31 * t564) * m(7) + t570 * t186 + ((t481 + (Ifges(5,5) * t324 + Ifges(5,6) * t326) * t332) * t586 + t551) * t440 + (t574 * (t259 * t328 + t260 * t321 - t251) + t575 * t251 + t542 * t260 + t543 * t259) * g(1) + (t574 * (t257 * t328 + t258 * t321 - t249) + t575 * t249 + t542 * t258 + t543 * t257) * g(2) + (Ifges(5,5) * t326 - Ifges(5,6) * t324) * t508 + (Ifges(5,1) * t326 - t489) * t510 + (-Ifges(5,2) * t324 + t488) * t511; (t171 * t326 - t172 * t324 + t241) * t415 + (t61 + (-t330 * t73 - t334 * t74) * qJD(6) + t554) * t266 + (t239 + t584) * t306 + t570 * t265 + t192 - t363 * t256 + t195 * t120 - t153 * t74 - t154 * t73 + t563 * t557 - (-g(1) * t259 - g(2) * t257 + g(3) * t457) * t595 + (-t13 * t153 - t14 * t154 + t256 * t375 + t266 * t339 - t31 * t557 + t500) * m(7) + (-t124 * t306 - t541) * m(6) + (-t163 * t306 + t370 * t415 + t374) * m(5) + (t203 * t306 + t208 * t415 + t146) * m(4) + t555; -t347 * t171 - t210 * t172 + t334 * t20 + t330 * t21 + t563 * t581 + t373 * qJD(6) + t363 * t582 + t112 + t30 + (t1 * t330 + t137 * t375 + t2 * t334 - t581 * t31 + t344) * m(7) + (t33 * t581 - t34 * t582 + t344 + t80) * m(6) + (-t100 * t347 - t210 * t99 + t107 + t344) * m(5); (-m(7) * t31 + t491 + t563) * t34 + (t176 * t567 + t177 * t402) * g(1) + (t402 * t219 - t567 * (-t318 * t327 - t319 * t457)) * g(3) + (-t178 * t567 - t180 * t402) * g(2) + ((-t433 + t467) * t14 + (-t432 + t466) * t13 + t553) * mrSges(7,3) + (m(7) * t339 - t432 * t74 - t433 * t73 + t554) * pkin(10) + t547 + (t376 * t517 + t377 * t521 + t379 * t519 - t360 + t544 - t585) * t582 + (-t540 + t580) * t581 + (-t120 + t492) * t33 + (t360 + t421) * qJD(6) + (t116 * t377 + t117 * t379 + t137 * t376) * qJD(6) / 0.2e1 + t45 * t405 + (t136 + t78) * t515 + t6 * t383 + (-pkin(5) * t6 - t13 * t22 - t14 * t23) * m(7) - t22 * t74 - t23 * t73 + t467 * t531 + (Ifges(7,2) * t334 + t485) * t534 + (Ifges(7,1) * t330 + t484) * t535 + t330 * t537 + t334 * t538 + (Ifges(7,5) * t330 + Ifges(7,6) * t334) * t528 + t466 * t530 - pkin(5) * t12 + t77 * t512 + (-t487 + t43) * t513; -t31 * (mrSges(7,1) * t117 + mrSges(7,2) * t116) + (Ifges(7,1) * t116 - t486) * t519 + t44 * t518 + (Ifges(7,5) * t116 - Ifges(7,6) * t117) * t517 - t13 * t73 + t14 * t74 - g(1) * (mrSges(7,1) * t129 - mrSges(7,2) * t130) - g(2) * ((t180 * t330 + t258 * t334) * mrSges(7,1) + (t180 * t334 - t258 * t330) * mrSges(7,2)) - g(3) * ((-t219 * t330 + t286) * mrSges(7,1) + (-t219 * t334 - t422) * mrSges(7,2)) + (t116 * t13 + t117 * t14) * mrSges(7,3) + t9 + (-Ifges(7,2) * t117 + t113 + t45) * t521 + t546;];
tau  = t26;
