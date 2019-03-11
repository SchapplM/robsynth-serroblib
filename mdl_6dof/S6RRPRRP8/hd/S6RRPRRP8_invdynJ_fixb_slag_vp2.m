% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:00
% EndTime: 2019-03-09 12:23:01
% DurationCPUTime: 36.29s
% Computational Cost: add. (16827->835), mult. (37701->1067), div. (0->0), fcn. (27936->14), ass. (0->381)
t332 = sin(pkin(10));
t333 = cos(pkin(10));
t336 = sin(qJ(4));
t339 = cos(qJ(4));
t278 = t332 * t339 + t333 * t336;
t340 = cos(qJ(2));
t361 = t278 * t340;
t229 = qJD(1) * t361;
t253 = t278 * qJD(4);
t606 = t229 - t253;
t374 = t332 * t336 - t333 * t339;
t360 = t374 * t340;
t230 = qJD(1) * t360;
t252 = t374 * qJD(4);
t557 = -t230 + t252;
t604 = mrSges(6,1) + mrSges(7,1);
t603 = mrSges(6,2) - mrSges(7,3);
t337 = sin(qJ(2));
t377 = pkin(2) * t337 - qJ(3) * t340;
t281 = t377 * qJD(1);
t432 = qJD(1) * t337;
t410 = t332 * t432;
t205 = pkin(7) * t410 + t333 * t281;
t444 = t333 * t340;
t372 = pkin(3) * t337 - pkin(8) * t444;
t170 = qJD(1) * t372 + t205;
t254 = t332 * t281;
t445 = t333 * t337;
t446 = t332 * t340;
t362 = -pkin(7) * t445 - pkin(8) * t446;
t192 = qJD(1) * t362 + t254;
t473 = pkin(8) + qJ(3);
t290 = t473 * t332;
t291 = t473 * t333;
t424 = qJD(4) * t339;
t427 = qJD(3) * t333;
t428 = qJD(3) * t332;
t561 = -t290 * t424 + (-t192 + t427) * t339 + (-qJD(4) * t291 - t170 - t428) * t336;
t203 = -t336 * t290 + t339 * t291;
t560 = -t278 * qJD(3) - qJD(4) * t203 - t339 * t170 + t192 * t336;
t609 = pkin(9) * t606 + t561;
t608 = -pkin(4) * t432 + pkin(9) * t557 + t560;
t575 = -Ifges(6,4) + Ifges(7,5);
t607 = t575 + Ifges(7,5);
t271 = qJD(2) * t333 - t410;
t408 = t333 * t432;
t272 = qJD(2) * t332 + t408;
t188 = t271 * t339 - t272 * t336;
t422 = qJD(1) * qJD(2);
t284 = qJDD(1) * t337 + t340 * t422;
t231 = qJDD(2) * t333 - t284 * t332;
t232 = qJDD(2) * t332 + t284 * t333;
t114 = qJD(4) * t188 + t231 * t336 + t232 * t339;
t189 = t271 * t336 + t272 * t339;
t115 = -qJD(4) * t189 + t231 * t339 - t232 * t336;
t335 = sin(qJ(5));
t488 = cos(qJ(5));
t125 = -t488 * t188 + t335 * t189;
t42 = -qJD(5) * t125 + t114 * t488 + t335 * t115;
t525 = t42 / 0.2e1;
t589 = t188 * t335 + t189 * t488;
t43 = qJD(5) * t589 + t335 * t114 - t115 * t488;
t523 = t43 / 0.2e1;
t605 = m(6) + m(7);
t324 = t340 * qJDD(1);
t403 = t337 * t422;
t283 = -t324 + t403;
t276 = qJDD(4) + t283;
t263 = qJDD(5) + t276;
t496 = t263 / 0.2e1;
t602 = -mrSges(6,3) - mrSges(7,2);
t576 = Ifges(6,1) + Ifges(7,1);
t574 = Ifges(7,4) + Ifges(6,5);
t573 = -Ifges(6,6) + Ifges(7,6);
t572 = Ifges(6,3) + Ifges(7,2);
t319 = pkin(7) * t432;
t286 = -qJD(2) * pkin(2) + qJD(3) + t319;
t204 = -t271 * pkin(3) + t286;
t138 = -pkin(4) * t188 + t204;
t431 = qJD(1) * t340;
t310 = qJD(4) - t431;
t304 = qJD(5) + t310;
t396 = -qJ(3) * t337 - pkin(1);
t289 = -pkin(2) * t340 + t396;
t259 = t289 * qJD(1);
t320 = pkin(7) * t431;
t293 = qJD(2) * qJ(3) + t320;
t194 = t333 * t259 - t293 * t332;
t144 = -pkin(3) * t431 - pkin(8) * t272 + t194;
t195 = t332 * t259 + t333 * t293;
t147 = pkin(8) * t271 + t195;
t88 = t336 * t144 + t339 * t147;
t79 = pkin(9) * t188 + t88;
t460 = t335 * t79;
t87 = t339 * t144 - t147 * t336;
t78 = -pkin(9) * t189 + t87;
t76 = pkin(4) * t310 + t78;
t27 = t488 * t76 - t460;
t597 = qJD(6) - t27;
t23 = -t304 * pkin(5) + t597;
t50 = t125 * pkin(5) - qJ(6) * t589 + t138;
t122 = Ifges(6,4) * t125;
t461 = Ifges(7,5) * t125;
t566 = t304 * t574 + t576 * t589 - t122 + t461;
t601 = mrSges(6,2) * t138 - mrSges(7,3) * t50 + mrSges(7,2) * t23 - mrSges(6,3) * t27 + t566 / 0.2e1;
t414 = t488 * t79;
t28 = t335 * t76 + t414;
t24 = t304 * qJ(6) + t28;
t121 = Ifges(7,5) * t589;
t64 = Ifges(7,6) * t304 + Ifges(7,3) * t125 + t121;
t462 = Ifges(6,4) * t589;
t67 = -Ifges(6,2) * t125 + Ifges(6,6) * t304 + t462;
t600 = -mrSges(7,2) * t24 - mrSges(6,3) * t28 + mrSges(6,1) * t138 + mrSges(7,1) * t50 + t64 / 0.2e1 - t67 / 0.2e1;
t331 = pkin(10) + qJ(4);
t325 = qJ(5) + t331;
t311 = sin(t325);
t312 = cos(t325);
t314 = pkin(3) * t333 + pkin(2);
t322 = sin(t331);
t323 = cos(t331);
t599 = m(5) * t314 + t323 * mrSges(5,1) - t322 * mrSges(5,2) - t603 * t311 + t312 * t604;
t463 = Ifges(5,4) * t189;
t117 = Ifges(5,2) * t188 + Ifges(5,6) * t310 + t463;
t513 = t117 / 0.2e1;
t202 = -t339 * t290 - t291 * t336;
t162 = -pkin(9) * t278 + t202;
t163 = -pkin(9) * t374 + t203;
t102 = t335 * t162 + t163 * t488;
t568 = -qJD(5) * t102 - t335 * t609 + t608 * t488;
t367 = t162 * t488 - t335 * t163;
t567 = qJD(5) * t367 + t608 * t335 + t488 * t609;
t598 = -m(7) * qJ(6) - mrSges(7,3);
t468 = mrSges(6,3) * t589;
t109 = mrSges(6,1) * t304 - t468;
t110 = -mrSges(7,1) * t304 + mrSges(7,2) * t589;
t596 = t110 - t109;
t409 = t332 * t431;
t264 = pkin(3) * t409 + t320;
t559 = -pkin(4) * t606 - t264;
t338 = sin(qJ(1));
t341 = cos(qJ(1));
t440 = t340 * t341;
t239 = -t322 * t440 + t338 * t323;
t595 = g(1) * t341 + g(2) * t338;
t426 = qJD(3) * t337;
t456 = qJDD(1) * pkin(1);
t180 = pkin(2) * t283 - qJ(3) * t284 - qJD(1) * t426 - t456;
t317 = pkin(7) * t324;
t236 = qJDD(2) * qJ(3) + t317 + (qJD(3) - t319) * qJD(2);
t136 = t332 * t180 + t333 * t236;
t106 = pkin(8) * t231 + t136;
t135 = t333 * t180 - t236 * t332;
t95 = pkin(3) * t283 - pkin(8) * t232 + t135;
t26 = -qJD(4) * t88 - t106 * t336 + t339 * t95;
t19 = pkin(4) * t276 - pkin(9) * t114 + t26;
t425 = qJD(4) * t336;
t25 = t339 * t106 + t144 * t424 - t147 * t425 + t336 * t95;
t21 = pkin(9) * t115 + t25;
t6 = -qJD(5) * t28 + t19 * t488 - t335 * t21;
t3 = -t263 * pkin(5) + qJDD(6) - t6;
t524 = -t43 / 0.2e1;
t270 = t284 * pkin(7);
t247 = -qJDD(2) * pkin(2) + qJDD(3) + t270;
t167 = -pkin(3) * t231 + t247;
t83 = -pkin(4) * t115 + t167;
t7 = pkin(5) * t43 - qJ(6) * t42 - qJD(6) * t589 + t83;
t593 = mrSges(6,2) * t83 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t7 + Ifges(6,4) * t524 + 0.2e1 * t496 * t574 + t523 * t607 + 0.2e1 * t525 * t576;
t492 = t304 / 0.2e1;
t506 = t589 / 0.2e1;
t509 = t125 / 0.2e1;
t510 = -t125 / 0.2e1;
t592 = -Ifges(6,2) * t510 + Ifges(7,3) * t509 + t492 * t573 + t506 * t575 + t600;
t493 = -t304 / 0.2e1;
t507 = -t589 / 0.2e1;
t591 = -Ifges(6,2) * t509 + Ifges(7,3) * t510 + t493 * t573 + t507 * t575 - t600;
t73 = pkin(5) * t589 + qJ(6) * t125;
t587 = Ifges(6,4) * t509 + Ifges(7,5) * t510 + t493 * t574 + t507 * t576 - t601;
t586 = Ifges(6,4) * t510 + Ifges(7,5) * t509 + t492 * t574 + t506 * t576 + t601;
t585 = -m(4) - m(5);
t515 = t114 / 0.2e1;
t514 = t115 / 0.2e1;
t583 = -t188 / 0.2e1;
t582 = t188 / 0.2e1;
t498 = t231 / 0.2e1;
t497 = t232 / 0.2e1;
t495 = t276 / 0.2e1;
t581 = -t283 / 0.2e1;
t494 = t283 / 0.2e1;
t580 = t284 / 0.2e1;
t577 = qJD(2) / 0.2e1;
t570 = -qJ(6) * t432 + t567;
t569 = pkin(5) * t432 - t568;
t365 = -t335 * t278 - t374 * t488;
t129 = qJD(5) * t365 - t252 * t488 - t335 * t253;
t191 = t278 * t488 - t335 * t374;
t130 = qJD(5) * t191 - t335 * t252 + t253 * t488;
t145 = t229 * t488 - t230 * t335;
t146 = -t335 * t229 - t230 * t488;
t565 = -qJD(6) * t191 + t559 + (-t129 + t146) * qJ(6) + (t130 - t145) * pkin(5);
t564 = t188 * Ifges(5,6);
t386 = -mrSges(4,1) * t333 + mrSges(4,2) * t332;
t359 = m(4) * pkin(2) - t386;
t562 = t340 * t359;
t268 = t333 * t289;
t193 = -pkin(8) * t445 + t268 + (-pkin(7) * t332 - pkin(3)) * t340;
t225 = pkin(7) * t444 + t332 * t289;
t447 = t332 * t337;
t201 = -pkin(8) * t447 + t225;
t133 = t336 * t193 + t339 * t201;
t556 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t271 + mrSges(4,2) * t272 + mrSges(3,3) * t432;
t464 = Ifges(4,4) * t333;
t380 = -Ifges(4,2) * t332 + t464;
t465 = Ifges(4,4) * t332;
t382 = Ifges(4,1) * t333 - t465;
t555 = t271 * (Ifges(4,6) * t337 + t340 * t380) + t272 * (Ifges(4,5) * t337 + t340 * t382);
t554 = t263 * t572 + t42 * t574 + t43 * t573;
t222 = t311 * t440 - t338 * t312;
t223 = t338 * t311 + t312 * t440;
t553 = t222 * t604 + t223 * t603;
t441 = t338 * t340;
t220 = t311 * t441 + t312 * t341;
t221 = -t311 * t341 + t312 * t441;
t552 = t220 * t604 + t221 * t603;
t269 = -pkin(7) * t403 + t317;
t551 = t269 * t340 + t270 * t337;
t550 = -t135 * t332 + t136 * t333;
t549 = t602 * t337;
t548 = -m(3) + t585;
t388 = mrSges(3,1) * t340 - mrSges(3,2) * t337;
t547 = t337 * mrSges(5,3) + mrSges(2,1) + t388;
t385 = mrSges(4,1) * t332 + t333 * mrSges(4,2);
t546 = -t286 * t340 * t385 - t195 * (-mrSges(4,2) * t337 - mrSges(4,3) * t446) - t194 * (mrSges(4,1) * t337 - mrSges(4,3) * t444);
t318 = Ifges(3,4) * t431;
t545 = t333 * (t272 * Ifges(4,1) + t271 * Ifges(4,4) - Ifges(4,5) * t431) + Ifges(3,1) * t432 + Ifges(3,5) * qJD(2) + t318;
t544 = t272 * Ifges(4,5) + Ifges(5,5) * t189 + t271 * Ifges(4,6) - Ifges(4,3) * t431 + Ifges(5,3) * t310 + t125 * t573 + t304 * t572 + t574 * t589 + t564;
t244 = t278 * t337;
t104 = -pkin(9) * t244 + t133;
t132 = t339 * t193 - t201 * t336;
t245 = t374 * t337;
t98 = -pkin(4) * t340 + pkin(9) * t245 + t132;
t472 = t488 * t104 + t335 * t98;
t164 = -qJD(2) * t360 - t253 * t337;
t430 = qJD(2) * t337;
t248 = qJD(2) * t377 - t426;
t416 = pkin(7) * t430;
t199 = t333 * t248 + t332 * t416;
t159 = qJD(2) * t372 + t199;
t234 = t332 * t248;
t171 = qJD(2) * t362 + t234;
t71 = -qJD(4) * t133 + t339 * t159 - t171 * t336;
t52 = pkin(4) * t430 - pkin(9) * t164 + t71;
t165 = -qJD(2) * t361 + t252 * t337;
t70 = t336 * t159 + t339 * t171 + t193 * t424 - t201 * t425;
t55 = pkin(9) * t165 + t70;
t11 = -qJD(5) * t472 - t335 * t55 + t488 * t52;
t543 = m(7) * pkin(5) + t604;
t542 = -mrSges(6,2) - t598;
t541 = -t26 * mrSges(5,1) + t25 * mrSges(5,2);
t483 = pkin(3) * t332;
t537 = -m(5) * t483 + mrSges(2,2) - mrSges(3,3) - t385;
t404 = qJD(5) * t488;
t423 = qJD(5) * t335;
t5 = t335 * t19 + t488 * t21 + t76 * t404 - t423 * t79;
t2 = qJ(6) * t263 + qJD(6) * t304 + t5;
t536 = -t6 * mrSges(6,1) + t3 * mrSges(7,1) + t5 * mrSges(6,2) - t2 * mrSges(7,3);
t467 = Ifges(3,4) * t337;
t381 = t340 * Ifges(3,2) + t467;
t533 = t24 * mrSges(7,3) + t27 * mrSges(6,1) + t87 * mrSges(5,1) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t381 / 0.2e1 + t564 / 0.2e1 - t23 * mrSges(7,1) - t28 * mrSges(6,2) - t88 * mrSges(5,2);
t529 = mrSges(6,1) * t83 + mrSges(7,1) * t7 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t523 - t42 * Ifges(6,4) / 0.2e1 - t263 * Ifges(6,6) / 0.2e1 + t607 * t525 + (t573 + Ifges(7,6)) * t496 + (-t524 + t523) * Ifges(6,2);
t522 = Ifges(5,4) * t515 + Ifges(5,2) * t514 + Ifges(5,6) * t495;
t521 = Ifges(5,1) * t515 + Ifges(5,4) * t514 + Ifges(5,5) * t495;
t186 = Ifges(5,4) * t188;
t118 = t189 * Ifges(5,1) + t310 * Ifges(5,5) + t186;
t512 = -t118 / 0.2e1;
t511 = t118 / 0.2e1;
t504 = Ifges(4,1) * t497 + Ifges(4,4) * t498 + Ifges(4,5) * t494;
t501 = -t189 / 0.2e1;
t500 = t189 / 0.2e1;
t491 = -t310 / 0.2e1;
t490 = t310 / 0.2e1;
t482 = pkin(4) * t189;
t480 = pkin(4) * t322;
t479 = pkin(4) * t335;
t476 = g(3) * t337;
t326 = t337 * pkin(7);
t471 = mrSges(6,2) * t312;
t470 = mrSges(5,3) * t189;
t469 = mrSges(6,3) * t125;
t466 = Ifges(3,4) * t340;
t452 = t247 * t337;
t330 = -pkin(9) - t473;
t448 = t330 * t337;
t443 = t337 * t341;
t280 = pkin(4) * t323 + t314;
t256 = t340 * t280;
t107 = -mrSges(7,2) * t125 + mrSges(7,3) * t304;
t108 = -mrSges(6,2) * t304 - t469;
t439 = -t107 - t108;
t429 = qJD(2) * t340;
t321 = pkin(7) * t429;
t407 = t332 * t429;
t265 = pkin(3) * t407 + t321;
t282 = pkin(3) * t447 + t326;
t433 = t341 * pkin(1) + t338 * pkin(7);
t418 = t488 * pkin(4);
t417 = pkin(4) * t423;
t412 = Ifges(5,5) * t114 + Ifges(5,6) * t115 + Ifges(5,3) * t276;
t406 = m(4) * qJ(3) + mrSges(4,3);
t405 = m(5) * t473 + mrSges(5,3);
t17 = t43 * mrSges(6,1) + t42 * mrSges(6,2);
t16 = t43 * mrSges(7,1) - t42 * mrSges(7,3);
t397 = t598 * t312 * t337;
t32 = -t263 * mrSges(7,1) + t42 * mrSges(7,2);
t394 = -t422 / 0.2e1;
t393 = t422 / 0.2e1;
t150 = -t231 * mrSges(4,1) + t232 * mrSges(4,2);
t63 = -t115 * mrSges(5,1) + t114 * mrSges(5,2);
t392 = -t220 * pkin(5) + t221 * qJ(6);
t391 = -t222 * pkin(5) + t223 * qJ(6);
t389 = pkin(4) * t404;
t141 = -pkin(4) * t165 + t265;
t197 = pkin(4) * t244 + t282;
t387 = mrSges(3,1) * t337 + mrSges(3,2) * t340;
t379 = Ifges(3,5) * t340 - Ifges(3,6) * t337;
t378 = Ifges(4,5) * t333 - Ifges(4,6) * t332;
t376 = pkin(5) * t312 + qJ(6) * t311;
t375 = t314 * t340 + t337 * t473;
t373 = t239 * pkin(4);
t233 = pkin(4) * t374 - t314;
t369 = pkin(1) * t387;
t61 = -t335 * t104 + t488 * t98;
t237 = t322 * t441 + t323 * t341;
t366 = -t244 * t488 + t335 * t245;
t158 = -t335 * t244 - t245 * t488;
t10 = -t104 * t423 + t335 * t52 + t98 * t404 + t488 * t55;
t364 = t337 * (Ifges(3,1) * t340 - t467);
t363 = t188 * mrSges(5,3);
t355 = t237 * pkin(4);
t350 = t340 * (Ifges(4,3) * t337 + t340 * t378);
t348 = -t536 + t554;
t347 = t337 * t406 + t562;
t328 = t341 * pkin(7);
t315 = -t418 - pkin(5);
t313 = qJ(6) + t479;
t305 = t389 + qJD(6);
t294 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t431;
t285 = t480 + t483;
t240 = t338 * t322 + t323 * t440;
t238 = t322 * t341 - t323 * t441;
t228 = -mrSges(4,1) * t431 - mrSges(4,3) * t272;
t227 = mrSges(4,2) * t431 + mrSges(4,3) * t271;
t224 = -pkin(7) * t446 + t268;
t206 = -pkin(7) * t408 + t254;
t200 = -t333 * t416 + t234;
t182 = t272 * Ifges(4,4) + t271 * Ifges(4,2) - Ifges(4,6) * t431;
t178 = mrSges(4,1) * t283 - mrSges(4,3) * t232;
t177 = -mrSges(4,2) * t283 + mrSges(4,3) * t231;
t155 = mrSges(5,1) * t310 - t470;
t154 = -t310 * mrSges(5,2) + t363;
t139 = t232 * Ifges(4,4) + t231 * Ifges(4,2) + t283 * Ifges(4,6);
t131 = -mrSges(5,1) * t188 + t189 * mrSges(5,2);
t100 = -pkin(5) * t365 - qJ(6) * t191 + t233;
t93 = -mrSges(5,2) * t276 + mrSges(5,3) * t115;
t92 = mrSges(5,1) * t276 - mrSges(5,3) * t114;
t84 = -pkin(5) * t366 - qJ(6) * t158 + t197;
t82 = qJD(5) * t158 + t335 * t164 - t165 * t488;
t81 = qJD(5) * t366 + t164 * t488 + t335 * t165;
t75 = mrSges(6,1) * t125 + mrSges(6,2) * t589;
t74 = mrSges(7,1) * t125 - mrSges(7,3) * t589;
t60 = t482 + t73;
t57 = t340 * pkin(5) - t61;
t56 = -qJ(6) * t340 + t472;
t34 = -mrSges(7,2) * t43 + mrSges(7,3) * t263;
t33 = -mrSges(6,2) * t263 - mrSges(6,3) * t43;
t31 = mrSges(6,1) * t263 - mrSges(6,3) * t42;
t30 = t488 * t78 - t460;
t29 = t335 * t78 + t414;
t22 = pkin(5) * t82 - qJ(6) * t81 - qJD(6) * t158 + t141;
t9 = -pkin(5) * t430 - t11;
t8 = qJ(6) * t430 - qJD(6) * t340 + t10;
t1 = [(-mrSges(3,1) * t326 + Ifges(3,5) * t337 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * t340) * qJDD(2) + (-Ifges(5,5) * t245 - Ifges(5,6) * t244) * t495 + t167 * (mrSges(5,1) * t244 - mrSges(5,2) * t245) + (-Ifges(5,1) * t245 - Ifges(5,4) * t244) * t515 + (-Ifges(5,4) * t245 - Ifges(5,2) * t244) * t514 + (-t164 * t87 + t165 * t88 - t244 * t25 + t245 * t26) * mrSges(5,3) + t592 * t82 + t593 * t158 + t586 * t81 + (Ifges(5,1) * t164 + Ifges(5,4) * t165) * t500 + (-t240 * mrSges(5,1) - t239 * mrSges(5,2) + t602 * t443 - t605 * (t280 * t440 + t338 * t285 - t330 * t443 + t433) - t543 * t223 - t542 * t222 + t548 * t433 + t537 * t338 + (-m(5) * t375 - t347 - t547) * t341) * g(2) + (-t238 * mrSges(5,1) - t237 * mrSges(5,2) - t605 * (t341 * t285 + t338 * t448 + t328) + t543 * t221 + t542 * t220 + t548 * t328 + t537 * t341 + (-m(4) * t396 + t337 * mrSges(4,3) + t562 + m(3) * pkin(1) - m(5) * (-pkin(1) - t375) - t605 * (-pkin(1) - t256) + t547 - t549) * t338) * g(1) + m(6) * (t10 * t28 + t11 * t27 + t138 * t141 + t197 * t83 + t472 * t5 + t6 * t61) + t472 * t33 + (Ifges(5,5) * t164 + Ifges(5,6) * t165) * t490 + t197 * t17 + t204 * (-mrSges(5,1) * t165 + mrSges(5,2) * t164) + (t379 * t577 - t546) * qJD(2) + (-pkin(7) * t283 * t340 + t284 * t326 + t551) * mrSges(3,3) + m(7) * (t2 * t56 + t22 * t50 + t23 * t9 + t24 * t8 + t3 * t57 + t7 * t84) + m(5) * (t132 * t26 + t133 * t25 + t167 * t282 + t204 * t265 + t70 * t88 + t71 * t87) + t224 * t178 + t225 * t177 + t200 * t227 + t199 * t228 + (Ifges(5,5) * t500 + Ifges(6,6) * t510 + Ifges(7,6) * t509 + Ifges(5,3) * t490 + t492 * t572 + t506 * t574 + t533 + t544 / 0.2e1) * t430 + t545 * t429 / 0.2e1 - t294 * t416 - t369 * t422 + m(4) * (t135 * t224 + t136 * t225 + t194 * t199 + t195 * t200 + (t286 * t429 + t452) * pkin(7)) + t466 * t580 + t381 * t581 + (Ifges(5,4) * t164 + Ifges(5,2) * t165) * t582 + (-t135 * t445 - t136 * t447) * mrSges(4,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t551) - (Ifges(4,5) * t232 + Ifges(4,6) * t231 + Ifges(4,3) * t283 + t412 + t554) * t340 / 0.2e1 + t556 * t321 - t139 * t447 / 0.2e1 - t182 * t407 / 0.2e1 - t529 * t366 + t388 * t456 + t70 * t154 + t71 * t155 + t141 * t75 + t132 * t92 + t133 * t93 + Ifges(2,3) * qJDD(1) + t9 * t110 + t8 * t107 + t10 * t108 + t11 * t109 + t84 * t16 + t22 * t74 + t61 * t31 + t56 * t34 + t57 * t32 - t245 * t521 - t244 * t522 + t164 * t511 + t165 * t513 + t445 * t504 + t265 * t131 + (-Ifges(4,5) * t497 - Ifges(4,6) * t498 + (-Ifges(3,2) * t337 + t466) * t393 - Ifges(7,6) * t523 - Ifges(6,6) * t524 - Ifges(5,6) * t514 - Ifges(5,5) * t515 - Ifges(4,3) * t494 - Ifges(5,3) * t495 + t536 - t574 * t525 - t572 * t496 + Ifges(3,4) * t580 + Ifges(3,2) * t581 + t136 * mrSges(4,2) - t135 * mrSges(4,1) + t541) * t340 + (Ifges(3,1) * t284 + Ifges(3,4) * t581 + t378 * t494 + t380 * t498 + t382 * t497) * t337 + t282 * t63 - pkin(1) * (mrSges(3,1) * t283 + mrSges(3,2) * t284) + t364 * t393 + t350 * t394 + t385 * t452 + t150 * t326 + t555 * t577; t591 * t145 + t592 * t130 + t593 * t191 - (t32 - t31) * t367 + (t102 * t5 + t138 * t559 + t233 * t83 + t27 * t568 + t28 * t567 + t367 * t6) * m(6) + (t100 * t7 + t102 * t2 + t23 * t569 + t24 * t570 - t3 * t367 + t50 * t565) * m(7) + t586 * t129 + (-t337 * t405 - t388 - t347 - t605 * (t256 - t448) + (-m(7) * t376 - t599) * t340 + t549) * g(3) + t568 * t109 + t569 * t110 + t570 * t107 + (Ifges(5,5) * t501 + Ifges(6,6) * t509 + Ifges(7,6) * t510 + Ifges(5,3) * t491 + t493 * t572 + t507 * t574 - t533) * t432 + t565 * t74 + t567 * t108 + t247 * t386 + t202 * t92 + t203 * t93 + (-t205 - t428) * t228 + (-Ifges(5,1) * t252 - Ifges(5,4) * t253) * t500 + (-Ifges(5,5) * t252 - Ifges(5,6) * t253) * t490 + (-t194 * t205 - t195 * t206 - t286 * t320 - pkin(2) * t247 + (-t194 * t332 + t195 * t333) * qJD(3) + t550 * qJ(3)) * m(4) + t550 * mrSges(4,3) + t559 * t75 + t560 * t155 + (-t167 * t314 + t202 * t26 + t203 * t25 - t204 * t264 + t560 * t87 + t561 * t88) * m(5) + t561 * t154 + t233 * t17 - t544 * t432 / 0.2e1 - (-Ifges(3,2) * t432 + t318 + t545) * t431 / 0.2e1 + (t177 * t333 - t178 * t332) * qJ(3) + (-t206 + t427) * t227 + (-Ifges(5,4) * t252 - Ifges(5,2) * t253) * t582 + (-mrSges(5,1) * t606 - mrSges(5,2) * t557) * t204 + (-t25 * t374 - t26 * t278 + t557 * t87 + t606 * t88) * mrSges(5,3) + t294 * t319 - t556 * t320 + t595 * (t387 + (t330 * t605 - t405 - t406 + t602) * t340 + (t359 + m(6) * t280 - m(7) * (-t280 - t376) + t599) * t337) + (t33 + t34) * t102 - t529 * t365 + (Ifges(5,4) * t278 - Ifges(5,2) * t374) * t514 + (Ifges(5,1) * t278 - Ifges(5,4) * t374) * t515 + (Ifges(5,5) * t278 - Ifges(5,6) * t374) * t495 + t167 * (mrSges(5,1) * t374 + mrSges(5,2) * t278) - t374 * t522 + (-Ifges(5,5) * t230 - Ifges(5,6) * t229) * t491 + (-Ifges(5,1) * t230 - Ifges(5,4) * t229) * t501 + (-Ifges(5,4) * t230 - Ifges(5,2) * t229) * t583 + t587 * t146 - pkin(2) * t150 + (t546 - t555 / 0.2e1 + (t369 + t350 / 0.2e1 - t364 / 0.2e1) * qJD(1)) * qJD(1) + t100 * t16 + t606 * t513 + t278 * t521 - t252 * t511 - t230 * t512 + t332 * t504 + (Ifges(4,5) * t332 + Ifges(4,6) * t333) * t494 + (Ifges(4,1) * t332 + t464) * t497 + (Ifges(4,2) * t333 + t465) * t498 + Ifges(3,3) * qJDD(2) - t264 * t131 - t269 * mrSges(3,2) - t270 * mrSges(3,1) + t182 * t409 / 0.2e1 - Ifges(3,6) * t283 + Ifges(3,5) * t284 - t314 * t63 + t379 * t394 + t333 * t139 / 0.2e1; -t596 * t589 - t439 * t125 - t188 * t154 + t189 * t155 - t271 * t227 + t272 * t228 + t150 + t16 + t17 + t63 + (t125 * t24 - t23 * t589 + t7) * m(7) + (t125 * t28 + t27 * t589 + t83) * m(6) + (-t188 * t88 + t189 * t87 + t167) * m(5) + (t194 * t272 - t195 * t271 + t247) * m(4) + (t340 * g(3) - t337 * t595) * (t605 - t585); t591 * t589 - t541 + t596 * (t417 - t29) + t412 + t348 + m(7) * (t2 * t313 + t23 * t417 + t24 * t305 + t3 * t315) + (-Ifges(5,2) * t189 + t186) * t583 + (-t30 + t305) * t107 + (t155 + t470) * t88 + (mrSges(5,1) * t237 - mrSges(5,2) * t238 - m(7) * (-t355 + t392) + m(6) * t355 + t552) * g(2) + (-mrSges(5,1) * t239 + mrSges(5,2) * t240 - m(7) * (t373 + t391) - m(6) * t373 + t553) * g(1) + ((t488 * t6 + t335 * t5 + (-t27 * t335 + t28 * t488) * qJD(5)) * pkin(4) + t480 * t476 - t138 * t482 + t27 * t29 - t28 * t30) * m(6) + (-t30 + t389) * t108 - t204 * (t189 * mrSges(5,1) + mrSges(5,2) * t188) + (mrSges(5,1) * t322 + mrSges(6,1) * t311 + mrSges(5,2) * t323 + t471) * t476 - t587 * t125 + (t363 - t154) * t87 - t60 * t74 + t188 * t512 + t117 * t500 + (Ifges(5,1) * t188 - t463) * t501 - g(3) * ((m(7) * (-pkin(5) * t311 - t480) - t311 * mrSges(7,1)) * t337 - t397) - t75 * t482 + t313 * t34 + t315 * t32 + t31 * t418 + t33 * t479 + (Ifges(5,5) * t188 - Ifges(5,6) * t189) * t491 - m(7) * (t23 * t29 + t24 * t30 + t50 * t60); (t125 * t23 + t24 * t589) * mrSges(7,2) + t348 + (-t596 + t468) * t28 + (t439 - t469) * t27 + t553 * g(1) + t552 * g(2) + ((t311 * t543 + t471) * t337 + t397) * g(3) - t138 * (mrSges(6,1) * t589 - mrSges(6,2) * t125) - t50 * (mrSges(7,1) * t589 + mrSges(7,3) * t125) + qJD(6) * t107 - t73 * t74 - pkin(5) * t32 + qJ(6) * t34 + (Ifges(7,3) * t589 - t461) * t510 + t67 * t506 + (-t125 * t574 + t573 * t589) * t493 + (-Ifges(6,2) * t589 - t122 + t566) * t509 + (-t125 * t576 + t121 - t462 + t64) * t507 + (-pkin(5) * t3 - t391 * g(1) - t392 * g(2) + qJ(6) * t2 - t23 * t28 + t24 * t597 - t50 * t73) * m(7); -t304 * t107 + t589 * t74 + (-g(1) * t222 - g(2) * t220 - t24 * t304 - t311 * t476 + t50 * t589 + t3) * m(7) + t32;];
tau  = t1;
