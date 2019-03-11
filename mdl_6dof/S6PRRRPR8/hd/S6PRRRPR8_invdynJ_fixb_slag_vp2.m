% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:15
% EndTime: 2019-03-08 23:49:16
% DurationCPUTime: 37.23s
% Computational Cost: add. (11509->903), mult. (28934->1214), div. (0->0), fcn. (24154->14), ass. (0->435)
t317 = cos(qJ(3));
t308 = sin(pkin(7));
t470 = qJD(2) * t308;
t445 = t317 * t470;
t653 = qJD(4) - t445;
t652 = Ifges(5,4) + Ifges(6,6);
t309 = sin(pkin(6));
t314 = sin(qJ(2));
t483 = t309 * t314;
t447 = qJD(1) * t483;
t267 = pkin(9) * t470 + t447;
t313 = sin(qJ(3));
t534 = cos(qJ(2));
t449 = t309 * t534;
t405 = qJD(1) * t449;
t279 = qJD(2) * pkin(2) + t405;
t310 = cos(pkin(6));
t488 = t308 * t310;
t500 = cos(pkin(7));
t641 = qJD(1) * t488 + t500 * t279;
t150 = t317 * t267 + t313 * t641;
t312 = sin(qJ(4));
t464 = qJD(4) * t312;
t651 = -qJD(5) * t312 - t150 + (-t312 * t445 + t464) * pkin(4);
t631 = m(7) + m(6);
t316 = cos(qJ(4));
t459 = qJD(2) * qJD(3);
t352 = qJDD(2) * t313 + t317 * t459;
t341 = t352 * t308;
t367 = qJD(2) * t500 + qJD(3);
t349 = qJD(4) * t367;
t363 = qJDD(2) * t500 + qJDD(3);
t486 = t308 * t313;
t455 = t312 * t486;
t406 = qJD(4) * t455;
t131 = qJD(2) * t406 - t312 * t363 + (-t341 - t349) * t316;
t554 = -t131 / 0.2e1;
t469 = qJD(2) * t313;
t443 = t316 * t469;
t132 = t312 * t349 - t316 * t363 + (qJD(4) * t443 + t312 * t352) * t308;
t552 = -t132 / 0.2e1;
t252 = (-qJDD(2) * t317 + t313 * t459) * t308;
t239 = qJDD(4) + t252;
t539 = t239 / 0.2e1;
t621 = Ifges(6,4) - Ifges(5,5);
t620 = Ifges(5,6) - Ifges(6,5);
t619 = Ifges(5,3) + Ifges(6,1);
t422 = t313 * t500;
t340 = -t314 * t422 + t317 * t534;
t223 = t340 * t309;
t198 = qJD(1) * t223;
t409 = t308 * t447;
t162 = t198 * t312 - t316 * t409;
t421 = t316 * t500;
t467 = qJD(3) * t317;
t441 = t308 * t467;
t188 = -qJD(4) * t421 - t316 * t441 + t406;
t468 = qJD(3) * t313;
t442 = t308 * t468;
t558 = pkin(4) + pkin(11);
t484 = t308 * t317;
t265 = pkin(2) * t422 + pkin(9) * t484;
t235 = pkin(10) * t500 + t265;
t388 = -pkin(3) * t317 - pkin(10) * t313;
t236 = (-pkin(2) + t388) * t308;
t361 = (pkin(3) * t313 - pkin(10) * t317) * t308;
t247 = qJD(3) * t361;
t420 = t317 * t500;
t264 = pkin(2) * t420 - pkin(9) * t486;
t248 = t264 * qJD(3);
t463 = qJD(4) * t316;
t77 = -t235 * t463 - t236 * t464 + t247 * t316 - t312 * t248;
t650 = -pkin(5) * t188 - t442 * t558 - t162 - t77;
t485 = t308 * t316;
t263 = t312 * t500 + t313 * t485;
t189 = qJD(4) * t263 + t312 * t441;
t339 = t313 * t534 + t314 * t420;
t222 = t339 * t309;
t197 = qJD(1) * t222;
t249 = t265 * qJD(3);
t325 = t188 * qJ(5) - t263 * qJD(5) + t249;
t649 = -t189 * t558 + t197 - t325;
t648 = t651 + t653 * (pkin(11) * t312 - qJ(5) * t316);
t557 = pkin(5) + pkin(10);
t294 = t557 * t316;
t472 = t316 * t317;
t149 = -t313 * t267 + t317 * t641;
t246 = qJD(2) * t361;
t96 = -t312 * t149 + t246 * t316;
t647 = qJD(4) * t294 - (pkin(5) * t472 - t313 * t558) * t470 + t96;
t127 = qJDD(6) - t131;
t446 = t308 * t469;
t220 = t312 * t446 - t316 * t367;
t311 = sin(qJ(6));
t315 = cos(qJ(6));
t167 = t220 * t315 - t311 * t653;
t49 = qJD(6) * t167 + t132 * t311 + t239 * t315;
t25 = mrSges(7,1) * t127 - mrSges(7,3) * t49;
t168 = t220 * t311 + t315 * t653;
t50 = -qJD(6) * t168 + t132 * t315 - t239 * t311;
t26 = -mrSges(7,2) * t127 + mrSges(7,3) * t50;
t221 = t308 * t443 + t312 * t367;
t209 = qJD(6) + t221;
t107 = -mrSges(7,2) * t209 + mrSges(7,3) * t167;
t108 = mrSges(7,1) * t209 - mrSges(7,3) * t168;
t369 = t315 * t107 - t311 * t108;
t646 = t369 * qJD(6) + t315 * t25 + t311 * t26;
t385 = mrSges(7,1) * t311 + mrSges(7,2) * t315;
t645 = -t239 / 0.2e1;
t207 = Ifges(5,4) * t220;
t612 = t221 * Ifges(5,1) + Ifges(5,5) * t653 + t168 * Ifges(7,5) + t167 * Ifges(7,6) + t209 * Ifges(7,3) - t207;
t423 = t310 * t500;
t389 = qJD(1) * t423;
t326 = t389 + (qJD(2) * t388 - t279) * t308;
t295 = qJD(2) * t405;
t254 = qJDD(1) * t483 + t295;
t633 = pkin(9) * qJDD(2) * t308 + qJD(3) * t641 + t254;
t444 = qJD(2) * t483;
t404 = qJD(1) * t444;
t253 = qJDD(1) * t449 - t404;
t232 = qJDD(2) * pkin(2) + t253;
t640 = qJDD(1) * t488 + t500 * t232;
t54 = -t267 * t468 + t313 * t640 + t317 * t633;
t642 = pkin(10) * t363 + qJD(4) * t326 + t54;
t134 = pkin(10) * t367 + t150;
t62 = t316 * t134 + t312 * t326;
t43 = -t220 * pkin(5) + t62;
t575 = -t385 + mrSges(5,2) - mrSges(6,3);
t587 = m(7) * pkin(11) + mrSges(5,1) - mrSges(6,2);
t639 = pkin(4) * t631 + t587;
t133 = -pkin(3) * t367 - t149;
t321 = -t221 * qJ(5) + t133;
t63 = t220 * pkin(4) + t321;
t638 = -t133 * mrSges(5,2) + t63 * mrSges(6,3);
t55 = -t267 * t467 - t313 * t633 + t317 * t640;
t637 = t55 * mrSges(4,1) - t54 * mrSges(4,2);
t61 = t134 * t312 - t316 * t326;
t362 = pkin(5) * t221 + t61;
t635 = t362 + qJD(5);
t10 = Ifges(7,5) * t49 + Ifges(7,6) * t50 + Ifges(7,3) * t127;
t553 = t131 / 0.2e1;
t555 = t127 / 0.2e1;
t561 = t50 / 0.2e1;
t562 = t49 / 0.2e1;
t29 = -t558 * t653 + t635;
t48 = t220 * t558 + t321;
t17 = t29 * t315 - t311 * t48;
t297 = qJDD(1) * t423;
t95 = t252 * pkin(3) + t297 + (-pkin(10) * t352 - t232) * t308;
t14 = -t134 * t463 - t312 * t642 + t316 * t95;
t353 = qJDD(5) - t14;
t5 = -pkin(5) * t131 - t239 * t558 + t353;
t47 = -pkin(3) * t363 - t55;
t320 = t131 * qJ(5) - t221 * qJD(5) + t47;
t9 = t132 * t558 + t320;
t1 = qJD(6) * t17 + t311 * t5 + t315 * t9;
t18 = t29 * t311 + t315 * t48;
t2 = -qJD(6) * t18 - t311 * t9 + t315 * t5;
t584 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t634 = 0.2e1 * Ifges(5,1) * t554 + Ifges(6,4) * t645 + Ifges(7,5) * t562 + Ifges(7,6) * t561 + Ifges(7,3) * t555 + t584 + (-t553 + t554) * Ifges(6,2) + t10 / 0.2e1 + t652 * t552 + (-t621 + Ifges(5,5)) * t539;
t551 = t132 / 0.2e1;
t632 = -Ifges(6,5) * t645 - 0.2e1 * Ifges(5,2) * t552 - t652 * t554 - (-t551 + t552) * Ifges(6,3) + (-t620 - Ifges(5,6)) * t539;
t630 = -t252 / 0.2e1;
t629 = t341 / 0.2e1;
t628 = t363 / 0.2e1;
t151 = -t312 * t235 + t236 * t316;
t143 = pkin(4) * t484 - t151;
t86 = pkin(5) * t263 + pkin(11) * t484 + t143;
t262 = -t421 + t455;
t234 = -pkin(3) * t500 - t264;
t330 = -t263 * qJ(5) + t234;
t94 = t262 * t558 + t330;
t28 = t311 * t86 + t315 * t94;
t627 = -qJD(6) * t28 + t649 * t311 + t315 * t650;
t27 = -t311 * t94 + t315 * t86;
t626 = qJD(6) * t27 + t311 * t650 - t649 * t315;
t82 = mrSges(5,1) * t239 + mrSges(5,3) * t131;
t85 = -t131 * mrSges(6,1) + t239 * mrSges(6,2);
t617 = -t82 + t85;
t83 = -mrSges(5,2) * t239 - mrSges(5,3) * t132;
t84 = mrSges(6,1) * t132 - mrSges(6,3) * t239;
t616 = -t84 + t83;
t496 = qJ(5) * t312;
t429 = -pkin(3) - t496;
t266 = -t316 * t558 + t429;
t293 = t557 * t312;
t182 = -t266 * t311 + t293 * t315;
t615 = qJD(6) * t182 + t311 * t647 + t315 * t648;
t183 = t266 * t315 + t293 * t311;
t614 = -qJD(6) * t183 - t311 * t648 + t315 * t647;
t522 = mrSges(6,1) * t220;
t173 = -mrSges(6,3) * t653 + t522;
t78 = -mrSges(7,1) * t167 + mrSges(7,2) * t168;
t611 = -t173 + t78;
t478 = t312 * t317;
t97 = t316 * t149 + t312 * t246;
t610 = -t557 * t464 - (-pkin(5) * t478 + qJ(5) * t313) * t470 - t97;
t498 = sin(pkin(12));
t391 = t498 * t534;
t499 = cos(pkin(12));
t417 = t499 * t314;
t260 = t310 * t417 + t391;
t392 = t499 * t534;
t416 = t498 * t314;
t333 = -t310 * t392 + t416;
t425 = t309 * t499;
t600 = t308 * t425 + t333 * t500;
t135 = t260 * t313 + t317 * t600;
t495 = t135 * t316;
t608 = -pkin(4) * t495 - t135 * t496;
t261 = -t310 * t416 + t392;
t334 = t310 * t391 + t417;
t424 = t309 * t498;
t599 = -t308 * t424 + t334 * t500;
t137 = t261 * t313 + t317 * t599;
t494 = t137 * t316;
t607 = -pkin(4) * t494 - t137 * t496;
t606 = -t620 * t220 - t621 * t221 + t619 * t653;
t393 = t500 * t534;
t597 = -t313 * t314 + t317 * t393;
t184 = -t309 * t597 - t310 * t484;
t493 = t184 * t316;
t605 = -pkin(4) * t493 - t184 * t496;
t474 = t315 * t317;
t202 = (-t311 * t313 + t312 * t474) * t470;
t461 = qJD(6) * t316;
t440 = t311 * t461;
t604 = -t315 * t464 + t202 - t440;
t203 = (t311 * t478 + t313 * t315) * t470;
t603 = -t311 * t464 + t315 * t461 + t203;
t602 = t248 - t198;
t407 = t316 * t445;
t601 = (t407 - t463) * qJ(5) + t651;
t598 = -t63 * (-mrSges(6,2) * t312 - mrSges(6,3) * t316) - t133 * (mrSges(5,1) * t312 + mrSges(5,2) * t316);
t596 = -t620 * t312 - t621 * t316;
t594 = t621 * t131 - t620 * t132 + t619 * t239;
t13 = -t134 * t464 + t312 * t95 + t316 * t642;
t592 = t13 * t316 - t14 * t312;
t7 = -qJ(5) * t239 - qJD(5) * t653 - t13;
t8 = -pkin(4) * t239 + t353;
t591 = t312 * t8 - t316 * t7;
t205 = -t308 * t279 + t389;
t590 = (t367 * (Ifges(4,5) * t317 - Ifges(4,6) * t313) / 0.2e1 + t205 * (mrSges(4,1) * t313 + mrSges(4,2) * t317)) * t308;
t489 = t653 * qJ(5);
t32 = t489 + t43;
t386 = t315 * mrSges(7,1) - t311 * mrSges(7,2);
t164 = Ifges(7,4) * t167;
t60 = Ifges(7,1) * t168 + Ifges(7,5) * t209 + t164;
t502 = t311 * t60;
t588 = t32 * t386 - t502 / 0.2e1;
t350 = t235 * t464 - t236 * t463 - t312 * t247 - t316 * t248;
t57 = -t308 * (qJ(5) * t468 - qJD(5) * t317) + t350;
t520 = mrSges(5,3) * t220;
t171 = -mrSges(5,2) * t653 - t520;
t56 = -t62 - t489;
t586 = -m(6) * t56 + t171 - t173;
t519 = mrSges(5,3) * t221;
t172 = mrSges(5,1) * t653 - t519;
t521 = mrSges(6,1) * t221;
t174 = mrSges(6,2) * t653 + t521;
t51 = -pkin(4) * t653 + qJD(5) + t61;
t585 = -m(6) * t51 + t172 - t174;
t414 = mrSges(4,3) * t446;
t583 = m(5) * t133 - mrSges(4,1) * t367 + mrSges(5,1) * t220 + mrSges(5,2) * t221 + t414;
t582 = t17 * mrSges(7,1) - t18 * mrSges(7,2);
t581 = -mrSges(7,3) - t587;
t580 = -t133 * mrSges(5,1) + t63 * mrSges(6,2);
t579 = -t631 * qJ(5) + t575;
t206 = Ifges(6,6) * t220;
t117 = Ifges(6,4) * t653 - t221 * Ifges(6,2) + t206;
t578 = -t117 / 0.2e1 + t582;
t577 = m(5) * t61 - t585;
t338 = t313 * t393 + t314 * t317;
t185 = t309 * t338 + t310 * t486;
t259 = -t308 * t449 + t423;
t140 = t185 * t312 - t259 * t316;
t136 = t260 * t317 - t313 * t600;
t322 = t308 * t333 - t425 * t500;
t71 = t136 * t312 - t316 * t322;
t138 = t261 * t317 - t313 * t599;
t323 = t308 * t334 + t424 * t500;
t73 = t138 * t312 - t316 * t323;
t576 = g(1) * t73 + g(2) * t71 + g(3) * t140;
t384 = mrSges(6,2) * t316 - mrSges(6,3) * t312;
t387 = -mrSges(5,1) * t316 + mrSges(5,2) * t312;
t573 = t312 * t385 + mrSges(4,1) - t384 - t387;
t572 = -t14 * mrSges(5,1) + t13 * mrSges(5,2) - t8 * mrSges(6,2) + t7 * mrSges(6,3);
t571 = m(5) * t62 + m(7) * t32 + t586 + t78;
t148 = -mrSges(6,2) * t220 - mrSges(6,3) * t221;
t570 = m(4) * t149 - m(6) * t63 - t148 - t583;
t569 = -m(7) * t557 - mrSges(6,1) + mrSges(4,2) - mrSges(5,3) - t386;
t545 = -t209 / 0.2e1;
t548 = -t168 / 0.2e1;
t550 = -t167 / 0.2e1;
t568 = -Ifges(7,5) * t548 - Ifges(7,6) * t550 - Ifges(7,3) * t545 + t582;
t566 = t308 ^ 2;
t319 = qJD(2) ^ 2;
t11 = Ifges(7,4) * t49 + Ifges(7,2) * t50 + Ifges(7,6) * t127;
t565 = -t11 / 0.2e1;
t12 = Ifges(7,1) * t49 + Ifges(7,4) * t50 + Ifges(7,5) * t127;
t564 = t12 / 0.2e1;
t513 = Ifges(7,4) * t168;
t59 = Ifges(7,2) * t167 + Ifges(7,6) * t209 + t513;
t560 = -t59 / 0.2e1;
t559 = t59 / 0.2e1;
t505 = t221 * Ifges(5,4);
t118 = -t220 * Ifges(5,2) + Ifges(5,6) * t653 + t505;
t556 = -t118 / 0.2e1;
t549 = t167 / 0.2e1;
t547 = t168 / 0.2e1;
t544 = t209 / 0.2e1;
t543 = -t220 / 0.2e1;
t542 = t220 / 0.2e1;
t541 = -t221 / 0.2e1;
t540 = t221 / 0.2e1;
t537 = t653 / 0.2e1;
t536 = -t653 / 0.2e1;
t533 = pkin(2) * t308;
t158 = t260 * t420 - t313 * t333;
t532 = pkin(10) * t158;
t160 = t261 * t420 - t313 * t334;
t531 = pkin(10) * t160;
t530 = pkin(10) * t222;
t527 = t2 * t315;
t518 = mrSges(7,3) * t311;
t517 = Ifges(4,4) * t313;
t516 = Ifges(4,4) * t317;
t515 = Ifges(5,4) * t312;
t514 = Ifges(5,4) * t316;
t512 = Ifges(7,4) * t311;
t511 = Ifges(7,4) * t315;
t510 = Ifges(6,6) * t312;
t509 = Ifges(6,6) * t316;
t506 = t17 * t311;
t504 = t221 * Ifges(6,6);
t497 = qJ(5) * t220;
t492 = t221 * t315;
t491 = t260 * t308;
t490 = t261 * t308;
t487 = t308 * t312;
t480 = t311 * t316;
t475 = t315 * t316;
t152 = t316 * t235 + t312 * t236;
t456 = t308 * t483;
t471 = pkin(2) * t449 + pkin(9) * t456;
t462 = qJD(6) * t315;
t458 = pkin(10) * t464;
t457 = pkin(10) * t463;
t452 = t223 * pkin(3) + t471;
t450 = Ifges(4,5) * t341 - Ifges(4,6) * t252 + Ifges(4,3) * t363;
t439 = t486 / 0.2e1;
t436 = -t470 / 0.2e1;
t435 = t470 / 0.2e1;
t430 = -t462 / 0.2e1;
t128 = t135 * pkin(3);
t428 = pkin(10) * t136 - t128;
t129 = t137 * pkin(3);
t427 = pkin(10) * t138 - t129;
t181 = t184 * pkin(3);
t426 = pkin(10) * t185 - t181;
t413 = mrSges(4,3) * t445;
t408 = t308 * t444;
t400 = t317 * t436;
t399 = t317 * t435;
t397 = -t333 * pkin(2) + pkin(9) * t491;
t396 = -t334 * pkin(2) + pkin(9) * t490;
t383 = Ifges(5,1) * t316 - t515;
t382 = Ifges(7,1) * t315 - t512;
t381 = Ifges(7,1) * t311 + t511;
t380 = -Ifges(5,2) * t312 + t514;
t378 = -Ifges(7,2) * t311 + t511;
t377 = Ifges(7,2) * t315 + t512;
t375 = Ifges(7,5) * t315 - Ifges(7,6) * t311;
t374 = Ifges(7,5) * t311 + Ifges(7,6) * t315;
t373 = -Ifges(6,2) * t316 + t510;
t372 = Ifges(6,3) * t312 - t509;
t370 = -t18 * t315 + t506;
t69 = t140 * t315 - t184 * t311;
t70 = t140 * t311 + t184 * t315;
t141 = t185 * t316 + t259 * t312;
t159 = -t260 * t422 - t317 * t333;
t366 = t159 * pkin(3) + t397;
t161 = -t261 * t422 - t317 * t334;
t365 = t161 * pkin(3) + t396;
t142 = qJ(5) * t484 - t152;
t360 = -t262 * t311 + t308 * t474;
t190 = t262 * t315 + t311 * t484;
t357 = (-mrSges(4,1) * t317 + mrSges(4,2) * t313) * t308;
t169 = t223 * t312 - t316 * t456;
t170 = t223 * t316 + t312 * t456;
t354 = t170 * pkin(4) + qJ(5) * t169 + t452;
t346 = t313 * t566 * (Ifges(4,1) * t317 - t517);
t103 = t159 * t312 - t260 * t485;
t104 = t159 * t316 + t260 * t487;
t337 = t104 * pkin(4) + qJ(5) * t103 + t366;
t105 = t161 * t312 - t261 * t485;
t106 = t161 * t316 + t261 * t487;
t336 = t106 * pkin(4) + qJ(5) * t105 + t365;
t332 = Ifges(4,6) * t500 + (t317 * Ifges(4,2) + t517) * t308;
t331 = -qJD(6) * t370 + t1 * t311 + t527;
t296 = Ifges(4,4) * t445;
t286 = -pkin(4) * t316 + t429;
t245 = qJD(2) * t357;
t244 = -mrSges(4,2) * t367 + t413;
t193 = Ifges(4,1) * t446 + Ifges(4,5) * t367 + t296;
t192 = Ifges(4,6) * qJD(3) + qJD(2) * t332;
t187 = mrSges(4,1) * t363 - mrSges(4,3) * t341;
t186 = -mrSges(4,2) * t363 - t252 * mrSges(4,3);
t180 = -t308 * t232 + t297;
t165 = t252 * mrSges(4,1) + mrSges(4,2) * t341;
t146 = pkin(4) * t221 + t497;
t139 = t262 * pkin(4) + t330;
t122 = t310 * t441 + (t340 * qJD(2) + qJD(3) * t597) * t309;
t121 = t310 * t442 + (qJD(2) * t339 + qJD(3) * t338) * t309;
t115 = Ifges(6,5) * t653 + t220 * Ifges(6,3) - t504;
t101 = t221 * t558 + t497;
t98 = -pkin(5) * t262 - t142;
t91 = qJD(6) * t360 + t189 * t315 - t311 * t442;
t90 = qJD(6) * t190 + t189 * t311 + t315 * t442;
t81 = -pkin(4) * t446 - t96;
t79 = -qJ(5) * t446 - t97;
t66 = t189 * pkin(4) + t325;
t64 = -pkin(4) * t442 - t77;
t53 = mrSges(5,1) * t132 - mrSges(5,2) * t131;
t52 = -mrSges(6,2) * t132 + mrSges(6,3) * t131;
t37 = qJD(4) * t141 + t122 * t312 - t316 * t408;
t31 = -pkin(5) * t189 - t57;
t22 = t101 * t315 + t311 * t43;
t21 = -t101 * t311 + t315 * t43;
t20 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t19 = t132 * pkin(4) + t320;
t16 = qJD(6) * t69 + t121 * t315 + t311 * t37;
t15 = -qJD(6) * t70 - t121 * t311 + t315 * t37;
t6 = -pkin(5) * t132 - t7;
t3 = [m(3) * (t310 ^ 2 * qJDD(1) + (t253 * t534 + t254 * t314) * t309) + m(4) * (t122 * t150 + t180 * t259 + t185 * t54 + t205 * t408) + t259 * t165 + t122 * t244 + t185 * t186 + t16 * t107 + t15 * t108 + t69 * t25 + t70 * t26 + m(7) * (t1 * t70 + t15 * t17 + t16 * t18 + t2 * t69) + t245 * t408 + m(2) * qJDD(1) + t577 * t37 + (-m(5) * t14 + m(6) * t8 + t617) * t140 + (-qJDD(2) * t483 - t319 * t449) * mrSges(3,2) + (qJDD(2) * t449 - t319 * t483) * mrSges(3,1) + t571 * (t312 * t408 - t185 * t464 + (qJD(4) * t259 + t122) * t316) + (-m(4) * t55 + m(5) * t47 + m(6) * t19 - t187 + t52 + t53) * t184 + (m(5) * t13 - m(6) * t7 + m(7) * t6 + t20 + t616) * t141 - t570 * t121 + (-m(2) - m(3) - m(4) - m(5) - t631) * g(3); m(5) * (t13 * t152 + t14 * t151 + t234 * t47 - t350 * t62 - t61 * t77) - t350 * t171 + (-Ifges(7,5) * t360 + Ifges(7,6) * t190) * t555 + (-Ifges(7,4) * t360 + Ifges(7,2) * t190) * t561 + (t1 * t190 - t17 * t90 + t18 * t91 + t2 * t360) * mrSges(7,3) + (-Ifges(7,1) * t360 + Ifges(7,4) * t190) * t562 + t6 * (-mrSges(7,1) * t190 - mrSges(7,2) * t360) - t360 * t564 + (t51 * mrSges(6,2) - t61 * mrSges(5,1) - t192 / 0.2e1 - t56 * mrSges(6,3) - t62 * mrSges(5,2) + Ifges(5,5) * t540 + Ifges(6,4) * t541 + Ifges(6,5) * t542 + Ifges(5,6) * t543 + t619 * t537) * t442 + (t56 * mrSges(6,1) - t62 * mrSges(5,3) + t556 - Ifges(5,4) * t540 + Ifges(6,6) * t541 + Ifges(6,3) * t542 - Ifges(5,2) * t543 + t115 / 0.2e1 - t620 * t537 - t580) * t189 + (Ifges(7,5) * t90 + Ifges(7,6) * t91) * t544 + (Ifges(4,1) * t341 - Ifges(4,4) * t252 + Ifges(4,5) * t363) * t439 + (Ifges(7,4) * t90 + Ifges(7,2) * t91) * t549 + (Ifges(7,1) * t90 + Ifges(7,4) * t91) * t547 + t626 * t107 + t627 * t108 + (t1 * t28 + t17 * t627 + t18 * t626 + t2 * t27 + t31 * t32 + t6 * t98) * m(7) + t583 * t249 + (-m(4) * t396 - t161 * mrSges(4,1) - mrSges(4,3) * t490 - m(7) * t336 + mrSges(3,1) * t334 + t261 * mrSges(3,2) - m(6) * (t336 + t531) - m(5) * (t365 + t531) + t575 * t105 + t569 * t160 + t581 * t106) * g(1) + (-m(4) * t471 - t223 * mrSges(4,1) - mrSges(4,3) * t456 - m(7) * t354 - (mrSges(3,1) * t534 - mrSges(3,2) * t314) * t309 - m(6) * (t354 + t530) - m(5) * (t452 + t530) + t575 * t169 + t569 * t222 + t581 * t170) * g(3) + (-m(4) * t397 - t159 * mrSges(4,1) - mrSges(4,3) * t491 - m(7) * t337 - m(6) * (t337 + t532) - m(5) * (t366 + t532) + mrSges(3,1) * t333 + t260 * mrSges(3,2) + t575 * t103 + t569 * t158 + t581 * t104) * g(2) + (t439 * t606 + t590) * qJD(3) - t577 * t162 + t570 * t197 - t571 * (t198 * t316 + t312 * t409) + (t566 * qJD(2) * (-Ifges(4,2) * t313 + t516) + t308 * t193) * t467 / 0.2e1 + (-t149 * t441 - t150 * t442 - t486 * t55) * mrSges(4,3) + (t54 * mrSges(4,3) + Ifges(4,4) * t629 - Ifges(6,4) * t553 - Ifges(5,5) * t554 - Ifges(6,5) * t551 + Ifges(4,2) * t630 + Ifges(4,6) * t628 - Ifges(5,6) * t552 - t539 * t619 + t572 - t594 / 0.2e1) * t484 + (-t254 + t295) * mrSges(3,2) + t602 * t244 + (-t149 * t249 + t150 * t602 - t180 * t533 - t205 * t409 + t264 * t55 + t265 * t54) * m(4) + (Ifges(4,5) * t313 + Ifges(4,6) * t317) * t308 * t628 + t346 * t459 / 0.2e1 + (t253 + t404) * mrSges(3,1) + t91 * t559 + (t313 * Ifges(4,1) + t516) * t308 * t629 + t180 * t357 - t165 * t533 + t264 * t187 + t265 * t186 + t234 * t53 + t190 * t11 / 0.2e1 + Ifges(3,3) * qJDD(2) + t77 * t172 + t57 * t173 + t64 * t174 + t151 * t82 + t152 * t83 + t142 * t84 + t143 * t85 + t66 * t148 + t139 * t52 + t98 * t20 + t90 * t60 / 0.2e1 + t32 * (-mrSges(7,1) * t91 + mrSges(7,2) * t90) - t245 * t409 + t31 * t78 + t27 * t25 + t28 * t26 + m(6) * (t139 * t19 + t142 * t7 + t143 * t8 + t51 * t64 + t56 * t57 + t63 * t66) + (t47 * mrSges(5,1) + t7 * mrSges(6,1) - t19 * mrSges(6,2) - t13 * mrSges(5,3) - Ifges(5,4) * t554 + Ifges(6,6) * t553 + t632) * t262 + (mrSges(6,1) * t8 + mrSges(5,2) * t47 - mrSges(5,3) * t14 - mrSges(6,3) * t19 + Ifges(5,4) * t552 - Ifges(6,6) * t551 + t634) * t263 + (t450 / 0.2e1 + Ifges(4,3) * t628 + Ifges(4,5) * t629 + t637) * t500 + (-t51 * mrSges(6,1) - t61 * mrSges(5,3) - Ifges(5,1) * t540 - Ifges(5,4) * t543 - Ifges(7,5) * t547 + Ifges(6,2) * t541 + Ifges(6,6) * t542 - Ifges(7,6) * t549 - Ifges(7,3) * t544 + t537 * t621 - t578 - t612 / 0.2e1 + t638) * t188 + t332 * t630; (t596 * t537 + (t383 / 0.2e1 - t373 / 0.2e1) * t221 + (-t380 / 0.2e1 + t372 / 0.2e1) * t220 + (Ifges(7,3) * t316 + t312 * t374) * t544 + (Ifges(7,5) * t316 + t312 * t381) * t547 + (Ifges(7,6) * t316 + t312 * t377) * t549 - t598) * qJD(4) - t509 * t553 + t514 * t554 + t614 * t108 + (t1 * t183 + t17 * t614 + t18 * t615 + t182 * t2 + t294 * t6 + t32 * t610) * m(7) + t615 * t107 + t610 * t78 + (Ifges(7,4) * t203 + Ifges(7,2) * t202) * t550 + (t61 * (mrSges(5,1) * t313 - mrSges(5,3) * t472) - t51 * (mrSges(6,1) * t472 + mrSges(6,2) * t313) - t62 * (-mrSges(5,2) * t313 - mrSges(5,3) * t478) - t56 * (mrSges(6,1) * t478 - mrSges(6,3) * t313)) * t470 + (t457 - t81) * t174 + (t458 - t79) * t173 + (-t457 - t96) * t172 + (t414 - t583) * t150 + (-t458 - t97) * t171 + (Ifges(7,5) * t203 + Ifges(7,6) * t202) * t545 + (t221 * (Ifges(5,5) * t313 + t317 * t383) + t220 * (Ifges(6,5) * t313 + t317 * t372) + t606 * t313) * t436 + (-t244 + t413) * t149 + t578 * t463 - t568 * t407 + (t315 * t59 + t115 + t502) * t464 / 0.2e1 + (-t375 * t544 - t378 * t549 - t382 * t547) * t461 + (t19 * t286 - t51 * t81 - t56 * t79 + t601 * t63) * m(6) + (-pkin(3) * t47 + t61 * t96 - t62 * t97) * m(5) - t510 * t551 + t515 * t552 + (Ifges(7,1) * t203 + Ifges(7,4) * t202) * t548 + (-Ifges(4,2) * t446 + t193 + t296) * t400 - t590 * qJD(2) + (t221 * (Ifges(6,4) * t313 + t317 * t373) + t220 * (Ifges(5,6) * t313 + t317 * t380) + t313 * t192 - (t313 * t619 + t317 * t596) * t653) * t435 + t601 * t148 + (mrSges(7,1) * t604 - mrSges(7,2) * t603) * t32 + (-t1 * t475 + t17 * t603 - t18 * t604 + t2 * t480) * mrSges(7,3) + (-m(7) * (-pkin(11) * t493 - t181 + t605) + mrSges(7,3) * t493 - m(6) * (t426 + t605) - m(5) * t426 + t569 * t185 + t573 * t184) * g(3) + (t463 * t51 + t464 * t56 + t591) * mrSges(6,1) + (t463 * t61 - t464 * t62 + t592) * mrSges(5,3) + t598 * t445 + t450 + (-m(7) * (-pkin(11) * t494 - t129 + t607) + mrSges(7,3) * t494 - m(6) * (t427 + t607) - m(5) * t427 + t569 * t138 + t573 * t137) * g(1) + (-m(7) * (-pkin(11) * t495 - t128 + t608) + mrSges(7,3) * t495 - m(6) * (t428 + t608) - m(5) * t428 + t569 * t136 + t573 * t135) * g(2) + t464 * t556 + t440 * t559 + t202 * t560 + t475 * t565 + t294 * t20 + t286 * t52 - t12 * t480 / 0.2e1 + t637 + (((t312 * t56 + t316 * t51) * qJD(4) + t591) * m(6) + ((-t312 * t62 + t316 * t61) * qJD(4) + t592) * m(5) + t616 * t316 + t617 * t312) * pkin(10) - t203 * t60 / 0.2e1 + t182 * t25 + t183 * t26 - pkin(3) * t53 + t19 * t384 + (t117 * t399 - t374 * t555 - t377 * t561 - t381 * t562 + t6 * t386 + t60 * t430 - t632) * t316 + t47 * t387 + (t115 * t400 + t118 * t399 + t634) * t312 + t612 * (t316 * t400 + t463 / 0.2e1) - t319 * t346 / 0.2e1; t362 * t78 + (Ifges(6,3) * t543 + t17 * t518 + t374 * t545 + t377 * t550 + t381 * t548 - t536 * t620 + t580 + t588) * t221 + t611 * qJD(5) + (-Ifges(5,2) * t221 - t207 + t612) * t542 - (t167 * t377 + t168 * t381 + t209 * t374) * qJD(6) / 0.2e1 + t594 + t588 * qJD(6) + (t519 + t585) * t62 + (t520 + t586) * t61 + (qJD(6) * t506 - t527 + (-t462 - t492) * t18 + t576) * mrSges(7,3) + (-t84 + t20) * qJ(5) + (-pkin(4) * t8 - qJ(5) * t7 - qJD(5) * t56 - t146 * t63) * m(6) - t572 + (-m(7) * t331 - t646) * t558 - t1 * t518 + t59 * t430 + (qJ(5) * t6 - t17 * t21 - t18 * t22 + t32 * t635) * m(7) + t375 * t555 + t492 * t560 + t378 * t561 + t382 * t562 + t315 * t564 + t311 * t565 - t56 * t521 - t146 * t148 - t22 * t107 - t21 * t108 - pkin(4) * t85 + (-t505 + t115) * t541 + (t504 + t118) * t540 + (t206 + t117) * t543 + t6 * t385 + (-Ifges(5,1) * t541 + Ifges(6,2) * t540 + t536 * t621 + t568 - t638) * t220 + (t579 * (t136 * t316 + t312 * t322) + t639 * t71) * g(2) + (t579 * (t138 * t316 + t312 * t323) + t639 * t73) * g(1) + (t140 * t639 + t141 * t579) * g(3) + t51 * t522; -t611 * t653 + (t148 + t369) * t221 + t85 + (-t221 * t370 - t32 * t653 + t331 - t576) * m(7) + (t221 * t63 + t56 * t653 - t576 + t8) * m(6) + t646; -t32 * (mrSges(7,1) * t168 + mrSges(7,2) * t167) + (Ifges(7,1) * t167 - t513) * t548 + t59 * t547 + (Ifges(7,5) * t167 - Ifges(7,6) * t168) * t545 - t17 * t107 + t18 * t108 - g(1) * ((-t137 * t311 + t315 * t73) * mrSges(7,1) + (-t137 * t315 - t311 * t73) * mrSges(7,2)) - g(2) * ((-t135 * t311 + t315 * t71) * mrSges(7,1) + (-t135 * t315 - t311 * t71) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t69 - mrSges(7,2) * t70) + (t167 * t17 + t168 * t18) * mrSges(7,3) + t10 + (-Ifges(7,2) * t168 + t164 + t60) * t550 + t584;];
tau  = t3;
