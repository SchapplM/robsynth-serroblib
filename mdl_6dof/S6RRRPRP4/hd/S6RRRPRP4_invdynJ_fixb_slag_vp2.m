% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:36
% EndTime: 2019-03-09 16:44:16
% DurationCPUTime: 24.56s
% Computational Cost: add. (10331->735), mult. (22109->878), div. (0->0), fcn. (14990->10), ass. (0->344)
t610 = m(6) + m(7);
t601 = -mrSges(6,3) - mrSges(7,2);
t609 = Ifges(4,5) - Ifges(5,4);
t608 = mrSges(4,1) - mrSges(5,2);
t573 = Ifges(6,1) + Ifges(7,1);
t607 = Ifges(4,4) + Ifges(5,6);
t606 = -Ifges(6,4) + Ifges(7,5);
t572 = Ifges(7,4) + Ifges(6,5);
t605 = Ifges(6,6) - Ifges(7,6);
t570 = Ifges(6,3) + Ifges(7,2);
t520 = pkin(3) + pkin(9);
t604 = t520 * t610 - t601;
t575 = mrSges(6,1) + mrSges(7,1);
t530 = -m(7) * pkin(5) - t575;
t318 = sin(qJ(2));
t502 = cos(qJ(2));
t400 = qJD(2) * t502;
t247 = qJD(1) * t400 + t318 * qJDD(1);
t317 = sin(qJ(3));
t501 = cos(qJ(3));
t242 = t317 * t502 + t318 * t501;
t332 = t242 * qJD(3);
t428 = qJD(1) * qJD(2);
t591 = t502 * qJDD(1) - t318 * t428;
t146 = qJD(1) * t332 + t317 * t247 - t501 * t591;
t313 = qJDD(2) + qJDD(3);
t316 = sin(qJ(5));
t320 = cos(qJ(5));
t383 = t501 * t502;
t436 = qJD(1) * t318;
t226 = -qJD(1) * t383 + t317 * t436;
t314 = qJD(2) + qJD(3);
t356 = t320 * t226 - t314 * t316;
t585 = qJD(5) * t356;
t80 = t146 * t316 + t313 * t320 + t585;
t523 = t80 / 0.2e1;
t199 = t226 * t316 + t314 * t320;
t81 = qJD(5) * t199 - t320 * t146 + t313 * t316;
t522 = -t81 / 0.2e1;
t145 = qJD(3) * t226 - t501 * t247 - t317 * t591;
t143 = qJDD(5) - t145;
t519 = t143 / 0.2e1;
t602 = -mrSges(5,1) - mrSges(4,3);
t574 = mrSges(6,2) - mrSges(7,3);
t315 = qJ(2) + qJ(3);
t308 = cos(t315);
t600 = t308 * t530;
t598 = -Ifges(4,6) + Ifges(5,5);
t597 = t143 * t572 + t573 * t80 + t606 * t81;
t196 = Ifges(6,4) * t356;
t227 = t242 * qJD(1);
t222 = qJD(5) + t227;
t479 = Ifges(7,5) * t356;
t567 = t199 * t573 + t222 * t572 + t196 - t479;
t362 = -pkin(5) * t320 - qJ(6) * t316;
t353 = -pkin(4) + t362;
t431 = qJD(5) * t320;
t432 = qJD(5) * t316;
t424 = t502 * pkin(7);
t263 = pkin(8) * t502 + t424;
t245 = t263 * qJD(1);
t446 = t317 * t245;
t596 = pkin(5) * t431 + qJ(6) * t432 - t320 * qJD(6) - t227 * t353 + qJD(4) + t446;
t322 = -pkin(8) - pkin(7);
t262 = t322 * t318;
t244 = qJD(1) * t262;
t411 = t501 * t244;
t187 = t411 - t446;
t399 = qJD(3) * t501;
t386 = pkin(2) * t399;
t279 = t386 + qJD(4);
t595 = -t279 + t187;
t307 = sin(t315);
t594 = t604 * t307;
t593 = t316 * t572 + t320 * t605;
t477 = Ifges(7,5) * t320;
t480 = Ifges(6,4) * t320;
t592 = t316 * t573 - t477 + t480;
t241 = t317 * t318 - t383;
t410 = t501 * t245;
t186 = t244 * t317 + t410;
t434 = qJD(3) * t317;
t421 = pkin(2) * t434;
t590 = t421 - t186;
t462 = t227 * t320;
t589 = -t431 - t462;
t463 = t227 * t316;
t588 = t432 + t463;
t237 = qJD(2) * pkin(2) + t244;
t212 = t501 * t237;
t183 = -t212 + t446;
t491 = t227 * pkin(4);
t349 = t183 + t491;
t584 = t349 + qJD(4);
t103 = -t314 * t520 + t584;
t311 = t502 * pkin(2);
t295 = t311 + pkin(1);
t261 = t295 * qJD(1);
t333 = -t227 * qJ(4) - t261;
t107 = t226 * t520 + t333;
t466 = qJDD(1) * pkin(1);
t218 = -pkin(2) * t591 - t466;
t325 = t145 * qJ(4) - t227 * qJD(4) + t218;
t22 = t146 * t520 + t325;
t240 = t247 * pkin(7);
t194 = qJDD(2) * pkin(2) - t247 * pkin(8) - t240;
t239 = t591 * pkin(7);
t197 = pkin(8) * t591 + t239;
t58 = t194 * t501 - t317 * t197 - t237 * t434 - t245 * t399;
t335 = qJDD(4) - t58;
t27 = -t145 * pkin(4) - t313 * t520 + t335;
t4 = t103 * t431 - t107 * t432 + t320 * t22 + t316 * t27;
t1 = qJ(6) * t143 + qJD(6) * t222 + t4;
t38 = t103 * t316 + t107 * t320;
t5 = -qJD(5) * t38 - t22 * t316 + t27 * t320;
t3 = -pkin(5) * t143 + qJDD(6) - t5;
t587 = -t1 * t316 + t3 * t320;
t586 = -t608 * t308 + (mrSges(4,2) - mrSges(5,3)) * t307;
t220 = Ifges(4,4) * t226;
t583 = t227 * Ifges(4,1) + t314 * Ifges(4,5) + t199 * t572 + t222 * t570 + t356 * t605 - t220;
t582 = Ifges(6,4) * t523 + Ifges(6,6) * t519 - t80 * Ifges(7,5) / 0.2e1 - t143 * Ifges(7,6) / 0.2e1 + (Ifges(6,2) + Ifges(7,3)) * t522;
t184 = t317 * t237 + t410;
t498 = pkin(4) * t226;
t138 = t184 - t498;
t306 = t314 * qJ(4);
t119 = t138 + t306;
t370 = t320 * mrSges(7,1) + t316 * mrSges(7,3);
t371 = mrSges(6,1) * t320 - mrSges(6,2) * t316;
t48 = -pkin(5) * t356 - qJ(6) * t199 + t119;
t472 = t199 * Ifges(6,4);
t97 = Ifges(6,2) * t356 + t222 * Ifges(6,6) + t472;
t581 = t119 * t371 + t370 * t48 - t320 * t97 / 0.2e1;
t154 = t226 * pkin(3) + t333;
t172 = -t306 - t184;
t580 = -t261 * mrSges(4,1) + t172 * mrSges(5,1) - t154 * mrSges(5,2) - t184 * mrSges(4,3);
t37 = t103 * t320 - t107 * t316;
t359 = t316 * t37 - t320 * t38;
t542 = t316 * t4 + t320 * t5;
t328 = -qJD(5) * t359 + t542;
t29 = -pkin(5) * t222 + qJD(6) - t37;
t30 = qJ(6) * t222 + t38;
t360 = t29 * t316 + t30 * t320;
t329 = qJD(5) * t360 - t587;
t32 = mrSges(6,1) * t143 - mrSges(6,3) * t80;
t33 = -t143 * mrSges(7,1) + t80 * mrSges(7,2);
t568 = t32 - t33;
t31 = -mrSges(7,2) * t81 + mrSges(7,3) * t143;
t34 = -mrSges(6,2) * t143 - mrSges(6,3) * t81;
t569 = t31 + t34;
t579 = m(6) * t328 + m(7) * t329 + t316 * t569 + t320 * t568;
t540 = -qJD(4) - t183;
t156 = -pkin(3) * t314 - t540;
t578 = t156 * mrSges(5,1) + t37 * mrSges(6,1) - t29 * mrSges(7,1) - t261 * mrSges(4,2) - t38 * mrSges(6,2) + t183 * mrSges(4,3) - t154 * mrSges(5,3) + t30 * mrSges(7,3);
t576 = t591 / 0.2e1;
t565 = t386 - t411 + t596;
t564 = -t212 + t596;
t122 = mrSges(5,1) * t146 - mrSges(5,3) * t313;
t24 = mrSges(6,1) * t81 + mrSges(6,2) * t80;
t563 = t24 - t122;
t118 = -mrSges(6,1) * t356 + mrSges(6,2) * t199;
t205 = mrSges(5,1) * t226 - mrSges(5,3) * t314;
t559 = t118 - t205;
t379 = -qJ(4) * t242 - t295;
t144 = t241 * t520 + t379;
t200 = -t501 * t262 + t263 * t317;
t170 = pkin(4) * t242 + t200;
t558 = t320 * t144 + t316 * t170;
t486 = mrSges(7,2) * t356;
t149 = mrSges(7,3) * t222 + t486;
t484 = mrSges(6,3) * t356;
t150 = -mrSges(6,2) * t222 + t484;
t557 = t149 + t150;
t483 = mrSges(6,3) * t199;
t151 = mrSges(6,1) * t222 - t483;
t485 = mrSges(7,2) * t199;
t152 = -mrSges(7,1) * t222 + t485;
t556 = -t151 + t152;
t203 = -mrSges(4,2) * t314 - mrSges(4,3) * t226;
t555 = t203 - t205;
t554 = t602 * t227 + t314 * t608;
t553 = t491 - t595;
t451 = t308 * t320;
t534 = -m(7) * qJ(6) - mrSges(7,3);
t550 = -t451 * t534 + t594;
t423 = t501 * pkin(2);
t294 = -t423 - pkin(3);
t283 = -pkin(9) + t294;
t549 = t283 * t431 + t316 * t421;
t548 = t283 * t432 - t320 * t421;
t547 = t502 * t239 + t240 * t318;
t546 = t143 * t570 + t572 * t80 - t605 * t81;
t401 = qJD(1) * t502;
t496 = pkin(7) * t318;
t545 = (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t401) * t496 + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t436) * t424;
t544 = t316 * pkin(5) - qJ(6) * t320;
t319 = sin(qJ(1));
t321 = cos(qJ(1));
t541 = g(1) * t321 + g(2) * t319;
t539 = m(5) + t610;
t348 = mrSges(3,1) * t318 + mrSges(3,2) * t502;
t482 = Ifges(3,4) * t318;
t538 = t318 * (Ifges(3,1) * t502 - t482) / 0.2e1 - pkin(1) * t348;
t450 = t308 * t321;
t266 = qJ(4) * t450;
t422 = mrSges(6,2) * t451;
t448 = t316 * t321;
t453 = t307 * t321;
t536 = -mrSges(5,2) * t453 - mrSges(5,3) * t450 - t266 * t610 - t321 * t422 + t448 * t600;
t452 = t308 * t319;
t264 = qJ(4) * t452;
t449 = t316 * t319;
t455 = t307 * t319;
t535 = -mrSges(5,2) * t455 - mrSges(5,3) * t452 - t264 * t610 - t319 * t422 + t449 * t600;
t531 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) + t602;
t191 = qJD(2) * t242 + t332;
t190 = t241 * t314;
t435 = qJD(2) * t318;
t301 = pkin(2) * t435;
t343 = qJ(4) * t190 - qJD(4) * t242 + t301;
t45 = t191 * t520 + t343;
t201 = t317 * t262 + t263 * t501;
t246 = qJD(2) * t262;
t334 = qJD(2) * t263;
t109 = qJD(3) * t201 + t317 * t246 + t501 * t334;
t74 = -t190 * pkin(4) + t109;
t11 = -qJD(5) * t558 - t316 * t45 + t320 * t74;
t260 = -mrSges(3,1) * t502 + t318 * mrSges(3,2);
t529 = -m(3) * pkin(1) - mrSges(2,1) + t260 + t586;
t454 = t307 * t320;
t456 = t307 * t316;
t528 = t308 * t601 - t454 * t574 - t456 * t575 + t586;
t527 = mrSges(6,2) + t534;
t499 = pkin(2) * t318;
t526 = -m(5) * (-pkin(3) * t307 - t499) - m(7) * (-qJ(6) * t451 - t499) + mrSges(7,3) * t451 + m(6) * t499 + t594;
t525 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t521 = t81 / 0.2e1;
t517 = t356 / 0.2e1;
t516 = -t356 / 0.2e1;
t515 = -t199 / 0.2e1;
t514 = t199 / 0.2e1;
t513 = -t222 / 0.2e1;
t512 = t222 / 0.2e1;
t511 = -t226 / 0.2e1;
t510 = t226 / 0.2e1;
t509 = -t227 / 0.2e1;
t508 = t227 / 0.2e1;
t505 = -t314 / 0.2e1;
t504 = t314 / 0.2e1;
t500 = pkin(2) * t317;
t497 = pkin(5) * t226;
t493 = g(3) * t308;
t291 = t308 * pkin(3);
t481 = Ifges(6,4) * t316;
t478 = Ifges(7,5) * t316;
t471 = t227 * Ifges(4,4);
t470 = t227 * Ifges(5,6);
t465 = t191 * t316;
t464 = t191 * t320;
t460 = t241 * t316;
t459 = t241 * t320;
t284 = t307 * qJ(4);
t444 = t319 * t320;
t443 = t320 * t321;
t178 = pkin(3) * t227 + qJ(4) * t226;
t299 = pkin(2) * t436;
t157 = t178 + t299;
t221 = t227 * pkin(9);
t115 = t157 + t221;
t147 = t186 - t498;
t47 = t320 * t115 + t316 * t147;
t130 = t178 + t221;
t51 = t320 * t130 + t316 * t138;
t437 = t291 + t284;
t430 = qJD(5) * t520;
t419 = Ifges(3,4) * t502;
t290 = t308 * pkin(9);
t414 = t290 + t437;
t413 = t311 + t437;
t406 = t316 * t430;
t405 = t320 * t430;
t397 = -t432 / 0.2e1;
t396 = t431 / 0.2e1;
t123 = -t145 * mrSges(5,1) + t313 * mrSges(5,2);
t391 = t321 * t295 - t319 * t322;
t256 = qJ(4) + t544;
t372 = mrSges(4,1) * t307 + mrSges(4,2) * t308;
t366 = Ifges(6,2) * t320 + t481;
t363 = -Ifges(7,3) * t320 + t478;
t46 = -t115 * t316 + t147 * t320;
t50 = -t130 * t316 + t138 * t320;
t63 = -t144 * t316 + t170 * t320;
t352 = pkin(3) * t450 + qJ(4) * t453 + t391;
t347 = t502 * Ifges(3,2) + t482;
t346 = Ifges(3,5) * t502 - Ifges(3,6) * t318;
t345 = t241 * t431 + t465;
t344 = t241 * t432 - t464;
t10 = -t144 * t432 + t170 * t431 + t316 * t74 + t320 * t45;
t57 = t317 * t194 + t501 * t197 + t237 * t399 - t245 * t434;
t108 = -t501 * t246 - t262 * t399 + t263 * t434 + t317 * t334;
t340 = pkin(5) * t456 - qJ(6) * t454 + t414;
t49 = -t313 * qJ(4) - t314 * qJD(4) - t57;
t28 = -pkin(4) * t146 - t49;
t164 = t314 * Ifges(5,5) + t226 * Ifges(5,3) - t470;
t219 = Ifges(5,6) * t226;
t165 = t314 * Ifges(5,4) - t227 * Ifges(5,2) + t219;
t166 = -t226 * Ifges(4,2) + t314 * Ifges(4,6) + t471;
t55 = -t313 * pkin(3) + t335;
t7 = pkin(5) * t81 - qJ(6) * t80 - qJD(6) * t199 + t28;
t195 = Ifges(7,5) * t199;
t94 = t222 * Ifges(7,6) - Ifges(7,3) * t356 + t195;
t324 = (-t471 + t164) * t509 + (t470 + t166) * t508 - (t199 * t592 + t222 * t593) * qJD(5) / 0.2e1 + (t219 + t165) * t511 + (-t366 / 0.2e1 + t363 / 0.2e1) * t585 + (-Ifges(4,1) * t509 + Ifges(5,2) * t508 - Ifges(6,6) * t516 - Ifges(7,6) * t517 - t505 * t609 - t513 * t570 - t515 * t572 + t578) * t226 - t609 * t145 + (t37 * t588 + t38 * t589 - t542) * mrSges(6,3) + (-t29 * t588 + t30 * t589 + t587) * mrSges(7,2) + (-t220 + t583) * t510 + (t28 * mrSges(6,1) + t7 * mrSges(7,1) - Ifges(6,2) * t522 + Ifges(7,3) * t521 - t519 * t605 - t582) * t316 + (Ifges(5,1) + Ifges(4,3)) * t313 + t581 * qJD(5) + t477 * t521 + t480 * t522 + (t462 / 0.2e1 + t396) * t94 + (t28 * mrSges(6,2) - t7 * mrSges(7,3) + t572 * t519 + t573 * t523 + t597 / 0.2e1) * t320 + (t478 - t481) * t523 + t598 * t146 + (-Ifges(4,2) * t510 + Ifges(5,3) * t511 + t363 * t517 + t366 * t516 + t505 * t598 + t513 * t593 + t515 * t592 - t580 + t581) * t227 - t49 * mrSges(5,3) + t55 * mrSges(5,2) - t57 * mrSges(4,2) + t58 * mrSges(4,1) + (t397 - t463 / 0.2e1) * t567;
t297 = Ifges(3,4) * t401;
t289 = qJ(4) + t500;
t236 = t256 + t500;
t225 = Ifges(3,1) * t436 + Ifges(3,5) * qJD(2) + t297;
t224 = Ifges(3,6) * qJD(2) + qJD(1) * t347;
t217 = -t307 * t449 + t443;
t216 = t307 * t444 + t448;
t215 = t307 * t448 + t444;
t214 = -t307 * t443 + t449;
t213 = t226 * qJ(6);
t182 = pkin(3) * t241 + t379;
t181 = -mrSges(5,2) * t226 - mrSges(5,3) * t227;
t180 = mrSges(4,1) * t226 + mrSges(4,2) * t227;
t171 = -t241 * pkin(4) + t201;
t121 = -mrSges(4,2) * t313 - mrSges(4,3) * t146;
t120 = mrSges(4,1) * t313 + mrSges(4,3) * t145;
t117 = -mrSges(7,1) * t356 - mrSges(7,3) * t199;
t116 = pkin(5) * t199 - qJ(6) * t356;
t100 = t241 * t353 + t201;
t87 = pkin(3) * t191 + t343;
t73 = -pkin(4) * t191 - t108;
t53 = -pkin(5) * t242 - t63;
t52 = qJ(6) * t242 + t558;
t42 = -t50 + t497;
t41 = -t213 + t51;
t40 = -t46 + t497;
t39 = -t213 + t47;
t36 = t146 * pkin(3) + t325;
t23 = mrSges(7,1) * t81 - mrSges(7,3) * t80;
t20 = (qJD(5) * t544 - qJD(6) * t316) * t241 + t353 * t191 - t108;
t9 = pkin(5) * t190 - t11;
t8 = -qJ(6) * t190 + qJD(6) * t242 + t10;
t2 = [(t218 * mrSges(4,1) + t49 * mrSges(5,1) - t36 * mrSges(5,2) - t57 * mrSges(4,3) - t28 * t371 + t363 * t521 + t366 * t522 - t7 * t370 + t97 * t397 + t592 * t523 + t593 * t519 + (Ifges(4,2) + Ifges(5,3)) * t146 + t607 * t145 + t598 * t313 + t567 * t396) * t241 + (-m(4) * t391 - m(5) * t352 + t601 * t450 - t610 * (t319 * pkin(4) + pkin(9) * t450 + t352) + t530 * t215 + t527 * t214 + t529 * t321 + t531 * t319) * g(2) + (t530 * t217 + t527 * t216 + (m(4) * t295 + m(5) * t291 - t539 * (-t295 - t284) + t604 * t308 - t529) * t319 + (-t610 * (pkin(4) - t322) + (m(4) + m(5)) * t322 + t531) * t321) * g(1) + t502 * (Ifges(3,4) * t247 + Ifges(3,2) * t591) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t591 + t247 * mrSges(3,2)) + (t247 * t496 + t424 * t591 + t547) * mrSges(3,3) + (-t583 / 0.2e1 - t578 - t570 * t512 - t572 * t514 - Ifges(4,1) * t508 + Ifges(5,2) * t509 + Ifges(5,6) * t510 - Ifges(4,4) * t511 - Ifges(7,6) * t516 - Ifges(6,6) * t517 - t609 * t504 + t165 / 0.2e1) * t190 + (t344 * t606 + t573 * t345) * t514 + t567 * t465 / 0.2e1 + (m(4) * t57 - m(5) * t49 + t121 - t122) * t201 + (-m(4) * t58 + m(5) * t55 - t120 + t123) * t200 + t347 * t576 + t582 * t459 + (-t344 * t605 + t345 * t572) * t512 + (-mrSges(3,1) * t496 - mrSges(3,2) * t424 + Ifges(3,5) * t318 + Ifges(3,6) * t502) * qJDD(2) + (Ifges(7,5) * t345 + Ifges(7,3) * t344) * t516 + (Ifges(6,4) * t345 - Ifges(6,2) * t344) * t517 + (t97 / 0.2e1 - t94 / 0.2e1) * t464 + (m(4) * t183 + m(5) * t156 - t554) * t109 + (-m(4) * t184 + m(5) * t172 - t555) * t108 - t260 * t466 - t224 * t435 / 0.2e1 + t247 * t419 / 0.2e1 + (t1 * t459 + t29 * t345 + t3 * t460 - t30 * t344) * mrSges(7,2) + m(6) * (t10 * t38 + t11 * t37 + t119 * t73 + t171 * t28 + t4 * t558 + t5 * t63) + t558 * t34 + m(7) * (t1 * t52 + t100 * t7 + t20 * t48 + t29 * t9 + t3 * t53 + t30 * t8) + m(5) * (t154 * t87 + t182 * t36) + t180 * t301 + (-t344 * t38 - t345 * t37 + t4 * t459 - t460 * t5) * mrSges(6,3) + (t346 * qJD(2) / 0.2e1 - t545) * qJD(2) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t547) + m(4) * (-t218 * t295 - t261 * t301) + (qJD(1) * (-Ifges(3,2) * t318 + t419) + t225) * t400 / 0.2e1 + Ifges(2,3) * qJDD(1) - t295 * (mrSges(4,1) * t146 - mrSges(4,2) * t145) + t87 * t181 + t182 * (-mrSges(5,2) * t146 + mrSges(5,3) * t145) + t171 * t24 + t9 * t152 + t8 * t149 + t10 * t150 + t11 * t151 + t20 * t117 + t73 * t118 + (qJD(5) * t94 + t597) * t460 / 0.2e1 + (-t166 / 0.2e1 + t164 / 0.2e1 - Ifges(4,4) * t508 + Ifges(5,6) * t509 + Ifges(5,3) * t510 - Ifges(4,2) * t511 + t598 * t504 + t580) * t191 + (t55 * mrSges(5,1) + t218 * mrSges(4,2) - t58 * mrSges(4,3) - t36 * mrSges(5,3) + Ifges(6,6) * t522 + Ifges(7,6) * t521 + t570 * t519 + t572 * t523 + t525 + t546 / 0.2e1 - t607 * t146 + (-Ifges(5,2) - Ifges(4,1)) * t145 + t609 * t313) * t242 + (Ifges(3,1) * t247 + Ifges(3,4) * t576) * t318 + t52 * t31 + t53 * t33 + t63 * t32 + t100 * t23 + t119 * (mrSges(6,1) * t344 + mrSges(6,2) * t345) + t48 * (mrSges(7,1) * t344 - mrSges(7,3) * t345) + t538 * t428; t579 * t283 + Ifges(3,6) * t591 + (-m(4) * t311 - m(6) * (t290 + t413) - m(5) * t413 + t260 - m(7) * (t311 + t340) + t528) * g(3) - t590 * t554 + t563 * t289 + t565 * t117 + (-t47 + t549) * t150 + (-t39 + t549) * t149 + t553 * t118 - t555 * t187 + t224 * t436 / 0.2e1 + (m(4) * t499 + t348 + t372) * t541 - t346 * t428 / 0.2e1 - t180 * t299 + t324 + (-qJD(1) * t538 + t545) * qJD(1) + (-t29 * t40 - t30 * t39 + t236 * t7 + (-t29 * t320 + t30 * t316) * t421 + t565 * t48) * m(7) + (-t37 * t46 - t38 * t47 + t28 * t289 + (t316 * t38 + t320 * t37) * t421 + t553 * t119) * m(6) + (-t154 * t157 + t156 * t590 + t172 * t595 - t289 * t49 + t294 * t55) * m(5) + (-t46 - t548) * t151 + (-t40 + t548) * t152 + (-t183 * t186 - t184 * t187 + t261 * t299 + (t501 * t58 + t317 * t57 + (t183 * t317 + t184 * t501) * qJD(3)) * pkin(2)) * m(4) - (-Ifges(3,2) * t436 + t225 + t297) * t401 / 0.2e1 + Ifges(3,3) * qJDD(2) + t294 * t123 - t279 * t205 + Ifges(3,5) * t247 - t239 * mrSges(3,2) - t240 * mrSges(3,1) + t236 * t23 - t157 * t181 + t203 * t386 + t121 * t500 + (-m(5) * t264 + t319 * t526 + t535) * g(2) + (-m(5) * t266 + t321 * t526 + t536) * g(1) + t120 * t423; -pkin(3) * t123 + t349 * t118 - t178 * t181 + t256 * t23 + t324 + t541 * t372 + t554 * t184 + t555 * t183 + (-t406 - t42) * t152 + (t406 - t50) * t151 + (-t405 - t51) * t150 + (-t405 - t41) * t149 + t564 * t117 + t559 * qJD(4) + t563 * qJ(4) + (t7 * t256 - t29 * t42 - t30 * t41 + t48 * t564) * m(7) + (t28 * qJ(4) + t119 * t584 - t37 * t50 - t38 * t51) * m(6) + (-pkin(3) * t55 - qJ(4) * t49 - t154 * t178 - t156 * t184 + t172 * t540) * m(5) + (-m(5) * (-pkin(3) * t455 + t264) + t550 * t319 + t535) * g(2) + (-m(5) * (-pkin(3) * t453 + t266) + t550 * t321 + t536) * g(1) + (-m(5) * t437 - m(6) * t414 - m(7) * t340 + t528) * g(3) - t579 * t520; t227 * t181 + (-t117 - t559) * t314 + (t222 * t557 + t568) * t320 + (t222 * t556 + t569) * t316 + t123 + (t227 * t360 - t314 * t48 + t329) * m(7) + (-t119 * t314 - t227 * t359 + t328) * m(6) + (t154 * t227 + t172 * t314 + t55) * m(5) + (-t307 * t541 + t493) * t539; t525 + (t370 + t371) * t493 + (-Ifges(6,2) * t199 + t196 + t567) * t516 + (-m(7) * t29 + t483 - t556) * t38 + t546 + (Ifges(7,3) * t199 + t479) * t517 + (-m(7) * t30 + t484 - t557) * t37 + (-t199 * t605 + t356 * t572) * t513 + (t356 * t573 + t195 - t472 + t94) * t515 - t48 * (mrSges(7,1) * t199 - mrSges(7,3) * t356) - t119 * (mrSges(6,1) * t199 + mrSges(6,2) * t356) + (-t216 * t575 - t217 * t574) * g(2) + qJD(6) * t149 - t116 * t117 + (-t362 * t493 - t116 * t48 - pkin(5) * t3 + qJ(6) * t1 + qJD(6) * t30 - g(1) * (-pkin(5) * t214 + qJ(6) * t215) - g(2) * (pkin(5) * t216 - qJ(6) * t217)) * m(7) + qJ(6) * t31 + (t214 * t575 + t215 * t574) * g(1) + t97 * t514 - pkin(5) * t33 + t30 * t485 - t29 * t486; t199 * t117 - t222 * t149 + (-g(1) * t214 + g(2) * t216 - g(3) * t451 + t48 * t199 - t30 * t222 + t3) * m(7) + t33;];
tau  = t2;
