% Calculate vector of inverse dynamics joint torques for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:49
% EndTime: 2019-12-31 22:39:55
% DurationCPUTime: 34.89s
% Computational Cost: add. (15827->887), mult. (38368->1221), div. (0->0), fcn. (30262->14), ass. (0->398)
t322 = cos(pkin(5));
t326 = sin(qJ(2));
t313 = t322 * t326 * pkin(1);
t325 = sin(qJ(3));
t330 = cos(qJ(3));
t370 = pkin(3) * t325 - pkin(9) * t330;
t321 = sin(pkin(5));
t331 = cos(qJ(2));
t433 = t321 * t331;
t585 = -(t313 + (pkin(7) + t370) * t433) * qJD(1) + t370 * qJD(3);
t324 = sin(qJ(4));
t329 = cos(qJ(4));
t419 = qJD(1) * t321;
t427 = t330 * t331;
t219 = (-t324 * t427 + t326 * t329) * t419;
t414 = qJD(3) * t330;
t584 = t324 * t414 + t219;
t411 = qJD(4) * t329;
t542 = t325 * t411 + t584;
t418 = qJD(1) * t322;
t305 = qJD(2) + t418;
t391 = t326 * t419;
t374 = t330 * t391;
t229 = t305 * t325 + t374;
t450 = t229 * Ifges(4,4);
t390 = t331 * t419;
t289 = qJD(3) - t390;
t551 = t289 * Ifges(4,6);
t228 = t305 * t330 - t325 * t391;
t575 = t228 * Ifges(4,2);
t132 = t450 + t551 + t575;
t223 = qJD(4) - t228;
t215 = qJD(5) + t223;
t487 = t215 / 0.2e1;
t183 = -t229 * t324 + t289 * t329;
t184 = t229 * t329 + t289 * t324;
t323 = sin(qJ(5));
t328 = cos(qJ(5));
t97 = t183 * t323 + t184 * t328;
t499 = t97 / 0.2e1;
t380 = t328 * t183 - t184 * t323;
t501 = t380 / 0.2e1;
t404 = pkin(1) * t418;
t248 = -pkin(7) * t391 + t331 * t404;
t211 = -pkin(2) * t305 - t248;
t111 = -pkin(3) * t228 - pkin(9) * t229 + t211;
t421 = pkin(7) * t433 + t313;
t251 = t421 * qJD(1);
t212 = pkin(8) * t305 + t251;
t221 = (-pkin(2) * t331 - pkin(8) * t326 - pkin(1)) * t419;
t140 = t212 * t330 + t221 * t325;
t114 = pkin(9) * t289 + t140;
t63 = t329 * t111 - t114 * t324;
t64 = t111 * t324 + t114 * t329;
t523 = -t63 * mrSges(5,1) + t64 * mrSges(5,2);
t553 = t184 * Ifges(5,5) + t97 * Ifges(6,5) + t183 * Ifges(5,6) + Ifges(6,6) * t380 + t223 * Ifges(5,3) + t215 * Ifges(6,3);
t583 = Ifges(6,5) * t499 + Ifges(6,6) * t501 + Ifges(6,3) * t487 - t132 / 0.2e1 + t553 / 0.2e1 - t523;
t348 = (pkin(2) * t326 - pkin(8) * t331) * t321;
t249 = qJD(1) * t348;
t179 = t330 * t248 + t325 * t249;
t159 = pkin(9) * t391 + t179;
t285 = -pkin(3) * t330 - pkin(9) * t325 - pkin(2);
t413 = qJD(4) * t324;
t415 = qJD(3) * t325;
t548 = -t329 * t159 + t285 * t411 + (-t329 * t415 - t330 * t413) * pkin(8) + t585 * t324;
t402 = pkin(8) * t415;
t582 = t585 * t329 + (t159 + t402) * t324;
t409 = qJD(1) * qJD(2);
t256 = (-qJDD(1) * t331 + t326 * t409) * t321;
t245 = qJDD(3) + t256;
t408 = qJDD(1) * t321;
t560 = pkin(7) * t408 + qJD(2) * t404;
t407 = qJDD(1) * t322;
t561 = -pkin(7) * t321 * t409 + pkin(1) * t407;
t185 = t326 * t561 + t331 * t560;
t304 = qJDD(2) + t407;
t167 = pkin(8) * t304 + t185;
t257 = (qJDD(1) * t326 + t331 * t409) * t321;
t175 = -pkin(1) * t408 + pkin(2) * t256 - pkin(8) * t257;
t61 = t330 * t167 + t325 * t175 - t212 * t415 + t221 * t414;
t52 = pkin(9) * t245 + t61;
t146 = qJD(3) * t228 + t257 * t330 + t304 * t325;
t147 = -qJD(3) * t374 - t325 * t257 + t304 * t330 - t305 * t415;
t186 = -t326 * t560 + t331 * t561;
t168 = -pkin(2) * t304 - t186;
t59 = -pkin(3) * t147 - pkin(9) * t146 + t168;
t13 = -qJD(4) * t64 - t324 * t52 + t329 * t59;
t556 = t13 * mrSges(5,1);
t12 = t111 * t411 - t114 * t413 + t324 * t59 + t329 * t52;
t557 = t12 * mrSges(5,2);
t581 = t556 - t557;
t367 = mrSges(4,1) * t330 - mrSges(4,2) * t325;
t332 = -pkin(10) - pkin(9);
t536 = -m(5) * pkin(9) + m(6) * t332 - mrSges(5,3) - mrSges(6,3);
t316 = pkin(4) * t329 + pkin(3);
t320 = qJ(4) + qJ(5);
t317 = sin(t320);
t318 = cos(t320);
t365 = -mrSges(5,1) * t329 + mrSges(5,2) * t324;
t538 = m(5) * pkin(3) + m(6) * t316 + mrSges(6,1) * t318 - mrSges(6,2) * t317 - t365;
t580 = t325 * t536 - t330 * t538 - t367;
t47 = -pkin(10) * t184 + t63;
t41 = pkin(4) * t223 + t47;
t48 = pkin(10) * t183 + t64;
t445 = t323 * t48;
t15 = t328 * t41 - t445;
t442 = t328 * t48;
t16 = t323 * t41 + t442;
t452 = t140 * mrSges(4,3);
t579 = t452 - t211 * mrSges(4,1) - t15 * mrSges(6,1) + t16 * mrSges(6,2) + t551 / 0.2e1;
t139 = -t325 * t212 + t221 * t330;
t453 = t139 * mrSges(4,3);
t552 = t289 * Ifges(4,5);
t578 = t211 * mrSges(4,2) - t453 + t552 / 0.2e1;
t577 = -m(4) - m(5);
t222 = Ifges(4,4) * t228;
t220 = (t324 * t326 + t329 * t427) * t419;
t428 = t329 * t330;
t314 = pkin(8) * t428;
t375 = t325 * t390;
t574 = -pkin(4) * t375 + pkin(10) * t220 + (pkin(4) * t325 - pkin(10) * t428) * qJD(3) + (-t314 + (pkin(10) * t325 - t285) * t324) * qJD(4) + t582;
t573 = -pkin(10) * t542 + t548;
t395 = qJD(4) * t332;
t438 = t228 * t324;
t172 = pkin(3) * t229 - pkin(9) * t228;
t81 = t329 * t139 + t324 * t172;
t572 = pkin(10) * t438 + t324 * t395 - t81;
t437 = t228 * t329;
t80 = -t139 * t324 + t329 * t172;
t571 = -pkin(4) * t229 + pkin(10) * t437 + t329 * t395 - t80;
t568 = t413 - t438;
t138 = qJDD(4) - t147;
t73 = qJD(4) * t183 + t146 * t329 + t245 * t324;
t10 = pkin(4) * t138 - pkin(10) * t73 + t13;
t74 = -qJD(4) * t184 - t146 * t324 + t245 * t329;
t11 = pkin(10) * t74 + t12;
t2 = qJD(5) * t15 + t10 * t323 + t11 * t328;
t3 = -qJD(5) * t16 + t10 * t328 - t11 * t323;
t567 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t481 = t245 / 0.2e1;
t494 = t147 / 0.2e1;
t495 = t146 / 0.2e1;
t509 = Ifges(4,1) * t495 + Ifges(4,4) * t494 + Ifges(4,5) * t481;
t127 = qJDD(5) + t138;
t498 = t127 / 0.2e1;
t25 = -qJD(5) * t97 - t323 * t73 + t328 * t74;
t516 = t25 / 0.2e1;
t24 = qJD(5) * t380 + t323 * t74 + t328 * t73;
t517 = t24 / 0.2e1;
t519 = Ifges(6,1) * t517 + Ifges(6,4) * t516 + Ifges(6,5) * t498;
t520 = Ifges(6,4) * t517 + Ifges(6,2) * t516 + Ifges(6,6) * t498;
t28 = Ifges(5,5) * t73 + Ifges(5,6) * t74 + Ifges(5,3) * t138;
t6 = Ifges(6,5) * t24 + Ifges(6,6) * t25 + Ifges(6,3) * t127;
t563 = t28 + t6;
t363 = t317 * mrSges(6,1) + t318 * mrSges(6,2);
t441 = t329 * mrSges(5,2);
t364 = t324 * mrSges(5,1) + t441;
t562 = -t363 - t364;
t178 = -t325 * t248 + t249 * t330;
t158 = -pkin(3) * t391 - t178;
t401 = pkin(8) * t414;
t559 = t401 - t158;
t410 = m(6) - t577;
t558 = pkin(2) * t410 + mrSges(3,1) - t580;
t508 = t73 / 0.2e1;
t507 = t74 / 0.2e1;
t496 = t138 / 0.2e1;
t465 = mrSges(3,2) - mrSges(4,3);
t275 = t329 * t285;
t431 = t325 * t329;
t191 = -pkin(10) * t431 + t275 + (-pkin(8) * t324 - pkin(4)) * t330;
t239 = t324 * t285 + t314;
t432 = t324 * t325;
t208 = -pkin(10) * t432 + t239;
t116 = t191 * t323 + t208 * t328;
t555 = -qJD(5) * t116 - t323 * t573 + t328 * t574;
t115 = t191 * t328 - t208 * t323;
t554 = qJD(5) * t115 + t323 * t574 + t328 * t573;
t291 = t332 * t324;
t292 = t332 * t329;
t213 = t291 * t328 + t292 * t323;
t550 = qJD(5) * t213 + t323 * t571 + t328 * t572;
t214 = t291 * t323 - t292 * t328;
t549 = -qJD(5) * t214 - t323 * t572 + t328 * t571;
t547 = -qJD(4) * t239 + t582;
t546 = pkin(4) * t568 - t140;
t518 = m(6) * pkin(4);
t545 = -t518 - mrSges(5,1);
t352 = t323 * t324 - t328 * t329;
t255 = t352 * t325;
t277 = t323 * t329 + t324 * t328;
t533 = qJD(4) + qJD(5);
t194 = t533 * t277;
t124 = -t194 * t325 - t352 * t414;
t151 = t219 * t323 + t220 * t328;
t544 = t124 - t151;
t125 = t255 * t533 - t277 * t414;
t150 = t219 * t328 - t220 * t323;
t543 = t125 - t150;
t436 = t321 * t326;
t306 = pkin(7) * t436;
t473 = pkin(1) * t331;
t240 = t306 + (-pkin(2) - t473) * t322;
t262 = -t322 * t330 + t325 * t436;
t434 = t321 * t330;
t263 = t322 * t325 + t326 * t434;
t155 = pkin(3) * t262 - pkin(9) * t263 + t240;
t241 = pkin(8) * t322 + t421;
t422 = pkin(2) * t433 + pkin(8) * t436;
t474 = pkin(1) * t321;
t242 = -t422 - t474;
t174 = t330 * t241 + t325 * t242;
t157 = -pkin(9) * t433 + t174;
t77 = t324 * t155 + t329 * t157;
t412 = qJD(4) * t325;
t541 = t324 * t412 - t329 * t414 + t220;
t540 = pkin(4) * t542 + t559;
t379 = mrSges(3,3) * t391;
t539 = -mrSges(3,1) * t305 - mrSges(4,1) * t228 + mrSges(4,2) * t229 + t379;
t271 = t322 * t473 - t306;
t537 = t375 - t415;
t62 = -t325 * t167 + t175 * t330 - t212 * t414 - t221 * t415;
t535 = -t325 * t62 + t330 * t61;
t534 = t12 * t329 - t13 * t324;
t462 = Ifges(3,4) * t326;
t532 = -t326 * (Ifges(3,1) * t331 - t462) / 0.2e1 + pkin(1) * (mrSges(3,1) * t326 + mrSges(3,2) * t331);
t531 = mrSges(4,1) + t538;
t339 = mrSges(4,2) + t536;
t113 = -pkin(3) * t289 - t139;
t82 = -pkin(4) * t183 + t113;
t530 = -t82 * mrSges(6,1) + mrSges(6,3) * t16;
t529 = t82 * mrSges(6,2) - mrSges(6,3) * t15;
t528 = t410 * pkin(8) - t465;
t525 = t62 * mrSges(4,1) - t61 * mrSges(4,2) + Ifges(4,5) * t146 + Ifges(4,6) * t147 + Ifges(4,3) * t245;
t524 = t186 * mrSges(3,1) - t185 * mrSges(3,2) + Ifges(3,5) * t257 - Ifges(3,6) * t256 + Ifges(3,3) * t304;
t470 = pkin(4) * t324;
t396 = pkin(8) + t470;
t522 = -m(6) * t396 + pkin(8) * t577 + t465 + t562;
t486 = -t223 / 0.2e1;
t488 = -t215 / 0.2e1;
t491 = -t184 / 0.2e1;
t493 = -t183 / 0.2e1;
t500 = -t97 / 0.2e1;
t502 = -t380 / 0.2e1;
t521 = Ifges(5,5) * t491 + Ifges(6,5) * t500 + Ifges(5,6) * t493 + Ifges(6,6) * t502 + Ifges(5,3) * t486 + Ifges(6,3) * t488 + t523;
t29 = t73 * Ifges(5,4) + t74 * Ifges(5,2) + t138 * Ifges(5,6);
t515 = t29 / 0.2e1;
t514 = Ifges(5,1) * t508 + Ifges(5,4) * t507 + Ifges(5,5) * t496;
t475 = Ifges(6,4) * t97;
t45 = Ifges(6,2) * t380 + Ifges(6,6) * t215 + t475;
t513 = -t45 / 0.2e1;
t512 = t45 / 0.2e1;
t93 = Ifges(6,4) * t380;
t46 = Ifges(6,1) * t97 + Ifges(6,5) * t215 + t93;
t511 = -t46 / 0.2e1;
t510 = t46 / 0.2e1;
t459 = Ifges(5,4) * t184;
t84 = Ifges(5,2) * t183 + Ifges(5,6) * t223 + t459;
t506 = -t84 / 0.2e1;
t505 = t84 / 0.2e1;
t181 = Ifges(5,4) * t183;
t85 = Ifges(5,1) * t184 + Ifges(5,5) * t223 + t181;
t504 = -t85 / 0.2e1;
t503 = t85 / 0.2e1;
t492 = t183 / 0.2e1;
t490 = t184 / 0.2e1;
t485 = t223 / 0.2e1;
t484 = -t228 / 0.2e1;
t482 = t229 / 0.2e1;
t478 = cos(qJ(1));
t472 = pkin(4) * t184;
t198 = -t263 * t324 - t329 * t433;
t471 = pkin(4) * t198;
t468 = pkin(8) * t330;
t464 = mrSges(5,3) * t183;
t463 = mrSges(5,3) * t184;
t461 = Ifges(4,4) * t325;
t460 = Ifges(4,4) * t330;
t458 = Ifges(5,4) * t324;
t457 = Ifges(5,4) * t329;
t451 = t228 * Ifges(4,6);
t449 = t229 * Ifges(4,5);
t448 = t289 * Ifges(4,3);
t447 = t305 * Ifges(3,5);
t446 = t305 * Ifges(3,6);
t53 = -pkin(3) * t245 - t62;
t444 = t325 * t53;
t327 = sin(qJ(1));
t435 = t321 * t327;
t430 = t326 * t327;
t429 = t327 * t331;
t393 = t478 * t326;
t265 = t322 * t393 + t429;
t394 = t321 * t478;
t201 = t265 * t330 - t325 * t394;
t392 = t478 * t331;
t264 = -t322 * t392 + t430;
t426 = (-t201 * t317 + t264 * t318) * mrSges(6,1) + (-t201 * t318 - t264 * t317) * mrSges(6,2);
t267 = -t322 * t430 + t392;
t205 = t267 * t330 + t325 * t435;
t266 = t322 * t429 + t393;
t152 = -t205 * t317 + t266 * t318;
t153 = t205 * t318 + t266 * t317;
t425 = t152 * mrSges(6,1) - t153 * mrSges(6,2);
t424 = (-t263 * t317 - t318 * t433) * mrSges(6,1) + (-t263 * t318 + t317 * t433) * mrSges(6,2);
t420 = t478 * pkin(1) + pkin(7) * t435;
t416 = qJD(2) * t321;
t406 = m(6) * t470;
t397 = t267 * pkin(2) + t420;
t389 = t326 * t416;
t388 = t331 * t416;
t383 = t414 / 0.2e1;
t382 = -t412 / 0.2e1;
t381 = -t327 * pkin(1) + pkin(7) * t394;
t76 = t329 * t155 - t157 * t324;
t173 = -t325 * t241 + t242 * t330;
t200 = -t265 * t325 - t330 * t394;
t378 = mrSges(3,3) * t390;
t156 = pkin(3) * t433 - t173;
t368 = mrSges(4,1) * t262 + mrSges(4,2) * t263;
t347 = -t263 * t329 + t324 * t433;
t366 = mrSges(5,1) * t198 + mrSges(5,2) * t347;
t362 = Ifges(4,1) * t330 - t461;
t361 = Ifges(5,1) * t329 - t458;
t360 = Ifges(5,1) * t324 + t457;
t359 = -Ifges(4,2) * t325 + t460;
t358 = -Ifges(5,2) * t324 + t457;
t357 = Ifges(5,2) * t329 + t458;
t356 = Ifges(4,5) * t330 - Ifges(4,6) * t325;
t355 = Ifges(5,5) * t329 - Ifges(5,6) * t324;
t354 = Ifges(5,5) * t324 + Ifges(5,6) * t329;
t58 = pkin(4) * t262 + pkin(10) * t347 + t76;
t65 = pkin(10) * t198 + t77;
t26 = -t323 * t65 + t328 * t58;
t27 = t323 * t58 + t328 * t65;
t119 = t198 * t328 + t323 * t347;
t120 = t198 * t323 - t328 * t347;
t160 = -t205 * t324 + t266 * t329;
t250 = qJD(2) * t348;
t252 = t271 * qJD(2);
t92 = -t241 * t414 - t242 * t415 + t250 * t330 - t325 * t252;
t350 = t567 + t6;
t346 = t113 * t364;
t196 = qJD(3) * t263 + t325 * t388;
t197 = -qJD(3) * t262 + t330 * t388;
t253 = t421 * qJD(2);
t105 = pkin(3) * t196 - pkin(9) * t197 + t253;
t91 = -t241 * t415 + t242 * t414 + t325 * t250 + t330 * t252;
t87 = pkin(9) * t389 + t91;
t34 = t324 * t105 + t155 * t411 - t157 * t413 + t329 * t87;
t88 = -pkin(3) * t389 - t92;
t35 = -qJD(4) * t77 + t329 * t105 - t324 * t87;
t301 = Ifges(3,4) * t390;
t282 = t396 * t325;
t268 = (-mrSges(3,1) * t331 + mrSges(3,2) * t326) * t321;
t254 = t277 * t325;
t247 = -mrSges(3,2) * t305 + t378;
t238 = -t324 * t468 + t275;
t210 = Ifges(3,1) * t391 + t301 + t447;
t209 = t446 + (t331 * Ifges(3,2) + t462) * t419;
t204 = t267 * t325 - t327 * t434;
t188 = mrSges(4,1) * t289 - mrSges(4,3) * t229;
t187 = -mrSges(4,2) * t289 + mrSges(4,3) * t228;
t161 = t205 * t329 + t266 * t324;
t149 = t352 * t228;
t148 = t277 * t228;
t133 = t229 * Ifges(4,1) + t222 + t552;
t131 = t448 + t449 + t451;
t118 = mrSges(5,1) * t223 - t463;
t117 = -mrSges(5,2) * t223 + t464;
t110 = qJD(4) * t198 + t197 * t329 + t324 * t389;
t109 = qJD(4) * t347 - t197 * t324 + t329 * t389;
t104 = -mrSges(4,2) * t245 + mrSges(4,3) * t147;
t103 = mrSges(4,1) * t245 - mrSges(4,3) * t146;
t102 = -mrSges(5,1) * t183 + mrSges(5,2) * t184;
t98 = t156 - t471;
t79 = mrSges(6,1) * t215 - mrSges(6,3) * t97;
t78 = -mrSges(6,2) * t215 + mrSges(6,3) * t380;
t75 = -mrSges(4,1) * t147 + mrSges(4,2) * t146;
t67 = t146 * Ifges(4,4) + t147 * Ifges(4,2) + t245 * Ifges(4,6);
t55 = -mrSges(6,1) * t380 + mrSges(6,2) * t97;
t54 = -pkin(4) * t109 + t88;
t50 = -mrSges(5,2) * t138 + mrSges(5,3) * t74;
t49 = mrSges(5,1) * t138 - mrSges(5,3) * t73;
t40 = -qJD(5) * t120 + t109 * t328 - t110 * t323;
t39 = qJD(5) * t119 + t109 * t323 + t110 * t328;
t38 = -mrSges(5,1) * t74 + mrSges(5,2) * t73;
t33 = -pkin(4) * t74 + t53;
t21 = pkin(10) * t109 + t34;
t20 = t328 * t47 - t445;
t19 = -t323 * t47 - t442;
t18 = -mrSges(6,2) * t127 + mrSges(6,3) * t25;
t17 = mrSges(6,1) * t127 - mrSges(6,3) * t24;
t14 = pkin(4) * t196 - pkin(10) * t110 + t35;
t9 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t5 = -qJD(5) * t27 + t14 * t328 - t21 * t323;
t4 = qJD(5) * t26 + t14 * t323 + t21 * t328;
t1 = [(t119 * t2 - t120 * t3 - t15 * t39 + t16 * t40) * mrSges(6,3) + (-Ifges(5,4) * t347 + Ifges(5,2) * t198) * t507 + (t109 * t64 - t110 * t63 + t12 * t198 + t13 * t347) * mrSges(5,3) + (Ifges(6,5) * t120 + Ifges(6,6) * t119) * t498 + (Ifges(6,5) * t39 + Ifges(6,6) * t40) * t487 + (Ifges(5,5) * t110 + Ifges(5,6) * t109) * t485 + (Ifges(5,4) * t110 + Ifges(5,2) * t109) * t492 + (Ifges(5,6) * t507 + Ifges(5,5) * t508 + Ifges(5,3) * t496 + t563 / 0.2e1 + Ifges(6,5) * t517 + Ifges(6,6) * t516 - t61 * mrSges(4,3) - Ifges(4,2) * t494 - Ifges(4,4) * t495 + Ifges(6,3) * t498 - Ifges(4,6) * t481 - t67 / 0.2e1 + t567 + t581) * t262 + (Ifges(6,4) * t120 + Ifges(6,2) * t119) * t516 + (Ifges(6,4) * t39 + Ifges(6,2) * t40) * t501 - t347 * t514 + (-Ifges(5,1) * t347 + Ifges(5,4) * t198) * t508 + (Ifges(5,1) * t110 + Ifges(5,4) * t109) * t490 + (Ifges(6,1) * t120 + Ifges(6,4) * t119) * t517 + (Ifges(6,1) * t39 + Ifges(6,4) * t40) * t499 + ((-g(1) * t478 - g(2) * t327) * mrSges(3,3) + (-mrSges(3,1) * t256 - mrSges(3,2) * t257 + (m(3) * t474 - t268) * qJDD(1)) * pkin(1) + (-mrSges(3,3) * t186 + Ifges(3,1) * t257 - Ifges(3,4) * t256 + Ifges(3,5) * t304) * t326 + (t185 * mrSges(3,3) + Ifges(3,4) * t257 - Ifges(3,2) * t256 + Ifges(3,6) * t304 - t525) * t331 + ((t447 / 0.2e1 + t210 / 0.2e1 - t248 * mrSges(3,3)) * t331 + (-t446 / 0.2e1 + t131 / 0.2e1 - t209 / 0.2e1 + t448 / 0.2e1 + t451 / 0.2e1 + t449 / 0.2e1 + t139 * mrSges(4,1) - t140 * mrSges(4,2) - t251 * mrSges(3,3)) * t326 + (t331 * (Ifges(3,4) * t331 - Ifges(3,2) * t326) / 0.2e1 - t532) * t419) * qJD(2)) * t321 + (-t478 * mrSges(2,1) + t327 * mrSges(2,2) - m(3) * t420 - t267 * mrSges(3,1) - m(4) * t397 - t205 * mrSges(4,1) - m(5) * (pkin(3) * t205 + t397) - t161 * mrSges(5,1) - t160 * mrSges(5,2) - m(6) * (t205 * t316 + t397) - t153 * mrSges(6,1) - t152 * mrSges(6,2) + (-t406 - t528) * t266 + t339 * t204) * g(2) + t539 * t253 + (-Ifges(5,5) * t347 + Ifges(5,6) * t198) * t496 + (t222 / 0.2e1 + Ifges(4,1) * t482 + t133 / 0.2e1 + t578) * t197 + t421 * (-mrSges(3,2) * t304 - mrSges(3,3) * t256) + m(3) * (t185 * t421 + t186 * t271 - t248 * t253 + t251 * t252) + t271 * (mrSges(3,1) * t304 - mrSges(3,3) * t257) + t168 * t368 - t53 * t366 + t524 * t322 + t120 * t519 + t119 * t520 + t110 * t503 + t109 * t505 + t39 * t510 + t40 * t512 + t198 * t515 + (-t62 * mrSges(4,3) + 0.2e1 * t509) * t263 + (t327 * mrSges(2,1) + t478 * mrSges(2,2) - m(3) * t381 + t265 * mrSges(3,1) + t531 * t201 + (-t324 * t545 + t363 + t441 + t528) * t264 + t339 * t200 + t410 * (t265 * pkin(2) - t381)) * g(1) + t252 * t247 + t240 * t75 + Ifges(2,3) * qJDD(1) + t91 * t187 + t92 * t188 + t173 * t103 + t174 * t104 + t156 * t38 + t35 * t118 + t33 * (-mrSges(6,1) * t119 + mrSges(6,2) * t120) + t113 * (-mrSges(5,1) * t109 + mrSges(5,2) * t110) + t34 * t117 + t88 * t102 + t98 * t9 + t82 * (-mrSges(6,1) * t40 + mrSges(6,2) * t39) + t76 * t49 + t77 * t50 + t4 * t78 + t5 * t79 + (Ifges(5,6) * t492 - t575 / 0.2e1 - Ifges(4,4) * t482 + Ifges(5,3) * t485 + Ifges(5,5) * t490 - t579 + t583) * t196 + m(6) * (t15 * t5 + t16 * t4 + t2 * t27 + t26 * t3 + t33 * t98 + t54 * t82) + m(5) * (t113 * t88 + t12 * t77 + t13 * t76 + t156 * t53 + t34 * t64 + t35 * t63) + m(4) * (t139 * t92 + t140 * t91 + t168 * t240 + t173 * t62 + t174 * t61 + t211 * t253) + t26 * t17 + t27 * t18 + t54 * t55; (-t140 * (-mrSges(4,3) * t325 * t331 - mrSges(4,2) * t326) - t139 * (mrSges(4,1) * t326 - mrSges(4,3) * t427)) * t419 + (-t401 - t178) * t188 + (-t402 - t179) * t187 + (t268 - t410 * t422 + (t580 * t331 + (-mrSges(4,3) - t406 + t562) * t326) * t321) * g(3) - ((-Ifges(3,2) * t391 + t330 * t133 + t325 * t553 + t210 + t301) * t331 + t228 * (Ifges(4,6) * t326 + t331 * t359) + t305 * (Ifges(3,5) * t331 - Ifges(3,6) * t326) + t326 * t131 + t289 * (Ifges(4,3) * t326 + t331 * t356) + t229 * (Ifges(4,5) * t326 + t331 * t362)) * t419 / 0.2e1 - t563 * t330 / 0.2e1 + (-t139 * t178 - t140 * t179 - pkin(2) * t168 + ((-t139 * t330 - t140 * t325) * qJD(3) + t535) * pkin(8)) * m(4) + t535 * mrSges(4,3) + t532 * qJD(1) ^ 2 * t321 ^ 2 + t289 * t211 * (mrSges(4,1) * t325 + mrSges(4,2) * t330) + (-m(4) * t211 + t379 - t539) * t251 + t540 * t55 + (mrSges(5,1) * t542 - mrSges(5,2) * t541) * t113 + (-t12 * t432 - t13 * t431 + t541 * t63 - t542 * t64) * mrSges(5,3) + (-t103 + t38) * pkin(8) * t325 + (t378 - t247) * t248 + (t382 * t84 + t383 * t85) * t329 + (mrSges(6,2) * t537 + mrSges(6,3) * t543) * t16 + (-mrSges(6,1) * t543 + mrSges(6,2) * t544) * t82 + (-mrSges(6,1) * t537 - mrSges(6,3) * t544) * t15 + t2 * (mrSges(6,2) * t330 - mrSges(6,3) * t254) + (-Ifges(6,4) * t255 - Ifges(6,2) * t254 - Ifges(6,6) * t330) * t516 + (-Ifges(6,1) * t255 - Ifges(6,4) * t254 - Ifges(6,5) * t330) * t517 + (-Ifges(6,5) * t255 - Ifges(6,6) * t254 - Ifges(6,3) * t330) * t498 + t33 * (mrSges(6,1) * t254 - mrSges(6,2) * t255) + t3 * (-mrSges(6,1) * t330 + mrSges(6,3) * t255) + (t228 * t359 + t229 * t362 + t289 * t356) * qJD(3) / 0.2e1 + (t132 / 0.2e1 + t521) * t375 + t524 + (Ifges(6,4) * t124 + Ifges(6,2) * t125) * t501 + (Ifges(6,4) * t151 + Ifges(6,2) * t150) * t502 + t282 * t9 + t584 * t506 - t168 * t367 + t209 * t391 / 0.2e1 - t414 * t453 + (Ifges(5,4) * t220 + Ifges(5,2) * t219) * t493 - t255 * t519 - t254 * t520 + t220 * t504 + (-Ifges(5,6) * t330 + t325 * t358) * t507 + (-Ifges(5,5) * t330 + t325 * t361) * t508 + t325 * t509 + t124 * t510 + t151 * t511 + t125 * t512 + t150 * t513 + t431 * t514 + t324 * t85 * t382 + t330 * t557 + (-t357 * t412 + (Ifges(5,6) * t325 + t330 * t358) * qJD(3)) * t492 + (Ifges(4,2) * t330 + t461) * t494 + (Ifges(4,1) * t325 + t460) * t495 + (-Ifges(5,3) * t330 + t325 * t355) * t496 - t29 * t432 / 0.2e1 + (Ifges(4,5) * t325 + Ifges(4,6) * t330) * t481 + (-t354 * t412 + (Ifges(5,3) * t325 + t330 * t355) * qJD(3)) * t485 + (-t360 * t412 + (Ifges(5,5) * t325 + t330 * t361) * qJD(3)) * t490 + t104 * t468 + (Ifges(6,1) * t124 + Ifges(6,4) * t125) * t499 + (Ifges(6,1) * t151 + Ifges(6,4) * t150) * t500 + (t266 * t558 + t267 * t522) * g(1) + (t264 * t558 + t265 * t522) * g(2) + t559 * t102 + (Ifges(5,5) * t220 + Ifges(5,6) * t219) * t486 + t547 * t118 + t548 * t117 + (-t113 * t158 + t12 * t239 + t13 * t238 + (t113 * t414 + t444) * pkin(8) + t548 * t64 + t547 * t63) * m(5) + (Ifges(6,5) * t124 + Ifges(6,6) * t125) * t487 + (Ifges(6,5) * t151 + Ifges(6,6) * t150) * t488 + (Ifges(5,1) * t220 + Ifges(5,4) * t219) * t491 + t238 * t49 + t239 * t50 + t115 * t17 + t116 * t18 - pkin(2) * t75 + t133 * t383 + (-t452 + t583) * t415 + t554 * t78 + t555 * t79 + (t115 * t3 + t116 * t2 + t15 * t555 + t16 * t554 + t282 * t33 + t540 * t82) * m(6) - t330 * t556 + t330 * t67 / 0.2e1 + t364 * t444; (-pkin(3) * t53 - t113 * t140 - t63 * t80 - t64 * t81) * m(5) + (m(5) * ((-t324 * t64 - t329 * t63) * qJD(4) + t534) - t118 * t411 - t117 * t413 - t324 * t49 + t329 * t50) * pkin(9) + (t204 * t531 + t205 * t339) * g(1) + (-t200 * t531 + t201 * t339) * g(2) - (Ifges(6,4) * t499 + Ifges(6,2) * t501 + Ifges(6,6) * t487 + t512 + t530) * t194 + (t148 * t16 - t149 * t15) * mrSges(6,3) + (t188 - t102) * t140 + (t262 * t538 + t263 * t536 + t368) * g(3) + (t355 * t486 + t358 * t493 + t361 * t491 - t346 - t578) * t228 + (-Ifges(4,2) * t484 + t521 + t579) * t229 + (-Ifges(6,4) * t149 - Ifges(6,2) * t148) * t502 + (-Ifges(6,1) * t149 - Ifges(6,4) * t148) * t500 + (-Ifges(6,5) * t149 - Ifges(6,6) * t148) * t488 - t82 * (mrSges(6,1) * t148 - mrSges(6,2) * t149) + (t183 * t358 + t184 * t361 + t223 * t355) * qJD(4) / 0.2e1 + t525 - t316 * t9 + t53 * t365 + (-t2 * mrSges(6,3) + t33 * mrSges(6,1) - 0.2e1 * t520 - (Ifges(6,1) * t499 + Ifges(6,4) * t501 + Ifges(6,5) * t487 + t510 + t529) * t533) * t352 + t411 * t503 + t437 * t504 + t438 * t505 + t413 * t506 + t357 * t507 + t360 * t508 - t149 * t511 - t148 * t513 + t324 * t514 + t329 * t515 + (t33 * mrSges(6,2) - t3 * mrSges(6,3) + 0.2e1 * t519) * t277 + t354 * t496 + (-t568 * t64 + (-t411 + t437) * t63 + t534) * mrSges(5,3) + t132 * t482 + t546 * t55 + t549 * t79 + t550 * t78 + (t15 * t549 + t16 * t550 + t2 * t214 + t213 * t3 - t316 * t33 + t546 * t82) * m(6) + (t222 + t133) * t484 + t213 * t17 + t214 * t18 - t139 * t187 - t80 * t118 - t81 * t117 - (Ifges(4,1) * t228 - t450 + t553) * t229 / 0.2e1 - pkin(3) * t38 + qJD(4) * t346; (Ifges(6,1) * t500 + Ifges(6,4) * t502 + Ifges(6,5) * t488 + t511 - t529) * t380 - (Ifges(6,4) * t500 + Ifges(6,2) * t502 + Ifges(6,6) * t488 + t513 - t530) * t97 + t581 + t350 + t28 + (t2 * t323 + t3 * t328 + (-t15 * t323 + t16 * t328) * qJD(5)) * t518 + (Ifges(5,1) * t183 - t459) * t491 - t55 * t472 - m(6) * (t15 * t19 + t16 * t20 + t472 * t82) + (Ifges(5,5) * t183 - Ifges(5,6) * t184) * t486 + t84 * t490 + (-t426 - (-t201 * t329 - t264 * t324) * mrSges(5,2) + t545 * (-t201 * t324 + t264 * t329)) * g(2) + (mrSges(5,2) * t161 + t160 * t545 - t425) * g(1) + ((-t323 * t79 + t328 * t78) * qJD(5) + t328 * t17 + t323 * t18) * pkin(4) - t113 * (mrSges(5,1) * t184 + mrSges(5,2) * t183) - t20 * t78 - t19 * t79 + (-m(6) * t471 - t366 - t424) * g(3) + (-Ifges(5,2) * t184 + t181 + t85) * t493 + (t464 - t117) * t63 + (t463 + t118) * t64; -t82 * (mrSges(6,1) * t97 + mrSges(6,2) * t380) + (Ifges(6,1) * t380 - t475) * t500 + t45 * t499 + (Ifges(6,5) * t380 - Ifges(6,6) * t97) * t488 - t15 * t78 + t16 * t79 - g(1) * t425 - g(2) * t426 - g(3) * t424 + (t15 * t380 + t16 * t97) * mrSges(6,3) + t350 + (-Ifges(6,2) * t97 + t46 + t93) * t502;];
tau = t1;
