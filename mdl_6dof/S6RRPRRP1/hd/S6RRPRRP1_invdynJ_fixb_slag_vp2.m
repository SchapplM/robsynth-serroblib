% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP1
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:39:00
% EndTime: 2019-03-09 11:39:43
% DurationCPUTime: 26.66s
% Computational Cost: add. (16359->736), mult. (38054->928), div. (0->0), fcn. (28878->14), ass. (0->353)
t323 = cos(qJ(5));
t600 = mrSges(7,1) + mrSges(6,1);
t601 = t323 * t600;
t599 = Ifges(6,4) + Ifges(7,4);
t321 = sin(qJ(2));
t414 = qJD(1) * qJD(2);
t389 = t321 * t414;
t325 = cos(qJ(2));
t413 = qJDD(1) * t325;
t265 = -t389 + t413;
t266 = qJDD(1) * t321 + t325 * t414;
t315 = sin(pkin(10));
t316 = cos(pkin(10));
t197 = t265 * t316 - t266 * t315;
t198 = t265 * t315 + t266 * t316;
t320 = sin(qJ(4));
t324 = cos(qJ(4));
t252 = -t315 * t321 + t316 * t325;
t236 = t252 * qJD(1);
t422 = qJD(1) * t325;
t423 = qJD(1) * t321;
t237 = -t315 * t422 - t316 * t423;
t375 = t324 * t236 + t237 * t320;
t109 = qJD(4) * t375 + t197 * t320 + t198 * t324;
t313 = qJD(2) + qJD(4);
t319 = sin(qJ(5));
t349 = t236 * t320 - t324 * t237;
t163 = t313 * t323 - t319 * t349;
t311 = qJDD(2) + qJDD(4);
t71 = qJD(5) * t163 + t109 * t323 + t311 * t319;
t516 = t71 / 0.2e1;
t164 = t313 * t319 + t323 * t349;
t72 = -qJD(5) * t164 - t109 * t319 + t311 * t323;
t515 = t72 / 0.2e1;
t110 = -qJD(4) * t349 + t197 * t324 - t198 * t320;
t108 = qJDD(5) - t110;
t514 = t108 / 0.2e1;
t594 = -mrSges(6,3) - mrSges(7,3);
t569 = Ifges(6,1) + Ifges(7,1);
t567 = Ifges(6,5) + Ifges(7,5);
t587 = Ifges(6,2) + Ifges(7,2);
t565 = Ifges(6,6) + Ifges(7,6);
t415 = qJD(5) * t323;
t546 = t375 * t323;
t598 = -t415 + t546;
t416 = qJD(5) * t319;
t547 = t375 * t319;
t597 = t416 - t547;
t175 = qJD(5) - t375;
t593 = t599 * t163;
t561 = t569 * t164 + t567 * t175 + t593;
t596 = t561 / 0.2e1;
t595 = t567 * t514 + t599 * t515 + t569 * t516;
t570 = mrSges(6,2) + mrSges(7,2);
t564 = Ifges(6,3) + Ifges(7,3);
t314 = qJ(2) + pkin(10);
t307 = qJ(4) + t314;
t295 = sin(t307);
t296 = cos(t307);
t491 = pkin(5) * t323;
t298 = pkin(4) + t491;
t317 = -qJ(6) - pkin(9);
t592 = -m(7) * (-t295 * t298 - t296 * t317) + t295 * t601;
t591 = t599 * t164;
t590 = t599 * t323;
t589 = t599 * t319;
t588 = t295 * t570;
t586 = t565 * t108 + t587 * t72 + t599 * t71;
t562 = t587 * t163 + t565 * t175 + t591;
t495 = pkin(2) * t316;
t297 = pkin(3) + t495;
t496 = pkin(2) * t315;
t229 = t297 * t324 - t320 * t496;
t214 = t229 * qJD(4);
t318 = -qJ(3) - pkin(7);
t280 = t318 * t321;
t263 = qJD(1) * t280;
t283 = t318 * t325;
t264 = qJD(1) * t283;
t433 = t316 * t264;
t195 = -t263 * t315 + t433;
t489 = pkin(8) * t236;
t165 = t195 - t489;
t240 = t315 * t264;
t196 = t316 * t263 + t240;
t488 = pkin(8) * t237;
t166 = t196 + t488;
t117 = t165 * t320 + t166 * t324;
t137 = pkin(4) * t349 - pkin(9) * t375;
t302 = pkin(2) * t423;
t205 = -pkin(3) * t237 + t302;
t118 = t137 + t205;
t51 = t323 * t117 + t319 * t118;
t584 = t214 * t323 - t51;
t230 = t320 * t297 + t324 * t496;
t542 = t230 * qJD(4) + t324 * t165 - t166 * t320;
t281 = -t325 * mrSges(3,1) + t321 * mrSges(3,2);
t305 = sin(t314);
t306 = cos(t314);
t583 = -t306 * mrSges(4,1) + t305 * mrSges(4,2) + t281;
t477 = mrSges(7,2) * t323;
t361 = mrSges(7,1) * t319 + t477;
t362 = mrSges(6,1) * t319 + mrSges(6,2) * t323;
t247 = qJD(2) * pkin(2) + t263;
t189 = t316 * t247 + t240;
t155 = qJD(2) * pkin(3) + t189 + t488;
t190 = t315 * t247 - t433;
t159 = t190 + t489;
t96 = t155 * t324 - t320 * t159;
t90 = -pkin(4) * t313 - t96;
t73 = -pkin(5) * t163 + qJD(6) + t90;
t582 = t73 * t361 + t90 * t362;
t538 = -t565 * t319 + t567 * t323;
t537 = -t587 * t319 + t590;
t536 = t569 * t323 - t589;
t310 = t325 * pkin(2);
t299 = t310 + pkin(1);
t268 = -qJD(1) * t299 + qJD(3);
t199 = -pkin(3) * t236 + t268;
t113 = -pkin(4) * t375 - pkin(9) * t349 + t199;
t251 = t266 * pkin(7);
t420 = qJD(3) * t321;
t187 = qJDD(2) * pkin(2) - qJ(3) * t266 - qJD(1) * t420 - t251;
t300 = pkin(7) * t413;
t421 = qJD(2) * t321;
t406 = pkin(7) * t421;
t419 = qJD(3) * t325;
t193 = qJ(3) * t265 + t300 + (-t406 + t419) * qJD(1);
t138 = t316 * t187 - t193 * t315;
t103 = qJDD(2) * pkin(3) - pkin(8) * t198 + t138;
t139 = t315 * t187 + t316 * t193;
t112 = pkin(8) * t197 + t139;
t417 = qJD(4) * t324;
t418 = qJD(4) * t320;
t28 = t320 * t103 + t324 * t112 + t155 * t417 - t159 * t418;
t23 = pkin(9) * t311 + t28;
t455 = qJDD(1) * pkin(1);
t221 = -pkin(2) * t265 + qJDD(3) - t455;
t162 = -pkin(3) * t197 + t221;
t42 = -pkin(4) * t110 - pkin(9) * t109 + t162;
t97 = t155 * t320 + t159 * t324;
t91 = pkin(9) * t313 + t97;
t5 = t113 * t415 + t323 * t23 + t319 * t42 - t416 * t91;
t47 = t113 * t319 + t323 * t91;
t6 = -qJD(5) * t47 - t23 * t319 + t323 * t42;
t581 = -t319 * t6 + t323 * t5;
t322 = sin(qJ(1));
t326 = cos(qJ(1));
t533 = g(1) * t326 + g(2) * t322;
t467 = Ifges(5,4) * t349;
t555 = t375 * Ifges(5,2);
t556 = t313 * Ifges(5,6);
t132 = t467 + t555 + t556;
t580 = t97 * mrSges(5,3) + t132 / 0.2e1;
t174 = Ifges(5,4) * t375;
t557 = t313 * Ifges(5,5);
t133 = Ifges(5,1) * t349 + t174 + t557;
t579 = t96 * mrSges(5,3) - t133 / 0.2e1;
t578 = -t296 * mrSges(5,1) + (mrSges(5,2) + t594) * t295;
t577 = t597 * pkin(5);
t576 = qJ(6) * t547 + t323 * qJD(6);
t539 = t296 * pkin(4) + t295 * pkin(9);
t541 = -t295 * t317 + t296 * t298;
t575 = -m(6) * t539 - m(7) * t541;
t1 = pkin(5) * t108 - qJ(6) * t71 - qJD(6) * t164 + t6;
t29 = t103 * t324 - t320 * t112 - t155 * t418 - t159 * t417;
t24 = -pkin(4) * t311 - t29;
t12 = -pkin(5) * t72 + qJDD(6) + t24;
t3 = qJ(6) * t72 + qJD(6) * t163 + t5;
t46 = t323 * t113 - t319 * t91;
t38 = -qJ(6) * t164 + t46;
t31 = pkin(5) * t175 + t38;
t384 = -t416 / 0.2e1;
t39 = qJ(6) * t163 + t47;
t431 = t322 * t319;
t432 = t319 * t326;
t435 = t296 * t326;
t437 = t296 * t322;
t574 = t24 * (-mrSges(6,1) * t323 + mrSges(6,2) * t319) + t12 * (-mrSges(7,1) * t323 + mrSges(7,2) * t319) + Ifges(5,3) * t311 + Ifges(5,5) * t109 + Ifges(5,6) * t110 + t29 * mrSges(5,1) - t28 * mrSges(5,2) + (t319 * t569 + t590) * t516 + (t323 * t587 + t589) * t515 + (t319 * t567 + t323 * t565) * t514 + t319 * t595 + t586 * t323 / 0.2e1 + t415 * t596 + t582 * qJD(5) + (t326 * t592 - t432 * t588 + t435 * t594) * g(1) + (t322 * t592 - t431 * t588 + t437 * t594) * g(2) + (t163 * t537 + t164 * t536 + t175 * t538) * qJD(5) / 0.2e1 - t561 * t546 / 0.2e1 + (t384 + t547 / 0.2e1) * t562 + (-t1 * t319 + t3 * t323 + t598 * t31 - t597 * t39) * mrSges(7,3) + (t598 * t46 - t597 * t47 + t581) * mrSges(6,3);
t517 = m(7) * pkin(5);
t573 = t265 / 0.2e1;
t572 = t266 / 0.2e1;
t504 = -t349 / 0.2e1;
t479 = qJD(2) / 0.2e1;
t563 = t163 * t565 + t164 * t567 + t175 * t564;
t253 = t315 * t325 + t316 * t321;
t560 = Ifges(4,5) * t253;
t559 = Ifges(4,6) * t252;
t558 = t199 * mrSges(5,2);
t456 = qJDD(2) / 0.2e1;
t227 = pkin(9) + t230;
t428 = -qJ(6) - t227;
t372 = qJD(5) * t428;
t554 = t319 * t372 + t576 + t584;
t50 = -t117 * t319 + t323 * t118;
t309 = t323 * qJ(6);
t526 = pkin(5) * t349 - t309 * t375;
t553 = (-qJD(6) - t214) * t319 + t323 * t372 - t50 - t526;
t552 = t577 + t542;
t380 = qJD(5) * t317;
t55 = t319 * t137 + t323 * t96;
t551 = t319 * t380 - t55 + t576;
t54 = t323 * t137 - t319 * t96;
t550 = -qJD(6) * t319 + t323 * t380 - t526 - t54;
t549 = t577 - t97;
t548 = t517 + mrSges(7,1);
t200 = t316 * t280 + t283 * t315;
t172 = -pkin(8) * t253 + t200;
t201 = t315 * t280 - t316 * t283;
t173 = pkin(8) * t252 + t201;
t131 = t172 * t320 + t173 * t324;
t128 = t323 * t131;
t192 = t252 * t320 + t253 * t324;
t210 = -pkin(3) * t252 - t299;
t348 = t324 * t252 - t253 * t320;
t129 = -pkin(4) * t348 - pkin(9) * t192 + t210;
t57 = t319 * t129 + t128;
t475 = mrSges(5,3) * t349;
t545 = mrSges(5,1) * t313 + mrSges(6,1) * t163 - mrSges(6,2) * t164 - t475;
t544 = t324 * t172 - t173 * t320;
t543 = t214 - t117;
t535 = t108 * t564 + t565 * t72 + t567 * t71;
t250 = -pkin(7) * t389 + t300;
t534 = t250 * t325 + t251 * t321;
t532 = m(7) + m(6) + m(5);
t531 = mrSges(6,1) + t548;
t529 = t578 + (t319 * t570 - t601) * t296;
t528 = -m(3) * pkin(7) + m(4) * t318 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t527 = -m(6) * t90 + t545;
t474 = mrSges(6,3) * t163;
t123 = -mrSges(6,2) * t175 + t474;
t473 = mrSges(6,3) * t164;
t125 = mrSges(6,1) * t175 - t473;
t34 = mrSges(6,1) * t108 - mrSges(6,3) * t71;
t525 = m(6) * ((-t319 * t47 - t323 * t46) * qJD(5) + t581) - t125 * t415 - t123 * t416 - t319 * t34;
t524 = m(3) * pkin(1) + m(4) * t299 + mrSges(2,1) - t578 - t583;
t523 = t6 * mrSges(6,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t521 = -t199 * mrSges(5,1) - t46 * mrSges(6,1) - t31 * mrSges(7,1) + t47 * mrSges(6,2) + t39 * mrSges(7,2);
t499 = -t313 / 0.2e1;
t508 = -t175 / 0.2e1;
t510 = -t164 / 0.2e1;
t512 = -t163 / 0.2e1;
t520 = Ifges(5,5) * t499 + t508 * t538 + t510 * t536 + t512 * t537 - t558 - t582;
t505 = -t375 / 0.2e1;
t519 = -Ifges(5,2) * t505 - Ifges(5,6) * t499 + t508 * t564 + t510 * t567 + t512 * t565 + t521;
t511 = t163 / 0.2e1;
t509 = t164 / 0.2e1;
t507 = t175 / 0.2e1;
t503 = t349 / 0.2e1;
t500 = -t237 / 0.2e1;
t494 = pkin(2) * t321;
t493 = pkin(4) * t295;
t490 = pkin(7) * t325;
t487 = pkin(9) * t323;
t476 = mrSges(5,3) * t375;
t472 = mrSges(7,3) * t163;
t471 = mrSges(7,3) * t164;
t470 = Ifges(3,4) * t321;
t469 = Ifges(3,4) * t325;
t468 = Ifges(4,4) * t237;
t460 = t189 * mrSges(4,3);
t459 = t190 * mrSges(4,3);
t238 = t253 * qJD(2);
t239 = t252 * qJD(2);
t140 = qJD(4) * t348 - t238 * t320 + t239 * t324;
t454 = t140 * t319;
t453 = t140 * t323;
t447 = t192 * t319;
t446 = t192 * t323;
t444 = t227 * t323;
t430 = t322 * t323;
t429 = t323 * t326;
t381 = qJD(2) * t318;
t233 = t321 * t381 + t419;
t234 = t325 * t381 - t420;
t171 = t316 * t233 + t315 * t234;
t424 = pkin(3) * t306 + t310;
t411 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t423) * t490;
t170 = -t233 * t315 + t316 * t234;
t149 = -pkin(8) * t239 + t170;
t150 = -pkin(8) * t238 + t171;
t52 = qJD(4) * t544 + t149 * t320 + t150 * t324;
t141 = qJD(4) * t192 + t324 * t238 + t239 * t320;
t303 = pkin(2) * t421;
t206 = pkin(3) * t238 + t303;
t65 = pkin(4) * t141 - pkin(9) * t140 + t206;
t410 = t129 * t415 + t319 * t65 + t323 * t52;
t399 = t192 * t415;
t25 = -t72 * mrSges(7,1) + t71 * mrSges(7,2);
t382 = -t319 * t52 + t323 * t65;
t378 = -t197 * mrSges(4,1) + t198 * mrSges(4,2);
t377 = -t110 * mrSges(5,1) + t109 * mrSges(5,2);
t56 = t323 * t129 - t131 * t319;
t369 = pkin(9) * t437 - t322 * t493;
t368 = pkin(9) * t435 - t326 * t493;
t226 = -pkin(4) - t229;
t366 = mrSges(3,1) * t321 + mrSges(3,2) * t325;
t363 = mrSges(5,1) * t295 + mrSges(5,2) * t296;
t358 = Ifges(3,2) * t325 + t470;
t355 = Ifges(3,5) * t325 - Ifges(3,6) * t321;
t351 = -t319 * t46 + t323 * t47;
t346 = -qJ(6) * t140 - qJD(6) * t192;
t343 = pkin(1) * t366;
t224 = -t296 * t432 + t430;
t222 = t296 * t431 + t429;
t340 = t399 + t454;
t339 = t192 * t416 - t453;
t338 = t321 * (Ifges(3,1) * t325 - t470);
t53 = qJD(4) * t131 - t324 * t149 + t150 * t320;
t312 = -pkin(8) + t318;
t301 = Ifges(3,4) * t422;
t282 = t309 + t487;
t279 = t317 * t319;
t278 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t422;
t267 = -pkin(3) * t305 - t494;
t262 = pkin(1) + t424;
t249 = t326 * t267;
t248 = t322 * t267;
t246 = Ifges(3,1) * t423 + Ifges(3,5) * qJD(2) + t301;
t245 = Ifges(3,6) * qJD(2) + qJD(1) * t358;
t228 = Ifges(4,4) * t236;
t225 = t296 * t429 + t431;
t223 = -t296 * t430 + t432;
t209 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t237;
t208 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t236;
t207 = t226 - t491;
t203 = t309 + t444;
t202 = t428 * t319;
t188 = -mrSges(4,1) * t236 - mrSges(4,2) * t237;
t185 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t198;
t184 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t197;
t179 = -t237 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t228;
t178 = t236 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t468;
t167 = -mrSges(5,2) * t313 + t476;
t136 = -mrSges(5,1) * t375 + mrSges(5,2) * t349;
t124 = mrSges(7,1) * t175 - t471;
t122 = -mrSges(7,2) * t175 + t472;
t119 = -mrSges(7,1) * t163 + mrSges(7,2) * t164;
t93 = -mrSges(5,2) * t311 + mrSges(5,3) * t110;
t92 = mrSges(5,1) * t311 - mrSges(5,3) * t109;
t87 = pkin(5) * t447 - t544;
t48 = -qJ(6) * t447 + t57;
t44 = -pkin(5) * t348 - t192 * t309 + t56;
t36 = -mrSges(6,2) * t108 + mrSges(6,3) * t72;
t35 = -mrSges(7,2) * t108 + mrSges(7,3) * t72;
t33 = mrSges(7,1) * t108 - mrSges(7,3) * t71;
t30 = pkin(5) * t340 + t53;
t26 = -mrSges(6,1) * t72 + mrSges(6,2) * t71;
t10 = -qJD(5) * t57 + t382;
t9 = -t131 * t416 + t410;
t8 = -qJ(6) * t399 + (-qJD(5) * t131 + t346) * t319 + t410;
t7 = pkin(5) * t141 + t346 * t323 + (-t128 + (qJ(6) * t192 - t129) * t319) * qJD(5) + t382;
t2 = [(-t454 / 0.2e1 - t399 / 0.2e1) * t562 + t446 * t595 + t453 * t596 + (Ifges(5,1) * t503 + t557 / 0.2e1 + t558 + t174 / 0.2e1 - t579) * t140 + (-Ifges(5,4) * t503 - t556 / 0.2e1 - t555 / 0.2e1 + t565 * t511 + t567 * t509 + t564 * t507 - t521 + t563 / 0.2e1 - t580) * t141 + (Ifges(4,5) * t239 - Ifges(4,6) * t238) * t479 + t268 * (mrSges(4,1) * t238 + mrSges(4,2) * t239) + t236 * (Ifges(4,4) * t239 - Ifges(4,2) * t238) / 0.2e1 + (Ifges(4,1) * t239 - Ifges(4,4) * t238) * t500 + (t325 * (-Ifges(3,2) * t321 + t469) + t338) * t414 / 0.2e1 + (-t1 * t446 - t3 * t447 + t31 * t339 - t340 * t39) * mrSges(7,3) + (-t339 * t599 - t340 * t587) * t511 + (-t339 * t569 - t340 * t599) * t509 + m(6) * (t10 * t46 + t47 * t9 + t5 * t57 + t56 * t6) + m(5) * (t131 * t28 + t162 * t210 + t199 * t206 + t52 * t97) + (t162 * mrSges(5,2) - t29 * mrSges(5,3) + Ifges(5,1) * t109 + Ifges(5,4) * t110 + Ifges(5,5) * t311 + t12 * t361 + t24 * t362 + t384 * t561 + t514 * t538 + t515 * t537 + t516 * t536) * t192 + (-t138 * t253 + t139 * t252) * mrSges(4,3) + (t339 * t46 - t340 * t47 - t446 * t6 - t447 * t5) * mrSges(6,3) - t245 * t421 / 0.2e1 + t469 * t572 + t358 * t573 + (t265 * t490 + t534) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t534) + (-m(5) * t96 - t527) * t53 + (t560 / 0.2e1 + Ifges(3,6) * t325 / 0.2e1 - mrSges(3,2) * t490 + t559 / 0.2e1) * qJDD(2) + (t559 + t560) * t456 + (-t339 * t567 - t340 * t565) * t507 + (t355 * t479 - t411) * qJD(2) - t343 * t414 - t238 * t459 - (-m(5) * t29 + m(6) * t24 + t26 - t92) * t544 - t239 * t460 - t281 * t455 + t188 * t303 + (Ifges(3,4) * t572 + Ifges(3,2) * t573 + Ifges(3,6) * t456 + t246 * t479) * t325 + (Ifges(3,1) * t266 + Ifges(3,4) * t573 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t266) + 0.2e1 * Ifges(3,5) * t456) * t321 + m(7) * (t1 * t44 + t12 * t87 + t3 * t48 + t30 * t73 + t31 * t7 + t39 * t8) + m(4) * (t138 * t200 + t139 * t201 + t170 * t189 + t171 * t190 - t221 * t299 + t268 * t303) - t278 * t406 + (-t431 * t517 - t532 * (t326 * t262 - t322 * t312) - t600 * t225 - t570 * t224 + t528 * t322 + (-t524 + t575) * t326) * g(2) + (-t600 * t223 - t570 * t222 + (t312 * t532 - t319 * t517 + t528) * t326 + (-m(7) * (-t262 - t541) - m(6) * (-t262 - t539) + m(5) * t262 + t524) * t322) * g(1) + t210 * t377 - t299 * t378 - pkin(1) * (-mrSges(3,1) * t265 + mrSges(3,2) * t266) + t221 * (-mrSges(4,1) * t252 + mrSges(4,2) * t253) + t197 * (Ifges(4,4) * t253 + Ifges(4,2) * t252) + t198 * (Ifges(4,1) * t253 + Ifges(4,4) * t252) - t238 * t178 / 0.2e1 + t239 * t179 / 0.2e1 + Ifges(2,3) * qJDD(1) + t206 * t136 + t171 * t208 + t170 * t209 + t200 * t185 + t201 * t184 + t52 * t167 + (-t535 / 0.2e1 - t162 * mrSges(5,1) - t1 * mrSges(7,1) + t28 * mrSges(5,3) + Ifges(5,4) * t109 + Ifges(5,2) * t110 + Ifges(5,6) * t311 - t514 * t564 - t515 * t565 - t516 * t567 - t523) * t348 + t131 * t93 + t30 * t119 + t8 * t122 + t9 * t123 + t7 * t124 + t10 * t125 + t73 * (mrSges(7,1) * t340 - mrSges(7,2) * t339) + t90 * (mrSges(6,1) * t340 - mrSges(6,2) * t339) + t87 * t25 + t56 * t34 + t57 * t36 + t48 * t35 + t44 * t33 - t586 * t447 / 0.2e1; t574 + t584 * t123 + (-m(7) * (t541 + t424) - m(5) * t424 - m(6) * (t424 + t539) - m(4) * t310 + t529 + t583) * g(3) - t542 * t545 + t533 * (m(4) * t494 - m(5) * t267 + mrSges(4,1) * t305 + mrSges(4,2) * t306 + t363 + t366) + (Ifges(5,1) * t504 + Ifges(5,4) * t505 + t520 + t579) * t375 + (-Ifges(5,4) * t504 + t519 + t580) * t349 - (Ifges(4,2) * t237 + t179 + t228) * t236 / 0.2e1 - (-Ifges(3,2) * t423 + t246 + t301) * t422 / 0.2e1 + t525 * t227 + t563 * t504 + t552 * t119 + t553 * t124 + (-g(1) * t249 - g(2) * t248 + t1 * t202 + t12 * t207 + t203 * t3 + t553 * t31 + t554 * t39 + t552 * t73) * m(7) + t554 * t122 + (t411 + (-t338 / 0.2e1 + t343) * qJD(1)) * qJD(1) - t355 * t414 / 0.2e1 + t236 * t460 - t237 * t459 + (-t46 * t50 - t47 * t51 + t214 * t351 + t226 * t24 - g(1) * (t249 + t368) - g(2) * (t248 + t369) + t542 * t90) * m(6) + (-t199 * t205 + t229 * t29 + t230 * t28 - t542 * t96 + t543 * t97) * m(5) + t543 * t167 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t237 * (Ifges(4,1) * t236 + t468) / 0.2e1 + t185 * t495 + t184 * t496 + t178 * t500 + (t245 / 0.2e1 + pkin(7) * t278) * t423 + Ifges(3,6) * t265 + Ifges(3,5) * t266 - t268 * (-mrSges(4,1) * t237 + mrSges(4,2) * t236) - t250 * mrSges(3,2) - t251 * mrSges(3,1) - qJD(2) * (Ifges(4,5) * t236 + Ifges(4,6) * t237) / 0.2e1 + t229 * t92 + t230 * t93 + t226 * t26 - t205 * t136 + t207 * t25 - t196 * t208 - t195 * t209 + Ifges(4,6) * t197 + Ifges(4,5) * t198 + t202 * t33 + t203 * t35 + t138 * mrSges(4,1) - t139 * mrSges(4,2) - t188 * t302 + ((t138 * t316 + t139 * t315) * pkin(2) - t189 * t195 - t190 * t196 - t268 * t302) * m(4) + (-t214 * t319 - t50) * t125 + t36 * t444; -t375 * t167 - t236 * t208 - t237 * t209 + (-t119 + t545) * t349 + (t33 + t34 + t175 * (t122 + t123)) * t323 + (t35 + t36 - t175 * (t124 + t125)) * t319 + t377 + t378 + (-g(1) * t322 + g(2) * t326) * (m(4) + t532) + (t1 * t323 - t349 * t73 + t3 * t319 + t175 * (-t31 * t319 + t323 * t39)) * m(7) + (t175 * t351 + t319 * t5 + t323 * t6 - t349 * t90) * m(6) + (t349 * t96 - t375 * t97 + t162) * m(5) + (-t189 * t237 - t190 * t236 + t221) * m(4); t519 * t349 + t520 * t375 + t574 + (t476 - t167) * t96 + t525 * pkin(9) + t533 * t363 + (t174 + t133) * t505 + (t475 + t527) * t97 + t549 * t119 + t550 * t124 + (t1 * t279 - t12 * t298 + t282 * t3 + t31 * t550 + t39 * t551 + t549 * t73) * m(7) + t551 * t122 + t132 * t503 + (-pkin(4) * t24 - g(1) * t368 - g(2) * t369 - t46 * t54 - t47 * t55) * m(6) + (Ifges(5,1) * t375 - t467 + t563) * t504 + t36 * t487 - t298 * t25 + t279 * t33 + t282 * t35 + (t529 + t575) * g(3) - t55 * t123 - t54 * t125 - pkin(4) * t26; (t163 * t569 - t591) * t510 + (t163 * t567 - t164 * t565) * t508 + (t222 * t531 - t223 * t570) * g(2) + (-t224 * t531 + t225 * t570) * g(1) + (-t164 * t587 + t561 + t593) * t512 + (t319 * t548 + t362 + t477) * g(3) * t295 + t548 * t1 + (t474 - t123) * t46 + (-m(7) * (-t31 + t38) + t471 + t124) * t39 + t31 * t472 + (t473 + t125) * t47 - t73 * (mrSges(7,1) * t164 + mrSges(7,2) * t163) - t90 * (mrSges(6,1) * t164 + mrSges(6,2) * t163) - t38 * t122 + t562 * t509 + t523 + t535 + ((-m(7) * t73 - t119) * t164 + t33) * pkin(5); -t163 * t122 + t164 * t124 + (g(3) * t296 - t39 * t163 + t31 * t164 - t295 * t533 + t12) * m(7) + t25;];
tau  = t2;
