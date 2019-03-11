% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:45
% EndTime: 2019-03-09 13:56:46
% DurationCPUTime: 37.94s
% Computational Cost: add. (15590->868), mult. (31885->1133), div. (0->0), fcn. (21646->12), ass. (0->404)
t331 = -pkin(10) - pkin(9);
t550 = -m(7) * t331 + mrSges(6,3) + mrSges(7,3);
t563 = -m(6) * pkin(9) - t550;
t636 = mrSges(5,2) + t563;
t325 = sin(qJ(1));
t323 = sin(qJ(4));
t324 = sin(qJ(2));
t328 = cos(qJ(4));
t329 = cos(qJ(2));
t562 = t329 * t323 - t324 * t328;
t196 = t562 * t325;
t229 = t323 * t324 + t328 * t329;
t197 = t229 * t325;
t322 = sin(qJ(5));
t327 = cos(qJ(5));
t258 = -mrSges(6,1) * t327 + mrSges(6,2) * t322;
t501 = t327 * pkin(5);
t291 = pkin(4) + t501;
t318 = qJ(5) + qJ(6);
t303 = sin(t318);
t304 = cos(t318);
t564 = m(6) * pkin(4) + m(7) * t291 + mrSges(7,1) * t304 - mrSges(7,2) * t303 - t258;
t632 = -mrSges(5,1) - t564;
t635 = -t196 * t632 + t636 * t197;
t330 = cos(qJ(1));
t447 = t329 * t330;
t450 = t324 * t330;
t199 = -t323 * t450 - t328 * t447;
t200 = t323 * t447 - t328 * t450;
t634 = -t636 * t199 - t200 * t632;
t580 = -m(7) - m(6);
t633 = -m(5) + t580;
t619 = Ifges(4,4) + Ifges(3,5);
t618 = Ifges(4,6) - Ifges(3,6);
t332 = -pkin(2) - pkin(3);
t249 = -t323 * qJ(3) + t328 * t332;
t189 = t328 * qJD(3) + qJD(4) * t249;
t440 = qJD(1) * t329;
t441 = qJD(1) * t324;
t213 = -t323 * t441 - t328 * t440;
t215 = -t323 * t440 + t328 * t441;
t151 = pkin(4) * t215 - pkin(9) * t213;
t287 = qJ(3) * t440;
t416 = t332 * t324;
t194 = qJD(1) * t416 + t287;
t108 = -t151 + t194;
t297 = pkin(7) * t441;
t241 = pkin(8) * t441 - t297;
t298 = pkin(7) * t440;
t242 = -pkin(8) * t440 + t298;
t161 = t241 * t328 + t242 * t323;
t68 = t322 * t108 + t327 * t161;
t631 = t189 * t327 - t68;
t67 = t327 * t108 - t161 * t322;
t630 = -t189 * t322 - t67;
t629 = -qJD(3) + t241;
t250 = t328 * qJ(3) + t323 * t332;
t240 = -pkin(9) + t250;
t500 = pkin(10) - t240;
t393 = qJD(5) * t500;
t461 = t213 * t322;
t425 = pkin(10) * t461;
t628 = t322 * t393 - t425 + t631;
t460 = t213 * t327;
t354 = -pkin(5) * t215 + pkin(10) * t460;
t627 = t327 * t393 - t354 + t630;
t411 = qJD(5) * t331;
t410 = t332 * qJD(2);
t188 = t410 - t629;
t317 = qJD(2) * qJ(3);
t216 = t242 + t317;
t134 = t188 * t328 - t323 * t216;
t77 = t327 * t134 + t322 * t151;
t626 = t322 * t411 + t425 - t77;
t76 = -t134 * t322 + t327 * t151;
t625 = t327 * t411 + t354 - t76;
t430 = qJD(1) * qJD(2);
t247 = -t329 * qJDD(1) + t324 * t430;
t248 = qJDD(1) * t324 + t329 * t430;
t344 = t229 * qJD(4);
t114 = -qJD(1) * t344 + t247 * t323 + t248 * t328;
t434 = qJD(4) * t328;
t435 = qJD(4) * t323;
t602 = -t324 * t434 + t329 * t435;
t115 = qJD(1) * t602 + t247 * t328 - t323 * t248;
t135 = t188 * t323 + t216 * t328;
t204 = Ifges(5,4) * t215;
t314 = -qJD(2) + qJD(4);
t615 = Ifges(5,6) * t314;
t616 = Ifges(5,2) * t213;
t136 = t204 + t615 + t616;
t203 = Ifges(5,4) * t213;
t614 = t314 * Ifges(5,5);
t137 = t215 * Ifges(5,1) + t203 + t614;
t326 = cos(qJ(6));
t206 = qJD(5) - t213;
t172 = t215 * t327 + t314 * t322;
t217 = -qJD(1) * pkin(1) - pkin(2) * t440 - qJ(3) * t441;
t183 = pkin(3) * t440 - t217;
t106 = -pkin(4) * t213 - pkin(9) * t215 + t183;
t122 = pkin(9) * t314 + t135;
t62 = t327 * t106 - t122 * t322;
t50 = -pkin(10) * t172 + t62;
t37 = pkin(5) * t206 + t50;
t321 = sin(qJ(6));
t171 = -t215 * t322 + t314 * t327;
t63 = t106 * t322 + t122 * t327;
t51 = pkin(10) * t171 + t63;
t475 = t321 * t51;
t16 = t326 * t37 - t475;
t431 = qJD(6) * t326;
t432 = qJD(5) * t327;
t454 = t321 * t322;
t554 = qJD(5) + qJD(6);
t162 = -t326 * t432 - t327 * t431 + t454 * t554;
t473 = t326 * t51;
t17 = t321 * t37 + t473;
t302 = t324 * qJD(3);
t319 = qJDD(1) * pkin(1);
t147 = t247 * pkin(2) - t248 * qJ(3) - qJD(1) * t302 - t319;
t112 = -pkin(3) * t247 - t147;
t42 = -pkin(4) * t115 - pkin(9) * t114 + t112;
t433 = qJD(5) * t322;
t313 = -qJDD(2) + qJDD(4);
t225 = t248 * pkin(7);
t397 = qJDD(3) + t225;
t149 = -pkin(8) * t248 + qJDD(2) * t332 + t397;
t224 = t247 * pkin(7);
t184 = qJDD(2) * qJ(3) + qJD(2) * qJD(3) - t224;
t153 = pkin(8) * t247 + t184;
t58 = t323 * t149 + t328 * t153 + t188 * t434 - t216 * t435;
t54 = pkin(9) * t313 + t58;
t11 = t106 * t432 - t122 * t433 + t322 * t42 + t327 * t54;
t74 = -qJD(5) * t172 - t114 * t322 + t313 * t327;
t10 = pkin(10) * t74 + t11;
t113 = qJDD(5) - t115;
t12 = -qJD(5) * t63 - t322 * t54 + t327 * t42;
t73 = qJD(5) * t171 + t114 * t327 + t313 * t322;
t8 = pkin(5) * t113 - pkin(10) * t73 + t12;
t2 = qJD(6) * t16 + t10 * t326 + t321 * t8;
t228 = -t326 * t327 + t454;
t231 = t321 * t327 + t322 * t326;
t27 = Ifges(6,4) * t73 + Ifges(6,2) * t74 + Ifges(6,6) * t113;
t3 = -qJD(6) * t17 - t10 * t321 + t326 * t8;
t59 = t149 * t328 - t323 * t153 - t188 * t435 - t216 * t434;
t55 = -pkin(4) * t313 - t59;
t33 = -pkin(5) * t74 + t55;
t366 = Ifges(6,5) * t327 - Ifges(6,6) * t322;
t490 = Ifges(6,4) * t327;
t371 = -Ifges(6,2) * t322 + t490;
t491 = Ifges(6,4) * t322;
t376 = Ifges(6,1) * t327 - t491;
t398 = -t433 / 0.2e1;
t513 = t327 / 0.2e1;
t168 = Ifges(6,4) * t171;
t84 = t172 * Ifges(6,1) + t206 * Ifges(6,5) + t168;
t413 = t84 * t513;
t478 = t134 * mrSges(5,3);
t497 = mrSges(5,3) * t215;
t517 = t215 / 0.2e1;
t520 = t206 / 0.2e1;
t521 = -t206 / 0.2e1;
t524 = t172 / 0.2e1;
t525 = -t172 / 0.2e1;
t526 = t171 / 0.2e1;
t527 = -t171 / 0.2e1;
t529 = t113 / 0.2e1;
t537 = t74 / 0.2e1;
t538 = t73 / 0.2e1;
t193 = qJD(6) + t206;
t389 = t326 * t171 - t172 * t321;
t87 = Ifges(7,4) * t389;
t94 = t171 * t321 + t172 * t326;
t47 = Ifges(7,1) * t94 + Ifges(7,5) * t193 + t87;
t539 = -t47 / 0.2e1;
t511 = Ifges(7,4) * t94;
t46 = Ifges(7,2) * t389 + Ifges(7,6) * t193 + t511;
t540 = -t46 / 0.2e1;
t541 = Ifges(6,1) * t538 + Ifges(6,4) * t537 + Ifges(6,5) * t529;
t110 = qJDD(6) + t113;
t530 = t110 / 0.2e1;
t24 = -qJD(6) * t94 - t321 * t73 + t326 * t74;
t542 = t24 / 0.2e1;
t23 = qJD(6) * t389 + t321 * t74 + t326 * t73;
t543 = t23 / 0.2e1;
t545 = Ifges(7,1) * t543 + Ifges(7,4) * t542 + Ifges(7,5) * t530;
t556 = mrSges(7,3) * t16 + t539;
t558 = -t17 * mrSges(7,3) + t540;
t560 = t11 * t327 - t12 * t322;
t572 = t231 * t213;
t573 = t228 * t213;
t620 = -t110 / 0.2e1;
t621 = -t24 / 0.2e1;
t622 = -t23 / 0.2e1;
t582 = Ifges(7,4) * t622 + Ifges(7,2) * t621 + Ifges(7,6) * t620;
t608 = t554 * t231;
t121 = -pkin(4) * t314 - t134;
t381 = t322 * mrSges(6,1) + t327 * mrSges(6,2);
t610 = t121 * t381;
t617 = Ifges(6,5) * t172 + Ifges(7,5) * t94 + Ifges(6,6) * t171 + Ifges(7,6) * t389 + Ifges(6,3) * t206 + Ifges(7,3) * t193;
t477 = t172 * Ifges(6,4);
t83 = t171 * Ifges(6,2) + t206 * Ifges(6,6) + t477;
t85 = -pkin(5) * t171 + t121;
t623 = (Ifges(6,5) * t215 + t213 * t376) * t525 + (Ifges(6,6) * t215 + t213 * t371) * t527 + (Ifges(6,3) * t215 + t213 * t366) * t521 + t558 * t608 - t572 * t540 - t573 * t539 - t213 * t610 + (t366 * t520 + t371 * t526 + t376 * t524 + t413 + t610) * qJD(5) + (Ifges(6,5) * t322 + Ifges(6,6) * t327) * t529 + (Ifges(6,2) * t327 + t491) * t537 + (Ifges(6,1) * t322 + t490) * t538 - (Ifges(5,1) * t213 - t204 + t617) * t215 / 0.2e1 - (-Ifges(5,2) * t215 + t137 + t203) * t213 / 0.2e1 - ((t162 - t573) * mrSges(7,2) + (t572 - t608) * mrSges(7,1)) * t85 - t16 * (mrSges(7,1) * t215 + mrSges(7,3) * t573) + t17 * (mrSges(7,2) * t215 + mrSges(7,3) * t572) + t228 * t582 + t63 * (mrSges(6,2) * t215 + mrSges(6,3) * t461) + t62 * (-mrSges(6,1) * t215 + mrSges(6,3) * t460) + t213 * t478 + t135 * t497 + (t461 / 0.2e1 + t398) * t83 + (-t2 * t228 - t231 * t3) * mrSges(7,3) - t84 * t460 / 0.2e1 + t556 * t162 + (-t432 * t62 - t433 * t63 + t560) * mrSges(6,3) + t231 * t545 + t322 * t541 + Ifges(5,3) * t313 - t314 * (Ifges(5,5) * t213 - Ifges(5,6) * t215) / 0.2e1 + t55 * t258 + t33 * (mrSges(7,1) * t228 + mrSges(7,2) * t231) - t183 * (mrSges(5,1) * t215 + mrSges(5,2) * t213) + t27 * t513 + t136 * t517 + Ifges(5,5) * t114 + Ifges(5,6) * t115 - t58 * mrSges(5,2) + t59 * mrSges(5,1);
t102 = mrSges(5,1) * t313 - mrSges(5,3) * t114;
t36 = -mrSges(6,1) * t74 + mrSges(6,2) * t73;
t613 = t102 - t36;
t437 = qJD(2) * t328;
t212 = -t322 * t437 + t327 * t441;
t214 = t322 * t441 + t327 * t437;
t612 = -t212 * t321 - t214 * t326 - t228 * t434 - t323 * t608;
t334 = t554 * t228;
t611 = -t212 * t326 + t214 * t321 - t231 * t434 + t323 * t334;
t571 = mrSges(5,1) * t314 + mrSges(6,1) * t171 - mrSges(6,2) * t172 - t497;
t607 = -qJD(4) * t250 - t328 * t242 + t323 * t629;
t295 = Ifges(3,4) * t440;
t488 = Ifges(4,5) * t329;
t377 = t324 * Ifges(4,1) - t488;
t606 = Ifges(3,1) * t441 + qJD(1) * t377 + qJD(2) * t619 + t295;
t415 = mrSges(4,2) * t441;
t605 = -mrSges(3,3) * t441 - t415 + (mrSges(3,1) + mrSges(4,1)) * qJD(2);
t604 = t324 * t618 + t329 * t619;
t385 = t329 * mrSges(4,1) + t324 * mrSges(4,3);
t387 = mrSges(3,1) * t329 - mrSges(3,2) * t324;
t603 = t385 + t387;
t165 = qJD(2) * t229 - t344;
t349 = t165 * t322 - t432 * t562;
t601 = (t433 - t461) * pkin(5);
t600 = -t224 * t329 + t225 * t324;
t192 = -qJDD(2) * pkin(2) + t397;
t599 = t184 * t329 + t192 * t324;
t598 = -mrSges(2,1) - t603;
t597 = -mrSges(3,3) - mrSges(4,2) + mrSges(5,3) + mrSges(2,2);
t390 = t229 * mrSges(5,1) - mrSges(5,2) * t562;
t592 = -t564 * t229 - t390;
t591 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t590 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t588 = Ifges(7,1) * t573 + Ifges(7,4) * t572 - Ifges(7,5) * t215;
t587 = Ifges(7,4) * t573 + Ifges(7,2) * t572 - Ifges(7,6) * t215;
t586 = Ifges(7,5) * t573 + Ifges(7,6) * t572 - Ifges(7,3) * t215;
t585 = m(5) * t134 - m(6) * t121 + t571;
t185 = t500 * t322;
t186 = t500 * t327;
t117 = t185 * t321 - t186 * t326;
t579 = -qJD(6) * t117 - t321 * t628 + t627 * t326;
t116 = t185 * t326 + t186 * t321;
t578 = qJD(6) * t116 + t627 * t321 + t326 * t628;
t263 = t331 * t322;
t265 = t331 * t327;
t175 = t263 * t321 - t265 * t326;
t577 = -qJD(6) * t175 - t321 * t626 + t326 * t625;
t173 = t263 * t326 + t265 * t321;
t576 = qJD(6) * t173 + t321 * t625 + t326 * t626;
t544 = m(7) * pkin(5);
t575 = mrSges(6,1) + t544;
t574 = -t135 + t601;
t570 = -t601 - t607;
t568 = -t161 + t189;
t306 = t324 * qJ(3);
t310 = t329 * pkin(2);
t443 = t310 + t306;
t412 = t329 * pkin(3) + t443;
t355 = pkin(9) * t562 + t412;
t132 = pkin(4) * t229 + pkin(1) + t355;
t531 = pkin(7) - pkin(8);
t264 = t531 * t324;
t266 = t531 * t329;
t176 = t264 * t323 + t266 * t328;
t167 = t327 * t176;
t81 = t322 * t132 + t167;
t566 = t328 * t264 - t266 * t323;
t311 = t330 * pkin(7);
t396 = -pkin(1) - t306;
t565 = t325 * (t329 * t332 + t396) - pkin(8) * t330 + t311;
t48 = mrSges(6,1) * t113 - mrSges(6,3) * t73;
t49 = -mrSges(6,2) * t113 + mrSges(6,3) * t74;
t561 = -t322 * t48 + t327 * t49;
t557 = g(1) * t330 + g(2) * t325;
t555 = qJD(3) + t297;
t496 = mrSges(6,3) * t171;
t119 = -mrSges(6,2) * t206 + t496;
t495 = mrSges(6,3) * t172;
t120 = mrSges(6,1) * t206 - t495;
t336 = (-t322 * t63 - t327 * t62) * qJD(5) + t560;
t548 = m(6) * t336 - t119 * t433 - t120 * t432 + t561;
t535 = -t389 / 0.2e1;
t534 = t389 / 0.2e1;
t533 = -t94 / 0.2e1;
t532 = t94 / 0.2e1;
t523 = -t193 / 0.2e1;
t522 = t193 / 0.2e1;
t514 = t324 / 0.2e1;
t510 = pkin(5) * t172;
t508 = pkin(5) * t322;
t507 = pkin(7) * t324;
t506 = pkin(7) * t329;
t503 = g(3) * t562;
t493 = Ifges(3,4) * t324;
t492 = Ifges(3,4) * t329;
t489 = Ifges(4,5) * t324;
t471 = t329 * mrSges(4,3);
t467 = t165 * t327;
t457 = t562 * t322;
t456 = t562 * t327;
t453 = t322 * t330;
t359 = t197 * t303 - t304 * t330;
t360 = -t197 * t304 - t303 * t330;
t446 = -t359 * mrSges(7,1) + t360 * mrSges(7,2);
t158 = t199 * t303 - t304 * t325;
t159 = -t199 * t304 - t303 * t325;
t445 = t158 * mrSges(7,1) - t159 * mrSges(7,2);
t436 = qJD(2) * t329;
t444 = qJ(3) * t436 + t302;
t442 = t330 * pkin(1) + t325 * pkin(7);
t439 = qJD(2) * t323;
t438 = qJD(2) * t324;
t424 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t110;
t423 = Ifges(6,5) * t73 + Ifges(6,6) * t74 + Ifges(6,3) * t113;
t418 = m(4) - t633;
t414 = mrSges(4,2) * t440;
t164 = t323 * t436 - t324 * t437 - t602;
t182 = t324 * t410 + t444;
t75 = pkin(4) * t164 - pkin(9) * t165 + t182;
t245 = t531 * t438;
t246 = qJD(2) * t266;
t98 = qJD(4) * t566 - t245 * t328 + t246 * t323;
t394 = -t322 * t98 + t327 * t75;
t392 = -t430 / 0.2e1;
t80 = t327 * t132 - t176 * t322;
t201 = -qJDD(2) * mrSges(4,1) + t248 * mrSges(4,2);
t388 = pkin(2) * t447 + qJ(3) * t450 + t442;
t239 = pkin(4) - t249;
t386 = mrSges(3,1) * t324 + mrSges(3,2) * t329;
t378 = -t303 * mrSges(7,1) - t304 * mrSges(7,2);
t375 = Ifges(7,1) * t162 + Ifges(7,4) * t608;
t374 = Ifges(7,1) * t231 - Ifges(7,4) * t228;
t373 = t329 * Ifges(3,2) + t493;
t370 = Ifges(7,4) * t162 + Ifges(7,2) * t608;
t369 = Ifges(7,4) * t231 - Ifges(7,2) * t228;
t365 = Ifges(7,5) * t162 + Ifges(7,6) * t608;
t364 = Ifges(7,5) * t231 - Ifges(7,6) * t228;
t64 = pkin(5) * t229 + pkin(10) * t456 + t80;
t66 = pkin(10) * t457 + t81;
t34 = -t321 * t66 + t326 * t64;
t35 = t321 * t64 + t326 * t66;
t363 = -t322 * t62 + t327 * t63;
t362 = pkin(3) * t447 + t388;
t358 = -t197 * t327 - t453;
t357 = t197 * t322 - t327 * t330;
t169 = t199 * t322 - t325 * t327;
t251 = -qJD(2) * pkin(2) + t555;
t255 = t298 + t317;
t356 = t251 * t329 - t255 * t324;
t353 = t424 + t591;
t350 = pkin(1) * t386;
t348 = -t433 * t562 - t467;
t347 = t217 * (t324 * mrSges(4,1) - t471);
t346 = t324 * (Ifges(3,1) * t329 - t493);
t345 = t329 * (Ifges(4,3) * t324 + t488);
t31 = t132 * t432 - t176 * t433 + t322 * t75 + t327 * t98;
t341 = -pkin(8) * t325 + t362;
t340 = t471 + (-m(4) * pkin(2) - mrSges(4,1)) * t324;
t99 = qJD(4) * t176 - t245 * t323 - t328 * t246;
t294 = Ifges(4,5) * t441;
t280 = qJ(3) * t447;
t279 = t329 * t325 * qJ(3);
t257 = qJD(2) * mrSges(4,3) + t414;
t256 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t440;
t252 = -pkin(1) - t443;
t234 = pkin(2) * t441 - t287;
t233 = t385 * qJD(1);
t221 = pkin(1) + t412;
t209 = Ifges(3,6) * qJD(2) + qJD(1) * t373;
t208 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t440 + t294;
t207 = t239 + t501;
t205 = pkin(2) * t438 - t444;
t202 = -mrSges(4,2) * t247 + qJDD(2) * mrSges(4,3);
t198 = t228 * t323;
t195 = t231 * t323;
t180 = -mrSges(5,2) * t314 + mrSges(5,3) * t213;
t170 = -t199 * t327 - t322 * t325;
t150 = t381 * t562;
t148 = -mrSges(5,1) * t213 + mrSges(5,2) * t215;
t144 = t228 * t562;
t143 = t231 * t562;
t133 = -pkin(5) * t457 - t566;
t103 = -mrSges(5,2) * t313 + mrSges(5,3) * t115;
t79 = mrSges(7,1) * t193 - mrSges(7,3) * t94;
t78 = -mrSges(7,2) * t193 + mrSges(7,3) * t389;
t60 = pkin(5) * t349 + t99;
t56 = -mrSges(7,1) * t389 + mrSges(7,2) * t94;
t41 = -t165 * t231 - t334 * t562;
t40 = -t165 * t228 + t562 * t608;
t32 = -qJD(5) * t81 + t394;
t20 = t326 * t50 - t475;
t19 = -t321 * t50 - t473;
t18 = -pkin(10) * t349 + t31;
t15 = -mrSges(7,2) * t110 + mrSges(7,3) * t24;
t14 = mrSges(7,1) * t110 - mrSges(7,3) * t23;
t13 = -pkin(10) * t467 + pkin(5) * t164 + (-t167 + (-pkin(10) * t562 - t132) * t322) * qJD(5) + t394;
t9 = -mrSges(7,1) * t24 + mrSges(7,2) * t23;
t5 = -qJD(6) * t35 + t13 * t326 - t18 * t321;
t4 = qJD(6) * t34 + t13 * t321 + t18 * t326;
t1 = [-(-m(5) * t59 + m(6) * t55 - t613) * t566 - (-t59 * mrSges(5,3) + Ifges(5,1) * t114 + Ifges(5,4) * t115 + Ifges(5,5) * t313 + t366 * t529 + t371 * t537 + t376 * t538 + t398 * t84) * t562 + (t324 * Ifges(3,1) + t377 + t492) * t248 / 0.2e1 + (t424 + t423) * t229 / 0.2e1 + t387 * t319 - t252 * mrSges(4,3) * t248 + (-t360 * mrSges(7,1) - t359 * mrSges(7,2) - (-pkin(5) * t453 - t197 * t291 + t565) * m(7) - t358 * mrSges(6,1) - t357 * mrSges(6,2) - (-pkin(4) * t197 + t565) * m(6) - m(5) * t565 + t197 * mrSges(5,1) + (-m(3) - m(4)) * t311 + (m(3) * pkin(1) - m(4) * (t396 - t310) - t598) * t325 - t636 * t196 + t597 * t330) * g(1) + (-m(4) * t388 - m(6) * (-pkin(4) * t199 + t341) - t170 * mrSges(6,1) - t169 * mrSges(6,2) - m(5) * t341 + t199 * mrSges(5,1) - m(7) * (-t199 * t291 + t362) - t159 * mrSges(7,1) - t158 * mrSges(7,2) - m(3) * t442 + t598 * t330 + t636 * t200 + (-m(7) * (-pkin(8) - t508) + t597) * t325) * g(2) + m(4) * (t147 * t252 + t205 * t217) + (Ifges(7,5) * t40 + Ifges(7,6) * t41) * t522 + (Ifges(7,5) * t144 + Ifges(7,6) * t143) * t530 + (-mrSges(5,3) * t58 - Ifges(5,4) * t114 + Ifges(6,5) * t538 + Ifges(7,5) * t543 - Ifges(5,2) * t115 - Ifges(5,6) * t313 + Ifges(6,6) * t537 + Ifges(7,6) * t542 + Ifges(6,3) * t529 + Ifges(7,3) * t530 + t590 + t591) * t229 - t585 * t99 + (-Ifges(6,5) * t348 - Ifges(6,6) * t349) * t520 + m(6) * (t11 * t81 + t12 * t80 + t31 * t63 + t32 * t62) + (t619 * t324 - t618 * t329) * qJDD(2) / 0.2e1 + ((Ifges(3,1) + Ifges(4,1)) * t248 + t619 * qJDD(2)) * t514 + (m(4) * (qJD(2) * t356 + t599) + m(3) * t600 - t605 * t436 + (-t257 - t256) * t438) * pkin(7) + (t251 * t436 - t255 * t438 + t599) * mrSges(4,2) + (t248 * t507 + t600) * mrSges(3,3) - t349 * t83 / 0.2e1 + (m(3) * pkin(1) ^ 2 + Ifges(2,3)) * qJDD(1) + t329 * (Ifges(3,4) * t248 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t543 + (Ifges(7,1) * t40 + Ifges(7,4) * t41) * t532 + (-Ifges(6,1) * t348 - Ifges(6,4) * t349) * t524 + (-qJDD(2) * mrSges(3,1) + t201) * t507 + (-t478 + t614 / 0.2e1 + t203 / 0.2e1 + Ifges(5,1) * t517 + t413 + t137 / 0.2e1 + t183 * mrSges(5,2)) * t165 - t143 * t582 + (t347 + t604 * qJD(2) / 0.2e1) * qJD(2) + (t11 * t457 + t12 * t456 + t348 * t62 - t349 * t63) * mrSges(6,3) + (t324 * (Ifges(4,1) * t329 + t489) + t329 * (-Ifges(3,2) * t324 + t492) + t346) * t430 / 0.2e1 - t329 * (Ifges(4,5) * t248 + Ifges(4,6) * qJDD(2)) / 0.2e1 + m(5) * (t112 * t221 + t135 * t98 + t176 * t58 + t182 * t183) + (-Ifges(6,4) * t348 - Ifges(6,2) * t349) * t526 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) - t615 / 0.2e1 - t616 / 0.2e1 - Ifges(5,4) * t517 + Ifges(6,3) * t520 + Ifges(7,3) * t522 + Ifges(6,5) * t524 + Ifges(6,6) * t526 + Ifges(7,5) * t532 + Ifges(7,6) * t534 + t16 * mrSges(7,1) - t17 * mrSges(7,2) - t136 / 0.2e1 + t183 * mrSges(5,1) - mrSges(5,3) * t135 + t617 / 0.2e1) * t164 + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t542 + (Ifges(7,4) * t40 + Ifges(7,2) * t41) * t534 + t202 * t506 + (t143 * t2 - t144 * t3 - t16 * t40 + t17 * t41) * mrSges(7,3) - t147 * t385 + t121 * (mrSges(6,1) * t349 - mrSges(6,2) * t348) + t606 * t436 / 0.2e1 + m(7) * (t133 * t33 + t16 * t5 + t17 * t4 + t2 * t35 + t3 * t34 + t60 * t85) + t34 * t14 + t35 * t15 + (-pkin(1) * t248 - qJDD(2) * t506) * mrSges(3,2) + (-t209 / 0.2e1 + t208 / 0.2e1) * t438 + t27 * t457 / 0.2e1 - t350 * t430 + t144 * t545 - t456 * t541 + t112 * t390 - t205 * t233 + t221 * (-mrSges(5,1) * t115 + mrSges(5,2) * t114) + t4 * t78 + t5 * t79 + t80 * t48 + t81 * t49 + t85 * (-mrSges(7,1) * t41 + mrSges(7,2) * t40) + (-mrSges(3,3) * t506 - t373 / 0.2e1 + t489 / 0.2e1 - pkin(1) * mrSges(3,1) + t252 * mrSges(4,1) + (Ifges(4,5) - Ifges(3,4)) * t514 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t329) * t247 + t31 * t119 + t32 * t120 + t133 * t9 + t33 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + t345 * t392 - t55 * t150 + t41 * t46 / 0.2e1 + t40 * t47 / 0.2e1 + t60 * t56 + t176 * t103 + t98 * t180 + t182 * t148; (-m(4) * t443 - m(6) * t355 + (-m(7) - m(5)) * t412 - t550 * t562 + t592 - t603) * g(3) + t364 * t620 + t369 * t621 + t374 * t622 - (Ifges(4,1) * t440 + t208 + t294) * t441 / 0.2e1 + (-pkin(2) * t192 + qJ(3) * t184 + qJD(3) * t255 - t217 * t234) * m(4) + (-pkin(7) * t356 * m(4) - t347 + (-t346 / 0.2e1 + t345 / 0.2e1 + t350) * qJD(1)) * qJD(1) + t586 * t523 + t587 * t535 + t588 * t533 + (t135 * t568 - t183 * t194 + t249 * t59 + t250 * t58) * m(5) + (t189 * t363 + t239 * t55 - t62 * t67 - t63 * t68) * m(6) + t618 * t247 + t619 * t248 + t630 * t120 + t631 * t119 + (-m(4) * t280 - t330 * t340 - m(5) * (t330 * t416 + t280) + t580 * (t332 * t450 + t280) - t634) * g(1) + (-m(4) * t279 - t325 * t340 + t633 * (t325 * t416 + t279) - t635) * g(2) + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) - t623 + t604 * t392 + t605 * t298 - (-Ifges(3,2) * t441 + t295 + t606) * t440 / 0.2e1 + t585 * t607 + t579 * t79 + (t116 * t3 + t117 * t2 + t16 * t579 + t17 * t578 + t207 * t33 + t570 * t85) * m(7) + t578 * t78 + t209 * t441 / 0.2e1 + t568 * t180 + t570 * t56 + t557 * t386 + t555 * t257 - t251 * t414 + t548 * t240 + t249 * t102 + t250 * t103 + t239 * t36 + t234 * t233 + t224 * mrSges(3,2) - t225 * mrSges(3,1) + t207 * t9 - pkin(2) * t201 + qJ(3) * t202 + t365 * t522 + t255 * t415 + t256 * t297 + t116 * t14 + t117 * t15 + t375 * t532 + t370 * t534 + t184 * mrSges(4,3) - t192 * mrSges(4,1) - t194 * t148; -qJD(2) * t257 - t214 * t119 - t212 * t120 - t195 * t14 - t198 * t15 + t611 * t79 + t612 * t78 + t418 * t329 * g(3) + (-qJD(2) * t180 - t9 + (t119 * t327 - t120 * t322 + t180) * qJD(4) + t613) * t328 + ((-t148 - t233) * qJD(1) - t557 * t418) * t324 + (t103 + (-t119 * t322 - t120 * t327) * qJD(5) + t314 * (t56 - t571) + t561) * t323 + t201 + (-t195 * t3 - t198 * t2 - t328 * t33 + (t435 - t439) * t85 + t612 * t17 + t611 * t16) * m(7) + ((qJD(4) * t363 - t55) * t328 + (qJD(4) * t121 + t336) * t323 - t121 * t439 - t212 * t62 - t214 * t63) * m(6) + (-t183 * t441 + t323 * t58 + t328 * t59 + t314 * (-t134 * t323 + t135 * t328)) * m(5) + (-qJD(2) * t255 + t217 * t441 + t192) * m(4); (-t562 * t563 - t592) * g(3) + (t365 - t586) * t523 + (t370 - t587) * t535 + (t375 - t588) * t533 + (-pkin(4) * t55 - t121 * t135 - t62 * t76 - t63 * t77) * m(6) + t634 * g(1) + t635 * g(2) + t623 - pkin(4) * t36 + t576 * t78 + t577 * t79 + (t16 * t577 + t17 * t576 + t173 * t3 + t175 * t2 - t291 * t33 + t574 * t85) * m(7) + t574 * t56 + t571 * t135 + t369 * t542 + t374 * t543 + t548 * pkin(9) - t291 * t9 - t77 * t119 - t76 * t120 + t364 * t530 + t173 * t14 + t175 * t15 - t134 * t180; (-Ifges(6,2) * t172 + t168 + t84) * t527 + (t496 - t119) * t62 + (t495 + t120) * t63 + t590 + t353 + t423 + ((-qJD(6) * t79 + t15) * t321 + t14 * t326 + t431 * t78) * pkin(5) + (-m(7) * t508 + t378) * t503 - t56 * t510 - m(7) * (t16 * t19 + t17 * t20 + t510 * t85) + (-mrSges(6,2) * t358 + t357 * t575 - t446) * g(2) + (mrSges(6,2) * t170 - t169 * t575 - t445) * g(1) + (-t85 * mrSges(7,2) + Ifges(7,1) * t533 + Ifges(7,4) * t535 + Ifges(7,5) * t523 + t556) * t389 - (t85 * mrSges(7,1) + Ifges(7,4) * t533 + Ifges(7,2) * t535 + Ifges(7,6) * t523 + t558) * t94 + (t2 * t321 + t3 * t326 + (-t16 * t321 + t17 * t326) * qJD(6)) * t544 + (Ifges(6,5) * t171 - Ifges(6,6) * t172) * t521 - t20 * t78 - t19 * t79 + t83 * t524 + (Ifges(6,1) * t171 - t477) * t525 - g(3) * t150 - t121 * (mrSges(6,1) * t172 + mrSges(6,2) * t171); -t85 * (mrSges(7,1) * t94 + mrSges(7,2) * t389) + (Ifges(7,1) * t389 - t511) * t533 + t46 * t532 + (Ifges(7,5) * t389 - Ifges(7,6) * t94) * t523 - t16 * t78 + t17 * t79 - g(1) * t445 - g(2) * t446 + t378 * t503 + (t16 * t389 + t17 * t94) * mrSges(7,3) + t353 + (-Ifges(7,2) * t94 + t47 + t87) * t535;];
tau  = t1;
