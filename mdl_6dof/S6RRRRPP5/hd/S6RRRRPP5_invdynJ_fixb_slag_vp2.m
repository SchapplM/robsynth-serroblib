% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:17
% EndTime: 2019-03-09 21:06:30
% DurationCPUTime: 44.13s
% Computational Cost: add. (11488->827), mult. (25056->1001), div. (0->0), fcn. (17051->10), ass. (0->395)
t624 = -mrSges(7,2) - mrSges(6,3);
t557 = -mrSges(5,2) - t624;
t322 = sin(qJ(2));
t485 = Ifges(3,4) * t322;
t326 = cos(qJ(2));
t506 = t326 / 0.2e1;
t630 = Ifges(3,2) * t506 + t485 / 0.2e1;
t625 = mrSges(5,1) + mrSges(6,1);
t589 = Ifges(7,4) + Ifges(6,5);
t555 = -Ifges(5,4) + t589;
t636 = t555 + t589;
t328 = -pkin(9) - pkin(8);
t383 = m(7) * (-qJ(6) - t328) - mrSges(7,3);
t623 = -mrSges(5,3) - mrSges(6,2);
t559 = t623 * t322;
t635 = -t383 * t322 + t559;
t440 = qJD(1) * t322;
t308 = pkin(7) * t440;
t271 = -qJD(2) * pkin(2) + t308;
t321 = sin(qJ(3));
t410 = t321 * t440;
t325 = cos(qJ(3));
t437 = qJD(2) * t325;
t359 = t410 - t437;
t182 = pkin(3) * t359 + t271;
t439 = qJD(1) * t326;
t292 = qJD(3) - t439;
t282 = -qJD(4) - t292;
t532 = pkin(4) + pkin(5);
t414 = t325 * t440;
t249 = qJD(2) * t321 + t414;
t320 = sin(qJ(4));
t324 = cos(qJ(4));
t338 = t324 * t249 - t320 * t359;
t573 = qJ(6) * t338;
t388 = pkin(2) * t326 + pkin(8) * t322;
t265 = -pkin(1) - t388;
t238 = t265 * qJD(1);
t309 = pkin(7) * t439;
t272 = qJD(2) * pkin(8) + t309;
t167 = t325 * t238 - t272 * t321;
t129 = -pkin(9) * t249 + t167;
t115 = pkin(3) * t292 + t129;
t168 = t321 * t238 + t325 * t272;
t130 = -pkin(9) * t359 + t168;
t463 = t320 * t130;
t60 = t324 * t115 - t463;
t32 = t60 + t573;
t571 = qJD(5) - t32;
t29 = t282 * t532 + t571;
t161 = t249 * t320 + t324 * t359;
t334 = -qJ(5) * t338 + t182;
t43 = -t161 * t532 + qJD(6) - t334;
t48 = pkin(4) * t282 + qJD(5) - t60;
t156 = Ifges(5,4) * t161;
t556 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t599 = t589 * t161;
t590 = Ifges(6,4) + Ifges(5,5);
t629 = -Ifges(7,5) + t590;
t552 = -t282 * t629 + t338 * t556 - t156 + t599;
t69 = t161 * pkin(4) + t334;
t634 = mrSges(5,2) * t182 + mrSges(6,2) * t48 + mrSges(7,2) * t43 - mrSges(5,3) * t60 - mrSges(6,3) * t69 - mrSges(7,3) * t29 + t552 / 0.2e1;
t264 = t282 * qJ(5);
t602 = qJ(6) * t161;
t453 = t324 * t130;
t61 = t320 * t115 + t453;
t33 = t61 + t602;
t30 = -t264 + t33;
t49 = -t264 + t61;
t480 = Ifges(5,4) * t338;
t85 = -Ifges(5,2) * t161 - Ifges(5,6) * t282 + t480;
t633 = -mrSges(5,1) * t182 - mrSges(6,1) * t69 + mrSges(7,1) * t43 + mrSges(6,2) * t49 + mrSges(5,3) * t61 - mrSges(7,3) * t30 + t85 / 0.2e1;
t430 = qJD(1) * qJD(2);
t260 = qJDD(1) * t322 + t326 * t430;
t346 = t359 * qJD(3);
t151 = qJDD(2) * t321 + t260 * t325 - t346;
t152 = -qJD(3) * t249 + qJDD(2) * t325 - t260 * t321;
t57 = -qJD(4) * t161 + t324 * t151 + t320 * t152;
t538 = t57 / 0.2e1;
t58 = qJD(4) * t338 + t320 * t151 - t324 * t152;
t536 = t58 / 0.2e1;
t526 = -t161 / 0.2e1;
t525 = t161 / 0.2e1;
t259 = t326 * qJDD(1) - t322 * t430;
t246 = qJDD(3) - t259;
t237 = qJDD(4) + t246;
t514 = t237 / 0.2e1;
t632 = -t249 / 0.2e1;
t508 = t282 / 0.2e1;
t631 = -t292 / 0.2e1;
t523 = -t338 / 0.2e1;
t522 = t338 / 0.2e1;
t593 = t359 / 0.2e1;
t250 = t320 * t321 - t324 * t325;
t558 = qJD(3) + qJD(4);
t172 = t558 * t250;
t356 = t250 * t326;
t199 = qJD(1) * t356;
t628 = t172 - t199;
t509 = -t282 / 0.2e1;
t627 = Ifges(5,4) * t526 + Ifges(7,5) * t508 + t509 * t590 + t522 * t556 + t525 * t589 + t634;
t585 = Ifges(6,6) - Ifges(7,6);
t587 = Ifges(7,2) + Ifges(6,3);
t601 = t589 * t338;
t579 = t161 * t587 - t282 * t585 + t601;
t586 = -Ifges(5,6) + Ifges(6,6);
t626 = -Ifges(5,2) * t526 + Ifges(7,6) * t508 + t509 * t586 + t522 * t555 + t525 * t587 + t579 / 0.2e1 - t633;
t244 = t259 * pkin(7);
t588 = Ifges(6,2) + Ifges(5,3);
t431 = qJD(4) * t324;
t67 = t324 * t129 - t463;
t622 = pkin(3) * t431 - t67;
t251 = t320 * t325 + t321 * t324;
t173 = t558 * t251;
t198 = t251 * t439;
t621 = t198 - t173;
t415 = t321 * t439;
t435 = qJD(3) * t321;
t566 = -t309 + (-t415 + t435) * pkin(3);
t620 = Ifges(5,4) * t525 + Ifges(7,5) * t509 + t590 * t508 + t556 * t523 + t589 * t526 - t634;
t619 = -Ifges(5,2) * t525 + Ifges(7,6) * t509 + t586 * t508 + t555 * t523 + t587 * t526 + t633;
t578 = mrSges(3,2) * t326;
t381 = mrSges(3,1) * t322 + t578;
t618 = pkin(1) * t381 - t322 * (Ifges(3,1) * t326 - t485) / 0.2e1;
t617 = mrSges(7,1) + t625;
t319 = qJ(3) + qJ(4);
t312 = sin(t319);
t313 = cos(t319);
t380 = -mrSges(4,1) * t325 + mrSges(4,2) * t321;
t615 = m(4) * pkin(2) + t312 * t557 + t625 * t313 - t380;
t468 = qJDD(1) * pkin(1);
t171 = -pkin(2) * t259 - pkin(8) * t260 - t468;
t216 = qJDD(2) * pkin(8) + t244;
t433 = qJD(3) * t325;
t73 = t321 * t171 + t325 * t216 + t238 * t433 - t272 * t435;
t74 = -qJD(3) * t168 + t325 * t171 - t216 * t321;
t614 = -t74 * mrSges(4,1) + t73 * mrSges(4,2);
t34 = pkin(3) * t246 - pkin(9) * t151 + t74;
t432 = qJD(4) * t320;
t46 = pkin(9) * t152 + t73;
t10 = -t115 * t432 - t130 * t431 - t320 * t46 + t324 * t34;
t357 = qJDD(5) - t10;
t1 = -qJ(6) * t57 - qJD(6) * t338 - t237 * t532 + t357;
t9 = t115 * t431 - t130 * t432 + t320 * t34 + t324 * t46;
t6 = t237 * qJ(5) - t282 * qJD(5) + t9;
t3 = qJ(6) * t58 + qJD(6) * t161 + t6;
t7 = -pkin(4) * t237 + t357;
t613 = -t10 * mrSges(5,1) + t7 * mrSges(6,1) + t1 * mrSges(7,1) + t9 * mrSges(5,2) - t3 * mrSges(7,2) - t6 * mrSges(6,3);
t610 = t29 * mrSges(7,1) + t48 * mrSges(6,1) + t61 * mrSges(5,2) + Ifges(4,5) * t632 + Ifges(4,6) * t593 + Ifges(4,3) * t631 + t586 * t526 + t590 * t523 + Ifges(7,5) * t522 + Ifges(3,6) * qJD(2) / 0.2e1 + Ifges(7,6) * t525 + qJD(1) * t630 - t30 * mrSges(7,2) - t49 * mrSges(6,3) - t60 * mrSges(5,1) + (t588 + Ifges(7,3)) * t508;
t245 = t260 * pkin(7);
t386 = qJDD(2) * pkin(2) - t245;
t354 = pkin(3) * t152 + t386;
t336 = qJ(5) * t57 + qJD(5) * t338 + t354;
t11 = pkin(4) * t58 - t336;
t4 = -t532 * t58 + qJDD(6) + t336;
t515 = -t237 / 0.2e1;
t537 = -t58 / 0.2e1;
t609 = -mrSges(5,2) * t354 + mrSges(6,2) * t7 + mrSges(7,2) * t4 - mrSges(5,3) * t10 - mrSges(6,3) * t11 - mrSges(7,3) * t1 + Ifges(5,4) * t537 + Ifges(7,5) * t515 + 0.2e1 * t538 * t556 + t636 * t536 + (t590 + t629) * t514;
t608 = -mrSges(5,1) * t354 + mrSges(6,1) * t11 - mrSges(7,1) * t4 - mrSges(6,2) * t6 - mrSges(5,3) * t9 + mrSges(7,3) * t3 + 0.2e1 * t536 * t587 - t57 * Ifges(5,4) / 0.2e1 + t636 * t538 + (Ifges(7,6) + Ifges(5,6)) * t515 + (t586 + t585) * t514 + (-t537 + t536) * Ifges(5,2);
t416 = qJD(3) * t328;
t256 = t321 * t416;
t257 = t325 * t416;
t273 = t328 * t321;
t274 = t328 * t325;
t111 = t324 * t256 + t320 * t257 + t273 * t431 + t274 * t432;
t387 = pkin(2) * t322 - pkin(8) * t326;
t255 = t387 * qJD(1);
t178 = pkin(7) * t410 + t325 * t255;
t451 = t325 * t326;
t365 = pkin(3) * t322 - pkin(9) * t451;
t144 = qJD(1) * t365 + t178;
t230 = t321 * t255;
t458 = t322 * t325;
t460 = t321 * t326;
t164 = t230 + (-pkin(7) * t458 - pkin(9) * t460) * qJD(1);
t90 = t320 * t144 + t324 * t164;
t77 = qJ(5) * t440 + t90;
t604 = t111 - t77;
t177 = t320 * t273 - t324 * t274;
t112 = qJD(4) * t177 + t256 * t320 - t324 * t257;
t89 = t144 * t324 - t320 * t164;
t603 = -t112 - t89;
t574 = qJD(5) + t622;
t470 = qJ(5) * t161;
t600 = qJ(5) * t628 - qJD(5) * t251 + t566;
t436 = qJD(2) * t326;
t409 = t321 * t436;
t350 = t322 * t433 + t409;
t323 = sin(qJ(1));
t327 = cos(qJ(1));
t450 = t326 * t327;
t227 = -t321 * t450 + t323 * t325;
t598 = g(1) * t327 + g(2) * t323;
t595 = -m(3) - m(4);
t594 = -m(6) - m(7);
t529 = t151 / 0.2e1;
t528 = t152 / 0.2e1;
t513 = t246 / 0.2e1;
t592 = pkin(4) * t338;
t591 = -mrSges(3,3) + mrSges(2,2);
t40 = mrSges(7,2) * t237 + mrSges(7,3) * t58;
t42 = -mrSges(6,2) * t58 + mrSges(6,3) * t237;
t583 = t42 + t40;
t582 = -qJ(6) * t621 + qJD(6) * t250 + t604;
t421 = t532 * t322;
t581 = qJ(6) * t628 + qJD(1) * t421 - qJD(6) * t251 - t603;
t580 = t532 * t621 - t600;
t577 = -pkin(4) * t621 + t600;
t575 = -t573 + t574;
t213 = t251 * t322;
t572 = t338 * t532;
t491 = mrSges(6,2) * t161;
t131 = -mrSges(6,3) * t282 - t491;
t487 = mrSges(7,3) * t161;
t132 = -mrSges(7,2) * t282 + t487;
t570 = t131 + t132;
t489 = mrSges(5,3) * t161;
t133 = mrSges(5,2) * t282 - t489;
t569 = t131 + t133;
t488 = mrSges(5,3) * t338;
t135 = -mrSges(5,1) * t282 - t488;
t490 = mrSges(6,2) * t338;
t136 = mrSges(6,1) * t282 + t490;
t568 = t135 - t136;
t242 = Ifges(4,4) * t359;
t143 = t249 * Ifges(4,1) + t292 * Ifges(4,5) - t242;
t307 = Ifges(3,4) * t439;
t567 = Ifges(3,1) * t440 + Ifges(3,5) * qJD(2) + t325 * t143 + t307;
t565 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t359 + t249 * mrSges(4,2) + mrSges(3,3) * t440;
t316 = t326 * pkin(4);
t469 = qJ(5) * t326;
t564 = t312 * t469 + t313 * t316;
t296 = pkin(7) * t451;
t185 = t321 * t265 + t296;
t391 = -m(7) * t532 - mrSges(7,1);
t464 = t313 * t322;
t465 = t312 * t322;
t563 = (mrSges(6,1) - t391) * t465 + t624 * t464;
t434 = qJD(3) * t322;
t351 = -t321 * t434 + t325 * t436;
t562 = t237 * t588 + t57 * t590 + t58 * t586;
t561 = t244 * t326 + t245 * t322;
t560 = -t321 * t74 + t325 * t73;
t382 = t326 * mrSges(3,1) - t322 * mrSges(3,2);
t550 = t322 * mrSges(4,3) + mrSges(2,1) + t382;
t549 = m(7) * pkin(5) + mrSges(7,1);
t449 = t327 * t312;
t210 = -t323 * t313 + t326 * t449;
t211 = t312 * t323 + t313 * t450;
t548 = t210 * t617 - t211 * t557;
t454 = t323 * t326;
t208 = t312 * t454 + t313 * t327;
t209 = t313 * t454 - t449;
t547 = t208 * t617 - t209 * t557;
t546 = -m(6) * t48 + t568;
t545 = t549 + t625;
t540 = m(5) * pkin(3);
t535 = Ifges(4,1) * t529 + Ifges(4,4) * t528 + Ifges(4,5) * t513;
t512 = t249 / 0.2e1;
t499 = pkin(3) * t249;
t498 = pkin(3) * t320;
t497 = pkin(3) * t324;
t494 = g(3) * t322;
t314 = t322 * pkin(7);
t493 = -qJD(1) / 0.2e1;
t492 = qJD(3) / 0.2e1;
t486 = mrSges(7,3) * t338;
t484 = Ifges(3,4) * t326;
t483 = Ifges(4,4) * t249;
t482 = Ifges(4,4) * t321;
t481 = Ifges(4,4) * t325;
t477 = t168 * mrSges(4,3);
t462 = t321 * t322;
t461 = t321 * t323;
t459 = t321 * t327;
t456 = t322 * t328;
t304 = pkin(3) * t325 + pkin(2);
t277 = t326 * t304;
t248 = t325 * t265;
t166 = -pkin(9) * t458 + t248 + (-pkin(7) * t321 - pkin(3)) * t326;
t175 = -pkin(9) * t462 + t185;
t102 = t320 * t166 + t324 * t175;
t258 = t387 * qJD(2);
t438 = qJD(2) * t322;
t423 = pkin(7) * t438;
t442 = t325 * t258 + t321 * t423;
t293 = pkin(3) * t462;
t261 = t314 + t293;
t441 = t327 * pkin(1) + t323 * pkin(7);
t310 = pkin(7) * t436;
t422 = m(4) * pkin(8) + mrSges(4,3);
t419 = t327 * t456;
t418 = Ifges(4,5) * t151 + Ifges(4,6) * t152 + Ifges(4,3) * t246;
t317 = t327 * pkin(7);
t417 = pkin(3) * t459 + t323 * t456 + t317;
t183 = pkin(3) * t350 + t310;
t303 = -pkin(4) - t497;
t142 = -Ifges(4,2) * t359 + Ifges(4,6) * t292 + t483;
t408 = -t321 * t142 / 0.2e1;
t21 = -t58 * mrSges(7,1) + t57 * mrSges(7,2);
t39 = -t237 * mrSges(6,1) + t57 * mrSges(6,2);
t37 = -t237 * mrSges(7,1) - t57 * mrSges(7,3);
t396 = -t208 * pkin(4) + qJ(5) * t209;
t395 = -t210 * pkin(4) + qJ(5) * t211;
t394 = -qJ(5) * t312 - t304;
t66 = t129 * t320 + t453;
t101 = t166 * t324 - t320 * t175;
t176 = -t324 * t273 - t274 * t320;
t393 = t277 - t456;
t392 = pkin(3) * t461 + t304 * t450 + t441;
t275 = qJ(5) * t464;
t390 = -pkin(4) * t465 + t275;
t93 = t102 - t469;
t214 = t250 * t322;
t384 = -qJ(5) * t214 - t261;
t94 = -t101 + t316;
t379 = mrSges(4,1) * t321 + mrSges(4,2) * t325;
t377 = -mrSges(5,1) * t312 - mrSges(5,2) * t313;
t375 = Ifges(4,1) * t325 - t482;
t374 = Ifges(4,1) * t321 + t481;
t372 = -Ifges(4,2) * t321 + t481;
t371 = Ifges(4,2) * t325 + t482;
t370 = Ifges(3,5) * t326 - Ifges(3,6) * t322;
t369 = Ifges(4,5) * t325 - t321 * Ifges(4,6);
t368 = Ifges(4,5) * t321 + t325 * Ifges(4,6);
t367 = -t499 - t470;
t366 = t227 * pkin(3);
t364 = qJ(5) * t251 + t304;
t100 = t365 * qJD(2) + (-t296 + (pkin(9) * t322 - t265) * t321) * qJD(3) + t442;
t127 = t321 * t258 + t265 * t433 + (-t322 * t437 - t326 * t435) * pkin(7);
t105 = -pkin(9) * t350 + t127;
t27 = t100 * t324 - t320 * t105 - t166 * t432 - t175 * t431;
t225 = t321 * t454 + t325 * t327;
t361 = t271 * t379;
t26 = t320 * t100 + t324 * t105 + t166 * t431 - t175 * t432;
t353 = Ifges(7,5) * t57 + Ifges(7,6) * t58 - Ifges(7,3) * t237;
t352 = t225 * pkin(3);
t349 = t211 * pkin(4) + qJ(5) * t210 + t392;
t348 = t359 * mrSges(4,3);
t347 = -g(1) * t210 - g(2) * t208 - g(3) * t465;
t106 = -qJD(2) * t356 - t213 * t558;
t344 = qJ(5) * t106 - qJD(5) * t214 - t183;
t342 = t366 + t395;
t341 = Ifges(4,5) * t322 + t326 * t375;
t340 = Ifges(4,6) * t322 + t326 * t372;
t339 = Ifges(4,3) * t322 + t326 * t369;
t24 = qJ(5) * t438 - qJD(5) * t326 + t26;
t337 = -t352 + t396;
t331 = -t353 + t562 - t613;
t299 = qJ(5) + t498;
t298 = -pkin(5) + t303;
t269 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t439;
t239 = t379 * t322;
t228 = t325 * t450 + t461;
t226 = -t323 * t451 + t459;
t203 = t210 * pkin(5);
t200 = t208 * pkin(5);
t184 = -pkin(7) * t460 + t248;
t181 = mrSges(4,1) * t292 - mrSges(4,3) * t249;
t180 = -t292 * mrSges(4,2) - t348;
t179 = -pkin(7) * t414 + t230;
t157 = pkin(4) * t250 - t364;
t139 = qJ(6) * t250 + t177;
t138 = -qJ(6) * t251 + t176;
t134 = mrSges(7,1) * t282 - t486;
t128 = -qJD(3) * t185 + t442;
t126 = -t250 * t532 + t364;
t119 = pkin(4) * t213 - t384;
t118 = -mrSges(4,2) * t246 + mrSges(4,3) * t152;
t117 = mrSges(4,1) * t246 - mrSges(4,3) * t151;
t107 = -t432 * t462 + (t458 * t558 + t409) * t324 + t351 * t320;
t99 = mrSges(5,1) * t161 + mrSges(5,2) * t338;
t98 = -mrSges(7,1) * t161 + mrSges(7,2) * t338;
t97 = mrSges(6,1) * t161 - mrSges(6,3) * t338;
t96 = t470 + t592;
t95 = -t213 * t532 + t384;
t91 = -mrSges(4,1) * t152 + mrSges(4,2) * t151;
t78 = -pkin(4) * t440 - t89;
t76 = -t367 + t592;
t71 = t151 * Ifges(4,4) + t152 * Ifges(4,2) + t246 * Ifges(4,6);
t70 = qJ(6) * t213 + t93;
t68 = pkin(5) * t326 + qJ(6) * t214 + t94;
t65 = -t470 - t572;
t47 = t367 - t572;
t41 = -mrSges(5,2) * t237 - mrSges(5,3) * t58;
t38 = mrSges(5,1) * t237 - mrSges(5,3) * t57;
t35 = t66 + t602;
t28 = pkin(4) * t107 - t344;
t25 = -pkin(4) * t438 - t27;
t23 = -t107 * t532 + t344;
t22 = mrSges(5,1) * t58 + mrSges(5,2) * t57;
t20 = mrSges(6,1) * t58 - mrSges(6,3) * t57;
t13 = qJ(6) * t107 + qJD(6) * t213 + t24;
t12 = -qJ(6) * t106 - qJD(2) * t421 + qJD(6) * t214 - t27;
t2 = [-(t142 * t325 + t143 * t321) * t434 / 0.2e1 + t565 * t310 + t567 * t436 / 0.2e1 + (-m(5) * (t392 - t419) - m(6) * (t349 - t419) - m(7) * t349 - t228 * mrSges(4,1) - t227 * mrSges(4,2) + t595 * t441 + t591 * t323 - t545 * t211 - t557 * t210 + (-m(4) * t388 - t550 + t635) * t327) * g(2) + (mrSges(3,3) * t314 + t322 * Ifges(3,1) + t484 / 0.2e1 - pkin(1) * mrSges(3,2) + Ifges(3,4) * t506) * t260 + t608 * t213 - t609 * t214 + (-t167 * t351 - t168 * t350 - t458 * t74 - t462 * t73) * mrSges(4,3) + (-m(5) * t417 - t226 * mrSges(4,1) - t225 * mrSges(4,2) + t594 * (-t209 * pkin(4) - qJ(5) * t208 + t417) + t591 * t327 + t595 * t317 + t545 * t209 + t557 * t208 + (m(3) * pkin(1) - m(4) * t265 + (-m(7) * qJ(6) - mrSges(7,3)) * t322 + (-m(5) + t594) * (-pkin(1) - t277) + t550 - t559) * t323) * g(1) + (t167 * mrSges(4,1) - t168 * mrSges(4,2) + Ifges(5,6) * t526 - Ifges(7,3) * t508 + t588 * t509 + t522 * t629 + t585 * t525 - t610) * t438 + (mrSges(3,3) * t244 + (-Ifges(3,2) * t322 + t484) * t430 / 0.2e1 - t629 * t538 - Ifges(4,6) * t528 - Ifges(4,5) * t529 - Ifges(5,6) * t537 + Ifges(7,3) * t515 - t585 * t536 - t588 * t514 - Ifges(4,3) * t513 + t613 - t418 / 0.2e1 - t562 / 0.2e1 + t614) * t326 + t322 * t372 * t528 + t322 * t375 * t529 + t322 * t369 * t513 - t618 * t430 + t626 * t107 + t627 * t106 - t71 * t462 / 0.2e1 + t561 * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t561) - t269 * t423 - t359 * (qJD(2) * t340 - t371 * t434) / 0.2e1 + m(7) * (t1 * t68 + t12 * t29 + t13 * t30 + t23 * t43 + t3 * t70 + t4 * t95) + m(6) * (t11 * t119 + t24 * t49 + t25 * t48 + t28 * t69 + t6 * t93 + t7 * t94) + t292 * (qJD(2) * t339 - t368 * t434) / 0.2e1 + t408 * t436 + (pkin(1) * mrSges(3,1) + 0.2e1 * t630) * t259 + Ifges(2,3) * qJDD(1) + t261 * t22 + t458 * t535 + t184 * t117 + t185 * t118 + t127 * t180 + t128 * t181 + t183 * t99 + t12 * t134 + t27 * t135 + t25 * t136 + t24 * t131 + t13 * t132 + t26 * t133 + t119 * t20 + qJD(2) ^ 2 * t370 / 0.2e1 + t101 * t38 + t102 * t41 + t93 * t42 + t94 * t39 + t95 * t21 + t28 * t97 + t23 * t98 + t68 * t37 + t70 * t40 + (-mrSges(3,1) * t314 + Ifges(3,5) * t322 + 0.2e1 * Ifges(3,6) * t506 - pkin(7) * t578) * qJDD(2) + t271 * (mrSges(4,1) * t350 + mrSges(4,2) * t351) + t91 * t314 + t353 * t506 + (qJD(2) * t341 - t374 * t434) * t512 + m(5) * (t10 * t101 + t102 * t9 + t182 * t183 + t26 * t61 - t261 * t354 + t27 * t60) - t386 * t239 + m(4) * (t127 * t168 + t128 * t167 + t184 * t74 + t185 * t73 + (t271 * t436 - t322 * t386) * pkin(7)) + t382 * t468; (t41 + t42) * t177 - t565 * t309 + t566 * t99 - t568 * t112 + t569 * t111 + (-t322 * t422 - t382 - m(5) * t393 - m(6) * (t393 + t564) - m(7) * (t277 + t564) + (-t313 * t549 - t615) * t326 + t635) * g(3) + t608 * t250 + t609 * t251 + (t339 * t493 + t369 * t492) * t292 + t598 * (t381 + (-t383 - t422 + (m(5) + m(6)) * t328 + t623) * t326 + (m(5) * t304 - m(6) * (-pkin(4) * t313 + t394) - m(7) * t394 - t313 * t391 + t615) * t322) + t577 * t97 + (Ifges(3,2) * t439 / 0.2e1 + t629 * t523 + t585 * t526 + t588 * t508 - Ifges(7,3) * t509 + Ifges(5,6) * t525 + t610) * t440 + (t340 * t593 - t167 * (mrSges(4,1) * t322 - mrSges(4,3) * t451) - t168 * (-mrSges(4,2) * t322 - mrSges(4,3) * t460) + t618 * qJD(1)) * qJD(1) + t580 * t98 + t581 * t134 + (t1 * t138 + t126 * t4 + t139 * t3 + t581 * t29 + t582 * t30 + t580 * t43) * m(7) + t582 * t132 + t143 * t433 / 0.2e1 + (-t579 / 0.2e1 + t619) * t198 - t620 * t199 + t142 * t415 / 0.2e1 + t626 * t173 - t627 * t172 + (t408 + t361) * qJD(3) - (t307 + t567) * t439 / 0.2e1 + (t341 * t493 + t375 * t492) * t249 + (t39 - t38) * t176 + (-t180 * t435 - t181 * t433 - t321 * t117 + m(4) * ((-t167 * t325 - t168 * t321) * qJD(3) + t560) + t325 * t118) * pkin(8) + (-t167 * t433 + t560) * mrSges(4,3) - t361 * t439 - t435 * t477 - t370 * t430 / 0.2e1 + (-t10 * t176 + t354 * t304 + t177 * t9 + (t111 - t90) * t61 + t603 * t60 + t566 * t182) * m(5) + (t11 * t157 + t176 * t7 + t177 * t6 + t577 * t69 + t604 * t49 + (t112 - t78) * t48) * m(6) + t325 * t71 / 0.2e1 - t304 * t22 + Ifges(3,5) * t260 + Ifges(3,6) * t259 + t321 * t535 + t371 * t528 + t374 * t529 - t244 * mrSges(3,2) - t245 * mrSges(3,1) + Ifges(3,3) * qJDD(2) - t179 * t180 - t178 * t181 + t157 * t20 - t89 * t135 - t78 * t136 + t138 * t37 + t139 * t40 + t126 * t21 - t77 * t131 - t90 * t133 - t372 * t346 / 0.2e1 - pkin(2) * t91 + t269 * t308 + t368 * t513 - t386 * t380 + (pkin(2) * t386 - t167 * t178 - t168 * t179 - t271 * t309) * m(4); (m(5) * t60 + t568) * t66 + (m(7) * t29 + t134 - t546) * pkin(3) * t432 + t574 * t131 + (t299 * t6 + t303 * t7 - t48 * t66 + t49 * t574 - t69 * t76) * m(6) + t575 * t132 + (t1 * t298 - t29 * t35 + t299 * t3 + t30 * t575 - t43 * t47) * m(7) + t583 * t299 + t619 * t338 - t620 * t161 - t99 * t499 + (-Ifges(4,2) * t249 + t143 - t242) * t593 - t614 - m(5) * (t182 * t499 + t61 * t67) + t418 + t622 * t133 + (mrSges(4,1) * t225 - mrSges(4,2) * t226 - m(7) * (-t200 + t337) - m(6) * t337 + m(5) * t352 + t547) * g(2) + (-mrSges(4,1) * t227 + mrSges(4,2) * t228 - m(7) * (-t203 + t342) - m(6) * t342 - m(5) * t366 + t548) * g(1) + (-m(7) * (t275 - t293) + t239 - m(6) * (-t293 + t390) + t563) * g(3) - (-t321 * t540 + t377) * t494 + t298 * t37 + t303 * t39 + (-Ifges(4,5) * t359 - Ifges(4,6) * t249) * t631 + (-Ifges(4,1) * t359 - t483) * t632 + (t10 * t324 + t320 * t9 + (-t320 * t60 + t324 * t61) * qJD(4)) * t540 + t168 * t181 - t35 * t134 - t76 * t97 - t47 * t98 + (-t348 - t180) * t167 + t331 + t249 * t477 + t38 * t497 + t41 * t498 + t579 * t523 + t142 * t512 - t271 * (t249 * mrSges(4,1) - mrSges(4,2) * t359); (-pkin(4) * t7 - g(3) * t390 + qJ(5) * t6 + qJD(5) * t49 - t69 * t96) * m(6) - t182 * (mrSges(5,1) * t338 - mrSges(5,2) * t161) - t69 * (mrSges(6,1) * t338 + mrSges(6,3) * t161) - t43 * (-mrSges(7,1) * t338 - mrSges(7,2) * t161) + (-Ifges(7,5) * t161 + Ifges(7,6) * t338) * t509 + t563 * g(3) + (-m(6) * t49 - t489 - t569) * t60 + t570 * qJD(5) + (t338 * t587 - t599) * t526 + (-t161 * t556 - t480 + t579 + t601) * t523 + t583 * qJ(5) + (t488 + t546) * t61 - t377 * t494 + (-m(7) * (-t200 + t396) - m(6) * t396 + t547) * g(2) + (-m(7) * (-t203 + t395) - m(6) * t395 + t548) * g(1) + (-Ifges(5,2) * t338 - t156 + t552) * t525 - t30 * t486 - t29 * t487 - t33 * t134 - t32 * t132 - t96 * t97 - t65 * t98 + (-t161 * t590 + t338 * t586) * t508 + t331 + t49 * t490 + t48 * t491 + t85 * t522 + (-g(3) * t275 + t3 * qJ(5) - t1 * t532 - t29 * t33 + t30 * t571 - t43 * t65) * m(7) - t532 * t37 - pkin(4) * t39; t570 * t282 + (t97 - t98) * t338 + t37 + t39 + (t282 * t30 - t338 * t43 + t1 + t347) * m(7) + (t282 * t49 + t338 * t69 + t347 + t7) * m(6); -t161 * t132 + t338 * t134 + (-g(3) * t326 - t161 * t30 + t29 * t338 + t322 * t598 + t4) * m(7) + t21;];
tau  = t2;
