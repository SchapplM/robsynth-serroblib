% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:45
% EndTime: 2019-03-09 04:22:33
% DurationCPUTime: 30.42s
% Computational Cost: add. (22386->892), mult. (71387->1164), div. (0->0), fcn. (61259->14), ass. (0->419)
t319 = sin(pkin(12));
t321 = sin(pkin(6));
t324 = sin(qJ(3));
t469 = cos(pkin(7));
t497 = cos(qJ(3));
t401 = t469 * t497;
t468 = cos(pkin(12));
t342 = t321 * (t319 * t401 + t324 * t468);
t246 = qJD(1) * t342;
t320 = sin(pkin(7));
t451 = qJD(3) * t324;
t428 = t320 * t451;
t599 = t428 - t246;
t587 = mrSges(5,2) - mrSges(4,1);
t586 = -mrSges(4,3) - mrSges(5,1);
t327 = cos(qJ(5));
t520 = pkin(3) + pkin(10);
t448 = qJD(5) * t520;
t395 = t469 * t468;
t367 = t324 * t395;
t470 = cos(pkin(6));
t413 = t470 * t320;
t402 = t324 * t413;
t227 = t402 + (t319 * t497 + t367) * t321;
t218 = t227 * qJD(1);
t368 = t321 * t395;
t352 = t497 * t368;
t345 = qJD(1) * t352;
t370 = t497 * t413;
t452 = qJD(1) * t321;
t430 = t319 * t452;
t215 = -qJD(1) * t370 + t324 * t430 - t345;
t467 = qJ(4) * t215;
t117 = t218 * t520 + t467;
t323 = sin(qJ(5));
t396 = t470 * t468;
t316 = pkin(1) * t396;
t309 = qJD(1) * t316;
t431 = t470 * pkin(2);
t461 = t319 * t321;
t343 = t431 + (-pkin(9) * t469 - qJ(2)) * t461;
t211 = qJD(1) * t343 + t309;
t415 = t324 * t469;
t197 = t211 * t415;
t418 = t319 * t470;
t315 = pkin(1) * t418;
t416 = t321 * t468;
t400 = qJD(1) * t416;
t264 = qJ(2) * t400 + qJD(1) * t315;
t341 = (t368 + t413) * pkin(9);
t202 = qJD(1) * t341 + t264;
t252 = (-pkin(9) * t319 * t320 - pkin(2) * t468 - pkin(1)) * t321;
t239 = qJD(1) * t252 + qJD(2);
t460 = t320 * t324;
t114 = t497 * t202 + t239 * t460 + t197;
t95 = -pkin(4) * t215 + t114;
t58 = t327 * t117 + t323 * t95;
t604 = -t327 * t448 - t58;
t432 = t320 * t497;
t113 = t202 * t324 - t211 * t401 - t239 * t432;
t554 = -qJD(4) - t113;
t603 = Ifges(5,1) + Ifges(4,3);
t585 = Ifges(5,4) - Ifges(4,5);
t584 = Ifges(5,5) - Ifges(4,6);
t399 = pkin(5) * t327 + pkin(11) * t323;
t602 = qJD(5) * t399 - (-pkin(4) - t399) * t218 - t554;
t601 = -pkin(11) * t215 - t604;
t349 = -t323 * t469 - t327 * t432;
t408 = t320 * t430;
t564 = qJD(5) * t349 + t323 * t599 - t327 * t408;
t371 = t497 * t416;
t417 = t321 * t469;
t405 = t319 * t417;
t247 = (-t324 * t405 + t371) * qJD(1);
t424 = qJD(3) * t497;
t406 = t320 * t424;
t600 = t406 - t247;
t210 = qJD(5) + t218;
t397 = t470 * t469;
t261 = -qJD(1) * t397 + t320 * t400 - qJD(3);
t358 = pkin(4) * t218 + t113;
t592 = qJD(4) + t358;
t70 = t261 * t520 + t592;
t161 = -t211 * t320 + t469 * t239;
t366 = -qJ(4) * t218 + t161;
t72 = t215 * t520 + t366;
t32 = -t323 * t72 + t327 * t70;
t28 = -pkin(5) * t210 - t32;
t167 = t215 * t323 - t261 * t327;
t486 = mrSges(6,3) * t167;
t122 = mrSges(6,1) * t210 - t486;
t322 = sin(qJ(6));
t326 = cos(qJ(6));
t119 = -t167 * t322 + t210 * t326;
t120 = t167 * t326 + t210 * t322;
t68 = -mrSges(7,1) * t119 + mrSges(7,2) * t120;
t576 = -t122 + t68;
t598 = -m(7) * t28 - t576;
t274 = qJ(2) * t416 + t315;
t222 = t341 + t274;
t230 = t316 + t343;
t369 = qJD(3) * t401;
t106 = -t324 * (qJD(2) * t405 + qJD(3) * t222) + qJD(2) * t371 + t230 * t369 + t252 * t406;
t269 = -t320 * t416 + t397;
t102 = -t269 * qJD(4) - t106;
t596 = -t586 * t218 - t261 * t587;
t325 = sin(qJ(1));
t328 = cos(qJ(1));
t346 = t328 * t319 + t325 * t396;
t234 = t346 * t320 + t325 * t417;
t439 = t324 * t461;
t595 = (-t370 + t439) * qJD(3);
t270 = t325 * t319 - t328 * t396;
t594 = -t270 * t320 + t328 * t417;
t388 = mrSges(7,1) * t322 + mrSges(7,2) * t326;
t360 = -m(7) * pkin(10) - t388;
t590 = -m(7) - m(6);
t593 = m(6) * pkin(10) - t360 + pkin(3) * (m(5) - t590) - t587;
t271 = t325 * t468 + t328 * t418;
t410 = t321 * t432;
t181 = t270 * t401 + t271 * t324 + t328 * t410;
t591 = t181 * t327 + t323 * t594;
t140 = t181 * t323 - t327 * t594;
t444 = qJDD(1) * t321;
t423 = t319 * t444;
t155 = qJD(1) * t595 - qJD(3) * t345 - qJDD(1) * t402 - t367 * t444 - t497 * t423;
t151 = qJDD(5) - t155;
t217 = t227 * qJD(3);
t562 = -t352 - t370;
t156 = qJD(1) * t217 + qJDD(1) * t562 + t324 * t423;
t166 = t215 * t327 + t261 * t323;
t260 = qJDD(1) * t269 + qJDD(3);
t88 = qJD(5) * t166 + t156 * t323 + t260 * t327;
t42 = qJD(6) * t119 + t151 * t322 + t326 * t88;
t530 = t42 / 0.2e1;
t43 = -qJD(6) * t120 + t151 * t326 - t322 * t88;
t529 = t43 / 0.2e1;
t89 = -qJD(5) * t167 + t156 * t327 - t260 * t323;
t83 = qJDD(6) - t89;
t524 = t83 / 0.2e1;
t254 = t261 * qJ(4);
t73 = -t254 + t95;
t589 = t73 * mrSges(6,2);
t588 = mrSges(4,2) - mrSges(5,3);
t14 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t63 = mrSges(6,1) * t151 - mrSges(6,3) * t88;
t583 = t14 - t63;
t164 = Ifges(6,4) * t166;
t582 = t166 * Ifges(6,2);
t581 = t166 * Ifges(6,6);
t580 = t210 * Ifges(6,5);
t579 = t210 * Ifges(6,6);
t578 = t210 * Ifges(6,3);
t389 = -mrSges(7,1) * t326 + mrSges(7,2) * t322;
t361 = m(7) * pkin(5) - t389;
t577 = mrSges(6,1) + t361;
t440 = m(7) * pkin(11) + mrSges(7,3);
t547 = mrSges(6,2) - t440;
t127 = mrSges(5,1) * t156 - mrSges(5,3) * t260;
t46 = -mrSges(6,1) * t89 + mrSges(6,2) * t88;
t575 = -t127 + t46;
t209 = Ifges(4,4) * t215;
t574 = t218 * Ifges(4,1) - t261 * Ifges(4,5) + t167 * Ifges(6,5) - t209 + t578 + t581;
t419 = -pkin(11) * t327 + qJ(4);
t298 = pkin(5) * t323 + t419;
t456 = t323 * t520;
t255 = t326 * t298 + t322 * t456;
t573 = qJD(6) * t255 + t322 * t602 - t326 * t601;
t256 = t322 * t298 - t326 * t456;
t572 = -qJD(6) * t256 + t322 * t601 + t326 * t602;
t108 = -mrSges(6,1) * t166 + mrSges(6,2) * t167;
t170 = mrSges(5,1) * t215 + mrSges(5,3) * t261;
t571 = t108 - t170;
t464 = t218 * t323;
t144 = -t215 * t326 - t322 * t464;
t450 = qJD(5) * t323;
t427 = t322 * t450;
t445 = qJD(6) * t327;
t570 = t326 * t445 + t144 - t427;
t145 = -t215 * t322 + t326 * t464;
t569 = t322 * t445 + t326 * t450 + t145;
t278 = -t323 * t432 + t327 * t469;
t241 = t278 * t326 + t322 * t460;
t568 = -qJD(6) * t241 - t322 * t564 + t326 * t600;
t240 = -t278 * t322 + t326 * t460;
t567 = qJD(6) * t240 + t322 * t600 + t326 * t564;
t169 = mrSges(4,2) * t261 - mrSges(4,3) * t215;
t566 = t169 - t170;
t391 = mrSges(6,1) * t323 + mrSges(6,2) * t327;
t561 = -m(6) * qJ(4) - m(7) * t419 + t327 * mrSges(7,3) - t323 * t361 - t391;
t559 = t585 * t155 + t584 * t156 + t260 * t603;
t449 = qJD(5) * t327;
t463 = t218 * t327;
t558 = t449 + t463;
t24 = mrSges(7,1) * t83 - mrSges(7,3) * t42;
t25 = -mrSges(7,2) * t83 + mrSges(7,3) * t43;
t557 = -t322 * t24 + t326 * t25;
t282 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t321;
t244 = qJDD(1) * t315 + t468 * t282;
t192 = qJDD(1) * t341 + t244;
t243 = qJDD(1) * t316 - t319 * t282;
t193 = (-pkin(9) * t405 + t431) * qJDD(1) + t243;
t235 = qJDD(1) * t252 + qJDD(2);
t62 = -qJD(3) * t197 - t324 * t192 + t193 * t401 - t202 * t424 + t235 * t432 - t239 * t428;
t333 = qJDD(4) - t62;
t35 = -t155 * pkin(4) - t260 * t520 + t333;
t143 = -t193 * t320 + t469 * t235;
t347 = qJ(4) * t155 - qJD(4) * t218 + t143;
t45 = t156 * t520 + t347;
t7 = t323 * t35 + t327 * t45 + t70 * t449 - t450 * t72;
t33 = t323 * t70 + t327 * t72;
t8 = -qJD(5) * t33 - t323 * t45 + t327 * t35;
t556 = -t7 * t323 - t8 * t327;
t29 = pkin(11) * t210 + t33;
t51 = -pkin(5) * t166 - pkin(11) * t167 + t73;
t12 = -t29 * t322 + t326 * t51;
t61 = t497 * t192 + t193 * t415 - t202 * t451 + t211 * t369 + t235 * t460 + t239 * t406;
t52 = -t260 * qJ(4) + t261 * qJD(4) - t61;
t36 = -pkin(4) * t156 - t52;
t15 = -pkin(5) * t89 - pkin(11) * t88 + t36;
t5 = pkin(11) * t151 + t7;
t1 = qJD(6) * t12 + t15 * t322 + t326 * t5;
t13 = t29 * t326 + t322 * t51;
t2 = -qJD(6) * t13 + t15 * t326 - t322 * t5;
t555 = t1 * t326 - t2 * t322;
t503 = t227 / 0.2e1;
t553 = -t227 / 0.2e1 - t503;
t502 = t260 / 0.2e1;
t552 = -t502 - t260 / 0.2e1;
t551 = -mrSges(6,3) + t587;
t458 = t321 * t328;
t459 = t321 * t325;
t549 = -g(1) * t459 + g(2) * t458 - g(3) * t470;
t338 = t324 * t222 - t230 * t401 - t252 * t432;
t78 = t227 * pkin(4) - t269 * t520 + t338;
t168 = -t230 * t320 + t469 * t252;
t226 = t439 + t562;
t224 = t226 * pkin(3);
t411 = qJ(4) * t227 - t224;
t110 = t168 - t411;
t496 = pkin(10) * t226;
t92 = t110 + t496;
t488 = t323 * t78 + t327 * t92;
t207 = t497 * t222;
t412 = t469 * t230;
t107 = qJD(2) * t342 + (t207 + (t252 * t320 + t412) * t324) * qJD(3);
t216 = -qJD(3) * t352 + t595;
t87 = -t216 * pkin(4) + t107;
t429 = qJD(2) * t461;
t291 = t320 * t429;
t364 = qJ(4) * t216 - qJD(4) * t227 + t291;
t96 = t217 * t520 + t364;
t19 = -qJD(5) * t488 - t323 * t96 + t327 * t87;
t546 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t545 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t498 = t327 / 0.2e1;
t476 = t167 * Ifges(6,4);
t85 = t476 + t579 + t582;
t523 = -t85 / 0.2e1;
t165 = qJD(6) - t166;
t54 = t120 * Ifges(7,5) + t119 * Ifges(7,6) + t165 * Ifges(7,3);
t544 = t73 * (mrSges(6,1) * t327 - mrSges(6,2) * t323) + t327 * t523 + t54 * t498;
t515 = t151 / 0.2e1;
t522 = t88 / 0.2e1;
t542 = Ifges(6,1) * t522 + Ifges(6,5) * t515;
t541 = -m(5) * qJ(4) + t561 + t588;
t540 = t73 * mrSges(6,1) + t12 * mrSges(7,1) - t13 * mrSges(7,2);
t272 = -t325 * t418 + t328 * t468;
t337 = t346 * t469;
t185 = t272 * t324 - t325 * t410 + t337 * t497;
t539 = g(1) * t185 + g(2) * t181 + g(3) * t226;
t103 = pkin(3) * t215 + t366;
t105 = t254 - t114;
t538 = mrSges(4,1) * t161 + t105 * mrSges(5,1) - mrSges(5,2) * t103 - t114 * mrSges(4,3);
t53 = -t260 * pkin(3) + t333;
t537 = t62 * mrSges(4,1) - t61 * mrSges(4,2) + t53 * mrSges(5,2) - t52 * mrSges(5,3);
t104 = pkin(3) * t261 - t554;
t536 = t104 * mrSges(5,1) + mrSges(6,1) * t32 + mrSges(4,2) * t161 - mrSges(6,2) * t33 + t113 * mrSges(4,3) - mrSges(5,3) * t103;
t521 = t89 / 0.2e1;
t9 = Ifges(7,5) * t42 + Ifges(7,6) * t43 + Ifges(7,3) * t83;
t535 = Ifges(7,5) * t530 + Ifges(7,6) * t529 + Ifges(7,3) * t524 - t88 * Ifges(6,4) / 0.2e1 + t9 / 0.2e1 + t545 + (-t515 - t151 / 0.2e1) * Ifges(6,6) + (-t521 - t89 / 0.2e1) * Ifges(6,2);
t10 = t42 * Ifges(7,4) + t43 * Ifges(7,2) + t83 * Ifges(7,6);
t533 = t10 / 0.2e1;
t532 = Ifges(7,1) * t530 + Ifges(7,4) * t529 + Ifges(7,5) * t524;
t483 = Ifges(7,4) * t120;
t55 = Ifges(7,2) * t119 + Ifges(7,6) * t165 + t483;
t528 = -t55 / 0.2e1;
t527 = t55 / 0.2e1;
t118 = Ifges(7,4) * t119;
t56 = Ifges(7,1) * t120 + Ifges(7,5) * t165 + t118;
t526 = -t56 / 0.2e1;
t525 = t56 / 0.2e1;
t519 = -t119 / 0.2e1;
t518 = t119 / 0.2e1;
t517 = -t120 / 0.2e1;
t516 = t120 / 0.2e1;
t514 = -t165 / 0.2e1;
t513 = t165 / 0.2e1;
t512 = -t166 / 0.2e1;
t511 = -t167 / 0.2e1;
t510 = t167 / 0.2e1;
t509 = -t210 / 0.2e1;
t508 = -t215 / 0.2e1;
t507 = t215 / 0.2e1;
t505 = -t218 / 0.2e1;
t504 = t218 / 0.2e1;
t501 = -t261 / 0.2e1;
t500 = t261 / 0.2e1;
t6 = -pkin(5) * t151 - t8;
t493 = t327 * t6;
t487 = mrSges(6,3) * t166;
t485 = Ifges(6,4) * t323;
t484 = Ifges(6,4) * t327;
t482 = Ifges(7,4) * t322;
t481 = Ifges(7,4) * t326;
t475 = t218 * Ifges(4,4);
t474 = t218 * Ifges(5,6);
t466 = t166 * t322;
t465 = t166 * t326;
t457 = t322 * t327;
t455 = t326 * t327;
t453 = t328 * pkin(1) + qJ(2) * t459;
t447 = qJD(6) * t322;
t446 = qJD(6) * t326;
t441 = Ifges(6,5) * t88 + Ifges(6,6) * t89 + Ifges(6,3) * t151;
t124 = t252 * t460 + t324 * t412 + t207;
t426 = t323 * t448;
t422 = -t450 / 0.2e1;
t421 = -t445 / 0.2e1;
t420 = -pkin(1) * t325 + qJ(2) * t458;
t128 = -t155 * mrSges(5,1) + t260 * mrSges(5,2);
t111 = -t269 * qJ(4) - t124;
t398 = qJDD(1) * t416;
t394 = mrSges(4,1) * t226 + mrSges(4,2) * t227;
t180 = t226 * t323 + t269 * t327;
t372 = t226 * t327 - t269 * t323;
t393 = -mrSges(6,1) * t372 + mrSges(6,2) * t180;
t135 = -t180 * t322 + t227 * t326;
t136 = t180 * t326 + t227 * t322;
t390 = mrSges(7,1) * t135 - mrSges(7,2) * t136;
t387 = -t226 * mrSges(5,2) - t227 * mrSges(5,3);
t386 = Ifges(6,1) * t323 + t484;
t385 = Ifges(7,1) * t326 - t482;
t384 = Ifges(7,1) * t322 + t481;
t383 = Ifges(6,2) * t327 + t485;
t382 = -Ifges(7,2) * t322 + t481;
t381 = Ifges(7,2) * t326 + t482;
t380 = Ifges(6,5) * t323 + Ifges(6,6) * t327;
t379 = Ifges(7,5) * t326 - Ifges(7,6) * t322;
t378 = Ifges(7,5) * t322 + Ifges(7,6) * t326;
t377 = t32 * t323 - t327 * t33;
t39 = pkin(11) * t227 + t488;
t97 = -pkin(4) * t226 - t111;
t60 = -pkin(5) * t372 - pkin(11) * t180 + t97;
t21 = t322 * t60 + t326 * t39;
t20 = -t322 * t39 + t326 * t60;
t47 = -t323 * t92 + t327 * t78;
t57 = -t117 * t323 + t327 * t95;
t18 = t323 * t87 + t327 * t96 + t78 * t449 - t450 * t92;
t365 = t28 * t388;
t363 = -(-qJ(2) * t430 + t309) * t319 + t264 * t468;
t357 = -mrSges(3,1) * t398 + mrSges(3,2) * t423;
t355 = mrSges(3,1) * t470 - mrSges(3,3) * t461;
t351 = -t271 * pkin(2) + pkin(9) * t594 + t420;
t350 = -mrSges(3,2) * t470 + mrSges(3,3) * t416;
t182 = -t270 * t415 + t271 * t497 - t458 * t460;
t340 = (-t12 * t326 - t13 * t322) * qJD(6) + t555;
t339 = -qJD(5) * t377 - t556;
t336 = -pkin(3) * t182 - qJ(4) * t181 + t351;
t334 = t272 * pkin(2) + pkin(9) * t234 + t453;
t71 = -t217 * pkin(4) - t102;
t186 = t272 * t497 + (t320 * t459 - t337) * t324;
t331 = t186 * pkin(3) + t185 * qJ(4) + t334;
t330 = t234 * pkin(4) + t186 * pkin(10) + t331;
t310 = -pkin(1) * t444 + qJDD(2);
t276 = t350 * qJD(1);
t275 = t355 * qJD(1);
t273 = -qJ(2) * t461 + t316;
t208 = Ifges(5,6) * t215;
t163 = -t218 * t455 - t261 * t322;
t162 = t218 * t457 - t261 * t326;
t154 = -mrSges(5,2) * t215 - mrSges(5,3) * t218;
t153 = mrSges(4,1) * t215 + mrSges(4,2) * t218;
t152 = pkin(3) * t218 + t467;
t142 = t185 * t323 + t234 * t327;
t141 = -t185 * t327 + t234 * t323;
t134 = qJD(5) * t180 - t217 * t327;
t133 = qJD(5) * t372 + t217 * t323;
t131 = -t215 * Ifges(4,2) - t261 * Ifges(4,6) + t475;
t130 = -t261 * Ifges(5,4) - t218 * Ifges(5,2) + t208;
t129 = -t261 * Ifges(5,5) + t215 * Ifges(5,3) - t474;
t126 = -mrSges(4,2) * t260 - mrSges(4,3) * t156;
t125 = mrSges(4,1) * t260 + mrSges(4,3) * t155;
t121 = -mrSges(6,2) * t210 + t487;
t115 = pkin(3) * t217 + t364;
t112 = -t269 * pkin(3) + t338;
t109 = pkin(5) * t167 - pkin(11) * t166;
t101 = t142 * t326 + t186 * t322;
t100 = -t142 * t322 + t186 * t326;
t99 = mrSges(4,1) * t156 - mrSges(4,2) * t155;
t98 = -mrSges(5,2) * t156 + mrSges(5,3) * t155;
t86 = t167 * Ifges(6,1) + t164 + t580;
t77 = mrSges(7,1) * t165 - mrSges(7,3) * t120;
t76 = -mrSges(7,2) * t165 + mrSges(7,3) * t119;
t67 = qJD(6) * t135 + t133 * t326 - t216 * t322;
t66 = -qJD(6) * t136 - t133 * t322 - t216 * t326;
t64 = -mrSges(6,2) * t151 + mrSges(6,3) * t89;
t59 = pkin(3) * t156 + t347;
t49 = pkin(5) * t215 - t57;
t38 = -pkin(5) * t227 - t47;
t37 = t134 * pkin(5) - t133 * pkin(11) + t71;
t31 = t88 * Ifges(6,1) + t89 * Ifges(6,4) + t151 * Ifges(6,5);
t27 = t109 * t322 + t32 * t326;
t26 = t109 * t326 - t32 * t322;
t17 = pkin(5) * t216 - t19;
t16 = -pkin(11) * t216 + t18;
t4 = -qJD(6) * t21 - t16 * t322 + t326 * t37;
t3 = qJD(6) * t20 + t16 * t326 + t322 * t37;
t11 = [qJD(2) * t276 * t416 + (Ifges(6,4) * t521 + t31 / 0.2e1 + t542) * t180 + m(7) * (t1 * t21 + t12 * t4 + t13 * t3 + t17 * t28 + t2 * t20 + t38 * t6) + (Ifges(7,1) * t136 + Ifges(7,4) * t135) * t530 + (Ifges(7,1) * t67 + Ifges(7,4) * t66) * t516 + (t553 * Ifges(4,1) - Ifges(5,2) * t227) * t155 + (t553 * Ifges(4,4) - Ifges(5,6) * t227) * t156 + (Ifges(7,4) * t136 + Ifges(7,2) * t135) * t529 + (Ifges(7,4) * t67 + Ifges(7,2) * t66) * t518 + m(3) * (t243 * t273 + t244 * t274) + (-t582 / 0.2e1 - Ifges(6,4) * t510 - t579 / 0.2e1 + t523 + Ifges(7,3) * t513 + Ifges(7,5) * t516 + Ifges(7,6) * t518 + t54 / 0.2e1 + t540) * t134 + (t470 * (Ifges(3,3) * t470 + (Ifges(3,5) * t319 + Ifges(3,6) * t468) * t321) + t274 * t350 + t273 * t355 + Ifges(2,3)) * qJDD(1) + (-t574 / 0.2e1 - Ifges(6,5) * t510 - Ifges(4,1) * t504 + Ifges(5,2) * t505 + Ifges(5,6) * t507 - Ifges(4,4) * t508 - t578 / 0.2e1 - t581 / 0.2e1 + t130 / 0.2e1 + t585 * t501 - t536) * t216 + m(6) * (t18 * t33 + t19 * t32 + t36 * t97 + t47 * t8 + t488 * t7 + t71 * t73) + t488 * t64 + (Ifges(7,5) * t136 + Ifges(7,6) * t135) * t524 + (Ifges(7,5) * t67 + Ifges(7,6) * t66) * t513 + (-m(3) * t420 - m(4) * t351 - m(5) * t336 + mrSges(2,1) * t325 + t271 * mrSges(3,1) + mrSges(2,2) * t328 - t270 * mrSges(3,2) - mrSges(3,3) * t458 + t590 * (pkin(4) * t594 - pkin(10) * t182 + t336) + t586 * t594 - t588 * t181 + t547 * t591 + t577 * t140 - (-t388 + t551) * t182) * g(1) + (m(4) * t113 + m(5) * t104 + t596) * t107 + (-m(3) * t453 - t272 * mrSges(3,1) + mrSges(3,2) * t346 - mrSges(3,3) * t459 - mrSges(2,1) * t328 + mrSges(2,2) * t325 - m(4) * t334 - m(5) * t331 - m(7) * (t142 * pkin(5) + t330) - t101 * mrSges(7,1) - t100 * mrSges(7,2) - m(6) * t330 - t142 * mrSges(6,1) + t586 * t234 + t588 * t185 + t547 * t141 + t551 * t186) * g(2) + (t164 / 0.2e1 + Ifges(6,1) * t510 + t580 / 0.2e1 + t86 / 0.2e1 + t589) * t133 + (-Ifges(4,4) * t504 + Ifges(5,6) * t505 + Ifges(5,3) * t507 - Ifges(4,2) * t508 + t129 / 0.2e1 - t131 / 0.2e1 + t584 * t501 + t538) * t217 + (t1 * t135 - t12 * t67 + t13 * t66 - t136 * t2) * mrSges(7,3) + (t53 * mrSges(5,1) - t62 * mrSges(4,3) + Ifges(5,4) * t552 + Ifges(4,5) * t502 + Ifges(6,5) * t522 + Ifges(6,6) * t521 + Ifges(6,3) * t515 + t546) * t227 + (Ifges(3,5) * t423 + Ifges(3,6) * t398) * t470 + (t559 / 0.2e1 + (Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1) * t155 + (Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1) * t156 + t603 * t502 + t537) * t269 + m(4) * (t106 * t114 + t124 * t61 + t143 * t168 + t161 * t291 - t338 * t62) - t338 * t125 + (-t133 * t32 - t134 * t33 - t180 * t8 + t372 * t7) * mrSges(6,3) - (-Ifges(6,4) * t522 + t535) * t372 - t6 * t390 + t36 * t393 + t143 * t394 + t59 * t387 + t38 * t14 + t243 * t355 + m(5) * (t102 * t105 + t103 * t115 + t110 * t59 + t111 * t52 + t112 * t53) + t153 * t291 + ((Ifges(3,1) * t319 + Ifges(3,4) * t468) * t423 + (Ifges(3,4) * t319 + Ifges(3,2) * t468) * t398 + t310 * (-mrSges(3,1) * t468 + mrSges(3,2) * t319) - pkin(1) * t357 + m(3) * (-pkin(1) * t310 + qJD(2) * t363)) * t321 - t275 * t429 + (t52 * mrSges(5,1) - t61 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t156 + (Ifges(5,6) + Ifges(4,4)) * t155 - t584 * t552) * t226 + t20 * t24 + t21 * t25 + t244 * t350 + t136 * t532 + t135 * t533 + t67 * t525 + t66 * t527 + t168 * t99 + t106 * t169 + t102 * t170 + t47 * t63 + t28 * (-mrSges(7,1) * t66 + mrSges(7,2) * t67) + t17 * t68 + t3 * t76 + t4 * t77 + t97 * t46 + t71 * t108 + t110 * t98 + t18 * t121 + t19 * t122 + t124 * t126 + t111 * t127 + t112 * t128 + (Ifges(4,5) * t260 + t441) * t503 + t115 * t154; (t126 + t575) * t460 + (t8 * t349 + t7 * t278 + (t324 * t36 + t424 * t73) * t320 - t247 * t73 + t564 * t33 + t549) * m(6) + t564 * t121 + (t98 + t99) * t469 + (-t128 + t125) * t432 + (-t113 * t246 - t114 * t247 - t161 * t408 + t143 * t469 + (t497 * t62 + t324 * t61 + (t113 * t324 + t114 * t497) * qJD(3)) * t320 + t549) * m(4) + (-t103 * t408 - t104 * t246 + t105 * t247 + t59 * t469 + (-t497 * t53 - t324 * t52 + (t104 * t324 - t105 * t497) * qJD(3)) * t320 + t549) * m(5) + (-t363 * t452 + t310 + t549) * m(3) - t583 * t349 + t567 * t76 + (t1 * t241 + t12 * t568 + t13 * t567 + t2 * t240 - t349 * t6 + t549) * m(7) + t568 * t77 + t275 * t430 - t276 * t400 + t357 + t278 * t64 + t240 * t24 + t241 * t25 + (-t154 - t153) * t408 + (-m(6) * t32 - t598) * (qJD(5) * t278 - t323 * t408 - t327 * t599) + t596 * t599 + t600 * (t108 + t566); t572 * t77 + t573 * t76 + (-t209 + t574) * t507 + t575 * qJ(4) + (t379 * t524 + t382 * t529 + t385 * t530 + t542) * t327 - t10 * t457 / 0.2e1 + (-m(6) * (-t224 - t496) - m(5) * t411 + t387 + t394 + m(7) * t224 - t226 * t360 + t561 * t227) * g(3) + t566 * t113 + (-t475 + t129) * t505 + t583 * t327 * t520 - (t166 * t383 + t167 * t386 + t210 * t380) * qJD(5) / 0.2e1 + t559 + (-t464 / 0.2e1 + t422) * t86 - t596 * t114 + t604 * t121 - t485 * t522 + (-Ifges(4,2) * t507 + Ifges(5,3) * t508 + t380 * t509 + t383 * t512 + t386 * t511 + t500 * t584 - t538 + t544) * t218 + (-Ifges(4,1) * t505 - Ifges(6,5) * t511 + Ifges(5,2) * t504 - Ifges(6,6) * t512 - Ifges(6,3) * t509 + t500 * t585 + t536) * t215 + (-t558 * t33 + (t450 + t464) * t32 + t539 + t556) * mrSges(6,3) + t484 * t521 + (mrSges(7,1) * t558 + mrSges(7,3) * t569) * t12 + (t570 * mrSges(7,1) - t569 * mrSges(7,2)) * t28 + (-mrSges(7,2) * t558 - mrSges(7,3) * t570) * t13 + t571 * qJD(4) + t537 + (t421 * t55 + t422 * t56) * t326 + (-pkin(3) * t53 - qJ(4) * t52 - t103 * t152 - t104 * t114 + t105 * t554) * m(5) + t535 * t323 + ((Ifges(7,3) * t327 - t323 * t379) * t513 + (Ifges(7,5) * t327 - t323 * t385) * t516 + (Ifges(7,6) * t327 - t323 * t382) * t518 + t544) * qJD(5) + t358 * t108 + t322 * t56 * t421 + t36 * t391 + (-t28 * t49 + t1 * t256 + t2 * t255 - (t28 * t450 - t493) * t520 + t573 * t13 + t572 * t12) * m(7) + (-t426 - t49) * t68 + (t426 - t57) * t122 + (t208 + t130) * t508 + (t474 + t131) * t504 + (-t378 * t513 - t381 * t518 - t384 * t516) * t445 + (t181 * t593 + t182 * t541) * g(2) + (t185 * t593 + t186 * t541) * g(1) + (t36 * qJ(4) - t32 * t57 - t33 * t58 - t520 * t339 + t592 * t73) * m(6) + (-t1 * t457 - t2 * t455) * mrSges(7,3) + t255 * t24 + t256 * t25 + t455 * t532 + t145 * t526 + t427 * t527 + t144 * t528 + (Ifges(7,5) * t145 + Ifges(7,6) * t144 - Ifges(7,3) * t463) * t514 + (Ifges(7,1) * t145 + Ifges(7,4) * t144 - Ifges(7,5) * t463) * t517 + (Ifges(7,4) * t145 + Ifges(7,2) * t144 - Ifges(7,6) * t463) * t519 + t31 * t498 + t388 * t493 - t64 * t456 - pkin(3) * t128 - t152 * t154; t218 * t154 - t162 * t77 - t163 * t76 + t571 * t261 + (t218 * t121 + (-t322 * t77 + t326 * t76 + t121) * qJD(5) - t583) * t327 + (t64 + (-t322 * t76 - t326 * t77) * qJD(6) + t210 * t576 + t557) * t323 + t128 + (-t539 + (-t6 + (-t12 * t322 + t13 * t326) * qJD(5)) * t327 + (qJD(5) * t28 + t340) * t323 - t12 * t162 - t13 * t163 + t28 * t464) * m(7) + (-t218 * t377 + t261 * t73 + t339 - t539) * m(6) + (t103 * t218 - t105 * t261 + t53 - t539) * m(5); (t141 * t577 + t142 * t547) * g(1) + (-t476 + t54) * t511 + (t486 + t598) * t33 + (t86 + t164) * t512 + t546 + (t119 * t382 + t120 * t385 + t165 * t379) * qJD(6) / 0.2e1 + (Ifges(7,5) * t517 - Ifges(6,2) * t512 - Ifges(6,6) * t509 + Ifges(7,6) * t519 + Ifges(7,3) * t514 - t540) * t167 + (-pkin(5) * t6 - t12 * t26 - t13 * t27) * m(7) + (Ifges(6,1) * t511 + Ifges(6,5) * t509 + t379 * t514 + t382 * t519 + t385 * t517 - t365 - t589) * t166 + ((-t447 + t466) * t13 + (-t446 + t465) * t12 + t555) * mrSges(7,3) + (m(7) * t340 - t446 * t77 - t447 * t76 + t557) * pkin(11) + (-t180 * t440 - t361 * t372 + t393) * g(3) + t441 + t6 * t389 - pkin(5) * t14 + (t547 * t140 - t577 * t591) * g(2) + qJD(6) * t365 + t322 * t532 + t326 * t533 + t378 * t524 + t446 * t525 + t465 * t526 + t466 * t527 + t447 * t528 + t381 * t529 + t384 * t530 + t85 * t510 - t27 * t76 - t26 * t77 + (t487 - t121) * t32; -t28 * (mrSges(7,1) * t120 + mrSges(7,2) * t119) + (Ifges(7,1) * t119 - t483) * t517 + t55 * t516 + (Ifges(7,5) * t119 - Ifges(7,6) * t120) * t514 - t12 * t76 + t13 * t77 - g(1) * (mrSges(7,1) * t100 - mrSges(7,2) * t101) - g(2) * ((-t140 * t322 + t182 * t326) * mrSges(7,1) + (-t140 * t326 - t182 * t322) * mrSges(7,2)) - g(3) * t390 + (t119 * t12 + t120 * t13) * mrSges(7,3) + t9 + (-Ifges(7,2) * t120 + t118 + t56) * t519 + t545;];
tau  = t11;
