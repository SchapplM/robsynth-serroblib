% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:26
% EndTime: 2019-03-09 05:16:15
% DurationCPUTime: 31.80s
% Computational Cost: add. (19537->844), mult. (47116->1104), div. (0->0), fcn. (36848->18), ass. (0->364)
t316 = sin(pkin(10));
t318 = cos(pkin(10));
t323 = sin(qJ(3));
t327 = cos(qJ(3));
t523 = -t316 * t323 + t327 * t318;
t266 = t523 * qJD(1);
t275 = t316 * t327 + t318 * t323;
t268 = t275 * qJD(1);
t200 = pkin(3) * t268 - pkin(8) * t266;
t427 = pkin(7) + qJ(2);
t286 = t427 * t316;
t279 = qJD(1) * t286;
t287 = t427 * t318;
t280 = qJD(1) * t287;
t207 = -t279 * t327 - t323 * t280;
t322 = sin(qJ(4));
t326 = cos(qJ(4));
t143 = t326 * t200 - t207 * t322;
t319 = -qJ(5) - pkin(8);
t365 = qJD(4) * t319;
t411 = t266 * t326;
t541 = -pkin(4) * t268 + qJ(5) * t411 - qJD(5) * t322 + t326 * t365 - t143;
t144 = t322 * t200 + t326 * t207;
t412 = t266 * t322;
t540 = -qJ(5) * t412 - qJD(5) * t326 - t322 * t365 + t144;
t315 = sin(pkin(11));
t317 = cos(pkin(11));
t274 = t315 * t326 + t317 * t322;
t173 = t274 * t266;
t265 = t274 * qJD(4);
t493 = t265 - t173;
t342 = t315 * t322 - t317 * t326;
t174 = t342 * t266;
t267 = t342 * qJD(4);
t539 = t267 - t174;
t526 = Ifges(6,3) + Ifges(5,3);
t508 = t540 * t315 + t317 * t541;
t507 = t315 * t541 - t540 * t317;
t383 = qJD(1) * qJD(2);
t294 = qJ(2) * qJDD(1) + t383;
t538 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t537 = -pkin(5) * t268 + pkin(9) * t539 + t508;
t536 = -pkin(9) * t493 + t507;
t389 = t316 ^ 2 + t318 ^ 2;
t464 = m(6) + m(7);
t520 = -t464 - m(4) - m(5);
t302 = pkin(2) * t318 + pkin(1);
t285 = -qJD(1) * t302 + qJD(2);
t503 = Ifges(4,5) * qJD(3);
t535 = -t285 * mrSges(4,2) + t207 * mrSges(4,3) - t503 / 0.2e1;
t270 = t275 * qJD(3);
t170 = -pkin(3) * t266 - pkin(8) * t268 + t285;
t208 = -t323 * t279 + t327 * t280;
t199 = qJD(3) * pkin(8) + t208;
t123 = t326 * t170 - t199 * t322;
t124 = t170 * t322 + t199 * t326;
t321 = sin(qJ(6));
t325 = cos(qJ(6));
t252 = qJD(4) - t266;
t220 = qJD(3) * t326 - t268 * t322;
t221 = qJD(3) * t322 + t268 * t326;
t162 = t220 * t315 + t221 * t317;
t527 = pkin(9) * t162;
t108 = qJ(5) * t220 + t124;
t100 = t315 * t108;
t107 = -qJ(5) * t221 + t123;
t93 = pkin(4) * t252 + t107;
t54 = t317 * t93 - t100;
t39 = pkin(5) * t252 - t527 + t54;
t360 = t317 * t220 - t221 * t315;
t511 = pkin(9) * t360;
t400 = t317 * t108;
t55 = t315 * t93 + t400;
t40 = t55 + t511;
t13 = -t321 * t40 + t325 * t39;
t14 = t321 * t39 + t325 * t40;
t502 = Ifges(4,6) * qJD(3);
t534 = -t285 * mrSges(4,1) - t123 * mrSges(5,1) - t54 * mrSges(6,1) - t13 * mrSges(7,1) + t124 * mrSges(5,2) + t55 * mrSges(6,2) + t14 * mrSges(7,2) + t208 * mrSges(4,3) + t502 / 0.2e1;
t314 = qJ(4) + pkin(11);
t308 = cos(t314);
t434 = pkin(4) * t326;
t284 = pkin(5) * t308 + t434;
t281 = pkin(3) + t284;
t309 = qJ(6) + t314;
t299 = sin(t309);
t300 = cos(t309);
t303 = pkin(3) + t434;
t306 = sin(t314);
t352 = -mrSges(5,1) * t326 + mrSges(5,2) * t322;
t533 = -m(5) * pkin(3) - m(6) * t303 - m(7) * t281 - mrSges(6,1) * t308 - mrSges(7,1) * t300 + mrSges(6,2) * t306 + mrSges(7,2) * t299 + t352;
t310 = -pkin(9) + t319;
t532 = -m(5) * pkin(8) + m(6) * t319 + m(7) * t310 - t538;
t484 = m(6) * pkin(4);
t435 = pkin(4) * t322;
t529 = m(6) * t435;
t283 = pkin(5) * t306 + t435;
t528 = m(7) * t283;
t525 = t220 * Ifges(5,6);
t386 = qJD(4) * t322;
t496 = -t208 + (t386 - t412) * pkin(4);
t494 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t220 - mrSges(5,2) * t221 - mrSges(4,3) * t268;
t269 = t523 * qJD(3);
t385 = qJD(4) * t326;
t371 = t275 * t385;
t339 = t269 * t322 + t371;
t203 = qJD(1) * t269 + qJDD(1) * t275;
t154 = qJD(4) * t220 + qJDD(3) * t322 + t203 * t326;
t155 = -qJD(4) * t221 + qJDD(3) * t326 - t203 * t322;
t380 = qJDD(1) * t318;
t381 = qJDD(1) * t316;
t204 = -qJD(1) * t270 - t323 * t381 + t327 * t380;
t197 = qJDD(4) - t204;
t83 = -t154 * t315 + t155 * t317;
t84 = t154 * t317 + t155 * t315;
t524 = Ifges(5,5) * t154 + Ifges(6,5) * t84 + Ifges(5,6) * t155 + Ifges(6,6) * t83 + t197 * t526;
t324 = sin(qJ(1));
t328 = cos(qJ(1));
t522 = g(1) * t328 + g(2) * t324;
t244 = qJD(6) + t252;
t517 = -t162 * t321 + t325 * t360;
t99 = t162 * t325 + t321 * t360;
t521 = t221 * Ifges(5,5) + t162 * Ifges(6,5) + t99 * Ifges(7,5) + Ifges(6,6) * t360 + Ifges(7,6) * t517 + t244 * Ifges(7,3) + t252 * t526 + t525;
t282 = -qJDD(1) * t302 + qJDD(2);
t136 = -pkin(3) * t204 - pkin(8) * t203 + t282;
t362 = pkin(7) * qJDD(1) + t294;
t246 = t362 * t316;
t247 = t362 * t318;
t387 = qJD(3) * t327;
t388 = qJD(3) * t323;
t141 = -t323 * t246 + t327 * t247 - t279 * t387 - t280 * t388;
t137 = qJDD(3) * pkin(8) + t141;
t53 = -t124 * qJD(4) + t326 * t136 - t137 * t322;
t34 = pkin(4) * t197 - qJ(5) * t154 - qJD(5) * t221 + t53;
t52 = t322 * t136 + t326 * t137 + t170 * t385 - t199 * t386;
t36 = qJ(5) * t155 + qJD(5) * t220 + t52;
t11 = -t315 * t36 + t317 * t34;
t6 = pkin(5) * t197 - pkin(9) * t84 + t11;
t12 = t315 * t34 + t317 * t36;
t7 = pkin(9) * t83 + t12;
t2 = qJD(6) * t13 + t321 * t6 + t325 * t7;
t3 = -qJD(6) * t14 - t321 * t7 + t325 * t6;
t516 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t369 = m(3) * qJ(2) + mrSges(3,3);
t515 = mrSges(2,2) - mrSges(4,3) - t369 - t528;
t313 = pkin(10) + qJ(3);
t305 = sin(t313);
t307 = cos(t313);
t354 = mrSges(4,1) * t307 - mrSges(4,2) * t305;
t355 = -mrSges(3,1) * t318 + mrSges(3,2) * t316;
t514 = m(3) * pkin(1) + t305 * t538 + mrSges(2,1) + t354 - t355;
t513 = t53 * mrSges(5,1) + t11 * mrSges(6,1) - t52 * mrSges(5,2) - t12 * mrSges(6,2);
t30 = qJD(6) * t517 + t321 * t83 + t325 * t84;
t483 = t30 / 0.2e1;
t31 = -qJD(6) * t99 - t321 * t84 + t325 * t83;
t482 = t31 / 0.2e1;
t474 = t83 / 0.2e1;
t473 = t84 / 0.2e1;
t463 = t154 / 0.2e1;
t462 = t155 / 0.2e1;
t194 = qJDD(6) + t197;
t457 = t194 / 0.2e1;
t456 = t197 / 0.2e1;
t288 = t319 * t322;
t289 = t319 * t326;
t213 = t317 * t288 + t289 * t315;
t178 = -pkin(9) * t274 + t213;
t214 = t315 * t288 - t317 * t289;
t179 = -pkin(9) * t342 + t214;
t122 = t178 * t321 + t179 * t325;
t510 = -qJD(6) * t122 - t536 * t321 + t325 * t537;
t121 = t178 * t325 - t179 * t321;
t509 = qJD(6) * t121 + t321 * t537 + t536 * t325;
t436 = pkin(4) * t317;
t301 = pkin(5) + t436;
t437 = pkin(4) * t315;
t258 = t301 * t325 - t321 * t437;
t59 = -t107 * t315 - t400;
t43 = t59 - t511;
t60 = t317 * t107 - t100;
t44 = t60 - t527;
t506 = t258 * qJD(6) - t321 * t43 - t325 * t44;
t259 = t301 * t321 + t325 * t437;
t505 = -t259 * qJD(6) + t321 * t44 - t325 * t43;
t504 = t484 + mrSges(5,1);
t501 = t305 * t522;
t117 = -t173 * t325 + t174 * t321;
t206 = t274 * t325 - t321 * t342;
t146 = -qJD(6) * t206 - t265 * t325 + t267 * t321;
t498 = -t117 + t146;
t118 = -t173 * t321 - t174 * t325;
t205 = -t274 * t321 - t325 * t342;
t145 = qJD(6) * t205 - t265 * t321 - t267 * t325;
t497 = -t118 + t145;
t202 = -pkin(3) * t523 - pkin(8) * t275 - t302;
t212 = -t286 * t323 + t287 * t327;
t209 = t326 * t212;
t153 = t322 * t202 + t209;
t495 = pkin(5) * t493 + t496;
t491 = -t327 * t286 - t287 * t323;
t490 = -t322 * t53 + t326 * t52;
t198 = -qJD(3) * pkin(3) - t207;
t156 = -pkin(4) * t220 + qJD(5) + t198;
t95 = -pkin(5) * t360 + t156;
t489 = -t95 * mrSges(7,1) + t14 * mrSges(7,3);
t488 = t95 * mrSges(7,2) - t13 * mrSges(7,3);
t486 = Ifges(7,4) * t483 + Ifges(7,2) * t482 + Ifges(7,6) * t457;
t485 = Ifges(7,1) * t483 + Ifges(7,4) * t482 + Ifges(7,5) * t457;
t481 = Ifges(6,4) * t473 + Ifges(6,2) * t474 + Ifges(6,6) * t456;
t480 = Ifges(6,1) * t473 + Ifges(6,4) * t474 + Ifges(6,5) * t456;
t439 = Ifges(7,4) * t99;
t49 = Ifges(7,2) * t517 + Ifges(7,6) * t244 + t439;
t479 = -t49 / 0.2e1;
t478 = t49 / 0.2e1;
t94 = Ifges(7,4) * t517;
t50 = Ifges(7,1) * t99 + Ifges(7,5) * t244 + t94;
t477 = -t50 / 0.2e1;
t476 = t50 / 0.2e1;
t475 = Ifges(5,1) * t463 + Ifges(5,4) * t462 + Ifges(5,5) * t456;
t91 = t162 * Ifges(6,4) + Ifges(6,2) * t360 + t252 * Ifges(6,6);
t472 = -t91 / 0.2e1;
t471 = t91 / 0.2e1;
t92 = t162 * Ifges(6,1) + Ifges(6,4) * t360 + t252 * Ifges(6,5);
t470 = -t92 / 0.2e1;
t469 = t92 / 0.2e1;
t468 = -t517 / 0.2e1;
t467 = t517 / 0.2e1;
t466 = -t99 / 0.2e1;
t465 = t99 / 0.2e1;
t461 = -t360 / 0.2e1;
t460 = t360 / 0.2e1;
t459 = -t162 / 0.2e1;
t458 = t162 / 0.2e1;
t455 = -t220 / 0.2e1;
t454 = -t221 / 0.2e1;
t453 = t221 / 0.2e1;
t452 = -t244 / 0.2e1;
t451 = t244 / 0.2e1;
t450 = -t252 / 0.2e1;
t449 = t252 / 0.2e1;
t447 = t266 / 0.2e1;
t445 = t268 / 0.2e1;
t441 = mrSges(6,3) * t54;
t440 = mrSges(6,3) * t55;
t438 = pkin(4) * t221;
t431 = g(3) * t305;
t341 = -qJ(5) * t269 - qJD(5) * t275;
t171 = qJD(2) * t523 + qJD(3) * t491;
t201 = pkin(3) * t270 - pkin(8) * t269;
t361 = -t171 * t322 + t326 * t201;
t63 = pkin(4) * t270 + t341 * t326 + (-t209 + (qJ(5) * t275 - t202) * t322) * qJD(4) + t361;
t374 = t326 * t171 + t322 * t201 + t202 * t385;
t67 = -qJ(5) * t371 + (-qJD(4) * t212 + t341) * t322 + t374;
t26 = t315 * t63 + t317 * t67;
t425 = Ifges(5,4) * t322;
t424 = Ifges(5,4) * t326;
t423 = t123 * mrSges(5,3);
t422 = t124 * mrSges(5,3);
t419 = t221 * Ifges(5,4);
t418 = t268 * Ifges(4,4);
t409 = t275 * t322;
t408 = t275 * t326;
t406 = t307 * t324;
t405 = t307 * t328;
t404 = t308 * t324;
t403 = t308 * t328;
t148 = t220 * Ifges(5,2) + t252 * Ifges(5,6) + t419;
t398 = t322 * t148;
t397 = t322 * t324;
t396 = t322 * t328;
t395 = t324 * t326;
t215 = Ifges(5,4) * t220;
t149 = t221 * Ifges(5,1) + t252 * Ifges(5,5) + t215;
t394 = t326 * t149;
t393 = t326 * t328;
t152 = t326 * t202 - t212 * t322;
t114 = -pkin(4) * t523 - qJ(5) * t408 + t152;
t125 = -qJ(5) * t409 + t153;
t71 = t315 * t114 + t317 * t125;
t222 = t299 * t406 + t300 * t328;
t223 = t299 * t328 - t300 * t406;
t391 = -t222 * mrSges(7,1) + t223 * mrSges(7,2);
t224 = -t299 * t405 + t300 * t324;
t225 = t299 * t324 + t300 * t405;
t390 = t224 * mrSges(7,1) - t225 * mrSges(7,2);
t379 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t194;
t370 = t394 / 0.2e1;
t42 = -t83 * mrSges(6,1) + t84 * mrSges(6,2);
t10 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t366 = -t386 / 0.2e1;
t25 = -t315 * t67 + t317 * t63;
t364 = -t204 * mrSges(4,1) + t203 * mrSges(4,2);
t70 = t317 * t114 - t125 * t315;
t358 = pkin(3) * t307 + pkin(8) * t305;
t175 = pkin(4) * t409 - t491;
t356 = -mrSges(3,1) * t380 + mrSges(3,2) * t381;
t351 = mrSges(5,1) * t322 + mrSges(5,2) * t326;
t350 = -mrSges(7,1) * t299 - mrSges(7,2) * t300;
t349 = Ifges(5,1) * t326 - t425;
t348 = -Ifges(5,2) * t322 + t424;
t347 = Ifges(5,5) * t326 - Ifges(5,6) * t322;
t182 = t342 * t275;
t56 = -pkin(5) * t523 + pkin(9) * t182 + t70;
t181 = t274 * t275;
t58 = -pkin(9) * t181 + t71;
t21 = -t321 * t58 + t325 * t56;
t22 = t321 * t56 + t325 * t58;
t168 = -mrSges(5,2) * t252 + mrSges(5,3) * t220;
t169 = mrSges(5,1) * t252 - mrSges(5,3) * t221;
t345 = t168 * t326 - t169 * t322;
t130 = -t181 * t325 + t182 * t321;
t131 = -t181 * t321 - t182 * t325;
t344 = t281 * t307 - t305 * t310;
t343 = t303 * t307 - t305 * t319;
t142 = -t246 * t327 - t323 * t247 + t279 * t388 - t280 * t387;
t340 = t379 + t516;
t255 = -t307 * t396 + t395;
t253 = t307 * t397 + t393;
t338 = -t269 * t326 + t275 * t386;
t337 = t198 * t351;
t138 = -qJDD(3) * pkin(3) - t142;
t79 = -pkin(4) * t155 + qJDD(5) + t138;
t172 = qJD(2) * t275 + qJD(3) * t212;
t135 = pkin(4) * t339 + t172;
t304 = -qJDD(1) * pkin(1) + qJDD(2);
t256 = t307 * t393 + t397;
t254 = -t307 * t395 + t396;
t248 = Ifges(4,4) * t266;
t236 = pkin(5) * t342 - t303;
t231 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t266;
t230 = t306 * t324 + t307 * t403;
t229 = -t306 * t405 + t404;
t228 = t306 * t328 - t307 * t404;
t227 = t306 * t406 + t403;
t184 = t268 * Ifges(4,1) + t248 + t503;
t183 = t266 * Ifges(4,2) + t418 + t502;
t140 = mrSges(6,1) * t252 - mrSges(6,3) * t162;
t139 = -mrSges(6,2) * t252 + mrSges(6,3) * t360;
t133 = pkin(5) * t162 + t438;
t132 = pkin(5) * t181 + t175;
t127 = -t265 * t275 - t269 * t342;
t126 = t267 * t275 - t269 * t274;
t111 = -mrSges(5,2) * t197 + mrSges(5,3) * t155;
t110 = mrSges(5,1) * t197 - mrSges(5,3) * t154;
t104 = -mrSges(6,1) * t360 + mrSges(6,2) * t162;
t87 = -mrSges(5,1) * t155 + mrSges(5,2) * t154;
t86 = mrSges(7,1) * t244 - mrSges(7,3) * t99;
t85 = -mrSges(7,2) * t244 + mrSges(7,3) * t517;
t78 = -qJD(4) * t153 + t361;
t77 = -t212 * t386 + t374;
t76 = -pkin(5) * t126 + t135;
t74 = t154 * Ifges(5,4) + t155 * Ifges(5,2) + t197 * Ifges(5,6);
t73 = mrSges(6,1) * t197 - mrSges(6,3) * t84;
t72 = -mrSges(6,2) * t197 + mrSges(6,3) * t83;
t57 = -mrSges(7,1) * t517 + mrSges(7,2) * t99;
t47 = -qJD(6) * t131 + t126 * t325 - t127 * t321;
t46 = qJD(6) * t130 + t126 * t321 + t127 * t325;
t41 = -pkin(5) * t83 + t79;
t24 = -mrSges(7,2) * t194 + mrSges(7,3) * t31;
t23 = mrSges(7,1) * t194 - mrSges(7,3) * t30;
t20 = pkin(9) * t126 + t26;
t17 = pkin(5) * t270 - pkin(9) * t127 + t25;
t5 = -qJD(6) * t22 + t17 * t325 - t20 * t321;
t4 = qJD(6) * t21 + t17 * t321 + t20 * t325;
t1 = [m(5) * (t123 * t78 + t124 * t77 + t152 * t53 + t153 * t52) + m(4) * (t141 * t212 + t171 * t208 - t282 * t302) + (-Ifges(6,4) * t182 - Ifges(6,2) * t181) * t474 + (-Ifges(6,5) * t182 - Ifges(6,6) * t181) * t456 + (t282 * mrSges(4,2) - t142 * mrSges(4,3) + t203 * Ifges(4,1) + t204 * Ifges(4,4) + Ifges(4,5) * qJDD(3) + t138 * t351 + t149 * t366 + t347 * t456 + t348 * t462 + t349 * t463) * t275 - (t379 + t524) * t523 / 0.2e1 - (t282 * mrSges(4,1) - t141 * mrSges(4,3) - Ifges(4,4) * t203 + Ifges(5,5) * t463 + Ifges(6,5) * t473 + Ifges(7,5) * t483 - Ifges(4,2) * t204 - Ifges(4,6) * qJDD(3) + Ifges(5,6) * t462 + Ifges(6,6) * t474 + Ifges(7,6) * t482 + Ifges(7,3) * t457 + t526 * t456 + t513 + t516) * t523 + (Ifges(7,5) * t46 + Ifges(7,6) * t47) * t451 + (Ifges(7,5) * t131 + Ifges(7,6) * t130) * t457 + (Ifges(7,4) * t46 + Ifges(7,2) * t47) * t467 + (Ifges(7,4) * t131 + Ifges(7,2) * t130) * t482 + (Ifges(6,4) * t127 + Ifges(6,2) * t126) * t460 + t220 * (-Ifges(5,4) * t338 - Ifges(5,2) * t339) / 0.2e1 + (-t13 * t46 + t130 * t2 - t131 * t3 + t14 * t47) * mrSges(7,3) + (Ifges(7,1) * t46 + Ifges(7,4) * t47) * t465 + (Ifges(7,1) * t131 + Ifges(7,4) * t130) * t483 + (Ifges(6,1) * t127 + Ifges(6,4) * t126) * t458 + (-Ifges(5,1) * t338 - Ifges(5,4) * t339) * t453 - (-m(4) * t142 + m(5) * t138 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t203 + t87) * t491 + (t123 * t338 - t124 * t339 - t408 * t53 - t409 * t52) * mrSges(5,3) + 0.2e1 * t389 * t294 * mrSges(3,3) + (-t254 * mrSges(5,1) - t228 * mrSges(6,1) - t223 * mrSges(7,1) - t253 * mrSges(5,2) - t227 * mrSges(6,2) - t222 * mrSges(7,2) + (t427 * t520 + t515 - t529) * t328 + (-m(6) * (-t302 - t343) - m(5) * (-t302 - t358) - m(7) * (-t302 - t344) + m(4) * t302 + t514) * t324) * g(1) + (t11 * t182 - t12 * t181 + t126 * t55 - t127 * t54) * mrSges(6,3) + (-m(4) * t207 + m(5) * t198 - t494) * t172 + (-Ifges(6,1) * t182 - Ifges(6,4) * t181) * t473 + (-t397 * t484 - t256 * mrSges(5,1) - t230 * mrSges(6,1) - t225 * mrSges(7,1) - t255 * mrSges(5,2) - t229 * mrSges(6,2) - t224 * mrSges(7,2) + t520 * (t328 * t302 + t324 * t427) + t515 * t324 + (-m(5) * t358 - m(6) * t343 - m(7) * t344 - t514) * t328) * g(2) + t79 * (mrSges(6,1) * t181 - mrSges(6,2) * t182) + (Ifges(3,4) * t316 + Ifges(3,2) * t318) * t380 + (Ifges(3,1) * t316 + Ifges(3,4) * t318) * t381 + m(7) * (t13 * t5 + t132 * t41 + t14 * t4 + t2 * t22 + t21 * t3 + t76 * t95) + m(6) * (t11 * t70 + t12 * t71 + t135 * t156 + t175 * t79 + t25 * t54 + t26 * t55) + m(3) * (-pkin(1) * t304 + (t294 + t383) * qJ(2) * t389) - t74 * t409 / 0.2e1 + Ifges(2,3) * qJDD(1) + (-Ifges(5,5) * t338 + Ifges(6,5) * t127 - Ifges(5,6) * t339 + Ifges(6,6) * t126) * t449 - t148 * t371 / 0.2e1 - t302 * t364 + t171 * t231 + t212 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t204) + t175 * t42 + t77 * t168 + t78 * t169 + t156 * (-mrSges(6,1) * t126 + mrSges(6,2) * t127) + t152 * t110 + t153 * t111 + t127 * t469 + t126 * t471 + t408 * t475 + t46 * t476 + t47 * t478 - t182 * t480 - t181 * t481 + t131 * t485 + t130 * t486 + t21 * t23 + t22 * t24 + t198 * (mrSges(5,1) * t339 - mrSges(5,2) * t338) + t71 * t72 + t70 * t73 + t76 * t57 + t4 * t85 + t5 * t86 + t95 * (-mrSges(7,1) * t47 + mrSges(7,2) * t46) + t304 * t355 - pkin(1) * t356 + (t521 / 0.2e1 - Ifges(4,4) * t445 - Ifges(4,2) * t447 + Ifges(7,3) * t451 + Ifges(5,5) * t453 + Ifges(6,5) * t458 + Ifges(6,6) * t460 + Ifges(7,5) * t465 + Ifges(7,6) * t467 + t525 / 0.2e1 + t526 * t449 - t183 / 0.2e1 - t534) * t270 + t41 * (-mrSges(7,1) * t130 + mrSges(7,2) * t131) + t132 * t10 + t135 * t104 + t26 * t139 + t25 * t140 - (-Ifges(4,1) * t445 - Ifges(4,4) * t447 + t398 / 0.2e1 - t184 / 0.2e1 - t370 + t535) * t269; t326 * t110 + t364 + (-t57 - t104 + t494) * t268 - t493 * t140 + t356 + m(3) * t304 - t539 * t139 + t498 * t86 + t497 * t85 + t322 * t111 + t274 * t72 - t342 * t73 + t206 * t24 + t205 * t23 + (-t231 - t345) * t266 + t345 * qJD(4) + (-g(1) * t324 + g(2) * t328) * (m(3) - t520) - t369 * t389 * qJD(1) ^ 2 + (t13 * t498 + t14 * t497 + t2 * t206 + t205 * t3 - t268 * t95) * m(7) + (-t11 * t342 + t12 * t274 - t156 * t268 - t493 * t54 - t539 * t55) * m(6) + (-t198 * t268 + t322 * t52 + t326 * t53 + t252 * (-t123 * t322 + t124 * t326)) * m(5) + (t207 * t268 - t208 * t266 + t282) * m(4); -(Ifges(6,4) * t458 + Ifges(6,2) * t460 + Ifges(6,6) * t449 + t440 + t471) * t265 - (Ifges(6,1) * t458 + Ifges(6,4) * t460 + Ifges(6,5) * t449 - t441 + t469) * t267 + t494 * t208 + t495 * t57 + t496 * t104 + (t337 + t370) * qJD(4) + (m(5) * ((-t123 * t326 - t124 * t322) * qJD(4) + t490) - t169 * t385 - t168 * t386 - t322 * t110 + t326 * t111) * pkin(8) + (t123 * t411 + t124 * t412 + t490) * mrSges(5,3) + (Ifges(7,5) * t118 + Ifges(7,6) * t117) * t452 + (-pkin(3) * t138 - t123 * t143 - t124 * t144 - t198 * t208) * m(5) + t507 * t139 + t508 * t140 + (t11 * t213 + t12 * t214 + t156 * t496 - t303 * t79 + t507 * t55 + t508 * t54) * m(6) + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t468 - (Ifges(4,1) * t266 - t418 + t521) * t268 / 0.2e1 + (Ifges(7,1) * t465 + Ifges(7,4) * t467 + Ifges(7,5) * t451 + t476 + t488) * t145 + (Ifges(7,4) * t465 + Ifges(7,2) * t467 + Ifges(7,6) * t451 + t478 + t489) * t146 + (Ifges(7,1) * t118 + Ifges(7,4) * t117) * t466 + (-t117 * t14 + t118 * t13 + t2 * t205 - t206 * t3) * mrSges(7,3) + t148 * t366 + (-Ifges(6,5) * t174 - Ifges(6,6) * t173) * t450 + (-Ifges(6,4) * t174 - Ifges(6,2) * t173) * t461 + (-Ifges(6,1) * t174 - Ifges(6,4) * t173) * t459 - (-Ifges(4,2) * t268 + t184 + t248 + t394) * t266 / 0.2e1 + t326 * t74 / 0.2e1 + (t220 * t348 + t221 * t349 + t252 * t347) * qJD(4) / 0.2e1 + (mrSges(6,1) * t493 - mrSges(6,2) * t539) * t156 + t509 * t85 + t510 * t86 + (t121 * t3 + t122 * t2 + t13 * t510 + t14 * t509 + t236 * t41 + t495 * t95) * m(7) - t386 * t422 - t385 * t423 - t303 * t42 - t207 * t231 + t236 * t10 + t213 * t73 + t214 * t72 + Ifges(4,3) * qJDD(3) + t41 * (-mrSges(7,1) * t205 + mrSges(7,2) * t206) + Ifges(4,5) * t203 + Ifges(4,6) * t204 - t144 * t168 - t143 * t169 + (-t11 * t274 - t12 * t342 + t173 * t55 - t174 * t54) * mrSges(6,3) + (Ifges(5,5) * t322 + Ifges(6,5) * t274 + Ifges(5,6) * t326 - Ifges(6,6) * t342) * t456 + t79 * (mrSges(6,1) * t342 + mrSges(6,2) * t274) + (Ifges(6,1) * t274 - Ifges(6,4) * t342) * t473 + (Ifges(6,4) * t274 - Ifges(6,2) * t342) * t474 - t342 * t481 + t183 * t445 + t398 * t447 + (Ifges(7,5) * t206 + Ifges(7,6) * t205) * t457 + (Ifges(5,2) * t326 + t425) * t462 + (Ifges(5,1) * t322 + t424) * t463 - t174 * t470 - t173 * t472 + t322 * t475 + t118 * t477 + t117 * t479 + t274 * t480 + (Ifges(7,4) * t206 + Ifges(7,2) * t205) * t482 + (Ifges(7,1) * t206 + Ifges(7,4) * t205) * t483 + t206 * t485 + t205 * t486 - pkin(3) * t87 + t138 * t352 - t95 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) - g(3) * t354 + t121 * t23 + t122 * t24 + (t522 * (mrSges(4,1) - t533) + t532 * g(3)) * t305 + (t522 * (mrSges(4,2) + t532) + t533 * g(3)) * t307 + (Ifges(5,5) * t454 + Ifges(6,5) * t459 + Ifges(7,5) * t466 + Ifges(5,6) * t455 + Ifges(6,6) * t461 + Ifges(7,6) * t468 + Ifges(7,3) * t452 + t526 * t450 + t534) * t268 - t141 * mrSges(4,2) + (t347 * t450 + t348 * t455 + t349 * t454 - t337 + t535) * t266 + t142 * mrSges(4,1); (-Ifges(5,2) * t221 + t149 + t215) * t455 + (t227 * mrSges(6,1) - t228 * mrSges(6,2) - m(7) * (-t283 * t406 - t284 * t328) - t391 - mrSges(5,2) * t254 + t504 * t253) * g(2) + (-t229 * mrSges(6,1) + t230 * mrSges(6,2) - m(7) * (-t283 * t405 + t284 * t324) - t390 + mrSges(5,2) * t256 - t504 * t255) * g(1) + t505 * t86 + t506 * t85 - (Ifges(7,4) * t466 + Ifges(7,2) * t468 + Ifges(7,6) * t452 + t479 - t489) * t99 + (Ifges(7,1) * t466 + Ifges(7,4) * t468 + Ifges(7,5) * t452 + t477 - t488) * t517 + (t13 * t505 - t133 * t95 + t14 * t506 + t2 * t259 + t258 * t3) * m(7) + (mrSges(6,1) * t306 + mrSges(6,2) * t308 - t350 + t351 + t528 + t529) * t431 + (-mrSges(6,2) * t156 + Ifges(6,1) * t459 + Ifges(6,4) * t461 + Ifges(6,5) * t450 + t441 + t470) * t360 + t340 + t524 + t513 - t104 * t438 - m(6) * (t156 * t438 + t54 * t59 + t55 * t60) + t258 * t23 + t259 * t24 - t198 * (mrSges(5,1) * t221 + mrSges(5,2) * t220) - t123 * t168 + t124 * t169 - (mrSges(6,1) * t156 + Ifges(6,4) * t459 + Ifges(6,2) * t461 + Ifges(6,6) * t450 - t440 + t472) * t162 + t221 * t422 + t220 * t423 + t73 * t436 + t72 * t437 + (Ifges(5,5) * t220 - Ifges(5,6) * t221) * t450 + t148 * t453 + (Ifges(5,1) * t220 - t419) * t454 + (t11 * t317 + t12 * t315) * t484 - t133 * t57 - t60 * t139 - t59 * t140; t464 * t307 * g(3) - t360 * t139 + t162 * t140 - t517 * t85 + t99 * t86 + t10 + t42 + (t13 * t99 - t14 * t517 + t41 - t501) * m(7) + (t162 * t54 - t360 * t55 - t501 + t79) * m(6); -t95 * (mrSges(7,1) * t99 + mrSges(7,2) * t517) + (Ifges(7,1) * t517 - t439) * t466 + t49 * t465 + (Ifges(7,5) * t517 - Ifges(7,6) * t99) * t452 - t13 * t85 + t14 * t86 - g(1) * t390 - g(2) * t391 - t350 * t431 + (t13 * t517 + t14 * t99) * mrSges(7,3) + t340 + (-Ifges(7,2) * t99 + t50 + t94) * t468;];
tau  = t1;
