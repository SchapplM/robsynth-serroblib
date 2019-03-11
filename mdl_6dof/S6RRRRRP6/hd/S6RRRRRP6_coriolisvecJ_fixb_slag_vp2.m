% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:35
% EndTime: 2019-03-10 01:27:26
% DurationCPUTime: 25.85s
% Computational Cost: add. (19408->774), mult. (47362->1027), div. (0->0), fcn. (33464->8), ass. (0->337)
t316 = sin(qJ(2));
t319 = cos(qJ(2));
t350 = pkin(2) * t316 - pkin(8) * t319;
t270 = t350 * qJD(1);
t318 = cos(qJ(3));
t315 = sin(qJ(3));
t388 = qJD(1) * t316;
t362 = t315 * t388;
t222 = pkin(7) * t362 + t318 * t270;
t392 = t318 * t319;
t339 = pkin(3) * t316 - pkin(9) * t392;
t465 = -pkin(9) - pkin(8);
t364 = qJD(3) * t465;
t536 = -qJD(1) * t339 + t318 * t364 - t222;
t249 = t315 * t270;
t393 = t316 * t318;
t394 = t315 * t319;
t535 = t249 + (-pkin(7) * t393 - pkin(9) * t394) * qJD(1) - t315 * t364;
t314 = sin(qJ(4));
t317 = cos(qJ(4));
t340 = t314 * t315 - t317 * t318;
t480 = qJD(3) + qJD(4);
t212 = t480 * t340;
t328 = t340 * t319;
t233 = qJD(1) * t328;
t534 = t212 - t233;
t266 = t314 * t318 + t315 * t317;
t213 = t480 * t266;
t329 = t266 * t319;
t232 = qJD(1) * t329;
t531 = t213 - t232;
t283 = t465 * t315;
t284 = t465 * t318;
t380 = qJD(4) * t317;
t381 = qJD(4) * t314;
t514 = t283 * t380 + t284 * t381 + t314 * t536 - t535 * t317;
t219 = t314 * t283 - t317 * t284;
t513 = -qJD(4) * t219 + t535 * t314 + t317 * t536;
t533 = -pkin(10) * t531 + t514;
t532 = -pkin(4) * t388 + pkin(10) * t534 + t513;
t385 = qJD(2) * t318;
t263 = -t362 + t385;
t361 = t318 * t388;
t264 = qJD(2) * t315 + t361;
t202 = t263 * t317 - t264 * t314;
t378 = qJD(2) * qJD(3);
t383 = qJD(3) * t315;
t384 = qJD(2) * t319;
t220 = t318 * t378 + (-t316 * t383 + t318 * t384) * qJD(1);
t382 = qJD(3) * t318;
t511 = t315 * t384 + t316 * t382;
t221 = -qJD(1) * t511 - t315 * t378;
t111 = qJD(4) * t202 + t220 * t317 + t221 * t314;
t203 = t263 * t314 + t264 * t317;
t112 = -qJD(4) * t203 - t220 * t314 + t221 * t317;
t313 = sin(qJ(5));
t433 = cos(qJ(5));
t135 = -t433 * t202 + t313 * t203;
t50 = -qJD(5) * t135 + t111 * t433 + t313 * t112;
t471 = t50 / 0.2e1;
t508 = t202 * t313 + t203 * t433;
t51 = qJD(5) * t508 + t313 * t111 - t112 * t433;
t470 = -t51 / 0.2e1;
t499 = Ifges(6,1) + Ifges(7,1);
t522 = Ifges(6,4) - Ifges(7,5);
t497 = Ifges(7,4) + Ifges(6,5);
t165 = t232 * t433 - t233 * t313;
t205 = t266 * t433 - t313 * t340;
t95 = qJD(5) * t205 - t313 * t212 + t213 * t433;
t402 = t165 - t95;
t166 = -t313 * t232 - t233 * t433;
t332 = -t313 * t266 - t340 * t433;
t94 = qJD(5) * t332 - t212 * t433 - t313 * t213;
t401 = t166 - t94;
t387 = qJD(1) * t319;
t309 = pkin(7) * t387;
t432 = pkin(3) * t315;
t258 = t387 * t432 + t309;
t530 = pkin(3) * t383 - t258;
t78 = pkin(5) * t508 + qJ(6) * t135;
t196 = Ifges(5,4) * t202;
t299 = qJD(3) - t387;
t286 = qJD(4) + t299;
t128 = t203 * Ifges(5,1) + t286 * Ifges(5,5) + t196;
t281 = -qJD(2) * pkin(2) + pkin(7) * t388;
t227 = -t263 * pkin(3) + t281;
t277 = qJD(5) + t286;
t386 = qJD(2) * t316;
t354 = qJD(1) * t386;
t276 = -pkin(2) * t319 - t316 * pkin(8) - pkin(1);
t255 = t276 * qJD(1);
t282 = qJD(2) * pkin(8) + t309;
t209 = t255 * t315 + t282 * t318;
t273 = t350 * qJD(2);
t256 = qJD(1) * t273;
t352 = pkin(7) * t354;
t149 = -qJD(3) * t209 + t318 * t256 + t315 * t352;
t104 = pkin(3) * t354 - pkin(9) * t220 + t149;
t148 = t255 * t382 + t315 * t256 - t282 * t383 - t318 * t352;
t117 = pkin(9) * t221 + t148;
t208 = t318 * t255 - t282 * t315;
t169 = -pkin(9) * t264 + t208;
t159 = pkin(3) * t299 + t169;
t170 = pkin(9) * t263 + t209;
t164 = t317 * t170;
t91 = t314 * t159 + t164;
t29 = -qJD(4) * t91 + t317 * t104 - t117 * t314;
t20 = pkin(4) * t354 - pkin(10) * t111 + t29;
t28 = t314 * t104 + t317 * t117 + t159 * t380 - t170 * t381;
t22 = pkin(10) * t112 + t28;
t355 = qJD(5) * t433;
t379 = qJD(5) * t313;
t523 = pkin(10) * t203;
t162 = t314 * t170;
t90 = t317 * t159 - t162;
t82 = t90 - t523;
t77 = pkin(4) * t286 + t82;
t428 = pkin(10) * t202;
t83 = t91 + t428;
t5 = t313 * t20 + t433 * t22 + t77 * t355 - t379 * t83;
t2 = qJ(6) * t354 + qJD(6) * t277 + t5;
t369 = t433 * t83;
t27 = t313 * t77 + t369;
t6 = -qJD(5) * t27 + t20 * t433 - t313 * t22;
t3 = -pkin(5) * t354 - t6;
t479 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t496 = -Ifges(6,6) + Ifges(7,6);
t526 = Ifges(7,2) + Ifges(6,3);
t483 = t354 * t526 + t496 * t51 + t497 * t50;
t325 = t479 + t483;
t366 = Ifges(5,5) * t111 + Ifges(5,6) * t112 + Ifges(5,3) * t354;
t417 = Ifges(5,4) * t203;
t439 = -t286 / 0.2e1;
t441 = -t277 / 0.2e1;
t451 = -t203 / 0.2e1;
t454 = -t508 / 0.2e1;
t457 = t135 / 0.2e1;
t458 = -t135 / 0.2e1;
t491 = t29 * mrSges(5,1) - t28 * mrSges(5,2);
t502 = -t202 / 0.2e1;
t157 = -pkin(4) * t202 + t227;
t24 = t277 * qJ(6) + t27;
t58 = t135 * pkin(5) - qJ(6) * t508 + t157;
t131 = Ifges(7,5) * t508;
t70 = Ifges(7,6) * t277 + Ifges(7,3) * t135 + t131;
t416 = Ifges(6,4) * t508;
t73 = -Ifges(6,2) * t135 + Ifges(6,6) * t277 + t416;
t507 = mrSges(6,1) * t157 + mrSges(7,1) * t58 - mrSges(7,2) * t24 - mrSges(6,3) * t27 + t70 / 0.2e1 - t73 / 0.2e1;
t403 = t313 * t83;
t26 = t433 * t77 - t403;
t485 = -t26 + qJD(6);
t23 = -t277 * pkin(5) + t485;
t132 = Ifges(6,4) * t135;
t412 = Ifges(7,5) * t135;
t495 = t277 * t497 + t499 * t508 - t132 + t412;
t525 = -mrSges(6,2) * t157 - mrSges(7,2) * t23 + mrSges(6,3) * t26 + mrSges(7,3) * t58 - t495 / 0.2e1;
t529 = -(Ifges(6,4) * t457 + Ifges(7,5) * t458 + t441 * t497 + t454 * t499 + t525) * t135 + t325 + t366 + t491 + (Ifges(5,5) * t202 - Ifges(5,6) * t203) * t439 - t227 * (t203 * mrSges(5,1) + mrSges(5,2) * t202) + (-Ifges(5,2) * t203 + t128 + t196) * t502 + (-Ifges(6,2) * t457 + Ifges(7,3) * t458 + t441 * t496 - t522 * t454 - t507) * t508 + (Ifges(5,1) * t202 - t417) * t451;
t358 = Ifges(3,5) * qJD(2) / 0.2e1;
t528 = t522 * t470 + t499 * t471 + t497 * t354 / 0.2e1;
t218 = t317 * t283 + t284 * t314;
t178 = -pkin(10) * t266 + t218;
t179 = -pkin(10) * t340 + t219;
t334 = t178 * t433 - t313 * t179;
t518 = qJD(5) * t334 + t532 * t313 + t533 * t433;
t123 = t313 * t178 + t179 * t433;
t517 = -qJD(5) * t123 - t533 * t313 + t532 * t433;
t512 = t531 * pkin(4) + t530;
t431 = pkin(4) * t203;
t520 = -pkin(5) * t402 + qJ(6) * t401 - qJD(6) * t205 + t512;
t519 = -qJ(6) * t388 + t518;
t516 = pkin(5) * t388 - t517;
t305 = pkin(3) * t317 + pkin(4);
t99 = -t169 * t314 - t164;
t338 = t99 - t428;
t377 = t313 * t314 * pkin(3);
t100 = t317 * t169 - t162;
t86 = t100 - t523;
t493 = -t313 * t338 + t305 * t355 + (-qJD(4) - qJD(5)) * t377 + (pkin(3) * t380 - t86) * t433;
t31 = t433 * t82 - t403;
t515 = pkin(4) * t355 - t31;
t422 = mrSges(6,3) * t508;
t120 = mrSges(6,1) * t277 - t422;
t121 = -mrSges(7,1) * t277 + mrSges(7,2) * t508;
t487 = -t120 + t121;
t307 = Ifges(3,4) * t387;
t406 = t264 * Ifges(4,4);
t181 = t263 * Ifges(4,2) + t299 * Ifges(4,6) + t406;
t259 = Ifges(4,4) * t263;
t182 = t264 * Ifges(4,1) + t299 * Ifges(4,5) + t259;
t341 = t208 * t318 + t209 * t315;
t418 = Ifges(4,4) * t318;
t345 = -Ifges(4,2) * t315 + t418;
t419 = Ifges(4,4) * t315;
t347 = Ifges(4,1) * t318 - t419;
t348 = mrSges(4,1) * t315 + mrSges(4,2) * t318;
t411 = Ifges(4,6) * t315;
t415 = Ifges(4,5) * t318;
t435 = t318 / 0.2e1;
t436 = -t315 / 0.2e1;
t442 = t264 / 0.2e1;
t321 = -t341 * mrSges(4,3) + t281 * t348 + t263 * t345 / 0.2e1 + t347 * t442 + t299 * (-t411 + t415) / 0.2e1 + t181 * t436 + t182 * t435;
t510 = t321 + Ifges(3,1) * t388 / 0.2e1 + t307 / 0.2e1 + t358;
t505 = -Ifges(6,6) / 0.2e1;
t504 = Ifges(7,6) / 0.2e1;
t469 = t51 / 0.2e1;
t127 = Ifges(5,2) * t202 + Ifges(5,6) * t286 + t417;
t503 = t127 / 0.2e1;
t501 = t202 / 0.2e1;
t440 = t277 / 0.2e1;
t453 = t508 / 0.2e1;
t357 = -Ifges(3,6) * qJD(2) / 0.2e1;
t494 = qJD(6) + t493;
t363 = t433 * t314;
t492 = -t313 * t86 + t338 * t433 + t305 * t379 + (t314 * t355 + (t313 * t317 + t363) * qJD(4)) * pkin(3);
t490 = qJD(6) + t515;
t241 = t340 * t316;
t262 = t318 * t276;
t429 = pkin(7) * t315;
t207 = -pkin(9) * t393 + t262 + (-pkin(3) - t429) * t319;
t302 = pkin(7) * t392;
t231 = t315 * t276 + t302;
t395 = t315 * t316;
t215 = -pkin(9) * t395 + t231;
t143 = t317 * t207 - t314 * t215;
t114 = -pkin(4) * t319 + t241 * pkin(10) + t143;
t144 = t314 * t207 + t317 * t215;
t240 = t266 * t316;
t124 = -pkin(10) * t240 + t144;
t488 = t313 * t114 + t433 * t124;
t486 = Ifges(4,5) * t220 + Ifges(4,6) * t221;
t482 = -t149 * mrSges(4,1) + t148 * mrSges(4,2);
t481 = pkin(1) * mrSges(3,2) * qJD(1);
t150 = -qJD(2) * t328 - t213 * t316;
t389 = t318 * t273 + t386 * t429;
t142 = t339 * qJD(2) + (-t302 + (pkin(9) * t316 - t276) * t315) * qJD(3) + t389;
t167 = t315 * t273 + t276 * t382 + (-t316 * t385 - t319 * t383) * pkin(7);
t147 = -pkin(9) * t511 + t167;
t57 = -qJD(4) * t144 + t317 * t142 - t147 * t314;
t36 = pkin(4) * t386 - pkin(10) * t150 + t57;
t151 = -qJD(2) * t329 + t241 * t480;
t56 = t314 * t142 + t317 * t147 + t207 * t380 - t215 * t381;
t44 = pkin(10) * t151 + t56;
t10 = -qJD(5) * t488 - t313 * t44 + t36 * t433;
t370 = t504 + t505;
t371 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t372 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t420 = Ifges(3,4) * t316;
t478 = t370 * t135 + t371 * t277 + t372 * t508 + t264 * Ifges(4,5) + t263 * Ifges(4,6) + t299 * Ifges(4,3) + Ifges(5,3) * t286 + Ifges(5,6) * t202 + Ifges(5,5) * t203 - t209 * mrSges(4,2) + t208 * mrSges(4,1) - t91 * mrSges(5,2) + t90 * mrSges(5,1) - t23 * mrSges(7,1) + t26 * mrSges(6,1) + t24 * mrSges(7,3) - t27 * mrSges(6,2) + Ifges(7,6) * t457 + Ifges(6,6) * t458 + t357 - (t319 * Ifges(3,2) + t420) * qJD(1) / 0.2e1 + t497 * t453 + t526 * t440;
t473 = Ifges(7,5) * t471 + Ifges(7,3) * t469 + t354 * t504;
t472 = -t50 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t469 + t354 * t505;
t464 = pkin(1) * mrSges(3,1);
t460 = t111 / 0.2e1;
t459 = t112 / 0.2e1;
t450 = t203 / 0.2e1;
t448 = t220 / 0.2e1;
t447 = t221 / 0.2e1;
t446 = -t240 / 0.2e1;
t445 = -t241 / 0.2e1;
t444 = -t263 / 0.2e1;
t443 = -t264 / 0.2e1;
t438 = t286 / 0.2e1;
t437 = -t299 / 0.2e1;
t430 = pkin(4) * t313;
t424 = mrSges(5,3) * t203;
t423 = mrSges(6,3) * t135;
t397 = qJD(2) * mrSges(3,2);
t244 = pkin(3) * t363 + t313 * t305;
t274 = pkin(3) * t395 + t316 * pkin(7);
t374 = t433 * pkin(4);
t365 = Ifges(4,3) * t354 + t486;
t228 = pkin(3) * t511 + pkin(7) * t384;
t306 = -pkin(3) * t318 - pkin(2);
t199 = -pkin(3) * t221 + qJD(2) * t309;
t351 = m(4) * t281 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t263 + mrSges(4,2) * t264 + mrSges(3,3) * t388;
t210 = pkin(4) * t240 + t274;
t161 = pkin(3) * t264 + t431;
t349 = mrSges(4,1) * t318 - mrSges(4,2) * t315;
t346 = Ifges(4,1) * t315 + t418;
t344 = Ifges(4,2) * t318 + t419;
t343 = Ifges(4,5) * t315 + Ifges(4,6) * t318;
t342 = t148 * t318 - t149 * t315;
t234 = pkin(4) * t340 + t306;
t125 = -pkin(4) * t151 + t228;
t42 = -mrSges(7,1) * t354 + t50 * mrSges(7,2);
t87 = -pkin(4) * t112 + t199;
t66 = t114 * t433 - t313 * t124;
t333 = -t240 * t433 + t313 * t241;
t174 = -t313 * t240 - t241 * t433;
t9 = t114 * t355 - t124 * t379 + t313 * t36 + t433 * t44;
t330 = t202 * mrSges(5,3);
t243 = t305 * t433 - t377;
t304 = -t374 - pkin(5);
t303 = qJ(6) + t430;
t279 = mrSges(3,3) * t387 - t397;
t242 = -pkin(5) - t243;
t239 = qJ(6) + t244;
t230 = -pkin(7) * t394 + t262;
t226 = mrSges(4,1) * t299 - mrSges(4,3) * t264;
t225 = -mrSges(4,2) * t299 + mrSges(4,3) * t263;
t223 = -pkin(7) * t361 + t249;
t198 = -mrSges(4,2) * t354 + mrSges(4,3) * t221;
t197 = mrSges(4,1) * t354 - mrSges(4,3) * t220;
t172 = mrSges(5,1) * t286 - t424;
t171 = -t286 * mrSges(5,2) + t330;
t168 = -qJD(3) * t231 + t389;
t160 = -mrSges(4,1) * t221 + mrSges(4,2) * t220;
t153 = t220 * Ifges(4,1) + t221 * Ifges(4,4) + Ifges(4,5) * t354;
t152 = t220 * Ifges(4,4) + t221 * Ifges(4,2) + Ifges(4,6) * t354;
t141 = -mrSges(5,1) * t202 + t203 * mrSges(5,2);
t119 = -mrSges(6,2) * t277 - t423;
t118 = -mrSges(7,2) * t135 + mrSges(7,3) * t277;
t116 = -pkin(5) * t332 - qJ(6) * t205 + t234;
t97 = -mrSges(5,2) * t354 + mrSges(5,3) * t112;
t96 = mrSges(5,1) * t354 - mrSges(5,3) * t111;
t88 = -pkin(5) * t333 - qJ(6) * t174 + t210;
t80 = mrSges(6,1) * t135 + mrSges(6,2) * t508;
t79 = mrSges(7,1) * t135 - mrSges(7,3) * t508;
t69 = qJD(5) * t174 + t313 * t150 - t151 * t433;
t68 = qJD(5) * t333 + t150 * t433 + t313 * t151;
t65 = t431 + t78;
t64 = -mrSges(5,1) * t112 + mrSges(5,2) * t111;
t63 = t319 * pkin(5) - t66;
t62 = -qJ(6) * t319 + t488;
t61 = t111 * Ifges(5,1) + t112 * Ifges(5,4) + Ifges(5,5) * t354;
t60 = t111 * Ifges(5,4) + t112 * Ifges(5,2) + Ifges(5,6) * t354;
t59 = t161 + t78;
t43 = -mrSges(6,2) * t354 - mrSges(6,3) * t51;
t41 = mrSges(6,1) * t354 - mrSges(6,3) * t50;
t40 = -mrSges(7,2) * t51 + mrSges(7,3) * t354;
t30 = t313 * t82 + t369;
t18 = pkin(5) * t69 - qJ(6) * t68 - qJD(6) * t174 + t125;
t17 = mrSges(6,1) * t51 + mrSges(6,2) * t50;
t16 = mrSges(7,1) * t51 - mrSges(7,3) * t50;
t11 = pkin(5) * t51 - qJ(6) * t50 - qJD(6) * t508 + t87;
t8 = -pkin(5) * t386 - t10;
t7 = qJ(6) * t386 - qJD(6) * t319 + t9;
t1 = [(t174 * t499 + t333 * t522) * t471 + (-Ifges(6,2) * t458 + Ifges(7,3) * t457 + t496 * t440 - t522 * t453 + t507) * t69 + m(6) * (t10 * t26 + t125 * t157 + t210 * t87 + t27 * t9 + t488 * t5 + t6 * t66) + t488 * t43 + (t347 * t448 + t345 * t447 + t153 * t435 + pkin(7) * t160 + t152 * t436 + (-t148 * t315 - t149 * t318) * mrSges(4,3) + (t182 * t436 - t318 * t181 / 0.2e1 + t346 * t443 + t343 * t437 + t281 * t349 + t344 * t444 + (t208 * t315 - t209 * t318) * mrSges(4,3)) * qJD(3) + (t357 - pkin(7) * t279 + (-0.2e1 * t464 + Ifges(5,5) * t445 + Ifges(5,6) * t446 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t415 / 0.2e1 - t411 / 0.2e1) * t316 + t372 * t174 - t370 * t333) * qJD(1) + t478) * qJD(2)) * t316 + (Ifges(7,5) * t174 - Ifges(7,3) * t333) * t469 + (-t174 * t6 + t333 * t5) * mrSges(6,3) + (Ifges(6,4) * t174 + Ifges(6,2) * t333) * t470 + (t174 * t3 + t2 * t333) * mrSges(7,2) + t11 * (-mrSges(7,1) * t333 - mrSges(7,3) * t174) + t87 * (-mrSges(6,1) * t333 + mrSges(6,2) * t174) - t333 * t472 - t333 * t473 + m(4) * (t148 * t231 + t149 * t230 + t209 * t167 + t208 * t168) + (-Ifges(5,1) * t241 - Ifges(5,4) * t240) * t460 + (-t150 * t90 + t151 * t91 - t240 * t28 + t241 * t29) * mrSges(5,3) + t199 * (mrSges(5,1) * t240 - mrSges(5,2) * t241) + (-Ifges(5,4) * t241 - Ifges(5,2) * t240) * t459 + m(7) * (t11 * t88 + t18 * t58 + t2 * t62 + t23 * t8 + t24 * t7 + t3 * t63) + m(5) * (t143 * t29 + t144 * t28 + t199 * t274 + t227 * t228 + t56 * t91 + t57 * t90) + (Ifges(6,4) * t458 + Ifges(7,5) * t457 + t497 * t440 + t499 * t453 - t525) * t68 + t274 * t64 + t230 * t197 + t231 * t198 + t167 * t225 + t168 * t226 + t227 * (-mrSges(5,1) * t151 + mrSges(5,2) * t150) + t228 * t141 + t210 * t17 + t56 * t171 + t57 * t172 + t150 * t128 / 0.2e1 + t174 * t528 + (-Ifges(5,5) * t460 - Ifges(5,6) * t459 - Ifges(6,6) * t470 - Ifges(7,6) * t469 - t497 * t471 + (0.3e1 / 0.2e1 * Ifges(3,4) * t384 + (-0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - Ifges(5,3) / 0.2e1 + (m(4) * pkin(7) + t348) * pkin(7) - t371) * t386) * qJD(1) - t479 + t482 - t491) * t319 + (Ifges(5,4) * t150 + Ifges(5,2) * t151) * t501 + t151 * t503 + t62 * t40 + t63 * t42 + t66 * t41 + (t351 * pkin(7) + t358 - 0.2e1 * t481 + t510) * t384 + t18 * t79 - (t366 + t365 + t483 + t486) * t319 / 0.2e1 + t88 * t16 + t7 * t118 + t9 * t119 + t10 * t120 + t8 * t121 + t125 * t80 + t143 * t96 + t144 * t97 + (Ifges(5,5) * t150 + Ifges(5,6) * t151) * t438 + t61 * t445 + t60 * t446 + (Ifges(5,1) * t150 + Ifges(5,4) * t151) * t450; (-Ifges(5,4) * t212 - Ifges(5,2) * t213) * t501 + (-Ifges(5,5) * t212 - Ifges(5,6) * t213) * t438 + (-Ifges(5,1) * t212 - Ifges(5,4) * t213) * t450 + (t232 / 0.2e1 - t213 / 0.2e1) * t127 + (t73 - t70) * (t165 / 0.2e1 - t95 / 0.2e1) - (t42 - t41) * t334 + (t2 * t332 + t205 * t3 - t23 * t401 + t24 * t402) * mrSges(7,2) + (-t205 * t6 + t26 * t401 + t27 * t402 + t332 * t5) * mrSges(6,3) + t11 * (-mrSges(7,1) * t332 - mrSges(7,3) * t205) + t87 * (-mrSges(6,1) * t332 + mrSges(6,2) * t205) + (Ifges(7,5) * t205 - Ifges(7,3) * t332) * t469 + (Ifges(6,4) * t205 + Ifges(6,2) * t332) * t470 - t332 * t472 - t332 * t473 + t199 * (mrSges(5,1) * t340 + mrSges(5,2) * t266) + (Ifges(5,4) * t266 - Ifges(5,2) * t340) * t459 + (Ifges(5,1) * t266 - Ifges(5,4) * t340) * t460 - t340 * t60 / 0.2e1 + t342 * mrSges(4,3) + (-mrSges(6,1) * t402 - mrSges(6,2) * t401) * t157 + (-mrSges(7,1) * t402 + mrSges(7,3) * t401) * t58 + (-t266 * t29 - t28 * t340 - t531 * t91 + t534 * t90) * mrSges(5,3) + (mrSges(5,1) * t531 - mrSges(5,2) * t534) * t227 - m(4) * (t208 * t222 + t209 * t223) + (Ifges(6,4) * t94 + Ifges(7,5) * t166 - Ifges(6,2) * t95 + Ifges(7,3) * t165) * t458 + ((-t225 * t315 - t226 * t318) * qJD(3) + m(4) * (-qJD(3) * t341 + t342) - t197 * t315 + t198 * t318) * pkin(8) + (t40 + t43) * t123 + t315 * t153 / 0.2e1 + t306 * t64 + t266 * t61 / 0.2e1 - t258 * t141 + t234 * t17 - t223 * t225 - t222 * t226 + t218 * t96 + t219 * t97 - pkin(2) * t160 + (Ifges(6,4) * t166 + Ifges(7,5) * t94 - Ifges(6,2) * t165 + Ifges(7,3) * t95) * t457 + t205 * t528 + (t233 / 0.2e1 - t212 / 0.2e1) * t128 + (-Ifges(5,4) * t233 - Ifges(5,2) * t232) * t502 + (-Ifges(5,5) * t233 - Ifges(5,6) * t232) * t439 + (-Ifges(5,1) * t233 - Ifges(5,4) * t232) * t451 + (t141 * t432 + t321) * qJD(3) + (t199 * t306 + t218 * t29 + t219 * t28 + t227 * t530 + t513 * t90 + t514 * t91) * m(5) + ((t358 - t307 / 0.2e1 + t481 + ((-m(4) * pkin(2) - mrSges(3,1) - t349) * qJD(2) - t351) * pkin(7) - t510) * t319 + ((t464 + t420 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t319) * qJD(1) + (t279 + t397) * pkin(7) + t357 - t478) * t316 + (Ifges(5,5) * t266 - Ifges(5,6) * t340 + t205 * t497 - t332 * t496 + t343) * t386 / 0.2e1) * qJD(1) + t512 * t80 + t513 * t172 + t514 * t171 + t116 * t16 + t495 * (-t166 / 0.2e1 + t94 / 0.2e1) + t516 * t121 + t517 * t120 + t518 * t119 + (t123 * t5 + t157 * t512 + t234 * t87 + t26 * t517 + t27 * t518 + t334 * t6) * m(6) + t519 * t118 + t520 * t79 + (t11 * t116 + t123 * t2 + t23 * t516 + t24 * t519 - t334 * t3 + t520 * t58) * m(7) + (t496 * t95 + t497 * t94) * t440 + (t165 * t496 + t166 * t497) * t441 + (t205 * t499 + t332 * t522) * t471 + (t499 * t94 - t522 * t95) * t453 + (-t165 * t522 + t166 * t499) * t454 + t152 * t435 + t344 * t447 + t346 * t448; (t202 * t90 + t203 * t91) * mrSges(5,3) + t203 * t503 + t487 * t492 + (-t157 * t161 + t243 * t6 + t244 * t5 - t26 * t492 + t27 * t493) * m(6) + t493 * t119 + (t2 * t239 + t23 * t492 + t24 * t494 + t242 * t3 - t58 * t59) * m(7) + t494 * t118 + t529 + t365 - t482 - m(5) * (t100 * t91 + t90 * t99) + (t208 * t263 + t209 * t264) * mrSges(4,3) + (-Ifges(4,2) * t264 + t182 + t259) * t444 + (-t264 * t141 + t314 * t97 + t317 * t96 + (t171 * t317 - t172 * t314) * qJD(4) + (0.2e1 * t227 * t443 + t28 * t314 + t29 * t317 + t380 * t91 - t381 * t90) * m(5)) * pkin(3) - t281 * (mrSges(4,1) * t264 + mrSges(4,2) * t263) + t242 * t42 + t243 * t41 + t244 * t43 + t239 * t40 - t208 * t225 + t209 * t226 - t100 * t171 - t99 * t172 - t161 * t80 - t59 * t79 + (Ifges(4,5) * t263 - Ifges(4,6) * t264) * t437 + t181 * t442 + (Ifges(4,1) * t263 - t406) * t443; (m(7) * t23 + t487) * (pkin(4) * t379 - t30) + t490 * t118 + (-t171 + t330) * t90 + t515 * t119 - t80 * t431 + (-t157 * t431 + t26 * t30 - t27 * t31 + (t433 * t6 + t313 * t5 + (-t26 * t313 + t27 * t433) * qJD(5)) * pkin(4)) * m(6) + (t2 * t303 + t24 * t490 + t3 * t304 - t58 * t65) * m(7) + t303 * t40 + t304 * t42 + (t172 + t424) * t91 - t65 * t79 + t41 * t374 + t43 * t430 + t127 * t450 + t529; t325 + (t422 - t487) * t27 + (-t118 - t119 - t423) * t26 + (t135 * t23 + t24 * t508) * mrSges(7,2) - t157 * (mrSges(6,1) * t508 - mrSges(6,2) * t135) + qJ(6) * t40 - pkin(5) * t42 - t78 * t79 + qJD(6) * t118 - t58 * (mrSges(7,1) * t508 + mrSges(7,3) * t135) + (Ifges(7,3) * t508 - t412) * t458 + t73 * t453 + (-t135 * t497 + t496 * t508) * t441 + (-pkin(5) * t3 + qJ(6) * t2 - t23 * t27 + t24 * t485 - t58 * t78) * m(7) + (-Ifges(6,2) * t508 - t132 + t495) * t457 + (-t135 * t499 + t131 - t416 + t70) * t454; -t277 * t118 + t508 * t79 + 0.2e1 * (t3 / 0.2e1 + t58 * t453 + t24 * t441) * m(7) + t42;];
tauc  = t1(:);
