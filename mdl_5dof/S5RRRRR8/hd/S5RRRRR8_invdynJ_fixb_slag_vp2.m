% Calculate vector of inverse dynamics joint torques for
% S5RRRRR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:24
% EndTime: 2019-12-31 22:24:55
% DurationCPUTime: 15.64s
% Computational Cost: add. (11734->696), mult. (25594->925), div. (0->0), fcn. (18190->14), ass. (0->343)
t551 = -mrSges(5,3) - mrSges(6,3);
t315 = sin(qJ(3));
t320 = cos(qJ(3));
t321 = cos(qJ(2));
t396 = qJD(1) * t321;
t316 = sin(qJ(2));
t397 = qJD(1) * t316;
t240 = -t315 * t397 + t320 * t396;
t256 = t315 * t321 + t316 * t320;
t241 = t256 * qJD(1);
t187 = pkin(3) * t241 - pkin(8) * t240;
t380 = pkin(2) * t397;
t167 = t187 + t380;
t324 = -pkin(7) - pkin(6);
t279 = t324 * t321;
t259 = qJD(1) * t279;
t242 = t315 * t259;
t277 = t324 * t316;
t258 = qJD(1) * t277;
t194 = t258 * t320 + t242;
t314 = sin(qJ(4));
t319 = cos(qJ(4));
t109 = t314 * t167 + t319 * t194;
t393 = qJD(3) * t320;
t377 = pkin(2) * t393;
t554 = t319 * t377 - t109;
t108 = t319 * t167 - t194 * t314;
t553 = -t314 * t377 - t108;
t392 = qJD(4) * t314;
t423 = t240 * t314;
t552 = t392 - t423;
t465 = pkin(2) * t315;
t297 = pkin(8) + t465;
t449 = -pkin(9) - t297;
t363 = qJD(4) * t449;
t387 = pkin(9) * t423;
t550 = t314 * t363 + t387 + t554;
t422 = t240 * t319;
t354 = t241 * pkin(4) - pkin(9) * t422;
t549 = t319 * t363 - t354 + t553;
t246 = qJD(2) * pkin(2) + t258;
t190 = t246 * t320 + t242;
t112 = t314 * t187 + t319 * t190;
t323 = -pkin(9) - pkin(8);
t371 = qJD(4) * t323;
t548 = t314 * t371 - t112 + t387;
t111 = t319 * t187 - t190 * t314;
t547 = t319 * t371 - t111 - t354;
t312 = qJ(2) + qJ(3);
t305 = sin(t312);
t311 = qJ(4) + qJ(5);
t304 = sin(t311);
t444 = mrSges(6,2) * t304;
t445 = mrSges(5,2) * t314;
t546 = (-t444 - t445) * t305;
t318 = cos(qJ(5));
t313 = sin(qJ(5));
t310 = qJD(2) + qJD(3);
t207 = -t241 * t314 + t310 * t319;
t462 = pkin(2) * t321;
t300 = pkin(1) + t462;
t275 = t300 * qJD(1);
t162 = -pkin(3) * t240 - pkin(8) * t241 - t275;
t243 = t320 * t259;
t191 = t246 * t315 - t243;
t172 = pkin(8) * t310 + t191;
t94 = t162 * t314 + t172 * t319;
t70 = pkin(9) * t207 + t94;
t429 = t313 * t70;
t236 = qJD(4) - t240;
t208 = t241 * t319 + t310 * t314;
t93 = t319 * t162 - t172 * t314;
t69 = -pkin(9) * t208 + t93;
t61 = pkin(4) * t236 + t69;
t23 = t318 * t61 - t429;
t427 = t318 * t70;
t24 = t313 * t61 + t427;
t526 = t310 * Ifges(4,6);
t545 = t275 * mrSges(4,1) - t93 * mrSges(5,1) - t23 * mrSges(6,1) + t94 * mrSges(5,2) + t24 * mrSges(6,2) + t526 / 0.2e1;
t527 = t310 * Ifges(4,5);
t544 = t275 * mrSges(4,2) + t190 * mrSges(4,3) - t527 / 0.2e1;
t495 = m(6) * pkin(4);
t130 = t207 * t313 + t208 * t318;
t228 = qJD(5) + t236;
t359 = t318 * t207 - t208 * t313;
t528 = t236 * Ifges(5,3);
t529 = t207 * Ifges(5,6);
t543 = t208 * Ifges(5,5) + t130 * Ifges(6,5) + Ifges(6,6) * t359 + t228 * Ifges(6,3) + t528 + t529;
t431 = t241 * mrSges(4,3);
t520 = mrSges(4,1) * t310 + mrSges(5,1) * t207 - mrSges(5,2) * t208 - t431;
t390 = qJD(1) * qJD(2);
t265 = qJDD(1) * t321 - t316 * t390;
t306 = cos(t311);
t447 = mrSges(6,1) * t306;
t448 = mrSges(5,1) * t319;
t542 = (t447 + t448) * t305;
t541 = t552 * pkin(4);
t193 = t258 * t315 - t243;
t394 = qJD(3) * t315;
t540 = pkin(2) * t394 - t193;
t391 = qJD(4) * t319;
t266 = qJDD(1) * t316 + t321 * t390;
t254 = t315 * t316 - t320 * t321;
t332 = t254 * qJD(3);
t152 = -qJD(1) * t332 + t265 * t315 + t266 * t320;
t153 = -qJD(3) * t241 + t265 * t320 - t315 * t266;
t426 = qJDD(1) * pkin(1);
t233 = -pkin(2) * t265 - t426;
t71 = -pkin(3) * t153 - pkin(8) * t152 + t233;
t309 = qJDD(2) + qJDD(3);
t250 = t266 * pkin(6);
t202 = qJDD(2) * pkin(2) - pkin(7) * t266 - t250;
t249 = t265 * pkin(6);
t206 = pkin(7) * t265 + t249;
t80 = t315 * t202 + t320 * t206 + t246 * t393 + t259 * t394;
t77 = pkin(8) * t309 + t80;
t15 = t162 * t391 - t172 * t392 + t314 * t71 + t319 * t77;
t16 = -qJD(4) * t94 - t314 * t77 + t319 * t71;
t539 = t15 * t319 - t16 * t314;
t538 = t445 - t448;
t307 = cos(t312);
t537 = -t307 * mrSges(4,1) + (mrSges(4,2) + t551) * t305;
t515 = t307 * pkin(3) + t305 * pkin(8);
t298 = pkin(4) * t319 + pkin(3);
t517 = t307 * t298 - t305 * t323;
t535 = -m(5) * t515 - m(6) * t517;
t97 = qJD(4) * t207 + t152 * t319 + t309 * t314;
t98 = -qJD(4) * t208 - t152 * t314 + t309 * t319;
t28 = qJD(5) * t359 + t313 * t98 + t318 * t97;
t493 = t28 / 0.2e1;
t29 = -qJD(5) * t130 - t313 * t97 + t318 * t98;
t492 = t29 / 0.2e1;
t486 = t97 / 0.2e1;
t485 = t98 / 0.2e1;
t151 = qJDD(4) - t153;
t146 = qJDD(5) + t151;
t480 = t146 / 0.2e1;
t479 = t151 / 0.2e1;
t533 = t265 / 0.2e1;
t466 = t321 / 0.2e1;
t532 = t314 * t495;
t247 = t449 * t314;
t308 = t319 * pkin(9);
t416 = t297 * t319;
t248 = t308 + t416;
t184 = t247 * t313 + t248 * t318;
t531 = -qJD(5) * t184 - t313 * t550 + t318 * t549;
t183 = t247 * t318 - t248 * t313;
t530 = qJD(5) * t183 + t313 * t549 + t318 * t550;
t525 = t321 * Ifges(3,2);
t524 = mrSges(5,1) + t495;
t276 = t323 * t314;
t456 = pkin(8) * t319;
t278 = t308 + t456;
t211 = t276 * t313 + t278 * t318;
t523 = -qJD(5) * t211 - t313 * t548 + t318 * t547;
t209 = t276 * t318 - t278 * t313;
t522 = qJD(5) * t209 + t313 * t547 + t318 * t548;
t171 = -pkin(3) * t310 - t190;
t349 = mrSges(5,1) * t314 + mrSges(5,2) * t319;
t521 = t171 * t349;
t340 = t313 * t314 - t318 * t319;
t177 = t340 * t256;
t189 = pkin(3) * t254 - pkin(8) * t256 - t300;
t212 = t277 * t315 - t279 * t320;
t204 = t319 * t212;
t121 = t314 * t189 + t204;
t519 = t320 * t277 + t279 * t315;
t518 = t540 + t541;
t341 = -t298 * t305 - t307 * t323;
t461 = pkin(3) * t305;
t464 = pkin(2) * t316;
t514 = -m(6) * (t341 - t464) - m(5) * (-t461 - t464) + t542;
t513 = -m(6) * t341 + t542;
t512 = -t191 + t541;
t457 = pkin(6) * t321;
t458 = pkin(6) * t316;
t510 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t397) * t457 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t396) * t458;
t322 = cos(qJ(1));
t413 = t307 * t322;
t509 = t546 * t322 + t413 * t551;
t317 = sin(qJ(1));
t414 = t307 * t317;
t508 = t546 * t317 + t414 * t551;
t507 = t249 * t321 + t250 * t316;
t506 = g(1) * t322 + g(2) * t317;
t505 = -m(5) - m(6) - m(4);
t504 = qJD(4) + qJD(5);
t503 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t502 = t537 + (t444 - t447 + t538) * t307;
t12 = pkin(4) * t151 - pkin(9) * t97 + t16;
t13 = pkin(9) * t98 + t15;
t3 = qJD(5) * t23 + t12 * t313 + t13 * t318;
t4 = -qJD(5) * t24 + t12 * t318 - t13 * t313;
t501 = t4 * mrSges(6,1) - t3 * mrSges(6,2);
t500 = t16 * mrSges(5,1) - t15 * mrSges(5,2);
t274 = -mrSges(3,1) * t321 + mrSges(3,2) * t316;
t499 = m(3) * pkin(1) + mrSges(2,1) - t274 - t537;
t443 = mrSges(5,3) * t207;
t155 = -mrSges(5,2) * t236 + t443;
t442 = mrSges(5,3) * t208;
t156 = mrSges(5,1) * t236 - t442;
t59 = mrSges(5,1) * t151 - mrSges(5,3) * t97;
t498 = m(5) * ((-t314 * t94 - t319 * t93) * qJD(4) + t539) - t156 * t391 - t155 * t392 - t314 * t59;
t496 = Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t480;
t494 = Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t480;
t491 = Ifges(5,1) * t486 + Ifges(5,4) * t485 + Ifges(5,5) * t479;
t436 = Ifges(6,4) * t130;
t57 = Ifges(6,2) * t359 + t228 * Ifges(6,6) + t436;
t490 = -t57 / 0.2e1;
t489 = t57 / 0.2e1;
t125 = Ifges(6,4) * t359;
t58 = t130 * Ifges(6,1) + t228 * Ifges(6,5) + t125;
t488 = -t58 / 0.2e1;
t487 = t58 / 0.2e1;
t484 = -t359 / 0.2e1;
t483 = t359 / 0.2e1;
t482 = -t130 / 0.2e1;
t481 = t130 / 0.2e1;
t477 = -t207 / 0.2e1;
t476 = -t208 / 0.2e1;
t475 = t208 / 0.2e1;
t474 = -t228 / 0.2e1;
t473 = t228 / 0.2e1;
t472 = -t236 / 0.2e1;
t470 = t240 / 0.2e1;
t468 = t241 / 0.2e1;
t463 = pkin(2) * t320;
t460 = pkin(4) * t208;
t453 = g(3) * t305;
t452 = t23 * mrSges(6,3);
t451 = t24 * mrSges(6,3);
t441 = Ifges(3,4) * t316;
t440 = Ifges(3,4) * t321;
t439 = Ifges(5,4) * t208;
t438 = Ifges(5,4) * t314;
t437 = Ifges(5,4) * t319;
t435 = pkin(4) * qJD(5);
t430 = t241 * Ifges(4,4);
t198 = -qJD(2) * t254 - t332;
t425 = t198 * t319;
t419 = t256 * t314;
t418 = t256 * t319;
t114 = t207 * Ifges(5,2) + t236 * Ifges(5,6) + t439;
t412 = t314 * t114;
t411 = t314 * t317;
t410 = t314 * t322;
t407 = t317 * t319;
t205 = Ifges(5,4) * t207;
t115 = Ifges(5,1) * t208 + Ifges(5,5) * t236 + t205;
t406 = t319 * t115;
t405 = t319 * t322;
t221 = t304 * t414 + t306 * t322;
t222 = t304 * t322 - t306 * t414;
t404 = -t221 * mrSges(6,1) + t222 * mrSges(6,2);
t223 = -t304 * t413 + t306 * t317;
t224 = t304 * t317 + t306 * t413;
t403 = t223 * mrSges(6,1) - t224 * mrSges(6,2);
t395 = qJD(2) * t316;
t384 = Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * t146;
t383 = Ifges(5,5) * t97 + Ifges(5,6) * t98 + Ifges(5,3) * t151;
t379 = pkin(2) * t395;
t372 = qJD(2) * t324;
t370 = t256 * t391;
t367 = t406 / 0.2e1;
t364 = -t392 / 0.2e1;
t362 = t390 / 0.2e1;
t199 = t310 * t256;
t118 = pkin(3) * t199 - pkin(8) * t198 + t379;
t263 = t316 * t372;
t264 = t321 * t372;
t133 = qJD(3) * t519 + t263 * t320 + t264 * t315;
t360 = t319 * t118 - t133 * t314;
t120 = t319 * t189 - t212 * t314;
t352 = mrSges(3,1) * t316 + mrSges(3,2) * t321;
t350 = mrSges(4,1) * t305 + mrSges(4,2) * t307;
t348 = -mrSges(6,1) * t304 - mrSges(6,2) * t306;
t347 = Ifges(5,1) * t319 - t438;
t346 = t441 + t525;
t345 = -Ifges(5,2) * t314 + t437;
t344 = Ifges(3,5) * t321 - Ifges(3,6) * t316;
t343 = Ifges(5,5) * t319 - Ifges(5,6) * t314;
t102 = -pkin(9) * t419 + t121;
t85 = pkin(4) * t254 - pkin(9) * t418 + t120;
t47 = t102 * t318 + t313 * t85;
t46 = -t102 * t313 + t318 * t85;
t255 = t313 * t319 + t314 * t318;
t81 = t202 * t320 - t315 * t206 - t246 * t394 + t259 * t393;
t337 = t384 + t501;
t336 = pkin(1) * t352;
t231 = -t307 * t410 + t407;
t229 = t307 * t411 + t405;
t335 = t198 * t314 + t370;
t334 = t256 * t392 - t425;
t333 = t316 * (Ifges(3,1) * t321 - t441);
t41 = t314 * t118 + t319 * t133 + t189 * t391 - t212 * t392;
t78 = -pkin(3) * t309 - t81;
t134 = qJD(3) * t212 + t263 * t315 - t320 * t264;
t197 = t504 * t255;
t119 = -pkin(4) * t207 + t171;
t160 = t255 * t240;
t161 = t340 * t240;
t168 = t240 * Ifges(4,2) + t430 + t526;
t234 = Ifges(4,4) * t240;
t169 = t241 * Ifges(4,1) + t234 + t527;
t196 = t504 * t340;
t34 = t97 * Ifges(5,4) + t98 * Ifges(5,2) + t151 * Ifges(5,6);
t43 = -pkin(4) * t98 + t78;
t325 = (t343 * t472 + t345 * t477 + t347 * t476 - t521 + t544) * t240 + (Ifges(5,5) * t476 + Ifges(6,5) * t482 + Ifges(5,6) * t477 + Ifges(6,6) * t484 + Ifges(5,3) * t472 + Ifges(6,3) * t474 + t545) * t241 + (t207 * t345 + t208 * t347 + t236 * t343) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t241 + t169 + t234 + t406) * t240 / 0.2e1 + (t160 * t24 - t161 * t23 - t255 * t4 - t3 * t340) * mrSges(6,3) + (Ifges(6,4) * t255 - Ifges(6,2) * t340) * t492 + (Ifges(6,1) * t255 - Ifges(6,4) * t340) * t493 + t43 * (mrSges(6,1) * t340 + mrSges(6,2) * t255) + (Ifges(6,5) * t255 - Ifges(6,6) * t340) * t480 - t340 * t496 + ((t161 - t196) * mrSges(6,2) + (-t160 + t197) * mrSges(6,1)) * t119 + (-Ifges(6,5) * t161 - Ifges(6,6) * t160) * t474 + (t521 + t367) * qJD(4) + (-Ifges(6,4) * t161 - Ifges(6,2) * t160) * t484 + (-Ifges(6,1) * t161 - Ifges(6,4) * t160) * t482 - (Ifges(4,1) * t240 - t430 + t543) * t241 / 0.2e1 + t78 * t538 + (-t552 * t94 + (-t391 + t422) * t93 + t539) * mrSges(5,3) - (Ifges(6,4) * t481 + Ifges(6,2) * t483 + Ifges(6,6) * t473 + t451 + t489) * t197 - (Ifges(6,1) * t481 + Ifges(6,4) * t483 + Ifges(6,5) * t473 - t452 + t487) * t196 + t319 * t34 / 0.2e1 + Ifges(4,3) * t309 + t191 * t431 + t114 * t364 + t168 * t468 + t412 * t470 + (Ifges(5,5) * t314 + Ifges(5,6) * t319) * t479 + (Ifges(5,2) * t319 + t438) * t485 + (Ifges(5,1) * t314 + t437) * t486 - t161 * t488 - t160 * t490 + t314 * t491 + t255 * t494 - t80 * mrSges(4,2) + t81 * mrSges(4,1) + Ifges(4,5) * t152 + Ifges(4,6) * t153;
t302 = Ifges(3,4) * t396;
t299 = -pkin(3) - t463;
t289 = pkin(8) * t413;
t288 = pkin(8) * t414;
t271 = -t298 - t463;
t239 = Ifges(3,1) * t397 + Ifges(3,5) * qJD(2) + t302;
t238 = Ifges(3,6) * qJD(2) + qJD(1) * t346;
t232 = t307 * t405 + t411;
t230 = -t307 * t407 + t410;
t215 = -mrSges(4,2) * t310 + mrSges(4,3) * t240;
t186 = -mrSges(4,1) * t240 + mrSges(4,2) * t241;
t176 = t255 * t256;
t163 = pkin(4) * t419 - t519;
t139 = -mrSges(4,2) * t309 + mrSges(4,3) * t153;
t138 = mrSges(4,1) * t309 - mrSges(4,3) * t152;
t105 = mrSges(6,1) * t228 - mrSges(6,3) * t130;
t104 = -mrSges(6,2) * t228 + mrSges(6,3) * t359;
t75 = pkin(4) * t335 + t134;
t66 = -mrSges(6,1) * t359 + mrSges(6,2) * t130;
t60 = -mrSges(5,2) * t151 + mrSges(5,3) * t98;
t55 = t177 * t504 - t255 * t198;
t54 = -t197 * t256 - t198 * t340;
t48 = -mrSges(5,1) * t98 + mrSges(5,2) * t97;
t42 = -qJD(4) * t121 + t360;
t31 = t318 * t69 - t429;
t30 = -t313 * t69 - t427;
t25 = -pkin(9) * t335 + t41;
t20 = -pkin(9) * t425 + pkin(4) * t199 + (-t204 + (pkin(9) * t256 - t189) * t314) * qJD(4) + t360;
t18 = -mrSges(6,2) * t146 + mrSges(6,3) * t29;
t17 = mrSges(6,1) * t146 - mrSges(6,3) * t28;
t11 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t6 = -qJD(5) * t47 + t20 * t318 - t25 * t313;
t5 = qJD(5) * t46 + t20 * t313 + t25 * t318;
t1 = [-(-m(4) * t81 + m(5) * t78 - t138 + t48) * t519 + (-t412 / 0.2e1 + t367 + Ifges(4,1) * t468 + Ifges(4,4) * t470 + t169 / 0.2e1 - t544) * t198 + (t528 / 0.2e1 + t529 / 0.2e1 - mrSges(4,3) * t191 - Ifges(4,4) * t468 - Ifges(4,2) * t470 + Ifges(6,3) * t473 + Ifges(5,5) * t475 + Ifges(6,5) * t481 + Ifges(6,6) * t483 - t168 / 0.2e1 + t543 / 0.2e1 - t545) * t199 + (-t230 * mrSges(5,1) - t222 * mrSges(6,1) - t229 * mrSges(5,2) - t221 * mrSges(6,2) + (-t324 * t505 + t503 - t532) * t322 + (m(4) * t300 - m(6) * (-t300 - t517) - m(5) * (-t300 - t515) + t499) * t317) * g(1) + (-t411 * t495 - t232 * mrSges(5,1) - t224 * mrSges(6,1) - t231 * mrSges(5,2) - t223 * mrSges(6,2) + t505 * (t322 * t300 - t317 * t324) + t503 * t317 + (-t499 + t535) * t322) * g(2) + t171 * (mrSges(5,1) * t335 - mrSges(5,2) * t334) + m(4) * (t133 * t191 + t212 * t80 - t233 * t300 - t275 * t379) + t266 * t440 / 0.2e1 + (Ifges(6,5) * t54 + Ifges(6,6) * t55) * t473 + (t321 * t440 + t333) * t362 + (Ifges(6,1) * t54 + Ifges(6,4) * t55) * t481 + m(5) * (t120 * t16 + t121 * t15 + t41 * t94 + t42 * t93) + (Ifges(3,4) * t266 + Ifges(3,2) * t265) * t466 + (-Ifges(5,1) * t334 - Ifges(5,4) * t335) * t475 + t346 * t533 + t236 * (-Ifges(5,5) * t334 - Ifges(5,6) * t335) / 0.2e1 + t186 * t379 + t207 * (-Ifges(5,4) * t334 - Ifges(5,2) * t335) / 0.2e1 + (-t15 * t419 - t16 * t418 + t334 * t93 - t335 * t94) * mrSges(5,3) + m(6) * (t119 * t75 + t163 * t43 + t23 * t6 + t24 * t5 + t3 * t47 + t4 * t46) - t300 * (-mrSges(4,1) * t153 + mrSges(4,2) * t152) + Ifges(2,3) * qJDD(1) - pkin(1) * (-mrSges(3,1) * t265 + mrSges(3,2) * t266) - t274 * t426 - t34 * t419 / 0.2e1 + t212 * t139 + t133 * t215 + (mrSges(4,1) * t233 - mrSges(4,3) * t80 - Ifges(4,4) * t152 + Ifges(5,5) * t486 + Ifges(6,5) * t493 - Ifges(4,2) * t153 - Ifges(4,6) * t309 + Ifges(5,6) * t485 + Ifges(6,6) * t492 + Ifges(5,3) * t479 + Ifges(6,3) * t480 + t500 + t501) * t254 + (Ifges(6,4) * t54 + Ifges(6,2) * t55) * t483 + (t384 + t383) * t254 / 0.2e1 + (-Ifges(6,5) * t177 - Ifges(6,6) * t176) * t480 + (-Ifges(6,1) * t177 - Ifges(6,4) * t176) * t493 + (-t176 * t3 + t177 * t4 - t23 * t54 + t24 * t55) * mrSges(6,3) + (-Ifges(6,4) * t177 - Ifges(6,2) * t176) * t492 + t43 * (mrSges(6,1) * t176 - mrSges(6,2) * t177) + (t265 * t457 + t266 * t458 + t507) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t507) + (t239 * t466 + t344 * qJD(2) / 0.2e1 - t510) * qJD(2) + (-m(4) * t190 + m(5) * t171 - t520) * t134 + t54 * t487 + t55 * t489 + t418 * t491 - t177 * t494 - t176 * t496 + (-mrSges(3,1) * t458 - mrSges(3,2) * t457 + 0.2e1 * Ifges(3,6) * t466) * qJDD(2) + (Ifges(3,1) * t266 + Ifges(3,4) * t533 + Ifges(3,5) * qJDD(2) - t362 * t525) * t316 + t46 * t17 + t47 * t18 - t114 * t370 / 0.2e1 + t75 * t66 + t5 * t104 + t6 * t105 + (t233 * mrSges(4,2) - t81 * mrSges(4,3) + Ifges(4,1) * t152 + Ifges(4,4) * t153 + Ifges(4,5) * t309 + t115 * t364 + t343 * t479 + t345 * t485 + t347 * t486 + t349 * t78) * t256 + t119 * (-mrSges(6,1) * t55 + mrSges(6,2) * t54) + t120 * t59 + t121 * t60 + t41 * t155 + t42 * t156 + t163 * t11 - t336 * t390 - t238 * t395 / 0.2e1; (-m(4) * t462 - m(6) * (t517 + t462) - m(5) * (t515 + t462) + t274 + t502) * g(3) - t520 * t540 + (m(4) * t464 + t350 + t352) * t506 + (t190 * t193 - t191 * t194 + t275 * t380 + (t315 * t80 + t320 * t81 + (-t190 * t315 + t191 * t320) * qJD(3)) * pkin(2)) * m(4) + (t510 + (-t333 / 0.2e1 + t336) * qJD(1)) * qJD(1) + t60 * t416 + (t299 * t78 + (t171 * t315 + (-t314 * t93 + t319 * t94) * t320) * qJD(3) * pkin(2) - t108 * t93 - t109 * t94 - t171 * t193 - g(1) * t289 - g(2) * t288) * m(5) + t553 * t156 + t554 * t155 + t325 + t299 * t48 + Ifges(3,6) * t265 + Ifges(3,5) * t266 + t271 * t11 - t249 * mrSges(3,2) - t250 * mrSges(3,1) + t498 * t297 - (-Ifges(3,2) * t397 + t239 + t302) * t396 / 0.2e1 + (t317 * t514 + t508) * g(2) + (t322 * t514 + t509) * g(1) + t518 * t66 + t138 * t463 + t139 * t465 + t530 * t104 + t531 * t105 + (t518 * t119 + t183 * t4 + t184 * t3 + t531 * t23 + t530 * t24 + t271 * t43) * m(6) + Ifges(3,3) * qJDD(2) + (t377 - t194) * t215 - t186 * t380 + t183 * t17 + t184 * t18 - t344 * t390 / 0.2e1 + t238 * t397 / 0.2e1; -pkin(3) * t48 - t298 * t11 - t111 * t156 - t112 * t155 + t209 * t17 + t211 * t18 - t190 * t215 + t60 * t456 + t325 + t512 * t66 + t506 * t350 + t520 * t191 + t523 * t105 + t522 * t104 + (t317 * t513 + t508) * g(2) + (t322 * t513 + t509) * g(1) + (t119 * t512 + t209 * t4 + t211 * t3 + t23 * t523 + t24 * t522 - t298 * t43) * m(6) + (-t111 * t93 - t112 * t94 - t171 * t191 - pkin(3) * t78 - g(1) * (-t322 * t461 + t289) - g(2) * (-t317 * t461 + t288)) * m(5) + (t502 + t535) * g(3) + t498 * pkin(8); (t443 - t155) * t93 + (t318 * t435 - t31) * t104 + t500 + t383 + t337 + (t318 * t17 + t313 * t18) * pkin(4) + (mrSges(5,2) * t232 - t231 * t524 - t403) * g(1) + (-t313 * t435 - t30) * t105 + (-mrSges(5,2) * t230 + t229 * t524 - t404) * g(2) + (-Ifges(5,2) * t208 + t115 + t205) * t477 + (t442 + t156) * t94 + (-t348 + t349 + t532) * t453 - t66 * t460 - m(6) * (t119 * t460 + t23 * t30 + t24 * t31) - t171 * (mrSges(5,1) * t208 + mrSges(5,2) * t207) - (t119 * mrSges(6,1) + Ifges(6,4) * t482 + Ifges(6,2) * t484 + Ifges(6,6) * t474 - t451 + t490) * t130 + (-t119 * mrSges(6,2) + Ifges(6,1) * t482 + Ifges(6,4) * t484 + Ifges(6,5) * t474 + t452 + t488) * t359 + (Ifges(5,5) * t207 - Ifges(5,6) * t208) * t472 + t114 * t475 + (Ifges(5,1) * t207 - t439) * t476 + (t3 * t313 + t318 * t4 + (-t23 * t313 + t24 * t318) * qJD(5)) * t495; -t119 * (mrSges(6,1) * t130 + mrSges(6,2) * t359) + (Ifges(6,1) * t359 - t436) * t482 + t57 * t481 + (Ifges(6,5) * t359 - Ifges(6,6) * t130) * t474 - t23 * t104 + t24 * t105 - g(1) * t403 - g(2) * t404 - t348 * t453 + (t130 * t24 + t23 * t359) * mrSges(6,3) + t337 + (-Ifges(6,2) * t130 + t125 + t58) * t484;];
tau = t1;
