% Calculate time derivative of joint inertia matrix for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:56
% EndTime: 2019-03-10 03:00:25
% DurationCPUTime: 13.63s
% Computational Cost: add. (21798->1013), mult. (65608->1396), div. (0->0), fcn. (65625->12), ass. (0->373)
t348 = sin(qJ(5));
t488 = (Ifges(7,4) + Ifges(6,5)) * t348;
t345 = sin(pkin(6));
t487 = 0.2e1 * t345;
t355 = cos(qJ(2));
t347 = cos(pkin(6));
t457 = pkin(1) * t347;
t335 = t355 * t457;
t323 = qJD(2) * t335;
t351 = sin(qJ(2));
t346 = cos(pkin(7));
t379 = t345 * (-pkin(10) * t346 - pkin(9));
t362 = t351 * t379;
t214 = pkin(2) * t347 + t335 + t362;
t441 = t214 * t346;
t485 = qJD(2) * t362 + qJD(3) * t441 + t323;
t349 = sin(qJ(4));
t353 = cos(qJ(4));
t420 = qJD(4) * t353;
t401 = t348 * t420;
t352 = cos(qJ(5));
t417 = qJD(5) * t352;
t357 = t349 * t417 + t401;
t354 = cos(qJ(3));
t430 = t351 * t354;
t350 = sin(qJ(3));
t431 = t350 * t355;
t360 = t346 * t431 + t430;
t344 = sin(pkin(7));
t439 = t344 * t350;
t207 = t345 * t360 + t347 * t439;
t437 = t345 * t355;
t269 = -t344 * t437 + t347 * t346;
t165 = t207 * t353 + t269 * t349;
t438 = t344 * t354;
t429 = t354 * t355;
t432 = t350 * t351;
t478 = t346 * t429 - t432;
t206 = -t345 * t478 - t347 * t438;
t122 = t165 * t352 + t206 * t348;
t423 = qJD(3) * t350;
t403 = t344 * t423;
t158 = t347 * t403 + (t360 * qJD(3) + (t346 * t430 + t431) * qJD(2)) * t345;
t422 = qJD(3) * t354;
t402 = t344 * t422;
t159 = t347 * t402 + (t478 * qJD(3) + (-t346 * t432 + t429) * qJD(2)) * t345;
t164 = t207 * t349 - t269 * t353;
t425 = qJD(2) * t345;
t405 = t351 * t425;
t382 = t344 * t405;
t97 = -qJD(4) * t164 + t159 * t353 + t349 * t382;
t42 = qJD(5) * t122 - t158 * t352 + t348 * t97;
t121 = t165 * t348 - t206 * t352;
t43 = -qJD(5) * t121 + t158 * t348 + t352 * t97;
t96 = qJD(4) * t165 + t159 * t349 - t353 * t382;
t7 = Ifges(6,5) * t43 - Ifges(6,6) * t42 + Ifges(6,3) * t96;
t8 = Ifges(7,4) * t43 + Ifges(7,2) * t96 + Ifges(7,6) * t42;
t484 = t7 + t8;
t270 = -t353 * t346 + t349 * t439;
t217 = -qJD(4) * t270 + t353 * t402;
t271 = t346 * t349 + t353 * t439;
t219 = t271 * t348 + t352 * t438;
t142 = -qJD(5) * t219 + t217 * t352 + t348 * t403;
t407 = t348 * t438;
t143 = -qJD(5) * t407 + t217 * t348 + t271 * t417 - t352 * t403;
t218 = qJD(4) * t271 + t349 * t402;
t71 = Ifges(6,5) * t142 - Ifges(6,6) * t143 + Ifges(6,3) * t218;
t72 = Ifges(7,4) * t142 + Ifges(7,2) * t218 + Ifges(7,6) * t143;
t482 = t71 + t72;
t456 = pkin(10) * t344;
t239 = (-pkin(2) * t355 - t351 * t456 - pkin(1)) * t345;
t153 = -t214 * t344 + t346 * t239;
t101 = pkin(3) * t206 - pkin(11) * t207 + t153;
t334 = t351 * t457;
t276 = pkin(9) * t437 + t334;
t203 = (t344 * t347 + t346 * t437) * pkin(10) + t276;
t436 = t346 * t350;
t118 = t354 * t203 + t214 * t436 + t239 * t439;
t107 = pkin(11) * t269 + t118;
t421 = qJD(4) * t349;
t216 = (t355 * t379 - t334) * qJD(2);
t243 = (pkin(2) * t351 - t355 * t456) * t425;
t64 = -t203 * t423 + t216 * t436 + t239 * t402 + t243 * t439 + t354 * t485;
t62 = pkin(11) * t382 + t64;
t160 = -t216 * t344 + t346 * t243;
t77 = pkin(3) * t158 - pkin(11) * t159 + t160;
t16 = -t101 * t421 - t107 * t420 - t349 * t62 + t353 * t77;
t14 = -pkin(4) * t158 - t16;
t18 = mrSges(6,1) * t42 + mrSges(6,2) * t43;
t481 = -m(6) * t14 - t18;
t275 = pkin(2) * t436 + pkin(10) * t438;
t250 = pkin(11) * t346 + t275;
t251 = (-pkin(3) * t354 - pkin(11) * t350 - pkin(2)) * t344;
t424 = qJD(3) * t344;
t263 = (pkin(3) * t350 - pkin(11) * t354) * t424;
t329 = pkin(10) * t439;
t435 = t346 * t354;
t273 = pkin(2) * t435 - t329;
t264 = t273 * qJD(3);
t124 = -t250 * t420 - t251 * t421 + t263 * t353 - t349 * t264;
t120 = -pkin(4) * t403 - t124;
t79 = mrSges(6,1) * t143 + mrSges(6,2) * t142;
t480 = -m(6) * t120 - t79;
t249 = t329 + (-pkin(2) * t354 - pkin(3)) * t346;
t166 = pkin(4) * t270 - pkin(12) * t271 + t249;
t180 = t353 * t250 + t349 * t251;
t168 = -pkin(12) * t438 + t180;
t479 = t348 * t166 + t352 * t168;
t442 = t349 * mrSges(5,2);
t477 = -m(5) * pkin(3) - t353 * mrSges(5,1) + t442;
t300 = (pkin(4) * t349 - pkin(12) * t353) * qJD(4);
t306 = -pkin(4) * t353 - pkin(12) * t349 - pkin(3);
t416 = qJD(5) * t353;
t419 = qJD(5) * t348;
t178 = pkin(11) * (t348 * t421 - t352 * t416) + t300 * t352 - t306 * t419;
t117 = -t350 * t203 + t354 * (t239 * t344 + t441);
t123 = -t250 * t421 + t251 * t420 + t349 * t263 + t353 * t264;
t119 = pkin(12) * t403 + t123;
t265 = t275 * qJD(3);
t136 = pkin(4) * t218 - pkin(12) * t217 + t265;
t36 = -qJD(5) * t479 - t119 * t348 + t136 * t352;
t15 = t101 * t420 - t107 * t421 + t349 * t77 + t353 * t62;
t13 = pkin(12) * t158 + t15;
t384 = -t203 * t422 - t239 * t403 - t350 * t485;
t63 = -t216 * t435 + (-pkin(3) * t405 - t243 * t354) * t344 - t384;
t24 = pkin(4) * t96 - pkin(12) * t97 + t63;
t50 = t349 * t101 + t353 * t107;
t46 = pkin(12) * t206 + t50;
t106 = -pkin(3) * t269 - t117;
t59 = pkin(4) * t164 - pkin(12) * t165 + t106;
t454 = t348 * t59 + t352 * t46;
t4 = -qJD(5) * t454 - t13 * t348 + t24 * t352;
t369 = pkin(5) * t352 + qJ(6) * t348;
t415 = qJD(6) * t352;
t476 = qJD(5) * t369 - t415;
t475 = 2 * m(4);
t474 = 0.2e1 * m(5);
t473 = 0.2e1 * m(6);
t472 = 2 * m(7);
t471 = 0.2e1 * pkin(11);
t470 = -2 * mrSges(3,3);
t469 = -2 * mrSges(4,3);
t34 = Ifges(5,1) * t97 - Ifges(5,4) * t96 + Ifges(5,5) * t158;
t467 = t34 / 0.2e1;
t86 = Ifges(5,1) * t165 - Ifges(5,4) * t164 + Ifges(5,5) * t206;
t466 = t86 / 0.2e1;
t148 = Ifges(5,1) * t217 - Ifges(5,4) * t218 + Ifges(5,5) * t403;
t464 = t148 / 0.2e1;
t183 = Ifges(5,1) * t271 - Ifges(5,4) * t270 - Ifges(5,5) * t438;
t463 = t183 / 0.2e1;
t451 = Ifges(5,4) * t349;
t295 = (Ifges(5,1) * t353 - t451) * qJD(4);
t462 = t295 / 0.2e1;
t461 = Ifges(5,5) * t349 / 0.2e1 + Ifges(5,6) * t353 / 0.2e1;
t450 = Ifges(5,4) * t353;
t318 = Ifges(5,1) * t349 + t450;
t460 = t318 / 0.2e1;
t459 = t346 / 0.2e1;
t455 = pkin(11) * t353;
t453 = Ifges(4,4) * t350;
t452 = Ifges(4,4) * t354;
t449 = Ifges(6,4) * t348;
t448 = Ifges(6,4) * t352;
t447 = Ifges(7,5) * t348;
t446 = Ifges(7,5) * t352;
t445 = Ifges(6,6) * t352;
t266 = -pkin(9) * t405 + t323;
t444 = t266 * mrSges(3,2);
t267 = t276 * qJD(2);
t443 = t267 * mrSges(3,1);
t440 = t306 * t352;
t434 = t348 * t349;
t433 = t349 * t352;
t370 = Ifges(7,3) * t348 + t446;
t252 = -Ifges(7,6) * t353 + t349 * t370;
t371 = -Ifges(6,2) * t348 + t448;
t255 = -Ifges(6,6) * t353 + t349 * t371;
t428 = t252 - t255;
t372 = Ifges(7,1) * t352 + t447;
t256 = -Ifges(7,4) * t353 + t349 * t372;
t373 = Ifges(6,1) * t352 - t449;
t257 = -Ifges(6,5) * t353 + t349 * t373;
t427 = t256 + t257;
t400 = t352 * t420;
t426 = Ifges(6,5) * t400 + Ifges(6,3) * t421;
t245 = t348 * t306 + t352 * t455;
t290 = Ifges(7,4) * t417 + Ifges(7,6) * t419;
t418 = qJD(5) * t349;
t6 = Ifges(7,5) * t43 + Ifges(7,6) * t96 + Ifges(7,3) * t42;
t9 = Ifges(6,4) * t43 - Ifges(6,2) * t42 + Ifges(6,6) * t96;
t414 = t6 / 0.2e1 - t9 / 0.2e1;
t32 = Ifges(5,5) * t97 - Ifges(5,6) * t96 + Ifges(5,3) * t158;
t10 = Ifges(7,1) * t43 + Ifges(7,4) * t96 + Ifges(7,5) * t42;
t11 = Ifges(6,1) * t43 - Ifges(6,4) * t42 + Ifges(6,5) * t96;
t412 = t10 / 0.2e1 + t11 / 0.2e1;
t51 = Ifges(7,5) * t122 + Ifges(7,6) * t164 + Ifges(7,3) * t121;
t54 = Ifges(6,4) * t122 - Ifges(6,2) * t121 + Ifges(6,6) * t164;
t411 = t51 / 0.2e1 - t54 / 0.2e1;
t55 = Ifges(7,1) * t122 + Ifges(7,4) * t164 + Ifges(7,5) * t121;
t56 = Ifges(6,1) * t122 - Ifges(6,4) * t121 + Ifges(6,5) * t164;
t410 = t55 / 0.2e1 + t56 / 0.2e1;
t70 = Ifges(7,5) * t142 + Ifges(7,6) * t218 + Ifges(7,3) * t143;
t73 = Ifges(6,4) * t142 - Ifges(6,2) * t143 + Ifges(6,6) * t218;
t409 = t70 / 0.2e1 - t73 / 0.2e1;
t74 = Ifges(7,1) * t142 + Ifges(7,4) * t218 + Ifges(7,5) * t143;
t75 = Ifges(6,1) * t142 - Ifges(6,4) * t143 + Ifges(6,5) * t218;
t408 = t74 / 0.2e1 + t75 / 0.2e1;
t93 = Ifges(4,5) * t159 - Ifges(4,6) * t158 + Ifges(4,3) * t382;
t146 = Ifges(5,5) * t217 - Ifges(5,6) * t218 + Ifges(5,3) * t403;
t399 = t348 * t418;
t220 = t271 * t352 - t407;
t130 = Ifges(7,5) * t220 + Ifges(7,6) * t270 + Ifges(7,3) * t219;
t133 = Ifges(6,4) * t220 - Ifges(6,2) * t219 + Ifges(6,6) * t270;
t397 = t130 / 0.2e1 - t133 / 0.2e1;
t134 = Ifges(7,1) * t220 + Ifges(7,4) * t270 + Ifges(7,5) * t219;
t135 = Ifges(6,1) * t220 - Ifges(6,4) * t219 + Ifges(6,5) * t270;
t396 = t135 / 0.2e1 + t134 / 0.2e1;
t310 = -Ifges(7,3) * t352 + t447;
t186 = -t310 * t418 + (Ifges(7,6) * t349 + t353 * t370) * qJD(4);
t314 = Ifges(6,2) * t352 + t449;
t189 = -t314 * t418 + (Ifges(6,6) * t349 + t353 * t371) * qJD(4);
t395 = t186 / 0.2e1 - t189 / 0.2e1;
t316 = Ifges(7,1) * t348 - t446;
t190 = -t316 * t418 + (Ifges(7,4) * t349 + t353 * t372) * qJD(4);
t317 = Ifges(6,1) * t348 + t448;
t191 = -t317 * t418 + (Ifges(6,5) * t349 + t353 * t373) * qJD(4);
t394 = t190 / 0.2e1 + t191 / 0.2e1;
t393 = t252 / 0.2e1 - t255 / 0.2e1;
t392 = t256 / 0.2e1 + t257 / 0.2e1;
t287 = t370 * qJD(5);
t291 = t371 * qJD(5);
t391 = t287 / 0.2e1 - t291 / 0.2e1;
t288 = Ifges(6,5) * t417 - Ifges(6,6) * t419;
t390 = t288 / 0.2e1 + t290 / 0.2e1;
t293 = t372 * qJD(5);
t294 = t373 * qJD(5);
t389 = t293 / 0.2e1 + t294 / 0.2e1;
t388 = -Ifges(7,6) * t352 / 0.2e1 + t445 / 0.2e1 + t488 / 0.2e1;
t387 = -t314 / 0.2e1 + t310 / 0.2e1;
t386 = t316 / 0.2e1 + t317 / 0.2e1;
t29 = -t96 * mrSges(7,1) + t43 * mrSges(7,2);
t113 = -t218 * mrSges(7,1) + t142 * mrSges(7,2);
t49 = t101 * t353 - t349 * t107;
t179 = -t349 * t250 + t251 * t353;
t33 = Ifges(5,4) * t97 - Ifges(5,2) * t96 + Ifges(5,6) * t158;
t385 = t7 / 0.2e1 + t8 / 0.2e1 - t33 / 0.2e1;
t383 = Ifges(7,4) * t400 + Ifges(7,2) * t421 + Ifges(7,6) * t357;
t52 = Ifges(6,5) * t122 - Ifges(6,6) * t121 + Ifges(6,3) * t164;
t53 = Ifges(7,4) * t122 + Ifges(7,2) * t164 + Ifges(7,6) * t121;
t85 = Ifges(5,4) * t165 - Ifges(5,2) * t164 + Ifges(5,6) * t206;
t381 = t52 / 0.2e1 + t53 / 0.2e1 - t85 / 0.2e1;
t147 = Ifges(5,4) * t217 - Ifges(5,2) * t218 + Ifges(5,6) * t403;
t380 = t71 / 0.2e1 + t72 / 0.2e1 - t147 / 0.2e1;
t167 = pkin(4) * t438 - t179;
t131 = Ifges(6,5) * t220 - Ifges(6,6) * t219 + Ifges(6,3) * t270;
t132 = Ifges(7,4) * t220 + Ifges(7,2) * t270 + Ifges(7,6) * t219;
t182 = Ifges(5,4) * t271 - Ifges(5,2) * t270 - Ifges(5,6) * t438;
t378 = -t182 / 0.2e1 + t131 / 0.2e1 + t132 / 0.2e1;
t187 = -Ifges(6,5) * t399 - Ifges(6,6) * t357 + t426;
t188 = -Ifges(7,4) * t399 + t383;
t292 = (-Ifges(5,2) * t349 + t450) * qJD(4);
t377 = -t292 / 0.2e1 + t187 / 0.2e1 + t188 / 0.2e1;
t253 = -Ifges(6,3) * t353 + (Ifges(6,5) * t352 - Ifges(6,6) * t348) * t349;
t254 = -Ifges(7,2) * t353 + (Ifges(7,4) * t352 + Ifges(7,6) * t348) * t349;
t315 = Ifges(5,2) * t353 + t451;
t376 = -t315 / 0.2e1 + t253 / 0.2e1 + t254 / 0.2e1;
t308 = -t352 * mrSges(6,1) + t348 * mrSges(6,2);
t375 = mrSges(6,1) * t348 + mrSges(6,2) * t352;
t307 = -t352 * mrSges(7,1) - t348 * mrSges(7,3);
t374 = mrSges(7,1) * t348 - mrSges(7,3) * t352;
t368 = pkin(5) * t348 - qJ(6) * t352;
t21 = -t348 * t46 + t352 * t59;
t108 = t166 * t352 - t168 * t348;
t361 = pkin(11) + t368;
t45 = -pkin(4) * t206 - t49;
t3 = t352 * t13 + t348 * t24 + t59 * t417 - t419 * t46;
t35 = t352 * t119 + t348 * t136 + t166 * t417 - t168 * t419;
t358 = -t399 + t400;
t233 = mrSges(7,2) * t400 + (-mrSges(7,1) * qJD(4) - mrSges(7,2) * t419) * t349;
t177 = t348 * t300 + t306 * t417 + (-t348 * t416 - t352 * t421) * pkin(11);
t342 = Ifges(5,5) * t420;
t322 = Ifges(3,5) * t355 * t425;
t321 = Ifges(4,5) * t402;
t303 = -pkin(4) - t369;
t299 = -mrSges(7,2) * t434 - mrSges(7,3) * t353;
t298 = mrSges(7,1) * t353 + mrSges(7,2) * t433;
t297 = -mrSges(6,1) * t353 - mrSges(6,3) * t433;
t296 = mrSges(6,2) * t353 - mrSges(6,3) * t434;
t289 = -Ifges(5,6) * t421 + t342;
t286 = (mrSges(5,1) * t349 + mrSges(5,2) * t353) * qJD(4);
t285 = t375 * qJD(5);
t284 = t374 * qJD(5);
t283 = -mrSges(4,2) * t346 + mrSges(4,3) * t438;
t282 = mrSges(4,1) * t346 - mrSges(4,3) * t439;
t278 = t375 * t349;
t277 = t374 * t349;
t274 = -pkin(9) * t345 * t351 + t335;
t268 = qJD(5) * t368 - qJD(6) * t348;
t262 = (Ifges(4,1) * t354 - t453) * t424;
t261 = (-Ifges(4,2) * t350 + t452) * t424;
t260 = -Ifges(4,6) * t403 + t321;
t259 = (mrSges(4,1) * t350 + mrSges(4,2) * t354) * t424;
t258 = t361 * t349;
t247 = Ifges(4,5) * t346 + (Ifges(4,1) * t350 + t452) * t344;
t246 = Ifges(4,6) * t346 + (Ifges(4,2) * t354 + t453) * t344;
t244 = -t348 * t455 + t440;
t235 = -mrSges(7,2) * t357 + mrSges(7,3) * t421;
t234 = -mrSges(6,2) * t421 - mrSges(6,3) * t357;
t232 = mrSges(6,1) * t421 - mrSges(6,3) * t358;
t227 = -t440 + (pkin(11) * t348 + pkin(5)) * t353;
t226 = -qJ(6) * t353 + t245;
t225 = -mrSges(5,1) * t438 - mrSges(5,3) * t271;
t224 = mrSges(5,2) * t438 - mrSges(5,3) * t270;
t201 = mrSges(6,1) * t357 + mrSges(6,2) * t358;
t200 = mrSges(7,1) * t357 - mrSges(7,3) * t358;
t195 = mrSges(5,1) * t270 + mrSges(5,2) * t271;
t185 = -mrSges(5,2) * t403 - mrSges(5,3) * t218;
t184 = mrSges(5,1) * t403 - mrSges(5,3) * t217;
t181 = Ifges(5,5) * t271 - Ifges(5,6) * t270 - Ifges(5,3) * t438;
t176 = t349 * t476 + t361 * t420;
t175 = -mrSges(7,1) * t270 + mrSges(7,2) * t220;
t174 = mrSges(6,1) * t270 - mrSges(6,3) * t220;
t173 = -mrSges(6,2) * t270 - mrSges(6,3) * t219;
t172 = -mrSges(7,2) * t219 + mrSges(7,3) * t270;
t171 = mrSges(4,1) * t269 - mrSges(4,3) * t207;
t170 = -mrSges(4,2) * t269 - mrSges(4,3) * t206;
t169 = -pkin(5) * t421 - t178;
t162 = qJ(6) * t421 - qJD(6) * t353 + t177;
t151 = mrSges(6,1) * t219 + mrSges(6,2) * t220;
t150 = mrSges(7,1) * t219 - mrSges(7,3) * t220;
t149 = mrSges(5,1) * t218 + mrSges(5,2) * t217;
t145 = mrSges(4,1) * t382 - mrSges(4,3) * t159;
t144 = -mrSges(4,2) * t382 - mrSges(4,3) * t158;
t128 = Ifges(4,1) * t207 - Ifges(4,4) * t206 + Ifges(4,5) * t269;
t127 = Ifges(4,4) * t207 - Ifges(4,2) * t206 + Ifges(4,6) * t269;
t126 = mrSges(5,1) * t206 - mrSges(5,3) * t165;
t125 = -mrSges(5,2) * t206 - mrSges(5,3) * t164;
t115 = -mrSges(7,2) * t143 + mrSges(7,3) * t218;
t114 = -mrSges(6,2) * t218 - mrSges(6,3) * t143;
t112 = mrSges(6,1) * t218 - mrSges(6,3) * t142;
t111 = mrSges(5,1) * t164 + mrSges(5,2) * t165;
t110 = pkin(5) * t219 - qJ(6) * t220 + t167;
t105 = mrSges(4,1) * t158 + mrSges(4,2) * t159;
t95 = Ifges(4,1) * t159 - Ifges(4,4) * t158 + Ifges(4,5) * t382;
t94 = Ifges(4,4) * t159 - Ifges(4,2) * t158 + Ifges(4,6) * t382;
t92 = -pkin(5) * t270 - t108;
t91 = qJ(6) * t270 + t479;
t84 = Ifges(5,5) * t165 - Ifges(5,6) * t164 + Ifges(5,3) * t206;
t83 = -mrSges(7,1) * t164 + mrSges(7,2) * t122;
t82 = mrSges(6,1) * t164 - mrSges(6,3) * t122;
t81 = -mrSges(6,2) * t164 - mrSges(6,3) * t121;
t80 = -mrSges(7,2) * t121 + mrSges(7,3) * t164;
t78 = mrSges(7,1) * t143 - mrSges(7,3) * t142;
t69 = mrSges(6,1) * t121 + mrSges(6,2) * t122;
t68 = mrSges(7,1) * t121 - mrSges(7,3) * t122;
t67 = mrSges(5,1) * t158 - mrSges(5,3) * t97;
t66 = -mrSges(5,2) * t158 - mrSges(5,3) * t96;
t65 = (t216 * t346 + t243 * t344) * t354 + t384;
t48 = mrSges(5,1) * t96 + mrSges(5,2) * t97;
t47 = pkin(5) * t143 - qJ(6) * t142 - qJD(6) * t220 + t120;
t31 = -pkin(5) * t218 - t36;
t30 = qJ(6) * t218 + qJD(6) * t270 + t35;
t28 = mrSges(6,1) * t96 - mrSges(6,3) * t43;
t27 = -mrSges(6,2) * t96 - mrSges(6,3) * t42;
t26 = -mrSges(7,2) * t42 + mrSges(7,3) * t96;
t25 = pkin(5) * t121 - qJ(6) * t122 + t45;
t20 = -pkin(5) * t164 - t21;
t19 = qJ(6) * t164 + t454;
t17 = mrSges(7,1) * t42 - mrSges(7,3) * t43;
t5 = pkin(5) * t42 - qJ(6) * t43 - qJD(6) * t122 + t14;
t2 = -pkin(5) * t96 - t4;
t1 = qJ(6) * t96 + qJD(6) * t164 + t3;
t12 = [(t14 * t45 + t21 * t4 + t3 * t454) * t473 + 0.2e1 * t454 * t27 + 0.2e1 * m(3) * (t266 * t276 - t267 * t274) + (t322 - 0.2e1 * t443 - 0.2e1 * t444) * t347 + (-t33 + t484) * t164 + (t1 * t19 + t2 * t20 + t25 * t5) * t472 + (t106 * t63 + t15 * t50 + t16 * t49) * t474 + (t117 * t65 + t118 * t64 + t153 * t160) * t475 + (t84 - t127) * t158 + (t6 - t9) * t121 + (t11 + t10) * t122 + t269 * t93 + t206 * t32 - t206 * t94 + 0.2e1 * t160 * (mrSges(4,1) * t206 + mrSges(4,2) * t207) + t207 * t95 + 0.2e1 * t64 * t170 + 0.2e1 * t65 * t171 + t165 * t34 + t159 * t128 + 0.2e1 * t153 * t105 + 0.2e1 * t118 * t144 + 0.2e1 * t117 * t145 + (t52 + t53 - t85) * t96 + 0.2e1 * t15 * t125 + 0.2e1 * t16 * t126 + 0.2e1 * t63 * t111 + 0.2e1 * t106 * t48 + t97 * t86 + 0.2e1 * t1 * t80 + 0.2e1 * t3 * t81 + 0.2e1 * t4 * t82 + 0.2e1 * t2 * t83 + 0.2e1 * t5 * t68 + 0.2e1 * t14 * t69 + 0.2e1 * t50 * t66 + 0.2e1 * t49 * t67 + 0.2e1 * t45 * t18 + 0.2e1 * t21 * t28 + 0.2e1 * t20 * t29 + 0.2e1 * t25 * t17 + 0.2e1 * t19 * t26 + (t55 + t56) * t43 + (t51 - t54) * t42 + (0.2e1 * (t266 * t355 + t267 * t351) * mrSges(3,3) + ((t274 * t470 + Ifges(3,5) * t347 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t355) * t487) * t355 + (t344 * (Ifges(4,5) * t207 - Ifges(4,6) * t206 + Ifges(4,3) * t269) + t276 * t470 - 0.2e1 * Ifges(3,6) * t347 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t351 + (Ifges(3,1) - Ifges(3,2)) * t355) * t487) * t351) * qJD(2)) * t345; t454 * t114 + m(7) * (t1 * t91 + t110 * t5 + t19 * t30 + t2 * t92 + t20 * t31 + t25 * t47) + m(5) * (t106 * t265 + t123 * t50 + t124 * t49 + t15 * t180 + t16 * t179 + t249 * t63) + (-t261 / 0.2e1 + t146 / 0.2e1) * t206 + ((t160 * mrSges(4,2) + t95 / 0.2e1) * t350 + (-t160 * mrSges(4,1) + t94 / 0.2e1 - t32 / 0.2e1) * t354 + (Ifges(4,3) * t459 + (Ifges(4,5) * t350 + Ifges(4,6) * t354) * t344 / 0.2e1) * t405 + (-m(4) * t160 - t105) * pkin(2) + ((-t117 * mrSges(4,3) + t128 / 0.2e1) * t354 + (-t118 * mrSges(4,3) - t127 / 0.2e1 + t84 / 0.2e1) * t350) * qJD(3)) * t344 + t414 * t219 + t409 * t121 + t410 * t142 + t411 * t143 + t412 * t220 + t408 * t122 - Ifges(3,6) * t405 + t397 * t42 + t396 * t43 + (-t246 / 0.2e1 + t181 / 0.2e1) * t158 + t385 * t270 + t322 + t381 * t218 + t380 * t164 + t378 * t96 + t97 * t463 + t165 * t464 + t217 * t466 + t271 * t467 + t93 * t459 + (t111 - t171) * t265 + m(4) * (-t117 * t265 + t118 * t264 + t273 * t65 + t275 * t64) - t443 - t444 + t65 * t282 + t64 * t283 + t273 * t145 + t275 * t144 + t269 * t260 / 0.2e1 + t153 * t259 + t207 * t262 / 0.2e1 + t264 * t170 + t159 * t247 / 0.2e1 + t249 * t48 + t15 * t224 + t16 * t225 + t63 * t195 + t179 * t67 + t180 * t66 + t49 * t184 + t50 * t185 + t1 * t172 + t3 * t173 + t4 * t174 + t2 * t175 + t167 * t18 + t5 * t150 + t14 * t151 + t106 * t149 + t123 * t125 + t124 * t126 + t120 * t69 + t110 * t17 + t21 * t112 + t20 * t113 + t19 * t115 + m(6) * (t108 * t4 + t120 * t45 + t14 * t167 + t21 * t36 + t3 * t479 + t35 * t454) + t479 * t27 + t108 * t28 + t91 * t26 + t92 * t29 + t45 * t79 + t30 * t80 + t35 * t81 + t36 * t82 + t31 * t83 + t25 * t78 + t47 * t68; (-t147 + t482) * t270 + (-0.2e1 * pkin(2) * t259 + t350 * t262 + (-t146 + t261) * t354 + ((t273 * t469 + t247) * t354 + (t275 * t469 + t181 - t246) * t350) * qJD(3)) * t344 + (t110 * t47 + t30 * t91 + t31 * t92) * t472 + (t123 * t180 + t124 * t179 + t249 * t265) * t474 + (t264 * t275 - t265 * t273) * t475 + 0.2e1 * (-t282 + t195) * t265 + t346 * t260 + 0.2e1 * t264 * t283 + t271 * t148 + 0.2e1 * t249 * t149 + 0.2e1 * t123 * t224 + 0.2e1 * t124 * t225 + t217 * t183 + 0.2e1 * t179 * t184 + 0.2e1 * t180 * t185 + 0.2e1 * t30 * t172 + 0.2e1 * t35 * t173 + 0.2e1 * t36 * t174 + 0.2e1 * t31 * t175 + 0.2e1 * t167 * t79 + 0.2e1 * t47 * t150 + 0.2e1 * t120 * t151 + 0.2e1 * t110 * t78 + 0.2e1 * t108 * t112 + 0.2e1 * t92 * t113 + 0.2e1 * t91 * t115 + (t108 * t36 + t120 * t167 + t35 * t479) * t473 + 0.2e1 * t479 * t114 + (t74 + t75) * t220 + (t70 - t73) * t219 + (t131 + t132 - t182) * t218 + (-t133 + t130) * t143 + (t134 + t135) * t142; t454 * t234 + m(6) * (t177 * t454 + t178 * t21 + t244 * t4 + t245 * t3) + m(7) * (t1 * t226 + t162 * t19 + t169 * t20 + t176 * t25 + t2 * t227 + t258 * t5) + (t467 - t16 * mrSges(5,3) + t412 * t352 + t414 * t348 + (-t348 * t410 + t352 * t411) * qJD(5) + (-t50 * mrSges(5,3) + t381) * qJD(4) + (-qJD(4) * t125 - t67 + m(5) * (-qJD(4) * t50 - t16) - t481) * pkin(11)) * t349 + t477 * t63 + (t15 * mrSges(5,3) + (m(5) * t15 + t66) * pkin(11) + (-t49 * mrSges(5,3) + t466 + t410 * t352 + t411 * t348 + (-m(5) * t49 + m(6) * t45 - t126 + t69) * pkin(11)) * qJD(4) - t385) * t353 + t392 * t43 + t393 * t42 + t394 * t122 + t395 * t121 + t376 * t96 + t377 * t164 + t97 * t460 + t158 * t461 + t165 * t462 + t3 * t296 + t4 * t297 + t2 * t298 + t1 * t299 + t5 * t277 + t14 * t278 + t106 * t286 + t206 * t289 / 0.2e1 + t258 * t17 + t21 * t232 + t20 * t233 + t19 * t235 + t244 * t28 + t245 * t27 + t226 * t26 + t227 * t29 + t25 * t200 + t45 * t201 + t169 * t83 + t176 * t68 + t177 * t81 + t178 * t82 + t162 * t80 - t64 * mrSges(4,2) + t65 * mrSges(4,1) - pkin(3) * t48 + t93; (t464 - t124 * mrSges(5,3) + t408 * t352 + t409 * t348 + (-t348 * t396 + t352 * t397) * qJD(5) + (-t180 * mrSges(5,3) + t378) * qJD(4) + (-qJD(4) * t224 - t184 + m(5) * (-qJD(4) * t180 - t124) - t480) * pkin(11)) * t349 + (-mrSges(4,1) + t477) * t265 + (t123 * mrSges(5,3) + (m(5) * t123 + t185) * pkin(11) + (-t179 * mrSges(5,3) + t463 + t396 * t352 + t397 * t348 + (-m(5) * t179 + m(6) * t167 + t151 - t225) * pkin(11)) * qJD(4) - t380) * t353 + (-t354 * t289 / 0.2e1 + (-Ifges(4,6) + t461) * t423) * t344 + m(7) * (t110 * t176 + t162 * t91 + t169 * t92 + t226 * t30 + t227 * t31 + t258 * t47) + t392 * t142 + t393 * t143 + t394 * t220 + t395 * t219 + t321 + t376 * t218 + t377 * t270 + t271 * t462 + t217 * t460 + t35 * t296 + t36 * t297 + t31 * t298 + t30 * t299 + t47 * t277 + t120 * t278 + t249 * t286 + t258 * t78 - t264 * mrSges(4,2) + t108 * t232 + t92 * t233 + t91 * t235 + t244 * t112 + t245 * t114 + t226 * t115 + t227 * t113 + t110 * t200 + t167 * t201 + t162 * t172 + t169 * t175 + t176 * t150 + t177 * t173 + t178 * t174 - pkin(3) * t149 + m(6) * (t108 * t178 + t177 * t479 + t244 * t36 + t245 * t35) + t479 * t234; -0.2e1 * pkin(3) * t286 + 0.2e1 * t162 * t299 + 0.2e1 * t169 * t298 + 0.2e1 * t176 * t277 + 0.2e1 * t177 * t296 + 0.2e1 * t178 * t297 + 0.2e1 * t258 * t200 + 0.2e1 * t226 * t235 + 0.2e1 * t227 * t233 + 0.2e1 * t244 * t232 + 0.2e1 * t245 * t234 + (t162 * t226 + t169 * t227 + t176 * t258) * t472 + (t177 * t245 + t178 * t244) * t473 + (-t187 - t188 + t292 + (t278 * t471 + t348 * t428 + t352 * t427 + t318) * qJD(4)) * t353 + (t201 * t471 + t295 + (t190 + t191) * t352 + (t186 - t189) * t348 + (-t348 * t427 + t352 * t428) * qJD(5) + (pkin(11) ^ 2 * t353 * t473 + t253 + t254 - t315) * qJD(4)) * t349; ((t20 * mrSges(7,2) - t21 * mrSges(6,3) + t410) * t352 + (-t19 * mrSges(7,2) - mrSges(6,3) * t454 + t411) * t348) * qJD(5) + (t2 * mrSges(7,2) - t4 * mrSges(6,3) + t412) * t348 + t387 * t42 + t388 * t96 + t389 * t122 + t390 * t164 + t391 * t121 + t386 * t43 + ((t26 + t27) * t352 + (-t28 + t29) * t348 + ((-t82 + t83) * t352 + (-t80 - t81) * t348) * qJD(5) + m(6) * (-t21 * t417 + t3 * t352 - t348 * t4 - t419 * t454) + m(7) * (t1 * t352 - t19 * t419 + t2 * t348 + t20 * t417)) * pkin(12) + t32 + t303 * t17 + t5 * t307 + t14 * t308 + t25 * t284 + t45 * t285 + t268 * t68 + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t414) * t352 + m(7) * (t25 * t268 + t303 * t5) - t15 * mrSges(5,2) + t16 * mrSges(5,1) + t481 * pkin(4); (t30 * mrSges(7,2) + t35 * mrSges(6,3) - t409) * t352 + (t31 * mrSges(7,2) - t36 * mrSges(6,3) + t408) * t348 + m(7) * (t110 * t268 + t303 * t47) + ((t92 * mrSges(7,2) - t108 * mrSges(6,3) + t396) * t352 + (-t91 * mrSges(7,2) - mrSges(6,3) * t479 + t397) * t348) * qJD(5) + t387 * t143 + t388 * t218 + t389 * t220 + t390 * t270 + t391 * t219 + t386 * t142 + ((t114 + t115) * t352 + (-t112 + t113) * t348 + ((-t174 + t175) * t352 + (-t172 - t173) * t348) * qJD(5) + m(6) * (-t108 * t417 - t348 * t36 + t35 * t352 - t419 * t479) + m(7) * (t30 * t352 + t31 * t348 + t417 * t92 - t419 * t91)) * pkin(12) + t146 + t303 * t78 + t47 * t307 + t120 * t308 + t110 * t284 + t167 * t285 + t268 * t150 - t123 * mrSges(5,2) + t124 * mrSges(5,1) + t480 * pkin(4); t342 + m(7) * (t176 * t303 + t258 * t268) + t303 * t200 + t176 * t307 + t268 * t277 + t258 * t284 - pkin(4) * t201 + t349 * pkin(11) * t285 - t390 * t353 + ((-Ifges(5,6) + t388) * t349 + (t442 + (-m(6) * pkin(4) - mrSges(5,1) + t308) * t353) * pkin(11)) * qJD(4) + (t162 * mrSges(7,2) + t177 * mrSges(6,3) + t389 * t349 + t386 * t420 + (t227 * mrSges(7,2) - t244 * mrSges(6,3) + t387 * t349 + t392) * qJD(5) + (t234 + t235 + (-t297 + t298) * qJD(5) + m(7) * (qJD(5) * t227 + t162) + m(6) * (-qJD(5) * t244 + t177)) * pkin(12) - t395) * t352 + (t169 * mrSges(7,2) - t178 * mrSges(6,3) + t391 * t349 + t387 * t420 + (-t226 * mrSges(7,2) - t245 * mrSges(6,3) - t386 * t349 + t393) * qJD(5) + (-t232 + t233 + (-t296 - t299) * qJD(5) + m(7) * (-qJD(5) * t226 + t169) + m(6) * (-qJD(5) * t245 - t178)) * pkin(12) + t394) * t348; -0.2e1 * pkin(4) * t285 + 0.2e1 * t284 * t303 + (-t287 + t291) * t352 + (t293 + t294) * t348 + 0.2e1 * (m(7) * t303 + t307) * t268 + ((t316 + t317) * t352 + (t310 - t314) * t348) * qJD(5); -pkin(5) * t29 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t19) + qJD(6) * t80 + qJ(6) * t26 + t1 * mrSges(7,3) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) + t484; t36 * mrSges(6,1) - pkin(5) * t113 + m(7) * (-pkin(5) * t31 + qJ(6) * t30 + qJD(6) * t91) + qJD(6) * t172 + qJ(6) * t115 + t30 * mrSges(7,3) - t35 * mrSges(6,2) - t31 * mrSges(7,1) + t482; -Ifges(6,6) * t401 - pkin(5) * t233 + m(7) * (-pkin(5) * t169 + qJ(6) * t162 + qJD(6) * t226) + qJD(6) * t299 + qJ(6) * t235 + t162 * mrSges(7,3) - t169 * mrSges(7,1) - t177 * mrSges(6,2) + t178 * mrSges(6,1) + (-t445 - t488) * t418 + t383 + t426; -t476 * mrSges(7,2) + (m(7) * t415 + (-m(7) * t369 + t307 + t308) * qJD(5)) * pkin(12) + t288 + t290; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t29; m(7) * t31 + t113; m(7) * t169 + t233; (m(7) * pkin(12) + mrSges(7,2)) * t417; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
