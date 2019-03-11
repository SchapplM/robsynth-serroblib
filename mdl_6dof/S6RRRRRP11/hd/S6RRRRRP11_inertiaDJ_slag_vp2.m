% Calculate time derivative of joint inertia matrix for
% S6RRRRRP11
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:46
% EndTime: 2019-03-10 02:33:20
% DurationCPUTime: 15.20s
% Computational Cost: add. (21695->1011), mult. (65376->1396), div. (0->0), fcn. (65546->12), ass. (0->379)
t359 = cos(qJ(5));
t489 = Ifges(6,6) + Ifges(7,6);
t493 = t489 * t359;
t355 = sin(qJ(5));
t492 = (Ifges(6,5) + Ifges(7,5)) * t355;
t352 = sin(pkin(6));
t491 = 0.2e1 * t352;
t362 = cos(qJ(2));
t354 = cos(pkin(6));
t456 = pkin(1) * t354;
t340 = t362 * t456;
t327 = qJD(2) * t340;
t358 = sin(qJ(2));
t353 = cos(pkin(7));
t377 = t352 * (-pkin(10) * t353 - pkin(9));
t367 = t358 * t377;
t216 = pkin(2) * t354 + t340 + t367;
t440 = t216 * t353;
t488 = qJD(2) * t367 + qJD(3) * t440 + t327;
t360 = cos(qJ(4));
t416 = qJD(4) * t360;
t396 = t359 * t416;
t356 = sin(qJ(4));
t414 = qJD(5) * t356;
t398 = t355 * t414;
t364 = t396 - t398;
t399 = t355 * t416;
t413 = qJD(5) * t359;
t363 = t356 * t413 + t399;
t361 = cos(qJ(3));
t430 = t358 * t361;
t357 = sin(qJ(3));
t431 = t357 * t362;
t365 = t353 * t431 + t430;
t351 = sin(pkin(7));
t439 = t351 * t357;
t209 = t352 * t365 + t354 * t439;
t437 = t352 * t362;
t271 = -t351 * t437 + t354 * t353;
t169 = t209 * t360 + t271 * t356;
t438 = t351 * t361;
t428 = t361 * t362;
t432 = t357 * t358;
t479 = t353 * t428 - t432;
t208 = -t352 * t479 - t354 * t438;
t123 = t169 * t359 + t208 * t355;
t419 = qJD(3) * t357;
t401 = t351 * t419;
t162 = t354 * t401 + (t365 * qJD(3) + (t353 * t430 + t431) * qJD(2)) * t352;
t418 = qJD(3) * t361;
t400 = t351 * t418;
t163 = t354 * t400 + (t479 * qJD(3) + (-t353 * t432 + t428) * qJD(2)) * t352;
t168 = t209 * t356 - t271 * t360;
t421 = qJD(2) * t352;
t403 = t358 * t421;
t380 = t351 * t403;
t99 = -qJD(4) * t168 + t163 * t360 + t356 * t380;
t43 = -qJD(5) * t123 + t162 * t359 - t355 * t99;
t122 = -t169 * t355 + t208 * t359;
t44 = qJD(5) * t122 + t162 * t355 + t359 * t99;
t98 = qJD(4) * t169 + t163 * t356 - t360 * t380;
t6 = Ifges(7,5) * t44 + Ifges(7,6) * t43 + Ifges(7,3) * t98;
t7 = Ifges(6,5) * t44 + Ifges(6,6) * t43 + Ifges(6,3) * t98;
t487 = t6 + t7;
t272 = -t360 * t353 + t356 * t439;
t219 = -qJD(4) * t272 + t360 * t400;
t273 = t353 * t356 + t360 * t439;
t221 = -t273 * t355 - t359 * t438;
t146 = qJD(5) * t221 + t219 * t359 + t355 * t401;
t366 = -t273 * t359 + t355 * t438;
t147 = qJD(5) * t366 - t219 * t355 + t359 * t401;
t220 = qJD(4) * t273 + t356 * t400;
t71 = Ifges(7,5) * t146 + Ifges(7,6) * t147 + Ifges(7,3) * t220;
t72 = Ifges(6,5) * t146 + Ifges(6,6) * t147 + Ifges(6,3) * t220;
t483 = t71 + t72;
t482 = -m(6) * pkin(12) - mrSges(6,3);
t436 = t353 * t357;
t278 = pkin(2) * t436 + pkin(10) * t438;
t252 = pkin(11) * t353 + t278;
t253 = (-pkin(3) * t361 - pkin(11) * t357 - pkin(2)) * t351;
t420 = qJD(3) * t351;
t264 = (pkin(3) * t357 - pkin(11) * t361) * t420;
t334 = pkin(10) * t439;
t435 = t353 * t361;
t276 = pkin(2) * t435 - t334;
t265 = t276 * qJD(3);
t417 = qJD(4) * t356;
t125 = -t252 * t416 - t253 * t417 + t264 * t360 - t356 * t265;
t121 = -pkin(4) * t401 - t125;
t82 = -mrSges(6,1) * t147 + mrSges(6,2) * t146;
t481 = -m(6) * t121 - t82;
t455 = pkin(10) * t351;
t241 = (-pkin(2) * t362 - t358 * t455 - pkin(1)) * t352;
t157 = -t216 * t351 + t353 * t241;
t103 = pkin(3) * t208 - pkin(11) * t209 + t157;
t339 = t358 * t456;
t279 = pkin(9) * t437 + t339;
t204 = (t351 * t354 + t353 * t437) * pkin(10) + t279;
t119 = t361 * t204 + t216 * t436 + t241 * t439;
t109 = pkin(11) * t271 + t119;
t218 = (t362 * t377 - t339) * qJD(2);
t245 = (pkin(2) * t358 - t362 * t455) * t421;
t65 = -t204 * t419 + t218 * t436 + t241 * t400 + t245 * t439 + t361 * t488;
t63 = pkin(11) * t380 + t65;
t164 = -t218 * t351 + t353 * t245;
t78 = pkin(3) * t162 - pkin(11) * t163 + t164;
t17 = -t103 * t417 - t109 * t416 - t356 * t63 + t360 * t78;
t14 = -pkin(4) * t162 - t17;
t20 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t480 = -m(6) * t14 - t20;
t251 = t334 + (-pkin(2) * t361 - pkin(3)) * t353;
t170 = pkin(4) * t272 - pkin(12) * t273 + t251;
t182 = t360 * t252 + t356 * t253;
t172 = -pkin(12) * t438 + t182;
t111 = t355 * t170 + t359 * t172;
t308 = -pkin(4) * t360 - pkin(12) * t356 - pkin(3);
t429 = t359 * t360;
t341 = pkin(11) * t429;
t247 = t355 * t308 + t341;
t441 = t356 * mrSges(5,2);
t478 = -m(5) * pkin(3) - t360 * mrSges(5,1) + t441;
t118 = -t357 * t204 + t361 * (t241 * t351 + t440);
t477 = 2 * m(4);
t476 = 0.2e1 * m(5);
t475 = 0.2e1 * m(6);
t474 = 2 * m(7);
t473 = 0.2e1 * pkin(11);
t472 = -2 * mrSges(3,3);
t471 = -2 * mrSges(4,3);
t470 = -2 * mrSges(7,3);
t467 = m(7) * pkin(5);
t35 = Ifges(5,1) * t99 - Ifges(5,4) * t98 + Ifges(5,5) * t162;
t466 = t35 / 0.2e1;
t90 = Ifges(5,1) * t169 - Ifges(5,4) * t168 + Ifges(5,5) * t208;
t465 = t90 / 0.2e1;
t153 = Ifges(5,1) * t219 - Ifges(5,4) * t220 + Ifges(5,5) * t401;
t463 = t153 / 0.2e1;
t185 = Ifges(5,1) * t273 - Ifges(5,4) * t272 - Ifges(5,5) * t438;
t462 = t185 / 0.2e1;
t449 = Ifges(5,4) * t356;
t299 = (Ifges(5,1) * t360 - t449) * qJD(4);
t461 = t299 / 0.2e1;
t460 = Ifges(5,5) * t356 / 0.2e1 + Ifges(5,6) * t360 / 0.2e1;
t448 = Ifges(5,4) * t360;
t322 = Ifges(5,1) * t356 + t448;
t459 = t322 / 0.2e1;
t458 = t353 / 0.2e1;
t454 = pkin(11) * t355;
t452 = -qJ(6) - pkin(12);
t50 = t356 * t103 + t360 * t109;
t47 = pkin(12) * t208 + t50;
t108 = -pkin(3) * t271 - t118;
t60 = pkin(4) * t168 - pkin(12) * t169 + t108;
t22 = t355 * t60 + t359 * t47;
t451 = Ifges(4,4) * t357;
t450 = Ifges(4,4) * t361;
t447 = Ifges(6,4) * t355;
t446 = Ifges(6,4) * t359;
t445 = Ifges(7,4) * t355;
t444 = Ifges(7,4) * t359;
t267 = -pkin(9) * t403 + t327;
t443 = t267 * mrSges(3,2);
t268 = t279 * qJD(2);
t442 = t268 * mrSges(3,1);
t434 = t355 * t356;
t433 = t356 * t359;
t369 = -Ifges(7,2) * t355 + t444;
t256 = -Ifges(7,6) * t360 + t356 * t369;
t370 = -Ifges(6,2) * t355 + t446;
t257 = -Ifges(6,6) * t360 + t356 * t370;
t427 = -t256 - t257;
t371 = Ifges(7,1) * t359 - t445;
t258 = -Ifges(7,5) * t360 + t356 * t371;
t372 = Ifges(6,1) * t359 - t447;
t259 = -Ifges(6,5) * t360 + t356 * t372;
t426 = t258 + t259;
t304 = (pkin(4) * t356 - pkin(12) * t360) * qJD(4);
t425 = t355 * t304 + t308 * t413;
t424 = t359 * t304 + t417 * t454;
t423 = Ifges(7,5) * t396 + Ifges(7,3) * t417;
t422 = Ifges(6,5) * t396 + Ifges(6,3) * t417;
t415 = qJD(5) * t355;
t288 = mrSges(7,1) * t415 + mrSges(7,2) * t413;
t412 = qJD(6) * t359;
t8 = Ifges(7,4) * t44 + Ifges(7,2) * t43 + Ifges(7,6) * t98;
t9 = Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t98;
t411 = t8 / 0.2e1 + t9 / 0.2e1;
t33 = Ifges(5,5) * t99 - Ifges(5,6) * t98 + Ifges(5,3) * t162;
t10 = Ifges(7,1) * t44 + Ifges(7,4) * t43 + Ifges(7,5) * t98;
t11 = Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t98;
t409 = t10 / 0.2e1 + t11 / 0.2e1;
t53 = Ifges(7,4) * t123 + Ifges(7,2) * t122 + Ifges(7,6) * t168;
t54 = Ifges(6,4) * t123 + Ifges(6,2) * t122 + Ifges(6,6) * t168;
t408 = t53 / 0.2e1 + t54 / 0.2e1;
t55 = Ifges(7,1) * t123 + Ifges(7,4) * t122 + Ifges(7,5) * t168;
t56 = Ifges(6,1) * t123 + Ifges(6,4) * t122 + Ifges(6,5) * t168;
t407 = t55 / 0.2e1 + t56 / 0.2e1;
t73 = Ifges(7,4) * t146 + Ifges(7,2) * t147 + Ifges(7,6) * t220;
t74 = Ifges(6,4) * t146 + Ifges(6,2) * t147 + Ifges(6,6) * t220;
t406 = t73 / 0.2e1 + t74 / 0.2e1;
t75 = Ifges(7,1) * t146 + Ifges(7,4) * t147 + Ifges(7,5) * t220;
t76 = Ifges(6,1) * t146 + Ifges(6,4) * t147 + Ifges(6,5) * t220;
t405 = t76 / 0.2e1 + t75 / 0.2e1;
t95 = Ifges(4,5) * t163 - Ifges(4,6) * t162 + Ifges(4,3) * t380;
t151 = Ifges(5,5) * t219 - Ifges(5,6) * t220 + Ifges(5,3) * t401;
t135 = -Ifges(7,4) * t366 + Ifges(7,2) * t221 + Ifges(7,6) * t272;
t136 = -Ifges(6,4) * t366 + Ifges(6,2) * t221 + Ifges(6,6) * t272;
t395 = t135 / 0.2e1 + t136 / 0.2e1;
t137 = -Ifges(7,1) * t366 + Ifges(7,4) * t221 + Ifges(7,5) * t272;
t138 = -Ifges(6,1) * t366 + Ifges(6,4) * t221 + Ifges(6,5) * t272;
t394 = t137 / 0.2e1 + t138 / 0.2e1;
t317 = Ifges(7,2) * t359 + t445;
t190 = -t317 * t414 + (Ifges(7,6) * t356 + t360 * t369) * qJD(4);
t318 = Ifges(6,2) * t359 + t447;
t191 = -t318 * t414 + (Ifges(6,6) * t356 + t360 * t370) * qJD(4);
t393 = t190 / 0.2e1 + t191 / 0.2e1;
t320 = Ifges(7,1) * t355 + t444;
t192 = -t320 * t414 + (Ifges(7,5) * t356 + t360 * t371) * qJD(4);
t321 = Ifges(6,1) * t355 + t446;
t193 = -t321 * t414 + (Ifges(6,5) * t356 + t360 * t372) * qJD(4);
t392 = t193 / 0.2e1 + t192 / 0.2e1;
t391 = t256 / 0.2e1 + t257 / 0.2e1;
t390 = t258 / 0.2e1 + t259 / 0.2e1;
t348 = Ifges(7,5) * t413;
t349 = Ifges(6,5) * t413;
t389 = t348 / 0.2e1 + t349 / 0.2e1 - t489 * t415 / 0.2e1;
t294 = t369 * qJD(5);
t295 = t370 * qJD(5);
t388 = -t295 / 0.2e1 - t294 / 0.2e1;
t297 = t371 * qJD(5);
t298 = t372 * qJD(5);
t387 = t297 / 0.2e1 + t298 / 0.2e1;
t386 = t492 / 0.2e1 + t493 / 0.2e1;
t385 = -t318 / 0.2e1 - t317 / 0.2e1;
t384 = t320 / 0.2e1 + t321 / 0.2e1;
t19 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t21 = -t355 * t47 + t359 * t60;
t383 = qJD(5) * t452;
t81 = -t147 * mrSges(7,1) + t146 * mrSges(7,2);
t49 = t103 * t360 - t356 * t109;
t110 = t359 * t170 - t172 * t355;
t181 = -t356 * t252 + t253 * t360;
t34 = Ifges(5,4) * t99 - Ifges(5,2) * t98 + Ifges(5,6) * t162;
t382 = t6 / 0.2e1 + t7 / 0.2e1 - t34 / 0.2e1;
t381 = -t204 * t418 - t241 * t401 - t357 * t488;
t51 = Ifges(7,5) * t123 + Ifges(7,6) * t122 + Ifges(7,3) * t168;
t52 = Ifges(6,5) * t123 + Ifges(6,6) * t122 + Ifges(6,3) * t168;
t89 = Ifges(5,4) * t169 - Ifges(5,2) * t168 + Ifges(5,6) * t208;
t379 = t51 / 0.2e1 + t52 / 0.2e1 - t89 / 0.2e1;
t152 = Ifges(5,4) * t219 - Ifges(5,2) * t220 + Ifges(5,6) * t401;
t378 = t71 / 0.2e1 + t72 / 0.2e1 - t152 / 0.2e1;
t171 = pkin(4) * t438 - t181;
t133 = -Ifges(7,5) * t366 + Ifges(7,6) * t221 + Ifges(7,3) * t272;
t134 = -Ifges(6,5) * t366 + Ifges(6,6) * t221 + Ifges(6,3) * t272;
t184 = Ifges(5,4) * t273 - Ifges(5,2) * t272 - Ifges(5,6) * t438;
t376 = t133 / 0.2e1 + t134 / 0.2e1 - t184 / 0.2e1;
t188 = -Ifges(7,5) * t398 - Ifges(7,6) * t363 + t423;
t189 = -Ifges(6,5) * t398 - Ifges(6,6) * t363 + t422;
t296 = (-Ifges(5,2) * t356 + t448) * qJD(4);
t375 = -t296 / 0.2e1 + t188 / 0.2e1 + t189 / 0.2e1;
t254 = -Ifges(7,3) * t360 + (Ifges(7,5) * t359 - Ifges(7,6) * t355) * t356;
t255 = -Ifges(6,3) * t360 + (Ifges(6,5) * t359 - Ifges(6,6) * t355) * t356;
t319 = Ifges(5,2) * t360 + t449;
t374 = -t319 / 0.2e1 + t254 / 0.2e1 + t255 / 0.2e1;
t373 = mrSges(6,1) * t355 + mrSges(6,2) * t359;
t46 = -pkin(4) * t208 - t49;
t16 = t103 * t416 - t109 * t417 + t356 * t78 + t360 * t63;
t13 = pkin(12) * t162 + t16;
t64 = -t218 * t435 + (-pkin(3) * t403 - t245 * t361) * t351 - t381;
t25 = pkin(4) * t98 - pkin(12) * t99 + t64;
t3 = t359 * t13 + t355 * t25 + t60 * t413 - t415 * t47;
t124 = -t252 * t417 + t253 * t416 + t356 * t264 + t360 * t265;
t120 = pkin(12) * t401 + t124;
t266 = t278 * qJD(3);
t139 = pkin(4) * t220 - pkin(12) * t219 + t266;
t36 = t359 * t120 + t355 * t139 + t170 * t413 - t172 * t415;
t202 = mrSges(7,1) * t363 + mrSges(7,2) * t364;
t4 = -qJD(5) * t22 - t13 * t355 + t359 * t25;
t37 = -qJD(5) * t111 - t120 * t355 + t359 * t139;
t350 = Ifges(5,5) * t416;
t343 = -pkin(5) * t359 - pkin(4);
t326 = Ifges(3,5) * t362 * t421;
t325 = Ifges(4,5) * t400;
t313 = t452 * t359;
t311 = -mrSges(6,1) * t359 + mrSges(6,2) * t355;
t310 = -mrSges(7,1) * t359 + mrSges(7,2) * t355;
t309 = t452 * t355;
t306 = (pkin(5) * t355 + pkin(11)) * t356;
t303 = -mrSges(6,1) * t360 - mrSges(6,3) * t433;
t302 = -mrSges(7,1) * t360 - mrSges(7,3) * t433;
t301 = mrSges(6,2) * t360 - mrSges(6,3) * t434;
t300 = mrSges(7,2) * t360 - mrSges(7,3) * t434;
t293 = -Ifges(5,6) * t417 + t350;
t290 = (mrSges(5,1) * t356 + mrSges(5,2) * t360) * qJD(4);
t289 = t373 * qJD(5);
t287 = -mrSges(4,2) * t353 + mrSges(4,3) * t438;
t286 = mrSges(4,1) * t353 - mrSges(4,3) * t439;
t285 = t359 * t308;
t281 = t373 * t356;
t280 = (mrSges(7,1) * t355 + mrSges(7,2) * t359) * t356;
t277 = -pkin(9) * t352 * t358 + t340;
t270 = -qJD(6) * t355 + t359 * t383;
t269 = t355 * t383 + t412;
t263 = (Ifges(4,1) * t361 - t451) * t420;
t262 = (-Ifges(4,2) * t357 + t450) * t420;
t261 = -Ifges(4,6) * t401 + t325;
t260 = (mrSges(4,1) * t357 + mrSges(4,2) * t361) * t420;
t249 = Ifges(4,5) * t353 + (Ifges(4,1) * t357 + t450) * t351;
t248 = Ifges(4,6) * t353 + (Ifges(4,2) * t361 + t451) * t351;
t246 = -t360 * t454 + t285;
t240 = pkin(5) * t363 + pkin(11) * t416;
t236 = -mrSges(6,2) * t417 - mrSges(6,3) * t363;
t235 = -mrSges(7,2) * t417 - mrSges(7,3) * t363;
t234 = mrSges(6,1) * t417 - mrSges(6,3) * t364;
t233 = mrSges(7,1) * t417 - mrSges(7,3) * t364;
t228 = -mrSges(5,1) * t438 - mrSges(5,3) * t273;
t227 = mrSges(5,2) * t438 - mrSges(5,3) * t272;
t224 = -qJ(6) * t434 + t247;
t207 = -qJ(6) * t433 + t285 + (-pkin(5) - t454) * t360;
t203 = mrSges(6,1) * t363 + mrSges(6,2) * t364;
t197 = mrSges(5,1) * t272 + mrSges(5,2) * t273;
t187 = -mrSges(5,2) * t401 - mrSges(5,3) * t220;
t186 = mrSges(5,1) * t401 - mrSges(5,3) * t219;
t183 = Ifges(5,5) * t273 - Ifges(5,6) * t272 - Ifges(5,3) * t438;
t180 = -qJD(5) * t247 + t424;
t179 = (-t359 * t417 - t360 * t415) * pkin(11) + t425;
t178 = mrSges(6,1) * t272 + mrSges(6,3) * t366;
t177 = mrSges(7,1) * t272 + mrSges(7,3) * t366;
t176 = -mrSges(6,2) * t272 + mrSges(6,3) * t221;
t175 = -mrSges(7,2) * t272 + mrSges(7,3) * t221;
t174 = mrSges(4,1) * t271 - mrSges(4,3) * t209;
t173 = -mrSges(4,2) * t271 - mrSges(4,3) * t208;
t156 = -mrSges(6,1) * t221 - mrSges(6,2) * t366;
t155 = -mrSges(7,1) * t221 - mrSges(7,2) * t366;
t154 = mrSges(5,1) * t220 + mrSges(5,2) * t219;
t150 = (-pkin(11) * qJD(4) - qJ(6) * qJD(5)) * t433 + (-qJD(6) * t356 + (-pkin(11) * qJD(5) - qJ(6) * qJD(4)) * t360) * t355 + t425;
t149 = mrSges(4,1) * t380 - mrSges(4,3) * t163;
t148 = -mrSges(4,2) * t380 - mrSges(4,3) * t162;
t140 = -t356 * t412 + (pkin(5) * t356 - qJ(6) * t429) * qJD(4) + (-t341 + (qJ(6) * t356 - t308) * t355) * qJD(5) + t424;
t130 = Ifges(4,1) * t209 - Ifges(4,4) * t208 + Ifges(4,5) * t271;
t129 = Ifges(4,4) * t209 - Ifges(4,2) * t208 + Ifges(4,6) * t271;
t128 = -pkin(5) * t221 + t171;
t127 = mrSges(5,1) * t208 - mrSges(5,3) * t169;
t126 = -mrSges(5,2) * t208 - mrSges(5,3) * t168;
t116 = -mrSges(6,2) * t220 + mrSges(6,3) * t147;
t115 = -mrSges(7,2) * t220 + mrSges(7,3) * t147;
t114 = mrSges(6,1) * t220 - mrSges(6,3) * t146;
t113 = mrSges(7,1) * t220 - mrSges(7,3) * t146;
t112 = mrSges(5,1) * t168 + mrSges(5,2) * t169;
t107 = mrSges(4,1) * t162 + mrSges(4,2) * t163;
t97 = Ifges(4,1) * t163 - Ifges(4,4) * t162 + Ifges(4,5) * t380;
t96 = Ifges(4,4) * t163 - Ifges(4,2) * t162 + Ifges(4,6) * t380;
t88 = Ifges(5,5) * t169 - Ifges(5,6) * t168 + Ifges(5,3) * t208;
t87 = qJ(6) * t221 + t111;
t86 = mrSges(6,1) * t168 - mrSges(6,3) * t123;
t85 = mrSges(7,1) * t168 - mrSges(7,3) * t123;
t84 = -mrSges(6,2) * t168 + mrSges(6,3) * t122;
t83 = -mrSges(7,2) * t168 + mrSges(7,3) * t122;
t80 = pkin(5) * t272 + qJ(6) * t366 + t110;
t79 = -pkin(5) * t147 + t121;
t70 = -mrSges(6,1) * t122 + mrSges(6,2) * t123;
t69 = -mrSges(7,1) * t122 + mrSges(7,2) * t123;
t68 = mrSges(5,1) * t162 - mrSges(5,3) * t99;
t67 = -mrSges(5,2) * t162 - mrSges(5,3) * t98;
t66 = (t218 * t353 + t245 * t351) * t361 + t381;
t48 = mrSges(5,1) * t98 + mrSges(5,2) * t99;
t32 = -pkin(5) * t122 + t46;
t31 = qJ(6) * t147 + qJD(6) * t221 + t36;
t30 = mrSges(6,1) * t98 - mrSges(6,3) * t44;
t29 = mrSges(7,1) * t98 - mrSges(7,3) * t44;
t28 = -mrSges(6,2) * t98 + mrSges(6,3) * t43;
t27 = -mrSges(7,2) * t98 + mrSges(7,3) * t43;
t26 = pkin(5) * t220 - qJ(6) * t146 + qJD(6) * t366 + t37;
t18 = qJ(6) * t122 + t22;
t15 = pkin(5) * t168 - qJ(6) * t123 + t21;
t5 = -pkin(5) * t43 + t14;
t2 = qJ(6) * t43 + qJD(6) * t122 + t3;
t1 = pkin(5) * t98 - qJ(6) * t44 - qJD(6) * t123 + t4;
t12 = [(0.2e1 * (t267 * t362 + t268 * t358) * mrSges(3,3) + ((t277 * t472 + Ifges(3,5) * t354 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t362) * t491) * t362 + (t351 * (Ifges(4,5) * t209 - Ifges(4,6) * t208 + Ifges(4,3) * t271) + t279 * t472 - 0.2e1 * Ifges(3,6) * t354 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t358 + (Ifges(3,1) - Ifges(3,2)) * t362) * t491) * t358) * qJD(2)) * t352 + (t10 + t11) * t123 + (t8 + t9) * t122 + (t88 - t129) * t162 + (-t34 + t487) * t168 + (t14 * t46 + t21 * t4 + t22 * t3) * t475 + (t108 * t64 + t16 * t50 + t17 * t49) * t476 + (t118 * t66 + t119 * t65 + t157 * t164) * t477 + 0.2e1 * m(3) * (t267 * t279 - t268 * t277) + (t326 - 0.2e1 * t442 - 0.2e1 * t443) * t354 + (t1 * t15 + t18 * t2 + t32 * t5) * t474 + t271 * t95 + 0.2e1 * t164 * (mrSges(4,1) * t208 + mrSges(4,2) * t209) + t209 * t97 + t208 * t33 - t208 * t96 + (t53 + t54) * t43 + (t55 + t56) * t44 + 0.2e1 * t18 * t27 + 0.2e1 * t22 * t28 + 0.2e1 * t15 * t29 + 0.2e1 * t21 * t30 + 0.2e1 * t32 * t19 + 0.2e1 * t46 * t20 + 0.2e1 * t50 * t67 + 0.2e1 * t49 * t68 + 0.2e1 * t5 * t69 + 0.2e1 * t14 * t70 + 0.2e1 * t2 * t83 + 0.2e1 * t3 * t84 + 0.2e1 * t1 * t85 + 0.2e1 * t4 * t86 + t99 * t90 + 0.2e1 * t108 * t48 + 0.2e1 * t64 * t112 + 0.2e1 * t16 * t126 + 0.2e1 * t17 * t127 + (t51 + t52 - t89) * t98 + 0.2e1 * t119 * t148 + 0.2e1 * t118 * t149 + 0.2e1 * t157 * t107 + t163 * t130 + t169 * t35 + 0.2e1 * t65 * t173 + 0.2e1 * t66 * t174; -t443 + m(6) * (t110 * t4 + t111 * t3 + t121 * t46 + t14 * t171 + t21 * t37 + t22 * t36) + m(7) * (t1 * t80 + t128 * t5 + t15 * t26 + t18 * t31 + t2 * t87 + t32 * t79) + m(5) * (t108 * t266 + t124 * t50 + t125 * t49 + t16 * t182 + t17 * t181 + t251 * t64) - t409 * t366 + ((t164 * mrSges(4,2) + t97 / 0.2e1) * t357 + (-t164 * mrSges(4,1) + t96 / 0.2e1 - t33 / 0.2e1) * t361 + (Ifges(4,3) * t458 + (Ifges(4,5) * t357 + Ifges(4,6) * t361) * t351 / 0.2e1) * t403 + (-m(4) * t164 - t107) * pkin(2) + ((-t118 * mrSges(4,3) + t130 / 0.2e1) * t361 + (-t119 * mrSges(4,3) + t88 / 0.2e1 - t129 / 0.2e1) * t357) * qJD(3)) * t351 + (-t262 / 0.2e1 + t151 / 0.2e1) * t208 + t411 * t221 + t405 * t123 + t406 * t122 + t407 * t146 + t408 * t147 - Ifges(3,6) * t403 + t395 * t43 + t394 * t44 + t382 * t272 + t95 * t458 + t99 * t462 + t169 * t463 + t219 * t465 + t273 * t466 + (-t248 / 0.2e1 + t183 / 0.2e1) * t162 + t379 * t220 + t378 * t168 + t376 * t98 + (-t174 + t112) * t266 + m(4) * (-t118 * t266 + t119 * t265 + t276 * t66 + t278 * t65) - t442 + t66 * t286 + t65 * t287 + t276 * t149 + t278 * t148 + t271 * t261 / 0.2e1 + t157 * t260 + t209 * t263 / 0.2e1 + t265 * t173 + t163 * t249 / 0.2e1 + t251 * t48 + t16 * t227 + t17 * t228 + t49 * t186 + t50 * t187 + t64 * t197 + t326 + t79 * t69 + t80 * t29 + t32 * t81 + t46 * t82 + t31 * t83 + t36 * t84 + t26 * t85 + t37 * t86 + t87 * t27 + t110 * t30 + t111 * t28 + t15 * t113 + t21 * t114 + t18 * t115 + t22 * t116 + t121 * t70 + t124 * t126 + t125 * t127 + t128 * t19 + t108 * t154 + t5 * t155 + t14 * t156 + t171 * t20 + t2 * t175 + t3 * t176 + t1 * t177 + t4 * t178 + t181 * t68 + t182 * t67; (-t152 + t483) * t272 + (t135 + t136) * t147 + (t137 + t138) * t146 - (t75 + t76) * t366 + (-0.2e1 * pkin(2) * t260 + t263 * t357 + (-t151 + t262) * t361 + ((t276 * t471 + t249) * t361 + (t278 * t471 + t183 - t248) * t357) * qJD(3)) * t351 + (t124 * t182 + t125 * t181 + t251 * t266) * t476 + (t265 * t278 - t266 * t276) * t477 + (t128 * t79 + t26 * t80 + t31 * t87) * t474 + (t110 * t37 + t111 * t36 + t121 * t171) * t475 + 0.2e1 * (-t286 + t197) * t266 + t353 * t261 + 0.2e1 * t265 * t287 + t273 * t153 + 0.2e1 * t251 * t154 + 0.2e1 * t124 * t227 + 0.2e1 * t125 * t228 + t219 * t185 + 0.2e1 * t181 * t186 + 0.2e1 * t182 * t187 + (t133 + t134 - t184) * t220 + (t73 + t74) * t221 + 0.2e1 * t80 * t113 + 0.2e1 * t110 * t114 + 0.2e1 * t87 * t115 + 0.2e1 * t111 * t116 + 0.2e1 * t128 * t81 + 0.2e1 * t79 * t155 + 0.2e1 * t121 * t156 + 0.2e1 * t171 * t82 + 0.2e1 * t31 * t175 + 0.2e1 * t36 * t176 + 0.2e1 * t26 * t177 + 0.2e1 * t37 * t178; t95 + m(7) * (t1 * t207 + t140 * t15 + t150 * t18 + t2 * t224 + t240 * t32 + t306 * t5) + (t16 * mrSges(5,3) + (m(5) * t16 + t67) * pkin(11) + (-t49 * mrSges(5,3) + t465 + t407 * t359 - t408 * t355 + (-m(5) * t49 + m(6) * t46 - t127 + t70) * pkin(11)) * qJD(4) - t382) * t360 + t478 * t64 + t390 * t44 + t391 * t43 + t392 * t123 + t393 * t122 + t99 * t459 + t162 * t460 + t169 * t461 + t374 * t98 + t375 * t168 + (t466 - t17 * mrSges(5,3) + t409 * t359 - t411 * t355 + (-t355 * t407 - t359 * t408) * qJD(5) + (-t50 * mrSges(5,3) + t379) * qJD(4) + (-qJD(4) * t126 - t68 + m(5) * (-qJD(4) * t50 - t17) - t480) * pkin(11)) * t356 + t2 * t300 + t3 * t301 + t1 * t302 + t4 * t303 + t306 * t19 + t108 * t290 + t208 * t293 / 0.2e1 + t14 * t281 + t5 * t280 + t22 * t236 + t240 * t69 + t246 * t30 + t247 * t28 + t224 * t27 + t15 * t233 + t21 * t234 + t18 * t235 + t207 * t29 + m(6) * (t179 * t22 + t180 * t21 + t246 * t4 + t247 * t3) + t32 * t202 + t46 * t203 - pkin(3) * t48 - t65 * mrSges(4,2) + t66 * mrSges(4,1) + t140 * t85 + t150 * t83 + t179 * t84 + t180 * t86; m(7) * (t128 * t240 + t140 * t80 + t150 * t87 + t207 * t26 + t224 * t31 + t306 * t79) - t392 * t366 + (t124 * mrSges(5,3) + (m(5) * t124 + t187) * pkin(11) + (-t181 * mrSges(5,3) + t462 + t394 * t359 - t395 * t355 + (-m(5) * t181 + m(6) * t171 + t156 - t228) * pkin(11)) * qJD(4) - t378) * t360 + (-t361 * t293 / 0.2e1 + (t460 - Ifges(4,6)) * t419) * t351 + (-mrSges(4,1) + t478) * t266 + t390 * t146 + t391 * t147 + t393 * t221 + t219 * t459 + t273 * t461 + t374 * t220 + t375 * t272 + (t463 - t125 * mrSges(5,3) + t405 * t359 - t406 * t355 + (-t355 * t394 - t359 * t395) * qJD(5) + (-t182 * mrSges(5,3) + t376) * qJD(4) + (-qJD(4) * t227 - t186 + m(5) * (-qJD(4) * t182 - t125) - t481) * pkin(11)) * t356 + m(6) * (t110 * t180 + t111 * t179 + t246 * t37 + t247 * t36) + t31 * t300 + t36 * t301 + t26 * t302 + t37 * t303 + t306 * t81 + t251 * t290 + t121 * t281 + t79 * t280 - t265 * mrSges(4,2) + t240 * t155 + t246 * t114 + t247 * t116 + t224 * t115 + t80 * t233 + t110 * t234 + t87 * t235 + t111 * t236 + t207 * t113 + t128 * t202 + t171 * t203 + t325 - pkin(3) * t154 + t150 * t175 + t140 * t177 + t179 * t176 + t180 * t178; -0.2e1 * pkin(3) * t290 + 0.2e1 * t140 * t302 + 0.2e1 * t150 * t300 + 0.2e1 * t179 * t301 + 0.2e1 * t180 * t303 + 0.2e1 * t306 * t202 + 0.2e1 * t207 * t233 + 0.2e1 * t224 * t235 + 0.2e1 * t246 * t234 + 0.2e1 * t247 * t236 + 0.2e1 * t240 * t280 + (t140 * t207 + t150 * t224 + t240 * t306) * t474 + (t179 * t247 + t180 * t246) * t475 + (-t188 - t189 + t296 + (t281 * t473 + t355 * t427 + t359 * t426 + t322) * qJD(4)) * t360 + (t203 * t473 + t299 + (t192 + t193) * t359 + (-t190 - t191) * t355 + (-t355 * t426 + t359 * t427) * qJD(5) + (pkin(11) ^ 2 * t360 * t475 + t254 + t255 - t319) * qJD(4)) * t356; m(7) * (t1 * t309 + t15 * t270 + t18 * t269 - t2 * t313 + t343 * t5) + t33 + (t2 * mrSges(7,3) + t3 * mrSges(6,3) + (-t21 * mrSges(6,3) - t15 * mrSges(7,3) + t407) * qJD(5) + (-qJD(5) * t86 + t28 + m(6) * (-qJD(5) * t21 + t3)) * pkin(12) + t411) * t359 + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) + (-m(6) * t4 - t30) * pkin(12) + (-t18 * mrSges(7,3) + pkin(5) * t69 - pkin(12) * t84 + t22 * t482 + t32 * t467 - t408) * qJD(5) + t409) * t355 - t385 * t43 + t386 * t98 + t387 * t123 - t388 * t122 + t389 * t168 + t384 * t44 + t343 * t19 + t309 * t29 + t5 * t310 + t14 * t311 - t313 * t27 + t32 * t288 + t46 * t289 + t269 * t83 + t270 * t85 - t16 * mrSges(5,2) + t17 * mrSges(5,1) + t480 * pkin(4); t151 + (-t37 * mrSges(6,3) - t26 * mrSges(7,3) + (-m(6) * t37 - t114) * pkin(12) + (-t87 * mrSges(7,3) + pkin(5) * t155 - pkin(12) * t176 + t111 * t482 + t128 * t467 - t395) * qJD(5) + t405) * t355 - t385 * t147 + t386 * t220 - t387 * t366 - t388 * t221 + t389 * t272 + (t36 * mrSges(6,3) + t31 * mrSges(7,3) + (-t110 * mrSges(6,3) - t80 * mrSges(7,3) + t394) * qJD(5) + (-qJD(5) * t178 + t116 + m(6) * (-qJD(5) * t110 + t36)) * pkin(12) + t406) * t359 + t384 * t146 + m(7) * (t26 * t309 + t269 * t87 + t270 * t80 - t31 * t313 + t343 * t79) + t343 * t81 + t309 * t113 + t79 * t310 + t121 * t311 - t313 * t115 + t128 * t288 + t171 * t289 + t269 * t175 + t270 * t177 - t124 * mrSges(5,2) + t125 * mrSges(5,1) + t481 * pkin(4); t343 * t202 + t269 * t300 + t270 * t302 + t306 * t288 + t309 * t233 + t240 * t310 - t313 * t235 - pkin(4) * t203 + m(7) * (t140 * t309 - t150 * t313 + t207 * t270 + t224 * t269 + t240 * t343) + t356 * pkin(11) * t289 + t350 - t389 * t360 + ((-Ifges(5,6) + t386) * t356 + (t441 + (-m(6) * pkin(4) - mrSges(5,1) + t311) * t360) * pkin(11)) * qJD(4) + (-t180 * mrSges(6,3) - t140 * mrSges(7,3) + t388 * t356 + t385 * t416 + (-m(6) * t180 - t234) * pkin(12) + (-t224 * mrSges(7,3) + pkin(5) * t280 - pkin(12) * t301 + t247 * t482 + t306 * t467 - t356 * t384 - t391) * qJD(5) + t392) * t355 + (t179 * mrSges(6,3) + t150 * mrSges(7,3) + t387 * t356 + t384 * t416 + (m(6) * t179 + t236) * pkin(12) + (-t207 * mrSges(7,3) - t246 * mrSges(6,3) + t385 * t356 + (-m(6) * t246 - t303) * pkin(12) + t390) * qJD(5) + t393) * t359; 0.2e1 * t343 * t288 + (-t269 * t313 + t270 * t309) * t474 - 0.2e1 * pkin(4) * t289 + (t270 * t470 + t297 + t298 + (-t313 * t470 - t317 - t318 + 0.2e1 * (m(7) * t343 + t310) * pkin(5)) * qJD(5)) * t355 + (0.2e1 * t269 * mrSges(7,3) + t294 + t295 + (t309 * t470 + t320 + t321) * qJD(5)) * t359; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + (m(7) * t1 + t29) * pkin(5) + t487; mrSges(6,1) * t37 + mrSges(7,1) * t26 - mrSges(6,2) * t36 - mrSges(7,2) * t31 + (m(7) * t26 + t113) * pkin(5) + t483; mrSges(6,1) * t180 + mrSges(7,1) * t140 - mrSges(6,2) * t179 - mrSges(7,2) * t150 - t489 * t399 + (m(7) * t140 + t233) * pkin(5) + (-t492 - t493) * t414 + t422 + t423; -mrSges(7,2) * t269 + t348 + t349 + (mrSges(7,1) + t467) * t270 + ((-mrSges(6,1) * pkin(12) - (mrSges(7,3) * pkin(5))) * t359 + (mrSges(6,2) * pkin(12) - t489) * t355) * qJD(5); 0; m(7) * t5 + t19; m(7) * t79 + t81; m(7) * t240 + t202; t415 * t467 + t288; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
