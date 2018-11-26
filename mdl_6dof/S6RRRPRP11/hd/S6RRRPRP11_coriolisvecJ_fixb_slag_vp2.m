% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:48:53
% EndTime: 2018-11-23 17:49:17
% DurationCPUTime: 24.62s
% Computational Cost: add. (10854->797), mult. (28262->1057), div. (0->0), fcn. (20832->8), ass. (0->315)
t453 = Ifges(7,4) + Ifges(6,4);
t281 = cos(qJ(2));
t387 = cos(pkin(6));
t338 = t387 * qJD(1);
t331 = pkin(1) * t338;
t278 = sin(qJ(2));
t275 = sin(pkin(6));
t373 = qJD(1) * t275;
t351 = t278 * t373;
t219 = -pkin(8) * t351 + t281 * t331;
t298 = (pkin(2) * t278 - pkin(9) * t281) * t275;
t220 = qJD(1) * t298;
t277 = sin(qJ(3));
t280 = cos(qJ(3));
t146 = -t277 * t219 + t220 * t280;
t274 = t280 * pkin(4);
t370 = qJD(3) * t280;
t426 = pkin(4) + pkin(9);
t427 = pkin(3) + pkin(10);
t489 = -(t274 * t281 - t278 * t427) * t373 + t146 + t426 * t370;
t455 = Ifges(7,1) + Ifges(6,1);
t451 = Ifges(7,5) + Ifges(6,5);
t450 = Ifges(7,2) + Ifges(6,2);
t449 = Ifges(7,6) + Ifges(6,6);
t276 = sin(qJ(5));
t279 = cos(qJ(5));
t379 = t279 * t281;
t188 = (-t276 * t278 + t277 * t379) * t373;
t369 = qJD(5) * t276;
t488 = -t280 * t369 + t188;
t304 = pkin(10) * t277 - qJ(4) * t280;
t350 = t281 * t373;
t222 = pkin(8) * t350 + t278 * t331;
t334 = t277 * t350;
t355 = pkin(3) * t334 + t222;
t135 = t304 * t350 + t355;
t371 = qJD(3) * t277;
t336 = pkin(3) * t371 - qJD(4) * t277;
t200 = qJD(3) * t304 + t336;
t340 = -qJ(4) * t277 - pkin(2);
t234 = -t280 * t427 + t340;
t256 = t426 * t277;
t368 = qJD(5) * t279;
t480 = -t234 * t369 + t256 * t368 + (-t135 + t200) * t279 + t489 * t276;
t487 = t135 * t276 + t279 * t489;
t301 = t338 + qJD(2);
t294 = t280 * t301;
t201 = t277 * t351 - t294;
t255 = -qJD(3) + t350;
t152 = t201 * t279 + t255 * t276;
t486 = t453 * t152;
t182 = -pkin(2) * t301 - t219;
t202 = t277 * t301 + t280 * t351;
t287 = -t202 * qJ(4) + t182;
t110 = t201 * pkin(3) + t287;
t183 = pkin(9) * t301 + t222;
t214 = (-pkin(2) * t281 - pkin(9) * t278 - pkin(1)) * t275;
t193 = qJD(1) * t214;
t129 = t280 * t183 + t277 * t193;
t243 = t255 * qJ(4);
t116 = t243 - t129;
t122 = -t255 * Ifges(5,5) - t202 * Ifges(5,6) + t201 * Ifges(5,3);
t125 = t202 * Ifges(4,4) - t201 * Ifges(4,2) - t255 * Ifges(4,6);
t360 = Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1;
t365 = Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1;
t485 = t110 * mrSges(5,2) + t129 * mrSges(4,3) - t122 / 0.2e1 + t125 / 0.2e1 - t116 * mrSges(5,1) - t182 * mrSges(4,1) + t365 * t202 - t360 * t255;
t153 = t201 * t276 - t255 * t279;
t484 = t453 * t153;
t128 = t183 * t277 - t280 * t193;
t299 = pkin(4) * t202 + t128;
t69 = t255 * t427 + qJD(4) + t299;
t74 = t201 * t427 + t287;
t20 = -t276 * t74 + t279 * t69;
t21 = t276 * t69 + t279 * t74;
t303 = t20 * t276 - t21 * t279;
t318 = mrSges(7,1) * t279 - mrSges(7,2) * t276;
t320 = mrSges(6,1) * t279 - mrSges(6,2) * t276;
t12 = qJ(6) * t152 + t21;
t11 = -qJ(6) * t153 + t20;
t196 = qJD(5) + t202;
t9 = pkin(5) * t196 + t11;
t321 = t12 * t279 - t276 * t9;
t405 = -t279 / 0.2e1;
t407 = -t276 / 0.2e1;
t98 = -pkin(4) * t201 + t129;
t82 = -t243 + t98;
t44 = -pkin(5) * t152 + qJD(6) + t82;
t445 = t455 * t153 + t451 * t196 + t486;
t446 = t450 * t152 + t449 * t196 + t484;
t483 = t303 * mrSges(6,3) - t321 * mrSges(7,3) + t44 * t318 + t82 * t320 + t405 * t446 + t407 * t445;
t456 = Ifges(5,1) + Ifges(4,3);
t454 = Ifges(5,4) - Ifges(4,5);
t452 = Ifges(5,5) - Ifges(4,6);
t381 = t276 * t281;
t189 = (t277 * t381 + t278 * t279) * t373;
t333 = t280 * t350;
t337 = qJ(6) * t280 - t234;
t482 = qJ(6) * t189 + t337 * t368 + (-qJ(6) * t371 - qJD(5) * t256 + qJD(6) * t280 - t200) * t276 + t487 + (-t333 + t370) * pkin(5);
t367 = t279 * qJD(6);
t481 = -t280 * t367 + t480 + (t279 * t371 - t488) * qJ(6);
t169 = t279 * t234 + t276 * t256;
t479 = -qJD(5) * t169 - t200 * t276 + t487;
t147 = t280 * t219 + t277 * t220;
t134 = -qJ(4) * t351 - t147;
t121 = -pkin(4) * t334 - t134;
t353 = -pkin(5) * t279 - pkin(4);
t478 = -t121 + (-pkin(9) + t353) * t371 + t488 * pkin(5);
t477 = -t426 * t371 - t121;
t476 = t453 * t279;
t475 = t453 * t276;
t463 = -qJD(4) - t128;
t410 = t196 / 0.2e1;
t418 = t153 / 0.2e1;
t420 = t152 / 0.2e1;
t437 = t455 * t276 + t476;
t438 = t450 * t279 + t475;
t439 = t276 * t451 + t279 * t449;
t474 = t410 * t439 + t418 * t437 + t420 * t438 - t483 - t485;
t473 = t446 / 0.2e1;
t366 = qJD(1) * qJD(2);
t343 = t281 * t366;
t329 = t275 * t343;
t383 = t275 * t278;
t357 = t277 * t383;
t332 = qJD(3) * t357;
t161 = qJD(1) * t332 - qJD(3) * t294 - t280 * t329;
t162 = qJD(3) * t202 + t277 * t329;
t342 = t278 * t366;
t330 = t275 * t342;
t80 = qJD(5) * t152 + t162 * t276 + t279 * t330;
t81 = -qJD(5) * t153 + t162 * t279 - t276 * t330;
t472 = -t449 * t161 + t450 * t81 + t453 * t80;
t471 = -t451 * t161 + t453 * t81 + t455 * t80;
t385 = qJ(4) * t201;
t115 = t202 * t427 + t385;
t378 = qJ(6) + t427;
t384 = qJ(6) * t202;
t91 = t279 * t98;
t470 = t369 * t378 - t367 + pkin(5) * t201 - t91 - (-t115 - t384) * t276;
t246 = t378 * t279;
t35 = t279 * t115 + t276 * t98;
t469 = -qJD(5) * t246 - t276 * qJD(6) - t279 * t384 - t35;
t468 = pkin(5) * t368 - t202 * t353 - t463;
t114 = pkin(3) * t255 - t463;
t195 = Ifges(4,4) * t201;
t127 = t202 * Ifges(4,1) - t255 * Ifges(4,5) - t195;
t358 = -Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1;
t359 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t363 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t462 = -Ifges(5,4) / 0.2e1;
t364 = Ifges(4,5) / 0.2e1 + t462;
t194 = Ifges(5,6) * t201;
t457 = -t255 / 0.2e1;
t458 = -t202 / 0.2e1;
t422 = Ifges(5,4) * t457 + Ifges(5,2) * t458 + t194 / 0.2e1;
t62 = t153 * Ifges(7,5) + t152 * Ifges(7,6) + t196 * Ifges(7,3);
t63 = t153 * Ifges(6,5) + t152 * Ifges(6,6) + t196 * Ifges(6,3);
t284 = t152 * t359 - t153 * t363 - t196 * t358 - t255 * t364 - t110 * mrSges(5,3) - t12 * mrSges(7,2) - t21 * mrSges(6,2) - t422 + t127 / 0.2e1 + t62 / 0.2e1 + t63 / 0.2e1 + t114 * mrSges(5,1) + t128 * mrSges(4,3) + t182 * mrSges(4,2) + t20 * mrSges(6,1) + t9 * mrSges(7,1);
t361 = Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1;
t467 = (-t201 * t365 + t202 * t361 + t284) * t280;
t466 = -t449 * t276 + t451 * t279;
t465 = -t450 * t276 + t476;
t464 = t455 * t279 - t475;
t417 = -t161 / 0.2e1;
t459 = -t162 / 0.2e1;
t448 = -Ifges(6,3) - Ifges(7,3);
t13 = Ifges(7,5) * t80 + Ifges(7,6) * t81 - Ifges(7,3) * t161;
t14 = Ifges(6,5) * t80 + Ifges(6,6) * t81 - Ifges(6,3) * t161;
t447 = t13 + t14;
t352 = pkin(1) * t387;
t232 = -pkin(8) * t383 + t281 * t352;
t223 = t232 * qJD(2);
t210 = qJD(1) * t223;
t444 = t210 * mrSges(3,2);
t339 = t387 * t280;
t228 = -t339 + t357;
t212 = -pkin(2) * t387 - t232;
t229 = t277 * t387 + t280 * t383;
t292 = -t229 * qJ(4) + t212;
t111 = t228 * t427 + t292;
t382 = t275 * t281;
t233 = pkin(8) * t382 + t278 * t352;
t213 = pkin(9) * t387 + t233;
t144 = -t277 * t213 + t214 * t280;
t133 = pkin(3) * t382 - t144;
t99 = pkin(4) * t229 + pkin(10) * t382 + t133;
t37 = t279 * t111 + t276 * t99;
t362 = Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t443 = t201 * t362;
t138 = -t161 * mrSges(5,1) + mrSges(5,2) * t330;
t221 = qJD(2) * t298;
t209 = qJD(1) * t221;
t58 = -t183 * t370 - t193 * t371 + t209 * t280 - t277 * t210;
t55 = -pkin(3) * t330 - t58;
t442 = m(5) * t55 + t138;
t441 = mrSges(3,1) * t301 - mrSges(4,1) * t201 - mrSges(4,2) * t202 - mrSges(3,3) * t351;
t435 = t454 * t161 + t452 * t162 + t330 * t456;
t372 = qJD(2) * t275;
t349 = t278 * t372;
t302 = t427 * t349;
t32 = -pkin(4) * t161 - qJD(1) * t302 - t58;
t224 = t233 * qJD(2);
t211 = qJD(1) * t224;
t286 = t161 * qJ(4) - t202 * qJD(4) + t211;
t38 = t162 * t427 + t286;
t5 = t276 * t32 + t279 * t38 + t69 * t368 - t369 * t74;
t6 = -qJD(5) * t21 - t276 * t38 + t279 * t32;
t434 = t276 * t5 + t279 * t6;
t430 = t80 / 0.2e1;
t429 = t81 / 0.2e1;
t428 = Ifges(5,2) * t417 + Ifges(5,6) * t459 + t330 * t462;
t421 = -t152 / 0.2e1;
t419 = -t153 / 0.2e1;
t411 = -t196 / 0.2e1;
t400 = Ifges(3,4) * t278;
t399 = Ifges(3,4) * t281;
t394 = Ifges(3,2) * t278;
t391 = t219 * mrSges(3,3);
t390 = t222 * mrSges(3,3);
t389 = t278 * Ifges(3,1);
t165 = mrSges(5,1) * t201 + mrSges(5,3) * t255;
t96 = -mrSges(6,1) * t152 + mrSges(6,2) * t153;
t388 = -t165 + t96;
t386 = Ifges(3,6) * qJD(2);
t380 = t279 * t280;
t163 = mrSges(4,2) * t255 - mrSges(4,3) * t201;
t375 = t163 - t165;
t164 = -mrSges(4,1) * t255 - mrSges(4,3) * t202;
t166 = mrSges(5,1) * t202 - mrSges(5,2) * t255;
t374 = -t166 + t164;
t145 = t280 * t213 + t277 * t214;
t257 = t280 * pkin(9) + t274;
t348 = t281 * t372;
t24 = -t81 * mrSges(7,1) + t80 * mrSges(7,2);
t345 = Ifges(3,5) * t387;
t344 = Ifges(3,6) * t387;
t36 = -t111 * t276 + t279 * t99;
t1 = -pkin(5) * t161 - qJ(6) * t80 - qJD(6) * t153 + t6;
t2 = qJ(6) * t81 + qJD(6) * t152 + t5;
t323 = -t1 * t279 - t2 * t276;
t319 = mrSges(6,1) * t276 + mrSges(6,2) * t279;
t317 = mrSges(7,1) * t276 + mrSges(7,2) * t279;
t132 = qJ(4) * t382 - t145;
t84 = -t213 * t370 - t214 * t371 + t221 * t280 - t277 * t223;
t296 = -t228 * t276 + t275 * t379;
t173 = t228 * t279 + t275 * t381;
t171 = -qJD(3) * t339 - t280 * t348 + t332;
t42 = -pkin(4) * t171 - t302 - t84;
t172 = qJD(3) * t229 + t277 * t348;
t288 = t171 * qJ(4) - t229 * qJD(4) + t224;
t53 = t172 * t427 + t288;
t7 = -t111 * t369 + t276 * t42 + t279 * t53 + t99 * t368;
t57 = -t183 * t371 + t193 * t370 + t277 * t209 + t280 * t210;
t83 = -t213 * t371 + t214 * t370 + t277 * t221 + t280 * t223;
t112 = -pkin(4) * t228 - t132;
t293 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t50 = -qJ(4) * t330 + t255 * qJD(4) - t57;
t8 = -qJD(5) * t37 - t276 * t53 + t279 * t42;
t30 = -pkin(4) * t162 - t50;
t59 = -qJ(4) * t349 + qJD(4) * t382 - t83;
t290 = (t344 + (Ifges(3,2) * t281 + t400) * t275) * qJD(1);
t43 = -pkin(4) * t172 - t59;
t271 = pkin(5) * t276 + qJ(4);
t260 = Ifges(3,4) * t350;
t253 = Ifges(3,5) * t329;
t248 = -pkin(3) * t280 + t340;
t245 = t378 * t276;
t237 = t279 * t256;
t227 = pkin(5) * t380 + t257;
t226 = -qJ(4) * t370 + t336;
t218 = -mrSges(3,2) * t301 + mrSges(3,3) * t350;
t179 = Ifges(3,1) * t351 + Ifges(3,5) * t301 + t260;
t178 = t290 + t386;
t168 = -t234 * t276 + t237;
t154 = -qJ(6) * t380 + t169;
t149 = -qJ(4) * t333 + t355;
t148 = pkin(5) * t277 + t276 * t337 + t237;
t143 = -mrSges(5,2) * t201 - mrSges(5,3) * t202;
t141 = pkin(3) * t202 + t385;
t140 = -mrSges(4,2) * t330 - mrSges(4,3) * t162;
t139 = mrSges(4,1) * t330 + mrSges(4,3) * t161;
t137 = mrSges(5,1) * t162 - mrSges(5,3) * t330;
t136 = -pkin(3) * t351 - t146;
t131 = t228 * pkin(3) + t292;
t126 = -t255 * Ifges(5,1) - t202 * Ifges(5,4) + t201 * Ifges(5,5);
t123 = t202 * Ifges(4,5) - t201 * Ifges(4,6) - t255 * Ifges(4,3);
t120 = mrSges(6,1) * t196 - mrSges(6,3) * t153;
t119 = mrSges(7,1) * t196 - mrSges(7,3) * t153;
t118 = -mrSges(6,2) * t196 + mrSges(6,3) * t152;
t117 = -mrSges(7,2) * t196 + mrSges(7,3) * t152;
t108 = qJD(5) * t296 + t172 * t279 - t276 * t349;
t107 = qJD(5) * t173 + t172 * t276 + t279 * t349;
t101 = mrSges(4,1) * t162 - mrSges(4,2) * t161;
t100 = -mrSges(5,2) * t162 + mrSges(5,3) * t161;
t95 = -mrSges(7,1) * t152 + mrSges(7,2) * t153;
t89 = -t161 * Ifges(4,1) - t162 * Ifges(4,4) + Ifges(4,5) * t330;
t88 = -t161 * Ifges(4,4) - t162 * Ifges(4,2) + Ifges(4,6) * t330;
t86 = Ifges(5,5) * t330 + t161 * Ifges(5,6) + t162 * Ifges(5,3);
t73 = t172 * pkin(3) + t288;
t71 = -pkin(3) * t349 - t84;
t70 = -pkin(5) * t173 + t112;
t54 = t162 * pkin(3) + t286;
t48 = mrSges(6,2) * t161 + mrSges(6,3) * t81;
t47 = mrSges(7,2) * t161 + mrSges(7,3) * t81;
t46 = -mrSges(6,1) * t161 - mrSges(6,3) * t80;
t45 = -mrSges(7,1) * t161 - mrSges(7,3) * t80;
t34 = -t115 * t276 + t91;
t26 = qJ(6) * t173 + t37;
t25 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t23 = pkin(5) * t229 + qJ(6) * t296 + t36;
t19 = -pkin(5) * t108 + t43;
t10 = -pkin(5) * t81 + t30;
t4 = qJ(6) * t108 + qJD(6) * t173 + t7;
t3 = -pkin(5) * t171 - qJ(6) * t107 + qJD(6) * t296 + t8;
t15 = [((t126 + t123) * t278 + t301 * (Ifges(3,5) * t281 - Ifges(3,6) * t278) + t281 * t179 + (t281 * (t345 + (t389 + t399) * t275) + (t228 * t452 - t229 * t454 - t382 * t456) * t278) * qJD(1)) * t372 / 0.2e1 + t1 * (mrSges(7,1) * t229 + mrSges(7,3) * t296) + t6 * (mrSges(6,1) * t229 + mrSges(6,3) * t296) + t10 * (-mrSges(7,1) * t173 - mrSges(7,2) * t296) + t30 * (-mrSges(6,1) * t173 - mrSges(6,2) * t296) + (t173 * t453 + t229 * t451 - t296 * t455) * t430 + (-Ifges(4,4) * t228 - Ifges(4,5) * t382 - t451 * t296 + t449 * t173 + (Ifges(4,1) - t448) * t229) * t417 + (t173 * t450 + t229 * t449 - t296 * t453) * t429 + m(4) * (-t128 * t84 + t129 * t83 + t144 * t58 + t145 * t57) - t435 * t382 / 0.2e1 + (-m(3) * t232 + m(4) * t212 - mrSges(3,1) * t387 + mrSges(4,1) * t228 + mrSges(4,2) * t229) * t211 - (t290 + t178) * t349 / 0.2e1 - (t127 + t63 + t62) * t171 / 0.2e1 + ((-t394 + t399) * t343 + (Ifges(3,1) * t281 - t400) * t342 - 0.2e1 * pkin(1) * (mrSges(3,1) * t278 + mrSges(3,2) * t281) * t366) * t275 ^ 2 + m(5) * (t110 * t73 + t114 * t71 + t116 * t59 + t131 * t54 + t132 * t50 + t133 * t55) + m(7) * (t1 * t23 + t10 * t70 + t12 * t4 + t19 * t44 + t2 * t26 + t3 * t9) + m(6) * (t112 * t30 + t20 * t8 + t21 * t7 + t36 * t6 + t37 * t5 + t43 * t82) + t108 * t473 + (t107 * t455 + t108 * t453 - t171 * t451) * t418 + (t210 * t382 + t211 * t383 - t232 * t329 - t233 * t330) * mrSges(3,3) + (-m(3) * t219 + m(4) * t182 - t441) * t224 + (t89 + t447) * t229 / 0.2e1 + (t107 * t451 + t108 * t449 + t171 * t448) * t410 + (t107 * t453 + t108 * t450 - t171 * t449) * t420 + t445 * t107 / 0.2e1 + m(3) * (t210 * t233 + t222 * t223) - t387 * t444 - t349 * t390 - t348 * t391 + t387 * (-Ifges(3,6) * t330 + t253) / 0.2e1 + t50 * (mrSges(5,1) * t228 + mrSges(5,3) * t382) + t57 * (mrSges(4,2) * t382 - mrSges(4,3) * t228) + t162 * (-Ifges(5,5) * t382 - Ifges(5,6) * t229 + Ifges(5,3) * t228) / 0.2e1 + t161 * (-Ifges(5,4) * t382 - Ifges(5,2) * t229 + Ifges(5,6) * t228) / 0.2e1 + t55 * (mrSges(5,1) * t229 - mrSges(5,2) * t382) + t58 * (-mrSges(4,1) * t382 - mrSges(4,3) * t229) + t116 * (mrSges(5,1) * t172 - mrSges(5,3) * t349) + t129 * (-mrSges(4,2) * t349 - mrSges(4,3) * t172) - t201 * (-Ifges(4,4) * t171 - Ifges(4,2) * t172 + Ifges(4,6) * t349) / 0.2e1 + t201 * (Ifges(5,5) * t349 + Ifges(5,6) * t171 + Ifges(5,3) * t172) / 0.2e1 + t202 * (-Ifges(4,1) * t171 - Ifges(4,4) * t172 + Ifges(4,5) * t349) / 0.2e1 + t114 * (-mrSges(5,1) * t171 + mrSges(5,2) * t349) - t128 * (mrSges(4,1) * t349 + mrSges(4,3) * t171) - t228 * t88 / 0.2e1 + t2 * (-mrSges(7,2) * t229 + mrSges(7,3) * t173) + t5 * (-mrSges(6,2) * t229 + mrSges(6,3) * t173) + t54 * (-mrSges(5,2) * t228 - mrSges(5,3) * t229) + t228 * t86 / 0.2e1 + t212 * t101 + t223 * t218 + t182 * (mrSges(4,1) * t172 - mrSges(4,2) * t171) + t110 * (-mrSges(5,2) * t172 + mrSges(5,3) * t171) + t172 * t122 / 0.2e1 - t172 * t125 / 0.2e1 + t21 * (mrSges(6,2) * t171 + mrSges(6,3) * t108) + t12 * (mrSges(7,2) * t171 + mrSges(7,3) * t108) + t20 * (-mrSges(6,1) * t171 - mrSges(6,3) * t107) + t9 * (-mrSges(7,1) * t171 - mrSges(7,3) * t107) + t71 * t166 + t83 * t163 + t84 * t164 + t59 * t165 + t73 * t143 + t144 * t139 + t145 * t140 + t132 * t137 + t133 * t138 + t131 * t100 + t3 * t119 + t8 * t120 + (t171 * t454 + t172 * t452 + t349 * t456) * t457 + (Ifges(5,4) * t349 + Ifges(5,2) * t171 + Ifges(5,6) * t172) * t458 + (Ifges(4,4) * t229 - Ifges(4,2) * t228 - Ifges(4,6) * t382) * t459 + t4 * t117 + t7 * t118 + t82 * (-mrSges(6,1) * t108 + mrSges(6,2) * t107) + t44 * (-mrSges(7,1) * t108 + mrSges(7,2) * t107) + t112 * t25 + t19 * t95 + t43 * t96 + t70 * t24 + t23 * t45 + t36 * t46 + t26 * t47 + t37 * t48 - t471 * t296 / 0.2e1 + t472 * t173 / 0.2e1 + t171 * t422 + t229 * t428; t481 * t117 + (t1 * t148 + t10 * t227 + t12 * t481 + t154 * t2 + t44 * t478 + t482 * t9) * m(7) + t482 * t119 - t445 * t189 / 0.2e1 - t446 * t188 / 0.2e1 + ((-t386 / 0.2e1 + t390 + t129 * mrSges(4,2) + t116 * mrSges(5,3) - t123 / 0.2e1 - t126 / 0.2e1 + t178 / 0.2e1 - t114 * mrSges(5,2) + t128 * mrSges(4,1) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t255 - t364 * t202 + t360 * t201 + (t344 / 0.2e1 + (pkin(1) * mrSges(3,1) + t400 / 0.2e1) * t275) * qJD(1) + (-t454 * t277 - t452 * t280) * qJD(2) / 0.2e1) * t278 + (t391 + (-t345 / 0.2e1 + (-t389 / 0.2e1 + pkin(1) * mrSges(3,2) + t394 / 0.2e1) * t275) * qJD(1) - qJD(2) * Ifges(3,5) / 0.2e1 - t179 / 0.2e1 - t260 / 0.2e1 + (-t443 + t485) * t277 - t467) * t281) * t373 + t478 * t95 + t479 * t120 + (t168 * t6 + t169 * t5 + t20 * t479 + t21 * t480 + t257 * t30 + t477 * t82) * m(6) + t480 * t118 + t477 * t96 + (t467 + (t443 + t474) * t277 + (-t374 * t280 - t375 * t277 + m(5) * (t114 * t280 + t116 * t277) + m(4) * (t128 * t280 - t129 * t277)) * pkin(9)) * qJD(3) + m(5) * (t110 * t226 + t248 * t54) + (t226 - t149) * t143 - t444 + (-m(4) * t211 - t101) * pkin(2) + (-t188 * t21 + t189 * t20) * mrSges(6,3) + t441 * t222 + (t55 * mrSges(5,1) - t58 * mrSges(4,3) + t211 * mrSges(4,2) + t89 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 - t54 * mrSges(5,3) + t428 + t359 * t81 - t363 * t80 - t365 * t162 + (t358 - t361) * t161 + (-m(4) * t58 - t139 + t442) * pkin(9) + t293) * t277 - m(4) * (-t128 * t146 + t129 * t147 + t182 * t222) - m(5) * (t110 * t149 + t114 * t136 + t116 * t134) + (-t12 * t188 + t189 * t9) * mrSges(7,3) + t257 * t25 + t248 * t100 + t227 * t24 - t219 * t218 - t211 * mrSges(3,1) - t82 * (-mrSges(6,1) * t188 + mrSges(6,2) * t189) - t44 * (-mrSges(7,1) * t188 + mrSges(7,2) * t189) - t136 * t166 + t168 * t46 + t169 * t48 - t147 * t163 - t146 * t164 - t134 * t165 + t154 * t47 + t148 * t45 + t253 + (t88 / 0.2e1 - t86 / 0.2e1 - t50 * mrSges(5,1) + t57 * mrSges(4,3) - t211 * mrSges(4,1) + t54 * mrSges(5,2) + t10 * t318 + t30 * t320 - t362 * t162 + (-t276 * t363 + t279 * t359 - t365) * t161 + (t1 * t276 - t2 * t279) * mrSges(7,3) + (t276 * t6 - t279 * t5) * mrSges(6,3) + (m(4) * t57 - m(5) * t50 - t137 + t140) * pkin(9) + (-t82 * t319 - t44 * t317 + (t12 * t276 + t279 * t9) * mrSges(7,3) + (t20 * t279 + t21 * t276) * mrSges(6,3) + t465 * t421 + t464 * t419 + t466 * t411 + t276 * t473) * qJD(5) - t437 * t80 / 0.2e1 - t438 * t81 / 0.2e1 + t471 * t407 + (qJD(5) * t445 + t472) * t405) * t280 + (Ifges(7,5) * t189 + Ifges(7,6) * t188) * t411 + (Ifges(6,5) * t189 + Ifges(6,6) * t188) * t411 + (Ifges(7,1) * t189 + Ifges(7,4) * t188) * t419 + (Ifges(6,1) * t189 + Ifges(6,4) * t188) * t419 + (Ifges(7,4) * t189 + Ifges(7,2) * t188) * t421 + (Ifges(6,4) * t189 + Ifges(6,2) * t188) * t421; (-(-m(6) * t303 + t279 * t118 - t276 * t120) * t427 + t438 * t421 + t437 * t419 + t439 * t411 + t483) * qJD(5) - t434 * mrSges(6,3) + ((t361 - t362) * t201 - t474) * t202 - m(6) * (t20 * t34 + t21 * t35 - t299 * t82) + t299 * t96 + m(6) * (t30 * qJ(4) + t82 * qJD(4) - t427 * t434) - (t276 * t48 + t279 * t46) * t427 + (-t137 + t25) * qJ(4) + (t284 - t194 / 0.2e1 - t195 / 0.2e1) * t201 + t388 * qJD(4) + t375 * t128 + t374 * t129 + t271 * t24 - t245 * t47 - t246 * t45 + t435 - t141 * t143 - pkin(3) * t138 - t34 * t120 - t35 * t118 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t55 * mrSges(5,2) - t50 * mrSges(5,3) + (-pkin(3) * t55 - qJ(4) * t50 - t110 * t141 - t114 * t129 + t116 * t463) * m(5) + t464 * t430 + t465 * t429 + t466 * t417 + t468 * t95 + t469 * t117 + t470 * t119 + (-t1 * t246 + t10 * t271 + t12 * t469 - t2 * t245 + t44 * t468 + t470 * t9) * m(7) + t471 * t279 / 0.2e1 + t472 * t407 + t10 * t317 + t30 * t319 + t323 * mrSges(7,3); t202 * t143 + (t95 + t388) * t255 + (t45 + t46 + t196 * (t117 + t118)) * t279 + (t47 + t48 + t196 * (-t119 - t120)) * t276 - m(5) * (-t110 * t202 + t116 * t255) + (t196 * t321 + t255 * t44 - t323) * m(7) + (-t196 * t303 + t255 * t82 + t434) * m(6) + t442; t293 + (-(t11 - t9) * t12 + (-t153 * t44 + t1) * pkin(5)) * m(7) + (t12 * t153 + t152 * t9) * mrSges(7,3) + (t152 * t20 + t153 * t21) * mrSges(6,3) + t447 + (-t153 * t95 + t45) * pkin(5) - t82 * (mrSges(6,1) * t153 + mrSges(6,2) * t152) - t44 * (mrSges(7,1) * t153 + mrSges(7,2) * t152) + t21 * t120 - t11 * t117 - t20 * t118 + t12 * t119 + (t152 * t455 - t484) * t419 + t446 * t418 + (t152 * t451 - t153 * t449) * t411 + (-t153 * t450 + t445 + t486) * t421; -t152 * t117 + t153 * t119 + 0.2e1 * (t10 / 0.2e1 + t12 * t421 + t9 * t418) * m(7) + t24;];
tauc  = t15(:);
