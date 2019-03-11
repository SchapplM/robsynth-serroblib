% Calculate vector of centrifugal and Coriolis load on the joints for
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:35
% EndTime: 2019-03-09 17:42:30
% DurationCPUTime: 28.45s
% Computational Cost: add. (10854->797), mult. (28262->1056), div. (0->0), fcn. (20832->8), ass. (0->314)
t451 = Ifges(7,4) + Ifges(6,4);
t280 = cos(qJ(2));
t385 = cos(pkin(6));
t337 = t385 * qJD(1);
t330 = pkin(1) * t337;
t277 = sin(qJ(2));
t274 = sin(pkin(6));
t371 = qJD(1) * t274;
t349 = t277 * t371;
t219 = -pkin(8) * t349 + t280 * t330;
t297 = t274 * (pkin(2) * t277 - pkin(9) * t280);
t220 = qJD(1) * t297;
t276 = sin(qJ(3));
t279 = cos(qJ(3));
t146 = -t276 * t219 + t220 * t279;
t273 = t279 * pkin(4);
t368 = qJD(3) * t279;
t424 = pkin(4) + pkin(9);
t425 = pkin(3) + pkin(10);
t489 = -(t273 * t280 - t277 * t425) * t371 + t146 + t424 * t368;
t453 = Ifges(7,1) + Ifges(6,1);
t449 = Ifges(7,5) + Ifges(6,5);
t448 = Ifges(7,2) + Ifges(6,2);
t447 = Ifges(7,6) + Ifges(6,6);
t275 = sin(qJ(5));
t278 = cos(qJ(5));
t377 = t278 * t280;
t188 = (-t275 * t277 + t276 * t377) * t371;
t367 = qJD(5) * t275;
t488 = -t279 * t367 + t188;
t303 = pkin(10) * t276 - qJ(4) * t279;
t348 = t280 * t371;
t222 = pkin(8) * t348 + t277 * t330;
t333 = t276 * t348;
t353 = pkin(3) * t333 + t222;
t135 = t303 * t348 + t353;
t369 = qJD(3) * t276;
t335 = pkin(3) * t369 - qJD(4) * t276;
t200 = qJD(3) * t303 + t335;
t338 = -qJ(4) * t276 - pkin(2);
t234 = -t279 * t425 + t338;
t256 = t424 * t276;
t366 = qJD(5) * t278;
t479 = -t234 * t367 + t256 * t366 + (-t135 + t200) * t278 + t489 * t275;
t487 = t135 * t275 + t489 * t278;
t300 = t337 + qJD(2);
t293 = t279 * t300;
t201 = t276 * t349 - t293;
t255 = -qJD(3) + t348;
t152 = t201 * t278 + t255 * t275;
t486 = t451 * t152;
t182 = -pkin(2) * t300 - t219;
t202 = t276 * t300 + t279 * t349;
t286 = -t202 * qJ(4) + t182;
t110 = t201 * pkin(3) + t286;
t183 = pkin(9) * t300 + t222;
t214 = (-pkin(2) * t280 - pkin(9) * t277 - pkin(1)) * t274;
t193 = qJD(1) * t214;
t129 = t279 * t183 + t276 * t193;
t243 = t255 * qJ(4);
t116 = t243 - t129;
t122 = -t255 * Ifges(5,5) - t202 * Ifges(5,6) + t201 * Ifges(5,3);
t125 = t202 * Ifges(4,4) - t201 * Ifges(4,2) - t255 * Ifges(4,6);
t358 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1;
t359 = Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1;
t485 = t110 * mrSges(5,2) + t129 * mrSges(4,3) - t122 / 0.2e1 + t125 / 0.2e1 - t116 * mrSges(5,1) - t182 * mrSges(4,1) + t358 * t202 - t359 * t255;
t153 = t201 * t275 - t255 * t278;
t484 = t451 * t153;
t128 = t183 * t276 - t279 * t193;
t298 = pkin(4) * t202 + t128;
t69 = t255 * t425 + qJD(4) + t298;
t74 = t201 * t425 + t286;
t20 = -t275 * t74 + t278 * t69;
t21 = t275 * t69 + t278 * t74;
t302 = t20 * t275 - t21 * t278;
t317 = mrSges(7,1) * t278 - mrSges(7,2) * t275;
t319 = mrSges(6,1) * t278 - mrSges(6,2) * t275;
t12 = qJ(6) * t152 + t21;
t11 = -qJ(6) * t153 + t20;
t196 = qJD(5) + t202;
t9 = pkin(5) * t196 + t11;
t320 = t12 * t278 - t275 * t9;
t403 = -t278 / 0.2e1;
t405 = -t275 / 0.2e1;
t98 = -pkin(4) * t201 + t129;
t82 = -t243 + t98;
t44 = -pkin(5) * t152 + qJD(6) + t82;
t443 = t153 * t453 + t449 * t196 + t486;
t444 = t152 * t448 + t196 * t447 + t484;
t483 = t302 * mrSges(6,3) - t320 * mrSges(7,3) + t44 * t317 + t82 * t319 + t403 * t444 + t405 * t443;
t454 = Ifges(5,1) + Ifges(4,3);
t482 = Ifges(5,4) - Ifges(4,5);
t450 = Ifges(5,5) - Ifges(4,6);
t379 = t275 * t280;
t189 = (t276 * t379 + t277 * t278) * t371;
t332 = t279 * t348;
t336 = qJ(6) * t279 - t234;
t481 = qJ(6) * t189 + t336 * t366 + (-qJ(6) * t369 - qJD(5) * t256 + qJD(6) * t279 - t200) * t275 + t487 + (-t332 + t368) * pkin(5);
t365 = t278 * qJD(6);
t480 = -t279 * t365 + t479 + (t278 * t369 - t488) * qJ(6);
t169 = t278 * t234 + t275 * t256;
t478 = -qJD(5) * t169 - t200 * t275 + t487;
t147 = t279 * t219 + t276 * t220;
t134 = -qJ(4) * t349 - t147;
t121 = -pkin(4) * t333 - t134;
t351 = -pkin(5) * t278 - pkin(4);
t477 = -t121 + (-pkin(9) + t351) * t369 + t488 * pkin(5);
t476 = -t424 * t369 - t121;
t475 = t451 * t278;
t474 = t451 * t275;
t461 = -qJD(4) - t128;
t408 = t196 / 0.2e1;
t416 = t153 / 0.2e1;
t418 = t152 / 0.2e1;
t435 = t275 * t453 + t475;
t436 = t278 * t448 + t474;
t437 = t275 * t449 + t278 * t447;
t473 = -t408 * t437 - t416 * t435 - t418 * t436 + t483 + t485;
t472 = -t201 / 0.2e1;
t471 = t444 / 0.2e1;
t364 = qJD(1) * qJD(2);
t341 = t280 * t364;
t328 = t274 * t341;
t381 = t274 * t277;
t355 = t276 * t381;
t331 = qJD(3) * t355;
t161 = qJD(1) * t331 - qJD(3) * t293 - t279 * t328;
t162 = qJD(3) * t202 + t276 * t328;
t340 = t277 * t364;
t329 = t274 * t340;
t80 = qJD(5) * t152 + t162 * t275 + t278 * t329;
t81 = -qJD(5) * t153 + t162 * t278 - t275 * t329;
t470 = -t161 * t447 + t448 * t81 + t451 * t80;
t469 = -t449 * t161 + t451 * t81 + t453 * t80;
t383 = qJ(4) * t201;
t115 = t202 * t425 + t383;
t376 = qJ(6) + t425;
t382 = qJ(6) * t202;
t91 = t278 * t98;
t468 = t367 * t376 - t365 + pkin(5) * t201 - t91 - (-t115 - t382) * t275;
t246 = t376 * t278;
t35 = t278 * t115 + t275 * t98;
t467 = -qJD(5) * t246 - t275 * qJD(6) - t278 * t382 - t35;
t466 = pkin(5) * t366 - t202 * t351 - t461;
t114 = pkin(3) * t255 - t461;
t195 = Ifges(4,4) * t201;
t127 = t202 * Ifges(4,1) - t255 * Ifges(4,5) - t195;
t356 = -Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1;
t357 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t362 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t460 = -Ifges(5,4) / 0.2e1;
t363 = Ifges(4,5) / 0.2e1 + t460;
t455 = t202 / 0.2e1;
t456 = Ifges(5,6) * t472;
t420 = t255 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t455 + t456;
t62 = t153 * Ifges(7,5) + t152 * Ifges(7,6) + t196 * Ifges(7,3);
t63 = t153 * Ifges(6,5) + t152 * Ifges(6,6) + t196 * Ifges(6,3);
t283 = -t152 * t357 - t153 * t362 - t196 * t356 - t255 * t363 + t114 * mrSges(5,1) + t128 * mrSges(4,3) + t182 * mrSges(4,2) + t20 * mrSges(6,1) + t9 * mrSges(7,1) + t420 + t127 / 0.2e1 + t62 / 0.2e1 + t63 / 0.2e1 - t110 * mrSges(5,3) - t12 * mrSges(7,2) - t21 * mrSges(6,2);
t360 = Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1;
t465 = (-t201 * t358 + t202 * t360 + t283) * t279;
t464 = -t275 * t447 + t278 * t449;
t463 = -t275 * t448 + t475;
t462 = t278 * t453 - t474;
t415 = -t161 / 0.2e1;
t457 = -t162 / 0.2e1;
t446 = Ifges(6,3) + Ifges(7,3);
t13 = Ifges(7,5) * t80 + Ifges(7,6) * t81 - Ifges(7,3) * t161;
t14 = Ifges(6,5) * t80 + Ifges(6,6) * t81 - Ifges(6,3) * t161;
t445 = t13 + t14;
t350 = pkin(1) * t385;
t232 = -pkin(8) * t381 + t280 * t350;
t223 = t232 * qJD(2);
t210 = qJD(1) * t223;
t442 = t210 * mrSges(3,2);
t228 = -t279 * t385 + t355;
t212 = -pkin(2) * t385 - t232;
t229 = t276 * t385 + t279 * t381;
t291 = -t229 * qJ(4) + t212;
t111 = t228 * t425 + t291;
t380 = t274 * t280;
t233 = pkin(8) * t380 + t277 * t350;
t213 = pkin(9) * t385 + t233;
t144 = -t276 * t213 + t214 * t279;
t133 = pkin(3) * t380 - t144;
t99 = pkin(4) * t229 + pkin(10) * t380 + t133;
t37 = t278 * t111 + t275 * t99;
t361 = Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t441 = t201 * t361;
t138 = -t161 * mrSges(5,1) + mrSges(5,2) * t329;
t221 = qJD(2) * t297;
t209 = qJD(1) * t221;
t58 = -t183 * t368 - t193 * t369 + t209 * t279 - t276 * t210;
t55 = -pkin(3) * t329 - t58;
t440 = m(5) * t55 + t138;
t439 = mrSges(3,1) * t300 - mrSges(4,1) * t201 - mrSges(4,2) * t202 - mrSges(3,3) * t349;
t433 = t482 * t161 + t450 * t162 + t329 * t454;
t370 = qJD(2) * t274;
t347 = t277 * t370;
t301 = t425 * t347;
t32 = -pkin(4) * t161 - qJD(1) * t301 - t58;
t224 = t233 * qJD(2);
t211 = qJD(1) * t224;
t285 = t161 * qJ(4) - t202 * qJD(4) + t211;
t38 = t162 * t425 + t285;
t5 = t275 * t32 + t278 * t38 + t69 * t366 - t367 * t74;
t6 = -qJD(5) * t21 - t275 * t38 + t278 * t32;
t432 = t275 * t5 + t278 * t6;
t428 = t80 / 0.2e1;
t427 = t81 / 0.2e1;
t426 = Ifges(5,2) * t415 + Ifges(5,6) * t457 + t329 * t460;
t419 = -t152 / 0.2e1;
t417 = -t153 / 0.2e1;
t409 = -t196 / 0.2e1;
t398 = Ifges(3,4) * t277;
t397 = Ifges(3,4) * t280;
t392 = Ifges(3,2) * t277;
t389 = t219 * mrSges(3,3);
t388 = t222 * mrSges(3,3);
t387 = t277 * Ifges(3,1);
t165 = mrSges(5,1) * t201 + mrSges(5,3) * t255;
t96 = -mrSges(6,1) * t152 + mrSges(6,2) * t153;
t386 = -t165 + t96;
t384 = Ifges(3,6) * qJD(2);
t378 = t278 * t279;
t163 = mrSges(4,2) * t255 - mrSges(4,3) * t201;
t373 = t163 - t165;
t164 = -mrSges(4,1) * t255 - mrSges(4,3) * t202;
t166 = mrSges(5,1) * t202 - mrSges(5,2) * t255;
t372 = -t166 + t164;
t145 = t279 * t213 + t276 * t214;
t257 = t279 * pkin(9) + t273;
t346 = t280 * t370;
t24 = -t81 * mrSges(7,1) + t80 * mrSges(7,2);
t343 = Ifges(3,5) * t385;
t342 = Ifges(3,6) * t385;
t36 = -t111 * t275 + t278 * t99;
t1 = -pkin(5) * t161 - qJ(6) * t80 - qJD(6) * t153 + t6;
t2 = qJ(6) * t81 + qJD(6) * t152 + t5;
t322 = -t1 * t278 - t2 * t275;
t318 = mrSges(6,1) * t275 + mrSges(6,2) * t278;
t316 = mrSges(7,1) * t275 + mrSges(7,2) * t278;
t132 = qJ(4) * t380 - t145;
t84 = -t213 * t368 - t214 * t369 + t221 * t279 - t276 * t223;
t295 = -t228 * t275 + t274 * t377;
t173 = t228 * t278 + t274 * t379;
t172 = -t331 + (qJD(3) * t385 + t346) * t279;
t42 = pkin(4) * t172 - t301 - t84;
t171 = qJD(3) * t229 + t276 * t346;
t287 = -t172 * qJ(4) - t229 * qJD(4) + t224;
t53 = t171 * t425 + t287;
t7 = -t111 * t367 + t275 * t42 + t278 * t53 + t99 * t366;
t57 = -t183 * t369 + t193 * t368 + t276 * t209 + t279 * t210;
t83 = -t213 * t369 + t214 * t368 + t276 * t221 + t279 * t223;
t112 = -pkin(4) * t228 - t132;
t292 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t50 = -qJ(4) * t329 + t255 * qJD(4) - t57;
t8 = -qJD(5) * t37 - t275 * t53 + t278 * t42;
t30 = -pkin(4) * t162 - t50;
t59 = -qJ(4) * t347 + qJD(4) * t380 - t83;
t289 = (t342 + (Ifges(3,2) * t280 + t398) * t274) * qJD(1);
t43 = -pkin(4) * t171 - t59;
t270 = pkin(5) * t275 + qJ(4);
t260 = Ifges(3,4) * t348;
t253 = Ifges(3,5) * t328;
t248 = -pkin(3) * t279 + t338;
t245 = t376 * t275;
t237 = t278 * t256;
t227 = pkin(5) * t378 + t257;
t226 = -qJ(4) * t368 + t335;
t218 = -mrSges(3,2) * t300 + mrSges(3,3) * t348;
t179 = Ifges(3,1) * t349 + Ifges(3,5) * t300 + t260;
t178 = t289 + t384;
t168 = -t234 * t275 + t237;
t154 = -qJ(6) * t378 + t169;
t149 = -qJ(4) * t332 + t353;
t148 = pkin(5) * t276 + t275 * t336 + t237;
t143 = -mrSges(5,2) * t201 - mrSges(5,3) * t202;
t141 = pkin(3) * t202 + t383;
t140 = -mrSges(4,2) * t329 - mrSges(4,3) * t162;
t139 = mrSges(4,1) * t329 + mrSges(4,3) * t161;
t137 = mrSges(5,1) * t162 - mrSges(5,3) * t329;
t136 = -pkin(3) * t349 - t146;
t131 = t228 * pkin(3) + t291;
t126 = -t255 * Ifges(5,1) - t202 * Ifges(5,4) + t201 * Ifges(5,5);
t123 = t202 * Ifges(4,5) - t201 * Ifges(4,6) - t255 * Ifges(4,3);
t120 = mrSges(6,1) * t196 - mrSges(6,3) * t153;
t119 = mrSges(7,1) * t196 - mrSges(7,3) * t153;
t118 = -mrSges(6,2) * t196 + mrSges(6,3) * t152;
t117 = -mrSges(7,2) * t196 + mrSges(7,3) * t152;
t108 = qJD(5) * t173 + t171 * t275 + t278 * t347;
t107 = qJD(5) * t295 + t171 * t278 - t275 * t347;
t101 = mrSges(4,1) * t162 - mrSges(4,2) * t161;
t100 = -mrSges(5,2) * t162 + mrSges(5,3) * t161;
t95 = -mrSges(7,1) * t152 + mrSges(7,2) * t153;
t89 = -t161 * Ifges(4,1) - t162 * Ifges(4,4) + Ifges(4,5) * t329;
t88 = -t161 * Ifges(4,4) - t162 * Ifges(4,2) + Ifges(4,6) * t329;
t86 = Ifges(5,5) * t329 + t161 * Ifges(5,6) + t162 * Ifges(5,3);
t73 = t171 * pkin(3) + t287;
t71 = -pkin(3) * t347 - t84;
t70 = -pkin(5) * t173 + t112;
t54 = t162 * pkin(3) + t285;
t48 = mrSges(6,2) * t161 + mrSges(6,3) * t81;
t47 = mrSges(7,2) * t161 + mrSges(7,3) * t81;
t46 = -mrSges(6,1) * t161 - mrSges(6,3) * t80;
t45 = -mrSges(7,1) * t161 - mrSges(7,3) * t80;
t34 = -t115 * t275 + t91;
t26 = qJ(6) * t173 + t37;
t25 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t23 = pkin(5) * t229 + qJ(6) * t295 + t36;
t19 = -pkin(5) * t107 + t43;
t10 = -pkin(5) * t81 + t30;
t4 = qJ(6) * t107 + qJD(6) * t173 + t7;
t3 = pkin(5) * t172 - qJ(6) * t108 + qJD(6) * t295 + t8;
t15 = [((t126 + t123) * t277 + t300 * (Ifges(3,5) * t280 - Ifges(3,6) * t277) + t280 * t179 + ((t228 * t450 - t229 * t482 - t380 * t454) * t277 + t280 * (t343 + (t387 + t397) * t274)) * qJD(1)) * t370 / 0.2e1 - (t171 * t450 - t172 * t482 + t347 * t454) * t255 / 0.2e1 + (-Ifges(4,4) * t228 - Ifges(4,5) * t380 - t449 * t295 + t447 * t173 + (Ifges(4,1) + t446) * t229) * t415 + (t173 * t448 + t229 * t447 - t295 * t451) * t427 + (t173 * t451 + t229 * t449 - t295 * t453) * t428 + t1 * (mrSges(7,1) * t229 + mrSges(7,3) * t295) + t6 * (mrSges(6,1) * t229 + mrSges(6,3) * t295) + t10 * (-mrSges(7,1) * t173 - mrSges(7,2) * t295) + t30 * (-mrSges(6,1) * t173 - mrSges(6,2) * t295) + (-m(3) * t232 + m(4) * t212 - mrSges(3,1) * t385 + mrSges(4,1) * t228 + mrSges(4,2) * t229) * t211 + (-m(3) * t219 + m(4) * t182 - t439) * t224 + (t210 * t380 + t211 * t381 - t232 * t328 - t233 * t329) * mrSges(3,3) - t433 * t380 / 0.2e1 - t385 * t442 + m(7) * (t1 * t23 + t10 * t70 + t12 * t4 + t19 * t44 + t2 * t26 + t3 * t9) + m(6) * (t112 * t30 + t20 * t8 + t21 * t7 + t36 * t6 + t37 * t5 + t43 * t82) + m(5) * (t110 * t73 + t114 * t71 + t116 * t59 + t131 * t54 + t132 * t50 + t133 * t55) + (t107 * t447 + t108 * t449 + t172 * t446) * t408 + (t107 * t448 + t108 * t451 + t172 * t447) * t418 + (t107 * t451 + t108 * t453 + t172 * t449) * t416 + m(4) * (-t128 * t84 + t129 * t83 + t144 * t58 + t145 * t57) + m(3) * (t210 * t233 + t222 * t223) - t346 * t389 - t347 * t388 + t385 * (-Ifges(3,6) * t329 + t253) / 0.2e1 + t162 * (-Ifges(5,5) * t380 - Ifges(5,6) * t229 + Ifges(5,3) * t228) / 0.2e1 + t161 * (-Ifges(5,4) * t380 - Ifges(5,2) * t229 + Ifges(5,6) * t228) / 0.2e1 + t55 * (mrSges(5,1) * t229 - mrSges(5,2) * t380) + t58 * (-mrSges(4,1) * t380 - mrSges(4,3) * t229) + t50 * (mrSges(5,1) * t228 + mrSges(5,3) * t380) + t57 * (mrSges(4,2) * t380 - mrSges(4,3) * t228) - t128 * (mrSges(4,1) * t347 - mrSges(4,3) * t172) + t116 * (mrSges(5,1) * t171 - mrSges(5,3) * t347) + t129 * (-mrSges(4,2) * t347 - mrSges(4,3) * t171) + t201 * (Ifges(5,5) * t347 - Ifges(5,6) * t172 + Ifges(5,3) * t171) / 0.2e1 - t202 * (Ifges(5,4) * t347 - Ifges(5,2) * t172 + Ifges(5,6) * t171) / 0.2e1 + t114 * (mrSges(5,1) * t172 + mrSges(5,2) * t347) - t228 * t88 / 0.2e1 + t2 * (-mrSges(7,2) * t229 + mrSges(7,3) * t173) + t5 * (-mrSges(6,2) * t229 + mrSges(6,3) * t173) + t54 * (-mrSges(5,2) * t228 - mrSges(5,3) * t229) + t228 * t86 / 0.2e1 + t223 * t218 + t212 * t101 + t182 * (mrSges(4,1) * t171 + mrSges(4,2) * t172) + t21 * (-mrSges(6,2) * t172 + mrSges(6,3) * t107) + t12 * (-mrSges(7,2) * t172 + mrSges(7,3) * t107) + t20 * (mrSges(6,1) * t172 - mrSges(6,3) * t108) + t9 * (mrSges(7,1) * t172 - mrSges(7,3) * t108) + t110 * (-mrSges(5,2) * t171 - mrSges(5,3) * t172) + t171 * t122 / 0.2e1 - t171 * t125 / 0.2e1 + t71 * t166 + t83 * t163 + t84 * t164 + t59 * t165 + t73 * t143 + t144 * t139 + t145 * t140 + t131 * t100 + t132 * t137 + t133 * t138 + t3 * t119 + t8 * t120 + t4 * t117 + t7 * t118 + t44 * (-mrSges(7,1) * t107 + mrSges(7,2) * t108) + t82 * (-mrSges(6,1) * t107 + mrSges(6,2) * t108) + t112 * t25 + t19 * t95 + t43 * t96 + t70 * t24 + t26 * t47 + t37 * t48 + t23 * t45 + t36 * t46 - t469 * t295 / 0.2e1 + t470 * t173 / 0.2e1 + (Ifges(4,1) * t172 - Ifges(4,4) * t171 + Ifges(4,5) * t347) * t455 + (Ifges(4,4) * t229 - Ifges(4,2) * t228 - Ifges(4,6) * t380) * t457 + t107 * t471 + (Ifges(4,4) * t172 - Ifges(4,2) * t171 + Ifges(4,6) * t347) * t472 - (t289 + t178) * t347 / 0.2e1 + (t127 + t63 + t62) * t172 / 0.2e1 + ((Ifges(3,1) * t280 - t398) * t340 + (-t392 + t397) * t341 - 0.2e1 * pkin(1) * (mrSges(3,1) * t277 + mrSges(3,2) * t280) * t364) * t274 ^ 2 + t443 * t108 / 0.2e1 + (t89 + t445) * t229 / 0.2e1 + t172 * t420 + t229 * t426; t476 * t96 - t443 * t189 / 0.2e1 - t444 * t188 / 0.2e1 + ((t388 - t384 / 0.2e1 + t128 * mrSges(4,1) - t114 * mrSges(5,2) + t129 * mrSges(4,2) + t116 * mrSges(5,3) - t123 / 0.2e1 - t126 / 0.2e1 + t178 / 0.2e1 + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t255 - t363 * t202 + t359 * t201 + (t342 / 0.2e1 + (pkin(1) * mrSges(3,1) + t398 / 0.2e1) * t274) * qJD(1) + (-t276 * t482 - t450 * t279) * qJD(2) / 0.2e1) * t277 + (t389 - qJD(2) * Ifges(3,5) / 0.2e1 + (-t343 / 0.2e1 + (pkin(1) * mrSges(3,2) - t387 / 0.2e1 + t392 / 0.2e1) * t274) * qJD(1) - t260 / 0.2e1 - t179 / 0.2e1 + (-t441 + t485) * t276 - t465) * t280) * t371 + t480 * t117 + (t1 * t148 + t10 * t227 + t12 * t480 + t154 * t2 + t44 * t477 + t481 * t9) * m(7) + t481 * t119 + (-m(4) * t211 - t101) * pkin(2) + t478 * t120 + (t168 * t6 + t169 * t5 + t20 * t478 + t21 * t479 + t257 * t30 + t476 * t82) * m(6) + t479 * t118 + (-t188 * t21 + t189 * t20) * mrSges(6,3) + t439 * t222 + (-t58 * mrSges(4,3) + t55 * mrSges(5,1) + t211 * mrSges(4,2) + t89 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 - t54 * mrSges(5,3) + t426 - t357 * t81 - t362 * t80 - t358 * t162 + (t356 - t360) * t161 + (-m(4) * t58 - t139 + t440) * pkin(9) + t292) * t276 + t477 * t95 - t442 + m(5) * (t110 * t226 + t248 * t54) + (t226 - t149) * t143 + (t465 + (t441 - t473) * t276 + (-t372 * t279 - t373 * t276 + m(5) * (t114 * t279 + t116 * t276) + m(4) * (t128 * t279 - t129 * t276)) * pkin(9)) * qJD(3) - m(4) * (-t128 * t146 + t129 * t147 + t182 * t222) - m(5) * (t110 * t149 + t114 * t136 + t116 * t134) + (-t12 * t188 + t189 * t9) * mrSges(7,3) + t257 * t25 + t248 * t100 + t227 * t24 - t219 * t218 - t211 * mrSges(3,1) - t82 * (-mrSges(6,1) * t188 + mrSges(6,2) * t189) - t44 * (-mrSges(7,1) * t188 + mrSges(7,2) * t189) - t136 * t166 + t168 * t46 + t169 * t48 - t147 * t163 - t146 * t164 - t134 * t165 + t154 * t47 + t148 * t45 + t253 + (t88 / 0.2e1 - t86 / 0.2e1 + t30 * t319 + t10 * t317 + t57 * mrSges(4,3) - t50 * mrSges(5,1) - t211 * mrSges(4,1) + t54 * mrSges(5,2) - t361 * t162 + (-t275 * t362 - t278 * t357 - t358) * t161 + (t1 * t275 - t2 * t278) * mrSges(7,3) + (t275 * t6 - t278 * t5) * mrSges(6,3) + (m(4) * t57 - m(5) * t50 - t137 + t140) * pkin(9) + (-t44 * t316 - t82 * t318 + (t12 * t275 + t278 * t9) * mrSges(7,3) + (t20 * t278 + t21 * t275) * mrSges(6,3) + t463 * t419 + t462 * t417 + t464 * t409 + t275 * t471) * qJD(5) - t435 * t80 / 0.2e1 - t436 * t81 / 0.2e1 + t469 * t405 + (qJD(5) * t443 + t470) * t403) * t279 + (Ifges(7,5) * t189 + Ifges(7,6) * t188) * t409 + (Ifges(6,5) * t189 + Ifges(6,6) * t188) * t409 + (Ifges(7,1) * t189 + Ifges(7,4) * t188) * t417 + (Ifges(6,1) * t189 + Ifges(6,4) * t188) * t417 + (Ifges(7,4) * t189 + Ifges(7,2) * t188) * t419 + (Ifges(6,4) * t189 + Ifges(6,2) * t188) * t419; -t432 * mrSges(6,3) + (-(-m(6) * t302 + t278 * t118 - t275 * t120) * t425 + t436 * t419 + t435 * t417 + t437 * t409 + t483) * qJD(5) + t298 * t96 - m(6) * (t20 * t34 + t21 * t35 - t298 * t82) + m(6) * (t30 * qJ(4) + t82 * qJD(4) - t425 * t432) - (t275 * t48 + t278 * t46) * t425 + t466 * t95 + t467 * t117 + (-t1 * t246 + t10 * t270 + t12 * t467 - t2 * t245 + t44 * t466 + t468 * t9) * m(7) + t468 * t119 + t462 * t428 + t463 * t427 + t464 * t415 + (-pkin(3) * t55 - qJ(4) * t50 - t110 * t141 - t114 * t129 + t116 * t461) * m(5) + t433 + ((t360 - t361) * t201 + t473) * t202 + (t283 - t195 / 0.2e1 + t456) * t201 + (-t137 + t25) * qJ(4) + t386 * qJD(4) + t372 * t129 + t373 * t128 + t270 * t24 - t245 * t47 - t246 * t45 - t141 * t143 - pkin(3) * t138 - t34 * t120 - t35 * t118 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t55 * mrSges(5,2) - t50 * mrSges(5,3) + t469 * t278 / 0.2e1 + t470 * t405 + t10 * t316 + t30 * t318 + t322 * mrSges(7,3); t202 * t143 + (t95 + t386) * t255 + (t45 + t46 + t196 * (t117 + t118)) * t278 + (t47 + t48 + t196 * (-t119 - t120)) * t275 - m(5) * (-t110 * t202 + t116 * t255) + (t196 * t320 + t255 * t44 - t322) * m(7) + (-t196 * t302 + t255 * t82 + t432) * m(6) + t440; (t12 * t153 + t152 * t9) * mrSges(7,3) + (t152 * t20 + t153 * t21) * mrSges(6,3) + t292 + t445 + (-t153 * t95 + t45) * pkin(5) + (-(t11 - t9) * t12 + (-t153 * t44 + t1) * pkin(5)) * m(7) - t82 * (mrSges(6,1) * t153 + mrSges(6,2) * t152) - t44 * (mrSges(7,1) * t153 + mrSges(7,2) * t152) + t21 * t120 - t11 * t117 - t20 * t118 + t12 * t119 + (t152 * t453 - t484) * t417 + t444 * t416 + (t152 * t449 - t153 * t447) * t409 + (-t153 * t448 + t443 + t486) * t419; -t152 * t117 + t153 * t119 + 0.2e1 * (t10 / 0.2e1 + t12 * t419 + t9 * t416) * m(7) + t24;];
tauc  = t15(:);
