% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2018-11-23 17:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:40:16
% EndTime: 2018-11-23 17:40:42
% DurationCPUTime: 27.26s
% Computational Cost: add. (13602->841), mult. (35440->1144), div. (0->0), fcn. (26955->10), ass. (0->363)
t326 = cos(qJ(2));
t318 = sin(pkin(6));
t400 = qJD(1) * t318;
t380 = t326 * t400;
t295 = -qJD(3) + t380;
t322 = sin(qJ(3));
t325 = cos(qJ(3));
t417 = cos(pkin(6));
t368 = t417 * qJD(1);
t340 = t368 + qJD(2);
t335 = t325 * t340;
t323 = sin(qJ(2));
t381 = t323 * t400;
t237 = t322 * t381 - t335;
t448 = -t237 / 0.2e1;
t361 = pkin(1) * t368;
t257 = -pkin(8) * t381 + t326 * t361;
t337 = (pkin(2) * t323 - pkin(9) * t326) * t318;
t258 = qJD(1) * t337;
t180 = -t322 * t257 + t325 * t258;
t316 = t325 * pkin(4);
t435 = pkin(3) + qJ(5);
t374 = t435 * t323;
t397 = qJD(3) * t325;
t469 = pkin(4) + pkin(9);
t542 = -(t316 * t326 - t374) * t400 + t180 + t469 * t397;
t398 = qJD(3) * t322;
t366 = pkin(3) * t398 - qJD(4) * t322;
t260 = pkin(8) * t380 + t323 * t361;
t364 = t322 * t380;
t385 = pkin(3) * t364 + t260;
t541 = qJD(5) * t325 - t366 + t385 + t295 * (-qJ(4) * t325 + qJ(5) * t322);
t540 = Ifges(5,6) * t448;
t238 = t322 * t340 + t325 * t381;
t445 = t238 / 0.2e1;
t539 = t295 / 0.2e1;
t220 = pkin(9) * t340 + t260;
t251 = (-pkin(2) * t326 - pkin(9) * t323 - pkin(1)) * t318;
t230 = qJD(1) * t251;
t145 = t220 * t322 - t325 * t230;
t485 = -qJD(4) - t145;
t131 = pkin(3) * t295 - t485;
t233 = qJD(6) + t238;
t449 = t233 / 0.2e1;
t317 = sin(pkin(11));
t319 = cos(pkin(11));
t188 = t237 * t317 - t295 * t319;
t321 = sin(qJ(6));
t324 = cos(qJ(6));
t367 = t319 * t237 + t295 * t317;
t111 = t188 * t324 + t321 * t367;
t462 = t111 / 0.2e1;
t515 = -t188 * t321 + t324 * t367;
t464 = t515 / 0.2e1;
t531 = Ifges(7,5) * t462 + Ifges(7,6) * t464 + Ifges(7,3) * t449;
t522 = Ifges(5,4) * t539 + Ifges(5,2) * t445 + t531 + t540;
t538 = -t131 * mrSges(5,1) - t145 * mrSges(4,3) - t522;
t503 = Ifges(6,3) + Ifges(4,1);
t453 = t188 / 0.2e1;
t455 = t367 / 0.2e1;
t532 = Ifges(6,5) * t453 + Ifges(6,6) * t455;
t533 = Ifges(4,4) * t448;
t537 = t503 * t445 + t532 + t533;
t528 = t317 * t541 + t319 * t542;
t527 = t317 * t542 - t319 * t541;
t407 = t322 * t326;
t223 = (t317 * t407 + t319 * t323) * t400;
t363 = t325 * t380;
t413 = t317 * t322;
t535 = -pkin(5) * t363 + t223 * pkin(10) + (pkin(5) * t325 - pkin(10) * t413) * qJD(3) + t528;
t222 = (-t317 * t323 + t319 * t407) * t400;
t437 = pkin(10) * t319;
t534 = -pkin(10) * t222 + t398 * t437 + t527;
t443 = -t295 / 0.2e1;
t506 = Ifges(5,1) + Ifges(4,3);
t505 = Ifges(5,4) - Ifges(4,5);
t504 = Ifges(5,5) - Ifges(4,6);
t370 = -qJ(4) * t322 - pkin(2);
t273 = -t325 * t435 + t370;
t296 = t469 * t322;
t277 = t319 * t296;
t179 = pkin(5) * t322 + t277 + (pkin(10) * t325 - t273) * t317;
t202 = t319 * t273 + t317 * t296;
t408 = t319 * t325;
t189 = -pkin(10) * t408 + t202;
t100 = t179 * t324 - t189 * t321;
t530 = qJD(6) * t100 + t535 * t321 + t534 * t324;
t101 = t179 * t321 + t189 * t324;
t529 = -qJD(6) * t101 - t534 * t321 + t535 * t324;
t181 = t325 * t257 + t322 * t258;
t162 = -qJ(4) * t381 - t181;
t137 = -pkin(4) * t364 - t162;
t383 = -pkin(5) * t319 - pkin(4);
t526 = t222 * pkin(5) - t137 + (-pkin(9) + t383) * t398;
t525 = -t469 * t398 - t137;
t396 = qJD(1) * qJD(2);
t373 = t326 * t396;
t359 = t318 * t373;
t196 = qJD(3) * t238 + t322 * t359;
t372 = t323 * t396;
t360 = t318 * t372;
t159 = t196 * t319 - t317 * t360;
t412 = t318 * t323;
t390 = t322 * t412;
t362 = qJD(3) * t390;
t195 = qJD(1) * t362 - qJD(3) * t335 - t325 * t359;
t382 = pkin(1) * t417;
t411 = t318 * t326;
t272 = pkin(8) * t411 + t323 * t382;
t262 = t272 * qJD(2);
t247 = qJD(1) * t262;
t329 = t195 * qJ(4) - t238 * qJD(4) + t247;
t50 = t237 * qJD(5) + t196 * t435 + t329;
t354 = qJD(2) * t374;
t259 = qJD(2) * t337;
t245 = qJD(1) * t259;
t271 = -pkin(8) * t412 + t326 * t382;
t261 = t271 * qJD(2);
t246 = qJD(1) * t261;
t83 = -t220 * t397 - t230 * t398 + t245 * t325 - t322 * t246;
t51 = -pkin(4) * t195 + qJD(5) * t295 - t354 * t400 - t83;
t19 = t317 * t51 + t319 * t50;
t11 = pkin(10) * t159 + t19;
t338 = pkin(4) * t238 + t145;
t517 = qJD(4) + t338;
t89 = t295 * t435 + t517;
t219 = -pkin(2) * t340 - t257;
t330 = -t238 * qJ(4) + t219;
t94 = t237 * t435 + t330;
t32 = -t317 * t94 + t319 * t89;
t25 = pkin(5) * t238 - pkin(10) * t188 + t32;
t33 = t317 * t89 + t319 * t94;
t26 = pkin(10) * t367 + t33;
t5 = t25 * t324 - t26 * t321;
t160 = t196 * t317 + t319 * t360;
t18 = -t317 * t50 + t319 * t51;
t7 = -pkin(5) * t195 - pkin(10) * t160 + t18;
t1 = qJD(6) * t5 + t11 * t324 + t321 * t7;
t6 = t25 * t321 + t26 * t324;
t2 = -qJD(6) * t6 - t11 * t321 + t324 * t7;
t342 = t324 * t317 + t321 * t319;
t266 = t342 * qJD(6);
t336 = t342 * t238;
t403 = t336 + t266;
t486 = -t317 * t321 + t319 * t324;
t265 = t486 * qJD(6);
t491 = t486 * t238;
t520 = t265 + t491;
t523 = -t1 * t342 - t2 * t486 + t403 * t5 - t520 * t6;
t126 = t237 * pkin(3) + t330;
t521 = t32 * mrSges(6,1) + t5 * mrSges(7,1) + t219 * mrSges(4,2) - t33 * mrSges(6,2) - t6 * mrSges(7,2) - t126 * mrSges(5,3) + Ifges(4,5) * t443 + t531 + t537;
t146 = t325 * t220 + t322 * t230;
t285 = t295 * qJ(4);
t132 = t285 - t146;
t518 = t132 * mrSges(5,1) - t146 * mrSges(4,3);
t355 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t393 = -Ifges(5,6) / 0.2e1 - Ifges(4,4) / 0.2e1;
t452 = -t195 / 0.2e1;
t483 = Ifges(5,2) * t452;
t511 = -t196 / 0.2e1;
t75 = t196 * pkin(3) + t329;
t76 = -pkin(3) * t360 - t83;
t516 = t355 + t76 * mrSges(5,1) + t18 * mrSges(6,1) - t19 * mrSges(6,2) - t83 * mrSges(4,3) - t75 * mrSges(5,3) - Ifges(5,4) * t360 / 0.2e1 + t483 + Ifges(5,6) * t511 + t393 * t196;
t513 = -t219 * mrSges(4,1) + t126 * mrSges(5,2) + Ifges(5,5) * t539 + Ifges(4,6) * t443 + (Ifges(5,3) + Ifges(4,2)) * t448 + (Ifges(5,6) + Ifges(4,4)) * t445;
t446 = -t238 / 0.2e1;
t447 = t237 / 0.2e1;
t512 = -Ifges(5,2) * t446 - Ifges(5,6) * t447 - t505 * t443 + t521 + t537;
t42 = qJD(6) * t515 + t159 * t321 + t160 * t324;
t473 = t42 / 0.2e1;
t43 = -qJD(6) * t111 + t159 * t324 - t160 * t321;
t472 = t43 / 0.2e1;
t507 = t83 * mrSges(4,1);
t496 = t246 * mrSges(3,2);
t434 = -pkin(10) - t435;
t282 = t434 * t317;
t283 = t434 * t319;
t205 = -t282 * t321 + t283 * t324;
t118 = -pkin(4) * t237 + t146;
t113 = t319 * t118;
t415 = qJ(4) * t237;
t130 = t238 * t435 + t415;
t34 = -pkin(5) * t237 + t113 + (-pkin(10) * t238 - t130) * t317;
t60 = t317 * t118 + t319 * t130;
t45 = t238 * t437 + t60;
t495 = -qJD(5) * t342 + qJD(6) * t205 - t321 * t34 - t324 * t45;
t206 = t282 * t324 + t283 * t321;
t494 = -qJD(5) * t486 - qJD(6) * t206 + t321 * t45 - t324 * t34;
t489 = mrSges(3,1) * t340 - mrSges(4,1) * t237 - mrSges(4,2) * t238 - mrSges(3,3) * t381;
t487 = t195 * t505 + t196 * t504 + t360 * t506;
t395 = -Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t482 = t395 * t295 + t521 + t532 - t538;
t480 = -t145 * mrSges(4,1) - t146 * mrSges(4,2) + t131 * mrSges(5,2) - t260 * mrSges(3,3) - t132 * mrSges(5,3);
t478 = -Ifges(4,4) * t445 - Ifges(4,2) * t448 + Ifges(5,6) * t446 + Ifges(5,3) * t447 + t443 * t504 - t513;
t475 = Ifges(7,4) * t473 + Ifges(7,2) * t472 + Ifges(7,6) * t452;
t474 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t452;
t87 = t188 * Ifges(6,4) + Ifges(6,2) * t367 + Ifges(6,6) * t238;
t471 = t87 / 0.2e1;
t88 = t188 * Ifges(6,1) + Ifges(6,4) * t367 + Ifges(6,5) * t238;
t470 = t88 / 0.2e1;
t465 = -t515 / 0.2e1;
t463 = -t111 / 0.2e1;
t458 = t159 / 0.2e1;
t457 = t160 / 0.2e1;
t456 = -t367 / 0.2e1;
t454 = -t188 / 0.2e1;
t450 = -t233 / 0.2e1;
t442 = -t317 / 0.2e1;
t441 = -t319 / 0.2e1;
t440 = t319 / 0.2e1;
t38 = Ifges(7,5) * t42;
t37 = Ifges(7,6) * t43;
t438 = pkin(9) * t322;
t315 = t325 * pkin(9);
t369 = t417 * t325;
t399 = qJD(2) * t318;
t378 = t326 * t399;
t210 = -qJD(3) * t369 - t325 * t378 + t362;
t250 = pkin(9) * t417 + t272;
t99 = -t250 * t397 - t251 * t398 + t325 * t259 - t322 * t261;
t65 = -t210 * pkin(4) + (qJD(5) * t326 - t354) * t318 - t99;
t268 = t322 * t417 + t325 * t412;
t211 = qJD(3) * t268 + t322 * t378;
t267 = -t369 + t390;
t331 = t210 * qJ(4) - t268 * qJD(4) + t262;
t69 = t267 * qJD(5) + t211 * t435 + t331;
t24 = t317 * t65 + t319 * t69;
t433 = Ifges(3,4) * t323;
t432 = Ifges(3,4) * t326;
t431 = Ifges(6,4) * t317;
t430 = Ifges(6,4) * t319;
t429 = Ifges(7,4) * t111;
t428 = Ifges(3,2) * t323;
t425 = t159 * Ifges(6,6);
t424 = t160 * Ifges(6,5);
t420 = t257 * mrSges(3,3);
t418 = t323 * Ifges(3,1);
t416 = Ifges(3,6) * qJD(2);
t410 = t319 * t322;
t175 = -t322 * t250 + t325 * t251;
t158 = pkin(3) * t411 - t175;
t119 = t268 * pkin(4) + qJ(5) * t411 + t158;
t249 = -pkin(2) * t417 - t271;
t333 = -t268 * qJ(4) + t249;
t127 = t267 * t435 + t333;
t59 = t317 * t119 + t319 * t127;
t147 = t222 * t324 - t223 * t321;
t183 = t266 * t325 + t398 * t486;
t406 = t147 - t183;
t148 = t222 * t321 + t223 * t324;
t252 = t486 * t325;
t182 = -qJD(6) * t252 + t342 * t398;
t405 = t148 - t182;
t197 = mrSges(4,2) * t295 - mrSges(4,3) * t237;
t199 = mrSges(5,1) * t237 + mrSges(5,3) * t295;
t402 = t197 - t199;
t198 = -mrSges(4,1) * t295 - mrSges(4,3) * t238;
t200 = mrSges(5,1) * t238 - mrSges(5,2) * t295;
t401 = t198 - t200;
t176 = t325 * t250 + t322 * t251;
t297 = t315 + t316;
t8 = -Ifges(7,3) * t195 + t37 + t38;
t394 = Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1;
t392 = -Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1;
t114 = -mrSges(6,1) * t367 + mrSges(6,2) * t188;
t52 = -mrSges(7,1) * t515 + mrSges(7,2) * t111;
t391 = t114 - t199 + t52;
t379 = t323 * t399;
t14 = -t43 * mrSges(7,1) + t42 * mrSges(7,2);
t376 = Ifges(3,5) * t417;
t375 = Ifges(3,6) * t417;
t371 = t399 / 0.2e1;
t23 = -t317 * t69 + t319 * t65;
t167 = -t195 * mrSges(5,1) + mrSges(5,2) * t360;
t79 = -t159 * mrSges(6,1) + t160 * mrSges(6,2);
t57 = t319 * t119 - t127 * t317;
t356 = -Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(4,1) / 0.2e1;
t353 = qJD(1) * t371;
t351 = mrSges(6,1) * t319 - mrSges(6,2) * t317;
t350 = Ifges(6,1) * t317 + t430;
t349 = Ifges(6,2) * t319 + t431;
t348 = t18 * t319 + t19 * t317;
t347 = t317 * t32 - t319 * t33;
t96 = mrSges(6,2) * t195 + mrSges(6,3) * t159;
t97 = -mrSges(6,1) * t195 - mrSges(6,3) * t160;
t346 = t317 * t96 + t319 * t97;
t209 = t317 * t267 - t319 * t411;
t35 = pkin(5) * t268 - pkin(10) * t209 + t57;
t208 = t319 * t267 + t317 * t411;
t44 = pkin(10) * t208 + t59;
t12 = -t321 * t44 + t324 * t35;
t13 = t321 * t35 + t324 * t44;
t344 = t131 * t325 + t132 * t322;
t343 = t145 * t325 - t146 * t322;
t133 = t208 * t324 - t209 * t321;
t134 = t208 * t321 + t209 * t324;
t339 = t323 * t353;
t157 = qJ(4) * t411 - t176;
t95 = qJD(5) + t118 - t285;
t82 = -t220 * t398 + t230 * t397 + t322 * t245 + t325 * t246;
t98 = -t250 * t398 + t251 * t397 + t322 * t259 + t325 * t261;
t129 = -t267 * pkin(4) - t157;
t74 = -qJ(4) * t360 + t295 * qJD(4) - t82;
t334 = Ifges(6,5) * t317 / 0.2e1 + Ifges(6,6) * t440 + t393;
t54 = -pkin(4) * t196 - t74;
t84 = -qJ(4) * t379 + qJD(4) * t411 - t98;
t332 = (t375 + (Ifges(3,2) * t326 + t433) * t318) * qJD(1);
t70 = -t211 * pkin(4) - t84;
t328 = -t394 * t295 + t513 - t518;
t310 = pkin(5) * t317 + qJ(4);
t300 = Ifges(3,4) * t380;
t293 = Ifges(3,5) * t359;
t288 = -pkin(3) * t325 + t370;
t264 = pkin(5) * t408 + t297;
t263 = -qJ(4) * t397 + t366;
t256 = -mrSges(3,2) * t340 + mrSges(3,3) * t380;
t253 = t342 * t325;
t216 = Ifges(3,1) * t381 + Ifges(3,5) * t340 + t300;
t215 = t332 + t416;
t201 = -t273 * t317 + t277;
t184 = -qJ(4) * t363 + t385;
t178 = t211 * t317 + t319 * t379;
t177 = t211 * t319 - t317 * t379;
t174 = -mrSges(5,2) * t237 - mrSges(5,3) * t238;
t172 = pkin(3) * t238 + t415;
t169 = -mrSges(4,2) * t360 - mrSges(4,3) * t196;
t168 = mrSges(4,1) * t360 + mrSges(4,3) * t195;
t166 = mrSges(5,1) * t196 - mrSges(5,3) * t360;
t163 = -pkin(3) * t381 - t180;
t155 = t267 * pkin(3) + t333;
t143 = -t295 * Ifges(5,1) - t238 * Ifges(5,4) + t237 * Ifges(5,5);
t140 = t238 * Ifges(4,5) - t237 * Ifges(4,6) - t295 * Ifges(4,3);
t136 = mrSges(6,1) * t238 - mrSges(6,3) * t188;
t135 = -mrSges(6,2) * t238 + mrSges(6,3) * t367;
t122 = mrSges(4,1) * t196 - mrSges(4,2) * t195;
t121 = -mrSges(5,2) * t196 + mrSges(5,3) * t195;
t107 = -t195 * Ifges(4,1) - t196 * Ifges(4,4) + Ifges(4,5) * t360;
t106 = -t195 * Ifges(4,4) - t196 * Ifges(4,2) + Ifges(4,6) * t360;
t104 = Ifges(5,5) * t360 + t195 * Ifges(5,6) + t196 * Ifges(5,3);
t103 = Ifges(7,4) * t515;
t93 = t211 * pkin(3) + t331;
t91 = -pkin(3) * t379 - t99;
t90 = -t208 * pkin(5) + t129;
t85 = t238 * t383 - t145;
t78 = mrSges(7,1) * t233 - mrSges(7,3) * t111;
t77 = -mrSges(7,2) * t233 + mrSges(7,3) * t515;
t71 = -pkin(5) * t367 + t95;
t68 = t160 * Ifges(6,1) + t159 * Ifges(6,4) - t195 * Ifges(6,5);
t67 = Ifges(6,4) * t160 + Ifges(6,2) * t159 - Ifges(6,6) * t195;
t66 = -t195 * Ifges(6,3) + t424 + t425;
t58 = -t130 * t317 + t113;
t56 = -qJD(6) * t134 + t177 * t324 - t178 * t321;
t55 = qJD(6) * t133 + t177 * t321 + t178 * t324;
t47 = -t177 * pkin(5) + t70;
t41 = t111 * Ifges(7,1) + t233 * Ifges(7,5) + t103;
t40 = Ifges(7,2) * t515 + t233 * Ifges(7,6) + t429;
t29 = -pkin(5) * t159 + t54;
t28 = mrSges(7,2) * t195 + mrSges(7,3) * t43;
t27 = -mrSges(7,1) * t195 - mrSges(7,3) * t42;
t20 = pkin(10) * t177 + t24;
t17 = -pkin(5) * t210 - pkin(10) * t178 + t23;
t4 = -qJD(6) * t13 + t17 * t324 - t20 * t321;
t3 = qJD(6) * t12 + t17 * t321 + t20 * t324;
t9 = [(-t512 + t538) * t210 + (-Ifges(4,4) * t267 - Ifges(4,5) * t411 + Ifges(6,5) * t209 + Ifges(7,5) * t134 + Ifges(6,6) * t208 + Ifges(7,6) * t133 + (Ifges(7,3) + t503) * t268) * t452 + (Ifges(5,4) * t446 + Ifges(4,5) * t445 + Ifges(5,5) * t447 + Ifges(4,6) * t448 + t443 * t506 + t480) * t379 + (t267 * t504 - t411 * t506) * t339 - t411 * t507 - t487 * t411 / 0.2e1 + t196 * (-Ifges(5,5) * t411 + Ifges(5,3) * t267) / 0.2e1 + t178 * t470 + t177 * t471 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t449 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t464 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t472 + m(3) * (t246 * t272 + t260 * t261) + m(4) * (-t145 * t99 + t146 * t98 + t175 * t83 + t176 * t82) + (Ifges(6,1) * t209 + Ifges(6,4) * t208) * t457 + (Ifges(6,1) * t178 + Ifges(6,4) * t177) * t453 - (t332 + t215) * t379 / 0.2e1 + (t107 + t66 + t8) * t268 / 0.2e1 + (-0.2e1 * pkin(1) * (mrSges(3,1) * t323 + mrSges(3,2) * t326) * t396 + (-t428 + t432) * t373 + (Ifges(3,1) * t326 - t433) * t372) * t318 ^ 2 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t462 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t473 + (Ifges(6,5) * t178 + Ifges(6,6) * t177) * t445 + t195 * (-Ifges(5,4) * t411 + Ifges(5,6) * t267) / 0.2e1 + (Ifges(6,4) * t209 + Ifges(6,2) * t208) * t458 + (Ifges(6,4) * t178 + Ifges(6,2) * t177) * t455 + m(7) * (t1 * t13 + t12 * t2 + t29 * t90 + t3 * t6 + t4 * t5 + t47 * t71) + m(6) * (t129 * t54 + t18 * t57 + t19 * t59 + t23 * t32 + t24 * t33 + t70 * t95) + m(5) * (t126 * t93 + t131 * t91 + t132 * t84 + t155 * t75 + t157 * t74 + t158 * t76) + t134 * t474 + t133 * t475 + t267 * t104 / 0.2e1 - t267 * t106 / 0.2e1 + t261 * t256 + t249 * t122 + t208 * t67 / 0.2e1 + t54 * (-mrSges(6,1) * t208 + mrSges(6,2) * t209) + t209 * t68 / 0.2e1 + ((t143 + t140) * t323 + t340 * (Ifges(3,5) * t326 - Ifges(3,6) * t323) + t326 * t216) * t371 + t98 * t197 + t99 * t198 + t84 * t199 + t91 * t200 + t175 * t168 + t176 * t169 + t95 * (-mrSges(6,1) * t177 + mrSges(6,2) * t178) + t93 * t174 + t157 * t166 + t158 * t167 + t155 * t121 + (t1 * t133 - t134 * t2 - t5 * t55 + t56 * t6) * mrSges(7,3) + t29 * (-mrSges(7,1) * t133 + mrSges(7,2) * t134) + t24 * t135 + t23 * t136 + t129 * t79 + t70 * t114 + t59 * t96 + t57 * t97 + (t478 + t518) * t211 + t90 * t14 + (Ifges(6,5) * t457 + Ifges(7,5) * t473 + Ifges(6,6) * t458 + Ifges(7,6) * t472 - t339 * t505 + t483 + t516) * t268 + t3 * t77 + t4 * t78 + t71 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t55 * t41 / 0.2e1 + t56 * t40 / 0.2e1 + t47 * t52 + (t177 * t33 - t178 * t32 - t18 * t209 + t19 * t208) * mrSges(6,3) + t12 * t27 + t13 * t28 - t378 * t420 + t417 * (-Ifges(3,6) * t360 + t293) / 0.2e1 + t74 * (t267 * mrSges(5,1) + mrSges(5,3) * t411) + t82 * (mrSges(4,2) * t411 - t267 * mrSges(4,3)) + (t246 * t411 + t247 * t412 - t271 * t359 - t272 * t360) * mrSges(3,3) - t417 * t496 + (-t267 * t75 - t411 * t76) * mrSges(5,2) + (-m(3) * t271 + m(4) * t249 - mrSges(3,1) * t417 + mrSges(4,1) * t267 + mrSges(4,2) * t268) * t247 + t326 * (t376 + (t418 + t432) * t318) * t353 + (-m(3) * t257 + m(4) * t219 - t489) * t262 + (-Ifges(4,2) * t267 - Ifges(4,6) * t411) * t511; ((-t416 / 0.2e1 - t140 / 0.2e1 - t143 / 0.2e1 + t215 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t295 + t395 * t238 + t394 * t237 + (t375 / 0.2e1 + (pkin(1) * mrSges(3,1) + t433 / 0.2e1) * t318) * qJD(1) + (-t505 * t322 - t504 * t325) * qJD(2) / 0.2e1 - t480) * t323 + (-t300 / 0.2e1 - t216 / 0.2e1 + ((pkin(1) * mrSges(3,2) + t428 / 0.2e1 - t418 / 0.2e1) * t318 - t376 / 0.2e1) * qJD(1) + t420 - qJD(2) * Ifges(3,5) / 0.2e1 + (t237 * t392 - t238 * t393 + t328) * t322 + (-t237 * t393 + t238 * t356 - t482) * t325) * t326) * t400 + (t263 - t184) * t174 + (Ifges(7,1) * t182 + Ifges(7,4) * t183) * t462 + (Ifges(7,1) * t148 + Ifges(7,4) * t147) * t463 + (Ifges(7,4) * t182 + Ifges(7,2) * t183) * t464 + (Ifges(7,4) * t148 + Ifges(7,2) * t147) * t465 + (-Ifges(7,4) * t253 - Ifges(7,2) * t252) * t472 + (-Ifges(7,5) * t253 - Ifges(7,6) * t252) * t452 + (-Ifges(7,1) * t253 - Ifges(7,4) * t252) * t473 + t29 * (mrSges(7,1) * t252 - mrSges(7,2) * t253) + (-t1 * t252 + t2 * t253 + t405 * t5 - t406 * t6) * mrSges(7,3) + (Ifges(6,5) * t223 + Ifges(6,6) * t222) * t446 + (Ifges(7,5) * t182 + Ifges(7,6) * t183) * t449 + (Ifges(7,5) * t148 + Ifges(7,6) * t147) * t450 + (Ifges(6,1) * t223 + Ifges(6,4) * t222) * t454 + (Ifges(6,4) * t223 + Ifges(6,2) * t222) * t456 + t529 * t78 + (t1 * t101 + t100 * t2 + t264 * t29 + t5 * t529 + t526 * t71 + t530 * t6) * m(7) + t530 * t77 - t496 + t293 + (t413 * t470 + t410 * t471 + (Ifges(6,5) * t413 + Ifges(6,6) * t410) * t445 + (Ifges(6,1) * t413 + Ifges(6,4) * t410) * t453 + (Ifges(6,4) * t413 + Ifges(6,2) * t410) * t455 + t344 * mrSges(5,1) + t343 * mrSges(4,3) + t95 * (-mrSges(6,1) * t410 + mrSges(6,2) * t413) + (m(4) * t343 + m(5) * t344) * pkin(9) + (-t32 * t413 + t33 * t410) * mrSges(6,3) + (-t402 * pkin(9) + t478) * t322 + (-t401 * pkin(9) + t512 + t522) * t325) * qJD(3) + (-t222 * t33 + t223 * t32) * mrSges(6,3) + t297 * t79 + t288 * t121 - t253 * t474 - t252 * t475 + (t183 / 0.2e1 - t147 / 0.2e1) * t40 + (t182 / 0.2e1 - t148 / 0.2e1) * t41 + t264 * t14 - t257 * t256 - t247 * mrSges(3,1) - t95 * (-mrSges(6,1) * t222 + mrSges(6,2) * t223) - t223 * t88 / 0.2e1 - t222 * t87 / 0.2e1 - t181 * t197 - t180 * t198 - t162 * t199 - t163 * t200 + t201 * t97 + t202 * t96 - pkin(2) * t122 + t100 * t27 + t101 * t28 - m(4) * (-t145 * t180 + t146 * t181 + t219 * t260) - m(5) * (t126 * t184 + t131 * t163 + t132 * t162) + (t38 / 0.2e1 + t37 / 0.2e1 + t107 / 0.2e1 + (-Ifges(7,3) / 0.2e1 + t356) * t195 + t66 / 0.2e1 + (t167 - t168) * pkin(9) + t425 / 0.2e1 + t424 / 0.2e1 + t8 / 0.2e1 + t247 * mrSges(4,2) + t516) * t322 + (t67 * t441 + t68 * t442 - t74 * mrSges(5,1) + t82 * mrSges(4,3) + t54 * t351 - t104 / 0.2e1 + t106 / 0.2e1 - t159 * t349 / 0.2e1 - t160 * t350 / 0.2e1 - t247 * mrSges(4,1) + t75 * mrSges(5,2) + t392 * t196 + t334 * t195 + (-t166 + t169) * pkin(9) + (t18 * t317 - t19 * t319) * mrSges(6,3)) * t325 + m(4) * (-pkin(2) * t247 + t315 * t82 - t438 * t83) + m(5) * (t126 * t263 + t288 * t75 - t315 * t74 + t438 * t76) + (mrSges(7,1) * t406 - mrSges(7,2) * t405) * t71 + t527 * t135 + t528 * t136 + (t18 * t201 + t19 * t202 + t297 * t54 + t32 * t528 + t33 * t527 + t525 * t95) * m(6) + t525 * t114 + t526 * t52 + t489 * t260; (t533 + t540 + t482) * t237 + (-t266 / 0.2e1 - t336 / 0.2e1) * t41 + (-t135 * t317 - t136 * t319) * qJD(5) + t54 * (mrSges(6,1) * t317 + mrSges(6,2) * t319) + (Ifges(6,1) * t319 - t431) * t457 + (-Ifges(6,2) * t317 + t430) * t458 + t486 * t474 + (Ifges(7,4) * t486 - Ifges(7,2) * t342) * t472 + (Ifges(7,1) * t486 - Ifges(7,4) * t342) * t473 + (Ifges(6,5) * t319 + Ifges(7,5) * t486 - Ifges(6,6) * t317 - Ifges(7,6) * t342) * t452 + t29 * (mrSges(7,1) * t342 + mrSges(7,2) * t486) + (mrSges(7,1) * t520 - mrSges(7,2) * t403) * t71 + (-pkin(3) * t76 - qJ(4) * t74 - t126 * t172 - t131 * t146 + t132 * t485) * m(5) + t310 * t14 + t494 * t78 + (t1 * t206 + t2 * t205 + t29 * t310 + (qJD(4) - t85) * t71 + t495 * t6 + t494 * t5) * m(7) + t495 * t77 + (-t265 / 0.2e1 - t491 / 0.2e1) * t40 + (Ifges(7,1) * t336 + Ifges(7,4) * t491) * t463 + (Ifges(7,4) * t336 + Ifges(7,2) * t491) * t465 + (Ifges(7,5) * t336 + Ifges(7,6) * t491) * t450 + t68 * t440 + t67 * t442 + t487 - t348 * mrSges(6,3) + t523 * mrSges(7,3) + t507 + (-Ifges(7,1) * t266 - Ifges(7,4) * t265) * t462 + (-Ifges(7,4) * t266 - Ifges(7,2) * t265) * t464 + (-Ifges(7,5) * t266 - Ifges(7,6) * t265) * t449 + (-t166 + t79) * qJ(4) + t338 * t114 - t342 * t475 - t346 * t435 + t205 * t27 + t206 * t28 - t172 * t174 - pkin(3) * t167 - t60 * t135 - t58 * t136 - t85 * t52 + (qJ(4) * t54 - t348 * t435 + (-t317 * t33 - t319 * t32) * qJD(5) - t32 * t58 - t33 * t60 + t517 * t95) * m(6) - t82 * mrSges(4,2) - t74 * mrSges(5,3) + t76 * mrSges(5,2) + (t349 * t456 + t350 * t454 + t95 * t351 + t87 * t441 + t88 * t442 - t334 * t238 + t347 * mrSges(6,3) + (-t356 + t392) * t237 + t328) * t238 + t401 * t146 + t402 * t145 + t391 * qJD(4); t486 * t27 + t342 * t28 - t403 * t78 + t520 * t77 + t391 * t295 + (t135 * t319 - t136 * t317 + t174) * t238 + t346 + t167 + (t295 * t71 - t523) * m(7) + (-t238 * t347 + t295 * t95 + t348) * m(6) + (t126 * t238 - t132 * t295 + t76) * m(5); t111 * t78 - t515 * t77 - t367 * t135 + t188 * t136 + t14 + t79 + (t111 * t5 - t515 * t6 + t29) * m(7) + (t188 * t32 - t33 * t367 + t54) * m(6); -t71 * (mrSges(7,1) * t111 + mrSges(7,2) * t515) + (Ifges(7,1) * t515 - t429) * t463 + t40 * t462 + (Ifges(7,5) * t515 - Ifges(7,6) * t111) * t450 - t5 * t77 + t6 * t78 + (t111 * t6 + t5 * t515) * mrSges(7,3) + t355 + t8 + (-Ifges(7,2) * t111 + t103 + t41) * t465;];
tauc  = t9(:);
