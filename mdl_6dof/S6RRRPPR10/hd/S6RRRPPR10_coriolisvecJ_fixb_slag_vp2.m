% Calculate vector of centrifugal and Coriolis load on the joints for
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:35
% EndTime: 2019-03-09 16:22:30
% DurationCPUTime: 31.63s
% Computational Cost: add. (13602->838), mult. (35440->1146), div. (0->0), fcn. (26955->10), ass. (0->359)
t325 = cos(qJ(2));
t317 = sin(pkin(6));
t398 = qJD(1) * t317;
t378 = t325 * t398;
t295 = -qJD(3) + t378;
t415 = cos(pkin(6));
t367 = t415 * qJD(1);
t360 = pkin(1) * t367;
t322 = sin(qJ(2));
t379 = t322 * t398;
t257 = -pkin(8) * t379 + t325 * t360;
t336 = (pkin(2) * t322 - pkin(9) * t325) * t317;
t258 = qJD(1) * t336;
t321 = sin(qJ(3));
t324 = cos(qJ(3));
t180 = -t321 * t257 + t324 * t258;
t315 = t324 * pkin(4);
t433 = pkin(3) + qJ(5);
t372 = t433 * t322;
t395 = qJD(3) * t324;
t466 = pkin(4) + pkin(9);
t546 = -(t315 * t325 - t372) * t398 + t180 + t466 * t395;
t396 = qJD(3) * t321;
t365 = pkin(3) * t396 - qJD(4) * t321;
t260 = pkin(8) * t378 + t322 * t360;
t363 = t321 * t378;
t383 = pkin(3) * t363 + t260;
t545 = qJD(5) * t324 - t365 + t383 + t295 * (-qJ(4) * t324 + qJ(5) * t321);
t339 = t367 + qJD(2);
t334 = t324 * t339;
t237 = t321 * t379 - t334;
t219 = -pkin(2) * t339 - t257;
t238 = t321 * t339 + t324 * t379;
t329 = -t238 * qJ(4) + t219;
t126 = t237 * pkin(3) + t329;
t316 = sin(pkin(11));
t318 = cos(pkin(11));
t220 = pkin(9) * t339 + t260;
t251 = (-pkin(2) * t325 - pkin(9) * t322 - pkin(1)) * t317;
t230 = qJD(1) * t251;
t145 = t220 * t321 - t324 * t230;
t337 = pkin(4) * t238 + t145;
t517 = qJD(4) + t337;
t89 = t295 * t433 + t517;
t94 = t237 * t433 + t329;
t32 = -t316 * t94 + t318 * t89;
t33 = t316 * t89 + t318 * t94;
t440 = -t295 / 0.2e1;
t442 = t238 / 0.2e1;
t233 = qJD(6) + t238;
t446 = t233 / 0.2e1;
t188 = t237 * t316 - t295 * t318;
t450 = t188 / 0.2e1;
t366 = t318 * t237 + t295 * t316;
t452 = t366 / 0.2e1;
t320 = sin(qJ(6));
t323 = cos(qJ(6));
t111 = t188 * t323 + t320 * t366;
t459 = t111 / 0.2e1;
t515 = -t188 * t320 + t323 * t366;
t461 = t515 / 0.2e1;
t25 = pkin(5) * t238 - pkin(10) * t188 + t32;
t26 = pkin(10) * t366 + t33;
t5 = t25 * t323 - t26 * t320;
t501 = Ifges(6,3) + Ifges(4,1);
t533 = t295 / 0.2e1;
t445 = -t237 / 0.2e1;
t541 = Ifges(4,4) * t445;
t538 = Ifges(5,6) * t445 + t541;
t6 = t25 * t320 + t26 * t323;
t544 = t5 * mrSges(7,1) - t6 * mrSges(7,2) + 0.2e1 * Ifges(7,5) * t459 + 0.2e1 * Ifges(7,6) * t461 + 0.2e1 * Ifges(7,3) * t446 + 0.2e1 * Ifges(6,5) * t450 + 0.2e1 * Ifges(6,6) * t452 + t32 * mrSges(6,1) + t219 * mrSges(4,2) - t33 * mrSges(6,2) - t126 * mrSges(5,3) + Ifges(5,4) * t533 + Ifges(4,5) * t440 + t538 + (Ifges(5,2) + t501) * t442;
t443 = -t238 / 0.2e1;
t444 = t237 / 0.2e1;
t532 = Ifges(5,4) - Ifges(4,5);
t543 = -Ifges(5,2) * t443 - Ifges(5,6) * t444 - t440 * t532 + t442 * t501 + t541 + t544;
t529 = t545 * t316 + t546 * t318;
t528 = t546 * t316 - t545 * t318;
t405 = t321 * t325;
t223 = (t316 * t405 + t318 * t322) * t398;
t362 = t324 * t378;
t411 = t316 * t321;
t540 = -pkin(5) * t362 + t223 * pkin(10) + (pkin(5) * t324 - pkin(10) * t411) * qJD(3) + t529;
t222 = (-t316 * t322 + t318 * t405) * t398;
t435 = pkin(10) * t318;
t539 = -pkin(10) * t222 + t396 * t435 + t528;
t504 = Ifges(5,1) + Ifges(4,3);
t502 = Ifges(5,5) - Ifges(4,6);
t368 = -qJ(4) * t321 - pkin(2);
t273 = -t324 * t433 + t368;
t296 = t466 * t321;
t277 = t318 * t296;
t179 = pkin(5) * t321 + t277 + (pkin(10) * t324 - t273) * t316;
t202 = t318 * t273 + t316 * t296;
t406 = t318 * t324;
t189 = -pkin(10) * t406 + t202;
t100 = t179 * t323 - t189 * t320;
t531 = qJD(6) * t100 + t320 * t540 + t323 * t539;
t101 = t179 * t320 + t189 * t323;
t530 = -qJD(6) * t101 - t320 * t539 + t323 * t540;
t181 = t324 * t257 + t321 * t258;
t162 = -qJ(4) * t379 - t181;
t137 = -pkin(4) * t363 - t162;
t381 = -pkin(5) * t318 - pkin(4);
t527 = t222 * pkin(5) - t137 + (-pkin(9) + t381) * t396;
t526 = -t466 * t396 - t137;
t483 = -qJD(4) - t145;
t131 = pkin(3) * t295 - t483;
t523 = t131 * mrSges(5,1) + t145 * mrSges(4,3);
t394 = qJD(1) * qJD(2);
t370 = t325 * t394;
t358 = t317 * t370;
t196 = qJD(3) * t238 + t321 * t358;
t371 = t322 * t394;
t359 = t317 * t371;
t159 = t196 * t318 - t316 * t359;
t410 = t317 * t322;
t388 = t321 * t410;
t361 = qJD(3) * t388;
t195 = qJD(1) * t361 - qJD(3) * t334 - t324 * t358;
t380 = pkin(1) * t415;
t409 = t317 * t325;
t272 = pkin(8) * t409 + t322 * t380;
t262 = t272 * qJD(2);
t247 = qJD(1) * t262;
t328 = t195 * qJ(4) - t238 * qJD(4) + t247;
t50 = t237 * qJD(5) + t196 * t433 + t328;
t353 = qJD(2) * t372;
t259 = qJD(2) * t336;
t245 = qJD(1) * t259;
t271 = -pkin(8) * t410 + t325 * t380;
t261 = t271 * qJD(2);
t246 = qJD(1) * t261;
t83 = -t220 * t395 - t230 * t396 + t245 * t324 - t321 * t246;
t51 = -pkin(4) * t195 + qJD(5) * t295 - t353 * t398 - t83;
t19 = t316 * t51 + t318 * t50;
t11 = pkin(10) * t159 + t19;
t160 = t196 * t316 + t318 * t359;
t18 = -t316 * t50 + t318 * t51;
t7 = -pkin(5) * t195 - pkin(10) * t160 + t18;
t1 = qJD(6) * t5 + t11 * t323 + t320 * t7;
t2 = -qJD(6) * t6 - t11 * t320 + t323 * t7;
t341 = t323 * t316 + t320 * t318;
t266 = t341 * qJD(6);
t335 = t341 * t238;
t401 = t335 + t266;
t484 = -t316 * t320 + t318 * t323;
t265 = t484 * qJD(6);
t489 = t484 * t238;
t520 = t265 + t489;
t522 = -t1 * t341 - t2 * t484 + t401 * t5 - t520 * t6;
t146 = t324 * t220 + t321 * t230;
t285 = t295 * qJ(4);
t132 = t285 - t146;
t518 = t132 * mrSges(5,1) - t146 * mrSges(4,3);
t354 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t391 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1;
t449 = -t195 / 0.2e1;
t481 = Ifges(5,2) * t449;
t510 = -t196 / 0.2e1;
t511 = -Ifges(5,4) / 0.2e1;
t75 = t196 * pkin(3) + t328;
t76 = -pkin(3) * t359 - t83;
t516 = t76 * mrSges(5,1) + t18 * mrSges(6,1) - t19 * mrSges(6,2) - t83 * mrSges(4,3) - t75 * mrSges(5,3) + Ifges(5,6) * t510 - t391 * t196 + t359 * t511 + t354 + t481;
t514 = -t219 * mrSges(4,1) + t126 * mrSges(5,2) + Ifges(5,5) * t533 + Ifges(4,6) * t440 + (Ifges(4,2) + Ifges(5,3)) * t445 + (Ifges(4,4) + Ifges(5,6)) * t442;
t42 = qJD(6) * t515 + t159 * t320 + t160 * t323;
t471 = t42 / 0.2e1;
t43 = -qJD(6) * t111 + t159 * t323 - t160 * t320;
t470 = t43 / 0.2e1;
t505 = t83 * mrSges(4,1);
t494 = t246 * mrSges(3,2);
t432 = -pkin(10) - t433;
t282 = t432 * t316;
t283 = t432 * t318;
t205 = -t282 * t320 + t283 * t323;
t118 = -pkin(4) * t237 + t146;
t113 = t318 * t118;
t413 = qJ(4) * t237;
t130 = t238 * t433 + t413;
t34 = -pkin(5) * t237 + t113 + (-pkin(10) * t238 - t130) * t316;
t60 = t316 * t118 + t318 * t130;
t45 = t238 * t435 + t60;
t493 = -qJD(5) * t341 + qJD(6) * t205 - t320 * t34 - t323 * t45;
t206 = t282 * t323 + t283 * t320;
t492 = -qJD(5) * t484 - qJD(6) * t206 + t320 * t45 - t323 * t34;
t487 = mrSges(3,1) * t339 - mrSges(4,1) * t237 - mrSges(4,2) * t238 - mrSges(3,3) * t379;
t485 = t195 * t532 + t196 * t502 + t359 * t504;
t393 = Ifges(4,5) / 0.2e1 + t511;
t480 = -t393 * t295 + t523 + t544;
t478 = -t145 * mrSges(4,1) - t146 * mrSges(4,2) + t131 * mrSges(5,2) - t260 * mrSges(3,3) - t132 * mrSges(5,3);
t476 = -Ifges(4,4) * t442 - Ifges(4,2) * t445 + Ifges(5,6) * t443 + Ifges(5,3) * t444 + t440 * t502 - t514;
t473 = Ifges(7,4) * t471 + Ifges(7,2) * t470 + Ifges(7,6) * t449;
t472 = Ifges(7,1) * t471 + Ifges(7,4) * t470 + Ifges(7,5) * t449;
t68 = t160 * Ifges(6,1) + t159 * Ifges(6,4) - t195 * Ifges(6,5);
t469 = t68 / 0.2e1;
t87 = t188 * Ifges(6,4) + Ifges(6,2) * t366 + Ifges(6,6) * t238;
t468 = t87 / 0.2e1;
t88 = t188 * Ifges(6,1) + Ifges(6,4) * t366 + Ifges(6,5) * t238;
t467 = t88 / 0.2e1;
t462 = -t515 / 0.2e1;
t460 = -t111 / 0.2e1;
t455 = t159 / 0.2e1;
t454 = t160 / 0.2e1;
t453 = -t366 / 0.2e1;
t451 = -t188 / 0.2e1;
t447 = -t233 / 0.2e1;
t439 = -t316 / 0.2e1;
t438 = -t318 / 0.2e1;
t38 = Ifges(7,5) * t42;
t37 = Ifges(7,6) * t43;
t436 = pkin(9) * t321;
t314 = t324 * pkin(9);
t397 = qJD(2) * t317;
t376 = t325 * t397;
t211 = -t361 + (qJD(3) * t415 + t376) * t324;
t250 = pkin(9) * t415 + t272;
t99 = -t250 * t395 - t251 * t396 + t324 * t259 - t321 * t261;
t65 = t211 * pkin(4) + (qJD(5) * t325 - t353) * t317 - t99;
t268 = t321 * t415 + t324 * t410;
t210 = qJD(3) * t268 + t321 * t376;
t267 = -t324 * t415 + t388;
t330 = -t211 * qJ(4) - t268 * qJD(4) + t262;
t69 = t267 * qJD(5) + t210 * t433 + t330;
t24 = t316 * t65 + t318 * t69;
t431 = Ifges(3,4) * t322;
t430 = Ifges(3,4) * t325;
t429 = Ifges(6,4) * t316;
t428 = Ifges(6,4) * t318;
t427 = Ifges(7,4) * t111;
t426 = Ifges(3,2) * t322;
t423 = t159 * Ifges(6,6);
t422 = t160 * Ifges(6,5);
t418 = t257 * mrSges(3,3);
t416 = t322 * Ifges(3,1);
t414 = Ifges(3,6) * qJD(2);
t408 = t318 * t321;
t175 = -t321 * t250 + t324 * t251;
t158 = pkin(3) * t409 - t175;
t119 = t268 * pkin(4) + qJ(5) * t409 + t158;
t249 = -pkin(2) * t415 - t271;
t332 = -t268 * qJ(4) + t249;
t127 = t267 * t433 + t332;
t59 = t316 * t119 + t318 * t127;
t147 = t222 * t323 - t223 * t320;
t183 = t266 * t324 + t396 * t484;
t404 = t147 - t183;
t148 = t222 * t320 + t223 * t323;
t252 = t484 * t324;
t182 = -qJD(6) * t252 + t341 * t396;
t403 = t148 - t182;
t197 = mrSges(4,2) * t295 - mrSges(4,3) * t237;
t199 = mrSges(5,1) * t237 + mrSges(5,3) * t295;
t400 = t197 - t199;
t198 = -mrSges(4,1) * t295 - mrSges(4,3) * t238;
t200 = mrSges(5,1) * t238 - mrSges(5,2) * t295;
t399 = t198 - t200;
t176 = t324 * t250 + t321 * t251;
t297 = t314 + t315;
t8 = -Ifges(7,3) * t195 + t37 + t38;
t392 = Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1;
t390 = -Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1;
t114 = -mrSges(6,1) * t366 + mrSges(6,2) * t188;
t52 = -mrSges(7,1) * t515 + mrSges(7,2) * t111;
t389 = t114 - t199 + t52;
t377 = t322 * t397;
t14 = -t43 * mrSges(7,1) + t42 * mrSges(7,2);
t374 = Ifges(3,5) * t415;
t373 = Ifges(3,6) * t415;
t369 = t397 / 0.2e1;
t23 = -t316 * t69 + t318 * t65;
t167 = -t195 * mrSges(5,1) + mrSges(5,2) * t359;
t79 = -t159 * mrSges(6,1) + t160 * mrSges(6,2);
t57 = t318 * t119 - t127 * t316;
t355 = -Ifges(5,2) / 0.2e1 - Ifges(4,1) / 0.2e1 - Ifges(6,3) / 0.2e1;
t352 = qJD(1) * t369;
t350 = mrSges(6,1) * t318 - mrSges(6,2) * t316;
t349 = Ifges(6,1) * t316 + t428;
t348 = Ifges(6,2) * t318 + t429;
t347 = t18 * t318 + t19 * t316;
t346 = t316 * t32 - t318 * t33;
t96 = mrSges(6,2) * t195 + mrSges(6,3) * t159;
t97 = -mrSges(6,1) * t195 - mrSges(6,3) * t160;
t345 = t316 * t96 + t318 * t97;
t209 = t267 * t316 - t318 * t409;
t35 = pkin(5) * t268 - pkin(10) * t209 + t57;
t208 = t267 * t318 + t316 * t409;
t44 = pkin(10) * t208 + t59;
t12 = -t320 * t44 + t323 * t35;
t13 = t320 * t35 + t323 * t44;
t343 = t131 * t324 + t132 * t321;
t342 = t145 * t324 - t146 * t321;
t133 = t208 * t323 - t209 * t320;
t134 = t208 * t320 + t209 * t323;
t338 = t322 * t352;
t157 = qJ(4) * t409 - t176;
t95 = qJD(5) + t118 - t285;
t82 = -t220 * t396 + t230 * t395 + t321 * t245 + t324 * t246;
t98 = -t250 * t396 + t251 * t395 + t321 * t259 + t324 * t261;
t129 = -t267 * pkin(4) - t157;
t74 = -qJ(4) * t359 + t295 * qJD(4) - t82;
t333 = Ifges(6,5) * t439 + Ifges(6,6) * t438 + t391;
t54 = -pkin(4) * t196 - t74;
t84 = -qJ(4) * t377 + qJD(4) * t409 - t98;
t331 = (t373 + (Ifges(3,2) * t325 + t431) * t317) * qJD(1);
t70 = -t210 * pkin(4) - t84;
t327 = -t392 * t295 + t514 - t518;
t309 = pkin(5) * t316 + qJ(4);
t300 = Ifges(3,4) * t378;
t293 = Ifges(3,5) * t358;
t288 = -pkin(3) * t324 + t368;
t264 = pkin(5) * t406 + t297;
t263 = -qJ(4) * t395 + t365;
t256 = -mrSges(3,2) * t339 + mrSges(3,3) * t378;
t253 = t341 * t324;
t216 = Ifges(3,1) * t379 + Ifges(3,5) * t339 + t300;
t215 = t331 + t414;
t201 = -t273 * t316 + t277;
t184 = -qJ(4) * t362 + t383;
t178 = t210 * t316 + t318 * t377;
t177 = t210 * t318 - t316 * t377;
t174 = -mrSges(5,2) * t237 - mrSges(5,3) * t238;
t172 = pkin(3) * t238 + t413;
t169 = -mrSges(4,2) * t359 - mrSges(4,3) * t196;
t168 = mrSges(4,1) * t359 + mrSges(4,3) * t195;
t166 = mrSges(5,1) * t196 - mrSges(5,3) * t359;
t163 = -pkin(3) * t379 - t180;
t155 = t267 * pkin(3) + t332;
t143 = -t295 * Ifges(5,1) - t238 * Ifges(5,4) + t237 * Ifges(5,5);
t140 = t238 * Ifges(4,5) - t237 * Ifges(4,6) - t295 * Ifges(4,3);
t136 = mrSges(6,1) * t238 - mrSges(6,3) * t188;
t135 = -mrSges(6,2) * t238 + mrSges(6,3) * t366;
t122 = mrSges(4,1) * t196 - mrSges(4,2) * t195;
t121 = -mrSges(5,2) * t196 + mrSges(5,3) * t195;
t107 = -t195 * Ifges(4,1) - t196 * Ifges(4,4) + Ifges(4,5) * t359;
t106 = -t195 * Ifges(4,4) - t196 * Ifges(4,2) + Ifges(4,6) * t359;
t104 = Ifges(5,5) * t359 + t195 * Ifges(5,6) + t196 * Ifges(5,3);
t103 = Ifges(7,4) * t515;
t93 = t210 * pkin(3) + t330;
t91 = -pkin(3) * t377 - t99;
t90 = -t208 * pkin(5) + t129;
t85 = t238 * t381 - t145;
t78 = mrSges(7,1) * t233 - mrSges(7,3) * t111;
t77 = -mrSges(7,2) * t233 + mrSges(7,3) * t515;
t71 = -pkin(5) * t366 + t95;
t67 = t160 * Ifges(6,4) + t159 * Ifges(6,2) - t195 * Ifges(6,6);
t66 = -t195 * Ifges(6,3) + t422 + t423;
t58 = -t130 * t316 + t113;
t56 = -qJD(6) * t134 + t177 * t323 - t178 * t320;
t55 = qJD(6) * t133 + t177 * t320 + t178 * t323;
t47 = -t177 * pkin(5) + t70;
t41 = t111 * Ifges(7,1) + t233 * Ifges(7,5) + t103;
t40 = Ifges(7,2) * t515 + t233 * Ifges(7,6) + t427;
t29 = -pkin(5) * t159 + t54;
t28 = mrSges(7,2) * t195 + mrSges(7,3) * t43;
t27 = -mrSges(7,1) * t195 - mrSges(7,3) * t42;
t20 = pkin(10) * t177 + t24;
t17 = pkin(5) * t211 - pkin(10) * t178 + t23;
t4 = -qJD(6) * t13 + t17 * t323 - t20 * t320;
t3 = qJD(6) * t12 + t17 * t320 + t20 * t323;
t9 = [m(3) * (t246 * t272 + t260 * t261) + m(4) * (-t145 * t99 + t146 * t98 + t175 * t83 + t176 * t82) + m(5) * (t126 * t93 + t131 * t91 + t132 * t84 + t155 * t75 + t157 * t74 + t158 * t76) + m(6) * (t129 * t54 + t18 * t57 + t19 * t59 + t23 * t32 + t24 * t33 + t70 * t95) + m(7) * (t1 * t13 + t12 * t2 + t29 * t90 + t3 * t6 + t4 * t5 + t47 * t71) + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t459 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t471 + (-t267 * t75 - t409 * t76) * mrSges(5,2) + (Ifges(6,5) * t178 + Ifges(6,6) * t177) * t442 + t74 * (t267 * mrSges(5,1) + mrSges(5,3) * t409) + t82 * (mrSges(4,2) * t409 - t267 * mrSges(4,3)) + t325 * (t374 + (t416 + t430) * t317) * t352 + (-m(3) * t271 + m(4) * t249 - mrSges(3,1) * t415 + mrSges(4,1) * t267 + mrSges(4,2) * t268) * t247 - (t331 + t215) * t377 / 0.2e1 + (t107 + t66 + t8) * t268 / 0.2e1 + ((-t426 + t430) * t370 + (Ifges(3,1) * t325 - t431) * t371 - 0.2e1 * pkin(1) * (mrSges(3,1) * t322 + mrSges(3,2) * t325) * t394) * t317 ^ 2 + (t246 * t409 + t247 * t410 - t271 * t358 - t272 * t359) * mrSges(3,3) - t415 * t494 + (t1 * t133 - t134 * t2 - t5 * t55 + t56 * t6) * mrSges(7,3) + (Ifges(6,1) * t178 + Ifges(6,4) * t177) * t450 + (Ifges(6,1) * t209 + Ifges(6,4) * t208) * t454 - t376 * t418 + (Ifges(5,4) * t443 + Ifges(4,5) * t442 + Ifges(5,5) * t444 + Ifges(4,6) * t445 + t440 * t504 + t478) * t377 + (t267 * t502 - t409 * t504) * t338 - t409 * t505 + (-Ifges(4,4) * t267 - Ifges(4,5) * t409 + Ifges(6,5) * t209 + Ifges(7,5) * t134 + Ifges(6,6) * t208 + Ifges(7,6) * t133 + (Ifges(7,3) + t501) * t268) * t449 + (t177 * t33 - t178 * t32 - t18 * t209 + t19 * t208) * mrSges(6,3) + t178 * t467 + t177 * t468 + t209 * t469 + t134 * t472 + t133 * t473 + (t523 + t543) * t211 + t415 * (-Ifges(3,6) * t359 + t293) / 0.2e1 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t446 + (Ifges(6,5) * t454 + Ifges(7,5) * t471 + Ifges(6,6) * t455 + Ifges(7,6) * t470 - t338 * t532 + t481 + t516) * t268 - t267 * t106 / 0.2e1 + t267 * t104 / 0.2e1 + t261 * t256 + t249 * t122 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t461 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t470 + t54 * (-mrSges(6,1) * t208 + mrSges(6,2) * t209) + t208 * t67 / 0.2e1 + t98 * t197 + t99 * t198 + t84 * t199 + t91 * t200 + t95 * (-mrSges(6,1) * t177 + mrSges(6,2) * t178) + t93 * t174 + t175 * t168 + t176 * t169 + t157 * t166 + t158 * t167 + t155 * t121 + t29 * (-mrSges(7,1) * t133 + mrSges(7,2) * t134) + t24 * t135 + t23 * t136 + t129 * t79 + (Ifges(6,4) * t178 + Ifges(6,2) * t177) * t452 + (Ifges(6,4) * t209 + Ifges(6,2) * t208) * t455 + t70 * t114 + t59 * t96 + t57 * t97 + t90 * t14 + t3 * t77 + t4 * t78 + t71 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t55 * t41 / 0.2e1 + t56 * t40 / 0.2e1 + t47 * t52 + t12 * t27 + t13 * t28 + t196 * (-Ifges(5,5) * t409 + Ifges(5,3) * t267) / 0.2e1 + (-Ifges(4,2) * t267 - Ifges(4,6) * t409) * t510 + (t339 * (Ifges(3,5) * t325 - Ifges(3,6) * t322) + t325 * t216 + (t143 + t140) * t322) * t369 - t485 * t409 / 0.2e1 + (-m(3) * t257 + m(4) * t219 - t487) * t262 + (t476 + t518) * t210 + t195 * (-Ifges(5,4) * t409 + Ifges(5,6) * t267) / 0.2e1; t526 * t114 + t527 * t52 + t293 + t528 * t135 + (t18 * t201 + t19 * t202 + t297 * t54 + t32 * t529 + t33 * t528 + t526 * t95) * m(6) + t529 * t136 + (t183 / 0.2e1 - t147 / 0.2e1) * t40 + m(4) * (-pkin(2) * t247 + t314 * t82 - t436 * t83) + m(5) * (t126 * t263 + t288 * t75 - t314 * t74 + t436 * t76) + (t54 * t350 - t104 / 0.2e1 + t106 / 0.2e1 + t75 * mrSges(5,2) - t159 * t348 / 0.2e1 - t160 * t349 / 0.2e1 - t247 * mrSges(4,1) + t67 * t438 + t68 * t439 + t82 * mrSges(4,3) - t74 * mrSges(5,1) + t390 * t196 - t333 * t195 + (-t166 + t169) * pkin(9) + (t18 * t316 - t19 * t318) * mrSges(6,3)) * t324 - t494 + (Ifges(7,5) * t182 + Ifges(7,6) * t183) * t446 + (Ifges(7,5) * t148 + Ifges(7,6) * t147) * t447 + (Ifges(6,1) * t223 + Ifges(6,4) * t222) * t451 + (Ifges(6,4) * t223 + Ifges(6,2) * t222) * t453 + (Ifges(7,1) * t182 + Ifges(7,4) * t183) * t459 + (Ifges(6,5) * t223 + Ifges(6,6) * t222) * t443 + (Ifges(7,1) * t148 + Ifges(7,4) * t147) * t460 + (Ifges(7,4) * t182 + Ifges(7,2) * t183) * t461 + (Ifges(7,4) * t148 + Ifges(7,2) * t147) * t462 - t253 * t472 - t252 * t473 + ((m(4) * t342 + m(5) * t343) * pkin(9) + t95 * (-mrSges(6,1) * t408 + mrSges(6,2) * t411) + (Ifges(6,1) * t411 + Ifges(6,4) * t408) * t450 + (Ifges(6,4) * t411 + Ifges(6,2) * t408) * t452 + (Ifges(6,5) * t411 + Ifges(6,6) * t408) * t442 + t342 * mrSges(4,3) + t343 * mrSges(5,1) + t411 * t467 + t408 * t468 + (-t32 * t411 + t33 * t408) * mrSges(6,3) + (-t400 * pkin(9) + t476) * t321 + (-t399 * pkin(9) + t543) * t324) * qJD(3) + ((-t414 / 0.2e1 - t140 / 0.2e1 - t143 / 0.2e1 + t215 / 0.2e1 + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t295 - t393 * t238 + t392 * t237 + (t373 / 0.2e1 + (pkin(1) * mrSges(3,1) + t431 / 0.2e1) * t317) * qJD(1) + (-t321 * t532 - t324 * t502) * qJD(2) / 0.2e1 - t478) * t322 + (-t216 / 0.2e1 - t300 / 0.2e1 + t418 - qJD(2) * Ifges(3,5) / 0.2e1 + (-t374 / 0.2e1 + (t426 / 0.2e1 + pkin(1) * mrSges(3,2) - t416 / 0.2e1) * t317) * qJD(1) + (t237 * t390 + t238 * t391 + t327) * t321 + (t237 * t391 + t238 * t355 - t480) * t324) * t325) * t398 + t297 * t79 + (-t222 * t33 + t223 * t32) * mrSges(6,3) + t288 * t121 + t264 * t14 - t257 * t256 - t247 * mrSges(3,1) - t222 * t87 / 0.2e1 - t95 * (-mrSges(6,1) * t222 + mrSges(6,2) * t223) - t223 * t88 / 0.2e1 - t181 * t197 - t180 * t198 - t162 * t199 - t163 * t200 + t201 * t97 + t202 * t96 - pkin(2) * t122 + t100 * t27 + t101 * t28 - m(4) * (-t145 * t180 + t146 * t181 + t219 * t260) - m(5) * (t126 * t184 + t131 * t163 + t132 * t162) + t487 * t260 + (-Ifges(7,5) * t253 - Ifges(7,6) * t252) * t449 + (-Ifges(7,4) * t253 - Ifges(7,2) * t252) * t470 + (-Ifges(7,1) * t253 - Ifges(7,4) * t252) * t471 + t29 * (mrSges(7,1) * t252 - mrSges(7,2) * t253) + (-t1 * t252 + t2 * t253 + t403 * t5 - t404 * t6) * mrSges(7,3) + t530 * t78 + (t1 * t101 + t100 * t2 + t264 * t29 + t5 * t530 + t527 * t71 + t531 * t6) * m(7) + t531 * t77 + (t247 * mrSges(4,2) + t38 / 0.2e1 + t37 / 0.2e1 + t107 / 0.2e1 + t423 / 0.2e1 + t66 / 0.2e1 + (t167 - t168) * pkin(9) + (-Ifges(7,3) / 0.2e1 + t355) * t195 + t422 / 0.2e1 + t8 / 0.2e1 + t516) * t321 + (mrSges(7,1) * t404 - mrSges(7,2) * t403) * t71 + (t263 - t184) * t174 + (t182 / 0.2e1 - t148 / 0.2e1) * t41; t389 * qJD(4) + t399 * t146 + t400 * t145 + (-pkin(3) * t76 - qJ(4) * t74 - t126 * t172 - t131 * t146 + t132 * t483) * m(5) + (Ifges(7,1) * t335 + Ifges(7,4) * t489) * t460 + (t480 + t538) * t237 + (Ifges(7,4) * t335 + Ifges(7,2) * t489) * t462 + (-t265 / 0.2e1 - t489 / 0.2e1) * t40 + (Ifges(7,5) * t335 + Ifges(7,6) * t489) * t447 - t341 * t473 - t345 * t433 + t337 * t114 + (-t166 + t79) * qJ(4) + (t87 * t438 + t88 * t439 + t348 * t453 + t349 * t451 + t95 * t350 + t333 * t238 + t346 * mrSges(6,3) + (-t355 + t390) * t237 + t327) * t238 + t492 * t78 + (t1 * t206 + t2 * t205 + t29 * t309 + (qJD(4) - t85) * t71 + t493 * t6 + t492 * t5) * m(7) + t493 * t77 + (-t135 * t316 - t136 * t318) * qJD(5) + (Ifges(6,1) * t318 - t429) * t454 + (-Ifges(6,2) * t316 + t428) * t455 - t347 * mrSges(6,3) + t505 + t67 * t439 + t522 * mrSges(7,3) + t54 * (mrSges(6,1) * t316 + mrSges(6,2) * t318) + t318 * t469 + t309 * t14 + t205 * t27 + t206 * t28 - t172 * t174 - pkin(3) * t167 - t60 * t135 - t58 * t136 - t82 * mrSges(4,2) - t85 * t52 - t74 * mrSges(5,3) + t76 * mrSges(5,2) + (Ifges(7,4) * t484 - Ifges(7,2) * t341) * t470 + (Ifges(7,1) * t484 - Ifges(7,4) * t341) * t471 + t29 * (mrSges(7,1) * t341 + mrSges(7,2) * t484) + (Ifges(6,5) * t318 + Ifges(7,5) * t484 - Ifges(6,6) * t316 - Ifges(7,6) * t341) * t449 + t484 * t472 + (mrSges(7,1) * t520 - mrSges(7,2) * t401) * t71 + (-t266 / 0.2e1 - t335 / 0.2e1) * t41 + (-Ifges(7,5) * t266 - Ifges(7,6) * t265) * t446 + (-Ifges(7,1) * t266 - Ifges(7,4) * t265) * t459 + (-Ifges(7,4) * t266 - Ifges(7,2) * t265) * t461 + t485 + (qJ(4) * t54 - t347 * t433 + (-t316 * t33 - t318 * t32) * qJD(5) - t32 * t58 - t33 * t60 + t517 * t95) * m(6); t484 * t27 + t341 * t28 - t401 * t78 + t520 * t77 + t389 * t295 + (t135 * t318 - t136 * t316 + t174) * t238 + t345 + t167 + (t295 * t71 - t522) * m(7) + (-t238 * t346 + t295 * t95 + t347) * m(6) + (t126 * t238 - t132 * t295 + t76) * m(5); t111 * t78 - t515 * t77 - t366 * t135 + t188 * t136 + t14 + t79 + (t111 * t5 - t515 * t6 + t29) * m(7) + (t188 * t32 - t33 * t366 + t54) * m(6); -t71 * (mrSges(7,1) * t111 + mrSges(7,2) * t515) + (Ifges(7,1) * t515 - t427) * t460 + t40 * t459 + (Ifges(7,5) * t515 - Ifges(7,6) * t111) * t447 - t5 * t77 + t6 * t78 + (t111 * t6 + t5 * t515) * mrSges(7,3) + t354 + t8 + (-Ifges(7,2) * t111 + t103 + t41) * t462;];
tauc  = t9(:);
