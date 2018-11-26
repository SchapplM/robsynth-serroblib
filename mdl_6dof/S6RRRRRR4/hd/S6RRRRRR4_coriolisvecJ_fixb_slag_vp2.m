% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:39:28
% EndTime: 2018-11-23 18:39:55
% DurationCPUTime: 27.17s
% Computational Cost: add. (34112->851), mult. (82478->1194), div. (0->0), fcn. (60474->10), ass. (0->384)
t377 = sin(qJ(2));
t382 = cos(qJ(2));
t403 = pkin(2) * t377 - pkin(8) * t382;
t334 = t403 * qJD(1);
t381 = cos(qJ(3));
t376 = sin(qJ(3));
t434 = qJD(1) * t377;
t415 = t376 * t434;
t287 = pkin(7) * t415 + t381 * t334;
t440 = t381 * t382;
t392 = pkin(3) * t377 - pkin(9) * t440;
t526 = -pkin(9) - pkin(8);
t416 = qJD(3) * t526;
t629 = -qJD(1) * t392 + t381 * t416 - t287;
t313 = t376 * t334;
t441 = t377 * t381;
t442 = t376 * t382;
t628 = t313 + (-pkin(7) * t441 - pkin(9) * t442) * qJD(1) - t376 * t416;
t431 = qJD(2) * t381;
t327 = -t415 + t431;
t414 = t381 * t434;
t328 = qJD(2) * t376 + t414;
t375 = sin(qJ(4));
t380 = cos(qJ(4));
t268 = t327 * t375 + t328 * t380;
t374 = sin(qJ(5));
t379 = cos(qJ(5));
t406 = t380 * t327 - t328 * t375;
t197 = t268 * t379 + t374 * t406;
t373 = sin(qJ(6));
t378 = cos(qJ(6));
t564 = -t268 * t374 + t379 * t406;
t599 = -t197 * t373 + t378 * t564;
t110 = Ifges(7,4) * t599;
t119 = t197 * t378 + t373 * t564;
t345 = -qJD(2) * pkin(2) + pkin(7) * t434;
t291 = -pkin(3) * t327 + t345;
t217 = -pkin(4) * t406 + t291;
t133 = -pkin(5) * t564 + t217;
t558 = pkin(11) * t564;
t340 = -pkin(2) * t382 - t377 * pkin(8) - pkin(1);
t319 = t340 * qJD(1);
t433 = qJD(1) * t382;
t369 = pkin(7) * t433;
t346 = qJD(2) * pkin(8) + t369;
t273 = t381 * t319 - t346 * t376;
t230 = -pkin(9) * t328 + t273;
t361 = qJD(3) - t433;
t219 = pkin(3) * t361 + t230;
t274 = t319 * t376 + t346 * t381;
t231 = pkin(9) * t327 + t274;
t223 = t375 * t231;
t143 = t380 * t219 - t223;
t581 = pkin(10) * t268;
t125 = t143 - t581;
t350 = qJD(4) + t361;
t120 = pkin(4) * t350 + t125;
t225 = t380 * t231;
t144 = t219 * t375 + t225;
t559 = pkin(10) * t406;
t126 = t144 + t559;
t124 = t379 * t126;
t64 = t120 * t374 + t124;
t53 = t64 + t558;
t457 = t373 * t53;
t341 = qJD(5) + t350;
t580 = pkin(11) * t197;
t122 = t374 * t126;
t63 = t379 * t120 - t122;
t52 = t63 - t580;
t49 = pkin(5) * t341 + t52;
t18 = t378 * t49 - t457;
t456 = t378 * t53;
t19 = t373 * t49 + t456;
t421 = qJD(2) * qJD(3);
t429 = qJD(3) * t376;
t430 = qJD(2) * t382;
t285 = t381 * t421 + (-t377 * t429 + t381 * t430) * qJD(1);
t428 = qJD(3) * t381;
t569 = t376 * t430 + t377 * t428;
t286 = -qJD(1) * t569 - t376 * t421;
t170 = qJD(4) * t406 + t285 * t380 + t286 * t375;
t171 = -qJD(4) * t268 - t285 * t375 + t286 * t380;
t85 = qJD(5) * t564 + t170 * t379 + t171 * t374;
t86 = -qJD(5) * t197 - t170 * t374 + t171 * t379;
t32 = qJD(6) * t599 + t373 * t86 + t378 * t85;
t33 = -qJD(6) * t119 - t373 * t85 + t378 * t86;
t432 = qJD(2) * t377;
t409 = qJD(1) * t432;
t420 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t409;
t475 = Ifges(7,4) * t119;
t339 = qJD(6) + t341;
t495 = -t339 / 0.2e1;
t517 = t119 / 0.2e1;
t518 = -t119 / 0.2e1;
t520 = -t599 / 0.2e1;
t424 = qJD(5) * t379;
t425 = qJD(5) * t374;
t337 = t403 * qJD(2);
t320 = qJD(1) * t337;
t405 = pkin(7) * t409;
t209 = -qJD(3) * t274 + t381 * t320 + t376 * t405;
t159 = pkin(3) * t409 - pkin(9) * t285 + t209;
t208 = t319 * t428 + t376 * t320 - t346 * t429 - t381 * t405;
t176 = pkin(9) * t286 + t208;
t66 = -qJD(4) * t144 + t380 * t159 - t176 * t375;
t45 = pkin(4) * t409 - pkin(10) * t170 + t66;
t426 = qJD(4) * t380;
t427 = qJD(4) * t375;
t65 = t375 * t159 + t380 * t176 + t219 * t426 - t231 * t427;
t48 = pkin(10) * t171 + t65;
t12 = t120 * t424 - t126 * t425 + t374 * t45 + t379 * t48;
t10 = pkin(11) * t86 + t12;
t13 = -qJD(5) * t64 - t374 * t48 + t379 * t45;
t9 = pkin(5) * t409 - pkin(11) * t85 + t13;
t2 = qJD(6) * t18 + t10 * t378 + t373 * t9;
t3 = -qJD(6) * t19 - t10 * t373 + t378 * t9;
t557 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t59 = Ifges(7,2) * t599 + Ifges(7,6) * t339 + t475;
t60 = Ifges(7,1) * t119 + Ifges(7,5) * t339 + t110;
t627 = t420 + t557 + t59 * t517 + (Ifges(7,1) * t599 - t475) * t518 + (Ifges(7,5) * t599 - Ifges(7,6) * t119) * t495 + (t119 * t19 + t18 * t599) * mrSges(7,3) - t133 * (mrSges(7,1) * t119 + mrSges(7,2) * t599) + (-Ifges(7,2) * t119 + t110 + t60) * t520;
t393 = t375 * t376 - t380 * t381;
t533 = qJD(3) + qJD(4);
t277 = t533 * t393;
t389 = t393 * t382;
t296 = qJD(1) * t389;
t626 = t277 - t296;
t330 = t375 * t381 + t376 * t380;
t278 = t533 * t330;
t390 = t330 * t382;
t295 = qJD(1) * t390;
t610 = t278 - t295;
t476 = Ifges(6,4) * t197;
t107 = Ifges(6,2) * t564 + Ifges(6,6) * t341 + t476;
t191 = Ifges(6,4) * t564;
t108 = Ifges(6,1) * t197 + Ifges(6,5) * t341 + t191;
t419 = Ifges(6,5) * t85 + Ifges(6,6) * t86 + Ifges(6,3) * t409;
t493 = -t341 / 0.2e1;
t509 = t197 / 0.2e1;
t510 = -t197 / 0.2e1;
t512 = -t564 / 0.2e1;
t545 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t625 = (Ifges(6,5) * t564 - Ifges(6,6) * t197) * t493 + (t197 * t64 + t564 * t63) * mrSges(6,3) - t217 * (mrSges(6,1) * t197 + mrSges(6,2) * t564) + (Ifges(6,1) * t564 - t476) * t510 + t107 * t509 + t419 + t545 + (-Ifges(6,2) * t197 + t108 + t191) * t512 + t627;
t347 = t526 * t376;
t348 = t526 * t381;
t284 = t375 * t347 - t380 * t348;
t572 = -qJD(4) * t284 + t628 * t375 + t380 * t629;
t571 = t347 * t426 + t348 * t427 + t375 * t629 - t628 * t380;
t261 = Ifges(5,4) * t406;
t186 = t268 * Ifges(5,1) + t350 * Ifges(5,5) + t261;
t418 = Ifges(5,5) * t170 + Ifges(5,6) * t171 + Ifges(5,3) * t409;
t477 = Ifges(5,4) * t268;
t491 = -t350 / 0.2e1;
t504 = -t268 / 0.2e1;
t506 = -t406 / 0.2e1;
t542 = t66 * mrSges(5,1) - t65 * mrSges(5,2);
t623 = t418 + t542 + (Ifges(5,5) * t406 - Ifges(5,6) * t268) * t491 + (t143 * t406 + t144 * t268) * mrSges(5,3) + (-Ifges(5,2) * t268 + t186 + t261) * t506 - t291 * (mrSges(5,1) * t268 + mrSges(5,2) * t406) + (Ifges(5,1) * t406 - t477) * t504 + t625;
t621 = -pkin(4) * t434 + pkin(10) * t626 + t572;
t620 = pkin(10) * t610 - t571;
t269 = -t330 * t374 - t379 * t393;
t149 = qJD(5) * t269 - t277 * t379 - t278 * t374;
t227 = -t295 * t374 - t296 * t379;
t439 = t149 - t227;
t270 = t330 * t379 - t374 * t393;
t150 = -qJD(5) * t270 + t277 * t374 - t278 * t379;
t226 = -t295 * t379 + t296 * t374;
t438 = t150 - t226;
t283 = t380 * t347 + t348 * t375;
t244 = -pkin(10) * t330 + t283;
t245 = -pkin(10) * t393 + t284;
t181 = t374 * t244 + t379 * t245;
t577 = -qJD(5) * t181 + t374 * t620 + t379 * t621;
t576 = t244 * t424 - t245 * t425 + t374 * t621 - t379 * t620;
t606 = pkin(5) * t197;
t605 = pkin(11) * t438 + t576;
t604 = -pkin(5) * t434 - pkin(11) * t439 + t577;
t154 = -t230 * t375 - t225;
t131 = t154 - t559;
t155 = t380 * t230 - t223;
t132 = t155 - t581;
t365 = pkin(3) * t380 + pkin(4);
t444 = t375 * t379;
t541 = -t379 * t131 + t132 * t374 - t365 * t425 + (-t375 * t424 + (-t374 * t380 - t444) * qJD(4)) * pkin(3);
t446 = t374 * t375;
t540 = -t374 * t131 - t379 * t132 + t365 * t424 + (-t375 * t425 + (t379 * t380 - t446) * qJD(4)) * pkin(3);
t485 = pkin(3) * t376;
t322 = t433 * t485 + t369;
t601 = pkin(3) * t429 - t322;
t411 = Ifges(3,5) * qJD(2) / 0.2e1;
t589 = t540 + t580;
t588 = t558 + t541;
t570 = pkin(4) * t610 + t601;
t582 = pkin(4) * t268;
t180 = t379 * t244 - t245 * t374;
t138 = -pkin(11) * t270 + t180;
t139 = pkin(11) * t269 + t181;
t90 = t138 * t373 + t139 * t378;
t579 = -qJD(6) * t90 - t373 * t605 + t378 * t604;
t89 = t138 * t378 - t139 * t373;
t578 = qJD(6) * t89 + t373 * t604 + t378 * t605;
t364 = pkin(4) * t379 + pkin(5);
t422 = qJD(6) * t378;
t423 = qJD(6) * t373;
t445 = t374 * t378;
t67 = -t125 * t374 - t124;
t54 = t67 - t558;
t68 = t379 * t125 - t122;
t55 = t68 - t580;
t575 = t373 * t55 - t378 * t54 - t364 * t423 + (-t374 * t422 + (-t373 * t379 - t445) * qJD(5)) * pkin(4);
t447 = t373 * t374;
t574 = -t373 * t54 - t378 * t55 + t364 * t422 + (-t374 * t423 + (t378 * t379 - t447) * qJD(5)) * pkin(4);
t573 = -pkin(5) * t438 + t570;
t367 = Ifges(3,4) * t433;
t463 = t328 * Ifges(4,4);
t247 = t327 * Ifges(4,2) + t361 * Ifges(4,6) + t463;
t323 = Ifges(4,4) * t327;
t248 = t328 * Ifges(4,1) + t361 * Ifges(4,5) + t323;
t394 = t273 * t381 + t274 * t376;
t478 = Ifges(4,4) * t381;
t398 = -Ifges(4,2) * t376 + t478;
t479 = Ifges(4,4) * t376;
t400 = Ifges(4,1) * t381 - t479;
t401 = mrSges(4,1) * t376 + mrSges(4,2) * t381;
t473 = Ifges(4,6) * t376;
t474 = Ifges(4,5) * t381;
t487 = t381 / 0.2e1;
t488 = -t376 / 0.2e1;
t496 = t328 / 0.2e1;
t384 = -t394 * mrSges(4,3) + t345 * t401 + t327 * t398 / 0.2e1 + t400 * t496 + t361 * (-t473 + t474) / 0.2e1 + t247 * t488 + t248 * t487;
t568 = t384 + Ifges(3,1) * t434 / 0.2e1 + t367 / 0.2e1 + t411;
t185 = Ifges(5,2) * t406 + t350 * Ifges(5,6) + t477;
t560 = t185 / 0.2e1;
t410 = -Ifges(3,6) * qJD(2) / 0.2e1;
t306 = -pkin(3) * t446 + t379 * t365;
t304 = pkin(5) + t306;
t308 = pkin(3) * t444 + t365 * t374;
t241 = t304 * t373 + t308 * t378;
t544 = -qJD(6) * t241 - t373 * t589 + t378 * t588;
t240 = t304 * t378 - t308 * t373;
t543 = qJD(6) * t240 + t373 * t588 + t378 * t589;
t303 = t393 * t377;
t326 = t381 * t340;
t484 = pkin(7) * t376;
t272 = -pkin(9) * t441 + t326 + (-pkin(3) - t484) * t382;
t363 = pkin(7) * t440;
t294 = t376 * t340 + t363;
t443 = t376 * t377;
t280 = -pkin(9) * t443 + t294;
t204 = t380 * t272 - t375 * t280;
t174 = -pkin(4) * t382 + t303 * pkin(10) + t204;
t205 = t375 * t272 + t380 * t280;
t302 = t330 * t377;
t182 = -pkin(10) * t302 + t205;
t101 = t374 * t174 + t379 * t182;
t536 = Ifges(4,5) * t285 + Ifges(4,6) * t286;
t535 = -t209 * mrSges(4,1) + t208 * mrSges(4,2);
t534 = qJD(1) * pkin(1) * mrSges(3,2);
t532 = t32 / 0.2e1;
t531 = t33 / 0.2e1;
t528 = t85 / 0.2e1;
t527 = t86 / 0.2e1;
t525 = pkin(1) * mrSges(3,1);
t519 = t599 / 0.2e1;
t235 = -t302 * t379 + t303 * t374;
t236 = -t302 * t374 - t303 * t379;
t165 = t235 * t378 - t236 * t373;
t516 = t165 / 0.2e1;
t166 = t235 * t373 + t236 * t378;
t515 = t166 / 0.2e1;
t514 = t170 / 0.2e1;
t513 = t171 / 0.2e1;
t511 = t564 / 0.2e1;
t508 = t235 / 0.2e1;
t507 = t236 / 0.2e1;
t505 = t406 / 0.2e1;
t503 = t268 / 0.2e1;
t502 = t285 / 0.2e1;
t501 = t286 / 0.2e1;
t500 = -t302 / 0.2e1;
t499 = -t303 / 0.2e1;
t498 = -t327 / 0.2e1;
t497 = -t328 / 0.2e1;
t494 = t339 / 0.2e1;
t492 = t341 / 0.2e1;
t490 = t350 / 0.2e1;
t489 = -t361 / 0.2e1;
t480 = Ifges(3,4) * t377;
t146 = t226 * t378 - t227 * t373;
t199 = t269 * t373 + t270 * t378;
t70 = -qJD(6) * t199 - t149 * t373 + t150 * t378;
t455 = t146 - t70;
t147 = t226 * t373 + t227 * t378;
t198 = t269 * t378 - t270 * t373;
t69 = qJD(6) * t198 + t149 * t378 + t150 * t373;
t454 = t147 - t69;
t451 = qJD(2) * mrSges(3,2);
t435 = t381 * t337 + t432 * t484;
t338 = pkin(3) * t443 + t377 * pkin(7);
t417 = Ifges(4,3) * t409 + t536;
t292 = pkin(3) * t569 + pkin(7) * t430;
t366 = -pkin(3) * t381 - pkin(2);
t264 = -pkin(3) * t286 + qJD(2) * t369;
t100 = t379 * t174 - t374 * t182;
t404 = m(4) * t345 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t327 + mrSges(4,2) * t328 + mrSges(3,3) * t434;
t275 = pkin(4) * t302 + t338;
t222 = pkin(3) * t328 + t582;
t402 = mrSges(4,1) * t381 - mrSges(4,2) * t376;
t399 = Ifges(4,1) * t376 + t478;
t397 = Ifges(4,2) * t381 + t479;
t396 = Ifges(4,5) * t376 + Ifges(4,6) * t381;
t91 = -pkin(5) * t382 - t236 * pkin(11) + t100;
t92 = pkin(11) * t235 + t101;
t39 = -t373 * t92 + t378 * t91;
t40 = t373 * t91 + t378 * t92;
t395 = t208 * t381 - t209 * t376;
t297 = pkin(4) * t393 + t366;
t211 = -qJD(2) * t390 + t303 * t533;
t183 = -pkin(4) * t211 + t292;
t134 = -pkin(4) * t171 + t264;
t210 = -qJD(2) * t389 - t278 * t377;
t203 = t392 * qJD(2) + (-t363 + (pkin(9) * t377 - t340) * t376) * qJD(3) + t435;
t228 = t376 * t337 + t340 * t428 + (-t377 * t431 - t382 * t429) * pkin(7);
t207 = -pkin(9) * t569 + t228;
t96 = -qJD(4) * t205 + t380 * t203 - t207 * t375;
t75 = pkin(4) * t432 - pkin(10) * t210 + t96;
t95 = t375 * t203 + t380 * t207 + t272 * t426 - t280 * t427;
t81 = pkin(10) * t211 + t95;
t22 = t174 * t424 - t182 * t425 + t374 * t75 + t379 * t81;
t23 = -qJD(5) * t101 - t374 * t81 + t379 * t75;
t383 = t410 - (t382 * Ifges(3,2) + t480) * qJD(1) / 0.2e1 + t328 * Ifges(4,5) + t327 * Ifges(4,6) + t341 * Ifges(6,3) + t339 * Ifges(7,3) + t361 * Ifges(4,3) + t350 * Ifges(5,3) + t197 * Ifges(6,5) + t268 * Ifges(5,5) + t406 * Ifges(5,6) + t273 * mrSges(4,1) - t274 * mrSges(4,2) + t63 * mrSges(6,1) - t64 * mrSges(6,2) + t18 * mrSges(7,1) - t19 * mrSges(7,2) + t143 * mrSges(5,1) - t144 * mrSges(5,2) + t564 * Ifges(6,6) + t119 * Ifges(7,5) + t599 * Ifges(7,6);
t343 = mrSges(3,3) * t433 - t451;
t307 = pkin(4) * t445 + t364 * t373;
t305 = -pkin(4) * t447 + t364 * t378;
t293 = -pkin(7) * t442 + t326;
t290 = mrSges(4,1) * t361 - mrSges(4,3) * t328;
t289 = -mrSges(4,2) * t361 + mrSges(4,3) * t327;
t288 = -pkin(7) * t414 + t313;
t263 = -mrSges(4,2) * t409 + mrSges(4,3) * t286;
t262 = mrSges(4,1) * t409 - mrSges(4,3) * t285;
t233 = mrSges(5,1) * t350 - mrSges(5,3) * t268;
t232 = -mrSges(5,2) * t350 + mrSges(5,3) * t406;
t229 = -qJD(3) * t294 + t435;
t221 = -mrSges(4,1) * t286 + mrSges(4,2) * t285;
t220 = -pkin(5) * t269 + t297;
t213 = t285 * Ifges(4,1) + t286 * Ifges(4,4) + Ifges(4,5) * t409;
t212 = t285 * Ifges(4,4) + t286 * Ifges(4,2) + Ifges(4,6) * t409;
t202 = -mrSges(5,1) * t406 + mrSges(5,2) * t268;
t190 = -pkin(5) * t235 + t275;
t179 = mrSges(6,1) * t341 - mrSges(6,3) * t197;
t178 = -mrSges(6,2) * t341 + mrSges(6,3) * t564;
t152 = -mrSges(5,2) * t409 + mrSges(5,3) * t171;
t151 = mrSges(5,1) * t409 - mrSges(5,3) * t170;
t140 = t582 + t606;
t135 = t222 + t606;
t121 = -mrSges(6,1) * t564 + mrSges(6,2) * t197;
t105 = mrSges(7,1) * t339 - mrSges(7,3) * t119;
t104 = -mrSges(7,2) * t339 + mrSges(7,3) * t599;
t103 = -qJD(5) * t236 - t210 * t374 + t211 * t379;
t102 = qJD(5) * t235 + t210 * t379 + t211 * t374;
t99 = -mrSges(5,1) * t171 + mrSges(5,2) * t170;
t98 = Ifges(5,1) * t170 + Ifges(5,4) * t171 + Ifges(5,5) * t409;
t97 = t170 * Ifges(5,4) + t171 * Ifges(5,2) + Ifges(5,6) * t409;
t88 = -pkin(5) * t103 + t183;
t80 = -mrSges(6,2) * t409 + mrSges(6,3) * t86;
t79 = mrSges(6,1) * t409 - mrSges(6,3) * t85;
t61 = -mrSges(7,1) * t599 + mrSges(7,2) * t119;
t47 = -pkin(5) * t86 + t134;
t42 = -qJD(6) * t166 - t102 * t373 + t103 * t378;
t41 = qJD(6) * t165 + t102 * t378 + t103 * t373;
t38 = -mrSges(6,1) * t86 + mrSges(6,2) * t85;
t35 = t85 * Ifges(6,1) + t86 * Ifges(6,4) + Ifges(6,5) * t409;
t34 = Ifges(6,4) * t85 + Ifges(6,2) * t86 + Ifges(6,6) * t409;
t29 = -mrSges(7,2) * t409 + mrSges(7,3) * t33;
t28 = mrSges(7,1) * t409 - mrSges(7,3) * t32;
t21 = t378 * t52 - t457;
t20 = -t373 * t52 - t456;
t17 = pkin(11) * t103 + t22;
t16 = pkin(5) * t432 - pkin(11) * t102 + t23;
t8 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t7 = t32 * Ifges(7,1) + t33 * Ifges(7,4) + Ifges(7,5) * t409;
t6 = t32 * Ifges(7,4) + t33 * Ifges(7,2) + Ifges(7,6) * t409;
t5 = -qJD(6) * t40 + t16 * t378 - t17 * t373;
t4 = qJD(6) * t39 + t16 * t373 + t17 * t378;
t1 = [(Ifges(7,1) * t166 + Ifges(7,4) * t165) * t532 + (Ifges(7,4) * t166 + Ifges(7,2) * t165) * t531 + (-Ifges(5,5) * t514 - Ifges(6,5) * t528 - Ifges(7,5) * t532 - Ifges(5,6) * t513 - Ifges(6,6) * t527 - Ifges(7,6) * t531 + (0.3e1 / 0.2e1 * Ifges(3,4) * t430 + (-Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1 - Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 + (m(4) * pkin(7) + t401) * pkin(7)) * t432) * qJD(1) + t535 - t542 - t545 - t557) * t382 + m(5) * (t143 * t96 + t144 * t95 + t204 * t66 + t205 * t65 + t264 * t338 + t291 * t292) + m(6) * (t100 * t13 + t101 * t12 + t134 * t275 + t183 * t217 + t22 * t64 + t23 * t63) + m(7) * (t133 * t88 + t18 * t5 + t19 * t4 + t190 * t47 + t2 * t40 + t3 * t39) + (-t143 * t210 + t144 * t211 - t302 * t65 + t303 * t66) * mrSges(5,3) + (-Ifges(5,1) * t303 - Ifges(5,4) * t302) * t514 + t100 * t79 + t101 * t80 + (Ifges(6,1) * t236 + Ifges(6,4) * t235) * t528 + t133 * (-mrSges(7,1) * t42 + mrSges(7,2) * t41) + (t165 * t2 - t166 * t3 - t18 * t41 + t19 * t42) * mrSges(7,3) + (t400 * t502 + t398 * t501 + t213 * t487 + t212 * t488 + pkin(7) * t221 + (-t208 * t376 - t209 * t381) * mrSges(4,3) + (t397 * t498 + t399 * t497 + t345 * t402 + t396 * t489 + t248 * t488 - t381 * t247 / 0.2e1 + (t273 * t376 - t274 * t381) * mrSges(4,3)) * qJD(3) + (t383 - pkin(7) * t343 + ((-0.3e1 / 0.2e1 * Ifges(3,4) + t474 / 0.2e1 - t473 / 0.2e1) * t377 + Ifges(7,5) * t515 + Ifges(7,6) * t516 + Ifges(6,5) * t507 + Ifges(6,6) * t508 + Ifges(5,5) * t499 + Ifges(5,6) * t500 - 0.2e1 * t525) * qJD(1) + t410) * qJD(2)) * t377 - (t420 + t419 + t418 + t417 + t536) * t382 / 0.2e1 + t5 * t105 + t103 * t107 / 0.2e1 + t102 * t108 / 0.2e1 + t4 * t104 + (-Ifges(5,4) * t303 - Ifges(5,2) * t302) * t513 + t264 * (mrSges(5,1) * t302 - mrSges(5,2) * t303) + m(4) * (t208 * t294 + t209 * t293 + t274 * t228 + t273 * t229) + t88 * t61 + t42 * t59 / 0.2e1 + t41 * t60 / 0.2e1 + t39 * t28 + t40 * t29 + (Ifges(6,4) * t102 + Ifges(6,2) * t103) * t511 + t7 * t515 + t6 * t516 + (Ifges(7,1) * t41 + Ifges(7,4) * t42) * t517 + (Ifges(7,4) * t41 + Ifges(7,2) * t42) * t519 + t98 * t499 + t97 * t500 + (Ifges(5,1) * t210 + Ifges(5,4) * t211) * t503 + (Ifges(5,4) * t210 + Ifges(5,2) * t211) * t505 + t35 * t507 + t34 * t508 + (Ifges(6,1) * t102 + Ifges(6,4) * t103) * t509 + (Ifges(5,5) * t210 + Ifges(5,6) * t211) * t490 + (Ifges(6,5) * t102 + Ifges(6,6) * t103) * t492 + (Ifges(7,5) * t41 + Ifges(7,6) * t42) * t494 + (Ifges(6,4) * t236 + Ifges(6,2) * t235) * t527 + t47 * (-mrSges(7,1) * t165 + mrSges(7,2) * t166) + t22 * t178 + t23 * t179 + t183 * t121 + t190 * t8 + t211 * t560 + t204 * t151 + t205 * t152 + t210 * t186 / 0.2e1 + (t404 * pkin(7) + t411 - 0.2e1 * t534 + t568) * t430 + t217 * (-mrSges(6,1) * t103 + mrSges(6,2) * t102) + t95 * t232 + t96 * t233 + t134 * (-mrSges(6,1) * t235 + mrSges(6,2) * t236) + t275 * t38 + t228 * t289 + t229 * t290 + t291 * (-mrSges(5,1) * t211 + mrSges(5,2) * t210) + t292 * t202 + t293 * t262 + t294 * t263 + (-t102 * t63 + t103 * t64 + t12 * t235 - t13 * t236) * mrSges(6,3) + t338 * t99; (t143 * t626 - t144 * t610 - t330 * t66 - t393 * t65) * mrSges(5,3) + (mrSges(5,1) * t610 - mrSges(5,2) * t626) * t291 + (t149 / 0.2e1 - t227 / 0.2e1) * t108 - m(4) * (t273 * t287 + t274 * t288) + (t150 / 0.2e1 - t226 / 0.2e1) * t107 + (t70 / 0.2e1 - t146 / 0.2e1) * t59 + t89 * t28 + t90 * t29 + (mrSges(7,1) * t455 - mrSges(7,2) * t454) * t133 + (t18 * t454 - t19 * t455 + t198 * t2 - t199 * t3) * mrSges(7,3) + (-mrSges(6,1) * t438 + mrSges(6,2) * t439) * t217 + (t12 * t269 - t13 * t270 + t438 * t64 - t439 * t63) * mrSges(6,3) + t395 * mrSges(4,3) + (-Ifges(5,1) * t296 - Ifges(5,4) * t295) * t504 + (-Ifges(5,4) * t296 - Ifges(5,2) * t295) * t506 + (-Ifges(5,5) * t296 - Ifges(5,6) * t295) * t491 + (-t277 / 0.2e1 + t296 / 0.2e1) * t186 + (-Ifges(5,1) * t277 - Ifges(5,4) * t278) * t503 + (-Ifges(5,4) * t277 - Ifges(5,2) * t278) * t505 + (-Ifges(5,5) * t277 - Ifges(5,6) * t278) * t490 + (-t278 / 0.2e1 + t295 / 0.2e1) * t185 + (Ifges(5,4) * t330 - Ifges(5,2) * t393) * t513 + (Ifges(5,1) * t330 - Ifges(5,4) * t393) * t514 + t264 * (mrSges(5,1) * t393 + mrSges(5,2) * t330) - t393 * t97 / 0.2e1 + (t202 * t485 + t384) * qJD(3) + ((-t289 * t376 - t290 * t381) * qJD(3) - t262 * t376 + t263 * t381 + m(4) * (-qJD(3) * t394 + t395)) * pkin(8) + (Ifges(6,4) * t270 + Ifges(6,2) * t269) * t527 + (Ifges(6,1) * t270 + Ifges(6,4) * t269) * t528 + (Ifges(7,4) * t199 + Ifges(7,2) * t198) * t531 + (Ifges(7,1) * t199 + Ifges(7,4) * t198) * t532 + (Ifges(6,1) * t227 + Ifges(6,4) * t226) * t510 + (Ifges(6,4) * t149 + Ifges(6,2) * t150) * t511 + (Ifges(6,4) * t227 + Ifges(6,2) * t226) * t512 + (Ifges(7,1) * t69 + Ifges(7,4) * t70) * t517 + (Ifges(7,1) * t147 + Ifges(7,4) * t146) * t518 + (Ifges(7,4) * t69 + Ifges(7,2) * t70) * t519 + (Ifges(7,4) * t147 + Ifges(7,2) * t146) * t520 + t397 * t501 + t399 * t502 + (Ifges(6,1) * t149 + Ifges(6,4) * t150) * t509 + t212 * t487 + (Ifges(6,5) * t149 + Ifges(6,6) * t150) * t492 + (Ifges(6,5) * t227 + Ifges(6,6) * t226) * t493 + (Ifges(7,5) * t69 + Ifges(7,6) * t70) * t494 + (Ifges(7,5) * t147 + Ifges(7,6) * t146) * t495 + (t572 * t143 + t571 * t144 + t264 * t366 + t283 * t66 + t284 * t65 + t291 * t601) * m(5) + t180 * t79 + t181 * t80 + t198 * t6 / 0.2e1 + t199 * t7 / 0.2e1 + t47 * (-mrSges(7,1) * t198 + mrSges(7,2) * t199) + ((t411 - t367 / 0.2e1 + t534 + ((-m(4) * pkin(2) - mrSges(3,1) - t402) * qJD(2) - t404) * pkin(7) - t568) * t382 + (-t383 + (t525 + t480 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t382) * qJD(1) + (t343 + t451) * pkin(7) + t410) * t377 + (Ifges(5,5) * t330 + Ifges(6,5) * t270 + Ifges(7,5) * t199 - Ifges(5,6) * t393 + Ifges(6,6) * t269 + Ifges(7,6) * t198 + t396) * t432 / 0.2e1) * qJD(1) + t220 * t8 - pkin(2) * t221 + t570 * t121 + t571 * t232 + t572 * t233 + t573 * t61 + t269 * t34 / 0.2e1 + t270 * t35 / 0.2e1 + t134 * (-mrSges(6,1) * t269 + mrSges(6,2) * t270) + t283 * t151 + t284 * t152 - t288 * t289 - t287 * t290 + t297 * t38 + t576 * t178 + t577 * t179 + (t12 * t181 + t13 * t180 + t134 * t297 + t217 * t570 + t576 * t64 + t577 * t63) * m(6) + t578 * t104 + t579 * t105 + (t133 * t573 + t18 * t579 + t19 * t578 + t2 * t90 + t220 * t47 + t3 * t89) * m(7) - t322 * t202 + (t69 / 0.2e1 - t147 / 0.2e1) * t60 + t330 * t98 / 0.2e1 + t366 * t99 + t376 * t213 / 0.2e1; (t12 * t308 + t13 * t306 - t217 * t222 + t540 * t64 + t541 * t63) * m(6) + t541 * t179 + t543 * t104 + (-t133 * t135 + t18 * t544 + t19 * t543 + t2 * t241 + t240 * t3) * m(7) + t544 * t105 - m(5) * (t143 * t154 + t144 * t155) - t535 + (-Ifges(4,2) * t328 + t248 + t323) * t498 + t417 + t540 * t178 - t135 * t61 + t623 + (t151 * t380 + t152 * t375 - t202 * t328 + (t232 * t380 - t233 * t375) * qJD(4) + (-t143 * t427 + t144 * t426 + 0.2e1 * t291 * t497 + t375 * t65 + t380 * t66) * m(5)) * pkin(3) + t268 * t560 + (Ifges(4,5) * t327 - Ifges(4,6) * t328) * t489 + t247 * t496 + (Ifges(4,1) * t327 - t463) * t497 + (t273 * t327 + t274 * t328) * mrSges(4,3) - t222 * t121 - t155 * t232 - t154 * t233 + t240 * t28 + t241 * t29 - t273 * t289 + t274 * t290 + t306 * t79 + t308 * t80 - t345 * (mrSges(4,1) * t328 + mrSges(4,2) * t327); -m(6) * (t63 * t67 + t64 * t68) + (-t268 * t121 + t374 * t80 + t379 * t79 + (t178 * t379 - t179 * t374) * qJD(5) + (t12 * t374 + t13 * t379 - t217 * t268 + t424 * t64 - t425 * t63) * m(6)) * pkin(4) + t185 * t503 - t140 * t61 - t68 * t178 - t67 * t179 - t143 * t232 + t144 * t233 + t574 * t104 + t575 * t105 + (-t133 * t140 + t18 * t575 + t19 * t574 + t2 * t307 + t3 * t305) * m(7) + t305 * t28 + t307 * t29 + t623; -m(7) * (t18 * t20 + t19 * t21) + (-t197 * t61 + t378 * t28 + t373 * t29 + (t104 * t378 - t105 * t373) * qJD(6) + (-t133 * t197 - t18 * t423 + t19 * t422 + t2 * t373 + t3 * t378) * m(7)) * pkin(5) - t21 * t104 - t20 * t105 - t63 * t178 + t64 * t179 + t625; -t18 * t104 + t19 * t105 + t627;];
tauc  = t1(:);
