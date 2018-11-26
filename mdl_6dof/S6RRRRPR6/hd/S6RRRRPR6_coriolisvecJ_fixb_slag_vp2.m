% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:15:11
% EndTime: 2018-11-23 18:15:36
% DurationCPUTime: 25.56s
% Computational Cost: add. (26441->819), mult. (65405->1143), div. (0->0), fcn. (47890->10), ass. (0->371)
t371 = sin(qJ(2));
t375 = cos(qJ(2));
t395 = pkin(2) * t371 - pkin(8) * t375;
t328 = t395 * qJD(1);
t374 = cos(qJ(3));
t370 = sin(qJ(3));
t422 = qJD(1) * t371;
t406 = t370 * t422;
t277 = pkin(7) * t406 + t374 * t328;
t428 = t374 * t375;
t384 = pkin(3) * t371 - pkin(9) * t428;
t507 = -pkin(9) - pkin(8);
t407 = qJD(3) * t507;
t596 = -qJD(1) * t384 + t374 * t407 - t277;
t307 = t370 * t328;
t429 = t371 * t374;
t430 = t370 * t375;
t595 = t307 + (-pkin(7) * t429 - pkin(9) * t430) * qJD(1) - t370 * t407;
t369 = sin(qJ(4));
t373 = cos(qJ(4));
t385 = t369 * t370 - t373 * t374;
t516 = qJD(3) + qJD(4);
t267 = t516 * t385;
t382 = t375 * t385;
t286 = qJD(1) * t382;
t594 = t267 - t286;
t324 = t369 * t374 + t370 * t373;
t268 = t516 * t324;
t381 = t324 * t375;
t285 = qJD(1) * t381;
t581 = t268 - t285;
t419 = qJD(2) * t374;
t321 = -t406 + t419;
t405 = t374 * t422;
t322 = qJD(2) * t370 + t405;
t260 = t321 * t369 + t322 * t373;
t366 = sin(pkin(11));
t367 = cos(pkin(11));
t398 = t373 * t321 - t322 * t369;
t191 = t260 * t367 + t366 * t398;
t368 = sin(qJ(6));
t372 = cos(qJ(6));
t541 = -t260 * t366 + t367 * t398;
t573 = -t191 * t368 + t372 * t541;
t106 = Ifges(7,4) * t573;
t114 = t191 * t372 + t368 * t541;
t339 = -qJD(2) * pkin(2) + pkin(7) * t422;
t281 = -pkin(3) * t321 + t339;
t209 = -pkin(4) * t398 + qJD(5) + t281;
t126 = -pkin(5) * t541 + t209;
t421 = qJD(1) * t375;
t354 = qJD(3) - t421;
t343 = qJD(4) + t354;
t556 = pkin(10) * t191;
t333 = -pkin(2) * t375 - t371 * pkin(8) - pkin(1);
t313 = t333 * qJD(1);
t362 = pkin(7) * t421;
t340 = qJD(2) * pkin(8) + t362;
t263 = t374 * t313 - t340 * t370;
t225 = -pkin(9) * t322 + t263;
t215 = pkin(3) * t354 + t225;
t264 = t313 * t370 + t340 * t374;
t226 = pkin(9) * t321 + t264;
t220 = t369 * t226;
t142 = t373 * t215 - t220;
t550 = qJ(5) * t260;
t121 = t142 - t550;
t115 = pkin(4) * t343 + t121;
t222 = t373 * t226;
t143 = t215 * t369 + t222;
t526 = qJ(5) * t398;
t122 = t143 + t526;
t117 = t366 * t122;
t58 = t367 * t115 - t117;
t44 = pkin(5) * t343 - t556 + t58;
t538 = pkin(10) * t541;
t433 = t367 * t122;
t59 = t366 * t115 + t433;
t46 = t59 + t538;
t16 = -t368 * t46 + t372 * t44;
t17 = t368 * t44 + t372 * t46;
t413 = qJD(2) * qJD(3);
t417 = qJD(3) * t370;
t418 = qJD(2) * t375;
t275 = t374 * t413 + (-t371 * t417 + t374 * t418) * qJD(1);
t416 = qJD(3) * t374;
t545 = t370 * t418 + t371 * t416;
t276 = -qJD(1) * t545 - t370 * t413;
t165 = qJD(4) * t398 + t275 * t373 + t276 * t369;
t166 = -qJD(4) * t260 - t275 * t369 + t276 * t373;
t96 = -t165 * t366 + t166 * t367;
t97 = t165 * t367 + t166 * t366;
t31 = qJD(6) * t573 + t368 * t96 + t372 * t97;
t32 = -qJD(6) * t114 - t368 * t97 + t372 * t96;
t420 = qJD(2) * t371;
t400 = qJD(1) * t420;
t412 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t400;
t457 = Ifges(7,4) * t114;
t334 = qJD(6) + t343;
t476 = -t334 / 0.2e1;
t498 = t114 / 0.2e1;
t499 = -t114 / 0.2e1;
t501 = -t573 / 0.2e1;
t331 = t395 * qJD(2);
t314 = qJD(1) * t331;
t397 = pkin(7) * t400;
t204 = -qJD(3) * t264 + t374 * t314 + t370 * t397;
t157 = pkin(3) * t400 - pkin(9) * t275 + t204;
t203 = t313 * t416 + t370 * t314 - t340 * t417 - t374 * t397;
t173 = pkin(9) * t276 + t203;
t64 = -qJD(4) * t143 + t373 * t157 - t173 * t369;
t39 = pkin(4) * t400 - qJ(5) * t165 - qJD(5) * t260 + t64;
t414 = qJD(4) * t373;
t415 = qJD(4) * t369;
t63 = t369 * t157 + t373 * t173 + t215 * t414 - t226 * t415;
t41 = qJ(5) * t166 + qJD(5) * t398 + t63;
t13 = t366 * t39 + t367 * t41;
t10 = pkin(10) * t96 + t13;
t12 = -t366 * t41 + t367 * t39;
t9 = pkin(5) * t400 - pkin(10) * t97 + t12;
t2 = qJD(6) * t16 + t10 * t372 + t368 * t9;
t3 = -qJD(6) * t17 - t10 * t368 + t372 * t9;
t537 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t54 = Ifges(7,2) * t573 + Ifges(7,6) * t334 + t457;
t55 = Ifges(7,1) * t114 + Ifges(7,5) * t334 + t106;
t593 = (Ifges(7,1) * t573 - t457) * t499 + (Ifges(7,5) * t573 - Ifges(7,6) * t114) * t476 + (t114 * t17 + t16 * t573) * mrSges(7,3) - t126 * (mrSges(7,1) * t114 + mrSges(7,2) * t573) + t54 * t498 + t412 + t537 + (-Ifges(7,2) * t114 + t106 + t55) * t501;
t341 = t507 * t370;
t342 = t507 * t374;
t548 = t341 * t414 + t342 * t415 + t369 * t596 - t595 * t373;
t274 = t369 * t341 - t373 * t342;
t547 = -qJD(4) * t274 + t595 * t369 + t373 * t596;
t551 = Ifges(6,4) * t191;
t104 = Ifges(6,2) * t541 + Ifges(6,6) * t343 + t551;
t534 = Ifges(6,4) * t541;
t105 = Ifges(6,1) * t191 + Ifges(6,5) * t343 + t534;
t251 = Ifges(5,4) * t398;
t181 = Ifges(5,1) * t260 + Ifges(5,5) * t343 + t251;
t458 = Ifges(5,4) * t260;
t474 = -t343 / 0.2e1;
t485 = -t260 / 0.2e1;
t487 = -t398 / 0.2e1;
t490 = t191 / 0.2e1;
t491 = -t191 / 0.2e1;
t493 = -t541 / 0.2e1;
t514 = t64 * mrSges(5,1) + t12 * mrSges(6,1) - t63 * mrSges(5,2) - t13 * mrSges(6,2);
t570 = Ifges(6,3) + Ifges(5,3);
t519 = Ifges(5,5) * t165 + Ifges(6,5) * t97 + Ifges(5,6) * t166 + Ifges(6,6) * t96 + t400 * t570;
t592 = t514 + t519 + (t142 * t398 + t143 * t260) * mrSges(5,3) + (-Ifges(5,2) * t260 + t181 + t251) * t487 - t281 * (mrSges(5,1) * t260 + mrSges(5,2) * t398) + (Ifges(5,1) * t398 - t458) * t485 + (Ifges(5,5) * t398 + Ifges(6,5) * t541 - Ifges(5,6) * t260 - Ifges(6,6) * t191) * t474 + (t191 * t59 + t541 * t58) * mrSges(6,3) + t104 * t490 - t209 * (mrSges(6,1) * t191 + mrSges(6,2) * t541) + (-Ifges(6,2) * t191 + t105 + t534) * t493 + (Ifges(6,1) * t541 - t551) * t491 + t593;
t590 = -pkin(4) * t422 + qJ(5) * t594 - qJD(5) * t324 + t547;
t589 = -qJ(5) * t581 - qJD(5) * t385 + t548;
t197 = t267 * t366 - t268 * t367;
t218 = -t285 * t367 + t286 * t366;
t427 = t197 - t218;
t198 = -t267 * t367 - t268 * t366;
t219 = -t285 * t366 - t286 * t367;
t426 = t198 - t219;
t553 = -t366 * t589 + t367 * t590;
t552 = t366 * t590 + t367 * t589;
t578 = pkin(5) * t191;
t577 = -pkin(5) * t422 - pkin(10) * t426 + t553;
t576 = pkin(10) * t427 + t552;
t153 = -t225 * t369 - t222;
t127 = t153 - t526;
t154 = t373 * t225 - t220;
t128 = t154 - t550;
t432 = t367 * t369;
t454 = pkin(3) * qJD(4);
t528 = -t367 * t127 + t128 * t366 + (-t366 * t373 - t432) * t454;
t434 = t366 * t369;
t527 = -t366 * t127 - t367 * t128 + (t367 * t373 - t434) * t454;
t468 = pkin(3) * t370;
t316 = t421 * t468 + t362;
t575 = pkin(3) * t417 - t316;
t402 = Ifges(3,5) * qJD(2) / 0.2e1;
t562 = t538 + t528;
t561 = t527 + t556;
t546 = pkin(4) * t581 + t575;
t467 = pkin(4) * t260;
t273 = t373 * t341 + t342 * t369;
t238 = -qJ(5) * t324 + t273;
t239 = -qJ(5) * t385 + t274;
t169 = t367 * t238 - t239 * t366;
t256 = t324 * t367 - t366 * t385;
t133 = -pkin(10) * t256 + t169;
t170 = t366 * t238 + t367 * t239;
t255 = -t324 * t366 - t367 * t385;
t134 = pkin(10) * t255 + t170;
t76 = t133 * t372 - t134 * t368;
t555 = qJD(6) * t76 + t368 * t577 + t372 * t576;
t77 = t133 * t368 + t134 * t372;
t554 = -qJD(6) * t77 - t368 * t576 + t372 * t577;
t549 = -pkin(5) * t427 + t546;
t360 = Ifges(3,4) * t421;
t444 = t322 * Ifges(4,4);
t241 = t321 * Ifges(4,2) + t354 * Ifges(4,6) + t444;
t317 = Ifges(4,4) * t321;
t242 = t322 * Ifges(4,1) + t354 * Ifges(4,5) + t317;
t386 = t263 * t374 + t264 * t370;
t459 = Ifges(4,4) * t374;
t390 = -Ifges(4,2) * t370 + t459;
t460 = Ifges(4,4) * t370;
t392 = Ifges(4,1) * t374 - t460;
t393 = mrSges(4,1) * t370 + mrSges(4,2) * t374;
t455 = Ifges(4,6) * t370;
t456 = Ifges(4,5) * t374;
t470 = t374 / 0.2e1;
t471 = -t370 / 0.2e1;
t477 = t322 / 0.2e1;
t377 = -t386 * mrSges(4,3) + t339 * t393 + t321 * t390 / 0.2e1 + t392 * t477 + t354 * (-t455 + t456) / 0.2e1 + t241 * t471 + t242 * t470;
t544 = t377 + Ifges(3,1) * t422 / 0.2e1 + t360 / 0.2e1 + t402;
t180 = Ifges(5,2) * t398 + Ifges(5,6) * t343 + t458;
t539 = t180 / 0.2e1;
t473 = t343 / 0.2e1;
t401 = -Ifges(3,6) * qJD(2) / 0.2e1;
t358 = pkin(3) * t373 + pkin(4);
t299 = -pkin(3) * t434 + t367 * t358;
t291 = pkin(5) + t299;
t301 = pkin(3) * t432 + t358 * t366;
t235 = t291 * t368 + t301 * t372;
t532 = -qJD(6) * t235 - t368 * t561 + t372 * t562;
t234 = t291 * t372 - t301 * t368;
t531 = qJD(6) * t234 + t368 * t562 + t372 * t561;
t357 = pkin(4) * t367 + pkin(5);
t466 = pkin(4) * t366;
t300 = t357 * t372 - t368 * t466;
t67 = -t121 * t366 - t433;
t47 = t67 - t538;
t68 = t367 * t121 - t117;
t48 = t68 - t556;
t530 = t300 * qJD(6) - t368 * t47 - t372 * t48;
t302 = t357 * t368 + t372 * t466;
t529 = -t302 * qJD(6) + t368 * t48 - t372 * t47;
t298 = t385 * t371;
t320 = t374 * t333;
t465 = pkin(7) * t370;
t262 = -pkin(9) * t429 + t320 + (-pkin(3) - t465) * t375;
t356 = pkin(7) * t428;
t284 = t370 * t333 + t356;
t431 = t370 * t371;
t270 = -pkin(9) * t431 + t284;
t200 = t369 * t262 + t373 * t270;
t520 = Ifges(4,5) * t275 + Ifges(4,6) * t276;
t518 = -t204 * mrSges(4,1) + t203 * mrSges(4,2);
t517 = pkin(1) * mrSges(3,2) * qJD(1);
t410 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t461 = Ifges(3,4) * t371;
t515 = t410 * t343 + t260 * Ifges(5,5) + t398 * Ifges(5,6) + t322 * Ifges(4,5) + t321 * Ifges(4,6) + t354 * Ifges(4,3) + t334 * Ifges(7,3) + t263 * mrSges(4,1) - t264 * mrSges(4,2) + t142 * mrSges(5,1) - t143 * mrSges(5,2) + t16 * mrSges(7,1) - t17 * mrSges(7,2) - t59 * mrSges(6,2) + t58 * mrSges(6,1) + t573 * Ifges(7,6) + t401 - (t375 * Ifges(3,2) + t461) * qJD(1) / 0.2e1 + t114 * Ifges(7,5) + t191 * Ifges(6,5) + t541 * Ifges(6,6) + t570 * t473;
t513 = t31 / 0.2e1;
t512 = t32 / 0.2e1;
t509 = t96 / 0.2e1;
t508 = t97 / 0.2e1;
t506 = pkin(1) * mrSges(3,1);
t500 = t573 / 0.2e1;
t297 = t324 * t371;
t229 = -t297 * t367 + t298 * t366;
t230 = -t297 * t366 - t298 * t367;
t158 = t229 * t372 - t230 * t368;
t497 = t158 / 0.2e1;
t159 = t229 * t368 + t230 * t372;
t496 = t159 / 0.2e1;
t495 = t165 / 0.2e1;
t494 = t166 / 0.2e1;
t492 = t541 / 0.2e1;
t489 = t229 / 0.2e1;
t488 = t230 / 0.2e1;
t486 = t398 / 0.2e1;
t484 = t260 / 0.2e1;
t483 = t275 / 0.2e1;
t482 = t276 / 0.2e1;
t481 = -t297 / 0.2e1;
t480 = -t298 / 0.2e1;
t479 = -t321 / 0.2e1;
t478 = -t322 / 0.2e1;
t475 = t334 / 0.2e1;
t472 = -t354 / 0.2e1;
t205 = -qJD(2) * t382 - t268 * t371;
t423 = t374 * t331 + t420 * t465;
t196 = t384 * qJD(2) + (-t356 + (pkin(9) * t371 - t333) * t370) * qJD(3) + t423;
t223 = t370 * t331 + t333 * t416 + (-t371 * t419 - t375 * t417) * pkin(7);
t202 = -pkin(9) * t545 + t223;
t88 = -qJD(4) * t200 + t373 * t196 - t202 * t369;
t65 = pkin(4) * t420 - qJ(5) * t205 + qJD(5) * t298 + t88;
t206 = -qJD(2) * t381 + t298 * t516;
t87 = t369 * t196 + t373 * t202 + t262 * t414 - t270 * t415;
t69 = qJ(5) * t206 - qJD(5) * t297 + t87;
t25 = t366 * t65 + t367 * t69;
t144 = t218 * t372 - t219 * t368;
t192 = t255 * t368 + t256 * t372;
t81 = -qJD(6) * t192 + t197 * t372 - t198 * t368;
t440 = t144 - t81;
t145 = t218 * t368 + t219 * t372;
t188 = t255 * t372 - t256 * t368;
t80 = qJD(6) * t188 + t197 * t368 + t198 * t372;
t439 = t145 - t80;
t436 = qJD(2) * mrSges(3,2);
t199 = t373 * t262 - t369 * t270;
t167 = -pkin(4) * t375 + t298 * qJ(5) + t199;
t175 = -qJ(5) * t297 + t200;
t100 = t366 * t167 + t367 * t175;
t332 = pkin(3) * t431 + t371 * pkin(7);
t408 = Ifges(4,3) * t400 + t520;
t282 = pkin(3) * t545 + pkin(7) * t418;
t359 = -pkin(3) * t374 - pkin(2);
t45 = -t96 * mrSges(6,1) + t97 * mrSges(6,2);
t8 = -t32 * mrSges(7,1) + t31 * mrSges(7,2);
t254 = -pkin(3) * t276 + qJD(2) * t362;
t24 = -t366 * t69 + t367 * t65;
t99 = t367 * t167 - t366 * t175;
t396 = m(4) * t339 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t321 + mrSges(4,2) * t322 + mrSges(3,3) * t422;
t265 = pkin(4) * t297 + t332;
t217 = pkin(3) * t322 + t467;
t394 = mrSges(4,1) * t374 - mrSges(4,2) * t370;
t391 = Ifges(4,1) * t370 + t459;
t389 = Ifges(4,2) * t374 + t460;
t388 = Ifges(4,5) * t370 + Ifges(4,6) * t374;
t78 = -pkin(5) * t375 - t230 * pkin(10) + t99;
t79 = pkin(10) * t229 + t100;
t35 = -t368 * t79 + t372 * t78;
t36 = t368 * t78 + t372 * t79;
t387 = t203 * t374 - t204 * t370;
t289 = pkin(4) * t385 + t359;
t178 = -pkin(4) * t206 + t282;
t131 = -pkin(4) * t166 + t254;
t336 = mrSges(3,3) * t421 - t436;
t283 = -pkin(7) * t430 + t320;
t280 = mrSges(4,1) * t354 - mrSges(4,3) * t322;
t279 = -mrSges(4,2) * t354 + mrSges(4,3) * t321;
t278 = -pkin(7) * t405 + t307;
t253 = -mrSges(4,2) * t400 + mrSges(4,3) * t276;
t252 = mrSges(4,1) * t400 - mrSges(4,3) * t275;
t228 = mrSges(5,1) * t343 - mrSges(5,3) * t260;
t227 = -mrSges(5,2) * t343 + mrSges(5,3) * t398;
t224 = -qJD(3) * t284 + t423;
t216 = -mrSges(4,1) * t276 + mrSges(4,2) * t275;
t214 = -pkin(5) * t255 + t289;
t208 = t275 * Ifges(4,1) + t276 * Ifges(4,4) + Ifges(4,5) * t400;
t207 = t275 * Ifges(4,4) + t276 * Ifges(4,2) + Ifges(4,6) * t400;
t195 = -mrSges(5,1) * t398 + mrSges(5,2) * t260;
t184 = -pkin(5) * t229 + t265;
t177 = mrSges(6,1) * t343 - mrSges(6,3) * t191;
t176 = -mrSges(6,2) * t343 + mrSges(6,3) * t541;
t151 = -mrSges(5,2) * t400 + mrSges(5,3) * t166;
t150 = mrSges(5,1) * t400 - mrSges(5,3) * t165;
t135 = t467 + t578;
t132 = t217 + t578;
t130 = t205 * t367 + t206 * t366;
t129 = -t205 * t366 + t206 * t367;
t116 = -mrSges(6,1) * t541 + mrSges(6,2) * t191;
t102 = mrSges(7,1) * t334 - mrSges(7,3) * t114;
t101 = -mrSges(7,2) * t334 + mrSges(7,3) * t573;
t98 = -mrSges(5,1) * t166 + mrSges(5,2) * t165;
t92 = t165 * Ifges(5,1) + t166 * Ifges(5,4) + Ifges(5,5) * t400;
t91 = t165 * Ifges(5,4) + t166 * Ifges(5,2) + Ifges(5,6) * t400;
t90 = mrSges(6,1) * t400 - mrSges(6,3) * t97;
t89 = -mrSges(6,2) * t400 + mrSges(6,3) * t96;
t84 = -pkin(5) * t129 + t178;
t57 = -mrSges(7,1) * t573 + mrSges(7,2) * t114;
t56 = -pkin(5) * t96 + t131;
t50 = -qJD(6) * t159 + t129 * t372 - t130 * t368;
t49 = qJD(6) * t158 + t129 * t368 + t130 * t372;
t43 = t97 * Ifges(6,1) + t96 * Ifges(6,4) + Ifges(6,5) * t400;
t42 = t97 * Ifges(6,4) + t96 * Ifges(6,2) + Ifges(6,6) * t400;
t27 = -mrSges(7,2) * t400 + mrSges(7,3) * t32;
t26 = mrSges(7,1) * t400 - mrSges(7,3) * t31;
t21 = pkin(10) * t129 + t25;
t20 = pkin(5) * t420 - pkin(10) * t130 + t24;
t7 = t31 * Ifges(7,1) + t32 * Ifges(7,4) + Ifges(7,5) * t400;
t6 = t31 * Ifges(7,4) + t32 * Ifges(7,2) + Ifges(7,6) * t400;
t5 = -qJD(6) * t36 + t20 * t372 - t21 * t368;
t4 = qJD(6) * t35 + t20 * t368 + t21 * t372;
t1 = [(t396 * pkin(7) + t402 - 0.2e1 * t517 + t544) * t418 + (Ifges(7,4) * t159 + Ifges(7,2) * t158) * t512 + (Ifges(6,4) * t230 + Ifges(6,2) * t229) * t509 + (Ifges(7,1) * t159 + Ifges(7,4) * t158) * t513 + (Ifges(6,1) * t230 + Ifges(6,4) * t229) * t508 + t254 * (mrSges(5,1) * t297 - mrSges(5,2) * t298) + (-Ifges(5,4) * t298 - Ifges(5,2) * t297) * t494 + (-Ifges(5,1) * t298 - Ifges(5,4) * t297) * t495 + (-t142 * t205 + t143 * t206 - t297 * t63 + t298 * t64) * mrSges(5,3) + (t208 * t470 + pkin(7) * t216 + t207 * t471 + t392 * t483 + t390 * t482 + (-t203 * t370 - t204 * t374) * mrSges(4,3) + (t389 * t479 + t391 * t478 + t388 * t472 + t339 * t394 + t242 * t471 - t374 * t241 / 0.2e1 + (t263 * t370 - t264 * t374) * mrSges(4,3)) * qJD(3) + (-pkin(7) * t336 + t401 + (-0.2e1 * t506 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t456 / 0.2e1 - t455 / 0.2e1) * t371 + Ifges(7,5) * t496 + Ifges(7,6) * t497 + Ifges(6,5) * t488 + Ifges(6,6) * t489 + Ifges(5,5) * t480 + Ifges(5,6) * t481) * qJD(1) + t515) * qJD(2)) * t371 + (-t12 * t230 + t129 * t59 + t13 * t229 - t130 * t58) * mrSges(6,3) + m(7) * (t126 * t84 + t16 * t5 + t17 * t4 + t184 * t56 + t2 * t36 + t3 * t35) + m(6) * (t100 * t13 + t12 * t99 + t131 * t265 + t178 * t209 + t24 * t58 + t25 * t59) + m(5) * (t142 * t88 + t143 * t87 + t199 * t64 + t200 * t63 + t254 * t332 + t281 * t282) + m(4) * (t203 * t284 + t204 * t283 + t264 * t223 + t263 * t224) + (Ifges(5,5) * t205 + Ifges(6,5) * t130 + Ifges(5,6) * t206 + Ifges(6,6) * t129) * t473 + (t158 * t2 - t159 * t3 - t16 * t49 + t17 * t50) * mrSges(7,3) + t43 * t488 + t42 * t489 + (Ifges(6,1) * t130 + Ifges(6,4) * t129) * t490 + (Ifges(6,4) * t130 + Ifges(6,2) * t129) * t492 + t7 * t496 + t6 * t497 + (Ifges(7,1) * t49 + Ifges(7,4) * t50) * t498 + (Ifges(7,5) * t49 + Ifges(7,6) * t50) * t475 + t92 * t480 + t91 * t481 + (Ifges(5,1) * t205 + Ifges(5,4) * t206) * t484 + (Ifges(5,4) * t205 + Ifges(5,2) * t206) * t486 + t332 * t98 + t283 * t252 + t284 * t253 + t223 * t279 + t224 * t280 + t281 * (-mrSges(5,1) * t206 + mrSges(5,2) * t205) + t282 * t195 + t265 * t45 + t206 * t539 + t131 * (-mrSges(6,1) * t229 + mrSges(6,2) * t230) + t88 * t228 + t87 * t227 + t209 * (-mrSges(6,1) * t129 + mrSges(6,2) * t130) + t205 * t181 / 0.2e1 + t199 * t150 + t200 * t151 + (Ifges(7,4) * t49 + Ifges(7,2) * t50) * t500 - (t412 + t408 + t519 + t520) * t375 / 0.2e1 + t35 * t26 + t36 * t27 + t50 * t54 / 0.2e1 + t49 * t55 / 0.2e1 + (-Ifges(5,5) * t495 - Ifges(6,5) * t508 - Ifges(7,5) * t513 - Ifges(5,6) * t494 - Ifges(6,6) * t509 - Ifges(7,6) * t512 + (0.3e1 / 0.2e1 * Ifges(3,4) * t418 + (-0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - Ifges(7,3) / 0.2e1 + (m(4) * pkin(7) + t393) * pkin(7) - t410) * t420) * qJD(1) - t514 + t518 - t537) * t375 + t84 * t57 + t99 * t90 + t100 * t89 + t4 * t101 + t5 * t102 + t126 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t129 * t104 / 0.2e1 + t130 * t105 / 0.2e1 + t56 * (-mrSges(7,1) * t158 + mrSges(7,2) * t159) + t25 * t176 + t24 * t177 + t178 * t116 + t184 * t8; t546 * t116 + t547 * t228 + t548 * t227 + t549 * t57 + ((t517 + t402 - t360 / 0.2e1 + ((-m(4) * pkin(2) - mrSges(3,1) - t394) * qJD(2) - t396) * pkin(7) - t544) * t375 + ((t336 + t436) * pkin(7) + (t506 + t461 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t375) * qJD(1) + t401 - t515) * t371 + (Ifges(5,5) * t324 + Ifges(6,5) * t256 + Ifges(7,5) * t192 - Ifges(5,6) * t385 + Ifges(6,6) * t255 + Ifges(7,6) * t188 + t388) * t420 / 0.2e1) * qJD(1) + (t547 * t142 + t548 * t143 + t254 * t359 + t273 * t64 + t274 * t63 + t281 * t575) * m(5) + (Ifges(5,4) * t324 - Ifges(5,2) * t385) * t494 + (Ifges(5,1) * t324 - Ifges(5,4) * t385) * t495 + t254 * (mrSges(5,1) * t385 + mrSges(5,2) * t324) - t385 * t91 / 0.2e1 + (-Ifges(5,5) * t267 + Ifges(6,5) * t198 - Ifges(5,6) * t268 + Ifges(6,6) * t197) * t473 + (-Ifges(5,1) * t267 - Ifges(5,4) * t268) * t484 + (-Ifges(5,4) * t267 - Ifges(5,2) * t268) * t486 + (t285 / 0.2e1 - t268 / 0.2e1) * t180 + (t286 / 0.2e1 - t267 / 0.2e1) * t181 + (-Ifges(5,5) * t286 + Ifges(6,5) * t219 - Ifges(5,6) * t285 + Ifges(6,6) * t218) * t474 + (-Ifges(5,4) * t286 - Ifges(5,2) * t285) * t487 + (-Ifges(5,1) * t286 - Ifges(5,4) * t285) * t485 + (t195 * t468 + t377) * qJD(3) + (t80 / 0.2e1 - t145 / 0.2e1) * t55 + (m(4) * (-qJD(3) * t386 + t387) - t252 * t370 + t253 * t374 + (-t279 * t370 - t280 * t374) * qJD(3)) * pkin(8) + (t81 / 0.2e1 - t144 / 0.2e1) * t54 - m(4) * (t263 * t277 + t264 * t278) + (-t218 / 0.2e1 + t197 / 0.2e1) * t104 + (-t219 / 0.2e1 + t198 / 0.2e1) * t105 + (t142 * t594 - t143 * t581 - t324 * t64 - t385 * t63) * mrSges(5,3) + (mrSges(5,1) * t581 - mrSges(5,2) * t594) * t281 + (Ifges(6,1) * t198 + Ifges(6,4) * t197) * t490 + (Ifges(6,1) * t219 + Ifges(6,4) * t218) * t491 + (Ifges(6,4) * t198 + Ifges(6,2) * t197) * t492 + (Ifges(6,4) * t219 + Ifges(6,2) * t218) * t493 + (Ifges(7,1) * t80 + Ifges(7,4) * t81) * t498 + (Ifges(7,1) * t145 + Ifges(7,4) * t144) * t499 + (Ifges(7,5) * t80 + Ifges(7,6) * t81) * t475 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t476 + t389 * t482 + t391 * t483 + t207 * t470 + t370 * t208 / 0.2e1 + t359 * t98 + t324 * t92 / 0.2e1 - t316 * t195 + (mrSges(7,1) * t440 - mrSges(7,2) * t439) * t126 + (t16 * t439 - t17 * t440 + t188 * t2 - t192 * t3) * mrSges(7,3) + (-t12 * t256 + t13 * t255 - t426 * t58 + t427 * t59) * mrSges(6,3) + (-mrSges(6,1) * t427 + mrSges(6,2) * t426) * t209 + t289 * t45 - t278 * t279 - t277 * t280 + t273 * t150 + t274 * t151 + t131 * (-mrSges(6,1) * t255 + mrSges(6,2) * t256) + t256 * t43 / 0.2e1 + t255 * t42 / 0.2e1 + t214 * t8 - pkin(2) * t216 + t56 * (-mrSges(7,1) * t188 + mrSges(7,2) * t192) + (Ifges(7,4) * t80 + Ifges(7,2) * t81) * t500 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t501 + t552 * t176 + t553 * t177 + (t12 * t169 + t13 * t170 + t131 * t289 + t209 * t546 + t552 * t59 + t553 * t58) * m(6) + t554 * t102 + t555 * t101 + (t126 * t549 + t16 * t554 + t17 * t555 + t2 * t77 + t214 * t56 + t3 * t76) * m(7) + (Ifges(6,1) * t256 + Ifges(6,4) * t255) * t508 + (Ifges(6,4) * t256 + Ifges(6,2) * t255) * t509 + (Ifges(7,4) * t192 + Ifges(7,2) * t188) * t512 + (Ifges(7,1) * t192 + Ifges(7,4) * t188) * t513 + t387 * mrSges(4,3) + t76 * t26 + t77 * t27 + t169 * t90 + t170 * t89 + t188 * t6 / 0.2e1 + t192 * t7 / 0.2e1; t260 * t539 + (t373 * t150 + t369 * t151 - t322 * t195 + (t227 * t373 - t228 * t369) * qJD(4) + (-t142 * t415 + t143 * t414 + 0.2e1 * t281 * t478 + t369 * t63 + t373 * t64) * m(5)) * pkin(3) + t592 + (-Ifges(4,2) * t322 + t242 + t317) * t479 + t408 - m(5) * (t142 * t153 + t143 * t154) + (t263 * t321 + t264 * t322) * mrSges(4,3) - t518 + t241 * t477 + (Ifges(4,1) * t321 - t444) * t478 + (Ifges(4,5) * t321 - Ifges(4,6) * t322) * t472 - t339 * (mrSges(4,1) * t322 + mrSges(4,2) * t321) + t299 * t90 + t301 * t89 - t263 * t279 + t264 * t280 + t234 * t26 + t235 * t27 - t153 * t228 - t154 * t227 - t217 * t116 + t527 * t176 + (t12 * t299 + t13 * t301 - t209 * t217 + t527 * t59 + t528 * t58) * m(6) + t528 * t177 + t531 * t101 + t532 * t102 + (-t126 * t132 + t16 * t532 + t17 * t531 + t2 * t235 + t234 * t3) * m(7) - t132 * t57; (-t116 * t260 + t366 * t89 + t367 * t90) * pkin(4) + (-t209 * t467 - t58 * t67 - t59 * t68 + (t12 * t367 + t13 * t366) * pkin(4)) * m(6) + t180 * t484 + t300 * t26 + t302 * t27 + t143 * t228 - t142 * t227 + t529 * t102 + t530 * t101 + (-t126 * t135 + t16 * t529 + t17 * t530 + t2 * t302 + t3 * t300) * m(7) - t135 * t57 - t68 * t176 - t67 * t177 + t592; -t573 * t101 + t114 * t102 - t541 * t176 + t191 * t177 + t45 + t8 + (t114 * t16 - t17 * t573 + t56) * m(7) + (t191 * t58 - t541 * t59 + t131) * m(6); -t16 * t101 + t17 * t102 + t593;];
tauc  = t1(:);
