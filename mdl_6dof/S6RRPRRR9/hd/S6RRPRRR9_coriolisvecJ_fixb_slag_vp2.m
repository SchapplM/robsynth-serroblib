% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:22
% EndTime: 2019-03-09 14:09:18
% DurationCPUTime: 27.82s
% Computational Cost: add. (31999->846), mult. (83657->1202), div. (0->0), fcn. (69362->12), ass. (0->387)
t341 = sin(pkin(12));
t343 = cos(pkin(12));
t347 = sin(qJ(4));
t351 = cos(qJ(4));
t312 = t341 * t351 + t343 * t347;
t342 = sin(pkin(6));
t352 = cos(qJ(2));
t434 = t342 * t352;
t364 = t312 * t434;
t265 = qJD(1) * t364;
t304 = t312 * qJD(4);
t425 = -t265 + t304;
t433 = t343 * t351;
t311 = -t341 * t347 + t433;
t303 = t311 * qJD(4);
t346 = sin(qJ(5));
t350 = cos(qJ(5));
t373 = t350 * t311 - t312 * t346;
t188 = qJD(5) * t373 + t303 * t350 - t304 * t346;
t363 = t311 * t434;
t266 = qJD(1) * t363;
t208 = -t265 * t346 + t266 * t350;
t428 = t188 - t208;
t245 = t311 * t346 + t312 * t350;
t189 = qJD(5) * t245 + t303 * t346 + t350 * t304;
t207 = t350 * t265 + t266 * t346;
t427 = t189 - t207;
t348 = sin(qJ(2));
t384 = pkin(2) * t348 - qJ(3) * t352;
t423 = qJD(1) * t342;
t295 = t384 * t423;
t410 = t348 * t423;
t344 = cos(pkin(6));
t422 = qJD(1) * t344;
t412 = pkin(1) * t422;
t296 = -pkin(8) * t410 + t352 * t412;
t234 = t343 * t295 - t341 * t296;
t432 = t343 * t352;
t365 = (pkin(3) * t348 - pkin(9) * t432) * t342;
t205 = qJD(1) * t365 + t234;
t235 = t341 * t295 + t343 * t296;
t409 = t352 * t423;
t400 = t341 * t409;
t218 = -pkin(9) * t400 + t235;
t135 = t351 * t205 - t218 * t347;
t125 = pkin(4) * t410 - pkin(10) * t266 + t135;
t136 = t347 * t205 + t351 * t218;
t127 = -pkin(10) * t265 + t136;
t472 = pkin(9) + qJ(3);
t319 = t472 * t341;
t320 = t472 * t343;
t419 = qJD(4) * t351;
t231 = -t319 * t419 + qJD(3) * t433 + (-qJD(3) * t341 - qJD(4) * t320) * t347;
t206 = -pkin(10) * t304 + t231;
t264 = -t347 * t319 + t351 * t320;
t232 = -t312 * qJD(3) - qJD(4) * t264;
t358 = -pkin(10) * t303 + t232;
t263 = -t351 * t319 - t320 * t347;
t237 = -pkin(10) * t312 + t263;
t238 = pkin(10) * t311 + t264;
t375 = t350 * t237 - t238 * t346;
t543 = qJD(5) * t375 + (-t127 + t206) * t350 + (-t125 + t358) * t346;
t297 = pkin(8) * t409 + t348 * t412;
t258 = pkin(3) * t400 + t297;
t541 = t425 * pkin(4) - t258;
t547 = -pkin(11) * t410 + t543;
t546 = t427 * pkin(5) - t428 * pkin(11) + t541;
t164 = t237 * t346 + t238 * t350;
t338 = -pkin(3) * t343 - pkin(2);
t286 = -pkin(4) * t311 + t338;
t165 = -pkin(5) * t373 - pkin(11) * t245 + t286;
t345 = sin(qJ(6));
t349 = cos(qJ(6));
t115 = t164 * t349 + t165 * t345;
t545 = -qJD(6) * t115 - t345 * t547 + t546 * t349;
t114 = -t164 * t345 + t165 * t349;
t544 = qJD(6) * t114 + t546 * t345 + t349 * t547;
t325 = qJD(4) - t409;
t318 = qJD(5) + t325;
t334 = qJD(2) + t422;
t281 = t334 * t343 - t341 * t410;
t282 = t334 * t341 + t343 * t410;
t227 = t281 * t347 + t282 * t351;
t401 = t351 * t281 - t282 * t347;
t537 = t350 * t227 + t346 * t401;
t133 = t318 * t349 - t345 * t537;
t413 = qJD(1) * qJD(2);
t405 = t342 * t413;
t398 = t348 * t405;
t360 = qJD(2) * t363;
t175 = qJD(1) * t360 + qJD(4) * t401;
t361 = qJD(2) * t364;
t176 = -qJD(1) * t361 - qJD(4) * t227;
t538 = -t227 * t346 + t350 * t401;
t89 = qJD(5) * t538 + t175 * t350 + t176 * t346;
t56 = qJD(6) * t133 + t345 * t398 + t349 * t89;
t134 = t318 * t345 + t349 * t537;
t57 = -qJD(6) * t134 - t345 * t89 + t349 * t398;
t21 = -mrSges(7,1) * t57 + mrSges(7,2) * t56;
t269 = qJ(3) * t334 + t297;
t290 = (-pkin(2) * t352 - qJ(3) * t348 - pkin(1)) * t342;
t274 = qJD(1) * t290;
t209 = -t341 * t269 + t343 * t274;
t166 = -pkin(3) * t409 - t282 * pkin(9) + t209;
t210 = t343 * t269 + t341 * t274;
t177 = pkin(9) * t281 + t210;
t121 = t166 * t347 + t177 * t351;
t105 = pkin(10) * t401 + t121;
t429 = t350 * t105;
t120 = t351 * t166 - t177 * t347;
t104 = -pkin(10) * t227 + t120;
t95 = pkin(4) * t325 + t104;
t45 = t346 * t95 + t429;
t272 = (qJD(2) * t384 - qJD(3) * t348) * t342;
t254 = qJD(1) * t272;
t485 = pkin(1) * t344;
t411 = qJD(2) * t485;
t331 = t352 * t411;
t287 = -pkin(8) * t398 + qJD(1) * t331;
t255 = qJD(3) * t334 + t287;
t197 = t343 * t254 - t341 * t255;
t362 = qJD(2) * t365;
t167 = qJD(1) * t362 + t197;
t198 = t341 * t254 + t343 * t255;
t397 = t352 * t405;
t372 = t341 * t397;
t178 = -pkin(9) * t372 + t198;
t72 = -qJD(4) * t121 + t351 * t167 - t178 * t347;
t49 = pkin(4) * t398 - pkin(10) * t175 + t72;
t420 = qJD(4) * t347;
t71 = t166 * t419 + t347 * t167 - t177 * t420 + t351 * t178;
t51 = pkin(10) * t176 + t71;
t11 = -qJD(5) * t45 - t346 * t51 + t350 * t49;
t8 = -pkin(5) * t398 - t11;
t542 = -m(7) * t8 - t21;
t540 = t231 - t136;
t539 = t232 - t135;
t145 = qJD(6) - t538;
t43 = pkin(11) * t318 + t45;
t260 = -t334 * pkin(2) + qJD(3) - t296;
t223 = -t281 * pkin(3) + t260;
t153 = -pkin(4) * t401 + t223;
t77 = -pkin(5) * t538 - pkin(11) * t537 + t153;
t22 = -t345 * t43 + t349 * t77;
t306 = pkin(8) * t434 + t348 * t485;
t288 = t306 * t413;
t248 = pkin(3) * t372 + t288;
t143 = -t176 * pkin(4) + t248;
t90 = qJD(5) * t537 + t175 * t346 - t350 * t176;
t30 = t90 * pkin(5) - t89 * pkin(11) + t143;
t417 = qJD(5) * t350;
t418 = qJD(5) * t346;
t10 = -t105 * t418 + t346 * t49 + t350 * t51 + t95 * t417;
t7 = pkin(11) * t398 + t10;
t2 = qJD(6) * t22 + t30 * t345 + t349 * t7;
t536 = t2 * mrSges(7,2);
t23 = t345 * t77 + t349 * t43;
t3 = -qJD(6) * t23 + t30 * t349 - t345 * t7;
t535 = t3 * mrSges(7,1);
t534 = -Ifges(3,6) * t334 / 0.2e1;
t289 = qJ(3) * t344 + t306;
t229 = -t341 * t289 + t343 * t290;
t435 = t342 * t348;
t301 = t341 * t344 + t343 * t435;
t190 = -pkin(3) * t434 - t301 * pkin(9) + t229;
t230 = t343 * t289 + t341 * t290;
t300 = -t341 * t435 + t343 * t344;
t202 = pkin(9) * t300 + t230;
t128 = t351 * t190 - t347 * t202;
t240 = t300 * t347 + t301 * t351;
t112 = -pkin(4) * t434 - t240 * pkin(10) + t128;
t129 = t347 * t190 + t351 * t202;
t239 = t300 * t351 - t301 * t347;
t116 = pkin(10) * t239 + t129;
t532 = t346 * t112 + t350 * t116;
t144 = Ifges(6,4) * t538;
t386 = Ifges(7,5) * t349 - Ifges(7,6) * t345;
t366 = t145 * t386;
t464 = Ifges(7,4) * t345;
t391 = Ifges(7,1) * t349 - t464;
t367 = t134 * t391;
t463 = Ifges(7,4) * t349;
t388 = -Ifges(7,2) * t345 + t463;
t368 = t133 * t388;
t393 = mrSges(7,1) * t345 + mrSges(7,2) * t349;
t438 = t105 * t346;
t44 = t350 * t95 - t438;
t42 = -pkin(5) * t318 - t44;
t369 = t42 * t393;
t446 = t318 * Ifges(6,5);
t487 = -t349 / 0.2e1;
t488 = t345 / 0.2e1;
t454 = t537 * Ifges(6,1);
t100 = t144 + t446 + t454;
t526 = t153 * mrSges(6,2) + t100 / 0.2e1 - t44 * mrSges(6,3);
t465 = Ifges(7,4) * t134;
t67 = Ifges(7,2) * t133 + Ifges(7,6) * t145 + t465;
t132 = Ifges(7,4) * t133;
t68 = Ifges(7,1) * t134 + Ifges(7,5) * t145 + t132;
t531 = -t368 / 0.2e1 - t367 / 0.2e1 - t366 / 0.2e1 - t446 / 0.2e1 + t67 * t488 + t68 * t487 - t369 - t526 - t144 / 0.2e1;
t530 = -t22 * t345 + t23 * t349;
t529 = -t11 * mrSges(6,1) + t10 * mrSges(6,2) - Ifges(6,5) * t89 + Ifges(6,6) * t90;
t326 = Ifges(3,4) * t409;
t468 = Ifges(4,4) * t343;
t389 = -Ifges(4,2) * t341 + t468;
t469 = Ifges(4,4) * t341;
t392 = Ifges(4,1) * t343 - t469;
t471 = mrSges(4,2) * t343;
t490 = t343 / 0.2e1;
t491 = -t341 / 0.2e1;
t528 = -(t209 * t343 + t210 * t341) * mrSges(4,3) + t260 * (mrSges(4,1) * t341 + t471) + Ifges(3,1) * t410 / 0.2e1 + t326 / 0.2e1 + t334 * Ifges(3,5) + t281 * t389 / 0.2e1 + t282 * t392 / 0.2e1 - t296 * mrSges(3,3) + (Ifges(4,4) * t282 + Ifges(4,2) * t281 - Ifges(4,6) * t409) * t491 + (Ifges(4,1) * t282 + Ifges(4,4) * t281 - Ifges(4,5) * t409) * t490;
t193 = qJD(4) * t239 + t360;
t421 = qJD(2) * t342;
t407 = t348 * t421;
t298 = -pkin(8) * t407 + t331;
t280 = qJD(3) * t344 + t298;
t216 = t343 * t272 - t341 * t280;
t183 = t216 + t362;
t217 = t341 * t272 + t343 * t280;
t408 = t352 * t421;
t399 = t341 * t408;
t196 = -pkin(9) * t399 + t217;
t81 = -qJD(4) * t129 + t351 * t183 - t196 * t347;
t65 = pkin(4) * t407 - pkin(10) * t193 + t81;
t194 = -qJD(4) * t240 - t361;
t80 = t347 * t183 + t190 * t419 + t351 * t196 - t202 * t420;
t70 = pkin(10) * t194 + t80;
t20 = -qJD(5) * t532 - t346 * t70 + t350 * t65;
t527 = -t72 * mrSges(5,1) + t71 * mrSges(5,2) - Ifges(5,5) * t175 - Ifges(5,6) * t176;
t109 = pkin(5) * t537 - pkin(11) * t538;
t457 = t145 * Ifges(7,3);
t458 = t134 * Ifges(7,5);
t459 = t133 * Ifges(7,6);
t66 = t457 + t458 + t459;
t445 = t318 * Ifges(6,6);
t456 = t538 * Ifges(6,2);
t466 = Ifges(6,4) * t537;
t99 = t445 + t456 + t466;
t525 = t23 * mrSges(7,2) + t45 * mrSges(6,3) - t66 / 0.2e1 + t99 / 0.2e1 - t153 * mrSges(6,1) - t22 * mrSges(7,1);
t523 = Ifges(6,2) / 0.2e1;
t54 = Ifges(7,6) * t57;
t55 = Ifges(7,5) * t56;
t14 = Ifges(7,3) * t90 + t54 + t55;
t522 = t14 / 0.2e1;
t521 = t56 / 0.2e1;
t520 = t57 / 0.2e1;
t519 = t68 / 0.2e1;
t518 = t90 / 0.2e1;
t517 = pkin(1) * mrSges(3,1);
t516 = pkin(1) * mrSges(3,2);
t515 = -t133 / 0.2e1;
t514 = t133 / 0.2e1;
t513 = -t134 / 0.2e1;
t512 = t134 / 0.2e1;
t510 = -t145 / 0.2e1;
t509 = t145 / 0.2e1;
t508 = t538 / 0.2e1;
t507 = t537 / 0.2e1;
t374 = t350 * t239 - t240 * t346;
t506 = t374 / 0.2e1;
t169 = t239 * t346 + t240 * t350;
t505 = t169 / 0.2e1;
t504 = t175 / 0.2e1;
t503 = t176 / 0.2e1;
t502 = -t401 / 0.2e1;
t501 = t401 / 0.2e1;
t500 = -t227 / 0.2e1;
t499 = t227 / 0.2e1;
t498 = t239 / 0.2e1;
t497 = t240 / 0.2e1;
t496 = t300 / 0.2e1;
t495 = t301 / 0.2e1;
t494 = t318 / 0.2e1;
t493 = -t325 / 0.2e1;
t492 = t325 / 0.2e1;
t489 = -t345 / 0.2e1;
t486 = t349 / 0.2e1;
t484 = pkin(1) * t352;
t481 = t2 * t349;
t480 = t3 * t345;
t477 = t89 * Ifges(6,1);
t476 = t89 * Ifges(6,4);
t475 = t90 * Ifges(6,4);
t470 = Ifges(3,4) * t348;
t467 = Ifges(5,4) * t227;
t451 = t22 * t349;
t447 = t287 * mrSges(3,2);
t441 = t341 * Ifges(4,6);
t440 = t343 * Ifges(4,5);
t139 = mrSges(6,1) * t318 - mrSges(6,3) * t537;
t92 = -mrSges(7,1) * t133 + mrSges(7,2) * t134;
t439 = t139 - t92;
t431 = t345 * t188;
t430 = t349 * t188;
t426 = -mrSges(3,1) * t334 - mrSges(4,1) * t281 + mrSges(4,2) * t282 + mrSges(3,3) * t410;
t424 = t266 - t303;
t270 = mrSges(4,1) * t372 + t397 * t471;
t299 = pkin(8) * t408 + t348 * t411;
t416 = qJD(6) * t345;
t415 = qJD(6) * t349;
t259 = pkin(3) * t399 + t299;
t37 = t90 * mrSges(6,1) + t89 * mrSges(6,2);
t123 = -t176 * mrSges(5,1) + t175 * mrSges(5,2);
t186 = -t208 * t345 + t349 * t410;
t404 = t186 + t431;
t187 = t208 * t349 + t345 * t410;
t403 = -t187 + t430;
t396 = t535 - t536;
t395 = -t2 * t345 - t3 * t349;
t394 = mrSges(7,1) * t349 - mrSges(7,2) * t345;
t390 = Ifges(7,1) * t345 + t463;
t387 = Ifges(7,2) * t349 + t464;
t385 = Ifges(7,5) * t345 + Ifges(7,6) * t349;
t383 = t23 * t345 + t451;
t59 = -pkin(11) * t434 + t532;
t335 = pkin(8) * t435;
t291 = t335 + (-pkin(2) - t484) * t344;
t241 = -t300 * pkin(3) + t291;
t181 = -t239 * pkin(4) + t241;
t93 = -pkin(5) * t374 - t169 * pkin(11) + t181;
t32 = t345 * t93 + t349 * t59;
t31 = -t345 * t59 + t349 * t93;
t96 = -mrSges(7,2) * t145 + mrSges(7,3) * t133;
t97 = mrSges(7,1) * t145 - mrSges(7,3) * t134;
t381 = -t345 * t97 + t349 * t96;
t60 = t350 * t112 - t346 * t116;
t78 = t125 * t350 - t127 * t346;
t156 = -pkin(4) * t194 + t259;
t138 = -mrSges(6,2) * t318 + mrSges(6,3) * t538;
t371 = -t138 - t381;
t154 = -t345 * t169 - t349 * t434;
t370 = -t349 * t169 + t345 * t434;
t19 = t112 * t417 - t116 * t418 + t346 * t65 + t350 * t70;
t15 = t56 * Ifges(7,4) + t57 * Ifges(7,2) + t90 * Ifges(7,6);
t16 = t56 * Ifges(7,1) + t57 * Ifges(7,4) + t90 * Ifges(7,5);
t321 = Ifges(6,3) * t398;
t359 = mrSges(7,3) * t481 + qJD(6) * t369 + t15 * t486 + t16 * t488 + t387 * t520 + t390 * t521 - t8 * t394 + t321 - t67 * t416 / 0.2e1 + t415 * t519 + t385 * t518 - t529 + (t368 + t367 + t366) * qJD(6) / 0.2e1;
t26 = mrSges(7,1) * t90 - mrSges(7,3) * t56;
t27 = -mrSges(7,2) * t90 + mrSges(7,3) * t57;
t357 = m(7) * (-t22 * t415 - t23 * t416 - t480 + t481) + t349 * t27 - t345 * t26 - t97 * t415 - t96 * t416;
t356 = -t458 / 0.2e1 - t457 / 0.2e1 - t459 / 0.2e1 + t445 / 0.2e1 + t466 / 0.2e1 + t525;
t353 = t120 * mrSges(5,1) + t209 * mrSges(4,1) + t44 * mrSges(6,1) + t325 * Ifges(5,3) + t227 * Ifges(5,5) + t401 * Ifges(5,6) - Ifges(4,3) * t409 / 0.2e1 + Ifges(4,6) * t281 + Ifges(4,5) * t282 + t534 - (t352 * Ifges(3,2) + t470) * t423 / 0.2e1 + t318 * Ifges(6,3) + t537 * Ifges(6,5) + t538 * Ifges(6,6) - t121 * mrSges(5,2) - t210 * mrSges(4,2) - t297 * mrSges(3,3) - t45 * mrSges(6,2);
t323 = Ifges(3,5) * t397;
t322 = Ifges(5,3) * t398;
t305 = t344 * t484 - t335;
t294 = -t334 * mrSges(3,2) + mrSges(3,3) * t409;
t276 = (mrSges(4,1) * t348 - mrSges(4,3) * t432) * t405;
t275 = (-mrSges(4,3) * t341 * t352 - mrSges(4,2) * t348) * t405;
t251 = -mrSges(4,1) * t409 - t282 * mrSges(4,3);
t250 = mrSges(4,2) * t409 + t281 * mrSges(4,3);
t243 = (Ifges(4,5) * t348 + t352 * t392) * t405;
t242 = (Ifges(4,6) * t348 + t352 * t389) * t405;
t222 = Ifges(5,4) * t401;
t204 = mrSges(5,1) * t325 - mrSges(5,3) * t227;
t203 = -mrSges(5,2) * t325 + mrSges(5,3) * t401;
t158 = -mrSges(5,2) * t398 + mrSges(5,3) * t176;
t157 = mrSges(5,1) * t398 - mrSges(5,3) * t175;
t152 = -mrSges(5,1) * t401 + mrSges(5,2) * t227;
t142 = t227 * Ifges(5,1) + t325 * Ifges(5,5) + t222;
t141 = Ifges(5,2) * t401 + t325 * Ifges(5,6) + t467;
t119 = t175 * Ifges(5,1) + t176 * Ifges(5,4) + Ifges(5,5) * t398;
t118 = t175 * Ifges(5,4) + t176 * Ifges(5,2) + Ifges(5,6) * t398;
t108 = -mrSges(6,1) * t538 + mrSges(6,2) * t537;
t107 = qJD(5) * t164 + t206 * t346 - t350 * t358;
t103 = qJD(5) * t169 + t193 * t346 - t350 * t194;
t102 = qJD(5) * t374 + t193 * t350 + t194 * t346;
t91 = pkin(4) * t227 + t109;
t83 = -mrSges(6,2) * t398 - mrSges(6,3) * t90;
t82 = mrSges(6,1) * t398 - mrSges(6,3) * t89;
t75 = -pkin(5) * t410 - t78;
t74 = qJD(6) * t370 - t345 * t102 + t349 * t407;
t73 = qJD(6) * t154 + t349 * t102 + t345 * t407;
t58 = pkin(5) * t434 - t60;
t47 = t104 * t350 - t438;
t46 = t104 * t346 + t429;
t38 = pkin(5) * t103 - pkin(11) * t102 + t156;
t36 = Ifges(6,5) * t398 - t475 + t477;
t35 = -t90 * Ifges(6,2) + Ifges(6,6) * t398 + t476;
t29 = t109 * t345 + t349 * t44;
t28 = t109 * t349 - t345 * t44;
t25 = t345 * t91 + t349 * t47;
t24 = -t345 * t47 + t349 * t91;
t18 = -pkin(5) * t407 - t20;
t17 = pkin(11) * t407 + t19;
t5 = -qJD(6) * t32 - t17 * t345 + t349 * t38;
t4 = qJD(6) * t31 + t17 * t349 + t345 * t38;
t1 = [(t10 * t374 - t11 * t169) * mrSges(6,3) + t143 * (-mrSges(6,1) * t374 + mrSges(6,2) * t169) + t89 * (Ifges(6,1) * t169 + Ifges(6,4) * t374) / 0.2e1 + t374 * t536 - t374 * t522 + m(6) * (t10 * t532 + t11 * t60 + t143 * t181 + t153 * t156 + t19 * t45 + t20 * t44) + t532 * t83 + (t154 * t2 - t22 * t73 + t23 * t74 + t3 * t370) * mrSges(7,3) + t8 * (-mrSges(7,1) * t154 - mrSges(7,2) * t370) - t370 * t16 / 0.2e1 - t374 * t535 + (-Ifges(7,5) * t370 + Ifges(7,6) * t154 - Ifges(7,3) * t374) * t518 + (-Ifges(7,4) * t370 + Ifges(7,2) * t154 - Ifges(7,6) * t374) * t520 + (-Ifges(7,1) * t370 + Ifges(7,4) * t154 - Ifges(7,5) * t374) * t521 - t90 * (Ifges(6,4) * t169 + Ifges(6,2) * t374) / 0.2e1 + (-t321 / 0.2e1 - t322 / 0.2e1 + t287 * mrSges(3,3) - t197 * mrSges(4,1) + t198 * mrSges(4,2) + t527 + t529) * t434 + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t509 + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t514 + (Ifges(6,1) * t507 + Ifges(6,4) * t508 + Ifges(6,5) * t494 + t526) * t102 + m(4) * (t197 * t229 + t198 * t230 + t209 * t216 + t210 * t217 + t260 * t299 + t288 * t291) + m(3) * (t287 * t306 - t288 * t305 - t296 * t299 + t297 * t298) + m(5) * (t120 * t81 + t121 * t80 + t128 * t72 + t129 * t71 + t223 * t259 + t241 * t248) + m(7) * (t18 * t42 + t2 * t32 + t22 * t5 + t23 * t4 + t3 * t31 + t58 * t8) + (-Ifges(6,4) * t507 + Ifges(7,5) * t512 - Ifges(6,2) * t508 - Ifges(6,6) * t494 + Ifges(7,6) * t514 + Ifges(7,3) * t509 - t525) * t103 + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t512 + (t528 * t352 + (t534 + t353) * t348) * t421 + ((-Ifges(3,6) * t344 + Ifges(5,5) * t497 + Ifges(5,6) * t498 + Ifges(6,5) * t505 + Ifges(6,6) * t506 + Ifges(4,5) * t495 + Ifges(4,6) * t496 - t306 * mrSges(3,3) + (-0.2e1 * t517 - 0.3e1 / 0.2e1 * t470) * t342) * t348 + (Ifges(3,5) * t344 / 0.2e1 - t305 * mrSges(3,3) + (Ifges(4,4) * t301 + Ifges(4,2) * t300) * t491 + (Ifges(4,1) * t301 + Ifges(4,4) * t300) * t490 + (-0.2e1 * t516 + (-0.3e1 / 0.2e1 * t440 + 0.3e1 / 0.2e1 * t441 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t352) * t342 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1) * t435) * t352) * t405 + (t323 / 0.2e1 - t447) * t344 + (-mrSges(3,1) * t344 - mrSges(4,1) * t300 + mrSges(4,2) * t301 + mrSges(3,3) * t435) * t288 + t426 * t299 + t298 * t294 + t291 * t270 + t230 * t275 + t229 * t276 + t248 * (-mrSges(5,1) * t239 + mrSges(5,2) * t240) + t217 * t250 + t216 * t251 + t259 * t152 + t241 * t123 + t223 * (-mrSges(5,1) * t194 + mrSges(5,2) * t193) + t80 * t203 + t81 * t204 + t194 * t141 / 0.2e1 + t193 * t142 / 0.2e1 + t181 * t37 + t156 * t108 + t128 * t157 + t129 * t158 + t154 * t15 / 0.2e1 + t19 * t138 + t20 * t139 + t4 * t96 + t5 * t97 + t18 * t92 + t60 * t82 + t42 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t74 * t67 / 0.2e1 + t58 * t21 + t31 * t26 + t32 * t27 + (-t120 * t193 + t121 * t194 + t239 * t71 - t240 * t72) * mrSges(5,3) + (Ifges(5,5) * t193 + Ifges(5,6) * t194) * t492 + t243 * t495 + t242 * t496 + t119 * t497 + t118 * t498 + (Ifges(5,1) * t193 + Ifges(5,4) * t194) * t499 + (Ifges(5,4) * t193 + Ifges(5,2) * t194) * t501 + (Ifges(5,4) * t240 + Ifges(5,2) * t239) * t503 + (Ifges(5,1) * t240 + Ifges(5,4) * t239) * t504 + t36 * t505 + t35 * t506 + t73 * t519 + (-t197 * t301 + t198 * t300) * mrSges(4,3); (t99 - t66) * (t207 / 0.2e1 - t189 / 0.2e1) - t537 * (Ifges(6,1) * t208 - Ifges(6,4) * t207) / 0.2e1 + (t16 * t486 + t15 * t489 + t386 * t518 + t388 * t520 + t391 * t521 + t8 * t393 + t143 * mrSges(6,2) - t475 / 0.2e1 + t477 / 0.2e1 + t36 / 0.2e1 - t11 * mrSges(6,3) + t395 * mrSges(7,3) + (-mrSges(7,3) * t530 + t385 * t510 + t387 * t515 + t390 * t513 + t394 * t42 + t487 * t67 + t489 * t68) * qJD(6)) * t245 - t538 * (Ifges(6,4) * t208 - Ifges(6,2) * t207) / 0.2e1 + ((-t326 / 0.2e1 + ((Ifges(4,2) * t343 + t469) * t491 + (Ifges(4,1) * t341 + t468) * t490) * qJD(2) + (t516 + (t440 / 0.2e1 - t441 / 0.2e1) * t352) * t423 - t528) * t352 + ((-qJD(2) + t334 / 0.2e1) * Ifges(3,6) - t353 + (t517 + t470 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t352) * t423 + (Ifges(4,5) * t341 + Ifges(5,5) * t312 + Ifges(6,5) * t245 + Ifges(4,6) * t343 + Ifges(5,6) * t311 + Ifges(6,6) * t373) * qJD(2) / 0.2e1) * t348) * t423 - (t54 / 0.2e1 + t55 / 0.2e1 + t522 - t35 / 0.2e1 + t143 * mrSges(6,1) - t476 / 0.2e1 - t10 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + t523) * t90 + t396) * t373 + (-t304 / 0.2e1 + t265 / 0.2e1) * t141 + (Ifges(5,5) * t266 - Ifges(5,6) * t265) * t493 + (Ifges(5,1) * t266 - Ifges(5,4) * t265) * t500 + (Ifges(5,4) * t266 - Ifges(5,2) * t265) * t502 - (t21 - t82) * t375 + (Ifges(5,5) * t303 - Ifges(5,6) * t304) * t492 + (Ifges(5,1) * t303 - Ifges(5,4) * t304) * t499 + (Ifges(5,4) * t303 - Ifges(5,2) * t304) * t501 + t544 * t96 + (t114 * t3 + t115 * t2 - t375 * t8 + (t107 - t75) * t42 + t544 * t23 + t545 * t22) * m(7) + t545 * t97 - t447 + t543 * t138 + (t10 * t164 + t11 * t375 + t143 * t286 + t543 * t45 + (-t107 - t78) * t44 + t541 * t153) * m(6) + (t288 * mrSges(4,2) + t243 / 0.2e1 - qJ(3) * t276 - qJD(3) * t251 - t197 * mrSges(4,3)) * t341 + (-t208 / 0.2e1 + t188 / 0.2e1) * t100 + (qJ(3) * t275 + qJD(3) * t250 + t198 * mrSges(4,3) - t288 * mrSges(4,1) + t242 / 0.2e1) * t343 - t318 * (Ifges(6,5) * t208 - Ifges(6,6) * t207) / 0.2e1 + (t303 / 0.2e1 - t266 / 0.2e1) * t142 - t439 * t107 + (-t431 / 0.2e1 - t186 / 0.2e1) * t67 + (t430 / 0.2e1 - t187 / 0.2e1) * t68 + (-t427 * t45 - t428 * t44) * mrSges(6,3) + t338 * t123 + (mrSges(5,1) * t425 - mrSges(5,2) * t424) * t223 + (t120 * t424 - t121 * t425 + t311 * t71 - t312 * t72) * mrSges(5,3) - t426 * t297 + (mrSges(7,1) * t427 - mrSges(7,3) * t403) * t22 + (-mrSges(7,2) * t427 - mrSges(7,3) * t404) * t23 + (mrSges(6,1) * t427 + mrSges(6,2) * t428) * t153 + t311 * t118 / 0.2e1 + t248 * (-mrSges(5,1) * t311 + mrSges(5,2) * t312) + t312 * t119 / 0.2e1 - t296 * t294 - t288 * mrSges(3,1) + t286 * t37 - pkin(2) * t270 + t263 * t157 + t264 * t158 - t235 * t250 - t234 * t251 - t258 * t152 + (mrSges(7,1) * t404 + mrSges(7,2) * t403) * t42 + t164 * t83 - t78 * t139 + t114 * t26 + t115 * t27 - t75 * t92 + t323 + (-t209 * t234 - t210 * t235 - t260 * t297 - pkin(2) * t288 + (-t209 * t341 + t210 * t343) * qJD(3) + (-t197 * t341 + t198 * t343) * qJ(3)) * m(4) + (Ifges(6,5) * t188 - Ifges(6,6) * t189) * t494 + (Ifges(5,4) * t312 + Ifges(5,2) * t311) * t503 + (Ifges(5,1) * t312 + Ifges(5,4) * t311) * t504 + (Ifges(6,1) * t188 - Ifges(6,4) * t189) * t507 + (Ifges(6,4) * t188 - Ifges(6,2) * t189) * t508 + (Ifges(7,5) * t430 - Ifges(7,6) * t431 + Ifges(7,3) * t189) * t509 + (Ifges(7,5) * t187 + Ifges(7,6) * t186 + Ifges(7,3) * t207) * t510 + (Ifges(7,1) * t430 - Ifges(7,4) * t431 + Ifges(7,5) * t189) * t512 + (Ifges(7,1) * t187 + Ifges(7,4) * t186 + Ifges(7,5) * t207) * t513 + (Ifges(7,4) * t430 - Ifges(7,2) * t431 + Ifges(7,6) * t189) * t514 + (Ifges(7,4) * t187 + Ifges(7,2) * t186 + Ifges(7,6) * t207) * t515 + t539 * t204 + t540 * t203 + (t120 * t539 + t121 * t540 - t223 * t258 + t248 * t338 + t263 * t72 + t264 * t71) * m(5) + t541 * t108; t381 * qJD(6) + t439 * t537 + t371 * t538 - t401 * t203 + t227 * t204 - t281 * t250 + t282 * t251 + t349 * t26 + t345 * t27 + t123 + t270 + t37 + (t145 * t530 - t537 * t42 - t395) * m(7) + (t44 * t537 - t45 * t538 + t143) * m(6) + (t120 * t227 - t121 * t401 + t248) * m(5) + (t209 * t282 - t210 * t281 + t288) * m(4); (-t454 / 0.2e1 + t531) * t538 + (-t227 * t108 + t346 * t83 + t350 * t82 + ((m(7) * t42 - t439) * t346 + (m(7) * t530 - t371) * t350) * qJD(5) + (t10 * t346 + t11 * t350 + 0.2e1 * t153 * t500 + (-t346 * t44 + t350 * t45) * qJD(5)) * m(6)) * pkin(4) + t359 + (-t145 * t451 + (-t145 * t23 - t3) * t345) * mrSges(7,3) + t357 * (pkin(4) * t346 + pkin(11)) - m(7) * (t22 * t24 + t23 * t25 + t42 * t46) - t527 + (t456 / 0.2e1 + t356) * t537 - t223 * (mrSges(5,1) * t227 + mrSges(5,2) * t401) + (t120 * t401 + t121 * t227) * mrSges(5,3) + (Ifges(5,5) * t401 - Ifges(5,6) * t227) * t493 + (Ifges(5,1) * t401 - t467) * t500 - m(6) * (-t44 * t46 + t45 * t47) + t439 * t46 - t120 * t203 + t121 * t204 - t47 * t138 - t25 * t96 - t24 * t97 + t322 + t141 * t499 + (-Ifges(5,2) * t227 + t142 + t222) * t502 - t542 * (-pkin(4) * t350 - pkin(5)); t359 + t356 * t537 + t357 * pkin(11) + ((t523 - Ifges(6,1) / 0.2e1) * t537 + t383 * mrSges(7,3) + t531) * t538 + t439 * t45 + (-qJD(6) * t383 - t480) * mrSges(7,3) - t44 * t138 - m(7) * (t22 * t28 + t23 * t29 + t42 * t45) - t29 * t96 - t28 * t97 + t542 * pkin(5); -t42 * (mrSges(7,1) * t134 + mrSges(7,2) * t133) + (Ifges(7,1) * t133 - t465) * t513 + t67 * t512 + (Ifges(7,5) * t133 - Ifges(7,6) * t134) * t510 - t22 * t96 + t23 * t97 + (t133 * t22 + t134 * t23) * mrSges(7,3) + t396 + t14 + (-Ifges(7,2) * t134 + t132 + t68) * t515;];
tauc  = t1(:);
