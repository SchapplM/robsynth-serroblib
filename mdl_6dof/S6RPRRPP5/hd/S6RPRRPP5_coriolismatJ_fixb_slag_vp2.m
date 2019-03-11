% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:52
% EndTime: 2019-03-09 04:43:02
% DurationCPUTime: 6.69s
% Computational Cost: add. (11616->570), mult. (23601->745), div. (0->0), fcn. (24487->6), ass. (0->279)
t478 = -m(7) / 0.2e1;
t497 = Ifges(6,4) + Ifges(5,5);
t318 = sin(qJ(4));
t320 = cos(qJ(4));
t316 = sin(pkin(9));
t317 = cos(pkin(9));
t319 = sin(qJ(3));
t453 = cos(qJ(3));
t275 = t316 * t453 + t319 * t317;
t411 = t275 * t318;
t251 = mrSges(7,3) * t411;
t273 = t316 * t319 - t317 * t453;
t174 = t273 * mrSges(7,2) + t251;
t410 = t275 * t320;
t179 = -t273 * mrSges(7,1) - mrSges(7,3) * t410;
t455 = t318 / 0.2e1;
t507 = t320 / 0.2e1;
t402 = t174 * t507 + t179 * t455;
t382 = -pkin(2) * t317 - pkin(1);
t169 = pkin(3) * t273 - pkin(8) * t275 + t382;
t443 = pkin(7) + qJ(2);
t279 = t443 * t317;
t370 = t443 * t316;
t193 = t279 * t453 - t319 * t370;
t397 = -t320 * t169 + t318 * t193;
t56 = -pkin(4) * t273 + t397;
t439 = t56 - t397;
t266 = t273 * qJ(5);
t74 = t169 * t318 + t193 * t320;
t55 = t266 + t74;
t440 = t55 - t74;
t409 = t318 * qJ(6);
t249 = t275 * t409;
t41 = t249 + t55;
t53 = t249 + t74;
t441 = t41 - t53;
t468 = pkin(4) + pkin(5);
t52 = qJ(6) * t410 - t397;
t34 = -t273 * t468 - t52;
t442 = t34 + t52;
t480 = -m(6) / 0.2e1;
t511 = (t318 * t439 + t320 * t440) * t480 + (t318 * t442 + t320 * t441) * t478 - t402;
t406 = t320 * qJ(5);
t407 = t318 * t468;
t277 = t406 - t407;
t280 = pkin(8) * t318 - t409;
t284 = (pkin(8) - qJ(6)) * t320;
t354 = -t318 * t34 - t320 * t41;
t510 = t402 + ((-t280 * t320 + t284 * t318) * t275 + t354) * t478;
t314 = t318 ^ 2;
t315 = t320 ^ 2;
t509 = (mrSges(6,2) + mrSges(5,3)) * (-t314 / 0.2e1 - t315 / 0.2e1);
t421 = t320 * mrSges(6,3);
t285 = t318 * mrSges(6,1) - t421;
t160 = t285 * t275;
t508 = t160 / 0.2e1;
t496 = Ifges(5,6) + Ifges(7,6);
t433 = Ifges(6,5) * t320;
t289 = Ifges(6,3) * t318 + t433;
t434 = Ifges(7,4) * t320;
t291 = Ifges(7,2) * t318 + t434;
t505 = -t291 / 0.2e1 - t289 / 0.2e1;
t418 = qJ(5) * t318;
t504 = t320 * t468 + t418;
t469 = m(6) + m(7);
t465 = t273 / 0.2e1;
t502 = -t275 / 0.2e1;
t501 = t275 / 0.2e1;
t498 = -mrSges(5,1) - mrSges(6,1);
t394 = t314 + t315;
t495 = mrSges(7,3) * t394;
t175 = -t273 * mrSges(5,2) - mrSges(5,3) * t411;
t388 = mrSges(6,2) * t411;
t182 = t273 * mrSges(6,3) - t388;
t399 = t175 + t182;
t304 = t320 * mrSges(7,2);
t425 = t318 * mrSges(7,1);
t286 = t304 - t425;
t306 = Ifges(6,5) * t318;
t494 = Ifges(6,1) * t320 + t306;
t493 = Ifges(6,6) * t318 + t320 * t497;
t308 = Ifges(7,4) * t318;
t492 = Ifges(7,1) * t320 + t308;
t310 = Ifges(5,4) * t320;
t491 = -Ifges(5,2) * t318 + t310;
t298 = Ifges(5,1) * t318 + t310;
t489 = t318 * pkin(4) - t406;
t488 = Ifges(6,3) * t320 - t306;
t487 = Ifges(7,2) * t320 - t308;
t282 = mrSges(7,1) * t320 + mrSges(7,2) * t318;
t360 = t320 * mrSges(6,1) + t318 * mrSges(6,3);
t462 = -t360 / 0.2e1;
t372 = t462 - t282 / 0.2e1;
t486 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t384 = -Ifges(5,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t366 = -Ifges(6,6) / 0.2e1 - t384;
t485 = Ifges(6,6) * t501 + t273 * t505 - t366 * t275 + t491 * t465 + t496 * t502;
t250 = t273 * t406;
t387 = -mrSges(7,1) / 0.2e1 - mrSges(6,1) / 0.2e1;
t414 = t273 * t318;
t390 = pkin(4) * t414;
t474 = -mrSges(6,3) / 0.2e1;
t475 = -mrSges(7,2) / 0.2e1;
t476 = mrSges(5,2) / 0.2e1;
t477 = m(7) / 0.2e1;
t479 = m(6) / 0.2e1;
t483 = (-t250 + t390) * t479 + (t273 * t407 - t250) * t477 + ((t474 + t475 + t476) * t320 + (mrSges(5,1) / 0.2e1 - t387) * t318) * t273;
t482 = 2 * qJD(4);
t481 = m(5) / 0.2e1;
t473 = mrSges(7,3) / 0.2e1;
t192 = t279 * t319 + t453 * t370;
t185 = t318 * t192;
t446 = pkin(8) * t273;
t447 = pkin(3) * t275;
t188 = t446 + t447;
t38 = -t185 + (qJ(6) * t273 - t188) * t320 - t468 * t275;
t472 = t38 / 0.2e1;
t471 = -t53 / 0.2e1;
t77 = t188 * t320 + t185;
t64 = -pkin(4) * t275 - t77;
t470 = t64 / 0.2e1;
t413 = t273 * t320;
t176 = -t275 * mrSges(7,1) + mrSges(7,3) * t413;
t466 = -t176 / 0.2e1;
t460 = t282 / 0.2e1;
t459 = t284 / 0.2e1;
t456 = -t318 / 0.2e1;
t452 = m(6) * t489;
t271 = -pkin(5) * t318 - t489;
t451 = m(7) * t271;
t450 = m(7) * t275;
t449 = m(7) * t277;
t448 = m(7) * t318;
t445 = -mrSges(6,2) + mrSges(7,3);
t444 = mrSges(7,2) + mrSges(6,3);
t435 = Ifges(5,4) * t318;
t422 = t320 * mrSges(5,2);
t287 = t318 * mrSges(5,1) + t422;
t162 = t287 * t275;
t180 = t273 * mrSges(5,1) - mrSges(5,3) * t410;
t181 = -t273 * mrSges(6,1) + mrSges(6,2) * t410;
t398 = -t180 + t181;
t161 = t286 * t275;
t400 = -t160 + t161;
t415 = t192 * t275;
t62 = t275 * t277 - t192;
t81 = t275 * t489 + t192;
t5 = (mrSges(4,3) * t275 + t162 - t400) * t275 + (mrSges(4,3) * t273 + (-t174 - t399) * t320 + (-t179 - t398) * t318) * t273 + m(7) * (t273 * t354 - t275 * t62) + m(6) * (t275 * t81 + (-t318 * t56 - t320 * t55) * t273) + m(5) * (t415 + (-t318 * t397 - t320 * t74) * t273) + m(4) * (-t193 * t273 + t415) + (m(3) * qJ(2) + mrSges(3,3)) * (t316 ^ 2 + t317 ^ 2);
t429 = qJD(1) * t5;
t268 = t275 * mrSges(4,1);
t173 = -t275 * mrSges(5,2) + mrSges(5,3) * t414;
t177 = t275 * mrSges(5,1) + mrSges(5,3) * t413;
t178 = -t275 * mrSges(6,1) - mrSges(6,2) * t413;
t172 = mrSges(7,2) * t275 - mrSges(7,3) * t414;
t183 = mrSges(6,2) * t414 + mrSges(6,3) * t275;
t375 = -t172 / 0.2e1 - t183 / 0.2e1;
t78 = t318 * t188 - t320 * t192;
t61 = t275 * qJ(5) + t78;
t44 = -t273 * t409 + t61;
t324 = (t173 / 0.2e1 - t375) * t318 + (t466 - t178 / 0.2e1 + t177 / 0.2e1) * t320 + (t318 * t78 + t320 * t77) * t481 + (t318 * t61 - t320 * t64) * t479 + (t318 * t44 - t320 * t38) * t477;
t270 = pkin(3) + t504;
t351 = -t280 * t318 - t284 * t320;
t369 = t394 * t446;
t355 = pkin(4) * t320 + t418;
t278 = -pkin(3) - t355;
t412 = t275 * t278;
t325 = -m(5) * (-t369 - t447) / 0.2e1 + (-t369 + t412) * t480 + (-t270 * t275 + t273 * t351) * t478;
t386 = t473 - mrSges(6,2) / 0.2e1;
t423 = t320 * mrSges(5,1);
t424 = t318 * mrSges(5,2);
t8 = t268 + (t423 / 0.2e1 - t424 / 0.2e1 - t372) * t275 + (-mrSges(4,2) + t394 * (mrSges(5,3) / 0.2e1 - t386)) * t273 + t324 + t325;
t428 = qJD(1) * t8;
t427 = t273 * Ifges(7,5);
t157 = t285 * t273;
t158 = t286 * t273;
t159 = t287 * t273;
t104 = t275 * t492 - t427;
t106 = t273 * Ifges(6,4) + t275 * t494;
t299 = Ifges(5,1) * t320 - t435;
t108 = t273 * Ifges(5,5) + t275 * t299;
t385 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t367 = -Ifges(7,5) / 0.2e1 + t385;
t336 = -t104 / 0.2e1 - t106 / 0.2e1 - t108 / 0.2e1 - t367 * t273;
t337 = Ifges(7,5) * t502 + t367 * t275 + t497 * t501 - (t299 + t492 + t494) * t273 / 0.2e1;
t256 = Ifges(7,4) * t410;
t100 = Ifges(7,2) * t411 - t273 * Ifges(7,6) + t256;
t102 = t273 * Ifges(5,6) + t275 * t491;
t255 = Ifges(6,5) * t410;
t98 = Ifges(6,6) * t273 + Ifges(6,3) * t411 + t255;
t364 = -t100 / 0.2e1 + t102 / 0.2e1 - t98 / 0.2e1;
t339 = -t250 - t193;
t63 = t414 * t468 + t339;
t82 = -t339 - t390;
t3 = -t81 * t157 - t62 * t158 + t82 * t160 + t63 * t161 + t41 * t172 + t74 * t173 + t44 * t174 + t78 * t175 + t34 * t176 - t397 * t177 + t56 * t178 + t38 * t179 + t77 * t180 + t64 * t181 + t61 * t182 + t55 * t183 - t192 * t159 + t193 * t162 + t382 * t268 + m(5) * (t192 * t193 - t397 * t77 + t74 * t78) + m(6) * (t55 * t61 + t56 * t64 + t81 * t82) + m(7) * (t34 * t38 + t41 * t44 + t62 * t63) + (Ifges(4,4) * t273 - t382 * mrSges(4,2) + t336 * t320 + (t273 * t366 + t364) * t318) * t273 + (-Ifges(4,4) * t275 + t337 * t320 + t485 * t318 + (-Ifges(4,1) + Ifges(4,2) + t486) * t273) * t275;
t426 = t3 * qJD(1);
t153 = t355 * t275;
t154 = t360 * t275;
t155 = -mrSges(7,1) * t410 - mrSges(7,2) * t411;
t361 = t423 - t424;
t156 = t361 * t275;
t163 = t488 * t275;
t164 = t487 * t275;
t292 = Ifges(5,2) * t320 + t435;
t165 = t275 * t292;
t166 = -Ifges(7,1) * t411 + t256;
t167 = -Ifges(6,1) * t411 + t255;
t168 = t275 * t298;
t254 = Ifges(6,6) * t410;
t94 = t504 * t275;
t4 = t81 * t154 + t62 * t155 + t153 * t160 - t94 * t161 + t52 * t174 + t53 * t179 + t192 * t156 + t254 * t465 + t398 * t74 - t399 * t397 + m(6) * (t153 * t81 - t397 * t55 + t56 * t74) + m(7) * (t34 * t53 + t41 * t52 - t62 * t94) + ((t166 / 0.2e1 + t167 / 0.2e1 - t168 / 0.2e1 - t74 * mrSges(5,3) + t41 * mrSges(7,3) - t55 * mrSges(6,2) + t384 * t273 - t364) * t320 + (t163 / 0.2e1 + t164 / 0.2e1 + t165 / 0.2e1 - t56 * mrSges(6,2) - t397 * mrSges(5,3) + t34 * mrSges(7,3) + t336) * t318) * t275;
t420 = t4 * qJD(1);
t373 = -t182 / 0.2e1 - t175 / 0.2e1;
t374 = -t181 / 0.2e1 + t180 / 0.2e1;
t6 = t374 * t318 + t373 * t320 + t483 + t511;
t419 = t6 * qJD(1);
t13 = t400 * t410 + (t174 + t182) * t273 + m(7) * (t273 * t41 + t410 * t62) + m(6) * (t273 * t55 - t410 * t81);
t417 = qJD(1) * t13;
t21 = (m(7) * (t41 * t318 - t34 * t320) + t318 * t174 - t320 * t179) * t275;
t416 = qJD(1) * t21;
t300 = t469 * t318;
t72 = t273 * t300;
t404 = t72 * qJD(1);
t86 = 0.2e1 * (t314 / 0.4e1 + t315 / 0.4e1 + 0.1e1 / 0.4e1) * t450;
t403 = t86 * qJD(1);
t393 = qJD(3) * t318;
t392 = qJD(4) * t318;
t391 = m(7) * t410;
t389 = t94 * t477;
t381 = -t414 / 0.2e1;
t379 = -t413 / 0.2e1;
t378 = t413 / 0.2e1;
t363 = t488 / 0.2e1 + t487 / 0.2e1 + t292 / 0.2e1;
t294 = Ifges(7,1) * t318 - t434;
t296 = Ifges(6,1) * t318 - t433;
t362 = -t296 / 0.2e1 - t298 / 0.2e1 - t294 / 0.2e1;
t322 = (t153 * t278 + t489 * t81) * t479 + (-t270 * t94 + t271 * t62 - t280 * t441 + t284 * t442) * t477 - pkin(3) * t156 / 0.2e1 + t153 * t462 + t192 * t287 / 0.2e1 + t270 * t155 / 0.2e1 + t271 * t161 / 0.2e1 + t278 * t154 / 0.2e1 - t280 * t174 / 0.2e1 + t489 * t508 + t179 * t459 + t62 * t286 / 0.2e1 + t81 * t285 / 0.2e1 - t94 * t460 + t493 * t273 / 0.4e1;
t323 = (-pkin(4) * t64 + qJ(5) * t61) * t480 + (qJ(5) * t44 - t38 * t468) * t478 + pkin(4) * t178 / 0.2e1 - t468 * t466 + mrSges(7,1) * t472 + t44 * t475 + t61 * t474 + mrSges(6,1) * t470 - t77 * mrSges(5,1) / 0.2e1 + t78 * t476;
t326 = (-t52 / 0.2e1 - t34 / 0.2e1) * mrSges(7,3) + (t56 / 0.2e1 - t397 / 0.2e1) * mrSges(6,2) + (t439 * t479 - t374) * pkin(8) + t104 / 0.4e1 + t106 / 0.4e1 + t108 / 0.4e1 - t163 / 0.4e1 - t164 / 0.4e1 - t165 / 0.4e1;
t327 = (t471 + t41 / 0.2e1) * mrSges(7,3) + (t74 / 0.2e1 - t55 / 0.2e1) * mrSges(6,2) + (t440 * t480 + t373) * pkin(8) + t100 / 0.4e1 - t102 / 0.4e1 + t166 / 0.4e1 + t167 / 0.4e1 - t168 / 0.4e1 + t98 / 0.4e1;
t333 = -t298 / 0.4e1 - t296 / 0.4e1 - t294 / 0.4e1 - t491 / 0.4e1 + t291 / 0.4e1 + t289 / 0.4e1 + t280 * t473;
t334 = -t292 / 0.4e1 - t487 / 0.4e1 - t488 / 0.4e1 + t299 / 0.4e1 + t494 / 0.4e1 + t492 / 0.4e1 + mrSges(7,3) * t459;
t335 = pkin(8) * t509;
t1 = t322 + ((-0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(7,6) + Ifges(6,6) / 0.2e1) * t273 + t333 * t275 + t327) * t318 + ((-0.3e1 / 0.4e1 * Ifges(7,5) + t385) * t273 + t334 * t275 + t326) * t320 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t335) * t275 + t375 * qJ(5) + t323;
t20 = -pkin(3) * t287 + t271 * t282 - t489 * t360 + (t285 + t452) * t278 + (t286 + t451) * t270 + (t491 / 0.2e1 - t362 + t505) * t320 + (t492 / 0.2e1 + t494 / 0.2e1 + t299 / 0.2e1 - t363) * t318;
t353 = t1 * qJD(1) + t20 * qJD(3);
t328 = (-t318 * t81 + (-t412 + t446) * t320) * t480 + (t270 * t410 + t273 * t284 + t318 * t62) * t478;
t347 = m(6) * t470 + m(7) * t472;
t10 = (t508 - t161 / 0.2e1) * t318 + t445 * t413 + (t320 * t372 + t387) * t275 + t328 + t347;
t83 = -t270 * t448 + (m(6) * t278 - t282 - t360) * t318;
t352 = -qJD(1) * t10 - qJD(3) * t83;
t114 = 0.2e1 * (-t277 / 0.4e1 - t271 / 0.4e1) * m(7) - t286;
t340 = t504 * t450;
t32 = -t389 - t340 / 0.2e1 + t155;
t350 = qJD(1) * t32 - qJD(3) * t114;
t126 = m(7) * t351 + t495;
t338 = (-t304 / 0.2e1 + t425 / 0.2e1) * t273 + t63 * t477;
t14 = t338 + t510;
t349 = -qJD(1) * t14 + qJD(3) * t126;
t341 = 0.2e1 * t266 + t74;
t329 = t444 * t273 + t341 * t479 + (t249 + t341) * t477;
t346 = m(7) * t471 + t480 * t74;
t19 = t329 + t346;
t301 = qJ(5) * t469 + t444;
t348 = qJD(1) * t19 + qJD(4) * t301;
t345 = -mrSges(6,2) * pkin(4) + mrSges(7,3) * t468 - Ifges(7,5);
t344 = qJ(5) * t445 - t496;
t209 = (qJD(1) * t410 + t393) * m(7);
t191 = m(7) * t284 + (m(6) * pkin(8) - t445) * t320;
t170 = -t449 / 0.2e1 + t451 / 0.2e1;
t85 = (-0.1e1 / 0.2e1 + t394 / 0.2e1) * t450;
t71 = (t477 + t479) * t414 + t469 * t381;
t43 = t340 / 0.2e1 - t389;
t17 = t251 + t329 - t346 - t388;
t15 = t338 - t510;
t11 = mrSges(6,2) * t378 + mrSges(7,3) * t379 + t160 * t456 + t161 * t455 + t387 * t275 + t386 * t413 - t328 + t347 + (t360 / 0.2e1 + t460) * t410;
t9 = t324 - t325 + t465 * t495 + (-t361 / 0.2e1 + t372) * t275 + t273 * t509;
t7 = t180 * t456 + t181 * t455 + t399 * t507 + t483 - t511;
t2 = t322 + Ifges(7,5) * t378 + Ifges(6,6) * t381 + (-t427 / 0.4e1 + t326) * t320 + ((-Ifges(5,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t273 + t327) * t318 + (t318 * t333 + t320 * t334 + t335) * t275 - t323 + (t172 + t183) * qJ(5) / 0.2e1 + t496 * t414 / 0.2e1 + t497 * t379 + t486 * t501;
t12 = [qJD(2) * t5 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t13 + qJD(6) * t21, qJD(3) * t9 + qJD(4) * t7 + qJD(5) * t71 + qJD(6) * t85 + t429, t426 + t9 * qJD(2) + t2 * qJD(4) + t11 * qJD(5) + t15 * qJD(6) + (-t77 * mrSges(5,3) + t64 * mrSges(6,2) - t38 * mrSges(7,3) + t193 * mrSges(5,2) + (-t177 + t178) * pkin(8) + t337) * t393 + (-t193 * mrSges(4,1) + t192 * mrSges(4,2) - Ifges(4,6) * t275 + pkin(3) * t159 - t278 * t157 - t270 * t158 + t284 * t172 + t280 * t176 - t82 * t360 + t63 * t282 + (t78 * mrSges(5,3) + t61 * mrSges(6,2) - t44 * mrSges(7,3) - t193 * mrSges(5,1) + (t173 + t183) * pkin(8) - t485) * t320 + 0.2e1 * (t278 * t82 + (t318 * t64 + t320 * t61) * pkin(8)) * t479 + 0.2e1 * (-pkin(3) * t193 + (-t318 * t77 + t320 * t78) * pkin(8)) * t481 + 0.2e1 * (t270 * t63 + t280 * t38 + t284 * t44) * t477 + (t318 * t363 + t320 * t362 - Ifges(4,5)) * t273) * qJD(3), t420 + t7 * qJD(2) + t2 * qJD(3) + t17 * qJD(5) + t43 * qJD(6) + ((qJ(5) * t52 - t468 * t53) * t477 + (-pkin(4) * t74 - qJ(5) * t397) * t479) * t482 + (-t53 * mrSges(7,1) + t52 * mrSges(7,2) + t254 + (t344 * t320 + (-t345 - t497) * t318) * t275 + t498 * t74 - (-mrSges(5,2) + mrSges(6,3)) * t397) * qJD(4), qJD(2) * t71 + qJD(3) * t11 + qJD(4) * t17 + t417, qJD(2) * t85 + qJD(3) * t15 + qJD(4) * t43 + t416; qJD(3) * t8 - qJD(4) * t6 + qJD(5) * t72 + qJD(6) * t86 - t429, 0, t428, -t419 + (t304 + t421 - t422) * qJD(4) + t300 * qJD(5) + (-mrSges(7,1) + t498) * t392 + (-t452 / 0.2e1 + t449 / 0.2e1) * t482, qJD(4) * t300 + t404, t403; -qJD(2) * t8 + qJD(4) * t1 - qJD(5) * t10 - qJD(6) * t14 - t426, -t428, qJD(4) * t20 - qJD(5) * t83 + qJD(6) * t126, t191 * qJD(5) + t170 * qJD(6) + t344 * t392 + t353 + (-t284 * mrSges(7,1) - t280 * mrSges(7,2) + m(7) * (-qJ(5) * t280 - t284 * t468) + t345 * t320 + (-m(6) * t355 - t360 - t361) * pkin(8) + t493) * qJD(4), qJD(4) * t191 + t352, qJD(4) * t170 + t349; qJD(2) * t6 - qJD(3) * t1 + qJD(5) * t19 - qJD(6) * t32 - t420, t419, t114 * qJD(6) - t353, t301 * qJD(5), t348, -t350; -qJD(2) * t72 + qJD(3) * t10 - qJD(4) * t19 - qJD(6) * t391 - t417, -t404, -qJD(6) * t448 - t352, -t348, 0, -t209; -qJD(2) * t86 + qJD(3) * t14 + qJD(4) * t32 + qJD(5) * t391 - t416, -t403, -qJD(4) * t114 + qJD(5) * t448 - t349, t350, t209, 0;];
Cq  = t12;
