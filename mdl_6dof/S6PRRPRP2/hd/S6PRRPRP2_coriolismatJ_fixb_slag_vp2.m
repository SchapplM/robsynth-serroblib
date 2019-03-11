% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:41
% EndTime: 2019-03-08 21:28:53
% DurationCPUTime: 8.17s
% Computational Cost: add. (11811->556), mult. (26884->762), div. (0->0), fcn. (29391->10), ass. (0->280)
t451 = sin(pkin(11));
t400 = t451 * pkin(3);
t297 = t400 + pkin(9);
t314 = sin(qJ(5));
t308 = t314 ^ 2;
t317 = cos(qJ(5));
t310 = t317 ^ 2;
t575 = t308 + t310;
t564 = t575 * t297;
t525 = -m(7) / 0.2e1;
t356 = t314 * pkin(5) - qJ(6) * t317;
t567 = m(7) * t356;
t574 = -t567 / 0.2e1;
t558 = mrSges(7,2) + mrSges(6,3);
t523 = m(5) * pkin(3);
t411 = t523 / 0.2e1;
t573 = t451 * t411;
t318 = cos(qJ(3));
t315 = sin(qJ(3));
t374 = t451 * t315;
t452 = cos(pkin(11));
t274 = -t318 * t452 + t374;
t375 = t452 * t315;
t276 = -t318 * t451 - t375;
t428 = t276 * t314;
t204 = -mrSges(6,2) * t274 + mrSges(6,3) * t428;
t467 = t274 * mrSges(7,3);
t209 = mrSges(7,2) * t428 + t467;
t551 = t204 + t209;
t572 = t558 * (-t310 / 0.2e1 - t308 / 0.2e1);
t401 = t452 * pkin(3);
t298 = -t401 - pkin(4);
t526 = m(6) / 0.2e1;
t571 = -t298 * t526 + t452 * t411;
t561 = m(7) + m(6);
t313 = cos(pkin(6));
t312 = sin(pkin(6));
t316 = sin(qJ(2));
t421 = t312 * t316;
t268 = t313 * t318 - t315 * t421;
t420 = t312 * t318;
t269 = t313 * t315 + t316 * t420;
t174 = -t452 * t268 + t269 * t451;
t570 = t174 / 0.2e1;
t569 = -t551 / 0.2e1;
t299 = -pkin(3) * t318 - pkin(2);
t568 = m(5) * t299;
t468 = t274 * mrSges(5,3);
t566 = -t451 * t523 + mrSges(5,2);
t319 = cos(qJ(2));
t419 = t312 * t319;
t239 = t274 * t419;
t188 = -t239 * t314 - t317 * t421;
t189 = -t239 * t317 + t314 * t421;
t519 = mrSges(7,3) / 0.2e1;
t520 = -mrSges(6,2) / 0.2e1;
t407 = t520 + t519;
t524 = m(7) / 0.2e1;
t565 = (-pkin(5) * t188 + qJ(6) * t189) * t524 + t407 * t189;
t535 = t524 + t526;
t332 = t268 * t451 + t269 * t452;
t129 = t314 * t332 + t317 * t419;
t130 = -t314 * t419 + t317 * t332;
t357 = -t317 * pkin(5) - t314 * qJ(6);
t194 = t357 * t276;
t427 = t276 * t317;
t207 = mrSges(6,1) * t274 + mrSges(6,3) * t427;
t208 = -mrSges(7,1) * t274 - mrSges(7,2) * t427;
t379 = t208 / 0.2e1 - t207 / 0.2e1;
t201 = pkin(4) * t274 + pkin(9) * t276 + t299;
t417 = t314 * t201;
t481 = -qJ(4) - pkin(8);
t284 = t481 * t318;
t549 = -t452 * t284 + t481 * t374;
t436 = t549 * t317;
t107 = t417 + t436;
t433 = t274 * qJ(6);
t76 = t107 + t433;
t453 = t107 - t76;
t106 = t201 * t317 - t314 * t549;
t77 = -pkin(5) * t274 - t106;
t454 = t106 + t77;
t562 = (t129 * t453 + t130 * t454 + t174 * t194) * t525 - t379 * t130;
t560 = m(4) * t319;
t559 = mrSges(6,1) + mrSges(7,1);
t557 = Ifges(7,4) + Ifges(6,5);
t556 = Ifges(6,6) - Ifges(7,6);
t431 = t274 * t314;
t203 = mrSges(6,2) * t276 + mrSges(6,3) * t431;
t210 = mrSges(7,2) * t431 - mrSges(7,3) * t276;
t552 = t203 + t210;
t550 = t207 - t208;
t281 = -mrSges(7,1) * t317 - t314 * mrSges(7,3);
t282 = -mrSges(6,1) * t317 + t314 * mrSges(6,2);
t548 = t281 + t282;
t475 = Ifges(7,5) * t317;
t288 = Ifges(7,1) * t314 - t475;
t303 = Ifges(6,4) * t317;
t289 = Ifges(6,1) * t314 + t303;
t547 = t289 + t288;
t302 = Ifges(7,5) * t314;
t286 = -Ifges(7,3) * t317 + t302;
t373 = Ifges(7,1) * t317 + t302;
t546 = t315 ^ 2 + t318 ^ 2;
t495 = t276 / 0.2e1;
t545 = t495 * t564;
t476 = Ifges(6,4) * t314;
t361 = Ifges(6,1) * t317 - t476;
t543 = t373 + t361;
t285 = t315 * mrSges(4,1) + t318 * mrSges(4,2);
t483 = pkin(3) * t315;
t202 = -pkin(4) * t276 + pkin(9) * t274 + t483;
t227 = -t284 * t451 - t481 * t375;
t110 = t314 * t202 - t227 * t317;
t78 = -qJ(6) * t276 + t110;
t109 = t202 * t317 + t227 * t314;
t80 = pkin(5) * t276 - t109;
t539 = t314 * t80 + t317 * t78;
t359 = Ifges(6,2) * t314 - t303;
t538 = -mrSges(4,1) * t318 + mrSges(4,2) * t315;
t457 = t317 * mrSges(7,3);
t462 = t314 * mrSges(7,1);
t362 = -t457 + t462;
t195 = t362 * t274;
t458 = t317 * mrSges(6,2);
t463 = t314 * mrSges(6,1);
t363 = t458 + t463;
t196 = t363 * t274;
t537 = -t195 / 0.2e1 - t196 / 0.2e1;
t534 = -t359 + t547;
t300 = m(7) * qJ(6) + mrSges(7,3);
t273 = t357 + t298;
t533 = m(7) * t273 + t281;
t138 = t274 * Ifges(7,4) - t276 * t373;
t140 = t274 * Ifges(6,5) - t276 * t361;
t406 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t532 = -t406 * t274 - t138 / 0.2e1 - t140 / 0.2e1;
t486 = t317 / 0.2e1;
t489 = t314 / 0.2e1;
t490 = -t314 / 0.2e1;
t531 = t207 * t490 + t208 * t489 + t486 * t551;
t530 = m(6) * t298 - t452 * t523 - mrSges(5,1) + t282;
t113 = -t274 * t356 + t549;
t114 = -t276 * t356 + t227;
t198 = t363 * t276;
t430 = t274 * t317;
t205 = -mrSges(6,1) * t276 + mrSges(6,3) * t430;
t410 = mrSges(7,2) * t430;
t465 = t276 * mrSges(7,1);
t206 = -t410 + t465;
t350 = t106 * t314 - t107 * t317;
t353 = t314 * t77 + t317 * t76;
t501 = t210 / 0.2e1;
t197 = t362 * t276;
t508 = -t197 / 0.2e1;
t527 = m(5) / 0.2e1;
t529 = -t419 * t483 * t527 + (-t109 * t129 + t110 * t130 + t227 * t332 + (t350 + t549) * t174) * t526 + (t114 * t332 + t129 * t80 + t130 * t78 + (t113 - t353) * t174) * t524 + (-t205 / 0.2e1 + t206 / 0.2e1) * t129 + (t203 / 0.2e1 + t501) * t130 + (t508 - t198 / 0.2e1) * t332;
t528 = 0.2e1 * m(7);
t522 = mrSges(6,1) / 0.2e1;
t521 = -mrSges(7,1) / 0.2e1;
t518 = m(7) * t80;
t517 = t107 / 0.2e1;
t514 = -t174 / 0.2e1;
t499 = t273 / 0.2e1;
t498 = -t274 / 0.2e1;
t494 = -t281 / 0.2e1;
t493 = -t282 / 0.2e1;
t487 = -t317 / 0.2e1;
t472 = t114 * mrSges(7,1);
t471 = t114 * mrSges(7,3);
t470 = t227 * mrSges(6,1);
t469 = t227 * mrSges(6,2);
t466 = t274 * Ifges(6,6);
t464 = t276 * mrSges(5,3);
t212 = -t276 * mrSges(5,1) - t274 * mrSges(5,2);
t414 = t564 * t274;
t429 = t276 * t281;
t434 = t273 * t276;
t322 = -t414 * t526 + (-t414 - t434) * t524 - t429 / 0.2e1 + (t493 + t571) * t276 + (t572 - t573) * t274;
t323 = (t109 * t317 + t110 * t314) * t526 + (t314 * t78 - t317 * t80) * t524 + t205 * t486 + t206 * t487 + t315 * t411 + t552 * t489;
t12 = -t212 + t322 - t323;
t448 = t12 * qJD(2);
t447 = t130 * t317;
t418 = t314 * t174;
t14 = t561 * (-t129 * t418 + (t332 - t447) * t174);
t446 = t14 * qJD(1);
t326 = (t314 * t454 - t317 * t453) * t524 + t531;
t330 = -t463 / 0.2e1 - t462 / 0.2e1 - t458 / 0.2e1 + t457 / 0.2e1;
t16 = (t574 + t330) * t274 + t326;
t445 = t16 * qJD(2);
t238 = t276 * t419;
t111 = t174 * t238;
t424 = t312 ^ 2 * t316;
t18 = (-t268 * t312 * t315 + t269 * t420 - t424) * t560 + t561 * (t129 * t188 + t130 * t189 - t111) + (-t239 * t332 - t319 * t424 - t111) * m(5);
t442 = t18 * qJD(1);
t441 = t188 * t314;
t440 = t189 * t317;
t438 = t227 * t238;
t437 = t227 * t276;
t432 = t274 * t130;
t127 = t276 * t174;
t426 = t297 * t314;
t425 = t297 * t317;
t416 = t314 * t286;
t413 = qJD(5) * t314;
t191 = m(7) * t431;
t412 = t191 * qJD(2);
t409 = m(7) * t514;
t408 = -mrSges(6,1) / 0.2e1 + t521;
t405 = -Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t404 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t397 = t431 / 0.2e1;
t396 = -t430 / 0.2e1;
t393 = -t427 / 0.4e1;
t392 = t427 / 0.2e1;
t388 = -t419 / 0.2e1;
t383 = t281 * t486;
t263 = Ifges(7,5) * t427;
t134 = Ifges(7,6) * t274 - Ifges(7,3) * t428 - t263;
t136 = t276 * t359 + t466;
t382 = -t134 / 0.2e1 + t136 / 0.2e1;
t380 = t570 + t514;
t378 = t209 / 0.2e1 + t204 / 0.2e1;
t376 = mrSges(7,2) * qJ(6) + Ifges(6,6);
t370 = -mrSges(7,2) * pkin(5) + t557;
t365 = t197 + t198 + t464;
t287 = Ifges(6,2) * t317 + t476;
t358 = Ifges(7,3) * t314 + t475;
t320 = t285 * t388 + t558 * (t441 / 0.2e1 + t440 / 0.2e1) + t535 * (t440 + t441) * t297 + (mrSges(5,2) / 0.2e1 - t573) * t239 + (-t273 * t524 + mrSges(5,1) / 0.2e1 - t548 / 0.2e1 + t571) * t238;
t1 = -t320 + (t212 + t285) * t388 + t529 + (t537 - t531) * t174;
t133 = -Ifges(7,6) * t276 - t274 * t358;
t135 = -Ifges(6,6) * t276 + t274 * t359;
t137 = -Ifges(7,4) * t276 - t274 * t373;
t139 = -Ifges(6,5) * t276 - t274 * t361;
t213 = mrSges(5,1) * t274 - mrSges(5,2) * t276;
t3 = -pkin(2) * t285 + t106 * t205 + t107 * t203 + t109 * t207 + t110 * t204 - t113 * t197 - t114 * t195 - t227 * t196 - t549 * t198 + t77 * t206 + t80 * t208 + t78 * t209 + t76 * t210 + t299 * t212 + (-Ifges(4,4) * t315 + pkin(3) * t213) * t315 + t483 * t568 + m(6) * (t106 * t109 + t107 * t110 + t549 * t227) + m(7) * (t113 * t114 + t76 * t78 + t77 * t80) + (Ifges(5,4) * t274 + t532 * t317 + (-t404 * t274 + t382) * t314) * t274 + (Ifges(4,4) * t318 + (Ifges(4,1) - Ifges(4,2)) * t315) * t318 + (-Ifges(5,4) * t276 + (-t137 / 0.2e1 - t139 / 0.2e1 + t406 * t276) * t317 + (-t133 / 0.2e1 + t135 / 0.2e1 + t404 * t276) * t314 + (-Ifges(7,2) - Ifges(6,3) - Ifges(5,2) + Ifges(5,1)) * t274) * t276;
t355 = t1 * qJD(1) + t3 * qJD(2);
t262 = Ifges(7,6) * t427;
t6 = m(7) * (t106 * t76 + t107 * t77 + t114 * t194) - t194 * t197 + t106 * t209 + t262 * t498 - t107 * t207 + t107 * t208 + t106 * t204 + ((-t472 + t76 * mrSges(7,2) + t466 / 0.2e1 + t107 * mrSges(6,3) - t470 + t263 / 0.2e1 - Ifges(6,4) * t427 / 0.2e1 + t382) * t317 + (-t471 - t106 * mrSges(6,3) + t77 * mrSges(7,2) + t469 + (-Ifges(7,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t428 + (Ifges(7,3) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1) * t427 - t532) * t314) * t276;
t7 = t408 * t188 + t378 * t129 + (t558 * (t129 * t490 - t447 / 0.2e1) + (t493 + t494) * t174) * t276 + t562 + t565;
t354 = -t7 * qJD(1) + t6 * qJD(2);
t325 = (-t274 * t332 - t127) * t527 + t535 * (-t129 * t431 - t130 * t430 - t127);
t328 = t421 * t527 + t535 * (-t188 * t317 + t189 * t314);
t20 = t325 - t328;
t9 = t365 * t276 + (t550 * t314 - t551 * t317 + t468) * t274 + m(7) * (-t114 * t276 - t274 * t353) + m(6) * (t274 * t350 - t437) + m(5) * (-t274 * t549 - t437);
t352 = qJD(1) * t20 + qJD(2) * t9;
t30 = m(7) * (t114 * t427 + t274 * t76) - t197 * t427 + t274 * t209;
t41 = (t188 / 0.4e1 + t174 * t393 - t432 / 0.4e1) * t528;
t351 = -qJD(1) * t41 + qJD(2) * t30;
t349 = t114 * t356 + t194 * t273;
t39 = t467 + (t433 / 0.2e1 + t417 / 0.4e1 + t436 / 0.4e1 - t107 / 0.4e1) * t528;
t346 = qJD(2) * t39 + qJD(5) * t300;
t342 = Ifges(7,3) / 0.4e1 + Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1;
t339 = t194 * t494 - t356 * t508;
t336 = t282 * t495;
t335 = t429 / 0.2e1;
t327 = (t356 * t525 + t330) * t174;
t10 = (t408 * t314 + t407 * t317 + t574) * t174 - t327;
t36 = t356 * t281 + t298 * t363 + t287 * t490 + t416 / 0.2e1 + t358 * t487 + (t362 + t567) * t273 + t543 * t489 + t534 * t486;
t324 = (-pkin(5) * t80 + qJ(6) * t78) * t524 - pkin(5) * t206 / 0.2e1 + qJ(6) * t501 + t109 * t522 + t110 * t520 + t78 * t519 + t80 * t521;
t4 = t349 * t525 + (-t140 / 0.4e1 - t138 / 0.4e1 + t471 / 0.2e1 - t469 / 0.2e1 + (-0.3e1 / 0.4e1 * Ifges(7,4) - 0.3e1 / 0.4e1 * Ifges(6,5)) * t274 + (-t106 / 0.2e1 - t77 / 0.2e1) * mrSges(7,2) + (t454 * t525 - t379) * t297) * t317 + (t263 / 0.4e1 - t472 / 0.2e1 - t470 / 0.2e1 + t136 / 0.4e1 - t134 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t274 + (-t107 / 0.2e1 + t76 / 0.2e1) * mrSges(7,2) + (t453 * t525 + t378) * t297) * t314 + ((t302 / 0.4e1 - t287 / 0.4e1 + t286 / 0.4e1 + t298 * t522 + mrSges(7,1) * t499 - t342 * t317) * t317 + (-t289 / 0.4e1 - t288 / 0.4e1 - t303 / 0.4e1 + t298 * t520 + mrSges(7,3) * t499 + t342 * t314 + (-0.3e1 / 0.4e1 * Ifges(6,4) + Ifges(7,5) / 0.2e1) * t317) * t314 + t405 + t297 * t572) * t276 + t324 + t339;
t334 = t10 * qJD(1) + t4 * qJD(2) - t36 * qJD(3);
t211 = t533 * t314;
t329 = (-t114 * t314 + (t274 * t297 + t434) * t317) * t524 - t197 * t490;
t27 = t410 + (t383 + t521) * t276 - t518 / 0.2e1 + t329;
t60 = t380 * t314 * m(7);
t333 = -qJD(1) * t60 - qJD(2) * t27 + qJD(3) * t211;
t270 = (m(7) * t297 + mrSges(7,2)) * t317;
t61 = t314 * t409 - t418 * t524;
t42 = (t174 * t427 + t188 + t432) * t524;
t35 = (t107 + 0.2e1 * t433) * t524 + m(7) * t517 + t209;
t26 = t276 * t383 + t518 / 0.2e1 + t465 / 0.2e1 + t329;
t19 = t325 + t328;
t15 = t274 * t567 / 0.2e1 + mrSges(7,3) * t396 + mrSges(6,2) * t430 / 0.2e1 + t326 + t559 * t397;
t13 = t322 + t323;
t11 = -t356 * t409 - t327 + (t362 + t363) * t570;
t8 = t129 * t569 - t559 * t188 / 0.2e1 + (t336 + t335) * t174 + t558 * (t130 * t392 + t129 * t428 / 0.2e1) - t562 + t565;
t5 = t324 + t298 * t336 + t273 * t335 - t339 - t358 * t428 / 0.4e1 + t405 * t276 + (-t314 * t404 - t317 * t406) * t274 - t314 * t136 / 0.4e1 + t287 * t392 + ((t314 * t453 + t317 * t454) * t297 + t349) * t524 + t114 * t362 / 0.2e1 + t227 * t363 / 0.2e1 + (-t556 * t314 + t557 * t317) * t274 / 0.4e1 + (Ifges(7,1) * t428 + t134 - t263) * t314 / 0.4e1 + (t140 + t138) * t317 / 0.4e1 + t379 * t425 + t426 * t569 + t545 * mrSges(6,3) + (t534 + t289) * t428 / 0.4e1 + (0.2e1 * t286 + t543) * t393 + ((-t76 / 0.2e1 + t517) * t314 + t454 * t486 + t545) * mrSges(7,2);
t2 = t320 + (-t285 / 0.2e1 - t212 / 0.2e1) * t419 + t380 * t468 + t529 + (-t314 * t379 - t317 * t378 + t537) * t174;
t17 = [qJD(2) * t18 + qJD(3) * t14, t2 * qJD(3) + t19 * qJD(4) + t8 * qJD(5) + t42 * qJD(6) + t442 + (-t550 * t188 + t551 * t189 + t239 * t468 + t365 * t238 + 0.2e1 * (-t114 * t238 + t188 * t77 + t189 * t76) * t524 + 0.2e1 * (-t106 * t188 + t107 * t189 - t438) * t526 + 0.2e1 * (-t239 * t549 - t438) * t527 + ((t546 * mrSges(4,3) - mrSges(3,2)) * t319 + t546 * pkin(8) * t560 + (-m(4) * pkin(2) - mrSges(3,1) + t213 + t538 + t568) * t316) * t312) * qJD(2), t2 * qJD(2) + t11 * qJD(5) + t61 * qJD(6) + t446 + (-t269 * mrSges(4,1) - t268 * mrSges(4,2) + (t530 + t533) * t332 + (-t558 * t575 - t561 * t564 + t566) * t174) * qJD(3), qJD(2) * t19, t8 * qJD(2) + t11 * qJD(3) + (-t559 * t130 + (mrSges(6,2) - mrSges(7,3)) * t129) * qJD(5) + ((-pkin(5) * t130 - qJ(6) * t129) * qJD(5) / 0.2e1 + t130 * qJD(6) / 0.2e1) * t528, m(7) * t130 * qJD(5) + t42 * qJD(2) + t61 * qJD(3); qJD(3) * t1 + qJD(4) * t20 - qJD(5) * t7 - qJD(6) * t41 - t442, qJD(3) * t3 + qJD(4) * t9 + qJD(5) * t6 + qJD(6) * t30 (-(t557 * t314 + t556 * t317) * t276 / 0.2e1 + t552 * t425 + t566 * t227 + t547 * t396 + t538 * pkin(8) + m(7) * (t113 * t273 + t539 * t297) + t539 * mrSges(7,2) + t530 * t549 + (t139 + t137) * t489 + (-t205 + t206) * t426 + Ifges(4,5) * t318 - Ifges(4,6) * t315 - t298 * t196 + t113 * t281 - Ifges(5,5) * t274 + Ifges(5,6) * t276 - t273 * t195 + t287 * t397 + t400 * t464 + t401 * t468 + t135 * t486 + t133 * t487 + t416 * t498 + (m(6) * t297 + mrSges(6,3)) * (-t109 * t314 + t110 * t317)) * qJD(3) + t13 * qJD(4) + t5 * qJD(5) + t26 * qJD(6) + t355, qJD(3) * t13 + qJD(5) * t15 + t352, t5 * qJD(3) + t15 * qJD(4) + t35 * qJD(6) + t354 + (-t262 + (-m(7) * pkin(5) - t559) * t107 + (-mrSges(6,2) + t300) * t106 + (t314 * t370 + t317 * t376) * t276) * qJD(5), qJD(3) * t26 + qJD(5) * t35 + t351; -qJD(2) * t1 - qJD(5) * t10 + qJD(6) * t60 - t446, qJD(4) * t12 - qJD(5) * t4 + qJD(6) * t27 - t355, qJD(5) * t36 - qJD(6) * t211, t448, t270 * qJD(6) + (Ifges(7,6) - t376) * t413 - t334 + (t370 * t317 + (m(7) * t357 + t548) * t297) * qJD(5), qJD(5) * t270 - t333; -qJD(2) * t20, -qJD(3) * t12 + qJD(5) * t16 + qJD(6) * t191 - t352, -t448, 0, t445 + (t457 - t458 - t567) * qJD(5) + (m(7) * qJD(6) - t559 * qJD(5)) * t314, m(7) * t413 + t412; qJD(2) * t7 + qJD(3) * t10, qJD(3) * t4 - qJD(4) * t16 + qJD(6) * t39 - t354, t334, -t445, t300 * qJD(6), t346; qJD(2) * t41 - qJD(3) * t60, -qJD(3) * t27 - qJD(4) * t191 - qJD(5) * t39 - t351, t333, -t412, -t346, 0;];
Cq  = t17;
