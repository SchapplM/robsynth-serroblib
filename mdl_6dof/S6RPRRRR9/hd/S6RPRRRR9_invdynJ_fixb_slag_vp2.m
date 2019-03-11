% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:09
% EndTime: 2019-03-09 07:23:59
% DurationCPUTime: 32.24s
% Computational Cost: add. (15858->859), mult. (31637->1151), div. (0->0), fcn. (21301->14), ass. (0->389)
t569 = -Ifges(4,2) / 0.2e1;
t316 = sin(qJ(3));
t321 = cos(qJ(3));
t365 = pkin(3) * t321 + pkin(8) * t316;
t268 = t365 * qJD(1);
t324 = -pkin(1) - pkin(7);
t289 = qJD(1) * t324 + qJD(2);
t320 = cos(qJ(4));
t315 = sin(qJ(4));
t425 = t315 * t321;
t177 = t320 * t268 - t289 * t425;
t323 = -pkin(9) - pkin(8);
t384 = qJD(4) * t323;
t422 = t316 * t320;
t392 = pkin(9) * t422;
t568 = -(pkin(4) * t321 + t392) * qJD(1) - t177 + t320 * t384;
t416 = t320 * t321;
t178 = t315 * t268 + t289 * t416;
t303 = t316 * qJD(1);
t383 = t315 * t303;
t567 = pkin(9) * t383 - t315 * t384 + t178;
t292 = t303 + qJD(4);
t287 = qJD(5) + t292;
t466 = -t287 / 0.2e1;
t279 = qJD(6) + t287;
t468 = -t279 / 0.2e1;
t406 = qJD(3) * t320;
t408 = qJD(1) * t321;
t257 = -t315 * t408 + t406;
t258 = qJD(3) * t315 + t320 * t408;
t314 = sin(qJ(5));
t319 = cos(qJ(5));
t169 = t257 * t314 + t258 * t319;
t475 = -t169 / 0.2e1;
t366 = t319 * t257 - t258 * t314;
t477 = -t366 / 0.2e1;
t313 = sin(qJ(6));
t318 = cos(qJ(6));
t100 = t169 * t318 + t313 * t366;
t481 = -t100 / 0.2e1;
t367 = -t169 * t313 + t318 * t366;
t484 = -t367 / 0.2e1;
t566 = Ifges(6,5) * t475 + Ifges(7,5) * t481 + Ifges(6,6) * t477 + Ifges(7,6) * t484 + Ifges(6,3) * t466 + Ifges(7,3) * t468;
t565 = -t257 / 0.2e1;
t564 = -t258 / 0.2e1;
t563 = -t292 / 0.2e1;
t445 = Ifges(4,4) * t321;
t562 = t316 * t569 + t445 / 0.2e1;
t262 = t314 * t320 + t315 * t319;
t512 = qJD(4) + qJD(5);
t176 = t512 * t262;
t232 = t262 * qJD(1);
t205 = t316 * t232;
t561 = t176 + t205;
t283 = t323 * t315;
t284 = t323 * t320;
t399 = qJD(5) * t319;
t400 = qJD(5) * t314;
t531 = t283 * t399 + t284 * t400 + t314 * t568 - t319 * t567;
t185 = t314 * t283 - t319 * t284;
t530 = -qJD(5) * t185 + t314 * t567 + t319 * t568;
t364 = pkin(3) * t316 - pkin(8) * t321;
t277 = qJ(2) + t364;
t233 = t277 * qJD(1);
t274 = t316 * t289;
t243 = qJD(3) * pkin(8) + t274;
t150 = t320 * t233 - t243 * t315;
t124 = -pkin(9) * t258 + t150;
t114 = pkin(4) * t292 + t124;
t151 = t233 * t315 + t243 * t320;
t125 = pkin(9) * t257 + t151;
t394 = qJD(1) * qJD(3);
t272 = qJDD(1) * t321 - t316 * t394;
t161 = qJD(4) * t257 + qJDD(3) * t315 + t272 * t320;
t273 = -t316 * qJDD(1) - t321 * t394;
t254 = qJDD(4) - t273;
t395 = qJD(1) * qJD(2);
t290 = qJDD(1) * qJ(2) + t395;
t164 = -pkin(3) * t273 - pkin(8) * t272 + t290;
t285 = qJDD(1) * t324 + qJDD(2);
t405 = qJD(3) * t321;
t187 = t316 * t285 + t289 * t405;
t183 = qJDD(3) * pkin(8) + t187;
t71 = -qJD(4) * t151 + t320 * t164 - t183 * t315;
t47 = pkin(4) * t254 - pkin(9) * t161 + t71;
t162 = -qJD(4) * t258 + qJDD(3) * t320 - t272 * t315;
t402 = qJD(4) * t320;
t403 = qJD(4) * t315;
t70 = t315 * t164 + t320 * t183 + t233 * t402 - t243 * t403;
t55 = pkin(9) * t162 + t70;
t12 = t114 * t399 - t125 * t400 + t314 * t47 + t319 * t55;
t66 = -qJD(5) * t169 - t161 * t314 + t162 * t319;
t10 = pkin(10) * t66 + t12;
t534 = pkin(10) * t366;
t119 = t319 * t125;
t61 = t114 * t314 + t119;
t52 = t61 + t534;
t437 = t313 * t52;
t548 = pkin(10) * t169;
t117 = t314 * t125;
t60 = t319 * t114 - t117;
t51 = t60 - t548;
t46 = pkin(5) * t287 + t51;
t16 = t318 * t46 - t437;
t13 = -qJD(5) * t61 - t314 * t55 + t319 * t47;
t247 = qJDD(5) + t254;
t65 = qJD(5) * t366 + t161 * t319 + t162 * t314;
t9 = pkin(5) * t247 - pkin(10) * t65 + t13;
t2 = qJD(6) * t16 + t10 * t318 + t313 * t9;
t435 = t318 * t52;
t17 = t313 * t46 + t435;
t3 = -qJD(6) * t17 - t10 * t313 + t318 * t9;
t560 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t312 = qJ(4) + qJ(5);
t305 = cos(t312);
t456 = pkin(4) * t320;
t276 = pkin(5) * t305 + t456;
t299 = pkin(3) + t456;
t559 = -m(6) * t299 - m(7) * (pkin(3) + t276);
t558 = t71 * mrSges(5,1) - t70 * mrSges(5,2);
t557 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t556 = -m(6) * t323 - m(7) * (-pkin(10) + t323) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t555 = t17 * mrSges(7,2) + t61 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t562 + Ifges(5,5) * t564 + Ifges(5,6) * t565 + Ifges(5,3) * t563 - t16 * mrSges(7,1) - t60 * mrSges(6,1) + t566;
t554 = t321 / 0.2e1;
t345 = t314 * t315 - t319 * t320;
t523 = t345 * t316;
t206 = qJD(1) * t523;
t544 = t512 * t345;
t553 = -pkin(5) * t408 + t530 + (t206 + t544) * pkin(10);
t552 = pkin(10) * t561 - t531;
t217 = t262 * t321;
t521 = t345 * qJD(1) - qJD(3) * t217 + t316 * t544;
t219 = t345 * t321;
t520 = -qJD(3) * t219 - t176 * t316 - t232;
t317 = sin(qJ(1));
t322 = cos(qJ(1));
t551 = g(1) * t317 - g(2) * t322;
t429 = t289 * t321;
t244 = -qJD(3) * pkin(3) - t429;
t181 = -pkin(4) * t257 + t244;
t108 = -pkin(5) * t366 + t181;
t22 = qJD(6) * t367 + t313 * t66 + t318 * t65;
t227 = qJDD(6) + t247;
t23 = -qJD(6) * t100 - t313 * t65 + t318 * t66;
t391 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t227;
t344 = t391 + t560;
t390 = Ifges(6,5) * t65 + Ifges(6,6) * t66 + Ifges(6,3) * t247;
t460 = mrSges(7,3) * t17;
t461 = mrSges(7,3) * t16;
t93 = Ifges(7,4) * t367;
t50 = Ifges(7,1) * t100 + Ifges(7,5) * t279 + t93;
t493 = -t50 / 0.2e1;
t441 = Ifges(7,4) * t100;
t49 = Ifges(7,2) * t367 + Ifges(7,6) * t279 + t441;
t495 = -t49 / 0.2e1;
t550 = t344 + t390 - (mrSges(7,1) * t108 + Ifges(7,4) * t481 + Ifges(7,2) * t484 + Ifges(7,6) * t468 - t460 + t495) * t100 + (-mrSges(7,2) * t108 + Ifges(7,1) * t481 + Ifges(7,4) * t484 + Ifges(7,5) * t468 + t461 + t493) * t367 + t557;
t501 = m(6) * pkin(4);
t458 = pkin(4) * t315;
t549 = m(6) * t458;
t455 = pkin(5) * t169;
t92 = -mrSges(5,1) * t162 + mrSges(5,2) * t161;
t547 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t272 - t92;
t519 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t257 - mrSges(5,2) * t258 - mrSges(4,3) * t408;
t407 = qJD(3) * t316;
t382 = t315 * t407;
t401 = qJD(4) * t321;
t543 = t320 * t401 - t382;
t517 = -t274 + (t383 + t403) * pkin(4);
t249 = Ifges(5,4) * t257;
t145 = t258 * Ifges(5,1) + t292 * Ifges(5,5) + t249;
t446 = Ifges(4,4) * t316;
t357 = t321 * Ifges(4,1) - t446;
t542 = Ifges(4,5) * qJD(3) + qJD(1) * t357 + t320 * t145;
t186 = t285 * t321 - t289 * t407;
t347 = t186 * t321 + t187 * t316;
t541 = m(5) + m(6) + m(7);
t363 = mrSges(4,1) * t321 - mrSges(4,2) * t316;
t540 = (-Ifges(4,1) * t316 - t445) * t554 + qJ(2) * t363;
t538 = -m(4) - t541;
t304 = sin(t312);
t454 = pkin(5) * t304;
t275 = t454 + t458;
t537 = m(7) * t275 + mrSges(2,1) - mrSges(3,2) + mrSges(4,3);
t362 = mrSges(4,1) * t316 + mrSges(4,2) * t321;
t536 = -m(5) * t364 + t316 * t559 + t321 * t556 + mrSges(2,2) - mrSges(3,3) - t362;
t499 = t22 / 0.2e1;
t498 = t23 / 0.2e1;
t491 = t65 / 0.2e1;
t490 = t66 / 0.2e1;
t479 = t161 / 0.2e1;
t478 = t162 / 0.2e1;
t473 = t227 / 0.2e1;
t472 = t247 / 0.2e1;
t471 = t254 / 0.2e1;
t184 = t319 * t283 + t284 * t314;
t140 = -pkin(10) * t262 + t184;
t141 = -pkin(10) * t345 + t185;
t84 = t140 * t318 - t141 * t313;
t533 = qJD(6) * t84 + t313 * t553 - t318 * t552;
t85 = t140 * t313 + t141 * t318;
t532 = -qJD(6) * t85 + t313 * t552 + t318 * t553;
t216 = t262 * t316;
t136 = -t216 * t313 - t318 * t523;
t529 = -qJD(6) * t136 - t313 * t520 + t318 * t521;
t134 = -t216 * t318 + t313 * t523;
t528 = qJD(6) * t134 + t313 * t521 + t318 * t520;
t457 = pkin(4) * t319;
t298 = pkin(5) + t457;
t397 = qJD(6) * t318;
t398 = qJD(6) * t313;
t428 = t313 * t314;
t68 = -t124 * t314 - t119;
t56 = t68 - t534;
t69 = t319 * t124 - t117;
t57 = t69 - t548;
t527 = -t313 * t56 - t318 * t57 + t298 * t397 + (-t314 * t398 + (t318 * t319 - t428) * qJD(5)) * pkin(4);
t427 = t314 * t318;
t526 = t313 * t57 - t318 * t56 - t298 * t398 + (-t314 * t397 + (-t313 * t319 - t427) * qJD(5)) * pkin(4);
t525 = -t501 - mrSges(5,1);
t522 = pkin(5) * t561 + t517;
t253 = t320 * t277;
t371 = -t315 * t324 + pkin(4);
t165 = -pkin(9) * t416 + t316 * t371 + t253;
t420 = t316 * t324;
t286 = t320 * t420;
t196 = t315 * t277 + t286;
t180 = -pkin(9) * t425 + t196;
t103 = t314 * t165 + t319 * t180;
t308 = qJ(6) + t312;
t296 = sin(t308);
t297 = cos(t308);
t358 = -mrSges(7,1) * t296 - mrSges(7,2) * t297;
t518 = mrSges(6,1) * t304 + mrSges(6,2) * t305 - t358;
t421 = t316 * t322;
t212 = t304 * t421 + t305 * t317;
t213 = -t304 * t317 + t305 * t421;
t199 = t296 * t421 + t297 * t317;
t200 = -t296 * t317 + t297 * t421;
t412 = t199 * mrSges(7,1) + t200 * mrSges(7,2);
t516 = -t212 * mrSges(6,1) - t213 * mrSges(6,2) - t412;
t423 = t316 * t317;
t210 = -t304 * t423 + t305 * t322;
t211 = t304 * t322 + t305 * t423;
t197 = -t296 * t423 + t297 * t322;
t198 = t296 * t322 + t297 * t423;
t413 = t197 * mrSges(7,1) - t198 * mrSges(7,2);
t515 = -t210 * mrSges(6,1) + t211 * mrSges(6,2) - t413;
t120 = mrSges(5,1) * t254 - mrSges(5,3) * t161;
t121 = -mrSges(5,2) * t254 + mrSges(5,3) * t162;
t514 = -t315 * t120 + t320 * t121;
t513 = -t315 * t71 + t320 * t70;
t361 = -mrSges(5,1) * t320 + mrSges(5,2) * t315;
t511 = m(5) * pkin(3) + mrSges(6,1) * t305 + mrSges(7,1) * t297 - mrSges(6,2) * t304 - mrSges(7,2) * t296 - t361 - t559;
t510 = -m(5) * pkin(8) - t556;
t509 = -mrSges(6,1) * t181 + mrSges(6,3) * t61;
t508 = mrSges(6,2) * t181 - mrSges(6,3) * t60;
t504 = qJD(1) ^ 2;
t503 = Ifges(7,4) * t499 + Ifges(7,2) * t498 + Ifges(7,6) * t473;
t502 = Ifges(7,1) * t499 + Ifges(7,4) * t498 + Ifges(7,5) * t473;
t500 = m(7) * pkin(5);
t497 = Ifges(6,4) * t491 + Ifges(6,2) * t490 + Ifges(6,6) * t472;
t496 = Ifges(6,1) * t491 + Ifges(6,4) * t490 + Ifges(6,5) * t472;
t494 = t49 / 0.2e1;
t492 = t50 / 0.2e1;
t489 = Ifges(5,1) * t479 + Ifges(5,4) * t478 + Ifges(5,5) * t471;
t442 = Ifges(6,4) * t169;
t88 = Ifges(6,2) * t366 + Ifges(6,6) * t287 + t442;
t488 = -t88 / 0.2e1;
t487 = t88 / 0.2e1;
t163 = Ifges(6,4) * t366;
t89 = Ifges(6,1) * t169 + Ifges(6,5) * t287 + t163;
t486 = -t89 / 0.2e1;
t485 = t89 / 0.2e1;
t483 = t367 / 0.2e1;
t482 = -m(3) - m(4);
t480 = t100 / 0.2e1;
t476 = t366 / 0.2e1;
t474 = t169 / 0.2e1;
t469 = t258 / 0.2e1;
t467 = t279 / 0.2e1;
t465 = t287 / 0.2e1;
t459 = pkin(4) * t258;
t451 = g(3) * t321;
t448 = mrSges(6,3) * t366;
t447 = mrSges(6,3) * t169;
t444 = Ifges(5,4) * t315;
t443 = Ifges(5,4) * t320;
t440 = t150 * mrSges(5,3);
t439 = t151 * mrSges(5,3);
t438 = t258 * Ifges(5,4);
t424 = t315 * t322;
t419 = t317 * t320;
t415 = t320 * t322;
t414 = t321 * t324;
t409 = t322 * pkin(1) + t317 * qJ(2);
t404 = qJD(3) * t324;
t396 = qJDD(1) * mrSges(3,2);
t387 = t315 * t420;
t386 = Ifges(5,5) * t161 + Ifges(5,6) * t162 + Ifges(5,3) * t254;
t291 = t316 * t404;
t381 = t321 * t404;
t373 = -t401 / 0.2e1;
t370 = -t394 / 0.2e1;
t368 = (t290 + t395) * qJ(2);
t102 = t319 * t165 - t180 * t314;
t256 = pkin(4) * t425 - t414;
t360 = mrSges(5,1) * t315 + mrSges(5,2) * t320;
t356 = Ifges(5,1) * t320 - t444;
t355 = Ifges(5,1) * t315 + t443;
t353 = -Ifges(5,2) * t315 + t443;
t352 = Ifges(5,2) * t320 + t444;
t351 = -Ifges(4,5) * t316 - Ifges(4,6) * t321;
t350 = Ifges(5,5) * t320 - Ifges(5,6) * t315;
t349 = Ifges(5,5) * t315 + Ifges(5,6) * t320;
t74 = pkin(5) * t316 + pkin(10) * t219 + t102;
t79 = -pkin(10) * t217 + t103;
t37 = -t313 * t79 + t318 * t74;
t38 = t313 * t74 + t318 * t79;
t348 = t150 * t320 + t151 * t315;
t189 = -mrSges(5,2) * t292 + mrSges(5,3) * t257;
t190 = mrSges(5,1) * t292 - mrSges(5,3) * t258;
t346 = -t315 * t189 - t320 * t190;
t135 = -t217 * t318 + t219 * t313;
t137 = -t217 * t313 - t219 * t318;
t172 = -t262 * t313 - t318 * t345;
t173 = t262 * t318 - t313 * t345;
t236 = t315 * t421 + t419;
t234 = -t315 * t423 + t415;
t182 = -qJDD(3) * pkin(3) - t186;
t342 = t316 * (-Ifges(4,2) * t321 - t446);
t188 = pkin(4) * t543 + t291;
t255 = qJD(3) * t365 + qJD(2);
t122 = -qJD(4) * t387 + t315 * t255 + t277 * t402 + t320 * t381;
t105 = -pkin(9) * t543 + t122;
t225 = t320 * t255;
t94 = t225 + (-t286 + (pkin(9) * t321 - t277) * t315) * qJD(4) + (t321 * t371 + t392) * qJD(3);
t35 = t319 * t105 + t165 * t399 - t180 * t400 + t314 * t94;
t336 = t315 * t401 + t316 * t406;
t106 = -pkin(4) * t162 + t182;
t331 = Ifges(5,5) * t321 - t316 * t356;
t330 = Ifges(5,6) * t321 - t316 * t353;
t329 = Ifges(5,3) * t321 - t316 * t350;
t36 = -qJD(5) * t103 - t105 * t314 + t319 * t94;
t327 = -qJD(4) * t348 + t513;
t300 = -pkin(1) * qJDD(1) + qJDD(2);
t280 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t303;
t266 = t362 * qJD(1);
t248 = t360 * t321;
t237 = -t315 * t317 + t316 * t415;
t235 = t316 * t419 + t424;
t223 = pkin(4) * t427 + t298 * t313;
t222 = -pkin(4) * t428 + t298 * t318;
t221 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t273;
t208 = pkin(5) * t345 - t299;
t195 = t253 - t387;
t171 = pkin(5) * t217 + t256;
t144 = t257 * Ifges(5,2) + t292 * Ifges(5,6) + t438;
t132 = mrSges(6,1) * t287 - t447;
t131 = -mrSges(6,2) * t287 + t448;
t128 = t205 * t313 + t206 * t318;
t127 = t205 * t318 - t206 * t313;
t126 = t459 + t455;
t123 = -qJD(4) * t196 - t315 * t381 + t225;
t113 = t262 * t407 + t321 * t544;
t111 = qJD(3) * t523 - t176 * t321;
t104 = -mrSges(6,1) * t366 + mrSges(6,2) * t169;
t86 = -pkin(5) * t113 + t188;
t83 = mrSges(7,1) * t279 - mrSges(7,3) * t100;
t82 = -mrSges(7,2) * t279 + mrSges(7,3) * t367;
t77 = t161 * Ifges(5,4) + t162 * Ifges(5,2) + t254 * Ifges(5,6);
t73 = -qJD(6) * t173 - t176 * t318 + t313 * t544;
t72 = qJD(6) * t172 - t176 * t313 - t318 * t544;
t59 = -mrSges(6,2) * t247 + mrSges(6,3) * t66;
t58 = mrSges(6,1) * t247 - mrSges(6,3) * t65;
t53 = -mrSges(7,1) * t367 + mrSges(7,2) * t100;
t43 = -qJD(6) * t137 - t111 * t313 + t113 * t318;
t41 = qJD(6) * t135 + t111 * t318 + t113 * t313;
t39 = -pkin(5) * t66 + t106;
t34 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t27 = pkin(10) * t113 + t35;
t26 = pkin(5) * t405 - pkin(10) * t111 + t36;
t19 = t318 * t51 - t437;
t18 = -t313 * t51 - t435;
t15 = -mrSges(7,2) * t227 + mrSges(7,3) * t23;
t14 = mrSges(7,1) * t227 - mrSges(7,3) * t22;
t8 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t5 = -qJD(6) * t38 + t26 * t318 - t27 * t313;
t4 = qJD(6) * t37 + t26 * t313 + t27 * t318;
t1 = [t540 * t394 - t347 * mrSges(4,3) - t542 * t407 / 0.2e1 - t519 * t291 + (Ifges(7,5) * t137 + Ifges(7,6) * t135) * t473 + (Ifges(7,5) * t41 + Ifges(7,6) * t43) * t467 + (Ifges(6,4) * t111 + Ifges(6,2) * t113) * t476 + (-Ifges(6,4) * t219 - Ifges(6,2) * t217) * t490 + (-t111 * t60 + t113 * t61 - t12 * t217 + t13 * t219) * mrSges(6,3) + (t150 * t336 - t151 * t543 - t416 * t71 - t425 * t70) * mrSges(5,3) + t244 * (mrSges(5,1) * t543 - mrSges(5,2) * t336) + (Ifges(7,1) * t41 + Ifges(7,4) * t43) * t480 + (Ifges(7,1) * t137 + Ifges(7,4) * t135) * t499 + (t382 / 0.2e1 + t320 * t373) * t144 + (Ifges(4,1) * t272 + Ifges(4,4) * t273) * t554 + (t135 * t2 - t137 * t3 - t16 * t41 + t17 * t43) * mrSges(7,3) + t273 * t562 + (t362 + 0.2e1 * mrSges(3,3)) * t290 + (-m(5) * t182 * t324 + Ifges(4,5) * qJDD(3) + t350 * t471 + t353 * t478 + t356 * t479) * t321 + (-Ifges(4,4) * t272 / 0.2e1 + t273 * t569 + Ifges(7,3) * t473 + Ifges(5,6) * t478 + Ifges(5,5) * t479 + Ifges(5,3) * t471 + t391 / 0.2e1 + t390 / 0.2e1 + t386 / 0.2e1 + Ifges(7,6) * t498 + Ifges(7,5) * t499 - Ifges(4,6) * qJDD(3) + Ifges(6,6) * t490 + Ifges(6,5) * t491 + Ifges(6,3) * t472 + t557 + t558 + t560) * t316 + (Ifges(6,5) * t111 + Ifges(6,6) * t113) * t465 + (-Ifges(6,5) * t219 - Ifges(6,6) * t217) * t472 + t111 * t485 + t113 * t487 + t416 * t489 + t41 * t492 + (qJD(3) * t331 - t355 * t401) * t469 - pkin(1) * t396 + (t150 * mrSges(5,1) - t151 * mrSges(5,2) + Ifges(6,5) * t474 + Ifges(7,5) * t480 + Ifges(6,6) * t476 + Ifges(7,6) * t483 + Ifges(6,3) * t465 + Ifges(7,3) * t467 - t555) * t405 + t221 * t420 + t547 * t414 + (Ifges(7,4) * t41 + Ifges(7,2) * t43) * t483 + (Ifges(7,4) * t137 + Ifges(7,2) * t135) * t498 + m(7) * (t108 * t86 + t16 * t5 + t17 * t4 + t171 * t39 + t2 * t38 + t3 * t37) + m(6) * (t102 * t13 + t103 * t12 + t106 * t256 + t181 * t188 + t35 * t61 + t36 * t60) - t77 * t425 / 0.2e1 + t342 * t370 + (Ifges(6,1) * t111 + Ifges(6,4) * t113) * t474 + (-Ifges(6,1) * t219 - Ifges(6,4) * t217) * t491 + t43 * t494 - t219 * t496 - t217 * t497 + t137 * t502 + t135 * t503 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t257 * (qJD(3) * t330 - t352 * t401) / 0.2e1 + t292 * (qJD(3) * t329 - t349 * t401) / 0.2e1 + m(4) * (t324 * t347 + t368) + m(3) * (-pkin(1) * t300 + t368) + t300 * mrSges(3,2) + qJD(2) * t266 + qJ(2) * (-mrSges(4,1) * t273 + mrSges(4,2) * t272) + t256 * t34 + t182 * t248 + t195 * t120 + t196 * t121 + t181 * (-mrSges(6,1) * t113 + mrSges(6,2) * t111) + t188 * t104 + t122 * t189 + t123 * t190 + t171 * t8 + t39 * (-mrSges(7,1) * t135 + mrSges(7,2) * t137) + t35 * t131 + t36 * t132 + t103 * t59 + t108 * (-mrSges(7,1) * t43 + mrSges(7,2) * t41) + t102 * t58 + qJD(3) ^ 2 * t351 / 0.2e1 + m(5) * (t244 * t324 * t407 + t122 * t151 + t123 * t150 + t195 * t71 + t196 * t70) + t315 * t145 * t373 + t280 * t381 + (-t424 * t501 - m(3) * t409 - t235 * mrSges(5,1) - t211 * mrSges(6,1) - t198 * mrSges(7,1) - t234 * mrSges(5,2) - t210 * mrSges(6,2) - t197 * mrSges(7,2) + t538 * (t322 * pkin(7) + t409) - t537 * t322 + t536 * t317) * g(2) + t272 * t357 / 0.2e1 + t37 * t14 + t38 * t15 + (-t237 * mrSges(5,1) - t213 * mrSges(6,1) - t200 * mrSges(7,1) + t236 * mrSges(5,2) + t212 * mrSges(6,2) + t199 * mrSges(7,2) + (m(3) * pkin(1) + t324 * t538 + t537 + t549) * t317 + ((-m(3) + t538) * qJ(2) + t536) * t322) * g(1) + t4 * t82 + t5 * t83 + t86 * t53 + t106 * (mrSges(6,1) * t217 - mrSges(6,2) * t219); t396 + t134 * t14 + t136 * t15 - t216 * t58 - t523 * t59 + t529 * t83 + t528 * t82 + t521 * t132 + t520 * t131 + (qJ(2) * t482 - mrSges(3,3)) * t504 + (-t266 + t346) * qJD(1) + (-t34 - t8 + (t189 * t320 - t190 * t315 + t280) * qJD(3) + t547) * t321 + (t221 + t346 * qJD(4) + (t104 + t53 - t519) * qJD(3) + t514) * t316 + m(3) * t300 + m(4) * t347 - t551 * (-t482 + t541) + (t108 * t407 + t134 * t3 + t136 * t2 + t16 * t529 + t17 * t528 - t321 * t39) * m(7) + (-t106 * t321 - t12 * t523 - t13 * t216 + t181 * t407 + t520 * t61 + t521 * t60) * m(6) + ((-t182 + (-t150 * t315 + t151 * t320) * qJD(3)) * t321 + (qJD(3) * t244 + t327) * t316 - t348 * qJD(1)) * m(5); (t342 / 0.2e1 - t540) * t504 + t292 * (t244 * t360 - t315 * t144 / 0.2e1) + t542 * t303 / 0.2e1 - (Ifges(6,1) * t474 + Ifges(6,4) * t476 + Ifges(6,5) * t465 + t485 + t508) * t544 + (-t151 * (mrSges(5,3) * t315 * t316 - mrSges(5,2) * t321) - t150 * (mrSges(5,1) * t321 + mrSges(5,3) * t422)) * qJD(1) + (Ifges(6,5) * t206 + Ifges(6,6) * t205) * t466 + (Ifges(7,5) * t128 + Ifges(7,6) * t127) * t468 - (Ifges(6,4) * t474 + Ifges(6,2) * t476 + Ifges(6,6) * t465 + t487 + t509) * t176 + (Ifges(7,4) * t480 + Ifges(7,2) * t483 + Ifges(7,6) * t467 + t460 + t494) * t73 + (t555 + t566) * t408 + (Ifges(6,4) * t206 + Ifges(6,2) * t205) * t477 + t206 * t486 + t205 * t488 + t315 * t489 + t128 * t493 + (Ifges(7,5) * t173 + Ifges(7,6) * t172) * t473 + t352 * t478 + t355 * t479 + t349 * t471 + (-t127 * t17 + t128 * t16 + t172 * t2 - t173 * t3) * mrSges(7,3) + (-t440 + t145 / 0.2e1) * t402 - t403 * t439 - t280 * t429 + (-pkin(3) * t182 - t150 * t177 - t151 * t178 - t244 * t274) * m(5) + (-t12 * t345 - t13 * t262 - t205 * t61 + t206 * t60) * mrSges(6,3) + (Ifges(6,4) * t262 - Ifges(6,2) * t345) * t490 + (Ifges(6,1) * t262 - Ifges(6,4) * t345) * t491 + (Ifges(6,5) * t262 - Ifges(6,6) * t345) * t472 + t106 * (mrSges(6,1) * t345 + mrSges(6,2) * t262) - t345 * t497 + (t316 * t511 + t321 * t510 + t362) * g(3) + (t257 * t353 + t258 * t356 + t292 * t350) * qJD(4) / 0.2e1 - (t257 * t330 + t258 * t331 + t292 * t329) * qJD(1) / 0.2e1 + t351 * t370 + ((-t128 + t72) * mrSges(7,2) + (t127 - t73) * mrSges(7,1)) * t108 + (Ifges(7,1) * t128 + Ifges(7,4) * t127) * t481 + t127 * t495 + t262 * t496 + (Ifges(7,4) * t173 + Ifges(7,2) * t172) * t498 + (Ifges(7,1) * t173 + Ifges(7,4) * t172) * t499 + t173 * t502 + t172 * t503 + t320 * t77 / 0.2e1 + Ifges(4,3) * qJDD(3) - t299 * t34 + Ifges(4,5) * t272 + Ifges(4,6) * t273 - t181 * (-mrSges(6,1) * t205 + mrSges(6,2) * t206) + t208 * t8 + t184 * t58 + t185 * t59 + t186 * mrSges(4,1) - t187 * mrSges(4,2) - t178 * t189 - t177 * t190 + t39 * (-mrSges(7,1) * t172 + mrSges(7,2) * t173) - pkin(3) * t92 + (Ifges(7,4) * t128 + Ifges(7,2) * t127) * t484 + t551 * (t316 * t510 - t321 * t511 - t363) + t513 * mrSges(5,3) + (m(5) * t327 - t189 * t403 - t190 * t402 + t514) * pkin(8) + t182 * t361 + t517 * t104 + t519 * t274 + t522 * t53 + t530 * t132 + t531 * t131 + (-t106 * t299 + t12 * t185 + t13 * t184 + t181 * t517 + t530 * t60 + t531 * t61) * m(6) + t532 * t83 + t533 * t82 + (t108 * t522 + t16 * t532 + t17 * t533 + t2 * t85 + t208 * t39 + t3 * t84) * m(7) + t84 * t14 + t85 * t15 + (Ifges(7,1) * t480 + Ifges(7,4) * t483 + Ifges(7,5) * t467 - t461 + t492) * t72 + (Ifges(6,1) * t206 + Ifges(6,4) * t205) * t475; (Ifges(6,1) * t475 + Ifges(6,4) * t477 + Ifges(6,5) * t466 + t486 - t508) * t366 + (t131 * t399 - t132 * t400 + t314 * t59) * pkin(4) + (Ifges(5,5) * t257 - Ifges(5,6) * t258) * t563 + (Ifges(5,1) * t257 - t438) * t564 + (-Ifges(5,2) * t258 + t145 + t249) * t565 + t58 * t457 + t144 * t469 + t258 * t439 + t257 * t440 + t550 + t386 - (Ifges(6,4) * t475 + Ifges(6,2) * t477 + Ifges(6,6) * t466 + t488 - t509) * t169 + t558 - t104 * t459 - m(6) * (t181 * t459 + t60 * t68 + t61 * t69) + (t12 * t314 + t13 * t319 + (-t314 * t60 + t319 * t61) * qJD(5)) * t501 - t244 * (mrSges(5,1) * t258 + mrSges(5,2) * t257) + g(3) * t248 + t222 * t14 + t223 * t15 - t150 * t189 + t151 * t190 - t69 * t131 - t68 * t132 - t126 * t53 + (t518 + t549) * t451 + (-m(7) * (-t275 * t423 + t276 * t322) + mrSges(5,2) * t235 + t525 * t234 + t515) * g(1) + (-m(7) * (t275 * t421 + t276 * t317) - mrSges(5,2) * t237 + t525 * t236 + t516) * g(2) + t526 * t83 + t527 * t82 + (-t108 * t126 + t16 * t526 + t17 * t527 + t2 * t223 + t222 * t3 + t275 * t451) * m(7); t88 * t474 + (Ifges(6,1) * t366 - t442) * t475 + (Ifges(6,5) * t366 - Ifges(6,6) * t169) * t466 - t53 * t455 - m(7) * (t108 * t455 + t16 * t18 + t17 * t19) + (t2 * t313 + t3 * t318 + (-t16 * t313 + t17 * t318) * qJD(6)) * t500 - t181 * (mrSges(6,1) * t169 + mrSges(6,2) * t366) - t19 * t82 - t18 * t83 + (t447 + t132) * t61 + (t448 - t131) * t60 + (-Ifges(6,2) * t169 + t163 + t89) * t477 + (m(7) * t454 + t518) * t451 + (-t212 * t500 + t516) * g(2) + (-t210 * t500 + t515) * g(1) + (t14 * t318 + t15 * t313 + t397 * t82 - t398 * t83) * pkin(5) + t550; -t108 * (mrSges(7,1) * t100 + mrSges(7,2) * t367) + (Ifges(7,1) * t367 - t441) * t481 + t49 * t480 + (Ifges(7,5) * t367 - Ifges(7,6) * t100) * t468 - t16 * t82 + t17 * t83 - g(1) * t413 - g(2) * t412 - t358 * t451 + (t100 * t17 + t16 * t367) * mrSges(7,3) + t344 + (-Ifges(7,2) * t100 + t50 + t93) * t484;];
tau  = t1;
