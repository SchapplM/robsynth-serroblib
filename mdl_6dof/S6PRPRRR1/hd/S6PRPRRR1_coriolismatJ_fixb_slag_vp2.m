% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:34
% EndTime: 2019-03-08 20:22:48
% DurationCPUTime: 7.64s
% Computational Cost: add. (15427->454), mult. (35985->631), div. (0->0), fcn. (40334->12), ass. (0->295)
t411 = qJD(4) + qJD(5);
t450 = sin(pkin(12));
t451 = sin(pkin(6));
t366 = t451 * t450;
t452 = cos(pkin(12));
t367 = t452 * t451;
t490 = sin(qJ(2));
t492 = cos(qJ(2));
t248 = t366 * t490 - t367 * t492;
t304 = sin(qJ(4));
t306 = cos(qJ(4));
t489 = sin(qJ(5));
t491 = cos(qJ(5));
t274 = -t304 * t491 - t306 * t489;
t163 = t274 * t248;
t273 = t304 * t489 - t306 * t491;
t164 = t273 * t248;
t303 = sin(qJ(6));
t305 = cos(qJ(6));
t277 = -mrSges(7,1) * t305 + mrSges(7,2) * t303;
t575 = t164 * mrSges(6,2) / 0.2e1 - (t277 / 0.2e1 - mrSges(6,1) / 0.2e1) * t163;
t322 = t366 * t492 + t367 * t490;
t117 = t305 * t164 + t303 * t322;
t447 = t117 * t305;
t116 = -t303 * t164 + t305 * t322;
t448 = t116 * t303;
t574 = (t448 / 0.2e1 - t447 / 0.2e1) * mrSges(7,3) + t575;
t453 = cos(pkin(6));
t230 = t304 * t453 + t306 * t322;
t318 = -t304 * t322 + t306 * t453;
t139 = t230 * t489 - t491 * t318;
t314 = t230 * t491 + t318 * t489;
t90 = t248 * t303 + t305 * t314;
t457 = t305 * t90;
t89 = t248 * t305 - t303 * t314;
t460 = t303 * t89;
t352 = t314 - t457 + t460;
t572 = m(7) * t352 * t139;
t573 = qJD(1) * t572;
t458 = t305 * mrSges(7,2);
t462 = t303 * mrSges(7,1);
t278 = t458 + t462;
t218 = t273 * t278;
t426 = t274 * t305;
t228 = mrSges(7,1) * t273 + mrSges(7,3) * t426;
t427 = t274 * t303;
t226 = -mrSges(7,2) * t273 + mrSges(7,3) * t427;
t493 = t305 / 0.2e1;
t381 = t226 * t493;
t496 = -t303 / 0.2e1;
t571 = -t228 * t496 - t381 - t218 / 0.2e1;
t296 = Ifges(7,4) * t305;
t365 = Ifges(7,2) * t303 - t296;
t175 = -Ifges(7,6) * t274 + t273 * t365;
t481 = Ifges(7,4) * t303;
t284 = Ifges(7,1) * t305 - t481;
t177 = -Ifges(7,5) * t274 - t273 * t284;
t281 = Ifges(7,2) * t305 + t481;
t283 = Ifges(7,1) * t303 + t296;
t364 = Ifges(7,5) * t303 + Ifges(7,6) * t305;
t344 = t274 * t364;
t429 = t273 * t305;
t387 = -t429 / 0.2e1;
t430 = t273 * t303;
t388 = t430 / 0.2e1;
t495 = t303 / 0.2e1;
t341 = t175 * t493 + t177 * t495 + t281 * t388 + t283 * t387 - t344 / 0.2e1 + Ifges(6,6) * t274 - Ifges(6,5) * t273;
t549 = pkin(2) * t450;
t369 = pkin(8) + t549;
t342 = t304 * (-pkin(9) - t369);
t285 = t306 * t369;
t414 = t306 * pkin(9) + t285;
t517 = t489 * t342 + t414 * t491;
t541 = t517 * t277;
t545 = t517 * mrSges(6,1);
t212 = -t491 * t342 + t414 * t489;
t564 = t212 * mrSges(6,2);
t569 = t341 + t541 - t545 + t564;
t348 = t462 / 0.2e1 + t458 / 0.2e1;
t337 = t139 * t348;
t441 = t139 * t278;
t461 = t303 * mrSges(7,3);
t397 = t461 / 0.2e1;
t398 = -t461 / 0.2e1;
t538 = t397 + t398;
t559 = t337 + t441 / 0.2e1 + t538 * t90;
t568 = t559 * qJD(6);
t567 = -t541 / 0.2e1 + t545 / 0.2e1 - t564 / 0.2e1;
t548 = pkin(2) * t452;
t289 = -pkin(3) - t548;
t276 = -t306 * pkin(4) + t289;
t565 = m(6) * t276;
t563 = t212 * t303;
t562 = t212 * t305;
t528 = t212 * t314;
t561 = t212 * t489;
t204 = t273 * pkin(5) + t274 * pkin(10) + t276;
t104 = t204 * t305 - t303 * t517;
t105 = t204 * t303 + t305 * t517;
t357 = t104 * t303 - t105 * t305;
t512 = m(7) / 0.2e1;
t320 = (t517 + t357) * t512 + t571;
t560 = t320 * t273;
t437 = t517 * t212;
t525 = t314 * t277;
t529 = t314 * mrSges(6,1);
t546 = t139 * mrSges(6,2);
t558 = t525 - t529 + t546;
t557 = -t525 / 0.2e1 + t529 / 0.2e1 - t546 / 0.2e1;
t299 = t303 ^ 2;
t301 = t305 ^ 2;
t413 = t299 + t301;
t542 = t413 * t139;
t554 = -pkin(5) * t314 - pkin(10) * t542;
t485 = pkin(5) * t274;
t233 = pkin(10) * t273 - t485;
t487 = pkin(4) * t304;
t222 = t233 + t487;
t118 = t222 * t305 + t563;
t119 = t222 * t303 - t562;
t432 = t248 * t304;
t514 = m(6) / 0.2e1;
t225 = mrSges(7,2) * t274 + mrSges(7,3) * t430;
t227 = -mrSges(7,1) * t274 + mrSges(7,3) * t429;
t219 = t278 * t274;
t501 = -t219 / 0.2e1;
t506 = t90 / 0.2e1;
t518 = t314 * t501 + t89 * t227 / 0.2e1 + t225 * t506;
t553 = pkin(4) * t432 * t514 + (t118 * t89 + t119 * t90 + t139 * t517 + t528) * t512 + t518 + (t357 * t512 + t571) * t139;
t231 = -t274 * mrSges(6,1) - t273 * mrSges(6,2);
t295 = Ifges(7,5) * t305;
t479 = Ifges(7,6) * t303;
t351 = t295 / 0.2e1 - t479 / 0.2e1;
t178 = t273 * Ifges(7,5) - t274 * t284;
t420 = t305 * t178;
t176 = t273 * Ifges(7,6) + t274 * t365;
t422 = t303 * t176;
t494 = -t305 / 0.2e1;
t524 = t351 * t273;
t552 = (Ifges(6,4) * t273 - t420 / 0.2e1 + t422 / 0.2e1 - t524) * t273 - t517 * t219 + (t177 * t494 + t175 * t495 + (-Ifges(6,4) + t351) * t274 + (-Ifges(7,3) + Ifges(6,1) - Ifges(6,2)) * t273) * t274 + t104 * t227 + t105 * t225 - t212 * t218 + t276 * t231;
t407 = t489 * pkin(4);
t293 = t407 + pkin(10);
t408 = t491 * pkin(4);
t294 = -t408 - pkin(5);
t551 = -t293 * t542 + t294 * t314;
t550 = -t278 / 0.2e1;
t547 = pkin(5) * t517;
t543 = t294 * t517;
t540 = -t365 + t283;
t463 = t301 * mrSges(7,3);
t464 = t299 * mrSges(7,3);
t537 = -t463 - t464;
t534 = t284 / 0.4e1 - t281 / 0.4e1;
t527 = t273 * t352;
t373 = t281 * t496 + t284 * t495 + t540 * t493;
t423 = t294 * t278;
t133 = t373 + t423;
t220 = t274 * t281;
t221 = t283 * t274;
t478 = Ifges(7,3) * t274;
t395 = Ifges(7,5) * t387 + Ifges(7,6) * t388 - t478 / 0.2e1;
t507 = -mrSges(7,2) / 0.2e1;
t508 = mrSges(7,1) / 0.2e1;
t329 = t118 * t508 + t119 * t507 + t395;
t379 = t299 / 0.2e1 + t301 / 0.2e1;
t368 = t379 * mrSges(7,3);
t428 = t274 * t277;
t425 = t294 * t428;
t497 = t293 / 0.2e1;
t280 = t295 - t479;
t519 = -t273 * t280 / 0.4e1 + t212 * t550;
t15 = -t425 / 0.2e1 + (-t178 / 0.4e1 - t220 / 0.4e1 + t228 * t497) * t305 + (-t221 / 0.4e1 + t176 / 0.4e1 + t226 * t497) * t303 + (t534 * t305 + (-t283 / 0.4e1 + t365 / 0.4e1) * t303 - t293 * t368) * t274 + t329 + t519;
t31 = -t441 / 0.2e1 + t337;
t336 = t550 + t348;
t154 = t336 * t273;
t412 = t154 * qJD(3);
t523 = -t31 * qJD(1) - t15 * qJD(2) + t133 * qJD(4) - t412;
t419 = t305 * t225;
t421 = t303 * t227;
t522 = -t419 / 0.2e1 + t421 / 0.2e1;
t346 = t226 * t496 + t228 * t494;
t484 = pkin(5) * t278;
t146 = t373 - t484;
t323 = -t422 / 0.4e1 + t420 / 0.4e1 + t303 * t221 / 0.4e1 + t305 * t220 / 0.4e1 - t519 - t534 * t426 + t540 * t427 / 0.4e1 + t538 * t105;
t483 = mrSges(7,3) * t274;
t503 = -t428 / 0.2e1;
t313 = t323 + (t379 * t483 + t346) * pkin(10) + pkin(5) * t503;
t122 = t233 * t305 + t563;
t123 = t233 * t303 - t562;
t509 = -mrSges(7,1) / 0.2e1;
t349 = t122 * t509 + t123 * mrSges(7,2) / 0.2e1;
t18 = t313 + t524 + t478 / 0.2e1 + t349;
t33 = t336 * t139;
t371 = t491 * t507;
t85 = (-t294 / 0.2e1 + pkin(5) / 0.2e1) * t278 + (pkin(4) * t371 - t283 / 0.2e1 + t365 / 0.2e1) * t305 + (t408 * t509 - t284 / 0.2e1 + t281 / 0.2e1) * t303;
t521 = t33 * qJD(1) - t18 * qJD(2) + t85 * qJD(4) - t146 * qJD(5) + t412;
t396 = mrSges(7,3) * t493;
t516 = t116 * t398 + t117 * t396 - t575;
t513 = -m(7) / 0.2e1;
t511 = m(6) * pkin(4);
t510 = m(7) * pkin(4);
t500 = t248 / 0.2e1;
t499 = t273 / 0.2e1;
t456 = t306 * mrSges(5,2);
t279 = t304 * mrSges(5,1) + t456;
t498 = t279 / 0.2e1;
t377 = t413 * t273;
t488 = m(7) * (-pkin(10) * t377 + t485);
t486 = pkin(5) * t218;
t466 = t273 * mrSges(6,3);
t465 = t274 * mrSges(6,3);
t459 = t304 * mrSges(5,2);
t372 = (-0.1e1 + t413) * t274;
t67 = m(7) * t273 * t372;
t417 = t67 * qJD(3);
t455 = -t154 * qJD(6) - t417;
t155 = t273 * t348 + t278 * t499;
t454 = t155 * qJD(6) + t417;
t356 = -t447 + t448;
t439 = t163 * t273;
t43 = (t274 * t356 + t439) * t512 + (-t164 * t274 + t439) * t514;
t449 = qJD(1) * t43;
t446 = t118 * t303;
t445 = t119 * t305;
t442 = t139 * t163;
t435 = t212 * t163;
t424 = t294 * t218;
t45 = t428 * t499 + (-t274 * t368 - t346) * t274;
t418 = t45 * qJD(2);
t333 = t219 / 0.2e1 + t522;
t355 = t445 - t446;
t20 = t560 + (t333 + (-t212 - t355) * t512) * t274;
t354 = -t122 * t303 + t123 * t305;
t23 = ((-t212 - t354) * t512 + t333) * t274 + t560;
t410 = t20 * qJD(4) + t23 * qJD(5) + t45 * qJD(6);
t406 = mrSges(7,3) * t445;
t405 = t488 / 0.2e1;
t404 = t293 * t421;
t403 = t293 * t419;
t400 = -t464 / 0.2e1;
t399 = -t463 / 0.2e1;
t394 = t303 * t491;
t393 = t305 * t491;
t392 = t489 * t273;
t389 = t231 * t500;
t370 = -t394 / 0.2e1;
t24 = ((-t139 + t542) * t274 + t527) * t512;
t362 = t24 * qJD(1) + t20 * qJD(2);
t25 = (t139 * t372 + t527) * t512;
t361 = t25 * qJD(1) + t23 * qJD(2);
t203 = t248 * t322;
t21 = m(7) * (t116 * t89 + t117 * t90 + t442) + m(6) * (t164 * t314 + t203 + t442) + m(5) * (t203 + (-t230 * t306 + t304 * t318) * t248);
t360 = t21 * qJD(1) + t43 * qJD(3);
t359 = t24 * qJD(3) + t573;
t358 = t25 * qJD(3) + t573;
t350 = t116 * t508 + t117 * t507;
t343 = t537 * t273 - t231 - t428;
t340 = t413 * t491;
t312 = (t163 * t294 - t293 * t356) * t512 + (-t163 * t491 + t164 * t489) * t511 / 0.2e1 + mrSges(5,1) * t432 / 0.2e1 + t456 * t500;
t2 = t248 * t498 - t312 + t389 + t553 + t574;
t232 = mrSges(6,1) * t273 - mrSges(6,2) * t274;
t6 = t118 * t228 + t119 * t226 + t289 * t279 + (-Ifges(5,4) * t304 + pkin(4) * t232) * t304 + t487 * t565 + m(7) * (t104 * t118 + t105 * t119 + t437) + (Ifges(5,4) * t306 + (Ifges(5,1) - Ifges(5,2)) * t304) * t306 + t552;
t339 = t2 * qJD(1) + t6 * qJD(2) + t20 * qJD(3);
t310 = t320 * t139 + (t122 * t89 + t123 * t90 + t528) * t512 + t389 + t518;
t324 = m(7) * (-pkin(5) * t163 - pkin(10) * t356);
t3 = -t324 / 0.2e1 + t310 + t574;
t8 = t123 * t226 + t122 * t228 + m(7) * (t104 * t122 + t105 * t123 + t437) + t552;
t338 = t3 * qJD(1) + t8 * qJD(2) + t23 * qJD(3);
t335 = -t294 * t274 - t293 * t377;
t328 = t139 * t503 - t89 * t226 / 0.2e1 + t228 * t506;
t12 = (t460 / 0.2e1 - t457 / 0.2e1) * t483 + t328 + t350;
t14 = t104 * t226 - t105 * t228 + t212 * t428 + (-t357 * mrSges(7,3) + t176 * t493 + t221 * t494 + t364 * t499 + (t178 + t220) * t495) * t274;
t334 = -t12 * qJD(1) + t14 * qJD(2) + t45 * qJD(3);
t317 = (t277 - mrSges(6,1)) * t407 + (mrSges(7,3) * t413 - mrSges(6,2)) * t408;
t108 = (t293 * t340 + t294 * t489) * t510 + t317;
t309 = ((t139 * t489 + t393 * t90 - t394 * t89) * pkin(4) + t551) * t512 + t139 * t400 + t139 * t399 - t557;
t315 = t554 * t513 - (t400 + t399) * t139 + t557;
t11 = t309 + t315;
t316 = m(7) * ((-t274 * t340 + t392) * pkin(4) + t335);
t55 = t405 - t316 / 0.2e1;
t307 = (t543 + t354 * t293 + (-t104 * t394 + t105 * t393 + t561) * pkin(4)) * t512 - t424 / 0.2e1 + t122 * t398 + t123 * t396 - t404 / 0.2e1 + t403 / 0.2e1 + t407 * t501 + pkin(4) * t228 * t370 + t381 * t408 - t567;
t311 = -t513 * t547 - t486 / 0.2e1 + t118 * t397 - t406 / 0.2e1 + (t355 * t513 + t522) * pkin(10) + t567;
t9 = t307 + t311;
t327 = t11 * qJD(1) + t9 * qJD(2) - t55 * qJD(3) + t108 * qJD(4);
t86 = t423 / 0.2e1 - t484 / 0.2e1 + (mrSges(7,1) * t370 + t305 * t371) * pkin(4) + t373;
t49 = t316 / 0.2e1 + t405 + t343;
t17 = t313 - t349 + t395;
t16 = t323 + t425 / 0.2e1 + t329 + t413 * t483 * t497 + t346 * t293;
t13 = -t328 + t350 + (t396 * t90 + t398 * t89) * t274;
t10 = t309 - t315;
t7 = -t311 + t307 + t341;
t5 = t43 * qJD(2) + t24 * qJD(4) + t25 * qJD(5);
t4 = t324 / 0.2e1 + t310 + t516;
t1 = (t231 / 0.2e1 + t498) * t248 + t312 + t516 + t553;
t19 = [t21 * qJD(2) + t411 * t572 (t116 * t228 + t117 * t226 + m(7) * (t104 * t116 + t105 * t117 + t435) - t164 * t466 + m(6) * (t164 * t517 + t435) + (-mrSges(3,1) * t490 - mrSges(3,2) * t492) * t451 + (-m(4) * t548 + m(5) * t289 - mrSges(5,1) * t306 - mrSges(4,1) + t232 + t459 + t565) * t322 + (-t219 - t465) * t163 + (-m(4) * t549 + mrSges(4,2) + (m(5) * t369 + mrSges(5,3)) * (-t304 ^ 2 - t306 ^ 2)) * t248) * qJD(2) + t1 * qJD(4) + t4 * qJD(5) + t13 * qJD(6) + t360, t5, t1 * qJD(2) + (m(7) * t551 - t314 * t491 * t511 - t318 * mrSges(5,2) - t230 * mrSges(5,1) - (t489 * t511 - t537) * t139 + t558) * qJD(4) + t10 * qJD(5) + t568 + t359, t4 * qJD(2) + t10 * qJD(4) + (m(7) * t554 - mrSges(7,3) * t542 + t558) * qJD(5) + t568 + t358, t13 * qJD(2) + (-mrSges(7,1) * t90 - mrSges(7,2) * t89) * qJD(6) + t411 * t559; qJD(4) * t2 + qJD(5) * t3 - qJD(6) * t12 - t360, qJD(4) * t6 + qJD(5) * t8 + qJD(6) * t14, t410 - t449 ((-t491 * t517 - t561) * t511 + Ifges(5,5) * t306 - Ifges(5,6) * t304 - t424 + m(7) * (t293 * t355 + t543) - mrSges(5,1) * t285 + t403 + t406 + t369 * t459 - t404 - mrSges(7,3) * t446 + t408 * t466 + t407 * t465 + t569) * qJD(4) + t7 * qJD(5) + t16 * qJD(6) + t339, t7 * qJD(4) + (-m(7) * t547 + t354 * mrSges(7,3) + t486 + (m(7) * t354 + t419 - t421) * pkin(10) + t569) * qJD(5) + t17 * qJD(6) + t338, t16 * qJD(4) + t17 * qJD(5) + (-mrSges(7,1) * t105 - mrSges(7,2) * t104 + t344) * qJD(6) + t334; t5, t410 + t449, t411 * t67 ((t274 * t491 - t392) * t511 + m(7) * t335 + t343 - t279) * qJD(4) + t49 * qJD(5) + t362 + t454, t49 * qJD(4) + (t343 + t488) * qJD(5) + t361 + t454, -qJD(6) * t428 + t155 * t411 + t418; -qJD(2) * t2 + qJD(5) * t11 - qJD(6) * t31 - t359, qJD(5) * t9 - qJD(6) * t15 - t339, -qJD(5) * t55 - t362 + t455, qJD(5) * t108 + qJD(6) * t133 ((-pkin(5) * t489 + pkin(10) * t340) * t510 + t317) * qJD(5) + t86 * qJD(6) + t327, t86 * qJD(5) + (t277 * t293 + t280) * qJD(6) + t523; -qJD(2) * t3 - qJD(4) * t11 - qJD(6) * t33 - t358, -qJD(4) * t9 + qJD(6) * t18 - t338, qJD(4) * t55 - t361 + t455, -qJD(6) * t85 - t327, t146 * qJD(6) (pkin(10) * t277 + t280) * qJD(6) - t521; t12 * qJD(2) + t31 * qJD(4) + t33 * qJD(5), qJD(4) * t15 - qJD(5) * t18 - t334, t154 * t411 - t418, qJD(5) * t85 - t523, t521, 0;];
Cq  = t19;
