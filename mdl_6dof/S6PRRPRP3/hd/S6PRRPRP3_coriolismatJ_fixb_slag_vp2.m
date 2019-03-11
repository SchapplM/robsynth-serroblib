% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRP3
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:21
% EndTime: 2019-03-08 21:34:34
% DurationCPUTime: 7.74s
% Computational Cost: add. (12778->620), mult. (31015->858), div. (0->0), fcn. (32645->10), ass. (0->316)
t548 = m(7) / 0.2e1;
t575 = Ifges(7,4) + Ifges(6,5);
t587 = -Ifges(6,6) + Ifges(7,6);
t586 = 0.2e1 * t548;
t585 = mrSges(7,2) + mrSges(6,3);
t381 = sin(qJ(3));
t383 = cos(qJ(3));
t348 = -pkin(3) * t383 - qJ(4) * t381 - pkin(2);
t380 = cos(pkin(11));
t329 = t380 * t348;
t378 = sin(pkin(11));
t466 = t380 * t381;
t229 = -pkin(9) * t466 + t329 + (-pkin(8) * t378 - pkin(4)) * t383;
t465 = t380 * t383;
t273 = pkin(8) * t465 + t348 * t378;
t471 = t378 * t381;
t250 = -pkin(9) * t471 + t273;
t513 = sin(qJ(5));
t514 = cos(qJ(5));
t111 = t229 * t514 - t250 * t513;
t99 = pkin(5) * t383 - t111;
t584 = t111 + t99;
t401 = -t378 * t513 + t380 * t514;
t484 = qJ(6) * t401;
t402 = t378 * t514 + t380 * t513;
t512 = pkin(5) * t402;
t422 = t484 - t512;
t583 = t422 * t548;
t366 = -pkin(4) * t380 - pkin(3);
t421 = -pkin(5) * t401 - qJ(6) * t402;
t213 = t366 + t421;
t239 = -mrSges(7,1) * t401 - mrSges(7,3) * t402;
t582 = m(7) * t213 + t239;
t379 = sin(pkin(6));
t382 = sin(qJ(2));
t384 = cos(qJ(2));
t462 = t383 * t384;
t270 = (-t378 * t462 + t380 * t382) * t379;
t271 = (t378 * t382 + t380 * t462) * t379;
t162 = -t270 * t514 + t271 * t513;
t163 = t270 * t513 + t271 * t514;
t546 = -mrSges(6,2) / 0.2e1;
t454 = mrSges(7,3) / 0.2e1 + t546;
t581 = t454 * t163 + (-pkin(5) * t162 + qJ(6) * t163) * t548;
t580 = -Ifges(6,6) / 0.2e1;
t579 = m(7) + m(6);
t306 = t402 * t383;
t526 = t306 / 0.2e1;
t307 = t401 * t383;
t578 = -t307 / 0.2e1;
t576 = mrSges(6,1) + mrSges(7,1);
t574 = Ifges(7,2) + Ifges(6,3);
t304 = t401 * t381;
t573 = t304 * mrSges(6,3);
t305 = t402 * t381;
t572 = t305 * mrSges(6,3);
t545 = mrSges(6,3) / 0.2e1;
t455 = t545 + mrSges(7,2) / 0.2e1;
t571 = t305 * t455;
t448 = t513 * t229;
t450 = t514 * t250;
t112 = t450 + t448;
t464 = t383 * qJ(6);
t96 = t112 - t464;
t240 = -mrSges(6,1) * t401 + mrSges(6,2) * t402;
t569 = t239 + t240;
t258 = -mrSges(7,2) * t306 + mrSges(7,3) * t381;
t261 = -mrSges(6,2) * t381 - mrSges(6,3) * t306;
t568 = t258 + t261;
t488 = t383 * mrSges(7,3);
t259 = -t305 * mrSges(7,2) - t488;
t260 = mrSges(6,2) * t383 - t572;
t567 = t260 + t259;
t262 = -mrSges(6,1) * t383 - t573;
t263 = mrSges(7,1) * t383 + t304 * mrSges(7,2);
t566 = t262 - t263;
t564 = t401 * t575 + t402 * t587;
t563 = t304 * t587 - t305 * t575;
t352 = pkin(3) * t381 - qJ(4) * t383;
t287 = pkin(8) * t471 + t352 * t380;
t288 = -pkin(8) * t466 + t352 * t378;
t562 = -t287 * t378 + t288 * t380;
t437 = -t262 / 0.2e1 + t263 / 0.2e1;
t550 = m(6) / 0.2e1;
t560 = t550 + t548;
t559 = m(5) + t579;
t199 = mrSges(7,1) * t305 - mrSges(7,3) * t304;
t200 = mrSges(6,1) * t305 + mrSges(6,2) * t304;
t492 = t380 * mrSges(5,2);
t493 = t378 * mrSges(5,1);
t425 = t492 + t493;
t318 = t425 * t381;
t558 = t199 + t200 + t318;
t370 = m(7) * qJ(6) + mrSges(7,3);
t557 = -m(7) * pkin(5) - t576;
t556 = -mrSges(6,2) + t370;
t555 = 0.2e1 * m(7);
t375 = t380 ^ 2;
t554 = 2 * qJD(3);
t553 = -m(5) / 0.2e1;
t552 = m(5) / 0.2e1;
t551 = -m(6) / 0.2e1;
t549 = -m(7) / 0.2e1;
t547 = -mrSges(6,1) / 0.2e1;
t469 = t379 * t382;
t485 = cos(pkin(6));
t317 = t381 * t485 + t383 * t469;
t468 = t379 * t384;
t249 = t317 * t380 - t378 * t468;
t404 = t317 * t378 + t380 * t468;
t121 = t249 * t514 - t404 * t513;
t543 = t121 / 0.2e1;
t541 = Ifges(6,4) * t578 + Ifges(6,2) * t526 + t381 * t580;
t316 = t381 * t469 - t383 * t485;
t190 = t402 * t316;
t540 = t190 / 0.2e1;
t191 = t401 * t316;
t539 = -t191 / 0.2e1;
t473 = t305 * qJ(6);
t511 = t304 * pkin(5);
t423 = -t473 - t511;
t538 = -t423 / 0.2e1;
t197 = mrSges(7,1) * t304 + mrSges(7,3) * t305;
t537 = t197 / 0.2e1;
t198 = mrSges(6,1) * t304 - mrSges(6,2) * t305;
t536 = t198 / 0.2e1;
t535 = -t422 / 0.2e1;
t237 = mrSges(7,1) * t402 - mrSges(7,3) * t401;
t534 = t237 / 0.2e1;
t238 = mrSges(6,1) * t402 + mrSges(6,2) * t401;
t533 = t238 / 0.2e1;
t509 = pkin(9) + qJ(4);
t350 = t509 * t380;
t436 = t509 * t378;
t256 = t350 * t513 + t436 * t514;
t532 = -t256 / 0.2e1;
t491 = t381 * mrSges(7,1);
t496 = t307 * mrSges(7,2);
t265 = -t491 + t496;
t529 = t265 / 0.2e1;
t528 = -t304 / 0.2e1;
t527 = -t306 / 0.2e1;
t525 = t307 / 0.2e1;
t523 = t317 / 0.2e1;
t522 = -t401 / 0.2e1;
t520 = -t402 / 0.2e1;
t519 = -t378 / 0.2e1;
t518 = t380 / 0.2e1;
t517 = t381 / 0.2e1;
t372 = t381 * pkin(8);
t373 = t383 * pkin(8);
t508 = Ifges(5,4) * t378;
t507 = Ifges(5,4) * t380;
t506 = Ifges(6,4) * t304;
t505 = Ifges(6,4) * t402;
t504 = Ifges(5,5) * t380;
t503 = Ifges(7,5) * t305;
t502 = Ifges(7,5) * t401;
t501 = Ifges(5,2) * t378;
t500 = Ifges(5,6) * t378;
t499 = t306 * mrSges(6,1);
t498 = t306 * mrSges(7,1);
t497 = t307 * mrSges(6,2);
t495 = t307 * mrSges(7,3);
t494 = t401 * mrSges(7,2);
t490 = t381 * mrSges(4,2);
t489 = t383 * mrSges(4,2);
t349 = -mrSges(5,1) * t380 + mrSges(5,2) * t378;
t487 = -mrSges(4,1) + t349;
t486 = t112 - t96;
t120 = t249 * t513 + t404 * t514;
t451 = t381 * t468;
t267 = t316 * t451;
t14 = m(5) * (t249 * t271 - t270 * t404 + t267) + m(4) * (t316 * t381 + t317 * t383 - t469) * t468 + t579 * (t120 * t162 + t121 * t163 + t267);
t483 = t14 * qJD(1);
t231 = t316 * t317;
t395 = t249 * t380 + t378 * t404;
t15 = m(5) * (-t316 * t395 + t231) + t579 * (-t120 * t190 - t121 * t191 + t231);
t482 = t15 * qJD(1);
t478 = t270 * t378;
t477 = t271 * t380;
t474 = t304 * t316;
t339 = -mrSges(5,1) * t383 - mrSges(5,3) * t466;
t472 = t378 * t339;
t470 = t378 * t383;
t337 = mrSges(5,2) * t383 - mrSges(5,3) * t471;
t467 = t380 * t337;
t463 = t383 * t121;
t232 = pkin(4) * t381 - pkin(9) * t465 + t287;
t253 = -pkin(9) * t470 + t288;
t115 = t232 * t513 + t253 * t514;
t341 = pkin(4) * t471 + t372;
t342 = pkin(4) * t470 + t373;
t457 = t378 ^ 2 + t375;
t456 = -mrSges(7,1) / 0.2e1 + t547;
t453 = -Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t452 = Ifges(7,6) / 0.2e1 + t580;
t444 = t316 * t520;
t443 = -t468 / 0.2e1;
t442 = t468 / 0.2e1;
t441 = t383 * t522;
t297 = Ifges(7,5) * t304;
t182 = -Ifges(7,6) * t383 + Ifges(7,3) * t305 + t297;
t184 = -Ifges(6,2) * t305 - Ifges(6,6) * t383 + t506;
t440 = t182 / 0.2e1 - t184 / 0.2e1;
t186 = Ifges(7,1) * t304 - Ifges(7,4) * t383 + t503;
t300 = Ifges(6,4) * t305;
t188 = Ifges(6,1) * t304 - Ifges(6,5) * t383 - t300;
t439 = t186 / 0.2e1 + t188 / 0.2e1;
t438 = -t259 / 0.2e1 - t260 / 0.2e1;
t257 = t350 * t514 - t436 * t513;
t434 = -t190 * t256 - t191 * t257;
t432 = t256 * t304 - t257 * t305;
t429 = t381 * t442;
t428 = t455 * t304;
t424 = pkin(5) * t305 - qJ(6) * t304;
t105 = qJ(6) * t381 + t115;
t114 = t232 * t514 - t253 * t513;
t106 = -pkin(5) * t381 - t114;
t151 = t424 + t341;
t152 = pkin(5) * t306 - qJ(6) * t307 + t342;
t201 = -t495 + t498;
t202 = t497 + t499;
t264 = mrSges(6,1) * t381 - mrSges(6,3) * t307;
t319 = t425 * t383;
t338 = -mrSges(5,2) * t381 - mrSges(5,3) * t470;
t340 = mrSges(5,1) * t381 - mrSges(5,3) * t465;
t353 = t381 * mrSges(4,1) + t489;
t272 = -pkin(8) * t470 + t329;
t412 = -t272 * t378 + t273 * t380;
t385 = (t288 * t249 - t287 * t404 + t317 * t372 + (-t412 + t373) * t316) * t552 + (-t112 * t191 + t115 * t121 + t316 * t342 + t317 * t341) * t550 + (t105 * t121 + t151 * t317 + t152 * t316 - t191 * t96) * t548 + t249 * t338 / 0.2e1 - t404 * t340 / 0.2e1 - t316 * t467 / 0.2e1 + t353 * t443 + t568 * t543 + t567 * t539 - (-t111 * t550 + t548 * t99 + t437) * t190 + (-t114 * t550 + t106 * t548 - t264 / 0.2e1 + t529) * t120 + t558 * t523 + (t201 + t202 + t319 + t472) * t316 / 0.2e1;
t413 = t162 * t256 + t163 * t257;
t387 = (-pkin(3) * t451 + (t477 - t478) * qJ(4)) * t553 + (t366 * t451 + t413) * t551 + (t213 * t451 + t413) * t549 + mrSges(4,1) * t429 + t442 * t489 + (t478 / 0.2e1 - t477 / 0.2e1) * mrSges(5,3) + (t349 + t569) * t381 * t443 + t585 * (t162 * t520 + t163 * t522);
t4 = t387 + t385;
t183 = Ifges(7,5) * t307 + Ifges(7,6) * t381 + Ifges(7,3) * t306;
t187 = Ifges(7,1) * t307 + Ifges(7,4) * t381 + Ifges(7,5) * t306;
t189 = Ifges(6,1) * t307 - Ifges(6,4) * t306 + Ifges(6,5) * t381;
t301 = Ifges(5,6) * t381 + (-t501 + t507) * t383;
t302 = Ifges(5,5) * t381 + (Ifges(5,1) * t380 - t508) * t383;
t396 = -t306 * t452 + t307 * t453;
t5 = (pkin(8) * t318 + (t375 * Ifges(5,1) / 0.2e1 + m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) + (-t507 + t501 / 0.2e1) * t378 - t574) * t381 + t396 + (Ifges(4,4) + t500 - t504) * t383) * t383 + m(5) * (t272 * t287 + t273 * t288) + t439 * t307 + t440 * t306 + m(6) * (t111 * t114 + t112 * t115 + t341 * t342) + m(7) * (t105 * t96 + t106 * t99 + t151 * t152) + (t302 * t518 + t301 * t519 + pkin(8) * t319 + (-Ifges(4,4) + t504 / 0.2e1 - t500 / 0.2e1) * t381 + t452 * t305 - t453 * t304) * t381 - (-t187 / 0.2e1 - t189 / 0.2e1) * t304 - pkin(2) * t353 + t288 * t337 + t273 * t338 + t287 * t339 + t272 * t340 + t341 * t202 + t342 * t200 + t111 * t264 + t99 * t265 + t96 * t258 + t105 * t259 + t115 * t260 + t112 * t261 + t114 * t262 + t106 * t263 + t152 * t199 + t151 * t201 + (t541 + t183 / 0.2e1) * t305;
t420 = t4 * qJD(1) + t5 * qJD(2);
t393 = (t537 + t536) * t316 + (-t428 + t437) * t121;
t397 = t438 - t571;
t407 = t121 * t584 - t316 * t423;
t7 = -t456 * t162 + t397 * t120 + (t120 * t486 + t407) * t548 + t393 - t581;
t203 = Ifges(7,3) * t304 - t503;
t204 = -Ifges(6,2) * t304 - t300;
t205 = -Ifges(7,1) * t305 + t297;
t206 = -Ifges(6,1) * t305 - t506;
t8 = -t423 * t199 + t341 * t198 + (-t99 * mrSges(7,2) + t203 / 0.2e1 - t204 / 0.2e1 - t439) * t305 - (t96 * mrSges(7,2) - t206 / 0.2e1 - t205 / 0.2e1 - t440) * t304 - t563 * t383 / 0.2e1 + (-m(7) * t423 + t197) * t151 + (m(7) * t99 - t566 - t573) * t112 + (m(7) * t96 + t567 + t572) * t111;
t419 = qJD(1) * t7 + qJD(2) * t8;
t17 = -t567 * t305 - t566 * t304 + m(7) * (t304 * t99 - t305 * t96) + m(6) * (-t111 * t304 - t112 * t305) + (m(5) * (-t272 * t380 - t273 * t378) - t378 * t337 - t380 * t339) * t381;
t389 = m(5) * (-t378 * t249 + t380 * t404) * t517 + t560 * (t120 * t304 - t121 * t305);
t25 = t429 * t559 - t389;
t418 = -qJD(1) * t25 + qJD(2) * t17;
t35 = t383 * t259 - m(7) * (-t151 * t304 - t383 * t96) + t304 * t199;
t46 = (t162 / 0.4e1 + t463 / 0.4e1 + t474 / 0.4e1) * t555;
t417 = -qJD(1) * t46 - qJD(2) * t35;
t43 = (-t511 / 0.4e1 - t473 / 0.4e1 + t423 / 0.4e1) * t555 - t197 - t198;
t53 = (-t512 / 0.4e1 + t484 / 0.4e1 + t422 / 0.4e1) * t555 - t237 - t238;
t416 = qJD(2) * t43 + qJD(3) * t53;
t38 = t488 + (t112 / 0.4e1 - t450 / 0.4e1 - t448 / 0.4e1 + t464 / 0.2e1) * t555;
t411 = qJD(2) * t38 - qJD(5) * t370;
t178 = m(7) * t304;
t215 = m(7) * t402;
t410 = -qJD(2) * t178 - qJD(3) * t215;
t409 = t152 * t549 + t342 * t551;
t408 = t106 * t549 + t491 / 0.2e1;
t405 = m(7) * (pkin(5) * t190 - qJ(6) * t191);
t391 = (t534 + t533 - t583) * t316;
t10 = t454 * t191 + t456 * t190 - t405 / 0.2e1 + t391;
t241 = Ifges(7,3) * t402 + t502;
t324 = Ifges(7,5) * t402;
t242 = -Ifges(7,3) * t401 + t324;
t327 = Ifges(6,4) * t401;
t243 = -Ifges(6,2) * t402 + t327;
t244 = Ifges(6,2) * t401 + t505;
t245 = Ifges(7,1) * t401 + t324;
t246 = Ifges(7,1) * t402 - t502;
t247 = Ifges(6,1) * t401 - t505;
t248 = Ifges(6,1) * t402 + t327;
t13 = t213 * t237 + t366 * t238 - (t244 / 0.2e1 - t245 / 0.2e1 - t242 / 0.2e1 - t247 / 0.2e1) * t402 - (-t243 / 0.2e1 - t246 / 0.2e1 + t241 / 0.2e1 - t248 / 0.2e1) * t401 - t582 * t422;
t386 = t438 * t256 - (t184 / 0.4e1 - t206 / 0.4e1 - t205 / 0.4e1 - t182 / 0.4e1) * t402 - (t203 / 0.4e1 - t204 / 0.4e1 - t188 / 0.4e1 - t186 / 0.4e1) * t401 + (t305 * t532 - (-t112 / 0.2e1 + t96 / 0.2e1) * t402 - (-t99 / 0.2e1 - t111 / 0.2e1) * t401) * mrSges(7,2) + (mrSges(6,3) * t532 - t248 / 0.4e1 - t246 / 0.4e1 - t243 / 0.4e1 + t241 / 0.4e1) * t305 - (-t247 / 0.4e1 - t245 / 0.4e1 - t242 / 0.4e1 + t244 / 0.4e1) * t304 + (-t151 * t422 - t213 * t423 + t256 * t486) * t548 + t151 * t534 + t239 * t538 + t213 * t537 + t199 * t535 + t341 * t533 + t366 * t536 - t564 * t383 / 0.4e1 + (t528 * mrSges(7,2) - t545 * t304 + t548 * t584 + t437) * t257;
t392 = (-pkin(5) * t106 + qJ(6) * t105) * t549 + pkin(5) * t529 - qJ(6) * t258 / 0.2e1 - t105 * mrSges(7,3) / 0.2e1 + t106 * mrSges(7,1) / 0.2e1 + t114 * t547 + t115 * mrSges(6,2) / 0.2e1;
t2 = t386 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t381 + t392 + t396;
t400 = t10 * qJD(1) + t2 * qJD(2) + t13 * qJD(3);
t388 = -(t438 + t571) * t401 - (-t428 - t437) * t402 + t412 * t552 + (-t111 * t402 + t112 * t401 + t432) * t550 + (t401 * t96 + t402 * t99 + t432) * t548 - t472 / 0.2e1 + t467 / 0.2e1;
t12 = t454 * t307 + t456 * t306 + (pkin(8) * t553 - t492 / 0.2e1 - t493 / 0.2e1) * t383 + t388 + t409;
t390 = t395 * t552 + t560 * (t120 * t402 + t121 * t401);
t23 = t523 * t559 - t390;
t32 = (t401 ^ 2 + t402 ^ 2) * t585 + (m(5) * qJ(4) + mrSges(5,3)) * t457 + t579 * (t256 * t402 + t257 * t401);
t399 = -qJD(1) * t23 + qJD(2) * t12 + qJD(3) * t32;
t394 = (-t151 * t402 - t213 * t304 - t257 * t383) * t548 + t239 * t528 + t199 * t520;
t30 = (t441 + t578) * mrSges(7,2) + t394 + t408;
t60 = t582 * t402;
t88 = (t444 + t540) * m(7);
t398 = -qJD(1) * t88 - qJD(2) * t30 + qJD(3) * t60;
t377 = t383 ^ 2;
t376 = t381 ^ 2;
t351 = t376 * pkin(8) * t468;
t122 = m(7) * t257 + t494;
t110 = m(7) * t535 + t583;
t87 = m(7) * t444 - t190 * t548;
t71 = m(7) * t538 + t423 * t548;
t47 = (-t463 - t474 + t162) * t548;
t39 = t586 * t96 + t259;
t29 = mrSges(7,2) * t441 + t496 / 0.2e1 + t394 - t408;
t26 = (t552 + t560) * t451 + t389;
t24 = 0.2e1 * (m(7) / 0.4e1 + m(6) / 0.4e1 + m(5) / 0.4e1) * t317 + t390;
t11 = -t495 / 0.2e1 + t498 / 0.2e1 + t497 / 0.2e1 + t499 / 0.2e1 + t373 * t552 + mrSges(5,2) * t465 / 0.2e1 + mrSges(5,1) * t470 / 0.2e1 + t388 - t409;
t9 = t405 / 0.2e1 - t191 * t546 + mrSges(7,3) * t539 + t391 + t576 * t540;
t6 = t407 * t548 + (t486 * t548 + t397) * t120 + t393 - t576 * t162 / 0.2e1 + t581;
t3 = -t387 + t385;
t1 = Ifges(6,6) * t527 + Ifges(7,6) * t526 + t517 * t574 + t525 * t575 + t386 - t392;
t16 = [qJD(2) * t14 + qJD(3) * t15, t3 * qJD(3) + t26 * qJD(4) + t6 * qJD(5) + t47 * qJD(6) + t483 + (-t162 * t262 + t162 * t263 + t163 * t259 + t163 * t260 + t270 * t339 + t271 * t337 + ((-mrSges(4,1) * t383 - mrSges(3,1) + t490) * t382 + (-mrSges(3,2) + (t376 + t377) * mrSges(4,3) + t558 * t381) * t384) * t379 + (t151 * t451 + t162 * t99 + t163 * t96) * t586 + 0.2e1 * (-t111 * t162 + t112 * t163 + t341 * t451) * t550 + 0.2e1 * (t270 * t272 + t271 * t273 + t351) * t552 + m(4) * (t351 + (pkin(8) * t377 * t384 - pkin(2) * t382) * t379)) * qJD(2), t482 + t3 * qJD(2) + t24 * qJD(4) + t9 * qJD(5) + t87 * qJD(6) + ((t213 * t317 + t434) * t548 + (t317 * t366 + t434) * t550 + (-qJ(4) * t316 * t457 - pkin(3) * t317) * t552) * t554 + ((-mrSges(5,3) * t457 + mrSges(4,2)) * t316 + (t487 + t569) * t317 + t585 * (-t190 * t402 - t191 * t401)) * qJD(3), qJD(2) * t26 + qJD(3) * t24, t6 * qJD(2) + t9 * qJD(3) + (-t576 * t121 + (mrSges(6,2) - mrSges(7,3)) * t120) * qJD(5) + ((-pkin(5) * t121 - qJ(6) * t120) * qJD(5) / 0.2e1 + qJD(6) * t543) * t555, m(7) * qJD(5) * t121 + qJD(2) * t47 + qJD(3) * t87; qJD(3) * t4 - qJD(4) * t25 + qJD(5) * t7 - qJD(6) * t46 - t483, qJD(3) * t5 + qJD(4) * t17 + qJD(5) * t8 - qJD(6) * t35, t11 * qJD(4) + t1 * qJD(5) + t29 * qJD(6) + ((-pkin(3) * t373 + qJ(4) * t562) * t552 + (-t114 * t256 + t115 * t257 + t342 * t366) * t550 + (t105 * t257 + t106 * t256 + t152 * t213) * t548) * t554 + t420 + (t105 * t494 + t242 * t526 + t244 * t527 - t401 * t541 + pkin(8) * t490 + t106 * t402 * mrSges(7,2) - Ifges(4,6) * t381 + t378 * t302 / 0.2e1 + t366 * t202 + t342 * t240 - pkin(3) * t319 + t152 * t239 + t213 * t201 + ((Ifges(5,2) * t380 + t508) * t519 + (Ifges(5,1) * t378 + t507) * t518 + Ifges(4,5) + t487 * pkin(8)) * t383 + t301 * t518 + t183 * t522 + (t248 + t246) * t525 + (t189 + t187) * t402 / 0.2e1 + t568 * t257 + (-t264 + t265) * t256 + (t338 * t380 - t340 * t378) * qJ(4) + (-t114 * t402 + t115 * t401) * mrSges(6,3) + t562 * mrSges(5,3) + (Ifges(5,5) * t378 + Ifges(5,6) * t380 - t401 * t587 + t575 * t402) * t517) * qJD(3), qJD(3) * t11 + qJD(5) * t71 + t418, t1 * qJD(3) + t71 * qJD(4) + (t424 * mrSges(7,2) + t111 * t556 + t112 * t557 + t563) * qJD(5) + t39 * qJD(6) + t419, qJD(3) * t29 + qJD(5) * t39 + t417; -qJD(2) * t4 - qJD(4) * t23 + qJD(5) * t10 + qJD(6) * t88 - t482, qJD(4) * t12 + qJD(5) * t2 + qJD(6) * t30 - t420, qJD(4) * t32 + qJD(5) * t13 - qJD(6) * t60, qJD(5) * t110 + t399, t110 * qJD(4) + (t421 * mrSges(7,2) - t256 * t556 + t257 * t557 + t564) * qJD(5) + t122 * qJD(6) + t400, qJD(5) * t122 - t398; qJD(2) * t25 + qJD(3) * t23, -qJD(3) * t12 - qJD(5) * t43 - qJD(6) * t178 - t418, -qJD(5) * t53 - qJD(6) * t215 - t399, 0, -t416, t410; -qJD(2) * t7 - qJD(3) * t10, -qJD(3) * t2 + qJD(4) * t43 - qJD(6) * t38 - t419, qJD(4) * t53 - t400, t416, t370 * qJD(6), -t411; qJD(2) * t46 - qJD(3) * t88, -qJD(3) * t30 + qJD(4) * t178 + qJD(5) * t38 - t417, qJD(4) * t215 + t398, -t410, t411, 0;];
Cq  = t16;
