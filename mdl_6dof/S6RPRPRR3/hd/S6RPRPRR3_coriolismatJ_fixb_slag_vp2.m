% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:45
% EndTime: 2019-03-09 03:41:05
% DurationCPUTime: 11.11s
% Computational Cost: add. (25827->630), mult. (53557->869), div. (0->0), fcn. (57602->10), ass. (0->312)
t380 = sin(pkin(11));
t381 = cos(pkin(11));
t553 = sin(qJ(5));
t554 = cos(qJ(5));
t345 = -t380 * t554 - t381 * t553;
t384 = sin(qJ(3));
t326 = t345 * t384;
t344 = -t553 * t380 + t554 * t381;
t328 = t344 * t384;
t383 = sin(qJ(6));
t385 = cos(qJ(6));
t247 = t326 * t383 + t328 * t385;
t438 = t385 * t326 - t328 * t383;
t635 = t247 * mrSges(7,1) + t438 * mrSges(7,2);
t604 = qJD(6) * t635;
t287 = t344 * t383 - t345 * t385;
t437 = t385 * t344 + t345 * t383;
t463 = Ifges(7,5) * t437 - Ifges(7,6) * t287;
t541 = pkin(8) + qJ(4);
t357 = t541 * t380;
t359 = t541 * t381;
t296 = -t357 * t553 + t359 * t554;
t252 = t344 * pkin(9) + t296;
t295 = -t554 * t357 - t359 * t553;
t402 = t345 * pkin(9) + t295;
t149 = t252 * t385 + t383 * t402;
t600 = -t252 * t383 + t385 * t402;
t645 = -t149 * mrSges(7,1) - t600 * mrSges(7,2);
t30 = t463 + t645;
t651 = t30 * qJD(6);
t465 = Ifges(7,5) * t438 - Ifges(7,6) * t247;
t386 = cos(qJ(3));
t327 = t386 * t345;
t370 = sin(pkin(10)) * pkin(1) + pkin(7);
t362 = t386 * t370;
t474 = t380 * t386;
t336 = pkin(4) * t474 + t362;
t265 = -pkin(5) * t327 + t336;
t329 = t386 * t344;
t248 = t327 * t383 + t329 * t385;
t511 = t248 * mrSges(7,2);
t245 = t327 * t385 - t329 * t383;
t515 = t245 * mrSges(7,1);
t467 = t515 / 0.2e1 - t511 / 0.2e1;
t502 = t329 * mrSges(6,2);
t504 = t327 * mrSges(6,1);
t436 = t467 - t502 / 0.2e1 + t504 / 0.2e1;
t594 = -m(7) / 0.2e1;
t596 = -m(6) / 0.2e1;
t650 = t265 * t594 + t336 * t596 + t436;
t592 = m(7) * pkin(5);
t649 = -(t245 * t385 + t248 * t383) * t592 / 0.2e1 - t436;
t530 = Ifges(7,4) * t247;
t132 = Ifges(7,2) * t438 - t386 * Ifges(7,6) + t530;
t196 = t386 * mrSges(7,2) + mrSges(7,3) * t438;
t371 = -pkin(4) * t381 - pkin(3);
t313 = -pkin(5) * t344 + t371;
t512 = t247 * mrSges(7,3);
t198 = -mrSges(7,1) * t386 - t512;
t584 = -t198 / 0.2e1;
t639 = Ifges(7,1) * t438 - t530;
t647 = t600 / 0.2e1;
t648 = (t639 / 0.4e1 - t132 / 0.4e1) * t287 + t149 * t584 + t196 * t647 + t313 * t635 / 0.2e1;
t361 = t384 * t370;
t475 = t380 * t384;
t335 = pkin(4) * t475 + t361;
t264 = -pkin(5) * t326 + t335;
t646 = t264 * t635;
t529 = Ifges(7,4) * t287;
t178 = Ifges(7,2) * t437 + t529;
t640 = Ifges(7,1) * t437 - t529;
t644 = t178 / 0.4e1 - t640 / 0.4e1;
t557 = t384 / 0.2e1;
t628 = t287 * mrSges(7,1);
t452 = t628 / 0.2e1;
t278 = Ifges(7,4) * t437;
t180 = Ifges(7,1) * t287 + t278;
t623 = -Ifges(7,2) * t287 + t278;
t638 = t623 + t180;
t229 = Ifges(7,4) * t438;
t134 = Ifges(7,1) * t247 - t386 * Ifges(7,5) + t229;
t624 = -Ifges(7,2) * t247 + t229;
t637 = t624 + t134;
t636 = t437 * mrSges(7,2) + t628;
t549 = t384 * pkin(3);
t360 = -qJ(4) * t386 + t549;
t297 = t381 * t360 + t380 * t361;
t298 = t380 * t360 - t361 * t381;
t414 = -t297 * t380 + t298 * t381;
t410 = -Ifges(7,5) * t248 / 0.2e1 - Ifges(7,6) * t245 / 0.2e1;
t471 = t381 * t386;
t263 = t384 * pkin(4) - pkin(8) * t471 + t297;
t288 = -pkin(8) * t474 + t298;
t167 = t554 * t263 - t288 * t553;
t117 = t384 * pkin(5) - t329 * pkin(9) + t167;
t168 = t553 * t263 + t554 * t288;
t118 = pkin(9) * t327 + t168;
t64 = t117 * t385 - t118 * t383;
t65 = t117 * t383 + t118 * t385;
t434 = Ifges(7,3) * t557 - t65 * mrSges(7,2) / 0.2e1 + t64 * mrSges(7,1) / 0.2e1 - t410;
t532 = Ifges(6,4) * t328;
t239 = Ifges(6,2) * t326 - t386 * Ifges(6,6) + t532;
t317 = Ifges(6,4) * t326;
t241 = Ifges(6,1) * t328 - t386 * Ifges(6,5) + t317;
t505 = t326 * mrSges(6,3);
t299 = mrSges(6,2) * t386 + t505;
t503 = t328 * mrSges(6,3);
t301 = -mrSges(6,1) * t386 - t503;
t443 = t328 * mrSges(6,1) + t326 * mrSges(6,2);
t562 = -t345 / 0.4e1;
t565 = -t344 / 0.4e1;
t572 = t296 / 0.2e1;
t291 = -t345 * mrSges(6,1) + t344 * mrSges(6,2);
t573 = -t291 / 0.2e1;
t586 = -t636 / 0.2e1;
t608 = Ifges(6,5) * t344 + Ifges(6,6) * t345 + t463;
t341 = Ifges(6,4) * t344;
t294 = -Ifges(6,1) * t345 + t341;
t609 = Ifges(6,2) * t345 + t294 + t341;
t622 = -(t299 / 0.2e1 - t505 / 0.2e1) * t295 + t264 * t586 + t301 * t572 + t335 * t573 + t241 * t565 + t239 * t562 - t371 * t443 / 0.2e1 - t609 * t326 / 0.4e1 + t608 * t386 / 0.4e1 - t637 * t437 / 0.4e1 - t638 * t438 / 0.4e1 + t644 * t247 - t648;
t621 = -t132 / 0.2e1;
t574 = t287 / 0.2e1;
t379 = t381 ^ 2;
t620 = t379 / 0.2e1;
t618 = t437 / 0.2e1;
t616 = t438 / 0.2e1;
t468 = t385 * t438;
t607 = t635 + t443;
t537 = mrSges(5,2) * t381;
t539 = mrSges(5,1) * t380;
t412 = t539 / 0.2e1 + t537 / 0.2e1;
t606 = Ifges(6,5) * t326 - Ifges(6,6) * t328 + t465;
t598 = 2 * qJD(3);
t597 = m(5) / 0.2e1;
t595 = m(6) / 0.2e1;
t593 = m(7) / 0.2e1;
t591 = -mrSges(7,3) / 0.2e1;
t588 = -t600 / 0.2e1;
t587 = -t149 / 0.2e1;
t583 = -t438 / 0.2e1;
t581 = t245 / 0.2e1;
t579 = t247 / 0.2e1;
t578 = t248 / 0.2e1;
t577 = -t287 / 0.2e1;
t570 = t326 / 0.2e1;
t569 = t327 / 0.2e1;
t568 = -t328 / 0.2e1;
t567 = t328 / 0.2e1;
t566 = t329 / 0.2e1;
t564 = t344 / 0.2e1;
t563 = -t345 / 0.2e1;
t561 = t345 / 0.2e1;
t560 = -t380 / 0.2e1;
t559 = t381 / 0.2e1;
t558 = t383 / 0.2e1;
t556 = -t386 / 0.2e1;
t552 = m(5) * t370;
t551 = pkin(5) * t328;
t550 = pkin(5) * t345;
t548 = t386 * pkin(5);
t451 = -cos(pkin(10)) * pkin(1) - pkin(2);
t342 = -pkin(3) * t386 - t384 * qJ(4) + t451;
t331 = t381 * t342;
t472 = t381 * t384;
t249 = -pkin(8) * t472 + t331 + (-t370 * t380 - pkin(4)) * t386;
t290 = t380 * t342 + t381 * t362;
t262 = -pkin(8) * t475 + t290;
t153 = t554 * t249 - t262 * t553;
t109 = -t328 * pkin(9) + t153;
t106 = t109 - t548;
t154 = t249 * t553 + t262 * t554;
t110 = t326 * pkin(9) + t154;
t493 = t110 * t383;
t56 = t106 * t385 - t493;
t547 = t56 * mrSges(7,2);
t492 = t110 * t385;
t57 = t106 * t383 + t492;
t546 = t57 * mrSges(7,1);
t60 = -t109 * t383 - t492;
t545 = t60 * mrSges(7,1);
t61 = t109 * t385 - t493;
t544 = t61 * mrSges(7,2);
t538 = mrSges(6,1) * t326;
t536 = mrSges(6,2) * t328;
t535 = mrSges(6,3) * t345;
t534 = Ifges(5,4) * t380;
t533 = Ifges(5,4) * t381;
t531 = Ifges(6,4) * t345;
t528 = Ifges(5,5) * t381;
t525 = Ifges(5,6) * t380;
t521 = pkin(5) * qJD(5);
t509 = t437 * mrSges(7,3);
t506 = t287 * mrSges(7,3);
t501 = t344 * mrSges(6,3);
t500 = t380 * Ifges(5,2);
t499 = t56 * t437;
t498 = t57 * t287;
t497 = t60 * t287;
t496 = t61 * t437;
t495 = -mrSges(5,1) * t381 + mrSges(5,2) * t380 - mrSges(4,1);
t491 = t149 * t247;
t16 = t247 * t584 + t635 * t556 + t196 * t616 + (-t438 ^ 2 / 0.2e1 - t247 ^ 2 / 0.2e1) * mrSges(7,3);
t490 = t16 * qJD(1);
t489 = t247 * t198;
t488 = t247 * t385;
t487 = t438 * t196;
t486 = t438 * t383;
t483 = t287 * t385;
t482 = t437 * t383;
t479 = t326 * t299;
t478 = t328 * t301;
t477 = t345 * t326;
t352 = t384 * mrSges(5,1) - mrSges(5,3) * t471;
t476 = t380 * t352;
t350 = -t384 * mrSges(5,2) - mrSges(5,3) * t474;
t473 = t381 * t350;
t470 = t384 * t386;
t469 = t385 * t196;
t378 = t380 ^ 2;
t461 = t378 + t379;
t460 = qJD(3) * t386;
t459 = m(5) * t557;
t458 = t551 / 0.2e1;
t457 = -t550 / 0.2e1;
t453 = mrSges(7,3) * t588;
t448 = t636 * t556;
t444 = t461 * mrSges(5,3);
t441 = -Ifges(6,2) * t328 + t317;
t439 = t461 * qJ(4);
t433 = Ifges(6,1) * t326 - t532;
t432 = Ifges(6,1) * t344 + t531;
t133 = Ifges(7,4) * t248 + Ifges(7,2) * t245 + Ifges(7,6) * t384;
t135 = Ifges(7,1) * t248 + Ifges(7,4) * t245 + Ifges(7,5) * t384;
t142 = -mrSges(7,1) * t438 + mrSges(7,2) * t247;
t143 = t511 - t515;
t197 = -mrSges(7,2) * t384 + mrSges(7,3) * t245;
t199 = mrSges(7,1) * t384 - mrSges(7,3) * t248;
t240 = Ifges(6,4) * t329 + Ifges(6,2) * t327 + Ifges(6,6) * t384;
t242 = Ifges(6,1) * t329 + Ifges(6,4) * t327 + Ifges(6,5) * t384;
t254 = t502 - t504;
t289 = -t362 * t380 + t331;
t300 = -t384 * mrSges(6,2) + mrSges(6,3) * t327;
t302 = t384 * mrSges(6,1) - mrSges(6,3) * t329;
t324 = Ifges(5,6) * t384 + (-t500 + t533) * t386;
t325 = Ifges(5,5) * t384 + (t381 * Ifges(5,1) - t534) * t386;
t339 = (t537 + t539) * t386;
t349 = t386 * mrSges(5,2) - mrSges(5,3) * t475;
t351 = -t386 * mrSges(5,1) - mrSges(5,3) * t472;
t411 = -Ifges(6,5) * t329 / 0.2e1 - Ifges(6,6) * t327 / 0.2e1;
t4 = m(5) * (t289 * t297 + t290 * t298) + m(6) * (t153 * t167 + t154 * t168 + t335 * t336) + m(7) * (t264 * t265 + t56 * t64 + t57 * t65) + t134 * t578 + t135 * t579 + t132 * t581 + t241 * t566 + t242 * t567 + t239 * t569 + t240 * t570 + t336 * (t536 - t538) + t298 * t349 + t290 * t350 + t297 * t351 + t289 * t352 + t335 * t254 + t168 * t299 + t154 * t300 + t167 * t301 + t153 * t302 + t264 * t143 + t265 * t142 + t65 * t196 + t57 * t197 + t64 * t198 + t56 * t199 + (t451 * mrSges(4,2) + (Ifges(4,1) - Ifges(4,2) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3) + Ifges(5,1) * t620 + (t537 + t552) * t370 + (t370 * mrSges(5,1) - t533 + t500 / 0.2e1) * t380) * t384 + t410 + t411 + (Ifges(4,4) + t525 - t528) * t386) * t386 + t133 * t616 + (t451 * mrSges(4,1) + Ifges(7,5) * t579 + Ifges(7,6) * t616 + Ifges(6,5) * t567 + Ifges(6,6) * t570 + t370 * t339 + t324 * t560 + t325 * t559 + (-Ifges(4,4) + t528 / 0.2e1 - t525 / 0.2e1) * t384) * t384;
t406 = t349 * t559 + t351 * t560;
t415 = -t289 * t380 + t290 * t381;
t9 = t198 * t581 + t199 * t616 + t196 * t578 + t197 * t579 + t301 * t569 + t302 * t570 + t299 * t566 + t300 * t567 + (t473 / 0.2e1 - t476 / 0.2e1 + t142 / 0.2e1 - t538 / 0.2e1 + t536 / 0.2e1 + t412 * t384) * t384 + (-t143 / 0.2e1 - t254 / 0.2e1 - t339 / 0.2e1 + t406) * t386 + (t56 * t245 + t65 * t247 + t57 * t248 + t264 * t384 - t265 * t386 + t438 * t64) * t593 + (t153 * t327 + t154 * t329 + t167 * t326 + t168 * t328 + t335 * t384 - t336 * t386) * t595 + ((t415 - t362) * t386 + (t414 + t361) * t384) * t597;
t425 = t4 * qJD(1) + t9 * qJD(2);
t424 = -t247 * t56 + t438 * t57;
t10 = -t489 / 0.2e1 + (t61 * t247 - t328 * t548 + t438 * t60 + t424) * t593 + t487 / 0.2e1 - t478 / 0.2e1 + t479 / 0.2e1 + (-t247 * t579 + t438 * t583) * mrSges(7,3) + (-t326 ^ 2 / 0.2e1 - t328 ^ 2 / 0.2e1) * mrSges(6,3) + t607 * t556;
t7 = t247 * t621 + t60 * t198 + m(7) * (t264 * t551 + t56 * t60 + t57 * t61) + t61 * t196 + t142 * t551 + t646 + t639 * t579 + t335 * t443 + t433 * t567 + t239 * t568 - t154 * t301 + t153 * t299 + (-t57 * t247 - t438 * t56) * mrSges(7,3) + (-t153 * t326 - t154 * t328) * mrSges(6,3) + (t241 + t441) * t570 + t606 * t556 + t637 * t616;
t423 = t7 * qJD(1) + t10 * qJD(2);
t8 = t465 * t556 - t57 * t198 + t56 * t196 + t646 + (t621 - t57 * mrSges(7,3) + t639 / 0.2e1) * t247 + (t624 / 0.2e1 - t56 * mrSges(7,3) + t134 / 0.2e1) * t438;
t422 = t8 * qJD(1) + t16 * qJD(2);
t39 = m(7) * (t245 * t438 + t247 * t248 - t470) + m(6) * (t326 * t327 + t328 * t329 - t470) + m(5) * (-0.1e1 + t461) * t470;
t421 = t9 * qJD(1) + t39 * qJD(2);
t420 = t10 * qJD(1);
t21 = t487 - t489 + t479 - t478 + m(7) * t424 + m(6) * (-t153 * t328 + t154 * t326) + (m(5) * (-t289 * t381 - t290 * t380) - t380 * t349 - t381 * t351) * t384;
t419 = qJD(1) * t21;
t53 = (t567 + t488 / 0.2e1 - t486 / 0.2e1) * t592 + t607;
t66 = (t563 + t483 / 0.2e1 - t482 / 0.2e1) * t592 + t636 + t291;
t418 = qJD(1) * t53 + qJD(3) * t66;
t75 = 0.2e1 * t618 * mrSges(7,2) + 0.2e1 * t452;
t417 = -qJD(1) * t635 - qJD(3) * t75;
t416 = -t264 * t345 + t313 * t328;
t407 = t247 * t574 + t438 * t618;
t401 = (t486 - t488) * t592;
t176 = -mrSges(7,1) * t437 + mrSges(7,2) * t287;
t293 = Ifges(6,2) * t344 - t531;
t391 = Ifges(6,3) * t557 + t167 * mrSges(6,1) / 0.2e1 - t168 * mrSges(6,2) / 0.2e1 - t411 + t434;
t393 = t385 * t199 / 0.2e1 + (t383 * t65 + t385 * t64) * t593 + t197 * t558;
t395 = (t57 + t60) * t600 + (-t56 + t61) * t149;
t1 = (t498 / 0.2e1 - t496 / 0.2e1 + t499 / 0.2e1 + t497 / 0.2e1 + t600 * t616 + t491 / 0.2e1) * mrSges(7,3) + t391 + t317 * t565 + (t142 * t561 + t176 * t568 + t416 * t594 + t393) * pkin(5) + t395 * t594 + (-t531 / 0.2e1 + t293 / 0.4e1 + mrSges(6,3) * t572 + (Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1) * t344) * t328 + Ifges(6,1) * t477 / 0.4e1 + t622;
t11 = -t176 * t550 + t178 * t577 + t371 * t291 + t293 * t561 + t640 * t574 + t432 * t563 + t609 * t564 + t638 * t618 + (-m(7) * t550 + t636) * t313;
t390 = (t586 + t573) * t386 + (t247 * t577 + t437 * t583 + t407) * mrSges(7,3) + t345 * t548 * t593;
t18 = t390 + t649;
t399 = -t1 * qJD(1) + t18 * qJD(2) + t11 * qJD(3);
t22 = t313 * t636 + (t640 / 0.2e1 - t178 / 0.2e1) * t287 + (t180 / 0.2e1 + t623 / 0.2e1) * t437;
t27 = t448 - t467;
t389 = (t134 / 0.4e1 + t624 / 0.4e1) * t437 + (t623 / 0.4e1 + t453 + t180 / 0.4e1) * t438 + (mrSges(7,3) * t587 - t644) * t247 + t264 * t636 / 0.2e1 - t386 * t463 / 0.4e1 + t648;
t6 = t389 - t434;
t398 = t6 * qJD(1) + t27 * qJD(2) + t22 * qJD(3);
t388 = t407 * mrSges(7,3) + (t326 * t564 + t328 * t563) * mrSges(6,3) + t415 * t597 + (t153 * t345 + t154 * t344 - t295 * t328 + t296 * t326) * t595 + (t149 * t438 - t247 * t600 - t287 * t56 + t437 * t57) * t593 + t198 * t577 + t196 * t618 + t299 * t564 + t301 * t561 + t406;
t13 = t388 + (-t552 / 0.2e1 - t412) * t386 + t650;
t33 = (t287 ^ 2 + t437 ^ 2) * mrSges(7,3) + (t344 ^ 2 + t345 ^ 2) * mrSges(6,3) + t444 + m(7) * (t149 * t437 - t287 * t600) + m(6) * (t295 * t345 + t296 * t344) + m(5) * t439;
t392 = (t247 * t437 - t287 * t438) * t593 + (t328 * t344 + t477) * t595;
t48 = (t594 + t596 + (t378 / 0.2e1 + t620 - 0.1e1 / 0.2e1) * m(5)) * t384 + t392;
t397 = qJD(1) * t13 + qJD(2) * t48 + qJD(3) * t33;
t14 = (-t61 / 0.2e1 + t56 / 0.2e1) * mrSges(7,2) + (t60 / 0.2e1 + t57 / 0.2e1) * mrSges(7,1) + (t198 * t558 - t469 / 0.2e1 + (t468 / 0.2e1 + t247 * t558) * mrSges(7,3)) * pkin(5);
t29 = (t588 + t647) * mrSges(7,2) + (t587 + t149 / 0.2e1) * mrSges(7,1);
t356 = (mrSges(7,1) * t383 + mrSges(7,2) * t385) * pkin(5);
t394 = qJD(1) * t14 - qJD(3) * t29 + qJD(5) * t356;
t346 = t356 * qJD(6);
t292 = -mrSges(6,1) * t344 - mrSges(6,2) * t345;
t119 = m(7) * t457 + (t482 - t483) * t592 / 0.2e1;
t95 = m(7) * t458 + t401 / 0.2e1;
t76 = t452 - t628 / 0.2e1;
t47 = t459 * t461 + t392 + t459 + (m(6) + m(7)) * t557;
t26 = t448 + t467;
t17 = t390 - t649;
t15 = -t546 / 0.2e1 - t547 / 0.2e1 + t545 / 0.2e1 - t544 / 0.2e1 + (-(t198 + t512) * t383 / 0.2e1 + t468 * t591 + t469 / 0.2e1) * pkin(5) + t465;
t12 = t362 * t597 + t412 * t386 + t388 - t650;
t5 = t389 + t434;
t3 = qJD(3) * t9 + qJD(5) * t10 + qJD(6) * t16;
t2 = t391 + t393 * pkin(5) + (pkin(5) * t416 + t395) * t593 + t433 * t562 + t344 * t441 / 0.4e1 + t142 * t457 - t296 * t503 / 0.2e1 + t176 * t458 + t438 * t453 + (t499 + t497) * t591 + (t432 / 0.4e1 - t293 / 0.4e1) * t328 + (-t498 + t496 - t491) * mrSges(7,3) / 0.2e1 - t622;
t19 = [qJD(3) * t4 + qJD(4) * t21 + qJD(5) * t7 + qJD(6) * t8, t3, t12 * qJD(4) + t2 * qJD(5) + t5 * qJD(6) + (Ifges(4,5) + (Ifges(5,2) * t381 + t534) * t560 + (Ifges(5,1) * t380 + t533) * t559 + t495 * t370) * t460 + ((-pkin(3) * t362 + qJ(4) * t414) * t597 + (t167 * t295 + t168 * t296 + t336 * t371) * t595 + (t149 * t65 + t265 * t313 + t600 * t64) * t593) * t598 + t425 + (t167 * t535 + t168 * t501 + t135 * t574 + t180 * t578 + t178 * t581 + t324 * t559 + t242 * t563 + t240 * t564 + t294 * t566 + t293 * t569 - t64 * t506 - Ifges(4,6) * t384 + t380 * t325 / 0.2e1 + t371 * t254 + t336 * t292 - pkin(3) * t339 + t313 * t143 + t296 * t300 + t295 * t302 + t265 * t176 + t149 * t197 + mrSges(4,2) * t361 + t65 * t509 + t414 * mrSges(5,3) + (-t476 + t473) * qJ(4) + (Ifges(5,5) * t380 - Ifges(6,5) * t345 + Ifges(7,5) * t287 + Ifges(5,6) * t381 + Ifges(6,6) * t344 + Ifges(7,6) * t437) * t557 + t600 * t199 + t133 * t618) * qJD(3), qJD(3) * t12 + qJD(5) * t95 + t419, t2 * qJD(3) + t95 * qJD(4) + (-t154 * mrSges(6,1) - t153 * mrSges(6,2) - t544 + t545 + t606) * qJD(5) + t15 * qJD(6) + (m(7) * (t383 * t61 + t385 * t60) + (-t247 * t383 - t468) * mrSges(7,3)) * t521 + t423, t5 * qJD(3) + t15 * qJD(5) + (t465 - t546 - t547) * qJD(6) + t422; t3, t39 * qJD(3), t47 * qJD(4) + t17 * qJD(5) + t26 * qJD(6) + (-mrSges(4,2) + t444) * t460 + ((t149 * t248 + t245 * t600 + t313 * t384) * t593 + (t295 * t327 + t296 * t329 + t371 * t384) * t595 + (t386 * t439 - t549) * t597) * t598 + t421 + (-t245 * t506 + t248 * t509 + t327 * t535 + t329 * t501 + (t176 + t292 + t495) * t384) * qJD(3), qJD(3) * t47, t17 * qJD(3) + (t401 - t607) * qJD(5) - t604 + t420, t26 * qJD(3) - qJD(5) * t635 + t490 - t604; qJD(4) * t13 - qJD(5) * t1 + qJD(6) * t6 - t425, qJD(4) * t48 + qJD(5) * t18 + qJD(6) * t27 - t421, qJD(4) * t33 + qJD(5) * t11 + qJD(6) * t22, qJD(5) * t119 + qJD(6) * t76 + t397, t119 * qJD(4) + (-t296 * mrSges(6,1) - t295 * mrSges(6,2) + t608 + t645) * qJD(5) + t651 + (m(7) * (-t149 * t385 + t383 * t600) + (-t287 * t383 - t385 * t437) * mrSges(7,3)) * t521 + t399, t76 * qJD(4) + t30 * qJD(5) + t398 + t651; -qJD(3) * t13 + qJD(5) * t53 - t419 + t604, -qJD(3) * t48, qJD(5) * t66 + qJD(6) * t75 - t397, 0, t418, -t417; qJD(3) * t1 - qJD(4) * t53 - qJD(6) * t14 - t423, -qJD(3) * t18 - t420, -qJD(4) * t66 + qJD(6) * t29 - t399, -t418, -t346, -t346 - t394; -qJD(3) * t6 - qJD(4) * t635 + qJD(5) * t14 - t422, -t27 * qJD(3) - t490, -qJD(4) * t75 - qJD(5) * t29 - t398, t417, t394, 0;];
Cq  = t19;
