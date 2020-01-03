% Calculate vector of inverse dynamics joint torques for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:00
% EndTime: 2019-12-31 22:33:01
% DurationCPUTime: 31.20s
% Computational Cost: add. (16198->821), mult. (39248->1140), div. (0->0), fcn. (31035->14), ass. (0->385)
t331 = cos(qJ(2));
t326 = sin(qJ(2));
t321 = sin(pkin(5));
t434 = qJD(1) * t321;
t398 = t326 * t434;
t322 = cos(pkin(5));
t433 = qJD(1) * t322;
t419 = pkin(1) * t433;
t256 = -pkin(7) * t398 + t331 * t419;
t350 = (pkin(2) * t326 - pkin(8) * t331) * t321;
t257 = qJD(1) * t350;
t325 = sin(qJ(3));
t330 = cos(qJ(3));
t167 = -t256 * t325 + t330 * t257;
t332 = -pkin(9) - pkin(8);
t402 = qJD(3) * t332;
t441 = t330 * t331;
t601 = -(pkin(3) * t326 - pkin(9) * t441) * t434 - t167 + t330 * t402;
t168 = t330 * t256 + t325 * t257;
t397 = t331 * t434;
t380 = t325 * t397;
t600 = -pkin(9) * t380 - t325 * t402 + t168;
t287 = qJD(3) - t397;
t279 = qJD(4) + t287;
t503 = t279 / 0.2e1;
t306 = qJD(2) + t433;
t227 = t306 * t330 - t325 * t398;
t228 = t306 * t325 + t330 * t398;
t324 = sin(qJ(4));
t329 = cos(qJ(4));
t357 = t227 * t324 + t329 * t228;
t509 = t357 / 0.2e1;
t599 = Ifges(5,4) * t509 + Ifges(5,6) * t503;
t387 = t329 * t227 - t228 * t324;
t151 = qJD(5) - t387;
t514 = -t151 / 0.2e1;
t323 = sin(qJ(5));
t328 = cos(qJ(5));
t123 = t279 * t323 + t328 * t357;
t520 = -t123 / 0.2e1;
t122 = t279 * t328 - t323 * t357;
t522 = -t122 / 0.2e1;
t598 = Ifges(6,5) * t520 + Ifges(6,6) * t522 + Ifges(6,3) * t514;
t594 = mrSges(5,2) - mrSges(6,3);
t511 = t387 / 0.2e1;
t597 = Ifges(5,2) * t511;
t372 = -mrSges(6,1) * t328 + mrSges(6,2) * t323;
t593 = mrSges(5,1) - t372;
t543 = m(6) * pkin(4);
t596 = -t543 - t593;
t542 = m(6) * pkin(10);
t385 = -t542 + t594;
t292 = t332 * t325;
t293 = t332 * t330;
t355 = t329 * t292 + t293 * t324;
t571 = qJD(4) * t355 + t601 * t324 - t600 * t329;
t303 = pkin(7) * t397;
t259 = t326 * t419 + t303;
t429 = qJD(3) * t325;
t559 = -t259 + (-t380 + t429) * pkin(3);
t207 = -pkin(2) * t306 - t256;
t157 = -pkin(3) * t227 + t207;
t208 = pkin(8) * t306 + t259;
t245 = (-pkin(2) * t331 - pkin(8) * t326 - pkin(1)) * t321;
t217 = qJD(1) * t245;
t140 = -t208 * t325 + t330 * t217;
t113 = -pkin(9) * t228 + t140;
t106 = pkin(3) * t287 + t113;
t141 = t208 * t330 + t217 * t325;
t114 = pkin(9) * t227 + t141;
t442 = t329 * t114;
t63 = t106 * t324 + t442;
t60 = pkin(10) * t279 + t63;
t76 = -pkin(4) * t387 - pkin(10) * t357 + t157;
t27 = -t323 * t60 + t328 * t76;
t28 = t323 * t76 + t328 * t60;
t595 = t598 + t599;
t590 = -t157 * mrSges(5,1) - t27 * mrSges(6,1) + t28 * mrSges(6,2) + t63 * mrSges(5,3) + t595 + t597;
t592 = -pkin(10) * t398 + t571;
t274 = t324 * t325 - t329 * t330;
t556 = qJD(3) + qJD(4);
t196 = t556 * t274;
t275 = t324 * t330 + t325 * t329;
t197 = t556 * t275;
t209 = t275 * t397;
t210 = t274 * t397;
t591 = t559 + (t196 - t210) * pkin(10) + (t197 - t209) * pkin(4);
t430 = qJD(2) * t331;
t263 = (qJD(1) * t430 + qJDD(1) * t326) * t321;
t422 = qJDD(1) * t322;
t305 = qJDD(2) + t422;
t144 = -qJD(3) * t228 - t263 * t325 + t305 * t330;
t450 = t321 * t326;
t307 = pkin(7) * t450;
t382 = qJD(2) * t419;
t413 = pkin(1) * t422;
t174 = -qJD(2) * t303 - qJDD(1) * t307 - t326 * t382 + t331 * t413;
t162 = -pkin(2) * t305 - t174;
t102 = -pkin(3) * t144 + t162;
t426 = qJD(4) * t329;
t427 = qJD(4) * t324;
t143 = qJD(3) * t227 + t263 * t330 + t305 * t325;
t431 = qJD(2) * t321;
t396 = t326 * t431;
t423 = qJDD(1) * t321;
t262 = -qJD(1) * t396 + t331 * t423;
t248 = qJDD(3) - t262;
t173 = pkin(7) * t262 + t326 * t413 + t331 * t382;
t161 = pkin(8) * t305 + t173;
t166 = -pkin(1) * t423 - pkin(2) * t262 - pkin(8) * t263;
t75 = -t141 * qJD(3) - t161 * t325 + t330 * t166;
t44 = pkin(3) * t248 - pkin(9) * t143 + t75;
t428 = qJD(3) * t330;
t74 = t330 * t161 + t325 * t166 - t208 * t429 + t217 * t428;
t48 = pkin(9) * t144 + t74;
t15 = t106 * t426 - t114 * t427 + t324 * t44 + t329 * t48;
t240 = qJDD(4) + t248;
t506 = t240 / 0.2e1;
t71 = -qJD(4) * t357 - t143 * t324 + t144 * t329;
t531 = t71 / 0.2e1;
t70 = qJD(4) * t387 + t143 * t329 + t144 * t324;
t532 = t70 / 0.2e1;
t69 = qJDD(5) - t71;
t533 = t69 / 0.2e1;
t41 = -qJD(5) * t123 + t240 * t328 - t323 * t70;
t537 = t41 / 0.2e1;
t40 = qJD(5) * t122 + t240 * t323 + t328 * t70;
t538 = t40 / 0.2e1;
t12 = pkin(10) * t240 + t15;
t18 = -pkin(4) * t71 - pkin(10) * t70 + t102;
t2 = qJD(5) * t27 + t12 * t328 + t18 * t323;
t3 = -qJD(5) * t28 - t12 * t323 + t18 * t328;
t551 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t8 = Ifges(6,5) * t40 + Ifges(6,6) * t41 + Ifges(6,3) * t69;
t589 = t551 + mrSges(5,1) * t102 - mrSges(5,3) * t15 + Ifges(6,5) * t538 + Ifges(6,6) * t537 + Ifges(6,3) * t533 + t8 / 0.2e1 + (-t506 - t240 / 0.2e1) * Ifges(5,6) + (-t531 - t71 / 0.2e1) * Ifges(5,2) + (-t532 - t70 / 0.2e1) * Ifges(5,4);
t455 = t114 * t324;
t62 = t106 * t329 - t455;
t554 = t157 * mrSges(5,2) - t62 * mrSges(5,3);
t150 = Ifges(5,4) * t387;
t476 = Ifges(5,5) * t279;
t98 = Ifges(5,1) * t357 + t150 + t476;
t588 = t554 + Ifges(5,1) * t509 + Ifges(5,4) * t511 + Ifges(5,5) * t503 + t98 / 0.2e1;
t513 = t151 / 0.2e1;
t519 = t123 / 0.2e1;
t521 = t122 / 0.2e1;
t586 = Ifges(6,5) * t519 + Ifges(6,6) * t521 + Ifges(6,3) * t513 - t590 - t597 - t599;
t264 = t322 * t330 - t325 * t450;
t354 = pkin(3) * t264;
t585 = t173 * mrSges(3,2);
t212 = t292 * t324 - t293 * t329;
t570 = -qJD(4) * t212 + t600 * t324 + t601 * t329;
t16 = -qJD(4) * t63 - t324 * t48 + t329 * t44;
t584 = t16 * mrSges(5,1) - t15 * mrSges(5,2);
t130 = mrSges(5,1) * t279 - mrSges(5,3) * t357;
t81 = -mrSges(6,1) * t122 + mrSges(6,2) * t123;
t583 = t81 - t130;
t320 = qJ(3) + qJ(4);
t317 = sin(t320);
t318 = cos(t320);
t374 = -mrSges(4,1) * t330 + mrSges(4,2) * t325;
t582 = -m(4) * pkin(2) + t317 * t385 + t318 * t596 + t374;
t500 = cos(qJ(1));
t399 = t500 * t331;
t327 = sin(qJ(1));
t444 = t326 * t327;
t269 = -t322 * t444 + t399;
t448 = t321 * t330;
t201 = -t269 * t325 + t327 * t448;
t581 = t75 * mrSges(4,1) - t74 * mrSges(4,2);
t504 = -t279 / 0.2e1;
t510 = -t357 / 0.2e1;
t512 = -t387 / 0.2e1;
t580 = -Ifges(5,4) * t510 - Ifges(5,2) * t512 - Ifges(5,6) * t504 + t590 + t598;
t578 = t102 * mrSges(5,2) - mrSges(5,3) * t16 + 0.2e1 * Ifges(5,1) * t532 + 0.2e1 * Ifges(5,4) * t531 + 0.2e1 * Ifges(5,5) * t506;
t575 = m(5) + m(6);
t516 = t143 / 0.2e1;
t515 = t144 / 0.2e1;
t505 = t248 / 0.2e1;
t316 = pkin(3) * t330 + pkin(2);
t180 = pkin(4) * t274 - pkin(10) * t275 - t316;
t125 = t180 * t328 - t212 * t323;
t574 = qJD(5) * t125 + t323 * t591 + t328 * t592;
t126 = t180 * t323 + t212 * t328;
t573 = -qJD(5) * t126 - t323 * t592 + t328 * t591;
t572 = -m(4) * pkin(8) - mrSges(4,3) - mrSges(5,3);
t569 = pkin(4) * t398 - t570;
t447 = t321 * t331;
t273 = t322 * t326 * pkin(1) + pkin(7) * t447;
t244 = pkin(8) * t322 + t273;
t164 = -t244 * t325 + t330 * t245;
t265 = t322 * t325 + t326 * t448;
t119 = -pkin(3) * t447 - pkin(9) * t265 + t164;
t165 = t330 * t244 + t325 * t245;
t131 = pkin(9) * t264 + t165;
t567 = t324 * t119 + t329 * t131;
t384 = mrSges(3,3) * t398;
t566 = -mrSges(3,1) * t306 - mrSges(4,1) * t227 + mrSges(4,2) * t228 + t384;
t169 = t210 * t323 + t328 * t398;
t424 = qJD(5) * t328;
t347 = -t323 * t196 + t275 * t424;
t565 = t169 + t347;
t170 = -t210 * t328 + t323 * t398;
t425 = qJD(5) * t323;
t346 = t328 * t196 + t275 * t425;
t564 = t170 + t346;
t499 = pkin(1) * t331;
t272 = t322 * t499 - t307;
t241 = -t317 * t450 + t318 * t322;
t242 = t317 * t322 + t318 * t450;
t562 = -t593 * t241 + t242 * t594;
t449 = t321 * t327;
t193 = t269 * t317 - t318 * t449;
t194 = t269 * t318 + t317 * t449;
t561 = t593 * t193 + t194 * t594;
t400 = t500 * t326;
t443 = t327 * t331;
t267 = t322 * t400 + t443;
t401 = t321 * t500;
t189 = -t267 * t317 - t318 * t401;
t190 = t267 * t318 - t317 * t401;
t560 = -t593 * t189 + t190 * t594;
t19 = mrSges(6,1) * t69 - mrSges(6,3) * t40;
t20 = -mrSges(6,2) * t69 + mrSges(6,3) * t41;
t558 = -t323 * t19 + t328 * t20;
t557 = -t325 * t75 + t330 * t74;
t371 = mrSges(6,1) * t323 + mrSges(6,2) * t328;
t555 = -t371 + t572;
t553 = mrSges(3,2) + t555;
t258 = qJD(2) * t350;
t260 = t272 * qJD(2);
t108 = -qJD(3) * t165 + t330 * t258 - t260 * t325;
t395 = t321 * t430;
t199 = qJD(3) * t264 + t330 * t395;
t87 = pkin(3) * t396 - pkin(9) * t199 + t108;
t107 = -t244 * t429 + t245 * t428 + t325 * t258 + t330 * t260;
t198 = -qJD(3) * t265 - t325 * t395;
t94 = pkin(9) * t198 + t107;
t24 = -qJD(4) * t567 - t324 * t94 + t329 * t87;
t552 = mrSges(3,1) - t582;
t104 = pkin(4) * t357 - pkin(10) * t387;
t59 = -pkin(4) * t279 - t62;
t348 = t59 * t371;
t120 = Ifges(6,4) * t122;
t58 = Ifges(6,1) * t123 + Ifges(6,5) * t151 + t120;
t457 = t328 * t58;
t501 = t323 / 0.2e1;
t526 = -t98 / 0.2e1;
t479 = Ifges(6,4) * t123;
t57 = Ifges(6,2) * t122 + Ifges(6,6) * t151 + t479;
t549 = -t348 - t457 / 0.2e1 + t57 * t501 + t526 - t554;
t546 = t321 ^ 2;
t9 = t40 * Ifges(6,4) + t41 * Ifges(6,2) + t69 * Ifges(6,6);
t544 = t9 / 0.2e1;
t10 = t40 * Ifges(6,1) + t41 * Ifges(6,4) + t69 * Ifges(6,5);
t541 = t10 / 0.2e1;
t530 = Ifges(4,4) * t516 + Ifges(4,2) * t515 + Ifges(4,6) * t505;
t529 = Ifges(4,1) * t516 + Ifges(4,4) * t515 + Ifges(4,5) * t505;
t524 = pkin(1) * mrSges(3,1);
t523 = pkin(1) * mrSges(3,2);
t466 = t228 * Ifges(4,4);
t136 = t227 * Ifges(4,2) + t287 * Ifges(4,6) + t466;
t518 = t136 / 0.2e1;
t218 = Ifges(4,4) * t227;
t137 = t228 * Ifges(4,1) + t287 * Ifges(4,5) + t218;
t517 = t137 / 0.2e1;
t507 = t228 / 0.2e1;
t502 = t322 / 0.2e1;
t498 = pkin(3) * t228;
t497 = pkin(3) * t324;
t496 = pkin(3) * t329;
t486 = mrSges(6,3) * t323;
t485 = mrSges(6,3) * t328;
t484 = Ifges(3,4) * t326;
t483 = Ifges(3,4) * t331;
t482 = Ifges(4,4) * t325;
t481 = Ifges(4,4) * t330;
t478 = Ifges(6,4) * t323;
t477 = Ifges(6,4) * t328;
t471 = t140 * mrSges(4,3);
t470 = t141 * mrSges(4,3);
t469 = t387 * Ifges(5,6);
t468 = t357 * Ifges(5,5);
t467 = t227 * Ifges(4,6);
t465 = t228 * Ifges(4,5);
t464 = t279 * Ifges(5,3);
t463 = t287 * Ifges(4,3);
t462 = t306 * Ifges(3,5);
t461 = t306 * Ifges(3,6);
t452 = t275 * t323;
t451 = t275 * t328;
t435 = t500 * pkin(1) + pkin(7) * t449;
t420 = Ifges(5,5) * t70 + Ifges(5,6) * t71 + Ifges(5,3) * t240;
t412 = t323 * t447;
t410 = t325 * t449;
t408 = t328 * t447;
t405 = t457 / 0.2e1;
t404 = Ifges(4,5) * t143 + Ifges(4,6) * t144 + Ifges(4,3) * t248;
t403 = Ifges(3,5) * t263 + Ifges(3,6) * t262 + Ifges(3,3) * t305;
t390 = -t425 / 0.2e1;
t389 = -pkin(1) * t327 + pkin(7) * t401;
t388 = -t193 * pkin(4) + pkin(10) * t194;
t296 = t325 * t401;
t386 = -t267 * t330 + t296;
t383 = mrSges(3,3) * t397;
t377 = t201 * pkin(3);
t376 = mrSges(3,2) + t572;
t375 = mrSges(4,1) * t264 - mrSges(4,2) * t265;
t370 = Ifges(4,1) * t330 - t482;
t369 = Ifges(6,1) * t328 - t478;
t368 = Ifges(3,2) * t331 + t484;
t367 = -Ifges(4,2) * t325 + t481;
t366 = -Ifges(6,2) * t323 + t477;
t365 = Ifges(4,5) * t330 - Ifges(4,6) * t325;
t364 = Ifges(6,5) * t328 - Ifges(6,6) * t323;
t363 = t27 * t328 + t28 * t323;
t78 = -pkin(10) * t447 + t567;
t172 = t264 * t324 + t265 * t329;
t243 = t307 + (-pkin(2) - t499) * t322;
t175 = t243 - t354;
t356 = t329 * t264 - t265 * t324;
t95 = -pkin(4) * t356 - pkin(10) * t172 + t175;
t37 = t323 * t95 + t328 * t78;
t36 = -t323 * t78 + t328 * t95;
t268 = t322 * t443 + t400;
t360 = pkin(3) * t410 - t268 * t332 + t269 * t316 + t435;
t79 = t119 * t329 - t131 * t324;
t148 = -t172 * t323 - t408;
t349 = -t172 * t328 + t412;
t344 = t122 * t366;
t343 = t123 * t369;
t342 = t151 * t364;
t23 = t119 * t426 - t131 * t427 + t324 * t87 + t329 * t94;
t338 = t267 * t325 + t330 * t401;
t261 = t273 * qJD(2);
t337 = t338 * pkin(3);
t336 = -qJD(5) * t363 - t3 * t323;
t160 = -pkin(3) * t198 + t261;
t335 = t2 * t328 + t336;
t13 = -pkin(4) * t240 - t16;
t333 = t10 * t501 + t13 * t372 + t2 * t485 + t328 * t544 + (Ifges(6,1) * t323 + t477) * t538 + (Ifges(6,2) * t328 + t478) * t537 + t420 + t57 * t390 + (Ifges(6,5) * t323 + Ifges(6,6) * t328) * t533 + (t348 + t405) * qJD(5) + (t344 + t343 + t342) * qJD(5) / 0.2e1 + t584;
t315 = -pkin(4) - t496;
t301 = Ifges(3,4) * t397;
t270 = (-mrSges(3,1) * t331 + mrSges(3,2) * t326) * t321;
t266 = -t322 * t399 + t444;
t255 = -mrSges(3,2) * t306 + t383;
t232 = t241 * pkin(4);
t205 = Ifges(3,1) * t398 + t301 + t462;
t204 = t368 * t434 + t461;
t202 = t269 * t330 + t410;
t187 = t189 * pkin(4);
t179 = mrSges(4,1) * t287 - mrSges(4,3) * t228;
t178 = -mrSges(4,2) * t287 + mrSges(4,3) * t227;
t146 = t194 * t328 + t268 * t323;
t145 = -t194 * t323 + t268 * t328;
t135 = t463 + t465 + t467;
t129 = -mrSges(5,2) * t279 + mrSges(5,3) * t387;
t116 = -mrSges(4,2) * t248 + mrSges(4,3) * t144;
t115 = mrSges(4,1) * t248 - mrSges(4,3) * t143;
t103 = -mrSges(5,1) * t387 + mrSges(5,2) * t357;
t101 = qJD(4) * t172 - t329 * t198 + t199 * t324;
t100 = qJD(4) * t356 + t198 * t324 + t199 * t329;
t99 = -mrSges(4,1) * t144 + mrSges(4,2) * t143;
t96 = t464 + t468 + t469;
t91 = mrSges(6,1) * t151 - mrSges(6,3) * t123;
t90 = -mrSges(6,2) * t151 + mrSges(6,3) * t122;
t88 = t104 + t498;
t77 = pkin(4) * t447 - t79;
t73 = t113 * t329 - t455;
t72 = t113 * t324 + t442;
t66 = qJD(5) * t349 - t100 * t323 + t328 * t396;
t65 = qJD(5) * t148 + t100 * t328 + t323 * t396;
t55 = -mrSges(5,2) * t240 + mrSges(5,3) * t71;
t54 = mrSges(5,1) * t240 - mrSges(5,3) * t70;
t42 = pkin(4) * t101 - pkin(10) * t100 + t160;
t33 = t104 * t323 + t328 * t62;
t32 = t104 * t328 - t323 * t62;
t31 = t323 * t88 + t328 * t73;
t30 = -t323 * t73 + t328 * t88;
t29 = -mrSges(5,1) * t71 + mrSges(5,2) * t70;
t22 = -pkin(4) * t396 - t24;
t21 = pkin(10) * t396 + t23;
t17 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t5 = -qJD(5) * t37 - t21 * t323 + t328 * t42;
t4 = qJD(5) * t36 + t21 * t328 + t323 * t42;
t1 = [(mrSges(3,3) * t173 - Ifges(4,5) * t516 - Ifges(5,5) * t532 - Ifges(4,6) * t515 - Ifges(5,6) * t531 - Ifges(4,3) * t505 - Ifges(5,3) * t506 - t581 - t584) * t447 + (Ifges(4,1) * t265 + Ifges(4,4) * t264) * t516 + t227 * (Ifges(4,4) * t199 + Ifges(4,2) * t198) / 0.2e1 + (-t500 * mrSges(2,1) - m(4) * (pkin(2) * t269 + t435) - t202 * mrSges(4,1) - t201 * mrSges(4,2) - m(5) * t360 - t194 * mrSges(5,1) - m(6) * (pkin(4) * t194 + t360) - t146 * mrSges(6,1) - t145 * mrSges(6,2) - m(3) * t435 - t269 * mrSges(3,1) + (-mrSges(3,3) * t321 + mrSges(2,2)) * t327 + t385 * t193 + t376 * t268) * g(2) + (t272 * mrSges(3,1) - t273 * mrSges(3,2) + Ifges(3,3) * t502 + (Ifges(3,5) * t326 + Ifges(3,6) * t331) * t321) * t305 + t174 * (mrSges(3,1) * t322 - mrSges(3,3) * t450) - (t404 + t420) * t447 / 0.2e1 + m(5) * (t102 * t175 + t15 * t567 + t157 * t160 + t16 * t79 + t23 * t63 + t24 * t62) + t567 * t55 + (-Ifges(6,1) * t349 + Ifges(6,4) * t148) * t538 + t13 * (-mrSges(6,1) * t148 - mrSges(6,2) * t349) + (-Ifges(6,5) * t349 + Ifges(6,6) * t148) * t533 + (-Ifges(6,4) * t349 + Ifges(6,2) * t148) * t537 + (t148 * t2 - t27 * t65 + t28 * t66 + t3 * t349) * mrSges(6,3) - t349 * t541 + (-t140 * t199 + t141 * t198 + t264 * t74 - t265 * t75) * mrSges(4,3) + t287 * (Ifges(4,5) * t199 + Ifges(4,6) * t198) / 0.2e1 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t546 + t173 * t273 + t174 * t272 - t256 * t261 + t259 * t260) + (-m(3) * t389 + t267 * mrSges(3,1) - mrSges(3,3) * t401 + t327 * mrSges(2,1) + t500 * mrSges(2,2) - m(4) * (-pkin(2) * t267 + t389) - t386 * mrSges(4,1) - t338 * mrSges(4,2) + t385 * t189 - t596 * t190 + (t371 - t376) * t266 + t575 * (-pkin(3) * t296 - t266 * t332 + t267 * t316 - t389)) * g(1) + (-t272 * mrSges(3,3) + Ifges(3,5) * t502 + (t326 * Ifges(3,1) + t483 - t523) * t321) * t263 + ((t462 / 0.2e1 + t205 / 0.2e1 - t256 * mrSges(3,3) + (-t523 + t483 / 0.2e1) * t434) * t331 + (t464 / 0.2e1 + t465 / 0.2e1 + t467 / 0.2e1 - t141 * mrSges(4,2) + t140 * mrSges(4,1) + t468 / 0.2e1 + t469 / 0.2e1 - t63 * mrSges(5,2) + t62 * mrSges(5,1) + t463 / 0.2e1 - t461 / 0.2e1 + t96 / 0.2e1 + t135 / 0.2e1 - t204 / 0.2e1 - t259 * mrSges(3,3) + (-t524 - t484 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t331) * t434) * t326) * t431 + (t273 * mrSges(3,3) + Ifges(3,6) * t502 + (t368 + t524) * t321) * t262 + t588 * t100 - t589 * t356 + t586 * t101 + m(6) * (t13 * t77 + t2 * t37 + t22 * t59 + t27 * t5 + t28 * t4 + t3 * t36) + m(4) * (t107 * t141 + t108 * t140 + t162 * t243 + t164 * t75 + t165 * t74 + t207 * t261) + t578 * t172 + (-pkin(1) * t270 * t321 + Ifges(2,3)) * qJDD(1) + t260 * t255 - t322 * t585 + (Ifges(6,1) * t65 + Ifges(6,4) * t66) * t519 + t243 * t99 + t207 * (-mrSges(4,1) * t198 + mrSges(4,2) * t199) - t162 * t375 + t108 * t179 + t175 * t29 + t107 * t178 + t164 * t115 + t165 * t116 + t160 * t103 + t23 * t129 + t24 * t130 + (Ifges(6,5) * t65 + Ifges(6,6) * t66) * t513 + t4 * t90 + t5 * t91 + t22 * t81 + t77 * t17 + t79 * t54 + t65 * t58 / 0.2e1 + t59 * (-mrSges(6,1) * t66 + mrSges(6,2) * t65) + t66 * t57 / 0.2e1 + t566 * t261 + t403 * t502 + (Ifges(4,1) * t199 + Ifges(4,4) * t198) * t507 + t199 * t517 + t198 * t518 + t265 * t529 + t264 * t530 + t148 * t544 + (Ifges(4,4) * t265 + Ifges(4,2) * t264) * t515 + (Ifges(6,4) * t65 + Ifges(6,2) * t66) * t521 + (Ifges(4,5) * t265 + Ifges(4,6) * t264) * t505 + t36 * t19 + t37 * t20; -t565 * t57 / 0.2e1 + (mrSges(6,1) * t565 - mrSges(6,2) * t564) * t59 + (-t2 * t452 + t27 * t564 - t28 * t565 - t3 * t451) * mrSges(6,3) + (m(4) * ((-t140 * t330 - t141 * t325) * qJD(3) + t557) - t179 * t428 - t178 * t429 - t325 * t115 + t330 * t116) * pkin(8) + t557 * mrSges(4,3) + t559 * t103 + t580 * t209 + (t270 - t575 * t316 * t447 + (t582 * t331 + (t575 * t332 + t555) * t326) * t321) * g(3) + (-Ifges(5,5) * t210 + Ifges(5,3) * t398) * t504 - t62 * (mrSges(5,1) * t398 + mrSges(5,3) * t210) + (t157 * t210 + t398 * t63) * mrSges(5,2) + (-Ifges(5,4) * t210 + Ifges(5,6) * t398) * t512 + (-Ifges(5,1) * t210 + Ifges(5,5) * t398) * t510 + (t227 * t367 + t228 * t370 + t287 * t365) * qJD(3) / 0.2e1 + (-t326 * (Ifges(3,1) * t331 - t484) / 0.2e1 + pkin(1) * (mrSges(3,1) * t326 + mrSges(3,2) * t331)) * qJD(1) ^ 2 * t546 + (-Ifges(6,4) * t346 - Ifges(6,2) * t347) * t521 + (Ifges(6,4) * t170 + Ifges(6,2) * t169) * t522 - t9 * t452 / 0.2e1 + (-Ifges(6,5) * t346 - Ifges(6,6) * t347) * t513 + (Ifges(6,5) * t170 + Ifges(6,6) * t169) * t514 - ((-Ifges(3,2) * t398 + t330 * t137 + t205 + t301) * t331 + (t135 + t96) * t326 + t227 * (Ifges(4,6) * t326 + t331 * t367) + t306 * (Ifges(3,5) * t331 - Ifges(3,6) * t326) + t228 * (Ifges(4,5) * t326 + t331 * t370) + t287 * (Ifges(4,3) * t326 + t331 * t365)) * t434 / 0.2e1 - (t17 - t54) * t355 + (-t102 * t316 + t15 * t212 + t157 * t559 + t16 * t355 + t570 * t62 + t571 * t63) * m(5) + (t125 * t3 + t126 * t2 - t13 * t355 + t27 * t573 + t28 * t574 + t569 * t59) * m(6) + (-pkin(2) * t162 - t140 * t167 - t141 * t168) * m(4) + t403 + (-t470 - t136 / 0.2e1) * t429 - (t405 + t588) * t196 + t589 * t274 + (-t141 * (-mrSges(4,3) * t325 * t331 - mrSges(4,2) * t326) - t140 * (mrSges(4,1) * t326 - mrSges(4,3) * t441)) * t434 + t586 * t197 + t204 * t398 / 0.2e1 + (t13 * t371 + t364 * t533 + t366 * t537 + t369 * t538 + t390 * t58 + t578) * t275 + (-t471 + t517) * t428 - t585 + (-Ifges(6,1) * t346 - Ifges(6,4) * t347) * t519 + (Ifges(6,1) * t170 + Ifges(6,4) * t169) * t520 + (-t255 + t383) * t256 + t212 * t55 + t162 * t374 - t167 * t179 + t174 * mrSges(3,1) - t168 * t178 - t170 * t58 / 0.2e1 + t125 * t19 + t126 * t20 - pkin(2) * t99 + t287 * t207 * (mrSges(4,1) * t325 + mrSges(4,2) * t330) + (-m(4) * t207 + t384 - t566) * t259 - t316 * t29 + t569 * t81 + t570 * t130 + t571 * t129 + t573 * t91 + t574 * t90 + (-t575 * (-t268 * t316 - t269 * t332) + t553 * t269 + t552 * t268) * g(1) + (-t575 * (-t266 * t316 - t267 * t332) + t553 * t267 + t552 * t266) * g(2) + (Ifges(4,5) * t325 + Ifges(4,6) * t330) * t505 + (Ifges(4,2) * t330 + t482) * t515 + (Ifges(4,1) * t325 + t481) * t516 + t380 * t518 - t210 * t526 + t325 * t529 + t330 * t530 + t451 * t541; (m(5) * t337 - m(6) * (t190 * pkin(10) + t187 - t337) + mrSges(4,1) * t338 - mrSges(4,2) * t386 + t560) * g(2) + (-m(5) * t377 - m(6) * (t377 + t388) - mrSges(4,1) * t201 + mrSges(4,2) * t202 + t561) * g(1) + (-m(5) * t354 - t375 - m(6) * (pkin(10) * t242 + t232 + t354) + t562) * g(3) + (m(6) * t335 - t424 * t91 - t425 * t90 + t558) * (pkin(10) + t497) + (-t323 * t91 + t328 * t90 + t129) * pkin(3) * t426 + t580 * t357 - (-Ifges(4,2) * t228 + t137 + t218) * t227 / 0.2e1 + (Ifges(5,1) * t510 + Ifges(5,4) * t512 + Ifges(5,5) * t504 + t27 * t485 + t28 * t486 + t364 * t514 + t366 * t522 + t369 * t520 + t549) * t387 + t583 * (pkin(3) * t427 - t72) + (-t27 * t424 - t28 * t425) * mrSges(6,3) + (-t27 * t30 - t28 * t31 - t59 * t72 + t13 * t315 + (t324 * t59 + (-t27 * t323 + t28 * t328) * t329) * qJD(4) * pkin(3)) * m(6) - t228 * (Ifges(4,1) * t227 - t466) / 0.2e1 + t581 - t103 * t498 + t404 + (-t157 * t498 + t62 * t72 - t63 * t73 + (t15 * t324 + t16 * t329 + (-t324 * t62 + t329 * t63) * qJD(4)) * pkin(3)) * m(5) + t333 - t207 * (mrSges(4,1) * t228 + mrSges(4,2) * t227) + t141 * t179 - t140 * t178 - t73 * t129 - t31 * t90 - t30 * t91 + t228 * t470 + t227 * t471 + t315 * t17 + t54 * t496 + t55 * t497 - t3 * t486 + t136 * t507 - t287 * (Ifges(4,5) * t227 - Ifges(4,6) * t228) / 0.2e1; -t583 * t63 - m(6) * (t27 * t32 + t28 * t33 + t59 * t63) + ((-t323 * t90 - t328 * t91) * qJD(5) + (-g(2) * t190 - g(3) * t242) * m(6) + t558) * pkin(10) + (t590 + t595) * t357 + t336 * mrSges(6,3) + (-m(6) * t187 + t560) * g(2) + (-m(6) * t388 + t561) * g(1) + (-m(6) * t232 + t562) * g(3) + t333 + (-t476 / 0.2e1 - t150 / 0.2e1 - t342 / 0.2e1 - t344 / 0.2e1 - t343 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t357 + t363 * mrSges(6,3) + t549) * t387 + t335 * t542 - t62 * t129 - t33 * t90 - t32 * t91 - t13 * t543 - pkin(4) * t17; -t59 * (mrSges(6,1) * t123 + mrSges(6,2) * t122) + (Ifges(6,1) * t122 - t479) * t520 + t57 * t519 + (Ifges(6,5) * t122 - Ifges(6,6) * t123) * t514 - t27 * t90 + t28 * t91 - g(1) * (mrSges(6,1) * t145 - mrSges(6,2) * t146) - g(2) * ((-t190 * t323 + t266 * t328) * mrSges(6,1) + (-t190 * t328 - t266 * t323) * mrSges(6,2)) - g(3) * ((-t242 * t323 - t408) * mrSges(6,1) + (-t242 * t328 + t412) * mrSges(6,2)) + (t122 * t27 + t123 * t28) * mrSges(6,3) + t8 + (-Ifges(6,2) * t123 + t120 + t58) * t522 + t551;];
tau = t1;
