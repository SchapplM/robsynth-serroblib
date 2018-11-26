% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 18:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:02:03
% EndTime: 2018-11-23 18:02:32
% DurationCPUTime: 30.14s
% Computational Cost: add. (18310->889), mult. (45925->1197), div. (0->0), fcn. (35078->10), ass. (0->394)
t335 = cos(qJ(2));
t327 = sin(pkin(6));
t421 = qJD(1) * t327;
t398 = t335 * t421;
t303 = -qJD(3) + t398;
t330 = sin(qJ(3));
t334 = cos(qJ(3));
t440 = cos(pkin(6));
t385 = t440 * qJD(1);
t352 = t385 + qJD(2);
t345 = t334 * t352;
t331 = sin(qJ(2));
t399 = t331 * t421;
t245 = t330 * t399 - t345;
t477 = -t245 / 0.2e1;
t377 = pkin(1) * t385;
t261 = -pkin(8) * t399 + t335 * t377;
t349 = t327 * (pkin(2) * t331 - pkin(9) * t335);
t262 = qJD(1) * t349;
t185 = -t261 * t330 + t262 * t334;
t326 = t334 * pkin(4);
t418 = qJD(3) * t334;
t499 = pkin(4) + pkin(9);
t500 = pkin(3) + pkin(10);
t592 = -(t326 * t335 - t331 * t500) * t421 + t185 + t499 * t418;
t419 = qJD(3) * t330;
t382 = pkin(3) * t419 - qJD(4) * t330;
t264 = pkin(8) * t398 + t331 * t377;
t380 = t330 * t398;
t403 = pkin(3) * t380 + t264;
t591 = -t382 + t403 + t303 * (pkin(10) * t330 - qJ(4) * t334);
t590 = Ifges(5,6) * t477;
t246 = t330 * t352 + t334 * t399;
t474 = t246 / 0.2e1;
t329 = sin(qJ(5));
t333 = cos(qJ(5));
t429 = t333 * t335;
t232 = (-t329 * t331 + t330 * t429) * t421;
t417 = qJD(5) * t329;
t589 = -t334 * t417 + t232;
t226 = pkin(9) * t352 + t264;
t258 = (-pkin(2) * t335 - pkin(9) * t331 - pkin(1)) * t327;
t237 = qJD(1) * t258;
t161 = t226 * t330 - t237 * t334;
t521 = -qJD(4) - t161;
t143 = pkin(3) * t303 - t521;
t240 = qJD(5) + t246;
t228 = qJD(6) + t240;
t480 = t228 / 0.2e1;
t191 = t245 * t333 + t303 * t329;
t192 = t245 * t329 - t303 * t333;
t328 = sin(qJ(6));
t332 = cos(qJ(6));
t119 = t191 * t328 + t192 * t332;
t491 = t119 / 0.2e1;
t384 = t191 * t332 - t192 * t328;
t493 = t384 / 0.2e1;
t580 = Ifges(7,5) * t491 + Ifges(7,6) * t493 + Ifges(7,3) * t480;
t567 = t303 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t474 + t590 + t580;
t588 = -t143 * mrSges(5,1) - t161 * mrSges(4,3) - t567;
t478 = t240 / 0.2e1;
t484 = t192 / 0.2e1;
t486 = t191 / 0.2e1;
t581 = Ifges(6,5) * t484 + Ifges(6,6) * t486 + Ifges(6,3) * t478;
t582 = Ifges(4,4) * t477;
t587 = Ifges(4,1) * t474 + t581 + t582;
t388 = -qJ(4) * t330 - pkin(2);
t278 = -t334 * t500 + t388;
t304 = t499 * t330;
t416 = qJD(5) * t333;
t576 = -t278 * t417 + t304 * t416 + t329 * t592 - t333 * t591;
t586 = t329 * t591 + t333 * t592;
t433 = t329 * t335;
t233 = (t330 * t433 + t331 * t333) * t421;
t284 = t329 * t304;
t379 = t334 * t398;
t387 = pkin(11) * t334 - t278;
t434 = t329 * t330;
t585 = -pkin(5) * t379 + pkin(11) * t233 + (pkin(5) * t334 - pkin(11) * t434) * qJD(3) + (t333 * t387 - t284) * qJD(5) + t586;
t583 = t576 + (t333 * t419 - t589) * pkin(11);
t472 = -t303 / 0.2e1;
t579 = Ifges(5,5) / 0.2e1;
t547 = Ifges(5,1) + Ifges(4,3);
t546 = -Ifges(4,5) + Ifges(5,4);
t545 = Ifges(5,5) - Ifges(4,6);
t285 = t333 * t304;
t187 = pkin(5) * t330 + t329 * t387 + t285;
t208 = t278 * t333 + t284;
t430 = t333 * t334;
t193 = -pkin(11) * t430 + t208;
t109 = t187 * t332 - t193 * t328;
t578 = qJD(6) * t109 + t328 * t585 + t332 * t583;
t110 = t187 * t328 + t193 * t332;
t577 = -qJD(6) * t110 - t328 * t583 + t332 * t585;
t575 = -qJD(5) * t208 + t586;
t162 = t226 * t334 + t237 * t330;
t127 = -pkin(4) * t245 + t162;
t121 = t333 * t127;
t438 = qJ(4) * t245;
t144 = t246 * t500 + t438;
t463 = pkin(11) + t500;
t467 = pkin(11) * t246;
t574 = t463 * t417 + pkin(5) * t245 - t121 - (-t144 - t467) * t329;
t294 = t463 * t333;
t64 = t127 * t329 + t144 * t333;
t573 = qJD(5) * t294 + t333 * t467 + t64;
t186 = t261 * t334 + t262 * t330;
t173 = -qJ(4) * t399 - t186;
t152 = -pkin(4) * t380 - t173;
t401 = -pkin(5) * t333 - pkin(4);
t572 = -t152 + (-pkin(9) + t401) * t419 + t589 * pkin(5);
t571 = -t419 * t499 - t152;
t414 = qJD(1) * qJD(2);
t390 = t335 * t414;
t375 = t327 * t390;
t201 = qJD(3) * t246 + t330 * t375;
t391 = t331 * t414;
t376 = t327 * t391;
t409 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1;
t436 = t327 * t331;
t408 = t330 * t436;
t378 = qJD(3) * t408;
t200 = qJD(1) * t378 - qJD(3) * t345 - t334 * t375;
t483 = -t200 / 0.2e1;
t518 = Ifges(5,2) * t483;
t420 = qJD(2) * t327;
t397 = t331 * t420;
t358 = t500 * t397;
t263 = qJD(2) * t349;
t253 = qJD(1) * t263;
t400 = pkin(1) * t440;
t276 = -pkin(8) * t436 + t335 * t400;
t265 = t276 * qJD(2);
t254 = qJD(1) * t265;
t87 = -t226 * t418 - t237 * t419 + t253 * t334 - t254 * t330;
t61 = -pkin(4) * t200 - qJD(1) * t358 - t87;
t435 = t327 * t335;
t277 = pkin(8) * t435 + t331 * t400;
t266 = t277 * qJD(2);
t255 = qJD(1) * t266;
t339 = qJ(4) * t200 - qJD(4) * t246 + t255;
t67 = t201 * t500 + t339;
t350 = pkin(4) * t246 + t161;
t559 = qJD(4) + t350;
t94 = t303 * t500 + t559;
t225 = -pkin(2) * t352 - t261;
t340 = -t246 * qJ(4) + t225;
t99 = t245 * t500 + t340;
t13 = t329 * t61 + t333 * t67 + t416 * t94 - t417 * t99;
t47 = t329 * t94 + t333 * t99;
t14 = -qJD(5) * t47 - t329 * t67 + t333 * t61;
t46 = -t329 * t99 + t333 * t94;
t38 = -pkin(11) * t192 + t46;
t36 = pkin(5) * t240 + t38;
t39 = pkin(11) * t191 + t47;
t443 = t328 * t39;
t15 = t332 * t36 - t443;
t102 = qJD(5) * t191 + t201 * t329 + t333 * t376;
t6 = -pkin(5) * t200 - pkin(11) * t102 + t14;
t103 = -qJD(5) * t192 + t201 * t333 - t329 * t376;
t7 = pkin(11) * t103 + t13;
t2 = qJD(6) * t15 + t328 * t6 + t332 * t7;
t441 = t332 * t39;
t16 = t328 * t36 + t441;
t528 = qJD(6) * t16;
t3 = -t328 * t7 + t332 * t6 - t528;
t564 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t544 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t564;
t552 = -t201 / 0.2e1;
t554 = -Ifges(5,4) / 0.2e1;
t81 = pkin(3) * t201 + t339;
t82 = -pkin(3) * t376 - t87;
t569 = t82 * mrSges(5,1) - t87 * mrSges(4,3) - t81 * mrSges(5,3) + Ifges(5,6) * t552 - t409 * t201 + t376 * t554 + t518 + t544;
t431 = t332 * t333;
t353 = t328 * t329 - t431;
t354 = t328 * t333 + t329 * t332;
t520 = qJD(5) + qJD(6);
t210 = t520 * t354;
t347 = t354 * t246;
t424 = t347 + t210;
t415 = qJD(6) * t328;
t209 = -t328 * t417 - t329 * t415 + t431 * t520;
t527 = t353 * t246;
t563 = t209 - t527;
t568 = t15 * t424 - t16 * t563 - t2 * t354 + t3 * t353;
t115 = Ifges(7,4) * t384;
t457 = Ifges(7,4) * t119;
t481 = -t228 / 0.2e1;
t492 = -t119 / 0.2e1;
t494 = -t384 / 0.2e1;
t52 = Ifges(7,1) * t119 + Ifges(7,5) * t228 + t115;
t291 = t303 * qJ(4);
t104 = t127 - t291;
t73 = -pkin(5) * t191 + t104;
t566 = (Ifges(7,5) * t384 - Ifges(7,6) * t119) * t481 + (t119 * t16 + t15 * t384) * mrSges(7,3) + (-Ifges(7,2) * t119 + t115 + t52) * t494 - t73 * (mrSges(7,1) * t119 + mrSges(7,2) * t384) + (Ifges(7,1) * t384 - t457) * t492;
t139 = t245 * pkin(3) + t340;
t565 = t46 * mrSges(6,1) + t15 * mrSges(7,1) + t225 * mrSges(4,2) - t47 * mrSges(6,2) - t16 * mrSges(7,2) - t139 * mrSges(5,3) + Ifges(4,5) * t472 + t580 + t587;
t145 = t291 - t162;
t560 = t145 * mrSges(5,1) - t162 * mrSges(4,3);
t557 = -t225 * mrSges(4,1) + t139 * mrSges(5,2) + Ifges(4,6) * t472 + t303 * t579 + (Ifges(4,2) + Ifges(5,3)) * t477 + (Ifges(4,4) + Ifges(5,6)) * t474;
t475 = -t246 / 0.2e1;
t476 = t245 / 0.2e1;
t556 = -Ifges(5,2) * t475 - Ifges(5,6) * t476 - t546 * t472 + t565 + t587;
t34 = qJD(6) * t384 + t102 * t332 + t103 * t328;
t505 = t34 / 0.2e1;
t35 = -qJD(6) * t119 - t102 * t328 + t103 * t332;
t504 = t35 / 0.2e1;
t51 = Ifges(7,2) * t384 + Ifges(7,6) * t228 + t457;
t553 = t51 / 0.2e1;
t100 = Ifges(6,6) * t103;
t101 = Ifges(6,5) * t102;
t40 = -Ifges(6,3) * t200 + t100 + t101;
t32 = Ifges(7,6) * t35;
t33 = Ifges(7,5) * t34;
t8 = -Ifges(7,3) * t200 + t32 + t33;
t551 = t40 + t8;
t548 = t87 * mrSges(4,1);
t532 = t254 * mrSges(3,2);
t293 = t463 * t329;
t216 = t293 * t328 - t294 * t332;
t531 = qJD(6) * t216 + t328 * t574 - t332 * t573;
t217 = -t293 * t332 - t294 * t328;
t530 = -qJD(6) * t217 + t328 * t573 + t332 * t574;
t529 = pkin(5) * t416 - t246 * t401 - t521;
t257 = pkin(9) * t440 + t277;
t183 = -t257 * t330 + t258 * t334;
t172 = pkin(3) * t435 - t183;
t273 = t330 * t440 + t334 * t436;
t128 = pkin(4) * t273 + pkin(10) * t435 + t172;
t386 = t440 * t334;
t272 = -t386 + t408;
t256 = -pkin(2) * t440 - t276;
t344 = -qJ(4) * t273 + t256;
t140 = t272 * t500 + t344;
t66 = t128 * t329 + t140 * t333;
t525 = mrSges(3,1) * t352 - mrSges(4,1) * t245 - mrSges(4,2) * t246 - mrSges(3,3) * t399;
t523 = t200 * t546 + t201 * t545 + t376 * t547;
t522 = t13 * t329 + t14 * t333;
t412 = Ifges(4,5) / 0.2e1 + t554;
t517 = -t412 * t303 + t565 + t581 - t588;
t516 = t334 * t520;
t512 = -t161 * mrSges(4,1) - t162 * mrSges(4,2) + t143 * mrSges(5,2) - t264 * mrSges(3,3) - t145 * mrSges(5,3);
t360 = t329 * t46 - t333 * t47;
t459 = Ifges(6,4) * t329;
t365 = Ifges(6,2) * t333 + t459;
t458 = Ifges(6,4) * t333;
t367 = Ifges(6,1) * t329 + t458;
t370 = mrSges(6,1) * t333 - mrSges(6,2) * t329;
t454 = Ifges(6,6) * t333;
t456 = Ifges(6,5) * t329;
t470 = -t333 / 0.2e1;
t471 = -t329 / 0.2e1;
t479 = -t240 / 0.2e1;
t485 = -t192 / 0.2e1;
t487 = -t191 / 0.2e1;
t460 = Ifges(6,4) * t192;
t92 = Ifges(6,2) * t191 + Ifges(6,6) * t240 + t460;
t190 = Ifges(6,4) * t191;
t93 = Ifges(6,1) * t192 + Ifges(6,5) * t240 + t190;
t511 = mrSges(6,3) * t360 + t104 * t370 + (t454 + t456) * t479 + t365 * t487 + t367 * t485 + t470 * t92 + t471 * t93;
t510 = -Ifges(4,4) * t474 - Ifges(4,2) * t477 + Ifges(5,6) * t475 + Ifges(5,3) * t476 + t472 * t545 - t557;
t507 = Ifges(7,4) * t505 + Ifges(7,2) * t504 + Ifges(7,6) * t483;
t506 = Ifges(7,1) * t505 + Ifges(7,4) * t504 + Ifges(7,5) * t483;
t42 = Ifges(6,1) * t102 + Ifges(6,4) * t103 - Ifges(6,5) * t200;
t503 = t42 / 0.2e1;
t502 = t92 / 0.2e1;
t501 = t93 / 0.2e1;
t497 = t102 / 0.2e1;
t496 = t103 / 0.2e1;
t468 = pkin(9) * t330;
t462 = Ifges(3,4) * t331;
t461 = Ifges(3,4) * t335;
t455 = Ifges(3,2) * t331;
t445 = t261 * mrSges(3,3);
t442 = t331 * Ifges(3,1);
t439 = Ifges(3,6) * qJD(2);
t432 = t330 * t333;
t125 = -mrSges(6,1) * t191 + mrSges(6,2) * t192;
t204 = mrSges(5,1) * t245 + mrSges(5,3) * t303;
t428 = t125 - t204;
t153 = t353 * t516 + t354 * t419;
t168 = t232 * t328 + t233 * t332;
t427 = t153 - t168;
t154 = -t353 * t419 + t354 * t516;
t167 = t232 * t332 - t233 * t328;
t426 = t154 - t167;
t202 = mrSges(4,2) * t303 - mrSges(4,3) * t245;
t423 = t202 - t204;
t203 = -mrSges(4,1) * t303 - mrSges(4,3) * t246;
t205 = mrSges(5,1) * t246 - mrSges(5,2) * t303;
t422 = t203 - t205;
t184 = t257 * t334 + t258 * t330;
t305 = pkin(9) * t334 + t326;
t413 = -Ifges(4,1) / 0.2e1 - Ifges(5,2) / 0.2e1;
t411 = t579 - Ifges(4,6) / 0.2e1;
t410 = -Ifges(4,2) / 0.2e1 - Ifges(5,3) / 0.2e1;
t396 = t335 * t420;
t393 = Ifges(3,5) * t440;
t392 = Ifges(3,6) * t440;
t389 = t420 / 0.2e1;
t177 = -mrSges(5,1) * t200 + mrSges(5,2) * t376;
t65 = t128 * t333 - t140 * t329;
t372 = qJD(1) * t389;
t369 = mrSges(6,1) * t329 + mrSges(6,2) * t333;
t368 = Ifges(6,1) * t333 - t459;
t366 = -Ifges(6,2) * t329 + t458;
t364 = Ifges(6,5) * t333 - Ifges(6,6) * t329;
t348 = -t272 * t329 + t327 * t429;
t49 = pkin(5) * t273 + pkin(11) * t348 + t65;
t214 = t272 * t333 + t327 * t433;
t54 = pkin(11) * t214 + t66;
t22 = -t328 * t54 + t332 * t49;
t23 = t328 * t49 + t332 * t54;
t74 = -mrSges(6,1) * t200 - mrSges(6,3) * t102;
t75 = mrSges(6,2) * t200 + mrSges(6,3) * t103;
t359 = t329 * t75 + t333 * t74;
t357 = t143 * t334 + t145 * t330;
t146 = -mrSges(6,2) * t240 + mrSges(6,3) * t191;
t147 = mrSges(6,1) * t240 - mrSges(6,3) * t192;
t356 = t146 * t333 - t147 * t329;
t355 = t161 * t334 - t162 * t330;
t150 = t214 * t332 + t328 * t348;
t151 = t214 * t328 - t332 * t348;
t351 = t331 * t372;
t171 = qJ(4) * t435 - t184;
t106 = -t257 * t418 - t258 * t419 + t263 * t334 - t265 * t330;
t212 = -qJD(3) * t386 - t334 * t396 + t378;
t71 = -pkin(4) * t212 - t106 - t358;
t213 = qJD(3) * t273 + t330 * t396;
t341 = qJ(4) * t212 - qJD(4) * t273 + t266;
t80 = t213 * t500 + t341;
t20 = t128 * t416 - t140 * t417 + t329 * t71 + t333 * t80;
t86 = -t226 * t419 + t237 * t418 + t253 * t330 + t254 * t334;
t105 = -t257 * t419 + t258 * t418 + t263 * t330 + t265 * t334;
t141 = -pkin(4) * t272 - t171;
t77 = -qJ(4) * t376 + qJD(4) * t303 - t86;
t21 = -qJD(5) * t66 - t329 * t80 + t333 * t71;
t59 = -pkin(4) * t201 - t77;
t88 = -qJ(4) * t397 + qJD(4) * t435 - t105;
t343 = (t392 + (Ifges(3,2) * t335 + t462) * t327) * qJD(1);
t72 = -pkin(4) * t213 - t88;
t338 = t246 * t409 + t303 * t411 + t557 - t560;
t321 = pkin(5) * t329 + qJ(4);
t308 = Ifges(3,4) * t398;
t301 = Ifges(3,5) * t375;
t296 = -pkin(3) * t334 + t388;
t271 = pkin(5) * t430 + t305;
t270 = -qJ(4) * t418 + t382;
t268 = t354 * t334;
t267 = t353 * t334;
t260 = -mrSges(3,2) * t352 + mrSges(3,3) * t398;
t222 = Ifges(3,1) * t399 + Ifges(3,5) * t352 + t308;
t221 = t343 + t439;
t207 = -t278 * t329 + t285;
t189 = -qJ(4) * t379 + t403;
t182 = -mrSges(5,2) * t245 - mrSges(5,3) * t246;
t180 = pkin(3) * t246 + t438;
t179 = -mrSges(4,2) * t376 - mrSges(4,3) * t201;
t178 = mrSges(4,1) * t376 + mrSges(4,3) * t200;
t176 = mrSges(5,1) * t201 - mrSges(5,3) * t376;
t175 = -pkin(3) * t399 - t185;
t170 = pkin(3) * t272 + t344;
t159 = -Ifges(5,1) * t303 - Ifges(5,4) * t246 + Ifges(5,5) * t245;
t156 = Ifges(4,5) * t246 - Ifges(4,6) * t245 - Ifges(4,3) * t303;
t137 = qJD(5) * t348 + t213 * t333 - t329 * t397;
t136 = qJD(5) * t214 + t213 * t329 + t333 * t397;
t130 = mrSges(4,1) * t201 - mrSges(4,2) * t200;
t129 = -mrSges(5,2) * t201 + mrSges(5,3) * t200;
t114 = -Ifges(4,1) * t200 - t201 * Ifges(4,4) + Ifges(4,5) * t376;
t113 = -t200 * Ifges(4,4) - t201 * Ifges(4,2) + Ifges(4,6) * t376;
t111 = Ifges(5,5) * t376 + t200 * Ifges(5,6) + t201 * Ifges(5,3);
t97 = pkin(3) * t213 + t341;
t96 = -pkin(3) * t397 - t106;
t95 = -pkin(5) * t214 + t141;
t84 = mrSges(7,1) * t228 - mrSges(7,3) * t119;
t83 = -mrSges(7,2) * t228 + mrSges(7,3) * t384;
t63 = -t144 * t329 + t121;
t56 = -mrSges(7,1) * t384 + mrSges(7,2) * t119;
t53 = -mrSges(6,1) * t103 + mrSges(6,2) * t102;
t45 = -pkin(5) * t137 + t72;
t44 = -qJD(6) * t151 - t136 * t328 + t137 * t332;
t43 = qJD(6) * t150 + t136 * t332 + t137 * t328;
t41 = Ifges(6,4) * t102 + Ifges(6,2) * t103 - Ifges(6,6) * t200;
t37 = -pkin(5) * t103 + t59;
t29 = mrSges(7,2) * t200 + mrSges(7,3) * t35;
t28 = -mrSges(7,1) * t200 - mrSges(7,3) * t34;
t19 = t332 * t38 - t443;
t18 = -t328 * t38 - t441;
t17 = pkin(11) * t137 + t20;
t12 = -pkin(5) * t212 - pkin(11) * t136 + t21;
t11 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t5 = -qJD(6) * t23 + t12 * t332 - t17 * t328;
t4 = qJD(6) * t22 + t12 * t328 + t17 * t332;
t1 = [(Ifges(7,4) * t43 + Ifges(7,2) * t44) * t493 + (Ifges(7,4) * t151 + Ifges(7,2) * t150) * t504 + t200 * (-Ifges(5,4) * t435 + Ifges(5,6) * t272) / 0.2e1 - (t343 + t221) * t397 / 0.2e1 + ((-t455 + t461) * t390 + (Ifges(3,1) * t335 - t462) * t391 - 0.2e1 * pkin(1) * (mrSges(3,1) * t331 + mrSges(3,2) * t335) * t414) * t327 ^ 2 - t396 * t445 + m(3) * (t254 * t277 + t264 * t265) + m(4) * (t105 * t162 - t106 * t161 + t183 * t87 + t184 * t86) + (Ifges(6,5) * t497 + Ifges(7,5) * t505 + Ifges(6,6) * t496 + Ifges(7,6) * t504 - t351 * t546 + t518 + t569) * t273 + m(5) * (t139 * t97 + t143 * t96 + t145 * t88 + t170 * t81 + t171 * t77 + t172 * t82) + m(6) * (t104 * t72 + t13 * t66 + t14 * t65 + t141 * t59 + t20 * t47 + t21 * t46) + m(7) * (t15 * t5 + t16 * t4 + t2 * t23 + t22 * t3 + t37 * t95 + t45 * t73) + (-t15 * t43 + t150 * t2 - t151 * t3 + t16 * t44) * mrSges(7,3) + t201 * (-Ifges(5,5) * t435 + Ifges(5,3) * t272) / 0.2e1 + (t254 * t435 + t255 * t436 - t276 * t375 - t277 * t376) * mrSges(3,3) + ((t156 + t159) * t331 + t352 * (Ifges(3,5) * t335 - Ifges(3,6) * t331) + t335 * t222) * t389 + (-m(3) * t276 + m(4) * t256 - mrSges(3,1) * t440 + mrSges(4,1) * t272 + mrSges(4,2) * t273) * t255 + (t510 + t560) * t213 + (Ifges(6,1) * t136 + Ifges(6,4) * t137) * t484 + (-t556 + t588) * t212 + (Ifges(6,4) * t136 + Ifges(6,2) * t137) * t486 + (Ifges(7,5) * t43 + Ifges(7,6) * t44) * t480 + (Ifges(6,5) * t136 + Ifges(6,6) * t137) * t478 + (t13 * t214 - t136 * t46 + t137 * t47 + t14 * t348) * mrSges(6,3) + (-Ifges(6,1) * t348 + Ifges(6,4) * t214) * t497 + (-Ifges(4,4) * t272 - Ifges(4,5) * t435 - Ifges(6,5) * t348 + Ifges(7,5) * t151 + Ifges(6,6) * t214 + Ifges(7,6) * t150 + (Ifges(4,1) + Ifges(6,3) + Ifges(7,3)) * t273) * t483 + (-Ifges(6,4) * t348 + Ifges(6,2) * t214) * t496 + t59 * (-mrSges(6,1) * t214 - mrSges(6,2) * t348) - t348 * t503 + (Ifges(7,1) * t43 + Ifges(7,4) * t44) * t491 + (Ifges(7,1) * t151 + Ifges(7,4) * t150) * t505 + t335 * (t393 + (t442 + t461) * t327) * t372 + (-t272 * t81 - t435 * t82) * mrSges(5,2) + t20 * t146 + t21 * t147 + t104 * (-mrSges(6,1) * t137 + mrSges(6,2) * t136) + t141 * t53 + t72 * t125 + t77 * (mrSges(5,1) * t272 + mrSges(5,3) * t435) + t86 * (mrSges(4,2) * t435 - mrSges(4,3) * t272) + t95 * t11 + t4 * t83 + t5 * t84 + t65 * t74 + t66 * t75 + t73 * (-mrSges(7,1) * t44 + mrSges(7,2) * t43) + t45 * t56 + t43 * t52 / 0.2e1 + t440 * (-Ifges(3,6) * t376 + t301) / 0.2e1 + t22 * t28 + t23 * t29 + t37 * (-mrSges(7,1) * t150 + mrSges(7,2) * t151) + t170 * t129 + t171 * t176 + t172 * t177 + t97 * t182 + t183 * t178 + t184 * t179 + t105 * t202 + t106 * t203 + t88 * t204 + t96 * t205 - t523 * t435 / 0.2e1 + (-m(3) * t261 + m(4) * t225 - t525) * t266 + (-Ifges(4,2) * t272 - Ifges(4,6) * t435) * t552 + t44 * t553 + t214 * t41 / 0.2e1 - t440 * t532 + t256 * t130 + t265 * t260 + t272 * t111 / 0.2e1 - t272 * t113 / 0.2e1 + (Ifges(5,4) * t475 + Ifges(4,5) * t474 + Ifges(5,5) * t476 + Ifges(4,6) * t477 + t472 * t547 + t512) * t397 + (t272 * t545 - t435 * t547) * t351 - t435 * t548 + (t114 + t551) * t273 / 0.2e1 + t136 * t501 + t137 * t502 + t151 * t506 + t150 * t507; (t153 / 0.2e1 - t168 / 0.2e1) * t52 + t575 * t147 + (t104 * t571 + t13 * t208 + t14 * t207 + t305 * t59 + t46 * t575 + t47 * t576) * m(6) + t576 * t146 + (t255 * mrSges(4,2) + t101 / 0.2e1 + t100 / 0.2e1 + t114 / 0.2e1 + t40 / 0.2e1 + t8 / 0.2e1 + t33 / 0.2e1 + t32 / 0.2e1 + (t177 - t178) * pkin(9) + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t413) * t200 + t569) * t330 + t571 * t125 + t572 * t56 + (t154 / 0.2e1 - t167 / 0.2e1) * t51 + m(4) * (-pkin(2) * t255 - t468 * t87) + m(5) * (t139 * t270 + t296 * t81 + t468 * t82) - m(4) * (-t161 * t185 + t162 * t186 + t225 * t264) - m(5) * (t139 * t189 + t143 * t175 + t145 * t173) + (-mrSges(7,1) * t426 + mrSges(7,2) * t427) * t73 + (-t15 * t427 + t16 * t426 + t2 * t267 + t268 * t3) * mrSges(7,3) + t37 * (-mrSges(7,1) * t267 - mrSges(7,2) * t268) + (-Ifges(7,5) * t268 + Ifges(7,6) * t267) * t483 + (-Ifges(7,4) * t268 + Ifges(7,2) * t267) * t504 + (-Ifges(7,1) * t268 + Ifges(7,4) * t267) * t505 + (-t189 + t270) * t182 + t577 * t84 + (t109 * t3 + t110 * t2 + t15 * t577 + t16 * t578 + t271 * t37 + t572 * t73) * m(7) + t578 * t83 + (-t232 * t47 + t233 * t46) * mrSges(6,3) + (t41 * t470 + t42 * t471 - t255 * mrSges(4,1) + t81 * mrSges(5,2) - t77 * mrSges(5,1) + t86 * mrSges(4,3) + t59 * t370 - t102 * t367 / 0.2e1 - t103 * t365 / 0.2e1 + t113 / 0.2e1 - t111 / 0.2e1 + t410 * t201 + (t456 / 0.2e1 + t454 / 0.2e1 - t409) * t200 + (-t13 * t333 + t14 * t329) * mrSges(6,3) + (m(4) * t86 - m(5) * t77 - t176 + t179) * pkin(9) + (-t104 * t369 + t329 * t502 + t93 * t470 + t366 * t487 + t368 * t485 + t364 * t479 + (t329 * t47 + t333 * t46) * mrSges(6,3)) * qJD(5)) * t334 - pkin(2) * t130 - t532 + ((m(4) * t355 + m(5) * t357) * pkin(9) + t355 * mrSges(4,3) + t357 * mrSges(5,1) + t104 * (-mrSges(6,1) * t432 + mrSges(6,2) * t434) + (Ifges(6,5) * t434 + Ifges(6,6) * t432) * t478 + (Ifges(6,1) * t434 + Ifges(6,4) * t432) * t484 + (Ifges(6,4) * t434 + Ifges(6,2) * t432) * t486 + t434 * t501 + t432 * t502 + (t432 * t47 - t434 * t46) * mrSges(6,3) + (-pkin(9) * t423 + t510) * t330 + (-t422 * pkin(9) + t556 + t567) * t334) * qJD(3) + t301 + t109 * t28 + t110 * t29 - t186 * t202 - t185 * t203 - t173 * t204 - t175 * t205 + t207 * t74 + t208 * t75 + t525 * t264 - t232 * t92 / 0.2e1 - t233 * t93 / 0.2e1 - t104 * (-mrSges(6,1) * t232 + mrSges(6,2) * t233) - t255 * mrSges(3,1) - t261 * t260 + t271 * t11 + t296 * t129 + ((-t439 / 0.2e1 + t221 / 0.2e1 - t156 / 0.2e1 - t159 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t303 - t412 * t246 - t411 * t245 + (t392 / 0.2e1 + (pkin(1) * mrSges(3,1) + t462 / 0.2e1) * t327) * qJD(1) + (-t546 * t330 - t545 * t334) * qJD(2) / 0.2e1 - t512) * t331 + (-qJD(2) * Ifges(3,5) / 0.2e1 + ((-t442 / 0.2e1 + pkin(1) * mrSges(3,2) + t455 / 0.2e1) * t327 - t393 / 0.2e1) * qJD(1) + t445 - t222 / 0.2e1 - t308 / 0.2e1 + (t245 * t410 + t338) * t330 + (t245 * t409 + t246 * t413 - t517) * t334) * t335) * t421 + t305 * t53 + (Ifges(6,5) * t233 + Ifges(6,6) * t232) * t479 + (Ifges(7,5) * t153 + Ifges(7,6) * t154) * t480 + (Ifges(7,5) * t168 + Ifges(7,6) * t167) * t481 + (Ifges(6,1) * t233 + Ifges(6,4) * t232) * t485 + (Ifges(6,4) * t233 + Ifges(6,2) * t232) * t487 + (Ifges(7,1) * t153 + Ifges(7,4) * t154) * t491 + (Ifges(7,1) * t168 + Ifges(7,4) * t167) * t492 + (Ifges(7,4) * t153 + Ifges(7,2) * t154) * t493 + (Ifges(7,4) * t168 + Ifges(7,2) * t167) * t494 - t268 * t506 + t267 * t507; t568 * mrSges(7,3) - t522 * mrSges(6,3) + (mrSges(7,1) * t563 - mrSges(7,2) * t424) * t73 + t428 * qJD(4) + t422 * t162 + t423 * t161 + (t59 * qJ(4) + t559 * t104 - t46 * t63 - t47 * t64) * m(6) + (t53 - t176) * qJ(4) + t350 * t125 + (Ifges(7,5) * t347 - Ifges(7,6) * t527) * t481 + (Ifges(7,1) * t347 - Ifges(7,4) * t527) * t492 + (Ifges(7,4) * t347 - Ifges(7,2) * t527) * t494 + (t527 / 0.2e1 - t209 / 0.2e1) * t51 - t354 * t507 + (-Ifges(7,5) * t353 - Ifges(7,6) * t354 + t364) * t483 + t37 * (mrSges(7,1) * t354 - mrSges(7,2) * t353) + (-Ifges(7,4) * t353 - Ifges(7,2) * t354) * t504 + (-Ifges(7,1) * t353 - Ifges(7,4) * t354) * t505 - t353 * t506 - ((-m(6) * t360 + t356) * qJD(5) + m(6) * t522 + t359) * t500 + (-t347 / 0.2e1 - t210 / 0.2e1) * t52 + (-Ifges(7,5) * t210 - Ifges(7,6) * t209) * t480 + (-Ifges(7,1) * t210 - Ifges(7,4) * t209) * t491 + (-Ifges(7,4) * t210 - Ifges(7,2) * t209) * t493 + (t338 + (t410 - t413) * t245 + t511) * t246 + t511 * qJD(5) + t523 - t64 * t146 - t63 * t147 - t86 * mrSges(4,2) + t82 * mrSges(5,2) - t77 * mrSges(5,3) + t59 * t369 + t548 - pkin(3) * t177 - t180 * t182 + (-pkin(3) * t82 - qJ(4) * t77 - t139 * t180 - t143 * t162 + t145 * t521) * m(5) + t216 * t28 + t217 * t29 + (t582 + t590 + t517) * t245 + t529 * t56 + t530 * t84 + t531 * t83 + (t15 * t530 + t16 * t531 + t2 * t217 + t216 * t3 + t321 * t37 + t529 * t73) * m(7) + t321 * t11 + t41 * t471 + t366 * t496 + t368 * t497 + t333 * t503; -t353 * t28 + t354 * t29 - t424 * t84 + t563 * t83 + t356 * qJD(5) + (t56 + t428) * t303 + (t182 + t356) * t246 + t359 + t177 + (t303 * t73 - t568) * m(7) + (t104 * t303 - t240 * t360 + t522) * m(6) + (t139 * t246 - t145 * t303 + t82) * m(5); -m(7) * (t15 * t18 + t16 * t19) + t119 * t553 - t46 * t146 + t47 * t147 + (-Ifges(6,2) * t192 + t190 + t93) * t487 + t551 - t19 * t83 - t18 * t84 + t544 - t104 * (mrSges(6,1) * t192 + mrSges(6,2) * t191) + (t191 * t46 + t192 * t47) * mrSges(6,3) + (-t192 * t56 + t332 * t28 + t328 * t29 + (-t328 * t84 + t332 * t83) * qJD(6) + (-t15 * t415 - t192 * t73 + t2 * t328 + (t3 + t528) * t332) * m(7)) * pkin(5) + (Ifges(6,5) * t191 - Ifges(6,6) * t192) * t479 + t92 * t484 + (Ifges(6,1) * t191 - t460) * t485 + t566; -t15 * t83 + t16 * t84 + t51 * t491 + t564 + t566 + t8;];
tauc  = t1(:);
