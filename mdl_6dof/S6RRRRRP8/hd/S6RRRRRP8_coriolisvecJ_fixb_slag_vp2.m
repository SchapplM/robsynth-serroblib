% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:31:50
% EndTime: 2018-11-23 18:32:23
% DurationCPUTime: 33.94s
% Computational Cost: add. (25317->913), mult. (65000->1246), div. (0->0), fcn. (51114->10), ass. (0->378)
t356 = cos(qJ(2));
t352 = sin(qJ(2));
t347 = sin(pkin(6));
t442 = qJD(1) * t347;
t419 = t352 * t442;
t348 = cos(pkin(6));
t441 = qJD(1) * t348;
t430 = pkin(1) * t441;
t302 = -pkin(8) * t419 + t356 * t430;
t375 = (pkin(2) * t352 - pkin(9) * t356) * t347;
t303 = qJD(1) * t375;
t351 = sin(qJ(3));
t355 = cos(qJ(3));
t239 = -t351 * t302 + t355 * t303;
t525 = -pkin(10) - pkin(9);
t420 = qJD(3) * t525;
t578 = -(-pkin(10) * t355 * t356 + pkin(3) * t352) * t442 - t239 + t355 * t420;
t240 = t355 * t302 + t351 * t303;
t415 = t356 * t442;
t409 = t351 * t415;
t577 = -pkin(10) * t409 - t351 * t420 + t240;
t350 = sin(qJ(4));
t354 = cos(qJ(4));
t317 = t350 * t355 + t351 * t354;
t432 = qJD(3) + qJD(4);
t258 = t432 * t317;
t269 = t317 * t415;
t573 = t258 - t269;
t331 = t525 * t351;
t332 = t525 * t355;
t535 = t354 * t331 + t332 * t350;
t564 = qJD(4) * t535 + t578 * t350 - t577 * t354;
t305 = pkin(8) * t415 + t352 * t430;
t259 = pkin(3) * t409 + t305;
t438 = qJD(3) * t351;
t576 = pkin(3) * t438 - t259;
t575 = pkin(11) * t419 - t564;
t316 = t350 * t351 - t354 * t355;
t257 = t432 * t316;
t270 = t316 * t415;
t574 = t576 + (t257 - t270) * pkin(11) + t573 * pkin(4);
t548 = Ifges(6,1) + Ifges(7,1);
t547 = Ifges(6,5) + Ifges(7,4);
t272 = t331 * t350 - t332 * t354;
t563 = -qJD(4) * t272 + t577 * t350 + t578 * t354;
t546 = Ifges(7,5) - Ifges(6,4);
t572 = Ifges(6,6) - Ifges(7,6);
t349 = sin(qJ(5));
t353 = cos(qJ(5));
t345 = -pkin(3) * t355 - pkin(2);
t255 = pkin(4) * t316 - pkin(11) * t317 + t345;
t536 = t349 * t255 + t353 * t272;
t567 = -qJD(5) * t536 + t349 * t575 + t574 * t353;
t433 = qJD(5) * t353;
t434 = qJD(5) * t349;
t566 = t255 * t433 - t272 * t434 + t574 * t349 - t353 * t575;
t562 = pkin(4) * t419 - t563;
t337 = qJD(2) + t441;
t268 = pkin(9) * t337 + t305;
t298 = (-pkin(2) * t356 - pkin(9) * t352 - pkin(1)) * t347;
t279 = qJD(1) * t298;
t219 = -t268 * t351 + t355 * t279;
t286 = t337 * t351 + t355 * t419;
t186 = -pkin(10) * t286 + t219;
t330 = qJD(3) - t415;
t170 = pkin(3) * t330 + t186;
t220 = t268 * t355 + t279 * t351;
t285 = t337 * t355 - t351 * t419;
t187 = pkin(10) * t285 + t220;
t178 = t350 * t187;
t107 = t354 * t170 - t178;
t376 = -qJD(4) - t330;
t103 = pkin(4) * t376 - t107;
t410 = t354 * t285 - t286 * t350;
t227 = qJD(5) - t410;
t179 = t354 * t187;
t108 = t350 * t170 + t179;
t104 = -pkin(11) * t376 + t108;
t267 = -t337 * pkin(2) - t302;
t232 = -t285 * pkin(3) + t267;
t378 = t285 * t350 + t354 * t286;
t117 = -pkin(4) * t410 - pkin(11) * t378 + t232;
t49 = t104 * t353 + t117 * t349;
t39 = qJ(6) * t227 + t49;
t197 = t349 * t378 + t353 * t376;
t198 = -t349 * t376 + t353 * t378;
t53 = t197 * pkin(5) - t198 * qJ(6) + t103;
t473 = Ifges(6,4) * t198;
t100 = -Ifges(6,2) * t197 + Ifges(6,6) * t227 + t473;
t194 = Ifges(7,5) * t198;
t97 = Ifges(7,6) * t227 + Ifges(7,3) * t197 + t194;
t555 = -t100 / 0.2e1 + t97 / 0.2e1;
t571 = t103 * mrSges(6,1) + t53 * mrSges(7,1) - t39 * mrSges(7,2) - t49 * mrSges(6,3) + t555;
t515 = -t197 / 0.2e1;
t195 = Ifges(6,4) * t197;
t467 = Ifges(7,5) * t197;
t538 = t198 * t548 + t547 * t227 - t195 + t467;
t570 = t538 / 0.2e1;
t569 = qJ(6) * t573 + qJD(6) * t316 + t566;
t568 = -pkin(5) * t573 - t567;
t241 = -t270 * t349 - t353 * t419;
t242 = -t270 * t353 + t349 * t419;
t388 = pkin(5) * t349 - qJ(6) * t353;
t389 = pkin(5) * t353 + qJ(6) * t349;
t565 = -pkin(5) * t241 + qJ(6) * t242 - t388 * t257 + (qJD(5) * t389 - qJD(6) * t353) * t317 + t562;
t439 = qJD(2) * t356;
t417 = t355 * t439;
t437 = qJD(3) * t355;
t251 = t337 * t437 + (-t352 * t438 + t417) * t442;
t418 = t351 * t439;
t252 = -t337 * t438 + (-t352 * t437 - t418) * t442;
t142 = qJD(4) * t378 + t251 * t350 - t354 * t252;
t141 = qJD(4) * t410 + t251 * t354 + t252 * t350;
t440 = qJD(2) * t347;
t413 = qJD(1) * t440;
t408 = t352 * t413;
t87 = -qJD(5) * t197 + t353 * t141 + t349 * t408;
t88 = qJD(5) * t198 + t349 * t141 - t353 * t408;
t560 = t547 * t142 + t546 * t88 + t548 * t87;
t304 = qJD(2) * t375;
t293 = qJD(1) * t304;
t450 = t347 * t352;
t338 = pkin(8) * t450;
t488 = pkin(1) * t356;
t312 = t348 * t488 - t338;
t306 = t312 * qJD(2);
t294 = qJD(1) * t306;
t162 = -qJD(3) * t220 + t355 * t293 - t294 * t351;
t116 = pkin(3) * t408 - pkin(10) * t251 + t162;
t161 = -t268 * t438 + t279 * t437 + t351 * t293 + t355 * t294;
t123 = pkin(10) * t252 + t161;
t435 = qJD(4) * t354;
t436 = qJD(4) * t350;
t31 = t116 * t354 - t350 * t123 - t170 * t436 - t187 * t435;
t28 = -pkin(4) * t408 - t31;
t33 = mrSges(6,1) * t88 + mrSges(6,2) * t87;
t559 = -m(6) * t28 - t33;
t558 = t349 * t547 + t353 * t572;
t465 = Ifges(7,5) * t353;
t471 = Ifges(6,4) * t353;
t557 = t349 * t548 - t465 + t471;
t30 = t350 * t116 + t354 * t123 + t170 * t435 - t187 * t436;
t27 = pkin(11) * t408 + t30;
t449 = t347 * t356;
t313 = t348 * t352 * pkin(1) + pkin(8) * t449;
t307 = t313 * qJD(2);
t295 = qJD(1) * t307;
t218 = -t252 * pkin(3) + t295;
t55 = t142 * pkin(4) - t141 * pkin(11) + t218;
t6 = -t104 * t434 + t117 * t433 + t353 * t27 + t349 * t55;
t7 = -qJD(5) * t49 - t27 * t349 + t353 * t55;
t556 = -t7 * t349 + t353 * t6;
t226 = Ifges(5,4) * t410;
t393 = Ifges(6,5) * t353 - Ifges(6,6) * t349;
t367 = t227 * t393;
t396 = Ifges(7,4) * t353 + Ifges(7,6) * t349;
t368 = t227 * t396;
t466 = Ifges(7,5) * t349;
t400 = Ifges(7,1) * t353 + t466;
t369 = t198 * t400;
t472 = Ifges(6,4) * t349;
t402 = Ifges(6,1) * t353 - t472;
t370 = t198 * t402;
t391 = Ifges(7,3) * t349 + t465;
t371 = t197 * t391;
t48 = -t104 * t349 + t117 * t353;
t539 = qJD(6) - t48;
t38 = -pkin(5) * t227 + t539;
t491 = -t353 / 0.2e1;
t492 = t349 / 0.2e1;
t493 = -t349 / 0.2e1;
t366 = Ifges(5,5) * t376;
t478 = Ifges(5,1) * t378;
t160 = t226 - t366 + t478;
t517 = -t160 / 0.2e1;
t398 = -Ifges(6,2) * t349 + t471;
t403 = mrSges(7,1) * t349 - mrSges(7,3) * t353;
t405 = mrSges(6,1) * t349 + mrSges(6,2) * t353;
t552 = -(t349 * t49 + t353 * t48) * mrSges(6,3) + t103 * t405 + t53 * t403 + t398 * t515;
t554 = (t349 * t39 - t353 * t38) * mrSges(7,2) + t107 * mrSges(5,3) + t517 - t371 / 0.2e1 - t370 / 0.2e1 - t369 / 0.2e1 - t368 / 0.2e1 - t367 / 0.2e1 - t232 * mrSges(5,2) + t100 * t492 + t97 * t493 + t538 * t491 - t552 - t226 / 0.2e1;
t2 = qJ(6) * t142 + qJD(6) * t227 + t6;
t4 = -pkin(5) * t142 - t7;
t553 = t2 * t353 + t349 * t4 + t38 * t433 - t39 * t434;
t551 = -t432 / 0.2e1;
t550 = t432 / 0.2e1;
t549 = pkin(5) * t378;
t545 = Ifges(7,2) + Ifges(6,3);
t20 = Ifges(6,5) * t87 - Ifges(6,6) * t88 + Ifges(6,3) * t142;
t21 = Ifges(7,4) * t87 + Ifges(7,2) * t142 + Ifges(7,6) * t88;
t543 = t21 + t20;
t308 = pkin(5) * t434 - qJ(6) * t433 - qJD(6) * t349;
t540 = t388 * t410;
t542 = t308 - t108 - t540;
t541 = qJ(6) * t378;
t297 = pkin(9) * t348 + t313;
t237 = -t351 * t297 + t355 * t298;
t310 = t348 * t351 + t355 * t450;
t193 = -pkin(3) * t449 - t310 * pkin(10) + t237;
t238 = t355 * t297 + t351 * t298;
t309 = t348 * t355 - t351 * t450;
t209 = pkin(10) * t309 + t238;
t129 = t350 * t193 + t354 * t209;
t120 = -pkin(11) * t449 + t129;
t245 = t309 * t350 + t310 * t354;
t296 = t338 + (-pkin(2) - t488) * t348;
t250 = -t309 * pkin(3) + t296;
t377 = t354 * t309 - t310 * t350;
t157 = -pkin(4) * t377 - t245 * pkin(11) + t250;
t537 = t353 * t120 + t349 * t157;
t534 = t366 / 0.2e1 + t554;
t457 = t286 * Ifges(4,4);
t216 = t285 * Ifges(4,2) + t330 * Ifges(4,6) + t457;
t280 = Ifges(4,4) * t285;
t217 = t286 * Ifges(4,1) + t330 * Ifges(4,5) + t280;
t380 = t219 * t355 + t220 * t351;
t475 = Ifges(4,4) * t355;
t476 = Ifges(4,4) * t351;
t489 = t355 / 0.2e1;
t494 = t330 / 0.2e1;
t497 = t286 / 0.2e1;
t499 = t285 / 0.2e1;
t533 = -t380 * mrSges(4,3) + t267 * (mrSges(4,1) * t351 + mrSges(4,2) * t355) + (-Ifges(4,2) * t351 + t475) * t499 + (Ifges(4,1) * t355 - t476) * t497 + (Ifges(4,5) * t355 - Ifges(4,6) * t351) * t494 - t351 * t216 / 0.2e1 + t217 * t489;
t172 = -qJD(3) * t238 + t355 * t304 - t306 * t351;
t260 = qJD(3) * t309 + t347 * t417;
t416 = t352 * t440;
t144 = pkin(3) * t416 - pkin(10) * t260 + t172;
t171 = -t297 * t438 + t298 * t437 + t351 * t304 + t355 * t306;
t261 = -qJD(3) * t310 - t347 * t418;
t155 = pkin(10) * t261 + t171;
t40 = t350 * t144 + t354 * t155 + t193 * t435 - t209 * t436;
t36 = pkin(11) * t416 + t40;
t164 = qJD(4) * t377 + t260 * t354 + t261 * t350;
t165 = qJD(4) * t245 + t260 * t350 - t354 * t261;
t235 = -t261 * pkin(3) + t307;
t70 = t165 * pkin(4) - t164 * pkin(11) + t235;
t13 = -qJD(5) * t537 - t349 * t36 + t353 * t70;
t167 = pkin(4) * t378 - pkin(11) * t410;
t426 = -Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1;
t427 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t428 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t474 = Ifges(5,4) * t378;
t365 = Ifges(5,6) * t376;
t464 = Ifges(5,2) * t410;
t159 = -t365 + t464 + t474;
t518 = t159 / 0.2e1;
t98 = t198 * Ifges(6,5) - t197 * Ifges(6,6) + t227 * Ifges(6,3);
t99 = t198 * Ifges(7,4) + t227 * Ifges(7,2) + t197 * Ifges(7,6);
t531 = t427 * t197 - t428 * t198 + t426 * t227 + t108 * mrSges(5,3) + t38 * mrSges(7,1) + t49 * mrSges(6,2) + t518 - t98 / 0.2e1 - t99 / 0.2e1 + t474 / 0.2e1 - t232 * mrSges(5,1) - t39 * mrSges(7,3) - t48 * mrSges(6,1);
t530 = Ifges(5,2) / 0.2e1;
t529 = t87 / 0.2e1;
t528 = -t88 / 0.2e1;
t527 = t88 / 0.2e1;
t523 = pkin(1) * mrSges(3,1);
t522 = pkin(1) * mrSges(3,2);
t519 = t142 / 0.2e1;
t514 = t197 / 0.2e1;
t513 = -t198 / 0.2e1;
t512 = t198 / 0.2e1;
t509 = -t227 / 0.2e1;
t508 = t227 / 0.2e1;
t506 = t377 / 0.2e1;
t504 = t245 / 0.2e1;
t503 = t251 / 0.2e1;
t502 = t252 / 0.2e1;
t501 = -t269 / 0.2e1;
t498 = -t286 / 0.2e1;
t496 = t309 / 0.2e1;
t495 = t310 / 0.2e1;
t490 = t353 / 0.2e1;
t487 = pkin(3) * t354;
t480 = mrSges(6,3) * t197;
t479 = mrSges(6,3) * t198;
t477 = Ifges(3,4) * t352;
t470 = Ifges(3,5) * t356;
t469 = Ifges(5,5) * t378;
t468 = Ifges(5,5) * t270;
t463 = Ifges(5,6) * t410;
t462 = Ifges(5,6) * t269;
t461 = t141 * Ifges(5,1);
t460 = t141 * Ifges(5,4);
t459 = t142 * Ifges(5,4);
t458 = t285 * Ifges(4,6);
t456 = t286 * Ifges(4,5);
t455 = t294 * mrSges(3,2);
t454 = t330 * Ifges(4,3);
t453 = t337 * Ifges(3,5);
t448 = t349 * t354;
t447 = t353 * t354;
t61 = t353 * t107 + t349 * t167;
t113 = t186 * t354 - t178;
t145 = pkin(3) * t286 + t167;
t59 = t353 * t113 + t349 * t145;
t148 = -mrSges(7,2) * t197 + mrSges(7,3) * t227;
t149 = -mrSges(6,2) * t227 - t480;
t446 = t148 + t149;
t150 = mrSges(6,1) * t227 - t479;
t151 = -mrSges(7,1) * t227 + mrSges(7,2) * t198;
t445 = -t150 + t151;
t132 = mrSges(6,1) * t197 + mrSges(6,2) * t198;
t208 = -mrSges(5,1) * t376 - mrSges(5,3) * t378;
t444 = t208 - t132;
t443 = -mrSges(3,1) * t337 - mrSges(4,1) * t285 + mrSges(4,2) * t286 + mrSges(3,3) * t419;
t423 = t349 * t449;
t422 = Ifges(5,5) * t141 - Ifges(5,6) * t142 + Ifges(5,3) * t408;
t421 = Ifges(4,5) * t251 + Ifges(4,6) * t252 + Ifges(4,3) * t408;
t45 = -t142 * mrSges(7,1) + t87 * mrSges(7,2);
t112 = t186 * t350 + t179;
t128 = t193 * t354 - t350 * t209;
t119 = pkin(4) * t449 - t128;
t406 = mrSges(6,1) * t353 - mrSges(6,2) * t349;
t404 = mrSges(7,1) * t353 + mrSges(7,3) * t349;
t397 = Ifges(6,2) * t353 + t472;
t394 = Ifges(5,5) * t164 - Ifges(5,6) * t165;
t390 = -Ifges(7,3) * t353 + t466;
t60 = -t107 * t349 + t167 * t353;
t58 = -t113 * t349 + t145 * t353;
t66 = -t120 * t349 + t157 * t353;
t382 = t161 * t355 - t162 * t351;
t202 = t255 * t353 - t272 * t349;
t323 = -pkin(4) - t389;
t41 = t144 * t354 - t350 * t155 - t193 * t436 - t209 * t435;
t224 = t349 * t245 + t353 * t449;
t12 = -t120 * t434 + t157 * t433 + t349 * t70 + t353 * t36;
t364 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t37 = -pkin(4) * t416 - t41;
t361 = -t365 / 0.2e1 + t531;
t44 = mrSges(6,1) * t142 - mrSges(6,3) * t87;
t46 = -mrSges(6,2) * t142 - mrSges(6,3) * t88;
t47 = -mrSges(7,2) * t88 + mrSges(7,3) * t142;
t360 = (t46 + t47) * t353 + (-t44 + t45) * t349 + (-t349 * t446 + t353 * t445) * qJD(5) + m(6) * (-t433 * t48 - t434 * t49 + t556) + m(7) * t553;
t11 = pkin(5) * t88 - qJ(6) * t87 - qJD(6) * t198 + t28;
t19 = t87 * Ifges(7,5) + t142 * Ifges(7,6) + t88 * Ifges(7,3);
t22 = t87 * Ifges(6,4) - t88 * Ifges(6,2) + t142 * Ifges(6,6);
t359 = t31 * mrSges(5,1) - t30 * mrSges(5,2) - t11 * t404 + t19 * t491 + t22 * t490 - t28 * t406 + t390 * t527 + t397 * t528 + t422 + t557 * t529 + t558 * t519 + t560 * t492 + t555 * t434 + t433 * t570 + t556 * mrSges(6,3) + t552 * qJD(5) + t553 * mrSges(7,2) + (t371 + t370 + t369 + t368 + t367) * qJD(5) / 0.2e1;
t333 = Ifges(3,4) * t415;
t329 = t413 * t470;
t314 = t323 - t487;
t301 = -t337 * mrSges(3,2) + mrSges(3,3) * t415;
t284 = pkin(3) * t436 + t308;
t265 = Ifges(3,1) * t419 + t333 + t453;
t264 = Ifges(3,6) * t337 + (Ifges(3,2) * t356 + t477) * t442;
t254 = mrSges(4,1) * t330 - mrSges(4,3) * t286;
t253 = -mrSges(4,2) * t330 + mrSges(4,3) * t285;
t234 = -mrSges(4,2) * t408 + mrSges(4,3) * t252;
t233 = mrSges(4,1) * t408 - mrSges(4,3) * t251;
t225 = t353 * t245 - t423;
t215 = t454 + t456 + t458;
t212 = t317 * t388 - t535;
t207 = mrSges(5,2) * t376 + mrSges(5,3) * t410;
t188 = -mrSges(4,1) * t252 + mrSges(4,2) * t251;
t185 = -pkin(5) * t316 - t202;
t184 = qJ(6) * t316 + t536;
t175 = t251 * Ifges(4,1) + t252 * Ifges(4,4) + Ifges(4,5) * t408;
t174 = t251 * Ifges(4,4) + t252 * Ifges(4,2) + Ifges(4,6) * t408;
t166 = -mrSges(5,1) * t410 + mrSges(5,2) * t378;
t158 = -Ifges(5,3) * t376 + t463 + t469;
t131 = mrSges(7,1) * t197 - mrSges(7,3) * t198;
t130 = pkin(5) * t198 + qJ(6) * t197;
t125 = -mrSges(5,2) * t408 - mrSges(5,3) * t142;
t124 = mrSges(5,1) * t408 - mrSges(5,3) * t141;
t110 = -qJD(5) * t423 + t164 * t349 + t245 * t433 - t353 * t416;
t109 = -qJD(5) * t224 + t353 * t164 + t349 * t416;
t71 = pkin(5) * t224 - qJ(6) * t225 + t119;
t69 = mrSges(5,1) * t142 + mrSges(5,2) * t141;
t65 = Ifges(5,5) * t408 - t459 + t461;
t64 = -t142 * Ifges(5,2) + Ifges(5,6) * t408 + t460;
t63 = t112 + t540;
t57 = pkin(5) * t377 - t66;
t56 = -qJ(6) * t377 + t537;
t52 = -t60 - t549;
t51 = t61 + t541;
t43 = -t58 - t549;
t42 = t59 + t541;
t32 = mrSges(7,1) * t88 - mrSges(7,3) * t87;
t14 = pkin(5) * t110 - qJ(6) * t109 - qJD(6) * t225 + t37;
t9 = -pkin(5) * t165 - t13;
t8 = qJ(6) * t165 - qJD(6) * t377 + t12;
t1 = [(t99 + t98) * t165 / 0.2e1 + t410 * (Ifges(5,4) * t164 - Ifges(5,2) * t165 + Ifges(5,6) * t416) / 0.2e1 + t378 * (Ifges(5,1) * t164 - Ifges(5,4) * t165 + Ifges(5,5) * t416) / 0.2e1 - (t421 + t422) * t449 / 0.2e1 + (t109 * t547 + t165 * t545) * t508 + (t109 * t548 + t165 * t547) * t512 + (Ifges(7,5) * t109 + Ifges(7,6) * t165) * t514 + (Ifges(6,4) * t109 + Ifges(6,6) * t165) * t515 + t31 * (-mrSges(5,1) * t449 - t245 * mrSges(5,3)) + t162 * (-mrSges(4,1) * t449 - t310 * mrSges(4,3)) + t161 * (mrSges(4,2) * t449 + t309 * mrSges(4,3)) + (-t356 * t394 / 0.2e1 + ((Ifges(3,5) * t348 / 0.2e1 - t312 * mrSges(3,3) + (-0.2e1 * t522 + 0.3e1 / 0.2e1 * Ifges(3,4) * t356) * t347) * t356 + (Ifges(4,5) * t495 + Ifges(4,6) * t496 + Ifges(5,5) * t504 + Ifges(5,6) * t506 - Ifges(3,6) * t348 - t313 * mrSges(3,3) + (-0.2e1 * t523 - 0.3e1 / 0.2e1 * t477) * t347 + (-Ifges(4,3) / 0.2e1 - Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t449) * t352) * qJD(2)) * t442 + (-Ifges(6,2) * t515 + Ifges(7,3) * t514 - t508 * t572 + t512 * t546 + t571) * t110 + (-t6 * mrSges(6,3) - t2 * mrSges(7,2) - t22 / 0.2e1 + t28 * mrSges(6,1) + t11 * mrSges(7,1) + t19 / 0.2e1 + Ifges(7,3) * t527 - Ifges(6,2) * t528 + t546 * t529 - t572 * t519) * t224 + t560 * t225 / 0.2e1 + (t329 / 0.2e1 - t455 - t295 * mrSges(3,1)) * t348 + t443 * t307 + m(5) * (t107 * t41 + t108 * t40 + t128 * t31 + t129 * t30 + t218 * t250 + t232 * t235) + m(7) * (t11 * t71 + t14 * t53 + t2 * t56 + t38 * t9 + t39 * t8 + t4 * t57) + m(4) * (t161 * t238 + t162 * t237 + t171 * t220 + t172 * t219 + t267 * t307 + t295 * t296) + m(3) * (t294 * t313 - t295 * t312 - t302 * t307 + t305 * t306) - t264 * t416 / 0.2e1 + t108 * (-mrSges(5,2) * t416 - mrSges(5,3) * t165) + t220 * (-mrSges(4,2) * t416 + mrSges(4,3) * t261) + t107 * (mrSges(5,1) * t416 - mrSges(5,3) * t164) + t219 * (mrSges(4,1) * t416 - mrSges(4,3) * t260) + ((t215 + t158) * t352 + t356 * t265 + t337 * (-Ifges(3,6) * t352 + t470)) * t440 / 0.2e1 + m(6) * (t103 * t37 + t119 * t28 + t12 * t49 + t13 * t48 + t537 * t6 + t66 * t7) + t537 * t46 + t7 * (-mrSges(6,1) * t377 - mrSges(6,3) * t225) + t218 * (-mrSges(5,1) * t377 + mrSges(5,2) * t245) + (t103 * t109 - t165 * t49 + t225 * t28 + t377 * t6) * mrSges(6,2) + (-t109 * t53 - t11 * t225 + t165 * t39 - t2 * t377) * mrSges(7,3) + (Ifges(6,4) * t225 - Ifges(6,6) * t377) * t528 + (Ifges(7,5) * t225 - Ifges(7,6) * t377) * t527 - t142 * (Ifges(5,4) * t245 + Ifges(5,2) * t377 - Ifges(5,6) * t449) / 0.2e1 + t141 * (Ifges(5,1) * t245 + Ifges(5,4) * t377 - Ifges(5,5) * t449) / 0.2e1 + t30 * (mrSges(5,2) * t449 + mrSges(5,3) * t377) + t4 * (mrSges(7,1) * t377 + mrSges(7,2) * t225) + (t225 * t548 - t377 * t547) * t529 - t543 * t377 / 0.2e1 + (t225 * t547 - t377 * t545) * t519 + t295 * (-mrSges(4,1) * t309 + mrSges(4,2) * t310) + t306 * t301 + t296 * t188 + t267 * (-mrSges(4,1) * t261 + mrSges(4,2) * t260) + t260 * t217 / 0.2e1 + t261 * t216 / 0.2e1 + t171 * t253 + t172 * t254 + t250 * t69 + t232 * (mrSges(5,1) * t165 + mrSges(5,2) * t164) + t235 * t166 + t237 * t233 + t238 * t234 + (t294 * t356 + t295 * t352 + (-t302 * t356 - t305 * t352) * qJD(2)) * t347 * mrSges(3,3) + t40 * t207 + t41 * t208 - t165 * t159 / 0.2e1 + t38 * (-mrSges(7,1) * t165 + mrSges(7,2) * t109) + t48 * (mrSges(6,1) * t165 - mrSges(6,3) * t109) + t56 * t47 + t57 * t45 + t66 * t44 + (Ifges(4,5) * t260 + Ifges(4,6) * t261 + Ifges(4,3) * t416) * t494 + t175 * t495 + t174 * t496 + (Ifges(4,1) * t260 + Ifges(4,4) * t261 + Ifges(4,5) * t416) * t497 + (Ifges(4,4) * t260 + Ifges(4,2) * t261 + Ifges(4,6) * t416) * t499 + (Ifges(4,4) * t310 + Ifges(4,2) * t309 - Ifges(4,6) * t449) * t502 + (Ifges(4,1) * t310 + Ifges(4,4) * t309 - Ifges(4,5) * t449) * t503 + t65 * t504 + t64 * t506 + t109 * t570 + (Ifges(5,3) * t416 + t394) * t550 + t71 * t32 + t119 * t33 + t128 * t124 + t129 * t125 + t14 * t131 + t37 * t132 + t8 * t148 + t12 * t149 + t13 * t150 + t9 * t151 + t164 * t160 / 0.2e1; (-t219 * t239 - t220 * t240 - t267 * t305 - pkin(2) * t295 + (-qJD(3) * t380 + t382) * pkin(9)) * m(4) + (t468 / 0.2e1 + t462 / 0.2e1 + t351 * pkin(3) * t166 + (-t253 * t351 - t254 * t355) * pkin(9) + t533) * qJD(3) - t410 * (-Ifges(5,4) * t270 - Ifges(5,2) * t269) / 0.2e1 - t378 * (-Ifges(5,1) * t270 - Ifges(5,4) * t269) / 0.2e1 + (-t107 * t270 + t108 * t269) * mrSges(5,3) - t232 * (mrSges(5,1) * t269 - mrSges(5,2) * t270) - qJD(4) * (-t462 - t468) / 0.2e1 - (t33 - t124) * t535 + (t242 * t53 - t269 * t39) * mrSges(7,3) + t562 * t132 + t563 * t208 + t564 * t207 - t538 * t242 / 0.2e1 + (t563 * t107 + t564 * t108 + t218 * t345 + t576 * t232 + t272 * t30 + t535 * t31) * m(5) + (Ifges(6,1) * t242 + Ifges(6,5) * t269) * t513 + (-t103 * t242 + t269 * t49) * mrSges(6,2) + (Ifges(7,4) * t242 + Ifges(7,2) * t269) * t509 + (Ifges(7,1) * t242 + Ifges(7,4) * t269) * t513 + (-t233 * t351 + t234 * t355) * pkin(9) - t455 + (Ifges(6,5) * t242 + Ifges(6,3) * t269) * t509 + (Ifges(5,6) * t551 - t464 / 0.2e1 - t531) * t258 + (-t30 * mrSges(5,3) - t460 / 0.2e1 + t20 / 0.2e1 + t21 / 0.2e1 - t64 / 0.2e1 + t218 * mrSges(5,1) - t427 * t88 + t428 * t87 + (t530 - t426) * t142 + t364) * t316 + (Ifges(6,4) * t242 + Ifges(6,6) * t269) * t514 + (Ifges(7,5) * t242 + Ifges(7,6) * t269) * t515 + (-mrSges(4,1) * t355 + mrSges(4,2) * t351 - mrSges(3,1)) * t295 + (-Ifges(6,2) * t514 + Ifges(7,3) * t515 - t509 * t572 + t513 * t546 - t571) * t241 + (t391 * t527 + t65 / 0.2e1 + t11 * t403 + t398 * t528 + t28 * t405 + t22 * t493 + t19 * t492 + t461 / 0.2e1 - t459 / 0.2e1 - t31 * mrSges(5,3) + t218 * mrSges(5,2) + (-t349 * t6 - t353 * t7) * mrSges(6,3) + (-t2 * t349 + t353 * t4) * mrSges(7,2) + (t397 * t514 + t100 * t491 + t390 * t515 + t53 * t404 + t103 * t406 + (t349 * t48 - t353 * t49) * mrSges(6,3) + (-t349 * t38 - t353 * t39) * mrSges(7,2) + t557 * t513 + t558 * t509 + t538 * t493) * qJD(5) + (t400 + t402) * t529 + (t393 + t396) * t519 + (qJD(5) * t97 + t560) * t490) * t317 - t443 * t305 - (Ifges(5,5) * t550 + t478 / 0.2e1 - t554) * t257 + t536 * t46 + ((-t265 / 0.2e1 - t333 / 0.2e1 - t453 / 0.2e1 + t302 * mrSges(3,3) + t442 * t522 + (t501 + t258 / 0.2e1) * Ifges(5,6) + (t257 / 0.2e1 - t270 / 0.2e1) * Ifges(5,5) - t533) * t356 + (-t219 * mrSges(4,1) - t458 / 0.2e1 - t456 / 0.2e1 - t454 / 0.2e1 + t220 * mrSges(4,2) - t215 / 0.2e1 - t158 / 0.2e1 + t264 / 0.2e1 + t305 * mrSges(3,3) - t469 / 0.2e1 - t463 / 0.2e1 + t108 * mrSges(5,2) - t107 * mrSges(5,1) + Ifges(5,3) * t551 + (t523 + t477 / 0.2e1 + (Ifges(5,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t356) * t442 + (t337 / 0.2e1 - qJD(2)) * Ifges(3,6) + (Ifges(4,5) * t351 + Ifges(5,5) * t317 + Ifges(4,6) * t355 - Ifges(5,6) * t316) * qJD(2) / 0.2e1) * t352) * t442 - t302 * t301 + t272 * t125 - t38 * (-mrSges(7,1) * t269 + mrSges(7,2) * t242) - t48 * (mrSges(6,1) * t269 - mrSges(6,3) * t242) - t259 * t166 - t240 * t253 - t239 * t254 + t212 * t32 + t382 * mrSges(4,3) + t202 * t44 - pkin(2) * t188 + t184 * t47 + t185 * t45 + t345 * t69 + t565 * t131 + t329 + t351 * t175 / 0.2e1 + t566 * t149 + t567 * t150 + (t103 * t562 + t202 * t7 - t535 * t28 + t48 * t567 + t49 * t566 + t536 * t6) * m(6) + t568 * t151 + t569 * t148 + (t11 * t212 + t184 * t2 + t185 * t4 + t38 * t568 + t39 * t569 + t53 * t565) * m(7) + t174 * t489 + t99 * t501 + t98 * t501 + (Ifges(4,2) * t355 + t476) * t502 + (Ifges(4,1) * t351 + t475) * t503 - t270 * t517 + t269 * t518; (t361 + t464 / 0.2e1) * t378 + (-t478 / 0.2e1 + t534) * t410 + t359 + (t219 * t285 + t220 * t286) * mrSges(4,3) + (t284 - t63) * t131 + t360 * (pkin(3) * t350 + pkin(11)) + m(7) * (t11 * t314 + t284 * t53) - m(5) * (-t107 * t112 + t108 * t113) - t559 * (-pkin(4) - t487) + t314 * t32 + t444 * t112 + (t354 * t124 + t350 * t125 - t286 * t166 + (-t444 * t350 + (t349 * t445 + t353 * t446 + t207) * t354 + m(7) * (t38 * t448 + t39 * t447) + m(6) * (t103 * t350 + t447 * t49 - t448 * t48)) * qJD(4) + (0.2e1 * t232 * t498 + t30 * t350 + t31 * t354 + (-t107 * t350 + t108 * t354) * qJD(4)) * m(5)) * pkin(3) - m(6) * (t103 * t112 + t48 * t58 + t49 * t59) - m(7) * (t38 * t43 + t39 * t42 + t53 * t63) - (-Ifges(4,2) * t286 + t217 + t280) * t285 / 0.2e1 - t267 * (mrSges(4,1) * t286 + mrSges(4,2) * t285) - t219 * t253 + t220 * t254 - t113 * t207 + t421 - t330 * (Ifges(4,5) * t285 - Ifges(4,6) * t286) / 0.2e1 + t216 * t497 + (Ifges(4,1) * t285 - t457) * t498 - t42 * t148 - t59 * t149 - t58 * t150 - t43 * t151 - t161 * mrSges(4,2) + t162 * mrSges(4,1); ((-Ifges(5,1) / 0.2e1 + t530) * t378 + t534) * t410 + t359 + t361 * t378 + t542 * t131 - m(6) * (t103 * t108 + t48 * t60 + t49 * t61) - t107 * t207 + t360 * pkin(11) + t323 * t32 + t444 * t108 - t51 * t148 - t61 * t149 - t60 * t150 - t52 * t151 + t559 * pkin(4) + (t11 * t323 - t38 * t52 - t39 * t51 + t53 * t542) * m(7); (-t445 + t479) * t49 + (-t446 - t480) * t48 + (t197 * t38 + t198 * t39) * mrSges(7,2) - t103 * (mrSges(6,1) * t198 - mrSges(6,2) * t197) + t100 * t512 - t53 * (mrSges(7,1) * t198 + mrSges(7,3) * t197) + (Ifges(7,3) * t198 - t467) * t515 + t364 - pkin(5) * t45 + qJ(6) * t47 - t130 * t131 + qJD(6) * t148 + (-t197 * t547 - t198 * t572) * t509 + (-pkin(5) * t4 + qJ(6) * t2 - t130 * t53 - t38 * t49 + t39 * t539) * m(7) + (-Ifges(6,2) * t198 - t195 + t538) * t514 + (-t197 * t548 + t194 - t473 + t97) * t513 + t543; t198 * t131 - t227 * t148 + 0.2e1 * (t4 / 0.2e1 + t53 * t512 + t39 * t509) * m(7) + t45;];
tauc  = t1(:);
