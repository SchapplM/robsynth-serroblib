% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:44:27
% EndTime: 2018-11-23 17:44:55
% DurationCPUTime: 28.87s
% Computational Cost: add. (18291->822), mult. (48565->1130), div. (0->0), fcn. (38194->10), ass. (0->343)
t318 = cos(qJ(2));
t315 = sin(qJ(2));
t310 = sin(pkin(6));
t387 = qJD(1) * t310;
t369 = t315 * t387;
t312 = cos(pkin(6));
t386 = qJD(1) * t312;
t378 = pkin(1) * t386;
t269 = -pkin(8) * t369 + t318 * t378;
t327 = (pkin(2) * t315 - pkin(9) * t318) * t310;
t270 = qJD(1) * t327;
t314 = sin(qJ(3));
t317 = cos(qJ(3));
t202 = -t314 * t269 + t317 * t270;
t425 = -qJ(4) - pkin(9);
t360 = qJD(3) * t425;
t522 = -(-qJ(4) * t317 * t318 + pkin(3) * t315) * t387 - t202 - qJD(4) * t314 + t317 * t360;
t203 = t317 * t269 + t314 * t270;
t368 = t318 * t387;
t356 = t314 * t368;
t521 = -qJ(4) * t356 - qJD(4) * t317 - t314 * t360 + t203;
t309 = sin(pkin(11));
t311 = cos(pkin(11));
t286 = t309 * t317 + t311 * t314;
t229 = t286 * t368;
t277 = t286 * qJD(3);
t520 = -t229 + t277;
t486 = Ifges(7,4) + Ifges(6,4);
t493 = t309 * t522 - t521 * t311;
t272 = pkin(8) * t368 + t315 * t378;
t225 = pkin(3) * t356 + t272;
t383 = qJD(3) * t314;
t519 = pkin(3) * t383 - t225;
t484 = Ifges(7,2) + Ifges(6,2);
t483 = Ifges(7,6) + Ifges(6,6);
t518 = pkin(10) * t369 - t493;
t285 = t309 * t314 - t311 * t317;
t230 = t285 * t368;
t278 = t285 * qJD(3);
t517 = t519 + (-t230 + t278) * pkin(10) + t520 * pkin(4);
t293 = qJD(3) - t368;
t313 = sin(qJ(5));
t316 = cos(qJ(5));
t300 = qJD(2) + t386;
t251 = t300 * t317 - t314 * t369;
t252 = t300 * t314 + t317 * t369;
t329 = t251 * t309 + t311 * t252;
t164 = t293 * t316 - t313 * t329;
t358 = t311 * t251 - t252 * t309;
t186 = qJD(5) - t358;
t165 = t293 * t313 + t316 * t329;
t512 = t486 * t165;
t478 = t164 * t484 + t186 * t483 + t512;
t516 = -t478 / 0.2e1;
t237 = pkin(9) * t300 + t272;
t265 = (-pkin(2) * t318 - pkin(9) * t315 - pkin(1)) * t310;
t244 = qJD(1) * t265;
t182 = t237 * t317 + t244 * t314;
t271 = qJD(2) * t327;
t259 = qJD(1) * t271;
t395 = t310 * t315;
t301 = pkin(8) * t395;
t430 = pkin(1) * t318;
t281 = t312 * t430 - t301;
t273 = t281 * qJD(2);
t260 = qJD(1) * t273;
t129 = -qJD(3) * t182 + t317 * t259 - t260 * t314;
t384 = qJD(2) * t318;
t366 = t317 * t384;
t382 = qJD(3) * t317;
t218 = t300 * t382 + (-t315 * t383 + t366) * t387;
t385 = qJD(2) * t310;
t362 = qJD(1) * t385;
t355 = t315 * t362;
t72 = pkin(3) * t355 - qJ(4) * t218 - qJD(4) * t252 + t129;
t128 = -t237 * t383 + t244 * t382 + t314 * t259 + t317 * t260;
t367 = t314 * t384;
t219 = -t300 * t383 + (-t315 * t382 - t367) * t387;
t81 = qJ(4) * t219 + qJD(4) * t251 + t128;
t27 = t309 * t72 + t311 * t81;
t515 = t27 * mrSges(5,3);
t487 = Ifges(7,1) + Ifges(6,1);
t485 = Ifges(7,5) + Ifges(6,5);
t394 = t310 * t318;
t282 = t312 * t315 * pkin(1) + pkin(8) * t394;
t274 = t282 * qJD(2);
t261 = qJD(1) * t274;
t178 = -t219 * pkin(3) + t261;
t514 = t178 * mrSges(5,1);
t495 = t521 * t309 + t311 * t522;
t207 = t230 * t313 + t316 * t369;
t380 = qJD(5) * t316;
t513 = -t286 * t380 - t207;
t181 = -t237 * t314 + t317 * t244;
t153 = -qJ(4) * t252 + t181;
t135 = pkin(3) * t293 + t153;
t154 = qJ(4) * t251 + t182;
t393 = t311 * t154;
t71 = t309 * t135 + t393;
t65 = pkin(10) * t293 + t71;
t236 = -t300 * pkin(2) - t269;
t187 = -t251 * pkin(3) + qJD(4) + t236;
t94 = -pkin(4) * t358 - pkin(10) * t329 + t187;
t28 = -t313 * t65 + t316 * t94;
t12 = -qJ(6) * t165 + t28;
t10 = pkin(5) * t186 + t12;
t400 = t293 * Ifges(5,6);
t407 = t329 * Ifges(5,4);
t410 = t358 * Ifges(5,2);
t126 = t400 + t407 + t410;
t29 = t313 * t94 + t316 * t65;
t13 = qJ(6) * t164 + t29;
t511 = t187 * mrSges(5,1) + t28 * mrSges(6,1) + t10 * mrSges(7,1) - t29 * mrSges(6,2) - t13 * mrSges(7,2) - t71 * mrSges(5,3) - t126 / 0.2e1;
t436 = t293 / 0.2e1;
t510 = t329 / 0.2e1;
t509 = t358 / 0.2e1;
t508 = Ifges(4,3) + Ifges(5,3);
t507 = t313 * t518 + t316 * t517;
t308 = -pkin(3) * t317 - pkin(2);
t220 = pkin(4) * t285 - pkin(10) * t286 + t308;
t506 = t220 * t380 + t313 * t517 - t316 * t518;
t505 = t486 * t164;
t494 = pkin(4) * t369 - t495;
t333 = t28 * t316 + t29 * t313;
t348 = mrSges(7,1) * t313 + mrSges(7,2) * t316;
t350 = mrSges(6,1) * t313 + mrSges(6,2) * t316;
t432 = t316 / 0.2e1;
t435 = -t313 / 0.2e1;
t144 = t309 * t154;
t70 = t135 * t311 - t144;
t64 = -pkin(4) * t293 - t70;
t44 = -pkin(5) * t164 + qJD(6) + t64;
t451 = t186 / 0.2e1;
t456 = t165 / 0.2e1;
t458 = t164 / 0.2e1;
t502 = t486 * t313;
t473 = t316 * t487 - t502;
t503 = t486 * t316;
t474 = -t313 * t484 + t503;
t475 = -t313 * t483 + t316 * t485;
t477 = t165 * t487 + t186 * t485 + t505;
t504 = t348 * t44 + t350 * t64 + t432 * t477 + t435 * t478 + t451 * t475 + t456 * t473 + t458 * t474 - t333 * mrSges(6,3) - (t10 * t316 + t13 * t313) * mrSges(7,3);
t482 = Ifges(7,3) + Ifges(6,3);
t208 = -t230 * t316 + t313 * t369;
t294 = t425 * t314;
t295 = t425 * t317;
t235 = t294 * t309 - t295 * t311;
t224 = t316 * t235;
t328 = qJ(6) * t278 - qJD(6) * t286;
t501 = qJ(6) * t208 + t328 * t316 + (-t224 + (qJ(6) * t286 - t220) * t313) * qJD(5) + t507 + t520 * pkin(5);
t500 = (-qJD(5) * t235 + t328) * t313 + t506 + t513 * qJ(6);
t58 = t165 * Ifges(7,5) + t164 * Ifges(7,6) + t186 * Ifges(7,3);
t59 = t165 * Ifges(6,5) + t164 * Ifges(6,6) + t186 * Ifges(6,3);
t499 = t59 + t58;
t381 = qJD(5) * t313;
t498 = -t235 * t381 + t506;
t163 = t313 * t220 + t224;
t497 = -qJD(5) * t163 + t507;
t496 = t494 + (-t278 * t313 - t513) * pkin(5);
t461 = Ifges(5,1) * t510 + Ifges(5,4) * t509 + Ifges(5,5) * t436;
t492 = t187 * mrSges(5,2) - t70 * mrSges(5,3) + 0.2e1 * t461 + t504;
t491 = t313 * t485 + t316 * t483;
t490 = t316 * t484 + t502;
t489 = t313 * t487 + t503;
t488 = qJD(5) * t286;
t155 = t218 * t309 - t311 * t219;
t156 = t218 * t311 + t219 * t309;
t87 = qJD(5) * t164 + t156 * t316 + t313 * t355;
t88 = -qJD(5) * t165 - t156 * t313 + t316 * t355;
t14 = Ifges(7,5) * t87 + Ifges(7,6) * t88 + Ifges(7,3) * t155;
t15 = Ifges(6,5) * t87 + Ifges(6,6) * t88 + Ifges(6,3) * t155;
t481 = t15 + t14;
t480 = t155 * t483 + t484 * t88 + t486 * t87;
t479 = t155 * t485 + t486 * t88 + t487 * t87;
t279 = t312 * t317 - t314 * t395;
t280 = t312 * t314 + t317 * t395;
t210 = -t311 * t279 + t280 * t309;
t211 = t279 * t309 + t280 * t311;
t263 = t301 + (-pkin(2) - t430) * t312;
t217 = -t279 * pkin(3) + t263;
t124 = t210 * pkin(4) - t211 * pkin(10) + t217;
t264 = pkin(9) * t312 + t282;
t196 = -t314 * t264 + t317 * t265;
t159 = -pkin(3) * t394 - t280 * qJ(4) + t196;
t197 = t317 * t264 + t314 * t265;
t173 = qJ(4) * t279 + t197;
t103 = t309 * t159 + t311 * t173;
t98 = -pkin(10) * t394 + t103;
t42 = t313 * t124 + t316 * t98;
t472 = Ifges(4,5) * t218 + Ifges(5,5) * t156 + Ifges(4,6) * t219 - Ifges(5,6) * t155 + t355 * t508;
t24 = pkin(10) * t355 + t27;
t55 = t155 * pkin(4) - t156 * pkin(10) + t178;
t5 = t316 * t24 + t313 * t55 + t94 * t380 - t381 * t65;
t6 = -qJD(5) * t29 - t24 * t313 + t316 * t55;
t471 = -t313 * t6 + t316 * t5;
t404 = t252 * Ifges(4,4);
t176 = t251 * Ifges(4,2) + t293 * Ifges(4,6) + t404;
t245 = Ifges(4,4) * t251;
t177 = t252 * Ifges(4,1) + t293 * Ifges(4,5) + t245;
t330 = t181 * t317 + t182 * t314;
t422 = Ifges(4,4) * t317;
t423 = Ifges(4,4) * t314;
t431 = t317 / 0.2e1;
t440 = t252 / 0.2e1;
t441 = t251 / 0.2e1;
t470 = -t330 * mrSges(4,3) + t236 * (mrSges(4,1) * t314 + mrSges(4,2) * t317) + (-Ifges(4,2) * t314 + t422) * t441 + (Ifges(4,1) * t317 - t423) * t440 + (Ifges(4,5) * t317 - Ifges(4,6) * t314) * t436 - t314 * t176 / 0.2e1 + t177 * t431;
t373 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t375 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t376 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t469 = t164 * t375 + t165 * t376 + t186 * t373 + t58 / 0.2e1 + t59 / 0.2e1 - t410 / 0.2e1 - t407 / 0.2e1 - t400 / 0.2e1 + t511;
t468 = t87 / 0.2e1;
t467 = t88 / 0.2e1;
t466 = pkin(1) * mrSges(3,1);
t465 = pkin(1) * mrSges(3,2);
t460 = t155 / 0.2e1;
t459 = -t164 / 0.2e1;
t457 = -t165 / 0.2e1;
t452 = -t186 / 0.2e1;
t448 = -t210 / 0.2e1;
t446 = t211 / 0.2e1;
t445 = t218 / 0.2e1;
t444 = t219 / 0.2e1;
t439 = t279 / 0.2e1;
t438 = t280 / 0.2e1;
t437 = -t293 / 0.2e1;
t429 = pkin(3) * t252;
t424 = Ifges(3,4) * t315;
t417 = Ifges(3,5) * t318;
t416 = t155 * Ifges(5,4);
t415 = t156 * Ifges(5,1);
t414 = t156 * Ifges(5,4);
t409 = t358 * Ifges(5,6);
t406 = t329 * Ifges(5,5);
t405 = t251 * Ifges(4,6);
t403 = t252 * Ifges(4,5);
t402 = t260 * mrSges(3,2);
t399 = t300 * Ifges(3,5);
t136 = -t264 * t383 + t265 * t382 + t314 * t271 + t317 * t273;
t227 = -qJD(3) * t280 - t310 * t367;
t104 = qJ(4) * t227 + qJD(4) * t279 + t136;
t137 = -qJD(3) * t197 + t317 * t271 - t273 * t314;
t226 = qJD(3) * t279 + t310 * t366;
t364 = t315 * t385;
t99 = pkin(3) * t364 - qJ(4) * t226 - qJD(4) * t280 + t137;
t40 = t311 * t104 + t309 * t99;
t114 = pkin(4) * t329 - pkin(10) * t358 + t429;
t80 = t153 * t311 - t144;
t35 = t313 * t114 + t316 * t80;
t398 = qJ(6) * t316;
t397 = t358 * t313;
t396 = t286 * t313;
t306 = pkin(3) * t309 + pkin(10);
t392 = qJ(6) + t306;
t106 = -mrSges(6,1) * t164 + mrSges(6,2) * t165;
t172 = mrSges(5,1) * t293 - mrSges(5,3) * t329;
t389 = t172 - t106;
t388 = -mrSges(3,1) * t300 - mrSges(4,1) * t251 + mrSges(4,2) * t252 + mrSges(3,3) * t369;
t374 = -Ifges(5,3) / 0.2e1 - Ifges(4,3) / 0.2e1;
t307 = -pkin(3) * t311 - pkin(4);
t30 = -t88 * mrSges(7,1) + t87 * mrSges(7,2);
t26 = -t309 * t81 + t311 * t72;
t89 = t155 * mrSges(5,1) + t156 * mrSges(5,2);
t39 = -t309 * t104 + t311 * t99;
t34 = t316 * t114 - t313 * t80;
t41 = t316 * t124 - t313 * t98;
t79 = t153 * t309 + t393;
t102 = t159 * t311 - t309 * t173;
t162 = t316 * t220 - t235 * t313;
t234 = -t311 * t294 - t295 * t309;
t357 = qJD(5) * t392;
t1 = pkin(5) * t155 - qJ(6) * t87 - qJD(6) * t165 + t6;
t2 = qJ(6) * t88 + qJD(6) * t164 + t5;
t353 = -t1 * t316 - t2 * t313;
t352 = -t313 * t5 - t316 * t6;
t97 = pkin(4) * t394 - t102;
t351 = mrSges(6,1) * t316 - mrSges(6,2) * t313;
t349 = mrSges(7,1) * t316 - mrSges(7,2) * t313;
t334 = t10 * t313 - t13 * t316;
t332 = t28 * t313 - t29 * t316;
t331 = t128 * t317 - t129 * t314;
t184 = -t313 * t211 - t316 * t394;
t326 = -t316 * t211 + t313 * t394;
t37 = pkin(10) * t364 + t40;
t169 = t226 * t309 - t311 * t227;
t170 = t226 * t311 + t227 * t309;
t194 = -t227 * pkin(3) + t274;
t68 = t169 * pkin(4) - t170 * pkin(10) + t194;
t7 = t124 * t380 + t313 * t68 + t316 * t37 - t381 * t98;
t36 = -pkin(4) * t364 - t39;
t323 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t23 = -pkin(4) * t355 - t26;
t8 = -qJD(5) * t42 - t313 * t37 + t316 * t68;
t296 = Ifges(3,4) * t368;
t292 = t362 * t417;
t289 = -pkin(5) * t316 + t307;
t284 = t392 * t316;
t283 = t392 * t313;
t268 = -t300 * mrSges(3,2) + mrSges(3,3) * t368;
t250 = -qJD(6) * t313 - t316 * t357;
t249 = qJD(6) * t316 - t313 * t357;
t232 = Ifges(3,1) * t369 + t296 + t399;
t231 = Ifges(3,6) * t300 + (Ifges(3,2) * t318 + t424) * t387;
t222 = mrSges(4,1) * t293 - mrSges(4,3) * t252;
t221 = -mrSges(4,2) * t293 + mrSges(4,3) * t251;
t201 = pkin(5) * t396 + t234;
t193 = -mrSges(4,2) * t355 + mrSges(4,3) * t219;
t192 = mrSges(4,1) * t355 - mrSges(4,3) * t218;
t175 = t293 * Ifges(4,3) + t403 + t405;
t171 = -mrSges(5,2) * t293 + mrSges(5,3) * t358;
t157 = -mrSges(4,1) * t219 + mrSges(4,2) * t218;
t143 = t218 * Ifges(4,1) + t219 * Ifges(4,4) + Ifges(4,5) * t355;
t142 = t218 * Ifges(4,4) + t219 * Ifges(4,2) + Ifges(4,6) * t355;
t140 = -qJ(6) * t396 + t163;
t134 = mrSges(5,1) * t355 - mrSges(5,3) * t156;
t133 = -mrSges(5,2) * t355 - mrSges(5,3) * t155;
t131 = pkin(5) * t285 - t286 * t398 + t162;
t130 = -mrSges(5,1) * t358 + mrSges(5,2) * t329;
t125 = t293 * Ifges(5,3) + t406 + t409;
t120 = mrSges(6,1) * t186 - mrSges(6,3) * t165;
t119 = mrSges(7,1) * t186 - mrSges(7,3) * t165;
t118 = -mrSges(6,2) * t186 + mrSges(6,3) * t164;
t117 = -mrSges(7,2) * t186 + mrSges(7,3) * t164;
t109 = qJD(5) * t326 - t313 * t170 + t316 * t364;
t108 = qJD(5) * t184 + t316 * t170 + t313 * t364;
t105 = -mrSges(7,1) * t164 + mrSges(7,2) * t165;
t75 = Ifges(5,5) * t355 + t415 - t416;
t74 = -t155 * Ifges(5,2) + Ifges(5,6) * t355 + t414;
t57 = -pkin(5) * t184 + t97;
t52 = pkin(5) * t397 + t79;
t48 = -mrSges(6,2) * t155 + mrSges(6,3) * t88;
t47 = -mrSges(7,2) * t155 + mrSges(7,3) * t88;
t46 = mrSges(6,1) * t155 - mrSges(6,3) * t87;
t45 = mrSges(7,1) * t155 - mrSges(7,3) * t87;
t32 = qJ(6) * t184 + t42;
t31 = -mrSges(6,1) * t88 + mrSges(6,2) * t87;
t25 = -qJ(6) * t397 + t35;
t21 = pkin(5) * t210 + qJ(6) * t326 + t41;
t20 = pkin(5) * t329 - t358 * t398 + t34;
t11 = -pkin(5) * t109 + t36;
t9 = -pkin(5) * t88 + t23;
t4 = qJ(6) * t109 + qJD(6) * t184 + t7;
t3 = pkin(5) * t169 - qJ(6) * t108 + qJD(6) * t326 + t8;
t16 = [(t260 * t318 + t261 * t315 + (-t269 * t318 - t272 * t315) * qJD(2)) * t310 * mrSges(3,3) - t210 * t515 + t143 * t438 + t142 * t439 + (Ifges(4,1) * t226 + Ifges(4,4) * t227 + Ifges(4,5) * t364) * t440 + (Ifges(4,4) * t226 + Ifges(4,2) * t227 + Ifges(4,6) * t364) * t441 + (Ifges(4,4) * t280 + Ifges(4,2) * t279 - Ifges(4,6) * t394) * t444 + (Ifges(4,1) * t280 + Ifges(4,4) * t279 - Ifges(4,5) * t394) * t445 + t75 * t446 + t74 * t448 - t472 * t394 / 0.2e1 + t236 * (-mrSges(4,1) * t227 + mrSges(4,2) * t226) + t227 * t176 / 0.2e1 + t136 * t221 + t137 * t222 + t226 * t177 / 0.2e1 + ((t175 + t125) * t315 + t318 * t232 + t300 * (-Ifges(3,6) * t315 + t417)) * t385 / 0.2e1 - t479 * t326 / 0.2e1 + (t184 * t483 + t210 * t482 - t326 * t485) * t460 + (t184 * t484 + t210 * t483 - t326 * t486) * t467 + (t184 * t486 + t210 * t485 - t326 * t487) * t468 + t6 * (mrSges(6,1) * t210 + mrSges(6,3) * t326) + t1 * (mrSges(7,1) * t210 + mrSges(7,3) * t326) + t23 * (-mrSges(6,1) * t184 - mrSges(6,2) * t326) + t9 * (-mrSges(7,1) * t184 - mrSges(7,2) * t326) + (Ifges(5,1) * t170 + Ifges(5,5) * t364) * t510 + (t478 / 0.2e1 + t483 * t451 + t484 * t458 + t486 * t456 + t29 * mrSges(6,3) + t13 * mrSges(7,3) - t64 * mrSges(6,1) - t44 * mrSges(7,1)) * t109 + (Ifges(5,4) * t170 + Ifges(5,6) * t364) * t509 + t480 * t184 / 0.2e1 + t481 * t210 / 0.2e1 + (Ifges(4,5) * t226 + Ifges(5,5) * t170 + Ifges(4,6) * t227 + t364 * t508) * t436 + t170 * t461 + (t170 * t187 + t178 * t211 + t27 * t394 - t364 * t71) * mrSges(5,2) + (t477 / 0.2e1 + t485 * t451 + t486 * t458 + t487 * t456 - t28 * mrSges(6,3) - t10 * mrSges(7,3) + t64 * mrSges(6,2) + t44 * mrSges(7,2)) * t108 + (t482 * t451 + t483 * t458 + t485 * t456 + t499 / 0.2e1 - Ifges(5,6) * t436 - Ifges(5,2) * t509 - Ifges(5,4) * t510 + t511) * t169 + t5 * (-mrSges(6,2) * t210 + mrSges(6,3) * t184) + t2 * (-mrSges(7,2) * t210 + mrSges(7,3) * t184) + t196 * t192 + t197 * t193 + t194 * t130 + t40 * t171 + t39 * t172 + t103 * t133 + t102 * t134 + t4 * t117 + t7 * t118 + t3 * t119 + t8 * t120 + t11 * t105 + t36 * t106 + t97 * t31 + t57 * t30 + t21 * t45 + t41 * t46 + t32 * t47 + t42 * t48 + ((-t281 * mrSges(3,3) + Ifges(3,5) * t312 / 0.2e1 + (-0.2e1 * t465 + 0.3e1 / 0.2e1 * Ifges(3,4) * t318) * t310) * t318 + (-t282 * mrSges(3,3) + Ifges(4,5) * t438 + Ifges(4,6) * t439 + Ifges(5,5) * t446 + Ifges(5,6) * t448 - Ifges(3,6) * t312 + (-0.2e1 * t466 - 0.3e1 / 0.2e1 * t424) * t310 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + t374) * t394) * t315) * t362 + (t292 / 0.2e1 - t261 * mrSges(3,1) - t402) * t312 - t155 * (Ifges(5,4) * t211 - Ifges(5,2) * t210 - Ifges(5,6) * t394) / 0.2e1 + t156 * (Ifges(5,1) * t211 - Ifges(5,4) * t210 - Ifges(5,5) * t394) / 0.2e1 + t26 * (-mrSges(5,1) * t394 - t211 * mrSges(5,3)) + t129 * (-mrSges(4,1) * t394 - t280 * mrSges(4,3)) + t128 * (mrSges(4,2) * t394 + t279 * mrSges(4,3)) + t388 * t274 + m(3) * (t260 * t282 - t261 * t281 - t269 * t274 + t272 * t273) + m(4) * (t128 * t197 + t129 * t196 + t136 * t182 + t137 * t181 + t236 * t274 + t261 * t263) + m(5) * (t102 * t26 + t103 * t27 + t178 * t217 + t187 * t194 + t39 * t70 + t40 * t71) + m(6) * (t23 * t97 + t28 * t8 + t29 * t7 + t36 * t64 + t41 * t6 + t42 * t5) + m(7) * (t1 * t21 + t10 * t3 + t11 * t44 + t13 * t4 + t2 * t32 + t57 * t9) - t231 * t364 / 0.2e1 + t182 * (-mrSges(4,2) * t364 + mrSges(4,3) * t227) + t70 * (mrSges(5,1) * t364 - mrSges(5,3) * t170) + t181 * (mrSges(4,1) * t364 - mrSges(4,3) * t226) + t263 * t157 + t273 * t268 + t261 * (-mrSges(4,1) * t279 + mrSges(4,2) * t280) + t210 * t514 + t217 * t89; t142 * t431 + (-mrSges(4,1) * t317 + mrSges(4,2) * t314 - mrSges(3,1)) * t261 + (t14 / 0.2e1 + t15 / 0.2e1 - t515 - t414 / 0.2e1 - t74 / 0.2e1 + t514 + t375 * t88 + t376 * t87 + (Ifges(5,2) / 0.2e1 + t373) * t155 + t323) * t285 + (Ifges(4,2) * t317 + t423) * t444 + (Ifges(4,1) * t314 + t422) * t445 + t235 * t133 + t229 * t126 / 0.2e1 + (t314 * pkin(3) * t130 + (-t221 * t314 - t222 * t317) * pkin(9) + t470) * qJD(3) + ((t269 * mrSges(3,3) + t387 * t465 - t399 / 0.2e1 - t232 / 0.2e1 - t296 / 0.2e1 - t470) * t318 + (t272 * mrSges(3,3) - t406 / 0.2e1 - t409 / 0.2e1 + t71 * mrSges(5,2) - t70 * mrSges(5,1) + t182 * mrSges(4,2) - t181 * mrSges(4,1) - t405 / 0.2e1 - t403 / 0.2e1 + t231 / 0.2e1 - t125 / 0.2e1 - t175 / 0.2e1 + t374 * t293 + (t466 + t424 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t318) * t387 + (-qJD(2) + t300 / 0.2e1) * Ifges(3,6) + (Ifges(4,5) * t314 + Ifges(5,5) * t286 + Ifges(4,6) * t317 - Ifges(5,6) * t285) * qJD(2) / 0.2e1) * t315) * t387 - t28 * (mrSges(6,1) * t229 - mrSges(6,3) * t208) - t10 * (mrSges(7,1) * t229 - mrSges(7,3) * t208) - t29 * (-mrSges(6,2) * t229 + mrSges(6,3) * t207) - t13 * (-mrSges(7,2) * t229 + mrSges(7,3) * t207) - t203 * t221 - t202 * t222 - t225 * t130 - t402 + t292 + (t75 / 0.2e1 - t416 / 0.2e1 + t415 / 0.2e1 - t26 * mrSges(5,3) + t178 * mrSges(5,2) + t23 * t350 + t9 * t348 + t353 * mrSges(7,3) + t352 * mrSges(6,3) + (mrSges(6,3) * t332 + mrSges(7,3) * t334 + t316 * t516 + t349 * t44 + t351 * t64) * qJD(5) + t473 * t468 + t474 * t467 + t475 * t460 + t479 * t432 + (qJD(5) * t477 + t480) * t435) * t286 + (-pkin(2) * t261 - t181 * t202 - t182 * t203 - t236 * t272 + (-qJD(3) * t330 + t331) * pkin(9)) * m(4) + (-t192 * t314 + t193 * t317) * pkin(9) + (t31 - t134) * t234 - t64 * (-mrSges(6,1) * t207 + mrSges(6,2) * t208) - t44 * (-mrSges(7,1) * t207 + mrSges(7,2) * t208) + t201 * t30 + t469 * t277 + t162 * t46 + t163 * t48 - pkin(2) * t157 + t140 * t47 + t131 * t45 + (t178 * t308 + t187 * t519 - t234 * t26 + t235 * t27 + t493 * t71 + t495 * t70) * m(5) - t388 * t272 + (t486 * t207 + t208 * t487 + t485 * t229 + t489 * t488) * t457 + (t207 * t483 + t208 * t485 + t229 * t482 + t488 * t491) * t452 + (t207 * t484 + t208 * t486 + t229 * t483 + t488 * t490) * t459 - t477 * t208 / 0.2e1 - t499 * t229 / 0.2e1 + t500 * t117 + (t1 * t131 + t10 * t501 + t13 * t500 + t140 * t2 + t201 * t9 + t44 * t496) * m(7) + t501 * t119 + t497 * t120 + t498 * t118 + (t162 * t6 + t163 * t5 + t23 * t234 + t28 * t497 + t29 * t498 + t494 * t64) * m(6) + t496 * t105 + (-Ifges(5,5) * t230 - Ifges(5,6) * t229) * t437 - t187 * (mrSges(5,1) * t229 - mrSges(5,2) * t230) - t358 * (-Ifges(5,4) * t230 - Ifges(5,2) * t229) / 0.2e1 - t329 * (-Ifges(5,1) * t230 - Ifges(5,4) * t229) / 0.2e1 + (t229 * t71 - t230 * t70) * mrSges(5,3) - t269 * t268 + t308 * t89 + t230 * t461 + t314 * t143 / 0.2e1 + t331 * mrSges(4,3) + t207 * t516 - t492 * t278 + t493 * t171 + t494 * t106 + t495 * t172; (Ifges(4,5) * t251 - Ifges(4,6) * t252) * t437 + t176 * t440 + (-t313 * t46 + t316 * t48 + m(6) * t471 + (-m(6) * t333 - t313 * t118 - t316 * t120) * qJD(5)) * t306 + t471 * mrSges(6,3) - t181 * t221 + t182 * t222 + (t181 * t251 + t182 * t252) * mrSges(4,3) + (t249 - t25) * t117 + (-t130 * t252 + t133 * t309 + t134 * t311) * pkin(3) - t236 * (mrSges(4,1) * t252 + mrSges(4,2) * t251) + (t250 - t20) * t119 + t479 * t313 / 0.2e1 + t480 * t432 + (t23 * t307 - t28 * t34 - t29 * t35 - t64 * t79) * m(6) + t472 - m(7) * (t10 * t20 + t13 * t25 + t44 * t52) + (-t1 * t313 + t2 * t316) * mrSges(7,3) - (-Ifges(4,2) * t252 + t177 + t245) * t251 / 0.2e1 - t469 * t329 - t80 * t171 - t128 * mrSges(4,2) + t129 * mrSges(4,1) - t35 * t118 - t34 * t120 - t52 * t105 + t26 * mrSges(5,1) - t27 * mrSges(5,2) - t252 * (Ifges(4,1) * t251 - t404) / 0.2e1 + t389 * t79 + ((m(7) * t44 + t105) * t313 * pkin(5) + t504) * qJD(5) + m(7) * (-t1 * t283 + t10 * t250 + t13 * t249 + t2 * t284 + t289 * t9) - t283 * t45 + t284 * t47 + t289 * t30 + t307 * t31 - t9 * t349 - t23 * t351 + (-t187 * t429 + t70 * t79 - t71 * t80 + (t26 * t311 + t27 * t309) * pkin(3)) * m(5) + t489 * t468 + t490 * t467 + t491 * t460 - t492 * t358; -t358 * t171 + (-t105 + t389) * t329 + (t45 + t46 + t186 * (t117 + t118)) * t316 + (t47 + t48 - t186 * (t119 + t120)) * t313 + t89 + (-t186 * t334 - t329 * t44 - t353) * m(7) + (-t186 * t332 - t329 * t64 - t352) * m(6) + (t329 * t70 - t358 * t71 + t178) * m(5); (-t105 * t165 + t45) * pkin(5) + t323 + t481 + (t10 * t164 + t13 * t165) * mrSges(7,3) + (t164 * t28 + t165 * t29) * mrSges(6,3) - t64 * (mrSges(6,1) * t165 + mrSges(6,2) * t164) - t44 * (mrSges(7,1) * t165 + mrSges(7,2) * t164) + (-(-t10 + t12) * t13 + (-t165 * t44 + t1) * pkin(5)) * m(7) - t28 * t118 + t13 * t119 + t29 * t120 - t12 * t117 + (t164 * t487 - t512) * t457 + t478 * t456 + (t164 * t485 - t165 * t483) * t452 + (-t165 * t484 + t477 + t505) * t459; -t164 * t117 + t165 * t119 + 0.2e1 * (t9 / 0.2e1 + t10 * t456 + t13 * t459) * m(7) + t30;];
tauc  = t16(:);
