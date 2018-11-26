% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2018-11-23 17:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:39:19
% EndTime: 2018-11-23 17:39:44
% DurationCPUTime: 25.54s
% Computational Cost: add. (13568->883), mult. (36232->1206), div. (0->0), fcn. (27765->10), ass. (0->367)
t318 = cos(qJ(2));
t405 = cos(pkin(6));
t374 = pkin(1) * t405;
t360 = t318 * t374;
t315 = sin(qJ(2));
t311 = sin(pkin(6));
t387 = qJD(1) * t311;
t372 = t315 * t387;
t267 = -pkin(8) * t372 + qJD(1) * t360;
t332 = (pkin(2) * t315 - pkin(9) * t318) * t311;
t268 = qJD(1) * t332;
t314 = sin(qJ(3));
t317 = cos(qJ(3));
t182 = t317 * t267 + t314 * t268;
t165 = qJ(4) * t372 + t182;
t306 = t315 * t374;
t342 = pkin(3) * t314 - qJ(4) * t317;
t399 = t311 * t318;
t187 = (t306 + (pkin(8) + t342) * t399) * qJD(1);
t310 = sin(pkin(11));
t312 = cos(pkin(11));
t106 = t312 * t165 + t310 * t187;
t371 = t318 * t387;
t359 = t314 * t371;
t384 = qJD(3) * t314;
t511 = -qJD(5) * t317 - t106 + (-t359 + t384) * qJ(5);
t105 = -t310 * t165 + t312 * t187;
t274 = qJD(3) * t342 - qJD(4) * t314;
t402 = t274 * t312;
t510 = t105 - t402;
t492 = Ifges(5,1) + Ifges(6,1);
t491 = Ifges(6,4) + Ifges(5,5);
t490 = Ifges(6,5) - Ifges(5,4);
t396 = t312 * t318;
t226 = (t310 * t315 + t317 * t396) * t387;
t375 = -pkin(9) * t310 - pkin(4);
t397 = t312 * t317;
t468 = pkin(4) + pkin(5);
t509 = t226 * pkin(10) + t359 * t468 + (-pkin(10) * t397 + (-pkin(5) + t375) * t314) * qJD(3) + t510;
t401 = t310 * t317;
t225 = -t312 * t372 + t371 * t401;
t257 = t310 * t274;
t398 = t312 * t314;
t508 = -pkin(10) * t225 + t257 + (-pkin(9) * t398 + pkin(10) * t401) * qJD(3) + t511;
t181 = -t314 * t267 + t317 * t268;
t331 = pkin(3) * t372 + t181;
t507 = qJ(5) * t226 - qJD(5) * t398 + t331;
t336 = qJD(1) * t405 + qJD(2);
t245 = t314 * t372 - t317 * t336;
t355 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1;
t298 = qJD(3) - t371;
t409 = t298 * Ifges(4,6);
t246 = t314 * t336 + t317 * t372;
t432 = Ifges(4,4) * t246;
t142 = -t245 * Ifges(4,2) + t409 + t432;
t388 = pkin(8) * t399 + t306;
t259 = t405 * pkin(9) + t388;
t221 = qJD(2) * pkin(9) + qJD(1) * t259;
t261 = (-pkin(2) * t318 - pkin(9) * t315 - pkin(1)) * t311;
t233 = qJD(1) * t261;
t146 = t221 * t317 + t233 * t314;
t189 = t246 * t310 - t312 * t298;
t220 = -pkin(2) * t336 - t267;
t338 = t312 * t246 + t298 * t310;
t377 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t378 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t238 = qJD(6) - t245;
t416 = t238 * Ifges(7,3);
t313 = sin(qJ(6));
t316 = cos(qJ(6));
t494 = t189 * t313 + t316 * t338;
t420 = t494 * Ifges(7,5);
t118 = t189 * t316 - t313 * t338;
t421 = t118 * Ifges(7,6);
t39 = t416 + t420 + t421;
t128 = t245 * pkin(3) - t246 * qJ(4) + t220;
t132 = qJ(4) * t298 + t146;
t57 = t128 * t312 - t310 * t132;
t351 = qJD(5) - t57;
t32 = -pkin(10) * t338 - t245 * t468 + t351;
t58 = t310 * t128 + t312 * t132;
t51 = t245 * qJ(5) + t58;
t34 = pkin(10) * t189 + t51;
t5 = -t313 * t34 + t316 * t32;
t50 = -pkin(4) * t245 + t351;
t6 = t313 * t32 + t316 * t34;
t98 = Ifges(5,5) * t338 - t189 * Ifges(5,6) + t245 * Ifges(5,3);
t99 = Ifges(6,4) * t338 + t245 * Ifges(6,2) + t189 * Ifges(6,6);
t475 = t377 * t189 + t378 * t338 - t146 * mrSges(4,3) - t5 * mrSges(7,1) - t50 * mrSges(6,1) - t58 * mrSges(5,2) - t142 / 0.2e1 - t39 / 0.2e1 + t98 / 0.2e1 + t99 / 0.2e1 - t432 / 0.2e1 - t421 / 0.2e1 - t420 / 0.2e1 + t220 * mrSges(4,1) - t416 / 0.2e1 - t409 / 0.2e1 + t51 * mrSges(6,3) + t57 * mrSges(5,1) + t6 * mrSges(7,2);
t506 = -t355 * t245 - t475;
t380 = qJD(1) * qJD(2);
t365 = t318 * t380;
t357 = t311 * t365;
t198 = -qJD(3) * t245 + t317 * t357;
t505 = t198 / 0.2e1;
t199 = qJD(3) * t246 + t314 * t357;
t454 = -t199 / 0.2e1;
t504 = Ifges(5,6) - Ifges(6,6);
t292 = -pkin(3) * t317 - qJ(4) * t314 - pkin(2);
t301 = pkin(9) * t401;
t309 = t317 * pkin(4);
t436 = pkin(10) * t314;
t197 = pkin(5) * t317 + t301 + t309 + (-t292 - t436) * t312;
t244 = pkin(9) * t397 + t310 * t292;
t222 = -qJ(5) * t317 + t244;
t203 = t310 * t436 + t222;
t124 = t197 * t316 - t203 * t313;
t503 = qJD(6) * t124 + t313 * t509 + t316 * t508;
t125 = t197 * t313 + t203 * t316;
t502 = -qJD(6) * t125 - t313 * t508 + t316 * t509;
t364 = t315 * t380;
t358 = t311 * t364;
t163 = t198 * t310 - t312 * t358;
t164 = t198 * t312 + t310 * t358;
t487 = t490 * t163 + t164 * t492 + t491 * t199;
t435 = -pkin(10) + qJ(4);
t293 = t435 * t310;
t294 = t435 * t312;
t214 = t293 * t316 - t294 * t313;
t337 = t310 * t313 + t312 * t316;
t145 = -t314 * t221 + t317 * t233;
t139 = t310 * t145;
t172 = pkin(3) * t246 + qJ(4) * t245;
t437 = pkin(10) * t245;
t47 = t139 + (-t172 + t437) * t312 - t468 * t246;
t86 = t312 * t145 + t310 * t172;
t68 = t246 * qJ(5) + t86;
t55 = -t310 * t437 + t68;
t501 = qJD(4) * t337 + qJD(6) * t214 - t313 * t47 - t316 * t55;
t215 = t293 * t313 + t294 * t316;
t285 = t310 * t316 - t312 * t313;
t500 = qJD(4) * t285 - qJD(6) * t215 + t313 * t55 - t316 * t47;
t379 = pkin(9) * t384;
t207 = -t312 * t379 + t257;
t499 = t207 + t511;
t403 = qJ(5) * t312;
t476 = t310 * t468 - t403;
t330 = -pkin(9) - t476;
t383 = qJD(3) * t317;
t498 = t225 * t468 + t330 * t383 - t507;
t497 = pkin(4) * t359 + t375 * t384 + t510;
t482 = t490 * t189 + t491 * t245 + t338 * t492;
t341 = -pkin(4) * t310 + t403;
t334 = pkin(9) - t341;
t496 = -pkin(4) * t225 + t334 * t383 + t507;
t237 = Ifges(4,4) * t245;
t410 = t298 * Ifges(4,5);
t414 = t246 * Ifges(4,1);
t143 = -t237 + t410 + t414;
t495 = t237 / 0.2e1 + t145 * mrSges(4,3) - t143 / 0.2e1 - t220 * mrSges(4,2) - t410 / 0.2e1;
t100 = Ifges(5,4) * t338 - t189 * Ifges(5,2) + t245 * Ifges(5,6);
t335 = pkin(3) * t298 - qJD(4) + t145;
t328 = qJ(5) * t338 + t335;
t56 = pkin(4) * t189 - t328;
t97 = Ifges(6,5) * t338 + t245 * Ifges(6,6) + t189 * Ifges(6,3);
t493 = t56 * mrSges(6,1) - t100 / 0.2e1 + t97 / 0.2e1 - t335 * mrSges(5,1) - t51 * mrSges(6,2) - t58 * mrSges(5,3);
t42 = qJD(6) * t118 + t163 * t313 + t164 * t316;
t470 = t42 / 0.2e1;
t43 = -qJD(6) * t494 + t163 * t316 - t164 * t313;
t469 = t43 / 0.2e1;
t489 = Ifges(6,2) + Ifges(5,3);
t400 = t311 * t315;
t280 = -pkin(8) * t400 + t360;
t271 = t280 * qJD(2);
t255 = qJD(1) * t271;
t486 = t255 * mrSges(3,2);
t481 = mrSges(3,1) * t336 - mrSges(4,1) * t245 - mrSges(4,2) * t246 - mrSges(3,3) * t372;
t422 = Ifges(6,6) * t310;
t423 = Ifges(5,6) * t310;
t427 = Ifges(5,5) * t312;
t429 = Ifges(6,4) * t312;
t479 = t422 + t429 - t423 + t427;
t426 = Ifges(6,5) * t310;
t431 = Ifges(5,4) * t310;
t478 = t312 * t492 + t426 - t431;
t477 = Ifges(4,1) * t505 + Ifges(4,4) * t454;
t269 = qJD(2) * t332;
t327 = qJD(1) * t269;
t95 = -qJD(3) * t146 - t255 * t314 + t317 * t327;
t425 = Ifges(6,5) * t312;
t343 = Ifges(6,3) * t310 + t425;
t430 = Ifges(5,4) * t312;
t346 = -Ifges(5,2) * t310 + t430;
t349 = mrSges(6,1) * t310 - mrSges(6,3) * t312;
t350 = mrSges(5,1) * t310 + mrSges(5,2) * t312;
t439 = t312 / 0.2e1;
t440 = t310 / 0.2e1;
t441 = -t310 / 0.2e1;
t455 = t338 / 0.2e1;
t457 = t189 / 0.2e1;
t458 = -t189 / 0.2e1;
t474 = -t335 * t350 + t343 * t457 + t346 * t458 + t56 * t349 + t455 * t478 + (-t310 * t51 + t312 * t50) * mrSges(6,2) + (-t310 * t58 - t312 * t57) * mrSges(5,3) + t100 * t441 + t97 * t440 + t482 * t439 - t495;
t472 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t454;
t471 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t454;
t467 = Ifges(4,5) * t358 / 0.2e1 + t477;
t466 = -t118 / 0.2e1;
t465 = t118 / 0.2e1;
t464 = -t494 / 0.2e1;
t463 = t494 / 0.2e1;
t462 = -t163 / 0.2e1;
t461 = t163 / 0.2e1;
t460 = t164 / 0.2e1;
t456 = -t338 / 0.2e1;
t453 = t199 / 0.2e1;
t447 = -t238 / 0.2e1;
t446 = t238 / 0.2e1;
t445 = -t245 / 0.2e1;
t444 = t245 / 0.2e1;
t38 = Ifges(7,5) * t42;
t37 = Ifges(7,6) * t43;
t438 = pkin(9) * t314;
t94 = -t221 * t384 + t233 * t383 + t317 * t255 + t314 * t327;
t75 = qJ(4) * t358 + qJD(4) * t298 + t94;
t272 = t388 * qJD(2);
t256 = qJD(1) * t272;
t78 = t199 * pkin(3) - t198 * qJ(4) - t246 * qJD(4) + t256;
t28 = t310 * t78 + t312 * t75;
t434 = Ifges(3,4) * t315;
t433 = Ifges(3,4) * t318;
t428 = Ifges(7,4) * t494;
t424 = Ifges(3,2) * t315;
t418 = t198 * Ifges(4,4);
t415 = t245 * Ifges(4,6);
t413 = t246 * Ifges(4,5);
t412 = t267 * mrSges(3,3);
t270 = t388 * qJD(1);
t411 = t270 * mrSges(3,3);
t408 = t298 * Ifges(4,3);
t407 = t315 * Ifges(3,1);
t122 = mrSges(6,1) * t189 - mrSges(6,3) * t338;
t49 = -mrSges(7,1) * t118 + mrSges(7,2) * t494;
t406 = t122 - t49;
t277 = t314 * t400 - t317 * t405;
t386 = qJD(2) * t311;
t369 = t318 * t386;
t211 = -qJD(3) * t277 + t317 * t369;
t278 = t314 * t405 + t317 * t400;
t212 = qJD(3) * t278 + t314 * t369;
t108 = t212 * pkin(3) - t211 * qJ(4) - t278 * qJD(4) + t272;
t113 = -t259 * t384 + t261 * t383 + t314 * t269 + t317 * t271;
t385 = qJD(2) * t315;
t96 = (qJ(4) * t385 - qJD(4) * t318) * t311 + t113;
t36 = t310 * t108 + t312 * t96;
t404 = Ifges(3,6) * qJD(2);
t135 = -mrSges(6,2) * t189 + mrSges(6,3) * t245;
t136 = -mrSges(5,2) * t245 - mrSges(5,3) * t189;
t395 = t135 + t136;
t137 = mrSges(5,1) * t245 - mrSges(5,3) * t338;
t138 = -mrSges(6,1) * t245 + mrSges(6,2) * t338;
t394 = t137 - t138;
t258 = -pkin(2) * t405 - t280;
t155 = t277 * pkin(3) - t278 * qJ(4) + t258;
t175 = t317 * t259 + t314 * t261;
t160 = -qJ(4) * t399 + t175;
t81 = t310 * t155 + t312 * t160;
t149 = t225 * t316 - t226 * t313;
t275 = t337 * qJD(6);
t184 = -t275 * t314 + t285 * t383;
t393 = t149 - t184;
t150 = t225 * t313 + t226 * t316;
t262 = t285 * t314;
t183 = qJD(6) * t262 + t337 * t383;
t392 = t150 - t183;
t151 = t285 * t245;
t276 = t285 * qJD(6);
t391 = -t151 + t276;
t152 = t337 * t245;
t390 = -t152 + t275;
t123 = mrSges(5,1) * t189 + mrSges(5,2) * t338;
t202 = mrSges(4,1) * t298 - mrSges(4,3) * t246;
t389 = t202 - t123;
t174 = -t314 * t259 + t317 * t261;
t382 = qJD(5) * t310;
t7 = -Ifges(7,3) * t199 + t37 + t38;
t71 = t277 * qJ(5) + t81;
t376 = Ifges(4,5) * t198 - Ifges(4,6) * t199 + Ifges(4,3) * t358;
t162 = pkin(3) * t399 - t174;
t370 = t311 * t385;
t367 = Ifges(3,5) * t405;
t366 = Ifges(3,6) * t405;
t362 = qJ(5) * t310 + pkin(3);
t27 = -t310 * t75 + t312 * t78;
t93 = t163 * mrSges(5,1) + t164 * mrSges(5,2);
t112 = -t199 * mrSges(6,1) + t164 * mrSges(6,2);
t92 = t163 * mrSges(6,1) - t164 * mrSges(6,3);
t35 = t108 * t312 - t310 * t96;
t85 = t172 * t312 - t139;
t80 = t155 * t312 - t310 * t160;
t243 = t292 * t312 - t301;
t19 = t199 * qJ(5) + t245 * qJD(5) + t28;
t26 = t212 * qJ(5) + t277 * qJD(5) + t36;
t114 = -t259 * t383 - t261 * t384 + t317 * t269 - t314 * t271;
t11 = -pkin(10) * t164 - t199 * t468 - t27;
t12 = pkin(10) * t163 + t19;
t1 = qJD(6) * t5 + t11 * t313 + t12 * t316;
t2 = -qJD(6) * t6 + t11 * t316 - t12 * t313;
t354 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t10 = -t43 * mrSges(7,1) + t42 * mrSges(7,2);
t340 = -t27 * t310 + t28 * t312;
t210 = t312 * t278 - t310 * t399;
t48 = -pkin(10) * t210 - t277 * t468 - t80;
t209 = t310 * t278 + t311 * t396;
t52 = pkin(10) * t209 + t71;
t13 = -t313 * t52 + t316 * t48;
t14 = t313 * t48 + t316 * t52;
t83 = -mrSges(7,2) * t238 + mrSges(7,3) * t118;
t84 = mrSges(7,1) * t238 - mrSges(7,3) * t494;
t339 = -t313 * t84 + t316 * t83;
t133 = t209 * t316 - t210 * t313;
t134 = t209 * t313 + t210 * t316;
t333 = qJ(5) * t210 - t162;
t326 = pkin(3) * t370 + t114;
t79 = -pkin(3) * t358 - t95;
t323 = (t366 + (Ifges(3,2) * t318 + t434) * t311) * qJD(1);
t180 = t211 * t312 + t310 * t370;
t322 = qJ(5) * t180 + qJD(5) * t210 + t326;
t25 = t163 * pkin(4) - t164 * qJ(5) - qJD(5) * t338 + t79;
t299 = Ifges(3,4) * t371;
t296 = Ifges(3,5) * t357;
t288 = -pkin(4) * t312 - t362;
t279 = t312 * t468 + t362;
t265 = -mrSges(3,2) * t336 + mrSges(3,3) * t371;
t263 = t337 * t314;
t260 = t334 * t314;
t223 = -t243 + t309;
t219 = t330 * t314;
t217 = Ifges(3,1) * t372 + Ifges(3,5) * t336 + t299;
t216 = t323 + t404;
t206 = t310 * t379 + t402;
t201 = -mrSges(4,2) * t298 - mrSges(4,3) * t245;
t179 = t211 * t310 - t312 * t370;
t171 = -mrSges(4,2) * t358 - mrSges(4,3) * t199;
t170 = mrSges(4,1) * t358 - mrSges(4,3) * t198;
t141 = t408 + t413 - t415;
t127 = mrSges(4,1) * t199 + mrSges(4,2) * t198;
t116 = -t199 * Ifges(4,2) + Ifges(4,6) * t358 + t418;
t115 = Ifges(7,4) * t118;
t111 = mrSges(5,1) * t199 - mrSges(5,3) * t164;
t110 = -mrSges(5,2) * t199 - mrSges(5,3) * t163;
t109 = -mrSges(6,2) * t163 + mrSges(6,3) * t199;
t91 = t245 * t341 + t146;
t82 = pkin(4) * t209 - t333;
t72 = -pkin(4) * t277 - t80;
t70 = -pkin(4) * t246 - t85;
t69 = t245 * t476 - t146;
t64 = t164 * Ifges(5,4) - t163 * Ifges(5,2) + t199 * Ifges(5,6);
t63 = t164 * Ifges(6,4) + t199 * Ifges(6,2) + t163 * Ifges(6,6);
t62 = t164 * Ifges(5,5) - t163 * Ifges(5,6) + t199 * Ifges(5,3);
t61 = t164 * Ifges(6,5) + t199 * Ifges(6,6) + t163 * Ifges(6,3);
t60 = -t209 * t468 + t333;
t54 = -qJD(6) * t134 + t179 * t316 - t180 * t313;
t53 = qJD(6) * t133 + t179 * t313 + t180 * t316;
t44 = -t189 * t468 + t328;
t41 = Ifges(7,1) * t494 + t238 * Ifges(7,5) + t115;
t40 = t118 * Ifges(7,2) + t238 * Ifges(7,6) + t428;
t33 = pkin(4) * t179 - t322;
t31 = mrSges(7,2) * t199 + mrSges(7,3) * t43;
t30 = -mrSges(7,1) * t199 - mrSges(7,3) * t42;
t29 = -pkin(4) * t212 - t35;
t24 = -t179 * t468 + t322;
t21 = -pkin(4) * t199 - t27;
t20 = pkin(10) * t179 + t26;
t18 = -pkin(10) * t180 - t212 * t468 - t35;
t17 = pkin(5) * t163 + t25;
t4 = -qJD(6) * t14 + t18 * t316 - t20 * t313;
t3 = qJD(6) * t13 + t18 * t313 + t20 * t316;
t8 = [(-m(3) * t280 + m(4) * t258 - mrSges(3,1) * t405 + mrSges(4,1) * t277 + mrSges(4,2) * t278) * t256 + m(4) * (t113 * t146 + t114 * t145 + t174 * t95 + t175 * t94) + (Ifges(4,1) * t278 - Ifges(4,4) * t277 - Ifges(4,5) * t399) * t505 - (t142 + t39) * t212 / 0.2e1 + ((-t424 + t433) * t365 + (Ifges(3,1) * t318 - t434) * t364 - 0.2e1 * pkin(1) * (mrSges(3,1) * t315 + mrSges(3,2) * t318) * t380) * t311 ^ 2 + m(3) * (t255 * t388 + t270 * t271) + (t255 * t399 + t256 * t400 - t280 * t357 - t358 * t388) * mrSges(3,3) + ((t315 * (Ifges(4,5) * t278 - Ifges(4,6) * t277 - Ifges(4,3) * t399) + t318 * (t367 + (t407 + t433) * t311)) * qJD(1) + t336 * (Ifges(3,5) * t318 - Ifges(3,6) * t315) + t318 * t217 + t315 * t141) * t386 / 0.2e1 + (Ifges(4,4) * t278 + Ifges(7,5) * t134 - Ifges(4,6) * t399 + Ifges(7,6) * t133 + (-Ifges(4,2) - Ifges(7,3)) * t277) * t454 - (t323 + t216) * t370 / 0.2e1 + (t63 + t62) * t277 / 0.2e1 - (t116 + t7) * t277 / 0.2e1 + (t99 + t98) * t212 / 0.2e1 - t376 * t399 / 0.2e1 + t95 * (-mrSges(4,1) * t399 - t278 * mrSges(4,3)) + t94 * (mrSges(4,2) * t399 - t277 * mrSges(4,3)) + t405 * (-Ifges(3,6) * t358 + t296) / 0.2e1 + (-m(3) * t267 + m(4) * t220 - t481) * t272 + t482 * t180 / 0.2e1 + t146 * (-mrSges(4,2) * t370 - mrSges(4,3) * t212) + t298 * (Ifges(4,5) * t211 - Ifges(4,6) * t212 + Ifges(4,3) * t370) / 0.2e1 + t246 * (Ifges(4,1) * t211 - Ifges(4,4) * t212 + Ifges(4,5) * t370) / 0.2e1 + t145 * (mrSges(4,1) * t370 - mrSges(4,3) * t211) - t370 * t411 - t369 * t412 + (t180 * t492 + t212 * t491) * t455 + (t210 * t492 + t277 * t491) * t460 + (Ifges(6,5) * t180 + Ifges(6,6) * t212) * t457 + (Ifges(6,5) * t210 + Ifges(6,6) * t277) * t461 + m(6) * (t19 * t71 + t21 * t72 + t25 * t82 + t26 * t51 + t29 * t50 + t33 * t56) + m(7) * (t1 * t14 + t13 * t2 - t17 * t60 + t24 * t44 + t3 * t6 + t4 * t5) + t487 * t210 / 0.2e1 + (t180 * t491 + t212 * t489) * t444 + (t210 * t491 + t277 * t489) * t453 + (-Ifges(5,2) * t458 + Ifges(6,3) * t457 - t444 * t504 + t490 * t455 + t493) * t179 + (-t28 * mrSges(5,3) - t19 * mrSges(6,2) - t64 / 0.2e1 + t79 * mrSges(5,1) + t25 * mrSges(6,1) + t61 / 0.2e1 + Ifges(6,3) * t461 - Ifges(5,2) * t462 + t490 * t460 - t504 * t453) * t209 - t405 * t486 - t326 * t123 + (-t180 * t335 + t210 * t79 - t212 * t58 - t277 * t28) * mrSges(5,2) + m(5) * (t162 * t79 + t27 * t80 + t28 * t81 + t326 * t335 + t35 * t57 + t36 * t58) + t2 * (-mrSges(7,1) * t277 - mrSges(7,3) * t134) + t27 * (mrSges(5,1) * t277 - mrSges(5,3) * t210) + t21 * (-mrSges(6,1) * t277 + mrSges(6,2) * t210) + t1 * (mrSges(7,2) * t277 + mrSges(7,3) * t133) + t271 * t265 + t258 * t127 + t5 * (-mrSges(7,1) * t212 - mrSges(7,3) * t53) + t6 * (mrSges(7,2) * t212 + mrSges(7,3) * t54) + t220 * (mrSges(4,1) * t212 + mrSges(4,2) * t211) + t211 * t143 / 0.2e1 + t57 * (mrSges(5,1) * t212 - mrSges(5,3) * t180) + t50 * (-mrSges(6,1) * t212 + mrSges(6,2) * t180) + t113 * t201 + t114 * t202 + t174 * t170 + t175 * t171 + (Ifges(5,4) * t180 + Ifges(5,6) * t212) * t458 + (Ifges(5,4) * t210 + Ifges(5,6) * t277) * t462 + t162 * t93 + t26 * t135 + t36 * t136 + t35 * t137 + t29 * t138 - t17 * (-mrSges(7,1) * t133 + mrSges(7,2) * t134) + t33 * t122 + t72 * t112 + t71 * t109 + t81 * t110 + t80 * t111 + t82 * t92 + t3 * t83 + t4 * t84 + t60 * t10 + t53 * t41 / 0.2e1 + t44 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t54 * t40 / 0.2e1 + t24 * t49 + t13 * t30 + t14 * t31 + (-t180 * t56 + t19 * t277 - t210 * t25 + t212 * t51) * mrSges(6,3) + (Ifges(4,4) * t211 - Ifges(4,2) * t212 + Ifges(4,6) * t370) * t445 + (Ifges(7,5) * t53 + Ifges(7,6) * t54 - Ifges(7,3) * t212) * t446 + (Ifges(7,1) * t53 + Ifges(7,4) * t54 - Ifges(7,5) * t212) * t463 + (Ifges(7,4) * t53 + Ifges(7,2) * t54 - Ifges(7,6) * t212) * t465 + t278 * t467 + (Ifges(7,4) * t134 + Ifges(7,2) * t133 - Ifges(7,6) * t277) * t469 + (Ifges(7,1) * t134 + Ifges(7,4) * t133 - Ifges(7,5) * t277) * t470 + t134 * t471 + t133 * t472; -t486 + (((-m(4) * t146 - t201) * pkin(9) - t506) * t314 + (t414 / 0.2e1 + (-m(4) * t145 - m(5) * t335 - t389) * pkin(9) + t479 * t444 + t474) * t317) * qJD(3) + ((t411 - t404 / 0.2e1 + t216 / 0.2e1 - t141 / 0.2e1 + qJD(2) * (Ifges(4,5) * t314 + Ifges(4,6) * t317) / 0.2e1 - t145 * mrSges(4,1) - t408 / 0.2e1 - t413 / 0.2e1 + t415 / 0.2e1 + t146 * mrSges(4,2) + (t366 / 0.2e1 + (pkin(1) * mrSges(3,1) + t434 / 0.2e1) * t311) * qJD(1)) * t315 + (t412 - Ifges(3,5) * qJD(2) / 0.2e1 - t299 / 0.2e1 - t217 / 0.2e1 + (-t367 / 0.2e1 + (pkin(1) * mrSges(3,2) - t407 / 0.2e1 + t424 / 0.2e1) * t311) * qJD(1) + (-t414 / 0.2e1 + t495) * t317 + t506 * t314) * t318) * t387 + t481 * t270 + (-t482 / 0.2e1 + mrSges(5,2) * t335 - mrSges(6,2) * t50 + mrSges(5,3) * t57 + mrSges(6,3) * t56 + Ifges(5,4) * t457 + Ifges(6,5) * t458 + t445 * t491 + t456 * t492) * t226 + (t206 - t105) * t137 + (mrSges(7,1) * t393 - mrSges(7,2) * t392) * t44 + (t1 * t262 - t2 * t263 + t392 * t5 - t393 * t6) * mrSges(7,3) - m(4) * (t145 * t181 + t146 * t182 + t220 * t270) + m(4) * (pkin(9) * t317 * t94 - pkin(2) * t256 - t438 * t95) + m(5) * (t206 * t57 + t207 * t58 + t243 * t27 + t244 * t28 + t438 * t79) + (t28 * mrSges(5,2) - t19 * mrSges(6,3) - t256 * mrSges(4,1) + t418 / 0.2e1 - t27 * mrSges(5,1) + t21 * mrSges(6,1) + t38 / 0.2e1 + t37 / 0.2e1 + t116 / 0.2e1 - t62 / 0.2e1 - t63 / 0.2e1 + t7 / 0.2e1 + pkin(9) * t171 + t94 * mrSges(4,3) - t378 * t164 - t377 * t163 + (-Ifges(7,3) / 0.2e1 - t355) * t199 + t354) * t317 + (t207 - t106) * t136 + (t256 * mrSges(4,2) + t343 * t461 + t346 * t462 + t467 + t79 * t350 + t25 * t349 - t95 * mrSges(4,3) + t61 * t440 + t64 * t441 + (-t170 + t93) * pkin(9) + (-t27 * t312 - t28 * t310) * mrSges(5,3) + (-t19 * t310 + t21 * t312) * mrSges(6,2) + t478 * t460 + t479 * t453 + t487 * t439 + t477) * t314 + (t183 / 0.2e1 - t150 / 0.2e1) * t41 - m(5) * (t105 * t57 + t106 * t58 + t331 * t335) + t331 * t123 + t296 + (t184 / 0.2e1 - t149 / 0.2e1) * t40 - t17 * (-mrSges(7,1) * t262 + mrSges(7,2) * t263) - t267 * t265 - t256 * mrSges(3,1) + t260 * t92 + t243 * t111 + t244 * t110 + t222 * t109 + t223 * t112 + t219 * t10 - t182 * t201 - t181 * t202 + t124 * t30 + t125 * t31 - pkin(2) * t127 + (Ifges(7,5) * t183 + Ifges(7,6) * t184) * t446 + (Ifges(7,5) * t150 + Ifges(7,6) * t149) * t447 + (Ifges(7,5) * t263 + Ifges(7,6) * t262) * t454 + t496 * t122 + t497 * t138 + t498 * t49 + t499 * t135 + (t19 * t222 + t21 * t223 + t25 * t260 + t496 * t56 + t497 * t50 + t499 * t51) * m(6) + t502 * t84 + t503 * t83 + (t1 * t125 + t124 * t2 - t17 * t219 + t44 * t498 + t5 * t502 + t503 * t6) * m(7) + (-Ifges(5,2) * t457 + Ifges(6,3) * t458 - t445 * t504 + t456 * t490 - t493) * t225 + (Ifges(7,1) * t183 + Ifges(7,4) * t184) * t463 + (Ifges(7,1) * t150 + Ifges(7,4) * t149) * t464 + (Ifges(7,4) * t183 + Ifges(7,2) * t184) * t465 + (Ifges(7,4) * t150 + Ifges(7,2) * t149) * t466 + (Ifges(7,4) * t263 + Ifges(7,2) * t262) * t469 + (Ifges(7,1) * t263 + Ifges(7,4) * t262) * t470 + t263 * t471 + t262 * t472; -t475 * t246 + (-t50 * t70 - t51 * t68 - t56 * t91 + t25 * t288 + (qJ(4) * t19 + qJD(4) * t51) * t312 + (qJ(4) * t21 + qJD(4) * t50 - qJD(5) * t56) * t310) * m(6) - t406 * t382 + (t19 * t312 + t21 * t310) * mrSges(6,2) + t79 * (-mrSges(5,1) * t312 + mrSges(5,2) * t310) + t25 * (-mrSges(6,1) * t312 - mrSges(6,3) * t310) - t312 * t61 / 0.2e1 + t376 + (mrSges(7,1) * t391 - mrSges(7,2) * t390) * t44 + (-t310 * t394 + t312 * t395) * qJD(4) + ((Ifges(4,1) / 0.2e1 - t355) * t246 + (t427 / 0.2e1 - t423 / 0.2e1 + t429 / 0.2e1 + t422 / 0.2e1) * t245 + t474) * t245 + t340 * mrSges(5,3) + (t335 * t146 - t57 * t85 - t58 * t86 - pkin(3) * t79 + (-t310 * t57 + t312 * t58) * qJD(4) + t340 * qJ(4)) * m(5) - t17 * (mrSges(7,1) * t337 + mrSges(7,2) * t285) + (Ifges(7,5) * t285 - Ifges(7,6) * t337) * t454 + (Ifges(7,4) * t285 - Ifges(7,2) * t337) * t469 + (Ifges(7,1) * t285 - Ifges(7,4) * t337) * t470 + (-t1 * t337 - t2 * t285 + t390 * t5 - t391 * t6) * mrSges(7,3) - t337 * t472 + (-Ifges(7,5) * t152 - Ifges(7,6) * t151) * t447 + (-Ifges(7,1) * t152 - Ifges(7,4) * t151) * t464 + (-Ifges(7,4) * t152 - Ifges(7,2) * t151) * t466 + (-t275 / 0.2e1 + t152 / 0.2e1) * t41 + (-t276 / 0.2e1 + t151 / 0.2e1) * t40 + (-Ifges(7,5) * t275 - Ifges(7,6) * t276) * t446 + (-Ifges(7,1) * t275 - Ifges(7,4) * t276) * t463 + (-Ifges(7,4) * t275 - Ifges(7,2) * t276) * t465 + t288 * t92 + t279 * t10 + t214 * t30 + t215 * t31 - t145 * t201 - t68 * t135 - t86 * t136 - t85 * t137 - t70 * t138 - t91 * t122 - pkin(3) * t93 - t94 * mrSges(4,2) + t95 * mrSges(4,1) - t69 * t49 + ((t109 + t110) * t312 + (-t111 + t112) * t310) * qJ(4) + t64 * t439 + t500 * t84 + t501 * t83 + (t1 * t215 - t17 * t279 + t2 * t214 + t501 * t6 + t500 * t5 + (t382 - t69) * t44) * m(7) + t487 * t440 + (t310 * t491 + t312 * t504) * t453 + t389 * t146 + (t310 * t492 - t425 + t430) * t460 + (-Ifges(6,3) * t312 + t426) * t461 + (Ifges(5,2) * t312 + t431) * t462 + t285 * t471; -t494 * t84 + t118 * t83 + t394 * t338 + t395 * t189 - t10 + t92 + t93 + (t118 * t6 - t494 * t5 + t17) * m(7) + (t189 * t51 - t338 * t50 + t25) * m(6) + (t189 * t58 + t338 * t57 + t79) * m(5); t316 * t30 + t313 * t31 + t406 * t338 + t339 * qJD(6) + (-t135 - t339) * t245 + t112 + (t1 * t313 - t338 * t44 + t2 * t316 + t238 * (-t313 * t5 + t316 * t6)) * m(7) + (-t245 * t51 + t338 * t56 + t21) * m(6); -t44 * (mrSges(7,1) * t494 + mrSges(7,2) * t118) + (Ifges(7,1) * t118 - t428) * t464 + t40 * t463 + (Ifges(7,5) * t118 - Ifges(7,6) * t494) * t447 - t5 * t83 + t6 * t84 + (t118 * t5 + t494 * t6) * mrSges(7,3) + t354 + t7 + (-Ifges(7,2) * t494 + t115 + t41) * t466;];
tauc  = t8(:);
