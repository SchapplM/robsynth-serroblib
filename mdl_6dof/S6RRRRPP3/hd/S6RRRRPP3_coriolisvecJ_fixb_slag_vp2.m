% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:05:39
% EndTime: 2018-11-23 18:05:53
% DurationCPUTime: 14.24s
% Computational Cost: add. (9215->643), mult. (22821->797), div. (0->0), fcn. (15402->6), ass. (0->295)
t508 = Ifges(5,1) + Ifges(7,3);
t478 = Ifges(7,4) + Ifges(6,5);
t477 = Ifges(5,5) + Ifges(7,5);
t507 = -Ifges(6,3) - Ifges(7,2);
t457 = Ifges(6,4) - t477;
t513 = Ifges(5,6) - t478;
t299 = sin(qJ(3));
t300 = sin(qJ(2));
t302 = cos(qJ(3));
t303 = cos(qJ(2));
t267 = t299 * t303 + t302 * t300;
t255 = t267 * qJD(1);
t301 = cos(qJ(4));
t266 = t299 * t300 - t302 * t303;
t313 = t266 * qJD(3);
t222 = -qJD(2) * t266 - t313;
t305 = qJD(1) * t222;
t368 = qJD(2) + qJD(3);
t339 = t301 * t368;
t298 = sin(qJ(4));
t371 = qJD(4) * t298;
t137 = -qJD(4) * t339 + t255 * t371 - t301 * t305;
t450 = -t137 / 0.2e1;
t449 = t137 / 0.2e1;
t231 = t301 * t255 + t298 * t368;
t372 = qJD(4) * t231;
t138 = t298 * t305 + t372;
t448 = -t138 / 0.2e1;
t447 = t138 / 0.2e1;
t223 = t368 * t267;
t205 = t223 * qJD(1);
t446 = t205 / 0.2e1;
t524 = -Ifges(5,4) + Ifges(7,6);
t523 = Ifges(6,6) - Ifges(7,6);
t370 = qJD(4) * t301;
t376 = qJD(1) * t303;
t377 = qJD(1) * t300;
t254 = -t299 * t377 + t302 * t376;
t391 = t254 * t301;
t522 = t370 - t391;
t290 = -pkin(2) * t303 - pkin(1);
t278 = qJD(1) * t290;
t504 = t368 * Ifges(4,5);
t521 = -t504 / 0.2e1 - t278 * mrSges(4,2);
t250 = qJD(4) - t254;
t425 = pkin(4) + qJ(6);
t451 = -pkin(8) - pkin(7);
t280 = t451 * t300;
t269 = qJD(1) * t280;
t261 = qJD(2) * pkin(2) + t269;
t282 = t451 * t303;
t270 = qJD(1) * t282;
t382 = t302 * t270;
t214 = t299 * t261 - t382;
t184 = pkin(9) * t368 + t214;
t312 = -t254 * pkin(3) - t255 * pkin(9) + t278;
t78 = t184 * t298 - t301 * t312;
t317 = pkin(5) * t231 + t78;
t491 = qJD(5) + t317;
t36 = -t250 * t425 + t491;
t503 = t368 * Ifges(4,6);
t465 = -qJD(5) - t78;
t52 = -pkin(4) * t250 - t465;
t79 = t301 * t184 + t298 * t312;
t53 = -t250 * qJ(5) - t79;
t520 = -t278 * mrSges(4,1) + t78 * mrSges(5,1) + t79 * mrSges(5,2) - t52 * mrSges(6,2) + t53 * mrSges(6,3) + t36 * mrSges(7,3) + t503 / 0.2e1;
t458 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t517 = -t205 / 0.2e1;
t516 = t36 * mrSges(7,1);
t515 = t78 * mrSges(5,3);
t226 = Ifges(7,6) * t231;
t230 = t255 * t298 - t339;
t404 = t231 * Ifges(6,6);
t502 = -t230 * t507 + t250 * t478 + t226 - t404;
t228 = Ifges(5,4) * t230;
t412 = Ifges(7,6) * t230;
t501 = t231 * t508 + t477 * t250 - t228 + t412;
t256 = t299 * t270;
t213 = t302 * t261 + t256;
t360 = qJD(2) * t451;
t338 = qJD(1) * t360;
t262 = t300 * t338;
t318 = t303 * t338;
t125 = qJD(3) * t213 + t302 * t262 + t299 * t318;
t500 = qJD(4) * t312 + t125;
t89 = t205 * pkin(3) + (pkin(9) * t313 + (t300 * pkin(2) + pkin(9) * t266) * qJD(2)) * qJD(1);
t13 = -t184 * t371 + t298 * t89 + t301 * t500;
t8 = -qJ(5) * t205 - qJD(5) * t250 - t13;
t4 = -pkin(5) * t138 - t8;
t512 = -t4 * mrSges(7,1) + Ifges(5,2) * t447 + Ifges(5,6) * t517 + t446 * t478 + t448 * t507 + (Ifges(5,4) + t523) * t449;
t14 = -t184 * t370 - t298 * t500 + t301 * t89;
t2 = -pkin(5) * t137 - qJD(6) * t250 - t205 * t425 - t14;
t511 = mrSges(7,1) * t2 + Ifges(6,4) * t517 + Ifges(6,6) * t448 + t477 * t446 + t524 * t447 + (Ifges(6,2) + t508) * t450;
t490 = t230 * pkin(5) - qJD(6);
t40 = -t53 - t490;
t510 = -t40 * mrSges(7,1) - t79 * mrSges(5,3);
t340 = pkin(4) * t371 - qJD(5) * t298;
t392 = t254 * t298;
t499 = (-qJ(5) * qJD(4) - qJD(6)) * t301 + t340 + (t371 - t392) * qJ(6);
t498 = pkin(5) * t522 + t255 * t425;
t411 = Ifges(7,6) * t298;
t416 = Ifges(5,4) * t298;
t497 = t301 * t508 + t411 - t416;
t410 = Ifges(7,6) * t301;
t413 = Ifges(6,6) * t301;
t496 = -t298 * t507 + t410 - t413;
t394 = t222 * t298;
t316 = t267 * t370 + t394;
t227 = Ifges(6,6) * t230;
t118 = Ifges(6,4) * t250 - Ifges(6,2) * t231 + t227;
t405 = t231 * Ifges(5,4);
t119 = -Ifges(5,2) * t230 + Ifges(5,6) * t250 + t405;
t493 = t301 * t118 + t298 * t119;
t488 = -t230 * t513 - t231 * t457 + t250 * t458;
t487 = -t298 * t513 - t301 * t457;
t11 = -pkin(4) * t205 - t14;
t486 = t11 * t298 + t52 * t370 + t53 * t371;
t183 = -pkin(3) * t368 - t213;
t329 = -mrSges(7,2) * t301 + mrSges(7,3) * t298;
t330 = -mrSges(6,2) * t298 - mrSges(6,3) * t301;
t331 = mrSges(5,1) * t298 + mrSges(5,2) * t301;
t311 = -t231 * qJ(5) + t183;
t44 = t230 * t425 + t311;
t70 = t230 * pkin(4) + t311;
t485 = t183 * t331 + t329 * t44 + t330 * t70;
t482 = t8 * mrSges(6,1);
t369 = qJD(1) * qJD(2);
t351 = t300 * t369;
t481 = pkin(2) * t351;
t480 = t14 * mrSges(5,3);
t208 = pkin(3) * t255 - pkin(9) * t254;
t177 = pkin(2) * t377 + t208;
t217 = t269 * t302 + t256;
t104 = t177 * t301 - t298 * t217;
t288 = pkin(2) * t299 + pkin(9);
t409 = pkin(2) * qJD(3);
t367 = t302 * t409;
t473 = t288 * t370 + t298 * t367 + t104 + t498;
t112 = t208 * t301 - t298 * t213;
t472 = pkin(9) * t370 + t112 + t498;
t126 = qJD(3) * t214 + t262 * t299 - t302 * t318;
t38 = mrSges(5,1) * t138 - mrSges(5,2) * t137;
t471 = m(5) * t126 + t38;
t105 = t298 * t177 + t301 * t217;
t247 = t255 * qJ(5);
t335 = -pkin(5) * t392 + t247;
t470 = -t335 - t105 + t301 * t367 + (-pkin(5) - t288) * t371;
t113 = t298 * t208 + t301 * t213;
t469 = -t335 - t113 + (-pkin(5) - pkin(9)) * t371;
t334 = pkin(4) * t392 - qJ(5) * t391;
t124 = t334 + t214;
t468 = -t124 + t499;
t216 = t269 * t299 - t382;
t128 = t216 + t334;
t294 = t299 * t409;
t467 = -t128 + t294 + t499;
t464 = -t79 + t490;
t249 = -qJ(5) * t370 + t340;
t463 = -t124 + t249;
t402 = t255 * mrSges(4,3);
t381 = -mrSges(4,1) * t368 + mrSges(5,1) * t230 + mrSges(5,2) * t231 + t402;
t462 = t302 * t280 + t282 * t299;
t406 = t14 * t298;
t407 = t13 * t301;
t426 = t301 * t8;
t461 = m(5) * (t370 * t78 - t371 * t79 - t406 + t407) + m(6) * (-t426 + t486);
t455 = t137 * t457 - t138 * t513 + t458 * t205;
t454 = t14 * mrSges(5,1) - t13 * mrSges(5,2) + t11 * mrSges(6,2) + t4 * mrSges(7,2) - t8 * mrSges(6,3) - t2 * mrSges(7,3);
t453 = m(4) / 0.2e1;
t443 = -t230 / 0.2e1;
t442 = t230 / 0.2e1;
t441 = -t231 / 0.2e1;
t440 = t231 / 0.2e1;
t439 = -t250 / 0.2e1;
t438 = t250 / 0.2e1;
t437 = -t254 / 0.2e1;
t436 = t254 / 0.2e1;
t434 = t255 / 0.2e1;
t429 = pkin(2) * t302;
t428 = pkin(4) * t255;
t56 = mrSges(6,1) * t138 - mrSges(6,3) * t205;
t60 = -mrSges(5,2) * t205 - mrSges(5,3) * t138;
t424 = -t56 + t60;
t58 = -t137 * mrSges(6,1) + t205 * mrSges(6,2);
t59 = mrSges(5,1) * t205 + mrSges(5,3) * t137;
t423 = t58 - t59;
t422 = mrSges(4,3) * t299;
t421 = mrSges(5,3) * t230;
t420 = mrSges(5,3) * t231;
t419 = Ifges(3,4) * t300;
t418 = Ifges(3,4) * t303;
t417 = Ifges(4,4) * t255;
t415 = Ifges(5,4) * t301;
t414 = Ifges(6,6) * t298;
t403 = t254 * mrSges(4,3);
t401 = Ifges(3,5) * qJD(2);
t400 = Ifges(3,6) * qJD(2);
t399 = qJ(5) * t230;
t398 = qJD(2) * mrSges(3,1);
t397 = qJD(2) * mrSges(3,2);
t396 = t126 * t462;
t395 = t126 * t267;
t393 = t222 * t301;
t390 = t267 * t298;
t389 = t267 * t301;
t386 = t298 * t302;
t384 = t301 * t302;
t159 = mrSges(6,1) * t230 - mrSges(6,3) * t250;
t160 = -mrSges(7,1) * t230 + mrSges(7,2) * t250;
t380 = t159 - t160;
t162 = -mrSges(5,2) * t250 - t421;
t379 = t159 - t162;
t161 = mrSges(6,1) * t231 + mrSges(6,2) * t250;
t163 = mrSges(5,1) * t250 - t420;
t378 = t161 - t163;
t212 = t266 * pkin(3) - t267 * pkin(9) + t290;
t233 = t280 * t299 - t282 * t302;
t140 = t298 * t212 + t301 * t233;
t375 = qJD(2) * t300;
t374 = qJD(2) * t303;
t366 = pkin(2) * t375;
t359 = t267 * t371;
t345 = -t371 / 0.2e1;
t344 = t371 / 0.2e1;
t343 = -t370 / 0.2e1;
t342 = t370 / 0.2e1;
t341 = -qJ(5) * t298 - pkin(3);
t57 = -t138 * mrSges(7,1) + t205 * mrSges(7,2);
t55 = -t137 * mrSges(7,1) - t205 * mrSges(7,3);
t224 = t298 * t233;
t139 = t212 * t301 - t224;
t90 = -qJ(5) * t266 - t140;
t333 = pkin(4) * t390 - t462;
t327 = -Ifges(5,2) * t298 + t415;
t323 = -Ifges(6,2) * t301 + t414;
t319 = -qJ(5) * t301 + qJ(6) * t298;
t275 = -pkin(4) * t301 + t341;
t136 = pkin(3) * t223 - pkin(9) * t222 + t366;
t273 = t300 * t360;
t274 = t303 * t360;
t144 = qJD(3) * t462 + t273 * t302 + t274 * t299;
t21 = t136 * t301 - t298 * t144 - t212 * t371 - t233 * t370;
t315 = t359 - t393;
t314 = (Ifges(3,2) * t303 + t419) * qJD(1);
t20 = t298 * t136 + t301 * t144 + t212 * t370 - t233 * t371;
t260 = -t301 * t425 + t341;
t145 = qJD(3) * t233 + t273 * t299 - t302 * t274;
t16 = -qJ(5) * t223 - qJD(5) * t266 - t20;
t307 = pkin(4) * t316 + qJ(5) * t359 + t145;
t306 = qJ(5) * t137 - qJD(5) * t231 + t126;
t178 = Ifges(4,2) * t254 + t417 + t503;
t248 = Ifges(4,4) * t254;
t179 = Ifges(4,1) * t255 + t248 + t504;
t18 = pkin(4) * t138 + t306;
t6 = qJD(6) * t230 + t138 * t425 + t306;
t304 = (-t391 * t52 - t392 * t53 + t486) * mrSges(6,1) - (Ifges(4,1) * t254 - t417 + t488) * t255 / 0.2e1 + (t254 * t487 + t255 * t458) * t439 + t118 * t343 + t119 * t345 + (-t414 + t411) * t447 + (t248 + t179) * t437 + (-t391 * t78 + t392 * t79 + t407) * mrSges(5,3) + (Ifges(6,4) * t440 - Ifges(4,2) * t437 + Ifges(5,6) * t442 + t477 * t441 + t478 * t443 + t520) * t255 + (t323 * t440 + t327 * t442 + t497 * t441 + t496 * t443 - t485 + t521) * t254 + t416 * t448 + (-t410 + t415) * t450 + (t438 * t487 + t485 + (t496 / 0.2e1 - t327 / 0.2e1) * t230) * qJD(4) - t413 * t449 + (mrSges(5,2) * t126 - t6 * mrSges(7,2) - t18 * mrSges(6,3) - Ifges(6,2) * t449 - t446 * t457 + t508 * t450 + t511) * t298 + (-t126 * mrSges(5,1) + t18 * mrSges(6,2) - t6 * mrSges(7,3) + Ifges(5,2) * t448 + t446 * t513 + t507 * t447 - t512) * t301 + (t497 / 0.2e1 - t323 / 0.2e1) * t372 + t493 * t436 + Ifges(4,5) * t305 + t522 * t516 + t213 * t403 + t178 * t434 + t501 * (t342 - t391 / 0.2e1) + t502 * (t344 - t392 / 0.2e1) - t125 * mrSges(4,2) - t126 * mrSges(4,1) - Ifges(4,6) * t205;
t296 = t301 * pkin(5);
t295 = t298 * pkin(5);
t291 = Ifges(3,4) * t376;
t281 = pkin(9) * t301 + t296;
t279 = pkin(9) * t298 + t295;
t277 = mrSges(3,3) * t376 - t397;
t276 = -mrSges(3,3) * t377 + t398;
t265 = t288 * t301 + t296;
t264 = t288 * t298 + t295;
t263 = t275 - t429;
t253 = Ifges(3,1) * t377 + t291 + t401;
t252 = t314 + t400;
t242 = t260 - t429;
t238 = t249 + t294;
t236 = -mrSges(4,2) * t368 + t403;
t207 = -mrSges(4,1) * t254 + mrSges(4,2) * t255;
t193 = -mrSges(7,1) * t392 + mrSges(7,2) * t255;
t158 = mrSges(7,1) * t231 - mrSges(7,3) * t250;
t151 = -mrSges(6,2) * t230 - mrSges(6,3) * t231;
t149 = pkin(4) * t231 + t399;
t148 = -mrSges(7,2) * t231 + mrSges(7,3) * t230;
t143 = -qJ(5) * t389 + t333;
t106 = t267 * t319 + t333;
t97 = -pkin(4) * t266 - t139;
t87 = t231 * t425 + t399;
t73 = -t112 - t428;
t72 = -t247 - t113;
t63 = -t104 - t428;
t62 = -t247 - t105;
t61 = -pkin(5) * t390 - t90;
t51 = t224 + (pkin(5) * t267 - t212) * t301 - t425 * t266;
t39 = mrSges(7,2) * t137 + mrSges(7,3) * t138;
t37 = -mrSges(6,2) * t138 + mrSges(6,3) * t137;
t28 = (-qJ(5) * t222 - qJD(5) * t267) * t301 + t307;
t19 = -pkin(4) * t223 - t21;
t15 = t319 * t222 + (qJD(6) * t298 + (qJ(6) * qJD(4) - qJD(5)) * t301) * t267 + t307;
t9 = -pkin(5) * t316 - t16;
t7 = -pkin(5) * t315 - qJD(6) * t266 - t223 * t425 - t21;
t1 = [(-t52 * mrSges(6,1) - t183 * mrSges(5,2) + t44 * mrSges(7,2) + t70 * mrSges(6,3) - Ifges(5,4) * t443 + Ifges(6,2) * t441 + t457 * t438 - t508 * t440 + t442 * t523 - t515 - t516) * t315 + (-t13 * mrSges(5,3) + t482 + t512) * t390 + m(5) * (t13 * t140 + t139 * t14 + t20 * t79 - t21 * t78 - t396) + m(4) * (t125 * t233 + t144 * t214 + 0.2e1 * t278 * t366 - t396) + t501 * t393 / 0.2e1 + t502 * t394 / 0.2e1 + (mrSges(4,2) * t481 + Ifges(4,1) * t305 - Ifges(4,4) * t205 + t118 * t344 + t119 * t343 + t18 * t330 + t323 * t449 + t327 * t448 + t329 * t6 + t342 * t502 + t345 * t501 + t446 * t487 + t447 * t496 + t450 * t497) * t267 + (t488 / 0.2e1 - t178 / 0.2e1 - Ifges(4,2) * t436 + t458 * t438 - t434 * Ifges(4,4) + t478 * t442 + t477 * t440 + Ifges(6,4) * t441 + Ifges(5,6) * t443 + t40 * mrSges(7,2) - t214 * mrSges(4,3) - t520) * t223 + (Ifges(4,1) * t434 - t213 * mrSges(4,3) + t436 * Ifges(4,4) - t493 / 0.2e1 + t179 / 0.2e1 - t521) * t222 + (mrSges(6,1) * t11 - t480 + t511) * t389 + (-m(4) * t213 + m(5) * t183 + t381) * t145 + (t303 * (-Ifges(3,2) * t300 + t418) - 0.2e1 * pkin(1) * (mrSges(3,1) * t300 + mrSges(3,2) * t303)) * t369 + (Ifges(3,1) * t303 - t419) * t351 + (qJD(1) * (Ifges(3,1) * t300 + t418) + t253) * t374 / 0.2e1 - (t314 + t252) * t375 / 0.2e1 + (-t125 * mrSges(4,3) - t305 * Ifges(4,4) + t455 / 0.2e1 + mrSges(4,1) * t481 + Ifges(6,4) * t449 + Ifges(4,2) * t205 + Ifges(5,6) * t448 + t458 * t446 + t478 * t447 + t477 * t450 + t454) * t266 + (-t276 * t374 - t277 * t375) * pkin(7) + t290 * (t205 * mrSges(4,1) + mrSges(4,2) * t305) + qJD(2) ^ 2 * (Ifges(3,5) * t303 - Ifges(3,6) * t300) / 0.2e1 + (t183 * mrSges(5,1) + t53 * mrSges(6,1) - t70 * mrSges(6,2) + t44 * mrSges(7,3) - Ifges(5,2) * t443 + Ifges(6,6) * t441 - t438 * t513 + t440 * t524 - t507 * t442 + t510) * t316 + t207 * t366 + m(6) * (t11 * t97 + t143 * t18 + t16 * t53 + t19 * t52 + t28 * t70 + t8 * t90) + m(7) * (t106 * t6 + t15 * t44 + t2 * t51 + t36 * t7 + t4 * t61 + t40 * t9) + t331 * t395 + t51 * t55 + t61 * t57 + t90 * t56 + (-t205 * t233 - t305 * t462 + t395) * mrSges(4,3) - t462 * t38 + t97 * t58 + t106 * t39 + t139 * t59 + t140 * t60 + t143 * t37 + t15 * t148 + t28 * t151 + t7 * t158 + t16 * t159 + t9 * t160 + t19 * t161 + t20 * t162 + t21 * t163 + t144 * t236; -m(4) * (-t213 * t216 + t214 * t217) + t304 + (qJD(4) * t515 - t379 * t367 - t482) * t301 + 0.2e1 * ((t125 * t299 - t126 * t302) * t453 + ((-t213 * t299 + t214 * t302) * t453 + m(6) * (-t384 * t53 + t386 * t52) / 0.2e1 + m(5) * (t183 * t299 + t384 * t79 + t386 * t78) / 0.2e1) * qJD(3)) * pkin(2) + ((t401 / 0.2e1 - t253 / 0.2e1 - t291 / 0.2e1 + pkin(1) * mrSges(3,2) * qJD(1) - pkin(2) * t368 * mrSges(4,3) * t302 ^ 2 + (t276 - t398) * pkin(7)) * t303 + (-t400 / 0.2e1 + t252 / 0.2e1 + (pkin(1) * mrSges(3,1) + t419 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t303) * qJD(1) + (t277 + t397) * pkin(7) + (t302 * t368 * t422 - m(4) * t278 - t207) * pkin(2)) * t300) * qJD(1) + (-t205 * t422 + (t236 * t302 + t299 * t381) * qJD(3)) * pkin(2) - m(5) * (-t104 * t78 + t105 * t79 + t183 * t216) - m(6) * (t128 * t70 + t52 * t63 + t53 * t62) + ((qJD(4) * t379 + t423) * t298 + (qJD(4) * t378 + t424) * t301 + t461) * t288 + (-t128 + t238) * t151 + (t510 * qJD(4) + t378 * t367 - t480) * t298 + m(6) * (t18 * t263 + t238 * t70) + t467 * t148 + t470 * t160 + t471 * (-pkin(3) - t429) + t473 * t158 + (t2 * t264 + t242 * t6 + t265 * t4 + t36 * t473 + t40 * t470 + t44 * t467) * m(7) - t381 * t216 + t214 * t402 - t62 * t159 - t63 * t161 - t105 * t162 - t104 * t163 - t40 * t193 - t217 * t236 + t242 * t39 + t263 * t37 + t264 * t55 + t265 * t57; -mrSges(6,1) * t426 + t304 + (t424 * t301 + t423 * t298 + (t298 * t379 + t301 * t378) * qJD(4) + t461) * pkin(9) + t469 * t160 - m(5) * (-t112 * t78 + t113 * t79 + t183 * t214) + t472 * t158 + t463 * t151 + t468 * t148 + (-t406 + (-t298 * t79 + t301 * t78) * qJD(4)) * mrSges(5,3) + (-t381 + t402) * t214 + (-mrSges(7,1) * t371 - t193) * t40 - t72 * t159 - t73 * t161 - t113 * t162 - t112 * t163 - t213 * t236 + t260 * t39 + t275 * t37 + t279 * t55 + t281 * t57 - t471 * pkin(3) + (t2 * t279 + t260 * t6 + t281 * t4 + t36 * t472 + t40 * t469 + t44 * t468) * m(7) + (t18 * t275 + t463 * t70 - t52 * t73 - t53 * t72) * m(6); (-Ifges(5,2) * t231 - t228 + t501) * t442 + t455 + t454 + (t230 * t36 + t231 * t40) * mrSges(7,1) + (-t230 * t508 + t226 - t405 + t502) * t441 + (t230 * t52 - t231 * t53) * mrSges(6,1) + (-t378 + t420) * t79 + (-t379 + t421) * t78 + (-t56 + t57) * qJ(5) + (t230 * t457 - t231 * t513) * t439 + (-t231 * t507 + t118 + t227 - t412) * t443 + (-pkin(4) * t11 - qJ(5) * t8 - t149 * t70 + t465 * t53 - t52 * t79) * m(6) + (qJ(5) * t4 - t2 * t425 + t464 * t36 + t40 * t491 - t44 * t87) * m(7) - t380 * qJD(5) + (Ifges(6,2) * t230 + t119 + t404) * t440 - pkin(4) * t58 + t464 * t158 - t87 * t148 - t149 * t151 + t317 * t160 - t70 * (-mrSges(6,2) * t231 + mrSges(6,3) * t230) - t44 * (mrSges(7,2) * t230 + mrSges(7,3) * t231) - t183 * (mrSges(5,1) * t231 - mrSges(5,2) * t230) - t425 * t55; t380 * t250 + (t148 + t151) * t231 + t55 + t58 + (t231 * t44 - t250 * t40 + t2) * m(7) + (t231 * t70 + t250 * t53 + t11) * m(6); -t230 * t148 + t250 * t158 + 0.2e1 * (t4 / 0.2e1 + t44 * t443 + t36 * t438) * m(7) + t57;];
tauc  = t1(:);
