% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:06:54
% EndTime: 2018-11-23 17:07:17
% DurationCPUTime: 24.34s
% Computational Cost: add. (14012->739), mult. (37442->996), div. (0->0), fcn. (29854->10), ass. (0->359)
t285 = cos(qJ(2));
t280 = cos(pkin(6));
t395 = qJD(1) * t280;
t388 = pkin(1) * t395;
t366 = t285 * t388;
t283 = sin(qJ(2));
t279 = sin(pkin(6));
t396 = qJD(1) * t279;
t383 = t283 * t396;
t234 = -pkin(8) * t383 + t366;
t365 = qJD(2) + t395;
t197 = -pkin(2) * t365 + qJD(3) - t234;
t278 = sin(pkin(11));
t408 = cos(pkin(11));
t295 = t278 * t383 - t365 * t408;
t160 = pkin(3) * t295 + t197;
t382 = t285 * t396;
t262 = -qJD(4) + t382;
t431 = -t262 / 0.2e1;
t235 = pkin(8) * t382 + t283 * t388;
t206 = qJ(3) * t365 + t235;
t228 = (-pkin(2) * t285 - qJ(3) * t283 - pkin(1)) * t279;
t213 = qJD(1) * t228;
t146 = t408 * t206 + t278 * t213;
t114 = -pkin(9) * t295 + t146;
t282 = sin(qJ(4));
t145 = -t206 * t278 + t408 * t213;
t370 = t279 * t408;
t351 = qJD(1) * t370;
t220 = t278 * t365 + t283 * t351;
t302 = -pkin(3) * t382 - pkin(9) * t220 + t145;
t428 = cos(qJ(4));
t98 = t428 * t302;
t53 = t114 * t282 - t98;
t491 = -qJD(5) - t53;
t47 = pkin(4) * t262 - t491;
t281 = sin(qJ(6));
t284 = cos(qJ(6));
t452 = pkin(4) + pkin(10);
t289 = t220 * t428 - t282 * t295;
t330 = pkin(5) * t289 + t53;
t511 = qJD(5) + t330;
t31 = t262 * t452 + t511;
t211 = t428 * t295;
t163 = t220 * t282 + t211;
t287 = -qJ(5) * t289 + t160;
t41 = t163 * t452 + t287;
t5 = -t281 * t41 + t284 * t31;
t438 = -t289 / 0.2e1;
t521 = Ifges(6,2) * t438;
t440 = -t163 / 0.2e1;
t522 = Ifges(6,6) * t440;
t6 = t281 * t31 + t284 * t41;
t65 = t163 * pkin(4) + t287;
t515 = t47 * mrSges(6,1) + t5 * mrSges(7,1) + t160 * mrSges(5,2) - t6 * mrSges(7,2) + t53 * mrSges(5,3) - t65 * mrSges(6,3) - Ifges(6,4) * t431 - t521 + t522;
t524 = t515 + t522;
t352 = t428 * t408;
t400 = t279 * t285;
t323 = t352 * t400;
t312 = qJD(1) * t323;
t361 = t278 * t382;
t203 = -t282 * t361 + t312;
t391 = qJD(4) * t282;
t239 = -qJD(4) * t352 + t278 * t391;
t523 = t239 + t203;
t247 = t278 * t428 + t282 * t408;
t303 = t247 * t400;
t202 = qJD(1) * t303;
t240 = t247 * qJD(4);
t517 = t240 - t202;
t342 = pkin(2) * t283 - qJ(3) * t285;
t233 = t342 * t396;
t172 = t408 * t233 - t278 * t234;
t369 = t285 * t408;
t309 = t279 * (pkin(3) * t283 - pkin(9) * t369);
t142 = qJD(1) * t309 + t172;
t173 = t278 * t233 + t408 * t234;
t154 = -pkin(9) * t361 + t173;
t255 = (-pkin(9) - qJ(3)) * t278;
t368 = t408 * qJ(3);
t256 = pkin(9) * t408 + t368;
t201 = t282 * t255 + t256 * t428;
t493 = qJD(3) * t247 + qJD(4) * t201 + t142 * t428 - t282 * t154;
t195 = pkin(3) * t361 + t235;
t520 = qJ(5) * t523 - qJD(5) * t247 - t195;
t389 = qJD(1) * qJD(2);
t373 = t285 * t389;
t358 = t279 * t373;
t337 = t278 * t358;
t112 = qJD(4) * t211 - qJD(2) * t312 + (qJD(4) * t220 + t337) * t282;
t451 = -t112 / 0.2e1;
t297 = qJD(2) * t303;
t113 = qJD(1) * t297 + qJD(4) * t289;
t449 = -t113 / 0.2e1;
t430 = t262 / 0.2e1;
t437 = t289 / 0.2e1;
t502 = Ifges(6,1) + Ifges(5,3);
t501 = Ifges(6,4) - Ifges(5,5);
t500 = Ifges(6,5) - Ifges(5,6);
t519 = t452 * t517 + t520;
t401 = t279 * t283;
t364 = t452 * t401;
t518 = -pkin(5) * t523 + qJD(1) * t364 + t493;
t375 = qJD(4) * t428;
t392 = qJD(3) * t278;
t496 = qJD(3) * t352 - t428 * t154 + t255 * t375 + (-qJD(4) * t256 - t142 - t392) * t282;
t448 = t113 / 0.2e1;
t516 = t448 - t449;
t296 = t282 * t302;
t54 = t114 * t428 + t296;
t48 = t262 * qJ(5) - t54;
t509 = -t160 * mrSges(5,1) + t65 * mrSges(6,2) + Ifges(6,5) * t430 + Ifges(5,6) * t431 + (Ifges(5,2) + Ifges(6,3)) * t440 + (Ifges(5,4) + Ifges(6,6)) * t437 - t48 * mrSges(6,1) + t54 * mrSges(5,3);
t394 = qJD(2) * t279;
t371 = t394 / 0.2e1;
t159 = qJD(6) + t289;
t410 = Ifges(7,3) * t159;
t131 = t163 * t284 + t262 * t281;
t411 = Ifges(7,6) * t131;
t132 = t163 * t281 - t262 * t284;
t412 = Ifges(7,5) * t132;
t43 = t410 + t411 + t412;
t158 = Ifges(5,4) * t163;
t85 = Ifges(5,1) * t289 - Ifges(5,5) * t262 - t158;
t514 = t85 + t43;
t495 = qJ(5) * t383 - t496;
t439 = t163 / 0.2e1;
t510 = -Ifges(5,4) * t437 - Ifges(5,2) * t440 + Ifges(6,6) * t438 + Ifges(6,3) * t439 + t431 * t500 - t509;
t441 = t159 / 0.2e1;
t443 = t132 / 0.2e1;
t445 = t131 / 0.2e1;
t508 = -Ifges(5,1) * t437 - Ifges(5,4) * t440 - Ifges(7,5) * t443 + Ifges(6,6) * t439 - Ifges(7,6) * t445 - Ifges(7,3) * t441 + t431 * t501 - t515 + t521;
t427 = pkin(1) * t280;
t242 = pkin(8) * t400 + t283 * t427;
t226 = t242 * t389;
t186 = pkin(3) * t337 + t226;
t210 = (qJD(2) * t342 - qJD(3) * t283) * t279;
t191 = qJD(1) * t210;
t374 = t279 * t389;
t359 = t283 * t374;
t225 = -pkin(8) * t359 + qJD(2) * t366;
t192 = qJD(3) * t365 + t225;
t134 = t408 * t191 - t278 * t192;
t307 = qJD(2) * t309;
t100 = qJD(1) * t307 + t134;
t135 = t278 * t191 + t408 * t192;
t115 = -pkin(9) * t337 + t135;
t322 = qJD(4) * t296 - t100 * t428 + t114 * t375 + t282 * t115;
t19 = -pkin(4) * t359 + t322;
t291 = qJ(5) * t112 - qJD(5) * t289 + t186;
t32 = pkin(4) * t113 + t291;
t353 = qJD(1) * t371;
t334 = t283 * t353;
t450 = t112 / 0.2e1;
t61 = -qJD(6) * t132 + t113 * t284 - t281 * t359;
t458 = t61 / 0.2e1;
t60 = qJD(6) * t131 + t113 * t281 + t284 * t359;
t459 = t60 / 0.2e1;
t23 = t113 * t452 + t291;
t338 = qJD(2) * t364;
t8 = -t112 * pkin(5) - qJD(1) * t338 + t322;
t1 = qJD(6) * t5 + t23 * t284 + t281 * t8;
t2 = -qJD(6) * t6 - t23 * t281 + t284 * t8;
t481 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t504 = -t359 / 0.2e1;
t9 = Ifges(7,5) * t60 + Ifges(7,6) * t61 - Ifges(7,3) * t112;
t507 = -t501 * t334 + t19 * mrSges(6,1) + t186 * mrSges(5,2) + t322 * mrSges(5,3) - t32 * mrSges(6,3) + 0.2e1 * Ifges(5,4) * t449 + Ifges(7,5) * t459 + Ifges(7,6) * t458 + Ifges(6,4) * t504 + t481 - t516 * Ifges(6,6) - Ifges(6,2) * t450 + Ifges(5,5) * t359 / 0.2e1 + t9 / 0.2e1 + (0.2e1 * Ifges(5,1) + Ifges(7,3) + Ifges(6,2)) * t451;
t506 = Ifges(6,5) / 0.2e1;
t416 = Ifges(4,4) * t278;
t300 = Ifges(4,5) * t283 + (Ifges(4,1) * t408 - t416) * t285;
t505 = t300 * t371;
t246 = t278 * t282 - t352;
t277 = -pkin(3) * t408 - pkin(2);
t317 = -t247 * qJ(5) + t277;
t166 = t246 * t452 + t317;
t200 = -t428 * t255 + t256 * t282;
t174 = pkin(5) * t247 + t200;
t91 = t166 * t284 + t174 * t281;
t499 = -qJD(6) * t91 - t281 * t519 + t284 * t518;
t90 = -t166 * t281 + t174 * t284;
t498 = qJD(6) * t90 + t281 * t518 + t284 * t519;
t497 = t225 * mrSges(3,2);
t494 = pkin(4) * t383 + t493;
t492 = -pkin(5) * t517 - t495;
t490 = pkin(4) * t517 + t520;
t176 = t202 * t284 - t281 * t383;
t390 = qJD(6) * t246;
t379 = t281 * t390;
t404 = t240 * t284;
t326 = t379 - t404;
t489 = t176 + t326;
t177 = t202 * t281 + t284 * t383;
t405 = t240 * t281;
t327 = t284 * t390 + t405;
t488 = t177 - t327;
t487 = -t134 * t278 + t408 * t135;
t486 = t112 * t501 + t113 * t500 + t359 * t502;
t485 = t1 * t281 + t2 * t284;
t393 = qJD(2) * t283;
t381 = t279 * t393;
t387 = qJD(2) * t427;
t236 = -pkin(8) * t381 + t285 * t387;
t219 = qJD(3) * t280 + t236;
t152 = t408 * t210 - t278 * t219;
t119 = t152 + t307;
t227 = qJ(3) * t280 + t242;
t168 = -t227 * t278 + t408 * t228;
t238 = t280 * t278 + t283 * t370;
t121 = -pkin(3) * t400 - pkin(9) * t238 + t168;
t153 = t278 * t210 + t408 * t219;
t380 = t285 * t394;
t360 = t278 * t380;
t133 = -pkin(9) * t360 + t153;
t169 = t408 * t227 + t278 * t228;
t386 = t278 * t401;
t313 = -t280 * t408 + t386;
t137 = -pkin(9) * t313 + t169;
t318 = -t282 * t119 - t121 * t375 - t428 * t133 + t137 * t391;
t24 = -t279 * (qJ(5) * t393 - qJD(5) * t285) + t318;
t376 = Ifges(4,4) * t408;
t299 = Ifges(4,6) * t283 + (-Ifges(4,2) * t278 + t376) * t285;
t483 = t197 * (mrSges(4,1) * t278 + mrSges(4,2) * t408) * t400 + (-t299 * t295 / 0.2e1 + t365 * (Ifges(3,5) * t285 - Ifges(3,6) * t283) / 0.2e1) * t279;
t482 = t283 * (Ifges(4,5) * t220 - Ifges(4,6) * t295 - Ifges(4,3) * t382 + t163 * t500 - t262 * t502 - t289 * t501);
t363 = mrSges(3,3) * t383;
t480 = -m(4) * t197 + mrSges(3,1) * t365 - mrSges(4,1) * t295 - t220 * mrSges(4,2) - t363;
t20 = qJD(4) * t98 + t282 * t100 - t114 * t391 + t428 * t115;
t16 = -qJ(5) * t359 + qJD(5) * t262 - t20;
t476 = -mrSges(5,1) * t322 - t20 * mrSges(5,2) + t19 * mrSges(6,2) - t16 * mrSges(6,3);
t475 = -t53 * mrSges(5,1) - t54 * mrSges(5,2) + t47 * mrSges(6,2) - t48 * mrSges(6,3);
t343 = Ifges(7,5) * t281 + Ifges(7,6) * t284;
t414 = Ifges(7,4) * t281;
t344 = Ifges(7,2) * t284 + t414;
t413 = Ifges(7,4) * t284;
t345 = Ifges(7,1) * t281 + t413;
t346 = mrSges(7,1) * t284 - mrSges(7,2) * t281;
t348 = t281 * t5 - t284 * t6;
t425 = t163 * pkin(5);
t35 = -t48 - t425;
t429 = -t281 / 0.2e1;
t442 = -t159 / 0.2e1;
t444 = -t132 / 0.2e1;
t446 = -t131 / 0.2e1;
t125 = Ifges(7,4) * t131;
t45 = Ifges(7,1) * t132 + Ifges(7,5) * t159 + t125;
t415 = Ifges(7,4) * t132;
t44 = Ifges(7,2) * t131 + Ifges(7,6) * t159 + t415;
t465 = -t44 / 0.2e1;
t473 = mrSges(7,3) * t348 + t284 * t465 + t343 * t442 + t344 * t446 + t345 * t444 + t35 * t346 + t429 * t45;
t472 = t186 * mrSges(5,1) + t16 * mrSges(6,1) - t32 * mrSges(6,2) - t20 * mrSges(5,3) + Ifges(5,6) * t504 + 0.2e1 * Ifges(6,3) * t448 + t359 * t506 + (Ifges(5,4) + 0.2e1 * Ifges(6,6)) * t450 + t516 * Ifges(5,2);
t468 = t279 ^ 2;
t466 = Ifges(7,1) * t459 + Ifges(7,4) * t458 + Ifges(7,5) * t451;
t464 = t44 / 0.2e1;
t463 = t45 / 0.2e1;
t435 = qJD(1) * t505;
t423 = t285 * pkin(1);
t418 = Ifges(3,4) * t283;
t417 = Ifges(3,4) * t285;
t138 = mrSges(6,1) * t163 + mrSges(6,3) * t262;
t70 = -mrSges(7,1) * t131 + mrSges(7,2) * t132;
t409 = -t138 + t70;
t407 = qJ(5) * t163;
t403 = t246 * t281;
t402 = t246 * t284;
t69 = t282 * t121 + t428 * t137;
t140 = mrSges(5,2) * t262 - mrSges(5,3) * t163;
t399 = -t138 + t140;
t139 = mrSges(6,1) * t289 - mrSges(6,2) * t262;
t141 = -mrSges(5,1) * t262 - mrSges(5,3) * t289;
t398 = -t139 + t141;
t332 = t369 * t394;
t207 = qJD(1) * mrSges(4,2) * t332 + mrSges(4,1) * t337;
t237 = pkin(8) * t380 + t283 * t387;
t196 = pkin(3) * t360 + t237;
t372 = -t396 / 0.2e1;
t93 = -t112 * mrSges(6,1) + mrSges(6,2) * t359;
t362 = mrSges(3,3) * t382;
t356 = t285 * t372;
t68 = t121 * t428 - t282 * t137;
t33 = -mrSges(7,1) * t112 - mrSges(7,3) * t60;
t34 = mrSges(7,2) * t112 + mrSges(7,3) * t61;
t341 = t281 * t34 + t284 * t33;
t179 = t238 * t428 - t282 * t313;
t67 = pkin(4) * t400 - t68;
t42 = t179 * pkin(5) + pkin(10) * t400 + t67;
t304 = t428 * t313;
t178 = t238 * t282 + t304;
t272 = pkin(8) * t401;
t180 = pkin(3) * t386 + t272 + (t277 - t423) * t280;
t290 = -t179 * qJ(5) + t180;
t63 = t178 * t452 + t290;
t17 = -t281 * t63 + t284 * t42;
t18 = t281 * t42 + t284 * t63;
t77 = -mrSges(7,2) * t159 + mrSges(7,3) * t131;
t78 = mrSges(7,1) * t159 - mrSges(7,3) * t132;
t340 = -t281 * t78 + t284 * t77;
t339 = -t281 * t77 - t284 * t78;
t336 = -t360 / 0.2e1;
t66 = qJ(5) * t400 - t69;
t329 = -mrSges(4,3) * t278 * t285 - mrSges(4,2) * t283;
t328 = -t178 * t281 + t284 * t400;
t155 = t178 * t284 + t281 * t400;
t325 = qJD(1) * t336;
t321 = t332 / 0.2e1;
t320 = -t119 * t428 + t121 * t391 + t282 * t133 + t137 * t375;
t319 = Ifges(3,5) * t358 - Ifges(3,6) * t359;
t316 = pkin(1) * t468 * (mrSges(3,1) * t283 + mrSges(3,2) * t285);
t315 = mrSges(4,1) * t283 - mrSges(4,3) * t369;
t314 = t283 * t468 * (Ifges(3,1) * t285 - t418);
t311 = qJD(1) * t321;
t126 = qJD(4) * t304 - qJD(2) * t323 + (qJD(4) * t238 + t360) * t282;
t310 = qJ(5) * t126 - qJD(5) * t179 + t196;
t308 = t329 * t394;
t305 = t315 * t394;
t298 = (Ifges(3,6) * t280 + (t285 * Ifges(3,2) + t418) * t279) * qJD(1);
t293 = t285 * t468 * (Ifges(4,3) * t283 + (Ifges(4,5) * t408 - Ifges(4,6) * t278) * t285);
t263 = Ifges(3,4) * t382;
t241 = t280 * t423 - t272;
t232 = -mrSges(3,2) * t365 + t362;
t230 = t272 + (-pkin(2) - t423) * t280;
t215 = qJD(1) * t305;
t214 = qJD(1) * t308;
t205 = Ifges(3,1) * t383 + Ifges(3,5) * t365 + t263;
t204 = Ifges(3,6) * qJD(2) + t298;
t189 = -mrSges(4,1) * t382 - mrSges(4,3) * t220;
t188 = mrSges(4,2) * t382 - mrSges(4,3) * t295;
t183 = t246 * pkin(4) + t317;
t181 = t299 * t374;
t175 = -t246 * pkin(5) + t201;
t151 = Ifges(4,1) * t220 - Ifges(4,4) * t295 - Ifges(4,5) * t382;
t150 = Ifges(4,4) * t220 - Ifges(4,2) * t295 - Ifges(4,6) * t382;
t127 = qJD(4) * t179 + t297;
t104 = t112 * mrSges(5,2);
t103 = t112 * mrSges(6,3);
t95 = -mrSges(5,2) * t359 - mrSges(5,3) * t113;
t94 = mrSges(5,1) * t359 + mrSges(5,3) * t112;
t92 = mrSges(6,1) * t113 - mrSges(6,3) * t359;
t89 = -mrSges(6,2) * t163 - mrSges(6,3) * t289;
t88 = mrSges(5,1) * t163 + mrSges(5,2) * t289;
t87 = pkin(4) * t289 + t407;
t79 = t178 * pkin(4) + t290;
t72 = qJD(6) * t328 + t127 * t284 - t281 * t381;
t71 = qJD(6) * t155 + t127 * t281 + t284 * t381;
t64 = t289 * t452 + t407;
t59 = t113 * mrSges(5,1) - t104;
t58 = -t113 * mrSges(6,2) + t103;
t46 = -pkin(5) * t178 - t66;
t40 = pkin(4) * t127 + t310;
t39 = t54 - t425;
t30 = t127 * t452 + t310;
t25 = -pkin(4) * t381 + t320;
t22 = -mrSges(7,1) * t61 + mrSges(7,2) * t60;
t15 = -pkin(5) * t127 - t24;
t14 = -t126 * pkin(5) + t320 - t338;
t13 = t281 * t39 + t284 * t64;
t12 = -t281 * t64 + t284 * t39;
t10 = t60 * Ifges(7,4) + t61 * Ifges(7,2) - t112 * Ifges(7,6);
t7 = -pkin(5) * t113 - t16;
t4 = -qJD(6) * t18 + t14 * t284 - t281 * t30;
t3 = qJD(6) * t17 + t14 * t281 + t284 * t30;
t11 = [(-m(3) * t234 - t480) * t237 + t483 * qJD(2) + (-Ifges(5,4) * t451 + t472) * t178 - (t298 + t204) * t381 / 0.2e1 + (Ifges(7,5) * t71 + Ifges(7,6) * t72) * t441 + (-t134 * mrSges(4,1) + t135 * mrSges(4,2) - Ifges(6,4) * t450 - Ifges(4,5) * t311 - Ifges(5,5) * t451 - Ifges(6,5) * t448 - Ifges(4,6) * t325 - Ifges(5,6) * t449 - t476 - t486 / 0.2e1 + t225 * mrSges(3,3) + (-Ifges(4,3) - t502) * t334) * t400 + m(4) * (t134 * t168 + t135 * t169 + t145 * t152 + t146 * t153) + m(3) * (t225 * t242 + t235 * t236) + (Ifges(7,1) * t71 + Ifges(7,4) * t72) * t443 + (-t497 + t319 / 0.2e1) * t280 + m(6) * (t16 * t66 + t19 * t67 + t24 * t48 + t25 * t47 + t32 * t79 + t40 * t65) + m(7) * (t1 * t18 + t15 * t35 + t17 * t2 + t3 * t6 + t4 * t5 + t46 * t7) + (Ifges(7,4) * t71 + Ifges(7,2) * t72) * t445 + (-0.2e1 * t316 - t293 + t314) * t389 + (Ifges(4,1) * t238 - Ifges(4,4) * t313) * t311 + (Ifges(4,4) * t238 - Ifges(4,2) * t313) * t325 - t318 * t140 + m(5) * (t160 * t196 + t180 * t186 + t20 * t69 - t318 * t54 + t320 * t53 - t322 * t68) - t320 * t141 + (-Ifges(7,5) * t328 + Ifges(7,6) * t155) * t451 + (-Ifges(7,1) * t328 + Ifges(7,4) * t155) * t459 + (-Ifges(7,4) * t328 + Ifges(7,2) * t155) * t458 + (t1 * t155 + t2 * t328 - t5 * t71 + t6 * t72) * mrSges(7,3) + t7 * (-mrSges(7,1) * t155 - mrSges(7,2) * t328) - t328 * t466 + (t205 * t371 + (Ifges(3,5) * t280 + (Ifges(3,1) * t283 + t417) * t279) * t353) * t285 + (-t514 / 0.2e1 + t508) * t126 + t510 * t127 + (-m(3) * t241 + m(4) * t230 - mrSges(3,1) * t280 + mrSges(4,1) * t313 + t238 * mrSges(4,2)) * t226 + t146 * t308 + t145 * t305 + (-t134 * t238 - t135 * t313) * mrSges(4,3) + t507 * t179 + (t226 * t401 - t234 * t380 - t235 * t381 - t241 * t358 - t242 * t359) * mrSges(3,3) - t313 * t181 / 0.2e1 + t71 * t463 + t72 * t464 + t238 * t435 + t468 * (-Ifges(3,2) * t283 + t417) * t373 + t236 * t232 + t230 * t207 + t169 * t214 + t168 * t215 + t153 * t188 + t152 * t189 + t196 * t88 + t180 * t59 + t151 * t321 + t150 * t336 + (Ifges(6,4) * t438 + Ifges(5,5) * t437 + Ifges(6,5) * t439 + Ifges(5,6) * t440 + t431 * t502 + t475) * t381 + (Ifges(4,5) * t238 - Ifges(4,6) * t313 + t500 * t178) * t334 + t371 * t482 + t220 * t505 + t17 * t33 + t18 * t34 + t46 * t22 + t15 * t70 + t35 * (-mrSges(7,1) * t72 + mrSges(7,2) * t71) + t3 * t77 + t4 * t78 + t79 * t58 + t40 * t89 + t66 * t92 + t67 * t93 + t68 * t94 + t69 * t95 + t24 * t138 + t25 * t139 + t155 * t10 / 0.2e1; (t363 + t480) * t235 + t492 * t70 - t493 * t141 + t494 * t139 + (-t16 * t201 + t183 * t32 + t19 * t200 + t47 * t494 + t48 * t495 + t490 * t65) * m(6) + t495 * t138 + t496 * t140 + (t95 - t92) * t201 + (-mrSges(4,1) * t408 + t278 * mrSges(4,2) - mrSges(3,1)) * t226 + (qJD(6) * t45 + t10) * t402 / 0.2e1 + (t220 * t300 + t482) * t372 + (-t392 - t172) * t189 + (Ifges(7,4) * t327 - Ifges(7,2) * t326) * t445 + (Ifges(7,4) * t177 + Ifges(7,2) * t176) * t446 + t514 * (-t239 / 0.2e1 - t203 / 0.2e1) + (t379 + t176) * t465 + (Ifges(7,5) * t327 - Ifges(7,6) * t326) * t441 + (Ifges(7,5) * t177 + Ifges(7,6) * t176) * t442 + (t205 + t263) * t356 + t376 * t311 + t416 * t325 + (Ifges(4,2) * t325 + t181 / 0.2e1 + Ifges(4,6) * t334) * t408 + t319 + (-t145 * t315 - t146 * t329) * t396 + (-t160 * t195 + t186 * t277 + t20 * t201 + t200 * t322 + t493 * t53 + t496 * t54) * m(5) + (qJD(3) * t408 - t173) * t188 + t508 * t239 + (-Ifges(5,4) * t438 - Ifges(5,2) * t439 + Ifges(6,6) * t437 + Ifges(6,3) * t440 + t430 * t500 + t509) * t202 + t510 * t240 + (-t483 + (t293 / 0.2e1 + t316 - t314 / 0.2e1) * qJD(1)) * qJD(1) + (t93 - t94) * t200 + (-t145 * t172 - t146 * t173 - t226 * pkin(2) + (-t145 * t278 + t146 * t408) * qJD(3) + t487 * qJ(3)) * m(4) + t487 * mrSges(4,3) + (mrSges(7,1) * t489 - mrSges(7,2) * t488) * t35 + (t1 * t402 - t2 * t403 + t488 * t5 - t489 * t6) * mrSges(7,3) + t490 * t89 + t507 * t247 + (t362 - t232) * t234 + (Ifges(5,1) * t438 + Ifges(5,4) * t439 + Ifges(7,5) * t444 - Ifges(6,2) * t437 + Ifges(7,6) * t446 + Ifges(7,3) * t442 - t430 * t501 - t524) * t203 + (Ifges(7,1) * t327 - Ifges(7,4) * t326) * t443 + (Ifges(7,1) * t177 + Ifges(7,4) * t176) * t444 + t405 * t463 + t404 * t464 + t403 * t466 + t277 * t59 + t150 * t361 / 0.2e1 - pkin(2) * t207 - t195 * t88 + t183 * t58 - t177 * t45 / 0.2e1 + t175 * t22 + (Ifges(4,1) * t311 + Ifges(4,5) * t334 - qJ(3) * t215 + t435) * t278 + t498 * t77 + t499 * t78 + (t1 * t91 + t175 * t7 + t2 * t90 + t35 * t492 + t498 * t6 + t499 * t5) * m(7) + (t344 * t458 + t345 * t459 - t346 * t7 + (-Ifges(5,4) + t343) * t451 + t500 * t334 + t472) * t246 + (t204 / 0.2e1 - Ifges(3,2) * t356 + Ifges(6,4) * t437 + Ifges(5,5) * t438 + Ifges(5,6) * t439 + Ifges(6,5) * t440 + t502 * t430 - t475) * t383 - t497 - t285 * t151 * t351 / 0.2e1 + t90 * t33 + t91 * t34 + t214 * t368; t295 * t188 - t104 + t103 + t284 * t34 - t281 * t33 + t220 * t189 + (-mrSges(6,2) + mrSges(5,1)) * t113 + t339 * qJD(6) - (-t70 - t399) * t163 + (t339 + t398) * t289 + t207 + (t1 * t284 + t163 * t35 - t2 * t281 - t159 * (t281 * t6 + t284 * t5)) * m(7) + (-t163 * t48 - t289 * t47 + t32) * m(6) + (t163 * t54 - t289 * t53 + t186) * m(5) + (t145 * t220 + t146 * t295 + t226) * m(4); (-t158 / 0.2e1 + t411 / 0.2e1 + t412 / 0.2e1 + t410 / 0.2e1 + t43 / 0.2e1 + t85 / 0.2e1 + (-Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t262 + t524) * t163 + ((Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t289 + (-Ifges(5,6) / 0.2e1 + t506) * t262 + (-Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(6,3) / 0.2e1) * t163 + t473 + t509) * t289 + t473 * qJD(6) + t476 + (t7 * qJ(5) - t12 * t5 - t13 * t6 + t511 * t35) * m(7) + t486 - t485 * mrSges(7,3) - (t341 + (-m(7) * t348 + t340) * qJD(6) + m(7) * t485) * t452 + (t22 - t92) * qJ(5) + t409 * qJD(5) + t398 * t54 + t399 * t53 + t284 * t466 + (Ifges(7,5) * t284 - Ifges(7,6) * t281) * t451 + (-Ifges(7,2) * t281 + t413) * t458 + (Ifges(7,1) * t284 - t414) * t459 + t10 * t429 + t7 * (mrSges(7,1) * t281 + mrSges(7,2) * t284) + (-pkin(4) * t19 - qJ(5) * t16 - t47 * t54 + t48 * t491 - t65 * t87) * m(6) + t330 * t70 - t13 * t77 - t12 * t78 - t87 * t89 - pkin(4) * t93; t409 * t262 + t340 * qJD(6) + (t340 + t89) * t289 + t341 + t93 + (-t159 * t348 + t262 * t35 + t485) * m(7) + (-t262 * t48 + t289 * t65 + t19) * m(6); -t35 * (mrSges(7,1) * t132 + mrSges(7,2) * t131) + (Ifges(7,1) * t131 - t415) * t444 + t44 * t443 + (Ifges(7,5) * t131 - Ifges(7,6) * t132) * t442 - t5 * t77 + t6 * t78 + (t131 * t5 + t132 * t6) * mrSges(7,3) + t9 + (-Ifges(7,2) * t132 + t125 + t45) * t446 + t481;];
tauc  = t11(:);
