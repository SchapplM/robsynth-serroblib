% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:23
% EndTime: 2019-03-09 07:01:13
% DurationCPUTime: 29.53s
% Computational Cost: add. (16461->868), mult. (34778->1166), div. (0->0), fcn. (23851->18), ass. (0->383)
t313 = sin(pkin(11));
t291 = pkin(1) * t313 + pkin(7);
t274 = t291 * qJD(1);
t318 = sin(qJ(3));
t323 = cos(qJ(3));
t216 = qJD(2) * t323 - t318 * t274;
t365 = pkin(3) * t318 - pkin(8) * t323;
t262 = t365 * qJD(1);
t317 = sin(qJ(4));
t322 = cos(qJ(4));
t155 = -t216 * t317 + t322 * t262;
t411 = t322 * t323;
t346 = pkin(4) * t318 - pkin(9) * t411;
t325 = -pkin(9) - pkin(8);
t383 = qJD(4) * t325;
t557 = -qJD(1) * t346 + t322 * t383 - t155;
t156 = t322 * t216 + t317 * t262;
t404 = qJD(1) * t323;
t382 = t317 * t404;
t556 = -pkin(9) * t382 - t317 * t383 + t156;
t316 = sin(qJ(5));
t321 = cos(qJ(5));
t348 = t316 * t317 - t321 * t322;
t512 = qJD(4) + qJD(5);
t177 = t512 * t348;
t340 = t348 * t323;
t209 = qJD(1) * t340;
t555 = t177 - t209;
t256 = t316 * t322 + t317 * t321;
t178 = t512 * t256;
t341 = t256 * t323;
t208 = qJD(1) * t341;
t551 = t178 - t208;
t280 = t325 * t317;
t281 = t325 * t322;
t396 = qJD(5) * t321;
t397 = qJD(5) * t316;
t529 = t280 * t396 + t281 * t397 + t316 * t557 - t556 * t321;
t191 = t316 * t280 - t321 * t281;
t528 = -qJD(5) * t191 + t556 * t316 + t321 * t557;
t554 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t405 = qJD(1) * t318;
t553 = -pkin(5) * t405 + pkin(10) * t555 + t528;
t552 = pkin(10) * t551 - t529;
t402 = qJD(3) * t322;
t253 = -t317 * t405 + t402;
t254 = qJD(3) * t317 + t322 * t405;
t172 = t253 * t316 + t254 * t321;
t315 = sin(qJ(6));
t320 = cos(qJ(6));
t367 = t321 * t253 - t254 * t316;
t105 = t172 * t320 + t315 * t367;
t205 = -qJD(3) * pkin(3) - t216;
t165 = -pkin(4) * t253 + t205;
t107 = -pkin(5) * t367 + t165;
t393 = qJD(1) * qJD(3);
t266 = qJDD(1) * t323 - t318 * t393;
t248 = qJDD(4) - t266;
t241 = qJDD(5) + t248;
t229 = qJDD(6) + t241;
t368 = -t172 * t315 + t320 * t367;
t267 = qJDD(1) * t318 + t323 * t393;
t166 = qJD(4) * t253 + qJDD(3) * t317 + t267 * t322;
t167 = -qJD(4) * t254 + qJDD(3) * t322 - t267 * t317;
t68 = qJD(5) * t367 + t166 * t321 + t167 * t316;
t69 = -qJD(5) * t172 - t166 * t316 + t167 * t321;
t26 = qJD(6) * t368 + t315 * t69 + t320 * t68;
t27 = -qJD(6) * t105 - t315 * t68 + t320 * t69;
t391 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t229;
t286 = qJD(4) - t404;
t282 = qJD(5) + t286;
t543 = pkin(10) * t172;
t217 = t318 * qJD(2) + t323 * t274;
t206 = qJD(3) * pkin(8) + t217;
t366 = pkin(3) * t323 + pkin(8) * t318;
t347 = -pkin(2) - t366;
t314 = cos(pkin(11));
t461 = pkin(1) * t314;
t243 = t347 - t461;
t207 = t243 * qJD(1);
t125 = -t206 * t317 + t322 * t207;
t112 = -pkin(9) * t254 + t125;
t106 = pkin(4) * t286 + t112;
t126 = t206 * t322 + t207 * t317;
t113 = pkin(9) * t253 + t126;
t109 = t316 * t113;
t56 = t321 * t106 - t109;
t45 = t56 - t543;
t43 = pkin(5) * t282 + t45;
t533 = pkin(10) * t367;
t111 = t321 * t113;
t57 = t106 * t316 + t111;
t46 = t57 + t533;
t437 = t315 * t46;
t14 = t320 * t43 - t437;
t403 = qJD(3) * t318;
t539 = qJD(2) * qJD(3) + t291 * qJDD(1);
t149 = t318 * qJDD(2) - t274 * t403 + t323 * t539;
t141 = qJDD(3) * pkin(8) + t149;
t292 = -pkin(2) - t461;
t272 = t292 * qJDD(1);
t160 = -pkin(3) * t266 - pkin(8) * t267 + t272;
t59 = -qJD(4) * t126 - t141 * t317 + t322 * t160;
t42 = pkin(4) * t248 - pkin(9) * t166 + t59;
t398 = qJD(4) * t322;
t400 = qJD(4) * t317;
t58 = t322 * t141 + t317 * t160 - t206 * t400 + t207 * t398;
t47 = pkin(9) * t167 + t58;
t13 = -qJD(5) * t57 - t316 * t47 + t321 * t42;
t6 = pkin(5) * t241 - pkin(10) * t68 + t13;
t12 = t106 * t396 - t113 * t397 + t316 * t42 + t321 * t47;
t9 = pkin(10) * t69 + t12;
t2 = qJD(6) * t14 + t315 * t6 + t320 * t9;
t432 = t320 * t46;
t15 = t315 * t43 + t432;
t3 = -qJD(6) * t15 - t315 * t9 + t320 * t6;
t538 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t345 = t391 - t538;
t390 = Ifges(6,5) * t68 + Ifges(6,6) * t69 + Ifges(6,3) * t241;
t462 = mrSges(7,3) * t15;
t463 = mrSges(7,3) * t14;
t273 = qJD(6) + t282;
t470 = -t273 / 0.2e1;
t483 = -t105 / 0.2e1;
t485 = -t368 / 0.2e1;
t99 = Ifges(7,4) * t368;
t54 = t105 * Ifges(7,1) + t273 * Ifges(7,5) + t99;
t494 = -t54 / 0.2e1;
t441 = Ifges(7,4) * t105;
t53 = Ifges(7,2) * t368 + t273 * Ifges(7,6) + t441;
t496 = -t53 / 0.2e1;
t536 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t550 = t345 + t390 - t536 + (-mrSges(7,2) * t107 + Ifges(7,1) * t483 + Ifges(7,4) * t485 + Ifges(7,5) * t470 + t463 + t494) * t368 - (mrSges(7,1) * t107 + Ifges(7,4) * t483 + Ifges(7,2) * t485 + Ifges(7,6) * t470 - t462 + t496) * t105;
t312 = qJ(4) + qJ(5);
t307 = cos(t312);
t456 = pkin(4) * t322;
t269 = pkin(5) * t307 + t456;
t261 = pkin(3) + t269;
t308 = qJ(6) + t312;
t297 = sin(t308);
t298 = cos(t308);
t300 = pkin(3) + t456;
t306 = sin(t312);
t362 = -mrSges(5,1) * t322 + mrSges(5,2) * t317;
t549 = m(5) * pkin(3) + m(6) * t300 + m(7) * t261 + mrSges(6,1) * t307 + mrSges(7,1) * t298 - mrSges(6,2) * t306 - mrSges(7,2) * t297 - t362;
t311 = -pkin(10) + t325;
t548 = -m(5) * pkin(8) + m(6) * t325 + m(7) * t311 - t554;
t547 = t266 / 0.2e1;
t546 = t267 / 0.2e1;
t545 = t393 / 0.2e1;
t544 = m(5) * t205;
t455 = pkin(5) * t172;
t542 = m(4) * qJD(3);
t401 = qJD(3) * t323;
t335 = t317 * t401 + t318 * t398;
t519 = -t217 + (-t382 + t400) * pkin(4);
t537 = -t59 * mrSges(5,1) + t58 * mrSges(5,2);
t446 = Ifges(4,4) * t318;
t356 = t323 * Ifges(4,2) + t446;
t535 = t15 * mrSges(7,2) + t57 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t356 / 0.2e1 - t14 * mrSges(7,1) - t56 * mrSges(6,1);
t502 = m(6) * pkin(4);
t500 = t26 / 0.2e1;
t499 = t27 / 0.2e1;
t492 = t68 / 0.2e1;
t491 = t69 / 0.2e1;
t481 = t166 / 0.2e1;
t480 = t167 / 0.2e1;
t475 = t229 / 0.2e1;
t474 = t241 / 0.2e1;
t473 = t248 / 0.2e1;
t532 = -qJD(1) / 0.2e1;
t190 = t321 * t280 + t281 * t316;
t147 = -pkin(10) * t256 + t190;
t148 = -pkin(10) * t348 + t191;
t86 = t147 * t315 + t148 * t320;
t531 = -qJD(6) * t86 + t552 * t315 + t320 * t553;
t85 = t147 * t320 - t148 * t315;
t530 = qJD(6) * t85 + t315 * t553 - t552 * t320;
t457 = pkin(4) * t321;
t299 = pkin(5) + t457;
t394 = qJD(6) * t320;
t395 = qJD(6) * t315;
t419 = t315 * t316;
t60 = -t112 * t316 - t111;
t50 = t60 - t533;
t61 = t321 * t112 - t109;
t51 = t61 - t543;
t527 = -t315 * t50 - t320 * t51 + t299 * t394 + (-t316 * t395 + (t320 * t321 - t419) * qJD(5)) * pkin(4);
t418 = t316 * t320;
t526 = t315 * t51 - t320 * t50 - t299 * t395 + (-t316 * t394 + (-t315 * t321 - t418) * qJD(5)) * pkin(4);
t525 = t502 + mrSges(5,1);
t98 = -mrSges(5,1) * t167 + mrSges(5,2) * t166;
t524 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t267 + t98;
t220 = t348 * t318;
t222 = t322 * t243;
t414 = t318 * t322;
t143 = -pkin(9) * t414 + t222 + (-t291 * t317 - pkin(4)) * t323;
t260 = t291 * t411;
t176 = t317 * t243 + t260;
t416 = t317 * t318;
t157 = -pkin(9) * t416 + t176;
t88 = t316 * t143 + t321 * t157;
t522 = pkin(5) * t551 + t519;
t387 = mrSges(4,3) * t405;
t521 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t253 + mrSges(5,2) * t254 + t387;
t359 = -mrSges(7,1) * t297 - mrSges(7,2) * t298;
t520 = mrSges(6,1) * t306 + mrSges(6,2) * t307 - t359;
t310 = qJ(1) + pkin(11);
t302 = sin(t310);
t303 = cos(t310);
t421 = t306 * t323;
t202 = t302 * t307 - t303 * t421;
t420 = t307 * t323;
t203 = t302 * t306 + t303 * t420;
t422 = t303 * t323;
t188 = -t297 * t422 + t298 * t302;
t189 = t297 * t302 + t298 * t422;
t409 = t188 * mrSges(7,1) - t189 * mrSges(7,2);
t518 = -t202 * mrSges(6,1) + t203 * mrSges(6,2) - t409;
t200 = t302 * t421 + t303 * t307;
t201 = -t302 * t420 + t303 * t306;
t424 = t302 * t323;
t186 = t297 * t424 + t298 * t303;
t187 = t297 * t303 - t298 * t424;
t410 = -t186 * mrSges(7,1) + t187 * mrSges(7,2);
t517 = t200 * mrSges(6,1) - t201 * mrSges(6,2) - t410;
t244 = Ifges(5,4) * t253;
t153 = t254 * Ifges(5,1) + t286 * Ifges(5,5) + t244;
t301 = Ifges(4,4) * t404;
t516 = Ifges(4,1) * t405 + Ifges(4,5) * qJD(3) + t322 * t153 + t301;
t122 = mrSges(5,1) * t248 - mrSges(5,3) * t166;
t123 = -mrSges(5,2) * t248 + mrSges(5,3) * t167;
t515 = -t317 * t122 + t322 * t123;
t150 = qJDD(2) * t323 - t274 * t401 - t318 * t539;
t514 = t149 * t323 - t150 * t318;
t513 = -t317 * t59 + t322 * t58;
t511 = t254 * Ifges(5,5) + t172 * Ifges(6,5) + t105 * Ifges(7,5) + t253 * Ifges(5,6) + Ifges(6,6) * t367 + Ifges(7,6) * t368 + t286 * Ifges(5,3) + t282 * Ifges(6,3) + t273 * Ifges(7,3);
t510 = -m(6) - m(5) - m(7) - m(4);
t454 = pkin(5) * t306;
t458 = pkin(4) * t317;
t268 = t454 + t458;
t508 = -m(7) * t268 + mrSges(3,2) - mrSges(4,3);
t364 = t323 * mrSges(4,1) - mrSges(4,2) * t318;
t507 = t318 * t554 + mrSges(3,1) + t364;
t504 = Ifges(7,4) * t500 + Ifges(7,2) * t499 + Ifges(7,6) * t475;
t503 = Ifges(7,1) * t500 + Ifges(7,4) * t499 + Ifges(7,5) * t475;
t501 = m(7) * pkin(5);
t498 = Ifges(6,4) * t492 + Ifges(6,2) * t491 + Ifges(6,6) * t474;
t497 = Ifges(6,1) * t492 + Ifges(6,4) * t491 + Ifges(6,5) * t474;
t495 = t53 / 0.2e1;
t493 = t54 / 0.2e1;
t490 = Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * t473;
t442 = Ifges(6,4) * t172;
t92 = Ifges(6,2) * t367 + t282 * Ifges(6,6) + t442;
t489 = -t92 / 0.2e1;
t488 = t92 / 0.2e1;
t168 = Ifges(6,4) * t367;
t93 = t172 * Ifges(6,1) + t282 * Ifges(6,5) + t168;
t487 = -t93 / 0.2e1;
t486 = t93 / 0.2e1;
t484 = t368 / 0.2e1;
t482 = t105 / 0.2e1;
t479 = -t367 / 0.2e1;
t478 = t367 / 0.2e1;
t477 = -t172 / 0.2e1;
t476 = t172 / 0.2e1;
t471 = t254 / 0.2e1;
t469 = t273 / 0.2e1;
t468 = -t282 / 0.2e1;
t467 = t282 / 0.2e1;
t465 = mrSges(6,3) * t56;
t464 = mrSges(6,3) * t57;
t319 = sin(qJ(1));
t460 = pkin(1) * t319;
t459 = pkin(4) * t254;
t451 = g(3) * t318;
t324 = cos(qJ(1));
t309 = t324 * pkin(1);
t448 = mrSges(6,3) * t367;
t447 = mrSges(6,3) * t172;
t445 = Ifges(4,4) * t323;
t444 = Ifges(5,4) * t317;
t443 = Ifges(5,4) * t322;
t440 = t125 * mrSges(5,3);
t439 = t126 * mrSges(5,3);
t438 = t254 * Ifges(5,4);
t425 = t302 * t317;
t423 = t303 * t317;
t415 = t317 * t323;
t279 = t318 * t291;
t265 = t365 * qJD(3);
t381 = t291 * t403;
t406 = t322 * t265 + t317 * t381;
t225 = pkin(4) * t416 + t279;
t399 = qJD(4) * t318;
t386 = mrSges(4,3) * t404;
t385 = Ifges(5,5) * t166 + Ifges(5,6) * t167 + Ifges(5,3) * t248;
t270 = t291 * t401;
t180 = pkin(4) * t335 + t270;
t152 = t253 * Ifges(5,2) + t286 * Ifges(5,6) + t438;
t378 = -t317 * t152 / 0.2e1;
t87 = t321 * t143 - t157 * t316;
t363 = mrSges(4,1) * t318 + mrSges(4,2) * t323;
t361 = mrSges(5,1) * t317 + mrSges(5,2) * t322;
t358 = Ifges(5,1) * t322 - t444;
t357 = Ifges(5,1) * t317 + t443;
t355 = -Ifges(5,2) * t317 + t443;
t354 = Ifges(5,2) * t322 + t444;
t353 = Ifges(4,5) * t323 - Ifges(4,6) * t318;
t352 = Ifges(5,5) * t322 - Ifges(5,6) * t317;
t351 = Ifges(5,5) * t317 + Ifges(5,6) * t322;
t70 = -pkin(5) * t323 + pkin(10) * t220 + t87;
t219 = t256 * t318;
t71 = -pkin(10) * t219 + t88;
t37 = -t315 * t71 + t320 * t70;
t38 = t315 * t70 + t320 * t71;
t139 = -t219 * t320 + t220 * t315;
t140 = -t219 * t315 - t220 * t320;
t173 = -t256 * t315 - t320 * t348;
t174 = t256 * t320 - t315 * t348;
t350 = t261 * t323 - t311 * t318;
t349 = t300 * t323 - t318 * t325;
t213 = t302 * t322 - t303 * t415;
t211 = t302 * t415 + t303 * t322;
t344 = t205 * t361;
t343 = t292 * qJD(1) * t363;
t342 = t318 * (Ifges(4,1) * t323 - t446);
t89 = t346 * qJD(3) + (-t260 + (pkin(9) * t318 - t243) * t317) * qJD(4) + t406;
t114 = t243 * t398 + t317 * t265 + (-t318 * t402 - t323 * t400) * t291;
t95 = -pkin(9) * t335 + t114;
t32 = t143 * t396 - t157 * t397 + t316 * t89 + t321 * t95;
t336 = -t317 * t399 + t322 * t401;
t142 = -qJDD(3) * pkin(3) - t150;
t334 = Ifges(5,5) * t318 + t323 * t358;
t333 = Ifges(5,6) * t318 + t323 * t355;
t332 = Ifges(5,3) * t318 + t323 * t352;
t33 = -qJD(5) * t88 - t316 * t95 + t321 * t89;
t94 = -pkin(4) * t167 + t142;
t330 = (-t125 * t322 - t126 * t317) * qJD(4) + t513;
t277 = -qJD(3) * mrSges(4,2) + t386;
t242 = t361 * t318;
t227 = pkin(4) * t418 + t299 * t315;
t226 = -pkin(4) * t419 + t299 * t320;
t223 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t266;
t214 = t303 * t411 + t425;
t212 = -t302 * t411 + t423;
t210 = pkin(5) * t348 - t300;
t199 = mrSges(5,1) * t286 - mrSges(5,3) * t254;
t198 = -mrSges(5,2) * t286 + mrSges(5,3) * t253;
t175 = -t291 * t415 + t222;
t161 = pkin(5) * t219 + t225;
t136 = mrSges(6,1) * t282 - t447;
t135 = -mrSges(6,2) * t282 + t448;
t131 = -t208 * t315 - t209 * t320;
t130 = -t208 * t320 + t209 * t315;
t127 = t459 + t455;
t117 = -qJD(3) * t341 + t220 * t512;
t116 = -qJD(3) * t340 - t178 * t318;
t115 = -qJD(4) * t176 + t406;
t108 = -mrSges(6,1) * t367 + mrSges(6,2) * t172;
t84 = mrSges(7,1) * t273 - mrSges(7,3) * t105;
t83 = -mrSges(7,2) * t273 + mrSges(7,3) * t368;
t80 = -pkin(5) * t117 + t180;
t78 = t166 * Ifges(5,4) + t167 * Ifges(5,2) + t248 * Ifges(5,6);
t75 = -qJD(6) * t174 + t177 * t315 - t178 * t320;
t74 = qJD(6) * t173 - t177 * t320 - t178 * t315;
t64 = -mrSges(6,2) * t241 + mrSges(6,3) * t69;
t63 = mrSges(6,1) * t241 - mrSges(6,3) * t68;
t55 = -mrSges(7,1) * t368 + mrSges(7,2) * t105;
t49 = -qJD(6) * t140 - t116 * t315 + t117 * t320;
t48 = qJD(6) * t139 + t116 * t320 + t117 * t315;
t39 = -pkin(5) * t69 + t94;
t36 = -mrSges(6,1) * t69 + mrSges(6,2) * t68;
t23 = pkin(10) * t117 + t32;
t22 = pkin(5) * t403 - pkin(10) * t116 + t33;
t21 = -mrSges(7,2) * t229 + mrSges(7,3) * t27;
t20 = mrSges(7,1) * t229 - mrSges(7,3) * t26;
t17 = t320 * t45 - t437;
t16 = -t315 * t45 - t432;
t10 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t5 = -qJD(6) * t38 + t22 * t320 - t23 * t315;
t4 = qJD(6) * t37 + t22 * t315 + t23 * t320;
t1 = [(Ifges(4,1) * t267 + Ifges(4,5) * qJDD(3) + t352 * t473 + t355 * t480 + t358 * t481 + (m(5) * t142 - t217 * t542) * t291) * t318 + (t511 / 0.2e1 - t217 * mrSges(4,3) + t125 * mrSges(5,1) - t126 * mrSges(5,2) + Ifges(6,3) * t467 + Ifges(7,3) * t469 + Ifges(6,5) * t476 + Ifges(6,6) * t478 + Ifges(7,5) * t482 + Ifges(7,6) * t484 - t535) * t403 + qJD(3) ^ 2 * t353 / 0.2e1 + (-t216 * mrSges(4,3) + t516 / 0.2e1 + t378 + t291 * t544) * t401 + t94 * (mrSges(6,1) * t219 - mrSges(6,2) * t220) - (t322 * t152 + t317 * t153) * t399 / 0.2e1 + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t469 + (Ifges(7,5) * t140 + Ifges(7,6) * t139) * t475 + (-Ifges(6,1) * t220 - Ifges(6,4) * t219) * t492 + (Ifges(6,1) * t116 + Ifges(6,4) * t117) * t476 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t314 - 0.2e1 * mrSges(3,2) * t313 + m(3) * (t313 ^ 2 + t314 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t205 * (mrSges(5,1) * t335 + mrSges(5,2) * t336) + (-t125 * t336 - t126 * t335 - t414 * t59 - t416 * t58) * mrSges(5,3) + qJD(3) * t343 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t482 + (Ifges(7,1) * t140 + Ifges(7,4) * t139) * t500 + t524 * t279 + t514 * mrSges(4,3) + t521 * t270 + m(7) * (t107 * t80 + t14 * t5 + t15 * t4 + t161 * t39 + t2 * t38 + t3 * t37) + m(6) * (t12 * t88 + t13 * t87 + t165 * t180 + t225 * t94 + t32 * t57 + t33 * t56) - t272 * t364 + (-Ifges(6,5) * t220 - Ifges(6,6) * t219) * t474 + (Ifges(6,5) * t116 + Ifges(6,6) * t117) * t467 + m(5) * (t114 * t126 + t115 * t125 + t175 * t59 + t176 * t58) + t342 * t545 + t445 * t546 + (t446 + t356) * t547 + (Ifges(7,4) * t48 + Ifges(7,2) * t49) * t484 + (Ifges(7,4) * t140 + Ifges(7,2) * t139) * t499 + (-t425 * t502 - m(3) * t309 - mrSges(2,1) * t324 - t214 * mrSges(5,1) - t203 * mrSges(6,1) - t189 * mrSges(7,1) + mrSges(2,2) * t319 - t213 * mrSges(5,2) - t202 * mrSges(6,2) - t188 * mrSges(7,2) + t510 * (t303 * pkin(2) + t302 * pkin(7) + t309) + t508 * t302 + (-m(5) * t366 - m(6) * t349 - m(7) * t350 - t507) * t303) * g(2) + (-t423 * t502 + m(3) * t460 + mrSges(2,1) * t319 - t212 * mrSges(5,1) - t201 * mrSges(6,1) - t187 * mrSges(7,1) + mrSges(2,2) * t324 - t211 * mrSges(5,2) - t200 * mrSges(6,2) - t186 * mrSges(7,2) + t510 * (t303 * pkin(7) - t460) + t508 * t303 + (-m(6) * (-pkin(2) - t349) - m(5) * t347 - m(7) * (-pkin(2) - t350) + m(4) * pkin(2) + t507) * t302) * g(1) + t292 * (-mrSges(4,1) * t266 + mrSges(4,2) * t267) + t142 * t242 + t225 * t36 + t114 * t198 + t115 * t199 + t180 * t108 + t175 * t122 + t176 * t123 + t165 * (-mrSges(6,1) * t117 + mrSges(6,2) * t116) + t161 * t10 + t39 * (-mrSges(7,1) * t139 + mrSges(7,2) * t140) + t32 * t135 + t33 * t136 + t107 * (-mrSges(7,1) * t49 + mrSges(7,2) * t48) + t4 * t83 + t5 * t84 + t87 * t63 + t88 * t64 + t80 * t55 - t78 * t416 / 0.2e1 + t37 * t20 + t38 * t21 + t253 * (qJD(3) * t333 - t354 * t399) / 0.2e1 + t286 * (qJD(3) * t332 - t351 * t399) / 0.2e1 - t277 * t381 + (t139 * t2 - t14 * t48 - t140 * t3 + t15 * t49) * mrSges(7,3) + m(4) * (t272 * t292 + t291 * t514) + (-t116 * t56 + t117 * t57 - t12 * t219 + t13 * t220) * mrSges(6,3) + (qJD(3) * t334 - t357 * t399) * t471 + t116 * t486 + t117 * t488 + t414 * t490 + t48 * t493 + t49 * t495 - t220 * t497 - t219 * t498 + t140 * t503 + t139 * t504 + (-Ifges(6,4) * t220 - Ifges(6,2) * t219) * t491 + (Ifges(6,4) * t116 + Ifges(6,2) * t117) * t478 + (-Ifges(5,3) * t473 - Ifges(7,3) * t475 - Ifges(5,6) * t480 - Ifges(5,5) * t481 - Ifges(7,6) * t499 - Ifges(7,5) * t500 + Ifges(4,4) * t546 + Ifges(4,2) * t547 - t390 / 0.2e1 - t391 / 0.2e1 - Ifges(6,3) * t474 - Ifges(6,6) * t491 - Ifges(6,5) * t492 + (-Ifges(4,2) * t318 + t445) * t545 - t385 / 0.2e1 + Ifges(4,6) * qJDD(3) + (-t216 * t542 + t223) * t291 + t536 + t537 + t538) * t323; m(3) * qJDD(2) + t116 * t135 + t117 * t136 + t139 * t20 + t140 * t21 - t219 * t63 - t220 * t64 + t48 * t83 + t49 * t84 + (-m(3) + t510) * g(3) + (-t10 - t36 + (t198 * t322 - t199 * t317 + t277) * qJD(3) - t524) * t323 + (t223 + (-t317 * t198 - t322 * t199) * qJD(4) + (t108 + t55 + t521) * qJD(3) + t515) * t318 + m(7) * (t107 * t403 + t139 * t3 + t14 * t49 + t140 * t2 + t15 * t48 - t323 * t39) + m(4) * (t149 * t318 + t150 * t323 + (-t216 * t318 + t217 * t323) * qJD(3)) + m(6) * (t116 * t57 + t117 * t56 - t12 * t220 - t13 * t219 + t165 * t403 - t323 * t94) + m(5) * ((-t142 + (-t125 * t317 + t126 * t322) * qJD(3)) * t323 + (qJD(3) * t205 + t330) * t318); (g(1) * t303 + g(2) * t302) * (t318 * t549 + t323 * t548 + t363) + (t318 * t548 - t323 * t549 - t364) * g(3) - t511 * t405 / 0.2e1 + (Ifges(6,5) * t477 + Ifges(7,5) * t483 + Ifges(6,6) * t479 + Ifges(7,6) * t485 + Ifges(6,3) * t468 + Ifges(7,3) * t470 + t535) * t405 + (t387 - t521 - t544) * t217 + (t253 * t355 + t254 * t358 + t286 * t352) * qJD(4) / 0.2e1 + t322 * t78 / 0.2e1 + (-Ifges(6,1) * t209 - Ifges(6,4) * t208) * t477 + (-Ifges(6,5) * t209 - Ifges(6,6) * t208) * t468 + (-t12 * t348 - t13 * t256 + t208 * t57 - t209 * t56) * mrSges(6,3) + (-Ifges(6,4) * t209 - Ifges(6,2) * t208) * t479 + (t344 + t378) * qJD(4) - (Ifges(6,4) * t476 + Ifges(6,2) * t478 + Ifges(6,6) * t467 + t464 + t488) * t178 + (Ifges(7,4) * t131 + Ifges(7,2) * t130) * t485 + (t253 * t333 + t254 * t334 + t286 * t332) * t532 + (Ifges(7,5) * t131 + Ifges(7,6) * t130) * t470 + t513 * mrSges(5,3) + (m(5) * t330 - t198 * t400 - t199 * t398 + t515) * pkin(8) - (-Ifges(4,2) * t405 + t301 + t516) * t404 / 0.2e1 + t519 * t108 + t142 * t362 + (-t440 + t153 / 0.2e1) * t398 - (Ifges(6,1) * t476 + Ifges(6,4) * t478 + Ifges(6,5) * t467 - t465 + t486) * t177 + t528 * t136 + (t12 * t191 + t13 * t190 + t165 * t519 - t300 * t94 + t528 * t56 + t529 * t57) * m(6) + t529 * t135 + t522 * t55 + (Ifges(6,5) * t256 - Ifges(6,6) * t348) * t474 + (Ifges(6,4) * t256 - Ifges(6,2) * t348) * t491 + (Ifges(6,1) * t256 - Ifges(6,4) * t348) * t492 + t94 * (mrSges(6,1) * t348 + mrSges(6,2) * t256) - t348 * t498 + (Ifges(7,1) * t131 + Ifges(7,4) * t130) * t483 + t530 * t83 + (t107 * t522 + t14 * t531 + t15 * t530 + t2 * t86 + t210 * t39 + t3 * t85) * m(7) + t531 * t84 + (-t343 - t126 * (-mrSges(5,2) * t318 - mrSges(5,3) * t415) - t125 * (mrSges(5,1) * t318 - mrSges(5,3) * t411) + t342 * t532) * qJD(1) - t300 * t36 + Ifges(4,6) * t266 + Ifges(4,5) * t267 + t210 * t10 - t156 * t198 - t155 * t199 + t190 * t63 + t191 * t64 + t39 * (-mrSges(7,1) * t173 + mrSges(7,2) * t174) + (Ifges(7,4) * t482 + Ifges(7,2) * t484 + Ifges(7,6) * t469 + t462 + t495) * t75 - t149 * mrSges(4,2) + t150 * mrSges(4,1) - t400 * t439 + (Ifges(7,1) * t482 + Ifges(7,4) * t484 + Ifges(7,5) * t469 - t463 + t493) * t74 + (-t130 * t15 + t131 * t14 + t173 * t2 - t174 * t3) * mrSges(7,3) - pkin(3) * t98 + t85 * t20 + t86 * t21 + Ifges(4,3) * qJDD(3) - t344 * t404 - t353 * t393 / 0.2e1 + t152 * t382 / 0.2e1 + (t386 - t277) * t216 + (t551 * mrSges(6,1) - mrSges(6,2) * t555) * t165 + t351 * t473 + (Ifges(7,5) * t174 + Ifges(7,6) * t173) * t475 + t354 * t480 + t357 * t481 - t209 * t487 - t208 * t489 + t317 * t490 + t131 * t494 + t130 * t496 + t256 * t497 + (Ifges(7,4) * t174 + Ifges(7,2) * t173) * t499 + (Ifges(7,1) * t174 + Ifges(7,4) * t173) * t500 + t174 * t503 + t173 * t504 + (-pkin(3) * t142 - t125 * t155 - t126 * t156) * m(5) + ((-t131 + t74) * mrSges(7,2) + (t130 - t75) * mrSges(7,1)) * t107; -(-Ifges(5,2) * t254 + t153 + t244) * t253 / 0.2e1 + (-mrSges(6,2) * t165 + Ifges(6,1) * t477 + Ifges(6,4) * t479 + Ifges(6,5) * t468 + t465 + t487) * t367 + (m(6) * t458 + t520) * t451 + t254 * t439 + t253 * t440 + (-mrSges(5,2) * t212 - m(7) * (-t268 * t424 - t269 * t303) + t525 * t211 + t517) * g(2) + (mrSges(5,2) * t214 - m(7) * (-t268 * t422 + t269 * t302) - t525 * t213 + t518) * g(1) + t526 * t84 + t527 * t83 + (-t107 * t127 + t14 * t526 + t15 * t527 + t2 * t227 + t226 * t3 + t268 * t451) * m(7) - (mrSges(6,1) * t165 + Ifges(6,4) * t477 + Ifges(6,2) * t479 + Ifges(6,6) * t468 - t464 + t489) * t172 + t385 - t286 * (Ifges(5,5) * t253 - Ifges(5,6) * t254) / 0.2e1 - t205 * (mrSges(5,1) * t254 + mrSges(5,2) * t253) + g(3) * t242 + (t135 * t396 - t136 * t397 + t316 * t64) * pkin(4) + t226 * t20 + t227 * t21 - t108 * t459 - m(6) * (t165 * t459 + t56 * t60 + t57 * t61) - t125 * t198 + t126 * t199 - t61 * t135 - t60 * t136 - t127 * t55 - t254 * (Ifges(5,1) * t253 - t438) / 0.2e1 + t63 * t457 + t152 * t471 + (t12 * t316 + t13 * t321 + (-t316 * t56 + t321 * t57) * qJD(5)) * t502 + t550 - t537; -t55 * t455 - m(7) * (t107 * t455 + t14 * t16 + t15 * t17) - t165 * (mrSges(6,1) * t172 + mrSges(6,2) * t367) - t16 * t84 - t17 * t83 + (Ifges(6,5) * t367 - Ifges(6,6) * t172) * t468 + t92 * t476 + (Ifges(6,1) * t367 - t442) * t477 + (t2 * t315 + t3 * t320 + (-t14 * t315 + t15 * t320) * qJD(6)) * t501 + (t136 + t447) * t57 + (-t135 + t448) * t56 + (-Ifges(6,2) * t172 + t168 + t93) * t479 + (m(7) * t454 + t520) * t451 + (t200 * t501 + t517) * g(2) + (-t202 * t501 + t518) * g(1) + (t20 * t320 + t21 * t315 + t394 * t83 - t395 * t84) * pkin(5) + t550; -t107 * (mrSges(7,1) * t105 + mrSges(7,2) * t368) + (Ifges(7,1) * t368 - t441) * t483 + t53 * t482 + (Ifges(7,5) * t368 - Ifges(7,6) * t105) * t470 - t14 * t83 + t15 * t84 - g(1) * t409 - g(2) * t410 - t359 * t451 + (t105 * t15 + t14 * t368) * mrSges(7,3) + t345 + (-Ifges(7,2) * t105 + t54 + t99) * t485;];
tau  = t1;
