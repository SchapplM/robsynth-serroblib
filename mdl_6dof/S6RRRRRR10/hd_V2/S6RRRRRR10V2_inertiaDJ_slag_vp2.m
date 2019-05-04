% Calculate time derivative of joint inertia matrix for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR10V2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:42
% EndTime: 2019-04-11 14:42:13
% DurationCPUTime: 9.81s
% Computational Cost: add. (9646->759), mult. (23489->1126), div. (0->0), fcn. (22717->10), ass. (0->354)
t521 = 2 * pkin(6);
t299 = sin(qJ(5));
t300 = sin(qJ(4));
t409 = qJD(5) * t300;
t384 = t299 * t409;
t304 = cos(qJ(5));
t305 = cos(qJ(4));
t411 = qJD(4) * t305;
t387 = t304 * t411;
t323 = -t384 + t387;
t413 = qJD(4) * t300;
t200 = mrSges(6,1) * t413 - mrSges(6,3) * t323;
t298 = sin(qJ(6));
t303 = cos(qJ(6));
t412 = qJD(4) * t304;
t360 = -qJD(6) + t412;
t361 = -qJD(6) * t304 + qJD(4);
t410 = qJD(5) * t299;
t419 = t303 * t305;
t142 = t360 * t419 + (t298 * t361 - t303 * t410) * t300;
t421 = t300 * t303;
t143 = t361 * t421 + (-t305 * t360 + t384) * t298;
t85 = -mrSges(7,1) * t143 + mrSges(7,2) * t142;
t512 = t85 - t200;
t520 = 0.2e1 * t512;
t295 = t300 ^ 2;
t297 = t305 ^ 2;
t519 = -mrSges(4,2) + (t295 + t297) * mrSges(5,3);
t408 = qJD(5) * t304;
t371 = -t408 / 0.2e1;
t374 = -t411 / 0.2e1;
t518 = t299 * t374 + t300 * t371;
t372 = -t410 / 0.2e1;
t373 = t411 / 0.2e1;
t517 = t300 * t372 + t304 * t373;
t370 = t408 / 0.2e1;
t516 = t299 * t373 + t300 * t370;
t405 = qJD(6) * t299;
t380 = t298 * t405;
t320 = t303 * t408 - t380;
t321 = t299 * t411 + t300 * t408;
t301 = sin(qJ(3));
t302 = sin(qJ(2));
t306 = cos(qJ(3));
t307 = cos(qJ(2));
t239 = t301 * t307 + t302 * t306;
t238 = t301 * t302 - t306 * t307;
t508 = qJD(2) + qJD(3);
t179 = t508 * t238;
t422 = t300 * t179;
t327 = t239 * t411 - t422;
t180 = t508 * t239;
t100 = pkin(2) * qJD(2) * t302 + pkin(3) * t180 + pkin(5) * t179;
t515 = 0.2e1 * t100;
t286 = -pkin(2) * t307 - pkin(1);
t171 = pkin(3) * t238 - pkin(5) * t239 + t286;
t514 = 0.2e1 * t171;
t463 = mrSges(5,2) * t300;
t513 = -mrSges(5,1) * t305 - mrSges(4,1) + t463;
t285 = pkin(2) * t301 + pkin(5);
t393 = -pkin(2) * t306 - pkin(3);
t418 = t304 * t305;
t206 = t285 * t418 + t299 * t393;
t451 = pkin(2) * qJD(3);
t399 = t306 * t451;
t363 = t305 * t399;
t390 = t299 * t413;
t400 = t301 * t451;
t117 = qJD(5) * t206 - t285 * t390 + t299 * t363 - t304 * t400;
t420 = t300 * t304;
t231 = -t298 * t420 - t419;
t232 = -t298 * t305 + t303 * t420;
t166 = -mrSges(7,1) * t231 + mrSges(7,2) * t232;
t254 = -mrSges(6,1) * t305 - mrSges(6,3) * t420;
t416 = -t254 + t166;
t511 = t416 * t117;
t510 = t295 - t297;
t170 = mrSges(6,1) * t321 + mrSges(6,2) * t323;
t462 = mrSges(6,2) * t304;
t345 = mrSges(6,1) * t299 + t462;
t234 = t345 * t300;
t509 = t285 * t170 + t234 * t399;
t506 = 2 * m(5);
t505 = 2 * m(6);
t504 = 2 * m(7);
t503 = 0.2e1 * pkin(5);
t104 = mrSges(7,1) * t321 - mrSges(7,3) * t142;
t502 = 0.2e1 * t104;
t105 = -mrSges(7,2) * t321 + mrSges(7,3) * t143;
t501 = 0.2e1 * t105;
t425 = t299 * t300;
t195 = -mrSges(7,2) * t425 + mrSges(7,3) * t231;
t500 = 0.2e1 * t195;
t196 = mrSges(7,1) * t425 - mrSges(7,3) * t232;
t499 = 0.2e1 * t196;
t201 = -mrSges(6,2) * t413 - mrSges(6,3) * t321;
t498 = 0.2e1 * t201;
t252 = mrSges(6,2) * t305 - mrSges(6,3) * t425;
t497 = 0.2e1 * t252;
t496 = 0.2e1 * t286;
t494 = 0.2e1 * t304;
t493 = m(6) / 0.2e1;
t492 = m(7) / 0.2e1;
t438 = t179 * t305;
t326 = t239 * t413 + t438;
t407 = qJD(5) * t305;
t60 = (-t239 * t407 + t180) * t299 + (qJD(5) * t238 - t326) * t304;
t491 = t60 / 0.2e1;
t71 = Ifges(7,1) * t142 + Ifges(7,4) * t143 + Ifges(7,5) * t321;
t490 = t71 / 0.2e1;
t149 = Ifges(7,5) * t232 + Ifges(7,6) * t231 + Ifges(7,3) * t425;
t489 = t149 / 0.2e1;
t151 = Ifges(7,1) * t232 + Ifges(7,4) * t231 + Ifges(7,5) * t425;
t488 = t151 / 0.2e1;
t457 = Ifges(7,4) * t298;
t266 = Ifges(7,2) * t303 + t457;
t456 = Ifges(7,4) * t303;
t338 = -Ifges(7,2) * t298 + t456;
t155 = -t266 * t405 + (Ifges(7,6) * t299 + t304 * t338) * qJD(5);
t487 = t155 / 0.2e1;
t269 = Ifges(7,1) * t298 + t456;
t341 = Ifges(7,1) * t303 - t457;
t157 = -t269 * t405 + (Ifges(7,5) * t299 + t304 * t341) * qJD(5);
t486 = t157 / 0.2e1;
t423 = t299 * t305;
t161 = -t304 * t238 + t239 * t423;
t485 = t161 / 0.2e1;
t162 = t238 * t299 + t239 * t418;
t484 = t162 / 0.2e1;
t213 = -Ifges(7,6) * t304 + t299 * t338;
t483 = t213 / 0.2e1;
t215 = -Ifges(7,5) * t304 + t299 * t341;
t482 = t215 / 0.2e1;
t481 = t231 / 0.2e1;
t480 = t232 / 0.2e1;
t245 = t338 * qJD(6);
t479 = t245 / 0.2e1;
t248 = t341 * qJD(6);
t478 = t248 / 0.2e1;
t264 = Ifges(7,5) * t298 + Ifges(7,6) * t303;
t477 = t264 / 0.2e1;
t476 = t266 / 0.2e1;
t475 = t269 / 0.2e1;
t474 = -t298 / 0.2e1;
t473 = t298 / 0.2e1;
t472 = -t299 / 0.2e1;
t471 = t299 / 0.2e1;
t470 = t303 / 0.2e1;
t469 = -t304 / 0.2e1;
t468 = -t305 / 0.2e1;
t467 = pkin(6) * t304;
t466 = pkin(6) * t305;
t317 = -qJD(6) * t162 + t327;
t433 = t239 * t300;
t348 = qJD(6) * t433 + t60;
t25 = t298 * t317 + t303 * t348;
t26 = -t298 * t348 + t303 * t317;
t465 = -mrSges(6,1) * t327 - mrSges(7,1) * t26 + mrSges(7,2) * t25 + mrSges(6,3) * t60;
t61 = qJD(5) * t162 - t179 * t423 - t304 * t180 - t239 * t390;
t464 = mrSges(5,1) * t180 - mrSges(6,1) * t61 - mrSges(6,2) * t60 + mrSges(5,3) * t326;
t461 = Ifges(5,4) * t300;
t460 = Ifges(5,4) * t305;
t459 = Ifges(6,4) * t299;
t458 = Ifges(6,4) * t304;
t455 = Ifges(5,5) * t180;
t454 = Ifges(5,5) * t238;
t453 = Ifges(5,6) * t238;
t452 = Ifges(7,6) * t298;
t450 = t180 * Ifges(5,6);
t165 = t171 ^ 2;
t294 = t299 ^ 2;
t388 = t300 * t411;
t444 = t100 * t171;
t397 = t295 * t444;
t449 = (t165 * t388 + t397) * t294;
t431 = t239 * t305;
t448 = mrSges(5,1) * t238 - mrSges(6,1) * t161 - mrSges(6,2) * t162 - mrSges(5,3) * t431;
t261 = -mrSges(7,1) * t303 + mrSges(7,2) * t298;
t447 = t261 - mrSges(6,1);
t446 = mrSges(6,1) * t304 - mrSges(6,2) * t299 + mrSges(5,1);
t108 = -t162 * t298 + t239 * t421;
t109 = t162 * t303 + t298 * t433;
t445 = -mrSges(6,1) * t433 - mrSges(7,1) * t108 + mrSges(7,2) * t109 + mrSges(6,3) * t162;
t443 = t100 * t300;
t322 = -t299 * t407 - t300 * t412;
t349 = t304 * t393;
t116 = qJD(5) * t349 + t285 * t322 + t299 * t400 + t304 * t363;
t442 = t116 * t304;
t205 = t285 * t423 - t349;
t441 = t117 * t205;
t440 = t117 * t299;
t202 = -pkin(3) * t408 + pkin(5) * t322;
t437 = t202 * t304;
t283 = pkin(5) * t418;
t203 = -pkin(3) * t410 - pkin(5) * t390 + qJD(5) * t283;
t255 = pkin(3) * t304 + pkin(5) * t423;
t436 = t203 * t255;
t435 = t203 * t299;
t434 = t205 * t299;
t432 = t239 * t301;
t430 = t255 * t299;
t428 = t297 * t306;
t427 = t298 * t299;
t426 = t299 * t254;
t424 = t299 * t303;
t417 = -Ifges(5,5) * t438 + Ifges(5,3) * t180;
t268 = Ifges(5,2) * t305 + t461;
t414 = qJD(4) * t268;
t406 = qJD(6) * t298;
t404 = qJD(6) * t303;
t265 = Ifges(6,5) * t299 + Ifges(6,6) * t304;
t403 = t265 / 0.2e1 - Ifges(5,6);
t402 = 0.2e1 * t307;
t4 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t61;
t401 = 0.2e1 * t411;
t398 = Ifges(5,6) * t413;
t396 = t205 * t423;
t395 = t255 * t423;
t394 = t254 * t423;
t392 = pkin(5) - t467;
t391 = t171 * t413;
t339 = -Ifges(6,2) * t299 + t458;
t214 = -Ifges(6,6) * t305 + t300 * t339;
t386 = t214 * t411;
t381 = t214 * t408;
t379 = -t425 / 0.2e1;
t378 = t425 / 0.2e1;
t377 = t420 / 0.2e1;
t376 = -t413 / 0.2e1;
t375 = t413 / 0.2e1;
t369 = -t405 / 0.2e1;
t256 = -pkin(3) * t299 + t283;
t368 = t285 - t467;
t367 = t305 * t448;
t366 = t416 * t203;
t365 = t285 * t401;
t17 = Ifges(6,5) * t60 - Ifges(6,6) * t61 + t327 * Ifges(6,3);
t364 = t295 * t399;
t359 = 0.2e1 * t445;
t69 = Ifges(7,5) * t142 + Ifges(7,6) * t143 + t321 * Ifges(7,3);
t340 = -Ifges(5,2) * t300 + t460;
t131 = t239 * t340 + t453;
t86 = Ifges(6,5) * t162 - Ifges(6,6) * t161 + Ifges(6,3) * t433;
t351 = -t131 + t86 - t453;
t174 = -mrSges(5,2) * t238 - mrSges(5,3) * t433;
t350 = -qJD(4) * t174 - t464;
t27 = -pkin(6) * t60 - t100 * t305 + t391;
t337 = pkin(6) * t239 + t171 * t304;
t28 = t337 * t411 + (-pkin(6) * t179 + t100 * t304 - t171 * t410) * t300;
t110 = t337 * t300;
t98 = -pkin(6) * t162 - t171 * t305;
t51 = t110 * t303 + t298 * t98;
t10 = -qJD(6) * t51 + t27 * t303 - t28 * t298;
t50 = -t110 * t298 + t303 * t98;
t9 = qJD(6) * t50 + t27 * t298 + t28 * t303;
t347 = -t10 * t298 + t9 * t303;
t296 = t304 ^ 2;
t346 = -mrSges(5,2) + (t294 + t296) * mrSges(6,3);
t344 = mrSges(7,1) * t298 + mrSges(7,2) * t303;
t343 = Ifges(5,1) * t305 - t461;
t342 = Ifges(6,1) * t304 - t459;
t270 = Ifges(6,1) * t299 + t458;
t267 = Ifges(6,2) * t304 + t459;
t293 = pkin(6) * t413;
t107 = t293 + t116;
t192 = t206 - t466;
t225 = t368 * t300;
t129 = -t192 * t298 + t225 * t303;
t282 = pkin(6) * t384;
t160 = t300 * t399 + t368 * t411 + t282;
t40 = qJD(6) * t129 + t107 * t303 + t160 * t298;
t130 = t192 * t303 + t225 * t298;
t41 = -qJD(6) * t130 - t107 * t298 + t160 * t303;
t336 = -t41 * t298 + t40 * t303;
t335 = -t298 * t51 - t303 * t50;
t226 = t256 - t466;
t257 = t392 * t300;
t163 = -t226 * t298 + t257 * t303;
t185 = t293 + t202;
t204 = t392 * t411 + t282;
t83 = qJD(6) * t163 + t185 * t303 + t204 * t298;
t164 = t226 * t303 + t257 * t298;
t84 = -qJD(6) * t164 - t185 * t298 + t204 * t303;
t334 = -t84 * t298 + t83 * t303;
t333 = -t100 * t234 - t171 * t170;
t332 = t117 * t255 + t203 * t205;
t331 = -t129 * t303 - t130 * t298;
t330 = -t163 * t303 - t164 * t298;
t329 = -t298 * t195 - t303 * t196;
t43 = Ifges(7,4) * t109 + Ifges(7,2) * t108 + Ifges(7,6) * t161;
t44 = Ifges(7,1) * t109 + Ifges(7,4) * t108 + Ifges(7,5) * t161;
t328 = t43 * t474 + t44 * t470;
t325 = t205 * t408 + t440;
t324 = t255 * t408 + t435;
t319 = t298 * t408 + t299 * t404;
t318 = -t303 * t104 - t298 * t105 - t195 * t404;
t120 = -mrSges(6,2) * t433 - mrSges(6,3) * t161;
t316 = t120 * t494 + t299 * t359 + 0.2e1 * t174;
t315 = -t299 * t200 + (-t252 * t299 - t254 * t304) * qJD(5);
t314 = -t298 * t104 - t195 * t406 - t196 * t404;
t150 = Ifges(7,4) * t232 + Ifges(7,2) * t231 + Ifges(7,6) * t425;
t158 = -t270 * t409 + (Ifges(6,5) * t300 + t305 * t342) * qJD(4);
t212 = -Ifges(6,3) * t305 + (Ifges(6,5) * t304 - Ifges(6,6) * t299) * t300;
t216 = -Ifges(6,5) * t305 + t300 * t342;
t247 = t340 * qJD(4);
t250 = t343 * qJD(4);
t271 = Ifges(5,1) * t300 + t460;
t70 = Ifges(7,4) * t142 + Ifges(7,2) * t143 + Ifges(7,6) * t321;
t313 = t142 * t151 + t143 * t150 + t321 * t149 + t158 * t420 + t212 * t413 + t216 * t387 + t231 * t70 + t232 * t71 + t305 * t247 + t300 * t250 + t271 * t411 + t69 * t425;
t242 = (mrSges(5,1) * t300 + mrSges(5,2) * t305) * qJD(4);
t312 = -0.2e1 * pkin(3) * t242 + t313;
t154 = t323 * Ifges(6,5) - Ifges(6,6) * t321 + Ifges(6,3) * t413;
t153 = t320 * Ifges(7,5) - Ifges(7,6) * t319 + Ifges(7,3) * t410;
t289 = Ifges(7,5) * t404;
t243 = -Ifges(7,6) * t406 + t289;
t311 = t303 * pkin(6) * t105 + t142 * t475 + t143 * t476 + t70 * t470 + t71 * t473 - t150 * t406 / 0.2e1 + t404 * t488 + t154 + t231 * t479 + t232 * t478 + t243 * t378 + t516 * t264;
t156 = -t267 * t409 + (Ifges(6,6) * t300 + t305 * t339) * qJD(4);
t211 = -Ifges(7,3) * t304 + (Ifges(7,5) * t303 - t452) * t299;
t290 = Ifges(6,5) * t408;
t244 = -Ifges(6,6) * t410 + t290;
t246 = t339 * qJD(5);
t249 = t342 * qJD(5);
t291 = Ifges(5,5) * t411;
t310 = -t70 * t427 / 0.2e1 + t249 * t377 + t153 * t378 + t246 * t379 + t214 * t372 + t265 * t375 + t216 * t370 + t304 * t156 / 0.2e1 + t291 + pkin(6) * t196 * t380 + t244 * t468 + t69 * t469 + t158 * t471 + t157 * t480 + t155 * t481 + t142 * t482 + t143 * t483 + t410 * t489 + t424 * t490 + t517 * t270 + t518 * t267 + t516 * t211 + (t298 * t369 + t303 * t370) * t151 + (t298 * t371 + t303 * t369) * t150;
t132 = t239 * t343 + t454;
t18 = Ifges(6,4) * t60 - Ifges(6,2) * t61 + Ifges(6,6) * t327;
t19 = Ifges(6,1) * t60 - Ifges(6,4) * t61 + Ifges(6,5) * t327;
t42 = Ifges(7,5) * t109 + Ifges(7,6) * t108 + Ifges(7,3) * t161;
t5 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t61;
t52 = -Ifges(5,4) * t326 - Ifges(5,2) * t327 + t450;
t53 = -Ifges(5,1) * t326 - Ifges(5,4) * t327 + t455;
t6 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t61;
t75 = mrSges(5,1) * t327 - mrSges(5,2) * t326;
t87 = Ifges(6,4) * t162 - Ifges(6,2) * t161 + Ifges(6,6) * t433;
t88 = Ifges(6,1) * t162 - Ifges(6,4) * t161 + Ifges(6,5) * t433;
t309 = t516 * t42 + t517 * t88 + t518 * t87 + (t212 * t373 + t268 * t374 + t271 * t376) * t239 + (t166 * t425 + t252 * t420) * t100 + (t321 * t166 + t201 * t420 + t252 * t387 + t425 * t85) * t171 - t271 * t438 / 0.2e1 + t250 * t431 / 0.2e1 + t131 * t376 + t19 * t377 + t4 * t378 + t18 * t379 - (-t268 / 0.2e1 + t212 / 0.2e1) * t422 + t132 * t373 + t86 * t375 + t238 * (t291 - t398) / 0.2e1 + t180 * (Ifges(5,5) * t300 + Ifges(5,6) * t305) / 0.2e1 + t305 * t52 / 0.2e1 + t300 * t53 / 0.2e1 + t234 * t391 + t9 * t195 + t10 * t196 - Ifges(4,5) * t179 - Ifges(4,6) * t180 - t161 * t156 / 0.2e1 + t142 * t44 / 0.2e1 + t143 * t43 / 0.2e1 + t26 * t150 / 0.2e1 + t50 * t104 + t51 * t105 + t108 * t70 / 0.2e1 - pkin(3) * t75 + (-t247 / 0.2e1 + t154 / 0.2e1) * t433 + (-t214 / 0.2e1 + t489) * t61 + t17 * t468 + t6 * t480 + t5 * t481 + t158 * t484 + t69 * t485 + t25 * t488 + t109 * t490 + t216 * t491;
t272 = pkin(5) * t364;
t253 = -mrSges(7,1) * t304 - mrSges(7,3) * t424;
t251 = mrSges(7,2) * t304 - mrSges(7,3) * t427;
t241 = t345 * qJD(5);
t240 = t344 * qJD(6);
t236 = t285 * t364;
t233 = t344 * t299;
t199 = -mrSges(7,2) * t410 - mrSges(7,3) * t319;
t198 = mrSges(7,1) * t410 - mrSges(7,3) * t320;
t169 = mrSges(7,1) * t319 + mrSges(7,2) * t320;
t93 = -mrSges(5,2) * t180 - mrSges(5,3) * t327;
t82 = mrSges(7,1) * t161 - mrSges(7,3) * t109;
t81 = -mrSges(7,2) * t161 + mrSges(7,3) * t108;
t74 = t297 * t444;
t34 = -mrSges(6,2) * t327 - mrSges(6,3) * t61;
t13 = -mrSges(7,2) * t61 + mrSges(7,3) * t26;
t12 = mrSges(7,1) * t61 - mrSges(7,3) * t25;
t1 = [(mrSges(4,1) * t180 - mrSges(4,2) * t179) * t496 - 0.2e1 * t179 * t239 * Ifges(4,1) + t162 * t19 + t109 * t6 + t108 * t5 + 0.2e1 * t9 * t81 + 0.2e1 * t10 * t82 + t60 * t88 + t26 * t43 + t25 * t44 + 0.2e1 * t50 * t12 + 0.2e1 * t51 * t13 + (-t87 + t42) * t61 + (t4 - t18) * t161 + (((2 * Ifges(4,2)) + Ifges(5,3)) * t180 + t417) * t238 + (t165 * t295 * t299 * t408 + t10 * t50 + t51 * t9 + t449) * t504 + (t74 + t397) * t506 + (t296 * t397 + t449 + t74) * t505 + 0.2e1 * (t179 * t238 - t180 * t239) * Ifges(4,4) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t307) * t402 + (0.2e1 * pkin(2) * (mrSges(4,1) * t238 + mrSges(4,2) * t239) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t496 - 0.2e1 * Ifges(3,4) * t302 + (Ifges(3,1) - Ifges(3,2)) * t402) * t302) * qJD(2) + (-t179 * t132 + (t53 + t455) * t239 + (t171 * t316 + t239 * t351) * qJD(4) + t448 * t515 + t464 * t514) * t305 + (m(6) * (t296 - 0.1e1) * t165 * t401 - t351 * t179 + (-t450 + t17 - t52 + (-t132 - t454) * qJD(4)) * t239 + t316 * t100 + (t34 * t494 + 0.2e1 * t93 + (-0.2e1 * t120 * t299 + t304 * t359) * qJD(5) - 0.2e1 * qJD(4) * t448 + 0.2e1 * t299 * t465) * t171) * t300; t309 + (-t306 * t75 + (mrSges(5,2) * t432 + t174 * t306) * t305 * qJD(3) + (t179 * t306 - t180 * t301 + (-t238 * t306 + t432) * qJD(3)) * mrSges(4,3)) * pkin(2) + (t285 * t93 + t333) * t305 + t465 * t205 + t445 * t117 + t206 * t34 + t129 * t12 + t130 * t13 + t116 * t120 + t40 * t81 + t41 * t82 + (t350 * t285 + (mrSges(5,1) * t432 - t306 * t448) * t451 + (-t426 + m(7) * t434 + m(6) * (t206 * t304 - t285 * t305 + t434)) * t100 + (m(7) * t325 + m(6) * (-t206 * t410 + t325 - t363 + t442) + t315) * t171) * t300 + (-t285 * t367 + (-t394 + m(7) * t396 + m(6) * (t206 * t418 + t285 * t510 + t396)) * t171) * qJD(4) + m(7) * (t10 * t129 + t130 * t9 + t40 * t51 + t41 * t50) + (Ifges(3,5) * t307 - Ifges(3,6) * t302) * qJD(2); -t156 * t425 - t268 * t413 + t313 + 0.2e1 * t393 * t242 - t299 * t386 - t216 * t384 - t305 * t154 + t234 * t365 + t116 * t497 + t206 * t498 + t41 * t499 + t40 * t500 + t130 * t501 + t129 * t502 + (t129 * t41 + t130 * t40 + t441) * t504 + (t285 ^ 2 * t388 + t116 * t206 + t236 + t441) * t505 + (t236 + (t285 * t428 + t301 * t393) * t451) * t506 + 0.2e1 * t519 * t399 + (-t381 + 0.2e1 * t509) * t300 + t205 * t520 + 0.2e1 * t400 * t513 + 0.2e1 * t511; t309 + t256 * t34 + t202 * t120 + t163 * t12 + t164 * t13 + t83 * t81 + t84 * t82 + (-t100 * t426 + t315 * t171 + t350 * pkin(5) + ((-pkin(5) * t305 + t256 * t304 + t430) * t493 + t430 * t492) * t515 + ((-t256 * t410 + t324 + t437) * t493 + t324 * t492) * t514) * t300 + m(7) * (t10 * t163 + t164 * t9 + t50 * t84 + t51 * t83) + (-pkin(5) * t367 + (-t394 + m(6) * (pkin(5) * t510 + t256 * t418 + t395) + m(7) * t395) * t171) * qJD(4) + (pkin(5) * t93 + t333) * t305 + t465 * t255 + t445 * t203; t312 + (-t154 + (-t214 * t299 + (pkin(5) + t285) * t234) * qJD(4)) * t305 + t366 + t511 + m(5) * t272 + (-t414 - t299 * t156 + (-t214 * t304 - t216 * t299) * qJD(5) + (m(6) * t365 + t170) * pkin(5) + t509) * t300 + ((m(5) * pkin(5) * t428 + (-m(5) * pkin(3) + t513) * t301) * qJD(3) + (t519 * qJD(3) - t242) * t306) * pkin(2) + m(6) * (t116 * t256 + t202 * t206 + t272 + t332) + m(7) * (t129 * t84 + t130 * t83 + t163 * t41 + t164 * t40 + t332) + (t202 + t116) * t252 + (t163 + t129) * t104 + (t164 + t130) * t105 + (t83 + t40) * t195 + (t84 + t41) * t196 + (t256 + t206) * t201 + t512 * (t255 + t205); t312 + t256 * t498 + t202 * t497 + t83 * t500 + t84 * t499 + t163 * t502 + t164 * t501 + (-t386 + (-qJD(5) * t216 - t156) * t300) * t299 + (qJD(4) * t234 * t503 - t154) * t305 + (t170 * t503 - t381 - t414) * t300 + (pkin(5) ^ 2 * t388 + t202 * t256 + t436) * t505 + (t163 * t84 + t164 * t83 + t436) * t504 + 0.2e1 * t366 + t255 * t520; t270 * t491 + t249 * t484 + t9 * t251 + t10 * t253 + t26 * t483 + t25 * t482 + t50 * t198 + t51 * t199 + t108 * t487 + t109 * t486 + (-t267 / 0.2e1 + t211 / 0.2e1) * t61 + (-t246 / 0.2e1 + t153 / 0.2e1) * t161 + (-t171 * t241 + t446 * t100 + (t171 * t346 + t239 * t403) * qJD(4)) * t305 + (t18 / 0.2e1 - t4 / 0.2e1 + (t88 / 0.2e1 + (m(7) * t335 - t298 * t81 - t303 * t82) * pkin(6) + t328) * qJD(5)) * t304 + (t239 * t244 / 0.2e1 + t171 * t233 * t408 - t403 * t179 + t346 * t100 + (-Ifges(5,5) * t239 - t171 * t446) * qJD(4)) * t300 + (t233 * t443 + t19 / 0.2e1 + t6 * t470 + t5 * t474 + (t169 * t300 + t233 * t411) * t171 + (-t303 * t43 / 0.2e1 + t44 * t474) * qJD(6) + (-t87 / 0.2e1 + t42 / 0.2e1) * qJD(5) + (-t303 * t12 + m(7) * (-t10 * t303 - t298 * t9 - t404 * t51 + t406 * t50) - t298 * t13 - t81 * t404 + t82 * t406) * pkin(6)) * t299 + t417; (t285 * t241 + (mrSges(5,2) * t285 - Ifges(5,6)) * qJD(4) - t446 * t399) * t300 + (t442 + t440 + (t205 * t304 - t206 * t299) * qJD(5)) * mrSges(6,3) + (-qJD(4) * t285 * t446 - mrSges(5,2) * t399) * t305 + t40 * t251 + t41 * t253 + t117 * t233 + t205 * t169 + t129 * t198 + t130 * t199 + ((m(7) * t331 + t329) * t408 + (m(7) * (t129 * t406 - t130 * t404 - t298 * t40 - t303 * t41) + t318) * t299) * pkin(6) + t310; t255 * t169 + t83 * t251 + t84 * t253 + t203 * t233 + t163 * t198 + t164 * t199 - t398 + (t437 + t435 + (t255 * t304 - t256 * t299) * qJD(5)) * mrSges(6,3) + (t300 * t241 + (-t305 * t446 + t463) * qJD(4)) * pkin(5) + ((m(7) * t330 + t329) * t408 + (m(7) * (t163 * t406 - t164 * t404 - t298 * t83 - t303 * t84) + t318) * t299) * pkin(6) + t310; (-t153 + t246 + (-t298 * t213 + t303 * t215 + t270 + (-t251 * t298 - t253 * t303) * t521) * qJD(5)) * t304 + (-t298 * t155 + t303 * t157 + t249 + (-t213 * t303 - t215 * t298) * qJD(6) + (t211 - t267 + m(7) * (t298 ^ 2 + t303 ^ 2) * (pkin(6) ^ 2) * t494) * qJD(5) + (-t303 * t198 - t298 * t199 + (-t251 * t303 + t253 * t298) * qJD(6)) * t521) * t299; t5 * t470 + t6 * t473 + t61 * t477 + t26 * t476 + t25 * t475 + t243 * t485 + t108 * t479 + t109 * t478 + t328 * qJD(6) + (t299 * t447 - t462) * t443 + (qJD(6) * t335 + t347) * mrSges(7,3) + ((-mrSges(6,2) * t411 + t409 * t447) * t304 + ((mrSges(6,2) * qJD(5) + t240) * t300 + t447 * t411) * t299) * t171 + (-t82 * t404 + m(7) * (-t404 * t50 - t406 * t51 + t347) + t303 * t13 - t298 * t12 - t81 * t406) * pkin(6) + t17; (qJD(6) * t331 + t336) * mrSges(7,3) + t447 * t117 + t205 * t240 - t116 * mrSges(6,2) + (m(7) * (-t129 * t404 - t130 * t406 + t336) + t314) * pkin(6) + t311; (qJD(6) * t330 + t334) * mrSges(7,3) + t447 * t203 + t255 * t240 - t202 * mrSges(6,2) + t311 + (m(7) * (-t163 * t404 - t164 * t406 + t334) + t314) * pkin(6); t243 * t469 + t290 + (t477 - Ifges(6,6)) * t410 + (pkin(6) * t199 + t487 + t269 * t370 + t248 * t471 + (-pkin(6) * t253 + t266 * t472 + t482) * qJD(6)) * t303 + (-pkin(6) * t198 + t486 + t266 * t371 + t245 * t472 + (-pkin(6) * t251 - t213 / 0.2e1 + t269 * t472) * qJD(6)) * t298; t245 * t303 + t248 * t298 + (-t266 * t298 + t269 * t303) * qJD(6); mrSges(7,1) * t10 - mrSges(7,2) * t9 + t4; mrSges(7,1) * t41 - mrSges(7,2) * t40 + t69; mrSges(7,1) * t84 - mrSges(7,2) * t83 + t69; (-mrSges(7,1) * t320 + mrSges(7,2) * t319) * pkin(6) + t153; t289 + (pkin(6) * t261 - t452) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
