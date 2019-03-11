% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:17
% EndTime: 2019-03-09 14:35:55
% DurationCPUTime: 19.69s
% Computational Cost: add. (18097->788), mult. (45072->1085), div. (0->0), fcn. (33798->10), ass. (0->377)
t300 = sin(qJ(2));
t295 = sin(pkin(6));
t391 = qJD(1) * t295;
t275 = t300 * t391;
t267 = pkin(2) * t275;
t303 = cos(qJ(2));
t335 = pkin(9) * t300 - qJ(3) * t303;
t213 = t335 * t391 + t267;
t364 = t303 * t391;
t296 = cos(pkin(6));
t390 = qJD(1) * t296;
t377 = pkin(1) * t390;
t240 = pkin(8) * t364 + t300 * t377;
t217 = pkin(3) * t364 + t240;
t299 = sin(qJ(4));
t302 = cos(qJ(4));
t139 = -t213 * t299 + t302 * t217;
t385 = qJD(4) * t299;
t304 = -pkin(2) - pkin(9);
t436 = pkin(10) - t304;
t533 = -(-pkin(10) * t299 * t300 + pkin(4) * t303) * t391 - t139 + t436 * t385;
t140 = t302 * t213 + t299 * t217;
t255 = t436 * t302;
t353 = t302 * t275;
t532 = pkin(10) * t353 + qJD(4) * t255 + t140;
t298 = sin(qJ(5));
t382 = qJD(5) * t298;
t450 = cos(qJ(5));
t358 = t450 * qJD(5);
t504 = t450 * qJD(4) + t358;
t193 = -t298 * t385 - t299 * t382 + t302 * t504;
t396 = t298 * t299;
t214 = -t275 * t396 + t353 * t450;
t395 = t193 + t214;
t384 = qJD(4) * t302;
t192 = -t298 * t384 - t299 * t504 - t302 * t382;
t252 = t298 * t302 + t299 * t450;
t215 = t252 * t275;
t394 = t215 - t192;
t272 = t303 * t377;
t531 = qJD(3) - t272;
t254 = t436 * t299;
t319 = t298 * t254 - t255 * t450;
t521 = qJD(5) * t319 + t533 * t298 - t532 * t450;
t476 = pkin(3) + pkin(8);
t520 = pkin(4) * t384 - (-pkin(4) * t302 - t476) * t275 + t531;
t530 = -pkin(11) * t364 + t521;
t529 = t395 * pkin(5) + t394 * pkin(11) + t520;
t204 = -t254 * t450 - t298 * t255;
t522 = -qJD(5) * t204 + t532 * t298 + t533 * t450;
t263 = t275 + qJD(4);
t253 = qJD(5) + t263;
t412 = t253 * Ifges(6,6);
t282 = qJD(2) + t390;
t225 = -t282 * t299 - t302 * t364;
t226 = t282 * t302 - t299 * t364;
t356 = t450 * t225 - t298 * t226;
t421 = t356 * Ifges(6,2);
t321 = t298 * t225 + t226 * t450;
t429 = Ifges(6,4) * t321;
t98 = t412 + t421 + t429;
t479 = -t98 / 0.2e1;
t148 = qJD(6) - t356;
t422 = t148 * Ifges(7,3);
t297 = sin(qJ(6));
t301 = cos(qJ(6));
t130 = t253 * t297 + t301 * t321;
t423 = t130 * Ifges(7,5);
t129 = t253 * t301 - t297 * t321;
t424 = t129 * Ifges(7,6);
t60 = t422 + t423 + t424;
t485 = t60 / 0.2e1;
t158 = t275 * t476 + t282 * t304 + t531;
t403 = qJ(3) * t300;
t207 = (t303 * t304 - pkin(1) - t403) * t295;
t186 = qJD(1) * t207;
t115 = t158 * t299 + t186 * t302;
t102 = pkin(10) * t225 + t115;
t365 = t450 * t102;
t114 = t302 * t158 - t186 * t299;
t101 = -pkin(10) * t226 + t114;
t96 = pkin(4) * t263 + t101;
t43 = t298 * t96 + t365;
t41 = t253 * pkin(11) + t43;
t265 = t282 * qJ(3);
t177 = t265 + t217;
t144 = -pkin(4) * t225 + t177;
t77 = -pkin(5) * t356 - pkin(11) * t321 + t144;
t22 = -t297 * t41 + t301 * t77;
t23 = t297 * t77 + t301 * t41;
t493 = t22 * mrSges(7,1) - t23 * mrSges(7,2);
t528 = t493 + t479 + t485;
t527 = mrSges(4,1) + mrSges(3,3);
t526 = mrSges(4,2) - mrSges(3,1);
t251 = -t302 * t450 + t396;
t287 = t299 * pkin(4) + qJ(3);
t183 = pkin(5) * t252 + pkin(11) * t251 + t287;
t127 = t183 * t301 - t204 * t297;
t525 = qJD(6) * t127 + t529 * t297 + t530 * t301;
t128 = t183 * t297 + t204 * t301;
t524 = -qJD(6) * t128 - t530 * t297 + t529 * t301;
t523 = pkin(5) * t364 - t522;
t135 = mrSges(6,1) * t253 - mrSges(6,3) * t321;
t83 = -mrSges(7,1) * t129 + mrSges(7,2) * t130;
t506 = t83 - t135;
t389 = qJD(2) * t295;
t357 = qJD(1) * t389;
t351 = t303 * t357;
t388 = qJD(2) * t300;
t361 = t299 * t388;
t383 = qJD(4) * t303;
t174 = -t282 * t385 + (-t302 * t383 + t361) * t391;
t360 = t302 * t388;
t175 = -t282 * t384 + (t299 * t383 + t360) * t391;
t89 = qJD(5) * t356 + t450 * t174 + t298 * t175;
t50 = qJD(6) * t129 + t297 * t351 + t301 * t89;
t90 = qJD(5) * t321 + t298 * t174 - t175 * t450;
t26 = mrSges(7,1) * t90 - mrSges(7,3) * t50;
t51 = -qJD(6) * t130 - t297 * t89 + t301 * t351;
t27 = -mrSges(7,2) * t90 + mrSges(7,3) * t51;
t332 = -t297 * t26 + t301 * t27;
t380 = qJD(6) * t301;
t381 = qJD(6) * t297;
t93 = -mrSges(7,2) * t148 + mrSges(7,3) * t129;
t94 = mrSges(7,1) * t148 - mrSges(7,3) * t130;
t519 = -t94 * t380 - t93 * t381 + t332;
t398 = t298 * t102;
t42 = t450 * t96 - t398;
t517 = t394 * t42 - t395 * t43;
t239 = pkin(8) * t275 - t272;
t516 = -qJD(3) - t239;
t346 = mrSges(7,1) * t297 + mrSges(7,2) * t301;
t40 = -t253 * pkin(5) - t42;
t323 = t40 * t346;
t338 = Ifges(7,5) * t301 - Ifges(7,6) * t297;
t426 = Ifges(7,4) * t301;
t341 = -Ifges(7,2) * t297 + t426;
t427 = Ifges(7,4) * t297;
t344 = Ifges(7,1) * t301 - t427;
t452 = t301 / 0.2e1;
t468 = t148 / 0.2e1;
t470 = t130 / 0.2e1;
t472 = t129 / 0.2e1;
t428 = Ifges(7,4) * t130;
t61 = Ifges(7,2) * t129 + Ifges(7,6) * t148 + t428;
t484 = -t61 / 0.2e1;
t126 = Ifges(7,4) * t129;
t62 = Ifges(7,1) * t130 + Ifges(7,5) * t148 + t126;
t515 = t297 * t484 + t338 * t468 + t341 * t472 + t344 * t470 + t62 * t452 + t323;
t456 = -t253 / 0.2e1;
t466 = -t321 / 0.2e1;
t467 = -t356 / 0.2e1;
t147 = Ifges(6,4) * t356;
t413 = t253 * Ifges(6,5);
t419 = t321 * Ifges(6,1);
t99 = t147 + t413 + t419;
t514 = Ifges(6,1) * t466 + Ifges(6,4) * t467 + Ifges(6,5) * t456 - t99 / 0.2e1;
t469 = -t148 / 0.2e1;
t471 = -t130 / 0.2e1;
t473 = -t129 / 0.2e1;
t513 = -Ifges(6,4) * t466 + Ifges(7,5) * t471 - Ifges(6,2) * t467 - Ifges(6,6) * t456 + Ifges(7,6) * t473 + Ifges(7,3) * t469;
t512 = -t282 / 0.2e1;
t511 = t282 / 0.2e1;
t510 = -t391 / 0.2e1;
t400 = t295 * t300;
t283 = pkin(8) * t400;
t449 = pkin(1) * t303;
t366 = -pkin(2) - t449;
t189 = pkin(3) * t400 + t283 + (-pkin(9) + t366) * t296;
t132 = t302 * t189 - t207 * t299;
t399 = t295 * t303;
t324 = -t296 * t302 + t299 * t399;
t110 = pkin(4) * t400 + pkin(10) * t324 + t132;
t133 = t299 * t189 + t302 * t207;
t244 = -t296 * t299 - t302 * t399;
t117 = pkin(10) * t244 + t133;
t505 = t298 * t110 + t450 * t117;
t352 = t300 * t357;
t261 = pkin(2) * t352;
t386 = qJD(3) * t300;
t311 = (qJD(2) * t335 - t386) * t295;
t165 = qJD(1) * t311 + t261;
t286 = t296 * t300 * pkin(1);
t218 = (t399 * t476 + t286) * qJD(2);
t190 = qJD(1) * t218;
t75 = t158 * t384 + t302 * t165 - t186 * t385 + t299 * t190;
t76 = -qJD(4) * t115 - t165 * t299 + t302 * t190;
t502 = t299 * t75 + t302 * t76;
t333 = -t22 * t297 + t23 * t301;
t134 = -mrSges(6,2) * t253 + mrSges(6,3) * t356;
t501 = -t297 * t94 + t301 * t93 + t134;
t262 = qJD(2) * t272;
t264 = t282 * qJD(3);
t354 = t476 * t400;
t330 = qJD(2) * t354;
t166 = -qJD(1) * t330 + t262 + t264;
t121 = -pkin(4) * t175 + t166;
t30 = pkin(5) * t90 - pkin(11) * t89 + t121;
t55 = pkin(4) * t351 - pkin(10) * t174 + t76;
t57 = pkin(10) * t175 + t75;
t10 = -t102 * t382 + t298 * t55 + t96 * t358 + t450 * t57;
t7 = pkin(11) * t351 + t10;
t3 = -qJD(6) * t23 - t297 * t7 + t30 * t301;
t444 = t297 * t3;
t500 = -t22 * t380 - t23 * t381 - t444;
t11 = -qJD(5) * t43 - t298 * t57 + t450 * t55;
t499 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + Ifges(6,5) * t89 - Ifges(6,6) * t90;
t498 = -t144 * mrSges(6,1) + t43 * mrSges(6,3);
t497 = -t144 * mrSges(6,2) + t42 * mrSges(6,3);
t205 = -pkin(2) * t282 - t516;
t336 = -pkin(2) * t303 - t403;
t230 = (-pkin(1) + t336) * t295;
t222 = qJD(1) * t230;
t266 = Ifges(3,4) * t364;
t425 = Ifges(4,6) * t303;
t496 = -(Ifges(4,4) / 0.2e1 - Ifges(3,5) / 0.2e1) * t282 + t114 * mrSges(5,1) + t205 * mrSges(4,1) + t239 * mrSges(3,3) + t42 * mrSges(6,1) + t263 * Ifges(5,3) + t226 * Ifges(5,5) + t225 * Ifges(5,6) + Ifges(3,1) * t275 / 0.2e1 + Ifges(3,5) * t511 + t266 / 0.2e1 + Ifges(4,4) * t512 + (-t300 * Ifges(4,2) - t425) * t510 + t253 * Ifges(6,3) + t321 * Ifges(6,5) + t356 * Ifges(6,6) - t115 * mrSges(5,2) - t222 * mrSges(4,3) - t43 * mrSges(6,2);
t212 = -t265 - t240;
t490 = Ifges(4,5) / 0.2e1;
t495 = (t490 - Ifges(3,6) / 0.2e1) * t282 + t212 * mrSges(4,1) + Ifges(3,6) * t512 + (Ifges(3,4) * t300 + Ifges(3,2) * t303) * t510 + Ifges(4,5) * t511 + (-Ifges(4,6) * t300 - Ifges(4,3) * t303) * t391 / 0.2e1 - t222 * mrSges(4,2) - t240 * mrSges(3,3);
t201 = qJD(4) * t244 + t295 * t361;
t387 = qJD(2) * t303;
t362 = t295 * t387;
t363 = t295 * t388;
t269 = pkin(2) * t363;
t182 = t269 + t311;
t85 = -qJD(4) * t133 - t182 * t299 + t302 * t218;
t68 = pkin(4) * t362 - pkin(10) * t201 + t85;
t202 = qJD(4) * t324 + t295 * t360;
t84 = t302 * t182 + t189 * t384 - t207 * t385 + t299 * t218;
t74 = pkin(10) * t202 + t84;
t20 = -qJD(5) * t505 - t298 * t74 + t450 * t68;
t494 = t76 * mrSges(5,1) - t75 * mrSges(5,2) + Ifges(5,5) * t174 + Ifges(5,6) * t175;
t106 = pkin(5) * t321 - pkin(11) * t356;
t334 = t22 * t301 + t23 * t297;
t477 = t99 / 0.2e1;
t492 = t477 + t413 / 0.2e1 + t147 / 0.2e1 - t334 * mrSges(7,3) + t515;
t432 = Ifges(5,4) * t226;
t142 = t225 * Ifges(5,2) + t263 * Ifges(5,6) + t432;
t223 = Ifges(5,4) * t225;
t143 = t226 * Ifges(5,1) + t263 * Ifges(5,5) + t223;
t329 = t114 * t299 - t115 * t302;
t430 = Ifges(5,4) * t302;
t431 = Ifges(5,4) * t299;
t453 = -t299 / 0.2e1;
t455 = -t263 / 0.2e1;
t460 = -t226 / 0.2e1;
t461 = -t225 / 0.2e1;
t491 = t177 * (mrSges(5,1) * t302 - mrSges(5,2) * t299) + (Ifges(5,5) * t299 + Ifges(5,6) * t302) * t455 + (Ifges(5,2) * t302 + t431) * t461 + (Ifges(5,1) * t299 + t430) * t460 + t329 * mrSges(5,3) + t143 * t453 - t302 * t142 / 0.2e1;
t489 = Ifges(6,2) / 0.2e1;
t48 = Ifges(7,6) * t51;
t49 = Ifges(7,5) * t50;
t14 = Ifges(7,3) * t90 + t48 + t49;
t488 = t14 / 0.2e1;
t487 = t50 / 0.2e1;
t486 = t51 / 0.2e1;
t483 = t61 / 0.2e1;
t482 = -t62 / 0.2e1;
t481 = -t90 / 0.2e1;
t480 = t90 / 0.2e1;
t475 = pkin(1) * mrSges(3,1);
t474 = pkin(1) * mrSges(3,2);
t320 = t244 * t450 + t298 * t324;
t465 = t320 / 0.2e1;
t169 = t298 * t244 - t324 * t450;
t464 = t169 / 0.2e1;
t463 = t174 / 0.2e1;
t462 = t175 / 0.2e1;
t459 = t226 / 0.2e1;
t458 = t244 / 0.2e1;
t457 = -t324 / 0.2e1;
t454 = t297 / 0.2e1;
t448 = pkin(4) * t226;
t447 = pkin(4) * t298;
t2 = qJD(6) * t22 + t297 * t30 + t301 * t7;
t445 = t2 * t301;
t439 = t89 * Ifges(6,1);
t438 = t89 * Ifges(6,4);
t437 = t90 * Ifges(6,4);
t21 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t80 = mrSges(6,1) * t351 - mrSges(6,3) * t89;
t435 = t21 - t80;
t434 = mrSges(5,3) * t225;
t433 = mrSges(5,3) * t226;
t402 = t356 * t297;
t401 = t356 * t301;
t157 = -mrSges(5,1) * t225 + mrSges(5,2) * t226;
t235 = -mrSges(4,1) * t364 - mrSges(4,3) * t282;
t393 = -t235 + t157;
t392 = t275 * t527 + t526 * t282;
t247 = pkin(8) * t399 + t286;
t379 = t296 * t449;
t378 = t450 * pkin(4);
t373 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t229 = -t296 * qJ(3) - t247;
t206 = pkin(3) * t399 - t229;
t350 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t347 = -mrSges(7,1) * t301 + mrSges(7,2) * t297;
t343 = Ifges(7,1) * t297 + t426;
t340 = Ifges(7,2) * t301 + t427;
t337 = Ifges(7,5) * t297 + Ifges(7,6) * t301;
t64 = pkin(11) * t400 + t505;
t159 = -pkin(4) * t244 + t206;
t92 = -pkin(5) * t320 - pkin(11) * t169 + t159;
t32 = t297 * t92 + t301 * t64;
t31 = -t297 * t64 + t301 * t92;
t154 = mrSges(5,1) * t351 - mrSges(5,3) * t174;
t155 = -mrSges(5,2) * t351 + mrSges(5,3) * t175;
t328 = t302 * t154 + t299 * t155;
t178 = -mrSges(5,2) * t263 + t434;
t179 = mrSges(5,1) * t263 - t433;
t327 = t302 * t178 - t299 * t179;
t273 = qJD(2) * t379;
t241 = -pkin(8) * t363 + t273;
t145 = -t169 * t297 + t301 * t400;
t146 = t169 * t301 + t297 * t400;
t65 = t110 * t450 - t298 * t117;
t19 = t110 * t358 - t117 * t382 + t298 * t68 + t450 * t74;
t227 = -pkin(8) * t352 + t262;
t292 = t296 * qJD(3);
t188 = t273 + t292 - t330;
t315 = (-qJ(3) * t387 - t386) * t295;
t242 = t247 * qJD(2);
t196 = -t227 - t264;
t228 = qJD(1) * t242;
t314 = -t227 * mrSges(3,2) - t196 * mrSges(4,3) + t228 * t526;
t313 = -qJD(6) * t334 - t444;
t138 = -pkin(4) * t202 + t188;
t310 = t313 + t445;
t15 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + t90 * Ifges(7,6);
t16 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + t90 * Ifges(7,5);
t257 = Ifges(6,3) * t351;
t8 = -pkin(5) * t351 - t11;
t308 = mrSges(7,3) * t445 + qJD(6) * t515 + t15 * t452 + t16 * t454 + t337 * t480 + t340 * t486 + t343 * t487 + t8 * t347 + t257 + t499;
t307 = -t429 / 0.2e1 + t424 / 0.2e1 + t423 / 0.2e1 + t422 / 0.2e1 - t412 / 0.2e1 + t528;
t290 = -t378 - pkin(5);
t260 = Ifges(3,5) * t351;
t259 = Ifges(4,5) * t352;
t258 = Ifges(5,3) * t351;
t246 = -t283 + t379;
t238 = -qJ(3) * t364 + t267;
t237 = (mrSges(4,2) * t303 - mrSges(4,3) * t300) * t391;
t234 = -mrSges(3,2) * t282 + mrSges(3,3) * t364;
t231 = t296 * t366 + t283;
t224 = -t241 - t292;
t219 = t269 + t315;
t216 = -qJD(1) * t354 + t272;
t194 = qJD(1) * t315 + t261;
t164 = t215 * t301 + t297 * t364;
t163 = -t215 * t297 + t301 * t364;
t162 = -t214 * t301 + t282 * t297;
t161 = t214 * t297 + t282 * t301;
t122 = -mrSges(5,1) * t175 + mrSges(5,2) * t174;
t113 = t174 * Ifges(5,1) + t175 * Ifges(5,4) + Ifges(5,5) * t351;
t112 = t174 * Ifges(5,4) + t175 * Ifges(5,2) + Ifges(5,6) * t351;
t105 = -mrSges(6,1) * t356 + mrSges(6,2) * t321;
t104 = qJD(5) * t169 + t298 * t201 - t202 * t450;
t103 = qJD(5) * t320 + t201 * t450 + t298 * t202;
t91 = t106 + t448;
t81 = -mrSges(6,2) * t351 - mrSges(6,3) * t90;
t70 = -qJD(6) * t146 - t103 * t297 + t301 * t362;
t69 = qJD(6) * t145 + t103 * t301 + t297 * t362;
t63 = -pkin(5) * t400 - t65;
t47 = t101 * t450 - t398;
t46 = t101 * t298 + t365;
t38 = mrSges(6,1) * t90 + mrSges(6,2) * t89;
t35 = Ifges(6,5) * t351 - t437 + t439;
t34 = -t90 * Ifges(6,2) + Ifges(6,6) * t351 + t438;
t33 = pkin(5) * t104 - pkin(11) * t103 + t138;
t29 = t106 * t297 + t301 * t42;
t28 = t106 * t301 - t297 * t42;
t25 = t297 * t91 + t301 * t47;
t24 = -t297 * t47 + t301 * t91;
t18 = -pkin(5) * t362 - t20;
t17 = pkin(11) * t362 + t19;
t5 = -qJD(6) * t32 - t17 * t297 + t301 * t33;
t4 = qJD(6) * t31 + t17 * t301 + t297 * t33;
t1 = [t505 * t81 + m(6) * (t10 * t505 + t11 * t65 + t121 * t159 + t138 * t144 + t19 * t43 + t20 * t42) + t121 * (-mrSges(6,1) * t320 + mrSges(6,2) * t169) + t3 * (-mrSges(7,1) * t320 - mrSges(7,3) * t146) + t2 * (mrSges(7,2) * t320 + mrSges(7,3) * t145) + (t10 * t320 - t103 * t42 - t104 * t43 - t11 * t169) * mrSges(6,3) + t89 * (Ifges(6,1) * t169 + Ifges(6,4) * t320) / 0.2e1 + (Ifges(7,5) * t146 + Ifges(7,6) * t145 - Ifges(7,3) * t320) * t480 + (Ifges(6,4) * t169 + Ifges(6,2) * t320) * t481 + (Ifges(7,4) * t146 + Ifges(7,2) * t145 - Ifges(7,6) * t320) * t486 + (Ifges(7,1) * t146 + Ifges(7,4) * t145 - Ifges(7,5) * t320) * t487 - t320 * t488 + (-t114 * t201 + t115 * t202 + t244 * t75 + t324 * t76) * mrSges(5,3) + t166 * (-mrSges(5,1) * t244 - mrSges(5,2) * t324) + (-Ifges(5,4) * t324 + Ifges(5,2) * t244) * t462 + (-Ifges(5,1) * t324 + Ifges(5,4) * t244) * t463 + (t300 * t495 + t303 * t496) * t389 + t225 * (Ifges(5,4) * t201 + Ifges(5,2) * t202) / 0.2e1 + (t260 / 0.2e1 + t259 / 0.2e1 + t314) * t296 + ((-mrSges(4,1) * t196 + mrSges(4,2) * t194 + mrSges(3,3) * t227) * t303 + (-t194 * mrSges(4,3) + t257 / 0.2e1 + t258 / 0.2e1 + t527 * t228 + t494 + t499) * t300) * t295 + t263 * (Ifges(5,5) * t201 + Ifges(5,6) * t202) / 0.2e1 + t241 * t234 + t219 * t237 + t224 * t235 + t206 * t122 + t201 * t143 / 0.2e1 + t177 * (-mrSges(5,1) * t202 + mrSges(5,2) * t201) + t202 * t142 / 0.2e1 + t188 * t157 + t84 * t178 + t85 * t179 + t132 * t154 + t133 * t155 + t159 * t38 + t8 * (-mrSges(7,1) * t145 + mrSges(7,2) * t146) + t146 * t16 / 0.2e1 + t144 * (mrSges(6,1) * t104 + mrSges(6,2) * t103) + t145 * t15 / 0.2e1 + t19 * t134 + t20 * t135 + t138 * t105 + t22 * (mrSges(7,1) * t104 - mrSges(7,3) * t69) + t23 * (-mrSges(7,2) * t104 + mrSges(7,3) * t70) + t4 * t93 + t5 * t94 + t18 * t83 + t65 * t80 + t40 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t69 * t62 / 0.2e1 + t63 * t21 + t31 * t26 + t32 * t27 + t392 * t242 + m(4) * (t194 * t230 + t196 * t229 + t205 * t242 + t212 * t224 + t219 * t222 + t228 * t231) + m(3) * (t227 * t247 - t228 * t246 + t239 * t242 + t240 * t241) + m(7) * (t18 * t40 + t2 * t32 + t22 * t5 + t23 * t4 + t3 * t31 + t63 * t8) + m(5) * (t114 * t85 + t115 * t84 + t132 * t76 + t133 * t75 + t166 * t206 + t177 * t188) + t253 * (Ifges(6,5) * t103 - Ifges(6,6) * t104) / 0.2e1 + t356 * (Ifges(6,4) * t103 - Ifges(6,2) * t104) / 0.2e1 + t321 * (Ifges(6,1) * t103 - Ifges(6,4) * t104) / 0.2e1 + ((-t230 * mrSges(4,3) + Ifges(5,5) * t457 + Ifges(5,6) * t458 + Ifges(6,5) * t464 + Ifges(6,6) * t465 + t231 * mrSges(4,1) - t246 * mrSges(3,3) + (-Ifges(4,4) + Ifges(3,5) / 0.2e1) * t296 + (t303 * t373 - 0.2e1 * t474) * t295) * t303 + (t229 * mrSges(4,1) - t230 * mrSges(4,2) - t247 * mrSges(3,3) + (-Ifges(3,6) + t490) * t296 + (-t300 * t373 - 0.2e1 * t475) * t295 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t399) * t300) * t357 + t113 * t457 + t112 * t458 + (Ifges(5,1) * t201 + Ifges(5,4) * t202) * t459 + t35 * t464 + t34 * t465 + (Ifges(7,5) * t69 + Ifges(7,6) * t70 + Ifges(7,3) * t104) * t468 + (Ifges(7,1) * t69 + Ifges(7,4) * t70 + Ifges(7,5) * t104) * t470 + (Ifges(7,4) * t69 + Ifges(7,2) * t70 + Ifges(7,6) * t104) * t472 + t103 * t477 + t104 * t479 + t70 * t483 + t104 * t485; (-t513 + t528) * t214 - t502 * mrSges(5,3) + t514 * t215 - t435 * t319 + (-t23 * t163 + t22 * t164 + (qJD(6) * t333 + t2 * t297 + t3 * t301) * t251) * mrSges(7,3) + ((-m(5) * t329 + t327) * t304 + t491) * qJD(4) + (t419 / 0.2e1 + t492) * t192 + ((-t266 / 0.2e1 + (t474 - t425 / 0.2e1) * t391 - t496) * t303 + (((t475 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t300) * t295 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t399) * qJD(1) + t491 - t495) * t300 + (t336 * mrSges(4,1) - Ifges(4,4) * t303 - Ifges(3,6) * t300 + (Ifges(5,5) * t302 - Ifges(6,5) * t251 - Ifges(5,6) * t299 - Ifges(6,6) * t252) * t303 / 0.2e1) * qJD(2)) * t391 + (Ifges(7,5) * t164 + Ifges(7,6) * t163) * t469 + (Ifges(7,4) * t164 + Ifges(7,2) * t163) * t473 + (Ifges(7,1) * t164 + Ifges(7,4) * t163) * t471 + (-t8 * t346 - t50 * t344 / 0.2e1 - t51 * t341 / 0.2e1 + t338 * t481 - t35 / 0.2e1 - t121 * mrSges(6,2) + t437 / 0.2e1 - t439 / 0.2e1 - t301 * t16 / 0.2e1 + t15 * t454 + t11 * mrSges(6,3) + (t337 * t468 + t340 * t472 + t343 * t470 + t347 * t40 + t452 * t61 + t454 * t62) * qJD(6)) * t251 + t521 * t134 + (t10 * t204 + t11 * t319 + t121 * t287 + t144 * t520 + t42 * t522 + t43 * t521) * m(6) + t259 + t260 - m(5) * (t114 * t139 + t115 * t140 + t177 * t216) + t328 * t304 + (-pkin(2) * t228 - qJ(3) * t196 - t205 * t240 + t212 * t516 - t222 * t238) * m(4) + t314 + (-t235 + t234) * t239 + t522 * t135 + t523 * t83 + t524 * t94 + t525 * t93 + (t127 * t3 + t128 * t2 + t22 * t524 + t23 * t525 - t319 * t8 + t40 * t523) * m(7) + t517 * mrSges(6,3) + t520 * t105 + t166 * (mrSges(5,1) * t299 + mrSges(5,2) * t302) + t302 * t113 / 0.2e1 + t287 * t38 - t238 * t237 - t216 * t157 + t204 * t81 - t140 * t178 - t139 * t179 - t40 * (-mrSges(7,1) * t163 + mrSges(7,2) * t164) + t127 * t26 + t128 * t27 + qJ(3) * t122 + m(5) * (t166 * qJ(3) + t177 * qJD(3) + t304 * t502) - t392 * t240 + t393 * qJD(3) + (mrSges(6,1) * t395 - mrSges(6,2) * t394) * t144 + (-t421 / 0.2e1 + t307) * t193 + (-t438 / 0.2e1 - t10 * mrSges(6,3) + t121 * mrSges(6,1) + t49 / 0.2e1 + t48 / 0.2e1 + t488 - t34 / 0.2e1 + (t489 + Ifges(7,3) / 0.2e1) * t90 + t350) * t252 + t112 * t453 + (-Ifges(5,2) * t299 + t430) * t462 + (Ifges(5,1) * t302 - t431) * t463 + t164 * t482 + t163 * t484; t214 * t134 - t161 * t94 - t162 * t93 + t435 * t251 + t327 * qJD(4) + (-t105 - t393) * t282 + t501 * t193 + (mrSges(4,1) * t387 + (t237 + t327) * t300) * t391 + (t81 + (-t297 * t93 - t301 * t94) * qJD(6) + t332) * t252 + t328 + t506 * t394 + (-t161 * t22 - t162 * t23 + t193 * t333 + t251 * t8 + t252 * t310 + t394 * t40) * m(7) + (t10 * t252 - t11 * t251 - t144 * t282 - t517) * m(6) + (-t177 * t282 - t263 * t329 + t502) * m(5) + (t212 * t282 + t222 * t275 + t228) * m(4); (t338 * t469 + t341 * t473 + t344 * t471 - t323 + t497 + t514) * t356 + (-t60 / 0.2e1 + t98 / 0.2e1 - t493 + t498 + t513) * t321 - t506 * (-pkin(4) * t382 + t46) + t258 + (-t178 + t434) * t114 + (t179 + t433) * t115 + t494 + (m(7) * t310 + t519) * (pkin(11) + t447) + (-Ifges(5,2) * t226 + t143 + t223) * t461 + t290 * t21 - t177 * (mrSges(5,1) * t226 + mrSges(5,2) * t225) - t47 * t134 - t25 * t93 - t24 * t94 + ((t450 * t11 + t10 * t298 + (-t298 * t42 + t43 * t450) * qJD(5)) * pkin(4) - t144 * t448 + t42 * t46 - t43 * t47) * m(6) + t501 * pkin(4) * t358 + (t22 * t401 + t23 * t402 + t500) * mrSges(7,3) + (-t22 * t24 - t23 * t25 - t40 * t46 + t8 * t290 + (t298 * t40 + t333 * t450) * qJD(5) * pkin(4)) * m(7) + t80 * t378 - t105 * t448 + t81 * t447 + (Ifges(5,5) * t225 - Ifges(5,6) * t226) * t455 + t142 * t459 + (Ifges(5,1) * t225 - t432) * t460 + t401 * t482 + t402 * t483 + t308; -t506 * t43 + ((-Ifges(6,1) / 0.2e1 + t489) * t321 - t492 + t497) * t356 + t313 * mrSges(7,3) + t519 * pkin(11) - t42 * t134 - t29 * t93 - t28 * t94 - pkin(5) * t21 + (-t307 + t498) * t321 + t308 + ((t445 + t500) * pkin(11) - t22 * t28 - t23 * t29 - t40 * t43 - pkin(5) * t8) * m(7); -t40 * (mrSges(7,1) * t130 + mrSges(7,2) * t129) + (Ifges(7,1) * t129 - t428) * t471 + t61 * t470 + (Ifges(7,5) * t129 - Ifges(7,6) * t130) * t469 - t22 * t93 + t23 * t94 + (t129 * t22 + t130 * t23) * mrSges(7,3) + t350 + t14 + (-Ifges(7,2) * t130 + t126 + t62) * t473;];
tauc  = t1(:);
