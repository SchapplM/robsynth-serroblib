% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:29:10
% EndTime: 2018-11-23 17:29:27
% DurationCPUTime: 17.81s
% Computational Cost: add. (18097->788), mult. (45072->1085), div. (0->0), fcn. (33798->10), ass. (0->375)
t300 = sin(qJ(2));
t295 = sin(pkin(6));
t389 = qJD(1) * t295;
t275 = t300 * t389;
t267 = pkin(2) * t275;
t303 = cos(qJ(2));
t335 = pkin(9) * t300 - qJ(3) * t303;
t213 = t335 * t389 + t267;
t362 = t303 * t389;
t296 = cos(pkin(6));
t388 = qJD(1) * t296;
t375 = pkin(1) * t388;
t240 = pkin(8) * t362 + t300 * t375;
t217 = pkin(3) * t362 + t240;
t299 = sin(qJ(4));
t302 = cos(qJ(4));
t139 = -t213 * t299 + t302 * t217;
t383 = qJD(4) * t299;
t304 = -pkin(2) - pkin(9);
t434 = pkin(10) - t304;
t531 = -(-pkin(10) * t299 * t300 + pkin(4) * t303) * t389 - t139 + t434 * t383;
t140 = t302 * t213 + t299 * t217;
t255 = t434 * t302;
t353 = t302 * t275;
t530 = pkin(10) * t353 + qJD(4) * t255 + t140;
t298 = sin(qJ(5));
t380 = qJD(5) * t298;
t448 = cos(qJ(5));
t358 = t448 * qJD(5);
t502 = t448 * qJD(4) + t358;
t193 = -t298 * t383 - t299 * t380 + t302 * t502;
t394 = t298 * t299;
t214 = -t275 * t394 + t353 * t448;
t393 = t193 + t214;
t382 = qJD(4) * t302;
t192 = -t298 * t382 - t299 * t502 - t302 * t380;
t252 = t298 * t302 + t299 * t448;
t215 = t252 * t275;
t392 = t215 - t192;
t272 = t303 * t375;
t529 = qJD(3) - t272;
t254 = t434 * t299;
t319 = t298 * t254 - t255 * t448;
t519 = qJD(5) * t319 + t531 * t298 - t530 * t448;
t474 = pkin(3) + pkin(8);
t518 = pkin(4) * t382 - (-pkin(4) * t302 - t474) * t275 + t529;
t528 = -pkin(11) * t362 + t519;
t527 = t393 * pkin(5) + t392 * pkin(11) + t518;
t204 = -t254 * t448 - t298 * t255;
t520 = -qJD(5) * t204 + t530 * t298 + t531 * t448;
t263 = t275 + qJD(4);
t253 = qJD(5) + t263;
t410 = t253 * Ifges(6,6);
t282 = qJD(2) + t388;
t225 = -t282 * t299 - t302 * t362;
t226 = t282 * t302 - t299 * t362;
t356 = t448 * t225 - t298 * t226;
t420 = t356 * Ifges(6,2);
t321 = t298 * t225 + t226 * t448;
t428 = Ifges(6,4) * t321;
t98 = t410 + t420 + t428;
t477 = -t98 / 0.2e1;
t148 = qJD(6) - t356;
t421 = t148 * Ifges(7,3);
t297 = sin(qJ(6));
t301 = cos(qJ(6));
t130 = t253 * t297 + t301 * t321;
t422 = t130 * Ifges(7,5);
t129 = t253 * t301 - t297 * t321;
t423 = t129 * Ifges(7,6);
t60 = t421 + t422 + t423;
t483 = t60 / 0.2e1;
t158 = t275 * t474 + t282 * t304 + t529;
t401 = qJ(3) * t300;
t207 = (t303 * t304 - pkin(1) - t401) * t295;
t186 = qJD(1) * t207;
t115 = t158 * t299 + t186 * t302;
t102 = pkin(10) * t225 + t115;
t363 = t448 * t102;
t114 = t302 * t158 - t186 * t299;
t101 = -pkin(10) * t226 + t114;
t96 = pkin(4) * t263 + t101;
t43 = t298 * t96 + t363;
t41 = t253 * pkin(11) + t43;
t265 = t282 * qJ(3);
t177 = t265 + t217;
t144 = -pkin(4) * t225 + t177;
t77 = -pkin(5) * t356 - pkin(11) * t321 + t144;
t22 = -t297 * t41 + t301 * t77;
t23 = t297 * t77 + t301 * t41;
t491 = t22 * mrSges(7,1) - t23 * mrSges(7,2);
t526 = t491 + t477 + t483;
t525 = mrSges(4,1) + mrSges(3,3);
t524 = mrSges(4,2) - mrSges(3,1);
t251 = -t302 * t448 + t394;
t287 = t299 * pkin(4) + qJ(3);
t183 = pkin(5) * t252 + pkin(11) * t251 + t287;
t127 = t183 * t301 - t204 * t297;
t523 = qJD(6) * t127 + t527 * t297 + t301 * t528;
t128 = t183 * t297 + t204 * t301;
t522 = -qJD(6) * t128 - t297 * t528 + t527 * t301;
t521 = pkin(5) * t362 - t520;
t135 = mrSges(6,1) * t253 - mrSges(6,3) * t321;
t83 = -mrSges(7,1) * t129 + mrSges(7,2) * t130;
t504 = t83 - t135;
t387 = qJD(2) * t295;
t357 = qJD(1) * t387;
t351 = t303 * t357;
t381 = qJD(4) * t303;
t386 = qJD(2) * t300;
t174 = -t282 * t383 + (t299 * t386 - t302 * t381) * t389;
t175 = -t282 * t382 + (t299 * t381 + t302 * t386) * t389;
t89 = qJD(5) * t356 + t448 * t174 + t298 * t175;
t50 = qJD(6) * t129 + t297 * t351 + t301 * t89;
t90 = qJD(5) * t321 + t298 * t174 - t175 * t448;
t26 = mrSges(7,1) * t90 - mrSges(7,3) * t50;
t51 = -qJD(6) * t130 - t297 * t89 + t301 * t351;
t27 = -mrSges(7,2) * t90 + mrSges(7,3) * t51;
t332 = -t297 * t26 + t301 * t27;
t378 = qJD(6) * t301;
t379 = qJD(6) * t297;
t93 = -mrSges(7,2) * t148 + mrSges(7,3) * t129;
t94 = mrSges(7,1) * t148 - mrSges(7,3) * t130;
t517 = -t94 * t378 - t93 * t379 + t332;
t396 = t298 * t102;
t42 = t448 * t96 - t396;
t515 = t392 * t42 - t393 * t43;
t239 = pkin(8) * t275 - t272;
t514 = -qJD(3) - t239;
t346 = mrSges(7,1) * t297 + mrSges(7,2) * t301;
t40 = -t253 * pkin(5) - t42;
t323 = t40 * t346;
t338 = Ifges(7,5) * t301 - Ifges(7,6) * t297;
t425 = Ifges(7,4) * t301;
t341 = -Ifges(7,2) * t297 + t425;
t426 = Ifges(7,4) * t297;
t344 = Ifges(7,1) * t301 - t426;
t450 = t301 / 0.2e1;
t466 = t148 / 0.2e1;
t468 = t130 / 0.2e1;
t470 = t129 / 0.2e1;
t427 = Ifges(7,4) * t130;
t61 = Ifges(7,2) * t129 + Ifges(7,6) * t148 + t427;
t482 = -t61 / 0.2e1;
t126 = Ifges(7,4) * t129;
t62 = Ifges(7,1) * t130 + Ifges(7,5) * t148 + t126;
t513 = t297 * t482 + t338 * t466 + t341 * t470 + t344 * t468 + t62 * t450 + t323;
t454 = -t253 / 0.2e1;
t464 = -t321 / 0.2e1;
t465 = -t356 / 0.2e1;
t147 = Ifges(6,4) * t356;
t411 = t253 * Ifges(6,5);
t418 = t321 * Ifges(6,1);
t99 = t147 + t411 + t418;
t512 = Ifges(6,1) * t464 + Ifges(6,4) * t465 + Ifges(6,5) * t454 - t99 / 0.2e1;
t467 = -t148 / 0.2e1;
t469 = -t130 / 0.2e1;
t471 = -t129 / 0.2e1;
t511 = -Ifges(6,4) * t464 + Ifges(7,5) * t469 - Ifges(6,2) * t465 - Ifges(6,6) * t454 + Ifges(7,6) * t471 + Ifges(7,3) * t467;
t510 = -t282 / 0.2e1;
t509 = t282 / 0.2e1;
t508 = -t389 / 0.2e1;
t398 = t295 * t300;
t283 = pkin(8) * t398;
t447 = pkin(1) * t303;
t364 = -pkin(2) - t447;
t189 = pkin(3) * t398 + t283 + (-pkin(9) + t364) * t296;
t132 = t302 * t189 - t207 * t299;
t397 = t295 * t303;
t324 = -t296 * t302 + t299 * t397;
t110 = pkin(4) * t398 + pkin(10) * t324 + t132;
t133 = t299 * t189 + t302 * t207;
t244 = -t296 * t299 - t302 * t397;
t117 = pkin(10) * t244 + t133;
t503 = t298 * t110 + t448 * t117;
t352 = t300 * t357;
t261 = pkin(2) * t352;
t384 = qJD(3) * t300;
t311 = (qJD(2) * t335 - t384) * t295;
t165 = qJD(1) * t311 + t261;
t286 = t296 * t300 * pkin(1);
t218 = (t397 * t474 + t286) * qJD(2);
t190 = qJD(1) * t218;
t75 = t158 * t382 + t302 * t165 - t186 * t383 + t299 * t190;
t76 = -qJD(4) * t115 - t165 * t299 + t302 * t190;
t500 = t299 * t75 + t302 * t76;
t333 = -t22 * t297 + t23 * t301;
t134 = -mrSges(6,2) * t253 + mrSges(6,3) * t356;
t499 = -t297 * t94 + t301 * t93 + t134;
t262 = qJD(2) * t272;
t264 = t282 * qJD(3);
t354 = t474 * t398;
t330 = qJD(2) * t354;
t166 = -qJD(1) * t330 + t262 + t264;
t121 = -pkin(4) * t175 + t166;
t30 = pkin(5) * t90 - pkin(11) * t89 + t121;
t55 = pkin(4) * t351 - pkin(10) * t174 + t76;
t57 = pkin(10) * t175 + t75;
t10 = -t102 * t380 + t298 * t55 + t96 * t358 + t448 * t57;
t7 = pkin(11) * t351 + t10;
t3 = -qJD(6) * t23 - t297 * t7 + t30 * t301;
t442 = t297 * t3;
t498 = -t22 * t378 - t23 * t379 - t442;
t11 = -qJD(5) * t43 - t298 * t57 + t448 * t55;
t497 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + Ifges(6,5) * t89 - Ifges(6,6) * t90;
t496 = -t144 * mrSges(6,1) + t43 * mrSges(6,3);
t495 = -t144 * mrSges(6,2) + t42 * mrSges(6,3);
t205 = -pkin(2) * t282 - t514;
t336 = -pkin(2) * t303 - t401;
t230 = (-pkin(1) + t336) * t295;
t222 = qJD(1) * t230;
t266 = Ifges(3,4) * t362;
t424 = Ifges(4,6) * t303;
t494 = -(Ifges(4,4) / 0.2e1 - Ifges(3,5) / 0.2e1) * t282 + t114 * mrSges(5,1) + t205 * mrSges(4,1) + t239 * mrSges(3,3) + t42 * mrSges(6,1) + t263 * Ifges(5,3) + t226 * Ifges(5,5) + t225 * Ifges(5,6) + Ifges(3,1) * t275 / 0.2e1 + Ifges(3,5) * t509 + t266 / 0.2e1 + Ifges(4,4) * t510 + (-Ifges(4,2) * t300 - t424) * t508 + t253 * Ifges(6,3) + t321 * Ifges(6,5) + t356 * Ifges(6,6) - t115 * mrSges(5,2) - t222 * mrSges(4,3) - t43 * mrSges(6,2);
t212 = -t265 - t240;
t488 = Ifges(4,5) / 0.2e1;
t493 = (t488 - Ifges(3,6) / 0.2e1) * t282 + t212 * mrSges(4,1) + Ifges(3,6) * t510 + (Ifges(3,4) * t300 + Ifges(3,2) * t303) * t508 + Ifges(4,5) * t509 + (-Ifges(4,6) * t300 - Ifges(4,3) * t303) * t389 / 0.2e1 - t222 * mrSges(4,2) - t240 * mrSges(3,3);
t361 = t295 * t386;
t201 = qJD(4) * t244 + t299 * t361;
t385 = qJD(2) * t303;
t360 = t295 * t385;
t269 = pkin(2) * t361;
t182 = t269 + t311;
t85 = -qJD(4) * t133 - t182 * t299 + t302 * t218;
t68 = pkin(4) * t360 - pkin(10) * t201 + t85;
t202 = qJD(4) * t324 + t302 * t361;
t84 = t302 * t182 + t189 * t382 - t207 * t383 + t299 * t218;
t74 = pkin(10) * t202 + t84;
t20 = -qJD(5) * t503 - t298 * t74 + t448 * t68;
t492 = t76 * mrSges(5,1) - t75 * mrSges(5,2) + Ifges(5,5) * t174 + Ifges(5,6) * t175;
t106 = pkin(5) * t321 - pkin(11) * t356;
t334 = t22 * t301 + t23 * t297;
t475 = t99 / 0.2e1;
t490 = t475 + t411 / 0.2e1 + t147 / 0.2e1 - t334 * mrSges(7,3) + t513;
t414 = t226 * Ifges(5,4);
t142 = t225 * Ifges(5,2) + t263 * Ifges(5,6) + t414;
t223 = Ifges(5,4) * t225;
t143 = t226 * Ifges(5,1) + t263 * Ifges(5,5) + t223;
t329 = t114 * t299 - t115 * t302;
t429 = Ifges(5,4) * t302;
t430 = Ifges(5,4) * t299;
t451 = -t299 / 0.2e1;
t453 = -t263 / 0.2e1;
t458 = -t226 / 0.2e1;
t459 = -t225 / 0.2e1;
t489 = t177 * (mrSges(5,1) * t302 - mrSges(5,2) * t299) + (Ifges(5,5) * t299 + Ifges(5,6) * t302) * t453 + (Ifges(5,2) * t302 + t430) * t459 + (Ifges(5,1) * t299 + t429) * t458 + t329 * mrSges(5,3) + t143 * t451 - t302 * t142 / 0.2e1;
t487 = Ifges(6,2) / 0.2e1;
t48 = Ifges(7,6) * t51;
t49 = Ifges(7,5) * t50;
t14 = Ifges(7,3) * t90 + t48 + t49;
t486 = t14 / 0.2e1;
t485 = t50 / 0.2e1;
t484 = t51 / 0.2e1;
t481 = t61 / 0.2e1;
t480 = -t62 / 0.2e1;
t479 = -t90 / 0.2e1;
t478 = t90 / 0.2e1;
t473 = pkin(1) * mrSges(3,1);
t472 = pkin(1) * mrSges(3,2);
t320 = t244 * t448 + t298 * t324;
t463 = t320 / 0.2e1;
t169 = t298 * t244 - t324 * t448;
t462 = t169 / 0.2e1;
t461 = t174 / 0.2e1;
t460 = t175 / 0.2e1;
t457 = t226 / 0.2e1;
t456 = t244 / 0.2e1;
t455 = -t324 / 0.2e1;
t452 = t297 / 0.2e1;
t446 = pkin(4) * t226;
t445 = pkin(4) * t298;
t2 = qJD(6) * t22 + t297 * t30 + t301 * t7;
t443 = t2 * t301;
t437 = t89 * Ifges(6,1);
t436 = t89 * Ifges(6,4);
t435 = t90 * Ifges(6,4);
t21 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t80 = mrSges(6,1) * t351 - mrSges(6,3) * t89;
t433 = t21 - t80;
t432 = mrSges(5,3) * t225;
t431 = mrSges(5,3) * t226;
t400 = t356 * t297;
t399 = t356 * t301;
t157 = -mrSges(5,1) * t225 + mrSges(5,2) * t226;
t235 = -mrSges(4,1) * t362 - mrSges(4,3) * t282;
t391 = -t235 + t157;
t390 = t275 * t525 + t524 * t282;
t247 = pkin(8) * t397 + t286;
t377 = t296 * t447;
t376 = t448 * pkin(4);
t371 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t229 = -t296 * qJ(3) - t247;
t206 = pkin(3) * t397 - t229;
t350 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t347 = -mrSges(7,1) * t301 + mrSges(7,2) * t297;
t343 = Ifges(7,1) * t297 + t425;
t340 = Ifges(7,2) * t301 + t426;
t337 = Ifges(7,5) * t297 + Ifges(7,6) * t301;
t64 = pkin(11) * t398 + t503;
t159 = -pkin(4) * t244 + t206;
t92 = -pkin(5) * t320 - pkin(11) * t169 + t159;
t32 = t297 * t92 + t301 * t64;
t31 = -t297 * t64 + t301 * t92;
t154 = mrSges(5,1) * t351 - mrSges(5,3) * t174;
t155 = -mrSges(5,2) * t351 + mrSges(5,3) * t175;
t328 = t302 * t154 + t299 * t155;
t178 = -mrSges(5,2) * t263 + t432;
t179 = mrSges(5,1) * t263 - t431;
t327 = t302 * t178 - t299 * t179;
t273 = qJD(2) * t377;
t241 = -pkin(8) * t361 + t273;
t145 = -t169 * t297 + t301 * t398;
t146 = t169 * t301 + t297 * t398;
t65 = t110 * t448 - t298 * t117;
t19 = t110 * t358 - t117 * t380 + t298 * t68 + t448 * t74;
t227 = -pkin(8) * t352 + t262;
t292 = t296 * qJD(3);
t188 = t273 + t292 - t330;
t315 = (-qJ(3) * t385 - t384) * t295;
t242 = t247 * qJD(2);
t196 = -t227 - t264;
t228 = qJD(1) * t242;
t314 = -t227 * mrSges(3,2) - t196 * mrSges(4,3) + t228 * t524;
t313 = -qJD(6) * t334 - t442;
t138 = -pkin(4) * t202 + t188;
t310 = t313 + t443;
t15 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + t90 * Ifges(7,6);
t16 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + t90 * Ifges(7,5);
t257 = Ifges(6,3) * t351;
t8 = -pkin(5) * t351 - t11;
t308 = mrSges(7,3) * t443 + qJD(6) * t513 + t15 * t450 + t16 * t452 + t337 * t478 + t340 * t484 + t343 * t485 + t8 * t347 + t257 + t497;
t307 = -t428 / 0.2e1 + t423 / 0.2e1 + t422 / 0.2e1 + t421 / 0.2e1 - t410 / 0.2e1 + t526;
t290 = -t376 - pkin(5);
t260 = Ifges(3,5) * t351;
t259 = Ifges(4,5) * t352;
t258 = Ifges(5,3) * t351;
t246 = -t283 + t377;
t238 = -qJ(3) * t362 + t267;
t237 = (mrSges(4,2) * t303 - mrSges(4,3) * t300) * t389;
t234 = -mrSges(3,2) * t282 + mrSges(3,3) * t362;
t231 = t296 * t364 + t283;
t224 = -t241 - t292;
t219 = t269 + t315;
t216 = -qJD(1) * t354 + t272;
t194 = qJD(1) * t315 + t261;
t164 = t215 * t301 + t297 * t362;
t163 = -t215 * t297 + t301 * t362;
t162 = -t214 * t301 + t282 * t297;
t161 = t214 * t297 + t282 * t301;
t122 = -mrSges(5,1) * t175 + mrSges(5,2) * t174;
t113 = t174 * Ifges(5,1) + t175 * Ifges(5,4) + Ifges(5,5) * t351;
t112 = t174 * Ifges(5,4) + t175 * Ifges(5,2) + Ifges(5,6) * t351;
t105 = -mrSges(6,1) * t356 + mrSges(6,2) * t321;
t104 = qJD(5) * t169 + t298 * t201 - t202 * t448;
t103 = qJD(5) * t320 + t201 * t448 + t298 * t202;
t91 = t106 + t446;
t81 = -mrSges(6,2) * t351 - mrSges(6,3) * t90;
t70 = -qJD(6) * t146 - t103 * t297 + t301 * t360;
t69 = qJD(6) * t145 + t103 * t301 + t297 * t360;
t63 = -pkin(5) * t398 - t65;
t47 = t101 * t448 - t396;
t46 = t101 * t298 + t363;
t38 = mrSges(6,1) * t90 + mrSges(6,2) * t89;
t35 = Ifges(6,5) * t351 - t435 + t437;
t34 = -t90 * Ifges(6,2) + Ifges(6,6) * t351 + t436;
t33 = pkin(5) * t104 - pkin(11) * t103 + t138;
t29 = t106 * t297 + t301 * t42;
t28 = t106 * t301 - t297 * t42;
t25 = t297 * t91 + t301 * t47;
t24 = -t297 * t47 + t301 * t91;
t18 = -pkin(5) * t360 - t20;
t17 = pkin(11) * t360 + t19;
t5 = -qJD(6) * t32 - t17 * t297 + t301 * t33;
t4 = qJD(6) * t31 + t17 * t301 + t297 * t33;
t1 = [t503 * t81 + m(6) * (t10 * t503 + t11 * t65 + t121 * t159 + t138 * t144 + t19 * t43 + t20 * t42) + (t10 * t320 - t103 * t42 - t104 * t43 - t11 * t169) * mrSges(6,3) + t89 * (Ifges(6,1) * t169 + Ifges(6,4) * t320) / 0.2e1 + t2 * (mrSges(7,2) * t320 + mrSges(7,3) * t145) + t3 * (-mrSges(7,1) * t320 - mrSges(7,3) * t146) + t121 * (-mrSges(6,1) * t320 + mrSges(6,2) * t169) + (Ifges(7,5) * t146 + Ifges(7,6) * t145 - Ifges(7,3) * t320) * t478 + (Ifges(6,4) * t169 + Ifges(6,2) * t320) * t479 + (Ifges(7,4) * t146 + Ifges(7,2) * t145 - Ifges(7,6) * t320) * t484 + (Ifges(7,1) * t146 + Ifges(7,4) * t145 - Ifges(7,5) * t320) * t485 - t320 * t486 + (-t114 * t201 + t115 * t202 + t244 * t75 + t324 * t76) * mrSges(5,3) + t166 * (-mrSges(5,1) * t244 - mrSges(5,2) * t324) + (-Ifges(5,4) * t324 + Ifges(5,2) * t244) * t460 + (-Ifges(5,1) * t324 + Ifges(5,4) * t244) * t461 + ((-mrSges(4,1) * t196 + mrSges(4,2) * t194 + mrSges(3,3) * t227) * t303 + (-t194 * mrSges(4,3) + t257 / 0.2e1 + t258 / 0.2e1 + t525 * t228 + t492 + t497) * t300) * t295 + t356 * (Ifges(6,4) * t103 - Ifges(6,2) * t104) / 0.2e1 + t263 * (Ifges(5,5) * t201 + Ifges(5,6) * t202) / 0.2e1 + (t260 / 0.2e1 + t259 / 0.2e1 + t314) * t296 + t225 * (Ifges(5,4) * t201 + Ifges(5,2) * t202) / 0.2e1 + t253 * (Ifges(6,5) * t103 - Ifges(6,6) * t104) / 0.2e1 + (t300 * t493 + t303 * t494) * t387 + t321 * (Ifges(6,1) * t103 - Ifges(6,4) * t104) / 0.2e1 + t241 * t234 + t219 * t237 + t224 * t235 + t206 * t122 + t201 * t143 / 0.2e1 + t177 * (-mrSges(5,1) * t202 + mrSges(5,2) * t201) + t202 * t142 / 0.2e1 + t188 * t157 + t84 * t178 + t85 * t179 + t132 * t154 + t133 * t155 + t159 * t38 + t145 * t15 / 0.2e1 + t8 * (-mrSges(7,1) * t145 + mrSges(7,2) * t146) + t146 * t16 / 0.2e1 + t138 * t105 + t144 * (mrSges(6,1) * t104 + mrSges(6,2) * t103) + t19 * t134 + t20 * t135 + t22 * (mrSges(7,1) * t104 - mrSges(7,3) * t69) + t23 * (-mrSges(7,2) * t104 + mrSges(7,3) * t70) + t4 * t93 + t5 * t94 + t18 * t83 + t65 * t80 + t40 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t69 * t62 / 0.2e1 + t63 * t21 + t32 * t27 + t31 * t26 + t390 * t242 + m(7) * (t18 * t40 + t2 * t32 + t22 * t5 + t23 * t4 + t3 * t31 + t63 * t8) + m(5) * (t114 * t85 + t115 * t84 + t132 * t76 + t133 * t75 + t166 * t206 + t177 * t188) + m(4) * (t194 * t230 + t196 * t229 + t205 * t242 + t212 * t224 + t219 * t222 + t228 * t231) + m(3) * (t227 * t247 - t228 * t246 + t239 * t242 + t240 * t241) + ((-t230 * mrSges(4,3) + Ifges(5,5) * t455 + Ifges(5,6) * t456 + Ifges(6,5) * t462 + Ifges(6,6) * t463 + t231 * mrSges(4,1) - t246 * mrSges(3,3) + (-Ifges(4,4) + Ifges(3,5) / 0.2e1) * t296 + (t303 * t371 - 0.2e1 * t472) * t295) * t303 + (t229 * mrSges(4,1) - t230 * mrSges(4,2) - t247 * mrSges(3,3) + (-Ifges(3,6) + t488) * t296 + (-t300 * t371 - 0.2e1 * t473) * t295 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t397) * t300) * t357 + t113 * t455 + t112 * t456 + (Ifges(5,1) * t201 + Ifges(5,4) * t202) * t457 + t35 * t462 + t34 * t463 + (Ifges(7,5) * t69 + Ifges(7,6) * t70 + Ifges(7,3) * t104) * t466 + (Ifges(7,1) * t69 + Ifges(7,4) * t70 + Ifges(7,5) * t104) * t468 + (Ifges(7,4) * t69 + Ifges(7,2) * t70 + Ifges(7,6) * t104) * t470 + t103 * t475 + t104 * t477 + t70 * t481 + t104 * t483; -t433 * t319 + t522 * t94 + (t127 * t3 + t128 * t2 + t22 * t522 + t23 * t523 - t319 * t8 + t40 * t521) * m(7) + t523 * t93 + t519 * t134 + (t10 * t204 + t11 * t319 + t121 * t287 + t144 * t518 + t42 * t520 + t43 * t519) * m(6) + t520 * t135 + t521 * t83 + t515 * mrSges(6,3) + t518 * t105 + (-pkin(2) * t228 - qJ(3) * t196 - t205 * t240 + t212 * t514 - t222 * t238) * m(4) + t512 * t215 + t328 * t304 + (t418 / 0.2e1 + t490) * t192 + (Ifges(7,4) * t164 + Ifges(7,2) * t163) * t471 + t260 + t259 + t314 + (Ifges(7,1) * t164 + Ifges(7,4) * t163) * t469 + m(5) * (t166 * qJ(3) + t177 * qJD(3) + t304 * t500) + ((-t266 / 0.2e1 + (t472 - t424 / 0.2e1) * t389 - t494) * t303 + (((t473 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t300) * t295 + (-Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t397) * qJD(1) + t489 - t493) * t300 + (t336 * mrSges(4,1) - Ifges(4,4) * t303 - Ifges(3,6) * t300 + (Ifges(5,5) * t302 - Ifges(6,5) * t251 - Ifges(5,6) * t299 - Ifges(6,6) * t252) * t303 / 0.2e1) * qJD(2)) * t389 + (Ifges(7,5) * t164 + Ifges(7,6) * t163) * t467 + ((-m(5) * t329 + t327) * t304 + t489) * qJD(4) + (-t235 + t234) * t239 - t500 * mrSges(5,3) + (-t8 * t346 + t338 * t479 - t50 * t344 / 0.2e1 - t51 * t341 / 0.2e1 - t35 / 0.2e1 - t121 * mrSges(6,2) - t437 / 0.2e1 + t435 / 0.2e1 - t301 * t16 / 0.2e1 + t15 * t452 + t11 * mrSges(6,3) + (t337 * t466 + t340 * t470 + t343 * t468 + t347 * t40 + t450 * t61 + t452 * t62) * qJD(6)) * t251 - m(5) * (t114 * t139 + t115 * t140 + t177 * t216) + t166 * (mrSges(5,1) * t299 + mrSges(5,2) * t302) + t302 * t113 / 0.2e1 + t287 * t38 - t238 * t237 - t216 * t157 + t204 * t81 - t140 * t178 - t139 * t179 - t40 * (-mrSges(7,1) * t163 + mrSges(7,2) * t164) + t128 * t27 + t127 * t26 + qJ(3) * t122 + (-t23 * t163 + t22 * t164 + (qJD(6) * t333 + t2 * t297 + t3 * t301) * t251) * mrSges(7,3) - t390 * t240 + t391 * qJD(3) + (mrSges(6,1) * t393 - mrSges(6,2) * t392) * t144 + (-t420 / 0.2e1 + t307) * t193 + (-t511 + t526) * t214 + (t49 / 0.2e1 + t48 / 0.2e1 + t486 + t121 * mrSges(6,1) - t436 / 0.2e1 - t34 / 0.2e1 - t10 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + t487) * t90 + t350) * t252 + t112 * t451 + (-Ifges(5,2) * t299 + t429) * t460 + (Ifges(5,1) * t302 - t430) * t461 + t164 * t480 + t163 * t482; t214 * t134 - t161 * t94 - t162 * t93 + t433 * t251 + t327 * qJD(4) + (-t105 - t391) * t282 + t499 * t193 + (mrSges(4,1) * t385 + (t237 + t327) * t300) * t389 + (t81 + (-t297 * t93 - t301 * t94) * qJD(6) + t332) * t252 + t328 + t504 * t392 + (-t161 * t22 - t162 * t23 + t193 * t333 + t251 * t8 + t252 * t310 + t392 * t40) * m(7) + (t10 * t252 - t11 * t251 - t144 * t282 - t515) * m(6) + (-t177 * t282 - t263 * t329 + t500) * m(5) + (t212 * t282 + t222 * t275 + t228) * m(4); -t504 * (-pkin(4) * t380 + t46) + (m(7) * t310 + t517) * (pkin(11) + t445) + (t338 * t467 + t341 * t471 + t344 * t469 - t323 + t495 + t512) * t356 + (-Ifges(5,2) * t226 + t143 + t223) * t459 + (-t60 / 0.2e1 + t98 / 0.2e1 - t491 + t496 + t511) * t321 + t308 + t499 * pkin(4) * t358 + t492 + (-t178 + t432) * t114 + t258 + ((t448 * t11 + t10 * t298 + (-t298 * t42 + t43 * t448) * qJD(5)) * pkin(4) - t144 * t446 + t42 * t46 - t43 * t47) * m(6) + (t179 + t431) * t115 + (t8 * t290 + (t298 * t40 + t333 * t448) * qJD(5) * pkin(4) - t22 * t24 - t23 * t25 - t40 * t46) * m(7) + (t22 * t399 + t23 * t400 + t498) * mrSges(7,3) + t290 * t21 - t177 * (mrSges(5,1) * t226 + mrSges(5,2) * t225) - t47 * t134 - t25 * t93 - t24 * t94 - t105 * t446 + t80 * t376 + t81 * t445 + (Ifges(5,5) * t225 - Ifges(5,6) * t226) * t453 + t142 * t457 + (Ifges(5,1) * t225 - t414) * t458 + t399 * t480 + t400 * t481; ((t487 - Ifges(6,1) / 0.2e1) * t321 - t490 + t495) * t356 + t308 + t313 * mrSges(7,3) - t504 * t43 - t42 * t134 - t29 * t93 - t28 * t94 - pkin(5) * t21 + t517 * pkin(11) + (-t307 + t496) * t321 + (-pkin(5) * t8 - t22 * t28 - t23 * t29 - t40 * t43 + (t443 + t498) * pkin(11)) * m(7); -t40 * (mrSges(7,1) * t130 + mrSges(7,2) * t129) + (Ifges(7,1) * t129 - t427) * t469 + t61 * t468 + (Ifges(7,5) * t129 - Ifges(7,6) * t130) * t467 - t22 * t93 + t23 * t94 + (t129 * t22 + t130 * t23) * mrSges(7,3) + t350 + t14 + (-Ifges(7,2) * t130 + t126 + t62) * t471;];
tauc  = t1(:);
