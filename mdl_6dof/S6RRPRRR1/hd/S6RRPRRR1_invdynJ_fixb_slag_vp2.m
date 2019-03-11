% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:38
% EndTime: 2019-03-09 13:13:12
% DurationCPUTime: 22.85s
% Computational Cost: add. (28115->787), mult. (67681->1039), div. (0->0), fcn. (53246->18), ass. (0->376)
t311 = qJD(2) + qJD(4);
t305 = qJD(5) + t311;
t322 = cos(qJ(5));
t317 = sin(qJ(5));
t315 = -qJ(3) - pkin(7);
t324 = cos(qJ(2));
t276 = t315 * t324;
t262 = qJD(1) * t276;
t313 = sin(pkin(11));
t244 = t313 * t262;
t319 = sin(qJ(2));
t274 = t315 * t319;
t261 = qJD(1) * t274;
t250 = qJD(2) * pkin(2) + t261;
t314 = cos(pkin(11));
t197 = t314 * t250 + t244;
t417 = qJD(1) * t324;
t418 = qJD(1) * t319;
t241 = -t313 * t417 - t314 * t418;
t477 = pkin(8) * t241;
t156 = qJD(2) * pkin(3) + t197 + t477;
t425 = t314 * t262;
t198 = t313 * t250 - t425;
t255 = -t313 * t319 + t314 * t324;
t240 = t255 * qJD(1);
t478 = pkin(8) * t240;
t161 = t198 + t478;
t318 = sin(qJ(4));
t323 = cos(qJ(4));
t105 = t156 * t318 + t161 * t323;
t372 = t323 * t240 + t241 * t318;
t538 = pkin(9) * t372;
t94 = t105 + t538;
t441 = t317 * t94;
t104 = t323 * t156 - t161 * t318;
t190 = t240 * t318 - t241 * t323;
t553 = pkin(9) * t190;
t93 = t104 - t553;
t91 = pkin(4) * t311 + t93;
t52 = t322 * t91 - t441;
t47 = -pkin(5) * t305 - t52;
t316 = sin(qJ(6));
t321 = cos(qJ(6));
t544 = t322 * t190 + t317 * t372;
t121 = t305 * t321 - t316 * t544;
t122 = t305 * t316 + t321 * t544;
t552 = -mrSges(6,1) * t305 - mrSges(7,1) * t121 + mrSges(7,2) * t122 + mrSges(6,3) * t544;
t557 = -m(6) * t52 + m(7) * t47 + t552;
t312 = qJ(2) + pkin(11);
t306 = qJ(4) + t312;
t293 = qJ(5) + t306;
t286 = sin(t293);
t287 = cos(t293);
t466 = mrSges(7,2) * t316;
t556 = t286 * t466 + t287 * (m(7) * pkin(10) + mrSges(7,3));
t545 = -t190 * t317 + t322 * t372;
t132 = -t545 + qJD(6);
t88 = pkin(5) * t544 - pkin(10) * t545;
t407 = qJD(1) * qJD(2);
t383 = t319 * t407;
t406 = qJDD(1) * t324;
t263 = -t383 + t406;
t264 = qJDD(1) * t319 + t324 * t407;
t205 = t263 * t314 - t264 * t313;
t206 = t263 * t313 + t264 * t314;
t113 = qJD(4) * t372 + t205 * t318 + t206 * t323;
t114 = -qJD(4) * t190 + t205 * t323 - t206 * t318;
t437 = t322 * t94;
t53 = t317 * t91 + t437;
t48 = pkin(10) * t305 + t53;
t308 = t324 * pkin(2);
t296 = t308 + pkin(1);
t271 = -qJD(1) * t296 + qJD(3);
t207 = -pkin(3) * t240 + t271;
t147 = -pkin(4) * t372 + t207;
t74 = -pkin(5) * t545 - pkin(10) * t544 + t147;
t22 = -t316 * t48 + t321 * t74;
t23 = t316 * t74 + t321 * t48;
t61 = qJD(5) * t545 + t113 * t322 + t114 * t317;
t62 = -qJD(5) * t544 - t113 * t317 + t114 * t322;
t434 = qJDD(1) * pkin(1);
t233 = -pkin(2) * t263 + qJDD(3) - t434;
t162 = -pkin(3) * t205 + t233;
t90 = -pkin(4) * t114 + t162;
t19 = -pkin(5) * t62 - pkin(10) * t61 + t90;
t309 = qJDD(2) + qJDD(4);
t254 = t264 * pkin(7);
t415 = qJD(3) * t319;
t195 = qJDD(2) * pkin(2) - qJ(3) * t264 - qJD(1) * t415 - t254;
t297 = pkin(7) * t406;
t416 = qJD(2) * t319;
t399 = pkin(7) * t416;
t414 = qJD(3) * t324;
t201 = qJ(3) * t263 + t297 + (-t399 + t414) * qJD(1);
t141 = t314 * t195 - t201 * t313;
t107 = qJDD(2) * pkin(3) - pkin(8) * t206 + t141;
t142 = t313 * t195 + t314 * t201;
t115 = pkin(8) * t205 + t142;
t45 = -qJD(4) * t105 + t323 * t107 - t115 * t318;
t34 = pkin(4) * t309 - pkin(9) * t113 + t45;
t412 = qJD(4) * t323;
t413 = qJD(4) * t318;
t44 = t318 * t107 + t323 * t115 + t156 * t412 - t161 * t413;
t36 = pkin(9) * t114 + t44;
t410 = qJD(5) * t322;
t411 = qJD(5) * t317;
t10 = t317 * t34 + t322 * t36 + t91 * t410 - t411 * t94;
t302 = qJDD(5) + t309;
t7 = pkin(10) * t302 + t10;
t3 = -qJD(6) * t23 + t19 * t321 - t316 * t7;
t11 = -qJD(5) * t53 - t317 * t36 + t322 * t34;
t41 = qJD(6) * t121 + t302 * t316 + t321 * t61;
t42 = -qJD(6) * t122 + t302 * t321 - t316 * t61;
t60 = qJDD(6) - t62;
t14 = Ifges(7,4) * t41 + Ifges(7,2) * t42 + Ifges(7,6) * t60;
t15 = Ifges(7,1) * t41 + Ifges(7,4) * t42 + Ifges(7,5) * t60;
t2 = qJD(6) * t22 + t19 * t316 + t321 * t7;
t358 = Ifges(7,5) * t321 - Ifges(7,6) * t316;
t342 = t132 * t358;
t458 = Ifges(7,4) * t316;
t362 = Ifges(7,1) * t321 - t458;
t343 = t122 * t362;
t457 = Ifges(7,4) * t321;
t360 = -Ifges(7,2) * t316 + t457;
t344 = t121 * t360;
t363 = mrSges(7,1) * t316 + mrSges(7,2) * t321;
t345 = t47 * t363;
t409 = qJD(6) * t316;
t381 = -t409 / 0.2e1;
t120 = Ifges(7,4) * t121;
t71 = t122 * Ifges(7,1) + t132 * Ifges(7,5) + t120;
t439 = t321 * t71;
t389 = t439 / 0.2e1;
t464 = mrSges(7,3) * t321;
t488 = t316 / 0.2e1;
t506 = t60 / 0.2e1;
t507 = t42 / 0.2e1;
t508 = t41 / 0.2e1;
t468 = mrSges(7,1) * t321;
t548 = t466 - t468;
t449 = t122 * Ifges(7,4);
t70 = t121 * Ifges(7,2) + t132 * Ifges(7,6) + t449;
t8 = -pkin(5) * t302 - t11;
t327 = -t10 * mrSges(6,2) + t2 * t464 + t15 * t488 + t321 * t14 / 0.2e1 + Ifges(6,3) * t302 + (Ifges(7,1) * t316 + t457) * t508 + (Ifges(7,2) * t321 + t458) * t507 + (Ifges(7,5) * t316 + Ifges(7,6) * t321) * t506 + Ifges(6,6) * t62 + Ifges(6,5) * t61 + t8 * t548 + t70 * t381 + t11 * mrSges(6,1) + (t345 + t389) * qJD(6) + (t344 + t343 + t342) * qJD(6) / 0.2e1;
t465 = mrSges(7,3) * t316;
t490 = -t305 / 0.2e1;
t496 = -t544 / 0.2e1;
t497 = -t545 / 0.2e1;
t498 = -t132 / 0.2e1;
t500 = -t122 / 0.2e1;
t501 = -t121 / 0.2e1;
t471 = t53 * mrSges(6,3);
t536 = t23 * mrSges(7,2);
t537 = t22 * mrSges(7,1);
t452 = Ifges(7,3) * t132;
t453 = Ifges(7,6) * t121;
t455 = Ifges(7,5) * t122;
t69 = t452 + t453 + t455;
t454 = Ifges(6,6) * t305;
t459 = Ifges(6,4) * t544;
t85 = Ifges(6,2) * t545 + t454 + t459;
t513 = -t147 * mrSges(6,1) + t471 + t85 / 0.2e1 - t69 / 0.2e1 + t536 - t537;
t472 = t52 * mrSges(6,3);
t131 = Ifges(6,4) * t545;
t456 = Ifges(6,5) * t305;
t86 = Ifges(6,1) * t544 + t131 + t456;
t514 = -t147 * mrSges(6,2) - t345 - t439 / 0.2e1 + t70 * t488 + t472 - t86 / 0.2e1;
t408 = qJD(6) * t321;
t525 = -t22 * t408 - t23 * t409;
t555 = t45 * mrSges(5,1) - t44 * mrSges(5,2) + mrSges(7,3) * t525 + Ifges(5,5) * t113 + Ifges(5,6) * t114 + Ifges(5,3) * t309 - t3 * t465 + (Ifges(6,1) * t496 + Ifges(6,4) * t497 + Ifges(6,5) * t490 + t22 * t464 + t23 * t465 + t358 * t498 + t360 * t501 + t362 * t500 + t514) * t545 + (-Ifges(6,4) * t496 + Ifges(7,5) * t500 - Ifges(6,2) * t497 - Ifges(6,6) * t490 + Ifges(7,6) * t501 + Ifges(7,3) * t498 + t513) * t544 + t327;
t527 = t287 * pkin(5) + t286 * pkin(10);
t554 = m(7) * t527;
t484 = pkin(4) * t190;
t275 = -t324 * mrSges(3,1) + t319 * mrSges(3,2);
t303 = sin(t312);
t304 = cos(t312);
t551 = -t304 * mrSges(4,1) + t303 * mrSges(4,2) + t275;
t550 = -t287 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t286;
t290 = sin(t306);
t291 = cos(t306);
t467 = mrSges(6,2) * t287;
t549 = mrSges(5,1) * t290 + mrSges(6,1) * t286 + mrSges(5,2) * t291 + t467;
t320 = sin(qJ(1));
t325 = cos(qJ(1));
t547 = g(1) * t325 + g(2) * t320;
t541 = m(6) + m(7);
t540 = t263 / 0.2e1;
t539 = t264 / 0.2e1;
t470 = qJD(2) / 0.2e1;
t256 = t313 * t324 + t314 * t319;
t535 = Ifges(4,5) * t256;
t534 = Ifges(4,6) * t255;
t436 = qJDD(2) / 0.2e1;
t203 = -t261 * t313 + t425;
t163 = t203 - t478;
t204 = t314 * t261 + t244;
t164 = t204 + t477;
t116 = t323 * t163 - t164 * t318;
t486 = pkin(2) * t314;
t292 = pkin(3) + t486;
t487 = pkin(2) * t313;
t237 = t292 * t318 + t323 * t487;
t223 = t237 * qJD(4);
t530 = -t116 - t223;
t117 = t318 * t163 + t323 * t164;
t236 = t323 * t292 - t318 * t487;
t222 = t236 * qJD(4);
t529 = -t117 + t222;
t123 = -mrSges(6,2) * t305 + mrSges(6,3) * t545;
t83 = -mrSges(7,2) * t132 + mrSges(7,3) * t121;
t84 = mrSges(7,1) * t132 - mrSges(7,3) * t122;
t355 = -t316 * t84 + t321 * t83;
t528 = -t123 - t355;
t208 = t314 * t274 + t276 * t313;
t180 = -pkin(8) * t256 + t208;
t209 = t313 * t274 - t314 * t276;
t181 = pkin(8) * t255 + t209;
t128 = t318 * t180 + t323 * t181;
t234 = pkin(4) + t236;
t174 = t317 * t234 + t322 * t237;
t376 = t291 * mrSges(5,1) - t290 * mrSges(5,2);
t526 = t287 * t548 + t550;
t199 = t255 * t323 - t256 * t318;
t200 = t255 * t318 + t256 * t323;
t146 = t199 * t317 + t200 * t322;
t242 = t256 * qJD(2);
t243 = t255 * qJD(2);
t143 = qJD(4) * t199 - t242 * t318 + t243 * t323;
t144 = -qJD(4) * t200 - t242 * t323 - t243 * t318;
t351 = t322 * t199 - t200 * t317;
t76 = qJD(5) * t351 + t143 * t322 + t144 * t317;
t347 = t146 * t408 + t316 * t76;
t253 = -pkin(7) * t383 + t297;
t524 = t253 * t324 + t254 * t319;
t521 = -t376 + t526;
t20 = mrSges(7,1) * t60 - mrSges(7,3) * t41;
t21 = -mrSges(7,2) * t60 + mrSges(7,3) * t42;
t519 = -t316 * t20 + t321 * t21 - t84 * t408 - t83 * t409;
t518 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t310 = -pkin(8) + t315;
t517 = -m(3) * pkin(7) + m(4) * t315 + m(5) * t310 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t419 = pkin(3) * t304 + t308;
t516 = mrSges(2,1) + m(5) * (pkin(1) + t419) + t376 + m(4) * t296 + m(3) * pkin(1) - t550 - t551;
t357 = t22 * t321 + t23 * t316;
t473 = t3 * t316;
t332 = -qJD(6) * t357 - t473;
t474 = t2 * t321;
t515 = m(7) * (t332 + t474) + t519;
t509 = m(7) * pkin(5);
t499 = t122 / 0.2e1;
t495 = -t372 / 0.2e1;
t494 = -t190 / 0.2e1;
t493 = t190 / 0.2e1;
t491 = -t241 / 0.2e1;
t489 = -t311 / 0.2e1;
t485 = pkin(2) * t319;
t483 = pkin(4) * t290;
t285 = pkin(4) * t291;
t482 = pkin(4) * t317;
t481 = pkin(4) * t322;
t480 = pkin(5) * t286;
t479 = pkin(7) * t324;
t463 = Ifges(3,4) * t319;
t462 = Ifges(3,4) * t324;
t461 = Ifges(4,4) * t241;
t460 = Ifges(5,4) * t190;
t451 = t104 * mrSges(5,3);
t450 = t105 * mrSges(5,3);
t448 = t197 * mrSges(4,3);
t447 = t198 * mrSges(4,3);
t433 = t146 * t316;
t432 = t146 * t321;
t424 = t316 * t325;
t423 = t320 * t316;
t422 = t320 * t321;
t421 = t321 * t325;
t380 = qJD(2) * t315;
t238 = t319 * t380 + t414;
t239 = t324 * t380 - t415;
t179 = t314 * t238 + t313 * t239;
t405 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t60;
t403 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t418) * t479;
t300 = pkin(2) * t416;
t398 = t286 * t468;
t299 = pkin(2) * t418;
t388 = t285 + t419;
t384 = -t62 * mrSges(6,1) + t61 * mrSges(6,2);
t210 = -pkin(3) * t241 + t299;
t211 = pkin(3) * t242 + t300;
t378 = -t205 * mrSges(4,1) + t206 * mrSges(4,2);
t377 = -t114 * mrSges(5,1) + t113 * mrSges(5,2);
t127 = t323 * t180 - t181 * t318;
t178 = -t238 * t313 + t314 * t239;
t217 = -pkin(3) * t255 - t296;
t371 = t556 * t320;
t370 = t556 * t325;
t270 = -pkin(3) * t303 - t485;
t368 = mrSges(3,1) * t319 + mrSges(3,2) * t324;
t361 = t324 * Ifges(3,2) + t463;
t359 = Ifges(3,5) * t324 - Ifges(3,6) * t319;
t356 = -t22 * t316 + t23 * t321;
t98 = -pkin(9) * t200 + t127;
t99 = pkin(9) * t199 + t128;
t73 = t317 * t98 + t322 * t99;
t157 = -pkin(4) * t199 + t217;
t81 = -pkin(5) * t351 - pkin(10) * t146 + t157;
t32 = t316 * t81 + t321 * t73;
t31 = -t316 * t73 + t321 * t81;
t72 = t317 * t99 - t322 * t98;
t173 = t234 * t322 - t237 * t317;
t119 = -pkin(4) * t144 + t211;
t148 = t210 + t484;
t349 = t116 - t538;
t348 = pkin(1) * t368;
t346 = t146 * t409 - t321 * t76;
t341 = t319 * (Ifges(3,1) * t324 - t463);
t153 = -pkin(8) * t243 + t178;
t154 = -pkin(8) * t242 + t179;
t78 = t318 * t153 + t323 * t154 + t180 * t412 - t181 * t413;
t230 = t270 - t483;
t333 = m(7) * (t230 - t480) - t398;
t331 = m(7) * (-t480 - t483) - t398;
t330 = t467 + (mrSges(6,1) + t468 + t509) * t286;
t79 = -qJD(4) * t128 + t323 * t153 - t154 * t318;
t328 = -pkin(9) * t143 + t79;
t307 = -pkin(9) + t310;
t298 = Ifges(3,4) * t417;
t295 = -pkin(5) - t481;
t273 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t417;
t249 = Ifges(3,1) * t418 + Ifges(3,5) * qJD(2) + t298;
t248 = Ifges(3,6) * qJD(2) + qJD(1) * t361;
t235 = Ifges(4,4) * t240;
t228 = t287 * t421 + t423;
t227 = -t287 * t424 + t422;
t226 = -t287 * t422 + t424;
t225 = t287 * t423 + t421;
t224 = pkin(1) + t388;
t215 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t241;
t214 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t240;
t196 = -mrSges(4,1) * t240 - mrSges(4,2) * t241;
t192 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t206;
t191 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t205;
t186 = -t241 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t235;
t185 = t240 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t461;
t182 = Ifges(5,4) * t372;
t171 = -pkin(5) - t173;
t166 = mrSges(5,1) * t311 - mrSges(5,3) * t190;
t165 = -mrSges(5,2) * t311 + mrSges(5,3) * t372;
t139 = -mrSges(5,1) * t372 + mrSges(5,2) * t190;
t130 = t190 * Ifges(5,1) + t311 * Ifges(5,5) + t182;
t129 = Ifges(5,2) * t372 + t311 * Ifges(5,6) + t460;
t101 = -mrSges(5,2) * t309 + mrSges(5,3) * t114;
t100 = mrSges(5,1) * t309 - mrSges(5,3) * t113;
t96 = t117 - t553;
t87 = -mrSges(6,1) * t545 + mrSges(6,2) * t544;
t80 = t484 + t88;
t77 = qJD(5) * t146 + t143 * t317 - t322 * t144;
t75 = t148 + t88;
t68 = pkin(9) * t144 + t78;
t64 = t317 * t349 + t322 * t96;
t56 = t322 * t93 - t441;
t55 = t317 * t93 + t437;
t51 = -mrSges(6,2) * t302 + mrSges(6,3) * t62;
t50 = mrSges(6,1) * t302 - mrSges(6,3) * t61;
t30 = t316 * t88 + t321 * t52;
t29 = -t316 * t52 + t321 * t88;
t28 = t316 * t80 + t321 * t56;
t27 = -t316 * t56 + t321 * t80;
t26 = t316 * t75 + t321 * t64;
t25 = -t316 * t64 + t321 * t75;
t24 = pkin(5) * t77 - pkin(10) * t76 + t119;
t17 = -qJD(5) * t72 + t317 * t328 + t322 * t68;
t16 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t5 = -qJD(6) * t32 - t17 * t316 + t24 * t321;
t4 = qJD(6) * t31 + t17 * t321 + t24 * t316;
t1 = [t557 * (qJD(5) * t73 + t317 * t68 - t322 * t328) - (Ifges(7,3) * t506 + Ifges(7,6) * t507 + Ifges(7,5) * t508 - Ifges(6,4) * t61 + t90 * mrSges(6,1) - Ifges(6,2) * t62 + t405 / 0.2e1 - Ifges(6,6) * t302 - t10 * mrSges(6,3) + t518) * t351 + (-t226 * mrSges(7,1) - t225 * mrSges(7,2) + (t307 * t541 + t517) * t325 + (-m(7) * (-t224 - t527) + m(6) * t224 + t516) * t320) * g(1) + (-m(6) * t11 + m(7) * t8 + t16 - t50) * t72 + m(7) * (t2 * t32 + t22 * t5 + t23 * t4 + t3 * t31) + m(6) * (t10 * t73 + t119 * t147 + t157 * t90 + t17 * t53) + (t359 * t470 - t403) * qJD(2) + (t534 / 0.2e1 - mrSges(3,2) * t479 + Ifges(3,6) * t324 / 0.2e1 + t535 / 0.2e1) * qJDD(2) + (t534 + t535) * t436 - t347 * t70 / 0.2e1 + (t263 * t479 + t524) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t524) + t77 * t537 + t462 * t539 + t361 * t540 + (-t141 * t256 + t142 * t255) * mrSges(4,3) + (mrSges(5,2) * t162 - mrSges(5,3) * t45 + Ifges(5,1) * t113 + Ifges(5,4) * t114 + Ifges(5,5) * t309) * t200 - t143 * t451 + m(5) * (t104 * t79 + t105 * t78 + t127 * t45 + t128 * t44 + t162 * t217 + t207 * t211) + t157 * t384 - t275 * t434 - t348 * t407 + t15 * t432 / 0.2e1 - t14 * t433 / 0.2e1 - t77 * t471 - t76 * t472 - t248 * t416 / 0.2e1 + t78 * t165 + t79 * t166 + t147 * (mrSges(6,1) * t77 + mrSges(6,2) * t76) + (-t2 * t433 + t22 * t346 - t23 * t347 - t3 * t432) * mrSges(7,3) - t243 * t448 + t144 * t129 / 0.2e1 + t217 * t377 - t296 * t378 + m(4) * (t141 * t208 + t142 * t209 + t178 * t197 + t179 * t198 - t233 * t296 + t271 * t300) - t273 * t399 + t76 * t389 + (-mrSges(5,1) * t162 + mrSges(5,3) * t44 + Ifges(5,4) * t113 + Ifges(5,2) * t114 + Ifges(5,6) * t309) * t199 + (t90 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t61 + Ifges(6,4) * t62 + Ifges(6,5) * t302 + t358 * t506 + t360 * t507 + t362 * t508 + t363 * t8 + t381 * t71) * t146 + t196 * t300 + (-t228 * mrSges(7,1) - t227 * mrSges(7,2) - t541 * (t325 * t224 - t320 * t307) + t517 * t320 + (-t516 - t554) * t325) * g(2) + t143 * t130 / 0.2e1 + t17 * t123 + t127 * t100 + t128 * t101 + t119 * t87 + Ifges(2,3) * qJDD(1) + t4 * t83 + t5 * t84 - t77 * t85 / 0.2e1 + t76 * t86 / 0.2e1 + t77 * t69 / 0.2e1 + t73 * t51 + t240 * (Ifges(4,4) * t243 - Ifges(4,2) * t242) / 0.2e1 + t271 * (mrSges(4,1) * t242 + mrSges(4,2) * t243) + (Ifges(4,5) * t243 - Ifges(4,6) * t242) * t470 + (Ifges(4,1) * t243 - Ifges(4,4) * t242) * t491 + t372 * (Ifges(5,4) * t143 + Ifges(5,2) * t144) / 0.2e1 + (t341 + t324 * (-Ifges(3,2) * t319 + t462)) * t407 / 0.2e1 + t31 * t20 + t32 * t21 + t47 * (mrSges(7,1) * t347 - mrSges(7,2) * t346) + t121 * (-Ifges(7,4) * t346 - Ifges(7,2) * t347 + Ifges(7,6) * t77) / 0.2e1 + t132 * (-Ifges(7,5) * t346 - Ifges(7,6) * t347 + Ifges(7,3) * t77) / 0.2e1 - t77 * t536 + (Ifges(3,4) * t539 + Ifges(3,2) * t540 + Ifges(3,6) * t436 + t249 * t470) * t324 + (Ifges(3,1) * t264 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t264) + Ifges(3,4) * t540 + 0.2e1 * Ifges(3,5) * t436) * t319 + t207 * (-mrSges(5,1) * t144 + mrSges(5,2) * t143) + t208 * t192 + t209 * t191 + t211 * t139 + t179 * t214 + t178 * t215 + t544 * (Ifges(6,1) * t76 - Ifges(6,4) * t77) / 0.2e1 + t545 * (Ifges(6,4) * t76 - Ifges(6,2) * t77) / 0.2e1 - t242 * t185 / 0.2e1 + t243 * t186 / 0.2e1 + t205 * (Ifges(4,4) * t256 + Ifges(4,2) * t255) + t206 * (Ifges(4,1) * t256 + Ifges(4,4) * t255) + t233 * (-mrSges(4,1) * t255 + mrSges(4,2) * t256) - pkin(1) * (-mrSges(3,1) * t263 + mrSges(3,2) * t264) - t242 * t447 + t144 * t450 + (Ifges(5,1) * t143 + Ifges(5,4) * t144) * t493 + (-Ifges(7,1) * t346 - Ifges(7,4) * t347 + Ifges(7,5) * t77) * t499 + t305 * (Ifges(6,5) * t76 - Ifges(6,6) * t77) / 0.2e1 + t311 * (Ifges(5,5) * t143 + Ifges(5,6) * t144) / 0.2e1; -t557 * (-qJD(5) * t174 + (-t223 - t349) * t322 + (-t222 + t96) * t317) + t515 * (pkin(10) + t174) - (-t129 / 0.2e1 - t450 + t207 * mrSges(5,1) + Ifges(5,6) * t489 + Ifges(5,4) * t494 + Ifges(5,2) * t495) * t190 + (t403 + (t348 - t341 / 0.2e1) * qJD(1)) * qJD(1) + (t10 * t174 + t11 * t173 - t147 * t148 - t53 * t64) * m(6) + (t171 * t8 - t22 * t25 - t23 * t26) * m(7) + (-m(5) * t419 - m(7) * (t388 + t527) - m(4) * t308 - m(6) * t388 + t521 + t551) * g(3) + t547 * (m(4) * t485 - m(5) * t270 - m(6) * t230 + mrSges(4,1) * t303 + mrSges(4,2) * t304 + t368 + t549) + t529 * t165 + (t104 * t530 + t105 * t529 - t207 * t210 + t236 * t45 + t237 * t44) * m(5) + t530 * t166 + (t248 / 0.2e1 + pkin(7) * t273) * t418 + t241 * (Ifges(4,1) * t240 + t461) / 0.2e1 - t359 * t407 / 0.2e1 + t174 * t51 + t171 * t16 + t173 * t50 - t148 * t87 - t196 * t299 - t241 * t447 + t555 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t141 * mrSges(4,1) - t142 * mrSges(4,2) - t64 * t123 - t26 * t83 - t25 * t84 + (-t130 / 0.2e1 - t207 * mrSges(5,2) + t451 + Ifges(5,5) * t489 + Ifges(5,1) * t494 + Ifges(5,4) * t495) * t372 - (-Ifges(3,2) * t418 + t249 + t298) * t417 / 0.2e1 - (Ifges(4,2) * t241 + t186 + t235) * t240 / 0.2e1 - g(1) * (t325 * t333 + t370) - g(2) * (t320 * t333 + t371) + Ifges(4,6) * t205 + Ifges(4,5) * t206 - t210 * t139 - t204 * t214 - t203 * t215 + (m(6) * t53 + m(7) * t356 - t528) * (qJD(5) * t173 + t222 * t322 - t223 * t317) + t236 * t100 + t237 * t101 - qJD(2) * (Ifges(4,5) * t240 + Ifges(4,6) * t241) / 0.2e1 - t253 * mrSges(3,2) - t254 * mrSges(3,1) + Ifges(3,6) * t263 + Ifges(3,5) * t264 - t271 * (-mrSges(4,1) * t241 + mrSges(4,2) * t240) + t240 * t448 + t192 * t486 + t191 * t487 + t185 * t491 + (-t197 * t203 - t198 * t204 - t271 * t299 + (t141 * t314 + t142 * t313) * pkin(2)) * m(4); t190 * t166 - t372 * t165 + t378 + t384 + t377 - t552 * t544 + t528 * t545 + t355 * qJD(6) - t240 * t214 - t241 * t215 + t316 * t21 + t321 * t20 + (-g(1) * t320 + g(2) * t325) * (m(4) + m(5) + t541) + (t132 * t356 + t2 * t316 + t3 * t321 - t544 * t47) * m(7) + (t52 * t544 - t53 * t545 + t90) * m(6) + (t104 * t190 - t105 * t372 + t162) * m(5) + (-t197 * t241 - t198 * t240 + t233) * m(4); t515 * (pkin(10) + t482) - t528 * pkin(4) * t410 + (-m(7) * (t285 + t527) - m(6) * t285 + t521) * g(3) + (m(6) * t483 + t549) * t547 + (-t22 * t27 - t23 * t28 - t47 * t55 + t295 * t8 + (t317 * t47 + t322 * t356) * qJD(5) * pkin(4)) * m(7) + ((t10 * t317 + t11 * t322 + (-t317 * t52 + t322 * t53) * qJD(5)) * pkin(4) - t147 * t484 + t52 * t55 - t53 * t56) * m(6) + (-Ifges(5,2) * t190 + t130 + t182) * t495 - t104 * t165 + t105 * t166 + t552 * (pkin(4) * t411 - t55) + t555 - t56 * t123 - t28 * t83 - t27 * t84 - t207 * (mrSges(5,1) * t190 + mrSges(5,2) * t372) + t372 * t451 + (Ifges(5,5) * t372 - Ifges(5,6) * t190) * t489 + (Ifges(5,1) * t372 - t460) * t494 - g(1) * (t325 * t331 + t370) - g(2) * (t320 * t331 + t371) + t190 * t450 + t50 * t481 + t51 * t482 + t129 * t493 + t295 * t16 - t87 * t484; (t526 - t554) * g(3) + (t325 * t330 - t370) * g(1) + (t320 * t330 - t371) * g(2) + t327 - m(7) * (t22 * t29 + t23 * t30 + t47 * t53) + (-t131 / 0.2e1 - t456 / 0.2e1 - t344 / 0.2e1 - t343 / 0.2e1 - t342 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t544 + t357 * mrSges(7,3) + t514) * t545 + (t454 / 0.2e1 - t453 / 0.2e1 - t455 / 0.2e1 - t452 / 0.2e1 + t459 / 0.2e1 + t513) * t544 - t552 * t53 + (m(7) * (-t473 + t474 + t525) + t519) * pkin(10) - t52 * t123 - t30 * t83 - t29 * t84 - pkin(5) * t16 - t8 * t509 + t332 * mrSges(7,3); -t47 * (mrSges(7,1) * t122 + mrSges(7,2) * t121) + (Ifges(7,1) * t121 - t449) * t500 + t70 * t499 + (Ifges(7,5) * t121 - Ifges(7,6) * t122) * t498 - t22 * t83 + t23 * t84 - g(1) * (mrSges(7,1) * t227 - mrSges(7,2) * t228) - g(2) * (-mrSges(7,1) * t225 + mrSges(7,2) * t226) + g(3) * t363 * t286 + (t121 * t22 + t122 * t23) * mrSges(7,3) + t405 + (-Ifges(7,2) * t122 + t120 + t71) * t501 + t518;];
tau  = t1;
