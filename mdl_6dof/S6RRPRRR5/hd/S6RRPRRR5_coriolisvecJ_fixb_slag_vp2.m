% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:15
% EndTime: 2019-03-09 13:39:18
% DurationCPUTime: 35.25s
% Computational Cost: add. (26169->914), mult. (76952->1306), div. (0->0), fcn. (62419->12), ass. (0->393)
t317 = sin(pkin(6));
t316 = sin(pkin(12));
t322 = sin(qJ(2));
t326 = cos(qJ(2));
t411 = cos(pkin(12));
t365 = t411 * t326;
t335 = -t316 * t322 + t365;
t331 = t317 * t335;
t264 = qJD(1) * t331;
t496 = -t264 + qJD(4);
t318 = cos(pkin(6));
t446 = pkin(1) * t318;
t311 = t326 * t446;
t306 = qJD(1) * t311;
t442 = pkin(8) + qJ(3);
t373 = t442 * t322;
t359 = t317 * t373;
t248 = -qJD(1) * t359 + t306;
t310 = t322 * t446;
t403 = t317 * t326;
t493 = t403 * t442 + t310;
t249 = t493 * qJD(1);
t367 = t411 * t249;
t193 = t248 * t316 + t367;
t321 = sin(qJ(4));
t325 = cos(qJ(4));
t535 = -t193 + t496 * (pkin(4) * t321 - pkin(10) * t325);
t366 = t411 * t322;
t392 = qJD(1) * t317;
t379 = t326 * t392;
t265 = -t316 * t379 - t366 * t392;
t324 = cos(qJ(5));
t320 = sin(qJ(5));
t401 = t320 * t325;
t215 = -t264 * t401 - t265 * t324;
t388 = qJD(4) * t325;
t534 = t320 * t388 + t215;
t385 = qJD(5) * t324;
t501 = t321 * t385 + t534;
t237 = t316 * t249;
t194 = t248 * t411 - t237;
t380 = t322 * t392;
t363 = pkin(2) * t380;
t210 = -pkin(3) * t265 - pkin(9) * t264 + t363;
t125 = t325 * t194 + t321 * t210;
t108 = -pkin(10) * t265 + t125;
t381 = t411 * pkin(2);
t314 = -t381 - pkin(3);
t289 = -t325 * pkin(4) - t321 * pkin(10) + t314;
t445 = pkin(2) * t316;
t313 = pkin(9) + t445;
t387 = qJD(5) * t320;
t389 = qJD(4) * t321;
t509 = -t324 * t108 + t289 * t385 + (-t324 * t389 - t325 * t387) * t313 + t535 * t320;
t378 = t313 * t389;
t533 = t535 * t324 + (t108 + t378) * t320;
t308 = qJD(1) * t318 + qJD(2);
t229 = pkin(2) * t308 + t248;
t167 = t316 * t229 + t367;
t159 = pkin(9) * t308 + t167;
t293 = (-pkin(2) * t326 - pkin(1)) * t317;
t284 = qJD(1) * t293 + qJD(3);
t187 = -t264 * pkin(3) + t265 * pkin(9) + t284;
t106 = t159 * t325 + t187 * t321;
t166 = t229 * t411 - t237;
t158 = -t308 * pkin(3) - t166;
t227 = t321 * t265 + t308 * t325;
t228 = -t265 * t325 + t308 * t321;
t102 = -t227 * pkin(4) - t228 * pkin(10) + t158;
t262 = -t335 * t392 + qJD(4);
t93 = pkin(10) * t262 + t106;
t48 = t324 * t102 - t320 * t93;
t49 = t102 * t320 + t324 * t93;
t532 = -t48 * mrSges(6,1) + t49 * mrSges(6,2) + t106 * mrSges(5,3);
t399 = t324 * t325;
t216 = t264 * t399 - t265 * t320;
t294 = t313 * t399;
t408 = t264 * t321;
t531 = -pkin(5) * t408 + pkin(11) * t216 + (pkin(5) * t321 - pkin(11) * t399) * qJD(4) + (-t294 + (pkin(11) * t321 - t289) * t320) * qJD(5) + t533;
t530 = -t501 * pkin(11) + t509;
t476 = -pkin(11) - pkin(10);
t382 = qJD(5) * t476;
t409 = t227 * t320;
t105 = -t321 * t159 + t187 * t325;
t164 = pkin(4) * t228 - pkin(10) * t227;
t77 = t324 * t105 + t320 * t164;
t529 = pkin(11) * t409 + t320 * t382 - t77;
t76 = -t105 * t320 + t324 * t164;
t528 = -pkin(5) * t228 - t76 + (pkin(11) * t227 + t382) * t324;
t419 = t262 * Ifges(5,6);
t438 = Ifges(5,4) * t228;
t139 = t227 * Ifges(5,2) + t419 + t438;
t527 = t139 / 0.2e1 + t532;
t177 = -t228 * t320 + t262 * t324;
t178 = t228 * t324 + t262 * t320;
t319 = sin(qJ(6));
t323 = cos(qJ(6));
t118 = t177 * t319 + t178 * t323;
t225 = qJD(5) - t227;
t42 = -pkin(11) * t178 + t48;
t35 = pkin(5) * t225 + t42;
t43 = pkin(11) * t177 + t49;
t417 = t319 * t43;
t16 = t323 * t35 - t417;
t414 = t323 * t43;
t17 = t319 * t35 + t414;
t271 = (t316 * t326 + t366) * t317;
t266 = qJD(2) * t271;
t259 = qJD(1) * t266;
t301 = qJD(2) * t306;
t330 = (-qJD(2) * t373 + qJD(3) * t326) * t317;
t220 = qJD(1) * t330 + t301;
t404 = t317 * t322;
t231 = -qJD(2) * t493 - qJD(3) * t404;
t328 = qJD(1) * t231;
t146 = t220 * t411 + t316 * t328;
t267 = qJD(2) * t331;
t260 = qJD(1) * t267;
t390 = qJD(2) * t317;
t372 = qJD(1) * t390;
t358 = t322 * t372;
t342 = pkin(2) * t358;
t188 = pkin(3) * t259 - pkin(9) * t260 + t342;
t58 = t325 * t146 - t159 * t389 + t187 * t388 + t321 * t188;
t50 = pkin(10) * t259 + t58;
t145 = t220 * t316 - t411 * t328;
t175 = t227 * qJD(4) + t260 * t325;
t176 = qJD(4) * t228 + t260 * t321;
t80 = pkin(4) * t176 - pkin(10) * t175 + t145;
t15 = -qJD(5) * t49 - t320 * t50 + t324 * t80;
t96 = qJD(5) * t177 + t175 * t324 + t259 * t320;
t6 = pkin(5) * t176 - pkin(11) * t96 + t15;
t14 = t102 * t385 + t320 * t80 + t324 * t50 - t387 * t93;
t97 = -qJD(5) * t178 - t175 * t320 + t259 * t324;
t7 = pkin(11) * t97 + t14;
t2 = qJD(6) * t16 + t319 * t6 + t323 * t7;
t3 = -qJD(6) * t17 - t319 * t7 + t323 * t6;
t364 = t323 * t177 - t178 * t319;
t432 = Ifges(7,4) * t118;
t221 = qJD(6) + t225;
t459 = -t221 / 0.2e1;
t472 = -t118 / 0.2e1;
t92 = -pkin(4) * t262 - t105;
t73 = -pkin(5) * t177 + t92;
t33 = qJD(6) * t364 + t319 * t97 + t323 * t96;
t34 = -qJD(6) * t118 - t319 * t96 + t323 * t97;
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t176;
t526 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t8 + (Ifges(7,5) * t364 - Ifges(7,6) * t118) * t459 + (t118 * t17 + t16 * t364) * mrSges(7,3) - t73 * (mrSges(7,1) * t118 + mrSges(7,2) * t364) + (Ifges(7,1) * t364 - t432) * t472;
t422 = t221 * Ifges(7,3);
t427 = t118 * Ifges(7,5);
t428 = t364 * Ifges(7,6);
t53 = t422 + t427 + t428;
t421 = t225 * Ifges(6,3);
t423 = t178 * Ifges(6,5);
t424 = t177 * Ifges(6,6);
t89 = t421 + t423 + t424;
t515 = t89 + t53;
t525 = t515 / 0.2e1;
t124 = -t321 * t194 + t210 * t325;
t107 = pkin(4) * t265 - t124;
t377 = t313 * t388;
t524 = -t107 + t377;
t114 = Ifges(7,4) * t364;
t523 = -Ifges(7,2) * t118 + t114;
t487 = t33 / 0.2e1;
t486 = t34 / 0.2e1;
t478 = t96 / 0.2e1;
t477 = t97 / 0.2e1;
t465 = t176 / 0.2e1;
t39 = Ifges(6,5) * t96 + Ifges(6,6) * t97 + Ifges(6,3) * t176;
t521 = t39 + t8;
t520 = t14 * mrSges(6,2);
t519 = t15 * mrSges(6,1);
t518 = Ifges(6,3) + Ifges(7,3);
t282 = t324 * t289;
t400 = t321 * t324;
t219 = -pkin(11) * t400 + t282 + (-t313 * t320 - pkin(5)) * t325;
t236 = t320 * t289 + t294;
t402 = t320 * t321;
t226 = -pkin(11) * t402 + t236;
t147 = t219 * t323 - t226 * t319;
t517 = qJD(6) * t147 + t319 * t531 + t530 * t323;
t148 = t219 * t319 + t226 * t323;
t516 = -qJD(6) * t148 - t530 * t319 + t323 * t531;
t141 = mrSges(5,1) * t259 - mrSges(5,3) * t175;
t46 = -mrSges(6,1) * t97 + mrSges(6,2) * t96;
t510 = -t141 + t46;
t508 = -qJD(5) * t236 + t533;
t303 = t476 * t320;
t304 = t476 * t324;
t256 = t303 * t323 + t304 * t319;
t507 = qJD(6) * t256 + t319 * t528 + t323 * t529;
t257 = t303 * t319 - t304 * t323;
t506 = -qJD(6) * t257 - t319 * t529 + t323 * t528;
t505 = pkin(5) * t501 + t524;
t63 = -mrSges(7,1) * t364 + mrSges(7,2) * t118;
t504 = m(7) * t73 + t63;
t339 = t319 * t320 - t323 * t324;
t280 = t339 * t321;
t245 = pkin(2) * t318 + t311 - t359;
t393 = pkin(8) * t403 + t310;
t263 = qJ(3) * t403 + t393;
t206 = t316 * t245 + t411 * t263;
t190 = pkin(9) * t318 + t206;
t270 = t316 * t404 - t317 * t365;
t214 = t270 * pkin(3) - t271 * pkin(9) + t293;
t128 = t325 * t190 + t321 * t214;
t113 = pkin(10) * t270 + t128;
t205 = t245 * t411 - t316 * t263;
t189 = -t318 * pkin(3) - t205;
t242 = t271 * t321 - t318 * t325;
t243 = t271 * t325 + t318 * t321;
t126 = t242 * pkin(4) - t243 * pkin(10) + t189;
t65 = t324 * t113 + t320 * t126;
t134 = t215 * t319 + t216 * t323;
t292 = t319 * t324 + t320 * t323;
t495 = qJD(5) + qJD(6);
t241 = t495 * t292;
t195 = -t241 * t321 - t339 * t388;
t503 = t195 - t134;
t133 = t215 * t323 - t216 * t319;
t196 = t280 * t495 - t292 * t388;
t502 = t196 - t133;
t386 = qJD(5) * t321;
t500 = t320 * t386 - t324 * t388 + t216;
t440 = mrSges(4,3) * t265;
t499 = mrSges(4,1) * t308 + mrSges(5,1) * t227 - mrSges(5,2) * t228 + t440;
t498 = t389 - t408;
t59 = -t321 * t146 - t159 * t388 - t187 * t389 + t188 * t325;
t497 = -t321 * t59 + t325 * t58;
t345 = t14 * t324 - t15 * t320;
t344 = t320 * t49 + t324 * t48;
t347 = Ifges(6,5) * t324 - Ifges(6,6) * t320;
t433 = Ifges(6,4) * t324;
t350 = -Ifges(6,2) * t320 + t433;
t434 = Ifges(6,4) * t320;
t353 = Ifges(6,1) * t324 - t434;
t355 = mrSges(6,1) * t320 + mrSges(6,2) * t324;
t172 = Ifges(6,4) * t177;
t91 = Ifges(6,1) * t178 + Ifges(6,5) * t225 + t172;
t413 = t324 * t91;
t456 = t225 / 0.2e1;
t461 = t178 / 0.2e1;
t463 = t177 / 0.2e1;
t435 = Ifges(6,4) * t178;
t90 = Ifges(6,2) * t177 + Ifges(6,6) * t225 + t435;
t479 = -t90 / 0.2e1;
t494 = -t344 * mrSges(6,3) + t320 * t479 + t413 / 0.2e1 + t347 * t456 + t350 * t463 + t353 * t461 + t92 * t355;
t490 = -0.2e1 * pkin(1);
t489 = Ifges(7,4) * t487 + Ifges(7,2) * t486 + Ifges(7,6) * t465;
t488 = Ifges(7,1) * t487 + Ifges(7,4) * t486 + Ifges(7,5) * t465;
t40 = Ifges(6,4) * t96 + Ifges(6,2) * t97 + Ifges(6,6) * t176;
t485 = t40 / 0.2e1;
t484 = Ifges(6,1) * t478 + Ifges(6,4) * t477 + Ifges(6,5) * t465;
t54 = Ifges(7,2) * t364 + Ifges(7,6) * t221 + t432;
t483 = -t54 / 0.2e1;
t482 = t54 / 0.2e1;
t55 = Ifges(7,1) * t118 + Ifges(7,5) * t221 + t114;
t481 = -t55 / 0.2e1;
t480 = t55 / 0.2e1;
t474 = -t364 / 0.2e1;
t473 = t364 / 0.2e1;
t471 = t118 / 0.2e1;
t470 = -t139 / 0.2e1;
t223 = Ifges(5,4) * t227;
t420 = t262 * Ifges(5,5);
t140 = t228 * Ifges(5,1) + t223 + t420;
t468 = -t140 / 0.2e1;
t467 = t175 / 0.2e1;
t466 = -t176 / 0.2e1;
t464 = -t177 / 0.2e1;
t462 = -t178 / 0.2e1;
t458 = t221 / 0.2e1;
t457 = -t225 / 0.2e1;
t455 = -t242 / 0.2e1;
t453 = t243 / 0.2e1;
t451 = -t265 / 0.2e1;
t449 = t318 / 0.2e1;
t448 = -t322 / 0.2e1;
t444 = pkin(5) * t320;
t441 = mrSges(4,3) * t264;
t439 = Ifges(3,4) * t322;
t437 = Ifges(5,4) * t321;
t436 = Ifges(5,4) * t325;
t431 = Ifges(3,5) * t326;
t430 = t105 * mrSges(5,3);
t418 = t265 * Ifges(4,4);
t51 = -pkin(4) * t259 - t59;
t416 = t321 * t51;
t407 = t264 * t325;
t149 = t292 * t227;
t398 = -t149 + t241;
t150 = t339 * t227;
t240 = t495 * t339;
t397 = -t150 + t240;
t120 = -mrSges(6,1) * t177 + mrSges(6,2) * t178;
t186 = mrSges(5,1) * t262 - mrSges(5,3) * t228;
t396 = t186 - t120;
t383 = Ifges(5,5) * t175 - Ifges(5,6) * t176 + Ifges(5,3) * t259;
t375 = t322 * t390;
t64 = -t113 * t320 + t324 * t126;
t127 = -t321 * t190 + t214 * t325;
t307 = qJD(2) * t311;
t230 = t307 + t330;
t162 = t230 * t316 - t411 * t231;
t362 = pkin(2) * t375;
t361 = mrSges(3,3) * t380;
t360 = mrSges(3,3) * t379;
t354 = Ifges(5,1) * t325 - t437;
t352 = Ifges(6,1) * t320 + t433;
t351 = -Ifges(5,2) * t321 + t436;
t349 = Ifges(6,2) * t324 + t434;
t348 = Ifges(5,5) * t325 - Ifges(5,6) * t321;
t346 = Ifges(6,5) * t320 + Ifges(6,6) * t324;
t213 = t243 * t324 + t270 * t320;
t44 = pkin(5) * t242 - pkin(11) * t213 + t64;
t212 = -t243 * t320 + t270 * t324;
t52 = pkin(11) * t212 + t65;
t22 = -t319 * t52 + t323 * t44;
t23 = t319 * t44 + t323 * t52;
t74 = mrSges(6,1) * t176 - mrSges(6,3) * t96;
t75 = -mrSges(6,2) * t176 + mrSges(6,3) * t97;
t343 = -t320 * t74 + t324 * t75;
t136 = -mrSges(6,2) * t225 + mrSges(6,3) * t177;
t137 = mrSges(6,1) * t225 - mrSges(6,3) * t178;
t340 = -t320 * t136 - t324 * t137;
t129 = t212 * t323 - t213 * t319;
t130 = t212 * t319 + t213 * t323;
t163 = t230 * t411 + t316 * t231;
t211 = pkin(3) * t266 - pkin(9) * t267 + t362;
t70 = -t321 * t163 - t190 * t388 + t211 * t325 - t214 * t389;
t112 = -pkin(4) * t270 - t127;
t336 = t308 * (-Ifges(3,6) * t322 + t431);
t69 = t325 * t163 - t190 * t389 + t321 * t211 + t214 * t388;
t61 = pkin(10) * t266 + t69;
t203 = qJD(4) * t243 + t267 * t321;
t204 = -qJD(4) * t242 + t267 * t325;
t85 = pkin(4) * t203 - pkin(10) * t204 + t162;
t20 = -t113 * t387 + t126 * t385 + t320 * t85 + t324 * t61;
t62 = -pkin(4) * t266 - t70;
t278 = t393 * qJD(2);
t268 = -pkin(8) * t358 + t301;
t269 = qJD(1) * t278;
t329 = -t269 * mrSges(3,1) - t145 * mrSges(4,1) - t268 * mrSges(3,2) - t146 * mrSges(4,2);
t21 = -qJD(5) * t65 - t320 * t61 + t324 * t85;
t315 = -pkin(5) * t324 - pkin(4);
t305 = Ifges(3,4) * t379;
t300 = t372 * t431;
t287 = -pkin(8) * t404 + t311;
t283 = (t313 + t444) * t321;
t279 = t292 * t321;
t277 = -pkin(8) * t375 + t307;
t276 = t393 * qJD(1);
t275 = -pkin(8) * t380 + t306;
t274 = -t308 * mrSges(3,2) + t360;
t273 = mrSges(3,1) * t308 - t361;
t261 = Ifges(4,4) * t264;
t253 = Ifges(4,5) * t260;
t252 = Ifges(4,6) * t259;
t250 = t260 * mrSges(4,2);
t247 = Ifges(3,1) * t380 + t308 * Ifges(3,5) + t305;
t246 = Ifges(3,6) * t308 + (Ifges(3,2) * t326 + t439) * t392;
t235 = -t313 * t401 + t282;
t232 = -mrSges(4,2) * t308 + t441;
t217 = -mrSges(4,1) * t264 - mrSges(4,2) * t265;
t209 = -t265 * Ifges(4,1) + t308 * Ifges(4,5) + t261;
t208 = t264 * Ifges(4,2) + t308 * Ifges(4,6) - t418;
t185 = -mrSges(5,2) * t262 + mrSges(5,3) * t227;
t142 = -mrSges(5,2) * t259 - mrSges(5,3) * t176;
t138 = t228 * Ifges(5,5) + t227 * Ifges(5,6) + t262 * Ifges(5,3);
t119 = mrSges(5,1) * t176 + mrSges(5,2) * t175;
t111 = qJD(5) * t212 + t204 * t324 + t266 * t320;
t110 = -qJD(5) * t213 - t204 * t320 + t266 * t324;
t99 = t175 * Ifges(5,1) - t176 * Ifges(5,4) + t259 * Ifges(5,5);
t98 = t175 * Ifges(5,4) - t176 * Ifges(5,2) + t259 * Ifges(5,6);
t88 = mrSges(7,1) * t221 - mrSges(7,3) * t118;
t87 = -mrSges(7,2) * t221 + mrSges(7,3) * t364;
t86 = pkin(5) * t409 + t106;
t82 = -pkin(5) * t212 + t112;
t38 = -qJD(6) * t130 + t110 * t323 - t111 * t319;
t37 = qJD(6) * t129 + t110 * t319 + t111 * t323;
t36 = -pkin(5) * t110 + t62;
t30 = -pkin(5) * t97 + t51;
t29 = -mrSges(7,2) * t176 + mrSges(7,3) * t34;
t28 = mrSges(7,1) * t176 - mrSges(7,3) * t33;
t19 = t323 * t42 - t417;
t18 = -t319 * t42 - t414;
t13 = pkin(11) * t110 + t20;
t12 = pkin(5) * t203 - pkin(11) * t111 + t21;
t11 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t5 = -qJD(6) * t23 + t12 * t323 - t13 * t319;
t4 = qJD(6) * t22 + t12 * t319 + t13 * t323;
t1 = [m(3) * (t268 * t393 - t269 * t287 - t275 * t278 + t276 * t277) - t499 * t162 + t203 * t525 + t3 * (mrSges(7,1) * t242 - mrSges(7,3) * t130) + t145 * (mrSges(5,1) * t242 + mrSges(5,2) * t243) + ((-t287 * mrSges(3,3) + Ifges(3,5) * t449 + (mrSges(3,2) * t490 + 0.3e1 / 0.2e1 * Ifges(3,4) * t326) * t317) * t326 + (-t393 * mrSges(3,3) - Ifges(3,6) * t318 + (mrSges(3,1) * t490 - 0.3e1 / 0.2e1 * t439 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t326) * t317 + (m(4) * t293 + mrSges(4,1) * t270 + mrSges(4,2) * t271) * pkin(2)) * t322) * t372 + m(5) * (t105 * t70 + t106 * t69 + t127 * t59 + t128 * t58 + t145 * t189 + t158 * t162) + m(6) * (t112 * t51 + t14 * t65 + t15 * t64 + t20 * t49 + t21 * t48 + t62 * t92) + m(7) * (t16 * t5 + t17 * t4 + t2 * t23 + t22 * t3 + t30 * t82 + t36 * t73) + t14 * (-mrSges(6,2) * t242 + mrSges(6,3) * t212) + t15 * (mrSges(6,1) * t242 - mrSges(6,3) * t213) + t2 * (-mrSges(7,2) * t242 + mrSges(7,3) * t129) + t163 * t232 + (Ifges(6,5) * t213 + Ifges(7,5) * t130 + Ifges(6,6) * t212 + Ifges(7,6) * t129 + t242 * t518) * t465 + (t322 * pkin(2) * t217 + t246 * t448 + t336 / 0.2e1 + t326 * t247 / 0.2e1 + (-t275 * t326 - t276 * t322) * mrSges(3,3)) * t390 + (-t205 * mrSges(4,3) + Ifges(4,1) * t271 - Ifges(4,4) * t270 + Ifges(4,5) * t449) * t260 + (-t206 * mrSges(4,3) + Ifges(5,5) * t453 + Ifges(5,6) * t455 - Ifges(4,4) * t271 + t293 * mrSges(4,1) - Ifges(4,6) * t318 / 0.2e1 + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t270) * t259 + t51 * (-mrSges(6,1) * t212 + mrSges(6,2) * t213) + t158 * (mrSges(5,1) * t203 + mrSges(5,2) * t204) + t204 * t140 / 0.2e1 + (t145 * t271 - t146 * t270 - t166 * t267 - t167 * t266) * mrSges(4,3) + t521 * t242 / 0.2e1 + t49 * (-mrSges(6,2) * t203 + mrSges(6,3) * t110) + t48 * (mrSges(6,1) * t203 - mrSges(6,3) * t111) + t16 * (mrSges(7,1) * t203 - mrSges(7,3) * t37) + t17 * (-mrSges(7,2) * t203 + mrSges(7,3) * t38) + t70 * t186 + t189 * t119 + t69 * t185 + t20 * t136 + t21 * t137 + t127 * t141 + t128 * t142 + t30 * (-mrSges(7,1) * t129 + mrSges(7,2) * t130) + t62 * t120 + t110 * t90 / 0.2e1 + t92 * (-mrSges(6,1) * t110 + mrSges(6,2) * t111) + t111 * t91 / 0.2e1 + t112 * t46 + t270 * t383 / 0.2e1 + t5 * t88 + t4 * t87 + t82 * t11 + t73 * (-mrSges(7,1) * t38 + mrSges(7,2) * t37) + t64 * t74 + t65 * t75 + t36 * t63 + m(4) * (-t145 * t205 + t146 * t206 - t162 * t166 + t163 * t167 + t284 * t362) + t22 * t28 + t23 * t29 + (t268 * t326 + t269 * t322) * t317 * mrSges(3,3) + (t300 / 0.2e1 + t253 / 0.2e1 - t252 / 0.2e1 + t329) * t318 + t227 * (Ifges(5,4) * t204 - Ifges(5,2) * t203 + Ifges(5,6) * t266) / 0.2e1 + t228 * (Ifges(5,1) * t204 - Ifges(5,4) * t203 + Ifges(5,5) * t266) / 0.2e1 + t106 * (-mrSges(5,2) * t266 - mrSges(5,3) * t203) + t105 * (mrSges(5,1) * t266 - mrSges(5,3) * t204) + t266 * t138 / 0.2e1 + t262 * (Ifges(5,5) * t204 - Ifges(5,6) * t203 + Ifges(5,3) * t266) / 0.2e1 - t266 * t208 / 0.2e1 + t267 * t209 / 0.2e1 + t264 * (Ifges(4,4) * t267 - Ifges(4,2) * t266) / 0.2e1 + t58 * (-mrSges(5,2) * t270 - mrSges(5,3) * t242) + t59 * (mrSges(5,1) * t270 - mrSges(5,3) * t243) + t277 * t274 - t278 * t273 + t284 * (mrSges(4,1) * t266 + mrSges(4,2) * t267) + t308 * (Ifges(4,5) * t267 - Ifges(4,6) * t266) / 0.2e1 + (Ifges(4,1) * t267 - Ifges(4,4) * t266) * t451 + t99 * t453 + t98 * t455 + (Ifges(6,5) * t111 + Ifges(6,6) * t110 + Ifges(6,3) * t203) * t456 + (Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t203) * t458 + (Ifges(6,1) * t111 + Ifges(6,4) * t110 + Ifges(6,5) * t203) * t461 + (Ifges(6,4) * t111 + Ifges(6,2) * t110 + Ifges(6,6) * t203) * t463 + (Ifges(5,4) * t243 - Ifges(5,2) * t242 + Ifges(5,6) * t270) * t466 + (Ifges(5,1) * t243 - Ifges(5,4) * t242 + Ifges(5,5) * t270) * t467 + t203 * t470 + (Ifges(7,1) * t37 + Ifges(7,4) * t38 + Ifges(7,5) * t203) * t471 + t293 * t250 + (Ifges(7,4) * t37 + Ifges(7,2) * t38 + Ifges(7,6) * t203) * t473 + (Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t242) * t477 + (Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t242) * t478 + t37 * t480 + t38 * t482 + t213 * t484 + t212 * t485 + (Ifges(7,4) * t130 + Ifges(7,2) * t129 + Ifges(7,6) * t242) * t486 + (Ifges(7,1) * t130 + Ifges(7,4) * t129 + Ifges(7,5) * t242) * t487 + t130 * t488 + t129 * t489; t534 * t479 + (Ifges(7,5) * t471 + Ifges(7,6) * t473 + Ifges(7,3) * t458 + t470 + t525 - t532) * t389 + (t105 * t407 + t497) * mrSges(5,3) - ((-Ifges(3,2) * t380 + t247 + t305) * t326 + t336) * t392 / 0.2e1 + t300 + (Ifges(6,5) * t216 + Ifges(6,6) * t215) * t457 + (-t125 - t378) * t185 + (-t107 * t92 + t14 * t236 + t15 * t235 + t508 * t48 + t509 * t49) * m(6) + (-t105 * t124 - t106 * t125 + t145 * t314 - t158 * t193) * m(5) + ((t388 * t92 + t416) * m(6) + ((-t105 * t325 - t106 * t321) * qJD(4) + t497) * m(5) + t510 * t321 + t325 * t142) * t313 + (Ifges(6,4) * t216 + Ifges(6,2) * t215) * t464 + (-t515 / 0.2e1 + Ifges(6,5) * t462 + Ifges(7,5) * t472 + Ifges(6,6) * t464 + Ifges(7,6) * t474 + Ifges(6,3) * t457 + Ifges(7,3) * t459 + t527) * t408 + t325 * t520 + t505 * t63 - t388 * t430 + (t361 + t273) * t276 + (t360 - t274) * t275 - t167 * t440 + t508 * t137 + t509 * t136 - (t320 * t91 + t324 * t90) * t386 / 0.2e1 + (t140 + t413) * t388 / 0.2e1 + (t227 * t351 + t228 * t354 + t262 * t348) * qJD(4) / 0.2e1 + (Ifges(4,1) * t264 + t138 + t418) * t265 / 0.2e1 - (Ifges(4,2) * t265 + t209 + t261) * t264 / 0.2e1 + (Ifges(7,4) * t195 + Ifges(7,2) * t196) * t473 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t474 + (Ifges(7,1) * t195 + Ifges(7,4) * t196) * t471 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t472 - t325 * t519 + t253 - t252 + (Ifges(7,5) * t195 + Ifges(7,6) * t196) * t458 + (Ifges(7,5) * t134 + Ifges(7,6) * t133) * t459 + t2 * (mrSges(7,2) * t325 - mrSges(7,3) * t279) + t30 * (mrSges(7,1) * t279 - mrSges(7,2) * t280) + (-Ifges(7,4) * t280 - Ifges(7,2) * t279 - Ifges(7,6) * t325) * t486 + (-t259 * t445 - t260 * t381) * mrSges(4,3) + (-t124 - t377) * t186 + t235 * t74 + t236 * t75 - t194 * t232 + t499 * t193 + (mrSges(6,1) * t501 - mrSges(6,2) * t500) * t92 + (-t14 * t402 - t15 * t400 + t48 * t500 - t49 * t501) * mrSges(6,3) + (-mrSges(7,2) * t498 + mrSges(7,3) * t502) * t17 + (-mrSges(7,1) * t502 + mrSges(7,2) * t503) * t73 + (mrSges(7,1) * t498 - mrSges(7,3) * t503) * t16 - t216 * t91 / 0.2e1 + (-Ifges(7,1) * t280 - Ifges(7,4) * t279 - Ifges(7,5) * t325) * t487 + t3 * (-mrSges(7,1) * t325 + mrSges(7,3) * t280) + t516 * t88 + (t147 * t3 + t148 * t2 + t16 * t516 + t17 * t517 + t283 * t30 + t505 * t73) * m(7) + t517 * t87 + (-Ifges(7,5) * t280 - Ifges(7,6) * t279 + t321 * t347 - t325 * t518) * t465 + t524 * t120 - t521 * t325 / 0.2e1 + t496 * t158 * (mrSges(5,1) * t321 + mrSges(5,2) * t325) + t329 - t40 * t402 / 0.2e1 + t147 * t28 + t148 * t29 + t246 * t380 / 0.2e1 - t217 * t363 - t106 * mrSges(5,2) * t265 + t105 * mrSges(5,1) * t265 - Ifges(3,6) * t358 + (t166 * t193 - t167 * t194 - t284 * t363 + (-t145 * t411 + t146 * t316) * pkin(2)) * m(4) + ((Ifges(3,1) * t326 - t439) * t448 + pkin(1) * (mrSges(3,1) * t322 + mrSges(3,2) * t326)) * qJD(1) ^ 2 * t317 ^ 2 + t283 * t11 - t284 * (-mrSges(4,1) * t265 + mrSges(4,2) * t264) - t308 * (Ifges(4,5) * t264 + Ifges(4,6) * t265) / 0.2e1 + t314 * t119 + t355 * t416 + t166 * t441 + t208 * t451 + (-t346 * t386 + (Ifges(6,3) * t321 + t325 * t347) * qJD(4)) * t456 + (-t352 * t386 + (Ifges(6,5) * t321 + t325 * t353) * qJD(4)) * t461 + (-t349 * t386 + (Ifges(6,6) * t321 + t325 * t350) * qJD(4)) * t463 + (Ifges(5,2) * t325 + t437) * t466 + (Ifges(5,1) * t321 + t436) * t467 + t407 * t468 + t321 * t99 / 0.2e1 + t325 * t98 / 0.2e1 + t259 * (Ifges(5,5) * t321 + Ifges(5,6) * t325) / 0.2e1 + t145 * (-mrSges(5,1) * t325 + mrSges(5,2) * t321) + (Ifges(6,1) * t216 + Ifges(6,4) * t215) * t462 + (-Ifges(6,6) * t325 + t321 * t350) * t477 + (-Ifges(6,5) * t325 + t321 * t353) * t478 + t195 * t480 + t134 * t481 + t196 * t482 + t133 * t483 + t400 * t484 - t280 * t488 - t279 * t489 - t262 * (-Ifges(5,3) * t265 + t264 * t348) / 0.2e1 - t227 * (-Ifges(5,6) * t265 + t264 * t351) / 0.2e1 - t228 * (-Ifges(5,5) * t265 + t264 * t354) / 0.2e1; t259 * mrSges(4,1) - t216 * t136 - t215 * t137 - t264 * t232 - t279 * t28 - t280 * t29 + t250 + t502 * t88 + t503 * t87 - t499 * t265 + (-t264 * t185 - t11 + (t136 * t324 - t137 * t320 + t185) * qJD(4) - t510) * t325 + (qJD(5) * t340 + t142 + t343 + t496 * (t63 - t396)) * t321 + (t16 * t502 + t17 * t503 - t2 * t280 - t279 * t3 - t30 * t325 + t498 * t73) * m(7) + (-t215 * t48 - t216 * t49 - t92 * t408 + (-t51 + (-t320 * t48 + t324 * t49) * qJD(4)) * t325 + (qJD(4) * t92 - qJD(5) * t344 + t345) * t321) * m(6) + (t158 * t265 + t321 * t58 + t325 * t59 + t496 * (-t105 * t321 + t106 * t325)) * m(5) + (-t166 * t265 - t167 * t264 + t342) * m(4); (t16 * t397 - t17 * t398 - t2 * t339 - t292 * t3) * mrSges(7,3) + t30 * (mrSges(7,1) * t339 + mrSges(7,2) * t292) + (Ifges(7,5) * t292 - Ifges(7,6) * t339 + t346) * t465 + (Ifges(7,4) * t292 - Ifges(7,2) * t339) * t486 + (Ifges(7,1) * t292 - Ifges(7,4) * t339) * t487 - t339 * t489 + (-t240 / 0.2e1 + t150 / 0.2e1) * t55 + (-Ifges(7,5) * t240 - Ifges(7,6) * t241) * t458 + (-Ifges(7,1) * t240 - Ifges(7,4) * t241) * t471 + (-Ifges(7,4) * t240 - Ifges(7,2) * t241) * t473 + (-t241 / 0.2e1 + t149 / 0.2e1) * t54 + (t419 / 0.2e1 - t421 / 0.2e1 - t424 / 0.2e1 - t423 / 0.2e1 - t422 / 0.2e1 - t428 / 0.2e1 - t427 / 0.2e1 - t89 / 0.2e1 - t53 / 0.2e1 + t17 * mrSges(7,2) - t16 * mrSges(7,1) - t158 * mrSges(5,1) + t438 / 0.2e1 + t527) * t228 + (-pkin(4) * t51 - t106 * t92 - t48 * t76 - t49 * t77) * m(6) + t256 * t28 + t257 * t29 + (t444 * t504 + t494) * qJD(5) + t506 * t88 + (t16 * t506 + t17 * t507 + t2 * t257 + t256 * t3 + t30 * t315 - t73 * t86) * m(7) + t507 * t87 + (-Ifges(7,5) * t150 - Ifges(7,6) * t149) * t459 + (-Ifges(7,1) * t150 - Ifges(7,4) * t149) * t472 + (-Ifges(7,4) * t150 - Ifges(7,2) * t149) * t474 + (m(6) * t345 + (-m(6) * t344 + t340) * qJD(5) + t343) * pkin(10) + t383 + (-t223 / 0.2e1 - t158 * mrSges(5,2) + t468 + t430 - t420 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t228 - t494) * t227 - t105 * t185 + (mrSges(7,1) * t398 - mrSges(7,2) * t397) * t73 + t396 * t106 - t76 * t137 - t77 * t136 - t86 * t63 - t58 * mrSges(5,2) + t59 * mrSges(5,1) - pkin(4) * t46 + t315 * t11 + t51 * (-mrSges(6,1) * t324 + mrSges(6,2) * t320) + t349 * t477 + t352 * t478 + t320 * t484 + t324 * t485 + t292 * t488 + t345 * mrSges(6,3); -t118 * t483 + (-Ifges(6,2) * t178 + t172 + t91) * t464 + t39 + (t323 * t28 + t319 * t29 + m(7) * (t2 * t319 + t3 * t323) - t504 * t178 + (-t319 * t88 + t323 * t87 + m(7) * (-t16 * t319 + t17 * t323)) * qJD(6)) * pkin(5) + t364 * t481 - m(7) * (t16 * t18 + t17 * t19) - t92 * (mrSges(6,1) * t178 + mrSges(6,2) * t177) + t49 * t137 - t48 * t136 - t18 * t88 - t19 * t87 + (t177 * t48 + t178 * t49) * mrSges(6,3) + t519 - t520 + (Ifges(6,5) * t177 - Ifges(6,6) * t178) * t457 + t90 * t461 + (Ifges(6,1) * t177 - t435) * t462 + t523 * t474 + t526; t54 * t471 - t16 * t87 + t17 * t88 + (t523 + t55) * t474 + t526;];
tauc  = t1(:);
