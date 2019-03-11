% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:58
% EndTime: 2019-03-09 15:25:28
% DurationCPUTime: 18.28s
% Computational Cost: add. (13495->732), mult. (31293->916), div. (0->0), fcn. (22918->14), ass. (0->351)
t533 = mrSges(6,2) - mrSges(5,1);
t311 = sin(qJ(6));
t315 = cos(qJ(6));
t523 = mrSges(7,1) * t311 + mrSges(7,2) * t315;
t312 = sin(qJ(3));
t313 = sin(qJ(2));
t316 = cos(qJ(3));
t317 = cos(qJ(2));
t247 = -t312 * t313 + t316 * t317;
t228 = t247 * qJD(1);
t248 = t312 * t317 + t313 * t316;
t229 = t248 * qJD(1);
t309 = sin(pkin(10));
t310 = cos(pkin(10));
t339 = t228 * t309 + t310 * t229;
t469 = t339 / 0.2e1;
t181 = t228 * t310 - t229 * t309;
t473 = t181 / 0.2e1;
t532 = -Ifges(6,2) * t469 - Ifges(6,6) * t473;
t307 = qJD(2) + qJD(3);
t464 = -t307 / 0.2e1;
t531 = -mrSges(6,1) - mrSges(5,3);
t472 = -t181 / 0.2e1;
t449 = pkin(5) * t181;
t529 = Ifges(5,5) - Ifges(6,4);
t528 = -Ifges(5,6) + Ifges(6,5);
t162 = -t181 * t311 + t307 * t315;
t174 = qJD(6) + t339;
t511 = Ifges(7,3) * t174;
t161 = -t181 * t315 - t307 * t311;
t512 = Ifges(7,6) * t161;
t527 = Ifges(5,1) * t339 + Ifges(5,4) * t181 + Ifges(5,5) * t307 + Ifges(7,5) * t162 + t511 + t512;
t319 = -pkin(8) - pkin(7);
t275 = t319 * t317;
t254 = qJD(1) * t275;
t233 = t316 * t254;
t274 = t319 * t313;
t253 = qJD(1) * t274;
t193 = -t253 * t312 + t233;
t424 = qJ(4) * t228;
t163 = t193 - t424;
t230 = t312 * t254;
t194 = t316 * t253 + t230;
t222 = t229 * qJ(4);
t164 = -t222 + t194;
t107 = t163 * t309 + t164 * t310;
t391 = qJD(3) * t316;
t379 = pkin(2) * t391;
t392 = qJD(3) * t312;
t380 = pkin(2) * t392;
t219 = -t309 * t380 + t310 * t379;
t208 = qJD(5) + t219;
t526 = -t208 + t107;
t308 = qJ(2) + qJ(3);
t301 = pkin(10) + t308;
t287 = sin(t301);
t288 = cos(t301);
t302 = sin(t308);
t303 = cos(t308);
t525 = mrSges(4,1) * t302 + mrSges(5,1) * t287 + mrSges(4,2) * t303 + mrSges(5,2) * t288;
t387 = qJD(1) * qJD(2);
t257 = qJDD(1) * t317 - t313 * t387;
t524 = t523 * t288;
t314 = sin(qJ(1));
t318 = cos(qJ(1));
t495 = g(1) * t318 + g(2) * t314;
t350 = mrSges(7,1) * t315 - mrSges(7,2) * t311;
t237 = qJD(2) * pkin(2) + t253;
t186 = t316 * t237 + t230;
t155 = t186 - t222;
t144 = pkin(3) * t307 + t155;
t187 = t237 * t312 - t233;
t156 = t187 + t424;
t411 = t310 * t156;
t87 = t309 * t144 + t411;
t84 = -qJ(5) * t307 - t87;
t58 = -t84 + t449;
t431 = Ifges(7,4) * t162;
t66 = Ifges(7,2) * t161 + Ifges(7,6) * t174 + t431;
t522 = t350 * t58 - t315 * t66 / 0.2e1;
t521 = -t303 * mrSges(4,1) + mrSges(4,2) * t302 + t533 * t288 + (mrSges(5,2) - mrSges(6,3)) * t287;
t304 = t317 * pkin(2);
t292 = t304 + pkin(1);
t273 = t292 * qJD(1);
t199 = -pkin(3) * t228 + qJD(4) - t273;
t323 = -qJ(5) * t339 + t199;
t95 = -pkin(4) * t181 + t323;
t520 = t199 * mrSges(5,1) + mrSges(6,1) * t84 - t95 * mrSges(6,2) - mrSges(5,3) * t87;
t147 = t309 * t156;
t86 = t144 * t310 - t147;
t355 = qJD(5) - t86;
t450 = pkin(5) * t339;
t480 = pkin(4) + pkin(9);
t54 = -t307 * t480 + t355 + t450;
t64 = -t181 * t480 + t323;
t21 = -t311 * t64 + t315 * t54;
t22 = t311 * t54 + t315 * t64;
t83 = -pkin(4) * t307 + t355;
t519 = mrSges(6,1) * t83 + t21 * mrSges(7,1) + t199 * mrSges(5,2) - t22 * mrSges(7,2) - mrSges(5,3) * t86 - t95 * mrSges(6,3) + Ifges(6,4) * t464 - t532;
t258 = qJDD(1) * t313 + t317 * t387;
t331 = t247 * qJD(3);
t165 = qJD(1) * t331 + t257 * t312 + t258 * t316;
t332 = t248 * qJD(3);
t166 = -qJD(1) * t332 + t257 * t316 - t258 * t312;
t108 = t165 * t309 - t310 * t166;
t305 = qJDD(2) + qJDD(3);
t55 = qJD(6) * t161 + t108 * t311 + t305 * t315;
t483 = t55 / 0.2e1;
t56 = -qJD(6) * t162 + t108 * t315 - t305 * t311;
t482 = t56 / 0.2e1;
t109 = t165 * t310 + t166 * t309;
t103 = qJDD(6) + t109;
t479 = t103 / 0.2e1;
t518 = t257 / 0.2e1;
t465 = t305 / 0.2e1;
t462 = t317 / 0.2e1;
t470 = -t339 / 0.2e1;
t19 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t90 = mrSges(6,1) * t108 - mrSges(6,3) * t305;
t517 = t19 - t90;
t516 = Ifges(5,4) * t339;
t515 = Ifges(4,5) * t248;
t514 = Ifges(4,6) * t247;
t513 = Ifges(6,6) * t339;
t510 = t317 * Ifges(3,2);
t508 = t450 - t526;
t507 = t287 * t495;
t111 = -mrSges(7,1) * t161 + mrSges(7,2) * t162;
t169 = -mrSges(6,1) * t181 - mrSges(6,3) * t307;
t506 = t111 - t169;
t167 = -mrSges(5,2) * t307 + mrSges(5,3) * t181;
t505 = -t167 + t169;
t504 = -t307 * t533 + t339 * t531;
t202 = t312 * t274 - t316 * t275;
t412 = t288 * t318;
t414 = t287 * t318;
t503 = pkin(4) * t412 + qJ(5) * t414;
t394 = qJD(1) * t317;
t395 = qJD(1) * t313;
t447 = pkin(7) * t317;
t448 = pkin(7) * t313;
t500 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t395) * t447 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t394) * t448;
t499 = -mrSges(6,2) * t414 - mrSges(6,3) * t412 - t318 * t524;
t413 = t288 * t314;
t415 = t287 * t314;
t498 = -mrSges(6,2) * t415 - mrSges(6,3) * t413 - t314 * t524;
t116 = -mrSges(7,2) * t174 + mrSges(7,3) * t161;
t117 = mrSges(7,1) * t174 - mrSges(7,3) * t162;
t340 = -t311 * t116 - t315 * t117;
t241 = t257 * pkin(7);
t242 = t258 * pkin(7);
t497 = t241 * t317 + t242 * t313;
t31 = mrSges(7,1) * t103 - mrSges(7,3) * t55;
t32 = -mrSges(7,2) * t103 + mrSges(7,3) * t56;
t496 = t315 * t31 + t311 * t32;
t493 = 0.2e1 * t465;
t491 = -t288 * mrSges(7,3) - t287 * t523 + t521;
t423 = qJDD(1) * pkin(1);
t223 = -pkin(2) * t257 - t423;
t134 = -pkin(3) * t166 + qJDD(4) + t223;
t321 = -qJ(5) * t109 - qJD(5) * t339 + t134;
t12 = t108 * t480 + t321;
t198 = qJDD(2) * pkin(2) - pkin(8) * t258 - t242;
t200 = pkin(8) * t257 + t241;
t114 = -qJD(3) * t187 + t316 * t198 - t200 * t312;
t47 = pkin(3) * t305 - qJ(4) * t165 - qJD(4) * t229 + t114;
t113 = t312 * t198 + t316 * t200 + t237 * t391 + t254 * t392;
t53 = qJ(4) * t166 + qJD(4) * t228 + t113;
t17 = -t309 * t53 + t310 * t47;
t356 = qJDD(5) - t17;
t6 = pkin(5) * t109 - t305 * t480 + t356;
t1 = qJD(6) * t21 + t12 * t315 + t311 * t6;
t2 = -qJD(6) * t22 - t12 * t311 + t315 * t6;
t490 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t272 = -mrSges(3,1) * t317 + mrSges(3,2) * t313;
t489 = -m(3) * pkin(1) - m(4) * t292 - mrSges(2,1) + t272 + t521;
t488 = -m(5) * t86 + m(6) * t83 - t504;
t342 = t21 * t311 - t22 * t315;
t443 = t1 * t311;
t322 = -qJD(6) * t342 + t2 * t315 + t443;
t389 = qJD(6) * t315;
t390 = qJD(6) * t311;
t487 = m(7) * t322 + t116 * t389 - t117 * t390 + t496;
t306 = -qJ(4) + t319;
t486 = -m(3) * pkin(7) + m(4) * t319 - m(7) * (pkin(5) - t306) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) + t531;
t484 = Ifges(7,1) * t483 + Ifges(7,4) * t482 + Ifges(7,5) * t479;
t481 = m(6) + m(7);
t477 = -t161 / 0.2e1;
t476 = -t162 / 0.2e1;
t475 = t162 / 0.2e1;
t474 = -t174 / 0.2e1;
t190 = t247 * t309 + t248 * t310;
t468 = t190 / 0.2e1;
t466 = t229 / 0.2e1;
t463 = t307 / 0.2e1;
t457 = pkin(2) * t312;
t456 = pkin(2) * t313;
t455 = pkin(2) * t316;
t454 = pkin(3) * t229;
t453 = pkin(3) * t302;
t290 = pkin(3) * t303;
t452 = pkin(3) * t309;
t451 = pkin(3) * t310;
t444 = g(3) * t288;
t282 = t288 * pkin(4);
t18 = t309 * t47 + t310 * t53;
t437 = mrSges(4,3) * t228;
t436 = mrSges(4,3) * t229;
t435 = mrSges(7,3) * t315;
t434 = Ifges(3,4) * t313;
t433 = Ifges(3,4) * t317;
t432 = Ifges(4,4) * t229;
t430 = Ifges(7,4) * t311;
t429 = Ifges(7,4) * t315;
t195 = qJD(2) * t247 + t331;
t196 = -qJD(2) * t248 - t332;
t135 = t195 * t309 - t310 * t196;
t422 = t135 * t311;
t421 = t135 * t315;
t420 = t339 * t311;
t189 = -t310 * t247 + t248 * t309;
t419 = t189 * t311;
t418 = t189 * t315;
t277 = t287 * qJ(5);
t410 = t310 * t312;
t408 = t311 * t314;
t407 = t311 * t318;
t406 = t314 * t315;
t404 = t315 * t318;
t259 = -t453 - t456;
t260 = qJ(5) * t413;
t402 = t314 * t259 + t260;
t262 = qJ(5) * t412;
t401 = t318 * t259 + t262;
t396 = t290 + t304;
t393 = qJD(2) * t313;
t385 = Ifges(7,5) * t55 + Ifges(7,6) * t56 + Ifges(7,3) * t103;
t296 = pkin(2) * t393;
t295 = pkin(2) * t395;
t375 = t282 + t277 + t290;
t289 = -pkin(4) - t451;
t374 = qJD(2) * t319;
t367 = -t390 / 0.2e1;
t183 = -pkin(3) * t196 + t296;
t366 = t387 / 0.2e1;
t91 = t109 * mrSges(6,1) + t305 * mrSges(6,2);
t255 = t313 * t374;
t256 = t317 * t374;
t140 = t316 * t255 + t312 * t256 + t274 * t391 + t275 * t392;
t104 = qJ(4) * t196 + qJD(4) * t247 + t140;
t141 = -qJD(3) * t202 - t255 * t312 + t316 * t256;
t105 = -qJ(4) * t195 - qJD(4) * t248 + t141;
t35 = t104 * t309 - t310 * t105;
t252 = pkin(1) + t396;
t362 = -t252 - t277;
t92 = t155 * t309 + t411;
t106 = -t310 * t163 + t164 * t309;
t201 = t316 * t274 + t275 * t312;
t177 = -qJ(4) * t248 + t201;
t178 = qJ(4) * t247 + t202;
t119 = -t310 * t177 + t178 * t309;
t236 = t318 * t252;
t361 = -t306 * t314 + t236;
t291 = pkin(3) + t455;
t220 = t291 * t310 - t309 * t457;
t360 = t304 + t375;
t206 = -pkin(3) * t247 - t292;
t359 = -m(7) * t480 - mrSges(7,3);
t216 = -pkin(4) - t220;
t358 = -pkin(4) * t287 - t453;
t354 = mrSges(3,1) * t313 + mrSges(3,2) * t317;
t348 = Ifges(7,1) * t311 + t429;
t347 = t434 + t510;
t346 = Ifges(7,2) * t315 + t430;
t345 = Ifges(3,5) * t317 - Ifges(3,6) * t313;
t344 = Ifges(7,5) * t311 + Ifges(7,6) * t315;
t343 = t21 * t315 + t22 * t311;
t333 = -qJ(5) * t190 + t206;
t76 = t189 * t480 + t333;
t77 = pkin(5) * t190 + t119;
t34 = t311 * t77 + t315 * t76;
t33 = -t311 * t76 + t315 * t77;
t36 = t104 * t310 + t105 * t309;
t341 = t116 * t315 - t117 * t311;
t93 = t155 * t310 - t147;
t120 = t177 * t309 + t178 * t310;
t338 = t359 * t287;
t221 = pkin(2) * t410 + t291 * t309;
t337 = pkin(1) * t354;
t336 = t189 * t389 + t422;
t335 = t189 * t390 - t421;
t334 = t313 * (Ifges(3,1) * t317 - t434);
t13 = -qJ(5) * t305 - qJD(5) * t307 - t18;
t115 = pkin(4) * t339 - qJ(5) * t181 + t454;
t110 = t115 + t295;
t136 = t195 * t310 + t196 * t309;
t325 = -qJ(5) * t136 - qJD(5) * t190 + t183;
t324 = -m(7) * t453 + t338;
t10 = Ifges(7,4) * t55 + Ifges(7,2) * t56 + Ifges(7,6) * t103;
t123 = t307 * Ifges(6,5) - Ifges(6,3) * t181 - t513;
t125 = Ifges(5,2) * t181 + t307 * Ifges(5,6) + t516;
t15 = -pkin(4) * t305 + t356;
t175 = Ifges(4,2) * t228 + Ifges(4,6) * t307 + t432;
t224 = Ifges(4,4) * t228;
t176 = t229 * Ifges(4,1) + t307 * Ifges(4,5) + t224;
t157 = Ifges(7,4) * t161;
t67 = Ifges(7,1) * t162 + Ifges(7,5) * t174 + t157;
t7 = -pkin(5) * t108 - t13;
t320 = (-Ifges(5,2) * t472 + Ifges(6,3) * t473 - t22 * t435 + t344 * t474 + t346 * t477 + t348 * t476 + t464 * t528 - t520 + t522) * t339 + t528 * t108 + t529 * t109 + (-t22 * t389 - t443 + (t390 + t420) * t21) * mrSges(7,3) + t187 * t436 + t186 * t437 + (t125 + t513) * t469 + (-t420 / 0.2e1 + t367) * t67 + (t123 - t516) * t470 - (t161 * t346 + t162 * t348 + t174 * t344) * qJD(6) / 0.2e1 - (-Ifges(4,2) * t229 + t176 + t224) * t228 / 0.2e1 + t522 * qJD(6) + t7 * t523 + t15 * mrSges(6,2) + t17 * mrSges(5,1) - t18 * mrSges(5,2) - t13 * mrSges(6,3) + (Ifges(6,1) + Ifges(4,3) + Ifges(5,3)) * t305 - t2 * t435 - t229 * (Ifges(4,1) * t228 - t432) / 0.2e1 + (Ifges(5,1) * t470 + Ifges(5,4) * t472 + Ifges(7,5) * t476 + Ifges(7,6) * t477 + Ifges(7,3) * t474 + t464 * t529 - t519 + t532) * t181 - t311 * t10 / 0.2e1 + (Ifges(4,5) * t228 - Ifges(4,6) * t229) * t464 + t175 * t466 + t527 * t472 + (Ifges(7,5) * t315 - Ifges(7,6) * t311) * t479 + (-Ifges(7,2) * t311 + t429) * t482 + (Ifges(7,1) * t315 - t430) * t483 + t315 * t484 - t113 * mrSges(4,2) + t114 * mrSges(4,1) + Ifges(4,5) * t165 + Ifges(4,6) * t166 + t273 * (mrSges(4,1) * t229 + mrSges(4,2) * t228);
t294 = Ifges(3,4) * t394;
t285 = qJ(5) + t452;
t281 = t288 * pkin(9);
t227 = Ifges(3,1) * t395 + Ifges(3,5) * qJD(2) + t294;
t226 = Ifges(3,6) * qJD(2) + qJD(1) * t347;
t215 = qJ(5) + t221;
t214 = -t287 * t408 + t404;
t213 = t287 * t406 + t407;
t212 = t287 * t407 + t406;
t211 = t287 * t404 - t408;
t205 = mrSges(4,1) * t307 - t436;
t204 = -mrSges(4,2) * t307 + t437;
t203 = t295 + t454;
t185 = -mrSges(4,1) * t228 + mrSges(4,2) * t229;
t173 = t339 * pkin(9);
t146 = -mrSges(4,2) * t305 + mrSges(4,3) * t166;
t145 = mrSges(4,1) * t305 - mrSges(4,3) * t165;
t133 = mrSges(6,2) * t181 - mrSges(6,3) * t339;
t132 = -mrSges(5,1) * t181 + mrSges(5,2) * t339;
t118 = pkin(4) * t189 + t333;
t98 = t109 * mrSges(6,3);
t97 = t109 * mrSges(5,2);
t89 = mrSges(5,1) * t305 - mrSges(5,3) * t109;
t88 = -mrSges(5,2) * t305 - mrSges(5,3) * t108;
t78 = -pkin(5) * t189 + t120;
t71 = t115 + t173;
t70 = t110 + t173;
t68 = t106 + t449;
t63 = t93 - t450;
t62 = t92 + t449;
t37 = pkin(4) * t135 + t325;
t30 = t135 * t480 + t325;
t29 = -pkin(5) * t135 + t36;
t28 = pkin(5) * t136 + t35;
t27 = t311 * t68 + t315 * t70;
t26 = -t311 * t70 + t315 * t68;
t25 = t311 * t62 + t315 * t71;
t24 = -t311 * t71 + t315 * t62;
t23 = pkin(4) * t108 + t321;
t4 = -qJD(6) * t34 + t28 * t315 - t30 * t311;
t3 = qJD(6) * t33 + t28 * t311 + t30 * t315;
t5 = [(t134 * mrSges(5,1) + t13 * mrSges(6,1) - t23 * mrSges(6,2) - t18 * mrSges(5,3) + t344 * t479 + t346 * t482 + t348 * t483 - t7 * t350 + t66 * t367 + (-Ifges(6,6) - Ifges(5,4)) * t109 + (Ifges(6,3) + Ifges(5,2)) * t108 + t528 * t493) * t189 + (-Ifges(5,4) * t469 + Ifges(6,6) * t470 + Ifges(6,3) * t472 - Ifges(5,2) * t473 + t123 / 0.2e1 - t125 / 0.2e1 + t528 * t463 + t520) * t135 + (t512 / 0.2e1 + t511 / 0.2e1 + Ifges(5,1) * t469 - Ifges(6,2) * t470 - Ifges(6,6) * t472 + Ifges(5,4) * t473 + Ifges(7,5) * t475 + t529 * t463 + t519 + t527 / 0.2e1) * t136 + t161 * (Ifges(7,4) * t336 - Ifges(7,2) * t335) / 0.2e1 + t33 * t31 + t34 * t32 + (Ifges(3,4) * t258 + Ifges(3,2) * t257) * t462 + (m(5) * t18 - m(6) * t13 + t88 - t90) * t120 + (mrSges(4,1) * t292 + Ifges(4,4) * t248 + Ifges(4,2) * t247) * t166 + t347 * t518 + (-m(5) * t17 + m(6) * t15 - t89 + t91) * t119 + (Ifges(7,1) * t336 - Ifges(7,4) * t335) * t475 + (t15 * mrSges(6,1) + t134 * mrSges(5,2) - t17 * mrSges(5,3) - t23 * mrSges(6,3) + Ifges(5,5) * t465 + Ifges(7,5) * t483 + Ifges(7,6) * t482 + Ifges(7,3) * t479 + (Ifges(6,2) + Ifges(5,1) / 0.2e1) * t109 + (-Ifges(6,6) - Ifges(5,4) / 0.2e1) * t108 - t493 * Ifges(6,4) + t490) * t190 + (-m(7) * (pkin(9) * t412 + t236 + t503) - t212 * mrSges(7,1) - t211 * mrSges(7,2) - mrSges(7,3) * t412 - m(6) * (t361 + t503) - m(5) * t361 + t489 * t318 + t486 * t314) * g(2) + (m(5) * t87 - m(6) * t84 - t505) * t36 + m(4) * (t113 * t202 + t114 * t201 + t140 * t187 + t141 * t186 - t223 * t292 - t273 * t296) + (-mrSges(3,1) * t448 - mrSges(3,2) * t447 + 0.2e1 * Ifges(3,6) * t462) * qJDD(2) + (Ifges(3,1) * t258 + Ifges(3,4) * t518 + Ifges(3,5) * qJDD(2) - t366 * t510) * t313 + t185 * t296 + (-mrSges(4,2) * t292 + Ifges(4,1) * t248 + Ifges(4,4) * t247) * t165 + (t113 * t247 - t114 * t248 - t186 * t195 + t187 * t196) * mrSges(4,3) + m(7) * (t1 * t34 + t2 * t33 + t21 * t4 + t22 * t3 + t29 * t58 + t7 * t78) + (-t214 * mrSges(7,1) + t213 * mrSges(7,2) + ((m(5) + m(6)) * t306 + t486) * t318 + (-m(6) * (t362 - t282) - m(7) * t362 - t288 * t359 + m(5) * t252 - t489) * t314) * g(1) + t258 * t433 / 0.2e1 + t174 * (Ifges(7,5) * t336 - Ifges(7,6) * t335) / 0.2e1 + m(6) * (t118 * t23 + t37 * t95) + m(5) * (t134 * t206 + t183 * t199) + (t514 / 0.2e1 + t515 / 0.2e1 + Ifges(5,5) * t468) * t305 + (t514 + t515) * t465 + t118 * (-t108 * mrSges(6,2) - t98) + (Ifges(5,1) * t109 - Ifges(5,4) * t108 + t385) * t468 + (t1 * t418 - t2 * t419 - t21 * t336 - t22 * t335) * mrSges(7,3) + (t257 * t447 + t258 * t448 + t497) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t497) + (t227 * t462 + t345 * qJD(2) / 0.2e1 - t500) * qJD(2) + Ifges(2,3) * qJDD(1) - t272 * t423 + t66 * t421 / 0.2e1 + t67 * t422 / 0.2e1 - t226 * t393 / 0.2e1 - t337 * t387 - t273 * (-mrSges(4,1) * t196 + mrSges(4,2) * t195) - pkin(1) * (-mrSges(3,1) * t257 + mrSges(3,2) * t258) + t223 * (-mrSges(4,1) * t247 + mrSges(4,2) * t248) + t228 * (Ifges(4,4) * t195 + Ifges(4,2) * t196) / 0.2e1 + t206 * (t108 * mrSges(5,1) + t97) + t201 * t145 + t202 * t146 + t140 * t204 + t141 * t205 + t195 * t176 / 0.2e1 + t196 * t175 / 0.2e1 + t183 * t132 + t78 * t19 + t58 * (mrSges(7,1) * t335 + mrSges(7,2) * t336) + (Ifges(4,5) * t195 + Ifges(4,6) * t196) * t463 + (Ifges(4,1) * t195 + Ifges(4,4) * t196) * t466 + t419 * t484 + t29 * t111 + t3 * t116 + t4 * t117 + t37 * t133 + (t317 * t433 + t334) * t366 + t488 * t35 + (qJD(6) * t67 + t10) * t418 / 0.2e1; t505 * t107 + (t500 + (t337 - t334 / 0.2e1) * qJD(1)) * qJD(1) + (t379 - t194) * t204 + (-t380 - t193) * t205 + t504 * t106 + ((t113 * t312 + t114 * t316 + (-t186 * t312 + t187 * t316) * qJD(3)) * pkin(2) - t186 * t193 - t187 * t194 + t273 * t295) * m(4) + t495 * (m(4) * t456 - m(5) * t259 + t354 + t525) + (-t13 * t215 + t15 * t216 - t106 * t83 - t110 * t95 - g(1) * (-pkin(4) * t414 + t401) - g(2) * (-pkin(4) * t415 + t402) + t526 * t84) * m(6) + (t106 * t86 + t17 * t220 + t18 * t221 - t199 * t203 + (t219 - t107) * t87) * m(5) + t517 * t215 + (m(7) * t343 - t340 + t488) * (t309 * t316 + t410) * qJD(3) * pkin(2) + (-g(1) * t401 - g(2) * t402 - t21 * t26 + t215 * t7 - t22 * t27 + t508 * t58) * m(7) + t508 * t111 + t320 + (-m(4) * t304 - m(5) * t396 - m(7) * (t281 + t360) - m(6) * t360 + t272 + t491) * g(3) + (-t314 * t338 + t498) * g(2) + (-t318 * t338 + t499) * g(1) + t226 * t395 / 0.2e1 - t345 * t387 / 0.2e1 - t185 * t295 + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t257 + Ifges(3,5) * t258 - t241 * mrSges(3,2) - t242 * mrSges(3,1) + t219 * t167 + t220 * t89 + t221 * t88 - t208 * t169 + t216 * t91 - t203 * t132 + t145 * t455 + t146 * t457 - t27 * t116 - t26 * t117 - t110 * t133 + t487 * (-pkin(9) + t216) - (-Ifges(3,2) * t395 + t227 + t294) * t394 / 0.2e1; -t63 * t111 - t115 * t133 - t25 * t116 - t24 * t117 - t132 * t454 - t186 * t204 + t187 * t205 + t289 * t91 + t89 * t451 + t88 * t452 + t320 + t505 * t93 + t504 * t92 + t517 * t285 + t506 * qJD(5) + (-t314 * t324 + t498) * g(2) + (-t318 * t324 + t499) * g(1) + (-g(1) * t262 - g(2) * t260 - t21 * t24 - t22 * t25 + t285 * t7 + (qJD(5) - t63) * t58) * m(7) + (-t13 * t285 + t15 * t289 - t115 * t95 - t83 * t92 - g(1) * (t318 * t358 + t262) - g(2) * (t314 * t358 + t260) + (-qJD(5) + t93) * t84) * m(6) + ((t17 * t310 + t18 * t309) * pkin(3) - t199 * t454 + t86 * t92 - t87 * t93) * m(5) + t487 * (-pkin(9) + t289) + (-m(5) * t290 - m(7) * (t281 + t375) - m(6) * t375 + t491) * g(3) + (m(5) * t453 + t525) * t495; -t311 * t31 + t315 * t32 + t97 - t98 - t533 * t108 + t340 * qJD(6) + (-t167 - t506) * t181 + (t340 + t504) * t339 + (-g(1) * t314 + g(2) * t318) * (m(5) + t481) + (t1 * t315 - t174 * t343 - t181 * t58 - t2 * t311) * m(7) + (t181 * t84 - t339 * t83 + t23) * m(6) + (-t181 * t87 + t339 * t86 + t134) * m(5); -t506 * t307 + t341 * qJD(6) + t481 * t444 + (t133 + t341) * t339 + t91 + (-t307 * t58 - t339 * t342 + t322 - t507) * m(7) + (t307 * t84 + t339 * t95 + t15 - t507) * m(6) + t496; -t58 * (mrSges(7,1) * t162 + mrSges(7,2) * t161) + (Ifges(7,1) * t161 - t431) * t476 + t66 * t475 + (Ifges(7,5) * t161 - Ifges(7,6) * t162) * t474 - t21 * t116 + t22 * t117 - g(1) * (mrSges(7,1) * t211 - mrSges(7,2) * t212) - g(2) * (mrSges(7,1) * t213 + mrSges(7,2) * t214) + t350 * t444 + (t161 * t21 + t162 * t22) * mrSges(7,3) + t385 + (-Ifges(7,2) * t162 + t157 + t67) * t477 + t490;];
tau  = t5;
