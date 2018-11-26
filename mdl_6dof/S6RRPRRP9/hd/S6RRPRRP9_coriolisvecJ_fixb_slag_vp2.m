% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:16:43
% EndTime: 2018-11-23 17:17:09
% DurationCPUTime: 26.34s
% Computational Cost: add. (17607->792), mult. (46911->1110), div. (0->0), fcn. (37752->10), ass. (0->344)
t303 = sin(pkin(11));
t305 = cos(pkin(11));
t308 = sin(qJ(4));
t311 = cos(qJ(4));
t276 = t303 * t311 + t305 * t308;
t304 = sin(pkin(6));
t312 = cos(qJ(2));
t385 = t304 * t312;
t321 = t276 * t385;
t231 = qJD(1) * t321;
t271 = t276 * qJD(4);
t513 = -t231 + t271;
t487 = Ifges(7,4) + Ifges(6,4);
t309 = sin(qJ(2));
t334 = pkin(2) * t309 - qJ(3) * t312;
t380 = qJD(1) * t304;
t261 = t334 * t380;
t366 = t309 * t380;
t306 = cos(pkin(6));
t379 = qJD(1) * t306;
t372 = pkin(1) * t379;
t262 = -pkin(8) * t366 + t312 * t372;
t196 = t305 * t261 - t303 * t262;
t384 = t305 * t312;
t322 = (pkin(3) * t309 - pkin(9) * t384) * t304;
t169 = qJD(1) * t322 + t196;
t197 = t303 * t261 + t305 * t262;
t365 = t312 * t380;
t358 = t303 * t365;
t179 = -pkin(9) * t358 + t197;
t275 = t303 * t308 - t311 * t305;
t424 = pkin(9) + qJ(3);
t282 = t424 * t303;
t283 = t424 * t305;
t477 = -t311 * t282 - t283 * t308;
t494 = -t275 * qJD(3) + qJD(4) * t477 - t308 * t169 - t311 * t179;
t485 = Ifges(7,2) + Ifges(6,2);
t484 = Ifges(7,6) + Ifges(6,6);
t512 = -pkin(10) * t366 + t494;
t263 = pkin(8) * t365 + t309 * t372;
t224 = pkin(3) * t358 + t263;
t320 = t275 * t385;
t232 = qJD(1) * t320;
t270 = t275 * qJD(4);
t511 = -t224 + (-t232 + t270) * pkin(10) + t513 * pkin(4);
t287 = qJD(4) - t365;
t307 = sin(qJ(5));
t310 = cos(qJ(5));
t295 = qJD(2) + t379;
t248 = t295 * t305 - t303 * t366;
t249 = t295 * t303 + t305 * t366;
t329 = t248 * t308 + t311 * t249;
t160 = t287 * t310 - t307 * t329;
t359 = t311 * t248 - t249 * t308;
t183 = qJD(5) - t359;
t161 = t287 * t307 + t310 * t329;
t508 = t487 * t161;
t479 = t160 * t485 + t183 * t484 + t508;
t510 = -t479 / 0.2e1;
t488 = Ifges(7,1) + Ifges(6,1);
t486 = Ifges(7,5) + Ifges(6,5);
t230 = -t282 * t308 + t283 * t311;
t493 = -qJD(3) * t276 - qJD(4) * t230 - t169 * t311 + t308 * t179;
t198 = t232 * t307 + t310 * t366;
t374 = qJD(5) * t310;
t509 = -t276 * t374 - t198;
t495 = pkin(4) * t366 - t493;
t507 = -t307 * t512 + t310 * t511;
t301 = -pkin(3) * t305 - pkin(2);
t211 = pkin(4) * t275 - pkin(10) * t276 + t301;
t506 = t211 * t374 + t307 * t511 + t310 * t512;
t505 = t487 * t160;
t504 = t487 * t310;
t503 = t487 * t307;
t182 = Ifges(5,4) * t359;
t226 = -t295 * pkin(2) + qJD(3) - t262;
t184 = -t248 * pkin(3) + t226;
t398 = t287 * Ifges(5,5);
t235 = qJ(3) * t295 + t263;
t256 = (-pkin(2) * t312 - qJ(3) * t309 - pkin(1)) * t304;
t241 = qJD(1) * t256;
t170 = -t303 * t235 + t305 * t241;
t123 = -pkin(3) * t365 - t249 * pkin(9) + t170;
t171 = t305 * t235 + t303 * t241;
t141 = pkin(9) * t248 + t171;
t70 = t123 * t308 + t141 * t311;
t65 = pkin(10) * t287 + t70;
t87 = -pkin(4) * t359 - pkin(10) * t329 + t184;
t27 = -t307 * t65 + t310 * t87;
t18 = -qJ(6) * t161 + t27;
t10 = pkin(5) * t183 + t18;
t28 = t307 * t87 + t310 * t65;
t19 = qJ(6) * t160 + t28;
t331 = t27 * t310 + t28 * t307;
t349 = mrSges(7,1) * t307 + mrSges(7,2) * t310;
t351 = mrSges(6,1) * t307 + mrSges(6,2) * t310;
t434 = t310 / 0.2e1;
t439 = -t307 / 0.2e1;
t69 = t123 * t311 - t308 * t141;
t64 = -pkin(4) * t287 - t69;
t44 = -pkin(5) * t160 + qJD(6) + t64;
t451 = t183 / 0.2e1;
t455 = t161 / 0.2e1;
t457 = t160 / 0.2e1;
t474 = t310 * t488 - t503;
t475 = -t307 * t485 + t504;
t476 = -t307 * t484 + t310 * t486;
t478 = t161 * t488 + t183 * t486 + t505;
t471 = t331 * mrSges(6,3) + (t10 * t310 + t19 * t307) * mrSges(7,3) - t349 * t44 - t351 * t64 - t434 * t478 - t439 * t479 - t451 * t476 - t455 * t474 - t457 * t475;
t502 = t184 * mrSges(5,2) + t182 / 0.2e1 + t398 / 0.2e1 - t69 * mrSges(5,3) - t471;
t483 = Ifges(7,3) + Ifges(6,3);
t199 = -t232 * t310 + t307 * t366;
t218 = t310 * t230;
t326 = qJ(6) * t270 - qJD(6) * t276;
t501 = qJ(6) * t199 + t326 * t310 + (-t218 + (qJ(6) * t276 - t211) * t307) * qJD(5) + t507 + t513 * pkin(5);
t500 = (-qJD(5) * t230 + t326) * t307 + t506 + t509 * qJ(6);
t58 = t161 * Ifges(7,5) + t160 * Ifges(7,6) + t183 * Ifges(7,3);
t59 = t161 * Ifges(6,5) + t160 * Ifges(6,6) + t183 * Ifges(6,3);
t499 = t59 + t58;
t375 = qJD(5) * t307;
t498 = -t230 * t375 + t506;
t152 = t307 * t211 + t218;
t497 = -qJD(5) * t152 + t507;
t496 = t495 + (-t270 * t307 - t509) * pkin(5);
t492 = t307 * t486 + t310 * t484;
t491 = t310 * t485 + t503;
t490 = t307 * t488 + t504;
t489 = qJD(5) * t276;
t317 = qJD(2) * t321;
t140 = qJD(1) * t317 + qJD(4) * t329;
t316 = qJD(2) * t320;
t139 = -qJD(1) * t316 + qJD(4) * t359;
t378 = qJD(2) * t304;
t362 = qJD(1) * t378;
t356 = t309 * t362;
t81 = qJD(5) * t160 + t139 * t310 + t307 * t356;
t82 = -qJD(5) * t161 - t139 * t307 + t310 * t356;
t12 = Ifges(7,5) * t81 + Ifges(7,6) * t82 + Ifges(7,3) * t140;
t13 = Ifges(6,5) * t81 + Ifges(6,6) * t82 + Ifges(6,3) * t140;
t482 = t13 + t12;
t481 = t140 * t484 + t485 * t82 + t487 * t81;
t480 = t486 * t140 + t487 * t82 + t488 * t81;
t386 = t304 * t309;
t268 = -t303 * t386 + t305 * t306;
t269 = t303 * t306 + t305 * t386;
t201 = t268 * t308 + t269 * t311;
t296 = pkin(8) * t386;
t432 = pkin(1) * t312;
t258 = t296 + (-pkin(2) - t432) * t306;
t207 = -t268 * pkin(3) + t258;
t328 = t311 * t268 - t269 * t308;
t109 = -pkin(4) * t328 - t201 * pkin(10) + t207;
t273 = t306 * t309 * pkin(1) + pkin(8) * t385;
t255 = qJ(3) * t306 + t273;
t191 = -t303 * t255 + t305 * t256;
t148 = -pkin(3) * t385 - t269 * pkin(9) + t191;
t192 = t305 * t255 + t303 * t256;
t166 = pkin(9) * t268 + t192;
t92 = t308 * t148 + t311 * t166;
t90 = -pkin(10) * t385 + t92;
t42 = t307 * t109 + t310 * t90;
t238 = (qJD(2) * t334 - qJD(3) * t309) * t304;
t220 = qJD(1) * t238;
t373 = t306 * t432;
t294 = qJD(2) * t373;
t253 = -pkin(8) * t356 + qJD(1) * t294;
t221 = qJD(3) * t295 + t253;
t163 = t305 * t220 - t303 * t221;
t319 = qJD(2) * t322;
t124 = qJD(1) * t319 + t163;
t164 = t303 * t220 + t305 * t221;
t355 = t312 * t362;
t327 = t303 * t355;
t142 = -pkin(9) * t327 + t164;
t376 = qJD(4) * t311;
t377 = qJD(4) * t308;
t25 = t123 * t376 + t308 * t124 - t141 * t377 + t311 * t142;
t23 = pkin(10) * t356 + t25;
t265 = t273 * qJD(2);
t254 = qJD(1) * t265;
t214 = pkin(3) * t327 + t254;
t57 = t140 * pkin(4) - t139 * pkin(10) + t214;
t5 = t310 * t23 + t307 * t57 + t87 * t374 - t375 * t65;
t6 = -qJD(5) * t28 - t23 * t307 + t310 * t57;
t473 = -t307 * t6 + t310 * t5;
t396 = t287 * Ifges(5,3);
t399 = t329 * Ifges(5,5);
t401 = t359 * Ifges(5,6);
t110 = t396 + t399 + t401;
t369 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t370 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t371 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t397 = t287 * Ifges(5,6);
t418 = Ifges(5,4) * t329;
t402 = t359 * Ifges(5,2);
t111 = t397 + t402 + t418;
t463 = t111 / 0.2e1;
t470 = t370 * t160 + t371 * t161 + t369 * t183 - t19 * mrSges(7,2) - t28 * mrSges(6,2) - t70 * mrSges(5,3) - t463 + t58 / 0.2e1 + t59 / 0.2e1 - t418 / 0.2e1 + t10 * mrSges(7,1) + t184 * mrSges(5,1) + t27 * mrSges(6,1) - t397 / 0.2e1;
t469 = -0.2e1 * pkin(1);
t468 = Ifges(5,2) / 0.2e1;
t467 = t81 / 0.2e1;
t466 = t82 / 0.2e1;
t400 = t329 * Ifges(5,1);
t112 = t182 + t398 + t400;
t462 = -t112 / 0.2e1;
t461 = t112 / 0.2e1;
t460 = t140 / 0.2e1;
t458 = -t160 / 0.2e1;
t456 = -t161 / 0.2e1;
t452 = -t183 / 0.2e1;
t448 = t328 / 0.2e1;
t446 = t201 / 0.2e1;
t444 = t268 / 0.2e1;
t443 = t269 / 0.2e1;
t442 = -t303 / 0.2e1;
t441 = t305 / 0.2e1;
t440 = t306 / 0.2e1;
t436 = t309 / 0.2e1;
t431 = t25 * mrSges(5,2);
t26 = -t123 * t377 + t124 * t311 - t141 * t376 - t308 * t142;
t430 = t26 * mrSges(5,1);
t427 = t69 * mrSges(5,1);
t425 = t70 * mrSges(5,2);
t423 = -qJ(6) - pkin(10);
t422 = mrSges(4,2) * t305;
t421 = Ifges(3,4) * t309;
t420 = Ifges(4,4) * t303;
t419 = Ifges(4,4) * t305;
t413 = Ifges(4,5) * t249;
t412 = Ifges(4,5) * t309;
t411 = Ifges(3,6) * t295;
t410 = Ifges(4,6) * t248;
t409 = Ifges(4,6) * t309;
t408 = t139 * Ifges(5,1);
t407 = t139 * Ifges(5,4);
t406 = t140 * Ifges(5,4);
t395 = t295 * Ifges(3,5);
t394 = t312 * Ifges(3,2);
t168 = mrSges(5,1) * t287 - mrSges(5,3) * t329;
t94 = -mrSges(6,1) * t160 + mrSges(6,2) * t161;
t393 = t168 - t94;
t116 = pkin(4) * t329 - pkin(10) * t359;
t40 = t307 * t116 + t310 * t69;
t392 = qJ(6) * t310;
t391 = t359 * t307;
t389 = t276 * t307;
t387 = t303 * t312;
t381 = -mrSges(3,1) * t295 - mrSges(4,1) * t248 + mrSges(4,2) * t249 + mrSges(3,3) * t366;
t363 = t309 * t378;
t264 = -pkin(8) * t363 + t294;
t247 = qJD(3) * t306 + t264;
t178 = t303 * t238 + t305 * t247;
t236 = mrSges(4,1) * t327 + t355 * t422;
t368 = Ifges(5,5) * t139 - Ifges(5,6) * t140 + Ifges(5,3) * t356;
t30 = -t82 * mrSges(7,1) + t81 * mrSges(7,2);
t361 = qJD(5) * t423;
t80 = t140 * mrSges(5,1) + t139 * mrSges(5,2);
t41 = t310 * t109 - t307 * t90;
t39 = t310 * t116 - t307 * t69;
t91 = t148 * t311 - t308 * t166;
t151 = t310 * t211 - t230 * t307;
t177 = t305 * t238 - t303 * t247;
t357 = t378 * t387;
t1 = pkin(5) * t140 - qJ(6) * t81 - qJD(6) * t161 + t6;
t2 = qJ(6) * t82 + qJD(6) * t160 + t5;
t354 = -t1 * t310 - t2 * t307;
t353 = -t5 * t307 - t6 * t310;
t89 = pkin(4) * t385 - t91;
t352 = mrSges(6,1) * t310 - mrSges(6,2) * t307;
t350 = mrSges(7,1) * t310 - mrSges(7,2) * t307;
t348 = Ifges(4,1) * t305 - t420;
t343 = -Ifges(4,2) * t303 + t419;
t332 = t10 * t307 - t19 * t310;
t330 = t27 * t307 - t28 * t310;
t146 = t177 + t319;
t162 = -pkin(9) * t357 + t178;
t38 = t146 * t311 - t148 * t377 - t308 * t162 - t166 * t376;
t325 = mrSges(4,1) * t309 - mrSges(4,3) * t384;
t324 = -mrSges(4,2) * t309 - mrSges(4,3) * t387;
t180 = -t307 * t201 - t310 * t385;
t323 = -t310 * t201 + t307 * t385;
t37 = t308 * t146 + t148 * t376 + t311 * t162 - t166 * t377;
t34 = pkin(10) * t363 + t37;
t156 = qJD(4) * t328 - t316;
t157 = qJD(4) * t201 + t317;
t225 = pkin(3) * t357 + t265;
t74 = t157 * pkin(4) - t156 * pkin(10) + t225;
t7 = t109 * t374 + t307 * t74 + t310 * t34 - t375 * t90;
t318 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t35 = -pkin(4) * t363 - t38;
t8 = -qJD(5) * t42 - t307 * t34 + t310 * t74;
t24 = -pkin(4) * t356 - t26;
t302 = -pkin(5) * t310 - pkin(4);
t290 = Ifges(3,4) * t365;
t289 = t423 * t310;
t288 = t423 * t307;
t285 = Ifges(3,5) * t355;
t272 = -t296 + t373;
t267 = -qJD(6) * t307 + t310 * t361;
t266 = qJD(6) * t310 + t307 * t361;
t260 = -t295 * mrSges(3,2) + mrSges(3,3) * t365;
t243 = t325 * t362;
t242 = t324 * t362;
t234 = Ifges(3,1) * t366 + t290 + t395;
t233 = t411 + (t394 + t421) * t380;
t217 = -mrSges(4,1) * t365 - t249 * mrSges(4,3);
t216 = mrSges(4,2) * t365 + t248 * mrSges(4,3);
t209 = (t312 * t348 + t412) * t362;
t208 = (t312 * t343 + t409) * t362;
t195 = pkin(5) * t389 - t477;
t176 = Ifges(4,1) * t249 + Ifges(4,4) * t248 - Ifges(4,5) * t365;
t175 = Ifges(4,4) * t249 + Ifges(4,2) * t248 - Ifges(4,6) * t365;
t174 = -Ifges(4,3) * t365 + t410 + t413;
t167 = -mrSges(5,2) * t287 + mrSges(5,3) * t359;
t127 = -qJ(6) * t389 + t152;
t119 = -mrSges(5,2) * t356 - mrSges(5,3) * t140;
t118 = mrSges(5,1) * t356 - mrSges(5,3) * t139;
t117 = pkin(5) * t275 - t276 * t392 + t151;
t115 = -mrSges(5,1) * t359 + mrSges(5,2) * t329;
t106 = mrSges(6,1) * t183 - mrSges(6,3) * t161;
t105 = mrSges(7,1) * t183 - mrSges(7,3) * t161;
t104 = -mrSges(6,2) * t183 + mrSges(6,3) * t160;
t103 = -mrSges(7,2) * t183 + mrSges(7,3) * t160;
t96 = qJD(5) * t323 - t307 * t156 + t310 * t363;
t95 = qJD(5) * t180 + t310 * t156 + t307 * t363;
t93 = -mrSges(7,1) * t160 + mrSges(7,2) * t161;
t68 = Ifges(5,5) * t356 - t406 + t408;
t67 = -t140 * Ifges(5,2) + Ifges(5,6) * t356 + t407;
t54 = -pkin(5) * t180 + t89;
t50 = pkin(5) * t391 + t70;
t48 = -mrSges(6,2) * t140 + mrSges(6,3) * t82;
t47 = -mrSges(7,2) * t140 + mrSges(7,3) * t82;
t46 = mrSges(6,1) * t140 - mrSges(6,3) * t81;
t45 = mrSges(7,1) * t140 - mrSges(7,3) * t81;
t32 = qJ(6) * t180 + t42;
t31 = -mrSges(6,1) * t82 + mrSges(6,2) * t81;
t29 = -qJ(6) * t391 + t40;
t21 = -pkin(5) * t328 + qJ(6) * t323 + t41;
t20 = pkin(5) * t329 - t359 * t392 + t39;
t11 = -pkin(5) * t96 + t35;
t9 = -pkin(5) * t82 + t24;
t4 = qJ(6) * t96 + qJD(6) * t180 + t7;
t3 = pkin(5) * t157 - qJ(6) * t95 + qJD(6) * t323 + t8;
t14 = [t478 * t95 / 0.2e1 + t479 * t96 / 0.2e1 + t481 * t180 / 0.2e1 + (t157 * t483 + t484 * t96 + t486 * t95) * t451 + (t157 * t484 + t485 * t96 + t487 * t95) * t457 + (t486 * t157 + t487 * t96 + t488 * t95) * t455 + ((-t273 * mrSges(3,3) - Ifges(3,6) * t306 + Ifges(5,5) * t446 + Ifges(5,6) * t448 + Ifges(4,5) * t443 + Ifges(4,6) * t444 + (mrSges(3,1) * t469 - 0.3e1 / 0.2e1 * t421) * t304) * t309 + (Ifges(3,5) * t440 + (Ifges(4,1) * t269 + Ifges(4,4) * t268) * t441 + (Ifges(4,4) * t269 + Ifges(4,2) * t268) * t442 - t272 * mrSges(3,3) + (mrSges(3,2) * t469 + (0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * t305 * Ifges(4,5) + 0.3e1 / 0.2e1 * t303 * Ifges(4,6)) * t312) * t304 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,3)) * t386) * t312) * t362 + t258 * t236 + t192 * t242 + t191 * t243 + ((-t262 * mrSges(3,3) + t176 * t441 + t175 * t442 + t226 * (mrSges(4,1) * t303 + t422) + t395 / 0.2e1 + t248 * t343 / 0.2e1 + t249 * t348 / 0.2e1 + t234 / 0.2e1 + (-t170 * t305 - t171 * t303) * mrSges(4,3)) * t312 + (-t263 * mrSges(3,3) - t233 / 0.2e1 + t174 / 0.2e1 + t110 / 0.2e1 + t399 / 0.2e1 + t401 / 0.2e1 - t425 + t427 + t396 / 0.2e1 - t411 / 0.2e1 - t171 * mrSges(4,2) + t170 * mrSges(4,1) + t410 / 0.2e1 + t413 / 0.2e1) * t309) * t378 + t381 * t265 + t163 * (-mrSges(4,1) * t385 - t269 * mrSges(4,3)) + t253 * (-t306 * mrSges(3,2) + mrSges(3,3) * t385) + t164 * (mrSges(4,2) * t385 + t268 * mrSges(4,3)) + (-mrSges(3,1) * t306 - mrSges(4,1) * t268 + mrSges(4,2) * t269 + mrSges(3,3) * t386) * t254 - t368 * t385 / 0.2e1 - t480 * t323 / 0.2e1 + t24 * (-mrSges(6,1) * t180 - mrSges(6,2) * t323) + t9 * (-mrSges(7,1) * t180 - mrSges(7,2) * t323) - t482 * t328 / 0.2e1 + (t180 * t484 - t323 * t486 - t328 * t483) * t460 + (t180 * t485 - t323 * t487 - t328 * t484) * t466 + (t487 * t180 - t323 * t488 - t328 * t486) * t467 + t6 * (-mrSges(6,1) * t328 + mrSges(6,3) * t323) + t1 * (-mrSges(7,1) * t328 + mrSges(7,3) * t323) + (-t156 * t69 - t157 * t70 - t201 * t26 + t25 * t328) * mrSges(5,3) - t140 * (Ifges(5,4) * t201 + Ifges(5,2) * t328 - Ifges(5,6) * t385) / 0.2e1 + t139 * (Ifges(5,1) * t201 + Ifges(5,4) * t328 - Ifges(5,5) * t385) / 0.2e1 + t214 * (-mrSges(5,1) * t328 + mrSges(5,2) * t201) + t5 * (mrSges(6,2) * t328 + mrSges(6,3) * t180) + t2 * (mrSges(7,2) * t328 + mrSges(7,3) * t180) - t385 * t430 + t287 * (Ifges(5,5) * t156 - Ifges(5,6) * t157) / 0.2e1 + t178 * t216 + t177 * t217 + t225 * t115 + m(7) * (t1 * t21 + t10 * t3 + t11 * t44 + t19 * t4 + t2 * t32 + t54 * t9) + m(6) * (t24 * t89 + t27 * t8 + t28 * t7 + t35 * t64 + t41 * t6 + t42 * t5) + m(5) * (t184 * t225 + t207 * t214 + t25 * t92 + t26 * t91 + t37 * t70 + t38 * t69) + m(4) * (t163 * t191 + t164 * t192 + t170 * t177 + t171 * t178 + t226 * t265 + t254 * t258) + m(3) * (t253 * t273 - t254 * t272 - t262 * t265 + t263 * t264) + t207 * t80 + t184 * (mrSges(5,1) * t157 + mrSges(5,2) * t156) + t37 * t167 + t38 * t168 + t19 * (-mrSges(7,2) * t157 + mrSges(7,3) * t96) + t28 * (-mrSges(6,2) * t157 + mrSges(6,3) * t96) + t27 * (mrSges(6,1) * t157 - mrSges(6,3) * t95) + t10 * (mrSges(7,1) * t157 - mrSges(7,3) * t95) - t157 * t111 / 0.2e1 + t92 * t119 + t91 * t118 + t8 * t106 + t4 * t103 + t7 * t104 + t3 * t105 + t64 * (-mrSges(6,1) * t96 + mrSges(6,2) * t95) + t44 * (-mrSges(7,1) * t96 + mrSges(7,2) * t95) + t11 * t93 + t35 * t94 + t89 * t31 + t54 * t30 + t21 * t45 + t41 * t46 + t32 * t47 + t42 * t48 + t156 * t461 + t209 * t443 + t208 * t444 + t68 * t446 + t67 * t448 + t285 * t440 + t385 * t431 + t264 * t260 + t499 * t157 / 0.2e1 + t359 * (Ifges(5,4) * t156 - Ifges(5,2) * t157) / 0.2e1 + t329 * (Ifges(5,1) * t156 - Ifges(5,4) * t157) / 0.2e1; ((t312 * (Ifges(4,5) * t384 - Ifges(4,6) * t387 + Ifges(4,3) * t309) / 0.2e1 + t394 * t436 + pkin(1) * (mrSges(3,1) * t309 + mrSges(3,2) * t312)) * t380 + t233 * t436 - t226 * (mrSges(4,1) * t387 + mrSges(4,2) * t384) - t295 * (Ifges(3,5) * t312 - Ifges(3,6) * t309) / 0.2e1 - t171 * t324 - t170 * t325 - t248 * (Ifges(4,4) * t384 - Ifges(4,2) * t387 + t409) / 0.2e1 - t249 * (Ifges(4,1) * t384 - Ifges(4,4) * t387 + t412) / 0.2e1 - t176 * t384 / 0.2e1 + t175 * t387 / 0.2e1 + t309 * t425 - t309 * t427 + (t262 * t312 + t263 * t309) * mrSges(3,3) + ((Ifges(4,5) * t303 / 0.2e1 + Ifges(4,6) * t441 + Ifges(5,5) * t276 / 0.2e1 - Ifges(5,6) * t275 / 0.2e1 - Ifges(3,6)) * t309 + ((Ifges(4,2) * t305 + t420) * t442 + (Ifges(4,1) * t303 + t419) * t441) * t312) * qJD(2) - (t234 + t290) * t312 / 0.2e1 - ((Ifges(3,1) * t312 - t421) * t380 + t174 + 0.2e1 * t110) * t309 / 0.2e1) * t380 - t329 * (-Ifges(5,1) * t232 - Ifges(5,4) * t231) / 0.2e1 - (t461 + t400 / 0.2e1 + t502) * t270 + (t214 * mrSges(5,1) - t67 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1 - t25 * mrSges(5,3) - t407 / 0.2e1 + t370 * t82 + t371 * t81 + (t468 + t369) * t140 + t318) * t275 - t262 * t260 + (-t163 * mrSges(4,3) - qJD(3) * t217 - qJ(3) * t243 + t209 / 0.2e1 + t254 * mrSges(4,2)) * t303 + (-pkin(2) * t254 + (-t170 * t303 + t171 * t305) * qJD(3) + (-t163 * t303 + t164 * t305) * qJ(3) - t170 * t196 - t171 * t197 - t226 * t263) * m(4) + (-t402 / 0.2e1 + t470) * t271 - t254 * mrSges(3,1) - t253 * mrSges(3,2) - pkin(2) * t236 + (t164 * mrSges(4,3) + qJD(3) * t216 + qJ(3) * t242 + t208 / 0.2e1 - t254 * mrSges(4,1)) * t305 - t381 * t263 + (t231 * t70 - t232 * t69) * mrSges(5,3) - t287 * (-Ifges(5,5) * t232 - Ifges(5,6) * t231) / 0.2e1 - t184 * (mrSges(5,1) * t231 - mrSges(5,2) * t232) - t359 * (-Ifges(5,4) * t232 - Ifges(5,2) * t231) / 0.2e1 + (t68 / 0.2e1 + t214 * mrSges(5,2) - t406 / 0.2e1 + t408 / 0.2e1 - t26 * mrSges(5,3) + t24 * t351 + t9 * t349 + t354 * mrSges(7,3) + t353 * mrSges(6,3) + (mrSges(6,3) * t330 + mrSges(7,3) * t332 + t310 * t510 + t350 * t44 + t352 * t64) * qJD(5) + t474 * t467 + t475 * t466 + t476 * t460 + t480 * t434 + (qJD(5) * t478 + t481) * t439) * t276 + t230 * t119 - t27 * (mrSges(6,1) * t231 - mrSges(6,3) * t199) - t10 * (mrSges(7,1) * t231 - mrSges(7,3) * t199) - t28 * (-mrSges(6,2) * t231 + mrSges(6,3) * t198) - t19 * (-mrSges(7,2) * t231 + mrSges(7,3) * t198) - t197 * t216 - t196 * t217 - t224 * t115 - t64 * (-mrSges(6,1) * t198 + mrSges(6,2) * t199) - t44 * (-mrSges(7,1) * t198 + mrSges(7,2) * t199) - (t31 - t118) * t477 + (-t184 * t224 + t214 * t301 + t230 * t25 + t26 * t477 + t493 * t69 + t494 * t70) * m(5) + (t151 * t6 + t152 * t5 - t24 * t477 + t27 * t497 + t28 * t498 + t495 * t64) * m(6) + t195 * t30 + t151 * t46 + t152 * t48 + t127 * t47 + t117 * t45 + t198 * t510 + t285 - t232 * t462 + t231 * t463 + t301 * t80 + t493 * t168 + t494 * t167 + t495 * t94 + t496 * t93 + t497 * t106 + t498 * t104 - t478 * t199 / 0.2e1 - t499 * t231 / 0.2e1 + t500 * t103 + t501 * t105 + (t1 * t117 + t10 * t501 + t127 * t2 + t19 * t500 + t195 * t9 + t44 * t496) * m(7) + (t198 * t484 + t199 * t486 + t231 * t483 + t489 * t492) * t452 + (t198 * t485 + t199 * t487 + t231 * t484 + t489 * t491) * t458 + (t487 * t198 + t199 * t488 + t486 * t231 + t490 * t489) * t456; -t359 * t167 - t248 * t216 + t249 * t217 + (-t93 + t393) * t329 + (t45 + t46 + t183 * (t103 + t104)) * t310 + (t47 + t48 - t183 * (t105 + t106)) * t307 + t80 + t236 + (-t183 * t332 - t329 * t44 - t354) * m(7) + (-t183 * t330 - t329 * t64 - t353) * m(6) + (t329 * t69 - t359 * t70 + t214) * m(5) + (t170 * t249 - t171 * t248 + t254) * m(4); t480 * t307 / 0.2e1 + t481 * t434 + (t462 + (t468 - Ifges(5,1) / 0.2e1) * t329 - t502) * t359 + (-t29 + t266) * t103 + t368 - t431 + ((-m(6) * t331 - t307 * t104 - t310 * t106) * qJD(5) + m(6) * t473 - t307 * t46 + t310 * t48) * pkin(10) + t473 * mrSges(6,3) - t9 * t350 - t24 * t352 + t430 - t470 * t329 + (-t1 * t307 + t2 * t310) * mrSges(7,3) + (-t20 + t267) * t105 + t393 * t70 + (-t471 + (m(7) * t44 + t93) * t307 * pkin(5)) * qJD(5) - t69 * t167 - t39 * t106 - t40 * t104 - t50 * t93 - pkin(4) * t31 + (-pkin(4) * t24 - t27 * t39 - t28 * t40 - t64 * t70) * m(6) + t288 * t45 - t289 * t47 + t302 * t30 + t490 * t467 + t491 * t466 + t492 * t460 - m(7) * (t10 * t20 + t19 * t29 + t44 * t50) + m(7) * (t1 * t288 + t10 * t267 + t19 * t266 - t2 * t289 + t302 * t9); (-(-t10 + t18) * t19 + (-t161 * t44 + t1) * pkin(5)) * m(7) + (-t161 * t93 + t45) * pkin(5) + t482 - t64 * (mrSges(6,1) * t161 + mrSges(6,2) * t160) - t44 * (mrSges(7,1) * t161 + mrSges(7,2) * t160) + t28 * t106 - t18 * t103 - t27 * t104 + t19 * t105 + (t10 * t160 + t161 * t19) * mrSges(7,3) + (t160 * t27 + t161 * t28) * mrSges(6,3) + t318 + (t160 * t488 - t508) * t456 + t479 * t455 + (t160 * t486 - t161 * t484) * t452 + (-t161 * t485 + t478 + t505) * t458; -t160 * t103 + t161 * t105 + 0.2e1 * (t9 / 0.2e1 + t10 * t455 + t19 * t458) * m(7) + t30;];
tauc  = t14(:);
