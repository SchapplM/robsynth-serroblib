% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPR6
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
% Datum: 2018-11-23 17:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:03:43
% EndTime: 2018-11-23 17:04:06
% DurationCPUTime: 24.17s
% Computational Cost: add. (12382->730), mult. (37316->969), div. (0->0), fcn. (29275->10), ass. (0->345)
t244 = sin(qJ(2));
t247 = cos(qJ(2));
t241 = sin(pkin(6));
t345 = qJD(1) * t241;
t360 = sin(pkin(11));
t299 = t360 * t345;
t361 = cos(pkin(11));
t300 = t361 * t345;
t198 = -t244 * t299 + t247 * t300;
t195 = qJD(4) - t198;
t402 = t195 / 0.2e1;
t362 = cos(pkin(6));
t333 = pkin(1) * t362;
t236 = t247 * t333;
t231 = qJD(1) * t236;
t387 = pkin(8) + qJ(3);
t321 = t387 * t244;
t306 = t241 * t321;
t264 = pkin(2) * t362 - t306;
t166 = qJD(2) * pkin(2) + qJD(1) * t264 + t231;
t235 = t244 * t333;
t353 = t241 * t247;
t446 = t353 * t387 + t235;
t182 = t446 * qJD(1);
t316 = t361 * t182;
t107 = t360 * t166 + t316;
t277 = qJD(1) * t362 + qJD(2);
t100 = pkin(9) * t277 + t107;
t243 = sin(qJ(4));
t246 = cos(qJ(4));
t199 = -t244 * t300 - t247 * t299;
t221 = (-pkin(2) * t247 - pkin(1)) * t241;
t346 = qJD(1) * t221;
t215 = qJD(3) + t346;
t256 = -pkin(3) * t198 + pkin(9) * t199 + t215;
t59 = t100 * t243 - t246 * t256;
t461 = -qJD(5) - t59;
t45 = -pkin(4) * t195 - t461;
t164 = -t199 * t243 - t246 * t277;
t409 = -t164 / 0.2e1;
t165 = -t246 * t199 + t243 * t277;
t163 = qJD(6) + t165;
t410 = t163 / 0.2e1;
t242 = sin(qJ(6));
t245 = cos(qJ(6));
t119 = t164 * t242 + t195 * t245;
t413 = t119 / 0.2e1;
t118 = t164 * t245 - t195 * t242;
t415 = t118 / 0.2e1;
t488 = Ifges(5,4) * t409 + Ifges(7,5) * t413 + Ifges(7,6) * t415 + Ifges(7,3) * t410;
t406 = t165 / 0.2e1;
t490 = Ifges(5,1) * t406;
t421 = pkin(4) + pkin(10);
t276 = pkin(5) * t165 + t59;
t479 = qJD(5) + t276;
t34 = -t195 * t421 + t479;
t172 = t360 * t182;
t106 = t166 * t361 - t172;
t99 = -pkin(3) * t277 - t106;
t250 = -t165 * qJ(5) + t99;
t39 = t164 * t421 + t250;
t5 = -t242 * t39 + t245 * t34;
t6 = t242 * t34 + t245 * t39;
t485 = t45 * mrSges(6,1) + t5 * mrSges(7,1) - t6 * mrSges(7,2) + t59 * mrSges(5,3) + Ifges(5,5) * t402 + t488 + t490;
t58 = t164 * pkin(4) + t250;
t491 = t99 * mrSges(5,2) - t58 * mrSges(6,3) + t485 + t488;
t472 = Ifges(5,6) - Ifges(6,5);
t487 = Ifges(6,6) + Ifges(5,4);
t181 = -qJD(1) * t306 + t231;
t134 = t181 * t360 + t316;
t341 = qJD(4) * t243;
t358 = t198 * t243;
t486 = -qJD(5) * t243 - t134 + (t341 - t358) * pkin(4);
t344 = qJD(2) * t241;
t452 = -t360 * t244 + t361 * t247;
t260 = t452 * t344;
t258 = qJD(1) * t260;
t457 = qJD(4) * t277 + t258;
t116 = -t199 * t341 - t246 * t457;
t420 = -t116 / 0.2e1;
t340 = qJD(4) * t246;
t117 = -t199 * t340 + t243 * t457;
t418 = -t117 / 0.2e1;
t204 = (t244 * t361 + t247 * t360) * t241;
t200 = qJD(2) * t204;
t193 = qJD(1) * t200;
t404 = t193 / 0.2e1;
t473 = Ifges(6,4) - Ifges(5,5);
t471 = Ifges(5,3) + Ifges(6,1);
t484 = t486 + t195 * (pkin(10) * t243 - qJ(5) * t246);
t135 = t181 * t361 - t172;
t125 = t243 * t135;
t330 = t244 * t345;
t310 = pkin(2) * t330;
t147 = -pkin(3) * t199 - pkin(9) * t198 + t310;
t331 = t360 * pkin(2);
t238 = t331 + pkin(9);
t388 = pkin(5) + t238;
t220 = t388 * t246;
t483 = qJD(4) * t220 - t125 - (pkin(5) * t198 - t147) * t246 - t421 * t199;
t60 = t246 * t100 + t243 * t256;
t48 = -t195 * qJ(5) - t60;
t448 = t48 * mrSges(6,1) - t60 * mrSges(5,3);
t481 = -t99 * mrSges(5,1) + t58 * mrSges(6,2) - t448;
t419 = t116 / 0.2e1;
t417 = t117 / 0.2e1;
t228 = qJD(2) * t231;
t263 = (-qJD(2) * t321 + qJD(3) * t247) * t241;
t158 = qJD(1) * t263 + t228;
t354 = t241 * t244;
t168 = -qJD(2) * t446 - qJD(3) * t354;
t249 = qJD(1) * t168;
t95 = t158 * t361 + t249 * t360;
t480 = qJD(4) * t256 + t95;
t94 = t158 * t360 - t361 * t249;
t262 = t116 * qJ(5) - t165 * qJD(5) + t94;
t31 = t117 * pkin(4) + t262;
t52 = -qJD(6) * t119 + t117 * t245 - t193 * t242;
t428 = t52 / 0.2e1;
t51 = qJD(6) * t118 + t117 * t242 + t193 * t245;
t429 = t51 / 0.2e1;
t21 = t117 * t421 + t262;
t328 = t244 * t344;
t305 = qJD(1) * t328;
t278 = pkin(2) * t305;
t131 = t193 * pkin(3) - pkin(9) * t258 + t278;
t24 = -t100 * t340 + t131 * t246 - t243 * t480;
t7 = -pkin(5) * t116 - t193 * t421 - t24;
t1 = qJD(6) * t5 + t21 * t245 + t242 * t7;
t2 = -qJD(6) * t6 - t21 * t242 + t245 * t7;
t445 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t475 = -Ifges(6,4) / 0.2e1;
t9 = Ifges(7,5) * t51 + Ifges(7,6) * t52 - Ifges(7,3) * t116;
t478 = mrSges(5,2) * t94 - mrSges(6,3) * t31 + Ifges(7,5) * t429 + Ifges(7,6) * t428 + t193 * t475 + t445 - Ifges(6,2) * t419 + t9 / 0.2e1 + (0.2e1 * Ifges(5,1) + Ifges(7,3) + Ifges(6,2)) * t420 + t487 * t418 + (-t473 + Ifges(5,5)) * t404;
t476 = -0.2e1 * Ifges(6,3) * t417 - t94 * mrSges(5,1) + t31 * mrSges(6,2) - t487 * t419 - (t417 - t418) * Ifges(5,2) + 0.2e1 * t404 * t472;
t332 = t361 * pkin(2);
t239 = -t332 - pkin(3);
t270 = -t243 * qJ(5) + t239;
t205 = -t246 * t421 + t270;
t219 = t388 * t243;
t157 = t205 * t245 + t219 * t242;
t470 = -qJD(6) * t157 - t242 * t484 + t245 * t483;
t156 = -t205 * t242 + t219 * t245;
t469 = qJD(6) * t156 + t242 * t483 + t245 * t484;
t468 = -t164 * t472 - t165 * t473 + t195 * t471;
t194 = Ifges(4,4) * t198;
t466 = Ifges(4,2) * t198;
t465 = t277 * Ifges(4,5);
t464 = t277 * Ifges(4,6);
t72 = t246 * t135 + t243 * t147;
t61 = qJ(5) * t199 - t72;
t463 = pkin(5) * t358 - t388 * t341 + t61;
t357 = t198 * t246;
t462 = (-t340 + t357) * qJ(5) + t486;
t369 = t199 * mrSges(4,3);
t460 = -mrSges(4,1) * t277 + mrSges(5,1) * t164 + mrSges(5,2) * t165 - t369;
t152 = t199 * t242 + t245 * t358;
t339 = qJD(6) * t246;
t324 = t242 * t339;
t459 = -t245 * t341 + t152 - t324;
t153 = -t199 * t245 + t242 * t358;
t458 = -t242 * t341 + t245 * t339 + t153;
t456 = qJD(4) * t362 + t260;
t338 = qJD(1) * qJD(2);
t320 = t247 * t338;
t304 = t241 * t320;
t455 = Ifges(3,5) * t304 + Ifges(4,5) * t258 - Ifges(3,6) * t305 - Ifges(4,6) * t193;
t454 = t99 * (mrSges(5,1) * t243 + mrSges(5,2) * t246) + t58 * (-mrSges(6,2) * t243 - mrSges(6,3) * t246);
t453 = -t243 * t472 - t246 * t473;
t451 = t116 * t473 - t117 * t472 + t193 * t471;
t23 = -t100 * t341 + t243 * t131 + t246 * t480;
t450 = t23 * t246 - t24 * t243;
t14 = -qJ(5) * t193 - qJD(5) * t195 - t23;
t20 = -pkin(4) * t193 - t24;
t449 = -t14 * t246 + t20 * t243;
t298 = t1 * t242 + t2 * t245;
t80 = t195 * Ifges(6,5) - t165 * Ifges(6,6) + t164 * Ifges(6,3);
t83 = t165 * Ifges(5,4) - t164 * Ifges(5,2) + t195 * Ifges(5,6);
t447 = t83 / 0.2e1 - t80 / 0.2e1;
t78 = -mrSges(7,2) * t163 + mrSges(7,3) * t118;
t79 = mrSges(7,1) * t163 - mrSges(7,3) * t119;
t281 = t242 * t79 - t245 * t78;
t32 = -mrSges(7,1) * t116 - mrSges(7,3) * t51;
t33 = mrSges(7,2) * t116 + mrSges(7,3) * t52;
t282 = t242 * t33 + t245 * t32;
t89 = -t116 * mrSges(6,1) + t193 * mrSges(6,2);
t442 = t281 * qJD(6) - t282 - t89;
t441 = t215 * mrSges(4,2) + t465 / 0.2e1;
t439 = t24 * mrSges(5,1) - t23 * mrSges(5,2) + t20 * mrSges(6,2) - t14 * mrSges(6,3);
t201 = -pkin(8) * t305 + t228;
t347 = pkin(8) * t353 + t235;
t211 = t347 * qJD(2);
t202 = qJD(1) * t211;
t438 = -t202 * mrSges(3,1) - t94 * mrSges(4,1) - t201 * mrSges(3,2) - t95 * mrSges(4,2);
t286 = Ifges(7,5) * t242 + Ifges(7,6) * t245;
t378 = Ifges(7,4) * t242;
t289 = Ifges(7,2) * t245 + t378;
t377 = Ifges(7,4) * t245;
t293 = Ifges(7,1) * t242 + t377;
t296 = mrSges(7,1) * t245 - mrSges(7,2) * t242;
t297 = t242 * t5 - t245 * t6;
t394 = t164 * pkin(5);
t36 = -t48 - t394;
t379 = Ifges(7,4) * t119;
t43 = Ifges(7,2) * t118 + Ifges(7,6) * t163 + t379;
t364 = t245 * t43;
t113 = Ifges(7,4) * t118;
t44 = Ifges(7,1) * t119 + Ifges(7,5) * t163 + t113;
t365 = t242 * t44;
t411 = -t163 / 0.2e1;
t414 = -t119 / 0.2e1;
t416 = -t118 / 0.2e1;
t437 = t286 * t411 + t289 * t416 + t293 * t414 + t36 * t296 + t297 * mrSges(7,3) - t365 / 0.2e1 - t364 / 0.2e1;
t435 = t215 * mrSges(4,1) + t45 * mrSges(6,2) - t464 / 0.2e1 - t48 * mrSges(6,3) - t59 * mrSges(5,1) - t60 * mrSges(5,2);
t433 = t241 ^ 2;
t11 = Ifges(7,1) * t51 + Ifges(7,4) * t52 - Ifges(7,5) * t116;
t431 = t11 / 0.2e1;
t430 = t43 / 0.2e1;
t160 = Ifges(6,6) * t164;
t82 = t195 * Ifges(6,4) - t165 * Ifges(6,2) + t160;
t425 = -t82 / 0.2e1;
t424 = t82 / 0.2e1;
t423 = -t83 / 0.2e1;
t408 = t164 / 0.2e1;
t407 = -t165 / 0.2e1;
t403 = -t195 / 0.2e1;
t401 = -t198 / 0.2e1;
t400 = -t199 / 0.2e1;
t19 = -mrSges(7,1) * t52 + mrSges(7,2) * t51;
t88 = mrSges(6,1) * t117 - mrSges(6,3) * t193;
t386 = t19 - t88;
t385 = mrSges(4,3) * t193;
t384 = Ifges(3,4) * t244;
t383 = Ifges(3,4) * t247;
t382 = Ifges(4,4) * t199;
t381 = Ifges(5,4) * t243;
t380 = Ifges(5,4) * t246;
t376 = Ifges(6,6) * t243;
t375 = Ifges(6,6) * t246;
t370 = t198 * mrSges(4,3);
t127 = mrSges(6,1) * t164 - mrSges(6,3) * t195;
t70 = -mrSges(7,1) * t118 + mrSges(7,2) * t119;
t363 = t70 - t127;
t359 = qJ(5) * t164;
t352 = t242 * t246;
t351 = t245 * t246;
t178 = t236 + t264;
t197 = qJ(3) * t353 + t347;
t143 = t360 * t178 + t361 * t197;
t133 = pkin(9) * t362 + t143;
t203 = t452 * t241;
t151 = -pkin(3) * t203 - pkin(9) * t204 + t221;
t75 = t246 * t133 + t243 * t151;
t129 = -mrSges(5,2) * t195 - mrSges(5,3) * t164;
t350 = -t127 + t129;
t128 = mrSges(6,1) * t165 + mrSges(6,2) * t195;
t130 = mrSges(5,1) * t195 - mrSges(5,3) * t165;
t349 = -t128 + t130;
t329 = t247 * t345;
t327 = t247 * t344;
t326 = t238 * t341;
t325 = t238 * t340;
t313 = t193 * mrSges(4,1) + mrSges(4,2) * t258;
t74 = -t243 * t133 + t151 * t246;
t71 = t147 * t246 - t125;
t309 = pkin(2) * t328;
t308 = mrSges(3,3) * t330;
t307 = mrSges(3,3) * t329;
t301 = t327 / 0.2e1;
t65 = qJ(5) * t203 - t75;
t295 = Ifges(5,1) * t246 - t381;
t294 = Ifges(7,1) * t245 - t378;
t292 = -Ifges(5,2) * t243 + t380;
t290 = -Ifges(7,2) * t242 + t377;
t287 = Ifges(7,5) * t245 - Ifges(7,6) * t242;
t285 = -Ifges(6,2) * t246 + t376;
t284 = Ifges(6,3) * t243 - t375;
t176 = t246 * t204 + t243 * t362;
t38 = pkin(5) * t176 + t203 * t421 - t74;
t175 = t204 * t243 - t246 * t362;
t142 = t178 * t361 - t360 * t197;
t132 = -pkin(3) * t362 - t142;
t259 = -t176 * qJ(5) + t132;
t57 = t175 * t421 + t259;
t15 = -t242 * t57 + t245 * t38;
t16 = t242 * t38 + t245 * t57;
t232 = qJD(2) * t236;
t167 = t232 + t263;
t104 = t167 * t360 - t361 * t168;
t149 = t175 * t245 + t203 * t242;
t150 = t175 * t242 - t203 * t245;
t105 = t167 * t361 + t168 * t360;
t148 = t200 * pkin(3) - pkin(9) * t260 + t309;
t30 = -t243 * t105 - t133 * t340 + t148 * t246 - t151 * t341;
t29 = t246 * t105 - t133 * t341 + t243 * t148 + t151 * t340;
t269 = pkin(1) * t433 * (mrSges(3,1) * t244 + mrSges(3,2) * t247);
t268 = t244 * t433 * (Ifges(3,1) * t247 - t384);
t140 = t204 * t341 - t246 * t456;
t261 = t140 * qJ(5) - t176 * qJD(5) + t104;
t257 = t260 / 0.2e1;
t22 = -qJ(5) * t200 + qJD(5) * t203 - t29;
t255 = t241 * t277 * (Ifges(3,5) * t247 - Ifges(3,6) * t244);
t254 = (Ifges(3,6) * t362 + (Ifges(3,2) * t247 + t384) * t241) * qJD(1);
t252 = mrSges(4,3) * t258;
t230 = Ifges(3,4) * t329;
t218 = -t246 * pkin(4) + t270;
t216 = -pkin(8) * t354 + t236;
t210 = -pkin(8) * t328 + t232;
t209 = t347 * qJD(1);
t208 = -pkin(8) * t330 + t231;
t207 = -mrSges(3,2) * t277 + t307;
t206 = mrSges(3,1) * t277 - t308;
t180 = Ifges(3,1) * t330 + Ifges(3,5) * t277 + t230;
t179 = Ifges(3,6) * qJD(2) + t254;
t169 = -mrSges(4,2) * t277 + t370;
t154 = -mrSges(4,1) * t198 - mrSges(4,2) * t199;
t146 = -Ifges(4,1) * t199 + t194 + t465;
t145 = -t382 + t464 + t466;
t141 = t204 * t340 + t243 * t456;
t103 = -mrSges(6,2) * t164 - mrSges(6,3) * t165;
t101 = pkin(4) * t165 + t359;
t87 = -mrSges(5,2) * t193 - mrSges(5,3) * t117;
t86 = mrSges(5,1) * t193 + mrSges(5,3) * t116;
t76 = t165 * t421 + t359;
t73 = t175 * pkin(4) + t259;
t69 = mrSges(5,1) * t117 - mrSges(5,2) * t116;
t68 = -mrSges(6,2) * t117 + mrSges(6,3) * t116;
t66 = pkin(4) * t203 - t74;
t64 = -qJD(6) * t150 + t141 * t245 - t200 * t242;
t63 = qJD(6) * t149 + t141 * t242 + t200 * t245;
t62 = pkin(4) * t199 - t71;
t46 = -pkin(5) * t175 - t65;
t41 = t60 - t394;
t35 = t141 * pkin(4) + t261;
t28 = t141 * t421 + t261;
t27 = -pkin(4) * t200 - t30;
t26 = t242 * t41 + t245 * t76;
t25 = -t242 * t76 + t245 * t41;
t13 = -pkin(5) * t141 - t22;
t12 = -pkin(5) * t140 - t200 * t421 - t30;
t10 = Ifges(7,4) * t51 + Ifges(7,2) * t52 - Ifges(7,6) * t116;
t8 = -pkin(5) * t117 - t14;
t4 = -qJD(6) * t16 + t12 * t245 - t242 * t28;
t3 = qJD(6) * t15 + t12 * t242 + t245 * t28;
t17 = [-(t254 + t179) * t328 / 0.2e1 + t438 * t362 + (-Ifges(4,4) * t204 - Ifges(4,2) * t203 - Ifges(4,6) * t362 / 0.2e1) * t193 + ((Ifges(3,5) * t362 + (t244 * Ifges(3,1) + t383) * t241) * t301 + (Ifges(4,1) * t204 + Ifges(4,4) * t203 + Ifges(4,5) * t362) * t257) * qJD(1) + (-t106 * t260 + t203 * t95 + t204 * t94) * mrSges(4,3) + (t468 / 0.2e1 - t466 / 0.2e1 - Ifges(4,4) * t400 + Ifges(5,5) * t406 + Ifges(6,4) * t407 + Ifges(6,5) * t408 + Ifges(5,6) * t409 - t145 / 0.2e1 + t471 * t402 + t435 - t107 * mrSges(4,3)) * t200 + (t201 * t353 + t202 * t354 - t208 * t327 - t209 * t328 - t216 * t304 - t305 * t347) * mrSges(3,3) + m(3) * (t201 * t347 - t202 * t216 - t208 * t211 + t209 * t210) + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t410 + (Ifges(7,5) * t150 + Ifges(7,6) * t149) * t420 - t142 * t252 + m(7) * (t1 * t16 + t13 * t36 + t15 * t2 + t3 * t6 + t4 * t5 + t46 * t8) + m(6) * (t14 * t65 + t20 * t66 + t22 * t48 + t27 * t45 + t31 * t73 + t35 * t58) + (t268 - 0.2e1 * t269) * t338 + t154 * t309 + t46 * t19 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t415 + (Ifges(7,4) * t150 + Ifges(7,2) * t149) * t428 + t146 * t257 + t180 * t301 + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t413 + (Ifges(7,1) * t150 + Ifges(7,4) * t149) * t429 + m(5) * (t132 * t94 + t23 * t75 + t24 * t74 + t29 * t60 - t30 * t59) + m(4) * (t105 * t107 - t142 * t94 + t143 * t95 + (t215 + t346) * t309) + (Ifges(6,2) * t407 + Ifges(6,6) * t408 + t473 * t402 + t424 - t490 - t491) * t140 + qJD(2) * t255 / 0.2e1 + (mrSges(4,2) * t278 + Ifges(4,1) * t258 / 0.2e1) * t204 + (-Ifges(5,4) * t406 + Ifges(6,6) * t407 + Ifges(6,3) * t408 - Ifges(5,2) * t409 + t80 / 0.2e1 + t423 - t472 * t402 - t481) * t141 + (mrSges(6,1) * t20 - mrSges(5,3) * t24 + Ifges(5,4) * t418 - Ifges(6,6) * t417 + t478) * t176 + (t14 * mrSges(6,1) - t23 * mrSges(5,3) - Ifges(5,4) * t420 + Ifges(6,6) * t419 - t476) * t175 + t15 * t32 + t16 * t33 - t143 * t385 + (t1 * t149 - t150 * t2 - t5 * t63 + t6 * t64) * mrSges(7,3) + t221 * t313 + t63 * t44 / 0.2e1 + t36 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t13 * t70 + t73 * t68 + t3 * t78 + t4 * t79 + t433 * (-Ifges(3,2) * t244 + t383) * t320 + t74 * t86 + t75 * t87 + t65 * t88 + t66 * t89 + t35 * t103 - t451 * t203 / 0.2e1 + t455 * t362 / 0.2e1 + (-m(4) * t106 + m(5) * t99 + t460) * t104 + t22 * t127 + t27 * t128 + t29 * t129 + t30 * t130 + t132 * t69 + t149 * t10 / 0.2e1 + t8 * (-mrSges(7,1) * t149 + mrSges(7,2) * t150) + (t194 / 0.2e1 + Ifges(4,1) * t400 + t441) * t260 - (mrSges(4,1) * t278 - Ifges(4,4) * t258 / 0.2e1 + Ifges(6,5) * t417 + Ifges(5,6) * t418 + Ifges(6,4) * t419 + Ifges(5,5) * t420 + t471 * t404 + t439) * t203 + t105 * t169 + t64 * t430 + t150 * t431 + t210 * t207 - t211 * t206; -(-Ifges(3,2) * t330 + t180 + t230) * t329 / 0.2e1 + (t80 + t364 + t365) * t341 / 0.2e1 - (qJD(6) * t44 + t10) * t351 / 0.2e1 + (-t134 * t99 + t238 * t450 + t239 * t94 + t59 * t71 - t60 * t72) * m(5) + (t218 * t31 + t238 * t449 - t45 * t62 + t462 * t58 - t48 * t61) * m(6) + (-t255 / 0.2e1 + (t269 - t268 / 0.2e1) * qJD(1)) * qJD(1) + (-t358 * t48 + t449) * mrSges(6,1) + (t358 * t60 + t450) * mrSges(5,3) + (-t325 - t71) * t130 + t438 + (t325 - t62) * t128 + (t308 + t206) * t209 + (Ifges(7,1) * t153 + Ifges(7,4) * t152) * t414 + (t146 + t194) * t401 + (t307 - t207) * t208 - t376 * t417 + t381 * t418 + (t453 * t402 + (t295 / 0.2e1 - t285 / 0.2e1) * t165 + (-t292 / 0.2e1 + t284 / 0.2e1) * t164 + (Ifges(7,3) * t246 + t243 * t286) * t410 + (Ifges(7,5) * t246 + t243 * t293) * t413 + (Ifges(7,6) * t246 + t243 * t289) * t415 + t454 + (m(5) * (-t60 * t243 + t59 * t246) + m(6) * (t243 * t48 + t246 * t45)) * t238) * qJD(4) + ((t360 * t95 - t361 * t94) * pkin(2) + t106 * t134 - t107 * t135 - t215 * t310) * m(4) - t375 * t419 + t380 * t420 + (-t326 - t72) * t129 + ((t89 - t86) * t238 + t478) * t243 + t455 + (Ifges(7,5) * t153 + Ifges(7,6) * t152) * t411 + ((t87 - t88) * t238 - t286 * t420 - t289 * t428 - t293 * t429 + t8 * t296 + t476) * t246 + (Ifges(7,4) * t153 + Ifges(7,2) * t152) * t416 + (Ifges(7,5) * t414 + Ifges(7,6) * t416 + Ifges(7,3) * t411 + t424 - t485) * t357 + (t425 + t485) * t340 + (t326 - t61) * t127 - t331 * t385 - t107 * t369 - t11 * t352 / 0.2e1 + t106 * t370 + t179 * t330 / 0.2e1 - t252 * t332 + (-t287 * t410 - t290 * t415 - t294 * t413) * t339 - t154 * t310 + t145 * t400 + t447 * t358 + (t423 + t448) * t341 + (t284 * t409 + t285 * t406 + t292 * t408 + t295 * t407 + t403 * t453 - t441 - t454) * t198 + (mrSges(7,1) * t459 - mrSges(7,2) * t458) * t36 + (-t1 * t351 + t2 * t352 + t458 * t5 - t459 * t6) * mrSges(7,3) - t460 * t134 + t462 * t103 + t463 * t70 - t152 * t43 / 0.2e1 - t153 * t44 / 0.2e1 + t156 * t32 + t157 * t33 + (Ifges(4,1) * t198 + t382 + t468) * t199 / 0.2e1 + t469 * t78 + t470 * t79 + (t1 * t157 + t156 * t2 + t220 * t8 + t36 * t463 + t469 * t6 + t470 * t5) * m(7) + (-Ifges(6,4) * t406 - Ifges(5,5) * t407 - Ifges(6,5) * t409 + Ifges(4,2) * t401 - Ifges(5,6) * t408 - t403 * t471 + t435) * t199 - t135 * t169 + t324 * t430 + t218 * t68 + t220 * t19 + t239 * t69; -t152 * t79 - t153 * t78 - t198 * t169 + (t103 + t460) * t199 + (t87 + t349 * t198 + (t242 * t78 + t245 * t79 - t349) * qJD(4) + t386) * t243 + (t86 - t195 * (-t70 - t350) + t442) * t246 + t313 + ((t8 + (t242 * t6 + t245 * t5) * qJD(4)) * t243 + (qJD(4) * t36 + qJD(6) * t297 - t298) * t246 - t152 * t5 - t153 * t6 - t36 * t357) * m(7) + (-t14 * t243 + t199 * t58 - t20 * t246 + t195 * (t243 * t45 - t246 * t48)) * m(6) + (t199 * t99 + t23 * t243 + t24 * t246 + t195 * (t243 * t59 + t246 * t60)) * m(5) + (-t106 * t199 - t107 * t198 + t278) * m(4); (t8 * qJ(5) - t25 * t5 - t26 * t6 + t479 * t36) * m(7) + t437 * qJD(6) + (-t160 / 0.2e1 + t425 + (t475 + Ifges(5,5) / 0.2e1) * t195 + t491) * t164 + t439 + (-pkin(4) * t20 - qJ(5) * t14 - t101 * t58 - t45 * t60 + t461 * t48) * m(6) + t451 + t386 * qJ(5) + t363 * qJD(5) + t349 * t60 + t350 * t59 + t276 * t70 - t298 * mrSges(7,3) - t26 * t78 - t25 * t79 - pkin(4) * t89 - t101 * t103 - (t282 + (-m(7) * t297 - t281) * qJD(6) + m(7) * t298) * t421 + ((Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t165 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t195 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t164 + t437 + t447 + t481) * t165 + t287 * t420 + t290 * t428 + t294 * t429 + t245 * t431 - t242 * t10 / 0.2e1 + t8 * (mrSges(7,1) * t242 + mrSges(7,2) * t245); -t363 * t195 + (t103 - t281) * t165 + (-t163 * t297 - t195 * t36 + t298) * m(7) + (t165 * t58 + t195 * t48 + t20) * m(6) - t442; -t36 * (mrSges(7,1) * t119 + mrSges(7,2) * t118) + (Ifges(7,1) * t118 - t379) * t414 + t43 * t413 + (Ifges(7,5) * t118 - Ifges(7,6) * t119) * t411 - t5 * t78 + t6 * t79 + (t118 * t5 + t119 * t6) * mrSges(7,3) + t9 + (-Ifges(7,2) * t119 + t113 + t44) * t416 + t445;];
tauc  = t17(:);
