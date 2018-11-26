% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:23:10
% EndTime: 2018-11-23 16:23:31
% DurationCPUTime: 22.08s
% Computational Cost: add. (19375->718), mult. (64249->951), div. (0->0), fcn. (54284->12), ass. (0->343)
t248 = sin(qJ(3));
t243 = sin(pkin(12));
t245 = sin(pkin(6));
t344 = qJD(1) * t245;
t330 = t243 * t344;
t251 = cos(qJ(3));
t357 = cos(pkin(12));
t358 = cos(pkin(7));
t304 = t358 * t357;
t285 = t251 * t304;
t244 = sin(pkin(7));
t359 = cos(pkin(6));
t316 = t359 * t244;
t308 = t251 * t316;
t481 = t245 * t285 + t308;
t190 = qJD(1) * t481 - t248 * t330;
t188 = qJD(4) - t190;
t397 = t188 / 0.2e1;
t247 = sin(qJ(4));
t250 = cos(qJ(4));
t331 = pkin(1) * t359;
t240 = t357 * t331;
t236 = qJD(1) * t240;
t353 = t243 * t245;
t269 = t359 * pkin(2) + (-pkin(9) * t358 - qJ(2)) * t353;
t189 = qJD(1) * t269 + t236;
t211 = (-pkin(9) * t243 * t244 - pkin(2) * t357 - pkin(1)) * t245;
t203 = qJD(1) * t211 + qJD(2);
t158 = -t244 * t189 + t203 * t358;
t196 = t248 * t316 + (t243 * t251 + t248 * t304) * t245;
t191 = t196 * qJD(1);
t270 = -t190 * pkin(3) - t191 * pkin(10) + t158;
t239 = t243 * t331;
t318 = t245 * t357;
t307 = qJD(1) * t318;
t216 = qJ(2) * t307 + qJD(1) * t239;
t268 = (t245 * t304 + t316) * pkin(9);
t182 = qJD(1) * t268 + t216;
t172 = t251 * t182;
t315 = t358 * t189;
t351 = t244 * t248;
t119 = t203 * t351 + t248 * t315 + t172;
t273 = t244 * t318 - t358 * t359;
t262 = -qJD(1) * t273 + qJD(3);
t99 = pkin(10) * t262 + t119;
t45 = t247 * t99 - t250 * t270;
t457 = -qJD(5) - t45;
t42 = -pkin(4) * t188 - t457;
t261 = t250 * t262;
t162 = t191 * t247 - t261;
t404 = -t162 / 0.2e1;
t347 = t250 * t191;
t163 = t247 * t262 + t347;
t161 = qJD(6) + t163;
t405 = t161 / 0.2e1;
t246 = sin(qJ(6));
t249 = cos(qJ(6));
t129 = t162 * t246 + t188 * t249;
t408 = t129 / 0.2e1;
t128 = t162 * t249 - t188 * t246;
t410 = t128 / 0.2e1;
t492 = Ifges(5,4) * t404 + Ifges(7,5) * t408 + Ifges(7,6) * t410 + Ifges(7,3) * t405;
t401 = t163 / 0.2e1;
t494 = Ifges(5,1) * t401;
t417 = pkin(4) + pkin(11);
t286 = pkin(5) * t163 + t45;
t479 = qJD(5) + t286;
t31 = -t188 * t417 + t479;
t171 = t248 * t182;
t350 = t244 * t251;
t98 = -pkin(3) * t262 - t203 * t350 - t251 * t315 + t171;
t54 = t162 * pkin(4) - t163 * qJ(5) + t98;
t41 = t162 * pkin(11) + t54;
t5 = -t246 * t41 + t249 * t31;
t6 = t246 * t31 + t249 * t41;
t489 = t42 * mrSges(6,1) + t5 * mrSges(7,1) - t6 * mrSges(7,2) + t45 * mrSges(5,3) + Ifges(5,5) * t397 + t492 + t494;
t495 = t98 * mrSges(5,2) - t54 * mrSges(6,3) + t489 + t492;
t319 = t243 * t358;
t272 = t245 * (-t248 * t319 + t251 * t357);
t209 = qJD(1) * t272;
t480 = qJD(3) * t350 - t209;
t491 = Ifges(6,6) + Ifges(5,4);
t339 = qJD(4) * t247;
t355 = t190 * t247;
t490 = -qJD(5) * t247 - t119 + (t339 - t355) * pkin(4);
t352 = t243 * t248;
t260 = qJD(3) * (t308 + (t285 - t352) * t245);
t258 = qJD(1) * t260;
t126 = -qJD(4) * t261 + t191 * t339 - t250 * t258;
t415 = -t126 / 0.2e1;
t259 = t247 * t260;
t265 = qJD(4) * t273;
t127 = (qJD(3) * t247 + t347) * qJD(4) + (-t247 * t265 + t259) * qJD(1);
t413 = -t127 / 0.2e1;
t192 = t196 * qJD(3);
t180 = qJD(1) * t192;
t399 = t180 / 0.2e1;
t473 = Ifges(6,4) - Ifges(5,5);
t472 = Ifges(5,6) - Ifges(6,5);
t471 = Ifges(5,3) + Ifges(6,1);
t488 = t490 + t188 * (pkin(11) * t247 - qJ(5) * t250);
t281 = t203 * t244 + t315;
t118 = t281 * t251 - t171;
t114 = t247 * t118;
t151 = pkin(3) * t191 - pkin(10) * t190;
t416 = pkin(5) + pkin(10);
t234 = t416 * t250;
t487 = qJD(4) * t234 - t114 - (pkin(5) * t190 - t151) * t250 + t417 * t191;
t221 = t247 * t358 + t250 * t351;
t311 = t244 * t330;
t486 = qJD(4) * t221 + t247 * t480 + t250 * t311;
t271 = t245 * (t248 * t357 + t251 * t319);
t208 = qJD(1) * t271;
t342 = qJD(3) * t248;
t328 = t244 * t342;
t485 = t328 - t208;
t46 = t247 * t270 + t250 * t99;
t43 = -t188 * qJ(5) - t46;
t444 = t43 * mrSges(6,1) - t46 * mrSges(5,3);
t483 = -t98 * mrSges(5,1) + t54 * mrSges(6,2) - t444;
t414 = t126 / 0.2e1;
t412 = t127 / 0.2e1;
t267 = qJD(2) * t272;
t92 = qJD(1) * t267 + qJD(3) * t118;
t482 = qJD(4) * t270 + t92;
t345 = qJ(2) * t318 + t239;
t194 = t268 + t345;
t197 = t240 + t269;
t314 = t358 * t197;
t280 = t211 * t244 + t314;
t134 = -t248 * t194 + t280 * t251;
t57 = qJD(6) * t128 + t127 * t246 + t180 * t249;
t58 = -qJD(6) * t129 + t127 * t249 - t180 * t246;
t13 = Ifges(7,5) * t57 + Ifges(7,6) * t58 - Ifges(7,3) * t126;
t266 = qJD(2) * t271;
t93 = qJD(1) * t266 + (t248 * t281 + t172) * qJD(3);
t253 = t126 * qJ(5) - t163 * qJD(5) + t93;
t32 = t127 * pkin(4) + t253;
t424 = t58 / 0.2e1;
t425 = t57 / 0.2e1;
t25 = t127 * t417 + t253;
t343 = qJD(2) * t245;
t329 = t243 * t343;
t309 = qJD(1) * t329;
t287 = t244 * t309;
t141 = t180 * pkin(3) - pkin(10) * t258 + t287;
t338 = qJD(4) * t250;
t24 = t141 * t250 - t247 * t482 - t99 * t338;
t7 = -pkin(5) * t126 - t180 * t417 - t24;
t1 = qJD(6) * t5 + t246 * t7 + t249 * t25;
t2 = -qJD(6) * t6 - t246 * t25 + t249 * t7;
t440 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t475 = -Ifges(6,4) / 0.2e1;
t478 = mrSges(5,2) * t93 - mrSges(6,3) * t32 + Ifges(7,5) * t425 + Ifges(7,6) * t424 + t180 * t475 + t440 - Ifges(6,2) * t414 + t13 / 0.2e1 + (0.2e1 * Ifges(5,1) + Ifges(7,3) + Ifges(6,2)) * t415 + t491 * t413 + (-t473 + Ifges(5,5)) * t399;
t476 = 0.2e1 * Ifges(6,3) * t412 + t93 * mrSges(5,1) - t180 * Ifges(5,6) / 0.2e1 - t32 * mrSges(6,2) + t491 * t414 + (t412 - t413) * Ifges(5,2) + (-t472 + Ifges(6,5)) * t399;
t474 = t93 * mrSges(4,1);
t469 = -t162 * t472 - t163 * t473 + t188 * t471;
t89 = -t126 * mrSges(6,1) + t180 * mrSges(6,2);
t90 = mrSges(5,1) * t180 + mrSges(5,3) * t126;
t467 = t89 - t90;
t88 = mrSges(6,1) * t127 - mrSges(6,3) * t180;
t91 = -mrSges(5,2) * t180 - mrSges(5,3) * t127;
t466 = t91 - t88;
t465 = Ifges(4,1) * t196;
t187 = Ifges(4,4) * t190;
t464 = Ifges(4,2) * t190;
t463 = t262 * Ifges(4,5);
t462 = t262 * Ifges(4,6);
t320 = -qJ(5) * t247 - pkin(3);
t222 = -t250 * t417 + t320;
t233 = t416 * t247;
t199 = t222 * t249 + t233 * t246;
t461 = -qJD(6) * t199 - t246 * t488 + t249 * t487;
t198 = -t222 * t246 + t233 * t249;
t460 = qJD(6) * t198 + t246 * t487 + t249 * t488;
t354 = t190 * t250;
t459 = (-t338 + t354) * qJ(5) + t490;
t75 = t250 * t118 + t247 * t151;
t65 = -qJ(5) * t191 - t75;
t458 = pkin(5) * t355 - t416 * t339 + t65;
t366 = t191 * mrSges(4,3);
t456 = -mrSges(4,1) * t262 + mrSges(5,1) * t162 + mrSges(5,2) * t163 + t366;
t130 = mrSges(6,1) * t162 - mrSges(6,3) * t188;
t132 = -mrSges(5,2) * t188 - mrSges(5,3) * t162;
t455 = -t130 + t132;
t131 = mrSges(6,1) * t163 + mrSges(6,2) * t188;
t133 = mrSges(5,1) * t188 - mrSges(5,3) * t163;
t454 = -t131 + t133;
t147 = -t191 * t246 + t249 * t355;
t337 = qJD(6) * t250;
t326 = t246 * t337;
t453 = -t249 * t339 + t147 - t326;
t148 = t191 * t249 + t246 * t355;
t452 = -t246 * t339 + t249 * t337 + t148;
t220 = t247 * t351 - t250 * t358;
t284 = -t220 * t246 + t249 * t350;
t451 = qJD(6) * t284 - t246 * t485 + t249 * t486;
t204 = t220 * t249 + t246 * t350;
t450 = qJD(6) * t204 + t246 * t486 + t249 * t485;
t449 = t98 * (mrSges(5,1) * t247 + mrSges(5,2) * t250) + t54 * (-mrSges(6,2) * t247 - mrSges(6,3) * t250);
t448 = -t247 * t472 - t250 * t473;
t447 = t126 * t473 - t127 * t472 + t180 * t471;
t23 = t247 * t141 + t250 * t482 - t339 * t99;
t446 = t23 * t250 - t24 * t247;
t16 = -qJ(5) * t180 - qJD(5) * t188 - t23;
t17 = -pkin(4) * t180 - t24;
t445 = -t16 * t250 + t17 * t247;
t443 = t1 * t246 + t2 * t249;
t80 = t188 * Ifges(6,5) - t163 * Ifges(6,6) + t162 * Ifges(6,3);
t83 = t163 * Ifges(5,4) - t162 * Ifges(5,2) + t188 * Ifges(5,6);
t442 = t83 / 0.2e1 - t80 / 0.2e1;
t441 = -m(5) * t98 - t456;
t437 = t158 * mrSges(4,2) + t463 / 0.2e1;
t435 = -m(4) * t118 - t441;
t433 = t24 * mrSges(5,1) - t23 * mrSges(5,2) + t17 * mrSges(6,2) - t16 * mrSges(6,3);
t293 = Ifges(7,5) * t246 + Ifges(7,6) * t249;
t377 = Ifges(7,4) * t246;
t296 = Ifges(7,2) * t249 + t377;
t376 = Ifges(7,4) * t249;
t300 = Ifges(7,1) * t246 + t376;
t303 = mrSges(7,1) * t249 - mrSges(7,2) * t246;
t305 = t246 * t5 - t249 * t6;
t387 = t162 * pkin(5);
t33 = -t43 - t387;
t372 = t129 * Ifges(7,4);
t51 = t128 * Ifges(7,2) + t161 * Ifges(7,6) + t372;
t362 = t249 * t51;
t125 = Ifges(7,4) * t128;
t52 = t129 * Ifges(7,1) + t161 * Ifges(7,5) + t125;
t363 = t246 * t52;
t406 = -t161 / 0.2e1;
t409 = -t129 / 0.2e1;
t411 = -t128 / 0.2e1;
t431 = t293 * t406 + t296 * t411 + t300 * t409 + t33 * t303 + t305 * mrSges(7,3) - t363 / 0.2e1 - t362 / 0.2e1;
t429 = t43 * mrSges(6,3) + t45 * mrSges(5,1) + t46 * mrSges(5,2) - t158 * mrSges(4,1) + t462 / 0.2e1 - t42 * mrSges(6,2);
t15 = Ifges(7,1) * t57 + Ifges(7,4) * t58 - Ifges(7,5) * t126;
t427 = t15 / 0.2e1;
t426 = t51 / 0.2e1;
t159 = Ifges(6,6) * t162;
t82 = t188 * Ifges(6,4) - t163 * Ifges(6,2) + t159;
t421 = -t82 / 0.2e1;
t420 = t82 / 0.2e1;
t419 = -t83 / 0.2e1;
t403 = t162 / 0.2e1;
t402 = -t163 / 0.2e1;
t398 = -t188 / 0.2e1;
t396 = -t190 / 0.2e1;
t394 = t191 / 0.2e1;
t381 = mrSges(4,3) * t180;
t380 = Ifges(4,4) * t191;
t379 = Ifges(5,4) * t247;
t378 = Ifges(5,4) * t250;
t375 = Ifges(6,6) * t247;
t374 = Ifges(6,6) * t250;
t367 = t190 * mrSges(4,3);
t361 = t251 * t93;
t73 = -mrSges(7,1) * t128 + mrSges(7,2) * t129;
t360 = t73 - t130;
t356 = qJ(5) * t162;
t349 = t246 * t250;
t348 = t249 * t250;
t186 = t251 * t194;
t164 = -t197 * t244 + t358 * t211;
t195 = t245 * t352 - t481;
t110 = pkin(3) * t195 - pkin(10) * t196 + t164;
t135 = t211 * t351 + t248 * t314 + t186;
t117 = -pkin(10) * t273 + t135;
t64 = t247 * t110 + t250 * t117;
t346 = Ifges(4,5) * t258 - Ifges(4,6) * t180;
t336 = pkin(10) * t339;
t335 = pkin(10) * t338;
t63 = t110 * t250 - t247 * t117;
t74 = t151 * t250 - t114;
t310 = t244 * t329;
t48 = -qJ(5) * t195 - t64;
t302 = Ifges(5,1) * t250 - t379;
t301 = Ifges(7,1) * t249 - t377;
t299 = -Ifges(5,2) * t247 + t378;
t297 = -Ifges(7,2) * t246 + t376;
t294 = Ifges(7,5) * t249 - Ifges(7,6) * t246;
t292 = -Ifges(6,2) * t250 + t375;
t291 = Ifges(6,3) * t247 - t374;
t168 = t250 * t196 - t247 * t273;
t35 = pkin(5) * t168 - t195 * t417 - t63;
t167 = t196 * t247 + t250 * t273;
t116 = pkin(3) * t273 - t134;
t255 = -t168 * qJ(5) + t116;
t47 = t167 * t417 + t255;
t9 = -t246 * t47 + t249 * t35;
t10 = t246 * t35 + t249 * t47;
t36 = -mrSges(7,1) * t126 - mrSges(7,3) * t57;
t37 = mrSges(7,2) * t126 + mrSges(7,3) * t58;
t289 = t246 * t37 + t249 * t36;
t78 = -mrSges(7,2) * t161 + mrSges(7,3) * t128;
t79 = mrSges(7,1) * t161 - mrSges(7,3) * t129;
t288 = -t246 * t79 + t249 * t78;
t142 = t167 * t249 - t195 * t246;
t143 = t167 * t246 + t195 * t249;
t102 = qJD(3) * t134 + t267;
t145 = t192 * pkin(3) - pkin(10) * t260 + t310;
t30 = -t247 * t102 - t110 * t339 - t117 * t338 + t145 * t250;
t279 = -(-qJ(2) * t330 + t236) * t243 + t216 * t357;
t29 = t250 * t102 + t110 * t338 - t117 * t339 + t247 * t145;
t218 = (mrSges(3,1) * t359 - mrSges(3,3) * t353) * qJD(1);
t219 = (-mrSges(3,2) * t359 + mrSges(3,3) * t318) * qJD(1);
t22 = -qJ(5) * t192 - qJD(5) * t195 - t29;
t257 = t260 / 0.2e1;
t256 = mrSges(4,3) * t258;
t103 = t266 + (t248 * t280 + t186) * qJD(3);
t139 = t196 * t339 + (-t260 + t265) * t250;
t254 = t139 * qJ(5) - t168 * qJD(5) + t103;
t229 = -pkin(4) * t250 + t320;
t165 = -mrSges(4,2) * t262 + t367;
t150 = -mrSges(4,1) * t190 + mrSges(4,2) * t191;
t146 = t180 * mrSges(4,1) + mrSges(4,2) * t258;
t140 = qJD(4) * t168 + t259;
t138 = Ifges(4,1) * t191 + t187 + t463;
t137 = t380 + t462 + t464;
t108 = -mrSges(6,2) * t162 - mrSges(6,3) * t163;
t106 = pkin(4) * t163 + t356;
t77 = t163 * t417 + t356;
t72 = mrSges(5,1) * t127 - mrSges(5,2) * t126;
t71 = -mrSges(6,2) * t127 + mrSges(6,3) * t126;
t70 = -qJD(6) * t143 + t140 * t249 - t192 * t246;
t69 = qJD(6) * t142 + t140 * t246 + t192 * t249;
t68 = t167 * pkin(4) + t255;
t67 = -pkin(4) * t191 - t74;
t49 = -pkin(4) * t195 - t63;
t40 = -pkin(5) * t167 - t48;
t39 = t46 - t387;
t34 = t140 * pkin(4) + t254;
t28 = t140 * t417 + t254;
t27 = -pkin(4) * t192 - t30;
t26 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t19 = t246 * t39 + t249 * t77;
t18 = -t246 * t77 + t249 * t39;
t14 = t57 * Ifges(7,4) + t58 * Ifges(7,2) - t126 * Ifges(7,6);
t12 = -pkin(5) * t140 - t22;
t11 = -pkin(5) * t139 - t192 * t417 - t30;
t8 = -pkin(5) * t127 - t16;
t4 = -qJD(6) * t10 + t11 * t249 - t246 * t28;
t3 = qJD(6) * t9 + t11 * t246 + t249 * t28;
t20 = [(t16 * mrSges(6,1) - t23 * mrSges(5,3) - Ifges(5,4) * t415 + Ifges(6,6) * t414 + t476) * t167 + (Ifges(7,5) * t69 + Ifges(7,6) * t70) * t405 + (Ifges(7,5) * t143 + Ifges(7,6) * t142) * t415 + (-Ifges(4,4) * t196 + Ifges(4,6) * t273 / 0.2e1) * t180 + (Ifges(4,2) * t180 + t447 / 0.2e1 - t92 * mrSges(4,3) + mrSges(4,1) * t287 + Ifges(6,4) * t414 + Ifges(5,5) * t415 + Ifges(6,5) * t412 + Ifges(5,6) * t413 + t471 * t399 + (-qJD(1) * t257 - t258 / 0.2e1) * Ifges(4,4) + t433) * t195 + m(5) * (t116 * t93 + t23 * t64 + t24 * t63 + t29 * t46 - t30 * t45) + m(4) * (t102 * t119 - t134 * t93 + t135 * t92 + (qJD(1) * t164 + t158) * t310) + (mrSges(6,1) * t17 - mrSges(5,3) * t24 + Ifges(5,4) * t413 - Ifges(6,6) * t412 + t478) * t168 + (Ifges(7,1) * t69 + Ifges(7,4) * t70) * t408 + (Ifges(7,1) * t143 + Ifges(7,4) * t142) * t425 - t134 * t256 + (-t118 * t260 + t196 * t93) * mrSges(4,3) + t150 * t310 + (Ifges(7,4) * t69 + Ifges(7,2) * t70) * t410 + (Ifges(7,4) * t143 + Ifges(7,2) * t142) * t424 + (t1 * t142 - t143 * t2 - t5 * t69 + t6 * t70) * mrSges(7,3) + t70 * t426 + t143 * t427 + (t196 * t287 + t273 * t92) * mrSges(4,2) + (-Ifges(5,4) * t401 + Ifges(6,6) * t402 + Ifges(6,3) * t403 - Ifges(5,2) * t404 + t419 + t80 / 0.2e1 - t472 * t397 - t483) * t140 + (-t119 * mrSges(4,3) + t469 / 0.2e1 + Ifges(5,5) * t401 + Ifges(6,4) * t402 + Ifges(6,5) * t403 + Ifges(5,6) * t404 - Ifges(4,4) * t394 - t464 / 0.2e1 - t137 / 0.2e1 + t471 * t397 - t429) * t192 + t435 * t103 - t135 * t381 - 0.2e1 * t218 * t329 + m(7) * (t1 * t10 + t12 * t33 + t2 * t9 + t3 * t6 + t4 * t5 + t40 * t8) + m(6) * (t16 * t48 + t17 * t49 + t22 * t43 + t27 * t42 + t32 * t68 + t34 * t54) + m(3) * ((t357 * t345 + (qJ(2) * t353 - t240) * t243) * qJD(1) + t279) * t343 + (Ifges(4,1) * t394 + t187 / 0.2e1 + t437) * t260 + t258 * t465 / 0.2e1 + (qJD(1) * (-Ifges(4,5) * t273 + t465) + t138) * t257 + (-t346 / 0.2e1 + t474) * t273 + (Ifges(6,2) * t402 + Ifges(6,6) * t403 + t397 * t473 + t420 - t494 - t495) * t139 + 0.2e1 * qJD(2) * t219 * t318 + t9 * t36 + t10 * t37 + t40 * t26 + t69 * t52 / 0.2e1 + t33 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t68 * t71 + t12 * t73 + t3 * t78 + t4 * t79 + t48 * t88 + t49 * t89 + t63 * t90 + t64 * t91 + t34 * t108 + t116 * t72 + t22 * t130 + t27 * t131 + t29 * t132 + t30 * t133 + t142 * t14 / 0.2e1 + t8 * (-mrSges(7,1) * t142 + mrSges(7,2) * t143) + t164 * t146 + t102 * t165; -m(3) * t279 * t344 + t358 * t146 - t150 * t311 + t204 * t36 - t284 * t37 + t218 * t330 - t219 * t307 - t351 * t381 + t451 * t79 + t450 * t78 + (m(5) * (t342 * t98 - t361) + m(6) * (-t251 * t32 + t342 * t54)) * t244 + (-m(5) * t24 + m(6) * t17 + t467) * t220 + t480 * t165 + (-t1 * t284 + t2 * t204 + t450 * t6 + t451 * t5) * m(7) + (-t119 * t209 - t158 * t311 + (t358 * t309 + t248 * t92 - t361 + (-t118 * t248 + t119 * t251) * qJD(3)) * t244) * m(4) + (-t256 - t72 - t71) * t350 + (t108 + t456) * t328 + (m(5) * t23 - m(6) * t16 + m(7) * t8 + t26 + t466) * t221 + (-m(6) * t54 - t108 - t435) * t208 + t486 * (m(5) * t45 + m(6) * t42 - t454) + (t220 * qJD(4) + t247 * t311 - t250 * t480) * (-m(5) * t46 + m(6) * t43 - m(7) * t33 - t455 - t73); ((Ifges(7,3) * t250 + t247 * t293) * t405 + (Ifges(7,5) * t250 + t247 * t300) * t408 + (Ifges(7,6) * t250 + t247 * t296) * t410 + t449 + (-t299 / 0.2e1 + t291 / 0.2e1) * t162 + t448 * t397 + (t302 / 0.2e1 - t292 / 0.2e1) * t163) * qJD(4) + (t336 - t65) * t130 + (Ifges(7,1) * t148 + Ifges(7,4) * t147) * t409 + (-t293 * t415 - t296 * t424 - t300 * t425 + t8 * t303 - t476) * t250 + (Ifges(7,5) * t148 + Ifges(7,6) * t147) * t406 + t442 * t355 + (t419 + t444) * t339 + (t291 * t404 + t292 * t401 + t299 * t403 + t302 * t402 + t448 * t398 - t437 - t449) * t190 + (mrSges(7,1) * t453 - mrSges(7,2) * t452) * t33 + (-t1 * t348 + t2 * t349 + t452 * t5 - t453 * t6) * mrSges(7,3) + (t367 - t165) * t118 + (-t355 * t43 + t445) * mrSges(6,1) + (t355 * t46 + t446) * mrSges(5,3) + (-t336 - t75) * t132 + (t366 + t441) * t119 - t474 + t346 + t137 * t394 + (t138 + t187) * t396 + (-pkin(3) * t93 + t45 * t74 - t46 * t75) * m(5) + (t229 * t32 - t42 * t67 - t43 * t65 + t459 * t54) * m(6) + (((-t247 * t46 + t250 * t45) * qJD(4) + t446) * m(5) + t466 * t250 + t467 * t247 + ((t247 * t43 + t250 * t42) * qJD(4) + t445) * m(6)) * pkin(10) + (t335 - t67) * t131 + t478 * t247 + (t362 + t363 + t80) * t339 / 0.2e1 - (qJD(6) * t52 + t14) * t348 / 0.2e1 + (-t294 * t405 - t297 * t410 - t301 * t408) * t337 + t326 * t426 + (-t335 - t74) * t133 - t374 * t414 + t378 * t415 + (Ifges(7,4) * t148 + Ifges(7,2) * t147) * t411 - t15 * t349 / 0.2e1 - t375 * t412 + t379 * t413 + t229 * t71 + t234 * t26 + t198 * t36 + t199 * t37 + t458 * t73 + t459 * t108 + t460 * t78 + t461 * t79 + (t1 * t199 + t198 * t2 + t234 * t8 + t33 * t458 + t460 * t6 + t461 * t5) * m(7) + (Ifges(7,5) * t409 + Ifges(7,6) * t411 + Ifges(7,3) * t406 + t420 - t489) * t354 + (t421 + t489) * t338 - (Ifges(4,1) * t190 - t380 + t469) * t191 / 0.2e1 + (Ifges(6,4) * t401 + Ifges(5,5) * t402 + Ifges(6,5) * t404 - Ifges(4,2) * t396 + Ifges(5,6) * t403 + t471 * t398 + t429) * t191 - pkin(3) * t72 - t92 * mrSges(4,2) - t147 * t51 / 0.2e1 - t148 * t52 / 0.2e1; (-t159 / 0.2e1 + t421 + (t475 + Ifges(5,5) / 0.2e1) * t188 + t495) * t162 + t433 - ((-m(7) * t305 + t288) * qJD(6) + m(7) * t443 + t289) * t417 + ((Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t163 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t188 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t162 + t431 + t442 + t483) * t163 + (-pkin(4) * t17 - qJ(5) * t16 - t106 * t54 - t42 * t46 + t43 * t457) * m(6) + t431 * qJD(6) + t447 + (t8 * qJ(5) - t18 * t5 - t19 * t6 + t479 * t33) * m(7) + t294 * t415 + t297 * t424 + t301 * t425 + t249 * t427 + (t26 - t88) * qJ(5) + t360 * qJD(5) + t8 * (mrSges(7,1) * t246 + mrSges(7,2) * t249) - t246 * t14 / 0.2e1 + t455 * t45 + t454 * t46 - t443 * mrSges(7,3) + t286 * t73 - t19 * t78 - t18 * t79 - pkin(4) * t89 - t106 * t108; -t360 * t188 + t288 * qJD(6) + (t108 + t288) * t163 + t289 + t89 + (-t161 * t305 - t188 * t33 + t443) * m(7) + (t163 * t54 + t188 * t43 + t17) * m(6); -t33 * (mrSges(7,1) * t129 + mrSges(7,2) * t128) + (Ifges(7,1) * t128 - t372) * t409 + t51 * t408 + (Ifges(7,5) * t128 - Ifges(7,6) * t129) * t406 - t5 * t78 + t6 * t79 + (t128 * t5 + t129 * t6) * mrSges(7,3) + t13 + (-Ifges(7,2) * t129 + t125 + t52) * t411 + t440;];
tauc  = t20(:);
