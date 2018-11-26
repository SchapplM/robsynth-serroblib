% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 17:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:50:36
% EndTime: 2018-11-23 17:50:45
% DurationCPUTime: 8.66s
% Computational Cost: add. (24500->615), mult. (64871->835), div. (0->0), fcn. (49846->10), ass. (0->299)
t295 = sin(qJ(3));
t296 = sin(qJ(2));
t299 = cos(qJ(3));
t300 = cos(qJ(2));
t265 = -t295 * t296 + t299 * t300;
t253 = t265 * qJD(1);
t266 = t295 * t300 + t299 * t296;
t254 = t266 * qJD(1);
t291 = sin(pkin(11));
t292 = cos(pkin(11));
t209 = t253 * t291 + t254 * t292;
t294 = sin(qJ(5));
t298 = cos(qJ(5));
t335 = t292 * t253 - t254 * t291;
t437 = t209 * t298 + t294 * t335;
t385 = Ifges(6,4) * t437;
t438 = -t209 * t294 + t298 * t335;
t142 = Ifges(6,4) * t438;
t143 = qJD(6) - t438;
t97 = pkin(5) * t437 - pkin(10) * t438;
t290 = qJD(2) + qJD(3);
t440 = t290 / 0.2e1;
t406 = t335 / 0.2e1;
t201 = pkin(9) * t209;
t439 = Ifges(5,4) * t209;
t414 = -pkin(8) - pkin(7);
t279 = t414 * t296;
t280 = t414 * t300;
t227 = t299 * t279 + t280 * t295;
t289 = qJD(5) + t290;
t293 = sin(qJ(6));
t297 = cos(qJ(6));
t131 = t289 * t297 - t293 * t437;
t329 = mrSges(7,1) * t293 + mrSges(7,2) * t297;
t249 = t254 * qJ(4);
t271 = qJD(1) * t280;
t255 = t295 * t271;
t270 = qJD(1) * t279;
t261 = qJD(2) * pkin(2) + t270;
t423 = t299 * t261 + t255;
t184 = t423 - t249;
t173 = pkin(3) * t290 + t184;
t258 = t299 * t271;
t219 = t261 * t295 - t258;
t362 = qJ(4) * t253;
t185 = t219 + t362;
t352 = t292 * t185;
t116 = t291 * t173 + t352;
t394 = pkin(9) * t335;
t100 = t116 + t394;
t174 = t291 * t185;
t115 = t292 * t173 - t174;
t98 = pkin(4) * t290 + t115 - t201;
t49 = -t100 * t294 + t298 * t98;
t47 = -pkin(5) * t289 - t49;
t311 = t47 * t329;
t324 = Ifges(7,5) * t297 - Ifges(7,6) * t293;
t383 = Ifges(7,4) * t297;
t326 = -Ifges(7,2) * t293 + t383;
t384 = Ifges(7,4) * t293;
t328 = Ifges(7,1) * t297 - t384;
t398 = t297 / 0.2e1;
t399 = -t293 / 0.2e1;
t132 = t289 * t293 + t297 * t437;
t408 = t132 / 0.2e1;
t375 = t132 * Ifges(7,4);
t60 = t131 * Ifges(7,2) + t143 * Ifges(7,6) + t375;
t130 = Ifges(7,4) * t131;
t61 = t132 * Ifges(7,1) + t143 * Ifges(7,5) + t130;
t436 = t131 * t326 / 0.2e1 + t328 * t408 + t143 * t324 / 0.2e1 + t311 + t60 * t399 + t61 * t398;
t435 = t439 / 0.2e1 + t406 * Ifges(5,2) + Ifges(5,6) * t440;
t433 = mrSges(5,3) * t335;
t432 = Ifges(5,4) * t335;
t120 = t292 * t184 - t174;
t102 = -t201 + t120;
t283 = pkin(3) * t292 + pkin(4);
t395 = pkin(3) * t291;
t248 = t294 * t283 + t298 * t395;
t119 = -t184 * t291 - t352;
t313 = t119 - t394;
t431 = -t248 * qJD(5) + t102 * t294 - t298 * t313;
t246 = t283 * t298 - t294 * t395;
t235 = t246 * qJD(5);
t52 = t298 * t102 + t294 * t313;
t430 = -t52 + t235;
t222 = -t270 * t295 + t258;
t187 = t222 - t362;
t223 = t299 * t270 + t255;
t188 = -t249 + t223;
t128 = t291 * t187 + t292 * t188;
t104 = -t201 + t128;
t284 = pkin(2) * t299 + pkin(3);
t353 = t291 * t295;
t245 = -pkin(2) * t353 + t292 * t284;
t239 = pkin(4) + t245;
t351 = t292 * t295;
t247 = pkin(2) * t351 + t284 * t291;
t198 = t294 * t239 + t298 * t247;
t376 = pkin(2) * qJD(3);
t243 = (-t291 * t299 - t351) * t376;
t244 = (t292 * t299 - t353) * t376;
t127 = t292 * t187 - t188 * t291;
t312 = t127 - t394;
t429 = -qJD(5) * t198 + (t243 - t312) * t298 + (t104 - t244) * t294;
t197 = t239 * t298 - t247 * t294;
t137 = qJD(5) * t197 + t243 * t294 + t244 * t298;
t54 = t298 * t104 + t294 * t312;
t428 = -t54 + t137;
t367 = -mrSges(6,1) * t289 - mrSges(7,1) * t131 + mrSges(7,2) * t132 + mrSges(6,3) * t437;
t425 = -t127 + t243;
t424 = -t128 + t244;
t228 = t295 * t279 - t299 * t280;
t224 = t290 * t265;
t214 = t224 * qJD(1);
t422 = -qJ(4) * t214 - qJD(4) * t254;
t225 = t290 * t266;
t215 = t225 * qJD(1);
t421 = -qJ(4) * t215 + qJD(4) * t253;
t50 = t100 * t298 + t294 * t98;
t48 = pkin(10) * t289 + t50;
t285 = -pkin(2) * t300 - pkin(1);
t278 = qJD(1) * t285;
t226 = -t253 * pkin(3) + qJD(4) + t278;
t164 = -pkin(4) * t335 + t226;
t78 = -pkin(5) * t438 - pkin(10) * t437 + t164;
t20 = -t293 * t48 + t297 * t78;
t21 = t293 * t78 + t297 * t48;
t420 = -t20 * t293 + t21 * t297;
t154 = -t214 * t291 - t215 * t292;
t350 = qJD(1) * t296;
t287 = pkin(2) * t350;
t191 = pkin(3) * t215 + qJD(2) * t287;
t110 = -pkin(4) * t154 + t191;
t155 = t214 * t292 - t215 * t291;
t74 = qJD(5) * t438 + t154 * t294 + t155 * t298;
t75 = qJD(5) * t437 - t298 * t154 + t155 * t294;
t17 = pkin(5) * t75 - pkin(10) * t74 + t110;
t341 = qJD(2) * t414;
t334 = qJD(1) * t341;
t262 = t296 * t334;
t263 = t300 * t334;
t315 = -t262 * t295 + t263 * t299;
t316 = t262 * t299 + t263 * t295;
t66 = -t291 * (t316 + t421) + t292 * (t315 + t422) + (-t219 * t292 - t291 * t423) * qJD(3);
t302 = -t155 * pkin(9) + t66;
t159 = qJD(3) * t423 + t316;
t160 = -qJD(3) * t219 + t315;
t67 = t292 * (t159 + t421) + t291 * (t160 + t422);
t45 = pkin(9) * t154 + t67;
t8 = qJD(5) * t49 + t294 * t302 + t298 * t45;
t2 = qJD(6) * t20 + t17 * t293 + t297 * t8;
t359 = qJD(6) * t21;
t3 = t17 * t297 - t293 * t8 - t359;
t41 = qJD(6) * t131 + t297 * t74;
t42 = -qJD(6) * t132 - t293 * t74;
t419 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t41 + Ifges(7,6) * t42;
t322 = t20 * t297 + t21 * t293;
t382 = Ifges(6,5) * t289;
t388 = Ifges(6,1) * t437;
t91 = t142 + t382 + t388;
t418 = t164 * mrSges(6,2) + t91 / 0.2e1 + t382 / 0.2e1 + t142 / 0.2e1 - t322 * mrSges(7,3) + t436;
t417 = t41 / 0.2e1;
t416 = t42 / 0.2e1;
t415 = t75 / 0.2e1;
t413 = pkin(1) * mrSges(3,1);
t412 = pkin(1) * mrSges(3,2);
t204 = -qJ(4) * t266 + t227;
t205 = qJ(4) * t265 + t228;
t135 = t292 * t204 - t205 * t291;
t221 = t265 * t291 + t266 * t292;
t111 = -pkin(9) * t221 + t135;
t136 = t291 * t204 + t292 * t205;
t220 = t265 * t292 - t266 * t291;
t112 = pkin(9) * t220 + t136;
t319 = t298 * t111 - t112 * t294;
t9 = qJD(5) * t50 + t294 * t45 - t298 * t302;
t411 = t319 * t9;
t410 = -t131 / 0.2e1;
t409 = -t132 / 0.2e1;
t407 = -t143 / 0.2e1;
t405 = t209 / 0.2e1;
t403 = t253 / 0.2e1;
t402 = t254 / 0.2e1;
t397 = m(4) * t278;
t396 = pkin(3) * t254;
t393 = t2 * t297;
t392 = t20 * mrSges(7,3);
t391 = t3 * t293;
t390 = t49 * mrSges(6,3);
t389 = mrSges(4,3) * t253;
t387 = Ifges(3,4) * t296;
t386 = Ifges(4,4) * t254;
t381 = Ifges(7,5) * t132;
t380 = Ifges(6,2) * t438;
t379 = Ifges(6,6) * t289;
t378 = Ifges(7,6) * t131;
t377 = Ifges(7,3) * t143;
t374 = t437 * t50;
t371 = t254 * mrSges(4,3);
t366 = Ifges(3,5) * qJD(2);
t365 = Ifges(3,6) * qJD(2);
t361 = qJD(2) * mrSges(3,1);
t360 = qJD(2) * mrSges(3,2);
t357 = t116 * t209;
t356 = t438 * t293;
t355 = t438 * t297;
t272 = t296 * t341;
t273 = t300 * t341;
t168 = qJD(3) * t227 + t299 * t272 + t295 * t273;
t125 = -qJ(4) * t225 + qJD(4) * t265 + t168;
t169 = -qJD(3) * t228 - t272 * t295 + t299 * t273;
t126 = -qJ(4) * t224 - qJD(4) * t266 + t169;
t84 = t292 * t125 + t291 * t126;
t349 = qJD(1) * t300;
t348 = qJD(2) * t296;
t345 = qJD(6) * t293;
t344 = qJD(6) * t297;
t340 = t366 / 0.2e1;
t339 = -t365 / 0.2e1;
t338 = t75 * mrSges(6,1) + t74 * mrSges(6,2);
t210 = pkin(2) * t348 + pkin(3) * t225;
t337 = -t154 * mrSges(5,1) + t155 * mrSges(5,2);
t18 = mrSges(7,1) * t75 - mrSges(7,3) * t41;
t87 = -mrSges(7,2) * t143 + mrSges(7,3) * t131;
t336 = -qJD(6) * t87 - t18;
t83 = -t125 * t291 + t292 * t126;
t170 = pkin(4) * t209 + t396;
t332 = -t2 * t293 - t297 * t3;
t331 = (-t3 - t359) * mrSges(7,3);
t330 = mrSges(7,1) * t297 - mrSges(7,2) * t293;
t327 = Ifges(7,1) * t293 + t383;
t325 = Ifges(7,2) * t297 + t384;
t323 = Ifges(7,5) * t293 + Ifges(7,6) * t297;
t77 = t111 * t294 + t112 * t298;
t158 = t220 * t294 + t221 * t298;
t237 = -t265 * pkin(3) + t285;
t178 = -t220 * pkin(4) + t237;
t318 = t298 * t220 - t221 * t294;
t85 = -pkin(5) * t318 - t158 * pkin(10) + t178;
t30 = t293 * t85 + t297 * t77;
t29 = -t293 * t77 + t297 * t85;
t88 = mrSges(7,1) * t143 - mrSges(7,3) * t132;
t320 = -t293 * t88 + t297 * t87;
t162 = -t224 * t291 - t225 * t292;
t124 = -pkin(4) * t162 + t210;
t163 = t224 * t292 - t225 * t291;
t314 = -pkin(9) * t163 + t83;
t80 = t170 + t97;
t307 = -qJD(6) * t322 - t391;
t306 = m(7) * (t307 + t393);
t12 = t41 * Ifges(7,4) + t42 * Ifges(7,2) + t75 * Ifges(7,6);
t13 = t41 * Ifges(7,1) + t42 * Ifges(7,4) + t75 * Ifges(7,5);
t305 = -t8 * mrSges(6,2) + mrSges(7,3) * t393 + t293 * t13 / 0.2e1 + t12 * t398 + t327 * t417 + t325 * t416 + t323 * t415 - Ifges(6,6) * t75 + Ifges(6,5) * t74 + (-mrSges(6,1) - t330) * t9 + t436 * qJD(6);
t59 = t377 + t378 + t381;
t90 = t379 + t380 + t385;
t304 = t164 * mrSges(6,1) + t20 * mrSges(7,1) + t59 / 0.2e1 - t90 / 0.2e1 - t385 / 0.2e1 + t381 / 0.2e1 - t379 / 0.2e1 + t378 / 0.2e1 + t377 / 0.2e1 - t21 * mrSges(7,2);
t141 = t209 * Ifges(5,1) + t290 * Ifges(5,5) + t432;
t202 = Ifges(4,2) * t253 + t290 * Ifges(4,6) + t386;
t250 = Ifges(4,4) * t253;
t203 = t254 * Ifges(4,1) + t290 * Ifges(4,5) + t250;
t301 = -t21 * (-mrSges(7,2) * t437 - mrSges(7,3) * t356) - t20 * (mrSges(7,1) * t437 - mrSges(7,3) * t355) + t437 * t90 / 0.2e1 - t164 * (mrSges(6,1) * t437 + mrSges(6,2) * t438) + (Ifges(7,3) * t437 + t324 * t438) * t407 + (Ifges(7,5) * t437 + t328 * t438) * t409 + (Ifges(7,6) * t437 + t326 * t438) * t410 - t289 * (Ifges(6,5) * t438 - Ifges(6,6) * t437) / 0.2e1 - t438 * t311 + t438 * t390 - (-Ifges(4,2) * t254 + t203 + t250) * t253 / 0.2e1 - (Ifges(4,5) * t253 + Ifges(5,5) * t335 - Ifges(4,6) * t254 - Ifges(5,6) * t209) * t290 / 0.2e1 - (-Ifges(5,2) * t209 + t141 + t432) * t335 / 0.2e1 - t226 * (mrSges(5,1) * t209 + mrSges(5,2) * t335) + t209 * t435 + t423 * t389 + t115 * t433 + t305 - t209 * (Ifges(5,1) * t335 - t439) / 0.2e1 - (-Ifges(6,2) * t437 + t142 + t91) * t438 / 0.2e1 - (Ifges(6,1) * t438 - t385 + t59) * t437 / 0.2e1 - t278 * (mrSges(4,1) * t254 + mrSges(4,2) * t253) + t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,6) * t154 + Ifges(5,5) * t155 - t159 * mrSges(4,2) + t160 * mrSges(4,1) - t61 * t355 / 0.2e1 + t60 * t356 / 0.2e1 + t202 * t402 + Ifges(4,5) * t214 - Ifges(4,6) * t215 - t254 * (Ifges(4,1) * t253 - t386) / 0.2e1;
t286 = Ifges(3,4) * t349;
t275 = mrSges(3,3) * t349 - t360;
t274 = -mrSges(3,3) * t350 + t361;
t252 = Ifges(3,1) * t350 + t286 + t366;
t251 = t365 + (Ifges(3,2) * t300 + t387) * qJD(1);
t241 = pkin(10) + t248;
t240 = -pkin(5) - t246;
t231 = mrSges(4,1) * t290 - t371;
t230 = -mrSges(4,2) * t290 + t389;
t229 = t287 + t396;
t217 = -mrSges(4,1) * t253 + mrSges(4,2) * t254;
t193 = pkin(10) + t198;
t192 = -pkin(5) - t197;
t190 = mrSges(5,1) * t290 - mrSges(5,3) * t209;
t189 = -mrSges(5,2) * t290 + t433;
t165 = t170 + t287;
t150 = -mrSges(5,1) * t335 + mrSges(5,2) * t209;
t133 = -mrSges(6,2) * t289 + mrSges(6,3) * t438;
t96 = -mrSges(6,1) * t438 + mrSges(6,2) * t437;
t82 = qJD(5) * t158 - t298 * t162 + t163 * t294;
t81 = qJD(5) * t318 + t162 * t294 + t163 * t298;
t79 = t287 + t80;
t71 = Ifges(7,3) * t75;
t58 = pkin(9) * t162 + t84;
t28 = t293 * t97 + t297 * t49;
t27 = -t293 * t49 + t297 * t97;
t26 = t293 * t79 + t297 * t54;
t25 = -t293 * t54 + t297 * t79;
t24 = t293 * t80 + t297 * t52;
t23 = -t293 * t52 + t297 * t80;
t22 = pkin(5) * t82 - pkin(10) * t81 + t124;
t19 = -mrSges(7,2) * t75 + mrSges(7,3) * t42;
t16 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t15 = qJD(5) * t77 + t294 * t58 - t298 * t314;
t14 = qJD(5) * t319 + t294 * t314 + t298 * t58;
t5 = -qJD(6) * t30 - t14 * t293 + t22 * t297;
t4 = qJD(6) * t29 + t14 * t297 + t22 * t293;
t1 = [(t13 * t398 + t12 * t399 + Ifges(6,1) * t74 - Ifges(6,4) * t75 + t110 * mrSges(6,2) + t324 * t415 + t328 * t417 + t326 * t416 + (mrSges(6,3) + t329) * t9 + t332 * mrSges(7,3) + (-t297 * t60 / 0.2e1 + t61 * t399 + t47 * t330 + t325 * t410 + t327 * t409 + t323 * t407 - t420 * mrSges(7,3)) * qJD(6)) * t158 + (t159 * t265 - t160 * t266 - t214 * t227 - t215 * t228 - t219 * t225 - t224 * t423) * mrSges(4,3) + (-t265 * t215 - t225 * t403) * Ifges(4,2) + (t265 * t214 - t266 * t215 + t224 * t403 - t225 * t402) * Ifges(4,4) + t285 * (mrSges(4,1) * t215 + mrSges(4,2) * t214) + (-t319 * t74 - t49 * t81 - t50 * t82 - t75 * t77) * mrSges(6,3) - t319 * t16 - (t71 / 0.2e1 - Ifges(6,4) * t74 + t110 * mrSges(6,1) - t8 * mrSges(6,3) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t75 + t419) * t318 + m(4) * (t159 * t228 + t160 * t227 + t168 * t219 + t169 * t423) + t278 * (mrSges(4,1) * t225 + mrSges(4,2) * t224) + t237 * t337 + t178 * t338 + m(7) * (t15 * t47 + t2 * t30 + t20 * t5 + t21 * t4 + t29 * t3 - t411) + m(6) * (t110 * t178 + t124 * t164 + t14 * t50 - t15 * t49 + t77 * t8 - t411) + m(5) * (t115 * t83 + t116 * t84 + t135 * t66 + t136 * t67 + t191 * t237 + t210 * t226) + (-t115 * t163 + t116 * t162 - t135 * t155 + t136 * t154 + t220 * t67 - t221 * t66) * mrSges(5,3) + (Ifges(4,5) * t224 + Ifges(5,5) * t163 - Ifges(4,6) * t225 + Ifges(5,6) * t162) * t440 + (t388 / 0.2e1 + t418) * t81 + (-pkin(7) * t274 + t252 / 0.2e1 + t340 + (-0.2e1 * t412 + 0.3e1 / 0.2e1 * Ifges(3,4) * t300) * qJD(1)) * t300 * qJD(2) + t4 * t87 + t5 * t88 + t29 * t18 + t30 * t19 + t124 * t96 + t14 * t133 + t163 * t141 / 0.2e1 + t367 * t15 + t84 * t189 + t83 * t190 + (-t380 / 0.2e1 + t304) * t82 + t210 * t150 + t191 * (-mrSges(5,1) * t220 + mrSges(5,2) * t221) + t224 * t203 / 0.2e1 - t225 * t202 / 0.2e1 + t226 * (-mrSges(5,1) * t162 + mrSges(5,2) * t163) + t168 * t230 + t169 * t231 + t162 * t435 + (t266 * t214 + t224 * t402) * Ifges(4,1) + (t221 * t155 + t163 * t405) * Ifges(5,1) + (t221 * t154 + t220 * t155 + t162 * t405 + t163 * t406) * Ifges(5,4) + (t220 * t154 + t162 * t406) * Ifges(5,2) + (-pkin(7) * t275 - t251 / 0.2e1 + t339 + (-0.2e1 * t413 - 0.3e1 / 0.2e1 * t387 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t300) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t265 + mrSges(4,2) * t266) + t217 + 0.2e1 * t397) * pkin(2)) * t348; ((t230 * t299 - t231 * t295) * qJD(3) + (-t214 * t299 - t215 * t295) * mrSges(4,3)) * pkin(2) + (-t137 * t88 + t193 * t336 + t331) * t293 + t301 + t425 * t190 + t424 * t189 - t26 * t87 - t25 * t88 + t193 * t306 + (t154 * t247 - t155 * t245 + t357) * mrSges(5,3) - t165 * t96 + t219 * t371 + t192 * t16 + (-t197 * t74 - t198 * t75 + t374) * mrSges(6,3) + (t137 * t87 + t193 * t19 + (-t193 * t88 - t392) * qJD(6)) * t297 - t229 * t150 - t223 * t230 - t222 * t231 + t428 * t133 + ((t340 - t286 / 0.2e1 - t252 / 0.2e1 + qJD(1) * t412 + (t274 - t361) * pkin(7)) * t300 + (t339 + t251 / 0.2e1 + (t413 + t387 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t300) * qJD(1) + (t275 + t360) * pkin(7) + (-t217 - t397) * pkin(2)) * t296) * qJD(1) - t367 * t429 + (t137 * t420 + t192 * t9 - t20 * t25 - t21 * t26 - t429 * t47) * m(7) + (-t164 * t165 - t197 * t9 + t198 * t8 + t428 * t50 + t429 * t49) * m(6) + (t115 * t425 + t116 * t424 - t226 * t229 + t245 * t66 + t247 * t67) * m(5) + (-t423 * t222 - t219 * t223 + (t159 * t295 + t160 * t299 + (t219 * t299 - t295 * t423) * qJD(3)) * pkin(2)) * m(4); (t231 + t371) * t219 + t301 + t430 * t133 + t241 * t306 + (-t235 * t88 + t241 * t336 + t331) * t293 + (t241 * t19 + t235 * t87 + (-t241 * t88 - t392) * qJD(6)) * t297 + (t357 + (t154 * t291 - t155 * t292) * pkin(3)) * mrSges(5,3) + (-t246 * t74 - t248 * t75 + t374) * mrSges(6,3) - t24 * t87 - t23 * t88 - t150 * t396 - t170 * t96 - t120 * t189 - t119 * t190 - t423 * t230 + t240 * t16 - t367 * t431 + (-t20 * t23 - t21 * t24 + t235 * t420 + t240 * t9 - t431 * t47) * m(7) + (-t164 * t170 - t246 * t9 + t248 * t8 + t430 * t50 + t431 * t49) * m(6) + (-t115 * t119 - t116 * t120 - t226 * t396 + (t291 * t67 + t292 * t66) * pkin(3)) * m(5); t297 * t18 - t335 * t189 + t293 * t19 + t209 * t190 - t367 * t437 + t320 * qJD(6) + (-t133 - t320) * t438 + t337 + t338 + (t143 * t420 - t437 * t47 - t332) * m(7) + (t437 * t49 - t438 * t50 + t110) * m(6) + (t115 * t209 - t116 * t335 + t191) * m(5); (-t293 * t18 + t297 * t19 - t88 * t344 - t87 * t345) * pkin(10) + t305 + (t50 * mrSges(6,3) - t304) * t437 + (t390 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t437 - t418) * t438 + t307 * mrSges(7,3) - t367 * t50 - t28 * t87 - t27 * t88 - pkin(5) * t16 - t49 * t133 + ((-t20 * t344 - t21 * t345 - t391 + t393) * pkin(10) - t20 * t27 - t21 * t28 - t47 * t50 - pkin(5) * t9) * m(7); t71 - t47 * (mrSges(7,1) * t132 + mrSges(7,2) * t131) + (Ifges(7,1) * t131 - t375) * t409 + t60 * t408 + (Ifges(7,5) * t131 - Ifges(7,6) * t132) * t407 - t20 * t87 + t21 * t88 + (t131 * t20 + t132 * t21) * mrSges(7,3) + (-Ifges(7,2) * t132 + t130 + t61) * t410 + t419;];
tauc  = t1(:);
