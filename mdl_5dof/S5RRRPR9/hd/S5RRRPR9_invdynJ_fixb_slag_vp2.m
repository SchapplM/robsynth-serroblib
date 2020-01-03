% Calculate vector of inverse dynamics joint torques for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:15
% EndTime: 2019-12-31 21:23:05
% DurationCPUTime: 28.12s
% Computational Cost: add. (9227->753), mult. (20879->1039), div. (0->0), fcn. (14672->14), ass. (0->334)
t287 = sin(qJ(2));
t291 = cos(qJ(2));
t328 = pkin(2) * t287 - pkin(7) * t291;
t236 = t328 * qJD(1);
t290 = cos(qJ(3));
t286 = sin(qJ(3));
t364 = qJD(1) * t287;
t346 = t286 * t364;
t162 = pkin(6) * t346 + t290 * t236;
t371 = t290 * t291;
t309 = pkin(3) * t287 - qJ(4) * t371;
t284 = -qJ(4) - pkin(7);
t332 = qJD(3) * t284;
t491 = -qJD(1) * t309 - qJD(4) * t286 + t290 * t332 - t162;
t211 = t286 * t236;
t356 = qJD(4) * t290;
t375 = t287 * t290;
t377 = t286 * t291;
t490 = t211 + (-pkin(6) * t375 - qJ(4) * t377) * qJD(1) - t286 * t332 - t356;
t282 = sin(pkin(9));
t283 = cos(pkin(9));
t227 = t282 * t290 + t283 * t286;
t305 = t227 * t291;
t176 = qJD(1) * t305;
t203 = t227 * qJD(3);
t489 = t176 - t203;
t311 = t282 * t286 - t283 * t290;
t304 = t311 * t291;
t177 = qJD(1) * t304;
t204 = t311 * qJD(3);
t488 = t177 - t204;
t465 = t490 * t282 + t283 * t491;
t464 = t282 * t491 - t490 * t283;
t487 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t281 = qJ(3) + pkin(9);
t274 = cos(t281);
t400 = pkin(3) * t290;
t242 = pkin(4) * t274 + t400;
t235 = pkin(2) + t242;
t275 = qJ(5) + t281;
t262 = sin(t275);
t263 = cos(t275);
t266 = pkin(2) + t400;
t273 = sin(t281);
t324 = -mrSges(4,1) * t290 + mrSges(4,2) * t286;
t486 = m(4) * pkin(2) + m(5) * t266 + m(6) * t235 + t274 * mrSges(5,1) + t263 * mrSges(6,1) - t273 * mrSges(5,2) - t262 * mrSges(6,2) - t324;
t485 = Ifges(4,3) + Ifges(5,3);
t484 = -pkin(4) * t364 - t488 * pkin(8) + t465;
t483 = t489 * pkin(8) + t464;
t325 = mrSges(3,1) * t287 + mrSges(3,2) * t291;
t394 = Ifges(3,4) * t287;
t481 = pkin(1) * t325 - t287 * (Ifges(3,1) * t291 - t394) / 0.2e1;
t280 = -pkin(8) + t284;
t480 = -m(4) * pkin(7) + m(5) * t284 + m(6) * t280 - t487;
t361 = qJD(2) * t290;
t233 = -t346 + t361;
t355 = qJD(1) * qJD(2);
t239 = qJDD(1) * t287 + t291 * t355;
t140 = qJD(3) * t233 + qJDD(2) * t286 + t239 * t290;
t238 = qJDD(1) * t291 - t287 * t355;
t225 = qJDD(3) - t238;
t344 = t290 * t364;
t234 = qJD(2) * t286 + t344;
t329 = pkin(2) * t291 + pkin(7) * t287;
t243 = -pkin(1) - t329;
t218 = t243 * qJD(1);
t363 = qJD(1) * t291;
t271 = pkin(6) * t363;
t250 = qJD(2) * pkin(7) + t271;
t155 = t218 * t286 + t250 * t290;
t382 = qJDD(1) * pkin(1);
t156 = -pkin(2) * t238 - pkin(7) * t239 - t382;
t223 = t238 * pkin(6);
t199 = qJDD(2) * pkin(7) + t223;
t68 = -t155 * qJD(3) + t290 * t156 - t199 * t286;
t42 = pkin(3) * t225 - qJ(4) * t140 - qJD(4) * t234 + t68;
t141 = -qJD(3) * t234 + qJDD(2) * t290 - t239 * t286;
t357 = qJD(3) * t290;
t359 = qJD(3) * t286;
t67 = t286 * t156 + t290 * t199 + t218 * t357 - t250 * t359;
t47 = qJ(4) * t141 + qJD(4) * t233 + t67;
t14 = t282 * t42 + t283 * t47;
t79 = -t140 * t282 + t141 * t283;
t10 = pkin(8) * t79 + t14;
t285 = sin(qJ(5));
t289 = cos(qJ(5));
t258 = qJD(3) - t363;
t149 = t233 * t282 + t234 * t283;
t476 = pkin(8) * t149;
t154 = t290 * t218 - t250 * t286;
t116 = -qJ(4) * t234 + t154;
t104 = pkin(3) * t258 + t116;
t117 = qJ(4) * t233 + t155;
t110 = t282 * t117;
t53 = t283 * t104 - t110;
t44 = pkin(4) * t258 - t476 + t53;
t330 = t283 * t233 - t234 * t282;
t471 = pkin(8) * t330;
t379 = t283 * t117;
t54 = t282 * t104 + t379;
t45 = t54 + t471;
t11 = -t285 * t45 + t289 * t44;
t13 = -t282 * t47 + t283 * t42;
t81 = t140 * t283 + t141 * t282;
t9 = pkin(4) * t225 - pkin(8) * t81 + t13;
t2 = qJD(5) * t11 + t10 * t289 + t285 * t9;
t12 = t285 * t44 + t289 * t45;
t3 = -qJD(5) * t12 - t10 * t285 + t289 * t9;
t479 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t478 = -t68 * mrSges(4,1) - t13 * mrSges(5,1) + t67 * mrSges(4,2) + t14 * mrSges(5,2);
t319 = t291 * Ifges(3,2) + t394;
t477 = t12 * mrSges(6,2) + t54 * mrSges(5,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t319 / 0.2e1 - t11 * mrSges(6,1) - t53 * mrSges(5,1);
t360 = qJD(2) * t291;
t299 = t286 * t360 + t287 * t357;
t345 = t286 * t363;
t457 = -t271 + (-t345 + t359) * pkin(3);
t288 = sin(qJ(1));
t292 = cos(qJ(1));
t475 = g(1) * t292 + g(2) * t288;
t474 = -t149 * t285 + t289 * t330;
t92 = t149 * t289 + t285 * t330;
t444 = m(5) * pkin(3);
t22 = qJD(5) * t474 + t285 * t79 + t289 * t81;
t443 = t22 / 0.2e1;
t23 = -qJD(5) * t92 - t285 * t81 + t289 * t79;
t442 = t23 / 0.2e1;
t430 = t79 / 0.2e1;
t429 = t81 / 0.2e1;
t473 = m(5) + m(6);
t424 = t140 / 0.2e1;
t423 = t141 / 0.2e1;
t217 = qJDD(5) + t225;
t418 = t217 / 0.2e1;
t417 = t225 / 0.2e1;
t246 = t284 * t286;
t248 = t284 * t290;
t159 = t283 * t246 + t248 * t282;
t128 = -pkin(8) * t227 + t159;
t160 = t282 * t246 - t283 * t248;
t129 = -pkin(8) * t311 + t160;
t62 = t128 * t289 - t129 * t285;
t470 = qJD(5) * t62 + t285 * t484 + t289 * t483;
t63 = t128 * t285 + t129 * t289;
t469 = -qJD(5) * t63 - t285 * t483 + t289 * t484;
t402 = pkin(3) * t283;
t264 = pkin(4) + t402;
t403 = pkin(3) * t282;
t196 = t264 * t285 + t289 * t403;
t57 = -t116 * t282 - t379;
t51 = t57 - t471;
t58 = t283 * t116 - t110;
t52 = t58 - t476;
t468 = -t196 * qJD(5) + t285 * t52 - t289 * t51;
t195 = t264 * t289 - t285 * t403;
t467 = t195 * qJD(5) - t285 * t51 - t289 * t52;
t466 = t444 + mrSges(4,1);
t461 = -pkin(4) * t489 + t457;
t222 = Ifges(4,4) * t233;
t133 = t234 * Ifges(4,1) + t258 * Ifges(4,5) + t222;
t269 = Ifges(3,4) * t363;
t460 = Ifges(3,1) * t364 + Ifges(3,5) * qJD(2) + t290 * t133 + t269;
t459 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t233 + mrSges(4,2) * t234 + mrSges(3,3) * t364;
t260 = pkin(6) * t371;
t179 = t286 * t243 + t260;
t458 = Ifges(4,5) * t140 + Ifges(5,5) * t81 + Ifges(4,6) * t141 + Ifges(5,6) * t79 + t225 * t485;
t224 = t239 * pkin(6);
t456 = t223 * t291 + t224 * t287;
t455 = -t286 * t68 + t290 * t67;
t251 = qJD(5) + t258;
t454 = t234 * Ifges(4,5) + t149 * Ifges(5,5) + t92 * Ifges(6,5) + t233 * Ifges(4,6) + Ifges(5,6) * t330 + t474 * Ifges(6,6) + t251 * Ifges(6,3) + t258 * t485;
t453 = -m(4) - m(3) - t473;
t270 = pkin(6) * t364;
t249 = -qJD(2) * pkin(2) + t270;
t161 = -pkin(3) * t233 + qJD(4) + t249;
t98 = -pkin(4) * t330 + t161;
t452 = -mrSges(6,1) * t98 + mrSges(6,3) * t12;
t451 = mrSges(6,2) * t98 - mrSges(6,3) * t11;
t401 = pkin(3) * t286;
t241 = pkin(4) * t273 + t401;
t449 = -m(6) * t241 + mrSges(2,2) - mrSges(3,3);
t326 = t291 * mrSges(3,1) - mrSges(3,2) * t287;
t448 = t287 * t487 + mrSges(2,1) + t326;
t446 = Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t418;
t445 = Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t418;
t441 = Ifges(5,4) * t429 + Ifges(5,2) * t430 + Ifges(5,6) * t417;
t440 = Ifges(5,1) * t429 + Ifges(5,4) * t430 + Ifges(5,5) * t417;
t405 = Ifges(6,4) * t92;
t35 = Ifges(6,2) * t474 + Ifges(6,6) * t251 + t405;
t439 = -t35 / 0.2e1;
t438 = t35 / 0.2e1;
t84 = Ifges(6,4) * t474;
t36 = Ifges(6,1) * t92 + Ifges(6,5) * t251 + t84;
t437 = -t36 / 0.2e1;
t436 = t36 / 0.2e1;
t435 = Ifges(4,1) * t424 + Ifges(4,4) * t423 + Ifges(4,5) * t417;
t77 = Ifges(5,4) * t149 + Ifges(5,2) * t330 + Ifges(5,6) * t258;
t434 = -t77 / 0.2e1;
t433 = t77 / 0.2e1;
t78 = Ifges(5,1) * t149 + Ifges(5,4) * t330 + Ifges(5,5) * t258;
t432 = -t78 / 0.2e1;
t431 = t78 / 0.2e1;
t428 = -t474 / 0.2e1;
t427 = t474 / 0.2e1;
t426 = -t92 / 0.2e1;
t425 = t92 / 0.2e1;
t422 = -t330 / 0.2e1;
t421 = t330 / 0.2e1;
t420 = -t149 / 0.2e1;
t419 = t149 / 0.2e1;
t415 = t234 / 0.2e1;
t414 = -t251 / 0.2e1;
t413 = t251 / 0.2e1;
t412 = -t258 / 0.2e1;
t411 = t258 / 0.2e1;
t409 = mrSges(5,3) * t53;
t408 = mrSges(5,3) * t54;
t404 = pkin(3) * t234;
t397 = g(3) * t287;
t276 = t287 * pkin(6);
t237 = t328 * qJD(2);
t362 = qJD(2) * t287;
t349 = pkin(6) * t362;
t366 = t290 * t237 + t286 * t349;
t83 = -t287 * t356 + t309 * qJD(2) + (-t260 + (qJ(4) * t287 - t243) * t286) * qJD(3) + t366;
t367 = t286 * t237 + t243 * t357;
t93 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t375 + (-qJD(4) * t287 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t291) * t286 + t367;
t41 = t282 * t83 + t283 * t93;
t393 = Ifges(3,4) * t291;
t392 = Ifges(4,4) * t286;
t391 = Ifges(4,4) * t290;
t390 = t154 * mrSges(4,3);
t389 = t155 * mrSges(4,3);
t388 = t234 * Ifges(4,4);
t378 = t286 * t287;
t376 = t286 * t292;
t374 = t288 * t286;
t373 = t288 * t291;
t370 = t291 * t292;
t229 = t290 * t243;
t151 = -qJ(4) * t375 + t229 + (-pkin(6) * t286 - pkin(3)) * t291;
t158 = -qJ(4) * t378 + t179;
t96 = t282 * t151 + t283 * t158;
t172 = t262 * t373 + t263 * t292;
t173 = t262 * t292 - t263 * t373;
t369 = -t172 * mrSges(6,1) + t173 * mrSges(6,2);
t174 = -t262 * t370 + t288 * t263;
t175 = t288 * t262 + t263 * t370;
t368 = t174 * mrSges(6,1) - t175 * mrSges(6,2);
t240 = pkin(3) * t378 + t276;
t358 = qJD(3) * t287;
t352 = Ifges(6,5) * t22 + Ifges(6,6) * t23 + Ifges(6,3) * t217;
t272 = pkin(6) * t360;
t171 = pkin(3) * t299 + t272;
t132 = t233 * Ifges(4,2) + t258 * Ifges(4,6) + t388;
t341 = -t286 * t132 / 0.2e1;
t37 = -t79 * mrSges(5,1) + t81 * mrSges(5,2);
t8 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t40 = -t282 * t93 + t283 * t83;
t95 = t283 * t151 - t158 * t282;
t200 = -qJDD(2) * pkin(2) + t224;
t323 = mrSges(4,1) * t286 + mrSges(4,2) * t290;
t322 = -mrSges(6,1) * t262 - mrSges(6,2) * t263;
t321 = Ifges(4,1) * t290 - t392;
t320 = Ifges(4,1) * t286 + t391;
t318 = -Ifges(4,2) * t286 + t391;
t317 = Ifges(4,2) * t290 + t392;
t316 = Ifges(3,5) * t291 - Ifges(3,6) * t287;
t315 = Ifges(4,5) * t290 - Ifges(4,6) * t286;
t314 = Ifges(4,5) * t286 + Ifges(4,6) * t290;
t193 = t311 * t287;
t61 = -pkin(4) * t291 + pkin(8) * t193 + t95;
t192 = t227 * t287;
t66 = -pkin(8) * t192 + t96;
t30 = -t285 * t66 + t289 * t61;
t31 = t285 * t61 + t289 * t66;
t120 = -t192 * t289 + t193 * t285;
t121 = -t192 * t285 - t193 * t289;
t142 = -t227 * t285 - t289 * t311;
t143 = t227 * t289 - t285 * t311;
t313 = t235 * t291 - t280 * t287;
t312 = t266 * t291 - t284 * t287;
t310 = t352 - t479;
t209 = -t286 * t370 + t288 * t290;
t207 = t286 * t373 + t290 * t292;
t307 = t249 * t323;
t300 = -t286 * t358 + t290 * t360;
t103 = -pkin(3) * t141 + qJDD(4) + t200;
t298 = Ifges(4,5) * t287 + t291 * t321;
t297 = Ifges(4,6) * t287 + t291 * t318;
t296 = Ifges(4,3) * t287 + t291 * t315;
t245 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t219 = t323 * t287;
t210 = t290 * t370 + t374;
t208 = -t288 * t371 + t376;
t190 = t288 * t273 + t274 * t370;
t189 = -t273 * t370 + t288 * t274;
t188 = t273 * t292 - t274 * t373;
t187 = t273 * t373 + t274 * t292;
t180 = pkin(4) * t311 - t266;
t178 = -pkin(6) * t377 + t229;
t170 = mrSges(4,1) * t258 - mrSges(4,3) * t234;
t169 = -mrSges(4,2) * t258 + mrSges(4,3) * t233;
t163 = -pkin(6) * t344 + t211;
t152 = pkin(4) * t192 + t240;
t123 = -qJD(2) * t304 - t203 * t287;
t122 = -qJD(2) * t305 + t311 * t358;
t119 = mrSges(5,1) * t258 - mrSges(5,3) * t149;
t118 = -mrSges(5,2) * t258 + mrSges(5,3) * t330;
t115 = -qJD(3) * t179 + t366;
t114 = (-t287 * t361 - t291 * t359) * pkin(6) + t367;
t109 = -t176 * t285 - t177 * t289;
t108 = -t176 * t289 + t177 * t285;
t107 = pkin(4) * t149 + t404;
t106 = -mrSges(4,2) * t225 + mrSges(4,3) * t141;
t105 = mrSges(4,1) * t225 - mrSges(4,3) * t140;
t97 = -pkin(4) * t122 + t171;
t94 = -mrSges(5,1) * t330 + mrSges(5,2) * t149;
t87 = -mrSges(4,1) * t141 + mrSges(4,2) * t140;
t70 = mrSges(6,1) * t251 - mrSges(6,3) * t92;
t69 = -mrSges(6,2) * t251 + mrSges(6,3) * t474;
t64 = t140 * Ifges(4,4) + t141 * Ifges(4,2) + t225 * Ifges(4,6);
t60 = mrSges(5,1) * t225 - mrSges(5,3) * t81;
t59 = -mrSges(5,2) * t225 + mrSges(5,3) * t79;
t50 = -pkin(4) * t79 + t103;
t49 = -qJD(5) * t121 + t122 * t289 - t123 * t285;
t48 = qJD(5) * t120 + t122 * t285 + t123 * t289;
t43 = -mrSges(6,1) * t474 + mrSges(6,2) * t92;
t29 = pkin(8) * t122 + t41;
t28 = pkin(4) * t362 - pkin(8) * t123 + t40;
t18 = -mrSges(6,2) * t217 + mrSges(6,3) * t23;
t17 = mrSges(6,1) * t217 - mrSges(6,3) * t22;
t5 = -qJD(5) * t31 + t28 * t289 - t285 * t29;
t4 = qJD(5) * t30 + t28 * t285 + t289 * t29;
t1 = [(t454 / 0.2e1 - t155 * mrSges(4,2) + t154 * mrSges(4,1) + Ifges(5,3) * t411 + Ifges(6,3) * t413 + Ifges(5,5) * t419 + Ifges(5,6) * t421 + Ifges(6,5) * t425 + Ifges(6,6) * t427 - t477) * t362 + (t276 * mrSges(3,3) + t287 * Ifges(3,1) + t393 / 0.2e1 - pkin(1) * mrSges(3,2)) * t239 + (Ifges(6,5) * t48 + Ifges(6,6) * t49) * t413 + (Ifges(6,5) * t121 + Ifges(6,6) * t120) * t418 + (Ifges(5,5) * t123 + Ifges(5,6) * t122 + qJD(2) * t296 - t314 * t358) * t411 + t456 * mrSges(3,3) + (Ifges(6,4) * t48 + Ifges(6,2) * t49) * t427 + (Ifges(6,4) * t121 + Ifges(6,2) * t120) * t442 + (-Ifges(5,4) * t193 - Ifges(5,2) * t192) * t430 + (Ifges(5,4) * t123 + Ifges(5,2) * t122) * t421 + (Ifges(6,1) * t48 + Ifges(6,4) * t49) * t425 + (Ifges(6,1) * t121 + Ifges(6,4) * t120) * t443 + (-Ifges(5,1) * t193 - Ifges(5,4) * t192) * t429 + (Ifges(5,1) * t123 + Ifges(5,4) * t122) * t419 + pkin(1) * mrSges(3,1) * t238 + t326 * t382 + t287 * t318 * t423 + t287 * t321 * t424 + (-mrSges(3,1) * t276 + Ifges(3,5) * t287) * qJDD(2) + (-Ifges(5,5) * t193 - Ifges(5,6) * t192 + t287 * t315) * t417 + m(5) * (t103 * t240 + t13 * t95 + t14 * t96 + t161 * t171 + t40 * t53 + t41 * t54) + m(6) * (t11 * t5 + t12 * t4 + t152 * t50 + t2 * t31 + t3 * t30 + t97 * t98) + t341 * t360 + t30 * t17 + t31 * t18 - t245 * t349 + t249 * (mrSges(4,1) * t299 + mrSges(4,2) * t300) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t456) + t459 * t272 + t460 * t360 / 0.2e1 + (t122 * t54 - t123 * t53 + t13 * t193 - t14 * t192) * mrSges(5,3) - t64 * t378 / 0.2e1 + (-t11 * t48 + t12 * t49 + t120 * t2 - t121 * t3) * mrSges(6,3) + (-t374 * t444 - t210 * mrSges(4,1) - t190 * mrSges(5,1) - t175 * mrSges(6,1) - t209 * mrSges(4,2) - t189 * mrSges(5,2) - t174 * mrSges(6,2) + t453 * (t292 * pkin(1) + t288 * pkin(6)) + t449 * t288 + (-m(4) * t329 - m(5) * t312 - m(6) * t313 - t448) * t292) * g(2) + m(4) * (t114 * t155 + t115 * t154 + t178 * t68 + t179 * t67 + (t200 * t287 + t249 * t360) * pkin(6)) + (-t154 * t300 - t155 * t299 - t375 * t68 - t378 * t67) * mrSges(4,3) + t233 * (qJD(2) * t297 - t317 * t358) / 0.2e1 + (t394 + t319) * t238 / 0.2e1 - (t290 * t132 + t286 * t133) * t358 / 0.2e1 + t103 * (mrSges(5,1) * t192 - mrSges(5,2) * t193) + (-t376 * t444 - t208 * mrSges(4,1) - t188 * mrSges(5,1) - t173 * mrSges(6,1) - t207 * mrSges(4,2) - t187 * mrSges(5,2) - t172 * mrSges(6,2) + (-m(5) * (-pkin(1) - t312) - m(4) * t243 - m(6) * (-pkin(1) - t313) + m(3) * pkin(1) + t448) * t288 + (t453 * pkin(6) + t449) * t292) * g(1) - t481 * t355 + t4 * t69 + t5 * t70 + t95 * t60 + t96 * t59 + t97 * t43 + t98 * (-mrSges(6,1) * t49 + mrSges(6,2) * t48) + ((-mrSges(3,2) * pkin(6) + Ifges(3,6)) * qJDD(2) + (-Ifges(3,2) * t287 + t393) * t355 / 0.2e1 - t485 * t417 - Ifges(5,5) * t429 - Ifges(5,6) * t430 - Ifges(6,3) * t418 - Ifges(4,6) * t423 - Ifges(4,5) * t424 - Ifges(6,6) * t442 - Ifges(6,5) * t443 - t352 / 0.2e1 - t458 / 0.2e1 + (pkin(6) * mrSges(3,3) + Ifges(3,2) / 0.2e1) * t238 + Ifges(3,4) * t239 / 0.2e1 + t478 + t479) * t291 + t41 * t118 + t40 * t119 + t50 * (-mrSges(6,1) * t120 + mrSges(6,2) * t121) + qJD(2) ^ 2 * t316 / 0.2e1 + t152 * t8 + t161 * (-mrSges(5,1) * t122 + mrSges(5,2) * t123) + t114 * t169 + t115 * t170 + t171 * t94 + t178 * t105 + t179 * t106 + t87 * t276 + (qJD(2) * t298 - t320 * t358) * t415 + t123 * t431 + t122 * t433 + t375 * t435 + t48 * t436 + t49 * t438 - t193 * t440 - t192 * t441 + t121 * t445 + t120 * t446 + t200 * t219 + t240 * t37 + Ifges(2,3) * qJDD(1); (Ifges(5,5) * t420 + Ifges(6,5) * t426 + Ifges(5,6) * t422 + Ifges(6,6) * t428 + Ifges(5,3) * t412 + Ifges(6,3) * t414 + t477) * t364 + (-mrSges(5,1) * t489 + t488 * mrSges(5,2)) * t161 + (Ifges(5,5) * t227 - Ifges(5,6) * t311 + t314) * t417 + (Ifges(5,1) * t227 - Ifges(5,4) * t311) * t429 + (Ifges(5,4) * t227 - Ifges(5,2) * t311) * t430 + t103 * (mrSges(5,1) * t311 + mrSges(5,2) * t227) - t311 * t441 + (-t13 * t227 - t14 * t311 + t176 * t54 - t177 * t53) * mrSges(5,3) + (-Ifges(5,4) * t177 - Ifges(5,2) * t176) * t422 + (-Ifges(5,5) * t177 - Ifges(5,6) * t176) * t412 + (-Ifges(5,1) * t177 - Ifges(5,4) * t176) * t420 + (-t390 + t133 / 0.2e1) * t357 + (Ifges(6,1) * t109 + Ifges(6,4) * t108) * t426 + (-t154 * (mrSges(4,1) * t287 - mrSges(4,3) * t371) - t155 * (-mrSges(4,2) * t287 - mrSges(4,3) * t377)) * qJD(1) + (Ifges(6,5) * t109 + Ifges(6,6) * t108) * t414 - t359 * t389 + t132 * t345 / 0.2e1 + (m(4) * ((-t154 * t290 - t155 * t286) * qJD(3) + t455) - t286 * t105 + t290 * t106 - t170 * t357 - t169 * t359) * pkin(7) + t455 * mrSges(4,3) + t457 * t94 - t459 * t271 - (-Ifges(3,2) * t364 + t269 + t460) * t363 / 0.2e1 + t461 * t43 + t469 * t70 + t470 * t69 + (t11 * t469 + t12 * t470 + t180 * t50 + t2 * t63 + t3 * t62 + t461 * t98) * m(6) + t464 * t118 + (-t103 * t266 + t13 * t159 + t14 * t160 + t161 * t457 + t464 * t54 + t465 * t53) * m(5) + t465 * t119 - (Ifges(5,4) * t419 + Ifges(5,2) * t421 + Ifges(5,6) * t411 + t408 + t433) * t203 - t454 * t364 / 0.2e1 - (Ifges(5,1) * t419 + Ifges(5,4) * t421 + Ifges(5,5) * t411 - t409 + t431) * t204 + (Ifges(6,4) * t109 + Ifges(6,2) * t108) * t428 + (Ifges(6,1) * t425 + Ifges(6,4) * t427 + Ifges(6,5) * t413 + t436 + t451) * (qJD(5) * t142 - t203 * t285 - t204 * t289) + (Ifges(6,4) * t425 + Ifges(6,2) * t427 + Ifges(6,6) * t413 + t438 + t452) * (-qJD(5) * t143 - t203 * t289 + t204 * t285) + t475 * (t287 * t486 + t325) + (-g(3) * t486 + t475 * t480) * t291 + (t307 + t341) * qJD(3) + (-t108 * t12 + t109 * t11 + t142 * t2 - t143 * t3) * mrSges(6,3) - t307 * t363 - t316 * t355 / 0.2e1 + t245 * t270 + (t233 * t318 + t234 * t321 + t258 * t315) * qJD(3) / 0.2e1 - (t233 * t297 + t234 * t298 + t258 * t296) * qJD(1) / 0.2e1 + (t287 * t480 - t326) * g(3) + t481 * qJD(1) ^ 2 + t62 * t17 + t63 * t18 - pkin(2) * t87 - t98 * (-mrSges(6,1) * t108 + mrSges(6,2) * t109) + t200 * t324 + t50 * (-mrSges(6,1) * t142 + mrSges(6,2) * t143) + t159 * t60 + t160 * t59 - t163 * t169 - t162 * t170 + t180 * t8 + (Ifges(6,5) * t143 + Ifges(6,6) * t142) * t418 + t317 * t423 + t320 * t424 - t177 * t432 - t176 * t434 + t286 * t435 + t109 * t437 + t108 * t439 + t227 * t440 + (Ifges(6,4) * t143 + Ifges(6,2) * t142) * t442 + (Ifges(6,1) * t143 + Ifges(6,4) * t142) * t443 + t143 * t445 + t142 * t446 + Ifges(3,3) * qJDD(2) - t223 * mrSges(3,2) - t224 * mrSges(3,1) + Ifges(3,6) * t238 + Ifges(3,5) * t239 - t266 * t37 + t290 * t64 / 0.2e1 + (-pkin(2) * t200 - t154 * t162 - t155 * t163 - t249 * t271) * m(4); -(Ifges(6,4) * t426 + Ifges(6,2) * t428 + Ifges(6,6) * t414 + t439 - t452) * t92 - (mrSges(5,1) * t161 + Ifges(5,4) * t420 + Ifges(5,2) * t422 + Ifges(5,6) * t412 - t408 + t434) * t149 + t310 - t94 * t404 - m(5) * (t161 * t404 + t53 * t57 + t54 * t58) - t234 * (Ifges(4,1) * t233 - t388) / 0.2e1 + (t187 * mrSges(5,1) - t188 * mrSges(5,2) - m(6) * (-t241 * t373 - t242 * t292) - t369 - mrSges(4,2) * t208 + t466 * t207) * g(2) + (-t189 * mrSges(5,1) + t190 * mrSges(5,2) - m(6) * (-t241 * t370 + t288 * t242) - t368 + mrSges(4,2) * t210 - t466 * t209) * g(1) + t467 * t69 + (-t107 * t98 + t11 * t468 + t12 * t467 + t195 * t3 + t196 * t2 + t241 * t397) * m(6) + t468 * t70 + (m(5) * t401 + mrSges(5,1) * t273 + mrSges(5,2) * t274 - t322) * t397 + t458 + (Ifges(6,1) * t426 + Ifges(6,4) * t428 + Ifges(6,5) * t414 + t437 - t451) * t474 - (-Ifges(4,2) * t234 + t133 + t222) * t233 / 0.2e1 + (-mrSges(5,2) * t161 + Ifges(5,1) * t420 + Ifges(5,4) * t422 + Ifges(5,5) * t412 + t409 + t432) * t330 - t478 - t107 * t43 - t58 * t118 - t57 * t119 - t154 * t169 + t155 * t170 + t234 * t389 + t233 * t390 + t60 * t402 + t59 * t403 + (Ifges(4,5) * t233 - Ifges(4,6) * t234) * t412 + t132 * t415 + (t13 * t283 + t14 * t282) * t444 + t195 * t17 + t196 * t18 + g(3) * t219 - t249 * (mrSges(4,1) * t234 + mrSges(4,2) * t233); -t330 * t118 + t149 * t119 - t474 * t69 + t92 * t70 + t37 + t8 + (t291 * g(3) - t287 * t475) * t473 + (t11 * t92 - t12 * t474 + t50) * m(6) + (t149 * t53 - t330 * t54 + t103) * m(5); -t98 * (mrSges(6,1) * t92 + mrSges(6,2) * t474) + (Ifges(6,1) * t474 - t405) * t426 + t35 * t425 + (Ifges(6,5) * t474 - Ifges(6,6) * t92) * t414 - t11 * t69 + t12 * t70 - g(1) * t368 - g(2) * t369 - t322 * t397 + (t11 * t474 + t12 * t92) * mrSges(6,3) + t310 + (-Ifges(6,2) * t92 + t36 + t84) * t428;];
tau = t1;
