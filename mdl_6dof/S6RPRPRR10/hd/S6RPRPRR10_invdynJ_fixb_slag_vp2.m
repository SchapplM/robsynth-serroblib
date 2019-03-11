% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:45
% EndTime: 2019-03-09 04:08:22
% DurationCPUTime: 27.10s
% Computational Cost: add. (11416->745), mult. (23184->1000), div. (0->0), fcn. (16133->14), ass. (0->326)
t279 = sin(qJ(3));
t283 = cos(qJ(3));
t345 = qJD(1) * qJD(3);
t240 = qJDD(1) * t279 + t283 * t345;
t406 = t240 / 0.2e1;
t485 = t279 / 0.2e1;
t278 = sin(qJ(5));
t282 = cos(qJ(5));
t349 = qJD(5) * t282;
t350 = qJD(5) * t278;
t239 = qJDD(1) * t283 - t279 * t345;
t274 = sin(pkin(10));
t275 = cos(pkin(10));
t187 = qJDD(3) * t274 + t239 * t275;
t321 = -qJD(4) * t283 + qJD(2);
t344 = qJDD(1) * qJ(2);
t130 = pkin(3) * t240 - qJ(4) * t239 + qJD(1) * t321 + t344;
t285 = -pkin(1) - pkin(7);
t251 = qJDD(1) * t285 + qJDD(2);
t253 = qJD(1) * t285 + qJD(2);
t354 = qJD(3) * t283;
t168 = t279 * t251 + t253 * t354;
t158 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t168;
t76 = t275 * t130 - t158 * t274;
t54 = pkin(4) * t240 - pkin(8) * t187 + t76;
t186 = qJDD(3) * t275 - t239 * t274;
t77 = t274 * t130 + t275 * t158;
t57 = pkin(8) * t186 + t77;
t243 = pkin(3) * t279 - qJ(4) * t283 + qJ(2);
t214 = t243 * qJD(1);
t242 = t279 * t253;
t217 = qJD(3) * qJ(4) + t242;
t135 = t275 * t214 - t217 * t274;
t356 = qJD(1) * t283;
t226 = qJD(3) * t274 + t275 * t356;
t267 = t279 * qJD(1);
t98 = pkin(4) * t267 - pkin(8) * t226 + t135;
t136 = t274 * t214 + t275 * t217;
t225 = qJD(3) * t275 - t274 * t356;
t99 = pkin(8) * t225 + t136;
t11 = t278 * t54 + t282 * t57 + t98 * t349 - t350 * t99;
t50 = t278 * t98 + t282 * t99;
t12 = -qJD(5) * t50 - t278 * t57 + t282 * t54;
t489 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t151 = t225 * t278 + t226 * t282;
t72 = -qJD(5) * t151 + t186 * t282 - t187 * t278;
t10 = pkin(9) * t72 + t11;
t281 = cos(qJ(6));
t277 = sin(qJ(6));
t322 = t282 * t225 - t226 * t278;
t42 = pkin(9) * t322 + t50;
t377 = t277 * t42;
t258 = t267 + qJD(5);
t49 = -t278 * t99 + t282 * t98;
t41 = -pkin(9) * t151 + t49;
t40 = pkin(5) * t258 + t41;
t13 = t281 * t40 - t377;
t230 = qJDD(5) + t240;
t71 = qJD(5) * t322 + t186 * t278 + t187 * t282;
t9 = pkin(5) * t230 - pkin(9) * t71 + t12;
t2 = qJD(6) * t13 + t10 * t281 + t277 * t9;
t376 = t281 * t42;
t14 = t277 * t40 + t376;
t3 = -qJD(6) * t14 - t10 * t277 + t281 * t9;
t488 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t319 = mrSges(4,1) * t279 + mrSges(4,2) * t283;
t388 = pkin(8) + qJ(4);
t476 = m(6) * t388 - m(7) * (-pkin(9) - t388) + mrSges(6,3) + mrSges(7,3) + m(5) * qJ(4) + mrSges(5,3);
t487 = -t283 * t476 + t319;
t402 = t258 / 0.2e1;
t252 = qJD(6) + t258;
t404 = t252 / 0.2e1;
t411 = t151 / 0.2e1;
t413 = t322 / 0.2e1;
t84 = t151 * t281 + t277 * t322;
t417 = t84 / 0.2e1;
t464 = -t151 * t277 + t281 * t322;
t419 = t464 / 0.2e1;
t486 = Ifges(6,5) * t411 + Ifges(7,5) * t417 + Ifges(6,6) * t413 + Ifges(7,6) * t419 + Ifges(6,3) * t402 + Ifges(7,3) * t404;
t329 = t267 / 0.2e1;
t484 = t283 / 0.2e1;
t384 = Ifges(4,4) * t283;
t483 = Ifges(4,2) * t485 - t384 / 0.2e1;
t232 = t274 * t282 + t275 * t278;
t209 = t232 * qJD(1);
t184 = t279 * t209;
t211 = t232 * qJD(5);
t449 = t184 + t211;
t305 = t274 * t278 - t275 * t282;
t451 = t279 * t305;
t185 = qJD(1) * t451;
t442 = qJD(5) * t305;
t448 = t185 + t442;
t309 = pkin(3) * t283 + qJ(4) * t279;
t236 = t309 * qJD(1);
t366 = t274 * t283;
t162 = t275 * t236 - t253 * t366;
t365 = t275 * t279;
t342 = pkin(8) * t365;
t113 = (pkin(4) * t283 + t342) * qJD(1) + t162;
t364 = t275 * t283;
t163 = t274 * t236 + t253 * t364;
t338 = t274 * t267;
t132 = pkin(8) * t338 + t163;
t244 = t388 * t274;
t245 = t388 * t275;
t351 = qJD(4) * t275;
t352 = qJD(4) * t274;
t456 = -t244 * t349 + (-t132 + t351) * t282 + (-qJD(5) * t245 - t113 - t352) * t278;
t166 = -t278 * t244 + t282 * t245;
t455 = -t232 * qJD(4) - qJD(5) * t166 - t282 * t113 + t132 * t278;
t310 = Ifges(5,5) * t275 - Ifges(5,6) * t274;
t320 = mrSges(4,1) * t283 - mrSges(4,2) * t279;
t478 = qJ(2) * t320 + (-Ifges(4,1) * t279 - t384) * t484 + (Ifges(5,3) * t283 - t279 * t310) * t485;
t261 = pkin(4) * t275 + pkin(3);
t273 = pkin(10) + qJ(5);
t264 = cos(t273);
t318 = -mrSges(5,1) * t275 + mrSges(5,2) * t274;
t477 = -m(6) * t261 - m(7) * (pkin(5) * t264 + t261) - m(5) * pkin(3) + t318;
t475 = t13 * mrSges(7,1) + t49 * mrSges(6,1) - Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t483 + t226 * Ifges(5,5) / 0.2e1 + t225 * Ifges(5,6) / 0.2e1 + Ifges(5,3) * t329 - t14 * mrSges(7,2) - t50 * mrSges(6,2) + t486;
t474 = -pkin(5) * t356 + pkin(9) * t448 + t455;
t473 = pkin(9) * t449 - t456;
t198 = t232 * t283;
t447 = t305 * qJD(1) - qJD(3) * t198 + t279 * t442;
t200 = t305 * t283;
t446 = -qJD(3) * t200 - t211 * t279 - t209;
t280 = sin(qJ(1));
t284 = cos(qJ(1));
t472 = g(1) * t280 - g(2) * t284;
t471 = qJD(3) / 0.2e1;
t470 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t225 + mrSges(5,2) * t226 + mrSges(4,3) * t356;
t382 = Ifges(5,4) * t275;
t312 = -Ifges(5,2) * t274 + t382;
t383 = Ifges(5,4) * t274;
t314 = Ifges(5,1) * t275 - t383;
t469 = t225 * (Ifges(5,6) * t283 - t279 * t312) + t226 * (Ifges(5,5) * t283 - t279 * t314);
t346 = qJD(1) * qJD(2);
t254 = t344 + t346;
t355 = qJD(3) * t279;
t167 = t251 * t283 - t253 * t355;
t306 = t167 * t283 + t168 * t279;
t347 = m(5) + m(6) + m(7);
t466 = -m(4) - t347;
t317 = t274 * mrSges(5,1) + t275 * mrSges(5,2);
t367 = t253 * t283;
t206 = -qJD(3) * pkin(3) + qJD(4) - t367;
t369 = t206 * t279;
t465 = t317 * t369 - t136 * (mrSges(5,3) * t274 * t279 - mrSges(5,2) * t283) - t135 * (mrSges(5,1) * t283 + mrSges(5,3) * t365);
t409 = t187 / 0.2e1;
t463 = Ifges(5,1) * t409 + Ifges(5,5) * t406;
t462 = t279 * t477 + mrSges(2,2) - mrSges(3,3) - t487;
t263 = sin(t273);
t395 = pkin(5) * t263;
t398 = pkin(4) * t274;
t460 = m(6) * t398 + m(7) * (t395 + t398) + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t317;
t22 = qJD(6) * t464 + t277 * t72 + t281 * t71;
t432 = t22 / 0.2e1;
t23 = -qJD(6) * t84 - t277 * t71 + t281 * t72;
t431 = t23 / 0.2e1;
t424 = t71 / 0.2e1;
t423 = t72 / 0.2e1;
t410 = t186 / 0.2e1;
t219 = qJDD(6) + t230;
t408 = t219 / 0.2e1;
t407 = t230 / 0.2e1;
t165 = -t282 * t244 - t245 * t278;
t123 = -pkin(9) * t232 + t165;
t124 = -pkin(9) * t305 + t166;
t62 = t123 * t281 - t124 * t277;
t459 = qJD(6) * t62 + t277 * t474 - t473 * t281;
t63 = t123 * t277 + t124 * t281;
t458 = -qJD(6) * t63 + t473 * t277 + t281 * t474;
t433 = m(7) * pkin(5);
t457 = -mrSges(6,1) - t433;
t197 = t232 * t279;
t117 = -t197 * t277 - t281 * t451;
t454 = -qJD(6) * t117 - t277 * t446 + t281 * t447;
t115 = -t197 * t281 + t277 * t451;
t453 = qJD(6) * t115 + t277 * t447 + t281 * t446;
t223 = t275 * t243;
t328 = -t274 * t285 + pkin(4);
t147 = -pkin(8) * t364 + t279 * t328 + t223;
t361 = t279 * t285;
t174 = t274 * t243 + t275 * t361;
t161 = -pkin(8) * t366 + t174;
t88 = t278 * t147 + t282 * t161;
t140 = -mrSges(5,2) * t240 + mrSges(5,3) * t186;
t141 = mrSges(5,1) * t240 - mrSges(5,3) * t187;
t445 = t140 * t275 - t141 * t274;
t308 = -t274 * t76 + t275 * t77;
t189 = -pkin(4) * t338 + t242;
t443 = pkin(5) * t449 - t189;
t108 = -t186 * mrSges(5,1) + t187 * mrSges(5,2);
t34 = -t72 * mrSges(6,1) + t71 * mrSges(6,2);
t8 = -t23 * mrSges(7,1) + t22 * mrSges(7,2);
t441 = -t108 - t34 - t8;
t265 = qJ(6) + t273;
t259 = sin(t265);
t260 = cos(t265);
t440 = mrSges(6,1) * t264 + mrSges(7,1) * t260 - mrSges(6,2) * t263 - mrSges(7,2) * t259 - t477;
t159 = -pkin(4) * t225 + t206;
t90 = -pkin(5) * t322 + t159;
t438 = -mrSges(7,1) * t90 + t14 * mrSges(7,3);
t437 = mrSges(7,2) * t90 - mrSges(7,3) * t13;
t436 = qJD(1) ^ 2;
t435 = Ifges(7,4) * t432 + Ifges(7,2) * t431 + Ifges(7,6) * t408;
t434 = Ifges(7,1) * t432 + Ifges(7,4) * t431 + Ifges(7,5) * t408;
t430 = Ifges(6,4) * t424 + Ifges(6,2) * t423 + Ifges(6,6) * t407;
t429 = Ifges(6,1) * t424 + Ifges(6,4) * t423 + Ifges(6,5) * t407;
t399 = Ifges(7,4) * t84;
t38 = Ifges(7,2) * t464 + Ifges(7,6) * t252 + t399;
t428 = -t38 / 0.2e1;
t427 = t38 / 0.2e1;
t78 = Ifges(7,4) * t464;
t39 = Ifges(7,1) * t84 + Ifges(7,5) * t252 + t78;
t426 = -t39 / 0.2e1;
t425 = t39 / 0.2e1;
t381 = Ifges(6,4) * t151;
t74 = Ifges(6,2) * t322 + t258 * Ifges(6,6) + t381;
t422 = t74 / 0.2e1;
t146 = Ifges(6,4) * t322;
t75 = t151 * Ifges(6,1) + t258 * Ifges(6,5) + t146;
t421 = t75 / 0.2e1;
t420 = -t464 / 0.2e1;
t418 = -t84 / 0.2e1;
t416 = Ifges(5,4) * t410 + t463;
t415 = -m(3) - m(4);
t414 = -t322 / 0.2e1;
t412 = -t151 / 0.2e1;
t405 = -t252 / 0.2e1;
t403 = -t258 / 0.2e1;
t397 = pkin(5) * t151;
t392 = g(3) * t283;
t387 = mrSges(6,3) * t322;
t386 = mrSges(6,3) * t151;
t385 = Ifges(4,4) * t279;
t164 = -qJDD(3) * pkin(3) + qJDD(4) - t167;
t372 = t164 * t283;
t363 = t279 * t280;
t362 = t279 * t284;
t360 = t283 * t285;
t176 = -t259 * t363 + t260 * t284;
t177 = t259 * t284 + t260 * t363;
t359 = t176 * mrSges(7,1) - t177 * mrSges(7,2);
t178 = t259 * t362 + t260 * t280;
t179 = -t259 * t280 + t260 * t362;
t358 = t178 * mrSges(7,1) + t179 * mrSges(7,2);
t201 = qJD(3) * t309 + t321;
t353 = qJD(3) * t285;
t336 = t283 * t353;
t157 = t274 * t201 + t275 * t336;
t357 = t284 * pkin(1) + t280 * qJ(2);
t348 = qJDD(1) * mrSges(3,2);
t341 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t219;
t340 = Ifges(6,5) * t71 + Ifges(6,6) * t72 + Ifges(6,3) * t230;
t337 = t274 * t355;
t257 = t279 * t353;
t327 = -t345 / 0.2e1;
t324 = (t254 + t346) * qJ(2);
t87 = t282 * t147 - t161 * t278;
t224 = pkin(4) * t366 - t360;
t316 = -mrSges(7,1) * t259 - mrSges(7,2) * t260;
t315 = t283 * Ifges(4,1) - t385;
t311 = -Ifges(4,5) * t279 - Ifges(4,6) * t283;
t60 = pkin(5) * t279 + pkin(9) * t200 + t87;
t61 = -pkin(9) * t198 + t88;
t30 = -t277 * t61 + t281 * t60;
t31 = t277 * t60 + t281 * t61;
t307 = -t135 * t274 + t136 * t275;
t116 = -t198 * t281 + t200 * t277;
t118 = -t198 * t277 - t200 * t281;
t153 = -t232 * t277 - t281 * t305;
t154 = t232 * t281 - t277 * t305;
t207 = -pkin(4) * t337 + t257;
t304 = t341 + t488;
t192 = t263 * t362 + t264 * t280;
t190 = -t263 * t363 + t264 * t284;
t302 = t279 * (-Ifges(4,2) * t283 - t385);
t181 = t275 * t201;
t107 = t181 + (t283 * t328 + t342) * qJD(3);
t129 = pkin(8) * t337 + t157;
t35 = t278 * t107 + t282 * t129 + t147 * t349 - t161 * t350;
t100 = -pkin(4) * t186 + t164;
t36 = -qJD(5) * t88 + t282 * t107 - t129 * t278;
t262 = -pkin(1) * qJDD(1) + qJDD(2);
t247 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t267;
t237 = t319 * qJD(1);
t216 = Ifges(4,5) * qJD(3) + qJD(1) * t315;
t203 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t240;
t202 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t239;
t193 = -t263 * t280 + t264 * t362;
t191 = t263 * t284 + t264 * t363;
t188 = pkin(5) * t305 - t261;
t183 = mrSges(5,1) * t267 - mrSges(5,3) * t226;
t182 = -mrSges(5,2) * t267 + mrSges(5,3) * t225;
t173 = -t274 * t361 + t223;
t156 = -t274 * t336 + t181;
t155 = pkin(5) * t198 + t224;
t144 = t226 * Ifges(5,1) + t225 * Ifges(5,4) + Ifges(5,5) * t267;
t143 = t226 * Ifges(5,4) + t225 * Ifges(5,2) + Ifges(5,6) * t267;
t128 = t232 * t355 + t283 * t442;
t126 = qJD(3) * t451 - t211 * t283;
t112 = mrSges(6,1) * t258 - t386;
t111 = -mrSges(6,2) * t258 + t387;
t104 = t184 * t277 + t185 * t281;
t103 = t184 * t281 - t185 * t277;
t97 = -pkin(5) * t128 + t207;
t95 = t187 * Ifges(5,4) + t186 * Ifges(5,2) + t240 * Ifges(5,6);
t89 = -mrSges(6,1) * t322 + mrSges(6,2) * t151;
t67 = mrSges(7,1) * t252 - mrSges(7,3) * t84;
t66 = -mrSges(7,2) * t252 + mrSges(7,3) * t464;
t59 = -mrSges(6,2) * t230 + mrSges(6,3) * t72;
t58 = mrSges(6,1) * t230 - mrSges(6,3) * t71;
t48 = -qJD(6) * t118 - t126 * t277 + t128 * t281;
t46 = qJD(6) * t116 + t126 * t281 + t128 * t277;
t44 = -pkin(5) * t72 + t100;
t43 = -mrSges(7,1) * t464 + mrSges(7,2) * t84;
t29 = pkin(9) * t128 + t35;
t28 = pkin(5) * t354 - pkin(9) * t126 + t36;
t18 = -mrSges(7,2) * t219 + mrSges(7,3) * t23;
t17 = mrSges(7,1) * t219 - mrSges(7,3) * t22;
t16 = t281 * t41 - t377;
t15 = -t277 * t41 - t376;
t5 = -qJD(6) * t31 - t277 * t29 + t28 * t281;
t4 = qJD(6) * t30 + t277 * t28 + t281 * t29;
t1 = [(-t193 * mrSges(6,1) - t179 * mrSges(7,1) + t192 * mrSges(6,2) + t178 * mrSges(7,2) + (m(3) * pkin(1) + t285 * t466 + t460) * t280 + ((-m(3) + t466) * qJ(2) + t462) * t284) * g(1) + (t311 * t471 - t465) * qJD(3) + t143 * t337 / 0.2e1 + m(4) * (t285 * t306 + t324) + m(3) * (-pkin(1) * t262 + t324) + (Ifges(6,4) * t126 + Ifges(6,2) * t128) * t413 + (-Ifges(6,4) * t200 - Ifges(6,2) * t198) * t423 + (Ifges(7,3) * t408 + Ifges(5,5) * t409 + Ifges(5,6) * t410 + Ifges(6,3) * t407 + Ifges(6,6) * t423 + Ifges(6,5) * t424 + t76 * mrSges(5,1) - t77 * mrSges(5,2) - Ifges(4,4) * t239 / 0.2e1 - Ifges(4,6) * qJDD(3) + Ifges(7,6) * t431 + Ifges(7,5) * t432 + t488 + t489 + (Ifges(5,3) + Ifges(4,2)) * t406) * t279 + t240 * t483 + (Ifges(4,1) * t239 - Ifges(4,4) * t240) * t484 + (Ifges(5,5) * t187 + Ifges(5,6) * t186 + Ifges(5,3) * t240 + t340 + t341) * t485 + t469 * t471 + (t475 + t486) * t354 + t364 * t416 + t317 * t372 + t203 * t361 + t239 * t315 / 0.2e1 + (t319 + 0.2e1 * mrSges(3,3)) * t254 + t100 * (mrSges(6,1) * t198 - mrSges(6,2) * t200) + t302 * t327 + m(5) * (t135 * t156 + t136 * t157 + t173 * t76 + t174 * t77 + (t206 * t355 - t372) * t285) + (Ifges(4,5) * qJDD(3) + t310 * t406 + t312 * t410 + t314 * t409) * t283 + m(7) * (t13 * t5 + t14 * t4 + t155 * t44 + t2 * t31 + t3 * t30 + t90 * t97) + m(6) * (t100 * t224 + t11 * t88 + t12 * t87 + t159 * t207 + t35 * t50 + t36 * t49) + (Ifges(7,1) * t46 + Ifges(7,4) * t48) * t417 + (Ifges(7,1) * t118 + Ifges(7,4) * t116) * t432 - pkin(1) * t348 - t95 * t366 / 0.2e1 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) - (t275 * t144 + t216) * t355 / 0.2e1 + (Ifges(7,5) * t118 + Ifges(7,6) * t116) * t408 + (Ifges(7,5) * t46 + Ifges(7,6) * t48) * t404 + (-m(3) * t357 - t191 * mrSges(6,1) - t177 * mrSges(7,1) - t190 * mrSges(6,2) - t176 * mrSges(7,2) + t466 * (t284 * pkin(7) + t357) - t460 * t284 + t462 * t280) * g(2) + (t202 - t108) * t360 + (Ifges(6,1) * t126 + Ifges(6,4) * t128) * t411 + (-Ifges(6,1) * t200 - Ifges(6,4) * t198) * t424 + (Ifges(7,4) * t46 + Ifges(7,2) * t48) * t419 + (Ifges(7,4) * t118 + Ifges(7,2) * t116) * t431 - t306 * mrSges(4,3) + t470 * t257 + (Ifges(6,5) * t126 + Ifges(6,6) * t128) * t402 + (-Ifges(6,5) * t200 - Ifges(6,6) * t198) * t407 + (t116 * t2 - t118 * t3 - t13 * t46 + t14 * t48) * mrSges(7,3) + t478 * t345 + (-t11 * t198 + t12 * t200 - t126 * t49 + t128 * t50) * mrSges(6,3) + (-t364 * t76 - t366 * t77) * mrSges(5,3) + t262 * mrSges(3,2) + qJD(2) * t237 + qJ(2) * (mrSges(4,1) * t240 + mrSges(4,2) * t239) + t224 * t34 + t207 * t89 + t157 * t182 + t156 * t183 + t173 * t141 + t174 * t140 + t155 * t8 + t159 * (-mrSges(6,1) * t128 + mrSges(6,2) * t126) + t44 * (-mrSges(7,1) * t116 + mrSges(7,2) * t118) + t35 * t111 + t36 * t112 + t97 * t43 + t87 * t58 + t88 * t59 + t90 * (-mrSges(7,1) * t48 + mrSges(7,2) * t46) + t4 * t66 + t5 * t67 + t31 * t18 + t30 * t17 + t247 * t336 + t126 * t421 + t128 * t422 + t46 * t425 + t48 * t427 - t200 * t429 - t198 * t430 + t118 * t434 + t116 * t435; t348 + t115 * t17 + t117 * t18 - t197 * t58 - t451 * t59 + t454 * t67 + t453 * t66 + t447 * t112 + t446 * t111 + (qJ(2) * t415 - mrSges(3,3)) * t436 + (t203 + t445) * t279 + (-t182 * t274 - t183 * t275 - t237) * qJD(1) + (t202 + t441) * t283 + ((t182 * t275 - t183 * t274 + t247) * t283 + (t43 + t89 + t470) * t279) * qJD(3) + m(3) * t262 + m(4) * t306 - t472 * (t347 - t415) + (t115 * t3 + t117 * t2 + t13 * t454 + t14 * t453 - t283 * t44 + t90 * t355) * m(7) + (-t100 * t283 - t11 * t451 - t12 * t197 + t159 * t355 + t446 * t50 + t447 * t49) * m(6) + (-t372 + t308 * t279 + (t283 * t307 + t369) * qJD(3) - (t135 * t275 + t136 * t274) * qJD(1)) * m(5); -t143 * t338 / 0.2e1 + (t351 - t163) * t182 + (t465 - t469 / 0.2e1) * qJD(1) + (Ifges(7,5) * t154 + Ifges(7,6) * t153) * t408 + t164 * t318 + t311 * t327 + t216 * t329 + (-t103 * t14 + t104 * t13 + t153 * t2 - t154 * t3) * mrSges(7,3) + t383 * t410 + (t144 * t329 + Ifges(5,6) * t406 + Ifges(5,2) * t410 + t95 / 0.2e1) * t275 + t382 * t409 - t470 * t242 + (Ifges(6,5) * t185 + Ifges(6,6) * t184) * t403 + (Ifges(7,5) * t104 + Ifges(7,6) * t103) * t405 - t247 * t367 - t442 * t421 + (-Ifges(6,1) * t442 - Ifges(6,4) * t211) * t411 + (-Ifges(6,4) * t442 - Ifges(6,2) * t211) * t413 + (-Ifges(6,5) * t442 - Ifges(6,6) * t211) * t402 + (Ifges(7,1) * t417 + Ifges(7,4) * t419 + Ifges(7,5) * t404 + t425 + t437) * (qJD(6) * t153 - t211 * t277 - t281 * t442) + (Ifges(7,4) * t417 + Ifges(7,2) * t419 + Ifges(7,6) * t404 + t427 + t438) * (-qJD(6) * t154 - t211 * t281 + t277 * t442) + (Ifges(6,5) * t232 - Ifges(6,6) * t305) * t407 + t100 * (mrSges(6,1) * t305 + mrSges(6,2) * t232) + (Ifges(6,4) * t232 - Ifges(6,2) * t305) * t423 + (Ifges(6,1) * t232 - Ifges(6,4) * t305) * t424 + (-t11 * t305 - t12 * t232 + t448 * t49 - t449 * t50) * mrSges(6,3) - t305 * t430 + (-pkin(3) * t164 + qJ(4) * t308 + qJD(4) * t307 - t135 * t162 - t136 * t163 - t206 * t242) * m(5) + t472 * (-t279 * t476 - t283 * t440 - t320) + (t416 + t463) * t274 + (t279 * t440 + t487) * g(3) + (Ifges(6,1) * t185 + Ifges(6,4) * t184) * t412 + (-t352 - t162) * t183 + (t302 / 0.2e1 - t478) * t436 + (Ifges(6,5) * t412 + Ifges(7,5) * t418 + Ifges(6,6) * t414 + Ifges(7,6) * t420 + Ifges(6,3) * t403 + Ifges(7,3) * t405 - t475) * t356 - t261 * t34 + Ifges(4,5) * t239 - Ifges(4,6) * t240 + t188 * t8 - t189 * t89 - t184 * t74 / 0.2e1 - t185 * t75 / 0.2e1 + t165 * t58 + t166 * t59 + t167 * mrSges(4,1) - t168 * mrSges(4,2) + t44 * (-mrSges(7,1) * t153 + mrSges(7,2) * t154) - pkin(3) * t108 - t90 * (-mrSges(7,1) * t103 + mrSges(7,2) * t104) + Ifges(4,3) * qJDD(3) + t62 * t17 + t63 * t18 + (Ifges(7,4) * t104 + Ifges(7,2) * t103) * t420 + (Ifges(7,1) * t104 + Ifges(7,4) * t103) * t418 + t443 * t43 + t308 * mrSges(5,3) + t445 * qJ(4) + (mrSges(6,1) * t449 - mrSges(6,2) * t448) * t159 + t455 * t112 + t456 * t111 + (-t100 * t261 + t11 * t166 + t12 * t165 - t159 * t189 + t455 * t49 + t456 * t50) * m(6) + t458 * t67 + t459 * t66 + (t13 * t458 + t14 * t459 + t188 * t44 + t2 * t63 + t3 * t62 + t443 * t90) * m(7) + (Ifges(6,4) * t185 + Ifges(6,2) * t184) * t414 - t211 * t422 + t104 * t426 + t103 * t428 + t232 * t429 + (Ifges(7,4) * t154 + Ifges(7,2) * t153) * t431 + (Ifges(7,1) * t154 + Ifges(7,4) * t153) * t432 + t154 * t434 + t153 * t435; -t322 * t111 + t151 * t112 - t225 * t182 + t226 * t183 - t464 * t66 + t84 * t67 + (t13 * t84 - t14 * t464 + t44) * m(7) + (t151 * t49 - t322 * t50 + t100) * m(6) + (t135 * t226 - t136 * t225 + t164) * m(5) + (-t279 * g(3) + t283 * t472) * t347 - t441; t74 * t411 + (t387 - t111) * t49 + (Ifges(6,1) * t322 - t381) * t412 + (Ifges(6,5) * t322 - Ifges(6,6) * t151) * t403 - t159 * (mrSges(6,1) * t151 + mrSges(6,2) * t322) - t43 * t397 - m(7) * (t13 * t15 + t14 * t16 + t397 * t90) + (-Ifges(6,2) * t151 + t146 + t75) * t414 + (Ifges(7,1) * t418 + Ifges(7,4) * t420 + Ifges(7,5) * t405 + t426 - t437) * t464 + (m(7) * t395 + mrSges(6,1) * t263 + mrSges(6,2) * t264 - t316) * t392 + t489 + t340 + (t386 + t112) * t50 + t304 + ((-t277 * t67 + t281 * t66) * qJD(6) + t17 * t281 + t18 * t277) * pkin(5) - t16 * t66 - t15 * t67 - (Ifges(7,4) * t418 + Ifges(7,2) * t420 + Ifges(7,6) * t405 + t428 - t438) * t84 + (-mrSges(6,2) * t193 + t192 * t457 - t358) * g(2) + (mrSges(6,2) * t191 + t190 * t457 - t359) * g(1) + (t2 * t277 + t281 * t3 + (-t13 * t277 + t14 * t281) * qJD(6)) * t433; -t90 * (mrSges(7,1) * t84 + mrSges(7,2) * t464) + (Ifges(7,1) * t464 - t399) * t418 + t38 * t417 + (Ifges(7,5) * t464 - Ifges(7,6) * t84) * t405 - t13 * t66 + t14 * t67 - g(1) * t359 - g(2) * t358 - t316 * t392 + (t13 * t464 + t14 * t84) * mrSges(7,3) + t304 + (-Ifges(7,2) * t84 + t39 + t78) * t420;];
tau  = t1;
