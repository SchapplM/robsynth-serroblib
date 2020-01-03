% Calculate vector of inverse dynamics joint torques for
% S5RRPRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:48
% EndTime: 2019-12-31 20:36:34
% DurationCPUTime: 24.87s
% Computational Cost: add. (11489->737), mult. (28564->1023), div. (0->0), fcn. (23021->14), ass. (0->340)
t272 = sin(pkin(10));
t274 = cos(pkin(10));
t278 = sin(qJ(4));
t282 = cos(qJ(4));
t232 = t272 * t282 + t274 * t278;
t273 = sin(pkin(5));
t283 = cos(qJ(2));
t371 = t273 * t283;
t287 = t232 * t371;
t173 = qJD(1) * t287;
t222 = t232 * qJD(4);
t486 = t173 - t222;
t231 = t272 * t278 - t282 * t274;
t286 = t231 * t371;
t174 = qJD(1) * t286;
t221 = t231 * qJD(4);
t485 = -t174 + t221;
t277 = sin(qJ(5));
t281 = cos(qJ(5));
t317 = -mrSges(6,1) * t281 + mrSges(6,2) * t277;
t466 = m(6) * pkin(4) + mrSges(5,1) - t317;
t330 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t279 = sin(qJ(2));
t310 = pkin(2) * t279 - qJ(3) * t283;
t362 = qJD(1) * t273;
t212 = t310 * t362;
t339 = t279 * t362;
t275 = cos(pkin(5));
t361 = qJD(1) * t275;
t350 = pkin(1) * t361;
t213 = -pkin(7) * t339 + t283 * t350;
t141 = t274 * t212 - t272 * t213;
t369 = t274 * t283;
t288 = (pkin(3) * t279 - pkin(8) * t369) * t273;
t107 = qJD(1) * t288 + t141;
t142 = t272 * t212 + t274 * t213;
t338 = t283 * t362;
t325 = t272 * t338;
t122 = -pkin(8) * t325 + t142;
t401 = pkin(8) + qJ(3);
t239 = t401 * t272;
t240 = t401 * t274;
t299 = -t282 * t239 - t240 * t278;
t464 = -qJD(3) * t231 + qJD(4) * t299 - t278 * t107 - t282 * t122;
t484 = -pkin(9) * t339 + t464;
t255 = pkin(7) * t338;
t214 = t279 * t350 + t255;
t166 = pkin(3) * t325 + t214;
t483 = -t486 * pkin(4) + t485 * pkin(9) - t166;
t246 = qJD(4) - t338;
t413 = t246 / 0.2e1;
t260 = qJD(2) + t361;
t191 = t260 * t274 - t272 * t339;
t192 = t260 * t272 + t274 * t339;
t301 = t191 * t278 + t282 * t192;
t420 = t301 / 0.2e1;
t331 = t282 * t191 - t192 * t278;
t422 = t331 / 0.2e1;
t129 = qJD(5) - t331;
t424 = t129 / 0.2e1;
t101 = t246 * t277 + t281 * t301;
t426 = t101 / 0.2e1;
t100 = t246 * t281 - t277 * t301;
t428 = t100 / 0.2e1;
t177 = qJ(3) * t260 + t214;
t200 = (-pkin(2) * t283 - qJ(3) * t279 - pkin(1)) * t273;
t182 = qJD(1) * t200;
t109 = -t272 * t177 + t274 * t182;
t81 = -pkin(3) * t338 - t192 * pkin(8) + t109;
t110 = t274 * t177 + t272 * t182;
t86 = pkin(8) * t191 + t110;
t41 = t278 * t81 + t282 * t86;
t39 = pkin(9) * t246 + t41;
t168 = -t260 * pkin(2) + qJD(3) - t213;
t130 = -t191 * pkin(3) + t168;
t47 = -pkin(4) * t331 - pkin(9) * t301 + t130;
t13 = -t277 * t39 + t281 * t47;
t14 = t277 * t47 + t281 * t39;
t448 = -t130 * mrSges(5,1) - t13 * mrSges(6,1) + t14 * mrSges(6,2);
t35 = t101 * Ifges(6,5) + t100 * Ifges(6,6) + t129 * Ifges(6,3);
t393 = Ifges(5,4) * t301;
t76 = Ifges(5,2) * t331 + t246 * Ifges(5,6) + t393;
t479 = -t35 / 0.2e1 + t76 / 0.2e1;
t482 = -Ifges(5,4) * t420 + Ifges(6,5) * t426 - Ifges(5,2) * t422 - Ifges(5,6) * t413 + Ifges(6,6) * t428 + Ifges(6,3) * t424 - t448 - t479;
t359 = qJD(2) * t273;
t337 = t279 * t359;
t477 = -qJD(1) * t337 + qJDD(1) * t371;
t205 = qJDD(4) - t477;
t417 = t205 / 0.2e1;
t358 = qJD(2) * t283;
t218 = (qJD(1) * t358 + qJDD(1) * t279) * t273;
t352 = qJDD(1) * t275;
t259 = qJDD(2) + t352;
t155 = -t218 * t272 + t259 * t274;
t156 = t218 * t274 + t259 * t272;
t67 = -qJD(4) * t301 + t155 * t282 - t156 * t278;
t432 = t67 / 0.2e1;
t66 = qJD(4) * t331 + t155 * t278 + t156 * t282;
t433 = t66 / 0.2e1;
t65 = qJDD(5) - t67;
t434 = t65 / 0.2e1;
t34 = -qJD(5) * t101 + t205 * t281 - t277 * t66;
t439 = t34 / 0.2e1;
t33 = qJD(5) * t100 + t205 * t277 + t281 * t66;
t440 = t33 / 0.2e1;
t373 = t273 * t279;
t261 = pkin(7) * t373;
t408 = pkin(1) * t275;
t349 = qJD(2) * t408;
t327 = qJD(1) * t349;
t347 = pkin(1) * t352;
t148 = -qJD(2) * t255 - qJDD(1) * t261 - t279 * t327 + t283 * t347;
t135 = -t259 * pkin(2) + qJDD(3) - t148;
t85 = -t155 * pkin(3) + t135;
t17 = -t67 * pkin(4) - t66 * pkin(9) + t85;
t355 = qJD(4) * t282;
t356 = qJD(4) * t278;
t147 = pkin(7) * t477 + t279 * t347 + t283 * t327;
t117 = qJ(3) * t259 + qJD(3) * t260 + t147;
t357 = qJD(3) * t279;
t125 = -pkin(2) * t477 - qJ(3) * t218 + (-pkin(1) * qJDD(1) - qJD(1) * t357) * t273;
t73 = -t117 * t272 + t274 * t125;
t43 = -pkin(3) * t477 - pkin(8) * t156 + t73;
t74 = t274 * t117 + t272 * t125;
t49 = pkin(8) * t155 + t74;
t10 = t278 * t43 + t282 * t49 + t81 * t355 - t356 * t86;
t8 = pkin(9) * t205 + t10;
t1 = qJD(5) * t13 + t17 * t277 + t281 * t8;
t2 = -qJD(5) * t14 + t17 * t281 - t277 * t8;
t450 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t5 = Ifges(6,5) * t33 + Ifges(6,6) * t34 + Ifges(6,3) * t65;
t481 = t450 + mrSges(5,1) * t85 + Ifges(6,5) * t440 + Ifges(6,6) * t439 + Ifges(6,3) * t434 + t5 / 0.2e1 + (-t417 - t205 / 0.2e1) * Ifges(5,6) + (-t432 - t67 / 0.2e1) * Ifges(5,2) + (-t433 - t66 / 0.2e1) * Ifges(5,4);
t467 = t130 * mrSges(5,2);
t128 = Ifges(5,4) * t331;
t77 = Ifges(5,1) * t301 + t246 * Ifges(5,5) + t128;
t480 = t467 + Ifges(5,1) * t420 + Ifges(5,4) * t422 + Ifges(5,5) * t413 + t77 / 0.2e1;
t478 = t147 * mrSges(3,2);
t172 = -t239 * t278 + t240 * t282;
t463 = -qJD(3) * t232 - qJD(4) * t172 - t107 * t282 + t122 * t278;
t271 = pkin(10) + qJ(4);
t268 = sin(t271);
t269 = cos(t271);
t320 = -mrSges(4,1) * t274 + mrSges(4,2) * t272;
t290 = m(4) * pkin(2) - t320;
t476 = t330 * t268 - t269 * t466 - t290;
t11 = -qJD(4) * t41 - t278 * t49 + t282 * t43;
t475 = -t11 * mrSges(5,1) + t10 * mrSges(5,2);
t473 = t85 * mrSges(5,2) + 0.2e1 * Ifges(5,1) * t433 + 0.2e1 * Ifges(5,4) * t432 + 0.2e1 * Ifges(5,5) * t417;
t470 = -m(6) - m(5);
t267 = pkin(3) * t274 + pkin(2);
t151 = pkin(4) * t231 - pkin(9) * t232 - t267;
t93 = t151 * t281 - t172 * t277;
t469 = qJD(5) * t93 + t483 * t277 + t281 * t484;
t94 = t151 * t277 + t172 * t281;
t468 = -qJD(5) * t94 - t277 * t484 + t483 * t281;
t465 = -m(4) * qJ(3) - mrSges(4,3) - mrSges(5,3);
t462 = pkin(4) * t339 - t463;
t398 = mrSges(5,3) * t301;
t106 = mrSges(5,1) * t246 - t398;
t56 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t461 = -t56 + t106;
t329 = mrSges(3,3) * t339;
t460 = -mrSges(3,1) * t260 - mrSges(4,1) * t191 + mrSges(4,2) * t192 + t329;
t143 = t174 * t277 + t281 * t339;
t353 = qJD(5) * t281;
t295 = -t277 * t221 + t232 * t353;
t459 = t143 + t295;
t144 = -t174 * t281 + t277 * t339;
t354 = qJD(5) * t277;
t294 = t281 * t221 + t232 * t354;
t458 = t144 + t294;
t456 = -t272 * t73 + t274 * t74;
t455 = t1 * t281 - t2 * t277;
t397 = Ifges(3,4) * t279;
t454 = pkin(1) * (mrSges(3,1) * t279 + mrSges(3,2) * t283) - t279 * (Ifges(3,1) * t283 - t397) / 0.2e1;
t316 = mrSges(6,1) * t277 + mrSges(6,2) * t281;
t453 = -t316 + t465;
t452 = mrSges(3,2) + t453;
t179 = (qJD(2) * t310 - t357) * t273;
t215 = -pkin(7) * t337 + t283 * t349;
t186 = qJD(3) * t275 + t215;
t119 = t272 * t179 + t274 * t186;
t324 = t272 * t273 * t358;
t102 = -pkin(8) * t324 + t119;
t229 = pkin(7) * t371 + t279 * t408;
t199 = qJ(3) * t275 + t229;
t138 = t274 * t199 + t272 * t200;
t219 = -t272 * t373 + t274 * t275;
t104 = pkin(8) * t219 + t138;
t137 = -t272 * t199 + t274 * t200;
t220 = t272 * t275 + t274 * t373;
t91 = -pkin(3) * t371 - t220 * pkin(8) + t137;
t400 = t282 * t104 + t278 * t91;
t118 = t274 * t179 - t272 * t186;
t89 = qJD(2) * t288 + t118;
t21 = -qJD(4) * t400 - t102 * t278 + t282 * t89;
t451 = mrSges(3,1) - t476;
t388 = Ifges(3,6) * t260;
t40 = -t278 * t86 + t282 * t81;
t449 = t41 * mrSges(5,2) + t388 / 0.2e1 + (t283 * Ifges(3,2) + t397) * t362 / 0.2e1 - t40 * mrSges(5,1);
t414 = -t246 / 0.2e1;
t421 = -t301 / 0.2e1;
t447 = Ifges(5,1) * t421 + Ifges(5,5) * t414 - t467;
t423 = -t331 / 0.2e1;
t425 = -t129 / 0.2e1;
t427 = -t101 / 0.2e1;
t429 = -t100 / 0.2e1;
t446 = Ifges(6,5) * t427 - Ifges(5,2) * t423 - Ifges(5,6) * t414 + Ifges(6,6) * t429 + Ifges(6,3) * t425 + t448;
t445 = t273 ^ 2;
t443 = Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * t434;
t392 = Ifges(6,4) * t101;
t36 = Ifges(6,2) * t100 + Ifges(6,6) * t129 + t392;
t436 = t36 / 0.2e1;
t96 = Ifges(6,4) * t100;
t37 = Ifges(6,1) * t101 + Ifges(6,5) * t129 + t96;
t435 = -t37 / 0.2e1;
t419 = t155 / 0.2e1;
t418 = t156 / 0.2e1;
t416 = t219 / 0.2e1;
t415 = t220 / 0.2e1;
t412 = t275 / 0.2e1;
t411 = t281 / 0.2e1;
t410 = t283 / 0.2e1;
t409 = cos(qJ(1));
t407 = pkin(1) * t283;
t399 = mrSges(5,3) * t331;
t396 = Ifges(3,4) * t283;
t395 = Ifges(4,4) * t272;
t394 = Ifges(4,4) * t274;
t391 = Ifges(6,4) * t277;
t390 = Ifges(6,4) * t281;
t389 = Ifges(4,5) * t192;
t387 = Ifges(4,6) * t191;
t386 = Ifges(4,3) * t279;
t385 = t331 * Ifges(5,6);
t384 = t301 * Ifges(5,5);
t383 = t246 * Ifges(5,3);
t382 = t260 * Ifges(3,5);
t379 = t331 * t277;
t378 = t331 * t281;
t376 = t232 * t277;
t375 = t232 * t281;
t374 = t272 * t283;
t280 = sin(qJ(1));
t372 = t273 * t280;
t370 = t274 * (Ifges(4,1) * t192 + Ifges(4,4) * t191 - Ifges(4,5) * t338);
t367 = t279 * t280;
t366 = t280 * t283;
t363 = t409 * pkin(1) + pkin(7) * t372;
t351 = Ifges(5,5) * t66 + Ifges(5,6) * t67 + Ifges(5,3) * t205;
t346 = t277 * t371;
t345 = t281 * t371;
t344 = t37 * t411;
t343 = Ifges(3,5) * t218 + Ifges(3,6) * t477 + Ifges(3,3) * t259;
t342 = t273 * t409;
t341 = t409 * t279;
t340 = t409 * t283;
t28 = -t67 * mrSges(5,1) + t66 * mrSges(5,2);
t333 = -t354 / 0.2e1;
t332 = -pkin(1) * t280 + pkin(7) * t342;
t95 = -t155 * mrSges(4,1) + t156 * mrSges(4,2);
t224 = t275 * t341 + t366;
t158 = t224 * t269 - t268 * t342;
t157 = -t224 * t268 - t269 * t342;
t328 = mrSges(3,3) * t338;
t326 = t272 * t342;
t319 = mrSges(4,1) * t272 + mrSges(4,2) * t274;
t315 = Ifges(4,1) * t274 - t395;
t314 = Ifges(6,1) * t281 - t391;
t313 = -Ifges(4,2) * t272 + t394;
t312 = -Ifges(6,2) * t277 + t390;
t311 = Ifges(6,5) * t281 - Ifges(6,6) * t277;
t51 = -pkin(9) * t371 + t400;
t146 = t219 * t278 + t220 * t282;
t202 = t261 + (-pkin(2) - t407) * t275;
t149 = -t219 * pkin(3) + t202;
t300 = t282 * t219 - t220 * t278;
t72 = -pkin(4) * t300 - t146 * pkin(9) + t149;
t25 = t277 * t72 + t281 * t51;
t24 = -t277 * t51 + t281 * t72;
t70 = -mrSges(6,2) * t129 + mrSges(6,3) * t100;
t71 = mrSges(6,1) * t129 - mrSges(6,3) * t101;
t308 = -t277 * t71 + t281 * t70;
t306 = mrSges(3,2) + t465;
t225 = t275 * t366 + t341;
t226 = -t275 * t367 + t340;
t305 = t272 * pkin(3) * t372 + t225 * t401 + t226 * t267 + t363;
t54 = -t278 * t104 + t282 * t91;
t126 = -t277 * t146 - t345;
t297 = -t281 * t146 + t346;
t38 = -pkin(4) * t246 - t40;
t296 = t38 * t316;
t20 = t282 * t102 - t104 * t356 + t278 * t89 + t91 * t355;
t292 = t168 * t319;
t223 = -t275 * t340 + t367;
t285 = -g(1) * t225 - g(2) * t223 + g(3) * t371;
t216 = t229 * qJD(2);
t167 = pkin(3) * t324 + t216;
t253 = Ifges(3,4) * t338;
t228 = t275 * t407 - t261;
t227 = (-mrSges(3,1) * t283 + mrSges(3,2) * t279) * t273;
t211 = -t260 * mrSges(3,2) + t328;
t198 = t268 * t275 + t269 * t373;
t176 = Ifges(3,1) * t339 + t253 + t382;
t162 = t226 * t269 + t268 * t372;
t161 = t226 * t268 - t269 * t372;
t154 = -mrSges(4,1) * t338 - t192 * mrSges(4,3);
t153 = mrSges(4,2) * t338 + t191 * mrSges(4,3);
t124 = -mrSges(4,1) * t477 - mrSges(4,3) * t156;
t123 = mrSges(4,2) * t477 + mrSges(4,3) * t155;
t121 = t162 * t281 + t225 * t277;
t120 = -t162 * t277 + t225 * t281;
t115 = Ifges(4,4) * t192 + Ifges(4,2) * t191 - Ifges(4,6) * t338;
t114 = -Ifges(4,3) * t338 + t387 + t389;
t105 = -mrSges(5,2) * t246 + t399;
t98 = qJD(2) * t287 + qJD(4) * t146;
t97 = -qJD(2) * t286 + qJD(4) * t300;
t83 = t156 * Ifges(4,1) + t155 * Ifges(4,4) - Ifges(4,5) * t477;
t82 = t156 * Ifges(4,4) + t155 * Ifges(4,2) - Ifges(4,6) * t477;
t79 = pkin(4) * t301 - pkin(9) * t331;
t78 = -mrSges(5,1) * t331 + mrSges(5,2) * t301;
t75 = t383 + t384 + t385;
t58 = qJD(5) * t297 - t277 * t97 + t281 * t337;
t57 = qJD(5) * t126 + t277 * t337 + t281 * t97;
t53 = -mrSges(5,2) * t205 + mrSges(5,3) * t67;
t52 = mrSges(5,1) * t205 - mrSges(5,3) * t66;
t50 = pkin(4) * t371 - t54;
t44 = t98 * pkin(4) - t97 * pkin(9) + t167;
t23 = t277 * t79 + t281 * t40;
t22 = -t277 * t40 + t281 * t79;
t19 = -pkin(4) * t337 - t21;
t18 = pkin(9) * t337 + t20;
t16 = -mrSges(6,2) * t65 + mrSges(6,3) * t34;
t15 = mrSges(6,1) * t65 - mrSges(6,3) * t33;
t12 = -mrSges(6,1) * t34 + mrSges(6,2) * t33;
t9 = -pkin(4) * t205 - t11;
t6 = t33 * Ifges(6,4) + t34 * Ifges(6,2) + t65 * Ifges(6,6);
t4 = -qJD(5) * t25 - t18 * t277 + t281 * t44;
t3 = qJD(5) * t24 + t18 * t281 + t277 * t44;
t7 = [t460 * t216 + (-pkin(1) * t227 * t273 + Ifges(2,3)) * qJDD(1) + (Ifges(6,1) * t57 + Ifges(6,4) * t58) * t426 + (-m(5) * t305 - t162 * mrSges(5,1) - m(6) * (pkin(4) * t162 + t305) - t121 * mrSges(6,1) - t120 * mrSges(6,2) - t409 * mrSges(2,1) + (-mrSges(3,1) - t290) * t226 + (mrSges(2,2) + (-mrSges(3,3) - t319) * t273) * t280 + t330 * t161 + t306 * t225 + (-m(3) - m(4)) * t363) * g(2) + (-Ifges(4,6) * t155 - t156 * Ifges(4,5) - t351 / 0.2e1 - t73 * mrSges(4,1) + t147 * mrSges(3,3) + t74 * mrSges(4,2) - Ifges(5,6) * t432 - Ifges(5,5) * t433 - Ifges(5,3) * t417 + t475) * t371 + t58 * t436 + t343 * t412 + t83 * t415 + t82 * t416 + (Ifges(4,1) * t220 + Ifges(4,4) * t219) * t418 + (Ifges(4,4) * t220 + Ifges(4,2) * t219) * t419 + (t1 * t126 - t13 * t57 + t14 * t58 + t2 * t297) * mrSges(6,3) - t297 * t443 + (-Ifges(6,5) * t297 + Ifges(6,6) * t126) * t434 + (-Ifges(6,4) * t297 + Ifges(6,2) * t126) * t439 + t9 * (-mrSges(6,1) * t126 - mrSges(6,2) * t297) + (-Ifges(6,1) * t297 + Ifges(6,4) * t126) * t440 + (-m(3) * t332 + t224 * mrSges(3,1) - mrSges(3,3) * t342 - m(4) * (-pkin(2) * t224 + t332) - (-t224 * t274 + t326) * mrSges(4,1) - (t224 * t272 + t274 * t342) * mrSges(4,2) + t280 * mrSges(2,1) + t409 * mrSges(2,2) + t330 * t157 + t466 * t158 + (-t306 + t316) * t223 - t470 * (-pkin(3) * t326 + t223 * t401 + t224 * t267 - t332)) * g(1) + (-t40 * mrSges(5,3) + t480) * t97 + (mrSges(5,3) * t10 - t481) * t300 + (-mrSges(5,3) * t41 + t482) * t98 + (Ifges(6,5) * t57 + Ifges(6,6) * t58) * t424 + m(6) * (t1 * t25 + t13 * t4 + t14 * t3 + t19 * t38 + t2 * t24 + t50 * t9) + m(4) * (t109 * t118 + t110 * t119 + t135 * t202 + t137 * t73 + t138 * t74 + t168 * t216) + t148 * (mrSges(3,1) * t275 - mrSges(3,3) * t373) - t73 * t220 * mrSges(4,3) + m(5) * (t10 * t400 + t11 * t54 + t130 * t167 + t149 * t85 + t20 * t41 + t21 * t40) + t400 * t53 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t445 + t147 * t229 + t148 * t228 - t213 * t216 + t214 * t215) + (-mrSges(5,3) * t11 + t473) * t146 + t74 * t219 * mrSges(4,3) + (-t229 * mrSges(3,2) + t228 * mrSges(3,1) + Ifges(3,3) * t412 + (Ifges(3,5) * t279 + Ifges(3,6) * t283) * t273) * t259 + (-t228 * mrSges(3,3) + Ifges(3,5) * t412 + (-pkin(1) * mrSges(3,2) + t279 * Ifges(3,1) + t396) * t273) * t218 + (Ifges(6,4) * t57 + Ifges(6,2) * t58) * t428 + t135 * (-mrSges(4,1) * t219 + mrSges(4,2) * t220) + t215 * t211 + t202 * t95 + t167 * t78 + t119 * t153 + t118 * t154 + t149 * t28 + t138 * t123 + t137 * t124 + t126 * t6 / 0.2e1 + t20 * t105 + t21 * t106 + t3 * t70 + t4 * t71 + t57 * t37 / 0.2e1 + t38 * (-mrSges(6,1) * t58 + mrSges(6,2) * t57) + t54 * t52 + t19 * t56 + t50 * t12 - t275 * t478 + ((t370 / 0.2e1 - t272 * t115 / 0.2e1 - t213 * mrSges(3,3) + t292 + t382 / 0.2e1 + t191 * t313 / 0.2e1 + t192 * t315 / 0.2e1 + t176 / 0.2e1 + (-t109 * t274 - t110 * t272) * mrSges(4,3)) * t283 + (-t214 * mrSges(3,3) + t383 / 0.2e1 + t384 / 0.2e1 + t385 / 0.2e1 - t388 / 0.2e1 - t110 * mrSges(4,2) + t387 / 0.2e1 + t389 / 0.2e1 + t109 * mrSges(4,1) + t114 / 0.2e1 + t75 / 0.2e1 - t449) * t279 + ((-Ifges(3,2) * t279 + t396) * t410 - t283 * (Ifges(4,5) * t369 - Ifges(4,6) * t374 + t386) / 0.2e1 - t454) * t362) * t359 - (-t229 * mrSges(3,3) + Ifges(4,5) * t415 + Ifges(4,6) * t416 - Ifges(3,6) * t275 / 0.2e1 + (-pkin(1) * mrSges(3,1) - t397 + (-Ifges(3,2) - Ifges(4,3)) * t283) * t273) * t477 + t24 * t15 + t25 * t16; -t459 * t36 / 0.2e1 + (mrSges(6,1) * t459 - mrSges(6,2) * t458) * t38 + (-t1 * t376 + t13 * t458 - t14 * t459 - t2 * t375) * mrSges(6,3) + (-m(4) * t168 + t329 - t460) * t214 + t375 * t443 + (t227 + t470 * t267 * t371 + (t476 * t283 + (t401 * t470 + t453) * t279) * t273) * g(3) + (t153 * t274 - t154 * t272) * qJD(3) + t144 * t435 + (Ifges(4,1) * t272 + t394) * t418 + (Ifges(4,2) * t274 + t395) * t419 - (t192 * (Ifges(4,5) * t279 + t283 * t315) + t191 * (Ifges(4,6) * t279 + t283 * t313) + t260 * (Ifges(3,5) * t283 - Ifges(3,6) * t279) + (-Ifges(3,2) * t339 + t176 + t253 + t370) * t283 + (t114 + t75) * t279) * t362 / 0.2e1 + (t1 * t94 + t13 * t468 + t14 * t469 + t2 * t93 - t299 * t9 + t38 * t462) * m(6) + (t10 * t172 + t11 * t299 - t130 * t166 - t267 * t85 + t40 * t463 + t41 * t464) * m(5) - (t12 - t52) * t299 + (t470 * (-t225 * t267 + t226 * t401) + t452 * t226 + t451 * t225) * g(1) + (t470 * (-t223 * t267 + t224 * t401) + t452 * t224 + t451 * t223) * g(2) + (t123 * t274 - t124 * t272) * qJ(3) + (-Ifges(5,4) * t421 + t446 + t479) * t173 - (t344 + t480) * t221 + t481 * t231 + t482 * t222 + (-t109 * (mrSges(4,1) * t279 - mrSges(4,3) * t369) - t110 * (-mrSges(4,2) * t279 - mrSges(4,3) * t374)) * t362 + t468 * t71 + t469 * t70 + (Ifges(5,5) * t421 + Ifges(5,6) * t423 + Ifges(5,3) * t414 + t449) * t339 - (Ifges(5,4) * t423 - t77 / 0.2e1 + t447) * t174 + (-t10 * t231 - t11 * t232 + t485 * t40 + t486 * t41) * mrSges(5,3) + t463 * t106 + t464 * t105 + t462 * t56 - t478 + (-Ifges(6,5) * t294 - Ifges(6,6) * t295) * t424 + (Ifges(6,5) * t144 + Ifges(6,6) * t143) * t425 - t6 * t376 / 0.2e1 + (-Ifges(6,4) * t294 - Ifges(6,2) * t295) * t428 + (Ifges(6,4) * t144 + Ifges(6,2) * t143) * t429 - t267 * t28 + t272 * t83 / 0.2e1 + (-Ifges(6,1) * t294 - Ifges(6,4) * t295) * t426 + (Ifges(6,1) * t144 + Ifges(6,4) * t143) * t427 + (t328 - t211) * t213 + ((t386 + (Ifges(4,5) * t274 - Ifges(4,6) * t272) * t283) * t410 + t454) * qJD(1) ^ 2 * t445 + (t311 * t434 + t312 * t439 + t314 * t440 + t316 * t9 + t333 * t37 + t473) * t232 + t343 + t172 * t53 - t166 * t78 - t292 * t338 - t142 * t153 - t141 * t154 + t148 * mrSges(3,1) + t93 * t15 + t94 * t16 - pkin(2) * t95 + t135 * t320 + t274 * t82 / 0.2e1 + (-t109 * t141 - t110 * t142 - pkin(2) * t135 + (-t109 * t272 + t110 * t274) * qJD(3) + t456 * qJ(3)) * m(4) + t456 * mrSges(4,3) + t115 * t325 / 0.2e1 - t477 * (Ifges(4,5) * t272 + Ifges(4,6) * t274) / 0.2e1; t281 * t15 - t191 * t153 + t192 * t154 + t277 * t16 + t461 * t301 + t308 * qJD(5) + (-t105 - t308) * t331 + t95 + t28 + (t1 * t277 - t301 * t38 + t2 * t281 + t285 + t129 * (-t13 * t277 + t14 * t281)) * m(6) + (t301 * t40 - t331 * t41 + t285 + t85) * m(5) + (t109 * t192 - t110 * t191 + t135 + t285) * m(4); (t344 + t296) * qJD(5) + (Ifges(6,1) * t277 + t390) * t440 + t277 * t443 + (Ifges(6,5) * t277 + Ifges(6,6) * t281) * t434 + t378 * t435 + t379 * t436 + (Ifges(6,2) * t281 + t391) * t439 + t6 * t411 + t76 * t420 + (t399 - t105) * t40 + (-pkin(4) * t9 - t13 * t22 - t14 * t23) * m(6) + t446 * t301 + (t311 * t425 + t312 * t429 + t314 * t427 - t296 + t447) * t331 + (-t393 + t35) * t421 + (t161 * t466 + t162 * t330) * g(1) + (t330 * t198 - t466 * (-t268 * t373 + t269 * t275)) * g(3) + (-t157 * t466 + t158 * t330) * g(2) + (-m(6) * t38 + t398 + t461) * t41 + t351 - t475 + (t100 * t312 + t101 * t314 + t129 * t311) * qJD(5) / 0.2e1 + (t77 + t128) * t423 + t36 * t333 - t23 * t70 - t22 * t71 + t9 * t317 + ((-t354 + t379) * t14 + (-t353 + t378) * t13 + t455) * mrSges(6,3) + (-t71 * t353 - t70 * t354 + m(6) * ((-t13 * t281 - t14 * t277) * qJD(5) + t455) - t277 * t15 + t281 * t16) * pkin(9) - pkin(4) * t12; -t38 * (mrSges(6,1) * t101 + mrSges(6,2) * t100) + (Ifges(6,1) * t100 - t392) * t427 + t36 * t426 + (Ifges(6,5) * t100 - Ifges(6,6) * t101) * t425 - t13 * t70 + t14 * t71 - g(1) * (mrSges(6,1) * t120 - mrSges(6,2) * t121) - g(2) * ((-t158 * t277 + t223 * t281) * mrSges(6,1) + (-t158 * t281 - t223 * t277) * mrSges(6,2)) - g(3) * ((-t198 * t277 - t345) * mrSges(6,1) + (-t198 * t281 + t346) * mrSges(6,2)) + (t100 * t13 + t101 * t14) * mrSges(6,3) + t5 + (-Ifges(6,2) * t101 + t37 + t96) * t429 + t450;];
tau = t7;
