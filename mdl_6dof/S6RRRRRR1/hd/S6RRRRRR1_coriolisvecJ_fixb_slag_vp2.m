% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:52
% EndTime: 2019-03-10 03:28:11
% DurationCPUTime: 11.22s
% Computational Cost: add. (32400->669), mult. (84439->919), div. (0->0), fcn. (64555->10), ass. (0->320)
t305 = sin(qJ(3));
t306 = sin(qJ(2));
t310 = cos(qJ(3));
t311 = cos(qJ(2));
t273 = -t305 * t306 + t310 * t311;
t263 = t273 * qJD(1);
t274 = t305 * t311 + t310 * t306;
t264 = t274 * qJD(1);
t304 = sin(qJ(4));
t309 = cos(qJ(4));
t220 = t263 * t304 + t264 * t309;
t303 = sin(qJ(5));
t308 = cos(qJ(5));
t349 = t309 * t263 - t264 * t304;
t469 = t308 * t220 + t303 * t349;
t474 = Ifges(6,4) * t469;
t470 = -t220 * t303 + t308 * t349;
t156 = Ifges(6,4) * t470;
t101 = pkin(5) * t469 - pkin(11) * t470;
t301 = qJD(2) + qJD(3);
t239 = t301 * t273;
t229 = t239 * qJD(1);
t240 = t301 * t274;
t230 = t240 * qJD(1);
t335 = -t229 * t309 + t230 * t304;
t125 = qJD(4) * t349 - t335;
t126 = -qJD(4) * t220 - t229 * t304 - t230 * t309;
t212 = Ifges(5,4) * t349;
t300 = qJD(4) + t301;
t153 = Ifges(5,1) * t220 + t300 * Ifges(5,5) + t212;
t296 = -pkin(2) * t311 - pkin(1);
t284 = qJD(1) * t296;
t241 = -t263 * pkin(3) + t284;
t291 = qJD(5) + t300;
t302 = sin(qJ(6));
t307 = cos(qJ(6));
t142 = t291 * t307 - t302 * t469;
t58 = qJD(5) * t470 + t125 * t308 + t126 * t303;
t39 = qJD(6) * t142 + t307 * t58;
t143 = t291 * t302 + t307 * t469;
t40 = -qJD(6) * t143 - t302 * t58;
t59 = qJD(5) * t469 + t125 * t303 - t308 * t126;
t12 = Ifges(7,4) * t39 + Ifges(7,2) * t40 + Ifges(7,6) * t59;
t13 = Ifges(7,1) * t39 + Ifges(7,4) * t40 + Ifges(7,5) * t59;
t339 = Ifges(7,5) * t302 + Ifges(7,6) * t307;
t403 = Ifges(7,4) * t302;
t341 = Ifges(7,2) * t307 + t403;
t402 = Ifges(7,4) * t307;
t343 = Ifges(7,1) * t302 + t402;
t346 = mrSges(7,1) * t307 - mrSges(7,2) * t302;
t372 = qJD(1) * t306;
t298 = pkin(2) * t372;
t197 = pkin(3) * t230 + qJD(2) * t298;
t102 = -pkin(4) * t126 + t197;
t17 = pkin(5) * t59 - pkin(11) * t58 + t102;
t256 = t264 * pkin(9);
t440 = -pkin(8) - pkin(7);
t454 = t440 * t311;
t279 = qJD(1) * t454;
t265 = t305 * t279;
t285 = t440 * t306;
t278 = qJD(1) * t285;
t271 = qJD(2) * pkin(2) + t278;
t451 = t310 * t271 + t265;
t191 = t451 - t256;
t180 = pkin(3) * t301 + t191;
t268 = t310 * t279;
t234 = t271 * t305 - t268;
t419 = pkin(9) * t263;
t192 = t234 + t419;
t181 = t304 * t192;
t127 = t309 * t180 - t181;
t215 = pkin(10) * t220;
t106 = t127 - t215;
t103 = pkin(4) * t300 + t106;
t183 = t309 * t192;
t128 = t180 * t304 + t183;
t418 = pkin(10) * t349;
t107 = t128 + t418;
t373 = t308 * t107;
t63 = t103 * t303 + t373;
t61 = pkin(11) * t291 + t63;
t171 = -pkin(4) * t349 + t241;
t85 = -pkin(5) * t470 - pkin(11) * t469 + t171;
t21 = -t302 * t61 + t307 * t85;
t358 = t310 * t440;
t359 = t305 * t440;
t360 = qJD(1) * qJD(2);
t317 = ((-t304 * t359 + t309 * t358) * t311 + (-t304 * t358 - t309 * t359) * t306) * t360;
t242 = t310 * t285 + t305 * t454;
t167 = t451 * qJD(3) + t242 * t360;
t414 = t230 * pkin(9);
t319 = t167 - t414;
t468 = -t305 * t285 + t310 * t454;
t168 = -t234 * qJD(3) + t360 * t468;
t415 = t229 * pkin(9);
t320 = t168 - t415;
t366 = qJD(4) * t309;
t367 = qJD(4) * t304;
t51 = t180 * t366 - t192 * t367 + t304 * t320 + t309 * t319;
t32 = pkin(10) * t126 + t51;
t325 = -t125 * pkin(10) - t180 * t367 - t192 * t366;
t368 = qJD(3) * t310;
t369 = qJD(3) * t305;
t383 = t107 * t303;
t62 = t103 * t308 - t383;
t8 = t62 * qJD(5) + t308 * t32 + (-t304 * (t271 * t368 + t279 * t369 - t414) + t309 * (-t271 * t369 + t279 * t368 - t415) + t325 + t317) * t303;
t2 = qJD(6) * t21 + t17 * t302 + t307 * t8;
t417 = t2 * t307;
t422 = t307 / 0.2e1;
t443 = t59 / 0.2e1;
t444 = t40 / 0.2e1;
t445 = t39 / 0.2e1;
t157 = qJD(6) - t470;
t345 = mrSges(7,1) * t302 + mrSges(7,2) * t307;
t60 = -pkin(5) * t291 - t62;
t330 = t60 * t345;
t340 = Ifges(7,5) * t307 - Ifges(7,6) * t302;
t342 = -Ifges(7,2) * t302 + t402;
t344 = Ifges(7,1) * t307 - t403;
t423 = -t302 / 0.2e1;
t434 = t143 / 0.2e1;
t395 = t143 * Ifges(7,4);
t75 = t142 * Ifges(7,2) + t157 * Ifges(7,6) + t395;
t141 = Ifges(7,4) * t142;
t76 = t143 * Ifges(7,1) + t157 * Ifges(7,5) + t141;
t467 = t142 * t342 / 0.2e1 + t344 * t434 + t157 * t340 / 0.2e1 + t330 + t75 * t423 + t76 * t422;
t9 = t303 * t32 - t308 * (-t304 * t319 + t309 * t320 + t325) + t63 * qJD(5);
t322 = -t8 * mrSges(6,2) + mrSges(7,3) * t417 + t302 * t13 / 0.2e1 + t12 * t422 + t343 * t445 + t341 * t444 + t339 * t443 - Ifges(6,6) * t59 + Ifges(6,5) * t58 + (-mrSges(6,1) - t346) * t9 + t467 * qJD(6);
t405 = Ifges(5,4) * t220;
t430 = -t220 / 0.2e1;
t52 = -t128 * qJD(4) + t335 * pkin(9) + (-t234 * t309 - t304 * t451) * qJD(3) + t317;
t473 = t52 * mrSges(5,1) - t51 * mrSges(5,2) + Ifges(5,5) * t125 + Ifges(5,6) * t126 + t322 - (Ifges(5,5) * t349 - Ifges(5,6) * t220) * t300 / 0.2e1 - (-Ifges(5,2) * t220 + t153 + t212) * t349 / 0.2e1 - t241 * (mrSges(5,1) * t220 + mrSges(5,2) * t349) + (Ifges(5,1) * t349 - t405) * t430;
t436 = -t142 / 0.2e1;
t435 = -t143 / 0.2e1;
t433 = -t157 / 0.2e1;
t472 = pkin(4) * t220;
t16 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t471 = m(7) * t9 + t16;
t382 = t128 * t220;
t152 = Ifges(5,2) * t349 + t300 * Ifges(5,6) + t405;
t465 = t152 / 0.2e1;
t22 = t302 * t85 + t307 * t61;
t392 = t21 * t307;
t338 = t22 * t302 + t392;
t461 = t338 * mrSges(7,3);
t295 = pkin(2) * t310 + pkin(3);
t376 = t304 * t305;
t258 = -pkin(2) * t376 + t309 * t295;
t253 = pkin(4) + t258;
t374 = t305 * t309;
t260 = pkin(2) * t374 + t295 * t304;
t204 = t253 * t308 - t260 * t303;
t226 = t295 * t366 + (-t305 * t367 + (t309 * t310 - t376) * qJD(3)) * pkin(2);
t227 = -t295 * t367 + (-t305 * t366 + (-t304 * t310 - t374) * qJD(3)) * pkin(2);
t116 = qJD(5) * t204 + t226 * t308 + t227 * t303;
t237 = -t278 * t305 + t268;
t193 = t237 - t419;
t238 = t310 * t278 + t265;
t194 = -t256 + t238;
t137 = t304 * t193 + t309 * t194;
t111 = -t215 + t137;
t136 = t309 * t193 - t194 * t304;
t331 = t136 - t418;
t69 = t308 * t111 + t303 * t331;
t460 = t116 - t69;
t205 = t303 * t253 + t308 * t260;
t459 = qJD(5) * t205 + (-t227 + t331) * t308 + (-t111 + t226) * t303;
t390 = mrSges(6,1) * t291 + mrSges(7,1) * t142 - mrSges(7,2) * t143 - mrSges(6,3) * t469;
t294 = pkin(3) * t309 + pkin(4);
t364 = qJD(5) * t308;
t365 = qJD(5) * t303;
t377 = t303 * t304;
t224 = t294 * t364 + (-t304 * t365 + (t308 * t309 - t377) * qJD(4)) * pkin(3);
t134 = t309 * t191 - t181;
t109 = -t215 + t134;
t133 = -t191 * t304 - t183;
t332 = t133 - t418;
t67 = t308 * t109 + t303 * t332;
t458 = t224 - t67;
t375 = t304 * t308;
t457 = -t109 * t303 + t308 * t332 + t294 * t365 + (t304 * t364 + (t303 * t309 + t375) * qJD(4)) * pkin(3);
t456 = t127 * t349;
t453 = -t136 + t227;
t452 = -t137 + t226;
t213 = -pkin(9) * t274 + t242;
t214 = pkin(9) * t273 - t468;
t155 = t304 * t213 + t309 * t214;
t450 = -t21 * t302 + t22 * t307;
t401 = Ifges(6,5) * t291;
t408 = Ifges(6,1) * t469;
t95 = t156 + t401 + t408;
t449 = t171 * mrSges(6,2) + t95 / 0.2e1 + t401 / 0.2e1 + t156 / 0.2e1 + t467;
t384 = qJD(6) * t22;
t3 = t17 * t307 - t302 * t8 - t384;
t448 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t39 + Ifges(7,6) * t40;
t447 = Ifges(6,6) * t291 / 0.2e1 + t474 / 0.2e1;
t442 = Ifges(7,5) * t435 + Ifges(7,6) * t436 + Ifges(7,3) * t433;
t399 = Ifges(6,2) * t470;
t441 = t399 / 0.2e1 + t447;
t439 = pkin(1) * mrSges(3,1);
t438 = pkin(1) * mrSges(3,2);
t154 = t309 * t213 - t214 * t304;
t236 = t273 * t304 + t274 * t309;
t119 = -pkin(10) * t236 + t154;
t235 = t273 * t309 - t274 * t304;
t120 = pkin(10) * t235 + t155;
t337 = t308 * t119 - t120 * t303;
t437 = t337 * t9;
t431 = t349 / 0.2e1;
t429 = t220 / 0.2e1;
t427 = t263 / 0.2e1;
t426 = -t264 / 0.2e1;
t425 = t264 / 0.2e1;
t421 = m(4) * t284;
t420 = pkin(3) * t264;
t416 = t21 * mrSges(7,3);
t413 = t3 * t302;
t412 = t62 * mrSges(6,3);
t411 = t63 * mrSges(6,3);
t410 = mrSges(4,3) * t263;
t409 = mrSges(4,3) * t264;
t407 = Ifges(3,4) * t306;
t406 = Ifges(4,4) * t264;
t394 = t469 * t63;
t388 = Ifges(3,5) * qJD(2);
t387 = Ifges(3,6) * qJD(2);
t386 = qJD(2) * mrSges(3,1);
t385 = qJD(2) * mrSges(3,2);
t381 = t470 * t302;
t380 = t470 * t307;
t259 = pkin(3) * t375 + t303 * t294;
t371 = qJD(1) * t311;
t370 = qJD(2) * t306;
t363 = qJD(6) * t302;
t362 = qJD(6) * t307;
t354 = qJD(2) * t440;
t353 = t388 / 0.2e1;
t352 = -t387 / 0.2e1;
t221 = pkin(2) * t370 + pkin(3) * t240;
t18 = mrSges(7,1) * t59 - mrSges(7,3) * t39;
t91 = -mrSges(7,2) * t157 + mrSges(7,3) * t142;
t351 = -qJD(6) * t91 - t18;
t179 = t420 + t472;
t347 = (-t3 - t384) * mrSges(7,3);
t84 = t119 * t303 + t120 * t308;
t170 = t235 * t303 + t236 * t308;
t248 = -t273 * pkin(3) + t296;
t184 = -t235 * pkin(4) + t248;
t334 = t308 * t235 - t236 * t303;
t89 = -pkin(5) * t334 - t170 * pkin(11) + t184;
t34 = t302 * t89 + t307 * t84;
t33 = -t302 * t84 + t307 * t89;
t139 = -qJD(4) * t236 - t239 * t304 - t240 * t309;
t112 = -pkin(4) * t139 + t221;
t257 = -pkin(3) * t377 + t294 * t308;
t280 = t306 * t354;
t281 = t311 * t354;
t175 = t310 * t280 + t305 * t281 + t285 * t368 + t369 * t454;
t149 = -pkin(9) * t240 + t175;
t176 = qJD(3) * t468 - t280 * t305 + t310 * t281;
t150 = -pkin(9) * t239 + t176;
t77 = t309 * t149 + t304 * t150 + t213 * t366 - t214 * t367;
t87 = t101 + t179;
t326 = -qJD(6) * t338 - t413;
t78 = -qJD(4) * t155 - t149 * t304 + t309 * t150;
t324 = m(7) * (t326 + t417);
t138 = qJD(4) * t235 + t239 * t309 - t240 * t304;
t323 = -pkin(10) * t138 + t78;
t318 = -t171 * mrSges(6,1) - t21 * mrSges(7,1) + t22 * mrSges(7,2) + t441 + 0.2e1 * t442 + t447;
t19 = -mrSges(7,2) * t59 + mrSges(7,3) * t40;
t92 = mrSges(7,1) * t157 - mrSges(7,3) * t143;
t316 = t307 * t19 + m(7) * (-t21 * t362 - t22 * t363 - t413 + t417) - t91 * t363 - t92 * t362 - t302 * t18;
t315 = -t399 / 0.2e1 - t318;
t313 = t408 / 0.2e1 + t449;
t210 = Ifges(4,2) * t263 + t301 * Ifges(4,6) + t406;
t255 = Ifges(4,4) * t263;
t211 = Ifges(4,1) * t264 + t301 * Ifges(4,5) + t255;
t312 = -(-Ifges(6,2) * t469 + t156 + t95) * t470 / 0.2e1 + t469 * t441 + t469 * t442 - t22 * (-mrSges(7,2) * t469 - mrSges(7,3) * t381) - t21 * (mrSges(7,1) * t469 - mrSges(7,3) * t380) + (Ifges(7,3) * t469 + t340 * t470) * t433 + (Ifges(7,5) * t469 + t344 * t470) * t435 + (Ifges(7,6) * t469 + t342 * t470) * t436 - t171 * (mrSges(6,1) * t469 + mrSges(6,2) * t470) - t291 * (Ifges(6,5) * t470 - Ifges(6,6) * t469) / 0.2e1 + t470 * t412 - t470 * t330 - t469 * (Ifges(6,1) * t470 - t474) / 0.2e1 + t473 + t168 * mrSges(4,1) - t167 * mrSges(4,2) + t210 * t425 + (Ifges(4,1) * t263 - t406) * t426 - (-Ifges(4,2) * t264 + t211 + t255) * t263 / 0.2e1 + t220 * t465 + t451 * t410 - t284 * (mrSges(4,1) * t264 + mrSges(4,2) * t263) + Ifges(4,5) * t229 - Ifges(4,6) * t230 - t301 * (Ifges(4,5) * t263 - Ifges(4,6) * t264) / 0.2e1 - t76 * t380 / 0.2e1 + t75 * t381 / 0.2e1 + mrSges(5,3) * t456;
t297 = Ifges(3,4) * t371;
t283 = mrSges(3,3) * t371 - t385;
t282 = -mrSges(3,3) * t372 + t386;
t262 = Ifges(3,1) * t372 + t297 + t388;
t261 = t387 + (t311 * Ifges(3,2) + t407) * qJD(1);
t254 = pkin(11) + t259;
t252 = -pkin(5) - t257;
t246 = mrSges(4,1) * t301 - t409;
t245 = -mrSges(4,2) * t301 + t410;
t244 = t298 + t420;
t232 = -mrSges(4,1) * t263 + mrSges(4,2) * t264;
t203 = pkin(11) + t205;
t202 = -pkin(5) - t204;
t196 = mrSges(5,1) * t300 - mrSges(5,3) * t220;
t195 = -mrSges(5,2) * t300 + mrSges(5,3) * t349;
t172 = t179 + t298;
t165 = -mrSges(5,1) * t349 + mrSges(5,2) * t220;
t144 = -mrSges(6,2) * t291 + mrSges(6,3) * t470;
t100 = -mrSges(6,1) * t470 + mrSges(6,2) * t469;
t88 = t101 + t472;
t86 = t298 + t87;
t73 = qJD(5) * t170 + t138 * t303 - t308 * t139;
t72 = qJD(5) * t334 + t138 * t308 + t139 * t303;
t65 = t106 * t308 - t383;
t64 = t106 * t303 + t373;
t55 = Ifges(7,3) * t59;
t43 = pkin(10) * t139 + t77;
t30 = t101 * t302 + t307 * t62;
t29 = t101 * t307 - t302 * t62;
t28 = t302 * t88 + t307 * t65;
t27 = -t302 * t65 + t307 * t88;
t26 = t302 * t86 + t307 * t69;
t25 = -t302 * t69 + t307 * t86;
t24 = t302 * t87 + t307 * t67;
t23 = -t302 * t67 + t307 * t87;
t20 = pkin(5) * t73 - pkin(11) * t72 + t112;
t15 = qJD(5) * t84 + t303 * t43 - t308 * t323;
t14 = qJD(5) * t337 + t303 * t323 + t308 * t43;
t5 = -qJD(6) * t34 - t14 * t302 + t20 * t307;
t4 = qJD(6) * t33 + t14 * t307 + t20 * t302;
t1 = [m(4) * (-t167 * t468 + t168 * t242 + t175 * t234 + t176 * t451) + (t167 * t273 - t168 * t274 - t229 * t242 + t230 * t468 - t234 * t240 - t239 * t451) * mrSges(4,3) + (-t273 * t230 - t240 * t427) * Ifges(4,2) + (t273 * t229 - t274 * t230 + t239 * t427 - t240 * t425) * Ifges(4,4) + t296 * (mrSges(4,1) * t230 + mrSges(4,2) * t229) + t284 * (mrSges(4,1) * t240 + mrSges(4,2) * t239) + t301 * (Ifges(4,5) * t239 - Ifges(4,6) * t240) / 0.2e1 + t221 * t165 + t139 * t465 + m(5) * (t127 * t78 + t128 * t77 + t154 * t52 + t155 * t51 + t197 * t248 + t221 * t241) + (-t337 * t58 - t59 * t84 - t62 * t72 - t63 * t73) * mrSges(6,3) - t337 * t16 - (t55 / 0.2e1 - Ifges(6,4) * t58 + t102 * mrSges(6,1) - t8 * mrSges(6,3) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t59 + t448) * t334 + t77 * t195 + t78 * t196 + m(7) * (t15 * t60 + t2 * t34 + t21 * t5 + t22 * t4 + t3 * t33 - t437) + m(6) * (t102 * t184 + t112 * t171 + t14 * t63 - t15 * t62 + t8 * t84 - t437) + t315 * t73 + t138 * t153 / 0.2e1 + t14 * t144 + t112 * t100 + t5 * t92 + t4 * t91 + (-t125 * t154 + t126 * t155 - t127 * t138 + t128 * t139 + t235 * t51 - t236 * t52) * mrSges(5,3) + t33 * t18 + t34 * t19 + (t262 / 0.2e1 - pkin(7) * t282 + t353 + (-0.2e1 * t438 + 0.3e1 / 0.2e1 * Ifges(3,4) * t311) * qJD(1)) * t311 * qJD(2) + t197 * (-mrSges(5,1) * t235 + mrSges(5,2) * t236) + t239 * t211 / 0.2e1 - t240 * t210 / 0.2e1 + t241 * (-mrSges(5,1) * t139 + mrSges(5,2) * t138) + t175 * t245 + t176 * t246 + (t102 * mrSges(6,2) - Ifges(6,4) * t59 + Ifges(6,1) * t58 + t13 * t422 + t340 * t443 + t344 * t445 + t342 * t444 + t12 * t423 + (mrSges(6,3) + t345) * t9 + (-t2 * t302 - t3 * t307) * mrSges(7,3) + (t76 * t423 + t339 * t433 + t343 * t435 + t341 * t436 + t60 * t346 - t307 * t75 / 0.2e1 - t450 * mrSges(7,3)) * qJD(6)) * t170 + t248 * (-mrSges(5,1) * t126 + mrSges(5,2) * t125) + t300 * (Ifges(5,5) * t138 + Ifges(5,6) * t139) / 0.2e1 + (t313 - t461) * t72 + t184 * (mrSges(6,1) * t59 + mrSges(6,2) * t58) - t390 * t15 + (t274 * t229 + t239 * t425) * Ifges(4,1) + (t236 * t125 + t138 * t429) * Ifges(5,1) + (t235 * t126 + t139 * t431) * Ifges(5,2) + (t235 * t125 + t236 * t126 + t138 * t431 + t139 * t429) * Ifges(5,4) + (-t261 / 0.2e1 - pkin(7) * t283 + t352 + (-0.2e1 * t439 - 0.3e1 / 0.2e1 * t407 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t311) * qJD(1) + (0.2e1 * t421 + t232 + qJD(1) * (-mrSges(4,1) * t273 + mrSges(4,2) * t274)) * pkin(2)) * t370; t453 * t196 + t203 * t324 + t202 * t16 + t452 * t195 + t312 - t25 * t92 - t26 * t91 + ((t245 * t310 - t246 * t305) * qJD(3) + (-t229 * t310 - t230 * t305) * mrSges(4,3)) * pkin(2) - t244 * t165 - t238 * t245 - t237 * t246 + (-t116 * t92 + t203 * t351 + t347) * t302 - t172 * t100 + (-t125 * t258 + t126 * t260 + t382) * mrSges(5,3) + t234 * t409 + t460 * t144 + (-t204 * t58 - t205 * t59 + t394) * mrSges(6,3) + (t116 * t91 + t203 * t19 + (-t203 * t92 - t416) * qJD(6)) * t307 + ((-t297 / 0.2e1 - t262 / 0.2e1 + t353 + qJD(1) * t438 + (t282 - t386) * pkin(7)) * t311 + (t261 / 0.2e1 + t352 + (t439 + t407 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t311) * qJD(1) + (t283 + t385) * pkin(7) + (-t232 - t421) * pkin(2)) * t306) * qJD(1) - t459 * t390 + (t116 * t450 + t202 * t9 - t21 * t25 - t22 * t26 + t459 * t60) * m(7) + (-t171 * t172 - t204 * t9 + t205 * t8 - t459 * t62 + t460 * t63) * m(6) + (t127 * t453 + t128 * t452 - t241 * t244 + t258 * t52 + t260 * t51) * m(5) + (-t451 * t237 - t234 * t238 + (t167 * t305 + t168 * t310 + (t234 * t310 - t305 * t451) * qJD(3)) * pkin(2)) * m(4); (t246 + t409) * t234 - m(5) * (t127 * t133 + t128 * t134) + (-t224 * t92 + t254 * t351 + t347) * t302 + t458 * t144 + (t254 * t19 + t224 * t91 + (-t254 * t92 - t416) * qJD(6)) * t307 - t134 * t195 - t133 * t196 + (-t257 * t58 - t259 * t59 + t394) * mrSges(6,3) + t254 * t324 + t312 - t23 * t92 - t24 * t91 - t451 * t245 + t252 * t16 - t179 * t100 + mrSges(5,3) * t382 - t390 * t457 + (-t21 * t23 - t22 * t24 + t224 * t450 + t252 * t9 + t457 * t60) * m(7) + (-t171 * t179 - t257 * t9 + t259 * t8 - t457 * t62 + t458 * t63) * m(6) + (-t264 * t165 + (t195 * t309 - t196 * t304) * qJD(4) + (-t125 * t309 + t126 * t304) * mrSges(5,3) + (-t127 * t367 + t128 * t366 + 0.2e1 * t241 * t426 + t304 * t51 + t309 * t52) * m(5)) * pkin(3); (-t157 * t392 + (-t157 * t22 - t3) * t302) * mrSges(7,3) + (-t315 + t411) * t469 + (-t313 + t412) * t470 - m(7) * (t21 * t27 + t22 * t28 + t60 * t64) + t152 * t429 + t471 * (-pkin(4) * t308 - pkin(5)) + t316 * (pkin(4) * t303 + pkin(11)) - t127 * t195 + t128 * t196 - m(6) * (-t62 * t64 + t63 * t65) - t65 * t144 - t28 * t91 - t27 * t92 + (-t220 * t100 + (-t303 * t59 - t308 * t58) * mrSges(6,3) + (-t390 * t303 + (-t302 * t92 + t307 * t91 + t144) * t308 + m(7) * (t303 * t60 + t308 * t450)) * qJD(5) + (t303 * t8 - t308 * t9 + 0.2e1 * t171 * t430 + (-t303 * t62 + t308 * t63) * qJD(5)) * m(6)) * pkin(4) + (t456 + t382) * mrSges(5,3) + t390 * t64 + t473; t316 * pkin(11) + (t318 + t411) * t469 - m(7) * (t21 * t29 + t22 * t30 + t60 * t63) + (t412 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t469 + t461 - t449) * t470 - t62 * t144 - t30 * t91 - t29 * t92 + t326 * mrSges(7,3) + t390 * t63 + t322 - t471 * pkin(5); t55 - t60 * (mrSges(7,1) * t143 + mrSges(7,2) * t142) + (Ifges(7,1) * t142 - t395) * t435 + t75 * t434 + (Ifges(7,5) * t142 - Ifges(7,6) * t143) * t433 - t21 * t91 + t22 * t92 + (t142 * t21 + t143 * t22) * mrSges(7,3) + (-Ifges(7,2) * t143 + t141 + t76) * t436 + t448;];
tauc  = t1(:);
