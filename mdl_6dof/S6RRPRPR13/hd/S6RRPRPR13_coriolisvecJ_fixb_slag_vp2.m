% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:09:12
% EndTime: 2018-11-23 17:09:26
% DurationCPUTime: 13.93s
% Computational Cost: add. (12894->822), mult. (33365->1139), div. (0->0), fcn. (24832->10), ass. (0->359)
t310 = cos(qJ(4));
t307 = sin(qJ(4));
t332 = pkin(4) * t310 + qJ(5) * t307;
t244 = qJD(4) * t332 - qJD(5) * t310 + qJD(3);
t302 = sin(pkin(11));
t304 = cos(pkin(11));
t312 = -pkin(2) - pkin(9);
t366 = qJD(4) * t310;
t350 = t312 * t366;
t187 = t302 * t244 + t304 * t350;
t308 = sin(qJ(2));
t303 = sin(pkin(6));
t372 = qJD(1) * t303;
t355 = t308 * t372;
t287 = pkin(2) * t355;
t311 = cos(qJ(2));
t331 = pkin(9) * t308 - qJ(3) * t311;
t216 = t331 * t372 + t287;
t354 = t311 * t372;
t305 = cos(pkin(6));
t371 = qJD(1) * t305;
t363 = pkin(1) * t371;
t252 = pkin(8) * t354 + t308 * t363;
t218 = pkin(3) * t354 + t252;
t125 = t310 * t216 + t307 * t218;
t113 = qJ(5) * t354 + t125;
t292 = t311 * t363;
t433 = pkin(3) + pkin(8);
t144 = t292 + (-t332 - t433) * t355;
t74 = t304 * t113 + t302 * t144;
t467 = t187 - t74;
t382 = t307 * t308;
t224 = (t302 * t311 + t304 * t382) * t372;
t232 = t304 * t244;
t342 = t310 * t355;
t347 = -t302 * t312 + pkin(5);
t414 = pkin(10) * t304;
t73 = -t113 * t302 + t304 * t144;
t466 = pkin(5) * t342 + pkin(10) * t224 + t232 + (t307 * t414 + t310 * t347) * qJD(4) - t73;
t223 = (-t302 * t382 + t304 * t311) * t372;
t367 = qJD(4) * t307;
t465 = t467 + (t302 * t367 - t223) * pkin(10);
t294 = qJD(2) + t371;
t461 = -t294 / 0.2e1;
t464 = -mrSges(3,1) + mrSges(4,2);
t463 = mrSges(3,3) + mrSges(4,1);
t341 = t307 * t354;
t234 = t294 * t310 - t341;
t282 = qJD(4) + t355;
t167 = t234 * t304 + t282 * t302;
t306 = sin(qJ(6));
t309 = cos(qJ(6));
t345 = -t234 * t302 + t304 * t282;
t462 = -t167 * t306 + t309 * t345;
t99 = t167 * t309 + t306 * t345;
t370 = qJD(2) * t308;
t351 = t307 * t370;
t174 = -t294 * t367 + (-t311 * t366 + t351) * t372;
t349 = qJD(2) * t372;
t339 = t311 * t349;
t138 = -t174 * t302 + t304 * t339;
t139 = t174 * t304 + t302 * t339;
t36 = qJD(6) * t462 + t138 * t306 + t139 * t309;
t440 = t36 / 0.2e1;
t37 = -qJD(6) * t99 + t138 * t309 - t139 * t306;
t439 = t37 / 0.2e1;
t340 = t308 * t349;
t175 = -qJD(4) * t341 + t294 * t366 - t310 * t340;
t426 = t175 / 0.2e1;
t460 = -t372 / 0.2e1;
t270 = pkin(4) * t307 - qJ(5) * t310 + qJ(3);
t263 = t304 * t270;
t413 = pkin(10) * t310;
t182 = -t304 * t413 + t307 * t347 + t263;
t381 = t307 * t312;
t227 = t302 * t270 + t304 * t381;
t192 = -t302 * t413 + t227;
t108 = t182 * t309 - t192 * t306;
t459 = qJD(6) * t108 + t306 * t466 + t465 * t309;
t109 = t182 * t306 + t192 * t309;
t458 = -qJD(6) * t109 - t465 * t306 + t309 * t466;
t410 = pkin(10) + qJ(5);
t274 = t410 * t302;
t275 = t410 * t304;
t207 = -t274 * t309 - t275 * t306;
t323 = t302 * t306 - t304 * t309;
t233 = t294 * t307 + t310 * t354;
t152 = t294 * t312 + t355 * t433 + qJD(3) - t292;
t348 = -qJ(3) * t308 - pkin(1);
t210 = (t311 * t312 + t348) * t303;
t188 = qJD(1) * t210;
t102 = t152 * t310 - t307 * t188;
t147 = pkin(4) * t234 + qJ(5) * t233;
t67 = -t102 * t302 + t304 * t147;
t41 = pkin(5) * t234 + t233 * t414 + t67;
t385 = t233 * t302;
t68 = t304 * t102 + t302 * t147;
t50 = pkin(10) * t385 + t68;
t457 = -qJD(5) * t323 + qJD(6) * t207 - t306 * t41 - t309 * t50;
t208 = -t274 * t306 + t275 * t309;
t267 = t302 * t309 + t304 * t306;
t456 = -qJD(5) * t267 - qJD(6) * t208 + t306 * t50 - t309 * t41;
t124 = -t307 * t216 + t218 * t310;
t114 = -pkin(4) * t354 - t124;
t346 = pkin(5) * t302 - t312;
t455 = pkin(5) * t223 - t346 * t367 - t114;
t221 = t294 * t304 + t302 * t342;
t222 = t294 * t302 - t304 * t342;
t243 = t323 * t310;
t257 = t267 * qJD(6);
t452 = -qJD(4) * t243 - t221 * t306 - t222 * t309 - t257 * t307;
t241 = t267 * t310;
t449 = qJD(6) * t323;
t451 = -qJD(4) * t241 - t221 * t309 + t222 * t306 + t307 * t449;
t228 = Ifges(5,4) * t233;
t394 = t234 * Ifges(5,1);
t450 = -t394 / 0.2e1 + t228 / 0.2e1;
t251 = pkin(8) * t355 - t292;
t448 = -qJD(3) - t251;
t280 = pkin(2) * t340;
t368 = qJD(3) * t308;
t317 = (qJD(2) * t331 - t368) * t303;
t162 = qJD(1) * t317 + t280;
t298 = t305 * t308 * pkin(1);
t383 = t303 * t311;
t219 = (t383 * t433 + t298) * qJD(2);
t193 = qJD(1) * t219;
t51 = t152 * t366 + t310 * t162 - t188 * t367 + t307 * t193;
t52 = -t152 * t367 - t307 * t162 - t188 * t366 + t193 * t310;
t447 = t52 * mrSges(5,1) - t51 * mrSges(5,2) + Ifges(5,5) * t174 - Ifges(5,6) * t175;
t446 = -Ifges(3,4) * t354 / 0.2e1 + Ifges(3,5) * t461;
t392 = t282 * Ifges(5,5);
t128 = -t228 + t392 + t394;
t284 = t294 * qJ(3);
t179 = t284 + t218;
t318 = t102 * mrSges(5,3) - t128 / 0.2e1 - t179 * mrSges(5,2) - t392 / 0.2e1;
t407 = Ifges(6,4) * t304;
t334 = -Ifges(6,2) * t302 + t407;
t408 = Ifges(6,4) * t302;
t335 = Ifges(6,1) * t304 - t408;
t336 = mrSges(6,1) * t302 + mrSges(6,2) * t304;
t417 = t304 / 0.2e1;
t418 = -t302 / 0.2e1;
t427 = t167 / 0.2e1;
t428 = t345 / 0.2e1;
t105 = pkin(4) * t233 - qJ(5) * t234 + t179;
t103 = t152 * t307 + t188 * t310;
t88 = qJ(5) * t282 + t103;
t43 = t304 * t105 - t302 * t88;
t44 = t302 * t105 + t304 * t88;
t83 = t167 * Ifges(6,4) + Ifges(6,2) * t345 + Ifges(6,6) * t233;
t84 = t167 * Ifges(6,1) + Ifges(6,4) * t345 + Ifges(6,5) * t233;
t87 = -pkin(4) * t282 + qJD(5) - t102;
t445 = -(t302 * t44 + t304 * t43) * mrSges(6,3) + t334 * t428 + t335 * t427 + t336 * t87 + t417 * t84 + t418 * t83 - t318;
t444 = Ifges(7,4) * t440 + Ifges(7,2) * t439 + Ifges(7,6) * t426;
t443 = Ifges(7,1) * t440 + Ifges(7,4) * t439 + Ifges(7,5) * t426;
t442 = Ifges(3,5) / 0.2e1;
t441 = Ifges(4,5) / 0.2e1;
t61 = t139 * Ifges(6,1) + t138 * Ifges(6,4) + t175 * Ifges(6,5);
t438 = t61 / 0.2e1;
t437 = -t462 / 0.2e1;
t436 = t462 / 0.2e1;
t435 = -t99 / 0.2e1;
t434 = t99 / 0.2e1;
t432 = pkin(1) * mrSges(3,1);
t431 = pkin(1) * mrSges(3,2);
t430 = t138 / 0.2e1;
t429 = t139 / 0.2e1;
t230 = qJD(6) + t233;
t424 = -t230 / 0.2e1;
t423 = t230 / 0.2e1;
t422 = -t233 / 0.2e1;
t258 = t305 * t307 + t310 * t383;
t421 = -t258 / 0.2e1;
t359 = t307 * t383;
t259 = t305 * t310 - t359;
t419 = t259 / 0.2e1;
t416 = Ifges(7,4) * t99;
t32 = Ifges(7,5) * t36;
t31 = Ifges(7,6) * t37;
t415 = pkin(1) * t311;
t412 = t462 * Ifges(7,6);
t411 = t99 * Ifges(7,5);
t42 = qJ(5) * t339 + qJD(5) * t282 + t51;
t281 = qJD(2) * t292;
t283 = t294 * qJD(3);
t384 = t303 * t308;
t344 = t433 * t384;
t326 = qJD(2) * t344;
t163 = -qJD(1) * t326 + t281 + t283;
t69 = pkin(4) * t175 - qJ(5) * t174 - qJD(5) * t234 + t163;
t19 = t302 * t69 + t304 * t42;
t369 = qJD(2) * t311;
t353 = t303 * t370;
t289 = pkin(2) * t353;
t184 = t289 + t317;
t295 = pkin(8) * t384;
t356 = -pkin(2) - t415;
t191 = pkin(3) * t384 + t295 + (-pkin(9) + t356) * t305;
t71 = t310 * t184 + t191 * t366 - t210 * t367 + t307 * t219;
t58 = (qJ(5) * t369 + qJD(5) * t308) * t303 + t71;
t365 = t305 * t415;
t293 = qJD(2) * t365;
t300 = t305 * qJD(3);
t190 = t293 + t300 - t326;
t204 = -qJD(4) * t258 + t303 * t351;
t205 = -qJD(4) * t359 + t305 * t366 - t310 * t353;
t80 = pkin(4) * t205 - qJ(5) * t204 - qJD(5) * t259 + t190;
t24 = t302 * t80 + t304 * t58;
t409 = Ifges(5,4) * t234;
t406 = Ifges(6,5) * t304;
t405 = Ifges(4,6) * t311;
t404 = Ifges(6,6) * t302;
t403 = t138 * Ifges(6,6);
t402 = t139 * Ifges(6,5);
t401 = t345 * Ifges(6,6);
t400 = t167 * Ifges(6,5);
t399 = t174 * Ifges(5,1);
t398 = t174 * Ifges(5,4);
t397 = t175 * Ifges(5,4);
t396 = t230 * Ifges(7,3);
t391 = t282 * Ifges(5,6);
t388 = t294 * Ifges(4,5);
t387 = t307 * t87;
t145 = mrSges(5,1) * t339 - mrSges(5,3) * t174;
t81 = -t138 * mrSges(6,1) + t139 * mrSges(6,2);
t386 = -t81 + t145;
t380 = t310 * t312;
t120 = t307 * t191 + t310 * t210;
t110 = qJ(5) * t384 + t120;
t261 = pkin(8) * t383 + t298;
t237 = -t305 * qJ(3) - t261;
t209 = pkin(3) * t383 - t237;
t123 = pkin(4) * t258 - qJ(5) * t259 + t209;
t64 = t304 * t110 + t302 * t123;
t131 = t223 * t309 - t224 * t306;
t160 = t267 * t367 + t310 * t449;
t379 = t131 - t160;
t132 = t223 * t306 + t224 * t309;
t158 = -t257 * t310 + t323 * t367;
t378 = t132 - t158;
t133 = t267 * t233;
t377 = t133 + t257;
t134 = t323 * t233;
t376 = t134 + t449;
t104 = -mrSges(6,1) * t345 + mrSges(6,2) * t167;
t181 = mrSges(5,1) * t282 - mrSges(5,3) * t234;
t375 = t181 - t104;
t148 = mrSges(5,1) * t233 + mrSges(5,2) * t234;
t247 = -mrSges(4,1) * t354 - mrSges(4,3) * t294;
t374 = -t247 + t148;
t373 = t294 * t464 + t355 * t463;
t364 = t294 / 0.2e1 - qJD(2);
t7 = Ifges(7,3) * t175 + t31 + t32;
t362 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t361 = Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1;
t352 = t303 * t369;
t11 = -t37 * mrSges(7,1) + t36 * mrSges(7,2);
t18 = -t302 * t42 + t304 * t69;
t23 = -t302 * t58 + t304 * t80;
t63 = -t110 * t302 + t304 * t123;
t119 = t191 * t310 - t307 * t210;
t343 = t307 * t355;
t10 = pkin(5) * t175 - pkin(10) * t139 + t18;
t12 = pkin(10) * t138 + t19;
t25 = pkin(5) * t233 - pkin(10) * t167 + t43;
t29 = pkin(10) * t345 + t44;
t5 = t25 * t309 - t29 * t306;
t1 = qJD(6) * t5 + t10 * t306 + t12 * t309;
t6 = t25 * t306 + t29 * t309;
t2 = -qJD(6) * t6 + t10 * t309 - t12 * t306;
t338 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t333 = -t404 + t406;
t330 = -t18 * t302 + t19 * t304;
t328 = -t302 * t43 + t304 * t44;
t90 = -mrSges(6,2) * t175 + mrSges(6,3) * t138;
t91 = mrSges(6,1) * t175 - mrSges(6,3) * t139;
t327 = -t302 * t91 + t304 * t90;
t203 = t259 * t304 + t302 * t384;
t38 = pkin(5) * t258 - pkin(10) * t203 + t63;
t202 = -t259 * t302 + t304 * t384;
t47 = pkin(10) * t202 + t64;
t13 = -t306 * t47 + t309 * t38;
t14 = t306 * t38 + t309 * t47;
t325 = t102 * t307 - t103 * t310;
t121 = -mrSges(6,2) * t233 + mrSges(6,3) * t345;
t122 = mrSges(6,1) * t233 - mrSges(6,3) * t167;
t324 = t121 * t304 - t122 * t302;
t115 = t202 * t309 - t203 * t306;
t116 = t202 * t306 + t203 * t309;
t253 = -pkin(8) * t353 + t293;
t72 = -t307 * t184 - t191 * t367 - t210 * t366 + t219 * t310;
t235 = -pkin(8) * t340 + t281;
t112 = -pkin(4) * t384 - t119;
t238 = (-pkin(2) * t311 + t348) * t303;
t321 = (-qJ(3) * t369 - t368) * t303;
t225 = qJD(1) * t238;
t320 = Ifges(3,6) * t461 + (Ifges(3,4) * t308 + Ifges(3,2) * t311) * t460 + t388 / 0.2e1 + (-Ifges(4,6) * t308 - Ifges(4,3) * t311) * t372 / 0.2e1 - t225 * mrSges(4,2) - t252 * mrSges(3,3);
t254 = t261 * qJD(2);
t199 = -t235 - t283;
t236 = qJD(1) * t254;
t319 = -t235 * mrSges(3,2) - t199 * mrSges(4,3) + t236 * t464;
t62 = -pkin(4) * t352 - t72;
t46 = -pkin(4) * t339 - t52;
t316 = t102 * mrSges(5,1) + t251 * mrSges(3,3) + t282 * Ifges(5,3) + t234 * Ifges(5,5) - t233 * Ifges(5,6) + Ifges(3,1) * t355 / 0.2e1 + Ifges(4,4) * t461 + (-t308 * Ifges(4,2) - t405) * t460 - t103 * mrSges(5,2) - t225 * mrSges(4,3) - t446;
t127 = -t233 * Ifges(5,2) + t391 + t409;
t33 = t396 + t411 + t412;
t82 = t233 * Ifges(6,3) + t400 + t401;
t314 = -t401 / 0.2e1 - t400 / 0.2e1 + t391 / 0.2e1 - t5 * mrSges(7,1) + t6 * mrSges(7,2) - t43 * mrSges(6,1) + t44 * mrSges(6,2) - t179 * mrSges(5,1) - t396 / 0.2e1 - t412 / 0.2e1 - t411 / 0.2e1 - t82 / 0.2e1 - t33 / 0.2e1 + t127 / 0.2e1 + t103 * mrSges(5,3) + t409 / 0.2e1;
t313 = (t233 * t361 - t314) * t310;
t299 = -pkin(5) * t304 - pkin(4);
t278 = Ifges(3,5) * t339;
t277 = Ifges(4,5) * t340;
t276 = Ifges(5,3) * t339;
t264 = t346 * t310;
t260 = -t295 + t365;
t250 = -qJ(3) * t354 + t287;
t249 = (mrSges(4,2) * t311 - mrSges(4,3) * t308) * t372;
t246 = -mrSges(3,2) * t294 + mrSges(3,3) * t354;
t242 = t323 * t307;
t240 = t267 * t307;
t239 = t305 * t356 + t295;
t229 = -t253 - t300;
t226 = -t302 * t381 + t263;
t220 = t289 + t321;
t217 = -qJD(1) * t344 + t292;
t215 = -t284 - t252;
t206 = -pkin(2) * t294 - t448;
t196 = qJD(1) * t321 + t280;
t186 = -t302 * t350 + t232;
t180 = -mrSges(5,2) * t282 - mrSges(5,3) * t233;
t154 = t204 * t304 + t302 * t352;
t153 = -t204 * t302 + t304 * t352;
t146 = -mrSges(5,2) * t339 - mrSges(5,3) * t175;
t107 = mrSges(5,1) * t175 + mrSges(5,2) * t174;
t95 = Ifges(5,5) * t339 - t397 + t399;
t94 = -t175 * Ifges(5,2) + Ifges(5,6) * t339 + t398;
t93 = Ifges(7,4) * t462;
t85 = -pkin(5) * t202 + t112;
t79 = -pkin(5) * t385 + t103;
t78 = mrSges(7,1) * t230 - mrSges(7,3) * t99;
t77 = -mrSges(7,2) * t230 + mrSges(7,3) * t462;
t70 = -pkin(5) * t345 + t87;
t60 = t139 * Ifges(6,4) + t138 * Ifges(6,2) + t175 * Ifges(6,6);
t59 = t175 * Ifges(6,3) + t402 + t403;
t54 = -qJD(6) * t116 + t153 * t309 - t154 * t306;
t53 = qJD(6) * t115 + t153 * t306 + t154 * t309;
t45 = -mrSges(7,1) * t462 + mrSges(7,2) * t99;
t40 = -pkin(5) * t153 + t62;
t35 = Ifges(7,1) * t99 + Ifges(7,5) * t230 + t93;
t34 = Ifges(7,2) * t462 + Ifges(7,6) * t230 + t416;
t28 = -pkin(5) * t138 + t46;
t27 = -mrSges(7,2) * t175 + mrSges(7,3) * t37;
t26 = mrSges(7,1) * t175 - mrSges(7,3) * t36;
t22 = pkin(10) * t153 + t24;
t17 = pkin(5) * t205 - pkin(10) * t154 + t23;
t4 = -qJD(6) * t14 + t17 * t309 - t22 * t306;
t3 = qJD(6) * t13 + t17 * t306 + t22 * t309;
t8 = [(Ifges(6,5) * t203 + Ifges(7,5) * t116 + Ifges(6,6) * t202 + Ifges(7,6) * t115 + (Ifges(6,3) + Ifges(7,3)) * t258) * t426 + (t59 + t7) * t258 / 0.2e1 + (t82 + t33) * t205 / 0.2e1 + m(4) * (t196 * t238 + t199 * t237 + t206 * t254 + t215 * t229 + t220 * t225 + t236 * t239) + m(3) * (t235 * t261 - t236 * t260 + t251 * t254 + t252 * t253) + m(5) * (t102 * t72 + t103 * t71 + t119 * t52 + t120 * t51 + t163 * t209 + t179 * t190) + m(7) * (t1 * t14 + t13 * t2 + t28 * t85 + t3 * t6 + t4 * t5 + t40 * t70) + m(6) * (t112 * t46 + t18 * t63 + t19 * t64 + t23 * t43 + t24 * t44 + t62 * t87) + (-t102 * t204 - t103 * t205 - t258 * t51 - t259 * t52) * mrSges(5,3) + t282 * (Ifges(5,5) * t204 - Ifges(5,6) * t205) / 0.2e1 + (t277 / 0.2e1 + t278 / 0.2e1 + t319) * t305 + ((-t238 * mrSges(4,3) + Ifges(5,5) * t419 + Ifges(5,6) * t421 - t260 * mrSges(3,3) + t239 * mrSges(4,1) + (-Ifges(4,4) + t442) * t305 + (t311 * t362 - 0.2e1 * t431) * t303) * t311 + (t237 * mrSges(4,1) - t238 * mrSges(4,2) - t261 * mrSges(3,3) + (-Ifges(3,6) + t441) * t305 + (-t308 * t362 - 0.2e1 * t432) * t303 + (0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + Ifges(5,3) / 0.2e1) * t383) * t308) * t349 + t163 * (mrSges(5,1) * t258 + mrSges(5,2) * t259) + t18 * (mrSges(6,1) * t258 - mrSges(6,3) * t203) + t19 * (-mrSges(6,2) * t258 + mrSges(6,3) * t202) + t1 * (-mrSges(7,2) * t258 + mrSges(7,3) * t115) + t2 * (mrSges(7,1) * t258 - mrSges(7,3) * t116) + t253 * t246 + t229 * t247 + t220 * t249 + t233 * (Ifges(6,5) * t154 + Ifges(6,6) * t153 + Ifges(6,3) * t205) / 0.2e1 + t179 * (mrSges(5,1) * t205 + mrSges(5,2) * t204) - t205 * t127 / 0.2e1 + t209 * t107 + t46 * (-mrSges(6,1) * t202 + mrSges(6,2) * t203) + t204 * t128 / 0.2e1 + t43 * (mrSges(6,1) * t205 - mrSges(6,3) * t154) + t44 * (-mrSges(6,2) * t205 + mrSges(6,3) * t153) + t6 * (-mrSges(7,2) * t205 + mrSges(7,3) * t54) + t5 * (mrSges(7,1) * t205 - mrSges(7,3) * t53) + t202 * t60 / 0.2e1 + t373 * t254 + t190 * t148 + t71 * t180 + t72 * t181 + t87 * (-mrSges(6,1) * t153 + mrSges(6,2) * t154) + t154 * t84 / 0.2e1 + t153 * t83 / 0.2e1 + t119 * t145 + t120 * t146 + t23 * t122 + t24 * t121 + t28 * (-mrSges(7,1) * t115 + mrSges(7,2) * t116) + t112 * t81 + t62 * t104 + t63 * t91 + t85 * t11 + t64 * t90 + t3 * t77 + t4 * t78 + t70 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t234 * (Ifges(5,1) * t204 - Ifges(5,4) * t205) / 0.2e1 + t174 * (Ifges(5,1) * t259 - Ifges(5,4) * t258) / 0.2e1 + t54 * t34 / 0.2e1 + t53 * t35 / 0.2e1 + t40 * t45 + t13 * t26 + t14 * t27 + t95 * t419 + t94 * t421 + (Ifges(5,4) * t204 - Ifges(5,2) * t205) * t422 + (Ifges(7,5) * t53 + Ifges(7,6) * t54 + Ifges(7,3) * t205) * t423 + (Ifges(6,1) * t154 + Ifges(6,4) * t153 + Ifges(6,5) * t205) * t427 + (Ifges(6,4) * t154 + Ifges(6,2) * t153 + Ifges(6,6) * t205) * t428 + (Ifges(6,1) * t203 + Ifges(6,4) * t202 + Ifges(6,5) * t258) * t429 + (Ifges(6,4) * t203 + Ifges(6,2) * t202 + Ifges(6,6) * t258) * t430 + (Ifges(7,1) * t53 + Ifges(7,4) * t54 + Ifges(7,5) * t205) * t434 + (Ifges(7,4) * t53 + Ifges(7,2) * t54 + Ifges(7,6) * t205) * t436 + t203 * t438 + (Ifges(7,4) * t116 + Ifges(7,2) * t115 + Ifges(7,6) * t258) * t439 + (Ifges(7,1) * t116 + Ifges(7,4) * t115 + Ifges(7,5) * t258) * t440 + t116 * t443 + t115 * t444 + ((-mrSges(4,1) * t199 + mrSges(4,2) * t196 + mrSges(3,3) * t235) * t311 + (t276 / 0.2e1 - t196 * mrSges(4,3) + t463 * t236 + t447) * t308 + ((t215 * mrSges(4,1) + (t441 - Ifges(3,6) / 0.2e1) * t294 + t320) * t308 + (t206 * mrSges(4,1) + (-Ifges(4,4) / 0.2e1 + t442) * t294 + t316) * t311) * qJD(2)) * t303 - t175 * (Ifges(5,4) * t259 - Ifges(5,2) * t258) / 0.2e1; -t345 * (Ifges(6,4) * t224 + Ifges(6,2) * t223) / 0.2e1 + t28 * (mrSges(7,1) * t241 - mrSges(7,2) * t243) + (-t1 * t241 + t2 * t243 + t378 * t5 - t379 * t6) * mrSges(7,3) + (-Ifges(7,5) * t243 - Ifges(7,6) * t241) * t426 + (-Ifges(7,4) * t243 - Ifges(7,2) * t241) * t439 + (-Ifges(7,1) * t243 - Ifges(7,4) * t241) * t440 + (t186 - t73) * t122 + (-pkin(2) * t236 - qJ(3) * t199 - t206 * t252 + t215 * t448 - t225 * t250) * m(4) + (t313 + (t333 * t422 - t445 + t450) * t307 + (-m(5) * t325 + m(6) * t387 + t310 * t180 - t307 * t375) * t312) * qJD(4) + (((t431 - t405 / 0.2e1) * t372 + t364 * Ifges(4,4) + (-pkin(2) * qJD(2) - t206) * mrSges(4,1) + qJD(2) * (Ifges(5,5) * t310 - Ifges(5,6) * t307) / 0.2e1 - t316 + t446) * t311 + (-t388 / 0.2e1 + ((t432 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t308) * t303 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t383) * qJD(1) + t364 * Ifges(3,6) + (-qJ(3) * qJD(2) - t215) * mrSges(4,1) + (t318 + t450) * t307 + t313 - t320) * t308) * t372 + (-t247 + t246) * t251 + t458 * t78 + (t1 * t109 + t108 * t2 + t264 * t28 + t455 * t70 + t458 * t5 + t459 * t6) * m(7) + t459 * t77 + (-t223 * t44 + t224 * t43) * mrSges(6,3) + t467 * t121 - m(5) * (t102 * t124 + t103 * t125 + t179 * t217) - m(6) * (t114 * t87 + t43 * t73 + t44 * t74) + t277 + t278 + t319 + t455 * t45 + (t333 * t426 + t334 * t430 + t335 * t429 + t46 * t336 + t95 / 0.2e1 - t397 / 0.2e1 + t399 / 0.2e1 + t163 * mrSges(5,2) + t61 * t417 + t60 * t418 - t52 * mrSges(5,3) + t386 * t312 + (-t18 * t304 - t19 * t302) * mrSges(6,3)) * t310 + (t18 * mrSges(6,1) + t403 / 0.2e1 + t402 / 0.2e1 - t19 * mrSges(6,2) - t398 / 0.2e1 + t163 * mrSges(5,1) - t94 / 0.2e1 + t59 / 0.2e1 + t7 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1 + t312 * t146 - t51 * mrSges(5,3) + (Ifges(7,3) / 0.2e1 + t361) * t175 + t338) * t307 + t264 * t11 - t250 * t249 + m(6) * (t18 * t226 + t43 * t186 + t44 * t187 + t19 * t227 - t380 * t46) + m(5) * (t163 * qJ(3) + t179 * qJD(3) + t380 * t52 + t381 * t51) + (mrSges(7,1) * t379 - mrSges(7,2) * t378) * t70 + t226 * t91 + t227 * t90 - t223 * t83 / 0.2e1 - t87 * (-mrSges(6,1) * t223 + mrSges(6,2) * t224) - t224 * t84 / 0.2e1 - t217 * t148 - t373 * t252 + t374 * qJD(3) - t125 * t180 - t124 * t181 - t167 * (Ifges(6,1) * t224 + Ifges(6,4) * t223) / 0.2e1 + t109 * t27 - t114 * t104 + qJ(3) * t107 + t108 * t26 + (t158 / 0.2e1 - t132 / 0.2e1) * t35 + (t160 / 0.2e1 - t131 / 0.2e1) * t34 + (Ifges(6,5) * t224 + Ifges(6,6) * t223) * t422 + (Ifges(7,5) * t158 + Ifges(7,6) * t160) * t423 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t424 + (Ifges(7,1) * t158 + Ifges(7,4) * t160) * t434 + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t435 + (Ifges(7,4) * t158 + Ifges(7,2) * t160) * t436 + (Ifges(7,4) * t132 + Ifges(7,2) * t131) * t437 - t243 * t443 - t241 * t444; -t222 * t121 - t221 * t122 - t240 * t26 - t242 * t27 + t451 * t78 + t452 * t77 - t374 * t294 + (mrSges(4,1) * t369 + t249 * t308) * t372 + (t180 * t355 - t11 + (t180 + t324) * qJD(4) + t386) * t310 + (t146 + t327 + t282 * (t45 - t375)) * t307 + (-t1 * t242 - t2 * t240 - t28 * t310 + (t343 + t367) * t70 + t452 * t6 + t451 * t5) * m(7) + (-t310 * t46 + t330 * t307 + (t310 * t328 + t387) * qJD(4) - t221 * t43 - t222 * t44 + t87 * t343) * m(6) + (-t179 * t294 - t282 * t325 + t307 * t51 + t310 * t52) * m(5) + (t215 * t294 + t225 * t355 + t236) * m(4); (-t449 / 0.2e1 - t134 / 0.2e1) * t35 + (-Ifges(7,5) * t449 - Ifges(7,6) * t257) * t423 + (-Ifges(7,1) * t449 - Ifges(7,4) * t257) * t434 + (-Ifges(7,4) * t449 - Ifges(7,2) * t257) * t436 + (Ifges(6,5) * t302 + Ifges(7,5) * t267 + Ifges(6,6) * t304 - Ifges(7,6) * t323) * t426 + t28 * (mrSges(7,1) * t323 + mrSges(7,2) * t267) + (-t1 * t323 - t2 * t267 + t376 * t5 - t377 * t6) * mrSges(7,3) + (Ifges(7,4) * t267 - Ifges(7,2) * t323) * t439 + (Ifges(7,1) * t267 - Ifges(7,4) * t323) * t440 - t323 * t444 + (-t257 / 0.2e1 - t133 / 0.2e1) * t34 + (-pkin(4) * t46 + qJ(5) * t330 + qJD(5) * t328 - t103 * t87 - t43 * t67 - t44 * t68) * m(6) + t314 * t234 + (-t228 / 0.2e1 + (t406 / 0.2e1 - t404 / 0.2e1) * t233 + (Ifges(5,1) / 0.2e1 - t361) * t234 + t445) * t233 + t447 + t276 + t456 * t78 + t457 * t77 + (t1 * t208 + t2 * t207 + t28 * t299 + t456 * t5 + t457 * t6 - t70 * t79) * m(7) + t46 * (-mrSges(6,1) * t304 + mrSges(6,2) * t302) + t299 * t11 + t207 * t26 + t208 * t27 + (mrSges(7,1) * t377 - mrSges(7,2) * t376) * t70 + t375 * t103 - t102 * t180 - t68 * t121 - t67 * t122 - t79 * t45 - pkin(4) * t81 + t60 * t417 + (Ifges(7,5) * t134 + Ifges(7,6) * t133) * t424 + (Ifges(6,1) * t302 + t407) * t429 + (Ifges(6,2) * t304 + t408) * t430 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t435 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t437 + t302 * t438 + t267 * t443 + t324 * qJD(5) + t327 * qJ(5) + t330 * mrSges(6,3); -t345 * t121 + t167 * t122 - t462 * t77 + t99 * t78 + t11 + t81 + (-t462 * t6 + t5 * t99 + t28) * m(7) + (t167 * t43 - t345 * t44 + t46) * m(6); -t70 * (mrSges(7,1) * t99 + mrSges(7,2) * t462) + (Ifges(7,1) * t462 - t416) * t435 + t34 * t434 + (Ifges(7,5) * t462 - Ifges(7,6) * t99) * t424 - t5 * t77 + t6 * t78 + (t462 * t5 + t6 * t99) * mrSges(7,3) + t338 + t7 + (-Ifges(7,2) * t99 + t35 + t93) * t437;];
tauc  = t8(:);
