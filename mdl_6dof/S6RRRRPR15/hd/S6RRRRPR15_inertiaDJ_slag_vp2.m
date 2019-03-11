% Calculate time derivative of joint inertia matrix for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR15_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:18
% EndTime: 2019-03-10 00:34:43
% DurationCPUTime: 10.68s
% Computational Cost: add. (17241->924), mult. (51811->1290), div. (0->0), fcn. (51573->12), ass. (0->358)
t312 = sin(pkin(6));
t452 = 0.2e1 * t312;
t322 = cos(qJ(2));
t314 = cos(pkin(6));
t417 = pkin(1) * t314;
t305 = t322 * t417;
t294 = qJD(2) * t305;
t318 = sin(qJ(2));
t313 = cos(pkin(7));
t353 = t312 * (-pkin(10) * t313 - pkin(9));
t336 = t318 * t353;
t194 = pkin(2) * t314 + t305 + t336;
t401 = t194 * t313;
t451 = qJD(2) * t336 + qJD(3) * t401 + t294;
t315 = sin(qJ(6));
t319 = cos(qJ(6));
t320 = cos(qJ(4));
t379 = qJD(6) * t320;
t316 = sin(qJ(4));
t381 = qJD(4) * t316;
t329 = t315 * t379 + t319 * t381;
t328 = t315 * t381 - t319 * t379;
t450 = -Ifges(6,4) / 0.2e1;
t449 = -t320 / 0.2e1;
t321 = cos(qJ(3));
t390 = t318 * t321;
t317 = sin(qJ(3));
t391 = t317 * t322;
t333 = t313 * t391 + t390;
t311 = sin(pkin(7));
t383 = qJD(3) * t317;
t370 = t311 * t383;
t142 = t314 * t370 + (t333 * qJD(3) + (t313 * t390 + t391) * qJD(2)) * t312;
t382 = qJD(3) * t321;
t369 = t311 * t382;
t386 = t321 * t322;
t392 = t317 * t318;
t445 = t313 * t386 - t392;
t143 = t314 * t369 + (t445 * qJD(3) + (-t313 * t392 + t386) * qJD(2)) * t312;
t400 = t311 * t317;
t186 = t312 * t333 + t314 * t400;
t398 = t312 * t322;
t242 = -t311 * t398 + t313 * t314;
t147 = t186 * t320 + t242 * t316;
t385 = qJD(2) * t312;
t372 = t318 * t385;
t356 = t311 * t372;
t81 = qJD(4) * t147 + t143 * t316 - t320 * t356;
t82 = t316 * t356 - t186 * t381 + (qJD(4) * t242 + t143) * t320;
t24 = Ifges(5,5) * t82 - Ifges(5,6) * t81 + Ifges(5,3) * t142;
t27 = Ifges(6,1) * t142 - Ifges(6,4) * t82 + Ifges(6,5) * t81;
t448 = t24 + t27;
t374 = t316 * t400;
t380 = qJD(4) * t320;
t197 = qJD(4) * t374 - t313 * t380 - t320 * t369;
t244 = t313 * t316 + t320 * t400;
t198 = qJD(4) * t244 + t316 * t369;
t129 = Ifges(6,1) * t370 + Ifges(6,4) * t197 + Ifges(6,5) * t198;
t130 = -Ifges(5,5) * t197 - Ifges(5,6) * t198 + Ifges(5,3) * t370;
t447 = -t129 - t130;
t275 = -t320 * mrSges(5,1) + t316 * mrSges(5,2);
t446 = -m(5) * pkin(3) + t275;
t299 = pkin(10) * t400;
t396 = t313 * t321;
t245 = pkin(2) * t396 - t299;
t304 = t318 * t417;
t248 = pkin(9) * t398 + t304;
t180 = (t311 * t314 + t313 * t398) * pkin(10) + t248;
t416 = pkin(10) * t311;
t218 = (-pkin(2) * t322 - t318 * t416 - pkin(1)) * t312;
t99 = -t317 * t180 + t321 * (t218 * t311 + t401);
t397 = t313 * t317;
t399 = t311 * t321;
t247 = pkin(2) * t397 + pkin(10) * t399;
t227 = pkin(11) * t313 + t247;
t228 = (-pkin(3) * t321 - pkin(11) * t317 - pkin(2)) * t311;
t384 = qJD(3) * t311;
t236 = (pkin(3) * t317 - pkin(11) * t321) * t384;
t237 = t245 * qJD(3);
t330 = t227 * t381 - t228 * t380 - t316 * t236 - t320 * t237;
t96 = -t311 * (qJ(5) * t383 - qJD(5) * t321) + t330;
t402 = qJ(5) * t316;
t435 = pkin(4) + pkin(12);
t250 = -t320 * t435 - pkin(3) - t402;
t434 = pkin(5) + pkin(11);
t287 = t434 * t316;
t183 = -t250 * t315 + t287 * t319;
t184 = t250 * t319 + t287 * t315;
t415 = mrSges(7,3) * t320;
t266 = mrSges(7,1) * t316 + t315 * t415;
t267 = -mrSges(7,2) * t316 - t319 * t415;
t444 = m(7) * (-t183 * t315 + t184 * t319) + t319 * t267 - t315 * t266;
t443 = 2 * m(4);
t442 = 0.2e1 * m(5);
t441 = 2 * m(6);
t440 = 0.2e1 * m(7);
t439 = -2 * mrSges(3,3);
t438 = -2 * mrSges(4,3);
t243 = -t320 * t313 + t374;
t199 = t243 * t319 + t315 * t399;
t121 = qJD(6) * t199 + t198 * t315 + t319 * t370;
t334 = -t243 * t315 + t319 * t399;
t122 = qJD(6) * t334 + t198 * t319 - t315 * t370;
t57 = Ifges(7,1) * t121 + Ifges(7,4) * t122 - Ifges(7,5) * t197;
t436 = t57 / 0.2e1;
t146 = t186 * t316 - t242 * t320;
t185 = -t312 * t445 - t314 * t399;
t103 = t146 * t319 - t185 * t315;
t433 = t103 / 0.2e1;
t410 = Ifges(7,4) * t315;
t284 = Ifges(7,1) * t319 - t410;
t409 = Ifges(7,4) * t319;
t347 = Ifges(7,1) * t315 + t409;
t168 = -t284 * t379 + (Ifges(7,5) * t320 + t316 * t347) * qJD(4);
t432 = t168 / 0.2e1;
t431 = t199 / 0.2e1;
t430 = -t334 / 0.2e1;
t346 = Ifges(7,2) * t319 + t410;
t230 = Ifges(7,6) * t316 - t320 * t346;
t429 = t230 / 0.2e1;
t231 = Ifges(7,5) * t316 - t320 * t347;
t428 = t231 / 0.2e1;
t345 = -Ifges(7,5) * t315 - Ifges(7,6) * t319;
t427 = t345 * qJD(6) / 0.2e1;
t261 = t346 * qJD(6);
t426 = -t261 / 0.2e1;
t264 = t347 * qJD(6);
t425 = -t264 / 0.2e1;
t279 = Ifges(7,5) * t319 - Ifges(7,6) * t315;
t424 = t279 / 0.2e1;
t281 = -Ifges(7,2) * t315 + t409;
t423 = t281 / 0.2e1;
t422 = t284 / 0.2e1;
t421 = t313 / 0.2e1;
t420 = -t315 / 0.2e1;
t419 = t315 / 0.2e1;
t418 = -t319 / 0.2e1;
t136 = -t194 * t311 + t313 * t218;
t87 = pkin(3) * t185 - pkin(11) * t186 + t136;
t100 = t321 * t180 + t194 * t397 + t218 * t400;
t93 = pkin(11) * t242 + t100;
t39 = t316 * t87 + t320 * t93;
t414 = Ifges(4,4) * t317;
t413 = Ifges(4,4) * t321;
t412 = Ifges(5,4) * t316;
t411 = Ifges(5,4) * t320;
t408 = Ifges(6,6) * t316;
t407 = Ifges(6,6) * t320;
t196 = (t322 * t353 - t304) * qJD(2);
t222 = (pkin(2) * t318 - t322 * t416) * t385;
t144 = -t196 * t311 + t313 * t222;
t406 = t144 * mrSges(4,1);
t405 = t144 * mrSges(4,2);
t239 = -pkin(9) * t372 + t294;
t404 = t239 * mrSges(3,2);
t240 = t248 * qJD(2);
t403 = t240 * mrSges(3,1);
t394 = t315 * t284;
t393 = t315 * t435;
t388 = t319 * t281;
t387 = t319 * t435;
t155 = t320 * t227 + t316 * t228;
t104 = t146 * t315 + t185 * t319;
t31 = -qJD(6) * t104 - t142 * t315 + t319 * t81;
t32 = qJD(6) * t103 + t142 * t319 + t315 * t81;
t6 = Ifges(7,5) * t32 + Ifges(7,6) * t31 + Ifges(7,3) * t82;
t23 = Ifges(6,5) * t142 - Ifges(6,6) * t82 + Ifges(6,3) * t81;
t26 = Ifges(5,4) * t82 - Ifges(5,2) * t81 + Ifges(5,6) * t142;
t377 = t23 / 0.2e1 - t26 / 0.2e1;
t65 = Ifges(6,5) * t185 - Ifges(6,6) * t147 + Ifges(6,3) * t146;
t68 = Ifges(5,4) * t147 - Ifges(5,2) * t146 + Ifges(5,6) * t185;
t376 = -t68 / 0.2e1 + t65 / 0.2e1;
t375 = m(6) * pkin(11) + mrSges(6,1);
t288 = t434 * t320;
t55 = Ifges(7,5) * t121 + Ifges(7,6) * t122 - Ifges(7,3) * t197;
t78 = Ifges(4,5) * t143 - Ifges(4,6) * t142 + Ifges(4,3) * t356;
t127 = Ifges(6,5) * t370 + Ifges(6,6) * t197 + Ifges(6,3) * t198;
t131 = -Ifges(5,4) * t197 - Ifges(5,2) * t198 + Ifges(5,6) * t370;
t364 = t127 / 0.2e1 - t131 / 0.2e1;
t156 = -Ifges(6,5) * t399 - Ifges(6,6) * t244 + Ifges(6,3) * t243;
t160 = Ifges(5,4) * t244 - Ifges(5,2) * t243 - Ifges(5,6) * t399;
t363 = t156 / 0.2e1 - t160 / 0.2e1;
t257 = (Ifges(6,3) * t316 - t407) * qJD(4);
t263 = (-Ifges(5,2) * t316 + t411) * qJD(4);
t362 = t257 / 0.2e1 - t263 / 0.2e1;
t308 = Ifges(6,5) * t381;
t309 = Ifges(5,5) * t380;
t361 = -Ifges(5,6) * t381 / 0.2e1 + t309 / 0.2e1 + t380 * t450 + t308 / 0.2e1;
t360 = Ifges(5,6) * t320 / 0.2e1 + Ifges(6,5) * t449 + (Ifges(5,5) / 0.2e1 + t450) * t316;
t277 = -Ifges(6,3) * t320 - t408;
t283 = Ifges(5,2) * t320 + t412;
t359 = -t283 / 0.2e1 + t277 / 0.2e1;
t51 = t82 * mrSges(6,1) + t142 * mrSges(6,2);
t38 = -t316 * t93 + t320 * t87;
t163 = -t197 * mrSges(6,1) + mrSges(6,2) * t370;
t154 = -t316 * t227 + t228 * t320;
t358 = pkin(4) * t381 - qJD(5) * t316;
t357 = -t180 * t382 - t218 * t370 - t451 * t317;
t25 = Ifges(6,4) * t142 - Ifges(6,2) * t82 + Ifges(6,6) * t81;
t28 = Ifges(5,1) * t82 - Ifges(5,4) * t81 + Ifges(5,5) * t142;
t355 = t28 / 0.2e1 + t6 / 0.2e1 - t25 / 0.2e1;
t40 = Ifges(7,5) * t104 + Ifges(7,6) * t103 + Ifges(7,3) * t147;
t67 = Ifges(6,4) * t185 - Ifges(6,2) * t147 + Ifges(6,6) * t146;
t70 = Ifges(5,1) * t147 - Ifges(5,4) * t146 + Ifges(5,5) * t185;
t354 = t67 / 0.2e1 - t70 / 0.2e1 - t40 / 0.2e1;
t34 = -qJ(5) * t185 - t39;
t128 = Ifges(6,4) * t370 + Ifges(6,2) * t197 + Ifges(6,6) * t198;
t132 = -Ifges(5,1) * t197 - Ifges(5,4) * t198 + Ifges(5,5) * t370;
t352 = t132 / 0.2e1 + t55 / 0.2e1 - t128 / 0.2e1;
t149 = pkin(4) * t399 - t154;
t113 = -Ifges(7,5) * t334 + Ifges(7,6) * t199 + Ifges(7,3) * t244;
t157 = -Ifges(6,4) * t399 - Ifges(6,2) * t244 + Ifges(6,6) * t243;
t161 = Ifges(5,1) * t244 - Ifges(5,4) * t243 - Ifges(5,5) * t399;
t351 = -t157 / 0.2e1 + t161 / 0.2e1 + t113 / 0.2e1;
t166 = t328 * Ifges(7,5) + t329 * Ifges(7,6) + Ifges(7,3) * t380;
t258 = (-Ifges(6,2) * t320 + t408) * qJD(4);
t265 = (Ifges(5,1) * t320 - t412) * qJD(4);
t350 = -t258 / 0.2e1 + t265 / 0.2e1 + t166 / 0.2e1;
t229 = Ifges(7,3) * t316 + t320 * t345;
t278 = -Ifges(6,2) * t316 - t407;
t285 = Ifges(5,1) * t316 + t411;
t349 = t285 / 0.2e1 - t278 / 0.2e1 + t229 / 0.2e1;
t348 = mrSges(7,1) * t319 - mrSges(7,2) * t315;
t276 = mrSges(7,1) * t315 + mrSges(7,2) * t319;
t274 = t320 * mrSges(6,2) - t316 * mrSges(6,3);
t344 = -pkin(4) * t320 - t402;
t21 = pkin(5) * t147 - t185 * t435 - t38;
t92 = -pkin(3) * t242 - t99;
t326 = -qJ(5) * t147 + t92;
t33 = t146 * t435 + t326;
t11 = t21 * t319 - t315 * t33;
t12 = t21 * t315 + t319 * t33;
t343 = t11 * t315 - t12 * t319;
t116 = pkin(5) * t244 + pkin(12) * t399 + t149;
t226 = t299 + (-pkin(2) * t321 - pkin(3)) * t313;
t327 = -qJ(5) * t244 + t226;
t125 = t243 * t435 + t327;
t60 = t116 * t319 - t125 * t315;
t61 = t116 * t315 + t125 * t319;
t342 = t315 * t60 - t319 * t61;
t63 = -mrSges(7,2) * t147 + mrSges(7,3) * t103;
t64 = mrSges(7,1) * t147 - mrSges(7,3) * t104;
t341 = -t315 * t64 + t319 * t63;
t48 = -t180 * t383 + t196 * t397 + t218 * t369 + t222 * t400 + t451 * t321;
t46 = pkin(11) * t356 + t48;
t59 = pkin(3) * t142 - pkin(11) * t143 + t144;
t14 = -t316 * t46 + t320 * t59 - t93 * t380 - t87 * t381;
t217 = (pkin(12) * t316 - qJ(5) * t320) * qJD(4) + t358;
t269 = qJD(4) * t288;
t119 = qJD(6) * t183 + t217 * t319 + t269 * t315;
t120 = -qJD(6) * t184 - t217 * t315 + t269 * t319;
t340 = -t119 * t315 - t120 * t319;
t152 = -mrSges(7,2) * t244 + mrSges(7,3) * t199;
t153 = mrSges(7,1) * t244 + mrSges(7,3) * t334;
t339 = t319 * t152 - t315 * t153;
t148 = qJ(5) * t399 - t155;
t106 = -t227 * t380 - t228 * t381 + t236 * t320 - t316 * t237;
t41 = Ifges(7,4) * t104 + Ifges(7,2) * t103 + Ifges(7,6) * t147;
t42 = Ifges(7,1) * t104 + Ifges(7,4) * t103 + Ifges(7,5) * t147;
t335 = t41 * t418 + t42 * t420;
t13 = t316 * t59 + t320 * t46 + t87 * t380 - t381 * t93;
t114 = -Ifges(7,4) * t334 + Ifges(7,2) * t199 + Ifges(7,6) * t244;
t115 = -Ifges(7,1) * t334 + Ifges(7,4) * t199 + Ifges(7,5) * t244;
t332 = t114 * t418 + t115 * t420;
t238 = t247 * qJD(3);
t5 = -qJ(5) * t142 - qJD(5) * t185 - t13;
t325 = qJ(5) * t197 - qJD(5) * t244 + t238;
t47 = -t196 * t396 + (-pkin(3) * t372 - t222 * t321) * t311 - t357;
t324 = -qJ(5) * t82 - qJD(5) * t147 + t47;
t292 = Ifges(3,5) * t322 * t385;
t291 = Ifges(4,5) * t369;
t272 = -pkin(3) + t344;
t268 = t434 * t381;
t256 = (mrSges(5,1) * t316 + mrSges(5,2) * t320) * qJD(4);
t255 = (-mrSges(6,2) * t316 - mrSges(6,3) * t320) * qJD(4);
t254 = t348 * qJD(6);
t253 = -mrSges(4,2) * t313 + mrSges(4,3) * t399;
t252 = mrSges(4,1) * t313 - mrSges(4,3) * t400;
t249 = t348 * t320;
t246 = -pkin(9) * t312 * t318 + t305;
t241 = -qJ(5) * t380 + t358;
t235 = (Ifges(4,1) * t321 - t414) * t384;
t234 = (-Ifges(4,2) * t317 + t413) * t384;
t233 = -Ifges(4,6) * t370 + t291;
t232 = (mrSges(4,1) * t317 + mrSges(4,2) * t321) * t384;
t224 = Ifges(4,5) * t313 + (Ifges(4,1) * t317 + t413) * t311;
t223 = Ifges(4,6) * t313 + (Ifges(4,2) * t321 + t414) * t311;
t213 = mrSges(7,1) * t380 - mrSges(7,3) * t328;
t212 = -mrSges(7,2) * t380 + mrSges(7,3) * t329;
t207 = -mrSges(5,1) * t399 - mrSges(5,3) * t244;
t206 = mrSges(5,2) * t399 - mrSges(5,3) * t243;
t205 = mrSges(6,1) * t244 - mrSges(6,2) * t399;
t204 = mrSges(6,1) * t243 + mrSges(6,3) * t399;
t179 = -mrSges(7,1) * t329 + mrSges(7,2) * t328;
t173 = -mrSges(6,2) * t243 - mrSges(6,3) * t244;
t172 = mrSges(5,1) * t243 + mrSges(5,2) * t244;
t167 = -t281 * t379 + (Ifges(7,6) * t320 + t316 * t346) * qJD(4);
t165 = -mrSges(5,2) * t370 - mrSges(5,3) * t198;
t164 = mrSges(5,1) * t370 + mrSges(5,3) * t197;
t162 = mrSges(6,1) * t198 - mrSges(6,3) * t370;
t159 = Ifges(5,5) * t244 - Ifges(5,6) * t243 - Ifges(5,3) * t399;
t158 = -Ifges(6,1) * t399 - Ifges(6,4) * t244 + Ifges(6,5) * t243;
t151 = mrSges(4,1) * t242 - mrSges(4,3) * t186;
t150 = -mrSges(4,2) * t242 - mrSges(4,3) * t185;
t145 = pkin(4) * t243 + t327;
t135 = -mrSges(7,1) * t199 - mrSges(7,2) * t334;
t134 = mrSges(5,1) * t198 - mrSges(5,2) * t197;
t133 = -mrSges(6,2) * t198 + mrSges(6,3) * t197;
t126 = -pkin(5) * t243 - t148;
t124 = mrSges(4,1) * t356 - mrSges(4,3) * t143;
t123 = -mrSges(4,2) * t356 - mrSges(4,3) * t142;
t112 = Ifges(4,1) * t186 - Ifges(4,4) * t185 + Ifges(4,5) * t242;
t111 = Ifges(4,4) * t186 - Ifges(4,2) * t185 + Ifges(4,6) * t242;
t110 = mrSges(5,1) * t185 - mrSges(5,3) * t147;
t109 = -mrSges(5,2) * t185 - mrSges(5,3) * t146;
t108 = mrSges(6,1) * t147 + mrSges(6,2) * t185;
t107 = mrSges(6,1) * t146 - mrSges(6,3) * t185;
t102 = pkin(4) * t198 + t325;
t101 = -pkin(4) * t370 - t106;
t98 = mrSges(7,2) * t197 + mrSges(7,3) * t122;
t97 = -mrSges(7,1) * t197 - mrSges(7,3) * t121;
t95 = -mrSges(6,2) * t146 - mrSges(6,3) * t147;
t94 = mrSges(5,1) * t146 + mrSges(5,2) * t147;
t91 = mrSges(4,1) * t142 + mrSges(4,2) * t143;
t86 = t198 * t435 + t325;
t80 = Ifges(4,1) * t143 - Ifges(4,4) * t142 + Ifges(4,5) * t356;
t79 = Ifges(4,4) * t143 - Ifges(4,2) * t142 + Ifges(4,6) * t356;
t72 = -pkin(5) * t198 - t96;
t71 = -pkin(5) * t197 - t370 * t435 - t106;
t69 = Ifges(6,1) * t185 - Ifges(6,4) * t147 + Ifges(6,5) * t146;
t66 = Ifges(5,5) * t147 - Ifges(5,6) * t146 + Ifges(5,3) * t185;
t62 = -mrSges(7,1) * t122 + mrSges(7,2) * t121;
t56 = Ifges(7,4) * t121 + Ifges(7,2) * t122 - Ifges(7,6) * t197;
t54 = -mrSges(7,1) * t103 + mrSges(7,2) * t104;
t53 = mrSges(5,1) * t142 - mrSges(5,3) * t82;
t52 = -mrSges(5,2) * t142 - mrSges(5,3) * t81;
t50 = mrSges(6,1) * t81 - mrSges(6,3) * t142;
t49 = (t196 * t313 + t222 * t311) * t321 + t357;
t43 = pkin(4) * t146 + t326;
t37 = -mrSges(6,2) * t81 - mrSges(6,3) * t82;
t36 = mrSges(5,1) * t81 + mrSges(5,2) * t82;
t35 = -pkin(4) * t185 - t38;
t22 = -pkin(5) * t146 - t34;
t20 = mrSges(7,1) * t82 - mrSges(7,3) * t32;
t19 = -mrSges(7,2) * t82 + mrSges(7,3) * t31;
t18 = -qJD(6) * t61 - t315 * t86 + t319 * t71;
t17 = qJD(6) * t60 + t315 * t71 + t319 * t86;
t16 = pkin(4) * t81 + t324;
t15 = -mrSges(7,1) * t31 + mrSges(7,2) * t32;
t10 = t435 * t81 + t324;
t9 = -pkin(4) * t142 - t14;
t8 = Ifges(7,1) * t32 + Ifges(7,4) * t31 + Ifges(7,5) * t82;
t7 = Ifges(7,4) * t32 + Ifges(7,2) * t31 + Ifges(7,6) * t82;
t4 = -pkin(5) * t81 - t5;
t3 = pkin(5) * t82 - t142 * t435 - t14;
t2 = -qJD(6) * t12 - t10 * t315 + t3 * t319;
t1 = qJD(6) * t11 + t10 * t319 + t3 * t315;
t29 = [0.2e1 * m(3) * (t239 * t248 - t240 * t246) + (t70 + t40 - t67) * t82 + (t65 - t68) * t81 + (t1 * t12 + t11 * t2 + t22 * t4) * t440 + (t16 * t43 + t34 * t5 + t35 * t9) * t441 + (t13 * t39 + t14 * t38 + t47 * t92) * t442 + (t100 * t48 + t136 * t144 + t49 * t99) * t443 + (t28 + t6 - t25) * t147 + (t23 - t26) * t146 + (t66 + t69 - t111) * t142 + (t292 - 0.2e1 * t403 - 0.2e1 * t404) * t314 + t242 * t78 + 0.2e1 * t48 * t150 + 0.2e1 * t49 * t151 + t143 * t112 + 0.2e1 * t136 * t91 + 0.2e1 * t100 * t123 + 0.2e1 * t99 * t124 + t104 * t8 + 0.2e1 * t5 * t107 + 0.2e1 * t9 * t108 + 0.2e1 * t13 * t109 + 0.2e1 * t14 * t110 + t103 * t7 + 0.2e1 * t92 * t36 + 0.2e1 * t47 * t94 + 0.2e1 * t16 * t95 + 0.2e1 * t1 * t63 + 0.2e1 * t2 * t64 + 0.2e1 * t38 * t53 + 0.2e1 * t4 * t54 + 0.2e1 * t34 * t50 + 0.2e1 * t35 * t51 + 0.2e1 * t39 * t52 + t31 * t41 + t32 * t42 + 0.2e1 * t43 * t37 + 0.2e1 * t22 * t15 + 0.2e1 * t12 * t19 + 0.2e1 * t11 * t20 + (t80 + 0.2e1 * t405) * t186 + (-t79 + 0.2e1 * t406 + t448) * t185 + (0.2e1 * (t239 * t322 + t240 * t318) * mrSges(3,3) + ((t246 * t439 + Ifges(3,5) * t314 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t322) * t452) * t322 + (t311 * (Ifges(4,5) * t186 - Ifges(4,6) * t185 + Ifges(4,3) * t242) + t248 * t439 - 0.2e1 * t314 * Ifges(3,6) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t318 + (Ifges(3,1) - Ifges(3,2)) * t322) * t452) * t318) * qJD(2)) * t312; t292 + (-t223 / 0.2e1 + t158 / 0.2e1 + t159 / 0.2e1) * t142 + (-t234 / 0.2e1 + t130 / 0.2e1 + t129 / 0.2e1) * t185 + (t94 - t151) * t238 + t78 * t421 + t8 * t430 + t7 * t431 + t56 * t433 + t104 * t436 + m(5) * (t106 * t38 + t13 * t155 + t14 * t154 + t226 * t47 + t238 * t92 - t330 * t39) - t330 * t109 + m(4) * (t100 * t237 - t238 * t99 + t245 * t49 + t247 * t48) - t403 - t404 + t245 * t124 + t247 * t123 + t49 * t252 + t48 * t253 + t242 * t233 / 0.2e1 + t237 * t150 + t143 * t224 / 0.2e1 + t226 * t36 + t136 * t232 + t186 * t235 / 0.2e1 + t5 * t204 + t9 * t205 + t13 * t206 + t14 * t207 + t39 * t165 + t47 * t172 + t16 * t173 + t1 * t152 + t2 * t153 + t154 * t53 + t155 * t52 + t34 * t162 + t35 * t163 + t38 * t164 + t148 * t50 + t149 * t51 + t145 * t37 + t43 * t133 + t92 * t134 + t4 * t135 + t121 * t42 / 0.2e1 + t122 * t41 / 0.2e1 + t126 * t15 + t31 * t114 / 0.2e1 + t32 * t115 / 0.2e1 + t96 * t107 + t101 * t108 + t106 * t110 + t102 * t95 + t11 * t97 + t12 * t98 + t72 * t54 + t17 * t63 + t18 * t64 + t60 * t20 + t61 * t19 + t22 * t62 + m(7) * (t1 * t61 + t11 * t18 + t12 * t17 + t126 * t4 + t2 * t60 + t22 * t72) + m(6) * (t101 * t35 + t102 * t43 + t145 * t16 + t148 * t5 + t149 * t9 + t34 * t96) + ((t405 + t80 / 0.2e1) * t317 + (-t406 - t24 / 0.2e1 - t27 / 0.2e1 + t79 / 0.2e1) * t321 + (Ifges(4,3) * t421 + (Ifges(4,5) * t317 + Ifges(4,6) * t321) * t311 / 0.2e1) * t372 + (-m(4) * t144 - t91) * pkin(2) + ((-t99 * mrSges(4,3) + t112 / 0.2e1) * t321 + (-t111 / 0.2e1 + t66 / 0.2e1 + t69 / 0.2e1 - t100 * mrSges(4,3)) * t317) * qJD(3)) * t311 + t351 * t82 + t352 * t147 + t354 * t197 + t355 * t244 + t363 * t81 + t364 * t146 - Ifges(3,6) * t372 + t376 * t198 + t377 * t243; 0.2e1 * (-t252 + t172) * t238 + (t126 * t72 + t17 * t61 + t18 * t60) * t440 + (t101 * t149 + t102 * t145 + t148 * t96) * t441 + (t237 * t247 - t238 * t245) * t443 + (t106 * t154 - t155 * t330 + t226 * t238) * t442 - 0.2e1 * t330 * t206 - t334 * t57 + t313 * t233 + 0.2e1 * t237 * t253 + 0.2e1 * t226 * t134 + 0.2e1 * t96 * t204 + 0.2e1 * t101 * t205 + 0.2e1 * t106 * t207 + t199 * t56 + 0.2e1 * t155 * t165 + 0.2e1 * t102 * t173 + 0.2e1 * t17 * t152 + 0.2e1 * t18 * t153 + 0.2e1 * t148 * t162 + 0.2e1 * t149 * t163 + 0.2e1 * t154 * t164 + 0.2e1 * t145 * t133 + 0.2e1 * t72 * t135 + t121 * t115 + t122 * t114 + 0.2e1 * t126 * t62 + 0.2e1 * t60 * t97 + 0.2e1 * t61 * t98 + (t132 + t55 - t128) * t244 + (t127 - t131) * t243 + (-t160 + t156) * t198 + (t157 - t161 - t113) * t197 + (-0.2e1 * pkin(2) * t232 + t317 * t235 + (t234 + t447) * t321 + ((t245 * t438 + t224) * t321 + (t247 * t438 + t158 + t159 - t223) * t317) * qJD(3)) * t311; t78 + m(7) * (t1 * t184 + t11 * t120 + t119 * t12 + t183 * t2 - t22 * t268 + t288 * t4) + t446 * t47 + t32 * t428 + t31 * t429 + t104 * t432 + t167 * t433 + m(6) * (t16 * t272 + t241 * t43) + t288 * t15 + t2 * t266 + t1 * t267 - t268 * t54 + t272 * t37 + t16 * t274 + t43 * t255 + t92 * t256 + t4 * t249 + t241 * t95 + t12 * t212 + t11 * t213 + t183 * t20 + t184 * t19 + t22 * t179 + t119 * t63 + t120 * t64 - t48 * mrSges(4,2) + t49 * mrSges(4,1) - pkin(3) * t36 + (t7 * t418 - t5 * mrSges(6,1) + t13 * mrSges(5,3) + t8 * t420 + (t41 * t419 + t418 * t42) * qJD(6) + (t35 * mrSges(6,1) - t38 * mrSges(5,3) - t354) * qJD(4) + (-t50 + t52 + (t108 - t110) * qJD(4) + m(5) * (-qJD(4) * t38 + t13) + m(6) * (qJD(4) * t35 - t5)) * pkin(11) - t377) * t320 + t349 * t82 + t350 * t147 + t359 * t81 + t360 * t142 + t361 * t185 + t362 * t146 + (t9 * mrSges(6,1) - t14 * mrSges(5,3) + (t34 * mrSges(6,1) - t39 * mrSges(5,3) - t335 + t376) * qJD(4) + (t51 - t53 + (t107 - t109) * qJD(4) + m(5) * (-qJD(4) * t39 - t14) + m(6) * (qJD(4) * t34 + t9)) * pkin(11) + t355) * t316; t291 + m(7) * (t119 * t61 + t120 * t60 - t126 * t268 + t17 * t184 + t18 * t183 + t288 * t72) + m(6) * (t102 * t272 + t145 * t241) + (-mrSges(4,1) + t446) * t238 + t121 * t428 + t122 * t429 + t168 * t430 + t167 * t431 + (t56 * t418 - t96 * mrSges(6,1) - t330 * mrSges(5,3) + t57 * t420 + (t114 * t419 + t115 * t418) * qJD(6) + (t149 * mrSges(6,1) - t154 * mrSges(5,3) + t351) * qJD(4) + (-t162 + t165 + (t205 - t207) * qJD(4) + m(5) * (-qJD(4) * t154 - t330) + m(6) * (qJD(4) * t149 - t96)) * pkin(11) - t364) * t320 + t288 * t62 + t18 * t266 + t17 * t267 - t268 * t135 + t272 * t133 + t102 * t274 + t145 * t255 + t226 * t256 + t72 * t249 + t241 * t173 - t237 * mrSges(4,2) + t61 * t212 + t60 * t213 + t183 * t97 + t184 * t98 + t126 * t179 + t119 * t152 + t120 * t153 - pkin(3) * t134 - t349 * t197 + t350 * t244 + t359 * t198 + t362 * t243 + (t101 * mrSges(6,1) - t106 * mrSges(5,3) + (t148 * mrSges(6,1) - t155 * mrSges(5,3) - t332 + t363) * qJD(4) + (t163 - t164 + (t204 - t206) * qJD(4) + m(5) * (-qJD(4) * t155 - t106) + m(6) * (qJD(4) * t148 + t101)) * pkin(11) + t352) * t316 + (-t361 * t321 + (-Ifges(4,6) + t360) * t383) * t311; 0.2e1 * t288 * t179 + 0.2e1 * t119 * t267 - 0.2e1 * t268 * t249 + 0.2e1 * t272 * t255 - 0.2e1 * pkin(3) * t256 + 0.2e1 * t120 * t266 + 0.2e1 * t184 * t212 + 0.2e1 * t183 * t213 + (t119 * t184 + t120 * t183 - t268 * t288) * t440 + 0.2e1 * (m(6) * t272 + t274) * t241 + (t166 - t258 + t265 + (t230 * t319 + t231 * t315 + t277 - t283) * qJD(4)) * t316 + (-t319 * t167 - t315 * t168 - t257 + t263 + (t230 * t315 - t231 * t319) * qJD(6) + (t229 - t278 + t285) * qJD(4)) * t320; t448 + (t8 / 0.2e1 - t435 * t20 - t2 * mrSges(7,3)) * t319 + (-t7 / 0.2e1 - t435 * t19 - t1 * mrSges(7,3)) * t315 + (t343 * mrSges(7,3) - (-m(7) * t343 + t341) * t435 + t335) * qJD(6) + t32 * t422 + t31 * t423 + t82 * t424 + t104 * t425 + t103 * t426 + t147 * t427 + (-t107 + t54) * qJD(5) + (-t50 + t15) * qJ(5) + m(6) * (-pkin(4) * t9 - qJ(5) * t5 - qJD(5) * t34) + m(7) * (qJ(5) * t4 + qJD(5) * t22 - t1 * t393 - t2 * t387) + t4 * t276 + t22 * t254 - pkin(4) * t51 + t14 * mrSges(5,1) - t13 * mrSges(5,2) + t9 * mrSges(6,2) - t5 * mrSges(6,3); (-t18 * mrSges(7,3) - t435 * t97 + t436) * t319 + (-t56 / 0.2e1 - t435 * t98 - t17 * mrSges(7,3)) * t315 + (t342 * mrSges(7,3) - (-m(7) * t342 + t339) * t435 + t332) * qJD(6) + (-t204 + t135) * qJD(5) + (-t162 + t62) * qJ(5) + m(6) * (-pkin(4) * t101 - qJ(5) * t96 - qJD(5) * t148) + t121 * t422 + t122 * t423 - t334 * t425 + t199 * t426 + t244 * t427 + m(7) * (qJ(5) * t72 + qJD(5) * t126 - t17 * t393 - t18 * t387) - t447 + t72 * t276 - t197 * t279 / 0.2e1 + t126 * t254 - pkin(4) * t163 + t330 * mrSges(5,2) + t106 * mrSges(5,1) + t101 * mrSges(6,2) - t96 * mrSges(6,3); t308 + t309 + t319 * t432 + t167 * t420 + t316 * t427 + t288 * t254 - t268 * t276 + qJD(5) * t249 + qJ(5) * t179 - t213 * t387 + m(7) * (-qJ(5) * t268 + qJD(5) * t288 - t119 * t393 - t120 * t387) - t212 * t393 + t340 * mrSges(7,3) + (qJD(5) * t375 - t261 * t418 - t264 * t420) * t320 + ((-t230 / 0.2e1 - t184 * mrSges(7,3) + t284 * t449) * t319 + (t183 * mrSges(7,3) + t320 * t423 - t231 / 0.2e1) * t315 - t444 * t435) * qJD(6) + ((-pkin(4) * mrSges(6,1) - Ifges(6,4) + t424) * t320 + (t394 / 0.2e1 - qJ(5) * mrSges(6,1) + t388 / 0.2e1 - Ifges(5,6)) * t316 + (m(6) * t344 + t274 + t275) * pkin(11)) * qJD(4); 0.2e1 * qJ(5) * t254 + t261 * t315 - t264 * t319 + (-t388 - t394) * qJD(6) + 0.2e1 * (mrSges(6,3) + t276 + (m(6) + m(7)) * qJ(5)) * qJD(5); t315 * t19 + t319 * t20 + t341 * qJD(6) + m(7) * (-qJD(6) * t343 + t1 * t315 + t2 * t319) + m(6) * t9 + t51; t315 * t98 + t319 * t97 + t339 * qJD(6) + m(7) * (-qJD(6) * t342 + t17 * t315 + t18 * t319) + m(6) * t101 + t163; -m(7) * t340 + qJD(6) * t444 + t315 * t212 + t319 * t213 + t375 * t380; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t18 - mrSges(7,2) * t17 + t55; mrSges(7,1) * t120 - mrSges(7,2) * t119 + t166; ((mrSges(7,2) * t435 - Ifges(7,6)) * t319 + (mrSges(7,1) * t435 - Ifges(7,5)) * t315) * qJD(6); -t276 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t29(1) t29(2) t29(4) t29(7) t29(11) t29(16); t29(2) t29(3) t29(5) t29(8) t29(12) t29(17); t29(4) t29(5) t29(6) t29(9) t29(13) t29(18); t29(7) t29(8) t29(9) t29(10) t29(14) t29(19); t29(11) t29(12) t29(13) t29(14) t29(15) t29(20); t29(16) t29(17) t29(18) t29(19) t29(20) t29(21);];
Mq  = res;
