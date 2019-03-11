% Calculate time derivative of joint inertia matrix for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:52:14
% EndTime: 2019-03-09 19:52:38
% DurationCPUTime: 11.57s
% Computational Cost: add. (26780->911), mult. (79785->1334), div. (0->0), fcn. (82669->14), ass. (0->373)
t330 = sin(pkin(6));
t466 = 0.2e1 * t330;
t328 = sin(pkin(13));
t331 = cos(pkin(13));
t335 = sin(qJ(5));
t421 = cos(qJ(5));
t343 = -t335 * t328 + t331 * t421;
t282 = t343 * qJD(5);
t290 = t328 * t421 + t335 * t331;
t338 = cos(qJ(6));
t334 = sin(qJ(6));
t384 = qJD(6) * t334;
t341 = -t282 * t338 + t290 * t384;
t427 = t331 / 0.2e1;
t424 = t334 / 0.2e1;
t422 = t338 / 0.2e1;
t332 = cos(pkin(7));
t336 = sin(qJ(3));
t398 = t332 * t336;
t329 = sin(pkin(7));
t339 = cos(qJ(3));
t401 = t329 * t339;
t286 = pkin(2) * t398 + pkin(10) * t401;
t263 = qJ(4) * t332 + t286;
t264 = (-pkin(3) * t339 - qJ(4) * t336 - pkin(2)) * t329;
t194 = -t263 * t328 + t331 * t264;
t402 = t329 * t336;
t281 = t328 * t332 + t331 * t402;
t156 = -pkin(4) * t401 - pkin(11) * t281 + t194;
t195 = t331 * t263 + t328 * t264;
t279 = -t328 * t402 + t331 * t332;
t170 = pkin(11) * t279 + t195;
t465 = t335 * t156 + t421 * t170;
t340 = cos(qJ(2));
t393 = t339 * t340;
t337 = sin(qJ(2));
t396 = t336 * t337;
t464 = t332 * t393 - t396;
t333 = cos(pkin(6));
t420 = pkin(1) * t333;
t321 = t337 * t420;
t400 = t330 * t340;
t287 = pkin(9) * t400 + t321;
t224 = (t329 * t333 + t332 * t400) * pkin(10) + t287;
t322 = t340 * t420;
t361 = t330 * (-pkin(10) * t332 - pkin(9));
t354 = t337 * t361;
t235 = pkin(2) * t333 + t322 + t354;
t419 = pkin(10) * t329;
t256 = (-pkin(2) * t340 - t337 * t419 - pkin(1)) * t330;
t130 = -t336 * t224 + t339 * (t235 * t332 + t256 * t329);
t249 = (-qJD(4) * t336 + (pkin(3) * t336 - qJ(4) * t339) * qJD(3)) * t329;
t388 = qJD(3) * t339;
t373 = t332 * t388;
t389 = qJD(3) * t336;
t375 = t329 * t389;
t271 = pkin(2) * t373 - pkin(10) * t375;
t255 = qJD(4) * t332 + t271;
t181 = t331 * t249 - t255 * t328;
t390 = qJD(3) * t329;
t399 = t331 * t339;
t152 = (pkin(4) * t336 - pkin(11) * t399) * t390 + t181;
t182 = t328 * t249 + t331 * t255;
t374 = t329 * t388;
t363 = t328 * t374;
t169 = -pkin(11) * t363 + t182;
t54 = -qJD(5) * t465 + t152 * t421 - t335 * t169;
t177 = t333 * t374 + (t464 * qJD(3) + (-t332 * t396 + t393) * qJD(2)) * t330;
t391 = qJD(2) * t330;
t376 = t337 * t391;
t364 = t329 * t376;
t154 = t177 * t331 + t328 * t364;
t394 = t337 * t339;
t395 = t336 * t340;
t348 = t332 * t395 + t394;
t176 = t333 * t375 + (t348 * qJD(3) + (t332 * t394 + t395) * qJD(2)) * t330;
t280 = -t329 * t400 + t333 * t332;
t316 = qJD(2) * t322;
t239 = qJD(2) * t354 + t316;
t240 = (t340 * t361 - t321) * qJD(2);
t259 = (pkin(2) * t337 - t340 * t419) * t391;
t84 = -t224 * t389 + t235 * t373 + t339 * t239 + t240 * t398 + t256 * t374 + t259 * t402;
t74 = qJ(4) * t364 + qJD(4) * t280 + t84;
t178 = -t240 * t329 + t332 * t259;
t228 = t330 * t348 + t333 * t402;
t83 = pkin(3) * t176 - qJ(4) * t177 - qJD(4) * t228 + t178;
t34 = -t328 * t74 + t331 * t83;
t26 = pkin(4) * t176 - pkin(11) * t154 + t34;
t153 = -t177 * t328 + t331 * t364;
t35 = t328 * t83 + t331 * t74;
t29 = pkin(11) * t153 + t35;
t184 = t228 * t331 + t280 * t328;
t227 = -t330 * t464 - t333 * t401;
t171 = -t235 * t329 + t332 * t256;
t120 = pkin(3) * t227 - qJ(4) * t228 + t171;
t214 = t235 * t398;
t131 = t339 * t224 + t256 * t402 + t214;
t124 = qJ(4) * t280 + t131;
t70 = t331 * t120 - t124 * t328;
t47 = pkin(4) * t227 - pkin(11) * t184 + t70;
t183 = -t228 * t328 + t280 * t331;
t71 = t328 * t120 + t331 * t124;
t50 = pkin(11) * t183 + t71;
t417 = t335 * t47 + t421 * t50;
t6 = -qJD(5) * t417 + t26 * t421 - t335 * t29;
t463 = 2 * m(4);
t462 = 2 * m(5);
t461 = 2 * m(6);
t460 = 2 * m(7);
t459 = -2 * mrSges(3,3);
t458 = -2 * mrSges(4,3);
t457 = -2 * mrSges(6,3);
t418 = pkin(11) + qJ(4);
t303 = t418 * t331;
t366 = t421 * qJD(4);
t367 = t335 * t418;
t368 = qJD(5) * t421;
t382 = t335 * qJD(4);
t197 = t303 * t368 + t331 * t382 + (-qJD(5) * t367 + t366) * t328;
t456 = 0.2e1 * t197;
t362 = t418 * t421;
t242 = t303 * t335 + t328 * t362;
t455 = 0.2e1 * t242;
t345 = t183 * t421 - t335 * t184;
t63 = qJD(5) * t345 + t335 * t153 + t154 * t421;
t127 = t335 * t183 + t184 * t421;
t64 = qJD(5) * t127 - t153 * t421 + t335 * t154;
t24 = Ifges(6,1) * t63 - Ifges(6,4) * t64 + Ifges(6,5) * t176;
t454 = t24 / 0.2e1;
t69 = Ifges(6,1) * t127 + Ifges(6,4) * t345 + Ifges(6,5) * t227;
t453 = t69 / 0.2e1;
t76 = Ifges(5,4) * t154 + Ifges(5,2) * t153 + Ifges(5,6) * t176;
t452 = t76 / 0.2e1;
t77 = Ifges(5,1) * t154 + Ifges(5,4) * t153 + Ifges(5,5) * t176;
t451 = t77 / 0.2e1;
t344 = t279 * t421 - t335 * t281;
t166 = qJD(5) * t344 + t343 * t374;
t202 = t335 * t279 + t281 * t421;
t167 = qJD(5) * t202 + t290 * t374;
t102 = Ifges(6,1) * t166 - Ifges(6,4) * t167 + Ifges(6,5) * t375;
t450 = t102 / 0.2e1;
t283 = t290 * qJD(5);
t383 = qJD(6) * t338;
t342 = t282 * t334 + t290 * t383;
t135 = -Ifges(7,4) * t341 - Ifges(7,2) * t342 + Ifges(7,6) * t283;
t449 = t135 / 0.2e1;
t136 = -Ifges(7,1) * t341 - Ifges(7,4) * t342 + Ifges(7,5) * t283;
t448 = t136 / 0.2e1;
t139 = Ifges(6,1) * t202 + Ifges(6,4) * t344 - Ifges(6,5) * t401;
t447 = t139 / 0.2e1;
t446 = t153 / 0.2e1;
t445 = t154 / 0.2e1;
t185 = -t202 * t334 - t338 * t401;
t444 = t185 / 0.2e1;
t349 = -t202 * t338 + t334 * t401;
t443 = -t349 / 0.2e1;
t410 = Ifges(7,4) * t338;
t357 = -Ifges(7,2) * t334 + t410;
t190 = -Ifges(7,6) * t343 + t290 * t357;
t442 = t190 / 0.2e1;
t411 = Ifges(7,4) * t334;
t358 = Ifges(7,1) * t338 - t411;
t191 = -Ifges(7,5) * t343 + t290 * t358;
t441 = t191 / 0.2e1;
t211 = Ifges(6,1) * t282 - Ifges(6,4) * t283;
t440 = t211 / 0.2e1;
t412 = Ifges(5,4) * t331;
t222 = (Ifges(5,6) * t336 + (-Ifges(5,2) * t328 + t412) * t339) * t390;
t439 = t222 / 0.2e1;
t413 = Ifges(5,4) * t328;
t223 = (Ifges(5,5) * t336 + (Ifges(5,1) * t331 - t413) * t339) * t390;
t438 = t223 / 0.2e1;
t232 = Ifges(6,1) * t290 + Ifges(6,4) * t343;
t437 = t232 / 0.2e1;
t436 = -t290 / 0.2e1;
t325 = Ifges(7,5) * t383;
t435 = -Ifges(7,6) * t384 / 0.2e1 + t325 / 0.2e1;
t295 = t357 * qJD(6);
t434 = t295 / 0.2e1;
t296 = t358 * qJD(6);
t433 = t296 / 0.2e1;
t432 = Ifges(7,5) * t424 + Ifges(7,6) * t422;
t309 = Ifges(7,2) * t338 + t411;
t431 = -t309 / 0.2e1;
t430 = t309 / 0.2e1;
t310 = Ifges(7,1) * t334 + t410;
t429 = t310 / 0.2e1;
t428 = -t328 / 0.2e1;
t426 = t332 / 0.2e1;
t425 = -t334 / 0.2e1;
t423 = -t338 / 0.2e1;
t416 = mrSges(7,3) * t290;
t415 = Ifges(4,4) * t336;
t414 = Ifges(4,4) * t339;
t409 = Ifges(7,6) * t334;
t408 = t178 * mrSges(4,1);
t407 = t178 * mrSges(4,2);
t273 = -pkin(9) * t376 + t316;
t406 = t273 * mrSges(3,2);
t274 = t287 * qJD(2);
t405 = t274 * mrSges(3,1);
t404 = t197 * t242;
t397 = t332 * t339;
t209 = Ifges(6,5) * t282 - Ifges(6,6) * t283;
t251 = t331 * mrSges(5,2) * t374 + mrSges(5,1) * t363;
t387 = qJD(5) * t335;
t324 = -pkin(4) * t331 - pkin(3);
t217 = -pkin(5) * t343 - pkin(12) * t290 + t324;
t243 = t303 * t421 - t328 * t367;
t164 = t217 * t338 - t243 * t334;
t386 = qJD(6) * t164;
t165 = t217 * t334 + t243 * t338;
t385 = qJD(6) * t165;
t91 = -t127 * t334 + t227 * t338;
t32 = qJD(6) * t91 + t176 * t334 + t338 * t63;
t92 = t127 * t338 + t227 * t334;
t33 = -qJD(6) * t92 + t176 * t338 - t334 * t63;
t7 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t64;
t23 = Ifges(6,4) * t63 - Ifges(6,2) * t64 + Ifges(6,6) * t176;
t381 = t7 / 0.2e1 - t23 / 0.2e1;
t22 = Ifges(6,5) * t63 - Ifges(6,6) * t64 + Ifges(6,3) * t176;
t36 = Ifges(7,5) * t92 + Ifges(7,6) * t91 - Ifges(7,3) * t345;
t68 = Ifges(6,4) * t127 + Ifges(6,2) * t345 + Ifges(6,6) * t227;
t380 = t36 / 0.2e1 - t68 / 0.2e1;
t101 = Ifges(6,4) * t166 - Ifges(6,2) * t167 + Ifges(6,6) * t375;
t111 = qJD(6) * t185 + t166 * t338 + t334 * t375;
t112 = qJD(6) * t349 - t166 * t334 + t338 * t375;
t40 = Ifges(7,5) * t111 + Ifges(7,6) * t112 + Ifges(7,3) * t167;
t378 = t40 / 0.2e1 - t101 / 0.2e1;
t138 = Ifges(6,4) * t202 + Ifges(6,2) * t344 - Ifges(6,6) * t401;
t95 = -Ifges(7,5) * t349 + Ifges(7,6) * t185 - Ifges(7,3) * t344;
t377 = t95 / 0.2e1 - t138 / 0.2e1;
t100 = Ifges(6,5) * t166 - Ifges(6,6) * t167 + Ifges(6,3) * t375;
t115 = Ifges(4,5) * t177 - Ifges(4,6) * t176 + Ifges(4,3) * t364;
t134 = -Ifges(7,5) * t341 - Ifges(7,6) * t342 + Ifges(7,3) * t283;
t210 = Ifges(6,4) * t282 - Ifges(6,2) * t283;
t371 = t134 / 0.2e1 - t210 / 0.2e1;
t189 = -Ifges(7,3) * t343 + (Ifges(7,5) * t338 - t409) * t290;
t231 = Ifges(6,4) * t290 + Ifges(6,2) * t343;
t370 = t189 / 0.2e1 - t231 / 0.2e1;
t369 = Ifges(6,5) * t290 / 0.2e1 + Ifges(6,6) * t343 / 0.2e1 + Ifges(5,5) * t328 / 0.2e1 + Ifges(5,6) * t427;
t27 = t64 * mrSges(6,1) + t63 * mrSges(6,2);
t90 = -t153 * mrSges(5,1) + t154 * mrSges(5,2);
t108 = t167 * mrSges(6,1) + t166 * mrSges(6,2);
t208 = t283 * mrSges(6,1) + t282 * mrSges(6,2);
t365 = -qJD(3) * t214 - t224 * t388 - t336 * t239 - t256 * t375;
t19 = pkin(12) * t227 + t417;
t125 = -pkin(3) * t280 - t130;
t86 = -pkin(4) * t183 + t125;
t39 = -pkin(5) * t345 - pkin(12) * t127 + t86;
t10 = -t19 * t334 + t338 * t39;
t82 = -t240 * t397 + (-pkin(3) * t376 - t259 * t339) * t329 - t365;
t55 = -pkin(4) * t153 + t82;
t13 = pkin(5) * t64 - pkin(12) * t63 + t55;
t5 = t335 * t26 + t421 * t29 + t47 * t368 - t387 * t50;
t3 = pkin(12) * t176 + t5;
t1 = qJD(6) * t10 + t13 * t334 + t3 * t338;
t11 = t19 * t338 + t334 * t39;
t2 = -qJD(6) * t11 + t13 * t338 - t3 * t334;
t360 = t1 * t338 - t2 * t334;
t307 = -mrSges(7,1) * t338 + mrSges(7,2) * t334;
t359 = mrSges(7,1) * t334 + mrSges(7,2) * t338;
t53 = t335 * t152 + t156 * t368 + t421 * t169 - t170 * t387;
t51 = pkin(12) * t375 + t53;
t317 = pkin(10) * t402;
t266 = t317 + (-pkin(2) * t339 - pkin(3)) * t332;
t203 = -pkin(4) * t279 + t266;
t122 = -pkin(5) * t344 - pkin(12) * t202 + t203;
t94 = -pkin(12) * t401 + t465;
t57 = t122 * t338 - t334 * t94;
t272 = t286 * qJD(3);
t238 = pkin(4) * t363 + t272;
t87 = pkin(5) * t167 - pkin(12) * t166 + t238;
t16 = qJD(6) * t57 + t334 * t87 + t338 * t51;
t58 = t122 * t334 + t338 * t94;
t17 = -qJD(6) * t58 - t334 * t51 + t338 * t87;
t356 = t16 * t338 - t17 * t334;
t37 = Ifges(7,4) * t92 + Ifges(7,2) * t91 - Ifges(7,6) * t345;
t38 = Ifges(7,1) * t92 + Ifges(7,4) * t91 - Ifges(7,5) * t345;
t353 = t37 * t425 + t38 * t422;
t96 = -Ifges(7,4) * t349 + Ifges(7,2) * t185 - Ifges(7,6) * t344;
t97 = -Ifges(7,1) * t349 + Ifges(7,4) * t185 - Ifges(7,5) * t344;
t352 = t422 * t97 + t425 * t96;
t20 = -t335 * t50 + t421 * t47;
t106 = t156 * t421 - t335 * t170;
t314 = Ifges(3,5) * t340 * t391;
t313 = Ifges(4,5) * t374;
t306 = Ifges(5,1) * t328 + t412;
t305 = Ifges(5,2) * t331 + t413;
t302 = -mrSges(5,1) * t331 + mrSges(5,2) * t328;
t293 = t359 * qJD(6);
t292 = -mrSges(4,2) * t332 + mrSges(4,3) * t401;
t291 = mrSges(4,1) * t332 - mrSges(4,3) * t402;
t285 = -pkin(9) * t330 * t337 + t322;
t284 = pkin(2) * t397 - t317;
t270 = (Ifges(4,1) * t339 - t415) * t390;
t269 = (-Ifges(4,2) * t336 + t414) * t390;
t268 = -Ifges(4,6) * t375 + t313;
t267 = (mrSges(4,1) * t336 + mrSges(4,2) * t339) * t390;
t262 = Ifges(4,5) * t332 + (Ifges(4,1) * t336 + t414) * t329;
t261 = Ifges(4,6) * t332 + (Ifges(4,2) * t339 + t415) * t329;
t258 = (mrSges(5,1) * t336 - mrSges(5,3) * t399) * t390;
t257 = (-mrSges(5,3) * t328 * t339 - mrSges(5,2) * t336) * t390;
t247 = -mrSges(5,1) * t401 - mrSges(5,3) * t281;
t246 = mrSges(5,2) * t401 + mrSges(5,3) * t279;
t229 = -mrSges(6,1) * t343 + mrSges(6,2) * t290;
t221 = (Ifges(5,3) * t336 + (Ifges(5,5) * t331 - Ifges(5,6) * t328) * t339) * t390;
t220 = -mrSges(7,1) * t343 - t338 * t416;
t219 = mrSges(7,2) * t343 - t334 * t416;
t213 = pkin(5) * t283 - pkin(12) * t282;
t212 = t359 * t290;
t204 = -mrSges(5,1) * t279 + mrSges(5,2) * t281;
t200 = Ifges(5,1) * t281 + Ifges(5,4) * t279 - Ifges(5,5) * t401;
t199 = Ifges(5,4) * t281 + Ifges(5,2) * t279 - Ifges(5,6) * t401;
t198 = Ifges(5,5) * t281 + Ifges(5,6) * t279 - Ifges(5,3) * t401;
t196 = t331 * t366 - t303 * t387 + (-qJD(5) * t362 - t382) * t328;
t193 = -mrSges(6,1) * t401 - mrSges(6,3) * t202;
t192 = mrSges(6,2) * t401 + mrSges(6,3) * t344;
t188 = mrSges(4,1) * t280 - mrSges(4,3) * t228;
t187 = -mrSges(4,2) * t280 - mrSges(4,3) * t227;
t180 = -mrSges(7,2) * t283 - mrSges(7,3) * t342;
t179 = mrSges(7,1) * t283 + mrSges(7,3) * t341;
t158 = mrSges(4,1) * t364 - mrSges(4,3) * t177;
t157 = -mrSges(4,2) * t364 - mrSges(4,3) * t176;
t149 = mrSges(7,1) * t342 - mrSges(7,2) * t341;
t146 = -mrSges(6,2) * t375 - mrSges(6,3) * t167;
t145 = mrSges(6,1) * t375 - mrSges(6,3) * t166;
t144 = Ifges(4,1) * t228 - Ifges(4,4) * t227 + Ifges(4,5) * t280;
t143 = Ifges(4,4) * t228 - Ifges(4,2) * t227 + Ifges(4,6) * t280;
t142 = -mrSges(6,1) * t344 + mrSges(6,2) * t202;
t141 = mrSges(5,1) * t227 - mrSges(5,3) * t184;
t140 = -mrSges(5,2) * t227 + mrSges(5,3) * t183;
t137 = Ifges(6,5) * t202 + Ifges(6,6) * t344 - Ifges(6,3) * t401;
t133 = -mrSges(7,1) * t344 + mrSges(7,3) * t349;
t132 = mrSges(7,2) * t344 + mrSges(7,3) * t185;
t129 = -mrSges(7,1) * t185 - mrSges(7,2) * t349;
t128 = -mrSges(5,1) * t183 + mrSges(5,2) * t184;
t123 = mrSges(4,1) * t176 + mrSges(4,2) * t177;
t117 = Ifges(4,1) * t177 - Ifges(4,4) * t176 + Ifges(4,5) * t364;
t116 = Ifges(4,4) * t177 - Ifges(4,2) * t176 + Ifges(4,6) * t364;
t114 = mrSges(5,1) * t176 - mrSges(5,3) * t154;
t113 = -mrSges(5,2) * t176 + mrSges(5,3) * t153;
t105 = Ifges(5,1) * t184 + Ifges(5,4) * t183 + Ifges(5,5) * t227;
t104 = Ifges(5,4) * t184 + Ifges(5,2) * t183 + Ifges(5,6) * t227;
t103 = Ifges(5,5) * t184 + Ifges(5,6) * t183 + Ifges(5,3) * t227;
t99 = mrSges(6,1) * t227 - mrSges(6,3) * t127;
t98 = -mrSges(6,2) * t227 + mrSges(6,3) * t345;
t93 = pkin(5) * t401 - t106;
t89 = -t196 * t334 + t213 * t338 - t385;
t88 = t196 * t338 + t213 * t334 + t386;
t85 = (t240 * t332 + t259 * t329) * t339 + t365;
t81 = -mrSges(7,2) * t167 + mrSges(7,3) * t112;
t80 = mrSges(7,1) * t167 - mrSges(7,3) * t111;
t75 = Ifges(5,5) * t154 + Ifges(5,6) * t153 + Ifges(5,3) * t176;
t72 = -mrSges(6,1) * t345 + mrSges(6,2) * t127;
t67 = Ifges(6,5) * t127 + Ifges(6,6) * t345 + Ifges(6,3) * t227;
t66 = -mrSges(7,1) * t345 - mrSges(7,3) * t92;
t65 = mrSges(7,2) * t345 + mrSges(7,3) * t91;
t56 = -mrSges(7,1) * t112 + mrSges(7,2) * t111;
t52 = -pkin(5) * t375 - t54;
t48 = -mrSges(7,1) * t91 + mrSges(7,2) * t92;
t46 = -mrSges(6,2) * t176 - mrSges(6,3) * t64;
t45 = mrSges(6,1) * t176 - mrSges(6,3) * t63;
t42 = Ifges(7,1) * t111 + Ifges(7,4) * t112 + Ifges(7,5) * t167;
t41 = Ifges(7,4) * t111 + Ifges(7,2) * t112 + Ifges(7,6) * t167;
t18 = -t227 * pkin(5) - t20;
t15 = -mrSges(7,2) * t64 + mrSges(7,3) * t33;
t14 = mrSges(7,1) * t64 - mrSges(7,3) * t32;
t12 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t9 = Ifges(7,1) * t32 + Ifges(7,4) * t33 + Ifges(7,5) * t64;
t8 = Ifges(7,4) * t32 + Ifges(7,2) * t33 + Ifges(7,6) * t64;
t4 = -t176 * pkin(5) - t6;
t21 = [-(t7 - t23) * t345 + (t314 - 0.2e1 * t405 - 0.2e1 * t406) * t333 + 0.2e1 * m(3) * (t273 * t287 - t274 * t285) + (t20 * t6 + t417 * t5 + t55 * t86) * t461 + 0.2e1 * t417 * t46 + (t117 + 0.2e1 * t407) * t228 + (-t116 + t22 + t75 + 0.2e1 * t408) * t227 + (0.2e1 * (t273 * t340 + t274 * t337) * mrSges(3,3) + ((t285 * t459 + Ifges(3,5) * t333 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t340) * t466) * t340 + (t329 * (Ifges(4,5) * t228 - Ifges(4,6) * t227 + Ifges(4,3) * t280) + t287 * t459 - 0.2e1 * Ifges(3,6) * t333 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t337 + (Ifges(3,1) - Ifges(3,2)) * t340) * t466) * t337) * qJD(2)) * t330 + (t36 - t68) * t64 + (t67 + t103 - t143) * t176 + (t1 * t11 + t10 * t2 + t18 * t4) * t460 + (t125 * t82 + t34 * t70 + t35 * t71) * t462 + (t130 * t85 + t131 * t84 + t171 * t178) * t463 + 0.2e1 * t10 * t14 + 0.2e1 * t11 * t15 + 0.2e1 * t18 * t12 + t33 * t37 + t32 * t38 + 0.2e1 * t20 * t45 + 0.2e1 * t4 * t48 + 0.2e1 * t1 * t65 + 0.2e1 * t2 * t66 + t63 * t69 + 0.2e1 * t55 * t72 + 0.2e1 * t86 * t27 + t91 * t8 + t92 * t9 + 0.2e1 * t5 * t98 + 0.2e1 * t6 * t99 + 0.2e1 * t71 * t113 + 0.2e1 * t70 * t114 + 0.2e1 * t125 * t90 + t127 * t24 + 0.2e1 * t82 * t128 + 0.2e1 * t35 * t140 + 0.2e1 * t34 * t141 + t153 * t104 + t154 * t105 + 0.2e1 * t131 * t157 + 0.2e1 * t130 * t158 + 0.2e1 * t171 * t123 + t177 * t144 + t183 * t76 + t184 * t77 + 0.2e1 * t84 * t187 + 0.2e1 * t85 * t188 + t280 * t115; m(6) * (t106 * t6 + t20 * t54 + t203 * t55 + t238 * t86 + t417 * t53 + t465 * t5) + t465 * t46 - t378 * t345 - t381 * t344 + t417 * t146 + (t128 - t188) * t272 + m(4) * (-t130 * t272 + t131 * t271 + t284 * t85 + t286 * t84) + ((t117 / 0.2e1 + t407) * t336 + (t116 / 0.2e1 - t75 / 0.2e1 - t22 / 0.2e1 - t408) * t339 + (Ifges(4,3) * t426 + (Ifges(4,5) * t336 + Ifges(4,6) * t339) * t329 / 0.2e1) * t376 + (-m(4) * t178 - t123) * pkin(2) + ((t105 * t427 + t104 * t428 - t130 * mrSges(4,3) + t144 / 0.2e1) * t339 + (-t131 * mrSges(4,3) + t67 / 0.2e1 + t103 / 0.2e1 - t143 / 0.2e1) * t336) * qJD(3)) * t329 + t380 * t167 + t377 * t64 - Ifges(3,6) * t376 + t314 - t405 - t406 + m(7) * (t1 * t58 + t10 * t17 + t11 * t16 + t18 * t52 + t2 * t57 + t4 * t93) + m(5) * (t125 * t272 + t181 * t70 + t182 * t71 + t194 * t34 + t195 * t35 + t266 * t82) + (t100 / 0.2e1 + t221 / 0.2e1 - t269 / 0.2e1) * t227 + t115 * t426 + t184 * t438 + t183 * t439 + t9 * t443 + t8 * t444 + t200 * t445 + t199 * t446 + t63 * t447 + t127 * t450 + t281 * t451 + t279 * t452 + t166 * t453 + t202 * t454 + t52 * t48 + t18 * t56 + t57 * t14 + t58 * t15 + t16 * t65 + t17 * t66 + (t137 / 0.2e1 + t198 / 0.2e1 - t261 / 0.2e1) * t176 + t10 * t80 + t11 * t81 + t91 * t41 / 0.2e1 + t92 * t42 / 0.2e1 + t93 * t12 + t33 * t96 / 0.2e1 + t32 * t97 / 0.2e1 + t53 * t98 + t54 * t99 + t106 * t45 + t86 * t108 + t111 * t38 / 0.2e1 + t112 * t37 / 0.2e1 + t4 * t129 + t1 * t132 + t2 * t133 + t55 * t142 + t20 * t145 + t181 * t141 + t182 * t140 + t5 * t192 + t6 * t193 + t194 * t114 + t195 * t113 + t203 * t27 + t82 * t204 + t238 * t72 + t35 * t246 + t34 * t247 + t125 * t251 + t71 * t257 + t70 * t258 + t177 * t262 / 0.2e1 + t266 * t90 + t171 * t267 + t228 * t270 / 0.2e1 + t271 * t187 + t280 * t268 / 0.2e1 + t284 * t158 + t286 * t157 + t85 * t291 + t84 * t292; (t106 * t54 + t203 * t238 + t465 * t53) * t461 + 0.2e1 * t465 * t146 - t349 * t42 - (t40 - t101) * t344 + 0.2e1 * (t204 - t291) * t272 + (-0.2e1 * pkin(2) * t267 + t336 * t270 + (-t100 - t221 + t269) * t339 + ((-t199 * t328 + t200 * t331 + t284 * t458 + t262) * t339 + (t286 * t458 + t137 + t198 - t261) * t336) * qJD(3)) * t329 + (t95 - t138) * t167 + (t16 * t58 + t17 * t57 + t52 * t93) * t460 + (t181 * t194 + t182 * t195 + t266 * t272) * t462 + (t271 * t286 - t272 * t284) * t463 + 0.2e1 * t57 * t80 + 0.2e1 * t58 * t81 + 0.2e1 * t93 * t56 + t111 * t97 + t112 * t96 + 0.2e1 * t52 * t129 + 0.2e1 * t16 * t132 + 0.2e1 * t17 * t133 + 0.2e1 * t106 * t145 + t166 * t139 + t185 * t41 + 0.2e1 * t53 * t192 + 0.2e1 * t54 * t193 + t202 * t102 + 0.2e1 * t203 * t108 + 0.2e1 * t238 * t142 + 0.2e1 * t182 * t246 + 0.2e1 * t181 * t247 + 0.2e1 * t195 * t257 + 0.2e1 * t194 * t258 + 0.2e1 * t266 * t251 + t279 * t222 + t281 * t223 + 0.2e1 * t271 * t292 + t332 * t268; -t371 * t345 - (-t5 * mrSges(6,3) + t381) * t343 + (-mrSges(6,3) * t417 + t380) * t283 + m(6) * (t196 * t417 - t197 * t20 - t242 * t6 + t243 * t5 + t324 * t55) + (t12 - t45) * t242 + (t48 - t99) * t197 + (-t34 * mrSges(5,3) - qJ(4) * t114 - qJD(4) * t141 + t451) * t328 + (t35 * mrSges(5,3) + qJ(4) * t113 + qJD(4) * t140 + t452) * t331 + (-t20 * mrSges(6,3) + t353 + t453) * t282 + (-t6 * mrSges(6,3) + t8 * t425 + t9 * t422 + t454 + (t37 * t423 + t38 * t425) * qJD(6)) * t290 + t369 * t176 + t370 * t64 + t115 + m(7) * (t1 * t165 + t10 * t89 + t11 * t88 + t164 * t2 + t18 * t197 + t242 * t4) + m(5) * (-pkin(3) * t82 + (-t328 * t70 + t331 * t71) * qJD(4) + (-t328 * t34 + t331 * t35) * qJ(4)) + t63 * t437 + t127 * t440 + t32 * t441 + t33 * t442 + t306 * t445 + t305 * t446 + t92 * t448 + t91 * t449 - t84 * mrSges(4,2) + t85 * mrSges(4,1) + t88 * t65 + t89 * t66 - pkin(3) * t90 + t18 * t149 + t164 * t14 + t165 * t15 + t10 * t179 + t11 * t180 + t196 * t98 + t86 * t208 + t4 * t212 + t1 * t219 + t2 * t220 + t227 * t209 / 0.2e1 + t55 * t229 + t243 * t46 + t82 * t302 + t324 * t27; (-mrSges(6,3) * t465 + t377) * t283 + m(6) * (-t106 * t197 + t196 * t465 + t238 * t324 - t242 * t54 + t243 * t53) - t371 * t344 - (-t53 * mrSges(6,3) + t378) * t343 + (t129 - t193) * t197 + (-mrSges(4,1) + t302) * t272 + (-t106 * mrSges(6,3) + t352 + t447) * t282 + (-t54 * mrSges(6,3) + t41 * t425 + t42 * t422 + t450 + (t423 * t96 + t425 * t97) * qJD(6)) * t290 + (-t181 * mrSges(5,3) - qJ(4) * t258 - qJD(4) * t247 + t438) * t328 + (t182 * mrSges(5,3) + qJ(4) * t257 + qJD(4) * t246 + t439) * t331 + (-t339 * t209 / 0.2e1 + ((t305 * t428 + t306 * t427) * t339 + (-Ifges(4,6) + t369) * t336) * qJD(3)) * t329 + (t56 - t145) * t242 + t370 * t167 + t313 + m(7) * (t16 * t165 + t164 * t17 + t197 * t93 + t242 * t52 + t57 * t89 + t58 * t88) + m(5) * (-pkin(3) * t272 + (-t194 * t328 + t195 * t331) * qJD(4) + (-t181 * t328 + t182 * t331) * qJ(4)) + t166 * t437 + t202 * t440 + t111 * t441 + t112 * t442 + t136 * t443 + t135 * t444 + t88 * t132 + t89 * t133 + t93 * t149 + t164 * t80 + t165 * t81 + t57 * t179 + t58 * t180 + t196 * t192 + t203 * t208 + t52 * t212 + t16 * t219 + t17 * t220 + t238 * t229 + t243 * t146 - pkin(3) * t251 - t271 * mrSges(4,2) + t324 * t108; t149 * t455 + 0.2e1 * t164 * t179 + 0.2e1 * t165 * t180 + t212 * t456 + 0.2e1 * t324 * t208 + 0.2e1 * t88 * t219 + 0.2e1 * t89 * t220 + (t164 * t89 + t165 * t88 + t404) * t460 + (t196 * t243 + t404) * t461 - (t196 * t457 + t134 - t210) * t343 + (t243 * t457 + t189 - t231) * t283 + (mrSges(6,3) * t455 - t190 * t334 + t191 * t338 + t232) * t282 + (mrSges(6,3) * t456 - t334 * t135 + t338 * t136 + t211 + (-t190 * t338 - t191 * t334) * qJD(6)) * t290 + (qJ(4) * t462 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t328 ^ 2 + t331 ^ 2); t338 * t14 + t334 * t15 + (-t334 * t66 + t338 * t65) * qJD(6) + m(7) * (t1 * t334 + t2 * t338 + (-t10 * t334 + t11 * t338) * qJD(6)) + m(6) * t55 + m(5) * t82 + t90 + t27; t334 * t81 + t338 * t80 + (t132 * t338 - t133 * t334) * qJD(6) + m(7) * (t16 * t334 + t17 * t338 + (-t334 * t57 + t338 * t58) * qJD(6)) + m(6) * t238 + m(5) * t272 + t108 + t251; m(7) * (t334 * t88 + t338 * t89 + (-t164 * t334 + t165 * t338) * qJD(6)) + t219 * t383 + t334 * t180 - t220 * t384 + t338 * t179 + t208; 0; -t5 * mrSges(6,2) + t6 * mrSges(6,1) + t18 * t293 - t345 * t435 + t91 * t434 + t92 * t433 + t4 * t307 + t64 * t432 + t33 * t430 + t32 * t429 + t9 * t424 + t8 * t422 + t353 * qJD(6) + (-m(7) * t4 - t12) * pkin(5) + ((-t10 * t338 - t11 * t334) * qJD(6) + t360) * mrSges(7,3) + (-t66 * t383 - t65 * t384 - t334 * t14 + t338 * t15 + m(7) * (-t10 * t383 - t11 * t384 + t360)) * pkin(12) + t22; -t53 * mrSges(6,2) + t54 * mrSges(6,1) + t93 * t293 - t344 * t435 + t185 * t434 - t349 * t433 + t52 * t307 + t167 * t432 + t112 * t430 + t111 * t429 + t42 * t424 + t41 * t422 + t352 * qJD(6) + (-m(7) * t52 - t56) * pkin(5) + ((-t334 * t58 - t338 * t57) * qJD(6) + t356) * mrSges(7,3) + (-t133 * t383 - t132 * t384 - t334 * t80 + t338 * t81 + m(7) * (-t383 * t57 - t384 * t58 + t356)) * pkin(12) + t100; -pkin(5) * t149 - t196 * mrSges(6,2) + t242 * t293 - t343 * t435 + t283 * t432 + (-m(7) * pkin(5) - mrSges(6,1) + t307) * t197 + (t88 * mrSges(7,3) + t290 * t433 + t282 * t429 + t449 + (-t164 * mrSges(7,3) + t290 * t431 + t441) * qJD(6) + (m(7) * (t88 - t386) + t180 - qJD(6) * t220) * pkin(12)) * t338 + (-t89 * mrSges(7,3) + t295 * t436 + t282 * t431 + t448 + (t310 * t436 - t165 * mrSges(7,3) - t190 / 0.2e1) * qJD(6) + (-qJD(6) * t219 - t179 + m(7) * (-t89 - t385)) * pkin(12)) * t334 + t209; 0; -0.2e1 * pkin(5) * t293 + t295 * t338 + t296 * t334 + (-t309 * t334 + t310 * t338) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t17 - mrSges(7,2) * t16 + t40; mrSges(7,1) * t89 - mrSges(7,2) * t88 + t134; -t293; t325 + (pkin(12) * t307 - t409) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
