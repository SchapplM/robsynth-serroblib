% Calculate time derivative of joint inertia matrix for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:29:06
% EndTime: 2019-03-10 04:29:29
% DurationCPUTime: 10.34s
% Computational Cost: add. (25078->881), mult. (64537->1283), div. (0->0), fcn. (64049->12), ass. (0->350)
t328 = sin(qJ(4));
t333 = cos(qJ(4));
t329 = sin(qJ(3));
t373 = qJD(4) * t329;
t334 = cos(qJ(3));
t375 = qJD(3) * t334;
t340 = -t328 * t373 + t333 * t375;
t372 = qJD(4) * t333;
t339 = t328 * t375 + t329 * t372;
t293 = -pkin(3) * t334 - pkin(10) * t329 - pkin(2);
t380 = t333 * t334;
t313 = pkin(9) * t380;
t245 = t328 * t293 + t313;
t441 = qJD(4) * t245;
t327 = sin(qJ(5));
t332 = cos(qJ(5));
t276 = t327 * t333 + t328 * t332;
t326 = sin(qJ(6));
t331 = cos(qJ(6));
t346 = t327 * t328 - t332 * t333;
t203 = -t276 * t326 - t331 * t346;
t417 = t203 / 0.2e1;
t204 = t276 * t331 - t326 * t346;
t416 = t204 / 0.2e1;
t406 = -t346 / 0.2e1;
t405 = t276 / 0.2e1;
t440 = t328 / 0.2e1;
t401 = t333 / 0.2e1;
t257 = t346 * t329;
t324 = sin(pkin(6));
t330 = sin(qJ(2));
t386 = t324 * t330;
t308 = pkin(8) * t386;
t325 = cos(pkin(6));
t335 = cos(qJ(2));
t399 = pkin(1) * t335;
t247 = t308 + (-pkin(2) - t399) * t325;
t260 = -t325 * t334 + t329 * t386;
t261 = t325 * t329 + t334 * t386;
t169 = t260 * pkin(3) - t261 * pkin(10) + t247;
t385 = t324 * t335;
t268 = t325 * t330 * pkin(1) + pkin(8) * t385;
t248 = pkin(9) * t325 + t268;
t249 = (-pkin(2) * t335 - pkin(9) * t330 - pkin(1)) * t324;
t179 = t334 * t248 + t329 * t249;
t171 = -pkin(10) * t385 + t179;
t99 = t328 * t169 + t333 * t171;
t274 = t333 * t293;
t381 = t329 * t333;
t398 = pkin(9) * t328;
t205 = -pkin(11) * t381 + t274 + (-pkin(4) - t398) * t334;
t382 = t328 * t329;
t225 = -pkin(11) * t382 + t245;
t147 = t327 * t205 + t332 * t225;
t427 = -pkin(11) - pkin(10);
t300 = t427 * t328;
t301 = t427 * t333;
t229 = t327 * t300 - t332 * t301;
t294 = -mrSges(5,1) * t333 + mrSges(5,2) * t328;
t439 = -m(5) * pkin(3) + t294;
t267 = t325 * t399 - t308;
t438 = qJD(4) + qJD(5);
t437 = 0.2e1 * m(5);
t436 = 2 * m(6);
t435 = 2 * m(7);
t434 = 0.2e1 * pkin(9);
t433 = -2 * mrSges(3,3);
t223 = -t261 * t328 - t333 * t385;
t342 = -t261 * t333 + t328 * t385;
t149 = t223 * t332 + t327 * t342;
t150 = t223 * t327 - t332 * t342;
t91 = t149 * t331 - t150 * t326;
t92 = t149 * t326 + t150 * t331;
t45 = Ifges(7,4) * t92 + Ifges(7,2) * t91 + Ifges(7,6) * t260;
t431 = t45 / 0.2e1;
t46 = Ifges(7,1) * t92 + Ifges(7,4) * t91 + Ifges(7,5) * t260;
t430 = t46 / 0.2e1;
t429 = t91 / 0.2e1;
t428 = t92 / 0.2e1;
t256 = t276 * t329;
t185 = -t256 * t331 + t257 * t326;
t186 = -t256 * t326 - t257 * t331;
t110 = Ifges(7,4) * t186 + Ifges(7,2) * t185 - Ifges(7,6) * t334;
t426 = t110 / 0.2e1;
t111 = Ifges(7,1) * t186 + Ifges(7,4) * t185 - Ifges(7,5) * t334;
t425 = t111 / 0.2e1;
t378 = qJD(2) * t324;
t361 = t335 * t378;
t222 = -qJD(3) * t260 + t334 * t361;
t377 = qJD(2) * t330;
t362 = t324 * t377;
t133 = qJD(4) * t342 - t222 * t328 + t333 * t362;
t424 = t133 / 0.2e1;
t141 = Ifges(7,4) * t204 + Ifges(7,2) * t203;
t423 = t141 / 0.2e1;
t142 = Ifges(7,1) * t204 + Ifges(7,4) * t203;
t422 = t142 / 0.2e1;
t421 = t149 / 0.2e1;
t420 = t150 / 0.2e1;
t419 = t185 / 0.2e1;
t418 = t186 / 0.2e1;
t215 = t438 * t346;
t415 = -t215 / 0.2e1;
t216 = t438 * t276;
t414 = -t216 / 0.2e1;
t219 = Ifges(6,4) * t276 - Ifges(6,2) * t346;
t413 = t219 / 0.2e1;
t220 = Ifges(6,1) * t276 - Ifges(6,4) * t346;
t412 = t220 / 0.2e1;
t411 = t223 / 0.2e1;
t410 = -t342 / 0.2e1;
t394 = Ifges(5,4) * t328;
t349 = Ifges(5,1) * t333 - t394;
t252 = -Ifges(5,5) * t334 + t329 * t349;
t409 = t252 / 0.2e1;
t408 = -t256 / 0.2e1;
t407 = -t257 / 0.2e1;
t393 = Ifges(5,4) * t333;
t298 = Ifges(5,1) * t328 + t393;
t404 = t298 / 0.2e1;
t403 = -t328 / 0.2e1;
t402 = -t333 / 0.2e1;
t400 = m(5) * t334;
t397 = pkin(9) * t334;
t323 = t329 * pkin(9);
t98 = t333 * t169 - t171 * t328;
t81 = pkin(4) * t260 + pkin(11) * t342 + t98;
t90 = pkin(11) * t223 + t99;
t41 = t327 * t81 + t332 * t90;
t396 = Ifges(4,4) * t329;
t395 = Ifges(4,4) * t334;
t392 = Ifges(5,6) * t328;
t315 = pkin(4) * t332 + pkin(5);
t368 = qJD(6) * t331;
t369 = qJD(6) * t326;
t384 = t326 * t327;
t198 = t315 * t368 + (-t327 * t369 + (t331 * t332 - t384) * qJD(5)) * pkin(4);
t391 = t198 * mrSges(7,2);
t254 = t267 * qJD(2);
t390 = t254 * mrSges(3,2);
t255 = t268 * qJD(2);
t389 = t255 * mrSges(3,1);
t388 = t255 * mrSges(4,1);
t387 = t255 * mrSges(4,2);
t383 = t327 * t331;
t105 = qJD(6) * t203 - t215 * t331 - t216 * t326;
t106 = -qJD(6) * t204 + t215 * t326 - t216 * t331;
t63 = Ifges(7,5) * t105 + Ifges(7,6) * t106;
t152 = -Ifges(6,5) * t215 - Ifges(6,6) * t216;
t291 = (pkin(3) * t329 - pkin(10) * t334) * qJD(3);
t376 = qJD(3) * t329;
t379 = t333 * t291 + t376 * t398;
t292 = pkin(4) * t382 + t323;
t374 = qJD(4) * t328;
t371 = qJD(5) * t327;
t370 = qJD(5) * t332;
t134 = qJD(4) * t223 + t222 * t333 + t328 * t362;
t57 = qJD(5) * t149 + t133 * t327 + t134 * t332;
t58 = -qJD(5) * t150 + t133 * t332 - t134 * t327;
t22 = qJD(6) * t91 + t326 * t58 + t331 * t57;
t221 = qJD(3) * t261 + t329 * t361;
t23 = -qJD(6) * t92 - t326 * t57 + t331 * t58;
t6 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t221;
t24 = Ifges(6,5) * t57 + Ifges(6,6) * t58 + Ifges(6,3) * t221;
t160 = -t216 * t329 - t346 * t375;
t161 = t257 * t438 - t276 * t375;
t76 = qJD(6) * t185 + t160 * t331 + t161 * t326;
t77 = -qJD(6) * t186 - t160 * t326 + t161 * t331;
t36 = Ifges(7,5) * t76 + Ifges(7,6) * t77 + Ifges(7,3) * t376;
t366 = pkin(4) * t374;
t365 = Ifges(4,6) * t385;
t68 = Ifges(5,5) * t134 + Ifges(5,6) * t133 + Ifges(5,3) * t221;
t94 = Ifges(6,5) * t160 + Ifges(6,6) * t161 + Ifges(6,3) * t376;
t364 = Ifges(4,5) * t222 - Ifges(4,6) * t221 + Ifges(4,3) * t362;
t240 = pkin(4) * t339 + pkin(9) * t375;
t316 = -pkin(4) * t333 - pkin(3);
t363 = qJD(4) * t427;
t129 = -Ifges(5,1) * t342 + Ifges(5,4) * t223 + Ifges(5,5) * t260;
t356 = t129 * t401;
t40 = -t327 * t90 + t332 * t81;
t199 = -t315 * t369 + (-t327 * t368 + (-t326 * t332 - t383) * qJD(5)) * pkin(4);
t196 = t199 * mrSges(7,1);
t355 = t196 - t391;
t146 = t332 * t205 - t225 * t327;
t178 = -t329 * t248 + t249 * t334;
t228 = t332 * t300 + t301 * t327;
t176 = t328 * t291 + t293 * t372 + (-t333 * t376 - t334 * t374) * pkin(9);
t244 = -t328 * t397 + t274;
t354 = -qJD(4) * t244 + t176;
t353 = t362 / 0.2e1;
t320 = Ifges(5,5) * t372;
t352 = -Ifges(5,6) * t374 / 0.2e1 + t320 / 0.2e1 + t152 / 0.2e1 + t63 / 0.2e1;
t170 = pkin(3) * t385 - t178;
t351 = Ifges(5,5) * t440 + Ifges(6,5) * t405 + Ifges(7,5) * t416 + Ifges(5,6) * t401 + Ifges(6,6) * t406 + Ifges(7,6) * t417;
t350 = mrSges(5,1) * t328 + mrSges(5,2) * t333;
t348 = -Ifges(5,2) * t328 + t393;
t296 = Ifges(5,2) * t333 + t394;
t31 = pkin(5) * t260 - pkin(12) * t150 + t40;
t33 = pkin(12) * t149 + t41;
t13 = t31 * t331 - t326 * t33;
t14 = t31 * t326 + t33 * t331;
t253 = (pkin(2) * t330 - pkin(9) * t335) * t378;
t112 = -t248 * t376 + t249 * t375 + t329 * t253 + t334 * t254;
t107 = pkin(10) * t362 + t112;
t130 = t221 * pkin(3) - t222 * pkin(10) + t255;
t42 = t333 * t107 + t328 * t130 + t169 * t372 - t171 * t374;
t43 = -qJD(4) * t99 - t107 * t328 + t333 * t130;
t347 = -t43 * t328 + t42 * t333;
t116 = -pkin(5) * t334 + pkin(12) * t257 + t146;
t120 = -pkin(12) * t256 + t147;
t66 = t116 * t331 - t120 * t326;
t67 = t116 * t326 + t120 * t331;
t193 = -pkin(12) * t276 + t228;
t194 = -pkin(12) * t346 + t229;
t118 = t193 * t331 - t194 * t326;
t119 = t193 * t326 + t194 * t331;
t289 = t328 * t363;
t290 = t333 * t363;
t164 = t332 * t289 + t327 * t290 + t300 * t370 + t301 * t371;
t114 = -pkin(12) * t216 + t164;
t165 = -qJD(5) * t229 - t289 * t327 + t332 * t290;
t115 = pkin(12) * t215 + t165;
t49 = qJD(6) * t118 + t114 * t331 + t115 * t326;
t50 = -qJD(6) * t119 - t114 * t326 + t115 * t331;
t345 = t50 * mrSges(7,1) - t49 * mrSges(7,2) + t63;
t113 = -t248 * t375 - t249 * t376 + t253 * t334 - t329 * t254;
t30 = pkin(4) * t221 - pkin(11) * t134 + t43;
t35 = pkin(11) * t133 + t42;
t12 = -qJD(5) * t41 + t332 * t30 - t327 * t35;
t4 = pkin(5) * t221 - pkin(12) * t57 + t12;
t11 = t327 * t30 + t332 * t35 + t81 * t370 - t371 * t90;
t5 = pkin(12) * t58 + t11;
t2 = qJD(6) * t13 + t326 * t4 + t331 * t5;
t3 = -qJD(6) * t14 - t326 * t5 + t331 * t4;
t344 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t6;
t143 = (pkin(4) * t329 - pkin(11) * t380) * qJD(3) + (-t313 + (pkin(11) * t329 - t293) * t328) * qJD(4) + t379;
t155 = -pkin(11) * t339 + t176;
t61 = -qJD(5) * t147 + t332 * t143 - t155 * t327;
t51 = pkin(5) * t376 - pkin(12) * t160 + t61;
t60 = t327 * t143 + t332 * t155 + t205 * t370 - t225 * t371;
t52 = pkin(12) * t161 + t60;
t16 = qJD(6) * t66 + t326 * t51 + t331 * t52;
t17 = -qJD(6) * t67 - t326 * t52 + t331 * t51;
t343 = t17 * mrSges(7,1) - t16 * mrSges(7,2) + t36;
t124 = -pkin(4) * t223 + t170;
t341 = (-mrSges(6,1) * t327 - mrSges(6,2) * t332) * qJD(5) * pkin(4);
t108 = -pkin(3) * t362 - t113;
t338 = t165 * mrSges(6,1) - t164 * mrSges(6,2) + t152 + t345;
t337 = t12 * mrSges(6,1) - t11 * mrSges(6,2) + t24 + t344;
t336 = t61 * mrSges(6,1) - t60 * mrSges(6,2) + t343 + t94;
t73 = -pkin(4) * t133 + t108;
t189 = Ifges(5,5) * t340 - Ifges(5,6) * t339 + Ifges(5,3) * t376;
t321 = Ifges(4,5) * t375;
t303 = Ifges(3,5) * t361;
t299 = Ifges(4,1) * t329 + t395;
t297 = Ifges(4,2) * t334 + t396;
t288 = -mrSges(5,1) * t334 - mrSges(5,3) * t381;
t287 = mrSges(5,2) * t334 - mrSges(5,3) * t382;
t286 = (Ifges(4,1) * t334 - t396) * qJD(3);
t285 = t349 * qJD(4);
t284 = (-Ifges(4,2) * t329 + t395) * qJD(3);
t283 = t348 * qJD(4);
t281 = (mrSges(4,1) * t329 + mrSges(4,2) * t334) * qJD(3);
t280 = t350 * qJD(4);
t271 = (-mrSges(7,1) * t326 - mrSges(7,2) * t331) * qJD(6) * pkin(5);
t269 = t350 * t329;
t259 = pkin(4) * t383 + t315 * t326;
t258 = -pkin(4) * t384 + t315 * t331;
t251 = -Ifges(5,6) * t334 + t329 * t348;
t250 = -Ifges(5,3) * t334 + (Ifges(5,5) * t333 - t392) * t329;
t246 = pkin(5) * t346 + t316;
t236 = -mrSges(5,2) * t376 - mrSges(5,3) * t339;
t235 = mrSges(5,1) * t376 - mrSges(5,3) * t340;
t231 = -mrSges(6,1) * t334 + mrSges(6,3) * t257;
t230 = mrSges(6,2) * t334 - mrSges(6,3) * t256;
t227 = -mrSges(4,1) * t385 - t261 * mrSges(4,3);
t226 = mrSges(4,2) * t385 - t260 * mrSges(4,3);
t217 = mrSges(6,1) * t346 + mrSges(6,2) * t276;
t208 = pkin(5) * t256 + t292;
t202 = mrSges(5,1) * t339 + mrSges(5,2) * t340;
t195 = pkin(5) * t216 + t366;
t192 = mrSges(6,1) * t256 - mrSges(6,2) * t257;
t191 = -t298 * t373 + (Ifges(5,5) * t329 + t334 * t349) * qJD(3);
t190 = -t296 * t373 + (Ifges(5,6) * t329 + t334 * t348) * qJD(3);
t188 = mrSges(4,1) * t362 - mrSges(4,3) * t222;
t187 = -mrSges(4,2) * t362 - mrSges(4,3) * t221;
t184 = Ifges(4,1) * t261 - Ifges(4,4) * t260 - Ifges(4,5) * t385;
t183 = Ifges(4,4) * t261 - Ifges(4,2) * t260 - t365;
t182 = -Ifges(6,1) * t257 - Ifges(6,4) * t256 - Ifges(6,5) * t334;
t181 = -Ifges(6,4) * t257 - Ifges(6,2) * t256 - Ifges(6,6) * t334;
t180 = -Ifges(6,5) * t257 - Ifges(6,6) * t256 - Ifges(6,3) * t334;
t177 = t379 - t441;
t175 = mrSges(5,1) * t260 + mrSges(5,3) * t342;
t174 = -mrSges(5,2) * t260 + mrSges(5,3) * t223;
t173 = -mrSges(7,1) * t334 - mrSges(7,3) * t186;
t172 = mrSges(7,2) * t334 + mrSges(7,3) * t185;
t159 = -mrSges(5,1) * t223 - mrSges(5,2) * t342;
t156 = mrSges(4,1) * t221 + mrSges(4,2) * t222;
t154 = -Ifges(6,1) * t215 - Ifges(6,4) * t216;
t153 = -Ifges(6,4) * t215 - Ifges(6,2) * t216;
t151 = mrSges(6,1) * t216 - mrSges(6,2) * t215;
t145 = -mrSges(6,2) * t376 + mrSges(6,3) * t161;
t144 = mrSges(6,1) * t376 - mrSges(6,3) * t160;
t139 = -mrSges(7,1) * t203 + mrSges(7,2) * t204;
t138 = Ifges(4,1) * t222 - Ifges(4,4) * t221 + Ifges(4,5) * t362;
t137 = Ifges(4,4) * t222 - Ifges(4,2) * t221 + Ifges(4,6) * t362;
t128 = -Ifges(5,4) * t342 + Ifges(5,2) * t223 + Ifges(5,6) * t260;
t127 = -Ifges(5,5) * t342 + Ifges(5,6) * t223 + Ifges(5,3) * t260;
t123 = -pkin(5) * t161 + t240;
t122 = mrSges(6,1) * t260 - mrSges(6,3) * t150;
t121 = -mrSges(6,2) * t260 + mrSges(6,3) * t149;
t117 = -mrSges(7,1) * t185 + mrSges(7,2) * t186;
t109 = Ifges(7,5) * t186 + Ifges(7,6) * t185 - Ifges(7,3) * t334;
t101 = mrSges(5,1) * t221 - mrSges(5,3) * t134;
t100 = -mrSges(5,2) * t221 + mrSges(5,3) * t133;
t97 = -mrSges(6,1) * t161 + mrSges(6,2) * t160;
t96 = Ifges(6,1) * t160 + Ifges(6,4) * t161 + Ifges(6,5) * t376;
t95 = Ifges(6,4) * t160 + Ifges(6,2) * t161 + Ifges(6,6) * t376;
t93 = -mrSges(6,1) * t149 + mrSges(6,2) * t150;
t88 = -pkin(5) * t149 + t124;
t87 = Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t260;
t86 = Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t260;
t85 = Ifges(6,5) * t150 + Ifges(6,6) * t149 + Ifges(6,3) * t260;
t84 = -mrSges(5,1) * t133 + mrSges(5,2) * t134;
t83 = mrSges(7,1) * t260 - mrSges(7,3) * t92;
t82 = -mrSges(7,2) * t260 + mrSges(7,3) * t91;
t72 = -mrSges(7,2) * t376 + mrSges(7,3) * t77;
t71 = mrSges(7,1) * t376 - mrSges(7,3) * t76;
t70 = Ifges(5,1) * t134 + Ifges(5,4) * t133 + Ifges(5,5) * t221;
t69 = Ifges(5,4) * t134 + Ifges(5,2) * t133 + Ifges(5,6) * t221;
t65 = Ifges(7,1) * t105 + Ifges(7,4) * t106;
t64 = Ifges(7,4) * t105 + Ifges(7,2) * t106;
t62 = -mrSges(7,1) * t106 + mrSges(7,2) * t105;
t54 = -mrSges(6,2) * t221 + mrSges(6,3) * t58;
t53 = mrSges(6,1) * t221 - mrSges(6,3) * t57;
t47 = -mrSges(7,1) * t91 + mrSges(7,2) * t92;
t44 = Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(7,3) * t260;
t39 = -mrSges(7,1) * t77 + mrSges(7,2) * t76;
t38 = Ifges(7,1) * t76 + Ifges(7,4) * t77 + Ifges(7,5) * t376;
t37 = Ifges(7,4) * t76 + Ifges(7,2) * t77 + Ifges(7,6) * t376;
t32 = -pkin(5) * t58 + t73;
t27 = -mrSges(6,1) * t58 + mrSges(6,2) * t57;
t26 = Ifges(6,1) * t57 + Ifges(6,4) * t58 + Ifges(6,5) * t221;
t25 = Ifges(6,4) * t57 + Ifges(6,2) * t58 + Ifges(6,6) * t221;
t19 = -mrSges(7,2) * t221 + mrSges(7,3) * t23;
t18 = mrSges(7,1) * t221 - mrSges(7,3) * t22;
t9 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t8 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + Ifges(7,5) * t221;
t7 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + Ifges(7,6) * t221;
t1 = [-t342 * t70 + (t13 * t3 + t14 * t2 + t32 * t88) * t435 + (t11 * t41 + t12 * t40 + t124 * t73) * t436 + (t108 * t170 + t42 * t99 + t43 * t98) * t437 + 0.2e1 * m(3) * (t254 * t268 - t255 * t267) + 0.2e1 * m(4) * (t112 * t179 + t113 * t178 + t247 * t255) + 0.2e1 * t247 * t156 + 0.2e1 * t112 * t226 + 0.2e1 * t113 * t227 + t222 * t184 + t223 * t69 + 0.2e1 * t179 * t187 + 0.2e1 * t178 * t188 + 0.2e1 * t170 * t84 + 0.2e1 * t42 * t174 + 0.2e1 * t43 * t175 + 0.2e1 * t108 * t159 + t149 * t25 + t150 * t26 + t133 * t128 + t134 * t129 + 0.2e1 * t11 * t121 + 0.2e1 * t12 * t122 + 0.2e1 * t124 * t27 + 0.2e1 * t99 * t100 + 0.2e1 * t98 * t101 + t91 * t7 + t92 * t8 + 0.2e1 * t73 * t93 + 0.2e1 * t88 * t9 + t58 * t86 + t57 * t87 + 0.2e1 * t2 * t82 + 0.2e1 * t3 * t83 + 0.2e1 * t40 * t53 + 0.2e1 * t41 * t54 + t23 * t45 + t22 * t46 + 0.2e1 * t32 * t47 + 0.2e1 * t13 * t18 + 0.2e1 * t14 * t19 + (t127 + t85 + t44 - t183) * t221 + (t303 - 0.2e1 * t389 - 0.2e1 * t390) * t325 + (t138 + 0.2e1 * t387) * t261 + (-t137 + t24 + t6 + t68 + 0.2e1 * t388) * t260 + (-t335 * t364 + 0.2e1 * (t254 * t335 + t255 * t330) * mrSges(3,3) + ((t267 * t433 + Ifges(3,5) * t325 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t335) * t324) * t335 + (t268 * t433 + Ifges(4,5) * t261 - 0.2e1 * Ifges(3,6) * t325 - Ifges(4,6) * t260 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t330 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t335) * t324) * t330) * qJD(2)) * t324; t8 * t418 + t7 * t419 + t96 * t420 + t95 * t421 + t251 * t424 + t22 * t425 + t23 * t426 + t38 * t428 + t37 * t429 + t76 * t430 + t77 * t431 + t26 * t407 + t25 * t408 + t134 * t409 + t191 * t410 + t190 * t411 + m(6) * (t11 * t147 + t12 * t146 + t124 * t240 + t292 * t73 + t40 * t61 + t41 * t60) + m(7) * (t123 * t88 + t13 * t17 + t14 * t16 + t2 * t67 + t208 * t32 + t3 * t66) + t303 + (-t284 / 0.2e1 + t189 / 0.2e1 + t94 / 0.2e1 + t36 / 0.2e1) * t260 + t292 * t27 + t222 * t299 / 0.2e1 + t247 * t281 + t261 * t286 / 0.2e1 + t42 * t287 + t43 * t288 + t108 * t269 + t240 * t93 + t244 * t101 + t245 * t100 + t11 * t230 + t12 * t231 + t98 * t235 + t99 * t236 + t208 * t9 + t170 * t202 + t73 * t192 + t58 * t181 / 0.2e1 + t57 * t182 / 0.2e1 + t2 * t172 + t3 * t173 + t176 * t174 + t177 * t175 + t160 * t87 / 0.2e1 + t161 * t86 / 0.2e1 + (-t297 / 0.2e1 + t250 / 0.2e1 + t180 / 0.2e1 + t109 / 0.2e1) * t221 - pkin(2) * t156 + t146 * t53 + t147 * t54 + t40 * t144 + t41 * t145 + t60 * t121 + t61 * t122 + t123 * t47 + t124 * t97 + t32 * t117 + t88 * t39 + t16 * t82 + t17 * t83 + t13 * t71 + t14 * t72 + t66 * t18 + t67 * t19 + ((t184 / 0.2e1 + t128 * t403 - t178 * mrSges(4,3) + t356) * t334 + (t365 / 0.2e1 - t183 / 0.2e1 + t127 / 0.2e1 + t85 / 0.2e1 + t44 / 0.2e1 - t179 * mrSges(4,3)) * t329 + (-t329 * t226 + (t159 - t227) * t334 + t170 * t400 + m(4) * (-t178 * t334 - t179 * t329)) * pkin(9)) * qJD(3) + (t387 + t138 / 0.2e1 + t70 * t401 + t69 * t403 - t113 * mrSges(4,3) + Ifges(4,5) * t353 + (t128 * t402 + t129 * t403) * qJD(4) + (-t188 + t84) * pkin(9)) * t329 + m(5) * (t108 * t323 + t176 * t99 + t177 * t98 + t244 * t43 + t245 * t42) + m(4) * (-pkin(2) * t255 + t112 * t397 - t113 * t323) + (t137 / 0.2e1 - t68 / 0.2e1 - t24 / 0.2e1 - t6 / 0.2e1 - t388 + pkin(9) * t187 + t112 * mrSges(4,3) + Ifges(4,6) * t353) * t334 - t389 - t390 + (-t335 * t321 / 0.2e1 - Ifges(3,6) * t377) * t324; (t123 * t208 + t16 * t67 + t17 * t66) * t435 + (t146 * t61 + t147 * t60 + t240 * t292) * t436 + (t176 * t245 + t177 * t244) * t437 + 0.2e1 * t292 * t97 - 0.2e1 * pkin(2) * t281 + 0.2e1 * t176 * t287 + 0.2e1 * t177 * t288 - t256 * t95 - t257 * t96 + 0.2e1 * t240 * t192 + 0.2e1 * t244 * t235 + 0.2e1 * t245 * t236 + 0.2e1 * t60 * t230 + 0.2e1 * t61 * t231 + 0.2e1 * t208 * t39 + t186 * t38 + t161 * t181 + t160 * t182 + t185 * t37 + 0.2e1 * t16 * t172 + 0.2e1 * t17 * t173 + 0.2e1 * t146 * t144 + 0.2e1 * t147 * t145 + 0.2e1 * t123 * t117 + t77 * t110 + t76 * t111 + 0.2e1 * t66 * t71 + 0.2e1 * t67 * t72 + (-t189 + t284 - t36 - t94 + (-t328 * t251 + t333 * t252 + t269 * t434 + t299) * qJD(3)) * t334 + (t202 * t434 - t328 * t190 + t333 * t191 + t286 + (-t251 * t333 - t252 * t328) * qJD(4) + (0.2e1 * pkin(9) ^ 2 * t400 + t109 + t180 + t250 - t297) * qJD(3)) * t329; (-t11 * t346 - t12 * t276 + t215 * t40 - t216 * t41) * mrSges(6,3) + t439 * t108 + t364 + (-t105 * t13 + t106 * t14 + t2 * t203 - t204 * t3) * mrSges(7,3) + t87 * t415 + t8 * t416 + t7 * t417 + t154 * t420 + t153 * t421 + t22 * t422 + t23 * t423 + t296 * t424 + t65 * t428 + t64 * t429 + t105 * t430 + t106 * t431 + t69 * t401 + t134 * t404 + t26 * t405 + t25 * t406 + t285 * t410 + t283 * t411 + t57 * t412 + t58 * t413 + t86 * t414 + t316 * t27 + t170 * t280 + t246 * t9 + t228 * t53 + t229 * t54 + t73 * t217 + t195 * t47 + t164 * t121 + t165 * t122 + t124 * t151 + t32 * t139 + t118 * t18 + t119 * t19 - t112 * mrSges(4,2) + t113 * mrSges(4,1) + t88 * t62 - pkin(3) * t84 + t49 * t82 + t50 * t83 + ((-t328 * t99 - t333 * t98) * qJD(4) + t347) * mrSges(5,3) + t351 * t221 + t352 * t260 + t70 * t440 + (t356 + (pkin(4) * t93 - t128 / 0.2e1) * t328) * qJD(4) + m(7) * (t118 * t3 + t119 * t2 + t13 * t50 + t14 * t49 + t195 * t88 + t246 * t32) + m(6) * (t11 * t229 + t12 * t228 + t124 * t366 + t164 * t41 + t165 * t40 + t316 * t73) + (-t174 * t374 + m(5) * (-t372 * t98 - t374 * t99 + t347) + t333 * t100 - t328 * t101 - t175 * t372) * pkin(10); (t146 * t215 - t147 * t216 - t276 * t61 - t346 * t60) * mrSges(6,3) + (-t105 * t66 + t106 * t67 + t16 * t203 - t17 * t204) * mrSges(7,3) + ((-mrSges(4,1) + t439) * qJD(3) * pkin(9) - t352) * t334 + t38 * t416 + t37 * t417 + t65 * t418 + t64 * t419 + t76 * t422 + t77 * t423 + t105 * t425 + t106 * t426 + t96 * t405 + t95 * t406 + t154 * t407 + t153 * t408 + t160 * t412 + t161 * t413 + t181 * t414 + t182 * t415 + (t191 / 0.2e1 - t177 * mrSges(5,3) - t296 * t375 / 0.2e1 + (-t245 * mrSges(5,3) - t251 / 0.2e1 + (m(6) * t292 + t192) * pkin(4)) * qJD(4) + (-qJD(4) * t287 + m(5) * (-t177 - t441) - t235) * pkin(10)) * t328 + t321 + t316 * t97 + t292 * t151 + t240 * t217 + t246 * t39 + t228 * t144 + t229 * t145 + t164 * t230 + t165 * t231 + t208 * t62 - pkin(3) * t202 + t195 * t117 + t49 * t172 + t50 * t173 + t123 * t139 + t118 * t71 + t119 * t72 + m(7) * (t118 * t17 + t119 * t16 + t123 * t246 + t195 * t208 + t49 * t67 + t50 * t66) + (t285 * t401 + t283 * t403 + pkin(9) * t280 + (t296 * t402 + t298 * t403) * qJD(4) + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t351) * qJD(3)) * t329 + (t190 / 0.2e1 + qJD(4) * t409 + t375 * t404 + t354 * mrSges(5,3) + (m(5) * t354 - qJD(4) * t288 + t236) * pkin(10)) * t333 + m(6) * (t146 * t165 + t147 * t164 + t228 * t61 + t229 * t60 + t240 * t316); -0.2e1 * pkin(3) * t280 + t105 * t142 + t106 * t141 + 0.2e1 * t195 * t139 + 0.2e1 * t316 * t151 - t346 * t153 + t276 * t154 + t203 * t64 + t204 * t65 - t215 * t220 - t216 * t219 + 0.2e1 * t246 * t62 + t333 * t283 + t328 * t285 + (t333 * t298 + (0.2e1 * pkin(4) * t217 - t296) * t328) * qJD(4) + (t164 * t229 + t165 * t228 + t316 * t366) * t436 + (t118 * t50 + t119 * t49 + t195 * t246) * t435 + 0.2e1 * (-t105 * t118 + t106 * t119 + t203 * t49 - t204 * t50) * mrSges(7,3) + 0.2e1 * (-t164 * t346 - t165 * t276 + t215 * t228 - t216 * t229) * mrSges(6,3); t337 + m(7) * (t13 * t199 + t14 * t198 + t2 * t259 + t258 * t3) + (-t122 * t371 + m(6) * (t11 * t327 + t12 * t332 + t370 * t41 - t371 * t40) + t332 * t53 + t327 * t54 + t121 * t370) * pkin(4) + t258 * t18 + t259 * t19 + t199 * t83 + t198 * t82 + t43 * mrSges(5,1) - t42 * mrSges(5,2) + t68; t336 + m(7) * (t16 * t259 + t17 * t258 + t198 * t67 + t199 * t66) + (-t231 * t371 + m(6) * (-t146 * t371 + t147 * t370 + t327 * t60 + t332 * t61) + t332 * t144 + t327 * t145 + t230 * t370) * pkin(4) + t189 + t258 * t71 + t259 * t72 + t199 * t173 + t198 * t172 - t176 * mrSges(5,2) + t177 * mrSges(5,1); m(7) * (t118 * t199 + t119 * t198 + t258 * t50 + t259 * t49) + t320 + (pkin(10) * t294 - t392) * qJD(4) + (-t105 * t258 + t106 * t259 + t198 * t203 - t199 * t204) * mrSges(7,3) + (m(6) * (t164 * t327 + t165 * t332 + (-t228 * t327 + t229 * t332) * qJD(5)) + (t332 * t215 - t327 * t216 + (t276 * t327 - t332 * t346) * qJD(5)) * mrSges(6,3)) * pkin(4) + t338; -0.2e1 * t391 + 0.2e1 * t196 + (t198 * t259 + t199 * t258) * t435 + 0.2e1 * t341; (m(7) * (-t13 * t369 + t14 * t368 + t2 * t326 + t3 * t331) + t82 * t368 + t326 * t19 - t83 * t369 + t331 * t18) * pkin(5) + t337; (-t173 * t369 + t331 * t71 + m(7) * (t16 * t326 + t17 * t331 + t368 * t67 - t369 * t66) + t172 * t368 + t326 * t72) * pkin(5) + t336; (m(7) * (t326 * t49 + t331 * t50 + (-t118 * t326 + t119 * t331) * qJD(6)) + (-t331 * t105 + t326 * t106 + (t203 * t331 + t204 * t326) * qJD(6)) * mrSges(7,3)) * pkin(5) + t338; t341 + (m(7) * (t198 * t326 + t199 * t331 - t258 * t369 + t259 * t368) - mrSges(7,2) * t368 - mrSges(7,1) * t369) * pkin(5) + t355; 0.2e1 * t271; t344; t343; t345; t355; t271; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
