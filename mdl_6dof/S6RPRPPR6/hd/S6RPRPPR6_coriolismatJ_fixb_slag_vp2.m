% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:20
% EndTime: 2019-03-09 02:53:30
% DurationCPUTime: 5.33s
% Computational Cost: add. (14809->436), mult. (27747->611), div. (0->0), fcn. (31348->8), ass. (0->233)
t254 = sin(pkin(9));
t256 = sin(qJ(3));
t340 = cos(pkin(9));
t370 = cos(qJ(3));
t233 = t254 * t370 + t256 * t340;
t253 = sin(pkin(10));
t255 = cos(pkin(10));
t368 = sin(qJ(6));
t369 = cos(qJ(6));
t401 = t368 * t253 - t369 * t255;
t409 = t401 * t233;
t382 = t409 / 0.2e1;
t232 = t254 * t256 - t340 * t370;
t410 = t401 * t232;
t381 = t410 / 0.2e1;
t323 = t253 * t233;
t257 = -pkin(1) - pkin(7);
t319 = t256 * t257;
t238 = -qJ(4) * t256 + t319;
t303 = t370 * t257;
t239 = -qJ(4) * t370 + t303;
t402 = t340 * t238 + t254 * t239;
t141 = -pkin(5) * t323 + t402;
t393 = -m(7) / 0.2e1;
t395 = -m(6) / 0.2e1;
t356 = t409 * mrSges(7,2);
t276 = t253 * t369 + t255 * t368;
t405 = t276 * t233;
t359 = t405 * mrSges(7,1);
t404 = t359 / 0.2e1 - t356 / 0.2e1;
t412 = t141 * t393 + t395 * t402 + t404;
t411 = t276 * t409 - t401 * t405;
t248 = t256 * pkin(3) + qJ(2);
t408 = m(5) * t248;
t272 = t276 * t232;
t383 = t272 / 0.2e1;
t252 = t255 ^ 2;
t312 = t253 ^ 2 + t252;
t406 = mrSges(6,3) * t312;
t400 = t232 * t254 + t340 * t233;
t306 = t370 * pkin(3);
t184 = -t232 * pkin(4) + t233 * qJ(5) + t306;
t199 = t238 * t254 - t340 * t239;
t103 = t255 * t184 + t199 * t253;
t104 = t253 * t184 - t199 * t255;
t289 = -t103 * t253 + t104 * t255;
t241 = -mrSges(6,1) * t255 + mrSges(6,2) * t253;
t247 = -pkin(3) * t340 - pkin(4);
t399 = m(6) * t247 + t241;
t206 = t232 ^ 2;
t398 = t233 ^ 2;
t397 = 0.2e1 * t233;
t396 = m(5) / 0.2e1;
t394 = m(6) / 0.2e1;
t392 = m(7) / 0.2e1;
t391 = m(5) * pkin(3);
t390 = mrSges(7,1) / 0.2e1;
t389 = -mrSges(7,2) / 0.2e1;
t387 = m(7) * (-t272 * t409 + t405 * t410);
t385 = -t405 / 0.2e1;
t384 = t405 / 0.2e1;
t242 = pkin(3) * t254 + qJ(5);
t367 = pkin(8) + t242;
t221 = t367 * t253;
t222 = t367 * t255;
t180 = -t221 * t368 + t222 * t369;
t380 = -t180 / 0.2e1;
t361 = Ifges(7,4) * t276;
t196 = -Ifges(7,2) * t401 + t361;
t379 = t196 / 0.2e1;
t228 = Ifges(7,4) * t401;
t198 = Ifges(7,1) * t276 - t228;
t378 = t198 / 0.2e1;
t377 = -t232 / 0.2e1;
t376 = t232 / 0.2e1;
t375 = -t401 / 0.2e1;
t374 = t276 / 0.2e1;
t373 = -t276 / 0.2e1;
t372 = t253 / 0.2e1;
t371 = t255 / 0.2e1;
t326 = t276 * t410;
t327 = t401 * t272;
t67 = t326 - t327;
t307 = t67 * t392;
t366 = qJD(3) * t307;
t365 = m(7) * qJD(4);
t364 = Ifges(6,4) * t253;
t363 = Ifges(6,4) * t255;
t362 = Ifges(7,4) * t410;
t183 = pkin(4) * t233 + qJ(5) * t232 + t248;
t101 = t255 * t183 - t253 * t402;
t102 = t253 * t183 + t255 * t402;
t116 = mrSges(7,2) * t232 + mrSges(7,3) * t405;
t117 = -mrSges(7,2) * t233 + mrSges(7,3) * t272;
t118 = -mrSges(7,1) * t232 - mrSges(7,3) * t409;
t119 = mrSges(7,1) * t233 - mrSges(7,3) * t410;
t346 = t253 * Ifges(6,2);
t138 = -Ifges(6,6) * t232 + (t346 - t363) * t233;
t139 = -Ifges(6,5) * t232 + (-t255 * Ifges(6,1) + t364) * t233;
t331 = t232 * t253;
t140 = -pkin(5) * t331 + t199;
t345 = t255 * mrSges(6,2);
t347 = t253 * mrSges(6,1);
t291 = -t345 - t347;
t181 = t291 * t233;
t182 = t291 * t232;
t186 = mrSges(6,2) * t232 + mrSges(6,3) * t323;
t187 = -t233 * mrSges(6,2) + mrSges(6,3) * t331;
t328 = t233 * t255;
t188 = -mrSges(6,1) * t232 + mrSges(6,3) * t328;
t329 = t232 * t255;
t189 = t233 * mrSges(6,1) + mrSges(6,3) * t329;
t286 = Ifges(7,5) * t382 + Ifges(7,6) * t384;
t294 = t233 * mrSges(5,1) - t232 * mrSges(5,2);
t295 = -t232 * mrSges(5,1) - t233 * mrSges(5,2);
t68 = pkin(5) * t233 + pkin(8) * t329 + t101;
t89 = pkin(8) * t331 + t102;
t42 = -t368 * t89 + t369 * t68;
t43 = t368 * t68 + t369 * t89;
t69 = -pkin(5) * t232 + pkin(8) * t328 + t103;
t90 = pkin(8) * t323 + t104;
t44 = -t368 * t90 + t369 * t69;
t45 = t368 * t69 + t369 * t90;
t81 = Ifges(7,4) * t409 + Ifges(7,2) * t405 - Ifges(7,6) * t232;
t82 = Ifges(7,2) * t272 + t233 * Ifges(7,6) + t362;
t83 = Ifges(7,1) * t409 + Ifges(7,4) * t405 - Ifges(7,5) * t232;
t159 = Ifges(7,4) * t272;
t84 = Ifges(7,1) * t410 + t233 * Ifges(7,5) + t159;
t92 = t356 - t359;
t354 = t410 * mrSges(7,2);
t358 = t272 * mrSges(7,1);
t93 = t354 - t358;
t1 = m(7) * (t140 * t141 + t42 * t44 + t43 * t45) + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t256) * t256 + (qJ(2) * mrSges(4,1) + (-Ifges(4,1) + Ifges(4,2)) * t256 - Ifges(4,4) * t370) * t370 + m(6) * (t101 * t103 + t102 * t104 + t199 * t402) + (-Ifges(6,5) * t329 + Ifges(7,5) * t410 + Ifges(6,6) * t331 + Ifges(7,6) * t272) * t377 - Ifges(5,4) * t206 - t139 * t329 / 0.2e1 + t138 * t331 / 0.2e1 + (t408 + t294) * t306 + t199 * t181 + t102 * t186 + t104 * t187 + t101 * t188 + t103 * t189 + t140 * t92 + t141 * t93 + t83 * t381 + t84 * t382 + t81 * t383 + t82 * t384 + t402 * t182 + ((t252 * Ifges(6,1) / 0.2e1 - Ifges(6,3) - Ifges(7,3) + Ifges(5,1) - Ifges(5,2) + (-t363 + t346 / 0.2e1) * t253) * t232 + t286 + (-Ifges(6,5) * t255 + Ifges(6,6) * t253 + Ifges(5,4)) * t233) * t233 + t248 * t295 + t43 * t116 + t45 * t117 + t42 * t118 + t44 * t119;
t360 = t1 * qJD(1);
t357 = t272 * mrSges(7,2);
t355 = t410 * mrSges(7,1);
t351 = t401 * mrSges(7,1);
t350 = t401 * mrSges(7,3);
t349 = t276 * mrSges(7,2);
t348 = t276 * mrSges(7,3);
t316 = Ifges(7,5) * t272 - Ifges(7,6) * t410;
t91 = t355 + t357;
t94 = -Ifges(7,2) * t410 + t159;
t95 = Ifges(7,1) * t272 - t362;
t6 = t42 * t117 - t43 * t119 + t233 * t316 / 0.2e1 + t140 * t91 + (-t43 * mrSges(7,3) + t95 / 0.2e1 - t82 / 0.2e1) * t410 + (-t42 * mrSges(7,3) + t84 / 0.2e1 + t94 / 0.2e1) * t272;
t344 = t6 * qJD(1);
t335 = t409 * t410;
t336 = t405 * t272;
t266 = (t335 / 0.2e1 + t336 / 0.2e1) * mrSges(7,3) + t117 * t385 + t119 * t382 + t91 * t376;
t284 = t351 / 0.2e1 + t349 / 0.2e1;
t7 = t266 + t284;
t343 = t7 * qJD(1);
t290 = t101 * t253 - t102 * t255;
t320 = t255 * t187;
t324 = t253 * t189;
t333 = t199 * t232;
t9 = t409 * t117 + t405 * t119 + (mrSges(5,3) * t233 - t320 + t324) * t233 + (t232 * mrSges(5,3) - t182 - t93) * t232 + m(7) * (-t140 * t232 + t405 * t42 + t409 * t43) + m(6) * (t233 * t290 - t333) + m(5) * (-t233 * t402 - t333);
t342 = t9 * qJD(1);
t265 = (-t405 ^ 2 - t409 ^ 2 - t206) * t392 + (-t312 * t398 - t206) * t394 + (-t206 - t398) * t396;
t297 = m(6) * t312;
t298 = t276 ^ 2 + t401 ^ 2;
t270 = -m(5) / 0.2e1 - t297 / 0.2e1 + t298 * t393;
t29 = t265 + t270;
t339 = qJD(1) * t29;
t268 = t253 * t187 + t255 * t189 + m(6) * (t101 * t255 + t102 * t253);
t271 = t256 * mrSges(4,1) + mrSges(4,2) * t370 + t294;
t21 = t276 * t117 - t401 * t119 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t276 * t43 - t401 * t42) + t408 + t268 + t271;
t332 = t21 * qJD(1);
t205 = t232 * t233;
t325 = t253 * t188;
t321 = t255 * t186;
t275 = (t272 * t276 + t401 * t410) * t392;
t310 = m(6) / 0.4e1 + m(7) / 0.4e1;
t50 = t275 + 0.2e1 * (t297 / 0.4e1 + t310) * t232;
t318 = t50 * qJD(1);
t313 = -Ifges(7,5) * t401 - Ifges(7,6) * t276;
t193 = mrSges(7,1) * t276 - mrSges(7,2) * t401;
t311 = t193 * qJD(6);
t309 = -t387 / 0.2e1;
t308 = t387 / 0.2e1;
t300 = t193 * t376;
t299 = -t328 / 0.2e1;
t293 = t312 * t242;
t292 = t297 / 0.2e1;
t54 = 0.2e1 * t381 * mrSges(7,1) + 0.2e1 * t383 * mrSges(7,2);
t288 = qJD(1) * t54 + qJD(3) * t193;
t285 = t358 / 0.2e1 - t354 / 0.2e1;
t283 = t117 * t375 + t119 * t373;
t282 = t324 / 0.2e1 - t320 / 0.2e1;
t281 = m(7) * (t276 * t405 + t401 * t409);
t179 = -t221 * t369 - t222 * t368;
t280 = m(7) * (-t179 * t401 + t180 * t276);
t194 = t349 + t351;
t240 = -t255 * pkin(5) + t247;
t260 = (-t232 * t247 - t233 * t293) * t394 + (t179 * t405 + t180 * t409 - t232 * t240) * t392 + (t232 * t340 - t233 * t254) * t391 / 0.2e1 - t348 * t384 - t350 * t382 + (t194 + t241) * t377 - t233 * t406 / 0.2e1;
t263 = (t103 * t255 + t104 * t253) * t394 + (t276 * t45 - t401 * t44) * t392 + t118 * t375 + t116 * t374 + t186 * t372 + t188 * t371 + t306 * t396;
t10 = t260 - t263 - t295;
t279 = -t10 * qJD(1) + qJD(2) * t307;
t264 = (-t326 / 0.2e1 + t327 / 0.2e1) * mrSges(7,3) + t283;
t14 = t264 - t404;
t278 = t14 * qJD(1);
t19 = t272 * t117 - t410 * t119 + m(7) * (t272 * t43 - t410 * t42) + t268 * t232;
t277 = -t19 * qJD(1) + qJD(2) * t309;
t262 = (t272 * t375 - t373 * t410) * mrSges(7,3) + t290 * t395 + (-t179 * t410 + t180 * t272 - t276 * t42 - t401 * t43) * t392 - t282 + t283;
t13 = (t345 / 0.2e1 + t347 / 0.2e1) * t233 + t262 + t412;
t47 = m(7) * (-t179 * t276 - t180 * t401) + mrSges(7,3) * t298 + t242 * t297 + t406;
t51 = -t281 / 0.2e1 + (-t297 / 0.4e1 + t310) * t397;
t274 = t13 * qJD(1) - t51 * qJD(2) + t47 * qJD(3);
t30 = m(6) * (-t312 + 0.1e1) * t205 + m(7) * (t205 - t335 - t336);
t259 = (t92 / 0.2e1 + t181 / 0.2e1 + t282) * t232 + (t93 / 0.2e1 - t325 / 0.2e1 + t321 / 0.2e1 + t182 / 0.2e1) * t233 + ((t199 + t289) * t233 + (t402 + t290) * t232) * t394 + (t140 * t233 + t141 * t232 + t272 * t42 - t405 * t44 - t409 * t45 + t410 * t43) * t392 + t118 * t385 + t119 * t383 - t409 * t116 / 0.2e1 + t117 * t381;
t5 = -t280 / 0.2e1 + t259;
t273 = -t5 * qJD(1) - t30 * qJD(2) - t365 * t67 / 0.2e1;
t195 = -Ifges(7,2) * t276 - t228;
t197 = -Ifges(7,1) * t401 - t361;
t261 = -(t84 / 0.4e1 + t94 / 0.4e1) * t401 - (-t95 / 0.4e1 + t82 / 0.4e1) * t276 + (t197 / 0.4e1 - t196 / 0.4e1 + mrSges(7,3) * t380) * t410 + (t198 / 0.4e1 + t195 / 0.4e1 - t179 * mrSges(7,3) / 0.2e1) * t272 + t140 * t193 / 0.2e1 + t179 * t117 / 0.2e1 + t119 * t380 + t233 * t313 / 0.4e1 + t240 * t91 / 0.2e1;
t267 = Ifges(7,3) * t377 + t389 * t45 + t390 * t44 + t286;
t2 = t261 - t267;
t23 = t300 - t285;
t25 = t240 * t193 - (-t197 / 0.2e1 + t379) * t276 - (t378 + t195 / 0.2e1) * t401;
t269 = -t2 * qJD(1) - t23 * qJD(2) - t25 * qJD(3);
t56 = qJD(5) * t308;
t55 = t357 / 0.2e1 + t355 / 0.2e1 + t272 * t389 - t410 * t390;
t52 = t233 * t292 + t281 / 0.2e1 + t310 * t397;
t49 = t232 * t292 + t275 + (m(6) + m(7)) * t377;
t28 = t265 - t270;
t24 = t300 + t285;
t15 = t264 + t404;
t12 = mrSges(6,2) * t299 - mrSges(6,1) * t323 / 0.2e1 + t262 - t412;
t11 = t260 + t263;
t8 = t266 - t284;
t4 = t280 / 0.2e1 + t259;
t3 = t261 + t267;
t16 = [qJD(2) * t21 + qJD(3) * t1 + qJD(4) * t9 + qJD(5) * t19 + qJD(6) * t6, -m(7) * t411 * qJD(2) + t4 * qJD(3) + t28 * qJD(4) + t8 * qJD(6) + t332 + t56, t360 + t4 * qJD(2) + ((m(6) * t289 + t321 - t325) * t242 + t289 * mrSges(6,3) + m(7) * (t141 * t240 + t179 * t44 + t180 * t45) + (Ifges(6,5) * t253 + Ifges(7,5) * t276 + Ifges(6,6) * t255 - Ifges(7,6) * t401) * t377 + t409 * t378 - Ifges(4,6) * t370 + (Ifges(6,1) * t253 + t363) * t299 + (Ifges(6,2) * t255 + t364) * t323 / 0.2e1 - mrSges(4,1) * t319 - Ifges(4,5) * t256 + t247 * t181 + t240 * t92 - Ifges(5,5) * t233 + Ifges(5,6) * t232 + t141 * t194 + t179 * t118 + t180 * t116 - mrSges(4,2) * t303 + t400 * mrSges(5,3) * pkin(3) - t44 * t348 - t45 * t350 + t138 * t371 + t139 * t372 + t83 * t374 + t81 * t375 + t405 * t379 + (-t340 * t391 - mrSges(5,1) + t399) * t402 - (t254 * t391 - mrSges(5,2)) * t199) * qJD(3) + t11 * qJD(4) + t12 * qJD(5) + t3 * qJD(6), t28 * qJD(2) + t11 * qJD(3) + t49 * qJD(5) + t15 * qJD(6) + t411 * t365 + t342, t12 * qJD(3) + t49 * qJD(4) + t55 * qJD(6) - t277, t344 + t8 * qJD(2) + t3 * qJD(3) + t15 * qJD(4) + t55 * qJD(5) + (-mrSges(7,1) * t43 - mrSges(7,2) * t42 + t316) * qJD(6); qJD(3) * t5 + qJD(4) * t29 + qJD(6) * t7 - t332 + t56, t30 * qJD(3) (m(7) * (t179 * t272 + t180 * t410) - t410 * t350 - t272 * t348 - t400 * t391 - t271 + (m(7) * t240 + t194 + t399) * t233 + (-m(6) * t293 - t406) * t232) * qJD(3) + t52 * qJD(5) + t24 * qJD(6) - t273, t339 + t366, qJD(1) * t308 + t52 * qJD(3), t343 + t24 * qJD(3) + (mrSges(7,1) * t409 + mrSges(7,2) * t405) * qJD(6); -qJD(2) * t5 + qJD(4) * t10 + qJD(5) * t13 + qJD(6) * t2 - t360, -t51 * qJD(5) + t23 * qJD(6) + t273, qJD(5) * t47 + qJD(6) * t25, -t279, t274 (-mrSges(7,1) * t180 - mrSges(7,2) * t179 + t313) * qJD(6) - t269; -qJD(2) * t29 - qJD(3) * t10 + qJD(5) * t50 + qJD(6) * t14 - t342, -t339 + t366, t279, 0, t318, t278 - t311; -t13 * qJD(3) - t50 * qJD(4) + t54 * qJD(6) + t277, qJD(1) * t309 + t51 * qJD(3), -t274 + t311, -t318, 0, t288; -qJD(2) * t7 - qJD(3) * t2 - qJD(4) * t14 - qJD(5) * t54 - t344, -t23 * qJD(3) - t343, -t193 * qJD(5) + t269, -t278, -t288, 0;];
Cq  = t16;
