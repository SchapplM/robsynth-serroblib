% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:58
% EndTime: 2019-03-09 03:31:27
% DurationCPUTime: 20.54s
% Computational Cost: add. (3993->604), mult. (7383->739), div. (0->0), fcn. (3799->6), ass. (0->278)
t422 = Ifges(6,1) + Ifges(7,1);
t411 = Ifges(6,5) + Ifges(7,4);
t420 = Ifges(6,6) - Ifges(7,6);
t170 = cos(qJ(3));
t423 = t170 / 0.2e1;
t421 = Ifges(5,4) - Ifges(4,5);
t417 = -Ifges(4,6) + Ifges(5,5);
t419 = Ifges(6,3) + Ifges(7,2);
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t403 = t166 * t411 + t169 * t420;
t323 = Ifges(7,5) * t169;
t326 = Ifges(6,4) * t169;
t399 = t166 * t422 - t323 + t326;
t167 = sin(qJ(3));
t233 = mrSges(4,1) * t170 - mrSges(4,2) * t167;
t329 = Ifges(4,4) * t170;
t418 = (-Ifges(4,1) * t167 - t329) * t423 + qJ(2) * t233;
t276 = qJD(1) * qJD(3);
t117 = qJDD(1) * t167 + t170 * t276;
t289 = qJD(3) * t166;
t291 = qJD(1) * t167;
t193 = t169 * t291 - t289;
t49 = qJD(5) * t193 + qJDD(3) * t169 + t117 * t166;
t363 = t49 / 0.2e1;
t287 = qJD(3) * t169;
t109 = t166 * t291 + t287;
t50 = qJD(5) * t109 + qJDD(3) * t166 - t169 * t117;
t362 = -t50 / 0.2e1;
t116 = -t170 * qJDD(1) + t167 * t276;
t107 = qJDD(5) - t116;
t358 = t107 / 0.2e1;
t320 = t167 * mrSges(5,2);
t226 = -t170 * mrSges(5,3) - t320;
t113 = t226 * qJD(1);
t290 = qJD(1) * t170;
t137 = qJD(5) + t290;
t331 = mrSges(6,3) * t109;
t67 = mrSges(6,1) * t137 - t331;
t333 = mrSges(7,2) * t109;
t68 = -mrSges(7,1) * t137 + t333;
t335 = t67 - t68;
t332 = mrSges(6,3) * t193;
t66 = -mrSges(6,2) * t137 + t332;
t334 = mrSges(7,2) * t193;
t69 = mrSges(7,3) * t137 + t334;
t336 = t66 + t69;
t187 = t166 * t335 - t169 * t336;
t416 = t113 - t187;
t25 = -mrSges(6,2) * t107 - mrSges(6,3) * t50;
t26 = -mrSges(7,2) * t50 + mrSges(7,3) * t107;
t337 = t25 + t26;
t23 = mrSges(6,1) * t107 - mrSges(6,3) * t49;
t24 = -t107 * mrSges(7,1) + t49 * mrSges(7,2);
t338 = t23 - t24;
t83 = -t116 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t415 = t187 * qJD(5) - t337 * t166 - t338 * t169 - t83;
t360 = -m(3) - m(4);
t383 = m(6) + m(7);
t414 = -t109 / 0.2e1;
t413 = -t137 / 0.2e1;
t382 = mrSges(6,1) + mrSges(7,1);
t381 = mrSges(6,2) - mrSges(7,3);
t412 = -mrSges(6,3) - mrSges(7,2);
t103 = Ifges(6,4) * t193;
t325 = Ifges(7,5) * t193;
t379 = t109 * t422 + t411 * t137 + t103 - t325;
t313 = qJDD(3) / 0.2e1;
t321 = Ifges(5,6) * t170;
t209 = t167 * Ifges(5,3) - t321;
t328 = Ifges(6,4) * t109;
t39 = Ifges(6,2) * t193 + t137 * Ifges(6,6) + t328;
t410 = Ifges(5,5) * qJD(3) + qJD(1) * t209 + t169 * t39;
t220 = -t167 * Ifges(4,2) + t329;
t102 = Ifges(7,5) * t109;
t36 = t137 * Ifges(7,6) - Ifges(7,3) * t193 + t102;
t409 = Ifges(4,6) * qJD(3) + qJD(1) * t220 + t169 * t36;
t173 = -pkin(1) - pkin(7);
t129 = qJDD(1) * t173 + qJDD(2);
t132 = qJD(1) * t173 + qJD(2);
t288 = qJD(3) * t167;
t61 = t129 * t170 - t132 * t288;
t202 = qJDD(4) - t61;
t55 = -qJDD(3) * pkin(3) + t202;
t118 = t167 * t132;
t164 = qJD(3) * qJ(4);
t91 = -t118 - t164;
t408 = -qJD(3) * t91 - t55;
t124 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t291;
t265 = mrSges(5,1) * t291;
t126 = -qJD(3) * mrSges(5,3) + t265;
t407 = -t126 + t124;
t264 = mrSges(5,1) * t290;
t375 = -mrSges(4,3) * t290 - t264 + (mrSges(4,1) - mrSges(5,2)) * qJD(3);
t406 = -t167 * t419 + t170 * t403;
t405 = -t167 * t411 + t170 * t399;
t330 = Ifges(4,4) * t167;
t404 = t167 * (-Ifges(4,2) * t170 - t330) + t170 * (Ifges(5,2) * t167 + t321);
t402 = -t166 * t420 + t169 * t411;
t401 = t167 * t421 + t170 * t417;
t324 = Ifges(7,5) * t166;
t327 = Ifges(6,4) * t166;
t400 = t169 * t422 + t324 - t327;
t228 = t166 * mrSges(7,1) - t169 * mrSges(7,3);
t230 = mrSges(6,1) * t166 + mrSges(6,2) * t169;
t398 = t228 + t230;
t397 = t107 * t419 + t411 * t49 - t420 * t50;
t275 = qJDD(1) * qJ(2);
t277 = qJD(1) * qJD(2);
t133 = t275 + t277;
t286 = qJD(3) * t170;
t62 = t167 * t129 + t132 * t286;
t54 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t62;
t396 = -t167 * t54 - t170 * t55;
t395 = -t167 * t62 - t170 * t61;
t359 = pkin(3) + pkin(8);
t119 = t170 * t132;
t79 = -pkin(4) * t290 + t119;
t368 = -t79 + qJD(4);
t63 = -qJD(3) * t359 + t368;
t155 = t170 * qJ(4);
t295 = pkin(3) * t291 + qJD(1) * qJ(2);
t349 = pkin(8) * t167;
t70 = (-t155 + t349) * qJD(1) + t295;
t21 = -t166 * t70 + t169 * t63;
t22 = t166 * t63 + t169 * t70;
t241 = -qJD(4) * t170 + qJD(2);
t177 = qJ(4) * t116 + qJD(1) * t241 + t275;
t20 = t117 * t359 + t177;
t281 = qJD(5) * t169;
t283 = qJD(5) * t166;
t31 = -pkin(4) * t116 - qJDD(3) * t359 + t202;
t3 = t166 * t31 + t169 * t20 + t63 * t281 - t283 * t70;
t4 = -qJD(5) * t22 - t166 * t20 + t169 * t31;
t234 = t166 * t3 + t169 * t4;
t394 = t21 * t283 - t22 * t281 - t234;
t16 = -pkin(5) * t137 + qJD(6) - t21;
t17 = qJ(6) * t137 + t22;
t1 = qJ(6) * t107 + qJD(6) * t137 + t3;
t2 = -pkin(5) * t107 + qJDD(6) - t4;
t235 = t1 * t166 - t169 * t2;
t393 = -t16 * t283 - t17 * t281 - t235;
t392 = qJD(4) - t119;
t391 = -t49 * Ifges(7,5) / 0.2e1 - t107 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t363 + Ifges(6,6) * t358 + (Ifges(7,3) + Ifges(6,2)) * t362;
t225 = t170 * Ifges(4,1) - t330;
t390 = Ifges(4,5) * qJD(3) + qJD(1) * t225 + t109 * t411 + t137 * t419 + t193 * t420;
t389 = 0.2e1 * t313;
t388 = -mrSges(5,1) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t232 = mrSges(4,1) * t167 + mrSges(4,2) * t170;
t387 = t320 - (-m(5) * qJ(4) - mrSges(5,3)) * t170 + mrSges(2,2) - mrSges(3,3) - t232;
t386 = m(7) * pkin(5) + t382;
t385 = m(7) * qJ(6) - t381;
t384 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t339 = pkin(4) - t173;
t208 = pkin(5) * t169 + qJ(6) * t166;
t199 = -pkin(4) - t208;
t191 = t199 * t170;
t378 = -qJD(1) * t191 + qJD(5) * t208 - qJD(6) * t169 + t392;
t15 = mrSges(6,1) * t50 + mrSges(6,2) * t49;
t82 = mrSges(5,1) * t117 - qJDD(3) * mrSges(5,3);
t377 = -t82 + t15;
t58 = -mrSges(6,1) * t193 + mrSges(6,2) * t109;
t376 = -t126 + t58;
t168 = sin(qJ(1));
t302 = t168 * t170;
t374 = pkin(3) * t286 + qJ(4) * t288;
t207 = pkin(5) * t166 - qJ(6) * t169;
t373 = -m(7) * t207 - t398;
t171 = cos(qJ(1));
t347 = g(2) * t171;
t348 = g(1) * t168;
t369 = t347 - t348;
t367 = g(1) * t302;
t105 = t339 * t288;
t294 = t167 * pkin(3) - t155;
t236 = -t294 - t349;
t100 = qJ(2) - t236;
t123 = t339 * t170;
t312 = t169 * t100 + t166 * t123;
t65 = qJD(2) + (qJD(3) * pkin(8) - qJD(4)) * t170 + t374;
t13 = -qJD(5) * t312 - t105 * t169 - t166 * t65;
t229 = t169 * mrSges(7,1) + t166 * mrSges(7,3);
t231 = mrSges(6,1) * t169 - mrSges(6,2) * t166;
t78 = -pkin(4) * t291 + t118;
t72 = t164 + t78;
t27 = -pkin(5) * t193 - qJ(6) * t109 + t72;
t366 = t27 * t229 + t72 * t231;
t364 = qJD(1) ^ 2;
t361 = t50 / 0.2e1;
t357 = t193 / 0.2e1;
t356 = -t193 / 0.2e1;
t354 = t109 / 0.2e1;
t346 = g(3) * t167;
t112 = pkin(3) * t290 + qJ(4) * t291;
t77 = pkin(8) * t290 + t112;
t34 = t166 * t78 + t169 * t77;
t322 = Ifges(5,6) * t167;
t311 = qJD(3) * t27;
t310 = qJD(3) * t72;
t308 = t166 * t167;
t307 = t166 * t170;
t305 = t167 * t168;
t304 = t167 * t169;
t303 = t167 * t171;
t147 = t167 * t173;
t301 = t169 * t170;
t299 = t170 * t171;
t297 = pkin(3) * t302 + qJ(4) * t305;
t156 = t171 * qJ(2);
t296 = pkin(3) * t303 + t156;
t293 = t171 * pkin(1) + t168 * qJ(2);
t285 = qJD(3) * t173;
t282 = qJD(5) * t167;
t280 = qJD(5) * t359;
t279 = qJDD(1) * mrSges(3,2);
t278 = m(5) + t383;
t57 = -mrSges(7,1) * t193 - mrSges(7,3) * t109;
t270 = -t57 - t376;
t262 = t171 * pkin(7) + t293;
t136 = t170 * t285;
t253 = -t290 / 0.2e1;
t248 = -t283 / 0.2e1;
t247 = t281 / 0.2e1;
t246 = -t276 / 0.2e1;
t243 = (t133 + t277) * qJ(2);
t240 = pkin(3) * t305 + t262;
t227 = -t170 * mrSges(5,2) + t167 * mrSges(5,3);
t218 = -Ifges(6,2) * t166 + t326;
t217 = Ifges(6,2) * t169 + t327;
t211 = Ifges(7,3) * t166 + t323;
t210 = -Ifges(7,3) * t169 + t324;
t206 = t16 * t166 + t169 * t17;
t204 = t166 * t21 - t169 * t22;
t33 = -t166 * t77 + t169 * t78;
t85 = -qJD(3) * pkin(3) + t392;
t203 = t167 * t85 - t170 * t91;
t52 = -t100 * t166 + t123 * t169;
t120 = qJ(4) + t207;
t86 = -qJ(4) * t290 + t295;
t198 = t86 * t227;
t12 = -t100 * t283 - t166 * t105 + t123 * t281 + t169 * t65;
t32 = -pkin(4) * t117 - t54;
t190 = -t166 * t282 + t169 * t286;
t189 = t166 * t286 + t167 * t281;
t182 = -Ifges(6,6) * t167 + t170 * t217;
t181 = -Ifges(7,6) * t167 + t170 * t210;
t176 = qJD(5) * t206 + t235;
t175 = -qJD(5) * t204 + t234;
t151 = -qJDD(1) * pkin(1) + qJDD(2);
t146 = Ifges(5,6) * t291;
t122 = -pkin(4) * t167 + t147;
t121 = qJ(2) + t294;
t114 = t232 * qJD(1);
t106 = -pkin(4) * t286 + t136;
t101 = t231 * t167;
t96 = -t166 * t302 + t169 * t171;
t95 = t166 * t171 + t168 * t301;
t94 = t166 * t299 + t168 * t169;
t93 = t166 * t168 - t169 * t299;
t90 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t290 + t146;
t81 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t117;
t80 = qJDD(3) * mrSges(4,1) + mrSges(4,3) * t116;
t75 = t241 + t374;
t64 = t167 * t199 + t147;
t56 = pkin(5) * t109 - qJ(6) * t193;
t48 = -pkin(5) * t170 - t52;
t42 = qJ(6) * t170 + t312;
t35 = pkin(3) * t117 + t177;
t30 = pkin(5) * t291 - t33;
t29 = -qJ(6) * t291 + t34;
t19 = t136 + (qJD(5) * t207 - qJD(6) * t166) * t167 + qJD(3) * t191;
t14 = mrSges(7,1) * t50 - mrSges(7,3) * t49;
t11 = pkin(5) * t288 - t13;
t10 = t49 * Ifges(6,1) - t50 * Ifges(6,4) + t107 * Ifges(6,5);
t9 = t49 * Ifges(7,1) + t107 * Ifges(7,4) + t50 * Ifges(7,5);
t6 = -qJ(6) * t288 + qJD(6) * t170 + t12;
t5 = pkin(5) * t50 - qJ(6) * t49 - qJD(6) * t109 + t32;
t7 = [(-t83 + t80) * t170 * t173 + (-m(3) * t293 - m(4) * t262 - m(5) * t240 - t383 * (t171 * pkin(4) + pkin(8) * t305 - t155 * t168 + t240) - t386 * t96 - t385 * t95 + t412 * t305 + t388 * t171 + t387 * t168) * g(2) + (-m(5) * t296 - t383 * (pkin(8) * t303 - qJ(4) * t299 - t339 * t168 + t296) + t386 * t94 + t385 * t93 + t412 * t303 + t360 * t156 + t387 * t171 + (m(3) * pkin(1) + (-m(4) - m(5)) * t173 - t388) * t168) * g(1) - t409 * t286 / 0.2e1 + (qJD(3) * t406 + t282 * t402) * t137 / 0.2e1 + t407 * t136 + m(5) * (t121 * t35 + t75 * t86 + (qJD(3) * t203 + t396) * t173) + t401 * qJD(3) ^ 2 / 0.2e1 + t404 * t246 + (qJD(3) * t405 + t282 * t400) * t354 + t395 * mrSges(4,3) + m(4) * (-t173 * t395 + t243) + (-t21 * t189 + t22 * t190 + t3 * t304 - t4 * t308) * mrSges(6,3) + (-t82 + t81) * t147 + (qJ(2) * mrSges(4,1) - t121 * mrSges(5,2) - t220 / 0.2e1 + t209 / 0.2e1) * t117 + (0.2e1 * mrSges(3,3) + t232) * t133 - pkin(1) * t279 + (t286 * t91 - t396) * mrSges(5,1) + (-t85 * mrSges(5,1) - t390 / 0.2e1 - t17 * mrSges(7,3) - t21 * mrSges(6,1) + t22 * mrSges(6,2) + t16 * mrSges(7,1) + t90 / 0.2e1) * t288 + t151 * mrSges(3,2) + t122 * t15 + t75 * t113 + qJD(2) * t114 + t391 * t304 + t418 * t276 + t106 * t58 - t32 * t101 + t64 * t14 + t12 * t66 + t13 * t67 + t11 * t68 + t6 * t69 + t19 * t57 + t52 * t23 + t42 * t26 + t48 * t24 + (t1 * t304 + t16 * t189 + t17 * t190 + t2 * t308) * mrSges(7,2) + (qJD(5) * t36 + t10 + t9) * t308 / 0.2e1 + (-qJ(2) * mrSges(4,2) + t121 * mrSges(5,3) + t322 / 0.2e1 - t225 / 0.2e1) * t116 + t312 * t25 + m(6) * (t106 * t72 + t12 * t22 + t122 * t32 + t13 * t21 + t3 * t312 + t4 * t52) + m(7) * (t1 * t42 + t11 * t16 + t17 * t6 + t19 * t27 + t2 * t48 + t5 * t64) + (-Ifges(5,2) * t116 - Ifges(5,6) * t117 / 0.2e1 + Ifges(4,5) * t313 + Ifges(7,6) * t361 + Ifges(6,6) * t362 + t411 * t363 + t419 * t358 - t389 * Ifges(5,4) + t384) * t170 + (-Ifges(4,1) * t116 - Ifges(4,4) * t117 + Ifges(4,5) * qJDD(3) + t397) * t423 + m(3) * (-pkin(1) * t151 + t243) + t35 * t226 + t72 * (-mrSges(6,1) * t190 + mrSges(6,2) * t189) + t27 * (-mrSges(7,1) * t190 - mrSges(7,3) * t189) + (-t375 * t285 + (Ifges(5,3) * t170 + t322) * t276 / 0.2e1 + t379 * t247 + t210 * t361 + t217 * t362 - t5 * t229 + t39 * t248 + t399 * t363 + t403 * t358 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t117 + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t116 + t417 * t389) * t167 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (t166 * t379 + t410) * t286 / 0.2e1 + (qJD(3) * t181 + t211 * t282) * t356 + (qJD(3) * t182 + t218 * t282) * t357 + qJD(3) * t198; m(3) * t151 + t279 + (qJ(2) * t360 - mrSges(3,3)) * t364 + (-m(5) * t86 + m(6) * t204 - m(7) * t206 - t114 - t416) * qJD(1) + (t14 + t81 + (t166 * t336 + t169 * t335 - t375) * qJD(3) + m(6) * (t21 * t287 + t22 * t289 + t32) + m(7) * (-t16 * t287 + t17 * t289 + t5) + m(5) * (qJD(3) * t85 - t54) + m(4) * t62 + t377) * t167 + (t80 + (t124 - t270) * qJD(3) + m(6) * (t310 + t394) + m(7) * (t311 + t393) + m(5) * t408 + m(4) * t61 + t415) * t170 + t369 * (t278 - t360); (-(t217 / 0.2e1 - t210 / 0.2e1) * t193 + t403 * t413 + t399 * t414 + t366) * qJD(5) + t409 * t290 / 0.2e1 + t410 * t253 - (t405 * t109 + t406 * t137) * qJD(1) / 0.2e1 - t407 * t119 + t400 * t363 + t401 * t246 + t402 * t358 + (-t367 + t393) * mrSges(7,2) + (-t367 + t394) * mrSges(6,3) + t390 * t291 / 0.2e1 + (t379 * t253 + t335 * t280 - t337 * t359 - t391) * t166 - t39 * t281 / 0.2e1 + t120 * t14 - t112 * t113 - (Ifges(5,3) * t290 + t146 + t90) * t291 / 0.2e1 + (-m(5) * t297 - t383 * (pkin(8) * t302 + t297) + (t167 * t373 - t227) * t168) * g(1) + (-pkin(3) * t55 - qJ(4) * t54 - qJD(4) * t91 - t112 * t86 - t132 * t203) * m(5) + (t404 / 0.2e1 - t418) * t364 - t79 * t58 - pkin(3) * t83 - t34 * t66 - t33 * t67 - t30 * t68 - t29 * t69 + t55 * mrSges(5,2) + t61 * mrSges(4,1) - t62 * mrSges(4,2) - t54 * mrSges(5,3) - t233 * t348 + t379 * t248 + (-(t182 / 0.2e1 - t181 / 0.2e1) * t193 - t198 - t22 * (mrSges(6,2) * t167 + mrSges(6,3) * t301) - t17 * (mrSges(7,2) * t301 - mrSges(7,3) * t167) - t21 * (-mrSges(6,1) * t167 - mrSges(6,3) * t307) - t16 * (mrSges(7,1) * t167 + mrSges(7,2) * t307)) * qJD(1) + (-t336 * t280 - t338 * t359 + t10 / 0.2e1 + t9 / 0.2e1) * t169 + (t120 * t5 - t16 * t30 - t17 * t29 - t176 * t359 + t27 * t378) * m(7) + (qJ(4) * t32 - t175 * t359 - t21 * t33 - t22 * t34 + t368 * t72) * m(6) + t377 * qJ(4) + t378 * t57 + t376 * qJD(4) - t91 * t264 + t375 * t118 + t85 * t265 + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + t421 * t116 + t366 * t290 + t5 * t228 + t32 * t230 + t417 * t117 + (t227 + t233 + (m(5) * pkin(3) + t359 * t383 - t412) * t170 + (m(7) * t120 + (m(5) + m(6)) * qJ(4) + t398) * t167) * t347 + (m(5) * t294 - t167 * t412 + t373 * t170 - t236 * t383 + t226 + t232) * g(3) + t36 * t247 + t211 * t361 + t218 * t362; t270 * qJD(3) - t278 * t346 + (qJD(1) * t416 - t369 * t278) * t170 + (t206 * t290 + t176 - t311) * m(7) + (-t204 * t290 + t175 - t310) * m(6) + (t290 * t86 - t408) * m(5) - t415; (t381 * t96 + t382 * t95) * g(1) + (t381 * t94 + t382 * t93) * g(2) - g(3) * t101 + (-t27 * t56 - g(1) * (-pkin(5) * t95 + qJ(6) * t96) - g(2) * (-pkin(5) * t93 + qJ(6) * t94) - t208 * t346 - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t17) * m(7) + qJD(6) * t69 - t56 * t57 + qJ(6) * t26 - pkin(5) * t24 + (-m(7) * t17 + t332 - t336) * t21 - t72 * (mrSges(6,1) * t109 + mrSges(6,2) * t193) + (-t109 * t420 + t411 * t193) * t413 + (t193 * t422 + t102 - t328 + t36) * t414 + (-Ifges(6,2) * t109 + t103 + t379) * t356 - t27 * (mrSges(7,1) * t109 - mrSges(7,3) * t193) + (-m(7) * t16 + t331 + t335) * t22 + (Ifges(7,3) * t109 + t325) * t357 - t229 * t346 + t384 + t397 + t17 * t333 - t16 * t334 + t39 * t354; t109 * t57 - t137 * t69 + (-g(1) * t95 - g(2) * t93 + g(3) * t304 + t109 * t27 - t137 * t17 + t2) * m(7) + t24;];
tau  = t7;
