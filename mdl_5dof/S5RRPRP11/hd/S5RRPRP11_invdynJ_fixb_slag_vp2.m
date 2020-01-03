% Calculate vector of inverse dynamics joint torques for
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:23
% EndTime: 2019-12-31 20:12:49
% DurationCPUTime: 16.19s
% Computational Cost: add. (2915->528), mult. (6131->690), div. (0->0), fcn. (3331->6), ass. (0->249)
t401 = m(5) + m(6);
t397 = Ifges(5,1) + Ifges(6,1);
t378 = Ifges(6,4) + Ifges(5,5);
t400 = Ifges(5,6) - Ifges(6,6);
t380 = -mrSges(5,3) - mrSges(6,2);
t154 = sin(qJ(4));
t158 = cos(qJ(2));
t267 = qJD(1) * t158;
t242 = t154 * t267;
t157 = cos(qJ(4));
t265 = qJD(2) * t157;
t100 = -t242 + t265;
t399 = -t100 / 0.2e1;
t155 = sin(qJ(2));
t268 = qJD(1) * t155;
t128 = qJD(4) + t268;
t398 = -t128 / 0.2e1;
t379 = -Ifges(4,4) + Ifges(3,5);
t377 = Ifges(4,5) - Ifges(3,6);
t396 = Ifges(6,2) + Ifges(5,3);
t192 = pkin(4) * t154 - qJ(5) * t157;
t212 = t154 * mrSges(6,1) - t157 * mrSges(6,3);
t214 = mrSges(5,1) * t154 + mrSges(5,2) * t157;
t394 = -m(6) * t192 - t212 - t214;
t363 = t378 * t154 + t400 * t157;
t298 = Ifges(6,5) * t157;
t300 = Ifges(5,4) * t157;
t360 = t397 * t154 - t298 + t300;
t326 = pkin(2) + pkin(7);
t393 = t401 * t326 - t380;
t304 = mrSges(5,3) * t100;
t61 = mrSges(5,1) * t128 - t304;
t305 = mrSges(6,2) * t100;
t62 = -mrSges(6,1) * t128 + t305;
t374 = -t61 + t62;
t99 = qJD(2) * t154 + t157 * t267;
t59 = -mrSges(5,2) * t128 - mrSges(5,3) * t99;
t60 = -mrSges(6,2) * t99 + mrSges(6,3) * t128;
t375 = t59 + t60;
t173 = t374 * t154 + t375 * t157;
t257 = qJD(1) * qJD(2);
t108 = -t158 * qJDD(1) + t155 * t257;
t284 = qJD(4) * t99;
t45 = qJDD(2) * t157 + t108 * t154 - t284;
t109 = qJDD(1) * t155 + t158 * t257;
t98 = qJDD(4) + t109;
t19 = mrSges(5,1) * t98 - mrSges(5,3) * t45;
t20 = -t98 * mrSges(6,1) + t45 * mrSges(6,2);
t261 = qJD(4) * t157;
t46 = qJD(2) * t261 - qJD(4) * t242 + qJDD(2) * t154 - t157 * t108;
t21 = -mrSges(5,2) * t98 - mrSges(5,3) * t46;
t22 = -mrSges(6,2) * t46 + mrSges(6,3) * t98;
t392 = t173 * qJD(4) + (t21 + t22) * t154 + (t19 - t20) * t157;
t156 = sin(qJ(1));
t391 = g(2) * t156;
t211 = t158 * mrSges(4,2) - t155 * mrSges(4,3);
t217 = mrSges(3,1) * t158 - mrSges(3,2) * t155;
t390 = t211 - t217;
t297 = Ifges(4,6) * t155;
t194 = -t158 * Ifges(4,3) - t297;
t295 = t100 * Ifges(5,4);
t33 = -t99 * Ifges(5,2) + t128 * Ifges(5,6) + t295;
t318 = Ifges(6,5) * t99;
t94 = Ifges(5,4) * t99;
t376 = t100 * t397 + t128 * t378 + t318 - t94;
t388 = Ifges(4,5) * qJD(2) + qJD(1) * t194 + t154 * t376 + t157 * t33;
t134 = pkin(6) * t268;
t387 = qJD(3) + t134;
t332 = t45 / 0.2e1;
t330 = t46 / 0.2e1;
t329 = t98 / 0.2e1;
t385 = t378 * t98 + (-Ifges(5,4) + Ifges(6,5)) * t46 + t397 * t45;
t263 = qJD(3) * t155;
t283 = qJDD(1) * pkin(1);
t172 = -qJ(3) * t109 - qJD(1) * t263 - t283;
t24 = t108 * t326 + t172;
t262 = qJD(4) * t154;
t97 = t109 * pkin(6);
t228 = qJDD(3) + t97;
t49 = pkin(3) * t109 - qJDD(2) * t326 + t228;
t140 = t155 * qJ(3);
t227 = -pkin(1) - t140;
t66 = (-t158 * t326 + t227) * qJD(1);
t104 = -pkin(3) * t268 - t134;
t349 = -t104 + qJD(3);
t68 = -qJD(2) * t326 + t349;
t3 = t154 * t49 + t157 * t24 + t68 * t261 - t262 * t66;
t1 = qJ(5) * t98 + qJD(5) * t128 + t3;
t384 = t1 * mrSges(6,3);
t383 = t3 * mrSges(5,2);
t382 = mrSges(5,1) + mrSges(6,1);
t381 = mrSges(5,2) - mrSges(6,3);
t193 = pkin(4) * t157 + qJ(5) * t154;
t183 = -pkin(3) - t193;
t373 = qJD(4) * t193 - qJD(5) * t157 - t183 * t268 + t387;
t247 = mrSges(4,1) * t267;
t116 = -qJD(2) * mrSges(4,3) - t247;
t55 = mrSges(5,1) * t99 + mrSges(5,2) * t100;
t372 = -t116 + t55;
t325 = pkin(3) + pkin(6);
t119 = t325 * t155;
t146 = t158 * pkin(2);
t271 = t146 + t140;
t246 = t158 * pkin(7) + t271;
t92 = -pkin(1) - t246;
t371 = t154 * t119 + t157 * t92;
t303 = Ifges(3,4) * t155;
t206 = t158 * Ifges(3,2) + t303;
t93 = Ifges(6,5) * t100;
t30 = t128 * Ifges(6,6) + t99 * Ifges(6,3) + t93;
t370 = Ifges(3,6) * qJD(2) + qJD(1) * t206 + t157 * t30;
t368 = qJD(2) * mrSges(3,2) - mrSges(3,3) * t267 + t116;
t248 = mrSges(4,1) * t268;
t367 = mrSges(3,3) * t268 + t248 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2);
t366 = t155 * t363 + t158 * t396;
t365 = t155 * t360 + t158 * t378;
t296 = Ifges(4,6) * t158;
t364 = t158 * (Ifges(4,3) * t155 - t296) + t155 * (-Ifges(4,2) * t158 + t297);
t362 = -t154 * t400 + t157 * t378;
t361 = t155 * t377 + t158 * t379;
t299 = Ifges(6,5) * t154;
t301 = Ifges(5,4) * t154;
t359 = t157 * t397 + t299 - t301;
t25 = -t154 * t66 + t157 * t68;
t358 = -t25 + qJD(5);
t357 = t378 * t45 + t396 * t98 - t400 * t46;
t96 = t108 * pkin(6);
t356 = t155 * t97 - t158 * t96;
t67 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t96;
t70 = -qJDD(2) * pkin(2) + t228;
t355 = t155 * t70 - t158 * t67;
t16 = -pkin(4) * t128 + t358;
t354 = -t16 * mrSges(6,2) + t25 * mrSges(5,3);
t26 = t154 * t68 + t157 * t66;
t4 = -qJD(4) * t26 - t154 * t24 + t157 * t49;
t353 = -t154 * t3 - t157 * t4;
t2 = -pkin(4) * t98 + qJDD(5) - t4;
t352 = -t1 * t154 + t157 * t2;
t159 = cos(qJ(1));
t351 = g(1) * t159 + t391;
t350 = -t45 * Ifges(5,4) / 0.2e1 - t98 * Ifges(5,6) / 0.2e1 + Ifges(6,5) * t332 + Ifges(6,6) * t329 + (Ifges(5,2) + Ifges(6,3)) * t330;
t258 = m(4) + t401;
t348 = -mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t133 = Ifges(3,4) * t267;
t347 = Ifges(3,1) * t268 + Ifges(3,5) * qJD(2) + t100 * t378 + t128 * t396 - t400 * t99 + t133;
t346 = -mrSges(2,1) + t390;
t290 = t158 * mrSges(4,3);
t345 = -t290 + t394 * t158 + (m(4) * pkin(2) - mrSges(4,2) + t393) * t155;
t136 = pkin(6) * t267;
t105 = pkin(3) * t267 + t136;
t153 = qJD(2) * qJ(3);
t80 = t153 + t105;
t27 = pkin(4) * t99 - qJ(5) * t100 + t80;
t344 = -t27 * t100 - t2;
t264 = qJD(2) * t158;
t107 = t325 * t264;
t286 = qJ(3) * t158;
t191 = pkin(7) * t155 - t286;
t266 = qJD(2) * t155;
t224 = pkin(2) * t266 - t263;
t64 = qJD(2) * t191 + t224;
t15 = -qJD(4) * t371 + t107 * t157 - t154 * t64;
t338 = -m(6) * pkin(4) - t382;
t337 = -m(6) * qJ(5) + t381;
t331 = -t46 / 0.2e1;
t328 = -t99 / 0.2e1;
t327 = t99 / 0.2e1;
t323 = t100 / 0.2e1;
t317 = pkin(6) * t155;
t314 = g(3) * t158;
t144 = t158 * pkin(6);
t135 = pkin(2) * t268;
t71 = qJD(1) * t191 + t135;
t43 = t154 * t105 + t157 * t71;
t302 = Ifges(3,4) * t158;
t282 = t154 * t155;
t281 = t154 * t156;
t280 = t154 * t158;
t278 = t155 * t157;
t277 = t155 * t159;
t276 = t156 * t157;
t275 = t157 * t158;
t274 = t157 * t159;
t272 = t158 * t159;
t120 = t158 * pkin(3) + t144;
t269 = t159 * pkin(1) + t156 * pkin(6);
t260 = qJD(4) * t158;
t245 = t154 * t260;
t226 = -t257 / 0.2e1;
t74 = t109 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t223 = pkin(2) * t272 + qJ(3) * t277 + t269;
t216 = mrSges(3,1) * t155 + mrSges(3,2) * t158;
t215 = mrSges(5,1) * t157 - mrSges(5,2) * t154;
t213 = t157 * mrSges(6,1) + t154 * mrSges(6,3);
t204 = -Ifges(5,2) * t154 + t300;
t203 = Ifges(5,2) * t157 + t301;
t197 = Ifges(6,3) * t154 + t298;
t196 = -Ifges(6,3) * t157 + t299;
t195 = -t155 * Ifges(4,2) - t296;
t17 = qJ(5) * t128 + t26;
t190 = t16 * t154 + t17 * t157;
t188 = t25 * t154 - t26 * t157;
t42 = t105 * t157 - t154 * t71;
t51 = t119 * t157 - t154 * t92;
t110 = -qJD(2) * pkin(2) + t387;
t114 = -t136 - t153;
t185 = t110 * t158 + t114 * t155;
t184 = t227 - t146;
t181 = pkin(1) * t216;
t14 = t154 * t107 + t119 * t261 + t157 * t64 - t262 * t92;
t81 = t184 * qJD(1);
t180 = t81 * (-mrSges(4,2) * t155 - t290);
t179 = t155 * (Ifges(3,1) * t158 - t303);
t175 = t155 * t265 + t245;
t174 = t154 * t266 - t157 * t260;
t168 = Ifges(5,6) * t158 + t155 * t203;
t167 = Ifges(6,6) * t158 + t155 * t196;
t50 = -pkin(3) * t108 - t67;
t164 = qJD(4) * t190 - t352;
t163 = -qJD(4) * t188 - t353;
t147 = t159 * pkin(6);
t112 = -pkin(1) - t271;
t111 = qJ(3) + t192;
t106 = t325 * t266;
t103 = -qJ(3) * t267 + t135;
t102 = t211 * qJD(1);
t91 = t215 * t158;
t85 = -t155 * t281 + t274;
t84 = t154 * t159 + t155 * t276;
t83 = t154 * t277 + t276;
t82 = -t155 * t274 + t281;
t79 = Ifges(4,4) * qJD(2) + qJD(1) * t195;
t75 = -qJ(3) * t264 + t224;
t73 = mrSges(4,1) * t108 - qJDD(2) * mrSges(4,3);
t63 = t158 * t193 + t120;
t54 = mrSges(6,1) * t99 - mrSges(6,3) * t100;
t53 = pkin(4) * t100 + qJ(5) * t99;
t48 = -pkin(4) * t155 - t51;
t47 = qJ(5) * t155 + t371;
t44 = pkin(2) * t108 + t172;
t29 = -pkin(4) * t267 - t42;
t28 = qJ(5) * t267 + t43;
t23 = (-qJD(4) * t192 + qJD(5) * t154) * t158 + (-pkin(6) + t183) * t266;
t13 = mrSges(5,1) * t46 + mrSges(5,2) * t45;
t12 = mrSges(6,1) * t46 - mrSges(6,3) * t45;
t11 = -pkin(4) * t264 - t15;
t10 = qJ(5) * t264 + qJD(5) * t155 + t14;
t5 = pkin(4) * t46 - qJ(5) * t45 - qJD(5) * t100 + t50;
t6 = [t355 * mrSges(4,1) + (t388 / 0.2e1 + t368 * pkin(6) - t370 / 0.2e1 + t114 * mrSges(4,1)) * t266 + (-Ifges(3,4) * t108 + Ifges(3,5) * qJDD(2) + t357) * t155 / 0.2e1 + t155 * t384 + m(6) * (t1 * t47 + t10 * t17 + t11 * t16 + t2 * t48 + t23 * t27 + t5 * t63) - (t154 * t30 + t157 * t376) * t260 / 0.2e1 + (-t174 * t25 + t175 * t26 - t275 * t3 + t280 * t4) * mrSges(5,3) + t158 * (Ifges(3,4) * t109 - Ifges(3,2) * t108 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t158 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t109 + Ifges(4,3) * t108) / 0.2e1 - t108 * t206 / 0.2e1 + t44 * t211 + t108 * t194 / 0.2e1 - t109 * t195 / 0.2e1 - t155 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t109 + Ifges(4,6) * t108) / 0.2e1 + t120 * t13 + t364 * t226 + (qJD(2) * t365 - t260 * t359) * t323 + (qJD(2) * t366 - t260 * t362) * t128 / 0.2e1 - t2 * mrSges(6,1) * t155 + (t155 * t396 - t363 * t158) * t329 + t350 * t275 + (-m(3) * t269 - m(4) * t223 - t401 * (t156 * pkin(3) + pkin(7) * t272 + t223) + t338 * t83 + t337 * t82 + t380 * t272 + t346 * t159 + t348 * t156) * g(2) + (-t401 * (t159 * pkin(3) + t147) + t338 * t85 + t337 * t84 + (-m(3) - m(4)) * t147 + t348 * t159 + (m(3) * pkin(1) - m(4) * t184 + t158 * t393 - t227 * t401 - t346) * t156) * g(1) - pkin(1) * (mrSges(3,1) * t108 + mrSges(3,2) * t109) + t112 * (-mrSges(4,2) * t108 - mrSges(4,3) * t109) + t75 * t102 - t106 * t55 + t50 * t91 + t14 * t59 + t10 * t60 + t15 * t61 + t11 * t62 + t63 * t12 + t48 * t20 + t51 * t19 + t23 * t54 + t47 * t22 + (-qJDD(2) * mrSges(3,1) + t74) * t317 + m(4) * (t112 * t44 + t75 * t81 + (qJD(2) * t185 + t355) * pkin(6)) + (-t108 * t144 + t109 * t317 + t356) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t356) + t361 * qJD(2) ^ 2 / 0.2e1 + t33 * t245 / 0.2e1 + (-t1 * t275 + t16 * t174 + t17 * t175 - t2 * t280) * mrSges(6,2) - t155 * t383 - t385 * t280 / 0.2e1 + (t179 + t158 * (-Ifges(3,2) * t155 + t302)) * t257 / 0.2e1 + t109 * t302 / 0.2e1 + (t155 * t378 - t158 * t360) * t332 + (t155 * t379 - t158 * t377) * qJDD(2) / 0.2e1 + (qJD(2) * t167 - t197 * t260) * t327 + (qJD(2) * t168 - t204 * t260) * t328 + (Ifges(6,6) * t155 - t158 * t196) * t330 + (Ifges(5,6) * t155 - t158 * t203) * t331 - t181 * t257 + (-qJDD(2) * mrSges(3,2) - t73) * t144 + t5 * t213 * t158 + t4 * mrSges(5,1) * t155 + (t367 * pkin(6) + t110 * mrSges(4,1) - t79 / 0.2e1 - t26 * mrSges(5,2) - t16 * mrSges(6,1) + t17 * mrSges(6,3) + t25 * mrSges(5,1) + t347 / 0.2e1) * t264 + t371 * t21 + m(5) * (-t106 * t80 + t120 * t50 + t14 * t26 + t15 * t25 + t3 * t371 + t4 * t51) + t217 * t283 + qJD(2) * t180 + Ifges(2,3) * qJDD(1) + t80 * (-mrSges(5,1) * t175 + mrSges(5,2) * t174) + t27 * (-mrSges(6,1) * t175 - mrSges(6,3) * t174) + t109 * Ifges(3,1) * t155; t350 * t154 + t351 * t216 + t352 * mrSges(6,2) + t353 * mrSges(5,3) + (t203 / 0.2e1 - t196 / 0.2e1) * t284 + t370 * t268 / 0.2e1 + t372 * qJD(3) + (t50 * qJ(3) - t25 * t42 - t26 * t43 + t349 * t80) * m(5) + (t5 * t111 - t16 * t29 - t17 * t28 + t27 * t373) * m(6) + t50 * t214 + t5 * t212 + (-m(5) * t163 - m(6) * t164 - t392) * t326 - t367 * t136 - t368 * t134 + t79 * t267 / 0.2e1 + (-m(4) * t271 + t155 * t394 + t380 * t158 - t246 * t401 + t390) * g(3) + t111 * t12 - t103 * t102 - t104 * t55 + t96 * mrSges(3,2) - t97 * mrSges(3,1) - t67 * mrSges(4,3) + t70 * mrSges(4,2) - pkin(2) * t74 - t43 * t59 - t28 * t60 - t42 * t61 - t29 * t62 - t110 * t247 - t114 * t248 + t359 * t332 + t361 * t226 + t362 * t329 - (t360 * t100 + t363 * t128) * qJD(4) / 0.2e1 - (-Ifges(3,2) * t268 + t133 + t347) * t267 / 0.2e1 + (-qJ(3) * t258 * t272 + t159 * t345) * g(1) + t385 * t157 / 0.2e1 + t128 * (t27 * t213 + t80 * t215) - t388 * t268 / 0.2e1 + t377 * t108 + t379 * t109 + (-pkin(2) * t70 - qJ(3) * t67 - qJD(3) * t114 - t103 * t81) * m(4) + t373 * t54 + t197 * t330 + t204 * t331 + (t354 - t376 / 0.2e1) * t262 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + (-t258 * t286 + t345) * t391 + (-t73 + t13) * qJ(3) + (-t33 / 0.2e1 + t30 / 0.2e1 - t26 * mrSges(5,3) - t17 * mrSges(6,2)) * t261 + (t365 * t399 + t366 * t398 + (t168 / 0.2e1 - t167 / 0.2e1) * t99 - m(4) * t185 * pkin(6) - t180 - t26 * (-t158 * mrSges(5,2) + mrSges(5,3) * t278) - t17 * (mrSges(6,2) * t278 + t158 * mrSges(6,3)) - t25 * (t158 * mrSges(5,1) - mrSges(5,3) * t282) - t16 * (-t158 * mrSges(6,1) + mrSges(6,2) * t282) + (t364 / 0.2e1 - t179 / 0.2e1 + t181) * qJD(1)) * qJD(1); (-t54 - t372) * qJD(2) + t258 * t314 + ((t102 + t173) * qJD(1) - t351 * t258) * t155 + t74 + (-qJD(2) * t27 + t190 * t268 + t164) * m(6) + (-qJD(2) * t80 - t188 * t268 + t163) * m(5) + (qJD(2) * t114 + t268 * t81 + t70) * m(4) + t392; (mrSges(5,2) * t80 - mrSges(6,3) * t27 - t354) * t99 - t375 * t25 + t213 * t314 + (-g(2) * (pkin(4) * t84 - qJ(5) * t85) - g(1) * (-pkin(4) * t82 + qJ(5) * t83) + t193 * t314 - t27 * t53 - pkin(4) * t2 + qJ(5) * t1 + t358 * t17) * m(6) + (t381 * t83 + t382 * t82) * g(1) + g(3) * t91 + qJD(5) * t60 - t53 * t54 + qJ(5) * t22 - pkin(4) * t20 - t383 + t384 + (-t381 * t85 - t382 * t84) * g(2) + t344 * mrSges(6,1) + (-m(6) * t16 + t304 - t374) * t26 + (-t100 * t400 - t378 * t99) * t398 + (-t397 * t99 - t295 + t30 + t93) * t399 + (-t100 * t80 + t4) * mrSges(5,1) + (-Ifges(5,2) * t100 + t376 - t94) * t327 + (Ifges(6,3) * t100 - t318) * t328 + t17 * t305 + t33 * t323 + t357; t100 * t54 - t128 * t60 + (-g(1) * t82 + g(2) * t84 - g(3) * t275 - t17 * t128 - t344) * m(6) + t20;];
tau = t6;
