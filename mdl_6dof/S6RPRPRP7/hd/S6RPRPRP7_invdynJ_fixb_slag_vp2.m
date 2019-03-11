% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:18
% EndTime: 2019-03-09 03:21:45
% DurationCPUTime: 18.11s
% Computational Cost: add. (6524->612), mult. (12997->786), div. (0->0), fcn. (8499->10), ass. (0->271)
t407 = Ifges(6,4) + Ifges(7,4);
t349 = cos(qJ(3));
t419 = t349 / 0.2e1;
t408 = Ifges(6,1) + Ifges(7,1);
t406 = -Ifges(7,5) - Ifges(6,5);
t405 = Ifges(6,2) + Ifges(7,2);
t404 = Ifges(6,6) + Ifges(7,6);
t403 = -Ifges(7,3) - Ifges(6,3);
t195 = cos(qJ(1));
t344 = g(2) * t195;
t188 = sin(pkin(9));
t309 = cos(pkin(9));
t236 = t309 * t349;
t192 = sin(qJ(3));
t288 = qJD(1) * t192;
t135 = -qJD(1) * t236 + t188 * t288;
t191 = sin(qJ(5));
t194 = cos(qJ(5));
t110 = qJD(3) * t194 + t135 * t191;
t418 = t407 * t110;
t196 = -pkin(1) - pkin(7);
t160 = qJDD(1) * t196 + qJDD(2);
t161 = qJD(1) * t196 + qJD(2);
t287 = qJD(3) * t192;
t108 = t349 * t160 - t161 * t287;
t280 = qJD(1) * qJD(3);
t150 = qJDD(1) * t349 - t192 * t280;
t253 = t349 * qJD(4);
t67 = qJDD(3) * pkin(3) - t150 * qJ(4) - qJD(1) * t253 + t108;
t254 = qJD(3) * t349;
t109 = t192 * t160 + t161 * t254;
t151 = -qJD(1) * t254 - t192 * qJDD(1);
t283 = t192 * qJD(4);
t79 = qJ(4) * t151 - qJD(1) * t283 + t109;
t32 = t188 * t67 + t309 * t79;
t24 = qJDD(3) * pkin(8) + t32;
t285 = qJD(5) * t194;
t286 = qJD(5) * t191;
t102 = -t188 * t150 + t151 * t309;
t103 = t150 * t309 + t188 * t151;
t281 = qJD(1) * qJD(2);
t162 = qJDD(1) * qJ(2) + t281;
t113 = -pkin(3) * t151 + qJDD(4) + t162;
t38 = -pkin(4) * t102 - pkin(8) * t103 + t113;
t152 = t349 * t161;
t255 = qJD(1) * t349;
t126 = -qJ(4) * t255 + t152;
t122 = qJD(3) * pkin(3) + t126;
t125 = -qJ(4) * t288 + t161 * t192;
t245 = t309 * t125;
t71 = t188 * t122 + t245;
t62 = qJD(3) * pkin(8) + t71;
t202 = -t188 * t349 - t192 * t309;
t136 = t202 * qJD(1);
t153 = pkin(3) * t288 + qJD(1) * qJ(2) + qJD(4);
t76 = -pkin(4) * t136 + pkin(8) * t135 + t153;
t3 = t191 * t38 + t194 * t24 + t76 * t285 - t286 * t62;
t30 = t191 * t76 + t194 * t62;
t4 = -qJD(5) * t30 - t191 * t24 + t194 * t38;
t235 = -t191 * t4 + t194 * t3;
t29 = -t191 * t62 + t194 * t76;
t417 = -t29 * t285 - t30 * t286 + t235;
t111 = qJD(3) * t191 - t135 * t194;
t416 = t407 * t111;
t415 = t407 * t194;
t414 = t407 * t191;
t413 = qJD(5) - t136;
t186 = qJ(3) + pkin(9);
t179 = sin(t186);
t180 = cos(t186);
t214 = t192 * mrSges(4,1) + mrSges(4,2) * t349;
t412 = mrSges(5,1) * t179 + t214 + (mrSges(5,2) - mrSges(6,3) - mrSges(7,3)) * t180;
t367 = m(7) * pkin(5);
t364 = -m(3) - m(4);
t101 = qJDD(5) - t102;
t54 = qJD(5) * t110 + qJDD(3) * t191 + t103 * t194;
t55 = -qJD(5) * t111 + qJDD(3) * t194 - t103 * t191;
t411 = t101 * t404 + t405 * t55 + t407 * t54;
t193 = sin(qJ(1));
t345 = g(1) * t193;
t410 = -mrSges(7,1) - mrSges(6,1);
t409 = mrSges(6,2) + mrSges(7,2);
t402 = -t101 * t406 + t407 * t55 + t408 * t54;
t16 = -mrSges(6,1) * t55 + mrSges(6,2) * t54;
t94 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t103;
t401 = t16 - t94;
t400 = t404 * t110 - t111 * t406 - t403 * t413;
t399 = t405 * t110 + t404 * t413 + t416;
t398 = t111 * t408 - t406 * t413 + t418;
t348 = pkin(3) * t188;
t170 = pkin(8) + t348;
t290 = qJ(6) + t170;
t240 = qJD(5) * t290;
t307 = t136 * t191;
t119 = t188 * t125;
t85 = t126 * t309 - t119;
t239 = pkin(3) * t255;
t87 = -t135 * pkin(4) - t136 * pkin(8) + t239;
t34 = t191 * t87 + t194 * t85;
t394 = qJ(6) * t307 + qJD(6) * t194 - t191 * t240 - t34;
t306 = t136 * t194;
t33 = -t191 * t85 + t194 * t87;
t393 = pkin(5) * t135 + qJ(6) * t306 - qJD(6) * t191 - t194 * t240 - t33;
t84 = t126 * t188 + t245;
t392 = -t84 + (t286 - t307) * pkin(5);
t391 = t367 + mrSges(7,1);
t330 = mrSges(5,3) * t135;
t390 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t110 + mrSges(6,2) * t111 - t330;
t389 = Ifges(5,5) * qJD(3);
t388 = Ifges(5,6) * qJD(3);
t387 = t179 * t193;
t339 = t194 * pkin(5);
t175 = pkin(4) + t339;
t231 = -mrSges(7,1) * t194 + mrSges(7,2) * t191;
t233 = -mrSges(6,1) * t194 + mrSges(6,2) * t191;
t386 = m(6) * pkin(4) + m(7) * t175 - t231 - t233;
t331 = mrSges(7,2) * t194;
t230 = mrSges(7,1) * t191 + t331;
t232 = mrSges(6,1) * t191 + mrSges(6,2) * t194;
t70 = t122 * t309 - t119;
t61 = -qJD(3) * pkin(4) - t70;
t40 = -t110 * pkin(5) + qJD(6) + t61;
t385 = t40 * t230 + t61 * t232;
t384 = -t191 * t404 - t194 * t406;
t383 = -t191 * t405 + t415;
t382 = t194 * t408 - t414;
t380 = -t101 * t403 + t404 * t55 - t406 * t54;
t379 = -t108 * t349 - t109 * t192;
t378 = t344 - t345;
t282 = -m(5) - m(6) - m(7);
t217 = mrSges(4,1) * t349 - mrSges(4,2) * t192;
t272 = Ifges(4,4) * t349;
t376 = (-Ifges(4,1) * t192 - t272) * t419 + qJ(2) * t217;
t375 = -mrSges(6,1) - t391;
t137 = -qJD(3) * t236 + t188 * t287;
t138 = t202 * qJD(3);
t144 = t188 * t192 - t236;
t31 = -t188 * t79 + t309 * t67;
t374 = -t137 * t71 + t138 * t70 - t144 * t31;
t373 = -m(4) * t379 + t192 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t151);
t372 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t1 = pkin(5) * t101 - qJ(6) * t54 - qJD(6) * t111 + t4;
t18 = -qJ(6) * t111 + t29;
t14 = pkin(5) * t413 + t18;
t19 = qJ(6) * t110 + t30;
t2 = qJ(6) * t55 + qJD(6) * t110 + t3;
t371 = -t1 * t191 - t14 * t285 - t19 * t286 + t194 * t2;
t370 = -g(1) * t387 + t179 * t344;
t189 = -qJ(6) - pkin(8);
t299 = t180 * t189;
t346 = pkin(8) * t180;
t369 = -m(7) * (t175 * t179 + t299) - mrSges(3,3) + mrSges(2,2) - m(6) * (pkin(4) * t179 - t346) - t412;
t368 = qJD(1) ^ 2;
t366 = t54 / 0.2e1;
t365 = t55 / 0.2e1;
t363 = t101 / 0.2e1;
t362 = -t110 / 0.2e1;
t360 = -t111 / 0.2e1;
t359 = t111 / 0.2e1;
t358 = -t413 / 0.2e1;
t356 = -t135 / 0.2e1;
t354 = -t136 / 0.2e1;
t23 = -qJDD(3) * pkin(4) - t31;
t7 = -t55 * pkin(5) + qJDD(6) + t23;
t341 = t144 * t7;
t184 = t192 * pkin(3);
t25 = mrSges(7,1) * t101 - mrSges(7,3) * t54;
t26 = mrSges(6,1) * t101 - mrSges(6,3) * t54;
t335 = -t25 - t26;
t27 = -mrSges(7,2) * t101 + mrSges(7,3) * t55;
t28 = -mrSges(6,2) * t101 + mrSges(6,3) * t55;
t334 = t27 + t28;
t326 = mrSges(7,3) * t110;
t72 = -mrSges(7,2) * t413 + t326;
t328 = mrSges(6,3) * t110;
t73 = -mrSges(6,2) * t413 + t328;
t333 = t72 + t73;
t325 = mrSges(7,3) * t111;
t74 = mrSges(7,1) * t413 - t325;
t327 = mrSges(6,3) * t111;
t75 = mrSges(6,1) * t413 - t327;
t332 = t74 + t75;
t172 = qJ(2) + t184;
t96 = -pkin(4) * t202 + pkin(8) * t144 + t172;
t291 = qJ(4) - t196;
t154 = t291 * t192;
t264 = t349 * t196;
t155 = -qJ(4) * t349 + t264;
t105 = -t154 * t309 + t188 * t155;
t99 = t194 * t105;
t42 = t191 * t96 + t99;
t329 = mrSges(5,3) * t136;
t324 = Ifges(4,4) * t192;
t323 = Ifges(5,4) * t135;
t314 = t144 * t23;
t305 = t137 * t191;
t304 = t137 * t194;
t303 = t138 * t191;
t302 = t144 * t191;
t301 = t144 * t194;
t298 = t191 * t193;
t297 = t191 * t195;
t158 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t255;
t295 = t192 * t158;
t294 = t193 * t194;
t293 = t194 * t138;
t292 = t194 * t195;
t289 = t195 * pkin(1) + t193 * qJ(2);
t284 = qJDD(1) * mrSges(3,2);
t163 = pkin(3) * t254 + qJD(2);
t123 = t287 * t291 - t253;
t124 = qJD(3) * t155 - t283;
t81 = t188 * t123 + t124 * t309;
t86 = -pkin(4) * t137 - pkin(8) * t138 + t163;
t276 = t191 * t86 + t194 * t81 + t96 * t285;
t57 = -mrSges(7,1) * t110 + mrSges(7,2) * t111;
t275 = t57 + t390;
t274 = t349 * pkin(3);
t263 = t309 * pkin(3);
t262 = t144 * t285;
t15 = -t55 * mrSges(7,1) + t54 * mrSges(7,2);
t248 = t285 / 0.2e1;
t183 = t195 * qJ(2);
t247 = -pkin(1) * t193 + t183;
t246 = -t191 * t81 + t194 * t86;
t244 = -t280 / 0.2e1;
t243 = -t102 * mrSges(5,1) + t103 * mrSges(5,2);
t41 = -t105 * t191 + t194 * t96;
t241 = (t162 + t281) * qJ(2);
t80 = -t309 * t123 + t124 * t188;
t104 = -t154 * t188 - t309 * t155;
t171 = -t263 - pkin(4);
t222 = t191 * t30 + t194 * t29;
t220 = -qJ(6) * t138 + qJD(6) * t144;
t216 = t349 * Ifges(4,1) - t324;
t215 = -Ifges(4,2) * t192 + t272;
t213 = -Ifges(4,5) * t192 - Ifges(4,6) * t349;
t131 = t179 * t297 + t294;
t129 = -t179 * t298 + t292;
t209 = t262 - t303;
t208 = t144 * t286 + t293;
t204 = t192 * (-Ifges(4,2) * t349 - t324);
t203 = -t191 * t333 - t194 * t332;
t190 = -qJ(4) - pkin(7);
t178 = -qJDD(1) * pkin(1) + qJDD(2);
t167 = t193 * t274;
t157 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t288;
t156 = t171 - t339;
t147 = t214 * qJD(1);
t142 = t290 * t194;
t141 = t290 * t191;
t140 = Ifges(4,5) * qJD(3) + qJD(1) * t216;
t139 = Ifges(4,6) * qJD(3) + qJD(1) * t215;
t133 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t150;
t132 = t179 * t292 - t298;
t130 = t179 * t294 + t297;
t127 = Ifges(5,4) * t136;
t117 = -qJD(3) * mrSges(5,2) + t329;
t95 = -mrSges(5,1) * t136 - mrSges(5,2) * t135;
t93 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t102;
t90 = -t135 * Ifges(5,1) + t127 + t389;
t89 = t136 * Ifges(5,2) - t323 + t388;
t68 = -pkin(5) * t302 + t104;
t39 = -pkin(5) * t209 + t80;
t37 = qJ(6) * t302 + t42;
t22 = -pkin(5) * t202 + qJ(6) * t301 + t41;
t13 = -qJD(5) * t42 + t246;
t12 = -t105 * t286 + t276;
t6 = qJ(6) * t262 + (-qJD(5) * t105 + t220) * t191 + t276;
t5 = -pkin(5) * t137 + t220 * t194 + (-t99 + (-qJ(6) * t144 - t96) * t191) * qJD(5) + t246;
t8 = [(qJD(5) * t398 + t411) * t302 / 0.2e1 + (-t144 * t382 + t202 * t406) * t366 + (-t297 * t367 + t364 * t289 + t282 * (t193 * t184 - t190 * t195 + t289) + t410 * t130 - t409 * t129 + (-m(4) * pkin(7) + t372) * t195 + t369 * t193) * g(2) + (t298 * t367 - m(3) * t247 - m(4) * t183 + t282 * (t195 * t184 + t193 * t190 + t247) + t410 * t132 + t409 * t131 + (-m(4) * t196 - t372) * t193 + t369 * t195) * g(1) + (-t137 * t404 + t208 * t407 + t209 * t405) * t110 / 0.2e1 + (t137 * t406 + t208 * t408 + t209 * t407) * t359 - t400 * t137 / 0.2e1 + (-m(5) * t31 + m(6) * t23 + t401) * t104 - t402 * t301 / 0.2e1 + t398 * t293 / 0.2e1 + m(4) * t241 + (-t303 / 0.2e1 + t144 * t248) * t399 + (Ifges(5,1) * t138 + Ifges(5,4) * t137) * t356 + t153 * (-mrSges(5,1) * t137 + mrSges(5,2) * t138) + t136 * (Ifges(5,4) * t138 + Ifges(5,2) * t137) / 0.2e1 + (0.2e1 * Ifges(4,5) * t419 - Ifges(5,5) * t144 - Ifges(4,6) * t192 + Ifges(5,6) * t202) * qJDD(3) - t230 * t341 + (-m(5) * t70 + m(6) * t61 + t390) * t80 - t232 * t314 - t140 * t287 / 0.2e1 - pkin(1) * t284 + (Ifges(4,1) * t150 + Ifges(4,4) * t151) * t419 + t133 * t264 + m(6) * (t12 * t30 + t13 * t29 + t3 * t42 + t4 * t41) + t379 * mrSges(4,3) + t204 * t244 + m(5) * (t105 * t32 + t113 * t172 + t153 * t163 + t71 * t81) + t376 * t280 + qJD(3) * (Ifges(5,5) * t138 + Ifges(5,6) * t137) / 0.2e1 + (t137 * t403 - t208 * t406 + t209 * t404) * t413 / 0.2e1 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (t157 * t254 - t158 * t287 + t373) * t196 + (-t144 * t384 + t202 * t403) * t363 + (-t144 * t383 - t202 * t404) * t365 - t380 * t202 / 0.2e1 + (t202 * t32 - t374) * mrSges(5,3) + t1 * (-mrSges(7,1) * t202 + mrSges(7,3) * t301) + t4 * (-mrSges(6,1) * t202 + mrSges(6,3) * t301) + t2 * (mrSges(7,2) * t202 + mrSges(7,3) * t302) + t3 * (mrSges(6,2) * t202 + mrSges(6,3) * t302) + t113 * (-mrSges(5,1) * t202 - mrSges(5,2) * t144) + t102 * (-Ifges(5,4) * t144 + Ifges(5,2) * t202) + t103 * (-Ifges(5,1) * t144 + Ifges(5,4) * t202) + (t214 + 0.2e1 * mrSges(3,3)) * t162 - t192 * (Ifges(4,4) * t150 + Ifges(4,2) * t151) / 0.2e1 - t139 * t254 / 0.2e1 + t172 * t243 + m(3) * (-pkin(1) * t178 + t241) + t178 * mrSges(3,2) + t163 * t95 + qJ(2) * (-mrSges(4,1) * t151 + mrSges(4,2) * t150) + qJD(2) * t147 + t138 * t90 / 0.2e1 + t137 * t89 / 0.2e1 + t81 * t117 + t105 * t93 + t13 * t75 + t68 * t15 + t6 * t72 + t12 * t73 + t5 * t74 + t39 * t57 + t41 * t26 + t42 * t28 + t37 * t27 + t22 * t25 + t14 * (-mrSges(7,1) * t137 - mrSges(7,3) * t208) + t29 * (-mrSges(6,1) * t137 - mrSges(6,3) * t208) + t30 * (mrSges(6,2) * t137 + mrSges(6,3) * t209) + t19 * (mrSges(7,2) * t137 + mrSges(7,3) * t209) + t40 * (-mrSges(7,1) * t209 + mrSges(7,2) * t208) + t61 * (-mrSges(6,1) * t209 + mrSges(6,2) * t208) + qJD(3) ^ 2 * t213 / 0.2e1 + t151 * t215 / 0.2e1 + t150 * t216 / 0.2e1 + m(7) * (t1 * t22 + t14 * t5 + t19 * t6 + t2 * t37 + t39 * t40 + t68 * t7); t284 + t349 * t133 + (t157 * t349 - t295) * qJD(3) + (qJ(2) * t364 - mrSges(3,3)) * t368 + (t15 + t401) * t144 - t275 * t138 + (t191 * t332 - t194 * t333 - t117) * t137 + m(7) * (-t138 * t40 + t14 * t305 - t19 * t304 + t341) + m(6) * (-t138 * t61 + t29 * t305 - t30 * t304 + t314) + m(5) * t374 + m(3) * t178 - (m(5) * t32 + m(6) * t417 + m(7) * t371 + t203 * qJD(5) + t335 * t191 + t334 * t194 + t93) * t202 + (-m(5) * t153 - t147 - t95 - m(7) * (t14 * t194 + t19 * t191) - m(6) * t222 + t203) * qJD(1) + t378 * (-t282 - t364) + t373; (-t153 * t239 + t70 * t84 - t71 * t85 + (t188 * t32 + t309 * t31) * pkin(3)) * m(5) + t411 * t194 / 0.2e1 + (-(-t175 * t180 + t179 * t189 - t274) * t344 - t1 * t141 + t142 * t2 + t156 * t7 - g(1) * (-t189 * t387 + t167) + t392 * t40 + t394 * t19 + t393 * t14) * m(7) + (-g(1) * (pkin(8) * t387 + t167) - (-pkin(4) * t180 - pkin(8) * t179 - t274) * t344 + t171 * t23 - t29 * t33 - t30 * t34 - t61 * t84) * m(6) + (-t191 * t406 + t194 * t404) * t363 + (Ifges(5,1) * t136 + t323 + t400) * t135 / 0.2e1 + t402 * t191 / 0.2e1 + (-t30 * mrSges(6,2) + t14 * mrSges(7,1) + t29 * mrSges(6,1) - t19 * mrSges(7,2) + Ifges(5,2) * t354 + t153 * mrSges(5,1) - t388 / 0.2e1 - t404 * t362 + t406 * t360 + t403 * t358) * t135 + (m(5) * t274 + mrSges(5,1) * t180 - mrSges(5,2) * t179 + t217) * t378 + (-t292 * t410 - t297 * t409) * t180 * g(2) + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (-m(6) * (-t184 + t346) + m(5) * t184 - m(7) * (-t184 - t299) + t386 * t179 + t412) * g(3) + (-t153 * mrSges(5,2) - t389 / 0.2e1 + t383 * t362 + t382 * t360 + t384 * t358 - t385) * t136 - t390 * t84 + t392 * t57 + t393 * t74 + t394 * t72 - t71 * t330 + t140 * t288 / 0.2e1 - t157 * t152 + t93 * t348 + t70 * t329 + t161 * t295 + t94 * t263 + t385 * qJD(5) + t213 * t244 + t89 * t356 + (t204 / 0.2e1 - t376) * t368 + (t110 * t383 + t111 * t382 + t384 * t413) * qJD(5) / 0.2e1 - t386 * t180 * t345 + (t194 * t28 - t191 * t26 + m(6) * (-qJD(5) * t222 + t235) - t75 * t285 - t73 * t286) * t170 + (t194 * t405 + t414) * t365 + (t191 * t408 + t415) * t366 + (t29 * t306 + t30 * t307 + t370 + t417) * mrSges(6,3) + t139 * t255 / 0.2e1 + (-t286 / 0.2e1 + t307 / 0.2e1) * t399 + (-t306 / 0.2e1 + t248) * t398 - t95 * t239 + (t14 * t306 + t19 * t307 + t370 + t371) * mrSges(7,3) + t171 * t16 + t156 * t15 + Ifges(4,5) * t150 + Ifges(4,6) * t151 - t141 * t25 + t142 * t27 - t85 * t117 + t108 * mrSges(4,1) - t109 * mrSges(4,2) + Ifges(5,6) * t102 + Ifges(5,5) * t103 - t33 * t75 - t34 * t73 + t31 * mrSges(5,1) - t32 * mrSges(5,2) + t7 * t231 + t23 * t233 + (t127 + t90) * t354; -t136 * t117 + t275 * t135 + (t333 * t413 - t335) * t194 + (-t332 * t413 + t334) * t191 + t243 + (g(1) * t195 + g(2) * t193) * t282 + (t1 * t194 + t135 * t40 + t2 * t191 + t413 * (-t14 * t191 + t19 * t194)) * m(7) + (t135 * t61 + t3 * t191 + t4 * t194 + t413 * (-t29 * t191 + t30 * t194)) * m(6) + (-t135 * t70 - t136 * t71 + t113) * m(5); (-t111 * t405 + t398 + t418) * t362 + (-t110 * t406 - t111 * t404) * t358 + (t110 * t408 - t416) * t360 + t391 * t1 + (t328 - t73) * t29 + t14 * t326 + (t191 * t391 + t232 + t331) * g(3) * t180 + t380 + (t131 * t375 - t132 * t409) * g(2) + (t129 * t375 + t130 * t409) * g(1) + t399 * t359 + (-m(7) * (-t14 + t18) + t325 + t74) * t19 - t61 * (mrSges(6,1) * t111 + mrSges(6,2) * t110) - t40 * (mrSges(7,1) * t111 + mrSges(7,2) * t110) - t18 * t72 - t2 * mrSges(7,2) - t3 * mrSges(6,2) + t4 * mrSges(6,1) + (t327 + t75) * t30 + ((-m(7) * t40 - t57) * t111 + t25) * pkin(5); -t110 * t72 + t111 * t74 + (-g(3) * t179 - t19 * t110 + t14 * t111 - t180 * t378 + t7) * m(7) + t15;];
tau  = t8;
