% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP8
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:26
% EndTime: 2019-03-09 03:24:52
% DurationCPUTime: 15.43s
% Computational Cost: add. (6587->638), mult. (12988->814), div. (0->0), fcn. (8445->10), ass. (0->283)
t399 = Ifges(6,1) + Ifges(7,1);
t411 = Ifges(6,4) - Ifges(7,5);
t398 = Ifges(7,4) + Ifges(6,5);
t397 = Ifges(6,6) - Ifges(7,6);
t396 = Ifges(6,3) + Ifges(7,2);
t187 = sin(pkin(9));
t190 = sin(qJ(3));
t312 = cos(pkin(9));
t348 = cos(qJ(3));
t201 = -t187 * t348 - t190 * t312;
t130 = t201 * qJD(1);
t405 = qJD(5) - t130;
t410 = t348 / 0.2e1;
t242 = t312 * t348;
t288 = qJD(1) * t190;
t129 = -qJD(1) * t242 + t187 * t288;
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t106 = qJD(3) * t189 - t129 * t192;
t279 = qJD(1) * qJD(3);
t143 = qJDD(1) * t348 - t190 * t279;
t255 = qJD(3) * t348;
t144 = -qJD(1) * t255 - t190 * qJDD(1);
t98 = t143 * t312 + t187 * t144;
t54 = qJD(5) * t106 - t192 * qJDD(3) + t189 * t98;
t364 = -t54 / 0.2e1;
t409 = Ifges(6,2) * t364;
t185 = qJ(3) + pkin(9);
t176 = sin(t185);
t177 = cos(t185);
t211 = t190 * mrSges(4,1) + mrSges(4,2) * t348;
t408 = -mrSges(5,1) * t176 - mrSges(5,2) * t177 - t211;
t194 = -pkin(1) - pkin(7);
t156 = qJDD(1) * t194 + qJDD(2);
t157 = qJD(1) * t194 + qJD(2);
t287 = qJD(3) * t190;
t103 = t348 * t156 - t157 * t287;
t254 = t348 * qJD(4);
t66 = qJDD(3) * pkin(3) - t143 * qJ(4) - qJD(1) * t254 + t103;
t104 = t190 * t156 + t157 * t255;
t282 = t190 * qJD(4);
t76 = qJ(4) * t144 - qJD(1) * t282 + t104;
t32 = t187 * t66 + t312 * t76;
t23 = qJDD(3) * pkin(8) + t32;
t284 = qJD(5) * t192;
t285 = qJD(5) * t189;
t280 = qJD(1) * qJD(2);
t158 = qJDD(1) * qJ(2) + t280;
t108 = -pkin(3) * t144 + qJDD(4) + t158;
t97 = -t187 * t143 + t144 * t312;
t36 = -pkin(4) * t97 - pkin(8) * t98 + t108;
t145 = t348 * t157;
t256 = qJD(1) * t348;
t119 = -qJ(4) * t256 + t145;
t115 = qJD(3) * pkin(3) + t119;
t118 = -qJ(4) * t288 + t157 * t190;
t248 = t312 * t118;
t69 = t187 * t115 + t248;
t62 = qJD(3) * pkin(8) + t69;
t146 = pkin(3) * t288 + qJD(1) * qJ(2) + qJD(4);
t74 = -pkin(4) * t130 + pkin(8) * t129 + t146;
t3 = t189 * t36 + t192 * t23 + t74 * t284 - t285 * t62;
t29 = t189 * t74 + t192 * t62;
t4 = -qJD(5) * t29 - t189 * t23 + t192 * t36;
t240 = -t189 * t4 + t192 * t3;
t28 = -t189 * t62 + t192 * t74;
t407 = -t28 * t284 - t29 * t285 + t240;
t17 = -pkin(5) * t405 + qJD(6) - t28;
t18 = qJ(6) * t405 + t29;
t96 = qJDD(5) - t97;
t1 = qJ(6) * t96 + qJD(6) * t405 + t3;
t2 = -pkin(5) * t96 + qJDD(6) - t4;
t241 = t1 * t192 + t189 * t2;
t406 = t17 * t284 - t18 * t285 + t241;
t26 = -mrSges(6,2) * t96 - mrSges(6,3) * t54;
t27 = -mrSges(7,2) * t54 + mrSges(7,3) * t96;
t336 = t26 + t27;
t219 = t192 * qJD(3) + t129 * t189;
t53 = qJD(5) * t219 + qJDD(3) * t189 + t192 * t98;
t24 = mrSges(6,1) * t96 - mrSges(6,3) * t53;
t25 = -t96 * mrSges(7,1) + t53 * mrSges(7,2);
t337 = -t24 + t25;
t404 = t337 * t189 + t336 * t192;
t365 = t53 / 0.2e1;
t362 = t96 / 0.2e1;
t361 = -m(3) - m(4);
t403 = -m(7) - m(6);
t193 = cos(qJ(1));
t345 = g(2) * t193;
t344 = g(3) * t177;
t402 = mrSges(6,1) + mrSges(7,1);
t401 = mrSges(6,2) - mrSges(7,3);
t400 = mrSges(6,3) + mrSges(7,2);
t395 = t398 * t96 + t399 * t53 - t411 * t54;
t16 = mrSges(6,1) * t54 + mrSges(6,2) * t53;
t89 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t98;
t394 = t16 - t89;
t393 = t106 * t398 + t219 * t397 + t396 * t405;
t102 = Ifges(6,4) * t219;
t322 = Ifges(7,5) * t219;
t392 = t106 * t399 + t398 * t405 + t102 - t322;
t332 = mrSges(7,2) * t219;
t70 = mrSges(7,3) * t405 + t332;
t328 = mrSges(6,3) * t219;
t71 = -mrSges(6,2) * t405 + t328;
t335 = t70 + t71;
t327 = mrSges(6,3) * t106;
t72 = mrSges(6,1) * t405 - t327;
t331 = mrSges(7,2) * t106;
t73 = -mrSges(7,1) * t405 + t331;
t334 = t72 - t73;
t227 = pkin(5) * t189 - qJ(6) * t192;
t80 = t119 * t187 + t248;
t391 = -qJD(6) * t189 + t227 * t405 - t80;
t330 = mrSges(5,3) * t129;
t390 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t219 + mrSges(6,2) * t106 - t330;
t389 = Ifges(5,5) * qJD(3);
t388 = Ifges(5,6) * qJD(3);
t191 = sin(qJ(1));
t301 = t176 * t191;
t235 = t189 * mrSges(7,1) - t192 * mrSges(7,3);
t237 = mrSges(6,1) * t189 + mrSges(6,2) * t192;
t112 = t187 * t118;
t68 = t115 * t312 - t112;
t61 = -qJD(3) * pkin(4) - t68;
t30 = -pkin(5) * t219 - t106 * qJ(6) + t61;
t387 = t30 * t235 + t61 * t237;
t214 = mrSges(4,1) * t348 - mrSges(4,2) * t190;
t274 = t348 * pkin(3);
t386 = -m(5) * t274 - mrSges(5,1) * t177 + mrSges(5,2) * t176 - t214;
t385 = -t189 * t397 + t192 * t398;
t321 = Ifges(7,5) * t189;
t324 = Ifges(6,4) * t189;
t384 = t192 * t399 + t321 - t324;
t381 = t396 * t96 - t397 * t54 + t398 * t53;
t380 = -t103 * t348 - t104 * t190;
t379 = Ifges(6,4) * t365 + t409 + Ifges(6,6) * t362 - t53 * Ifges(7,5) / 0.2e1 - t96 * Ifges(7,6) / 0.2e1 + Ifges(7,3) * t364;
t273 = Ifges(4,4) * t348;
t377 = (-Ifges(4,1) * t190 - t273) * t410 + qJ(2) * t214;
t313 = qJDD(3) / 0.2e1;
t376 = -t313 - qJDD(3) / 0.2e1;
t131 = -qJD(3) * t242 + t187 * t287;
t132 = t201 * qJD(3);
t137 = t187 * t190 - t242;
t31 = -t187 * t76 + t312 * t66;
t375 = t131 * t69 - t132 * t68 + t137 * t31;
t374 = -m(4) * t380 + t190 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t144);
t373 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t372 = -mrSges(3,3) + mrSges(2,2) + t408;
t291 = qJ(4) - t194;
t147 = t291 * t190;
t264 = t348 * t194;
t148 = -qJ(4) * t348 + t264;
t100 = -t147 * t312 + t187 * t148;
t183 = t190 * pkin(3);
t169 = qJ(2) + t183;
t91 = -pkin(4) * t201 + pkin(8) * t137 + t169;
t333 = t192 * t100 + t189 * t91;
t116 = t287 * t291 - t254;
t117 = qJD(3) * t148 - t282;
t78 = t187 * t116 + t117 * t312;
t159 = pkin(3) * t255 + qJD(2);
t82 = -pkin(4) * t131 - pkin(8) * t132 + t159;
t13 = -qJD(5) * t333 - t189 * t78 + t192 * t82;
t371 = -m(7) * pkin(5) - t402;
t370 = -m(7) * qJ(6) + t401;
t369 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t368 = qJD(1) ^ 2;
t363 = t54 / 0.2e1;
t360 = t219 / 0.2e1;
t359 = -t219 / 0.2e1;
t358 = -t106 / 0.2e1;
t357 = t106 / 0.2e1;
t356 = -t405 / 0.2e1;
t354 = -t129 / 0.2e1;
t353 = t129 / 0.2e1;
t352 = -t130 / 0.2e1;
t349 = t189 / 0.2e1;
t347 = pkin(3) * t187;
t346 = g(1) * t191;
t22 = -qJDD(3) * pkin(4) - t31;
t5 = t54 * pkin(5) - t53 * qJ(6) - t106 * qJD(6) + t22;
t342 = t137 * t5;
t81 = t119 * t312 - t112;
t244 = pkin(3) * t256;
t83 = -t129 * pkin(4) - t130 * pkin(8) + t244;
t34 = t189 * t83 + t192 * t81;
t329 = mrSges(5,3) * t130;
t326 = Ifges(4,4) * t190;
t325 = Ifges(5,4) * t129;
t323 = Ifges(6,4) * t192;
t320 = Ifges(7,5) * t192;
t319 = t106 * Ifges(6,4);
t316 = t137 * t22;
t314 = t176 * mrSges(7,2);
t310 = t130 * t189;
t309 = t130 * t192;
t308 = t131 * t189;
t307 = t131 * t192;
t306 = t132 * t192;
t304 = t137 * t192;
t300 = t176 * t193;
t299 = t177 * t191;
t298 = t177 * t193;
t297 = t189 * t191;
t296 = t189 * t193;
t150 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t256;
t294 = t190 * t150;
t293 = t191 * t192;
t292 = t193 * t192;
t289 = t193 * pkin(1) + t191 * qJ(2);
t283 = qJDD(1) * mrSges(3,2);
t281 = -m(5) + t403;
t57 = -mrSges(7,1) * t219 - mrSges(7,3) * t106;
t275 = t57 + t390;
t101 = Ifges(7,5) * t106;
t42 = Ifges(7,6) * t405 - Ifges(7,3) * t219 + t101;
t268 = -t189 * t42 / 0.2e1;
t45 = Ifges(6,2) * t219 + Ifges(6,6) * t405 + t319;
t267 = t45 * t349;
t266 = pkin(4) * t299 + pkin(8) * t301 + t191 * t274;
t263 = t312 * pkin(3);
t257 = -t97 * mrSges(5,1) + t98 * mrSges(5,2);
t251 = t285 / 0.2e1;
t250 = t284 / 0.2e1;
t182 = t193 * qJ(2);
t249 = -pkin(1) * t191 + t182;
t247 = -t279 / 0.2e1;
t245 = (t158 + t280) * qJ(2);
t77 = -t312 * t116 + t117 * t187;
t99 = -t147 * t187 - t312 * t148;
t238 = mrSges(6,1) * t192 - mrSges(6,2) * t189;
t236 = t192 * mrSges(7,1) + t189 * mrSges(7,3);
t232 = -Ifges(6,2) * t189 + t323;
t229 = Ifges(7,3) * t189 + t320;
t228 = t192 * pkin(5) + t189 * qJ(6);
t226 = t17 * t192 - t18 * t189;
t223 = t189 * t29 + t192 * t28;
t33 = -t189 * t81 + t192 * t83;
t40 = -t100 * t189 + t192 * t91;
t188 = -qJ(4) - pkin(7);
t218 = t193 * t183 + t191 * t188 + t249;
t217 = -pkin(4) - t228;
t216 = t191 * t183 - t188 * t193 + t289;
t215 = -pkin(8) * t176 - t274;
t213 = t348 * Ifges(4,1) - t326;
t212 = -Ifges(4,2) * t190 + t273;
t210 = -Ifges(4,5) * t190 - Ifges(4,6) * t348;
t206 = -t132 * t189 + t137 * t284;
t205 = t137 * t285 + t306;
t12 = -t100 * t285 + t189 * t82 + t192 * t78 + t91 * t284;
t203 = t190 * (-Ifges(4,2) * t348 - t326);
t202 = -t189 * t335 - t192 * t334;
t196 = m(7) * t217 - t236;
t175 = -qJDD(1) * pkin(1) + qJDD(2);
t168 = -t263 - pkin(4);
t149 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t288;
t140 = t211 * qJD(1);
t135 = -t263 + t217;
t134 = Ifges(4,5) * qJD(3) + qJD(1) * t213;
t133 = Ifges(4,6) * qJD(3) + qJD(1) * t212;
t126 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t143;
t125 = t176 * t292 - t297;
t124 = t176 * t296 + t293;
t123 = t176 * t293 + t296;
t122 = t176 * t297 - t292;
t120 = Ifges(5,4) * t130;
t110 = -qJD(3) * mrSges(5,2) + t329;
t90 = -mrSges(5,1) * t130 - mrSges(5,2) * t129;
t88 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t97;
t86 = -t129 * Ifges(5,1) + t120 + t389;
t85 = t130 * Ifges(5,2) - t325 + t388;
t56 = pkin(5) * t106 - qJ(6) * t219;
t55 = -t137 * t227 + t99;
t38 = pkin(5) * t201 - t40;
t37 = -qJ(6) * t201 + t333;
t20 = pkin(5) * t129 - t33;
t19 = -qJ(6) * t129 + t34;
t15 = mrSges(7,1) * t54 - mrSges(7,3) * t53;
t14 = t227 * t132 + (-qJD(5) * t228 + qJD(6) * t192) * t137 + t77;
t7 = pkin(5) * t131 - t13;
t6 = -qJ(6) * t131 - qJD(6) * t201 + t12;
t8 = [(t392 * t251 + (mrSges(7,2) * t1 + mrSges(6,3) * t3 + t379) * t189 - t108 * mrSges(5,2) - Ifges(5,1) * t98 - Ifges(5,4) * t97 + Ifges(5,5) * t376 - t229 * t363 - t232 * t364 + t45 * t250 - t362 * t385 - t365 * t384) * t137 + (-mrSges(7,2) * t2 + mrSges(6,3) * t4) * t304 + m(5) * (t100 * t32 + t108 * t169 + t146 * t159 + t69 * t78) + t392 * t306 / 0.2e1 - t393 * t131 / 0.2e1 + (-m(5) * t31 + m(6) * t22 + t394) * t99 - (qJD(5) * t42 + t395) * t304 / 0.2e1 + (-m(5) * t68 + m(6) * t61 + t390) * t77 + t380 * mrSges(4,3) + t377 * t279 + t375 * mrSges(5,3) + (t149 * t255 - t150 * t287 + t374) * t194 + (Ifges(4,5) * t348 - Ifges(4,6) * t190) * t313 + m(7) * (t1 * t37 + t14 * t30 + t17 * t7 + t18 * t6 + t2 * t38 + t5 * t55) - t235 * t342 - t237 * t316 - t134 * t287 / 0.2e1 - pkin(1) * t283 + t203 * t247 + m(4) * t245 + (Ifges(6,4) * t205 + Ifges(6,2) * t206 - Ifges(6,6) * t131) * t360 + (Ifges(7,5) * t205 - Ifges(7,6) * t131 - Ifges(7,3) * t206) * t359 + (Ifges(4,1) * t143 + Ifges(4,4) * t144 + Ifges(4,5) * qJDD(3)) * t410 + (-t398 * t131 + t399 * t205 + t206 * t411) * t357 - t132 * t267 - t132 * t268 + m(6) * (t12 * t29 + t13 * t28 + t3 * t333 + t4 * t40) + t333 * t26 + (-t131 * t396 + t205 * t398 + t206 * t397) * t405 / 0.2e1 + (Ifges(5,1) * t132 + Ifges(5,4) * t131) * t354 + t146 * (-mrSges(5,1) * t131 + mrSges(5,2) * t132) + t130 * (Ifges(5,4) * t132 + Ifges(5,2) * t131) / 0.2e1 + qJD(3) * (Ifges(5,5) * t132 + Ifges(5,6) * t131) / 0.2e1 + m(3) * (-pkin(1) * t175 + t245) + t126 * t264 - t190 * (Ifges(4,4) * t143 + Ifges(4,2) * t144 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t175 * mrSges(3,2) + t159 * t90 + qJD(2) * t140 + qJ(2) * (-mrSges(4,1) * t144 + mrSges(4,2) * t143) + t132 * t86 / 0.2e1 + t131 * t85 / 0.2e1 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + t78 * t110 + t100 * t88 + t13 * t72 + t7 * t73 + t6 * t70 + t12 * t71 + t55 * t15 + (-t108 * mrSges(5,1) + t32 * mrSges(5,3) + Ifges(5,4) * t98 + Ifges(5,2) * t97 - Ifges(5,6) * t376 - Ifges(6,6) * t364 - Ifges(7,6) * t363 - t362 * t396 - t365 * t398 - t369 - t381 / 0.2e1) * t201 + (t211 + 0.2e1 * mrSges(3,3)) * t158 + t14 * t57 + t40 * t24 + t37 * t27 + t38 * t25 + (-m(3) * t249 - m(4) * t182 - m(5) * t218 + t400 * t298 + t403 * (pkin(4) * t300 - pkin(8) * t298 + t218) + t371 * t125 + t370 * t124 + t372 * t193 + (-m(4) * t194 - t373) * t191) * g(1) + (-m(5) * t216 + t400 * t299 + t361 * t289 + t403 * (pkin(4) * t301 - pkin(8) * t299 + t216) + t371 * t123 + t370 * t122 + (-m(4) * pkin(7) + t373) * t193 + t372 * t191) * g(2) + t17 * (mrSges(7,1) * t131 + mrSges(7,2) * t205) + t28 * (-mrSges(6,1) * t131 - mrSges(6,3) * t205) + t18 * (mrSges(7,2) * t206 - mrSges(7,3) * t131) + t29 * (mrSges(6,2) * t131 + mrSges(6,3) * t206) + t61 * (-mrSges(6,1) * t206 + mrSges(6,2) * t205) + t30 * (-mrSges(7,1) * t206 - mrSges(7,3) * t205) + qJD(3) ^ 2 * t210 / 0.2e1 + t144 * t212 / 0.2e1 + t143 * t213 / 0.2e1 - t133 * t255 / 0.2e1 + t169 * t257; t283 + t348 * t126 + (t149 * t348 - t294) * qJD(3) + (qJ(2) * t361 - mrSges(3,3)) * t368 + (t15 + t394) * t137 - t275 * t132 + (t189 * t334 - t192 * t335 - t110) * t131 + m(6) * (-t132 * t61 + t28 * t308 - t29 * t307 + t316) + m(7) * (-t132 * t30 - t17 * t308 - t18 * t307 + t342) - m(5) * t375 + m(3) * t175 - (m(5) * t32 + m(6) * t407 + m(7) * t406 + t202 * qJD(5) + t404 + t88) * t201 + (-m(5) * t146 - m(6) * t223 + m(7) * t226 - t140 + t202 - t90) * qJD(1) + (t345 - t346) * (-t281 - t361) + t374; (-t320 + t323) * t365 + (-t17 * t309 + t18 * t310 - t344 + t406) * mrSges(7,2) + (-g(1) * t301 + g(2) * t300 + t28 * t309 + t29 * t310 + t407) * mrSges(6,3) + (m(5) * t183 - t177 * mrSges(6,3) + t403 * (t177 * pkin(8) - t183) + (m(6) * pkin(4) - t196 + t238) * t176 - t408) * g(3) + (t28 * mrSges(6,1) - t17 * mrSges(7,1) - t29 * mrSges(6,2) + t18 * mrSges(7,3) - Ifges(7,6) * t360 + Ifges(5,2) * t352 - Ifges(6,6) * t359 + t146 * mrSges(5,1) - t388 / 0.2e1 - t398 * t358 - t396 * t356) * t129 + (t362 * t398 + t365 * t399) * t189 + (t325 + t393) * t353 + t395 * t349 + (t229 * t360 + Ifges(5,1) * t353 + t232 * t359 + t267 + t268 - t146 * mrSges(5,2) - t389 / 0.2e1 + t384 * t358 + t385 * t356 - t387) * t130 - t390 * t80 + t391 * t57 + (-m(7) * t215 - t177 * t196 + t314 - t386) * t345 + t386 * t346 + (t203 / 0.2e1 - t377) * t368 + (-m(7) * t266 + (-t314 + (-m(7) * t228 - t236 - t238) * t177) * t191) * g(1) + (-t334 * t284 - t335 * t285 + (qJD(5) * t226 + t241) * m(7) + (-qJD(5) * t223 + t240) * m(6) + t404) * (pkin(8) + t347) - t69 * t330 + (t86 + t120) * t352 + t134 * t288 / 0.2e1 - t45 * t285 / 0.2e1 - t149 * t145 + t42 * t251 + t210 * t247 + t85 * t354 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (-Ifges(7,3) * t363 + t362 * t397 + t379 + t409) * t192 + (t135 * t5 - t17 * t20 - t18 * t19 + t391 * t30) * m(7) + (-(-pkin(4) * t177 + t215) * t345 - g(1) * t266 - t28 * t33 - t29 * t34 - t61 * t80 + t168 * t22) * m(6) + (t106 * t384 + t385 * t405) * qJD(5) / 0.2e1 + t321 * t363 + t324 * t364 + (-t146 * t244 + t68 * t80 - t69 * t81 + (t187 * t32 + t31 * t312) * pkin(3)) * m(5) - t90 * t244 + t88 * t347 + t68 * t329 + t157 * t294 + t168 * t16 + Ifges(4,5) * t143 + Ifges(4,6) * t144 + t135 * t15 - t81 * t110 + t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(5,5) * t98 + Ifges(5,6) * t97 - t5 * t236 - t22 * t238 - t33 * t72 - t20 * t73 - t19 * t70 - t34 * t71 + (t387 - (-t232 / 0.2e1 + t229 / 0.2e1) * t219) * qJD(5) + t31 * mrSges(5,1) - t32 * mrSges(5,2) + (t250 - t309 / 0.2e1) * t392 - g(2) * (-mrSges(6,1) * t292 + mrSges(6,2) * t296) * t177 + t89 * t263 + t133 * t256 / 0.2e1; -t130 * t110 + t275 * t129 + (t335 * t405 - t337) * t192 + (-t334 * t405 + t336) * t189 + t257 + (g(1) * t193 + g(2) * t191) * t281 + (t1 * t189 + t129 * t30 - t2 * t192 + t405 * (t17 * t189 + t18 * t192)) * m(7) + (t129 * t61 + t3 * t189 + t4 * t192 + t405 * (-t28 * t189 + t29 * t192)) * m(6) + (-t129 * t68 - t130 * t69 + t108) * m(5); t381 + (-Ifges(6,2) * t106 + t102 + t392) * t359 + (-t402 * t124 - t401 * t125) * g(2) + t45 * t357 + (-t106 * t397 + t219 * t398) * t356 + (t219 * t399 + t101 - t319 + t42) * t358 - t30 * (mrSges(7,1) * t106 - mrSges(7,3) * t219) - t61 * (mrSges(6,1) * t106 + mrSges(6,2) * t219) + (-m(7) * t17 + t327 + t334) * t29 + (t235 + t237) * t344 + (Ifges(7,3) * t106 + t322) * t360 + (-m(7) * t18 + t328 - t335) * t28 + (t122 * t402 + t123 * t401) * g(1) + t369 + t18 * t331 - t17 * t332 + (t227 * t344 - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t18 - t30 * t56 - g(2) * (pkin(5) * t124 - qJ(6) * t125) - g(1) * (-pkin(5) * t122 + qJ(6) * t123)) * m(7) + qJD(6) * t70 - t56 * t57 - pkin(5) * t25 + qJ(6) * t27; t106 * t57 - t405 * t70 + (-g(1) * t122 + g(2) * t124 + t30 * t106 - t18 * t405 - t189 * t344 + t2) * m(7) + t25;];
tau  = t8;
