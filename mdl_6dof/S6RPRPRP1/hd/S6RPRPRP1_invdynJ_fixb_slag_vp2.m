% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:21
% EndTime: 2019-03-09 03:01:46
% DurationCPUTime: 16.26s
% Computational Cost: add. (6654->599), mult. (14218->754), div. (0->0), fcn. (9587->14), ass. (0->273)
t406 = Ifges(6,4) + Ifges(7,4);
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t272 = qJD(1) * qJD(3);
t152 = qJDD(1) * t194 - t191 * t272;
t153 = qJDD(1) * t191 + t194 * t272;
t295 = sin(pkin(10));
t296 = cos(pkin(10));
t107 = t152 * t295 + t153 * t296;
t150 = t296 * t191 + t295 * t194;
t140 = t150 * qJD(1);
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t117 = qJD(3) * t193 - t140 * t190;
t55 = qJD(5) * t117 + qJDD(3) * t190 + t107 * t193;
t350 = t55 / 0.2e1;
t118 = qJD(3) * t190 + t140 * t193;
t56 = -qJD(5) * t118 + qJDD(3) * t193 - t107 * t190;
t349 = t56 / 0.2e1;
t106 = t152 * t296 - t153 * t295;
t105 = qJDD(5) - t106;
t348 = t105 / 0.2e1;
t390 = Ifges(6,1) + Ifges(7,1);
t388 = Ifges(6,5) + Ifges(7,5);
t387 = Ifges(6,2) + Ifges(7,2);
t386 = Ifges(7,6) + Ifges(6,6);
t385 = Ifges(7,3) + Ifges(6,3);
t325 = t193 * pkin(5);
t172 = pkin(4) + t325;
t223 = -mrSges(7,1) * t193 + mrSges(7,2) * t190;
t225 = -mrSges(6,1) * t193 + mrSges(6,2) * t190;
t405 = -m(6) * pkin(4) - m(7) * t172 + t223 + t225;
t273 = m(5) + m(6) + m(7);
t404 = -mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t403 = t388 * t348 + t406 * t349 + t390 * t350;
t402 = t406 * t117;
t186 = sin(pkin(9));
t168 = pkin(1) * t186 + pkin(7);
t157 = t168 * qJDD(1);
t401 = qJD(2) * qJD(3) + t157;
t163 = -mrSges(4,1) * t194 + mrSges(4,2) * t191;
t184 = qJ(3) + pkin(10);
t175 = sin(t184);
t177 = cos(t184);
t400 = -t177 * mrSges(5,1) - t175 * t404 + t163;
t159 = t168 * qJD(1);
t180 = t194 * qJDD(2);
t271 = qJD(1) * qJD(4);
t276 = qJD(3) * t194;
t62 = -t159 * t276 + qJDD(3) * pkin(3) - qJ(4) * t153 + t180 + (-t271 - t401) * t191;
t277 = qJD(3) * t191;
t95 = t191 * qJDD(2) - t159 * t277 + t194 * t401;
t66 = qJ(4) * t152 + t194 * t271 + t95;
t24 = t295 * t62 + t296 * t66;
t22 = qJDD(3) * pkin(8) + t24;
t274 = qJD(5) * t193;
t275 = qJD(5) * t190;
t187 = cos(pkin(9));
t170 = -pkin(1) * t187 - pkin(2);
t158 = t170 * qJDD(1);
t116 = -pkin(3) * t152 + qJDD(4) + t158;
t39 = -pkin(4) * t106 - pkin(8) * t107 + t116;
t181 = t194 * qJD(2);
t233 = qJ(4) * qJD(1) + t159;
t119 = -t191 * t233 + t181;
t112 = qJD(3) * pkin(3) + t119;
t278 = qJD(2) * t191;
t120 = t194 * t233 + t278;
t240 = t296 * t120;
t65 = t295 * t112 + t240;
t61 = qJD(3) * pkin(8) + t65;
t182 = t194 * pkin(3);
t156 = t170 - t182;
t138 = qJD(1) * t156 + qJD(4);
t149 = t191 * t295 - t194 * t296;
t139 = t149 * qJD(1);
t77 = pkin(4) * t139 - pkin(8) * t140 + t138;
t3 = t190 * t39 + t193 * t22 + t77 * t274 - t275 * t61;
t27 = t190 * t77 + t193 * t61;
t4 = -qJD(5) * t27 - t190 * t22 + t193 * t39;
t229 = -t190 * t4 + t193 * t3;
t26 = -t190 * t61 + t193 * t77;
t399 = -t26 * t274 - t27 * t275 + t229;
t398 = t406 * t118;
t397 = t406 * t193;
t396 = t406 * t190;
t185 = qJ(1) + pkin(9);
t176 = sin(t185);
t178 = cos(t185);
t395 = g(1) * t178 + g(2) * t176;
t394 = qJD(5) + t139;
t351 = m(7) * pkin(5);
t393 = m(3) + m(4);
t392 = -mrSges(6,1) - mrSges(7,1);
t391 = mrSges(6,2) + mrSges(7,2);
t384 = t105 * t386 + t387 * t56 + t406 * t55;
t382 = t117 * t386 + t118 * t388 + t385 * t394;
t381 = t117 * t387 + t386 * t394 + t398;
t380 = t118 * t390 + t388 * t394 + t402;
t379 = Ifges(5,4) * t139;
t378 = t139 * Ifges(5,2);
t297 = qJDD(3) / 0.2e1;
t255 = t295 * pkin(3);
t167 = t255 + pkin(8);
t281 = qJ(6) + t167;
t234 = qJD(5) * t281;
t293 = t139 * t190;
t109 = t295 * t120;
t68 = t119 * t296 - t109;
t280 = qJD(1) * t191;
t265 = pkin(3) * t280;
t89 = pkin(4) * t140 + pkin(8) * t139 + t265;
t30 = t190 * t89 + t193 * t68;
t377 = -qJ(6) * t293 + qJD(6) * t193 - t190 * t234 - t30;
t29 = -t190 * t68 + t193 * t89;
t292 = t139 * t193;
t376 = -pkin(5) * t140 - qJ(6) * t292 - qJD(6) * t190 - t193 * t234 - t29;
t100 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t107;
t19 = -mrSges(6,1) * t56 + mrSges(6,2) * t55;
t375 = t19 - t100;
t67 = t119 * t295 + t240;
t374 = -t67 + (t275 + t293) * pkin(5);
t373 = t351 + mrSges(7,1);
t303 = t140 * mrSges(5,3);
t372 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t117 + mrSges(6,2) * t118 + t303;
t371 = Ifges(5,5) * qJD(3);
t370 = Ifges(5,6) * qJD(3);
t317 = mrSges(7,2) * t193;
t222 = mrSges(7,1) * t190 + t317;
t224 = mrSges(6,1) * t190 + mrSges(6,2) * t193;
t64 = t112 * t296 - t109;
t60 = -qJD(3) * pkin(4) - t64;
t40 = -t117 * pkin(5) + qJD(6) + t60;
t369 = t40 * t222 + t60 * t224;
t368 = -t190 * t386 + t193 * t388;
t367 = -t190 * t387 + t397;
t366 = t193 * t390 - t396;
t364 = t105 * t385 + t386 * t56 + t388 * t55;
t132 = t159 * t194 + t278;
t96 = -qJD(3) * t132 - t157 * t191 + t180;
t363 = -t191 * t96 + t194 * t95;
t361 = mrSges(6,1) + t373;
t360 = 0.2e1 * t297;
t1 = pkin(5) * t105 - qJ(6) * t55 - qJD(6) * t118 + t4;
t15 = -qJ(6) * t118 + t26;
t14 = pkin(5) * t394 + t15;
t16 = qJ(6) * t117 + t27;
t2 = qJ(6) * t56 + qJD(6) * t117 + t3;
t357 = -t1 * t190 - t14 * t274 - t16 * t275 + t193 * t2;
t356 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t355 = m(4) * pkin(2) + mrSges(3,1) - t400;
t354 = t138 * mrSges(5,2) + t371 / 0.2e1;
t353 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t352 = t138 * mrSges(5,1) + t14 * mrSges(7,1) + t26 * mrSges(6,1) - t370 / 0.2e1 - t16 * mrSges(7,2) - t27 * mrSges(6,2);
t347 = -t117 / 0.2e1;
t346 = t117 / 0.2e1;
t345 = -t118 / 0.2e1;
t344 = t118 / 0.2e1;
t343 = -t394 / 0.2e1;
t342 = t394 / 0.2e1;
t341 = t139 / 0.2e1;
t340 = t140 / 0.2e1;
t339 = -t140 / 0.2e1;
t192 = sin(qJ(1));
t334 = pkin(1) * t192;
t331 = pkin(8) * t175;
t195 = cos(qJ(1));
t183 = t195 * pkin(1);
t31 = mrSges(7,1) * t105 - mrSges(7,3) * t55;
t32 = mrSges(6,1) * t105 - mrSges(6,3) * t55;
t321 = -t31 - t32;
t33 = -mrSges(7,2) * t105 + mrSges(7,3) * t56;
t34 = -mrSges(6,2) * t105 + mrSges(6,3) * t56;
t320 = t33 + t34;
t313 = mrSges(7,3) * t117;
t78 = -mrSges(7,2) * t394 + t313;
t315 = mrSges(6,3) * t117;
t79 = -mrSges(6,2) * t394 + t315;
t319 = t78 + t79;
t312 = mrSges(7,3) * t118;
t80 = mrSges(7,1) * t394 - t312;
t314 = mrSges(6,3) * t118;
t81 = mrSges(6,1) * t394 - t314;
t318 = t80 + t81;
t282 = qJ(4) + t168;
t146 = t282 * t191;
t148 = t282 * t194;
t98 = -t146 * t295 + t148 * t296;
t91 = t193 * t98;
t92 = pkin(4) * t149 - pkin(8) * t150 + t156;
t42 = t190 * t92 + t91;
t316 = mrSges(5,3) * t139;
t311 = Ifges(4,4) * t191;
t310 = Ifges(4,4) * t194;
t302 = t140 * Ifges(5,4);
t142 = t149 * qJD(3);
t291 = t142 * t190;
t290 = t142 * t193;
t289 = t150 * t190;
t288 = t150 * t193;
t188 = -qJ(6) - pkin(8);
t287 = t175 * t188;
t286 = t176 * t190;
t285 = t176 * t193;
t284 = t178 * t190;
t283 = t178 * t193;
t279 = qJD(1) * t194;
t235 = qJD(3) * t282;
t123 = qJD(4) * t194 - t191 * t235;
t124 = -qJD(4) * t191 - t194 * t235;
t76 = t123 * t296 + t124 * t295;
t141 = t150 * qJD(3);
t264 = pkin(3) * t277;
t90 = pkin(4) * t141 + pkin(8) * t142 + t264;
t269 = t190 * t90 + t193 * t76 + t92 * t274;
t73 = -mrSges(7,1) * t117 + mrSges(7,2) * t118;
t268 = t73 + t372;
t262 = mrSges(4,3) * t280;
t261 = mrSges(4,3) * t279;
t256 = t296 * pkin(3);
t254 = t150 * t274;
t18 = -t56 * mrSges(7,1) + t55 * mrSges(7,2);
t243 = -t275 / 0.2e1;
t241 = -t190 * t76 + t193 * t90;
t41 = -t190 * t98 + t193 * t92;
t236 = -t106 * mrSges(5,1) + t107 * mrSges(5,2);
t169 = -t256 - pkin(4);
t230 = pkin(4) * t177 + t331;
t23 = -t295 * t66 + t296 * t62;
t227 = mrSges(4,1) * t191 + mrSges(4,2) * t194;
t219 = t194 * Ifges(4,2) + t311;
t216 = Ifges(4,5) * t194 - Ifges(4,6) * t191;
t75 = t123 * t295 - t296 * t124;
t97 = t296 * t146 + t148 * t295;
t210 = t172 * t177 - t287;
t209 = qJ(6) * t142 - qJD(6) * t150;
t127 = -t177 * t284 + t285;
t125 = t177 * t286 + t283;
t206 = t254 - t291;
t205 = t150 * t275 + t290;
t204 = t170 * qJD(1) * t227;
t203 = t191 * (Ifges(4,1) * t194 - t311);
t21 = -qJDD(3) * pkin(4) - t23;
t189 = -qJ(4) - pkin(7);
t174 = Ifges(4,4) * t279;
t173 = t182 + pkin(2);
t162 = -qJD(3) * mrSges(4,2) + t261;
t160 = qJD(3) * mrSges(4,1) - t262;
t155 = t169 - t325;
t147 = t281 * t193;
t145 = t281 * t190;
t144 = Ifges(4,1) * t280 + Ifges(4,5) * qJD(3) + t174;
t143 = Ifges(4,6) * qJD(3) + qJD(1) * t219;
t137 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t153;
t136 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t152;
t131 = -t159 * t191 + t181;
t129 = -qJD(3) * mrSges(5,2) - t316;
t128 = t177 * t283 + t286;
t126 = -t177 * t285 + t284;
t101 = mrSges(5,1) * t139 + mrSges(5,2) * t140;
t99 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t106;
t94 = t140 * Ifges(5,1) + t371 - t379;
t93 = t302 + t370 - t378;
t72 = pkin(5) * t289 + t97;
t38 = pkin(5) * t206 + t75;
t35 = -qJ(6) * t289 + t42;
t28 = pkin(5) * t149 - qJ(6) * t288 + t41;
t9 = -qJD(5) * t42 + t241;
t8 = -t275 * t98 + t269;
t7 = -t56 * pkin(5) + qJDD(6) + t21;
t6 = -qJ(6) * t254 + (-qJD(5) * t98 + t209) * t190 + t269;
t5 = pkin(5) * t141 + t209 * t193 + (-t91 + (qJ(6) * t150 - t92) * t190) * qJD(5) + t241;
t10 = [t153 * t310 / 0.2e1 + t191 * (Ifges(4,4) * t152 + Ifges(4,5) * qJDD(3)) / 0.2e1 + t153 * t191 * Ifges(4,1) + (-t1 * t288 + t14 * t205 - t16 * t206 - t2 * t289) * mrSges(7,3) + (t205 * t26 - t206 * t27 - t288 * t4 - t289 * t3) * mrSges(6,3) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t187 - 0.2e1 * mrSges(3,2) * t186 + m(3) * (t186 ^ 2 + t187 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) - t384 * t289 / 0.2e1 + (-t284 * t351 + mrSges(2,1) * t192 + mrSges(2,2) * t195 + t393 * t334 - t273 * (-t178 * t189 - t334) + t392 * t126 - t391 * t125 + t356 * t178 + (m(5) * t173 - m(6) * (-t173 - t230) - m(7) * (-t173 - t210) + t355) * t176) * g(1) + (t203 + t194 * (-Ifges(4,2) * t191 + t310)) * t272 / 0.2e1 + m(6) * (t26 * t9 + t27 * t8 + t3 * t42 + t4 * t41) + m(5) * (t116 * t156 + t138 * t264 + t24 * t98 + t65 * t76) + t288 * t403 + (-t205 * t390 - t206 * t406) * t344 + (-t205 * t406 - t206 * t387) * t346 + m(7) * (t1 * t28 + t14 * t5 + t16 * t6 + t2 * t35 + t38 * t40 + t7 * t72) + (-t205 * t388 - t206 * t386) * t342 + (t204 + t216 * qJD(3) / 0.2e1) * qJD(3) + (t64 * mrSges(5,3) - t94 / 0.2e1 + t379 / 0.2e1 - Ifges(5,1) * t340 - t354) * t142 + (t116 * mrSges(5,2) - t23 * mrSges(5,3) + Ifges(5,1) * t107 + Ifges(5,4) * t106 + Ifges(5,5) * t360 + t21 * t224 + t7 * t222 + t243 * t380 + t348 * t368 + t349 * t367 + t350 * t366) * t150 + t152 * t219 / 0.2e1 + (-m(5) * t64 + m(6) * t60 + t372) * t75 + (-m(5) * t23 + m(6) * t21 + t375) * t97 + (-t160 * t276 + m(4) * ((-t131 * t194 - t132 * t191) * qJD(3) + t363) - t162 * t277 - t191 * t137 + t194 * t136) * t168 + (-t131 * t276 - t132 * t277 + t363) * mrSges(4,3) + (-t286 * t351 - mrSges(2,1) * t195 + mrSges(2,2) * t192 - t273 * (t178 * t173 - t176 * t189 + t183) - t393 * t183 + t392 * t128 - t391 * t127 + t356 * t176 + (-m(6) * t230 - m(7) * t210 - t355) * t178) * g(2) + (t116 * mrSges(5,1) + t1 * mrSges(7,1) - t24 * mrSges(5,3) - Ifges(5,4) * t107 - Ifges(5,2) * t106 - t360 * Ifges(5,6) + t385 * t348 + t386 * t349 + t388 * t350 + t353 + t364 / 0.2e1) * t149 + (-t254 / 0.2e1 + t291 / 0.2e1) * t381 + t101 * t264 - t380 * t290 / 0.2e1 + t40 * (mrSges(7,1) * t206 - mrSges(7,2) * t205) + t60 * (mrSges(6,1) * t206 - mrSges(6,2) * t205) - t143 * t277 / 0.2e1 + t144 * t276 / 0.2e1 + (m(4) * t170 + t163) * t158 + t194 * (Ifges(4,4) * t153 + Ifges(4,2) * t152 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t156 * t236 + t170 * (-mrSges(4,1) * t152 + mrSges(4,2) * t153) + t76 * t129 + t98 * t99 + t5 * t80 + t9 * t81 + t72 * t18 + t38 * t73 + t6 * t78 + t8 * t79 + t42 * t34 + t41 * t32 + t28 * t31 + t35 * t33 + (-t93 / 0.2e1 + t378 / 0.2e1 - Ifges(5,4) * t340 + t386 * t346 + t388 * t344 + t385 * t342 + t352 - t65 * mrSges(5,3) + t382 / 0.2e1) * t141 + (Ifges(4,5) * t191 + Ifges(4,6) * t194) * t297; m(3) * qJDD(2) + t191 * t136 + t194 * t137 + (-t160 * t191 + t162 * t194) * qJD(3) + (t18 + t375) * t149 + t268 * t141 - (-t190 * t318 + t193 * t319 + t129) * t142 + (-t273 - t393) * g(3) + m(5) * (-t141 * t64 - t142 * t65 - t149 * t23) + m(7) * (t14 * t291 + t141 * t40 + t149 * t7 - t16 * t290) + m(6) * (t141 * t60 + t149 * t21 + t26 * t291 - t27 * t290) + m(4) * (t191 * t95 + t194 * t96 + (-t131 * t191 + t132 * t194) * qJD(3)) + (t99 + t320 * t193 + t321 * t190 + (-t190 * t319 - t193 * t318) * qJD(5) + m(5) * t24 + m(7) * t357 + m(6) * t399) * t150; (-t81 * t274 - t79 * t275 + m(6) * ((-t190 * t27 - t193 * t26) * qJD(5) + t229) - t190 * t32 + t193 * t34) * t167 + (-t293 / 0.2e1 + t243) * t381 + (-t302 + t382) * t339 + t384 * t193 / 0.2e1 + (t190 * t388 + t193 * t386) * t348 - (-Ifges(4,2) * t280 + t144 + t174) * t279 / 0.2e1 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t193 * t387 + t396) * t349 + (t190 * t390 + t397) * t350 + (-t26 * t292 - t27 * t293 + t399) * mrSges(6,3) + t190 * t403 + (-t138 * t265 + t64 * t67 - t65 * t68 + (t23 * t296 + t24 * t295) * pkin(3)) * m(5) + (-t379 + t94) * t341 + (t117 * t367 + t118 * t366 + t368 * t394) * qJD(5) / 0.2e1 + t395 * (t227 + t273 * pkin(3) * t191 + (-m(6) * pkin(8) + m(7) * t188 - t404) * t177 + (mrSges(5,1) - t405) * t175) + (-m(6) * (t182 + t331) - m(5) * t182 - m(7) * (t182 - t287) + t405 * t177 + t400) * g(3) + (t261 - t162) * t131 + (t169 * t21 - t26 * t29 - t27 * t30 - t60 * t67) * m(6) + (t262 + t160) * t132 - (Ifges(5,2) * t341 - t343 * t385 - t345 * t388 - t347 * t386 + t352) * t140 + (t274 / 0.2e1 + t292 / 0.2e1) * t380 + t7 * t223 + t21 * t225 - t372 * t67 + t374 * t73 + t376 * t80 + t377 * t78 + (-t1 * t145 + t14 * t376 + t147 * t2 + t155 * t7 + t16 * t377 + t374 * t40) * m(7) - (Ifges(5,1) * t339 + t343 * t368 + t345 * t366 + t347 * t367 - t354 - t369) * t139 + t369 * qJD(5) + (-t204 - t203 * qJD(1) / 0.2e1) * qJD(1) + t99 * t255 + t100 * t256 + (-t14 * t292 - t16 * t293 + t357) * mrSges(7,3) + t65 * t303 + t143 * t280 / 0.2e1 - t216 * t272 / 0.2e1 - t101 * t265 + t169 * t19 + Ifges(4,6) * t152 + Ifges(4,5) * t153 + t155 * t18 - t145 * t31 + t147 * t33 - t68 * t129 + Ifges(5,6) * t106 + Ifges(5,5) * t107 - t29 * t81 - t95 * mrSges(4,2) + t96 * mrSges(4,1) - t30 * t79 + t23 * mrSges(5,1) - t24 * mrSges(5,2) - t64 * t316 + t93 * t340; t139 * t129 - t268 * t140 + (t319 * t394 - t321) * t193 + (-t318 * t394 + t320) * t190 + t236 + (-g(1) * t176 + g(2) * t178) * t273 + (t1 * t193 - t140 * t40 + t190 * t2 + t394 * (-t14 * t190 + t16 * t193)) * m(7) + (-t140 * t60 + t190 * t3 + t193 * t4 + t394 * (-t190 * t26 + t193 * t27)) * m(6) + (t139 * t65 + t140 * t64 + t116) * m(5); t353 + (-t118 * t387 + t380 + t402) * t347 + (-t127 * t361 + t128 * t391) * g(1) + (t125 * t361 - t126 * t391) * g(2) + (t117 * t388 - t118 * t386) * t343 + (t117 * t390 - t398) * t345 + (t190 * t373 + t224 + t317) * g(3) * t175 + t381 * t344 + t373 * t1 + (t81 + t314) * t27 + (-m(7) * (-t14 + t15) + t80 + t312) * t16 + (-t79 + t315) * t26 - t60 * (mrSges(6,1) * t118 + mrSges(6,2) * t117) - t40 * (mrSges(7,1) * t118 + mrSges(7,2) * t117) - t15 * t78 + t14 * t313 + t364 + ((-m(7) * t40 - t73) * t118 + t31) * pkin(5); -t117 * t78 + t118 * t80 + (g(3) * t177 - t16 * t117 + t14 * t118 - t175 * t395 + t7) * m(7) + t18;];
tau  = t10;
