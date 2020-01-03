% Calculate vector of inverse dynamics joint torques for
% S5RRPRP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:56
% EndTime: 2019-12-31 20:00:25
% DurationCPUTime: 16.45s
% Computational Cost: add. (4725->557), mult. (10768->721), div. (0->0), fcn. (7350->10), ass. (0->248)
t362 = m(5) + m(6);
t366 = m(4) + t362;
t369 = -mrSges(5,3) - mrSges(6,2);
t357 = Ifges(5,1) + Ifges(6,1);
t368 = -Ifges(5,4) + Ifges(6,5);
t356 = Ifges(6,4) + Ifges(5,5);
t355 = Ifges(5,6) - Ifges(6,6);
t354 = -Ifges(5,3) - Ifges(6,2);
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t204 = t171 * mrSges(6,1) + t168 * mrSges(6,3);
t206 = mrSges(5,1) * t171 - mrSges(5,2) * t168;
t367 = -t204 - t206;
t169 = sin(qJ(2));
t172 = cos(qJ(2));
t269 = sin(pkin(8));
t270 = cos(pkin(8));
t132 = t169 * t269 - t172 * t270;
t119 = t132 * qJD(1);
t363 = qJD(4) + t119;
t133 = t270 * t169 + t269 * t172;
t120 = t133 * qJD(1);
t187 = t171 * qJD(2) - t120 * t168;
t241 = qJD(1) * qJD(2);
t222 = t169 * t241;
t240 = qJDD(1) * t172;
t137 = -t222 + t240;
t138 = qJDD(1) * t169 + t172 * t241;
t97 = t137 * t269 + t138 * t270;
t53 = qJD(4) * t187 + qJDD(2) * t168 + t171 * t97;
t319 = t53 / 0.2e1;
t104 = qJD(2) * t168 + t120 * t171;
t54 = qJD(4) * t104 - t171 * qJDD(2) + t168 * t97;
t317 = t54 / 0.2e1;
t96 = t137 * t270 - t138 * t269;
t93 = qJDD(4) - t96;
t316 = t93 / 0.2e1;
t170 = sin(qJ(1));
t299 = g(2) * t170;
t143 = -mrSges(3,1) * t172 + mrSges(3,2) * t169;
t365 = -m(3) * pkin(1) - mrSges(2,1) + t143;
t361 = t137 / 0.2e1;
t360 = qJD(2) / 0.2e1;
t359 = mrSges(5,1) + mrSges(6,1);
t358 = mrSges(5,2) - mrSges(6,3);
t353 = t356 * t93 + t357 * t53 + t368 * t54;
t19 = mrSges(5,1) * t93 - mrSges(5,3) * t53;
t20 = -t93 * mrSges(6,1) + t53 * mrSges(6,2);
t352 = t20 - t19;
t21 = -mrSges(5,2) * t93 - mrSges(5,3) * t54;
t22 = -mrSges(6,2) * t54 + mrSges(6,3) * t93;
t351 = t22 + t21;
t350 = t104 * t356 + t187 * t355 - t354 * t363;
t102 = Ifges(5,4) * t187;
t279 = Ifges(6,5) * t187;
t349 = t104 * t357 + t356 * t363 + t102 - t279;
t288 = mrSges(6,2) * t187;
t60 = mrSges(6,3) * t363 + t288;
t286 = mrSges(5,3) * t187;
t61 = -mrSges(5,2) * t363 + t286;
t348 = t60 + t61;
t285 = mrSges(5,3) * t104;
t62 = mrSges(5,1) * t363 - t285;
t287 = mrSges(6,2) * t104;
t63 = -mrSges(6,1) * t363 + t287;
t347 = t62 - t63;
t271 = qJDD(2) / 0.2e1;
t193 = pkin(4) * t168 - qJ(5) * t171;
t167 = -qJ(3) - pkin(6);
t142 = t167 * t169;
t135 = qJD(1) * t142;
t144 = t167 * t172;
t136 = qJD(1) * t144;
t215 = t270 * t136;
t94 = t135 * t269 - t215;
t346 = -qJD(5) * t168 + t193 * t363 - t94;
t275 = t120 * mrSges(4,3);
t345 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t187 + mrSges(5,2) * t104 + t275;
t344 = Ifges(4,5) * qJD(2);
t343 = Ifges(4,6) * qJD(2);
t166 = qJ(2) + pkin(8);
t162 = sin(t166);
t163 = cos(t166);
t342 = t163 * pkin(3) + t162 * pkin(7);
t203 = t168 * mrSges(6,1) - t171 * mrSges(6,3);
t205 = mrSges(5,1) * t168 + mrSges(5,2) * t171;
t123 = t269 * t136;
t128 = qJD(2) * pkin(2) + t135;
t87 = t128 * t270 + t123;
t77 = -qJD(2) * pkin(3) - t87;
t27 = -pkin(4) * t187 - t104 * qJ(5) + t77;
t340 = t27 * t203 + t77 * t205;
t339 = -t168 * t355 + t171 * t356;
t278 = Ifges(6,5) * t168;
t281 = Ifges(5,4) * t168;
t338 = t171 * t357 + t278 - t281;
t122 = t132 * qJD(2);
t243 = qJD(4) * t171;
t180 = -t168 * t122 + t133 * t243;
t337 = -t354 * t93 - t355 * t54 + t356 * t53;
t266 = t119 * t171;
t336 = t243 + t266;
t244 = qJD(4) * t168;
t267 = t119 * t168;
t335 = -t244 - t267;
t160 = pkin(6) * t240;
t130 = -pkin(6) * t222 + t160;
t131 = t138 * pkin(6);
t334 = t130 * t172 + t131 * t169;
t268 = qJDD(1) * pkin(1);
t109 = -pkin(2) * t137 + qJDD(3) - t268;
t30 = -pkin(3) * t96 - pkin(7) * t97 + t109;
t247 = qJD(3) * t169;
t83 = qJDD(2) * pkin(2) - qJ(3) * t138 - qJD(1) * t247 - t131;
t248 = qJD(2) * t169;
t234 = pkin(6) * t248;
t246 = qJD(3) * t172;
t91 = qJ(3) * t137 + t160 + (-t234 + t246) * qJD(1);
t38 = t269 * t83 + t270 * t91;
t36 = qJDD(2) * pkin(7) + t38;
t302 = pkin(2) * t172;
t159 = pkin(1) + t302;
t139 = -qJD(1) * t159 + qJD(3);
t64 = pkin(3) * t119 - pkin(7) * t120 + t139;
t88 = t269 * t128 - t215;
t78 = qJD(2) * pkin(7) + t88;
t3 = t168 * t30 + t171 * t36 + t64 * t243 - t244 * t78;
t26 = t168 * t64 + t171 * t78;
t4 = -qJD(4) * t26 - t168 * t36 + t171 * t30;
t333 = -t4 * t168 + t3 * t171;
t1 = qJ(5) * t93 + qJD(5) * t363 + t3;
t2 = -pkin(4) * t93 + qJDD(5) - t4;
t332 = t1 * t171 + t2 * t168;
t329 = Ifges(6,5) * t319 + Ifges(6,6) * t316 - t53 * Ifges(5,4) / 0.2e1 - t93 * Ifges(5,6) / 0.2e1 + (Ifges(6,3) + Ifges(5,2)) * t317;
t328 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t327 = 0.2e1 * t271;
t207 = t163 * mrSges(4,1) - t162 * mrSges(4,2);
t326 = -t162 * t369 + t207;
t100 = t142 * t269 - t144 * t270;
t86 = pkin(3) * t132 - pkin(7) * t133 - t159;
t289 = t171 * t100 + t168 * t86;
t216 = qJD(2) * t167;
t117 = t169 * t216 + t246;
t118 = t172 * t216 - t247;
t69 = t117 * t270 + t118 * t269;
t121 = t133 * qJD(2);
t235 = pkin(2) * t248;
t71 = pkin(3) * t121 + pkin(7) * t122 + t235;
t13 = -qJD(4) * t289 - t168 * t69 + t171 * t71;
t325 = m(6) * pkin(4) + t359;
t324 = m(6) * qJ(5) - t358;
t173 = cos(qJ(1));
t257 = t163 * t173;
t323 = (-g(1) * t257 - t163 * t299) * pkin(7);
t322 = t4 * mrSges(5,1) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t1 * mrSges(6,3);
t318 = -t54 / 0.2e1;
t315 = t187 / 0.2e1;
t314 = -t187 / 0.2e1;
t313 = -t104 / 0.2e1;
t312 = t104 / 0.2e1;
t311 = -t363 / 0.2e1;
t309 = t119 / 0.2e1;
t308 = t120 / 0.2e1;
t307 = -t120 / 0.2e1;
t304 = t168 / 0.2e1;
t301 = pkin(6) * t172;
t298 = g(3) * t162;
t250 = qJD(1) * t169;
t236 = pkin(2) * t250;
t70 = pkin(3) * t120 + pkin(7) * t119 + t236;
t95 = t135 * t270 + t123;
t34 = t168 * t70 + t171 * t95;
t284 = Ifges(3,4) * t169;
t283 = Ifges(3,4) * t172;
t282 = Ifges(5,4) * t104;
t280 = Ifges(5,4) * t171;
t277 = Ifges(6,5) * t171;
t276 = t119 * mrSges(4,3);
t274 = t120 * Ifges(4,4);
t265 = t122 * t171;
t258 = t162 * t173;
t254 = t170 * t168;
t253 = t170 * t171;
t252 = t171 * t173;
t251 = t173 * t168;
t249 = qJD(1) * t172;
t237 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t250) * t301;
t101 = Ifges(6,5) * t104;
t40 = Ifges(6,6) * t363 - Ifges(6,3) * t187 + t101;
t232 = t40 * t304;
t231 = t270 * pkin(2);
t230 = t269 * pkin(2);
t223 = -t96 * mrSges(4,1) + t97 * mrSges(4,2);
t218 = -t244 / 0.2e1;
t217 = t243 / 0.2e1;
t211 = t173 * t159 - t167 * t170;
t37 = -t269 * t91 + t270 * t83;
t208 = mrSges(3,1) * t169 + mrSges(3,2) * t172;
t200 = t172 * Ifges(3,2) + t284;
t199 = -Ifges(5,2) * t168 + t280;
t197 = Ifges(3,5) * t172 - Ifges(3,6) * t169;
t195 = Ifges(6,3) * t168 + t277;
t194 = t171 * pkin(4) + t168 * qJ(5);
t25 = -t168 * t78 + t171 * t64;
t33 = -t168 * t95 + t171 * t70;
t68 = t117 * t269 - t270 * t118;
t99 = -t270 * t142 - t144 * t269;
t46 = -t100 * t168 + t171 * t86;
t186 = -pkin(3) - t194;
t183 = pkin(1) * t208;
t179 = t133 * t244 + t265;
t12 = -t100 * t244 + t168 * t71 + t171 * t69 + t86 * t243;
t178 = t169 * (Ifges(3,1) * t172 - t284);
t35 = -qJDD(2) * pkin(3) - t37;
t161 = Ifges(3,4) * t249;
t156 = -t231 - pkin(3);
t141 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t249;
t129 = -t231 + t186;
t127 = Ifges(3,1) * t250 + Ifges(3,5) * qJD(2) + t161;
t126 = Ifges(3,6) * qJD(2) + qJD(1) * t200;
t115 = t163 * t252 + t254;
t114 = t163 * t251 - t253;
t113 = t163 * t253 - t251;
t112 = t163 * t254 + t252;
t110 = Ifges(4,4) * t119;
t105 = -qJD(2) * mrSges(4,2) - t276;
t84 = mrSges(4,1) * t119 + mrSges(4,2) * t120;
t80 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t97;
t79 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t96;
t75 = t120 * Ifges(4,1) - t110 + t344;
t74 = -t119 * Ifges(4,2) + t274 + t343;
t57 = -mrSges(6,1) * t187 - mrSges(6,3) * t104;
t56 = pkin(4) * t104 - qJ(5) * t187;
t55 = t133 * t193 + t99;
t43 = Ifges(5,2) * t187 + Ifges(5,6) * t363 + t282;
t31 = -pkin(4) * t132 - t46;
t29 = qJ(5) * t132 + t289;
t24 = -pkin(4) * t120 - t33;
t23 = qJ(5) * t120 + t34;
t18 = qJ(5) * t363 + t26;
t17 = -pkin(4) * t363 + qJD(5) - t25;
t16 = mrSges(5,1) * t54 + mrSges(5,2) * t53;
t15 = mrSges(6,1) * t54 - mrSges(6,3) * t53;
t14 = -t193 * t122 + (qJD(4) * t194 - qJD(5) * t171) * t133 + t68;
t7 = -pkin(4) * t121 - t13;
t6 = qJ(5) * t121 + qJD(5) * t132 + t12;
t5 = t54 * pkin(4) - t53 * qJ(5) - t104 * qJD(5) + t35;
t8 = [(t356 * t121 - t357 * t179 + t180 * t368) * t312 + (-m(4) * t211 + t369 * t258 - t362 * (pkin(3) * t257 + pkin(7) * t258 + t211) - t325 * t115 - t324 * t114 + (-t207 + t365) * t173 + t328 * t170) * g(2) + ((mrSges(6,2) * t2 - mrSges(5,3) * t4 + t353 / 0.2e1) * t171 + (-t1 * mrSges(6,2) - t3 * mrSges(5,3) + t329) * t168 + t109 * mrSges(4,2) - t37 * mrSges(4,3) + Ifges(4,1) * t97 + Ifges(4,4) * t96 + Ifges(4,5) * t327 + t195 * t317 + t199 * t318 + t5 * t203 + t35 * t205 + t40 * t217 + t316 * t339 + t319 * t338 + t349 * t218) * t133 + (t172 * (-Ifges(3,2) * t169 + t283) + t178) * t241 / 0.2e1 + m(5) * (t12 * t26 + t13 * t25 + t289 * t3 + t4 * t46) + t289 * t21 + (-t121 * t354 - t179 * t356 - t180 * t355) * t363 / 0.2e1 - t141 * t234 - t143 * t268 + (t137 * t301 + t334) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t334) - t180 * t43 / 0.2e1 + (-Ifges(6,5) * t179 + Ifges(6,6) * t121 + Ifges(6,3) * t180) * t314 + (-Ifges(5,4) * t179 - Ifges(5,2) * t180 + Ifges(5,6) * t121) * t315 + m(4) * (t100 * t38 - t109 * t159 + t139 * t235 + t69 * t88) - t122 * t232 + t84 * t235 + t138 * t283 / 0.2e1 + t26 * (-mrSges(5,2) * t121 - mrSges(5,3) * t180) + (-Ifges(4,5) * t122 - Ifges(4,6) * t121 + t172 * t127) * t360 + t200 * t361 + (-Ifges(4,1) * t122 - Ifges(4,4) * t121) * t308 + (-t121 * t88 + t122 * t87) * mrSges(4,3) - t119 * (-Ifges(4,4) * t122 - Ifges(4,2) * t121) / 0.2e1 + t139 * (mrSges(4,1) * t121 - mrSges(4,2) * t122) + t25 * (mrSges(5,1) * t121 + mrSges(5,3) * t179) + t17 * (-mrSges(6,1) * t121 - mrSges(6,2) * t179) + t27 * (mrSges(6,1) * t180 + mrSges(6,3) * t179) + t77 * (mrSges(5,1) * t180 - mrSges(5,2) * t179) + t18 * (-mrSges(6,2) * t180 + mrSges(6,3) * t121) + Ifges(3,6) * t172 * t271 - t126 * t248 / 0.2e1 - t183 * t241 + (t325 * t113 + t324 * t112 + (m(4) * t159 - t362 * (-t159 - t342) + t326 - t365) * t170 + (t167 * t366 + t328) * t173) * g(1) - t159 * t223 - qJDD(2) * mrSges(3,2) * t301 + (-m(4) * t87 + m(5) * t77 + t345) * t68 + t29 * t22 + t31 * t20 + t46 * t19 + t55 * t15 - t349 * t265 / 0.2e1 + t14 * t57 + t350 * t121 / 0.2e1 + t6 * t60 + t12 * t61 + t13 * t62 + t7 * t63 + (t197 * t360 - t237) * qJD(2) + (Ifges(3,1) * t138 - pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t138) + Ifges(3,4) * t361 + t327 * Ifges(3,5)) * t169 + t100 * t79 + t69 * t105 - t121 * t74 / 0.2e1 - t122 * t75 / 0.2e1 - pkin(1) * (-mrSges(3,1) * t137 + mrSges(3,2) * t138) + m(6) * (t1 * t29 + t14 * t27 + t17 * t7 + t18 * t6 + t2 * t31 + t5 * t55) + t172 * (Ifges(3,4) * t138 + Ifges(3,2) * t137 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (-m(4) * t37 + m(5) * t35 + t16 - t80) * t99 + Ifges(2,3) * qJDD(1) + (t337 / 0.2e1 + t109 * mrSges(4,1) - t38 * mrSges(4,3) - Ifges(4,4) * t97 - Ifges(4,2) * t96 - Ifges(4,6) * t327 + Ifges(5,6) * t318 + Ifges(6,6) * t317 - t316 * t354 + t319 * t356 + t322) * t132; (-m(4) * t302 + t143 - t362 * (t302 + t342) + (-m(6) * t194 + t367) * t163 - t326) * g(3) + t74 * t308 + (-t139 * t236 + t87 * t94 - t88 * t95 + (t269 * t38 + t270 * t37) * pkin(2)) * m(4) - (-Ifges(3,2) * t250 + t127 + t161) * t249 / 0.2e1 - t5 * t204 - t35 * t206 + (t104 * t338 + t339 * t363) * qJD(4) / 0.2e1 - t84 * t236 + t278 * t317 + t281 * t318 + (pkin(6) * t141 + t126 / 0.2e1) * t250 + (t17 * t336 + t18 * t335 + t332) * mrSges(6,2) + (-t25 * t336 + t26 * t335 + t333) * mrSges(5,3) - (t43 / 0.2e1 - t40 / 0.2e1) * t267 + (-t277 + t280) * t319 + t43 * t218 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t208 + t366 * pkin(2) * t169 + (mrSges(4,2) + t369) * t163 + (m(5) * pkin(3) - m(6) * t186 + mrSges(4,1) - t367) * t162) * (g(1) * t173 + t299) + (t237 + (t183 - t178 / 0.2e1) * qJD(1)) * qJD(1) + (t129 * t5 - t17 * t24 - t18 * t23 + t346 * t27 + t323) * m(6) + (t156 * t35 - t25 * t33 - t26 * t34 - t77 * t94 + t323) * m(5) - t87 * t276 + t79 * t230 + t80 * t231 + (-t274 + t350) * t307 + (t217 + t266 / 0.2e1) * t349 + (-t347 * t243 - t348 * t244 + t351 * t171 + t352 * t168 + ((-t18 * t168 + t17 * t171) * qJD(4) + t332) * m(6) + ((-t26 * t168 - t25 * t171) * qJD(4) + t333) * m(5)) * (t230 + pkin(7)) + (-(-t199 / 0.2e1 + t195 / 0.2e1) * t187 + t232 + t340) * qJD(4) + (t75 - t110) * t309 - t197 * t241 / 0.2e1 - (Ifges(4,1) * t307 + t199 * t314 + t195 * t315 - t344 / 0.2e1 - t139 * mrSges(4,2) + t338 * t313 + t339 * t311 - t340) * t119 - t345 * t94 + t346 * t57 + t88 * t275 + t37 * mrSges(4,1) - t38 * mrSges(4,2) + t353 * t304 - t23 * t60 - t34 * t61 - t33 * t62 - t24 * t63 + (Ifges(5,2) * t318 - Ifges(6,3) * t317 + t316 * t355 - t329) * t171 - (Ifges(4,2) * t309 + t18 * mrSges(6,3) - Ifges(5,6) * t314 - Ifges(6,6) * t315 + t25 * mrSges(5,1) - t17 * mrSges(6,1) - t26 * mrSges(5,2) - t343 / 0.2e1 + t139 * mrSges(4,1) - t356 * t313 + t354 * t311) * t120 + (t316 * t356 + t319 * t357) * t168 + Ifges(4,6) * t96 + Ifges(4,5) * t97 - t95 * t105 + t129 * t15 - t130 * mrSges(3,2) - t131 * mrSges(3,1) + Ifges(3,6) * t137 + Ifges(3,5) * t138 + t156 * t16; t119 * t105 - (t57 + t345) * t120 + (t348 * t363 - t352) * t171 + (-t347 * t363 + t351) * t168 + t223 + (-g(1) * t170 + g(2) * t173) * t366 + (t1 * t168 - t120 * t27 - t171 * t2 + t363 * (t168 * t17 + t171 * t18)) * m(6) + (-t120 * t77 + t168 * t3 + t171 * t4 + t363 * (-t168 * t25 + t171 * t26)) * m(5) + (t119 * t88 + t120 * t87 + t109) * m(4); t18 * t287 - t17 * t288 + t43 * t312 + (-m(6) * t18 + t286 - t348) * t25 + (t114 * t359 + t115 * t358) * g(1) + (t112 * t359 + t113 * t358) * g(2) + (-m(6) * t17 + t285 + t347) * t26 - pkin(4) * t20 + qJ(5) * t22 + (t187 * t357 + t101 - t282 + t40) * t313 + (-t104 * t355 + t187 * t356) * t311 - t27 * (mrSges(6,1) * t104 - mrSges(6,3) * t187) - t77 * (mrSges(5,1) * t104 + mrSges(5,2) * t187) + (Ifges(6,3) * t104 + t279) * t315 + (t205 + t203) * t298 + (-Ifges(5,2) * t104 + t102 + t349) * t314 + t337 + (-t27 * t56 - pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t18 + t193 * t298 - g(2) * (-pkin(4) * t112 + qJ(5) * t113) - g(1) * (-pkin(4) * t114 + qJ(5) * t115)) * m(6) - t56 * t57 + qJD(5) * t60 + t322; t104 * t57 - t363 * t60 + (-g(1) * t114 - g(2) * t112 + t27 * t104 - t168 * t298 - t18 * t363 + t2) * m(6) + t20;];
tau = t8;
