% Calculate vector of inverse dynamics joint torques for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:31
% EndTime: 2019-12-05 18:05:57
% DurationCPUTime: 10.01s
% Computational Cost: add. (3827->428), mult. (9236->564), div. (0->0), fcn. (6284->10), ass. (0->228)
t334 = Ifges(5,4) + Ifges(6,4);
t335 = Ifges(5,1) + Ifges(6,1);
t333 = Ifges(5,5) + Ifges(6,5);
t332 = Ifges(5,2) + Ifges(6,2);
t331 = Ifges(5,6) + Ifges(6,6);
t194 = sin(pkin(8));
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t251 = qJD(1) * qJD(3);
t136 = (qJDD(1) * t200 - t197 * t251) * t194;
t137 = (-qJDD(1) * t197 - t200 * t251) * t194;
t196 = sin(qJ(4));
t199 = cos(qJ(4));
t218 = t196 * t197 - t199 * t200;
t214 = t218 * qJD(4);
t258 = qJD(1) * t194;
t52 = -t136 * t196 + t137 * t199 + t214 * t258;
t358 = t331 * t52;
t152 = t196 * t200 + t197 * t199;
t215 = qJD(1) * t152;
t109 = t194 * t215;
t51 = -qJD(4) * t109 + t136 * t199 + t137 * t196;
t357 = t333 * t51;
t195 = cos(pkin(8));
t248 = qJDD(1) * t195;
t173 = qJDD(3) - t248;
t163 = qJDD(4) + t173;
t356 = (Ifges(5,3) + Ifges(6,3)) * t163;
t250 = qJDD(1) * qJ(2);
t252 = qJD(1) * qJD(2);
t170 = t250 + t252;
t355 = t170 + t252;
t337 = mrSges(6,1) + mrSges(5,1);
t336 = -mrSges(6,2) - mrSges(5,2);
t353 = t334 * t109;
t235 = t200 * t258;
t236 = t197 * t258;
t111 = -t196 * t236 + t199 * t235;
t352 = t334 * t111;
t257 = qJD(1) * t195;
t175 = qJD(3) - t257;
t164 = qJD(4) + t175;
t351 = -t332 * t109 + t331 * t164 + t352;
t350 = t335 * t111 + t333 * t164 - t353;
t349 = Ifges(4,5) * t136;
t348 = Ifges(4,6) * t137;
t347 = Ifges(4,3) * t173;
t191 = t194 ^ 2;
t346 = t355 * t191;
t193 = qJ(3) + qJ(4);
t186 = sin(t193);
t295 = pkin(3) * t197;
t156 = pkin(4) * t186 + t295;
t318 = m(4) + m(6) + m(5) + m(3);
t345 = -m(6) * t156 - qJ(2) * t318 + mrSges(2,2) - mrSges(3,3);
t160 = -pkin(2) * t195 - pkin(6) * t194 - pkin(1);
t135 = qJDD(1) * t160 + qJDD(2);
t274 = t195 * t197;
t138 = qJD(1) * t160 + qJD(2);
t238 = qJ(2) * t257;
t97 = t138 * t197 + t200 * t238;
t50 = -qJD(3) * t97 + t200 * t135 - t170 * t274;
t21 = pkin(3) * t173 - pkin(7) * t136 + t50;
t123 = t200 * t138;
t240 = qJ(2) * t274;
t275 = t194 * t200;
t247 = pkin(7) * t275;
t213 = -t240 - t247;
t79 = qJD(1) * t213 + t123;
t62 = pkin(3) * t175 + t79;
t80 = -pkin(7) * t236 + t97;
t74 = t199 * t80;
t24 = t196 * t62 + t74;
t225 = qJD(3) * t240;
t255 = qJD(3) * t200;
t272 = t195 * t200;
t49 = -qJD(1) * t225 + t197 * t135 + t138 * t255 + t170 * t272;
t27 = pkin(7) * t137 + t49;
t6 = -qJD(4) * t24 - t196 * t27 + t199 * t21;
t2 = pkin(4) * t163 - qJ(5) * t51 - qJD(5) * t111 + t6;
t253 = qJD(4) * t199;
t254 = qJD(4) * t196;
t5 = t196 * t21 + t199 * t27 + t62 * t253 - t254 * t80;
t3 = qJ(5) * t52 - qJD(5) * t109 + t5;
t344 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2);
t343 = -t50 * mrSges(4,1) + t49 * mrSges(4,2);
t312 = m(5) * pkin(3);
t338 = m(5) * t295;
t36 = -mrSges(6,2) * t163 + mrSges(6,3) * t52;
t37 = -mrSges(5,2) * t163 + mrSges(5,3) * t52;
t330 = t37 + t36;
t286 = mrSges(6,3) * t109;
t83 = -mrSges(6,2) * t164 - t286;
t288 = mrSges(5,3) * t109;
t84 = -mrSges(5,2) * t164 - t288;
t329 = t83 + t84;
t285 = mrSges(6,3) * t111;
t85 = mrSges(6,1) * t164 - t285;
t287 = mrSges(5,3) * t111;
t86 = mrSges(5,1) * t164 - t287;
t328 = t86 + t85;
t327 = mrSges(4,1) + t312;
t103 = t111 * qJ(5);
t72 = t196 * t80;
t23 = t199 * t62 - t72;
t13 = -t103 + t23;
t319 = -qJD(3) - qJD(4);
t326 = (t257 + t319) * t218;
t95 = t319 * t152;
t325 = t195 * t215 + t95;
t174 = qJ(2) * t272;
t107 = t197 * t160 + t174;
t192 = t195 ^ 2;
t324 = t192 + t191;
t323 = t356 + t357 + t358;
t201 = cos(qJ(1));
t265 = t201 * t186;
t187 = cos(t193);
t198 = sin(qJ(1));
t268 = t198 * t187;
t126 = t195 * t265 - t268;
t271 = t195 * t201;
t127 = -t186 * t198 - t187 * t271;
t322 = t126 * t337 + t336 * t127;
t273 = t195 * t198;
t124 = t186 * t273 + t187 * t201;
t125 = t195 * t268 - t265;
t321 = -t124 * t337 + t336 * t125;
t320 = mrSges(5,1) * t186 - t187 * t336;
t317 = t192 * t355;
t283 = Ifges(4,4) * t200;
t284 = Ifges(4,4) * t197;
t315 = (qJ(2) * (mrSges(4,1) * t200 - mrSges(4,2) * t197) - t197 * (-Ifges(4,2) * t200 - t284) / 0.2e1 + t200 * (-Ifges(4,1) * t197 - t283) / 0.2e1) * t191;
t189 = t200 * pkin(3);
t157 = pkin(4) * t187 + t189;
t202 = -pkin(7) - pkin(6);
t223 = -mrSges(3,1) * t195 + mrSges(3,2) * t194;
t314 = -m(4) * t160 - m(6) * (-(pkin(2) + t157) * t195 - pkin(1)) - m(5) * (-(t189 + pkin(2)) * t195 - pkin(1)) + m(3) * pkin(1) - t223 + mrSges(2,1) + (mrSges(4,3) - m(6) * (-qJ(5) + t202) + mrSges(6,3) - m(5) * t202 + mrSges(5,3)) * t194;
t216 = (-t197 * Ifges(4,2) + t283) * t194;
t217 = (t200 * Ifges(4,1) - t284) * t194;
t313 = t200 * (t175 * Ifges(4,6) + qJD(1) * t216) + t197 * (t175 * Ifges(4,5) + qJD(1) * t217);
t203 = qJD(1) ^ 2;
t311 = m(6) * pkin(4);
t303 = t111 / 0.2e1;
t294 = pkin(3) * t199;
t293 = pkin(4) * t111;
t292 = g(1) * t194;
t29 = t199 * t79 - t72;
t149 = t200 * t160;
t87 = -t247 + t149 + (-qJ(2) * t197 - pkin(3)) * t195;
t276 = t194 * t197;
t98 = -pkin(7) * t276 + t107;
t39 = t196 * t87 + t199 * t98;
t279 = qJ(2) * t203;
t278 = qJ(5) * t109;
t185 = t194 * qJ(2);
t155 = t194 * t170;
t270 = t197 * t198;
t269 = t197 * t201;
t267 = t198 * t200;
t266 = t200 * t201;
t256 = qJD(2) * t195;
t260 = t160 * t255 + t200 * t256;
t259 = t346 * qJ(2);
t139 = pkin(3) * t236 + qJ(2) * t258;
t184 = t194 * qJD(2);
t233 = t194 * t255;
t146 = pkin(3) * t233 + t184;
t150 = pkin(3) * t276 + t185;
t249 = qJDD(1) * t194;
t244 = mrSges(4,3) * t276;
t239 = t347 + t348 + t349;
t234 = t197 * t256;
t232 = -t52 * mrSges(6,1) + t51 * mrSges(6,2);
t93 = -pkin(3) * t137 + t155;
t28 = -t196 * t79 - t74;
t38 = -t196 * t98 + t199 * t87;
t228 = pkin(3) * t235;
t227 = mrSges(4,3) * t236;
t226 = mrSges(4,3) * t235;
t224 = -mrSges(3,1) * t248 + mrSges(3,2) * t249;
t222 = mrSges(4,1) * t197 + mrSges(4,2) * t200;
t132 = -mrSges(4,2) * t175 - t227;
t133 = mrSges(4,1) * t175 - t226;
t219 = t132 * t200 - t133 * t197;
t142 = t195 * t269 - t267;
t140 = t195 * t270 + t266;
t75 = qJD(3) * t213 + t260;
t76 = -t234 + (-t174 + (pkin(7) * t194 - t160) * t197) * qJD(3);
t9 = t196 * t76 + t199 * t75 + t87 * t253 - t254 * t98;
t210 = t175 * t194 * (-Ifges(4,5) * t197 - Ifges(4,6) * t200);
t10 = -qJD(4) * t39 - t196 * t75 + t199 * t76;
t12 = pkin(4) * t164 + t13;
t14 = t24 - t278;
t77 = pkin(4) * t109 + qJD(5) + t139;
t204 = -t77 * (mrSges(6,1) * t111 - mrSges(6,2) * t109) - t139 * (mrSges(5,1) * t111 - mrSges(5,2) * t109) + t14 * t285 + t24 * t287 - t12 * t286 - t23 * t288 + t323 - (-t335 * t109 - t352) * t111 / 0.2e1 + t351 * t303 - (-t109 * t333 - t111 * t331) * t164 / 0.2e1 + (-t332 * t111 + t350 - t353) * t109 / 0.2e1 + t344;
t183 = -qJDD(1) * pkin(1) + qJDD(2);
t181 = pkin(4) + t294;
t177 = t191 * t279;
t144 = t222 * t194;
t143 = -t195 * t266 - t270;
t141 = t195 * t267 - t269;
t134 = t222 * t258;
t131 = t218 * t194;
t130 = t152 * t194;
t106 = t149 - t240;
t100 = -mrSges(4,2) * t173 + mrSges(4,3) * t137;
t99 = mrSges(4,1) * t173 - mrSges(4,3) * t136;
t96 = -t197 * t238 + t123;
t92 = t228 + t293;
t90 = -qJD(3) * t107 - t234;
t89 = -t225 + t260;
t88 = pkin(4) * t130 + t150;
t71 = (qJD(3) * t218 + t214) * t194;
t70 = t95 * t194;
t66 = mrSges(5,1) * t109 + mrSges(5,2) * t111;
t65 = mrSges(6,1) * t109 + mrSges(6,2) * t111;
t53 = -pkin(4) * t71 + t146;
t35 = mrSges(5,1) * t163 - mrSges(5,3) * t51;
t34 = mrSges(6,1) * t163 - mrSges(6,3) * t51;
t25 = -qJ(5) * t130 + t39;
t22 = -pkin(4) * t195 + qJ(5) * t131 + t38;
t18 = -pkin(4) * t52 + qJDD(5) + t93;
t17 = -t103 + t29;
t16 = t28 + t278;
t8 = -qJ(5) * t70 + qJD(5) * t131 + t10;
t7 = qJ(5) * t71 - qJD(5) * t130 + t9;
t1 = [(-t143 * mrSges(4,1) - t142 * mrSges(4,2) - t337 * t127 + t336 * t126 + t314 * t201 + (t338 - t345) * t198) * g(2) + (t141 * mrSges(4,1) - t140 * mrSges(4,2) + t124 * t336 + t125 * t337 + t198 * t314 + t201 * t345 - t269 * t312) * g(3) + t350 * t70 / 0.2e1 + t351 * t71 / 0.2e1 + (mrSges(5,1) * t93 + mrSges(6,1) * t18 - mrSges(5,3) * t5 - mrSges(6,3) * t3 - t163 * t331 - t332 * t52 - t334 * t51) * t130 + (Ifges(3,2) * t248 + Ifges(3,4) * t249 - t239 / 0.2e1 - t323 / 0.2e1 - t347 / 0.2e1 - t356 / 0.2e1 - t358 / 0.2e1 - t357 / 0.2e1 - t349 / 0.2e1 - t348 / 0.2e1 + t343 - t344) * t195 - pkin(1) * t224 + t183 * t223 + (-mrSges(4,1) * t137 + mrSges(4,2) * t136) * t185 + t134 * t184 + (t324 * t250 + t317 + t346) * mrSges(3,3) + (-t97 * t233 - t50 * t275) * mrSges(4,3) + t136 * t217 / 0.2e1 - (Ifges(4,4) * t136 + Ifges(4,2) * t137 + Ifges(4,6) * t173) * t276 / 0.2e1 + t144 * t155 + (Ifges(4,1) * t136 + Ifges(4,4) * t137 + Ifges(4,5) * t173) * t275 / 0.2e1 + m(4) * (t106 * t50 + t107 * t49 + t89 * t97 + t90 * t96 + t259) + t150 * (-mrSges(5,1) * t52 + mrSges(5,2) * t51) + t139 * (-mrSges(5,1) * t71 + mrSges(5,2) * t70) + t146 * t66 + t89 * t132 + t90 * t133 + t106 * t99 + t107 * t100 + t10 * t86 + t77 * (-mrSges(6,1) * t71 + mrSges(6,2) * t70) + t7 * t83 + t9 * t84 + t8 * t85 + (-t313 * qJD(3) / 0.2e1 + Ifges(3,4) * t248 + Ifges(3,1) * t249 + t173 * (Ifges(4,5) * t200 - Ifges(4,6) * t197) / 0.2e1) * t194 + t53 * t65 + t38 * t35 + t39 * t37 + t22 * t34 + t25 * t36 + (-mrSges(5,2) * t93 - mrSges(6,2) * t18 + mrSges(5,3) * t6 + mrSges(6,3) * t2 - t163 * t333 - t334 * t52 - t335 * t51) * t131 + m(5) * (t10 * t23 + t139 * t146 + t150 * t93 + t24 * t9 + t38 * t6 + t39 * t5) + m(6) * (t12 * t8 + t14 * t7 + t18 * t88 + t2 * t22 + t25 * t3 + t53 * t77) + (t96 * t244 + t210 / 0.2e1) * qJD(3) + t315 * t251 + m(3) * (-pkin(1) * t183 + qJ(2) * t317 + t259) + t137 * t216 / 0.2e1 + (t331 * t71 + t333 * t70) * t164 / 0.2e1 - (t332 * t71 + t334 * t70) * t109 / 0.2e1 + (t334 * t71 + t335 * t70) * t303 - t49 * t244 + (-t12 * t70 + t14 * t71) * mrSges(6,3) + (-t23 * t70 + t24 * t71) * mrSges(5,3) + t88 * t232 + Ifges(2,3) * qJDD(1); t224 + t197 * t100 - t324 * t203 * mrSges(3,3) + t330 * t152 - (t34 + t35) * t218 + t219 * qJD(3) + (-t219 * t195 + (-t134 - t65 - t66) * t194) * qJD(1) + t200 * t99 + t326 * t329 + t325 * t328 - (g(2) * t201 + g(3) * t198) * t318 + (t12 * t325 + t14 * t326 + t152 * t3 - t2 * t218 - t258 * t77) * m(6) + (-t139 * t258 + t152 * t5 - t218 * t6 + t23 * t325 + t24 * t326) * m(5) + (t197 * t49 + t200 * t50 - t177 + t175 * (-t96 * t197 + t97 * t200)) * m(4) + (-t192 * t279 - t177 + t183) * m(3); t204 + (-t12 * t16 - t14 * t17 + t156 * t292 + t181 * t2 - t77 * t92) * m(6) + ((t196 * t3 + (-t12 * t196 + t14 * t199) * qJD(4)) * m(6) - t328 * t254 + t329 * t253 + t330 * t196) * pkin(3) + (-t132 - t227) * t96 + t35 * t294 + (-m(6) * (-t156 * t271 + t157 * t198) - mrSges(4,2) * t143 + t327 * t142 + t322) * g(3) + (-m(6) * (t156 * t273 + t157 * t201) - mrSges(4,2) * t141 - t327 * t140 + t321) * g(2) - t343 + (mrSges(6,1) * t186 + t320 + t338) * t292 + (t133 + t226) * t97 + t239 + t181 * t34 + g(1) * t144 - t92 * t65 - t17 * t83 - t29 * t84 - t16 * t85 - t28 * t86 - t315 * t203 + t313 * t258 / 0.2e1 - t66 * t228 - m(5) * (t139 * t228 + t23 * t28 + t24 * t29) - qJD(1) * t210 / 0.2e1 + (t196 * t5 + t199 * t6 + (-t196 * t23 + t199 * t24) * qJD(4)) * t312; t204 - t13 * t83 - t23 * t84 + t14 * t85 + t24 * t86 - t65 * t293 + t2 * t311 + pkin(4) * t34 - m(6) * (t77 * t293 + (-t12 + t13) * t14) + (-(-mrSges(6,1) - t311) * t186 + t320) * t292 + (t126 * t311 + t322) * g(3) + (-t124 * t311 + t321) * g(2); t109 * t83 + t111 * t85 + (g(1) * t195 + t109 * t14 + t111 * t12 + t18 + (g(2) * t198 - g(3) * t201) * t194) * m(6) + t232;];
tau = t1;
