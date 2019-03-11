% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:47
% EndTime: 2019-03-09 02:59:06
% DurationCPUTime: 11.69s
% Computational Cost: add. (3022->556), mult. (5442->669), div. (0->0), fcn. (2571->6), ass. (0->257)
t353 = Ifges(4,1) + Ifges(5,1);
t352 = -Ifges(4,4) + Ifges(5,5);
t343 = -Ifges(6,5) + Ifges(5,6);
t351 = -Ifges(4,6) + t343;
t141 = sin(qJ(3));
t144 = cos(qJ(3));
t242 = qJD(1) * qJD(3);
t91 = -t144 * qJDD(1) + t141 * t242;
t348 = t91 / 0.2e1;
t92 = qJDD(1) * t141 + t144 * t242;
t311 = t92 / 0.2e1;
t344 = Ifges(5,4) + Ifges(4,5);
t279 = qJDD(3) / 0.2e1;
t290 = Ifges(5,5) * t141;
t296 = Ifges(4,4) * t141;
t350 = t353 * t144 + t290 - t296;
t140 = sin(qJ(6));
t143 = cos(qJ(6));
t131 = t144 * pkin(5);
t309 = -pkin(4) - pkin(8);
t136 = -pkin(3) + t309;
t173 = t136 * t141 - qJ(2);
t247 = t144 * qJD(1);
t117 = qJ(4) * t247;
t245 = qJD(5) + t117;
t31 = (t173 + t131) * qJD(1) + t245;
t147 = -pkin(1) - pkin(7);
t105 = qJD(1) * t147 + qJD(2);
t51 = (qJ(5) * qJD(1) + t105) * t144;
t162 = qJD(4) - t51;
t36 = qJD(3) * t136 + t162;
t10 = t140 * t31 + t143 * t36;
t104 = qJDD(1) * t147 + qJDD(2);
t241 = qJD(1) * qJD(5);
t256 = qJD(3) * t141;
t84 = t105 * t256;
t151 = qJ(5) * t91 + qJDD(4) + t84 + (-t104 - t241) * t144;
t11 = qJDD(3) * t136 + t151;
t243 = qJD(1) * qJD(2);
t106 = qJDD(1) * qJ(2) + t243;
t124 = t144 * qJD(4);
t22 = t92 * pkin(3) + t91 * qJ(4) - qJD(1) * t124 + t106;
t152 = qJDD(5) - t22;
t7 = -pkin(5) * t91 + t309 * t92 + t152;
t9 = -t140 * t36 + t143 * t31;
t1 = qJD(6) * t9 + t11 * t143 + t140 * t7;
t2 = -qJD(6) * t10 - t11 * t140 + t143 * t7;
t201 = t1 * t143 - t140 * t2;
t250 = qJD(6) * t143;
t252 = qJD(6) * t140;
t349 = t10 * t252 + t9 * t250 - t201;
t347 = -t92 / 0.2e1;
t346 = -m(6) - m(7);
t56 = -qJDD(3) * mrSges(6,1) - mrSges(6,3) * t92;
t255 = qJD(3) * t143;
t257 = qJD(1) * t141;
t82 = -t140 * t257 - t255;
t28 = qJD(6) * t82 - qJDD(3) * t140 + t143 * t92;
t248 = t140 * qJD(3);
t83 = t143 * t257 - t248;
t29 = -qJD(6) * t83 - qJDD(3) * t143 - t140 * t92;
t8 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t345 = t8 - t56;
t231 = mrSges(6,3) * t257;
t297 = -qJD(3) * mrSges(6,1) + mrSges(7,1) * t82 - mrSges(7,2) * t83 - t231;
t58 = -mrSges(5,2) * t92 + qJDD(3) * mrSges(5,3);
t342 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t92 + t58;
t192 = mrSges(7,1) * t140 + mrSges(7,2) * t143;
t116 = qJ(5) * t257;
t93 = t141 * t105;
t67 = qJD(3) * qJ(4) + t93;
t46 = -t116 - t67;
t40 = qJD(3) * pkin(5) - t46;
t341 = t192 * t40;
t228 = mrSges(5,2) * t257;
t102 = qJD(3) * mrSges(5,3) - t228;
t98 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t257;
t339 = t102 + t98;
t146 = -pkin(3) - pkin(4);
t16 = qJDD(3) * t146 + t151;
t338 = -qJD(3) * t46 - t16;
t38 = t104 * t144 - t84;
t33 = -qJDD(3) * pkin(3) + qJDD(4) - t38;
t337 = qJD(3) * t67 - t33;
t295 = Ifges(4,4) * t144;
t186 = -t141 * Ifges(4,2) + t295;
t108 = qJD(6) + t247;
t305 = Ifges(7,4) * t83;
t24 = Ifges(7,2) * t82 + Ifges(7,6) * t108 + t305;
t336 = Ifges(4,6) * qJD(3) + qJD(1) * t186 + t140 * t24;
t232 = mrSges(5,2) * t247;
t335 = -mrSges(4,3) * t247 - t232 + (mrSges(4,1) + mrSges(5,1)) * qJD(3);
t127 = t144 * qJ(4);
t260 = t127 + t131;
t129 = t141 * mrSges(6,2);
t259 = t144 * mrSges(6,1) + t129;
t334 = t141 * (Ifges(5,3) * t144 - t290) + t144 * (-Ifges(4,1) * t141 - t295);
t293 = Ifges(6,4) * t144;
t333 = t141 * (-Ifges(4,2) * t144 - t296) + t144 * (Ifges(6,2) * t141 + t293);
t332 = (-Ifges(4,6) + Ifges(5,6)) * t144 - t344 * t141;
t254 = qJD(3) * t144;
t39 = t141 * t104 + t105 * t254;
t32 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t39;
t330 = t141 * t32 - t144 * t33;
t80 = qJDD(6) - t91;
t17 = mrSges(7,1) * t80 - mrSges(7,3) * t28;
t18 = -mrSges(7,2) * t80 + mrSges(7,3) * t29;
t177 = t140 * t17 - t143 * t18;
t329 = -t141 * t39 - t144 * t38;
t145 = cos(qJ(1));
t301 = g(2) * t145;
t142 = sin(qJ(1));
t302 = g(1) * t142;
t328 = t301 - t302;
t327 = t83 * Ifges(7,5) + t82 * Ifges(7,6) + t108 * Ifges(7,3) + qJD(1) * t350 + qJD(3) * t344;
t118 = Ifges(5,5) * t247;
t189 = t141 * Ifges(6,1) - t293;
t75 = Ifges(7,4) * t82;
t25 = Ifges(7,1) * t83 + Ifges(7,5) * t108 + t75;
t326 = Ifges(5,3) * t257 + qJD(1) * t189 + qJD(3) * t343 + t143 * t25 + t118;
t325 = -mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t324 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t42 = -mrSges(7,2) * t108 + mrSges(7,3) * t82;
t43 = mrSges(7,1) * t108 - mrSges(7,3) * t83;
t175 = t140 * t42 + t143 * t43;
t53 = qJDD(3) * mrSges(6,2) + t91 * mrSges(6,3);
t55 = -qJDD(3) * mrSges(5,1) - t91 * mrSges(5,2);
t323 = t175 * qJD(6) + t177 - t53 - t55;
t197 = mrSges(4,1) * t141 + mrSges(4,2) * t144;
t222 = -m(6) * qJ(4) - mrSges(6,1);
t286 = t141 * mrSges(5,1);
t322 = mrSges(2,2) - mrSges(3,3) - m(7) * (pkin(8) * t141 - t260) - t141 * mrSges(7,3) + t129 - t286 - t197 + (m(5) * qJ(4) + mrSges(5,3) - t222) * t144;
t321 = -m(6) * t46 + m(7) * t40 - t297;
t320 = qJD(1) ^ 2;
t318 = t28 / 0.2e1;
t317 = t29 / 0.2e1;
t316 = t80 / 0.2e1;
t315 = -t82 / 0.2e1;
t314 = -t83 / 0.2e1;
t313 = t83 / 0.2e1;
t310 = -m(3) - m(4);
t308 = -t108 / 0.2e1;
t304 = pkin(1) * t142;
t303 = pkin(3) * t141;
t139 = qJ(4) + pkin(5);
t294 = Ifges(6,4) * t141;
t292 = Ifges(7,4) * t140;
t291 = Ifges(7,4) * t143;
t289 = Ifges(5,5) * t144;
t278 = qJ(4) * t141;
t277 = qJD(3) * t40;
t274 = t105 * t144;
t273 = t140 * t141;
t272 = t140 * t144;
t271 = t141 * t142;
t270 = t141 * t143;
t269 = t141 * t145;
t267 = t142 * t144;
t266 = t143 * t144;
t265 = t144 * t145;
t263 = qJ(5) + t147;
t262 = pkin(3) * t267 + qJ(4) * t271;
t128 = t145 * qJ(2);
t261 = pkin(3) * t269 + t128;
t258 = t145 * pkin(1) + t142 * qJ(2);
t253 = qJD(3) * t147;
t251 = qJD(6) * t141;
t249 = qJDD(1) * mrSges(3,2);
t246 = -qJD(2) + t124;
t244 = m(5) - t346;
t236 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t80;
t235 = -t102 + t297;
t233 = t146 * t141;
t230 = mrSges(6,3) * t247;
t226 = t145 * pkin(7) + t258;
t221 = -t91 * mrSges(6,1) + t92 * mrSges(6,2);
t215 = -t250 / 0.2e1;
t212 = t127 - t303;
t211 = -t242 / 0.2e1;
t210 = t242 / 0.2e1;
t208 = (t106 + t243) * qJ(2);
t207 = qJD(3) * t263;
t205 = pkin(3) * t271 + t226;
t204 = pkin(4) * t269 + t142 * qJ(5) + t261;
t203 = t127 + t233;
t202 = g(1) * t145 + g(2) * t142;
t200 = t10 * t140 + t9 * t143;
t199 = m(7) * t136 - mrSges(7,3);
t198 = mrSges(4,1) * t144 - mrSges(4,2) * t141;
t196 = t144 * mrSges(5,1) + t141 * mrSges(5,3);
t195 = -t144 * mrSges(5,3) + t286;
t194 = t141 * mrSges(6,1) - t144 * mrSges(6,2);
t193 = t143 * mrSges(7,1) - t140 * mrSges(7,2);
t188 = Ifges(7,1) * t143 - t292;
t187 = -Ifges(7,1) * t140 - t291;
t184 = -Ifges(7,2) * t140 + t291;
t183 = -Ifges(7,2) * t143 - t292;
t181 = Ifges(6,5) * t144 + Ifges(6,6) * t141;
t180 = Ifges(7,5) * t143 - Ifges(7,6) * t140;
t179 = -Ifges(7,5) * t140 - Ifges(7,6) * t143;
t178 = pkin(3) * t144 + t278;
t176 = -t140 * t43 + t143 * t42;
t44 = t173 + t260;
t95 = t263 * t144;
t20 = t140 * t95 + t143 * t44;
t21 = t140 * t44 - t143 * t95;
t59 = -qJD(3) * pkin(3) + qJD(4) - t274;
t174 = t59 * t141 + t67 * t144;
t99 = qJD(3) * mrSges(6,2) - t230;
t172 = t176 + t99;
t171 = t144 * t146 - t278;
t45 = (-qJ(2) + t233) * qJD(1) + t245;
t170 = t45 * t194;
t60 = -t117 + (qJ(2) + t303) * qJD(1);
t169 = t60 * t196;
t168 = qJ(2) * t198;
t87 = t259 * qJD(1);
t88 = t195 * qJD(1);
t163 = t175 + t87 - t88;
t161 = -t140 * t251 + t143 * t254;
t160 = t141 * t250 + t144 * t248;
t156 = t136 * t144 - t139 * t141;
t155 = -Ifges(7,5) * t141 + t144 * t188;
t154 = -Ifges(7,6) * t141 + t144 * t184;
t153 = -Ifges(7,3) * t141 + t144 * t180;
t19 = -qJ(5) * t92 - t141 * t241 - t32;
t149 = -qJD(6) * t200 + t201;
t120 = -qJDD(1) * pkin(1) + qJDD(2);
t119 = Ifges(6,4) * t257;
t96 = qJ(2) - t212;
t89 = t197 * qJD(1);
t86 = t178 * qJD(1);
t74 = t192 * t141;
t73 = -qJ(2) + t203;
t71 = t140 * t142 - t143 * t265;
t70 = t140 * t265 + t142 * t143;
t69 = t140 * t145 + t142 * t266;
t68 = t140 * t267 - t143 * t145;
t62 = -Ifges(6,2) * t247 - Ifges(6,6) * qJD(3) + t119;
t54 = qJDD(3) * mrSges(4,1) + mrSges(4,3) * t91;
t52 = t171 * qJD(1);
t50 = t93 + t116;
t49 = qJD(3) * t178 - t246;
t47 = -qJD(5) * t144 + t141 * t207;
t41 = qJD(3) * t171 + t246;
t37 = qJD(3) * t146 + t162;
t35 = t156 * qJD(1);
t30 = qJD(3) * t156 + t246;
t15 = t140 * t35 + t143 * t50;
t14 = -t140 * t50 + t143 * t35;
t13 = qJDD(3) * pkin(5) - t19;
t12 = -pkin(4) * t92 + t152;
t6 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + t80 * Ifges(7,5);
t5 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + t80 * Ifges(7,6);
t4 = -qJD(6) * t21 - t140 * t47 + t143 * t30;
t3 = qJD(6) * t20 + t140 * t30 + t143 * t47;
t23 = [(t188 * t318 + t180 * t316 + t184 * t317 - t19 * mrSges(6,3) + t24 * t215 - t335 * t253 + t342 * t147 + (-m(6) * t19 + m(7) * t13 + t345) * t263 + t321 * qJD(5) + ((-Ifges(5,1) + Ifges(6,1)) * t144 + t294) * t210 + (Ifges(6,1) + 0.2e1 * Ifges(5,3) + Ifges(4,2)) * t311 + (Ifges(6,4) - t352) * t348 - (Ifges(6,5) + Ifges(4,6)) * qJDD(3) / 0.2e1 + (Ifges(5,6) + t351) * t279) * t141 + m(6) * (t12 * t73 - t16 * t95 + t37 * t47 + t41 * t45) - (qJD(6) * t25 + t5) * t273 / 0.2e1 + m(7) * (t1 * t21 + t10 * t3 + t2 * t20 + t4 * t9) + (t332 / 0.2e1 - t181 / 0.2e1) * qJD(3) ^ 2 + (t197 + 0.2e1 * mrSges(3,3)) * t106 + (t37 * mrSges(6,3) - t59 * mrSges(5,2) - t327 / 0.2e1 + t62 / 0.2e1 - t9 * mrSges(7,1) + t10 * mrSges(7,2)) * t256 + t186 * t347 + t294 * t348 - t330 * mrSges(5,2) + (-t46 * mrSges(6,3) - t336 / 0.2e1 - t67 * mrSges(5,2) + t326 / 0.2e1) * t254 + (t289 + t189) * t311 + t329 * mrSges(4,3) + m(4) * (-t329 * t147 + t208) + m(5) * (t22 * t96 + t49 * t60 + (qJD(3) * t174 + t330) * t147) + t333 * t211 + t334 * t210 + (-m(3) * t258 - m(4) * t226 - m(5) * t205 + t69 * mrSges(7,1) - t68 * mrSges(7,2) + t346 * (pkin(4) * t271 - qJ(5) * t145 + t205) + t325 * t145 + t322 * t142) * g(2) + (-m(5) * t261 - m(6) * t204 - m(4) * t128 - m(3) * (t128 - t304) - m(7) * (t204 - t304) - t71 * mrSges(7,1) - t70 * mrSges(7,2) + (m(7) * pkin(7) + (-m(4) - m(5) - m(6)) * t147 - t325) * t142 + t322 * t145) * g(1) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + qJD(3) * t169 - qJD(3) * t170 + t6 * t270 / 0.2e1 + t12 * t259 - pkin(1) * t249 + t82 * (qJD(3) * t154 + t183 * t251) / 0.2e1 + t108 * (qJD(3) * t153 + t179 * t251) / 0.2e1 + t73 * t221 + m(3) * (-pkin(1) * t120 + t208) + t120 * mrSges(3,2) + t47 * t99 + qJ(2) * (mrSges(4,1) * t92 - mrSges(4,2) * t91) - t95 * t53 + t96 * (mrSges(5,1) * t92 + mrSges(5,3) * t91) + t41 * t87 + t49 * t88 + qJD(2) * t89 + t13 * t74 + ((-t55 + t54) * t147 - t16 * mrSges(6,3) + Ifges(6,4) * t347 + Ifges(7,5) * t318 - Ifges(6,2) * t91 + Ifges(7,6) * t317 + Ifges(7,3) * t316 + t289 * t210 + t344 * t279 + t324 + t236 / 0.2e1 + t352 * t311 - t353 * t348 + t339 * t253 + t321 * t207 + (Ifges(6,6) + t344 / 0.2e1) * qJDD(3)) * t144 + (-t1 * t273 - t10 * t160 - t161 * t9 - t2 * t270) * mrSges(7,3) + t3 * t42 + t4 * t43 + t20 * t17 + t21 * t18 + (qJD(3) * t155 + t187 * t251) * t313 + t22 * t195 - t350 * t91 / 0.2e1 + t40 * (mrSges(7,1) * t160 + mrSges(7,2) * t161) + t168 * t242; m(3) * t120 + t249 + (qJ(2) * t310 - mrSges(3,3)) * t320 + ((t172 - t335) * qJD(3) + m(7) * (t10 * t255 - t248 * t9 + t13) + m(5) * (qJD(3) * t59 + t32) + m(6) * (qJD(3) * t37 - t19) + m(4) * t39 + t342 + t345) * t141 + (t54 + (t98 - t235) * qJD(3) + m(7) * (t277 + t349) + m(5) * t337 + m(6) * t338 + m(4) * t38 + t323) * t144 + (-m(5) * t60 + m(6) * t45 + m(7) * t200 + t163 - t89) * qJD(1) + t328 * (t244 - t310); t351 * t92 + (-m(5) * t212 - m(6) * t203 - m(7) * t260 - t141 * t199 - t144 * t193 + t195 + t197 - t259) * g(3) + (-t9 * (-mrSges(7,1) * t141 - mrSges(7,3) * t266) - t10 * (mrSges(7,2) * t141 - mrSges(7,3) * t272) - t169 + t170) * qJD(1) - (Ifges(6,1) * t247 + t119 + t62) * t257 / 0.2e1 - (t108 * t153 + t154 * t82 + t155 * t83) * qJD(1) / 0.2e1 + t349 * mrSges(7,3) + (-t56 + t58) * qJ(4) - t247 * t341 + (t180 * t308 + t184 * t315 + t188 * t314 - t341) * qJD(6) + t336 * t247 / 0.2e1 - t339 * t274 + t335 * t93 + (m(7) * t149 - t250 * t43 - t252 * t42 - t177) * t136 + t332 * t211 - (-Ifges(5,1) * t257 + t118 + t326) * t247 / 0.2e1 + t327 * t257 / 0.2e1 + (-m(5) * t262 + t346 * (pkin(4) * t267 + t262) + (-t196 - t194 - (m(7) * pkin(8) + mrSges(7,3)) * t144 - (m(7) * pkin(5) + t193) * t141) * t142) * g(1) + t297 * t51 + (-Ifges(6,6) - t344) * t91 + (-t10 * t15 + t13 * t139 - t14 * t9 - t40 * t51) * m(7) + (-pkin(3) * t33 + qJ(4) * t32 - t105 * t174 - t60 * t86) * m(5) - t198 * t302 + t25 * t215 + t181 * t210 + t24 * t252 / 0.2e1 - t37 * t231 + t146 * t53 - t143 * t5 / 0.2e1 - t140 * t6 / 0.2e1 + t139 * t8 - t50 * t99 - t52 * t87 - t86 * t88 + t59 * t228 + t46 * t230 + t67 * t232 + (-qJ(4) * t19 + t146 * t16 - t37 * t50 - t45 * t52 + t46 * t51) * m(6) + (Ifges(6,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(3) + (t333 / 0.2e1 - t334 / 0.2e1 - t168) * t320 - pkin(3) * t55 + t38 * mrSges(4,1) - t39 * mrSges(4,2) - t15 * t42 - t14 * t43 + t32 * mrSges(5,3) - t33 * mrSges(5,1) + t16 * mrSges(6,2) - t19 * mrSges(6,1) + (m(5) * t178 + t196 + t198 + (-m(6) * t146 - mrSges(6,2) - t199) * t144 + (m(7) * t139 + t193 - t222) * t141) * t301 + (m(5) * t67 + t102 + t321) * qJD(4) + t187 * t318 + t179 * t316 + t183 * t317 + t13 * t193; -t163 * t247 + t235 * qJD(3) + (-t141 * g(3) - t144 * t328) * t244 + (-t200 * t247 + t149 - t277) * m(7) + (-t247 * t45 - t338) * m(6) + (t247 * t60 - t337) * m(5) - t323; t140 * t18 + t143 * t17 + t176 * qJD(6) + (t1 * t140 + t2 * t143 + (t10 * t143 - t9 * t140) * qJD(6) + t202) * m(7) + (t12 + t202) * m(6) + (t297 * t141 + t172 * t144 - m(6) * (-t141 * t46 - t144 * t37) - m(7) * (-t10 * t266 + t141 * t40 + t272 * t9)) * qJD(1) + t221; -t40 * (mrSges(7,1) * t83 + mrSges(7,2) * t82) + (Ifges(7,1) * t82 - t305) * t314 + t24 * t313 + (Ifges(7,5) * t82 - Ifges(7,6) * t83) * t308 - t9 * t42 + t10 * t43 - g(1) * (mrSges(7,1) * t68 + mrSges(7,2) * t69) - g(2) * (-mrSges(7,1) * t70 + mrSges(7,2) * t71) + g(3) * t74 + (t10 * t83 + t82 * t9) * mrSges(7,3) + t236 + (-Ifges(7,2) * t83 + t25 + t75) * t315 + t324;];
tau  = t23;
