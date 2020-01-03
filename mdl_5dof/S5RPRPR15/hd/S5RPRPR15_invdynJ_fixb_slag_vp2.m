% Calculate vector of inverse dynamics joint torques for
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR15_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:38
% EndTime: 2019-12-31 18:36:59
% DurationCPUTime: 11.95s
% Computational Cost: add. (3272->470), mult. (6609->652), div. (0->0), fcn. (4114->10), ass. (0->221)
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t221 = qJD(1) * qJD(3);
t132 = qJDD(1) * t160 + t163 * t221;
t263 = t132 / 0.2e1;
t156 = sin(pkin(8));
t157 = cos(pkin(8));
t131 = qJDD(1) * t163 - t160 * t221;
t198 = -qJD(4) * t163 + qJD(2);
t220 = qJDD(1) * qJ(2);
t48 = pkin(3) * t132 - qJ(4) * t131 + qJD(1) * t198 + t220;
t165 = -pkin(1) - pkin(6);
t142 = qJDD(1) * t165 + qJDD(2);
t143 = qJD(1) * t165 + qJD(2);
t225 = qJD(3) * t163;
t76 = t160 * t142 + t143 * t225;
t66 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t76;
t24 = -t156 * t66 + t157 * t48;
t87 = qJDD(3) * t156 + t131 * t157;
t10 = pkin(4) * t132 - pkin(7) * t87 + t24;
t25 = t156 * t48 + t157 * t66;
t86 = qJDD(3) * t157 - t131 * t156;
t11 = pkin(7) * t86 + t25;
t159 = sin(qJ(5));
t162 = cos(qJ(5));
t227 = qJD(1) * t163;
t122 = qJD(3) * t156 + t157 * t227;
t228 = qJD(1) * t160;
t134 = pkin(3) * t160 - qJ(4) * t163 + qJ(2);
t112 = t134 * qJD(1);
t133 = t160 * t143;
t115 = qJD(3) * qJ(4) + t133;
t50 = t157 * t112 - t115 * t156;
t31 = pkin(4) * t228 - pkin(7) * t122 + t50;
t121 = qJD(3) * t157 - t156 * t227;
t51 = t156 * t112 + t157 * t115;
t32 = pkin(7) * t121 + t51;
t8 = -t159 * t32 + t162 * t31;
t1 = qJD(5) * t8 + t10 * t159 + t11 * t162;
t9 = t159 * t31 + t162 * t32;
t2 = -qJD(5) * t9 + t10 * t162 - t11 * t159;
t320 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t195 = mrSges(4,1) * t160 + mrSges(4,2) * t163;
t253 = pkin(7) + qJ(4);
t311 = m(5) * qJ(4) + m(6) * t253 + mrSges(5,3) + mrSges(6,3);
t319 = -t163 * t311 + t195;
t145 = qJD(5) + t228;
t262 = -t145 / 0.2e1;
t63 = t121 * t159 + t122 * t162;
t270 = -t63 / 0.2e1;
t199 = t162 * t121 - t122 * t159;
t272 = -t199 / 0.2e1;
t318 = Ifges(6,5) * t270 + Ifges(6,6) * t272 + Ifges(6,3) * t262;
t317 = t160 / 0.2e1;
t316 = t163 / 0.2e1;
t308 = qJD(3) / 0.2e1;
t187 = Ifges(5,5) * t157 - Ifges(5,6) * t156;
t196 = mrSges(4,1) * t163 - mrSges(4,2) * t160;
t251 = Ifges(4,4) * t163;
t313 = (-Ifges(4,1) * t160 - t251) * t316 + (Ifges(5,3) * t163 - t160 * t187) * t317 + qJ(2) * t196;
t146 = pkin(4) * t157 + pkin(3);
t194 = -mrSges(5,1) * t157 + mrSges(5,2) * t156;
t312 = -m(5) * pkin(3) - m(6) * t146 + t194;
t190 = -t160 * Ifges(4,2) + t251;
t310 = t9 * mrSges(6,2) + Ifges(4,6) * t308 + qJD(1) * t190 / 0.2e1 - t122 * Ifges(5,5) / 0.2e1 - t121 * Ifges(5,6) / 0.2e1 - Ifges(5,3) * t228 / 0.2e1 - t8 * mrSges(6,1) + t318;
t161 = sin(qJ(1));
t164 = cos(qJ(1));
t309 = g(1) * t161 - g(2) * t164;
t265 = -m(5) - m(6);
t306 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t121 + mrSges(5,2) * t122 + mrSges(4,3) * t227;
t249 = Ifges(5,4) * t157;
t189 = -Ifges(5,2) * t156 + t249;
t250 = Ifges(5,4) * t156;
t191 = Ifges(5,1) * t157 - t250;
t305 = t121 * (Ifges(5,6) * t163 - t160 * t189) + t122 * (Ifges(5,5) * t163 - t160 * t191);
t222 = qJD(1) * qJD(2);
t144 = t220 + t222;
t226 = qJD(3) * t160;
t75 = t142 * t163 - t143 * t226;
t183 = t160 * t76 + t163 * t75;
t303 = -m(4) + t265;
t193 = t156 * mrSges(5,1) + t157 * mrSges(5,2);
t235 = t157 * t160;
t237 = t143 * t163;
t104 = -qJD(3) * pkin(3) + qJD(4) - t237;
t239 = t104 * t160;
t302 = t193 * t239 - t51 * (mrSges(5,3) * t156 * t160 - mrSges(5,2) * t163) - t50 * (mrSges(5,1) * t163 + mrSges(5,3) * t235);
t301 = mrSges(2,1) + mrSges(4,3) + t193 - mrSges(3,2);
t300 = t160 * t312 + mrSges(2,2) - mrSges(3,3) - t319;
t267 = t87 / 0.2e1;
t299 = Ifges(5,1) * t267 + Ifges(5,5) * t263;
t19 = qJD(5) * t199 + t159 * t86 + t162 * t87;
t277 = t19 / 0.2e1;
t20 = -qJD(5) * t63 - t159 * t87 + t162 * t86;
t276 = t20 / 0.2e1;
t268 = t86 / 0.2e1;
t123 = qJDD(5) + t132;
t264 = t123 / 0.2e1;
t35 = -t86 * mrSges(5,1) + t87 * mrSges(5,2);
t5 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t297 = -t35 - t5;
t182 = t156 * t159 - t157 * t162;
t218 = pkin(7) * t235;
t186 = pkin(3) * t163 + qJ(4) * t160;
t128 = t186 * qJD(1);
t236 = t156 * t163;
t70 = t157 * t128 - t143 * t236;
t40 = (pkin(4) * t163 + t218) * qJD(1) + t70;
t214 = t156 * t228;
t234 = t157 * t163;
t71 = t156 * t128 + t143 * t234;
t49 = pkin(7) * t214 + t71;
t135 = t253 * t156;
t136 = t253 * t157;
t73 = -t135 * t162 - t136 * t159;
t296 = -qJD(4) * t182 + qJD(5) * t73 - t159 * t40 - t162 * t49;
t125 = t156 * t162 + t157 * t159;
t74 = -t135 * t159 + t136 * t162;
t295 = -qJD(4) * t125 - qJD(5) * t74 + t159 * t49 - t162 * t40;
t281 = qJD(5) * t182;
t97 = t125 * t163;
t294 = t182 * qJD(1) - qJD(3) * t97 + t160 * t281;
t107 = t125 * qJD(1);
t109 = t125 * qJD(5);
t99 = t182 * t163;
t293 = -qJD(3) * t99 - t109 * t160 - t107;
t288 = t182 * t160;
t85 = qJD(1) * t288;
t292 = t281 + t85;
t84 = t160 * t107;
t291 = t109 + t84;
t289 = t163 * t309;
t155 = pkin(8) + qJ(5);
t148 = sin(t155);
t149 = cos(t155);
t287 = mrSges(6,1) * t149 - mrSges(6,2) * t148 - t312;
t82 = -mrSges(5,2) * t228 + mrSges(5,3) * t121;
t83 = mrSges(5,1) * t228 - mrSges(5,3) * t122;
t285 = -t156 * t83 + t157 * t82;
t52 = -mrSges(5,2) * t132 + mrSges(5,3) * t86;
t53 = mrSges(5,1) * t132 - mrSges(5,3) * t87;
t284 = -t156 * t53 + t157 * t52;
t185 = -t156 * t24 + t157 * t25;
t258 = pkin(4) * t156;
t282 = -t258 + t165;
t280 = qJD(1) ^ 2;
t279 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t264;
t278 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t264;
t259 = Ifges(6,4) * t63;
t22 = Ifges(6,2) * t199 + Ifges(6,6) * t145 + t259;
t275 = t22 / 0.2e1;
t58 = Ifges(6,4) * t199;
t23 = Ifges(6,1) * t63 + Ifges(6,5) * t145 + t58;
t274 = t23 / 0.2e1;
t273 = Ifges(5,4) * t268 + t299;
t271 = t199 / 0.2e1;
t269 = t63 / 0.2e1;
t266 = -m(3) - m(4);
t261 = t145 / 0.2e1;
t252 = Ifges(4,4) * t160;
t72 = -qJDD(3) * pkin(3) + qJDD(4) - t75;
t241 = t163 * t72;
t100 = qJD(3) * t186 + t198;
t224 = qJD(3) * t165;
t211 = t163 * t224;
t65 = t156 * t100 + t157 * t211;
t233 = t160 * t161;
t232 = t160 * t164;
t231 = t160 * t165;
t78 = t156 * t134 + t157 * t231;
t229 = t164 * pkin(1) + t161 * qJ(2);
t223 = qJDD(1) * mrSges(3,2);
t217 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t123;
t213 = t156 * t226;
t208 = t228 / 0.2e1;
t204 = -t156 * t165 + pkin(4);
t202 = -t221 / 0.2e1;
t200 = (t144 + t222) * qJ(2);
t192 = t163 * Ifges(4,1) - t252;
t188 = -Ifges(4,5) * t160 - Ifges(4,6) * t163;
t184 = -t156 * t50 + t157 * t51;
t119 = t157 * t134;
t59 = -pkin(7) * t234 + t160 * t204 + t119;
t69 = -pkin(7) * t236 + t78;
t26 = -t159 * t69 + t162 * t59;
t27 = t159 * t59 + t162 * t69;
t180 = t160 * (-Ifges(4,2) * t163 - t252);
t147 = -qJDD(1) * pkin(1) + qJDD(2);
t138 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t228;
t129 = t195 * qJD(1);
t120 = t282 * t163;
t114 = Ifges(4,5) * qJD(3) + qJD(1) * t192;
t105 = t282 * t226;
t102 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t132;
t101 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t131;
t96 = t125 * t160;
t92 = -t148 * t161 + t149 * t232;
t91 = t148 * t232 + t149 * t161;
t90 = t148 * t164 + t149 * t233;
t89 = -t148 * t233 + t149 * t164;
t88 = -pkin(4) * t214 + t133;
t81 = t157 * t100;
t77 = -t156 * t231 + t119;
t67 = -pkin(4) * t121 + t104;
t64 = -t156 * t211 + t81;
t56 = t122 * Ifges(5,1) + t121 * Ifges(5,4) + Ifges(5,5) * t228;
t55 = t122 * Ifges(5,4) + t121 * Ifges(5,2) + Ifges(5,6) * t228;
t47 = pkin(7) * t213 + t65;
t46 = t125 * t226 + t163 * t281;
t44 = qJD(3) * t288 - t109 * t163;
t39 = mrSges(6,1) * t145 - mrSges(6,3) * t63;
t38 = -mrSges(6,2) * t145 + mrSges(6,3) * t199;
t34 = t81 + (t163 * t204 + t218) * qJD(3);
t33 = -pkin(4) * t86 + t72;
t29 = t87 * Ifges(5,4) + t86 * Ifges(5,2) + t132 * Ifges(5,6);
t28 = -mrSges(6,1) * t199 + mrSges(6,2) * t63;
t13 = -mrSges(6,2) * t123 + mrSges(6,3) * t20;
t12 = mrSges(6,1) * t123 - mrSges(6,3) * t19;
t7 = -qJD(5) * t27 - t159 * t47 + t162 * t34;
t6 = qJD(5) * t26 + t159 * t34 + t162 * t47;
t3 = [(-Ifges(6,5) * t99 - Ifges(6,6) * t97) * t264 + t33 * (mrSges(6,1) * t97 - mrSges(6,2) * t99) + (-t92 * mrSges(6,1) + t91 * mrSges(6,2) + (m(3) * pkin(1) - t282 * m(6) + (-m(4) - m(5)) * t165 + t301) * t161 + ((-m(3) + t303) * qJ(2) + t300) * t164) * g(1) + (t188 * t308 - t302) * qJD(3) + (-m(3) * t229 - t90 * mrSges(6,1) - t89 * mrSges(6,2) + t303 * (t164 * pkin(6) + t229) + (-m(6) * t258 - t301) * t164 + t300 * t161) * g(2) - t183 * mrSges(4,3) - (t157 * t56 + t114) * t226 / 0.2e1 + (-t234 * t24 - t236 * t25) * mrSges(5,3) + (Ifges(6,1) * t44 + Ifges(6,4) * t46) * t269 + (-t1 * t97 + t2 * t99 - t44 * t8 + t46 * t9) * mrSges(6,3) + t147 * mrSges(3,2) + qJD(2) * t129 + qJ(2) * (mrSges(4,1) * t132 + mrSges(4,2) * t131) - t120 * t5 + t105 * t28 + t65 * t82 + t64 * t83 + t67 * (-mrSges(6,1) * t46 + mrSges(6,2) * t44) + t77 * t53 + t78 * t52 + t6 * t38 + t7 * t39 + t313 * t221 + t138 * t211 + t26 * t12 + t27 * t13 + t180 * t202 + (-Ifges(4,4) * t131 / 0.2e1 + Ifges(6,6) * t276 + Ifges(6,5) * t277 + Ifges(6,3) * t264 + t306 * t224 - Ifges(4,6) * qJDD(3) + Ifges(5,5) * t267 + Ifges(5,6) * t268 + t24 * mrSges(5,1) - t25 * mrSges(5,2) + t320 + (Ifges(4,2) + Ifges(5,3)) * t263) * t160 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-Ifges(6,1) * t99 - Ifges(6,4) * t97) * t277 + (Ifges(6,5) * t269 + Ifges(6,6) * t271 + Ifges(6,3) * t261 - t310) * t225 + (Ifges(6,4) * t44 + Ifges(6,2) * t46) * t271 + (-Ifges(6,4) * t99 - Ifges(6,2) * t97) * t276 + t234 * t273 + t44 * t274 + t46 * t275 - t99 * t278 - t97 * t279 + t193 * t241 + t102 * t231 + (Ifges(4,5) * qJDD(3) + (t101 - t35) * t165 + t187 * t263 + t191 * t267 + t189 * t268) * t163 + (Ifges(6,5) * t44 + Ifges(6,6) * t46) * t261 + m(6) * (t1 * t27 + t105 * t67 - t120 * t33 + t2 * t26 + t6 * t9 + t7 * t8) + (Ifges(4,1) * t131 - Ifges(4,4) * t132) * t316 + (Ifges(5,5) * t87 + Ifges(5,6) * t86 + Ifges(5,3) * t132 + t217) * t317 + t305 * t308 - t132 * t190 / 0.2e1 + t131 * t192 / 0.2e1 + (t195 + 0.2e1 * mrSges(3,3)) * t144 + m(4) * (t165 * t183 + t200) + m(3) * (-pkin(1) * t147 + t200) + t55 * t213 / 0.2e1 - pkin(1) * t223 - t29 * t236 / 0.2e1 + m(5) * (t24 * t77 + t25 * t78 + t50 * t64 + t51 * t65 + (t104 * t226 - t241) * t165); t223 - t96 * t12 - t288 * t13 + t294 * t39 + t293 * t38 + (qJ(2) * t266 - mrSges(3,3)) * t280 + (t101 + t297) * t163 + (t102 + t284) * t160 + (-t156 * t82 - t157 * t83 - t129) * qJD(1) + ((t138 + t285) * t163 + (t28 + t306) * t160) * qJD(3) + m(4) * t183 + m(3) * t147 - t309 * (-t265 - t266) + (-t1 * t288 - t163 * t33 - t2 * t96 + t67 * t226 + t293 * t9 + t294 * t8) * m(6) + (-t241 + t185 * t160 + (t163 * t184 + t239) * qJD(3) - (t156 * t51 + t157 * t50) * qJD(1)) * m(5); (t310 + t318) * t227 + (t273 + t299) * t156 + (-Ifges(6,4) * t281 - Ifges(6,2) * t109) * t271 + (-Ifges(6,5) * t281 - Ifges(6,6) * t109) * t261 + (-Ifges(6,1) * t281 - Ifges(6,4) * t109) * t269 - t281 * t274 + t33 * (mrSges(6,1) * t182 + mrSges(6,2) * t125) + (Ifges(6,4) * t125 - Ifges(6,2) * t182) * t276 + (Ifges(6,1) * t125 - Ifges(6,4) * t182) * t277 + (Ifges(6,5) * t125 - Ifges(6,6) * t182) * t264 + (-t1 * t182 - t125 * t2 - t291 * t9 + t292 * t8) * mrSges(6,3) - t182 * t279 - t146 * t5 + Ifges(4,5) * t131 - Ifges(4,6) * t132 - t88 * t28 - t71 * t82 - t70 * t83 - t84 * t22 / 0.2e1 - t85 * t23 / 0.2e1 + t73 * t12 + t74 * t13 + t75 * mrSges(4,1) - t76 * mrSges(4,2) - pkin(3) * t35 + (mrSges(6,1) * t291 - mrSges(6,2) * t292) * t67 + (t180 / 0.2e1 - t313) * t280 + (Ifges(6,5) * t85 + Ifges(6,6) * t84) * t262 + t188 * t202 + t114 * t208 + (-t305 / 0.2e1 + t302) * qJD(1) + t309 * (-t160 * t311 - t163 * t287 - t196) + t295 * t39 + (t1 * t74 - t146 * t33 + t2 * t73 + t295 * t8 + t296 * t9 - t67 * t88) * m(6) + t296 * t38 + (Ifges(6,1) * t85 + Ifges(6,4) * t84) * t270 + t185 * mrSges(5,3) + t284 * qJ(4) + t285 * qJD(4) + (Ifges(6,4) * t85 + Ifges(6,2) * t84) * t272 - t109 * t275 + t125 * t278 + (-pkin(3) * t72 + qJ(4) * t185 + qJD(4) * t184 - t104 * t133 - t50 * t70 - t51 * t71) * m(5) + t250 * t268 + (t29 / 0.2e1 + t56 * t208 + Ifges(5,6) * t263 + Ifges(5,2) * t268) * t157 + t249 * t267 + (t160 * t287 + t319) * g(3) + t72 * t194 + Ifges(4,3) * qJDD(3) - t306 * t133 - t55 * t214 / 0.2e1 - t138 * t237; t265 * t160 * g(3) - t121 * t82 + t122 * t83 - t199 * t38 + t63 * t39 + (-t199 * t9 + t63 * t8 + t289 + t33) * m(6) + (-t121 * t51 + t122 * t50 + t289 + t72) * m(5) - t297; -t67 * (mrSges(6,1) * t63 + mrSges(6,2) * t199) + (Ifges(6,1) * t199 - t259) * t270 + t22 * t269 + (Ifges(6,5) * t199 - Ifges(6,6) * t63) * t262 - t8 * t38 + t9 * t39 - g(1) * (mrSges(6,1) * t89 - mrSges(6,2) * t90) - g(2) * (mrSges(6,1) * t91 + mrSges(6,2) * t92) - g(3) * (-mrSges(6,1) * t148 - mrSges(6,2) * t149) * t163 + (t199 * t8 + t63 * t9) * mrSges(6,3) + t217 + (-Ifges(6,2) * t63 + t23 + t58) * t272 + t320;];
tau = t3;
