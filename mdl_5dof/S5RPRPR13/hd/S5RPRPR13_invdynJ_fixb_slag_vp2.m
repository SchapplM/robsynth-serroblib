% Calculate vector of inverse dynamics joint torques for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:57
% EndTime: 2019-12-31 18:32:15
% DurationCPUTime: 9.86s
% Computational Cost: add. (3406->458), mult. (8085->568), div. (0->0), fcn. (5613->10), ass. (0->215)
t283 = -mrSges(5,1) - mrSges(4,3);
t282 = mrSges(4,2) - mrSges(5,3);
t281 = -mrSges(5,2) + mrSges(4,1);
t133 = sin(pkin(8));
t134 = cos(pkin(8));
t199 = t133 ^ 2 + t134 ^ 2;
t194 = qJD(1) * qJD(2);
t117 = qJDD(1) * qJ(2) + t194;
t265 = Ifges(5,5) - Ifges(4,6);
t138 = sin(qJ(1));
t234 = g(2) * t138;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t238 = cos(qJ(3));
t185 = t238 * t134;
t171 = qJD(1) * t185;
t137 = sin(qJ(3));
t205 = t133 * t137;
t94 = qJD(1) * t205 - t171;
t72 = -qJD(3) * t136 + t139 * t94;
t280 = t72 * Ifges(6,6);
t101 = t133 * t238 + t137 * t134;
t95 = t101 * qJD(1);
t88 = qJD(5) + t95;
t279 = t88 * Ifges(6,3);
t266 = Ifges(5,4) - Ifges(4,5);
t73 = qJD(3) * t139 + t136 * t94;
t84 = Ifges(4,4) * t94;
t278 = t95 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t73 * Ifges(6,5) + t279 + t280 - t84;
t74 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t94;
t76 = mrSges(5,1) * t94 - qJD(3) * mrSges(5,3);
t277 = t74 - t76;
t264 = t281 * qJD(3) + t283 * t95;
t211 = qJDD(3) / 0.2e1;
t131 = pkin(8) + qJ(3);
t127 = sin(t131);
t128 = cos(t131);
t276 = t282 * t127 - t281 * t128;
t140 = cos(qJ(1));
t259 = g(1) * t140 + t234;
t221 = pkin(6) + qJ(2);
t106 = t221 * t133;
t102 = qJD(1) * t106;
t107 = t221 * t134;
t103 = qJD(1) * t107;
t207 = t103 * t137;
t67 = t238 * t102 + t207;
t152 = pkin(4) * t95 + t67;
t275 = t152 + qJD(4);
t274 = 0.2e1 * t211;
t166 = -mrSges(3,1) * t134 + mrSges(3,2) * t133;
t272 = -m(3) * pkin(1) - mrSges(2,1) + t166 + t276;
t179 = qJD(3) * t238;
t198 = qJD(3) * t137;
t176 = pkin(6) * qJDD(1) + t117;
t81 = t176 * t133;
t82 = t176 * t134;
t23 = t102 * t198 - t103 * t179 - t137 * t82 - t238 * t81;
t147 = qJDD(4) - t23;
t241 = pkin(3) + pkin(7);
t177 = qJDD(1) * t238;
t182 = t133 * t198;
t190 = qJDD(1) * t137;
t65 = qJD(1) * t182 - qJD(3) * t171 - t133 * t177 - t134 * t190;
t11 = -t65 * pkin(4) - qJDD(3) * t241 + t147;
t125 = pkin(2) * t134 + pkin(1);
t104 = -qJDD(1) * t125 + qJDD(2);
t144 = qJ(4) * t65 - qJD(4) * t95 + t104;
t97 = t101 * qJD(3);
t66 = qJD(1) * t97 + t133 * t190 - t134 * t177;
t8 = t241 * t66 + t144;
t105 = -qJD(1) * t125 + qJD(2);
t146 = -qJ(4) * t95 + t105;
t32 = t241 * t94 + t146;
t35 = -qJD(3) * t241 + t275;
t9 = -t136 * t32 + t139 * t35;
t1 = qJD(5) * t9 + t11 * t136 + t139 * t8;
t10 = t136 * t35 + t139 * t32;
t2 = -qJD(5) * t10 + t11 * t139 - t136 * t8;
t271 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t63 = qJDD(5) - t65;
t252 = t63 / 0.2e1;
t29 = qJD(5) * t72 + qJDD(3) * t139 + t136 * t66;
t254 = t29 / 0.2e1;
t270 = Ifges(6,1) * t254 + Ifges(6,5) * t252;
t180 = m(3) * qJ(2) + mrSges(3,3);
t269 = -m(6) * (pkin(4) + t221) + mrSges(2,2) - t180 + t283;
t43 = pkin(3) * t94 + t146;
t68 = -t137 * t102 + t103 * t238;
t58 = -qJD(3) * qJ(4) - t68;
t268 = t105 * mrSges(4,1) + t58 * mrSges(5,1) - t43 * mrSges(5,2) - t68 * mrSges(4,3);
t258 = -t67 - qJD(4);
t57 = -qJD(3) * pkin(3) - t258;
t267 = -t57 * mrSges(5,1) - t9 * mrSges(6,1) - t105 * mrSges(4,2) + t10 * mrSges(6,2) - t67 * mrSges(4,3) + t43 * mrSges(5,3);
t30 = -qJD(5) * t73 - qJDD(3) * t136 + t139 * t66;
t253 = t30 / 0.2e1;
t233 = g(3) * t128;
t38 = -mrSges(6,1) * t72 + mrSges(6,2) * t73;
t263 = -t76 + t38;
t262 = t127 * t259;
t121 = t127 * qJ(4);
t206 = t128 * t140;
t261 = pkin(3) * t206 + t140 * t121;
t163 = mrSges(6,1) * t139 - mrSges(6,2) * t136;
t236 = Ifges(6,4) * t73;
t25 = t72 * Ifges(6,2) + t88 * Ifges(6,6) + t236;
t240 = t94 * pkin(4);
t36 = -t58 - t240;
t257 = -t139 * t25 / 0.2e1 + t36 * t163;
t46 = t137 * (qJD(2) * t133 + qJD(3) * t107) - qJD(2) * t185 + t106 * t179;
t188 = -t102 * t179 - t137 * t81 + t238 * t82;
t20 = -qJDD(3) * qJ(4) + qJD(3) * (-qJD(4) + t207) - t188;
t256 = (-g(1) * t206 - t128 * t234) * qJ(4);
t255 = Ifges(6,4) * t253 + t270;
t251 = -t72 / 0.2e1;
t250 = -t73 / 0.2e1;
t249 = t73 / 0.2e1;
t248 = -t88 / 0.2e1;
t247 = -t94 / 0.2e1;
t246 = t94 / 0.2e1;
t245 = -t95 / 0.2e1;
t244 = t95 / 0.2e1;
t242 = m(5) + m(6);
t232 = t1 * t136;
t122 = t128 * pkin(3);
t227 = t95 * Ifges(4,4);
t226 = t95 * Ifges(5,6);
t225 = -qJD(3) / 0.2e1;
t224 = qJD(3) / 0.2e1;
t219 = mrSges(6,3) * t139;
t218 = Ifges(6,4) * t136;
t217 = Ifges(6,4) * t139;
t216 = qJ(4) * t94;
t214 = t136 * t95;
t213 = t136 * t97;
t212 = t139 * t97;
t100 = -t185 + t205;
t209 = t100 * t136;
t208 = t100 * t139;
t204 = t136 * t138;
t203 = t136 * t140;
t202 = t138 * t139;
t201 = t139 * t140;
t200 = t122 + t121;
t197 = qJD(5) * t136;
t196 = qJD(5) * t139;
t195 = qJD(5) * t241;
t192 = qJDD(1) * t133;
t191 = qJDD(1) * t134;
t189 = Ifges(6,5) * t29 + Ifges(6,6) * t30 + Ifges(6,3) * t63;
t178 = -t197 / 0.2e1;
t51 = -t65 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t69 = t238 * t106 + t107 * t137;
t175 = -t125 - t121;
t111 = t140 * t125;
t174 = t138 * t221 + t111;
t172 = -m(6) * t241 - mrSges(6,3);
t169 = -mrSges(3,1) * t191 + mrSges(3,2) * t192;
t168 = t10 * t139 - t136 * t9;
t162 = mrSges(6,1) * t136 + mrSges(6,2) * t139;
t160 = Ifges(6,1) * t136 + t217;
t159 = Ifges(6,2) * t139 + t218;
t158 = Ifges(6,5) * t136 + Ifges(6,6) * t139;
t153 = -qJ(4) * t101 - t125;
t40 = t100 * t241 + t153;
t48 = pkin(4) * t101 + t69;
t19 = t136 * t48 + t139 * t40;
t18 = -t136 * t40 + t139 * t48;
t44 = -mrSges(6,2) * t88 + mrSges(6,3) * t72;
t45 = mrSges(6,1) * t88 - mrSges(6,3) * t73;
t157 = -t136 * t44 - t139 * t45;
t96 = -t134 * t179 + t182;
t156 = qJ(4) * t96 - qJD(4) * t101;
t151 = t100 * t196 + t213;
t150 = t100 * t197 - t212;
t70 = -t137 * t106 + t107 * t238;
t145 = qJD(5) * t168 + t139 * t2 + t232;
t47 = qJD(2) * t101 + qJD(3) * t70;
t126 = -qJDD(1) * pkin(1) + qJDD(2);
t92 = -t127 * t204 + t201;
t91 = t127 * t202 + t203;
t90 = t127 * t203 + t202;
t89 = t127 * t201 - t204;
t83 = Ifges(5,6) * t94;
t71 = Ifges(6,4) * t72;
t64 = pkin(3) * t100 + t153;
t62 = t65 * mrSges(4,2);
t61 = t65 * mrSges(5,3);
t60 = -mrSges(5,2) * t94 - mrSges(5,3) * t95;
t59 = pkin(3) * t95 + t216;
t54 = -t94 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t227;
t53 = Ifges(5,4) * qJD(3) - t95 * Ifges(5,2) + t83;
t52 = Ifges(5,5) * qJD(3) + t94 * Ifges(5,3) - t226;
t50 = mrSges(5,1) * t66 - qJDD(3) * mrSges(5,3);
t49 = -t100 * pkin(4) + t70;
t42 = t68 - t240;
t39 = pkin(3) * t97 + t156;
t37 = t241 * t95 + t216;
t34 = -t96 * pkin(4) + t47;
t33 = -pkin(4) * t97 - t46;
t31 = t241 * t97 + t156;
t26 = t73 * Ifges(6,1) + t88 * Ifges(6,5) + t71;
t22 = -t103 * t198 + t188;
t21 = -qJDD(3) * pkin(3) + t147;
t17 = pkin(3) * t66 + t144;
t16 = -mrSges(6,2) * t63 + mrSges(6,3) * t30;
t15 = mrSges(6,1) * t63 - mrSges(6,3) * t29;
t14 = t136 * t42 + t139 * t37;
t13 = -t136 * t37 + t139 * t42;
t12 = -pkin(4) * t66 - t20;
t7 = -mrSges(6,1) * t30 + mrSges(6,2) * t29;
t5 = t29 * Ifges(6,4) + t30 * Ifges(6,2) + t63 * Ifges(6,6);
t4 = -qJD(5) * t19 - t136 * t31 + t139 * t34;
t3 = qJD(5) * t18 + t136 * t34 + t139 * t31;
t6 = [(-m(4) * t104 - t66 * mrSges(4,1) + t62) * t125 + 0.2e1 * t199 * t117 * mrSges(3,3) + (t53 / 0.2e1 - t280 / 0.2e1 - t279 / 0.2e1 - Ifges(4,1) * t244 + Ifges(5,2) * t245 + Ifges(5,6) * t246 - Ifges(4,4) * t247 - Ifges(6,5) * t249 + t266 * t224 + t267 - t278 / 0.2e1) * t96 + t88 * (Ifges(6,5) * t151 - Ifges(6,6) * t150) / 0.2e1 + t72 * (Ifges(6,4) * t151 - Ifges(6,2) * t150) / 0.2e1 + (Ifges(6,1) * t151 - Ifges(6,4) * t150) * t249 + m(5) * (t17 * t64 + t39 * t43) + (-Ifges(4,4) * t66 + Ifges(4,5) * qJDD(3) + t189) * t101 / 0.2e1 + (t21 * mrSges(5,1) + t104 * mrSges(4,2) - t23 * mrSges(4,3) - t17 * mrSges(5,3) + Ifges(4,5) * t211 + Ifges(6,5) * t254 + Ifges(6,6) * t253 + Ifges(6,3) * t252 + (-Ifges(5,6) - Ifges(4,4) / 0.2e1) * t66 - Ifges(5,2) * t65 - t274 * Ifges(5,4) + t271) * t101 + (m(4) * t22 - m(5) * t20 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t66 - t50) * t70 + (-m(4) * t23 + m(5) * t21 - qJDD(3) * mrSges(4,1) - mrSges(4,3) * t65 + t51) * t69 + (qJD(5) * t26 + t5) * t208 / 0.2e1 + (t1 * t208 - t10 * t150 - t9 * t151 - t2 * t209) * mrSges(6,3) + t26 * t213 / 0.2e1 + t25 * t212 / 0.2e1 + t126 * t166 - pkin(1) * t169 + (m(4) * t67 + m(5) * t57 - t264) * t47 + (-m(4) * t68 + m(5) * t58 - t277) * t46 + (t52 / 0.2e1 - t54 / 0.2e1 - Ifges(4,4) * t244 + Ifges(5,6) * t245 + Ifges(5,3) * t246 - Ifges(4,2) * t247 + t265 * t224 + t268) * t97 + t36 * (mrSges(6,1) * t150 + mrSges(6,2) * t151) - t65 * Ifges(4,1) * t101 + (-t92 * mrSges(6,1) + t91 * mrSges(6,2) + ((-m(5) - m(4)) * t221 + t269) * t140 + (-m(5) * (t175 - t122) - m(6) * t175 - t128 * t172 + m(4) * t125 - t272) * t138) * g(1) + (-m(4) * t174 - m(6) * (pkin(7) * t206 + t111 + t261) - t90 * mrSges(6,1) - t89 * mrSges(6,2) - mrSges(6,3) * t206 - m(5) * (t174 + t261) + t272 * t140 + t269 * t138) * g(2) + t64 * (-t66 * mrSges(5,2) + t61) + t39 * t60 + t3 * t44 + t4 * t45 + t49 * t7 + t33 * t38 + t18 * t15 + t19 * t16 + (t104 * mrSges(4,1) + t20 * mrSges(5,1) - t17 * mrSges(5,2) - t22 * mrSges(4,3) - t12 * t163 + t158 * t252 + t159 * t253 + t160 * t254 + t25 * t178 + (Ifges(5,3) + Ifges(4,2)) * t66 + (Ifges(5,6) + Ifges(4,4)) * t65 + t265 * t274) * t100 + (Ifges(3,4) * t133 + Ifges(3,2) * t134) * t191 + (Ifges(3,1) * t133 + Ifges(3,4) * t134) * t192 + t209 * t255 + m(6) * (t1 * t19 + t10 * t3 + t12 * t49 + t18 * t2 + t33 * t36 + t4 * t9) + Ifges(2,3) * qJDD(1) + m(3) * (-pkin(1) * t126 + (t117 + t194) * qJ(2) * t199); -t136 * t15 + t139 * t16 + t61 - t62 + t281 * t66 + t157 * qJD(5) + (t74 + t263) * t94 + (t157 + t264) * t95 + m(3) * t126 + t169 + (-g(1) * t138 + g(2) * t140) * (m(3) + m(4) + t242) - t180 * t199 * qJD(1) ^ 2 + (t1 * t139 - t136 * t2 + t36 * t94 - t88 * (t10 * t136 + t139 * t9)) * m(6) + (-t57 * t95 - t58 * t94 + t17) * m(5) + (-t67 * t95 + t68 * t94 + t104) * m(4); (-t227 + t52) * t245 + (t226 + t54) * t244 - (t158 * t88 + t159 * t72 + t160 * t73) * qJD(5) / 0.2e1 + ((-t162 + t282) * t128 + (m(5) * pkin(3) - t172 + t281) * t127) * t259 + (t83 + t53) * t247 + (-t214 / 0.2e1 + t178) * t26 + t152 * t38 + (-t241 * t16 - t5 / 0.2e1 + t45 * t195 - Ifges(6,6) * t252 - Ifges(6,2) * t253) * t136 + t263 * qJD(4) + t264 * t68 + t265 * t66 + t266 * t65 + (-t233 - t10 * t196 - t232 + (t197 + t214) * t9) * mrSges(6,3) - t2 * t219 + t12 * t162 + t277 * t67 + (-t84 + t278) * t246 + (-m(5) * t200 - m(6) * (pkin(7) * t128 + t200) - t162 * t127 + t276) * g(3) + (qJ(4) * t12 - t10 * t14 - t13 * t9 - t241 * t145 + t275 * t36 + t256) * m(6) + (-t15 * t241 - t195 * t44 + t255 + t270) * t139 - t59 * t60 - t14 * t44 - t13 * t45 - pkin(3) * t51 + (-Ifges(4,2) * t246 + Ifges(5,3) * t247 - t10 * t219 + t158 * t248 + t159 * t251 + t160 * t250 + t225 * t265 + t257 - t268) * t95 + (-Ifges(4,1) * t245 - Ifges(6,5) * t250 + Ifges(5,2) * t244 - Ifges(6,6) * t251 - Ifges(6,3) * t248 + t225 * t266 - t267) * t94 - t20 * mrSges(5,3) + t21 * mrSges(5,2) - t22 * mrSges(4,2) + t23 * mrSges(4,1) + t257 * qJD(5) + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + (-pkin(3) * t21 - qJ(4) * t20 + t258 * t58 - t43 * t59 - t57 * t68 + t256) * m(5) + t217 * t253 + (-t50 + t7) * qJ(4) - t218 * t254; t95 * t60 - t263 * qJD(3) + t242 * t233 + (t44 * t88 + t15) * t139 + (-t45 * t88 + t16) * t136 + t51 + (-qJD(3) * t36 + t168 * t95 + t145 - t262) * m(6) + (qJD(3) * t58 + t43 * t95 + t21 - t262) * m(5); -t36 * (mrSges(6,1) * t73 + mrSges(6,2) * t72) + (Ifges(6,1) * t72 - t236) * t250 + t25 * t249 + (Ifges(6,5) * t72 - Ifges(6,6) * t73) * t248 - t9 * t44 + t10 * t45 - g(1) * (mrSges(6,1) * t89 - mrSges(6,2) * t90) - g(2) * (mrSges(6,1) * t91 + mrSges(6,2) * t92) + t163 * t233 + (t10 * t73 + t72 * t9) * mrSges(6,3) + t189 + (-Ifges(6,2) * t73 + t26 + t71) * t251 + t271;];
tau = t6;
