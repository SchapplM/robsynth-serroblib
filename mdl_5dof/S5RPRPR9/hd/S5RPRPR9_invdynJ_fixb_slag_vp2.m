% Calculate vector of inverse dynamics joint torques for
% S5RPRPR9
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:23:56
% DurationCPUTime: 7.36s
% Computational Cost: add. (2132->424), mult. (4271->569), div. (0->0), fcn. (2265->10), ass. (0->203)
t279 = Ifges(4,5) - Ifges(5,4);
t278 = Ifges(5,5) - Ifges(4,6);
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t202 = qJD(3) * t116;
t120 = cos(qJ(3));
t203 = qJD(1) * t120;
t72 = -t119 * t203 - t202;
t117 = sin(qJ(3));
t193 = qJD(1) * qJD(3);
t76 = -t120 * qJDD(1) + t117 * t193;
t29 = qJD(5) * t72 + qJDD(3) * t119 + t116 * t76;
t77 = qJDD(1) * t117 + t120 * t193;
t69 = qJDD(5) + t77;
t13 = mrSges(6,1) * t69 - mrSges(6,3) * t29;
t200 = qJD(3) * t119;
t134 = t116 * t203 - t200;
t30 = qJD(5) * t134 - qJDD(3) * t116 + t119 * t76;
t14 = -mrSges(6,2) * t69 + mrSges(6,3) * t30;
t204 = qJD(1) * t117;
t92 = qJD(5) + t204;
t37 = -mrSges(6,2) * t92 + mrSges(6,3) * t72;
t38 = mrSges(6,1) * t92 + mrSges(6,3) * t134;
t264 = -t116 * t38 + t119 * t37;
t288 = -qJD(5) * t264 - t116 * t14 - t119 * t13;
t199 = qJD(3) * t120;
t114 = sin(pkin(8));
t96 = pkin(1) * t114 + pkin(6);
t283 = qJD(2) * qJD(3) + t96 * qJDD(1);
t80 = t96 * qJD(1);
t23 = qJDD(2) * t120 - t283 * t117 - t80 * t199;
t133 = qJDD(4) - t23;
t249 = pkin(3) + pkin(7);
t10 = pkin(4) * t77 - qJDD(3) * t249 + t133;
t198 = qJD(4) * t117;
t115 = cos(pkin(8));
t245 = pkin(1) * t115;
t97 = -pkin(2) - t245;
t79 = t97 * qJDD(1);
t125 = -qJ(4) * t77 - qJD(1) * t198 + t79;
t12 = t249 * t76 + t125;
t108 = t120 * qJD(2);
t35 = -t117 * (pkin(4) * qJD(1) + t80) + t108;
t272 = qJD(4) - t35;
t31 = -qJD(3) * t249 + t272;
t109 = t117 * qJ(4);
t171 = -pkin(2) - t109;
t34 = (-t120 * t249 + t171 - t245) * qJD(1);
t8 = -t116 * t34 + t119 * t31;
t1 = qJD(5) * t8 + t10 * t116 + t119 * t12;
t287 = t1 * mrSges(6,2);
t9 = t116 * t31 + t119 * t34;
t2 = -qJD(5) * t9 + t10 * t119 - t116 * t12;
t286 = t2 * mrSges(6,1);
t112 = qJ(1) + pkin(8);
t104 = sin(t112);
t241 = g(2) * t104;
t250 = m(5) + m(6);
t285 = -m(4) - t250;
t18 = -qJDD(3) * pkin(3) + t133;
t49 = qJD(2) * t117 + t120 * t80;
t41 = -qJD(3) * qJ(4) - t49;
t284 = -qJD(3) * t41 - t18;
t156 = t120 * mrSges(5,2) - t117 * mrSges(5,3);
t160 = mrSges(4,1) * t120 - mrSges(4,2) * t117;
t282 = t156 - t160;
t105 = cos(t112);
t263 = g(1) * t105 + t241;
t162 = t1 * t116 + t119 * t2;
t196 = qJD(5) * t119;
t197 = qJD(5) * t116;
t281 = -t9 * t196 + t8 * t197 - t162;
t280 = g(3) * t120;
t53 = mrSges(5,1) * t76 - qJDD(3) * mrSges(5,3);
t277 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t76 - t53;
t100 = Ifges(4,4) * t203;
t276 = Ifges(4,1) * t204 + Ifges(4,5) * qJD(3) - Ifges(6,5) * t134 + t72 * Ifges(6,6) + t92 * Ifges(6,3) + t100;
t32 = -mrSges(6,1) * t72 - mrSges(6,2) * t134;
t183 = mrSges(5,1) * t203;
t84 = -qJD(3) * mrSges(5,3) - t183;
t275 = -t84 + t32;
t181 = mrSges(4,3) * t203;
t83 = -qJD(3) * mrSges(4,2) + t181;
t274 = t84 - t83;
t182 = mrSges(4,3) * t204;
t184 = mrSges(5,1) * t204;
t273 = t182 + t184 + (-mrSges(4,1) + mrSges(5,2)) * qJD(3);
t230 = Ifges(5,6) * t120;
t231 = Ifges(5,6) * t117;
t270 = t117 * (-Ifges(5,2) * t120 + t231) + t120 * (Ifges(5,3) * t117 - t230);
t269 = t278 * t117 + t279 * t120;
t188 = t117 * qJDD(2) + t283 * t120;
t201 = qJD(3) * t117;
t22 = -t201 * t80 + t188;
t266 = -t117 * t23 + t120 * t22;
t224 = t117 * t80;
t17 = -qJDD(3) * qJ(4) + qJD(3) * (-qJD(4) + t224) - t188;
t265 = t117 * t18 - t120 * t17;
t262 = -mrSges(3,1) + t282;
t145 = -t120 * Ifges(5,3) - t231;
t246 = Ifges(6,4) * t134;
t25 = Ifges(6,2) * t72 + Ifges(6,6) * t92 - t246;
t65 = Ifges(6,4) * t72;
t26 = -Ifges(6,1) * t134 + Ifges(6,5) * t92 + t65;
t261 = Ifges(5,5) * qJD(3) + qJD(1) * t145 + t116 * t26 + t119 * t25;
t210 = t105 * t120;
t213 = qJ(4) * t120;
t260 = -g(1) * qJ(4) * t210 - t213 * t241;
t54 = t77 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t259 = t288 - t54;
t258 = -m(6) * pkin(4) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t256 = -t29 * Ifges(6,4) / 0.2e1 - t30 * Ifges(6,2) / 0.2e1 - t69 * Ifges(6,6) / 0.2e1;
t255 = t29 / 0.2e1;
t254 = t30 / 0.2e1;
t253 = t69 / 0.2e1;
t251 = -t134 / 0.2e1;
t247 = pkin(4) + t96;
t118 = sin(qJ(1));
t244 = pkin(1) * t118;
t243 = pkin(7) * t120;
t110 = t120 * pkin(3);
t121 = cos(qJ(1));
t111 = t121 * pkin(1);
t235 = Ifges(4,4) * t117;
t234 = Ifges(4,4) * t120;
t233 = Ifges(6,4) * t116;
t232 = Ifges(6,4) * t119;
t219 = t120 * mrSges(5,3);
t102 = pkin(4) * t203;
t33 = t102 - t41;
t212 = qJD(3) * t33;
t209 = t116 * t117;
t208 = t116 * t120;
t207 = t117 * t119;
t206 = t119 * t120;
t205 = t110 + t109;
t195 = qJD(5) * t120;
t189 = Ifges(6,5) * t29 + Ifges(6,6) * t30 + Ifges(6,3) * t69;
t187 = t105 * pkin(2) + t104 * pkin(6) + t111;
t67 = t247 * t120;
t178 = t116 * t195;
t172 = -t196 / 0.2e1;
t169 = -t193 / 0.2e1;
t48 = -t108 + t224;
t166 = pkin(3) * t201 - t198;
t165 = -m(6) * t249 - mrSges(6,3);
t164 = pkin(3) * t210 + t105 * t109 + t187;
t161 = t116 * t8 - t119 * t9;
t159 = mrSges(4,1) * t117 + mrSges(4,2) * t120;
t158 = mrSges(6,1) * t119 - mrSges(6,2) * t116;
t157 = mrSges(6,1) * t116 + mrSges(6,2) * t119;
t155 = Ifges(6,1) * t119 - t233;
t154 = Ifges(6,1) * t116 + t232;
t153 = t120 * Ifges(4,2) + t235;
t151 = -Ifges(6,2) * t116 + t232;
t150 = Ifges(6,2) * t119 + t233;
t148 = Ifges(6,5) * t119 - Ifges(6,6) * t116;
t147 = Ifges(6,5) * t116 + Ifges(6,6) * t119;
t146 = -t117 * Ifges(5,2) - t230;
t144 = pkin(7) * t117 - t213;
t63 = t97 - t205;
t47 = t63 - t243;
t66 = t247 * t117;
t21 = t116 * t66 + t119 * t47;
t20 = -t116 * t47 + t119 * t66;
t141 = t171 - t110;
t42 = (t141 - t245) * qJD(1);
t139 = t42 * (-mrSges(5,2) * t117 - t219);
t138 = t97 * qJD(1) * t159;
t137 = t117 * (Ifges(4,1) * t120 - t235);
t131 = t117 * t200 + t178;
t130 = t116 * t201 - t119 * t195;
t129 = Ifges(6,5) * t120 + t117 * t154;
t128 = Ifges(6,6) * t120 + t117 * t150;
t127 = Ifges(6,3) * t120 + t117 * t147;
t126 = -qJD(5) * t161 + t162;
t101 = pkin(3) * t204;
t75 = -qJ(4) * t203 + t101;
t74 = t156 * qJD(1);
t64 = t158 * t120;
t61 = Ifges(5,4) * qJD(3) + qJD(1) * t146;
t58 = Ifges(4,6) * qJD(3) + qJD(1) * t153;
t57 = -qJ(4) * t199 + t166;
t56 = qJD(3) * t67;
t55 = t247 * t201;
t52 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t77;
t50 = qJD(1) * t144 + t101;
t46 = -t104 * t209 + t105 * t119;
t45 = t104 * t207 + t105 * t116;
t44 = t104 * t119 + t105 * t209;
t43 = -t104 * t116 + t105 * t207;
t40 = qJD(3) * t144 + t166;
t39 = -qJD(3) * pkin(3) + qJD(4) + t48;
t36 = t102 + t49;
t19 = pkin(3) * t76 + t125;
t16 = t116 * t36 + t119 * t50;
t15 = -t116 * t50 + t119 * t36;
t11 = -pkin(4) * t76 - t17;
t7 = -mrSges(6,1) * t30 + mrSges(6,2) * t29;
t6 = -qJD(5) * t21 - t116 * t40 + t119 * t56;
t5 = qJD(5) * t20 + t116 * t56 + t119 * t40;
t4 = t29 * Ifges(6,1) + t30 * Ifges(6,4) + t69 * Ifges(6,5);
t3 = [(m(5) * ((t117 * t41 + t120 * t39) * qJD(3) + t265) + m(4) * ((-t117 * t49 + t120 * t48) * qJD(3) + t266) + (-t52 + t54) * t117 + t273 * t199 + t274 * t201 + t277 * t120) * t96 + t266 * mrSges(4,3) + t265 * mrSges(5,1) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t115 - 0.2e1 * mrSges(3,2) * t114 + m(3) * (t114 ^ 2 + t115 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (m(4) * t97 - t160) * t79 - t117 * t287 + (m(3) * t244 + mrSges(2,1) * t118 - t46 * mrSges(6,1) + mrSges(2,2) * t121 + t45 * mrSges(6,2) + t285 * (t105 * pkin(6) - t244) + t258 * t105 + (m(4) * pkin(2) - m(5) * t141 - m(6) * t171 - t120 * t165 - t262) * t104) * g(1) + (-Ifges(4,4) * t76 + Ifges(4,5) * qJDD(3) + t189) * t117 / 0.2e1 + (-m(3) * t111 - m(6) * (pkin(7) * t210 + t164) - t44 * mrSges(6,1) - t43 * mrSges(6,2) - mrSges(2,1) * t121 + mrSges(2,2) * t118 - m(4) * t187 - m(5) * t164 + t262 * t105 + t258 * t104) * g(2) + (-g(2) * t210 - t1 * t206 - t130 * t8 + t131 * t9 + t2 * t208) * mrSges(6,3) + t269 * qJD(3) ^ 2 / 0.2e1 + t270 * t169 + (t120 * (-Ifges(4,2) * t117 + t234) + t137) * t193 / 0.2e1 + t120 * t26 * t172 + t77 * t234 / 0.2e1 - t4 * t208 / 0.2e1 + t72 * (qJD(3) * t128 - t151 * t195) / 0.2e1 + t92 * (qJD(3) * t127 - t148 * t195) / 0.2e1 + t120 * (Ifges(4,4) * t77 - Ifges(4,2) * t76 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t120 * (Ifges(5,5) * qJDD(3) - Ifges(5,6) * t77 + Ifges(5,3) * t76) / 0.2e1 - t117 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t77 + Ifges(5,6) * t76) / 0.2e1 + t97 * (mrSges(4,1) * t76 + mrSges(4,2) * t77) + t63 * (-mrSges(5,2) * t76 - mrSges(5,3) * t77) + t57 * t74 + t77 * t117 * Ifges(4,1) + m(5) * (t19 * t63 + t42 * t57) - t55 * t32 + t11 * t64 + t67 * t7 + t25 * t178 / 0.2e1 + t5 * t37 + t6 * t38 + t20 * t13 + t21 * t14 + m(6) * (t1 * t21 + t11 * t67 + t2 * t20 - t33 * t55 + t5 * t9 + t6 * t8) + (t39 * mrSges(5,1) + t48 * mrSges(4,3) - t61 / 0.2e1 - t9 * mrSges(6,2) + t8 * mrSges(6,1) + t276 / 0.2e1) * t199 + t117 * t286 + (qJD(3) * t129 - t155 * t195) * t251 + (Ifges(6,3) * t117 - t120 * t147) * t253 + (Ifges(6,6) * t117 - t120 * t150) * t254 + (Ifges(6,5) * t117 - t120 * t154) * t255 + t206 * t256 + (t41 * mrSges(5,1) - t49 * mrSges(4,3) + t261 / 0.2e1 - t58 / 0.2e1) * t201 + (t279 * t117 - t278 * t120) * qJDD(3) / 0.2e1 + t33 * (-mrSges(6,1) * t131 + mrSges(6,2) * t130) + qJD(3) * t138 + qJD(3) * t139 + t76 * t145 / 0.2e1 - t77 * t146 / 0.2e1 - t76 * t153 / 0.2e1 + t19 * t156; m(3) * qJDD(2) + (-m(3) + t285) * g(3) + (t7 + (t116 * t37 + t119 * t38 + t273) * qJD(3) + m(4) * (qJD(3) * t48 + t22) + m(5) * (qJD(3) * t39 - t17) + m(6) * (t200 * t8 + t202 * t9 + t11) + t277) * t117 + (t52 + (t83 + t275) * qJD(3) + m(4) * (qJD(3) * t49 + t23) + m(5) * t284 + m(6) * (t212 + t281) + t259) * t120; (t270 / 0.2e1 - t137 / 0.2e1) * qJD(1) ^ 2 - (-t134 * t154 + t147 * t92 + t150 * t72) * qJD(5) / 0.2e1 - (t92 * t127 + t72 * t128 - t129 * t134) * qJD(1) / 0.2e1 + (-t280 + t281) * mrSges(6,3) + (-m(5) * t205 - m(6) * (t205 + t243) - t157 * t117 + t282) * g(3) + (-t120 * t157 + t159 - t219 + (m(5) * pkin(3) - mrSges(5,2) - t165) * t117) * t263 + t269 * t169 - t261 * t204 / 0.2e1 + (-t53 + t7) * qJ(4) + t92 * t33 * t158 + (-m(6) * t126 + t288) * t249 + (-pkin(3) * t18 - qJ(4) * t17 - qJD(4) * t41 - t42 * t75 + t260) * m(5) + (qJ(4) * t11 - t15 * t8 - t16 * t9 + t272 * t33 + t260) * m(6) + t61 * t203 / 0.2e1 + t58 * t204 / 0.2e1 - t26 * t197 / 0.2e1 + t119 * t4 / 0.2e1 - t39 * t183 - t41 * t184 - t75 * t74 - pkin(3) * t54 - t16 * t37 - t15 * t38 - t35 * t32 + t18 * mrSges(5,2) - t22 * mrSges(4,2) + t23 * mrSges(4,1) - t17 * mrSges(5,3) + t148 * t253 + t151 * t254 + t155 * t255 + t116 * t256 + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + (-t9 * (-mrSges(6,2) * t120 + mrSges(6,3) * t207) - t8 * (mrSges(6,1) * t120 - mrSges(6,3) * t209) - t138 - t139) * qJD(1) + (-m(5) * t39 + t182 - t273) * t49 + (-m(5) * t41 - t181 - t274) * t48 + t275 * qJD(4) - (-Ifges(4,2) * t204 + t100 + t276) * t203 / 0.2e1 + t278 * t76 + t279 * t77 + t25 * t172 + t11 * t157; -t275 * qJD(3) + t250 * t280 + ((t264 + t74) * qJD(1) - t250 * t263) * t117 - t259 + (-t161 * t204 + t126 - t212) * m(6) + (t204 * t42 - t284) * m(5); -t287 + t286 - t33 * (-mrSges(6,1) * t134 + mrSges(6,2) * t72) + t134 * (Ifges(6,1) * t72 + t246) / 0.2e1 + t25 * t251 - t92 * (Ifges(6,5) * t72 + Ifges(6,6) * t134) / 0.2e1 - t8 * t37 + t9 * t38 - g(1) * (mrSges(6,1) * t43 - mrSges(6,2) * t44) - g(2) * (mrSges(6,1) * t45 + mrSges(6,2) * t46) + g(3) * t64 + (-t134 * t9 + t72 * t8) * mrSges(6,3) + t189 - (Ifges(6,2) * t134 + t26 + t65) * t72 / 0.2e1;];
tau = t3;
