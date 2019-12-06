% Calculate vector of inverse dynamics joint torques for
% S5RPRRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:17
% EndTime: 2019-12-05 18:11:44
% DurationCPUTime: 11.64s
% Computational Cost: add. (7952->437), mult. (20248->589), div. (0->0), fcn. (15697->16), ass. (0->190)
t192 = sin(qJ(5));
t196 = cos(qJ(5));
t189 = sin(pkin(9));
t194 = sin(qJ(3));
t190 = cos(pkin(9));
t198 = cos(qJ(3));
t240 = t190 * t198;
t151 = -t189 * t194 + t240;
t142 = t151 * qJD(1);
t152 = t189 * t198 + t190 * t194;
t143 = t152 * qJD(1);
t193 = sin(qJ(4));
t197 = cos(qJ(4));
t114 = t142 * t193 + t143 * t197;
t214 = t197 * t142 - t143 * t193;
t288 = -t114 * t192 + t196 * t214;
t144 = t151 * qJD(3);
t118 = qJD(1) * t144 + qJDD(1) * t152;
t145 = t152 * qJD(3);
t119 = -qJD(1) * t145 + qJDD(1) * t151;
t52 = qJD(4) * t214 + t118 * t197 + t119 * t193;
t53 = -qJD(4) * t114 - t118 * t193 + t119 * t197;
t16 = qJD(5) * t288 + t192 * t53 + t196 * t52;
t72 = t114 * t196 + t192 * t214;
t17 = -qJD(5) * t72 - t192 * t52 + t196 * t53;
t183 = qJDD(3) + qJDD(4);
t177 = qJDD(5) + t183;
t282 = pkin(8) * t214;
t251 = pkin(6) + qJ(2);
t158 = t251 * t189;
t153 = qJD(1) * t158;
t159 = t251 * t190;
t154 = qJD(1) * t159;
t123 = -t153 * t194 + t154 * t198;
t95 = pkin(7) * t142 + t123;
t92 = t197 * t95;
t241 = t154 * t194;
t122 = -t198 * t153 - t241;
t94 = -pkin(7) * t143 + t122;
t93 = qJD(3) * pkin(3) + t94;
t56 = t193 * t93 + t92;
t38 = t56 + t282;
t245 = t192 * t38;
t188 = qJD(3) + qJD(4);
t291 = pkin(8) * t114;
t90 = t193 * t95;
t55 = t197 * t93 - t90;
t37 = t55 - t291;
t36 = pkin(4) * t188 + t37;
t18 = t196 * t36 - t245;
t180 = qJD(5) + t188;
t244 = t196 * t38;
t19 = t192 * t36 + t244;
t230 = qJD(1) * qJD(2);
t162 = qJ(2) * qJDD(1) + t230;
t216 = pkin(6) * qJDD(1) + t162;
t135 = t216 * t189;
t136 = t216 * t190;
t75 = -qJD(3) * t123 - t198 * t135 - t136 * t194;
t51 = qJDD(3) * pkin(3) - pkin(7) * t118 + t75;
t236 = qJD(3) * t198;
t74 = -qJD(3) * t241 - t194 * t135 + t198 * t136 - t153 * t236;
t54 = pkin(7) * t119 + t74;
t10 = -qJD(4) * t56 - t193 * t54 + t197 * t51;
t6 = pkin(4) * t183 - pkin(8) * t52 + t10;
t234 = qJD(4) * t197;
t235 = qJD(4) * t193;
t9 = t193 * t51 + t197 * t54 + t93 * t234 - t235 * t95;
t7 = pkin(8) * t53 + t9;
t2 = qJD(5) * t18 + t192 * t6 + t196 * t7;
t255 = Ifges(6,4) * t72;
t264 = t72 / 0.2e1;
t3 = -qJD(5) * t19 - t192 * t7 + t196 * t6;
t32 = Ifges(6,2) * t288 + Ifges(6,6) * t180 + t255;
t66 = Ifges(6,4) * t288;
t33 = Ifges(6,1) * t72 + Ifges(6,5) * t180 + t66;
t172 = pkin(2) * t190 + pkin(1);
t157 = -qJD(1) * t172 + qJD(2);
t124 = -pkin(3) * t142 + t157;
t80 = -pkin(4) * t214 + t124;
t304 = -(Ifges(6,1) * t288 - t255) * t72 / 0.2e1 - (Ifges(6,5) * t288 - Ifges(6,6) * t72) * t180 / 0.2e1 + (t18 * t288 + t19 * t72) * mrSges(6,3) - t80 * (mrSges(6,1) * t72 + mrSges(6,2) * t288) + t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t16 + Ifges(6,6) * t17 + Ifges(6,3) * t177 + t32 * t264 - (-Ifges(6,2) * t72 + t33 + t66) * t288 / 0.2e1;
t108 = Ifges(5,4) * t214;
t248 = Ifges(5,4) * t114;
t65 = Ifges(5,1) * t114 + Ifges(5,5) * t188 + t108;
t303 = t10 * mrSges(5,1) - t9 * mrSges(5,2) + Ifges(5,5) * t52 + Ifges(5,6) * t53 + Ifges(5,3) * t183 - (Ifges(5,5) * t214 - Ifges(5,6) * t114) * t188 / 0.2e1 + (t114 * t56 + t214 * t55) * mrSges(5,3) - (-Ifges(5,2) * t114 + t108 + t65) * t214 / 0.2e1 - t124 * (mrSges(5,1) * t114 + mrSges(5,2) * t214) - (Ifges(5,1) * t214 - t248) * t114 / 0.2e1 + t304;
t237 = t189 ^ 2 + t190 ^ 2;
t187 = pkin(9) + qJ(3);
t178 = sin(t187);
t179 = cos(t187);
t181 = qJ(4) + t187;
t170 = sin(t181);
t171 = cos(t181);
t173 = qJ(5) + t181;
t165 = sin(t173);
t166 = cos(t173);
t217 = t166 * mrSges(6,1) - t165 * mrSges(6,2);
t207 = -t171 * mrSges(5,1) + mrSges(5,2) * t170 - t217;
t294 = -mrSges(4,1) * t179 + mrSges(4,2) * t178 + t207;
t168 = pkin(3) * t179;
t155 = t168 + t172;
t164 = pkin(4) * t171;
t212 = -mrSges(3,1) * t190 + mrSges(3,2) * t189;
t287 = mrSges(2,1) + m(6) * (t155 + t164) + m(5) * t155 + m(4) * t172 + m(3) * pkin(1) - t212 - t294;
t184 = -pkin(7) - t251;
t224 = m(3) * qJ(2) + mrSges(3,3);
t286 = mrSges(2,2) + m(6) * (-pkin(8) + t184) - mrSges(6,3) + m(5) * t184 - mrSges(5,3) - m(4) * t251 - mrSges(4,3) - t224;
t64 = Ifges(5,2) * t214 + Ifges(5,6) * t188 + t248;
t283 = t64 / 0.2e1;
t174 = pkin(3) * t197 + pkin(4);
t232 = qJD(5) * t196;
t233 = qJD(5) * t192;
t238 = t193 * t196;
t57 = -t193 * t94 - t92;
t39 = t57 - t282;
t58 = t197 * t94 - t90;
t40 = t58 - t291;
t273 = t192 * t40 - t196 * t39 - t174 * t233 + (-t193 * t232 + (-t192 * t197 - t238) * qJD(4)) * pkin(3);
t239 = t192 * t193;
t272 = -t192 * t39 - t196 * t40 + t174 * t232 + (-t193 * t233 + (t196 * t197 - t239) * qJD(4)) * pkin(3);
t125 = -t198 * t158 - t159 * t194;
t106 = -pkin(7) * t152 + t125;
t126 = -t194 * t158 + t198 * t159;
t107 = pkin(7) * t151 + t126;
t63 = t193 * t106 + t197 * t107;
t195 = sin(qJ(1));
t199 = cos(qJ(1));
t269 = g(1) * t199 + g(2) * t195;
t261 = t114 / 0.2e1;
t259 = t143 / 0.2e1;
t254 = pkin(3) * t145;
t249 = Ifges(4,4) * t143;
t229 = qJDD(1) * t189;
t228 = qJDD(1) * t190;
t220 = -t53 * mrSges(5,1) + t52 * mrSges(5,2);
t219 = -t17 * mrSges(6,1) + t16 * mrSges(6,2);
t218 = -t119 * mrSges(4,1) + t118 * mrSges(4,2);
t62 = t197 * t106 - t107 * t193;
t213 = -mrSges(3,1) * t228 + mrSges(3,2) * t229;
t210 = mrSges(6,1) * t165 + mrSges(6,2) * t166;
t121 = t151 * t193 + t152 * t197;
t41 = -pkin(8) * t121 + t62;
t120 = t151 * t197 - t152 * t193;
t42 = pkin(8) * t120 + t63;
t26 = -t192 * t42 + t196 * t41;
t27 = t192 * t41 + t196 * t42;
t78 = t120 * t196 - t121 * t192;
t79 = t120 * t192 + t121 * t196;
t132 = -pkin(3) * t151 - t172;
t100 = -t158 * t236 + qJD(2) * t240 + (-qJD(2) * t189 - qJD(3) * t159) * t194;
t84 = -pkin(7) * t145 + t100;
t101 = -t152 * qJD(2) - qJD(3) * t126;
t85 = -pkin(7) * t144 + t101;
t30 = t106 * t234 - t107 * t235 + t193 * t85 + t197 * t84;
t156 = -qJDD(1) * t172 + qJDD(2);
t96 = -pkin(3) * t119 + t156;
t203 = mrSges(5,1) * t170 + mrSges(5,2) * t171 + t210;
t31 = -qJD(4) * t63 - t193 * t84 + t197 * t85;
t176 = -qJDD(1) * pkin(1) + qJDD(2);
t140 = pkin(3) * t238 + t174 * t192;
t139 = -pkin(3) * t239 + t174 * t196;
t137 = Ifges(4,4) * t142;
t128 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t143;
t127 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t142;
t110 = t143 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t137;
t109 = t142 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t249;
t98 = mrSges(5,1) * t188 - mrSges(5,3) * t114;
t97 = -mrSges(5,2) * t188 + mrSges(5,3) * t214;
t89 = -pkin(4) * t120 + t132;
t86 = pkin(3) * t143 + pkin(4) * t114;
t77 = -qJD(4) * t121 - t144 * t193 - t145 * t197;
t76 = qJD(4) * t120 + t144 * t197 - t145 * t193;
t73 = -mrSges(5,1) * t214 + mrSges(5,2) * t114;
t61 = mrSges(6,1) * t180 - mrSges(6,3) * t72;
t60 = -mrSges(6,2) * t180 + mrSges(6,3) * t288;
t59 = -pkin(4) * t77 + t254;
t44 = -mrSges(5,2) * t183 + mrSges(5,3) * t53;
t43 = mrSges(5,1) * t183 - mrSges(5,3) * t52;
t35 = -pkin(4) * t53 + t96;
t34 = -mrSges(6,1) * t288 + mrSges(6,2) * t72;
t29 = -qJD(5) * t79 - t192 * t76 + t196 * t77;
t28 = qJD(5) * t78 + t192 * t77 + t196 * t76;
t25 = -pkin(8) * t76 + t31;
t24 = pkin(8) * t77 + t30;
t21 = t196 * t37 - t245;
t20 = -t192 * t37 - t244;
t12 = -mrSges(6,2) * t177 + mrSges(6,3) * t17;
t11 = mrSges(6,1) * t177 - mrSges(6,3) * t16;
t5 = -qJD(5) * t27 - t192 * t24 + t196 * t25;
t4 = qJD(5) * t26 + t192 * t25 + t196 * t24;
t1 = [m(6) * (t18 * t5 + t19 * t4 + t2 * t27 + t26 * t3 + t35 * t89 + t59 * t80) + m(4) * (t100 * t123 + t101 * t122 + t125 * t75 + t126 * t74 - t156 * t172) + (-mrSges(4,1) * t156 + mrSges(4,3) * t74 + Ifges(4,4) * t118 + Ifges(4,2) * t119 + Ifges(4,6) * qJDD(3)) * t151 + t77 * t283 + (t125 * mrSges(4,1) - t126 * mrSges(4,2)) * qJDD(3) + (-t18 * t28 + t19 * t29) * mrSges(6,3) + (mrSges(4,2) * t156 - mrSges(4,3) * t75 + Ifges(4,1) * t118 + Ifges(4,4) * t119 + Ifges(4,5) * qJDD(3)) * t152 + (Ifges(4,1) * t144 - Ifges(4,4) * t145) * t259 + qJD(3) * (Ifges(4,5) * t144 - Ifges(4,6) * t145) / 0.2e1 + t142 * (Ifges(4,4) * t144 - Ifges(4,2) * t145) / 0.2e1 + t157 * (mrSges(4,1) * t145 + mrSges(4,2) * t144) + m(5) * (t10 * t62 + t124 * t254 + t132 * t96 + t30 * t56 + t31 * t55 + t63 * t9) + (-mrSges(5,1) * t96 + mrSges(5,3) * t9 + Ifges(5,4) * t52 + Ifges(5,2) * t53 + Ifges(5,6) * t183) * t120 + (-t55 * t76 + t56 * t77) * mrSges(5,3) + (Ifges(5,1) * t76 + Ifges(5,4) * t77) * t261 + (Ifges(6,1) * t28 + Ifges(6,4) * t29) * t264 + t29 * t32 / 0.2e1 + t28 * t33 / 0.2e1 + t26 * t11 + t27 * t12 + (t195 * t287 + t199 * t286) * g(1) + (t195 * t286 - t199 * t287) * g(2) + t214 * (Ifges(5,4) * t76 + Ifges(5,2) * t77) / 0.2e1 + (mrSges(6,2) * t35 - mrSges(6,3) * t3 + Ifges(6,1) * t16 + Ifges(6,4) * t17 + Ifges(6,5) * t177) * t79 + 0.2e1 * t237 * t162 * mrSges(3,3) + (mrSges(5,2) * t96 - mrSges(5,3) * t10 + Ifges(5,1) * t52 + Ifges(5,4) * t53 + Ifges(5,5) * t183) * t121 + t288 * (Ifges(6,4) * t28 + Ifges(6,2) * t29) / 0.2e1 + t176 * t212 - pkin(1) * t213 + t59 * t34 + t4 * t60 + t5 * t61 + t62 * t43 + t63 * t44 - t172 * t218 + t89 * t219 + t132 * t220 + t76 * t65 / 0.2e1 + t80 * (-mrSges(6,1) * t29 + mrSges(6,2) * t28) + t30 * t97 + t31 * t98 + t124 * (-mrSges(5,1) * t77 + mrSges(5,2) * t76) + t100 * t127 + t101 * t128 + t144 * t110 / 0.2e1 - t145 * t109 / 0.2e1 + t73 * t254 + t180 * (Ifges(6,5) * t28 + Ifges(6,6) * t29) / 0.2e1 + t188 * (Ifges(5,5) * t76 + Ifges(5,6) * t77) / 0.2e1 + m(3) * (-pkin(1) * t176 + (t162 + t230) * qJ(2) * t237) + (-mrSges(6,1) * t35 + mrSges(6,3) * t2 + Ifges(6,4) * t16 + Ifges(6,2) * t17 + Ifges(6,6) * t177) * t78 + (Ifges(3,4) * t189 + Ifges(3,2) * t190) * t228 + (Ifges(3,1) * t189 + Ifges(3,4) * t190) * t229 + Ifges(2,3) * qJDD(1) + (-t118 * t125 + t119 * t126 - t122 * t144 - t123 * t145) * mrSges(4,3); m(3) * t176 + t114 * t98 - t214 * t97 - t142 * t127 + t143 * t128 - t288 * t60 + t72 * t61 + t213 + t218 + t219 + t220 + (-g(1) * t195 + g(2) * t199) * (m(3) + m(4) + m(5) + m(6)) - t224 * t237 * qJD(1) ^ 2 + (t18 * t72 - t19 * t288 + t35) * m(6) + (t114 * t55 - t214 * t56 + t96) * m(5) + (t122 * t143 - t123 * t142 + t156) * m(4); t303 + t272 * t60 + (-t80 * t86 + t139 * t3 + t140 * t2 + (-t164 - t168) * g(3) + t272 * t19 + t273 * t18) * m(6) + t273 * t61 + t269 * (-m(6) * (-pkin(3) * t178 - pkin(4) * t170) + mrSges(4,2) * t179 + (m(5) * pkin(3) + mrSges(4,1)) * t178 + t203) + t109 * t259 - m(5) * (t55 * t57 + t56 * t58) + t114 * t283 - (-Ifges(4,2) * t143 + t110 + t137) * t142 / 0.2e1 + t294 * g(3) - t74 * mrSges(4,2) + t75 * mrSges(4,1) + (t122 * t142 + t123 * t143) * mrSges(4,3) - t86 * t34 - t58 * t97 - t57 * t98 + (-t143 * t73 + t193 * t44 + t197 * t43 + (-t193 * t98 + t197 * t97) * qJD(4) + (-g(3) * t179 + t10 * t197 - t124 * t143 + t193 * t9 + t234 * t56 - t235 * t55) * m(5)) * pkin(3) + Ifges(4,5) * t118 + Ifges(4,6) * t119 - t122 * t127 + t123 * t128 + t139 * t11 + t140 * t12 - t143 * (Ifges(4,1) * t142 - t249) / 0.2e1 - qJD(3) * (Ifges(4,5) * t142 - Ifges(4,6) * t143) / 0.2e1 - t157 * (mrSges(4,1) * t143 + mrSges(4,2) * t142) + Ifges(4,3) * qJDD(3); t207 * g(3) + (t196 * t11 - t114 * t34 + t192 * t12 + (-t192 * t61 + t196 * t60) * qJD(5) + (-g(3) * t171 - t114 * t80 + t170 * t269 - t18 * t233 + t19 * t232 + t192 * t2 + t196 * t3) * m(6)) * pkin(4) + t269 * t203 + t64 * t261 - m(6) * (t18 * t20 + t19 * t21) - t21 * t60 - t20 * t61 - t55 * t97 + t56 * t98 + t303; -g(3) * t217 - t18 * t60 + t19 * t61 + t269 * t210 + t304;];
tau = t1;
