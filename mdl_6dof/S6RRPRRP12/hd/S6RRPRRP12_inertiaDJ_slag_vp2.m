% Calculate time derivative of joint inertia matrix for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:26
% EndTime: 2019-03-09 12:50:36
% DurationCPUTime: 4.09s
% Computational Cost: add. (4104->429), mult. (8726->602), div. (0->0), fcn. (6955->6), ass. (0->190)
t262 = Ifges(7,4) + Ifges(6,5);
t260 = Ifges(6,6) - Ifges(7,6);
t156 = cos(qJ(4));
t155 = sin(qJ(2));
t209 = qJD(2) * t155;
t194 = t156 * t209;
t154 = sin(qJ(4));
t157 = cos(qJ(2));
t205 = qJD(4) * t157;
t197 = t154 * t205;
t166 = t194 + t197;
t265 = -2 * mrSges(6,3) - 2 * mrSges(7,2);
t264 = mrSges(6,1) + mrSges(7,1);
t263 = -mrSges(6,2) + mrSges(7,3);
t261 = Ifges(7,2) + Ifges(6,3);
t153 = sin(qJ(5));
t204 = qJD(5) * t153;
t206 = qJD(4) * t156;
t236 = cos(qJ(5));
t192 = qJD(5) * t236;
t254 = t236 * qJD(4) + t192;
t74 = -t153 * t206 - t154 * t254 - t156 * t204;
t207 = qJD(4) * t154;
t75 = -t153 * t207 - t154 * t204 + t156 * t254;
t259 = -t260 * t75 + t262 * t74;
t258 = m(4) * pkin(7);
t257 = m(6) + m(7);
t171 = -t153 * t156 - t154 * t236;
t96 = t171 * t157;
t88 = mrSges(6,1) * t155 - mrSges(6,3) * t96;
t89 = -mrSges(7,1) * t155 + mrSges(7,2) * t96;
t227 = t88 - t89;
t198 = t236 * t156;
t108 = t153 * t154 - t198;
t255 = -t108 * t74 - t171 * t75;
t253 = -mrSges(5,1) * t207 - mrSges(5,2) * t206;
t158 = -pkin(2) - pkin(8);
t215 = qJ(3) * t155;
t104 = t157 * t158 - pkin(1) - t215;
t208 = qJD(2) * t157;
t240 = pkin(3) + pkin(7);
t118 = t240 * t208;
t125 = t240 * t155;
t186 = pkin(2) * t209 - qJD(3) * t155;
t92 = (pkin(8) * t155 - qJ(3) * t157) * qJD(2) + t186;
t25 = -t104 * t207 + t154 * t118 + t125 * t206 + t156 * t92;
t187 = t156 * t118 - t154 * t92;
t110 = t154 * t125;
t73 = t156 * t104 + t110;
t26 = -qJD(4) * t73 + t187;
t252 = t154 * t25 + t156 * t26;
t251 = -t153 * t75 - t236 * t74;
t250 = qJD(4) + qJD(5);
t234 = pkin(9) - t158;
t191 = t234 * t156;
t107 = qJD(4) * t191;
t183 = t234 * t207;
t119 = t234 * t154;
t82 = -t153 * t119 + t191 * t236;
t33 = -qJD(5) * t82 - t107 * t236 + t153 * t183;
t83 = -t119 * t236 - t153 * t191;
t34 = qJD(5) * t83 - t153 * t107 - t183 * t236;
t249 = t108 * t34 - t171 * t33 - t74 * t82 + t75 * t83;
t189 = pkin(9) * t157 - t104;
t15 = (-pkin(9) * t154 * t155 + pkin(4) * t157) * qJD(2) + (t156 * t189 - t110) * qJD(4) + t187;
t19 = pkin(9) * t166 + t25;
t111 = t156 * t125;
t57 = pkin(4) * t155 + t154 * t189 + t111;
t210 = t156 * t157;
t63 = -pkin(9) * t210 + t73;
t231 = t153 * t57 + t236 * t63;
t6 = -qJD(5) * t231 + t15 * t236 - t153 * t19;
t195 = t154 * t209;
t248 = Ifges(5,5) * t195 + t166 * Ifges(5,6) + Ifges(5,3) * t208;
t184 = pkin(4) * t192;
t131 = t184 + qJD(6);
t235 = pkin(4) * t153;
t140 = qJ(6) + t235;
t203 = t236 * pkin(4);
t143 = -t203 - pkin(5);
t202 = pkin(4) * t204;
t185 = t108 * t202;
t247 = -t131 * t171 + t140 * t75 - t143 * t74 + t185;
t246 = 2 * m(6);
t245 = 2 * m(7);
t244 = -0.2e1 * pkin(1);
t98 = -qJ(3) * t208 + t186;
t243 = 0.2e1 * t98;
t177 = -pkin(2) * t157 - t215;
t120 = -pkin(1) + t177;
t242 = -0.2e1 * t120;
t241 = m(6) * pkin(4);
t237 = -t157 / 0.2e1;
t150 = t157 * pkin(7);
t47 = -t153 * t195 + t236 * t194 - t250 * t96;
t27 = mrSges(7,2) * t47 + mrSges(7,3) * t208;
t30 = -mrSges(6,2) * t208 + mrSges(6,3) * t47;
t233 = t27 + t30;
t46 = t108 * t157 * t250 - t171 * t209;
t28 = mrSges(6,1) * t208 - mrSges(6,3) * t46;
t29 = -mrSges(7,1) * t208 + t46 * mrSges(7,2);
t232 = -t28 + t29;
t212 = t154 * t157;
t95 = t153 * t212 - t157 * t198;
t86 = mrSges(7,2) * t95 + mrSges(7,3) * t155;
t87 = -mrSges(6,2) * t155 + mrSges(6,3) * t95;
t228 = t86 + t87;
t226 = Ifges(5,4) * t154;
t225 = Ifges(5,4) * t156;
t180 = Ifges(5,1) * t154 + t225;
t218 = t155 * Ifges(5,5);
t94 = -t157 * t180 + t218;
t219 = t154 * t94;
t179 = Ifges(5,2) * t156 + t226;
t93 = t155 * Ifges(5,6) - t157 * t179;
t216 = t156 * t93;
t124 = Ifges(5,1) * t156 - t226;
t213 = t154 * t124;
t123 = -Ifges(5,2) * t154 + t225;
t211 = t156 * t123;
t141 = t154 * pkin(4) + qJ(3);
t126 = t157 * pkin(3) + t150;
t132 = pkin(4) * t206 + qJD(3);
t100 = pkin(4) * t210 + t126;
t196 = t156 * t205;
t193 = t83 * t33 + t34 * t82;
t181 = mrSges(5,1) * t156 - mrSges(5,2) * t154;
t178 = -Ifges(5,5) * t154 - Ifges(5,6) * t156;
t72 = -t104 * t154 + t111;
t176 = t72 * t154 - t73 * t156;
t115 = mrSges(5,1) * t155 + mrSges(5,3) * t212;
t116 = -mrSges(5,2) * t155 - mrSges(5,3) * t210;
t175 = -t154 * t115 + t156 * t116;
t174 = t261 * t208 + t260 * t47 + t262 * t46;
t21 = -t153 * t63 + t236 * t57;
t5 = t153 * t15 + t236 * t19 + t57 * t192 - t204 * t63;
t170 = t263 * t75 + t264 * t74;
t167 = pkin(5) * t74 + qJ(6) * t75 - qJD(6) * t171;
t165 = t195 - t196;
t18 = qJ(6) * t155 + t231;
t2 = qJ(6) * t208 + qJD(6) * t155 + t5;
t20 = -t155 * pkin(5) - t21;
t3 = -pkin(5) * t208 - t6;
t164 = t108 * t3 - t171 * t2 + t18 * t75 - t20 * t74;
t163 = t108 * t6 + t171 * t5 - t21 * t74 - t231 * t75;
t162 = t263 * t33 - t264 * t34 + t259;
t161 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t174;
t84 = -pkin(4) * t197 + (-pkin(4) * t156 - t240) * t209;
t159 = -mrSges(6,2) * t184 + t131 * mrSges(7,3) - t264 * t202;
t152 = qJD(6) * mrSges(7,3);
t122 = mrSges(5,1) * t154 + mrSges(5,2) * t156;
t117 = t240 * t209;
t114 = t180 * qJD(4);
t113 = t179 * qJD(4);
t112 = t181 * qJD(4);
t103 = t181 * t157;
t91 = mrSges(5,1) * t208 - mrSges(5,3) * t165;
t90 = -mrSges(5,2) * t208 + mrSges(5,3) * t166;
t81 = -Ifges(6,1) * t108 + Ifges(6,4) * t171;
t80 = -Ifges(7,1) * t108 - Ifges(7,5) * t171;
t79 = -Ifges(6,4) * t108 + Ifges(6,2) * t171;
t78 = -Ifges(7,5) * t108 - Ifges(7,3) * t171;
t77 = -mrSges(6,1) * t171 - mrSges(6,2) * t108;
t76 = -mrSges(7,1) * t171 + mrSges(7,3) * t108;
t65 = -mrSges(5,1) * t166 + mrSges(5,2) * t165;
t64 = -pkin(5) * t171 + qJ(6) * t108 + t141;
t62 = -mrSges(6,1) * t95 + mrSges(6,2) * t96;
t61 = -mrSges(7,1) * t95 - mrSges(7,3) * t96;
t60 = -t124 * t205 + (t157 * Ifges(5,5) + t155 * t180) * qJD(2);
t59 = -t123 * t205 + (t157 * Ifges(5,6) + t155 * t179) * qJD(2);
t55 = Ifges(6,1) * t96 + Ifges(6,4) * t95 + Ifges(6,5) * t155;
t54 = Ifges(7,1) * t96 + Ifges(7,4) * t155 - Ifges(7,5) * t95;
t53 = Ifges(6,4) * t96 + Ifges(6,2) * t95 + Ifges(6,6) * t155;
t52 = Ifges(7,5) * t96 + Ifges(7,6) * t155 - Ifges(7,3) * t95;
t49 = -pkin(5) * t95 - qJ(6) * t96 + t100;
t40 = Ifges(6,1) * t74 - Ifges(6,4) * t75;
t39 = Ifges(7,1) * t74 + Ifges(7,5) * t75;
t38 = Ifges(6,4) * t74 - Ifges(6,2) * t75;
t37 = Ifges(7,5) * t74 + Ifges(7,3) * t75;
t36 = mrSges(6,1) * t75 + mrSges(6,2) * t74;
t35 = mrSges(7,1) * t75 - mrSges(7,3) * t74;
t23 = pkin(5) * t75 - qJ(6) * t74 + qJD(6) * t108 + t132;
t13 = -mrSges(6,1) * t47 + mrSges(6,2) * t46;
t12 = -mrSges(7,1) * t47 - mrSges(7,3) * t46;
t11 = Ifges(6,1) * t46 + Ifges(6,4) * t47 + Ifges(6,5) * t208;
t10 = Ifges(7,1) * t46 + Ifges(7,4) * t208 - Ifges(7,5) * t47;
t9 = Ifges(6,4) * t46 + Ifges(6,2) * t47 + Ifges(6,6) * t208;
t8 = Ifges(7,5) * t46 + Ifges(7,6) * t208 - Ifges(7,3) * t47;
t7 = -pkin(5) * t47 - qJ(6) * t46 - qJD(6) * t96 + t84;
t1 = [0.2e1 * m(5) * (-t117 * t126 + t25 * t73 + t26 * t72) + (t100 * t84 + t21 * t6 + t231 * t5) * t246 + 0.2e1 * t231 * t30 + (t10 + t11) * t96 + (t9 - t8) * t95 + (-0.2e1 * t98 * mrSges(4,3) + t174 + t248) * t155 + m(4) * t120 * t243 + (-t52 + t53) * t47 + (t54 + t55) * t46 + (t18 * t2 + t20 * t3 + t49 * t7) * t245 + 0.2e1 * t126 * t65 + 0.2e1 * t26 * t115 + 0.2e1 * t25 * t116 - 0.2e1 * t117 * t103 + 0.2e1 * t100 * t13 + 0.2e1 * t2 * t86 + 0.2e1 * t5 * t87 + 0.2e1 * t6 * t88 + 0.2e1 * t3 * t89 + 0.2e1 * t73 * t90 + 0.2e1 * t72 * t91 + 0.2e1 * t84 * t62 + 0.2e1 * t49 * t12 + 0.2e1 * t7 * t61 + 0.2e1 * t18 * t27 + 0.2e1 * t21 * t28 + 0.2e1 * t20 * t29 + ((mrSges(3,1) * t244 + mrSges(4,2) * t242 + t219 + t216 + 0.2e1 * (-Ifges(3,4) - Ifges(4,6)) * t155) * t155 + (mrSges(3,2) * t244 + mrSges(4,3) * t242 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t178) * t157 + t262 * t96 + t260 * t95 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3) + t261) * t155) * t157) * qJD(2) + (mrSges(4,2) * t243 - t154 * t60 - t156 * t59 + (t154 * t93 + (-t94 - t218) * t156) * qJD(4)) * t157; m(5) * (-qJ(3) * t117 + t158 * t252) + (t60 / 0.2e1 - t26 * mrSges(5,3) - t113 * t237 + t158 * t91) * t156 + (-t59 / 0.2e1 - t25 * mrSges(5,3) - t114 * t237 + t158 * t90) * t154 + m(6) * (t100 * t132 + t141 * t84 - t21 * t34 + t231 * t33 + t5 * t83 - t6 * t82) + (t178 * qJD(4) + t259) * t155 / 0.2e1 + (-t10 / 0.2e1 - t11 / 0.2e1) * t108 + (m(4) * t150 + m(5) * t126 + t157 * mrSges(4,1) + t103) * qJD(3) - (t8 / 0.2e1 - t9 / 0.2e1) * t171 + (t177 * t258 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) + Ifges(3,5) + Ifges(5,5) * t156 / 0.2e1 - Ifges(5,6) * t154 / 0.2e1 - (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t171 + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t108 + (-mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t157 + (t211 / 0.2e1 - qJ(3) * mrSges(4,1) + t213 / 0.2e1 + Ifges(4,5) - Ifges(3,6) + (mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t155) * qJD(2) - t164 * mrSges(7,2) + t163 * mrSges(6,3) + (-t219 / 0.2e1 - t216 / 0.2e1 + (-t156 * t124 / 0.2e1 + t154 * t123 / 0.2e1) * t157 + t176 * mrSges(5,3) + (-m(5) * t176 + t175) * t158) * qJD(4) + m(7) * (t18 * t33 + t2 * t83 + t20 * t34 + t23 * t49 + t3 * t82 + t64 * t7) + (t39 / 0.2e1 + t40 / 0.2e1) * t96 + (t38 / 0.2e1 - t37 / 0.2e1) * t95 + t132 * t62 + t141 * t13 - t117 * t122 + t126 * t112 + t100 * t36 + t7 * t76 + t84 * t77 + t64 * t12 + qJ(3) * t65 + t49 * t35 + t23 * t61 + (-t78 / 0.2e1 + t79 / 0.2e1) * t47 + (t80 / 0.2e1 + t81 / 0.2e1) * t46 + (t54 / 0.2e1 + t55 / 0.2e1) * t74 + (t52 / 0.2e1 - t53 / 0.2e1) * t75 - t227 * t34 + t228 * t33 + t232 * t82 + t233 * t83; 0.2e1 * qJ(3) * t112 + t154 * t113 - t156 * t114 + 0.2e1 * t132 * t77 + 0.2e1 * t141 * t36 + 0.2e1 * t23 * t76 + 0.2e1 * t64 * t35 + (t78 - t79) * t75 + (t80 + t81) * t74 - (t37 - t38) * t171 + (-t39 - t40) * t108 + (-t211 - t213) * qJD(4) + (t23 * t64 + t193) * t245 + (t132 * t141 + t193) * t246 + 0.2e1 * (mrSges(4,3) + t122 + (m(4) + m(5)) * qJ(3)) * qJD(3) + t249 * t265; t154 * t90 + t156 * t91 + t228 * t75 + t227 * t74 - t233 * t171 + t232 * t108 + t175 * qJD(4) + (mrSges(4,1) + t258) * t208 + m(7) * t164 - m(6) * t163 + m(5) * (-qJD(4) * t176 + t252); t249 * t257 + t255 * t265; 0.2e1 * t257 * t255; m(7) * (t131 * t18 + t140 * t2 + t143 * t3) + (t236 * t6 + t153 * t5 + (-t153 * t21 + t231 * t236) * qJD(5)) * t241 - Ifges(5,5) * t196 + t87 * t184 + t131 * t86 + t140 * t27 + t143 * t29 - t25 * mrSges(5,2) + t26 * mrSges(5,1) + t161 + t28 * t203 + t30 * t235 + (m(7) * t20 - t227) * t202 + t248; m(7) * (t131 * t83 + t140 * t33 + t143 * t34 + t202 * t82) + (-t236 * t34 + t153 * t33 + (t153 * t82 + t236 * t83) * qJD(5)) * t241 + t162 - Ifges(5,5) * t207 - Ifges(5,6) * t206 + t253 * t158 - t247 * mrSges(7,2) + (pkin(4) * t251 + t171 * t184 - t185) * mrSges(6,3); ((t108 * t153 - t171 * t236) * qJD(5) - t251) * t241 + m(7) * t247 + t170 + t253; 0.2e1 * m(7) * (t131 * t140 + t143 * t202) + 0.2e1 * t159; m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t18) + qJD(6) * t86 + qJ(6) * t27 - pkin(5) * t29 + t161; m(7) * (-pkin(5) * t34 + qJ(6) * t33 + qJD(6) * t83) - t167 * mrSges(7,2) + t162; m(7) * t167 + t170; t152 + m(7) * (-pkin(5) * t202 + qJ(6) * t131 + qJD(6) * t140) + t159; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t152; m(7) * t3 + t29; m(7) * t34 + t74 * mrSges(7,2); -m(7) * t74; m(7) * t202; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
