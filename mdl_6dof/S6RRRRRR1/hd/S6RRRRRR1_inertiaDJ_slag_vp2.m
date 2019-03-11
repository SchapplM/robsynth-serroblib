% Calculate time derivative of joint inertia matrix for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:55
% EndTime: 2019-03-10 03:28:08
% DurationCPUTime: 5.92s
% Computational Cost: add. (18997->428), mult. (40576->647), div. (0->0), fcn. (41709->10), ass. (0->198)
t257 = qJD(2) + qJD(3);
t144 = sin(qJ(6));
t142 = t144 ^ 2;
t149 = cos(qJ(6));
t143 = t149 ^ 2;
t205 = t142 + t143;
t197 = qJD(6) * t149;
t145 = sin(qJ(5));
t150 = cos(qJ(5));
t146 = sin(qJ(4));
t151 = cos(qJ(4));
t147 = sin(qJ(3));
t148 = sin(qJ(2));
t152 = cos(qJ(3));
t153 = cos(qJ(2));
t117 = -t147 * t148 + t152 * t153;
t98 = t257 * t117;
t118 = t147 * t153 + t152 * t148;
t99 = t257 * t118;
t175 = t146 * t99 - t151 * t98;
t93 = t117 * t151 - t118 * t146;
t55 = qJD(4) * t93 - t175;
t94 = t117 * t146 + t118 * t151;
t56 = -qJD(4) * t94 - t146 * t98 - t151 * t99;
t66 = t145 * t94 - t150 * t93;
t27 = -qJD(5) * t66 + t145 * t56 + t150 * t55;
t222 = t144 * t27;
t67 = t145 * t93 + t150 * t94;
t171 = t67 * t197 + t222;
t245 = -pkin(8) - pkin(7);
t130 = t245 * t148;
t261 = t245 * t153;
t263 = -t147 * t130 + t152 * t261;
t260 = t152 * t130 + t147 * t261;
t172 = -pkin(9) * t118 + t260;
t83 = pkin(9) * t117 - t263;
t63 = -t146 * t83 + t151 * t172;
t164 = -t94 * pkin(10) + t63;
t64 = t146 * t172 + t151 * t83;
t51 = pkin(10) * t93 + t64;
t34 = t145 * t164 + t150 * t51;
t139 = -pkin(2) * t153 - pkin(1);
t103 = -t117 * pkin(3) + t139;
t72 = -t93 * pkin(4) + t103;
t39 = t66 * pkin(5) - t67 * pkin(11) + t72;
t17 = t144 * t39 + t149 * t34;
t206 = t17 * qJD(6);
t193 = t152 * t245;
t194 = t147 * t245;
t157 = ((-t146 * t194 + t151 * t193) * t153 + (-t146 * t193 - t151 * t194) * t148) * qJD(2);
t166 = qJD(4) * t172;
t201 = qJD(4) * t151;
t161 = -t55 * pkin(10) - t146 * t166 - t201 * t83;
t243 = t99 * pkin(9);
t69 = t257 * t260;
t158 = t69 - t243;
t244 = t98 * pkin(9);
t70 = t257 * t263;
t159 = t70 - t244;
t202 = qJD(4) * t146;
t30 = t146 * t159 - t202 * t83 + (t158 + t166) * t151;
t19 = pkin(10) * t56 + t30;
t203 = qJD(3) * t152;
t204 = qJD(3) * t147;
t33 = t145 * t51 - t150 * t164;
t10 = -t33 * qJD(5) + t150 * t19 + (-t146 * (t130 * t203 + t204 * t261 - t243) + t151 * (-t130 * t204 + t203 * t261 - t244) + t161 + t157) * t145;
t28 = qJD(5) * t67 + t145 * t55 - t150 * t56;
t85 = qJD(2) * t148 * pkin(2) + pkin(3) * t99;
t43 = -pkin(4) * t56 + t85;
t13 = pkin(5) * t28 - pkin(11) * t27 + t43;
t3 = -t10 * t144 + t13 * t149 - t206;
t262 = -t3 - t206;
t138 = pkin(2) * t152 + pkin(3);
t209 = t146 * t147;
t110 = -pkin(2) * t209 + t151 * t138;
t107 = pkin(4) + t110;
t207 = t147 * t151;
t112 = pkin(2) * t207 + t138 * t146;
t80 = t145 * t107 + t150 * t112;
t259 = t205 * t150;
t135 = pkin(4) * t145 + pkin(11);
t258 = t205 * t135;
t31 = t175 * pkin(9) - t64 * qJD(4) + (-t146 * t260 + t151 * t263) * qJD(3) + t157;
t256 = t31 * mrSges(5,1) - t30 * mrSges(5,2) + Ifges(5,5) * t55 + Ifges(5,6) * t56;
t255 = t70 * mrSges(4,1) - t69 * mrSges(4,2) + Ifges(4,5) * t98 - Ifges(4,6) * t99;
t180 = mrSges(7,1) * t144 + mrSges(7,2) * t149;
t122 = t180 * qJD(6);
t136 = -pkin(4) * t150 - pkin(5);
t105 = t136 * t122;
t125 = -mrSges(7,1) * t149 + mrSges(7,2) * t144;
t200 = qJD(5) * t145;
t115 = pkin(4) * t125 * t200;
t199 = qJD(5) * t150;
t195 = pkin(4) * t199;
t182 = mrSges(7,3) * t195;
t128 = t142 * t182;
t129 = t143 * t182;
t254 = t105 + t115 + t128 + t129;
t137 = pkin(3) * t151 + pkin(4);
t208 = t146 * t150;
t89 = t137 * t200 + (t146 * t199 + (t145 * t151 + t208) * qJD(4)) * pkin(3);
t73 = t89 * t125;
t229 = mrSges(7,3) * t142;
t210 = t145 * t146;
t88 = t137 * t199 + (-t146 * t200 + (t150 * t151 - t210) * qJD(4)) * pkin(3);
t81 = t88 * t229;
t228 = mrSges(7,3) * t143;
t82 = t88 * t228;
t86 = t89 * mrSges(6,1);
t109 = -pkin(3) * t210 + t137 * t150;
t106 = -pkin(5) - t109;
t97 = t106 * t122;
t253 = t73 + t81 + t82 + t97 - t86;
t16 = -t144 * t34 + t149 * t39;
t219 = t149 * t17;
t221 = t144 * t67;
t41 = -mrSges(7,2) * t66 - mrSges(7,3) * t221;
t217 = t149 * t67;
t42 = mrSges(7,1) * t66 - mrSges(7,3) * t217;
t252 = -t66 * mrSges(6,3) - t144 * t42 + t149 * t41 + m(7) * (-t144 * t16 + t219) + m(6) * t34;
t251 = 2 * m(5);
t250 = 0.2e1 * m(6);
t249 = 0.2e1 * m(7);
t11 = t145 * t19 - t150 * (-t146 * t158 + t151 * t159 + t161) + t34 * qJD(5);
t248 = 0.2e1 * t11;
t247 = 0.2e1 * t43;
t246 = 0.2e1 * t139;
t241 = pkin(5) * t122;
t240 = t11 * t33;
t2 = qJD(6) * t16 + t10 * t149 + t13 * t144;
t239 = t149 * t2;
t238 = t3 * t144;
t90 = t138 * t201 + (-t147 * t202 + (t151 * t152 - t209) * qJD(3)) * pkin(2);
t91 = -t138 * t202 + (-t147 * t201 + (-t146 * t152 - t207) * qJD(3)) * pkin(2);
t50 = qJD(5) * t80 + t145 * t90 - t150 * t91;
t236 = t33 * t50;
t235 = t33 * t89;
t79 = t107 * t150 - t112 * t145;
t49 = qJD(5) * t79 + t145 * t91 + t150 * t90;
t234 = t49 * mrSges(6,2);
t232 = t88 * mrSges(6,2);
t231 = t90 * mrSges(5,2);
t218 = t149 * t27;
t230 = Ifges(7,5) * t218 + Ifges(7,3) * t28;
t227 = Ifges(7,4) * t144;
t226 = Ifges(7,4) * t149;
t225 = Ifges(7,6) * t144;
t224 = pkin(4) * qJD(5);
t78 = pkin(11) + t80;
t216 = t149 * t78;
t111 = pkin(3) * t208 + t145 * t137;
t198 = qJD(6) * t144;
t196 = 0.2e1 * t153;
t191 = t67 * t198;
t189 = t16 * t197;
t170 = t191 - t218;
t12 = mrSges(7,1) * t171 - mrSges(7,2) * t170;
t188 = m(7) * t11 + t12;
t187 = pkin(11) * t205;
t186 = -t198 / 0.2e1;
t185 = -(2 * Ifges(6,4)) - t225;
t184 = t205 * t78;
t108 = pkin(11) + t111;
t183 = t108 * t205;
t181 = -t145 * mrSges(6,1) - t150 * mrSges(6,2);
t179 = Ifges(7,1) * t149 - t227;
t178 = -Ifges(7,2) * t144 + t226;
t177 = Ifges(7,5) * t144 + Ifges(7,6) * t149;
t123 = t178 * qJD(6);
t124 = t179 * qJD(6);
t126 = Ifges(7,2) * t149 + t227;
t127 = Ifges(7,1) * t144 + t226;
t169 = t149 * t123 + t144 * t124 - t126 * t198 + t127 * t197;
t168 = (-mrSges(4,1) * t147 - mrSges(4,2) * t152) * qJD(3) * pkin(2);
t167 = (-mrSges(5,1) * t146 - mrSges(5,2) * t151) * qJD(4) * pkin(3);
t44 = t50 * t125;
t45 = t49 * t229;
t46 = t49 * t228;
t47 = t50 * mrSges(6,1);
t77 = -pkin(5) - t79;
t71 = t77 * t122;
t165 = t169 + t44 + t45 + t46 - t47 + t71;
t163 = t169 - t232 + t253;
t87 = t91 * mrSges(5,1);
t162 = t165 + t87 - t231;
t140 = Ifges(7,5) * t197;
t37 = Ifges(7,6) * t66 + t178 * t67;
t38 = Ifges(7,5) * t66 + t179 * t67;
t7 = -Ifges(7,4) * t170 - Ifges(7,2) * t171 + Ifges(7,6) * t28;
t8 = -Ifges(7,1) * t170 - Ifges(7,4) * t171 + Ifges(7,5) * t28;
t160 = -t10 * mrSges(6,2) + mrSges(7,3) * t239 + Ifges(6,5) * t27 + t33 * t122 + t37 * t186 + t38 * t197 / 0.2e1 + t144 * t8 / 0.2e1 + t149 * t7 / 0.2e1 - t123 * t221 / 0.2e1 + t124 * t217 / 0.2e1 + t66 * (-Ifges(7,6) * t198 + t140) / 0.2e1 + (t177 / 0.2e1 - Ifges(6,6)) * t28 - t171 * t126 / 0.2e1 + (t218 / 0.2e1 + t67 * t186) * t127 + (t125 - mrSges(6,1)) * t11;
t14 = mrSges(7,1) * t28 + mrSges(7,3) * t170;
t15 = -mrSges(7,2) * t28 - mrSges(7,3) * t171;
t156 = -t42 * t197 - t41 * t198 + m(7) * (-t17 * t198 - t189 - t238 + t239) + t149 * t15 - t144 * t14;
t155 = t160 + (-t238 + (-t144 * t17 - t149 * t16) * qJD(6)) * mrSges(7,3);
t154 = t155 + t256;
t40 = t180 * t67;
t1 = [0.2e1 * m(4) * (t260 * t70 - t263 * t69) + 0.2e1 * (t117 * t69 - t118 * t70 - t260 * t98 + t263 * t99) * mrSges(4,3) + t40 * t248 + (t16 * t3 + t17 * t2 + t240) * t249 + (t10 * t34 + t43 * t72 + t240) * t250 + (t103 * t85 + t30 * t64 + t31 * t63) * t251 + t38 * t218 + (mrSges(4,1) * t99 + mrSges(4,2) * t98) * t246 + 0.2e1 * (t117 * t98 - t118 * t99) * Ifges(4,4) + 0.2e1 * t103 * (-mrSges(5,1) * t56 + mrSges(5,2) * t55) + 0.2e1 * t85 * (-mrSges(5,1) * t93 + mrSges(5,2) * t94) + 0.2e1 * t72 * (mrSges(6,1) * t28 + mrSges(6,2) * t27) + 0.2e1 * t2 * t41 + 0.2e1 * t3 * t42 + 0.2e1 * t33 * t12 + 0.2e1 * t16 * t14 + 0.2e1 * t17 * t15 + 0.2e1 * t98 * t118 * Ifges(4,1) + 0.2e1 * t55 * t94 * Ifges(5,1) - 0.2e1 * t117 * Ifges(4,2) * t99 + 0.2e1 * t93 * Ifges(5,2) * t56 + 0.2e1 * (t27 * t33 - t28 * t34) * mrSges(6,3) + 0.2e1 * (t55 * t93 + t56 * t94) * Ifges(5,4) + 0.2e1 * (t30 * t93 - t31 * t94 - t55 * t63 + t56 * t64) * mrSges(5,3) - t37 * t222 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t153) * t196 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t117 + mrSges(4,2) * t118) + m(4) * pkin(2) * t246 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t148 + (Ifges(3,1) - Ifges(3,2)) * t196) * t148) * qJD(2) + (mrSges(6,1) * t247 - 0.2e1 * mrSges(6,3) * t10 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t28 + t185 * t27 + t230) * t66 + (mrSges(6,2) * t247 + mrSges(6,3) * t248 + 0.2e1 * Ifges(6,1) * t27 - t144 * t7 + t149 * t8 + (Ifges(7,5) * t149 + t185) * t28 + (-t144 * t38 - t149 * t37 - t177 * t66) * qJD(6)) * t67; (-t110 * t55 + t112 * t56 + t90 * t93 - t91 * t94) * mrSges(5,3) + t255 + (-t79 * t27 - t80 * t28 - t49 * t66 + t50 * t67) * mrSges(6,3) + m(5) * (t110 * t31 + t112 * t30 + t63 * t91 + t64 * t90) + (m(4) * (t147 * t69 + t152 * t70 + (-t147 * t260 - t152 * t263) * qJD(3)) + (-t147 * t99 - t152 * t98 + (t117 * t152 + t118 * t147) * qJD(3)) * mrSges(4,3)) * pkin(2) + t160 + (Ifges(3,5) * t153 - Ifges(3,6) * t148 + (-mrSges(3,1) * t153 + mrSges(3,2) * t148) * pkin(7)) * qJD(2) + (t78 * t15 + t49 * t41 + (-mrSges(7,3) * t16 - t42 * t78) * qJD(6)) * t149 + t77 * t12 + t50 * t40 + (t262 * mrSges(7,3) + (-m(7) * t16 - t42) * t49 + (m(7) * t262 - qJD(6) * t41 - t14) * t78) * t144 + m(6) * (t10 * t80 - t11 * t79 + t34 * t49 + t236) + m(7) * (t11 * t77 - t78 * t189 + t2 * t216 + t49 * t219 + t236) + t256; -0.2e1 * t231 - 0.2e1 * t234 + 0.2e1 * t44 + 0.2e1 * t45 + 0.2e1 * t46 - 0.2e1 * t47 + 0.2e1 * t71 + 0.2e1 * t87 + 0.2e1 * t168 + (t184 * t49 + t50 * t77) * t249 + (t49 * t80 - t50 * t79) * t250 + (t110 * t91 + t112 * t90) * t251 + t169; t156 * t108 + (m(5) * (t146 * t30 + t151 * t31 + (-t146 * t63 + t151 * t64) * qJD(4)) + (t146 * t56 - t151 * t55 + (t146 * t94 + t151 * t93) * qJD(4)) * mrSges(5,3)) * pkin(3) + m(6) * (t10 * t111 - t109 * t11 + t235) + t154 + t106 * t12 + t89 * t40 + m(7) * (t106 * t11 + t235) + (-t109 * t27 - t111 * t28 + t67 * t89) * mrSges(6,3) + t252 * t88 + t255; (-mrSges(5,2) * t201 - mrSges(5,1) * t202 + m(5) * (-t110 * t202 + t112 * t201 + t146 * t90 + t151 * t91)) * pkin(3) + m(7) * (t106 * t50 + t183 * t49 + t184 * t88 + t77 * t89) + t162 + t168 + m(6) * (-t109 * t50 + t111 * t49 - t79 * t89 + t80 * t88) + (-t88 - t49) * mrSges(6,2) + t253; -0.2e1 * t232 + 0.2e1 * t73 + 0.2e1 * t81 + 0.2e1 * t82 - 0.2e1 * t86 + 0.2e1 * t97 + 0.2e1 * t167 + (t106 * t89 + t183 * t88) * t249 + (-t109 * t89 + t111 * t88) * t250 + t169; (m(6) * (t10 * t145 - t11 * t150) + (-t145 * t28 - t150 * t27) * mrSges(6,3) + (t252 * t150 + (t67 * mrSges(6,3) + t40 + (m(7) + m(6)) * t33) * t145) * qJD(5)) * pkin(4) + t156 * t135 + t188 * t136 + t154; m(7) * (t136 * t50 + t258 * t49) + (m(6) * (t145 * t49 - t150 * t50) + (m(7) * (t145 * t77 + t259 * t78) + m(6) * (-t145 * t79 + t150 * t80) + t181) * qJD(5)) * pkin(4) + t162 - t234 + t254; m(7) * (t136 * t89 + t258 * t88) + (m(6) * (t145 * t88 - t150 * t89) + (m(7) * (t106 * t145 + t108 * t259) + m(6) * (-t109 * t145 + t111 * t150) + t181) * qJD(5)) * pkin(4) + t163 + t167 + t254; 0.2e1 * t105 + 0.2e1 * t115 + 0.2e1 * t128 + 0.2e1 * t129 + 0.2e1 * (m(7) * (t135 * t259 + t136 * t145) + t181) * t224 + t169; -pkin(5) * t188 + pkin(11) * t156 + t155; -t241 + m(7) * (-pkin(5) * t50 + t187 * t49) - t234 + t165; -t241 + m(7) * (-pkin(5) * t89 + t187 * t88) + t163; -t241 + (m(7) * (-pkin(5) * t145 + pkin(11) * t259) + t181) * t224 + t169 + t254; t169 - 0.2e1 * t241; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t191 - Ifges(7,6) * t171 + t230; t140 - t180 * t49 + (-mrSges(7,1) * t216 + (mrSges(7,2) * t78 - Ifges(7,6)) * t144) * qJD(6); t140 - t180 * t88 + (t108 * t125 - t225) * qJD(6); t140 - t180 * t195 + (t125 * t135 - t225) * qJD(6); t140 + (pkin(11) * t125 - t225) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
