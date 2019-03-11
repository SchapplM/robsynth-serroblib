% Calculate time derivative of joint inertia matrix for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:09:05
% EndTime: 2019-03-09 05:09:09
% DurationCPUTime: 2.54s
% Computational Cost: add. (8419->339), mult. (18403->506), div. (0->0), fcn. (19383->10), ass. (0->142)
t147 = sin(pkin(11));
t149 = cos(pkin(11));
t206 = -mrSges(6,1) * t149 + mrSges(6,2) * t147 - mrSges(5,1);
t152 = sin(qJ(4));
t155 = cos(qJ(4));
t148 = sin(pkin(10));
t191 = pkin(7) + qJ(2);
t134 = t191 * t148;
t156 = cos(qJ(3));
t125 = t156 * t134;
t150 = cos(pkin(10));
t136 = t191 * t150;
t153 = sin(qJ(3));
t108 = -t136 * t153 - t125;
t130 = t148 * t156 + t150 * t153;
t95 = -pkin(8) * t130 + t108;
t110 = -t134 * t153 + t136 * t156;
t174 = t150 * t156;
t128 = -t148 * t153 + t174;
t96 = pkin(8) * t128 + t110;
t204 = -t152 * t96 + t155 * t95;
t103 = t128 * t152 + t130 * t155;
t168 = -pkin(2) * t150 - pkin(1);
t111 = -pkin(3) * t128 + t168;
t160 = t128 * t155 - t130 * t152;
t62 = -pkin(4) * t160 - qJ(5) * t103 + t111;
t64 = t152 * t95 + t155 * t96;
t46 = -t147 * t64 + t149 * t62;
t47 = t147 * t62 + t149 * t64;
t203 = -t147 * t46 + t149 * t47;
t151 = sin(qJ(6));
t154 = cos(qJ(6));
t159 = t147 * t151 - t149 * t154;
t117 = t159 * qJD(6);
t129 = t147 * t154 + t149 * t151;
t104 = mrSges(7,1) * t159 + mrSges(7,2) * t129;
t186 = pkin(3) * qJD(4);
t202 = (-mrSges(5,2) * t155 + (t104 + t206) * t152) * t186;
t201 = 2 * m(6);
t200 = 2 * m(7);
t199 = -0.2e1 * t204;
t118 = t129 * qJD(6);
t99 = mrSges(7,1) * t118 - mrSges(7,2) * t117;
t198 = 0.2e1 * t99;
t145 = t149 ^ 2;
t197 = 0.2e1 * t111;
t195 = pkin(3) * t155;
t120 = t130 * qJD(3);
t194 = t120 * pkin(3);
t119 = t128 * qJD(3);
t93 = -qJD(2) * t130 - qJD(3) * t110;
t158 = -t119 * pkin(8) + t93;
t92 = -qJD(3) * t125 + qJD(2) * t174 + (-qJD(2) * t148 - qJD(3) * t136) * t153;
t85 = -pkin(8) * t120 + t92;
t37 = qJD(4) * t204 + t152 * t158 + t155 * t85;
t77 = qJD(4) * t160 + t119 * t155 - t120 * t152;
t78 = qJD(4) * t103 + t119 * t152 + t120 * t155;
t43 = pkin(4) * t78 - qJ(5) * t77 - qJD(5) * t103 + t194;
t15 = -t147 * t37 + t149 * t43;
t193 = t15 * mrSges(6,3);
t38 = qJD(4) * t64 + t152 * t85 - t155 * t158;
t192 = t38 * t204;
t16 = t147 * t43 + t149 * t37;
t180 = t149 * t77;
t183 = t147 * t77;
t50 = mrSges(6,1) * t183 + mrSges(6,2) * t180;
t189 = Ifges(6,4) * t147;
t188 = Ifges(6,4) * t149;
t187 = Ifges(6,2) * t147;
t185 = t147 * t15;
t182 = t149 * t16;
t179 = t152 * t204;
t177 = t103 * t147;
t176 = t103 * t149;
t139 = pkin(3) * t152 + qJ(5);
t175 = t139 * t149;
t173 = -Ifges(7,5) * t117 - Ifges(7,6) * t118;
t172 = t147 ^ 2 + t145;
t171 = 2 * mrSges(7,3);
t27 = -t103 * t118 - t159 * t77;
t28 = t103 * t117 - t129 * t77;
t170 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t78;
t169 = t152 * t186;
t140 = -pkin(5) * t149 - pkin(4);
t13 = -t28 * mrSges(7,1) + mrSges(7,2) * t27;
t100 = -Ifges(7,4) * t117 - Ifges(7,2) * t118;
t101 = -Ifges(7,1) * t117 - Ifges(7,4) * t118;
t105 = Ifges(7,4) * t129 - Ifges(7,2) * t159;
t106 = Ifges(7,1) * t129 - Ifges(7,4) * t159;
t167 = -t100 * t159 + t101 * t129 - t105 * t118 - t106 * t117;
t137 = t155 * t186 + qJD(5);
t166 = t172 * t137;
t165 = t172 * qJD(5);
t162 = Ifges(6,5) * t149 - Ifges(6,6) * t147;
t23 = -pkin(5) * t160 - pkin(9) * t176 + t46;
t35 = -pkin(9) * t177 + t47;
t11 = -t151 * t35 + t154 * t23;
t12 = t151 * t23 + t154 * t35;
t161 = 0.2e1 * t172 * mrSges(6,3);
t121 = (-pkin(9) - t139) * t147;
t142 = t149 * pkin(9);
t122 = t142 + t175;
t97 = t121 * t154 - t122 * t151;
t98 = t121 * t151 + t122 * t154;
t132 = (-pkin(9) - qJ(5)) * t147;
t135 = qJ(5) * t149 + t142;
t107 = t132 * t154 - t135 * t151;
t109 = t132 * t151 + t135 * t154;
t4 = pkin(5) * t78 - pkin(9) * t180 + t15;
t9 = -pkin(9) * t183 + t16;
t2 = qJD(6) * t11 + t151 * t4 + t154 * t9;
t20 = pkin(5) * t183 + t38;
t3 = -qJD(6) * t12 - t151 * t9 + t154 * t4;
t32 = Ifges(6,6) * t78 + (-t187 + t188) * t77;
t33 = Ifges(6,5) * t78 + (Ifges(6,1) * t149 - t189) * t77;
t68 = t129 * t103;
t69 = t159 * t103;
t44 = -Ifges(7,4) * t69 - Ifges(7,2) * t68 - Ifges(7,6) * t160;
t45 = -Ifges(7,1) * t69 - Ifges(7,4) * t68 - Ifges(7,5) * t160;
t53 = pkin(5) * t177 - t204;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t78;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t78;
t157 = -(Ifges(6,2) * t149 + t189) * t183 / 0.2e1 + (Ifges(6,1) * t147 + t188) * t180 / 0.2e1 - t160 * t173 / 0.2e1 + mrSges(6,3) * t182 - t37 * mrSges(5,2) + Ifges(5,5) * t77 - Ifges(5,6) * t78 + t53 * t99 - t68 * t100 / 0.2e1 - t69 * t101 / 0.2e1 + t20 * t104 + t28 * t105 / 0.2e1 + t27 * t106 / 0.2e1 - t117 * t45 / 0.2e1 - t118 * t44 / 0.2e1 - t159 * t7 / 0.2e1 + t129 * t8 / 0.2e1 + t147 * t33 / 0.2e1 + t149 * t32 / 0.2e1 + t206 * t38 + (Ifges(6,5) * t147 + Ifges(7,5) * t129 + Ifges(6,6) * t149 - Ifges(7,6) * t159) * t78 / 0.2e1 + (t11 * t117 - t118 * t12 - t129 * t3 - t159 * t2) * mrSges(7,3);
t141 = -pkin(4) - t195;
t131 = t140 - t195;
t113 = t119 * mrSges(4,2);
t91 = -qJD(5) * t129 - qJD(6) * t109;
t90 = -qJD(5) * t159 + qJD(6) * t107;
t83 = -qJD(6) * t98 - t129 * t137;
t82 = qJD(6) * t97 - t137 * t159;
t80 = -mrSges(6,1) * t160 - mrSges(6,3) * t176;
t79 = mrSges(6,2) * t160 - mrSges(6,3) * t177;
t73 = t77 * mrSges(5,2);
t72 = (mrSges(6,1) * t147 + mrSges(6,2) * t149) * t103;
t57 = -mrSges(7,1) * t160 + mrSges(7,3) * t69;
t56 = mrSges(7,2) * t160 - mrSges(7,3) * t68;
t52 = mrSges(6,1) * t78 - mrSges(6,3) * t180;
t51 = -mrSges(6,2) * t78 - mrSges(6,3) * t183;
t49 = mrSges(7,1) * t68 - mrSges(7,2) * t69;
t19 = -mrSges(7,2) * t78 + mrSges(7,3) * t28;
t18 = mrSges(7,1) * t78 - mrSges(7,3) * t27;
t1 = [-0.2e1 * t128 * Ifges(4,2) * t120 + 0.2e1 * m(4) * (t108 * t93 + t110 * t92) + 0.2e1 * t130 * t119 * Ifges(4,1) - t32 * t177 + (mrSges(5,3) * t199 + (Ifges(6,1) * t145 + (2 * Ifges(5,1)) + (t187 - 0.2e1 * t188) * t147) * t103 - 0.2e1 * (-Ifges(5,4) + t162) * t160) * t77 - t160 * t170 + 0.2e1 * (t103 * t38 + t160 * t37) * mrSges(5,3) + 0.2e1 * (-mrSges(5,1) * t160 + mrSges(5,2) * t103) * t194 + (mrSges(5,1) * t197 - 0.2e1 * mrSges(5,3) * t64 - Ifges(7,5) * t69 - Ifges(7,6) * t68 + (-0.2e1 * Ifges(5,4) + t162) * t103 - ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t160) * t78 + (t15 * t46 + t16 * t47 - t192) * t201 + 0.2e1 * m(5) * (t111 * t194 + t37 * t64 - t192) + t33 * t176 + t73 * t197 + t50 * t199 + (t11 * t3 + t12 * t2 + t20 * t53) * t200 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t148 ^ 2 + t150 ^ 2) * qJD(2) + 0.2e1 * t11 * t18 + 0.2e1 * t12 * t19 + t28 * t44 + t27 * t45 + 0.2e1 * t20 * t49 + 0.2e1 * t47 * t51 + 0.2e1 * t46 * t52 + 0.2e1 * t53 * t13 + 0.2e1 * t2 * t56 + 0.2e1 * t3 * t57 - t68 * t7 - t69 * t8 + 0.2e1 * t38 * t72 + 0.2e1 * t16 * t79 + 0.2e1 * t15 * t80 + 0.2e1 * t168 * (t120 * mrSges(4,1) + t113) + 0.2e1 * (-t108 * t119 - t110 * t120 + t128 * t92 - t130 * t93) * mrSges(4,3) + 0.2e1 * (t119 * t128 - t120 * t130) * Ifges(4,4); t78 * mrSges(5,1) - t117 * t56 - t118 * t57 - t159 * t18 + t129 * t19 + t147 * t51 + t149 * t52 + t113 + t73 - (-m(5) * pkin(3) - mrSges(4,1)) * t120 + m(7) * (-t11 * t118 - t117 * t12 + t129 * t2 - t159 * t3) + m(6) * (t147 * t16 + t149 * t15); (-t117 * t129 + t118 * t159) * t200; t157 + (m(5) * (t152 * t37 - t155 * t38) + (-t152 * t78 - t155 * t77) * mrSges(5,3) + (t155 * t160 * mrSges(5,3) - m(6) * t179 + m(5) * (t155 * t64 - t179) + (m(7) * t53 + t103 * mrSges(5,3) + t49 + t72) * t152) * qJD(4)) * pkin(3) + (t137 * t79 + t139 * t51) * t149 + (-t137 * t80 - t139 * t52 - t193) * t147 + m(6) * (t137 * t203 - t139 * t185 + t141 * t38 + t16 * t175) + m(7) * (t11 * t83 + t12 * t82 + t131 * t20 + t2 * t98 + t3 * t97) + t82 * t56 + t83 * t57 - t92 * mrSges(4,2) + t93 * mrSges(4,1) + t97 * t18 + t98 * t19 + Ifges(4,5) * t119 - Ifges(4,6) * t120 + t131 * t13 + t141 * t50; m(7) * (-t117 * t98 - t118 * t97 + t129 * t82 - t159 * t83); t131 * t198 + t137 * t161 + 0.2e1 * t202 + (t131 * t169 + t82 * t98 + t83 * t97) * t200 + (t139 * t166 + t141 * t169) * t201 + (t97 * t117 - t98 * t118 - t83 * t129 - t159 * t82) * t171 + t167; t157 + (qJ(5) * t51 + qJD(5) * t79) * t149 + (-qJ(5) * t52 - qJD(5) * t80 - t193) * t147 + m(7) * (t107 * t3 + t109 * t2 + t11 * t91 + t12 * t90 + t140 * t20) + m(6) * (-pkin(4) * t38 + t203 * qJD(5) + (t182 - t185) * qJ(5)) - pkin(4) * t50 + t90 * t56 + t91 * t57 + t107 * t18 + t109 * t19 + t140 * t13; m(7) * (-t107 * t118 - t109 * t117 + t129 * t90 - t159 * t91); (t131 + t140) * t99 + t202 + m(7) * (t107 * t83 + t109 * t82 + t140 * t169 + t90 * t98 + t91 * t97) + m(6) * (-pkin(4) * t169 + qJ(5) * t166 + t139 * t165) + (t166 + t165) * mrSges(6,3) + ((-t83 - t91) * t129 - (t82 + t90) * t159 - (t109 + t98) * t118 - (-t107 - t97) * t117) * mrSges(7,3) + t167; t140 * t198 + (t107 * t91 + t109 * t90) * t200 + (qJ(5) * t172 * t201 + t161) * qJD(5) + (t107 * t117 - t109 * t118 - t91 * t129 - t159 * t90) * t171 + t167; m(6) * t38 + m(7) * t20 + t13 + t50; 0; (m(6) + m(7)) * t169 + t99; t99; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t170; -t99; mrSges(7,1) * t83 - mrSges(7,2) * t82 + t173; mrSges(7,1) * t91 - mrSges(7,2) * t90 + t173; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
