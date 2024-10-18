% Calculate time derivative of joint inertia matrix for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:24
% EndTime: 2024-09-27 17:30:27
% DurationCPUTime: 2.05s
% Computational Cost: add. (3587->277), mult. (9126->407), div. (0->0), fcn. (7474->10), ass. (0->171)
t133 = sin(qJ(3));
t172 = qJD(3) * t133;
t163 = pkin(2) * t172;
t129 = sin(pkin(5));
t136 = cos(qJ(4));
t169 = qJD(4) * t136;
t158 = t129 * t169;
t205 = Ifges(5,5) * t158;
t201 = qJD(4) + qJD(5);
t131 = sin(qJ(5));
t132 = sin(qJ(4));
t135 = cos(qJ(5));
t91 = (-t131 * t132 + t135 * t136) * t129;
t59 = t201 * t91;
t92 = (t131 * t136 + t132 * t135) * t129;
t60 = t201 * t92;
t186 = Ifges(6,5) * t59 - Ifges(6,6) * t60;
t128 = t129 ^ 2;
t153 = t128 * t163;
t190 = pkin(10) * t129;
t138 = cos(qJ(2));
t125 = pkin(1) * t138 + pkin(2);
t134 = sin(qJ(2));
t137 = cos(qJ(3));
t174 = t134 * t137;
t100 = pkin(1) * t174 + t125 * t133;
t126 = t129 * pkin(9);
t89 = t126 + t100;
t156 = -t89 - t190;
t149 = t156 * t132;
t130 = cos(pkin(5));
t157 = t130 * t169;
t177 = t130 * t132;
t171 = qJD(3) * t137;
t175 = t133 * t134;
t65 = t125 * t171 + (-t134 * t172 + (t137 * t138 - t175) * qJD(2)) * pkin(1);
t66 = -t125 * t172 + (-t134 * t171 + (-t133 * t138 - t174) * qJD(2)) * pkin(1);
t99 = -pkin(1) * t175 + t125 * t137;
t97 = pkin(3) + t99;
t164 = t136 * t65 + t157 * t97 + t177 * t66;
t16 = qJD(4) * t149 + t164;
t176 = t130 * t136;
t154 = -t132 * t65 + t176 * t66;
t85 = t97 * t177;
t17 = (t136 * t156 - t85) * qJD(4) + t154;
t127 = t130 * pkin(4);
t86 = t97 * t176;
t41 = t127 + t86 + t149;
t178 = t129 * t136;
t119 = pkin(10) * t178;
t50 = t136 * t89 + t85;
t47 = t119 + t50;
t20 = -t131 * t47 + t135 * t41;
t4 = qJD(5) * t20 + t131 * t17 + t135 * t16;
t21 = t131 * t41 + t135 * t47;
t5 = -qJD(5) * t21 - t131 * t16 + t135 * t17;
t204 = mrSges(6,1) * t5 - t4 * mrSges(6,2);
t124 = pkin(2) * t137 + pkin(3);
t111 = t124 * t176;
t112 = pkin(2) * t133 + t126;
t155 = -t112 - t190;
t64 = t132 * t155 + t111 + t127;
t110 = t124 * t177;
t79 = t136 * t112 + t110;
t67 = t119 + t79;
t38 = t131 * t64 + t135 * t67;
t152 = t130 * t163;
t162 = pkin(2) * t171;
t173 = t124 * t157 + t136 * t162;
t43 = (qJD(4) * t155 - t152) * t132 + t173;
t183 = pkin(2) * qJD(3);
t141 = (-t132 * t137 - t133 * t176) * t183;
t44 = t141 + (t136 * t155 - t110) * qJD(4);
t10 = -qJD(5) * t38 - t131 * t43 + t135 * t44;
t37 = -t131 * t67 + t135 * t64;
t9 = qJD(5) * t37 + t131 * t44 + t135 * t43;
t203 = t10 * mrSges(6,1) - t9 * mrSges(6,2);
t122 = pkin(3) * t176;
t161 = t129 * (-pkin(9) - pkin(10));
t150 = t132 * t161;
t73 = t122 + t127 + t150;
t121 = pkin(3) * t177;
t103 = pkin(9) * t178 + t121;
t87 = t119 + t103;
t45 = -t131 * t87 + t135 * t73;
t117 = pkin(3) * t157;
t74 = qJD(4) * t150 + t117;
t75 = (t136 * t161 - t121) * qJD(4);
t24 = qJD(5) * t45 + t131 * t75 + t135 * t74;
t46 = t131 * t73 + t135 * t87;
t25 = -qJD(5) * t46 - t131 * t74 + t135 * t75;
t202 = t25 * mrSges(6,1) - t24 * mrSges(6,2);
t76 = -mrSges(6,2) * t130 + mrSges(6,3) * t91;
t14 = t24 * t76;
t77 = mrSges(6,1) * t130 - mrSges(6,3) * t92;
t15 = t25 * t77;
t191 = pkin(4) * t136;
t109 = (-pkin(3) - t191) * t129;
t36 = mrSges(6,1) * t60 + mrSges(6,2) * t59;
t30 = t109 * t36;
t192 = mrSges(6,3) * t60;
t31 = t46 * t192;
t170 = qJD(4) * t132;
t159 = t129 * t170;
t116 = pkin(4) * t159;
t61 = -mrSges(6,1) * t91 + mrSges(6,2) * t92;
t48 = t61 * t116;
t107 = -mrSges(5,2) * t130 + mrSges(5,3) * t178;
t95 = -pkin(9) * t159 + t117;
t69 = t95 * t107;
t179 = t129 * t132;
t106 = mrSges(5,1) * t130 - mrSges(5,3) * t179;
t96 = t103 * qJD(4);
t70 = t96 * t106;
t200 = t14 + t15 + t30 - t31 + t48 + t69 - t70;
t1 = t4 * t76;
t11 = t21 * t192;
t18 = -t170 * t89 + t164;
t12 = t18 * t107;
t19 = -qJD(4) * t50 + t154;
t13 = t19 * t106;
t2 = t5 * t77;
t71 = (-t97 - t191) * t129;
t26 = t71 * t36;
t51 = -t129 * t66 + t116;
t33 = t51 * t61;
t63 = t66 * mrSges(4,1);
t199 = t1 - t11 + t12 + t13 + t2 + t26 + t33 + t63;
t198 = 2 * m(5);
t197 = 2 * m(6);
t93 = (mrSges(5,1) * t132 + mrSges(5,2) * t136) * t129 * qJD(4);
t196 = -0.2e1 * t93;
t195 = m(6) * pkin(4);
t188 = t59 * mrSges(6,3);
t187 = t65 * mrSges(4,2);
t185 = Ifges(5,4) * t132;
t184 = Ifges(5,4) * t136;
t182 = t128 * t66;
t90 = Ifges(5,6) * t130 + (Ifges(5,2) * t136 + t185) * t129;
t181 = t132 * t90;
t101 = (-mrSges(5,1) * t136 + mrSges(5,2) * t132) * t129;
t180 = t66 * t101;
t168 = qJD(5) * t131;
t167 = qJD(5) * t135;
t166 = -0.2e1 * t188;
t165 = 0.2e1 * mrSges(5,3);
t151 = t129 * t163;
t146 = -Ifges(5,6) * t159 + t205;
t145 = -t135 * t188 - t168 * t77;
t144 = (-mrSges(3,1) * t134 - mrSges(3,2) * t138) * qJD(2) * pkin(1);
t143 = (-mrSges(4,1) * t133 - mrSges(4,2) * t137) * t183;
t142 = -0.2e1 * t91 * Ifges(6,2) * t60 + 0.2e1 * t92 * t59 * Ifges(6,1) + (Ifges(5,1) * t132 + t184) * t129 * t158 + 0.2e1 * (t59 * t91 - t60 * t92) * Ifges(6,4) + (t146 + 0.2e1 * t186 + t205) * t130 + ((Ifges(5,1) * t136 - t185) * t170 + (-Ifges(5,2) * t132 + t184) * t169) * t128;
t140 = t146 + t186 + (-t131 * t192 + t167 * t76) * pkin(4);
t22 = t38 * t192;
t98 = (-t124 - t191) * t129;
t29 = t98 * t36;
t52 = (-qJD(4) * t112 - t152) * t132 + t173;
t39 = t52 * t107;
t53 = -qJD(4) * t79 + t141;
t40 = t53 * t106;
t94 = t116 + t151;
t42 = t94 * t61;
t6 = t9 * t76;
t7 = t10 * t77;
t84 = t101 * t151;
t139 = t29 + t142 - t22 + t7 + t6 + t84 + t42 + t40 + t39;
t104 = (-mrSges(6,1) * t131 - mrSges(6,2) * t135) * qJD(5) * pkin(4);
t102 = -pkin(9) * t179 + t122;
t78 = -t112 * t132 + t111;
t49 = -t132 * t89 + t86;
t3 = [-0.2e1 * t187 + (-0.2e1 * t180 + t97 * t196 + (-t181 + (-t132 * t50 - t136 * t49) * t165) * qJD(4)) * t129 + t20 * t166 + 0.2e1 * t144 + 0.2e1 * t1 + 0.2e1 * t33 + 0.2e1 * t26 + t142 + 0.2e1 * t2 + 0.2e1 * t13 + 0.2e1 * t12 - 0.2e1 * t11 + 0.2e1 * t63 + (t18 * t50 + t182 * t97 + t19 * t49) * t198 + 0.2e1 * m(4) * (t100 * t65 + t66 * t99) + (t20 * t5 + t21 * t4 + t51 * t71) * t197; m(4) * (t100 * t171 + t133 * t65 + t137 * t66 - t172 * t99) * pkin(2) + m(6) * (t10 * t20 + t21 * t9 + t37 * t5 + t38 * t4 + t51 * t98 + t71 * t94) + (-t180 + (-t124 - t97) * t93 + (-t181 + ((-t49 - t78) * t136 + (-t50 - t79) * t132) * mrSges(5,3)) * qJD(4)) * t129 + (-t65 - t162) * mrSges(4,2) - mrSges(4,1) * t163 + (-t20 - t37) * t188 + t139 + t144 + (t124 * t182 - t153 * t97 + t79 * t18 + t78 * t19 + t53 * t49 + t52 * t50) * m(5) + t199; (t124 * t196 + (-t181 + (-t132 * t79 - t136 * t78) * t165) * qJD(4)) * t129 + t37 * t166 + 0.2e1 * t143 + 0.2e1 * t29 + t142 - 0.2e1 * t22 + 0.2e1 * t7 + 0.2e1 * t6 + 0.2e1 * t84 + 0.2e1 * t42 + (-t124 * t153 + t52 * t79 + t53 * t78) * t198 + (t10 * t37 + t38 * t9 + t94 * t98) * t197 + 0.2e1 * t40 + 0.2e1 * t39; (-t180 + (-pkin(3) - t97) * t93 + ((t71 * t195 - t90) * t132 + ((-t102 - t49) * t136 + (-t103 - t50) * t132) * mrSges(5,3)) * qJD(4)) * t129 + (-t20 - t45) * t188 + t142 - t187 + m(6) * (t109 * t51 + t20 * t25 + t21 * t24 + t4 * t46 + t45 * t5) + m(5) * (pkin(3) * t182 + t102 * t19 + t103 * t18 - t49 * t96 + t50 * t95) + t199 + t200; ((-pkin(3) - t124) * t93 + ((t98 * t195 - t90) * t132 + ((-t102 - t78) * t136 + (-t103 - t79) * t132) * mrSges(5,3)) * qJD(4)) * t129 + t139 + m(6) * (t10 * t45 + t109 * t94 + t24 * t38 + t25 * t37 + t46 * t9) + m(5) * (-pkin(3) * t153 + t102 * t53 + t103 * t52 - t78 * t96 + t79 * t95) + (-t37 - t45) * t188 + t143 + t200; (t109 * t116 + t24 * t46 + t25 * t45) * t197 + (-t102 * t96 + t103 * t95) * t198 + (pkin(3) * t196 + (-t181 + (-t102 * t136 - t103 * t132) * t165) * qJD(4)) * t129 + t45 * t166 - 0.2e1 * t31 + 0.2e1 * t30 + t142 + 0.2e1 * t15 + 0.2e1 * t14 + 0.2e1 * t69 - 0.2e1 * t70 + 0.2e1 * t48; t19 * mrSges(5,1) - t18 * mrSges(5,2) + (m(6) * (t131 * t4 + t135 * t5 + t167 * t21 - t168 * t20) + t145) * pkin(4) + t140 + t204; t53 * mrSges(5,1) - t52 * mrSges(5,2) + (m(6) * (t10 * t135 + t131 * t9 + t167 * t38 - t168 * t37) + t145) * pkin(4) + t140 + t203; -t96 * mrSges(5,1) - t95 * mrSges(5,2) + (m(6) * (t131 * t24 + t135 * t25 + t167 * t46 - t168 * t45) + t145) * pkin(4) + t140 + t202; 0.2e1 * t104; t186 + t204; t186 + t203; t186 + t202; t104; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
