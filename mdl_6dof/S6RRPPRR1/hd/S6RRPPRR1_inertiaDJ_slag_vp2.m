% Calculate time derivative of joint inertia matrix for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:54
% EndTime: 2019-03-09 08:45:59
% DurationCPUTime: 2.81s
% Computational Cost: add. (4575->311), mult. (9619->441), div. (0->0), fcn. (9414->8), ass. (0->148)
t142 = (mrSges(5,2) + mrSges(4,3));
t199 = 2 * t142;
t90 = sin(qJ(6));
t93 = cos(qJ(6));
t139 = t90 ^ 2 + t93 ^ 2;
t121 = t139 * mrSges(7,3);
t94 = cos(qJ(5));
t198 = (mrSges(6,2) - t121) * t94;
t71 = -t93 * mrSges(7,1) + mrSges(7,2) * t90;
t191 = t71 - mrSges(6,1);
t183 = -m(7) * pkin(5) + t191;
t136 = cos(pkin(10));
t89 = sin(pkin(10));
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t102 = t136 * t95 - t89 * t92;
t61 = t136 * t92 + t89 * t95;
t91 = sin(qJ(5));
t108 = -t102 * t94 - t61 * t91;
t154 = t90 * mrSges(7,3);
t37 = -t102 * t91 + t61 * t94;
t24 = mrSges(7,2) * t108 - t154 * t37;
t149 = t93 * mrSges(7,3);
t25 = -mrSges(7,1) * t108 - t149 * t37;
t109 = -t93 * t24 + t90 * t25;
t81 = -t95 * pkin(2) - pkin(1);
t35 = -pkin(3) * t102 - t61 * qJ(4) + t81;
t30 = pkin(4) * t102 - t35;
t13 = -pkin(5) * t108 - pkin(9) * t37 + t30;
t141 = -qJ(3) - pkin(7);
t72 = t141 * t92;
t73 = t141 * t95;
t40 = -t136 * t72 - t73 * t89;
t105 = -pkin(8) * t61 + t40;
t41 = -t136 * t73 + t89 * t72;
t34 = -pkin(8) * t102 + t41;
t17 = t105 * t91 + t94 * t34;
t7 = t13 * t93 - t17 * t90;
t8 = t13 * t90 + t17 * t93;
t187 = -t7 * t90 + t8 * t93;
t197 = m(7) * t187 + t108 * mrSges(6,3) - t109;
t193 = m(7) * pkin(9);
t119 = t136 * pkin(2);
t80 = -t119 - pkin(3);
t77 = -pkin(4) + t80;
t172 = pkin(2) * t89;
t78 = qJ(4) + t172;
t44 = t77 * t94 - t78 * t91;
t38 = qJD(4) * t94 + qJD(5) * t44;
t192 = t38 * mrSges(6,2);
t45 = t91 * t77 + t94 * t78;
t128 = qJD(6) * t93;
t129 = qJD(6) * t90;
t190 = t7 * t128 + t8 * t129;
t123 = t37 * t129;
t55 = t61 * qJD(2);
t56 = t102 * qJD(2);
t21 = qJD(5) * t108 + t55 * t91 + t56 * t94;
t147 = t93 * t21;
t103 = t123 - t147;
t131 = qJD(5) * t37;
t22 = -t94 * t55 + t56 * t91 + t131;
t11 = mrSges(7,1) * t22 + mrSges(7,3) * t103;
t122 = t37 * t128;
t152 = t90 * t21;
t104 = t122 + t152;
t12 = -mrSges(7,2) * t22 - mrSges(7,3) * t104;
t189 = -t90 * t11 + t93 * t12;
t188 = t139 * t94;
t74 = Ifges(7,5) * t90 + Ifges(7,6) * t93;
t186 = t74 / 0.2e1 - Ifges(6,6);
t4 = -Ifges(7,1) * t103 - Ifges(7,4) * t104 + Ifges(7,5) * t22;
t67 = Ifges(7,4) * t128 - Ifges(7,2) * t129;
t185 = t37 * t67 / 0.2e1 - t4 / 0.2e1;
t184 = m(5) * t78 + mrSges(5,3);
t168 = Ifges(7,4) * t93;
t76 = Ifges(7,1) * t90 + t168;
t145 = t93 * t76;
t169 = Ifges(7,4) * t90;
t75 = Ifges(7,2) * t93 + t169;
t150 = t90 * t75;
t68 = Ifges(7,1) * t128 - Ifges(7,4) * t129;
t97 = -(-t145 + t150) * qJD(6) + t93 * t67 + t90 * t68;
t182 = 2 * m(6);
t181 = 0.2e1 * m(7);
t132 = qJD(5) * t17;
t115 = qJD(2) * t141;
t100 = qJD(3) * t95 + t115 * t92;
t99 = -t92 * qJD(3) + t115 * t95;
t33 = t136 * t100 + t89 * t99;
t28 = pkin(8) * t55 + t33;
t32 = t100 * t89 - t136 * t99;
t96 = -t56 * pkin(8) + t32;
t10 = t91 * t28 - t94 * t96 + t132;
t180 = 0.2e1 * t10;
t135 = qJD(2) * t92;
t126 = pkin(2) * t135;
t29 = t55 * pkin(3) - t56 * qJ(4) - t61 * qJD(4) + t126;
t23 = pkin(4) * t55 + t29;
t179 = -0.2e1 * t23;
t110 = mrSges(7,1) * t90 + mrSges(7,2) * t93;
t65 = t110 * qJD(6);
t178 = -0.2e1 * t65;
t177 = 0.2e1 * t81;
t176 = m(4) * pkin(2);
t175 = t21 / 0.2e1;
t174 = -t75 / 0.2e1;
t167 = Ifges(7,5) * t93;
t16 = -t105 * t94 + t91 * t34;
t166 = t10 * t16;
t39 = qJD(4) * t91 + t45 * qJD(5);
t165 = t16 * t39;
t164 = t21 * mrSges(6,3);
t161 = t37 * t68;
t160 = t37 * t76;
t159 = t38 * t91;
t158 = t39 * t94;
t157 = t56 * mrSges(5,2);
t144 = t94 * t65;
t140 = Ifges(7,5) * t147 + Ifges(7,3) * t22;
t138 = t7 * qJD(6);
t137 = t8 * qJD(6);
t134 = qJD(2) * t95;
t133 = qJD(5) * t16;
t130 = qJD(5) * t94;
t127 = 0.2e1 * t95;
t120 = t139 * t38;
t117 = -Ifges(7,6) * t90 - (2 * Ifges(6,4));
t116 = t32 * t40 + t41 * t33;
t114 = -2 * Ifges(4,4) + 2 * Ifges(5,5);
t6 = pkin(5) * t22 - pkin(9) * t21 - t23;
t9 = t94 * t28 + t91 * t96 - t133;
t1 = t6 * t90 + t9 * t93 + t138;
t2 = t6 * t93 - t9 * t90 - t137;
t113 = t1 * t93 - t2 * t90;
t111 = t22 * mrSges(6,1) + t21 * mrSges(6,2);
t82 = Ifges(7,6) * t129;
t98 = Ifges(6,5) * t21 + t16 * t65 - t108 * (Ifges(7,5) * t128 - t82) / 0.2e1 - t9 * mrSges(6,2);
t52 = t56 * mrSges(4,2);
t51 = t55 * mrSges(5,1);
t43 = -pkin(9) + t45;
t42 = pkin(5) - t44;
t19 = t110 * t37;
t15 = -Ifges(7,5) * t108 + (Ifges(7,1) * t93 - t169) * t37;
t14 = -Ifges(7,6) * t108 + (-Ifges(7,2) * t90 + t168) * t37;
t5 = mrSges(7,1) * t104 - mrSges(7,2) * t103;
t3 = -Ifges(7,4) * t103 - Ifges(7,2) * t104 + Ifges(7,6) * t22;
t18 = [0.2e1 * t29 * (-mrSges(5,1) * t102 - mrSges(5,3) * t61) + 0.2e1 * t35 * t51 + 0.2e1 * t30 * t111 + 0.2e1 * t1 * t24 + 0.2e1 * t2 * t25 + 0.2e1 * t16 * t5 + t19 * t180 + 0.2e1 * t7 * t11 + 0.2e1 * t8 * t12 - t14 * t152 + t15 * t147 + t52 * t177 + 0.2e1 * (t16 * t21 - t17 * t22) * mrSges(6,3) + (-0.2e1 * t35 * mrSges(5,3) - t102 * t114 + 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t61 + t40 * t199) * t56 + (mrSges(4,1) * t177 + t61 * t114 - 0.2e1 * (Ifges(4,2) + Ifges(5,3)) * t102 - 0.2e1 * t142 * t41) * t55 + 0.2e1 * m(4) * t116 + (t1 * t8 + t2 * t7 + t166) * t181 + (t17 * t9 - t23 * t30 + t166) * t182 + 0.2e1 * m(5) * (t29 * t35 + t116) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t95) * t127 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t102 + mrSges(4,2) * t61) - 0.2e1 * pkin(1) * mrSges(3,1) + t176 * t177 - 0.2e1 * Ifges(3,4) * t92 + (Ifges(3,1) - Ifges(3,2)) * t127) * t92) * qJD(2) - (mrSges(6,1) * t179 - 0.2e1 * mrSges(6,3) * t9 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t22 + t117 * t21 + t140) * t108 + (mrSges(6,2) * t179 + mrSges(6,3) * t180 + 0.2e1 * Ifges(6,1) * t21 - t90 * t3 + t93 * t4 + (t117 + t167) * t22 + (t108 * t74 - t93 * t14 - t90 * t15) * qJD(6)) * t37 + (t102 * t33 + t32 * t61) * t199; (-mrSges(3,1) * t134 + mrSges(3,2) * t135) * pkin(7) - t1 * t149 - t21 * t145 / 0.2e1 - Ifges(3,6) * t135 + t150 * t175 + t80 * t157 - t98 + (t37 * mrSges(6,3) + t19) * t39 - t15 * t128 / 0.2e1 - t44 * t164 + (-mrSges(5,2) * t78 - mrSges(4,3) * t172 - Ifges(4,6) + Ifges(5,6)) * t55 + t75 * t122 / 0.2e1 + (-m(6) * t44 + m(7) * t42 - t191) * t10 + (m(7) * ((-t7 * t93 - t8 * t90) * qJD(6) + t113) - t25 * t128 - t24 * t129 + t189) * t43 + t190 * mrSges(7,3) + (-t45 * mrSges(6,3) - t186) * t22 + (t176 * t89 - mrSges(4,2) + t184) * t33 + t185 * t90 + Ifges(3,5) * t134 + t2 * t154 + (-mrSges(4,3) * t119 + Ifges(5,4) + Ifges(4,5)) * t56 + (m(5) * t41 + mrSges(5,2) * t102) * qJD(4) + m(7) * t165 + m(6) * (t45 * t9 + t165) + (m(5) * t80 - t136 * t176 - mrSges(4,1) - mrSges(5,1)) * t32 + t42 * t5 + (m(6) * t17 + t197) * t38 + (t14 + t160) * t129 / 0.2e1 - (t3 + t161) * t93 / 0.2e1; t42 * t178 + 0.2e1 * t192 + (t120 * t43 + t39 * t42) * t181 + (t38 * t45 - t39 * t44) * t182 + t97 - 0.2e1 * t191 * t39 + 0.2e1 * t184 * qJD(4) - 0.2e1 * t121 * t38; m(4) * t126 + t55 * mrSges(4,1) - t56 * mrSges(5,3) - t93 * t11 - t90 * t12 + t51 + t52 + t109 * qJD(6) + m(7) * (-qJD(6) * t187 - t1 * t90 - t2 * t93) + m(6) * t23 + m(5) * t29 - t111; 0; 0; t157 + m(5) * t32 + (-t164 - t5 - m(7) * t10 + m(6) * (-t10 + t132) + t197 * qJD(5)) * t94 + (qJD(5) * t19 + (-t90 * t24 - t93 * t25) * qJD(6) + (-t22 + t131) * mrSges(6,3) + m(7) * (t113 + t133 - t190) + m(6) * (t9 + t133) + t189) * t91; t144 + m(7) * (t139 * t159 - t158) + m(6) * (-t158 + t159) + (-t191 * t91 + t198 + m(7) * (t188 * t43 + t42 * t91) + m(6) * (-t44 * t91 + t45 * t94)) * qJD(5); 0; (-0.1e1 + t139) * t91 * t130 * t181; -pkin(5) * t5 + t186 * t22 + t183 * t10 + (t3 / 0.2e1 + t1 * mrSges(7,3) + t161 / 0.2e1 + t76 * t175 + (t37 * t174 - t7 * mrSges(7,3) + t15 / 0.2e1) * qJD(6) + (m(7) * (t1 - t138) + t12 - qJD(6) * t25) * pkin(9)) * t93 + (-t2 * mrSges(7,3) + t21 * t174 + (-t160 / 0.2e1 - t8 * mrSges(7,3) - t14 / 0.2e1) * qJD(6) + (m(7) * (-t2 - t137) - qJD(6) * t24 - t11) * pkin(9) - t185) * t90 + t98; -t192 + (t42 + pkin(5)) * t65 + (mrSges(7,3) + t193) * t120 + t183 * t39 - t97; 0; -t144 + (t183 * t91 + t188 * t193 - t198) * qJD(5); pkin(5) * t178 + t97; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,5) * t123 - Ifges(7,6) * t104 + t140; t82 - t110 * t38 + (t43 * t71 - t167) * qJD(6); t65; (t129 * t91 - t130 * t93) * mrSges(7,2) + (-t128 * t91 - t130 * t90) * mrSges(7,1); -t82 + (pkin(9) * t71 + t167) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
