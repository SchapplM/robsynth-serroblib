% Calculate joint inertia matrix for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:09:23
% EndTime: 2018-11-23 16:09:25
% DurationCPUTime: 1.64s
% Computational Cost: add. (5702->352), mult. (14985->527), div. (0->0), fcn. (17387->14), ass. (0->143)
t139 = sin(qJ(5));
t174 = cos(qJ(5));
t130 = sin(pkin(13));
t134 = cos(pkin(13));
t131 = sin(pkin(12));
t133 = sin(pkin(6));
t135 = cos(pkin(12));
t137 = cos(pkin(6));
t142 = cos(qJ(3));
t136 = cos(pkin(7));
t140 = sin(qJ(3));
t159 = t136 * t140;
t132 = sin(pkin(7));
t162 = t132 * t140;
t79 = t137 * t162 + (t131 * t142 + t135 * t159) * t133;
t160 = t133 * t135;
t93 = -t132 * t160 + t136 * t137;
t55 = -t130 * t79 + t134 * t93;
t56 = t130 * t93 + t134 * t79;
t42 = t139 * t56 - t174 * t55;
t43 = t139 * t55 + t174 * t56;
t22 = t42 * mrSges(6,1) + t43 * mrSges(6,2);
t44 = -t55 * mrSges(5,1) + t56 * mrSges(5,2);
t191 = -t22 - t44;
t106 = -t134 * mrSges(5,1) + t130 * mrSges(5,2);
t103 = t130 * t139 - t134 * t174;
t104 = t130 * t174 + t139 * t134;
t81 = t103 * mrSges(6,1) + t104 * mrSges(6,2);
t190 = -t106 - t81;
t189 = m(7) * pkin(11) + mrSges(7,3);
t153 = t136 * t160;
t173 = pkin(1) * t137;
t97 = qJ(2) * t160 + t131 * t173;
t72 = (t132 * t137 + t153) * pkin(9) + t97;
t118 = t135 * t173;
t163 = t131 * t133;
t80 = pkin(2) * t137 + t118 + (-pkin(9) * t136 - qJ(2)) * t163;
t89 = (-pkin(9) * t131 * t132 - pkin(2) * t135 - pkin(1)) * t133;
t45 = -t140 * t72 + (t132 * t89 + t136 * t80) * t142;
t138 = sin(qJ(6));
t141 = cos(qJ(6));
t110 = -mrSges(7,1) * t141 + mrSges(7,2) * t138;
t188 = -m(7) * pkin(5) - mrSges(6,1) + t110;
t94 = -t130 * t162 + t134 * t136;
t95 = t130 * t136 + t134 * t162;
t64 = t139 * t95 - t174 * t94;
t187 = t64 ^ 2;
t169 = pkin(10) + qJ(4);
t107 = t169 * t134;
t152 = t169 * t130;
t85 = t107 * t139 + t152 * t174;
t186 = t85 ^ 2;
t185 = 2 * mrSges(3,1);
t184 = 0.2e1 * t85;
t183 = 0.2e1 * t137;
t161 = t132 * t142;
t78 = -t137 * t161 + t140 * t163 - t142 * t153;
t26 = -t138 * t43 + t141 * t78;
t182 = t26 / 0.2e1;
t27 = t138 * t78 + t141 * t43;
t181 = t27 / 0.2e1;
t179 = -m(5) - m(6);
t111 = Ifges(7,5) * t138 + Ifges(7,6) * t141;
t178 = t111 / 0.2e1;
t177 = -t138 / 0.2e1;
t176 = t138 / 0.2e1;
t175 = t141 / 0.2e1;
t172 = pkin(11) * t138;
t171 = pkin(11) * t141;
t170 = t64 * t85;
t13 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t29 = mrSges(6,1) * t78 - mrSges(6,3) * t43;
t168 = t13 - t29;
t52 = -t132 * t80 + t136 * t89;
t34 = pkin(3) * t78 - qJ(4) * t79 + t52;
t46 = t142 * t72 + t80 * t159 + t89 * t162;
t36 = qJ(4) * t93 + t46;
t20 = -t130 * t36 + t134 * t34;
t12 = pkin(4) * t78 - pkin(10) * t56 + t20;
t21 = t130 * t34 + t134 * t36;
t15 = pkin(10) * t55 + t21;
t6 = t139 * t12 + t174 * t15;
t167 = Ifges(7,4) * t138;
t166 = Ifges(7,4) * t141;
t165 = t104 * t138;
t164 = t104 * t141;
t158 = Ifges(6,5) * t104 - Ifges(6,6) * t103;
t157 = t130 ^ 2 + t134 ^ 2;
t156 = t138 ^ 2 + t141 ^ 2;
t7 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t42;
t155 = Ifges(6,5) * t43 - Ifges(6,6) * t42 + Ifges(6,3) * t78;
t154 = Ifges(4,5) * t79 - Ifges(4,6) * t78 + Ifges(4,3) * t93;
t121 = -pkin(4) * t134 - pkin(3);
t37 = -pkin(3) * t93 - t45;
t23 = -pkin(4) * t55 + t37;
t10 = pkin(5) * t42 - pkin(11) * t43 + t23;
t4 = pkin(11) * t78 + t6;
t1 = t10 * t141 - t138 * t4;
t2 = t10 * t138 + t141 * t4;
t151 = -t1 * t138 + t141 * t2;
t150 = mrSges(7,1) * t138 + mrSges(7,2) * t141;
t149 = -t20 * t130 + t21 * t134;
t148 = -t130 * t94 + t134 * t95;
t61 = Ifges(7,5) * t164 - Ifges(7,6) * t165 + Ifges(7,3) * t103;
t5 = t12 * t174 - t139 * t15;
t126 = t132 ^ 2;
t119 = t126 * t142 ^ 2;
t116 = mrSges(3,2) * t163;
t113 = Ifges(7,1) * t138 + t166;
t112 = Ifges(7,2) * t141 + t167;
t109 = Ifges(5,1) * t130 + Ifges(5,4) * t134;
t108 = Ifges(5,4) * t130 + Ifges(5,2) * t134;
t96 = -qJ(2) * t163 + t118;
t87 = t107 * t174 - t139 * t152;
t83 = Ifges(6,1) * t104 - Ifges(6,4) * t103;
t82 = Ifges(6,4) * t104 - Ifges(6,2) * t103;
t74 = mrSges(7,1) * t103 - mrSges(7,3) * t164;
t73 = -mrSges(7,2) * t103 - mrSges(7,3) * t165;
t71 = pkin(5) * t103 - pkin(11) * t104 + t121;
t70 = t150 * t104;
t66 = t139 * t94 + t174 * t95;
t63 = Ifges(7,5) * t103 + (Ifges(7,1) * t141 - t167) * t104;
t62 = Ifges(7,6) * t103 + (-Ifges(7,2) * t138 + t166) * t104;
t60 = mrSges(4,1) * t93 - mrSges(4,3) * t79;
t59 = -mrSges(4,2) * t93 - mrSges(4,3) * t78;
t58 = -t138 * t161 + t141 * t66;
t57 = -t138 * t66 - t141 * t161;
t51 = mrSges(4,1) * t78 + mrSges(4,2) * t79;
t50 = t138 * t71 + t141 * t87;
t49 = -t138 * t87 + t141 * t71;
t48 = mrSges(5,1) * t78 - mrSges(5,3) * t56;
t47 = -mrSges(5,2) * t78 + mrSges(5,3) * t55;
t31 = Ifges(5,1) * t56 + Ifges(5,4) * t55 + t78 * Ifges(5,5);
t30 = Ifges(5,4) * t56 + Ifges(5,2) * t55 + t78 * Ifges(5,6);
t28 = -mrSges(6,2) * t78 - mrSges(6,3) * t42;
t19 = Ifges(6,1) * t43 - Ifges(6,4) * t42 + Ifges(6,5) * t78;
t18 = Ifges(6,4) * t43 - Ifges(6,2) * t42 + Ifges(6,6) * t78;
t17 = mrSges(7,1) * t42 - mrSges(7,3) * t27;
t16 = -mrSges(7,2) * t42 + mrSges(7,3) * t26;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t42;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t42;
t3 = -t78 * pkin(5) - t5;
t11 = [(t7 - t18) * t42 + (-0.2e1 * t97 * mrSges(3,2) + Ifges(3,3) * t137 + t185 * t96) * t137 + (-0.2e1 * pkin(1) * t116 + (-0.2e1 * t96 * mrSges(3,3) + Ifges(3,1) * t163 + Ifges(3,5) * t183) * t131 + (0.2e1 * t97 * mrSges(3,3) + Ifges(3,6) * t183 + (0.2e1 * Ifges(3,4) * t131 + Ifges(3,2) * t135 + pkin(1) * t185) * t133) * t135) * t133 + t93 * t154 + (-0.2e1 * Ifges(4,4) * t79 + Ifges(5,5) * t56 - Ifges(4,6) * t93 + Ifges(5,6) * t55 + (Ifges(5,3) + Ifges(4,2)) * t78 + t155) * t78 + Ifges(2,3) + m(4) * (t45 ^ 2 + t46 ^ 2 + t52 ^ 2) + m(6) * (t23 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2 + t37 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(3) * (pkin(1) ^ 2 * t133 ^ 2 + t96 ^ 2 + t97 ^ 2) + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t16 + 0.2e1 * t1 * t17 + 0.2e1 * t23 * t22 + t26 * t8 + t27 * t9 + 0.2e1 * t6 * t28 + 0.2e1 * t5 * t29 + t43 * t19 + 0.2e1 * t37 * t44 + 0.2e1 * t21 * t47 + 0.2e1 * t20 * t48 + 0.2e1 * t52 * t51 + t55 * t30 + t56 * t31 + 0.2e1 * t46 * t59 + 0.2e1 * t45 * t60 + t79 * (Ifges(4,1) * t79 + Ifges(4,5) * t93); t136 * t51 + t58 * t16 + t57 * t17 + t66 * t28 + t95 * t47 + t94 * t48 + t116 + t168 * t64 + (-m(3) * pkin(1) - mrSges(3,1) * t135) * t133 + (t140 * t59 + (t60 + t191) * t142) * t132 + m(7) * (t1 * t57 + t2 * t58 + t3 * t64) + m(6) * (-t161 * t23 - t5 * t64 + t6 * t66) + m(5) * (-t161 * t37 + t20 * t94 + t21 * t95) + m(4) * (t136 * t52 + (t140 * t46 + t142 * t45) * t132); m(3) + m(7) * (t57 ^ 2 + t58 ^ 2 + t187) + m(6) * (t66 ^ 2 + t119 + t187) + m(5) * (t94 ^ 2 + t95 ^ 2 + t119) + m(4) * (t126 * t140 ^ 2 + t136 ^ 2 + t119); t168 * t85 + m(7) * (t1 * t49 + t2 * t50 + t3 * t85) + m(6) * (t121 * t23 - t5 * t85 + t6 * t87) + t63 * t181 + t62 * t182 + (t61 / 0.2e1 - t82 / 0.2e1) * t42 + (Ifges(5,5) * t130 + Ifges(5,6) * t134 + t158) * t78 / 0.2e1 + (-t130 * t48 + t134 * t47) * qJ(4) + t154 + (-t6 * mrSges(6,3) + t7 / 0.2e1 - t18 / 0.2e1) * t103 + (-t5 * mrSges(6,3) + t8 * t177 + t9 * t175 + t19 / 0.2e1) * t104 + m(5) * (-pkin(3) * t37 + qJ(4) * t149) + t149 * mrSges(5,3) - pkin(3) * t44 + t45 * mrSges(4,1) - t46 * mrSges(4,2) + t49 * t17 + t50 * t16 + t3 * t70 + t2 * t73 + t1 * t74 + t23 * t81 + t43 * t83 / 0.2e1 + t87 * t28 + t37 * t106 + t55 * t108 / 0.2e1 + t56 * t109 / 0.2e1 + t121 * t22 + t130 * t31 / 0.2e1 + t134 * t30 / 0.2e1; t57 * t74 + t58 * t73 + t64 * t70 + (-t103 * t66 + t104 * t64) * mrSges(6,3) + t148 * mrSges(5,3) + (-t140 * mrSges(4,2) + (mrSges(4,1) + t190) * t142) * t132 + m(7) * (t49 * t57 + t50 * t58 + t170) + m(6) * (-t121 * t161 + t66 * t87 + t170) + m(5) * (pkin(3) * t161 + qJ(4) * t148); -0.2e1 * pkin(3) * t106 + t134 * t108 + t130 * t109 + 0.2e1 * t121 * t81 + 0.2e1 * t49 * t74 + 0.2e1 * t50 * t73 + t70 * t184 + Ifges(4,3) + 0.2e1 * t157 * qJ(4) * mrSges(5,3) + (-0.2e1 * t87 * mrSges(6,3) + t61 - t82) * t103 + m(7) * (t49 ^ 2 + t50 ^ 2 + t186) + m(6) * (t121 ^ 2 + t87 ^ 2 + t186) + m(5) * (qJ(4) ^ 2 * t157 + pkin(3) ^ 2) + (mrSges(6,3) * t184 - t138 * t62 + t141 * t63 + t83) * t104; t138 * t16 + t141 * t17 + m(7) * (t1 * t141 + t138 * t2) + m(6) * t23 + m(5) * t37 - t191; m(7) * (t138 * t58 + t141 * t57) + t179 * t161; -m(5) * pkin(3) + t138 * t73 + t141 * t74 + m(7) * (t138 * t50 + t141 * t49) + m(6) * t121 - t190; m(7) * t156 - t179; t113 * t181 + m(7) * (-pkin(5) * t3 + pkin(11) * t151) + t16 * t171 - t17 * t172 + t9 * t176 + t8 * t175 - pkin(5) * t13 + t3 * t110 + t112 * t182 + t42 * t178 - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t151 * mrSges(7,3) + t155; -t66 * mrSges(6,2) + t189 * (-t138 * t57 + t141 * t58) + t188 * t64; t73 * t171 - pkin(5) * t70 - t74 * t172 + t63 * t176 + t62 * t175 + t103 * t178 - t87 * mrSges(6,2) + (t112 * t177 + t113 * t175) * t104 + t158 + t189 * (-t138 * t49 + t141 * t50) + t188 * t85; 0; Ifges(6,3) - 0.2e1 * pkin(5) * t110 + t138 * t113 + t141 * t112 + m(7) * (pkin(11) ^ 2 * t156 + pkin(5) ^ 2) + 0.2e1 * t156 * pkin(11) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t57 - mrSges(7,2) * t58; mrSges(7,1) * t49 - mrSges(7,2) * t50 + t61; -t110; -pkin(11) * t150 + t111; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
