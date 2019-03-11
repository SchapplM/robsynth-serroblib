% Calculate joint inertia matrix for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:40
% EndTime: 2019-03-09 10:55:44
% DurationCPUTime: 1.72s
% Computational Cost: add. (4272->391), mult. (9553->569), div. (0->0), fcn. (10979->12), ass. (0->143)
t153 = sin(pkin(6));
t161 = cos(qJ(2));
t171 = t153 * t161;
t155 = cos(pkin(11));
t181 = pkin(9) + qJ(3);
t129 = t181 * t155;
t158 = sin(qJ(4));
t152 = sin(pkin(11));
t165 = t181 * t152;
t183 = cos(qJ(4));
t94 = t129 * t158 + t183 * t165;
t194 = t94 ^ 2;
t193 = 0.2e1 * t94;
t151 = sin(pkin(12));
t154 = cos(pkin(12));
t156 = cos(pkin(6));
t159 = sin(qJ(2));
t172 = t153 * t159;
t110 = -t152 * t172 + t155 * t156;
t111 = t152 * t156 + t155 * t172;
t74 = t158 * t110 + t111 * t183;
t54 = -t151 * t74 - t154 * t171;
t55 = -t151 * t171 + t154 * t74;
t73 = -t110 * t183 + t111 * t158;
t19 = Ifges(6,1) * t55 + Ifges(6,4) * t54 + Ifges(6,5) * t73;
t192 = t19 / 0.2e1;
t122 = t152 * t158 - t155 * t183;
t124 = t152 * t183 + t158 * t155;
t177 = Ifges(6,4) * t154;
t57 = Ifges(6,6) * t122 + (-Ifges(6,2) * t151 + t177) * t124;
t191 = t57 / 0.2e1;
t178 = Ifges(6,4) * t151;
t58 = Ifges(6,5) * t122 + (Ifges(6,1) * t154 - t178) * t124;
t190 = t58 / 0.2e1;
t157 = sin(qJ(6));
t160 = cos(qJ(6));
t121 = -t151 * t157 + t154 * t160;
t123 = t151 * t160 + t154 * t157;
t89 = Ifges(7,4) * t123 + Ifges(7,2) * t121;
t189 = t89 / 0.2e1;
t91 = Ifges(7,1) * t123 + Ifges(7,4) * t121;
t188 = t91 / 0.2e1;
t187 = t121 / 0.2e1;
t186 = t123 / 0.2e1;
t133 = Ifges(6,1) * t151 + t177;
t185 = t133 / 0.2e1;
t182 = pkin(1) * t161;
t180 = pkin(10) + qJ(5);
t113 = t156 * t159 * pkin(1) + pkin(8) * t171;
t103 = qJ(3) * t156 + t113;
t104 = (-pkin(2) * t161 - qJ(3) * t159 - pkin(1)) * t153;
t61 = -t103 * t152 + t155 * t104;
t45 = -pkin(3) * t171 - pkin(9) * t111 + t61;
t62 = t155 * t103 + t152 * t104;
t50 = pkin(9) * t110 + t62;
t23 = t158 * t45 + t183 * t50;
t20 = -qJ(5) * t171 + t23;
t137 = pkin(8) * t172;
t106 = t137 + (-pkin(2) - t182) * t156;
t79 = -pkin(3) * t110 + t106;
t26 = pkin(4) * t73 - qJ(5) * t74 + t79;
t6 = t151 * t26 + t154 * t20;
t179 = -Ifges(5,5) * t74 + Ifges(5,6) * t73;
t143 = -pkin(3) * t155 - pkin(2);
t82 = pkin(4) * t122 - qJ(5) * t124 + t143;
t97 = t129 * t183 - t158 * t165;
t47 = t151 * t82 + t154 * t97;
t112 = t156 * t182 - t137;
t176 = t112 * mrSges(3,1);
t175 = t113 * mrSges(3,2);
t174 = t124 * t151;
t173 = t124 * t154;
t81 = mrSges(6,1) * t174 + mrSges(6,2) * t173;
t88 = Ifges(7,5) * t123 + Ifges(7,6) * t121;
t170 = Ifges(5,5) * t124 - Ifges(5,6) * t122;
t169 = t151 ^ 2 + t154 ^ 2;
t168 = t152 ^ 2 + t155 ^ 2;
t30 = -t157 * t55 + t160 * t54;
t31 = t157 * t54 + t160 * t55;
t7 = Ifges(7,5) * t31 + Ifges(7,6) * t30 + Ifges(7,3) * t73;
t75 = t123 * t124;
t76 = t121 * t124;
t34 = Ifges(7,5) * t76 - Ifges(7,6) * t75 + Ifges(7,3) * t122;
t167 = t88 / 0.2e1 + Ifges(6,5) * t151 / 0.2e1 + Ifges(6,6) * t154 / 0.2e1;
t166 = Ifges(3,5) * t172 + Ifges(3,6) * t171 + Ifges(3,3) * t156;
t42 = t73 * mrSges(5,1) + t74 * mrSges(5,2);
t32 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t10 = -t30 * mrSges(7,1) + t31 * mrSges(7,2);
t43 = t75 * mrSges(7,1) + t76 * mrSges(7,2);
t5 = -t151 * t20 + t154 * t26;
t46 = -t151 * t97 + t154 * t82;
t80 = -t110 * mrSges(4,1) + t111 * mrSges(4,2);
t127 = -t155 * mrSges(4,1) + t152 * mrSges(4,2);
t87 = t122 * mrSges(5,1) + t124 * mrSges(5,2);
t126 = -t154 * mrSges(6,1) + t151 * mrSges(6,2);
t22 = -t158 * t50 + t183 * t45;
t86 = -t121 * mrSges(7,1) + t123 * mrSges(7,2);
t164 = -t61 * t152 + t62 * t155;
t21 = pkin(4) * t171 - t22;
t142 = -pkin(5) * t154 - pkin(4);
t134 = Ifges(4,1) * t152 + Ifges(4,4) * t155;
t132 = Ifges(4,4) * t152 + Ifges(4,2) * t155;
t131 = Ifges(6,2) * t154 + t178;
t128 = t180 * t154;
t125 = t180 * t151;
t99 = -mrSges(4,1) * t171 - mrSges(4,3) * t111;
t98 = mrSges(4,2) * t171 + mrSges(4,3) * t110;
t96 = -t125 * t157 + t128 * t160;
t93 = -t125 * t160 - t128 * t157;
t92 = Ifges(5,1) * t124 - Ifges(5,4) * t122;
t90 = Ifges(5,4) * t124 - Ifges(5,2) * t122;
t84 = mrSges(6,1) * t122 - mrSges(6,3) * t173;
t83 = -mrSges(6,2) * t122 - mrSges(6,3) * t174;
t65 = Ifges(4,1) * t111 + Ifges(4,4) * t110 - Ifges(4,5) * t171;
t64 = Ifges(4,4) * t111 + Ifges(4,2) * t110 - Ifges(4,6) * t171;
t63 = pkin(5) * t174 + t94;
t60 = -mrSges(5,1) * t171 - mrSges(5,3) * t74;
t59 = mrSges(5,2) * t171 - mrSges(5,3) * t73;
t56 = Ifges(6,3) * t122 + (Ifges(6,5) * t154 - Ifges(6,6) * t151) * t124;
t52 = mrSges(7,1) * t122 - mrSges(7,3) * t76;
t51 = -mrSges(7,2) * t122 - mrSges(7,3) * t75;
t41 = -pkin(10) * t174 + t47;
t40 = Ifges(5,1) * t74 - Ifges(5,4) * t73 - Ifges(5,5) * t171;
t39 = Ifges(5,4) * t74 - Ifges(5,2) * t73 - Ifges(5,6) * t171;
t38 = mrSges(6,1) * t73 - mrSges(6,3) * t55;
t37 = -mrSges(6,2) * t73 + mrSges(6,3) * t54;
t36 = Ifges(7,1) * t76 - Ifges(7,4) * t75 + Ifges(7,5) * t122;
t35 = Ifges(7,4) * t76 - Ifges(7,2) * t75 + Ifges(7,6) * t122;
t33 = pkin(5) * t122 - pkin(10) * t173 + t46;
t18 = Ifges(6,4) * t55 + Ifges(6,2) * t54 + Ifges(6,6) * t73;
t17 = Ifges(6,5) * t55 + Ifges(6,6) * t54 + Ifges(6,3) * t73;
t15 = mrSges(7,1) * t73 - mrSges(7,3) * t31;
t14 = -mrSges(7,2) * t73 + mrSges(7,3) * t30;
t13 = t157 * t33 + t160 * t41;
t12 = -t157 * t41 + t160 * t33;
t11 = -t54 * pkin(5) + t21;
t9 = Ifges(7,1) * t31 + Ifges(7,4) * t30 + Ifges(7,5) * t73;
t8 = Ifges(7,4) * t31 + Ifges(7,2) * t30 + Ifges(7,6) * t73;
t4 = pkin(10) * t54 + t6;
t3 = pkin(5) * t73 - pkin(10) * t55 + t5;
t2 = t157 * t3 + t160 * t4;
t1 = -t157 * t4 + t160 * t3;
t16 = [0.2e1 * t2 * t14 + 0.2e1 * t1 * t15 + 0.2e1 * t11 * t10 + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(6) * (t21 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2 + t79 ^ 2) + m(4) * (t106 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(3) * (pkin(1) ^ 2 * t153 ^ 2 + t112 ^ 2 + t113 ^ 2) + (t17 + t7 - t39) * t73 + t30 * t8 + t31 * t9 + 0.2e1 * t21 * t32 + 0.2e1 * t6 * t37 + 0.2e1 * t5 * t38 + (t166 - 0.2e1 * t175 + 0.2e1 * t176) * t156 + t54 * t18 + Ifges(2,3) + t55 * t19 + 0.2e1 * t23 * t59 + 0.2e1 * t22 * t60 + ((-0.2e1 * t112 * mrSges(3,3) + Ifges(3,5) * t156 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t159) * t153) * t159 + (0.2e1 * t113 * mrSges(3,3) - Ifges(4,5) * t111 + Ifges(3,6) * t156 - Ifges(4,6) * t110 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t159 + (Ifges(5,3) + Ifges(4,3) + Ifges(3,2)) * t161) * t153 + t179) * t161) * t153 + t74 * t40 + 0.2e1 * t79 * t42 + 0.2e1 * t62 * t98 + 0.2e1 * t61 * t99 + 0.2e1 * t106 * t80 + t110 * t64 + t111 * t65; t13 * t14 + t12 * t15 + (-t23 * mrSges(5,3) + t17 / 0.2e1 + t7 / 0.2e1 - t39 / 0.2e1) * t122 - t175 + t176 + t55 * t190 + t54 * t191 + t164 * mrSges(4,3) + m(4) * (-pkin(2) * t106 + qJ(3) * t164) - (Ifges(4,5) * t152 + Ifges(4,6) * t155 + t170) * t171 / 0.2e1 + t166 + t30 * t35 / 0.2e1 + t31 * t36 / 0.2e1 + t11 * t43 + t46 * t38 + t47 * t37 + t2 * t51 + t1 * t52 + m(7) * (t1 * t12 + t11 * t63 + t13 * t2) + m(6) * (t21 * t94 + t46 * t5 + t47 * t6) + m(5) * (t143 * t79 - t22 * t94 + t23 * t97) + t63 * t10 - t75 * t8 / 0.2e1 + t76 * t9 / 0.2e1 - pkin(2) * t80 + t21 * t81 + t6 * t83 + t5 * t84 + t79 * t87 + t74 * t92 / 0.2e1 + t97 * t59 + (-t22 * mrSges(5,3) - t151 * t18 / 0.2e1 + t154 * t192 + t40 / 0.2e1) * t124 + t106 * t127 + t110 * t132 / 0.2e1 + t111 * t134 / 0.2e1 + t143 * t42 + t152 * t65 / 0.2e1 + t155 * t64 / 0.2e1 + (-t152 * t99 + t155 * t98) * qJ(3) + (t34 / 0.2e1 + t56 / 0.2e1 - t90 / 0.2e1) * t73 + (t32 - t60) * t94; -0.2e1 * pkin(2) * t127 + 0.2e1 * t12 * t52 + 0.2e1 * t13 * t51 + t155 * t132 + t152 * t134 + 0.2e1 * t143 * t87 - t75 * t35 + t76 * t36 + 0.2e1 * t63 * t43 + 0.2e1 * t46 * t84 + 0.2e1 * t47 * t83 + t81 * t193 + Ifges(3,3) + 0.2e1 * t168 * qJ(3) * mrSges(4,3) + (mrSges(5,3) * t193 - t151 * t57 + t154 * t58 + t92) * t124 + (-0.2e1 * mrSges(5,3) * t97 + t34 + t56 - t90) * t122 + m(7) * (t12 ^ 2 + t13 ^ 2 + t63 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2 + t194) + m(5) * (t143 ^ 2 + t97 ^ 2 + t194) + m(4) * (qJ(3) ^ 2 * t168 + pkin(2) ^ 2); t121 * t15 + t123 * t14 + t151 * t37 + t154 * t38 + m(7) * (t1 * t121 + t123 * t2) + m(6) * (t151 * t6 + t154 * t5) + m(5) * t79 + m(4) * t106 + t80 + t42; -m(4) * pkin(2) + t121 * t52 + t123 * t51 + t151 * t83 + t154 * t84 + m(7) * (t12 * t121 + t123 * t13) + m(6) * (t151 * t47 + t154 * t46) + m(5) * t143 + t87 + t127; m(4) + m(5) + m(6) * t169 + m(7) * (t121 ^ 2 + t123 ^ 2); m(7) * (t1 * t93 + t11 * t142 + t2 * t96) + (-t1 * t123 + t121 * t2) * mrSges(7,3) - t179 + t167 * t73 + m(6) * (-pkin(4) * t21 + (-t151 * t5 + t154 * t6) * qJ(5)) + (qJ(5) * t37 + t6 * mrSges(6,3) + t18 / 0.2e1) * t154 + (-t5 * mrSges(6,3) - qJ(5) * t38 + t192) * t151 + t22 * mrSges(5,1) - t23 * mrSges(5,2) - Ifges(5,3) * t171 - pkin(4) * t32 + t11 * t86 + t30 * t189 + t31 * t188 + t93 * t15 + t96 * t14 + t8 * t187 + t9 * t186 + t21 * t126 + t54 * t131 / 0.2e1 + t55 * t185 + t142 * t10; -pkin(4) * t81 + t63 * t86 - t75 * t189 + t76 * t188 + t93 * t52 + t96 * t51 - t97 * mrSges(5,2) + t35 * t187 + t36 * t186 + t142 * t43 + (-mrSges(5,1) + t126) * t94 + t167 * t122 + (-t12 * t123 + t121 * t13) * mrSges(7,3) + (t47 * mrSges(6,3) + qJ(5) * t83 + t124 * t185 + t191) * t154 + (-t46 * mrSges(6,3) - qJ(5) * t84 - t124 * t131 / 0.2e1 + t190) * t151 + m(7) * (t12 * t93 + t13 * t96 + t142 * t63) + m(6) * (-pkin(4) * t94 + (-t151 * t46 + t154 * t47) * qJ(5)) + t170; m(7) * (t121 * t93 + t123 * t96); -0.2e1 * pkin(4) * t126 + t121 * t89 + t123 * t91 + t154 * t131 + t151 * t133 + 0.2e1 * t142 * t86 + Ifges(5,3) + m(7) * (t142 ^ 2 + t93 ^ 2 + t96 ^ 2) + m(6) * (qJ(5) ^ 2 * t169 + pkin(4) ^ 2) + 0.2e1 * (t121 * t96 - t123 * t93) * mrSges(7,3) + 0.2e1 * t169 * qJ(5) * mrSges(6,3); m(6) * t21 + m(7) * t11 + t10 + t32; m(6) * t94 + m(7) * t63 + t43 + t81; 0; -m(6) * pkin(4) + m(7) * t142 + t126 + t86; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t12 - mrSges(7,2) * t13 + t34; -t86; mrSges(7,1) * t93 - t96 * mrSges(7,2) + t88; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
