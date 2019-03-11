% Calculate joint inertia matrix for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:26
% EndTime: 2019-03-09 21:57:29
% DurationCPUTime: 1.19s
% Computational Cost: add. (3602->302), mult. (6744->432), div. (0->0), fcn. (7732->10), ass. (0->134)
t130 = sin(pkin(11));
t131 = cos(pkin(11));
t132 = sin(qJ(6));
t136 = cos(qJ(6));
t98 = -t132 * t130 + t136 * t131;
t180 = mrSges(7,3) * t98;
t138 = cos(qJ(3));
t121 = t138 * pkin(2) + pkin(3);
t133 = sin(qJ(4));
t137 = cos(qJ(4));
t134 = sin(qJ(3));
t179 = pkin(2) * t134;
t86 = t133 * t121 + t137 * t179;
t83 = qJ(5) + t86;
t73 = (-pkin(10) - t83) * t130;
t125 = t131 * pkin(10);
t163 = t131 * t83;
t74 = t125 + t163;
t44 = t132 * t73 + t136 * t74;
t38 = t44 * t180;
t99 = t136 * t130 + t132 * t131;
t66 = -t98 * mrSges(7,1) + t99 * mrSges(7,2);
t177 = t131 * pkin(5);
t85 = t137 * t121 - t133 * t179;
t84 = -pkin(4) - t85;
t79 = t84 - t177;
t45 = t79 * t66;
t108 = -t131 * mrSges(6,1) + t130 * mrSges(6,2);
t69 = t84 * t108;
t126 = t130 ^ 2;
t171 = mrSges(6,3) * t126;
t77 = t83 * t171;
t127 = t131 ^ 2;
t170 = mrSges(6,3) * t127;
t78 = t83 * t170;
t82 = t85 * mrSges(5,1);
t187 = t38 + t45 + t69 + t77 + t78 + t82;
t176 = t133 * pkin(3);
t116 = qJ(5) + t176;
t104 = t116 * t171;
t105 = t116 * t170;
t175 = t137 * pkin(3);
t123 = mrSges(5,1) * t175;
t89 = (-pkin(10) - t116) * t130;
t90 = t131 * t116 + t125;
t58 = t132 * t89 + t136 * t90;
t46 = t58 * t180;
t117 = -pkin(4) - t177;
t106 = t117 - t175;
t54 = t106 * t66;
t120 = -pkin(4) - t175;
t87 = t120 * t108;
t186 = t104 + t105 + t123 + t46 + t54 + t87;
t135 = sin(qJ(2));
t139 = cos(qJ(2));
t101 = t134 * t139 + t138 * t135;
t181 = -pkin(8) - pkin(7);
t155 = t181 * t135;
t156 = t181 * t139;
t75 = t134 * t156 + t138 * t155;
t146 = -t101 * pkin(9) + t75;
t100 = -t134 * t135 + t138 * t139;
t76 = t134 * t155 - t138 * t156;
t55 = t100 * pkin(9) + t76;
t30 = t133 * t55 - t137 * t146;
t185 = t30 ^ 2;
t184 = 0.2e1 * t30;
t122 = -t139 * pkin(2) - pkin(1);
t80 = -t100 * pkin(3) + t122;
t183 = 0.2e1 * t80;
t178 = pkin(4) * t108;
t174 = t86 * mrSges(5,2);
t173 = t99 * mrSges(7,3);
t64 = -t137 * t100 + t133 * t101;
t65 = t133 * t100 + t137 * t101;
t28 = t64 * pkin(4) - t65 * qJ(5) + t80;
t32 = t133 * t146 + t137 * t55;
t12 = t130 * t28 + t131 * t32;
t164 = t131 * t65;
t165 = t130 * t65;
t39 = mrSges(6,1) * t165 + mrSges(6,2) * t164;
t172 = Ifges(7,5) * t99 + Ifges(7,6) * t98;
t169 = Ifges(6,4) * t130;
t168 = Ifges(6,4) * t131;
t11 = -t130 * t32 + t131 * t28;
t167 = t11 * t130;
t166 = t12 * t131;
t162 = t126 + t127;
t161 = t135 ^ 2 + t139 ^ 2;
t160 = -0.2e1 * t173;
t159 = mrSges(5,2) * t176;
t36 = t99 * t65;
t37 = t98 * t65;
t158 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t64;
t157 = mrSges(6,3) * t167;
t13 = t36 * mrSges(7,1) + t37 * mrSges(7,2);
t154 = t162 * t116;
t110 = Ifges(6,2) * t131 + t169;
t111 = Ifges(6,1) * t130 + t168;
t67 = Ifges(7,4) * t99 + Ifges(7,2) * t98;
t68 = Ifges(7,1) * t99 + Ifges(7,4) * t98;
t153 = t131 * t110 + t130 * t111 + t98 * t67 + t99 * t68 + Ifges(5,3);
t152 = t166 - t167;
t40 = -t64 * mrSges(6,2) - mrSges(6,3) * t165;
t41 = t64 * mrSges(6,1) - mrSges(6,3) * t164;
t151 = -t130 * t41 + t131 * t40;
t150 = Ifges(4,3) + t153;
t149 = (t138 * mrSges(4,1) - t134 * mrSges(4,2)) * pkin(2);
t148 = t108 + t66;
t118 = qJ(5) * t171;
t119 = qJ(5) * t170;
t107 = (-pkin(10) - qJ(5)) * t130;
t109 = t131 * qJ(5) + t125;
t71 = t132 * t107 + t136 * t109;
t53 = t71 * t180;
t56 = t117 * t66;
t147 = t118 + t119 + t153 + t53 + t56 - t178;
t15 = pkin(5) * t165 + t30;
t4 = t64 * pkin(5) - pkin(10) * t164 + t11;
t5 = -pkin(10) * t165 + t12;
t2 = -t132 * t5 + t136 * t4;
t22 = Ifges(6,6) * t64 + (-Ifges(6,2) * t130 + t168) * t65;
t23 = Ifges(6,5) * t64 + (Ifges(6,1) * t131 - t169) * t65;
t3 = t132 * t4 + t136 * t5;
t8 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t64;
t9 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t64;
t145 = -t32 * mrSges(5,2) - t2 * t173 + t3 * t180 + mrSges(6,3) * t166 + t15 * t66 + t130 * t23 / 0.2e1 + t131 * t22 / 0.2e1 - t36 * t67 / 0.2e1 + t37 * t68 / 0.2e1 - t110 * t165 / 0.2e1 + t111 * t164 / 0.2e1 + t98 * t8 / 0.2e1 - Ifges(5,6) * t64 + Ifges(5,5) * t65 + t99 * t9 / 0.2e1 + (t108 - mrSges(5,1)) * t30 + (Ifges(6,5) * t130 + Ifges(6,6) * t131 + t172) * t64 / 0.2e1;
t144 = t75 * mrSges(4,1) - t76 * mrSges(4,2) + Ifges(4,5) * t101 + Ifges(4,6) * t100 + t145;
t70 = t136 * t107 - t132 * t109;
t57 = -t132 * t90 + t136 * t89;
t43 = -t132 * t74 + t136 * t73;
t17 = mrSges(7,1) * t64 - mrSges(7,3) * t37;
t16 = -mrSges(7,2) * t64 - mrSges(7,3) * t36;
t1 = [-t36 * t8 + t37 * t9 + t39 * t184 + 0.2e1 * t12 * t40 + 0.2e1 * t11 * t41 + 0.2e1 * t2 * t17 + 0.2e1 * t15 * t13 + 0.2e1 * t3 * t16 + Ifges(2,3) - 0.2e1 * pkin(1) * (-t139 * mrSges(3,1) + t135 * mrSges(3,2)) + t135 * (Ifges(3,1) * t135 + Ifges(3,4) * t139) + t139 * (Ifges(3,4) * t135 + Ifges(3,2) * t139) + 0.2e1 * t122 * (-t100 * mrSges(4,1) + t101 * mrSges(4,2)) + t100 * (Ifges(4,4) * t101 + Ifges(4,2) * t100) + t101 * (Ifges(4,1) * t101 + Ifges(4,4) * t100) + (mrSges(5,2) * t183 + mrSges(5,3) * t184 + Ifges(5,1) * t65 - t130 * t22 + t131 * t23) * t65 + (mrSges(5,1) * t183 - 0.2e1 * t32 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t64 + (Ifges(6,5) * t131 - Ifges(6,6) * t130 - (2 * Ifges(5,4))) * t65 + t158) * t64 + m(3) * (pkin(7) ^ 2 * t161 + pkin(1) ^ 2) + m(4) * (t122 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t32 ^ 2 + t80 ^ 2 + t185) + m(7) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2 + t185) + 0.2e1 * (t76 * t100 - t75 * t101) * mrSges(4,3) + 0.2e1 * t161 * pkin(7) * mrSges(3,3); t40 * t163 + m(6) * (t152 * t83 + t84 * t30) + t43 * t17 + t44 * t16 + (-t11 * mrSges(6,3) - t83 * t41) * t130 + m(5) * (-t85 * t30 + t86 * t32) + m(7) * (t79 * t15 + t43 * t2 + t44 * t3) + (-t135 * mrSges(3,1) - t139 * mrSges(3,2)) * pkin(7) + (-t86 * t64 - t85 * t65) * mrSges(5,3) + t144 + Ifges(3,6) * t139 + Ifges(3,5) * t135 + t79 * t13 + t84 * t39 + (m(4) * (t134 * t76 + t138 * t75) + (t100 * t134 - t101 * t138) * mrSges(4,3)) * pkin(2); 0.2e1 * t45 + 0.2e1 * t38 + t150 + t43 * t160 + m(4) * (t134 ^ 2 + t138 ^ 2) * pkin(2) ^ 2 + m(6) * (t162 * t83 ^ 2 + t84 ^ 2) + 0.2e1 * t149 + 0.2e1 * t82 + Ifges(3,3) + 0.2e1 * t69 + m(7) * (t43 ^ 2 + t44 ^ 2 + t79 ^ 2) + m(5) * (t85 ^ 2 + t86 ^ 2) + 0.2e1 * t77 + 0.2e1 * t78 - 0.2e1 * t174; -t157 + t151 * t116 + (m(5) * (t133 * t32 - t137 * t30) + (-t133 * t64 - t137 * t65) * mrSges(5,3)) * pkin(3) + m(7) * (t106 * t15 + t57 * t2 + t58 * t3) + m(6) * (t116 * t152 + t120 * t30) + t57 * t17 + t58 * t16 + t144 + t120 * t39 + t106 * t13; t150 + (-t86 - t176) * mrSges(5,2) + m(6) * (t120 * t84 + t154 * t83) + t149 + (-t43 - t57) * t173 + m(7) * (t106 * t79 + t57 * t43 + t58 * t44) + t186 + m(5) * (t133 * t86 + t137 * t85) * pkin(3) + t187; -0.2e1 * t159 + t57 * t160 + 0.2e1 * t104 + 0.2e1 * t105 + 0.2e1 * t123 + 0.2e1 * t46 + 0.2e1 * t54 + 0.2e1 * t87 + m(7) * (t106 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t116 ^ 2 * t162 + t120 ^ 2) + m(5) * (t133 ^ 2 + t137 ^ 2) * pkin(3) ^ 2 + t150; -t157 + m(6) * (-pkin(4) * t30 + qJ(5) * t152) + t145 - pkin(4) * t39 + t151 * qJ(5) + m(7) * (t117 * t15 + t70 * t2 + t71 * t3) + t117 * t13 + t70 * t17 + t71 * t16; m(6) * (qJ(5) * t162 * t83 - pkin(4) * t84) + t147 + (-t43 - t70) * t173 + m(7) * (t117 * t79 + t70 * t43 + t71 * t44) - t174 + t187; -t159 + m(6) * (-pkin(4) * t120 + qJ(5) * t154) + t147 + (-t57 - t70) * t173 + m(7) * (t117 * t106 + t70 * t57 + t71 * t58) + t186; t70 * t160 - 0.2e1 * t178 + 0.2e1 * t118 + 0.2e1 * t119 + 0.2e1 * t53 + 0.2e1 * t56 + m(7) * (t117 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (qJ(5) ^ 2 * t162 + pkin(4) ^ 2) + t153; m(6) * t30 + m(7) * t15 + t13 + t39; m(6) * t84 + m(7) * t79 + t148; m(6) * t120 + m(7) * t106 + t148; -m(6) * pkin(4) + m(7) * t117 + t148; m(6) + m(7); mrSges(7,1) * t2 - mrSges(7,2) * t3 + t158; t43 * mrSges(7,1) - t44 * mrSges(7,2) + t172; t57 * mrSges(7,1) - t58 * mrSges(7,2) + t172; t70 * mrSges(7,1) - t71 * mrSges(7,2) + t172; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
