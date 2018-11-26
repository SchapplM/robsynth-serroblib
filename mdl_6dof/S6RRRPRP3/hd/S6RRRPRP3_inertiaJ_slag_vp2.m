% Calculate joint inertia matrix for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:42:13
% EndTime: 2018-11-23 17:42:14
% DurationCPUTime: 1.39s
% Computational Cost: add. (2358->280), mult. (4506->373), div. (0->0), fcn. (4855->8), ass. (0->108)
t195 = Ifges(6,1) + Ifges(7,1);
t194 = Ifges(6,4) - Ifges(7,5);
t193 = Ifges(7,4) + Ifges(6,5);
t192 = Ifges(6,6) - Ifges(7,6);
t186 = (mrSges(6,3) + mrSges(7,2));
t191 = 2 * t186;
t190 = Ifges(7,2) + Ifges(6,3);
t133 = sin(qJ(3));
t134 = sin(qJ(2));
t135 = cos(qJ(3));
t136 = cos(qJ(2));
t108 = t133 * t134 - t135 * t136;
t130 = sin(pkin(10));
t131 = cos(pkin(10));
t132 = sin(qJ(5));
t173 = cos(qJ(5));
t107 = t173 * t130 + t132 * t131;
t109 = t133 * t136 + t134 * t135;
t55 = t107 * t109;
t147 = -t132 * t130 + t173 * t131;
t56 = t147 * t109;
t189 = t193 * t108 - t194 * t55 + t195 * t56;
t188 = t195 * t107 + t194 * t147;
t123 = -pkin(2) * t136 - pkin(1);
t66 = pkin(3) * t108 - qJ(4) * t109 + t123;
t176 = -pkin(8) - pkin(7);
t116 = t176 * t136;
t157 = t176 * t134;
t85 = -t135 * t116 + t133 * t157;
t29 = -t130 * t85 + t131 * t66;
t30 = t130 * t66 + t131 * t85;
t152 = -t130 * t29 + t131 * t30;
t187 = t193 * t107 + t192 * t147;
t185 = -m(7) * pkin(5) - mrSges(7,1);
t184 = -mrSges(6,1) + t185;
t183 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t83 = -t116 * t133 - t135 * t157;
t182 = t83 ^ 2;
t69 = -mrSges(7,1) * t147 - t107 * mrSges(7,3);
t181 = 0.2e1 * t69;
t70 = -mrSges(6,1) * t147 + t107 * mrSges(6,2);
t180 = 0.2e1 * t70;
t179 = 0.2e1 * t83;
t178 = 0.2e1 * t123;
t172 = pkin(2) * t135;
t162 = t109 * t131;
t14 = pkin(4) * t108 - pkin(9) * t162 + t29;
t163 = t109 * t130;
t20 = -pkin(9) * t163 + t30;
t6 = t132 * t14 + t173 * t20;
t33 = mrSges(6,1) * t108 - mrSges(6,3) * t56;
t34 = -t108 * mrSges(7,1) + t56 * mrSges(7,2);
t171 = t34 - t33;
t32 = -mrSges(6,2) * t108 - mrSges(6,3) * t55;
t35 = -mrSges(7,2) * t55 + mrSges(7,3) * t108;
t170 = t35 + t32;
t65 = mrSges(5,1) * t163 + mrSges(5,2) * t162;
t167 = Ifges(5,4) * t130;
t166 = Ifges(5,4) * t131;
t93 = t107 * mrSges(7,2);
t161 = t130 ^ 2 + t131 ^ 2;
t160 = t134 ^ 2 + t136 ^ 2;
t119 = pkin(2) * t133 + qJ(4);
t154 = (-pkin(9) - t119) * t130;
t125 = t131 * pkin(9);
t91 = t119 * t131 + t125;
t60 = t132 * t91 - t173 * t154;
t62 = t132 * t154 + t173 * t91;
t159 = t60 ^ 2 + t62 ^ 2;
t113 = qJ(4) * t131 + t125;
t155 = (-pkin(9) - qJ(4)) * t130;
t79 = t113 * t132 - t173 * t155;
t81 = t173 * t113 + t132 * t155;
t158 = t79 ^ 2 + t81 ^ 2;
t120 = -pkin(4) * t131 - pkin(3);
t22 = t55 * mrSges(6,1) + t56 * mrSges(6,2);
t21 = t55 * mrSges(7,1) - t56 * mrSges(7,3);
t156 = t60 * t79 + t81 * t62;
t112 = -t131 * mrSges(5,1) + t130 * mrSges(5,2);
t153 = t161 * t119;
t46 = pkin(4) * t163 + t83;
t67 = -mrSges(5,2) * t108 - mrSges(5,3) * t163;
t68 = mrSges(5,1) * t108 - mrSges(5,3) * t162;
t151 = -t130 * t68 + t131 * t67;
t150 = 0.2e1 * t161 * mrSges(5,3);
t149 = t190 * t108 - t192 * t55 + t193 * t56;
t5 = -t132 * t20 + t173 * t14;
t148 = (t135 * mrSges(4,1) - t133 * mrSges(4,2)) * pkin(2);
t114 = Ifges(5,2) * t131 + t167;
t115 = Ifges(5,1) * t130 + t166;
t71 = Ifges(7,5) * t107 - Ifges(7,3) * t147;
t72 = Ifges(6,4) * t107 + Ifges(6,2) * t147;
t146 = t131 * t114 + t130 * t115 + Ifges(4,3) + t188 * t107 + (-t71 + t72) * t147;
t64 = -pkin(5) * t147 - qJ(6) * t107 + t120;
t142 = t112 + t69 + t70;
t141 = (-pkin(5) * t107 + qJ(6) * t147) * mrSges(7,2) + t187;
t15 = Ifges(7,5) * t56 + Ifges(7,6) * t108 + Ifges(7,3) * t55;
t16 = Ifges(6,4) * t56 - Ifges(6,2) * t55 + Ifges(6,6) * t108;
t3 = qJ(6) * t108 + t6;
t38 = Ifges(5,6) * t108 + (-Ifges(5,2) * t130 + t166) * t109;
t39 = Ifges(5,5) * t108 + (Ifges(5,1) * t131 - t167) * t109;
t4 = -t108 * pkin(5) - t5;
t8 = pkin(5) * t55 - qJ(6) * t56 + t46;
t140 = t131 * t38 / 0.2e1 + t130 * t39 / 0.2e1 - Ifges(4,6) * t108 + Ifges(4,5) * t109 - t85 * mrSges(4,2) + t8 * t69 + t46 * t70 + t115 * t162 / 0.2e1 - t114 * t163 / 0.2e1 + t4 * t93 + (t112 - mrSges(4,1)) * t83 + (t71 / 0.2e1 - t72 / 0.2e1) * t55 + t188 * t56 / 0.2e1 + t152 * mrSges(5,3) + (Ifges(5,5) * t130 + Ifges(5,6) * t131 + t187) * t108 / 0.2e1 + (-t15 / 0.2e1 + t16 / 0.2e1 + t6 * mrSges(6,3) + t3 * mrSges(7,2)) * t147 + (-t5 * mrSges(6,3) + t189 / 0.2e1) * t107;
t122 = -pkin(3) - t172;
t111 = t120 - t172;
t47 = t64 - t172;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t136 + mrSges(3,2) * t134) + t134 * (Ifges(3,1) * t134 + Ifges(3,4) * t136) + t136 * (Ifges(3,4) * t134 + Ifges(3,2) * t136) + t65 * t179 + 0.2e1 * t30 * t67 + 0.2e1 * t29 * t68 + 0.2e1 * t3 * t35 + 0.2e1 * t46 * t22 + 0.2e1 * t8 * t21 + 0.2e1 * t6 * t32 + 0.2e1 * t5 * t33 + 0.2e1 * t4 * t34 + Ifges(2,3) + t189 * t56 + (t15 - t16) * t55 + 0.2e1 * t160 * pkin(7) * mrSges(3,3) + (mrSges(4,2) * t178 + mrSges(4,3) * t179 + Ifges(4,1) * t109 - t130 * t38 + t131 * t39) * t109 + (mrSges(4,1) * t178 - 0.2e1 * t85 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t108 + (Ifges(5,5) * t131 - Ifges(5,6) * t130 - (2 * Ifges(4,4))) * t109 + t149) * t108 + m(3) * (t160 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t123 ^ 2 + t85 ^ 2 + t182) + m(5) * (t29 ^ 2 + t30 ^ 2 + t182) + m(6) * (t46 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2); t140 + m(5) * (t152 * t119 + t122 * t83) + t170 * t62 + t171 * t60 + t151 * t119 + Ifges(3,5) * t134 + Ifges(3,6) * t136 + t122 * t65 + t111 * t22 + t47 * t21 + m(6) * (t111 * t46 - t5 * t60 + t6 * t62) + m(7) * (t3 * t62 + t4 * t60 + t47 * t8) + (-t134 * mrSges(3,1) - t136 * mrSges(3,2)) * pkin(7) + (m(4) * (t133 * t85 - t135 * t83) + (-t108 * t133 - t109 * t135) * mrSges(4,3)) * pkin(2); t146 + m(5) * (t161 * t119 ^ 2 + t122 ^ 2) + m(6) * (t111 ^ 2 + t159) + m(7) * (t47 ^ 2 + t159) + 0.2e1 * t122 * t112 + t111 * t180 + t47 * t181 + t119 * t150 + 0.2e1 * t148 + m(4) * (t133 ^ 2 + t135 ^ 2) * pkin(2) ^ 2 + Ifges(3,3) + (t60 * t107 + t147 * t62) * t191; t140 + m(5) * (-pkin(3) * t83 + t152 * qJ(4)) + t120 * t22 + t64 * t21 - pkin(3) * t65 + t170 * t81 + t171 * t79 + t151 * qJ(4) + m(6) * (t120 * t46 - t5 * t79 + t6 * t81) + m(7) * (t3 * t81 + t4 * t79 + t64 * t8); t146 + (t111 + t120) * t70 + (t64 + t47) * t69 + (t122 - pkin(3)) * t112 + t148 + m(5) * (-pkin(3) * t122 + qJ(4) * t153) + (t161 * qJ(4) + t153) * mrSges(5,3) + m(6) * (t111 * t120 + t156) + m(7) * (t47 * t64 + t156) + t186 * ((t60 + t79) * t107 - (-t62 - t81) * t147); -0.2e1 * pkin(3) * t112 + t120 * t180 + t64 * t181 + qJ(4) * t150 + m(7) * (t64 ^ 2 + t158) + m(6) * (t120 ^ 2 + t158) + m(5) * (t161 * qJ(4) ^ 2 + pkin(3) ^ 2) + t146 + (t79 * t107 + t147 * t81) * t191; m(5) * t83 + m(6) * t46 + m(7) * t8 + t21 + t22 + t65; m(5) * t122 + m(6) * t111 + m(7) * t47 + t142; -m(5) * pkin(3) + m(6) * t120 + m(7) * t64 + t142; m(5) + m(6) + m(7); -pkin(5) * t34 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) + qJ(6) * t35 + t3 * mrSges(7,3) - t6 * mrSges(6,2) + t5 * mrSges(6,1) - t4 * mrSges(7,1) + t149; t183 * t62 + t184 * t60 + t141; t183 * t81 + t184 * t79 + t141; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t190; m(7) * t4 + t34; m(7) * t60 + t93; m(7) * t79 + t93; 0; t185; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
