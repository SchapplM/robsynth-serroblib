% Calculate joint inertia matrix for
% S6RRRPRP1
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:55
% EndTime: 2019-03-09 16:30:59
% DurationCPUTime: 1.58s
% Computational Cost: add. (2233->259), mult. (4166->355), div. (0->0), fcn. (4634->8), ass. (0->98)
t173 = Ifges(6,4) + Ifges(7,4);
t168 = Ifges(6,6) + Ifges(7,6);
t102 = sin(pkin(10));
t103 = cos(pkin(10));
t105 = sin(qJ(3));
t106 = sin(qJ(2));
t108 = cos(qJ(3));
t109 = cos(qJ(2));
t72 = -t105 * t106 + t108 * t109;
t73 = t105 * t109 + t106 * t108;
t49 = t102 * t73 - t103 * t72;
t172 = t168 * t49;
t171 = Ifges(6,1) + Ifges(7,1);
t170 = Ifges(6,5) + Ifges(7,5);
t169 = Ifges(6,2) + Ifges(7,2);
t107 = cos(qJ(5));
t167 = t173 * t107;
t104 = sin(qJ(5));
t166 = t173 * t104;
t165 = Ifges(6,3) + Ifges(7,3);
t50 = t102 * t72 + t103 * t73;
t164 = (-t169 * t104 + t167) * t50 + t172;
t163 = (t171 * t107 - t166) * t50 + t170 * t49;
t162 = t169 * t107 + t166;
t161 = t171 * t104 + t167;
t122 = t170 * t104 + t168 * t107;
t129 = t104 ^ 2 + t107 ^ 2;
t160 = 0.2e1 * t129;
t152 = -pkin(8) - pkin(7);
t126 = t152 * t106;
t127 = t152 * t109;
t54 = t105 * t127 + t108 * t126;
t114 = -t73 * qJ(4) + t54;
t55 = t105 * t126 - t108 * t127;
t39 = qJ(4) * t72 + t55;
t23 = t102 * t39 - t103 * t114;
t159 = t23 ^ 2;
t158 = 0.2e1 * t23;
t91 = -pkin(2) * t109 - pkin(1);
t57 = -pkin(3) * t72 + t91;
t157 = 0.2e1 * t57;
t93 = t104 * mrSges(7,2);
t77 = -t107 * mrSges(7,1) + t93;
t156 = 0.2e1 * t77;
t140 = mrSges(6,2) * t104;
t78 = -mrSges(6,1) * t107 + t140;
t155 = 0.2e1 * t78;
t154 = m(7) * pkin(5);
t149 = pkin(2) * t105;
t148 = pkin(5) * t107;
t22 = pkin(4) * t49 - pkin(9) * t50 + t57;
t25 = t102 * t114 + t103 * t39;
t6 = t104 * t22 + t107 * t25;
t147 = t107 * t6;
t90 = pkin(2) * t108 + pkin(3);
t61 = -t102 * t149 + t103 * t90;
t146 = t61 * mrSges(5,1);
t62 = t102 * t90 + t103 * t149;
t145 = t62 * mrSges(5,2);
t133 = t107 * t50;
t134 = t104 * t50;
t26 = mrSges(7,1) * t134 + mrSges(7,2) * t133;
t135 = t104 * mrSges(7,3);
t60 = pkin(9) + t62;
t132 = t107 * t60;
t88 = pkin(3) * t102 + pkin(9);
t131 = t107 * t88;
t130 = t106 ^ 2 + t109 ^ 2;
t92 = t107 * qJ(6);
t128 = 0.2e1 * mrSges(7,3);
t89 = -pkin(3) * t103 - pkin(4);
t5 = -t104 * t25 + t107 * t22;
t123 = t170 * t133 + t165 * t49;
t59 = -pkin(4) - t61;
t1 = pkin(5) * t49 - t50 * t92 + t5;
t121 = -t5 * mrSges(6,3) - t1 * mrSges(7,3);
t120 = -t104 * t5 + t147;
t119 = mrSges(6,3) * t160;
t118 = t103 * mrSges(5,1) - t102 * mrSges(5,2);
t117 = mrSges(6,1) * t104 + mrSges(6,2) * t107;
t116 = t104 * t161 + t107 * t162 + Ifges(4,3) + Ifges(5,3);
t115 = (mrSges(4,1) * t108 - mrSges(4,2) * t105) * pkin(2);
t3 = -qJ(6) * t134 + t6;
t8 = pkin(5) * t134 + t23;
t113 = t54 * mrSges(4,1) - t55 * mrSges(4,2) - t25 * mrSges(5,2) + mrSges(6,3) * t147 + Ifges(4,5) * t73 + Ifges(5,5) * t50 + Ifges(4,6) * t72 + t8 * t77 + (-mrSges(5,1) + t78) * t23 + t163 * t104 / 0.2e1 - t162 * t134 / 0.2e1 + t161 * t133 / 0.2e1 + (-Ifges(5,6) + t122 / 0.2e1) * t49 + (t3 * mrSges(7,3) + t164 / 0.2e1) * t107;
t76 = t89 - t148;
t65 = t92 + t131;
t64 = (-qJ(6) - t88) * t104;
t56 = t59 - t148;
t53 = t92 + t132;
t52 = (-qJ(6) - t60) * t104;
t44 = t50 * mrSges(5,2);
t31 = mrSges(6,1) * t49 - mrSges(6,3) * t133;
t30 = mrSges(7,1) * t49 - mrSges(7,3) * t133;
t29 = -mrSges(6,2) * t49 - mrSges(6,3) * t134;
t28 = -mrSges(7,2) * t49 - mrSges(7,3) * t134;
t27 = t117 * t50;
t2 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t109 + mrSges(3,2) * t106) + t106 * (Ifges(3,1) * t106 + Ifges(3,4) * t109) + t109 * (Ifges(3,4) * t106 + Ifges(3,2) * t109) + 0.2e1 * t91 * (-mrSges(4,1) * t72 + mrSges(4,2) * t73) + t73 * (Ifges(4,1) * t73 + Ifges(4,4) * t72) + t72 * (Ifges(4,4) * t73 + Ifges(4,2) * t72) + t44 * t157 + 0.2e1 * t8 * t26 + t27 * t158 + 0.2e1 * t3 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t1 * t30 + 0.2e1 * t5 * t31 + Ifges(2,3) + (mrSges(5,1) * t157 - 0.2e1 * t25 * mrSges(5,3) + Ifges(5,2) * t49 + t123) * t49 + (mrSges(5,3) * t158 + Ifges(5,1) * t50 - 0.2e1 * Ifges(5,4) * t49 + t163 * t107 + (-t164 - t172) * t104) * t50 + m(3) * (t130 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t54 ^ 2 + t55 ^ 2 + t91 ^ 2) + m(5) * (t25 ^ 2 + t57 ^ 2 + t159) + m(6) * (t5 ^ 2 + t6 ^ 2 + t159) + m(7) * (t1 ^ 2 + t3 ^ 2 + t8 ^ 2) + 0.2e1 * (-t54 * t73 + t55 * t72) * mrSges(4,3) + 0.2e1 * t130 * pkin(7) * mrSges(3,3); t113 + m(6) * (t120 * t60 + t23 * t59) + t29 * t132 + (-t60 * t31 + t121) * t104 + m(7) * (t1 * t52 + t3 * t53 + t56 * t8) + m(5) * (-t23 * t61 + t25 * t62) + Ifges(3,6) * t109 + Ifges(3,5) * t106 + t56 * t26 + t59 * t27 + t52 * t30 + t53 * t28 + (-t106 * mrSges(3,1) - t109 * mrSges(3,2)) * pkin(7) + (-t62 * t49 - t61 * t50) * mrSges(5,3) + (m(4) * (t105 * t55 + t108 * t54) + (t105 * t72 - t108 * t73) * mrSges(4,3)) * pkin(2); t116 + (-t52 * t104 + t53 * t107) * t128 + m(4) * (t105 ^ 2 + t108 ^ 2) * pkin(2) ^ 2 + m(7) * (t52 ^ 2 + t53 ^ 2 + t56 ^ 2) + m(5) * (t61 ^ 2 + t62 ^ 2) + m(6) * (t129 * t60 ^ 2 + t59 ^ 2) + t56 * t156 + t59 * t155 + 0.2e1 * t146 - 0.2e1 * t145 + t60 * t119 + 0.2e1 * t115 + Ifges(3,3); t113 + (m(5) * (t102 * t25 - t103 * t23) + (-t102 * t49 - t103 * t50) * mrSges(5,3)) * pkin(3) + m(6) * (t120 * t88 + t23 * t89) + t29 * t131 + t89 * t27 + t76 * t26 + t64 * t30 + t65 * t28 + (-t88 * t31 + t121) * t104 + m(7) * (t1 * t64 + t3 * t65 + t76 * t8); t146 - t145 + (t89 + t59) * t78 + (t76 + t56) * t77 + t115 + m(6) * (t129 * t88 * t60 + t59 * t89) + m(7) * (t52 * t64 + t53 * t65 + t56 * t76) + (m(5) * (t102 * t62 + t103 * t61) + t118) * pkin(3) + ((t53 + t65) * t107 + (-t52 - t64) * t104) * mrSges(7,3) + t116 + t129 * mrSges(6,3) * (t60 + t88); t76 * t156 + t89 * t155 + (-t64 * t104 + t65 * t107) * t128 + t88 * t119 + m(7) * (t64 ^ 2 + t65 ^ 2 + t76 ^ 2) + m(6) * (t129 * t88 ^ 2 + t89 ^ 2) + t116 + (0.2e1 * t118 + m(5) * (t102 ^ 2 + t103 ^ 2) * pkin(3)) * pkin(3); t49 * mrSges(5,1) + t44 + (t30 + t31) * t107 + (t28 + t29) * t104 + m(6) * (t104 * t6 + t107 * t5) + m(7) * (t1 * t107 + t104 * t3) + m(5) * t57; m(7) * (t104 * t53 + t107 * t52); m(7) * (t104 * t65 + t107 * t64); m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t160; mrSges(6,1) * t5 + mrSges(7,1) * t1 - mrSges(6,2) * t6 - mrSges(7,2) * t3 - t168 * t134 + (m(7) * t1 + t30) * pkin(5) + t123; mrSges(7,1) * t52 - mrSges(7,2) * t53 - t117 * t60 + (m(7) * t52 - t135) * pkin(5) + t122; mrSges(7,1) * t64 - mrSges(7,2) * t65 - t117 * t88 + (m(7) * t64 - t135) * pkin(5) + t122; -t140 - t93 + (mrSges(6,1) + mrSges(7,1) + t154) * t107; (0.2e1 * mrSges(7,1) + t154) * pkin(5) + t165; m(7) * t8 + t26; m(7) * t56 + t77; m(7) * t76 + t77; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
