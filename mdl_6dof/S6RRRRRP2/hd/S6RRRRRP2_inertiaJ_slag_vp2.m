% Calculate joint inertia matrix for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:43
% EndTime: 2019-03-10 01:00:46
% DurationCPUTime: 1.41s
% Computational Cost: add. (2694->268), mult. (5044->375), div. (0->0), fcn. (5489->8), ass. (0->131)
t205 = Ifges(6,1) + Ifges(7,1);
t204 = Ifges(7,4) + Ifges(6,5);
t203 = Ifges(7,2) + Ifges(6,3);
t130 = cos(qJ(5));
t126 = sin(qJ(5));
t171 = Ifges(7,5) * t126;
t173 = Ifges(6,4) * t126;
t127 = sin(qJ(4));
t131 = cos(qJ(4));
t128 = sin(qJ(3));
t129 = sin(qJ(2));
t132 = cos(qJ(3));
t133 = cos(qJ(2));
t83 = -t128 * t129 + t132 * t133;
t84 = t128 * t133 + t129 * t132;
t49 = t127 * t84 - t131 * t83;
t50 = t127 * t83 + t131 * t84;
t202 = (t205 * t130 + t171 - t173) * t50 + t204 * t49;
t170 = Ifges(7,5) * t130;
t172 = Ifges(6,4) * t130;
t201 = t205 * t126 - t170 + t172;
t122 = t126 ^ 2;
t124 = t130 ^ 2;
t200 = t122 + t124;
t199 = (Ifges(6,6) - Ifges(7,6)) * t130 + t204 * t126;
t185 = pkin(3) * t127;
t106 = pkin(10) + t185;
t163 = t106 * t124;
t164 = t106 * t122;
t198 = t163 + t164;
t109 = -pkin(2) * t133 - pkin(1);
t68 = -pkin(3) * t83 + t109;
t22 = pkin(4) * t49 - pkin(10) * t50 + t68;
t191 = -pkin(8) - pkin(7);
t156 = t191 * t129;
t157 = t191 * t133;
t56 = t128 * t157 + t132 * t156;
t142 = -t84 * pkin(9) + t56;
t57 = t128 * t156 - t132 * t157;
t39 = pkin(9) * t83 + t57;
t25 = t127 * t142 + t131 * t39;
t6 = -t126 * t25 + t130 * t22;
t7 = t126 * t22 + t130 * t25;
t151 = -t126 * t6 + t130 * t7;
t108 = pkin(2) * t132 + pkin(3);
t186 = pkin(2) * t128;
t73 = t108 * t131 - t127 * t186;
t87 = -pkin(5) * t130 - qJ(6) * t126 - pkin(4);
t55 = -t73 + t87;
t94 = -mrSges(7,1) * t130 - mrSges(7,3) * t126;
t44 = t55 * t94;
t71 = -pkin(4) - t73;
t95 = -mrSges(6,1) * t130 + mrSges(6,2) * t126;
t51 = t71 * t95;
t74 = t127 * t108 + t131 * t186;
t72 = pkin(10) + t74;
t169 = t122 * t72;
t59 = mrSges(6,3) * t169;
t60 = mrSges(7,2) * t169;
t168 = t124 * t72;
t61 = mrSges(6,3) * t168;
t62 = mrSges(7,2) * t168;
t69 = t73 * mrSges(5,1);
t197 = t44 + t51 + t59 + t60 + t61 + t62 + t69;
t184 = pkin(3) * t131;
t115 = mrSges(5,1) * t184;
t77 = t87 - t184;
t58 = t77 * t94;
t107 = -pkin(4) - t184;
t76 = t107 * t95;
t88 = mrSges(6,3) * t164;
t89 = mrSges(7,2) * t164;
t90 = mrSges(6,3) * t163;
t91 = mrSges(7,2) * t163;
t196 = t115 + t58 + t76 + t88 + t89 + t90 + t91;
t23 = t127 * t39 - t131 * t142;
t195 = t23 ^ 2;
t194 = 0.2e1 * t23;
t193 = 0.2e1 * t68;
t190 = pkin(4) * t95;
t188 = m(7) * t126;
t187 = Ifges(6,6) * t49;
t183 = pkin(10) * t122;
t182 = pkin(10) * t124;
t3 = qJ(6) * t49 + t7;
t180 = t130 * t3;
t178 = t74 * mrSges(5,2);
t177 = t198 * t72;
t176 = (t168 + t169) * pkin(10);
t175 = t200 * t72 ^ 2;
t174 = t198 * pkin(10);
t116 = t126 * mrSges(7,2);
t167 = t126 * t50;
t166 = t130 * t50;
t165 = qJ(6) * t130;
t162 = t200 * t106 ^ 2;
t160 = t200 * pkin(10) ^ 2;
t159 = t129 ^ 2 + t133 ^ 2;
t158 = mrSges(5,2) * t185;
t30 = -t49 * mrSges(7,1) + mrSges(7,2) * t166;
t153 = Ifges(7,6) * t167 + t204 * t166 + t203 * t49;
t4 = -pkin(5) * t49 - t6;
t152 = t126 * t4 + t180;
t150 = t126 * mrSges(6,1) + t130 * mrSges(6,2);
t149 = t126 * mrSges(7,1) - t130 * mrSges(7,3);
t148 = pkin(5) * t126 - t165;
t96 = -Ifges(7,3) * t130 + t171;
t97 = Ifges(6,2) * t130 + t173;
t147 = Ifges(5,3) + (-t96 + t97) * t130 + t201 * t126;
t146 = (mrSges(4,1) * t132 - mrSges(4,2) * t128) * pkin(2);
t145 = Ifges(4,3) + t147;
t28 = -mrSges(6,2) * t49 - mrSges(6,3) * t167;
t29 = mrSges(6,1) * t49 - mrSges(6,3) * t166;
t31 = -mrSges(7,2) * t167 + mrSges(7,3) * t49;
t144 = (t28 + t31) * t130 + (-t29 + t30) * t126;
t143 = mrSges(7,2) * t165 - pkin(5) * t116 + t199;
t111 = mrSges(6,3) * t183;
t112 = mrSges(7,2) * t183;
t113 = mrSges(6,3) * t182;
t114 = mrSges(7,2) * t182;
t67 = t87 * t94;
t141 = t111 + t112 + t113 + t114 + t147 + t67 - t190;
t140 = -m(7) * t148 - t149 - t150;
t14 = Ifges(7,6) * t49 + (Ifges(7,3) * t126 + t170) * t50;
t15 = t187 + (-Ifges(6,2) * t126 + t172) * t50;
t9 = t148 * t50 + t23;
t139 = -t25 * mrSges(5,2) + mrSges(7,2) * t180 + Ifges(5,5) * t50 + t4 * t116 + t9 * t94 + (t95 - mrSges(5,1)) * t23 + t202 * t126 / 0.2e1 + (t96 / 0.2e1 - t97 / 0.2e1) * t167 + t201 * t166 / 0.2e1 + (-t14 / 0.2e1 + t15 / 0.2e1) * t130 + t151 * mrSges(6,3) + (-Ifges(5,6) + t199 / 0.2e1) * t49;
t138 = t56 * mrSges(4,1) - t57 * mrSges(4,2) + Ifges(4,5) * t84 + Ifges(4,6) * t83 + t139;
t27 = t150 * t50;
t26 = t149 * t50;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t133 + mrSges(3,2) * t129) + t133 * (Ifges(3,4) * t129 + Ifges(3,2) * t133) + t129 * (Ifges(3,1) * t129 + Ifges(3,4) * t133) + 0.2e1 * t109 * (-mrSges(4,1) * t83 + mrSges(4,2) * t84) + t84 * (Ifges(4,1) * t84 + Ifges(4,4) * t83) + t83 * (Ifges(4,4) * t84 + Ifges(4,2) * t83) + 0.2e1 * t9 * t26 + t27 * t194 + 0.2e1 * t7 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t4 * t30 + 0.2e1 * t3 * t31 + Ifges(2,3) + (mrSges(5,1) * t193 - 0.2e1 * t25 * mrSges(5,3) + Ifges(5,2) * t49 + t153) * t49 + (mrSges(5,2) * t193 + mrSges(5,3) * t194 + Ifges(5,1) * t50 - 0.2e1 * Ifges(5,4) * t49 + t202 * t130 + (t14 - t15 - t187) * t126) * t50 + m(3) * (t159 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t109 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t25 ^ 2 + t68 ^ 2 + t195) + m(6) * (t6 ^ 2 + t7 ^ 2 + t195) + m(7) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + 0.2e1 * (-t56 * t84 + t57 * t83) * mrSges(4,3) + 0.2e1 * t159 * pkin(7) * mrSges(3,3); t138 + Ifges(3,6) * t133 + Ifges(3,5) * t129 + t71 * t27 + t55 * t26 + m(6) * (t151 * t72 + t23 * t71) + m(7) * (t152 * t72 + t55 * t9) + t144 * t72 + m(5) * (-t23 * t73 + t25 * t74) + (-t129 * mrSges(3,1) - t133 * mrSges(3,2)) * pkin(7) + (-t74 * t49 - t73 * t50) * mrSges(5,3) + (m(4) * (t128 * t57 + t132 * t56) + (t128 * t83 - t132 * t84) * mrSges(4,3)) * pkin(2); m(7) * (t55 ^ 2 + t175) + m(6) * (t71 ^ 2 + t175) + m(5) * (t73 ^ 2 + t74 ^ 2) + m(4) * (t128 ^ 2 + t132 ^ 2) * pkin(2) ^ 2 + t145 + 0.2e1 * t69 + 0.2e1 * t62 + 0.2e1 * t60 + 0.2e1 * t61 + 0.2e1 * t59 + 0.2e1 * t51 + 0.2e1 * t44 - 0.2e1 * t178 + 0.2e1 * t146 + Ifges(3,3); t138 + t107 * t27 + t77 * t26 + (m(5) * (t127 * t25 - t131 * t23) + (-t127 * t49 - t131 * t50) * mrSges(5,3)) * pkin(3) + m(7) * (t152 * t106 + t77 * t9) + m(6) * (t151 * t106 + t107 * t23) + t144 * t106; m(5) * (t127 * t74 + t131 * t73) * pkin(3) + m(6) * (t107 * t71 + t177) + m(7) * (t55 * t77 + t177) + t146 + t145 + t196 + (-t74 - t185) * mrSges(5,2) + t197; -0.2e1 * t158 + 0.2e1 * t115 + 0.2e1 * t58 + 0.2e1 * t76 + 0.2e1 * t88 + 0.2e1 * t89 + 0.2e1 * t90 + 0.2e1 * t91 + m(7) * (t77 ^ 2 + t162) + m(6) * (t107 ^ 2 + t162) + m(5) * (t127 ^ 2 + t131 ^ 2) * pkin(3) ^ 2 + t145; t139 + t87 * t26 - pkin(4) * t27 + m(7) * (t152 * pkin(10) + t87 * t9) + m(6) * (-pkin(4) * t23 + t151 * pkin(10)) + t144 * pkin(10); m(7) * (t55 * t87 + t176) + m(6) * (-pkin(4) * t71 + t176) + t141 - t178 + t197; m(7) * (t77 * t87 + t174) + m(6) * (-pkin(4) * t107 + t174) + t141 - t158 + t196; -0.2e1 * t190 + 0.2e1 * t111 + 0.2e1 * t112 + 0.2e1 * t113 + 0.2e1 * t114 + 0.2e1 * t67 + m(7) * (t87 ^ 2 + t160) + m(6) * (pkin(4) ^ 2 + t160) + t147; qJ(6) * t31 + t3 * mrSges(7,3) - t7 * mrSges(6,2) - t4 * mrSges(7,1) - pkin(5) * t30 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) - Ifges(6,6) * t167 + t6 * mrSges(6,1) + t153; t140 * t72 + t143; t140 * t106 + t143; t140 * pkin(10) + t143; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t203; m(7) * t4 + t30; t72 * t188 + t116; t106 * t188 + t116; pkin(10) * t188 + t116; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
