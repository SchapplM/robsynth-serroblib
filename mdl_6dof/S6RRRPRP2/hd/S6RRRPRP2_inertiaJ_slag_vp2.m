% Calculate joint inertia matrix for
% S6RRRPRP2
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

function Mq = S6RRRPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:41:47
% EndTime: 2018-11-23 17:41:48
% DurationCPUTime: 1.20s
% Computational Cost: add. (2246->245), mult. (4195->335), div. (0->0), fcn. (4657->8), ass. (0->103)
t177 = Ifges(6,1) + Ifges(7,1);
t176 = Ifges(7,4) + Ifges(6,5);
t175 = Ifges(7,2) + Ifges(6,3);
t107 = cos(qJ(5));
t104 = sin(qJ(5));
t142 = Ifges(7,5) * t104;
t144 = Ifges(6,4) * t104;
t102 = sin(pkin(10));
t103 = cos(pkin(10));
t105 = sin(qJ(3));
t106 = sin(qJ(2));
t108 = cos(qJ(3));
t109 = cos(qJ(2));
t72 = -t105 * t106 + t108 * t109;
t73 = t105 * t109 + t106 * t108;
t49 = t102 * t73 - t103 * t72;
t50 = t102 * t72 + t103 * t73;
t174 = (t177 * t107 + t142 - t144) * t50 + t176 * t49;
t141 = Ifges(7,5) * t107;
t143 = Ifges(6,4) * t107;
t173 = t177 * t104 - t141 + t143;
t137 = t104 ^ 2 + t107 ^ 2;
t172 = (Ifges(6,6) - Ifges(7,6)) * t107 + t176 * t104;
t92 = -pkin(2) * t109 - pkin(1);
t59 = -pkin(3) * t72 + t92;
t22 = pkin(4) * t49 - pkin(9) * t50 + t59;
t163 = -pkin(8) - pkin(7);
t134 = t163 * t106;
t135 = t163 * t109;
t55 = t105 * t135 + t108 * t134;
t115 = -t73 * qJ(4) + t55;
t56 = t105 * t134 - t108 * t135;
t39 = qJ(4) * t72 + t56;
t25 = t102 * t115 + t103 * t39;
t6 = -t104 * t25 + t107 * t22;
t7 = t104 * t22 + t107 * t25;
t128 = -t104 * t6 + t107 * t7;
t171 = 0.2e1 * t137;
t170 = mrSges(6,3) + mrSges(7,2);
t23 = t102 * t39 - t103 * t115;
t169 = t23 ^ 2;
t168 = 0.2e1 * t23;
t167 = 0.2e1 * t59;
t76 = -t107 * mrSges(7,1) - t104 * mrSges(7,3);
t166 = 0.2e1 * t76;
t77 = -t107 * mrSges(6,1) + t104 * mrSges(6,2);
t165 = 0.2e1 * t77;
t161 = m(7) * t104;
t160 = Ifges(6,6) * t49;
t159 = pkin(2) * t105;
t158 = pkin(3) * t103;
t3 = qJ(6) * t49 + t7;
t156 = t107 * t3;
t91 = pkin(2) * t108 + pkin(3);
t63 = -t102 * t159 + t103 * t91;
t153 = t63 * mrSges(5,1);
t64 = t102 * t91 + t103 * t159;
t152 = t64 * mrSges(5,2);
t140 = t104 * t50;
t28 = -mrSges(6,2) * t49 - mrSges(6,3) * t140;
t31 = -mrSges(7,2) * t140 + mrSges(7,3) * t49;
t151 = t28 + t31;
t139 = t107 * t50;
t29 = mrSges(6,1) * t49 - mrSges(6,3) * t139;
t30 = -t49 * mrSges(7,1) + mrSges(7,2) * t139;
t150 = t29 - t30;
t62 = pkin(9) + t64;
t89 = pkin(3) * t102 + pkin(9);
t149 = t137 * t62 * t89;
t148 = t137 * t62 ^ 2;
t146 = t137 * t89 ^ 2;
t94 = t104 * mrSges(7,2);
t138 = t106 ^ 2 + t109 ^ 2;
t136 = qJ(6) * t107;
t130 = Ifges(7,6) * t140 + t176 * t139 + t175 * t49;
t4 = -pkin(5) * t49 - t6;
t129 = t104 * t4 + t156;
t127 = t103 * mrSges(5,1) - t102 * mrSges(5,2);
t126 = t104 * mrSges(6,1) + t107 * mrSges(6,2);
t125 = t104 * mrSges(7,1) - t107 * mrSges(7,3);
t124 = pkin(5) * t107 + qJ(6) * t104;
t123 = pkin(5) * t104 - t136;
t122 = -pkin(4) - t124;
t121 = (mrSges(4,1) * t108 - mrSges(4,2) * t105) * pkin(2);
t78 = -Ifges(7,3) * t107 + t142;
t79 = Ifges(6,2) * t107 + t144;
t120 = Ifges(4,3) + Ifges(5,3) + (-t78 + t79) * t107 + t173 * t104;
t118 = -t150 * t104 + t151 * t107;
t117 = mrSges(7,2) * t136 - pkin(5) * t94 + t172;
t116 = t170 * t171;
t114 = -m(7) * t123 - t125 - t126;
t14 = Ifges(7,6) * t49 + (Ifges(7,3) * t104 + t141) * t50;
t15 = t160 + (-Ifges(6,2) * t104 + t143) * t50;
t9 = t123 * t50 + t23;
t113 = t55 * mrSges(4,1) - t56 * mrSges(4,2) - t25 * mrSges(5,2) + mrSges(7,2) * t156 + Ifges(4,5) * t73 + Ifges(5,5) * t50 + Ifges(4,6) * t72 + t4 * t94 + t9 * t76 + (t77 - mrSges(5,1)) * t23 + t174 * t104 / 0.2e1 + (t78 / 0.2e1 - t79 / 0.2e1) * t140 + t173 * t139 / 0.2e1 + (-t14 / 0.2e1 + t15 / 0.2e1) * t107 + t128 * mrSges(6,3) + (-Ifges(5,6) + t172 / 0.2e1) * t49;
t90 = -pkin(4) - t158;
t66 = t122 - t158;
t61 = -pkin(4) - t63;
t54 = t122 - t63;
t44 = t50 * mrSges(5,2);
t27 = t126 * t50;
t26 = t125 * t50;
t1 = [t106 * (Ifges(3,1) * t106 + Ifges(3,4) * t109) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t109 + mrSges(3,2) * t106) + t109 * (Ifges(3,4) * t106 + Ifges(3,2) * t109) + 0.2e1 * t92 * (-mrSges(4,1) * t72 + mrSges(4,2) * t73) + t73 * (Ifges(4,1) * t73 + Ifges(4,4) * t72) + t72 * (Ifges(4,4) * t73 + Ifges(4,2) * t72) + t44 * t167 + 0.2e1 * t4 * t30 + 0.2e1 * t3 * t31 + 0.2e1 * t9 * t26 + t27 * t168 + 0.2e1 * t7 * t28 + 0.2e1 * t6 * t29 + Ifges(2,3) + (mrSges(5,1) * t167 - 0.2e1 * t25 * mrSges(5,3) + Ifges(5,2) * t49 + t130) * t49 + (mrSges(5,3) * t168 + Ifges(5,1) * t50 - 0.2e1 * Ifges(5,4) * t49 + t174 * t107 + (t14 - t15 - t160) * t104) * t50 + m(3) * (t138 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t55 ^ 2 + t56 ^ 2 + t92 ^ 2) + m(5) * (t25 ^ 2 + t59 ^ 2 + t169) + m(7) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2 + t169) + 0.2e1 * (-t55 * t73 + t56 * t72) * mrSges(4,3) + 0.2e1 * t138 * pkin(7) * mrSges(3,3); m(7) * (t129 * t62 + t54 * t9) + m(6) * (t128 * t62 + t23 * t61) + t113 + t118 * t62 + Ifges(3,5) * t106 + Ifges(3,6) * t109 + t61 * t27 + t54 * t26 + m(5) * (-t23 * t63 + t25 * t64) + (-t106 * mrSges(3,1) - t109 * mrSges(3,2)) * pkin(7) + (-t64 * t49 - t63 * t50) * mrSges(5,3) + (m(4) * (t105 * t56 + t108 * t55) + (t105 * t72 - t108 * t73) * mrSges(4,3)) * pkin(2); 0.2e1 * t153 - 0.2e1 * t152 + t54 * t166 + t61 * t165 + Ifges(3,3) + 0.2e1 * t121 + t116 * t62 + m(7) * (t54 ^ 2 + t148) + m(6) * (t61 ^ 2 + t148) + m(5) * (t63 ^ 2 + t64 ^ 2) + m(4) * (t105 ^ 2 + t108 ^ 2) * pkin(2) ^ 2 + t120; (m(5) * (t102 * t25 - t103 * t23) + (-t102 * t49 - t103 * t50) * mrSges(5,3)) * pkin(3) + m(7) * (t129 * t89 + t66 * t9) + m(6) * (t128 * t89 + t23 * t90) + t113 + t90 * t27 + t66 * t26 + t118 * t89; t153 - t152 + (t90 + t61) * t77 + (t66 + t54) * t76 + t121 + m(6) * (t61 * t90 + t149) + m(7) * (t54 * t66 + t149) + (m(5) * (t102 * t64 + t103 * t63) + t127) * pkin(3) + t120 + t170 * t137 * (t62 + t89); t66 * t166 + t90 * t165 + m(7) * (t66 ^ 2 + t146) + m(6) * (t90 ^ 2 + t146) + t116 * t89 + t120 + (0.2e1 * t127 + m(5) * (t102 ^ 2 + t103 ^ 2) * pkin(3)) * pkin(3); t49 * mrSges(5,1) + t44 + t150 * t107 + t151 * t104 + m(6) * (t104 * t7 + t107 * t6) + m(7) * (t104 * t3 - t107 * t4) + m(5) * t59; 0; 0; m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t171; -Ifges(6,6) * t140 - pkin(5) * t30 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) + qJ(6) * t31 + t3 * mrSges(7,3) - t7 * mrSges(6,2) - t4 * mrSges(7,1) + t6 * mrSges(6,1) + t130; t114 * t62 + t117; t114 * t89 + t117; m(7) * t124 - t76 - t77; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t175; m(7) * t4 + t30; t62 * t161 + t94; t89 * t161 + t94; -m(7) * t107; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
