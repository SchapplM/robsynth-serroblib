% Calculate joint inertia matrix for
% S6RRRRRP1
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:37
% EndTime: 2019-03-10 00:55:41
% DurationCPUTime: 1.82s
% Computational Cost: add. (2680->289), mult. (4993->397), div. (0->0), fcn. (5464->8), ass. (0->123)
t196 = Ifges(6,4) + Ifges(7,4);
t191 = Ifges(6,6) + Ifges(7,6);
t120 = sin(qJ(4));
t124 = cos(qJ(4));
t121 = sin(qJ(3));
t122 = sin(qJ(2));
t125 = cos(qJ(3));
t126 = cos(qJ(2));
t82 = -t121 * t122 + t125 * t126;
t83 = t121 * t126 + t122 * t125;
t48 = t120 * t83 - t124 * t82;
t195 = t191 * t48;
t194 = Ifges(6,1) + Ifges(7,1);
t193 = Ifges(6,5) + Ifges(7,5);
t192 = Ifges(6,2) + Ifges(7,2);
t123 = cos(qJ(5));
t190 = t196 * t123;
t119 = sin(qJ(5));
t189 = t196 * t119;
t188 = Ifges(6,3) + Ifges(7,3);
t49 = t120 * t82 + t124 * t83;
t187 = (-t192 * t119 + t190) * t49 + t195;
t186 = (t194 * t123 - t189) * t49 + t193 * t48;
t185 = t192 * t123 + t189;
t184 = t194 * t119 + t190;
t141 = t193 * t119 + t191 * t123;
t170 = pkin(5) * t123;
t103 = pkin(2) * t125 + pkin(3);
t173 = pkin(2) * t121;
t68 = t103 * t124 - t120 * t173;
t66 = -pkin(4) - t68;
t60 = t66 - t170;
t90 = -t123 * mrSges(7,1) + t119 * mrSges(7,2);
t50 = t60 * t90;
t163 = mrSges(7,3) * t123;
t109 = t123 * qJ(6);
t69 = t120 * t103 + t124 * t173;
t67 = pkin(10) + t69;
t155 = t123 * t67;
t55 = t109 + t155;
t51 = t55 * t163;
t91 = -mrSges(6,1) * t123 + mrSges(6,2) * t119;
t52 = t66 * t91;
t115 = t119 ^ 2;
t165 = mrSges(6,3) * t115;
t58 = t67 * t165;
t117 = t123 ^ 2;
t164 = mrSges(6,3) * t117;
t59 = t67 * t164;
t64 = t68 * mrSges(5,1);
t183 = t50 + t51 + t52 + t58 + t59 + t64;
t171 = pkin(3) * t124;
t108 = mrSges(5,1) * t171;
t104 = -pkin(4) - t170;
t88 = t104 - t171;
t61 = t88 * t90;
t172 = pkin(3) * t120;
t101 = pkin(10) + t172;
t154 = t101 * t123;
t74 = t109 + t154;
t63 = t74 * t163;
t102 = -pkin(4) - t171;
t71 = t102 * t91;
t86 = t101 * t165;
t87 = t101 * t164;
t182 = t108 + t61 + t63 + t71 + t86 + t87;
t177 = -pkin(8) - pkin(7);
t146 = t177 * t122;
t147 = t177 * t126;
t56 = t121 * t147 + t125 * t146;
t133 = -t83 * pkin(9) + t56;
t57 = t121 * t146 - t125 * t147;
t39 = pkin(9) * t82 + t57;
t23 = t120 * t39 - t124 * t133;
t181 = t23 ^ 2;
t180 = 0.2e1 * t23;
t105 = -pkin(2) * t126 - pkin(1);
t62 = -pkin(3) * t82 + t105;
t179 = 0.2e1 * t62;
t176 = pkin(4) * t91;
t169 = pkin(10) * t123;
t22 = pkin(4) * t48 - pkin(10) * t49 + t62;
t25 = t120 * t133 + t124 * t39;
t6 = t119 * t22 + t123 * t25;
t168 = t123 * t6;
t167 = t69 * mrSges(5,2);
t156 = t123 * t49;
t157 = t119 * t49;
t26 = mrSges(7,1) * t157 + mrSges(7,2) * t156;
t158 = t119 * mrSges(7,3);
t151 = t115 + t117;
t150 = t122 ^ 2 + t126 ^ 2;
t149 = -0.2e1 * t158;
t148 = mrSges(5,2) * t172;
t5 = -t119 * t25 + t123 * t22;
t143 = t193 * t156 + t188 * t48;
t142 = t151 * t101;
t140 = t184 * t119 + t123 * t185 + Ifges(5,3);
t1 = pkin(5) * t48 - t49 * t109 + t5;
t139 = -t5 * mrSges(6,3) - t1 * mrSges(7,3);
t138 = -t119 * t5 + t168;
t137 = mrSges(6,1) * t119 + mrSges(6,2) * t123;
t136 = Ifges(4,3) + t140;
t135 = (mrSges(4,1) * t125 - mrSges(4,2) * t121) * pkin(2);
t106 = pkin(10) * t165;
t107 = pkin(10) * t164;
t72 = t104 * t90;
t92 = t109 + t169;
t77 = t92 * t163;
t134 = t106 + t107 + t140 + t72 + t77 - t176;
t3 = -qJ(6) * t157 + t6;
t8 = pkin(5) * t157 + t23;
t132 = -t25 * mrSges(5,2) + mrSges(6,3) * t168 + Ifges(5,5) * t49 + t3 * t163 + t8 * t90 + (-mrSges(5,1) + t91) * t23 + t186 * t119 / 0.2e1 + t187 * t123 / 0.2e1 - t185 * t157 / 0.2e1 + t184 * t156 / 0.2e1 + (-Ifges(5,6) + t141 / 0.2e1) * t48;
t131 = t56 * mrSges(4,1) - t57 * mrSges(4,2) + Ifges(4,5) * t83 + Ifges(4,6) * t82 + t132;
t89 = (-qJ(6) - pkin(10)) * t119;
t73 = (-qJ(6) - t101) * t119;
t54 = (-qJ(6) - t67) * t119;
t31 = mrSges(6,1) * t48 - mrSges(6,3) * t156;
t30 = mrSges(7,1) * t48 - mrSges(7,3) * t156;
t29 = -mrSges(6,2) * t48 - mrSges(6,3) * t157;
t28 = -mrSges(7,2) * t48 - mrSges(7,3) * t157;
t27 = t137 * t49;
t2 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t126 + mrSges(3,2) * t122) + t122 * (Ifges(3,1) * t122 + Ifges(3,4) * t126) + t126 * (Ifges(3,4) * t122 + Ifges(3,2) * t126) + 0.2e1 * t105 * (-mrSges(4,1) * t82 + mrSges(4,2) * t83) + t83 * (Ifges(4,1) * t83 + Ifges(4,4) * t82) + t82 * (Ifges(4,4) * t83 + Ifges(4,2) * t82) + 0.2e1 * t8 * t26 + t27 * t180 + 0.2e1 * t3 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t1 * t30 + 0.2e1 * t5 * t31 + Ifges(2,3) + (mrSges(5,1) * t179 - 0.2e1 * t25 * mrSges(5,3) + Ifges(5,2) * t48 + t143) * t48 + (mrSges(5,2) * t179 + mrSges(5,3) * t180 + Ifges(5,1) * t49 - 0.2e1 * Ifges(5,4) * t48 + t186 * t123 + (-t187 - t195) * t119) * t49 + m(3) * (t150 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t105 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t25 ^ 2 + t62 ^ 2 + t181) + m(6) * (t5 ^ 2 + t6 ^ 2 + t181) + m(7) * (t1 ^ 2 + t3 ^ 2 + t8 ^ 2) + 0.2e1 * (-t56 * t83 + t57 * t82) * mrSges(4,3) + 0.2e1 * t150 * pkin(7) * mrSges(3,3); t131 + t29 * t155 + Ifges(3,6) * t126 + Ifges(3,5) * t122 + t60 * t26 + t66 * t27 + t54 * t30 + t55 * t28 + m(6) * (t138 * t67 + t23 * t66) + (-t67 * t31 + t139) * t119 + m(5) * (-t23 * t68 + t25 * t69) + m(7) * (t1 * t54 + t3 * t55 + t60 * t8) + (-t122 * mrSges(3,1) - t126 * mrSges(3,2)) * pkin(7) + (-t69 * t48 - t68 * t49) * mrSges(5,3) + (m(4) * (t121 * t57 + t125 * t56) + (t121 * t82 - t125 * t83) * mrSges(4,3)) * pkin(2); t136 + m(4) * (t121 ^ 2 + t125 ^ 2) * pkin(2) ^ 2 + m(7) * (t54 ^ 2 + t55 ^ 2 + t60 ^ 2) + m(5) * (t68 ^ 2 + t69 ^ 2) + m(6) * (t151 * t67 ^ 2 + t66 ^ 2) + 0.2e1 * t64 + 0.2e1 * t58 + 0.2e1 * t59 + 0.2e1 * t51 + 0.2e1 * t52 + 0.2e1 * t50 - 0.2e1 * t167 + t54 * t149 + 0.2e1 * t135 + Ifges(3,3); t131 + (m(5) * (t120 * t25 - t124 * t23) + (-t120 * t48 - t124 * t49) * mrSges(5,3)) * pkin(3) + t29 * t154 + t102 * t27 + t74 * t28 + t88 * t26 + t73 * t30 + m(6) * (t138 * t101 + t102 * t23) + (-t101 * t31 + t139) * t119 + m(7) * (t1 * t73 + t3 * t74 + t8 * t88); t182 + t136 + t135 + m(7) * (t54 * t73 + t55 * t74 + t60 * t88) + m(5) * (t120 * t69 + t124 * t68) * pkin(3) + (-t69 - t172) * mrSges(5,2) + m(6) * (t102 * t66 + t67 * t142) + (-t54 - t73) * t158 + t183; -0.2e1 * t148 + t73 * t149 + 0.2e1 * t108 + 0.2e1 * t61 + 0.2e1 * t63 + 0.2e1 * t71 + 0.2e1 * t86 + 0.2e1 * t87 + m(7) * (t73 ^ 2 + t74 ^ 2 + t88 ^ 2) + m(6) * (t151 * t101 ^ 2 + t102 ^ 2) + m(5) * (t120 ^ 2 + t124 ^ 2) * pkin(3) ^ 2 + t136; t132 + (-pkin(10) * t31 + t139) * t119 + m(7) * (t1 * t89 + t104 * t8 + t3 * t92) + t29 * t169 + t89 * t30 + t92 * t28 + t104 * t26 - pkin(4) * t27 + m(6) * (-pkin(4) * t23 + t138 * pkin(10)); t134 + m(7) * (t104 * t60 + t54 * t89 + t55 * t92) + (-t54 - t89) * t158 - t167 + m(6) * (t151 * t67 * pkin(10) - pkin(4) * t66) + t183; t134 + m(7) * (t104 * t88 + t73 * t89 + t74 * t92) + (-t73 - t89) * t158 + m(6) * (-pkin(4) * t102 + pkin(10) * t142) - t148 + t182; t89 * t149 - 0.2e1 * t176 + 0.2e1 * t106 + 0.2e1 * t107 + 0.2e1 * t72 + 0.2e1 * t77 + m(7) * (t104 ^ 2 + t89 ^ 2 + t92 ^ 2) + m(6) * (t151 * pkin(10) ^ 2 + pkin(4) ^ 2) + t140; mrSges(6,1) * t5 + mrSges(7,1) * t1 - mrSges(6,2) * t6 - mrSges(7,2) * t3 - t191 * t157 + (m(7) * t1 + t30) * pkin(5) + t143; mrSges(7,1) * t54 - mrSges(7,2) * t55 - t137 * t67 + (m(7) * t54 - t158) * pkin(5) + t141; mrSges(7,1) * t73 - mrSges(7,2) * t74 - t137 * t101 + (m(7) * t73 - t158) * pkin(5) + t141; mrSges(7,1) * t89 - mrSges(7,2) * t92 - t137 * pkin(10) + (m(7) * t89 - t158) * pkin(5) + t141; (m(7) * pkin(5) + 0.2e1 * mrSges(7,1)) * pkin(5) + t188; m(7) * t8 + t26; m(7) * t60 + t90; m(7) * t88 + t90; m(7) * t104 + t90; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
