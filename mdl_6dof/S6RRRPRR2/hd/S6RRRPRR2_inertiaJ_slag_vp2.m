% Calculate joint inertia matrix for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:20
% EndTime: 2019-03-09 18:06:22
% DurationCPUTime: 1.13s
% Computational Cost: add. (3347->290), mult. (6236->422), div. (0->0), fcn. (7245->10), ass. (0->113)
t121 = cos(qJ(5));
t114 = sin(pkin(11));
t115 = cos(pkin(11));
t118 = sin(qJ(3));
t119 = sin(qJ(2));
t122 = cos(qJ(3));
t123 = cos(qJ(2));
t90 = -t118 * t119 + t122 * t123;
t92 = t118 * t123 + t119 * t122;
t62 = t114 * t90 + t115 * t92;
t150 = t121 * t62;
t61 = t114 * t92 - t115 * t90;
t172 = Ifges(6,5) * t150 + Ifges(6,3) * t61;
t116 = sin(qJ(6));
t117 = sin(qJ(5));
t148 = Ifges(6,5) * t117 + Ifges(6,6) * t121;
t120 = cos(qJ(6));
t89 = -t116 * t117 + t120 * t121;
t158 = t89 * mrSges(7,3);
t171 = t116 * pkin(5) * t158 + t148;
t164 = -pkin(8) - pkin(7);
t141 = t164 * t119;
t142 = t164 * t123;
t69 = t118 * t142 + t122 * t141;
t129 = -t92 * qJ(4) + t69;
t70 = t118 * t141 - t122 * t142;
t51 = qJ(4) * t90 + t70;
t31 = t114 * t51 - t115 * t129;
t170 = t31 ^ 2;
t169 = 0.2e1 * t31;
t91 = t116 * t121 + t117 * t120;
t63 = -t89 * mrSges(7,1) + t91 * mrSges(7,2);
t168 = 0.2e1 * t63;
t106 = -pkin(2) * t123 - pkin(1);
t72 = -pkin(3) * t90 + t106;
t167 = 0.2e1 * t72;
t96 = -t121 * mrSges(6,1) + t117 * mrSges(6,2);
t166 = 0.2e1 * t96;
t163 = Ifges(6,6) * t61;
t162 = pkin(2) * t118;
t161 = pkin(5) * t121;
t105 = pkin(2) * t122 + pkin(3);
t77 = t105 * t115 - t114 * t162;
t160 = t77 * mrSges(5,1);
t78 = t114 * t105 + t115 * t162;
t159 = t78 * mrSges(5,2);
t157 = t91 * mrSges(7,3);
t30 = pkin(4) * t61 - pkin(9) * t62 + t72;
t33 = t114 * t129 + t115 * t51;
t13 = t117 * t30 + t121 * t33;
t156 = Ifges(7,5) * t91 + Ifges(7,6) * t89;
t155 = Ifges(6,4) * t117;
t154 = Ifges(6,4) * t121;
t12 = -t117 * t33 + t121 * t30;
t153 = t117 * t12;
t152 = t117 * t62;
t151 = t121 * t13;
t75 = pkin(9) + t78;
t149 = t121 * t75;
t147 = t117 ^ 2 + t121 ^ 2;
t146 = t119 ^ 2 + t123 ^ 2;
t145 = 0.2e1 * mrSges(7,3);
t36 = t91 * t62;
t37 = t89 * t62;
t144 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t61;
t143 = t120 * t157;
t104 = -t115 * pkin(3) - pkin(4);
t103 = t114 * pkin(3) + pkin(9);
t140 = t147 * t103;
t74 = -pkin(4) - t77;
t139 = t115 * mrSges(5,1) - t114 * mrSges(5,2);
t138 = t117 * mrSges(6,1) + t121 * mrSges(6,2);
t137 = t151 - t153;
t67 = (-pkin(10) - t75) * t117;
t109 = t121 * pkin(10);
t68 = t109 + t149;
t43 = -t116 * t68 + t120 * t67;
t44 = t116 * t67 + t120 * t68;
t136 = t43 * mrSges(7,1) - t44 * mrSges(7,2) + t156;
t80 = (-pkin(10) - t103) * t117;
t81 = t103 * t121 + t109;
t54 = -t116 * t81 + t120 * t80;
t55 = t116 * t80 + t120 * t81;
t135 = t54 * mrSges(7,1) - t55 * mrSges(7,2) + t156;
t134 = 0.2e1 * mrSges(6,3) * t147;
t64 = Ifges(7,4) * t91 + Ifges(7,2) * t89;
t65 = Ifges(7,1) * t91 + Ifges(7,4) * t89;
t97 = Ifges(6,2) * t121 + t155;
t98 = Ifges(6,1) * t117 + t154;
t133 = t117 * t98 + t121 * t97 + t89 * t64 + t91 * t65 + Ifges(4,3) + Ifges(5,3);
t5 = pkin(5) * t61 - pkin(10) * t150 + t12;
t6 = -pkin(10) * t152 + t13;
t3 = -t116 * t6 + t120 * t5;
t4 = t116 * t5 + t120 * t6;
t132 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t144;
t131 = (t122 * mrSges(4,1) - t118 * mrSges(4,2)) * pkin(2);
t130 = (mrSges(7,1) * t120 - mrSges(7,2) * t116) * pkin(5);
t10 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t61;
t16 = pkin(5) * t152 + t31;
t21 = t163 + (-Ifges(6,2) * t117 + t154) * t62;
t22 = Ifges(6,5) * t61 + (Ifges(6,1) * t121 - t155) * t62;
t9 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t61;
t128 = -t70 * mrSges(4,2) - t33 * mrSges(5,2) - t3 * t157 + t4 * t158 + mrSges(6,3) * t151 + t16 * t63 + t117 * t22 / 0.2e1 + t121 * t21 / 0.2e1 - t36 * t64 / 0.2e1 + t37 * t65 / 0.2e1 - t97 * t152 / 0.2e1 + t98 * t150 / 0.2e1 - Ifges(5,6) * t61 + Ifges(5,5) * t62 + t69 * mrSges(4,1) + t89 * t9 / 0.2e1 + t91 * t10 / 0.2e1 + Ifges(4,6) * t90 + Ifges(4,5) * t92 + (t96 - mrSges(5,1)) * t31 + (t156 + t148) * t61 / 0.2e1;
t95 = t104 - t161;
t71 = t74 - t161;
t56 = t62 * mrSges(5,2);
t40 = mrSges(6,1) * t61 - mrSges(6,3) * t150;
t39 = -mrSges(6,2) * t61 - mrSges(6,3) * t152;
t38 = t138 * t62;
t18 = mrSges(7,1) * t61 - mrSges(7,3) * t37;
t17 = -mrSges(7,2) * t61 - mrSges(7,3) * t36;
t14 = mrSges(7,1) * t36 + mrSges(7,2) * t37;
t1 = [t119 * (Ifges(3,1) * t119 + Ifges(3,4) * t123) + t123 * (Ifges(3,4) * t119 + Ifges(3,2) * t123) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t123 + mrSges(3,2) * t119) + 0.2e1 * t106 * (-mrSges(4,1) * t90 + mrSges(4,2) * t92) + (mrSges(5,1) * t167 - 0.2e1 * t33 * mrSges(5,3) + Ifges(5,2) * t61 + t144 + t172) * t61 + t92 * (Ifges(4,1) * t92 + Ifges(4,4) * t90) + t90 * (Ifges(4,4) * t92 + Ifges(4,2) * t90) + 0.2e1 * (-t69 * t92 + t70 * t90) * mrSges(4,3) + 0.2e1 * t146 * pkin(7) * mrSges(3,3) - t36 * t9 + t37 * t10 + 0.2e1 * t13 * t39 + 0.2e1 * t12 * t40 + 0.2e1 * t4 * t17 + 0.2e1 * t3 * t18 + 0.2e1 * t16 * t14 + m(4) * (t106 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(7) * (t16 ^ 2 + t3 ^ 2 + t4 ^ 2) + (mrSges(5,3) * t169 + Ifges(5,1) * t62 - 0.2e1 * Ifges(5,4) * t61 + t121 * t22 + (-t21 - t163) * t117) * t62 + m(5) * (t33 ^ 2 + t72 ^ 2 + t170) + m(6) * (t12 ^ 2 + t13 ^ 2 + t170) + m(3) * (pkin(7) ^ 2 * t146 + pkin(1) ^ 2) + Ifges(2,3) + t56 * t167 + t38 * t169; t39 * t149 + Ifges(3,6) * t123 + Ifges(3,5) * t119 + t71 * t14 + t74 * t38 + t128 + m(6) * (t137 * t75 + t31 * t74) + (-t12 * mrSges(6,3) - t75 * t40) * t117 + m(7) * (t16 * t71 + t3 * t43 + t4 * t44) + m(5) * (-t31 * t77 + t33 * t78) + t43 * t18 + t44 * t17 + (-t119 * mrSges(3,1) - t123 * mrSges(3,2)) * pkin(7) + (-t78 * t61 - t77 * t62) * mrSges(5,3) + (m(4) * (t118 * t70 + t122 * t69) + (t118 * t90 - t122 * t92) * mrSges(4,3)) * pkin(2); m(4) * (t118 ^ 2 + t122 ^ 2) * pkin(2) ^ 2 + t74 * t166 + t71 * t168 + 0.2e1 * t160 - 0.2e1 * t159 + m(6) * (t147 * t75 ^ 2 + t74 ^ 2) + (-t43 * t91 + t44 * t89) * t145 + t133 + t75 * t134 + 0.2e1 * t131 + Ifges(3,3) + m(7) * (t43 ^ 2 + t44 ^ 2 + t71 ^ 2) + m(5) * (t77 ^ 2 + t78 ^ 2); -mrSges(6,3) * t153 + t95 * t14 + t104 * t38 + t128 + (m(5) * (t114 * t33 - t115 * t31) + (-t114 * t61 - t115 * t62) * mrSges(5,3)) * pkin(3) + m(6) * (t103 * t137 + t104 * t31) + t54 * t18 + t55 * t17 + (-t117 * t40 + t121 * t39) * t103 + m(7) * (t16 * t95 + t3 * t54 + t4 * t55); t160 - t159 + (t74 + t104) * t96 + (t95 + t71) * t63 + t131 + m(6) * (t104 * t74 + t140 * t75) + m(7) * (t43 * t54 + t44 * t55 + t71 * t95) + (m(5) * (t114 * t78 + t115 * t77) + t139) * pkin(3) + ((-t43 - t54) * t91 + (t44 + t55) * t89) * mrSges(7,3) + (t147 * t75 + t140) * mrSges(6,3) + t133; t104 * t166 + t95 * t168 + (-t54 * t91 + t55 * t89) * t145 + t103 * t134 + m(7) * (t54 ^ 2 + t55 ^ 2 + t95 ^ 2) + m(6) * (t103 ^ 2 * t147 + t104 ^ 2) + t133 + (0.2e1 * t139 + m(5) * (t114 ^ 2 + t115 ^ 2) * pkin(3)) * pkin(3); t61 * mrSges(5,1) + t117 * t39 + t121 * t40 + t91 * t17 + t89 * t18 + t56 + m(7) * (t3 * t89 + t4 * t91) + m(6) * (t117 * t13 + t12 * t121) + m(5) * t72; m(7) * (t43 * t89 + t44 * t91); m(7) * (t54 * t89 + t55 * t91); m(5) + m(6) * t147 + m(7) * (t89 ^ 2 + t91 ^ 2); -Ifges(6,6) * t152 + t12 * mrSges(6,1) - t13 * mrSges(6,2) + (t120 * t18 + m(7) * (t116 * t4 + t120 * t3) + t116 * t17) * pkin(5) + t132 + t172; -t138 * t75 + (m(7) * (t116 * t44 + t120 * t43) - t143) * pkin(5) + t136 + t171; -t138 * t103 + (m(7) * (t116 * t55 + t120 * t54) - t143) * pkin(5) + t135 + t171; m(7) * (t116 * t91 + t120 * t89) * pkin(5) - t96 - t63; Ifges(6,3) + Ifges(7,3) + m(7) * (t116 ^ 2 + t120 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t130; t132; t136; t135; -t63; Ifges(7,3) + t130; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
