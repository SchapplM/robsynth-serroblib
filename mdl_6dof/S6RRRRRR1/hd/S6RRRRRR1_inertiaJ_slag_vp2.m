% Calculate joint inertia matrix for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:54
% EndTime: 2019-03-10 03:27:56
% DurationCPUTime: 1.19s
% Computational Cost: add. (4061->276), mult. (7543->398), div. (0->0), fcn. (8750->10), ass. (0->123)
t103 = cos(qJ(6));
t98 = sin(qJ(6));
t144 = t98 * mrSges(7,3);
t104 = cos(qJ(5));
t100 = sin(qJ(4));
t105 = cos(qJ(4));
t101 = sin(qJ(3));
t102 = sin(qJ(2));
t106 = cos(qJ(3));
t107 = cos(qJ(2));
t72 = -t101 * t102 + t106 * t107;
t73 = t101 * t107 + t102 * t106;
t46 = -t100 * t73 + t105 * t72;
t47 = t100 * t72 + t105 * t73;
t99 = sin(qJ(5));
t29 = -t104 * t46 + t47 * t99;
t30 = t104 * t47 + t46 * t99;
t15 = -mrSges(7,2) * t29 - t144 * t30;
t135 = t103 * t30;
t16 = mrSges(7,1) * t29 - mrSges(7,3) * t135;
t121 = t103 * t15 - t98 * t16;
t76 = -mrSges(7,1) * t103 + mrSges(7,2) * t98;
t156 = pkin(5) * t76;
t94 = t98 ^ 2;
t155 = mrSges(7,3) * t94;
t88 = pkin(11) * t155;
t96 = t103 ^ 2;
t154 = mrSges(7,3) * t96;
t89 = pkin(11) * t154;
t166 = t88 + t89 - t156;
t122 = mrSges(7,1) * t98 + mrSges(7,2) * t103;
t14 = t122 * t30;
t133 = t101 * t105;
t134 = t100 * t101;
t157 = -pkin(8) - pkin(7);
t20 = -t47 * pkin(9) + ((t100 * t106 + t133) * t107 + (t105 * t106 - t134) * t102) * t157;
t113 = -t47 * pkin(10) + t20;
t128 = t157 * t102;
t129 = t157 * t107;
t50 = t101 * t129 + t106 * t128;
t51 = t101 * t128 - t106 * t129;
t21 = t105 * (t72 * pkin(9) + t51) + t100 * (-t73 * pkin(9) + t50);
t18 = pkin(10) * t46 + t21;
t6 = -t104 * t113 + t99 * t18;
t165 = m(7) * t6 + t14;
t86 = pkin(2) * t106 + pkin(3);
t64 = -pkin(2) * t134 + t105 * t86;
t61 = pkin(4) + t64;
t66 = pkin(2) * t133 + t100 * t86;
t41 = t104 * t61 - t66 * t99;
t39 = -pkin(5) - t41;
t32 = t39 * t76;
t42 = t104 * t66 + t99 * t61;
t40 = pkin(11) + t42;
t33 = t40 * t155;
t34 = t40 * t154;
t35 = t41 * mrSges(6,1);
t164 = t32 + t33 + t34 + t35;
t150 = pkin(4) * t104;
t84 = -pkin(5) - t150;
t67 = t84 * t76;
t83 = pkin(4) * t99 + pkin(11);
t74 = t83 * t155;
t75 = t83 * t154;
t90 = mrSges(6,1) * t150;
t163 = t67 + t74 + t75 + t90;
t162 = Ifges(5,3) + t163;
t87 = -pkin(2) * t107 - pkin(1);
t55 = -pkin(3) * t72 + t87;
t31 = -pkin(4) * t46 + t55;
t13 = pkin(5) * t29 - pkin(11) * t30 + t31;
t8 = t104 * t18 + t113 * t99;
t3 = t103 * t8 + t13 * t98;
t149 = t103 * t3;
t2 = t103 * t13 - t8 * t98;
t123 = -t2 * t98 + t149;
t161 = m(7) * t123 + t121;
t160 = t6 ^ 2;
t159 = 0.2e1 * t6;
t158 = 0.2e1 * t31;
t153 = Ifges(7,4) * t98;
t152 = pkin(3) * t100;
t151 = pkin(3) * t105;
t148 = t30 * t98;
t147 = t42 * mrSges(6,2);
t85 = pkin(4) + t151;
t65 = t104 * t152 + t99 * t85;
t146 = t65 * mrSges(6,2);
t145 = t66 * mrSges(5,2);
t142 = t99 * mrSges(6,2);
t141 = Ifges(7,5) * t135 + Ifges(7,3) * t29;
t140 = Ifges(7,5) * t98 + Ifges(7,6) * t103;
t139 = t94 + t96;
t138 = t102 ^ 2 + t107 ^ 2;
t137 = Ifges(7,4) * t103;
t132 = pkin(4) * t142;
t131 = mrSges(5,2) * t152;
t77 = Ifges(7,2) * t103 + t153;
t78 = Ifges(7,1) * t98 + t137;
t130 = t103 * t77 + t98 * t78 + Ifges(6,3);
t62 = pkin(11) + t65;
t127 = t139 * t62;
t126 = t139 * t83;
t125 = Ifges(5,3) + t130;
t124 = Ifges(4,3) + t125;
t63 = t104 * t85 - t152 * t99;
t60 = -pkin(5) - t63;
t48 = t60 * t76;
t52 = t62 * t155;
t53 = t62 * t154;
t57 = t63 * mrSges(6,1);
t120 = t48 + t52 + t53 + t57 + t130;
t119 = (mrSges(4,1) * t106 - mrSges(4,2) * t101) * pkin(2);
t118 = t130 - t147 + t164;
t117 = t120 - t146;
t11 = Ifges(7,6) * t29 + (-Ifges(7,2) * t98 + t137) * t30;
t12 = Ifges(7,5) * t29 + (Ifges(7,1) * t103 - t153) * t30;
t116 = -t8 * mrSges(6,2) + mrSges(7,3) * t149 - t144 * t2 + t103 * t11 / 0.2e1 - t77 * t148 / 0.2e1 + t78 * t135 / 0.2e1 + Ifges(6,5) * t30 + t98 * t12 / 0.2e1 + (-mrSges(6,1) + t76) * t6 + (t140 / 0.2e1 - Ifges(6,6)) * t29;
t115 = t20 * mrSges(5,1) - t21 * mrSges(5,2) + Ifges(5,5) * t47 + Ifges(5,6) * t46 + t116;
t114 = t50 * mrSges(4,1) - t51 * mrSges(4,2) + Ifges(4,5) * t73 + Ifges(4,6) * t72 + t115;
t91 = mrSges(5,1) * t151;
t58 = t64 * mrSges(5,1);
t1 = [-0.2e1 * pkin(1) * (-t107 * mrSges(3,1) + t102 * mrSges(3,2)) + t102 * (Ifges(3,1) * t102 + Ifges(3,4) * t107) + t107 * (Ifges(3,4) * t102 + Ifges(3,2) * t107) + 0.2e1 * t87 * (-t72 * mrSges(4,1) + t73 * mrSges(4,2)) + t73 * (Ifges(4,1) * t73 + Ifges(4,4) * t72) + t72 * (Ifges(4,4) * t73 + Ifges(4,2) * t72) + t47 * (Ifges(5,1) * t47 + Ifges(5,4) * t46) + t46 * (Ifges(5,4) * t47 + Ifges(5,2) * t46) + 0.2e1 * t55 * (-t46 * mrSges(5,1) + t47 * mrSges(5,2)) + t14 * t159 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + (mrSges(6,1) * t158 - 0.2e1 * mrSges(6,3) * t8 + Ifges(6,2) * t29 + t141) * t29 + (mrSges(6,2) * t158 + mrSges(6,3) * t159 + Ifges(6,1) * t30 + t103 * t12 - t98 * t11 + (-Ifges(7,6) * t98 - (2 * Ifges(6,4))) * t29) * t30 + m(7) * (t2 ^ 2 + t3 ^ 2 + t160) + m(6) * (t31 ^ 2 + t8 ^ 2 + t160) + m(5) * (t20 ^ 2 + t21 ^ 2 + t55 ^ 2) + m(4) * (t50 ^ 2 + t51 ^ 2 + t87 ^ 2) + m(3) * (pkin(7) ^ 2 * t138 + pkin(1) ^ 2) + 0.2e1 * (-t20 * t47 + t21 * t46) * mrSges(5,3) + 0.2e1 * (-t50 * t73 + t51 * t72) * mrSges(4,3) + 0.2e1 * t138 * pkin(7) * mrSges(3,3); (-mrSges(3,1) * t102 - mrSges(3,2) * t107) * pkin(7) + (-t29 * t42 - t30 * t41) * mrSges(6,3) + (t46 * t66 - t47 * t64) * mrSges(5,3) + Ifges(3,6) * t107 + Ifges(3,5) * t102 + t39 * t14 + t121 * t40 + m(6) * (-t41 * t6 + t42 * t8) + m(5) * (t20 * t64 + t21 * t66) + t114 + m(7) * (t123 * t40 + t39 * t6) + ((t101 * t72 - t106 * t73) * mrSges(4,3) + m(4) * (t101 * t51 + t106 * t50)) * pkin(2); -0.2e1 * t145 - 0.2e1 * t147 + Ifges(3,3) + 0.2e1 * t32 + 0.2e1 * t33 + 0.2e1 * t34 + 0.2e1 * t35 + 0.2e1 * t58 + 0.2e1 * t119 + m(7) * (t139 * t40 ^ 2 + t39 ^ 2) + m(6) * (t41 ^ 2 + t42 ^ 2) + m(5) * (t64 ^ 2 + t66 ^ 2) + m(4) * (t101 ^ 2 + t106 ^ 2) * pkin(2) ^ 2 + t124; t60 * t14 + t121 * t62 + m(6) * (-t6 * t63 + t65 * t8) + (-t29 * t65 - t30 * t63) * mrSges(6,3) + t114 + m(7) * (t123 * t62 + t60 * t6) + (m(5) * (t100 * t21 + t105 * t20) + (t100 * t46 - t105 * t47) * mrSges(5,3)) * pkin(3); m(7) * (t127 * t40 + t60 * t39) + t120 + t119 + t91 + m(6) * (t41 * t63 + t42 * t65) + (-t42 - t65) * mrSges(6,2) + (-t66 - t152) * mrSges(5,2) + m(5) * (t100 * t66 + t105 * t64) * pkin(3) + t58 + Ifges(4,3) + Ifges(5,3) + t164; -0.2e1 * t131 - 0.2e1 * t146 + 0.2e1 * t48 + 0.2e1 * t52 + 0.2e1 * t53 + 0.2e1 * t57 + 0.2e1 * t91 + m(7) * (t139 * t62 ^ 2 + t60 ^ 2) + m(6) * (t63 ^ 2 + t65 ^ 2) + m(5) * (t100 ^ 2 + t105 ^ 2) * pkin(3) ^ 2 + t124; t115 + (m(6) * (-t104 * t6 + t8 * t99) + (-t104 * t30 - t29 * t99) * mrSges(6,3)) * pkin(4) + t165 * t84 + t161 * t83; -t145 + t118 + t58 + m(7) * (t126 * t40 + t84 * t39) + (m(6) * (t104 * t41 + t42 * t99) - t142) * pkin(4) + t162; t91 + t117 - t131 + m(7) * (t126 * t62 + t84 * t60) + (-t142 + m(6) * (t104 * t63 + t65 * t99)) * pkin(4) + t162; -0.2e1 * t132 + 0.2e1 * t67 + 0.2e1 * t74 + 0.2e1 * t75 + 0.2e1 * t90 + m(7) * (t139 * t83 ^ 2 + t84 ^ 2) + m(6) * (t104 ^ 2 + t99 ^ 2) * pkin(4) ^ 2 + t125; -t165 * pkin(5) + t161 * pkin(11) + t116; m(7) * (pkin(11) * t139 * t40 - pkin(5) * t39) + t118 + t166; m(7) * (-pkin(5) * t60 + pkin(11) * t127) + t117 + t166; m(7) * (-pkin(5) * t84 + pkin(11) * t126) - t132 + t130 + t163 + t166; -0.2e1 * t156 + m(7) * (pkin(11) ^ 2 * t139 + pkin(5) ^ 2) + 0.2e1 * t88 + 0.2e1 * t89 + t130; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t148 + t141; -t122 * t40 + t140; -t122 * t62 + t140; -t122 * t83 + t140; -pkin(11) * t122 + t140; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
