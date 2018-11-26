% Calculate joint inertia matrix for
% S6RRRPRR6
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:54:04
% EndTime: 2018-11-23 17:54:05
% DurationCPUTime: 1.35s
% Computational Cost: add. (3714->338), mult. (7520->489), div. (0->0), fcn. (8303->10), ass. (0->135)
t172 = 2 * pkin(7);
t137 = sin(qJ(5));
t141 = cos(qJ(5));
t134 = sin(pkin(11));
t135 = cos(pkin(11));
t138 = sin(qJ(3));
t142 = cos(qJ(3));
t104 = t134 * t142 + t135 * t138;
t139 = sin(qJ(2));
t94 = t104 * t139;
t103 = -t134 * t138 + t135 * t142;
t95 = t103 * t139;
t61 = -t137 * t95 - t141 * t94;
t62 = -t137 * t94 + t141 * t95;
t171 = -Ifges(6,5) * t62 - Ifges(6,6) * t61;
t157 = t139 * t142;
t170 = -Ifges(4,5) * t157 - Ifges(5,5) * t95 + Ifges(5,6) * t94;
t169 = pkin(3) * t134;
t143 = cos(qJ(2));
t168 = pkin(7) * t143;
t129 = t139 * pkin(7);
t136 = sin(qJ(6));
t140 = cos(qJ(6));
t123 = pkin(3) * t135 + pkin(4);
t97 = t141 * t123 - t137 * t169;
t96 = pkin(5) + t97;
t98 = t123 * t137 + t141 * t169;
t64 = t136 * t96 + t140 * t98;
t167 = t64 * mrSges(7,2);
t166 = t97 * mrSges(6,1);
t165 = t98 * mrSges(6,2);
t164 = -Ifges(7,3) - Ifges(6,3);
t163 = -qJ(4) - pkin(8);
t27 = -t136 * t62 + t140 * t61;
t28 = t136 * t61 + t140 * t62;
t162 = -Ifges(7,5) * t28 - Ifges(7,6) * t27;
t113 = -pkin(2) * t143 - pkin(8) * t139 - pkin(1);
t106 = t142 * t113;
t75 = -qJ(4) * t157 + t106 + (-pkin(7) * t138 - pkin(3)) * t143;
t158 = t138 * t139;
t87 = t138 * t113 + t142 * t168;
t81 = -qJ(4) * t158 + t87;
t47 = -t134 * t81 + t135 * t75;
t33 = -pkin(4) * t143 - pkin(9) * t95 + t47;
t48 = t134 * t75 + t135 * t81;
t37 = -pkin(9) * t94 + t48;
t14 = t137 * t33 + t141 * t37;
t114 = t163 * t138;
t116 = t163 * t142;
t82 = t135 * t114 + t116 * t134;
t65 = -pkin(9) * t104 + t82;
t83 = t134 * t114 - t135 * t116;
t66 = pkin(9) * t103 + t83;
t36 = t137 * t65 + t141 * t66;
t161 = Ifges(4,4) * t138;
t160 = Ifges(4,4) * t142;
t159 = t136 * mrSges(7,2);
t112 = pkin(3) * t158 + t129;
t156 = t138 ^ 2 + t142 ^ 2;
t155 = pkin(5) * t159;
t124 = -pkin(3) * t142 - pkin(2);
t67 = t94 * mrSges(5,1) + t95 * mrSges(5,2);
t31 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t73 = t103 * t141 - t104 * t137;
t74 = t103 * t137 + t104 * t141;
t44 = -t73 * mrSges(6,1) + t74 * mrSges(6,2);
t11 = -t27 * mrSges(7,1) + t28 * mrSges(7,2);
t42 = -t136 * t74 + t140 * t73;
t43 = t136 * t73 + t140 * t74;
t15 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t78 = -t103 * mrSges(5,1) + t104 * mrSges(5,2);
t13 = -t137 * t37 + t141 * t33;
t35 = -t137 * t66 + t141 * t65;
t154 = Ifges(5,3) + Ifges(4,3) - t164;
t63 = -t136 * t98 + t140 * t96;
t60 = t63 * mrSges(7,1);
t153 = Ifges(7,3) + t60 - t167;
t76 = pkin(4) * t94 + t112;
t7 = -pkin(5) * t143 - pkin(10) * t62 + t13;
t8 = pkin(10) * t61 + t14;
t2 = -t136 * t8 + t140 * t7;
t3 = t136 * t7 + t140 * t8;
t152 = t2 * mrSges(7,1) - t3 * mrSges(7,2) - t162;
t40 = Ifges(7,6) * t42;
t41 = Ifges(7,5) * t43;
t18 = -pkin(10) * t74 + t35;
t19 = pkin(10) * t73 + t36;
t5 = -t136 * t19 + t140 * t18;
t6 = t136 * t18 + t140 * t19;
t151 = t5 * mrSges(7,1) - t6 * mrSges(7,2) + t40 + t41;
t150 = mrSges(4,1) * t138 + mrSges(4,2) * t142;
t88 = -pkin(4) * t103 + t124;
t149 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t152 - t171;
t71 = Ifges(6,6) * t73;
t72 = Ifges(6,5) * t74;
t148 = t35 * mrSges(6,1) - t36 * mrSges(6,2) + t151 + t71 + t72;
t145 = pkin(7) ^ 2;
t133 = t143 ^ 2;
t131 = t139 ^ 2;
t128 = t131 * t145;
t127 = Ifges(4,5) * t138;
t126 = Ifges(4,6) * t142;
t125 = t140 * pkin(5) * mrSges(7,1);
t118 = Ifges(4,1) * t138 + t160;
t117 = Ifges(4,2) * t142 + t161;
t115 = -mrSges(4,1) * t142 + mrSges(4,2) * t138;
t111 = -mrSges(4,1) * t143 - mrSges(4,3) * t157;
t110 = mrSges(4,2) * t143 - mrSges(4,3) * t158;
t102 = t150 * t139;
t101 = Ifges(5,5) * t104;
t100 = Ifges(5,6) * t103;
t93 = -Ifges(4,5) * t143 + (Ifges(4,1) * t142 - t161) * t139;
t92 = -Ifges(4,6) * t143 + (-Ifges(4,2) * t138 + t160) * t139;
t86 = -t138 * t168 + t106;
t85 = -mrSges(5,1) * t143 - mrSges(5,3) * t95;
t84 = mrSges(5,2) * t143 - mrSges(5,3) * t94;
t80 = Ifges(5,1) * t104 + Ifges(5,4) * t103;
t79 = Ifges(5,4) * t104 + Ifges(5,2) * t103;
t56 = Ifges(5,1) * t95 - Ifges(5,4) * t94 - Ifges(5,5) * t143;
t55 = Ifges(5,4) * t95 - Ifges(5,2) * t94 - Ifges(5,6) * t143;
t51 = -mrSges(6,1) * t143 - mrSges(6,3) * t62;
t50 = mrSges(6,2) * t143 + mrSges(6,3) * t61;
t49 = -pkin(5) * t73 + t88;
t46 = Ifges(6,1) * t74 + Ifges(6,4) * t73;
t45 = Ifges(6,4) * t74 + Ifges(6,2) * t73;
t38 = -pkin(5) * t61 + t76;
t26 = Ifges(6,1) * t62 + Ifges(6,4) * t61 - Ifges(6,5) * t143;
t25 = Ifges(6,4) * t62 + Ifges(6,2) * t61 - Ifges(6,6) * t143;
t21 = -mrSges(7,1) * t143 - mrSges(7,3) * t28;
t20 = mrSges(7,2) * t143 + mrSges(7,3) * t27;
t17 = Ifges(7,1) * t43 + Ifges(7,4) * t42;
t16 = Ifges(7,4) * t43 + Ifges(7,2) * t42;
t10 = Ifges(7,1) * t28 + Ifges(7,4) * t27 - Ifges(7,5) * t143;
t9 = Ifges(7,4) * t28 + Ifges(7,2) * t27 - Ifges(7,6) * t143;
t1 = [0.2e1 * t87 * t110 + 0.2e1 * t86 * t111 + 0.2e1 * t112 * t67 - t94 * t55 + t95 * t56 + 0.2e1 * t48 * t84 + 0.2e1 * t47 * t85 + 0.2e1 * t76 * t31 + t61 * t25 + t62 * t26 + 0.2e1 * t14 * t50 + 0.2e1 * t13 * t51 + 0.2e1 * t38 * t11 + t27 * t9 + t28 * t10 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + m(7) * (t2 ^ 2 + t3 ^ 2 + t38 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2 + t76 ^ 2) + m(5) * (t112 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(4) * (t86 ^ 2 + t87 ^ 2 + t128) + m(3) * (pkin(1) ^ 2 + t133 * t145 + t128) + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t139 + t102 * t172 - t138 * t92 + t142 * t93) * t139 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t154) * t143 + (Ifges(4,6) * t138 + (2 * Ifges(3,4))) * t139 + t162 + t170 + t171) * t143 + Ifges(2,3) + (t131 + t133) * mrSges(3,3) * t172; (-t13 * t74 + t14 * t73) * mrSges(6,3) + (t103 * t48 - t104 * t47) * mrSges(5,3) + t124 * t67 + t103 * t55 / 0.2e1 + t104 * t56 / 0.2e1 + t112 * t78 - t94 * t79 / 0.2e1 + t95 * t80 / 0.2e1 - pkin(2) * t102 + t83 * t84 + t82 * t85 + t88 * t31 + t73 * t25 / 0.2e1 + t74 * t26 / 0.2e1 + t76 * t44 + t61 * t45 / 0.2e1 + t62 * t46 / 0.2e1 + t49 * t11 + t36 * t50 + t35 * t51 + t38 * t15 + t42 * t9 / 0.2e1 + t43 * t10 / 0.2e1 + t27 * t16 / 0.2e1 + t28 * t17 / 0.2e1 + t6 * t20 + t5 * t21 + (t92 / 0.2e1 + t87 * mrSges(4,3) + pkin(8) * t110) * t142 + (t93 / 0.2e1 - t86 * mrSges(4,3) - pkin(8) * t111) * t138 + m(7) * (t2 * t5 + t3 * t6 + t38 * t49) + m(6) * (t13 * t35 + t14 * t36 + t76 * t88) + m(5) * (t112 * t124 + t47 * t82 + t48 * t83) + m(4) * (-pkin(2) * t129 + (-t138 * t86 + t142 * t87) * pkin(8)) + (-t127 / 0.2e1 - t126 / 0.2e1 - t101 / 0.2e1 - t100 / 0.2e1 - t72 / 0.2e1 - t71 / 0.2e1 - t41 / 0.2e1 - t40 / 0.2e1 + Ifges(3,6) - pkin(7) * mrSges(3,2)) * t143 + (-t2 * t43 + t3 * t42) * mrSges(7,3) + (Ifges(3,5) - t138 * t117 / 0.2e1 + t142 * t118 / 0.2e1 + (-mrSges(3,1) + t115) * pkin(7)) * t139; -0.2e1 * pkin(2) * t115 + t103 * t79 + t104 * t80 + t142 * t117 + t138 * t118 + 0.2e1 * t124 * t78 + 0.2e1 * t49 * t15 + t42 * t16 + t43 * t17 + 0.2e1 * t88 * t44 + t73 * t45 + t74 * t46 + Ifges(3,3) + m(7) * (t49 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2 + t88 ^ 2) + m(5) * (t124 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(4) * (pkin(8) ^ 2 * t156 + pkin(2) ^ 2) + 0.2e1 * (t42 * t6 - t43 * t5) * mrSges(7,3) + 0.2e1 * (-t35 * t74 + t36 * t73) * mrSges(6,3) + 0.2e1 * (t103 * t83 - t104 * t82) * mrSges(5,3) + 0.2e1 * t156 * pkin(8) * mrSges(4,3); t97 * t51 + t98 * t50 + t86 * mrSges(4,1) - t87 * mrSges(4,2) + t63 * t21 + (t134 * t84 + t135 * t85 + m(5) * (t134 * t48 + t135 * t47)) * pkin(3) + t64 * t20 - t48 * mrSges(5,2) + t47 * mrSges(5,1) + t149 + m(7) * (t2 * t63 + t3 * t64) + m(6) * (t13 * t97 + t14 * t98) - t154 * t143 - Ifges(4,6) * t158 - t170; t82 * mrSges(5,1) - t83 * mrSges(5,2) + (m(5) * (t134 * t83 + t135 * t82) + (t103 * t134 - t104 * t135) * mrSges(5,3)) * pkin(3) + t148 + t126 + t127 + t100 + t101 + (t73 * t98 - t74 * t97) * mrSges(6,3) + (t42 * t64 - t43 * t63) * mrSges(7,3) - t150 * pkin(8) + m(7) * (t5 * t63 + t6 * t64) + m(6) * (t35 * t97 + t36 * t98); 0.2e1 * t166 - 0.2e1 * t165 - 0.2e1 * t167 + 0.2e1 * t60 + m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t97 ^ 2 + t98 ^ 2) + t154 + (0.2e1 * mrSges(5,1) * t135 - 0.2e1 * mrSges(5,2) * t134 + m(5) * (t134 ^ 2 + t135 ^ 2) * pkin(3)) * pkin(3); m(5) * t112 + m(6) * t76 + m(7) * t38 + t11 + t31 + t67; m(5) * t124 + m(6) * t88 + m(7) * t49 + t15 + t44 + t78; 0; m(5) + m(6) + m(7); t164 * t143 + (m(7) * (t136 * t3 + t140 * t2) + t140 * t21 + t136 * t20) * pkin(5) + t149; (m(7) * (t136 * t6 + t140 * t5) + (t136 * t42 - t140 * t43) * mrSges(7,3)) * pkin(5) + t148; t166 - t165 + Ifges(6,3) + t125 + (m(7) * (t136 * t64 + t140 * t63) - t159) * pkin(5) + t153; 0; -0.2e1 * t155 + 0.2e1 * t125 + m(7) * (t136 ^ 2 + t140 ^ 2) * pkin(5) ^ 2 - t164; -Ifges(7,3) * t143 + t152; t151; t153; 0; Ifges(7,3) + t125 - t155; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
