% Calculate joint inertia matrix for
% S6RRRRPR6
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:15:11
% EndTime: 2018-11-23 18:15:13
% DurationCPUTime: 1.40s
% Computational Cost: add. (3888->354), mult. (7847->513), div. (0->0), fcn. (8618->10), ass. (0->139)
t176 = 2 * pkin(7);
t139 = sin(qJ(6));
t143 = cos(qJ(6));
t137 = sin(pkin(11));
t138 = cos(pkin(11));
t146 = cos(qJ(2));
t140 = sin(qJ(4));
t144 = cos(qJ(4));
t142 = sin(qJ(2));
t116 = -pkin(2) * t146 - pkin(8) * t142 - pkin(1);
t145 = cos(qJ(3));
t107 = t145 * t116;
t141 = sin(qJ(3));
t161 = t142 * t145;
t76 = -pkin(9) * t161 + t107 + (-pkin(7) * t141 - pkin(3)) * t146;
t162 = t141 * t142;
t170 = pkin(7) * t146;
t89 = t141 * t116 + t145 * t170;
t82 = -pkin(9) * t162 + t89;
t48 = -t140 * t82 + t144 * t76;
t108 = -t140 * t141 + t144 * t145;
t98 = t108 * t142;
t32 = -pkin(4) * t146 - qJ(5) * t98 + t48;
t49 = t140 * t76 + t144 * t82;
t109 = t140 * t145 + t141 * t144;
t97 = t109 * t142;
t37 = -qJ(5) * t97 + t49;
t13 = -t137 * t37 + t138 * t32;
t60 = -t137 * t97 + t138 * t98;
t7 = -pkin(5) * t146 - pkin(10) * t60 + t13;
t14 = t137 * t32 + t138 * t37;
t59 = -t137 * t98 - t138 * t97;
t8 = pkin(10) * t59 + t14;
t2 = -t139 * t8 + t143 * t7;
t3 = t139 * t7 + t143 * t8;
t175 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t174 = -pkin(9) - pkin(8);
t172 = pkin(3) * t140;
t171 = pkin(4) * t137;
t132 = t142 * pkin(7);
t127 = pkin(3) * t144 + pkin(4);
t101 = t127 * t137 + t138 * t172;
t99 = t138 * t127 - t137 * t172;
t96 = pkin(5) + t99;
t65 = t101 * t143 + t139 * t96;
t169 = t65 * mrSges(7,2);
t168 = t99 * mrSges(6,1);
t27 = -t139 * t60 + t143 * t59;
t28 = t139 * t59 + t143 * t60;
t167 = -Ifges(7,5) * t28 - Ifges(7,6) * t27;
t121 = t174 * t141;
t122 = t174 * t145;
t84 = t144 * t121 + t122 * t140;
t66 = -qJ(5) * t109 + t84;
t85 = t140 * t121 - t144 * t122;
t67 = qJ(5) * t108 + t85;
t36 = t137 * t66 + t138 * t67;
t166 = Ifges(4,4) * t141;
t165 = Ifges(4,4) * t145;
t164 = t101 * mrSges(6,2);
t126 = pkin(4) * t138 + pkin(5);
t102 = t126 * t139 + t143 * t171;
t163 = t102 * mrSges(7,2);
t115 = pkin(3) * t162 + t132;
t160 = t141 ^ 2 + t145 ^ 2;
t159 = -Ifges(7,3) - Ifges(5,3) - Ifges(6,3);
t128 = -pkin(3) * t145 - pkin(2);
t31 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t74 = t108 * t138 - t109 * t137;
t75 = t108 * t137 + t109 * t138;
t44 = -t74 * mrSges(6,1) + t75 * mrSges(6,2);
t11 = -t27 * mrSges(7,1) + t28 * mrSges(7,2);
t42 = -t139 * t75 + t143 * t74;
t43 = t139 * t74 + t143 * t75;
t15 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t35 = -t137 * t67 + t138 * t66;
t158 = Ifges(4,3) - t159;
t77 = pkin(4) * t97 + t115;
t40 = Ifges(7,6) * t42;
t41 = Ifges(7,5) * t43;
t18 = -pkin(10) * t75 + t35;
t19 = pkin(10) * t74 + t36;
t5 = -t139 * t19 + t143 * t18;
t6 = t139 * t18 + t143 * t19;
t157 = t5 * mrSges(7,1) - t6 * mrSges(7,2) + t40 + t41;
t156 = mrSges(4,1) * t141 + mrSges(4,2) * t145;
t155 = t138 * mrSges(6,1) - t137 * mrSges(6,2);
t154 = -Ifges(5,5) * t98 - Ifges(6,5) * t60 + Ifges(5,6) * t97 - Ifges(6,6) * t59 + t167;
t90 = -pkin(4) * t108 + t128;
t153 = (mrSges(5,1) * t144 - mrSges(5,2) * t140) * pkin(3);
t152 = t48 * mrSges(5,1) + t13 * mrSges(6,1) - t49 * mrSges(5,2) - t14 * mrSges(6,2) - t154 + t175;
t104 = Ifges(5,6) * t108;
t105 = Ifges(5,5) * t109;
t72 = Ifges(6,6) * t74;
t73 = Ifges(6,5) * t75;
t151 = t84 * mrSges(5,1) + t35 * mrSges(6,1) - t85 * mrSges(5,2) - t36 * mrSges(6,2) + t104 + t105 + t157 + t72 + t73;
t148 = pkin(7) ^ 2;
t136 = t146 ^ 2;
t134 = t142 ^ 2;
t131 = t134 * t148;
t130 = Ifges(4,5) * t141;
t129 = Ifges(4,6) * t145;
t123 = Ifges(4,5) * t161;
t119 = Ifges(4,1) * t141 + t165;
t118 = Ifges(4,2) * t145 + t166;
t117 = -mrSges(4,1) * t145 + mrSges(4,2) * t141;
t114 = -mrSges(4,1) * t146 - mrSges(4,3) * t161;
t113 = mrSges(4,2) * t146 - mrSges(4,3) * t162;
t103 = t156 * t142;
t100 = t126 * t143 - t139 * t171;
t95 = t100 * mrSges(7,1);
t94 = -Ifges(4,5) * t146 + (Ifges(4,1) * t145 - t166) * t142;
t93 = -Ifges(4,6) * t146 + (-Ifges(4,2) * t141 + t165) * t142;
t88 = -t141 * t170 + t107;
t87 = -mrSges(5,1) * t146 - mrSges(5,3) * t98;
t86 = mrSges(5,2) * t146 - mrSges(5,3) * t97;
t81 = Ifges(5,1) * t109 + Ifges(5,4) * t108;
t80 = Ifges(5,4) * t109 + Ifges(5,2) * t108;
t79 = -mrSges(5,1) * t108 + mrSges(5,2) * t109;
t68 = mrSges(5,1) * t97 + mrSges(5,2) * t98;
t64 = -t101 * t139 + t143 * t96;
t58 = Ifges(5,1) * t98 - Ifges(5,4) * t97 - Ifges(5,5) * t146;
t57 = Ifges(5,4) * t98 - Ifges(5,2) * t97 - Ifges(5,6) * t146;
t56 = t64 * mrSges(7,1);
t52 = -mrSges(6,1) * t146 - mrSges(6,3) * t60;
t51 = mrSges(6,2) * t146 + mrSges(6,3) * t59;
t50 = -pkin(5) * t74 + t90;
t46 = Ifges(6,1) * t75 + Ifges(6,4) * t74;
t45 = Ifges(6,4) * t75 + Ifges(6,2) * t74;
t38 = -pkin(5) * t59 + t77;
t26 = Ifges(6,1) * t60 + Ifges(6,4) * t59 - Ifges(6,5) * t146;
t25 = Ifges(6,4) * t60 + Ifges(6,2) * t59 - Ifges(6,6) * t146;
t21 = -mrSges(7,1) * t146 - mrSges(7,3) * t28;
t20 = mrSges(7,2) * t146 + mrSges(7,3) * t27;
t17 = Ifges(7,1) * t43 + Ifges(7,4) * t42;
t16 = Ifges(7,4) * t43 + Ifges(7,2) * t42;
t10 = Ifges(7,1) * t28 + Ifges(7,4) * t27 - Ifges(7,5) * t146;
t9 = Ifges(7,4) * t28 + Ifges(7,2) * t27 - Ifges(7,6) * t146;
t1 = [0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + (0.2e1 * pkin(1) * mrSges(3,1) - t123 + (Ifges(3,2) + t158) * t146 + (Ifges(4,6) * t141 + (2 * Ifges(3,4))) * t142 + t154) * t146 + m(3) * (pkin(1) ^ 2 + t136 * t148 + t131) + m(4) * (t88 ^ 2 + t89 ^ 2 + t131) + m(5) * (t115 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2 + t77 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t38 ^ 2) + t27 * t9 + t28 * t10 + 0.2e1 * t38 * t11 + 0.2e1 * t14 * t51 + 0.2e1 * t13 * t52 + t59 * t25 + t60 * t26 + Ifges(2,3) + 0.2e1 * t77 * t31 + 0.2e1 * t49 * t86 + 0.2e1 * t48 * t87 - t97 * t57 + t98 * t58 + 0.2e1 * t89 * t113 + 0.2e1 * t88 * t114 + 0.2e1 * t115 * t68 + (t134 + t136) * mrSges(3,3) * t176 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t142 + t103 * t176 - t141 * t93 + t145 * t94) * t142; m(5) * (t115 * t128 + t48 * t84 + t49 * t85) + m(6) * (t13 * t35 + t14 * t36 + t77 * t90) + m(7) * (t2 * t5 + t3 * t6 + t38 * t50) + t6 * t20 + t5 * t21 + (pkin(8) * t113 + t89 * mrSges(4,3) + t93 / 0.2e1) * t145 + (-pkin(8) * t114 - t88 * mrSges(4,3) + t94 / 0.2e1) * t141 + (-pkin(7) * mrSges(3,2) - t41 / 0.2e1 - t40 / 0.2e1 - t73 / 0.2e1 - t72 / 0.2e1 - t105 / 0.2e1 - t104 / 0.2e1 - t130 / 0.2e1 - t129 / 0.2e1 + Ifges(3,6)) * t146 + (-t2 * t43 + t3 * t42) * mrSges(7,3) + (-t13 * t75 + t14 * t74) * mrSges(6,3) + (t108 * t49 - t109 * t48) * mrSges(5,3) + m(4) * (-pkin(2) * t132 + (-t141 * t88 + t145 * t89) * pkin(8)) + t27 * t16 / 0.2e1 + t28 * t17 / 0.2e1 + t38 * t15 + t42 * t9 / 0.2e1 + t43 * t10 / 0.2e1 + t50 * t11 + t36 * t51 + t35 * t52 + t59 * t45 / 0.2e1 + t60 * t46 / 0.2e1 + t74 * t25 / 0.2e1 + t75 * t26 / 0.2e1 + t77 * t44 + t85 * t86 + t84 * t87 + t90 * t31 - t97 * t80 / 0.2e1 + t98 * t81 / 0.2e1 - pkin(2) * t103 + t108 * t57 / 0.2e1 + t109 * t58 / 0.2e1 + t115 * t79 + t128 * t68 + (t145 * t119 / 0.2e1 - t141 * t118 / 0.2e1 + Ifges(3,5) + (t117 - mrSges(3,1)) * pkin(7)) * t142; -0.2e1 * pkin(2) * t117 + t108 * t80 + t109 * t81 + t145 * t118 + t141 * t119 + 0.2e1 * t128 * t79 + 0.2e1 * t50 * t15 + t42 * t16 + t43 * t17 + 0.2e1 * t90 * t44 + t74 * t45 + t75 * t46 + Ifges(3,3) + m(4) * (t160 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t128 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2 + t90 ^ 2) + m(7) * (t5 ^ 2 + t50 ^ 2 + t6 ^ 2) + 0.2e1 * (t42 * t6 - t43 * t5) * mrSges(7,3) + 0.2e1 * (-t35 * t75 + t36 * t74) * mrSges(6,3) + 0.2e1 * (t108 * t85 - t109 * t84) * mrSges(5,3) + 0.2e1 * t160 * pkin(8) * mrSges(4,3); -Ifges(4,6) * t162 + t152 + m(6) * (t101 * t14 + t13 * t99) + m(7) * (t2 * t64 + t3 * t65) - t158 * t146 + t123 + (m(5) * (t140 * t49 + t144 * t48) + t144 * t87 + t140 * t86) * pkin(3) + t64 * t21 + t65 * t20 + t88 * mrSges(4,1) - t89 * mrSges(4,2) + t99 * t52 + t101 * t51; t151 + m(6) * (t101 * t36 + t35 * t99) + m(7) * (t5 * t64 + t6 * t65) + t129 + t130 - t156 * pkin(8) + (t42 * t65 - t43 * t64) * mrSges(7,3) + (t101 * t74 - t75 * t99) * mrSges(6,3) + (m(5) * (t140 * t85 + t144 * t84) + (t108 * t140 - t109 * t144) * mrSges(5,3)) * pkin(3); 0.2e1 * t168 - 0.2e1 * t164 - 0.2e1 * t169 + 0.2e1 * t56 + 0.2e1 * t153 + m(6) * (t101 ^ 2 + t99 ^ 2) + m(7) * (t64 ^ 2 + t65 ^ 2) + m(5) * (t140 ^ 2 + t144 ^ 2) * pkin(3) ^ 2 + t158; t152 + m(7) * (t100 * t2 + t102 * t3) + (t137 * t51 + t138 * t52 + m(6) * (t13 * t138 + t137 * t14)) * pkin(4) + t159 * t146 + t100 * t21 + t102 * t20; t151 + (m(6) * (t137 * t36 + t138 * t35) + (t137 * t74 - t138 * t75) * mrSges(6,3)) * pkin(4) + m(7) * (t100 * t5 + t102 * t6) + (-t100 * t43 + t102 * t42) * mrSges(7,3); t95 + t56 + m(7) * (t100 * t64 + t102 * t65) - t164 + t168 + t153 + (-t65 - t102) * mrSges(7,2) + (m(6) * (t101 * t137 + t138 * t99) + t155) * pkin(4) - t159; -0.2e1 * t163 + 0.2e1 * t95 + m(7) * (t100 ^ 2 + t102 ^ 2) - t159 + (0.2e1 * t155 + m(6) * (t137 ^ 2 + t138 ^ 2) * pkin(4)) * pkin(4); m(6) * t77 + m(7) * t38 + t11 + t31; m(6) * t90 + m(7) * t50 + t15 + t44; 0; 0; m(6) + m(7); -Ifges(7,3) * t146 - t167 + t175; t157; Ifges(7,3) + t56 - t169; Ifges(7,3) + t95 - t163; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
