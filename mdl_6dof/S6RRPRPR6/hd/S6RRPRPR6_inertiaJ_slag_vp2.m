% Calculate joint inertia matrix for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:13
% EndTime: 2019-03-09 10:38:18
% DurationCPUTime: 1.53s
% Computational Cost: add. (2290->330), mult. (5406->453), div. (0->0), fcn. (5733->10), ass. (0->130)
t116 = sin(qJ(6));
t119 = cos(qJ(6));
t136 = t116 ^ 2 + t119 ^ 2;
t134 = m(7) * t136;
t171 = m(6) + t134;
t170 = Ifges(6,1) + Ifges(5,3);
t169 = Ifges(3,3) + Ifges(4,3);
t117 = sin(qJ(4));
t108 = t117 ^ 2;
t120 = cos(qJ(4));
t110 = t120 ^ 2;
t168 = t108 + t110;
t167 = 0.2e1 * t168;
t166 = m(6) * pkin(4);
t158 = -t116 / 0.2e1;
t112 = sin(pkin(11));
t114 = cos(pkin(11));
t113 = sin(pkin(6));
t121 = cos(qJ(2));
t141 = t113 * t121;
t118 = sin(qJ(2));
t142 = t113 * t118;
t61 = t112 * t142 - t114 * t141;
t62 = (t112 * t121 + t114 * t118) * t113;
t165 = Ifges(4,5) * t62 - Ifges(4,6) * t61;
t164 = mrSges(7,3) * t136 - mrSges(6,2);
t73 = (-pkin(2) * t121 - pkin(1)) * t113;
t163 = 0.2e1 * t73;
t115 = cos(pkin(6));
t47 = -t115 * t120 + t117 * t62;
t32 = -t116 * t61 + t119 * t47;
t162 = t32 / 0.2e1;
t144 = Ifges(7,4) * t119;
t66 = Ifges(7,5) * t117 + (-Ifges(7,1) * t116 - t144) * t120;
t161 = t66 / 0.2e1;
t105 = Ifges(7,5) * t119;
t160 = Ifges(7,6) * t158 + t105 / 0.2e1;
t159 = pkin(4) + pkin(10);
t157 = -t119 / 0.2e1;
t98 = pkin(2) * t112 + pkin(9);
t156 = pkin(5) + t98;
t155 = pkin(1) * t115;
t154 = pkin(4) * t120;
t153 = t61 * Ifges(6,4);
t152 = t61 * Ifges(6,5);
t94 = t121 * t155;
t67 = -pkin(8) * t142 + t94;
t151 = t67 * mrSges(3,1);
t68 = pkin(8) * t141 + t118 * t155;
t150 = t68 * mrSges(3,2);
t49 = pkin(2) * t115 + t94 + (-pkin(8) - qJ(3)) * t142;
t54 = qJ(3) * t141 + t68;
t30 = t112 * t49 + t114 * t54;
t26 = pkin(9) * t115 + t30;
t34 = pkin(3) * t61 - pkin(9) * t62 + t73;
t13 = t117 * t34 + t120 * t26;
t35 = mrSges(6,1) * t47 - mrSges(6,3) * t61;
t37 = -mrSges(5,2) * t61 - mrSges(5,3) * t47;
t149 = -t35 + t37;
t48 = t115 * t117 + t120 * t62;
t36 = t48 * mrSges(6,1) + t61 * mrSges(6,2);
t38 = mrSges(5,1) * t61 - mrSges(5,3) * t48;
t148 = t36 - t38;
t78 = mrSges(7,1) * t116 + mrSges(7,2) * t119;
t147 = t78 + mrSges(6,3);
t146 = t168 * t98 ^ 2;
t145 = Ifges(7,4) * t116;
t139 = t119 * t120;
t75 = -mrSges(7,2) * t117 - mrSges(7,3) * t139;
t143 = t116 * t75;
t140 = t116 * t120;
t138 = t119 * t159;
t137 = Ifges(5,5) * t117 + Ifges(5,6) * t120;
t33 = t116 * t47 + t119 * t61;
t6 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t48;
t99 = -pkin(2) * t114 - pkin(3);
t29 = -t112 * t54 + t114 * t49;
t12 = -t117 * t26 + t120 * t34;
t133 = t136 * t159;
t102 = t117 * qJ(5);
t132 = -t102 + t99;
t9 = -qJ(5) * t61 - t13;
t3 = pkin(5) * t48 - t159 * t61 - t12;
t25 = -pkin(3) * t115 - t29;
t125 = -qJ(5) * t48 + t25;
t5 = t159 * t47 + t125;
t1 = -t116 * t5 + t119 * t3;
t2 = t116 * t3 + t119 * t5;
t130 = t1 * t119 + t2 * t116;
t15 = -mrSges(7,2) * t48 + mrSges(7,3) * t32;
t16 = mrSges(7,1) * t48 - mrSges(7,3) * t33;
t129 = t116 * t15 + t119 * t16;
t63 = -t120 * t159 + t132;
t71 = t156 * t117;
t39 = -t116 * t63 + t119 * t71;
t40 = t116 * t71 + t119 * t63;
t128 = t116 * t40 + t119 * t39;
t127 = t170 * t61 + (-Ifges(6,4) + Ifges(5,5)) * t48 + (Ifges(6,5) - Ifges(5,6)) * t47;
t126 = Ifges(3,5) * t142 + Ifges(3,6) * t141 + t169 * t115 + t165;
t64 = Ifges(7,3) * t117 + (-Ifges(7,5) * t116 - Ifges(7,6) * t119) * t120;
t123 = qJ(5) ^ 2;
t85 = Ifges(5,1) * t117 + Ifges(5,4) * t120;
t84 = Ifges(7,1) * t119 - t145;
t83 = Ifges(5,4) * t117 + Ifges(5,2) * t120;
t82 = -Ifges(7,2) * t116 + t144;
t80 = -Ifges(6,2) * t117 - Ifges(6,6) * t120;
t79 = -Ifges(6,6) * t117 - Ifges(6,3) * t120;
t77 = -t120 * mrSges(5,1) + t117 * mrSges(5,2);
t76 = t120 * mrSges(6,2) - t117 * mrSges(6,3);
t74 = mrSges(7,1) * t117 + mrSges(7,3) * t140;
t72 = t156 * t120;
t70 = -mrSges(7,1) * t139 + mrSges(7,2) * t140;
t69 = t132 - t154;
t65 = Ifges(7,6) * t117 + (-Ifges(7,2) * t119 - t145) * t120;
t56 = t62 * mrSges(4,2);
t51 = mrSges(4,1) * t115 - mrSges(4,3) * t62;
t50 = -mrSges(4,2) * t115 - mrSges(4,3) * t61;
t24 = -mrSges(6,2) * t47 - mrSges(6,3) * t48;
t23 = mrSges(5,1) * t47 + mrSges(5,2) * t48;
t20 = Ifges(5,1) * t48 - Ifges(5,4) * t47 + Ifges(5,5) * t61;
t19 = Ifges(5,4) * t48 - Ifges(5,2) * t47 + Ifges(5,6) * t61;
t18 = -Ifges(6,2) * t48 + Ifges(6,6) * t47 + t153;
t17 = -Ifges(6,6) * t48 + Ifges(6,3) * t47 + t152;
t14 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t11 = pkin(4) * t47 + t125;
t10 = -pkin(4) * t61 - t12;
t8 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t48;
t7 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t48;
t4 = -pkin(5) * t47 - t9;
t21 = [m(4) * (t29 ^ 2 + t30 ^ 2 + t73 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t25 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t9 ^ 2) + (t20 + t6 - t18) * t48 + (mrSges(4,1) * t163 - 0.2e1 * Ifges(4,4) * t62 + Ifges(4,2) * t61 + t127) * t61 + (-t19 + t17) * t47 + ((t118 * Ifges(3,5) + t121 * Ifges(3,6)) * t115 + 0.2e1 * (-t67 * t118 + t68 * t121) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t121 + mrSges(3,2) * t118) + t118 * (Ifges(3,1) * t118 + Ifges(3,4) * t121) + t121 * (Ifges(3,4) * t118 + Ifges(3,2) * t121) + m(3) * pkin(1) ^ 2) * t113) * t113 + 0.2e1 * t30 * t50 + 0.2e1 * t29 * t51 + t32 * t7 + t33 * t8 + 0.2e1 * t9 * t35 + 0.2e1 * t10 * t36 + 0.2e1 * t13 * t37 + 0.2e1 * t12 * t38 + 0.2e1 * t2 * t15 + 0.2e1 * t1 * t16 + 0.2e1 * t11 * t24 + 0.2e1 * t25 * t23 + 0.2e1 * t4 * t14 + m(3) * (t67 ^ 2 + t68 ^ 2) + t56 * t163 + (t126 - 0.2e1 * t150 + 0.2e1 * t151 + t165) * t115 + Ifges(4,1) * t62 ^ 2 + Ifges(2,3); m(5) * (t25 * t99 + (-t12 * t117 + t13 * t120) * t98) + m(6) * (t11 * t69 + (t10 * t117 - t9 * t120) * t98) + (t7 * t157 + t8 * t158 + t13 * mrSges(5,3) - t9 * mrSges(6,1) - t152 / 0.2e1 - t17 / 0.2e1 + t19 / 0.2e1 + t149 * t98) * t120 + (-t12 * mrSges(5,3) + t10 * mrSges(6,1) - t153 / 0.2e1 - t18 / 0.2e1 + t20 / 0.2e1 + t6 / 0.2e1 + t148 * t98) * t117 + t61 * t137 / 0.2e1 + t99 * t23 + t1 * t74 + t2 * t75 + t11 * t76 + t25 * t77 + t69 * t24 - t4 * t70 + t72 * t14 + t39 * t16 + t40 * t15 + t29 * mrSges(4,1) - t30 * mrSges(4,2) + (m(4) * (t112 * t30 + t114 * t29) + t112 * t50 + t114 * t51) * pkin(2) + m(7) * (t1 * t39 + t2 * t40 + t4 * t72) + (-t80 / 0.2e1 + t85 / 0.2e1 + t64 / 0.2e1) * t48 + (t79 / 0.2e1 - t83 / 0.2e1) * t47 + t126 + t33 * t161 + t65 * t162 + t151 - t150; 0.2e1 * t39 * t74 + 0.2e1 * t40 * t75 + 0.2e1 * t69 * t76 - 0.2e1 * t72 * t70 + 0.2e1 * t99 * t77 + (t64 + t85 - t80) * t117 + (-t116 * t66 - t119 * t65 - t79 + t83) * t120 + m(7) * (t39 ^ 2 + t40 ^ 2 + t72 ^ 2) + m(6) * (t69 ^ 2 + t146) + m(5) * (t99 ^ 2 + t146) + (mrSges(6,1) + mrSges(5,3)) * t98 * t167 + (0.2e1 * t114 * mrSges(4,1) - 0.2e1 * t112 * mrSges(4,2) + m(4) * (t112 ^ 2 + t114 ^ 2) * pkin(2)) * pkin(2) + t169; t61 * mrSges(4,1) + t56 + (t14 + t149) * t117 + (-t129 - t148) * t120 + m(7) * (t117 * t4 - t120 * t130) + m(6) * (-t10 * t120 - t117 * t9) + m(5) * (t117 * t13 + t12 * t120) + m(4) * t73; m(7) * (t117 * t72 - t120 * t128) - t75 * t140 - t74 * t139 - t117 * t70; m(4) + m(7) * (t110 * t136 + t108) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t167; m(7) * (qJ(5) * t4 - t130 * t159) + t4 * t78 + t48 * t160 + t82 * t162 + t33 * t84 / 0.2e1 - pkin(4) * t36 - t9 * mrSges(6,3) + t10 * mrSges(6,2) + t12 * mrSges(5,1) - t13 * mrSges(5,2) + t127 + m(6) * (-pkin(4) * t10 - qJ(5) * t9) + (-t35 + t14) * qJ(5) + (-t159 * t15 - t2 * mrSges(7,3) - t7 / 0.2e1) * t116 + (-t159 * t16 - t1 * mrSges(7,3) + t8 / 0.2e1) * t119; -t159 * t143 - t74 * t138 + m(7) * (qJ(5) * t72 - t128 * t159) + t119 * t161 + t65 * t158 + t72 * t78 - qJ(5) * t70 - t128 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) - Ifges(6,4) + t160) * t117 + (qJ(5) * mrSges(6,1) + t82 * t157 + t84 * t158 - Ifges(6,5)) * t120 + ((m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t120 + (-mrSges(5,1) + mrSges(6,2) - t166) * t117) * t98 + t137; (mrSges(5,1) + t164) * t120 + m(6) * (t102 + t154) + m(7) * (t120 * t133 + t102) + (-mrSges(5,2) + t147) * t117; -0.2e1 * pkin(4) * mrSges(6,2) - t116 * t82 + t119 * t84 + m(6) * (pkin(4) ^ 2 + t123) + m(7) * (t136 * t159 ^ 2 + t123) + 0.2e1 * mrSges(7,3) * t133 + 0.2e1 * qJ(5) * t147 + t170; m(6) * t10 + m(7) * t130 + t129 + t36; m(7) * t128 + t143 + t119 * t74 + (m(6) * t98 + mrSges(6,1)) * t117; -t171 * t120; -t134 * t159 - t164 - t166; t171; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t6; mrSges(7,1) * t39 - mrSges(7,2) * t40 + t64; t70; -mrSges(7,1) * t138 + t105 + (mrSges(7,2) * t159 - Ifges(7,6)) * t116; mrSges(7,1) * t119 - t116 * mrSges(7,2); Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
