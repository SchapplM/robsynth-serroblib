% Calculate joint inertia matrix for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:58
% EndTime: 2019-03-08 22:23:01
% DurationCPUTime: 1.32s
% Computational Cost: add. (2478->330), mult. (5950->489), div. (0->0), fcn. (6774->14), ass. (0->129)
t125 = sin(qJ(5));
t158 = cos(qJ(5));
t118 = sin(pkin(13));
t121 = cos(pkin(13));
t122 = cos(pkin(7));
t119 = sin(pkin(7));
t126 = sin(qJ(3));
t145 = t119 * t126;
t77 = -t118 * t145 + t121 * t122;
t79 = t118 * t122 + t121 * t145;
t48 = t125 * t79 - t158 * t77;
t49 = t125 * t77 + t158 * t79;
t21 = t48 * mrSges(6,1) + t49 * mrSges(6,2);
t51 = -t77 * mrSges(5,1) + t79 * mrSges(5,2);
t172 = t21 + t51;
t88 = t118 * t125 - t158 * t121;
t89 = t158 * t118 + t125 * t121;
t60 = t88 * mrSges(6,1) + t89 * mrSges(6,2);
t93 = -t121 * mrSges(5,1) + t118 * mrSges(5,2);
t171 = t60 + t93;
t129 = cos(qJ(3));
t144 = t119 * t129;
t170 = m(7) * pkin(11) + mrSges(7,3);
t124 = sin(qJ(6));
t128 = cos(qJ(6));
t97 = -mrSges(7,1) * t128 + mrSges(7,2) * t124;
t169 = -m(7) * pkin(5) - mrSges(6,1) + t97;
t120 = sin(pkin(6));
t123 = cos(pkin(6));
t127 = sin(qJ(2));
t130 = cos(qJ(2));
t143 = t122 * t130;
t59 = t123 * t145 + (t126 * t143 + t127 * t129) * t120;
t78 = -t119 * t120 * t130 + t122 * t123;
t29 = -t118 * t59 + t121 * t78;
t30 = t118 * t78 + t121 * t59;
t13 = t125 * t30 - t158 * t29;
t168 = t13 ^ 2;
t57 = -t123 * t144 + (t126 * t127 - t129 * t143) * t120;
t56 = t57 ^ 2;
t153 = pkin(10) + qJ(4);
t139 = t153 * t118;
t94 = t153 * t121;
t63 = t125 * t94 + t158 * t139;
t167 = t63 ^ 2;
t166 = 0.2e1 * t63;
t33 = -t124 * t49 - t128 * t144;
t165 = t33 / 0.2e1;
t34 = -t124 * t144 + t128 * t49;
t164 = t34 / 0.2e1;
t98 = Ifges(7,5) * t124 + Ifges(7,6) * t128;
t163 = t98 / 0.2e1;
t162 = -t124 / 0.2e1;
t161 = t124 / 0.2e1;
t160 = t128 / 0.2e1;
t157 = pkin(2) * t129;
t156 = pkin(11) * t124;
t155 = pkin(11) * t128;
t154 = t13 * t63;
t16 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t39 = -mrSges(6,1) * t144 - mrSges(6,3) * t49;
t152 = t16 - t39;
t82 = t122 * t126 * pkin(2) + pkin(9) * t144;
t71 = qJ(4) * t122 + t82;
t72 = (-pkin(3) * t129 - qJ(4) * t126 - pkin(2)) * t119;
t40 = -t118 * t71 + t121 * t72;
t23 = -pkin(4) * t144 - pkin(10) * t79 + t40;
t41 = t118 * t72 + t121 * t71;
t27 = pkin(10) * t77 + t41;
t11 = t125 * t23 + t158 * t27;
t151 = -Ifges(6,5) * t49 + Ifges(6,6) * t48;
t150 = Ifges(6,5) * t89 - Ifges(6,6) * t88;
t149 = Ifges(7,4) * t124;
t148 = Ifges(7,4) * t128;
t147 = t124 * t89;
t146 = t128 * t89;
t142 = t118 ^ 2 + t121 ^ 2;
t141 = t124 ^ 2 + t128 ^ 2;
t7 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t48;
t140 = Ifges(4,5) * t145 + Ifges(4,6) * t144 + Ifges(4,3) * t122;
t108 = -pkin(4) * t121 - pkin(3);
t104 = pkin(9) * t145;
t74 = t104 + (-pkin(3) - t157) * t122;
t50 = -pkin(4) * t77 + t74;
t12 = pkin(5) * t48 - pkin(11) * t49 + t50;
t6 = -pkin(11) * t144 + t11;
t1 = t12 * t128 - t124 * t6;
t2 = t12 * t124 + t128 * t6;
t138 = -t1 * t124 + t128 * t2;
t136 = mrSges(7,1) * t124 + mrSges(7,2) * t128;
t135 = -t118 * t29 + t121 * t30;
t134 = -t40 * t118 + t41 * t121;
t35 = Ifges(7,5) * t146 - Ifges(7,6) * t147 + Ifges(7,3) * t88;
t10 = -t125 * t27 + t158 * t23;
t100 = Ifges(7,1) * t124 + t148;
t99 = Ifges(7,2) * t128 + t149;
t96 = Ifges(5,1) * t118 + Ifges(5,4) * t121;
t95 = Ifges(5,4) * t118 + Ifges(5,2) * t121;
t91 = -mrSges(4,2) * t122 + mrSges(4,3) * t144;
t90 = mrSges(4,1) * t122 - mrSges(4,3) * t145;
t81 = t122 * t157 - t104;
t80 = (-mrSges(4,1) * t129 + mrSges(4,2) * t126) * t119;
t67 = -mrSges(5,1) * t144 - mrSges(5,3) * t79;
t66 = mrSges(5,2) * t144 + mrSges(5,3) * t77;
t65 = -t125 * t139 + t158 * t94;
t62 = Ifges(6,1) * t89 - Ifges(6,4) * t88;
t61 = Ifges(6,4) * t89 - Ifges(6,2) * t88;
t55 = mrSges(7,1) * t88 - mrSges(7,3) * t146;
t54 = -mrSges(7,2) * t88 - mrSges(7,3) * t147;
t53 = pkin(5) * t88 - pkin(11) * t89 + t108;
t52 = t136 * t89;
t43 = Ifges(5,1) * t79 + Ifges(5,4) * t77 - Ifges(5,5) * t144;
t42 = Ifges(5,4) * t79 + Ifges(5,2) * t77 - Ifges(5,6) * t144;
t38 = mrSges(6,2) * t144 - mrSges(6,3) * t48;
t37 = Ifges(7,5) * t88 + (Ifges(7,1) * t128 - t149) * t89;
t36 = Ifges(7,6) * t88 + (-Ifges(7,2) * t124 + t148) * t89;
t25 = t124 * t53 + t128 * t65;
t24 = -t124 * t65 + t128 * t53;
t20 = Ifges(6,1) * t49 - Ifges(6,4) * t48 - Ifges(6,5) * t144;
t19 = Ifges(6,4) * t49 - Ifges(6,2) * t48 - Ifges(6,6) * t144;
t18 = mrSges(7,1) * t48 - mrSges(7,3) * t34;
t17 = -mrSges(7,2) * t48 + mrSges(7,3) * t33;
t15 = t125 * t29 + t158 * t30;
t9 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t48;
t8 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t48;
t5 = pkin(5) * t144 - t10;
t4 = t124 * t57 + t128 * t15;
t3 = -t124 * t15 + t128 * t57;
t14 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t168) + m(6) * (t15 ^ 2 + t168 + t56) + m(5) * (t29 ^ 2 + t30 ^ 2 + t56) + m(4) * (t59 ^ 2 + t78 ^ 2 + t56) + m(3) * (t123 ^ 2 + (t127 ^ 2 + t130 ^ 2) * t120 ^ 2); t15 * t38 + t4 * t17 + t3 * t18 + t29 * t67 + t30 * t66 + t59 * t91 + t78 * t80 + t152 * t13 + (mrSges(3,1) * t130 - mrSges(3,2) * t127) * t120 + (-t90 + t172) * t57 + m(7) * (t1 * t3 + t13 * t5 + t2 * t4) + m(6) * (-t10 * t13 + t11 * t15 + t50 * t57) + m(5) * (t29 * t40 + t30 * t41 + t57 * t74) + m(4) * (-pkin(2) * t119 * t78 - t57 * t81 + t59 * t82); t122 * t140 + 0.2e1 * t81 * t90 + 0.2e1 * t82 * t91 + t77 * t42 + t79 * t43 + 0.2e1 * t41 * t66 + 0.2e1 * t40 * t67 + 0.2e1 * t74 * t51 + t49 * t20 + 0.2e1 * t50 * t21 + t33 * t8 + t34 * t9 + 0.2e1 * t11 * t38 + 0.2e1 * t10 * t39 + 0.2e1 * t5 * t16 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + Ifges(3,3) + (t7 - t19) * t48 + (-0.2e1 * pkin(2) * t80 + (Ifges(4,1) * t145 + Ifges(4,5) * t122) * t126 + (-Ifges(5,5) * t79 + Ifges(4,6) * t122 - Ifges(5,6) * t77 + (0.2e1 * Ifges(4,4) * t126 + (Ifges(4,2) + Ifges(5,3) + Ifges(6,3)) * t129) * t119 + t151) * t129) * t119 + m(4) * (pkin(2) ^ 2 * t119 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (t40 ^ 2 + t41 ^ 2 + t74 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t50 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); -t59 * mrSges(4,2) + t13 * t52 + t3 * t55 + t4 * t54 + (t13 * t89 - t15 * t88) * mrSges(6,3) + t135 * mrSges(5,3) + (-mrSges(4,1) + t171) * t57 + m(7) * (t24 * t3 + t25 * t4 + t154) + m(6) * (t108 * t57 + t15 * t65 + t154) + m(5) * (-pkin(3) * t57 + t135 * qJ(4)); -(Ifges(5,5) * t118 + Ifges(5,6) * t121 + t150) * t144 / 0.2e1 + (t20 / 0.2e1 + t9 * t160 + t8 * t162 - t10 * mrSges(6,3)) * t89 + t152 * t63 + t118 * t43 / 0.2e1 + t121 * t42 / 0.2e1 + t108 * t21 + t74 * t93 + t77 * t95 / 0.2e1 + t79 * t96 / 0.2e1 + t81 * mrSges(4,1) - t82 * mrSges(4,2) + t65 * t38 - pkin(3) * t51 + t5 * t52 + t2 * t54 + t1 * t55 + t50 * t60 + t49 * t62 / 0.2e1 + t134 * mrSges(5,3) + m(5) * (-pkin(3) * t74 + t134 * qJ(4)) + t24 * t18 + t25 * t17 + t140 + (-t118 * t67 + t121 * t66) * qJ(4) + m(6) * (-t10 * t63 + t108 * t50 + t11 * t65) + m(7) * (t1 * t24 + t2 * t25 + t5 * t63) + (t7 / 0.2e1 - t19 / 0.2e1 - t11 * mrSges(6,3)) * t88 + (-t61 / 0.2e1 + t35 / 0.2e1) * t48 + t37 * t164 + t36 * t165; -0.2e1 * pkin(3) * t93 + 0.2e1 * t108 * t60 + t118 * t96 + t121 * t95 + 0.2e1 * t24 * t55 + 0.2e1 * t25 * t54 + t52 * t166 + Ifges(4,3) + 0.2e1 * t142 * qJ(4) * mrSges(5,3) + (-0.2e1 * t65 * mrSges(6,3) + t35 - t61) * t88 + m(7) * (t24 ^ 2 + t25 ^ 2 + t167) + m(6) * (t108 ^ 2 + t65 ^ 2 + t167) + m(5) * (t142 * qJ(4) ^ 2 + pkin(3) ^ 2) + (mrSges(6,3) * t166 - t124 * t36 + t128 * t37 + t62) * t89; m(7) * (t124 * t4 + t128 * t3) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t57; t124 * t17 + t128 * t18 + m(7) * (t1 * t128 + t124 * t2) + m(6) * t50 + m(5) * t74 + t172; -m(5) * pkin(3) + t124 * t54 + t128 * t55 + m(7) * (t124 * t25 + t128 * t24) + m(6) * t108 + t171; m(7) * t141 + m(5) + m(6); -t15 * mrSges(6,2) + t170 * (-t124 * t3 + t128 * t4) + t169 * t13; m(7) * (-pkin(5) * t5 + t138 * pkin(11)) + t17 * t155 - t18 * t156 + t8 * t160 + t9 * t161 - pkin(5) * t16 + t5 * t97 + t100 * t164 + t99 * t165 + t48 * t163 - Ifges(6,3) * t144 - t11 * mrSges(6,2) + t10 * mrSges(6,1) + t138 * mrSges(7,3) - t151; -pkin(5) * t52 + t54 * t155 - t55 * t156 + t37 * t161 + t36 * t160 + t88 * t163 - t65 * mrSges(6,2) + (t100 * t160 + t99 * t162) * t89 + t150 + t170 * (-t124 * t24 + t128 * t25) + t169 * t63; 0; Ifges(6,3) + m(7) * (t141 * pkin(11) ^ 2 + pkin(5) ^ 2) + t124 * t100 + t128 * t99 - 0.2e1 * pkin(5) * t97 + 0.2e1 * t141 * pkin(11) * mrSges(7,3); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t24 - mrSges(7,2) * t25 + t35; -t97; -t136 * pkin(11) + t98; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
