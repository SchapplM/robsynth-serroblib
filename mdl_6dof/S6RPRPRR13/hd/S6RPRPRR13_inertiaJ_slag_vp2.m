% Calculate joint inertia matrix for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:10:47
% EndTime: 2018-11-23 16:10:49
% DurationCPUTime: 1.54s
% Computational Cost: add. (3416->327), mult. (8990->471), div. (0->0), fcn. (10033->12), ass. (0->138)
t178 = Ifges(5,1) + Ifges(4,3);
t177 = -Ifges(5,4) + Ifges(4,5);
t176 = Ifges(5,5) - Ifges(4,6);
t165 = (pkin(3) + pkin(10));
t175 = 2 * t165;
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t116 = cos(pkin(6));
t119 = sin(qJ(3));
t122 = cos(qJ(3));
t115 = cos(pkin(7));
t113 = sin(pkin(6));
t114 = cos(pkin(12));
t148 = t113 * t114;
t139 = t115 * t148;
t112 = sin(pkin(7));
t149 = t112 * t122;
t111 = sin(pkin(12));
t151 = t111 * t113;
t55 = -t116 * t149 + t119 * t151 - t122 * t139;
t73 = -t112 * t148 + t115 * t116;
t40 = t118 * t55 + t121 * t73;
t147 = t115 * t119;
t150 = t112 * t119;
t56 = t116 * t150 + (t111 * t122 + t114 * t147) * t113;
t29 = -t117 * t40 + t120 * t56;
t39 = t118 * t73 - t121 * t55;
t17 = -mrSges(7,2) * t39 + mrSges(7,3) * t29;
t30 = t117 * t56 + t120 * t40;
t18 = mrSges(7,1) * t39 - mrSges(7,3) * t30;
t174 = -t117 * t18 + t120 * t17;
t173 = m(7) * pkin(11) + mrSges(7,3);
t160 = pkin(1) * t116;
t75 = qJ(2) * t148 + t111 * t160;
t49 = (t112 * t116 + t139) * pkin(9) + t75;
t95 = t114 * t160;
t57 = t116 * pkin(2) + t95 + (-pkin(9) * t115 - qJ(2)) * t151;
t62 = (-pkin(9) * t111 * t112 - pkin(2) * t114 - pkin(1)) * t113;
t25 = -t119 * t49 + (t112 * t62 + t115 * t57) * t122;
t84 = -mrSges(7,1) * t120 + mrSges(7,2) * t117;
t172 = m(7) * pkin(5) + mrSges(6,1) - t84;
t76 = t115 * t118 + t121 * t149;
t171 = t76 ^ 2;
t170 = 2 * mrSges(3,1);
t169 = 0.2e1 * t116;
t168 = t29 / 0.2e1;
t167 = t30 / 0.2e1;
t86 = Ifges(7,5) * t117 + Ifges(7,6) * t120;
t166 = t86 / 0.2e1;
t164 = -t117 / 0.2e1;
t163 = t117 / 0.2e1;
t162 = t120 / 0.2e1;
t161 = Ifges(6,6) * t56;
t11 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t32 = mrSges(6,1) * t56 - mrSges(6,3) * t40;
t159 = t11 - t32;
t13 = t56 * pkin(4) - t165 * t73 - t25;
t35 = -t112 * t57 + t115 * t62;
t127 = -t56 * qJ(4) + t35;
t15 = t165 * t55 + t127;
t6 = t118 * t13 + t121 * t15;
t24 = mrSges(6,1) * t39 + mrSges(6,2) * t40;
t43 = mrSges(5,1) * t55 - mrSges(5,3) * t73;
t158 = -t43 + t24;
t85 = mrSges(6,1) * t118 + mrSges(6,2) * t121;
t157 = t85 + mrSges(5,3);
t156 = Ifges(7,4) * t117;
t155 = Ifges(7,4) * t120;
t152 = t121 * t76;
t146 = t117 * t121;
t145 = t118 * t165;
t144 = t120 * t121;
t143 = t121 * t165;
t142 = t117 ^ 2 + t120 ^ 2;
t107 = t118 ^ 2;
t109 = t121 ^ 2;
t141 = t107 + t109;
t8 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t39;
t140 = Ifges(6,5) * t40 - Ifges(6,6) * t39 + Ifges(6,3) * t56;
t26 = t122 * t49 + t57 * t147 + t62 * t150;
t44 = t56 * mrSges(5,1) + t73 * mrSges(5,2);
t138 = t141 * mrSges(6,3);
t22 = -t73 * qJ(4) - t26;
t4 = pkin(11) * t56 + t6;
t16 = -pkin(4) * t55 - t22;
t7 = pkin(5) * t39 - pkin(11) * t40 + t16;
t1 = -t117 * t4 + t120 * t7;
t2 = t117 * t7 + t120 * t4;
t136 = -t1 * t117 + t120 * t2;
t5 = -t118 * t15 + t121 * t13;
t135 = t118 * t6 + t121 * t5;
t134 = mrSges(7,1) * t117 + mrSges(7,2) * t120;
t78 = t115 * t121 - t118 * t149;
t59 = -t117 * t78 + t120 * t150;
t60 = t117 * t150 + t120 * t78;
t132 = -t117 * t59 + t120 * t60;
t83 = pkin(5) * t118 - pkin(11) * t121 + qJ(4);
t63 = t117 * t145 + t120 * t83;
t64 = t117 * t83 - t120 * t145;
t131 = -t117 * t63 + t120 * t64;
t81 = -mrSges(7,2) * t118 - mrSges(7,3) * t146;
t82 = mrSges(7,1) * t118 - mrSges(7,3) * t144;
t130 = -t117 * t82 + t120 * t81;
t129 = t118 * t78 - t152;
t128 = t176 * t55 + t177 * t56 + t178 * t73;
t65 = Ifges(7,5) * t144 - Ifges(7,6) * t146 + Ifges(7,3) * t118;
t124 = qJ(4) ^ 2;
t110 = t165 ^ 2;
t105 = t112 ^ 2;
t104 = Ifges(6,5) * t121;
t100 = t109 * t165;
t99 = t109 * t110;
t98 = t105 * t119 ^ 2;
t96 = qJ(4) * t150;
t93 = mrSges(3,2) * t151;
t90 = Ifges(6,1) * t121 - Ifges(6,4) * t118;
t89 = Ifges(7,1) * t117 + t155;
t88 = Ifges(6,4) * t121 - Ifges(6,2) * t118;
t87 = Ifges(7,2) * t120 + t156;
t79 = t134 * t121;
t74 = -qJ(2) * t151 + t95;
t67 = Ifges(7,5) * t118 + (Ifges(7,1) * t120 - t156) * t121;
t66 = Ifges(7,6) * t118 + (-Ifges(7,2) * t117 + t155) * t121;
t42 = mrSges(4,1) * t73 - mrSges(4,3) * t56;
t41 = -mrSges(4,2) * t73 - mrSges(4,3) * t55;
t34 = -mrSges(5,2) * t55 - mrSges(5,3) * t56;
t33 = mrSges(4,1) * t55 + mrSges(4,2) * t56;
t31 = -mrSges(6,2) * t56 - mrSges(6,3) * t39;
t23 = -t73 * pkin(3) - t25;
t21 = pkin(3) * t55 + t127;
t20 = Ifges(6,1) * t40 - Ifges(6,4) * t39 + Ifges(6,5) * t56;
t19 = Ifges(6,4) * t40 - Ifges(6,2) * t39 + t161;
t10 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t39;
t9 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t39;
t3 = -pkin(5) * t56 - t5;
t12 = [t40 * t20 + 0.2e1 * t26 * t41 + 0.2e1 * t25 * t42 + 0.2e1 * t22 * t43 + 0.2e1 * t23 * t44 + 0.2e1 * t21 * t34 + 0.2e1 * t35 * t33 + t29 * t9 + t30 * t10 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t16 * t24 + 0.2e1 * t3 * t11 + m(3) * (pkin(1) ^ 2 * t113 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t25 ^ 2 + t26 ^ 2 + t35 ^ 2) + m(6) * (t16 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (t8 - t19) * t39 + t128 * t73 + (-0.2e1 * mrSges(3,2) * t75 + Ifges(3,3) * t116 + t74 * t170) * t116 + (-0.2e1 * pkin(1) * t93 + (-0.2e1 * mrSges(3,3) * t74 + Ifges(3,1) * t151 + Ifges(3,5) * t169) * t111 + (0.2e1 * t75 * mrSges(3,3) + Ifges(3,6) * t169 + (0.2e1 * Ifges(3,4) * t111 + Ifges(3,2) * t114 + pkin(1) * t170) * t113) * t114) * t113 + Ifges(2,3) + ((Ifges(4,2) + Ifges(5,3)) * t55 + t176 * t73) * t55 + ((Ifges(4,1) + Ifges(5,2)) * t56 + t177 * t73 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t55 + t140) * t56; t60 * t17 + t59 * t18 + t78 * t31 + t93 + t159 * t76 + (t33 + t34) * t115 + (-m(3) * pkin(1) - mrSges(3,1) * t114) * t113 + ((t42 - t44) * t122 + (t41 + t158) * t119) * t112 + m(7) * (t1 * t59 + t2 * t60 + t3 * t76) + m(6) * (t16 * t150 - t5 * t76 + t6 * t78) + m(5) * (t115 * t21 + (-t119 * t22 - t122 * t23) * t112) + m(4) * (t115 * t35 + (t119 * t26 + t122 * t25) * t112); m(3) + m(7) * (t59 ^ 2 + t60 ^ 2 + t171) + m(6) * (t78 ^ 2 + t171 + t98) + 0.2e1 * (m(5) / 0.2e1 + m(4) / 0.2e1) * (t105 * t122 ^ 2 + t115 ^ 2 + t98); (t10 * t162 + t9 * t164 - t5 * mrSges(6,3) + t20 / 0.2e1 + t159 * t165) * t121 + t3 * t79 + t2 * t81 + t1 * t82 + t16 * t85 + t40 * t90 / 0.2e1 + t63 * t18 + t64 * t17 + t66 * t168 + t67 * t167 - pkin(3) * t44 + t25 * mrSges(4,1) - t26 * mrSges(4,2) - t22 * mrSges(5,3) + t23 * mrSges(5,2) + (-t161 / 0.2e1 - t19 / 0.2e1 + t8 / 0.2e1 - t6 * mrSges(6,3) - t165 * t31) * t118 + t158 * qJ(4) + m(5) * (-pkin(3) * t23 - qJ(4) * t22) + m(7) * (t1 * t63 + t3 * t143 + t2 * t64) + (-t88 / 0.2e1 + t65 / 0.2e1) * t39 + t56 * t104 / 0.2e1 + m(6) * (qJ(4) * t16 - t135 * t165) + t128; t59 * t82 + t60 * t81 + t76 * t79 - t129 * mrSges(6,3) + ((mrSges(4,1) - mrSges(5,2)) * t122 + (-mrSges(4,2) + t157) * t119) * t112 + m(7) * (t76 * t143 + t59 * t63 + t60 * t64) + m(6) * (-t129 * t165 + t96) + m(5) * (pkin(3) * t149 + t96); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t63 * t82 + 0.2e1 * t64 * t81 + (t65 - t88) * t118 + m(7) * (t63 ^ 2 + t64 ^ 2 + t99) + m(6) * (t107 * t110 + t124 + t99) + m(5) * (pkin(3) ^ 2 + t124) + (-t117 * t66 + t120 * t67 + t79 * t175 + t90) * t121 + 0.2e1 * t157 * qJ(4) + t138 * t175 + t178; -t159 * t121 + (t31 + t174) * t118 + m(7) * (t136 * t118 - t121 * t3) + m(6) * t135 + m(5) * t23 + t44; -m(5) * t149 + m(7) * (t132 * t118 - t152) + m(6) * t129; -m(5) * pkin(3) - t121 * t79 + mrSges(5,2) + t130 * t118 - t138 + m(7) * (t131 * t118 - t100) + m(6) * (-t107 * t165 - t100); m(5) + m(6) * t141 + m(7) * (t142 * t107 + t109); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t136 * mrSges(7,3) + t10 * t163 + t9 * t162 + t39 * t166 + t89 * t167 + t87 * t168 + t3 * t84 + t140 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t136 + t174) * pkin(11); -t78 * mrSges(6,2) + t173 * t132 - t172 * t76; -pkin(5) * t79 + t67 * t163 + t66 * t162 + t104 + (m(7) * t131 + t130) * pkin(11) + (t89 * t162 + t87 * t164 - t165 * t172) * t121 + t131 * mrSges(7,3) + (mrSges(6,2) * t165 - Ifges(6,6) + t166) * t118; t172 * t121 + (t142 * t173 - mrSges(6,2)) * t118; Ifges(6,3) - 0.2e1 * pkin(5) * t84 + t117 * t89 + t120 * t87 + m(7) * (t142 * pkin(11) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t142 * pkin(11) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t59 - mrSges(7,2) * t60; mrSges(7,1) * t63 - mrSges(7,2) * t64 + t65; -t134 * t118; -t134 * pkin(11) + t86; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
