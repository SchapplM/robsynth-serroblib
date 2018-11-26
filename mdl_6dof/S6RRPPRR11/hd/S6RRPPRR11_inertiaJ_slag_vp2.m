% Calculate joint inertia matrix for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:56:02
% EndTime: 2018-11-23 16:56:03
% DurationCPUTime: 1.34s
% Computational Cost: add. (2412->320), mult. (5123->449), div. (0->0), fcn. (5480->10), ass. (0->123)
t172 = Ifges(4,1) + Ifges(3,3);
t121 = sin(qJ(6));
t124 = cos(qJ(6));
t116 = sin(pkin(11));
t118 = cos(pkin(11));
t122 = sin(qJ(5));
t159 = cos(qJ(5));
t79 = t116 * t122 - t118 * t159;
t147 = t121 * t79;
t80 = t116 * t159 + t122 * t118;
t48 = -mrSges(7,2) * t80 + mrSges(7,3) * t147;
t144 = t124 * t79;
t49 = mrSges(7,1) * t80 + mrSges(7,3) * t144;
t171 = -t121 * t49 + t124 * t48;
t117 = sin(pkin(6));
t123 = sin(qJ(2));
t143 = t117 * t123;
t119 = cos(pkin(6));
t125 = cos(qJ(2));
t142 = t117 * t125;
t68 = -t116 * t119 - t118 * t142;
t69 = -t116 * t142 + t118 * t119;
t42 = t122 * t68 + t159 * t69;
t27 = -t121 * t42 + t124 * t143;
t41 = t122 * t69 - t159 * t68;
t12 = -mrSges(7,2) * t41 + mrSges(7,3) * t27;
t28 = t121 * t143 + t124 * t42;
t13 = mrSges(7,1) * t41 - mrSges(7,3) * t28;
t170 = t124 * t12 - t121 * t13;
t88 = -mrSges(7,1) * t124 + mrSges(7,2) * t121;
t169 = -m(7) * pkin(5) - mrSges(6,1) + t88;
t120 = -pkin(2) - qJ(4);
t154 = -pkin(9) + t120;
t136 = t154 * t118;
t84 = t154 * t116;
t52 = t122 * t84 - t136 * t159;
t168 = t52 ^ 2;
t78 = t79 ^ 2;
t167 = -2 * mrSges(6,3);
t166 = t27 / 0.2e1;
t165 = t28 / 0.2e1;
t89 = Ifges(7,5) * t121 + Ifges(7,6) * t124;
t164 = t89 / 0.2e1;
t163 = t121 / 0.2e1;
t161 = -t124 / 0.2e1;
t160 = t124 / 0.2e1;
t158 = pkin(1) * t125;
t157 = t52 * t79;
t96 = pkin(8) * t143;
t70 = t119 * t158 - t96;
t156 = t70 * mrSges(3,1);
t71 = t119 * t123 * pkin(1) + pkin(8) * t142;
t155 = t71 * mrSges(3,2);
t11 = -mrSges(7,1) * t27 + mrSges(7,2) * t28;
t33 = mrSges(6,1) * t143 - mrSges(6,3) * t42;
t153 = t11 - t33;
t138 = -pkin(2) - t158;
t50 = pkin(3) * t143 + t96 + (-qJ(4) + t138) * t119;
t135 = -qJ(3) * t123 - pkin(1);
t58 = (t120 * t125 + t135) * t117;
t23 = -t116 * t58 + t118 * t50;
t17 = pkin(4) * t143 - pkin(9) * t69 + t23;
t24 = t116 * t50 + t118 * t58;
t19 = pkin(9) * t68 + t24;
t6 = t122 * t17 + t159 * t19;
t152 = -Ifges(6,5) * t79 - Ifges(6,6) * t80;
t151 = Ifges(7,4) * t121;
t150 = Ifges(7,4) * t124;
t83 = mrSges(4,1) * t143 + t119 * mrSges(4,2);
t100 = t116 * pkin(4) + qJ(3);
t85 = t116 * mrSges(5,1) + t118 * mrSges(5,2);
t141 = -t116 ^ 2 - t118 ^ 2;
t140 = t121 ^ 2 + t124 ^ 2;
t7 = Ifges(7,5) * t28 + Ifges(7,6) * t27 + Ifges(7,3) * t41;
t139 = Ifges(6,5) * t42 - Ifges(6,6) * t41 + Ifges(6,3) * t143;
t62 = -t119 * qJ(3) - t71;
t43 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t20 = t41 * mrSges(6,1) + t42 * mrSges(6,2);
t55 = t80 * mrSges(6,1) - t79 * mrSges(6,2);
t137 = m(5) * t141;
t134 = t141 * mrSges(5,3);
t59 = pkin(3) * t142 - t62;
t132 = Ifges(3,5) * t143 + Ifges(3,6) * t142 + t172 * t119;
t34 = -pkin(4) * t68 + t59;
t10 = pkin(5) * t41 - pkin(10) * t42 + t34;
t4 = pkin(10) * t143 + t6;
t1 = t10 * t124 - t121 * t4;
t2 = t10 * t121 + t124 * t4;
t131 = -t1 * t121 + t124 * t2;
t130 = -mrSges(7,1) * t121 - mrSges(7,2) * t124;
t129 = t24 * t116 + t23 * t118;
t44 = pkin(5) * t80 + pkin(10) * t79 + t100;
t54 = t122 * t136 + t159 * t84;
t21 = -t121 * t54 + t124 * t44;
t22 = t121 * t44 + t124 * t54;
t128 = -t121 * t21 + t124 * t22;
t29 = -Ifges(7,5) * t144 + Ifges(7,6) * t147 + Ifges(7,3) * t80;
t5 = -t122 * t19 + t159 * t17;
t126 = qJ(3) ^ 2;
t91 = Ifges(7,1) * t121 + t150;
t90 = Ifges(7,2) * t124 + t151;
t87 = Ifges(5,1) * t118 - Ifges(5,4) * t116;
t86 = Ifges(5,4) * t118 - Ifges(5,2) * t116;
t82 = -mrSges(4,1) * t142 - mrSges(4,3) * t119;
t77 = t80 ^ 2;
t65 = t119 * t138 + t96;
t63 = (-pkin(2) * t125 + t135) * t117;
t61 = mrSges(5,1) * t143 - mrSges(5,3) * t69;
t60 = -mrSges(5,2) * t143 + mrSges(5,3) * t68;
t57 = -Ifges(6,1) * t79 - Ifges(6,4) * t80;
t56 = -Ifges(6,4) * t79 - Ifges(6,2) * t80;
t45 = t130 * t79;
t36 = Ifges(5,1) * t69 + Ifges(5,4) * t68 + Ifges(5,5) * t143;
t35 = Ifges(5,4) * t69 + Ifges(5,2) * t68 + Ifges(5,6) * t143;
t32 = -mrSges(6,2) * t143 - mrSges(6,3) * t41;
t31 = Ifges(7,5) * t80 + (-Ifges(7,1) * t124 + t151) * t79;
t30 = Ifges(7,6) * t80 + (Ifges(7,2) * t121 - t150) * t79;
t16 = Ifges(6,1) * t42 - Ifges(6,4) * t41 + Ifges(6,5) * t143;
t15 = Ifges(6,4) * t42 - Ifges(6,2) * t41 + Ifges(6,6) * t143;
t9 = Ifges(7,1) * t28 + Ifges(7,4) * t27 + Ifges(7,5) * t41;
t8 = Ifges(7,4) * t28 + Ifges(7,2) * t27 + Ifges(7,6) * t41;
t3 = -pkin(5) * t143 - t5;
t14 = [(m(3) * pkin(1) ^ 2 * t117 + (0.2e1 * t63 * mrSges(4,2) + 0.2e1 * t71 * mrSges(3,3) + (-(2 * Ifges(4,5)) + Ifges(3,6)) * t119 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t125) * t117) * t125 + (-0.2e1 * t70 * mrSges(3,3) - 0.2e1 * t63 * mrSges(4,3) + Ifges(5,5) * t69 + Ifges(5,6) * t68 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t119 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1) + Ifges(5,3)) * t123 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t125) * t117 + t139) * t123) * t117 + m(3) * (t70 ^ 2 + t71 ^ 2) + (t132 - 0.2e1 * t155 + 0.2e1 * t156) * t119 + m(5) * (t23 ^ 2 + t24 ^ 2 + t59 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2 + t65 ^ 2) + m(6) * (t34 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + 0.2e1 * t62 * t82 + 0.2e1 * t65 * t83 + t68 * t35 + t69 * t36 + 0.2e1 * t59 * t43 + 0.2e1 * t24 * t60 + 0.2e1 * t23 * t61 + t42 * t16 + t27 * t8 + t28 * t9 + 0.2e1 * t6 * t32 + 0.2e1 * t5 * t33 + 0.2e1 * t34 * t20 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t1 * t13 + (t7 - t15) * t41 + Ifges(2,3); (-Ifges(4,4) * t123 - Ifges(4,5) * t125 + (Ifges(5,5) * t118 - Ifges(5,6) * t116 + t152) * t123 / 0.2e1) * t117 + (-t16 / 0.2e1 + t9 * t161 + t8 * t163 + t5 * mrSges(6,3)) * t79 + t153 * t52 + m(5) * (qJ(3) * t59 + t120 * t129) + t59 * t85 + t68 * t86 / 0.2e1 + t69 * t87 / 0.2e1 + t100 * t20 - pkin(2) * t83 + t42 * t57 / 0.2e1 - t62 * mrSges(4,3) + t65 * mrSges(4,2) + t3 * t45 + t2 * t48 + t1 * t49 + t54 * t32 + t34 * t55 + t21 * t13 + t22 * t12 + (t36 / 0.2e1 - t23 * mrSges(5,3) + t120 * t61) * t118 + (-t35 / 0.2e1 + t120 * t60 - t24 * mrSges(5,3)) * t116 + (-t82 + t43) * qJ(3) + m(6) * (t100 * t34 - t5 * t52 + t54 * t6) + m(4) * (-pkin(2) * t65 - qJ(3) * t62) + m(7) * (t1 * t21 + t2 * t22 + t3 * t52) + (t7 / 0.2e1 - t15 / 0.2e1 - t6 * mrSges(6,3)) * t80 + (-t56 / 0.2e1 + t29 / 0.2e1) * t41 + t31 * t165 + t30 * t166 + t132 - t155 + t156; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t100 * t55 - t116 * t86 + t118 * t87 + 0.2e1 * t21 * t49 + 0.2e1 * t22 * t48 + 0.2e1 * t52 * t45 + (t167 * t54 + t29 - t56) * t80 + (t121 * t30 - t124 * t31 + t167 * t52 - t57) * t79 + m(7) * (t21 ^ 2 + t22 ^ 2 + t168) + m(6) * (t100 ^ 2 + t54 ^ 2 + t168) + m(5) * (-t120 ^ 2 * t141 + t126) + m(4) * (pkin(2) ^ 2 + t126) + 0.2e1 * (mrSges(4,3) + t85) * qJ(3) + 0.2e1 * t120 * t134 + t172; t116 * t60 + t118 * t61 + t153 * t79 + (t32 + t170) * t80 + m(7) * (t131 * t80 + t3 * t79) + m(6) * (-t5 * t79 + t6 * t80) + m(5) * t129 + m(4) * t65 + t83; -m(4) * pkin(2) - t78 * mrSges(6,3) + t79 * t45 + mrSges(4,2) + t134 + (-mrSges(6,3) * t80 + t171) * t80 + m(7) * (t128 * t80 + t157) + m(6) * (t54 * t80 + t157) - t120 * t137; m(4) - t137 + m(6) * (t77 + t78) + m(7) * (t140 * t77 + t78); t121 * t12 + t124 * t13 + m(7) * (t1 * t124 + t121 * t2) + m(6) * t34 + m(5) * t59 + t20 + t43; m(5) * qJ(3) + t121 * t48 + t124 * t49 + m(7) * (t121 * t22 + t124 * t21) + m(6) * t100 + t55 + t85; 0; m(7) * t140 + m(5) + m(6); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t131 * mrSges(7,3) + t8 * t160 + t9 * t163 + t41 * t164 + t91 * t165 + t90 * t166 + t3 * t88 + t139 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t131 + t170) * pkin(10); t31 * t163 + t30 * t160 + t80 * t164 - pkin(5) * t45 - t54 * mrSges(6,2) + (t161 * t91 + t163 * t90) * t79 + t128 * mrSges(7,3) + t152 + t169 * t52 + (m(7) * t128 + t171) * pkin(10); t169 * t79 + (-mrSges(6,2) + (m(7) * pkin(10) + mrSges(7,3)) * t140) * t80; 0; Ifges(6,3) - 0.2e1 * pkin(5) * t88 + m(7) * (pkin(10) ^ 2 * t140 + pkin(5) ^ 2) + t121 * t91 + t124 * t90 + 0.2e1 * t140 * pkin(10) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t21 - mrSges(7,2) * t22 + t29; t130 * t80; -t88; pkin(10) * t130 + t89; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
