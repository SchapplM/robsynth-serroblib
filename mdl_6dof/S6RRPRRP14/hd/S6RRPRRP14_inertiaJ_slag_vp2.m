% Calculate joint inertia matrix for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:20:17
% EndTime: 2018-11-23 17:20:19
% DurationCPUTime: 1.56s
% Computational Cost: add. (1891->370), mult. (4121->496), div. (0->0), fcn. (4021->8), ass. (0->137)
t183 = Ifges(4,1) + Ifges(3,3);
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t158 = t125 ^ 2 + t128 ^ 2;
t123 = sin(pkin(6));
t127 = sin(qJ(2));
t165 = t123 * t127;
t124 = cos(pkin(6));
t126 = sin(qJ(4));
t129 = cos(qJ(4));
t130 = cos(qJ(2));
t164 = t123 * t130;
t63 = t124 * t129 - t126 * t164;
t37 = t125 * t63 - t128 * t165;
t38 = t125 * t165 + t128 * t63;
t62 = t124 * t126 + t129 * t164;
t7 = Ifges(6,5) * t38 - Ifges(6,6) * t37 + Ifges(6,3) * t62;
t8 = Ifges(7,4) * t38 + Ifges(7,2) * t62 + Ifges(7,6) * t37;
t182 = t7 + t8;
t161 = t128 * t129;
t163 = t125 * t129;
t52 = Ifges(6,5) * t161 - Ifges(6,6) * t163 + Ifges(6,3) * t126;
t53 = Ifges(7,4) * t161 + Ifges(7,2) * t126 + Ifges(7,6) * t163;
t181 = t52 + t53;
t180 = (mrSges(7,2) + mrSges(6,3)) * t158;
t131 = -pkin(2) - pkin(9);
t179 = pkin(1) * t130;
t94 = pkin(8) * t165;
t64 = t124 * t179 - t94;
t177 = t64 * mrSges(3,1);
t65 = t124 * t127 * pkin(1) + pkin(8) * t164;
t176 = t65 * mrSges(3,2);
t78 = -mrSges(6,1) * t128 + mrSges(6,2) * t125;
t175 = mrSges(5,1) - t78;
t148 = -pkin(2) - t179;
t29 = pkin(3) * t165 + t94 + (-pkin(9) + t148) * t124;
t147 = -qJ(3) * t127 - pkin(1);
t41 = (t130 * t131 + t147) * t123;
t16 = t126 * t29 + t129 * t41;
t14 = pkin(10) * t165 + t16;
t48 = -t124 * qJ(3) - t65;
t40 = pkin(3) * t164 - t48;
t18 = pkin(4) * t62 - pkin(10) * t63 + t40;
t4 = t125 * t18 + t128 * t14;
t20 = mrSges(6,1) * t37 + mrSges(6,2) * t38;
t44 = mrSges(5,1) * t165 - mrSges(5,3) * t63;
t174 = -t20 + t44;
t21 = -mrSges(6,2) * t62 - mrSges(6,3) * t37;
t24 = -mrSges(7,2) * t37 + mrSges(7,3) * t62;
t173 = t21 + t24;
t22 = mrSges(6,1) * t62 - mrSges(6,3) * t38;
t23 = -t62 * mrSges(7,1) + t38 * mrSges(7,2);
t172 = -t22 + t23;
t162 = t126 * t131;
t75 = pkin(4) * t126 - pkin(10) * t129 + qJ(3);
t46 = t125 * t75 + t128 * t162;
t171 = Ifges(6,4) * t125;
t170 = Ifges(6,4) * t128;
t169 = Ifges(7,5) * t125;
t168 = Ifges(7,5) * t128;
t167 = t128 * t75;
t15 = -t126 * t41 + t129 * t29;
t13 = -pkin(4) * t165 - t15;
t166 = t129 * t13;
t70 = mrSges(4,1) * t165 + t124 * mrSges(4,2);
t160 = t158 * pkin(10) * t126;
t81 = Ifges(6,5) * t125 + Ifges(6,6) * t128;
t159 = t158 * pkin(10) ^ 2;
t119 = t126 ^ 2;
t121 = t129 ^ 2;
t157 = -t121 - t119;
t6 = Ifges(7,5) * t38 + Ifges(7,6) * t62 + Ifges(7,3) * t37;
t9 = Ifges(6,4) * t38 - Ifges(6,2) * t37 + Ifges(6,6) * t62;
t156 = t6 / 0.2e1 - t9 / 0.2e1;
t155 = Ifges(5,5) * t63 - Ifges(5,6) * t62 + Ifges(5,3) * t165;
t10 = Ifges(7,1) * t38 + Ifges(7,4) * t62 + Ifges(7,5) * t37;
t11 = Ifges(6,1) * t38 - Ifges(6,4) * t37 + Ifges(6,5) * t62;
t154 = t10 / 0.2e1 + t11 / 0.2e1;
t51 = Ifges(7,6) * t126 + (Ifges(7,3) * t125 + t168) * t129;
t54 = Ifges(6,6) * t126 + (-Ifges(6,2) * t125 + t170) * t129;
t153 = t51 / 0.2e1 - t54 / 0.2e1;
t55 = Ifges(7,4) * t126 + (Ifges(7,1) * t128 + t169) * t129;
t56 = Ifges(6,5) * t126 + (Ifges(6,1) * t128 - t171) * t129;
t152 = t55 / 0.2e1 + t56 / 0.2e1;
t80 = -Ifges(7,3) * t128 + t169;
t83 = Ifges(6,2) * t128 + t171;
t151 = t80 / 0.2e1 - t83 / 0.2e1;
t82 = Ifges(7,4) * t125 - Ifges(7,6) * t128;
t150 = t81 / 0.2e1 + t82 / 0.2e1;
t85 = Ifges(7,1) * t125 - t168;
t86 = Ifges(6,1) * t125 + t170;
t149 = t85 / 0.2e1 + t86 / 0.2e1;
t73 = -t126 * mrSges(7,1) + mrSges(7,2) * t161;
t146 = t157 * mrSges(5,3);
t145 = Ifges(3,5) * t165 + Ifges(3,6) * t164 + t183 * t124;
t1 = qJ(6) * t62 + t4;
t3 = -t125 * t14 + t128 * t18;
t2 = -pkin(5) * t62 - t3;
t143 = t1 * t128 + t2 * t125;
t142 = -t3 * t125 + t4 * t128;
t141 = t125 * mrSges(6,1) + t128 * mrSges(6,2);
t140 = t125 * mrSges(7,1) - t128 * mrSges(7,3);
t139 = -pkin(5) * t125 + qJ(6) * t128;
t39 = qJ(6) * t126 + t46;
t42 = -t167 + (t125 * t131 - pkin(5)) * t126;
t138 = t125 * t42 + t128 * t39;
t45 = -t125 * t162 + t167;
t137 = -t125 * t45 + t128 * t46;
t136 = t16 * t126 + t15 * t129;
t71 = -t126 * mrSges(6,2) - mrSges(6,3) * t163;
t72 = t126 * mrSges(6,1) - mrSges(6,3) * t161;
t74 = -mrSges(7,2) * t163 + t126 * mrSges(7,3);
t135 = (t71 + t74) * t128 + (-t72 + t73) * t125;
t134 = m(7) * t139 - t140 - t141;
t132 = qJ(3) ^ 2;
t122 = t131 ^ 2;
t114 = Ifges(5,5) * t129;
t105 = t121 * t131;
t104 = t121 * t122;
t87 = Ifges(5,1) * t129 - Ifges(5,4) * t126;
t84 = Ifges(5,4) * t129 - Ifges(5,2) * t126;
t79 = mrSges(5,1) * t126 + mrSges(5,2) * t129;
t77 = -mrSges(7,1) * t128 - mrSges(7,3) * t125;
t76 = -pkin(5) * t128 - qJ(6) * t125 - pkin(4);
t69 = -mrSges(4,1) * t164 - mrSges(4,3) * t124;
t67 = t141 * t129;
t66 = t140 * t129;
t50 = t124 * t148 + t94;
t49 = (-pkin(2) * t130 + t147) * t123;
t47 = (-t131 - t139) * t129;
t43 = -mrSges(5,2) * t165 - mrSges(5,3) * t62;
t27 = mrSges(5,1) * t62 + mrSges(5,2) * t63;
t26 = Ifges(5,1) * t63 - Ifges(5,4) * t62 + Ifges(5,5) * t165;
t25 = Ifges(5,4) * t63 - Ifges(5,2) * t62 + Ifges(5,6) * t165;
t19 = mrSges(7,1) * t37 - mrSges(7,3) * t38;
t5 = pkin(5) * t37 - qJ(6) * t38 + t13;
t12 = [0.2e1 * t1 * t24 + 0.2e1 * t13 * t20 + 0.2e1 * t15 * t44 + 0.2e1 * t16 * t43 + 0.2e1 * t5 * t19 + 0.2e1 * t2 * t23 + 0.2e1 * t4 * t21 + 0.2e1 * t3 * t22 + t63 * t26 + 0.2e1 * t40 * t27 + 0.2e1 * t48 * t69 + 0.2e1 * t50 * t70 + Ifges(2,3) + (t10 + t11) * t38 + (t6 - t9) * t37 + (-t25 + t182) * t62 + (t145 - 0.2e1 * t176 + 0.2e1 * t177) * t124 + m(4) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t40 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(3) * (t64 ^ 2 + t65 ^ 2) + (0.2e1 * t49 * (mrSges(4,2) * t130 - mrSges(4,3) * t127) + t127 * t155 + (m(3) * pkin(1) ^ 2 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t130) * t130 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1)) * t127 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t130) * t127) * t123 + 0.2e1 * (-t64 * t127 + t65 * t130) * mrSges(3,3) + ((-(2 * Ifges(4,5)) + Ifges(3,6)) * t130 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t127) * t124) * t123; (t26 / 0.2e1 - t15 * mrSges(5,3) + t174 * t131 + t154 * t128 + t156 * t125) * t129 + (-t25 / 0.2e1 + t7 / 0.2e1 + t8 / 0.2e1 + t131 * t43 - t16 * mrSges(5,3)) * t126 + m(6) * (-t131 * t166 + t3 * t45 + t4 * t46) + t152 * t38 + t153 * t37 + (-t69 + t27) * qJ(3) + m(4) * (-pkin(2) * t50 - qJ(3) * t48) + m(7) * (t1 * t39 + t2 * t42 + t47 * t5) + (-t84 / 0.2e1 + t52 / 0.2e1 + t53 / 0.2e1) * t62 + m(5) * (qJ(3) * t40 + t131 * t136) + t177 - t176 + t2 * t73 + t1 * t74 + t40 * t79 + t63 * t87 / 0.2e1 + t5 * t66 + t13 * t67 - pkin(2) * t70 + t4 * t71 + t3 * t72 + t39 * t24 + t42 * t23 + t45 * t22 + t46 * t21 + t47 * t19 - t48 * mrSges(4,3) + t50 * mrSges(4,2) + (-Ifges(4,5) * t130 + t127 * (-Ifges(5,6) * t126 + t114) / 0.2e1 - Ifges(4,4) * t127) * t123 + t145; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t39 * t74 + 0.2e1 * t42 * t73 + 0.2e1 * t45 * t72 + 0.2e1 * t46 * t71 + 0.2e1 * t47 * t66 + (-t84 + t181) * t126 + m(5) * (t119 * t122 + t104 + t132) + m(4) * (pkin(2) ^ 2 + t132) + m(6) * (t45 ^ 2 + t46 ^ 2 + t104) + m(7) * (t39 ^ 2 + t42 ^ 2 + t47 ^ 2) + (-0.2e1 * t131 * t67 + t87 + (t55 + t56) * t128 + (t51 - t54) * t125) * t129 + 0.2e1 * (t79 + mrSges(4,3)) * qJ(3) + 0.2e1 * t131 * t146 + t183; (-t19 + t174) * t129 + (t125 * t172 + t128 * t173 + t43) * t126 + m(7) * (t126 * t143 - t129 * t5) + m(6) * (t126 * t142 - t166) + m(5) * t136 + m(4) * t50 + t70; -m(4) * pkin(2) + mrSges(4,2) + (-t66 - t67) * t129 + t146 + t135 * t126 + m(6) * (t126 * t137 + t105) + m(7) * (t126 * t138 - t129 * t47) + m(5) * (t119 * t131 + t105); m(4) - m(5) * t157 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t119 * t158 + t121); t15 * mrSges(5,1) - t16 * mrSges(5,2) - pkin(4) * t20 + t13 * t78 + t76 * t19 + t5 * t77 + t150 * t62 + t149 * t38 + t151 * t37 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) + pkin(10) * t173 - t156) * t128 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + pkin(10) * t172 + t154) * t125 + m(6) * (-pkin(4) * t13 + pkin(10) * t142) + m(7) * (pkin(10) * t143 + t5 * t76) + t155; -pkin(4) * t67 + t76 * t66 + t114 + (m(7) * t76 + t77) * t47 + (t39 * mrSges(7,2) + t46 * mrSges(6,3) - t153) * t128 + (t42 * mrSges(7,2) - t45 * mrSges(6,3) + t152) * t125 + (m(6) * t137 + m(7) * t138 + t135) * pkin(10) + (t149 * t128 + t151 * t125 + (m(6) * pkin(4) + t175) * t131) * t129 + (-t131 * mrSges(5,2) - Ifges(5,6) + t150) * t126; (-t77 + t175) * t129 + m(6) * (pkin(4) * t129 + t160) + m(7) * (-t129 * t76 + t160) + (-mrSges(5,2) + t180) * t126; -0.2e1 * pkin(4) * t78 + 0.2e1 * t76 * t77 + Ifges(5,3) + (-t80 + t83) * t128 + (t85 + t86) * t125 + m(7) * (t76 ^ 2 + t159) + m(6) * (pkin(4) ^ 2 + t159) + 0.2e1 * pkin(10) * t180; -pkin(5) * t23 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t24 + t1 * mrSges(7,3) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t182; -pkin(5) * t73 + m(7) * (-pkin(5) * t42 + qJ(6) * t39) + qJ(6) * t74 + t39 * mrSges(7,3) + t45 * mrSges(6,1) - t42 * mrSges(7,1) - t46 * mrSges(6,2) + t181; t134 * t126; mrSges(7,2) * t139 + pkin(10) * t134 + t81 + t82; Ifges(7,2) + Ifges(6,3) + 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2); m(7) * t2 + t23; m(7) * t42 + t73; m(7) * t125 * t126; (m(7) * pkin(10) + mrSges(7,2)) * t125; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
