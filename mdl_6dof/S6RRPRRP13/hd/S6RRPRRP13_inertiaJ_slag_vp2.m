% Calculate joint inertia matrix for
% S6RRPRRP13
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

function Mq = S6RRPRRP13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:19:48
% EndTime: 2018-11-23 17:19:49
% DurationCPUTime: 1.60s
% Computational Cost: add. (1872->385), mult. (4080->508), div. (0->0), fcn. (4014->8), ass. (0->132)
t174 = Ifges(4,1) + Ifges(3,3);
t173 = -2 * mrSges(7,3);
t126 = sin(qJ(5));
t129 = cos(qJ(5));
t124 = sin(pkin(6));
t128 = sin(qJ(2));
t159 = t124 * t128;
t125 = cos(pkin(6));
t127 = sin(qJ(4));
t130 = cos(qJ(4));
t131 = cos(qJ(2));
t158 = t124 * t131;
t63 = t125 * t130 - t127 * t158;
t40 = -t126 * t63 + t129 * t159;
t41 = t126 * t159 + t129 * t63;
t62 = t125 * t127 + t130 * t158;
t6 = Ifges(7,5) * t41 + Ifges(7,6) * t40 + Ifges(7,3) * t62;
t7 = Ifges(6,5) * t41 + Ifges(6,6) * t40 + Ifges(6,3) * t62;
t172 = t6 + t7;
t171 = m(6) * pkin(4);
t170 = m(7) * pkin(5);
t132 = -pkin(2) - pkin(9);
t169 = pkin(1) * t131;
t96 = pkin(8) * t159;
t64 = t125 * t169 - t96;
t168 = t64 * mrSges(3,1);
t65 = t125 * t128 * pkin(1) + pkin(8) * t158;
t167 = t65 * mrSges(3,2);
t80 = -mrSges(6,1) * t129 + mrSges(6,2) * t126;
t166 = mrSges(5,1) - t80;
t165 = -qJ(6) - pkin(10);
t141 = -pkin(2) - t169;
t31 = pkin(3) * t159 + t96 + (-pkin(9) + t141) * t125;
t140 = -qJ(3) * t128 - pkin(1);
t43 = (t132 * t131 + t140) * t124;
t16 = t127 * t31 + t130 * t43;
t14 = pkin(10) * t159 + t16;
t48 = -t125 * qJ(3) - t65;
t42 = pkin(3) * t158 - t48;
t19 = pkin(4) * t62 - pkin(10) * t63 + t42;
t4 = t126 * t19 + t129 * t14;
t21 = -mrSges(6,1) * t40 + mrSges(6,2) * t41;
t45 = mrSges(5,1) * t159 - mrSges(5,3) * t63;
t164 = t45 - t21;
t156 = t127 * t132;
t77 = pkin(4) * t127 - pkin(10) * t130 + qJ(3);
t47 = t126 * t77 + t129 * t156;
t163 = Ifges(6,4) * t126;
t162 = Ifges(6,4) * t129;
t161 = Ifges(7,4) * t126;
t160 = Ifges(7,4) * t129;
t72 = mrSges(4,1) * t159 + t125 * mrSges(4,2);
t157 = t126 * t130;
t155 = t129 * t130;
t154 = t130 * t132;
t66 = mrSges(7,1) * t157 + mrSges(7,2) * t155;
t153 = Ifges(7,5) * t155 + Ifges(7,3) * t127;
t152 = Ifges(6,5) * t155 + Ifges(6,3) * t127;
t83 = Ifges(7,5) * t126 + Ifges(7,6) * t129;
t84 = Ifges(6,5) * t126 + Ifges(6,6) * t129;
t151 = t126 ^ 2 + t129 ^ 2;
t120 = t127 ^ 2;
t122 = t130 ^ 2;
t150 = -t122 - t120;
t8 = Ifges(7,4) * t41 + Ifges(7,2) * t40 + Ifges(7,6) * t62;
t9 = Ifges(6,4) * t41 + Ifges(6,2) * t40 + Ifges(6,6) * t62;
t149 = -t8 / 0.2e1 - t9 / 0.2e1;
t148 = Ifges(5,5) * t63 - Ifges(5,6) * t62 + Ifges(5,3) * t159;
t10 = Ifges(7,1) * t41 + Ifges(7,4) * t40 + Ifges(7,5) * t62;
t11 = Ifges(6,1) * t41 + Ifges(6,4) * t40 + Ifges(6,5) * t62;
t147 = t10 / 0.2e1 + t11 / 0.2e1;
t53 = Ifges(7,6) * t127 + (-Ifges(7,2) * t126 + t160) * t130;
t54 = Ifges(6,6) * t127 + (-Ifges(6,2) * t126 + t162) * t130;
t146 = t53 / 0.2e1 + t54 / 0.2e1;
t55 = Ifges(7,5) * t127 + (Ifges(7,1) * t129 - t161) * t130;
t56 = Ifges(6,5) * t127 + (Ifges(6,1) * t129 - t163) * t130;
t145 = t55 / 0.2e1 + t56 / 0.2e1;
t144 = t83 / 0.2e1 + t84 / 0.2e1;
t85 = Ifges(7,2) * t129 + t161;
t86 = Ifges(6,2) * t129 + t163;
t143 = t86 / 0.2e1 + t85 / 0.2e1;
t88 = Ifges(7,1) * t126 + t160;
t89 = Ifges(6,1) * t126 + t162;
t142 = t88 / 0.2e1 + t89 / 0.2e1;
t20 = -t40 * mrSges(7,1) + t41 * mrSges(7,2);
t3 = -t126 * t14 + t129 * t19;
t15 = -t127 * t43 + t130 * t31;
t139 = t151 * mrSges(6,3);
t138 = t150 * mrSges(5,3);
t79 = -t129 * mrSges(7,1) + t126 * mrSges(7,2);
t137 = Ifges(3,5) * t159 + Ifges(3,6) * t158 + t174 * t125;
t136 = mrSges(6,1) * t126 + mrSges(6,2) * t129;
t135 = t16 * t127 + t15 * t130;
t13 = -pkin(4) * t159 - t15;
t133 = qJ(3) ^ 2;
t123 = t132 ^ 2;
t118 = Ifges(5,5) * t130;
t106 = t122 * t132;
t105 = t122 * t123;
t104 = -pkin(5) * t129 - pkin(4);
t90 = Ifges(5,1) * t130 - Ifges(5,4) * t127;
t87 = Ifges(5,4) * t130 - Ifges(5,2) * t127;
t82 = mrSges(5,1) * t127 + mrSges(5,2) * t130;
t81 = t165 * t129;
t78 = t165 * t126;
t76 = mrSges(6,1) * t127 - mrSges(6,3) * t155;
t75 = mrSges(7,1) * t127 - mrSges(7,3) * t155;
t74 = -mrSges(6,2) * t127 - mrSges(6,3) * t157;
t73 = -mrSges(7,2) * t127 - mrSges(7,3) * t157;
t71 = -mrSges(4,1) * t158 - mrSges(4,3) * t125;
t70 = (pkin(5) * t126 - t132) * t130;
t69 = t129 * t77;
t67 = t136 * t130;
t52 = -Ifges(6,6) * t157 + t152;
t51 = -Ifges(7,6) * t157 + t153;
t50 = t125 * t141 + t96;
t49 = (-pkin(2) * t131 + t140) * t124;
t46 = -t126 * t156 + t69;
t44 = -mrSges(5,2) * t159 - mrSges(5,3) * t62;
t39 = -qJ(6) * t157 + t47;
t30 = -qJ(6) * t155 + t69 + (-t126 * t132 + pkin(5)) * t127;
t28 = mrSges(5,1) * t62 + mrSges(5,2) * t63;
t27 = Ifges(5,1) * t63 - Ifges(5,4) * t62 + Ifges(5,5) * t159;
t26 = Ifges(5,4) * t63 - Ifges(5,2) * t62 + Ifges(5,6) * t159;
t25 = mrSges(6,1) * t62 - mrSges(6,3) * t41;
t24 = mrSges(7,1) * t62 - mrSges(7,3) * t41;
t23 = -mrSges(6,2) * t62 + mrSges(6,3) * t40;
t22 = -mrSges(7,2) * t62 + mrSges(7,3) * t40;
t5 = -pkin(5) * t40 + t13;
t2 = qJ(6) * t40 + t4;
t1 = pkin(5) * t62 - qJ(6) * t41 + t3;
t12 = [0.2e1 * t1 * t24 + 0.2e1 * t13 * t21 + 0.2e1 * t15 * t45 + 0.2e1 * t16 * t44 + 0.2e1 * t2 * t22 + 0.2e1 * t5 * t20 + 0.2e1 * t4 * t23 + 0.2e1 * t3 * t25 + t63 * t27 + 0.2e1 * t42 * t28 + 0.2e1 * t48 * t71 + 0.2e1 * t50 * t72 + Ifges(2,3) + (t10 + t11) * t41 + (t8 + t9) * t40 + (-t26 + t172) * t62 + (t137 - 0.2e1 * t167 + 0.2e1 * t168) * t125 + m(4) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t42 ^ 2) + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (t64 ^ 2 + t65 ^ 2) + (t128 * t148 + 0.2e1 * t49 * (mrSges(4,2) * t131 - mrSges(4,3) * t128) + (m(3) * pkin(1) ^ 2 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t131) * t131 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(3,1) + Ifges(4,2)) * t128 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t131) * t128) * t124 + 0.2e1 * (-t64 * t128 + t65 * t131) * mrSges(3,3) + ((-(2 * Ifges(4,5)) + Ifges(3,6)) * t131 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t128) * t125) * t124; m(6) * (-t13 * t154 + t3 * t46 + t4 * t47) + t145 * t41 + t146 * t40 + m(5) * (qJ(3) * t42 + t132 * t135) + (t132 * t44 - t16 * mrSges(5,3) - t26 / 0.2e1 + t6 / 0.2e1 + t7 / 0.2e1) * t127 + (-t71 + t28) * qJ(3) + m(7) * (t1 * t30 + t2 * t39 + t5 * t70) + m(4) * (-pkin(2) * t50 - qJ(3) * t48) + (-t87 / 0.2e1 + t51 / 0.2e1 + t52 / 0.2e1) * t62 + t63 * t90 / 0.2e1 + t70 * t20 - pkin(2) * t72 + t2 * t73 + t4 * t74 + t1 * t75 + t3 * t76 + t42 * t82 + t5 * t66 + t13 * t67 - t48 * mrSges(4,3) + t50 * mrSges(4,2) + t46 * t25 + t47 * t23 + t30 * t24 + t39 * t22 + t168 - t167 + (-Ifges(4,5) * t131 + t128 * (-Ifges(5,6) * t127 + t118) / 0.2e1 - Ifges(4,4) * t128) * t124 + t137 + (-t15 * mrSges(5,3) + t27 / 0.2e1 + t164 * t132 + t147 * t129 + t149 * t126) * t130; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t30 * t75 + 0.2e1 * t39 * t73 + 0.2e1 * t46 * t76 + 0.2e1 * t47 * t74 + 0.2e1 * t70 * t66 + (-t87 + t51 + t52) * t127 + m(5) * (t120 * t123 + t105 + t133) + m(4) * (pkin(2) ^ 2 + t133) + m(6) * (t46 ^ 2 + t47 ^ 2 + t105) + m(7) * (t30 ^ 2 + t39 ^ 2 + t70 ^ 2) + (-0.2e1 * t132 * t67 + t90 + (t55 + t56) * t129 + (-t53 - t54) * t126) * t130 + 0.2e1 * (t82 + mrSges(4,3)) * qJ(3) + 0.2e1 * t132 * t138 + t174; (-t20 + t164) * t130 + (t44 + (t22 + t23) * t129 + (-t24 - t25) * t126) * t127 + m(7) * (-t130 * t5 + (-t1 * t126 + t129 * t2) * t127) + m(6) * (-t13 * t130 + (-t126 * t3 + t129 * t4) * t127) + m(5) * t135 + m(4) * t50 + t72; -m(4) * pkin(2) + mrSges(4,2) + (-t66 - t67) * t130 + t138 + ((t73 + t74) * t129 + (-t75 - t76) * t126) * t127 + m(6) * (t106 + (-t126 * t46 + t129 * t47) * t127) + m(7) * (-t130 * t70 + (-t126 * t30 + t129 * t39) * t127) + m(5) * (t120 * t132 + t106); m(4) - m(5) * t150 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t120 * t151 + t122); t15 * mrSges(5,1) - t16 * mrSges(5,2) - pkin(4) * t21 + t104 * t20 - t81 * t22 + t78 * t24 + t5 * t79 + t144 * t62 + t142 * t41 + t143 * t40 + m(7) * (t1 * t78 + t104 * t5 - t2 * t81) + (t4 * mrSges(6,3) + t2 * mrSges(7,3) + (m(6) * t4 + t23) * pkin(10) - t149) * t129 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + (-m(6) * t3 - t25) * pkin(10) + t147) * t126 + t148 + (t80 - t171) * t13; m(7) * (t104 * t70 + t30 * t78 - t39 * t81) + t118 + t104 * t66 + t78 * t75 + t70 * t79 - t81 * t73 - pkin(4) * t67 + (t166 + t171) * t154 + (-t132 * mrSges(5,2) - Ifges(5,6) + t144) * t127 + (t47 * mrSges(6,3) + t39 * mrSges(7,3) + t142 * t130 + (m(6) * t47 + t74) * pkin(10) + t146) * t129 + (-t46 * mrSges(6,3) - t30 * mrSges(7,3) - t143 * t130 + (-m(6) * t46 - t76) * pkin(10) + t145) * t126; (-t79 + t166) * t130 + (mrSges(7,3) * t151 - mrSges(5,2) + t139) * t127 + m(6) * (pkin(10) * t127 * t151 + pkin(4) * t130) + m(7) * (-t104 * t130 + (-t126 * t78 - t129 * t81) * t127); -0.2e1 * pkin(4) * t80 + 0.2e1 * t104 * t79 + Ifges(5,3) + 0.2e1 * pkin(10) * t139 + m(7) * (t104 ^ 2 + t78 ^ 2 + t81 ^ 2) + m(6) * (pkin(10) ^ 2 * t151 + pkin(4) ^ 2) + (t81 * t173 + t85 + t86) * t129 + (t78 * t173 + t88 + t89) * t126; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t24) * pkin(5) + t172; mrSges(6,1) * t46 + mrSges(7,1) * t30 - mrSges(6,2) * t47 - mrSges(7,2) * t39 + (-Ifges(6,6) - Ifges(7,6)) * t157 + (m(7) * t30 + t75) * pkin(5) + t152 + t153; ((-mrSges(6,2) - mrSges(7,2)) * t129 + (-mrSges(6,1) - mrSges(7,1) - t170) * t126) * t127; mrSges(7,1) * t78 + mrSges(7,2) * t81 - t136 * pkin(10) + (m(7) * t78 - t126 * mrSges(7,3)) * pkin(5) + t84 + t83; Ifges(6,3) + Ifges(7,3) + (0.2e1 * mrSges(7,1) + t170) * pkin(5); m(7) * t5 + t20; m(7) * t70 + t66; -m(7) * t130; m(7) * t104 + t79; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
