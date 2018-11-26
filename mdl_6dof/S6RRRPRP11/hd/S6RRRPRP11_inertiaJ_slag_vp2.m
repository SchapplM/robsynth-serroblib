% Calculate joint inertia matrix for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:48:53
% EndTime: 2018-11-23 17:48:55
% DurationCPUTime: 1.93s
% Computational Cost: add. (1987->400), mult. (4331->523), div. (0->0), fcn. (4347->8), ass. (0->135)
t182 = pkin(4) + pkin(9);
t181 = Ifges(6,6) + Ifges(7,6);
t133 = sin(qJ(3));
t136 = cos(qJ(3));
t180 = t133 ^ 2 + t136 ^ 2;
t132 = sin(qJ(5));
t135 = cos(qJ(5));
t130 = sin(pkin(6));
t137 = cos(qJ(2));
t162 = t130 * t137;
t131 = cos(pkin(6));
t134 = sin(qJ(2));
t163 = t130 * t134;
t70 = -t131 * t136 + t133 * t163;
t44 = t132 * t162 + t135 * t70;
t45 = t132 * t70 - t135 * t162;
t71 = t131 * t133 + t136 * t163;
t8 = Ifges(7,5) * t45 + Ifges(7,6) * t44 + Ifges(7,3) * t71;
t9 = Ifges(6,5) * t45 + Ifges(6,6) * t44 + Ifges(6,3) * t71;
t179 = t8 + t9;
t177 = -m(5) * pkin(3) + mrSges(5,2);
t176 = -2 * mrSges(7,3);
t175 = m(7) * pkin(5);
t174 = pkin(3) + pkin(10);
t73 = t131 * t134 * pkin(1) + pkin(8) * t162;
t54 = pkin(9) * t131 + t73;
t55 = (-pkin(2) * t137 - pkin(9) * t134 - pkin(1)) * t130;
t27 = -t133 * t54 + t136 * t55;
t22 = pkin(3) * t162 - t27;
t14 = pkin(4) * t71 + pkin(10) * t162 + t22;
t105 = pkin(8) * t163;
t173 = pkin(1) * t137;
t53 = t105 + (-pkin(2) - t173) * t131;
t141 = -qJ(4) * t71 + t53;
t16 = t174 * t70 + t141;
t4 = t132 * t14 + t135 * t16;
t72 = t131 * t173 - t105;
t172 = t72 * mrSges(3,1);
t171 = t73 * mrSges(3,2);
t170 = Ifges(5,1) + Ifges(4,3);
t28 = t133 * t55 + t136 * t54;
t100 = t182 * t133;
t146 = -qJ(4) * t133 - pkin(2);
t76 = -t174 * t136 + t146;
t38 = t132 * t100 + t135 * t76;
t168 = mrSges(6,1) * t135;
t167 = Ifges(6,4) * t132;
t166 = Ifges(6,4) * t135;
t165 = Ifges(7,4) * t132;
t164 = Ifges(7,4) * t135;
t161 = t132 * t136;
t160 = t135 * t136;
t159 = -qJ(6) - t174;
t88 = t132 * mrSges(7,1) + t135 * mrSges(7,2);
t158 = Ifges(4,5) * t133 + Ifges(4,6) * t136;
t157 = t180 * pkin(9) ^ 2;
t101 = t182 * t136;
t156 = t132 ^ 2 + t135 ^ 2;
t10 = Ifges(7,4) * t45 + Ifges(7,2) * t44 + Ifges(7,6) * t71;
t11 = Ifges(6,4) * t45 + Ifges(6,2) * t44 + Ifges(6,6) * t71;
t155 = -t10 / 0.2e1 - t11 / 0.2e1;
t12 = Ifges(7,1) * t45 + Ifges(7,4) * t44 + Ifges(7,5) * t71;
t13 = Ifges(6,1) * t45 + Ifges(6,4) * t44 + Ifges(6,5) * t71;
t154 = -t12 / 0.2e1 - t13 / 0.2e1;
t58 = Ifges(7,6) * t133 + (-Ifges(7,2) * t135 - t165) * t136;
t59 = Ifges(6,6) * t133 + (-Ifges(6,2) * t135 - t167) * t136;
t153 = t58 / 0.2e1 + t59 / 0.2e1;
t60 = Ifges(7,5) * t133 + (-Ifges(7,1) * t132 - t164) * t136;
t61 = Ifges(6,5) * t133 + (-Ifges(6,1) * t132 - t166) * t136;
t152 = t60 / 0.2e1 + t61 / 0.2e1;
t118 = Ifges(7,5) * t135;
t119 = Ifges(6,5) * t135;
t151 = t118 / 0.2e1 + t119 / 0.2e1 - t181 * t132 / 0.2e1;
t94 = -Ifges(7,2) * t132 + t164;
t95 = -Ifges(6,2) * t132 + t166;
t150 = t94 / 0.2e1 + t95 / 0.2e1;
t97 = Ifges(7,1) * t135 - t165;
t98 = Ifges(6,1) * t135 - t167;
t149 = t97 / 0.2e1 + t98 / 0.2e1;
t148 = Ifges(3,5) * t163 + Ifges(3,6) * t162 + Ifges(3,3) * t131;
t147 = t162 / 0.2e1;
t18 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t3 = -t132 * t16 + t135 * t14;
t145 = t156 * mrSges(6,3);
t144 = (Ifges(5,4) - Ifges(4,5)) * t71 + (-Ifges(5,5) + Ifges(4,6)) * t70;
t74 = mrSges(7,1) * t160 - mrSges(7,2) * t161;
t142 = t132 * t4 + t135 * t3;
t21 = qJ(4) * t162 - t28;
t47 = t71 * mrSges(5,1) - mrSges(5,2) * t162;
t17 = -pkin(4) * t70 - t21;
t139 = qJ(4) ^ 2;
t116 = Ifges(6,3) * t133;
t115 = Ifges(7,3) * t133;
t111 = pkin(5) * t132 + qJ(4);
t99 = Ifges(4,1) * t133 + Ifges(4,4) * t136;
t96 = Ifges(4,4) * t133 + Ifges(4,2) * t136;
t91 = -Ifges(5,2) * t133 - Ifges(5,6) * t136;
t90 = -Ifges(5,6) * t133 - Ifges(5,3) * t136;
t89 = mrSges(6,1) * t132 + mrSges(6,2) * t135;
t87 = -mrSges(4,1) * t136 + mrSges(4,2) * t133;
t86 = mrSges(5,2) * t136 - mrSges(5,3) * t133;
t85 = -pkin(3) * t136 + t146;
t84 = t159 * t135;
t83 = t159 * t132;
t82 = -mrSges(6,2) * t133 - mrSges(6,3) * t160;
t81 = -mrSges(7,2) * t133 - mrSges(7,3) * t160;
t80 = mrSges(6,1) * t133 + mrSges(6,3) * t161;
t79 = mrSges(7,1) * t133 + mrSges(7,3) * t161;
t78 = t135 * t100;
t75 = (-mrSges(6,2) * t132 + t168) * t136;
t69 = pkin(5) * t160 + t101;
t57 = t116 + (-Ifges(6,5) * t132 - Ifges(6,6) * t135) * t136;
t56 = t115 + (-Ifges(7,5) * t132 - Ifges(7,6) * t135) * t136;
t49 = -mrSges(4,1) * t162 - mrSges(4,3) * t71;
t48 = mrSges(4,2) * t162 - mrSges(4,3) * t70;
t46 = mrSges(5,1) * t70 + mrSges(5,3) * t162;
t37 = -t132 * t76 + t78;
t36 = -mrSges(5,2) * t70 - mrSges(5,3) * t71;
t35 = mrSges(4,1) * t70 + mrSges(4,2) * t71;
t34 = -qJ(6) * t160 + t38;
t33 = pkin(5) * t133 + t78 + (qJ(6) * t136 - t76) * t132;
t32 = Ifges(4,1) * t71 - Ifges(4,4) * t70 - Ifges(4,5) * t162;
t31 = Ifges(4,4) * t71 - Ifges(4,2) * t70 - Ifges(4,6) * t162;
t30 = -Ifges(5,4) * t162 - Ifges(5,2) * t71 + Ifges(5,6) * t70;
t29 = -Ifges(5,5) * t162 - Ifges(5,6) * t71 + Ifges(5,3) * t70;
t26 = mrSges(6,1) * t71 - mrSges(6,3) * t45;
t25 = mrSges(7,1) * t71 - mrSges(7,3) * t45;
t24 = -mrSges(6,2) * t71 + mrSges(6,3) * t44;
t23 = -mrSges(7,2) * t71 + mrSges(7,3) * t44;
t20 = pkin(3) * t70 + t141;
t19 = -mrSges(6,1) * t44 + mrSges(6,2) * t45;
t5 = -pkin(5) * t44 + t17;
t2 = qJ(6) * t44 + t4;
t1 = pkin(5) * t71 - qJ(6) * t45 + t3;
t6 = [0.2e1 * t1 * t25 + 0.2e1 * t17 * t19 + 0.2e1 * t5 * t18 + 0.2e1 * t2 * t23 + 0.2e1 * t20 * t36 + 0.2e1 * t21 * t46 + 0.2e1 * t22 * t47 + 0.2e1 * t4 * t24 + 0.2e1 * t3 * t26 + 0.2e1 * t27 * t49 + 0.2e1 * t28 * t48 + 0.2e1 * t53 * t35 + Ifges(2,3) + (t29 - t31) * t70 + (t12 + t13) * t45 + (t10 + t11) * t44 + (t148 - 0.2e1 * t171 + 0.2e1 * t172) * t131 + (-t30 + t32 + t179) * t71 + ((-0.2e1 * t72 * mrSges(3,3) + Ifges(3,5) * t131 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t134) * t130) * t134 + (0.2e1 * t73 * mrSges(3,3) + Ifges(3,6) * t131 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t134 + (Ifges(3,2) + t170) * t137) * t130 + t144) * t137) * t130 + m(3) * (pkin(1) ^ 2 * t130 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(4) * (t27 ^ 2 + t28 ^ 2 + t53 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t17 ^ 2 + t3 ^ 2 + t4 ^ 2); -t158 * t162 / 0.2e1 + (-t91 / 0.2e1 + t99 / 0.2e1 + t56 / 0.2e1 + t57 / 0.2e1) * t71 + (t90 / 0.2e1 - t96 / 0.2e1) * t70 + t148 + t85 * t36 + t20 * t86 + t53 * t87 + t101 * t19 + t5 * t74 + t17 * t75 + t1 * t79 + t3 * t80 + t2 * t81 + t4 * t82 + m(5) * (t20 * t85 + (t133 * t22 - t136 * t21) * pkin(9)) + m(4) * (-pkin(2) * t53 + (-t133 * t27 + t136 * t28) * pkin(9)) + t69 * t18 + t34 * t23 - pkin(2) * t35 + t37 * t26 + t38 * t24 + t33 * t25 + t152 * t45 + t153 * t44 + (Ifges(5,5) * t147 - t21 * mrSges(5,1) + t28 * mrSges(4,3) - t29 / 0.2e1 + t31 / 0.2e1 + t155 * t135 + t154 * t132 + (t48 - t46) * pkin(9)) * t136 + (Ifges(5,4) * t147 + t22 * mrSges(5,1) - t27 * mrSges(4,3) - t30 / 0.2e1 + t32 / 0.2e1 + t8 / 0.2e1 + t9 / 0.2e1 + (-t49 + t47) * pkin(9)) * t133 - t171 + t172 + m(6) * (t101 * t17 + t3 * t37 + t38 * t4) + m(7) * (t1 * t33 + t2 * t34 + t5 * t69); -0.2e1 * pkin(2) * t87 + 0.2e1 * t101 * t75 + 0.2e1 * t33 * t79 + 0.2e1 * t34 * t81 + 0.2e1 * t37 * t80 + 0.2e1 * t38 * t82 + 0.2e1 * t69 * t74 + 0.2e1 * t85 * t86 + Ifges(3,3) + (-t91 + t99 + t56 + t57) * t133 + m(5) * (t85 ^ 2 + t157) + m(4) * (pkin(2) ^ 2 + t157) + m(6) * (t101 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(7) * (t33 ^ 2 + t34 ^ 2 + t69 ^ 2) + (-t90 + t96 + (-t58 - t59) * t135 + (-t60 - t61) * t132) * t136 + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * pkin(9) * t180; -t144 - t170 * t162 + t111 * t18 + t83 * t23 + t84 * t25 + t5 * t88 + t17 * t89 + (-t3 * mrSges(6,3) - t1 * mrSges(7,3) - t174 * t26 - t154) * t135 + (-t4 * mrSges(6,3) - t2 * mrSges(7,3) - t174 * t24 + t155) * t132 - pkin(3) * t47 + t27 * mrSges(4,1) - t28 * mrSges(4,2) - t21 * mrSges(5,3) + t22 * mrSges(5,2) + (-t46 + t19) * qJ(4) + m(7) * (t1 * t84 + t111 * t5 + t2 * t83) + m(5) * (-pkin(3) * t22 - qJ(4) * t21) + t151 * t71 + t149 * t45 + t150 * t44 + m(6) * (qJ(4) * t17 - t142 * t174); qJ(4) * t75 + t111 * t74 + t69 * t88 + t84 * t79 + t83 * t81 + m(7) * (t111 * t69 + t33 * t84 + t34 * t83) + (-t37 * mrSges(6,3) - t33 * mrSges(7,3) - (m(6) * t37 + t80) * t174 + t152) * t135 + (-t34 * mrSges(7,3) - t38 * mrSges(6,3) - (m(6) * t38 + t82) * t174 - t153) * t132 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t151) * t133 + (qJ(4) * mrSges(5,1) - t149 * t132 - t150 * t135 - Ifges(5,5)) * t136 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t136 + (-mrSges(4,1) + t177) * t133) * pkin(9) + t158 + (m(6) * qJ(4) + t89) * t101; -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t111 * t88 + (t84 * t176 + t97 + t98) * t135 + (t83 * t176 - t94 - t95) * t132 + m(7) * (t111 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(6) * (t156 * t174 ^ 2 + t139) + m(5) * (pkin(3) ^ 2 + t139) + t170 + 0.2e1 * (t89 + mrSges(5,3)) * qJ(4) + 0.2e1 * t174 * t145; (t25 + t26) * t135 + (t23 + t24) * t132 + m(7) * (t1 * t135 + t132 * t2) + m(6) * t142 + m(5) * t22 + t47; (t79 + t80) * t135 + (m(5) * pkin(9) + mrSges(5,1)) * t133 + (t81 + t82) * t132 + m(7) * (t132 * t34 + t135 * t33) + m(6) * (t132 * t38 + t135 * t37); -t145 + m(7) * (t132 * t83 + t135 * t84) + (-m(6) * t174 - mrSges(7,3)) * t156 + t177; m(5) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t156; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t179; mrSges(6,1) * t37 + mrSges(7,1) * t33 - mrSges(6,2) * t38 - mrSges(7,2) * t34 + t115 + t116 + (m(7) * t33 + t79) * pkin(5) + (-t181 * t135 + (-Ifges(6,5) - Ifges(7,5)) * t132) * t136; -t174 * t168 + mrSges(7,1) * t84 - mrSges(7,2) * t83 + t118 + t119 + (m(7) * t84 - t135 * mrSges(7,3)) * pkin(5) + (mrSges(6,2) * t174 - t181) * t132; (-mrSges(6,2) - mrSges(7,2)) * t132 + (mrSges(6,1) + mrSges(7,1) + t175) * t135; Ifges(6,3) + Ifges(7,3) + (0.2e1 * mrSges(7,1) + t175) * pkin(5); m(7) * t5 + t18; m(7) * t69 + t74; m(7) * t111 + t88; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
