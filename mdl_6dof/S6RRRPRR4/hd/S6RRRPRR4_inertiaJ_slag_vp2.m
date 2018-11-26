% Calculate joint inertia matrix for
% S6RRRPRR4
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
% Datum: 2018-11-23 17:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:52:36
% EndTime: 2018-11-23 17:52:37
% DurationCPUTime: 1.32s
% Computational Cost: add. (3593->319), mult. (6857->446), div. (0->0), fcn. (7805->10), ass. (0->123)
t148 = sin(pkin(11));
t149 = cos(pkin(11));
t152 = sin(qJ(3));
t153 = sin(qJ(2));
t156 = cos(qJ(3));
t157 = cos(qJ(2));
t126 = t152 * t153 - t156 * t157;
t127 = t152 * t157 + t153 * t156;
t141 = -pkin(2) * t157 - pkin(1);
t83 = pkin(3) * t126 - qJ(4) * t127 + t141;
t194 = -pkin(8) - pkin(7);
t135 = t194 * t157;
t174 = t194 * t153;
t99 = -t156 * t135 + t152 * t174;
t52 = -t148 * t99 + t149 * t83;
t53 = t148 * t83 + t149 * t99;
t172 = -t148 * t52 + t149 * t53;
t150 = sin(qJ(6));
t151 = sin(qJ(5));
t155 = cos(qJ(5));
t124 = -t148 * t151 + t149 * t155;
t125 = t148 * t155 + t149 * t151;
t181 = Ifges(6,5) * t125 + Ifges(6,6) * t124;
t154 = cos(qJ(6));
t84 = t124 * t154 - t125 * t150;
t190 = t84 * mrSges(7,3);
t201 = t150 * pkin(5) * t190 + t181;
t72 = t125 * t127;
t73 = t124 * t127;
t200 = Ifges(6,5) * t73 - Ifges(6,6) * t72 + Ifges(6,3) * t126;
t97 = -t135 * t152 - t156 * t174;
t199 = t97 ^ 2;
t85 = t124 * t150 + t125 * t154;
t47 = -t84 * mrSges(7,1) + t85 * mrSges(7,2);
t198 = 0.2e1 * t47;
t88 = -t124 * mrSges(6,1) + t125 * mrSges(6,2);
t197 = 0.2e1 * t88;
t196 = 0.2e1 * t97;
t195 = 0.2e1 * t141;
t192 = pkin(2) * t156;
t191 = pkin(10) * t125;
t189 = t85 * mrSges(7,3);
t182 = t127 * t149;
t29 = pkin(4) * t126 - pkin(9) * t182 + t52;
t183 = t127 * t148;
t40 = -pkin(9) * t183 + t53;
t11 = t151 * t29 + t155 * t40;
t188 = Ifges(7,5) * t85 + Ifges(7,6) * t84;
t187 = Ifges(5,4) * t148;
t186 = Ifges(5,4) * t149;
t137 = pkin(2) * t152 + qJ(4);
t109 = (-pkin(9) - t137) * t148;
t143 = t149 * pkin(9);
t110 = t137 * t149 + t143;
t77 = t151 * t109 + t155 * t110;
t82 = mrSges(5,1) * t183 + mrSges(5,2) * t182;
t130 = (-pkin(9) - qJ(4)) * t148;
t132 = qJ(4) * t149 + t143;
t95 = t151 * t130 + t155 * t132;
t180 = t148 ^ 2 + t149 ^ 2;
t179 = t153 ^ 2 + t157 ^ 2;
t178 = 2 * mrSges(6,3);
t177 = 0.2e1 * mrSges(7,3);
t176 = t154 * t189;
t41 = -t150 * t73 - t154 * t72;
t42 = -t150 * t72 + t154 * t73;
t175 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t126;
t138 = -pkin(4) * t149 - pkin(3);
t44 = t72 * mrSges(6,1) + t73 * mrSges(6,2);
t14 = -t41 * mrSges(7,1) + t42 * mrSges(7,2);
t10 = -t151 * t40 + t155 * t29;
t131 = -t149 * mrSges(5,1) + t148 * mrSges(5,2);
t173 = t180 * t137;
t76 = t155 * t109 - t110 * t151;
t94 = t155 * t130 - t132 * t151;
t65 = pkin(4) * t183 + t97;
t86 = -mrSges(5,2) * t126 - mrSges(5,3) * t183;
t87 = mrSges(5,1) * t126 - mrSges(5,3) * t182;
t171 = -t148 * t87 + t149 * t86;
t58 = t76 - t191;
t114 = t124 * pkin(10);
t59 = t114 + t77;
t19 = -t150 * t59 + t154 * t58;
t20 = t150 * t58 + t154 * t59;
t170 = t19 * mrSges(7,1) - t20 * mrSges(7,2) + t188;
t66 = t94 - t191;
t67 = t114 + t95;
t33 = -t150 * t67 + t154 * t66;
t34 = t150 * t66 + t154 * t67;
t169 = t33 * mrSges(7,1) - t34 * mrSges(7,2) + t188;
t168 = 0.2e1 * mrSges(5,3) * t180;
t102 = -pkin(5) * t124 + t138;
t5 = pkin(5) * t126 - pkin(10) * t73 + t10;
t6 = -pkin(10) * t72 + t11;
t3 = -t150 * t6 + t154 * t5;
t4 = t150 * t5 + t154 * t6;
t167 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t175;
t166 = (mrSges(4,1) * t156 - mrSges(4,2) * t152) * pkin(2);
t165 = (mrSges(7,1) * t154 - mrSges(7,2) * t150) * pkin(5);
t133 = Ifges(5,2) * t149 + t187;
t134 = Ifges(5,1) * t148 + t186;
t48 = Ifges(7,4) * t85 + Ifges(7,2) * t84;
t49 = Ifges(7,1) * t85 + Ifges(7,4) * t84;
t89 = Ifges(6,4) * t125 + Ifges(6,2) * t124;
t90 = Ifges(6,1) * t125 + Ifges(6,4) * t124;
t164 = t124 * t89 + t125 * t90 + t149 * t133 + t148 * t134 + t84 * t48 + t85 * t49 + Ifges(4,3);
t163 = t131 + t47 + t88;
t12 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t126;
t13 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t126;
t30 = Ifges(6,4) * t73 - Ifges(6,2) * t72 + Ifges(6,6) * t126;
t31 = Ifges(6,1) * t73 - Ifges(6,4) * t72 + Ifges(6,5) * t126;
t43 = pkin(5) * t72 + t65;
t60 = Ifges(5,6) * t126 + (-Ifges(5,2) * t148 + t186) * t127;
t61 = Ifges(5,5) * t126 + (Ifges(5,1) * t149 - t187) * t127;
t162 = t149 * t60 / 0.2e1 + t148 * t61 / 0.2e1 + t124 * t30 / 0.2e1 + t125 * t31 / 0.2e1 - Ifges(4,6) * t126 + Ifges(4,5) * t127 + t65 * t88 - t72 * t89 / 0.2e1 + t73 * t90 / 0.2e1 - t99 * mrSges(4,2) + t84 * t12 / 0.2e1 + t85 * t13 / 0.2e1 + t43 * t47 + t41 * t48 / 0.2e1 + t42 * t49 / 0.2e1 - t133 * t183 / 0.2e1 + t134 * t182 / 0.2e1 - t3 * t189 + t4 * t190 + (t131 - mrSges(4,1)) * t97 + (-t10 * t125 + t11 * t124) * mrSges(6,3) + t172 * mrSges(5,3) + (Ifges(5,5) * t148 + Ifges(5,6) * t149 + t181 + t188) * t126 / 0.2e1;
t140 = -pkin(3) - t192;
t129 = t138 - t192;
t100 = t102 - t192;
t55 = mrSges(6,1) * t126 - mrSges(6,3) * t73;
t54 = -mrSges(6,2) * t126 - mrSges(6,3) * t72;
t28 = mrSges(7,1) * t126 - mrSges(7,3) * t42;
t27 = -mrSges(7,2) * t126 + mrSges(7,3) * t41;
t1 = [-0.2e1 * pkin(1) * (-t157 * mrSges(3,1) + t153 * mrSges(3,2)) + t157 * (Ifges(3,4) * t153 + Ifges(3,2) * t157) + t153 * (Ifges(3,1) * t153 + Ifges(3,4) * t157) + 0.2e1 * t52 * t87 - t72 * t30 + t73 * t31 + 0.2e1 * t53 * t86 + 0.2e1 * t65 * t44 + 0.2e1 * t43 * t14 + 0.2e1 * t11 * t54 + 0.2e1 * t10 * t55 + t41 * t12 + t42 * t13 + 0.2e1 * t4 * t27 + 0.2e1 * t3 * t28 + m(3) * (pkin(7) ^ 2 * t179 + pkin(1) ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t65 ^ 2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t43 ^ 2) + m(4) * (t141 ^ 2 + t99 ^ 2 + t199) + m(5) * (t52 ^ 2 + t53 ^ 2 + t199) + (mrSges(4,2) * t195 + mrSges(4,3) * t196 + Ifges(4,1) * t127 - t148 * t60 + t149 * t61) * t127 + 0.2e1 * t179 * pkin(7) * mrSges(3,3) + (mrSges(4,1) * t195 - 0.2e1 * t99 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t126 + (Ifges(5,5) * t149 - Ifges(5,6) * t148 - (2 * Ifges(4,4))) * t127 + t175 + t200) * t126 + Ifges(2,3) + t82 * t196; Ifges(3,6) * t157 + Ifges(3,5) * t153 + t140 * t82 + t129 * t44 + t100 * t14 + t76 * t55 + t77 * t54 + t20 * t27 + t19 * t28 + m(5) * (t137 * t172 + t140 * t97) + t162 + t171 * t137 + m(6) * (t10 * t76 + t11 * t77 + t129 * t65) + m(7) * (t100 * t43 + t19 * t3 + t20 * t4) + (-mrSges(3,1) * t153 - mrSges(3,2) * t157) * pkin(7) + (m(4) * (t152 * t99 - t156 * t97) + (-t126 * t152 - t127 * t156) * mrSges(4,3)) * pkin(2); 0.2e1 * t140 * t131 + t129 * t197 + t100 * t198 + m(5) * (t137 ^ 2 * t180 + t140 ^ 2) + (-t19 * t85 + t20 * t84) * t177 + 0.2e1 * t166 + (t77 * t124 - t76 * t125) * t178 + t164 + m(4) * (t152 ^ 2 + t156 ^ 2) * pkin(2) ^ 2 + t137 * t168 + m(7) * (t100 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(6) * (t129 ^ 2 + t76 ^ 2 + t77 ^ 2) + Ifges(3,3); t138 * t44 + t94 * t55 + t95 * t54 + t102 * t14 - pkin(3) * t82 + t33 * t28 + t34 * t27 + m(5) * (-pkin(3) * t97 + qJ(4) * t172) + t162 + t171 * qJ(4) + m(6) * (t10 * t94 + t11 * t95 + t138 * t65) + m(7) * (t102 * t43 + t3 * t33 + t34 * t4); m(5) * (-pkin(3) * t140 + qJ(4) * t173) + t164 + (qJ(4) * t180 + t173) * mrSges(5,3) + ((-t76 - t94) * t125 + (t77 + t95) * t124) * mrSges(6,3) + ((-t19 - t33) * t85 + (t20 + t34) * t84) * mrSges(7,3) + m(6) * (t129 * t138 + t76 * t94 + t77 * t95) + m(7) * (t100 * t102 + t19 * t33 + t20 * t34) + (t129 + t138) * t88 + (t100 + t102) * t47 + (t140 - pkin(3)) * t131 + t166; -0.2e1 * pkin(3) * t131 + t102 * t198 + t138 * t197 + (-t33 * t85 + t34 * t84) * t177 + (t95 * t124 - t94 * t125) * t178 + qJ(4) * t168 + m(7) * (t102 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t138 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(5) * (qJ(4) ^ 2 * t180 + pkin(3) ^ 2) + t164; m(5) * t97 + m(6) * t65 + m(7) * t43 + t14 + t44 + t82; m(5) * t140 + m(6) * t129 + m(7) * t100 + t163; -m(5) * pkin(3) + m(6) * t138 + m(7) * t102 + t163; m(5) + m(6) + m(7); t10 * mrSges(6,1) - t11 * mrSges(6,2) + (m(7) * (t150 * t4 + t154 * t3) + t154 * t28 + t150 * t27) * pkin(5) + t167 + t200; t76 * mrSges(6,1) - t77 * mrSges(6,2) + (m(7) * (t150 * t20 + t154 * t19) - t176) * pkin(5) + t170 + t201; t94 * mrSges(6,1) - t95 * mrSges(6,2) + (m(7) * (t150 * t34 + t154 * t33) - t176) * pkin(5) + t169 + t201; 0; Ifges(6,3) + Ifges(7,3) + m(7) * (t150 ^ 2 + t154 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t165; t167; t170; t169; 0; Ifges(7,3) + t165; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
