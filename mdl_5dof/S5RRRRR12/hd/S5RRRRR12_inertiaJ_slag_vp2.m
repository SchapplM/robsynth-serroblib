% Calculate joint inertia matrix for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:39
% EndTime: 2019-12-31 22:46:46
% DurationCPUTime: 1.77s
% Computational Cost: add. (3449->383), mult. (9014->574), div. (0->0), fcn. (9937->12), ass. (0->153)
t198 = 2 * pkin(10);
t144 = sin(pkin(6));
t146 = cos(pkin(6));
t150 = sin(qJ(3));
t154 = cos(qJ(3));
t151 = sin(qJ(2));
t145 = sin(pkin(5));
t155 = cos(qJ(2));
t173 = t145 * t155;
t147 = cos(pkin(5));
t186 = pkin(1) * t147;
t104 = pkin(8) * t173 + t151 * t186;
t171 = t146 * t155;
t165 = t145 * t171;
t63 = (t144 * t147 + t165) * pkin(9) + t104;
t130 = t155 * t186;
t174 = t145 * t151;
t71 = t147 * pkin(2) + t130 + (-pkin(9) * t146 - pkin(8)) * t174;
t80 = (-pkin(9) * t144 * t151 - pkin(2) * t155 - pkin(1)) * t145;
t25 = -t150 * t63 + (t144 * t80 + t146 * t71) * t154;
t149 = sin(qJ(4));
t153 = cos(qJ(4));
t176 = t144 * t150;
t68 = t147 * t176 + (t150 * t171 + t151 * t154) * t145;
t97 = -t144 * t173 + t147 * t146;
t45 = t68 * t149 - t97 * t153;
t46 = t97 * t149 + t68 * t153;
t175 = t144 * t154;
t67 = -t147 * t175 + t150 * t174 - t154 * t165;
t16 = Ifges(5,1) * t46 - Ifges(5,4) * t45 + Ifges(5,5) * t67;
t197 = t16 / 0.2e1;
t98 = -t153 * t146 + t149 * t176;
t99 = t149 * t146 + t153 * t176;
t58 = Ifges(5,1) * t99 - Ifges(5,4) * t98 - Ifges(5,5) * t175;
t196 = t58 / 0.2e1;
t148 = sin(qJ(5));
t152 = cos(qJ(5));
t179 = Ifges(6,4) * t152;
t91 = -Ifges(6,6) * t153 + (-Ifges(6,2) * t148 + t179) * t149;
t195 = t91 / 0.2e1;
t180 = Ifges(6,4) * t148;
t92 = -Ifges(6,5) * t153 + (Ifges(6,1) * t152 - t180) * t149;
t194 = t92 / 0.2e1;
t114 = Ifges(6,5) * t148 + Ifges(6,6) * t152;
t193 = t114 / 0.2e1;
t116 = Ifges(6,2) * t152 + t180;
t192 = t116 / 0.2e1;
t118 = Ifges(6,1) * t148 + t179;
t191 = t118 / 0.2e1;
t119 = Ifges(5,1) * t149 + Ifges(5,4) * t153;
t190 = t119 / 0.2e1;
t189 = -t148 / 0.2e1;
t188 = t148 / 0.2e1;
t187 = t152 / 0.2e1;
t185 = pkin(2) * t154;
t184 = pkin(10) * t149;
t183 = pkin(10) * t153;
t182 = pkin(11) * t148;
t181 = pkin(11) * t152;
t41 = -t144 * t71 + t146 * t80;
t18 = t67 * pkin(3) - t68 * pkin(10) + t41;
t172 = t146 * t150;
t26 = t154 * t63 + t71 * t172 + t80 * t176;
t21 = t97 * pkin(10) + t26;
t6 = t149 * t18 + t153 * t21;
t103 = pkin(2) * t172 + pkin(9) * t175;
t88 = t146 * pkin(10) + t103;
t89 = (-pkin(3) * t154 - pkin(10) * t150 - pkin(2)) * t144;
t55 = t149 * t89 + t153 * t88;
t102 = -pkin(8) * t174 + t130;
t178 = t102 * mrSges(3,1);
t177 = t104 * mrSges(3,2);
t170 = t148 * t149;
t169 = t149 * t152;
t115 = Ifges(5,5) * t149 + Ifges(5,6) * t153;
t168 = t148 ^ 2 + t152 ^ 2;
t29 = -t46 * t148 + t67 * t152;
t30 = t67 * t148 + t46 * t152;
t7 = Ifges(6,5) * t30 + Ifges(6,6) * t29 + Ifges(6,3) * t45;
t14 = Ifges(5,5) * t46 - Ifges(5,6) * t45 + Ifges(5,3) * t67;
t33 = Ifges(4,5) * t68 - Ifges(4,6) * t67 + Ifges(4,3) * t97;
t72 = -t148 * t99 - t152 * t175;
t73 = -t148 * t175 + t152 * t99;
t36 = Ifges(6,5) * t73 + Ifges(6,6) * t72 + Ifges(6,3) * t98;
t15 = Ifges(5,4) * t46 - Ifges(5,2) * t45 + Ifges(5,6) * t67;
t167 = t7 / 0.2e1 - t15 / 0.2e1;
t57 = Ifges(5,4) * t99 - Ifges(5,2) * t98 - Ifges(5,6) * t175;
t166 = -t57 / 0.2e1 + t36 / 0.2e1;
t117 = Ifges(5,4) * t149 + Ifges(5,2) * t153;
t90 = Ifges(6,5) * t169 - Ifges(6,6) * t170 - Ifges(6,3) * t153;
t164 = t90 / 0.2e1 - t117 / 0.2e1;
t83 = Ifges(4,5) * t176 + Ifges(4,6) * t175 + Ifges(4,3) * t146;
t163 = Ifges(3,5) * t174 + Ifges(3,6) * t173 + Ifges(3,3) * t147;
t20 = -t97 * pkin(3) - t25;
t10 = t45 * pkin(4) - t46 * pkin(11) + t20;
t4 = t67 * pkin(11) + t6;
t1 = t152 * t10 - t148 * t4;
t2 = t148 * t10 + t152 * t4;
t162 = -t1 * t148 + t2 * t152;
t161 = mrSges(6,1) * t148 + mrSges(6,2) * t152;
t125 = pkin(9) * t176;
t87 = t125 + (-pkin(3) - t185) * t146;
t47 = t98 * pkin(4) - t99 * pkin(11) + t87;
t49 = -pkin(11) * t175 + t55;
t22 = -t148 * t49 + t152 * t47;
t23 = t148 * t47 + t152 * t49;
t159 = -t22 * t148 + t23 * t152;
t111 = -t153 * pkin(4) - t149 * pkin(11) - pkin(3);
t81 = t152 * t111 - t148 * t183;
t82 = t148 * t111 + t152 * t183;
t158 = -t81 * t148 + t82 * t152;
t5 = -t149 * t21 + t153 * t18;
t54 = -t149 * t88 + t153 * t89;
t56 = Ifges(5,5) * t99 - Ifges(5,6) * t98 - Ifges(5,3) * t175;
t157 = pkin(10) ^ 2;
t143 = t153 ^ 2;
t141 = t149 ^ 2;
t139 = t141 * t157;
t113 = -t153 * mrSges(5,1) + t149 * mrSges(5,2);
t112 = -t152 * mrSges(6,1) + t148 * mrSges(6,2);
t109 = -t153 * mrSges(6,1) - mrSges(6,3) * t169;
t108 = t153 * mrSges(6,2) - mrSges(6,3) * t170;
t107 = -t146 * mrSges(4,2) + mrSges(4,3) * t175;
t106 = t146 * mrSges(4,1) - mrSges(4,3) * t176;
t105 = t161 * t149;
t101 = t146 * t185 - t125;
t100 = (-mrSges(4,1) * t154 + mrSges(4,2) * t150) * t144;
t85 = Ifges(4,5) * t146 + (Ifges(4,1) * t150 + Ifges(4,4) * t154) * t144;
t84 = Ifges(4,6) * t146 + (Ifges(4,4) * t150 + Ifges(4,2) * t154) * t144;
t77 = -mrSges(5,1) * t175 - t99 * mrSges(5,3);
t76 = mrSges(5,2) * t175 - t98 * mrSges(5,3);
t59 = t98 * mrSges(5,1) + t99 * mrSges(5,2);
t53 = t98 * mrSges(6,1) - t73 * mrSges(6,3);
t52 = -t98 * mrSges(6,2) + t72 * mrSges(6,3);
t51 = t97 * mrSges(4,1) - t68 * mrSges(4,3);
t50 = -t97 * mrSges(4,2) - t67 * mrSges(4,3);
t48 = pkin(4) * t175 - t54;
t40 = -t72 * mrSges(6,1) + t73 * mrSges(6,2);
t39 = t67 * mrSges(4,1) + t68 * mrSges(4,2);
t38 = Ifges(6,1) * t73 + Ifges(6,4) * t72 + Ifges(6,5) * t98;
t37 = Ifges(6,4) * t73 + Ifges(6,2) * t72 + Ifges(6,6) * t98;
t35 = Ifges(4,1) * t68 - Ifges(4,4) * t67 + Ifges(4,5) * t97;
t34 = Ifges(4,4) * t68 - Ifges(4,2) * t67 + Ifges(4,6) * t97;
t32 = t67 * mrSges(5,1) - t46 * mrSges(5,3);
t31 = -t67 * mrSges(5,2) - t45 * mrSges(5,3);
t24 = t45 * mrSges(5,1) + t46 * mrSges(5,2);
t13 = t45 * mrSges(6,1) - t30 * mrSges(6,3);
t12 = -t45 * mrSges(6,2) + t29 * mrSges(6,3);
t11 = -t29 * mrSges(6,1) + t30 * mrSges(6,2);
t9 = Ifges(6,1) * t30 + Ifges(6,4) * t29 + Ifges(6,5) * t45;
t8 = Ifges(6,4) * t30 + Ifges(6,2) * t29 + Ifges(6,6) * t45;
t3 = -t67 * pkin(4) - t5;
t17 = [0.2e1 * t1 * t13 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + t46 * t16 + 0.2e1 * t20 * t24 + 0.2e1 * t25 * t51 + 0.2e1 * t26 * t50 + t29 * t8 + t30 * t9 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + t97 * t33 + t68 * t35 + 0.2e1 * t41 * t39 + Ifges(2,3) + (t14 - t34) * t67 + (t7 - t15) * t45 + (t163 - 0.2e1 * t177 + 0.2e1 * t178) * t147 + m(4) * (t25 ^ 2 + t26 ^ 2 + t41 ^ 2) + m(5) * (t20 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(3) * (t102 ^ 2 + t104 ^ 2) + ((t151 * Ifges(3,5) + t155 * Ifges(3,6)) * t147 + 0.2e1 * (-t102 * t151 + t104 * t155) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t155 + mrSges(3,2) * t151) + t155 * (Ifges(3,4) * t151 + Ifges(3,2) * t155) + t151 * (Ifges(3,1) * t151 + Ifges(3,4) * t155) + m(3) * pkin(1) ^ 2) * t145) * t145; -t177 + t178 + t163 + t72 * t8 / 0.2e1 + t73 * t9 / 0.2e1 + t6 * t76 + t5 * t77 + t68 * t85 / 0.2e1 + t2 * t52 + t1 * t53 + t54 * t32 + t55 * t31 + t20 * t59 + t48 * t11 + t29 * t37 / 0.2e1 + t30 * t38 / 0.2e1 + t3 * t40 + t22 * t13 + t23 * t12 + t167 * t98 + t166 * t45 + t46 * t196 + t99 * t197 + (-pkin(2) * t39 + t150 * t35 / 0.2e1 + (t34 / 0.2e1 - t14 / 0.2e1) * t154) * t144 + m(4) * (-t144 * pkin(2) * t41 + t101 * t25 + t103 * t26) + m(5) * (t87 * t20 + t54 * t5 + t55 * t6) + m(6) * (t22 * t1 + t23 * t2 + t48 * t3) + (-t84 / 0.2e1 + t56 / 0.2e1) * t67 + t87 * t24 + t97 * t83 / 0.2e1 + t41 * t100 + t101 * t51 + t103 * t50 + t25 * t106 + t26 * t107 + t146 * t33 / 0.2e1; 0.2e1 * t101 * t106 + 0.2e1 * t103 * t107 + t146 * t83 + 0.2e1 * t22 * t53 + 0.2e1 * t23 * t52 + t72 * t37 + t73 * t38 + 0.2e1 * t48 * t40 + 0.2e1 * t54 * t77 + 0.2e1 * t55 * t76 + t99 * t58 + 0.2e1 * t87 * t59 + Ifges(3,3) + (t36 - t57) * t98 + (-0.2e1 * pkin(2) * t100 + t150 * t85 + (-t56 + t84) * t154) * t144 + m(6) * (t22 ^ 2 + t23 ^ 2 + t48 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2 + t87 ^ 2) + m(4) * (t144 ^ 2 * pkin(2) ^ 2 + t101 ^ 2 + t103 ^ 2); t81 * t13 + t82 * t12 + t25 * mrSges(4,1) - t26 * mrSges(4,2) - pkin(3) * t24 + (-t5 * mrSges(5,3) + t8 * t189 + t9 * t187 + t197 + (-t32 + t11) * pkin(10)) * t149 + (t6 * mrSges(5,3) + pkin(10) * t31 - t167) * t153 + t164 * t45 + m(5) * (-pkin(3) * t20 + (-t5 * t149 + t6 * t153) * pkin(10)) + m(6) * (t81 * t1 + t3 * t184 + t82 * t2) + t29 * t195 + t30 * t194 + t3 * t105 + t2 * t108 + t1 * t109 + t20 * t113 + t67 * t115 / 0.2e1 + t46 * t190 + t33; t81 * t53 + t82 * t52 - pkin(3) * t59 + (-t54 * mrSges(5,3) + t37 * t189 + t38 * t187 + t196 + (-t77 + t40) * pkin(10)) * t149 + (t55 * mrSges(5,3) + pkin(10) * t76 - t166) * t153 + m(6) * (t48 * t184 + t81 * t22 + t82 * t23) + t164 * t98 + m(5) * (-pkin(3) * t87 + (-t54 * t149 + t55 * t153) * pkin(10)) + t83 - t115 * t175 / 0.2e1 + t72 * t195 + t73 * t194 + t101 * mrSges(4,1) - t103 * mrSges(4,2) + t48 * t105 + t23 * t108 + t22 * t109 + t87 * t113 + t99 * t190; -0.2e1 * pkin(3) * t113 + 0.2e1 * t82 * t108 + 0.2e1 * t81 * t109 + Ifges(4,3) + (-t90 + t117) * t153 + (t141 + t143) * mrSges(5,3) * t198 + m(6) * (t81 ^ 2 + t82 ^ 2 + t139) + m(5) * (pkin(3) ^ 2 + t143 * t157 + t139) + (t105 * t198 - t148 * t91 + t152 * t92 + t119) * t149; t9 * t188 + t8 * t187 + t3 * t112 + m(6) * (-pkin(4) * t3 + t162 * pkin(11)) + t12 * t181 - t13 * t182 - pkin(4) * t11 + t30 * t191 + t29 * t192 + t45 * t193 - t6 * mrSges(5,2) + t5 * mrSges(5,1) + t162 * mrSges(6,3) + t14; t98 * t193 - pkin(4) * t40 - t53 * t182 + t38 * t188 + t37 * t187 + t48 * t112 + t52 * t181 + m(6) * (-pkin(4) * t48 + t159 * pkin(11)) + t73 * t191 + t72 * t192 - t55 * mrSges(5,2) + t54 * mrSges(5,1) + t159 * mrSges(6,3) + t56; t92 * t188 + t91 * t187 - pkin(4) * t105 + (m(6) * t158 + t152 * t108 - t148 * t109) * pkin(11) + (t118 * t187 + t116 * t189 + (-m(6) * pkin(4) - mrSges(5,1) + t112) * pkin(10)) * t149 + (-t114 / 0.2e1 - pkin(10) * mrSges(5,2)) * t153 + t158 * mrSges(6,3) + t115; Ifges(5,3) - 0.2e1 * pkin(4) * t112 + m(6) * (t168 * pkin(11) ^ 2 + pkin(4) ^ 2) + t148 * t118 + t152 * t116 + 0.2e1 * t168 * pkin(11) * mrSges(6,3); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t7; t22 * mrSges(6,1) - t23 * mrSges(6,2) + t36; t81 * mrSges(6,1) - t82 * mrSges(6,2) + t90; -t161 * pkin(11) + t114; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
