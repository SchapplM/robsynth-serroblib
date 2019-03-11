% Calculate joint inertia matrix for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:55
% EndTime: 2019-03-09 21:42:00
% DurationCPUTime: 2.08s
% Computational Cost: add. (2229->442), mult. (5007->575), div. (0->0), fcn. (5042->8), ass. (0->152)
t202 = Ifges(7,4) + Ifges(6,5);
t149 = sin(qJ(4));
t152 = cos(qJ(4));
t201 = t149 ^ 2 + t152 ^ 2;
t200 = 2 * pkin(9);
t153 = cos(qJ(3));
t141 = t153 * pkin(4);
t192 = pkin(9) * t153;
t127 = t149 * t192;
t150 = sin(qJ(3));
t95 = -pkin(3) * t153 - pkin(10) * t150 - pkin(2);
t57 = t152 * t95 - t127;
t53 = t141 - t57;
t176 = t150 * t152;
t93 = mrSges(6,1) * t176 - t153 * mrSges(6,2);
t198 = m(6) * t53 + t93;
t146 = sin(pkin(6));
t154 = cos(qJ(2));
t178 = t146 * t154;
t147 = cos(pkin(6));
t151 = sin(qJ(2));
t179 = t146 * t151;
t80 = t147 * t150 + t153 * t179;
t47 = t149 * t80 + t152 * t178;
t48 = -t149 * t178 + t152 * t80;
t79 = -t147 * t153 + t150 * t179;
t11 = Ifges(5,5) * t48 - Ifges(5,6) * t47 + Ifges(5,3) * t79;
t15 = Ifges(7,1) * t79 + Ifges(7,4) * t47 + Ifges(7,5) * t48;
t16 = Ifges(6,1) * t79 - Ifges(6,4) * t48 + Ifges(6,5) * t47;
t197 = t11 + t15 + t16;
t196 = 2 * mrSges(7,1);
t195 = pkin(5) + pkin(10);
t193 = pkin(1) * t154;
t139 = t150 * pkin(9);
t116 = pkin(8) * t179;
t81 = t147 * t193 - t116;
t191 = t81 * mrSges(3,1);
t82 = pkin(1) * t147 * t151 + pkin(8) * t178;
t190 = t82 * mrSges(3,2);
t189 = mrSges(6,1) + mrSges(7,1);
t188 = pkin(4) + qJ(6);
t59 = t116 + (-pkin(2) - t193) * t147;
t24 = pkin(3) * t79 - pkin(10) * t80 + t59;
t60 = pkin(9) * t147 + t82;
t61 = (-pkin(2) * t154 - pkin(9) * t151 - pkin(1)) * t146;
t34 = t150 * t61 + t153 * t60;
t26 = -pkin(10) * t178 + t34;
t7 = t149 * t24 + t152 * t26;
t187 = -Ifges(4,5) * t80 + Ifges(4,6) * t79;
t186 = Ifges(5,4) * t149;
t185 = Ifges(5,4) * t152;
t184 = Ifges(6,6) * t149;
t183 = Ifges(6,6) * t152;
t182 = Ifges(7,6) * t149;
t181 = Ifges(7,6) * t152;
t58 = t149 * t95 + t152 * t192;
t180 = qJ(5) * t152;
t177 = t149 * t150;
t175 = Ifges(7,4) * t177 + Ifges(7,5) * t176;
t174 = pkin(4) * t177 + t139;
t90 = mrSges(7,1) * t176 + mrSges(7,3) * t153;
t104 = Ifges(5,5) * t149 + Ifges(5,6) * t152;
t173 = Ifges(4,5) * t150 + Ifges(4,6) * t153;
t172 = t201 * pkin(10) ^ 2;
t171 = -Ifges(5,3) - Ifges(7,1) - Ifges(6,1);
t170 = Ifges(3,5) * t179 + Ifges(3,6) * t178 + Ifges(3,3) * t147;
t30 = t48 * mrSges(6,1) + mrSges(6,2) * t79;
t29 = -mrSges(7,1) * t47 + mrSges(7,2) * t79;
t27 = mrSges(7,1) * t48 - t79 * mrSges(7,3);
t169 = -qJ(5) * t149 - pkin(3);
t6 = -t149 * t26 + t152 * t24;
t33 = -t150 * t60 + t153 * t61;
t13 = Ifges(6,4) * t79 - Ifges(6,2) * t48 + Ifges(6,6) * t47;
t17 = Ifges(5,1) * t48 - Ifges(5,4) * t47 + Ifges(5,5) * t79;
t9 = Ifges(7,5) * t79 + Ifges(7,6) * t47 + Ifges(7,3) * t48;
t167 = t9 / 0.2e1 - t13 / 0.2e1 + t17 / 0.2e1;
t4 = -qJ(5) * t79 - t7;
t10 = Ifges(6,5) * t79 - Ifges(6,6) * t48 + Ifges(6,3) * t47;
t12 = Ifges(7,4) * t79 + Ifges(7,2) * t47 + Ifges(7,6) * t48;
t14 = Ifges(5,4) * t48 - Ifges(5,2) * t47 + Ifges(5,6) * t79;
t166 = t10 / 0.2e1 + t12 / 0.2e1 - t14 / 0.2e1;
t63 = -Ifges(5,6) * t153 + (-Ifges(5,2) * t149 + t185) * t150;
t66 = -Ifges(6,5) * t153 + (Ifges(6,3) * t149 - t183) * t150;
t67 = -Ifges(7,4) * t153 + (Ifges(7,2) * t149 + t181) * t150;
t165 = -t63 / 0.2e1 + t66 / 0.2e1 + t67 / 0.2e1;
t64 = -Ifges(5,5) * t153 + (Ifges(5,1) * t152 - t186) * t150;
t65 = -Ifges(7,5) * t153 + (Ifges(7,3) * t152 + t182) * t150;
t68 = -Ifges(6,4) * t153 + (-Ifges(6,2) * t152 + t184) * t150;
t164 = t64 / 0.2e1 + t65 / 0.2e1 - t68 / 0.2e1;
t25 = pkin(3) * t178 - t33;
t52 = qJ(5) * t153 - t58;
t100 = Ifges(7,3) * t149 - t181;
t103 = -Ifges(6,2) * t149 - t183;
t109 = Ifges(5,1) * t149 + t185;
t163 = -t103 / 0.2e1 + t109 / 0.2e1 + t100 / 0.2e1;
t134 = Ifges(7,5) * t149;
t162 = t104 / 0.2e1 + t134 / 0.2e1 - Ifges(6,4) * t149 / 0.2e1 - t202 * t152 / 0.2e1;
t101 = -Ifges(6,3) * t152 - t184;
t102 = -Ifges(7,2) * t152 + t182;
t107 = Ifges(5,2) * t152 + t186;
t161 = -t107 / 0.2e1 + t101 / 0.2e1 + t102 / 0.2e1;
t160 = t149 * mrSges(5,1) + t152 * mrSges(5,2);
t159 = -t149 * mrSges(6,2) - t152 * mrSges(6,3);
t92 = -mrSges(7,1) * t177 - mrSges(7,2) * t153;
t158 = -qJ(5) * t48 + t25;
t157 = pkin(9) ^ 2;
t155 = qJ(5) ^ 2;
t145 = t153 ^ 2;
t143 = t150 ^ 2;
t138 = t143 * t157;
t124 = Ifges(5,5) * t176;
t123 = Ifges(6,5) * t177;
t112 = t195 * t152;
t111 = t195 * t149;
t110 = Ifges(4,1) * t150 + Ifges(4,4) * t153;
t108 = Ifges(4,4) * t150 + Ifges(4,2) * t153;
t99 = -mrSges(7,2) * t149 - mrSges(7,3) * t152;
t98 = -mrSges(4,1) * t153 + mrSges(4,2) * t150;
t97 = -mrSges(5,1) * t152 + mrSges(5,2) * t149;
t96 = mrSges(6,2) * t152 - mrSges(6,3) * t149;
t94 = -pkin(4) * t152 + t169;
t91 = mrSges(6,1) * t177 + mrSges(6,3) * t153;
t89 = -mrSges(5,1) * t153 - mrSges(5,3) * t176;
t88 = mrSges(5,2) * t153 - mrSges(5,3) * t177;
t86 = t160 * t150;
t85 = t159 * t150;
t84 = (-mrSges(7,2) * t152 + mrSges(7,3) * t149) * t150;
t83 = -t152 * t188 + t169;
t71 = -qJ(5) * t176 + t174;
t70 = -Ifges(6,1) * t153 - Ifges(6,4) * t176 + t123;
t69 = -Ifges(7,1) * t153 + t175;
t62 = -Ifges(5,6) * t177 - Ifges(5,3) * t153 + t124;
t51 = -mrSges(4,1) * t178 - mrSges(4,3) * t80;
t50 = mrSges(4,2) * t178 - mrSges(4,3) * t79;
t49 = (qJ(6) * t149 - t180) * t150 + t174;
t39 = -pkin(5) * t177 - t52;
t38 = mrSges(4,1) * t79 + mrSges(4,2) * t80;
t37 = qJ(6) * t153 + t127 + t141 + (pkin(5) * t150 - t95) * t152;
t36 = Ifges(4,1) * t80 - Ifges(4,4) * t79 - Ifges(4,5) * t178;
t35 = Ifges(4,4) * t80 - Ifges(4,2) * t79 - Ifges(4,6) * t178;
t32 = mrSges(5,1) * t79 - mrSges(5,3) * t48;
t31 = -mrSges(5,2) * t79 - mrSges(5,3) * t47;
t28 = mrSges(6,1) * t47 - mrSges(6,3) * t79;
t20 = -mrSges(6,2) * t47 - mrSges(6,3) * t48;
t19 = mrSges(5,1) * t47 + mrSges(5,2) * t48;
t18 = -mrSges(7,2) * t48 + mrSges(7,3) * t47;
t8 = pkin(4) * t47 + t158;
t5 = -pkin(4) * t79 - t6;
t3 = t188 * t47 + t158;
t2 = -pkin(5) * t47 - t4;
t1 = pkin(5) * t48 - t188 * t79 - t6;
t21 = [0.2e1 * t1 * t27 + 0.2e1 * t3 * t18 + 0.2e1 * t25 * t19 + 0.2e1 * t2 * t29 + 0.2e1 * t8 * t20 + 0.2e1 * t4 * t28 + 0.2e1 * t5 * t30 + 0.2e1 * t7 * t31 + 0.2e1 * t6 * t32 + 0.2e1 * t33 * t51 + 0.2e1 * t34 * t50 + t80 * t36 + 0.2e1 * t59 * t38 + Ifges(2,3) + (t9 + t17 - t13) * t48 + (t10 + t12 - t14) * t47 + (t170 - 0.2e1 * t190 + 0.2e1 * t191) * t147 + (-t35 + t197) * t79 + ((-0.2e1 * t81 * mrSges(3,3) + Ifges(3,5) * t147 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t151) * t146) * t151 + (0.2e1 * t82 * mrSges(3,3) + Ifges(3,6) * t147 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t151 + (Ifges(3,2) + Ifges(4,3)) * t154) * t146 + t187) * t154) * t146 + m(3) * (pkin(1) ^ 2 * t146 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2 + t59 ^ 2) + m(5) * (t25 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2); t59 * t98 + t80 * t110 / 0.2e1 + t170 + (t35 / 0.2e1 - t11 / 0.2e1 - t15 / 0.2e1 - t16 / 0.2e1 + pkin(9) * t50 + t34 * mrSges(4,3)) * t153 + m(6) * (t4 * t52 + t5 * t53 + t71 * t8) + m(7) * (t1 * t37 + t2 * t39 + t3 * t49) + (-t108 / 0.2e1 + t62 / 0.2e1 + t69 / 0.2e1 + t70 / 0.2e1) * t79 + t3 * t84 + t8 * t85 + t25 * t86 + t7 * t88 + t6 * t89 + t1 * t90 + t4 * t91 + t2 * t92 + t5 * t93 + t57 * t32 + t58 * t31 + t71 * t20 + t49 * t18 + t52 * t28 + t53 * t30 + t37 * t27 - pkin(2) * t38 + t39 * t29 + m(4) * (-pkin(2) * t59 + (-t150 * t33 + t153 * t34) * pkin(9)) - t173 * t178 / 0.2e1 + t165 * t47 + (t36 / 0.2e1 - t33 * mrSges(4,3) + (-t51 + t19) * pkin(9) + t167 * t152 + t166 * t149) * t150 + t164 * t48 - t190 + t191 + m(5) * (t139 * t25 + t57 * t6 + t58 * t7); -0.2e1 * pkin(2) * t98 + 0.2e1 * t37 * t90 + 0.2e1 * t39 * t92 + 0.2e1 * t49 * t84 + 0.2e1 * t52 * t91 + 0.2e1 * t53 * t93 + 0.2e1 * t57 * t89 + 0.2e1 * t58 * t88 + 0.2e1 * t71 * t85 + Ifges(3,3) + (t143 + t145) * mrSges(4,3) * t200 + (t108 - t62 - t69 - t70) * t153 + m(4) * (pkin(2) ^ 2 + t145 * t157 + t138) + m(5) * (t57 ^ 2 + t58 ^ 2 + t138) + m(6) * (t52 ^ 2 + t53 ^ 2 + t71 ^ 2) + m(7) * (t37 ^ 2 + t39 ^ 2 + t49 ^ 2) + (t86 * t200 + t110 + (t64 + t65 - t68) * t152 + (-t63 + t66 + t67) * t149) * t150; -Ifges(4,3) * t178 + t33 * mrSges(4,1) - t34 * mrSges(4,2) - pkin(3) * t19 + t111 * t27 + t112 * t29 + t83 * t18 + t94 * t20 + t25 * t97 + t3 * t99 + t8 * t96 + t162 * t79 + t163 * t48 + t161 * t47 + (t7 * mrSges(5,3) + t2 * mrSges(7,1) - t4 * mrSges(6,1) + (-t28 + t31) * pkin(10) - t166) * t152 + (-t6 * mrSges(5,3) + t1 * mrSges(7,1) + t5 * mrSges(6,1) + (t30 - t32) * pkin(10) + t167) * t149 + m(6) * (t8 * t94 + (t149 * t5 - t152 * t4) * pkin(10)) + m(5) * (-pkin(3) * t25 + (-t149 * t6 + t152 * t7) * pkin(10)) + m(7) * (t1 * t111 + t112 * t2 + t3 * t83) - t187; -pkin(3) * t86 + t111 * t90 + t112 * t92 + t49 * t99 + t83 * t84 + t94 * t85 + m(7) * (t111 * t37 + t112 * t39 + t49 * t83) - t162 * t153 + (-t153 * mrSges(4,2) + (-m(5) * pkin(3) - mrSges(4,1) + t97) * t150) * pkin(9) + (t58 * mrSges(5,3) - t52 * mrSges(6,1) + t39 * mrSges(7,1) + t163 * t150 + (m(5) * t58 - m(6) * t52 + t88 - t91) * pkin(10) - t165) * t152 + (t53 * mrSges(6,1) - t57 * mrSges(5,3) + t37 * mrSges(7,1) + t161 * t150 + (-m(5) * t57 + t198 - t89) * pkin(10) + t164) * t149 + t173 + (m(6) * t94 + t96) * t71; -0.2e1 * pkin(3) * t97 + 0.2e1 * t83 * t99 + 0.2e1 * t94 * t96 + Ifges(4,3) + m(6) * (t94 ^ 2 + t172) + m(7) * (t111 ^ 2 + t112 ^ 2 + t83 ^ 2) + m(5) * (pkin(3) ^ 2 + t172) + (t112 * t196 - t101 - t102 + t107) * t152 + (t111 * t196 + t100 - t103 + t109) * t149 + 0.2e1 * (mrSges(6,1) + mrSges(5,3)) * pkin(10) * t201; t6 * mrSges(5,1) - t7 * mrSges(5,2) + t5 * mrSges(6,2) + t2 * mrSges(7,2) - t4 * mrSges(6,3) - t1 * mrSges(7,3) - pkin(4) * t30 - t188 * t27 + (-t28 + t29) * qJ(5) + m(6) * (-pkin(4) * t5 - qJ(5) * t4) + m(7) * (qJ(5) * t2 - t1 * t188) + t197; t57 * mrSges(5,1) - t58 * mrSges(5,2) + t53 * mrSges(6,2) + t39 * mrSges(7,2) - t52 * mrSges(6,3) - t37 * mrSges(7,3) - pkin(4) * t93 - t188 * t90 + t123 + t124 + (-Ifges(6,4) * t152 - Ifges(5,6) * t149) * t150 + (-t91 + t92) * qJ(5) + m(6) * (-pkin(4) * t53 - qJ(5) * t52) + m(7) * (qJ(5) * t39 - t188 * t37) + t171 * t153 + t175; m(7) * (qJ(5) * t112 - t111 * t188) + t134 + t112 * mrSges(7,2) - t111 * mrSges(7,3) + (-mrSges(6,1) * pkin(4) - mrSges(7,1) * t188 - Ifges(6,4)) * t149 + (qJ(5) * t189 - t202) * t152 + (m(6) * (-pkin(4) * t149 + t180) - t159 - t160) * pkin(10) + t104; -0.2e1 * pkin(4) * mrSges(6,2) + 0.2e1 * t188 * mrSges(7,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJ(5) + m(6) * (pkin(4) ^ 2 + t155) + m(7) * (t188 ^ 2 + t155) - t171; m(6) * t5 + m(7) * t1 + t27 + t30; m(7) * t37 + t198 + t90; m(7) * t111 + (m(6) * pkin(10) + t189) * t149; -m(6) * pkin(4) - m(7) * t188 + mrSges(6,2) - mrSges(7,3); m(6) + m(7); m(7) * t2 + t29; m(7) * t39 + t92; m(7) * t112 + t152 * mrSges(7,1); m(7) * qJ(5) + mrSges(7,2); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
