% Calculate joint inertia matrix for
% S6RRRRPP8
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:22
% EndTime: 2019-03-09 21:30:27
% DurationCPUTime: 2.07s
% Computational Cost: add. (2227->440), mult. (5010->574), div. (0->0), fcn. (5044->8), ass. (0->148)
t145 = sin(qJ(4));
t148 = cos(qJ(4));
t198 = t145 ^ 2 + t148 ^ 2;
t197 = 2 * pkin(9);
t149 = cos(qJ(3));
t139 = t149 * pkin(4);
t190 = pkin(9) * t149;
t123 = t145 * t190;
t146 = sin(qJ(3));
t94 = -pkin(3) * t149 - pkin(10) * t146 - pkin(2);
t57 = t148 * t94 - t123;
t52 = t139 - t57;
t175 = t146 * t148;
t91 = t149 * mrSges(6,1) + mrSges(6,2) * t175;
t196 = m(6) * t52 + t91;
t144 = sin(pkin(6));
t150 = cos(qJ(2));
t177 = t144 * t150;
t147 = sin(qJ(2));
t178 = t144 * t147;
t180 = cos(pkin(6));
t79 = t146 * t180 + t149 * t178;
t46 = t145 * t79 + t148 * t177;
t47 = -t145 * t177 + t148 * t79;
t78 = t146 * t178 - t149 * t180;
t11 = Ifges(5,5) * t47 - Ifges(5,6) * t46 + Ifges(5,3) * t78;
t13 = Ifges(6,4) * t47 + Ifges(6,2) * t78 + Ifges(6,6) * t46;
t9 = Ifges(7,5) * t47 + Ifges(7,6) * t46 - Ifges(7,3) * t78;
t195 = t11 + t13 - t9;
t194 = -2 * mrSges(7,3);
t193 = t144 ^ 2;
t192 = pkin(4) + pkin(5);
t189 = mrSges(6,2) - mrSges(7,3);
t188 = pkin(10) - qJ(6);
t169 = pkin(1) * t180;
t80 = -pkin(8) * t178 + t150 * t169;
t59 = -pkin(2) * t180 - t80;
t24 = t78 * pkin(3) - t79 * pkin(10) + t59;
t81 = pkin(8) * t177 + t147 * t169;
t60 = pkin(9) * t180 + t81;
t61 = (-pkin(2) * t150 - pkin(9) * t147 - pkin(1)) * t144;
t34 = t146 * t61 + t149 * t60;
t26 = -pkin(10) * t177 + t34;
t7 = t145 * t24 + t148 * t26;
t33 = -t146 * t60 + t149 * t61;
t187 = Ifges(5,4) * t145;
t186 = Ifges(5,4) * t148;
t185 = Ifges(7,4) * t145;
t184 = Ifges(7,4) * t148;
t183 = Ifges(6,5) * t145;
t182 = Ifges(6,5) * t148;
t181 = Ifges(7,5) * t148;
t58 = t145 * t94 + t148 * t190;
t179 = qJ(5) * t148;
t176 = t145 * t146;
t174 = Ifges(6,4) * t175 + Ifges(6,6) * t176;
t97 = t148 * mrSges(7,1) + t145 * mrSges(7,2);
t103 = Ifges(5,5) * t145 + Ifges(5,6) * t148;
t173 = Ifges(4,5) * t146 + Ifges(4,6) * t149;
t172 = t198 * pkin(10) ^ 2;
t171 = -Ifges(5,3) - Ifges(7,3) - Ifges(6,2);
t4 = t78 * qJ(5) + t7;
t25 = pkin(3) * t177 - t33;
t170 = Ifges(3,5) * t178 + Ifges(3,6) * t177 + Ifges(3,3) * t180;
t32 = -t78 * mrSges(6,1) + t47 * mrSges(6,2);
t19 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t30 = -t78 * mrSges(7,1) - t47 * mrSges(7,3);
t168 = qJ(5) * t145 + pkin(3);
t6 = -t145 * t26 + t148 * t24;
t10 = Ifges(6,5) * t47 + Ifges(6,6) * t78 + Ifges(6,3) * t46;
t12 = Ifges(7,4) * t47 + Ifges(7,2) * t46 - Ifges(7,6) * t78;
t14 = Ifges(5,4) * t47 - Ifges(5,2) * t46 + Ifges(5,6) * t78;
t166 = t10 / 0.2e1 + t12 / 0.2e1 - t14 / 0.2e1;
t15 = Ifges(7,1) * t47 + Ifges(7,4) * t46 - Ifges(7,5) * t78;
t16 = Ifges(6,1) * t47 + Ifges(6,4) * t78 + Ifges(6,5) * t46;
t17 = Ifges(5,1) * t47 - Ifges(5,4) * t46 + Ifges(5,5) * t78;
t165 = t15 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1;
t63 = -Ifges(6,6) * t149 + (Ifges(6,3) * t145 + t182) * t146;
t65 = Ifges(7,6) * t149 + (Ifges(7,2) * t145 + t184) * t146;
t67 = -Ifges(5,6) * t149 + (-Ifges(5,2) * t145 + t186) * t146;
t164 = t63 / 0.2e1 + t65 / 0.2e1 - t67 / 0.2e1;
t68 = Ifges(7,5) * t149 + (Ifges(7,1) * t148 + t185) * t146;
t69 = -Ifges(6,4) * t149 + (Ifges(6,1) * t148 + t183) * t146;
t70 = -Ifges(5,5) * t149 + (Ifges(5,1) * t148 - t187) * t146;
t163 = t68 / 0.2e1 + t69 / 0.2e1 + t70 / 0.2e1;
t83 = -mrSges(7,1) * t176 + mrSges(7,2) * t175;
t89 = mrSges(7,1) * t149 - mrSges(7,3) * t175;
t51 = -qJ(5) * t149 + t58;
t135 = Ifges(6,4) * t145;
t162 = Ifges(7,5) * t145 / 0.2e1 - t103 / 0.2e1 - t135 / 0.2e1 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t148;
t102 = -Ifges(6,3) * t148 + t183;
t104 = -Ifges(7,2) * t148 + t185;
t106 = Ifges(5,2) * t148 + t187;
t161 = t102 / 0.2e1 + t104 / 0.2e1 - t106 / 0.2e1;
t108 = Ifges(7,1) * t145 - t184;
t109 = Ifges(6,1) * t145 - t182;
t110 = Ifges(5,1) * t145 + t186;
t160 = t108 / 0.2e1 + t109 / 0.2e1 + t110 / 0.2e1;
t159 = t145 * mrSges(5,1) + t148 * mrSges(5,2);
t158 = t145 * mrSges(6,1) - t148 * mrSges(6,3);
t157 = -pkin(4) * t145 + t179;
t156 = Ifges(4,5) * t79 - Ifges(4,6) * t78 - Ifges(4,3) * t177;
t155 = qJ(5) * t47 - t25;
t154 = pkin(9) ^ 2;
t152 = qJ(5) ^ 2;
t143 = t149 ^ 2;
t141 = t146 ^ 2;
t137 = t141 * t154;
t121 = Ifges(5,5) * t175;
t111 = Ifges(4,1) * t146 + Ifges(4,4) * t149;
t107 = Ifges(4,4) * t146 + Ifges(4,2) * t149;
t100 = t188 * t148;
t99 = -mrSges(4,1) * t149 + mrSges(4,2) * t146;
t98 = -mrSges(5,1) * t148 + mrSges(5,2) * t145;
t96 = -mrSges(6,1) * t148 - mrSges(6,3) * t145;
t95 = t188 * t145;
t93 = -pkin(4) * t148 - t168;
t92 = -mrSges(6,2) * t176 - mrSges(6,3) * t149;
t90 = -mrSges(5,1) * t149 - mrSges(5,3) * t175;
t88 = mrSges(5,2) * t149 - mrSges(5,3) * t176;
t87 = -mrSges(7,2) * t149 + mrSges(7,3) * t176;
t85 = t148 * t192 + t168;
t84 = t159 * t146;
t82 = t158 * t146;
t71 = (pkin(9) - t157) * t146;
t66 = -Ifges(6,2) * t149 + t174;
t64 = -Ifges(5,6) * t176 - Ifges(5,3) * t149 + t121;
t62 = Ifges(7,3) * t149 + (Ifges(7,6) * t145 + t181) * t146;
t50 = -mrSges(4,1) * t177 - mrSges(4,3) * t79;
t49 = mrSges(4,2) * t177 - mrSges(4,3) * t78;
t48 = (-t145 * t192 - pkin(9) + t179) * t146;
t39 = qJ(6) * t176 + t51;
t38 = mrSges(4,1) * t78 + mrSges(4,2) * t79;
t37 = pkin(5) * t149 + t123 + t139 + (-qJ(6) * t146 - t94) * t148;
t36 = Ifges(4,1) * t79 - Ifges(4,4) * t78 - Ifges(4,5) * t177;
t35 = Ifges(4,4) * t79 - Ifges(4,2) * t78 - Ifges(4,6) * t177;
t31 = mrSges(5,1) * t78 - mrSges(5,3) * t47;
t29 = -mrSges(5,2) * t78 - mrSges(5,3) * t46;
t28 = mrSges(7,2) * t78 + mrSges(7,3) * t46;
t27 = -mrSges(6,2) * t46 + mrSges(6,3) * t78;
t20 = mrSges(5,1) * t46 + mrSges(5,2) * t47;
t18 = mrSges(6,1) * t46 - mrSges(6,3) * t47;
t8 = pkin(4) * t46 - t155;
t5 = -pkin(4) * t78 - t6;
t3 = -t192 * t46 + t155;
t2 = qJ(6) * t46 + t4;
t1 = -qJ(6) * t47 - t192 * t78 - t6;
t21 = [m(4) * (t33 ^ 2 + t34 ^ 2 + t59 ^ 2) + m(5) * (t25 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (t15 + t16 + t17) * t47 + (t10 + t12 - t14) * t46 + t79 * t36 + 0.2e1 * t59 * t38 + 0.2e1 * t34 * t49 + 0.2e1 * t33 * t50 + 0.2e1 * t25 * t20 + 0.2e1 * t4 * t27 + 0.2e1 * t2 * t28 + 0.2e1 * t7 * t29 + 0.2e1 * t1 * t30 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + 0.2e1 * t8 * t18 + 0.2e1 * t3 * t19 + (-t35 + t195) * t78 - 0.2e1 * t193 * pkin(1) * (-mrSges(3,1) * t150 + mrSges(3,2) * t147) + Ifges(2,3) - t156 * t177 + (Ifges(3,6) * t180 + (Ifges(3,4) * t147 + Ifges(3,2) * t150) * t144) * t177 + (Ifges(3,5) * t180 + (Ifges(3,1) * t147 + Ifges(3,4) * t150) * t144) * t178 + t180 * t170 + 0.2e1 * t81 * (-mrSges(3,2) * t180 + mrSges(3,3) * t177) + 0.2e1 * t80 * (mrSges(3,1) * t180 - mrSges(3,3) * t178) + m(3) * (pkin(1) ^ 2 * t193 + t80 ^ 2 + t81 ^ 2); t59 * t99 + t79 * t111 / 0.2e1 + t8 * t82 + t3 * t83 + t25 * t84 + t2 * t87 + t7 * t88 + t1 * t89 + t6 * t90 + t5 * t91 + t4 * t92 + t80 * mrSges(3,1) - t81 * mrSges(3,2) + t170 + (pkin(9) * t49 + t34 * mrSges(4,3) + t35 / 0.2e1 + t9 / 0.2e1 - t11 / 0.2e1 - t13 / 0.2e1) * t149 + m(5) * (pkin(9) * t146 * t25 + t57 * t6 + t58 * t7) + m(6) * (t4 * t51 + t5 * t52 + t71 * t8) + m(7) * (t1 * t37 + t2 * t39 + t3 * t48) + (-t107 / 0.2e1 + t64 / 0.2e1 + t66 / 0.2e1 - t62 / 0.2e1) * t78 + m(4) * (-pkin(2) * t59 + (-t146 * t33 + t149 * t34) * pkin(9)) + t58 * t29 + t71 * t18 + t48 * t19 + t51 * t27 + t52 * t32 + t57 * t31 + t37 * t30 - pkin(2) * t38 + t39 * t28 + t163 * t47 + t164 * t46 + (t36 / 0.2e1 - t33 * mrSges(4,3) + (-t50 + t20) * pkin(9) + t165 * t148 + t166 * t145) * t146 - t173 * t177 / 0.2e1; -0.2e1 * pkin(2) * t99 + 0.2e1 * t37 * t89 + 0.2e1 * t39 * t87 + 0.2e1 * t48 * t83 + 0.2e1 * t51 * t92 + 0.2e1 * t52 * t91 + 0.2e1 * t57 * t90 + 0.2e1 * t58 * t88 + 0.2e1 * t71 * t82 + Ifges(3,3) + (t141 + t143) * mrSges(4,3) * t197 + (-t66 + t107 + t62 - t64) * t149 + m(4) * (pkin(2) ^ 2 + t143 * t154 + t137) + m(5) * (t57 ^ 2 + t58 ^ 2 + t137) + m(6) * (t51 ^ 2 + t52 ^ 2 + t71 ^ 2) + m(7) * (t37 ^ 2 + t39 ^ 2 + t48 ^ 2) + (t84 * t197 + t111 + (t68 + t69 + t70) * t148 + (t63 + t65 - t67) * t145) * t146; m(7) * (t1 * t95 + t100 * t2 + t3 * t85) - t162 * t78 + t160 * t47 + t161 * t46 + t95 * t30 + t8 * t96 + t3 * t97 + t25 * t98 + t100 * t28 + t85 * t19 + t93 * t18 + m(6) * (t8 * t93 + (t145 * t5 + t148 * t4) * pkin(10)) + m(5) * (-pkin(3) * t25 + (-t145 * t6 + t148 * t7) * pkin(10)) + t33 * mrSges(4,1) - t34 * mrSges(4,2) - pkin(3) * t20 + t156 + (t4 * mrSges(6,2) + t7 * mrSges(5,3) - t2 * mrSges(7,3) + (t29 + t27) * pkin(10) - t166) * t148 + (t5 * mrSges(6,2) - t6 * mrSges(5,3) - t1 * mrSges(7,3) + (-t31 + t32) * pkin(10) + t165) * t145; -pkin(3) * t84 + t100 * t87 + t48 * t97 + t93 * t82 + t85 * t83 + t95 * t89 + m(7) * (t100 * t39 + t37 * t95 + t48 * t85) + t162 * t149 + (-t149 * mrSges(4,2) + (-m(5) * pkin(3) - mrSges(4,1) + t98) * t146) * pkin(9) + (t58 * mrSges(5,3) + t51 * mrSges(6,2) - t39 * mrSges(7,3) + t160 * t146 + (m(5) * t58 + m(6) * t51 + t88 + t92) * pkin(10) - t164) * t148 + (t52 * mrSges(6,2) - t37 * mrSges(7,3) - t57 * mrSges(5,3) + t161 * t146 + (-m(5) * t57 + t196 - t90) * pkin(10) + t163) * t145 + t173 + (m(6) * t93 + t96) * t71; -0.2e1 * pkin(3) * t98 + 0.2e1 * t85 * t97 + 0.2e1 * t93 * t96 + Ifges(4,3) + m(6) * (t93 ^ 2 + t172) + m(7) * (t100 ^ 2 + t85 ^ 2 + t95 ^ 2) + m(5) * (pkin(3) ^ 2 + t172) + (t100 * t194 - t102 - t104 + t106) * t148 + (t194 * t95 + t108 + t109 + t110) * t145 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * pkin(10) * t198; t6 * mrSges(5,1) - t5 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3) - pkin(4) * t32 - t192 * t30 + (t27 + t28) * qJ(5) + m(6) * (-pkin(4) * t5 + qJ(5) * t4) + m(7) * (qJ(5) * t2 - t1 * t192) + t195; t57 * mrSges(5,1) - t52 * mrSges(6,1) - t37 * mrSges(7,1) - t58 * mrSges(5,2) + t39 * mrSges(7,2) + t51 * mrSges(6,3) - pkin(4) * t91 - t192 * t89 + t121 + (t87 + t92) * qJ(5) + m(7) * (qJ(5) * t39 - t192 * t37) + m(6) * (-pkin(4) * t52 + qJ(5) * t51) + t171 * t149 + (-t181 + (-Ifges(5,6) - Ifges(7,6)) * t145) * t146 + t174; t135 + m(7) * (qJ(5) * t100 - t192 * t95) - t95 * mrSges(7,1) + t100 * mrSges(7,2) + (-mrSges(6,2) * pkin(4) + mrSges(7,3) * t192 - Ifges(7,5)) * t145 + (qJ(5) * t189 - Ifges(6,6) + Ifges(7,6)) * t148 + (m(6) * t157 - t158 - t159) * pkin(10) + t103; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * t192 * mrSges(7,1) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJ(5) + m(6) * (pkin(4) ^ 2 + t152) + m(7) * (t192 ^ 2 + t152) - t171; m(6) * t5 + m(7) * t1 + t30 + t32; m(7) * t37 + t196 + t89; m(7) * t95 + (m(6) * pkin(10) + t189) * t145; -m(6) * pkin(4) - m(7) * t192 - mrSges(6,1) - mrSges(7,1); m(6) + m(7); m(7) * t3 + t19; m(7) * t48 + t83; m(7) * t85 + t97; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
