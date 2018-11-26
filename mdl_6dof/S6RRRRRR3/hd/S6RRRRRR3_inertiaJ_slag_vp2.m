% Calculate joint inertia matrix for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:38:49
% EndTime: 2018-11-23 18:38:51
% DurationCPUTime: 1.54s
% Computational Cost: add. (4211->356), mult. (7974->500), div. (0->0), fcn. (9028->10), ass. (0->145)
t157 = sin(qJ(4));
t162 = cos(qJ(4));
t164 = cos(qJ(2));
t218 = -pkin(8) - pkin(7);
t139 = t218 * t164;
t158 = sin(qJ(3));
t163 = cos(qJ(3));
t159 = sin(qJ(2));
t185 = t218 * t159;
t100 = -t163 * t139 + t158 * t185;
t126 = t158 * t159 - t163 * t164;
t128 = t158 * t164 + t159 * t163;
t146 = -pkin(2) * t164 - pkin(1);
t83 = pkin(3) * t126 - pkin(9) * t128 + t146;
t52 = -t100 * t157 + t162 * t83;
t53 = t162 * t100 + t157 * t83;
t181 = -t157 * t52 + t162 * t53;
t156 = sin(qJ(5));
t161 = cos(qJ(5));
t125 = -t156 * t157 + t161 * t162;
t127 = t156 * t162 + t157 * t161;
t195 = Ifges(6,5) * t127 + Ifges(6,6) * t125;
t155 = sin(qJ(6));
t160 = cos(qJ(6));
t84 = t125 * t160 - t127 * t155;
t85 = t125 * t155 + t127 * t160;
t206 = Ifges(7,5) * t85 + Ifges(7,6) * t84;
t231 = Ifges(5,5) * t157 + Ifges(5,6) * t162 + t195 + t206;
t210 = t84 * mrSges(7,3);
t230 = t155 * pkin(5) * t210 + t195;
t214 = pkin(11) * t127;
t142 = pkin(2) * t158 + pkin(9);
t121 = (-pkin(10) - t142) * t157;
t150 = t162 * pkin(10);
t122 = t142 * t162 + t150;
t78 = t161 * t121 - t122 * t156;
t59 = t78 - t214;
t120 = t125 * pkin(11);
t79 = t156 * t121 + t161 * t122;
t60 = t120 + t79;
t22 = -t155 * t60 + t160 * t59;
t23 = t155 * t59 + t160 * t60;
t229 = t22 * mrSges(7,1) - t23 * mrSges(7,2);
t137 = (-pkin(10) - pkin(9)) * t157;
t138 = pkin(9) * t162 + t150;
t96 = t161 * t137 - t138 * t156;
t67 = t96 - t214;
t99 = t156 * t137 + t161 * t138;
t68 = t120 + t99;
t35 = -t155 * t68 + t160 * t67;
t36 = t155 * t67 + t160 * t68;
t228 = t35 * mrSges(7,1) - t36 * mrSges(7,2);
t227 = t78 * mrSges(6,1) - t79 * mrSges(6,2);
t226 = t96 * mrSges(6,1) - t99 * mrSges(6,2);
t196 = t128 * t162;
t225 = Ifges(5,5) * t196 + Ifges(5,3) * t126;
t71 = t127 * t128;
t72 = t125 * t128;
t224 = Ifges(6,5) * t72 - Ifges(6,6) * t71 + Ifges(6,3) * t126;
t97 = -t139 * t158 - t163 * t185;
t223 = t97 ^ 2;
t47 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t222 = 0.2e1 * t47;
t88 = -mrSges(6,1) * t125 + mrSges(6,2) * t127;
t221 = 0.2e1 * t88;
t220 = 0.2e1 * t97;
t219 = 0.2e1 * t146;
t216 = pkin(2) * t163;
t215 = pkin(4) * t156;
t209 = t85 * mrSges(7,3);
t207 = Ifges(6,3) + Ifges(7,3);
t30 = pkin(4) * t126 - pkin(10) * t196 + t52;
t197 = t128 * t157;
t40 = -pkin(10) * t197 + t53;
t12 = t156 * t30 + t161 * t40;
t205 = Ifges(5,4) * t157;
t204 = Ifges(5,4) * t162;
t143 = pkin(4) * t161 + pkin(5);
t112 = t143 * t155 + t160 * t215;
t203 = t112 * mrSges(7,2);
t202 = t125 * mrSges(6,3);
t201 = t127 * mrSges(6,3);
t200 = t155 * mrSges(7,2);
t193 = t157 ^ 2 + t162 ^ 2;
t192 = t159 ^ 2 + t164 ^ 2;
t191 = 0.2e1 * mrSges(6,3);
t190 = 0.2e1 * mrSges(7,3);
t189 = pkin(5) * t200;
t188 = t160 * t209;
t41 = -t155 * t72 - t160 * t71;
t42 = -t155 * t71 + t160 * t72;
t187 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t126;
t186 = t161 * t201;
t145 = -pkin(4) * t162 - pkin(3);
t11 = -t156 * t40 + t161 * t30;
t184 = t193 * t142;
t111 = t143 * t160 - t155 * t215;
t104 = t111 * mrSges(7,1);
t183 = Ifges(7,3) + t104 - t203;
t66 = pkin(4) * t197 + t97;
t182 = mrSges(5,1) * t157 + mrSges(5,2) * t162;
t86 = -mrSges(5,2) * t126 - mrSges(5,3) * t197;
t87 = mrSges(5,1) * t126 - mrSges(5,3) * t196;
t180 = -t157 * t87 + t162 * t86;
t179 = t206 + t229;
t178 = t206 + t228;
t177 = 0.2e1 * t193 * mrSges(5,3);
t103 = -pkin(5) * t125 + t145;
t5 = pkin(5) * t126 - pkin(11) * t72 + t11;
t6 = -pkin(11) * t71 + t12;
t3 = -t155 * t6 + t160 * t5;
t4 = t155 * t5 + t160 * t6;
t176 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t187;
t175 = (mrSges(4,1) * t163 - mrSges(4,2) * t158) * pkin(2);
t174 = (mrSges(6,1) * t161 - mrSges(6,2) * t156) * pkin(4);
t135 = Ifges(5,2) * t162 + t205;
t136 = Ifges(5,1) * t157 + t204;
t48 = Ifges(7,4) * t85 + Ifges(7,2) * t84;
t49 = Ifges(7,1) * t85 + Ifges(7,4) * t84;
t89 = Ifges(6,4) * t127 + Ifges(6,2) * t125;
t90 = Ifges(6,1) * t127 + Ifges(6,4) * t125;
t173 = t125 * t89 + t127 * t90 + t162 * t135 + t157 * t136 + t84 * t48 + t85 * t49 + Ifges(4,3);
t172 = -t111 * t209 + t112 * t210 + t202 * t215 + t231;
t171 = t11 * mrSges(6,1) - t12 * mrSges(6,2) + t176 + t224;
t13 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t126;
t134 = -mrSges(5,1) * t162 + mrSges(5,2) * t157;
t14 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t126;
t31 = Ifges(6,4) * t72 - Ifges(6,2) * t71 + Ifges(6,6) * t126;
t32 = Ifges(6,1) * t72 - Ifges(6,4) * t71 + Ifges(6,5) * t126;
t43 = pkin(5) * t71 + t66;
t61 = Ifges(5,6) * t126 + (-Ifges(5,2) * t157 + t204) * t128;
t62 = Ifges(5,5) * t126 + (Ifges(5,1) * t162 - t205) * t128;
t170 = -t100 * mrSges(4,2) + t84 * t13 / 0.2e1 + t85 * t14 / 0.2e1 + t66 * t88 - t71 * t89 / 0.2e1 + t72 * t90 / 0.2e1 + t43 * t47 + t41 * t48 / 0.2e1 + t42 * t49 / 0.2e1 + t12 * t202 + t4 * t210 - t3 * t209 - t11 * t201 - t135 * t197 / 0.2e1 + t136 * t196 / 0.2e1 + t125 * t31 / 0.2e1 + t127 * t32 / 0.2e1 + Ifges(4,5) * t128 + t157 * t62 / 0.2e1 + t162 * t61 / 0.2e1 + (-mrSges(4,1) + t134) * t97 + t181 * mrSges(5,3) + (-Ifges(4,6) + t231 / 0.2e1) * t126;
t147 = t160 * pkin(5) * mrSges(7,1);
t144 = -pkin(3) - t216;
t133 = t145 - t216;
t101 = t103 - t216;
t80 = t182 * t128;
t55 = mrSges(6,1) * t126 - mrSges(6,3) * t72;
t54 = -mrSges(6,2) * t126 - mrSges(6,3) * t71;
t44 = mrSges(6,1) * t71 + mrSges(6,2) * t72;
t29 = mrSges(7,1) * t126 - mrSges(7,3) * t42;
t28 = -mrSges(7,2) * t126 + mrSges(7,3) * t41;
t15 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t1 = [0.2e1 * t53 * t86 + 0.2e1 * t52 * t87 - t71 * t31 + t72 * t32 + 0.2e1 * t12 * t54 + 0.2e1 * t11 * t55 + 0.2e1 * t66 * t44 + t41 * t13 + t42 * t14 + 0.2e1 * t43 * t15 + 0.2e1 * t4 * t28 + 0.2e1 * t3 * t29 + m(3) * (t192 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2 + t66 ^ 2) + 0.2e1 * t192 * pkin(7) * mrSges(3,3) + m(7) * (t3 ^ 2 + t4 ^ 2 + t43 ^ 2) + t80 * t220 + (mrSges(4,2) * t219 + mrSges(4,3) * t220 + Ifges(4,1) * t128 - t157 * t61 + t162 * t62) * t128 + m(5) * (t52 ^ 2 + t53 ^ 2 + t223) + m(4) * (t100 ^ 2 + t146 ^ 2 + t223) + (mrSges(4,1) * t219 - 0.2e1 * t100 * mrSges(4,3) + Ifges(4,2) * t126 + (-Ifges(5,6) * t157 - (2 * Ifges(4,4))) * t128 + t187 + t224 + t225) * t126 + Ifges(2,3) + t159 * (Ifges(3,1) * t159 + Ifges(3,4) * t164) + t164 * (Ifges(3,4) * t159 + Ifges(3,2) * t164) - 0.2e1 * pkin(1) * (-t164 * mrSges(3,1) + t159 * mrSges(3,2)); t180 * t142 + m(7) * (t101 * t43 + t22 * t3 + t23 * t4) + m(6) * (t11 * t78 + t12 * t79 + t133 * t66) + t101 * t15 + t79 * t54 + t78 * t55 + t23 * t28 + t22 * t29 + t170 + (-mrSges(3,1) * t159 - mrSges(3,2) * t164) * pkin(7) + t133 * t44 + t144 * t80 + Ifges(3,5) * t159 + Ifges(3,6) * t164 + m(5) * (t181 * t142 + t144 * t97) + (m(4) * (t100 * t158 - t163 * t97) + (-t126 * t158 - t128 * t163) * mrSges(4,3)) * pkin(2); t101 * t222 + t173 + t142 * t177 + m(4) * (t158 ^ 2 + t163 ^ 2) * pkin(2) ^ 2 + m(7) * (t101 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(6) * (t133 ^ 2 + t78 ^ 2 + t79 ^ 2) + 0.2e1 * t175 + (t79 * t125 - t78 * t127) * t191 + (-t22 * t85 + t23 * t84) * t190 + t133 * t221 + 0.2e1 * t144 * t134 + Ifges(3,3) + m(5) * (t193 * t142 ^ 2 + t144 ^ 2); t96 * t55 + t99 * t54 + t103 * t15 - pkin(3) * t80 + t35 * t29 + t36 * t28 + t170 + t180 * pkin(9) + m(7) * (t103 * t43 + t3 * t35 + t36 * t4) + m(6) * (t11 * t96 + t12 * t99 + t145 * t66) + t145 * t44 + m(5) * (-pkin(3) * t97 + t181 * pkin(9)); t173 + m(7) * (t101 * t103 + t22 * t35 + t23 * t36) + m(6) * (t133 * t145 + t78 * t96 + t79 * t99) + (t133 + t145) * t88 + (t101 + t103) * t47 + (-pkin(3) + t144) * t134 + t175 + (t193 * pkin(9) + t184) * mrSges(5,3) + ((-t78 - t96) * t127 + (t79 + t99) * t125) * mrSges(6,3) + ((-t22 - t35) * t85 + (t23 + t36) * t84) * mrSges(7,3) + m(5) * (-pkin(3) * t144 + pkin(9) * t184); -0.2e1 * pkin(3) * t134 + t103 * t222 + t145 * t221 + (-t35 * t85 + t36 * t84) * t190 + (t99 * t125 - t96 * t127) * t191 + pkin(9) * t177 + m(7) * (t103 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t145 ^ 2 + t96 ^ 2 + t99 ^ 2) + m(5) * (t193 * pkin(9) ^ 2 + pkin(3) ^ 2) + t173; t111 * t29 + t112 * t28 + t52 * mrSges(5,1) - t53 * mrSges(5,2) + t171 + m(7) * (t111 * t3 + t112 * t4) - Ifges(5,6) * t197 + (m(6) * (t11 * t161 + t12 * t156) + t156 * t54 + t161 * t55) * pkin(4) + t225; t172 + (m(6) * (t156 * t79 + t161 * t78) - t186) * pkin(4) - t182 * t142 + m(7) * (t111 * t22 + t112 * t23) + t227 + t229; t172 - t182 * pkin(9) + m(7) * (t111 * t35 + t112 * t36) + (m(6) * (t156 * t99 + t161 * t96) - t186) * pkin(4) + t226 + t228; -0.2e1 * t203 + Ifges(5,3) + 0.2e1 * t104 + 0.2e1 * t174 + m(7) * (t111 ^ 2 + t112 ^ 2) + m(6) * (t156 ^ 2 + t161 ^ 2) * pkin(4) ^ 2 + t207; (m(7) * (t155 * t4 + t160 * t3) + t155 * t28 + t160 * t29) * pkin(5) + t171; (m(7) * (t155 * t23 + t160 * t22) - t188) * pkin(5) + t179 + t227 + t230; (m(7) * (t155 * t36 + t160 * t35) - t188) * pkin(5) + t178 + t226 + t230; Ifges(6,3) + t147 + t174 + (m(7) * (t111 * t160 + t112 * t155) - t200) * pkin(5) + t183; -0.2e1 * t189 + 0.2e1 * t147 + m(7) * (t155 ^ 2 + t160 ^ 2) * pkin(5) ^ 2 + t207; t176; t179; t178; t183; Ifges(7,3) + t147 - t189; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
