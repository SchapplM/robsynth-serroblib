% Calculate joint inertia matrix for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:20:47
% EndTime: 2018-11-23 16:20:49
% DurationCPUTime: 1.94s
% Computational Cost: add. (5951->375), mult. (15609->565), div. (0->0), fcn. (18065->14), ass. (0->142)
t194 = Ifges(5,3) + Ifges(6,3);
t136 = sin(pkin(13));
t140 = cos(pkin(13));
t145 = sin(qJ(4));
t148 = cos(qJ(4));
t104 = t136 * t145 - t140 * t148;
t105 = t136 * t148 + t140 * t145;
t193 = Ifges(5,5) * t145 + Ifges(6,5) * t105 + Ifges(5,6) * t148 - Ifges(6,6) * t104;
t138 = sin(pkin(7));
t142 = cos(pkin(7));
t146 = sin(qJ(3));
t149 = cos(qJ(3));
t143 = cos(pkin(6));
t139 = sin(pkin(6));
t141 = cos(pkin(12));
t167 = t139 * t141;
t160 = t142 * t167;
t137 = sin(pkin(12));
t180 = pkin(1) * t143;
t97 = qJ(2) * t167 + t137 * t180;
t73 = (t138 * t143 + t160) * pkin(9) + t97;
t120 = t141 * t180;
t170 = t137 * t139;
t82 = pkin(2) * t143 + t120 + (-pkin(9) * t142 - qJ(2)) * t170;
t91 = (-pkin(9) * t137 * t138 - pkin(2) * t141 - pkin(1)) * t139;
t45 = -t146 * t73 + (t138 * t91 + t142 * t82) * t149;
t169 = t138 * t146;
t98 = t142 * t148 - t145 * t169;
t99 = t142 * t145 + t148 * t169;
t65 = t136 * t99 - t140 * t98;
t192 = t65 ^ 2;
t178 = -qJ(5) - pkin(10);
t110 = t178 * t148;
t159 = t178 * t145;
t87 = -t110 * t136 - t140 * t159;
t191 = t87 ^ 2;
t190 = 2 * mrSges(3,1);
t189 = 0.2e1 * t87;
t188 = 0.2e1 * t143;
t144 = sin(qJ(6));
t147 = cos(qJ(6));
t166 = t142 * t146;
t81 = t143 * t169 + (t137 * t149 + t141 * t166) * t139;
t95 = -t138 * t167 + t142 * t143;
t56 = -t145 * t81 + t148 * t95;
t57 = t145 * t95 + t148 * t81;
t43 = t136 * t56 + t140 * t57;
t168 = t138 * t149;
t80 = -t143 * t168 + t146 * t170 - t149 * t160;
t26 = -t144 * t43 + t147 * t80;
t187 = t26 / 0.2e1;
t27 = t144 * t80 + t147 * t43;
t186 = t27 / 0.2e1;
t111 = Ifges(7,5) * t144 + Ifges(7,6) * t147;
t184 = t111 / 0.2e1;
t183 = -t144 / 0.2e1;
t182 = t144 / 0.2e1;
t181 = t147 / 0.2e1;
t179 = t65 * t87;
t13 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t29 = mrSges(6,1) * t80 - mrSges(6,3) * t43;
t177 = t13 - t29;
t52 = -t138 * t82 + t142 * t91;
t34 = pkin(3) * t80 - pkin(10) * t81 + t52;
t46 = t149 * t73 + t82 * t166 + t91 * t169;
t37 = pkin(10) * t95 + t46;
t20 = -t145 * t37 + t148 * t34;
t12 = pkin(4) * t80 - qJ(5) * t57 + t20;
t21 = t145 * t34 + t148 * t37;
t15 = qJ(5) * t56 + t21;
t6 = t136 * t12 + t140 * t15;
t176 = Ifges(7,4) * t144;
t175 = Ifges(7,4) * t147;
t174 = t105 * t144;
t173 = t105 * t147;
t124 = pkin(4) * t136 + pkin(11);
t172 = t124 * t144;
t171 = t124 * t147;
t163 = t144 ^ 2 + t147 ^ 2;
t162 = t145 ^ 2 + t148 ^ 2;
t42 = t136 * t57 - t140 * t56;
t7 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t42;
t161 = Ifges(4,5) * t81 - Ifges(4,6) * t80 + Ifges(4,3) * t95;
t126 = -pkin(4) * t148 - pkin(3);
t22 = t42 * mrSges(6,1) + t43 * mrSges(6,2);
t83 = t104 * mrSges(6,1) + t105 * mrSges(6,2);
t36 = -pkin(3) * t95 - t45;
t23 = -pkin(4) * t56 + t36;
t10 = pkin(5) * t42 - pkin(11) * t43 + t23;
t4 = pkin(11) * t80 + t6;
t1 = t10 * t147 - t144 * t4;
t2 = t10 * t144 + t147 * t4;
t158 = -t1 * t144 + t2 * t147;
t157 = mrSges(7,1) * t144 + mrSges(7,2) * t147;
t5 = t12 * t140 - t136 * t15;
t72 = pkin(5) * t104 - pkin(11) * t105 + t126;
t89 = -t140 * t110 + t136 * t159;
t49 = -t144 * t89 + t147 * t72;
t50 = t144 * t72 + t147 * t89;
t155 = -t49 * t144 + t50 * t147;
t67 = t136 * t98 + t140 * t99;
t58 = -t144 * t67 - t147 * t168;
t59 = -t144 * t168 + t147 * t67;
t154 = -t144 * t58 + t147 * t59;
t153 = -t145 * t98 + t148 * t99;
t152 = Ifges(5,5) * t57 + Ifges(6,5) * t43 + Ifges(5,6) * t56 - Ifges(6,6) * t42 + t194 * t80;
t62 = Ifges(7,5) * t173 - Ifges(7,6) * t174 + Ifges(7,3) * t104;
t131 = t138 ^ 2;
t125 = -pkin(4) * t140 - pkin(5);
t122 = t131 * t149 ^ 2;
t118 = mrSges(3,2) * t170;
t115 = Ifges(5,1) * t145 + Ifges(5,4) * t148;
t114 = Ifges(7,1) * t144 + t175;
t113 = Ifges(5,4) * t145 + Ifges(5,2) * t148;
t112 = Ifges(7,2) * t147 + t176;
t109 = -mrSges(5,1) * t148 + mrSges(5,2) * t145;
t108 = -mrSges(7,1) * t147 + mrSges(7,2) * t144;
t96 = -qJ(2) * t170 + t120;
t85 = Ifges(6,1) * t105 - Ifges(6,4) * t104;
t84 = Ifges(6,4) * t105 - Ifges(6,2) * t104;
t75 = mrSges(7,1) * t104 - mrSges(7,3) * t173;
t74 = -mrSges(7,2) * t104 - mrSges(7,3) * t174;
t71 = t157 * t105;
t64 = Ifges(7,5) * t104 + (Ifges(7,1) * t147 - t176) * t105;
t63 = Ifges(7,6) * t104 + (-Ifges(7,2) * t144 + t175) * t105;
t61 = mrSges(4,1) * t95 - mrSges(4,3) * t81;
t60 = -mrSges(4,2) * t95 - mrSges(4,3) * t80;
t51 = mrSges(4,1) * t80 + mrSges(4,2) * t81;
t48 = mrSges(5,1) * t80 - mrSges(5,3) * t57;
t47 = -mrSges(5,2) * t80 + mrSges(5,3) * t56;
t44 = -mrSges(5,1) * t56 + mrSges(5,2) * t57;
t31 = Ifges(5,1) * t57 + Ifges(5,4) * t56 + Ifges(5,5) * t80;
t30 = Ifges(5,4) * t57 + Ifges(5,2) * t56 + Ifges(5,6) * t80;
t28 = -mrSges(6,2) * t80 - mrSges(6,3) * t42;
t19 = Ifges(6,1) * t43 - Ifges(6,4) * t42 + Ifges(6,5) * t80;
t18 = Ifges(6,4) * t43 - Ifges(6,2) * t42 + Ifges(6,6) * t80;
t17 = mrSges(7,1) * t42 - mrSges(7,3) * t27;
t16 = -mrSges(7,2) * t42 + mrSges(7,3) * t26;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t42;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t42;
t3 = -pkin(5) * t80 - t5;
t11 = [(t7 - t18) * t42 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t23 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2 + t36 ^ 2) + m(4) * (t45 ^ 2 + t46 ^ 2 + t52 ^ 2) + m(3) * (pkin(1) ^ 2 * t139 ^ 2 + t96 ^ 2 + t97 ^ 2) + t81 * (Ifges(4,1) * t81 + Ifges(4,5) * t95) + (-0.2e1 * Ifges(4,4) * t81 + Ifges(4,2) * t80 - Ifges(4,6) * t95 + t152) * t80 + (-0.2e1 * t97 * mrSges(3,2) + Ifges(3,3) * t143 + t96 * t190) * t143 + (-0.2e1 * pkin(1) * t118 + (-0.2e1 * t96 * mrSges(3,3) + Ifges(3,1) * t170 + Ifges(3,5) * t188) * t137 + (0.2e1 * t97 * mrSges(3,3) + Ifges(3,6) * t188 + (0.2e1 * Ifges(3,4) * t137 + Ifges(3,2) * t141 + pkin(1) * t190) * t139) * t141) * t139 + t95 * t161 + Ifges(2,3) + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t16 + 0.2e1 * t1 * t17 + 0.2e1 * t23 * t22 + t26 * t8 + t27 * t9 + 0.2e1 * t6 * t28 + 0.2e1 * t5 * t29 + t43 * t19 + 0.2e1 * t36 * t44 + 0.2e1 * t21 * t47 + 0.2e1 * t20 * t48 + 0.2e1 * t52 * t51 + t56 * t30 + t57 * t31 + 0.2e1 * t46 * t60 + 0.2e1 * t45 * t61; t142 * t51 + t59 * t16 + t58 * t17 + t67 * t28 + t99 * t47 + t98 * t48 + t118 + t177 * t65 + (-m(3) * pkin(1) - mrSges(3,1) * t141) * t139 + (t146 * t60 + (-t22 - t44 + t61) * t149) * t138 + m(7) * (t1 * t58 + t2 * t59 + t3 * t65) + m(6) * (-t23 * t168 - t5 * t65 + t6 * t67) + m(5) * (-t36 * t168 + t20 * t98 + t21 * t99) + m(4) * (t142 * t52 + (t146 * t46 + t149 * t45) * t138); m(3) + m(6) * (t67 ^ 2 + t122 + t192) + m(7) * (t58 ^ 2 + t59 ^ 2 + t192) + m(5) * (t98 ^ 2 + t99 ^ 2 + t122) + m(4) * (t131 * t146 ^ 2 + t142 ^ 2 + t122); t126 * t22 + (t30 / 0.2e1 + t21 * mrSges(5,3) + pkin(10) * t47) * t148 + (t31 / 0.2e1 - t20 * mrSges(5,3) - pkin(10) * t48) * t145 + m(5) * (-pkin(3) * t36 + (-t20 * t145 + t21 * t148) * pkin(10)) + t193 * t80 / 0.2e1 + m(7) * (t1 * t49 + t2 * t50 + t3 * t87) + m(6) * (t126 * t23 - t5 * t87 + t6 * t89) + (t62 / 0.2e1 - t84 / 0.2e1) * t42 + t161 + (-t6 * mrSges(6,3) + t7 / 0.2e1 - t18 / 0.2e1) * t104 + (-t5 * mrSges(6,3) + t8 * t183 + t9 * t181 + t19 / 0.2e1) * t105 + t177 * t87 + t64 * t186 + t63 * t187 - pkin(3) * t44 + t45 * mrSges(4,1) - t46 * mrSges(4,2) + t49 * t17 + t50 * t16 + t3 * t71 + t2 * t74 + t1 * t75 + t23 * t83 + t43 * t85 / 0.2e1 + t89 * t28 + t36 * t109 + t56 * t113 / 0.2e1 + t57 * t115 / 0.2e1; t58 * t75 + t59 * t74 + t65 * t71 + (-t104 * t67 + t105 * t65) * mrSges(6,3) + t153 * mrSges(5,3) + (-t146 * mrSges(4,2) + (mrSges(4,1) - t109 - t83) * t149) * t138 + m(6) * (-t126 * t168 + t67 * t89 + t179) + m(7) * (t49 * t58 + t50 * t59 + t179) + m(5) * (pkin(3) * t168 + t153 * pkin(10)); -0.2e1 * pkin(3) * t109 + t148 * t113 + t145 * t115 + 0.2e1 * t126 * t83 + 0.2e1 * t49 * t75 + 0.2e1 * t50 * t74 + t71 * t189 + Ifges(4,3) + 0.2e1 * t162 * pkin(10) * mrSges(5,3) + (-0.2e1 * t89 * mrSges(6,3) + t62 - t84) * t104 + m(7) * (t49 ^ 2 + t50 ^ 2 + t191) + m(6) * (t126 ^ 2 + t89 ^ 2 + t191) + m(5) * (t162 * pkin(10) ^ 2 + pkin(3) ^ 2) + (mrSges(6,3) * t189 - t144 * t63 + t147 * t64 + t85) * t105; t9 * t182 + t8 * t181 + t125 * t13 + t152 + t158 * mrSges(7,3) + m(7) * (t158 * t124 + t125 * t3) + (t136 * t28 + t140 * t29 + m(6) * (t136 * t6 + t140 * t5)) * pkin(4) - t17 * t172 + t16 * t171 + t5 * mrSges(6,1) - t6 * mrSges(6,2) + t20 * mrSges(5,1) - t21 * mrSges(5,2) + t3 * t108 + t42 * t184 + t112 * t187 + t114 * t186; t98 * mrSges(5,1) - t99 * mrSges(5,2) - t67 * mrSges(6,2) + (-mrSges(6,1) + t108) * t65 + t154 * mrSges(7,3) + m(7) * (t154 * t124 + t125 * t65) + m(6) * (t136 * t67 - t140 * t65) * pkin(4); t74 * t171 - t75 * t172 + t104 * t184 + t64 * t182 + t63 * t181 + t125 * t71 + m(7) * (t155 * t124 + t125 * t87) + t87 * t108 - t87 * mrSges(6,1) - t89 * mrSges(6,2) + (t112 * t183 + t114 * t181) * t105 + (-t145 * mrSges(5,1) - t148 * mrSges(5,2)) * pkin(10) + t155 * mrSges(7,3) + (m(6) * (t136 * t89 - t140 * t87) + (-t136 * t104 - t140 * t105) * mrSges(6,3)) * pkin(4) + t193; 0.2e1 * t125 * t108 + t147 * t112 + t144 * t114 + m(7) * (t163 * t124 ^ 2 + t125 ^ 2) + m(6) * (t136 ^ 2 + t140 ^ 2) * pkin(4) ^ 2 + 0.2e1 * (mrSges(6,1) * t140 - mrSges(6,2) * t136) * pkin(4) + 0.2e1 * t163 * t124 * mrSges(7,3) + t194; t144 * t16 + t147 * t17 + m(7) * (t1 * t147 + t144 * t2) + m(6) * t23 + t22; -m(6) * t168 + m(7) * (t144 * t59 + t147 * t58); t144 * t74 + t147 * t75 + m(7) * (t144 * t50 + t147 * t49) + m(6) * t126 + t83; 0; m(7) * t163 + m(6); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t58 - mrSges(7,2) * t59; mrSges(7,1) * t49 - mrSges(7,2) * t50 + t62; -t157 * t124 + t111; -t108; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
