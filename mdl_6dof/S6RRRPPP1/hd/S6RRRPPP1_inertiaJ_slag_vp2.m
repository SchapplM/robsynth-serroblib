% Calculate joint inertia matrix for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2018-11-23 17:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:31:56
% EndTime: 2018-11-23 17:31:58
% DurationCPUTime: 1.94s
% Computational Cost: add. (2521->436), mult. (6030->588), div. (0->0), fcn. (6130->8), ass. (0->165)
t203 = 2 * pkin(8);
t150 = sin(pkin(6));
t153 = sin(qJ(3));
t155 = cos(qJ(3));
t108 = -qJ(4) * t150 * t153 - pkin(3) * t155 - pkin(2);
t152 = cos(pkin(6));
t191 = qJ(4) * t152;
t175 = pkin(9) + t191;
t117 = t175 * t153;
t151 = cos(pkin(10));
t202 = t151 * (t108 * t150 - t117 * t152);
t156 = cos(qJ(2));
t154 = sin(qJ(2));
t183 = t153 * t154;
t160 = t150 * t156 + t152 * t183;
t149 = sin(pkin(10));
t122 = -pkin(2) * t156 - pkin(9) * t154 - pkin(1);
t197 = pkin(8) * t156;
t83 = t153 * t122 + t155 * t197;
t48 = -qJ(4) * t160 + t83;
t110 = t155 * t122;
t182 = t154 * t155;
t52 = -t182 * t191 + t110 + (-pkin(8) * t153 - pkin(3)) * t156;
t187 = t150 * t155;
t71 = (pkin(3) * t153 - qJ(4) * t187 + pkin(8)) * t154;
t7 = -t149 * t48 + (t150 * t71 + t152 * t52) * t151;
t201 = m(5) * pkin(3);
t200 = m(6) + m(7);
t25 = -t150 * t52 + t152 * t71;
t199 = m(5) * t25;
t53 = t152 * t108 + t117 * t150;
t198 = m(5) * t53;
t196 = -mrSges(6,2) + mrSges(5,1);
t195 = pkin(4) + qJ(6);
t194 = Ifges(4,4) * t153;
t193 = Ifges(4,4) * t155;
t192 = Ifges(4,6) * t156;
t184 = t152 * t155;
t98 = t149 * t184 + t151 * t153;
t64 = t98 * mrSges(7,1) + mrSges(7,3) * t187;
t190 = t149 * t150;
t189 = t149 * t152;
t188 = t150 * t151;
t185 = t151 * t152;
t118 = t175 * t155;
t104 = t149 * t118;
t181 = pkin(4) * t187 + t104;
t103 = pkin(3) * t189 + qJ(4) * t188;
t115 = mrSges(7,1) * t188 + t152 * mrSges(7,2);
t116 = mrSges(6,1) * t190 + t152 * mrSges(6,2);
t180 = Ifges(4,5) * t153 + Ifges(4,6) * t155;
t179 = t153 ^ 2 + t155 ^ 2;
t8 = t151 * t48 + t52 * t189 + t71 * t190;
t29 = t108 * t190 - t117 * t189 + t151 * t118;
t177 = -pkin(3) * t151 - pkin(4);
t106 = t150 * t183 - t152 * t156;
t60 = -t160 * t149 + t151 * t182;
t35 = t60 * mrSges(6,1) + t106 * mrSges(6,2);
t59 = t149 * t182 + t151 * t160;
t34 = -t59 * mrSges(7,1) + t106 * mrSges(7,2);
t22 = -t60 * mrSges(7,2) + t59 * mrSges(7,3);
t97 = t149 * t153 - t151 * t184;
t49 = -t98 * mrSges(7,2) + t97 * mrSges(7,3);
t32 = t60 * mrSges(7,1) - t106 * mrSges(7,3);
t176 = -qJ(5) * t149 - pkin(3);
t113 = mrSges(7,1) * t190 - t152 * mrSges(7,3);
t11 = Ifges(7,5) * t106 + Ifges(7,6) * t59 + Ifges(7,3) * t60;
t15 = Ifges(6,4) * t106 - Ifges(6,2) * t60 + Ifges(6,6) * t59;
t19 = Ifges(5,1) * t60 - Ifges(5,4) * t59 + Ifges(5,5) * t106;
t174 = t11 / 0.2e1 - t15 / 0.2e1 + t19 / 0.2e1;
t12 = Ifges(6,5) * t106 - Ifges(6,6) * t60 + Ifges(6,3) * t59;
t14 = Ifges(7,4) * t106 + Ifges(7,2) * t59 + Ifges(7,6) * t60;
t16 = Ifges(5,4) * t60 - Ifges(5,2) * t59 + Ifges(5,6) * t106;
t173 = t12 / 0.2e1 + t14 / 0.2e1 - t16 / 0.2e1;
t13 = Ifges(5,5) * t60 - Ifges(5,6) * t59 + Ifges(5,3) * t106;
t17 = Ifges(7,1) * t106 + Ifges(7,4) * t59 + Ifges(7,5) * t60;
t18 = Ifges(6,1) * t106 - Ifges(6,4) * t60 + Ifges(6,5) * t59;
t172 = t13 / 0.2e1 + t17 / 0.2e1 + t18 / 0.2e1;
t36 = -Ifges(7,5) * t187 + Ifges(7,6) * t97 + Ifges(7,3) * t98;
t39 = -Ifges(6,4) * t187 - Ifges(6,2) * t98 + Ifges(6,6) * t97;
t44 = Ifges(5,1) * t98 - Ifges(5,4) * t97 - Ifges(5,5) * t187;
t171 = -t39 / 0.2e1 + t44 / 0.2e1 + t36 / 0.2e1;
t40 = -Ifges(7,1) * t187 + Ifges(7,4) * t97 + Ifges(7,5) * t98;
t41 = -Ifges(6,1) * t187 - Ifges(6,4) * t98 + Ifges(6,5) * t97;
t42 = Ifges(5,5) * t98 - Ifges(5,6) * t97 - Ifges(5,3) * t187;
t170 = t40 / 0.2e1 + t41 / 0.2e1 + t42 / 0.2e1;
t37 = -Ifges(6,5) * t187 - Ifges(6,6) * t98 + Ifges(6,3) * t97;
t38 = -Ifges(7,4) * t187 + Ifges(7,2) * t97 + Ifges(7,6) * t98;
t43 = Ifges(5,4) * t98 - Ifges(5,2) * t97 - Ifges(5,6) * t187;
t169 = t43 / 0.2e1 - t37 / 0.2e1 - t38 / 0.2e1;
t73 = Ifges(5,3) * t152 + (Ifges(5,5) * t149 + Ifges(5,6) * t151) * t150;
t80 = Ifges(7,1) * t152 + (-Ifges(7,4) * t151 + Ifges(7,5) * t149) * t150;
t81 = Ifges(6,1) * t152 + (-Ifges(6,4) * t149 - Ifges(6,5) * t151) * t150;
t168 = t73 / 0.2e1 + t80 / 0.2e1 + t81 / 0.2e1;
t75 = Ifges(5,5) * t152 + (Ifges(5,1) * t149 + Ifges(5,4) * t151) * t150;
t76 = Ifges(7,5) * t152 + (-Ifges(7,6) * t151 + Ifges(7,3) * t149) * t150;
t79 = Ifges(6,4) * t152 + (-Ifges(6,2) * t149 - Ifges(6,6) * t151) * t150;
t167 = t75 / 0.2e1 + t76 / 0.2e1 - t79 / 0.2e1;
t74 = Ifges(5,6) * t152 + (Ifges(5,4) * t149 + Ifges(5,2) * t151) * t150;
t77 = Ifges(6,5) * t152 + (-Ifges(6,6) * t149 - Ifges(6,3) * t151) * t150;
t78 = Ifges(7,4) * t152 + (-Ifges(7,2) * t151 + Ifges(7,6) * t149) * t150;
t166 = t77 / 0.2e1 + t78 / 0.2e1 - t74 / 0.2e1;
t72 = -qJ(5) * t152 - t103;
t165 = mrSges(4,1) * t153 + mrSges(4,2) * t155;
t4 = -qJ(5) * t106 - t8;
t67 = t98 * mrSges(6,1) - mrSges(6,2) * t187;
t66 = -t97 * mrSges(7,1) - mrSges(7,2) * t187;
t162 = -qJ(5) * t60 + t25;
t161 = -qJ(5) * t98 + t53;
t26 = qJ(5) * t187 - t29;
t158 = pkin(8) ^ 2;
t148 = t156 ^ 2;
t146 = t154 ^ 2;
t144 = t146 * t158;
t138 = Ifges(4,5) * t182;
t133 = mrSges(6,2) * t188;
t130 = mrSges(5,2) * t190;
t128 = qJ(4) * t190;
t125 = Ifges(4,1) * t153 + t193;
t124 = Ifges(4,2) * t155 + t194;
t123 = -mrSges(4,1) * t155 + mrSges(4,2) * t153;
t120 = -mrSges(4,1) * t156 - mrSges(4,3) * t182;
t119 = mrSges(4,2) * t156 - mrSges(4,3) * t183;
t114 = -mrSges(6,1) * t188 - mrSges(6,3) * t152;
t112 = -mrSges(5,2) * t152 + mrSges(5,3) * t188;
t111 = mrSges(5,1) * t152 - mrSges(5,3) * t190;
t107 = t165 * t154;
t102 = (-mrSges(7,2) * t149 - t151 * mrSges(7,3)) * t150;
t101 = -mrSges(5,1) * t188 + t130;
t100 = -mrSges(6,3) * t190 + t133;
t99 = pkin(3) * t185 - t128;
t93 = t98 * mrSges(6,3);
t91 = t98 * mrSges(5,2);
t88 = -Ifges(4,5) * t156 + (Ifges(4,1) * t155 - t194) * t154;
t87 = -t192 + (-Ifges(4,2) * t153 + t193) * t154;
t85 = (-pkin(4) * t151 + t176) * t150;
t84 = t152 * t177 + t128;
t82 = -t153 * t197 + t110;
t69 = -mrSges(5,1) * t187 - mrSges(5,3) * t98;
t68 = mrSges(5,2) * t187 - mrSges(5,3) * t97;
t65 = mrSges(6,1) * t97 + mrSges(6,3) * t187;
t62 = (-t151 * t195 + t176) * t150;
t61 = pkin(5) * t188 - t72;
t58 = t60 * mrSges(6,3);
t56 = t60 * mrSges(5,2);
t54 = pkin(5) * t190 + t128 + (-qJ(6) + t177) * t152;
t51 = -t97 * mrSges(6,2) - t93;
t50 = t97 * mrSges(5,1) + t91;
t33 = mrSges(6,1) * t59 - mrSges(6,3) * t106;
t31 = mrSges(5,1) * t106 - mrSges(5,3) * t60;
t30 = -mrSges(5,2) * t106 - mrSges(5,3) * t59;
t28 = -t104 + t202;
t27 = t181 - t202;
t24 = -t59 * mrSges(6,2) - t58;
t23 = t59 * mrSges(5,1) + t56;
t21 = pkin(4) * t97 + t161;
t20 = -pkin(5) * t97 - t26;
t10 = t195 * t97 + t161;
t9 = t117 * t185 + pkin(5) * t98 + (qJ(6) * t155 - t108 * t151) * t150 + t181;
t6 = pkin(4) * t59 + t162;
t5 = -pkin(4) * t106 - t7;
t3 = t195 * t59 + t162;
t2 = -pkin(5) * t59 - t4;
t1 = pkin(5) * t60 - t106 * t195 - t7;
t45 = [0.2e1 * t1 * t32 + 0.2e1 * t83 * t119 + 0.2e1 * t82 * t120 + 0.2e1 * t2 * t34 + 0.2e1 * t3 * t22 + 0.2e1 * t25 * t23 + 0.2e1 * t6 * t24 + 0.2e1 * t8 * t30 + 0.2e1 * t7 * t31 + 0.2e1 * t4 * t33 + 0.2e1 * t5 * t35 + Ifges(2,3) + (t146 + t148) * mrSges(3,3) * t203 + (t19 - t15 + t11) * t60 + (t12 + t14 - t16) * t59 + (t13 + t17 + t18) * t106 + (0.2e1 * pkin(1) * mrSges(3,1) - t138 + (Ifges(4,3) + Ifges(3,2)) * t156) * t156 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t154 + 0.2e1 * Ifges(3,4) * t156 + t107 * t203 + t155 * t88 + (-t87 + t192) * t153) * t154 + m(3) * (pkin(1) ^ 2 + t148 * t158 + t144) + m(4) * (t82 ^ 2 + t83 ^ 2 + t144) + m(5) * (t25 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2); (Ifges(3,6) - pkin(8) * mrSges(3,2) - t180 / 0.2e1) * t156 + m(4) * (-t153 * t82 + t155 * t83) * pkin(9) - pkin(2) * t107 + t1 * t64 + t4 * t65 + t2 * t66 + t5 * t67 + t8 * t68 + t7 * t69 + (-pkin(9) * t120 - t82 * mrSges(4,3) + t88 / 0.2e1) * t153 + (pkin(9) * t119 + t83 * mrSges(4,3) + t87 / 0.2e1 - t172 * t150) * t155 + t3 * t49 + t25 * t50 + t6 * t51 + t53 * t23 + t29 * t30 + t28 * t31 + t9 * t32 + t26 * t33 + t20 * t34 + t27 * t35 + t10 * t22 + t21 * t24 + m(5) * (t25 * t53 + t28 * t7 + t29 * t8) + m(6) * (t21 * t6 + t26 * t4 + t27 * t5) + m(7) * (t1 * t9 + t10 * t3 + t2 * t20) + (Ifges(3,5) - t153 * t124 / 0.2e1 + t155 * t125 / 0.2e1 + (-m(4) * pkin(2) - mrSges(3,1) + t123) * pkin(8)) * t154 - t169 * t59 + t170 * t106 + t171 * t60 + t173 * t97 + t174 * t98; -0.2e1 * pkin(2) * t123 + 0.2e1 * t10 * t49 + t153 * t125 + 0.2e1 * t20 * t66 + 0.2e1 * t21 * t51 + 0.2e1 * t26 * t65 + 0.2e1 * t27 * t67 + 0.2e1 * t28 * t69 + 0.2e1 * t29 * t68 + 0.2e1 * t53 * t50 + 0.2e1 * t9 * t64 + Ifges(3,3) + 0.2e1 * t179 * pkin(9) * mrSges(4,3) + (t36 + t44 - t39) * t98 + (t37 + t38 - t43) * t97 + (t124 + (-t40 - t41 - t42) * t150) * t155 + m(4) * (pkin(9) ^ 2 * t179 + pkin(2) ^ 2) + m(5) * (t28 ^ 2 + t29 ^ 2 + t53 ^ 2) + m(6) * (t21 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(7) * (t10 ^ 2 + t20 ^ 2 + t9 ^ 2); m(5) * (t103 * t8 + t7 * t99) + t138 - Ifges(4,3) * t156 + t7 * t111 + t8 * t112 + t1 * t113 + t4 * t114 + t2 * t115 + t5 * t116 + t99 * t31 + t6 * t100 + t25 * t101 + t3 * t102 + t103 * t30 + t72 * t33 + t82 * mrSges(4,1) - t83 * mrSges(4,2) + t84 * t35 + t85 * t24 + t61 * t34 + t62 * t22 + m(6) * (t4 * t72 + t5 * t84 + t6 * t85) + m(7) * (t1 * t54 + t2 * t61 + t3 * t62) + t54 * t32 + t166 * t59 + t167 * t60 + t168 * t106 + t172 * t152 - Ifges(4,6) * t183 + ((-t23 - t199) * pkin(3) - t173 * t151 + t174 * t149) * t150; t180 + t28 * t111 + t29 * t112 + t9 * t113 + t26 * t114 + t20 * t115 + t27 * t116 + t99 * t69 + t21 * t100 + t53 * t101 + t10 * t102 + t103 * t68 + t72 * t65 + t84 * t67 + t85 * t51 + t62 * t49 + t54 * t64 + t61 * t66 + t166 * t97 + t167 * t98 + t170 * t152 - t165 * pkin(9) + m(5) * (t103 * t29 + t28 * t99) + m(6) * (t21 * t85 + t26 * t72 + t27 * t84) + m(7) * (t10 * t62 + t20 * t61 + t54 * t9) + ((-t50 - t198) * pkin(3) - t168 * t155 + t169 * t151 + t171 * t149) * t150; 0.2e1 * t85 * t100 + 0.2e1 * t62 * t102 + 0.2e1 * t103 * t112 + 0.2e1 * t99 * t111 + 0.2e1 * t54 * t113 + 0.2e1 * t72 * t114 + 0.2e1 * t61 * t115 + 0.2e1 * t84 * t116 + Ifges(4,3) + (t73 + t80 + t81) * t152 + m(6) * (t72 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(7) * (t54 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t103 ^ 2 + t99 ^ 2) + ((t150 * t201 - 0.2e1 * t101) * pkin(3) + (t74 - t77 - t78) * t151 + (t75 + t76 - t79) * t149) * t150; m(6) * t6 + m(7) * t3 + t196 * t59 + t199 + t22 + t56 - t58; m(6) * t21 + m(7) * t10 + t196 * t97 + t198 + t49 + t91 - t93; t130 + t133 + m(6) * t85 + m(7) * t62 + (-t201 + (-mrSges(5,1) - mrSges(7,3)) * t151 + (-mrSges(7,2) - mrSges(6,3)) * t149) * t150; m(5) + t200; m(6) * t5 + m(7) * t1 + t32 + t35; m(6) * t27 + m(7) * t9 + t64 + t67; m(6) * t84 + m(7) * t54 + t113 + t116; 0; t200; m(7) * t2 + t34; m(7) * t20 + t66; m(7) * t61 + t115; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t45(1) t45(2) t45(4) t45(7) t45(11) t45(16); t45(2) t45(3) t45(5) t45(8) t45(12) t45(17); t45(4) t45(5) t45(6) t45(9) t45(13) t45(18); t45(7) t45(8) t45(9) t45(10) t45(14) t45(19); t45(11) t45(12) t45(13) t45(14) t45(15) t45(20); t45(16) t45(17) t45(18) t45(19) t45(20) t45(21);];
Mq  = res;
