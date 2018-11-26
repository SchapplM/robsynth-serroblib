% Calculate joint inertia matrix for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:44:27
% EndTime: 2018-11-23 17:44:29
% DurationCPUTime: 1.71s
% Computational Cost: add. (3349->402), mult. (7443->555), div. (0->0), fcn. (8277->10), ass. (0->143)
t202 = 2 * mrSges(7,3);
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t173 = t153 ^ 2 + t156 ^ 2;
t201 = 0.2e1 * t173;
t150 = sin(pkin(6));
t158 = cos(qJ(2));
t177 = t150 * t158;
t149 = sin(pkin(11));
t151 = cos(pkin(11));
t152 = cos(pkin(6));
t154 = sin(qJ(3));
t157 = cos(qJ(3));
t155 = sin(qJ(2));
t178 = t150 * t155;
t98 = t152 * t157 - t154 * t178;
t99 = t152 * t154 + t157 * t178;
t64 = t149 * t98 + t151 * t99;
t43 = -t153 * t64 - t156 * t177;
t44 = -t153 * t177 + t156 * t64;
t63 = t149 * t99 - t151 * t98;
t7 = Ifges(7,5) * t44 + Ifges(7,6) * t43 + Ifges(7,3) * t63;
t8 = Ifges(6,5) * t44 + Ifges(6,6) * t43 + Ifges(6,3) * t63;
t200 = t7 + t8;
t109 = t149 * t154 - t151 * t157;
t110 = t149 * t157 + t151 * t154;
t199 = Ifges(4,5) * t154 + Ifges(5,5) * t110 + Ifges(4,6) * t157 - Ifges(5,6) * t109;
t190 = -qJ(4) - pkin(9);
t116 = t190 * t157;
t163 = t190 * t154;
t80 = -t116 * t149 - t151 * t163;
t198 = t80 ^ 2;
t197 = 0.2e1 * t80;
t196 = m(7) * pkin(5);
t101 = t152 * t155 * pkin(1) + pkin(8) * t177;
t89 = pkin(9) * t152 + t101;
t90 = (-pkin(2) * t158 - pkin(9) * t155 - pkin(1)) * t150;
t53 = -t154 * t89 + t157 * t90;
t32 = -pkin(3) * t177 - qJ(4) * t99 + t53;
t54 = t154 * t90 + t157 * t89;
t37 = qJ(4) * t98 + t54;
t16 = t149 * t32 + t151 * t37;
t14 = -pkin(10) * t177 + t16;
t128 = pkin(8) * t178;
t194 = pkin(1) * t158;
t88 = t128 + (-pkin(2) - t194) * t152;
t67 = -pkin(3) * t98 + t88;
t19 = pkin(4) * t63 - pkin(10) * t64 + t67;
t4 = t156 * t14 + t153 * t19;
t193 = pkin(3) * t149;
t192 = pkin(3) * t151;
t191 = Ifges(4,3) + Ifges(5,3);
t136 = -pkin(3) * t157 - pkin(2);
t71 = pkin(4) * t109 - pkin(10) * t110 + t136;
t82 = -t151 * t116 + t149 * t163;
t34 = t153 * t71 + t156 * t82;
t179 = t110 * t156;
t180 = t110 * t153;
t68 = mrSges(7,1) * t180 + mrSges(7,2) * t179;
t189 = mrSges(6,2) * t153;
t188 = Ifges(6,4) * t153;
t187 = Ifges(6,4) * t156;
t186 = Ifges(7,4) * t153;
t185 = Ifges(7,4) * t156;
t100 = t152 * t194 - t128;
t184 = t100 * mrSges(3,1);
t183 = t101 * mrSges(3,2);
t182 = Ifges(7,5) * t179 + Ifges(7,3) * t109;
t181 = Ifges(6,5) * t179 + Ifges(6,3) * t109;
t134 = pkin(10) + t193;
t176 = qJ(6) + t134;
t117 = Ifges(7,5) * t153 + Ifges(7,6) * t156;
t118 = Ifges(6,5) * t153 + Ifges(6,6) * t156;
t172 = t154 ^ 2 + t157 ^ 2;
t10 = Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t63;
t9 = Ifges(7,4) * t44 + Ifges(7,2) * t43 + Ifges(7,6) * t63;
t171 = -t9 / 0.2e1 - t10 / 0.2e1;
t11 = Ifges(7,1) * t44 + Ifges(7,4) * t43 + Ifges(7,5) * t63;
t12 = Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t63;
t170 = t11 / 0.2e1 + t12 / 0.2e1;
t47 = Ifges(7,6) * t109 + (-Ifges(7,2) * t153 + t185) * t110;
t48 = Ifges(6,6) * t109 + (-Ifges(6,2) * t153 + t187) * t110;
t169 = t47 / 0.2e1 + t48 / 0.2e1;
t49 = Ifges(7,5) * t109 + (Ifges(7,1) * t156 - t186) * t110;
t50 = Ifges(6,5) * t109 + (Ifges(6,1) * t156 - t188) * t110;
t168 = t49 / 0.2e1 + t50 / 0.2e1;
t167 = Ifges(3,5) * t178 + Ifges(3,6) * t177 + Ifges(3,3) * t152;
t135 = -pkin(4) - t192;
t166 = t117 / 0.2e1 + t118 / 0.2e1;
t119 = Ifges(7,2) * t156 + t186;
t120 = Ifges(6,2) * t156 + t188;
t165 = t119 / 0.2e1 + t120 / 0.2e1;
t122 = Ifges(7,1) * t153 + t185;
t123 = Ifges(6,1) * t153 + t187;
t164 = t122 / 0.2e1 + t123 / 0.2e1;
t30 = t63 * mrSges(5,1) + t64 * mrSges(5,2);
t20 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t3 = -t14 * t153 + t156 * t19;
t15 = -t149 * t37 + t151 * t32;
t33 = -t153 * t82 + t156 * t71;
t162 = -Ifges(4,5) * t99 - Ifges(5,5) * t64 - Ifges(4,6) * t98 + Ifges(5,6) * t63;
t76 = t109 * mrSges(5,1) + t110 * mrSges(5,2);
t138 = t153 * mrSges(7,2);
t113 = -t156 * mrSges(7,1) + t138;
t13 = pkin(4) * t177 - t15;
t161 = mrSges(6,1) * t153 + mrSges(6,2) * t156;
t124 = Ifges(4,1) * t154 + Ifges(4,4) * t157;
t121 = Ifges(4,4) * t154 + Ifges(4,2) * t157;
t115 = -mrSges(4,1) * t157 + mrSges(4,2) * t154;
t114 = -mrSges(6,1) * t156 + t189;
t112 = -pkin(5) * t156 + t135;
t108 = t176 * t156;
t107 = t176 * t153;
t84 = -mrSges(4,1) * t177 - mrSges(4,3) * t99;
t83 = mrSges(4,2) * t177 + mrSges(4,3) * t98;
t78 = Ifges(5,1) * t110 - Ifges(5,4) * t109;
t77 = Ifges(5,4) * t110 - Ifges(5,2) * t109;
t75 = mrSges(6,1) * t109 - mrSges(6,3) * t179;
t74 = mrSges(7,1) * t109 - mrSges(7,3) * t179;
t73 = -mrSges(6,2) * t109 - mrSges(6,3) * t180;
t72 = -mrSges(7,2) * t109 - mrSges(7,3) * t180;
t70 = -mrSges(4,1) * t98 + mrSges(4,2) * t99;
t69 = t161 * t110;
t57 = Ifges(4,1) * t99 + Ifges(4,4) * t98 - Ifges(4,5) * t177;
t56 = Ifges(4,4) * t99 + Ifges(4,2) * t98 - Ifges(4,6) * t177;
t55 = pkin(5) * t180 + t80;
t52 = -mrSges(5,1) * t177 - mrSges(5,3) * t64;
t51 = mrSges(5,2) * t177 - mrSges(5,3) * t63;
t46 = -Ifges(6,6) * t180 + t181;
t45 = -Ifges(7,6) * t180 + t182;
t29 = -qJ(6) * t180 + t34;
t28 = Ifges(5,1) * t64 - Ifges(5,4) * t63 - Ifges(5,5) * t177;
t27 = Ifges(5,4) * t64 - Ifges(5,2) * t63 - Ifges(5,6) * t177;
t26 = mrSges(6,1) * t63 - mrSges(6,3) * t44;
t25 = mrSges(7,1) * t63 - mrSges(7,3) * t44;
t24 = -mrSges(6,2) * t63 + mrSges(6,3) * t43;
t23 = -mrSges(7,2) * t63 + mrSges(7,3) * t43;
t22 = pkin(5) * t109 - qJ(6) * t179 + t33;
t21 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t5 = -pkin(5) * t43 + t13;
t2 = qJ(6) * t43 + t4;
t1 = pkin(5) * t63 - qJ(6) * t44 + t3;
t6 = [(-t27 + t200) * t63 + ((-0.2e1 * t100 * mrSges(3,3) + Ifges(3,5) * t152 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t155) * t150) * t155 + (0.2e1 * t101 * mrSges(3,3) + Ifges(3,6) * t152 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t155 + (Ifges(3,2) + t191) * t158) * t150 + t162) * t158) * t150 + (t167 - 0.2e1 * t183 + 0.2e1 * t184) * t152 + t98 * t56 + t99 * t57 + 0.2e1 * t54 * t83 + 0.2e1 * t53 * t84 + 0.2e1 * t88 * t70 + 0.2e1 * t67 * t30 + t64 * t28 + 0.2e1 * t16 * t51 + 0.2e1 * t15 * t52 + 0.2e1 * t5 * t20 + 0.2e1 * t13 * t21 + 0.2e1 * t2 * t23 + 0.2e1 * t4 * t24 + 0.2e1 * t1 * t25 + 0.2e1 * t3 * t26 + m(4) * (t53 ^ 2 + t54 ^ 2 + t88 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t67 ^ 2) + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (pkin(1) ^ 2 * t150 ^ 2 + t100 ^ 2 + t101 ^ 2) + (t11 + t12) * t44 + (t10 + t9) * t43 + Ifges(2,3); -t183 + t184 - t199 * t177 / 0.2e1 + t167 + t168 * t44 + t169 * t43 + (t28 / 0.2e1 - t15 * mrSges(5,3) + t170 * t156 + t171 * t153) * t110 + t136 * t30 + t88 * t115 + t98 * t121 / 0.2e1 + t99 * t124 / 0.2e1 + t82 * t51 + t5 * t68 + t13 * t69 - pkin(2) * t70 + t2 * t72 + t4 * t73 + t1 * t74 + t3 * t75 + t67 * t76 + t64 * t78 / 0.2e1 + t55 * t20 + t29 * t23 + t33 * t26 + t34 * t24 + t22 * t25 + (-t77 / 0.2e1 + t45 / 0.2e1 + t46 / 0.2e1) * t63 + (t7 / 0.2e1 + t8 / 0.2e1 - t27 / 0.2e1 - t16 * mrSges(5,3)) * t109 + (t21 - t52) * t80 + (t56 / 0.2e1 + pkin(9) * t83 + t54 * mrSges(4,3)) * t157 + (t57 / 0.2e1 - pkin(9) * t84 - t53 * mrSges(4,3)) * t154 + m(4) * (-pkin(2) * t88 + (-t53 * t154 + t54 * t157) * pkin(9)) + m(6) * (t13 * t80 + t3 * t33 + t34 * t4) + m(7) * (t1 * t22 + t2 * t29 + t5 * t55) + m(5) * (t136 * t67 - t15 * t80 + t16 * t82); -0.2e1 * pkin(2) * t115 + t157 * t121 + t154 * t124 + 0.2e1 * t136 * t76 + 0.2e1 * t22 * t74 + 0.2e1 * t29 * t72 + 0.2e1 * t33 * t75 + 0.2e1 * t34 * t73 + 0.2e1 * t55 * t68 + t69 * t197 + Ifges(3,3) + 0.2e1 * t172 * pkin(9) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t82 + t45 + t46 - t77) * t109 + m(6) * (t33 ^ 2 + t34 ^ 2 + t198) + m(7) * (t22 ^ 2 + t29 ^ 2 + t55 ^ 2) + m(4) * (pkin(9) ^ 2 * t172 + pkin(2) ^ 2) + m(5) * (t136 ^ 2 + t82 ^ 2 + t198) + (mrSges(5,3) * t197 + t78 + (t49 + t50) * t156 + (-t47 - t48) * t153) * t110; -t162 + t135 * t21 + t112 * t20 + t5 * t113 - t107 * t25 + t108 * t23 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + t15 * mrSges(5,1) - t16 * mrSges(5,2) + (t2 * mrSges(7,3) + t4 * mrSges(6,3) + (m(6) * t4 + t24) * t134 - t171) * t156 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + (-m(6) * t3 - t26) * t134 + t170) * t153 + (t151 * t52 + m(5) * (t149 * t16 + t15 * t151) + t149 * t51) * pkin(3) - t191 * t177 + m(7) * (-t1 * t107 + t108 * t2 + t112 * t5) + t166 * t63 + t164 * t44 + t165 * t43 + (m(6) * t135 + t114) * t13; -t82 * mrSges(5,2) - t107 * t74 + t108 * t72 + t112 * t68 + t55 * t113 + t135 * t69 + (t114 - mrSges(5,1)) * t80 + (-mrSges(4,1) * t154 - mrSges(4,2) * t157) * pkin(9) + (-mrSges(5,3) * t193 + t166) * t109 + (t34 * mrSges(6,3) + t29 * mrSges(7,3) + t134 * t73 + t169) * t156 + (-t33 * mrSges(6,3) - t22 * mrSges(7,3) - t134 * t75 + t168) * t153 + m(6) * (t135 * t80 + (-t153 * t33 + t156 * t34) * t134) + m(7) * (-t107 * t22 + t108 * t29 + t112 * t55) + m(5) * (t149 * t82 - t151 * t80) * pkin(3) + (-mrSges(5,3) * t192 - t153 * t165 + t156 * t164) * t110 + t199; 0.2e1 * t112 * t113 + 0.2e1 * t135 * t114 + (t108 * t202 + t119 + t120) * t156 + (t107 * t202 + t122 + t123) * t153 + m(7) * (t107 ^ 2 + t108 ^ 2 + t112 ^ 2) + m(6) * (t134 ^ 2 * t173 + t135 ^ 2) + m(5) * (t149 ^ 2 + t151 ^ 2) * pkin(3) ^ 2 + t191 + 0.2e1 * (mrSges(5,1) * t151 - mrSges(5,2) * t149) * pkin(3) + t134 * mrSges(6,3) * t201; (t25 + t26) * t156 + (t23 + t24) * t153 + m(7) * (t1 * t156 + t153 * t2) + m(6) * (t153 * t4 + t156 * t3) + m(5) * t67 + t30; (t74 + t75) * t156 + (t72 + t73) * t153 + m(7) * (t153 * t29 + t156 * t22) + m(6) * (t153 * t34 + t156 * t33) + m(5) * t136 + t76; m(7) * (-t107 * t156 + t108 * t153); m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t201; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t200; mrSges(6,1) * t33 + mrSges(7,1) * t22 - mrSges(6,2) * t34 - mrSges(7,2) * t29 + (-Ifges(6,6) - Ifges(7,6)) * t180 + (m(7) * t22 + t74) * pkin(5) + t181 + t182; -mrSges(7,1) * t107 - mrSges(7,2) * t108 - t161 * t134 + (-m(7) * t107 - t153 * mrSges(7,3)) * pkin(5) + t118 + t117; -t189 - t138 + (mrSges(6,1) + mrSges(7,1) + t196) * t156; Ifges(6,3) + Ifges(7,3) + (0.2e1 * mrSges(7,1) + t196) * pkin(5); m(7) * t5 + t20; m(7) * t55 + t68; m(7) * t112 + t113; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
