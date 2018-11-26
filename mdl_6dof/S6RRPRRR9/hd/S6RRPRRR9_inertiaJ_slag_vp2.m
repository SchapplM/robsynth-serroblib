% Calculate joint inertia matrix for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:26:39
% EndTime: 2018-11-23 17:26:41
% DurationCPUTime: 1.88s
% Computational Cost: add. (5029->366), mult. (11146->529), div. (0->0), fcn. (12974->12), ass. (0->141)
t155 = sin(qJ(6));
t159 = cos(qJ(6));
t156 = sin(qJ(5));
t160 = cos(qJ(5));
t157 = sin(qJ(4));
t161 = cos(qJ(4));
t153 = cos(pkin(12));
t196 = pkin(9) + qJ(3);
t175 = t196 * t153;
t151 = sin(pkin(12));
t176 = t196 * t151;
t101 = -t157 * t175 - t161 * t176;
t124 = t151 * t161 + t153 * t157;
t168 = -t124 * pkin(10) + t101;
t102 = -t157 * t176 + t161 * t175;
t123 = -t151 * t157 + t153 * t161;
t83 = pkin(10) * t123 + t102;
t51 = t156 * t168 + t160 * t83;
t140 = -pkin(3) * t153 - pkin(2);
t108 = -pkin(4) * t123 + t140;
t96 = -t160 * t123 + t124 * t156;
t97 = t123 * t156 + t124 * t160;
t52 = pkin(5) * t96 - pkin(11) * t97 + t108;
t23 = -t155 * t51 + t159 * t52;
t24 = t155 * t52 + t159 * t51;
t171 = -t155 * t23 + t159 * t24;
t154 = cos(pkin(6));
t152 = sin(pkin(6));
t158 = sin(qJ(2));
t183 = t152 * t158;
t113 = -t151 * t183 + t153 * t154;
t114 = t151 * t154 + t153 * t183;
t87 = t113 * t161 - t114 * t157;
t88 = t113 * t157 + t114 * t161;
t57 = t156 * t88 - t160 * t87;
t58 = t156 * t87 + t160 * t88;
t136 = pkin(8) * t183;
t162 = cos(qJ(2));
t202 = pkin(1) * t162;
t111 = t136 + (-pkin(2) - t202) * t154;
t90 = -pkin(3) * t113 + t111;
t60 = -pkin(4) * t87 + t90;
t15 = pkin(5) * t57 - pkin(11) * t58 + t60;
t182 = t152 * t162;
t117 = t154 * t158 * pkin(1) + pkin(8) * t182;
t109 = qJ(3) * t154 + t117;
t110 = (-pkin(2) * t162 - qJ(3) * t158 - pkin(1)) * t152;
t76 = -t109 * t151 + t153 * t110;
t66 = -pkin(3) * t182 - pkin(9) * t114 + t76;
t77 = t153 * t109 + t151 * t110;
t71 = pkin(9) * t113 + t77;
t29 = -t157 * t71 + t161 * t66;
t18 = -pkin(4) * t182 - pkin(10) * t88 + t29;
t30 = t157 * t66 + t161 * t71;
t25 = pkin(10) * t87 + t30;
t9 = t156 * t18 + t160 * t25;
t6 = -pkin(11) * t182 + t9;
t2 = t15 * t159 - t155 * t6;
t3 = t15 * t155 + t159 * t6;
t173 = -t155 * t2 + t159 * t3;
t211 = -Ifges(5,5) * t88 - Ifges(5,6) * t87;
t49 = t156 * t83 - t160 * t168;
t210 = t49 ^ 2;
t209 = 0.2e1 * t49;
t41 = -t155 * t58 - t159 * t182;
t208 = t41 / 0.2e1;
t131 = Ifges(7,5) * t155 + Ifges(7,6) * t159;
t207 = t131 / 0.2e1;
t192 = Ifges(7,4) * t159;
t133 = Ifges(7,1) * t155 + t192;
t206 = t133 / 0.2e1;
t205 = t155 / 0.2e1;
t204 = t159 / 0.2e1;
t201 = pkin(11) * t155;
t200 = pkin(11) * t159;
t197 = Ifges(5,3) + Ifges(6,3);
t195 = -Ifges(6,5) * t58 + Ifges(6,6) * t57;
t194 = Ifges(6,5) * t97 - Ifges(6,6) * t96;
t193 = Ifges(7,4) * t155;
t116 = t154 * t202 - t136;
t191 = t116 * mrSges(3,1);
t190 = t117 * mrSges(3,2);
t188 = t155 * t97;
t186 = t159 * t97;
t141 = pkin(4) * t156 + pkin(11);
t185 = t141 * t155;
t184 = t141 * t159;
t181 = Ifges(5,5) * t124 + Ifges(5,6) * t123;
t180 = t151 ^ 2 + t153 ^ 2;
t179 = t155 ^ 2 + t159 ^ 2;
t42 = -t155 * t182 + t159 * t58;
t12 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t57;
t132 = Ifges(7,2) * t159 + t193;
t178 = t159 * t132 + t155 * t133 + Ifges(6,3);
t177 = Ifges(3,5) * t183 + Ifges(3,6) * t182 + Ifges(3,3) * t154;
t59 = -t87 * mrSges(5,1) + t88 * mrSges(5,2);
t28 = t57 * mrSges(6,1) + t58 * mrSges(6,2);
t67 = t96 * mrSges(6,1) + t97 * mrSges(6,2);
t91 = -t113 * mrSges(4,1) + t114 * mrSges(4,2);
t127 = -t153 * mrSges(4,1) + t151 * mrSges(4,2);
t98 = -t123 * mrSges(5,1) + t124 * mrSges(5,2);
t174 = t179 * t141;
t172 = mrSges(7,1) * t155 + mrSges(7,2) * t159;
t8 = -t156 * t25 + t160 * t18;
t170 = 0.2e1 * mrSges(7,3) * t179;
t33 = Ifges(7,5) * t186 - Ifges(7,6) * t188 + Ifges(7,3) * t96;
t169 = (mrSges(6,1) * t160 - mrSges(6,2) * t156) * pkin(4);
t13 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t57;
t130 = -mrSges(7,1) * t159 + mrSges(7,2) * t155;
t14 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t57;
t5 = pkin(5) * t182 - t8;
t167 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + t173 * mrSges(7,3) + t13 * t204 + t5 * t130 + t132 * t208 + t14 * t205 + t42 * t206 + t57 * t207 - t195;
t34 = Ifges(7,6) * t96 + (-Ifges(7,2) * t155 + t192) * t97;
t35 = Ifges(7,5) * t96 + (Ifges(7,1) * t159 - t193) * t97;
t166 = -t51 * mrSges(6,2) + t34 * t204 + t35 * t205 + t194 - t132 * t188 / 0.2e1 + t186 * t206 + t96 * t207 + (-mrSges(6,1) + t130) * t49 + t171 * mrSges(7,3);
t142 = -pkin(4) * t160 - pkin(5);
t129 = Ifges(4,1) * t151 + Ifges(4,4) * t153;
t128 = Ifges(4,4) * t151 + Ifges(4,2) * t153;
t104 = -mrSges(4,1) * t182 - mrSges(4,3) * t114;
t103 = mrSges(4,2) * t182 + mrSges(4,3) * t113;
t100 = Ifges(5,1) * t124 + Ifges(5,4) * t123;
t99 = Ifges(5,4) * t124 + Ifges(5,2) * t123;
t81 = Ifges(4,1) * t114 + Ifges(4,4) * t113 - Ifges(4,5) * t182;
t80 = Ifges(4,4) * t114 + Ifges(4,2) * t113 - Ifges(4,6) * t182;
t73 = -mrSges(5,1) * t182 - mrSges(5,3) * t88;
t72 = mrSges(5,2) * t182 + mrSges(5,3) * t87;
t69 = Ifges(6,1) * t97 - Ifges(6,4) * t96;
t68 = Ifges(6,4) * t97 - Ifges(6,2) * t96;
t65 = mrSges(7,1) * t96 - mrSges(7,3) * t186;
t64 = -mrSges(7,2) * t96 - mrSges(7,3) * t188;
t61 = t172 * t97;
t48 = Ifges(5,1) * t88 + Ifges(5,4) * t87 - Ifges(5,5) * t182;
t47 = Ifges(5,4) * t88 + Ifges(5,2) * t87 - Ifges(5,6) * t182;
t45 = -mrSges(6,1) * t182 - mrSges(6,3) * t58;
t44 = mrSges(6,2) * t182 - mrSges(6,3) * t57;
t27 = Ifges(6,1) * t58 - Ifges(6,4) * t57 - Ifges(6,5) * t182;
t26 = Ifges(6,4) * t58 - Ifges(6,2) * t57 - Ifges(6,6) * t182;
t21 = mrSges(7,1) * t57 - mrSges(7,3) * t42;
t20 = -mrSges(7,2) * t57 + mrSges(7,3) * t41;
t16 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t1 = [m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + m(6) * (t60 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2 + t90 ^ 2) + m(4) * (t111 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (pkin(1) ^ 2 * t152 ^ 2 + t116 ^ 2 + t117 ^ 2) + (t12 - t26) * t57 + ((-0.2e1 * t116 * mrSges(3,3) + Ifges(3,5) * t154 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t158) * t152) * t158 + (0.2e1 * t117 * mrSges(3,3) - Ifges(4,5) * t114 + Ifges(3,6) * t154 - Ifges(4,6) * t113 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t158 + (Ifges(4,3) + Ifges(3,2) + t197) * t162) * t152 + t195 + t211) * t162) * t152 + (t177 - 0.2e1 * t190 + 0.2e1 * t191) * t154 + Ifges(2,3) + 0.2e1 * t5 * t16 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + t41 * t13 + t42 * t14 + 0.2e1 * t9 * t44 + 0.2e1 * t8 * t45 + t58 * t27 + 0.2e1 * t60 * t28 + 0.2e1 * t30 * t72 + 0.2e1 * t29 * t73 + t87 * t47 + t88 * t48 + 0.2e1 * t90 * t59 + 0.2e1 * t77 * t103 + 0.2e1 * t76 * t104 + 0.2e1 * t111 * t91 + t113 * t80 + t114 * t81; (t16 - t45) * t49 - t190 + t191 + (t77 * mrSges(4,3) + qJ(3) * t103 + t80 / 0.2e1) * t153 + (-t76 * mrSges(4,3) - qJ(3) * t104 + t81 / 0.2e1) * t151 + m(4) * (-pkin(2) * t111 + (-t76 * t151 + t77 * t153) * qJ(3)) + t34 * t208 + m(7) * (t2 * t23 + t24 * t3 + t49 * t5) + m(6) * (t108 * t60 - t49 * t8 + t51 * t9) + m(5) * (t101 * t29 + t102 * t30 + t140 * t90) - (Ifges(4,5) * t151 + Ifges(4,6) * t153 + t181 + t194) * t182 / 0.2e1 + (-t8 * mrSges(6,3) - t155 * t13 / 0.2e1 + t14 * t204 + t27 / 0.2e1) * t97 + (-t9 * mrSges(6,3) + t12 / 0.2e1 - t26 / 0.2e1) * t96 + t177 + (t33 / 0.2e1 - t68 / 0.2e1) * t57 + t23 * t21 + t24 * t20 + t42 * t35 / 0.2e1 + t51 * t44 + t5 * t61 + t3 * t64 + t2 * t65 + t60 * t67 + t58 * t69 / 0.2e1 - pkin(2) * t91 + t90 * t98 + t87 * t99 / 0.2e1 + t88 * t100 / 0.2e1 + t101 * t73 + t102 * t72 + t108 * t28 + t123 * t47 / 0.2e1 + t124 * t48 / 0.2e1 + t111 * t127 + t113 * t128 / 0.2e1 + t114 * t129 / 0.2e1 + t140 * t59 + (t30 * t123 - t29 * t124) * mrSges(5,3); -0.2e1 * pkin(2) * t127 + t124 * t100 + 0.2e1 * t108 * t67 + t123 * t99 + t153 * t128 + t151 * t129 + 0.2e1 * t140 * t98 + 0.2e1 * t23 * t65 + 0.2e1 * t24 * t64 + t61 * t209 + Ifges(3,3) + (-0.2e1 * mrSges(6,3) * t51 + t33 - t68) * t96 + (mrSges(6,3) * t209 - t155 * t34 + t159 * t35 + t69) * t97 + m(7) * (t23 ^ 2 + t24 ^ 2 + t210) + m(6) * (t108 ^ 2 + t51 ^ 2 + t210) + m(5) * (t101 ^ 2 + t102 ^ 2 + t140 ^ 2) + m(4) * (qJ(3) ^ 2 * t180 + pkin(2) ^ 2) + 0.2e1 * (-t101 * t124 + t102 * t123) * mrSges(5,3) + 0.2e1 * t180 * qJ(3) * mrSges(4,3); t155 * t20 + t159 * t21 + m(7) * (t155 * t3 + t159 * t2) + m(6) * t60 + m(5) * t90 + m(4) * t111 + t91 + t28 + t59; -m(4) * pkin(2) + t155 * t64 + t159 * t65 + m(7) * (t155 * t24 + t159 * t23) + m(6) * t108 + m(5) * t140 + t98 + t127 + t67; m(7) * t179 + m(4) + m(5) + m(6); t167 + m(7) * (t141 * t173 + t142 * t5) + (t156 * t44 + t160 * t45 + m(6) * (t156 * t9 + t160 * t8)) * pkin(4) - t197 * t182 - t21 * t185 + t20 * t184 + t29 * mrSges(5,1) - t30 * mrSges(5,2) + t142 * t16 - t211; t166 + m(7) * (t141 * t171 + t142 * t49) + (m(6) * (t156 * t51 - t160 * t49) + (-t156 * t96 - t160 * t97) * mrSges(6,3)) * pkin(4) - t65 * t185 + t64 * t184 + t101 * mrSges(5,1) - t102 * mrSges(5,2) + t142 * t61 + t181; 0; 0.2e1 * t142 * t130 + Ifges(5,3) + 0.2e1 * t169 + t141 * t170 + m(7) * (t141 ^ 2 * t179 + t142 ^ 2) + m(6) * (t156 ^ 2 + t160 ^ 2) * pkin(4) ^ 2 + t178; t167 + m(7) * (-pkin(5) * t5 + pkin(11) * t173) - t21 * t201 + t20 * t200 - Ifges(6,3) * t182 - pkin(5) * t16; t166 + m(7) * (-pkin(5) * t49 + pkin(11) * t171) - t65 * t201 + t64 * t200 - pkin(5) * t61; 0; m(7) * (-pkin(5) * t142 + pkin(11) * t174) + (t142 - pkin(5)) * t130 + t169 + (pkin(11) * t179 + t174) * mrSges(7,3) + t178; -0.2e1 * pkin(5) * t130 + m(7) * (pkin(11) ^ 2 * t179 + pkin(5) ^ 2) + pkin(11) * t170 + t178; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t23 - mrSges(7,2) * t24 + t33; -t130; -t141 * t172 + t131; -pkin(11) * t172 + t131; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
