% Calculate joint inertia matrix for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:30:29
% EndTime: 2018-11-23 16:30:32
% DurationCPUTime: 2.17s
% Computational Cost: add. (4489->403), mult. (11885->559), div. (0->0), fcn. (13414->12), ass. (0->149)
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t103 = -mrSges(6,1) * t149 + mrSges(6,2) * t146;
t198 = -m(6) * pkin(4) + t103;
t197 = -mrSges(5,1) + t198;
t196 = 2 * pkin(10);
t195 = -2 * mrSges(7,3);
t147 = sin(qJ(4));
t150 = cos(qJ(4));
t140 = sin(pkin(12));
t142 = sin(pkin(6));
t143 = cos(pkin(12));
t145 = cos(pkin(6));
t151 = cos(qJ(3));
t144 = cos(pkin(7));
t148 = sin(qJ(3));
t171 = t144 * t148;
t141 = sin(pkin(7));
t174 = t141 * t148;
t66 = t145 * t174 + (t140 * t151 + t143 * t171) * t142;
t172 = t142 * t143;
t84 = -t141 * t172 + t144 * t145;
t54 = t147 * t84 + t150 * t66;
t160 = t144 * t172;
t173 = t141 * t151;
t175 = t140 * t142;
t65 = -t145 * t173 + t148 * t175 - t151 * t160;
t42 = -t146 * t54 + t149 * t65;
t43 = t146 * t65 + t149 * t54;
t53 = t147 * t66 - t150 * t84;
t11 = Ifges(7,5) * t43 + Ifges(7,6) * t42 + Ifges(7,3) * t53;
t12 = Ifges(6,5) * t43 + Ifges(6,6) * t42 + Ifges(6,3) * t53;
t194 = t11 + t12;
t102 = -t149 * mrSges(7,1) + t146 * mrSges(7,2);
t126 = -pkin(5) * t149 - pkin(4);
t193 = m(7) * t126 + t102;
t187 = pkin(1) * t145;
t86 = qJ(2) * t172 + t140 * t187;
t60 = (t141 * t145 + t160) * pkin(9) + t86;
t118 = t143 * t187;
t67 = pkin(2) * t145 + t118 + (-pkin(9) * t144 - qJ(2)) * t175;
t73 = (-pkin(9) * t140 * t141 - pkin(2) * t143 - pkin(1)) * t142;
t35 = -t148 * t60 + (t141 * t73 + t144 * t67) * t151;
t87 = -t150 * t144 + t147 * t174;
t192 = t87 ^ 2;
t191 = 2 * mrSges(3,1);
t190 = 0.2e1 * t145;
t188 = m(7) * pkin(5);
t32 = -pkin(3) * t84 - t35;
t19 = pkin(4) * t53 - pkin(11) * t54 + t32;
t48 = -t141 * t67 + t144 * t73;
t29 = pkin(3) * t65 - pkin(10) * t66 + t48;
t36 = t151 * t60 + t67 * t171 + t73 * t174;
t33 = pkin(10) * t84 + t36;
t10 = t147 * t29 + t150 * t33;
t8 = pkin(11) * t65 + t10;
t4 = t146 * t19 + t149 * t8;
t186 = pkin(10) * t150;
t185 = -Ifges(6,3) - Ifges(7,3);
t184 = -qJ(6) - pkin(11);
t21 = -mrSges(6,1) * t42 + mrSges(6,2) * t43;
t45 = mrSges(5,1) * t65 - mrSges(5,3) * t54;
t183 = t21 - t45;
t182 = Ifges(6,4) * t146;
t181 = Ifges(6,4) * t149;
t180 = Ifges(7,4) * t146;
t179 = Ifges(7,4) * t149;
t178 = t147 * t87;
t89 = t144 * t147 + t150 * t174;
t177 = t150 * t89;
t100 = -pkin(4) * t150 - pkin(11) * t147 - pkin(3);
t75 = t146 * t100 + t149 * t186;
t170 = t146 * t147;
t169 = t147 * t149;
t90 = mrSges(7,1) * t170 + mrSges(7,2) * t169;
t106 = Ifges(7,5) * t146 + Ifges(7,6) * t149;
t107 = Ifges(6,5) * t146 + Ifges(6,6) * t149;
t168 = Ifges(5,5) * t147 + Ifges(5,6) * t150;
t167 = t146 ^ 2 + t149 ^ 2;
t166 = Ifges(5,5) * t54 - Ifges(5,6) * t53 + Ifges(5,3) * t65;
t165 = Ifges(4,5) * t66 - Ifges(4,6) * t65 + Ifges(4,3) * t84;
t13 = Ifges(7,4) * t43 + Ifges(7,2) * t42 + Ifges(7,6) * t53;
t14 = Ifges(6,4) * t43 + Ifges(6,2) * t42 + Ifges(6,6) * t53;
t164 = -t13 / 0.2e1 - t14 / 0.2e1;
t15 = Ifges(7,1) * t43 + Ifges(7,4) * t42 + Ifges(7,5) * t53;
t16 = Ifges(6,1) * t43 + Ifges(6,4) * t42 + Ifges(6,5) * t53;
t163 = t15 / 0.2e1 + t16 / 0.2e1;
t79 = -Ifges(7,6) * t150 + (-Ifges(7,2) * t146 + t179) * t147;
t80 = -Ifges(6,6) * t150 + (-Ifges(6,2) * t146 + t181) * t147;
t162 = t79 / 0.2e1 + t80 / 0.2e1;
t81 = -Ifges(7,5) * t150 + (Ifges(7,1) * t149 - t180) * t147;
t82 = -Ifges(6,5) * t150 + (Ifges(6,1) * t149 - t182) * t147;
t161 = t81 / 0.2e1 + t82 / 0.2e1;
t159 = t106 / 0.2e1 + t107 / 0.2e1;
t108 = Ifges(7,2) * t149 + t180;
t109 = Ifges(6,2) * t149 + t182;
t158 = t108 / 0.2e1 + t109 / 0.2e1;
t111 = Ifges(7,1) * t146 + t179;
t112 = Ifges(6,1) * t146 + t181;
t157 = t111 / 0.2e1 + t112 / 0.2e1;
t20 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t3 = -t146 * t8 + t149 * t19;
t9 = -t147 * t33 + t150 * t29;
t156 = mrSges(6,1) * t146 + mrSges(6,2) * t149;
t7 = -pkin(4) * t65 - t9;
t153 = pkin(10) ^ 2;
t139 = t150 ^ 2;
t137 = t147 ^ 2;
t135 = t141 ^ 2;
t134 = t137 * t153;
t124 = t135 * t151 ^ 2;
t122 = Ifges(6,5) * t169;
t121 = Ifges(7,5) * t169;
t116 = mrSges(3,2) * t175;
t113 = Ifges(5,1) * t147 + Ifges(5,4) * t150;
t110 = Ifges(5,4) * t147 + Ifges(5,2) * t150;
t105 = t184 * t149;
t104 = -mrSges(5,1) * t150 + mrSges(5,2) * t147;
t101 = t184 * t146;
t99 = (pkin(5) * t146 + pkin(10)) * t147;
t98 = -mrSges(6,1) * t150 - mrSges(6,3) * t169;
t97 = -mrSges(7,1) * t150 - mrSges(7,3) * t169;
t96 = mrSges(6,2) * t150 - mrSges(6,3) * t170;
t95 = mrSges(7,2) * t150 - mrSges(7,3) * t170;
t93 = t149 * t100;
t91 = t156 * t147;
t85 = -qJ(2) * t175 + t118;
t78 = -Ifges(6,6) * t170 - Ifges(6,3) * t150 + t122;
t77 = -Ifges(7,6) * t170 - Ifges(7,3) * t150 + t121;
t74 = -t146 * t186 + t93;
t71 = -qJ(6) * t170 + t75;
t70 = -t146 * t173 + t149 * t89;
t69 = -t146 * t89 - t149 * t173;
t64 = -qJ(6) * t169 + t93 + (-pkin(10) * t146 - pkin(5)) * t150;
t56 = mrSges(4,1) * t84 - mrSges(4,3) * t66;
t55 = -mrSges(4,2) * t84 - mrSges(4,3) * t65;
t47 = mrSges(4,1) * t65 + mrSges(4,2) * t66;
t44 = -mrSges(5,2) * t65 - mrSges(5,3) * t53;
t34 = mrSges(5,1) * t53 + mrSges(5,2) * t54;
t27 = Ifges(5,1) * t54 - Ifges(5,4) * t53 + Ifges(5,5) * t65;
t26 = Ifges(5,4) * t54 - Ifges(5,2) * t53 + Ifges(5,6) * t65;
t25 = mrSges(6,1) * t53 - mrSges(6,3) * t43;
t24 = mrSges(7,1) * t53 - mrSges(7,3) * t43;
t23 = -mrSges(6,2) * t53 + mrSges(6,3) * t42;
t22 = -mrSges(7,2) * t53 + mrSges(7,3) * t42;
t5 = -pkin(5) * t42 + t7;
t2 = qJ(6) * t42 + t4;
t1 = pkin(5) * t53 - qJ(6) * t43 + t3;
t6 = [t66 * (Ifges(4,1) * t66 + Ifges(4,5) * t84) + (-0.2e1 * t86 * mrSges(3,2) + Ifges(3,3) * t145 + t191 * t85) * t145 + (-0.2e1 * pkin(1) * t116 + (-0.2e1 * t85 * mrSges(3,3) + Ifges(3,1) * t175 + Ifges(3,5) * t190) * t140 + (0.2e1 * t86 * mrSges(3,3) + Ifges(3,6) * t190 + (0.2e1 * Ifges(3,4) * t140 + Ifges(3,2) * t143 + pkin(1) * t191) * t142) * t143) * t142 + t84 * t165 + (-0.2e1 * Ifges(4,4) * t66 + Ifges(4,2) * t65 - Ifges(4,6) * t84 + t166) * t65 + (-t26 + t194) * t53 + (t13 + t14) * t42 + Ifges(2,3) + (t15 + t16) * t43 + 0.2e1 * t5 * t20 + 0.2e1 * t7 * t21 + 0.2e1 * t2 * t22 + 0.2e1 * t4 * t23 + 0.2e1 * t1 * t24 + 0.2e1 * t3 * t25 + 0.2e1 * t32 * t34 + 0.2e1 * t10 * t44 + 0.2e1 * t9 * t45 + 0.2e1 * t48 * t47 + t54 * t27 + 0.2e1 * t36 * t55 + 0.2e1 * t35 * t56 + m(3) * (pkin(1) ^ 2 * t142 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) + m(5) * (t10 ^ 2 + t32 ^ 2 + t9 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2 + t48 ^ 2); t144 * t47 + t89 * t44 + t116 + (t22 + t23) * t70 + (t24 + t25) * t69 + (-m(3) * pkin(1) - mrSges(3,1) * t143) * t142 + (t20 + t183) * t87 + (t148 * t55 + (-t34 + t56) * t151) * t141 + m(7) * (t1 * t69 + t2 * t70 + t5 * t87) + m(6) * (t3 * t69 + t4 * t70 + t7 * t87) + m(5) * (t10 * t89 - t173 * t32 - t87 * t9) + m(4) * (t144 * t48 + (t148 * t36 + t151 * t35) * t141); m(3) + m(5) * (t89 ^ 2 + t124 + t192) + m(4) * (t135 * t148 ^ 2 + t144 ^ 2 + t124) + 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t69 ^ 2 + t70 ^ 2 + t192); (t77 / 0.2e1 + t78 / 0.2e1 - t110 / 0.2e1) * t53 + (-t9 * mrSges(5,3) + t27 / 0.2e1 + t163 * t149 + t164 * t146 + t183 * pkin(10)) * t147 + t65 * t168 / 0.2e1 + t161 * t43 + t162 * t42 + (t10 * mrSges(5,3) + pkin(10) * t44 - t11 / 0.2e1 - t12 / 0.2e1 + t26 / 0.2e1) * t150 + m(7) * (t1 * t64 + t2 * t71 + t5 * t99) + m(6) * (pkin(10) * t147 * t7 + t3 * t74 + t4 * t75) + t165 + m(5) * (-pkin(3) * t32 + (t10 * t150 - t9 * t147) * pkin(10)) - pkin(3) * t34 + t35 * mrSges(4,1) - t36 * mrSges(4,2) + t64 * t24 + t71 * t22 + t74 * t25 + t75 * t23 + t5 * t90 + t7 * t91 + t2 * t95 + t4 * t96 + t1 * t97 + t3 * t98 + t99 * t20 + t32 * t104 + t54 * t113 / 0.2e1; mrSges(5,3) * t177 + (t95 + t96) * t70 + (t97 + t98) * t69 + (t147 * mrSges(5,3) + t90 + t91) * t87 + (-t148 * mrSges(4,2) + (mrSges(4,1) - t104) * t151) * t141 + m(7) * (t64 * t69 + t70 * t71 + t87 * t99) + m(6) * (pkin(10) * t178 + t69 * t74 + t70 * t75) + m(5) * (pkin(3) * t173 + (t177 + t178) * pkin(10)); -0.2e1 * pkin(3) * t104 + 0.2e1 * t64 * t97 + 0.2e1 * t71 * t95 + 0.2e1 * t74 * t98 + 0.2e1 * t75 * t96 + 0.2e1 * t99 * t90 + Ifges(4,3) + (t137 + t139) * mrSges(5,3) * t196 + (-t78 - t77 + t110) * t150 + m(7) * (t64 ^ 2 + t71 ^ 2 + t99 ^ 2) + m(6) * (t74 ^ 2 + t75 ^ 2 + t134) + m(5) * (pkin(3) ^ 2 + t139 * t153 + t134) + (t91 * t196 + t113 + (t81 + t82) * t149 + (-t79 - t80) * t146) * t147; t9 * mrSges(5,1) - t10 * mrSges(5,2) - pkin(4) * t21 + t101 * t24 + t5 * t102 - t105 * t22 + t126 * t20 + t159 * t53 + t157 * t43 + t158 * t42 + m(7) * (t1 * t101 - t105 * t2 + t126 * t5) + (t2 * mrSges(7,3) + t4 * mrSges(6,3) + (m(6) * t4 + t23) * pkin(11) - t164) * t149 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + (-m(6) * t3 - t25) * pkin(11) + t163) * t146 + t166 + t198 * t7; -t89 * mrSges(5,2) + m(7) * (t101 * t69 - t105 * t70) + (t193 + t197) * t87 + (m(6) * pkin(11) + mrSges(6,3) + mrSges(7,3)) * (-t146 * t69 + t149 * t70); m(7) * (t101 * t64 - t105 * t71 + t126 * t99) - pkin(4) * t91 + t101 * t97 + t99 * t102 - t105 * t95 + t126 * t90 - t159 * t150 + (-t150 * mrSges(5,2) + t197 * t147) * pkin(10) + (t71 * mrSges(7,3) + t75 * mrSges(6,3) + t157 * t147 + (m(6) * t75 + t96) * pkin(11) + t162) * t149 + (-t64 * mrSges(7,3) - t74 * mrSges(6,3) - t158 * t147 + (-m(6) * t74 - t98) * pkin(11) + t161) * t146 + t168; -0.2e1 * pkin(4) * t103 + 0.2e1 * t126 * t102 + Ifges(5,3) + 0.2e1 * t167 * pkin(11) * mrSges(6,3) + m(7) * (t101 ^ 2 + t105 ^ 2 + t126 ^ 2) + m(6) * (pkin(11) ^ 2 * t167 + pkin(4) ^ 2) + (t105 * t195 + t108 + t109) * t149 + (t101 * t195 + t111 + t112) * t146; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t24) * pkin(5) + t194; (-mrSges(6,2) - mrSges(7,2)) * t70 + (mrSges(6,1) + mrSges(7,1) + t188) * t69; mrSges(6,1) * t74 + mrSges(7,1) * t64 - mrSges(6,2) * t75 - mrSges(7,2) * t71 + t121 + t122 + t185 * t150 + (-Ifges(6,6) - Ifges(7,6)) * t170 + (m(7) * t64 + t97) * pkin(5); mrSges(7,1) * t101 + mrSges(7,2) * t105 - t156 * pkin(11) + (m(7) * t101 - t146 * mrSges(7,3)) * pkin(5) + t107 + t106; (0.2e1 * mrSges(7,1) + t188) * pkin(5) - t185; m(7) * t5 + t20; m(7) * t87; m(7) * t99 + t90; t193; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
