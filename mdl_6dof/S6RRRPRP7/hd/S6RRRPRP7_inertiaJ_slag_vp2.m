% Calculate joint inertia matrix for
% S6RRRPRP7
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:59
% EndTime: 2019-03-09 17:04:03
% DurationCPUTime: 2.04s
% Computational Cost: add. (3389->397), mult. (7520->555), div. (0->0), fcn. (8342->10), ass. (0->144)
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t204 = t149 ^ 2 + t152 ^ 2;
t203 = 0.2e1 * t204;
t202 = m(5) * pkin(3);
t146 = sin(pkin(6));
t154 = cos(qJ(2));
t177 = t146 * t154;
t145 = sin(pkin(11));
t147 = cos(pkin(11));
t148 = cos(pkin(6));
t150 = sin(qJ(3));
t153 = cos(qJ(3));
t151 = sin(qJ(2));
t178 = t146 * t151;
t96 = t148 * t153 - t150 * t178;
t97 = t148 * t150 + t153 * t178;
t63 = t145 * t96 + t147 * t97;
t43 = t149 * t63 + t152 * t177;
t44 = -t149 * t177 + t152 * t63;
t62 = t145 * t97 - t147 * t96;
t8 = Ifges(6,5) * t44 - Ifges(6,6) * t43 + Ifges(6,3) * t62;
t9 = Ifges(7,4) * t44 + Ifges(7,2) * t62 + Ifges(7,6) * t43;
t201 = t8 + t9;
t106 = t145 * t150 - t147 * t153;
t107 = t145 * t153 + t147 * t150;
t179 = t107 * t152;
t180 = t107 * t149;
t46 = Ifges(6,5) * t179 - Ifges(6,6) * t180 + Ifges(6,3) * t106;
t47 = Ifges(7,4) * t179 + Ifges(7,2) * t106 + Ifges(7,6) * t180;
t200 = t46 + t47;
t199 = Ifges(4,5) * t150 + Ifges(5,5) * t107 + Ifges(4,6) * t153 - Ifges(5,6) * t106;
t189 = -qJ(4) - pkin(9);
t112 = t189 * t153;
t163 = t189 * t150;
t78 = -t112 * t145 - t147 * t163;
t198 = t78 ^ 2;
t197 = 0.2e1 * t78;
t99 = t148 * t151 * pkin(1) + pkin(8) * t177;
t87 = pkin(9) * t148 + t99;
t88 = (-pkin(2) * t154 - pkin(9) * t151 - pkin(1)) * t146;
t53 = -t150 * t87 + t153 * t88;
t31 = -pkin(3) * t177 - qJ(4) * t97 + t53;
t54 = t150 * t88 + t153 * t87;
t36 = qJ(4) * t96 + t54;
t16 = t145 * t31 + t147 * t36;
t14 = -pkin(10) * t177 + t16;
t126 = pkin(8) * t178;
t195 = pkin(1) * t154;
t86 = t126 + (-pkin(2) - t195) * t148;
t65 = -pkin(3) * t96 + t86;
t18 = pkin(4) * t62 - pkin(10) * t63 + t65;
t4 = t152 * t14 + t149 * t18;
t194 = pkin(3) * t145;
t193 = pkin(3) * t147;
t98 = t148 * t195 - t126;
t192 = t98 * mrSges(3,1);
t191 = t99 * mrSges(3,2);
t190 = Ifges(4,3) + Ifges(5,3);
t21 = -mrSges(7,2) * t43 + mrSges(7,3) * t62;
t22 = -mrSges(6,2) * t62 - mrSges(6,3) * t43;
t188 = t21 + t22;
t23 = mrSges(6,1) * t62 - mrSges(6,3) * t44;
t24 = -t62 * mrSges(7,1) + t44 * mrSges(7,2);
t187 = t23 - t24;
t70 = -mrSges(6,2) * t106 - mrSges(6,3) * t180;
t73 = -mrSges(7,2) * t180 + mrSges(7,3) * t106;
t186 = t70 + t73;
t71 = mrSges(6,1) * t106 - mrSges(6,3) * t179;
t72 = -t106 * mrSges(7,1) + mrSges(7,2) * t179;
t185 = t71 - t72;
t134 = -pkin(3) * t153 - pkin(2);
t69 = pkin(4) * t106 - pkin(10) * t107 + t134;
t80 = -t147 * t112 + t145 * t163;
t33 = t149 * t69 + t152 * t80;
t184 = Ifges(6,4) * t149;
t183 = Ifges(6,4) * t152;
t182 = Ifges(7,5) * t149;
t181 = Ifges(7,5) * t152;
t132 = pkin(10) + t194;
t175 = t204 * t132 ^ 2;
t114 = Ifges(6,5) * t149 + Ifges(6,6) * t152;
t172 = t150 ^ 2 + t153 ^ 2;
t10 = Ifges(6,4) * t44 - Ifges(6,2) * t43 + Ifges(6,6) * t62;
t7 = Ifges(7,5) * t44 + Ifges(7,6) * t62 + Ifges(7,3) * t43;
t171 = t7 / 0.2e1 - t10 / 0.2e1;
t11 = Ifges(7,1) * t44 + Ifges(7,4) * t62 + Ifges(7,5) * t43;
t12 = Ifges(6,1) * t44 - Ifges(6,4) * t43 + Ifges(6,5) * t62;
t170 = t11 / 0.2e1 + t12 / 0.2e1;
t45 = Ifges(7,6) * t106 + (Ifges(7,3) * t149 + t181) * t107;
t48 = Ifges(6,6) * t106 + (-Ifges(6,2) * t149 + t183) * t107;
t169 = t45 / 0.2e1 - t48 / 0.2e1;
t49 = Ifges(7,4) * t106 + (Ifges(7,1) * t152 + t182) * t107;
t50 = Ifges(6,5) * t106 + (Ifges(6,1) * t152 - t184) * t107;
t168 = t49 / 0.2e1 + t50 / 0.2e1;
t167 = Ifges(3,5) * t178 + Ifges(3,6) * t177 + Ifges(3,3) * t148;
t133 = -pkin(4) - t193;
t113 = -Ifges(7,3) * t152 + t182;
t116 = Ifges(6,2) * t152 + t184;
t166 = t113 / 0.2e1 - t116 / 0.2e1;
t115 = Ifges(7,4) * t149 - Ifges(7,6) * t152;
t165 = t114 / 0.2e1 + t115 / 0.2e1;
t118 = Ifges(7,1) * t149 - t181;
t119 = Ifges(6,1) * t149 + t183;
t164 = t118 / 0.2e1 + t119 / 0.2e1;
t29 = t62 * mrSges(5,1) + t63 * mrSges(5,2);
t15 = -t145 * t36 + t147 * t31;
t162 = -Ifges(4,5) * t97 - Ifges(5,5) * t63 - Ifges(4,6) * t96 + Ifges(5,6) * t62;
t74 = t106 * mrSges(5,1) + t107 * mrSges(5,2);
t13 = pkin(4) * t177 - t15;
t110 = -t152 * mrSges(6,1) + t149 * mrSges(6,2);
t160 = t149 * mrSges(6,1) + t152 * mrSges(6,2);
t109 = -t152 * mrSges(7,1) - t149 * mrSges(7,3);
t159 = t149 * mrSges(7,1) - t152 * mrSges(7,3);
t158 = pkin(5) * t152 + qJ(6) * t149;
t157 = pkin(5) * t149 - qJ(6) * t152;
t3 = -t14 * t149 + t152 * t18;
t32 = -t149 * t80 + t152 * t69;
t120 = Ifges(4,1) * t150 + Ifges(4,4) * t153;
t117 = Ifges(4,4) * t150 + Ifges(4,2) * t153;
t111 = -mrSges(4,1) * t153 + mrSges(4,2) * t150;
t105 = t133 - t158;
t82 = -mrSges(4,1) * t177 - mrSges(4,3) * t97;
t81 = mrSges(4,2) * t177 + mrSges(4,3) * t96;
t76 = Ifges(5,1) * t107 - Ifges(5,4) * t106;
t75 = Ifges(5,4) * t107 - Ifges(5,2) * t106;
t68 = -mrSges(4,1) * t96 + mrSges(4,2) * t97;
t67 = t160 * t107;
t66 = t159 * t107;
t56 = Ifges(4,1) * t97 + Ifges(4,4) * t96 - Ifges(4,5) * t177;
t55 = Ifges(4,4) * t97 + Ifges(4,2) * t96 - Ifges(4,6) * t177;
t52 = -mrSges(5,1) * t177 - mrSges(5,3) * t63;
t51 = mrSges(5,2) * t177 - mrSges(5,3) * t62;
t37 = t107 * t157 + t78;
t28 = -pkin(5) * t106 - t32;
t27 = qJ(6) * t106 + t33;
t26 = Ifges(5,1) * t63 - Ifges(5,4) * t62 - Ifges(5,5) * t177;
t25 = Ifges(5,4) * t63 - Ifges(5,2) * t62 - Ifges(5,6) * t177;
t20 = mrSges(6,1) * t43 + mrSges(6,2) * t44;
t19 = mrSges(7,1) * t43 - mrSges(7,3) * t44;
t5 = pkin(5) * t43 - qJ(6) * t44 + t13;
t2 = -pkin(5) * t62 - t3;
t1 = qJ(6) * t62 + t4;
t6 = [(-t25 + t201) * t62 + ((-0.2e1 * t98 * mrSges(3,3) + Ifges(3,5) * t148 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t151) * t146) * t151 + (0.2e1 * t99 * mrSges(3,3) + Ifges(3,6) * t148 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t151 + (Ifges(3,2) + t190) * t154) * t146 + t162) * t154) * t146 + (t167 - 0.2e1 * t191 + 0.2e1 * t192) * t148 + t96 * t55 + t97 * t56 + 0.2e1 * t54 * t81 + 0.2e1 * t53 * t82 + 0.2e1 * t86 * t68 + t63 * t26 + 0.2e1 * t65 * t29 + 0.2e1 * t16 * t51 + 0.2e1 * t15 * t52 + 0.2e1 * t2 * t24 + 0.2e1 * t5 * t19 + 0.2e1 * t13 * t20 + 0.2e1 * t1 * t21 + 0.2e1 * t4 * t22 + 0.2e1 * t3 * t23 + (t11 + t12) * t44 + (t7 - t10) * t43 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t65 ^ 2) + m(4) * (t53 ^ 2 + t54 ^ 2 + t86 ^ 2) + m(3) * (pkin(1) ^ 2 * t146 ^ 2 + t98 ^ 2 + t99 ^ 2) + Ifges(2,3); -t191 + t192 + t167 - t199 * t177 / 0.2e1 + t168 * t44 + t169 * t43 + (t26 / 0.2e1 - t15 * mrSges(5,3) + t170 * t152 + t171 * t149) * t107 + t80 * t51 + t5 * t66 + t13 * t67 - pkin(2) * t68 + t4 * t70 + t3 * t71 + t2 * t72 + t1 * t73 + t65 * t74 + t63 * t76 / 0.2e1 + t37 * t19 + t27 * t21 + t28 * t24 + t32 * t23 + t33 * t22 + m(4) * (-pkin(2) * t86 + (-t53 * t150 + t54 * t153) * pkin(9)) + (t8 / 0.2e1 + t9 / 0.2e1 - t25 / 0.2e1 - t16 * mrSges(5,3)) * t106 + (t20 - t52) * t78 + (t55 / 0.2e1 + t54 * mrSges(4,3) + pkin(9) * t81) * t153 + (t56 / 0.2e1 - t53 * mrSges(4,3) - pkin(9) * t82) * t150 + m(7) * (t1 * t27 + t2 * t28 + t37 * t5) + m(6) * (t13 * t78 + t3 * t32 + t33 * t4) + m(5) * (t134 * t65 - t15 * t78 + t16 * t80) + (-t75 / 0.2e1 + t46 / 0.2e1 + t47 / 0.2e1) * t62 + t86 * t111 + t96 * t117 / 0.2e1 + t97 * t120 / 0.2e1 + t134 * t29; -0.2e1 * pkin(2) * t111 + t153 * t117 + t150 * t120 + 0.2e1 * t134 * t74 + 0.2e1 * t27 * t73 + 0.2e1 * t28 * t72 + 0.2e1 * t32 * t71 + 0.2e1 * t33 * t70 + 0.2e1 * t37 * t66 + t67 * t197 + Ifges(3,3) + 0.2e1 * t172 * pkin(9) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t80 + t200 - t75) * t106 + m(7) * (t27 ^ 2 + t28 ^ 2 + t37 ^ 2) + m(6) * (t32 ^ 2 + t33 ^ 2 + t198) + m(5) * (t134 ^ 2 + t80 ^ 2 + t198) + m(4) * (pkin(9) ^ 2 * t172 + pkin(2) ^ 2) + (mrSges(5,3) * t197 + t76 + (t49 + t50) * t152 + (t45 - t48) * t149) * t107; -t162 + (t145 * t51 + t147 * t52 + m(5) * (t145 * t16 + t147 * t15)) * pkin(3) + t105 * t19 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + t15 * mrSges(5,1) - t16 * mrSges(5,2) + t165 * t62 + m(6) * (t13 * t133 + (-t149 * t3 + t152 * t4) * t132) + t164 * t44 + m(7) * (t105 * t5 + (t1 * t152 + t149 * t2) * t132) + t166 * t43 - t190 * t177 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) - t132 * t187 + t170) * t149 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) + t132 * t188 - t171) * t152 + t5 * t109 + t13 * t110 + t133 * t20; -t80 * mrSges(5,2) + t105 * t66 + t37 * t109 + t133 * t67 + (-mrSges(5,1) + t110) * t78 + (-mrSges(4,1) * t150 - mrSges(4,2) * t153) * pkin(9) + (-mrSges(5,3) * t194 + t165) * t106 + (t27 * mrSges(7,2) + t33 * mrSges(6,3) + t132 * t186 - t169) * t152 + (t28 * mrSges(7,2) - t32 * mrSges(6,3) - t132 * t185 + t168) * t149 + m(7) * (t105 * t37 + (t149 * t28 + t152 * t27) * t132) + m(6) * (t133 * t78 + (-t149 * t32 + t152 * t33) * t132) + (t145 * t80 - t147 * t78) * t202 + (-mrSges(5,3) * t193 + t149 * t166 + t152 * t164) * t107 + t199; 0.2e1 * t105 * t109 + 0.2e1 * t133 * t110 + (-t113 + t116) * t152 + (t118 + t119) * t149 + m(7) * (t105 ^ 2 + t175) + m(6) * (t133 ^ 2 + t175) + t190 + (mrSges(7,2) + mrSges(6,3)) * t132 * t203 + (0.2e1 * mrSges(5,1) * t147 - 0.2e1 * mrSges(5,2) * t145 + (t145 ^ 2 + t147 ^ 2) * t202) * pkin(3); t187 * t152 + t188 * t149 + m(6) * (t149 * t4 + t152 * t3) + m(7) * (t1 * t149 - t152 * t2) + m(5) * t65 + t29; t185 * t152 + t186 * t149 + m(6) * (t149 * t33 + t152 * t32) + m(7) * (t149 * t27 - t152 * t28) + m(5) * t134 + t74; 0; m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t203; -pkin(5) * t24 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t21 + t1 * mrSges(7,3) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t201; -pkin(5) * t72 + m(7) * (-pkin(5) * t28 + qJ(6) * t27) + qJ(6) * t73 + t27 * mrSges(7,3) - t33 * mrSges(6,2) + t32 * mrSges(6,1) - t28 * mrSges(7,1) + t200; -t157 * mrSges(7,2) + (-m(7) * t157 - t159 - t160) * t132 + t115 + t114; m(7) * t158 - t109 - t110; Ifges(7,2) + Ifges(6,3) + 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2); m(7) * t2 + t24; m(7) * t28 + t72; (m(7) * t132 + mrSges(7,2)) * t149; -m(7) * t152; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
