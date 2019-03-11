% Calculate joint inertia matrix for
% S6RRPRRR5
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:07
% EndTime: 2019-03-09 13:38:12
% DurationCPUTime: 2.04s
% Computational Cost: add. (4194->406), mult. (9974->585), div. (0->0), fcn. (11063->12), ass. (0->160)
t207 = Ifges(3,3) + Ifges(4,3);
t150 = sin(qJ(4));
t149 = sin(qJ(5));
t153 = cos(qJ(5));
t165 = mrSges(6,1) * t149 + mrSges(6,2) * t153;
t102 = t165 * t150;
t148 = sin(qJ(6));
t152 = cos(qJ(6));
t107 = t148 * t153 + t149 * t152;
t95 = t107 * t150;
t106 = -t148 * t149 + t152 * t153;
t96 = t106 * t150;
t57 = t95 * mrSges(7,1) + t96 * mrSges(7,2);
t206 = -t102 - t57;
t144 = sin(pkin(12));
t131 = pkin(2) * t144 + pkin(9);
t205 = 0.2e1 * t131;
t146 = cos(pkin(12));
t145 = sin(pkin(6));
t155 = cos(qJ(2));
t176 = t145 * t155;
t151 = sin(qJ(2));
t177 = t145 * t151;
t90 = t144 * t177 - t146 * t176;
t91 = (t144 * t155 + t146 * t151) * t145;
t204 = Ifges(4,5) * t91 - Ifges(4,6) * t90;
t108 = (-pkin(2) * t155 - pkin(1)) * t145;
t203 = 0.2e1 * t108;
t147 = cos(pkin(6));
t154 = cos(qJ(4));
t67 = t147 * t150 + t154 * t91;
t48 = -t149 * t67 + t153 * t90;
t49 = t149 * t90 + t153 * t67;
t66 = -t147 * t154 + t150 * t91;
t17 = Ifges(6,1) * t49 + Ifges(6,4) * t48 + Ifges(6,5) * t66;
t202 = t17 / 0.2e1;
t28 = -t148 * t49 + t152 * t48;
t201 = t28 / 0.2e1;
t29 = t148 * t48 + t152 * t49;
t200 = t29 / 0.2e1;
t184 = Ifges(6,4) * t153;
t93 = -Ifges(6,6) * t154 + (-Ifges(6,2) * t149 + t184) * t150;
t199 = t93 / 0.2e1;
t185 = Ifges(6,4) * t149;
t94 = -Ifges(6,5) * t154 + (Ifges(6,1) * t153 - t185) * t150;
t198 = t94 / 0.2e1;
t197 = -t95 / 0.2e1;
t196 = t96 / 0.2e1;
t195 = -pkin(11) - pkin(10);
t194 = t106 / 0.2e1;
t193 = t107 / 0.2e1;
t117 = Ifges(6,1) * t149 + t184;
t192 = t117 / 0.2e1;
t191 = pkin(1) * t147;
t190 = pkin(4) * t154;
t189 = pkin(10) * t150;
t188 = -Ifges(7,3) - Ifges(6,3);
t127 = t155 * t191;
t73 = pkin(2) * t147 + t127 + (-pkin(8) - qJ(3)) * t177;
t101 = pkin(8) * t176 + t151 * t191;
t82 = qJ(3) * t176 + t101;
t46 = t144 * t73 + t146 * t82;
t42 = pkin(9) * t147 + t46;
t50 = pkin(3) * t90 - pkin(9) * t91 + t108;
t25 = t150 * t50 + t154 * t42;
t20 = pkin(10) * t90 + t25;
t45 = -t144 * t82 + t146 * t73;
t41 = -pkin(3) * t147 - t45;
t23 = pkin(4) * t66 - pkin(10) * t67 + t41;
t7 = t149 * t23 + t153 * t20;
t30 = -mrSges(6,1) * t48 + mrSges(6,2) * t49;
t52 = mrSges(5,1) * t90 - mrSges(5,3) * t67;
t187 = t30 - t52;
t186 = Ifges(7,5) * t96 - Ifges(7,6) * t95;
t183 = Ifges(7,3) * t154;
t100 = -pkin(8) * t177 + t127;
t182 = t100 * mrSges(3,1);
t181 = t101 * mrSges(3,2);
t112 = -mrSges(6,1) * t153 + mrSges(6,2) * t149;
t180 = -mrSges(5,1) + t112;
t132 = -pkin(2) * t146 - pkin(3);
t103 = t132 - t189 - t190;
t178 = t131 * t154;
t65 = t149 * t103 + t153 * t178;
t179 = t131 * t150;
t175 = t149 * t150;
t174 = t150 * t153;
t69 = Ifges(7,5) * t107 + Ifges(7,6) * t106;
t114 = Ifges(6,5) * t149 + Ifges(6,6) * t153;
t173 = Ifges(5,5) * t150 + Ifges(5,6) * t154;
t172 = t149 ^ 2 + t153 ^ 2;
t141 = t150 ^ 2;
t143 = t154 ^ 2;
t171 = t141 + t143;
t8 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t66;
t15 = Ifges(6,5) * t49 + Ifges(6,6) * t48 + Ifges(6,3) * t66;
t170 = Ifges(5,5) * t67 - Ifges(5,6) * t66 + Ifges(5,3) * t90;
t169 = t69 / 0.2e1 + t114 / 0.2e1;
t6 = -t149 * t20 + t153 * t23;
t24 = -t150 * t42 + t154 * t50;
t168 = t172 * mrSges(6,3);
t167 = Ifges(6,5) * t174 - Ifges(6,6) * t175;
t166 = -t149 * t6 + t153 * t7;
t98 = t153 * t103;
t64 = -t149 * t178 + t98;
t164 = -t149 * t64 + t153 * t65;
t56 = -pkin(11) * t174 + t98 + (-t131 * t149 - pkin(5)) * t154;
t58 = -pkin(11) * t175 + t65;
t34 = -t148 * t58 + t152 * t56;
t35 = t148 * t56 + t152 * t58;
t163 = t34 * mrSges(7,1) - t35 * mrSges(7,2) + t186;
t121 = t195 * t149;
t122 = t195 * t153;
t78 = t121 * t152 + t122 * t148;
t79 = t121 * t148 - t122 * t152;
t162 = t78 * mrSges(7,1) - t79 * mrSges(7,2) + t69;
t19 = -pkin(4) * t90 - t24;
t4 = pkin(5) * t66 - pkin(11) * t49 + t6;
t5 = pkin(11) * t48 + t7;
t2 = -t148 * t5 + t152 * t4;
t3 = t148 * t4 + t152 * t5;
t161 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t160 = Ifges(3,5) * t177 + Ifges(3,6) * t176 + t207 * t147 + t204;
t159 = (mrSges(7,1) * t152 - mrSges(7,2) * t148) * pkin(5);
t133 = -pkin(5) * t153 - pkin(4);
t129 = t131 ^ 2;
t119 = t141 * t129;
t118 = Ifges(5,1) * t150 + Ifges(5,4) * t154;
t116 = Ifges(5,4) * t150 + Ifges(5,2) * t154;
t115 = Ifges(6,2) * t153 + t185;
t113 = -t154 * mrSges(5,1) + t150 * mrSges(5,2);
t110 = -mrSges(6,1) * t154 - mrSges(6,3) * t174;
t109 = mrSges(6,2) * t154 - mrSges(6,3) * t175;
t99 = (pkin(5) * t149 + t131) * t150;
t92 = -Ifges(6,3) * t154 + t167;
t83 = t91 * mrSges(4,2);
t81 = -mrSges(7,1) * t154 - mrSges(7,3) * t96;
t80 = mrSges(7,2) * t154 - mrSges(7,3) * t95;
t75 = mrSges(4,1) * t147 - mrSges(4,3) * t91;
t74 = -mrSges(4,2) * t147 - mrSges(4,3) * t90;
t71 = Ifges(7,1) * t107 + Ifges(7,4) * t106;
t70 = Ifges(7,4) * t107 + Ifges(7,2) * t106;
t68 = -mrSges(7,1) * t106 + mrSges(7,2) * t107;
t55 = Ifges(7,1) * t96 - Ifges(7,4) * t95 - Ifges(7,5) * t154;
t54 = Ifges(7,4) * t96 - Ifges(7,2) * t95 - Ifges(7,6) * t154;
t53 = -t183 + t186;
t51 = -mrSges(5,2) * t90 - mrSges(5,3) * t66;
t40 = mrSges(5,1) * t66 + mrSges(5,2) * t67;
t37 = Ifges(5,1) * t67 - Ifges(5,4) * t66 + Ifges(5,5) * t90;
t36 = Ifges(5,4) * t67 - Ifges(5,2) * t66 + Ifges(5,6) * t90;
t32 = mrSges(6,1) * t66 - mrSges(6,3) * t49;
t31 = -mrSges(6,2) * t66 + mrSges(6,3) * t48;
t16 = Ifges(6,4) * t49 + Ifges(6,2) * t48 + Ifges(6,6) * t66;
t14 = mrSges(7,1) * t66 - mrSges(7,3) * t29;
t13 = -mrSges(7,2) * t66 + mrSges(7,3) * t28;
t12 = -pkin(5) * t48 + t19;
t11 = -mrSges(7,1) * t28 + mrSges(7,2) * t29;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t66;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t66;
t1 = [0.2e1 * t12 * t11 + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t14 + Ifges(4,1) * t91 ^ 2 + m(3) * (t100 ^ 2 + t101 ^ 2) + (t160 - 0.2e1 * t181 + 0.2e1 * t182 + t204) * t147 + ((Ifges(3,5) * t151 + Ifges(3,6) * t155) * t147 + 0.2e1 * (-t100 * t151 + t101 * t155) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t155 + mrSges(3,2) * t151) + t151 * (Ifges(3,1) * t151 + Ifges(3,4) * t155) + t155 * (t151 * Ifges(3,4) + Ifges(3,2) * t155) + m(3) * pkin(1) ^ 2) * t145) * t145 + t28 * t9 + t29 * t10 + 0.2e1 * t19 * t30 + 0.2e1 * t7 * t31 + 0.2e1 * t6 * t32 + 0.2e1 * t41 * t40 + t48 * t16 + t49 * t17 + 0.2e1 * t25 * t51 + 0.2e1 * t24 * t52 + Ifges(2,3) + t67 * t37 + 0.2e1 * t46 * t74 + 0.2e1 * t45 * t75 + (mrSges(4,1) * t203 - 0.2e1 * Ifges(4,4) * t91 + Ifges(4,2) * t90 + t170) * t90 + t83 * t203 + (t15 + t8 - t36) * t66 + m(7) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t19 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t41 ^ 2) + m(4) * (t108 ^ 2 + t45 ^ 2 + t46 ^ 2); m(5) * (t132 * t41 + (-t24 * t150 + t25 * t154) * t131) + (t25 * mrSges(5,3) + t131 * t51 - t15 / 0.2e1 - t8 / 0.2e1 + t36 / 0.2e1) * t154 - t181 + t182 + t90 * t173 / 0.2e1 + t34 * t14 + t35 * t13 + m(6) * (t179 * t19 + t6 * t64 + t65 * t7) + t45 * mrSges(4,1) - t46 * mrSges(4,2) + t12 * t57 + t64 * t32 + t65 * t31 + t160 + t3 * t80 + t2 * t81 + t99 * t11 + t19 * t102 + (-t24 * mrSges(5,3) - t149 * t16 / 0.2e1 + t153 * t202 + t37 / 0.2e1 + t187 * t131) * t150 + t7 * t109 + t6 * t110 + t41 * t113 + t67 * t118 / 0.2e1 + t132 * t40 + t10 * t196 + t9 * t197 + t49 * t198 + t48 * t199 + t55 * t200 + t54 * t201 + (m(4) * (t144 * t46 + t146 * t45) + t144 * t74 + t146 * t75) * pkin(2) + (t53 / 0.2e1 + t92 / 0.2e1 - t116 / 0.2e1) * t66 + m(7) * (t12 * t99 + t2 * t34 + t3 * t35); 0.2e1 * t65 * t109 + 0.2e1 * t64 * t110 + 0.2e1 * t132 * t113 + 0.2e1 * t34 * t81 + 0.2e1 * t35 * t80 - t95 * t54 + t96 * t55 + 0.2e1 * t99 * t57 + (-t53 - t92 + t116) * t154 + (t102 * t205 - t149 * t93 + t153 * t94 + t118) * t150 + m(7) * (t34 ^ 2 + t35 ^ 2 + t99 ^ 2) + m(6) * (t64 ^ 2 + t65 ^ 2 + t119) + m(5) * (t129 * t143 + t132 ^ 2 + t119) + m(4) * (t144 ^ 2 + t146 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t146 - mrSges(4,2) * t144) * pkin(2) + t171 * mrSges(5,3) * t205 + t207; t90 * mrSges(4,1) + t96 * t13 - t95 * t14 + t83 + (-t11 - t187) * t154 + (-t149 * t32 + t153 * t31 + t51) * t150 + m(7) * (-t12 * t154 - t2 * t95 + t3 * t96) + m(6) * (t150 * t166 - t154 * t19) + m(5) * (t150 * t25 + t154 * t24) + m(4) * t108; m(7) * (-t34 * t95 + t35 * t96) + t96 * t80 - t95 * t81 + (m(6) * (t164 - t178) + t153 * t109 - t149 * t110) * t150 + (-m(7) * t99 + t206) * t154; m(4) + m(5) * t171 + m(6) * (t141 * t172 + t143) + m(7) * (t95 ^ 2 + t96 ^ 2 + t143); t24 * mrSges(5,1) - t25 * mrSges(5,2) + (t7 * mrSges(6,3) + pkin(10) * t31 + t16 / 0.2e1) * t153 + (-t6 * mrSges(6,3) - pkin(10) * t32 + t202) * t149 + m(7) * (t12 * t133 + t2 * t78 + t3 * t79) + t169 * t66 + m(6) * (-pkin(4) * t19 + pkin(10) * t166) + (t106 * t3 - t107 * t2) * mrSges(7,3) + t170 - pkin(4) * t30 + t12 * t68 + t70 * t201 + t71 * t200 + t78 * t14 + t79 * t13 + t9 * t194 + t10 * t193 + t19 * t112 + t48 * t115 / 0.2e1 + t49 * t192 + t133 * t11; t79 * t80 + t78 * t81 + t70 * t197 + t71 * t196 + t99 * t68 - pkin(4) * t102 + t54 * t194 + t55 * t193 + t133 * t57 + t180 * t179 + (t106 * t35 - t107 * t34) * mrSges(7,3) + (t65 * mrSges(6,3) + pkin(10) * t109 + t150 * t192 + t199) * t153 + (-t64 * mrSges(6,3) - pkin(10) * t110 - t150 * t115 / 0.2e1 + t198) * t149 + m(7) * (t133 * t99 + t34 * t78 + t35 * t79) + m(6) * (-pkin(4) * t179 + pkin(10) * t164) + (-t131 * mrSges(5,2) - t169) * t154 + t173; (t106 * t96 + t107 * t95) * mrSges(7,3) + (-t68 - t180) * t154 + (-mrSges(5,2) + t168) * t150 + m(6) * (t172 * t189 + t190) + m(7) * (-t133 * t154 - t78 * t95 + t79 * t96); -0.2e1 * pkin(4) * t112 + t106 * t70 + t107 * t71 + t153 * t115 + t149 * t117 + 0.2e1 * t133 * t68 + Ifges(5,3) + m(7) * (t133 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(6) * (pkin(10) ^ 2 * t172 + pkin(4) ^ 2) + 0.2e1 * (t106 * t79 - t107 * t78) * mrSges(7,3) + 0.2e1 * pkin(10) * t168; t6 * mrSges(6,1) - t7 * mrSges(6,2) + (m(7) * (t148 * t3 + t152 * t2) + t148 * t13 + t152 * t14) * pkin(5) + t161 + t15; t64 * mrSges(6,1) - t65 * mrSges(6,2) + t188 * t154 + (m(7) * (t148 * t35 + t152 * t34) + t148 * t80 + t152 * t81) * pkin(5) + t163 + t167; m(7) * (t148 * t96 - t152 * t95) * pkin(5) + t206; -t165 * pkin(10) + (m(7) * (t148 * t79 + t152 * t78) + (t106 * t148 - t107 * t152) * mrSges(7,3)) * pkin(5) + t162 + t114; m(7) * (t148 ^ 2 + t152 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t159 - t188; t161; t163 - t183; -t57; t162; Ifges(7,3) + t159; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
