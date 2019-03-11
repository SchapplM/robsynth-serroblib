% Calculate joint inertia matrix for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:26
% EndTime: 2019-03-09 20:12:31
% DurationCPUTime: 1.99s
% Computational Cost: add. (2985->416), mult. (6432->579), div. (0->0), fcn. (6810->10), ass. (0->157)
t212 = pkin(4) + pkin(9);
t152 = sin(qJ(3));
t156 = cos(qJ(3));
t211 = t152 ^ 2 + t156 ^ 2;
t151 = sin(qJ(5));
t200 = -t151 / 0.2e1;
t210 = -m(5) * pkin(3) + mrSges(5,2);
t150 = sin(qJ(6));
t154 = cos(qJ(6));
t155 = cos(qJ(5));
t148 = sin(pkin(6));
t157 = cos(qJ(2));
t189 = t148 * t157;
t149 = cos(pkin(6));
t153 = sin(qJ(2));
t190 = t148 * t153;
t91 = -t149 * t156 + t152 * t190;
t60 = t151 * t189 + t155 * t91;
t61 = t151 * t91 - t155 * t189;
t29 = -t150 * t61 + t154 * t60;
t209 = t29 / 0.2e1;
t30 = t150 * t60 + t154 * t61;
t208 = t30 / 0.2e1;
t207 = t60 / 0.2e1;
t192 = Ifges(6,4) * t155;
t80 = Ifges(6,5) * t152 + (-Ifges(6,1) * t151 - t192) * t156;
t206 = t80 / 0.2e1;
t166 = t150 * t151 - t154 * t155;
t81 = t166 * t156;
t205 = t81 / 0.2e1;
t101 = -t150 * t155 - t154 * t151;
t82 = t101 * t156;
t204 = t82 / 0.2e1;
t203 = pkin(3) + pkin(10);
t202 = t101 / 0.2e1;
t201 = -t166 / 0.2e1;
t199 = -t155 / 0.2e1;
t198 = pkin(1) * t157;
t125 = pkin(8) * t190;
t93 = t149 * t198 - t125;
t197 = t93 * mrSges(3,1);
t94 = t149 * t153 * pkin(1) + pkin(8) * t189;
t196 = t94 * mrSges(3,2);
t195 = Ifges(5,1) + Ifges(4,3);
t194 = -pkin(11) - t203;
t74 = pkin(9) * t149 + t94;
t75 = (-pkin(2) * t157 - pkin(9) * t153 - pkin(1)) * t148;
t37 = -t152 * t74 + t156 * t75;
t34 = pkin(3) * t189 - t37;
t92 = t149 * t152 + t156 * t190;
t23 = pkin(4) * t92 + pkin(10) * t189 + t34;
t73 = t125 + (-pkin(2) - t198) * t149;
t162 = -qJ(4) * t92 + t73;
t25 = t203 * t91 + t162;
t7 = t151 * t23 + t155 * t25;
t38 = t152 * t75 + t156 * t74;
t57 = -Ifges(7,5) * t166 + Ifges(7,6) * t101;
t193 = Ifges(6,4) * t151;
t191 = t155 * mrSges(6,1);
t121 = t212 * t152;
t178 = -qJ(4) * t152 - pkin(2);
t96 = -t156 * t203 + t178;
t52 = t151 * t121 + t155 * t96;
t188 = t155 * t156;
t187 = Ifges(4,5) * t152 + Ifges(4,6) * t156;
t186 = t211 * pkin(9) ^ 2;
t122 = t212 * t156;
t185 = t151 ^ 2 + t155 ^ 2;
t8 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t92;
t20 = Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * t92;
t184 = t101 ^ 2 + t166 ^ 2;
t39 = Ifges(7,5) * t82 + Ifges(7,6) * t81 + Ifges(7,3) * t152;
t137 = Ifges(6,5) * t155;
t183 = Ifges(6,6) * t200 + t137 / 0.2e1 + t57 / 0.2e1;
t182 = Ifges(3,5) * t190 + Ifges(3,6) * t189 + Ifges(3,3) * t149;
t181 = t189 / 0.2e1;
t180 = m(6) * t185;
t179 = -mrSges(7,1) * t166 + t101 * mrSges(7,2);
t6 = -t151 * t25 + t155 * t23;
t177 = t185 * mrSges(6,3);
t176 = (Ifges(5,4) - Ifges(4,5)) * t92 + (-Ifges(5,5) + Ifges(4,6)) * t91;
t4 = pkin(5) * t92 - pkin(11) * t61 + t6;
t5 = pkin(11) * t60 + t7;
t2 = -t150 * t5 + t154 * t4;
t3 = t150 * t4 + t154 * t5;
t174 = t101 * t3 + t166 * t2;
t173 = t151 * t7 + t155 * t6;
t172 = -t151 * mrSges(6,2) + t191;
t105 = t155 * t121;
t46 = pkin(5) * t152 + t105 + (pkin(11) * t156 - t96) * t151;
t48 = -pkin(11) * t188 + t52;
t14 = -t150 * t48 + t154 * t46;
t15 = t150 * t46 + t154 * t48;
t171 = t101 * t15 + t14 * t166;
t108 = t194 * t151;
t109 = t194 * t155;
t62 = -t108 * t150 + t109 * t154;
t63 = t108 * t154 + t109 * t150;
t170 = t101 * t63 + t166 * t62;
t51 = -t151 * t96 + t105;
t169 = t151 * t52 + t155 * t51;
t168 = t62 * mrSges(7,1) - t63 * mrSges(7,2) + t57;
t167 = t101 * t150 + t154 * t166;
t33 = qJ(4) * t189 - t38;
t165 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t65 = t92 * mrSges(5,1) - mrSges(5,2) * t189;
t164 = t14 * mrSges(7,1) - t15 * mrSges(7,2) + t39;
t163 = (mrSges(7,1) * t154 - mrSges(7,2) * t150) * pkin(5);
t26 = -pkin(4) * t91 - t33;
t78 = Ifges(6,3) * t152 + (-Ifges(6,5) * t151 - Ifges(6,6) * t155) * t156;
t159 = qJ(4) ^ 2;
t130 = pkin(5) * t151 + qJ(4);
t120 = Ifges(4,1) * t152 + Ifges(4,4) * t156;
t119 = Ifges(6,1) * t155 - t193;
t118 = Ifges(4,4) * t152 + Ifges(4,2) * t156;
t117 = -Ifges(6,2) * t151 + t192;
t115 = -Ifges(5,2) * t152 - Ifges(5,6) * t156;
t114 = -Ifges(5,6) * t152 - Ifges(5,3) * t156;
t113 = mrSges(6,1) * t151 + mrSges(6,2) * t155;
t112 = -mrSges(4,1) * t156 + mrSges(4,2) * t152;
t111 = mrSges(5,2) * t156 - mrSges(5,3) * t152;
t110 = -pkin(3) * t156 + t178;
t107 = -mrSges(6,2) * t152 - mrSges(6,3) * t188;
t106 = mrSges(6,3) * t151 * t156 + mrSges(6,1) * t152;
t95 = t172 * t156;
t90 = pkin(5) * t188 + t122;
t79 = Ifges(6,6) * t152 + (-Ifges(6,2) * t155 - t193) * t156;
t69 = mrSges(7,1) * t152 - mrSges(7,3) * t82;
t68 = -mrSges(7,2) * t152 + mrSges(7,3) * t81;
t67 = -mrSges(4,1) * t189 - mrSges(4,3) * t92;
t66 = mrSges(4,2) * t189 - mrSges(4,3) * t91;
t64 = mrSges(5,1) * t91 + mrSges(5,3) * t189;
t59 = -Ifges(7,1) * t166 + Ifges(7,4) * t101;
t58 = -Ifges(7,4) * t166 + Ifges(7,2) * t101;
t56 = -mrSges(7,1) * t101 - mrSges(7,2) * t166;
t50 = -mrSges(5,2) * t91 - mrSges(5,3) * t92;
t49 = mrSges(4,1) * t91 + mrSges(4,2) * t92;
t47 = -mrSges(7,1) * t81 + mrSges(7,2) * t82;
t45 = Ifges(4,1) * t92 - Ifges(4,4) * t91 - Ifges(4,5) * t189;
t44 = Ifges(4,4) * t92 - Ifges(4,2) * t91 - Ifges(4,6) * t189;
t43 = -Ifges(5,4) * t189 - Ifges(5,2) * t92 + Ifges(5,6) * t91;
t42 = -Ifges(5,5) * t189 - Ifges(5,6) * t92 + Ifges(5,3) * t91;
t41 = Ifges(7,1) * t82 + Ifges(7,4) * t81 + Ifges(7,5) * t152;
t40 = Ifges(7,4) * t82 + Ifges(7,2) * t81 + Ifges(7,6) * t152;
t36 = mrSges(6,1) * t92 - mrSges(6,3) * t61;
t35 = -mrSges(6,2) * t92 + mrSges(6,3) * t60;
t32 = pkin(3) * t91 + t162;
t31 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t22 = Ifges(6,1) * t61 + Ifges(6,4) * t60 + Ifges(6,5) * t92;
t21 = Ifges(6,4) * t61 + Ifges(6,2) * t60 + Ifges(6,6) * t92;
t17 = mrSges(7,1) * t92 - mrSges(7,3) * t30;
t16 = -mrSges(7,2) * t92 + mrSges(7,3) * t29;
t12 = -pkin(5) * t60 + t26;
t11 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t10 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t92;
t9 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t92;
t1 = [0.2e1 * t37 * t67 + 0.2e1 * t73 * t49 + t60 * t21 + t61 * t22 + 0.2e1 * t33 * t64 + 0.2e1 * t34 * t65 + 0.2e1 * t38 * t66 + m(3) * (pkin(1) ^ 2 * t148 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2 + t73 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(7) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t26 ^ 2 + t6 ^ 2 + t7 ^ 2) + 0.2e1 * t32 * t50 + 0.2e1 * t7 * t35 + 0.2e1 * t6 * t36 + t29 * t9 + t30 * t10 + 0.2e1 * t26 * t31 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + 0.2e1 * t12 * t11 + (t45 + t20 + t8 - t43) * t92 + (-t44 + t42) * t91 + Ifges(2,3) + ((-0.2e1 * t93 * mrSges(3,3) + Ifges(3,5) * t149 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t153) * t148) * t153 + (0.2e1 * t94 * mrSges(3,3) + Ifges(3,6) * t149 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t153 + (Ifges(3,2) + t195) * t157) * t148 + t176) * t157) * t148 + (t182 - 0.2e1 * t196 + 0.2e1 * t197) * t149; -t196 + m(6) * (t122 * t26 + t51 * t6 + t52 * t7) + m(7) * (t12 * t90 + t14 * t2 + t15 * t3) + (-t115 / 0.2e1 + t120 / 0.2e1 + t78 / 0.2e1 + t39 / 0.2e1) * t92 + (t114 / 0.2e1 - t118 / 0.2e1) * t91 + t6 * t106 + t7 * t107 + t110 * t50 + t32 * t111 + t73 * t112 + t122 * t31 + t26 * t95 + t90 * t11 + t3 * t68 + t2 * t69 + t182 + t12 * t47 - pkin(2) * t49 + t51 * t36 + t52 * t35 + t15 * t16 + t14 * t17 + t10 * t204 + t9 * t205 + t61 * t206 + t79 * t207 + t41 * t208 + t40 * t209 + t197 + (-t43 / 0.2e1 + t45 / 0.2e1 + t20 / 0.2e1 + t8 / 0.2e1 + Ifges(5,4) * t181 + t34 * mrSges(5,1) - t37 * mrSges(4,3) + (t65 - t67) * pkin(9)) * t152 + m(4) * (-pkin(2) * t73 + (-t152 * t37 + t156 * t38) * pkin(9)) - t187 * t189 / 0.2e1 + m(5) * (t110 * t32 + (t152 * t34 - t156 * t33) * pkin(9)) + (-t42 / 0.2e1 + t44 / 0.2e1 + Ifges(5,5) * t181 + t21 * t199 + t22 * t200 - t33 * mrSges(5,1) + t38 * mrSges(4,3) + (-t64 + t66) * pkin(9)) * t156; -0.2e1 * pkin(2) * t112 + 0.2e1 * t51 * t106 + 0.2e1 * t52 * t107 + 0.2e1 * t110 * t111 + 0.2e1 * t122 * t95 + 0.2e1 * t14 * t69 + 0.2e1 * t15 * t68 + t81 * t40 + t82 * t41 + 0.2e1 * t90 * t47 + Ifges(3,3) + (-t151 * t80 - t155 * t79 - t114 + t118) * t156 + (-t115 + t120 + t78 + t39) * t152 + m(5) * (t110 ^ 2 + t186) + m(4) * (pkin(2) ^ 2 + t186) + m(6) * (t122 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(7) * (t14 ^ 2 + t15 ^ 2 + t90 ^ 2) + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * pkin(9) * t211; t130 * t11 + t26 * t113 + t61 * t119 / 0.2e1 + t62 * t17 + t63 * t16 - pkin(3) * t65 - t176 + t12 * t56 + t37 * mrSges(4,1) - t38 * mrSges(4,2) - t33 * mrSges(5,3) + t34 * mrSges(5,2) + (t22 / 0.2e1 - t203 * t36 - t6 * mrSges(6,3)) * t155 + (-t21 / 0.2e1 - t203 * t35 - t7 * mrSges(6,3)) * t151 + m(6) * (qJ(4) * t26 - t173 * t203) + (-t64 + t31) * qJ(4) + m(7) * (t12 * t130 + t2 * t62 + t3 * t63) + m(5) * (-pkin(3) * t34 - qJ(4) * t33) + t10 * t201 + t9 * t202 + t117 * t207 + t59 * t208 + t58 * t209 + t174 * mrSges(7,3) + t183 * t92 - t195 * t189; t130 * t47 + t122 * t113 + t40 * t202 + t41 * t201 + qJ(4) * t95 + t59 * t204 + t90 * t56 + t63 * t68 + t62 * t69 + t58 * t205 + t171 * mrSges(7,3) + (-t51 * mrSges(6,3) - t106 * t203 + t206) * t155 + (-t79 / 0.2e1 - t203 * t107 - t52 * mrSges(6,3)) * t151 + m(6) * (qJ(4) * t122 - t169 * t203) + m(7) * (t130 * t90 + t14 * t62 + t15 * t63) + (qJ(4) * mrSges(5,1) + t117 * t199 + t119 * t200 - Ifges(5,5)) * t156 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t183) * t152 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t156 + (-mrSges(4,1) + t210) * t152) * pkin(9) + t187; -0.2e1 * pkin(3) * mrSges(5,2) + t101 * t58 - t166 * t59 - t151 * t117 + t155 * t119 + 0.2e1 * t130 * t56 + m(7) * (t130 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t185 * t203 ^ 2 + t159) + m(5) * (pkin(3) ^ 2 + t159) + t195 + 0.2e1 * (t113 + mrSges(5,3)) * qJ(4) + 0.2e1 * t170 * mrSges(7,3) + 0.2e1 * t203 * t177; m(5) * t34 + m(6) * t173 - m(7) * t174 - t101 * t16 + t151 * t35 + t155 * t36 - t166 * t17 + t65; -t101 * t68 - t166 * t69 + t155 * t106 + t151 * t107 + (m(5) * pkin(9) + mrSges(5,1)) * t152 - m(7) * t171 + m(6) * t169; -m(7) * t170 - mrSges(7,3) * t184 - t180 * t203 - t177 + t210; m(7) * t184 + m(5) + t180; t6 * mrSges(6,1) - t7 * mrSges(6,2) + (m(7) * (t150 * t3 + t154 * t2) + t154 * t17 + t150 * t16) * pkin(5) + t165 + t20; t51 * mrSges(6,1) - t52 * mrSges(6,2) + (t150 * t68 + t154 * t69 + m(7) * (t14 * t154 + t15 * t150)) * pkin(5) + t78 + t164; -t203 * t191 + t137 + (mrSges(6,2) * t203 - Ifges(6,6)) * t151 + (m(7) * (t150 * t63 + t154 * t62) + t167 * mrSges(7,3)) * pkin(5) + t168; -m(7) * pkin(5) * t167 + t172 + t179; Ifges(6,3) + Ifges(7,3) + m(7) * (t150 ^ 2 + t154 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t163; t165; t164; t168; t179; Ifges(7,3) + t163; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
