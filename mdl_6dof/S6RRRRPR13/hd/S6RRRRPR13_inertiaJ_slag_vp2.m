% Calculate joint inertia matrix for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:10
% EndTime: 2019-03-09 23:51:15
% DurationCPUTime: 2.40s
% Computational Cost: add. (3115->466), mult. (6902->642), div. (0->0), fcn. (7284->10), ass. (0->171)
t162 = sin(qJ(4));
t166 = cos(qJ(4));
t223 = t162 ^ 2 + t166 ^ 2;
t222 = 2 * pkin(9);
t160 = sin(pkin(6));
t168 = cos(qJ(2));
t199 = t160 * t168;
t163 = sin(qJ(3));
t167 = cos(qJ(3));
t164 = sin(qJ(2));
t200 = t160 * t164;
t202 = cos(pkin(6));
t99 = t202 * t163 + t167 * t200;
t61 = t162 * t99 + t166 * t199;
t62 = -t162 * t199 + t166 * t99;
t98 = t163 * t200 - t202 * t167;
t18 = Ifges(5,5) * t62 - Ifges(5,6) * t61 + Ifges(5,3) * t98;
t19 = Ifges(6,4) * t62 + Ifges(6,2) * t98 + Ifges(6,6) * t61;
t221 = t18 + t19;
t220 = t160 ^ 2;
t161 = sin(qJ(6));
t165 = cos(qJ(6));
t27 = -t161 * t62 + t165 * t61;
t219 = t27 / 0.2e1;
t28 = t161 * t61 + t165 * t62;
t218 = t28 / 0.2e1;
t109 = -t161 * t166 + t162 * t165;
t90 = t109 * t163;
t217 = t90 / 0.2e1;
t177 = t161 * t162 + t165 * t166;
t91 = t177 * t163;
t216 = t91 / 0.2e1;
t215 = pkin(4) + pkin(5);
t214 = pkin(10) - pkin(11);
t213 = -t177 / 0.2e1;
t212 = t109 / 0.2e1;
t211 = pkin(9) * t163;
t210 = pkin(9) * t167;
t209 = -Ifges(5,3) - Ifges(6,2);
t188 = pkin(1) * t202;
t100 = -pkin(8) * t200 + t168 * t188;
t78 = -t202 * pkin(2) - t100;
t34 = t98 * pkin(3) - t99 * pkin(10) + t78;
t101 = pkin(8) * t199 + t164 * t188;
t79 = t202 * pkin(9) + t101;
t80 = (-pkin(2) * t168 - pkin(9) * t164 - pkin(1)) * t160;
t42 = t163 * t80 + t167 * t79;
t36 = -pkin(10) * t199 + t42;
t13 = t162 * t34 + t166 * t36;
t41 = -t163 * t79 + t167 * t80;
t208 = Ifges(5,4) * t162;
t207 = Ifges(5,4) * t166;
t206 = Ifges(6,5) * t162;
t205 = Ifges(6,5) * t166;
t114 = -qJ(5) * t161 - t165 * t215;
t204 = t114 * mrSges(7,1);
t115 = qJ(5) * t165 - t161 * t215;
t203 = t115 * mrSges(7,2);
t201 = qJ(5) * t166;
t198 = t162 * t163;
t197 = t163 * t166;
t58 = Ifges(7,5) * t109 - Ifges(7,6) * t177;
t112 = t167 * mrSges(6,1) + mrSges(6,2) * t197;
t196 = Ifges(6,4) * t197 + Ifges(6,6) * t198;
t117 = -pkin(3) * t167 - pkin(10) * t163 - pkin(2);
t77 = t162 * t117 + t166 * t210;
t122 = Ifges(5,5) * t162 + Ifges(5,6) * t166;
t195 = Ifges(4,5) * t163 + Ifges(4,6) * t167;
t194 = t223 * pkin(10) ^ 2;
t5 = Ifges(7,5) * t28 + Ifges(7,6) * t27 - Ifges(7,3) * t98;
t10 = t98 * qJ(5) + t13;
t35 = pkin(3) * t199 - t41;
t43 = Ifges(7,5) * t91 + Ifges(7,6) * t90 + Ifges(7,3) * t167;
t17 = Ifges(6,5) * t62 + Ifges(6,6) * t98 + Ifges(6,3) * t61;
t20 = Ifges(5,4) * t62 - Ifges(5,2) * t61 + Ifges(5,6) * t98;
t193 = -t20 / 0.2e1 + t17 / 0.2e1;
t21 = Ifges(6,1) * t62 + Ifges(6,4) * t98 + Ifges(6,5) * t61;
t22 = Ifges(5,1) * t62 - Ifges(5,4) * t61 + Ifges(5,5) * t98;
t192 = t21 / 0.2e1 + t22 / 0.2e1;
t83 = -Ifges(6,6) * t167 + (Ifges(6,3) * t162 + t205) * t163;
t86 = -Ifges(5,6) * t167 + (-Ifges(5,2) * t162 + t207) * t163;
t191 = t83 / 0.2e1 - t86 / 0.2e1;
t87 = -Ifges(6,4) * t167 + (Ifges(6,1) * t166 + t206) * t163;
t88 = -Ifges(5,5) * t167 + (Ifges(5,1) * t166 - t208) * t163;
t190 = t87 / 0.2e1 + t88 / 0.2e1;
t189 = Ifges(3,5) * t200 + Ifges(3,6) * t199 + Ifges(3,3) * t202;
t121 = -Ifges(6,3) * t166 + t206;
t124 = Ifges(5,2) * t166 + t208;
t187 = t121 / 0.2e1 - t124 / 0.2e1;
t126 = Ifges(6,1) * t162 - t205;
t127 = Ifges(5,1) * t162 + t207;
t186 = t126 / 0.2e1 + t127 / 0.2e1;
t40 = -t98 * mrSges(6,1) + t62 * mrSges(6,2);
t185 = qJ(5) * t162 + pkin(3);
t12 = -t162 * t36 + t166 * t34;
t123 = Ifges(6,4) * t162 - Ifges(6,6) * t166;
t141 = t162 * t210;
t76 = t117 * t166 - t141;
t183 = Ifges(5,5) * t197 - Ifges(5,6) * t198;
t182 = t122 / 0.2e1 + t123 / 0.2e1 - t58 / 0.2e1;
t68 = -qJ(5) * t167 + t77;
t181 = t162 * mrSges(5,1) + t166 * mrSges(5,2);
t180 = t162 * mrSges(6,1) - t166 * mrSges(6,3);
t179 = t165 * mrSges(7,1) - t161 * mrSges(7,2);
t178 = -pkin(4) * t162 + t201;
t176 = Ifges(4,5) * t99 - Ifges(4,6) * t98 - Ifges(4,3) * t199;
t175 = qJ(5) * t62 - t35;
t129 = t214 * t162;
t130 = t214 * t166;
t66 = t129 * t165 - t130 * t161;
t67 = t129 * t161 + t130 * t165;
t174 = t66 * mrSges(7,1) - t67 * mrSges(7,2) + t58;
t3 = -pkin(11) * t62 - t215 * t98 - t12;
t4 = pkin(11) * t61 + t10;
t1 = -t161 * t4 + t165 * t3;
t2 = t161 * t3 + t165 * t4;
t173 = t1 * mrSges(7,1) - t2 * mrSges(7,2) + t5;
t155 = t167 * pkin(4);
t49 = pkin(5) * t167 + t141 + t155 + (-pkin(11) * t163 - t117) * t166;
t51 = pkin(11) * t198 + t68;
t23 = -t161 * t51 + t165 * t49;
t24 = t161 * t49 + t165 * t51;
t172 = t23 * mrSges(7,1) - t24 * mrSges(7,2) + t43;
t171 = pkin(9) ^ 2;
t159 = t167 ^ 2;
t157 = t163 ^ 2;
t153 = t157 * t171;
t128 = Ifges(4,1) * t163 + Ifges(4,4) * t167;
t125 = Ifges(4,4) * t163 + Ifges(4,2) * t167;
t120 = -mrSges(4,1) * t167 + mrSges(4,2) * t163;
t119 = -mrSges(5,1) * t166 + mrSges(5,2) * t162;
t118 = -mrSges(6,1) * t166 - mrSges(6,3) * t162;
t116 = -pkin(4) * t166 - t185;
t113 = -mrSges(6,2) * t198 - mrSges(6,3) * t167;
t111 = -mrSges(5,1) * t167 - mrSges(5,3) * t197;
t110 = mrSges(5,2) * t167 - mrSges(5,3) * t198;
t104 = t215 * t166 + t185;
t103 = t181 * t163;
t102 = t180 * t163;
t89 = (pkin(9) - t178) * t163;
t85 = -Ifges(6,2) * t167 + t196;
t84 = -Ifges(5,3) * t167 + t183;
t71 = mrSges(7,1) * t167 - mrSges(7,3) * t91;
t70 = -mrSges(7,2) * t167 + mrSges(7,3) * t90;
t69 = t155 - t76;
t65 = -mrSges(4,1) * t199 - mrSges(4,3) * t99;
t64 = mrSges(4,2) * t199 - mrSges(4,3) * t98;
t63 = (-t215 * t162 - pkin(9) + t201) * t163;
t60 = Ifges(7,1) * t109 - Ifges(7,4) * t177;
t59 = Ifges(7,4) * t109 - Ifges(7,2) * t177;
t57 = mrSges(7,1) * t177 + mrSges(7,2) * t109;
t50 = mrSges(4,1) * t98 + mrSges(4,2) * t99;
t48 = -mrSges(7,1) * t90 + mrSges(7,2) * t91;
t47 = Ifges(4,1) * t99 - Ifges(4,4) * t98 - Ifges(4,5) * t199;
t46 = Ifges(4,4) * t99 - Ifges(4,2) * t98 - Ifges(4,6) * t199;
t45 = Ifges(7,1) * t91 + Ifges(7,4) * t90 + Ifges(7,5) * t167;
t44 = Ifges(7,4) * t91 + Ifges(7,2) * t90 + Ifges(7,6) * t167;
t39 = mrSges(5,1) * t98 - mrSges(5,3) * t62;
t38 = -mrSges(5,2) * t98 - mrSges(5,3) * t61;
t37 = -mrSges(6,2) * t61 + mrSges(6,3) * t98;
t30 = mrSges(5,1) * t61 + mrSges(5,2) * t62;
t29 = mrSges(6,1) * t61 - mrSges(6,3) * t62;
t16 = -mrSges(7,1) * t98 - mrSges(7,3) * t28;
t15 = mrSges(7,2) * t98 + mrSges(7,3) * t27;
t14 = pkin(4) * t61 - t175;
t11 = -pkin(4) * t98 - t12;
t9 = -mrSges(7,1) * t27 + mrSges(7,2) * t28;
t8 = -t215 * t61 + t175;
t7 = Ifges(7,1) * t28 + Ifges(7,4) * t27 - Ifges(7,5) * t98;
t6 = Ifges(7,4) * t28 + Ifges(7,2) * t27 - Ifges(7,6) * t98;
t25 = [t99 * t47 + 0.2e1 * t78 * t50 + 0.2e1 * t42 * t64 + 0.2e1 * t41 * t65 + 0.2e1 * t35 * t30 + 0.2e1 * t10 * t37 + 0.2e1 * t13 * t38 + 0.2e1 * t12 * t39 + 0.2e1 * t11 * t40 + t27 * t6 + t28 * t7 + 0.2e1 * t14 * t29 + 0.2e1 * t2 * t15 + 0.2e1 * t1 * t16 + 0.2e1 * t8 * t9 - 0.2e1 * t220 * pkin(1) * (-mrSges(3,1) * t168 + mrSges(3,2) * t164) + m(4) * (t41 ^ 2 + t42 ^ 2 + t78 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t35 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t14 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + (t21 + t22) * t62 + (t17 - t20) * t61 + m(3) * (pkin(1) ^ 2 * t220 + t100 ^ 2 + t101 ^ 2) + (Ifges(3,6) * t202 + (Ifges(3,4) * t164 + Ifges(3,2) * t168) * t160) * t199 + (Ifges(3,5) * t202 + (Ifges(3,1) * t164 + Ifges(3,4) * t168) * t160) * t200 + t202 * t189 + 0.2e1 * t100 * (t202 * mrSges(3,1) - mrSges(3,3) * t200) + 0.2e1 * t101 * (-t202 * mrSges(3,2) + mrSges(3,3) * t199) - t176 * t199 + Ifges(2,3) + (-t46 - t5 + t221) * t98; t78 * t120 + t99 * t128 / 0.2e1 + t13 * t110 + t12 * t111 + t11 * t112 + t10 * t113 + t100 * mrSges(3,1) - t101 * mrSges(3,2) + t14 * t102 + t35 * t103 + t76 * t39 + t77 * t38 + t89 * t29 + t7 * t216 + t6 * t217 + t45 * t218 + t44 * t219 + t63 * t9 + t68 * t37 + t69 * t40 + t2 * t70 + t1 * t71 + t8 * t48 - pkin(2) * t50 + t23 * t16 + t24 * t15 + m(4) * (-pkin(2) * t78 + (-t41 * t163 + t42 * t167) * pkin(9)) + t190 * t62 + t191 * t61 + (t47 / 0.2e1 - t41 * mrSges(4,3) + t192 * t166 + t193 * t162 + (-t65 + t30) * pkin(9)) * t163 + t189 + m(5) * (t12 * t76 + t13 * t77 + t35 * t211) - t195 * t199 / 0.2e1 + (-t125 / 0.2e1 - t43 / 0.2e1 + t84 / 0.2e1 + t85 / 0.2e1) * t98 + m(6) * (t10 * t68 + t11 * t69 + t14 * t89) + m(7) * (t1 * t23 + t2 * t24 + t63 * t8) + (pkin(9) * t64 + t42 * mrSges(4,3) + t46 / 0.2e1 - t18 / 0.2e1 - t19 / 0.2e1 + t5 / 0.2e1) * t167; -0.2e1 * pkin(2) * t120 + 0.2e1 * t89 * t102 + 0.2e1 * t77 * t110 + 0.2e1 * t76 * t111 + 0.2e1 * t69 * t112 + 0.2e1 * t68 * t113 + 0.2e1 * t23 * t71 + 0.2e1 * t24 * t70 + t90 * t44 + t91 * t45 + 0.2e1 * t63 * t48 + Ifges(3,3) + (t157 + t159) * mrSges(4,3) * t222 + (t125 - t84 - t85 + t43) * t167 + m(4) * (pkin(2) ^ 2 + t159 * t171 + t153) + m(5) * (t76 ^ 2 + t77 ^ 2 + t153) + m(6) * (t68 ^ 2 + t69 ^ 2 + t89 ^ 2) + m(7) * (t23 ^ 2 + t24 ^ 2 + t63 ^ 2) + (t103 * t222 + t128 + (t87 + t88) * t166 + (t83 - t86) * t162) * t163; t14 * t118 + t35 * t119 + t6 * t213 + t7 * t212 + t116 * t29 + t104 * t9 + m(7) * (t1 * t66 + t104 * t8 + t2 * t67) + t182 * t98 + t186 * t62 + t187 * t61 + (-t1 * t109 - t177 * t2) * mrSges(7,3) + t66 * t16 + t67 * t15 + t8 * t57 + t59 * t219 + t60 * t218 + t41 * mrSges(4,1) - t42 * mrSges(4,2) - pkin(3) * t30 + (t11 * mrSges(6,2) - t12 * mrSges(5,3) + (-t39 + t40) * pkin(10) + t192) * t162 + (t10 * mrSges(6,2) + t13 * mrSges(5,3) + (t37 + t38) * pkin(10) - t193) * t166 + m(6) * (t116 * t14 + (t10 * t166 + t11 * t162) * pkin(10)) + m(5) * (-pkin(3) * t35 + (-t12 * t162 + t13 * t166) * pkin(10)) + t176; t89 * t118 + t44 * t213 + t45 * t212 + t116 * t102 - pkin(3) * t103 + t104 * t48 + t60 * t216 + t59 * t217 + t63 * t57 + t67 * t70 + t66 * t71 + (-mrSges(4,1) + t119) * t211 + (-t109 * t23 - t177 * t24) * mrSges(7,3) + (t77 * mrSges(5,3) + t68 * mrSges(6,2) + t186 * t163 + (t110 + t113) * pkin(10) - t191) * t166 + (t69 * mrSges(6,2) - t76 * mrSges(5,3) + t187 * t163 + (-t111 + t112) * pkin(10) + t190) * t162 + m(5) * (-pkin(3) * t211 + (-t162 * t76 + t166 * t77) * pkin(10)) + m(6) * (t116 * t89 + (t162 * t69 + t166 * t68) * pkin(10)) + m(7) * (t104 * t63 + t23 * t66 + t24 * t67) + (-pkin(9) * mrSges(4,2) - t182) * t167 + t195; -0.2e1 * pkin(3) * t119 + 0.2e1 * t104 * t57 - t177 * t59 + t109 * t60 + 0.2e1 * t116 * t118 + Ifges(4,3) + (-t121 + t124) * t166 + (t126 + t127) * t162 + 0.2e1 * (-t109 * t66 - t177 * t67) * mrSges(7,3) + m(7) * (t104 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t116 ^ 2 + t194) + m(5) * (pkin(3) ^ 2 + t194) + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * pkin(10) * t223; m(7) * (t1 * t114 + t115 * t2) + m(6) * (-pkin(4) * t11 + qJ(5) * t10) + t114 * t16 + t115 * t15 + qJ(5) * t37 - pkin(4) * t40 + t10 * mrSges(6,3) - t11 * mrSges(6,1) + t12 * mrSges(5,1) - t13 * mrSges(5,2) - t173 + t221; -pkin(4) * t112 + qJ(5) * t113 + t114 * t71 + t115 * t70 + t76 * mrSges(5,1) - t77 * mrSges(5,2) + t68 * mrSges(6,3) - t69 * mrSges(6,1) + m(7) * (t114 * t23 + t115 * t24) + m(6) * (-pkin(4) * t69 + qJ(5) * t68) + t209 * t167 - t172 + t183 + t196; m(7) * (t114 * t66 + t115 * t67) + (-t109 * t114 - t115 * t177) * mrSges(7,3) + t178 * mrSges(6,2) + (m(6) * t178 - t180 - t181) * pkin(10) - t174 + t123 + t122; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t204 + 0.2e1 * t203 + 0.2e1 * qJ(5) * mrSges(6,3) + Ifges(7,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + m(7) * (t114 ^ 2 + t115 ^ 2) - t209; t161 * t15 + t165 * t16 + m(7) * (t1 * t165 + t161 * t2) + m(6) * t11 + t40; t161 * t70 + t165 * t71 + m(7) * (t161 * t24 + t165 * t23) + m(6) * t69 + t112; m(7) * (t161 * t67 + t165 * t66) + (m(6) * pkin(10) + mrSges(6,2)) * t162 + (-t109 * t165 - t161 * t177) * mrSges(7,3); -mrSges(6,1) - m(6) * pkin(4) + m(7) * (t114 * t165 + t115 * t161) - t179; m(6) + m(7) * (t161 ^ 2 + t165 ^ 2); t173; t172; t174; -Ifges(7,3) - t203 + t204; t179; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
