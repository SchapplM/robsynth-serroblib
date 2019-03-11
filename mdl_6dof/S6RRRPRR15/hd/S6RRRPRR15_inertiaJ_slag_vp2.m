% Calculate joint inertia matrix for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR15_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:39
% EndTime: 2019-03-09 20:24:46
% DurationCPUTime: 2.73s
% Computational Cost: add. (4851->494), mult. (12401->697), div. (0->0), fcn. (13675->12), ass. (0->189)
t244 = (pkin(3) + pkin(11));
t250 = 2 * t244;
t178 = sin(qJ(6));
t182 = cos(qJ(6));
t176 = cos(pkin(7));
t179 = sin(qJ(5));
t183 = cos(qJ(5));
t174 = sin(pkin(7));
t184 = cos(qJ(3));
t218 = t174 * t184;
t118 = t176 * t179 + t183 * t218;
t119 = t176 * t183 - t179 * t218;
t180 = sin(qJ(3));
t219 = t174 * t180;
t87 = -t119 * t178 + t182 * t219;
t63 = -mrSges(7,2) * t118 + mrSges(7,3) * t87;
t88 = t119 * t182 + t178 * t219;
t64 = mrSges(7,1) * t118 - mrSges(7,3) * t88;
t249 = -t178 * t64 + t182 * t63;
t177 = cos(pkin(6));
t175 = sin(pkin(6));
t185 = cos(qJ(2));
t216 = t175 * t185;
t117 = -t174 * t216 + t176 * t177;
t214 = t176 * t185;
t205 = t175 * t214;
t181 = sin(qJ(2));
t217 = t175 * t181;
t81 = -t177 * t218 + t180 * t217 - t184 * t205;
t58 = t117 * t183 + t179 * t81;
t82 = t177 * t219 + (t180 * t214 + t181 * t184) * t175;
t32 = -t178 * t58 + t182 * t82;
t57 = t117 * t179 - t81 * t183;
t17 = -mrSges(7,2) * t57 + mrSges(7,3) * t32;
t33 = t178 * t82 + t182 * t58;
t18 = mrSges(7,1) * t57 - mrSges(7,3) * t33;
t248 = t182 * t17 - t178 * t18;
t233 = pkin(1) * t177;
t125 = pkin(9) * t216 + t181 * t233;
t74 = (t174 * t177 + t205) * pkin(10) + t125;
t156 = t185 * t233;
t86 = pkin(2) * t177 + t156 + (-pkin(10) * t176 - pkin(9)) * t217;
t97 = (-pkin(10) * t174 * t181 - pkin(2) * t185 - pkin(1)) * t175;
t28 = -t180 * t74 + (t174 * t97 + t176 * t86) * t184;
t135 = -mrSges(7,1) * t182 + mrSges(7,2) * t178;
t247 = m(7) * pkin(5) + mrSges(6,1) - t135;
t21 = Ifges(6,1) * t58 - Ifges(6,4) * t57 + Ifges(6,5) * t82;
t246 = t21 / 0.2e1;
t67 = Ifges(6,1) * t119 - Ifges(6,4) * t118 + Ifges(6,5) * t219;
t245 = t67 / 0.2e1;
t227 = Ifges(7,4) * t182;
t108 = Ifges(7,6) * t179 + (-Ifges(7,2) * t178 + t227) * t183;
t243 = t108 / 0.2e1;
t228 = Ifges(7,4) * t178;
t109 = Ifges(7,5) * t179 + (Ifges(7,1) * t182 - t228) * t183;
t242 = t109 / 0.2e1;
t137 = Ifges(7,5) * t178 + Ifges(7,6) * t182;
t241 = t137 / 0.2e1;
t168 = Ifges(6,5) * t183;
t240 = -Ifges(6,6) * t179 / 0.2e1 + t168 / 0.2e1;
t139 = Ifges(7,2) * t182 + t228;
t239 = t139 / 0.2e1;
t141 = Ifges(7,1) * t178 + t227;
t238 = t141 / 0.2e1;
t142 = Ifges(6,1) * t183 - Ifges(6,4) * t179;
t237 = t142 / 0.2e1;
t236 = -t178 / 0.2e1;
t235 = t178 / 0.2e1;
t234 = t182 / 0.2e1;
t232 = pkin(2) * t184;
t11 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t35 = mrSges(6,1) * t82 - mrSges(6,3) * t58;
t231 = -t11 + t35;
t13 = pkin(4) * t82 - t244 * t117 - t28;
t53 = -t174 * t86 + t176 * t97;
t190 = -qJ(4) * t82 + t53;
t15 = t244 * t81 + t190;
t6 = t179 * t13 + t183 * t15;
t37 = Ifges(4,5) * t82 - Ifges(4,6) * t81 + Ifges(4,3) * t117;
t40 = Ifges(5,1) * t117 - Ifges(5,4) * t82 + Ifges(5,5) * t81;
t230 = t37 + t40;
t52 = -mrSges(7,1) * t87 + mrSges(7,2) * t88;
t94 = mrSges(6,1) * t219 - mrSges(6,3) * t119;
t229 = -t52 + t94;
t150 = pkin(10) * t219;
t203 = -pkin(3) - t232;
t80 = pkin(4) * t219 + t150 + (-pkin(11) + t203) * t176;
t201 = -qJ(4) * t180 - pkin(2);
t91 = (-t244 * t184 + t201) * t174;
t48 = t179 * t80 + t183 * t91;
t226 = Ifges(5,5) * t184;
t123 = -pkin(9) * t217 + t156;
t225 = t123 * mrSges(3,1);
t224 = t125 * mrSges(3,2);
t215 = t176 * t180;
t213 = t178 * t183;
t212 = t179 * t244;
t211 = t182 * t183;
t210 = t183 * t244;
t124 = pkin(2) * t215 + pkin(10) * t218;
t130 = mrSges(5,1) * t219 + t176 * mrSges(5,2);
t209 = t178 ^ 2 + t182 ^ 2;
t170 = t179 ^ 2;
t172 = t183 ^ 2;
t208 = t170 + t172;
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t57;
t19 = Ifges(6,5) * t58 - Ifges(6,6) * t57 + Ifges(6,3) * t82;
t29 = t184 * t74 + t86 * t215 + t97 * t219;
t20 = Ifges(6,4) * t58 - Ifges(6,2) * t57 + Ifges(6,6) * t82;
t207 = t8 / 0.2e1 - t20 / 0.2e1;
t42 = Ifges(7,5) * t88 + Ifges(7,6) * t87 + Ifges(7,3) * t118;
t66 = Ifges(6,4) * t119 - Ifges(6,2) * t118 + Ifges(6,6) * t219;
t206 = t42 / 0.2e1 - t66 / 0.2e1;
t65 = Ifges(6,5) * t119 - Ifges(6,6) * t118 + Ifges(6,3) * t219;
t98 = Ifges(4,5) * t219 + Ifges(4,6) * t218 + Ifges(4,3) * t176;
t204 = Ifges(3,5) * t217 + Ifges(3,6) * t216 + Ifges(3,3) * t177;
t104 = -t176 * qJ(4) - t124;
t107 = Ifges(7,5) * t211 - Ifges(7,6) * t213 + Ifges(7,3) * t179;
t140 = Ifges(6,4) * t183 - Ifges(6,2) * t179;
t202 = t107 / 0.2e1 - t140 / 0.2e1;
t62 = t82 * mrSges(5,1) + t117 * mrSges(5,2);
t200 = t208 * mrSges(6,3);
t25 = -t117 * qJ(4) - t29;
t90 = pkin(4) * t218 - t104;
t4 = pkin(12) * t82 + t6;
t16 = -pkin(4) * t81 - t25;
t7 = pkin(5) * t57 - pkin(12) * t58 + t16;
t1 = -t178 * t4 + t182 * t7;
t2 = t178 * t7 + t182 * t4;
t198 = -t1 * t178 + t182 * t2;
t5 = t13 * t183 - t15 * t179;
t197 = t6 * t179 + t5 * t183;
t196 = mrSges(7,1) * t178 + mrSges(7,2) * t182;
t46 = pkin(12) * t219 + t48;
t51 = pkin(5) * t118 - pkin(12) * t119 + t90;
t22 = -t178 * t46 + t182 * t51;
t23 = t178 * t51 + t182 * t46;
t194 = -t178 * t22 + t182 * t23;
t134 = pkin(5) * t179 - pkin(12) * t183 + qJ(4);
t95 = t134 * t182 + t178 * t212;
t96 = t134 * t178 - t182 * t212;
t193 = -t178 * t95 + t182 * t96;
t47 = -t179 * t91 + t183 * t80;
t192 = t48 * t179 + t47 * t183;
t131 = -mrSges(7,2) * t179 - mrSges(7,3) * t213;
t132 = mrSges(7,1) * t179 - mrSges(7,3) * t211;
t191 = t182 * t131 - t178 * t132;
t187 = qJ(4) ^ 2;
t173 = t244 ^ 2;
t164 = Ifges(5,1) * t176;
t159 = t172 * t244;
t158 = t172 * t173;
t136 = mrSges(6,1) * t179 + mrSges(6,2) * t183;
t129 = -mrSges(5,1) * t218 - mrSges(5,3) * t176;
t128 = -mrSges(4,2) * t176 + mrSges(4,3) * t218;
t127 = mrSges(4,1) * t176 - mrSges(4,3) * t219;
t126 = t196 * t183;
t122 = t176 * t232 - t150;
t121 = (-mrSges(4,1) * t184 + mrSges(4,2) * t180) * t174;
t120 = (mrSges(5,2) * t184 - mrSges(5,3) * t180) * t174;
t106 = t203 * t176 + t150;
t105 = (-pkin(3) * t184 + t201) * t174;
t103 = t164 + (-Ifges(5,4) * t180 - t226) * t174;
t102 = Ifges(5,4) * t176 + (-Ifges(5,2) * t180 - Ifges(5,6) * t184) * t174;
t101 = Ifges(5,5) * t176 + (-Ifges(5,6) * t180 - Ifges(5,3) * t184) * t174;
t100 = Ifges(4,5) * t176 + (Ifges(4,1) * t180 + Ifges(4,4) * t184) * t174;
t99 = Ifges(4,6) * t176 + (Ifges(4,4) * t180 + Ifges(4,2) * t184) * t174;
t93 = -mrSges(6,2) * t219 - mrSges(6,3) * t118;
t68 = mrSges(6,1) * t118 + mrSges(6,2) * t119;
t61 = mrSges(5,1) * t81 - mrSges(5,3) * t117;
t60 = mrSges(4,1) * t117 - mrSges(4,3) * t82;
t59 = -mrSges(4,2) * t117 - mrSges(4,3) * t81;
t50 = -mrSges(5,2) * t81 - mrSges(5,3) * t82;
t49 = mrSges(4,1) * t81 + mrSges(4,2) * t82;
t45 = -pkin(5) * t219 - t47;
t44 = Ifges(7,1) * t88 + Ifges(7,4) * t87 + Ifges(7,5) * t118;
t43 = Ifges(7,4) * t88 + Ifges(7,2) * t87 + Ifges(7,6) * t118;
t41 = Ifges(4,1) * t82 - Ifges(4,4) * t81 + Ifges(4,5) * t117;
t39 = Ifges(4,4) * t82 - Ifges(4,2) * t81 + Ifges(4,6) * t117;
t38 = Ifges(5,4) * t117 - Ifges(5,2) * t82 + Ifges(5,6) * t81;
t36 = Ifges(5,5) * t117 - Ifges(5,6) * t82 + Ifges(5,3) * t81;
t34 = -mrSges(6,2) * t82 - mrSges(6,3) * t57;
t27 = mrSges(6,1) * t57 + mrSges(6,2) * t58;
t26 = -pkin(3) * t117 - t28;
t24 = pkin(3) * t81 + t190;
t10 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t57;
t9 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t57;
t3 = -pkin(5) * t82 - t5;
t12 = [m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t16 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2 + t53 ^ 2) + (t19 + t41 - t38) * t82 + (t36 - t39) * t81 + (t8 - t20) * t57 + m(3) * (t123 ^ 2 + t125 ^ 2) + ((t181 * Ifges(3,5) + t185 * Ifges(3,6)) * t177 + 0.2e1 * (-t123 * t181 + t125 * t185) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t185 + mrSges(3,2) * t181) + t181 * (Ifges(3,1) * t181 + Ifges(3,4) * t185) + t185 * (Ifges(3,4) * t181 + Ifges(3,2) * t185)) * t175) * t175 + Ifges(2,3) + (t204 - 0.2e1 * t224 + 0.2e1 * t225) * t177 + t230 * t117 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t16 * t27 + t32 * t9 + t33 * t10 + 0.2e1 * t6 * t34 + 0.2e1 * t5 * t35 + 0.2e1 * t24 * t50 + 0.2e1 * t53 * t49 + t58 * t21 + 0.2e1 * t29 * t59 + 0.2e1 * t28 * t60 + 0.2e1 * t25 * t61 + 0.2e1 * t26 * t62; (-pkin(2) * t49 + (-t36 / 0.2e1 + t39 / 0.2e1) * t184 + (t19 / 0.2e1 - t38 / 0.2e1 + t41 / 0.2e1) * t180) * t174 + (t37 / 0.2e1 + t40 / 0.2e1) * t176 + (t98 / 0.2e1 + t103 / 0.2e1) * t117 + t225 - t224 + m(7) * (t1 * t22 + t2 * t23 + t3 * t45) + m(6) * (t16 * t90 + t47 * t5 + t48 * t6) + m(5) * (t104 * t25 + t105 * t24 + t106 * t26) + m(4) * (-pkin(2) * t174 * t53 + t122 * t28 + t124 * t29) + (t65 / 0.2e1 + t100 / 0.2e1 - t102 / 0.2e1) * t82 + (-t99 / 0.2e1 + t101 / 0.2e1) * t81 + t204 + t206 * t57 + t207 * t118 + t58 * t245 + t119 * t246 + t22 * t18 + t23 * t17 + t32 * t43 / 0.2e1 + t33 * t44 / 0.2e1 + t45 * t11 + t47 * t35 + t48 * t34 + t3 * t52 + t2 * t63 + t1 * t64 + t16 * t68 + t87 * t9 / 0.2e1 + t88 * t10 / 0.2e1 + t90 * t27 + t6 * t93 + t5 * t94 + t104 * t61 + t105 * t50 + t106 * t62 + t24 * t120 + t53 * t121 + t122 * t60 + t124 * t59 + t28 * t127 + t29 * t128 + t25 * t129 + t26 * t130; 0.2e1 * t104 * t129 + 0.2e1 * t105 * t120 + 0.2e1 * t106 * t130 + t119 * t67 + 0.2e1 * t122 * t127 + 0.2e1 * t124 * t128 + 0.2e1 * t22 * t64 + 0.2e1 * t23 * t63 + t87 * t43 + t88 * t44 + 0.2e1 * t45 * t52 + 0.2e1 * t47 * t94 + 0.2e1 * t48 * t93 + 0.2e1 * t90 * t68 + Ifges(3,3) + (t98 + t103) * t176 + (t42 - t66) * t118 + (-0.2e1 * pkin(2) * t121 + (-t101 + t99) * t184 + (t100 - t102 + t65) * t180) * t174 + m(7) * (t22 ^ 2 + t23 ^ 2 + t45 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2 + t90 ^ 2) + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) + m(4) * (pkin(2) ^ 2 * t174 ^ 2 + t122 ^ 2 + t124 ^ 2); m(6) * (qJ(4) * t16 - t197 * t244) + (-t5 * mrSges(6,3) + t10 * t234 - t231 * t244 + t9 * t236 + t246) * t183 + t230 + (t27 - t61) * qJ(4) + m(5) * (-pkin(3) * t26 - qJ(4) * t25) + t202 * t57 + (-t6 * mrSges(6,3) - t244 * t34 + t207) * t179 + m(7) * (t1 * t95 + t2 * t96 + t3 * t210) - t25 * mrSges(5,3) + t26 * mrSges(5,2) + t28 * mrSges(4,1) - t29 * mrSges(4,2) - pkin(3) * t62 + t95 * t18 + t96 * t17 + t32 * t243 + t33 * t242 + t3 * t126 + t2 * t131 + t1 * t132 + t16 * t136 + t82 * t240 + t58 * t237; (-t47 * mrSges(6,3) - t229 * t244 + t44 * t234 + t43 * t236 + t245) * t183 + m(6) * (qJ(4) * t90 - t192 * t244) + (-t226 + (-Ifges(5,4) + t240) * t180) * t174 + (-t48 * mrSges(6,3) - t244 * t93 + t206) * t179 + (t68 - t129) * qJ(4) + t164 + m(5) * (-pkin(3) * t106 - qJ(4) * t104) + m(7) * (t45 * t210 + t22 * t95 + t23 * t96) + t98 + t202 * t118 + t95 * t64 + t96 * t63 - t104 * mrSges(5,3) + t106 * mrSges(5,2) + t87 * t243 + t88 * t242 + t122 * mrSges(4,1) - t124 * mrSges(4,2) + t45 * t126 - pkin(3) * t130 + t23 * t131 + t22 * t132 + t90 * t136 + t119 * t237; -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t96 * t131 + 0.2e1 * t95 * t132 + Ifges(5,1) + Ifges(4,3) + (t107 - t140) * t179 + m(7) * (t95 ^ 2 + t96 ^ 2 + t158) + m(6) * (t170 * t173 + t158 + t187) + m(5) * (pkin(3) ^ 2 + t187) + (-t108 * t178 + t109 * t182 + t126 * t250 + t142) * t183 + 0.2e1 * (t136 + mrSges(5,3)) * qJ(4) + t200 * t250; t231 * t183 + (t34 + t248) * t179 + m(7) * (t198 * t179 - t183 * t3) + m(6) * t197 + m(5) * t26 + t62; t229 * t183 + (t93 + t249) * t179 + m(7) * (t194 * t179 - t183 * t45) + m(6) * t192 + m(5) * t106 + t130; -m(5) * pkin(3) - t183 * t126 + mrSges(5,2) + t191 * t179 - t200 + m(7) * (t193 * t179 - t159) + m(6) * (-t170 * t244 - t159); m(5) + m(6) * t208 + m(7) * (t209 * t170 + t172); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t198 * mrSges(7,3) + t10 * t235 + t3 * t135 + t9 * t234 + t33 * t238 + t32 * t239 + t57 * t241 + t19 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t198 + t248) * pkin(12); t47 * mrSges(6,1) - t48 * mrSges(6,2) + t194 * mrSges(7,3) + t118 * t241 + t45 * t135 + t43 * t234 + t44 * t235 + t88 * t238 + t87 * t239 + t65 + (-m(7) * t45 - t52) * pkin(5) + (m(7) * t194 + t249) * pkin(12); t109 * t235 + t108 * t234 - pkin(5) * t126 + t168 + (m(7) * t193 + t191) * pkin(12) + (t139 * t236 + t141 * t234 - t244 * t247) * t183 + t193 * mrSges(7,3) + (mrSges(6,2) * t244 - Ifges(6,6) + t241) * t179; t247 * t183 + (-mrSges(6,2) + (m(7) * pkin(12) + mrSges(7,3)) * t209) * t179; Ifges(6,3) + m(7) * (t209 * pkin(12) ^ 2 + pkin(5) ^ 2) + t178 * t141 + t182 * t139 - 0.2e1 * pkin(5) * t135 + 0.2e1 * t209 * pkin(12) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t22 - mrSges(7,2) * t23 + t42; mrSges(7,1) * t95 - mrSges(7,2) * t96 + t107; -t196 * t179; -pkin(12) * t196 + t137; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
