% Calculate joint inertia matrix for
% S6RRRPRP10
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:42
% EndTime: 2019-03-09 17:29:47
% DurationCPUTime: 1.95s
% Computational Cost: add. (3387->436), mult. (7700->592), div. (0->0), fcn. (8372->10), ass. (0->153)
t214 = 2 * pkin(9);
t170 = cos(pkin(6));
t172 = sin(qJ(3));
t174 = cos(qJ(3));
t168 = sin(pkin(6));
t173 = sin(qJ(2));
t195 = t168 * t173;
t115 = -t170 * t174 + t172 * t195;
t171 = sin(qJ(5));
t207 = cos(qJ(5));
t116 = t170 * t172 + t174 * t195;
t167 = sin(pkin(11));
t169 = cos(pkin(11));
t175 = cos(qJ(2));
t194 = t168 * t175;
t78 = -t116 * t167 - t169 * t194;
t79 = t116 * t169 - t167 * t194;
t39 = t171 * t79 - t207 * t78;
t40 = t171 * t78 + t207 * t79;
t10 = Ifges(7,4) * t40 + Ifges(7,2) * t115 + Ifges(7,6) * t39;
t9 = Ifges(6,5) * t40 - Ifges(6,6) * t39 + Ifges(6,3) * t115;
t213 = t10 + t9;
t212 = -m(7) * pkin(5) - mrSges(7,1);
t27 = Ifges(5,1) * t79 + Ifges(5,4) * t78 + Ifges(5,5) * t115;
t211 = t27 / 0.2e1;
t200 = Ifges(5,4) * t169;
t104 = -Ifges(5,6) * t174 + (-Ifges(5,2) * t167 + t200) * t172;
t210 = t104 / 0.2e1;
t201 = Ifges(5,4) * t167;
t105 = -Ifges(5,5) * t174 + (Ifges(5,1) * t169 - t201) * t172;
t209 = t105 / 0.2e1;
t139 = Ifges(5,1) * t167 + t200;
t208 = t139 / 0.2e1;
t148 = pkin(8) * t195;
t206 = pkin(1) * t175;
t106 = t148 + (-pkin(2) - t206) * t170;
t45 = pkin(3) * t115 - qJ(4) * t116 + t106;
t119 = t170 * t173 * pkin(1) + pkin(8) * t194;
t107 = pkin(9) * t170 + t119;
t108 = (-pkin(2) * t175 - pkin(9) * t173 - pkin(1)) * t168;
t52 = t174 * t107 + t172 * t108;
t46 = -qJ(4) * t194 + t52;
t19 = t167 * t45 + t169 * t46;
t15 = pkin(10) * t78 + t19;
t18 = -t167 * t46 + t169 * t45;
t7 = pkin(4) * t115 - pkin(10) * t79 + t18;
t4 = t207 * t15 + t171 * t7;
t205 = pkin(9) * t174;
t162 = t172 * pkin(9);
t203 = -Ifges(6,3) - Ifges(7,2);
t202 = pkin(10) + qJ(4);
t134 = -pkin(3) * t174 - qJ(4) * t172 - pkin(2);
t127 = t169 * t134;
t193 = t169 * t172;
t66 = -pkin(10) * t193 + t127 + (-pkin(9) * t167 - pkin(4)) * t174;
t196 = t167 * t172;
t95 = t167 * t134 + t169 * t205;
t80 = -pkin(10) * t196 + t95;
t31 = t171 * t66 + t207 * t80;
t118 = t170 * t206 - t148;
t199 = t118 * mrSges(3,1);
t198 = t119 * mrSges(3,2);
t130 = t167 * t207 + t171 * t169;
t109 = t130 * t172;
t178 = -t171 * t167 + t169 * t207;
t110 = t178 * t172;
t197 = Ifges(7,4) * t110 + Ifges(7,6) * t109;
t90 = t174 * mrSges(7,1) + t110 * mrSges(7,2);
t192 = Ifges(6,5) * t110 - Ifges(6,6) * t109;
t191 = -Ifges(4,5) * t116 + Ifges(4,6) * t115;
t73 = Ifges(6,5) * t130 + Ifges(6,6) * t178;
t74 = Ifges(7,4) * t130 - Ifges(7,6) * t178;
t117 = mrSges(5,1) * t196 + mrSges(5,2) * t193;
t190 = Ifges(4,5) * t172 + Ifges(4,6) * t174;
t133 = pkin(4) * t196 + t162;
t189 = t167 ^ 2 + t169 ^ 2;
t136 = t202 * t169;
t180 = t202 * t167;
t82 = t136 * t171 + t180 * t207;
t84 = t136 * t207 - t171 * t180;
t188 = t82 ^ 2 + t84 ^ 2;
t11 = Ifges(6,4) * t40 - Ifges(6,2) * t39 + Ifges(6,6) * t115;
t8 = Ifges(7,5) * t40 + Ifges(7,6) * t115 + Ifges(7,3) * t39;
t187 = t8 / 0.2e1 - t11 / 0.2e1;
t12 = Ifges(7,1) * t40 + Ifges(7,4) * t115 + Ifges(7,5) * t39;
t13 = Ifges(6,1) * t40 - Ifges(6,4) * t39 + Ifges(6,5) * t115;
t186 = t12 / 0.2e1 + t13 / 0.2e1;
t53 = Ifges(7,5) * t110 - Ifges(7,6) * t174 + Ifges(7,3) * t109;
t56 = Ifges(6,4) * t110 - Ifges(6,2) * t109 - Ifges(6,6) * t174;
t185 = t53 / 0.2e1 - t56 / 0.2e1;
t57 = Ifges(7,1) * t110 - Ifges(7,4) * t174 + Ifges(7,5) * t109;
t58 = Ifges(6,1) * t110 - Ifges(6,4) * t109 - Ifges(6,5) * t174;
t184 = t57 / 0.2e1 + t58 / 0.2e1;
t72 = Ifges(7,5) * t130 - Ifges(7,3) * t178;
t75 = Ifges(6,4) * t130 + Ifges(6,2) * t178;
t183 = t72 / 0.2e1 - t75 / 0.2e1;
t76 = Ifges(7,1) * t130 - Ifges(7,5) * t178;
t77 = Ifges(6,1) * t130 + Ifges(6,4) * t178;
t182 = t76 / 0.2e1 + t77 / 0.2e1;
t181 = Ifges(3,5) * t195 + Ifges(3,6) * t194 + Ifges(3,3) * t170;
t155 = -pkin(4) * t169 - pkin(3);
t41 = -t78 * mrSges(5,1) + t79 * mrSges(5,2);
t17 = t39 * mrSges(6,1) + t40 * mrSges(6,2);
t16 = t39 * mrSges(7,1) - t40 * mrSges(7,3);
t62 = t109 * mrSges(6,1) + t110 * mrSges(6,2);
t23 = -t115 * mrSges(7,1) + t40 * mrSges(7,2);
t61 = t109 * mrSges(7,1) - t110 * mrSges(7,3);
t135 = -t169 * mrSges(5,1) + t167 * mrSges(5,2);
t71 = -mrSges(6,1) * t178 + t130 * mrSges(6,2);
t70 = -mrSges(7,1) * t178 - t130 * mrSges(7,3);
t51 = -t172 * t107 + t108 * t174;
t179 = Ifges(5,5) * t167 / 0.2e1 + Ifges(5,6) * t169 / 0.2e1 + t73 / 0.2e1 + t74 / 0.2e1;
t47 = pkin(3) * t194 - t51;
t3 = -t171 * t15 + t207 * t7;
t30 = -t171 * t80 + t207 * t66;
t24 = -pkin(4) * t78 + t47;
t177 = pkin(9) ^ 2;
t166 = t174 ^ 2;
t165 = t172 ^ 2;
t161 = t165 * t177;
t142 = Ifges(4,1) * t172 + Ifges(4,4) * t174;
t141 = Ifges(4,4) * t172 + Ifges(4,2) * t174;
t140 = -mrSges(4,1) * t174 + mrSges(4,2) * t172;
t138 = Ifges(5,2) * t169 + t201;
t132 = -mrSges(5,1) * t174 - mrSges(5,3) * t193;
t131 = mrSges(5,2) * t174 - mrSges(5,3) * t196;
t103 = -Ifges(5,3) * t174 + (Ifges(5,5) * t169 - Ifges(5,6) * t167) * t172;
t94 = -t167 * t205 + t127;
t89 = -mrSges(6,1) * t174 - mrSges(6,3) * t110;
t88 = mrSges(6,2) * t174 - mrSges(6,3) * t109;
t87 = -mrSges(7,2) * t109 - mrSges(7,3) * t174;
t86 = -mrSges(4,1) * t194 - mrSges(4,3) * t116;
t85 = mrSges(4,2) * t194 - mrSges(4,3) * t115;
t65 = mrSges(4,1) * t115 + mrSges(4,2) * t116;
t64 = -pkin(5) * t178 - qJ(6) * t130 + t155;
t60 = Ifges(4,1) * t116 - Ifges(4,4) * t115 - Ifges(4,5) * t194;
t59 = Ifges(4,4) * t116 - Ifges(4,2) * t115 - Ifges(4,6) * t194;
t55 = -Ifges(7,2) * t174 + t197;
t54 = -Ifges(6,3) * t174 + t192;
t50 = mrSges(5,1) * t115 - mrSges(5,3) * t79;
t49 = -mrSges(5,2) * t115 + mrSges(5,3) * t78;
t48 = pkin(5) * t109 - qJ(6) * t110 + t133;
t29 = t174 * pkin(5) - t30;
t28 = -qJ(6) * t174 + t31;
t26 = Ifges(5,4) * t79 + Ifges(5,2) * t78 + Ifges(5,6) * t115;
t25 = Ifges(5,5) * t79 + Ifges(5,6) * t78 + Ifges(5,3) * t115;
t22 = mrSges(6,1) * t115 - mrSges(6,3) * t40;
t21 = -mrSges(6,2) * t115 - mrSges(6,3) * t39;
t20 = -mrSges(7,2) * t39 + mrSges(7,3) * t115;
t5 = pkin(5) * t39 - qJ(6) * t40 + t24;
t2 = -t115 * pkin(5) - t3;
t1 = qJ(6) * t115 + t4;
t6 = [(-t59 + t25 + t213) * t115 + ((-0.2e1 * t118 * mrSges(3,3) + Ifges(3,5) * t170 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t173) * t168) * t173 + (0.2e1 * t119 * mrSges(3,3) + Ifges(3,6) * t170 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t173 + (Ifges(3,2) + Ifges(4,3)) * t175) * t168 + t191) * t175) * t168 + 0.2e1 * t19 * t49 + 0.2e1 * t18 * t50 + 0.2e1 * t47 * t41 + 0.2e1 * t1 * t20 + 0.2e1 * t4 * t21 + 0.2e1 * t3 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t24 * t17 + 0.2e1 * t5 * t16 + m(3) * (pkin(1) ^ 2 * t168 ^ 2 + t118 ^ 2 + t119 ^ 2) + m(4) * (t106 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2 + t47 ^ 2) + m(6) * (t24 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (t12 + t13) * t40 + (t8 - t11) * t39 + (t181 - 0.2e1 * t198 + 0.2e1 * t199) * t170 + t116 * t60 + Ifges(2,3) + 0.2e1 * t106 * t65 + 0.2e1 * t52 * t85 + 0.2e1 * t51 * t86 + t78 * t26 + t79 * t27; t184 * t40 + t185 * t39 + t186 * t110 + t187 * t109 - t190 * t194 / 0.2e1 + t5 * t61 + t24 * t62 - pkin(2) * t65 + t48 * t16 + t30 * t22 + t31 * t21 + t28 * t20 + t29 * t23 + m(6) * (t133 * t24 + t3 * t30 + t31 * t4) + m(7) * (t1 * t28 + t2 * t29 + t48 * t5) + m(4) * (-pkin(2) * t106 + (-t51 * t172 + t52 * t174) * pkin(9)) + (pkin(9) * t85 + t52 * mrSges(4,3) + t59 / 0.2e1 - t25 / 0.2e1 - t9 / 0.2e1 - t10 / 0.2e1) * t174 + (-t141 / 0.2e1 + t54 / 0.2e1 + t55 / 0.2e1 + t103 / 0.2e1) * t115 + m(5) * (t47 * t162 + t18 * t94 + t19 * t95) + (-t167 * t26 / 0.2e1 + t169 * t211 - t51 * mrSges(4,3) + t60 / 0.2e1 + (-t86 + t41) * pkin(9)) * t172 + t106 * t140 + t116 * t142 / 0.2e1 + t19 * t131 + t18 * t132 + t133 * t17 + t47 * t117 + t79 * t209 + t78 * t210 - t198 + t199 + t181 + t1 * t87 + t4 * t88 + t3 * t89 + t2 * t90 + t94 * t50 + t95 * t49; -0.2e1 * pkin(2) * t140 + 0.2e1 * t95 * t131 + 0.2e1 * t94 * t132 + 0.2e1 * t133 * t62 + 0.2e1 * t28 * t87 + 0.2e1 * t29 * t90 + 0.2e1 * t30 * t89 + 0.2e1 * t31 * t88 + 0.2e1 * t48 * t61 + Ifges(3,3) + (t57 + t58) * t110 + (t53 - t56) * t109 + (t165 + t166) * mrSges(4,3) * t214 + (-t103 - t54 - t55 + t141) * t174 + (-t104 * t167 + t105 * t169 + t117 * t214 + t142) * t172 + m(4) * (pkin(2) ^ 2 + t166 * t177 + t161) + m(5) * (t94 ^ 2 + t95 ^ 2 + t161) + m(6) * (t133 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(7) * (t28 ^ 2 + t29 ^ 2 + t48 ^ 2); -t191 - Ifges(4,3) * t194 + t182 * t40 + t183 * t39 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t186) * t130 - (-t1 * mrSges(7,2) - t4 * mrSges(6,3) + t187) * t178 + t64 * t16 + t5 * t70 + t51 * mrSges(4,1) - t52 * mrSges(4,2) - pkin(3) * t41 + m(6) * (t155 * t24 - t3 * t82 + t4 * t84) + m(7) * (t1 * t84 + t2 * t82 + t5 * t64) + t179 * t115 + m(5) * (-pkin(3) * t47 + (-t18 * t167 + t19 * t169) * qJ(4)) + t79 * t208 + t155 * t17 + t47 * t135 + t78 * t138 / 0.2e1 + (-t18 * mrSges(5,3) - qJ(4) * t50 + t211) * t167 + (qJ(4) * t49 + t19 * mrSges(5,3) + t26 / 0.2e1) * t169 + t24 * t71 + (t23 - t22) * t82 + (t20 + t21) * t84; -pkin(3) * t117 + t133 * t71 + t155 * t62 + t48 * t70 + t64 * t61 + (t87 + t88) * t84 + (-t89 + t90) * t82 + t182 * t110 + t183 * t109 + (-mrSges(4,1) + t135) * t162 + (t95 * mrSges(5,3) + qJ(4) * t131 + t172 * t208 + t210) * t169 + (-t172 * t138 / 0.2e1 - qJ(4) * t132 - t94 * mrSges(5,3) + t209) * t167 + m(5) * (-pkin(3) * t162 + (-t167 * t94 + t169 * t95) * qJ(4)) + m(6) * (t133 * t155 - t30 * t82 + t31 * t84) + m(7) * (t28 * t84 + t29 * t82 + t48 * t64) + (-pkin(9) * mrSges(4,2) - t179) * t174 + (t29 * mrSges(7,2) - t30 * mrSges(6,3) + t184) * t130 - (-t28 * mrSges(7,2) - t31 * mrSges(6,3) + t185) * t178 + t190; -0.2e1 * pkin(3) * t135 + t169 * t138 + t167 * t139 + 0.2e1 * t155 * t71 + 0.2e1 * t64 * t70 + Ifges(4,3) + m(6) * (t155 ^ 2 + t188) + m(7) * (t64 ^ 2 + t188) + m(5) * (t189 * qJ(4) ^ 2 + pkin(3) ^ 2) + (t76 + t77) * t130 - (t72 - t75) * t178 + 0.2e1 * t189 * qJ(4) * mrSges(5,3) + 0.2e1 * (t130 * t82 + t178 * t84) * (mrSges(7,2) + mrSges(6,3)); m(5) * t47 + m(6) * t24 + m(7) * t5 + t16 + t17 + t41; m(5) * t162 + m(6) * t133 + m(7) * t48 + t117 + t61 + t62; -m(5) * pkin(3) + m(6) * t155 + m(7) * t64 + t135 + t70 + t71; m(5) + m(6) + m(7); m(7) * (-pkin(5) * t2 + qJ(6) * t1) - pkin(5) * t23 + t1 * mrSges(7,3) + qJ(6) * t20 - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t213; -t31 * mrSges(6,2) + t30 * mrSges(6,1) - t29 * mrSges(7,1) - pkin(5) * t90 + t28 * mrSges(7,3) + qJ(6) * t87 + m(7) * (-pkin(5) * t29 + qJ(6) * t28) + t203 * t174 + t192 + t197; (-pkin(5) * t130 + qJ(6) * t178) * mrSges(7,2) + t74 + t73 + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t84 + (-mrSges(6,1) + t212) * t82; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t203; m(7) * t2 + t23; m(7) * t29 + t90; m(7) * t82 + t130 * mrSges(7,2); 0; t212; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
