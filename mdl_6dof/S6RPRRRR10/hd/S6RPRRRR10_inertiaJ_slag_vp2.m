% Calculate joint inertia matrix for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:28:01
% EndTime: 2019-03-09 07:28:07
% DurationCPUTime: 2.20s
% Computational Cost: add. (6672->387), mult. (17463->580), div. (0->0), fcn. (20184->14), ass. (0->152)
t161 = sin(qJ(6));
t165 = cos(qJ(6));
t127 = -mrSges(7,1) * t165 + mrSges(7,2) * t161;
t228 = -mrSges(6,1) + t127;
t156 = sin(pkin(7));
t168 = cos(qJ(3));
t195 = t156 * t168;
t159 = cos(pkin(7));
t163 = sin(qJ(4));
t167 = cos(qJ(4));
t164 = sin(qJ(3));
t196 = t156 * t164;
t116 = t159 * t167 - t163 * t196;
t117 = t159 * t163 + t167 * t196;
t162 = sin(qJ(5));
t166 = cos(qJ(5));
t80 = t116 * t162 + t117 * t166;
t71 = -t161 * t80 - t165 * t195;
t72 = -t161 * t195 + t165 * t80;
t179 = -t161 * t71 + t165 * t72;
t219 = -pkin(11) - pkin(10);
t135 = t219 * t167;
t185 = t219 * t163;
t107 = -t166 * t135 + t162 * t185;
t123 = t162 * t163 - t166 * t167;
t124 = t162 * t167 + t163 * t166;
t145 = -pkin(4) * t167 - pkin(3);
t87 = pkin(5) * t123 - pkin(12) * t124 + t145;
t58 = -t107 * t161 + t165 * t87;
t59 = t107 * t165 + t161 * t87;
t180 = -t161 * t58 + t165 * t59;
t160 = cos(pkin(6));
t157 = sin(pkin(6));
t158 = cos(pkin(13));
t194 = t157 * t158;
t112 = -t156 * t194 + t159 * t160;
t155 = sin(pkin(13));
t109 = (-pkin(9) * t155 * t156 - pkin(2) * t158 - pkin(1)) * t157;
t214 = pkin(1) * t160;
t115 = qJ(2) * t194 + t155 * t214;
t184 = t159 * t194;
t86 = (t156 * t160 + t184) * pkin(9) + t115;
t139 = t158 * t214;
t197 = t155 * t157;
t96 = t160 * pkin(2) + t139 + (-pkin(9) * t159 - qJ(2)) * t197;
t52 = -t164 * t86 + (t109 * t156 + t159 * t96) * t168;
t44 = -t112 * pkin(3) - t52;
t193 = t159 * t164;
t95 = t160 * t196 + (t155 * t168 + t158 * t193) * t157;
t65 = t112 * t167 - t163 * t95;
t28 = -pkin(4) * t65 + t44;
t66 = t112 * t163 + t167 * t95;
t49 = t162 * t66 - t166 * t65;
t50 = t162 * t65 + t166 * t66;
t15 = pkin(5) * t49 - pkin(12) * t50 + t28;
t60 = t159 * t109 - t156 * t96;
t94 = -t160 * t195 + t164 * t197 - t168 * t184;
t42 = pkin(3) * t94 - pkin(10) * t95 + t60;
t53 = t109 * t196 + t168 * t86 + t96 * t193;
t45 = pkin(10) * t112 + t53;
t25 = -t163 * t45 + t167 * t42;
t17 = pkin(4) * t94 - pkin(11) * t66 + t25;
t26 = t163 * t42 + t167 * t45;
t20 = pkin(11) * t65 + t26;
t9 = t162 * t17 + t166 * t20;
t6 = pkin(12) * t94 + t9;
t2 = t15 * t165 - t161 * t6;
t3 = t15 * t161 + t165 * t6;
t182 = -t161 * t2 + t165 * t3;
t227 = Ifges(5,5) * t66 + Ifges(5,6) * t65 + Ifges(5,3) * t94;
t78 = -t166 * t116 + t117 * t162;
t226 = t78 ^ 2;
t225 = 2 * mrSges(3,1);
t105 = -t135 * t162 - t166 * t185;
t224 = t105 ^ 2;
t223 = 0.2e1 * t105;
t222 = 0.2e1 * t160;
t33 = -t161 * t50 + t165 * t94;
t221 = t33 / 0.2e1;
t129 = Ifges(7,5) * t161 + Ifges(7,6) * t165;
t218 = t129 / 0.2e1;
t207 = Ifges(7,4) * t165;
t132 = Ifges(7,1) * t161 + t207;
t217 = t132 / 0.2e1;
t216 = t161 / 0.2e1;
t215 = t165 / 0.2e1;
t213 = pkin(12) * t161;
t212 = pkin(12) * t165;
t34 = t161 * t94 + t165 * t50;
t19 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t36 = mrSges(6,1) * t94 - mrSges(6,3) * t50;
t209 = t19 - t36;
t208 = Ifges(7,4) * t161;
t206 = t105 * t78;
t201 = t124 * t161;
t200 = t124 * t165;
t143 = pkin(4) * t162 + pkin(12);
t199 = t143 * t161;
t198 = t143 * t165;
t192 = Ifges(6,5) * t124 - Ifges(6,6) * t123;
t191 = Ifges(5,5) * t163 + Ifges(5,6) * t167;
t190 = t161 ^ 2 + t165 ^ 2;
t189 = t163 ^ 2 + t167 ^ 2;
t12 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t49;
t188 = Ifges(6,5) * t50 - Ifges(6,6) * t49 + Ifges(6,3) * t94;
t187 = Ifges(4,5) * t95 - Ifges(4,6) * t94 + Ifges(4,3) * t112;
t130 = Ifges(7,2) * t165 + t208;
t186 = t165 * t130 + t161 * t132 + Ifges(6,3);
t183 = t190 * t143;
t181 = mrSges(7,1) * t161 + mrSges(7,2) * t165;
t8 = -t162 * t20 + t166 * t17;
t178 = 0.2e1 * mrSges(7,3) * t190;
t176 = -t116 * t163 + t117 * t167;
t74 = Ifges(7,5) * t200 - Ifges(7,6) * t201 + Ifges(7,3) * t123;
t175 = (mrSges(6,1) * t166 - mrSges(6,2) * t162) * pkin(4);
t174 = -t80 * mrSges(6,2) + t179 * mrSges(7,3) + t228 * t78;
t13 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t49;
t14 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t49;
t5 = -pkin(5) * t94 - t8;
t173 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + t182 * mrSges(7,3) + t5 * t127 + t13 * t215 + t130 * t221 + t14 * t216 + t34 * t217 + t49 * t218 + t188;
t75 = Ifges(7,6) * t123 + (-Ifges(7,2) * t161 + t207) * t124;
t76 = Ifges(7,5) * t123 + (Ifges(7,1) * t165 - t208) * t124;
t172 = -t107 * mrSges(6,2) + t75 * t215 + t76 * t216 - t130 * t201 / 0.2e1 + t200 * t217 + t123 * t218 + t192 + t228 * t105 + t180 * mrSges(7,3);
t150 = t156 ^ 2;
t144 = -pkin(4) * t166 - pkin(5);
t141 = t150 * t168 ^ 2;
t137 = mrSges(3,2) * t197;
t133 = Ifges(5,1) * t163 + Ifges(5,4) * t167;
t131 = Ifges(5,4) * t163 + Ifges(5,2) * t167;
t128 = -mrSges(5,1) * t167 + mrSges(5,2) * t163;
t114 = -qJ(2) * t197 + t139;
t100 = Ifges(6,1) * t124 - Ifges(6,4) * t123;
t99 = Ifges(6,4) * t124 - Ifges(6,2) * t123;
t98 = mrSges(6,1) * t123 + mrSges(6,2) * t124;
t93 = mrSges(7,1) * t123 - mrSges(7,3) * t200;
t92 = -mrSges(7,2) * t123 - mrSges(7,3) * t201;
t85 = t181 * t124;
t70 = mrSges(4,1) * t112 - mrSges(4,3) * t95;
t69 = -mrSges(4,2) * t112 - mrSges(4,3) * t94;
t57 = mrSges(4,1) * t94 + mrSges(4,2) * t95;
t55 = mrSges(5,1) * t94 - mrSges(5,3) * t66;
t54 = -mrSges(5,2) * t94 + mrSges(5,3) * t65;
t51 = -mrSges(5,1) * t65 + mrSges(5,2) * t66;
t38 = Ifges(5,1) * t66 + Ifges(5,4) * t65 + Ifges(5,5) * t94;
t37 = Ifges(5,4) * t66 + Ifges(5,2) * t65 + Ifges(5,6) * t94;
t35 = -mrSges(6,2) * t94 - mrSges(6,3) * t49;
t27 = mrSges(6,1) * t49 + mrSges(6,2) * t50;
t24 = Ifges(6,1) * t50 - Ifges(6,4) * t49 + Ifges(6,5) * t94;
t23 = Ifges(6,4) * t50 - Ifges(6,2) * t49 + Ifges(6,6) * t94;
t22 = mrSges(7,1) * t49 - mrSges(7,3) * t34;
t21 = -mrSges(7,2) * t49 + mrSges(7,3) * t33;
t1 = [t33 * t13 + t34 * t14 + 0.2e1 * t9 * t35 + 0.2e1 * t8 * t36 + 0.2e1 * t28 * t27 + 0.2e1 * t3 * t21 + 0.2e1 * t2 * t22 + 0.2e1 * t5 * t19 + t65 * t37 + t66 * t38 + 0.2e1 * t53 * t69 + 0.2e1 * t52 * t70 + 0.2e1 * t60 * t57 + t50 * t24 + 0.2e1 * t44 * t51 + 0.2e1 * t26 * t54 + 0.2e1 * t25 * t55 + m(3) * (pkin(1) ^ 2 * t157 ^ 2 + t114 ^ 2 + t115 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + m(6) * (t28 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2 + t44 ^ 2) + m(4) * (t52 ^ 2 + t53 ^ 2 + t60 ^ 2) + (t12 - t23) * t49 + (-0.2e1 * Ifges(4,4) * t95 + Ifges(4,2) * t94 - Ifges(4,6) * t112 + t188 + t227) * t94 + t95 * (Ifges(4,1) * t95 + Ifges(4,5) * t112) + Ifges(2,3) + t112 * t187 + (-0.2e1 * mrSges(3,2) * t115 + Ifges(3,3) * t160 + t114 * t225) * t160 + (-0.2e1 * pkin(1) * t137 + (-0.2e1 * mrSges(3,3) * t114 + Ifges(3,1) * t197 + Ifges(3,5) * t222) * t155 + (0.2e1 * t115 * mrSges(3,3) + Ifges(3,6) * t222 + (0.2e1 * Ifges(3,4) * t155 + Ifges(3,2) * t158 + pkin(1) * t225) * t157) * t158) * t157; t116 * t55 + t117 * t54 + t159 * t57 + t72 * t21 + t71 * t22 + t80 * t35 + t137 + t209 * t78 + (-m(3) * pkin(1) - mrSges(3,1) * t158) * t157 + (t164 * t69 + (-t27 - t51 + t70) * t168) * t156 + m(7) * (t2 * t71 + t3 * t72 + t5 * t78) + m(6) * (-t195 * t28 - t78 * t8 + t80 * t9) + m(5) * (t116 * t25 + t117 * t26 - t195 * t44) + m(4) * (t159 * t60 + (t164 * t53 + t168 * t52) * t156); m(3) + m(7) * (t71 ^ 2 + t72 ^ 2 + t226) + m(6) * (t80 ^ 2 + t141 + t226) + m(5) * (t116 ^ 2 + t117 ^ 2 + t141) + m(4) * (t150 * t164 ^ 2 + t159 ^ 2 + t141); t145 * t27 + t44 * t128 + t65 * t131 / 0.2e1 + t66 * t133 / 0.2e1 + t28 * t98 + t50 * t100 / 0.2e1 + t107 * t35 + t3 * t92 + t2 * t93 + t5 * t85 + t34 * t76 / 0.2e1 + t58 * t22 + t59 * t21 - pkin(3) * t51 + t52 * mrSges(4,1) - t53 * mrSges(4,2) + t187 + (t192 + t191) * t94 / 0.2e1 + (-t99 / 0.2e1 + t74 / 0.2e1) * t49 + (t37 / 0.2e1 + t26 * mrSges(5,3) + pkin(10) * t54) * t167 + (t38 / 0.2e1 - t25 * mrSges(5,3) - pkin(10) * t55) * t163 + m(5) * (-pkin(3) * t44 + (-t163 * t25 + t167 * t26) * pkin(10)) + m(7) * (t105 * t5 + t2 * t58 + t3 * t59) + m(6) * (-t105 * t8 + t107 * t9 + t145 * t28) + (t12 / 0.2e1 - t23 / 0.2e1 - t9 * mrSges(6,3)) * t123 + t75 * t221 + t209 * t105 + (t24 / 0.2e1 + t14 * t215 - t161 * t13 / 0.2e1 - t8 * mrSges(6,3)) * t124; t71 * t93 + t72 * t92 + t78 * t85 + (-t123 * t80 + t124 * t78) * mrSges(6,3) + t176 * mrSges(5,3) + (-t164 * mrSges(4,2) + (mrSges(4,1) - t128 - t98) * t168) * t156 + m(7) * (t58 * t71 + t59 * t72 + t206) + m(6) * (t107 * t80 - t145 * t195 + t206) + m(5) * (pkin(3) * t195 + pkin(10) * t176); -0.2e1 * pkin(3) * t128 + t85 * t223 + t167 * t131 + t163 * t133 + 0.2e1 * t145 * t98 + 0.2e1 * t58 * t93 + 0.2e1 * t59 * t92 + Ifges(4,3) + 0.2e1 * t189 * pkin(10) * mrSges(5,3) + (-0.2e1 * mrSges(6,3) * t107 + t74 - t99) * t123 + m(7) * (t58 ^ 2 + t59 ^ 2 + t224) + m(6) * (t107 ^ 2 + t145 ^ 2 + t224) + m(5) * (pkin(10) ^ 2 * t189 + pkin(3) ^ 2) + (mrSges(6,3) * t223 - t161 * t75 + t165 * t76 + t100) * t124; t25 * mrSges(5,1) - t26 * mrSges(5,2) + t144 * t19 + t173 + t21 * t198 - t22 * t199 + m(7) * (t143 * t182 + t144 * t5) + (m(6) * (t162 * t9 + t166 * t8) + t166 * t36 + t162 * t35) * pkin(4) + t227; t116 * mrSges(5,1) - t117 * mrSges(5,2) + m(7) * (t143 * t179 + t144 * t78) + m(6) * (t162 * t80 - t166 * t78) * pkin(4) + t174; t144 * t85 + t172 + (-mrSges(5,1) * t163 - mrSges(5,2) * t167) * pkin(10) + t92 * t198 - t93 * t199 + m(7) * (t144 * t105 + t143 * t180) + (m(6) * (-t105 * t166 + t107 * t162) + (-t123 * t162 - t124 * t166) * mrSges(6,3)) * pkin(4) + t191; 0.2e1 * t144 * t127 + Ifges(5,3) + 0.2e1 * t175 + t143 * t178 + m(7) * (t143 ^ 2 * t190 + t144 ^ 2) + m(6) * (t162 ^ 2 + t166 ^ 2) * pkin(4) ^ 2 + t186; -pkin(5) * t19 + t173 + t21 * t212 - t22 * t213 + m(7) * (-pkin(5) * t5 + pkin(12) * t182); m(7) * (-pkin(5) * t78 + pkin(12) * t179) + t174; -pkin(5) * t85 + t172 + t92 * t212 - t93 * t213 + m(7) * (-pkin(5) * t105 + pkin(12) * t180); m(7) * (-pkin(5) * t144 + pkin(12) * t183) + (-pkin(5) + t144) * t127 + t175 + (pkin(12) * t190 + t183) * mrSges(7,3) + t186; -0.2e1 * pkin(5) * t127 + m(7) * (pkin(12) ^ 2 * t190 + pkin(5) ^ 2) + pkin(12) * t178 + t186; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t71 - mrSges(7,2) * t72; mrSges(7,1) * t58 - mrSges(7,2) * t59 + t74; -t143 * t181 + t129; -pkin(12) * t181 + t129; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
