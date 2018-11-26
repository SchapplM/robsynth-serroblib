% Calculate joint inertia matrix for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:14:36
% EndTime: 2018-11-23 18:14:38
% DurationCPUTime: 1.57s
% Computational Cost: add. (2158->310), mult. (4016->418), div. (0->0), fcn. (4153->8), ass. (0->126)
t217 = Ifges(5,1) + Ifges(6,1);
t216 = Ifges(6,4) + Ifges(5,5);
t215 = Ifges(6,2) + Ifges(5,3);
t136 = cos(qJ(4));
t132 = sin(qJ(4));
t181 = Ifges(6,5) * t132;
t183 = Ifges(5,4) * t132;
t133 = sin(qJ(3));
t134 = sin(qJ(2));
t137 = cos(qJ(3));
t138 = cos(qJ(2));
t92 = t133 * t134 - t137 * t138;
t94 = t133 * t138 + t134 * t137;
t214 = (t136 * t217 + t181 - t183) * t94 + t216 * t92;
t180 = Ifges(6,5) * t136;
t182 = Ifges(5,4) * t136;
t213 = t132 * t217 - t180 + t182;
t212 = t132 ^ 2 + t136 ^ 2;
t211 = (Ifges(5,6) - Ifges(6,6)) * t136 + t216 * t132;
t116 = -pkin(2) * t138 - pkin(1);
t47 = pkin(3) * t92 - pkin(9) * t94 + t116;
t199 = -pkin(8) - pkin(7);
t109 = t199 * t138;
t163 = t199 * t134;
t66 = -t137 * t109 + t133 * t163;
t55 = t132 * t66;
t21 = t136 * t47 - t55;
t22 = t132 * t47 + t136 * t66;
t154 = -t132 * t21 + t136 * t22;
t13 = t92 * qJ(5) + t22;
t14 = -pkin(4) * t92 - t21;
t155 = t13 * t136 + t132 * t14;
t210 = t136 * pkin(4) + t132 * qJ(5);
t198 = pkin(9) - pkin(10);
t107 = t198 * t132;
t108 = t198 * t136;
t131 = sin(qJ(6));
t135 = cos(qJ(6));
t62 = t107 * t135 - t108 * t131;
t65 = t107 * t131 + t108 * t135;
t209 = t62 * mrSges(7,1) - t65 * mrSges(7,2);
t114 = pkin(2) * t133 + pkin(9);
t185 = -pkin(10) + t114;
t86 = t185 * t132;
t87 = t185 * t136;
t43 = -t131 * t87 + t135 * t86;
t44 = t131 * t86 + t135 * t87;
t208 = t43 * mrSges(7,1) - t44 * mrSges(7,2);
t207 = (mrSges(5,3) + mrSges(6,2)) * t212;
t174 = t136 * t94;
t176 = t132 * t94;
t206 = Ifges(6,6) * t176 + t174 * t216 + t215 * t92;
t63 = -t109 * t133 - t137 * t163;
t205 = t63 ^ 2;
t153 = t131 * t132 + t135 * t136;
t93 = -t131 * t136 + t132 * t135;
t52 = mrSges(7,1) * t153 + mrSges(7,2) * t93;
t204 = 0.2e1 * t52;
t203 = 0.2e1 * t63;
t101 = -mrSges(6,1) * t136 - mrSges(6,3) * t132;
t202 = 0.2e1 * t101;
t201 = 0.2e1 * t116;
t139 = -pkin(4) - pkin(5);
t196 = m(6) * t132;
t195 = Ifges(5,6) * t92;
t189 = t153 * mrSges(7,3);
t188 = t93 * mrSges(7,3);
t96 = -qJ(5) * t131 + t135 * t139;
t187 = t96 * mrSges(7,1);
t97 = qJ(5) * t135 + t131 * t139;
t186 = t97 * mrSges(7,2);
t184 = Ifges(7,5) * t93 - Ifges(7,6) * t153;
t119 = t132 * mrSges(6,2);
t173 = t212 * pkin(9) * t114;
t172 = qJ(5) * t136;
t171 = t212 * t114 ^ 2;
t169 = t212 * pkin(9) ^ 2;
t167 = t134 ^ 2 + t138 ^ 2;
t166 = 0.2e1 * mrSges(7,3);
t39 = t93 * t94;
t40 = t153 * t94;
t165 = Ifges(7,5) * t40 + Ifges(7,6) * t39 - Ifges(7,3) * t92;
t164 = pkin(3) + t210;
t115 = -pkin(2) * t137 - pkin(3);
t50 = -t92 * mrSges(6,1) + mrSges(6,2) * t174;
t159 = t132 * mrSges(5,1) + t136 * mrSges(5,2);
t158 = t132 * mrSges(6,1) - t136 * mrSges(6,3);
t157 = t135 * mrSges(7,1) - t131 * mrSges(7,2);
t156 = pkin(4) * t132 - t172;
t76 = t115 - t210;
t152 = -t131 * t189 - t135 * t188 + t119;
t151 = (t137 * mrSges(4,1) - t133 * mrSges(4,2)) * pkin(2);
t4 = t55 + t139 * t92 + (-pkin(10) * t94 - t47) * t136;
t8 = pkin(10) * t176 + t13;
t2 = -t131 * t8 + t135 * t4;
t3 = t131 * t4 + t135 * t8;
t150 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t165;
t48 = -mrSges(5,2) * t92 - mrSges(5,3) * t176;
t49 = mrSges(5,1) * t92 - mrSges(5,3) * t174;
t51 = -mrSges(6,2) * t176 + mrSges(6,3) * t92;
t149 = (t48 + t51) * t136 + (-t49 + t50) * t132;
t103 = -Ifges(6,3) * t136 + t181;
t104 = Ifges(5,2) * t136 + t183;
t53 = Ifges(7,4) * t93 - Ifges(7,2) * t153;
t54 = Ifges(7,1) * t93 - Ifges(7,4) * t153;
t148 = -t153 * t53 + t93 * t54 + Ifges(4,3) + (-t103 + t104) * t136 + t213 * t132;
t146 = 0.2e1 * t207;
t145 = mrSges(6,2) * t172 - pkin(4) * t119 - t96 * t188 - t97 * t189 - t184 + t211;
t144 = -m(6) * t156 - t158 - t159;
t10 = Ifges(7,1) * t40 + Ifges(7,4) * t39 - Ifges(7,5) * t92;
t102 = -mrSges(5,1) * t136 + mrSges(5,2) * t132;
t16 = (t139 * t132 + t172) * t94 - t63;
t23 = t156 * t94 + t63;
t30 = Ifges(6,6) * t92 + (Ifges(6,3) * t132 + t180) * t94;
t31 = t195 + (-Ifges(5,2) * t132 + t182) * t94;
t9 = Ifges(7,4) * t40 + Ifges(7,2) * t39 - Ifges(7,6) * t92;
t143 = -t2 * t188 + t23 * t101 - t153 * t9 / 0.2e1 + t93 * t10 / 0.2e1 + Ifges(4,5) * t94 + t39 * t53 / 0.2e1 + t40 * t54 / 0.2e1 - t66 * mrSges(4,2) + t16 * t52 - t3 * t189 + (t102 - mrSges(4,1)) * t63 + t214 * t132 / 0.2e1 + (-t104 / 0.2e1 + t103 / 0.2e1) * t176 + t213 * t174 / 0.2e1 + (t31 / 0.2e1 - t30 / 0.2e1) * t136 + t154 * mrSges(5,3) + t155 * mrSges(6,2) + (-t184 / 0.2e1 - Ifges(4,6) + t211 / 0.2e1) * t92;
t125 = t136 * pkin(5);
t77 = t125 + t164;
t70 = t125 - t76;
t46 = t159 * t94;
t45 = t158 * t94;
t25 = -mrSges(7,1) * t92 - mrSges(7,3) * t40;
t24 = mrSges(7,2) * t92 + mrSges(7,3) * t39;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t40;
t1 = [t134 * (Ifges(3,1) * t134 + Ifges(3,4) * t138) + t138 * (Ifges(3,4) * t134 + Ifges(3,2) * t138) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t138 + mrSges(3,2) * t134) + t46 * t203 + t40 * t10 + 0.2e1 * t23 * t45 + 0.2e1 * t22 * t48 + 0.2e1 * t21 * t49 + 0.2e1 * t14 * t50 + 0.2e1 * t13 * t51 + 0.2e1 * t3 * t24 + 0.2e1 * t2 * t25 + t39 * t9 + 0.2e1 * t16 * t15 + Ifges(2,3) + 0.2e1 * t167 * pkin(7) * mrSges(3,3) + (mrSges(4,1) * t201 - 0.2e1 * t66 * mrSges(4,3) + Ifges(4,2) * t92 - t165 + t206) * t92 + (mrSges(4,2) * t201 + mrSges(4,3) * t203 + Ifges(4,1) * t94 - 0.2e1 * Ifges(4,4) * t92 + t214 * t136 + (t30 - t31 - t195) * t132) * t94 + m(3) * (t167 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t116 ^ 2 + t66 ^ 2 + t205) + m(5) * (t21 ^ 2 + t22 ^ 2 + t205) + m(6) * (t13 ^ 2 + t14 ^ 2 + t23 ^ 2) + m(7) * (t16 ^ 2 + t2 ^ 2 + t3 ^ 2); t143 + m(6) * (t114 * t155 + t23 * t76) + m(5) * (t114 * t154 + t115 * t63) + Ifges(3,6) * t138 + Ifges(3,5) * t134 + t115 * t46 + t70 * t15 + t76 * t45 + t43 * t25 + t44 * t24 + t149 * t114 + m(7) * (t16 * t70 + t2 * t43 + t3 * t44) + (-t134 * mrSges(3,1) - t138 * mrSges(3,2)) * pkin(7) + (m(4) * (t133 * t66 - t137 * t63) + (-t133 * t92 - t137 * t94) * mrSges(4,3)) * pkin(2); t148 + m(5) * (t115 ^ 2 + t171) + m(6) * (t76 ^ 2 + t171) + m(7) * (t43 ^ 2 + t44 ^ 2 + t70 ^ 2) + m(4) * (t133 ^ 2 + t137 ^ 2) * pkin(2) ^ 2 + 0.2e1 * t115 * t102 + t76 * t202 + t70 * t204 + t146 * t114 + 0.2e1 * t151 + (-t153 * t44 - t43 * t93) * t166 + Ifges(3,3); t143 + m(6) * (pkin(9) * t155 - t164 * t23) + m(5) * (-pkin(3) * t63 + pkin(9) * t154) - t164 * t45 + t77 * t15 + t62 * t25 + t65 * t24 - pkin(3) * t46 + t149 * pkin(9) + m(7) * (t16 * t77 + t2 * t62 + t3 * t65); t148 + m(6) * (-t164 * t76 + t173) + m(5) * (-pkin(3) * t115 + t173) + m(7) * (t43 * t62 + t44 * t65 + t70 * t77) + (t77 + t70) * t52 + (t115 - pkin(3)) * t102 + (-t164 + t76) * t101 + t151 + ((-t43 - t62) * t93 - (t44 + t65) * t153) * mrSges(7,3) + (pkin(9) + t114) * t207; -0.2e1 * pkin(3) * t102 - t164 * t202 + t77 * t204 + (-t153 * t65 - t62 * t93) * t166 + m(7) * (t62 ^ 2 + t65 ^ 2 + t77 ^ 2) + m(6) * (t164 ^ 2 + t169) + m(5) * (pkin(3) ^ 2 + t169) + t146 * pkin(9) + t148; -t150 + m(7) * (t2 * t96 + t3 * t97) + m(6) * (-pkin(4) * t14 + qJ(5) * t13) + t96 * t25 + t97 * t24 - pkin(4) * t50 + qJ(5) * t51 + t13 * mrSges(6,3) - t14 * mrSges(6,1) + t21 * mrSges(5,1) - t22 * mrSges(5,2) - Ifges(5,6) * t176 + t206; t145 + m(7) * (t43 * t96 + t44 * t97) + t144 * t114 - t208; t145 + m(7) * (t62 * t96 + t65 * t97) + t144 * pkin(9) - t209; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t187 + 0.2e1 * t186 + 0.2e1 * qJ(5) * mrSges(6,3) + Ifges(7,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + m(7) * (t96 ^ 2 + t97 ^ 2) + t215; t131 * t24 + t135 * t25 + m(7) * (t131 * t3 + t135 * t2) + m(6) * t14 + t50; m(7) * (t131 * t44 + t135 * t43) + t114 * t196 + t152; m(7) * (t131 * t65 + t135 * t62) + pkin(9) * t196 + t152; -mrSges(6,1) - m(6) * pkin(4) + m(7) * (t131 * t97 + t135 * t96) - t157; m(6) + m(7) * (t131 ^ 2 + t135 ^ 2); t150; t184 + t208; t184 + t209; -Ifges(7,3) - t186 + t187; t157; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
