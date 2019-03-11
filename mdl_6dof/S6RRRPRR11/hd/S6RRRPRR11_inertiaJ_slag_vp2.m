% Calculate joint inertia matrix for
% S6RRRPRR11
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:54
% EndTime: 2019-03-09 19:24:57
% DurationCPUTime: 1.95s
% Computational Cost: add. (2861->400), mult. (6196->551), div. (0->0), fcn. (6561->10), ass. (0->146)
t139 = sin(qJ(3));
t143 = cos(qJ(3));
t199 = t139 ^ 2 + t143 ^ 2;
t137 = sin(qJ(6));
t141 = cos(qJ(6));
t95 = -t141 * mrSges(7,1) + mrSges(7,2) * t137;
t177 = -t95 + mrSges(6,1);
t164 = t137 ^ 2 + t141 ^ 2;
t159 = t164 * mrSges(7,3);
t138 = sin(qJ(5));
t142 = cos(qJ(5));
t87 = -t138 * t143 + t139 * t142;
t172 = t137 * t87;
t86 = t138 * t139 + t142 * t143;
t51 = -mrSges(7,2) * t86 - mrSges(7,3) * t172;
t169 = t141 * t87;
t52 = mrSges(7,1) * t86 - mrSges(7,3) * t169;
t198 = -t137 * t52 + t141 * t51;
t135 = sin(pkin(6));
t144 = cos(qJ(2));
t166 = t135 * t144;
t136 = cos(pkin(6));
t140 = sin(qJ(2));
t167 = t135 * t140;
t78 = -t136 * t143 + t139 * t167;
t79 = t136 * t139 + t143 * t167;
t47 = t138 * t78 + t142 * t79;
t29 = -t137 * t47 + t141 * t166;
t46 = t138 * t79 - t78 * t142;
t12 = -mrSges(7,2) * t46 + mrSges(7,3) * t29;
t30 = t137 * t166 + t141 * t47;
t13 = mrSges(7,1) * t46 - mrSges(7,3) * t30;
t197 = t141 * t12 - t137 * t13;
t196 = (-Ifges(5,4) - Ifges(4,5)) * t79 + (Ifges(4,6) - Ifges(5,6)) * t78;
t195 = m(7) * pkin(5) + t177;
t188 = pkin(9) - pkin(10);
t105 = t188 * t143;
t161 = t188 * t139;
t60 = t105 * t138 - t142 * t161;
t194 = t60 ^ 2;
t193 = 0.2e1 * t60;
t192 = -0.2e1 * t95;
t191 = t29 / 0.2e1;
t190 = t30 / 0.2e1;
t98 = Ifges(7,5) * t137 + Ifges(7,6) * t141;
t189 = t98 / 0.2e1;
t187 = -t137 / 0.2e1;
t186 = t137 / 0.2e1;
t185 = t141 / 0.2e1;
t184 = pkin(1) * t136;
t80 = -pkin(8) * t167 + t144 * t184;
t183 = t80 * mrSges(3,1);
t81 = pkin(8) * t166 + t140 * t184;
t182 = t81 * mrSges(3,2);
t145 = -pkin(3) - pkin(4);
t92 = -qJ(4) * t138 + t142 * t145;
t181 = t92 * mrSges(6,1);
t93 = t142 * qJ(4) + t138 * t145;
t180 = t93 * mrSges(6,2);
t179 = Ifges(5,2) + Ifges(4,3);
t11 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t35 = mrSges(6,1) * t166 - mrSges(6,3) * t47;
t178 = t11 - t35;
t67 = pkin(9) * t136 + t81;
t68 = (-pkin(2) * t144 - pkin(9) * t140 - pkin(1)) * t135;
t36 = -t139 * t67 + t143 * t68;
t28 = pkin(3) * t166 - t36;
t18 = pkin(4) * t166 - pkin(10) * t79 + t28;
t37 = t139 * t68 + t143 * t67;
t25 = -qJ(4) * t166 + t37;
t23 = pkin(10) * t78 + t25;
t6 = t138 * t18 + t142 * t23;
t176 = Ifges(7,4) * t137;
t175 = Ifges(7,4) * t141;
t168 = t142 * t60;
t59 = mrSges(5,1) * t166 + t79 * mrSges(5,2);
t165 = t199 * pkin(9) ^ 2;
t8 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t46;
t94 = -t143 * pkin(3) - t139 * qJ(4) - pkin(2);
t163 = Ifges(6,5) * t47 - Ifges(6,6) * t46 + Ifges(6,3) * t166;
t162 = Ifges(5,6) * t166;
t160 = Ifges(3,5) * t167 + Ifges(3,6) * t166 + Ifges(3,3) * t136;
t66 = -t136 * pkin(2) - t80;
t91 = -pkin(11) + t93;
t158 = t164 * t91;
t157 = t164 * t138;
t82 = t143 * pkin(4) - t94;
t62 = t142 * t105 + t138 * t161;
t84 = Ifges(6,6) * t86;
t85 = Ifges(6,5) * t87;
t156 = -t62 * mrSges(6,2) - t84 + t85;
t4 = pkin(11) * t166 + t6;
t24 = t78 * pkin(3) - t79 * qJ(4) + t66;
t22 = -pkin(4) * t78 - t24;
t7 = pkin(5) * t46 - pkin(11) * t47 + t22;
t1 = -t137 * t4 + t141 * t7;
t2 = t137 * t7 + t141 * t4;
t154 = -t1 * t137 + t141 * t2;
t153 = mrSges(7,1) * t137 + mrSges(7,2) * t141;
t152 = -pkin(3) * t139 + qJ(4) * t143;
t42 = pkin(5) * t86 - pkin(11) * t87 + t82;
t19 = -t137 * t62 + t141 * t42;
t20 = t137 * t42 + t141 * t62;
t151 = -t137 * t19 + t141 * t20;
t5 = -t138 * t23 + t142 * t18;
t31 = Ifges(7,5) * t169 - Ifges(7,6) * t172 + Ifges(7,3) * t86;
t100 = Ifges(7,2) * t141 + t176;
t102 = Ifges(7,1) * t137 + t175;
t150 = t141 * t100 + t137 * t102 + Ifges(6,3);
t149 = t100 * t187 + t102 * t185;
t3 = -pkin(5) * t166 - t5;
t148 = t5 * mrSges(6,1) - t6 * mrSges(6,2) + t100 * t191 + t102 * t190 + t46 * t189 + t3 * t95 + t163;
t133 = t142 ^ 2;
t130 = t138 ^ 2;
t124 = Ifges(5,4) * t139;
t123 = Ifges(4,5) * t139;
t122 = Ifges(4,6) * t143;
t104 = Ifges(4,1) * t139 + Ifges(4,4) * t143;
t103 = Ifges(5,1) * t139 - Ifges(5,5) * t143;
t101 = Ifges(4,4) * t139 + Ifges(4,2) * t143;
t99 = Ifges(5,5) * t139 - Ifges(5,3) * t143;
t97 = -mrSges(4,1) * t143 + mrSges(4,2) * t139;
t96 = -mrSges(5,1) * t143 - mrSges(5,3) * t139;
t90 = pkin(5) - t92;
t58 = -mrSges(4,1) * t166 - mrSges(4,3) * t79;
t57 = mrSges(4,2) * t166 - mrSges(4,3) * t78;
t56 = -mrSges(5,2) * t78 - mrSges(5,3) * t166;
t55 = Ifges(6,1) * t87 - Ifges(6,4) * t86;
t54 = Ifges(6,4) * t87 - Ifges(6,2) * t86;
t53 = mrSges(6,1) * t86 + mrSges(6,2) * t87;
t50 = t153 * t87;
t49 = mrSges(4,1) * t78 + mrSges(4,2) * t79;
t48 = mrSges(5,1) * t78 - mrSges(5,3) * t79;
t41 = Ifges(4,1) * t79 - Ifges(4,4) * t78 - Ifges(4,5) * t166;
t40 = Ifges(5,1) * t79 - Ifges(5,4) * t166 + Ifges(5,5) * t78;
t39 = Ifges(4,4) * t79 - Ifges(4,2) * t78 - Ifges(4,6) * t166;
t38 = Ifges(5,5) * t79 + Ifges(5,3) * t78 - t162;
t34 = -mrSges(6,2) * t166 - mrSges(6,3) * t46;
t33 = Ifges(7,5) * t86 + (Ifges(7,1) * t141 - t176) * t87;
t32 = Ifges(7,6) * t86 + (-Ifges(7,2) * t137 + t175) * t87;
t16 = mrSges(6,1) * t46 + mrSges(6,2) * t47;
t15 = Ifges(6,1) * t47 - Ifges(6,4) * t46 + Ifges(6,5) * t166;
t14 = Ifges(6,4) * t47 - Ifges(6,2) * t46 + Ifges(6,6) * t166;
t10 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t46;
t9 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t46;
t17 = [((-0.2e1 * t80 * mrSges(3,3) + Ifges(3,5) * t136 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t140) * t135) * t140 + (0.2e1 * t81 * mrSges(3,3) + Ifges(3,6) * t136 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t140 + (Ifges(3,2) + t179) * t144) * t135 + t163 + t196) * t144) * t135 + m(3) * (pkin(1) ^ 2 * t135 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2 + t66 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t28 ^ 2) + m(6) * (t22 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (t40 + t41) * t79 + (-t39 + t38) * t78 + (t8 - t14) * t46 + 0.2e1 * t66 * t49 + 0.2e1 * t25 * t56 + 0.2e1 * t37 * t57 + 0.2e1 * t36 * t58 + 0.2e1 * t28 * t59 + t47 * t15 + 0.2e1 * t24 * t48 + 0.2e1 * t6 * t34 + 0.2e1 * t5 * t35 + t29 * t9 + t30 * t10 + 0.2e1 * t2 * t12 + 0.2e1 * t1 * t13 + 0.2e1 * t22 * t16 + 0.2e1 * t3 * t11 + Ifges(2,3) + (t160 - 0.2e1 * t182 + 0.2e1 * t183) * t136; t24 * t96 + t66 * t97 + t94 * t48 + t82 * t16 + t160 + (t40 / 0.2e1 + t41 / 0.2e1 + t28 * mrSges(5,2) - t36 * mrSges(4,3) + (-t58 + t59) * pkin(9)) * t139 + m(6) * (t22 * t82 - t5 * t60 + t6 * t62) + m(7) * (t1 * t19 + t2 * t20 + t3 * t60) + (t8 / 0.2e1 - t14 / 0.2e1 - t6 * mrSges(6,3)) * t86 + (t103 / 0.2e1 + t104 / 0.2e1) * t79 + (t99 / 0.2e1 - t101 / 0.2e1) * t78 + (-t54 / 0.2e1 + t31 / 0.2e1) * t46 + t62 * t34 - pkin(2) * t49 + t3 * t50 + t2 * t51 + t1 * t52 + t22 * t53 + t47 * t55 / 0.2e1 + t19 * t13 + t20 * t12 - t182 + t183 + (-t124 / 0.2e1 - t123 / 0.2e1 - t122 / 0.2e1 + t85 / 0.2e1 - t84 / 0.2e1) * t166 + (-t38 / 0.2e1 + t39 / 0.2e1 + t25 * mrSges(5,2) + t37 * mrSges(4,3) + t162 / 0.2e1 + (t56 + t57) * pkin(9)) * t143 + t33 * t190 + t32 * t191 + t178 * t60 + (t15 / 0.2e1 + t10 * t185 + t9 * t187 - t5 * mrSges(6,3)) * t87 + m(4) * (-pkin(2) * t66 + (-t139 * t36 + t143 * t37) * pkin(9)) + m(5) * (t24 * t94 + (t139 * t28 + t143 * t25) * pkin(9)); -0.2e1 * pkin(2) * t97 + 0.2e1 * t19 * t52 + 0.2e1 * t20 * t51 + t50 * t193 + 0.2e1 * t82 * t53 + 0.2e1 * t94 * t96 + Ifges(3,3) + (-t99 + t101) * t143 + (t103 + t104) * t139 + (-0.2e1 * mrSges(6,3) * t62 + t31 - t54) * t86 + (mrSges(6,3) * t193 - t137 * t32 + t141 * t33 + t55) * t87 + m(5) * (t94 ^ 2 + t165) + m(4) * (pkin(2) ^ 2 + t165) + m(6) * (t62 ^ 2 + t82 ^ 2 + t194) + m(7) * (t19 ^ 2 + t20 ^ 2 + t194) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(9) * t199; -t196 + t90 * t11 + t92 * t35 + t93 * t34 + qJ(4) * t56 - pkin(3) * t59 + t36 * mrSges(4,1) - t37 * mrSges(4,2) + t25 * mrSges(5,3) - t28 * mrSges(5,1) + (-t9 / 0.2e1 - t2 * mrSges(7,3) + t91 * t12) * t141 + (-t10 / 0.2e1 - t91 * t13 + t1 * mrSges(7,3)) * t137 + m(6) * (t5 * t92 + t6 * t93) + m(5) * (-pkin(3) * t28 + qJ(4) * t25) + m(7) * (t154 * t91 + t3 * t90) - t179 * t166 - t148; t90 * t50 + t122 + t123 + t124 + (-t98 / 0.2e1 - t93 * mrSges(6,3)) * t86 + t177 * t60 + t152 * mrSges(5,2) + (-t32 / 0.2e1 - t20 * mrSges(7,3) + t91 * t51) * t141 + (-t33 / 0.2e1 - t91 * t52 + t19 * mrSges(7,3)) * t137 + m(7) * (t151 * t91 + t60 * t90) + m(6) * (-t60 * t92 + t62 * t93) + (-t92 * mrSges(6,3) - t149) * t87 + (m(5) * t152 + (-mrSges(4,1) - mrSges(5,1)) * t139) * pkin(9) - t156 + (-Ifges(5,6) + (-mrSges(4,2) + mrSges(5,3)) * pkin(9)) * t143; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t181 + 0.2e1 * t180 + 0.2e1 * qJ(4) * mrSges(5,3) + t90 * t192 - 0.2e1 * t91 * t159 + m(7) * (t164 * t91 ^ 2 + t90 ^ 2) + m(6) * (t92 ^ 2 + t93 ^ 2) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t150 + t179; -t178 * t142 + (t34 + t197) * t138 + m(7) * (t154 * t138 - t142 * t3) + m(6) * (t138 * t6 + t142 * t5) + m(5) * t28 + t59; (-t87 * mrSges(6,3) - t50) * t142 + (m(5) * pkin(9) + mrSges(5,2)) * t139 + (-t86 * mrSges(6,3) + t198) * t138 + m(7) * (t151 * t138 - t168) + m(6) * (t138 * t62 - t168); -m(5) * pkin(3) - mrSges(5,1) - t177 * t142 + (mrSges(6,2) - t159) * t138 + m(7) * (-t142 * t90 + t91 * t157) + m(6) * (t138 * t93 + t142 * t92); m(5) + m(6) * (t130 + t133) + m(7) * (t164 * t130 + t133); t154 * mrSges(7,3) + t10 * t186 + t9 * t185 + t148 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t154 + t197) * pkin(11); t151 * mrSges(7,3) - pkin(5) * t50 + t149 * t87 + t32 * t185 + t33 * t186 + t86 * t189 + t156 - t195 * t60 + (m(7) * t151 + t198) * pkin(11); m(7) * (-pkin(5) * t90 + pkin(11) * t158) - t180 + t181 + (pkin(5) + t90) * t95 + (-t164 * pkin(11) + t158) * mrSges(7,3) - t150; -t138 * mrSges(6,2) + (m(7) * pkin(11) + mrSges(7,3)) * t157 + t195 * t142; m(7) * (t164 * pkin(11) ^ 2 + pkin(5) ^ 2) + pkin(5) * t192 + 0.2e1 * pkin(11) * t159 + t150; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t19 - mrSges(7,2) * t20 + t31; (-mrSges(7,2) * t91 - Ifges(7,6)) * t141 + (-mrSges(7,1) * t91 - Ifges(7,5)) * t137; -t153 * t138; -t153 * pkin(11) + t98; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
