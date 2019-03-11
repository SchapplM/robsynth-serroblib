% Calculate joint inertia matrix for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:43
% EndTime: 2019-03-10 03:56:48
% DurationCPUTime: 2.40s
% Computational Cost: add. (6010->409), mult. (13217->599), div. (0->0), fcn. (15201->12), ass. (0->174)
t171 = sin(qJ(6));
t176 = cos(qJ(6));
t173 = sin(qJ(4));
t174 = sin(qJ(3));
t178 = cos(qJ(4));
t179 = cos(qJ(3));
t132 = -t173 * t174 + t178 * t179;
t156 = -pkin(3) * t179 - pkin(2);
t110 = -pkin(4) * t132 + t156;
t133 = t173 * t179 + t174 * t178;
t172 = sin(qJ(5));
t177 = cos(qJ(5));
t94 = -t132 * t177 + t133 * t172;
t95 = t132 * t172 + t133 * t177;
t50 = pkin(5) * t94 - pkin(12) * t95 + t110;
t241 = -pkin(10) - pkin(9);
t198 = t241 * t174;
t199 = t179 * t241;
t103 = t173 * t199 + t178 * t198;
t189 = -t133 * pkin(11) + t103;
t104 = t173 * t198 - t178 * t199;
t85 = pkin(11) * t132 + t104;
t56 = t172 * t189 + t177 * t85;
t26 = -t171 * t56 + t176 * t50;
t27 = t171 * t50 + t176 * t56;
t192 = -t171 * t26 + t176 * t27;
t170 = cos(pkin(6));
t169 = sin(pkin(6));
t175 = sin(qJ(2));
t208 = t169 * t175;
t122 = t170 * t179 - t174 * t208;
t123 = t170 * t174 + t179 * t208;
t86 = t122 * t178 - t123 * t173;
t87 = t122 * t173 + t123 * t178;
t57 = t172 * t87 - t177 * t86;
t58 = t172 * t86 + t177 * t87;
t148 = pkin(8) * t208;
t180 = cos(qJ(2));
t235 = pkin(1) * t180;
t111 = t148 + (-pkin(2) - t235) * t170;
t88 = -pkin(3) * t122 + t111;
t60 = -pkin(4) * t86 + t88;
t15 = pkin(5) * t57 - pkin(12) * t58 + t60;
t207 = t169 * t180;
t127 = pkin(1) * t170 * t175 + pkin(8) * t207;
t112 = pkin(9) * t170 + t127;
t113 = (-pkin(2) * t180 - pkin(9) * t175 - pkin(1)) * t169;
t76 = -t112 * t174 + t113 * t179;
t66 = -pkin(3) * t207 - pkin(10) * t123 + t76;
t77 = t112 * t179 + t113 * t174;
t71 = pkin(10) * t122 + t77;
t30 = -t173 * t71 + t178 * t66;
t18 = -pkin(4) * t207 - pkin(11) * t87 + t30;
t31 = t173 * t66 + t178 * t71;
t23 = pkin(11) * t86 + t31;
t9 = t172 * t18 + t177 * t23;
t6 = -pkin(12) * t207 + t9;
t2 = t15 * t176 - t171 * t6;
t3 = t15 * t171 + t176 * t6;
t194 = -t171 * t2 + t176 * t3;
t250 = -Ifges(5,5) * t87 - Ifges(5,6) * t86;
t249 = -Ifges(4,5) * t123 - Ifges(4,6) * t122;
t216 = t171 * t95;
t64 = -mrSges(7,2) * t94 - mrSges(7,3) * t216;
t211 = t176 * t95;
t65 = mrSges(7,1) * t94 - mrSges(7,3) * t211;
t248 = -t171 * t65 + t176 * t64;
t41 = -t171 * t58 - t176 * t207;
t19 = -mrSges(7,2) * t57 + mrSges(7,3) * t41;
t42 = -t171 * t207 + t176 * t58;
t20 = mrSges(7,1) * t57 - mrSges(7,3) * t42;
t247 = -t171 * t20 + t176 * t19;
t165 = t171 ^ 2;
t226 = mrSges(7,3) * t165;
t157 = pkin(12) * t226;
t167 = t176 ^ 2;
t225 = mrSges(7,3) * t167;
t158 = pkin(12) * t225;
t138 = -mrSges(7,1) * t176 + mrSges(7,2) * t171;
t232 = pkin(5) * t138;
t246 = t157 + t158 - t232;
t233 = pkin(4) * t177;
t154 = -pkin(5) - t233;
t125 = t154 * t138;
t153 = pkin(4) * t172 + pkin(12);
t136 = t153 * t226;
t137 = t153 * t225;
t159 = mrSges(6,1) * t233;
t245 = t125 + t136 + t137 + t159;
t54 = t172 * t85 - t177 * t189;
t244 = t54 ^ 2;
t243 = 0.2e1 * t54;
t242 = t41 / 0.2e1;
t140 = Ifges(7,5) * t171 + Ifges(7,6) * t176;
t240 = t140 / 0.2e1;
t223 = Ifges(7,4) * t176;
t143 = Ifges(7,1) * t171 + t223;
t239 = t143 / 0.2e1;
t238 = t171 / 0.2e1;
t237 = t176 / 0.2e1;
t234 = pkin(3) * t173;
t229 = Ifges(5,3) + Ifges(6,3);
t228 = -Ifges(6,5) * t58 + Ifges(6,6) * t57;
t227 = Ifges(6,5) * t95 - Ifges(6,6) * t94;
t224 = Ifges(7,4) * t171;
t155 = pkin(3) * t178 + pkin(4);
t121 = t155 * t172 + t177 * t234;
t222 = t121 * mrSges(6,2);
t126 = t170 * t235 - t148;
t221 = t126 * mrSges(3,1);
t220 = t127 * mrSges(3,2);
t215 = t172 * mrSges(6,2);
t210 = t153 * t171;
t209 = t153 * t176;
t206 = Ifges(5,5) * t133 + Ifges(5,6) * t132;
t205 = Ifges(4,5) * t174 + Ifges(4,6) * t179;
t204 = t165 + t167;
t203 = t174 ^ 2 + t179 ^ 2;
t202 = pkin(4) * t215;
t201 = Ifges(4,3) + t229;
t12 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t57;
t141 = Ifges(7,2) * t176 + t224;
t200 = t141 * t176 + t143 * t171 + Ifges(6,3);
t197 = Ifges(3,5) * t208 + Ifges(3,6) * t207 + Ifges(3,3) * t170;
t196 = t204 * t153;
t195 = Ifges(5,3) + t200;
t193 = mrSges(7,1) * t171 + mrSges(7,2) * t176;
t8 = -t172 * t23 + t177 * t18;
t36 = Ifges(7,5) * t211 - Ifges(7,6) * t216 + Ifges(7,3) * t94;
t120 = t155 * t177 - t172 * t234;
t191 = (mrSges(5,1) * t178 - mrSges(5,2) * t173) * pkin(3);
t119 = pkin(12) + t121;
t108 = t119 * t226;
t109 = t119 * t225;
t114 = t120 * mrSges(6,1);
t118 = -pkin(5) - t120;
t99 = t118 * t138;
t190 = t108 + t109 + t114 + t200 + t99 - t222;
t13 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t57;
t14 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t57;
t5 = pkin(5) * t207 - t8;
t188 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + mrSges(7,3) * t194 + t13 * t237 + t5 * t138 + t14 * t238 + t141 * t242 + t42 * t239 + t57 * t240 - t228;
t37 = Ifges(7,6) * t94 + (-Ifges(7,2) * t171 + t223) * t95;
t38 = Ifges(7,5) * t94 + (Ifges(7,1) * t176 - t224) * t95;
t187 = -t56 * mrSges(6,2) + t37 * t237 + t38 * t238 + t227 - t141 * t216 / 0.2e1 + t211 * t239 + t94 * t240 + (-mrSges(6,1) + t138) * t54 + t192 * mrSges(7,3);
t186 = mrSges(5,1) * t30 - t31 * mrSges(5,2) + t188 - t250;
t185 = mrSges(5,1) * t103 - t104 * mrSges(5,2) + t187 + t206;
t144 = Ifges(4,1) * t174 + Ifges(4,4) * t179;
t142 = Ifges(4,4) * t174 + Ifges(4,2) * t179;
t139 = -mrSges(4,1) * t179 + mrSges(4,2) * t174;
t102 = -mrSges(4,1) * t207 - mrSges(4,3) * t123;
t101 = mrSges(4,2) * t207 + mrSges(4,3) * t122;
t98 = Ifges(5,1) * t133 + Ifges(5,4) * t132;
t97 = Ifges(5,4) * t133 + Ifges(5,2) * t132;
t96 = -mrSges(5,1) * t132 + mrSges(5,2) * t133;
t90 = -mrSges(4,1) * t122 + mrSges(4,2) * t123;
t80 = Ifges(4,1) * t123 + Ifges(4,4) * t122 - Ifges(4,5) * t207;
t79 = Ifges(4,4) * t123 + Ifges(4,2) * t122 - Ifges(4,6) * t207;
t73 = -mrSges(5,1) * t207 - mrSges(5,3) * t87;
t72 = mrSges(5,2) * t207 + mrSges(5,3) * t86;
t69 = Ifges(6,1) * t95 - Ifges(6,4) * t94;
t68 = Ifges(6,4) * t95 - Ifges(6,2) * t94;
t67 = mrSges(6,1) * t94 + mrSges(6,2) * t95;
t61 = t193 * t95;
t59 = -mrSges(5,1) * t86 + mrSges(5,2) * t87;
t48 = Ifges(5,1) * t87 + Ifges(5,4) * t86 - Ifges(5,5) * t207;
t47 = Ifges(5,4) * t87 + Ifges(5,2) * t86 - Ifges(5,6) * t207;
t46 = -mrSges(6,1) * t207 - mrSges(6,3) * t58;
t45 = mrSges(6,2) * t207 - mrSges(6,3) * t57;
t28 = mrSges(6,1) * t57 + mrSges(6,2) * t58;
t25 = Ifges(6,1) * t58 - Ifges(6,4) * t57 - Ifges(6,5) * t207;
t24 = Ifges(6,4) * t58 - Ifges(6,2) * t57 - Ifges(6,6) * t207;
t16 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t1 = [m(3) * (pkin(1) ^ 2 * t169 ^ 2 + t126 ^ 2 + t127 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + m(6) * (t60 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t88 ^ 2) + m(4) * (t111 ^ 2 + t76 ^ 2 + t77 ^ 2) + (t12 - t24) * t57 + ((-0.2e1 * t126 * mrSges(3,3) + Ifges(3,5) * t170 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t175) * t169) * t175 + (0.2e1 * t127 * mrSges(3,3) + Ifges(3,6) * t170 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t175 + (Ifges(3,2) + t201) * t180) * t169 + t228 + t249 + t250) * t180) * t169 + Ifges(2,3) + (t197 - 0.2e1 * t220 + 0.2e1 * t221) * t170 + 0.2e1 * t5 * t16 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + t41 * t13 + t42 * t14 + 0.2e1 * t9 * t45 + 0.2e1 * t8 * t46 + t58 * t25 + 0.2e1 * t60 * t28 + 0.2e1 * t31 * t72 + 0.2e1 * t30 * t73 + t86 * t47 + t87 * t48 + 0.2e1 * t88 * t59 + 0.2e1 * t77 * t101 + 0.2e1 * t76 * t102 + 0.2e1 * t111 * t90 + t122 * t79 + t123 * t80; (t132 * t31 - t133 * t30) * mrSges(5,3) + (t16 - t46) * t54 + (t77 * mrSges(4,3) + pkin(9) * t101 + t79 / 0.2e1) * t179 + (-t76 * mrSges(4,3) - pkin(9) * t102 + t80 / 0.2e1) * t174 + m(7) * (t2 * t26 + t27 * t3 + t5 * t54) + m(6) * (t110 * t60 - t54 * t8 + t56 * t9) + m(5) * (t103 * t30 + t104 * t31 + t156 * t88) + (-t9 * mrSges(6,3) + t12 / 0.2e1 - t24 / 0.2e1) * t94 + (t36 / 0.2e1 - t68 / 0.2e1) * t57 - t220 + t221 - (t205 + t206 + t227) * t207 / 0.2e1 + m(4) * (-pkin(2) * t111 + (-t174 * t76 + t179 * t77) * pkin(9)) + (-t8 * mrSges(6,3) - t171 * t13 / 0.2e1 + t14 * t237 + t25 / 0.2e1) * t95 + t37 * t242 + t197 + t26 * t20 + t27 * t19 + t42 * t38 / 0.2e1 + t56 * t45 + t5 * t61 + t3 * t64 + t2 * t65 + t60 * t67 + t58 * t69 / 0.2e1 - pkin(2) * t90 + t88 * t96 + t86 * t97 / 0.2e1 + t87 * t98 / 0.2e1 + t103 * t73 + t104 * t72 + t110 * t28 + t132 * t47 / 0.2e1 + t133 * t48 / 0.2e1 + t111 * t139 + t122 * t142 / 0.2e1 + t123 * t144 / 0.2e1 + t156 * t59; -0.2e1 * pkin(2) * t139 + 0.2e1 * t110 * t67 + t132 * t97 + t133 * t98 + t179 * t142 + t174 * t144 + 0.2e1 * t156 * t96 + 0.2e1 * t26 * t65 + 0.2e1 * t27 * t64 + t61 * t243 + Ifges(3,3) + (-0.2e1 * mrSges(6,3) * t56 + t36 - t68) * t94 + (mrSges(6,3) * t243 - t171 * t37 + t176 * t38 + t69) * t95 + m(7) * (t26 ^ 2 + t27 ^ 2 + t244) + m(6) * (t110 ^ 2 + t56 ^ 2 + t244) + m(5) * (t103 ^ 2 + t104 ^ 2 + t156 ^ 2) + m(4) * (pkin(9) ^ 2 * t203 + pkin(2) ^ 2) + 0.2e1 * (-t103 * t133 + t104 * t132) * mrSges(5,3) + 0.2e1 * t203 * pkin(9) * mrSges(4,3); t186 + t247 * t119 + m(7) * (t118 * t5 + t119 * t194) + m(6) * (t120 * t8 + t121 * t9) + (t173 * t72 + t178 * t73 + m(5) * (t173 * t31 + t178 * t30)) * pkin(3) + t76 * mrSges(4,1) - t77 * mrSges(4,2) + t118 * t16 + t120 * t46 + t121 * t45 - t201 * t207 - t249; (-t120 * t95 - t121 * t94) * mrSges(6,3) + t248 * t119 + t185 + m(6) * (-t120 * t54 + t121 * t56) + (-mrSges(4,1) * t174 - mrSges(4,2) * t179) * pkin(9) + m(7) * (t118 * t54 + t119 * t192) + (m(5) * (t103 * t178 + t104 * t173) + (t132 * t173 - t133 * t178) * mrSges(5,3)) * pkin(3) + t118 * t61 + t205; -0.2e1 * t222 + Ifges(4,3) + 0.2e1 * t108 + 0.2e1 * t109 + 0.2e1 * t114 + 0.2e1 * t99 + 0.2e1 * t191 + m(7) * (t119 ^ 2 * t204 + t118 ^ 2) + m(6) * (t120 ^ 2 + t121 ^ 2) + m(5) * (t173 ^ 2 + t178 ^ 2) * pkin(3) ^ 2 + t195; t186 - t229 * t207 + m(7) * (t153 * t194 + t154 * t5) + (t172 * t45 + t177 * t46 + m(6) * (t172 * t9 + t177 * t8)) * pkin(4) - t20 * t210 + t19 * t209 + t154 * t16; t185 + m(7) * (t153 * t192 + t154 * t54) + (m(6) * (t172 * t56 - t177 * t54) + (-t172 * t94 - t177 * t95) * mrSges(6,3)) * pkin(4) - t65 * t210 + t64 * t209 + t154 * t61; t190 + t191 + m(7) * (t118 * t154 + t119 * t196) + (m(6) * (t120 * t177 + t121 * t172) - t215) * pkin(4) + Ifges(5,3) + t245; -0.2e1 * t202 + 0.2e1 * t125 + 0.2e1 * t136 + 0.2e1 * t137 + 0.2e1 * t159 + m(7) * (t153 ^ 2 * t204 + t154 ^ 2) + m(6) * (t172 ^ 2 + t177 ^ 2) * pkin(4) ^ 2 + t195; -Ifges(6,3) * t207 + t188 + (-m(7) * t5 - t16) * pkin(5) + (m(7) * t194 + t247) * pkin(12); t187 + (-m(7) * t54 - t61) * pkin(5) + (m(7) * t192 + t248) * pkin(12); m(7) * (pkin(12) * t119 * t204 - pkin(5) * t118) + t190 + t246; -t202 + m(7) * (-pkin(5) * t154 + pkin(12) * t196) + t200 + t245 + t246; 0.2e1 * t158 + 0.2e1 * t157 - 0.2e1 * t232 + m(7) * (pkin(12) ^ 2 * t204 + pkin(5) ^ 2) + t200; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t26 - mrSges(7,2) * t27 + t36; -t119 * t193 + t140; -t153 * t193 + t140; -pkin(12) * t193 + t140; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
