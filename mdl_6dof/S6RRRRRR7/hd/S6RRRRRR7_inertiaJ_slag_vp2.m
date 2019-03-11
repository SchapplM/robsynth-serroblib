% Calculate joint inertia matrix for
% S6RRRRRR7
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:29:05
% EndTime: 2019-03-10 04:29:10
% DurationCPUTime: 2.42s
% Computational Cost: add. (5919->493), mult. (13299->711), div. (0->0), fcn. (14933->12), ass. (0->194)
t253 = 2 * pkin(9);
t188 = sin(qJ(6));
t193 = cos(qJ(6));
t187 = cos(pkin(6));
t191 = sin(qJ(3));
t196 = cos(qJ(3));
t186 = sin(pkin(6));
t192 = sin(qJ(2));
t223 = t186 * t192;
t136 = t187 * t191 + t196 * t223;
t190 = sin(qJ(4));
t195 = cos(qJ(4));
t197 = cos(qJ(2));
t222 = t186 * t197;
t101 = -t136 * t190 - t195 * t222;
t102 = t136 * t195 - t190 * t222;
t189 = sin(qJ(5));
t194 = cos(qJ(5));
t58 = t101 * t194 - t102 * t189;
t59 = t101 * t189 + t102 * t194;
t29 = -t188 * t59 + t193 * t58;
t252 = t29 / 0.2e1;
t30 = t188 * t58 + t193 * t59;
t251 = t30 / 0.2e1;
t135 = -t187 * t196 + t191 * t223;
t48 = Ifges(5,1) * t102 + Ifges(5,4) * t101 + Ifges(5,5) * t135;
t250 = t48 / 0.2e1;
t249 = t58 / 0.2e1;
t248 = t59 / 0.2e1;
t145 = t189 * t195 + t190 * t194;
t126 = t145 * t191;
t144 = -t189 * t190 + t194 * t195;
t127 = t144 * t191;
t80 = -t126 * t193 - t127 * t188;
t247 = t80 / 0.2e1;
t81 = -t126 * t188 + t127 * t193;
t246 = t81 / 0.2e1;
t90 = t144 * t193 - t145 * t188;
t245 = t90 / 0.2e1;
t91 = t144 * t188 + t145 * t193;
t244 = t91 / 0.2e1;
t243 = -pkin(11) - pkin(10);
t229 = Ifges(5,4) * t195;
t123 = -Ifges(5,6) * t196 + (-Ifges(5,2) * t190 + t229) * t191;
t242 = t123 / 0.2e1;
t230 = Ifges(5,4) * t190;
t124 = -Ifges(5,5) * t196 + (Ifges(5,1) * t195 - t230) * t191;
t241 = t124 / 0.2e1;
t240 = -t126 / 0.2e1;
t239 = t127 / 0.2e1;
t238 = t144 / 0.2e1;
t237 = t145 / 0.2e1;
t158 = Ifges(5,1) * t190 + t229;
t236 = t158 / 0.2e1;
t235 = pkin(1) * t197;
t234 = pkin(4) * t189;
t233 = pkin(9) * t196;
t181 = t191 * pkin(9);
t232 = -Ifges(7,3) - Ifges(6,3);
t164 = pkin(8) * t223;
t117 = t164 + (-pkin(2) - t235) * t187;
t64 = pkin(3) * t135 - pkin(10) * t136 + t117;
t138 = t187 * t192 * pkin(1) + pkin(8) * t222;
t118 = pkin(9) * t187 + t138;
t119 = (-pkin(2) * t197 - pkin(9) * t192 - pkin(1)) * t186;
t72 = t196 * t118 + t191 * t119;
t66 = -pkin(10) * t222 + t72;
t32 = -t190 * t66 + t195 * t64;
t18 = pkin(4) * t135 - pkin(11) * t102 + t32;
t33 = t190 * t64 + t195 * t66;
t26 = pkin(11) * t101 + t33;
t8 = t189 * t18 + t194 * t26;
t231 = Ifges(7,5) * t81 + Ifges(7,6) * t80;
t50 = Ifges(7,5) * t91 + Ifges(7,6) * t90;
t152 = -pkin(3) * t196 - pkin(10) * t191 - pkin(2);
t115 = t190 * t152 + t195 * t233;
t221 = t190 * t191;
t103 = -pkin(11) * t221 + t115;
t143 = t195 * t152;
t220 = t191 * t195;
t92 = -pkin(11) * t220 + t143 + (-pkin(9) * t190 - pkin(4)) * t196;
t55 = t194 * t103 + t189 * t92;
t228 = Ifges(7,3) * t196;
t172 = pkin(4) * t194 + pkin(5);
t134 = t172 * t188 + t193 * t234;
t227 = t134 * mrSges(7,2);
t137 = t187 * t235 - t164;
t226 = t137 * mrSges(3,1);
t225 = t138 * mrSges(3,2);
t224 = t188 * mrSges(7,2);
t219 = Ifges(6,5) * t127 - Ifges(6,6) * t126;
t218 = -Ifges(4,5) * t136 + Ifges(4,6) * t135;
t98 = Ifges(6,5) * t145 + Ifges(6,6) * t144;
t160 = t243 * t190;
t161 = t243 * t195;
t108 = t189 * t160 - t194 * t161;
t155 = Ifges(5,5) * t190 + Ifges(5,6) * t195;
t217 = Ifges(4,5) * t191 + Ifges(4,6) * t196;
t151 = pkin(4) * t221 + t181;
t216 = t190 ^ 2 + t195 ^ 2;
t215 = pkin(5) * t224;
t214 = -Ifges(5,3) + t232;
t9 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t135;
t21 = Ifges(6,5) * t59 + Ifges(6,6) * t58 + Ifges(6,3) * t135;
t46 = Ifges(5,5) * t102 + Ifges(5,6) * t101 + Ifges(5,3) * t135;
t213 = Ifges(3,5) * t223 + Ifges(3,6) * t222 + Ifges(3,3) * t187;
t173 = -pkin(4) * t195 - pkin(3);
t7 = t194 * t18 - t189 * t26;
t54 = -t103 * t189 + t194 * t92;
t71 = -t191 * t118 + t119 * t196;
t107 = t194 * t160 + t161 * t189;
t133 = t172 * t193 - t188 * t234;
t125 = t133 * mrSges(7,1);
t212 = Ifges(7,3) + t125 - t227;
t211 = t50 / 0.2e1 + t98 / 0.2e1 + t155 / 0.2e1;
t210 = Ifges(5,5) * t220 - Ifges(5,6) * t221;
t65 = pkin(3) * t222 - t71;
t209 = mrSges(5,1) * t190 + mrSges(5,2) * t195;
t37 = -pkin(5) * t196 - pkin(12) * t127 + t54;
t42 = -pkin(12) * t126 + t55;
t14 = -t188 * t42 + t193 * t37;
t15 = t188 * t37 + t193 * t42;
t208 = t14 * mrSges(7,1) - t15 * mrSges(7,2) + t231;
t83 = -pkin(12) * t145 + t107;
t84 = pkin(12) * t144 + t108;
t40 = -t188 * t84 + t193 * t83;
t41 = t188 * t83 + t193 * t84;
t207 = t40 * mrSges(7,1) - t41 * mrSges(7,2) + t50;
t4 = pkin(5) * t135 - pkin(12) * t59 + t7;
t5 = pkin(12) * t58 + t8;
t2 = -t188 * t5 + t193 * t4;
t3 = t188 * t4 + t193 * t5;
t206 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t9;
t205 = (mrSges(6,1) * t194 - mrSges(6,2) * t189) * pkin(4);
t45 = -pkin(4) * t101 + t65;
t204 = t54 * mrSges(6,1) - t55 * mrSges(6,2) + t208 + t219;
t203 = t107 * mrSges(6,1) - t108 * mrSges(6,2) + t207 + t98;
t202 = t7 * mrSges(6,1) - t8 * mrSges(6,2) + t206 + t21;
t199 = pkin(9) ^ 2;
t185 = t196 ^ 2;
t183 = t191 ^ 2;
t180 = t183 * t199;
t174 = t193 * pkin(5) * mrSges(7,1);
t159 = Ifges(4,1) * t191 + Ifges(4,4) * t196;
t157 = Ifges(4,4) * t191 + Ifges(4,2) * t196;
t156 = Ifges(5,2) * t195 + t230;
t154 = -mrSges(4,1) * t196 + mrSges(4,2) * t191;
t153 = -mrSges(5,1) * t195 + mrSges(5,2) * t190;
t150 = -mrSges(5,1) * t196 - mrSges(5,3) * t220;
t149 = mrSges(5,2) * t196 - mrSges(5,3) * t221;
t139 = t209 * t191;
t122 = -Ifges(5,3) * t196 + t210;
t116 = -pkin(5) * t144 + t173;
t114 = -t190 * t233 + t143;
t110 = -mrSges(6,1) * t196 - mrSges(6,3) * t127;
t109 = mrSges(6,2) * t196 - mrSges(6,3) * t126;
t106 = -mrSges(4,1) * t222 - mrSges(4,3) * t136;
t105 = mrSges(4,2) * t222 - mrSges(4,3) * t135;
t100 = Ifges(6,1) * t145 + Ifges(6,4) * t144;
t99 = Ifges(6,4) * t145 + Ifges(6,2) * t144;
t97 = -mrSges(6,1) * t144 + mrSges(6,2) * t145;
t93 = pkin(5) * t126 + t151;
t85 = mrSges(4,1) * t135 + mrSges(4,2) * t136;
t82 = mrSges(6,1) * t126 + mrSges(6,2) * t127;
t79 = Ifges(4,1) * t136 - Ifges(4,4) * t135 - Ifges(4,5) * t222;
t78 = Ifges(4,4) * t136 - Ifges(4,2) * t135 - Ifges(4,6) * t222;
t77 = Ifges(6,1) * t127 - Ifges(6,4) * t126 - Ifges(6,5) * t196;
t76 = Ifges(6,4) * t127 - Ifges(6,2) * t126 - Ifges(6,6) * t196;
t75 = -Ifges(6,3) * t196 + t219;
t70 = mrSges(5,1) * t135 - mrSges(5,3) * t102;
t69 = -mrSges(5,2) * t135 + mrSges(5,3) * t101;
t68 = -mrSges(7,1) * t196 - mrSges(7,3) * t81;
t67 = mrSges(7,2) * t196 + mrSges(7,3) * t80;
t60 = -mrSges(5,1) * t101 + mrSges(5,2) * t102;
t52 = Ifges(7,1) * t91 + Ifges(7,4) * t90;
t51 = Ifges(7,4) * t91 + Ifges(7,2) * t90;
t49 = -mrSges(7,1) * t90 + mrSges(7,2) * t91;
t47 = Ifges(5,4) * t102 + Ifges(5,2) * t101 + Ifges(5,6) * t135;
t44 = mrSges(6,1) * t135 - mrSges(6,3) * t59;
t43 = -mrSges(6,2) * t135 + mrSges(6,3) * t58;
t38 = -mrSges(7,1) * t80 + mrSges(7,2) * t81;
t36 = Ifges(7,1) * t81 + Ifges(7,4) * t80 - Ifges(7,5) * t196;
t35 = Ifges(7,4) * t81 + Ifges(7,2) * t80 - Ifges(7,6) * t196;
t34 = -t228 + t231;
t31 = -mrSges(6,1) * t58 + mrSges(6,2) * t59;
t24 = -pkin(5) * t58 + t45;
t23 = Ifges(6,1) * t59 + Ifges(6,4) * t58 + Ifges(6,5) * t135;
t22 = Ifges(6,4) * t59 + Ifges(6,2) * t58 + Ifges(6,6) * t135;
t20 = mrSges(7,1) * t135 - mrSges(7,3) * t30;
t19 = -mrSges(7,2) * t135 + mrSges(7,3) * t29;
t12 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t11 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t135;
t10 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t135;
t1 = [(t213 - 0.2e1 * t225 + 0.2e1 * t226) * t187 + ((-0.2e1 * t137 * mrSges(3,3) + Ifges(3,5) * t187 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t192) * t186) * t192 + (0.2e1 * t138 * mrSges(3,3) + Ifges(3,6) * t187 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t192 + (Ifges(4,3) + Ifges(3,2)) * t197) * t186 + t218) * t197) * t186 + (t9 + t21 + t46 - t78) * t135 + Ifges(2,3) + m(7) * (t2 ^ 2 + t24 ^ 2 + t3 ^ 2) + m(6) * (t45 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2 + t65 ^ 2) + m(4) * (t117 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(3) * (pkin(1) ^ 2 * t186 ^ 2 + t137 ^ 2 + t138 ^ 2) + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t24 * t12 + t29 * t10 + t30 * t11 + 0.2e1 * t8 * t43 + 0.2e1 * t7 * t44 + 0.2e1 * t45 * t31 + t58 * t22 + t59 * t23 + 0.2e1 * t65 * t60 + 0.2e1 * t33 * t69 + 0.2e1 * t32 * t70 + t101 * t47 + t102 * t48 + 0.2e1 * t72 * t105 + 0.2e1 * t71 * t106 + 0.2e1 * t117 * t85 + t136 * t79; m(4) * (-pkin(2) * t117 + (-t191 * t71 + t196 * t72) * pkin(9)) + t23 * t239 + t22 * t240 + t102 * t241 + t101 * t242 + t11 * t246 + t10 * t247 + t77 * t248 + t76 * t249 + t36 * t251 + t35 * t252 - t225 + t226 + t213 + (-t71 * mrSges(4,3) - t190 * t47 / 0.2e1 + t195 * t250 + t79 / 0.2e1 + (t60 - t106) * pkin(9)) * t191 + m(5) * (t114 * t32 + t115 * t33 + t181 * t65) + (t72 * mrSges(4,3) + pkin(9) * t105 - t9 / 0.2e1 - t21 / 0.2e1 - t46 / 0.2e1 + t78 / 0.2e1) * t196 + m(7) * (t14 * t2 + t15 * t3 + t24 * t93) + m(6) * (t151 * t45 + t54 * t7 + t55 * t8) + (t34 / 0.2e1 + t75 / 0.2e1 + t122 / 0.2e1 - t157 / 0.2e1) * t135 - t217 * t222 / 0.2e1 + t15 * t19 + t14 * t20 + t24 * t38 + t54 * t44 + t55 * t43 + t3 * t67 + t2 * t68 + t45 * t82 - pkin(2) * t85 + t93 * t12 + t8 * t109 + t7 * t110 + t114 * t70 + t115 * t69 + t65 * t139 + t33 * t149 + t32 * t150 + t151 * t31 + t117 * t154 + t136 * t159 / 0.2e1; -0.2e1 * pkin(2) * t154 + 0.2e1 * t55 * t109 + 0.2e1 * t54 * t110 + 0.2e1 * t114 * t150 + 0.2e1 * t115 * t149 - t126 * t76 + t127 * t77 + 0.2e1 * t14 * t68 + 0.2e1 * t15 * t67 + 0.2e1 * t151 * t82 + t80 * t35 + t81 * t36 + 0.2e1 * t93 * t38 + Ifges(3,3) + (t183 + t185) * mrSges(4,3) * t253 + (-t34 - t75 - t122 + t157) * t196 + (-t123 * t190 + t124 * t195 + t139 * t253 + t159) * t191 + m(7) * (t14 ^ 2 + t15 ^ 2 + t93 ^ 2) + m(6) * (t151 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t114 ^ 2 + t115 ^ 2 + t180) + m(4) * (pkin(2) ^ 2 + t185 * t199 + t180); m(5) * (-pkin(3) * t65 + (-t190 * t32 + t195 * t33) * pkin(10)) + t11 * t244 + t10 * t245 + t100 * t248 + t99 * t249 + t52 * t251 + t51 * t252 + t102 * t236 + t23 * t237 + t22 * t238 + (t33 * mrSges(5,3) + pkin(10) * t69 + t47 / 0.2e1) * t195 + m(7) * (t116 * t24 + t2 * t40 + t3 * t41) + m(6) * (t107 * t7 + t108 * t8 + t173 * t45) - t218 + (-t32 * mrSges(5,3) - pkin(10) * t70 + t250) * t190 + (-t2 * t91 + t3 * t90) * mrSges(7,3) + (t144 * t8 - t145 * t7) * mrSges(6,3) - Ifges(4,3) * t222 + t211 * t135 + t40 * t20 + t41 * t19 + t24 * t49 - pkin(3) * t60 + t71 * mrSges(4,1) - t72 * mrSges(4,2) + t45 * t97 + t107 * t44 + t108 * t43 + t116 * t12 + t65 * t153 + t101 * t156 / 0.2e1 + t173 * t31; t99 * t240 + t36 * t244 + t35 * t245 + t52 * t246 + t51 * t247 + (-mrSges(4,1) + t153) * t181 + t77 * t237 + t76 * t238 + t100 * t239 + (-t14 * t91 + t15 * t90) * mrSges(7,3) + t217 + (t144 * t55 - t145 * t54) * mrSges(6,3) + (-t114 * mrSges(5,3) - pkin(10) * t150 - t191 * t156 / 0.2e1 + t241) * t190 + (t115 * mrSges(5,3) + pkin(10) * t149 + t191 * t236 + t242) * t195 + m(5) * (-pkin(3) * t181 + (-t114 * t190 + t115 * t195) * pkin(10)) + m(7) * (t116 * t93 + t14 * t40 + t15 * t41) + m(6) * (t107 * t54 + t108 * t55 + t151 * t173) + (-pkin(9) * mrSges(4,2) - t211) * t196 + t41 * t67 + t40 * t68 + t93 * t49 + t108 * t109 + t107 * t110 + t116 * t38 - pkin(3) * t139 + t151 * t97 + t173 * t82; -0.2e1 * pkin(3) * t153 + t145 * t100 + 0.2e1 * t116 * t49 + t144 * t99 + t195 * t156 + t190 * t158 + 0.2e1 * t173 * t97 + t90 * t51 + t91 * t52 + Ifges(4,3) + m(7) * (t116 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t107 ^ 2 + t108 ^ 2 + t173 ^ 2) + m(5) * (t216 * pkin(10) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t40 * t91 + t41 * t90) * mrSges(7,3) + 0.2e1 * (-t107 * t145 + t108 * t144) * mrSges(6,3) + 0.2e1 * t216 * pkin(10) * mrSges(5,3); t202 + (m(6) * (t189 * t8 + t194 * t7) + t189 * t43 + t194 * t44) * pkin(4) + m(7) * (t133 * t2 + t134 * t3) + t32 * mrSges(5,1) - t33 * mrSges(5,2) + t133 * t20 + t134 * t19 + t46; t204 + m(7) * (t133 * t14 + t134 * t15) + t214 * t196 + (m(6) * (t189 * t55 + t194 * t54) + t189 * t109 + t194 * t110) * pkin(4) + t114 * mrSges(5,1) - t115 * mrSges(5,2) + t133 * t68 + t134 * t67 + t210; m(7) * (t133 * t40 + t134 * t41) - t209 * pkin(10) + (-t133 * t91 + t134 * t90) * mrSges(7,3) + (m(6) * (t107 * t194 + t108 * t189) + (t144 * t189 - t145 * t194) * mrSges(6,3)) * pkin(4) + t203 + t155; -0.2e1 * t227 + 0.2e1 * t125 + 0.2e1 * t205 + m(7) * (t133 ^ 2 + t134 ^ 2) + m(6) * (t189 ^ 2 + t194 ^ 2) * pkin(4) ^ 2 - t214; (m(7) * (t188 * t3 + t193 * t2) + t188 * t19 + t193 * t20) * pkin(5) + t202; t232 * t196 + (t193 * t68 + m(7) * (t14 * t193 + t15 * t188) + t188 * t67) * pkin(5) + t204; (m(7) * (t188 * t41 + t193 * t40) + (t188 * t90 - t193 * t91) * mrSges(7,3)) * pkin(5) + t203; Ifges(6,3) + t174 + t205 + (m(7) * (t133 * t193 + t134 * t188) - t224) * pkin(5) + t212; -0.2e1 * t215 + 0.2e1 * t174 + m(7) * (t188 ^ 2 + t193 ^ 2) * pkin(5) ^ 2 - t232; t206; t208 - t228; t207; t212; Ifges(7,3) + t174 - t215; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
