% Calculate joint inertia matrix for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:07
% EndTime: 2019-12-31 21:16:13
% DurationCPUTime: 2.05s
% Computational Cost: add. (5557->348), mult. (5791->511), div. (0->0), fcn. (6124->10), ass. (0->171)
t162 = cos(pkin(9));
t147 = pkin(4) * t162 + pkin(3);
t163 = -pkin(8) - qJ(4);
t161 = sin(pkin(9));
t165 = sin(qJ(1));
t206 = t161 * t165;
t160 = qJ(2) + qJ(3);
t153 = cos(t160);
t167 = cos(qJ(1));
t207 = t153 * t167;
t152 = sin(t160);
t209 = t152 * t167;
t157 = pkin(9) + qJ(5);
t150 = sin(t157);
t151 = cos(t157);
t111 = -t150 * t207 + t151 * t165;
t112 = t150 * t165 + t151 * t207;
t61 = t112 * rSges(6,1) + t111 * rSges(6,2) + rSges(6,3) * t209;
t233 = pkin(4) * t206 + t147 * t207 - t163 * t209 + t61;
t201 = qJ(4) + t163;
t223 = -pkin(3) + t147;
t85 = -rSges(6,3) * t153 + (rSges(6,1) * t151 - rSges(6,2) * t150) * t152;
t232 = -t152 * t223 - t153 * t201 - t85;
t203 = t162 * t167;
t121 = -t153 * t206 - t203;
t204 = t162 * t165;
t205 = t161 * t167;
t122 = t153 * t204 - t205;
t159 = t167 ^ 2;
t211 = Icges(4,4) * t153;
t180 = -Icges(4,2) * t152 + t211;
t101 = Icges(4,6) * t165 + t167 * t180;
t212 = Icges(4,4) * t152;
t182 = Icges(4,1) * t153 - t212;
t103 = Icges(4,5) * t165 + t167 * t182;
t176 = -t101 * t152 + t103 * t153;
t100 = -Icges(4,6) * t167 + t165 * t180;
t102 = -Icges(4,5) * t167 + t165 * t182;
t177 = t100 * t152 - t102 * t153;
t210 = t152 * t165;
t69 = Icges(5,5) * t122 + Icges(5,6) * t121 + Icges(5,3) * t210;
t123 = -t153 * t205 + t204;
t124 = t153 * t203 + t206;
t70 = Icges(5,5) * t124 + Icges(5,6) * t123 + Icges(5,3) * t209;
t71 = Icges(5,4) * t122 + Icges(5,2) * t121 + Icges(5,6) * t210;
t72 = Icges(5,4) * t124 + Icges(5,2) * t123 + Icges(5,6) * t209;
t73 = Icges(5,1) * t122 + Icges(5,4) * t121 + Icges(5,5) * t210;
t74 = Icges(5,1) * t124 + Icges(5,4) * t123 + Icges(5,5) * t209;
t208 = t153 * t165;
t109 = -t150 * t208 - t151 * t167;
t110 = -t150 * t167 + t151 * t208;
t54 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t210;
t56 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t210;
t58 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t210;
t15 = t109 * t56 + t110 * t58 + t210 * t54;
t55 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t209;
t57 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t209;
t59 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t209;
t16 = t109 * t57 + t110 * t59 + t210 * t55;
t8 = -t15 * t167 + t16 * t165;
t178 = Icges(4,5) * t153 - Icges(4,6) * t152;
t98 = -Icges(4,3) * t167 + t165 * t178;
t99 = Icges(4,3) * t165 + t167 * t178;
t231 = -t8 + (t121 * t71 + t122 * t73 + t69 * t210) * t167 - t159 * t98 + (-t121 * t72 - t122 * t74 - t70 * t210 - (t177 - t99) * t167 - t176 * t165) * t165;
t230 = t153 ^ 2;
t158 = t165 ^ 2;
t229 = m(5) / 0.2e1;
t228 = m(6) / 0.2e1;
t227 = t165 / 0.2e1;
t226 = -t167 / 0.2e1;
t164 = sin(qJ(2));
t225 = pkin(2) * t164;
t224 = pkin(3) * t153;
t166 = cos(qJ(2));
t148 = pkin(2) * t166 + pkin(1);
t140 = t167 * t148;
t156 = t167 * pkin(6);
t168 = -pkin(7) - pkin(6);
t202 = t167 * t168;
t222 = t165 * (t202 + t156 + (-pkin(1) + t148) * t165) + t167 * (-t167 * pkin(1) + t140 + (-pkin(6) - t168) * t165);
t172 = rSges(4,1) * t207 - rSges(4,2) * t209 + t165 * rSges(4,3);
t187 = rSges(4,1) * t153 - rSges(4,2) * t152;
t64 = t165 * (-t167 * rSges(4,3) + t165 * t187) + t167 * t172;
t221 = rSges(3,1) * t166;
t220 = rSges(3,2) * t164;
t83 = -Icges(6,6) * t153 + (Icges(6,4) * t151 - Icges(6,2) * t150) * t152;
t219 = t150 * t83;
t218 = t167 * rSges(3,3);
t22 = -t153 * t54 + (-t150 * t56 + t151 * t58) * t152;
t217 = t22 * t167;
t23 = -t153 * t55 + (-t150 * t57 + t151 * t59) * t152;
t216 = t23 * t165;
t130 = pkin(3) * t152 - qJ(4) * t153;
t95 = -rSges(5,3) * t153 + (rSges(5,1) * t162 - rSges(5,2) * t161) * t152;
t215 = -t130 - t95;
t214 = Icges(3,4) * t164;
t213 = Icges(3,4) * t166;
t199 = pkin(3) * t207 + qJ(4) * t209;
t200 = t158 * (qJ(4) * t152 + t224) + t167 * t199;
t198 = t165 * rSges(3,3) + t167 * t221;
t197 = t158 + t159;
t17 = t111 * t56 + t112 * t58 + t209 * t54;
t18 = t111 * t57 + t112 * t59 + t209 * t55;
t9 = t165 * t18 - t167 * t17;
t196 = (t9 + t158 * t99 + (t123 * t72 + t124 * t74 + t70 * t209) * t165 + (-t123 * t71 - t124 * t73 - t69 * t209 + (t176 - t98) * t165 + t177 * t167) * t167) * t165;
t195 = -t130 + t232;
t194 = t124 * rSges(5,1) + t123 * rSges(5,2) + rSges(5,3) * t209;
t193 = -t130 - t225;
t131 = rSges(4,1) * t152 + rSges(4,2) * t153;
t192 = -t131 - t225;
t191 = -t165 * t168 + t140;
t186 = -t122 * rSges(5,1) - t121 * rSges(5,2);
t28 = t165 * (rSges(5,3) * t210 - t186) + t167 * t194 + t200;
t82 = -Icges(6,3) * t153 + (Icges(6,5) * t151 - Icges(6,6) * t150) * t152;
t84 = -Icges(6,5) * t153 + (Icges(6,1) * t151 - Icges(6,4) * t150) * t152;
t31 = t109 * t83 + t110 * t84 + t210 * t82;
t3 = -t153 * t31 + (t15 * t165 + t16 * t167) * t152;
t32 = t111 * t83 + t112 * t84 + t209 * t82;
t4 = -t153 * t32 + (t165 * t17 + t167 * t18) * t152;
t190 = t8 * t210 / 0.2e1 + t3 * t226 + t4 * t227 - t153 * (t216 - t217) / 0.2e1 + t9 * t209 / 0.2e1;
t189 = t193 - t95;
t188 = -t220 + t221;
t185 = -t110 * rSges(6,1) - t109 * rSges(6,2);
t184 = t193 + t232;
t183 = Icges(3,1) * t166 - t214;
t181 = -Icges(3,2) * t164 + t213;
t179 = Icges(3,5) * t166 - Icges(3,6) * t164;
t128 = Icges(4,2) * t153 + t212;
t129 = Icges(4,1) * t152 + t211;
t173 = -t128 * t152 + t129 * t153;
t60 = rSges(6,3) * t210 - t185;
t14 = t200 + (-t199 + t233) * t167 + (-pkin(4) * t205 + t60 + (-t152 * t201 + t153 * t223) * t165) * t165;
t170 = t167 * t231 + t196;
t127 = Icges(4,5) * t152 + Icges(4,6) * t153;
t92 = -Icges(5,3) * t153 + (Icges(5,5) * t162 - Icges(5,6) * t161) * t152;
t93 = -Icges(5,6) * t153 + (Icges(5,4) * t162 - Icges(5,2) * t161) * t152;
t94 = -Icges(5,5) * t153 + (Icges(5,1) * t162 - Icges(5,4) * t161) * t152;
t169 = -t217 / 0.2e1 + t216 / 0.2e1 + (t123 * t93 + t124 * t94 + t127 * t165 + t173 * t167 + t92 * t209 + t32 + (t101 - t70) * t153 + (-t161 * t72 + t162 * t74 + t103) * t152) * t227 + (t121 * t93 + t122 * t94 - t127 * t167 + t173 * t165 + t92 * t210 + t31 + (t100 - t69) * t153 + (-t161 * t71 + t162 * t73 + t102) * t152) * t226;
t138 = rSges(2,1) * t167 - rSges(2,2) * t165;
t137 = -rSges(2,1) * t165 - rSges(2,2) * t167;
t136 = rSges(3,1) * t164 + rSges(3,2) * t166;
t114 = Icges(3,3) * t165 + t167 * t179;
t113 = -Icges(3,3) * t167 + t165 * t179;
t97 = t192 * t167;
t96 = t192 * t165;
t87 = t165 * pkin(6) + (pkin(1) - t220) * t167 + t198;
t86 = t218 + t156 + (-pkin(1) - t188) * t165;
t80 = t172 + t191;
t79 = (rSges(4,3) - t168) * t167 + (-t148 - t187) * t165;
t78 = t152 * t151 * t84;
t77 = t167 * (-t167 * t220 + t198) + (t165 * t188 - t218) * t165;
t76 = t215 * t167;
t75 = t215 * t165;
t66 = t189 * t167;
t65 = t189 * t165;
t47 = t191 + t194 + t199;
t46 = -t202 + (-t224 - t148 + (-rSges(5,3) - qJ(4)) * t152) * t165 + t186;
t45 = t195 * t167;
t44 = t195 * t165;
t43 = t184 * t167;
t42 = t184 * t165;
t41 = -t153 * t61 - t209 * t85;
t40 = t153 * t60 + t210 * t85;
t39 = t191 + t233;
t38 = (pkin(4) * t161 - t168) * t167 + (-t147 * t153 - t148 + (-rSges(6,3) + t163) * t152) * t165 + t185;
t37 = t64 + t222;
t36 = -t152 * t219 - t153 * t82 + t78;
t33 = (-t165 * t61 + t167 * t60) * t152;
t19 = t28 + t222;
t13 = t14 + t222;
t1 = [t166 * (Icges(3,2) * t166 + t214) + t164 * (Icges(3,1) * t164 + t213) + Icges(2,3) + t78 + (-t82 + t128 - t92) * t153 + (-t161 * t93 + t162 * t94 + t129 - t219) * t152 + m(6) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t46 ^ 2 + t47 ^ 2) + m(4) * (t79 ^ 2 + t80 ^ 2) + m(3) * (t86 ^ 2 + t87 ^ 2) + m(2) * (t137 ^ 2 + t138 ^ 2); ((-Icges(3,6) * t167 + t165 * t181) * t166 + (-Icges(3,5) * t167 + t165 * t183) * t164) * t226 + ((Icges(3,6) * t165 + t167 * t181) * t166 + (Icges(3,5) * t165 + t167 * t183) * t164) * t227 + t169 + m(3) * (-t165 * t87 - t167 * t86) * t136 + (t158 / 0.2e1 + t159 / 0.2e1) * (Icges(3,5) * t164 + Icges(3,6) * t166) + m(6) * (t38 * t43 + t39 * t42) + m(5) * (t46 * t66 + t47 * t65) + m(4) * (t79 * t97 + t80 * t96); m(6) * (t13 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t19 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(4) * (t37 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(3) * (t136 ^ 2 * t197 + t77 ^ 2) + t165 * t114 * t158 + t196 + (-t113 * t159 + (-t113 * t165 + t114 * t167) * t165 + t231) * t167; t169 + m(6) * (t38 * t45 + t39 * t44) + m(5) * (t46 * t76 + t47 * t75) + m(4) * (-t165 * t80 - t167 * t79) * t131; m(6) * (t13 * t14 + t42 * t44 + t43 * t45) + m(5) * (t19 * t28 + t65 * t75 + t66 * t76) + m(4) * (t37 * t64 + (-t165 * t96 - t167 * t97) * t131) + t170; m(6) * (t14 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(5) * (t28 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t131 ^ 2 * t197 + t64 ^ 2) + t170; 0.2e1 * ((t165 * t39 + t167 * t38) * t228 + (t165 * t47 + t167 * t46) * t229) * t152; m(6) * (-t13 * t153 + (t165 * t42 + t167 * t43) * t152) + m(5) * (-t153 * t19 + (t165 * t65 + t167 * t66) * t152); m(6) * (-t14 * t153 + (t165 * t44 + t167 * t45) * t152) + m(5) * (-t153 * t28 + (t165 * t75 + t167 * t76) * t152); 0.2e1 * (t229 + t228) * (t152 ^ 2 * t197 + t230); m(6) * (t38 * t40 + t39 * t41) - t36 * t153 + ((t32 / 0.2e1 + t23 / 0.2e1) * t167 + (t31 / 0.2e1 + t22 / 0.2e1) * t165) * t152; m(6) * (t13 * t33 + t40 * t43 + t41 * t42) + t190; m(6) * (t14 * t33 + t40 * t45 + t41 * t44) + t190; m(6) * (-t153 * t33 + (t165 * t41 + t167 * t40) * t152); m(6) * (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) + t230 * t36 + (t167 * t4 + t165 * t3 - t153 * (t165 * t22 + t167 * t23)) * t152;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
