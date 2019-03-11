% Calculate joint inertia matrix for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:58
% EndTime: 2019-03-09 07:17:03
% DurationCPUTime: 2.34s
% Computational Cost: add. (7478->364), mult. (7165->522), div. (0->0), fcn. (7418->10), ass. (0->183)
t172 = qJ(3) + qJ(4);
t164 = qJ(5) + t172;
t157 = sin(t164);
t158 = cos(t164);
t173 = sin(qJ(6));
t176 = cos(qJ(6));
t83 = t157 * rSges(7,3) + (rSges(7,1) * t176 - rSges(7,2) * t173) * t158;
t260 = t158 * pkin(5) + t157 * pkin(10) + t83;
t175 = sin(qJ(1));
t178 = cos(qJ(1));
t259 = t175 * t178;
t258 = (rSges(6,1) * t157 + rSges(6,2) * t158) * t178;
t161 = cos(t172);
t160 = sin(t172);
t245 = rSges(5,1) * t160;
t257 = (rSges(5,2) * t161 + t245) * t178;
t174 = sin(qJ(3));
t177 = cos(qJ(3));
t256 = (rSges(4,1) * t174 + rSges(4,2) * t177) * t178;
t170 = t175 ^ 2;
t171 = t178 ^ 2;
t179 = -pkin(8) - pkin(7);
t255 = t175 / 0.2e1;
t254 = t178 / 0.2e1;
t253 = pkin(3) * t177;
t252 = pkin(4) * t161;
t251 = pkin(5) * t157;
t250 = t174 * pkin(3);
t225 = t178 * t173;
t227 = t175 * t176;
t118 = t157 * t225 + t227;
t224 = t178 * t176;
t228 = t175 * t173;
t119 = -t157 * t224 + t228;
t206 = -t119 * rSges(7,1) - t118 * rSges(7,2);
t231 = t158 * t178;
t66 = rSges(7,3) * t231 - t206;
t249 = t178 * t66 + t171 * (pkin(10) * t158 - t251);
t79 = Icges(7,3) * t157 + (Icges(7,5) * t176 - Icges(7,6) * t173) * t158;
t81 = Icges(7,5) * t157 + (Icges(7,1) * t176 - Icges(7,4) * t173) * t158;
t248 = t158 * t176 * t81 + t157 * t79;
t221 = t175 * t179 + t178 * t250;
t136 = pkin(4) * t160 + t250;
t169 = -pkin(9) + t179;
t222 = -t178 * t136 - t175 * t169;
t71 = t178 * (t221 + t222);
t82 = t178 * (t175 * rSges(6,3) - t258);
t247 = t71 + t82;
t129 = t175 * t136;
t229 = t174 * t175;
t153 = pkin(3) * t229;
t75 = t129 - t153 + (-t169 + t179) * t178;
t232 = t158 * t175;
t233 = t157 * t175;
t93 = rSges(6,1) * t233 + rSges(6,2) * t232 + t178 * rSges(6,3);
t246 = -t75 - t93;
t80 = Icges(7,6) * t157 + (Icges(7,4) * t176 - Icges(7,2) * t173) * t158;
t244 = t173 * t80;
t116 = -t157 * t228 + t224;
t117 = t157 * t227 + t225;
t57 = Icges(7,5) * t117 + Icges(7,6) * t116 - Icges(7,3) * t232;
t59 = Icges(7,4) * t117 + Icges(7,2) * t116 - Icges(7,6) * t232;
t61 = Icges(7,1) * t117 + Icges(7,4) * t116 - Icges(7,5) * t232;
t24 = t157 * t57 + (-t173 * t59 + t176 * t61) * t158;
t243 = t24 * t178;
t58 = Icges(7,5) * t119 + Icges(7,6) * t118 + Icges(7,3) * t231;
t60 = Icges(7,4) * t119 + Icges(7,2) * t118 + Icges(7,6) * t231;
t62 = Icges(7,1) * t119 + Icges(7,4) * t118 + Icges(7,5) * t231;
t25 = t157 * t58 + (-t173 * t60 + t176 * t62) * t158;
t242 = t25 * t175;
t140 = pkin(5) * t233;
t223 = t117 * rSges(7,1) + t116 * rSges(7,2);
t65 = -rSges(7,3) * t232 + t223;
t241 = pkin(10) * t232 - t140 - t65;
t63 = t260 * t175;
t239 = Icges(4,4) * t174;
t238 = Icges(4,4) * t177;
t237 = Icges(5,4) * t160;
t236 = Icges(5,4) * t161;
t235 = Icges(6,4) * t157;
t234 = Icges(6,4) * t158;
t230 = t161 * t175;
t226 = t175 * t177;
t127 = t158 * rSges(6,1) - t157 * rSges(6,2);
t149 = pkin(4) * t230;
t84 = t175 * t127 + t149;
t220 = t178 * pkin(1) + t175 * qJ(2);
t219 = t170 + t171;
t218 = t71 + t249;
t217 = -t75 + t241;
t50 = t149 + t63;
t163 = t178 * qJ(2);
t216 = t163 - t222;
t103 = rSges(5,2) * t230 + t178 * rSges(5,3) + t175 * t245;
t215 = rSges(4,1) * t229 + rSges(4,2) * t226 + t178 * rSges(4,3);
t214 = (-rSges(7,3) - pkin(10)) * t158;
t18 = t118 * t59 + t119 * t61 + t57 * t231;
t19 = t118 * t60 + t119 * t62 + t58 * t231;
t10 = t19 * t175 + t18 * t178;
t191 = Icges(6,5) * t157 + Icges(6,6) * t158;
t87 = Icges(6,3) * t178 + t191 * t175;
t88 = Icges(6,3) * t175 - t191 * t178;
t16 = t116 * t59 + t117 * t61 - t57 * t232;
t17 = t116 * t60 + t117 * t62 - t58 * t232;
t9 = t16 * t178 + t17 * t175;
t213 = (t171 * t87 + t9) * t178 + (t170 * t88 + t10 + (t175 * t87 + t178 * t88) * t178) * t175;
t30 = t116 * t80 + t117 * t81 - t79 * t232;
t3 = t30 * t157 + (-t16 * t175 + t17 * t178) * t158;
t31 = t118 * t80 + t119 * t81 + t79 * t231;
t4 = t31 * t157 + (-t175 * t18 + t178 * t19) * t158;
t212 = t3 * t254 + t4 * t255 + t157 * (t242 + t243) / 0.2e1 - t9 * t232 / 0.2e1 + t10 * t231 / 0.2e1;
t211 = -t252 - t253;
t192 = Icges(5,5) * t160 + Icges(5,6) * t161;
t97 = Icges(5,3) * t178 + t192 * t175;
t98 = Icges(5,3) * t175 - t192 * t178;
t210 = t178 * (t171 * t97 + t98 * t259) + t175 * (t170 * t98 + t97 * t259) + t213;
t154 = pkin(3) * t226;
t73 = t154 + t84;
t74 = (-t127 + t211) * t178;
t203 = t73 * t175 - t74 * t178;
t85 = (-t127 - t252) * t178;
t202 = t84 * t175 - t85 * t178;
t135 = t161 * rSges(5,1) - t160 * rSges(5,2);
t94 = t175 * t135 + t154;
t95 = (-t135 - t253) * t178;
t201 = t94 * t175 - t95 * t178;
t199 = Icges(4,1) * t174 + t238;
t198 = Icges(5,1) * t160 + t236;
t197 = Icges(6,1) * t157 + t234;
t196 = Icges(4,2) * t177 + t239;
t195 = Icges(5,2) * t161 + t237;
t194 = Icges(6,2) * t158 + t235;
t193 = Icges(4,5) * t174 + Icges(4,6) * t177;
t125 = -Icges(6,2) * t157 + t234;
t126 = Icges(6,1) * t158 - t235;
t187 = t125 * t158 + t126 * t157;
t133 = -Icges(5,2) * t160 + t236;
t134 = Icges(5,1) * t161 - t237;
t186 = t133 * t161 + t134 * t160;
t185 = -t178 * t169 + t129 + t220;
t76 = t163 + t256 + (-rSges(4,3) - pkin(1) - pkin(7)) * t175;
t77 = t178 * pkin(7) + t215 + t220;
t184 = m(4) * (t175 * t76 - t178 * t77);
t68 = t163 + t257 + (-rSges(5,3) - pkin(1)) * t175 + t221;
t69 = -t178 * t179 + t103 + t153 + t220;
t183 = m(5) * (t175 * t68 - t178 * t69);
t53 = t258 + (-rSges(6,3) - pkin(1)) * t175 + t216;
t54 = t185 + t93;
t182 = m(6) * (t175 * t53 - t178 * t54);
t124 = Icges(6,5) * t158 - Icges(6,6) * t157;
t181 = t242 / 0.2e1 + t243 / 0.2e1 + (t175 * t124 - t157 * (Icges(6,6) * t175 - t194 * t178) + t158 * (Icges(6,5) * t175 - t197 * t178) - t187 * t178 + t31) * t255 + (t178 * t124 - t157 * (Icges(6,6) * t178 + t194 * t175) + t158 * (Icges(6,5) * t178 + t197 * t175) + t187 * t175 + t30) * t254;
t132 = Icges(5,5) * t161 - Icges(5,6) * t160;
t180 = t181 + (-t160 * (Icges(5,6) * t175 - t195 * t178) + t161 * (Icges(5,5) * t175 - t198 * t178) + t175 * t132 - t186 * t178) * t255 + (t161 * (Icges(5,5) * t178 + t198 * t175) + t178 * t132 - t160 * (Icges(5,6) * t178 + t195 * t175) + t186 * t175) * t254;
t146 = t178 * rSges(2,1) - t175 * rSges(2,2);
t145 = t177 * rSges(4,1) - t174 * rSges(4,2);
t144 = -t175 * rSges(2,1) - t178 * rSges(2,2);
t122 = t153 + (-pkin(7) - t179) * t178;
t121 = -t178 * rSges(3,2) + t175 * rSges(3,3) + t220;
t120 = t178 * rSges(3,3) + t163 + (rSges(3,2) - pkin(1)) * t175;
t109 = Icges(4,3) * t175 - t193 * t178;
t108 = Icges(4,3) * t178 + t193 * t175;
t104 = t178 * (-t175 * pkin(7) - t221);
t86 = t178 * (t175 * rSges(5,3) - t257);
t67 = -t175 * t215 + (t175 * rSges(4,3) - t256) * t178;
t64 = t260 * t178;
t56 = -t175 * t103 + t86;
t52 = -t175 * t93 + t82;
t51 = (-t260 - t252) * t178;
t45 = (t211 - t260) * t178;
t44 = t154 + t50;
t39 = t104 + t86 + (-t103 - t122) * t175;
t38 = -t157 * t66 + t83 * t231;
t37 = t157 * t65 + t83 * t232;
t36 = t175 * t214 + t140 + t185 + t223;
t35 = -t175 * pkin(1) + (t214 + t251) * t178 + t206 + t216;
t34 = (-t158 * t244 + t248) * t157;
t33 = (-t175 * t66 - t178 * t65) * t158;
t32 = t246 * t175 + t247;
t27 = t241 * t175 + t249;
t26 = t104 + (-t122 + t246) * t175 + t247;
t13 = t217 * t175 + t218;
t12 = t104 + (-t122 + t217) * t175 + t218;
t1 = [-t157 * t125 - t160 * t133 + t161 * t134 - t174 * (-Icges(4,2) * t174 + t238) + t177 * (Icges(4,1) * t177 - t239) + Icges(3,1) + Icges(2,3) + (t126 - t244) * t158 + m(7) * (t35 ^ 2 + t36 ^ 2) + m(6) * (t53 ^ 2 + t54 ^ 2) + m(5) * (t68 ^ 2 + t69 ^ 2) + m(4) * (t76 ^ 2 + t77 ^ 2) + m(2) * (t144 ^ 2 + t146 ^ 2) + m(3) * (t120 ^ 2 + t121 ^ 2) + t248; m(7) * (t175 * t35 - t178 * t36) + t182 + t183 + t184 + m(3) * (t175 * t120 - t178 * t121); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t219; t180 + (t170 / 0.2e1 + t171 / 0.2e1) * (Icges(4,5) * t177 - Icges(4,6) * t174) + t145 * t184 + m(7) * (t44 * t35 + t45 * t36) + m(6) * (t73 * t53 + t74 * t54) + m(5) * (t94 * t68 + t95 * t69) + (-t174 * (Icges(4,6) * t175 - t196 * t178) + t177 * (Icges(4,5) * t175 - t199 * t178)) * t255 + (-t174 * (Icges(4,6) * t178 + t196 * t175) + t177 * (Icges(4,5) * t178 + t199 * t175)) * t254; m(5) * t201 + m(6) * t203 + m(7) * (t44 * t175 - t45 * t178) + m(4) * t219 * t145; m(7) * (t12 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t26 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t39 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(4) * (t219 * t145 ^ 2 + t67 ^ 2) + t175 * (t108 * t259 + t170 * t109) + t178 * (t171 * t108 + t109 * t259) + t210; t180 + m(7) * (t50 * t35 + t51 * t36) + m(6) * (t84 * t53 + t85 * t54) + t135 * t183; m(6) * t202 + m(7) * (t50 * t175 - t51 * t178) + m(5) * t219 * t135; m(7) * (t13 * t12 + t50 * t44 + t51 * t45) + m(6) * (t32 * t26 + t84 * t73 + t85 * t74) + m(5) * (t201 * t135 + t56 * t39) + t210; m(7) * (t13 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t32 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(5) * (t219 * t135 ^ 2 + t56 ^ 2) + t210; m(7) * (t63 * t35 - t64 * t36) + t127 * t182 + t181; m(7) * (t63 * t175 + t64 * t178) + m(6) * t219 * t127; m(7) * (t27 * t12 + t63 * t44 - t64 * t45) + m(6) * (t203 * t127 + t52 * t26) + t213; m(7) * (t27 * t13 + t63 * t50 - t64 * t51) + m(6) * (t202 * t127 + t52 * t32) + t213; m(6) * (t219 * t127 ^ 2 + t52 ^ 2) + m(7) * (t27 ^ 2 + t63 ^ 2 + t64 ^ 2) + t213; m(7) * (t38 * t35 + t37 * t36) + t34 + ((t25 / 0.2e1 + t31 / 0.2e1) * t178 + (-t24 / 0.2e1 - t30 / 0.2e1) * t175) * t158; m(7) * (t38 * t175 - t37 * t178); m(7) * (t33 * t12 + t37 * t45 + t38 * t44) + t212; m(7) * (t33 * t13 + t37 * t51 + t38 * t50) + t212; m(7) * (t33 * t27 - t37 * t64 + t38 * t63) + t212; t157 * t34 + m(7) * (t33 ^ 2 + t37 ^ 2 + t38 ^ 2) + (-t175 * t3 + t178 * t4 + t157 * (-t175 * t24 + t178 * t25)) * t158;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
