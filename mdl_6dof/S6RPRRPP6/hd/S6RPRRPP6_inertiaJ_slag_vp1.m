% Calculate joint inertia matrix for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:23
% EndTime: 2019-03-09 04:46:31
% DurationCPUTime: 3.13s
% Computational Cost: add. (5007->397), mult. (7732->570), div. (0->0), fcn. (8245->8), ass. (0->198)
t186 = cos(qJ(3));
t270 = -0.2e1 * t186;
t269 = Icges(4,5) * t186;
t183 = sin(qJ(3));
t268 = Icges(4,6) * t183;
t267 = t269 / 0.2e1 - t268 / 0.2e1;
t177 = qJ(4) + pkin(9);
t170 = sin(t177);
t187 = cos(qJ(1));
t229 = t187 * t170;
t171 = cos(t177);
t184 = sin(qJ(1));
t234 = t184 * t171;
t123 = t183 * t229 + t234;
t228 = t187 * t171;
t235 = t184 * t170;
t125 = -t183 * t228 + t235;
t230 = t186 * t187;
t265 = rSges(7,3) + qJ(6);
t266 = rSges(7,1) + pkin(5);
t247 = rSges(7,2) * t230 - t265 * t123 + t266 * t125;
t107 = Icges(7,6) * t183 + (Icges(7,5) * t171 + Icges(7,3) * t170) * t186;
t108 = Icges(6,3) * t183 + (Icges(6,5) * t171 - Icges(6,6) * t170) * t186;
t109 = Icges(7,2) * t183 + (Icges(7,4) * t171 + Icges(7,6) * t170) * t186;
t111 = Icges(7,4) * t183 + (Icges(7,1) * t171 + Icges(7,5) * t170) * t186;
t112 = Icges(6,5) * t183 + (Icges(6,1) * t171 - Icges(6,4) * t170) * t186;
t182 = sin(qJ(4));
t185 = cos(qJ(4));
t126 = Icges(5,3) * t183 + (Icges(5,5) * t185 - Icges(5,6) * t182) * t186;
t132 = Icges(5,5) * t183 + (Icges(5,1) * t185 - Icges(5,4) * t182) * t186;
t240 = t170 * t186;
t264 = t107 * t240 + (t108 + t109 + t126) * t183 + (t185 * t132 + (t111 + t112) * t171) * t186;
t173 = t187 * qJ(2);
t169 = t185 * pkin(4) + pkin(3);
t181 = -qJ(5) - pkin(8);
t236 = t183 * t187;
t220 = t169 * t236 + t181 * t230;
t256 = -pkin(1) - pkin(7);
t188 = t173 + (-pkin(4) * t182 + t256) * t184 + t220;
t77 = t125 * rSges(6,1) + t123 * rSges(6,2) + rSges(6,3) * t230;
t51 = t188 - t77;
t219 = t187 * pkin(1) + t184 * qJ(2);
t209 = t187 * pkin(7) + t219;
t227 = t187 * t182;
t231 = t184 * t186;
t237 = t183 * t184;
t211 = pkin(4) * t227 + t169 * t237 + t181 * t231;
t191 = t209 + t211;
t121 = t183 * t235 - t228;
t122 = t183 * t234 + t229;
t75 = t122 * rSges(6,1) - t121 * rSges(6,2) - rSges(6,3) * t231;
t52 = t191 + t75;
t263 = m(6) * (t184 * t51 - t187 * t52);
t39 = t188 - t247;
t249 = -rSges(7,2) * t231 + t265 * t121 + t266 * t122;
t40 = t191 + t249;
t262 = m(7) * (t184 * t39 - t187 * t40);
t261 = (rSges(4,1) * t183 + rSges(4,2) * t186) * t187;
t62 = Icges(7,5) * t122 - Icges(7,6) * t231 + Icges(7,3) * t121;
t66 = Icges(7,4) * t122 - Icges(7,2) * t231 + Icges(7,6) * t121;
t70 = Icges(7,1) * t122 - Icges(7,4) * t231 + Icges(7,5) * t121;
t19 = t121 * t62 + t122 * t70 - t231 * t66;
t63 = Icges(7,5) * t125 + Icges(7,6) * t230 - Icges(7,3) * t123;
t67 = Icges(7,4) * t125 + Icges(7,2) * t230 - Icges(7,6) * t123;
t71 = Icges(7,1) * t125 + Icges(7,4) * t230 - Icges(7,5) * t123;
t20 = t121 * t63 + t122 * t71 - t231 * t67;
t64 = Icges(6,5) * t122 - Icges(6,6) * t121 - Icges(6,3) * t231;
t68 = Icges(6,4) * t122 - Icges(6,2) * t121 - Icges(6,6) * t231;
t72 = Icges(6,1) * t122 - Icges(6,4) * t121 - Icges(6,5) * t231;
t21 = -t121 * t68 + t122 * t72 - t231 * t64;
t65 = Icges(6,5) * t125 + Icges(6,6) * t123 + Icges(6,3) * t230;
t69 = Icges(6,4) * t125 + Icges(6,2) * t123 + Icges(6,6) * t230;
t73 = Icges(6,1) * t125 + Icges(6,4) * t123 + Icges(6,5) * t230;
t22 = -t121 * t69 + t122 * t73 - t231 * t65;
t226 = t187 * t185;
t233 = t184 * t182;
t143 = -t183 * t233 + t226;
t232 = t184 * t185;
t144 = t183 * t232 + t227;
t83 = Icges(5,5) * t144 + Icges(5,6) * t143 - Icges(5,3) * t231;
t85 = Icges(5,4) * t144 + Icges(5,2) * t143 - Icges(5,6) * t231;
t87 = Icges(5,1) * t144 + Icges(5,4) * t143 - Icges(5,5) * t231;
t31 = t143 * t85 + t144 * t87 - t231 * t83;
t145 = t183 * t227 + t232;
t146 = -t183 * t226 + t233;
t84 = Icges(5,5) * t146 + Icges(5,6) * t145 + Icges(5,3) * t230;
t86 = Icges(5,4) * t146 + Icges(5,2) * t145 + Icges(5,6) * t230;
t88 = Icges(5,1) * t146 + Icges(5,4) * t145 + Icges(5,5) * t230;
t32 = t143 * t86 + t144 * t88 - t231 * t84;
t42 = t121 * t107 - t109 * t231 + t122 * t111;
t110 = Icges(6,6) * t183 + (Icges(6,4) * t171 - Icges(6,2) * t170) * t186;
t43 = -t108 * t231 - t121 * t110 + t122 * t112;
t129 = Icges(5,6) * t183 + (Icges(5,4) * t185 - Icges(5,2) * t182) * t186;
t46 = -t126 * t231 + t143 * t129 + t144 * t132;
t260 = ((t20 + t22 + t32) * t187 + (-t19 - t21 - t31) * t184) * t186 + (t42 + t43 + t46) * t183;
t23 = -t123 * t62 + t125 * t70 + t230 * t66;
t24 = -t123 * t63 + t125 * t71 + t230 * t67;
t25 = t123 * t68 + t125 * t72 + t230 * t64;
t26 = t123 * t69 + t125 * t73 + t230 * t65;
t33 = t145 * t85 + t146 * t87 + t230 * t83;
t34 = t145 * t86 + t146 * t88 + t230 * t84;
t44 = -t123 * t107 + t109 * t230 + t125 * t111;
t45 = t108 * t230 + t123 * t110 + t125 * t112;
t47 = t126 * t230 + t145 * t129 + t146 * t132;
t259 = ((t24 + t26 + t34) * t187 + (-t23 - t25 - t33) * t184) * t186 + (t44 + t45 + t47) * t183;
t27 = t183 * t66 + (t170 * t62 + t171 * t70) * t186;
t29 = t183 * t64 + (-t170 * t68 + t171 * t72) * t186;
t37 = t183 * t83 + (-t182 * t85 + t185 * t87) * t186;
t258 = t27 + t29 + t37;
t28 = t183 * t67 + (t170 * t63 + t171 * t71) * t186;
t30 = t183 * t65 + (-t170 * t69 + t171 * t73) * t186;
t38 = t183 * t84 + (-t182 * t86 + t185 * t88) * t186;
t257 = t28 + t30 + t38;
t178 = t184 ^ 2;
t180 = t187 ^ 2;
t253 = t184 / 0.2e1;
t251 = t187 / 0.2e1;
t157 = t186 * rSges(4,1) - t183 * rSges(4,2);
t250 = m(4) * t157;
t166 = pkin(3) * t237;
t147 = -pkin(8) * t231 + t166;
t89 = -t147 + t211;
t248 = -t75 - t89;
t168 = pkin(3) * t236;
t206 = pkin(8) * t230 - t168;
t90 = pkin(4) * t233 - t206 - t220;
t246 = -t77 - t90;
t106 = (-pkin(3) + t169) * t186 + (-pkin(8) - t181) * t183;
t245 = t106 * t231 + t183 * t89;
t139 = t187 * t206;
t243 = t187 * t90 + t139;
t238 = t182 * t129;
t159 = t186 * pkin(3) + t183 * pkin(8);
t149 = t184 * t159;
t225 = t184 * t106 + t149;
t223 = -t106 - t159;
t222 = t183 * rSges(7,2) + (t265 * t170 + t266 * t171) * t186;
t221 = t144 * rSges(5,1) + t143 * rSges(5,2);
t218 = t178 + t180;
t217 = m(6) / 0.2e1 + m(7) / 0.2e1;
t216 = (-t110 * t240 - t186 * t238 + t264) * t183;
t215 = -t89 - t249;
t214 = -t90 - t247;
t210 = rSges(4,1) * t237 + rSges(4,2) * t231 + t187 * rSges(4,3);
t208 = (-rSges(5,3) - pkin(8)) * t186;
t207 = t222 * t186;
t204 = -t146 * rSges(5,1) - t145 * rSges(5,2);
t16 = t183 * t249 + t184 * t207 + t245;
t101 = t106 * t230;
t17 = t183 * t214 + t187 * t207 + t101;
t203 = t16 * t187 - t17 * t184;
t114 = t183 * rSges(6,3) + (rSges(6,1) * t171 - rSges(6,2) * t170) * t186;
t35 = t114 * t231 + t183 * t75 + t245;
t36 = t114 * t230 + t183 * t246 + t101;
t202 = t36 * t184 - t35 * t187;
t53 = t184 * t222 + t225;
t54 = (-t222 + t223) * t187;
t199 = t53 * t184 - t54 * t187;
t60 = t184 * t114 + t225;
t61 = (-t114 + t223) * t187;
t198 = t60 * t184 - t61 * t187;
t195 = Icges(4,5) * t183 + Icges(4,6) * t186;
t194 = t121 * t184 + t123 * t187;
t190 = t44 / 0.2e1 + t28 / 0.2e1 + t38 / 0.2e1 + t30 / 0.2e1 + t47 / 0.2e1 + t45 / 0.2e1;
t189 = -t46 / 0.2e1 - t43 / 0.2e1 - t27 / 0.2e1 - t37 / 0.2e1 - t29 / 0.2e1 - t42 / 0.2e1;
t179 = t186 ^ 2;
t158 = t187 * rSges(2,1) - t184 * rSges(2,2);
t156 = -t184 * rSges(2,1) - t187 * rSges(2,2);
t153 = -t268 + t269;
t138 = -t187 * rSges(3,2) + t184 * rSges(3,3) + t219;
t137 = t187 * rSges(3,3) + t173 + (rSges(3,2) - pkin(1)) * t184;
t136 = t183 * rSges(5,3) + (rSges(5,1) * t185 - rSges(5,2) * t182) * t186;
t128 = Icges(4,3) * t184 - t187 * t195;
t127 = Icges(4,3) * t187 + t184 * t195;
t99 = t209 + t210;
t98 = t173 + t261 + (-rSges(4,3) + t256) * t184;
t94 = (-t136 - t159) * t187;
t93 = t184 * t136 + t149;
t92 = rSges(5,3) * t230 - t204;
t91 = -rSges(5,3) * t231 + t221;
t78 = -t184 * t210 + (t184 * rSges(4,3) - t261) * t187;
t59 = t184 * t208 + t166 + t209 + t221;
t58 = t184 * t256 + t187 * t208 + t168 + t173 + t204;
t57 = t136 * t230 - t183 * t92;
t56 = t136 * t231 + t183 * t91;
t50 = (-t184 * t92 - t187 * t91) * t186;
t41 = t187 * t92 + t139 + (-t147 - t91) * t184;
t18 = (t184 * t246 + t187 * t248) * t186;
t15 = t187 * t77 + (-t147 + t248) * t184 + t243;
t14 = (t184 * t214 + t187 * t215) * t186;
t13 = t247 * t187 + (-t147 + t215) * t184 + t243;
t12 = t34 * t184 + t33 * t187;
t11 = t32 * t184 + t31 * t187;
t10 = t26 * t184 + t25 * t187;
t9 = t24 * t184 + t23 * t187;
t8 = t22 * t184 + t21 * t187;
t7 = t20 * t184 + t19 * t187;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t186 - t170 * t110 - t238) * t186 + m(7) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(6) * (t51 ^ 2 + t52 ^ 2) + m(4) * (t98 ^ 2 + t99 ^ 2) + m(3) * (t137 ^ 2 + t138 ^ 2) + m(2) * (t156 ^ 2 + t158 ^ 2) + t264 + (Icges(4,4) * t270 + Icges(4,2) * t183) * t183; t262 + m(5) * (t184 * t58 - t187 * t59) + t263 + m(4) * (t184 * t98 - t187 * t99) + m(3) * (t184 * t137 - t187 * t138); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t217) * t218; m(7) * (t53 * t39 + t54 * t40) + m(5) * (t93 * t58 + t94 * t59) + m(6) * (t60 * t51 + t61 * t52) + (t153 * t251 + t187 * t267 - t250 * t99 - t189) * t187 + (t153 * t253 + t184 * t267 + t250 * t98 + t190) * t184; m(5) * (t93 * t184 - t94 * t187) + m(6) * t198 + m(7) * t199 + t218 * t250; m(6) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(7) * (t13 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t41 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(4) * (t157 ^ 2 * t218 + t78 ^ 2) + (t180 * t127 + t11 + t7 + t8) * t187 + (t178 * t128 + t10 + t12 + t9 + (t184 * t127 + t187 * t128) * t187) * t184; m(7) * (t16 * t40 + t17 * t39) + m(5) * (t56 * t59 + t57 * t58) + m(6) * (t35 * t52 + t36 * t51) + (t184 * t189 + t187 * t190) * t186 + t216; m(5) * (t57 * t184 - t56 * t187) + m(6) * t202 - m(7) * t203; m(6) * (t18 * t15 + t35 * t61 + t36 * t60) + m(7) * (t14 * t13 + t16 * t54 + t17 * t53) + m(5) * (t50 * t41 + t56 * t94 + t57 * t93) + ((t9 / 0.2e1 + t12 / 0.2e1 + t10 / 0.2e1) * t187 + (-t11 / 0.2e1 - t8 / 0.2e1 - t7 / 0.2e1) * t184) * t186 + (t257 * t184 + t258 * t187) * t183 / 0.2e1 + t259 * t253 + t260 * t251; t216 * t183 + m(7) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(6) * (t18 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t50 ^ 2 + t56 ^ 2 + t57 ^ 2) + (t259 * t187 - t260 * t184 + (-t258 * t184 + t257 * t187) * t183) * t186; 0.2e1 * (-t262 / 0.2e1 - t263 / 0.2e1) * t186; t217 * t218 * t270; m(6) * (t183 * t15 - t186 * t198) + m(7) * (t183 * t13 - t186 * t199); m(7) * (t183 * t14 + t186 * t203) + m(6) * (t183 * t18 - t186 * t202); 0.2e1 * t217 * (t179 * t218 + t183 ^ 2); m(7) * (t121 * t39 - t123 * t40); m(7) * t194; m(7) * (t121 * t53 - t123 * t54 + t13 * t240); m(7) * (t121 * t17 - t123 * t16 + t14 * t240); m(7) * (t170 * t183 - t194) * t186; m(7) * (t179 * t170 ^ 2 + t121 ^ 2 + t123 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
