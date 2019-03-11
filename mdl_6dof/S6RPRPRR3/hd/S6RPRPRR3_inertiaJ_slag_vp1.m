% Calculate joint inertia matrix for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:39
% EndTime: 2019-03-09 03:40:44
% DurationCPUTime: 2.62s
% Computational Cost: add. (9761->411), mult. (7578->595), div. (0->0), fcn. (8094->12), ass. (0->201)
t180 = qJ(1) + pkin(10);
t172 = sin(t180);
t251 = t172 / 0.2e1;
t174 = cos(t180);
t249 = t174 / 0.2e1;
t183 = cos(pkin(11));
t168 = t183 * pkin(4) + pkin(3);
t179 = pkin(11) + qJ(5);
t173 = cos(t179);
t143 = pkin(5) * t173 + t168;
t171 = sin(t179);
t182 = sin(pkin(11));
t144 = pkin(4) * t182 + pkin(5) * t171;
t187 = cos(qJ(3));
t228 = t174 * t187;
t175 = qJ(6) + t179;
t166 = sin(t175);
t167 = cos(t175);
t102 = -t166 * t228 + t167 * t172;
t103 = t166 * t172 + t167 * t228;
t185 = sin(qJ(3));
t229 = t174 * t185;
t72 = t103 * rSges(7,1) + t102 * rSges(7,2) + rSges(7,3) * t229;
t256 = t143 * t228 + t172 * t144 + t72;
t169 = t172 ^ 2;
t170 = t174 ^ 2;
t255 = m(5) / 0.2e1;
t254 = m(6) / 0.2e1;
t253 = m(7) / 0.2e1;
t233 = t172 * t185;
t232 = t172 * t187;
t100 = -t166 * t232 - t167 * t174;
t101 = -t166 * t174 + t167 * t232;
t65 = Icges(7,5) * t101 + Icges(7,6) * t100 + Icges(7,3) * t233;
t67 = Icges(7,4) * t101 + Icges(7,2) * t100 + Icges(7,6) * t233;
t69 = Icges(7,1) * t101 + Icges(7,4) * t100 + Icges(7,5) * t233;
t20 = t100 * t67 + t101 * t69 + t65 * t233;
t66 = Icges(7,5) * t103 + Icges(7,6) * t102 + Icges(7,3) * t229;
t68 = Icges(7,4) * t103 + Icges(7,2) * t102 + Icges(7,6) * t229;
t70 = Icges(7,1) * t103 + Icges(7,4) * t102 + Icges(7,5) * t229;
t21 = t100 * t68 + t101 * t70 + t66 * t233;
t104 = -Icges(7,3) * t187 + (Icges(7,5) * t167 - Icges(7,6) * t166) * t185;
t105 = -Icges(7,6) * t187 + (Icges(7,4) * t167 - Icges(7,2) * t166) * t185;
t106 = -Icges(7,5) * t187 + (Icges(7,1) * t167 - Icges(7,4) * t166) * t185;
t40 = t100 * t105 + t101 * t106 + t104 * t233;
t5 = -t40 * t187 + (t172 * t20 + t174 * t21) * t185;
t22 = t102 * t67 + t103 * t69 + t65 * t229;
t23 = t102 * t68 + t103 * t70 + t66 * t229;
t41 = t102 * t105 + t103 * t106 + t104 * t229;
t6 = -t41 * t187 + (t172 * t22 + t174 * t23) * t185;
t252 = t6 * t229 + t5 * t233;
t250 = -t174 / 0.2e1;
t248 = -t187 / 0.2e1;
t151 = rSges(4,1) * t185 + rSges(4,2) * t187;
t247 = m(4) * t151;
t246 = pkin(3) * t187;
t186 = sin(qJ(1));
t245 = t186 * pkin(1);
t244 = -pkin(3) + t168;
t184 = -pkin(8) - qJ(4);
t178 = -pkin(9) + t184;
t214 = t178 - t184;
t234 = t172 * t182;
t219 = -pkin(4) * t234 - t168 * t228;
t243 = -t214 * t229 + t219 + t256;
t113 = -t187 * rSges(7,3) + (rSges(7,1) * t167 - rSges(7,2) * t166) * t185;
t196 = -t101 * rSges(7,1) - t100 * rSges(7,2);
t71 = rSges(7,3) * t233 - t196;
t51 = t113 * t233 + t187 * t71;
t242 = t174 * rSges(4,3);
t236 = t105 * t166;
t95 = t185 * t167 * t106;
t53 = -t187 * t104 - t185 * t236 + t95;
t241 = t53 * t187;
t218 = t143 - t168;
t94 = t218 * t185 + t214 * t187;
t240 = -t113 - t94;
t239 = Icges(4,4) * t185;
t238 = Icges(4,4) * t187;
t237 = qJ(4) * t185;
t124 = -Icges(6,6) * t187 + (Icges(6,4) * t173 - Icges(6,2) * t171) * t185;
t235 = t124 * t171;
t231 = t174 * t144;
t230 = t174 * t182;
t227 = t178 * t185;
t226 = t182 * t187;
t225 = t183 * t187;
t224 = t184 * t185;
t150 = t185 * pkin(3) - t187 * qJ(4);
t223 = -(qJ(4) + t184) * t187 - t244 * t185 - t150;
t216 = pkin(3) * t228 + qJ(4) * t229;
t222 = t169 * (t237 + t246) + t174 * t216;
t220 = t187 * rSges(5,3) - (rSges(5,1) * t183 - rSges(5,2) * t182) * t185 - t150;
t217 = pkin(4) * t230 + t172 * t224;
t215 = t169 + t170;
t213 = -m(5) - m(6) - m(7);
t119 = -t171 * t232 - t173 * t174;
t120 = -t171 * t174 + t173 * t232;
t74 = Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t233;
t76 = Icges(6,4) * t120 + Icges(6,2) * t119 + Icges(6,6) * t233;
t78 = Icges(6,1) * t120 + Icges(6,4) * t119 + Icges(6,5) * t233;
t34 = -t187 * t74 + (-t171 * t76 + t173 * t78) * t185;
t123 = -Icges(6,3) * t187 + (Icges(6,5) * t173 - Icges(6,6) * t171) * t185;
t125 = -Icges(6,5) * t187 + (Icges(6,1) * t173 - Icges(6,4) * t171) * t185;
t45 = t119 * t124 + t120 * t125 + t123 * t233;
t212 = t45 / 0.2e1 + t34 / 0.2e1;
t121 = -t171 * t228 + t172 * t173;
t122 = t171 * t172 + t173 * t228;
t75 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t229;
t77 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t229;
t79 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t229;
t35 = -t187 * t75 + (-t171 * t77 + t173 * t79) * t185;
t46 = t121 * t124 + t122 * t125 + t123 * t229;
t211 = t46 / 0.2e1 + t35 / 0.2e1;
t81 = t122 * rSges(6,1) + t121 * rSges(6,2) + rSges(6,3) * t229;
t126 = -t187 * rSges(6,3) + (rSges(6,1) * t173 - rSges(6,2) * t171) * t185;
t210 = -t126 + t223;
t133 = t172 * t183 - t174 * t226;
t134 = t174 * t225 + t234;
t209 = t134 * rSges(5,1) + t133 * rSges(5,2) + rSges(5,3) * t229;
t188 = cos(qJ(1));
t177 = t188 * pkin(1);
t208 = t174 * pkin(2) + t172 * pkin(7) + t177;
t207 = t233 / 0.2e1;
t206 = t229 / 0.2e1;
t205 = t174 * pkin(7) - t245;
t32 = -t187 * t65 + (-t166 * t67 + t167 * t69) * t185;
t33 = -t187 * t66 + (-t166 * t68 + t167 * t70) * t185;
t204 = (t32 + t40) * t207 + (t33 + t41) * t206;
t189 = -t174 * t224 - t219;
t203 = t172 * ((t244 * t187 - t237) * t172 - t217) + t174 * (t189 - t216) + t222;
t9 = -t241 + (t172 * t32 + t174 * t33) * t185;
t202 = -t187 * t9 + t252;
t201 = t223 + t240;
t12 = t172 * t21 - t174 * t20;
t13 = t172 * t23 - t174 * t22;
t200 = t12 * t207 + t13 * t206 + t5 * t250 + t6 * t251 + (t33 * t172 - t32 * t174) * t248;
t199 = rSges(4,1) * t187 - rSges(4,2) * t185;
t131 = -t172 * t226 - t174 * t183;
t132 = t172 * t225 - t230;
t198 = -t132 * rSges(5,1) - t131 * rSges(5,2);
t197 = -t120 * rSges(6,1) - t119 * rSges(6,2);
t195 = Icges(4,1) * t187 - t239;
t194 = -Icges(4,2) * t185 + t238;
t193 = Icges(4,5) * t187 - Icges(4,6) * t185;
t190 = rSges(4,1) * t228 - rSges(4,2) * t229 + t172 * rSges(4,3);
t153 = rSges(2,1) * t188 - t186 * rSges(2,2);
t152 = -t186 * rSges(2,1) - rSges(2,2) * t188;
t147 = Icges(4,5) * t185 + Icges(4,6) * t187;
t140 = rSges(3,1) * t174 - rSges(3,2) * t172 + t177;
t139 = -rSges(3,1) * t172 - rSges(3,2) * t174 - t245;
t137 = -Icges(5,5) * t187 + (Icges(5,1) * t183 - Icges(5,4) * t182) * t185;
t136 = -Icges(5,6) * t187 + (Icges(5,4) * t183 - Icges(5,2) * t182) * t185;
t108 = Icges(4,3) * t172 + t193 * t174;
t107 = -Icges(4,3) * t174 + t193 * t172;
t97 = t185 * t173 * t125;
t93 = t220 * t174;
t92 = t220 * t172;
t91 = t190 + t208;
t90 = t242 + (-pkin(2) - t199) * t172 + t205;
t89 = Icges(5,1) * t134 + Icges(5,4) * t133 + Icges(5,5) * t229;
t88 = Icges(5,1) * t132 + Icges(5,4) * t131 + Icges(5,5) * t233;
t87 = Icges(5,4) * t134 + Icges(5,2) * t133 + Icges(5,6) * t229;
t86 = Icges(5,4) * t132 + Icges(5,2) * t131 + Icges(5,6) * t233;
t85 = Icges(5,5) * t134 + Icges(5,6) * t133 + Icges(5,3) * t229;
t84 = Icges(5,5) * t132 + Icges(5,6) * t131 + Icges(5,3) * t233;
t80 = rSges(6,3) * t233 - t197;
t73 = t174 * t190 + (t199 * t172 - t242) * t172;
t63 = t71 * t229;
t62 = t210 * t174;
t61 = t210 * t172;
t59 = -t231 + (t218 * t187 - t227) * t172 + t217;
t58 = t208 + t209 + t216;
t57 = (-t246 - pkin(2) + (-rSges(5,3) - qJ(4)) * t185) * t172 + t198 + t205;
t56 = -t187 * t123 - t185 * t235 + t97;
t55 = -t126 * t229 - t187 * t81;
t54 = t126 * t233 + t187 * t80;
t52 = -t113 * t229 - t187 * t72;
t50 = t189 + t208 + t81;
t49 = (-rSges(6,3) * t185 - t168 * t187 - pkin(2)) * t172 + t197 + t205 + t217;
t48 = t201 * t174;
t47 = t201 * t172;
t44 = -t174 * t227 + t208 + t256;
t43 = t231 + (-t143 * t187 - pkin(2) + (-rSges(7,3) + t178) * t185) * t172 + t196 + t205;
t42 = (-t172 * t81 + t174 * t80) * t185;
t39 = -t72 * t233 + t63;
t36 = t172 * (rSges(5,3) * t233 - t198) + t174 * t209 + t222;
t29 = -t243 * t187 + t240 * t229;
t28 = t187 * t59 + t94 * t233 + t51;
t27 = t121 * t77 + t122 * t79 + t75 * t229;
t26 = t121 * t76 + t122 * t78 + t74 * t229;
t25 = t119 * t77 + t120 * t79 + t75 * t233;
t24 = t119 * t76 + t120 * t78 + t74 * t233;
t19 = t172 * t80 + t174 * t81 + t203;
t18 = t63 + (-t243 * t172 + t174 * t59) * t185;
t16 = t243 * t174 + (t59 + t71) * t172 + t203;
t15 = t172 * t27 - t174 * t26;
t14 = t172 * t25 - t174 * t24;
t8 = -t46 * t187 + (t172 * t26 + t174 * t27) * t185;
t7 = -t45 * t187 + (t172 * t24 + t174 * t25) * t185;
t1 = [Icges(2,3) + Icges(3,3) + t95 + t97 + (-t104 - t123 + t239 - (Icges(5,5) * t183 - Icges(5,6) * t182) * t185 + (Icges(4,2) + Icges(5,3)) * t187) * t187 + (Icges(4,1) * t185 - t136 * t182 + t137 * t183 - t235 - t236 + t238) * t185 + m(7) * (t43 ^ 2 + t44 ^ 2) + m(6) * (t49 ^ 2 + t50 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t90 ^ 2 + t91 ^ 2) + m(2) * (t152 ^ 2 + t153 ^ 2) + m(3) * (t139 ^ 2 + t140 ^ 2); 0; m(3) + m(4) - t213; m(7) * (t43 * t48 + t44 * t47) + m(6) * (t49 * t62 + t50 * t61) + m(5) * (t57 * t93 + t58 * t92) + (-t40 / 0.2e1 - t32 / 0.2e1 - t90 * t247 - t131 * t136 / 0.2e1 - t132 * t137 / 0.2e1 + t147 * t249 + (Icges(4,6) * t249 - t194 * t172 / 0.2e1 + t84 / 0.2e1) * t187 - t212) * t174 + (t33 / 0.2e1 + t41 / 0.2e1 - t91 * t247 + t133 * t136 / 0.2e1 + t134 * t137 / 0.2e1 + t147 * t251 + (Icges(4,6) * t251 + t194 * t249 - t85 / 0.2e1) * t187 + t211) * t172 + ((Icges(4,5) * t172 + t195 * t174 - t182 * t87 + t183 * t89) * t251 + (-Icges(4,5) * t174 + t195 * t172 - t182 * t86 + t183 * t88) * t250) * t185; m(4) * t73 + m(5) * t36 + m(6) * t19 + m(7) * t16; m(7) * (t16 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t19 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t36 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(4) * (t215 * t151 ^ 2 + t73 ^ 2) + (-t170 * t107 - t12 - t14 + (t131 * t86 + t132 * t88 + t84 * t233) * t174) * t174 + (t13 + t15 + t169 * t108 + (t133 * t87 + t134 * t89 + t85 * t229) * t172 + (-t172 * t107 + t174 * t108 - t131 * t87 - t132 * t89 - t133 * t86 - t134 * t88 - t84 * t229 - t85 * t233) * t174) * t172; 0.2e1 * ((t172 * t44 + t174 * t43) * t253 + (t172 * t50 + t174 * t49) * t254 + (t172 * t58 + t174 * t57) * t255) * t185; t213 * t187; m(7) * (-t187 * t16 + (t172 * t47 + t174 * t48) * t185) + m(6) * (-t187 * t19 + (t172 * t61 + t174 * t62) * t185) + m(5) * (-t187 * t36 + (t172 * t92 + t174 * t93) * t185); 0.2e1 * (t255 + t254 + t253) * (t215 * t185 ^ 2 + t187 ^ 2); (-t53 - t56) * t187 + m(7) * (t28 * t43 + t29 * t44) + m(6) * (t49 * t54 + t50 * t55) + (t212 * t172 + t211 * t174) * t185 + t204; m(6) * t42 + m(7) * t18; t8 * t251 + t7 * t250 + (t35 * t172 - t34 * t174) * t248 + (t14 * t251 + t15 * t249) * t185 + m(7) * (t16 * t18 + t28 * t48 + t29 * t47) + m(6) * (t19 * t42 + t54 * t62 + t55 * t61) + t200; m(6) * (-t42 * t187 + (t172 * t55 + t174 * t54) * t185) + m(7) * (-t18 * t187 + (t172 * t29 + t174 * t28) * t185); (t56 * t187 - t9) * t187 + m(7) * (t18 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t42 ^ 2 + t54 ^ 2 + t55 ^ 2) + (t174 * t8 + t172 * t7 - t187 * (t172 * t34 + t174 * t35)) * t185 + t252; -t241 + m(7) * (t43 * t51 + t44 * t52) + t204; m(7) * t39; m(7) * (t16 * t39 + t47 * t52 + t48 * t51) + t200; m(7) * (-t39 * t187 + (t172 * t52 + t174 * t51) * t185); m(7) * (t18 * t39 + t28 * t51 + t29 * t52) + t202; m(7) * (t39 ^ 2 + t51 ^ 2 + t52 ^ 2) + t202;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
