% Calculate joint inertia matrix for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:55
% EndTime: 2019-03-08 20:05:03
% DurationCPUTime: 4.49s
% Computational Cost: add. (21357->534), mult. (35424->777), div. (0->0), fcn. (44912->12), ass. (0->230)
t257 = rSges(7,3) + qJ(6);
t197 = sin(pkin(6));
t256 = t197 ^ 2;
t195 = sin(pkin(11));
t198 = cos(pkin(11));
t200 = cos(pkin(6));
t204 = sin(qJ(2));
t242 = t197 * t204;
t182 = -t195 * t242 + t198 * t200;
t245 = t195 * t200;
t183 = t198 * t242 + t245;
t206 = cos(qJ(2));
t241 = t197 * t206;
t146 = Icges(4,5) * t183 + Icges(4,6) * t182 - Icges(4,3) * t241;
t177 = Icges(3,6) * t200 + (Icges(3,4) * t204 + Icges(3,2) * t206) * t197;
t255 = -t177 + t146;
t254 = pkin(3) * t198;
t205 = cos(qJ(5));
t253 = pkin(5) * t205;
t196 = sin(pkin(10));
t199 = cos(pkin(10));
t240 = t200 * t204;
t185 = t196 * t206 + t199 * t240;
t232 = pkin(11) + qJ(4);
t194 = sin(t232);
t214 = cos(t232);
t243 = t197 * t199;
t167 = t185 * t214 - t194 * t243;
t239 = t200 * t206;
t184 = t196 * t204 - t199 * t239;
t203 = sin(qJ(5));
t138 = -t167 * t203 + t184 * t205;
t247 = t184 * t203;
t139 = t167 * t205 + t247;
t210 = t197 * t214;
t166 = t185 * t194 + t199 * t210;
t251 = rSges(7,1) * t139 + rSges(7,2) * t138 + pkin(5) * t247 + t257 * t166 + t167 * t253;
t187 = -t196 * t240 + t199 * t206;
t244 = t196 * t197;
t169 = t187 * t214 + t194 * t244;
t186 = t196 * t239 + t199 * t204;
t140 = -t169 * t203 + t186 * t205;
t246 = t186 * t203;
t141 = t169 * t205 + t246;
t168 = t187 * t194 - t196 * t210;
t250 = rSges(7,1) * t141 + rSges(7,2) * t140 + pkin(5) * t246 + t257 * t168 + t169 * t253;
t180 = t200 * t194 + t204 * t210;
t170 = -t180 * t203 - t205 * t241;
t218 = t203 * t241;
t171 = t180 * t205 - t218;
t179 = t194 * t242 - t200 * t214;
t249 = rSges(7,1) * t171 + rSges(7,2) * t170 - pkin(5) * t218 + t257 * t179 + t180 * t253;
t136 = pkin(4) * t169 + pkin(9) * t168;
t96 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t168;
t248 = -t136 - t96;
t220 = t195 * t244;
t119 = pkin(3) * t220 + pkin(8) * t186 + t187 * t254;
t165 = pkin(2) * t187 + qJ(3) * t186;
t163 = t200 * t165;
t237 = t200 * t119 + t163;
t219 = t195 * t243;
t118 = -pkin(3) * t219 + pkin(8) * t184 + t185 * t254;
t164 = pkin(2) * t185 + qJ(3) * t184;
t236 = -t118 - t164;
t135 = pkin(4) * t167 + pkin(9) * t166;
t156 = t180 * pkin(4) + t179 * pkin(9);
t235 = t135 * t241 + t184 * t156;
t188 = (pkin(2) * t204 - qJ(3) * t206) * t197;
t234 = -pkin(3) * t245 - (-pkin(8) * t206 + t204 * t254) * t197 - t188;
t233 = t164 * t244 + t165 * t243;
t81 = Icges(7,5) * t139 + Icges(7,6) * t138 + Icges(7,3) * t166;
t85 = Icges(7,4) * t139 + Icges(7,2) * t138 + Icges(7,6) * t166;
t89 = Icges(7,1) * t139 + Icges(7,4) * t138 + Icges(7,5) * t166;
t35 = t138 * t85 + t139 * t89 + t166 * t81;
t82 = Icges(7,5) * t141 + Icges(7,6) * t140 + Icges(7,3) * t168;
t86 = Icges(7,4) * t141 + Icges(7,2) * t140 + Icges(7,6) * t168;
t90 = Icges(7,1) * t141 + Icges(7,4) * t140 + Icges(7,5) * t168;
t36 = t138 * t86 + t139 * t90 + t166 * t82;
t101 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t179;
t103 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t179;
t105 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t179;
t51 = t101 * t166 + t103 * t138 + t105 * t139;
t1 = t166 * t35 + t168 * t36 + t179 * t51;
t83 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t166;
t87 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t166;
t91 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t166;
t37 = t138 * t87 + t139 * t91 + t166 * t83;
t84 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t168;
t88 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t168;
t92 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t168;
t38 = t138 * t88 + t139 * t92 + t166 * t84;
t102 = Icges(6,5) * t171 + Icges(6,6) * t170 + Icges(6,3) * t179;
t104 = Icges(6,4) * t171 + Icges(6,2) * t170 + Icges(6,6) * t179;
t106 = Icges(6,1) * t171 + Icges(6,4) * t170 + Icges(6,5) * t179;
t52 = t102 * t166 + t104 * t138 + t106 * t139;
t2 = t166 * t37 + t168 * t38 + t179 * t52;
t231 = t2 / 0.2e1 + t1 / 0.2e1;
t39 = t140 * t85 + t141 * t89 + t168 * t81;
t40 = t140 * t86 + t141 * t90 + t168 * t82;
t53 = t101 * t168 + t103 * t140 + t105 * t141;
t3 = t166 * t39 + t168 * t40 + t179 * t53;
t41 = t140 * t87 + t141 * t91 + t168 * t83;
t42 = t140 * t88 + t141 * t92 + t168 * t84;
t54 = t102 * t168 + t104 * t140 + t106 * t141;
t4 = t166 * t41 + t168 * t42 + t179 * t54;
t230 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t35 * t184 + t36 * t186 - t241 * t51;
t6 = t37 * t184 + t38 * t186 - t241 * t52;
t229 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t39 * t184 + t40 * t186 - t241 * t53;
t8 = t41 * t184 + t42 * t186 - t241 * t54;
t228 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t200 * t52 + (t196 * t38 - t199 * t37) * t197;
t9 = t200 * t51 + (t196 * t36 - t199 * t35) * t197;
t227 = t10 / 0.2e1 + t9 / 0.2e1;
t226 = -t136 - t250;
t11 = t200 * t53 + (t196 * t40 - t199 * t39) * t197;
t12 = t200 * t54 + (t196 * t42 - t199 * t41) * t197;
t225 = t11 / 0.2e1 + t12 / 0.2e1;
t45 = t170 * t85 + t171 * t89 + t179 * t81;
t46 = t170 * t86 + t171 * t90 + t179 * t82;
t58 = t101 * t179 + t103 * t170 + t105 * t171;
t13 = t166 * t45 + t168 * t46 + t179 * t58;
t47 = t170 * t87 + t171 * t91 + t179 * t83;
t48 = t170 * t88 + t171 * t92 + t179 * t84;
t59 = t102 * t179 + t104 * t170 + t106 * t171;
t14 = t166 * t47 + t168 * t48 + t179 * t59;
t224 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t45 * t184 + t46 * t186 - t241 * t58;
t16 = t47 * t184 + t48 * t186 - t241 * t59;
t223 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t200 * t58 + (t196 * t46 - t199 * t45) * t197;
t18 = t200 * t59 + (t196 * t48 - t199 * t47) * t197;
t222 = t18 / 0.2e1 + t17 / 0.2e1;
t221 = -m(4) - m(5) - m(6) - m(7);
t217 = t200 * t136 + t237;
t216 = -t135 + t236;
t215 = -t156 + t234;
t213 = (-t183 * rSges(4,1) - t182 * rSges(4,2) + rSges(4,3) * t241 - t188) * t197;
t212 = t118 * t244 + t119 * t243 + t233;
t145 = t180 * rSges(5,1) - t179 * rSges(5,2) - rSges(5,3) * t241;
t211 = (-t145 + t234) * t197;
t109 = rSges(6,1) * t171 + rSges(6,2) * t170 + rSges(6,3) * t179;
t209 = (-t109 + t215) * t197;
t208 = t135 * t244 + t136 * t243 + t212;
t207 = (t215 - t249) * t197;
t181 = t200 * rSges(3,3) + (rSges(3,1) * t204 + rSges(3,2) * t206) * t197;
t178 = Icges(3,5) * t200 + (Icges(3,1) * t204 + Icges(3,4) * t206) * t197;
t176 = Icges(3,3) * t200 + (Icges(3,5) * t204 + Icges(3,6) * t206) * t197;
t175 = t187 * t198 + t220;
t174 = -t187 * t195 + t198 * t244;
t173 = t185 * t198 - t219;
t172 = -t185 * t195 - t198 * t243;
t159 = rSges(3,1) * t187 - rSges(3,2) * t186 + rSges(3,3) * t244;
t158 = rSges(3,1) * t185 - rSges(3,2) * t184 - rSges(3,3) * t243;
t154 = Icges(3,1) * t187 - Icges(3,4) * t186 + Icges(3,5) * t244;
t153 = Icges(3,1) * t185 - Icges(3,4) * t184 - Icges(3,5) * t243;
t152 = Icges(3,4) * t187 - Icges(3,2) * t186 + Icges(3,6) * t244;
t151 = Icges(3,4) * t185 - Icges(3,2) * t184 - Icges(3,6) * t243;
t150 = Icges(3,5) * t187 - Icges(3,6) * t186 + Icges(3,3) * t244;
t149 = Icges(3,5) * t185 - Icges(3,6) * t184 - Icges(3,3) * t243;
t148 = Icges(4,1) * t183 + Icges(4,4) * t182 - Icges(4,5) * t241;
t147 = Icges(4,4) * t183 + Icges(4,2) * t182 - Icges(4,6) * t241;
t144 = Icges(5,1) * t180 - Icges(5,4) * t179 - Icges(5,5) * t241;
t143 = Icges(5,4) * t180 - Icges(5,2) * t179 - Icges(5,6) * t241;
t142 = Icges(5,5) * t180 - Icges(5,6) * t179 - Icges(5,3) * t241;
t134 = -t158 * t200 - t181 * t243;
t133 = t159 * t200 - t181 * t244;
t128 = rSges(4,1) * t175 + rSges(4,2) * t174 + rSges(4,3) * t186;
t127 = rSges(4,1) * t173 + rSges(4,2) * t172 + rSges(4,3) * t184;
t126 = Icges(4,1) * t175 + Icges(4,4) * t174 + Icges(4,5) * t186;
t125 = Icges(4,1) * t173 + Icges(4,4) * t172 + Icges(4,5) * t184;
t124 = Icges(4,4) * t175 + Icges(4,2) * t174 + Icges(4,6) * t186;
t123 = Icges(4,4) * t173 + Icges(4,2) * t172 + Icges(4,6) * t184;
t122 = Icges(4,5) * t175 + Icges(4,6) * t174 + Icges(4,3) * t186;
t121 = Icges(4,5) * t173 + Icges(4,6) * t172 + Icges(4,3) * t184;
t120 = t186 * t135;
t117 = rSges(5,1) * t169 - rSges(5,2) * t168 + rSges(5,3) * t186;
t116 = rSges(5,1) * t167 - rSges(5,2) * t166 + rSges(5,3) * t184;
t115 = Icges(5,1) * t169 - Icges(5,4) * t168 + Icges(5,5) * t186;
t114 = Icges(5,1) * t167 - Icges(5,4) * t166 + Icges(5,5) * t184;
t113 = Icges(5,4) * t169 - Icges(5,2) * t168 + Icges(5,6) * t186;
t112 = Icges(5,4) * t167 - Icges(5,2) * t166 + Icges(5,6) * t184;
t111 = Icges(5,5) * t169 - Icges(5,6) * t168 + Icges(5,3) * t186;
t110 = Icges(5,5) * t167 - Icges(5,6) * t166 + Icges(5,3) * t184;
t98 = (t158 * t196 + t159 * t199) * t197;
t94 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t166;
t80 = -t117 * t241 - t186 * t145;
t79 = t116 * t241 + t184 * t145;
t76 = (-t127 - t164) * t200 + t199 * t213;
t75 = t128 * t200 + t196 * t213 + t163;
t74 = -t142 * t241 - t179 * t143 + t180 * t144;
t73 = t116 * t186 - t117 * t184;
t72 = t142 * t186 - t143 * t168 + t144 * t169;
t71 = t142 * t184 - t143 * t166 + t144 * t167;
t70 = (t127 * t196 + t128 * t199) * t197 + t233;
t69 = -t109 * t168 + t179 * t96;
t68 = t109 * t166 - t179 * t94;
t67 = -t111 * t241 - t179 * t113 + t180 * t115;
t66 = -t110 * t241 - t179 * t112 + t180 * t114;
t65 = t111 * t186 - t113 * t168 + t115 * t169;
t64 = t110 * t186 - t112 * t168 + t114 * t169;
t63 = t111 * t184 - t113 * t166 + t115 * t167;
t62 = t110 * t184 - t112 * t166 + t114 * t167;
t61 = (-t116 + t236) * t200 + t199 * t211;
t60 = t117 * t200 + t196 * t211 + t237;
t57 = -t166 * t96 + t168 * t94;
t56 = t248 * t241 + (-t109 - t156) * t186;
t55 = t184 * t109 + t241 * t94 + t235;
t50 = (t116 * t196 + t117 * t199) * t197 + t212;
t49 = t184 * t248 + t186 * t94 + t120;
t44 = (-t94 + t216) * t200 + t199 * t209;
t43 = t196 * t209 + t200 * t96 + t217;
t34 = -t168 * t249 + t179 * t250;
t33 = t166 * t249 - t179 * t251;
t32 = t226 * t241 + (-t156 - t249) * t186;
t31 = t184 * t249 + t241 * t251 + t235;
t30 = (t196 * t94 + t199 * t96) * t197 + t208;
t29 = -t166 * t250 + t168 * t251;
t28 = t200 * t74 + (t196 * t67 - t199 * t66) * t197;
t27 = (t216 - t251) * t200 + t199 * t207;
t26 = t196 * t207 + t200 * t250 + t217;
t25 = t66 * t184 + t67 * t186 - t241 * t74;
t24 = t184 * t226 + t186 * t251 + t120;
t23 = t200 * t72 + (t196 * t65 - t199 * t64) * t197;
t22 = t200 * t71 + (t196 * t63 - t199 * t62) * t197;
t21 = t64 * t184 + t65 * t186 - t241 * t72;
t20 = t62 * t184 + t63 * t186 - t241 * t71;
t19 = (t196 * t251 + t199 * t250) * t197 + t208;
t77 = [m(2) + m(3) - t221; m(3) * t98 + m(4) * t70 + m(5) * t50 + m(6) * t30 + m(7) * t19; m(7) * (t19 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t30 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t50 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t70 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(3) * (t133 ^ 2 + t134 ^ 2 + t98 ^ 2) + (t11 + t12 + t23 + ((t122 * t186 + t174 * t124 + t175 * t126) * t196 - (t186 * t121 + t174 * t123 + t175 * t125) * t199) * t197 + (t150 * t244 - t152 * t186 + t187 * t154) * t244) * t244 + (-t10 - t9 - t22 - ((t184 * t122 + t172 * t124 + t173 * t126) * t196 - (t121 * t184 + t123 * t172 + t125 * t173) * t199) * t197 + (-t149 * t243 - t151 * t184 + t153 * t185) * t243 + (-t149 * t244 + t150 * t243 + t186 * t151 + t184 * t152 - t187 * t153 - t185 * t154) * t244) * t243 + (t17 + t18 + t28 + ((t152 * t206 + t154 * t204) * t196 - (t151 * t206 + t153 * t204) * t199) * t256 + ((-t149 * t199 + t150 * t196 + t177 * t206 + t178 * t204) * t197 - t146 * t241 + t182 * t147 + t183 * t148 + t200 * t176) * t200 + (-t122 * t241 + t182 * t124 + t183 * t126 + t147 * t174 + t148 * t175 + t176 * t244 + t178 * t187 + t255 * t186) * t244 + (t121 * t241 - t182 * t123 - t183 * t125 - t147 * t172 - t148 * t173 + t176 * t243 - t178 * t185 - t255 * t184) * t243) * t200; t221 * t241; m(7) * (t184 * t26 + t186 * t27 - t19 * t241) + m(6) * (t184 * t43 + t186 * t44 - t241 * t30) + m(5) * (t184 * t60 + t186 * t61 - t241 * t50) + m(4) * (t184 * t75 + t186 * t76 - t241 * t70); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t256 * t206 ^ 2 + t184 ^ 2 + t186 ^ 2); m(5) * t73 + m(6) * t49 + m(7) * t24; (t25 / 0.2e1 + t223) * t200 + (t23 / 0.2e1 + t225) * t186 + (t22 / 0.2e1 + t227) * t184 + m(7) * (t19 * t24 + t26 * t32 + t27 * t31) + m(6) * (t30 * t49 + t43 * t56 + t44 * t55) + m(5) * (t50 * t73 + t60 * t80 + t61 * t79) + ((-t28 / 0.2e1 - t222) * t206 + (-t20 / 0.2e1 - t229) * t199 + (t21 / 0.2e1 + t228) * t196) * t197; m(5) * (t80 * t184 + t79 * t186 - t241 * t73) + m(6) * (t56 * t184 + t55 * t186 - t241 * t49) + m(7) * (t32 * t184 + t31 * t186 - t24 * t241); (-t15 - t16 - t25) * t241 + (t8 + t7 + t21) * t186 + (t6 + t5 + t20) * t184 + m(7) * (t24 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t49 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t73 ^ 2 + t79 ^ 2 + t80 ^ 2); m(6) * t57 + m(7) * t29; t224 * t200 + t222 * t179 + t225 * t168 + t227 * t166 + m(7) * (t19 * t29 + t26 * t34 + t27 * t33) + m(6) * (t30 * t57 + t43 * t69 + t44 * t68) + (t196 * t230 - t199 * t231) * t197; m(6) * (t69 * t184 + t68 * t186 - t241 * t57) + m(7) * (t34 * t184 + t33 * t186 - t241 * t29); -t224 * t241 + t230 * t186 + t231 * t184 + t223 * t179 + t228 * t168 + t229 * t166 + m(7) * (t24 * t29 + t31 * t33 + t32 * t34) + m(6) * (t49 * t57 + t55 * t68 + t56 * t69); (t13 + t14) * t179 + (t3 + t4) * t168 + (t2 + t1) * t166 + m(7) * (t29 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t57 ^ 2 + t68 ^ 2 + t69 ^ 2); m(7) * t179; m(7) * (t166 * t26 + t168 * t27 + t179 * t19); m(7) * (t166 * t184 + t168 * t186 - t179 * t241); m(7) * (t166 * t32 + t168 * t31 + t179 * t24); m(7) * (t166 * t34 + t168 * t33 + t179 * t29); m(7) * (t166 ^ 2 + t168 ^ 2 + t179 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t77(1) t77(2) t77(4) t77(7) t77(11) t77(16); t77(2) t77(3) t77(5) t77(8) t77(12) t77(17); t77(4) t77(5) t77(6) t77(9) t77(13) t77(18); t77(7) t77(8) t77(9) t77(10) t77(14) t77(19); t77(11) t77(12) t77(13) t77(14) t77(15) t77(20); t77(16) t77(17) t77(18) t77(19) t77(20) t77(21);];
Mq  = res;
