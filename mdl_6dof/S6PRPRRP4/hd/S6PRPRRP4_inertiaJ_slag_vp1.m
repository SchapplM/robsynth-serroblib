% Calculate joint inertia matrix for
% S6PRPRRP4
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:18
% EndTime: 2019-03-08 20:09:28
% DurationCPUTime: 4.42s
% Computational Cost: add. (20453->529), mult. (34652->770), div. (0->0), fcn. (44100->12), ass. (0->227)
t253 = rSges(7,1) + pkin(5);
t252 = rSges(7,3) + qJ(6);
t198 = sin(pkin(6));
t251 = t198 ^ 2;
t196 = sin(pkin(11));
t199 = cos(pkin(11));
t201 = cos(pkin(6));
t204 = sin(qJ(2));
t241 = t198 * t204;
t184 = -t196 * t241 + t199 * t201;
t244 = t196 * t201;
t185 = t199 * t241 + t244;
t205 = cos(qJ(2));
t240 = t198 * t205;
t146 = Icges(4,5) * t185 + Icges(4,6) * t184 - Icges(4,3) * t240;
t177 = Icges(3,6) * t201 + (Icges(3,4) * t204 + Icges(3,2) * t205) * t198;
t250 = -t177 + t146;
t249 = cos(qJ(5));
t248 = pkin(3) * t199;
t197 = sin(pkin(10));
t200 = cos(pkin(10));
t239 = t201 * t204;
t187 = t197 * t205 + t200 * t239;
t230 = pkin(11) + qJ(4);
t195 = sin(t230);
t213 = cos(t230);
t242 = t198 * t200;
t167 = t187 * t213 - t195 * t242;
t238 = t201 * t205;
t186 = t197 * t204 - t200 * t238;
t203 = sin(qJ(5));
t138 = t167 * t203 - t186 * t249;
t139 = t167 * t249 + t186 * t203;
t209 = t198 * t213;
t166 = t187 * t195 + t200 * t209;
t247 = rSges(7,2) * t166 + t252 * t138 + t253 * t139;
t189 = -t197 * t239 + t200 * t205;
t243 = t197 * t198;
t169 = t189 * t213 + t195 * t243;
t188 = t197 * t238 + t200 * t204;
t140 = t169 * t203 - t188 * t249;
t141 = t169 * t249 + t188 * t203;
t168 = t189 * t195 - t197 * t209;
t246 = rSges(7,2) * t168 + t252 * t140 + t253 * t141;
t135 = pkin(4) * t169 + pkin(9) * t168;
t94 = rSges(6,1) * t141 - rSges(6,2) * t140 + rSges(6,3) * t168;
t245 = -t135 - t94;
t218 = t196 * t243;
t118 = pkin(3) * t218 + pkin(8) * t188 + t189 * t248;
t165 = pkin(2) * t189 + qJ(3) * t188;
t163 = t201 * t165;
t236 = t201 * t118 + t163;
t180 = t201 * t195 + t204 * t209;
t170 = t180 * t203 + t240 * t249;
t171 = t180 * t249 - t203 * t240;
t179 = t195 * t241 - t201 * t213;
t235 = rSges(7,2) * t179 + t252 * t170 + t253 * t171;
t217 = t196 * t242;
t117 = -pkin(3) * t217 + pkin(8) * t186 + t187 * t248;
t164 = pkin(2) * t187 + qJ(3) * t186;
t234 = -t117 - t164;
t134 = pkin(4) * t167 + pkin(9) * t166;
t156 = pkin(4) * t180 + pkin(9) * t179;
t233 = t134 * t240 + t186 * t156;
t190 = (pkin(2) * t204 - qJ(3) * t205) * t198;
t232 = -pkin(3) * t244 - (-pkin(8) * t205 + t204 * t248) * t198 - t190;
t231 = t164 * t243 + t165 * t242;
t79 = Icges(7,5) * t139 + Icges(7,6) * t166 + Icges(7,3) * t138;
t83 = Icges(7,4) * t139 + Icges(7,2) * t166 + Icges(7,6) * t138;
t87 = Icges(7,1) * t139 + Icges(7,4) * t166 + Icges(7,5) * t138;
t33 = t138 * t79 + t139 * t87 + t166 * t83;
t80 = Icges(7,5) * t141 + Icges(7,6) * t168 + Icges(7,3) * t140;
t84 = Icges(7,4) * t141 + Icges(7,2) * t168 + Icges(7,6) * t140;
t88 = Icges(7,1) * t141 + Icges(7,4) * t168 + Icges(7,5) * t140;
t34 = t138 * t80 + t139 * t88 + t166 * t84;
t100 = Icges(7,5) * t171 + Icges(7,6) * t179 + Icges(7,3) * t170;
t102 = Icges(7,4) * t171 + Icges(7,2) * t179 + Icges(7,6) * t170;
t104 = Icges(7,1) * t171 + Icges(7,4) * t179 + Icges(7,5) * t170;
t51 = t100 * t138 + t102 * t166 + t104 * t139;
t1 = t166 * t33 + t168 * t34 + t179 * t51;
t81 = Icges(6,5) * t139 - Icges(6,6) * t138 + Icges(6,3) * t166;
t85 = Icges(6,4) * t139 - Icges(6,2) * t138 + Icges(6,6) * t166;
t89 = Icges(6,1) * t139 - Icges(6,4) * t138 + Icges(6,5) * t166;
t35 = -t138 * t85 + t139 * t89 + t166 * t81;
t82 = Icges(6,5) * t141 - Icges(6,6) * t140 + Icges(6,3) * t168;
t86 = Icges(6,4) * t141 - Icges(6,2) * t140 + Icges(6,6) * t168;
t90 = Icges(6,1) * t141 - Icges(6,4) * t140 + Icges(6,5) * t168;
t36 = -t138 * t86 + t139 * t90 + t166 * t82;
t101 = Icges(6,5) * t171 - Icges(6,6) * t170 + Icges(6,3) * t179;
t103 = Icges(6,4) * t171 - Icges(6,2) * t170 + Icges(6,6) * t179;
t105 = Icges(6,1) * t171 - Icges(6,4) * t170 + Icges(6,5) * t179;
t52 = t101 * t166 - t103 * t138 + t105 * t139;
t2 = t166 * t35 + t168 * t36 + t179 * t52;
t229 = t2 / 0.2e1 + t1 / 0.2e1;
t37 = t140 * t79 + t141 * t87 + t168 * t83;
t38 = t140 * t80 + t141 * t88 + t168 * t84;
t53 = t100 * t140 + t102 * t168 + t104 * t141;
t3 = t166 * t37 + t168 * t38 + t179 * t53;
t39 = -t140 * t85 + t141 * t89 + t168 * t81;
t40 = -t140 * t86 + t141 * t90 + t168 * t82;
t54 = t101 * t168 - t103 * t140 + t105 * t141;
t4 = t166 * t39 + t168 * t40 + t179 * t54;
t228 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t33 * t186 + t34 * t188 - t240 * t51;
t6 = t35 * t186 + t36 * t188 - t240 * t52;
t227 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t37 * t186 + t38 * t188 - t240 * t53;
t8 = t39 * t186 + t40 * t188 - t240 * t54;
t226 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t201 * t52 + (t197 * t36 - t200 * t35) * t198;
t9 = t201 * t51 + (t197 * t34 - t200 * t33) * t198;
t225 = t9 / 0.2e1 + t10 / 0.2e1;
t224 = -t135 - t246;
t11 = t201 * t53 + (t197 * t38 - t200 * t37) * t198;
t12 = t201 * t54 + (t197 * t40 - t200 * t39) * t198;
t223 = t12 / 0.2e1 + t11 / 0.2e1;
t43 = t170 * t79 + t171 * t87 + t179 * t83;
t44 = t170 * t80 + t171 * t88 + t179 * t84;
t58 = t100 * t170 + t102 * t179 + t104 * t171;
t13 = t166 * t43 + t168 * t44 + t179 * t58;
t45 = -t170 * t85 + t171 * t89 + t179 * t81;
t46 = -t170 * t86 + t171 * t90 + t179 * t82;
t59 = t101 * t179 - t103 * t170 + t105 * t171;
t14 = t166 * t45 + t168 * t46 + t179 * t59;
t222 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t43 * t186 + t44 * t188 - t240 * t58;
t16 = t45 * t186 + t46 * t188 - t240 * t59;
t221 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t201 * t58 + (t197 * t44 - t200 * t43) * t198;
t18 = t201 * t59 + (t197 * t46 - t200 * t45) * t198;
t220 = t18 / 0.2e1 + t17 / 0.2e1;
t219 = -m(4) - m(5) - m(6) - m(7);
t216 = t201 * t135 + t236;
t215 = -t134 + t234;
t214 = -t156 + t232;
t212 = (-t185 * rSges(4,1) - t184 * rSges(4,2) + rSges(4,3) * t240 - t190) * t198;
t211 = t117 * t243 + t118 * t242 + t231;
t145 = t180 * rSges(5,1) - t179 * rSges(5,2) - rSges(5,3) * t240;
t210 = (-t145 + t232) * t198;
t108 = rSges(6,1) * t171 - rSges(6,2) * t170 + rSges(6,3) * t179;
t208 = (-t108 + t214) * t198;
t207 = t134 * t243 + t135 * t242 + t211;
t206 = (t214 - t235) * t198;
t181 = t201 * rSges(3,3) + (rSges(3,1) * t204 + rSges(3,2) * t205) * t198;
t178 = Icges(3,5) * t201 + (Icges(3,1) * t204 + Icges(3,4) * t205) * t198;
t176 = Icges(3,3) * t201 + (Icges(3,5) * t204 + Icges(3,6) * t205) * t198;
t175 = t189 * t199 + t218;
t174 = -t189 * t196 + t199 * t243;
t173 = t187 * t199 - t217;
t172 = -t187 * t196 - t199 * t242;
t159 = rSges(3,1) * t189 - rSges(3,2) * t188 + rSges(3,3) * t243;
t158 = rSges(3,1) * t187 - rSges(3,2) * t186 - rSges(3,3) * t242;
t154 = Icges(3,1) * t189 - Icges(3,4) * t188 + Icges(3,5) * t243;
t153 = Icges(3,1) * t187 - Icges(3,4) * t186 - Icges(3,5) * t242;
t152 = Icges(3,4) * t189 - Icges(3,2) * t188 + Icges(3,6) * t243;
t151 = Icges(3,4) * t187 - Icges(3,2) * t186 - Icges(3,6) * t242;
t150 = Icges(3,5) * t189 - Icges(3,6) * t188 + Icges(3,3) * t243;
t149 = Icges(3,5) * t187 - Icges(3,6) * t186 - Icges(3,3) * t242;
t148 = Icges(4,1) * t185 + Icges(4,4) * t184 - Icges(4,5) * t240;
t147 = Icges(4,4) * t185 + Icges(4,2) * t184 - Icges(4,6) * t240;
t144 = Icges(5,1) * t180 - Icges(5,4) * t179 - Icges(5,5) * t240;
t143 = Icges(5,4) * t180 - Icges(5,2) * t179 - Icges(5,6) * t240;
t142 = Icges(5,5) * t180 - Icges(5,6) * t179 - Icges(5,3) * t240;
t133 = -t158 * t201 - t181 * t242;
t132 = t159 * t201 - t181 * t243;
t127 = rSges(4,1) * t175 + rSges(4,2) * t174 + rSges(4,3) * t188;
t126 = rSges(4,1) * t173 + rSges(4,2) * t172 + rSges(4,3) * t186;
t125 = Icges(4,1) * t175 + Icges(4,4) * t174 + Icges(4,5) * t188;
t124 = Icges(4,1) * t173 + Icges(4,4) * t172 + Icges(4,5) * t186;
t123 = Icges(4,4) * t175 + Icges(4,2) * t174 + Icges(4,6) * t188;
t122 = Icges(4,4) * t173 + Icges(4,2) * t172 + Icges(4,6) * t186;
t121 = Icges(4,5) * t175 + Icges(4,6) * t174 + Icges(4,3) * t188;
t120 = Icges(4,5) * t173 + Icges(4,6) * t172 + Icges(4,3) * t186;
t119 = t188 * t134;
t116 = rSges(5,1) * t169 - rSges(5,2) * t168 + rSges(5,3) * t188;
t115 = rSges(5,1) * t167 - rSges(5,2) * t166 + rSges(5,3) * t186;
t114 = Icges(5,1) * t169 - Icges(5,4) * t168 + Icges(5,5) * t188;
t113 = Icges(5,1) * t167 - Icges(5,4) * t166 + Icges(5,5) * t186;
t112 = Icges(5,4) * t169 - Icges(5,2) * t168 + Icges(5,6) * t188;
t111 = Icges(5,4) * t167 - Icges(5,2) * t166 + Icges(5,6) * t186;
t110 = Icges(5,5) * t169 - Icges(5,6) * t168 + Icges(5,3) * t188;
t109 = Icges(5,5) * t167 - Icges(5,6) * t166 + Icges(5,3) * t186;
t97 = (t158 * t197 + t159 * t200) * t198;
t92 = rSges(6,1) * t139 - rSges(6,2) * t138 + rSges(6,3) * t166;
t78 = -t116 * t240 - t188 * t145;
t77 = t115 * t240 + t186 * t145;
t76 = (-t126 - t164) * t201 + t200 * t212;
t75 = t127 * t201 + t197 * t212 + t163;
t74 = -t142 * t240 - t179 * t143 + t180 * t144;
t73 = t115 * t188 - t116 * t186;
t72 = t142 * t188 - t143 * t168 + t144 * t169;
t71 = t142 * t186 - t143 * t166 + t144 * t167;
t70 = (t126 * t197 + t127 * t200) * t198 + t231;
t69 = -t108 * t168 + t179 * t94;
t68 = t108 * t166 - t179 * t92;
t67 = -t110 * t240 - t179 * t112 + t180 * t114;
t66 = -t109 * t240 - t179 * t111 + t180 * t113;
t65 = t110 * t188 - t112 * t168 + t114 * t169;
t64 = t109 * t188 - t111 * t168 + t113 * t169;
t63 = t110 * t186 - t112 * t166 + t114 * t167;
t62 = t109 * t186 - t111 * t166 + t113 * t167;
t61 = (-t115 + t234) * t201 + t200 * t210;
t60 = t116 * t201 + t197 * t210 + t236;
t57 = -t166 * t94 + t168 * t92;
t56 = t245 * t240 + (-t108 - t156) * t188;
t55 = t186 * t108 + t240 * t92 + t233;
t50 = (t115 * t197 + t116 * t200) * t198 + t211;
t49 = t186 * t245 + t188 * t92 + t119;
t48 = -t168 * t235 + t179 * t246;
t47 = t166 * t235 - t179 * t247;
t42 = (-t92 + t215) * t201 + t200 * t208;
t41 = t197 * t208 + t201 * t94 + t216;
t32 = t224 * t240 + (-t156 - t235) * t188;
t31 = t186 * t235 + t240 * t247 + t233;
t30 = -t166 * t246 + t168 * t247;
t29 = (t197 * t92 + t200 * t94) * t198 + t207;
t28 = (t215 - t247) * t201 + t200 * t206;
t27 = t197 * t206 + t201 * t246 + t216;
t26 = t186 * t224 + t188 * t247 + t119;
t25 = t201 * t74 + (t197 * t67 - t200 * t66) * t198;
t24 = t66 * t186 + t67 * t188 - t240 * t74;
t23 = t201 * t72 + (t197 * t65 - t200 * t64) * t198;
t22 = t201 * t71 + (t197 * t63 - t200 * t62) * t198;
t21 = t64 * t186 + t65 * t188 - t240 * t72;
t20 = t62 * t186 + t63 * t188 - t240 * t71;
t19 = (t197 * t247 + t200 * t246) * t198 + t207;
t91 = [m(2) + m(3) - t219; m(3) * t97 + m(4) * t70 + m(5) * t50 + m(6) * t29 + m(7) * t19; m(7) * (t19 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t29 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t50 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t70 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(3) * (t132 ^ 2 + t133 ^ 2 + t97 ^ 2) + (t12 + t11 + t23 + ((t121 * t188 + t174 * t123 + t175 * t125) * t197 - (t188 * t120 + t174 * t122 + t175 * t124) * t200) * t198 + (t150 * t243 - t152 * t188 + t189 * t154) * t243) * t243 + (-t9 - t10 - t22 - ((t186 * t121 + t172 * t123 + t173 * t125) * t197 - (t120 * t186 + t122 * t172 + t124 * t173) * t200) * t198 + (-t149 * t242 - t151 * t186 + t153 * t187) * t242 + (-t149 * t243 + t150 * t242 + t188 * t151 + t186 * t152 - t189 * t153 - t187 * t154) * t243) * t242 + (t18 + t17 + t25 + ((t152 * t205 + t154 * t204) * t197 - (t151 * t205 + t153 * t204) * t200) * t251 + ((-t149 * t200 + t150 * t197 + t177 * t205 + t178 * t204) * t198 - t146 * t240 + t184 * t147 + t185 * t148 + t201 * t176) * t201 + (-t121 * t240 + t184 * t123 + t185 * t125 + t174 * t147 + t175 * t148 + t176 * t243 + t189 * t178 + t188 * t250) * t243 + (t120 * t240 - t184 * t122 - t185 * t124 - t172 * t147 - t173 * t148 + t176 * t242 - t187 * t178 - t186 * t250) * t242) * t201; t219 * t240; m(7) * (t186 * t27 + t188 * t28 - t19 * t240) + m(6) * (t186 * t41 + t188 * t42 - t240 * t29) + m(5) * (t186 * t60 + t188 * t61 - t240 * t50) + m(4) * (t186 * t75 + t188 * t76 - t240 * t70); 0.2e1 * (m(7) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t251 * t205 ^ 2 + t186 ^ 2 + t188 ^ 2); m(5) * t73 + m(6) * t49 + m(7) * t26; (t24 / 0.2e1 + t221) * t201 + (t23 / 0.2e1 + t223) * t188 + (t22 / 0.2e1 + t225) * t186 + m(7) * (t19 * t26 + t27 * t32 + t28 * t31) + m(6) * (t29 * t49 + t41 * t56 + t42 * t55) + m(5) * (t50 * t73 + t60 * t78 + t61 * t77) + ((-t25 / 0.2e1 - t220) * t205 + (-t20 / 0.2e1 - t227) * t200 + (t21 / 0.2e1 + t226) * t197) * t198; m(7) * (t32 * t186 + t31 * t188 - t240 * t26) + m(5) * (t78 * t186 + t77 * t188 - t240 * t73) + m(6) * (t56 * t186 + t55 * t188 - t240 * t49); (-t15 - t16 - t24) * t240 + (t8 + t7 + t21) * t188 + (t6 + t5 + t20) * t186 + m(7) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t49 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t73 ^ 2 + t77 ^ 2 + t78 ^ 2); m(6) * t57 + m(7) * t30; t222 * t201 + t220 * t179 + t223 * t168 + t225 * t166 + m(7) * (t19 * t30 + t27 * t48 + t28 * t47) + m(6) * (t29 * t57 + t41 * t69 + t42 * t68) + (t197 * t228 - t200 * t229) * t198; m(7) * (t48 * t186 + t47 * t188 - t240 * t30) + m(6) * (t69 * t186 + t68 * t188 - t240 * t57); -t222 * t240 + t228 * t188 + t229 * t186 + t221 * t179 + t226 * t168 + t227 * t166 + m(7) * (t26 * t30 + t31 * t47 + t32 * t48) + m(6) * (t49 * t57 + t55 * t68 + t56 * t69); (t13 + t14) * t179 + (t4 + t3) * t168 + (t1 + t2) * t166 + m(7) * (t30 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t57 ^ 2 + t68 ^ 2 + t69 ^ 2); m(7) * t170; m(7) * (t138 * t27 + t140 * t28 + t170 * t19); m(7) * (t138 * t186 + t140 * t188 - t170 * t240); m(7) * (t138 * t32 + t140 * t31 + t170 * t26); m(7) * (t138 * t48 + t140 * t47 + t170 * t30); m(7) * (t138 ^ 2 + t140 ^ 2 + t170 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t91(1) t91(2) t91(4) t91(7) t91(11) t91(16); t91(2) t91(3) t91(5) t91(8) t91(12) t91(17); t91(4) t91(5) t91(6) t91(9) t91(13) t91(18); t91(7) t91(8) t91(9) t91(10) t91(14) t91(19); t91(11) t91(12) t91(13) t91(14) t91(15) t91(20); t91(16) t91(17) t91(18) t91(19) t91(20) t91(21);];
Mq  = res;
