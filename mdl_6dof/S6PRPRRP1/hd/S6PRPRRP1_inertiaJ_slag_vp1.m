% Calculate joint inertia matrix for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:18
% EndTime: 2019-03-08 19:55:25
% DurationCPUTime: 5.28s
% Computational Cost: add. (23803->531), mult. (60698->778), div. (0->0), fcn. (79954->12), ass. (0->228)
t259 = rSges(7,3) + qJ(6);
t245 = sin(pkin(11));
t246 = cos(pkin(11));
t253 = sin(qJ(2));
t255 = cos(qJ(2));
t192 = -t245 * t253 + t246 * t255;
t201 = sin(pkin(6));
t182 = t192 * t201;
t209 = t245 * t255 + t246 * t253;
t183 = t209 * t201;
t203 = cos(pkin(6));
t258 = Icges(4,5) * t183 + Icges(4,6) * t182 + (Icges(3,5) * t253 + Icges(3,6) * t255) * t201 + (Icges(4,3) + Icges(3,3)) * t203;
t200 = sin(pkin(10));
t202 = cos(pkin(10));
t208 = t203 * t192;
t171 = -t200 * t209 + t202 * t208;
t184 = t209 * t203;
t172 = t184 * t202 + t192 * t200;
t240 = t201 * t202;
t127 = Icges(4,5) * t172 + Icges(4,6) * t171 - Icges(4,3) * t240;
t220 = t203 * t255;
t187 = -t200 * t253 + t202 * t220;
t219 = t203 * t253;
t188 = t200 * t255 + t202 * t219;
t160 = Icges(3,5) * t188 + Icges(3,6) * t187 - Icges(3,3) * t240;
t257 = -t160 - t127;
t173 = -t200 * t208 - t202 * t209;
t174 = -t184 * t200 + t192 * t202;
t241 = t200 * t201;
t128 = Icges(4,5) * t174 + Icges(4,6) * t173 + Icges(4,3) * t241;
t189 = -t200 * t220 - t202 * t253;
t190 = -t200 * t219 + t202 * t255;
t161 = Icges(3,5) * t190 + Icges(3,6) * t189 + Icges(3,3) * t241;
t256 = t161 + t128;
t254 = cos(qJ(4));
t252 = t255 * pkin(2);
t207 = cos(qJ(5));
t251 = pkin(5) * t207;
t206 = sin(qJ(4));
t239 = t201 * t206;
t150 = t172 * t254 - t202 * t239;
t205 = sin(qJ(5));
t123 = -t150 * t205 - t171 * t207;
t244 = t171 * t205;
t124 = t150 * t207 - t244;
t221 = t201 * t254;
t149 = t172 * t206 + t202 * t221;
t249 = rSges(7,1) * t124 + rSges(7,2) * t123 - pkin(5) * t244 + t149 * t259 + t251 * t150;
t152 = t174 * t254 + t200 * t239;
t125 = -t152 * t205 - t173 * t207;
t243 = t173 * t205;
t126 = t152 * t207 - t243;
t151 = t174 * t206 - t200 * t221;
t248 = rSges(7,1) * t126 + rSges(7,2) * t125 - pkin(5) * t243 + t151 * t259 + t251 * t152;
t176 = t183 * t254 + t203 * t206;
t147 = -t176 * t205 - t182 * t207;
t242 = t182 * t205;
t148 = t176 * t207 - t242;
t175 = t183 * t206 - t203 * t254;
t247 = rSges(7,1) * t148 + rSges(7,2) * t147 - pkin(5) * t242 + t175 * t259 + t251 * t176;
t139 = pkin(3) * t174 - pkin(8) * t173;
t218 = pkin(2) * t219 - qJ(3) * t201;
t169 = -t200 * t218 + t202 * t252;
t159 = t203 * t169;
t238 = t139 * t203 + t159;
t138 = pkin(3) * t172 - pkin(8) * t171;
t168 = t200 * t252 + t202 * t218;
t237 = -t138 - t168;
t236 = t168 * t241 + t169 * t240;
t193 = pkin(2) * t201 * t253 + qJ(3) * t203;
t235 = -pkin(3) * t183 + pkin(8) * t182 - t193;
t79 = Icges(7,5) * t124 + Icges(7,6) * t123 + Icges(7,3) * t149;
t83 = Icges(7,4) * t124 + Icges(7,2) * t123 + Icges(7,6) * t149;
t87 = Icges(7,1) * t124 + Icges(7,4) * t123 + Icges(7,5) * t149;
t33 = t123 * t83 + t124 * t87 + t149 * t79;
t80 = Icges(7,5) * t126 + Icges(7,6) * t125 + Icges(7,3) * t151;
t84 = Icges(7,4) * t126 + Icges(7,2) * t125 + Icges(7,6) * t151;
t88 = Icges(7,1) * t126 + Icges(7,4) * t125 + Icges(7,5) * t151;
t34 = t123 * t84 + t124 * t88 + t149 * t80;
t107 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t175;
t109 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t175;
t111 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t175;
t52 = t107 * t149 + t109 * t123 + t111 * t124;
t1 = t149 * t33 + t151 * t34 + t175 * t52;
t81 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t149;
t85 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t149;
t89 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t149;
t35 = t123 * t85 + t124 * t89 + t149 * t81;
t82 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t151;
t86 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t151;
t90 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t151;
t36 = t123 * t86 + t124 * t90 + t149 * t82;
t108 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t175;
t110 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t175;
t112 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t175;
t53 = t108 * t149 + t110 * t123 + t112 * t124;
t2 = t149 * t35 + t151 * t36 + t175 * t53;
t234 = -t2 / 0.2e1 - t1 / 0.2e1;
t37 = t125 * t83 + t126 * t87 + t151 * t79;
t38 = t125 * t84 + t126 * t88 + t151 * t80;
t54 = t107 * t151 + t109 * t125 + t111 * t126;
t3 = t149 * t37 + t151 * t38 + t175 * t54;
t39 = t125 * t85 + t126 * t89 + t151 * t81;
t40 = t125 * t86 + t126 * t90 + t151 * t82;
t55 = t108 * t151 + t110 * t125 + t112 * t126;
t4 = t149 * t39 + t151 * t40 + t175 * t55;
t233 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = -t171 * t33 - t173 * t34 - t182 * t52;
t6 = -t171 * t35 - t173 * t36 - t182 * t53;
t232 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = -t171 * t37 - t173 * t38 - t182 * t54;
t8 = -t171 * t39 - t173 * t40 - t182 * t55;
t231 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t203 * t53 + (t200 * t36 - t202 * t35) * t201;
t9 = t203 * t52 + (t200 * t34 - t202 * t33) * t201;
t230 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t203 * t54 + (t200 * t38 - t202 * t37) * t201;
t12 = t203 * t55 + (t200 * t40 - t202 * t39) * t201;
t229 = t12 / 0.2e1 + t11 / 0.2e1;
t43 = t147 * t83 + t148 * t87 + t175 * t79;
t44 = t147 * t84 + t148 * t88 + t175 * t80;
t62 = t107 * t175 + t109 * t147 + t111 * t148;
t13 = t149 * t43 + t151 * t44 + t175 * t62;
t45 = t147 * t85 + t148 * t89 + t175 * t81;
t46 = t147 * t86 + t148 * t90 + t175 * t82;
t63 = t108 * t175 + t110 * t147 + t112 * t148;
t14 = t149 * t45 + t151 * t46 + t175 * t63;
t228 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = -t171 * t43 - t173 * t44 - t182 * t62;
t16 = -t171 * t45 - t173 * t46 - t182 * t63;
t227 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t203 * t62 + (t200 * t44 - t202 * t43) * t201;
t18 = t203 * t63 + (t200 * t46 - t202 * t45) * t201;
t226 = t18 / 0.2e1 + t17 / 0.2e1;
t225 = m(4) + m(5) + m(6) + m(7);
t120 = pkin(4) * t152 + pkin(9) * t151;
t224 = t120 * t203 + t238;
t119 = pkin(4) * t150 + pkin(9) * t149;
t223 = -t119 + t237;
t146 = pkin(4) * t176 + pkin(9) * t175;
t222 = -t146 + t235;
t217 = (-rSges(4,1) * t183 - rSges(4,2) * t182 - rSges(4,3) * t203 - t193) * t201;
t216 = t138 * t241 + t139 * t240 + t236;
t143 = rSges(5,1) * t176 - rSges(5,2) * t175 - rSges(5,3) * t182;
t215 = (-t143 + t235) * t201;
t114 = rSges(6,1) * t148 + rSges(6,2) * t147 + rSges(6,3) * t175;
t212 = (-t114 + t222) * t201;
t211 = t119 * t241 + t120 * t240 + t216;
t210 = (t222 - t247) * t201;
t181 = t203 * rSges(3,3) + (rSges(3,1) * t253 + rSges(3,2) * t255) * t201;
t180 = Icges(3,5) * t203 + (Icges(3,1) * t253 + Icges(3,4) * t255) * t201;
t179 = Icges(3,6) * t203 + (Icges(3,4) * t253 + Icges(3,2) * t255) * t201;
t167 = rSges(3,1) * t190 + rSges(3,2) * t189 + rSges(3,3) * t241;
t166 = rSges(3,1) * t188 + rSges(3,2) * t187 - rSges(3,3) * t240;
t165 = Icges(3,1) * t190 + Icges(3,4) * t189 + Icges(3,5) * t241;
t164 = Icges(3,1) * t188 + Icges(3,4) * t187 - Icges(3,5) * t240;
t163 = Icges(3,4) * t190 + Icges(3,2) * t189 + Icges(3,6) * t241;
t162 = Icges(3,4) * t188 + Icges(3,2) * t187 - Icges(3,6) * t240;
t157 = Icges(4,1) * t183 + Icges(4,4) * t182 + Icges(4,5) * t203;
t156 = Icges(4,4) * t183 + Icges(4,2) * t182 + Icges(4,6) * t203;
t145 = -t166 * t203 - t181 * t240;
t144 = t167 * t203 - t181 * t241;
t142 = Icges(5,1) * t176 - Icges(5,4) * t175 - Icges(5,5) * t182;
t141 = Icges(5,4) * t176 - Icges(5,2) * t175 - Icges(5,6) * t182;
t140 = Icges(5,5) * t176 - Icges(5,6) * t175 - Icges(5,3) * t182;
t134 = rSges(4,1) * t174 + rSges(4,2) * t173 + rSges(4,3) * t241;
t133 = rSges(4,1) * t172 + rSges(4,2) * t171 - rSges(4,3) * t240;
t132 = Icges(4,1) * t174 + Icges(4,4) * t173 + Icges(4,5) * t241;
t131 = Icges(4,1) * t172 + Icges(4,4) * t171 - Icges(4,5) * t240;
t130 = Icges(4,4) * t174 + Icges(4,2) * t173 + Icges(4,6) * t241;
t129 = Icges(4,4) * t172 + Icges(4,2) * t171 - Icges(4,6) * t240;
t122 = (t166 * t200 + t167 * t202) * t201;
t121 = t171 * t146;
t115 = t182 * t120;
t106 = t173 * t119;
t105 = rSges(5,1) * t152 - rSges(5,2) * t151 - rSges(5,3) * t173;
t104 = rSges(5,1) * t150 - rSges(5,2) * t149 - rSges(5,3) * t171;
t103 = Icges(5,1) * t152 - Icges(5,4) * t151 - Icges(5,5) * t173;
t102 = Icges(5,1) * t150 - Icges(5,4) * t149 - Icges(5,5) * t171;
t101 = Icges(5,4) * t152 - Icges(5,2) * t151 - Icges(5,6) * t173;
t100 = Icges(5,4) * t150 - Icges(5,2) * t149 - Icges(5,6) * t171;
t99 = Icges(5,5) * t152 - Icges(5,6) * t151 - Icges(5,3) * t173;
t98 = Icges(5,5) * t150 - Icges(5,6) * t149 - Icges(5,3) * t171;
t96 = (-t133 - t168) * t203 + t202 * t217;
t95 = t134 * t203 + t200 * t217 + t159;
t94 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t151;
t92 = rSges(6,1) * t124 + rSges(6,2) * t123 + rSges(6,3) * t149;
t76 = (t133 * t200 + t134 * t202) * t201 + t236;
t75 = -t105 * t182 + t143 * t173;
t74 = t104 * t182 - t143 * t171;
t73 = -t140 * t182 - t141 * t175 + t142 * t176;
t72 = -t104 * t173 + t105 * t171;
t71 = -t140 * t173 - t141 * t151 + t142 * t152;
t70 = -t140 * t171 - t141 * t149 + t142 * t150;
t69 = (-t104 + t237) * t203 + t202 * t215;
t68 = t105 * t203 + t200 * t215 + t238;
t67 = -t114 * t151 + t175 * t94;
t66 = t114 * t149 - t175 * t92;
t65 = -t101 * t175 + t103 * t176 - t182 * t99;
t64 = -t100 * t175 + t102 * t176 - t182 * t98;
t61 = (t104 * t200 + t105 * t202) * t201 + t216;
t60 = -t101 * t151 + t103 * t152 - t173 * t99;
t59 = -t100 * t151 + t102 * t152 - t173 * t98;
t58 = -t101 * t149 + t103 * t150 - t171 * t99;
t57 = -t100 * t149 + t102 * t150 - t171 * t98;
t56 = -t149 * t94 + t151 * t92;
t51 = -t182 * t94 - t115 + (t114 + t146) * t173;
t50 = -t114 * t171 - t121 - (-t119 - t92) * t182;
t49 = (-t92 + t223) * t203 + t202 * t212;
t48 = t200 * t212 + t203 * t94 + t224;
t47 = -t173 * t92 - t106 + (t120 + t94) * t171;
t42 = -t151 * t247 + t175 * t248;
t41 = t149 * t247 - t175 * t249;
t32 = (t200 * t92 + t202 * t94) * t201 + t211;
t31 = -t115 - t248 * t182 + (t146 + t247) * t173;
t30 = -t121 - t247 * t171 - (-t119 - t249) * t182;
t29 = (t223 - t249) * t203 + t202 * t210;
t28 = t200 * t210 + t203 * t248 + t224;
t27 = -t149 * t248 + t151 * t249;
t26 = t203 * t73 + (t200 * t65 - t202 * t64) * t201;
t25 = -t171 * t64 - t173 * t65 - t182 * t73;
t24 = -t106 - t249 * t173 + (t120 + t248) * t171;
t23 = (t200 * t249 + t202 * t248) * t201 + t211;
t22 = t203 * t71 + (t200 * t60 - t202 * t59) * t201;
t21 = t203 * t70 + (t200 * t58 - t202 * t57) * t201;
t20 = -t171 * t59 - t173 * t60 - t182 * t71;
t19 = -t171 * t57 - t173 * t58 - t182 * t70;
t77 = [m(2) + m(3) + t225; m(3) * t122 + m(4) * t76 + m(5) * t61 + m(6) * t32 + m(7) * t23; m(7) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t32 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t61 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t76 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(3) * (t122 ^ 2 + t144 ^ 2 + t145 ^ 2) + (t17 + t18 + t26 + (t182 * t156 + t183 * t157 + t258 * t203) * t203 + ((t130 * t182 + t132 * t183) * t200 - (t129 * t182 + t131 * t183) * t202 + (-t127 * t202 + t128 * t200 + t179 * t255 + t180 * t253) * t203) * t201) * t203 + (t12 + t11 + t22 + (t173 * t130 + t174 * t132 + t189 * t163 + t190 * t165 + t241 * t256) * t241 + ((t163 * t255 + t165 * t253) * t201 + t179 * t189 + t180 * t190 + t156 * t173 + t157 * t174 + t258 * t241 + t203 * t161) * t203) * t241 + (-t9 - t10 - t21 + (t171 * t129 + t172 * t131 + t187 * t162 + t188 * t164 + t240 * t257) * t240 + (-(t162 * t255 + t164 * t253) * t201 - t179 * t187 - t180 * t188 - t156 * t171 - t157 * t172 + t258 * t240 - t203 * t160) * t203 + (-t129 * t173 - t130 * t171 - t131 * t174 - t132 * t172 - t162 * t189 - t163 * t187 - t164 * t190 - t165 * t188 + t240 * t256 + t241 * t257) * t241) * t240; t225 * t203; m(7) * (t203 * t23 + (t200 * t29 - t202 * t28) * t201) + m(6) * (t203 * t32 + (t200 * t49 - t202 * t48) * t201) + m(5) * (t203 * t61 + (t200 * t69 - t202 * t68) * t201) + m(4) * (t203 * t76 + (t200 * t96 - t202 * t95) * t201); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t203 ^ 2 + (t200 ^ 2 + t202 ^ 2) * t201 ^ 2); m(5) * t72 + m(6) * t47 + m(7) * t24; (t25 / 0.2e1 + t227) * t203 - (t26 / 0.2e1 + t226) * t182 + (-t22 / 0.2e1 - t229) * t173 + (-t21 / 0.2e1 - t230) * t171 + m(7) * (t23 * t24 + t28 * t31 + t29 * t30) + m(6) * (t32 * t47 + t48 * t51 + t49 * t50) + m(5) * (t61 * t72 + t68 * t75 + t69 * t74) + ((-t19 / 0.2e1 - t232) * t202 + (t20 / 0.2e1 + t231) * t200) * t201; m(5) * (t203 * t72 + (t200 * t74 - t202 * t75) * t201) + m(6) * (t203 * t47 + (t200 * t50 - t202 * t51) * t201) + m(7) * (t203 * t24 + (t200 * t30 - t202 * t31) * t201); -(t15 + t16 + t25) * t182 + (-t8 - t7 - t20) * t173 + (-t5 - t6 - t19) * t171 + m(7) * (t24 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t72 ^ 2 + t74 ^ 2 + t75 ^ 2); m(6) * t56 + m(7) * t27; t228 * t203 + t226 * t175 + t229 * t151 + t230 * t149 + m(7) * (t23 * t27 + t28 * t42 + t29 * t41) + m(6) * (t32 * t56 + t48 * t67 + t49 * t66) + (t200 * t233 + t202 * t234) * t201; m(6) * (t203 * t56 + (t200 * t66 - t202 * t67) * t201) + m(7) * (t203 * t27 + (t200 * t41 - t202 * t42) * t201); -t228 * t182 + t227 * t175 - t233 * t173 + t234 * t171 + t231 * t151 + t232 * t149 + m(7) * (t24 * t27 + t30 * t41 + t31 * t42) + m(6) * (t47 * t56 + t50 * t66 + t51 * t67); (t14 + t13) * t175 + (t4 + t3) * t151 + (t2 + t1) * t149 + m(7) * (t27 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * t175; m(7) * (t149 * t28 + t151 * t29 + t175 * t23); m(7) * (t175 * t203 + (-t149 * t202 + t151 * t200) * t201); m(7) * (t149 * t31 + t151 * t30 + t175 * t24); m(7) * (t149 * t42 + t151 * t41 + t175 * t27); m(7) * (t149 ^ 2 + t151 ^ 2 + t175 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t77(1) t77(2) t77(4) t77(7) t77(11) t77(16); t77(2) t77(3) t77(5) t77(8) t77(12) t77(17); t77(4) t77(5) t77(6) t77(9) t77(13) t77(18); t77(7) t77(8) t77(9) t77(10) t77(14) t77(19); t77(11) t77(12) t77(13) t77(14) t77(15) t77(20); t77(16) t77(17) t77(18) t77(19) t77(20) t77(21);];
Mq  = res;
