% Calculate joint inertia matrix for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:05
% EndTime: 2019-03-08 19:39:13
% DurationCPUTime: 3.75s
% Computational Cost: add. (18586->521), mult. (28395->771), div. (0->0), fcn. (35807->14), ass. (0->232)
t203 = sin(pkin(6));
t258 = t203 ^ 2;
t228 = m(6) / 0.2e1 + m(7) / 0.2e1;
t257 = 0.2e1 * t228;
t201 = sin(pkin(11));
t205 = cos(pkin(11));
t207 = cos(pkin(6));
t210 = sin(qJ(2));
t240 = t203 * t210;
t184 = -t201 * t240 + t205 * t207;
t243 = t201 * t207;
t185 = t205 * t240 + t243;
t211 = cos(qJ(2));
t239 = t203 * t211;
t146 = Icges(4,5) * t185 + Icges(4,6) * t184 - Icges(4,3) * t239;
t179 = Icges(3,6) * t207 + (Icges(3,4) * t210 + Icges(3,2) * t211) * t203;
t256 = -t179 + t146;
t202 = sin(pkin(10));
t206 = cos(pkin(10));
t238 = t207 * t210;
t187 = t202 * t211 + t206 * t238;
t229 = pkin(11) + qJ(4);
t197 = sin(t229);
t219 = cos(t229);
t215 = t203 * t219;
t170 = t187 * t197 + t206 * t215;
t189 = -t202 * t238 + t206 * t211;
t172 = t189 * t197 - t202 * t215;
t181 = t197 * t240 - t207 * t219;
t242 = t202 * t203;
t173 = t189 * t219 + t197 * t242;
t237 = t207 * t211;
t188 = t202 * t237 + t206 * t210;
t199 = pkin(12) + qJ(6);
t196 = sin(t199);
t198 = cos(t199);
t135 = -t173 * t196 + t188 * t198;
t136 = t173 * t198 + t188 * t196;
t241 = t203 * t206;
t171 = t187 * t219 - t197 * t241;
t186 = t202 * t210 - t206 * t237;
t133 = -t171 * t196 + t186 * t198;
t134 = t171 * t198 + t186 * t196;
t75 = Icges(7,5) * t134 + Icges(7,6) * t133 + Icges(7,3) * t170;
t77 = Icges(7,4) * t134 + Icges(7,2) * t133 + Icges(7,6) * t170;
t79 = Icges(7,1) * t134 + Icges(7,4) * t133 + Icges(7,5) * t170;
t31 = t135 * t77 + t136 * t79 + t172 * t75;
t76 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t172;
t78 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t172;
t80 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t172;
t32 = t135 * t78 + t136 * t80 + t172 * t76;
t182 = t207 * t197 + t210 * t215;
t164 = -t182 * t196 - t198 * t239;
t165 = t182 * t198 - t196 * t239;
t93 = Icges(7,5) * t165 + Icges(7,6) * t164 + Icges(7,3) * t181;
t94 = Icges(7,4) * t165 + Icges(7,2) * t164 + Icges(7,6) * t181;
t95 = Icges(7,1) * t165 + Icges(7,4) * t164 + Icges(7,5) * t181;
t46 = t135 * t94 + t136 * t95 + t172 * t93;
t2 = t170 * t31 + t172 * t32 + t181 * t46;
t255 = t2 / 0.2e1;
t254 = t170 / 0.2e1;
t253 = t172 / 0.2e1;
t252 = t181 / 0.2e1;
t251 = pkin(3) * t205;
t204 = cos(pkin(12));
t250 = pkin(5) * t204;
t200 = sin(pkin(12));
t245 = t186 * t200;
t81 = rSges(7,1) * t134 + rSges(7,2) * t133 + rSges(7,3) * t170;
t249 = pkin(5) * t245 + pkin(9) * t170 + t250 * t171 + t81;
t244 = t188 * t200;
t82 = rSges(7,1) * t136 + rSges(7,2) * t135 + rSges(7,3) * t172;
t248 = pkin(5) * t244 + pkin(9) * t172 + t250 * t173 + t82;
t225 = t200 * t239;
t97 = rSges(7,1) * t165 + rSges(7,2) * t164 + rSges(7,3) * t181;
t247 = -pkin(5) * t225 + pkin(9) * t181 + t250 * t182 + t97;
t132 = pkin(4) * t173 + qJ(5) * t172;
t140 = -t173 * t200 + t188 * t204;
t141 = t173 * t204 + t244;
t90 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t172;
t246 = -t132 - t90;
t224 = t201 * t242;
t115 = pkin(3) * t224 + pkin(8) * t188 + t251 * t189;
t167 = pkin(2) * t189 + qJ(3) * t188;
t163 = t207 * t167;
t234 = t207 * t115 + t163;
t223 = t201 * t241;
t114 = -pkin(3) * t223 + pkin(8) * t186 + t251 * t187;
t166 = pkin(2) * t187 + qJ(3) * t186;
t233 = -t114 - t166;
t131 = pkin(4) * t171 + qJ(5) * t170;
t149 = t182 * pkin(4) + t181 * qJ(5);
t232 = t131 * t239 + t186 * t149;
t190 = (pkin(2) * t210 - qJ(3) * t211) * t203;
t231 = -pkin(3) * t243 - (-pkin(8) * t211 + t251 * t210) * t203 - t190;
t230 = t166 * t242 + t167 * t241;
t227 = -t132 - t248;
t226 = -m(4) - m(5) - m(6) - m(7);
t222 = t207 * t132 + t234;
t221 = -t131 + t233;
t220 = -t149 + t231;
t218 = (-t185 * rSges(4,1) - t184 * rSges(4,2) + rSges(4,3) * t239 - t190) * t203;
t217 = t114 * t242 + t115 * t241 + t230;
t145 = t182 * rSges(5,1) - t181 * rSges(5,2) - rSges(5,3) * t239;
t216 = (-t145 + t231) * t203;
t168 = -t182 * t200 - t204 * t239;
t169 = t182 * t204 - t225;
t103 = rSges(6,1) * t169 + rSges(6,2) * t168 + rSges(6,3) * t181;
t214 = (-t103 + t220) * t203;
t213 = t131 * t242 + t132 * t241 + t217;
t212 = (t220 - t247) * t203;
t183 = t207 * rSges(3,3) + (rSges(3,1) * t210 + rSges(3,2) * t211) * t203;
t180 = Icges(3,5) * t207 + (Icges(3,1) * t210 + Icges(3,4) * t211) * t203;
t178 = Icges(3,3) * t207 + (Icges(3,5) * t210 + Icges(3,6) * t211) * t203;
t177 = t189 * t205 + t224;
t176 = -t189 * t201 + t205 * t242;
t175 = t187 * t205 - t223;
t174 = -t187 * t201 - t205 * t241;
t159 = rSges(3,1) * t189 - rSges(3,2) * t188 + rSges(3,3) * t242;
t158 = rSges(3,1) * t187 - rSges(3,2) * t186 - rSges(3,3) * t241;
t155 = Icges(3,1) * t189 - Icges(3,4) * t188 + Icges(3,5) * t242;
t154 = Icges(3,1) * t187 - Icges(3,4) * t186 - Icges(3,5) * t241;
t153 = Icges(3,4) * t189 - Icges(3,2) * t188 + Icges(3,6) * t242;
t152 = Icges(3,4) * t187 - Icges(3,2) * t186 - Icges(3,6) * t241;
t151 = Icges(3,5) * t189 - Icges(3,6) * t188 + Icges(3,3) * t242;
t150 = Icges(3,5) * t187 - Icges(3,6) * t186 - Icges(3,3) * t241;
t148 = Icges(4,1) * t185 + Icges(4,4) * t184 - Icges(4,5) * t239;
t147 = Icges(4,4) * t185 + Icges(4,2) * t184 - Icges(4,6) * t239;
t144 = Icges(5,1) * t182 - Icges(5,4) * t181 - Icges(5,5) * t239;
t143 = Icges(5,4) * t182 - Icges(5,2) * t181 - Icges(5,6) * t239;
t142 = Icges(5,5) * t182 - Icges(5,6) * t181 - Icges(5,3) * t239;
t139 = t171 * t204 + t245;
t138 = -t171 * t200 + t186 * t204;
t130 = -t158 * t207 - t183 * t241;
t129 = t159 * t207 - t183 * t242;
t124 = rSges(4,1) * t177 + rSges(4,2) * t176 + rSges(4,3) * t188;
t123 = rSges(4,1) * t175 + rSges(4,2) * t174 + rSges(4,3) * t186;
t122 = Icges(4,1) * t177 + Icges(4,4) * t176 + Icges(4,5) * t188;
t121 = Icges(4,1) * t175 + Icges(4,4) * t174 + Icges(4,5) * t186;
t120 = Icges(4,4) * t177 + Icges(4,2) * t176 + Icges(4,6) * t188;
t119 = Icges(4,4) * t175 + Icges(4,2) * t174 + Icges(4,6) * t186;
t118 = Icges(4,5) * t177 + Icges(4,6) * t176 + Icges(4,3) * t188;
t117 = Icges(4,5) * t175 + Icges(4,6) * t174 + Icges(4,3) * t186;
t116 = t188 * t131;
t113 = rSges(5,1) * t173 - rSges(5,2) * t172 + rSges(5,3) * t188;
t112 = rSges(5,1) * t171 - rSges(5,2) * t170 + rSges(5,3) * t186;
t111 = Icges(5,1) * t173 - Icges(5,4) * t172 + Icges(5,5) * t188;
t110 = Icges(5,1) * t171 - Icges(5,4) * t170 + Icges(5,5) * t186;
t109 = Icges(5,4) * t173 - Icges(5,2) * t172 + Icges(5,6) * t188;
t108 = Icges(5,4) * t171 - Icges(5,2) * t170 + Icges(5,6) * t186;
t107 = Icges(5,5) * t173 - Icges(5,6) * t172 + Icges(5,3) * t188;
t106 = Icges(5,5) * t171 - Icges(5,6) * t170 + Icges(5,3) * t186;
t102 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t181;
t101 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t181;
t100 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t181;
t96 = (t158 * t202 + t159 * t206) * t203;
t89 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t170;
t88 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t172;
t87 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t170;
t86 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t172;
t85 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t170;
t84 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t172;
t83 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t170;
t74 = -t113 * t239 - t188 * t145;
t73 = t112 * t239 + t186 * t145;
t70 = (-t123 - t166) * t207 + t206 * t218;
t69 = t124 * t207 + t202 * t218 + t163;
t68 = -t142 * t239 - t181 * t143 + t182 * t144;
t67 = t112 * t188 - t113 * t186;
t66 = t142 * t188 - t143 * t172 + t144 * t173;
t65 = t142 * t186 - t143 * t170 + t144 * t171;
t64 = (t123 * t202 + t124 * t206) * t203 + t230;
t63 = -t107 * t239 - t181 * t109 + t182 * t111;
t62 = -t106 * t239 - t181 * t108 + t182 * t110;
t61 = -t172 * t97 + t181 * t82;
t60 = t170 * t97 - t181 * t81;
t59 = t107 * t188 - t109 * t172 + t111 * t173;
t58 = t106 * t188 - t108 * t172 + t110 * t173;
t57 = t107 * t186 - t109 * t170 + t111 * t171;
t56 = t106 * t186 - t108 * t170 + t110 * t171;
t55 = (-t112 + t233) * t207 + t206 * t216;
t54 = t113 * t207 + t202 * t216 + t234;
t53 = t100 * t181 + t101 * t168 + t102 * t169;
t52 = -t170 * t82 + t172 * t81;
t51 = t164 * t94 + t165 * t95 + t181 * t93;
t50 = t246 * t239 + (-t103 - t149) * t188;
t49 = t186 * t103 + t239 * t89 + t232;
t48 = t100 * t172 + t101 * t140 + t102 * t141;
t47 = t100 * t170 + t101 * t138 + t102 * t139;
t45 = t133 * t94 + t134 * t95 + t170 * t93;
t44 = (t112 * t202 + t113 * t206) * t203 + t217;
t43 = t186 * t246 + t188 * t89 + t116;
t42 = t168 * t86 + t169 * t88 + t181 * t84;
t41 = t168 * t85 + t169 * t87 + t181 * t83;
t40 = t164 * t78 + t165 * t80 + t181 * t76;
t39 = t164 * t77 + t165 * t79 + t181 * t75;
t38 = (-t89 + t221) * t207 + t206 * t214;
t37 = t202 * t214 + t207 * t90 + t222;
t36 = t140 * t86 + t141 * t88 + t172 * t84;
t35 = t140 * t85 + t141 * t87 + t172 * t83;
t34 = t138 * t86 + t139 * t88 + t170 * t84;
t33 = t138 * t85 + t139 * t87 + t170 * t83;
t30 = t133 * t78 + t134 * t80 + t170 * t76;
t29 = t133 * t77 + t134 * t79 + t170 * t75;
t28 = t227 * t239 + (-t149 - t247) * t188;
t27 = t186 * t247 + t239 * t249 + t232;
t26 = (t202 * t89 + t206 * t90) * t203 + t213;
t25 = t207 * t68 + (t202 * t63 - t206 * t62) * t203;
t24 = t62 * t186 + t63 * t188 - t239 * t68;
t23 = (t221 - t249) * t207 + t206 * t212;
t22 = t202 * t212 + t207 * t248 + t222;
t21 = t186 * t227 + t188 * t249 + t116;
t20 = t207 * t66 + (t202 * t59 - t206 * t58) * t203;
t19 = t207 * t65 + (t202 * t57 - t206 * t56) * t203;
t18 = t58 * t186 + t59 * t188 - t239 * t66;
t17 = t56 * t186 + t57 * t188 - t239 * t65;
t16 = (t202 * t249 + t206 * t248) * t203 + t213;
t15 = t207 * t53 + (t202 * t42 - t206 * t41) * t203;
t14 = t41 * t186 + t42 * t188 - t239 * t53;
t13 = t207 * t51 + (t202 * t40 - t206 * t39) * t203;
t12 = t39 * t186 + t40 * t188 - t239 * t51;
t11 = t170 * t39 + t172 * t40 + t181 * t51;
t10 = t207 * t48 + (t202 * t36 - t206 * t35) * t203;
t9 = t207 * t47 + (t202 * t34 - t206 * t33) * t203;
t8 = t35 * t186 + t36 * t188 - t239 * t48;
t7 = t33 * t186 + t34 * t188 - t239 * t47;
t6 = t207 * t46 + (t202 * t32 - t206 * t31) * t203;
t5 = t207 * t45 + (t202 * t30 - t206 * t29) * t203;
t4 = t31 * t186 + t32 * t188 - t239 * t46;
t3 = t29 * t186 + t30 * t188 - t239 * t45;
t1 = t170 * t29 + t172 * t30 + t181 * t45;
t71 = [m(2) + m(3) - t226; m(3) * t96 + m(4) * t64 + m(5) * t44 + m(6) * t26 + m(7) * t16; m(7) * (t16 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(6) * (t26 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t44 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t64 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(3) * (t129 ^ 2 + t130 ^ 2 + t96 ^ 2) + (t6 + t20 + t10 + ((t118 * t188 + t176 * t120 + t177 * t122) * t202 - (t188 * t117 + t176 * t119 + t177 * t121) * t206) * t203 + (t151 * t242 - t153 * t188 + t189 * t155) * t242) * t242 + (-t5 - t19 - t9 - ((t186 * t118 + t174 * t120 + t175 * t122) * t202 - (t117 * t186 + t119 * t174 + t121 * t175) * t206) * t203 + (-t150 * t241 - t152 * t186 + t154 * t187) * t241 + (-t150 * t242 + t151 * t241 + t188 * t152 + t186 * t153 - t189 * t154 - t187 * t155) * t242) * t241 + (t13 + t25 + t15 + ((t153 * t211 + t155 * t210) * t202 - (t152 * t211 + t154 * t210) * t206) * t258 + ((-t150 * t206 + t151 * t202 + t179 * t211 + t180 * t210) * t203 - t146 * t239 + t184 * t147 + t185 * t148 + t207 * t178) * t207 + (-t118 * t239 + t184 * t120 + t185 * t122 + t147 * t176 + t148 * t177 + t178 * t242 + t180 * t189 + t188 * t256) * t242 + (t117 * t239 - t184 * t119 - t185 * t121 - t147 * t174 - t148 * t175 + t178 * t241 - t180 * t187 - t186 * t256) * t241) * t207; t226 * t239; m(7) * (-t16 * t239 + t186 * t22 + t188 * t23) + m(6) * (t186 * t37 + t188 * t38 - t239 * t26) + m(5) * (t186 * t54 + t188 * t55 - t239 * t44) + m(4) * (t186 * t69 + t188 * t70 - t239 * t64); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t228) * (t211 ^ 2 * t258 + t186 ^ 2 + t188 ^ 2); m(5) * t67 + m(6) * t43 + m(7) * t21; (t12 / 0.2e1 + t24 / 0.2e1 + t14 / 0.2e1) * t207 + (t6 / 0.2e1 + t20 / 0.2e1 + t10 / 0.2e1) * t188 + (t5 / 0.2e1 + t19 / 0.2e1 + t9 / 0.2e1) * t186 + m(7) * (t16 * t21 + t22 * t28 + t23 * t27) + m(6) * (t26 * t43 + t37 * t50 + t38 * t49) + m(5) * (t44 * t67 + t54 * t74 + t55 * t73) + ((-t13 / 0.2e1 - t25 / 0.2e1 - t15 / 0.2e1) * t211 + (-t3 / 0.2e1 - t17 / 0.2e1 - t7 / 0.2e1) * t206 + (t4 / 0.2e1 + t18 / 0.2e1 + t8 / 0.2e1) * t202) * t203; m(5) * (t74 * t186 + t73 * t188 - t239 * t67) + m(6) * (t50 * t186 + t49 * t188 - t239 * t43) + m(7) * (t28 * t186 + t27 * t188 - t21 * t239); (-t12 - t14 - t24) * t239 + (t4 + t8 + t18) * t188 + (t3 + t7 + t17) * t186 + m(7) * (t21 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t43 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t67 ^ 2 + t73 ^ 2 + t74 ^ 2); t181 * t257; m(7) * (t16 * t181 + t170 * t22 + t172 * t23) + m(6) * (t170 * t37 + t172 * t38 + t181 * t26); (t170 * t186 + t172 * t188 - t181 * t239) * t257; m(7) * (t170 * t28 + t172 * t27 + t181 * t21) + m(6) * (t170 * t50 + t172 * t49 + t181 * t43); (t170 ^ 2 + t172 ^ 2 + t181 ^ 2) * t257; m(7) * t52; t207 * t11 / 0.2e1 + t13 * t252 + m(7) * (t16 * t52 + t22 * t61 + t23 * t60) + t5 * t254 + t6 * t253 + (t202 * t255 - t206 * t1 / 0.2e1) * t203; m(7) * (t61 * t186 + t60 * t188 - t239 * t52); m(7) * (t21 * t52 + t27 * t60 + t28 * t61) + t4 * t253 + t186 * t1 / 0.2e1 + t12 * t252 + t188 * t255 + t3 * t254 - t11 * t239 / 0.2e1; m(7) * (t170 * t61 + t172 * t60 + t181 * t52); t172 * t2 + t170 * t1 + t181 * t11 + m(7) * (t52 ^ 2 + t60 ^ 2 + t61 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t71(1) t71(2) t71(4) t71(7) t71(11) t71(16); t71(2) t71(3) t71(5) t71(8) t71(12) t71(17); t71(4) t71(5) t71(6) t71(9) t71(13) t71(18); t71(7) t71(8) t71(9) t71(10) t71(14) t71(19); t71(11) t71(12) t71(13) t71(14) t71(15) t71(20); t71(16) t71(17) t71(18) t71(19) t71(20) t71(21);];
Mq  = res;
