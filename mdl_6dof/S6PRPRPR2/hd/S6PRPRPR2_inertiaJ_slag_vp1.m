% Calculate joint inertia matrix for
% S6PRPRPR2
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:37
% EndTime: 2019-03-08 19:29:47
% DurationCPUTime: 5.15s
% Computational Cost: add. (20539->518), mult. (48092->770), div. (0->0), fcn. (63125->14), ass. (0->231)
t231 = m(6) / 0.2e1 + m(7) / 0.2e1;
t261 = 0.2e1 * t231;
t206 = sin(pkin(10));
t209 = cos(pkin(10));
t243 = sin(pkin(11));
t244 = cos(pkin(11));
t250 = sin(qJ(2));
t252 = cos(qJ(2));
t194 = -t250 * t243 + t252 * t244;
t210 = cos(pkin(6));
t213 = t210 * t194;
t214 = t243 * t252 + t244 * t250;
t175 = -t206 * t213 - t209 * t214;
t186 = t214 * t210;
t176 = -t186 * t206 + t194 * t209;
t207 = sin(pkin(6));
t239 = t206 * t207;
t128 = Icges(4,5) * t176 + Icges(4,6) * t175 + Icges(4,3) * t239;
t225 = t210 * t252;
t191 = -t206 * t225 - t209 * t250;
t224 = t210 * t250;
t192 = -t206 * t224 + t209 * t252;
t163 = Icges(3,5) * t192 + Icges(3,6) * t191 + Icges(3,3) * t239;
t260 = t128 + t163;
t184 = t194 * t207;
t185 = t214 * t207;
t259 = Icges(4,5) * t185 + Icges(4,6) * t184 + (Icges(3,5) * t250 + Icges(3,6) * t252) * t207 + (Icges(4,3) + Icges(3,3)) * t210;
t173 = -t206 * t214 + t209 * t213;
t174 = t186 * t209 + t194 * t206;
t238 = t207 * t209;
t127 = Icges(4,5) * t174 + Icges(4,6) * t173 - Icges(4,3) * t238;
t189 = -t206 * t250 + t209 * t225;
t190 = t206 * t252 + t209 * t224;
t162 = Icges(3,5) * t190 + Icges(3,6) * t189 - Icges(3,3) * t238;
t258 = -t162 - t127;
t212 = sin(qJ(4));
t251 = cos(qJ(4));
t226 = t207 * t251;
t151 = t174 * t212 + t209 * t226;
t153 = t176 * t212 - t206 * t226;
t177 = t185 * t212 - t210 * t251;
t237 = t207 * t212;
t152 = t174 * t251 - t209 * t237;
t204 = pkin(12) + qJ(6);
t201 = sin(t204);
t202 = cos(t204);
t117 = -t152 * t201 - t173 * t202;
t118 = t152 * t202 - t173 * t201;
t73 = Icges(7,5) * t118 + Icges(7,6) * t117 + Icges(7,3) * t151;
t75 = Icges(7,4) * t118 + Icges(7,2) * t117 + Icges(7,6) * t151;
t77 = Icges(7,1) * t118 + Icges(7,4) * t117 + Icges(7,5) * t151;
t28 = t117 * t75 + t118 * t77 + t151 * t73;
t154 = t176 * t251 + t206 * t237;
t119 = -t154 * t201 - t175 * t202;
t120 = t154 * t202 - t175 * t201;
t74 = Icges(7,5) * t120 + Icges(7,6) * t119 + Icges(7,3) * t153;
t76 = Icges(7,4) * t120 + Icges(7,2) * t119 + Icges(7,6) * t153;
t78 = Icges(7,1) * t120 + Icges(7,4) * t119 + Icges(7,5) * t153;
t29 = t117 * t76 + t118 * t78 + t151 * t74;
t178 = t185 * t251 + t210 * t212;
t147 = -t178 * t201 - t184 * t202;
t148 = t178 * t202 - t184 * t201;
t91 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t177;
t92 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t177;
t93 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t177;
t44 = t117 * t92 + t118 * t93 + t151 * t91;
t1 = t151 * t28 + t153 * t29 + t177 * t44;
t257 = -t1 / 0.2e1;
t37 = t147 * t75 + t148 * t77 + t177 * t73;
t38 = t147 * t76 + t148 * t78 + t177 * t74;
t51 = t147 * t92 + t148 * t93 + t177 * t91;
t11 = t151 * t37 + t153 * t38 + t177 * t51;
t256 = t11 / 0.2e1;
t255 = t151 / 0.2e1;
t254 = t153 / 0.2e1;
t253 = t177 / 0.2e1;
t249 = pkin(2) * t252;
t208 = cos(pkin(12));
t248 = pkin(5) * t208;
t205 = sin(pkin(12));
t242 = t173 * t205;
t79 = rSges(7,1) * t118 + rSges(7,2) * t117 + rSges(7,3) * t151;
t247 = -pkin(5) * t242 + pkin(9) * t151 + t152 * t248 + t79;
t241 = t175 * t205;
t80 = rSges(7,1) * t120 + rSges(7,2) * t119 + rSges(7,3) * t153;
t246 = -pkin(5) * t241 + pkin(9) * t153 + t154 * t248 + t80;
t240 = t184 * t205;
t94 = rSges(7,1) * t148 + rSges(7,2) * t147 + rSges(7,3) * t177;
t245 = pkin(5) * t240 - pkin(9) * t177 - t178 * t248 - t94;
t139 = pkin(3) * t176 - pkin(8) * t175;
t223 = pkin(2) * t224 - qJ(3) * t207;
t171 = -t206 * t223 + t209 * t249;
t161 = t210 * t171;
t235 = t210 * t139 + t161;
t138 = pkin(3) * t174 - pkin(8) * t173;
t170 = t206 * t249 + t209 * t223;
t234 = -t138 - t170;
t233 = t170 * t239 + t171 * t238;
t195 = pkin(2) * t207 * t250 + t210 * qJ(3);
t232 = -pkin(3) * t185 + pkin(8) * t184 - t195;
t230 = m(4) + m(5) + m(6) + m(7);
t116 = pkin(4) * t154 + qJ(5) * t153;
t229 = t210 * t116 + t235;
t115 = pkin(4) * t152 + qJ(5) * t151;
t228 = -t115 + t234;
t146 = pkin(4) * t178 + qJ(5) * t177;
t227 = -t146 + t232;
t222 = (-rSges(4,1) * t185 - rSges(4,2) * t184 - rSges(4,3) * t210 - t195) * t207;
t221 = t138 * t239 + t139 * t238 + t233;
t143 = rSges(5,1) * t178 - rSges(5,2) * t177 - rSges(5,3) * t184;
t220 = (-t143 + t232) * t207;
t149 = -t178 * t205 - t184 * t208;
t150 = t178 * t208 - t240;
t109 = rSges(6,1) * t150 + rSges(6,2) * t149 + rSges(6,3) * t177;
t217 = (-t109 + t227) * t207;
t216 = t115 * t239 + t116 * t238 + t221;
t215 = (t227 + t245) * t207;
t183 = t210 * rSges(3,3) + (rSges(3,1) * t250 + rSges(3,2) * t252) * t207;
t182 = Icges(3,5) * t210 + (Icges(3,1) * t250 + Icges(3,4) * t252) * t207;
t181 = Icges(3,6) * t210 + (Icges(3,4) * t250 + Icges(3,2) * t252) * t207;
t169 = rSges(3,1) * t192 + rSges(3,2) * t191 + rSges(3,3) * t239;
t168 = rSges(3,1) * t190 + rSges(3,2) * t189 - rSges(3,3) * t238;
t167 = Icges(3,1) * t192 + Icges(3,4) * t191 + Icges(3,5) * t239;
t166 = Icges(3,1) * t190 + Icges(3,4) * t189 - Icges(3,5) * t238;
t165 = Icges(3,4) * t192 + Icges(3,2) * t191 + Icges(3,6) * t239;
t164 = Icges(3,4) * t190 + Icges(3,2) * t189 - Icges(3,6) * t238;
t159 = Icges(4,1) * t185 + Icges(4,4) * t184 + Icges(4,5) * t210;
t158 = Icges(4,4) * t185 + Icges(4,2) * t184 + Icges(4,6) * t210;
t145 = -t168 * t210 - t183 * t238;
t144 = t169 * t210 - t183 * t239;
t142 = Icges(5,1) * t178 - Icges(5,4) * t177 - Icges(5,5) * t184;
t141 = Icges(5,4) * t178 - Icges(5,2) * t177 - Icges(5,6) * t184;
t140 = Icges(5,5) * t178 - Icges(5,6) * t177 - Icges(5,3) * t184;
t134 = rSges(4,1) * t176 + rSges(4,2) * t175 + rSges(4,3) * t239;
t133 = rSges(4,1) * t174 + rSges(4,2) * t173 - rSges(4,3) * t238;
t132 = Icges(4,1) * t176 + Icges(4,4) * t175 + Icges(4,5) * t239;
t131 = Icges(4,1) * t174 + Icges(4,4) * t173 - Icges(4,5) * t238;
t130 = Icges(4,4) * t176 + Icges(4,2) * t175 + Icges(4,6) * t239;
t129 = Icges(4,4) * t174 + Icges(4,2) * t173 - Icges(4,6) * t238;
t126 = t154 * t208 - t241;
t125 = -t154 * t205 - t175 * t208;
t124 = t152 * t208 - t242;
t123 = -t152 * t205 - t173 * t208;
t122 = (t168 * t206 + t169 * t209) * t207;
t121 = t173 * t146;
t111 = t184 * t116;
t108 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t177;
t107 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t177;
t106 = Icges(6,5) * t150 + Icges(6,6) * t149 + Icges(6,3) * t177;
t105 = t175 * t115;
t104 = rSges(5,1) * t154 - rSges(5,2) * t153 - rSges(5,3) * t175;
t103 = rSges(5,1) * t152 - rSges(5,2) * t151 - rSges(5,3) * t173;
t102 = Icges(5,1) * t154 - Icges(5,4) * t153 - Icges(5,5) * t175;
t101 = Icges(5,1) * t152 - Icges(5,4) * t151 - Icges(5,5) * t173;
t100 = Icges(5,4) * t154 - Icges(5,2) * t153 - Icges(5,6) * t175;
t99 = Icges(5,4) * t152 - Icges(5,2) * t151 - Icges(5,6) * t173;
t98 = Icges(5,5) * t154 - Icges(5,6) * t153 - Icges(5,3) * t175;
t97 = Icges(5,5) * t152 - Icges(5,6) * t151 - Icges(5,3) * t173;
t90 = (-t133 - t170) * t210 + t209 * t222;
t89 = t134 * t210 + t206 * t222 + t161;
t88 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t153;
t87 = rSges(6,1) * t124 + rSges(6,2) * t123 + rSges(6,3) * t151;
t86 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t153;
t85 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t151;
t84 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t153;
t83 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t151;
t82 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t153;
t81 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t151;
t70 = (t133 * t206 + t134 * t209) * t207 + t233;
t69 = -t104 * t184 + t143 * t175;
t68 = t103 * t184 - t143 * t173;
t67 = -t140 * t184 - t141 * t177 + t142 * t178;
t66 = -t103 * t175 + t104 * t173;
t65 = -t140 * t175 - t141 * t153 + t142 * t154;
t64 = -t140 * t173 - t141 * t151 + t142 * t152;
t63 = (-t103 + t234) * t210 + t209 * t220;
t62 = t104 * t210 + t206 * t220 + t235;
t61 = -t153 * t94 + t177 * t80;
t60 = t151 * t94 - t177 * t79;
t59 = -t100 * t177 + t102 * t178 - t184 * t98;
t58 = t101 * t178 - t177 * t99 - t184 * t97;
t57 = t106 * t177 + t107 * t149 + t108 * t150;
t56 = (t103 * t206 + t104 * t209) * t207 + t221;
t55 = -t100 * t153 + t102 * t154 - t175 * t98;
t54 = t101 * t154 - t153 * t99 - t175 * t97;
t53 = -t100 * t151 + t102 * t152 - t173 * t98;
t52 = t101 * t152 - t151 * t99 - t173 * t97;
t50 = -t151 * t80 + t153 * t79;
t49 = t106 * t153 + t107 * t125 + t108 * t126;
t48 = t106 * t151 + t107 * t123 + t108 * t124;
t47 = -t184 * t88 - t111 + (t109 + t146) * t175;
t46 = -t109 * t173 - t121 - (-t115 - t87) * t184;
t45 = t119 * t92 + t120 * t93 + t153 * t91;
t43 = (-t87 + t228) * t210 + t209 * t217;
t42 = t206 * t217 + t210 * t88 + t229;
t41 = -t175 * t87 - t105 + (t116 + t88) * t173;
t40 = t149 * t84 + t150 * t86 + t177 * t82;
t39 = t149 * t83 + t150 * t85 + t177 * t81;
t36 = t125 * t84 + t126 * t86 + t153 * t82;
t35 = t125 * t83 + t126 * t85 + t153 * t81;
t34 = t123 * t84 + t124 * t86 + t151 * t82;
t33 = t123 * t83 + t124 * t85 + t151 * t81;
t32 = (t206 * t87 + t209 * t88) * t207 + t216;
t31 = t119 * t76 + t120 * t78 + t153 * t74;
t30 = t119 * t75 + t120 * t77 + t153 * t73;
t27 = -t111 - t246 * t184 + (t146 - t245) * t175;
t26 = -t121 + t245 * t173 - (-t115 - t247) * t184;
t25 = (t228 - t247) * t210 + t209 * t215;
t24 = t206 * t215 + t210 * t246 + t229;
t23 = t210 * t67 + (t206 * t59 - t209 * t58) * t207;
t22 = -t173 * t58 - t175 * t59 - t184 * t67;
t21 = -t105 - t247 * t175 + (t116 + t246) * t173;
t20 = (t206 * t247 + t209 * t246) * t207 + t216;
t19 = t210 * t65 + (t206 * t55 - t209 * t54) * t207;
t18 = t210 * t64 + (t206 * t53 - t209 * t52) * t207;
t17 = -t173 * t54 - t175 * t55 - t184 * t65;
t16 = -t173 * t52 - t175 * t53 - t184 * t64;
t15 = t210 * t57 + (t206 * t40 - t209 * t39) * t207;
t14 = -t173 * t39 - t175 * t40 - t184 * t57;
t13 = t210 * t51 + (t206 * t38 - t209 * t37) * t207;
t12 = -t173 * t37 - t175 * t38 - t184 * t51;
t10 = t210 * t49 + (t206 * t36 - t209 * t35) * t207;
t9 = t210 * t48 + (t206 * t34 - t209 * t33) * t207;
t8 = -t173 * t35 - t175 * t36 - t184 * t49;
t7 = -t173 * t33 - t175 * t34 - t184 * t48;
t6 = t210 * t45 + (t206 * t31 - t209 * t30) * t207;
t5 = t210 * t44 + (t206 * t29 - t209 * t28) * t207;
t4 = -t173 * t30 - t175 * t31 - t184 * t45;
t3 = -t173 * t28 - t175 * t29 - t184 * t44;
t2 = t151 * t30 + t153 * t31 + t177 * t45;
t71 = [m(3) + m(2) + t230; m(3) * t122 + m(4) * t70 + m(5) * t56 + m(6) * t32 + m(7) * t20; m(7) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t56 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(4) * (t70 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(3) * (t122 ^ 2 + t144 ^ 2 + t145 ^ 2) + (t13 + t23 + t15 + (t158 * t184 + t159 * t185 + t259 * t210) * t210 + ((t130 * t184 + t132 * t185) * t206 - (t129 * t184 + t131 * t185) * t209 + (-t127 * t209 + t128 * t206 + t181 * t252 + t182 * t250) * t210) * t207) * t210 + (t6 + t19 + t10 + (t130 * t175 + t132 * t176 + t165 * t191 + t167 * t192 + t260 * t239) * t239 + ((t165 * t252 + t167 * t250) * t207 + t158 * t175 + t159 * t176 + t181 * t191 + t182 * t192 + t259 * t239 + t210 * t163) * t210) * t239 + (-t5 - t18 - t9 + (t129 * t173 + t131 * t174 + t164 * t189 + t166 * t190 + t258 * t238) * t238 + (-(t164 * t252 + t166 * t250) * t207 - t181 * t189 - t182 * t190 - t158 * t173 - t159 * t174 + t259 * t238 - t210 * t162) * t210 + (-t129 * t175 - t130 * t173 - t131 * t176 - t132 * t174 - t164 * t191 - t165 * t189 - t166 * t192 - t167 * t190 + t260 * t238 + t258 * t239) * t239) * t238; t230 * t210; m(7) * (t20 * t210 + (t206 * t25 - t209 * t24) * t207) + m(6) * (t210 * t32 + (t206 * t43 - t209 * t42) * t207) + m(5) * (t210 * t56 + (t206 * t63 - t209 * t62) * t207) + m(4) * (t210 * t70 + (t206 * t90 - t209 * t89) * t207); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t231) * (t210 ^ 2 + (t206 ^ 2 + t209 ^ 2) * t207 ^ 2); m(5) * t66 + m(6) * t41 + m(7) * t21; (t12 / 0.2e1 + t22 / 0.2e1 + t14 / 0.2e1) * t210 - (t13 / 0.2e1 + t23 / 0.2e1 + t15 / 0.2e1) * t184 + (-t6 / 0.2e1 - t19 / 0.2e1 - t10 / 0.2e1) * t175 + (-t5 / 0.2e1 - t18 / 0.2e1 - t9 / 0.2e1) * t173 + m(7) * (t20 * t21 + t24 * t27 + t25 * t26) + m(6) * (t32 * t41 + t42 * t47 + t43 * t46) + m(5) * (t56 * t66 + t62 * t69 + t63 * t68) + ((-t3 / 0.2e1 - t7 / 0.2e1 - t16 / 0.2e1) * t209 + (t4 / 0.2e1 + t8 / 0.2e1 + t17 / 0.2e1) * t206) * t207; m(5) * (t210 * t66 + (t206 * t68 - t209 * t69) * t207) + m(6) * (t210 * t41 + (t206 * t46 - t209 * t47) * t207) + m(7) * (t21 * t210 + (t206 * t26 - t209 * t27) * t207); -(t12 + t14 + t22) * t184 + (-t4 - t8 - t17) * t175 + (-t3 - t7 - t16) * t173 + m(7) * (t21 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t41 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t66 ^ 2 + t68 ^ 2 + t69 ^ 2); t177 * t261; m(7) * (t151 * t24 + t153 * t25 + t177 * t20) + m(6) * (t151 * t42 + t153 * t43 + t177 * t32); (t177 * t210 + (-t151 * t209 + t153 * t206) * t207) * t261; m(7) * (t151 * t27 + t153 * t26 + t177 * t21) + m(6) * (t151 * t47 + t153 * t46 + t177 * t41); (t151 ^ 2 + t153 ^ 2 + t177 ^ 2) * t261; m(7) * t50; t6 * t254 + t13 * t253 + m(7) * (t20 * t50 + t24 * t61 + t25 * t60) + t210 * t256 + t5 * t255 + (t206 * t2 / 0.2e1 + t209 * t257) * t207; m(7) * (t210 * t50 + (t206 * t60 - t209 * t61) * t207); m(7) * (t21 * t50 + t26 * t60 + t27 * t61) + t173 * t257 + t3 * t255 - t184 * t256 + t4 * t254 - t175 * t2 / 0.2e1 + t12 * t253; m(7) * (t151 * t61 + t153 * t60 + t177 * t50); t153 * t2 + t151 * t1 + t177 * t11 + m(7) * (t50 ^ 2 + t60 ^ 2 + t61 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t71(1) t71(2) t71(4) t71(7) t71(11) t71(16); t71(2) t71(3) t71(5) t71(8) t71(12) t71(17); t71(4) t71(5) t71(6) t71(9) t71(13) t71(18); t71(7) t71(8) t71(9) t71(10) t71(14) t71(19); t71(11) t71(12) t71(13) t71(14) t71(15) t71(20); t71(16) t71(17) t71(18) t71(19) t71(20) t71(21);];
Mq  = res;
