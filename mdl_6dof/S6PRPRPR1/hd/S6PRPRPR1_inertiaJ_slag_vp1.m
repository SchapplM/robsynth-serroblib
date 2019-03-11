% Calculate joint inertia matrix for
% S6PRPRPR1
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:45
% EndTime: 2019-03-08 19:24:55
% DurationCPUTime: 5.23s
% Computational Cost: add. (21014->519), mult. (43728->771), div. (0->0), fcn. (56960->14), ass. (0->232)
t235 = m(6) / 0.2e1 + m(7) / 0.2e1;
t264 = 0.2e1 * t235;
t205 = sin(pkin(10));
t207 = cos(pkin(10));
t246 = sin(pkin(11));
t247 = cos(pkin(11));
t254 = sin(qJ(2));
t255 = cos(qJ(2));
t196 = -t254 * t246 + t255 * t247;
t208 = cos(pkin(6));
t214 = t208 * t196;
t215 = t255 * t246 + t254 * t247;
t171 = -t205 * t215 + t207 * t214;
t188 = t215 * t208;
t173 = t188 * t207 + t196 * t205;
t206 = sin(pkin(6));
t244 = t206 * t207;
t122 = Icges(4,5) * t173 + Icges(4,6) * t171 - Icges(4,3) * t244;
t228 = t208 * t255;
t191 = -t205 * t254 + t207 * t228;
t227 = t208 * t254;
t192 = t205 * t255 + t207 * t227;
t160 = Icges(3,5) * t192 + Icges(3,6) * t191 - Icges(3,3) * t244;
t263 = -t160 - t122;
t174 = -t205 * t214 - t207 * t215;
t176 = -t188 * t205 + t196 * t207;
t245 = t205 * t206;
t123 = Icges(4,5) * t176 + Icges(4,6) * t174 + Icges(4,3) * t245;
t193 = -t205 * t228 - t207 * t254;
t194 = -t205 * t227 + t207 * t255;
t161 = Icges(3,5) * t194 + Icges(3,6) * t193 + Icges(3,3) * t245;
t262 = t161 + t123;
t186 = t196 * t206;
t187 = t215 * t206;
t261 = Icges(4,5) * t187 + Icges(4,6) * t186 + (t254 * Icges(3,5) + t255 * Icges(3,6)) * t206 + (Icges(4,3) + Icges(3,3)) * t208;
t236 = qJ(4) + pkin(12);
t203 = sin(t236);
t225 = cos(t236);
t219 = t206 * t225;
t145 = t173 * t203 + t207 * t219;
t147 = t176 * t203 - t205 * t219;
t177 = t187 * t203 - t208 * t225;
t146 = t173 * t225 - t203 * t244;
t210 = sin(qJ(6));
t212 = cos(qJ(6));
t112 = -t146 * t210 - t171 * t212;
t113 = t146 * t212 - t171 * t210;
t71 = Icges(7,5) * t113 + Icges(7,6) * t112 + Icges(7,3) * t145;
t73 = Icges(7,4) * t113 + Icges(7,2) * t112 + Icges(7,6) * t145;
t75 = Icges(7,1) * t113 + Icges(7,4) * t112 + Icges(7,5) * t145;
t26 = t112 * t73 + t113 * t75 + t145 * t71;
t148 = t176 * t225 + t203 * t245;
t114 = -t148 * t210 - t174 * t212;
t115 = t148 * t212 - t174 * t210;
t72 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t147;
t74 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t147;
t76 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t147;
t27 = t112 * t74 + t113 * t76 + t145 * t72;
t178 = t187 * t225 + t208 * t203;
t143 = -t178 * t210 - t186 * t212;
t144 = t178 * t212 - t186 * t210;
t96 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t177;
t97 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t177;
t98 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t177;
t38 = t112 * t97 + t113 * t98 + t145 * t96;
t1 = t145 * t26 + t147 * t27 + t177 * t38;
t260 = -t1 / 0.2e1;
t33 = t143 * t73 + t144 * t75 + t177 * t71;
t34 = t143 * t74 + t144 * t76 + t177 * t72;
t47 = t143 * t97 + t144 * t98 + t177 * t96;
t7 = t145 * t33 + t147 * t34 + t177 * t47;
t259 = t7 / 0.2e1;
t258 = t145 / 0.2e1;
t257 = t147 / 0.2e1;
t256 = t177 / 0.2e1;
t253 = t255 * pkin(2);
t213 = cos(qJ(4));
t252 = pkin(4) * t213;
t77 = rSges(7,1) * t113 + rSges(7,2) * t112 + rSges(7,3) * t145;
t250 = pkin(5) * t146 + pkin(9) * t145 + t77;
t78 = rSges(7,1) * t115 + rSges(7,2) * t114 + rSges(7,3) * t147;
t249 = pkin(5) * t148 + pkin(9) * t147 + t78;
t99 = rSges(7,1) * t144 + rSges(7,2) * t143 + rSges(7,3) * t177;
t248 = -pkin(5) * t178 - pkin(9) * t177 - t99;
t211 = sin(qJ(4));
t243 = t206 * t211;
t242 = t206 * t213;
t241 = t208 * t211;
t135 = pkin(3) * t176 - pkin(8) * t174;
t226 = pkin(2) * t227 - qJ(3) * t206;
t169 = -t226 * t205 + t253 * t207;
t159 = t208 * t169;
t240 = t208 * t135 + t159;
t134 = pkin(3) * t173 - pkin(8) * t171;
t168 = t253 * t205 + t226 * t207;
t239 = -t134 - t168;
t238 = t168 * t245 + t169 * t244;
t197 = t206 * t254 * pkin(2) + t208 * qJ(3);
t237 = -pkin(3) * t187 + pkin(8) * t186 - t197;
t231 = t205 * t243;
t87 = pkin(4) * t231 - qJ(5) * t174 + t252 * t176;
t234 = t208 * t87 + t240;
t230 = t207 * t243;
t86 = -pkin(4) * t230 - qJ(5) * t171 + t252 * t173;
t233 = -t86 + t239;
t232 = m(4) + m(5) + m(6) + m(7);
t118 = pkin(4) * t241 - qJ(5) * t186 + t252 * t187;
t229 = -t118 + t237;
t224 = (-rSges(4,1) * t187 - rSges(4,2) * t186 - rSges(4,3) * t208 - t197) * t206;
t223 = t134 * t245 + t135 * t244 + t238;
t179 = -t187 * t211 + t208 * t213;
t180 = t187 * t213 + t241;
t139 = rSges(5,1) * t180 + rSges(5,2) * t179 - rSges(5,3) * t186;
t222 = (-t139 + t237) * t206;
t128 = rSges(6,1) * t178 - rSges(6,2) * t177 - rSges(6,3) * t186;
t218 = (-t128 + t229) * t206;
t217 = t87 * t244 + t86 * t245 + t223;
t216 = (t229 + t248) * t206;
t185 = t208 * rSges(3,3) + (t254 * rSges(3,1) + t255 * rSges(3,2)) * t206;
t184 = Icges(3,5) * t208 + (t254 * Icges(3,1) + t255 * Icges(3,4)) * t206;
t183 = Icges(3,6) * t208 + (t254 * Icges(3,4) + t255 * Icges(3,2)) * t206;
t167 = rSges(3,1) * t194 + rSges(3,2) * t193 + rSges(3,3) * t245;
t166 = rSges(3,1) * t192 + rSges(3,2) * t191 - rSges(3,3) * t244;
t165 = Icges(3,1) * t194 + Icges(3,4) * t193 + Icges(3,5) * t245;
t164 = Icges(3,1) * t192 + Icges(3,4) * t191 - Icges(3,5) * t244;
t163 = Icges(3,4) * t194 + Icges(3,2) * t193 + Icges(3,6) * t245;
t162 = Icges(3,4) * t192 + Icges(3,2) * t191 - Icges(3,6) * t244;
t157 = Icges(4,1) * t187 + Icges(4,4) * t186 + Icges(4,5) * t208;
t156 = Icges(4,4) * t187 + Icges(4,2) * t186 + Icges(4,6) * t208;
t152 = t176 * t213 + t231;
t151 = -t176 * t211 + t205 * t242;
t150 = t173 * t213 - t230;
t149 = -t173 * t211 - t207 * t242;
t141 = -t166 * t208 - t185 * t244;
t140 = t167 * t208 - t185 * t245;
t138 = Icges(5,1) * t180 + Icges(5,4) * t179 - Icges(5,5) * t186;
t137 = Icges(5,4) * t180 + Icges(5,2) * t179 - Icges(5,6) * t186;
t136 = Icges(5,5) * t180 + Icges(5,6) * t179 - Icges(5,3) * t186;
t130 = rSges(4,1) * t176 + rSges(4,2) * t174 + rSges(4,3) * t245;
t129 = rSges(4,1) * t173 + rSges(4,2) * t171 - rSges(4,3) * t244;
t127 = Icges(4,1) * t176 + Icges(4,4) * t174 + Icges(4,5) * t245;
t126 = Icges(4,1) * t173 + Icges(4,4) * t171 - Icges(4,5) * t244;
t125 = Icges(4,4) * t176 + Icges(4,2) * t174 + Icges(4,6) * t245;
t124 = Icges(4,4) * t173 + Icges(4,2) * t171 - Icges(4,6) * t244;
t121 = Icges(6,1) * t178 - Icges(6,4) * t177 - Icges(6,5) * t186;
t120 = Icges(6,4) * t178 - Icges(6,2) * t177 - Icges(6,6) * t186;
t119 = Icges(6,5) * t178 - Icges(6,6) * t177 - Icges(6,3) * t186;
t117 = (t166 * t205 + t167 * t207) * t206;
t108 = t171 * t118;
t107 = rSges(5,1) * t152 + rSges(5,2) * t151 - rSges(5,3) * t174;
t106 = rSges(5,1) * t150 + rSges(5,2) * t149 - rSges(5,3) * t171;
t105 = Icges(5,1) * t152 + Icges(5,4) * t151 - Icges(5,5) * t174;
t104 = Icges(5,1) * t150 + Icges(5,4) * t149 - Icges(5,5) * t171;
t103 = Icges(5,4) * t152 + Icges(5,2) * t151 - Icges(5,6) * t174;
t102 = Icges(5,4) * t150 + Icges(5,2) * t149 - Icges(5,6) * t171;
t101 = Icges(5,5) * t152 + Icges(5,6) * t151 - Icges(5,3) * t174;
t100 = Icges(5,5) * t150 + Icges(5,6) * t149 - Icges(5,3) * t171;
t95 = rSges(6,1) * t148 - rSges(6,2) * t147 - rSges(6,3) * t174;
t94 = rSges(6,1) * t146 - rSges(6,2) * t145 - rSges(6,3) * t171;
t93 = Icges(6,1) * t148 - Icges(6,4) * t147 - Icges(6,5) * t174;
t92 = Icges(6,1) * t146 - Icges(6,4) * t145 - Icges(6,5) * t171;
t91 = Icges(6,4) * t148 - Icges(6,2) * t147 - Icges(6,6) * t174;
t90 = Icges(6,4) * t146 - Icges(6,2) * t145 - Icges(6,6) * t171;
t89 = Icges(6,5) * t148 - Icges(6,6) * t147 - Icges(6,3) * t174;
t88 = Icges(6,5) * t146 - Icges(6,6) * t145 - Icges(6,3) * t171;
t82 = t186 * t87;
t81 = (-t129 - t168) * t208 + t207 * t224;
t80 = t130 * t208 + t205 * t224 + t159;
t79 = t174 * t86;
t70 = (t129 * t205 + t130 * t207) * t206 + t238;
t69 = -t107 * t186 + t139 * t174;
t68 = t106 * t186 - t139 * t171;
t67 = -t136 * t186 + t137 * t179 + t138 * t180;
t66 = -t119 * t186 - t120 * t177 + t121 * t178;
t65 = -t106 * t174 + t107 * t171;
t64 = -t136 * t174 + t137 * t151 + t138 * t152;
t63 = -t136 * t171 + t137 * t149 + t138 * t150;
t62 = -t119 * t174 - t120 * t147 + t121 * t148;
t61 = -t119 * t171 - t120 * t145 + t121 * t146;
t60 = (-t106 + t239) * t208 + t207 * t222;
t59 = t107 * t208 + t205 * t222 + t240;
t58 = -t147 * t99 + t177 * t78;
t57 = t145 * t99 - t177 * t77;
t56 = -t101 * t186 + t103 * t179 + t105 * t180;
t55 = -t100 * t186 + t102 * t179 + t104 * t180;
t54 = -t177 * t91 + t178 * t93 - t186 * t89;
t53 = -t177 * t90 + t178 * t92 - t186 * t88;
t52 = (t106 * t205 + t107 * t207) * t206 + t223;
t51 = -t101 * t174 + t103 * t151 + t105 * t152;
t50 = -t100 * t174 + t102 * t151 + t104 * t152;
t49 = -t101 * t171 + t103 * t149 + t105 * t150;
t48 = -t100 * t171 + t102 * t149 + t104 * t150;
t46 = -t147 * t91 + t148 * t93 - t174 * t89;
t45 = -t147 * t90 + t148 * t92 - t174 * t88;
t44 = -t145 * t91 + t146 * t93 - t171 * t89;
t43 = -t145 * t90 + t146 * t92 - t171 * t88;
t42 = -t186 * t95 - t82 + (t118 + t128) * t174;
t41 = -t128 * t171 - t108 - (-t86 - t94) * t186;
t40 = -t145 * t78 + t147 * t77;
t39 = t114 * t97 + t115 * t98 + t147 * t96;
t37 = (-t94 + t233) * t208 + t207 * t218;
t36 = t205 * t218 + t208 * t95 + t234;
t35 = -t174 * t94 - t79 + (t87 + t95) * t171;
t32 = (t205 * t94 + t207 * t95) * t206 + t217;
t31 = -t82 - t249 * t186 + (t118 - t248) * t174;
t30 = -t108 + t248 * t171 - (-t86 - t250) * t186;
t29 = t114 * t74 + t115 * t76 + t147 * t72;
t28 = t114 * t73 + t115 * t75 + t147 * t71;
t25 = (t233 - t250) * t208 + t207 * t216;
t24 = t205 * t216 + t249 * t208 + t234;
t23 = -t79 - t250 * t174 + (t87 + t249) * t171;
t22 = t208 * t67 + (t205 * t56 - t207 * t55) * t206;
t21 = (t250 * t205 + t249 * t207) * t206 + t217;
t20 = -t171 * t55 - t174 * t56 - t186 * t67;
t19 = t208 * t66 + (t205 * t54 - t207 * t53) * t206;
t18 = -t171 * t53 - t174 * t54 - t186 * t66;
t17 = t208 * t64 + (t205 * t51 - t207 * t50) * t206;
t16 = t208 * t63 + (t205 * t49 - t207 * t48) * t206;
t15 = -t171 * t50 - t174 * t51 - t186 * t64;
t14 = -t171 * t48 - t174 * t49 - t186 * t63;
t13 = t208 * t62 + (t205 * t46 - t207 * t45) * t206;
t12 = t208 * t61 + (t205 * t44 - t207 * t43) * t206;
t11 = -t171 * t45 - t174 * t46 - t186 * t62;
t10 = -t171 * t43 - t174 * t44 - t186 * t61;
t9 = t208 * t47 + (t205 * t34 - t207 * t33) * t206;
t8 = -t171 * t33 - t174 * t34 - t186 * t47;
t6 = t208 * t39 + (t205 * t29 - t207 * t28) * t206;
t5 = t208 * t38 + (t205 * t27 - t207 * t26) * t206;
t4 = -t171 * t28 - t174 * t29 - t186 * t39;
t3 = -t171 * t26 - t174 * t27 - t186 * t38;
t2 = t145 * t28 + t147 * t29 + t177 * t39;
t83 = [m(2) + m(3) + t232; m(3) * t117 + m(4) * t70 + m(5) * t52 + m(6) * t32 + m(7) * t21; m(7) * (t21 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t52 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t70 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(3) * (t117 ^ 2 + t140 ^ 2 + t141 ^ 2) + (t9 + t22 + t19 + (t186 * t156 + t187 * t157 + t261 * t208) * t208 + ((t186 * t125 + t187 * t127) * t205 - (t186 * t124 + t187 * t126) * t207 + (-t122 * t207 + t123 * t205 + t255 * t183 + t254 * t184) * t208) * t206) * t208 + (t6 + t17 + t13 + (t174 * t125 + t176 * t127 + t193 * t163 + t194 * t165 + t262 * t245) * t245 + ((t255 * t163 + t254 * t165) * t206 + t193 * t183 + t194 * t184 + t174 * t156 + t176 * t157 + t261 * t245 + t208 * t161) * t208) * t245 + (-t5 - t16 - t12 + (t171 * t124 + t173 * t126 + t191 * t162 + t192 * t164 + t263 * t244) * t244 + (-(t255 * t162 + t254 * t164) * t206 - t191 * t183 - t192 * t184 - t171 * t156 - t173 * t157 + t261 * t244 - t208 * t160) * t208 + (-t174 * t124 - t171 * t125 - t176 * t126 - t173 * t127 - t193 * t162 - t191 * t163 - t194 * t164 - t192 * t165 + t262 * t244 + t263 * t245) * t245) * t244; t232 * t208; m(7) * (t208 * t21 + (t205 * t25 - t207 * t24) * t206) + m(6) * (t208 * t32 + (t205 * t37 - t207 * t36) * t206) + m(5) * (t208 * t52 + (t205 * t60 - t207 * t59) * t206) + m(4) * (t208 * t70 + (t205 * t81 - t207 * t80) * t206); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t235) * (t208 ^ 2 + (t205 ^ 2 + t207 ^ 2) * t206 ^ 2); m(5) * t65 + m(6) * t35 + m(7) * t23; (t8 / 0.2e1 + t20 / 0.2e1 + t18 / 0.2e1) * t208 - (t9 / 0.2e1 + t22 / 0.2e1 + t19 / 0.2e1) * t186 + (-t6 / 0.2e1 - t17 / 0.2e1 - t13 / 0.2e1) * t174 + (-t5 / 0.2e1 - t16 / 0.2e1 - t12 / 0.2e1) * t171 + m(7) * (t21 * t23 + t24 * t31 + t25 * t30) + m(6) * (t32 * t35 + t36 * t42 + t37 * t41) + m(5) * (t52 * t65 + t59 * t69 + t60 * t68) + ((-t3 / 0.2e1 - t14 / 0.2e1 - t10 / 0.2e1) * t207 + (t4 / 0.2e1 + t15 / 0.2e1 + t11 / 0.2e1) * t205) * t206; m(6) * (t208 * t35 + (t205 * t41 - t207 * t42) * t206) + m(7) * (t208 * t23 + (t205 * t30 - t207 * t31) * t206) + m(5) * (t208 * t65 + (t205 * t68 - t207 * t69) * t206); -(t8 + t18 + t20) * t186 + (-t4 - t15 - t11) * t174 + (-t3 - t14 - t10) * t171 + m(7) * (t23 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t35 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t65 ^ 2 + t68 ^ 2 + t69 ^ 2); -t186 * t264; m(7) * (-t171 * t24 - t174 * t25 - t186 * t21) + m(6) * (-t171 * t36 - t174 * t37 - t186 * t32); (-t186 * t208 + (t171 * t207 - t174 * t205) * t206) * t264; m(7) * (-t171 * t31 - t174 * t30 - t186 * t23) + m(6) * (-t171 * t42 - t174 * t41 - t186 * t35); (t171 ^ 2 + t174 ^ 2 + t186 ^ 2) * t264; m(7) * t40; m(7) * (t21 * t40 + t24 * t58 + t25 * t57) + t9 * t256 + t5 * t258 + t6 * t257 + t208 * t259 + (t207 * t260 + t205 * t2 / 0.2e1) * t206; m(7) * (t208 * t40 + (t205 * t57 - t207 * t58) * t206); t8 * t256 + t3 * t258 + t4 * t257 + m(7) * (t23 * t40 + t30 * t57 + t31 * t58) + t171 * t260 - t186 * t259 - t174 * t2 / 0.2e1; m(7) * (-t171 * t58 - t174 * t57 - t186 * t40); m(7) * (t40 ^ 2 + t57 ^ 2 + t58 ^ 2) + t147 * t2 + t145 * t1 + t177 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t83(1) t83(2) t83(4) t83(7) t83(11) t83(16); t83(2) t83(3) t83(5) t83(8) t83(12) t83(17); t83(4) t83(5) t83(6) t83(9) t83(13) t83(18); t83(7) t83(8) t83(9) t83(10) t83(14) t83(19); t83(11) t83(12) t83(13) t83(14) t83(15) t83(20); t83(16) t83(17) t83(18) t83(19) t83(20) t83(21);];
Mq  = res;
