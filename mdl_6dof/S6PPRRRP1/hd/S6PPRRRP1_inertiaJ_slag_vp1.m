% Calculate joint inertia matrix for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:24
% EndTime: 2019-03-08 18:52:36
% DurationCPUTime: 5.06s
% Computational Cost: add. (35990->476), mult. (99327->684), div. (0->0), fcn. (131388->14), ass. (0->211)
t232 = rSges(7,3) + qJ(6);
t231 = cos(qJ(3));
t230 = cos(qJ(4));
t187 = cos(qJ(5));
t229 = pkin(5) * t187;
t179 = sin(pkin(11));
t181 = cos(pkin(11));
t182 = cos(pkin(6));
t219 = sin(pkin(12));
t196 = t182 * t219;
t221 = cos(pkin(12));
t172 = t179 * t221 + t181 * t196;
t186 = sin(qJ(3));
t197 = t182 * t221;
t191 = t179 * t219 - t181 * t197;
t222 = cos(pkin(7));
t189 = t191 * t222;
t180 = sin(pkin(6));
t220 = sin(pkin(7));
t198 = t180 * t220;
t156 = t172 * t231 + (-t181 * t198 - t189) * t186;
t199 = t180 * t222;
t165 = -t181 * t199 + t191 * t220;
t185 = sin(qJ(4));
t145 = t156 * t230 + t165 * t185;
t194 = t231 * t220;
t192 = t180 * t194;
t155 = t172 * t186 + t181 * t192 + t189 * t231;
t184 = sin(qJ(5));
t121 = -t145 * t184 + t155 * t187;
t218 = t155 * t184;
t122 = t145 * t187 + t218;
t144 = t156 * t185 - t165 * t230;
t227 = rSges(7,1) * t122 + rSges(7,2) * t121 + pkin(5) * t218 + t232 * t144 + t145 * t229;
t173 = -t179 * t196 + t181 * t221;
t190 = t179 * t197 + t181 * t219;
t188 = t190 * t222;
t158 = t173 * t231 + (t179 * t198 - t188) * t186;
t166 = t179 * t199 + t190 * t220;
t147 = t158 * t230 + t166 * t185;
t157 = t173 * t186 - t179 * t192 + t188 * t231;
t123 = -t147 * t184 + t157 * t187;
t217 = t157 * t184;
t124 = t147 * t187 + t217;
t146 = t158 * t185 - t166 * t230;
t226 = rSges(7,1) * t124 + rSges(7,2) * t123 + pkin(5) * t217 + t232 * t146 + t147 * t229;
t193 = t222 * t221;
t164 = t182 * t220 * t186 + (t186 * t193 + t219 * t231) * t180;
t171 = t182 * t222 - t198 * t221;
t160 = t164 * t230 + t171 * t185;
t163 = -t182 * t194 + (t186 * t219 - t193 * t231) * t180;
t148 = -t160 * t184 + t163 * t187;
t216 = t163 * t184;
t149 = t160 * t187 + t216;
t159 = t164 * t185 - t171 * t230;
t225 = rSges(7,1) * t149 + rSges(7,2) * t148 + pkin(5) * t216 + t232 * t159 + t160 * t229;
t118 = pkin(4) * t145 + pkin(10) * t144;
t91 = rSges(6,1) * t122 + rSges(6,2) * t121 + rSges(6,3) * t144;
t224 = -t118 - t91;
t119 = pkin(4) * t147 + pkin(10) * t146;
t93 = rSges(6,1) * t124 + rSges(6,2) * t123 + rSges(6,3) * t146;
t223 = -t119 - t93;
t113 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t159;
t142 = pkin(4) * t160 + pkin(10) * t159;
t215 = -t113 - t142;
t140 = pkin(3) * t156 + pkin(9) * t155;
t137 = t166 * t140;
t214 = t166 * t118 + t137;
t141 = pkin(3) * t158 + pkin(9) * t157;
t139 = t171 * t141;
t213 = t171 * t119 + t139;
t154 = pkin(3) * t164 + pkin(9) * t163;
t143 = t165 * t154;
t212 = t165 * t142 + t143;
t78 = Icges(7,5) * t122 + Icges(7,6) * t121 + Icges(7,3) * t144;
t82 = Icges(7,4) * t122 + Icges(7,2) * t121 + Icges(7,6) * t144;
t86 = Icges(7,1) * t122 + Icges(7,4) * t121 + Icges(7,5) * t144;
t34 = t121 * t82 + t122 * t86 + t144 * t78;
t79 = Icges(7,5) * t124 + Icges(7,6) * t123 + Icges(7,3) * t146;
t83 = Icges(7,4) * t124 + Icges(7,2) * t123 + Icges(7,6) * t146;
t87 = Icges(7,1) * t124 + Icges(7,4) * t123 + Icges(7,5) * t146;
t35 = t121 * t83 + t122 * t87 + t144 * t79;
t106 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t159;
t108 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t159;
t110 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t159;
t50 = t106 * t144 + t108 * t121 + t110 * t122;
t1 = t144 * t34 + t146 * t35 + t159 * t50;
t80 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t144;
t84 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t144;
t88 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t144;
t36 = t121 * t84 + t122 * t88 + t144 * t80;
t81 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t146;
t85 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t146;
t89 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t146;
t37 = t121 * t85 + t122 * t89 + t144 * t81;
t107 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t159;
t109 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t159;
t111 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t159;
t51 = t107 * t144 + t109 * t121 + t111 * t122;
t2 = t144 * t36 + t146 * t37 + t159 * t51;
t211 = t1 / 0.2e1 + t2 / 0.2e1;
t38 = t123 * t82 + t124 * t86 + t146 * t78;
t39 = t123 * t83 + t124 * t87 + t146 * t79;
t52 = t106 * t146 + t108 * t123 + t110 * t124;
t3 = t144 * t38 + t146 * t39 + t159 * t52;
t40 = t123 * t84 + t124 * t88 + t146 * t80;
t41 = t123 * t85 + t124 * t89 + t146 * t81;
t53 = t107 * t146 + t109 * t123 + t111 * t124;
t4 = t144 * t40 + t146 * t41 + t159 * t53;
t210 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t155 * t34 + t157 * t35 + t163 * t50;
t6 = t155 * t36 + t157 * t37 + t163 * t51;
t209 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t155 * t38 + t157 * t39 + t163 * t52;
t8 = t155 * t40 + t157 * t41 + t163 * t53;
t208 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t165 * t36 + t166 * t37 + t171 * t51;
t9 = t165 * t34 + t166 * t35 + t171 * t50;
t207 = t10 / 0.2e1 + t9 / 0.2e1;
t206 = -t118 - t227;
t205 = -t119 - t226;
t11 = t165 * t38 + t166 * t39 + t171 * t52;
t12 = t165 * t40 + t166 * t41 + t171 * t53;
t204 = t11 / 0.2e1 + t12 / 0.2e1;
t43 = t148 * t82 + t149 * t86 + t159 * t78;
t44 = t148 * t83 + t149 * t87 + t159 * t79;
t61 = t106 * t159 + t108 * t148 + t110 * t149;
t13 = t144 * t43 + t146 * t44 + t159 * t61;
t45 = t148 * t84 + t149 * t88 + t159 * t80;
t46 = t148 * t85 + t149 * t89 + t159 * t81;
t62 = t107 * t159 + t109 * t148 + t111 * t149;
t14 = t144 * t45 + t146 * t46 + t159 * t62;
t203 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t155 * t43 + t157 * t44 + t163 * t61;
t16 = t155 * t45 + t157 * t46 + t163 * t62;
t202 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t165 * t43 + t166 * t44 + t171 * t61;
t18 = t165 * t45 + t166 * t46 + t171 * t62;
t201 = t18 / 0.2e1 + t17 / 0.2e1;
t200 = -t142 - t225;
t195 = m(3) + m(4) + m(5) + m(6) + m(7);
t153 = rSges(4,1) * t164 - rSges(4,2) * t163 + rSges(4,3) * t171;
t152 = Icges(4,1) * t164 - Icges(4,4) * t163 + Icges(4,5) * t171;
t151 = Icges(4,4) * t164 - Icges(4,2) * t163 + Icges(4,6) * t171;
t150 = Icges(4,5) * t164 - Icges(4,6) * t163 + Icges(4,3) * t171;
t136 = rSges(5,1) * t160 - rSges(5,2) * t159 + rSges(5,3) * t163;
t135 = Icges(5,1) * t160 - Icges(5,4) * t159 + Icges(5,5) * t163;
t134 = Icges(5,4) * t160 - Icges(5,2) * t159 + Icges(5,6) * t163;
t133 = Icges(5,5) * t160 - Icges(5,6) * t159 + Icges(5,3) * t163;
t132 = rSges(4,1) * t158 - rSges(4,2) * t157 + rSges(4,3) * t166;
t131 = rSges(4,1) * t156 - rSges(4,2) * t155 + rSges(4,3) * t165;
t130 = Icges(4,1) * t158 - Icges(4,4) * t157 + Icges(4,5) * t166;
t129 = Icges(4,1) * t156 - Icges(4,4) * t155 + Icges(4,5) * t165;
t128 = Icges(4,4) * t158 - Icges(4,2) * t157 + Icges(4,6) * t166;
t127 = Icges(4,4) * t156 - Icges(4,2) * t155 + Icges(4,6) * t165;
t126 = Icges(4,5) * t158 - Icges(4,6) * t157 + Icges(4,3) * t166;
t125 = Icges(4,5) * t156 - Icges(4,6) * t155 + Icges(4,3) * t165;
t120 = t155 * t142;
t115 = t163 * t119;
t114 = t157 * t118;
t105 = rSges(5,1) * t147 - rSges(5,2) * t146 + rSges(5,3) * t157;
t104 = rSges(5,1) * t145 - rSges(5,2) * t144 + rSges(5,3) * t155;
t103 = Icges(5,1) * t147 - Icges(5,4) * t146 + Icges(5,5) * t157;
t102 = Icges(5,1) * t145 - Icges(5,4) * t144 + Icges(5,5) * t155;
t101 = Icges(5,4) * t147 - Icges(5,2) * t146 + Icges(5,6) * t157;
t100 = Icges(5,4) * t145 - Icges(5,2) * t144 + Icges(5,6) * t155;
t99 = Icges(5,5) * t147 - Icges(5,6) * t146 + Icges(5,3) * t157;
t98 = Icges(5,5) * t145 - Icges(5,6) * t144 + Icges(5,3) * t155;
t96 = t132 * t171 - t153 * t166;
t95 = -t131 * t171 + t153 * t165;
t94 = t131 * t166 - t132 * t165;
t75 = t105 * t163 - t136 * t157;
t74 = -t104 * t163 + t136 * t155;
t73 = t133 * t163 - t134 * t159 + t135 * t160;
t72 = t104 * t157 - t105 * t155;
t71 = t105 * t171 + t139 + (-t136 - t154) * t166;
t70 = t136 * t165 + t143 + (-t104 - t140) * t171;
t69 = t133 * t157 - t134 * t146 + t135 * t147;
t68 = t133 * t155 - t134 * t144 + t135 * t145;
t67 = -t113 * t146 + t159 * t93;
t66 = t113 * t144 - t159 * t91;
t65 = t104 * t166 + t137 + (-t105 - t141) * t165;
t64 = -t101 * t159 + t103 * t160 + t163 * t99;
t63 = -t100 * t159 + t102 * t160 + t163 * t98;
t60 = -t101 * t146 + t103 * t147 + t157 * t99;
t59 = -t100 * t146 + t102 * t147 + t157 * t98;
t58 = -t101 * t144 + t103 * t145 + t155 * t99;
t57 = -t100 * t144 + t102 * t145 + t155 * t98;
t56 = -t144 * t93 + t146 * t91;
t55 = t157 * t215 + t163 * t93 + t115;
t54 = t113 * t155 + t163 * t224 + t120;
t49 = t171 * t93 + (-t154 + t215) * t166 + t213;
t48 = t113 * t165 + (-t140 + t224) * t171 + t212;
t47 = t155 * t223 + t157 * t91 + t114;
t42 = t166 * t91 + (-t141 + t223) * t165 + t214;
t33 = -t146 * t225 + t159 * t226;
t32 = t144 * t225 - t159 * t227;
t31 = t157 * t200 + t163 * t226 + t115;
t30 = t155 * t225 + t163 * t206 + t120;
t29 = t226 * t171 + (-t154 + t200) * t166 + t213;
t28 = t225 * t165 + (-t140 + t206) * t171 + t212;
t27 = -t144 * t226 + t146 * t227;
t26 = t155 * t205 + t157 * t227 + t114;
t25 = t165 * t63 + t166 * t64 + t171 * t73;
t24 = t155 * t63 + t157 * t64 + t163 * t73;
t23 = t227 * t166 + (-t141 + t205) * t165 + t214;
t22 = t165 * t59 + t166 * t60 + t171 * t69;
t21 = t165 * t57 + t166 * t58 + t171 * t68;
t20 = t155 * t59 + t157 * t60 + t163 * t69;
t19 = t155 * t57 + t157 * t58 + t163 * t68;
t76 = [m(2) + t195; t195 * t182; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t182 ^ 2 + (t179 ^ 2 + t181 ^ 2) * t180 ^ 2); m(4) * t94 + m(5) * t65 + m(6) * t42 + m(7) * t23; m(4) * (t182 * t94 + (t179 * t95 - t181 * t96) * t180) + m(5) * (t182 * t65 + (t179 * t70 - t181 * t71) * t180) + m(6) * (t182 * t42 + (t179 * t48 - t181 * t49) * t180) + m(7) * (t182 * t23 + (t179 * t28 - t181 * t29) * t180); (t18 + t17 + t25 + (t150 * t171 - t151 * t163 + t152 * t164) * t171) * t171 + (t22 + t11 + t12 + (t126 * t166 - t128 * t157 + t130 * t158) * t166 + (t126 * t171 - t128 * t163 + t130 * t164 + t150 * t166 - t151 * t157 + t152 * t158) * t171) * t166 + (t10 + t9 + t21 + (t125 * t165 - t127 * t155 + t129 * t156) * t165 + (t125 * t171 - t127 * t163 + t129 * t164 + t150 * t165 - t151 * t155 + t152 * t156) * t171 + (t125 * t166 + t126 * t165 - t127 * t157 - t128 * t155 + t129 * t158 + t130 * t156) * t166) * t165 + m(5) * (t65 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(7) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t42 ^ 2 + t48 ^ 2 + t49 ^ 2); m(5) * t72 + m(6) * t47 + m(7) * t26; m(5) * (t182 * t72 + (t179 * t74 - t181 * t75) * t180) + m(6) * (t182 * t47 + (t179 * t54 - t181 * t55) * t180) + m(7) * (t182 * t26 + (t179 * t30 - t181 * t31) * t180); (t24 / 0.2e1 + t202) * t171 + (t20 / 0.2e1 + t208) * t166 + (t19 / 0.2e1 + t209) * t165 + (t25 / 0.2e1 + t201) * t163 + (t22 / 0.2e1 + t204) * t157 + (t21 / 0.2e1 + t207) * t155 + m(7) * (t23 * t26 + t28 * t30 + t29 * t31) + m(6) * (t42 * t47 + t48 * t54 + t49 * t55) + m(5) * (t65 * t72 + t70 * t74 + t71 * t75); (t15 + t16 + t24) * t163 + (t7 + t8 + t20) * t157 + (t5 + t6 + t19) * t155 + m(7) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t72 ^ 2 + t74 ^ 2 + t75 ^ 2); m(6) * t56 + m(7) * t27; m(6) * (t182 * t56 + (t179 * t66 - t181 * t67) * t180) + m(7) * (t182 * t27 + (t179 * t32 - t181 * t33) * t180); t203 * t171 + t210 * t166 + t211 * t165 + t201 * t159 + t204 * t146 + t207 * t144 + m(7) * (t23 * t27 + t28 * t32 + t29 * t33) + m(6) * (t42 * t56 + t48 * t66 + t49 * t67); t203 * t163 + t202 * t159 + t210 * t157 + t211 * t155 + t208 * t146 + t209 * t144 + m(7) * (t26 * t27 + t30 * t32 + t31 * t33) + m(6) * (t47 * t56 + t54 * t66 + t55 * t67); (t14 + t13) * t159 + (t3 + t4) * t146 + (t1 + t2) * t144 + m(7) * (t27 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * t159; m(7) * (t159 * t182 + (-t144 * t181 + t146 * t179) * t180); m(7) * (t144 * t29 + t146 * t28 + t159 * t23); m(7) * (t144 * t31 + t146 * t30 + t159 * t26); m(7) * (t144 * t33 + t146 * t32 + t159 * t27); m(7) * (t144 ^ 2 + t146 ^ 2 + t159 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t76(1) t76(2) t76(4) t76(7) t76(11) t76(16); t76(2) t76(3) t76(5) t76(8) t76(12) t76(17); t76(4) t76(5) t76(6) t76(9) t76(13) t76(18); t76(7) t76(8) t76(9) t76(10) t76(14) t76(19); t76(11) t76(12) t76(13) t76(14) t76(15) t76(20); t76(16) t76(17) t76(18) t76(19) t76(20) t76(21);];
Mq  = res;
