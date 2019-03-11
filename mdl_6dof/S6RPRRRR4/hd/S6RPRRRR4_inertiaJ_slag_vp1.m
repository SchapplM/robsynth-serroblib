% Calculate joint inertia matrix for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:48
% EndTime: 2019-03-09 07:04:53
% DurationCPUTime: 2.23s
% Computational Cost: add. (10279->362), mult. (7227->543), div. (0->0), fcn. (7498->12), ass. (0->183)
t171 = sin(qJ(1));
t165 = t171 ^ 2;
t173 = cos(qJ(1));
t164 = pkin(11) + qJ(3);
t157 = qJ(4) + t164;
t151 = sin(t157);
t152 = cos(t157);
t190 = Icges(5,5) * t152 - Icges(5,6) * t151;
t101 = -Icges(5,3) * t173 + t190 * t171;
t102 = Icges(5,3) * t171 + t190 * t173;
t166 = t173 ^ 2;
t225 = Icges(5,4) * t152;
t193 = -Icges(5,2) * t151 + t225;
t104 = Icges(5,6) * t171 + t193 * t173;
t226 = Icges(5,4) * t151;
t196 = Icges(5,1) * t152 - t226;
t106 = Icges(5,5) * t171 + t196 * t173;
t187 = -t104 * t151 + t106 * t152;
t103 = -Icges(5,6) * t173 + t193 * t171;
t105 = -Icges(5,5) * t173 + t196 * t171;
t188 = t103 * t151 - t105 * t152;
t154 = qJ(5) + t157;
t147 = sin(t154);
t148 = cos(t154);
t223 = Icges(6,4) * t148;
t192 = -Icges(6,2) * t147 + t223;
t94 = Icges(6,6) * t171 + t192 * t173;
t224 = Icges(6,4) * t147;
t195 = Icges(6,1) * t148 - t224;
t96 = Icges(6,5) * t171 + t195 * t173;
t198 = -t147 * t94 + t148 * t96;
t93 = -Icges(6,6) * t173 + t192 * t171;
t95 = -Icges(6,5) * t173 + t195 * t171;
t199 = t147 * t93 - t148 * t95;
t172 = cos(qJ(6));
t216 = t173 * t172;
t170 = sin(qJ(6));
t219 = t171 * t170;
t117 = -t148 * t219 - t216;
t217 = t173 * t170;
t218 = t171 * t172;
t118 = t148 * t218 - t217;
t222 = t147 * t171;
t61 = Icges(7,5) * t118 + Icges(7,6) * t117 + Icges(7,3) * t222;
t63 = Icges(7,4) * t118 + Icges(7,2) * t117 + Icges(7,6) * t222;
t65 = Icges(7,1) * t118 + Icges(7,4) * t117 + Icges(7,5) * t222;
t16 = t117 * t63 + t118 * t65 + t61 * t222;
t119 = -t148 * t217 + t218;
t120 = t148 * t216 + t219;
t221 = t147 * t173;
t62 = Icges(7,5) * t120 + Icges(7,6) * t119 + Icges(7,3) * t221;
t64 = Icges(7,4) * t120 + Icges(7,2) * t119 + Icges(7,6) * t221;
t66 = Icges(7,1) * t120 + Icges(7,4) * t119 + Icges(7,5) * t221;
t17 = t117 * t64 + t118 * t66 + t62 * t222;
t8 = -t16 * t173 + t17 * t171;
t189 = Icges(6,5) * t148 - Icges(6,6) * t147;
t91 = -Icges(6,3) * t173 + t189 * t171;
t92 = Icges(6,3) * t171 + t189 * t173;
t243 = -t166 * t91 - (t198 * t171 + (t199 - t92) * t173) * t171 - t8;
t248 = -t166 * t101 - (t187 * t171 + (-t102 + t188) * t173) * t171 + t243;
t220 = t148 * t173;
t68 = t120 * rSges(7,1) + t119 * rSges(7,2) + rSges(7,3) * t221;
t247 = pkin(5) * t220 + pkin(10) * t221 + t68;
t202 = rSges(5,1) * t152 - rSges(5,2) * t151;
t246 = t171 / 0.2e1;
t245 = -t173 / 0.2e1;
t18 = t119 * t63 + t120 * t65 + t61 * t221;
t19 = t119 * t64 + t120 * t66 + t62 * t221;
t9 = t19 * t171 - t18 * t173;
t244 = (t165 * t92 + t9 + (t199 * t173 + (t198 - t91) * t171) * t173) * t171;
t155 = sin(t164);
t242 = pkin(3) * t155;
t241 = pkin(4) * t151;
t240 = pkin(5) * t148;
t169 = -pkin(7) - qJ(2);
t168 = cos(pkin(11));
t153 = t168 * pkin(2) + pkin(1);
t156 = cos(t164);
t138 = pkin(3) * t156 + t153;
t116 = pkin(4) * t152 + t138;
t113 = t173 * t116;
t132 = t173 * t138;
t239 = t173 * (t113 - t132) + (t116 - t138) * t165;
t238 = t173 * (-t173 * t153 + t132) + (t138 - t153) * t165;
t180 = rSges(6,1) * t220 - rSges(6,2) * t221 + t171 * rSges(6,3);
t201 = rSges(6,1) * t148 - rSges(6,2) * t147;
t51 = t171 * (-t173 * rSges(6,3) + t201 * t171) + t173 * t180;
t181 = t171 * rSges(5,3) + t202 * t173;
t54 = t171 * (-t173 * rSges(5,3) + t202 * t171) + t173 * t181;
t237 = rSges(4,1) * t156;
t235 = rSges(4,2) * t155;
t80 = -Icges(7,6) * t148 + (Icges(7,4) * t172 - Icges(7,2) * t170) * t147;
t233 = t170 * t80;
t24 = -t148 * t61 + (-t170 * t63 + t172 * t65) * t147;
t232 = t24 * t173;
t25 = -t148 * t62 + (-t170 * t64 + t172 * t66) * t147;
t231 = t25 * t171;
t230 = rSges(3,3) + qJ(2);
t82 = -t148 * rSges(7,3) + (rSges(7,1) * t172 - rSges(7,2) * t170) * t147;
t229 = -t147 * pkin(5) + t148 * pkin(10) - t82;
t228 = Icges(4,4) * t155;
t227 = Icges(4,4) * t156;
t214 = t171 * rSges(4,3) + t173 * t237;
t211 = t165 + t166;
t163 = -pkin(8) + t169;
t210 = t171 * (t165 * t102 + (t188 * t173 + (-t101 + t187) * t171) * t173) + t244;
t131 = t151 * rSges(5,1) + t152 * rSges(5,2);
t209 = -t131 - t242;
t125 = t147 * rSges(6,1) + t148 * rSges(6,2);
t208 = -t125 - t241;
t26 = t51 + t239;
t200 = -t118 * rSges(7,1) - t117 * rSges(7,2);
t67 = rSges(7,3) * t222 - t200;
t27 = t171 * t67 + t165 * (pkin(10) * t147 + t240) + t247 * t173;
t158 = -pkin(9) + t163;
t207 = -t171 * t158 + t113;
t79 = -Icges(7,3) * t148 + (Icges(7,5) * t172 - Icges(7,6) * t170) * t147;
t81 = -Icges(7,5) * t148 + (Icges(7,1) * t172 - Icges(7,4) * t170) * t147;
t30 = t117 * t80 + t118 * t81 + t79 * t222;
t3 = -t30 * t148 + (t16 * t171 + t17 * t173) * t147;
t31 = t119 * t80 + t120 * t81 + t79 * t221;
t4 = -t31 * t148 + (t171 * t18 + t173 * t19) * t147;
t206 = t3 * t245 + t4 * t246 - t148 * (t231 - t232) / 0.2e1 + t8 * t222 / 0.2e1 + t9 * t221 / 0.2e1;
t205 = t229 - t241;
t204 = -t241 - t242;
t203 = -t235 + t237;
t197 = Icges(4,1) * t156 - t228;
t194 = -Icges(4,2) * t155 + t227;
t191 = Icges(4,5) * t156 - Icges(4,6) * t155;
t123 = Icges(6,2) * t148 + t224;
t124 = Icges(6,1) * t147 + t223;
t184 = -t123 * t147 + t124 * t148;
t129 = Icges(5,2) * t152 + t226;
t130 = Icges(5,1) * t151 + t225;
t183 = -t129 * t151 + t130 * t152;
t182 = t243 * t173 + t244;
t12 = t27 + t239;
t167 = sin(pkin(11));
t179 = rSges(3,1) * t168 - rSges(3,2) * t167 + pkin(1);
t178 = -t125 + t204;
t177 = t204 + t229;
t122 = Icges(6,5) * t147 + Icges(6,6) * t148;
t176 = -t232 / 0.2e1 + t231 / 0.2e1 + (t171 * t122 + t147 * t96 + t148 * t94 + t184 * t173 + t31) * t246 + (-t173 * t122 + t147 * t95 + t148 * t93 + t184 * t171 + t30) * t245;
t175 = t248 * t173 + t210;
t128 = Icges(5,5) * t151 + Icges(5,6) * t152;
t174 = t176 + (t152 * t104 + t151 * t106 + t171 * t128 + t183 * t173) * t246 + (t152 * t103 + t151 * t105 - t173 * t128 + t183 * t171) * t245;
t145 = t173 * rSges(2,1) - t171 * rSges(2,2);
t144 = -t171 * rSges(2,1) - t173 * rSges(2,2);
t137 = t155 * rSges(4,1) + t156 * rSges(4,2);
t108 = Icges(4,3) * t171 + t191 * t173;
t107 = -Icges(4,3) * t173 + t191 * t171;
t98 = t230 * t171 + t179 * t173;
t97 = -t179 * t171 + t230 * t173;
t88 = t209 * t173;
t87 = t209 * t171;
t84 = t208 * t173;
t83 = t208 * t171;
t78 = -t171 * t169 + (t153 - t235) * t173 + t214;
t77 = (rSges(4,3) - t169) * t173 + (-t153 - t203) * t171;
t76 = t178 * t173;
t75 = t178 * t171;
t72 = t147 * t172 * t81;
t71 = -t171 * t163 + t132 + t181;
t70 = (rSges(5,3) - t163) * t173 + (-t138 - t202) * t171;
t69 = t173 * (-t173 * t235 + t214) + (-t173 * rSges(4,3) + t203 * t171) * t171;
t56 = t229 * t173;
t55 = t229 * t171;
t53 = t180 + t207;
t52 = (rSges(6,3) - t158) * t173 + (-t116 - t201) * t171;
t50 = t205 * t173;
t49 = t205 * t171;
t44 = t177 * t173;
t43 = t177 * t171;
t38 = -t148 * t68 - t82 * t221;
t37 = t148 * t67 + t82 * t222;
t36 = t207 + t247;
t35 = -t173 * t158 + (-t240 - t116 + (-rSges(7,3) - pkin(10)) * t147) * t171 + t200;
t34 = -t147 * t233 - t148 * t79 + t72;
t33 = (-t171 * t68 + t173 * t67) * t147;
t32 = t54 + t238;
t13 = t26 + t238;
t11 = t12 + t238;
t1 = [Icges(3,2) * t168 ^ 2 + t152 * t129 + t151 * t130 + t156 * (Icges(4,2) * t156 + t228) + t155 * (Icges(4,1) * t155 + t227) + Icges(2,3) + t72 + (Icges(3,1) * t167 + 0.2e1 * Icges(3,4) * t168) * t167 + (-t79 + t123) * t148 + (t124 - t233) * t147 + m(7) * (t35 ^ 2 + t36 ^ 2) + m(6) * (t52 ^ 2 + t53 ^ 2) + m(5) * (t70 ^ 2 + t71 ^ 2) + m(4) * (t77 ^ 2 + t78 ^ 2) + m(3) * (t97 ^ 2 + t98 ^ 2) + m(2) * (t144 ^ 2 + t145 ^ 2); m(7) * (t171 * t35 - t173 * t36) + m(6) * (t171 * t52 - t173 * t53) + m(5) * (t171 * t70 - t173 * t71) + m(4) * (t171 * t77 - t173 * t78) + m(3) * (t171 * t97 - t173 * t98); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t211; m(4) * (-t171 * t78 - t173 * t77) * t137 + m(7) * (t44 * t35 + t43 * t36) + m(6) * (t76 * t52 + t75 * t53) + m(5) * (t88 * t70 + t87 * t71) + t174 + (t166 / 0.2e1 + t165 / 0.2e1) * (Icges(4,5) * t155 + Icges(4,6) * t156) + (t156 * (Icges(4,6) * t171 + t194 * t173) + t155 * (Icges(4,5) * t171 + t197 * t173)) * t246 + (t156 * (-Icges(4,6) * t173 + t194 * t171) + t155 * (-Icges(4,5) * t173 + t197 * t171)) * t245; m(5) * (t88 * t171 - t87 * t173) + m(6) * (t76 * t171 - t75 * t173) + m(7) * (t44 * t171 - t43 * t173); m(7) * (t11 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t13 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t32 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t137 ^ 2 * t211 + t69 ^ 2) + t171 * t165 * t108 + t210 + (-t166 * t107 + (-t171 * t107 + t173 * t108) * t171 + t248) * t173; m(7) * (t50 * t35 + t49 * t36) + m(6) * (t84 * t52 + t83 * t53) + m(5) * (-t171 * t71 - t173 * t70) * t131 + t174; m(6) * (t84 * t171 - t83 * t173) + m(7) * (t50 * t171 - t49 * t173); m(7) * (t12 * t11 + t49 * t43 + t50 * t44) + m(6) * (t26 * t13 + t83 * t75 + t84 * t76) + m(5) * (t54 * t32 + (-t171 * t87 - t173 * t88) * t131) + t175; m(7) * (t12 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t26 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(5) * (t131 ^ 2 * t211 + t54 ^ 2) + t175; m(7) * (t56 * t35 + t55 * t36) + m(6) * (-t171 * t53 - t173 * t52) * t125 + t176; m(7) * (t56 * t171 - t55 * t173); m(7) * (t27 * t11 + t55 * t43 + t56 * t44) + m(6) * (t51 * t13 + (-t171 * t75 - t173 * t76) * t125) + t182; m(7) * (t27 * t12 + t55 * t49 + t56 * t50) + m(6) * (t51 * t26 + (-t171 * t83 - t173 * t84) * t125) + t182; m(6) * (t125 ^ 2 * t211 + t51 ^ 2) + m(7) * (t27 ^ 2 + t55 ^ 2 + t56 ^ 2) + t182; m(7) * (t37 * t35 + t38 * t36) - t34 * t148 + ((t25 / 0.2e1 + t31 / 0.2e1) * t173 + (t24 / 0.2e1 + t30 / 0.2e1) * t171) * t147; m(7) * (t37 * t171 - t38 * t173); m(7) * (t33 * t11 + t37 * t44 + t38 * t43) + t206; m(7) * (t33 * t12 + t37 * t50 + t38 * t49) + t206; m(7) * (t33 * t27 + t37 * t56 + t38 * t55) + t206; t148 ^ 2 * t34 + m(7) * (t33 ^ 2 + t37 ^ 2 + t38 ^ 2) + (t173 * t4 + t171 * t3 - t148 * (t171 * t24 + t173 * t25)) * t147;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
