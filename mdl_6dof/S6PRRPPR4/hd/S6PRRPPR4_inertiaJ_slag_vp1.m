% Calculate joint inertia matrix for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:34
% EndTime: 2019-03-08 21:12:43
% DurationCPUTime: 5.46s
% Computational Cost: add. (16680->536), mult. (43273->778), div. (0->0), fcn. (56310->12), ass. (0->239)
t263 = m(6) + m(7);
t264 = m(5) + t263;
t215 = cos(pkin(6));
t212 = sin(pkin(10));
t220 = cos(qJ(2));
t248 = t220 * t212;
t214 = cos(pkin(10));
t218 = sin(qJ(2));
t250 = t214 * t218;
t202 = t215 * t250 + t248;
t217 = sin(qJ(3));
t213 = sin(pkin(6));
t258 = cos(qJ(3));
t231 = t213 * t258;
t191 = t202 * t217 + t214 * t231;
t247 = t220 * t214;
t249 = t218 * t212;
t204 = -t215 * t249 + t247;
t193 = t204 * t217 - t212 * t231;
t252 = t213 * t217;
t205 = -t215 * t258 + t218 * t252;
t194 = t204 * t258 + t212 * t252;
t203 = t215 * t248 + t250;
t211 = sin(pkin(11));
t255 = cos(pkin(11));
t167 = t194 * t211 - t203 * t255;
t168 = t194 * t255 + t203 * t211;
t216 = sin(qJ(6));
t219 = cos(qJ(6));
t128 = t167 * t219 - t168 * t216;
t129 = t167 * t216 + t168 * t219;
t192 = t202 * t258 - t214 * t252;
t201 = -t215 * t247 + t249;
t165 = t192 * t211 - t201 * t255;
t166 = t192 * t255 + t201 * t211;
t126 = t165 * t219 - t166 * t216;
t127 = t165 * t216 + t166 * t219;
t82 = Icges(7,5) * t127 + Icges(7,6) * t126 - Icges(7,3) * t191;
t84 = Icges(7,4) * t127 + Icges(7,2) * t126 - Icges(7,6) * t191;
t86 = Icges(7,1) * t127 + Icges(7,4) * t126 - Icges(7,5) * t191;
t26 = t128 * t84 + t129 * t86 - t193 * t82;
t83 = Icges(7,5) * t129 + Icges(7,6) * t128 - Icges(7,3) * t193;
t85 = Icges(7,4) * t129 + Icges(7,2) * t128 - Icges(7,6) * t193;
t87 = Icges(7,1) * t129 + Icges(7,4) * t128 - Icges(7,5) * t193;
t27 = t128 * t85 + t129 * t87 - t193 * t83;
t206 = t215 * t217 + t218 * t231;
t251 = t213 * t220;
t189 = t206 * t211 + t251 * t255;
t190 = t206 * t255 - t211 * t251;
t160 = t189 * t219 - t190 * t216;
t161 = t189 * t216 + t190 * t219;
t100 = Icges(7,1) * t161 + Icges(7,4) * t160 - Icges(7,5) * t205;
t98 = Icges(7,5) * t161 + Icges(7,6) * t160 - Icges(7,3) * t205;
t99 = Icges(7,4) * t161 + Icges(7,2) * t160 - Icges(7,6) * t205;
t43 = t100 * t129 + t128 * t99 - t193 * t98;
t2 = -t191 * t26 - t193 * t27 - t205 * t43;
t262 = t2 / 0.2e1;
t261 = -t191 / 0.2e1;
t260 = -t193 / 0.2e1;
t259 = -t205 / 0.2e1;
t88 = rSges(7,1) * t127 + rSges(7,2) * t126 - rSges(7,3) * t191;
t257 = pkin(5) * t166 - pkin(9) * t191 + t88;
t89 = rSges(7,1) * t129 + rSges(7,2) * t128 - rSges(7,3) * t193;
t256 = pkin(5) * t168 - pkin(9) * t193 + t89;
t254 = t212 * t213;
t253 = t213 * t214;
t101 = rSges(7,1) * t161 + rSges(7,2) * t160 - rSges(7,3) * t205;
t246 = pkin(5) * t190 - pkin(9) * t205 + t101;
t119 = rSges(5,1) * t168 - rSges(5,2) * t167 + rSges(5,3) * t193;
t164 = pkin(3) * t194 + qJ(4) * t193;
t245 = -t119 - t164;
t130 = pkin(4) * t166 + qJ(5) * t165;
t163 = pkin(3) * t192 + qJ(4) * t191;
t152 = t203 * t163;
t244 = t203 * t130 + t152;
t131 = pkin(4) * t168 + qJ(5) * t167;
t243 = -t131 - t164;
t150 = rSges(5,1) * t190 - rSges(5,2) * t189 + rSges(5,3) * t205;
t188 = pkin(3) * t206 + qJ(4) * t205;
t242 = -t150 - t188;
t241 = t163 * t251 + t201 * t188;
t187 = pkin(2) * t204 + pkin(8) * t203;
t185 = t215 * t187;
t240 = t215 * t164 + t185;
t162 = pkin(4) * t190 + qJ(5) * t189;
t239 = -t162 - t188;
t186 = pkin(2) * t202 + pkin(8) * t201;
t238 = -t163 - t186;
t237 = t186 * t254 + t187 * t253;
t118 = rSges(6,1) * t168 + rSges(6,2) * t193 + rSges(6,3) * t167;
t235 = -t118 + t243;
t234 = t215 * t131 + t240;
t233 = -t130 + t238;
t149 = rSges(6,1) * t190 + rSges(6,2) * t205 + rSges(6,3) * t189;
t232 = -t149 + t239;
t182 = t206 * rSges(4,1) - t205 * rSges(4,2) - rSges(4,3) * t251;
t207 = (pkin(2) * t218 - pkin(8) * t220) * t213;
t230 = (-t182 - t207) * t213;
t228 = t243 - t256;
t227 = t239 - t246;
t226 = t130 * t251 + t201 * t162 + t241;
t225 = t163 * t254 + t164 * t253 + t237;
t224 = (-t207 + t242) * t213;
t223 = (-t207 + t232) * t213;
t222 = t130 * t254 + t131 * t253 + t225;
t221 = (-t207 + t227) * t213;
t198 = t215 * rSges(3,3) + (rSges(3,1) * t218 + rSges(3,2) * t220) * t213;
t197 = Icges(3,5) * t215 + (Icges(3,1) * t218 + Icges(3,4) * t220) * t213;
t196 = Icges(3,6) * t215 + (Icges(3,4) * t218 + Icges(3,2) * t220) * t213;
t195 = Icges(3,3) * t215 + (Icges(3,5) * t218 + Icges(3,6) * t220) * t213;
t181 = Icges(4,1) * t206 - Icges(4,4) * t205 - Icges(4,5) * t251;
t180 = Icges(4,4) * t206 - Icges(4,2) * t205 - Icges(4,6) * t251;
t179 = Icges(4,5) * t206 - Icges(4,6) * t205 - Icges(4,3) * t251;
t178 = rSges(3,1) * t204 - rSges(3,2) * t203 + rSges(3,3) * t254;
t177 = rSges(3,1) * t202 - rSges(3,2) * t201 - rSges(3,3) * t253;
t176 = Icges(3,1) * t204 - Icges(3,4) * t203 + Icges(3,5) * t254;
t175 = Icges(3,1) * t202 - Icges(3,4) * t201 - Icges(3,5) * t253;
t174 = Icges(3,4) * t204 - Icges(3,2) * t203 + Icges(3,6) * t254;
t173 = Icges(3,4) * t202 - Icges(3,2) * t201 - Icges(3,6) * t253;
t172 = Icges(3,5) * t204 - Icges(3,6) * t203 + Icges(3,3) * t254;
t171 = Icges(3,5) * t202 - Icges(3,6) * t201 - Icges(3,3) * t253;
t155 = -t177 * t215 - t198 * t253;
t154 = t178 * t215 - t198 * t254;
t148 = rSges(4,1) * t194 - rSges(4,2) * t193 + rSges(4,3) * t203;
t147 = rSges(4,1) * t192 - rSges(4,2) * t191 + rSges(4,3) * t201;
t146 = Icges(5,1) * t190 - Icges(5,4) * t189 + Icges(5,5) * t205;
t145 = Icges(6,1) * t190 + Icges(6,4) * t205 + Icges(6,5) * t189;
t144 = Icges(5,4) * t190 - Icges(5,2) * t189 + Icges(5,6) * t205;
t143 = Icges(6,4) * t190 + Icges(6,2) * t205 + Icges(6,6) * t189;
t142 = Icges(5,5) * t190 - Icges(5,6) * t189 + Icges(5,3) * t205;
t141 = Icges(6,5) * t190 + Icges(6,6) * t205 + Icges(6,3) * t189;
t140 = Icges(4,1) * t194 - Icges(4,4) * t193 + Icges(4,5) * t203;
t139 = Icges(4,1) * t192 - Icges(4,4) * t191 + Icges(4,5) * t201;
t138 = Icges(4,4) * t194 - Icges(4,2) * t193 + Icges(4,6) * t203;
t137 = Icges(4,4) * t192 - Icges(4,2) * t191 + Icges(4,6) * t201;
t136 = Icges(4,5) * t194 - Icges(4,6) * t193 + Icges(4,3) * t203;
t135 = Icges(4,5) * t192 - Icges(4,6) * t191 + Icges(4,3) * t201;
t132 = (t177 * t212 + t178 * t214) * t213;
t117 = rSges(5,1) * t166 - rSges(5,2) * t165 + rSges(5,3) * t191;
t116 = rSges(6,1) * t166 + rSges(6,2) * t191 + rSges(6,3) * t165;
t115 = Icges(5,1) * t168 - Icges(5,4) * t167 + Icges(5,5) * t193;
t114 = Icges(5,1) * t166 - Icges(5,4) * t165 + Icges(5,5) * t191;
t113 = Icges(6,1) * t168 + Icges(6,4) * t193 + Icges(6,5) * t167;
t112 = Icges(6,1) * t166 + Icges(6,4) * t191 + Icges(6,5) * t165;
t111 = Icges(5,4) * t168 - Icges(5,2) * t167 + Icges(5,6) * t193;
t110 = Icges(5,4) * t166 - Icges(5,2) * t165 + Icges(5,6) * t191;
t109 = Icges(6,4) * t168 + Icges(6,2) * t193 + Icges(6,6) * t167;
t108 = Icges(6,4) * t166 + Icges(6,2) * t191 + Icges(6,6) * t165;
t107 = Icges(5,5) * t168 - Icges(5,6) * t167 + Icges(5,3) * t193;
t106 = Icges(5,5) * t166 - Icges(5,6) * t165 + Icges(5,3) * t191;
t105 = Icges(6,5) * t168 + Icges(6,6) * t193 + Icges(6,3) * t167;
t104 = Icges(6,5) * t166 + Icges(6,6) * t191 + Icges(6,3) * t165;
t103 = -t148 * t251 - t203 * t182;
t102 = t147 * t251 + t201 * t182;
t97 = -t179 * t251 - t205 * t180 + t206 * t181;
t96 = t147 * t203 - t148 * t201;
t95 = (-t147 - t186) * t215 + t214 * t230;
t94 = t148 * t215 + t212 * t230 + t185;
t93 = t179 * t203 - t180 * t193 + t181 * t194;
t92 = t179 * t201 - t180 * t191 + t181 * t192;
t90 = (t147 * t212 + t148 * t214) * t213 + t237;
t81 = -t136 * t251 - t205 * t138 + t206 * t140;
t80 = -t135 * t251 - t205 * t137 + t206 * t139;
t79 = t142 * t205 - t144 * t189 + t146 * t190;
t78 = t141 * t189 + t143 * t205 + t145 * t190;
t77 = t136 * t203 - t138 * t193 + t140 * t194;
t76 = t135 * t203 - t137 * t193 + t139 * t194;
t75 = t136 * t201 - t138 * t191 + t140 * t192;
t74 = t135 * t201 - t137 * t191 + t139 * t192;
t73 = t203 * t242 + t245 * t251;
t72 = t117 * t251 + t201 * t150 + t241;
t71 = t142 * t193 - t144 * t167 + t146 * t168;
t70 = t141 * t167 + t143 * t193 + t145 * t168;
t69 = t142 * t191 - t144 * t165 + t146 * t166;
t68 = t141 * t165 + t143 * t191 + t145 * t166;
t67 = (-t117 + t238) * t215 + t214 * t224;
t66 = t119 * t215 + t212 * t224 + t240;
t65 = t101 * t193 - t205 * t89;
t64 = -t101 * t191 + t205 * t88;
t63 = t117 * t203 + t201 * t245 + t152;
t62 = t107 * t205 - t111 * t189 + t115 * t190;
t61 = t106 * t205 - t110 * t189 + t114 * t190;
t60 = t105 * t189 + t109 * t205 + t113 * t190;
t59 = t104 * t189 + t108 * t205 + t112 * t190;
t58 = (t117 * t212 + t119 * t214) * t213 + t225;
t57 = t203 * t232 + t235 * t251;
t56 = t116 * t251 + t201 * t149 + t226;
t55 = t107 * t193 - t111 * t167 + t115 * t168;
t54 = t106 * t193 - t110 * t167 + t114 * t168;
t53 = t105 * t167 + t109 * t193 + t113 * t168;
t52 = t104 * t167 + t108 * t193 + t112 * t168;
t51 = t107 * t191 - t111 * t165 + t115 * t166;
t50 = t106 * t191 - t110 * t165 + t114 * t166;
t49 = t105 * t165 + t109 * t191 + t113 * t166;
t48 = t104 * t165 + t108 * t191 + t112 * t166;
t47 = (-t116 + t233) * t215 + t214 * t223;
t46 = t118 * t215 + t212 * t223 + t234;
t45 = t191 * t89 - t193 * t88;
t44 = t100 * t161 + t160 * t99 - t205 * t98;
t42 = t100 * t127 + t126 * t99 - t191 * t98;
t41 = t116 * t203 + t201 * t235 + t244;
t40 = t215 * t97 + (t212 * t81 - t214 * t80) * t213;
t39 = (t116 * t212 + t118 * t214) * t213 + t222;
t38 = t80 * t201 + t81 * t203 - t251 * t97;
t37 = t215 * t93 + (t212 * t77 - t214 * t76) * t213;
t36 = t215 * t92 + (t212 * t75 - t214 * t74) * t213;
t35 = t76 * t201 + t77 * t203 - t251 * t93;
t34 = t74 * t201 + t75 * t203 - t251 * t92;
t33 = t203 * t227 + t228 * t251;
t32 = t201 * t246 + t251 * t257 + t226;
t31 = (t233 - t257) * t215 + t214 * t221;
t30 = t212 * t221 + t215 * t256 + t234;
t29 = t160 * t85 + t161 * t87 - t205 * t83;
t28 = t160 * t84 + t161 * t86 - t205 * t82;
t25 = t126 * t85 + t127 * t87 - t191 * t83;
t24 = t126 * t84 + t127 * t86 - t191 * t82;
t23 = t201 * t228 + t203 * t257 + t244;
t22 = (t212 * t257 + t214 * t256) * t213 + t222;
t21 = t215 * t79 + (t212 * t62 - t214 * t61) * t213;
t20 = t215 * t78 + (t212 * t60 - t214 * t59) * t213;
t19 = t61 * t201 + t62 * t203 - t251 * t79;
t18 = t59 * t201 + t60 * t203 - t251 * t78;
t17 = t215 * t71 + (t212 * t55 - t214 * t54) * t213;
t16 = t215 * t70 + (t212 * t53 - t214 * t52) * t213;
t15 = t215 * t69 + (t212 * t51 - t214 * t50) * t213;
t14 = t215 * t68 + (t212 * t49 - t214 * t48) * t213;
t13 = t54 * t201 + t55 * t203 - t251 * t71;
t12 = t52 * t201 + t53 * t203 - t251 * t70;
t11 = t50 * t201 + t51 * t203 - t251 * t69;
t10 = t48 * t201 + t49 * t203 - t251 * t68;
t9 = t215 * t44 + (t212 * t29 - t214 * t28) * t213;
t8 = t28 * t201 + t29 * t203 - t251 * t44;
t7 = -t191 * t28 - t193 * t29 - t205 * t44;
t6 = t215 * t43 + (t212 * t27 - t214 * t26) * t213;
t5 = t215 * t42 + (t212 * t25 - t214 * t24) * t213;
t4 = t26 * t201 + t27 * t203 - t251 * t43;
t3 = t24 * t201 + t25 * t203 - t251 * t42;
t1 = -t191 * t24 - t193 * t25 - t205 * t42;
t91 = [m(4) + m(2) + m(3) + t264; m(3) * t132 + m(4) * t90 + m(5) * t58 + m(6) * t39 + m(7) * t22; m(7) * (t22 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t39 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t58 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(4) * (t90 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(3) * (t132 ^ 2 + t154 ^ 2 + t155 ^ 2) + (t6 + t37 + t17 + t16 + (t172 * t254 - t174 * t203 + t176 * t204) * t254) * t254 + (-t5 - t36 - t15 - t14 + (-t171 * t253 - t173 * t201 + t175 * t202) * t253 + (-t171 * t254 + t172 * t253 + t173 * t203 + t174 * t201 - t175 * t204 - t176 * t202) * t254) * t253 + (-(-t195 * t253 - t196 * t201 + t197 * t202) * t253 + (t195 * t254 - t196 * t203 + t197 * t204) * t254 + t9 + t40 + t21 + t20 + ((t174 * t220 + t176 * t218) * t212 - (t173 * t220 + t175 * t218) * t214) * t213 ^ 2 + ((-t171 * t214 + t172 * t212 + t196 * t220 + t197 * t218) * t213 + t215 * t195) * t215) * t215; m(4) * t96 + m(5) * t63 + m(6) * t41 + m(7) * t23; (t8 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1 + t38 / 0.2e1) * t215 + (t6 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1 + t37 / 0.2e1) * t203 + (t5 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1 + t36 / 0.2e1) * t201 + m(7) * (t22 * t23 + t30 * t33 + t31 * t32) + m(6) * (t39 * t41 + t46 * t57 + t47 * t56) + m(5) * (t58 * t63 + t66 * t73 + t67 * t72) + m(4) * (t102 * t95 + t103 * t94 + t90 * t96) + ((-t20 / 0.2e1 - t21 / 0.2e1 - t40 / 0.2e1 - t9 / 0.2e1) * t220 + (-t3 / 0.2e1 - t34 / 0.2e1 - t10 / 0.2e1 - t11 / 0.2e1) * t214 + (t4 / 0.2e1 + t35 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1) * t212) * t213; (-t18 - t19 - t38 - t8) * t251 + (t4 + t13 + t12 + t35) * t203 + (t3 + t11 + t10 + t34) * t201 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t41 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t63 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2 + t96 ^ 2); t205 * t264; m(7) * (t191 * t30 + t193 * t31 + t205 * t22) + m(6) * (t191 * t46 + t193 * t47 + t205 * t39) + m(5) * (t191 * t66 + t193 * t67 + t205 * t58); m(7) * (t191 * t33 + t193 * t32 + t205 * t23) + m(6) * (t191 * t57 + t193 * t56 + t205 * t41) + m(5) * (t191 * t73 + t193 * t72 + t205 * t63); (t191 ^ 2 + t193 ^ 2 + t205 ^ 2) * t264; t189 * t263; m(7) * (t165 * t30 + t167 * t31 + t189 * t22) + m(6) * (t165 * t46 + t167 * t47 + t189 * t39); m(7) * (t165 * t33 + t167 * t32 + t189 * t23) + m(6) * (t165 * t57 + t167 * t56 + t189 * t41); (t165 * t191 + t167 * t193 + t189 * t205) * t263; (t165 ^ 2 + t167 ^ 2 + t189 ^ 2) * t263; m(7) * t45; t9 * t259 + t6 * t260 + t5 * t261 + t215 * t7 / 0.2e1 + m(7) * (t22 * t45 + t30 * t65 + t31 * t64) + (t212 * t262 - t214 * t1 / 0.2e1) * t213; -t7 * t251 / 0.2e1 + m(7) * (t23 * t45 + t32 * t64 + t33 * t65) + t3 * t261 + t4 * t260 + t8 * t259 + t203 * t262 + t201 * t1 / 0.2e1; m(7) * (t191 * t65 + t193 * t64 + t205 * t45); m(7) * (t165 * t65 + t167 * t64 + t189 * t45); -t193 * t2 - t191 * t1 - t205 * t7 + m(7) * (t45 ^ 2 + t64 ^ 2 + t65 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t91(1) t91(2) t91(4) t91(7) t91(11) t91(16); t91(2) t91(3) t91(5) t91(8) t91(12) t91(17); t91(4) t91(5) t91(6) t91(9) t91(13) t91(18); t91(7) t91(8) t91(9) t91(10) t91(14) t91(19); t91(11) t91(12) t91(13) t91(14) t91(15) t91(20); t91(16) t91(17) t91(18) t91(19) t91(20) t91(21);];
Mq  = res;
