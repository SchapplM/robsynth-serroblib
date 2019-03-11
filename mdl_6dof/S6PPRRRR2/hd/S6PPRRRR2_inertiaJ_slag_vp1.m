% Calculate joint inertia matrix for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:19
% EndTime: 2019-03-08 19:03:30
% DurationCPUTime: 5.32s
% Computational Cost: add. (46298->504), mult. (120527->723), div. (0->0), fcn. (159630->16), ass. (0->234)
t209 = sin(pkin(12));
t211 = cos(pkin(12));
t212 = cos(pkin(6));
t246 = sin(pkin(13));
t228 = t212 * t246;
t248 = cos(pkin(13));
t199 = t209 * t248 + t211 * t228;
t215 = sin(qJ(3));
t229 = t212 * t248;
t221 = t209 * t246 - t211 * t229;
t249 = cos(pkin(7));
t219 = t221 * t249;
t210 = sin(pkin(6));
t247 = sin(pkin(7));
t230 = t210 * t247;
t255 = cos(qJ(3));
t183 = t199 * t255 + (-t211 * t230 - t219) * t215;
t231 = t210 * t249;
t192 = -t211 * t231 + t221 * t247;
t214 = sin(qJ(4));
t254 = cos(qJ(4));
t171 = t183 * t214 - t192 * t254;
t264 = t171 / 0.2e1;
t200 = -t209 * t228 + t211 * t248;
t220 = t209 * t229 + t211 * t246;
t218 = t220 * t249;
t185 = t200 * t255 + (t209 * t230 - t218) * t215;
t193 = t209 * t231 + t220 * t247;
t173 = t185 * t214 - t193 * t254;
t263 = t173 / 0.2e1;
t226 = t255 * t247;
t222 = t210 * t226;
t182 = t199 * t215 + t211 * t222 + t255 * t219;
t262 = t182 / 0.2e1;
t184 = t200 * t215 - t209 * t222 + t255 * t218;
t261 = t184 / 0.2e1;
t224 = t249 * t248;
t191 = t212 * t247 * t215 + (t215 * t224 + t246 * t255) * t210;
t198 = t212 * t249 - t248 * t230;
t186 = t191 * t214 - t198 * t254;
t260 = t186 / 0.2e1;
t190 = -t212 * t226 + (t215 * t246 - t224 * t255) * t210;
t259 = t190 / 0.2e1;
t258 = t192 / 0.2e1;
t257 = t193 / 0.2e1;
t256 = t198 / 0.2e1;
t216 = cos(qJ(5));
t253 = pkin(5) * t216;
t172 = t183 * t254 + t192 * t214;
t208 = qJ(5) + qJ(6);
t205 = sin(t208);
t206 = cos(t208);
t142 = -t172 * t205 + t182 * t206;
t143 = t172 * t206 + t182 * t205;
t104 = rSges(7,1) * t143 + rSges(7,2) * t142 + rSges(7,3) * t171;
t213 = sin(qJ(5));
t245 = t182 * t213;
t96 = pkin(5) * t245 + pkin(11) * t171 + t253 * t172;
t251 = t104 + t96;
t174 = t185 * t254 + t193 * t214;
t144 = -t174 * t205 + t184 * t206;
t145 = t174 * t206 + t184 * t205;
t105 = rSges(7,1) * t145 + rSges(7,2) * t144 + rSges(7,3) * t173;
t244 = t184 * t213;
t97 = pkin(5) * t244 + pkin(11) * t173 + t253 * t174;
t250 = t105 + t97;
t243 = t190 * t213;
t146 = -t172 * t213 + t182 * t216;
t147 = t172 * t216 + t245;
t112 = rSges(6,1) * t147 + rSges(6,2) * t146 + rSges(6,3) * t171;
t139 = t172 * pkin(4) + t171 * pkin(10);
t242 = -t112 - t139;
t148 = -t174 * t213 + t184 * t216;
t149 = t174 * t216 + t244;
t113 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t173;
t140 = t174 * pkin(4) + t173 * pkin(10);
t241 = -t113 - t140;
t187 = t191 * t254 + t198 * t214;
t118 = pkin(5) * t243 + pkin(11) * t186 + t253 * t187;
t168 = -t187 * t205 + t190 * t206;
t169 = t187 * t206 + t190 * t205;
t122 = rSges(7,1) * t169 + rSges(7,2) * t168 + rSges(7,3) * t186;
t240 = t118 + t122;
t175 = -t187 * t213 + t190 * t216;
t176 = t187 * t216 + t243;
t134 = rSges(6,1) * t176 + rSges(6,2) * t175 + rSges(6,3) * t186;
t167 = t187 * pkin(4) + t186 * pkin(10);
t239 = -t134 - t167;
t165 = pkin(3) * t183 + pkin(9) * t182;
t162 = t193 * t165;
t238 = t193 * t139 + t162;
t166 = pkin(3) * t185 + pkin(9) * t184;
t164 = t198 * t166;
t237 = t198 * t140 + t164;
t181 = pkin(3) * t191 + pkin(9) * t190;
t170 = t192 * t181;
t236 = t192 * t167 + t170;
t100 = Icges(7,4) * t143 + Icges(7,2) * t142 + Icges(7,6) * t171;
t102 = Icges(7,1) * t143 + Icges(7,4) * t142 + Icges(7,5) * t171;
t98 = Icges(7,5) * t143 + Icges(7,6) * t142 + Icges(7,3) * t171;
t58 = t100 * t168 + t102 * t169 + t186 * t98;
t101 = Icges(7,4) * t145 + Icges(7,2) * t144 + Icges(7,6) * t173;
t103 = Icges(7,1) * t145 + Icges(7,4) * t144 + Icges(7,5) * t173;
t99 = Icges(7,5) * t145 + Icges(7,6) * t144 + Icges(7,3) * t173;
t59 = t101 * t168 + t103 * t169 + t186 * t99;
t119 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t186;
t120 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t186;
t121 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t186;
t72 = t119 * t186 + t120 * t168 + t121 * t169;
t26 = t171 * t58 + t173 * t59 + t186 * t72;
t49 = t100 * t142 + t102 * t143 + t171 * t98;
t50 = t101 * t142 + t103 * t143 + t171 * t99;
t65 = t119 * t171 + t120 * t142 + t121 * t143;
t7 = t171 * t49 + t173 * t50 + t186 * t65;
t51 = t100 * t144 + t102 * t145 + t173 * t98;
t52 = t101 * t144 + t103 * t145 + t173 * t99;
t66 = t119 * t173 + t120 * t144 + t121 * t145;
t8 = t171 * t51 + t173 * t52 + t186 * t66;
t235 = t171 * t7 + t173 * t8 + t186 * t26;
t234 = -t139 - t251;
t233 = -t140 - t250;
t232 = -t167 - t240;
t227 = m(3) + m(4) + m(5) + m(6) + m(7);
t15 = t182 * t49 + t184 * t50 + t190 * t65;
t16 = t182 * t51 + t184 * t52 + t190 * t66;
t29 = t182 * t58 + t184 * t59 + t190 * t72;
t225 = t15 * t264 + t16 * t263 + t26 * t259 + t29 * t260 + t8 * t261 + t7 * t262;
t17 = t192 * t49 + t193 * t50 + t198 * t65;
t18 = t192 * t51 + t193 * t52 + t198 * t66;
t31 = t192 * t58 + t193 * t59 + t198 * t72;
t223 = t17 * t264 + t18 * t263 + t26 * t256 + t8 * t257 + t7 * t258 + t31 * t260;
t180 = rSges(4,1) * t191 - rSges(4,2) * t190 + rSges(4,3) * t198;
t179 = Icges(4,1) * t191 - Icges(4,4) * t190 + Icges(4,5) * t198;
t178 = Icges(4,4) * t191 - Icges(4,2) * t190 + Icges(4,6) * t198;
t177 = Icges(4,5) * t191 - Icges(4,6) * t190 + Icges(4,3) * t198;
t161 = rSges(5,1) * t187 - rSges(5,2) * t186 + rSges(5,3) * t190;
t160 = Icges(5,1) * t187 - Icges(5,4) * t186 + Icges(5,5) * t190;
t159 = Icges(5,4) * t187 - Icges(5,2) * t186 + Icges(5,6) * t190;
t158 = Icges(5,5) * t187 - Icges(5,6) * t186 + Icges(5,3) * t190;
t157 = rSges(4,1) * t185 - rSges(4,2) * t184 + rSges(4,3) * t193;
t156 = rSges(4,1) * t183 - rSges(4,2) * t182 + rSges(4,3) * t192;
t155 = Icges(4,1) * t185 - Icges(4,4) * t184 + Icges(4,5) * t193;
t154 = Icges(4,1) * t183 - Icges(4,4) * t182 + Icges(4,5) * t192;
t153 = Icges(4,4) * t185 - Icges(4,2) * t184 + Icges(4,6) * t193;
t152 = Icges(4,4) * t183 - Icges(4,2) * t182 + Icges(4,6) * t192;
t151 = Icges(4,5) * t185 - Icges(4,6) * t184 + Icges(4,3) * t193;
t150 = Icges(4,5) * t183 - Icges(4,6) * t182 + Icges(4,3) * t192;
t141 = t182 * t167;
t136 = t190 * t140;
t135 = t184 * t139;
t133 = Icges(6,1) * t176 + Icges(6,4) * t175 + Icges(6,5) * t186;
t132 = Icges(6,4) * t176 + Icges(6,2) * t175 + Icges(6,6) * t186;
t131 = Icges(6,5) * t176 + Icges(6,6) * t175 + Icges(6,3) * t186;
t130 = rSges(5,1) * t174 - rSges(5,2) * t173 + rSges(5,3) * t184;
t129 = rSges(5,1) * t172 - rSges(5,2) * t171 + rSges(5,3) * t182;
t128 = Icges(5,1) * t174 - Icges(5,4) * t173 + Icges(5,5) * t184;
t127 = Icges(5,1) * t172 - Icges(5,4) * t171 + Icges(5,5) * t182;
t126 = Icges(5,4) * t174 - Icges(5,2) * t173 + Icges(5,6) * t184;
t125 = Icges(5,4) * t172 - Icges(5,2) * t171 + Icges(5,6) * t182;
t124 = Icges(5,5) * t174 - Icges(5,6) * t173 + Icges(5,3) * t184;
t123 = Icges(5,5) * t172 - Icges(5,6) * t171 + Icges(5,3) * t182;
t117 = t157 * t198 - t180 * t193;
t116 = -t156 * t198 + t180 * t192;
t115 = t171 * t122;
t114 = t156 * t193 - t157 * t192;
t111 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t173;
t110 = Icges(6,1) * t147 + Icges(6,4) * t146 + Icges(6,5) * t171;
t109 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t173;
t108 = Icges(6,4) * t147 + Icges(6,2) * t146 + Icges(6,6) * t171;
t107 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t173;
t106 = Icges(6,5) * t147 + Icges(6,6) * t146 + Icges(6,3) * t171;
t95 = t186 * t105;
t94 = t130 * t190 - t161 * t184;
t93 = -t129 * t190 + t161 * t182;
t92 = t173 * t104;
t91 = t158 * t190 - t159 * t186 + t160 * t187;
t90 = t129 * t184 - t130 * t182;
t89 = t130 * t198 + t164 + (-t161 - t181) * t193;
t88 = t161 * t192 + t170 + (-t129 - t165) * t198;
t87 = t158 * t184 - t159 * t173 + t160 * t174;
t86 = t158 * t182 - t159 * t171 + t160 * t172;
t85 = t113 * t186 - t134 * t173;
t84 = -t112 * t186 + t134 * t171;
t83 = -t122 * t173 + t95;
t82 = -t104 * t186 + t115;
t81 = t129 * t193 + t162 + (-t130 - t166) * t192;
t80 = t124 * t190 - t126 * t186 + t128 * t187;
t79 = t123 * t190 - t125 * t186 + t127 * t187;
t78 = t131 * t186 + t132 * t175 + t133 * t176;
t77 = t124 * t184 - t126 * t173 + t128 * t174;
t76 = t123 * t184 - t125 * t173 + t127 * t174;
t75 = t124 * t182 - t126 * t171 + t128 * t172;
t74 = t123 * t182 - t125 * t171 + t127 * t172;
t73 = t112 * t173 - t113 * t171;
t71 = -t105 * t171 + t92;
t70 = t113 * t190 + t239 * t184 + t136;
t69 = t134 * t182 + t242 * t190 + t141;
t68 = t131 * t173 + t132 * t148 + t133 * t149;
t67 = t131 * t171 + t132 * t146 + t133 * t147;
t64 = t113 * t198 + (-t181 + t239) * t193 + t237;
t63 = t134 * t192 + (-t165 + t242) * t198 + t236;
t62 = t112 * t184 + t241 * t182 + t135;
t61 = t107 * t186 + t109 * t175 + t111 * t176;
t60 = t106 * t186 + t108 * t175 + t110 * t176;
t57 = t112 * t193 + (-t166 + t241) * t192 + t238;
t56 = t107 * t173 + t109 * t148 + t111 * t149;
t55 = t106 * t173 + t108 * t148 + t110 * t149;
t54 = t107 * t171 + t109 * t146 + t111 * t147;
t53 = t106 * t171 + t108 * t146 + t110 * t147;
t48 = -t240 * t173 + t186 * t97 + t95;
t47 = t118 * t171 - t251 * t186 + t115;
t46 = t232 * t184 + t250 * t190 + t136;
t45 = t240 * t182 + t234 * t190 + t141;
t44 = t250 * t198 + (-t181 + t232) * t193 + t237;
t43 = t240 * t192 + (-t165 + t234) * t198 + t236;
t42 = -t250 * t171 + t173 * t96 + t92;
t41 = t192 * t79 + t193 * t80 + t198 * t91;
t40 = t233 * t182 + t251 * t184 + t135;
t39 = t182 * t79 + t184 * t80 + t190 * t91;
t38 = t251 * t193 + (-t166 + t233) * t192 + t238;
t37 = t192 * t76 + t193 * t77 + t198 * t87;
t36 = t192 * t74 + t193 * t75 + t198 * t86;
t35 = t182 * t76 + t184 * t77 + t190 * t87;
t34 = t182 * t74 + t184 * t75 + t190 * t86;
t33 = t192 * t60 + t193 * t61 + t198 * t78;
t32 = t182 * t60 + t184 * t61 + t190 * t78;
t28 = t171 * t60 + t173 * t61 + t186 * t78;
t22 = t192 * t55 + t193 * t56 + t198 * t68;
t21 = t192 * t53 + t193 * t54 + t198 * t67;
t20 = t182 * t55 + t184 * t56 + t190 * t68;
t19 = t182 * t53 + t184 * t54 + t190 * t67;
t14 = t171 * t55 + t173 * t56 + t186 * t68;
t13 = t171 * t53 + t173 * t54 + t186 * t67;
t1 = [m(2) + t227; t227 * t212; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t212 ^ 2 + (t209 ^ 2 + t211 ^ 2) * t210 ^ 2); m(4) * t114 + m(5) * t81 + m(6) * t57 + m(7) * t38; m(4) * (t114 * t212 + (t116 * t209 - t117 * t211) * t210) + m(5) * (t212 * t81 + (t209 * t88 - t211 * t89) * t210) + m(6) * (t212 * t57 + (t209 * t63 - t211 * t64) * t210) + m(7) * (t212 * t38 + (t209 * t43 - t211 * t44) * t210); (t31 + t33 + t41 + (t177 * t198 - t178 * t190 + t179 * t191) * t198) * t198 + (t18 + t22 + t37 + (t151 * t193 - t153 * t184 + t155 * t185) * t193 + (t151 * t198 - t153 * t190 + t155 * t191 + t177 * t193 - t178 * t184 + t179 * t185) * t198) * t193 + (t17 + t21 + t36 + (t150 * t192 - t152 * t182 + t154 * t183) * t192 + (t150 * t198 - t152 * t190 + t154 * t191 + t177 * t192 - t178 * t182 + t179 * t183) * t198 + (t150 * t193 + t151 * t192 - t152 * t184 - t153 * t182 + t154 * t185 + t155 * t183) * t193) * t192 + m(7) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t57 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t81 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t114 ^ 2 + t116 ^ 2 + t117 ^ 2); m(5) * t90 + m(6) * t62 + m(7) * t40; m(5) * (t212 * t90 + (t209 * t93 - t211 * t94) * t210) + m(6) * (t212 * t62 + (t209 * t69 - t211 * t70) * t210) + m(7) * (t212 * t40 + (t209 * t45 - t211 * t46) * t210); (t29 / 0.2e1 + t32 / 0.2e1 + t39 / 0.2e1) * t198 + (t16 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1) * t193 + (t15 / 0.2e1 + t19 / 0.2e1 + t34 / 0.2e1) * t192 + (t31 / 0.2e1 + t33 / 0.2e1 + t41 / 0.2e1) * t190 + (t18 / 0.2e1 + t22 / 0.2e1 + t37 / 0.2e1) * t184 + (t17 / 0.2e1 + t21 / 0.2e1 + t36 / 0.2e1) * t182 + m(7) * (t38 * t40 + t43 * t45 + t44 * t46) + m(6) * (t57 * t62 + t63 * t69 + t64 * t70) + m(5) * (t81 * t90 + t88 * t93 + t89 * t94); (t29 + t32 + t39) * t190 + (t16 + t20 + t35) * t184 + (t15 + t19 + t34) * t182 + m(7) * (t40 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t62 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t90 ^ 2 + t93 ^ 2 + t94 ^ 2); m(6) * t73 + m(7) * t42; m(6) * (t212 * t73 + (t209 * t84 - t211 * t85) * t210) + m(7) * (t212 * t42 + (t209 * t47 - t211 * t48) * t210); t13 * t258 + t28 * t256 + t21 * t264 + t22 * t263 + t33 * t260 + t14 * t257 + m(7) * (t38 * t42 + t43 * t47 + t44 * t48) + m(6) * (t57 * t73 + t63 * t84 + t64 * t85) + t223; t19 * t264 + t20 * t263 + t14 * t261 + t13 * t262 + t28 * t259 + t32 * t260 + m(7) * (t40 * t42 + t45 * t47 + t46 * t48) + m(6) * (t62 * t73 + t69 * t84 + t70 * t85) + t225; t171 * t13 + t173 * t14 + t186 * t28 + m(7) * (t42 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t73 ^ 2 + t84 ^ 2 + t85 ^ 2) + t235; m(7) * t71; m(7) * (t212 * t71 + (t209 * t82 - t211 * t83) * t210); m(7) * (t38 * t71 + t43 * t82 + t44 * t83) + t223; m(7) * (t40 * t71 + t45 * t82 + t46 * t83) + t225; m(7) * (t42 * t71 + t47 * t82 + t48 * t83) + t235; m(7) * (t71 ^ 2 + t82 ^ 2 + t83 ^ 2) + t235;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
