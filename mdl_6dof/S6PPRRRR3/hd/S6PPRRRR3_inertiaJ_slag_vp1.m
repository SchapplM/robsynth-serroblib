% Calculate joint inertia matrix for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_inertiaJ_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:24
% EndTime: 2019-03-08 19:07:44
% DurationCPUTime: 8.62s
% Computational Cost: add. (91927->528), mult. (261479->760), div. (0->0), fcn. (349198->18), ass. (0->245)
t269 = m(7) / 0.2e1;
t273 = m(6) / 0.2e1 + t269;
t202 = sin(pkin(14));
t203 = sin(pkin(13));
t206 = cos(pkin(14));
t207 = cos(pkin(13));
t209 = cos(pkin(6));
t250 = t207 * t209;
t197 = -t202 * t203 + t206 * t250;
t204 = sin(pkin(7));
t205 = sin(pkin(6));
t208 = cos(pkin(7));
t252 = t205 * t208;
t192 = -t197 * t204 - t207 * t252;
t272 = -0.2e1 * t192;
t271 = m(5) / 0.2e1;
t198 = t202 * t250 + t203 * t206;
t213 = sin(qJ(3));
t215 = cos(qJ(3));
t254 = t204 * t205;
t221 = t197 * t208 - t207 * t254;
t182 = -t198 * t213 + t221 * t215;
t183 = t198 * t215 + t221 * t213;
t212 = sin(qJ(4));
t256 = sin(pkin(8));
t257 = cos(pkin(8));
t262 = cos(qJ(4));
t152 = t183 * t262 + (t257 * t182 + t256 * t192) * t212;
t172 = -t182 * t256 + t192 * t257;
t211 = sin(qJ(5));
t261 = cos(qJ(5));
t137 = t152 * t211 - t172 * t261;
t255 = t203 * t209;
t200 = -t202 * t255 + t206 * t207;
t199 = -t202 * t207 - t206 * t255;
t220 = t199 * t208 + t203 * t254;
t184 = -t200 * t213 + t220 * t215;
t185 = t200 * t215 + t220 * t213;
t193 = -t199 * t204 + t203 * t252;
t154 = t185 * t262 + (t257 * t184 + t256 * t193) * t212;
t173 = -t184 * t256 + t193 * t257;
t139 = t154 * t211 - t173 * t261;
t251 = t206 * t208;
t253 = t204 * t209;
t190 = t215 * t253 + (-t202 * t213 + t215 * t251) * t205;
t191 = t213 * t253 + (t202 * t215 + t213 * t251) * t205;
t196 = -t206 * t254 + t208 * t209;
t171 = t191 * t262 + (t257 * t190 + t256 * t196) * t212;
t186 = -t190 * t256 + t196 * t257;
t155 = t171 * t211 - t186 * t261;
t138 = t152 * t261 + t172 * t211;
t222 = t262 * t256;
t223 = t257 * t262;
t151 = -t182 * t223 + t183 * t212 - t192 * t222;
t210 = sin(qJ(6));
t214 = cos(qJ(6));
t104 = -t138 * t210 + t151 * t214;
t105 = t138 * t214 + t151 * t210;
t76 = Icges(7,5) * t105 + Icges(7,6) * t104 + Icges(7,3) * t137;
t78 = Icges(7,4) * t105 + Icges(7,2) * t104 + Icges(7,6) * t137;
t80 = Icges(7,1) * t105 + Icges(7,4) * t104 + Icges(7,5) * t137;
t28 = t104 * t78 + t105 * t80 + t137 * t76;
t140 = t154 * t261 + t173 * t211;
t153 = -t184 * t223 + t185 * t212 - t193 * t222;
t106 = -t140 * t210 + t153 * t214;
t107 = t140 * t214 + t153 * t210;
t77 = Icges(7,5) * t107 + Icges(7,6) * t106 + Icges(7,3) * t139;
t79 = Icges(7,4) * t107 + Icges(7,2) * t106 + Icges(7,6) * t139;
t81 = Icges(7,1) * t107 + Icges(7,4) * t106 + Icges(7,5) * t139;
t29 = t104 * t79 + t105 * t81 + t137 * t77;
t156 = t171 * t261 + t186 * t211;
t170 = -t190 * t223 + t191 * t212 - t196 * t222;
t135 = -t156 * t210 + t170 * t214;
t136 = t156 * t214 + t170 * t210;
t90 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t155;
t91 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t155;
t92 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t155;
t41 = t104 * t91 + t105 * t92 + t137 * t90;
t1 = t137 * t28 + t139 * t29 + t155 * t41;
t268 = t1 / 0.2e1;
t30 = t106 * t78 + t107 * t80 + t139 * t76;
t31 = t106 * t79 + t107 * t81 + t139 * t77;
t42 = t106 * t91 + t107 * t92 + t139 * t90;
t2 = t137 * t30 + t139 * t31 + t155 * t42;
t267 = t2 / 0.2e1;
t34 = t135 * t78 + t136 * t80 + t155 * t76;
t35 = t135 * t79 + t136 * t81 + t155 * t77;
t46 = t135 * t91 + t136 * t92 + t155 * t90;
t9 = t137 * t34 + t139 * t35 + t155 * t46;
t266 = t9 / 0.2e1;
t265 = t137 / 0.2e1;
t264 = t139 / 0.2e1;
t263 = t155 / 0.2e1;
t82 = rSges(7,1) * t105 + rSges(7,2) * t104 + rSges(7,3) * t137;
t260 = pkin(5) * t138 + pkin(12) * t137 + t82;
t83 = rSges(7,1) * t107 + rSges(7,2) * t106 + rSges(7,3) * t139;
t259 = pkin(5) * t140 + pkin(12) * t139 + t83;
t93 = rSges(7,1) * t136 + rSges(7,2) * t135 + rSges(7,3) * t155;
t258 = pkin(5) * t156 + pkin(12) * t155 + t93;
t100 = rSges(6,1) * t138 - rSges(6,2) * t137 + rSges(6,3) * t151;
t129 = pkin(4) * t152 + pkin(11) * t151;
t249 = -t100 - t129;
t101 = rSges(6,1) * t140 - rSges(6,2) * t139 + rSges(6,3) * t153;
t130 = pkin(4) * t154 + pkin(11) * t153;
t248 = -t101 - t130;
t111 = rSges(6,1) * t156 - rSges(6,2) * t155 + rSges(6,3) * t170;
t146 = pkin(4) * t171 + pkin(11) * t170;
t247 = -t111 - t146;
t127 = t193 * t129;
t157 = t183 * pkin(3) + pkin(10) * t172;
t149 = t193 * t157;
t245 = t127 + t149;
t158 = t185 * pkin(3) + pkin(10) * t173;
t150 = t196 * t158;
t244 = t196 * t130 + t150;
t174 = t191 * pkin(3) + pkin(10) * t186;
t167 = t192 * t174;
t243 = t192 * t146 + t167;
t242 = t158 * t272 + 0.2e1 * t149;
t94 = Icges(6,5) * t138 - Icges(6,6) * t137 + Icges(6,3) * t151;
t96 = Icges(6,4) * t138 - Icges(6,2) * t137 + Icges(6,6) * t151;
t98 = Icges(6,1) * t138 - Icges(6,4) * t137 + Icges(6,5) * t151;
t47 = -t137 * t96 + t138 * t98 + t151 * t94;
t95 = Icges(6,5) * t140 - Icges(6,6) * t139 + Icges(6,3) * t153;
t97 = Icges(6,4) * t140 - Icges(6,2) * t139 + Icges(6,6) * t153;
t99 = Icges(6,1) * t140 - Icges(6,4) * t139 + Icges(6,5) * t153;
t48 = -t137 * t97 + t138 * t99 + t151 * t95;
t108 = Icges(6,5) * t156 - Icges(6,6) * t155 + Icges(6,3) * t170;
t109 = Icges(6,4) * t156 - Icges(6,2) * t155 + Icges(6,6) * t170;
t110 = Icges(6,1) * t156 - Icges(6,4) * t155 + Icges(6,5) * t170;
t59 = t108 * t151 - t109 * t137 + t110 * t138;
t13 = t151 * t47 + t153 * t48 + t170 * t59;
t3 = t151 * t28 + t153 * t29 + t170 * t41;
t240 = t3 / 0.2e1 + t13 / 0.2e1;
t49 = -t139 * t96 + t140 * t98 + t153 * t94;
t50 = -t139 * t97 + t140 * t99 + t153 * t95;
t60 = t108 * t153 - t109 * t139 + t110 * t140;
t14 = t151 * t49 + t153 * t50 + t170 * t60;
t4 = t151 * t30 + t153 * t31 + t170 * t42;
t239 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t172 * t47 + t173 * t48 + t186 * t59;
t5 = t172 * t28 + t173 * t29 + t186 * t41;
t238 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t172 * t49 + t173 * t50 + t186 * t60;
t6 = t172 * t30 + t173 * t31 + t186 * t42;
t237 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t192 * t47 + t193 * t48 + t196 * t59;
t7 = t192 * t28 + t193 * t29 + t196 * t41;
t236 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t192 * t49 + t193 * t50 + t196 * t60;
t8 = t192 * t30 + t193 * t31 + t196 * t42;
t235 = t8 / 0.2e1 + t18 / 0.2e1;
t10 = t151 * t34 + t153 * t35 + t170 * t46;
t52 = -t155 * t96 + t156 * t98 + t170 * t94;
t53 = -t155 * t97 + t156 * t99 + t170 * t95;
t64 = t108 * t170 - t109 * t155 + t110 * t156;
t19 = t151 * t52 + t153 * t53 + t170 * t64;
t234 = t10 / 0.2e1 + t19 / 0.2e1;
t11 = t172 * t34 + t173 * t35 + t186 * t46;
t20 = t172 * t52 + t173 * t53 + t186 * t64;
t233 = t11 / 0.2e1 + t20 / 0.2e1;
t12 = t192 * t34 + t193 * t35 + t196 * t46;
t21 = t192 * t52 + t193 * t53 + t196 * t64;
t232 = t12 / 0.2e1 + t21 / 0.2e1;
t231 = -t129 - t260;
t230 = -t130 - t259;
t229 = -t146 - t258;
t225 = m(3) + m(4) + m(5) + m(6) + m(7);
t45 = -t137 * t83 + t139 * t82;
t219 = m(6) * t100 + m(7) * t260;
t218 = -m(6) * t101 - m(7) * t259;
t118 = rSges(5,1) * t152 - rSges(5,2) * t151 + rSges(5,3) * t172;
t217 = m(5) * t118 + t219;
t119 = rSges(5,1) * t154 - rSges(5,2) * t153 + rSges(5,3) * t173;
t216 = -m(5) * t119 + t218;
t178 = rSges(4,1) * t191 + rSges(4,2) * t190 + rSges(4,3) * t196;
t177 = Icges(4,1) * t191 + Icges(4,4) * t190 + Icges(4,5) * t196;
t176 = Icges(4,4) * t191 + Icges(4,2) * t190 + Icges(4,6) * t196;
t175 = Icges(4,5) * t191 + Icges(4,6) * t190 + Icges(4,3) * t196;
t166 = rSges(4,1) * t185 + rSges(4,2) * t184 + rSges(4,3) * t193;
t165 = rSges(4,1) * t183 + rSges(4,2) * t182 + rSges(4,3) * t192;
t164 = Icges(4,1) * t185 + Icges(4,4) * t184 + Icges(4,5) * t193;
t163 = Icges(4,1) * t183 + Icges(4,4) * t182 + Icges(4,5) * t192;
t162 = Icges(4,4) * t185 + Icges(4,2) * t184 + Icges(4,6) * t193;
t161 = Icges(4,4) * t183 + Icges(4,2) * t182 + Icges(4,6) * t192;
t160 = Icges(4,5) * t185 + Icges(4,6) * t184 + Icges(4,3) * t193;
t159 = Icges(4,5) * t183 + Icges(4,6) * t182 + Icges(4,3) * t192;
t144 = rSges(5,1) * t171 - rSges(5,2) * t170 + rSges(5,3) * t186;
t143 = Icges(5,1) * t171 - Icges(5,4) * t170 + Icges(5,5) * t186;
t142 = Icges(5,4) * t171 - Icges(5,2) * t170 + Icges(5,6) * t186;
t141 = Icges(5,5) * t171 - Icges(5,6) * t170 + Icges(5,3) * t186;
t134 = t166 * t196 - t178 * t193;
t133 = -t165 * t196 + t178 * t192;
t132 = t172 * t146;
t124 = t165 * t193 - t166 * t192;
t123 = t186 * t130;
t122 = t173 * t129;
t117 = Icges(5,1) * t154 - Icges(5,4) * t153 + Icges(5,5) * t173;
t116 = Icges(5,1) * t152 - Icges(5,4) * t151 + Icges(5,5) * t172;
t115 = Icges(5,4) * t154 - Icges(5,2) * t153 + Icges(5,6) * t173;
t114 = Icges(5,4) * t152 - Icges(5,2) * t151 + Icges(5,6) * t172;
t113 = Icges(5,5) * t154 - Icges(5,6) * t153 + Icges(5,3) * t173;
t112 = Icges(5,5) * t152 - Icges(5,6) * t151 + Icges(5,3) * t172;
t89 = t119 * t186 - t144 * t173;
t88 = -t118 * t186 + t144 * t172;
t87 = t141 * t186 - t142 * t170 + t143 * t171;
t86 = t196 * t119 + t150 + (-t144 - t174) * t193;
t85 = t192 * t144 + t167 + (-t118 - t157) * t196;
t84 = t118 * t173 - t119 * t172;
t75 = t141 * t173 - t142 * t153 + t143 * t154;
t74 = t141 * t172 - t142 * t151 + t143 * t152;
t73 = t193 * t118 + t149 + (-t119 - t158) * t192;
t72 = t101 * t170 - t111 * t153;
t71 = -t100 * t170 + t111 * t151;
t70 = t113 * t186 - t115 * t170 + t117 * t171;
t69 = t112 * t186 - t114 * t170 + t116 * t171;
t68 = t113 * t173 - t115 * t153 + t117 * t154;
t67 = t112 * t173 - t114 * t153 + t116 * t154;
t66 = t113 * t172 - t115 * t151 + t117 * t152;
t65 = t112 * t172 - t114 * t151 + t116 * t152;
t63 = t100 * t153 - t101 * t151;
t62 = t186 * t101 + t247 * t173 + t123;
t61 = t172 * t111 + t249 * t186 + t132;
t58 = t196 * t101 + (-t174 + t247) * t193 + t244;
t57 = t192 * t111 + (-t157 + t249) * t196 + t243;
t56 = -t139 * t93 + t155 * t83;
t55 = t137 * t93 - t155 * t82;
t54 = t173 * t100 + t248 * t172 + t122;
t51 = t193 * t100 + (-t158 + t248) * t192 + t245;
t44 = -t258 * t153 + t259 * t170;
t43 = t258 * t151 - t260 * t170;
t40 = t259 * t196 + (-t174 + t229) * t193 + t244;
t39 = t258 * t192 + (-t157 + t231) * t196 + t243;
t38 = t229 * t173 + t259 * t186 + t123;
t37 = t258 * t172 + t231 * t186 + t132;
t36 = -t259 * t151 + t260 * t153;
t33 = t192 * t69 + t193 * t70 + t196 * t87;
t32 = t172 * t69 + t173 * t70 + t186 * t87;
t27 = t260 * t193 + (-t158 + t230) * t192 + t245;
t26 = t230 * t172 + t260 * t173 + t122;
t25 = t192 * t67 + t193 * t68 + t196 * t75;
t24 = t192 * t65 + t193 * t66 + t196 * t74;
t23 = t172 * t67 + t173 * t68 + t186 * t75;
t22 = t172 * t65 + t173 * t66 + t186 * t74;
t102 = [m(2) + t225; t225 * t209; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t271 + t273) * (t209 ^ 2 + (t203 ^ 2 + t207 ^ 2) * t205 ^ 2); t242 * t271 + (m(4) * t165 + t217) * t193 + (-m(4) * t166 + t216) * t192 + t273 * (t130 * t272 + 0.2e1 * t127 + t242); (m(4) * t124 + m(5) * t73 + m(6) * t51 + m(7) * t27) * t209 + ((-m(4) * t134 - m(5) * t86 - m(6) * t58 - m(7) * t40) * t207 + (m(4) * t133 + m(5) * t85 + m(6) * t57 + m(7) * t39) * t203) * t205; (t27 ^ 2 + t39 ^ 2 + t40 ^ 2) * m(7) + (t51 ^ 2 + t57 ^ 2 + t58 ^ 2) * m(6) + (t73 ^ 2 + t85 ^ 2 + t86 ^ 2) * m(5) + m(4) * (t124 ^ 2 + t133 ^ 2 + t134 ^ 2) + (t12 + t21 + t33 + (t175 * t196 + t176 * t190 + t177 * t191) * t196) * t196 + (t8 + t18 + t25 + (t160 * t193 + t162 * t184 + t164 * t185) * t193 + (t160 * t196 + t162 * t190 + t164 * t191 + t175 * t193 + t176 * t184 + t177 * t185) * t196) * t193 + (t7 + t17 + t24 + (t159 * t192 + t161 * t182 + t163 * t183) * t192 + (t159 * t196 + t161 * t190 + t163 * t191 + t175 * t192 + t176 * t182 + t177 * t183) * t196 + (t159 * t193 + t160 * t192 + t161 * t184 + t162 * t182 + t163 * t185 + t164 * t183) * t193) * t192; t216 * t172 + t217 * t173 + 0.2e1 * t273 * (-t172 * t130 + t122); (m(5) * t84 + m(6) * t54 + m(7) * t26) * t209 + ((-m(5) * t89 - m(6) * t62 - m(7) * t38) * t207 + (m(5) * t88 + m(6) * t61 + m(7) * t37) * t203) * t205; (t26 * t27 + t37 * t39 + t38 * t40) * m(7) + (t51 * t54 + t57 * t61 + t58 * t62) * m(6) + (t73 * t84 + t85 * t88 + t86 * t89) * m(5) + (t32 / 0.2e1 + t233) * t196 + (t23 / 0.2e1 + t237) * t193 + (t22 / 0.2e1 + t238) * t192 + (t33 / 0.2e1 + t232) * t186 + (t25 / 0.2e1 + t235) * t173 + (t24 / 0.2e1 + t236) * t172; (t26 ^ 2 + t37 ^ 2 + t38 ^ 2) * m(7) + (t54 ^ 2 + t61 ^ 2 + t62 ^ 2) * m(6) + (t84 ^ 2 + t88 ^ 2 + t89 ^ 2) * m(5) + (t11 + t20 + t32) * t186 + (t6 + t16 + t23) * t173 + (t5 + t15 + t22) * t172; t151 * t218 + t153 * t219; (m(6) * t63 + m(7) * t36) * t209 + ((-m(6) * t72 - m(7) * t44) * t207 + (m(6) * t71 + m(7) * t43) * t203) * t205; (t27 * t36 + t39 * t43 + t40 * t44) * m(7) + (t51 * t63 + t57 * t71 + t58 * t72) * m(6) + t234 * t196 + t239 * t193 + t240 * t192 + t232 * t170 + t235 * t153 + t236 * t151; (t26 * t36 + t37 * t43 + t38 * t44) * m(7) + (t54 * t63 + t61 * t71 + t62 * t72) * m(6) + t234 * t186 + t239 * t173 + t240 * t172 + t233 * t170 + t237 * t153 + t238 * t151; (t36 ^ 2 + t43 ^ 2 + t44 ^ 2) * m(7) + (t63 ^ 2 + t71 ^ 2 + t72 ^ 2) * m(6) + (t10 + t19) * t170 + (t4 + t14) * t153 + (t3 + t13) * t151; t45 * m(7); 0.2e1 * (t45 * t209 + (t203 * t55 - t207 * t56) * t205) * t269; t192 * t268 + t12 * t263 + t8 * t264 + (t27 * t45 + t39 * t55 + t40 * t56) * m(7) + t7 * t265 + t193 * t267 + t196 * t266; (t26 * t45 + t37 * t55 + t38 * t56) * m(7) + t6 * t264 + t186 * t266 + t173 * t267 + t11 * t263 + t5 * t265 + t172 * t268; t153 * t267 + t3 * t265 + t151 * t268 + (t36 * t45 + t43 * t55 + t44 * t56) * m(7) + t4 * t264 + t10 * t263 + t170 * t266; t137 * t1 + t155 * t9 + t139 * t2 + (t45 ^ 2 + t55 ^ 2 + t56 ^ 2) * m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t102(1) t102(2) t102(4) t102(7) t102(11) t102(16); t102(2) t102(3) t102(5) t102(8) t102(12) t102(17); t102(4) t102(5) t102(6) t102(9) t102(13) t102(18); t102(7) t102(8) t102(9) t102(10) t102(14) t102(19); t102(11) t102(12) t102(13) t102(14) t102(15) t102(20); t102(16) t102(17) t102(18) t102(19) t102(20) t102(21);];
Mq  = res;
