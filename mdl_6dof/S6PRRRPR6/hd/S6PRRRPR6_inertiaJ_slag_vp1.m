% Calculate joint inertia matrix for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:14
% EndTime: 2019-03-08 23:31:27
% DurationCPUTime: 6.02s
% Computational Cost: add. (27675->608), mult. (71734->862), div. (0->0), fcn. (93927->12), ass. (0->269)
t230 = sin(pkin(11));
t232 = cos(pkin(11));
t239 = cos(qJ(2));
t233 = cos(pkin(6));
t237 = sin(qJ(2));
t279 = t233 * t237;
t221 = t230 * t239 + t232 * t279;
t236 = sin(qJ(3));
t231 = sin(pkin(6));
t285 = cos(qJ(3));
t258 = t231 * t285;
t208 = t221 * t236 + t232 * t258;
t223 = -t230 * t279 + t232 * t239;
t210 = t223 * t236 - t230 * t258;
t281 = t231 * t236;
t224 = -t233 * t285 + t237 * t281;
t209 = t221 * t285 - t232 * t281;
t278 = t233 * t239;
t220 = t230 * t237 - t232 * t278;
t235 = sin(qJ(4));
t284 = cos(qJ(4));
t184 = t209 * t235 - t220 * t284;
t185 = t209 * t284 + t220 * t235;
t234 = sin(qJ(6));
t238 = cos(qJ(6));
t145 = t184 * t238 - t185 * t234;
t146 = t184 * t234 + t185 * t238;
t100 = Icges(7,5) * t146 + Icges(7,6) * t145 - Icges(7,3) * t208;
t102 = Icges(7,4) * t146 + Icges(7,2) * t145 - Icges(7,6) * t208;
t104 = Icges(7,1) * t146 + Icges(7,4) * t145 - Icges(7,5) * t208;
t34 = -t100 * t208 + t102 * t145 + t104 * t146;
t211 = t223 * t285 + t230 * t281;
t222 = t230 * t278 + t232 * t237;
t186 = t211 * t235 - t222 * t284;
t187 = t211 * t284 + t222 * t235;
t147 = t186 * t238 - t187 * t234;
t148 = t186 * t234 + t187 * t238;
t101 = Icges(7,5) * t148 + Icges(7,6) * t147 - Icges(7,3) * t210;
t103 = Icges(7,4) * t148 + Icges(7,2) * t147 - Icges(7,6) * t210;
t105 = Icges(7,1) * t148 + Icges(7,4) * t147 - Icges(7,5) * t210;
t35 = -t101 * t208 + t103 * t145 + t105 * t146;
t225 = t233 * t236 + t237 * t258;
t280 = t231 * t239;
t212 = t225 * t235 + t280 * t284;
t213 = t225 * t284 - t235 * t280;
t179 = t212 * t238 - t213 * t234;
t180 = t212 * t234 + t213 * t238;
t115 = Icges(7,5) * t180 + Icges(7,6) * t179 - Icges(7,3) * t224;
t116 = Icges(7,4) * t180 + Icges(7,2) * t179 - Icges(7,6) * t224;
t117 = Icges(7,1) * t180 + Icges(7,4) * t179 - Icges(7,5) * t224;
t54 = -t115 * t208 + t116 * t145 + t117 * t146;
t1 = t208 * t34 + t210 * t35 + t224 * t54;
t291 = -t1 / 0.2e1;
t294 = t1 / 0.2e1;
t36 = -t100 * t210 + t102 * t147 + t104 * t148;
t37 = -t101 * t210 + t103 * t147 + t105 * t148;
t55 = -t115 * t210 + t116 * t147 + t117 * t148;
t3 = t208 * t36 + t210 * t37 + t224 * t55;
t290 = -t3 / 0.2e1;
t38 = -t100 * t224 + t102 * t179 + t104 * t180;
t39 = -t101 * t224 + t103 * t179 + t105 * t180;
t57 = -t115 * t224 + t116 * t179 + t117 * t180;
t9 = t208 * t38 + t210 * t39 + t224 * t57;
t289 = -t9 / 0.2e1;
t293 = t9 / 0.2e1;
t292 = m(6) + m(7);
t288 = -t208 / 0.2e1;
t287 = -t210 / 0.2e1;
t286 = -t224 / 0.2e1;
t283 = t230 * t231;
t282 = t231 * t232;
t106 = rSges(7,1) * t146 + rSges(7,2) * t145 - rSges(7,3) * t208;
t277 = pkin(5) * t185 - pkin(10) * t208 + t106;
t107 = rSges(7,1) * t148 + rSges(7,2) * t147 - rSges(7,3) * t210;
t276 = pkin(5) * t187 - pkin(10) * t210 + t107;
t118 = rSges(7,1) * t180 + rSges(7,2) * t179 - rSges(7,3) * t224;
t275 = pkin(5) * t213 - pkin(10) * t224 + t118;
t135 = rSges(6,1) * t187 + rSges(6,2) * t210 + rSges(6,3) * t186;
t150 = pkin(4) * t187 + qJ(5) * t186;
t274 = -t135 - t150;
t136 = rSges(5,1) * t187 - rSges(5,2) * t186 + rSges(5,3) * t210;
t182 = pkin(3) * t211 + pkin(9) * t210;
t273 = -t136 - t182;
t149 = pkin(4) * t185 + qJ(5) * t184;
t181 = pkin(3) * t209 + pkin(9) * t208;
t169 = t222 * t181;
t272 = t222 * t149 + t169;
t170 = rSges(6,1) * t213 + rSges(6,2) * t224 + rSges(6,3) * t212;
t183 = pkin(4) * t213 + qJ(5) * t212;
t271 = -t170 - t183;
t171 = rSges(5,1) * t213 - rSges(5,2) * t212 + rSges(5,3) * t224;
t207 = pkin(3) * t225 + pkin(9) * t224;
t270 = -t171 - t207;
t269 = t181 * t280 + t220 * t207;
t206 = pkin(2) * t223 + pkin(8) * t222;
t204 = t233 * t206;
t268 = t233 * t182 + t204;
t205 = pkin(2) * t221 + pkin(8) * t220;
t267 = -t181 - t205;
t266 = t205 * t283 + t206 * t282;
t264 = -t150 - t276;
t263 = -t183 - t275;
t262 = -t182 + t274;
t261 = t233 * t150 + t268;
t260 = -t149 + t267;
t259 = -t207 + t271;
t201 = t225 * rSges(4,1) - t224 * rSges(4,2) - rSges(4,3) * t280;
t226 = (pkin(2) * t237 - pkin(8) * t239) * t231;
t257 = (-t201 - t226) * t231;
t256 = -t182 + t264;
t255 = -t207 + t263;
t254 = t149 * t280 + t220 * t183 + t269;
t253 = t181 * t283 + t182 * t282 + t266;
t121 = Icges(6,5) * t185 + Icges(6,6) * t208 + Icges(6,3) * t184;
t125 = Icges(6,4) * t185 + Icges(6,2) * t208 + Icges(6,6) * t184;
t129 = Icges(6,1) * t185 + Icges(6,4) * t208 + Icges(6,5) * t184;
t61 = t121 * t184 + t125 * t208 + t129 * t185;
t122 = Icges(6,5) * t187 + Icges(6,6) * t210 + Icges(6,3) * t186;
t126 = Icges(6,4) * t187 + Icges(6,2) * t210 + Icges(6,6) * t186;
t130 = Icges(6,1) * t187 + Icges(6,4) * t210 + Icges(6,5) * t186;
t62 = t122 * t184 + t126 * t208 + t130 * t185;
t163 = Icges(6,5) * t213 + Icges(6,6) * t224 + Icges(6,3) * t212;
t165 = Icges(6,4) * t213 + Icges(6,2) * t224 + Icges(6,6) * t212;
t167 = Icges(6,1) * t213 + Icges(6,4) * t224 + Icges(6,5) * t212;
t83 = t163 * t184 + t165 * t208 + t167 * t185;
t13 = t208 * t61 + t210 * t62 + t224 * t83;
t123 = Icges(5,5) * t185 - Icges(5,6) * t184 + Icges(5,3) * t208;
t127 = Icges(5,4) * t185 - Icges(5,2) * t184 + Icges(5,6) * t208;
t131 = Icges(5,1) * t185 - Icges(5,4) * t184 + Icges(5,5) * t208;
t63 = t123 * t208 - t127 * t184 + t131 * t185;
t124 = Icges(5,5) * t187 - Icges(5,6) * t186 + Icges(5,3) * t210;
t128 = Icges(5,4) * t187 - Icges(5,2) * t186 + Icges(5,6) * t210;
t132 = Icges(5,1) * t187 - Icges(5,4) * t186 + Icges(5,5) * t210;
t64 = t124 * t208 - t128 * t184 + t132 * t185;
t164 = Icges(5,5) * t213 - Icges(5,6) * t212 + Icges(5,3) * t224;
t166 = Icges(5,4) * t213 - Icges(5,2) * t212 + Icges(5,6) * t224;
t168 = Icges(5,1) * t213 - Icges(5,4) * t212 + Icges(5,5) * t224;
t84 = t164 * t208 - t166 * t184 + t168 * t185;
t14 = t208 * t63 + t210 * t64 + t224 * t84;
t252 = t294 + t13 / 0.2e1 + t14 / 0.2e1;
t65 = t121 * t186 + t125 * t210 + t129 * t187;
t66 = t122 * t186 + t126 * t210 + t130 * t187;
t85 = t163 * t186 + t165 * t210 + t167 * t187;
t15 = t208 * t65 + t210 * t66 + t224 * t85;
t67 = t123 * t210 - t127 * t186 + t131 * t187;
t68 = t124 * t210 - t128 * t186 + t132 * t187;
t86 = t164 * t210 - t166 * t186 + t168 * t187;
t16 = t208 * t67 + t210 * t68 + t224 * t86;
t251 = t3 / 0.2e1 + t15 / 0.2e1 + t16 / 0.2e1;
t17 = t61 * t220 + t62 * t222 - t280 * t83;
t18 = t63 * t220 + t64 * t222 - t280 * t84;
t5 = t34 * t220 + t35 * t222 - t280 * t54;
t250 = t5 / 0.2e1 + t17 / 0.2e1 + t18 / 0.2e1;
t19 = t65 * t220 + t66 * t222 - t280 * t85;
t20 = t67 * t220 + t68 * t222 - t280 * t86;
t6 = t36 * t220 + t37 * t222 - t280 * t55;
t249 = t6 / 0.2e1 + t20 / 0.2e1 + t19 / 0.2e1;
t21 = t233 * t83 + (t230 * t62 - t232 * t61) * t231;
t22 = t233 * t84 + (t230 * t64 - t232 * t63) * t231;
t7 = t233 * t54 + (t230 * t35 - t232 * t34) * t231;
t248 = t7 / 0.2e1 + t22 / 0.2e1 + t21 / 0.2e1;
t23 = t233 * t85 + (t230 * t66 - t232 * t65) * t231;
t24 = t233 * t86 + (t230 * t68 - t232 * t67) * t231;
t8 = t233 * t55 + (t230 * t37 - t232 * t36) * t231;
t247 = t8 / 0.2e1 + t24 / 0.2e1 + t23 / 0.2e1;
t72 = t121 * t212 + t125 * t224 + t129 * t213;
t73 = t122 * t212 + t126 * t224 + t130 * t213;
t94 = t163 * t212 + t165 * t224 + t167 * t213;
t25 = t208 * t72 + t210 * t73 + t224 * t94;
t74 = t123 * t224 - t127 * t212 + t131 * t213;
t75 = t124 * t224 - t128 * t212 + t132 * t213;
t95 = t164 * t224 - t166 * t212 + t168 * t213;
t26 = t208 * t74 + t210 * t75 + t224 * t95;
t246 = t293 + t25 / 0.2e1 + t26 / 0.2e1;
t11 = t38 * t220 + t39 * t222 - t280 * t57;
t27 = t72 * t220 + t73 * t222 - t280 * t94;
t28 = t74 * t220 + t75 * t222 - t280 * t95;
t245 = t11 / 0.2e1 + t28 / 0.2e1 + t27 / 0.2e1;
t12 = t233 * t57 + (t230 * t39 - t232 * t38) * t231;
t29 = t233 * t94 + (t230 * t73 - t232 * t72) * t231;
t30 = t233 * t95 + (t230 * t75 - t232 * t74) * t231;
t244 = t12 / 0.2e1 + t29 / 0.2e1 + t30 / 0.2e1;
t243 = (-t226 + t270) * t231;
t242 = (-t226 + t259) * t231;
t241 = t149 * t283 + t150 * t282 + t253;
t240 = (-t226 + t255) * t231;
t217 = t233 * rSges(3,3) + (rSges(3,1) * t237 + rSges(3,2) * t239) * t231;
t216 = Icges(3,5) * t233 + (Icges(3,1) * t237 + Icges(3,4) * t239) * t231;
t215 = Icges(3,6) * t233 + (Icges(3,4) * t237 + Icges(3,2) * t239) * t231;
t214 = Icges(3,3) * t233 + (Icges(3,5) * t237 + Icges(3,6) * t239) * t231;
t200 = Icges(4,1) * t225 - Icges(4,4) * t224 - Icges(4,5) * t280;
t199 = Icges(4,4) * t225 - Icges(4,2) * t224 - Icges(4,6) * t280;
t198 = Icges(4,5) * t225 - Icges(4,6) * t224 - Icges(4,3) * t280;
t197 = rSges(3,1) * t223 - rSges(3,2) * t222 + rSges(3,3) * t283;
t196 = rSges(3,1) * t221 - rSges(3,2) * t220 - rSges(3,3) * t282;
t195 = Icges(3,1) * t223 - Icges(3,4) * t222 + Icges(3,5) * t283;
t194 = Icges(3,1) * t221 - Icges(3,4) * t220 - Icges(3,5) * t282;
t193 = Icges(3,4) * t223 - Icges(3,2) * t222 + Icges(3,6) * t283;
t192 = Icges(3,4) * t221 - Icges(3,2) * t220 - Icges(3,6) * t282;
t191 = Icges(3,5) * t223 - Icges(3,6) * t222 + Icges(3,3) * t283;
t190 = Icges(3,5) * t221 - Icges(3,6) * t220 - Icges(3,3) * t282;
t174 = -t196 * t233 - t217 * t282;
t173 = t197 * t233 - t217 * t283;
t162 = rSges(4,1) * t211 - rSges(4,2) * t210 + rSges(4,3) * t222;
t161 = rSges(4,1) * t209 - rSges(4,2) * t208 + rSges(4,3) * t220;
t160 = Icges(4,1) * t211 - Icges(4,4) * t210 + Icges(4,5) * t222;
t159 = Icges(4,1) * t209 - Icges(4,4) * t208 + Icges(4,5) * t220;
t158 = Icges(4,4) * t211 - Icges(4,2) * t210 + Icges(4,6) * t222;
t157 = Icges(4,4) * t209 - Icges(4,2) * t208 + Icges(4,6) * t220;
t156 = Icges(4,5) * t211 - Icges(4,6) * t210 + Icges(4,3) * t222;
t155 = Icges(4,5) * t209 - Icges(4,6) * t208 + Icges(4,3) * t220;
t152 = t208 * t183;
t151 = (t196 * t230 + t197 * t232) * t231;
t140 = t224 * t150;
t137 = t210 * t149;
t134 = rSges(5,1) * t185 - rSges(5,2) * t184 + rSges(5,3) * t208;
t133 = rSges(6,1) * t185 + rSges(6,2) * t208 + rSges(6,3) * t184;
t120 = -t162 * t280 - t222 * t201;
t119 = t161 * t280 + t220 * t201;
t114 = -t198 * t280 - t224 * t199 + t225 * t200;
t113 = t161 * t222 - t162 * t220;
t112 = (-t161 - t205) * t233 + t232 * t257;
t111 = t162 * t233 + t230 * t257 + t204;
t110 = t198 * t222 - t199 * t210 + t200 * t211;
t109 = t198 * t220 - t199 * t208 + t200 * t209;
t108 = (t161 * t230 + t162 * t232) * t231 + t266;
t99 = t136 * t224 - t171 * t210;
t98 = -t134 * t224 + t171 * t208;
t97 = -t156 * t280 - t224 * t158 + t225 * t160;
t96 = -t155 * t280 - t224 * t157 + t225 * t159;
t93 = t156 * t222 - t158 * t210 + t160 * t211;
t92 = t155 * t222 - t157 * t210 + t159 * t211;
t91 = t156 * t220 - t158 * t208 + t160 * t209;
t90 = t155 * t220 - t157 * t208 + t159 * t209;
t89 = t134 * t210 - t136 * t208;
t88 = t222 * t270 + t273 * t280;
t87 = t134 * t280 + t220 * t171 + t269;
t82 = (-t134 + t267) * t233 + t232 * t243;
t81 = t136 * t233 + t230 * t243 + t268;
t80 = -t107 * t224 + t118 * t210;
t79 = t106 * t224 - t118 * t208;
t78 = t134 * t222 + t220 * t273 + t169;
t77 = t135 * t224 + t210 * t271 + t140;
t76 = t170 * t208 + t152 + (-t133 - t149) * t224;
t71 = (t134 * t230 + t136 * t232) * t231 + t253;
t70 = t222 * t259 + t262 * t280;
t69 = t133 * t280 + t220 * t170 + t254;
t60 = (-t133 + t260) * t233 + t232 * t242;
t59 = t135 * t233 + t230 * t242 + t261;
t58 = -t106 * t210 + t107 * t208;
t56 = t133 * t210 + t208 * t274 + t137;
t53 = t133 * t222 + t220 * t262 + t272;
t52 = (t133 * t230 + t135 * t232) * t231 + t241;
t51 = t114 * t233 + (t230 * t97 - t232 * t96) * t231;
t50 = -t114 * t280 + t96 * t220 + t97 * t222;
t49 = t210 * t263 + t224 * t276 + t140;
t48 = t152 + t275 * t208 + (-t149 - t277) * t224;
t47 = t110 * t233 + (t230 * t93 - t232 * t92) * t231;
t46 = t109 * t233 + (t230 * t91 - t232 * t90) * t231;
t45 = -t110 * t280 + t92 * t220 + t93 * t222;
t44 = -t109 * t280 + t90 * t220 + t91 * t222;
t43 = t222 * t255 + t256 * t280;
t42 = t220 * t275 + t277 * t280 + t254;
t41 = (t260 - t277) * t233 + t232 * t240;
t40 = t230 * t240 + t233 * t276 + t261;
t33 = t208 * t264 + t210 * t277 + t137;
t32 = t220 * t256 + t222 * t277 + t272;
t31 = (t230 * t277 + t232 * t276) * t231 + t241;
t2 = [m(3) + m(4) + m(5) + m(2) + t292; m(3) * t151 + m(4) * t108 + m(5) * t71 + m(6) * t52 + m(7) * t31; m(7) * (t31 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t52 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t71 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t108 ^ 2 + t111 ^ 2 + t112 ^ 2) + m(3) * (t151 ^ 2 + t173 ^ 2 + t174 ^ 2) + (t8 + t23 + t24 + t47 + (t191 * t283 - t193 * t222 + t195 * t223) * t283) * t283 + (-t7 - t21 - t22 - t46 + (-t190 * t282 - t192 * t220 + t194 * t221) * t282 + (-t190 * t283 + t191 * t282 + t192 * t222 + t193 * t220 - t194 * t223 - t195 * t221) * t283) * t282 + ((t214 * t283 - t222 * t215 + t223 * t216) * t283 - (-t214 * t282 - t220 * t215 + t221 * t216) * t282 + t12 + t29 + t30 + t51 + ((t193 * t239 + t195 * t237) * t230 - (t192 * t239 + t194 * t237) * t232) * t231 ^ 2 + ((-t190 * t232 + t191 * t230 + t215 * t239 + t216 * t237) * t231 + t233 * t214) * t233) * t233; m(4) * t113 + m(5) * t78 + m(6) * t53 + m(7) * t32; (t50 / 0.2e1 + t245) * t233 + (t47 / 0.2e1 + t247) * t222 + (t46 / 0.2e1 + t248) * t220 + m(6) * (t52 * t53 + t59 * t70 + t60 * t69) + m(4) * (t108 * t113 + t111 * t120 + t112 * t119) + m(5) * (t71 * t78 + t81 * t88 + t82 * t87) + m(7) * (t31 * t32 + t40 * t43 + t41 * t42) + ((-t51 / 0.2e1 - t244) * t239 + (-t44 / 0.2e1 - t250) * t232 + (t45 / 0.2e1 + t249) * t230) * t231; (-t11 - t27 - t28 - t50) * t280 + (t6 + t19 + t20 + t45) * t222 + (t5 + t17 + t18 + t44) * t220 + m(7) * (t32 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t53 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t78 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t113 ^ 2 + t119 ^ 2 + t120 ^ 2); m(5) * t89 + m(6) * t56 + m(7) * t33; t246 * t233 + t244 * t224 + t247 * t210 + t248 * t208 + m(7) * (t31 * t33 + t40 * t49 + t41 * t48) + m(6) * (t52 * t56 + t59 * t77 + t60 * t76) + m(5) * (t71 * t89 + t81 * t99 + t82 * t98) + (t230 * t251 - t232 * t252) * t231; -t246 * t280 + t245 * t224 + t251 * t222 + t252 * t220 + t249 * t210 + t250 * t208 + m(7) * (t32 * t33 + t42 * t48 + t43 * t49) + m(6) * (t53 * t56 + t69 * t76 + t70 * t77) + m(5) * (t78 * t89 + t87 * t98 + t88 * t99); (t9 + t26 + t25) * t224 + (t3 + t15 + t16) * t210 + (t1 + t13 + t14) * t208 + m(7) * (t33 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t56 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t89 ^ 2 + t98 ^ 2 + t99 ^ 2); t212 * t292; m(7) * (t184 * t40 + t186 * t41 + t212 * t31) + m(6) * (t184 * t59 + t186 * t60 + t212 * t52); m(7) * (t184 * t43 + t186 * t42 + t212 * t32) + m(6) * (t184 * t70 + t186 * t69 + t212 * t53); m(7) * (t184 * t49 + t186 * t48 + t212 * t33) + m(6) * (t184 * t77 + t186 * t76 + t212 * t56); (t184 ^ 2 + t186 ^ 2 + t212 ^ 2) * t292; m(7) * t58; m(7) * (t31 * t58 + t40 * t80 + t41 * t79) + t7 * t288 + t8 * t287 + t12 * t286 + t233 * t289 + (t230 * t290 + t232 * t294) * t231; t280 * t293 + m(7) * (t32 * t58 + t42 * t79 + t43 * t80) + t220 * t291 + t11 * t286 + t6 * t287 + t5 * t288 + t222 * t290; m(7) * (t33 * t58 + t48 * t79 + t49 * t80) + 0.2e1 * t289 * t224 + 0.2e1 * t290 * t210 + 0.2e1 * t291 * t208; m(7) * (t184 * t80 + t186 * t79 + t212 * t58); t210 * t3 + t208 * t1 + t224 * t9 + m(7) * (t58 ^ 2 + t79 ^ 2 + t80 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
