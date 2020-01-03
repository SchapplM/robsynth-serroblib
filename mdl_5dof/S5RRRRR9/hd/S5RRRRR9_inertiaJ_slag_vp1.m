% Calculate joint inertia matrix for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:45
% EndTime: 2019-12-31 22:27:55
% DurationCPUTime: 3.54s
% Computational Cost: add. (10887->442), mult. (12524->615), div. (0->0), fcn. (13670->10), ass. (0->229)
t221 = sin(qJ(2));
t299 = Icges(3,5) * t221;
t298 = t299 / 0.2e1;
t219 = qJ(3) + qJ(4);
t212 = qJ(5) + t219;
t207 = sin(t212);
t208 = cos(t212);
t222 = sin(qJ(1));
t224 = cos(qJ(2));
t225 = cos(qJ(1));
t267 = t224 * t225;
t159 = -t207 * t267 + t208 * t222;
t160 = t207 * t222 + t208 * t267;
t270 = t221 * t225;
t109 = t160 * rSges(6,1) + t159 * rSges(6,2) + rSges(6,3) * t270;
t223 = cos(qJ(3));
t209 = t223 * pkin(3) + pkin(2);
t211 = cos(t219);
t187 = pkin(4) * t211 + t209;
t210 = sin(t219);
t220 = sin(qJ(3));
t189 = pkin(3) * t220 + pkin(4) * t210;
t297 = t187 * t267 + t222 * t189 + t109;
t296 = t222 ^ 2;
t295 = t225 ^ 2;
t226 = -pkin(8) - pkin(7);
t271 = t221 * t222;
t268 = t222 * t224;
t157 = -t207 * t268 - t208 * t225;
t158 = -t207 * t225 + t208 * t268;
t102 = Icges(6,5) * t158 + Icges(6,6) * t157 + Icges(6,3) * t271;
t104 = Icges(6,4) * t158 + Icges(6,2) * t157 + Icges(6,6) * t271;
t106 = Icges(6,1) * t158 + Icges(6,4) * t157 + Icges(6,5) * t271;
t35 = t102 * t271 + t104 * t157 + t106 * t158;
t103 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t270;
t105 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t270;
t107 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t270;
t36 = t103 * t271 + t105 * t157 + t107 * t158;
t145 = -Icges(6,3) * t224 + (Icges(6,5) * t208 - Icges(6,6) * t207) * t221;
t146 = -Icges(6,6) * t224 + (Icges(6,4) * t208 - Icges(6,2) * t207) * t221;
t147 = -Icges(6,5) * t224 + (Icges(6,1) * t208 - Icges(6,4) * t207) * t221;
t64 = t145 * t271 + t146 * t157 + t147 * t158;
t5 = -t224 * t64 + (t222 * t35 + t225 * t36) * t221;
t37 = t102 * t270 + t104 * t159 + t106 * t160;
t38 = t103 * t270 + t105 * t159 + t107 * t160;
t65 = t145 * t270 + t146 * t159 + t147 * t160;
t6 = -t224 * t65 + (t222 * t37 + t225 * t38) * t221;
t294 = t6 * t270 + t5 * t271;
t293 = t222 / 0.2e1;
t292 = -t224 / 0.2e1;
t291 = -t225 / 0.2e1;
t290 = t225 / 0.2e1;
t194 = rSges(3,1) * t221 + rSges(3,2) * t224;
t289 = m(3) * t194;
t288 = pkin(2) * t224;
t287 = pkin(7) * t221;
t286 = -pkin(2) + t209;
t138 = t221 * t208 * t147;
t276 = t207 * t146;
t77 = -t224 * t145 - t221 * t276 + t138;
t281 = t77 * t224;
t48 = -t102 * t224 + (-t104 * t207 + t106 * t208) * t221;
t49 = -t103 * t224 + (-t105 * t207 + t107 * t208) * t221;
t13 = -t281 + (t222 * t48 + t225 * t49) * t221;
t173 = -t210 * t268 - t211 * t225;
t174 = -t210 * t225 + t211 * t268;
t111 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t271;
t113 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t271;
t115 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t271;
t52 = -t111 * t224 + (-t113 * t210 + t115 * t211) * t221;
t175 = -t210 * t267 + t211 * t222;
t176 = t210 * t222 + t211 * t267;
t112 = Icges(5,5) * t176 + Icges(5,6) * t175 + Icges(5,3) * t270;
t114 = Icges(5,4) * t176 + Icges(5,2) * t175 + Icges(5,6) * t270;
t116 = Icges(5,1) * t176 + Icges(5,4) * t175 + Icges(5,5) * t270;
t53 = -t112 * t224 + (-t114 * t210 + t116 * t211) * t221;
t153 = -Icges(5,5) * t224 + (Icges(5,1) * t211 - Icges(5,4) * t210) * t221;
t139 = t221 * t211 * t153;
t151 = -Icges(5,3) * t224 + (Icges(5,5) * t211 - Icges(5,6) * t210) * t221;
t152 = -Icges(5,6) * t224 + (Icges(5,4) * t211 - Icges(5,2) * t210) * t221;
t275 = t210 * t152;
t81 = -t224 * t151 - t221 * t275 + t139;
t285 = -t13 + t81 * t224 - (t222 * t52 + t225 * t53) * t221;
t284 = -t77 - t81;
t218 = -pkin(9) + t226;
t269 = t221 * t226;
t272 = t220 * t225;
t257 = pkin(3) * t272 + t222 * t269;
t259 = t187 - t209;
t277 = t189 * t225;
t93 = -t277 + (-t218 * t221 + t259 * t224) * t222 + t257;
t237 = -rSges(6,1) * t158 - rSges(6,2) * t157;
t108 = rSges(6,3) * t271 - t237;
t99 = t108 * t270;
t283 = t93 * t270 + t99;
t282 = t225 * rSges(3,3);
t254 = t218 - t226;
t273 = t220 * t222;
t258 = -pkin(3) * t273 - t209 * t267;
t280 = -t254 * t270 + t258 + t297;
t278 = Icges(3,4) * t224;
t166 = -Icges(4,6) * t224 + (Icges(4,4) * t223 - Icges(4,2) * t220) * t221;
t274 = t220 * t166;
t118 = t176 * rSges(5,1) + t175 * rSges(5,2) + rSges(5,3) * t270;
t230 = -t225 * t269 - t258;
t256 = pkin(2) * t267 + pkin(7) * t270;
t134 = t230 - t256;
t266 = -t118 - t134;
t133 = (t286 * t224 - t287) * t222 - t257;
t150 = (pkin(7) + t226) * t224 + t286 * t221;
t265 = t224 * t133 + t150 * t271;
t135 = t259 * t221 + t254 * t224;
t149 = -rSges(6,3) * t224 + (rSges(6,1) * t208 - rSges(6,2) * t207) * t221;
t264 = -t135 - t149;
t82 = t224 * t108 + t149 * t271;
t238 = -rSges(5,1) * t174 - rSges(5,2) * t173;
t117 = rSges(5,3) * t271 - t238;
t156 = -rSges(5,3) * t224 + (rSges(5,1) * t211 - rSges(5,2) * t210) * t221;
t86 = t224 * t117 + t156 * t271;
t263 = -t150 - t156;
t172 = -rSges(4,3) * t224 + (rSges(4,1) * t223 - rSges(4,2) * t220) * t221;
t197 = t221 * pkin(2) - t224 * pkin(7);
t262 = -t172 - t197;
t260 = t296 * (t287 + t288) + t225 * t256;
t255 = t225 * pkin(1) + t222 * pkin(6);
t182 = -t220 * t268 - t223 * t225;
t183 = t223 * t268 - t272;
t125 = Icges(4,5) * t183 + Icges(4,6) * t182 + Icges(4,3) * t271;
t127 = Icges(4,4) * t183 + Icges(4,2) * t182 + Icges(4,6) * t271;
t129 = Icges(4,1) * t183 + Icges(4,4) * t182 + Icges(4,5) * t271;
t60 = -t125 * t224 + (-t127 * t220 + t129 * t223) * t221;
t163 = -Icges(4,3) * t224 + (Icges(4,5) * t223 - Icges(4,6) * t220) * t221;
t169 = -Icges(4,5) * t224 + (Icges(4,1) * t223 - Icges(4,4) * t220) * t221;
t78 = t163 * t271 + t166 * t182 + t169 * t183;
t253 = t60 / 0.2e1 + t78 / 0.2e1;
t184 = -t220 * t267 + t222 * t223;
t185 = t223 * t267 + t273;
t126 = Icges(4,5) * t185 + Icges(4,6) * t184 + Icges(4,3) * t270;
t128 = Icges(4,4) * t185 + Icges(4,2) * t184 + Icges(4,6) * t270;
t130 = Icges(4,1) * t185 + Icges(4,4) * t184 + Icges(4,5) * t270;
t61 = -t126 * t224 + (-t128 * t220 + t130 * t223) * t221;
t79 = t163 * t270 + t166 * t184 + t169 * t185;
t252 = t79 / 0.2e1 + t61 / 0.2e1;
t251 = -t134 - t280;
t42 = t111 * t271 + t113 * t173 + t115 * t174;
t43 = t112 * t271 + t114 * t173 + t116 * t174;
t70 = t151 * t271 + t152 * t173 + t153 * t174;
t11 = -t224 * t70 + (t222 * t42 + t225 * t43) * t221;
t44 = t111 * t270 + t113 * t175 + t115 * t176;
t45 = t112 * t270 + t114 * t175 + t116 * t176;
t71 = t151 * t270 + t152 * t175 + t153 * t176;
t12 = -t224 * t71 + (t222 * t44 + t225 * t45) * t221;
t250 = t11 * t271 + t12 * t270 + t294;
t249 = -t150 + t264;
t248 = -t197 + t263;
t132 = t185 * rSges(4,1) + t184 * rSges(4,2) + rSges(4,3) * t270;
t247 = t271 / 0.2e1;
t246 = t270 / 0.2e1;
t245 = (t48 + t64) * t247 + (t49 + t65) * t246;
t40 = t135 * t271 + t224 * t93 + t82;
t244 = -t224 * t13 + t294;
t243 = t222 * t133 + t225 * t134 + t260;
t242 = -t197 + t249;
t19 = t222 * t36 - t225 * t35;
t20 = t222 * t38 - t225 * t37;
t241 = t19 * t247 + t20 * t246 + t5 * t291 + t6 * t293 + (t49 * t222 - t48 * t225) * t292;
t240 = rSges(3,1) * t224 - rSges(3,2) * t221;
t239 = -rSges(4,1) * t183 - rSges(4,2) * t182;
t235 = -Icges(3,2) * t221 + t278;
t234 = Icges(3,5) * t224 - Icges(3,6) * t221;
t231 = rSges(3,1) * t267 - rSges(3,2) * t270 + t222 * rSges(3,3);
t229 = t285 * t224 + t250;
t228 = t245 + (t52 + t70) * t247 + (t53 + t71) * t246;
t23 = t222 * t43 - t225 * t42;
t24 = t222 * t45 - t225 * t44;
t227 = t11 * t291 + t12 * t293 + t23 * t247 + t24 * t246 + t241 + (t53 * t222 - t52 * t225) * t292;
t216 = t225 * pkin(6);
t196 = rSges(2,1) * t225 - rSges(2,2) * t222;
t195 = -rSges(2,1) * t222 - rSges(2,2) * t225;
t191 = Icges(3,6) * t224 + t299;
t165 = Icges(3,3) * t222 + t225 * t234;
t164 = -Icges(3,3) * t225 + t222 * t234;
t148 = t221 * t223 * t169;
t144 = t231 + t255;
t143 = t282 + t216 + (-pkin(1) - t240) * t222;
t137 = t262 * t225;
t136 = t262 * t222;
t131 = rSges(4,3) * t271 - t239;
t120 = t225 * t231 + (t240 * t222 - t282) * t222;
t119 = t133 * t270;
t101 = t117 * t270;
t98 = t132 + t255 + t256;
t97 = t216 + (-t288 - pkin(1) + (-rSges(4,3) - pkin(7)) * t221) * t222 + t239;
t96 = t248 * t225;
t95 = t248 * t222;
t92 = -t132 * t224 - t172 * t270;
t91 = t131 * t224 + t172 * t271;
t89 = -t224 * t163 - t221 * t274 + t148;
t87 = -t118 * t224 - t156 * t270;
t85 = t230 + t118 + t255;
t84 = t216 + (-rSges(5,3) * t221 - t209 * t224 - pkin(1)) * t222 + t238 + t257;
t83 = -t109 * t224 - t149 * t270;
t80 = (t131 * t225 - t132 * t222) * t221;
t76 = -t218 * t270 + t255 + t297;
t75 = t277 + t216 + (-t187 * t224 - pkin(1) + (-rSges(6,3) + t218) * t221) * t222 + t237;
t74 = t242 * t225;
t73 = t242 * t222;
t72 = -t118 * t271 + t101;
t69 = -t109 * t271 + t99;
t66 = t131 * t222 + t132 * t225 + t260;
t59 = t266 * t224 + t263 * t270;
t58 = t86 + t265;
t57 = t126 * t270 + t128 * t184 + t130 * t185;
t56 = t125 * t270 + t127 * t184 + t129 * t185;
t55 = t126 * t271 + t128 * t182 + t130 * t183;
t54 = t125 * t271 + t127 * t182 + t129 * t183;
t41 = -t280 * t224 + t264 * t270;
t39 = t266 * t271 + t101 + t119;
t34 = t117 * t222 + t118 * t225 + t243;
t33 = -t280 * t271 + t283;
t32 = t251 * t224 + t249 * t270;
t31 = t40 + t265;
t30 = t222 * t57 - t225 * t56;
t29 = t222 * t55 - t225 * t54;
t27 = t251 * t271 + t119 + t283;
t25 = t280 * t225 + (t108 + t93) * t222 + t243;
t18 = -t224 * t79 + (t222 * t56 + t225 * t57) * t221;
t17 = -t224 * t78 + (t222 * t54 + t225 * t55) * t221;
t1 = [Icges(2,3) + t138 + t139 + t148 + (Icges(3,4) * t221 + Icges(3,2) * t224 - t145 - t151 - t163) * t224 + (Icges(3,1) * t221 - t274 - t275 - t276 + t278) * t221 + m(6) * (t75 ^ 2 + t76 ^ 2) + m(5) * (t84 ^ 2 + t85 ^ 2) + m(4) * (t97 ^ 2 + t98 ^ 2) + m(3) * (t143 ^ 2 + t144 ^ 2) + m(2) * (t195 ^ 2 + t196 ^ 2); m(6) * (t73 * t76 + t74 * t75) + m(5) * (t84 * t96 + t85 * t95) + m(4) * (t136 * t98 + t137 * t97) + (-t64 / 0.2e1 - t48 / 0.2e1 - t70 / 0.2e1 - t52 / 0.2e1 + (-Icges(3,6) * t225 + t222 * t235) * t292 + t225 * t298 - t143 * t289 + t191 * t290 - t253) * t225 + (t49 / 0.2e1 + t65 / 0.2e1 + t53 / 0.2e1 + t71 / 0.2e1 + (Icges(3,6) * t222 + t225 * t235) * t224 / 0.2e1 + t222 * t298 - t144 * t289 + t191 * t293 + t252) * t222; m(6) * (t25 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t34 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(4) * (t136 ^ 2 + t137 ^ 2 + t66 ^ 2) + m(3) * (t120 ^ 2 + (t295 + t296) * t194 ^ 2) + (-t295 * t164 - t19 - t23 - t29) * t225 + (t296 * t165 + t20 + t24 + t30 + (-t222 * t164 + t225 * t165) * t225) * t222; (-t89 + t284) * t224 + m(6) * (t31 * t75 + t32 * t76) + m(5) * (t58 * t84 + t59 * t85) + m(4) * (t91 * t97 + t92 * t98) + (t253 * t222 + t252 * t225) * t221 + t228; m(6) * (t25 * t27 + t31 * t74 + t32 * t73) + m(5) * (t39 * t34 + t58 * t96 + t59 * t95) + m(4) * (t136 * t92 + t137 * t91 + t66 * t80) + (t29 * t293 + t30 * t290) * t221 + t17 * t291 + (t61 * t222 - t60 * t225) * t292 + t18 * t293 + t227; (t89 * t224 + t285) * t224 + (t225 * t18 + t222 * t17 - t224 * (t60 * t222 + t61 * t225)) * t221 + m(6) * (t27 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t39 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t80 ^ 2 + t91 ^ 2 + t92 ^ 2) + t250; t284 * t224 + m(6) * (t40 * t75 + t41 * t76) + m(5) * (t84 * t86 + t85 * t87) + t228; m(6) * (t25 * t33 + t40 * t74 + t41 * t73) + m(5) * (t72 * t34 + t86 * t96 + t87 * t95) + t227; m(6) * (t27 * t33 + t31 * t40 + t32 * t41) + m(5) * (t72 * t39 + t58 * t86 + t59 * t87) + t229; m(6) * (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t72 ^ 2 + t86 ^ 2 + t87 ^ 2) + t229; -t281 + m(6) * (t75 * t82 + t76 * t83) + t245; m(6) * (t69 * t25 + t73 * t83 + t74 * t82) + t241; m(6) * (t69 * t27 + t31 * t82 + t32 * t83) + t244; m(6) * (t69 * t33 + t40 * t82 + t41 * t83) + t244; m(6) * (t69 ^ 2 + t82 ^ 2 + t83 ^ 2) + t244;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
