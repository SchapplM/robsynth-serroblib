% Calculate joint inertia matrix for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:36:04
% EndTime: 2019-03-08 20:36:14
% DurationCPUTime: 4.72s
% Computational Cost: add. (28397->562), mult. (42522->818), div. (0->0), fcn. (54168->14), ass. (0->251)
t227 = sin(pkin(6));
t287 = t227 ^ 2;
t225 = sin(pkin(12));
t228 = cos(pkin(12));
t230 = cos(pkin(6));
t233 = sin(qJ(2));
t269 = t227 * t233;
t209 = -t225 * t269 + t228 * t230;
t272 = t225 * t230;
t210 = t228 * t269 + t272;
t235 = cos(qJ(2));
t268 = t227 * t235;
t171 = Icges(4,5) * t210 + Icges(4,6) * t209 - Icges(4,3) * t268;
t204 = Icges(3,6) * t230 + (Icges(3,4) * t233 + Icges(3,2) * t235) * t227;
t286 = -t204 + t171;
t226 = sin(pkin(11));
t229 = cos(pkin(11));
t267 = t230 * t233;
t212 = t226 * t235 + t229 * t267;
t257 = pkin(12) + qJ(4);
t221 = sin(t257);
t246 = cos(t257);
t240 = t227 * t246;
t193 = t212 * t221 + t229 * t240;
t285 = t193 / 0.2e1;
t214 = -t226 * t267 + t229 * t235;
t195 = t214 * t221 - t226 * t240;
t284 = t195 / 0.2e1;
t206 = t221 * t269 - t230 * t246;
t283 = t206 / 0.2e1;
t266 = t230 * t235;
t211 = t226 * t233 - t229 * t266;
t282 = t211 / 0.2e1;
t213 = t226 * t266 + t229 * t233;
t281 = t213 / 0.2e1;
t280 = t230 / 0.2e1;
t279 = pkin(3) * t228;
t234 = cos(qJ(5));
t278 = pkin(5) * t234;
t270 = t227 * t229;
t194 = t212 * t246 - t221 * t270;
t224 = qJ(5) + qJ(6);
t222 = sin(t224);
t223 = cos(t224);
t159 = -t194 * t222 + t211 * t223;
t160 = t194 * t223 + t211 * t222;
t107 = rSges(7,1) * t160 + rSges(7,2) * t159 + rSges(7,3) * t193;
t232 = sin(qJ(5));
t274 = t211 * t232;
t97 = pkin(5) * t274 + pkin(10) * t193 + t194 * t278;
t276 = t107 + t97;
t271 = t226 * t227;
t196 = t214 * t246 + t221 * t271;
t161 = -t196 * t222 + t213 * t223;
t162 = t196 * t223 + t213 * t222;
t108 = rSges(7,1) * t162 + rSges(7,2) * t161 + rSges(7,3) * t195;
t273 = t213 * t232;
t98 = pkin(5) * t273 + pkin(10) * t195 + t196 * t278;
t275 = t108 + t98;
t165 = -t196 * t232 + t213 * t234;
t166 = t196 * t234 + t273;
t116 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t195;
t157 = t196 * pkin(4) + t195 * pkin(9);
t264 = -t116 - t157;
t207 = t230 * t221 + t233 * t240;
t251 = t232 * t268;
t118 = -pkin(5) * t251 + pkin(10) * t206 + t207 * t278;
t189 = -t207 * t222 - t223 * t268;
t190 = t207 * t223 - t222 * t268;
t123 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t206;
t263 = t118 + t123;
t253 = t225 * t271;
t140 = pkin(3) * t253 + pkin(8) * t213 + t214 * t279;
t192 = pkin(2) * t214 + qJ(3) * t213;
t188 = t230 * t192;
t262 = t230 * t140 + t188;
t252 = t225 * t270;
t139 = -pkin(3) * t252 + pkin(8) * t211 + t212 * t279;
t191 = pkin(2) * t212 + qJ(3) * t211;
t261 = -t139 - t191;
t156 = t194 * pkin(4) + t193 * pkin(9);
t181 = t207 * pkin(4) + t206 * pkin(9);
t260 = t156 * t268 + t211 * t181;
t215 = (pkin(2) * t233 - qJ(3) * t235) * t227;
t259 = -pkin(3) * t272 - (-pkin(8) * t235 + t233 * t279) * t227 - t215;
t258 = t191 * t271 + t192 * t270;
t101 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t193;
t103 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t193;
t105 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t193;
t60 = t101 * t206 + t103 * t189 + t105 * t190;
t102 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t195;
t104 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t195;
t106 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t195;
t61 = t102 * t206 + t104 * t189 + t106 * t190;
t120 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t206;
t121 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t206;
t122 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t206;
t72 = t120 * t206 + t121 * t189 + t122 * t190;
t26 = t193 * t60 + t195 * t61 + t206 * t72;
t50 = t101 * t193 + t103 * t159 + t105 * t160;
t51 = t102 * t193 + t104 * t159 + t106 * t160;
t66 = t120 * t193 + t121 * t159 + t122 * t160;
t7 = t193 * t50 + t195 * t51 + t206 * t66;
t52 = t101 * t195 + t103 * t161 + t105 * t162;
t53 = t102 * t195 + t104 * t161 + t106 * t162;
t67 = t120 * t195 + t121 * t161 + t122 * t162;
t8 = t193 * t52 + t195 * t53 + t206 * t67;
t256 = t193 * t7 + t195 * t8 + t206 * t26;
t255 = -t157 - t275;
t254 = -m(4) - m(5) - m(6) - m(7);
t250 = t230 * t157 + t262;
t249 = -t156 + t261;
t248 = -t181 + t259;
t247 = -t268 / 0.2e1;
t245 = (-rSges(4,1) * t210 - rSges(4,2) * t209 + rSges(4,3) * t268 - t215) * t227;
t244 = t139 * t271 + t140 * t270 + t258;
t170 = rSges(5,1) * t207 - rSges(5,2) * t206 - rSges(5,3) * t268;
t243 = (-t170 + t259) * t227;
t13 = t211 * t50 + t213 * t51 - t268 * t66;
t14 = t211 * t52 + t213 * t53 - t268 * t67;
t28 = t211 * t60 + t213 * t61 - t268 * t72;
t242 = t13 * t285 + t14 * t284 + t26 * t247 + t28 * t283 + t8 * t281 + t7 * t282;
t15 = t230 * t66 + (t226 * t51 - t229 * t50) * t227;
t16 = t230 * t67 + (t226 * t53 - t229 * t52) * t227;
t30 = t230 * t72 + (t226 * t61 - t229 * t60) * t227;
t241 = t15 * t285 + t16 * t284 + t26 * t280 + t30 * t283 + t8 * t271 / 0.2e1 - t7 * t270 / 0.2e1;
t197 = -t207 * t232 - t234 * t268;
t198 = t207 * t234 - t251;
t130 = rSges(6,1) * t198 + rSges(6,2) * t197 + rSges(6,3) * t206;
t239 = (-t130 + t248) * t227;
t238 = t156 * t271 + t157 * t270 + t244;
t237 = (t248 - t263) * t227;
t208 = rSges(3,3) * t230 + (rSges(3,1) * t233 + rSges(3,2) * t235) * t227;
t205 = Icges(3,5) * t230 + (Icges(3,1) * t233 + Icges(3,4) * t235) * t227;
t203 = Icges(3,3) * t230 + (Icges(3,5) * t233 + Icges(3,6) * t235) * t227;
t202 = t214 * t228 + t253;
t201 = -t214 * t225 + t228 * t271;
t200 = t212 * t228 - t252;
t199 = -t212 * t225 - t228 * t270;
t184 = rSges(3,1) * t214 - rSges(3,2) * t213 + rSges(3,3) * t271;
t183 = rSges(3,1) * t212 - rSges(3,2) * t211 - rSges(3,3) * t270;
t179 = Icges(3,1) * t214 - Icges(3,4) * t213 + Icges(3,5) * t271;
t178 = Icges(3,1) * t212 - Icges(3,4) * t211 - Icges(3,5) * t270;
t177 = Icges(3,4) * t214 - Icges(3,2) * t213 + Icges(3,6) * t271;
t176 = Icges(3,4) * t212 - Icges(3,2) * t211 - Icges(3,6) * t270;
t175 = Icges(3,5) * t214 - Icges(3,6) * t213 + Icges(3,3) * t271;
t174 = Icges(3,5) * t212 - Icges(3,6) * t211 - Icges(3,3) * t270;
t173 = Icges(4,1) * t210 + Icges(4,4) * t209 - Icges(4,5) * t268;
t172 = Icges(4,4) * t210 + Icges(4,2) * t209 - Icges(4,6) * t268;
t169 = Icges(5,1) * t207 - Icges(5,4) * t206 - Icges(5,5) * t268;
t168 = Icges(5,4) * t207 - Icges(5,2) * t206 - Icges(5,6) * t268;
t167 = Icges(5,5) * t207 - Icges(5,6) * t206 - Icges(5,3) * t268;
t164 = t194 * t234 + t274;
t163 = -t194 * t232 + t211 * t234;
t155 = -t183 * t230 - t208 * t270;
t154 = t184 * t230 - t208 * t271;
t149 = rSges(4,1) * t202 + rSges(4,2) * t201 + rSges(4,3) * t213;
t148 = rSges(4,1) * t200 + rSges(4,2) * t199 + rSges(4,3) * t211;
t147 = Icges(4,1) * t202 + Icges(4,4) * t201 + Icges(4,5) * t213;
t146 = Icges(4,1) * t200 + Icges(4,4) * t199 + Icges(4,5) * t211;
t145 = Icges(4,4) * t202 + Icges(4,2) * t201 + Icges(4,6) * t213;
t144 = Icges(4,4) * t200 + Icges(4,2) * t199 + Icges(4,6) * t211;
t143 = Icges(4,5) * t202 + Icges(4,6) * t201 + Icges(4,3) * t213;
t142 = Icges(4,5) * t200 + Icges(4,6) * t199 + Icges(4,3) * t211;
t141 = t213 * t156;
t138 = rSges(5,1) * t196 - rSges(5,2) * t195 + rSges(5,3) * t213;
t137 = rSges(5,1) * t194 - rSges(5,2) * t193 + rSges(5,3) * t211;
t136 = Icges(5,1) * t196 - Icges(5,4) * t195 + Icges(5,5) * t213;
t135 = Icges(5,1) * t194 - Icges(5,4) * t193 + Icges(5,5) * t211;
t134 = Icges(5,4) * t196 - Icges(5,2) * t195 + Icges(5,6) * t213;
t133 = Icges(5,4) * t194 - Icges(5,2) * t193 + Icges(5,6) * t211;
t132 = Icges(5,5) * t196 - Icges(5,6) * t195 + Icges(5,3) * t213;
t131 = Icges(5,5) * t194 - Icges(5,6) * t193 + Icges(5,3) * t211;
t128 = Icges(6,1) * t198 + Icges(6,4) * t197 + Icges(6,5) * t206;
t127 = Icges(6,4) * t198 + Icges(6,2) * t197 + Icges(6,6) * t206;
t126 = Icges(6,5) * t198 + Icges(6,6) * t197 + Icges(6,3) * t206;
t119 = (t183 * t226 + t184 * t229) * t227;
t117 = t193 * t123;
t115 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t193;
t114 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t195;
t113 = Icges(6,1) * t164 + Icges(6,4) * t163 + Icges(6,5) * t193;
t112 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t195;
t111 = Icges(6,4) * t164 + Icges(6,2) * t163 + Icges(6,6) * t193;
t110 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t195;
t109 = Icges(6,5) * t164 + Icges(6,6) * t163 + Icges(6,3) * t193;
t100 = -t138 * t268 - t170 * t213;
t99 = t137 * t268 + t170 * t211;
t96 = t206 * t108;
t95 = t195 * t107;
t94 = (-t148 - t191) * t230 + t229 * t245;
t93 = t149 * t230 + t226 * t245 + t188;
t92 = -t167 * t268 - t168 * t206 + t169 * t207;
t91 = t137 * t213 - t138 * t211;
t90 = t167 * t213 - t168 * t195 + t169 * t196;
t89 = t167 * t211 - t168 * t193 + t169 * t194;
t88 = (t148 * t226 + t149 * t229) * t227 + t258;
t87 = t116 * t206 - t130 * t195;
t86 = -t115 * t206 + t130 * t193;
t85 = -t123 * t195 + t96;
t84 = -t107 * t206 + t117;
t83 = -t132 * t268 - t134 * t206 + t136 * t207;
t82 = -t131 * t268 - t133 * t206 + t135 * t207;
t81 = t132 * t213 - t134 * t195 + t136 * t196;
t80 = t131 * t213 - t133 * t195 + t135 * t196;
t79 = t132 * t211 - t134 * t193 + t136 * t194;
t78 = t131 * t211 - t133 * t193 + t135 * t194;
t77 = (-t137 + t261) * t230 + t229 * t243;
t76 = t138 * t230 + t226 * t243 + t262;
t75 = t126 * t206 + t127 * t197 + t128 * t198;
t74 = t115 * t195 - t116 * t193;
t73 = -t108 * t193 + t95;
t71 = t264 * t268 + (-t130 - t181) * t213;
t70 = t115 * t268 + t130 * t211 + t260;
t69 = t126 * t195 + t127 * t165 + t128 * t166;
t68 = t126 * t193 + t127 * t163 + t128 * t164;
t65 = (t137 * t226 + t138 * t229) * t227 + t244;
t64 = t115 * t213 + t211 * t264 + t141;
t63 = t110 * t206 + t112 * t197 + t114 * t198;
t62 = t109 * t206 + t111 * t197 + t113 * t198;
t59 = (-t115 + t249) * t230 + t229 * t239;
t58 = t116 * t230 + t226 * t239 + t250;
t57 = t110 * t195 + t112 * t165 + t114 * t166;
t56 = t109 * t195 + t111 * t165 + t113 * t166;
t55 = t110 * t193 + t112 * t163 + t114 * t164;
t54 = t109 * t193 + t111 * t163 + t113 * t164;
t49 = -t195 * t263 + t206 * t98 + t96;
t48 = t118 * t193 - t206 * t276 + t117;
t47 = t255 * t268 + (-t181 - t263) * t213;
t46 = t211 * t263 + t268 * t276 + t260;
t45 = (t115 * t226 + t116 * t229) * t227 + t238;
t44 = -t193 * t275 + t195 * t97 + t95;
t43 = t230 * t92 + (t226 * t83 - t229 * t82) * t227;
t42 = t211 * t82 + t213 * t83 - t268 * t92;
t41 = (t249 - t276) * t230 + t229 * t237;
t40 = t226 * t237 + t230 * t275 + t250;
t39 = t211 * t255 + t213 * t276 + t141;
t38 = t230 * t90 + (t226 * t81 - t229 * t80) * t227;
t37 = t230 * t89 + (t226 * t79 - t229 * t78) * t227;
t36 = t211 * t80 + t213 * t81 - t268 * t90;
t35 = t211 * t78 + t213 * t79 - t268 * t89;
t34 = (t226 * t276 + t229 * t275) * t227 + t238;
t33 = t230 * t75 + (t226 * t63 - t229 * t62) * t227;
t32 = t211 * t62 + t213 * t63 - t268 * t75;
t31 = t193 * t62 + t195 * t63 + t206 * t75;
t23 = t230 * t69 + (t226 * t57 - t229 * t56) * t227;
t22 = t230 * t68 + (t226 * t55 - t229 * t54) * t227;
t20 = t211 * t56 + t213 * t57 - t268 * t69;
t19 = t211 * t54 + t213 * t55 - t268 * t68;
t18 = t193 * t56 + t195 * t57 + t206 * t69;
t17 = t193 * t54 + t195 * t55 + t206 * t68;
t1 = [m(2) + m(3) - t254; m(3) * t119 + m(4) * t88 + m(5) * t65 + m(6) * t45 + m(7) * t34; m(7) * (t34 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t45 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t65 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(4) * (t88 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(3) * (t119 ^ 2 + t154 ^ 2 + t155 ^ 2) + (t16 + t23 + t38 + ((t213 * t143 + t201 * t145 + t202 * t147) * t226 - (t213 * t142 + t201 * t144 + t202 * t146) * t229) * t227 + (t175 * t271 - t213 * t177 + t214 * t179) * t271) * t271 + (-t15 - t22 - t37 - ((t211 * t143 + t199 * t145 + t200 * t147) * t226 - (t211 * t142 + t199 * t144 + t200 * t146) * t229) * t227 + (-t174 * t270 - t211 * t176 + t212 * t178) * t270 + (-t174 * t271 + t175 * t270 + t213 * t176 + t211 * t177 - t214 * t178 - t212 * t179) * t271) * t270 + (t30 + t33 + t43 + ((t177 * t235 + t179 * t233) * t226 - (t176 * t235 + t178 * t233) * t229) * t287 + ((-t174 * t229 + t175 * t226 + t204 * t235 + t205 * t233) * t227 - t171 * t268 + t209 * t172 + t210 * t173 + t230 * t203) * t230 + (-t143 * t268 + t145 * t209 + t210 * t147 + t201 * t172 + t202 * t173 + t203 * t271 + t214 * t205 + t213 * t286) * t271 + (t142 * t268 - t209 * t144 - t210 * t146 - t199 * t172 - t200 * t173 + t203 * t270 - t212 * t205 - t211 * t286) * t270) * t230; t254 * t268; m(7) * (t211 * t40 + t213 * t41 - t268 * t34) + m(6) * (t211 * t58 + t213 * t59 - t268 * t45) + m(5) * (t211 * t76 + t213 * t77 - t268 * t65) + m(4) * (t211 * t93 + t213 * t94 - t268 * t88); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1 + m(4) / 0.2e1) * (t235 ^ 2 * t287 + t211 ^ 2 + t213 ^ 2); m(5) * t91 + m(6) * t64 + m(7) * t39; (t28 / 0.2e1 + t32 / 0.2e1 + t42 / 0.2e1) * t230 + (t16 / 0.2e1 + t23 / 0.2e1 + t38 / 0.2e1) * t213 + (t15 / 0.2e1 + t22 / 0.2e1 + t37 / 0.2e1) * t211 + m(7) * (t34 * t39 + t40 * t47 + t41 * t46) + m(6) * (t45 * t64 + t58 * t71 + t59 * t70) + m(5) * (t100 * t76 + t65 * t91 + t77 * t99) + ((-t30 / 0.2e1 - t33 / 0.2e1 - t43 / 0.2e1) * t235 + (-t13 / 0.2e1 - t19 / 0.2e1 - t35 / 0.2e1) * t229 + (t14 / 0.2e1 + t20 / 0.2e1 + t36 / 0.2e1) * t226) * t227; m(5) * (t100 * t211 + t213 * t99 - t268 * t91) + m(6) * (t211 * t71 + t213 * t70 - t268 * t64) + m(7) * (t211 * t47 + t213 * t46 - t268 * t39); (-t28 - t32 - t42) * t268 + (t14 + t20 + t36) * t213 + (t13 + t19 + t35) * t211 + m(7) * (t39 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t64 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t100 ^ 2 + t91 ^ 2 + t99 ^ 2); m(6) * t74 + m(7) * t44; t23 * t284 + t33 * t283 + t31 * t280 + t22 * t285 + (-t229 * t17 / 0.2e1 + t226 * t18 / 0.2e1) * t227 + m(7) * (t34 * t44 + t40 * t49 + t41 * t48) + m(6) * (t45 * t74 + t58 * t87 + t59 * t86) + t241; m(6) * (t211 * t87 + t213 * t86 - t268 * t74) + m(7) * (t211 * t49 + t213 * t48 - t268 * t44); t32 * t283 + t20 * t284 + t18 * t281 + t17 * t282 + t19 * t285 + t31 * t247 + m(7) * (t39 * t44 + t46 * t48 + t47 * t49) + m(6) * (t64 * t74 + t70 * t86 + t71 * t87) + t242; t193 * t17 + t195 * t18 + t206 * t31 + m(7) * (t44 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t74 ^ 2 + t86 ^ 2 + t87 ^ 2) + t256; m(7) * t73; m(7) * (t34 * t73 + t40 * t85 + t41 * t84) + t241; m(7) * (t211 * t85 + t213 * t84 - t268 * t73); m(7) * (t39 * t73 + t46 * t84 + t47 * t85) + t242; m(7) * (t44 * t73 + t48 * t84 + t49 * t85) + t256; m(7) * (t73 ^ 2 + t84 ^ 2 + t85 ^ 2) + t256;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
