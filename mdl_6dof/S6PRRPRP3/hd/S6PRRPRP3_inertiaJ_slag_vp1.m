% Calculate joint inertia matrix for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:14
% EndTime: 2019-03-08 21:34:24
% DurationCPUTime: 4.47s
% Computational Cost: add. (20098->557), mult. (41723->798), div. (0->0), fcn. (53333->12), ass. (0->253)
t281 = rSges(7,1) + pkin(5);
t280 = rSges(7,3) + qJ(6);
t279 = m(5) + m(6) + m(7);
t278 = cos(qJ(3));
t220 = cos(pkin(11));
t277 = pkin(4) * t220;
t222 = cos(pkin(6));
t218 = sin(pkin(10));
t226 = cos(qJ(2));
t267 = t226 * t218;
t221 = cos(pkin(10));
t225 = sin(qJ(2));
t269 = t221 * t225;
t206 = t222 * t269 + t267;
t224 = sin(qJ(3));
t219 = sin(pkin(6));
t238 = t219 * t278;
t195 = t206 * t224 + t221 * t238;
t271 = t219 * t224;
t196 = t206 * t278 - t221 * t271;
t266 = t226 * t221;
t268 = t225 * t218;
t205 = -t222 * t266 + t268;
t217 = sin(pkin(11));
t275 = t205 * t217;
t100 = pkin(4) * t275 + pkin(9) * t195 + t196 * t277;
t162 = pkin(3) * t196 + qJ(4) * t195;
t207 = t222 * t267 + t269;
t153 = t207 * t162;
t276 = t207 * t100 + t153;
t274 = t207 * t217;
t273 = t218 * t219;
t272 = t219 * t221;
t270 = t219 * t226;
t208 = -t222 * t268 + t266;
t197 = t208 * t224 - t218 * t238;
t198 = t208 * t278 + t218 * t271;
t101 = pkin(4) * t274 + pkin(9) * t197 + t198 * t277;
t163 = pkin(3) * t198 + qJ(4) * t197;
t264 = -t101 - t163;
t253 = pkin(11) + qJ(5);
t216 = sin(t253);
t237 = cos(t253);
t164 = t196 * t216 - t205 * t237;
t165 = t196 * t237 + t205 * t216;
t263 = rSges(7,2) * t195 + t280 * t164 + t281 * t165;
t166 = t198 * t216 - t207 * t237;
t167 = t198 * t237 + t207 * t216;
t262 = rSges(7,2) * t197 + t280 * t166 + t281 * t167;
t170 = -t198 * t217 + t207 * t220;
t171 = t198 * t220 + t274;
t127 = rSges(5,1) * t171 + rSges(5,2) * t170 + rSges(5,3) * t197;
t261 = -t127 - t163;
t210 = t222 * t224 + t225 * t238;
t191 = t210 * t216 + t237 * t270;
t192 = t210 * t237 - t216 * t270;
t209 = -t222 * t278 + t225 * t271;
t260 = rSges(7,2) * t209 + t280 * t191 + t281 * t192;
t242 = t217 * t270;
t140 = -pkin(4) * t242 + pkin(9) * t209 + t210 * t277;
t190 = t210 * pkin(3) + t209 * qJ(4);
t259 = -t140 - t190;
t193 = -t210 * t217 - t220 * t270;
t194 = t210 * t220 - t242;
t152 = rSges(5,1) * t194 + rSges(5,2) * t193 + rSges(5,3) * t209;
t258 = -t152 - t190;
t257 = t162 * t270 + t205 * t190;
t189 = pkin(2) * t208 + pkin(8) * t207;
t187 = t222 * t189;
t256 = t222 * t163 + t187;
t188 = pkin(2) * t206 + pkin(8) * t205;
t255 = -t162 - t188;
t254 = t188 * t273 + t189 * t272;
t102 = Icges(7,5) * t165 + Icges(7,6) * t195 + Icges(7,3) * t164;
t106 = Icges(7,4) * t165 + Icges(7,2) * t195 + Icges(7,6) * t164;
t110 = Icges(7,1) * t165 + Icges(7,4) * t195 + Icges(7,5) * t164;
t44 = t102 * t164 + t106 * t195 + t110 * t165;
t103 = Icges(7,5) * t167 + Icges(7,6) * t197 + Icges(7,3) * t166;
t107 = Icges(7,4) * t167 + Icges(7,2) * t197 + Icges(7,6) * t166;
t111 = Icges(7,1) * t167 + Icges(7,4) * t197 + Icges(7,5) * t166;
t45 = t103 * t164 + t107 * t195 + t111 * t165;
t132 = Icges(7,5) * t192 + Icges(7,6) * t209 + Icges(7,3) * t191;
t134 = Icges(7,4) * t192 + Icges(7,2) * t209 + Icges(7,6) * t191;
t136 = Icges(7,1) * t192 + Icges(7,4) * t209 + Icges(7,5) * t191;
t66 = t132 * t164 + t134 * t195 + t136 * t165;
t1 = t195 * t44 + t197 * t45 + t209 * t66;
t104 = Icges(6,5) * t165 - Icges(6,6) * t164 + Icges(6,3) * t195;
t108 = Icges(6,4) * t165 - Icges(6,2) * t164 + Icges(6,6) * t195;
t112 = Icges(6,1) * t165 - Icges(6,4) * t164 + Icges(6,5) * t195;
t46 = t104 * t195 - t108 * t164 + t112 * t165;
t105 = Icges(6,5) * t167 - Icges(6,6) * t166 + Icges(6,3) * t197;
t109 = Icges(6,4) * t167 - Icges(6,2) * t166 + Icges(6,6) * t197;
t113 = Icges(6,1) * t167 - Icges(6,4) * t166 + Icges(6,5) * t197;
t47 = t105 * t195 - t109 * t164 + t113 * t165;
t133 = Icges(6,5) * t192 - Icges(6,6) * t191 + Icges(6,3) * t209;
t135 = Icges(6,4) * t192 - Icges(6,2) * t191 + Icges(6,6) * t209;
t137 = Icges(6,1) * t192 - Icges(6,4) * t191 + Icges(6,5) * t209;
t67 = t133 * t195 - t135 * t164 + t137 * t165;
t2 = t195 * t46 + t197 * t47 + t209 * t67;
t252 = -t2 / 0.2e1 - t1 / 0.2e1;
t48 = t102 * t166 + t106 * t197 + t110 * t167;
t49 = t103 * t166 + t107 * t197 + t111 * t167;
t68 = t132 * t166 + t134 * t197 + t136 * t167;
t3 = t195 * t48 + t197 * t49 + t209 * t68;
t50 = t104 * t197 - t108 * t166 + t112 * t167;
t51 = t105 * t197 - t109 * t166 + t113 * t167;
t69 = t133 * t197 - t135 * t166 + t137 * t167;
t4 = t195 * t50 + t197 * t51 + t209 * t69;
t251 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t44 * t205 + t45 * t207 - t270 * t66;
t6 = t46 * t205 + t47 * t207 - t270 * t67;
t250 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t48 * t205 + t49 * t207 - t270 * t68;
t8 = t50 * t205 + t51 * t207 - t270 * t69;
t249 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t222 * t67 + (t218 * t47 - t221 * t46) * t219;
t9 = t222 * t66 + (t218 * t45 - t221 * t44) * t219;
t248 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t222 * t68 + (t218 * t49 - t221 * t48) * t219;
t12 = t222 * t69 + (t218 * t51 - t221 * t50) * t219;
t247 = t11 / 0.2e1 + t12 / 0.2e1;
t56 = t102 * t191 + t106 * t209 + t110 * t192;
t57 = t103 * t191 + t107 * t209 + t111 * t192;
t77 = t132 * t191 + t134 * t209 + t136 * t192;
t17 = t195 * t56 + t197 * t57 + t209 * t77;
t58 = t104 * t209 - t108 * t191 + t112 * t192;
t59 = t105 * t209 - t109 * t191 + t113 * t192;
t78 = t133 * t209 - t135 * t191 + t137 * t192;
t18 = t195 * t58 + t197 * t59 + t209 * t78;
t246 = t18 / 0.2e1 + t17 / 0.2e1;
t19 = t56 * t205 + t57 * t207 - t270 * t77;
t20 = t58 * t205 + t59 * t207 - t270 * t78;
t245 = t19 / 0.2e1 + t20 / 0.2e1;
t21 = t222 * t77 + (t218 * t57 - t221 * t56) * t219;
t22 = t222 * t78 + (t218 * t59 - t221 * t58) * t219;
t244 = t22 / 0.2e1 + t21 / 0.2e1;
t243 = t222 * t101 + t256;
t241 = -t100 + t255;
t117 = rSges(6,1) * t167 - rSges(6,2) * t166 + rSges(6,3) * t197;
t240 = -t117 + t264;
t139 = rSges(6,1) * t192 - rSges(6,2) * t191 + rSges(6,3) * t209;
t239 = -t139 + t259;
t184 = t210 * rSges(4,1) - t209 * rSges(4,2) - rSges(4,3) * t270;
t211 = (pkin(2) * t225 - pkin(8) * t226) * t219;
t236 = (-t184 - t211) * t219;
t234 = t100 * t270 + t205 * t140 + t257;
t233 = -t262 + t264;
t232 = t259 - t260;
t231 = t162 * t273 + t163 * t272 + t254;
t230 = (-t211 + t258) * t219;
t229 = (-t211 + t239) * t219;
t228 = t100 * t273 + t101 * t272 + t231;
t227 = (-t211 + t232) * t219;
t204 = t222 * rSges(3,3) + (rSges(3,1) * t225 + rSges(3,2) * t226) * t219;
t203 = Icges(3,5) * t222 + (Icges(3,1) * t225 + Icges(3,4) * t226) * t219;
t202 = Icges(3,6) * t222 + (Icges(3,4) * t225 + Icges(3,2) * t226) * t219;
t201 = Icges(3,3) * t222 + (Icges(3,5) * t225 + Icges(3,6) * t226) * t219;
t183 = Icges(4,1) * t210 - Icges(4,4) * t209 - Icges(4,5) * t270;
t182 = Icges(4,4) * t210 - Icges(4,2) * t209 - Icges(4,6) * t270;
t181 = Icges(4,5) * t210 - Icges(4,6) * t209 - Icges(4,3) * t270;
t180 = rSges(3,1) * t208 - rSges(3,2) * t207 + rSges(3,3) * t273;
t179 = rSges(3,1) * t206 - rSges(3,2) * t205 - rSges(3,3) * t272;
t178 = Icges(3,1) * t208 - Icges(3,4) * t207 + Icges(3,5) * t273;
t177 = Icges(3,1) * t206 - Icges(3,4) * t205 - Icges(3,5) * t272;
t176 = Icges(3,4) * t208 - Icges(3,2) * t207 + Icges(3,6) * t273;
t175 = Icges(3,4) * t206 - Icges(3,2) * t205 - Icges(3,6) * t272;
t174 = Icges(3,5) * t208 - Icges(3,6) * t207 + Icges(3,3) * t273;
t173 = Icges(3,5) * t206 - Icges(3,6) * t205 - Icges(3,3) * t272;
t169 = t196 * t220 + t275;
t168 = -t196 * t217 + t205 * t220;
t156 = -t179 * t222 - t204 * t272;
t155 = t180 * t222 - t204 * t273;
t151 = rSges(4,1) * t198 - rSges(4,2) * t197 + rSges(4,3) * t207;
t150 = rSges(4,1) * t196 - rSges(4,2) * t195 + rSges(4,3) * t205;
t149 = Icges(5,1) * t194 + Icges(5,4) * t193 + Icges(5,5) * t209;
t148 = Icges(5,4) * t194 + Icges(5,2) * t193 + Icges(5,6) * t209;
t147 = Icges(5,5) * t194 + Icges(5,6) * t193 + Icges(5,3) * t209;
t146 = Icges(4,1) * t198 - Icges(4,4) * t197 + Icges(4,5) * t207;
t145 = Icges(4,1) * t196 - Icges(4,4) * t195 + Icges(4,5) * t205;
t144 = Icges(4,4) * t198 - Icges(4,2) * t197 + Icges(4,6) * t207;
t143 = Icges(4,4) * t196 - Icges(4,2) * t195 + Icges(4,6) * t205;
t142 = Icges(4,5) * t198 - Icges(4,6) * t197 + Icges(4,3) * t207;
t141 = Icges(4,5) * t196 - Icges(4,6) * t195 + Icges(4,3) * t205;
t131 = (t179 * t218 + t180 * t221) * t219;
t126 = rSges(5,1) * t169 + rSges(5,2) * t168 + rSges(5,3) * t195;
t125 = Icges(5,1) * t171 + Icges(5,4) * t170 + Icges(5,5) * t197;
t124 = Icges(5,1) * t169 + Icges(5,4) * t168 + Icges(5,5) * t195;
t123 = Icges(5,4) * t171 + Icges(5,2) * t170 + Icges(5,6) * t197;
t122 = Icges(5,4) * t169 + Icges(5,2) * t168 + Icges(5,6) * t195;
t121 = Icges(5,5) * t171 + Icges(5,6) * t170 + Icges(5,3) * t197;
t120 = Icges(5,5) * t169 + Icges(5,6) * t168 + Icges(5,3) * t195;
t119 = -t151 * t270 - t207 * t184;
t118 = t150 * t270 + t205 * t184;
t115 = rSges(6,1) * t165 - rSges(6,2) * t164 + rSges(6,3) * t195;
t94 = -t181 * t270 - t209 * t182 + t210 * t183;
t93 = t150 * t207 - t151 * t205;
t92 = (-t150 - t188) * t222 + t221 * t236;
t91 = t151 * t222 + t218 * t236 + t187;
t90 = t181 * t207 - t182 * t197 + t183 * t198;
t89 = t181 * t205 - t182 * t195 + t183 * t196;
t88 = (t150 * t218 + t151 * t221) * t219 + t254;
t87 = -t142 * t270 - t209 * t144 + t210 * t146;
t86 = -t141 * t270 - t209 * t143 + t210 * t145;
t85 = t117 * t209 - t139 * t197;
t84 = -t115 * t209 + t139 * t195;
t83 = t147 * t209 + t148 * t193 + t149 * t194;
t82 = t142 * t207 - t144 * t197 + t146 * t198;
t81 = t141 * t207 - t143 * t197 + t145 * t198;
t80 = t142 * t205 - t144 * t195 + t146 * t196;
t79 = t141 * t205 - t143 * t195 + t145 * t196;
t76 = t115 * t197 - t117 * t195;
t75 = t207 * t258 + t261 * t270;
t74 = t126 * t270 + t205 * t152 + t257;
t73 = t147 * t197 + t148 * t170 + t149 * t171;
t72 = t147 * t195 + t148 * t168 + t149 * t169;
t71 = (-t126 + t255) * t222 + t221 * t230;
t70 = t127 * t222 + t218 * t230 + t256;
t65 = t126 * t207 + t205 * t261 + t153;
t64 = t121 * t209 + t123 * t193 + t125 * t194;
t63 = t120 * t209 + t122 * t193 + t124 * t194;
t62 = -t197 * t260 + t209 * t262;
t61 = t195 * t260 - t209 * t263;
t60 = (t126 * t218 + t127 * t221) * t219 + t231;
t55 = t121 * t197 + t123 * t170 + t125 * t171;
t54 = t120 * t197 + t122 * t170 + t124 * t171;
t53 = t121 * t195 + t123 * t168 + t125 * t169;
t52 = t120 * t195 + t122 * t168 + t124 * t169;
t43 = t207 * t239 + t240 * t270;
t42 = t115 * t270 + t205 * t139 + t234;
t41 = -t195 * t262 + t197 * t263;
t40 = (-t115 + t241) * t222 + t221 * t229;
t39 = t117 * t222 + t218 * t229 + t243;
t38 = t222 * t94 + (t218 * t87 - t221 * t86) * t219;
t37 = t86 * t205 + t87 * t207 - t270 * t94;
t36 = t115 * t207 + t205 * t240 + t276;
t35 = (t115 * t218 + t117 * t221) * t219 + t228;
t34 = t207 * t232 + t233 * t270;
t33 = t205 * t260 + t263 * t270 + t234;
t32 = t222 * t90 + (t218 * t82 - t221 * t81) * t219;
t31 = t222 * t89 + (t218 * t80 - t221 * t79) * t219;
t30 = (t241 - t263) * t222 + t221 * t227;
t29 = t218 * t227 + t222 * t262 + t243;
t28 = t81 * t205 + t82 * t207 - t270 * t90;
t27 = t79 * t205 + t80 * t207 - t270 * t89;
t26 = t205 * t233 + t207 * t263 + t276;
t25 = (t218 * t263 + t221 * t262) * t219 + t228;
t24 = t222 * t83 + (t218 * t64 - t221 * t63) * t219;
t23 = t63 * t205 + t64 * t207 - t270 * t83;
t16 = t222 * t73 + (t218 * t55 - t221 * t54) * t219;
t15 = t222 * t72 + (t218 * t53 - t221 * t52) * t219;
t14 = t54 * t205 + t55 * t207 - t270 * t73;
t13 = t52 * t205 + t53 * t207 - t270 * t72;
t95 = [m(2) + m(3) + m(4) + t279; m(3) * t131 + m(4) * t88 + m(5) * t60 + m(6) * t35 + m(7) * t25; m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t35 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t60 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(4) * (t88 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(3) * (t131 ^ 2 + t155 ^ 2 + t156 ^ 2) + (t11 + t12 + t16 + t32 + (t174 * t273 - t176 * t207 + t178 * t208) * t273) * t273 + (-t10 - t9 - t31 - t15 + (-t173 * t272 - t175 * t205 + t177 * t206) * t272 + (-t173 * t273 + t174 * t272 + t175 * t207 + t176 * t205 - t177 * t208 - t178 * t206) * t273) * t272 + ((t201 * t273 - t202 * t207 + t203 * t208) * t273 - (-t201 * t272 - t202 * t205 + t203 * t206) * t272 + t22 + t21 + t38 + t24 + ((t176 * t226 + t178 * t225) * t218 - (t175 * t226 + t177 * t225) * t221) * t219 ^ 2 + ((-t173 * t221 + t174 * t218 + t202 * t226 + t203 * t225) * t219 + t222 * t201) * t222) * t222; m(4) * t93 + m(5) * t65 + m(6) * t36 + m(7) * t26; (t23 / 0.2e1 + t37 / 0.2e1 + t245) * t222 + (t16 / 0.2e1 + t32 / 0.2e1 + t247) * t207 + (t15 / 0.2e1 + t31 / 0.2e1 + t248) * t205 + m(4) * (t118 * t92 + t119 * t91 + t88 * t93) + m(7) * (t25 * t26 + t29 * t34 + t30 * t33) + m(6) * (t35 * t36 + t39 * t43 + t40 * t42) + m(5) * (t60 * t65 + t70 * t75 + t71 * t74) + ((-t24 / 0.2e1 - t38 / 0.2e1 - t244) * t226 + (-t13 / 0.2e1 - t27 / 0.2e1 - t250) * t221 + (t14 / 0.2e1 + t28 / 0.2e1 + t249) * t218) * t219; (-t19 - t20 - t23 - t37) * t270 + (t7 + t8 + t28 + t14) * t207 + (t5 + t6 + t27 + t13) * t205 + m(7) * (t26 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t36 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t65 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t118 ^ 2 + t119 ^ 2 + t93 ^ 2); t209 * t279; m(7) * (t195 * t29 + t197 * t30 + t209 * t25) + m(6) * (t195 * t39 + t197 * t40 + t209 * t35) + m(5) * (t195 * t70 + t197 * t71 + t209 * t60); m(7) * (t195 * t34 + t197 * t33 + t209 * t26) + m(6) * (t195 * t43 + t197 * t42 + t209 * t36) + m(5) * (t195 * t75 + t197 * t74 + t209 * t65); (t195 ^ 2 + t197 ^ 2 + t209 ^ 2) * t279; m(6) * t76 + m(7) * t41; t246 * t222 + t244 * t209 + t247 * t197 + t248 * t195 + m(7) * (t25 * t41 + t29 * t62 + t30 * t61) + m(6) * (t35 * t76 + t39 * t85 + t40 * t84) + (t218 * t251 + t221 * t252) * t219; -t246 * t270 + t245 * t209 + t251 * t207 - t252 * t205 + t249 * t197 + t250 * t195 + m(7) * (t26 * t41 + t33 * t61 + t34 * t62) + m(6) * (t36 * t76 + t42 * t84 + t43 * t85); m(6) * (t195 * t85 + t197 * t84 + t209 * t76) + m(7) * (t195 * t62 + t197 * t61 + t209 * t41); (t17 + t18) * t209 + (t4 + t3) * t197 + (t1 + t2) * t195 + m(7) * (t41 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t76 ^ 2 + t84 ^ 2 + t85 ^ 2); m(7) * t191; m(7) * (t164 * t29 + t166 * t30 + t191 * t25); m(7) * (t164 * t34 + t166 * t33 + t191 * t26); m(7) * (t164 * t195 + t166 * t197 + t191 * t209); m(7) * (t164 * t62 + t166 * t61 + t191 * t41); m(7) * (t164 ^ 2 + t166 ^ 2 + t191 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t95(1) t95(2) t95(4) t95(7) t95(11) t95(16); t95(2) t95(3) t95(5) t95(8) t95(12) t95(17); t95(4) t95(5) t95(6) t95(9) t95(13) t95(18); t95(7) t95(8) t95(9) t95(10) t95(14) t95(19); t95(11) t95(12) t95(13) t95(14) t95(15) t95(20); t95(16) t95(17) t95(18) t95(19) t95(20) t95(21);];
Mq  = res;
