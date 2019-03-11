% Calculate joint inertia matrix for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:30
% EndTime: 2019-03-09 15:14:42
% DurationCPUTime: 4.17s
% Computational Cost: add. (9858->503), mult. (26752->701), div. (0->0), fcn. (33134->10), ass. (0->233)
t229 = sin(qJ(2));
t318 = Icges(3,5) * t229;
t317 = t318 / 0.2e1;
t228 = sin(qJ(3));
t231 = cos(qJ(3));
t233 = cos(qJ(1));
t230 = sin(qJ(1));
t232 = cos(qJ(2));
t286 = t230 * t232;
t203 = -t228 * t286 - t231 * t233;
t204 = -t228 * t233 + t231 * t286;
t294 = sin(pkin(6));
t295 = cos(pkin(10));
t254 = t295 * t294;
t242 = t229 * t254;
t296 = cos(pkin(6));
t256 = t296 * t295;
t293 = sin(pkin(10));
t150 = -t203 * t256 + t204 * t293 - t230 * t242;
t266 = t229 * t296;
t258 = t230 * t266;
t268 = t203 * t294;
t178 = -t268 + t258;
t315 = rSges(7,1) + pkin(5);
t316 = -t150 * rSges(7,2) - t315 * t178;
t314 = rSges(7,3) + qJ(6);
t285 = t232 * t233;
t205 = -t228 * t285 + t230 * t231;
t179 = -t205 * t294 + t233 * t266;
t274 = m(6) / 0.2e1 + m(7) / 0.2e1;
t313 = 0.2e1 * t274;
t187 = -Icges(4,3) * t232 + (Icges(4,5) * t231 - Icges(4,6) * t228) * t229;
t190 = -Icges(4,6) * t232 + (Icges(4,4) * t231 - Icges(4,2) * t228) * t229;
t193 = -Icges(4,5) * t232 + (Icges(4,1) * t231 - Icges(4,4) * t228) * t229;
t289 = t229 * t230;
t116 = t187 * t289 + t190 * t203 + t193 * t204;
t253 = t294 * t293;
t241 = t229 * t253;
t255 = t296 * t293;
t151 = t203 * t255 + t204 * t295 + t230 * t241;
t77 = Icges(5,5) * t151 - Icges(5,6) * t150 + Icges(5,3) * t178;
t83 = Icges(5,4) * t151 - Icges(5,2) * t150 + Icges(5,6) * t178;
t89 = Icges(5,1) * t151 - Icges(5,4) * t150 + Icges(5,5) * t178;
t19 = -t150 * t83 + t151 * t89 + t178 * t77;
t206 = t230 * t228 + t231 * t285;
t152 = -t205 * t256 + t206 * t293 - t233 * t242;
t153 = t205 * t255 + t206 * t295 + t233 * t241;
t78 = Icges(5,5) * t153 - Icges(5,6) * t152 + Icges(5,3) * t179;
t84 = Icges(5,4) * t153 - Icges(5,2) * t152 + Icges(5,6) * t179;
t90 = Icges(5,1) * t153 - Icges(5,4) * t152 + Icges(5,5) * t179;
t20 = -t150 * t84 + t151 * t90 + t178 * t78;
t73 = Icges(7,5) * t178 + Icges(7,6) * t150 + Icges(7,3) * t151;
t79 = Icges(7,4) * t178 + Icges(7,2) * t150 + Icges(7,6) * t151;
t85 = Icges(7,1) * t178 + Icges(7,4) * t150 + Icges(7,5) * t151;
t23 = t150 * t79 + t151 * t73 + t178 * t85;
t74 = Icges(7,5) * t179 + Icges(7,6) * t152 + Icges(7,3) * t153;
t80 = Icges(7,4) * t179 + Icges(7,2) * t152 + Icges(7,6) * t153;
t86 = Icges(7,1) * t179 + Icges(7,4) * t152 + Icges(7,5) * t153;
t24 = t150 * t80 + t151 * t74 + t178 * t86;
t75 = Icges(6,5) * t178 - Icges(6,6) * t151 + Icges(6,3) * t150;
t81 = Icges(6,4) * t178 - Icges(6,2) * t151 + Icges(6,6) * t150;
t87 = Icges(6,1) * t178 - Icges(6,4) * t151 + Icges(6,5) * t150;
t25 = t150 * t75 - t151 * t81 + t178 * t87;
t76 = Icges(6,5) * t179 - Icges(6,6) * t153 + Icges(6,3) * t152;
t82 = Icges(6,4) * t179 - Icges(6,2) * t153 + Icges(6,6) * t152;
t88 = Icges(6,1) * t179 - Icges(6,4) * t153 + Icges(6,5) * t152;
t26 = t150 * t76 - t151 * t82 + t178 * t88;
t172 = t232 * t254 + (t228 * t256 + t231 * t293) * t229;
t288 = t229 * t231;
t290 = t228 * t229;
t173 = -t232 * t253 - t255 * t290 + t288 * t295;
t202 = -t232 * t296 + t290 * t294;
t123 = Icges(5,5) * t173 - Icges(5,6) * t172 + Icges(5,3) * t202;
t126 = Icges(5,4) * t173 - Icges(5,2) * t172 + Icges(5,6) * t202;
t129 = Icges(5,1) * t173 - Icges(5,4) * t172 + Icges(5,5) * t202;
t45 = t123 * t178 - t126 * t150 + t129 * t151;
t121 = Icges(7,5) * t202 + Icges(7,6) * t172 + Icges(7,3) * t173;
t124 = Icges(7,4) * t202 + Icges(7,2) * t172 + Icges(7,6) * t173;
t127 = Icges(7,1) * t202 + Icges(7,4) * t172 + Icges(7,5) * t173;
t47 = t121 * t151 + t124 * t150 + t127 * t178;
t122 = Icges(6,5) * t202 - Icges(6,6) * t173 + Icges(6,3) * t172;
t125 = Icges(6,4) * t202 - Icges(6,2) * t173 + Icges(6,6) * t172;
t128 = Icges(6,1) * t202 - Icges(6,4) * t173 + Icges(6,5) * t172;
t48 = t122 * t150 - t125 * t151 + t128 * t178;
t161 = Icges(4,5) * t204 + Icges(4,6) * t203 + Icges(4,3) * t289;
t163 = Icges(4,4) * t204 + Icges(4,2) * t203 + Icges(4,6) * t289;
t165 = Icges(4,1) * t204 + Icges(4,4) * t203 + Icges(4,5) * t289;
t67 = t161 * t289 + t163 * t203 + t165 * t204;
t287 = t229 * t233;
t162 = Icges(4,5) * t206 + Icges(4,6) * t205 + Icges(4,3) * t287;
t164 = Icges(4,4) * t206 + Icges(4,2) * t205 + Icges(4,6) * t287;
t166 = Icges(4,1) * t206 + Icges(4,4) * t205 + Icges(4,5) * t287;
t68 = t162 * t289 + t164 * t203 + t166 * t204;
t312 = (-t45 - t47 - t116 - t48) * t232 + ((t20 + t24 + t26 + t68) * t233 + (t19 + t23 + t25 + t67) * t230) * t229;
t117 = t187 * t287 + t205 * t190 + t206 * t193;
t21 = -t152 * t83 + t153 * t89 + t179 * t77;
t22 = -t152 * t84 + t153 * t90 + t179 * t78;
t27 = t152 * t79 + t153 * t73 + t179 * t85;
t28 = t152 * t80 + t153 * t74 + t179 * t86;
t29 = t152 * t75 - t153 * t81 + t179 * t87;
t30 = t152 * t76 - t153 * t82 + t179 * t88;
t46 = t123 * t179 - t126 * t152 + t129 * t153;
t49 = t121 * t153 + t124 * t152 + t127 * t179;
t50 = t122 * t152 - t125 * t153 + t128 * t179;
t69 = t161 * t287 + t205 * t163 + t206 * t165;
t70 = t162 * t287 + t205 * t164 + t206 * t166;
t311 = (-t46 - t117 - t49 - t50) * t232 + ((t22 + t28 + t30 + t70) * t233 + (t21 + t27 + t29 + t69) * t230) * t229;
t31 = -t172 * t83 + t173 * t89 + t202 * t77;
t33 = t172 * t79 + t173 * t73 + t202 * t85;
t35 = t172 * t75 - t173 * t81 + t202 * t87;
t71 = -t232 * t161 + (-t163 * t228 + t165 * t231) * t229;
t310 = -t31 - t33 - t35 - t71;
t32 = -t172 * t84 + t173 * t90 + t202 * t78;
t34 = t172 * t80 + t173 * t74 + t202 * t86;
t36 = t172 * t76 - t173 * t82 + t202 * t88;
t72 = -t232 * t162 + (-t164 * t228 + t166 * t231) * t229;
t309 = t32 + t34 + t36 + t72;
t308 = t193 * t288 + (t123 + t127 + t128) * t202 + (t121 - t125 + t129) * t173 + (t122 + t124 - t126) * t172;
t307 = t230 ^ 2;
t306 = t233 ^ 2;
t305 = t230 / 0.2e1;
t304 = -t232 / 0.2e1;
t302 = pkin(2) * t232;
t300 = t233 * rSges(3,3);
t299 = t314 * t151 - t316;
t298 = t152 * rSges(7,2) + t314 * t153 + t315 * t179;
t159 = t206 * pkin(3) + t179 * qJ(4);
t92 = t153 * rSges(5,1) - t152 * rSges(5,2) + t179 * rSges(5,3);
t297 = -t159 - t92;
t291 = Icges(3,4) * t232;
t141 = t150 * qJ(5);
t111 = t151 * pkin(4) + t141;
t269 = t204 * pkin(3) - qJ(4) * t268;
t158 = qJ(4) * t258 + t269;
t149 = t158 * t287;
t284 = t111 * t287 + t149;
t112 = t153 * pkin(4) + t152 * qJ(5);
t283 = -t112 - t159;
t130 = rSges(5,1) * t173 - rSges(5,2) * t172 + rSges(5,3) * t202;
t183 = pkin(3) * t288 + qJ(4) * t202;
t282 = -t130 - t183;
t281 = rSges(7,2) * t172 + t314 * t173 + t315 * t202;
t138 = pkin(4) * t173 + qJ(5) * t172;
t280 = -t138 - t183;
t279 = t232 * t158 + t183 * t289;
t196 = -t232 * rSges(4,3) + (rSges(4,1) * t231 - rSges(4,2) * t228) * t229;
t219 = pkin(2) * t229 - pkin(9) * t232;
t278 = -t196 - t219;
t276 = pkin(2) * t285 + pkin(9) * t287;
t277 = t307 * (pkin(9) * t229 + t302) + t233 * t276;
t275 = t233 * pkin(1) + t230 * pkin(8);
t96 = t179 * rSges(6,1) - t153 * rSges(6,2) + t152 * rSges(6,3);
t273 = -t96 + t283;
t272 = -t219 + t282;
t132 = rSges(6,1) * t202 - rSges(6,2) * t173 + rSges(6,3) * t172;
t271 = -t132 + t280;
t168 = t206 * rSges(4,1) + t205 * rSges(4,2) + rSges(4,3) * t287;
t270 = -pkin(1) - t302;
t265 = -t232 * t187 - t190 * t290 + t308;
t264 = t283 - t298;
t263 = t232 * t111 + t138 * t289 + t279;
t262 = t280 - t281;
t261 = -t219 + t271;
t260 = t230 * t158 + t233 * t159 + t277;
t259 = t275 + t276;
t252 = rSges(3,1) * t232 - rSges(3,2) * t229;
t251 = -t204 * rSges(4,1) - t203 * rSges(4,2);
t250 = -t178 * rSges(6,1) - t150 * rSges(6,3);
t249 = -t219 + t262;
t247 = -Icges(3,2) * t229 + t291;
t246 = Icges(3,5) * t232 - Icges(3,6) * t229;
t243 = rSges(3,1) * t285 - rSges(3,2) * t287 + t230 * rSges(3,3);
t240 = t230 * t111 + t233 * t112 + t260;
t91 = t151 * rSges(5,1) - t150 * rSges(5,2) + t178 * rSges(5,3);
t239 = t159 + t259;
t238 = t71 / 0.2e1 + t31 / 0.2e1 + t35 / 0.2e1 + t33 / 0.2e1 + t116 / 0.2e1 + t45 / 0.2e1 + t48 / 0.2e1 + t47 / 0.2e1;
t237 = t72 / 0.2e1 + t32 / 0.2e1 + t36 / 0.2e1 + t34 / 0.2e1 + t117 / 0.2e1 + t46 / 0.2e1 + t50 / 0.2e1 + t49 / 0.2e1;
t236 = t112 + t239;
t226 = t233 * pkin(8);
t235 = t226 + ((-qJ(4) * t296 - pkin(9)) * t229 + t270) * t230 - t269;
t234 = -t141 + t235;
t217 = rSges(2,1) * t233 - t230 * rSges(2,2);
t216 = -t230 * rSges(2,1) - rSges(2,2) * t233;
t215 = rSges(3,1) * t229 + rSges(3,2) * t232;
t212 = Icges(3,6) * t232 + t318;
t189 = Icges(3,3) * t230 + t233 * t246;
t188 = -Icges(3,3) * t233 + t230 * t246;
t181 = t243 + t275;
t180 = t300 + t226 + (-pkin(1) - t252) * t230;
t170 = t278 * t233;
t169 = t278 * t230;
t167 = rSges(4,3) * t289 - t251;
t160 = t233 * t243 + (t230 * t252 - t300) * t230;
t140 = t259 + t168;
t139 = t226 + ((-rSges(4,3) - pkin(9)) * t229 + t270) * t230 + t251;
t137 = -t232 * t168 - t196 * t287;
t136 = t167 * t232 + t196 * t289;
t118 = (t167 * t233 - t168 * t230) * t229;
t100 = t230 * t167 + t168 * t233 + t277;
t99 = t272 * t233;
t98 = t272 * t230;
t94 = -t151 * rSges(6,2) - t250;
t66 = t261 * t233;
t65 = t261 * t230;
t64 = t239 + t92;
t63 = t235 - t91;
t61 = t249 * t233;
t60 = t249 * t230;
t59 = t232 * t297 + t282 * t287;
t58 = t130 * t289 + t232 * t91 + t279;
t57 = t236 + t96;
t56 = (rSges(6,2) - pkin(4)) * t151 + t234 + t250;
t52 = t70 * t230 - t233 * t69;
t51 = t68 * t230 - t233 * t67;
t44 = t149 + (t230 * t297 + t233 * t91) * t229;
t43 = t236 + t298;
t42 = (-pkin(4) - t314) * t151 + t234 + t316;
t41 = t230 * t91 + t233 * t92 + t260;
t38 = t232 * t273 + t271 * t287;
t37 = t132 * t289 + t232 * t94 + t263;
t18 = (t230 * t273 + t233 * t94) * t229 + t284;
t17 = t232 * t264 + t262 * t287;
t16 = t232 * t299 + t281 * t289 + t263;
t15 = t230 * t94 + t233 * t96 + t240;
t14 = (t230 * t264 + t233 * t299) * t229 + t284;
t13 = t230 * t299 + t233 * t298 + t240;
t12 = t30 * t230 - t233 * t29;
t11 = t28 * t230 - t233 * t27;
t10 = t26 * t230 - t233 * t25;
t9 = -t23 * t233 + t24 * t230;
t8 = -t21 * t233 + t22 * t230;
t7 = -t19 * t233 + t20 * t230;
t1 = [Icges(2,3) + (Icges(3,4) * t229 + Icges(3,2) * t232 - t187) * t232 + (Icges(3,1) * t229 - t190 * t228 + t291) * t229 + m(6) * (t56 ^ 2 + t57 ^ 2) + m(7) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t63 ^ 2 + t64 ^ 2) + m(4) * (t139 ^ 2 + t140 ^ 2) + m(3) * (t180 ^ 2 + t181 ^ 2) + m(2) * (t216 ^ 2 + t217 ^ 2) + t308; (t230 * t247 * t304 - t238 + (-Icges(3,6) * t304 + t317 + t212 / 0.2e1) * t233) * t233 + (t232 * (Icges(3,6) * t230 + t233 * t247) / 0.2e1 + t230 * t317 + t212 * t305 + t237) * t230 + m(7) * (t42 * t61 + t43 * t60) + m(6) * (t56 * t66 + t57 * t65) + m(5) * (t63 * t99 + t64 * t98) + m(4) * (t139 * t170 + t140 * t169) + m(3) * (-t180 * t233 - t181 * t230) * t215; m(3) * (t160 ^ 2 + (t306 + t307) * t215 ^ 2) + m(7) * (t13 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t15 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t41 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(4) * (t100 ^ 2 + t169 ^ 2 + t170 ^ 2) + (-t306 * t188 - t10 - t51 - t7 - t9) * t233 + (t307 * t189 + t11 + t12 + t52 + t8 + (-t230 * t188 + t233 * t189) * t233) * t230; -t265 * t232 + m(6) * (t37 * t56 + t38 * t57) + m(7) * (t16 * t42 + t17 * t43) + m(5) * (t58 * t63 + t59 * t64) + m(4) * (t136 * t139 + t137 * t140) + (t230 * t238 + t233 * t237) * t229; m(7) * (t13 * t14 + t16 * t61 + t17 * t60) + m(6) * (t15 * t18 + t37 * t66 + t38 * t65) + m(5) * (t44 * t41 + t58 * t99 + t59 * t98) + m(4) * (t100 * t118 + t136 * t170 + t137 * t169) + ((t11 / 0.2e1 + t52 / 0.2e1 + t8 / 0.2e1 + t12 / 0.2e1) * t233 + (t10 / 0.2e1 + t9 / 0.2e1 + t7 / 0.2e1 + t51 / 0.2e1) * t230) * t229 + t311 * t305 + (t309 * t230 + t310 * t233) * t304 - t312 * t233 / 0.2e1; m(6) * (t18 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(4) * (t118 ^ 2 + t136 ^ 2 + t137 ^ 2) + m(7) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(5) * (t44 ^ 2 + t58 ^ 2 + t59 ^ 2) + t265 * t232 ^ 2 + ((-t309 * t232 + t311) * t233 + (t310 * t232 + t312) * t230) * t229; m(6) * (t178 * t57 + t179 * t56) + m(7) * (t178 * t43 + t179 * t42) + m(5) * (t178 * t64 + t179 * t63); m(7) * (t13 * t202 + t178 * t60 + t179 * t61) + m(6) * (t15 * t202 + t178 * t65 + t179 * t66) + m(5) * (t178 * t98 + t179 * t99 + t202 * t41); m(7) * (t14 * t202 + t16 * t179 + t17 * t178) + m(5) * (t178 * t59 + t179 * t58 + t202 * t44) + m(6) * (t178 * t38 + t179 * t37 + t18 * t202); 0.2e1 * (m(5) / 0.2e1 + t274) * (t178 ^ 2 + t179 ^ 2 + t202 ^ 2); m(6) * (t150 * t57 + t152 * t56) + m(7) * (t150 * t43 + t152 * t42); m(7) * (t13 * t172 + t150 * t60 + t152 * t61) + m(6) * (t15 * t172 + t150 * t65 + t152 * t66); m(7) * (t14 * t172 + t150 * t17 + t152 * t16) + m(6) * (t150 * t38 + t152 * t37 + t172 * t18); (t150 * t178 + t152 * t179 + t172 * t202) * t313; (t150 ^ 2 + t152 ^ 2 + t172 ^ 2) * t313; m(7) * (t151 * t43 + t153 * t42); m(7) * (t13 * t173 + t151 * t60 + t153 * t61); m(7) * (t14 * t173 + t151 * t17 + t153 * t16); m(7) * (t151 * t178 + t153 * t179 + t173 * t202); m(7) * (t150 * t151 + t152 * t153 + t172 * t173); m(7) * (t151 ^ 2 + t153 ^ 2 + t173 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
