% Calculate joint inertia matrix for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:21
% EndTime: 2019-03-09 04:01:31
% DurationCPUTime: 7.19s
% Computational Cost: add. (37837->540), mult. (102653->767), div. (0->0), fcn. (136255->16), ass. (0->251)
t283 = sin(pkin(13));
t284 = sin(pkin(7));
t251 = t284 * t283;
t285 = cos(pkin(13));
t252 = t285 * t284;
t293 = sin(qJ(3));
t295 = cos(qJ(3));
t207 = t251 * t295 + t252 * t293;
t231 = cos(pkin(7));
t241 = t283 * t295 + t285 * t293;
t209 = t241 * t231;
t230 = cos(pkin(12));
t228 = sin(pkin(12));
t235 = sin(qJ(1));
t279 = t235 * t228;
t232 = cos(pkin(6));
t237 = cos(qJ(1));
t280 = t232 * t237;
t214 = t230 * t280 - t279;
t278 = t235 * t230;
t215 = t228 * t280 + t278;
t219 = -t283 * t293 + t285 * t295;
t229 = sin(pkin(6));
t281 = t229 * t237;
t173 = -t207 * t281 + t214 * t209 + t215 * t219;
t234 = sin(qJ(5));
t294 = cos(qJ(5));
t309 = t214 * t284 + t231 * t281;
t156 = t173 * t234 + t294 * t309;
t216 = -t228 * t237 - t232 * t278;
t217 = t230 * t237 - t232 * t279;
t282 = t229 * t235;
t175 = t207 * t282 + t209 * t216 + t217 * t219;
t244 = t216 * t284 - t231 * t282;
t158 = t175 * t234 + t244 * t294;
t157 = t173 * t294 - t234 * t309;
t208 = t219 * t231;
t239 = -t251 * t293 + t252 * t295;
t238 = t229 * t239;
t172 = t208 * t214 - t215 * t241 - t237 * t238;
t233 = sin(qJ(6));
t236 = cos(qJ(6));
t112 = -t157 * t233 - t172 * t236;
t113 = t157 * t236 - t172 * t233;
t68 = Icges(7,5) * t113 + Icges(7,6) * t112 + Icges(7,3) * t156;
t70 = Icges(7,4) * t113 + Icges(7,2) * t112 + Icges(7,6) * t156;
t72 = Icges(7,1) * t113 + Icges(7,4) * t112 + Icges(7,5) * t156;
t17 = t112 * t70 + t113 * t72 + t156 * t68;
t183 = t232 * t207 + (t209 * t230 + t219 * t228) * t229;
t213 = -t229 * t230 * t284 + t232 * t231;
t170 = t183 * t234 - t213 * t294;
t159 = t175 * t294 - t234 * t244;
t174 = t216 * t208 - t217 * t241 + t235 * t238;
t114 = -t159 * t233 - t174 * t236;
t115 = t159 * t236 - t174 * t233;
t69 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t158;
t71 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t158;
t73 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t158;
t18 = t112 * t71 + t113 * t73 + t156 * t69;
t171 = t183 * t294 + t213 * t234;
t182 = t232 * t239 + (t208 * t230 - t228 * t241) * t229;
t133 = -t171 * t233 - t182 * t236;
t134 = t171 * t236 - t182 * t233;
t79 = Icges(7,5) * t134 + Icges(7,6) * t133 + Icges(7,3) * t170;
t80 = Icges(7,4) * t134 + Icges(7,2) * t133 + Icges(7,6) * t170;
t81 = Icges(7,1) * t134 + Icges(7,4) * t133 + Icges(7,5) * t170;
t26 = t112 * t80 + t113 * t81 + t156 * t79;
t1 = t156 * t17 + t158 * t18 + t170 * t26;
t310 = -t1 / 0.2e1;
t135 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t213;
t136 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t213;
t137 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t213;
t254 = t295 * t284;
t261 = t231 * t295;
t197 = t232 * t254 + (-t228 * t293 + t230 * t261) * t229;
t253 = t284 * t293;
t260 = t231 * t293;
t198 = t232 * t253 + (t228 * t295 + t230 * t260) * t229;
t176 = Icges(4,5) * t198 + Icges(4,6) * t197 + Icges(4,3) * t213;
t177 = Icges(4,4) * t198 + Icges(4,2) * t197 + Icges(4,6) * t213;
t178 = Icges(4,1) * t198 + Icges(4,4) * t197 + Icges(4,5) * t213;
t308 = t182 * t136 + t183 * t137 + t197 * t177 + t198 * t178 + (t135 + t176) * t213;
t299 = m(7) / 0.2e1;
t300 = m(6) / 0.2e1;
t301 = m(5) / 0.2e1;
t255 = t301 + t300 + t299;
t307 = 0.2e1 * t255;
t19 = t114 * t70 + t115 * t72 + t158 * t68;
t20 = t114 * t71 + t115 * t73 + t158 * t69;
t27 = t114 * t80 + t115 * t81 + t158 * t79;
t2 = t156 * t19 + t158 * t20 + t170 * t27;
t306 = -t2 / 0.2e1;
t305 = -t217 * pkin(2) + pkin(9) * t244;
t304 = t232 ^ 2;
t303 = m(3) / 0.2e1;
t302 = m(4) / 0.2e1;
t298 = t156 / 0.2e1;
t297 = t158 / 0.2e1;
t296 = t170 / 0.2e1;
t292 = pkin(5) * t157;
t225 = pkin(3) * t295 + pkin(2);
t290 = -pkin(2) + t225;
t289 = t308 * t213;
t250 = -rSges(7,1) * t113 - rSges(7,2) * t112;
t74 = rSges(7,3) * t156 - t250;
t288 = pkin(11) * t156 + t292 + t74;
t75 = t115 * rSges(7,1) + t114 * rSges(7,2) + t158 * rSges(7,3);
t287 = t159 * pkin(5) + pkin(11) * t158 + t75;
t82 = rSges(7,1) * t134 + rSges(7,2) * t133 + rSges(7,3) * t170;
t286 = pkin(5) * t171 + pkin(11) * t170 + t82;
t129 = pkin(4) * t173 - t172 * pkin(10);
t210 = pkin(3) * t253 + (pkin(9) + qJ(4)) * t231;
t259 = t284 * pkin(9);
t211 = pkin(3) * t260 - qJ(4) * t284 - t259;
t256 = -t210 * t281 + t211 * t214;
t269 = t309 * pkin(9);
t154 = t215 * t290 + t256 + t269;
t139 = t244 * t154;
t277 = -t129 * t244 - t139;
t130 = t175 * pkin(4) - pkin(10) * t174;
t262 = t210 * t282 + t216 * t211 + t217 * t225;
t155 = t262 + t305;
t140 = t213 * t155;
t276 = t213 * t130 + t140;
t275 = -t129 - t154;
t274 = -t130 - t155;
t149 = pkin(4) * t183 - pkin(10) * t182;
t180 = (-t231 * pkin(9) + t210) * t232 + ((t259 + t211) * t230 + t290 * t228) * t229;
t162 = t309 * t180;
t273 = -t149 * t309 - t162;
t116 = Icges(5,5) * t173 + Icges(5,6) * t172 - Icges(5,3) * t309;
t248 = t229 * t254;
t188 = t214 * t261 - t215 * t293 - t237 * t248;
t247 = t229 * t253;
t189 = t214 * t260 + t215 * t295 - t237 * t247;
t141 = Icges(4,5) * t189 + Icges(4,6) * t188 - Icges(4,3) * t309;
t272 = t141 + t116;
t117 = Icges(5,5) * t175 + Icges(5,6) * t174 - Icges(5,3) * t244;
t190 = t216 * t261 - t217 * t293 + t235 * t248;
t191 = t216 * t260 + t217 * t295 + t235 * t247;
t142 = Icges(4,5) * t191 + Icges(4,6) * t190 - Icges(4,3) * t244;
t271 = t142 + t117;
t270 = -t149 - t180;
t268 = t237 * pkin(1) + qJ(2) * t282;
t33 = t133 * t80 + t134 * t81 + t170 * t79;
t104 = Icges(6,5) * t171 - Icges(6,6) * t170 - Icges(6,3) * t182;
t105 = Icges(6,4) * t171 - Icges(6,2) * t170 - Icges(6,6) * t182;
t106 = Icges(6,1) * t171 - Icges(6,4) * t170 - Icges(6,5) * t182;
t52 = -t182 * t104 - t170 * t105 + t171 * t106;
t23 = t133 * t70 + t134 * t72 + t170 * t68;
t267 = t26 / 0.2e1 + t23 / 0.2e1;
t24 = t133 * t71 + t134 * t73 + t170 * t69;
t266 = t27 / 0.2e1 + t24 / 0.2e1;
t93 = t159 * rSges(6,1) - t158 * rSges(6,2) - t174 * rSges(6,3);
t123 = t175 * rSges(5,1) + t174 * rSges(5,2) - rSges(5,3) * t244;
t148 = t191 * rSges(4,1) + t190 * rSges(4,2) - rSges(4,3) * t244;
t258 = -t235 * pkin(1) + qJ(2) * t281;
t249 = t262 + t268;
t86 = Icges(6,5) * t157 - Icges(6,6) * t156 - Icges(6,3) * t172;
t88 = Icges(6,4) * t157 - Icges(6,2) * t156 - Icges(6,6) * t172;
t90 = Icges(6,1) * t157 - Icges(6,4) * t156 - Icges(6,5) * t172;
t40 = -t170 * t88 + t171 * t90 - t182 * t86;
t46 = -t104 * t172 - t105 * t156 + t106 * t157;
t246 = -t40 / 0.2e1 - t46 / 0.2e1 - t267;
t87 = Icges(6,5) * t159 - Icges(6,6) * t158 - Icges(6,3) * t174;
t89 = Icges(6,4) * t159 - Icges(6,2) * t158 - Icges(6,6) * t174;
t91 = Icges(6,1) * t159 - Icges(6,4) * t158 - Icges(6,5) * t174;
t41 = -t170 * t89 + t171 * t91 - t182 * t87;
t47 = -t104 * t174 - t105 * t158 + t106 * t159;
t245 = -t47 / 0.2e1 - t41 / 0.2e1 - t266;
t147 = rSges(4,1) * t189 + rSges(4,2) * t188 - rSges(4,3) * t309;
t122 = rSges(5,1) * t173 + rSges(5,2) * t172 - rSges(5,3) * t309;
t92 = rSges(6,1) * t157 - rSges(6,2) * t156 - rSges(6,3) * t172;
t243 = t130 + t249;
t242 = -t215 * t225 - t256 + t258;
t240 = -t129 + t242;
t222 = rSges(2,1) * t237 - t235 * rSges(2,2);
t221 = -t235 * rSges(2,1) - rSges(2,2) * t237;
t187 = rSges(3,1) * t217 + rSges(3,2) * t216 + rSges(3,3) * t282 + t268;
t186 = -t215 * rSges(3,1) - t214 * rSges(3,2) + rSges(3,3) * t281 + t258;
t179 = rSges(4,1) * t198 + rSges(4,2) * t197 + rSges(4,3) * t213;
t146 = Icges(4,1) * t191 + Icges(4,4) * t190 - Icges(4,5) * t244;
t145 = Icges(4,1) * t189 + Icges(4,4) * t188 - Icges(4,5) * t309;
t144 = Icges(4,4) * t191 + Icges(4,2) * t190 - Icges(4,6) * t244;
t143 = Icges(4,4) * t189 + Icges(4,2) * t188 - Icges(4,6) * t309;
t138 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t213;
t127 = t148 + t268 - t305;
t126 = -pkin(2) * t215 - t147 + t258 + t269;
t121 = Icges(5,1) * t175 + Icges(5,4) * t174 - Icges(5,5) * t244;
t120 = Icges(5,1) * t173 + Icges(5,4) * t172 - Icges(5,5) * t309;
t119 = Icges(5,4) * t175 + Icges(5,2) * t174 - Icges(5,6) * t244;
t118 = Icges(5,4) * t173 + Icges(5,2) * t172 - Icges(5,6) * t309;
t107 = rSges(6,1) * t171 - rSges(6,2) * t170 - rSges(6,3) * t182;
t101 = t148 * t213 + t179 * t244;
t100 = -t147 * t213 - t179 * t309;
t99 = t249 + t123;
t98 = -t122 + t242;
t96 = -t147 * t244 + t148 * t309;
t84 = -t176 * t244 + t177 * t190 + t178 * t191;
t83 = -t176 * t309 + t177 * t188 + t178 * t189;
t78 = t142 * t213 + t144 * t197 + t146 * t198;
t77 = t141 * t213 + t143 * t197 + t145 * t198;
t65 = t243 + t93;
t64 = t240 - t92;
t62 = t123 * t213 + t140 - (-t138 - t180) * t244;
t61 = -t138 * t309 - t162 + (-t122 - t154) * t213;
t60 = -t135 * t244 + t136 * t174 + t137 * t175;
t59 = -t135 * t309 + t136 * t172 + t137 * t173;
t58 = t107 * t174 - t182 * t93;
t57 = -t107 * t172 + t182 * t92;
t56 = t117 * t213 + t119 * t182 + t121 * t183;
t55 = t116 * t213 + t118 * t182 + t120 * t183;
t54 = -t244 * t122 - t139 - (-t123 - t155) * t309;
t53 = t172 * t93 - t174 * t92;
t51 = t52 * t213;
t50 = t52 * t182;
t49 = t243 + t287;
t48 = -t292 + (-rSges(7,3) - pkin(11)) * t156 + t240 + t250;
t45 = -t158 * t82 + t170 * t75;
t44 = t156 * t82 - t170 * t74;
t43 = t213 * t93 - (-t107 + t270) * t244 + t276;
t42 = -t107 * t309 + (-t92 + t275) * t213 + t273;
t39 = -t158 * t89 + t159 * t91 - t174 * t87;
t38 = -t158 * t88 + t159 * t90 - t174 * t86;
t37 = -t156 * t89 + t157 * t91 - t172 * t87;
t36 = -t156 * t88 + t157 * t90 - t172 * t86;
t35 = -t156 * t75 + t158 * t74;
t34 = -t244 * t92 - (-t93 + t274) * t309 + t277;
t32 = t33 * t213;
t31 = t174 * t286 - t182 * t287;
t30 = -t172 * t286 + t182 * t288;
t29 = t33 * t182;
t28 = t33 * t170;
t25 = t172 * t287 - t174 * t288;
t22 = t287 * t213 - (t270 - t286) * t244 + t276;
t21 = -t286 * t309 + (t275 - t288) * t213 + t273;
t16 = -t288 * t244 - (t274 - t287) * t309 + t277;
t15 = -t244 * t41 - t309 * t40 + t51;
t14 = -t40 * t172 - t41 * t174 - t50;
t13 = t213 * t47 - t244 * t39 - t309 * t38;
t12 = t213 * t46 - t244 * t37 - t309 * t36;
t11 = -t172 * t38 - t174 * t39 - t182 * t47;
t10 = -t172 * t36 - t174 * t37 - t182 * t46;
t9 = -t23 * t309 - t24 * t244 + t32;
t8 = -t23 * t172 - t24 * t174 - t29;
t7 = t23 * t156 + t24 * t158 + t28;
t6 = -t19 * t309 - t20 * t244 + t213 * t27;
t5 = -t17 * t309 - t18 * t244 + t213 * t26;
t4 = -t172 * t19 - t174 * t20 - t182 * t27;
t3 = -t17 * t172 - t174 * t18 - t182 * t26;
t63 = [m(7) * (t48 ^ 2 + t49 ^ 2) + m(6) * (t64 ^ 2 + t65 ^ 2) + m(5) * (t98 ^ 2 + t99 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t186 ^ 2 + t187 ^ 2) + m(2) * (t221 ^ 2 + t222 ^ 2) + Icges(3,3) * t304 + Icges(2,3) + ((Icges(3,2) * t230 ^ 2 + (Icges(3,1) * t228 + 0.2e1 * Icges(3,4) * t230) * t228) * t229 + 0.2e1 * t232 * (Icges(3,5) * t228 + Icges(3,6) * t230)) * t229 + t52 + t33 + t308; 0.2e1 * ((t235 * t48 - t237 * t49) * t299 + (t235 * t64 - t237 * t65) * t300 + (t235 * t98 - t237 * t99) * t301 + (t126 * t235 - t127 * t237) * t302 + (t186 * t235 - t187 * t237) * t303) * t229; 0.2e1 * (t303 + t302 + t255) * (t304 + (t235 ^ 2 + t237 ^ 2) * t229 ^ 2); t32 + t51 + m(7) * (t21 * t48 + t22 * t49) + m(6) * (t42 * t64 + t43 * t65) + m(5) * (t61 * t98 + t62 * t99) + m(4) * (t100 * t126 + t101 * t127) - (t84 / 0.2e1 + t60 / 0.2e1 + t78 / 0.2e1 + t56 / 0.2e1 - t245) * t244 - (t77 / 0.2e1 + t55 / 0.2e1 + t83 / 0.2e1 + t59 / 0.2e1 - t246) * t309 + t289; m(4) * (t96 * t232 + (t100 * t235 - t101 * t237) * t229) + m(5) * (t54 * t232 + (t235 * t61 - t237 * t62) * t229) + m(6) * (t34 * t232 + (t235 * t42 - t237 * t43) * t229) + m(7) * (t16 * t232 + (t21 * t235 - t22 * t237) * t229); (t9 + t15 + t289) * t213 + m(7) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(6) * (t34 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t54 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(4) * (t100 ^ 2 + t101 ^ 2 + t96 ^ 2) - (t13 + t6 - (t174 * t119 + t175 * t121 + t190 * t144 + t191 * t146 - t244 * t271) * t244 + (t56 + t60 + t78 + t84) * t213) * t244 - (t5 + t12 - (t172 * t118 + t173 * t120 + t188 * t143 + t189 * t145 - t272 * t309) * t309 + (t77 + t55 + t83 + t59) * t213 - (t174 * t118 + t172 * t119 + t175 * t120 + t173 * t121 + t190 * t143 + t188 * t144 + t191 * t145 + t189 * t146 - t244 * t272 - t271 * t309) * t244) * t309; m(7) * (-t244 * t48 - t309 * t49) + m(6) * (-t244 * t64 - t309 * t65) + m(5) * (-t244 * t98 - t309 * t99); (t213 * t232 + (-t235 * t244 + t237 * t309) * t229) * t307; m(7) * (t16 * t213 - t21 * t244 - t22 * t309) + m(6) * (t213 * t34 - t244 * t42 - t309 * t43) + m(5) * (t213 * t54 - t244 * t61 - t309 * t62); (t213 ^ 2 + t244 ^ 2 + t309 ^ 2) * t307; -t29 - t50 + m(7) * (t30 * t48 + t31 * t49) + m(6) * (t57 * t64 + t58 * t65) + t245 * t174 + t246 * t172; m(6) * (t53 * t232 + (t235 * t57 - t237 * t58) * t229) + m(7) * (t25 * t232 + (t235 * t30 - t237 * t31) * t229); (t8 / 0.2e1 + t14 / 0.2e1) * t213 - (t4 / 0.2e1 + t11 / 0.2e1) * t244 - (t3 / 0.2e1 + t10 / 0.2e1) * t309 + (-t9 / 0.2e1 - t15 / 0.2e1) * t182 + (-t6 / 0.2e1 - t13 / 0.2e1) * t174 + (-t5 / 0.2e1 - t12 / 0.2e1) * t172 + m(7) * (t16 * t25 + t21 * t30 + t22 * t31) + m(6) * (t34 * t53 + t42 * t57 + t43 * t58); m(6) * (t213 * t53 - t244 * t57 - t309 * t58) + m(7) * (t213 * t25 - t244 * t30 - t309 * t31); (-t8 - t14) * t182 + (-t4 - t11) * t174 + (-t3 - t10) * t172 + m(7) * (t25 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t53 ^ 2 + t57 ^ 2 + t58 ^ 2); t28 + m(7) * (t44 * t48 + t45 * t49) + t266 * t158 + t267 * t156; m(7) * (t35 * t232 + (t235 * t44 - t237 * t45) * t229); m(7) * (t16 * t35 + t21 * t44 + t22 * t45) + t309 * t310 + t6 * t297 + t5 * t298 + t9 * t296 + t213 * t7 / 0.2e1 + t244 * t306; m(7) * (t213 * t35 - t244 * t44 - t309 * t45); m(7) * (t25 * t35 + t30 * t44 + t31 * t45) + t172 * t310 + t4 * t297 + t8 * t296 + t3 * t298 - t182 * t7 / 0.2e1 + t174 * t306; t158 * t2 + t156 * t1 + t170 * t7 + m(7) * (t35 ^ 2 + t44 ^ 2 + t45 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t63(1) t63(2) t63(4) t63(7) t63(11) t63(16); t63(2) t63(3) t63(5) t63(8) t63(12) t63(17); t63(4) t63(5) t63(6) t63(9) t63(13) t63(18); t63(7) t63(8) t63(9) t63(10) t63(14) t63(19); t63(11) t63(12) t63(13) t63(14) t63(15) t63(20); t63(16) t63(17) t63(18) t63(19) t63(20) t63(21);];
Mq  = res;
