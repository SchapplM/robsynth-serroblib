% Calculate joint inertia matrix for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP14_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP14_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:53
% EndTime: 2019-03-09 13:04:06
% DurationCPUTime: 5.54s
% Computational Cost: add. (15652->605), mult. (40123->833), div. (0->0), fcn. (50962->10), ass. (0->265)
t259 = cos(pkin(6));
t264 = cos(qJ(2));
t265 = cos(qJ(1));
t311 = t264 * t265;
t262 = sin(qJ(2));
t263 = sin(qJ(1));
t314 = t262 * t263;
t240 = -t259 * t311 + t314;
t261 = sin(qJ(4));
t258 = sin(pkin(6));
t322 = cos(qJ(4));
t286 = t258 * t322;
t212 = t240 * t261 - t265 * t286;
t312 = t263 * t264;
t313 = t262 * t265;
t241 = t259 * t313 + t312;
t260 = sin(qJ(5));
t321 = cos(qJ(5));
t168 = t212 * t260 - t241 * t321;
t169 = t212 * t321 + t241 * t260;
t315 = t258 * t265;
t211 = t240 * t322 + t261 * t315;
t329 = rSges(7,3) + qJ(6);
t330 = rSges(7,1) + pkin(5);
t309 = -rSges(7,2) * t211 + t329 * t168 + t169 * t330;
t219 = Icges(3,3) * t259 + (Icges(3,5) * t262 + Icges(3,6) * t264) * t258;
t220 = Icges(3,6) * t259 + (Icges(3,4) * t262 + Icges(3,2) * t264) * t258;
t221 = Icges(3,5) * t259 + (Icges(3,1) * t262 + Icges(3,4) * t264) * t258;
t222 = Icges(4,5) * t259 + (-Icges(4,6) * t262 - Icges(4,3) * t264) * t258;
t223 = Icges(4,4) * t259 + (-Icges(4,2) * t262 - Icges(4,6) * t264) * t258;
t224 = Icges(4,1) * t259 + (-Icges(4,4) * t262 - Icges(4,5) * t264) * t258;
t316 = t258 * t264;
t318 = t258 * t262;
t328 = (-t222 * t264 - t223 * t262) * t258 + t220 * t316 + t221 * t318 + (t224 + t219) * t259;
t176 = -Icges(4,5) * t315 - Icges(4,6) * t241 + Icges(4,3) * t240;
t183 = Icges(3,4) * t241 - Icges(3,2) * t240 - Icges(3,6) * t315;
t327 = t176 - t183;
t178 = -Icges(4,4) * t315 - Icges(4,2) * t241 + Icges(4,6) * t240;
t185 = Icges(3,1) * t241 - Icges(3,4) * t240 - Icges(3,5) * t315;
t326 = t178 - t185;
t242 = t259 * t312 + t313;
t243 = -t259 * t314 + t311;
t317 = t258 * t263;
t175 = Icges(4,5) * t317 - Icges(4,6) * t243 + Icges(4,3) * t242;
t184 = Icges(3,4) * t243 - Icges(3,2) * t242 + Icges(3,6) * t317;
t325 = -t184 + t175;
t177 = Icges(4,4) * t317 - Icges(4,2) * t243 + Icges(4,6) * t242;
t186 = Icges(3,1) * t243 - Icges(3,4) * t242 + Icges(3,5) * t317;
t324 = t186 - t177;
t323 = t258 ^ 2;
t320 = t241 * pkin(2);
t210 = t242 * t261 + t263 * t286;
t166 = t210 * t260 - t243 * t321;
t167 = t210 * t321 + t243 * t260;
t209 = -t242 * t322 + t261 * t317;
t310 = t209 * rSges(7,2) + t329 * t166 + t167 * t330;
t308 = t328 * t259;
t239 = t259 * t322 - t261 * t316;
t206 = t239 * t260 - t318 * t321;
t207 = t239 * t321 + t260 * t318;
t238 = t259 * t261 + t264 * t286;
t307 = rSges(7,2) * t238 + t329 * t206 + t207 * t330;
t180 = -Icges(4,1) * t315 - Icges(4,4) * t241 + Icges(4,5) * t240;
t181 = Icges(3,5) * t241 - Icges(3,6) * t240 - Icges(3,3) * t315;
t306 = -t181 - t180;
t179 = Icges(4,1) * t317 - Icges(4,4) * t243 + Icges(4,5) * t242;
t182 = Icges(3,5) * t243 - Icges(3,6) * t242 + Icges(3,3) * t317;
t305 = t182 + t179;
t230 = t240 * qJ(3);
t196 = t230 + t320;
t197 = t243 * pkin(2) + qJ(3) * t242;
t304 = t196 * t317 + t197 * t315;
t194 = t259 * t197;
t215 = pkin(3) * t317 + pkin(9) * t243;
t303 = t259 * t215 + t194;
t216 = -pkin(3) * t315 + t241 * pkin(9);
t302 = -t196 - t216;
t244 = (pkin(2) * t262 - qJ(3) * t264) * t258;
t301 = -pkin(3) * t259 - pkin(9) * t318 - t244;
t300 = t265 * pkin(1) + pkin(8) * t317;
t103 = Icges(7,1) * t167 + Icges(7,4) * t209 + Icges(7,5) * t166;
t95 = Icges(7,5) * t167 + Icges(7,6) * t209 + Icges(7,3) * t166;
t99 = Icges(7,4) * t167 + Icges(7,2) * t209 + Icges(7,6) * t166;
t33 = t103 * t167 + t166 * t95 + t209 * t99;
t100 = Icges(7,4) * t169 - Icges(7,2) * t211 + Icges(7,6) * t168;
t104 = Icges(7,1) * t169 - Icges(7,4) * t211 + Icges(7,5) * t168;
t96 = Icges(7,5) * t169 - Icges(7,6) * t211 + Icges(7,3) * t168;
t34 = t100 * t209 + t104 * t167 + t166 * t96;
t130 = Icges(7,5) * t207 + Icges(7,6) * t238 + Icges(7,3) * t206;
t132 = Icges(7,4) * t207 + Icges(7,2) * t238 + Icges(7,6) * t206;
t134 = Icges(7,1) * t207 + Icges(7,4) * t238 + Icges(7,5) * t206;
t50 = t130 * t166 + t132 * t209 + t134 * t167;
t1 = t209 * t33 - t211 * t34 + t238 * t50;
t101 = Icges(6,4) * t167 - Icges(6,2) * t166 + Icges(6,6) * t209;
t105 = Icges(6,1) * t167 - Icges(6,4) * t166 + Icges(6,5) * t209;
t97 = Icges(6,5) * t167 - Icges(6,6) * t166 + Icges(6,3) * t209;
t35 = -t101 * t166 + t105 * t167 + t209 * t97;
t102 = Icges(6,4) * t169 - Icges(6,2) * t168 - Icges(6,6) * t211;
t106 = Icges(6,1) * t169 - Icges(6,4) * t168 - Icges(6,5) * t211;
t98 = Icges(6,5) * t169 - Icges(6,6) * t168 - Icges(6,3) * t211;
t36 = -t102 * t166 + t106 * t167 + t209 * t98;
t131 = Icges(6,5) * t207 - Icges(6,6) * t206 + Icges(6,3) * t238;
t133 = Icges(6,4) * t207 - Icges(6,2) * t206 + Icges(6,6) * t238;
t135 = Icges(6,1) * t207 - Icges(6,4) * t206 + Icges(6,5) * t238;
t51 = t131 * t209 - t133 * t166 + t135 * t167;
t2 = t209 * t35 - t211 * t36 + t238 * t51;
t299 = t2 / 0.2e1 + t1 / 0.2e1;
t37 = t103 * t169 + t168 * t95 - t211 * t99;
t38 = -t100 * t211 + t104 * t169 + t168 * t96;
t52 = t130 * t168 - t132 * t211 + t134 * t169;
t3 = t209 * t37 - t211 * t38 + t238 * t52;
t39 = -t101 * t168 + t105 * t169 - t211 * t97;
t40 = -t102 * t168 + t106 * t169 - t211 * t98;
t53 = -t131 * t211 - t133 * t168 + t135 * t169;
t4 = t209 * t39 - t211 * t40 + t238 * t53;
t298 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t241 * t34 + t243 * t33 + t318 * t50;
t6 = t241 * t36 + t243 * t35 + t318 * t51;
t297 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t241 * t38 + t243 * t37 + t318 * t52;
t8 = t241 * t40 + t243 * t39 + t318 * t53;
t296 = -t7 / 0.2e1 - t8 / 0.2e1;
t10 = t51 * t259 + (t263 * t35 - t265 * t36) * t258;
t9 = t50 * t259 + (t263 * t33 - t265 * t34) * t258;
t295 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t52 * t259 + (t263 * t37 - t265 * t38) * t258;
t12 = t53 * t259 + (t263 * t39 - t265 * t40) * t258;
t294 = -t11 / 0.2e1 - t12 / 0.2e1;
t41 = t103 * t207 + t206 * t95 + t238 * t99;
t42 = t100 * t238 + t104 * t207 + t206 * t96;
t65 = t206 * t130 + t238 * t132 + t207 * t134;
t59 = t65 * t238;
t13 = t41 * t209 - t42 * t211 + t59;
t43 = -t101 * t206 + t105 * t207 + t238 * t97;
t44 = -t102 * t206 + t106 * t207 + t238 * t98;
t66 = t238 * t131 - t206 * t133 + t207 * t135;
t60 = t66 * t238;
t14 = t43 * t209 - t44 * t211 + t60;
t293 = t13 / 0.2e1 + t14 / 0.2e1;
t61 = t65 * t318;
t15 = t42 * t241 + t41 * t243 + t61;
t62 = t66 * t318;
t16 = t44 * t241 + t43 * t243 + t62;
t292 = t15 / 0.2e1 + t16 / 0.2e1;
t63 = t65 * t259;
t17 = t63 + (t41 * t263 - t42 * t265) * t258;
t64 = t66 * t259;
t18 = t64 + (t43 * t263 - t44 * t265) * t258;
t291 = t18 / 0.2e1 + t17 / 0.2e1;
t154 = t210 * pkin(4) + pkin(10) * t209;
t290 = t259 * t154 + t303;
t155 = pkin(4) * t212 - t211 * pkin(10);
t289 = -t155 + t302;
t171 = Icges(5,5) * t239 - Icges(5,6) * t238 + Icges(5,3) * t318;
t172 = Icges(5,4) * t239 - Icges(5,2) * t238 + Icges(5,6) * t318;
t173 = Icges(5,1) * t239 - Icges(5,4) * t238 + Icges(5,5) * t318;
t86 = t171 * t318 - t238 * t172 + t239 * t173;
t108 = t167 * rSges(6,1) - t166 * rSges(6,2) + t209 * rSges(6,3);
t195 = pkin(4) * t239 + pkin(10) * t238;
t288 = -t195 + t301;
t144 = t210 * rSges(5,1) - t209 * rSges(5,2) + t243 * rSges(5,3);
t190 = t243 * rSges(3,1) - t242 * rSges(3,2) + rSges(3,3) * t317;
t187 = rSges(4,1) * t317 - t243 * rSges(4,2) + t242 * rSges(4,3);
t285 = -t263 * pkin(1) + pkin(8) * t315;
t284 = t258 * (-rSges(4,1) * t259 - (-rSges(4,2) * t262 - rSges(4,3) * t264) * t258 - t244);
t283 = t215 * t315 + t216 * t317 + t304;
t282 = -t230 + t285;
t174 = rSges(5,1) * t239 - rSges(5,2) * t238 + rSges(5,3) * t318;
t281 = t258 * (-t174 + t301);
t280 = -rSges(5,1) * t212 - rSges(5,2) * t211;
t137 = rSges(6,1) * t207 - rSges(6,2) * t206 + rSges(6,3) * t238;
t279 = t258 * (-t137 + t288);
t278 = t197 + t300;
t277 = -t42 / 0.2e1 - t44 / 0.2e1 - t53 / 0.2e1 - t52 / 0.2e1;
t276 = t43 / 0.2e1 + t50 / 0.2e1 + t51 / 0.2e1 + t41 / 0.2e1;
t275 = rSges(4,1) * t315 - t240 * rSges(4,3);
t274 = t154 * t315 + t155 * t317 + t283;
t273 = t282 - t216;
t272 = t258 * (t288 - t307);
t110 = rSges(6,1) * t169 - rSges(6,2) * t168 - rSges(6,3) * t211;
t189 = t241 * rSges(3,1) - t240 * rSges(3,2) - rSges(3,3) * t315;
t139 = Icges(5,5) * t212 + Icges(5,6) * t211 + Icges(5,3) * t241;
t141 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t241;
t143 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t241;
t73 = t139 * t318 - t141 * t238 + t143 * t239;
t81 = t171 * t241 + t172 * t211 + t173 * t212;
t270 = t73 / 0.2e1 + t81 / 0.2e1 - t277;
t138 = Icges(5,5) * t210 - Icges(5,6) * t209 + Icges(5,3) * t243;
t140 = Icges(5,4) * t210 - Icges(5,2) * t209 + Icges(5,6) * t243;
t142 = Icges(5,1) * t210 - Icges(5,4) * t209 + Icges(5,5) * t243;
t72 = t138 * t318 - t140 * t238 + t142 * t239;
t80 = t171 * t243 - t172 * t209 + t173 * t210;
t269 = t72 / 0.2e1 + t80 / 0.2e1 + t276;
t268 = t215 + t278;
t267 = -t155 + t273 - t320;
t266 = t154 + t268;
t247 = rSges(2,1) * t265 - t263 * rSges(2,2);
t246 = -t263 * rSges(2,1) - rSges(2,2) * t265;
t225 = rSges(3,3) * t259 + (rSges(3,1) * t262 + rSges(3,2) * t264) * t258;
t188 = -t241 * rSges(4,2) - t275;
t165 = t241 * t195;
t159 = t190 + t300;
t158 = -t189 + t285;
t150 = t154 * t318;
t148 = -t259 * t189 - t225 * t315;
t147 = t190 * t259 - t225 * t317;
t146 = t243 * t155;
t145 = rSges(5,3) * t241 - t280;
t127 = t278 + t187;
t126 = (rSges(4,2) - pkin(2)) * t241 + t275 + t282;
t123 = (t189 * t263 + t190 * t265) * t258;
t122 = t219 * t317 - t220 * t242 + t221 * t243;
t121 = -t219 * t315 - t240 * t220 + t241 * t221;
t120 = t240 * t222 - t241 * t223 - t224 * t315;
t119 = t222 * t242 - t223 * t243 + t224 * t317;
t112 = (-t188 - t196) * t259 + t265 * t284;
t111 = t187 * t259 + t263 * t284 + t194;
t94 = t144 * t318 - t174 * t243;
t93 = -t145 * t318 + t174 * t241;
t92 = t180 * t259 + (-t176 * t264 - t178 * t262) * t258;
t91 = t179 * t259 + (-t175 * t264 - t177 * t262) * t258;
t90 = t182 * t259 + (t184 * t264 + t186 * t262) * t258;
t89 = t181 * t259 + (t183 * t264 + t185 * t262) * t258;
t88 = t268 + t144;
t87 = (-rSges(5,3) - pkin(2)) * t241 + t273 + t280;
t85 = t86 * t259;
t84 = t86 * t318;
t83 = (t187 * t265 + t188 * t263) * t258 + t304;
t82 = -t144 * t241 + t145 * t243;
t79 = (-t145 + t302) * t259 + t265 * t281;
t78 = t144 * t259 + t263 * t281 + t303;
t77 = -t110 * t238 - t137 * t211;
t76 = t108 * t238 - t137 * t209;
t75 = t266 + t108;
t74 = -t110 + t267;
t71 = t139 * t241 + t141 * t211 + t143 * t212;
t70 = t138 * t241 + t140 * t211 + t142 * t212;
t69 = t139 * t243 - t141 * t209 + t143 * t210;
t68 = t138 * t243 - t140 * t209 + t142 * t210;
t67 = (t144 * t265 + t145 * t263) * t258 + t283;
t58 = t108 * t211 + t110 * t209;
t57 = t108 * t318 + t150 + (-t137 - t195) * t243;
t56 = t137 * t241 + t165 + (-t110 - t155) * t318;
t55 = t266 + t310;
t54 = t267 - t309;
t49 = (-t110 + t289) * t259 + t265 * t279;
t48 = t108 * t259 + t263 * t279 + t290;
t47 = t110 * t243 + t146 + (-t108 - t154) * t241;
t46 = -t211 * t307 - t238 * t309;
t45 = -t209 * t307 + t238 * t310;
t32 = t150 + t310 * t318 + (-t195 - t307) * t243;
t31 = t165 + t307 * t241 + (-t155 - t309) * t318;
t30 = (t108 * t265 + t110 * t263) * t258 + t274;
t29 = t209 * t309 + t211 * t310;
t28 = (t289 - t309) * t259 + t265 * t272;
t27 = t259 * t310 + t263 * t272 + t290;
t26 = t146 + t309 * t243 + (-t154 - t310) * t241;
t25 = t85 + (t72 * t263 - t73 * t265) * t258;
t24 = t73 * t241 + t72 * t243 + t84;
t23 = (t263 * t309 + t265 * t310) * t258 + t274;
t22 = t81 * t259 + (t263 * t70 - t265 * t71) * t258;
t21 = t80 * t259 + (t263 * t68 - t265 * t69) * t258;
t20 = t241 * t71 + t243 * t70 + t318 * t81;
t19 = t241 * t69 + t243 * t68 + t318 * t80;
t107 = [Icges(2,3) + m(7) * (t54 ^ 2 + t55 ^ 2) + m(6) * (t74 ^ 2 + t75 ^ 2) + m(5) * (t87 ^ 2 + t88 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t158 ^ 2 + t159 ^ 2) + m(2) * (t246 ^ 2 + t247 ^ 2) + t86 + t66 + t65 + t328; t64 + t63 + t85 + m(7) * (t27 * t55 + t28 * t54) + m(6) * (t48 * t75 + t49 * t74) + m(5) * (t78 * t88 + t79 * t87) + m(4) * (t111 * t127 + t112 * t126) + m(3) * (t147 * t159 + t148 * t158) + ((-t92 / 0.2e1 - t89 / 0.2e1 - t120 / 0.2e1 - t121 / 0.2e1 - t270) * t265 + (t91 / 0.2e1 + t90 / 0.2e1 + t119 / 0.2e1 + t122 / 0.2e1 + t269) * t263) * t258 + t308; (t18 + t17 + t25 + t308) * t259 + m(7) * (t23 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t30 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t67 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2 + t83 ^ 2) + m(3) * (t123 ^ 2 + t147 ^ 2 + t148 ^ 2) + ((-t11 - t12 - t22 + ((t327 * t240 - t326 * t241) * t258 + t306 * t323 * t265) * t265 + (-t120 - t121 - t89 - t92) * t259) * t265 + (t10 + t9 + t21 + ((t325 * t242 + t324 * t243) * t258 + t305 * t323 * t263) * t263 + (t122 + t119 + t91 + t90) * t259 + ((t263 * t306 + t265 * t305) * t258 + t326 * t243 - t327 * t242 - t324 * t241 - t325 * t240) * t315) * t263) * t258; m(7) * (t240 * t55 + t242 * t54) + m(6) * (t240 * t75 + t242 * t74) + m(5) * (t240 * t88 + t242 * t87) + m(4) * (t126 * t242 + t127 * t240); m(7) * (-t23 * t316 + t240 * t27 + t242 * t28) + m(6) * (t240 * t48 + t242 * t49 - t30 * t316) + m(5) * (t240 * t78 + t242 * t79 - t316 * t67) + m(4) * (t111 * t240 + t112 * t242 - t316 * t83); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t264 ^ 2 * t323 + t240 ^ 2 + t242 ^ 2); t62 + t61 + t84 + m(7) * (t31 * t54 + t32 * t55) + m(6) * (t56 * t74 + t57 * t75) + m(5) * (t87 * t93 + t88 * t94) + t269 * t243 + t270 * t241; (t24 / 0.2e1 + t292) * t259 + (t21 / 0.2e1 + t295) * t243 + (t22 / 0.2e1 - t294) * t241 + m(7) * (t23 * t26 + t27 * t32 + t28 * t31) + m(6) * (t30 * t47 + t48 * t57 + t49 * t56) + m(5) * (t67 * t82 + t78 * t94 + t79 * t93) + ((-t20 / 0.2e1 + t296) * t265 + (t19 / 0.2e1 + t297) * t263 + (t25 / 0.2e1 + t291) * t262) * t258; m(5) * (t240 * t94 + t242 * t93 - t316 * t82) + m(6) * (t240 * t57 + t242 * t56 - t316 * t47) + m(7) * (t240 * t32 + t242 * t31 - t26 * t316); (t15 + t16 + t24) * t318 + (t6 + t5 + t19) * t243 + (t8 + t7 + t20) * t241 + m(7) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t47 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t82 ^ 2 + t93 ^ 2 + t94 ^ 2); t60 + t59 + m(7) * (t45 * t55 + t46 * t54) + m(6) * (t74 * t77 + t75 * t76) + t277 * t211 + t276 * t209; t293 * t259 + t291 * t238 + t294 * t211 + t295 * t209 + m(7) * (t23 * t29 + t27 * t45 + t28 * t46) + m(6) * (t58 * t30 + t48 * t76 + t49 * t77) + (t263 * t299 - t265 * t298) * t258; m(6) * (t240 * t76 + t242 * t77 - t316 * t58) + m(7) * (t240 * t45 + t242 * t46 - t29 * t316); t293 * t318 + t299 * t243 + t298 * t241 + t292 * t238 + t296 * t211 + t297 * t209 + m(7) * (t26 * t29 + t31 * t46 + t32 * t45) + m(6) * (t58 * t47 + t56 * t77 + t57 * t76); (t13 + t14) * t238 + (-t4 - t3) * t211 + (t1 + t2) * t209 + m(7) * (t29 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t58 ^ 2 + t76 ^ 2 + t77 ^ 2); m(7) * (t166 * t54 + t168 * t55); m(7) * (t166 * t28 + t168 * t27 + t206 * t23); m(7) * (t166 * t242 + t168 * t240 - t206 * t316); m(7) * (t166 * t31 + t168 * t32 + t206 * t26); m(7) * (t166 * t46 + t168 * t45 + t206 * t29); m(7) * (t166 ^ 2 + t168 ^ 2 + t206 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t107(1) t107(2) t107(4) t107(7) t107(11) t107(16); t107(2) t107(3) t107(5) t107(8) t107(12) t107(17); t107(4) t107(5) t107(6) t107(9) t107(13) t107(18); t107(7) t107(8) t107(9) t107(10) t107(14) t107(19); t107(11) t107(12) t107(13) t107(14) t107(15) t107(20); t107(16) t107(17) t107(18) t107(19) t107(20) t107(21);];
Mq  = res;
