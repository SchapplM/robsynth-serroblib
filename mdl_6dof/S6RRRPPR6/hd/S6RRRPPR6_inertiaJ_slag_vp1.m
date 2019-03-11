% Calculate joint inertia matrix for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:48
% EndTime: 2019-03-09 15:47:58
% DurationCPUTime: 4.84s
% Computational Cost: add. (19862->645), mult. (34483->888), div. (0->0), fcn. (42823->12), ass. (0->294)
t315 = m(6) / 0.2e1 + m(7) / 0.2e1;
t352 = 0.2e1 * t315;
t272 = sin(pkin(6));
t282 = cos(qJ(1));
t339 = t272 * t282;
t351 = t272 ^ 2;
t273 = cos(pkin(6));
t278 = sin(qJ(1));
t281 = cos(qJ(2));
t335 = t278 * t281;
t277 = sin(qJ(2));
t336 = t277 * t282;
t255 = t273 * t336 + t335;
t316 = qJ(3) + pkin(11);
t302 = sin(t316);
t292 = t272 * t302;
t303 = cos(t316);
t223 = t255 * t303 - t282 * t292;
t334 = t281 * t282;
t337 = t277 * t278;
t257 = -t273 * t337 + t334;
t225 = t257 * t303 + t278 * t292;
t293 = t272 * t303;
t241 = t273 * t302 + t277 * t293;
t222 = t255 * t302 + t282 * t293;
t254 = -t273 * t334 + t337;
t275 = sin(qJ(6));
t279 = cos(qJ(6));
t177 = t222 * t279 - t254 * t275;
t178 = t222 * t275 + t254 * t279;
t100 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t223;
t102 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t223;
t224 = t257 * t302 - t278 * t293;
t256 = t273 * t335 + t336;
t179 = t224 * t279 - t256 * t275;
t180 = t224 * t275 + t256 * t279;
t98 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t223;
t37 = t100 * t179 + t102 * t180 + t225 * t98;
t101 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t225;
t103 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t225;
t99 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t225;
t38 = t101 * t179 + t103 * t180 + t225 * t99;
t240 = -t273 * t303 + t277 * t292;
t340 = t272 * t281;
t220 = t240 * t279 + t275 * t340;
t221 = t240 * t275 - t279 * t340;
t122 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t241;
t123 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t241;
t124 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t241;
t48 = t122 * t225 + t123 * t179 + t124 * t180;
t2 = t223 * t37 + t225 * t38 + t241 * t48;
t350 = t2 / 0.2e1;
t349 = t223 / 0.2e1;
t348 = t225 / 0.2e1;
t347 = t241 / 0.2e1;
t280 = cos(qJ(3));
t270 = pkin(3) * t280 + pkin(2);
t346 = -pkin(2) + t270;
t345 = t222 * rSges(6,3);
t199 = Icges(3,5) * t255 - Icges(3,6) * t254 - Icges(3,3) * t339;
t344 = t199 * t282;
t343 = t272 * t277;
t342 = t272 * t278;
t341 = t272 * t280;
t276 = sin(qJ(3));
t338 = t273 * t276;
t294 = -t178 * rSges(7,1) - t177 * rSges(7,2);
t104 = t223 * rSges(7,3) - t294;
t333 = t254 * pkin(5) + t223 * pkin(10) + t104;
t105 = t180 * rSges(7,1) + t179 * rSges(7,2) + t225 * rSges(7,3);
t332 = t256 * pkin(5) + pkin(10) * t225 + t105;
t249 = t254 * pkin(9);
t311 = t276 * t339;
t262 = pkin(3) * t311;
t274 = -qJ(4) - pkin(9);
t148 = -t254 * t274 + t255 * t346 - t249 - t262;
t118 = t256 * t148;
t213 = t222 * qJ(5);
t169 = t223 * pkin(4) + t213;
t331 = t256 * t169 + t118;
t125 = rSges(7,1) * t221 + rSges(7,2) * t220 + rSges(7,3) * t241;
t330 = -pkin(5) * t340 + pkin(10) * t241 + t125;
t197 = pkin(3) * t338 + ((pkin(9) + t274) * t281 + t346 * t277) * t272;
t329 = t148 * t340 + t254 * t197;
t212 = t257 * pkin(2) + pkin(9) * t256;
t312 = t276 * t342;
t305 = pkin(3) * t312 - t256 * t274 + t257 * t270;
t149 = -t212 + t305;
t210 = t273 * t212;
t328 = t273 * t149 + t210;
t146 = t225 * rSges(5,1) - t224 * rSges(5,2) + t256 * rSges(5,3);
t327 = -t146 - t149;
t211 = pkin(2) * t255 + t249;
t326 = -t148 - t211;
t170 = t225 * pkin(4) + qJ(5) * t224;
t325 = -t149 - t170;
t229 = -t255 * t276 - t280 * t339;
t230 = t255 * t280 - t311;
t157 = rSges(4,1) * t230 + rSges(4,2) * t229 + rSges(4,3) * t254;
t324 = -t157 - t211;
t183 = -Icges(6,5) * t340 - Icges(6,6) * t241 + Icges(6,3) * t240;
t184 = -Icges(6,4) * t340 - Icges(6,2) * t241 + Icges(6,6) * t240;
t323 = t240 * t183 - t241 * t184;
t187 = Icges(5,4) * t241 - Icges(5,2) * t240 - Icges(5,6) * t340;
t188 = Icges(5,1) * t241 - Icges(5,4) * t240 - Icges(5,5) * t340;
t322 = -t240 * t187 + t241 * t188;
t252 = t273 * t280 - t276 * t343;
t253 = t277 * t341 + t338;
t195 = Icges(4,4) * t253 + Icges(4,2) * t252 - Icges(4,6) * t340;
t196 = Icges(4,1) * t253 + Icges(4,4) * t252 - Icges(4,5) * t340;
t321 = t252 * t195 + t253 * t196;
t190 = rSges(5,1) * t241 - rSges(5,2) * t240 - rSges(5,3) * t340;
t320 = -t190 - t197;
t193 = pkin(4) * t241 + qJ(5) * t240;
t319 = -t193 - t197;
t318 = t211 * t342 + t212 * t339;
t317 = t282 * pkin(1) + pkin(8) * t342;
t40 = t100 * t220 + t102 * t221 + t241 * t98;
t47 = t122 * t223 + t123 * t177 + t124 * t178;
t314 = t40 / 0.2e1 + t47 / 0.2e1;
t41 = t101 * t220 + t103 * t221 + t241 * t99;
t313 = t41 / 0.2e1 + t48 / 0.2e1;
t52 = t241 * t122 + t220 * t123 + t221 * t124;
t310 = t273 * t170 + t328;
t144 = t256 * rSges(6,1) - t225 * rSges(6,2) + t224 * rSges(6,3);
t309 = -t144 + t325;
t308 = -t169 + t326;
t189 = -rSges(6,1) * t340 - rSges(6,2) * t241 + rSges(6,3) * t240;
t307 = -t189 + t319;
t231 = -t257 * t276 + t278 * t341;
t232 = t257 * t280 + t312;
t158 = t232 * rSges(4,1) + t231 * rSges(4,2) + t256 * rSges(4,3);
t237 = Icges(3,3) * t273 + (Icges(3,5) * t277 + Icges(3,6) * t281) * t272;
t238 = Icges(3,6) * t273 + (Icges(3,4) * t277 + Icges(3,2) * t281) * t272;
t239 = Icges(3,5) * t273 + (Icges(3,1) * t277 + Icges(3,4) * t281) * t272;
t306 = t273 * t237 + t238 * t340 + t239 * t343;
t206 = t257 * rSges(3,1) - t256 * rSges(3,2) + rSges(3,3) * t342;
t304 = -t278 * pkin(1) + pkin(8) * t339;
t198 = rSges(4,1) * t253 + rSges(4,2) * t252 - rSges(4,3) * t340;
t258 = (pkin(2) * t277 - pkin(9) * t281) * t272;
t301 = t272 * (-t198 - t258);
t300 = t325 - t332;
t299 = t319 - t330;
t298 = t148 * t342 + t149 * t339 + t318;
t297 = t169 * t340 + t254 * t193 + t329;
t296 = t272 * (-t258 + t320);
t295 = -t223 * rSges(5,1) + t222 * rSges(5,2);
t291 = t305 + t317;
t290 = t272 * (-t258 + t307);
t289 = t169 * t342 + t170 * t339 + t298;
t288 = t272 * (-t258 + t299);
t287 = -t255 * t270 + t262 + t304;
t286 = -t213 + t287;
t205 = t255 * rSges(3,1) - t254 * rSges(3,2) - rSges(3,3) * t339;
t285 = t170 + t291;
t131 = Icges(6,5) * t254 - Icges(6,6) * t223 + Icges(6,3) * t222;
t135 = Icges(6,4) * t254 - Icges(6,2) * t223 + Icges(6,6) * t222;
t139 = Icges(6,1) * t254 - Icges(6,4) * t223 + Icges(6,5) * t222;
t68 = t131 * t240 - t135 * t241 - t139 * t340;
t133 = Icges(5,5) * t223 - Icges(5,6) * t222 + Icges(5,3) * t254;
t137 = Icges(5,4) * t223 - Icges(5,2) * t222 + Icges(5,6) * t254;
t141 = Icges(5,1) * t223 - Icges(5,4) * t222 + Icges(5,5) * t254;
t70 = -t133 * t340 - t137 * t240 + t141 * t241;
t151 = Icges(4,5) * t230 + Icges(4,6) * t229 + Icges(4,3) * t254;
t153 = Icges(4,4) * t230 + Icges(4,2) * t229 + Icges(4,6) * t254;
t155 = Icges(4,1) * t230 + Icges(4,4) * t229 + Icges(4,5) * t254;
t78 = -t151 * t340 + t153 * t252 + t155 * t253;
t185 = -Icges(6,1) * t340 - Icges(6,4) * t241 + Icges(6,5) * t240;
t81 = t183 * t222 - t184 * t223 + t185 * t254;
t186 = Icges(5,5) * t241 - Icges(5,6) * t240 - Icges(5,3) * t340;
t83 = t186 * t254 - t187 * t222 + t188 * t223;
t194 = Icges(4,5) * t253 + Icges(4,6) * t252 - Icges(4,3) * t340;
t91 = t194 * t254 + t195 * t229 + t196 * t230;
t284 = t70 / 0.2e1 + t68 / 0.2e1 + t91 / 0.2e1 + t83 / 0.2e1 + t81 / 0.2e1 + t78 / 0.2e1 + t314;
t132 = Icges(6,5) * t256 - Icges(6,6) * t225 + Icges(6,3) * t224;
t136 = Icges(6,4) * t256 - Icges(6,2) * t225 + Icges(6,6) * t224;
t140 = Icges(6,1) * t256 - Icges(6,4) * t225 + Icges(6,5) * t224;
t69 = t132 * t240 - t136 * t241 - t140 * t340;
t134 = Icges(5,5) * t225 - Icges(5,6) * t224 + Icges(5,3) * t256;
t138 = Icges(5,4) * t225 - Icges(5,2) * t224 + Icges(5,6) * t256;
t142 = Icges(5,1) * t225 - Icges(5,4) * t224 + Icges(5,5) * t256;
t71 = -t134 * t340 - t138 * t240 + t142 * t241;
t152 = Icges(4,5) * t232 + Icges(4,6) * t231 + Icges(4,3) * t256;
t154 = Icges(4,4) * t232 + Icges(4,2) * t231 + Icges(4,6) * t256;
t156 = Icges(4,1) * t232 + Icges(4,4) * t231 + Icges(4,5) * t256;
t79 = -t152 * t340 + t154 * t252 + t156 * t253;
t82 = t183 * t224 - t184 * t225 + t185 * t256;
t84 = t186 * t256 - t187 * t224 + t188 * t225;
t92 = t194 * t256 + t195 * t231 + t196 * t232;
t283 = t69 / 0.2e1 + t92 / 0.2e1 + t84 / 0.2e1 + t82 / 0.2e1 + t79 / 0.2e1 + t71 / 0.2e1 + t313;
t264 = rSges(2,1) * t282 - t278 * rSges(2,2);
t263 = -t278 * rSges(2,1) - rSges(2,2) * t282;
t242 = t273 * rSges(3,3) + (rSges(3,1) * t277 + rSges(3,2) * t281) * t272;
t204 = Icges(3,1) * t257 - Icges(3,4) * t256 + Icges(3,5) * t342;
t203 = Icges(3,1) * t255 - Icges(3,4) * t254 - Icges(3,5) * t339;
t202 = Icges(3,4) * t257 - Icges(3,2) * t256 + Icges(3,6) * t342;
t201 = Icges(3,4) * t255 - Icges(3,2) * t254 - Icges(3,6) * t339;
t200 = Icges(3,5) * t257 - Icges(3,6) * t256 + Icges(3,3) * t342;
t182 = t206 + t317;
t181 = -t205 + t304;
t168 = -t273 * t205 - t242 * t339;
t167 = t206 * t273 - t242 * t342;
t150 = t306 * t273;
t145 = rSges(5,3) * t254 - t295;
t143 = t254 * rSges(6,1) - t223 * rSges(6,2) + t345;
t121 = (t205 * t278 + t206 * t282) * t272;
t120 = t237 * t342 - t238 * t256 + t239 * t257;
t119 = -t237 * t339 - t254 * t238 + t255 * t239;
t116 = t212 + t158 + t317;
t115 = t304 + t324;
t111 = -t158 * t340 - t198 * t256;
t110 = t157 * t340 + t198 * t254;
t109 = t273 * t200 + (t202 * t281 + t204 * t277) * t272;
t108 = t273 * t199 + (t201 * t281 + t203 * t277) * t272;
t107 = t291 + t146;
t106 = (-rSges(5,3) + t274) * t254 + t287 + t295;
t97 = -t194 * t340 + t321;
t96 = t97 * t273;
t95 = t157 * t256 - t158 * t254;
t94 = t273 * t324 + t282 * t301;
t93 = t273 * t158 + t278 * t301 + t210;
t90 = -t186 * t340 + t322;
t89 = -t185 * t340 + t323;
t88 = t90 * t273;
t87 = t89 * t273;
t86 = t285 + t144;
t85 = -t345 + (-rSges(6,1) + t274) * t254 + (rSges(6,2) - pkin(4)) * t223 + t286;
t80 = (t157 * t278 + t158 * t282) * t272 + t318;
t77 = t105 * t241 - t125 * t225;
t76 = -t104 * t241 + t125 * t223;
t75 = t152 * t256 + t154 * t231 + t156 * t232;
t74 = t151 * t256 + t153 * t231 + t155 * t232;
t73 = t152 * t254 + t154 * t229 + t156 * t230;
t72 = t151 * t254 + t153 * t229 + t155 * t230;
t67 = t256 * t320 + t327 * t340;
t66 = t145 * t340 + t190 * t254 + t329;
t65 = t285 + t332;
t64 = (-pkin(5) + t274) * t254 + (-rSges(7,3) - pkin(4) - pkin(10)) * t223 + t286 + t294;
t63 = t134 * t256 - t138 * t224 + t142 * t225;
t62 = t133 * t256 - t137 * t224 + t141 * t225;
t61 = t134 * t254 - t138 * t222 + t142 * t223;
t60 = t133 * t254 - t137 * t222 + t141 * t223;
t59 = t132 * t224 - t136 * t225 + t140 * t256;
t58 = t131 * t224 - t135 * t225 + t139 * t256;
t57 = t132 * t222 - t136 * t223 + t140 * t254;
t56 = t131 * t222 - t135 * t223 + t139 * t254;
t55 = (-t145 + t326) * t273 + t282 * t296;
t54 = t273 * t146 + t278 * t296 + t328;
t53 = t104 * t225 - t105 * t223;
t51 = t52 * t273;
t50 = t52 * t241;
t49 = t256 * t145 + t254 * t327 + t118;
t46 = t256 * t307 + t309 * t340;
t45 = t143 * t340 + t189 * t254 + t297;
t44 = (t145 * t278 + t146 * t282) * t272 + t298;
t43 = (-t143 + t308) * t273 + t282 * t290;
t42 = t273 * t144 + t278 * t290 + t310;
t39 = t256 * t143 + t254 * t309 + t331;
t36 = t101 * t177 + t103 * t178 + t223 * t99;
t35 = t100 * t177 + t102 * t178 + t223 * t98;
t34 = (t143 * t278 + t144 * t282) * t272 + t289;
t33 = t256 * t299 + t300 * t340;
t32 = t254 * t330 + t333 * t340 + t297;
t31 = (t308 - t333) * t273 + t282 * t288;
t30 = t273 * t332 + t278 * t288 + t310;
t29 = t96 + (t79 * t278 - t78 * t282) * t272;
t28 = t78 * t254 + t79 * t256 - t340 * t97;
t27 = t92 * t273 + (t278 * t75 - t282 * t74) * t272;
t26 = t91 * t273 + (t278 * t73 - t282 * t72) * t272;
t25 = t254 * t300 + t256 * t333 + t331;
t24 = t88 + (t71 * t278 - t70 * t282) * t272;
t23 = t87 + (t69 * t278 - t68 * t282) * t272;
t22 = (t278 * t333 + t282 * t332) * t272 + t289;
t21 = t254 * t74 + t256 * t75 - t340 * t92;
t20 = t254 * t72 + t256 * t73 - t340 * t91;
t19 = t70 * t254 + t71 * t256 - t340 * t90;
t18 = t68 * t254 + t69 * t256 - t340 * t89;
t17 = t84 * t273 + (t278 * t63 - t282 * t62) * t272;
t16 = t83 * t273 + (t278 * t61 - t282 * t60) * t272;
t15 = t82 * t273 + (t278 * t59 - t282 * t58) * t272;
t14 = t81 * t273 + (t278 * t57 - t282 * t56) * t272;
t13 = t254 * t62 + t256 * t63 - t340 * t84;
t12 = t254 * t60 + t256 * t61 - t340 * t83;
t11 = t254 * t58 + t256 * t59 - t340 * t82;
t10 = t254 * t56 + t256 * t57 - t340 * t81;
t9 = t51 + (t41 * t278 - t40 * t282) * t272;
t8 = t40 * t254 + t41 * t256 - t340 * t52;
t7 = t40 * t223 + t41 * t225 + t50;
t6 = t48 * t273 + (t278 * t38 - t282 * t37) * t272;
t5 = t47 * t273 + (t278 * t36 - t282 * t35) * t272;
t4 = t254 * t37 + t256 * t38 - t340 * t48;
t3 = t254 * t35 + t256 * t36 - t340 * t47;
t1 = t223 * t35 + t225 * t36 + t241 * t47;
t112 = [Icges(2,3) + (-t185 - t186 - t194) * t340 + m(7) * (t64 ^ 2 + t65 ^ 2) + m(5) * (t106 ^ 2 + t107 ^ 2) + m(6) * (t85 ^ 2 + t86 ^ 2) + m(4) * (t115 ^ 2 + t116 ^ 2) + m(3) * (t181 ^ 2 + t182 ^ 2) + m(2) * (t263 ^ 2 + t264 ^ 2) + t306 + t52 + t321 + t322 + t323; t96 + t51 + t88 + t87 + t150 + m(7) * (t30 * t65 + t31 * t64) + m(6) * (t42 * t86 + t43 * t85) + m(5) * (t106 * t55 + t107 * t54) + m(4) * (t115 * t94 + t116 * t93) + m(3) * (t167 * t182 + t168 * t181) + ((-t108 / 0.2e1 - t119 / 0.2e1 - t284) * t282 + (t109 / 0.2e1 + t120 / 0.2e1 + t283) * t278) * t272; (t9 + t24 + t23 + t29 + t150) * t273 + m(7) * (t22 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t34 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t44 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t80 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(3) * (t121 ^ 2 + t167 ^ 2 + t168 ^ 2) + ((-t5 - t26 - t16 - t14 + ((-t254 * t201 + t255 * t203) * t272 - t351 * t344) * t282) * t282 + (t6 + t17 + t15 + t27 + ((-t202 * t256 + t204 * t257 + (t200 * t278 - t344) * t272) * t278 + (t200 * t339 + t201 * t256 + t254 * t202 - t203 * t257 - t255 * t204) * t282) * t272) * t278 + ((-t108 - t119) * t282 + (t109 + t120) * t278) * t273) * t272; (-t52 - t89 - t90 - t97) * t340 + m(7) * (t32 * t64 + t33 * t65) + m(5) * (t106 * t66 + t107 * t67) + m(6) * (t45 * t85 + t46 * t86) + m(4) * (t110 * t115 + t111 * t116) + t283 * t256 + t284 * t254; (t8 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1 + t28 / 0.2e1) * t273 + (t6 / 0.2e1 + t15 / 0.2e1 + t17 / 0.2e1 + t27 / 0.2e1) * t256 + (t5 / 0.2e1 + t14 / 0.2e1 + t16 / 0.2e1 + t26 / 0.2e1) * t254 + m(7) * (t22 * t25 + t30 * t33 + t31 * t32) + m(6) * (t34 * t39 + t42 * t46 + t43 * t45) + m(5) * (t49 * t44 + t54 * t67 + t55 * t66) + m(4) * (t110 * t94 + t111 * t93 + t80 * t95) + ((-t3 / 0.2e1 - t10 / 0.2e1 - t12 / 0.2e1 - t20 / 0.2e1) * t282 + (-t9 / 0.2e1 - t23 / 0.2e1 - t24 / 0.2e1 - t29 / 0.2e1) * t281 + (t4 / 0.2e1 + t11 / 0.2e1 + t13 / 0.2e1 + t21 / 0.2e1) * t278) * t272; (-t18 - t19 - t28 - t8) * t340 + (t4 + t13 + t11 + t21) * t256 + (t3 + t12 + t20 + t10) * t254 + m(7) * (t25 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t49 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(4) * (t110 ^ 2 + t111 ^ 2 + t95 ^ 2); m(7) * (t254 * t65 + t256 * t64) + m(5) * (t106 * t256 + t107 * t254) + m(6) * (t254 * t86 + t256 * t85); m(7) * (-t22 * t340 + t254 * t30 + t256 * t31) + m(6) * (t254 * t42 + t256 * t43 - t34 * t340) + m(5) * (t254 * t54 + t256 * t55 - t340 * t44); m(7) * (-t25 * t340 + t254 * t33 + t256 * t32) + m(6) * (t254 * t46 + t256 * t45 - t340 * t39) + m(5) * (t254 * t67 + t256 * t66 - t340 * t49); 0.2e1 * (m(5) / 0.2e1 + t315) * (t281 ^ 2 * t351 + t254 ^ 2 + t256 ^ 2); m(7) * (t222 * t65 + t224 * t64) + m(6) * (t222 * t86 + t224 * t85); m(7) * (t22 * t240 + t222 * t30 + t224 * t31) + m(6) * (t222 * t42 + t224 * t43 + t240 * t34); m(7) * (t222 * t33 + t224 * t32 + t240 * t25) + m(6) * (t222 * t46 + t224 * t45 + t240 * t39); (t222 * t254 + t224 * t256 - t240 * t340) * t352; (t222 ^ 2 + t224 ^ 2 + t240 ^ 2) * t352; m(7) * (t64 * t76 + t65 * t77) + t50 + t313 * t225 + t314 * t223; m(7) * (t53 * t22 + t30 * t77 + t31 * t76) + t5 * t349 + t9 * t347 + t6 * t348 + t273 * t7 / 0.2e1 + (t278 * t350 - t282 * t1 / 0.2e1) * t272; -t7 * t340 / 0.2e1 + m(7) * (t53 * t25 + t32 * t76 + t33 * t77) + t8 * t347 + t256 * t350 + t254 * t1 / 0.2e1 + t4 * t348 + t3 * t349; m(7) * (t254 * t77 + t256 * t76 - t340 * t53); m(7) * (t222 * t77 + t224 * t76 + t240 * t53); t225 * t2 + t223 * t1 + t241 * t7 + m(7) * (t53 ^ 2 + t76 ^ 2 + t77 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t112(1) t112(2) t112(4) t112(7) t112(11) t112(16); t112(2) t112(3) t112(5) t112(8) t112(12) t112(17); t112(4) t112(5) t112(6) t112(9) t112(13) t112(18); t112(7) t112(8) t112(9) t112(10) t112(14) t112(19); t112(11) t112(12) t112(13) t112(14) t112(15) t112(20); t112(16) t112(17) t112(18) t112(19) t112(20) t112(21);];
Mq  = res;
