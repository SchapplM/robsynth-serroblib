% Calculate joint inertia matrix for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:52
% EndTime: 2019-03-09 15:32:01
% DurationCPUTime: 4.70s
% Computational Cost: add. (11654->536), mult. (16446->769), div. (0->0), fcn. (18966->10), ass. (0->261)
t238 = sin(qJ(2));
t334 = Icges(3,5) * t238;
t333 = t334 / 0.2e1;
t231 = qJ(3) + pkin(10);
t225 = sin(t231);
t226 = cos(t231);
t242 = cos(qJ(2));
t164 = -Icges(5,3) * t242 + (Icges(5,5) * t226 - Icges(5,6) * t225) * t238;
t165 = -Icges(6,2) * t242 + (Icges(6,4) * t226 + Icges(6,6) * t225) * t238;
t237 = sin(qJ(3));
t241 = cos(qJ(3));
t184 = -Icges(4,3) * t242 + (Icges(4,5) * t241 - Icges(4,6) * t237) * t238;
t332 = -t164 - t165 - t184;
t166 = -Icges(5,6) * t242 + (Icges(5,4) * t226 - Icges(5,2) * t225) * t238;
t302 = t225 * t238;
t187 = -Icges(4,6) * t242 + (Icges(4,4) * t241 - Icges(4,2) * t237) * t238;
t303 = t187 * t237;
t163 = -Icges(6,6) * t242 + (Icges(6,5) * t226 + Icges(6,3) * t225) * t238;
t167 = -Icges(6,4) * t242 + (Icges(6,1) * t226 + Icges(6,5) * t225) * t238;
t168 = -Icges(5,5) * t242 + (Icges(5,1) * t226 - Icges(5,4) * t225) * t238;
t190 = -Icges(4,5) * t242 + (Icges(4,1) * t241 - Icges(4,4) * t237) * t238;
t301 = t226 * t238;
t326 = t238 * t241 * t190 + t163 * t302 + (t167 + t168) * t301;
t331 = (-t166 * t302 - t238 * t303 + t242 * t332 + t326) * t242;
t239 = sin(qJ(1));
t243 = cos(qJ(1));
t296 = t239 * t242;
t180 = t225 * t296 + t226 * t243;
t181 = -t225 * t243 + t226 * t296;
t299 = t238 * t239;
t106 = Icges(6,5) * t181 + Icges(6,6) * t299 + Icges(6,3) * t180;
t110 = Icges(6,4) * t181 + Icges(6,2) * t299 + Icges(6,6) * t180;
t114 = Icges(6,1) * t181 + Icges(6,4) * t299 + Icges(6,5) * t180;
t43 = t106 * t180 + t110 * t299 + t114 * t181;
t295 = t242 * t243;
t182 = t225 * t295 - t226 * t239;
t183 = t225 * t239 + t226 * t295;
t298 = t238 * t243;
t107 = Icges(6,5) * t183 + Icges(6,6) * t298 + Icges(6,3) * t182;
t111 = Icges(6,4) * t183 + Icges(6,2) * t298 + Icges(6,6) * t182;
t115 = Icges(6,1) * t183 + Icges(6,4) * t298 + Icges(6,5) * t182;
t44 = t107 * t180 + t111 * t299 + t115 * t181;
t108 = Icges(5,5) * t181 - Icges(5,6) * t180 + Icges(5,3) * t299;
t112 = Icges(5,4) * t181 - Icges(5,2) * t180 + Icges(5,6) * t299;
t116 = Icges(5,1) * t181 - Icges(5,4) * t180 + Icges(5,5) * t299;
t45 = t108 * t299 - t112 * t180 + t116 * t181;
t109 = Icges(5,5) * t183 - Icges(5,6) * t182 + Icges(5,3) * t298;
t113 = Icges(5,4) * t183 - Icges(5,2) * t182 + Icges(5,6) * t298;
t117 = Icges(5,1) * t183 - Icges(5,4) * t182 + Icges(5,5) * t298;
t46 = t109 * t299 - t113 * t180 + t117 * t181;
t200 = -t237 * t296 - t241 * t243;
t300 = t237 * t243;
t201 = t241 * t296 - t300;
t140 = Icges(4,5) * t201 + Icges(4,6) * t200 + Icges(4,3) * t299;
t142 = Icges(4,4) * t201 + Icges(4,2) * t200 + Icges(4,6) * t299;
t144 = Icges(4,1) * t201 + Icges(4,4) * t200 + Icges(4,5) * t299;
t57 = t140 * t299 + t142 * t200 + t144 * t201;
t202 = -t237 * t295 + t239 * t241;
t297 = t239 * t237;
t203 = t241 * t295 + t297;
t141 = Icges(4,5) * t203 + Icges(4,6) * t202 + Icges(4,3) * t298;
t143 = Icges(4,4) * t203 + Icges(4,2) * t202 + Icges(4,6) * t298;
t145 = Icges(4,1) * t203 + Icges(4,4) * t202 + Icges(4,5) * t298;
t58 = t141 * t299 + t143 * t200 + t145 * t201;
t70 = t163 * t180 + t165 * t299 + t167 * t181;
t71 = t164 * t299 - t166 * t180 + t168 * t181;
t82 = t184 * t299 + t187 * t200 + t190 * t201;
t330 = (-t71 - t82 - t70) * t242 + ((t44 + t46 + t58) * t243 + (t43 + t45 + t57) * t239) * t238;
t47 = t106 * t182 + t110 * t298 + t114 * t183;
t48 = t107 * t182 + t111 * t298 + t115 * t183;
t49 = t108 * t298 - t112 * t182 + t116 * t183;
t50 = t109 * t298 - t113 * t182 + t117 * t183;
t59 = t140 * t298 + t142 * t202 + t144 * t203;
t60 = t141 * t298 + t143 * t202 + t145 * t203;
t72 = t163 * t182 + t165 * t298 + t167 * t183;
t73 = t164 * t298 - t166 * t182 + t168 * t183;
t83 = t184 * t298 + t187 * t202 + t190 * t203;
t329 = (-t83 - t73 - t72) * t242 + ((t48 + t50 + t60) * t243 + (t47 + t49 + t59) * t239) * t238;
t51 = -t110 * t242 + (t106 * t225 + t114 * t226) * t238;
t53 = -t108 * t242 + (-t112 * t225 + t116 * t226) * t238;
t63 = -t140 * t242 + (-t142 * t237 + t144 * t241) * t238;
t328 = -t51 - t53 - t63;
t52 = -t111 * t242 + (t107 * t225 + t115 * t226) * t238;
t54 = -t109 * t242 + (-t113 * t225 + t117 * t226) * t238;
t64 = -t141 * t242 + (-t143 * t237 + t145 * t241) * t238;
t327 = t52 + t54 + t64;
t236 = sin(qJ(6));
t240 = cos(qJ(6));
t130 = t180 * t240 - t181 * t236;
t131 = t180 * t236 + t181 * t240;
t74 = Icges(7,5) * t131 + Icges(7,6) * t130 - Icges(7,3) * t299;
t76 = Icges(7,4) * t131 + Icges(7,2) * t130 - Icges(7,6) * t299;
t78 = Icges(7,1) * t131 + Icges(7,4) * t130 - Icges(7,5) * t299;
t15 = t130 * t76 + t131 * t78 - t299 * t74;
t132 = t182 * t240 - t183 * t236;
t133 = t182 * t236 + t183 * t240;
t75 = Icges(7,5) * t133 + Icges(7,6) * t132 - Icges(7,3) * t298;
t77 = Icges(7,4) * t133 + Icges(7,2) * t132 - Icges(7,6) * t298;
t79 = Icges(7,1) * t133 + Icges(7,4) * t132 - Icges(7,5) * t298;
t16 = t130 * t77 + t131 * t79 - t299 * t75;
t171 = (t225 * t240 - t226 * t236) * t238;
t172 = (t225 * t236 + t226 * t240) * t238;
t102 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t242;
t103 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t242;
t104 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t242;
t33 = -t102 * t299 + t103 * t130 + t104 * t131;
t1 = -t33 * t242 + (t15 * t239 + t16 * t243) * t238;
t17 = t132 * t76 + t133 * t78 - t298 * t74;
t18 = t132 * t77 + t133 * t79 - t298 * t75;
t34 = -t102 * t298 + t103 * t132 + t104 * t133;
t2 = -t34 * t242 + (t17 * t239 + t18 * t243) * t238;
t25 = t171 * t76 + t172 * t78 + t242 * t74;
t26 = t171 * t77 + t172 * t79 + t242 * t75;
t276 = t102 * t242 + t103 * t171 + t104 * t172;
t39 = t276 * t242;
t3 = -t39 + (t25 * t239 + t26 * t243) * t238;
t325 = (t239 * t1 + t243 * t2) * t238 - t242 * t3;
t233 = t239 ^ 2;
t234 = t243 ^ 2;
t324 = 0.2e1 * t238;
t323 = m(5) / 0.2e1;
t322 = m(6) / 0.2e1;
t321 = m(7) / 0.2e1;
t320 = -t239 / 0.2e1;
t319 = t239 / 0.2e1;
t318 = -t242 / 0.2e1;
t317 = t242 / 0.2e1;
t316 = -t243 / 0.2e1;
t315 = t243 / 0.2e1;
t314 = rSges(7,3) + pkin(9);
t313 = pkin(2) * t242;
t312 = pkin(8) * t238;
t224 = pkin(3) * t241 + pkin(2);
t310 = -pkin(2) + t224;
t309 = rSges(6,3) * t180;
t308 = t243 * rSges(3,3);
t257 = -rSges(7,1) * t131 - rSges(7,2) * t130;
t80 = -rSges(7,3) * t299 - t257;
t307 = pkin(5) * t181 - pkin(9) * t299 + t80;
t178 = t183 * pkin(5);
t291 = rSges(7,1) * t133 + rSges(7,2) * t132;
t81 = -rSges(7,3) * t298 + t291;
t306 = -pkin(9) * t298 + t178 + t81;
t304 = Icges(3,4) * t242;
t105 = rSges(7,1) * t172 + rSges(7,2) * t171 + rSges(7,3) * t242;
t294 = pkin(5) * t301 + pkin(9) * t242 + t105;
t121 = rSges(5,1) * t183 - rSges(5,2) * t182 + rSges(5,3) * t298;
t235 = -qJ(4) - pkin(8);
t273 = t235 * t298;
t283 = pkin(3) * t297 + t224 * t295;
t249 = -t273 + t283;
t281 = pkin(2) * t295 + pkin(8) * t298;
t147 = t249 - t281;
t293 = -t121 - t147;
t282 = -pkin(3) * t300 - t235 * t299;
t146 = (t242 * t310 - t312) * t239 + t282;
t126 = t146 * t298;
t173 = t180 * qJ(5);
t138 = pkin(4) * t181 + t173;
t292 = t138 * t298 + t126;
t162 = (pkin(8) + t235) * t242 + t310 * t238;
t290 = t146 * t242 + t162 * t299;
t139 = pkin(4) * t183 + t182 * qJ(5);
t289 = -t139 - t147;
t170 = -t242 * rSges(5,3) + (rSges(5,1) * t226 - rSges(5,2) * t225) * t238;
t287 = -t162 - t170;
t193 = (pkin(4) * t226 + qJ(5) * t225) * t238;
t286 = -t162 - t193;
t194 = -t242 * rSges(4,3) + (rSges(4,1) * t241 - rSges(4,2) * t237) * t238;
t213 = pkin(2) * t238 - pkin(8) * t242;
t285 = -t194 - t213;
t284 = t233 * (t312 + t313) + t243 * t281;
t280 = pkin(1) * t243 + pkin(7) * t239;
t279 = t233 + t234;
t278 = t322 + t321;
t275 = -t25 / 0.2e1 - t33 / 0.2e1;
t274 = -t26 / 0.2e1 - t34 / 0.2e1;
t120 = rSges(6,1) * t183 + rSges(6,2) * t298 + rSges(6,3) * t182;
t272 = -t120 + t289;
t169 = -t242 * rSges(6,2) + (rSges(6,1) * t226 + rSges(6,3) * t225) * t238;
t271 = -t169 + t286;
t270 = -t213 + t287;
t149 = rSges(4,1) * t203 + rSges(4,2) * t202 + rSges(4,3) * t298;
t229 = t243 * pkin(7);
t269 = t229 - t282;
t268 = -t224 * t242 - pkin(1);
t267 = t289 - t306;
t266 = t286 - t294;
t265 = t138 * t242 + t193 * t299 + t290;
t264 = t146 * t239 + t147 * t243 + t284;
t263 = -t213 + t271;
t262 = -t173 + t269;
t260 = rSges(3,1) * t242 - rSges(3,2) * t238;
t259 = -rSges(4,1) * t201 - rSges(4,2) * t200;
t258 = -rSges(5,1) * t181 + rSges(5,2) * t180;
t256 = -t213 + t266;
t254 = -Icges(3,2) * t238 + t304;
t253 = Icges(3,5) * t242 - Icges(3,6) * t238;
t250 = rSges(3,1) * t295 - rSges(3,2) * t298 + rSges(3,3) * t239;
t248 = t138 * t239 + t139 * t243 + t264;
t247 = t139 + t280 + t283;
t246 = t1 * t315 + t2 * t320 + (t26 * t239 - t25 * t243) * t317;
t245 = t51 / 0.2e1 + t82 / 0.2e1 + t71 / 0.2e1 + t70 / 0.2e1 + t63 / 0.2e1 + t53 / 0.2e1 - t275;
t244 = t64 / 0.2e1 + t54 / 0.2e1 + t52 / 0.2e1 + t83 / 0.2e1 + t73 / 0.2e1 + t72 / 0.2e1 - t274;
t232 = t238 ^ 2;
t212 = rSges(2,1) * t243 - rSges(2,2) * t239;
t211 = -rSges(2,1) * t239 - rSges(2,2) * t243;
t210 = rSges(3,1) * t238 + rSges(3,2) * t242;
t207 = Icges(3,6) * t242 + t334;
t186 = Icges(3,3) * t239 + t243 * t253;
t185 = -Icges(3,3) * t243 + t239 * t253;
t159 = t250 + t280;
t158 = t308 + t229 + (-pkin(1) - t260) * t239;
t151 = t285 * t243;
t150 = t285 * t239;
t148 = rSges(4,3) * t299 - t259;
t134 = t243 * t250 + (t239 * t260 - t308) * t239;
t119 = rSges(5,3) * t299 - t258;
t118 = rSges(6,1) * t181 + rSges(6,2) * t299 + t309;
t99 = t149 + t280 + t281;
t98 = t229 + (-t313 - pkin(1) + (-rSges(4,3) - pkin(8)) * t238) * t239 + t259;
t97 = t270 * t243;
t96 = t270 * t239;
t95 = -t149 * t242 - t194 * t298;
t94 = t148 * t242 + t194 * t299;
t90 = t249 + t121 + t280;
t89 = (-rSges(5,3) * t238 + t268) * t239 + t258 + t269;
t88 = t263 * t243;
t87 = t263 * t239;
t84 = (t148 * t243 - t149 * t239) * t238;
t69 = t148 * t239 + t149 * t243 + t284;
t68 = t120 + t247 - t273;
t67 = -t309 + (-rSges(6,1) - pkin(4)) * t181 + (-rSges(6,2) * t238 + t268) * t239 + t262;
t66 = t256 * t243;
t65 = t256 * t239;
t62 = t242 * t293 + t287 * t298;
t61 = t119 * t242 + t170 * t299 + t290;
t56 = t105 * t298 + t242 * t81;
t55 = -t105 * t299 - t242 * t80;
t42 = t126 + (t119 * t243 + t239 * t293) * t238;
t41 = t178 + (-t235 - t314) * t298 + t247 + t291;
t40 = (-pkin(4) - pkin(5)) * t181 + (t238 * t314 + t268) * t239 + t257 + t262;
t38 = (t239 * t81 - t243 * t80) * t238;
t37 = t242 * t272 + t271 * t298;
t36 = t118 * t242 + t169 * t299 + t265;
t35 = t119 * t239 + t121 * t243 + t264;
t30 = (t118 * t243 + t239 * t272) * t238 + t292;
t29 = t118 * t239 + t120 * t243 + t248;
t28 = t239 * t60 - t243 * t59;
t27 = t239 * t58 - t243 * t57;
t24 = t242 * t267 + t266 * t298;
t23 = t242 * t307 + t294 * t299 + t265;
t22 = t239 * t50 - t243 * t49;
t21 = t239 * t48 - t243 * t47;
t20 = t239 * t46 - t243 * t45;
t19 = t239 * t44 - t243 * t43;
t12 = (t239 * t267 + t243 * t307) * t238 + t292;
t11 = t239 * t307 + t243 * t306 + t248;
t5 = -t17 * t243 + t18 * t239;
t4 = -t15 * t243 + t16 * t239;
t6 = [Icges(2,3) + (Icges(3,1) * t238 - t166 * t225 - t303 + t304) * t238 + (Icges(3,4) * t238 + Icges(3,2) * t242 + t332) * t242 + m(7) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t89 ^ 2 + t90 ^ 2) + m(6) * (t67 ^ 2 + t68 ^ 2) + m(4) * (t98 ^ 2 + t99 ^ 2) + m(3) * (t158 ^ 2 + t159 ^ 2) + m(2) * (t211 ^ 2 + t212 ^ 2) + t276 + t326; ((-Icges(3,6) * t243 + t239 * t254) * t318 + t243 * t333 + t207 * t315 - t245) * t243 + ((Icges(3,6) * t239 + t243 * t254) * t317 + t239 * t333 + t207 * t319 + t244) * t239 + m(6) * (t67 * t88 + t68 * t87) + m(7) * (t40 * t66 + t41 * t65) + m(5) * (t89 * t97 + t90 * t96) + m(4) * (t150 * t99 + t151 * t98) + m(3) * (-t158 * t243 - t159 * t239) * t210; m(3) * (t210 ^ 2 * t279 + t134 ^ 2) + m(7) * (t11 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t29 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t35 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(4) * (t150 ^ 2 + t151 ^ 2 + t69 ^ 2) + (-t234 * t185 - t19 - t20 - t27 - t4) * t243 + (t233 * t186 + t21 + t22 + t28 + t5 + (-t239 * t185 + t243 * t186) * t243) * t239; -t39 - t331 + m(7) * (t23 * t40 + t24 * t41) + m(5) * (t61 * t89 + t62 * t90) + m(6) * (t36 * t67 + t37 * t68) + m(4) * (t94 * t98 + t95 * t99) + (t239 * t245 + t243 * t244) * t238; m(7) * (t11 * t12 + t23 * t66 + t24 * t65) + m(6) * (t29 * t30 + t36 * t88 + t37 * t87) + m(5) * (t35 * t42 + t61 * t97 + t62 * t96) + m(4) * (t150 * t95 + t151 * t94 + t69 * t84) + ((t5 / 0.2e1 + t28 / 0.2e1 + t22 / 0.2e1 + t21 / 0.2e1) * t243 + (t4 / 0.2e1 + t20 / 0.2e1 + t19 / 0.2e1 + t27 / 0.2e1) * t239) * t238 - t246 + t329 * t319 + (t239 * t327 + t243 * t328) * t318 + t330 * t316; m(4) * (t84 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(7) * (t12 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t30 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t42 ^ 2 + t61 ^ 2 + t62 ^ 2) + (-t3 + t331) * t242 + ((-t242 * t327 + t2 + t329) * t243 + (t242 * t328 + t1 + t330) * t239) * t238; ((t239 * t41 + t243 * t40) * t321 + (t239 * t90 + t243 * t89) * t323 + (t239 * t68 + t243 * t67) * t322) * t324; m(7) * (-t242 * t11 + (t239 * t65 + t243 * t66) * t238) + m(6) * (-t242 * t29 + (t239 * t87 + t243 * t88) * t238) + m(5) * (-t242 * t35 + (t239 * t96 + t243 * t97) * t238); m(7) * (-t242 * t12 + (t23 * t243 + t239 * t24) * t238) + m(6) * (-t242 * t30 + (t239 * t37 + t243 * t36) * t238) + m(5) * (-t242 * t42 + (t239 * t62 + t243 * t61) * t238); 0.2e1 * (t323 + t278) * (t232 * t279 + t242 ^ 2); m(7) * (t180 * t41 + t182 * t40) + m(6) * (t180 * t68 + t182 * t67); m(7) * (t11 * t302 + t180 * t65 + t182 * t66) + m(6) * (t180 * t87 + t182 * t88 + t29 * t302); m(7) * (t12 * t302 + t180 * t24 + t182 * t23) + m(6) * (t180 * t37 + t182 * t36 + t30 * t302); t278 * (t180 * t239 + t182 * t243 - t225 * t242) * t324; 0.2e1 * t278 * (t225 ^ 2 * t232 + t180 ^ 2 + t182 ^ 2); t39 + m(7) * (t40 * t55 + t41 * t56) + (t239 * t275 + t243 * t274) * t238; m(7) * (t11 * t38 + t55 * t66 + t56 * t65) + (t316 * t5 + t320 * t4) * t238 + t246; m(7) * (t12 * t38 + t23 * t55 + t24 * t56) - t325; m(7) * (-t38 * t242 + (t239 * t56 + t243 * t55) * t238); m(7) * (t180 * t56 + t182 * t55 + t302 * t38); m(7) * (t38 ^ 2 + t55 ^ 2 + t56 ^ 2) + t325;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
