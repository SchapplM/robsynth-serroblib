% Calculate joint inertia matrix for
% S6RRPRRP13
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP13_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP13_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:42
% EndTime: 2019-03-09 12:55:55
% DurationCPUTime: 5.70s
% Computational Cost: add. (16127->616), mult. (40708->843), div. (0->0), fcn. (51500->10), ass. (0->272)
t337 = rSges(7,3) + qJ(6) + pkin(10);
t260 = cos(pkin(6));
t267 = cos(qJ(2));
t268 = cos(qJ(1));
t315 = t267 * t268;
t264 = sin(qJ(2));
t265 = sin(qJ(1));
t318 = t264 * t265;
t241 = -t260 * t315 + t318;
t263 = sin(qJ(4));
t259 = sin(pkin(6));
t329 = cos(qJ(4));
t288 = t259 * t329;
t214 = t241 * t263 - t268 * t288;
t316 = t265 * t267;
t317 = t264 * t268;
t242 = t260 * t317 + t316;
t262 = sin(qJ(5));
t266 = cos(qJ(5));
t167 = -t214 * t262 + t242 * t266;
t324 = t242 * t262;
t168 = t214 * t266 + t324;
t319 = t259 * t268;
t212 = t241 * t329 + t263 * t319;
t338 = -rSges(7,1) * t168 - rSges(7,2) * t167 + t337 * t212;
t221 = Icges(3,3) * t260 + (Icges(3,5) * t264 + Icges(3,6) * t267) * t259;
t222 = Icges(3,6) * t260 + (Icges(3,4) * t264 + Icges(3,2) * t267) * t259;
t223 = Icges(3,5) * t260 + (Icges(3,1) * t264 + Icges(3,4) * t267) * t259;
t224 = Icges(4,5) * t260 + (-Icges(4,6) * t264 - Icges(4,3) * t267) * t259;
t225 = Icges(4,4) * t260 + (-Icges(4,2) * t264 - Icges(4,6) * t267) * t259;
t226 = Icges(4,1) * t260 + (-Icges(4,4) * t264 - Icges(4,5) * t267) * t259;
t320 = t259 * t267;
t322 = t259 * t264;
t336 = (-t224 * t267 - t225 * t264) * t259 + t222 * t320 + t223 * t322 + (t226 + t221) * t260;
t243 = t260 * t316 + t317;
t211 = t243 * t263 + t265 * t288;
t244 = -t260 * t318 + t315;
t165 = -t211 * t262 + t244 * t266;
t323 = t244 * t262;
t166 = t211 * t266 + t323;
t321 = t259 * t265;
t210 = -t243 * t329 + t263 * t321;
t257 = pkin(5) * t266 + pkin(4);
t335 = t166 * rSges(7,1) + t165 * rSges(7,2) + pkin(5) * t323 + t337 * t210 + t211 * t257;
t175 = -Icges(4,5) * t319 - Icges(4,6) * t242 + Icges(4,3) * t241;
t182 = Icges(3,4) * t242 - Icges(3,2) * t241 - Icges(3,6) * t319;
t334 = t175 - t182;
t177 = -Icges(4,4) * t319 - Icges(4,2) * t242 + Icges(4,6) * t241;
t184 = Icges(3,1) * t242 - Icges(3,4) * t241 - Icges(3,5) * t319;
t333 = t177 - t184;
t174 = Icges(4,5) * t321 - Icges(4,6) * t244 + Icges(4,3) * t243;
t183 = Icges(3,4) * t244 - Icges(3,2) * t243 + Icges(3,6) * t321;
t332 = -t183 + t174;
t176 = Icges(4,4) * t321 - Icges(4,2) * t244 + Icges(4,6) * t243;
t185 = Icges(3,1) * t244 - Icges(3,4) * t243 + Icges(3,5) * t321;
t331 = t185 - t176;
t330 = t259 ^ 2;
t328 = t242 * pkin(2);
t327 = -pkin(4) + t257;
t154 = t211 * pkin(4) + pkin(10) * t210;
t326 = -t154 + t335;
t206 = t212 * pkin(10);
t325 = pkin(5) * t324 + t214 * t327 + t206 - t338;
t110 = t166 * rSges(6,1) + t165 * rSges(6,2) + t210 * rSges(6,3);
t314 = -t110 - t154;
t112 = rSges(6,1) * t168 + rSges(6,2) * t167 - rSges(6,3) * t212;
t155 = pkin(4) * t214 - t206;
t313 = -t112 - t155;
t240 = t260 * t329 - t263 * t320;
t207 = -t240 * t262 + t266 * t322;
t294 = t262 * t322;
t208 = t240 * t266 + t294;
t239 = t260 * t263 + t267 * t288;
t312 = rSges(7,1) * t208 + rSges(7,2) * t207 + pkin(5) * t294 + t327 * t240 + (-pkin(10) + t337) * t239;
t311 = t336 * t260;
t179 = -Icges(4,1) * t319 - Icges(4,4) * t242 + Icges(4,5) * t241;
t180 = Icges(3,5) * t242 - Icges(3,6) * t241 - Icges(3,3) * t319;
t310 = -t180 - t179;
t178 = Icges(4,1) * t321 - Icges(4,4) * t244 + Icges(4,5) * t243;
t181 = Icges(3,5) * t244 - Icges(3,6) * t243 + Icges(3,3) * t321;
t309 = t181 + t178;
t231 = t241 * qJ(3);
t195 = t231 + t328;
t196 = t244 * pkin(2) + qJ(3) * t243;
t308 = t195 * t321 + t196 * t319;
t193 = t260 * t196;
t217 = pkin(3) * t321 + pkin(9) * t244;
t307 = t260 * t217 + t193;
t218 = -pkin(3) * t319 + t242 * pkin(9);
t306 = -t195 - t218;
t245 = (pkin(2) * t264 - qJ(3) * t267) * t259;
t305 = -pkin(3) * t260 - pkin(9) * t322 - t245;
t304 = t268 * pkin(1) + pkin(8) * t321;
t101 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t210;
t105 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t210;
t97 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t210;
t33 = t101 * t165 + t105 * t166 + t210 * t97;
t102 = Icges(7,4) * t168 + Icges(7,2) * t167 - Icges(7,6) * t212;
t106 = Icges(7,1) * t168 + Icges(7,4) * t167 - Icges(7,5) * t212;
t98 = Icges(7,5) * t168 + Icges(7,6) * t167 - Icges(7,3) * t212;
t34 = t102 * t165 + t106 * t166 + t210 * t98;
t131 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t239;
t133 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t239;
t135 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t239;
t50 = t131 * t210 + t133 * t165 + t135 * t166;
t1 = t210 * t33 - t212 * t34 + t239 * t50;
t103 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t210;
t107 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t210;
t99 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t210;
t35 = t103 * t165 + t107 * t166 + t210 * t99;
t100 = Icges(6,5) * t168 + Icges(6,6) * t167 - Icges(6,3) * t212;
t104 = Icges(6,4) * t168 + Icges(6,2) * t167 - Icges(6,6) * t212;
t108 = Icges(6,1) * t168 + Icges(6,4) * t167 - Icges(6,5) * t212;
t36 = t100 * t210 + t104 * t165 + t108 * t166;
t132 = Icges(6,5) * t208 + Icges(6,6) * t207 + Icges(6,3) * t239;
t134 = Icges(6,4) * t208 + Icges(6,2) * t207 + Icges(6,6) * t239;
t136 = Icges(6,1) * t208 + Icges(6,4) * t207 + Icges(6,5) * t239;
t51 = t132 * t210 + t134 * t165 + t136 * t166;
t2 = t210 * t35 - t212 * t36 + t239 * t51;
t303 = t1 / 0.2e1 + t2 / 0.2e1;
t37 = t101 * t167 + t105 * t168 - t212 * t97;
t38 = t102 * t167 + t106 * t168 - t212 * t98;
t52 = -t131 * t212 + t133 * t167 + t135 * t168;
t3 = t210 * t37 - t212 * t38 + t239 * t52;
t39 = t103 * t167 + t107 * t168 - t212 * t99;
t40 = -t100 * t212 + t104 * t167 + t108 * t168;
t53 = -t132 * t212 + t134 * t167 + t136 * t168;
t4 = t210 * t39 - t212 * t40 + t239 * t53;
t302 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t242 * t34 + t244 * t33 + t322 * t50;
t6 = t242 * t36 + t244 * t35 + t322 * t51;
t301 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t242 * t38 + t244 * t37 + t322 * t52;
t8 = t242 * t40 + t244 * t39 + t322 * t53;
t300 = -t8 / 0.2e1 - t7 / 0.2e1;
t10 = t51 * t260 + (t265 * t35 - t268 * t36) * t259;
t9 = t50 * t260 + (t265 * t33 - t268 * t34) * t259;
t299 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t52 * t260 + (t265 * t37 - t268 * t38) * t259;
t12 = t53 * t260 + (t265 * t39 - t268 * t40) * t259;
t298 = -t11 / 0.2e1 - t12 / 0.2e1;
t43 = t101 * t207 + t105 * t208 + t239 * t97;
t44 = t102 * t207 + t106 * t208 + t239 * t98;
t63 = t239 * t131 + t207 * t133 + t208 * t135;
t57 = t63 * t239;
t13 = t43 * t210 - t44 * t212 + t57;
t45 = t103 * t207 + t107 * t208 + t239 * t99;
t46 = t100 * t239 + t104 * t207 + t108 * t208;
t64 = t239 * t132 + t207 * t134 + t208 * t136;
t58 = t64 * t239;
t14 = t45 * t210 - t46 * t212 + t58;
t297 = t14 / 0.2e1 + t13 / 0.2e1;
t59 = t63 * t322;
t15 = t44 * t242 + t43 * t244 + t59;
t60 = t64 * t322;
t16 = t46 * t242 + t45 * t244 + t60;
t296 = t15 / 0.2e1 + t16 / 0.2e1;
t61 = t63 * t260;
t17 = t61 + (t43 * t265 - t44 * t268) * t259;
t62 = t64 * t260;
t18 = t62 + (t45 * t265 - t46 * t268) * t259;
t295 = t17 / 0.2e1 + t18 / 0.2e1;
t293 = t260 * t154 + t307;
t292 = -t155 + t306;
t170 = Icges(5,5) * t240 - Icges(5,6) * t239 + Icges(5,3) * t322;
t171 = Icges(5,4) * t240 - Icges(5,2) * t239 + Icges(5,6) * t322;
t172 = Icges(5,1) * t240 - Icges(5,4) * t239 + Icges(5,5) * t322;
t86 = t170 * t322 - t239 * t171 + t240 * t172;
t194 = pkin(4) * t240 + pkin(10) * t239;
t291 = -t194 + t305;
t145 = t211 * rSges(5,1) - t210 * rSges(5,2) + t244 * rSges(5,3);
t189 = t244 * rSges(3,1) - t243 * rSges(3,2) + rSges(3,3) * t321;
t186 = rSges(4,1) * t321 - t244 * rSges(4,2) + t243 * rSges(4,3);
t287 = -t265 * pkin(1) + pkin(8) * t319;
t286 = t259 * (-rSges(4,1) * t260 - (-rSges(4,2) * t264 - rSges(4,3) * t267) * t259 - t245);
t285 = t217 * t319 + t218 * t321 + t308;
t284 = -t231 + t287;
t173 = rSges(5,1) * t240 - rSges(5,2) * t239 + rSges(5,3) * t322;
t283 = t259 * (-t173 + t305);
t282 = -rSges(5,1) * t214 - rSges(5,2) * t212;
t138 = rSges(6,1) * t208 + rSges(6,2) * t207 + rSges(6,3) * t239;
t280 = t259 * (-t138 + t291);
t279 = t196 + t304;
t278 = t45 / 0.2e1 + t43 / 0.2e1 + t51 / 0.2e1 + t50 / 0.2e1;
t277 = -t46 / 0.2e1 - t44 / 0.2e1 - t53 / 0.2e1 - t52 / 0.2e1;
t276 = rSges(4,1) * t319 - t241 * rSges(4,3);
t275 = t154 * t319 + t155 * t321 + t285;
t274 = t284 - t218;
t273 = t259 * (t291 - t312);
t188 = t242 * rSges(3,1) - t241 * rSges(3,2) - rSges(3,3) * t319;
t139 = Icges(5,5) * t211 - Icges(5,6) * t210 + Icges(5,3) * t244;
t141 = Icges(5,4) * t211 - Icges(5,2) * t210 + Icges(5,6) * t244;
t143 = Icges(5,1) * t211 - Icges(5,4) * t210 + Icges(5,5) * t244;
t72 = t139 * t322 - t141 * t239 + t143 * t240;
t80 = t170 * t244 - t171 * t210 + t172 * t211;
t271 = t72 / 0.2e1 + t80 / 0.2e1 + t278;
t140 = Icges(5,5) * t214 + Icges(5,6) * t212 + Icges(5,3) * t242;
t142 = Icges(5,4) * t214 + Icges(5,2) * t212 + Icges(5,6) * t242;
t144 = Icges(5,1) * t214 + Icges(5,4) * t212 + Icges(5,5) * t242;
t73 = t140 * t322 - t142 * t239 + t144 * t240;
t81 = t170 * t242 + t171 * t212 + t172 * t214;
t270 = t81 / 0.2e1 + t73 / 0.2e1 - t277;
t269 = t217 + t279;
t248 = rSges(2,1) * t268 - t265 * rSges(2,2);
t247 = -t265 * rSges(2,1) - rSges(2,2) * t268;
t227 = rSges(3,3) * t260 + (rSges(3,1) * t264 + rSges(3,2) * t267) * t259;
t187 = -t242 * rSges(4,2) - t276;
t164 = t242 * t194;
t159 = t189 + t304;
t158 = -t188 + t287;
t151 = t154 * t322;
t149 = -t260 * t188 - t227 * t319;
t148 = t189 * t260 - t227 * t321;
t147 = t244 * t155;
t146 = rSges(5,3) * t242 - t282;
t127 = t279 + t186;
t126 = (rSges(4,2) - pkin(2)) * t242 + t276 + t284;
t123 = (t188 * t265 + t189 * t268) * t259;
t122 = t221 * t321 - t222 * t243 + t223 * t244;
t121 = -t221 * t319 - t241 * t222 + t242 * t223;
t120 = t241 * t224 - t242 * t225 - t226 * t319;
t119 = t224 * t243 - t225 * t244 + t226 * t321;
t114 = (-t187 - t195) * t260 + t268 * t286;
t113 = t186 * t260 + t265 * t286 + t193;
t96 = t145 * t322 - t173 * t244;
t95 = -t146 * t322 + t173 * t242;
t94 = t179 * t260 + (-t175 * t267 - t177 * t264) * t259;
t93 = t178 * t260 + (-t174 * t267 - t176 * t264) * t259;
t92 = t181 * t260 + (t183 * t267 + t185 * t264) * t259;
t91 = t180 * t260 + (t182 * t267 + t184 * t264) * t259;
t88 = t269 + t145;
t87 = (-rSges(5,3) - pkin(2)) * t242 + t274 + t282;
t85 = t86 * t260;
t84 = t86 * t322;
t83 = (t186 * t268 + t187 * t265) * t259 + t308;
t82 = -t145 * t242 + t146 * t244;
t79 = (-t146 + t306) * t260 + t268 * t283;
t78 = t145 * t260 + t265 * t283 + t307;
t77 = -t112 * t239 - t138 * t212;
t76 = t110 * t239 - t138 * t210;
t75 = t269 - t314;
t74 = t274 + t313 - t328;
t71 = t269 + t335;
t70 = -t214 * t257 + (-pkin(5) * t262 - pkin(2)) * t242 + t274 + t338;
t69 = t140 * t242 + t142 * t212 + t144 * t214;
t68 = t139 * t242 + t141 * t212 + t143 * t214;
t67 = t140 * t244 - t142 * t210 + t144 * t211;
t66 = t139 * t244 - t141 * t210 + t143 * t211;
t65 = (t145 * t268 + t146 * t265) * t259 + t285;
t56 = t110 * t212 + t112 * t210;
t55 = t110 * t322 + t151 + (-t138 - t194) * t244;
t54 = t138 * t242 + t313 * t322 + t164;
t49 = (-t112 + t292) * t260 + t268 * t280;
t48 = t110 * t260 + t265 * t280 + t293;
t47 = t112 * t244 + t242 * t314 + t147;
t42 = -t212 * t312 - t239 * t325;
t41 = -t210 * t312 + t239 * t326;
t32 = (t110 * t268 + t112 * t265) * t259 + t275;
t31 = t151 + t326 * t322 + (-t194 - t312) * t244;
t30 = t164 + t312 * t242 + (-t155 - t325) * t322;
t29 = (t292 - t325) * t260 + t268 * t273;
t28 = t260 * t326 + t265 * t273 + t293;
t27 = t210 * t325 + t212 * t326;
t26 = t85 + (t72 * t265 - t73 * t268) * t259;
t25 = t73 * t242 + t72 * t244 + t84;
t24 = t147 + t325 * t244 + (-t154 - t326) * t242;
t23 = t81 * t260 + (t265 * t68 - t268 * t69) * t259;
t22 = t80 * t260 + (t265 * t66 - t268 * t67) * t259;
t21 = t242 * t69 + t244 * t68 + t322 * t81;
t20 = t242 * t67 + t244 * t66 + t322 * t80;
t19 = (t265 * t325 + t268 * t326) * t259 + t275;
t89 = [Icges(2,3) + m(7) * (t70 ^ 2 + t71 ^ 2) + m(6) * (t74 ^ 2 + t75 ^ 2) + m(5) * (t87 ^ 2 + t88 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t158 ^ 2 + t159 ^ 2) + m(2) * (t247 ^ 2 + t248 ^ 2) + t86 + t64 + t63 + t336; t62 + t61 + t85 + m(7) * (t28 * t71 + t29 * t70) + m(6) * (t48 * t75 + t49 * t74) + m(5) * (t78 * t88 + t79 * t87) + m(4) * (t113 * t127 + t114 * t126) + m(3) * (t148 * t159 + t149 * t158) + ((-t91 / 0.2e1 - t94 / 0.2e1 - t120 / 0.2e1 - t121 / 0.2e1 - t270) * t268 + (t92 / 0.2e1 + t93 / 0.2e1 + t119 / 0.2e1 + t122 / 0.2e1 + t271) * t265) * t259 + t311; (t17 + t18 + t26 + t311) * t260 + m(7) * (t19 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t32 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t65 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t113 ^ 2 + t114 ^ 2 + t83 ^ 2) + m(3) * (t123 ^ 2 + t148 ^ 2 + t149 ^ 2) + ((-t11 - t12 - t23 + ((t241 * t334 - t242 * t333) * t259 + t310 * t330 * t268) * t268 + (-t120 - t121 - t91 - t94) * t260) * t268 + (t9 + t10 + t22 + ((t243 * t332 + t244 * t331) * t259 + t309 * t330 * t265) * t265 + (t122 + t119 + t92 + t93) * t260 + ((t265 * t310 + t268 * t309) * t259 + t333 * t244 - t334 * t243 - t331 * t242 - t332 * t241) * t319) * t265) * t259; m(7) * (t241 * t71 + t243 * t70) + m(6) * (t241 * t75 + t243 * t74) + m(5) * (t241 * t88 + t243 * t87) + m(4) * (t126 * t243 + t127 * t241); m(7) * (-t19 * t320 + t241 * t28 + t243 * t29) + m(6) * (t241 * t48 + t243 * t49 - t32 * t320) + m(5) * (t241 * t78 + t243 * t79 - t320 * t65) + m(4) * (t113 * t241 + t114 * t243 - t320 * t83); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t267 ^ 2 * t330 + t241 ^ 2 + t243 ^ 2); t60 + t59 + t84 + m(7) * (t30 * t70 + t31 * t71) + m(6) * (t54 * t74 + t55 * t75) + m(5) * (t87 * t95 + t88 * t96) + t271 * t244 + t270 * t242; (t25 / 0.2e1 + t296) * t260 + (t22 / 0.2e1 + t299) * t244 + (t23 / 0.2e1 - t298) * t242 + m(7) * (t19 * t24 + t28 * t31 + t29 * t30) + m(6) * (t32 * t47 + t48 * t55 + t49 * t54) + m(5) * (t65 * t82 + t78 * t96 + t79 * t95) + ((-t21 / 0.2e1 + t300) * t268 + (t20 / 0.2e1 + t301) * t265 + (t26 / 0.2e1 + t295) * t264) * t259; m(5) * (t241 * t96 + t243 * t95 - t320 * t82) + m(6) * (t241 * t55 + t243 * t54 - t320 * t47) + m(7) * (-t24 * t320 + t241 * t31 + t243 * t30); (t15 + t16 + t25) * t322 + (t6 + t5 + t20) * t244 + (t7 + t8 + t21) * t242 + m(7) * (t24 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t82 ^ 2 + t95 ^ 2 + t96 ^ 2); t58 + t57 + m(7) * (t41 * t71 + t42 * t70) + m(6) * (t74 * t77 + t75 * t76) + t277 * t212 + t278 * t210; t297 * t260 + t295 * t239 + t298 * t212 + t299 * t210 + m(7) * (t19 * t27 + t28 * t41 + t29 * t42) + m(6) * (t56 * t32 + t48 * t76 + t49 * t77) + (t265 * t303 - t268 * t302) * t259; m(6) * (t241 * t76 + t243 * t77 - t320 * t56) + m(7) * (t241 * t41 + t243 * t42 - t27 * t320); t297 * t322 + t303 * t244 + t302 * t242 + t296 * t239 + t300 * t212 + t301 * t210 + m(7) * (t24 * t27 + t30 * t42 + t31 * t41) + m(6) * (t56 * t47 + t54 * t77 + t55 * t76); (t14 + t13) * t239 + (-t4 - t3) * t212 + (t2 + t1) * t210 + m(7) * (t27 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t56 ^ 2 + t76 ^ 2 + t77 ^ 2); m(7) * (t210 * t70 - t212 * t71); m(7) * (t19 * t239 + t210 * t29 - t212 * t28); m(7) * (t210 * t243 - t212 * t241 - t239 * t320); m(7) * (t210 * t30 - t212 * t31 + t239 * t24); m(7) * (t210 * t42 - t212 * t41 + t239 * t27); m(7) * (t210 ^ 2 + t212 ^ 2 + t239 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t89(1) t89(2) t89(4) t89(7) t89(11) t89(16); t89(2) t89(3) t89(5) t89(8) t89(12) t89(17); t89(4) t89(5) t89(6) t89(9) t89(13) t89(18); t89(7) t89(8) t89(9) t89(10) t89(14) t89(19); t89(11) t89(12) t89(13) t89(14) t89(15) t89(20); t89(16) t89(17) t89(18) t89(19) t89(20) t89(21);];
Mq  = res;
