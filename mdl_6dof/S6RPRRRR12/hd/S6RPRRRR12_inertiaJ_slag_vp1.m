% Calculate joint inertia matrix for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_inertiaJ_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:35
% EndTime: 2019-03-09 07:51:02
% DurationCPUTime: 11.07s
% Computational Cost: add. (106758->618), mult. (303670->850), div. (0->0), fcn. (405083->18), ass. (0->280)
t262 = cos(pkin(14));
t264 = cos(pkin(6));
t272 = cos(qJ(1));
t259 = sin(pkin(14));
t269 = sin(qJ(1));
t324 = t269 * t259;
t250 = t262 * t272 - t264 * t324;
t268 = sin(qJ(3));
t271 = cos(qJ(3));
t323 = t269 * t262;
t249 = -t259 * t272 - t264 * t323;
t260 = sin(pkin(7));
t263 = cos(pkin(7));
t261 = sin(pkin(6));
t328 = t261 * t269;
t281 = t249 * t263 + t260 * t328;
t231 = -t250 * t268 + t271 * t281;
t241 = -t249 * t260 + t263 * t328;
t332 = sin(pkin(8));
t333 = cos(pkin(8));
t213 = -t231 * t332 + t241 * t333;
t325 = t264 * t272;
t248 = t259 * t325 + t323;
t247 = t262 * t325 - t324;
t327 = t261 * t272;
t282 = t247 * t263 - t260 * t327;
t229 = -t248 * t268 + t271 * t282;
t347 = t247 * t260 + t263 * t327;
t212 = -t229 * t332 - t333 * t347;
t346 = t264 ^ 2;
t345 = m(7) / 0.2e1;
t230 = t248 * t271 + t268 * t282;
t267 = sin(qJ(4));
t338 = cos(qJ(4));
t190 = t230 * t338 + (t229 * t333 - t332 * t347) * t267;
t266 = sin(qJ(5));
t337 = cos(qJ(5));
t171 = t190 * t266 - t212 * t337;
t232 = t250 * t271 + t268 * t281;
t192 = t232 * t338 + (t231 * t333 + t241 * t332) * t267;
t173 = t192 * t266 - t213 * t337;
t326 = t262 * t263;
t329 = t260 * t264;
t238 = t271 * t329 + (-t259 * t268 + t271 * t326) * t261;
t239 = t268 * t329 + (t259 * t271 + t268 * t326) * t261;
t246 = -t260 * t261 * t262 + t263 * t264;
t208 = t239 * t338 + (t238 * t333 + t246 * t332) * t267;
t223 = -t238 * t332 + t246 * t333;
t182 = t208 * t266 - t223 * t337;
t172 = t190 * t337 + t212 * t266;
t286 = t338 * t332;
t287 = t333 * t338;
t189 = -t229 * t287 + t230 * t267 + t286 * t347;
t265 = sin(qJ(6));
t270 = cos(qJ(6));
t133 = -t172 * t265 + t189 * t270;
t134 = t172 * t270 + t189 * t265;
t87 = Icges(7,5) * t134 + Icges(7,6) * t133 + Icges(7,3) * t171;
t89 = Icges(7,4) * t134 + Icges(7,2) * t133 + Icges(7,6) * t171;
t91 = Icges(7,1) * t134 + Icges(7,4) * t133 + Icges(7,5) * t171;
t28 = t133 * t89 + t134 * t91 + t171 * t87;
t174 = t192 * t337 + t213 * t266;
t191 = -t231 * t287 + t232 * t267 - t241 * t286;
t135 = -t174 * t265 + t191 * t270;
t136 = t174 * t270 + t191 * t265;
t88 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t173;
t90 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t173;
t92 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t173;
t29 = t133 * t90 + t134 * t92 + t171 * t88;
t183 = t208 * t337 + t223 * t266;
t207 = -t238 * t287 + t239 * t267 - t246 * t286;
t162 = -t183 * t265 + t207 * t270;
t163 = t183 * t270 + t207 * t265;
t106 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t182;
t107 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t182;
t108 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t182;
t41 = t106 * t171 + t107 * t133 + t108 * t134;
t1 = t171 * t28 + t173 * t29 + t182 * t41;
t344 = t1 / 0.2e1;
t30 = t135 * t89 + t136 * t91 + t173 * t87;
t31 = t135 * t90 + t136 * t92 + t173 * t88;
t42 = t106 * t173 + t107 * t135 + t108 * t136;
t2 = t171 * t30 + t173 * t31 + t182 * t42;
t343 = t2 / 0.2e1;
t34 = t162 * t89 + t163 * t91 + t182 * t87;
t35 = t162 * t90 + t163 * t92 + t182 * t88;
t49 = t182 * t106 + t162 * t107 + t163 * t108;
t45 = t49 * t182;
t9 = t34 * t171 + t35 * t173 + t45;
t342 = t9 / 0.2e1;
t341 = t171 / 0.2e1;
t340 = t173 / 0.2e1;
t339 = t182 / 0.2e1;
t336 = t172 * pkin(5);
t285 = -t134 * rSges(7,1) - t133 * rSges(7,2);
t93 = t171 * rSges(7,3) - t285;
t335 = t171 * pkin(13) + t336 + t93;
t94 = t136 * rSges(7,1) + t135 * rSges(7,2) + t173 * rSges(7,3);
t334 = t174 * pkin(5) + t173 * pkin(13) + t94;
t109 = rSges(7,1) * t163 + rSges(7,2) * t162 + rSges(7,3) * t182;
t322 = pkin(5) * t183 + pkin(13) * t182 + t109;
t116 = t172 * rSges(6,1) - t171 * rSges(6,2) + t189 * rSges(6,3);
t156 = t190 * pkin(4) + t189 * pkin(12);
t321 = -t116 - t156;
t117 = t174 * rSges(6,1) - t173 * rSges(6,2) + t191 * rSges(6,3);
t157 = t192 * pkin(4) + t191 * pkin(12);
t320 = -t117 - t157;
t130 = rSges(6,1) * t183 - rSges(6,2) * t182 + rSges(6,3) * t207;
t179 = pkin(4) * t208 + pkin(12) * t207;
t319 = -t130 - t179;
t193 = t230 * pkin(3) + pkin(11) * t212;
t180 = t241 * t193;
t318 = t241 * t156 + t180;
t194 = t232 * pkin(3) + pkin(11) * t213;
t181 = t246 * t194;
t317 = t246 * t157 + t181;
t214 = t239 * pkin(3) + pkin(11) * t223;
t203 = t347 * t214;
t316 = -t179 * t347 - t203;
t315 = t272 * pkin(1) + qJ(2) * t328;
t110 = Icges(6,5) * t172 - Icges(6,6) * t171 + Icges(6,3) * t189;
t112 = Icges(6,4) * t172 - Icges(6,2) * t171 + Icges(6,6) * t189;
t114 = Icges(6,1) * t172 - Icges(6,4) * t171 + Icges(6,5) * t189;
t51 = t110 * t189 - t112 * t171 + t114 * t172;
t111 = Icges(6,5) * t174 - Icges(6,6) * t173 + Icges(6,3) * t191;
t113 = Icges(6,4) * t174 - Icges(6,2) * t173 + Icges(6,6) * t191;
t115 = Icges(6,1) * t174 - Icges(6,4) * t173 + Icges(6,5) * t191;
t52 = t111 * t189 - t113 * t171 + t115 * t172;
t127 = Icges(6,5) * t183 - Icges(6,6) * t182 + Icges(6,3) * t207;
t128 = Icges(6,4) * t183 - Icges(6,2) * t182 + Icges(6,6) * t207;
t129 = Icges(6,1) * t183 - Icges(6,4) * t182 + Icges(6,5) * t207;
t65 = t127 * t189 - t128 * t171 + t129 * t172;
t13 = t189 * t51 + t191 * t52 + t207 * t65;
t3 = t189 * t28 + t191 * t29 + t207 * t41;
t313 = t3 / 0.2e1 + t13 / 0.2e1;
t53 = t110 * t191 - t112 * t173 + t114 * t174;
t54 = t111 * t191 - t113 * t173 + t115 * t174;
t66 = t127 * t191 - t128 * t173 + t129 * t174;
t14 = t189 * t53 + t191 * t54 + t207 * t66;
t4 = t189 * t30 + t191 * t31 + t207 * t42;
t312 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t212 * t51 + t213 * t52 + t223 * t65;
t5 = t212 * t28 + t213 * t29 + t223 * t41;
t311 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t212 * t53 + t213 * t54 + t223 * t66;
t6 = t212 * t30 + t213 * t31 + t223 * t42;
t310 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t241 * t52 + t246 * t65 - t347 * t51;
t7 = t241 * t29 + t246 * t41 - t28 * t347;
t309 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t241 * t54 + t246 * t66 - t347 * t53;
t8 = t241 * t31 + t246 * t42 - t30 * t347;
t308 = t8 / 0.2e1 + t18 / 0.2e1;
t46 = t49 * t207;
t10 = t34 * t189 + t35 * t191 + t46;
t56 = t110 * t207 - t112 * t182 + t114 * t183;
t57 = t111 * t207 - t113 * t182 + t115 * t183;
t72 = t207 * t127 - t182 * t128 + t183 * t129;
t69 = t72 * t207;
t19 = t56 * t189 + t57 * t191 + t69;
t307 = t10 / 0.2e1 + t19 / 0.2e1;
t47 = t49 * t223;
t11 = t34 * t212 + t35 * t213 + t47;
t70 = t72 * t223;
t20 = t56 * t212 + t57 * t213 + t70;
t306 = t11 / 0.2e1 + t20 / 0.2e1;
t48 = t49 * t246;
t12 = t35 * t241 - t34 * t347 + t48;
t71 = t72 * t246;
t21 = t57 * t241 - t347 * t56 + t71;
t305 = t12 / 0.2e1 + t21 / 0.2e1;
t304 = t34 / 0.2e1 + t41 / 0.2e1;
t303 = t35 / 0.2e1 + t42 / 0.2e1;
t302 = -t156 - t335;
t301 = -t157 - t334;
t298 = -t179 - t322;
t168 = Icges(5,5) * t208 - Icges(5,6) * t207 + Icges(5,3) * t223;
t169 = Icges(5,4) * t208 - Icges(5,2) * t207 + Icges(5,6) * t223;
t170 = Icges(5,1) * t208 - Icges(5,4) * t207 + Icges(5,5) * t223;
t97 = t223 * t168 - t207 * t169 + t208 * t170;
t144 = t192 * rSges(5,1) - t191 * rSges(5,2) + t213 * rSges(5,3);
t215 = Icges(4,5) * t239 + Icges(4,6) * t238 + Icges(4,3) * t246;
t216 = Icges(4,4) * t239 + Icges(4,2) * t238 + Icges(4,6) * t246;
t217 = Icges(4,1) * t239 + Icges(4,4) * t238 + Icges(4,5) * t246;
t297 = t246 * t215 + t238 * t216 + t239 * t217;
t202 = t232 * rSges(4,1) + t231 * rSges(4,2) + t241 * rSges(4,3);
t296 = -t269 * pkin(1) + qJ(2) * t327;
t284 = t56 / 0.2e1 + t65 / 0.2e1 + t304;
t283 = t57 / 0.2e1 + t66 / 0.2e1 + t303;
t201 = rSges(4,1) * t230 + rSges(4,2) * t229 - rSges(4,3) * t347;
t143 = t190 * rSges(5,1) - t189 * rSges(5,2) + t212 * rSges(5,3);
t280 = t250 * pkin(2) + pkin(10) * t241 + t315;
t279 = -t248 * pkin(2) + pkin(10) * t347 + t296;
t137 = Icges(5,5) * t190 - Icges(5,6) * t189 + Icges(5,3) * t212;
t139 = Icges(5,4) * t190 - Icges(5,2) * t189 + Icges(5,6) * t212;
t141 = Icges(5,1) * t190 - Icges(5,4) * t189 + Icges(5,5) * t212;
t78 = t137 * t223 - t139 * t207 + t141 * t208;
t85 = t168 * t212 - t169 * t189 + t170 * t190;
t278 = t78 / 0.2e1 + t85 / 0.2e1 + t284;
t138 = Icges(5,5) * t192 - Icges(5,6) * t191 + Icges(5,3) * t213;
t140 = Icges(5,4) * t192 - Icges(5,2) * t191 + Icges(5,6) * t213;
t142 = Icges(5,1) * t192 - Icges(5,4) * t191 + Icges(5,5) * t213;
t79 = t138 * t223 - t140 * t207 + t142 * t208;
t86 = t168 * t213 - t169 * t191 + t170 * t192;
t277 = t79 / 0.2e1 + t86 / 0.2e1 + t283;
t276 = t280 + t194;
t275 = -t193 + t279;
t274 = t157 + t276;
t273 = -t156 + t275;
t254 = rSges(2,1) * t272 - t269 * rSges(2,2);
t253 = -t269 * rSges(2,1) - rSges(2,2) * t272;
t228 = rSges(3,1) * t250 + rSges(3,2) * t249 + rSges(3,3) * t328 + t315;
t227 = -t248 * rSges(3,1) - t247 * rSges(3,2) + rSges(3,3) * t327 + t296;
t218 = rSges(4,1) * t239 + rSges(4,2) * t238 + rSges(4,3) * t246;
t200 = Icges(4,1) * t232 + Icges(4,4) * t231 + Icges(4,5) * t241;
t199 = Icges(4,1) * t230 + Icges(4,4) * t229 - Icges(4,5) * t347;
t198 = Icges(4,4) * t232 + Icges(4,2) * t231 + Icges(4,6) * t241;
t197 = Icges(4,4) * t230 + Icges(4,2) * t229 - Icges(4,6) * t347;
t196 = Icges(4,5) * t232 + Icges(4,6) * t231 + Icges(4,3) * t241;
t195 = Icges(4,5) * t230 + Icges(4,6) * t229 - Icges(4,3) * t347;
t178 = t280 + t202;
t177 = -t201 + t279;
t175 = rSges(5,1) * t208 - rSges(5,2) * t207 + rSges(5,3) * t223;
t161 = t202 * t246 - t218 * t241;
t160 = -t201 * t246 - t218 * t347;
t159 = t212 * t179;
t152 = t201 * t241 + t202 * t347;
t149 = t297 * t246;
t148 = t223 * t157;
t147 = t215 * t241 + t216 * t231 + t217 * t232;
t146 = -t215 * t347 + t216 * t229 + t217 * t230;
t145 = t213 * t156;
t124 = t196 * t246 + t198 * t238 + t200 * t239;
t123 = t195 * t246 + t197 * t238 + t199 * t239;
t119 = t276 + t144;
t118 = -t143 + t275;
t105 = t144 * t223 - t175 * t213;
t104 = -t143 * t223 + t175 * t212;
t100 = t144 * t246 + t181 + (-t175 - t214) * t241;
t99 = -t175 * t347 - t203 + (-t143 - t193) * t246;
t98 = t143 * t213 - t144 * t212;
t96 = t97 * t246;
t95 = t97 * t223;
t84 = t274 + t117;
t83 = -t116 + t273;
t82 = t143 * t241 + t180 - (-t144 - t194) * t347;
t81 = t117 * t207 - t130 * t191;
t80 = -t116 * t207 + t130 * t189;
t77 = t138 * t213 - t140 * t191 + t142 * t192;
t76 = t137 * t213 - t139 * t191 + t141 * t192;
t75 = t138 * t212 - t140 * t189 + t142 * t190;
t74 = t137 * t212 - t139 * t189 + t141 * t190;
t73 = t116 * t191 - t117 * t189;
t68 = t117 * t223 + t213 * t319 + t148;
t67 = t130 * t212 + t223 * t321 + t159;
t64 = t117 * t246 + (-t214 + t319) * t241 + t317;
t63 = -t130 * t347 + (-t193 + t321) * t246 + t316;
t62 = t274 + t334;
t61 = (-rSges(7,3) - pkin(13)) * t171 + t273 - t336 + t285;
t60 = -t109 * t173 + t182 * t94;
t59 = t109 * t171 - t182 * t93;
t58 = t116 * t213 + t212 * t320 + t145;
t55 = t116 * t241 - (-t194 + t320) * t347 + t318;
t50 = -t171 * t94 + t173 * t93;
t44 = -t191 * t322 + t207 * t334;
t43 = t189 * t322 - t207 * t335;
t40 = t334 * t246 + (-t214 + t298) * t241 + t317;
t39 = -t322 * t347 + (-t193 + t302) * t246 + t316;
t38 = t213 * t298 + t223 * t334 + t148;
t37 = t212 * t322 + t223 * t302 + t159;
t36 = -t189 * t334 + t191 * t335;
t33 = t79 * t241 - t347 * t78 + t96;
t32 = t78 * t212 + t79 * t213 + t95;
t27 = t335 * t241 - (-t194 + t301) * t347 + t318;
t26 = t212 * t301 + t213 * t335 + t145;
t25 = t241 * t77 + t246 * t86 - t347 * t76;
t24 = t241 * t75 + t246 * t85 - t347 * t74;
t23 = t212 * t76 + t213 * t77 + t223 * t86;
t22 = t212 * t74 + t213 * t75 + t223 * t85;
t101 = [(t83 ^ 2 + t84 ^ 2) * m(6) + (t61 ^ 2 + t62 ^ 2) * m(7) + ((Icges(3,2) * t262 ^ 2 + (Icges(3,1) * t259 + 0.2e1 * Icges(3,4) * t262) * t259) * t261 + 0.2e1 * t264 * (Icges(3,5) * t259 + Icges(3,6) * t262)) * t261 + t297 + (t118 ^ 2 + t119 ^ 2) * m(5) + (t177 ^ 2 + t178 ^ 2) * m(4) + Icges(2,3) + Icges(3,3) * t346 + t97 + t72 + t49 + m(3) * (t227 ^ 2 + t228 ^ 2) + m(2) * (t253 ^ 2 + t254 ^ 2); (-m(3) * t228 - m(4) * t178 - m(5) * t119 - m(6) * t84 - m(7) * t62) * t327 + (m(3) * t227 + m(4) * t177 + m(5) * t118 + m(6) * t83 + m(7) * t61) * t328; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + t345) * (t346 + (t269 ^ 2 + t272 ^ 2) * t261 ^ 2); (t39 * t61 + t40 * t62) * m(7) + t48 + (t63 * t83 + t64 * t84) * m(6) + t71 + (t100 * t119 + t118 * t99) * m(5) + t96 + (t160 * t177 + t161 * t178) * m(4) + t149 + (t147 / 0.2e1 + t124 / 0.2e1 + t277) * t241 - (t123 / 0.2e1 + t146 / 0.2e1 + t278) * t347; (m(4) * t152 + m(5) * t82 + m(6) * t55 + m(7) * t27) * t264 + ((-m(4) * t161 - m(5) * t100 - m(6) * t64 - m(7) * t40) * t272 + (m(4) * t160 + m(5) * t99 + m(6) * t63 + m(7) * t39) * t269) * t261; (t27 ^ 2 + t39 ^ 2 + t40 ^ 2) * m(7) + (t55 ^ 2 + t63 ^ 2 + t64 ^ 2) * m(6) + (t100 ^ 2 + t82 ^ 2 + t99 ^ 2) * m(5) + (t152 ^ 2 + t160 ^ 2 + t161 ^ 2) * m(4) + (t12 + t21 + t33 + t149) * t246 + (t18 + t25 + t8 + (t241 * t196 + t198 * t231 + t200 * t232) * t241 + (t124 + t147) * t246) * t241 - (t7 + t17 + t24 - (-t195 * t347 + t229 * t197 + t230 * t199) * t347 + (t123 + t146) * t246 + (t241 * t195 - t196 * t347 + t231 * t197 + t229 * t198 + t232 * t199 + t230 * t200) * t241) * t347; (t37 * t61 + t38 * t62) * m(7) + t47 + (t67 * t83 + t68 * t84) * m(6) + t70 + (t104 * t118 + t105 * t119) * m(5) + t95 + t277 * t213 + t278 * t212; (m(5) * t98 + m(6) * t58 + m(7) * t26) * t264 + ((-m(5) * t105 - m(6) * t68 - m(7) * t38) * t272 + (m(5) * t104 + m(6) * t67 + m(7) * t37) * t269) * t261; (t26 * t27 + t37 * t39 + t38 * t40) * m(7) + (t58 * t55 + t67 * t63 + t68 * t64) * m(6) + (t105 * t100 + t104 * t99 + t98 * t82) * m(5) + (t32 / 0.2e1 + t306) * t246 + (t23 / 0.2e1 + t310) * t241 - (t22 / 0.2e1 + t311) * t347 + (t33 / 0.2e1 + t305) * t223 + (t25 / 0.2e1 + t308) * t213 + (t24 / 0.2e1 + t309) * t212; (t26 ^ 2 + t37 ^ 2 + t38 ^ 2) * m(7) + (t58 ^ 2 + t67 ^ 2 + t68 ^ 2) * m(6) + (t104 ^ 2 + t105 ^ 2 + t98 ^ 2) * m(5) + (t11 + t20 + t32) * t223 + (t6 + t16 + t23) * t213 + (t5 + t15 + t22) * t212; (t43 * t61 + t44 * t62) * m(7) + t46 + (t80 * t83 + t81 * t84) * m(6) + t69 + t283 * t191 + t284 * t189; (m(6) * t73 + m(7) * t36) * t264 + ((-m(6) * t81 - m(7) * t44) * t272 + (m(6) * t80 + m(7) * t43) * t269) * t261; (t27 * t36 + t39 * t43 + t40 * t44) * m(7) + (t55 * t73 + t63 * t80 + t64 * t81) * m(6) + t307 * t246 + t312 * t241 - t313 * t347 + t305 * t207 + t308 * t191 + t309 * t189; (t26 * t36 + t37 * t43 + t38 * t44) * m(7) + (t58 * t73 + t67 * t80 + t68 * t81) * m(6) + t307 * t223 + t312 * t213 + t313 * t212 + t306 * t207 + t310 * t191 + t311 * t189; (t36 ^ 2 + t43 ^ 2 + t44 ^ 2) * m(7) + (t73 ^ 2 + t80 ^ 2 + t81 ^ 2) * m(6) + (t10 + t19) * t207 + (t4 + t14) * t191 + (t3 + t13) * t189; (t59 * t61 + t60 * t62) * m(7) + t45 + t303 * t173 + t304 * t171; 0.2e1 * (t50 * t264 + (t269 * t59 - t272 * t60) * t261) * t345; (t27 * t50 + t39 * t59 + t40 * t60) * m(7) + t246 * t342 - t347 * t344 + t7 * t341 + t8 * t340 + t241 * t343 + t12 * t339; (t26 * t50 + t37 * t59 + t38 * t60) * m(7) + t11 * t339 + t223 * t342 + t212 * t344 + t6 * t340 + t5 * t341 + t213 * t343; t4 * t340 + t191 * t343 + t189 * t344 + (t36 * t50 + t43 * t59 + t44 * t60) * m(7) + t207 * t342 + t3 * t341 + t10 * t339; t173 * t2 + t171 * t1 + t182 * t9 + (t50 ^ 2 + t59 ^ 2 + t60 ^ 2) * m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t101(1) t101(2) t101(4) t101(7) t101(11) t101(16); t101(2) t101(3) t101(5) t101(8) t101(12) t101(17); t101(4) t101(5) t101(6) t101(9) t101(13) t101(18); t101(7) t101(8) t101(9) t101(10) t101(14) t101(19); t101(11) t101(12) t101(13) t101(14) t101(15) t101(20); t101(16) t101(17) t101(18) t101(19) t101(20) t101(21);];
Mq  = res;
