% Calculate joint inertia matrix for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:36
% EndTime: 2019-03-09 08:54:53
% DurationCPUTime: 7.52s
% Computational Cost: add. (22565->583), mult. (45684->828), div. (0->0), fcn. (59352->14), ass. (0->268)
t256 = sin(pkin(6));
t314 = sin(pkin(11));
t315 = cos(pkin(11));
t322 = sin(qJ(2));
t323 = cos(qJ(2));
t265 = t314 * t323 + t315 * t322;
t227 = t265 * t256;
t255 = sin(pkin(12));
t257 = cos(pkin(12));
t258 = cos(pkin(6));
t216 = -t227 * t255 + t257 * t258;
t311 = t255 * t258;
t217 = t227 * t257 + t311;
t238 = -t322 * t314 + t323 * t315;
t226 = t238 * t256;
t151 = Icges(5,5) * t217 + Icges(5,6) * t216 - Icges(5,3) * t226;
t152 = Icges(5,4) * t217 + Icges(5,2) * t216 - Icges(5,6) * t226;
t153 = Icges(5,1) * t217 + Icges(5,4) * t216 - Icges(5,5) * t226;
t184 = Icges(4,5) * t227 + Icges(4,6) * t226 + Icges(4,3) * t258;
t185 = Icges(4,4) * t227 + Icges(4,2) * t226 + Icges(4,6) * t258;
t186 = Icges(4,1) * t227 + Icges(4,4) * t226 + Icges(4,5) * t258;
t222 = Icges(3,3) * t258 + (Icges(3,5) * t322 + Icges(3,6) * t323) * t256;
t223 = Icges(3,6) * t258 + (Icges(3,4) * t322 + Icges(3,2) * t323) * t256;
t224 = Icges(3,5) * t258 + (Icges(3,1) * t322 + Icges(3,4) * t323) * t256;
t291 = t256 * t322;
t337 = t256 * t323 * t223 + t216 * t152 + t217 * t153 + t227 * t186 + t224 * t291 + (t184 + t222) * t258 + (-t151 + t185) * t226;
t329 = m(7) / 0.2e1;
t330 = m(6) / 0.2e1;
t331 = m(5) / 0.2e1;
t282 = t330 + t329 + t331;
t336 = 0.2e1 * t282;
t261 = sin(qJ(1));
t263 = cos(qJ(1));
t264 = t258 * t238;
t208 = -t261 * t265 + t263 * t264;
t228 = t265 * t258;
t210 = t228 * t263 + t261 * t238;
t309 = t256 * t263;
t140 = Icges(4,5) * t210 + Icges(4,6) * t208 - Icges(4,3) * t309;
t288 = t263 * t323;
t289 = t261 * t322;
t233 = t258 * t288 - t289;
t287 = t263 * t322;
t290 = t261 * t323;
t234 = t258 * t287 + t290;
t191 = Icges(3,5) * t234 + Icges(3,6) * t233 - Icges(3,3) * t309;
t335 = -t140 - t191;
t211 = -t261 * t264 - t263 * t265;
t213 = -t261 * t228 + t238 * t263;
t310 = t256 * t261;
t141 = Icges(4,5) * t213 + Icges(4,6) * t211 + Icges(4,3) * t310;
t235 = -t258 * t290 - t287;
t236 = -t258 * t289 + t288;
t192 = Icges(3,5) * t236 + Icges(3,6) * t235 + Icges(3,3) * t310;
t334 = t192 + t141;
t333 = t337 * t258;
t332 = m(4) / 0.2e1;
t303 = pkin(12) + qJ(5);
t253 = sin(t303);
t286 = cos(t303);
t274 = t256 * t286;
t169 = t210 * t253 + t263 * t274;
t171 = t213 * t253 - t261 * t274;
t170 = t210 * t286 - t253 * t309;
t260 = sin(qJ(6));
t262 = cos(qJ(6));
t127 = -t170 * t260 - t208 * t262;
t128 = t170 * t262 - t208 * t260;
t69 = Icges(7,5) * t128 + Icges(7,6) * t127 + Icges(7,3) * t169;
t71 = Icges(7,4) * t128 + Icges(7,2) * t127 + Icges(7,6) * t169;
t73 = Icges(7,1) * t128 + Icges(7,4) * t127 + Icges(7,5) * t169;
t19 = t127 * t71 + t128 * t73 + t169 * t69;
t172 = t213 * t286 + t253 * t310;
t129 = -t172 * t260 - t211 * t262;
t130 = t172 * t262 - t211 * t260;
t70 = Icges(7,5) * t130 + Icges(7,6) * t129 + Icges(7,3) * t171;
t72 = Icges(7,4) * t130 + Icges(7,2) * t129 + Icges(7,6) * t171;
t74 = Icges(7,1) * t130 + Icges(7,4) * t129 + Icges(7,5) * t171;
t20 = t127 * t72 + t128 * t74 + t169 * t70;
t214 = t227 * t253 - t258 * t286;
t215 = t227 * t286 + t258 * t253;
t163 = -t215 * t260 - t226 * t262;
t164 = t215 * t262 - t226 * t260;
t92 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t214;
t93 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t214;
t94 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t214;
t29 = t127 * t93 + t128 * t94 + t169 * t92;
t1 = t169 * t19 + t171 * t20 + t214 * t29;
t328 = -t1 / 0.2e1;
t24 = t163 * t71 + t164 * t73 + t214 * t69;
t25 = t163 * t72 + t164 * t74 + t214 * t70;
t37 = t163 * t93 + t164 * t94 + t214 * t92;
t33 = t37 * t214;
t7 = t24 * t169 + t25 * t171 + t33;
t327 = t7 / 0.2e1;
t326 = t169 / 0.2e1;
t325 = t171 / 0.2e1;
t324 = t214 / 0.2e1;
t321 = pkin(1) * t263;
t320 = t170 * pkin(5);
t251 = pkin(4) * t257 + pkin(3);
t319 = -pkin(3) + t251;
t275 = -t128 * rSges(7,1) - t127 * rSges(7,2);
t75 = t169 * rSges(7,3) - t275;
t318 = t169 * pkin(10) + t320 + t75;
t76 = t130 * rSges(7,1) + t129 * rSges(7,2) + t171 * rSges(7,3);
t317 = t172 * pkin(5) + pkin(10) * t171 + t76;
t95 = rSges(7,1) * t164 + rSges(7,2) * t163 + rSges(7,3) * t214;
t316 = pkin(5) * t215 + pkin(10) * t214 + t95;
t259 = -pkin(9) - qJ(4);
t313 = t208 * t259;
t229 = t258 * t322 * pkin(2) + (-pkin(8) - qJ(3)) * t256;
t312 = t229 * t263;
t252 = pkin(2) * t323 + pkin(1);
t308 = t261 * t252;
t246 = t263 * t252;
t201 = -t321 + t246 + (-t256 * pkin(8) - t229) * t261;
t189 = t258 * t201;
t284 = -t213 * pkin(3) + qJ(4) * t211;
t307 = -t258 * t284 + t189;
t202 = t208 * qJ(4);
t155 = t210 * pkin(3) - t202;
t250 = pkin(8) * t309;
t200 = t312 + t250 + (-pkin(1) + t252) * t261;
t306 = -t155 - t200;
t305 = t200 * t310 + t201 * t309;
t239 = pkin(2) * t291 + t258 * qJ(3);
t304 = -pkin(3) * t227 + qJ(4) * t226 - t239;
t30 = t129 * t93 + t130 * t94 + t171 * t92;
t302 = t25 / 0.2e1 + t30 / 0.2e1;
t301 = t29 / 0.2e1 + t24 / 0.2e1;
t298 = t255 * t310;
t293 = pkin(4) * t298 + t211 * t259 + t213 * t251;
t91 = t284 + t293;
t300 = t258 * t91 + t307;
t297 = t255 * t309;
t242 = pkin(4) * t297;
t90 = t210 * t319 + t202 - t242 + t313;
t299 = -t90 + t306;
t136 = Icges(6,5) * t215 - Icges(6,6) * t214 - Icges(6,3) * t226;
t137 = Icges(6,4) * t215 - Icges(6,2) * t214 - Icges(6,6) * t226;
t138 = Icges(6,1) * t215 - Icges(6,4) * t214 - Icges(6,5) * t226;
t60 = -t226 * t136 - t214 * t137 + t215 * t138;
t295 = -pkin(4) * t311 - t319 * t227 + (-qJ(4) - t259) * t226 + t304;
t103 = t172 * rSges(6,1) - t171 * rSges(6,2) - t211 * rSges(6,3);
t179 = -t213 * t255 + t257 * t310;
t180 = t213 * t257 + t298;
t111 = t180 * rSges(5,1) + t179 * rSges(5,2) - t211 * rSges(5,3);
t147 = t213 * rSges(4,1) + t211 * rSges(4,2) + rSges(4,3) * t310;
t198 = t236 * rSges(3,1) + t235 * rSges(3,2) + rSges(3,3) * t310;
t285 = t256 * (-rSges(4,1) * t227 - rSges(4,2) * t226 - rSges(4,3) * t258 - t239);
t283 = -t261 * t229 + t246;
t281 = t155 * t310 - t284 * t309 + t305;
t280 = t256 * (-rSges(5,1) * t217 - rSges(5,2) * t216 + rSges(5,3) * t226 + t304);
t277 = -t210 * rSges(4,1) - t208 * rSges(4,2);
t276 = -t170 * rSges(6,1) + t169 * rSges(6,2);
t273 = -t308 - t312;
t139 = rSges(6,1) * t215 - rSges(6,2) * t214 - rSges(6,3) * t226;
t272 = t256 * (-t139 + t295);
t101 = Icges(6,1) * t172 - Icges(6,4) * t171 - Icges(6,5) * t211;
t97 = Icges(6,5) * t172 - Icges(6,6) * t171 - Icges(6,3) * t211;
t99 = Icges(6,4) * t172 - Icges(6,2) * t171 - Icges(6,6) * t211;
t44 = t101 * t215 - t214 * t99 - t226 * t97;
t54 = -t136 * t211 - t137 * t171 + t138 * t172;
t271 = -t54 / 0.2e1 - t44 / 0.2e1 - t302;
t100 = Icges(6,1) * t170 - Icges(6,4) * t169 - Icges(6,5) * t208;
t96 = Icges(6,5) * t170 - Icges(6,6) * t169 - Icges(6,3) * t208;
t98 = Icges(6,4) * t170 - Icges(6,2) * t169 - Icges(6,6) * t208;
t43 = t100 * t215 - t214 * t98 - t226 * t96;
t53 = -t136 * t208 - t137 * t169 + t138 * t170;
t270 = -t43 / 0.2e1 - t53 / 0.2e1 - t301;
t269 = t91 * t309 + t90 * t310 + t281;
t268 = t256 * (t295 - t316);
t267 = t283 + t293;
t177 = -t210 * t255 - t257 * t309;
t178 = t210 * t257 - t297;
t110 = t178 * rSges(5,1) + t177 * rSges(5,2) - t208 * rSges(5,3);
t197 = t234 * rSges(3,1) + t233 * rSges(3,2) - rSges(3,3) * t309;
t266 = -t210 * t251 + t242 + t273;
t244 = rSges(2,1) * t263 - t261 * rSges(2,2);
t243 = -t261 * rSges(2,1) - rSges(2,2) * t263;
t225 = t258 * rSges(3,3) + (rSges(3,1) * t322 + rSges(3,2) * t323) * t256;
t196 = Icges(3,1) * t236 + Icges(3,4) * t235 + Icges(3,5) * t310;
t195 = Icges(3,1) * t234 + Icges(3,4) * t233 - Icges(3,5) * t309;
t194 = Icges(3,4) * t236 + Icges(3,2) * t235 + Icges(3,6) * t310;
t193 = Icges(3,4) * t234 + Icges(3,2) * t233 - Icges(3,6) * t309;
t176 = pkin(8) * t310 + t198 + t321;
t175 = -t261 * pkin(1) - t197 + t250;
t159 = -t258 * t197 - t225 * t309;
t158 = t198 * t258 - t225 * t310;
t146 = -rSges(4,3) * t309 - t277;
t145 = Icges(4,1) * t213 + Icges(4,4) * t211 + Icges(4,5) * t310;
t144 = Icges(4,1) * t210 + Icges(4,4) * t208 - Icges(4,5) * t309;
t143 = Icges(4,4) * t213 + Icges(4,2) * t211 + Icges(4,6) * t310;
t142 = Icges(4,4) * t210 + Icges(4,2) * t208 - Icges(4,6) * t309;
t135 = (t197 * t261 + t198 * t263) * t256;
t133 = t222 * t310 + t223 * t235 + t224 * t236;
t132 = -t222 * t309 + t233 * t223 + t234 * t224;
t117 = t283 + t147;
t116 = -t308 + (rSges(4,3) * t256 - t229) * t263 + t277;
t113 = t258 * t192 + (t194 * t323 + t196 * t322) * t256;
t112 = t258 * t191 + (t193 * t323 + t195 * t322) * t256;
t109 = Icges(5,1) * t180 + Icges(5,4) * t179 - Icges(5,5) * t211;
t108 = Icges(5,1) * t178 + Icges(5,4) * t177 - Icges(5,5) * t208;
t107 = Icges(5,4) * t180 + Icges(5,2) * t179 - Icges(5,6) * t211;
t106 = Icges(5,4) * t178 + Icges(5,2) * t177 - Icges(5,6) * t208;
t105 = Icges(5,5) * t180 + Icges(5,6) * t179 - Icges(5,3) * t211;
t104 = Icges(5,5) * t178 + Icges(5,6) * t177 - Icges(5,3) * t208;
t102 = -t208 * rSges(6,3) - t276;
t84 = (-t146 - t200) * t258 + t263 * t285;
t83 = t147 * t258 + t261 * t285 + t189;
t80 = t184 * t310 + t185 * t211 + t186 * t213;
t79 = -t184 * t309 + t208 * t185 + t210 * t186;
t78 = t283 - t284 + t111;
t77 = -t110 - t155 + t273;
t68 = t267 + t103;
t67 = (rSges(6,3) - t259) * t208 + t266 + t276;
t66 = (t146 * t261 + t147 * t263) * t256 + t305;
t65 = -t103 * t226 + t139 * t211;
t64 = t102 * t226 - t139 * t208;
t63 = t141 * t258 + t143 * t226 + t145 * t227;
t62 = t140 * t258 + t142 * t226 + t144 * t227;
t59 = t60 * t258;
t58 = t60 * t226;
t57 = -t102 * t211 + t103 * t208;
t56 = -t151 * t211 + t152 * t179 + t153 * t180;
t55 = -t151 * t208 + t152 * t177 + t153 * t178;
t52 = (-t110 + t306) * t258 + t263 * t280;
t51 = t111 * t258 + t261 * t280 + t307;
t50 = t267 + t317;
t49 = -t320 - t313 + (-rSges(7,3) - pkin(10)) * t169 + t266 + t275;
t48 = -t171 * t95 + t214 * t76;
t47 = t169 * t95 - t214 * t75;
t46 = -t105 * t226 + t107 * t216 + t109 * t217;
t45 = -t104 * t226 + t106 * t216 + t108 * t217;
t42 = (t110 * t261 + t111 * t263) * t256 + t281;
t41 = t101 * t172 - t171 * t99 - t211 * t97;
t40 = t100 * t172 - t171 * t98 - t211 * t96;
t39 = t101 * t170 - t169 * t99 - t208 * t97;
t38 = t100 * t170 - t169 * t98 - t208 * t96;
t36 = t37 * t258;
t35 = -t169 * t76 + t171 * t75;
t34 = t37 * t226;
t32 = t211 * t316 - t226 * t317;
t31 = -t208 * t316 + t226 * t318;
t28 = (-t102 + t299) * t258 + t263 * t272;
t27 = t103 * t258 + t261 * t272 + t300;
t26 = t208 * t317 - t211 * t318;
t23 = (t102 * t261 + t103 * t263) * t256 + t269;
t22 = t129 * t72 + t130 * t74 + t171 * t70;
t21 = t129 * t71 + t130 * t73 + t171 * t69;
t18 = (t299 - t318) * t258 + t263 * t268;
t17 = t258 * t317 + t261 * t268 + t300;
t16 = (t261 * t318 + t263 * t317) * t256 + t269;
t15 = t59 + (t44 * t261 - t43 * t263) * t256;
t14 = -t43 * t208 - t44 * t211 - t58;
t13 = t54 * t258 + (t261 * t41 - t263 * t40) * t256;
t12 = t53 * t258 + (t261 * t39 - t263 * t38) * t256;
t11 = -t208 * t40 - t211 * t41 - t226 * t54;
t10 = -t208 * t38 - t211 * t39 - t226 * t53;
t9 = t36 + (-t24 * t263 + t25 * t261) * t256;
t8 = -t24 * t208 - t25 * t211 - t34;
t6 = t30 * t258 + (-t21 * t263 + t22 * t261) * t256;
t5 = t29 * t258 + (-t19 * t263 + t20 * t261) * t256;
t4 = -t208 * t21 - t211 * t22 - t226 * t30;
t3 = -t19 * t208 - t20 * t211 - t226 * t29;
t2 = t169 * t21 + t171 * t22 + t214 * t30;
t61 = [m(7) * (t49 ^ 2 + t50 ^ 2) + m(6) * (t67 ^ 2 + t68 ^ 2) + m(5) * (t77 ^ 2 + t78 ^ 2) + m(4) * (t116 ^ 2 + t117 ^ 2) + m(3) * (t175 ^ 2 + t176 ^ 2) + m(2) * (t243 ^ 2 + t244 ^ 2) + Icges(2,3) + t60 + t37 + t337; t59 + t36 + m(7) * (t17 * t50 + t18 * t49) + m(6) * (t27 * t68 + t28 * t67) + m(5) * (t51 * t78 + t52 * t77) + m(4) * (t116 * t84 + t117 * t83) + m(3) * (t158 * t176 + t159 * t175) + ((-t112 / 0.2e1 - t62 / 0.2e1 - t45 / 0.2e1 - t55 / 0.2e1 - t79 / 0.2e1 - t132 / 0.2e1 + t270) * t263 + (t113 / 0.2e1 + t63 / 0.2e1 + t46 / 0.2e1 + t56 / 0.2e1 + t80 / 0.2e1 + t133 / 0.2e1 - t271) * t261) * t256 + t333; m(7) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(6) * (t23 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(5) * (t42 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(4) * (t66 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(3) * (t135 ^ 2 + t158 ^ 2 + t159 ^ 2) + (t9 + t15 + ((-t112 - t45 - t62) * t263 + (t113 + t46 + t63) * t261) * t256 + t333) * t258 + (t6 + t13 + ((-t105 * t211 + t107 * t179 + t109 * t180) * t261 - (-t104 * t211 + t106 * t179 + t108 * t180) * t263) * t256 + (t143 * t211 + t145 * t213 + t194 * t235 + t196 * t236 + t310 * t334) * t310 + (t56 + t133 + t80) * t258) * t310 + (-t5 - t12 - ((-t105 * t208 + t107 * t177 + t109 * t178) * t261 - (-t104 * t208 + t106 * t177 + t108 * t178) * t263) * t256 + (t208 * t142 + t210 * t144 + t233 * t193 + t234 * t195 + t309 * t335) * t309 + (-t79 - t55 - t132) * t258 + (-t142 * t211 - t208 * t143 - t144 * t213 - t210 * t145 - t193 * t235 - t233 * t194 - t195 * t236 - t234 * t196 + t309 * t334 + t310 * t335) * t310) * t309; 0.2e1 * ((t261 * t49 - t263 * t50) * t329 + (t261 * t67 - t263 * t68) * t330 + (t261 * t77 - t263 * t78) * t331 + (t116 * t261 - t117 * t263) * t332) * t256; m(7) * (t258 * t16 + (-t17 * t263 + t18 * t261) * t256) + m(6) * (t258 * t23 + (t261 * t28 - t263 * t27) * t256) + m(5) * (t258 * t42 + (t261 * t52 - t263 * t51) * t256) + m(4) * (t258 * t66 + (t261 * t84 - t263 * t83) * t256); 0.2e1 * (t332 + t282) * (t258 ^ 2 + (t261 ^ 2 + t263 ^ 2) * t256 ^ 2); m(7) * (-t208 * t50 - t211 * t49) + m(6) * (-t208 * t68 - t211 * t67) + m(5) * (-t208 * t78 - t211 * t77); m(7) * (-t16 * t226 - t17 * t208 - t18 * t211) + m(6) * (-t208 * t27 - t211 * t28 - t226 * t23) + m(5) * (-t208 * t51 - t211 * t52 - t226 * t42); (-t226 * t258 + (t208 * t263 - t211 * t261) * t256) * t336; (t208 ^ 2 + t211 ^ 2 + t226 ^ 2) * t336; -t34 - t58 + m(7) * (t31 * t49 + t32 * t50) + m(6) * (t64 * t67 + t65 * t68) + t271 * t211 + t270 * t208; (t8 / 0.2e1 + t14 / 0.2e1) * t258 - (t9 / 0.2e1 + t15 / 0.2e1) * t226 + (-t6 / 0.2e1 - t13 / 0.2e1) * t211 + (-t5 / 0.2e1 - t12 / 0.2e1) * t208 + m(7) * (t16 * t26 + t17 * t32 + t18 * t31) + m(6) * (t23 * t57 + t27 * t65 + t28 * t64) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t263 + (t4 / 0.2e1 + t11 / 0.2e1) * t261) * t256; m(6) * (t57 * t258 + (t261 * t64 - t263 * t65) * t256) + m(7) * (t26 * t258 + (t261 * t31 - t263 * t32) * t256); m(6) * (-t208 * t65 - t211 * t64 - t226 * t57) + m(7) * (-t208 * t32 - t211 * t31 - t226 * t26); -(t8 + t14) * t226 + (-t4 - t11) * t211 + (-t3 - t10) * t208 + m(7) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t57 ^ 2 + t64 ^ 2 + t65 ^ 2); t33 + m(7) * (t47 * t49 + t48 * t50) + t302 * t171 + t301 * t169; m(7) * (t16 * t35 + t17 * t48 + t18 * t47) + t258 * t327 + t9 * t324 + t5 * t326 + t6 * t325 + (t261 * t2 / 0.2e1 + t263 * t328) * t256; m(7) * (t35 * t258 + (t261 * t47 - t263 * t48) * t256); m(7) * (-t208 * t48 - t211 * t47 - t226 * t35); t208 * t328 - t226 * t327 + t3 * t326 + m(7) * (t26 * t35 + t31 * t47 + t32 * t48) - t211 * t2 / 0.2e1 + t8 * t324 + t4 * t325; m(7) * (t35 ^ 2 + t47 ^ 2 + t48 ^ 2) + t171 * t2 + t169 * t1 + t214 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t61(1) t61(2) t61(4) t61(7) t61(11) t61(16); t61(2) t61(3) t61(5) t61(8) t61(12) t61(17); t61(4) t61(5) t61(6) t61(9) t61(13) t61(18); t61(7) t61(8) t61(9) t61(10) t61(14) t61(19); t61(11) t61(12) t61(13) t61(14) t61(15) t61(20); t61(16) t61(17) t61(18) t61(19) t61(20) t61(21);];
Mq  = res;
