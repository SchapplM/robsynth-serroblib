% Calculate joint inertia matrix for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:14
% EndTime: 2019-03-09 10:38:30
% DurationCPUTime: 6.48s
% Computational Cost: add. (20125->607), mult. (51403->855), div. (0->0), fcn. (66973->12), ass. (0->274)
t315 = sin(pkin(11));
t316 = cos(pkin(11));
t322 = sin(qJ(2));
t324 = cos(qJ(2));
t248 = -t322 * t315 + t324 * t316;
t262 = sin(pkin(6));
t236 = t248 * t262;
t270 = t324 * t315 + t322 * t316;
t237 = t270 * t262;
t263 = cos(pkin(6));
t199 = Icges(4,5) * t237 + Icges(4,6) * t236 + Icges(4,3) * t263;
t200 = Icges(4,4) * t237 + Icges(4,2) * t236 + Icges(4,6) * t263;
t201 = Icges(4,1) * t237 + Icges(4,4) * t236 + Icges(4,5) * t263;
t232 = Icges(3,3) * t263 + (t322 * Icges(3,5) + t324 * Icges(3,6)) * t262;
t233 = Icges(3,6) * t263 + (t322 * Icges(3,4) + t324 * Icges(3,2)) * t262;
t234 = Icges(3,5) * t263 + (t322 * Icges(3,1) + t324 * Icges(3,4)) * t262;
t293 = t262 * t322;
t336 = t262 * t324 * t233 + t236 * t200 + t237 * t201 + t234 * t293 + (t199 + t232) * t263;
t330 = m(7) / 0.2e1;
t331 = m(6) / 0.2e1;
t302 = t331 + t330;
t335 = 0.2e1 * t302;
t334 = t262 ^ 2;
t333 = m(4) / 0.2e1;
t332 = m(5) / 0.2e1;
t238 = t270 * t263;
t265 = sin(qJ(1));
t267 = cos(qJ(1));
t223 = t238 * t267 + t265 * t248;
t321 = sin(qJ(4));
t292 = t262 * t321;
t323 = cos(qJ(4));
t192 = t223 * t323 - t267 * t292;
t225 = -t265 * t238 + t248 * t267;
t194 = t225 * t323 + t265 * t292;
t227 = t237 * t323 + t263 * t321;
t294 = t262 * t323;
t191 = t223 * t321 + t267 * t294;
t268 = t263 * t248;
t222 = -t265 * t270 + t267 * t268;
t264 = sin(qJ(6));
t266 = cos(qJ(6));
t146 = t191 * t266 + t222 * t264;
t147 = t191 * t264 - t222 * t266;
t86 = Icges(7,5) * t147 + Icges(7,6) * t146 + Icges(7,3) * t192;
t88 = Icges(7,4) * t147 + Icges(7,2) * t146 + Icges(7,6) * t192;
t90 = Icges(7,1) * t147 + Icges(7,4) * t146 + Icges(7,5) * t192;
t24 = t146 * t88 + t147 * t90 + t192 * t86;
t193 = t225 * t321 - t265 * t294;
t224 = -t265 * t268 - t267 * t270;
t148 = t193 * t266 + t224 * t264;
t149 = t193 * t264 - t224 * t266;
t87 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t194;
t89 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t194;
t91 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t194;
t25 = t146 * t89 + t147 * t91 + t192 * t87;
t226 = t237 * t321 - t263 * t323;
t180 = t226 * t266 + t236 * t264;
t181 = t226 * t264 - t236 * t266;
t120 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t227;
t121 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t227;
t122 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t227;
t35 = t120 * t192 + t121 * t146 + t122 * t147;
t1 = t192 * t24 + t194 * t25 + t227 * t35;
t329 = -t1 / 0.2e1;
t30 = t180 * t88 + t181 * t90 + t227 * t86;
t31 = t180 * t89 + t181 * t91 + t227 * t87;
t53 = t227 * t120 + t180 * t121 + t181 * t122;
t41 = t53 * t227;
t7 = t30 * t192 + t31 * t194 + t41;
t328 = t7 / 0.2e1;
t327 = t192 / 0.2e1;
t326 = t194 / 0.2e1;
t325 = t227 / 0.2e1;
t320 = t222 * pkin(5);
t319 = t267 * pkin(1);
t279 = -t147 * rSges(7,1) - t146 * rSges(7,2);
t92 = t192 * rSges(7,3) - t279;
t318 = t192 * pkin(10) - t320 + t92;
t93 = t149 * rSges(7,1) + t148 * rSges(7,2) + t194 * rSges(7,3);
t317 = -t224 * pkin(5) + pkin(10) * t194 + t93;
t314 = t262 * t265;
t313 = t262 * t267;
t260 = t324 * pkin(2) + pkin(1);
t312 = t265 * t260;
t239 = t263 * t322 * pkin(2) + (-pkin(8) - qJ(3)) * t262;
t311 = t267 * t239;
t310 = t336 * t263;
t123 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t227;
t309 = pkin(5) * t236 - pkin(10) * t227 - t123;
t173 = t225 * pkin(3) - pkin(9) * t224;
t253 = t267 * t260;
t214 = -t319 + t253 + (-t262 * pkin(8) - t239) * t265;
t203 = t263 * t214;
t308 = t263 * t173 + t203;
t172 = t223 * pkin(3) - t222 * pkin(9);
t258 = pkin(8) * t313;
t213 = t311 + t258 + (-pkin(1) + t260) * t265;
t307 = -t172 - t213;
t306 = t213 * t314 + t214 * t313;
t153 = Icges(4,5) * t223 + Icges(4,6) * t222 - Icges(4,3) * t313;
t289 = t267 * t324;
t290 = t265 * t322;
t243 = t263 * t289 - t290;
t288 = t267 * t322;
t291 = t265 * t324;
t244 = t263 * t288 + t291;
t204 = Icges(3,5) * t244 + Icges(3,6) * t243 - Icges(3,3) * t313;
t305 = -t204 - t153;
t154 = Icges(4,5) * t225 + Icges(4,6) * t224 + Icges(4,3) * t314;
t245 = -t263 * t291 - t288;
t246 = -t263 * t290 + t289;
t205 = Icges(3,5) * t246 + Icges(3,6) * t245 + Icges(3,3) * t314;
t304 = t205 + t154;
t249 = pkin(2) * t293 + t263 * qJ(3);
t303 = -pkin(3) * t237 + pkin(9) * t236 - t249;
t301 = t30 / 0.2e1 + t35 / 0.2e1;
t36 = t120 * t194 + t121 * t148 + t122 * t149;
t300 = t31 / 0.2e1 + t36 / 0.2e1;
t165 = Icges(5,5) * t227 - Icges(5,6) * t226 - Icges(5,3) * t236;
t167 = Icges(5,4) * t227 - Icges(5,2) * t226 - Icges(5,6) * t236;
t169 = Icges(5,1) * t227 - Icges(5,4) * t226 - Icges(5,5) * t236;
t78 = -t236 * t165 - t226 * t167 + t227 * t169;
t138 = t194 * pkin(4) + qJ(5) * t193;
t299 = t263 * t138 + t308;
t182 = t191 * qJ(5);
t137 = t192 * pkin(4) + t182;
t298 = -t137 + t307;
t164 = -Icges(6,5) * t236 - Icges(6,6) * t227 + Icges(6,3) * t226;
t166 = -Icges(6,4) * t236 - Icges(6,2) * t227 + Icges(6,6) * t226;
t168 = -Icges(6,1) * t236 - Icges(6,4) * t227 + Icges(6,5) * t226;
t77 = t226 * t164 - t227 * t166 - t236 * t168;
t179 = pkin(4) * t227 + qJ(5) * t226;
t296 = -t179 + t303;
t118 = t194 * rSges(5,1) - t193 * rSges(5,2) - t224 * rSges(5,3);
t116 = -t224 * rSges(6,1) - t194 * rSges(6,2) + t193 * rSges(6,3);
t160 = t225 * rSges(4,1) + t224 * rSges(4,2) + rSges(4,3) * t314;
t211 = t246 * rSges(3,1) + t245 * rSges(3,2) + rSges(3,3) * t314;
t287 = t262 * (-rSges(4,1) * t237 - rSges(4,2) * t236 - rSges(4,3) * t263 - t249);
t286 = -t265 * t239 + t253;
t285 = t172 * t314 + t173 * t313 + t306;
t171 = rSges(5,1) * t227 - rSges(5,2) * t226 - rSges(5,3) * t236;
t284 = t262 * (-t171 + t303);
t281 = -t223 * rSges(4,1) - t222 * rSges(4,2);
t280 = t222 * rSges(6,1) - t191 * rSges(6,3);
t170 = -rSges(6,1) * t236 - rSges(6,2) * t227 + rSges(6,3) * t226;
t278 = t262 * (-t170 + t296);
t277 = t137 * t314 + t138 * t313 + t285;
t276 = t262 * (t296 + t309);
t275 = t173 + t286;
t117 = t192 * rSges(5,1) - t191 * rSges(5,2) - t222 * rSges(5,3);
t210 = t244 * rSges(3,1) + t243 * rSges(3,2) - rSges(3,3) * t313;
t274 = -t172 - t311 - t312;
t103 = -Icges(6,5) * t222 - Icges(6,6) * t192 + Icges(6,3) * t191;
t107 = -Icges(6,4) * t222 - Icges(6,2) * t192 + Icges(6,6) * t191;
t111 = -Icges(6,1) * t222 - Icges(6,4) * t192 + Icges(6,5) * t191;
t56 = t103 * t226 - t107 * t227 - t111 * t236;
t105 = Icges(5,5) * t192 - Icges(5,6) * t191 - Icges(5,3) * t222;
t109 = Icges(5,4) * t192 - Icges(5,2) * t191 - Icges(5,6) * t222;
t113 = Icges(5,1) * t192 - Icges(5,4) * t191 - Icges(5,5) * t222;
t58 = -t105 * t236 - t109 * t226 + t113 * t227;
t66 = t164 * t191 - t166 * t192 - t168 * t222;
t68 = -t165 * t222 - t167 * t191 + t169 * t192;
t273 = -t58 / 0.2e1 - t56 / 0.2e1 - t68 / 0.2e1 - t66 / 0.2e1 - t301;
t104 = -Icges(6,5) * t224 - Icges(6,6) * t194 + Icges(6,3) * t193;
t108 = -Icges(6,4) * t224 - Icges(6,2) * t194 + Icges(6,6) * t193;
t112 = -Icges(6,1) * t224 - Icges(6,4) * t194 + Icges(6,5) * t193;
t57 = t104 * t226 - t108 * t227 - t112 * t236;
t106 = Icges(5,5) * t194 - Icges(5,6) * t193 - Icges(5,3) * t224;
t110 = Icges(5,4) * t194 - Icges(5,2) * t193 - Icges(5,6) * t224;
t114 = Icges(5,1) * t194 - Icges(5,4) * t193 - Icges(5,5) * t224;
t59 = -t106 * t236 - t110 * t226 + t114 * t227;
t67 = t164 * t193 - t166 * t194 - t168 * t224;
t69 = -t165 * t224 - t167 * t193 + t169 * t194;
t272 = -t69 / 0.2e1 - t67 / 0.2e1 - t59 / 0.2e1 - t57 / 0.2e1 - t300;
t271 = -t182 + t274;
t269 = t138 + t275;
t251 = rSges(2,1) * t267 - t265 * rSges(2,2);
t250 = -t265 * rSges(2,1) - rSges(2,2) * t267;
t235 = t263 * rSges(3,3) + (t322 * rSges(3,1) + t324 * rSges(3,2)) * t262;
t209 = Icges(3,1) * t246 + Icges(3,4) * t245 + Icges(3,5) * t314;
t208 = Icges(3,1) * t244 + Icges(3,4) * t243 - Icges(3,5) * t313;
t207 = Icges(3,4) * t246 + Icges(3,2) * t245 + Icges(3,6) * t314;
t206 = Icges(3,4) * t244 + Icges(3,2) * t243 - Icges(3,6) * t313;
t190 = pkin(8) * t314 + t211 + t319;
t189 = -t265 * pkin(1) - t210 + t258;
t176 = -t263 * t210 - t235 * t313;
t175 = t211 * t263 - t235 * t314;
t159 = -rSges(4,3) * t313 - t281;
t158 = Icges(4,1) * t225 + Icges(4,4) * t224 + Icges(4,5) * t314;
t157 = Icges(4,1) * t223 + Icges(4,4) * t222 - Icges(4,5) * t313;
t156 = Icges(4,4) * t225 + Icges(4,2) * t224 + Icges(4,6) * t314;
t155 = Icges(4,4) * t223 + Icges(4,2) * t222 - Icges(4,6) * t313;
t150 = (t210 * t265 + t211 * t267) * t262;
t145 = t232 * t314 + t233 * t245 + t234 * t246;
t144 = -t232 * t313 + t243 * t233 + t244 * t234;
t141 = t222 * t179;
t129 = t286 + t160;
t128 = -t312 + (rSges(4,3) * t262 - t239) * t267 + t281;
t127 = t236 * t138;
t125 = t263 * t205 + (t324 * t207 + t322 * t209) * t262;
t124 = t263 * t204 + (t324 * t206 + t322 * t208) * t262;
t119 = t224 * t137;
t115 = -t192 * rSges(6,2) - t280;
t97 = (-t159 - t213) * t263 + t267 * t287;
t96 = t263 * t160 + t265 * t287 + t203;
t95 = t199 * t314 + t200 * t224 + t201 * t225;
t94 = -t199 * t313 + t222 * t200 + t223 * t201;
t85 = t275 + t118;
t84 = -t117 + t274;
t83 = (t159 * t265 + t160 * t267) * t262 + t306;
t82 = -t118 * t236 + t171 * t224;
t81 = t117 * t236 - t171 * t222;
t80 = t154 * t263 + t156 * t236 + t158 * t237;
t79 = t153 * t263 + t155 * t236 + t157 * t237;
t76 = t78 * t263;
t75 = t77 * t263;
t74 = t78 * t236;
t73 = t77 * t236;
t72 = t269 + t116;
t71 = (rSges(6,2) - pkin(4)) * t192 + t271 + t280;
t70 = -t117 * t224 + t118 * t222;
t65 = (-t117 + t307) * t263 + t267 * t284;
t64 = t263 * t118 + t265 * t284 + t308;
t63 = -t123 * t194 + t227 * t93;
t62 = t123 * t192 - t227 * t92;
t61 = -t236 * t116 - t127 + (t170 + t179) * t224;
t60 = -t222 * t170 - t141 - (-t115 - t137) * t236;
t55 = t269 + t317;
t54 = t320 + (-rSges(7,3) - pkin(4) - pkin(10)) * t192 + t271 + t279;
t52 = (t117 * t265 + t118 * t267) * t262 + t285;
t51 = t53 * t263;
t50 = -t106 * t224 - t110 * t193 + t114 * t194;
t49 = -t105 * t224 - t109 * t193 + t113 * t194;
t48 = -t106 * t222 - t110 * t191 + t114 * t192;
t47 = -t105 * t222 - t109 * t191 + t113 * t192;
t46 = t104 * t193 - t108 * t194 - t112 * t224;
t45 = t103 * t193 - t107 * t194 - t111 * t224;
t44 = t104 * t191 - t108 * t192 - t112 * t222;
t43 = t103 * t191 - t107 * t192 - t111 * t222;
t42 = t53 * t236;
t40 = -t192 * t93 + t194 * t92;
t39 = (-t115 + t298) * t263 + t267 * t278;
t38 = t263 * t116 + t265 * t278 + t299;
t37 = -t224 * t115 - t119 + (t116 + t138) * t222;
t34 = (t115 * t265 + t116 * t267) * t262 + t277;
t33 = -t127 - t317 * t236 + (t179 - t309) * t224;
t32 = -t141 + t309 * t222 - (-t137 - t318) * t236;
t29 = (t298 - t318) * t263 + t267 * t276;
t28 = t263 * t317 + t265 * t276 + t299;
t27 = t148 * t89 + t149 * t91 + t194 * t87;
t26 = t148 * t88 + t149 * t90 + t194 * t86;
t23 = -t119 - t318 * t224 + (t138 + t317) * t222;
t22 = (t265 * t318 + t267 * t317) * t262 + t277;
t21 = t76 + (t59 * t265 - t58 * t267) * t262;
t20 = t75 + (t57 * t265 - t56 * t267) * t262;
t19 = -t58 * t222 - t59 * t224 - t74;
t18 = -t56 * t222 - t57 * t224 - t73;
t17 = t69 * t263 + (t265 * t50 - t267 * t49) * t262;
t16 = t68 * t263 + (t265 * t48 - t267 * t47) * t262;
t15 = t67 * t263 + (t265 * t46 - t267 * t45) * t262;
t14 = t66 * t263 + (t265 * t44 - t267 * t43) * t262;
t13 = -t222 * t49 - t224 * t50 - t236 * t69;
t12 = -t222 * t47 - t224 * t48 - t236 * t68;
t11 = -t222 * t45 - t224 * t46 - t236 * t67;
t10 = -t222 * t43 - t224 * t44 - t236 * t66;
t9 = t51 + (t31 * t265 - t30 * t267) * t262;
t8 = -t30 * t222 - t31 * t224 - t42;
t6 = t36 * t263 + (-t26 * t267 + t265 * t27) * t262;
t5 = t35 * t263 + (-t24 * t267 + t25 * t265) * t262;
t4 = -t222 * t26 - t224 * t27 - t236 * t36;
t3 = -t222 * t24 - t224 * t25 - t236 * t35;
t2 = t192 * t26 + t194 * t27 + t227 * t36;
t98 = [m(7) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t84 ^ 2 + t85 ^ 2) + m(6) * (t71 ^ 2 + t72 ^ 2) + m(4) * (t128 ^ 2 + t129 ^ 2) + m(3) * (t189 ^ 2 + t190 ^ 2) + m(2) * (t250 ^ 2 + t251 ^ 2) + Icges(2,3) + t77 + t78 + t53 + t336; t76 + t75 + t51 + m(7) * (t28 * t55 + t29 * t54) + m(6) * (t38 * t72 + t39 * t71) + m(5) * (t64 * t85 + t65 * t84) + m(4) * (t128 * t97 + t129 * t96) + m(3) * (t175 * t190 + t176 * t189) + ((-t94 / 0.2e1 - t144 / 0.2e1 - t124 / 0.2e1 - t79 / 0.2e1 + t273) * t267 + (t95 / 0.2e1 + t145 / 0.2e1 + t125 / 0.2e1 + t80 / 0.2e1 - t272) * t265) * t262 + t310; (t9 + t20 + t21 + t310) * t263 + m(7) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t34 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t52 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t83 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(3) * (t150 ^ 2 + t175 ^ 2 + t176 ^ 2) + ((-t14 - t16 - t5 + ((t222 * t155 + t223 * t157 + t243 * t206 + t244 * t208) * t262 + t305 * t334 * t267) * t267 + (-t124 - t144 - t79 - t94) * t263) * t267 + (t6 + t17 + t15 + ((t224 * t156 + t225 * t158 + t245 * t207 + t246 * t209) * t262 + t304 * t334 * t265) * t265 + (t145 + t95 + t125 + t80) * t263 + (-t224 * t155 - t222 * t156 - t225 * t157 - t223 * t158 - t245 * t206 - t243 * t207 - t246 * t208 - t244 * t209 + (t265 * t305 + t267 * t304) * t262) * t313) * t265) * t262; 0.2e1 * ((t265 * t54 - t267 * t55) * t330 + (t265 * t84 - t267 * t85) * t332 + (t265 * t71 - t267 * t72) * t331 + (t128 * t265 - t129 * t267) * t333) * t262; m(7) * (t263 * t22 + (t265 * t29 - t267 * t28) * t262) + m(6) * (t263 * t34 + (t265 * t39 - t267 * t38) * t262) + m(5) * (t263 * t52 + (t265 * t65 - t267 * t64) * t262) + m(4) * (t263 * t83 + (t265 * t97 - t267 * t96) * t262); 0.2e1 * (t333 + t332 + t302) * (t263 ^ 2 + (t265 ^ 2 + t267 ^ 2) * t334); -t42 - t74 - t73 + m(7) * (t32 * t54 + t33 * t55) + m(5) * (t81 * t84 + t82 * t85) + m(6) * (t60 * t71 + t61 * t72) + t272 * t224 + t273 * t222; (t8 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1) * t263 - (t9 / 0.2e1 + t20 / 0.2e1 + t21 / 0.2e1) * t236 + (-t6 / 0.2e1 - t15 / 0.2e1 - t17 / 0.2e1) * t224 + (-t5 / 0.2e1 - t14 / 0.2e1 - t16 / 0.2e1) * t222 + m(7) * (t22 * t23 + t28 * t33 + t29 * t32) + m(6) * (t34 * t37 + t38 * t61 + t39 * t60) + m(5) * (t52 * t70 + t64 * t82 + t65 * t81) + ((-t3 / 0.2e1 - t12 / 0.2e1 - t10 / 0.2e1) * t267 + (t4 / 0.2e1 + t13 / 0.2e1 + t11 / 0.2e1) * t265) * t262; m(5) * (t70 * t263 + (t265 * t81 - t267 * t82) * t262) + m(6) * (t37 * t263 + (t265 * t60 - t267 * t61) * t262) + m(7) * (t23 * t263 + (t265 * t32 - t267 * t33) * t262); -(t8 + t18 + t19) * t236 + (-t4 - t11 - t13) * t224 + (-t3 - t10 - t12) * t222 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t37 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t70 ^ 2 + t81 ^ 2 + t82 ^ 2); m(7) * (t191 * t55 + t193 * t54) + m(6) * (t191 * t72 + t193 * t71); m(7) * (t191 * t28 + t193 * t29 + t22 * t226) + m(6) * (t191 * t38 + t193 * t39 + t226 * t34); (t226 * t263 + (-t191 * t267 + t193 * t265) * t262) * t335; m(7) * (t191 * t33 + t193 * t32 + t226 * t23) + m(6) * (t191 * t61 + t193 * t60 + t226 * t37); (t191 ^ 2 + t193 ^ 2 + t226 ^ 2) * t335; m(7) * (t54 * t62 + t55 * t63) + t41 + t300 * t194 + t301 * t192; m(7) * (t22 * t40 + t28 * t63 + t29 * t62) + t5 * t327 + t6 * t326 + t9 * t325 + t263 * t328 + (t265 * t2 / 0.2e1 + t267 * t329) * t262; m(7) * (t40 * t263 + (t265 * t62 - t267 * t63) * t262); m(7) * (t23 * t40 + t32 * t62 + t33 * t63) + t4 * t326 + t3 * t327 + t8 * t325 - t224 * t2 / 0.2e1 - t236 * t328 + t222 * t329; m(7) * (t191 * t63 + t193 * t62 + t226 * t40); t194 * t2 + t192 * t1 + t227 * t7 + m(7) * (t40 ^ 2 + t62 ^ 2 + t63 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t98(1) t98(2) t98(4) t98(7) t98(11) t98(16); t98(2) t98(3) t98(5) t98(8) t98(12) t98(17); t98(4) t98(5) t98(6) t98(9) t98(13) t98(18); t98(7) t98(8) t98(9) t98(10) t98(14) t98(19); t98(11) t98(12) t98(13) t98(14) t98(15) t98(20); t98(16) t98(17) t98(18) t98(19) t98(20) t98(21);];
Mq  = res;
