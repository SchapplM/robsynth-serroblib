% Calculate joint inertia matrix for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:01
% EndTime: 2019-12-31 22:32:14
% DurationCPUTime: 4.88s
% Computational Cost: add. (21147->510), mult. (36316->714), div. (0->0), fcn. (45791->12), ass. (0->261)
t263 = cos(pkin(5));
t270 = cos(qJ(2));
t271 = cos(qJ(1));
t317 = t271 * t270;
t266 = sin(qJ(2));
t267 = sin(qJ(1));
t321 = t267 * t266;
t244 = -t263 * t317 + t321;
t318 = t271 * t266;
t320 = t267 * t270;
t246 = t263 * t320 + t318;
t262 = sin(pkin(5));
t324 = t262 * t270;
t245 = t263 * t318 + t320;
t316 = qJ(3) + qJ(4);
t260 = sin(t316);
t287 = cos(t316);
t323 = t262 * t271;
t216 = t245 * t287 - t260 * t323;
t264 = sin(qJ(5));
t268 = cos(qJ(5));
t180 = -t216 * t264 + t244 * t268;
t181 = t216 * t268 + t244 * t264;
t281 = t262 * t287;
t215 = t245 * t260 + t271 * t281;
t112 = Icges(6,5) * t181 + Icges(6,6) * t180 + Icges(6,3) * t215;
t114 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t215;
t116 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t215;
t46 = t112 * t215 + t114 * t180 + t116 * t181;
t247 = -t263 * t321 + t317;
t326 = t262 * t267;
t218 = t247 * t287 + t260 * t326;
t182 = -t218 * t264 + t246 * t268;
t183 = t218 * t268 + t246 * t264;
t217 = t247 * t260 - t267 * t281;
t113 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t217;
t115 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t217;
t117 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t217;
t47 = t113 * t215 + t115 * t180 + t117 * t181;
t234 = t260 * t263 + t266 * t281;
t213 = -t234 * t264 - t268 * t324;
t214 = t234 * t268 - t264 * t324;
t327 = t262 * t266;
t233 = t260 * t327 - t263 * t287;
t141 = Icges(6,5) * t214 + Icges(6,6) * t213 + Icges(6,3) * t233;
t142 = Icges(6,4) * t214 + Icges(6,2) * t213 + Icges(6,6) * t233;
t143 = Icges(6,1) * t214 + Icges(6,4) * t213 + Icges(6,5) * t233;
t62 = t141 * t215 + t142 * t180 + t143 * t181;
t11 = t244 * t46 + t246 * t47 - t324 * t62;
t146 = Icges(5,5) * t216 - Icges(5,6) * t215 + Icges(5,3) * t244;
t148 = Icges(5,4) * t216 - Icges(5,2) * t215 + Icges(5,6) * t244;
t150 = Icges(5,1) * t216 - Icges(5,4) * t215 + Icges(5,5) * t244;
t73 = t146 * t244 - t148 * t215 + t150 * t216;
t147 = Icges(5,5) * t218 - Icges(5,6) * t217 + Icges(5,3) * t246;
t149 = Icges(5,4) * t218 - Icges(5,2) * t217 + Icges(5,6) * t246;
t151 = Icges(5,1) * t218 - Icges(5,4) * t217 + Icges(5,5) * t246;
t74 = t147 * t244 - t149 * t215 + t151 * t216;
t186 = Icges(5,5) * t234 - Icges(5,6) * t233 - Icges(5,3) * t324;
t187 = Icges(5,4) * t234 - Icges(5,2) * t233 - Icges(5,6) * t324;
t188 = Icges(5,1) * t234 - Icges(5,4) * t233 - Icges(5,5) * t324;
t96 = t186 * t244 - t187 * t215 + t188 * t216;
t349 = t244 * t73 + t246 * t74 - t324 * t96 + t11;
t48 = t112 * t217 + t114 * t182 + t116 * t183;
t49 = t113 * t217 + t115 * t182 + t117 * t183;
t63 = t141 * t217 + t142 * t182 + t143 * t183;
t12 = t244 * t48 + t246 * t49 - t324 * t63;
t75 = t146 * t246 - t148 * t217 + t150 * t218;
t76 = t147 * t246 - t149 * t217 + t151 * t218;
t97 = t186 * t246 - t187 * t217 + t188 * t218;
t348 = t244 * t75 + t246 * t76 - t324 * t97 + t12;
t306 = -t187 * t233 + t188 * t234;
t102 = -t186 * t324 + t306;
t57 = t113 * t233 + t115 * t213 + t117 * t214;
t332 = t57 * t246;
t56 = t112 * t233 + t114 * t213 + t116 * t214;
t333 = t56 * t244;
t70 = t141 * t233 + t142 * t213 + t143 * t214;
t21 = -t324 * t70 + t332 + t333;
t86 = -t147 * t324 - t149 * t233 + t151 * t234;
t330 = t86 * t246;
t85 = -t146 * t324 - t148 * t233 + t150 * t234;
t331 = t85 * t244;
t347 = -t102 * t324 + t21 + t330 + t331;
t100 = t102 * t263;
t68 = t70 * t263;
t23 = t68 + (t57 * t267 - t56 * t271) * t262;
t346 = t23 + t100 + (t86 * t267 - t85 * t271) * t262;
t279 = -t181 * rSges(6,1) - t180 * rSges(6,2);
t118 = rSges(6,3) * t215 - t279;
t336 = t216 * pkin(4);
t314 = pkin(10) * t215 + t118 + t336;
t144 = rSges(6,1) * t214 + rSges(6,2) * t213 + rSges(6,3) * t233;
t345 = pkin(4) * t234 + pkin(10) * t233 + t144;
t344 = t215 / 0.2e1;
t343 = t217 / 0.2e1;
t342 = t233 / 0.2e1;
t341 = t244 / 0.2e1;
t340 = t246 / 0.2e1;
t339 = t263 / 0.2e1;
t338 = t267 / 0.2e1;
t337 = -t271 / 0.2e1;
t269 = cos(qJ(3));
t259 = pkin(3) * t269 + pkin(2);
t335 = -pkin(2) + t259;
t329 = -t102 - t70;
t272 = -pkin(9) - pkin(8);
t328 = t244 * t272;
t325 = t262 * t269;
t265 = sin(qJ(3));
t322 = t263 * t265;
t319 = t267 * t271;
t315 = t314 * t246;
t119 = rSges(6,1) * t183 + rSges(6,2) * t182 + rSges(6,3) * t217;
t313 = pkin(4) * t218 + pkin(10) * t217 + t119;
t240 = t244 * pkin(8);
t296 = t265 * t323;
t251 = pkin(3) * t296;
t154 = t245 * t335 - t240 - t251 - t328;
t201 = pkin(3) * t322 + ((pkin(8) + t272) * t270 + t335 * t266) * t262;
t312 = t154 * t324 + t201 * t244;
t208 = pkin(2) * t247 + pkin(8) * t246;
t297 = t265 * t326;
t292 = pkin(3) * t297 - t246 * t272 + t247 * t259;
t155 = -t208 + t292;
t206 = t263 * t208;
t310 = t155 * t263 + t206;
t153 = rSges(5,1) * t218 - rSges(5,2) * t217 + rSges(5,3) * t246;
t309 = -t153 - t155;
t207 = pkin(2) * t245 + t240;
t308 = -t154 - t207;
t221 = -t245 * t265 - t269 * t323;
t222 = t245 * t269 - t296;
t164 = rSges(4,1) * t222 + rSges(4,2) * t221 + rSges(4,3) * t244;
t307 = -t164 - t207;
t280 = -t216 * rSges(5,1) + t215 * rSges(5,2);
t152 = rSges(5,3) * t244 - t280;
t189 = rSges(5,1) * t234 - rSges(5,2) * t233 - rSges(5,3) * t324;
t109 = t152 * t324 + t189 * t244;
t242 = t263 * t269 - t265 * t327;
t243 = t266 * t325 + t322;
t191 = Icges(4,4) * t243 + Icges(4,2) * t242 - Icges(4,6) * t324;
t192 = Icges(4,1) * t243 + Icges(4,4) * t242 - Icges(4,5) * t324;
t305 = t191 * t242 + t192 * t243;
t304 = -t189 - t201;
t303 = t207 * t326 + t208 * t323;
t302 = pkin(1) * t271 + pkin(7) * t326;
t301 = t62 / 0.2e1 + t56 / 0.2e1;
t300 = t63 / 0.2e1 + t57 / 0.2e1;
t158 = Icges(4,5) * t222 + Icges(4,6) * t221 + Icges(4,3) * t244;
t160 = Icges(4,4) * t222 + Icges(4,2) * t221 + Icges(4,6) * t244;
t162 = Icges(4,1) * t222 + Icges(4,4) * t221 + Icges(4,5) * t244;
t91 = -t158 * t324 + t160 * t242 + t162 * t243;
t190 = Icges(4,5) * t243 + Icges(4,6) * t242 - Icges(4,3) * t324;
t98 = t190 * t244 + t191 * t221 + t192 * t222;
t299 = t98 / 0.2e1 + t91 / 0.2e1;
t223 = -t247 * t265 + t267 * t325;
t224 = t247 * t269 + t297;
t159 = Icges(4,5) * t224 + Icges(4,6) * t223 + Icges(4,3) * t246;
t161 = Icges(4,4) * t224 + Icges(4,2) * t223 + Icges(4,6) * t246;
t163 = Icges(4,1) * t224 + Icges(4,4) * t223 + Icges(4,5) * t246;
t92 = -t159 * t324 + t161 * t242 + t163 * t243;
t99 = t190 * t246 + t191 * t223 + t192 * t224;
t298 = t99 / 0.2e1 + t92 / 0.2e1;
t295 = -t155 - t313;
t294 = -t201 - t345;
t165 = rSges(4,1) * t224 + rSges(4,2) * t223 + rSges(4,3) * t246;
t229 = Icges(3,3) * t263 + (Icges(3,5) * t266 + Icges(3,6) * t270) * t262;
t230 = Icges(3,6) * t263 + (Icges(3,4) * t266 + Icges(3,2) * t270) * t262;
t231 = Icges(3,5) * t263 + (Icges(3,1) * t266 + Icges(3,4) * t270) * t262;
t293 = t229 * t263 + t230 * t324 + t231 * t327;
t203 = rSges(3,1) * t247 - rSges(3,2) * t246 + rSges(3,3) * t326;
t290 = -t324 / 0.2e1;
t288 = t244 * t349 + t246 * t348;
t286 = -t267 * pkin(1) + pkin(7) * t323;
t193 = rSges(4,1) * t243 + rSges(4,2) * t242 - rSges(4,3) * t324;
t248 = (pkin(2) * t266 - pkin(8) * t270) * t262;
t285 = t262 * (-t193 - t248);
t284 = t154 * t326 + t155 * t323 + t303;
t65 = t244 * t345 + t314 * t324;
t283 = t262 * (-t248 + t304);
t67 = t70 * t233;
t18 = t56 * t215 + t57 * t217 + t67;
t3 = t215 * t46 + t217 * t47 + t233 * t62;
t4 = t215 * t48 + t217 * t49 + t233 * t63;
t282 = t11 * t344 + t12 * t343 + t18 * t290 + t21 * t342 + t3 * t341 + t340 * t4;
t278 = t292 + t302;
t277 = t262 * (-t248 + t294);
t276 = t333 / 0.2e1 + t332 / 0.2e1 + t331 / 0.2e1 + t330 / 0.2e1 + (t62 + t96) * t341 + (t63 + t97) * t340;
t275 = -t245 * t259 + t251 + t286;
t274 = -t324 * t347 + t288;
t202 = rSges(3,1) * t245 - rSges(3,2) * t244 - rSges(3,3) * t323;
t15 = t62 * t263 + (t267 * t47 - t271 * t46) * t262;
t16 = t63 * t263 + (t267 * t49 - t271 * t48) * t262;
t32 = t96 * t263 + (t267 * t74 - t271 * t73) * t262;
t33 = t97 * t263 + (t267 * t76 - t271 * t75) * t262;
t273 = (t15 + t32) * t341 + (t16 + t33) * t340 + t347 * t339 + t348 * t326 / 0.2e1 + t346 * t290 - t349 * t323 / 0.2e1;
t253 = rSges(2,1) * t271 - rSges(2,2) * t267;
t252 = -rSges(2,1) * t267 - rSges(2,2) * t271;
t232 = t263 * rSges(3,3) + (rSges(3,1) * t266 + rSges(3,2) * t270) * t262;
t200 = Icges(3,1) * t247 - Icges(3,4) * t246 + Icges(3,5) * t326;
t199 = Icges(3,1) * t245 - Icges(3,4) * t244 - Icges(3,5) * t323;
t198 = Icges(3,4) * t247 - Icges(3,2) * t246 + Icges(3,6) * t326;
t197 = Icges(3,4) * t245 - Icges(3,2) * t244 - Icges(3,6) * t323;
t196 = Icges(3,5) * t247 - Icges(3,6) * t246 + Icges(3,3) * t326;
t195 = Icges(3,5) * t245 - Icges(3,6) * t244 - Icges(3,3) * t323;
t185 = t203 + t302;
t184 = -t202 + t286;
t168 = -t202 * t263 - t232 * t323;
t167 = t203 * t263 - t232 * t326;
t157 = t293 * t263;
t136 = (t202 * t267 + t203 * t271) * t262;
t135 = t229 * t326 - t230 * t246 + t231 * t247;
t134 = -t229 * t323 - t230 * t244 + t231 * t245;
t133 = t246 * t154;
t132 = t246 * t152;
t129 = t208 + t165 + t302;
t128 = t286 + t307;
t125 = -t165 * t324 - t193 * t246;
t124 = t164 * t324 + t193 * t244;
t123 = t263 * t196 + (t198 * t270 + t200 * t266) * t262;
t122 = t263 * t195 + (t197 * t270 + t199 * t266) * t262;
t121 = t278 + t153;
t120 = (-rSges(5,3) + t272) * t244 + t275 + t280;
t110 = -t153 * t324 - t246 * t189;
t107 = -t190 * t324 + t305;
t106 = t107 * t263;
t105 = t164 * t246 - t165 * t244;
t104 = t263 * t307 + t271 * t285;
t103 = t263 * t165 + t267 * t285 + t206;
t101 = -t153 * t244 + t132;
t95 = (t164 * t267 + t165 * t271) * t262 + t303;
t90 = t278 + t313;
t89 = -t336 + t328 + (-rSges(6,3) - pkin(10)) * t215 + t275 + t279;
t88 = t119 * t233 - t144 * t217;
t87 = -t118 * t233 + t144 * t215;
t84 = t246 * t304 + t309 * t324;
t83 = t109 + t312;
t82 = t159 * t246 + t161 * t223 + t163 * t224;
t81 = t158 * t246 + t160 * t223 + t162 * t224;
t80 = t159 * t244 + t161 * t221 + t163 * t222;
t79 = t158 * t244 + t160 * t221 + t162 * t222;
t72 = (-t152 + t308) * t263 + t271 * t283;
t71 = t263 * t153 + t267 * t283 + t310;
t69 = t118 * t217 - t119 * t215;
t66 = -t246 * t345 - t313 * t324;
t64 = t244 * t309 + t132 + t133;
t59 = (t152 * t267 + t153 * t271) * t262 + t284;
t58 = -t244 * t313 + t315;
t55 = t246 * t294 + t295 * t324;
t54 = t65 + t312;
t51 = (t308 - t314) * t263 + t271 * t277;
t50 = t263 * t313 + t267 * t277 + t310;
t45 = t244 * t295 + t133 + t315;
t44 = (t267 * t314 + t271 * t313) * t262 + t284;
t43 = t106 + (t92 * t267 - t91 * t271) * t262;
t42 = -t107 * t324 + t91 * t244 + t92 * t246;
t39 = t99 * t263 + (t267 * t82 - t271 * t81) * t262;
t38 = t98 * t263 + (t267 * t80 - t271 * t79) * t262;
t35 = t244 * t81 + t246 * t82 - t324 * t99;
t34 = t244 * t79 + t246 * t80 - t324 * t98;
t1 = [Icges(2,3) + (-t186 - t190) * t324 + m(6) * (t89 ^ 2 + t90 ^ 2) + m(5) * (t120 ^ 2 + t121 ^ 2) + m(4) * (t128 ^ 2 + t129 ^ 2) + m(3) * (t184 ^ 2 + t185 ^ 2) + m(2) * (t252 ^ 2 + t253 ^ 2) + t293 + t70 + t305 + t306; t68 + t100 + t106 + t157 + m(6) * (t50 * t90 + t51 * t89) + m(5) * (t120 * t72 + t121 * t71) + m(4) * (t103 * t129 + t104 * t128) + m(3) * (t167 * t185 + t168 * t184) + ((-t85 / 0.2e1 - t96 / 0.2e1 - t122 / 0.2e1 - t134 / 0.2e1 - t299 - t301) * t271 + (t97 / 0.2e1 + t86 / 0.2e1 + t123 / 0.2e1 + t135 / 0.2e1 + t298 + t300) * t267) * t262; (t43 + t157 + t346) * t263 + m(6) * (t44 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t59 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t103 ^ 2 + t104 ^ 2 + t95 ^ 2) + m(3) * (t136 ^ 2 + t167 ^ 2 + t168 ^ 2) + (-t271 * t15 + t267 * t16 - t271 * t32 + t267 * t33 - t271 * t38 + t267 * t39 + (t267 * ((-t198 * t246 + t200 * t247) * t267 - (-t197 * t246 + t199 * t247) * t271) - t271 * ((-t198 * t244 + t200 * t245) * t267 - (-t197 * t244 + t199 * t245) * t271) + (t267 * (t196 * t267 ^ 2 - t195 * t319) - t271 * (t195 * t271 ^ 2 - t196 * t319)) * t262) * t262 + ((-t122 - t134) * t271 + (t123 + t135) * t267) * t263) * t262; t298 * t246 + t299 * t244 + (-t107 + t329) * t324 + m(6) * (t54 * t89 + t55 * t90) + m(5) * (t120 * t83 + t121 * t84) + m(4) * (t124 * t128 + t125 * t129) + t276; t273 + m(6) * (t44 * t45 + t50 * t55 + t51 * t54) + m(5) * (t64 * t59 + t71 * t84 + t72 * t83) + m(4) * (t103 * t125 + t104 * t124 + t105 * t95) + (-t270 * t43 / 0.2e1 + t34 * t337 + t35 * t338) * t262 + t38 * t341 + t39 * t340 + t42 * t339; t244 * t34 + t246 * t35 + (-t42 - t347) * t324 + m(6) * (t45 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t64 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t105 ^ 2 + t124 ^ 2 + t125 ^ 2) + t288; t329 * t324 + m(6) * (t65 * t89 + t66 * t90) + m(5) * (t109 * t120 + t110 * t121) + t276; m(6) * (t58 * t44 + t50 * t66 + t65 * t51) + m(5) * (t101 * t59 + t109 * t72 + t110 * t71) + t273; m(6) * (t58 * t45 + t65 * t54 + t55 * t66) + m(5) * (t101 * t64 + t109 * t83 + t110 * t84) + t274; m(6) * (t58 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t101 ^ 2 + t109 ^ 2 + t110 ^ 2) + t274; m(6) * (t87 * t89 + t88 * t90) + t67 + t300 * t217 + t301 * t215; m(6) * (t44 * t69 + t50 * t88 + t51 * t87) + t15 * t344 + t18 * t339 + t23 * t342 + t16 * t343 + (t3 * t337 + t338 * t4) * t262; m(6) * (t45 * t69 + t54 * t87 + t55 * t88) + t282; m(6) * (t58 * t69 + t65 * t87 + t66 * t88) + t282; m(6) * (t69 ^ 2 + t87 ^ 2 + t88 ^ 2) + t217 * t4 + t215 * t3 + t233 * t18;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
