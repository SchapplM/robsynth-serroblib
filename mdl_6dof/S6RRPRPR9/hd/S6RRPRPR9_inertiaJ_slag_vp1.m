% Calculate joint inertia matrix for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:39
% EndTime: 2019-03-09 10:55:54
% DurationCPUTime: 6.04s
% Computational Cost: add. (22736->623), mult. (34588->867), div. (0->0), fcn. (43483->14), ass. (0->280)
t269 = sin(pkin(11));
t272 = cos(pkin(11));
t273 = cos(pkin(6));
t270 = sin(pkin(6));
t276 = sin(qJ(2));
t327 = t270 * t276;
t244 = -t269 * t327 + t272 * t273;
t328 = t269 * t273;
t245 = t272 * t327 + t328;
t278 = cos(qJ(2));
t325 = t270 * t278;
t185 = Icges(4,4) * t245 + Icges(4,2) * t244 - Icges(4,6) * t325;
t186 = Icges(4,1) * t245 + Icges(4,4) * t244 - Icges(4,5) * t325;
t230 = Icges(3,3) * t273 + (Icges(3,5) * t276 + Icges(3,6) * t278) * t270;
t231 = Icges(3,6) * t273 + (Icges(3,4) * t276 + Icges(3,2) * t278) * t270;
t232 = Icges(3,5) * t273 + (Icges(3,1) * t276 + Icges(3,4) * t278) * t270;
t344 = t244 * t185 + t245 * t186 + t273 * t230 + t231 * t325 + t232 * t327;
t307 = m(6) / 0.2e1 + m(7) / 0.2e1;
t343 = 0.2e1 * t307;
t308 = pkin(11) + qJ(4);
t264 = sin(t308);
t293 = cos(t308);
t287 = t270 * t293;
t234 = t273 * t264 + t276 * t287;
t267 = pkin(12) + qJ(6);
t263 = sin(t267);
t265 = cos(t267);
t203 = -t234 * t263 - t265 * t325;
t204 = t234 * t265 - t263 * t325;
t233 = t264 * t327 - t273 * t293;
t118 = Icges(7,5) * t204 + Icges(7,6) * t203 + Icges(7,3) * t233;
t119 = Icges(7,4) * t204 + Icges(7,2) * t203 + Icges(7,6) * t233;
t120 = Icges(7,1) * t204 + Icges(7,4) * t203 + Icges(7,5) * t233;
t53 = t233 * t118 + t203 * t119 + t204 * t120;
t268 = sin(pkin(12));
t271 = cos(pkin(12));
t214 = -t234 * t268 - t271 * t325;
t303 = t268 * t325;
t215 = t234 * t271 - t303;
t125 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t233;
t126 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t233;
t127 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t233;
t56 = t233 * t125 + t214 * t126 + t215 * t127;
t342 = -t53 - t56;
t184 = Icges(4,5) * t245 + Icges(4,6) * t244 - Icges(4,3) * t325;
t341 = (-t184 * t325 + t344) * t273;
t279 = cos(qJ(1));
t320 = t278 * t279;
t277 = sin(qJ(1));
t323 = t276 * t277;
t249 = -t273 * t323 + t320;
t218 = t249 * t264 - t277 * t287;
t326 = t270 * t277;
t219 = t249 * t293 + t264 * t326;
t261 = pkin(5) * t271 + pkin(4);
t274 = -pkin(10) - qJ(5);
t321 = t277 * t278;
t322 = t276 * t279;
t248 = t273 * t321 + t322;
t329 = t248 * t268;
t170 = -t219 * t263 + t248 * t265;
t171 = t219 * t265 + t248 * t263;
t95 = t171 * rSges(7,1) + t170 * rSges(7,2) + t218 * rSges(7,3);
t340 = pkin(5) * t329 - t218 * t274 + t219 * t261 + t95;
t247 = t273 * t322 + t321;
t216 = t247 * t264 + t279 * t287;
t324 = t270 * t279;
t217 = t247 * t293 - t264 * t324;
t246 = -t273 * t320 + t323;
t168 = -t217 * t263 + t246 * t265;
t169 = t217 * t265 + t246 * t263;
t88 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t216;
t90 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t216;
t92 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t216;
t31 = t170 * t90 + t171 * t92 + t218 * t88;
t89 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t218;
t91 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t218;
t93 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t218;
t32 = t170 * t91 + t171 * t93 + t218 * t89;
t45 = t118 * t218 + t119 * t170 + t120 * t171;
t2 = t216 * t31 + t218 * t32 + t233 * t45;
t339 = t2 / 0.2e1;
t338 = t216 / 0.2e1;
t337 = t218 / 0.2e1;
t336 = t233 / 0.2e1;
t262 = pkin(3) * t272 + pkin(2);
t335 = -pkin(2) + t262;
t334 = -pkin(4) + t261;
t208 = t216 * qJ(5);
t331 = t246 * t268;
t288 = -t169 * rSges(7,1) - t168 * rSges(7,2);
t94 = t216 * rSges(7,3) - t288;
t333 = pkin(5) * t331 - t216 * t274 + t217 * t334 - t208 + t94;
t162 = t219 * pkin(4) + qJ(5) * t218;
t332 = -t162 + t340;
t275 = -pkin(9) - qJ(3);
t330 = t246 * t275;
t176 = -t219 * t268 + t248 * t271;
t177 = t219 * t271 + t329;
t103 = t177 * rSges(6,1) + t176 * rSges(6,2) + t218 * rSges(6,3);
t319 = -t103 - t162;
t121 = rSges(7,1) * t204 + rSges(7,2) * t203 + rSges(7,3) * t233;
t318 = -pkin(5) * t303 + t334 * t234 + (-qJ(5) - t274) * t233 + t121;
t207 = t249 * pkin(2) + qJ(3) * t248;
t302 = t269 * t326;
t295 = pkin(3) * t302 - t248 * t275 + t249 * t262;
t142 = -t207 + t295;
t202 = t273 * t207;
t317 = t273 * t142 + t202;
t238 = t246 * qJ(3);
t301 = t269 * t324;
t253 = pkin(3) * t301;
t141 = t247 * t335 - t238 - t253 - t330;
t206 = pkin(2) * t247 + t238;
t316 = -t141 - t206;
t222 = -t247 * t269 - t272 * t324;
t223 = t247 * t272 - t301;
t151 = rSges(4,1) * t223 + rSges(4,2) * t222 + rSges(4,3) * t246;
t315 = -t151 - t206;
t161 = pkin(4) * t217 + t208;
t187 = t234 * pkin(4) + t233 * qJ(5);
t314 = t161 * t325 + t246 * t187;
t181 = Icges(5,4) * t234 - Icges(5,2) * t233 - Icges(5,6) * t325;
t182 = Icges(5,1) * t234 - Icges(5,4) * t233 - Icges(5,5) * t325;
t313 = -t233 * t181 + t234 * t182;
t250 = (pkin(2) * t276 - qJ(3) * t278) * t270;
t311 = -pkin(3) * t328 - ((qJ(3) + t275) * t278 + t335 * t276) * t270 - t250;
t310 = t206 * t326 + t207 * t324;
t309 = t279 * pkin(1) + pkin(8) * t326;
t306 = -t162 - t332;
t39 = t203 * t90 + t204 * t92 + t233 * t88;
t44 = t118 * t216 + t119 * t168 + t120 * t169;
t305 = t39 / 0.2e1 + t44 / 0.2e1;
t40 = t203 * t91 + t204 * t93 + t233 * t89;
t304 = t40 / 0.2e1 + t45 / 0.2e1;
t300 = t273 * t162 + t317;
t299 = -t161 + t316;
t298 = -t187 + t311;
t140 = t219 * rSges(5,1) - t218 * rSges(5,2) + t248 * rSges(5,3);
t224 = -t249 * t269 + t272 * t326;
t225 = t249 * t272 + t302;
t152 = t225 * rSges(4,1) + t224 * rSges(4,2) + t248 * rSges(4,3);
t197 = t249 * rSges(3,1) - t248 * rSges(3,2) + rSges(3,3) * t326;
t294 = -t277 * pkin(1) + pkin(8) * t324;
t292 = t270 * (-rSges(4,1) * t245 - rSges(4,2) * t244 + rSges(4,3) * t325 - t250);
t291 = t141 * t326 + t142 * t324 + t310;
t183 = rSges(5,1) * t234 - rSges(5,2) * t233 - rSges(5,3) * t325;
t290 = t270 * (-t183 + t311);
t289 = -t217 * rSges(5,1) + t216 * rSges(5,2);
t286 = t295 + t309;
t128 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t233;
t285 = t270 * (-t128 + t298);
t284 = t161 * t326 + t162 * t324 + t291;
t283 = t270 * (t298 - t318);
t282 = -t247 * t262 + t253 + t294;
t174 = -t217 * t268 + t246 * t271;
t175 = t217 * t271 + t331;
t102 = rSges(6,1) * t175 + rSges(6,2) * t174 + rSges(6,3) * t216;
t196 = t247 * rSges(3,1) - t246 * rSges(3,2) - rSges(3,3) * t324;
t100 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t216;
t96 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t216;
t98 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t216;
t41 = t100 * t215 + t214 * t98 + t233 * t96;
t47 = t125 * t216 + t126 * t174 + t127 * t175;
t133 = Icges(5,5) * t217 - Icges(5,6) * t216 + Icges(5,3) * t246;
t135 = Icges(5,4) * t217 - Icges(5,2) * t216 + Icges(5,6) * t246;
t137 = Icges(5,1) * t217 - Icges(5,4) * t216 + Icges(5,5) * t246;
t67 = -t133 * t325 - t135 * t233 + t137 * t234;
t180 = Icges(5,5) * t234 - Icges(5,6) * t233 - Icges(5,3) * t325;
t74 = t180 * t246 - t181 * t216 + t182 * t217;
t281 = t74 / 0.2e1 + t47 / 0.2e1 + t67 / 0.2e1 + t41 / 0.2e1 + t305;
t101 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t218;
t97 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t218;
t99 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t218;
t42 = t101 * t215 + t214 * t99 + t233 * t97;
t48 = t125 * t218 + t126 * t176 + t127 * t177;
t134 = Icges(5,5) * t219 - Icges(5,6) * t218 + Icges(5,3) * t248;
t136 = Icges(5,4) * t219 - Icges(5,2) * t218 + Icges(5,6) * t248;
t138 = Icges(5,1) * t219 - Icges(5,4) * t218 + Icges(5,5) * t248;
t68 = -t134 * t325 - t136 * t233 + t138 * t234;
t75 = t180 * t248 - t181 * t218 + t182 * t219;
t280 = t68 / 0.2e1 + t42 / 0.2e1 + t48 / 0.2e1 + t75 / 0.2e1 + t304;
t255 = rSges(2,1) * t279 - t277 * rSges(2,2);
t254 = -t277 * rSges(2,1) - rSges(2,2) * t279;
t235 = t273 * rSges(3,3) + (rSges(3,1) * t276 + rSges(3,2) * t278) * t270;
t195 = Icges(3,1) * t249 - Icges(3,4) * t248 + Icges(3,5) * t326;
t194 = Icges(3,1) * t247 - Icges(3,4) * t246 - Icges(3,5) * t324;
t193 = Icges(3,4) * t249 - Icges(3,2) * t248 + Icges(3,6) * t326;
t192 = Icges(3,4) * t247 - Icges(3,2) * t246 - Icges(3,6) * t324;
t191 = Icges(3,5) * t249 - Icges(3,6) * t248 + Icges(3,3) * t326;
t190 = Icges(3,5) * t247 - Icges(3,6) * t246 - Icges(3,3) * t324;
t179 = t197 + t309;
t178 = -t196 + t294;
t160 = -t273 * t196 - t235 * t324;
t159 = t197 * t273 - t235 * t326;
t150 = Icges(4,1) * t225 + Icges(4,4) * t224 + Icges(4,5) * t248;
t149 = Icges(4,1) * t223 + Icges(4,4) * t222 + Icges(4,5) * t246;
t148 = Icges(4,4) * t225 + Icges(4,2) * t224 + Icges(4,6) * t248;
t147 = Icges(4,4) * t223 + Icges(4,2) * t222 + Icges(4,6) * t246;
t146 = Icges(4,5) * t225 + Icges(4,6) * t224 + Icges(4,3) * t248;
t145 = Icges(4,5) * t223 + Icges(4,6) * t222 + Icges(4,3) * t246;
t143 = t248 * t161;
t139 = rSges(5,3) * t246 - t289;
t124 = (t196 * t277 + t197 * t279) * t270;
t123 = t230 * t326 - t231 * t248 + t232 * t249;
t122 = -t230 * t324 - t246 * t231 + t247 * t232;
t114 = t207 + t152 + t309;
t113 = t294 + t315;
t107 = t273 * t191 + (t193 * t278 + t195 * t276) * t270;
t106 = t273 * t190 + (t192 * t278 + t194 * t276) * t270;
t105 = t286 + t140;
t104 = (-rSges(5,3) + t275) * t246 + t282 + t289;
t87 = -t140 * t325 - t183 * t248;
t86 = t139 * t325 + t183 * t246;
t82 = t273 * t315 + t279 * t292;
t81 = t273 * t152 + t277 * t292 + t202;
t80 = -t180 * t325 + t313;
t79 = t139 * t248 - t140 * t246;
t78 = t80 * t273;
t77 = t184 * t248 + t185 * t224 + t186 * t225;
t76 = t184 * t246 + t185 * t222 + t186 * t223;
t73 = (t151 * t277 + t152 * t279) * t270 + t310;
t72 = -t146 * t325 + t148 * t244 + t150 * t245;
t71 = -t145 * t325 + t147 * t244 + t149 * t245;
t70 = t286 - t319;
t69 = -t102 - t161 + t282 + t330;
t66 = -t121 * t218 + t233 * t95;
t65 = t121 * t216 - t233 * t94;
t64 = t286 + t340;
t63 = -t217 * t261 + (-pkin(5) * t268 + t275) * t246 + (-rSges(7,3) + t274) * t216 + t282 + t288;
t62 = t134 * t248 - t136 * t218 + t138 * t219;
t61 = t133 * t248 - t135 * t218 + t137 * t219;
t60 = t134 * t246 - t136 * t216 + t138 * t217;
t59 = t133 * t246 - t135 * t216 + t137 * t217;
t58 = (-t139 + t316) * t273 + t279 * t290;
t57 = t273 * t140 + t277 * t290 + t317;
t55 = t56 * t273;
t54 = -t216 * t95 + t218 * t94;
t52 = t319 * t325 + (-t128 - t187) * t248;
t51 = t102 * t325 + t128 * t246 + t314;
t50 = t53 * t273;
t49 = t53 * t233;
t46 = (t139 * t277 + t140 * t279) * t270 + t291;
t43 = t248 * t102 + t246 * t319 + t143;
t38 = (-t102 + t299) * t273 + t279 * t285;
t37 = t273 * t103 + t277 * t285 + t300;
t36 = t101 * t177 + t176 * t99 + t218 * t97;
t35 = t100 * t177 + t176 * t98 + t218 * t96;
t34 = t101 * t175 + t174 * t99 + t216 * t97;
t33 = t100 * t175 + t174 * t98 + t216 * t96;
t30 = t168 * t91 + t169 * t93 + t216 * t89;
t29 = t168 * t90 + t169 * t92 + t216 * t88;
t28 = t306 * t325 + (-t187 - t318) * t248;
t27 = t246 * t318 + t325 * t333 + t314;
t26 = (t102 * t277 + t103 * t279) * t270 + t284;
t25 = t78 + (t68 * t277 - t67 * t279) * t270;
t24 = t67 * t246 + t68 * t248 - t325 * t80;
t23 = (t299 - t333) * t273 + t279 * t283;
t22 = t273 * t332 + t277 * t283 + t300;
t21 = t246 * t306 + t248 * t333 + t143;
t20 = t75 * t273 + (t277 * t62 - t279 * t61) * t270;
t19 = t74 * t273 + (t277 * t60 - t279 * t59) * t270;
t18 = t246 * t61 + t248 * t62 - t325 * t75;
t17 = t246 * t59 + t248 * t60 - t325 * t74;
t16 = (t277 * t333 + t279 * t332) * t270 + t284;
t15 = t55 + (t42 * t277 - t41 * t279) * t270;
t14 = t41 * t246 + t42 * t248 - t325 * t56;
t13 = t50 + (t40 * t277 - t39 * t279) * t270;
t12 = t39 * t246 + t40 * t248 - t325 * t53;
t11 = t39 * t216 + t40 * t218 + t49;
t10 = t48 * t273 + (t277 * t36 - t279 * t35) * t270;
t9 = t47 * t273 + (t277 * t34 - t279 * t33) * t270;
t8 = t246 * t35 + t248 * t36 - t325 * t48;
t7 = t246 * t33 + t248 * t34 - t325 * t47;
t6 = t45 * t273 + (t277 * t32 - t279 * t31) * t270;
t5 = t44 * t273 + (t277 * t30 - t279 * t29) * t270;
t4 = t246 * t31 + t248 * t32 - t325 * t45;
t3 = t246 * t29 + t248 * t30 - t325 * t44;
t1 = t216 * t29 + t218 * t30 + t233 * t44;
t83 = [m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t104 ^ 2 + t105 ^ 2) + m(4) * (t113 ^ 2 + t114 ^ 2) + m(3) * (t178 ^ 2 + t179 ^ 2) + m(2) * (t254 ^ 2 + t255 ^ 2) + (-t180 - t184) * t325 + Icges(2,3) + t313 - t342 + t344; t50 + t78 + t55 + m(7) * (t22 * t64 + t23 * t63) + m(6) * (t37 * t70 + t38 * t69) + m(5) * (t104 * t58 + t105 * t57) + m(4) * (t113 * t82 + t114 * t81) + m(3) * (t159 * t179 + t160 * t178) + ((-t71 / 0.2e1 - t106 / 0.2e1 - t76 / 0.2e1 - t122 / 0.2e1 - t281) * t279 + (t72 / 0.2e1 + t107 / 0.2e1 + t77 / 0.2e1 + t123 / 0.2e1 + t280) * t277) * t270 + t341; m(7) * (t16 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(6) * (t26 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t46 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(4) * (t73 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(3) * (t124 ^ 2 + t159 ^ 2 + t160 ^ 2) + (t6 + t10 + t20 + ((t146 * t248 + t224 * t148 + t225 * t150) * t277 - (t248 * t145 + t224 * t147 + t225 * t149) * t279) * t270 + (t191 * t326 - t193 * t248 + t249 * t195) * t326) * t326 + (-t5 - t9 - t19 - ((t246 * t146 + t222 * t148 + t223 * t150) * t277 - (t145 * t246 + t147 * t222 + t149 * t223) * t279) * t270 + (-t190 * t324 - t246 * t192 + t247 * t194) * t324 + (-t190 * t326 + t191 * t324 + t248 * t192 + t246 * t193 - t249 * t194 - t247 * t195) * t326) * t324 + (t13 + t15 + t25 + (t123 + t77) * t326 + (-t122 - t76) * t324 + ((-t106 - t71) * t279 + (t107 + t72) * t277) * t270 + t341) * t273; m(7) * (t246 * t64 + t248 * t63) + m(6) * (t246 * t70 + t248 * t69) + m(5) * (t104 * t248 + t105 * t246) + m(4) * (t113 * t248 + t114 * t246); m(7) * (-t16 * t325 + t22 * t246 + t23 * t248) + m(6) * (t246 * t37 + t248 * t38 - t26 * t325) + m(5) * (t246 * t57 + t248 * t58 - t325 * t46) + m(4) * (t246 * t81 + t248 * t82 - t325 * t73); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t307) * (t270 ^ 2 * t278 ^ 2 + t246 ^ 2 + t248 ^ 2); (-t80 + t342) * t325 + m(7) * (t27 * t63 + t28 * t64) + m(6) * (t51 * t69 + t52 * t70) + m(5) * (t104 * t86 + t105 * t87) + t280 * t248 + t281 * t246; (t12 / 0.2e1 + t14 / 0.2e1 + t24 / 0.2e1) * t273 + (t6 / 0.2e1 + t10 / 0.2e1 + t20 / 0.2e1) * t248 + (t5 / 0.2e1 + t19 / 0.2e1 + t9 / 0.2e1) * t246 + m(7) * (t16 * t21 + t22 * t28 + t23 * t27) + m(6) * (t26 * t43 + t37 * t52 + t38 * t51) + m(5) * (t46 * t79 + t57 * t87 + t58 * t86) + ((-t3 / 0.2e1 - t7 / 0.2e1 - t17 / 0.2e1) * t279 + (-t13 / 0.2e1 - t15 / 0.2e1 - t25 / 0.2e1) * t278 + (t4 / 0.2e1 + t8 / 0.2e1 + t18 / 0.2e1) * t277) * t270; m(5) * (t246 * t87 + t248 * t86 - t325 * t79) + m(6) * (t246 * t52 + t248 * t51 - t325 * t43) + m(7) * (-t21 * t325 + t246 * t28 + t248 * t27); (-t12 - t14 - t24) * t325 + (t4 + t18 + t8) * t248 + (t3 + t7 + t17) * t246 + m(7) * (t21 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t43 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t79 ^ 2 + t86 ^ 2 + t87 ^ 2); m(7) * (t216 * t64 + t218 * t63) + m(6) * (t216 * t70 + t218 * t69); m(7) * (t16 * t233 + t216 * t22 + t218 * t23) + m(6) * (t216 * t37 + t218 * t38 + t233 * t26); (t216 * t246 + t218 * t248 - t233 * t325) * t343; m(7) * (t21 * t233 + t216 * t28 + t218 * t27) + m(6) * (t216 * t52 + t218 * t51 + t233 * t43); (t216 ^ 2 + t218 ^ 2 + t233 ^ 2) * t343; m(7) * (t63 * t65 + t64 * t66) + t49 + t304 * t218 + t305 * t216; t6 * t337 + m(7) * (t16 * t54 + t22 * t66 + t23 * t65) + t273 * t11 / 0.2e1 + t5 * t338 + t13 * t336 + (t277 * t339 - t279 * t1 / 0.2e1) * t270; m(7) * (t246 * t66 + t248 * t65 - t325 * t54); -t11 * t325 / 0.2e1 + m(7) * (t21 * t54 + t27 * t65 + t28 * t66) + t12 * t336 + t246 * t1 / 0.2e1 + t248 * t339 + t3 * t338 + t4 * t337; m(7) * (t216 * t66 + t218 * t65 + t233 * t54); t218 * t2 + t216 * t1 + t233 * t11 + m(7) * (t54 ^ 2 + t65 ^ 2 + t66 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t83(1) t83(2) t83(4) t83(7) t83(11) t83(16); t83(2) t83(3) t83(5) t83(8) t83(12) t83(17); t83(4) t83(5) t83(6) t83(9) t83(13) t83(18); t83(7) t83(8) t83(9) t83(10) t83(14) t83(19); t83(11) t83(12) t83(13) t83(14) t83(15) t83(20); t83(16) t83(17) t83(18) t83(19) t83(20) t83(21);];
Mq  = res;
