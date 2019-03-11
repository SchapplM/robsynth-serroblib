% Calculate joint inertia matrix for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR13_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR13_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:12
% EndTime: 2019-03-09 11:24:23
% DurationCPUTime: 5.37s
% Computational Cost: add. (15106->606), mult. (33098->838), div. (0->0), fcn. (41649->12), ass. (0->275)
t264 = sin(pkin(6));
t266 = cos(pkin(6));
t269 = sin(qJ(2));
t271 = cos(qJ(2));
t222 = Icges(3,3) * t266 + (Icges(3,5) * t269 + Icges(3,6) * t271) * t264;
t223 = Icges(3,6) * t266 + (Icges(3,4) * t269 + Icges(3,2) * t271) * t264;
t224 = Icges(3,5) * t266 + (Icges(3,1) * t269 + Icges(3,4) * t271) * t264;
t225 = Icges(4,5) * t266 + (-Icges(4,6) * t269 - Icges(4,3) * t271) * t264;
t226 = Icges(4,4) * t266 + (-Icges(4,2) * t269 - Icges(4,6) * t271) * t264;
t227 = Icges(4,1) * t266 + (-Icges(4,4) * t269 - Icges(4,5) * t271) * t264;
t316 = t264 * t271;
t318 = t264 * t269;
t338 = (-t271 * t225 - t269 * t226) * t264 + t223 * t316 + t224 * t318 + (t227 + t222) * t266;
t299 = m(6) / 0.2e1 + m(7) / 0.2e1;
t337 = 0.2e1 * t299;
t270 = sin(qJ(1));
t312 = t270 * t271;
t272 = cos(qJ(1));
t313 = t269 * t272;
t244 = t266 * t312 + t313;
t268 = sin(qJ(4));
t317 = t264 * t270;
t325 = cos(qJ(4));
t211 = -t244 * t325 + t268 * t317;
t290 = t264 * t325;
t212 = t244 * t268 + t270 * t290;
t265 = cos(pkin(11));
t258 = pkin(5) * t265 + pkin(4);
t267 = -pkin(10) - qJ(5);
t311 = t271 * t272;
t314 = t269 * t270;
t245 = -t266 * t314 + t311;
t263 = sin(pkin(11));
t319 = t245 * t263;
t262 = pkin(11) + qJ(6);
t259 = sin(t262);
t260 = cos(t262);
t155 = -t212 * t259 + t245 * t260;
t156 = t212 * t260 + t245 * t259;
t94 = t156 * rSges(7,1) + t155 * rSges(7,2) + t211 * rSges(7,3);
t336 = pkin(5) * t319 - t211 * t267 + t212 * t258 + t94;
t242 = -t266 * t311 + t314;
t243 = t266 * t313 + t312;
t315 = t264 * t272;
t174 = -Icges(4,5) * t315 - Icges(4,6) * t243 + Icges(4,3) * t242;
t181 = Icges(3,4) * t243 - Icges(3,2) * t242 - Icges(3,6) * t315;
t335 = t174 - t181;
t176 = -Icges(4,4) * t315 - Icges(4,2) * t243 + Icges(4,6) * t242;
t183 = Icges(3,1) * t243 - Icges(3,4) * t242 - Icges(3,5) * t315;
t334 = t176 - t183;
t173 = Icges(4,5) * t317 - Icges(4,6) * t245 + Icges(4,3) * t244;
t182 = Icges(3,4) * t245 - Icges(3,2) * t244 + Icges(3,6) * t317;
t333 = -t182 + t173;
t175 = Icges(4,4) * t317 - Icges(4,2) * t245 + Icges(4,6) * t244;
t184 = Icges(3,1) * t245 - Icges(3,4) * t244 + Icges(3,5) * t317;
t332 = t184 - t175;
t331 = t264 ^ 2;
t213 = t242 * t325 + t268 * t315;
t240 = t266 * t268 + t271 * t290;
t84 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t211;
t86 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t211;
t88 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t211;
t28 = t155 * t86 + t156 * t88 + t211 * t84;
t215 = t242 * t268 - t272 * t290;
t157 = -t215 * t259 + t243 * t260;
t158 = t215 * t260 + t243 * t259;
t85 = Icges(7,5) * t158 + Icges(7,6) * t157 - Icges(7,3) * t213;
t87 = Icges(7,4) * t158 + Icges(7,2) * t157 - Icges(7,6) * t213;
t89 = Icges(7,1) * t158 + Icges(7,4) * t157 - Icges(7,5) * t213;
t29 = t155 * t87 + t156 * t89 + t211 * t85;
t241 = t266 * t325 - t268 * t316;
t197 = -t241 * t259 + t260 * t318;
t198 = t241 * t260 + t259 * t318;
t120 = Icges(7,5) * t198 + Icges(7,6) * t197 + Icges(7,3) * t240;
t121 = Icges(7,4) * t198 + Icges(7,2) * t197 + Icges(7,6) * t240;
t122 = Icges(7,1) * t198 + Icges(7,4) * t197 + Icges(7,5) * t240;
t44 = t120 * t211 + t121 * t155 + t122 * t156;
t1 = t211 * t28 - t213 * t29 + t240 * t44;
t330 = t1 / 0.2e1;
t37 = t197 * t86 + t198 * t88 + t240 * t84;
t38 = t197 * t87 + t198 * t89 + t240 * t85;
t53 = t240 * t120 + t197 * t121 + t198 * t122;
t48 = t53 * t240;
t11 = t37 * t211 - t38 * t213 + t48;
t329 = t11 / 0.2e1;
t328 = t211 / 0.2e1;
t327 = -t213 / 0.2e1;
t326 = t240 / 0.2e1;
t324 = t243 * pkin(2);
t323 = -pkin(4) + t258;
t149 = t212 * pkin(4) + qJ(5) * t211;
t322 = -t149 + t336;
t202 = t213 * qJ(5);
t320 = t243 * t263;
t283 = -rSges(7,1) * t158 - rSges(7,2) * t157;
t95 = -rSges(7,3) * t213 - t283;
t321 = pkin(5) * t320 + t213 * t267 + t215 * t323 + t202 + t95;
t164 = -t212 * t263 + t245 * t265;
t165 = t212 * t265 + t319;
t104 = t165 * rSges(6,1) + t164 * rSges(6,2) + t211 * rSges(6,3);
t310 = -t104 - t149;
t166 = -t215 * t263 + t243 * t265;
t167 = t215 * t265 + t320;
t105 = rSges(6,1) * t167 + rSges(6,2) * t166 - rSges(6,3) * t213;
t150 = pkin(4) * t215 - t202;
t309 = -t105 - t150;
t125 = rSges(7,1) * t198 + rSges(7,2) * t197 + rSges(7,3) * t240;
t296 = t263 * t318;
t308 = t125 + pkin(5) * t296 + t323 * t241 + (-qJ(5) - t267) * t240;
t307 = t338 * t266;
t178 = -Icges(4,1) * t315 - Icges(4,4) * t243 + Icges(4,5) * t242;
t179 = Icges(3,5) * t243 - Icges(3,6) * t242 - Icges(3,3) * t315;
t306 = -t179 - t178;
t177 = Icges(4,1) * t317 - Icges(4,4) * t245 + Icges(4,5) * t244;
t180 = Icges(3,5) * t245 - Icges(3,6) * t244 + Icges(3,3) * t317;
t305 = t180 + t177;
t232 = t242 * qJ(3);
t194 = t232 + t324;
t195 = t245 * pkin(2) + qJ(3) * t244;
t304 = t194 * t317 + t195 * t315;
t192 = t266 * t195;
t218 = pkin(3) * t317 + pkin(9) * t245;
t303 = t266 * t218 + t192;
t219 = -pkin(3) * t315 + t243 * pkin(9);
t302 = -t194 - t219;
t246 = (pkin(2) * t269 - qJ(3) * t271) * t264;
t301 = -pkin(3) * t266 - pkin(9) * t318 - t246;
t300 = t272 * pkin(1) + pkin(8) * t317;
t298 = t37 / 0.2e1 + t44 / 0.2e1;
t45 = -t120 * t213 + t121 * t157 + t122 * t158;
t297 = -t38 / 0.2e1 - t45 / 0.2e1;
t208 = -t241 * t263 + t265 * t318;
t209 = t241 * t265 + t296;
t129 = Icges(6,5) * t209 + Icges(6,6) * t208 + Icges(6,3) * t240;
t130 = Icges(6,4) * t209 + Icges(6,2) * t208 + Icges(6,6) * t240;
t131 = Icges(6,1) * t209 + Icges(6,4) * t208 + Icges(6,5) * t240;
t57 = t240 * t129 + t208 * t130 + t209 * t131;
t295 = t266 * t149 + t303;
t294 = -t150 + t302;
t169 = Icges(5,5) * t241 - Icges(5,6) * t240 + Icges(5,3) * t318;
t170 = Icges(5,4) * t241 - Icges(5,2) * t240 + Icges(5,6) * t318;
t171 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t318;
t79 = t169 * t318 - t240 * t170 + t241 * t171;
t193 = pkin(4) * t241 + qJ(5) * t240;
t293 = -t193 + t301;
t139 = t212 * rSges(5,1) - t211 * rSges(5,2) + t245 * rSges(5,3);
t188 = t245 * rSges(3,1) - t244 * rSges(3,2) + rSges(3,3) * t317;
t185 = rSges(4,1) * t317 - t245 * rSges(4,2) + t244 * rSges(4,3);
t289 = -t270 * pkin(1) + pkin(8) * t315;
t288 = t264 * (-rSges(4,1) * t266 - (-rSges(4,2) * t269 - rSges(4,3) * t271) * t264 - t246);
t287 = t218 * t315 + t219 * t317 + t304;
t286 = -t232 + t289;
t172 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t318;
t285 = t264 * (-t172 + t301);
t284 = -rSges(5,1) * t215 - rSges(5,2) * t213;
t132 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t240;
t282 = t264 * (-t132 + t293);
t281 = t195 + t300;
t280 = rSges(4,1) * t315 - t242 * rSges(4,3);
t279 = t149 * t315 + t150 * t317 + t287;
t278 = t286 - t219;
t277 = t264 * (t293 - t308);
t187 = t243 * rSges(3,1) - t242 * rSges(3,2) - rSges(3,3) * t315;
t100 = Icges(6,4) * t165 + Icges(6,2) * t164 + Icges(6,6) * t211;
t102 = Icges(6,1) * t165 + Icges(6,4) * t164 + Icges(6,5) * t211;
t98 = Icges(6,5) * t165 + Icges(6,6) * t164 + Icges(6,3) * t211;
t39 = t100 * t208 + t102 * t209 + t240 * t98;
t46 = t129 * t211 + t130 * t164 + t131 * t165;
t133 = Icges(5,5) * t212 - Icges(5,6) * t211 + Icges(5,3) * t245;
t135 = Icges(5,4) * t212 - Icges(5,2) * t211 + Icges(5,6) * t245;
t137 = Icges(5,1) * t212 - Icges(5,4) * t211 + Icges(5,5) * t245;
t67 = t133 * t318 - t135 * t240 + t137 * t241;
t73 = t169 * t245 - t170 * t211 + t171 * t212;
t275 = t73 / 0.2e1 + t46 / 0.2e1 + t67 / 0.2e1 + t39 / 0.2e1 + t298;
t101 = Icges(6,4) * t167 + Icges(6,2) * t166 - Icges(6,6) * t213;
t103 = Icges(6,1) * t167 + Icges(6,4) * t166 - Icges(6,5) * t213;
t99 = Icges(6,5) * t167 + Icges(6,6) * t166 - Icges(6,3) * t213;
t40 = t101 * t208 + t103 * t209 + t240 * t99;
t47 = -t129 * t213 + t130 * t166 + t131 * t167;
t134 = Icges(5,5) * t215 + Icges(5,6) * t213 + Icges(5,3) * t243;
t136 = Icges(5,4) * t215 + Icges(5,2) * t213 + Icges(5,6) * t243;
t138 = Icges(5,1) * t215 + Icges(5,4) * t213 + Icges(5,5) * t243;
t68 = t134 * t318 - t136 * t240 + t138 * t241;
t74 = t169 * t243 + t170 * t213 + t171 * t215;
t274 = t68 / 0.2e1 + t40 / 0.2e1 + t74 / 0.2e1 + t47 / 0.2e1 - t297;
t273 = t218 + t281;
t249 = rSges(2,1) * t272 - t270 * rSges(2,2);
t248 = -t270 * rSges(2,1) - rSges(2,2) * t272;
t228 = rSges(3,3) * t266 + (rSges(3,1) * t269 + rSges(3,2) * t271) * t264;
t186 = -t243 * rSges(4,2) - t280;
t163 = t243 * t193;
t160 = t188 + t300;
t159 = -t187 + t289;
t146 = t149 * t318;
t144 = -t266 * t187 - t228 * t315;
t143 = t188 * t266 - t228 * t317;
t141 = t245 * t150;
t140 = rSges(5,3) * t243 - t284;
t124 = t281 + t185;
t123 = (rSges(4,2) - pkin(2)) * t243 + t280 + t286;
t118 = (t187 * t270 + t188 * t272) * t264;
t117 = t222 * t317 - t223 * t244 + t224 * t245;
t116 = -t222 * t315 - t242 * t223 + t243 * t224;
t115 = t242 * t225 - t243 * t226 - t227 * t315;
t114 = t225 * t244 - t226 * t245 + t227 * t317;
t107 = (-t186 - t194) * t266 + t272 * t288;
t106 = t185 * t266 + t270 * t288 + t192;
t97 = t139 * t318 - t172 * t245;
t96 = -t140 * t318 + t172 * t243;
t93 = t178 * t266 + (-t174 * t271 - t176 * t269) * t264;
t92 = t177 * t266 + (-t173 * t271 - t175 * t269) * t264;
t91 = t180 * t266 + (t182 * t271 + t184 * t269) * t264;
t90 = t179 * t266 + (t181 * t271 + t183 * t269) * t264;
t81 = t273 + t139;
t80 = (-rSges(5,3) - pkin(2)) * t243 + t278 + t284;
t78 = t79 * t266;
t77 = t79 * t318;
t76 = (t185 * t272 + t186 * t270) * t264 + t304;
t75 = -t139 * t243 + t140 * t245;
t72 = (-t140 + t302) * t266 + t272 * t285;
t71 = t139 * t266 + t270 * t285 + t303;
t70 = t273 - t310;
t69 = t278 + t309 - t324;
t66 = -t125 * t213 - t240 * t95;
t65 = -t125 * t211 + t240 * t94;
t64 = t273 + t336;
t63 = -t215 * t258 + (-pkin(5) * t263 - pkin(2)) * t243 + (rSges(7,3) - t267) * t213 + t278 + t283;
t62 = t134 * t243 + t136 * t213 + t138 * t215;
t61 = t133 * t243 + t135 * t213 + t137 * t215;
t60 = t134 * t245 - t136 * t211 + t138 * t212;
t59 = t133 * t245 - t135 * t211 + t137 * t212;
t58 = (t139 * t272 + t140 * t270) * t264 + t287;
t56 = t57 * t266;
t55 = t57 * t318;
t54 = t211 * t95 + t213 * t94;
t52 = t53 * t266;
t51 = t53 * t318;
t50 = t104 * t318 + t146 + (-t132 - t193) * t245;
t49 = t132 * t243 + t309 * t318 + t163;
t43 = (-t105 + t294) * t266 + t272 * t282;
t42 = t104 * t266 + t270 * t282 + t295;
t41 = t105 * t245 + t310 * t243 + t141;
t36 = t101 * t166 + t103 * t167 - t213 * t99;
t35 = t100 * t166 + t102 * t167 - t213 * t98;
t34 = t101 * t164 + t103 * t165 + t211 * t99;
t33 = t100 * t164 + t102 * t165 + t211 * t98;
t32 = (t104 * t272 + t105 * t270) * t264 + t279;
t31 = t157 * t87 + t158 * t89 - t213 * t85;
t30 = t157 * t86 + t158 * t88 - t213 * t84;
t27 = t146 + t322 * t318 + (-t193 - t308) * t245;
t26 = t163 + t308 * t243 + (-t150 - t321) * t318;
t25 = (t294 - t321) * t266 + t272 * t277;
t24 = t266 * t322 + t270 * t277 + t295;
t23 = t78 + (t67 * t270 - t68 * t272) * t264;
t22 = t68 * t243 + t67 * t245 + t77;
t21 = t141 + t321 * t245 + (-t149 - t322) * t243;
t20 = t74 * t266 + (t270 * t61 - t272 * t62) * t264;
t19 = t73 * t266 + (t270 * t59 - t272 * t60) * t264;
t18 = t243 * t62 + t245 * t61 + t318 * t74;
t17 = t243 * t60 + t245 * t59 + t318 * t73;
t16 = (t270 * t321 + t272 * t322) * t264 + t279;
t15 = t56 + (t39 * t270 - t40 * t272) * t264;
t14 = t40 * t243 + t39 * t245 + t55;
t13 = t52 + (t37 * t270 - t38 * t272) * t264;
t12 = t38 * t243 + t37 * t245 + t51;
t10 = t47 * t266 + (t270 * t35 - t272 * t36) * t264;
t9 = t46 * t266 + (t270 * t33 - t272 * t34) * t264;
t8 = t243 * t36 + t245 * t35 + t318 * t47;
t7 = t243 * t34 + t245 * t33 + t318 * t46;
t6 = t45 * t266 + (t270 * t30 - t272 * t31) * t264;
t5 = t44 * t266 + (t270 * t28 - t272 * t29) * t264;
t4 = t243 * t31 + t245 * t30 + t318 * t45;
t3 = t243 * t29 + t245 * t28 + t318 * t44;
t2 = t211 * t30 - t213 * t31 + t240 * t45;
t82 = [Icges(2,3) + m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(4) * (t123 ^ 2 + t124 ^ 2) + m(3) * (t159 ^ 2 + t160 ^ 2) + m(2) * (t248 ^ 2 + t249 ^ 2) + t79 + t57 + t53 + t338; t52 + t78 + t56 + m(7) * (t24 * t64 + t25 * t63) + m(6) * (t42 * t70 + t43 * t69) + m(5) * (t71 * t81 + t72 * t80) + m(4) * (t106 * t124 + t107 * t123) + m(3) * (t143 * t160 + t144 * t159) + ((-t93 / 0.2e1 - t90 / 0.2e1 - t115 / 0.2e1 - t116 / 0.2e1 - t274) * t272 + (t92 / 0.2e1 + t91 / 0.2e1 + t114 / 0.2e1 + t117 / 0.2e1 + t275) * t270) * t264 + t307; (t13 + t23 + t15 + t307) * t266 + m(7) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t106 ^ 2 + t107 ^ 2 + t76 ^ 2) + m(3) * (t118 ^ 2 + t143 ^ 2 + t144 ^ 2) + ((-t10 - t20 - t6 + ((t242 * t335 - t243 * t334) * t264 + t306 * t331 * t272) * t272 + (-t115 - t116 - t90 - t93) * t266) * t272 + (t5 + t19 + t9 + ((t244 * t333 + t245 * t332) * t264 + t305 * t331 * t270) * t270 + (t117 + t114 + t92 + t91) * t266 + ((t270 * t306 + t305 * t272) * t264 + t334 * t245 - t335 * t244 - t332 * t243 - t333 * t242) * t315) * t270) * t264; m(7) * (t242 * t64 + t244 * t63) + m(6) * (t242 * t70 + t244 * t69) + m(5) * (t242 * t81 + t244 * t80) + m(4) * (t123 * t244 + t124 * t242); m(7) * (-t16 * t316 + t24 * t242 + t244 * t25) + m(6) * (t242 * t42 + t244 * t43 - t316 * t32) + m(5) * (t242 * t71 + t244 * t72 - t316 * t58) + m(4) * (t106 * t242 + t107 * t244 - t316 * t76); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t299) * (t271 ^ 2 * t331 + t242 ^ 2 + t244 ^ 2); t51 + t77 + t55 + m(7) * (t26 * t63 + t27 * t64) + m(6) * (t49 * t69 + t50 * t70) + m(5) * (t80 * t96 + t81 * t97) + t275 * t245 + t274 * t243; (t12 / 0.2e1 + t22 / 0.2e1 + t14 / 0.2e1) * t266 + (t5 / 0.2e1 + t19 / 0.2e1 + t9 / 0.2e1) * t245 + (t6 / 0.2e1 + t10 / 0.2e1 + t20 / 0.2e1) * t243 + m(7) * (t16 * t21 + t24 * t27 + t25 * t26) + m(6) * (t32 * t41 + t42 * t50 + t43 * t49) + m(5) * (t58 * t75 + t71 * t97 + t72 * t96) + ((-t4 / 0.2e1 - t18 / 0.2e1 - t8 / 0.2e1) * t272 + (t3 / 0.2e1 + t17 / 0.2e1 + t7 / 0.2e1) * t270 + (t13 / 0.2e1 + t23 / 0.2e1 + t15 / 0.2e1) * t269) * t264; m(5) * (t242 * t97 + t244 * t96 - t316 * t75) + m(6) * (t242 * t50 + t244 * t49 - t316 * t41) + m(7) * (-t21 * t316 + t242 * t27 + t244 * t26); (t12 + t14 + t22) * t318 + (t3 + t7 + t17) * t245 + (t4 + t18 + t8) * t243 + m(7) * (t21 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t41 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t75 ^ 2 + t96 ^ 2 + t97 ^ 2); m(7) * (t211 * t63 - t213 * t64) + m(6) * (t211 * t69 - t213 * t70); m(7) * (t16 * t240 + t211 * t25 - t213 * t24) + m(6) * (t211 * t43 - t213 * t42 + t240 * t32); (t211 * t244 - t213 * t242 - t240 * t316) * t337; m(7) * (t21 * t240 + t211 * t26 - t213 * t27) + m(6) * (t211 * t49 - t213 * t50 + t240 * t41); (t211 ^ 2 + t213 ^ 2 + t240 ^ 2) * t337; m(7) * (t63 * t66 + t64 * t65) + t48 + t297 * t213 + t298 * t211; m(7) * (t54 * t16 + t24 * t65 + t25 * t66) + t13 * t326 + t266 * t329 + t6 * t327 + t5 * t328 + (t270 * t330 - t272 * t2 / 0.2e1) * t264; m(7) * (t242 * t65 + t244 * t66 - t316 * t54); t318 * t329 + m(7) * (t54 * t21 + t26 * t66 + t27 * t65) + t245 * t330 + t243 * t2 / 0.2e1 + t3 * t328 + t12 * t326 + t4 * t327; m(7) * (t211 * t66 - t213 * t65 + t240 * t54); t211 * t1 - t213 * t2 + t240 * t11 + m(7) * (t54 ^ 2 + t65 ^ 2 + t66 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t82(1) t82(2) t82(4) t82(7) t82(11) t82(16); t82(2) t82(3) t82(5) t82(8) t82(12) t82(17); t82(4) t82(5) t82(6) t82(9) t82(13) t82(18); t82(7) t82(8) t82(9) t82(10) t82(14) t82(19); t82(11) t82(12) t82(13) t82(14) t82(15) t82(20); t82(16) t82(17) t82(18) t82(19) t82(20) t82(21);];
Mq  = res;
