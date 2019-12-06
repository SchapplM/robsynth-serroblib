% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:47
% EndTime: 2019-12-05 15:05:10
% DurationCPUTime: 14.28s
% Computational Cost: add. (14313->730), mult. (23213->1108), div. (0->0), fcn. (24727->10), ass. (0->314)
t257 = qJ(3) + pkin(9);
t254 = sin(t257);
t259 = sin(pkin(7));
t260 = cos(pkin(8));
t350 = t259 * t260;
t255 = cos(t257);
t261 = cos(pkin(7));
t355 = t255 * t261;
t206 = t254 * t350 + t355;
t258 = sin(pkin(8));
t317 = qJD(3) * t258;
t300 = t259 * t317;
t183 = qJD(5) * t206 + t300;
t378 = t183 / 0.2e1;
t347 = t261 * t254;
t208 = -t259 * t255 + t260 * t347;
t314 = qJD(3) * t261;
t299 = t258 * t314;
t184 = qJD(5) * t208 + t299;
t376 = t184 / 0.2e1;
t308 = qJD(5) * t258;
t316 = qJD(3) * t260;
t241 = t254 * t308 - t316;
t374 = t241 / 0.2e1;
t389 = qJD(5) / 0.2e1;
t192 = t206 * qJD(3);
t303 = t255 * t350;
t193 = qJD(3) * t303 - t254 * t314;
t108 = -Icges(5,5) * t192 - Icges(5,6) * t193;
t266 = cos(qJ(3));
t345 = t261 * t266;
t264 = sin(qJ(3));
t349 = t259 * t264;
t234 = -t260 * t349 - t345;
t220 = t234 * qJD(3);
t346 = t261 * t264;
t348 = t259 * t266;
t235 = t260 * t348 - t346;
t221 = t235 * qJD(3);
t145 = Icges(4,5) * t220 - Icges(4,6) * t221;
t388 = t108 + t145;
t194 = t208 * qJD(3);
t209 = t254 * t259 + t260 * t355;
t195 = t209 * qJD(3);
t109 = -Icges(5,5) * t194 - Icges(5,6) * t195;
t236 = -t260 * t346 + t348;
t222 = t236 * qJD(3);
t237 = t260 * t345 + t349;
t223 = t237 * qJD(3);
t146 = Icges(4,5) * t222 - Icges(4,6) * t223;
t387 = t109 + t146;
t386 = Icges(4,5) * t236 - Icges(5,5) * t208 - Icges(4,6) * t237 - Icges(5,6) * t209;
t207 = t303 - t347;
t385 = Icges(4,5) * t234 - Icges(5,5) * t206 - Icges(4,6) * t235 - Icges(5,6) * t207;
t212 = (-Icges(5,5) * t254 - Icges(5,6) * t255) * t258;
t238 = (-Icges(4,5) * t264 - Icges(4,6) * t266) * t258;
t384 = t212 + t238;
t196 = qJD(3) * t212;
t224 = qJD(3) * t238;
t383 = -t224 - t196;
t353 = t258 * t261;
t354 = t258 * t259;
t382 = -t260 * t384 + t353 * t386 + t354 * t385;
t265 = cos(qJ(5));
t263 = sin(qJ(5));
t352 = t258 * t263;
t228 = -t255 * t352 - t260 * t265;
t351 = t258 * t265;
t229 = t255 * t351 - t260 * t263;
t357 = t254 * t258;
t119 = Icges(6,5) * t229 + Icges(6,6) * t228 + Icges(6,3) * t357;
t358 = Icges(6,4) * t229;
t120 = Icges(6,2) * t228 + Icges(6,6) * t357 + t358;
t216 = Icges(6,4) * t228;
t121 = Icges(6,1) * t229 + Icges(6,5) * t357 + t216;
t173 = -t207 * t263 + t259 * t351;
t174 = t207 * t265 + t259 * t352;
t302 = t254 * t317;
t179 = -qJD(5) * t229 + t263 * t302;
t180 = qJD(5) * t228 - t265 * t302;
t301 = t255 * t317;
t75 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t301;
t76 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t301;
t77 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t301;
t97 = -qJD(5) * t174 + t192 * t263;
t98 = qJD(5) * t173 - t192 * t265;
t14 = t119 * t193 + t120 * t97 + t121 * t98 + t173 * t76 + t174 * t77 + t206 * t75;
t65 = Icges(6,5) * t174 + Icges(6,6) * t173 + Icges(6,3) * t206;
t360 = Icges(6,4) * t174;
t67 = Icges(6,2) * t173 + Icges(6,6) * t206 + t360;
t166 = Icges(6,4) * t173;
t69 = Icges(6,1) * t174 + Icges(6,5) * t206 + t166;
t20 = t173 * t67 + t174 * t69 + t206 * t65;
t175 = -t209 * t263 + t261 * t351;
t176 = t209 * t265 + t261 * t352;
t66 = Icges(6,5) * t176 + Icges(6,6) * t175 + Icges(6,3) * t208;
t359 = Icges(6,4) * t176;
t68 = Icges(6,2) * t175 + Icges(6,6) * t208 + t359;
t167 = Icges(6,4) * t175;
t70 = Icges(6,1) * t176 + Icges(6,5) * t208 + t167;
t21 = t173 * t68 + t174 * t70 + t206 * t66;
t287 = t193 * t20 + t195 * t21;
t40 = t119 * t206 + t120 * t173 + t121 * t174;
t52 = Icges(6,5) * t98 + Icges(6,6) * t97 + Icges(6,3) * t193;
t54 = Icges(6,4) * t98 + Icges(6,2) * t97 + Icges(6,6) * t193;
t56 = Icges(6,1) * t98 + Icges(6,4) * t97 + Icges(6,5) * t193;
t6 = t173 * t54 + t174 * t56 + t193 * t65 + t206 * t52 + t67 * t97 + t69 * t98;
t100 = qJD(5) * t175 - t194 * t265;
t99 = -qJD(5) * t176 + t194 * t263;
t53 = Icges(6,5) * t100 + Icges(6,6) * t99 + Icges(6,3) * t195;
t55 = Icges(6,4) * t100 + Icges(6,2) * t99 + Icges(6,6) * t195;
t57 = Icges(6,1) * t100 + Icges(6,4) * t99 + Icges(6,5) * t195;
t7 = t173 * t55 + t174 * t57 + t193 * t66 + t206 * t53 + t68 * t97 + t70 * t98;
t381 = t14 * t374 + t6 * t378 + t7 * t376 + (t301 * t40 + t287) * t389;
t15 = t100 * t121 + t119 * t195 + t120 * t99 + t175 * t76 + t176 * t77 + t208 * t75;
t22 = t175 * t67 + t176 * t69 + t208 * t65;
t23 = t175 * t68 + t176 * t70 + t208 * t66;
t286 = t22 * t193 + t23 * t195;
t41 = t119 * t208 + t120 * t175 + t121 * t176;
t8 = t100 * t69 + t175 * t54 + t176 * t56 + t195 * t65 + t208 * t52 + t67 * t99;
t9 = t100 * t70 + t175 * t55 + t176 * t57 + t195 * t66 + t208 * t53 + t68 * t99;
t380 = t15 * t374 + t8 * t378 + t9 * t376 + (t301 * t41 + t286) * t389;
t379 = -t183 / 0.2e1;
t377 = -t184 / 0.2e1;
t375 = -t241 / 0.2e1;
t372 = pkin(3) * t266;
t117 = -pkin(4) * t192 + pkin(6) * t193;
t58 = rSges(6,1) * t98 + rSges(6,2) * t97 + rSges(6,3) * t193;
t370 = t117 + t58;
t134 = pkin(4) * t207 + pkin(6) * t206;
t71 = rSges(6,1) * t174 + rSges(6,2) * t173 + rSges(6,3) * t206;
t369 = t134 + t71;
t368 = Icges(4,4) * t235;
t367 = Icges(4,4) * t237;
t366 = Icges(4,4) * t264;
t365 = Icges(4,4) * t266;
t364 = Icges(5,4) * t207;
t363 = Icges(5,4) * t209;
t362 = Icges(5,4) * t254;
t361 = Icges(5,4) * t255;
t356 = t255 * t258;
t101 = -Icges(5,2) * t206 + Icges(5,6) * t354 + t364;
t344 = -Icges(5,1) * t206 - t101 - t364;
t102 = -Icges(5,2) * t208 + Icges(5,6) * t353 + t363;
t343 = -Icges(5,1) * t208 - t102 - t363;
t201 = Icges(5,4) * t206;
t103 = Icges(5,1) * t207 + Icges(5,5) * t354 - t201;
t342 = Icges(5,2) * t207 - t103 + t201;
t202 = Icges(5,4) * t208;
t104 = Icges(5,1) * t209 + Icges(5,5) * t353 - t202;
t341 = Icges(5,2) * t209 - t104 + t202;
t270 = qJ(4) * t258 + t260 * t372;
t307 = pkin(3) * t349;
t142 = t261 * t270 + t307;
t340 = -rSges(5,1) * t209 + rSges(5,2) * t208 - rSges(5,3) * t353 - t142;
t313 = qJD(4) * t258;
t249 = t261 * t313;
t186 = pkin(3) * t222 + t249;
t339 = rSges(5,1) * t194 + rSges(5,2) * t195 - t186;
t338 = pkin(4) * t194 - pkin(6) * t195 - t186;
t141 = -pkin(3) * t346 + t259 * t270;
t188 = -qJ(4) * t260 + t258 * t372;
t337 = t260 * t141 + t188 * t354;
t233 = t236 * pkin(3);
t336 = rSges(5,1) * t208 + rSges(5,2) * t209 - t233;
t335 = pkin(4) * t208 - pkin(6) * t209 - t233;
t334 = -pkin(4) * t209 - pkin(6) * t208 - t142;
t137 = Icges(4,2) * t234 + Icges(4,6) * t354 + t368;
t333 = Icges(4,1) * t234 - t137 - t368;
t138 = Icges(4,2) * t236 + Icges(4,6) * t353 + t367;
t332 = Icges(4,1) * t236 - t138 - t367;
t230 = Icges(4,4) * t234;
t139 = Icges(4,1) * t235 + Icges(4,5) * t354 + t230;
t331 = -Icges(4,2) * t235 + t139 + t230;
t231 = Icges(4,4) * t236;
t140 = Icges(4,1) * t237 + Icges(4,5) * t353 + t231;
t330 = -Icges(4,2) * t237 + t140 + t231;
t248 = t259 * t313;
t185 = pkin(3) * t220 + t248;
t306 = pkin(3) * qJD(3) * t264;
t312 = qJD(4) * t260;
t243 = -t258 * t306 - t312;
t329 = t185 * t316 + t243 * t300;
t328 = t260 * t185 + t243 * t354;
t191 = -rSges(5,3) * t260 + (rSges(5,1) * t255 - rSges(5,2) * t254) * t258;
t327 = -t188 - t191;
t219 = (pkin(4) * t255 + pkin(6) * t254) * t258;
t326 = -t188 - t219;
t189 = -Icges(5,6) * t260 + (-Icges(5,2) * t254 + t361) * t258;
t214 = (-Icges(5,1) * t254 - t361) * t258;
t325 = t189 - t214;
t190 = -Icges(5,5) * t260 + (Icges(5,1) * t255 - t362) * t258;
t213 = (-Icges(5,2) * t255 - t362) * t258;
t324 = t190 + t213;
t272 = (-rSges(5,1) * t254 - rSges(5,2) * t255) * t258;
t200 = qJD(3) * t272;
t323 = -t200 - t243;
t274 = (-pkin(4) * t254 + pkin(6) * t255) * t258;
t203 = qJD(3) * t274;
t322 = -t203 - t243;
t210 = -Icges(4,6) * t260 + (-Icges(4,2) * t264 + t365) * t258;
t240 = (-Icges(4,1) * t264 - t365) * t258;
t321 = t210 - t240;
t211 = -Icges(4,5) * t260 + (Icges(4,1) * t266 - t366) * t258;
t239 = (-Icges(4,2) * t266 - t366) * t258;
t320 = t211 + t239;
t319 = qJD(2) * t261;
t318 = qJD(3) * t255;
t311 = qJD(5) * t193;
t310 = qJD(5) * t195;
t309 = qJD(5) * t255;
t59 = rSges(6,1) * t100 + rSges(6,2) * t99 + rSges(6,3) * t195;
t305 = -t59 + t338;
t72 = rSges(6,1) * t176 + rSges(6,2) * t175 + rSges(6,3) * t208;
t304 = -t72 + t334;
t297 = t311 / 0.2e1;
t296 = t310 / 0.2e1;
t256 = t258 ^ 2;
t294 = t256 * t307;
t293 = t248 - t319;
t253 = qJD(2) * t259;
t292 = t141 * t316 + t188 * t300 + t249 + t253;
t291 = t301 / 0.2e1;
t288 = -rSges(6,1) * t265 + rSges(6,2) * t263;
t30 = t228 * t67 + t229 * t69 + t357 * t65;
t31 = t228 * t68 + t229 * t70 + t357 * t66;
t285 = t193 * t30 + t195 * t31;
t284 = -t193 * t72 + t195 * t71;
t143 = rSges(4,1) * t235 + rSges(4,2) * t234 + rSges(4,3) * t354;
t215 = -t260 * rSges(4,3) + (rSges(4,1) * t266 - rSges(4,2) * t264) * t258;
t73 = t253 + (t143 * t260 + t215 * t354) * qJD(3);
t144 = rSges(4,1) * t237 + rSges(4,2) * t236 + rSges(4,3) * t353;
t74 = -t319 + (-t144 * t260 - t215 * t353) * qJD(3);
t283 = t73 * t259 - t74 * t261;
t151 = rSges(4,1) * t220 - rSges(4,2) * t221;
t273 = (-rSges(4,1) * t264 - rSges(4,2) * t266) * t258;
t227 = qJD(3) * t273;
t79 = (t151 * t260 + t227 * t354) * qJD(3);
t152 = rSges(4,1) * t222 - rSges(4,2) * t223;
t80 = (-t152 * t260 - t227 * t353) * qJD(3);
t282 = t259 * t79 - t261 * t80;
t281 = t141 * t299 + qJD(1) - t312;
t280 = -Icges(6,1) * t265 + Icges(6,4) * t263;
t279 = -Icges(6,4) * t265 + Icges(6,2) * t263;
t278 = -Icges(6,5) * t265 + Icges(6,6) * t263;
t277 = t143 * t261 - t144 * t259;
t276 = t151 * t261 - t152 * t259;
t275 = qJD(5) * t291;
t271 = (Icges(6,5) * t228 - Icges(6,6) * t229) * t241 + t183 * (Icges(6,5) * t173 - Icges(6,6) * t174) + t184 * (Icges(6,5) * t175 - Icges(6,6) * t176);
t269 = (Icges(6,1) * t175 - t359 - t68) * t184 + (Icges(6,1) * t173 - t360 - t67) * t183 + (Icges(6,1) * t228 - t120 - t358) * t241;
t268 = (-Icges(6,2) * t176 + t167 + t70) * t184 + (-Icges(6,2) * t174 + t166 + t69) * t183 + (-Icges(6,2) * t229 + t121 + t216) * t241;
t267 = (Icges(6,3) * t209 + t208 * t278 + t263 * t68 - t265 * t70) * t184 + (Icges(6,3) * t207 + t206 * t278 + t263 * t67 - t265 * t69) * t183 + (t120 * t263 - t121 * t265 + (Icges(6,3) * t255 + t254 * t278) * t258) * t241;
t244 = t261 * t256 * t306;
t232 = t234 * pkin(3);
t226 = qJD(3) * t240;
t225 = qJD(3) * t239;
t199 = t232 * t316;
t198 = qJD(3) * t214;
t197 = qJD(3) * t213;
t187 = t232 * t299;
t172 = (rSges(6,3) * t255 + t254 * t288) * t258;
t171 = (Icges(6,5) * t255 + t254 * t280) * t258;
t170 = (Icges(6,6) * t255 + t254 * t279) * t258;
t168 = t185 * t353;
t165 = rSges(4,1) * t236 - rSges(4,2) * t237;
t164 = rSges(4,1) * t234 - rSges(4,2) * t235;
t163 = t185 * t299;
t156 = rSges(6,1) * t228 - rSges(6,2) * t229;
t150 = Icges(4,1) * t222 - Icges(4,4) * t223;
t149 = Icges(4,1) * t220 - Icges(4,4) * t221;
t148 = Icges(4,4) * t222 - Icges(4,2) * t223;
t147 = Icges(4,4) * t220 - Icges(4,2) * t221;
t133 = -pkin(4) * t206 + pkin(6) * t207;
t131 = -rSges(5,1) * t206 - rSges(5,2) * t207;
t122 = rSges(6,1) * t229 + rSges(6,2) * t228 + rSges(6,3) * t357;
t116 = t141 * t353;
t114 = -rSges(5,1) * t192 - rSges(5,2) * t193;
t113 = -Icges(5,1) * t194 - Icges(5,4) * t195;
t112 = -Icges(5,1) * t192 - Icges(5,4) * t193;
t111 = -Icges(5,4) * t194 - Icges(5,2) * t195;
t110 = -Icges(5,4) * t192 - Icges(5,2) * t193;
t105 = rSges(5,1) * t207 - rSges(5,2) * t206 + rSges(5,3) * t354;
t96 = rSges(6,1) * t175 - rSges(6,2) * t176;
t95 = rSges(6,1) * t173 - rSges(6,2) * t174;
t88 = rSges(6,3) * t209 + t208 * t288;
t87 = rSges(6,3) * t207 + t206 * t288;
t86 = Icges(6,5) * t209 + t208 * t280;
t85 = Icges(6,5) * t207 + t206 * t280;
t84 = Icges(6,6) * t209 + t208 * t279;
t83 = Icges(6,6) * t207 + t206 * t279;
t78 = rSges(6,1) * t180 + rSges(6,2) * t179 + rSges(6,3) * t301;
t64 = t276 * t317;
t63 = t277 * t317 + qJD(1);
t61 = (t339 * t260 + t323 * t353) * qJD(3);
t60 = (t114 * t260 + t200 * t354) * qJD(3) + t329;
t48 = (t340 * t260 + t327 * t353) * qJD(3) + t293;
t47 = (t105 * t260 + t191 * t354) * qJD(3) + t292;
t46 = t163 + (t114 * t261 + t339 * t259) * t317;
t45 = t119 * t357 + t120 * t228 + t121 * t229;
t42 = (t105 * t261 + t340 * t259) * t317 + t281;
t25 = -t122 * t184 + t241 * t72 + (t260 * t334 + t326 * t353) * qJD(3) + t293;
t24 = t122 * t183 - t241 * t71 + (t134 * t260 + t219 * t354) * qJD(3) + t292;
t19 = -t183 * t72 + t184 * t71 + (t134 * t261 + t259 * t334) * t317 + t281;
t18 = t120 * t179 + t121 * t180 + t228 * t76 + t229 * t77 + (t119 * t318 + t254 * t75) * t258;
t17 = -t122 * t310 - t184 * t78 + t241 * t59 + (t338 * t260 + (t261 * t322 + t309 * t72) * t258) * qJD(3);
t16 = t122 * t311 + t183 * t78 - t241 * t58 + (t117 * t260 + (t203 * t259 - t309 * t71) * t258) * qJD(3) + t329;
t13 = -t183 * t59 + t184 * t58 + t163 + t284 * qJD(5) + (t117 * t261 + t259 * t338) * t317;
t12 = t179 * t68 + t180 * t70 + t228 * t55 + t229 * t57 + (t254 * t53 + t318 * t66) * t258;
t11 = t179 * t67 + t180 * t69 + t228 * t54 + t229 * t56 + (t254 * t52 + t318 * t65) * t258;
t10 = t183 * t30 + t184 * t31 + t241 * t45;
t5 = t183 * t22 + t23 * t184 + t241 * t41;
t4 = t20 * t183 + t184 * t21 + t241 * t40;
t3 = t11 * t183 + t12 * t184 + t241 * t18 + (t301 * t45 + t285) * qJD(5);
t1 = [m(4) * t64 + m(5) * t46 + m(6) * t13; m(4) * t282 + m(5) * (t259 * t60 - t261 * t61) + m(6) * (t16 * t259 - t17 * t261); (-t18 * t260 + (t11 * t259 + t12 * t261) * t258) * t374 + ((t228 * t84 + t229 * t86) * t184 + (t228 * t83 + t229 * t85) * t183 + (t170 * t228 + t171 * t229) * t241 + (t207 * t30 + t209 * t31) * qJD(5) + ((qJD(5) * t45 + t119 * t241 + t183 * t65 + t184 * t66) * t255 + t267 * t254) * t258) * t375 + (-t15 * t260 + (t259 * t8 + t261 * t9) * t258) * t376 + ((t175 * t84 + t176 * t86 + t209 * t66) * t184 + (t175 * t83 + t176 * t85 + t209 * t65) * t183 + (t119 * t209 + t170 * t175 + t171 * t176) * t241 + (t207 * t22 + t209 * t23 + t356 * t41) * qJD(5) + t267 * t208) * t377 + (-t14 * t260 + (t259 * t6 + t261 * t7) * t258) * t378 + ((t173 * t84 + t174 * t86 + t207 * t66) * t184 + (t173 * t83 + t174 * t85 + t207 * t65) * t183 + (t119 * t207 + t170 * t173 + t171 * t174) * t241 + (t20 * t207 + t209 * t21 + t356 * t40) * qJD(5) + t267 * t206) * t379 + t353 * t380 + t354 * t381 + (-t260 * t41 + (t22 * t259 + t23 * t261) * t258) * t296 + (-t260 * t40 + (t20 * t259 + t21 * t261) * t258) * t297 - t255 * t10 * t308 / 0.2e1 + (-t260 * t45 + (t259 * t30 + t261 * t31) * t258) * t275 - t260 * t3 / 0.2e1 - (t207 * t4 + t209 * t5) * qJD(5) / 0.2e1 + (t260 * (-t260 * t224 + (-t225 * t264 + t226 * t266 + (-t210 * t266 - t211 * t264) * qJD(3)) * t258) - (t259 * (-t260 * t145 + (-t147 * t264 + t149 * t266 + (-t137 * t266 - t139 * t264) * qJD(3)) * t258) + t261 * (-t260 * t146 + (-t148 * t264 + t150 * t266 + (-t138 * t266 - t140 * t264) * qJD(3)) * t258)) * t258 + t260 * (-t196 * t260 + (-t197 * t254 + t198 * t255 + (-t189 * t255 - t190 * t254) * qJD(3)) * t258) - (t259 * (-t108 * t260 + (-t110 * t254 + t112 * t255 + (-t101 * t255 - t103 * t254) * qJD(3)) * t258) + t261 * (-t109 * t260 + (-t111 * t254 + t113 * t255 + (-t102 * t255 - t104 * t254) * qJD(3)) * t258)) * t258) * t316 + ((t189 * t193 + t190 * t192 + t197 * t206 - t198 * t207 + t210 * t221 - t211 * t220 - t225 * t234 - t226 * t235 + t354 * t383) * t260 + ((-t102 * t193 - t104 * t192 - t111 * t206 + t113 * t207 - t138 * t221 + t140 * t220 + t148 * t234 + t150 * t235 + t354 * t387) * t261 + (-t101 * t193 - t103 * t192 - t110 * t206 + t112 * t207 - t137 * t221 + t139 * t220 + t147 * t234 + t149 * t235 + t354 * t388) * t259) * t258) * t300 + ((t189 * t195 + t190 * t194 + t197 * t208 - t198 * t209 + t210 * t223 - t211 * t222 - t225 * t236 - t226 * t237 + t353 * t383) * t260 + ((-t102 * t195 - t104 * t194 - t111 * t208 + t113 * t209 - t138 * t223 + t140 * t222 + t148 * t236 + t150 * t237 + t353 * t387) * t261 + (-t101 * t195 - t103 * t194 - t110 * t208 + t112 * t209 - t137 * t223 + t139 * t222 + t147 * t236 + t149 * t237 + t353 * t388) * t259) * t258) * t299 + (((t254 * t324 + t255 * t325 - t385 * t259 - t261 * t386 + t264 * t320 + t266 * t321) * t260 + ((t259 * t344 + t261 * t343) * t255 + (t259 * t342 + t261 * t341) * t254 + (t259 * t333 + t261 * t332) * t266 + (-t259 * t331 - t261 * t330) * t264) * t258) * t317 + t384 * qJD(3) * t260 ^ 2) * t316 / 0.2e1 + (t16 * t337 + t24 * t328 + t13 * t116 + t19 * t168 + (t16 * t369 + t17 * t304 + t24 * t370 + t25 * t305) * t260 + ((t17 * (-t122 + t326) + t25 * (-t78 + t322) + t13 * t369 + t19 * t370) * t261 + (t16 * (t122 + t219) + t24 * (t203 + t78) + t13 * t304 + t19 * t305) * t259) * t258 - t24 * (t172 * t183 - t241 * t87 + t199) - t25 * (-t172 * t184 + t241 * t88 + t244) - t19 * (-t183 * t88 + t184 * t87 + t187) - (t24 * (t122 * t207 - t356 * t71) + t25 * (-t122 * t209 + t356 * t72) + t19 * (-t207 * t72 + t209 * t71)) * qJD(5) - (-t24 * t294 + (t24 * t133 + t25 * t335) * t260 + (t19 * (t133 * t261 + t259 * t335) + (t24 * t259 - t25 * t261) * t274) * t258) * qJD(3)) * m(6) + (t60 * t337 + t47 * t328 + t46 * t116 + t42 * t168 + (t60 * t105 + t47 * t114 + t339 * t48 + t340 * t61) * t260 + ((t46 * t105 + t42 * t114 + t323 * t48 + t327 * t61) * t261 + (t60 * t191 + t47 * t200 + t339 * t42 + t340 * t46) * t259) * t258 - t42 * t187 - t47 * t199 - t48 * t244 - (-t47 * t294 + (t47 * t131 + t336 * t48) * t260 + (t42 * (t131 * t261 + t259 * t336) + (t47 * t259 - t48 * t261) * t272) * t258) * qJD(3)) * m(5) + (-((t164 * t73 - t165 * t74) * t260 + (t63 * (t164 * t261 - t165 * t259) + t283 * t273) * t258) * qJD(3) + (t143 * t79 - t144 * t80 + t151 * t73 - t152 * t74) * t260 + (t215 * t282 + t227 * t283 + t276 * t63 + t277 * t64) * t258) * m(4) - (((t206 * t341 + t207 * t343 + t234 * t330 + t235 * t332) * t353 + (t206 * t324 + t207 * t325 - t234 * t320 + t235 * t321) * t260 + (t206 * t342 + t207 * t344 + t234 * t331 + t235 * t333 + t382) * t354) * t259 + ((t208 * t342 + t209 * t344 + t236 * t331 + t237 * t333) * t354 + (t208 * t324 + t209 * t325 - t236 * t320 + t237 * t321) * t260 + (t208 * t341 + t209 * t343 + t236 * t330 + t237 * t332 + t382) * t353) * t261) * t258 * qJD(3) ^ 2 / 0.2e1; m(5) * (-t260 * t46 + (t259 * t61 + t261 * t60) * t258) + m(6) * (-t13 * t260 + (t16 * t261 + t17 * t259) * t258); t195 * t5 / 0.2e1 + t208 * t380 + (t206 * t22 + t208 * t23 + t357 * t41) * t296 + (t8 * t206 + t9 * t208 + (t15 * t254 + t318 * t41) * t258 + t286) * t376 + t193 * t4 / 0.2e1 + t206 * t381 + (t20 * t206 + t208 * t21 + t357 * t40) * t297 + (t206 * t6 + t208 * t7 + (t14 * t254 + t318 * t40) * t258 + t287) * t378 + t10 * t291 + t3 * t357 / 0.2e1 + (t206 * t30 + t208 * t31 + t357 * t45) * t275 + (t11 * t206 + t12 * t208 + (t18 * t254 + t318 * t45) * t258 + t285) * t374 + (t175 * t268 + t176 * t269 + t208 * t271) * t377 + (t173 * t268 + t174 * t269 + t206 * t271) * t379 + (t228 * t268 + t229 * t269 + t271 * t357) * t375 + (t13 * (-t206 * t72 + t208 * t71) + (t206 * t24 - t208 * t25) * t78 + (t16 * t206 - t17 * t208 + t193 * t24 - t195 * t25) * t122 + ((-t24 * t71 + t25 * t72) * t318 + (-t16 * t71 + t17 * t72 - t24 * t58 + t25 * t59) * t254) * t258 - t24 * (t156 * t183 - t241 * t95) - t25 * (-t156 * t184 + t241 * t96) + (t183 * t96 - t184 * t95 - t206 * t59 + t208 * t58 + t284) * t19) * m(6);];
tauc = t1(:);
